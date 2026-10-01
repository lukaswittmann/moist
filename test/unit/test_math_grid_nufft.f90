!> Transform-level tests for the atom-centered molecular 3D grid (FINUFFT NUFFT)
!>
!> Analytic checks of the pieces a 3D-RISM UV k-space solve is built from
!>
!> - Round-trip diagonal, forward vs analytic FT, k-space convolution
!> - Quadrature refinement, type-1/2 vs type-3 routes, stale and recreated trafo
!> - Tolerances assume a nucleus-localised field, the cloud resolves `2*pi/k` only near the nuclei
module test_math_grid_nufft
   use, intrinsic :: iso_fortran_env, only: output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_env, only: wp, mctc_error_type => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use moist_math_grid_3d_base, only: moist_math_grid_3d_trafo_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, &
      & new_cartesian_grid_3d
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type, &
      & moist_math_grid_atomic_recipe_override_type
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, &
      & new_molecular_grid, molecular_grid_set_kgrid, &
      & moist_math_grid_3d_molecular_trafo_type, new_molecular_grid_trafo, &
      & default_nufft_tol
   use test_helpers, only: get_uniform_recipe

   implicit none(type, external)
   private

   public :: collect_math_grid_nufft

   !> Reciprocal resolution (bohr), k_max = pi/dr, as in the 3D-UV molecular runs
   real(wp), parameter :: probe_dr = 0.70_wp
   !> Outer radial clamp of the probe grid (bohr)
   real(wp), parameter :: probe_rmax = 16.0_wp
   !> Radial shells per atom of the probe grid
   integer, parameter :: probe_nrad = 50
   !> Raw Lebedev point count per shell of the probe grid
   integer, parameter :: probe_nang = 110
   !> Exponent of the nucleus-localised probe Gaussian (bohr^-2), 1/e width 1 bohr
   real(wp), parameter :: probe_alpha = 1.0_wp
   !> Acceptance bound for the type-1/2 vs type-3 comparison
   !>
   !> - 100 times `default_nufft_tol`, headroom for two independent approximations
   !> - A convention bug lands at O(1), caught at any tolerance
   real(wp), parameter :: route_agreement_tol = 100.0_wp*default_nufft_tol

contains

   !> Collect the molecular-NUFFT transform tests
   !>
   !> @param[out] testsuite  registered unit tests
   subroutine collect_math_grid_nufft(testsuite)
      !> Registered unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("nufft_roundtrip_diagonal", test_roundtrip_diagonal), &
                  new_unittest("nufft_forward_matches_analytic_ft", test_forward_analytic), &
                  new_unittest("nufft_convolution_matches_analytic", test_convolution), &
                  new_unittest("nufft_quadrature_refines", test_quadrature_refines), &
                  new_unittest("nufft_type12_matches_type3", test_type12_vs_type3_molecular), &
                  new_unittest("nufft_type12_matches_type3_uniform", test_type12_vs_type3_uniform), &
                  new_unittest("nufft_trafo_stale_after_update", test_trafo_stale_after_update), &
                  new_unittest("nufft_trafo_stale_after_reconstruction", test_trafo_stale_after_reconstruction), &
                  new_unittest("nufft_trafo_destroy_recreate", test_trafo_destroy_recreate) &
                  ]
   end subroutine collect_math_grid_nufft

   !> Water geometry (bohr) as carrier solute of every probe grid
   !>
   !> @param[out] mol  three-site water structure
   subroutine probe_structure(mol)
      !> Three-site water structure
      type(structure_type), intent(out) :: mol

      call new(mol, [8, 1, 1], reshape([ &
                                       0.0_wp, 0.0_wp, 0.0_wp, &
                                       1.43_wp, 0.0_wp, 1.11_wp, &
                                       -1.43_wp, 0.0_wp, 1.11_wp], [3, 3]))
   end subroutine probe_structure

   !> Build the shared probe grid and its prepared single-column transform
   !>
   !> - `ready` false with `error` set on failure, callers release the grid either way
   !>
   !> @param[out]   mol    carrier structure
   !> @param[out]   mg     molecular grid with its reciprocal grid configured
   !> @param[out]   trafo  transform prepared for one column
   !> @param[out]   ready  whether a usable transform was produced
   !> @param[out]   error  propagated test failure
   subroutine setup_probe(mol, mg, trafo, ready, error)
      !> Carrier structure
      type(structure_type), intent(out) :: mol
      !> Molecular grid
      type(moist_math_grid_3d_molecular_type), intent(out), target :: mg
      !> Prepared transform
      type(moist_math_grid_3d_molecular_trafo_type), intent(out) :: trafo
      !> Whether a usable transform was produced
      logical, intent(out) :: ready
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error_type), allocatable :: merr
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)

      ready = .false.
      call probe_structure(mol)
      call get_uniform_recipe(recipe, overrides, probe_nrad, probe_nang, merr, rmax=probe_rmax)
      if (.not. allocated(merr)) call new_molecular_grid(mg, merr, recipe=recipe, overrides=overrides, &
                                                         reciprocal=.false.)
      if (.not. allocated(merr)) call mg%update(mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call molecular_grid_set_kgrid(mg, probe_dr, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call new_molecular_grid_trafo(trafo, mg)
      call trafo%prepare(1, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      ready = .true.
   end subroutine setup_probe

   !> Forward transform through the abstract engine, failing the test on error
   !>
   !> @param[out]    error  test error
   !> @param[in,out] trafo  transform engine
   !> @param[in,out] f_r    real-space block
   !> @param[out]    f_k    reciprocal-space block
   subroutine forward(error, trafo, f_r, f_k)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Transform engine
      class(moist_math_grid_3d_trafo_type), intent(inout) :: trafo
      !> Real-space block
      real(wp), intent(inout), contiguous, target :: f_r(:, :)
      !> Reciprocal-space block
      complex(wp), intent(out), contiguous, target :: f_k(:, :)

      type(mctc_error_type), allocatable :: merr

      call trafo%fft_r2k(f_r, f_k, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine forward

   !> Backward transform through the abstract engine, failing the test on error
   !>
   !> @param[out]    error  test error
   !> @param[in,out] trafo  transform engine
   !> @param[in,out] f_k    reciprocal-space block (destroyed)
   !> @param[out]    f_r    real-space block
   subroutine backward(error, trafo, f_k, f_r)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Transform engine
      class(moist_math_grid_3d_trafo_type), intent(inout) :: trafo
      !> Reciprocal-space block (destroyed)
      complex(wp), intent(inout), contiguous, target :: f_k(:, :)
      !> Real-space block
      real(wp), intent(out), contiguous, target :: f_r(:, :)

      type(mctc_error_type), allocatable :: merr

      call trafo%fft_k2r(f_k, f_r, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine backward

   !> Check forward-then-backward transform of a unit impulse against `w_i/dr^3`
   !>
   !> - Composite operator `g_j = (dkx*dky*dkz/(2*pi)^3) * sum_k sum_i w_i f_i exp(i k.(r_j - r_i))`
   !> - At `j = i` the k sum is `nkx*nky*nkz`
   !> - Hence `g_i = w_i/dr^3 * f_i`, no other point contributes
   !>
   !> @param[out] error  propagated test failure
   subroutine test_roundtrip_diagonal(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      real(wp), allocatable :: f(:, :), g(:, :)
      complex(wp), allocatable :: fk(:, :)
      real(wp) :: dv, expected
      logical :: ready
      integer :: isel, ip, probes(3)

      call setup_probe(mol, mg, trafo, ready, error)
      if (allocated(error) .or. .not. ready) then
         call mg%destroy()
         return
      end if

      dv = probe_dr**3
      allocate (f(mg%ngrid, 1), g(mg%ngrid, 1), fk(mg%npts_k, 1))

      ! Sample the extremes and the middle of the weight distribution
      probes = [1, mg%ngrid/2, maxloc(mg%w, dim=1)]
      do ip = 1, size(probes)
         isel = probes(ip)
         f = 0.0_wp
         f(isel, 1) = 1.0_wp
         call forward(error, trafo, f, fk)
         if (allocated(error)) exit
         call backward(error, trafo, fk, g)
         if (allocated(error)) exit
         expected = mg%w(isel)/dv
         call check(error, abs(g(isel, 1) - expected) <= 1.0e-8_wp*abs(expected) + 1.0e-12_wp, &
                    "molecular NUFFT round-trip diagonal is not w_i/dr^3")
         if (allocated(error)) exit
      end do

      call trafo%destroy()
      call mg%destroy()
   end subroutine test_roundtrip_diagonal

   !> Check forward transform of a nucleus-localised Gaussian against the analytic FT
   !>
   !> - `f(r) = exp(-a |r - c|^2)` with `r0 = point(1)`
   !> - `F(k) = (pi/a)^(3/2) exp(-k^2/(4a)) exp(-i k.(c - r0))`
   !> - Pins measure, sign, phase reference and k ordering of both backends
   !>
   !> @param[out] error  propagated test failure
   subroutine test_forward_analytic(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      type(moist_math_grid_3d_cartesian_type), target :: cg
      class(moist_math_grid_3d_trafo_type), allocatable :: ctrafo
      type(mctc_error_type), allocatable :: merr
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)
      real(wp) :: cen(3), emax
      logical :: ready
      integer :: j

      call setup_probe(mol, mg, trafo, ready, error)
      if (allocated(error) .or. .not. ready) then
         call mg%destroy()
         return
      end if
      cen = mol%xyz(:, 1)

      allocate (f(mg%ngrid, 1), fk(mg%npts_k, 1))
      do j = 1, mg%ngrid
         f(j, 1) = gaussian(mg%xyz(:, j), cen, probe_alpha)
      end do
      call forward(error, trafo, f, fk)
      if (allocated(error)) then
         call trafo%destroy()
         call mg%destroy()
         return
      end if
      ! Phase against mg%kref (point(1) directly, solute centroid after a domain update)
      call max_ft_error(fk(:, 1), mg%npts_k, mg%kref, cen, mg, .true., cg, emax)
      ! Quadrature-limited: 3.0e-4 at (50, 110), 2.7e-7 from Lebedev order 302
      call check(error, emax <= 1.0e-3_wp, &
                 "molecular forward NUFFT disagrees with the analytic Gaussian FT")
      call trafo%destroy()
      call mg%destroy()
      deallocate (f, fk)
      if (allocated(error)) return

      call new_cartesian_grid_3d(cg, 48, 48, 48, probe_dr, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call cg%new_trafo(ctrafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      allocate (f(cg%ngrid, 1), fk(cg%npts_k, 1))
      do j = 1, cg%ngrid
         f(j, 1) = gaussian(cg%point(j), cen, probe_alpha)
      end do
      call forward(error, ctrafo, f, fk)
      if (allocated(error)) then
         call ctrafo%destroy()
         call cg%destroy()
         return
      end if
      call max_ft_error(fk(:, 1), cg%npts_k, cg%point(1), cen, mg, .false., cg, emax)
      ! Midpoint-rule limited near Nyquist, 6.5e-3 measured
      call check(error, emax <= 2.0e-2_wp, &
                 "Cartesian forward FFT disagrees with the analytic Gaussian FT")

      call ctrafo%destroy()
      call cg%destroy()
   end subroutine test_forward_analytic

   !> Find the largest deviation `|F_num(k) - F_exact(k)| / F_exact(0)` in the band `|k| <= pi/dr`
   !>
   !> - Only the inscribed ball of the cube grid is resolved
   !> - Full cube gives 3.7e-3 instead of 3.0e-4
   !> - `use_molecular` selects the grid providing k-points, the other is ignored
   !>
   !> @param[in]  fk             numerical transform, length nk
   !> @param[in]  nk             number of reciprocal points
   !> @param[in]  kref           phase reference of the transform (bohr)
   !> @param[in]  cen            gaussian center (bohr)
   !> @param[in]  mg             molecular grid supplying kpoint (if selected)
   !> @param[in]  use_molecular  whether to read k-points from `mg`
   !> @param[in]  cg             cartesian grid supplying kpoint (if selected)
   !> @param[out] emax           largest relative deviation
   subroutine max_ft_error(fk, nk, kref, cen, mg, use_molecular, cg, emax)
      !> Numerical transform
      complex(wp), intent(in) :: fk(:)
      !> Number of reciprocal points
      integer, intent(in) :: nk
      !> Phase reference of the transform
      real(wp), intent(in) :: kref(3)
      !> Gaussian center
      real(wp), intent(in) :: cen(3)
      !> Molecular grid
      type(moist_math_grid_3d_molecular_type), intent(in) :: mg
      !> Whether to read k-points from `mg`
      logical, intent(in) :: use_molecular
      !> Cartesian grid
      type(moist_math_grid_3d_cartesian_type), intent(in) :: cg
      !> Largest relative deviation
      real(wp), intent(out) :: emax

      real(wp) :: kvec(3), f0, phase, kband, dev
      complex(wp) :: fex
      integer :: j

      f0 = (pi/probe_alpha)**1.5_wp
      kband = pi/probe_dr
      emax = 0.0_wp
      do j = 1, nk
         if (use_molecular) then
            kvec = mg%kpoint(j)
         else
            kvec = cg%kpoint(j)
         end if
         if (sqrt(sum(kvec**2)) > kband) cycle
         phase = dot_product(kvec, cen - kref)
         fex = f0*exp(-sum(kvec**2)/(4.0_wp*probe_alpha)) &
               *cmplx(cos(phase), -sin(phase), wp)
         dev = abs(fk(j) - fex)/f0
         ! Fail on a non-finite deviation, NaN does not propagate reliably through max()
         if (.not. ieee_is_finite(dev)) then
            emax = huge(1.0_wp)
            return
         end if
         emax = max(emax, dev)
      end do
   end subroutine max_ft_error

   !> Check the k-space convolution route of a 3D-RISM UV step against an exact answer
   !>
   !> - Forward transform `exp(-a1 |r - c|^2)`, multiply by the FT of `exp(-a2 |r|^2)`
   !> - Backward transform gives the convolution
   !> - Result `(pi/(a1+a2))^(3/2) exp(-a1 a2/(a1+a2)|r-c|^2)`
   !>
   !> @param[out] error  propagated test failure
   subroutine test_convolution(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Width parameter of the convolution kernel (bohr^-2)
      real(wp), parameter :: a2 = 0.5_wp
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      real(wp), allocatable :: f(:, :), g(:, :)
      complex(wp), allocatable :: fk(:, :)
      real(wp) :: cen(3), kvec(3), ab, peak, emax, exact, dev
      logical :: ready
      integer :: j

      call setup_probe(mol, mg, trafo, ready, error)
      if (allocated(error) .or. .not. ready) then
         call mg%destroy()
         return
      end if
      cen = mol%xyz(:, 1)
      ab = probe_alpha*a2/(probe_alpha + a2)
      peak = (pi/(probe_alpha + a2))**1.5_wp

      allocate (f(mg%ngrid, 1), g(mg%ngrid, 1), fk(mg%npts_k, 1))
      do j = 1, mg%ngrid
         f(j, 1) = gaussian(mg%xyz(:, j), cen, probe_alpha)
      end do
      call forward(error, trafo, f, fk)
      if (allocated(error)) then
         call trafo%destroy()
         call mg%destroy()
         return
      end if
      do j = 1, mg%npts_k
         kvec = mg%kpoint(j)
         fk(j, 1) = fk(j, 1)*(pi/a2)**1.5_wp*exp(-sum(kvec**2)/(4.0_wp*a2))
      end do
      call backward(error, trafo, fk, g)
      if (allocated(error)) then
         call trafo%destroy()
         call mg%destroy()
         return
      end if

      emax = 0.0_wp
      do j = 1, mg%ngrid
         exact = peak*exp(-ab*sum((mg%xyz(:, j) - cen)**2))
         dev = abs(g(j, 1) - exact)
         ! Fail on a non-finite deviation, NaN does not propagate reliably through max()
         if (.not. ieee_is_finite(dev)) then
            emax = huge(1.0_wp)
            exit
         end if
         emax = max(emax, dev)
      end do
      ! Measured 2.5e-5 against a peak of 3.03
      call check(error, emax <= 1.0e-3_wp*peak, &
                 "molecular NUFFT convolution disagrees with the analytic Gaussian")

      call trafo%destroy()
      call mg%destroy()
   end subroutine test_convolution

   !> Check quadrature convergence for a diffuse, shell-shaped field
   !>
   !> - `integrate_field` of `exp(-(|r - c| - R0)^2)`, `R0 = 6` bohr
   !> - Exact value `4*pi*integral r^2 exp(-(r-R0)^2) dr`
   !> - Sparse outer region, so a transform failure is not the quadrature rule
   !>
   !> @param[out] error  propagated test failure
   subroutine test_quadrature_refines(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Radius of the probe shell (bohr)
      real(wp), parameter :: r_shell = 6.0_wp
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(mctc_error_type), allocatable :: merr
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      real(wp), allocatable :: f(:)
      real(wp) :: reference, coarse, fine, cen(3), val
      integer :: j

      call probe_structure(mol)
      cen = mol%xyz(:, 1)
      reference = shell_integral(r_shell)
      coarse = 0.0_wp
      fine = 0.0_wp

      do j = 1, 2
         if (j == 1) then
            call get_uniform_recipe(recipe, overrides, probe_nrad, probe_nang, merr, rmax=probe_rmax)
         else
            call get_uniform_recipe(recipe, overrides, 100, 590, merr, rmax=probe_rmax)
         end if
         if (.not. allocated(merr)) call new_molecular_grid(mg, merr, recipe=recipe, overrides=overrides, &
                                                            reciprocal=.false.)
         if (.not. allocated(merr)) call mg%update(mol, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         allocate (f(mg%ngrid))
         call shell_field(mg, cen, r_shell, f)
         call mg%integrate_field(f, val)
         if (j == 1) then
            coarse = abs(val - reference)/reference
         else
            fine = abs(val - reference)/reference
         end if
         deallocate (f)
         call mg%destroy()
      end do

      !> Measured 6.1e-4 -> 3.1e-6 between (50, 110) and (100, 590)
      call check(error, coarse <= 5.0e-3_wp, &
                 "molecular quadrature of a diffuse shell field is off at (50, 110)")
      if (allocated(error)) return
      call check(error, fine <= 1.0e-4_wp, &
                 "molecular quadrature of a diffuse shell field is off at (100, 590)")
      if (allocated(error)) return
      call check(error, fine <= 0.1_wp*coarse, &
                 "molecular quadrature of a diffuse shell field does not refine")
   end subroutine test_quadrature_refines

   !> Tabulate `exp(-(|r - cen| - r0)^2)` on the grid points
   !>
   !> @param[in]  mg   molecular grid
   !> @param[in]  cen  shell center (bohr)
   !> @param[in]  r0   shell radius (bohr)
   !> @param[out] f    field values, length npts
   pure subroutine shell_field(mg, cen, r0, f)
      !> Molecular grid
      type(moist_math_grid_3d_molecular_type), intent(in) :: mg
      !> Shell center
      real(wp), intent(in) :: cen(3)
      !> Shell radius
      real(wp), intent(in) :: r0
      !> Field values
      real(wp), intent(out) :: f(:)

      integer :: j

      do j = 1, mg%ngrid
         f(j) = exp(-(norm2(mg%xyz(:, j) - cen) - r0)**2)
      end do
   end subroutine shell_field

   !> `4*pi*integral_0^20 r^2 exp(-(r - r0)^2) dr` by composite Simpson
   !>
   !> @param[in]  r0   shell radius (bohr)
   pure function shell_integral(r0) result(val)
      !> Shell radius
      real(wp), intent(in) :: r0
      !> Exact volume integral of the shell field
      real(wp) :: val

      integer, parameter :: nq = 200000
      real(wp), parameter :: rtop = 20.0_wp
      real(wp) :: h, s, r, y
      integer :: j

      h = rtop/real(nq, wp)
      s = 0.0_wp
      do j = 0, nq
         r = real(j, wp)*h
         y = r*r*exp(-(r - r0)**2)
         if (j == 0 .or. j == nq) then
            s = s + y
         else if (mod(j, 2) == 1) then
            s = s + 4.0_wp*y
         else
            s = s + 2.0_wp*y
         end if
      end do
      val = 4.0_wp*pi*s*h/3.0_wp
   end function shell_integral

   !> `exp(-a |r - c|^2)`
   !>
   !> @param[in]  r    evaluation point (bohr)
   !> @param[in]  c    gaussian center (bohr)
   !> @param[in]  a    exponent (bohr^-2)
   pure function gaussian(r, c, a) result(val)
      !> Evaluation point
      real(wp), intent(in) :: r(3)
      !> Gaussian center
      real(wp), intent(in) :: c(3)
      !> Exponent
      real(wp), intent(in) :: a
      !> Function value
      real(wp) :: val

      val = exp(-a*sum((r - c)**2))
   end function gaussian

   !> Run both transform routes on one grid, report their largest disagreement
   !>
   !> - Both routes evaluate the same two sums, so they agree to the FINUFFT tolerance
   !> - Larger deviation is a convention bug
   !> - Mode or wraparound errors give O(1), a wrong phase reference grows with `|k|`
   !>
   !> @param[in,out] mg     grid with its reciprocal grid already configured
   !> @param[in]    nv      batch width to prepare and transform
   !> @param[out]   dev_f   max|F12 - F3| / max|F3| (forward)
   !> @param[out]   dev_b   max|g12 - g3| / max|g3| (backward)
   !> @param[out]   ready   whether a usable transform was produced
   !> @param[out]   error   propagated test failure
   subroutine compare_transform_routes(mg, nv, dev_f, dev_b, ready, error)
      !> Grid under test
      type(moist_math_grid_3d_molecular_type), intent(inout), target :: mg
      !> Batch width
      integer, intent(in) :: nv
      !> Forward deviation (relative to the reference maximum)
      real(wp), intent(out) :: dev_f
      !> Backward deviation (relative to the reference maximum)
      real(wp), intent(out) :: dev_b
      !> Whether a usable transform was produced
      logical, intent(out) :: ready
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_molecular_trafo_type) :: t12, t3
      type(mctc_error_type), allocatable :: merr
      real(wp), allocatable :: f(:, :), g12(:, :), g3(:, :)
      complex(wp), allocatable :: fk12(:, :), fk3(:, :), fin(:, :)
      real(wp) :: scale_f, scale_b, r(3)
      integer :: i, iv

      ready = .false.
      dev_f = 0.0_wp
      dev_b = 0.0_wp

      call new_molecular_grid_trafo(t12, mg, use_type12=.true.)
      call t12%prepare(nv, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call new_molecular_grid_trafo(t3, mg, use_type12=.false.)
      call t3%prepare(nv, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call t12%destroy()
         return
      end if
      if (.not. t12%is_type12 .or. t3%is_type12) then
         call test_failed(error, "transform-route switch did not take effect")
         call t12%destroy(); call t3%destroy()
         return
      end if

      allocate (f(mg%ngrid, nv), g12(mg%ngrid, nv), g3(mg%ngrid, nv))
      allocate (fk12(mg%npts_k, nv), fk3(mg%npts_k, nv), fin(mg%npts_k, nv))

      !> Off-center Gaussians per column, nontrivial phase against `kref = point(1)`
      do iv = 1, nv
         do i = 1, mg%ngrid
            r = mg%point(i)
            f(i, iv) = gaussian(r, [0.4_wp*iv, -0.3_wp*iv, 0.7_wp*iv], probe_alpha)
         end do
      end do

      call forward(error, t12, f, fk12)
      if (.not. allocated(error)) call forward(error, t3, f, fk3)
      if (allocated(error)) then
         call t12%destroy(); call t3%destroy()
         return
      end if
      ! Fail on a zero or NaN type-3 reference, else dev_f stays at its initial value
      scale_f = maxval(abs(fk3))
      if (.not. ieee_is_finite(scale_f) .or. scale_f <= 0.0_wp) then
         call test_failed(error, "molecular NUFFT: type-3 forward reference has a non-finite or zero norm")
         call t12%destroy(); call t3%destroy()
         return
      end if
      dev_f = maxval(abs(fk12 - fk3))/scale_f

      !> Bit-identical input in both directions, so only the transform differs
      fin = fk3
      call backward(error, t12, fin, g12)
      fin = fk3
      if (.not. allocated(error)) call backward(error, t3, fin, g3)
      if (allocated(error)) then
         call t12%destroy(); call t3%destroy()
         return
      end if
      scale_b = maxval(abs(g3))
      if (.not. ieee_is_finite(scale_b) .or. scale_b <= 0.0_wp) then
         call test_failed(error, "molecular NUFFT: type-3 backward reference has a non-finite or zero norm")
         call t12%destroy(); call t3%destroy()
         return
      end if
      dev_b = maxval(abs(g12 - g3))/scale_b

      call t12%destroy()
      call t3%destroy()
      ready = .true.
   end subroutine compare_transform_routes

   !> Check type-1/2 against type-3 elementwise on the molecular probe grid
   !>
   !> @param[out] error  propagated test failure
   subroutine test_type12_vs_type3_molecular(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(mctc_error_type), allocatable :: merr
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      real(wp) :: dev_f, dev_b
      logical :: ready

      call probe_structure(mol)
      call get_uniform_recipe(recipe, overrides, probe_nrad, probe_nang, merr, rmax=probe_rmax)
      if (.not. allocated(merr)) call new_molecular_grid(mg, merr, recipe=recipe, overrides=overrides, &
                                                         reciprocal=.false.)
      if (.not. allocated(merr)) call mg%update(mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call molecular_grid_set_kgrid(mg, probe_dr, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call mg%destroy()
         return
      end if

      call compare_transform_routes(mg, 2, dev_f, dev_b, ready, error)
      if (allocated(error) .or. .not. ready) then
         call mg%destroy()
         return
      end if
      write (output_unit, "(a,es10.3,a,es10.3)") &
         "    molecular grid: type-1/2 vs type-3 max dev  forward ", dev_f, &
         "   backward ", dev_b
      call check(error, dev_f <= route_agreement_tol, &
                 "type-1 forward disagrees with type-3 beyond the NUFFT tolerance")
      if (.not. allocated(error)) then
         call check(error, dev_b <= route_agreement_tol, &
                    "type-2 backward disagrees with type-3 beyond the NUFFT tolerance")
      end if
      call mg%destroy()
   end subroutine test_type12_vs_type3_molecular

   !> Check type-1/2 against type-3 on a near-uniform point cloud
   !>
   !> - Opposite extreme to the clustered molecular cloud, every fine-grid cell equally loaded
   !>
   !> @param[out] error  propagated test failure
   subroutine test_type12_vs_type3_uniform(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Points per axis of the near-uniform probe lattice
      integer, parameter :: nlat = 12
      !> Lattice spacing (bohr)
      real(wp), parameter :: alat = 0.8_wp

      type(moist_math_grid_3d_molecular_type), target :: mg
      type(mctc_error_type), allocatable :: merr
      real(wp) :: dev_f, dev_b, o
      logical :: ready
      integer :: ix, iy, iz, i

      mg%ngrid = nlat**3
      allocate (mg%xyz(3, mg%ngrid), mg%w(mg%ngrid))
      allocate (mg%atom_offset(2), mg%nrad_per_atom(1), mg%nang_per_atom(1))
      mg%atom_offset = [1, mg%ngrid + 1]
      mg%nrad_per_atom = nlat
      mg%nang_per_atom = nlat*nlat
      i = 0
      do iz = 1, nlat
         do iy = 1, nlat
            do ix = 1, nlat
               i = i + 1
               !> Deterministic sub-spacing jitter
               o = 0.05_wp*alat*sin(real(7*i, wp))
               mg%xyz(:, i) = ([real(ix, wp), real(iy, wp), real(iz, wp)] &
                               - 0.5_wp*real(nlat + 1, wp))*alat + o
               mg%w(i) = alat**3
            end do
         end do
      end do

      call molecular_grid_set_kgrid(mg, alat, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call mg%destroy()
         return
      end if

      call compare_transform_routes(mg, 1, dev_f, dev_b, ready, error)
      if (allocated(error) .or. .not. ready) then
         call mg%destroy()
         return
      end if
      write (output_unit, "(a,es10.3,a,es10.3)") &
         "    uniform lattice: type-1/2 vs type-3 max dev forward ", dev_f, &
         "   backward ", dev_b
      call check(error, dev_f <= route_agreement_tol, &
                 "type-1 forward disagrees with type-3 on a near-uniform cloud")
      if (.not. allocated(error)) then
         call check(error, dev_b <= route_agreement_tol, &
                    "type-2 backward disagrees with type-3 on a near-uniform cloud")
      end if
      call mg%destroy()
   end subroutine test_type12_vs_type3_uniform

   !> Check that a stale trafo refuses to transform after a grid update
   !>
   !> - Holds for a pure translation (`ngrid` unchanged), only the generation guard catches it
   !>
   !> @param[out] error  propagated test failure
   subroutine test_trafo_stale_after_update(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      type(mctc_error_type), allocatable :: merr
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)
      logical :: ready

      call setup_probe(mol, mg, trafo, ready, error)
      if (allocated(error) .or. .not. ready) then
         call mg%destroy()
         return
      end if

      mol%xyz(1, 1) = mol%xyz(1, 1) + 0.05_wp
      call mg%update(mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call trafo%destroy()
         call mg%destroy()
         return
      end if

      allocate (f(mg%ngrid, 1), fk(mg%npts_k, 1))
      f = 0.0_wp
      call trafo%fft_r2k(f, fk, merr)
      call check(error, allocated(merr), "a trafo must refuse to transform after its grid was updated")
      if (.not. allocated(error)) then
         call check(error, index(merr%message, "geometry changed") > 0, &
                    "the stale-trafo error must name the geometry-generation mismatch")
      end if
      call trafo%destroy()
      call mg%destroy()
   end subroutine test_trafo_stale_after_update

   !> Check that a stale trafo refuses to transform after its grid was constructed again
   !>
   !> - Same recipe, same sequence of update and `molecular_grid_set_kgrid`, one
   !>   atom moved: the point count is unchanged, so only a generation
   !>   counter that survives `new_molecular_grid` catches the stale plans
   !>
   !> @param[out] error  propagated test failure
   subroutine test_trafo_stale_after_reconstruction(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      type(mctc_error_type), allocatable :: merr
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)
      logical :: ready

      call setup_probe(mol, mg, trafo, ready, error)
      if (allocated(error) .or. .not. ready) then
         call mg%destroy()
         return
      end if

      mol%xyz(1, 1) = mol%xyz(1, 1) + 0.05_wp
      call get_uniform_recipe(recipe, overrides, probe_nrad, probe_nang, merr, rmax=probe_rmax)
      if (.not. allocated(merr)) call new_molecular_grid(mg, merr, recipe=recipe, overrides=overrides, &
                                                         reciprocal=.false.)
      if (.not. allocated(merr)) call mg%update(mol, merr)
      if (.not. allocated(merr)) call molecular_grid_set_kgrid(mg, probe_dr, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call trafo%destroy()
         call mg%destroy()
         return
      end if

      allocate (f(mg%ngrid, 1), fk(mg%npts_k, 1))
      f = 0.0_wp
      call trafo%fft_r2k(f, fk, merr)
      call check(error, allocated(merr), &
                 "a trafo must refuse to transform after its grid was constructed again")
      if (.not. allocated(error)) then
         call check(error, index(merr%message, "geometry changed") > 0, &
                    "the stale-trafo error must name the geometry-generation mismatch")
      end if
      call trafo%destroy()
      call mg%destroy()
   end subroutine test_trafo_stale_after_reconstruction

   !> Check destroy of a stale trafo and fresh prepare on the updated grid
   !>
   !> @param[out] error  propagated test failure
   subroutine test_trafo_destroy_recreate(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      type(mctc_error_type), allocatable :: merr
      type(moist_math_grid_3d_cartesian_type) :: cg
      real(wp), allocatable :: f(:, :), g(:, :)
      complex(wp), allocatable :: fk(:, :)
      real(wp) :: cen(3), emax, dev
      logical :: ready
      integer :: j

      call setup_probe(mol, mg, trafo, ready, error)
      if (allocated(error) .or. .not. ready) then
         call mg%destroy()
         return
      end if

      mol%xyz(1, 1) = mol%xyz(1, 1) + 0.05_wp
      call mg%update(mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call trafo%destroy()
         call mg%destroy()
         return
      end if

      call trafo%destroy()
      call new_molecular_grid_trafo(trafo, mg)
      call trafo%prepare(1, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call mg%destroy()
         return
      end if

      cen = mol%xyz(:, 1)
      allocate (f(mg%ngrid, 1), g(mg%ngrid, 1), fk(mg%npts_k, 1))
      do j = 1, mg%ngrid
         f(j, 1) = gaussian(mg%xyz(:, j), cen, probe_alpha)
      end do
      call forward(error, trafo, f, fk)
      if (.not. allocated(error)) then
         ! Same bound as test_forward_analytic; 3.1e-4 measured, 4.3e-2 with the pre-update center
         call max_ft_error(fk(:, 1), mg%npts_k, mg%kref, cen, mg, .true., cg, emax)
         call check(error, emax <= 1.0e-3_wp, &
                    "recreated molecular NUFFT disagrees with the analytic Gaussian FT")
      end if
      if (.not. allocated(error)) call backward(error, trafo, fk, g)
      if (.not. allocated(error)) then
         ! Band-limited to ~6.5e-3 at k_max = pi/dr
         ! max |g - f| = 4.8e-3 measured, ~6x headroom
         emax = 0.0_wp
         do j = 1, mg%ngrid
            dev = abs(g(j, 1) - f(j, 1))
            if (.not. ieee_is_finite(dev)) then
               emax = huge(1.0_wp)
               exit
            end if
            emax = max(emax, dev)
         end do
         call check(error, emax <= 3.0e-2_wp, &
                    "recreated molecular NUFFT round trip does not reproduce the Gaussian")
      end if
      call trafo%destroy()
      call mg%destroy()
   end subroutine test_trafo_destroy_recreate

   ! No "failed prepare" test: FINUFFT aborts on bad input instead of returning `ier /= 0`

end module test_math_grid_nufft
