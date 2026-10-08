!> Transform-level tests for the atom-centered molecular 3D grid (FINUFFT NUFFT)
!>
!> Analytic checks of the pieces a 3D-RISM UV k-space solve is built from
!>
!> - Round-trip diagonal, forward vs analytic FT, k-space convolution
!> - Type-1/2 vs type-3 routes
!> - Cross-validation: molecular NUFFT vs Cartesian FFT of an LJ + Yukawa field
!> - Tolerances assume a nucleus-localised field, the cloud resolves `2*pi/k` only near the nuclei
module test_math_grid_3d_nufft
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_env, only: wp, mctc_error_type => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use mstore, only: get_structure
   use moist_math_grid_3d_base, only: moist_math_grid_3d_trafo_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type, &
      & moist_math_grid_atomic_recipe_override_type
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, &
      & new_molecular_point_grid, molecular_grid_set_kgrid, &
      & moist_math_grid_3d_molecular_trafo_type, new_molecular_grid_trafo, &
      & default_nufft_tol
   use test_helpers, only: get_uniform_recipe, get_cartesian_gaussian_grid, center_at_origin

   implicit none(type, external)
   private

   public :: collect_math_grid_3d_nufft

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
   subroutine collect_math_grid_3d_nufft(testsuite)
      !> Registered unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("roundtrip_diagonal", test_roundtrip_diagonal), &
                  new_unittest("forward_matches_analytic_ft", test_forward_analytic), &
                  new_unittest("convolution_matches_analytic", test_convolution), &
                  new_unittest("type12_matches_type3", test_type12_vs_type3_molecular), &
                  new_unittest("ft_lj_coulomb_nufft_vs_fft", test_ft_lj_coulomb) &
                  ]
   end subroutine collect_math_grid_3d_nufft

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
      if (.not. allocated(merr)) call new_molecular_point_grid(mg, merr, recipe=recipe, overrides=overrides, &
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
      ! The analytic probes cover the production route
      if (.not. trafo%uses_type12()) then
         call test_failed(error, "default molecular transform must use NUFFT type 1/2")
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
   !> - Pins measure, sign, phase reference and k ordering of the NUFFT
   subroutine test_forward_analytic(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
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
      call max_ft_error(fk(:, 1), mg%npts_k, mg%get_kref(), cen, mg, emax)
      ! Quadrature-limited: 3.0e-4 at (50, 110), 2.7e-7 from Lebedev order 302
      call check(error, emax <= 1.0e-3_wp, &
                 "molecular forward NUFFT disagrees with the analytic Gaussian FT")
      call trafo%destroy()
      call mg%destroy()
   end subroutine test_forward_analytic

   !> Find the largest deviation `|F_num(k) - F_exact(k)| / F_exact(0)` in the band `|k| <= pi/dr`
   !>
   !> - Only the inscribed ball of the cube grid is resolved
   !> - Full cube gives 3.7e-3 instead of 3.0e-4
   !>
   !> @param[in]  fk    numerical transform, length nk
   !> @param[in]  nk    number of reciprocal points
   !> @param[in]  kref  phase reference of the transform (bohr)
   !> @param[in]  cen   gaussian center (bohr)
   !> @param[in]  mg    molecular grid supplying kpoint
   !> @param[out] emax  largest relative deviation
   subroutine max_ft_error(fk, nk, kref, cen, mg, emax)
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
      !> Largest relative deviation
      real(wp), intent(out) :: emax

      real(wp) :: kvec(3), f0, phase, kband, dev
      complex(wp) :: fex
      integer :: j

      f0 = (pi/probe_alpha)**1.5_wp
      kband = pi/probe_dr
      emax = 0.0_wp
      do j = 1, nk
         kvec = mg%kpoint(j)
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

   !> Run both transform routes on one grid and compare them elementwise
   !>
   !> - Both routes evaluate the same two sums, so they agree to the FINUFFT tolerance
   !> - Larger deviation is a convention bug
   !> - Mode or wraparound errors give O(1), a wrong phase reference grows with `|k|`
   !>
   !> @param[in,out] mg             grid with its reciprocal grid already configured
   !> @param[in]     nv             batch width to prepare and transform
   !> @param[in]     forward_what   failure context of the forward comparison
   !> @param[in]     backward_what  failure context of the backward comparison
   !> @param[out]    error          propagated test failure
   subroutine compare_transform_routes(mg, nv, forward_what, backward_what, error)
      !> Grid under test
      type(moist_math_grid_3d_molecular_type), intent(inout), target :: mg
      !> Batch width
      integer, intent(in) :: nv
      !> Failure context of the forward comparison
      character(len=*), intent(in) :: forward_what
      !> Failure context of the backward comparison
      character(len=*), intent(in) :: backward_what
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_molecular_trafo_type) :: t12, t3
      type(mctc_error_type), allocatable :: merr
      real(wp), allocatable :: f(:, :), g12(:, :), g3(:, :)
      complex(wp), allocatable :: fk12(:, :), fk3(:, :), fin(:, :)
      real(wp) :: scale_f, scale_b, r(3)
      integer :: i, iv

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
         return
      end if
      ! Two engines on one backend would pass by comparing it against itself
      if (.not. t12%uses_type12() .or. t3%uses_type12()) then
         call test_failed(error, "transform-route switch did not take effect")
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
      if (allocated(error)) return
      ! Fail on a zero or NaN type-3 reference, else the scaled bound is meaningless
      scale_f = maxval(abs(fk3))
      if (.not. ieee_is_finite(scale_f) .or. scale_f <= 0.0_wp) then
         call test_failed(error, "molecular NUFFT: type-3 forward reference has a non-finite or zero norm")
         return
      end if
      forward_check: do iv = 1, nv
         do i = 1, mg%npts_k
            call check(error, ieee_is_finite(real(fk3(i, iv))) .and. ieee_is_finite(aimag(fk3(i, iv))), &
               & "type-3 forward reference is not finite")
            if (allocated(error)) exit forward_check
            call check(error, fk12(i, iv), fk3(i, iv), thr=route_agreement_tol*scale_f, more=forward_what)
            if (allocated(error)) exit forward_check
         end do
      end do forward_check
      if (allocated(error)) return

      !> Bit-identical input in both directions, so only the transform differs
      fin = fk3
      call backward(error, t12, fin, g12)
      fin = fk3
      if (.not. allocated(error)) call backward(error, t3, fin, g3)
      if (allocated(error)) return
      scale_b = maxval(abs(g3))
      if (.not. ieee_is_finite(scale_b) .or. scale_b <= 0.0_wp) then
         call test_failed(error, "molecular NUFFT: type-3 backward reference has a non-finite or zero norm")
         return
      end if
      backward_check: do iv = 1, nv
         do i = 1, mg%ngrid
            call check(error, ieee_is_finite(g3(i, iv)), "type-3 backward reference is not finite")
            if (allocated(error)) exit backward_check
            call check(error, g12(i, iv), g3(i, iv), thr=route_agreement_tol*scale_b, more=backward_what)
            if (allocated(error)) exit backward_check
         end do
      end do backward_check
      if (allocated(error)) return

      call t12%destroy()
      call t3%destroy()
   end subroutine compare_transform_routes

   !> Check type-1/2 against type-3 elementwise on the molecular probe grid
   subroutine test_type12_vs_type3_molecular(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(mctc_error_type), allocatable :: merr
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)

      call probe_structure(mol)
      call get_uniform_recipe(recipe, overrides, probe_nrad, probe_nang, merr, rmax=probe_rmax)
      if (.not. allocated(merr)) call new_molecular_point_grid(mg, merr, recipe=recipe, overrides=overrides, &
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

      call compare_transform_routes(mg, 2, "type-1 forward disagrees with type-3 beyond the NUFFT tolerance", &
         & "type-2 backward disagrees with type-3 beyond the NUFFT tolerance", error)
      call mg%destroy()
   end subroutine test_type12_vs_type3_molecular

   !> Lennard-Jones + screened-Coulomb (Yukawa) field at `r`, summed over sources
   !>
   !> - Sources at `coord(:,a)` with per-atom well depth `eps`, size `sig`, charge `q`
   !> - Source distance softened to `dsoft`, so the field is bounded at the grid
   !>   points nearest a core
   !> - Coulomb part screened by `kappa`, so the field is compact
   !> - Together they keep the transform well resolved on both the uniform
   !>   Cartesian and the atom-centered grid, which makes a cross-grid comparison
   !>   meaningful
   !>
   !> @param[in] r      field point in bohr
   !> @param[in] coord  source coordinates, shape (3, nat)
   !> @param[in] eps    per-source LJ well depth
   !> @param[in] sig    per-source LJ size
   !> @param[in] q      per-source charge
   !> @param[in] dsoft  distance softening in bohr
   !> @param[in] kappa  Coulomb screening in 1/bohr
   pure function lj_screened_coulomb(r, coord, eps, sig, q, dsoft, kappa) result(val)
      !> Field point (bohr)
      real(wp), intent(in) :: r(3)
      !> Source coordinates, shape (3, nat)
      real(wp), intent(in) :: coord(:, :)
      !> Per-source LJ well depth / size and charge
      real(wp), intent(in) :: eps(:), sig(:), q(:)
      !> Distance softening and Coulomb screening (bohr, 1/bohr)
      real(wp), intent(in) :: dsoft, kappa
      !> Field value at r
      real(wp) :: val

      integer :: a
      real(wp) :: d, sr6

      val = 0.0_wp
      do a = 1, size(eps)
         d = sqrt((r(1) - coord(1, a))**2 + (r(2) - coord(2, a))**2 &
                  + (r(3) - coord(3, a))**2)
         d = max(d, dsoft)
         sr6 = (sig(a)/d)**6
         val = val + 4.0_wp*eps(a)*(sr6*sr6 - sr6) + q(a)*exp(-kappa*d)/d
      end do
   end function lj_screened_coulomb

   !> Finiteness of a complex value, true iff both its parts are finite
   !>
   !> - Guards the `max(maxerr, abs(...))` accumulators below
   !> - A NaN sample would otherwise vanish silently, an IEEE comparison against
   !>   NaN is always false so gfortran's `max` can keep the running finite
   !>   accumulator
   !>
   !> @param[in] z  value to test
   pure function complex_is_finite(z) result(ok)
      !> Value to test
      complex(wp), intent(in) :: z
      !> True if both the real and imaginary parts are finite
      logical :: ok

      ok = ieee_is_finite(real(z, wp)) .and. ieee_is_finite(aimag(z))
   end function complex_is_finite

   !> Cross-validate the molecular-grid NUFFT against the Cartesian-grid FFT
   !>
   !> Common Lennard-Jones + screened-Coulomb field sourced at the atoms
   !> - Molecular NUFFT forward transform reproduces the brute-force direct DFT
   !>   `sum_j w_j f(r_j) exp(-i k.(r_j-kref))` at sampled k-points (tight,
   !>   FINUFFT correctness check independent of grid accuracy)
   !> - DC mode (k=0) agrees between the two grids, both equal the grid
   !>   quadrature of the field, `integral f dV`
   !> - At the lowest k along an axis the Cartesian FFT magnitude matches the
   !>   molecular grid's transform magnitude (same Fourier coefficient from two
   !>   independent grids and transforms, to grid accuracy)
   subroutine test_ft_lj_coulomb(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(mctc_error_type), allocatable :: merr
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: mtrafo
      type(moist_math_grid_3d_cartesian_type), target :: cgrid
      class(moist_math_grid_3d_trafo_type), allocatable :: ctrafo
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)

      !> Softened LJ + screened Coulomb model parameters (see field helper)
      real(wp), parameter :: dsoft = 3.5_wp, kappa = 0.5_wp
      integer, parameter :: nsamp = 8
      real(wp), allocatable :: eps(:), sig(:), q(:), coord(:, :)
      real(wp), allocatable :: f_mol(:, :), f_cart(:, :)
      complex(wp), allocatable :: fk_mol(:, :), fk_cart(:, :)
      real(wp) :: kvec(3), ph, fscale, maxerr, dc_mol, dc_cart, dcref
      complex(wp) :: acc
      integer :: nat, a, j, s, jk, k0_mol, mx

      ! --- Solute structure + synthetic LJ/charge parameters ---
      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      nat = mol%nat
      allocate (coord(3, nat), eps(nat), sig(nat), q(nat))
      coord = mol%xyz
      do a = 1, nat
         eps(a) = 0.01_wp
         sig(a) = 5.0_wp
         q(a) = real(1 - 2*mod(a, 2), wp)*0.5_wp   ! alternating +/- 0.5 charges
      end do

      ! --- Molecular grid + NUFFT of the field ---
      call get_uniform_recipe(recipe, overrides, 40, 110, merr, rmax=8.0_wp)
      if (.not. allocated(merr)) call new_molecular_point_grid(mgrid, merr, recipe=recipe, &
         & overrides=overrides, reciprocal=.false.)
      if (.not. allocated(merr)) call mgrid%update(mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call molecular_grid_set_kgrid(mgrid, 0.5_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call new_molecular_grid_trafo(mtrafo, mgrid)
      call mtrafo%prepare(1, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if

      allocate (f_mol(mgrid%ngrid, 1), fk_mol(mgrid%npts_k, 1))
      do j = 1, mgrid%ngrid
         f_mol(j, 1) = lj_screened_coulomb(mgrid%xyz(:, j), coord, eps, sig, q, dsoft, kappa)
      end do
      call forward(error, mtrafo, f_mol, fk_mol)
      if (allocated(error)) return
      fscale = maxval(abs(fk_mol(:, 1)))

      ! --- Cartesian grid (centered, enclosing the field) + FFT of same field ---
      call get_cartesian_gaussian_grid(cgrid, 64, 64, 64, 0.3_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call cgrid%new_trafo(ctrafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      allocate (f_cart(cgrid%ngrid, 1), fk_cart(cgrid%npts_k, 1))
      do j = 1, cgrid%ngrid
         f_cart(j, 1) = lj_screened_coulomb(cgrid%point(j), coord, eps, sig, q, dsoft, kappa)
      end do
      call forward(error, ctrafo, f_cart, fk_cart)
      if (allocated(error)) return

      ! --- Check 1: molecular NUFFT == brute-force direct DFT (tight) ---
      maxerr = 0.0_wp
      do s = 1, nsamp
         jk = 1 + (s - 1)*(mgrid%npts_k/nsamp)
         kvec = mgrid%kpoint(jk)
         acc = (0.0_wp, 0.0_wp)
         do j = 1, mgrid%ngrid
            ph = -dot_product(kvec, mgrid%xyz(:, j) - mgrid%get_kref())
            acc = acc + mgrid%w(j)*f_mol(j, 1)*cmplx(cos(ph), sin(ph), wp)
         end do
         if (.not. complex_is_finite(fk_mol(jk, 1)) .or. .not. complex_is_finite(acc)) then
            call test_failed(error, "molecular NUFFT check produced a non-finite value")
            return
         end if
         maxerr = max(maxerr, abs(fk_mol(jk, 1) - acc))
      end do
      call check(error, maxerr <= 1.0e-6_wp*fscale, &
         & "molecular NUFFT deviates from the direct DFT of the LJ+Coulomb field")

      ! --- Check 2: DC mode (k=0) agrees between the two grids ---
      if (.not. allocated(error)) then
         k0_mol = 0
         do j = 1, mgrid%npts_k
            if (maxval(abs(mgrid%kpoint(j))) < 1.0e-12_wp) then
               k0_mol = j; exit
            end if
         end do
         call check(error, k0_mol > 0, "molecular k-grid has no DC mode")
         if (.not. allocated(error)) then
            dc_mol = real(fk_mol(k0_mol, 1), wp)
            dc_cart = real(fk_cart(1, 1), wp)          ! Cartesian kpoint(1) = (0,0,0)
            call check(error, ieee_is_finite(dc_mol) .and. ieee_is_finite(dc_cart), &
               & "DC-mode comparison produced a non-finite value")
            if (.not. allocated(error)) then
               dcref = max(abs(dc_cart), 1.0e-30_wp)
               call check(error, abs(dc_mol - dc_cart) <= 2.0e-2_wp*dcref, &
                  & "molecular and Cartesian grids disagree on the field integral (DC)")
            end if
         end if
      end if

      ! --- Check 3: low-k Fourier magnitudes agree across the two grids ---
      ! Cartesian FFT value at q=(mx*dkx,0,0) (R2C layout, kx fastest) vs the
      ! molecular grid's transform at the same physical q
      ! Magnitudes are phase-reference independent; with check 1 this ties
      ! the NUFFT and FFT spectra of the same field together at low k
      if (.not. allocated(error)) then
         maxerr = 0.0_wp
         do mx = 1, 3
            jk = mx + 1
            kvec = cgrid%kpoint(jk)               ! = (mx*dkx, 0, 0)
            acc = (0.0_wp, 0.0_wp)
            do j = 1, mgrid%ngrid
               ph = -dot_product(kvec, mgrid%xyz(:, j))
               acc = acc + mgrid%w(j)*f_mol(j, 1)*cmplx(cos(ph), sin(ph), wp)
            end do
            if (.not. complex_is_finite(fk_cart(jk, 1)) .or. .not. complex_is_finite(acc)) then
               call test_failed(error, "low-k cross-grid check produced a non-finite value")
               return
            end if
            maxerr = max(maxerr, abs(abs(fk_cart(jk, 1)) - abs(acc)))
         end do
         call check(error, maxerr <= 3.0e-2_wp*abs(fk_cart(1, 1)), &
            & "Cartesian FFT and molecular transform disagree at low k")
      end if

      call mtrafo%destroy()
      call mgrid%destroy()
      call ctrafo%destroy()
      call cgrid%destroy()
   end subroutine test_ft_lj_coulomb

end module test_math_grid_3d_nufft
