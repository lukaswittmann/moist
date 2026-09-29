!> Test suite for the moist math/grid submodule
!>
!> The central goal is to verify that every concrete grid type integrates
!> analytically known functions to its expected accuracy; all four grid
!> types are driven through the abstract bases ([[moist_math_grid_1d_base]],
!> [[moist_math_grid_3d_base]]) via shared checker routines, so the same
!> assertions exercise polymorphic dispatch of `integrate` / `integrate_field`
!> / `measure` / `point` for:
!>
!>   1D radial:  uniform (equidistant) and Chebyshev-2 grids
!>   3D volume:  Cartesian (uniform box) and atom-centered molecular grids
!>
!> Reference integrals (all over R^3, lengths in bohr):
!>   integral exp(-r^2)            dV = pi^(3/2)
!>   integral x^2 exp(-r^2)        dV = pi^(3/2)/2
!>
!> The radial grids integrate the spherically symmetric Gaussian through the
!> 4*pi*r^2 dr volume measure; the 3D grids sample the full Cartesian field
!>
!> A handful of structural tests (Chebyshev FBT adjoint identity, idempotent
!> destroy, molecular-grid pruning) round out the suite
module test_math_grid
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: center_at_origin
   use moist_math_grid, only: moist_math_grid_3d_molecular_type, moist_math_grid_3d_type, &
      & new_molecular_grid, new_molecular_grid_uniform, new_molecular_grid_uniform_qc_handymod, &
      & molecular_grid_set_kgrid, moist_math_grid_3d_molecular_trafo_type, new_molecular_grid_trafo
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, &
      & new_cartesian_grid_3d
   use moist_math_grid_1d_base, only: moist_math_grid_1d_type, &
      & moist_math_grid_1d_trafo_type
   use moist_math_grid_1d_uniform, only: moist_math_grid_1d_uniform_type, &
      & new_uniform_radial_grid
   use moist_math_grid_1d_chebyshev, only: moist_math_grid_1d_chebyshev_type, &
      & new_chebyshev_radial_grid
   implicit none(type, external)
   private

   public :: collect_math_grid

contains

   !> Collect all math_grid tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_grid(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("grid_uniform_radial_gaussian", test_grid_uniform_radial), &
         new_unittest("grid_chebyshev_radial_gaussian", test_grid_chebyshev_radial), &
         new_unittest("grid_cartesian_3d_gaussian", test_grid_cartesian_3d), &
         new_unittest("grid_cartesian_3d_anisotropic", test_grid_cartesian_aniso), &
         new_unittest("grid_cartesian_3d_modulated", test_grid_cartesian_modulated), &
         new_unittest("grid_molecular_3d_gaussian", test_grid_molecular_3d), &
         new_unittest("grid_molecular_uniform_gaussian", test_grid_molecular_uniform), &
         new_unittest("grid_molecular_modulated_gaussian", test_grid_molecular_modulated), &
         new_unittest("grid_molecular_nufft_dc_mode", test_grid_molecular_nufft), &
         new_unittest("chebyshev_radial_fbt_adjoint", test_cheb_radial_fbt_adjoint), &
         new_unittest("chebyshev_radial_fbt_batched", test_cheb_radial_fbt_batched), &
         new_unittest("chebyshev_radial_destroy_idempotent", test_cheb_radial_destroy), &
         new_unittest("molecular_grid_integrate_constant", test_mol_grid_integrate_const), &
         new_unittest("molecular_grid_qc_handymod_smoke", test_mol_grid_qc_handymod_smoke), &
         new_unittest("molecular_grid_uniform_regression", test_mol_grid_uniform_regression), &
         new_unittest("molecular_grid_destroy_idempotent", test_mol_grid_destroy), &
         new_unittest("molecular_grid_pruning_threshold", test_mol_grid_pruning) &
         ]
   end subroutine collect_math_grid

   ! --------------------------------------------------------------------------
   ! Analytic integrands (module-scope so they match the procedure interfaces)
   ! --------------------------------------------------------------------------

   !> Spherically symmetric unit Gaussian f(r) = exp(-r^2)
   !>
   !> @param[in] r  Radius (bohr)
   pure function gaussian_unit_1d(r) result(val)
      !> Radius (bohr)
      real(wp), intent(in) :: r
      !> exp(-r^2)
      real(wp) :: val
      val = exp(-r*r)
   end function gaussian_unit_1d

   !> 3D unit Gaussian f(r) = exp(-(x^2+y^2+z^2))
   !>
   !> @param[in] r  Point (bohr)
   pure function gaussian_unit_3d(r) result(val)
      !> Point (bohr)
      real(wp), intent(in) :: r(3)
      !> exp(-|r|^2)
      real(wp) :: val
      val = exp(-(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))
   end function gaussian_unit_3d

   !> Anisotropic integrand f(r) = x^2 * exp(-|r|^2)
   !>
   !> @param[in] r  Point (bohr)
   pure function x2_gaussian_3d(r) result(val)
      !> Point (bohr)
      real(wp), intent(in) :: r(3)
      !> x^2 exp(-|r|^2)
      real(wp) :: val
      val = r(1)*r(1)*exp(-(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))
   end function x2_gaussian_3d

   !> Higher-frequency, origin-centered integrand: a unit Gaussian modulated by cos(2x)
   !>
   !> Over R^3 (each Cartesian factor separable):
   !>   integral exp(-|r|^2) cos(2x) dV = pi^(3/2) * exp(-1)
   !> The cos(2x) factor adds ~2 oscillations across the Gaussian's support,
   !> probing how finely the grid resolves an oscillatory field that the
   !> smooth Gaussian alone cannot reveal
   !>
   !> @param[in] r  Point (bohr)
   pure function gaussian_cos_3d(r) result(val)
      !> Point (bohr)
      real(wp), intent(in) :: r(3)
      !> exp(-|r|^2) cos(2x)
      real(wp) :: val
      val = exp(-(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))*cos(2.0_wp*r(1))
   end function gaussian_cos_3d

   !> Trivial integrand f == 1 used by the constant-integration test
   !>
   !> @param[in] r  Point (bohr)
   pure function one_3d(r) result(val)
      !> Point (bohr)
      real(wp), intent(in) :: r(3)
      !> 1
      real(wp) :: val
      val = 1.0_wp + 0.0_wp*r(1)   ! reference r to silence unused warning
   end function one_3d

   ! --------------------------------------------------------------------------
   ! Shared, fully type-agnostic checkers (operate on the abstract bases)
   ! --------------------------------------------------------------------------

   !> Integrate exp(-r^2) over R^3 through the abstract 1D radial grid
   !>
   !> Compares to 4*pi*integral_0^inf r^2 exp(-r^2) dr = pi^(3/2); both the
   !> analytic-sampling `integrate` and the tabulated-field `integrate_field`
   !> paths must hit the reference to the same tolerance
   !>
   !> @param[out] error  Test error
   !> @param[in]  grid   Any concrete 1D radial grid, via its abstract base
   !> @param[in]  thr    Absolute tolerance on the quadrature
   !> @param[in]  tag    Short grid label used in failure messages
   subroutine check_radial_gaussian(error, grid, thr, tag)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Grid under test (abstract base)
      class(moist_math_grid_1d_type), intent(in) :: grid
      !> Absolute tolerance
      real(wp), intent(in) :: thr
      !> Grid label
      character(*), intent(in) :: tag

      real(wp) :: res, field_res
      real(wp), allocatable :: f(:)
      integer :: i

      call grid%integrate(gaussian_unit_1d, res)
      call check(error, res, pi*sqrt(pi), thr=thr, &
         & more=tag//": radial Gaussian integral deviates from pi^(3/2)")
      if (allocated(error)) return

      allocate (f(grid%npts))
      do i = 1, grid%npts
         f(i) = gaussian_unit_1d(grid%r(i))
      end do
      call grid%integrate_field(f, field_res)
      call check(error, field_res, pi*sqrt(pi), thr=thr, &
         & more=tag//": integrate_field of the Gaussian deviates from pi^(3/2)")
   end subroutine check_radial_gaussian

   !> Integrate exp(-|r|^2) over R^3 through the abstract 3D grid
   !>
   !> Compares to pi^(3/2); both the analytic-sampling `integrate` and the
   !> tabulated-field `integrate_field` paths must hit the reference to the
   !> same tolerance
   !>
   !> @param[out] error  Test error
   !> @param[in]  grid   Any concrete 3D grid, via its abstract base
   !> @param[in]  thr    Absolute tolerance on the quadrature
   !> @param[in]  tag    Short grid label used in failure messages
   subroutine check_3d_gaussian(error, grid, thr, tag)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Grid under test (abstract base)
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Absolute tolerance
      real(wp), intent(in) :: thr
      !> Grid label
      character(*), intent(in) :: tag

      real(wp) :: res, field_res
      real(wp), allocatable :: f(:)
      integer :: i

      call grid%integrate(gaussian_unit_3d, res)
      call check(error, res, pi*sqrt(pi), thr=thr, &
         & more=tag//": 3D Gaussian integral deviates from pi^(3/2)")
      if (allocated(error)) return

      allocate (f(grid%ngrid))
      do i = 1, grid%ngrid
         f(i) = gaussian_unit_3d(grid%point(i))
      end do
      call grid%integrate_field(f, field_res)
      call check(error, field_res, pi*sqrt(pi), thr=thr, &
         & more=tag//": integrate_field of the Gaussian deviates from pi^(3/2)")
   end subroutine check_3d_gaussian

   !> Integrate the cos(2x)-modulated unit Gaussian over R^3, abstract 3D grid
   !>
   !> Compares to pi^(3/2)*exp(-1); same two-path contract as
   !> [[check_3d_gaussian]], but the integrand carries more spatial structure,
   !> stressing how finely the grid resolves an oscillatory, origin-centered
   !> field rather than a single smooth bump
   !>
   !> @param[out] error  Test error
   !> @param[in]  grid   Any concrete 3D grid, via its abstract base
   !> @param[in]  thr    Absolute tolerance on the quadrature
   !> @param[in]  tag    Short grid label used in failure messages
   subroutine check_3d_modulated(error, grid, thr, tag)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Grid under test (abstract base)
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Absolute tolerance
      real(wp), intent(in) :: thr
      !> Grid label
      character(*), intent(in) :: tag

      real(wp) :: res, field_res, ref
      real(wp), allocatable :: f(:)
      integer :: i

      ref = pi*sqrt(pi)*exp(-1.0_wp)

      call grid%integrate(gaussian_cos_3d, res)
      call check(error, res, ref, thr=thr, &
         & more=tag//": modulated Gaussian integral deviates from pi^(3/2)/e")
      if (allocated(error)) return

      allocate (f(grid%ngrid))
      do i = 1, grid%ngrid
         f(i) = gaussian_cos_3d(grid%point(i))
      end do
      call grid%integrate_field(f, field_res)
      call check(error, field_res, ref, thr=thr, &
         & more=tag//": integrate_field of the modulated Gaussian deviates from pi^(3/2)/e")
   end subroutine check_3d_modulated

   ! --------------------------------------------------------------------------
   ! 1D radial grids
   ! --------------------------------------------------------------------------

   !> Uniform (equidistant) radial grid integrates a Gaussian on R^3
   !>
   !> The midpoint rule reaches machine precision once the box covers the
   !> Gaussian's support
   !>
   !> @param[out] error  Test failure
   subroutine test_grid_uniform_radial(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_1d_uniform_type), target :: ugrid
      class(moist_math_grid_1d_type), pointer :: grid
      type(mctc_error), allocatable :: merr

      ! 2000 nodes at dr = 0.01 bohr -> r_max = 20 bohr; the even-Gaussian
      ! midpoint rule is then round-off-limited (~1e-15) -- few enough terms
      ! that the summation round-off stays far below the threshold
      call new_uniform_radial_grid(ugrid, 2000, 0.01_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      grid => ugrid
      call check_radial_gaussian(error, grid, 1.0e-13_wp, "uniform")
      call grid%destroy()
   end subroutine test_grid_uniform_radial

   !> Chebyshev-2 radial grid integrates the same Gaussian
   !>
   !> The non-uniform nodes are driven through the abstract base
   !>
   !> @param[out] error  Test failure
   subroutine test_grid_chebyshev_radial(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_1d_chebyshev_type), target :: cgrid
      class(moist_math_grid_1d_type), pointer :: grid
      type(mctc_error), allocatable :: merr

      ! 400 r- and k-nodes; the Chebyshev-2 rule then resolves the Gaussian to
      ! round-off (~1e-14)
      call new_chebyshev_radial_grid(cgrid, 400, 1.0_wp, 400, 1.0_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      grid => cgrid
      call check_radial_gaussian(error, grid, 1.0e-13_wp, "chebyshev")
      call grid%destroy()
   end subroutine test_grid_chebyshev_radial

   ! --------------------------------------------------------------------------
   ! 3D volume grids
   ! --------------------------------------------------------------------------

   !> Cartesian uniform box grid: integrates a centered 3D Gaussian
   !>
   !> @param[out] error  Test failure
   subroutine test_grid_cartesian_3d(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: cgrid
      class(moist_math_grid_3d_type), pointer :: grid
      type(mctc_error), allocatable :: merr

      ! 128^3 box, dr = 0.11 bohr -> half-width ~7 bohr (Gaussian fully
      ! enclosed); the midpoint rule converges exponentially, so the residual
      ! is pure summation round-off over ~2.1M terms (~4e-12) -- hence the 5e-12
      ! floor rather than a machine-precision threshold
      call new_cartesian_grid_3d(cgrid, 128, 128, 128, 0.11_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      grid => cgrid
      call check_3d_gaussian(error, grid, 5.0e-12_wp, "cartesian")
      call grid%destroy()
   end subroutine test_grid_cartesian_3d

   !> Cartesian grid integrates the anisotropic x^2 exp(-|r|^2) to pi^(3/2)/2
   !>
   !> A check the spherically symmetric Gaussian cannot catch
   !>
   !> @param[out] error  Test failure
   subroutine test_grid_cartesian_aniso(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: cgrid
      class(moist_math_grid_3d_type), pointer :: grid
      real(wp) :: res
      type(mctc_error), allocatable :: merr

      ! Same 128^3 box, so also summation-round-off-limited (~2e-12)
      call new_cartesian_grid_3d(cgrid, 128, 128, 128, 0.11_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      grid => cgrid
      call grid%integrate(x2_gaussian_3d, res)
      call check(error, res, pi*sqrt(pi)/2.0_wp, thr=5.0e-12_wp, &
         & more="cartesian: x^2 exp(-|r|^2) integral deviates from pi^(3/2)/2")
      call grid%destroy()
   end subroutine test_grid_cartesian_aniso

   !> Cartesian grid integrates a higher-frequency, origin-centered field
   !>
   !> Unit Gaussian modulated by cos(2x); the midpoint rule still converges
   !> exponentially, so the result is summation-round-off-limited (~1e-13)
   !>
   !> @param[out] error  Test failure
   subroutine test_grid_cartesian_modulated(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: cgrid
      class(moist_math_grid_3d_type), pointer :: grid
      type(mctc_error), allocatable :: merr

      call new_cartesian_grid_3d(cgrid, 128, 128, 128, 0.11_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      grid => cgrid
      call check_3d_modulated(error, grid, 1.0e-10_wp, "cartesian")
      call grid%destroy()
   end subroutine test_grid_cartesian_modulated

   !> Atom-centered molecular grid (default per-element sizes) integrates a centered Gaussian
   !>
   !> Reference pi^(3/2) is integrand-only, so it does not depend on the
   !> carrier molecule (centered MB16-43/H2 here)
   !>
   !> @param[out] error  Test failure
   subroutine test_grid_molecular_3d(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(mctc_error), allocatable :: merr
      class(moist_math_grid_3d_type), pointer :: grid

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid(mgrid, mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      ! Default per-element sizes are coarse (H -> (50, 302)); the atom-centered
      ! product grid then integrates the off-center Gaussian to ~3e-8
      grid => mgrid
      call check_3d_gaussian(error, grid, 1.0e-7_wp, "molecular")
      call grid%destroy()
   end subroutine test_grid_molecular_3d

   !> Same integrand on the molecular grid built with the uniform constructor
   !>
   !> At very fine sizes (nrad=300, nang=1202), where the atom-centered
   !> product grid resolves the smooth Gaussian down to the summation
   !> round-off floor (~4e-11)
   !>
   !> @param[out] error  Test failure
   subroutine test_grid_molecular_uniform(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(mctc_error), allocatable :: merr
      class(moist_math_grid_3d_type), pointer :: grid

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid_uniform(mgrid, mol, nrad=300, nang=1202, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      grid => mgrid
      call check_3d_gaussian(error, grid, 1.0e-9_wp, "molecular-uniform")
      call grid%destroy()
   end subroutine test_grid_molecular_uniform

   !> Atom-centered molecular grid (very fine uniform sizes), same field
   !>
   !> The same higher-frequency, origin-centered field as
   !> [[test_grid_cartesian_modulated]]; unlike the uniform Cartesian box, the
   !> Becke-partitioned product grid must resolve the cos(2x) oscillation from
   !> two off-origin atomic grids, making this the stricter resolution test
   !>
   !> @param[out] error  Test failure
   subroutine test_grid_molecular_modulated(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(mctc_error), allocatable :: merr
      class(moist_math_grid_3d_type), pointer :: grid

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid_uniform(mgrid, mol, nrad=300, nang=1202, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      ! The fine grid resolves the oscillation to the summation round-off floor
      ! (~5e-12; cos(2x) sign cancellation makes this tighter than the plain
      ! Gaussian on the same grid)
      grid => mgrid
      call check_3d_modulated(error, grid, 1.0e-10_wp, "molecular-modulated")
      call grid%destroy()
   end subroutine test_grid_molecular_modulated

   !> Molecular-grid NUFFT wiring check
   !>
   !> The DC (k = 0) mode of the forward transform must equal the grid
   !> quadrature of the field, sum_j w_j f_j, to the FINUFFT tolerance,
   !> independent of grid resolution, because exp(-i*0.r) = 1; a nonzero
   !> mode is also checked against direct summation so the test validates
   !> the type-3 phase convention, explicit reciprocal targets, baked-in
   !> quadrature weights, k-grid layout, and the plan/execute chain
   !>
   !> @param[out] error  Test failure
   subroutine test_grid_molecular_nufft(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: field(:, :)
      complex(wp), allocatable :: fk(:, :)
      real(wp) :: kvec(3), dc_ref, phase
      complex(wp) :: fk_ref
      integer :: j, npts, k0, ktest

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      !> Clamp the radial extent so the auto-sized reciprocal grid stays small:
      !> a molecular grid otherwise keeps far, low-weight shells that bloat the
      !> bounding box (and the NUFFT grid) far beyond what this test needs
      call new_molecular_grid_uniform(mgrid, mol, nrad=40, nang=110, error=merr, &
         & rmax=6.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      !> Auto-size the reciprocal grid from the point cloud, dr = 0.5 bohr
      call molecular_grid_set_kgrid(mgrid, 0.5_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, mgrid%has_kgrid)
      if (allocated(error)) return
      call check(error, mgrid%npts_k == mgrid%nkx*mgrid%nky*mgrid%nkz)
      if (allocated(error)) return

      call new_molecular_grid_trafo(trafo, mgrid)
      call trafo%prepare(1, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call mgrid%destroy()
         return
      end if

      npts = mgrid%ngrid
      allocate (field(npts, 1), fk(mgrid%npts_k, 1))
      do j = 1, npts
         field(j, 1) = exp(-0.3_wp*sum(mgrid%xyz(:, j)**2)) + 0.2_wp
      end do

      !> Reference DC value = quadrature of the field over the grid
      call mgrid%integrate_field(field(:, 1), dc_ref)

      call trafo%fft_r2k(field, fk, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call mgrid%destroy()
         return
      end if

      !> Locate the k = 0 mode in the CMCL layout and compare
      k0 = 0
      ktest = 0
      do j = 1, mgrid%npts_k
         kvec = mgrid%kpoint(j)
         if (maxval(abs(kvec)) < 1.0e-12_wp) then
            k0 = j
         else if (ktest == 0 .and. kvec(1) > 0.0_wp .and. &
            & max(abs(kvec(2)), abs(kvec(3))) < 1.0e-12_wp) then
            ktest = j
         end if
      end do
      call check(error, k0 > 0, "molecular k-grid has no k = 0 mode")
      if (.not. allocated(error)) &
         call check(error, abs(real(fk(k0, 1), wp) - dc_ref) <= 1.0e-7_wp*abs(dc_ref) + 1.0e-9_wp)
      if (.not. allocated(error)) &
         call check(error, abs(aimag(fk(k0, 1))) <= 1.0e-7_wp*abs(dc_ref) + 1.0e-9_wp)

      !> A nonzero mode catches the type-3 source/target scaling and phase
      call check(error, ktest > 0, "molecular k-grid has no positive x-axis mode")
      if (.not. allocated(error)) then
         kvec = mgrid%kpoint(ktest)
         fk_ref = (0.0_wp, 0.0_wp)
         do j = 1, npts
            phase = dot_product(kvec, mgrid%xyz(:, j) - mgrid%kref)
            fk_ref = fk_ref + mgrid%w(j)*field(j, 1)*cmplx(cos(phase), -sin(phase), wp)
         end do
         call check(error, abs(fk(ktest, 1) - fk_ref) <= 1.0e-7_wp*abs(fk_ref) + 1.0e-9_wp)
      end if

      call trafo%destroy()
      call mgrid%destroy()
      deallocate (field, fk)
   end subroutine test_grid_molecular_nufft

   ! --------------------------------------------------------------------------
   ! Structural / contract tests
   ! --------------------------------------------------------------------------

   !> The Chebyshev FBT adjoints must be exact transposes of the forward/backward transforms
   !>
   !> In the Euclidean inner product: <F a, b> = <a, F^T b>; driven through
   !> the abstract trafo base
   !>
   !> @param[out] error  Test failure
   subroutine test_cheb_radial_fbt_adjoint(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nr = 64, nk = 80
      type(moist_math_grid_1d_chebyshev_type), target :: cgrid
      class(moist_math_grid_1d_trafo_type), allocatable :: trafo
      real(wp) :: ar(nr), bk(nk), fwd(nk), adj(nr)
      real(wp) :: br(nr), ak(nk), bwd(nr), adjk(nk)
      real(wp) :: lhs, rhs
      integer :: i
      type(mctc_error), allocatable :: merr

      call new_chebyshev_radial_grid(cgrid, nr, 1.2_wp, nk, 1.5_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call cgrid%new_trafo(trafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if

      do i = 1, nr
         ar(i) = sin(0.3_wp*real(i, wp))
         br(i) = 1.0_wp/(1.0_wp + 0.1_wp*real(i, wp))
      end do
      do i = 1, nk
         bk(i) = cos(0.2_wp*real(i, wp))
         ak(i) = exp(-0.05_wp*real(i, wp))
      end do

      ! Forward M (r->k) vs its adjoint M^T (k->r): <M ar, bk> = <ar, M^T bk>
      call trafo%fbt_r2k(ar, fwd, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call trafo%fbt_r2k_adj(bk, adj, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      lhs = sum(fwd*bk)
      rhs = sum(ar*adj)
      call check(error, abs(lhs - rhs) <= 1.0e-10_wp*(abs(lhs) + abs(rhs) + 1.0_wp), &
         & "Chebyshev fbt_r2k adjoint identity failed")
      if (allocated(error)) then
         call trafo%destroy()
         call cgrid%destroy()
         return
      end if

      ! Backward B (k->r) vs its adjoint B^T (r->k): <B ak, br> = <ak, B^T br>
      call trafo%fbt_k2r(ak, bwd, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call trafo%fbt_k2r_adj(br, adjk, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      lhs = sum(bwd*br)
      rhs = sum(ak*adjk)
      call check(error, abs(lhs - rhs) <= 1.0e-10_wp*(abs(lhs) + abs(rhs) + 1.0_wp), &
         & "Chebyshev fbt_k2r adjoint identity failed")

      call trafo%destroy()
      call cgrid%destroy()
   end subroutine test_cheb_radial_fbt_adjoint

   !> The batched (`fbt_*_all`) Chebyshev transforms must reproduce the scalar loop
   !>
   !> A column-by-column loop of the scalar transforms to machine precision,
   !> for all four directions; driven through the abstract trafo base so the
   !> test also exercises polymorphic dispatch onto the concrete batched
   !> overrides
   !>
   !> @param[out] error  Test failure
   subroutine test_cheb_radial_fbt_batched(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nr = 64, nk = 80, nb = 4
      type(moist_math_grid_1d_chebyshev_type), target :: cgrid
      class(moist_math_grid_1d_trafo_type), allocatable :: trafo
      real(wp) :: fr(nr, nb), fk(nk, nb)
      real(wp) :: loop_k(nk, nb), batch_k(nk, nb)
      real(wp) :: loop_r(nr, nb), batch_r(nr, nb)
      real(wp) :: relerr
      integer :: i, j
      type(mctc_error), allocatable :: merr

      call new_chebyshev_radial_grid(cgrid, nr, 1.2_wp, nk, 1.5_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call cgrid%new_trafo(trafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if

      do j = 1, nb
         do i = 1, nr
            fr(i, j) = sin(0.3_wp*real(i, wp) + 0.2_wp*real(j, wp))
         end do
         do i = 1, nk
            fk(i, j) = cos(0.2_wp*real(i, wp) - 0.1_wp*real(j, wp))
         end do
      end do

      relerr = 0.0_wp

      ! Forward r -> k
      do j = 1, nb
         call trafo%fbt_r2k(fr(:, j), loop_k(:, j), merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message); return
         end if
      end do
      call trafo%fbt_r2k_all(fr, batch_k, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      relerr = max(relerr, batched_relerr(loop_k, batch_k))

      ! Backward k -> r
      do j = 1, nb
         call trafo%fbt_k2r(fk(:, j), loop_r(:, j), merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message); return
         end if
      end do
      call trafo%fbt_k2r_all(fk, batch_r, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      relerr = max(relerr, batched_relerr(loop_r, batch_r))

      ! Forward adjoint k -> r
      do j = 1, nb
         call trafo%fbt_r2k_adj(fk(:, j), loop_r(:, j), merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message); return
         end if
      end do
      call trafo%fbt_r2k_adj_all(fk, batch_r, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      relerr = max(relerr, batched_relerr(loop_r, batch_r))

      ! Backward adjoint r -> k
      do j = 1, nb
         call trafo%fbt_k2r_adj(fr(:, j), loop_k(:, j), merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message); return
         end if
      end do
      call trafo%fbt_k2r_adj_all(fr, batch_k, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      relerr = max(relerr, batched_relerr(loop_k, batch_k))

      call trafo%destroy()
      call cgrid%destroy()

      call check(error, relerr <= 1.0e-12_wp, &
         & "Chebyshev batched transform disagrees with the scalar loop")
   end subroutine test_cheb_radial_fbt_batched

   !> Largest column entry where the batched result deviates from the scalar loop
   !>
   !> Normalised by the loop magnitude (so GEMM-vs-GEMV rounding stays
   !> O(epsilon)); a non-finite entry returns `huge(1.0_wp)` instead of
   !> feeding NaN into `maxval`: a NaN comparison is always false, so
   !> gfortran's pairwise reduction can silently keep the running finite
   !> extremum and hide the corrupted value from the `relerr <= thr` check
   !> at the call site
   !>
   !> @param[in] loop   Reference result from looping the scalar transform per column
   !> @param[in] batch  Result from the batched transform
   pure function batched_relerr(loop, batch) result(relerr)
      !> Reference result from looping the scalar transform per column
      real(wp), intent(in) :: loop(:, :)
      !> Result from the batched transform
      real(wp), intent(in) :: batch(:, :)
      !> Relative deviation
      real(wp) :: relerr

      if (.not. all(ieee_is_finite(loop)) .or. .not. all(ieee_is_finite(batch))) then
         relerr = huge(1.0_wp)
         return
      end if
      relerr = maxval(abs(loop - batch))/max(1.0_wp, maxval(abs(loop)))
   end function batched_relerr

   !> Destroying the Chebyshev grid (and its trafo) twice must be a safe no-op
   !>
   !> @param[out] error  Test failure
   subroutine test_cheb_radial_destroy(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_1d_chebyshev_type), target :: grid
      class(moist_math_grid_1d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr

      call new_chebyshev_radial_grid(grid, 32, 1.0_wp, 32, 1.0_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call grid%new_trafo(trafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if

      call trafo%destroy()
      call trafo%destroy()
      call grid%destroy()
      call grid%destroy()
      call check(error, grid%npts == 0 .and. grid%nk == 0, &
         & "Chebyshev grid not reset after destroy")
   end subroutine test_cheb_radial_destroy

   !> %integrate(f=1, result) should return sum(weights) on the molecular grid
   !>
   !> @param[out] error  Test failure
   subroutine test_mol_grid_integrate_const(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: result, expected

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid(grid, mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call grid%integrate(one_3d, result)
      expected = sum(grid%w)

      call check(error, abs(result - expected) < 1.0e-10_wp, &
         & "integrate(f=1) != sum(weights)")
      call grid%destroy()
   end subroutine test_mol_grid_integrate_const

   !> HandyMod molecular grid should construct finite data
   !>
   !> @param[out] error  Test failure
   subroutine test_mol_grid_qc_handymod_smoke(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: total_weight

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid_uniform_qc_handymod(grid, mol, nrad=2000, &
         & lebedev_degree=77, error=merr, rmin=0.0_wp, rmax=6.0_wp, m=0.1_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call check(error, grid%ngrid > 0, "qc-HandyMod molecular grid has no points")
      if (allocated(error)) return
      call check(error, all(ieee_is_finite(grid%xyz)), &
         & "qc-HandyMod molecular grid has non-finite coordinates")
      if (allocated(error)) return
      call check(error, all(ieee_is_finite(grid%w)), &
         & "qc-HandyMod molecular grid has non-finite weights")
      if (allocated(error)) return
      call check(error, minval(grid%w) >= -1.0e-12_wp, &
         & "qc-HandyMod molecular grid has strongly negative weights")
      if (allocated(error)) return

      total_weight = sum(grid%w)
      call check(error, ieee_is_finite(total_weight) .and. total_weight > 0.0_wp, &
         & "qc-HandyMod molecular grid constant integral is not finite and positive")
      call grid%destroy()
   end subroutine test_mol_grid_qc_handymod_smoke

   !> The new qc-HandyMod path must not perturb existing uniform Chebyshev grids
   !>
   !> @param[out] error  Test failure
   subroutine test_mol_grid_uniform_regression(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type) :: before_grid, qc_grid, after_grid
      type(mctc_error), allocatable :: merr
      integer :: npts_before
      real(wp) :: weight_before

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)

      call new_molecular_grid_uniform(before_grid, mol, nrad=12, nang=26, error=merr, rmax=3.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      npts_before = before_grid%ngrid
      weight_before = sum(before_grid%w)

      call new_molecular_grid_uniform_qc_handymod(qc_grid, mol, nrad=16, &
         & lebedev_degree=7, error=merr, rmin=0.0_wp, rmax=3.0_wp, m=0.1_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call before_grid%destroy()
         return
      end if

      call new_molecular_grid_uniform(after_grid, mol, nrad=12, nang=26, error=merr, rmax=3.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call before_grid%destroy()
         call qc_grid%destroy()
         return
      end if

      call check(error, after_grid%ngrid == npts_before, &
         & "existing uniform molecular grid point count changed after qc-HandyMod construction")
      if (allocated(error)) return
      call check(error, abs(sum(after_grid%w) - weight_before) <= 1.0e-12_wp, &
         & "existing uniform molecular grid total weight changed after qc-HandyMod construction")

      call before_grid%destroy()
      call qc_grid%destroy()
      call after_grid%destroy()
   end subroutine test_mol_grid_uniform_regression

   !> destroy() must be safely callable twice on the molecular grid
   !>
   !> @param[out] error  Test failure
   subroutine test_mol_grid_destroy(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type) :: grid
      type(mctc_error), allocatable :: merr

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid(grid, mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call grid%destroy()
      call grid%destroy()

      call check(error, grid%ngrid == 0, "npts should be 0 after destroy")
      if (allocated(error)) return
      call check(error, .not. allocated(grid%xyz), "xyz should be deallocated")
   end subroutine test_mol_grid_destroy

   !> All retained weights should have |w| >= default threshold (1e-14)
   !>
   !> @param[out] error  Test failure
   subroutine test_mol_grid_pruning(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type) :: grid
      type(mctc_error), allocatable :: merr
      integer :: i

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid(grid, mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      do i = 1, grid%ngrid
         call check(error, abs(grid%w(i)) >= 1.0e-14_wp, &
            & "Point with |weight| below threshold survived pruning")
         if (allocated(error)) exit
      end do
      call grid%destroy()
   end subroutine test_mol_grid_pruning

end module test_math_grid
