!> Analytical volume integrals on all four 3D grid variants
!>
!> - Cartesian point, Cartesian Gaussian, molecular point and molecular Gaussian
!> - Both integrand callbacks and tabulated fields against independent references
!> - Whole-space references; unbounded molecular radii and a large Cartesian box
!> - Smooth screened Lennard-Jones and Yukawa terms with a regularized distance
!> - Cases and grid settings shared with the rigid-motion invariance suite
!> - Two-center Gaussians under every partition scheme and the per-element defaults
!> - Rotated water integrals converging with the angular degree
!> - Cartesian box volume and cell volume of a tabulated field
!> - Gaussian width law on the Cartesian and molecular Gaussian grids
module test_math_grid_3d_integration
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: check_moist_error, get_cartesian_gaussian_grid, get_becke_recipe
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type, integrand_3d
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, &
      & new_cartesian_point_grid, new_cartesian_gaussian_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, &
      & new_molecular_point_grid, new_molecular_gaussian_grid, partition_becke, partition_ssf, &
      & partition_pvoronoi
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type, &
      & moist_math_grid_atomic_shell_constant_type, new_constant_shell_policy, &
      & moist_math_grid_atomic_recipe_override_type, default_element_recipes
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_gauss_legendre_type, &
      & new_gauss_legendre_rule
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_becke_type, new_becke_mapping
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_generator_lebedev_type, new_lebedev_generator
   implicit none(type, external)
   private

   public :: collect_math_grid_3d_integration
   public :: integration_case_type, make_cases, make_grid, grid_names, integral_tol

   !> Absolute and relative integral tolerances, each below 1e-12
   real(wp), parameter :: integral_tol = 5.0e-13_wp
   !> Integral of exp(-|r|**2) over all space
   real(wp), parameter :: gaussian_volume = pi*sqrt(pi)
   !> Off-axis function centers, independent of the carrier molecule (bohr)
   real(wp), parameter :: centers(3, 3) = reshape([ &
      & 0.37_wp, -0.29_wp, 0.43_wp, -0.81_wp, 0.53_wp, -0.17_wp, &
      & 0.63_wp, 0.71_wp, -0.59_wp], [3, 3])
   !> Principal Gaussian exponents (bohr**-2)
   real(wp), parameter :: axes(3) = [0.4_wp, 1.1_wp, 2.3_wp]
   !> Wave vector in the Gaussian principal axes (bohr**-1)
   real(wp), parameter :: wave(3) = [2.1_wp, -1.7_wp, 1.3_wp]
   !> Phase of the anisotropic oscillatory field (radians)
   real(wp), parameter :: phase = 0.4_wp
   !> Exponents of the signed mixture (bohr**-2)
   real(wp), parameter :: exponents(3) = [0.35_wp, 1.2_wp, 3.0_wp]
   !> Positive and negative mixture coefficients
   real(wp), parameter :: coefficients(3) = [0.9_wp, -0.6_wp, 0.4_wp]
   !> Radius of the diffuse Gaussian shell (bohr)
   real(wp), parameter :: shell_radius = 6.0_wp
   !> Smooth interaction core radius (bohr)
   real(wp), parameter :: core = 3.5_wp
   !> Lennard-Jones size (bohr)
   real(wp), parameter :: sigma = 5.0_wp
   !> Lennard-Jones well depth
   real(wp), parameter :: well_depth = 0.01_wp
   !> Exponential screening of both interaction terms (bohr**-1)
   real(wp), parameter :: screening = 2.0_wp
   !> Per-center Lennard-Jones coefficients
   real(wp), parameter :: lj_coefficients(3) = [0.6_wp, 1.1_wp, 0.8_wp]
   !> Per-center signed Yukawa charges
   real(wp), parameter :: charges(3) = [0.7_wp, -0.2_wp, 0.4_wp]
   !> Cartesian point count per axis
   integer, parameter :: cartesian_points = 160
   !> Cartesian spacing (bohr); box half-width 24 bohr
   real(wp), parameter :: cartesian_spacing = 0.3_wp
   !> Gauss-Legendre radial point count
   integer, parameter :: radial_points = 240
   !> Lebedev polynomial exactness degree
   integer, parameter :: angular_degree = 131
   !> Scale of the Becke map r = scale*(1+x)/(1-x) onto [0, infinity) (bohr)
   real(wp), parameter :: radial_scale = 1.0_wp
   !> Labels distinguish point and Gaussian variants in assertion failures
   character(len=20), parameter :: grid_names(4) = [character(len=20) :: &
      & "Cartesian point", "Cartesian Gaussian", "molecular point", "molecular Gaussian"]

   !> Integrand and its independently calculated whole-space integral
   type :: integration_case_type
      !> Integrand name used in assertion failures
      character(len=32) :: name
      !> Analytical function sampled by the grid
      procedure(integrand_3d), pointer, nopass :: f => null()
      !> Exact or independently converged integral
      real(wp) :: reference
   end type integration_case_type

contains

   !> Register geometries; each exercises every integrand on all four grids
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_math_grid_3d_integration(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("single_atom", test_single_atom), &
         & new_unittest("asymmetric_trimer", test_asymmetric_trimer), &
         & new_unittest("partition_integrals", test_partition_integrals), &
         & new_unittest("multicenter_integrals", test_multicenter_integrals), &
         & new_unittest("rotational_convergence", test_rotational_convergence), &
         & new_unittest("integrate_constant_is_box_volume", test_integrate_const), &
         & new_unittest("gaussian_width_law", test_gaussian_width_law)]
   end subroutine collect_math_grid_3d_integration

   !> Translated single-center quadrature with off-center integrands
   subroutine test_single_atom(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [8], reshape([0.19_wp, -0.13_wp, 0.23_wp], [3, 1]))
      call check_integrals(mol, error)
   end subroutine test_single_atom

   !> Heteronuclear partition on an asymmetric, off-axis three-center geometry
   subroutine test_asymmetric_trimer(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [8, 1, 6], reshape([-0.6_wp, -0.2_wp, 0.1_wp, &
         & 0.7_wp, 0.3_wp, -0.25_wp, 0.1_wp, 1.2_wp, 0.4_wp], [3, 3]))
      call check_integrals(mol, error)
   end subroutine test_asymmetric_trimer

   !> Integrate analytic two-center Gaussians with every partition scheme
   subroutine test_partition_integrals(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      integer :: scheme, ischeme, a, b, level, nlevels
      integer, parameter :: schemes(3) = [partition_becke, partition_ssf, partition_pvoronoi]
      real(wp) :: exact, approx, max_error(2)
      character(len=120) :: msg

      call make_water(mol)
      do ischeme = 1, size(schemes)
         scheme = schemes(ischeme)
         nlevels = merge(2, 1, scheme == partition_pvoronoi)
         max_error = 0.0_wp
         do level = 1, nlevels
            if (level == 1) then
               call get_becke_recipe(recipe, 100, 0.5_wp, 83, merr)
            else
               call get_becke_recipe(recipe, 160, 0.5_wp, 131, merr)
            end if
            if (.not. allocated(merr)) call new_molecular_point_grid(grid, merr, recipe=recipe, &
               & partition=scheme, reciprocal=.false., weight_threshold=0.0_wp)
            if (.not. allocated(merr)) call grid%update(mol, merr)
            call check_moist_error(error, merr, "integration grid")
            if (allocated(error)) return
            do a = 1, mol%nat
               do b = a, mol%nat
                  exact = gauss_pair_exact(mol, a, b, 0.8_wp)
                  approx = gauss_pair_integral(grid, mol, a, b, 0.8_wp)
                  max_error(level) = max(max_error(level), abs(approx - exact)/exact)
                  if (level /= nlevels) cycle
                  write (msg, "(a,i0,a,i0,a,i0,a,es10.3)") "scheme ", scheme, ", Gaussian (", a, ",", b, &
                     & ") relative error ", abs(approx - exact)/exact
                  call check(error, abs(approx - exact) < 1.0e-5_wp*exact, trim(msg))
                  if (allocated(error)) return
               end do
            end do
         end do
         if (nlevels == 2) then
            call check(error, max_error(2) < 0.5_wp*max_error(1), "power partition quadrature must converge with refinement")
            if (allocated(error)) return
         end if
      end do
   end subroutine test_partition_integrals

   !> One- and two-center Gaussian integrals over a water grid with the per-element defaults
   subroutine test_multicenter_integrals(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: alpha = 0.8_wp
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      character(len=120) :: msg
      real(wp) :: exact, approx
      integer :: a, b

      call make_water(mol)
      call default_element_recipes(recipe, overrides, merr)
      if (.not. allocated(merr)) call new_molecular_point_grid(grid, merr, recipe=recipe, overrides=overrides, &
         & reciprocal=.false.)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      call check_moist_error(error, merr, "water")
      if (allocated(error)) return
      do a = 1, mol%nat
         do b = a, mol%nat
            exact = gauss_pair_exact(mol, a, b, alpha)
            approx = gauss_pair_integral(grid, mol, a, b, alpha)
            write (msg, "(a,i0,a,i0,a,es10.3)") "Gaussian (", a, ", ", b, ") relative error ", &
               & abs(approx - exact)/exact
            call check(error, abs(approx - exact) <= 2.0e-6_wp*exact, trim(msg))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_multicenter_integrals

   !> Rigidly rotating molecule and integrand changes the integral less at a higher angular degree
   subroutine test_rotational_convergence(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: alpha = 1.2_wp
      integer, parameter :: degrees(2) = [11, 29]
      real(wp), parameter :: angles(3, 3) = reshape([0.3_wp, 1.1_wp, -0.7_wp, &
         & 2.0_wp, -0.4_wp, 0.9_wp, -1.3_wp, 0.6_wp, 2.4_wp], [3, 3])
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol, rotated
      type(mctc_error), allocatable :: merr
      real(wp) :: spread_deg(2), reference, value, rot(3, 3), c(3)
      character(len=120) :: msg
      integer :: id, ir, iat

      call make_water(mol)
      c = centroid(mol)
      do id = 1, 2
         call get_becke_recipe(recipe, 40, 0.5_wp, degrees(id), merr)
         if (.not. allocated(merr)) call new_molecular_point_grid(grid, merr, recipe=recipe, reciprocal=.false.)
         if (.not. allocated(merr)) call grid%update(mol, merr)
         call check_moist_error(error, merr, "reference orientation")
         if (allocated(error)) return
         reference = gauss_pair_integral(grid, mol, 1, 2, alpha) + gauss_pair_integral(grid, mol, 2, 3, alpha)
         spread_deg(id) = 0.0_wp
         do ir = 1, size(angles, 2)
            rot = rotation(angles(:, ir))
            rotated = mol
            do iat = 1, mol%nat
               rotated%xyz(:, iat) = c + matmul(rot, mol%xyz(:, iat) - c)
            end do
            call grid%update(rotated, merr)
            call check_moist_error(error, merr, "rotated orientation")
            if (allocated(error)) return
            value = gauss_pair_integral(grid, rotated, 1, 2, alpha) &
               & + gauss_pair_integral(grid, rotated, 2, 3, alpha)
            spread_deg(id) = max(spread_deg(id), abs(value - reference)/reference)
         end do
      end do
      write (msg, "(a,es10.3,a,es10.3)") "rotational spread, degree 11: ", spread_deg(1), &
         & ", degree 29: ", spread_deg(2)
      call check(error, spread_deg(2) < 0.1_wp*spread_deg(1) .and. spread_deg(2) < 1.0e-6_wp, trim(msg))

   contains

      !> Rotation matrix from three Euler angles (z-y-z)
      !>
      !> @param[in] ang  Euler angles (rad)
      pure function rotation(ang) result(r)
         !> Euler angles
         real(wp), intent(in) :: ang(3)
         !> Rotation matrix
         real(wp) :: r(3, 3)

         real(wp) :: rz1(3, 3), ry(3, 3), rz2(3, 3)

         rz1 = reshape([cos(ang(1)), sin(ang(1)), 0.0_wp, -sin(ang(1)), cos(ang(1)), 0.0_wp, &
            & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
         ry = reshape([cos(ang(2)), 0.0_wp, -sin(ang(2)), 0.0_wp, 1.0_wp, 0.0_wp, &
            & sin(ang(2)), 0.0_wp, cos(ang(2))], [3, 3])
         rz2 = reshape([cos(ang(3)), sin(ang(3)), 0.0_wp, -sin(ang(3)), cos(ang(3)), 0.0_wp, &
            & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
         r = matmul(rz2, matmul(ry, rz1))
      end function rotation
   end subroutine test_rotational_convergence

   !> Integral of 1 dV over the box equals the box volume, dispatched polymorphically
   !>
   !> - Field integral of 1 + x^2 over the 8^3 box of spacing 0.25 bohr, nodes
   !>   at the cell centers +-0.125 .. +-0.875 bohr: box volume 8 times
   !>   (1 + mean x^2) = 8*(1 + 0.328125) = 10.625
   subroutine test_integrate_const(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_type), pointer :: g
      real(wp) :: res
      type(mctc_error), allocatable :: merr

      call get_cartesian_gaussian_grid(grid, 8, 8, 8, 0.25_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      g => grid
      call g%integrate(one, res)
      call check(error, abs(res - grid%get_box_volume()) < 1.0e-10_wp*grid%get_box_volume())
      if (allocated(error)) return
      call g%integrate_field(1.0_wp + grid%xyz(1, :)**2, res, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, abs(res - 10.625_wp) < 1.0e-12_wp, &
         & "Cartesian field integral must use the cell volume")
      call grid%destroy()
   end subroutine test_integrate_const

   !> Gaussian widths follow the spacing and the weights at every geometry
   !>
   !> - Cartesian: xi0 = xi0_factor/dr at every node
   !> - Molecular: xi0**3*w = xi0_factor**3 at every retained node
   subroutine test_gaussian_width_law(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Width scale of both grids
      real(wp), parameter :: width_factor = 1.2_wp
      !> Cartesian spacing (bohr)
      real(wp), parameter :: spacing = 0.25_wp
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(structure_type) :: mols(2)
      type(mctc_error), allocatable :: merr
      character(len=32) :: label
      integer :: step, i

      call make_water(mols(1))
      call displace(mols(1), mols(2))
      call new_cartesian_gaussian_grid(cart, merr, nx=16, ny=16, nz=16, dr=spacing, xi0_factor=width_factor, &
         & margin=0.0_wp)
      if (.not. allocated(merr)) call get_becke_recipe(recipe, 10, 0.5_wp, 7, merr, rcut_upper=5.0_wp)
      if (.not. allocated(merr)) call new_molecular_gaussian_grid(grid, merr, recipe=recipe, &
         & xi0_factor=width_factor, ssf_a=0.64_wp, reciprocal=.false.)
      call check_moist_error(error, merr, "construction")
      if (allocated(error)) return
      do step = 1, size(mols)
         write (label, "(a,i0)") "geometry ", step
         call cart%update(mols(step), merr)
         if (.not. allocated(merr)) call grid%update(mols(step), merr)
         call check_moist_error(error, merr, trim(label)//": update")
         if (allocated(error)) return
         call check(error, allocated(cart%xi0), trim(label)//": Cartesian widths")
         if (allocated(error)) return
         call check(error, size(cart%xi0) == cart%ngrid, trim(label)//": Cartesian width count")
         if (allocated(error)) return
         do i = 1, cart%ngrid
            call check(error, cart%xi0(i), width_factor/spacing, thr=1.0e-13_wp, &
               & more=trim(label)//": Cartesian width law")
            if (allocated(error)) return
         end do
         call check(error, allocated(grid%xi0), trim(label)//": molecular widths")
         if (allocated(error)) return
         call check(error, size(grid%xi0) == grid%ngrid, trim(label)//": molecular width count")
         if (allocated(error)) return
         do i = 1, grid%ngrid
            call check(error, grid%xi0(i)**3*grid%w(i), width_factor**3, thr=1.0e-13_wp, &
               & more=trim(label)//": molecular width law")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_gaussian_width_law

   !> Shared function catalogue for every geometry and grid variant
   !>
   !> @param[out] cases Functions and whole-space reference values
   subroutine make_cases(cases)
      !> Functions and independently calculated integrals
      type(integration_case_type), allocatable, intent(out) :: cases(:)
      real(wp) :: anisotropic_volume, pair_volume, shell_volume

      anisotropic_volume = gaussian_volume/sqrt(product(axes))
      pair_volume = (pi/2.0_wp)**1.5_wp*exp(-0.48_wp*sum((centers(:, 1) - centers(:, 2))**2))
      shell_volume = 4.0_wp*pi*(sqrt(pi)*(0.25_wp + 0.5_wp*shell_radius**2)*(1.0_wp + erf(shell_radius)) &
         & + 0.5_wp*shell_radius*exp(-shell_radius**2))
      cases = [ &
         & integration_case_type("Gaussian", gaussian, gaussian_volume), &
         & integration_case_type("translated Gaussian", shifted_gaussian, gaussian_volume/1.3_wp**1.5_wp), &
         & integration_case_type("two-center Gaussian product", gaussian_product, pair_volume), &
         & integration_case_type("quadratic Gaussian moment", quadratic_gaussian, gaussian_volume/2.0_wp), &
         & integration_case_type("mixed sixth-order moment", mixed_gaussian, gaussian_volume/8.0_wp), &
         & integration_case_type("cosine-modulated Gaussian", modulated_gaussian, gaussian_volume*exp(-1.0_wp)), &
         & integration_case_type("rotated anisotropic Gaussian", anisotropic_gaussian, anisotropic_volume), &
         & integration_case_type("anisotropic oscillatory Gaussian", oscillatory_gaussian, &
         & anisotropic_volume*exp(-sum(wave**2/(4.0_wp*axes)))*cos(phase)), &
         & integration_case_type("odd Gaussian moment", odd_gaussian, 0.0_wp), &
         & integration_case_type("signed multiscale mixture", gaussian_mixture, &
         & gaussian_volume*sum(coefficients/exponents**1.5_wp)), &
         & integration_case_type("diffuse Gaussian shell", gaussian_shell, shell_volume), &
         & integration_case_type("Gaussian-screened Coulomb", gaussian_coulomb, &
         & 2.0_wp*pi*0.7_wp/(0.8_wp*sqrt(0.8_wp + 0.7_wp**2))), &
         & integration_case_type("smooth Lennard-Jones + Yukawa", smooth_interaction, interaction_reference()) &
         & ]
   end subroutine make_cases

   !> Check both integration entry points against each independent reference
   !>
   !> A field of the wrong length is rejected, and the next valid call succeeds
   !>
   !> @param[in] mol Carrier molecule
   !> @param[out] error Test failure
   subroutine check_integrals(mol, error)
      !> Geometry supplying grid centers
      type(structure_type), intent(in) :: mol
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      class(moist_math_grid_3d_type), allocatable :: grid
      type(integration_case_type), allocatable :: cases(:)
      real(wp), allocatable :: values(:)
      real(wp) :: callback_integral, field_integral
      type(mctc_error), allocatable :: merr
      character(len=:), allocatable :: label
      integer :: kind, icase, i

      call make_cases(cases)
      do kind = 1, 4
         call make_grid(kind, mol, grid, error)
         if (allocated(error)) return
         allocate (values(grid%ngrid), source=1.0_wp)
         ! A short field is rejected; the valid calls below still succeed
         call grid%integrate_field(values(2:), field_integral, merr)
         call check(error, allocated(merr), trim(grid_names(kind))//": integrate_field accepted a short field")
         if (allocated(error)) return
         do icase = 1, size(cases)
            label = trim(grid_names(kind))//": "//trim(cases(icase)%name)
            call check(error, ieee_is_finite(cases(icase)%reference), label//": non-finite reference")
            if (allocated(error)) return
            do i = 1, grid%ngrid
               values(i) = cases(icase)%f(grid%point(i))
            end do
            call check(error, all(ieee_is_finite(values)), label//": non-finite samples")
            if (allocated(error)) return
            call grid%integrate(cases(icase)%f, callback_integral)
            call grid%integrate_field(values, field_integral, merr)
            if (allocated(merr)) then
               call test_failed(error, merr%message)
               return
            end if
            call check(error, callback_integral, cases(icase)%reference, thr_abs=integral_tol, &
               & thr_rel=integral_tol, more=label//": integrate")
            if (allocated(error)) return
            call check(error, field_integral, cases(icase)%reference, thr_abs=integral_tol, &
               & thr_rel=integral_tol, more=label//": integrate_field")
            if (allocated(error)) return
         end do
         deallocate (values)
      end do
   end subroutine check_integrals

   !> Build a point or Gaussian variant with common quadrature settings
   !>
   !> @param[in] kind Grid variant, 1..4
   !> @param[in] mol Carrier geometry
   !> @param[out] grid Updated grid
   !> @param[out] error Test failure
   subroutine make_grid(kind, mol, grid, error)
      !> Grid variant
      integer, intent(in) :: kind
      !> Geometry supplying grid centers
      type(structure_type), intent(in) :: mol
      !> Updated grid
      class(moist_math_grid_3d_type), allocatable, intent(out) :: grid
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_cartesian_type), allocatable :: cart
      type(moist_math_grid_3d_molecular_type), allocatable :: molecular
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(mctc_error), allocatable :: merr

      select case (kind)
      case (1, 2)
         allocate (cart)
         if (kind == 1) then
            call new_cartesian_point_grid(cart, merr, nx=cartesian_points, ny=cartesian_points, &
               & nz=cartesian_points, dr=cartesian_spacing)
         else
            call new_cartesian_gaussian_grid(cart, merr, nx=cartesian_points, ny=cartesian_points, &
               & nz=cartesian_points, dr=cartesian_spacing, xi0_factor=1.2_wp)
         end if
         if (.not. allocated(merr)) call cart%update(mol, merr)
         call check_moist_error(error, merr, trim(grid_names(kind))//" setup")
         if (allocated(error)) return
         call move_alloc(cart, grid)
      case default
         call make_recipe(recipe, merr)
         allocate (molecular)
         if (.not. allocated(merr)) then
            if (kind == 3) then
               call new_molecular_point_grid(molecular, merr, recipe=recipe, reciprocal=.false., &
                  & weight_threshold=0.0_wp, pruning_threshold=0.0_wp)
            else
               call new_molecular_gaussian_grid(molecular, merr, recipe=recipe, reciprocal=.false., &
                  & weight_threshold=0.0_wp, pruning_threshold=0.0_wp, xi0_factor=1.2_wp)
            end if
         end if
         molecular%nthreads = 2
         if (.not. allocated(merr)) call molecular%update(mol, merr)
         call check_moist_error(error, merr, trim(grid_names(kind))//" setup")
         if (allocated(error)) return
         call move_alloc(molecular, grid)
      end select
   end subroutine make_grid

   !> Positive Gauss-Legendre times Lebedev rule on every atom
   !>
   !> Rational radial mapping resolves both the near-nucleus partition and tails
   !>
   !> @param[out] recipe Atomic quadrature recipe
   !> @param[out] error Invalid mapping or shell settings
   subroutine make_recipe(recipe, error)
      !> Shared recipe for every element
      type(moist_math_grid_atomic_recipe_type), intent(out) :: recipe
      !> Construction failure
      type(mctc_error), allocatable, intent(out) :: error
      type(moist_math_grid_radial_rule_gauss_legendre_type) :: rule
      type(moist_math_grid_radial_mapping_becke_type) :: mapping
      type(moist_math_grid_atomic_shell_constant_type) :: shells
      type(moist_math_grid_angular_generator_lebedev_type) :: angular

      call new_gauss_legendre_rule(rule)
      call new_becke_mapping(mapping, error, scale=radial_scale)
      if (allocated(error)) return
      call new_constant_shell_policy(shells, angular_degree, error)
      if (allocated(error)) return
      call new_lebedev_generator(angular, positive_weights_only=.true.)
      recipe%radial%npts = radial_points
      allocate (recipe%radial%rule, source=rule)
      allocate (recipe%radial%mapping, source=mapping)
      allocate (recipe%angular, source=angular)
      allocate (recipe%shells, source=shells)
   end subroutine make_recipe

   !> Water, oxygen near the origin
   !>
   !> @param[out] mol Structure
   subroutine make_water(mol)
      !> Structure
      type(structure_type), intent(out) :: mol

      call new(mol, [8, 1, 1], reshape([0.05_wp, -0.02_wp, 0.12_wp, &
         & 0.0_wp, 1.43_wp, -0.98_wp, &
         & 0.0_wp, -1.43_wp, -0.98_wp], [3, 3]))
   end subroutine make_water

   !> Deterministic small displacement of every nucleus
   !>
   !> @param[in] mol Structure
   !> @param[out] moved Structure with displaced nuclei (at most 0.06 bohr per axis)
   subroutine displace(mol, moved)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Displaced structure
      type(structure_type), intent(out) :: moved

      integer :: iat, k

      moved = mol
      do iat = 1, mol%nat
         do k = 1, 3
            moved%xyz(k, iat) = mol%xyz(k, iat) + 0.06_wp*sin(real(7*iat + 3*k, wp))
         end do
      end do
   end subroutine displace

   !> Centroid of the nuclei
   !>
   !> @param[in] mol Structure
   pure function centroid(mol) result(c)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Centroid (bohr)
      real(wp) :: c(3)

      c = sum(mol%xyz, dim=2)/real(mol%nat, wp)
   end function centroid

   !> Sum of w*f over the grid for a Gaussian product centered on two nuclei
   !>
   !> `f = exp(-alpha |r - R_a|**2 - alpha |r - R_b|**2)`, or a single
   !> Gaussian on `R_a` when `a == b`
   !>
   !> @param[in] grid Realized grid
   !> @param[in] mol Structure
   !> @param[in] a First nucleus
   !> @param[in] b Second nucleus
   !> @param[in] alpha Exponent (1/bohr^2)
   function gauss_pair_integral(grid, mol, a, b, alpha) result(total)
      !> Realized grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Structure
      type(structure_type), intent(in) :: mol
      !> First nucleus
      integer, intent(in) :: a
      !> Second nucleus
      integer, intent(in) :: b
      !> Exponent
      real(wp), intent(in) :: alpha
      !> Quadrature sum
      real(wp) :: total

      integer :: i
      real(wp) :: ra, rb

      total = 0.0_wp
      do i = 1, grid%ngrid
         ra = sum((grid%xyz(:, i) - mol%xyz(:, a))**2)
         rb = sum((grid%xyz(:, i) - mol%xyz(:, b))**2)
         if (a == b) then
            total = total + grid%w(i)*exp(-alpha*ra)
         else
            total = total + grid%w(i)*exp(-alpha*(ra + rb))
         end if
      end do
   end function gauss_pair_integral

   !> Exact integral of the Gaussian product of `gauss_pair_integral`
   !>
   !> @param[in] mol Structure
   !> @param[in] a First nucleus
   !> @param[in] b Second nucleus
   !> @param[in] alpha Exponent (1/bohr^2)
   pure function gauss_pair_exact(mol, a, b, alpha) result(exact)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> First nucleus
      integer, intent(in) :: a
      !> Second nucleus
      integer, intent(in) :: b
      !> Exponent
      real(wp), intent(in) :: alpha
      !> Exact integral
      real(wp) :: exact

      if (a == b) then
         exact = (pi/alpha)**1.5_wp
      else
         exact = (pi/(2.0_wp*alpha))**1.5_wp*exp(-0.5_wp*alpha*sum((mol%xyz(:, a) - mol%xyz(:, b))**2))
      end if
   end function gauss_pair_exact

   !> Constant integrand, a volume probe
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function one(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value

      value = 1.0_wp + 0.0_wp*r(1)
   end function one

   !> Unit Gaussian
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function gaussian(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value

      value = exp(-sum(r**2))
   end function gaussian

   !> Off-center Gaussian with a nonunit exponent
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function shifted_gaussian(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value

      value = exp(-1.3_wp*sum((r - centers(:, 1))**2))
   end function shifted_gaussian

   !> Product of two Gaussians with distinct centers and exponents
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function gaussian_product(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value

      value = exp(-0.8_wp*sum((r - centers(:, 1))**2) - 1.2_wp*sum((r - centers(:, 2))**2))
   end function gaussian_product

   !> Second Cartesian moment of a Gaussian
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function quadratic_gaussian(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value

      value = r(1)**2*exp(-sum(r**2))
   end function quadratic_gaussian

   !> Mixed sixth-order Gaussian moment, x**2*y**2*z**2
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function mixed_gaussian(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value

      value = product(r**2)*exp(-sum(r**2))
   end function mixed_gaussian

   !> Gaussian modulated by cos(2*x)
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function modulated_gaussian(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value

      value = exp(-sum(r**2))*cos(2.0_wp*r(1))
   end function modulated_gaussian

   !> Orthonormal coordinates about an off-axis center
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function principal_coordinates(r) result(u)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Rotated displacement (bohr)
      real(wp) :: u(3)
      real(wp) :: x(3)

      x = r - centers(:, 1)
      u = [2.0_wp*x(1) + x(2) + 2.0_wp*x(3), &
         & 2.0_wp*x(1) - 2.0_wp*x(2) - x(3), x(1) + 2.0_wp*x(2) - 2.0_wp*x(3)]/3.0_wp
   end function principal_coordinates

   !> Tilted Gaussian with three distinct principal widths
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function anisotropic_gaussian(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value

      value = exp(-sum(axes*principal_coordinates(r)**2))
   end function anisotropic_gaussian

   !> Tilted anisotropic Gaussian with a nonzero oscillation phase
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function oscillatory_gaussian(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value
      real(wp) :: u(3)

      u = principal_coordinates(r)
      value = exp(-sum(axes*u**2))*cos(dot_product(wave, u) + phase)
   end function oscillatory_gaussian

   !> Signed odd moment with an exactly zero integral
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function odd_gaussian(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value
      real(wp) :: u(3)

      u = principal_coordinates(r)
      value = (u(1) + 2.0_wp*u(2) - u(3))*exp(-sum(axes*u**2))
   end function odd_gaussian

   !> Signed mixture of narrow and diffuse Gaussians at three centers
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function gaussian_mixture(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value
      integer :: i

      value = 0.0_wp
      do i = 1, size(exponents)
         value = value + coefficients(i)*exp(-exponents(i)*sum((r - centers(:, i))**2))
      end do
   end function gaussian_mixture

   !> Diffuse shell centered away from the carrier nuclei
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function gaussian_shell(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value

      value = exp(-(sqrt(sum((r - centers(:, 1))**2)) - shell_radius)**2)
   end function gaussian_shell

   !> Gaussian times the finite Coulomb potential erf(beta*r)/r
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function gaussian_coulomb(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value
      real(wp) :: radius2, radius, potential

      radius2 = sum((r - centers(:, 2))**2)
      radius = sqrt(radius2)
      if (radius < 1.0e-8_wp) then
         potential = 1.4_wp/sqrt(pi)*(1.0_wp - 0.7_wp**2*radius2/3.0_wp)
      else
         potential = erf(0.7_wp*radius)/radius
      end if
      value = exp(-0.8_wp*radius2)*potential
   end function gaussian_coulomb

   !> Smooth screened Lennard-Jones plus Yukawa field from three centers
   !>
   !> - d = sqrt(r**2 + core**2) removes the hard distance-clamp kink
   !> - Screening both terms removes the algebraic Lennard-Jones tail
   !> - exp(-screening*(d-core)) keeps the screening factor one at the center
   !>
   !> @param[in] r Position (bohr), shape (3)
   pure function smooth_interaction(r) result(value)
      !> Position
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: value
      real(wp) :: distance, sr6
      integer :: i

      value = 0.0_wp
      do i = 1, size(charges)
         distance = sqrt(sum((r - centers(:, i))**2) + core**2)
         sr6 = (sigma/distance)**6
         value = value + exp(-screening*(distance - core)) &
                & *(lj_coefficients(i)*4.0_wp*well_depth*(sr6**2 - sr6) + charges(i)/distance)
      end do
   end function smooth_interaction

   !> Independent radial reference for the smooth interaction field
   !>
   !> 4*pi*integral_0^infinity r**2*f(r) dr, computed at 80 decimal digits
   !> by both tanh-sinh and Gauss-Legendre quadrature; translation leaves
   !> each center's whole-space integral unchanged
   !> Parameters: core=3.5, sigma=5, well_depth=0.01, screening=2
   !> LJ reference includes 4*well_depth*((sigma/d)**12-(sigma/d)**6)
   !> Yukawa reference has unit charge; both include exp(-2*(d-3.5))
   pure function interaction_reference() result(value)
      !> Whole-space integral
      real(wp) :: value
      !> Independent integral of one screened Lennard-Jones center
      real(wp), parameter :: lj_integral = 21.7305618221666_wp
      !> Independent integral of one screened, unit-charge Yukawa center
      real(wp), parameter :: yukawa_integral = 10.953166008972923_wp

      value = sum(lj_coefficients)*lj_integral + sum(charges)*yukawa_integral
   end function interaction_reference

end module test_math_grid_3d_integration
