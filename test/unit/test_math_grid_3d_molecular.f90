!> Test suite for the molecular grid built from atomic recipes
!>
!> - Lifecycle: construction errors, update before construction, private
!>   configuration, independent copies, destroy, accessors
!> - Assembly: translated, multicenter, and rotated integrals, reciprocal
!>   grid on and off, element overrides, per-shell results
!> - Updates: fixed local nodes and weights, repartitioning and pruning
!>   against an independent reconstruction, period guard, Gaussian widths
module test_math_grid_3d_molecular
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_chebyshev2_type, &
      & new_chebyshev2_rule, moist_math_grid_radial_rule_midpoint_type, new_midpoint_rule
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_type, &
      & moist_math_grid_radial_mapping_becke_type, new_becke_mapping, &
      & moist_math_grid_radial_mapping_handymod_type, new_handymod_mapping
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_generator_lebedev_type, &
      & new_lebedev_generator
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_shell_type, &
      & moist_math_grid_atomic_shell_constant_type, new_constant_shell_policy, &
      & moist_math_grid_atomic_recipe_type, moist_math_grid_atomic_recipe_override_type, &
      & default_element_recipes, element_override_index
   use moist_math_grid_atomic_grid, only: moist_math_grid_atomic_type, new_atomic_grid
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_partition, only: becke_partition_weights, ssf_partition_weights, &
      & pvoronoi_partition_weights, &
      & partition_becke, partition_ssf, partition_pvoronoi
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_grid, &
      & molecular_grid_set_kgrid
   implicit none(type, external)
   private

   public :: collect_math_grid_3d_molecular

   !> Rule selectors for `make_recipe`
   integer, parameter :: rule_chebyshev2 = 1, rule_midpoint = 2

   !> Weight-pruning threshold of the grid (bohr^3)
   real(wp), parameter :: wthr = 1.0e-14_wp

contains

   !> Collect all math_grid_3d_molecular tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_grid_3d_molecular(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("constructor_errors", test_constructor_errors), &
         new_unittest("partition_options", test_partition_options), &
         new_unittest("partition_assembly", test_partition_assembly), &
         new_unittest("partition_integrals", test_partition_integrals), &
         new_unittest("update_before_construction", test_update_before_construction), &
         new_unittest("update_geometry_errors", test_update_geometry_errors), &
         new_unittest("private_configuration", test_private_configuration), &
         new_unittest("independent_copies", test_independent_copies), &
         new_unittest("accessors_and_defaults", test_accessors), &
         new_unittest("translated_integral", test_translated_integral), &
         new_unittest("multicenter_integrals", test_multicenter_integrals), &
         new_unittest("rotational_convergence", test_rotational_convergence), &
         new_unittest("reciprocal_on", test_reciprocal_on), &
         new_unittest("reciprocal_off", test_reciprocal_off), &
         new_unittest("element_overrides", test_element_overrides), &
         new_unittest("fixed_local_quadrature", test_fixed_local_quadrature), &
         new_unittest("update_changes_atom_count", test_update_atom_count), &
         new_unittest("period_guard", test_period_guard), &
         new_unittest("gaussian_widths", test_gaussian_widths), &
         new_unittest("destroy_keeps_configuration", test_destroy) &
         ]
   end subroutine collect_math_grid_3d_molecular

   !> Integrate analytic two-center Gaussians with every partition scheme
   !>
   !> @param[out] error Test failure
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
               call becke_recipe(recipe, 100, 0.5_wp, 83, merr)
            else
               call becke_recipe(recipe, 160, 0.5_wp, 131, merr)
            end if
            if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, &
               & partition=scheme, reciprocal=.false., weight_threshold=0.0_wp)
            if (.not. allocated(merr)) call grid%update(mol, merr)
            call require_ok(error, merr, "integration grid")
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

   !> Validate scheme-specific options and exact-zero pruning settings
   !>
   !> @param[out] error Test failure
   subroutine test_partition_options(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_molecular_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: nan

      nan = ieee_value(0.0_wp, ieee_quiet_nan)
      call new_molecular_grid(grid, merr, partition=4)
      call require_error(error, merr, "unknown partition", "scheme identifier above the supported range")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, becke_k=0)
      call require_error(error, merr, "becke_k must be >= 1", "unsmoothed Becke weights")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, partition=-1)
      call require_error(error, merr, "unknown partition", "unknown scheme")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, partition=partition_ssf, becke_k=3)
      call require_error(error, merr, "becke_k requires", "conflicting Becke option")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, partition=partition_pvoronoi, ssf_a=0.64_wp)
      call require_error(error, merr, "ssf_a requires", "conflicting SSF option")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, power_width=1.0_wp)
      call require_error(error, merr, "power_width requires", "power option without power partition")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, partition=partition_pvoronoi, power_width=-1.0_wp)
      call require_error(error, merr, "must be positive", "negative power width")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, partition=partition_pvoronoi, power_width=nan)
      call require_error(error, merr, "must be finite", "NaN power width")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, weight_threshold=-1.0_wp)
      call require_error(error, merr, "must be nonnegative", "negative weight threshold")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, weight_threshold=nan)
      call require_error(error, merr, "must be finite", "NaN weight threshold")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, partition=partition_pvoronoi, power_width=2.5_wp, weight_threshold=0.0_wp)
      call require_ok(error, merr, "power options")
      if (allocated(error)) return
      call check(error, grid%get_partition() == partition_pvoronoi .and. grid%get_power_width() == 2.5_wp &
         & .and. grid%get_weight_threshold() == 0.0_wp, "power options must be stored")
      if (allocated(error)) return
   end subroutine test_partition_options

   !> All schemes select, prune and repartition the molecular grid correctly
   !>
   !> @param[out] error Test failure
   subroutine test_partition_assembly(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_type) :: local(2)
      type(moist_math_grid_3d_molecular_type) :: grid, pruned
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: points(:, :), bw(:)
      real(wp) :: weight
      integer :: scheme, ischeme, step, iat, i, j, total, zeros
      integer, parameter :: schemes(3) = [partition_becke, partition_ssf, partition_pvoronoi]
      integer, parameter :: numbers(3) = [8, 1, 1]

      call becke_recipe(recipe, 14, 0.5_wp, 11, merr, rcut_upper=5.0_wp)
      if (.not. allocated(merr)) call new_atomic_grid(local(1), recipe, 8, merr)
      if (.not. allocated(merr)) call new_atomic_grid(local(2), recipe, 1, merr)
      call require_ok(error, merr, "local grids")
      if (allocated(error)) return
      total = local(1)%npts + 2*local(2)%npts
      do ischeme = 1, size(schemes)
         scheme = schemes(ischeme)
         call make_water(mol)
         call new_molecular_grid(grid, merr, recipe=recipe, reciprocal=.false., partition=scheme, weight_threshold=0.0_wp)
         call require_ok(error, merr, "partition construction")
         if (allocated(error)) return
         do step = 1, 2
            call grid%update(mol, merr)
            call require_ok(error, merr, "partition update")
            if (allocated(error)) return
            j = 0
            zeros = 0
            do iat = 1, 3
               call check(error, grid%atom_offset(iat) == j + 1, "pruned owner offset")
               if (allocated(error)) return
               associate (atom => local(merge(1, 2, iat == 1)))
                  allocate (points(3, atom%npts), bw(atom%npts))
                  do i = 1, atom%npts
                     points(:, i) = atom%shell_r(atom%shell(i))*atom%u(:, i) + mol%xyz(:, iat)
                  end do
                  select case (scheme)
                  case (partition_becke)
                     call becke_partition_weights(iat, points, mol%xyz, numbers, bw)
                  case (partition_ssf)
                     call ssf_partition_weights(iat, points, mol%xyz, numbers, bw)
                  case (partition_pvoronoi)
                     call pvoronoi_partition_weights(iat, points, mol%xyz, numbers, bw)
                  case default
                     call test_failed(error, "unknown partition in assembly test")
                     return
                  end select
                  do i = 1, atom%npts
                     weight = atom%w(i)*bw(i)
                     if (weight == 0.0_wp) then
                        zeros = zeros + 1
                        cycle
                     end if
                     j = j + 1
                     call check(error, j <= grid%ngrid, "nonzero point must be retained")
                     if (allocated(error)) return
                     call check(error, grid%owner(j) == iat .and. all(grid%xyz(:, j) == points(:, i)) &
                        & .and. abs(grid%w(j) - weight) <= 1.0e-14_wp*abs(weight), "partitioned point or weight mismatch")
                     if (allocated(error)) return
                  end do
                  deallocate (points, bw)
               end associate
            end do
            call check(error, j == grid%ngrid .and. j + zeros == total .and. grid%atom_offset(4) == j + 1, &
               & "exact-zero compaction must preserve every nonzero point")
            if (allocated(error)) return
            if (scheme == partition_ssf .or. scheme == partition_pvoronoi) then
               call check(error, zeros > 0, "compact partition must remove redundant points")
               if (allocated(error)) return
            end if
            if (step == 1) mol%xyz(:, 2) = mol%xyz(:, 2) + [7.0_wp, 0.35_wp, -0.2_wp]
         end do
      end do
      ! A nonzero threshold remains an explicit, separate pruning option
      call new_molecular_grid(pruned, merr, recipe=recipe, reciprocal=.false., &
         & partition=partition_pvoronoi, weight_threshold=1.0e-6_wp)
      if (.not. allocated(merr)) call pruned%update(mol, merr)
      call require_ok(error, merr, "threshold pruning")
      if (allocated(error)) return
      call check(error, pruned%ngrid < grid%ngrid, "positive threshold must prune additional nonzero weights")
   end subroutine test_partition_assembly

   !* ================================================================================= *!
   !*                                      Helpers                                      *!
   !* ================================================================================= *!

   !> Forward a library error into a test failure
   !>
   !> @param[out] error  Test failure
   !> @param[in]  merr   Library error, ignored if unallocated
   !> @param[in]  label  Case description
   subroutine require_ok(error, merr, label)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(mctc_error), allocatable, intent(in) :: merr
      !> Case description
      character(len=*), intent(in) :: label

      if (allocated(merr)) call test_failed(error, label//": "//merr%message)
   end subroutine require_ok

   !> Require a library error whose message contains a substring
   !>
   !> @param[out] error     Test failure
   !> @param[in]  merr      Library error
   !> @param[in]  expected  Substring of the expected message
   !> @param[in]  label     Case description
   subroutine require_error(error, merr, expected, label)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(mctc_error), allocatable, intent(in) :: merr
      !> Substring of the expected message
      character(len=*), intent(in) :: expected
      !> Case description
      character(len=*), intent(in) :: label

      if (.not. allocated(merr)) then
         call test_failed(error, label//": expected an error containing '"//expected//"'")
         return
      end if
      if (index(merr%message, expected) == 0) then
         call test_failed(error, label//": message '"//merr%message//"' lacks '"//expected//"'")
      end if
   end subroutine require_error

   !> Build a rule of the requested kind behind the abstract type
   !>
   !> @param[in]     kind    Rule selector
   !> @param[in,out] recipe  Recipe whose radial rule is set
   subroutine set_rule(kind, recipe)
      !> Rule selector
      integer, intent(in) :: kind
      !> Recipe whose radial rule is set
      type(moist_math_grid_atomic_recipe_type), intent(inout) :: recipe

      type(moist_math_grid_radial_rule_chebyshev2_type) :: cheb
      type(moist_math_grid_radial_rule_midpoint_type) :: mid

      if (allocated(recipe%radial%rule)) deallocate (recipe%radial%rule)
      select case (kind)
      case (rule_chebyshev2)
         call new_chebyshev2_rule(cheb)
         allocate (recipe%radial%rule, source=cheb)
      case default
         call new_midpoint_rule(mid)
         allocate (recipe%radial%rule, source=mid)
      end select
   end subroutine set_rule

   !> Assemble an atomic recipe from its parts
   !>
   !> @param[out] recipe      Assembled recipe
   !> @param[in]  kind        Rule selector
   !> @param[in]  npts        Radial node count
   !> @param[in]  mapping     Radial mapping
   !> @param[in]  shells      Shell policy
   !> @param[in]  positive    Lebedev generator admits positive-weight rules only
   !> @param[in]  rcut_lower  Optional lower radial cutoff (bohr)
   !> @param[in]  rcut_upper  Optional upper radial cutoff (bohr)
   subroutine make_recipe(recipe, kind, npts, mapping, shells, positive, rcut_lower, rcut_upper)
      !> Assembled recipe
      type(moist_math_grid_atomic_recipe_type), intent(out) :: recipe
      !> Rule selector
      integer, intent(in) :: kind
      !> Radial node count
      integer, intent(in) :: npts
      !> Radial mapping
      class(moist_math_grid_radial_mapping_type), intent(in) :: mapping
      !> Shell policy
      class(moist_math_grid_atomic_shell_type), intent(in) :: shells
      !> Positive-weight rules only
      logical, intent(in) :: positive
      !> Lower radial cutoff
      real(wp), intent(in), optional :: rcut_lower
      !> Upper radial cutoff
      real(wp), intent(in), optional :: rcut_upper

      type(moist_math_grid_angular_generator_lebedev_type) :: generator

      call set_rule(kind, recipe)
      recipe%radial%npts = npts
      allocate (recipe%radial%mapping, source=mapping)
      if (present(rcut_lower)) recipe%radial%rcut_lower = rcut_lower
      if (present(rcut_upper)) recipe%radial%rcut_upper = rcut_upper
      call new_lebedev_generator(generator, positive_weights_only=positive)
      allocate (recipe%angular, source=generator)
      allocate (recipe%shells, source=shells)
   end subroutine make_recipe

   !> HandyMod recipe with an arbitrary rule and shell policy
   !>
   !> @param[out] recipe    Recipe
   !> @param[in]  kind      Rule selector
   !> @param[in]  nrad      Radial node count
   !> @param[in]  rmin      HandyMod inner radius (bohr)
   !> @param[in]  rmax      HandyMod outer radius (bohr)
   !> @param[in]  m         HandyMod parameter
   !> @param[in]  shells    Shell policy
   !> @param[in]  positive  Positive-weight Lebedev rules only
   !> @param[out] merr      Construction error
   subroutine handymod_recipe(recipe, kind, nrad, rmin, rmax, m, shells, positive, merr)
      !> Recipe
      type(moist_math_grid_atomic_recipe_type), intent(out) :: recipe
      !> Rule selector
      integer, intent(in) :: kind
      !> Radial node count
      integer, intent(in) :: nrad
      !> HandyMod inner radius
      real(wp), intent(in) :: rmin
      !> HandyMod outer radius
      real(wp), intent(in) :: rmax
      !> HandyMod parameter
      real(wp), intent(in) :: m
      !> Shell policy
      class(moist_math_grid_atomic_shell_type), intent(in) :: shells
      !> Positive-weight rules only
      logical, intent(in) :: positive
      !> Construction error
      type(mctc_error), allocatable, intent(out) :: merr

      type(moist_math_grid_radial_mapping_handymod_type) :: handymod

      call new_handymod_mapping(handymod, rmin, rmax, m, merr)
      if (allocated(merr)) return
      call make_recipe(recipe, kind, nrad, handymod, shells, positive)
   end subroutine handymod_recipe

   !> Chebyshev-II x Becke recipe with a constant minimum degree
   !>
   !> @param[out] recipe         Recipe
   !> @param[in]  nrad           Radial node count
   !> @param[in]  radius_factor  Becke scale as a multiple of the covalent radius
   !> @param[in]  degree         Constant minimum Lebedev degree
   !> @param[out] merr           Construction error
   !> @param[in]  rcut_upper     Optional upper radial cutoff (bohr)
   subroutine becke_recipe(recipe, nrad, radius_factor, degree, merr, rcut_upper)
      !> Recipe
      type(moist_math_grid_atomic_recipe_type), intent(out) :: recipe
      !> Radial node count
      integer, intent(in) :: nrad
      !> Becke scale factor
      real(wp), intent(in) :: radius_factor
      !> Constant minimum Lebedev degree
      integer, intent(in) :: degree
      !> Construction error
      type(mctc_error), allocatable, intent(out) :: merr
      !> Upper radial cutoff
      real(wp), intent(in), optional :: rcut_upper

      type(moist_math_grid_atomic_shell_constant_type) :: shells
      type(moist_math_grid_radial_mapping_becke_type) :: becke

      call new_constant_shell_policy(shells, degree, merr)
      if (allocated(merr)) return
      call new_becke_mapping(becke, merr, radius_factor=radius_factor)
      if (allocated(merr)) return
      call make_recipe(recipe, rule_chebyshev2, nrad, becke, shells, .true., rcut_upper=rcut_upper)
   end subroutine becke_recipe

   !> Hydrogen molecule on the x axis, 1.4 bohr bond
   !>
   !> @param[out] mol  Structure
   subroutine make_h2(mol)
      !> Structure
      type(structure_type), intent(out) :: mol

      call new(mol, [1, 1], reshape([-0.7_wp, 0.0_wp, 0.0_wp, 0.7_wp, 0.0_wp, 0.0_wp], [3, 2]))
   end subroutine make_h2

   !> Water, oxygen near the origin
   !>
   !> @param[out] mol  Structure
   subroutine make_water(mol)
      !> Structure
      type(structure_type), intent(out) :: mol

      call new(mol, [8, 1, 1], reshape([0.05_wp, -0.02_wp, 0.12_wp, &
         & 0.0_wp, 1.43_wp, -0.98_wp, &
         & 0.0_wp, -1.43_wp, -0.98_wp], [3, 3]))
   end subroutine make_water

   !> Deterministic small displacement of every nucleus
   !>
   !> @param[in]  mol    Structure
   !> @param[out] moved  Structure with displaced nuclei (at most 0.06 bohr per axis)
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
   !> @param[in] mol  Structure
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
   !> @param[in] grid   Realized grid
   !> @param[in] mol    Structure
   !> @param[in] a      First nucleus
   !> @param[in] b      Second nucleus
   !> @param[in] alpha  Exponent (1/bohr^2)
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
   !> @param[in] mol    Structure
   !> @param[in] a      First nucleus
   !> @param[in] b      Second nucleus
   !> @param[in] alpha  Exponent (1/bohr^2)
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

   !* ================================================================================= *!
   !*                                     Lifecycle                                     *!
   !* ================================================================================= *!

   !> Invalid options and incomplete recipes fail construction; the grid stays unconstructed
   !>
   !> @param[out] error  Test failure
   subroutine test_constructor_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_atomic_recipe_type) :: recipe, broken
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      real(wp) :: nan

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      call make_h2(mol)
      call becke_recipe(recipe, 8, 0.5_wp, 5, merr)
      call require_ok(error, merr, "recipe")
      if (allocated(error)) return

      call new_molecular_grid(grid, merr, becke_k=-1)
      call require_error(error, merr, "becke_k must be >= 1", "becke_k")
      if (allocated(error)) return
      call grid%update(mol, merr)
      call require_error(error, merr, "new_molecular_grid", "update after a failed construction")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, ssf_a=0.0_wp)
      call require_error(error, merr, "ssf_a must be in (0, 1]", "ssf_a = 0")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, ssf_a=1.5_wp)
      call require_error(error, merr, "ssf_a must be in (0, 1]", "ssf_a = 1.5")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, ssf_a=nan)
      call require_error(error, merr, "ssf_a must be finite", "ssf_a = NaN")
      if (allocated(error)) return
      ! SSF replaces the Becke polynomial, so its iteration count would be ignored
      call new_molecular_grid(grid, merr, becke_k=3, ssf_a=0.64_wp)
      call require_error(error, merr, "either becke_k or ssf_a", "becke_k with ssf_a")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, pruning_threshold=1.0_wp)
      call require_error(error, merr, "pruning_threshold must be in [0, 1)", "pruning_threshold = 1")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, pruning_threshold=-0.1_wp)
      call require_error(error, merr, "pruning_threshold must be in [0, 1)", "pruning_threshold < 0")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, dr=nan)
      call require_error(error, merr, "settings must be finite", "dr = NaN")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, dr=0.0_wp)
      call require_error(error, merr, "invalid dr, kbuffer or nufft_tol", "dr = 0")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, kbuffer=-1.0_wp)
      call require_error(error, merr, "invalid dr, kbuffer or nufft_tol", "kbuffer < 0")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, nufft_tol=0.0_wp)
      call require_error(error, merr, "invalid dr, kbuffer or nufft_tol", "nufft_tol = 0")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, xi0_factor=0.0_wp)
      call require_error(error, merr, "xi0_factor must be positive", "xi0_factor = 0")
      if (allocated(error)) return

      broken = recipe
      deallocate (broken%shells)
      call new_molecular_grid(grid, merr, recipe=broken)
      call require_error(error, merr, "needs a radial rule", "recipe without shell policy")
      if (allocated(error)) return
      broken = recipe
      broken%radial%npts = 0
      call new_molecular_grid(grid, merr, recipe=broken)
      call require_error(error, merr, "needs at least one radial node", "npts = 0")
      if (allocated(error)) return
      broken = recipe
      broken%radial%rcut_upper = nan
      call new_molecular_grid(grid, merr, recipe=broken)
      call require_error(error, merr, "settings must be finite", "rcut_upper = NaN")
      if (allocated(error)) return
      broken%radial%rcut_upper = ieee_value(1.0_wp, ieee_positive_inf)
      call new_molecular_grid(grid, merr, recipe=broken)
      call require_error(error, merr, "settings must be finite", "rcut_upper = Inf")
      if (allocated(error)) return
      broken = recipe
      broken%radial%rcut_lower = 3.0_wp
      broken%radial%rcut_upper = 2.0_wp
      call new_molecular_grid(grid, merr, recipe=broken)
      call require_error(error, merr, "rcut_lower < rcut_upper", "crossed cutoffs")
      if (allocated(error)) return
      broken = recipe
      select type (shells => broken%shells)
      type is (moist_math_grid_atomic_shell_constant_type)
         shells%degree = -1
      end select
      call new_molecular_grid(grid, merr, recipe=broken)
      call require_error(error, merr, "Constant shell policy", "invalid shell policy")
      if (allocated(error)) return

      allocate (overrides(1))
      overrides(1)%recipe = recipe
      call new_molecular_grid(grid, merr, recipe=recipe, overrides=overrides)
      call require_error(error, merr, "override 1 lists no elements", "override without elements")
      if (allocated(error)) return
      overrides(1)%elements = [1]
      deallocate (overrides(1)%recipe%angular)
      call new_molecular_grid(grid, merr, recipe=recipe, overrides=overrides)
      call require_error(error, merr, "override 1 needs a radial rule", "incomplete override recipe")
      if (allocated(error)) return

      ! A generator request that no rule satisfies fails at the first update
      broken = recipe
      select type (shells => broken%shells)
      type is (moist_math_grid_atomic_shell_constant_type)
         shells%degree = 200
      end select
      call new_molecular_grid(grid, merr, recipe=broken)
      call require_ok(error, merr, "unreachable degree construction")
      if (allocated(error)) return
      call grid%update(mol, merr)
      call require_error(error, merr, "molecular domain: Atomic grid for z = 1", "unreachable degree")
   end subroutine test_constructor_errors

   !> A never-constructed grid refuses update, validate, and rebuild
   !>
   !> @param[out] error  Test failure
   subroutine test_update_before_construction(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_molecular_type) :: grid
      class(moist_math_grid_3d_type), allocatable :: copy
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr

      call make_h2(mol)
      call grid%update(mol, merr)
      call require_error(error, merr, "construct the grid with new_molecular_grid before update", "update")
      if (allocated(error)) return
      call check(error, grid%ngrid == 0 .and. .not. allocated(grid%xyz), "failed update left points")
      if (allocated(error)) return
      call grid%validate(merr)
      call require_error(error, merr, "new_molecular_grid", "validate")
      if (allocated(error)) return
      call grid%rebuild(merr)
      call require_error(error, merr, "update before rebuild", "rebuild")
      if (allocated(error)) return
      allocate (copy, source=grid)
      call copy%update(mol, merr)
      call require_error(error, merr, "new_molecular_grid", "update of a copy")
   end subroutine test_update_before_construction

   !> Geometry errors of a constructed grid
   !>
   !> @param[out] error  Test failure
   subroutine test_update_geometry_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol, empty
      type(mctc_error), allocatable :: merr

      call new_molecular_grid(grid, merr)
      call require_ok(error, merr, "construction")
      if (allocated(error)) return
      call grid%update(empty, merr)
      call require_error(error, merr, "at least one solute atom", "no atoms")
      if (allocated(error)) return
      call make_h2(mol)
      mol%xyz(2, 1) = ieee_value(1.0_wp, ieee_quiet_nan)
      call grid%update(mol, merr)
      call require_error(error, merr, "coordinates must be finite", "NaN coordinate")
      if (allocated(error)) return
      ! A dummy atom has no covalent radius for the partition
      call make_h2(mol)
      mol%num(1) = 0
      call grid%update(mol, merr)
      call require_error(error, merr, "atomic numbers must be between", "dummy atom")
      if (allocated(error)) return
      call new_molecular_grid(grid, merr, dr=1.0e-12_wp)
      call require_ok(error, merr, "tiny dr construction")
      if (allocated(error)) return
      call make_h2(mol)
      call grid%update(mol, merr)
      call require_error(error, merr, "reciprocal grid is too large", "tiny dr")
   end subroutine test_update_geometry_errors

   !> Later changes to the caller's recipe objects do not reach the grid
   !>
   !> The grid built from mutated sources equals a grid built from pristine
   !> copies, point for point
   !>
   !> @param[out] error  Test failure
   subroutine test_private_configuration(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe, pristine
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:), pristine_overrides(:)
      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(moist_math_grid_atomic_shell_constant_type) :: shells
      type(moist_math_grid_3d_molecular_type) :: grid, reference
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr

      call make_water(mol)
      call becke_recipe(recipe, 12, 0.5_wp, 7, merr)
      call require_ok(error, merr, "recipe")
      if (allocated(error)) return
      allocate (overrides(1))
      overrides(1)%elements = [1]
      call becke_recipe(overrides(1)%recipe, 9, 1.0_wp, 5, merr)
      call require_ok(error, merr, "override recipe")
      if (allocated(error)) return
      pristine = recipe
      pristine_overrides = overrides

      call new_molecular_grid(grid, merr, recipe=recipe, overrides=overrides, becke_k=2)
      call require_ok(error, merr, "construction")
      if (allocated(error)) return

      ! Mutate every level of the sources after construction
      recipe%radial%npts = 30
      deallocate (recipe%radial%mapping)
      call new_becke_mapping(becke, merr, radius_factor=2.0_wp)
      allocate (recipe%radial%mapping, source=becke)
      call new_constant_shell_policy(shells, 29, merr)
      deallocate (recipe%shells)
      allocate (recipe%shells, source=shells)
      overrides(1)%elements = [8]
      overrides(1)%recipe%radial%npts = 3
      select type (policy => overrides(1)%recipe%shells)
      type is (moist_math_grid_atomic_shell_constant_type)
         policy%degree = 41
      end select
      deallocate (overrides)

      call grid%update(mol, merr)
      call require_ok(error, merr, "update")
      if (allocated(error)) return
      call new_molecular_grid(reference, merr, recipe=pristine, overrides=pristine_overrides, becke_k=2)
      if (.not. allocated(merr)) call reference%update(mol, merr)
      call require_ok(error, merr, "reference")
      if (allocated(error)) return
      call check(error, grid%ngrid == reference%ngrid .and. all(grid%nrad_per_atom == [12, 9, 9]) &
         & .and. all(grid%nang_per_atom == [26, 14, 14]), "grid followed changes to its source recipes")
      if (allocated(error)) return
      call check(error, all(grid%xyz == reference%xyz) .and. all(grid%w == reference%w), &
         & "grid differs from one built from pristine recipes")
   end subroutine test_private_configuration

   !> `allocate(copy, source=grid)` gives an independent grid, before and after the first update
   !>
   !> @param[out] error  Test failure
   subroutine test_independent_copies(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: template, grid, fresh
      class(moist_math_grid_3d_type), allocatable :: before, after
      type(structure_type) :: mol, moved
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: xyz(:, :), w(:)

      call make_water(mol)
      moved = mol
      moved%xyz = mol%xyz + spread([0.3_wp, -0.2_wp, 0.1_wp], 2, mol%nat)
      call becke_recipe(recipe, 10, 0.5_wp, 7, merr, rcut_upper=5.0_wp)
      if (.not. allocated(merr)) call new_molecular_grid(template, merr, recipe=recipe, dr=0.8_wp)
      call require_ok(error, merr, "construction")
      if (allocated(error)) return

      ! Copy of the configured template, as MOZ 3D takes it
      allocate (before, source=template)
      grid = template
      call grid%update(mol, merr)
      if (.not. allocated(merr)) call before%update(moved, merr)
      call require_ok(error, merr, "updates of template copies")
      if (allocated(error)) return
      call check(error, template%ngrid == 0 .and. .not. allocated(template%xyz), &
         & "updating a copy realized the template")
      if (allocated(error)) return
      call new_molecular_grid(fresh, merr, recipe=recipe, dr=0.8_wp)
      if (.not. allocated(merr)) call fresh%update(moved, merr)
      call require_ok(error, merr, "fresh grid")
      if (allocated(error)) return
      call check(error, before%ngrid == fresh%ngrid, "template copy differs from a fresh grid")
      if (allocated(error)) return
      call check(error, all(before%xyz == fresh%xyz) .and. all(before%w == fresh%w), &
         & "template copy differs from a fresh grid")
      if (allocated(error)) return

      ! Copy of a realized grid: the source moves, the copy keeps its points
      xyz = grid%xyz
      w = grid%w
      allocate (after, source=grid)
      call grid%update(moved, merr)
      call require_ok(error, merr, "update of the source")
      if (allocated(error)) return
      call check(error, all(after%xyz == xyz) .and. all(after%w == w), "copy followed its source")
      if (allocated(error)) return
      call check(error, all(grid%xyz == fresh%xyz) .and. all(grid%w == fresh%w), &
         & "moved source differs from a fresh grid")
      if (allocated(error)) return
      select type (after)
      type is (moist_math_grid_3d_molecular_type)
         call check(error, all(after%kref == centroid(mol)), "copy kref followed its source")
         if (allocated(error)) return
         call after%update(moved, merr)
         call require_ok(error, merr, "update of the copy")
         if (allocated(error)) return
         call check(error, all(after%xyz == fresh%xyz) .and. all(after%w == fresh%w) &
            & .and. after%nkx == grid%nkx .and. all(after%kref == grid%kref), &
            & "moved copy differs from its moved source")
      class default
         call test_failed(error, "copy lost its dynamic type")
      end select
   end subroutine test_independent_copies

   !> Accessors return the constructor values and the documented defaults
   !>
   !> @param[out] error  Test failure
   subroutine test_accessors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr

      call new_molecular_grid(grid, merr)
      call require_ok(error, merr, "defaults")
      if (allocated(error)) return
      call check(error, grid%get_dr() == 0.5_wp .and. grid%get_kbuffer() == 2.0_wp &
         & .and. grid%get_xi0_factor() == 1.0_wp .and. grid%get_becke_k() == 3 &
         & .and. grid%get_pruning_threshold() == 0.0_wp .and. grid%get_nufft_tol() == 1.0e-10_wp &
         & .and. .not. grid%has_geometry_dependent_xi0(), "default settings")
      if (allocated(error)) return
      call make_h2(mol)
      call grid%update(mol, merr)
      call require_ok(error, merr, "default update")
      if (allocated(error)) return
      ! Default recipe: 50 midpoint shells, 110-point rule, reciprocal grid on
      call check(error, all(grid%nrad_per_atom == 50) .and. all(grid%nang_per_atom == 110) &
         & .and. all(grid%shell_npts == 110) .and. all(grid%shell_degree == 17) &
         & .and. grid%has_kgrid .and. maxval(grid%shell_r) < 10.0_wp, "default recipe")
      if (allocated(error)) return

      call new_molecular_grid(grid, merr, becke_k=5, pruning_threshold=0.25_wp, dr=0.7_wp, &
         & kbuffer=1.5_wp, nufft_tol=1.0e-7_wp, gaussian=.true., xi0_factor=1.3_wp)
      call require_ok(error, merr, "options")
      if (allocated(error)) return
      call check(error, grid%get_dr() == 0.7_wp .and. grid%get_kbuffer() == 1.5_wp &
         & .and. grid%get_xi0_factor() == 1.3_wp .and. grid%get_becke_k() == 5 &
         & .and. grid%get_pruning_threshold() == 0.25_wp .and. grid%get_nufft_tol() == 1.0e-7_wp &
         & .and. grid%has_geometry_dependent_xi0(), "constructor settings")
   end subroutine test_accessors

   !* ================================================================================= *!
   !*                                     Assembly                                      *!
   !* ================================================================================= *!

   !> A translated atom keeps its weights and moves its points; its Gaussian integral stays exact
   !>
   !> @param[out] error  Test failure
   subroutine test_translated_integral(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: alpha = 0.7_wp
      real(wp), parameter :: shift(3) = [0.37_wp, -1.21_wp, 2.05_wp]
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: xyz(:, :), w(:)
      real(wp) :: exact, before, after
      character(len=120) :: msg

      call new(mol, [6], reshape([0.1_wp, 0.2_wp, -0.3_wp], [3, 1]))
      call becke_recipe(recipe, 60, 0.5_wp, 29, merr)
      if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, reciprocal=.false.)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      call require_ok(error, merr, "carbon atom")
      if (allocated(error)) return
      exact = gauss_pair_exact(mol, 1, 1, alpha)
      before = gauss_pair_integral(grid, mol, 1, 1, alpha)
      xyz = grid%xyz
      w = grid%w

      mol%xyz(:, 1) = mol%xyz(:, 1) + shift
      call grid%update(mol, merr)
      call require_ok(error, merr, "translated atom")
      if (allocated(error)) return
      after = gauss_pair_integral(grid, mol, 1, 1, alpha)
      write (msg, "(a,es10.3,a,es10.3)") "Gaussian integral errors ", abs(before - exact)/exact, &
         & " and ", abs(after - exact)/exact
      call check(error, abs(before - exact) <= 1.0e-10_wp*exact .and. abs(after - exact) <= 1.0e-10_wp*exact, &
         & trim(msg))
      if (allocated(error)) return
      call check(error, size(grid%w) == size(w), "translation changed the point count")
      if (allocated(error)) return
      ! One atom: partition weight exactly 1, local weights fixed
      call check(error, all(grid%w == w), "translation changed single-atom weights")
      if (allocated(error)) return
      call check(error, maxval(abs(grid%xyz - xyz - spread(shift, 2, grid%ngrid))) &
         & <= 8.0_wp*epsilon(1.0_wp)*(1.0_wp + maxval(abs(xyz))), "points did not move with the atom")
   end subroutine test_translated_integral

   !> One- and two-center Gaussian integrals over a water grid with the per-element defaults
   !>
   !> @param[out] error  Test failure
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
      if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, overrides=overrides, &
         & reciprocal=.false.)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      call require_ok(error, merr, "water")
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
   !>
   !> @param[out] error  Test failure
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
         call becke_recipe(recipe, 40, 0.5_wp, degrees(id), merr)
         if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, reciprocal=.false.)
         if (.not. allocated(merr)) call grid%update(mol, merr)
         call require_ok(error, merr, "reference orientation")
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
            call require_ok(error, merr, "rotated orientation")
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

   !> `reciprocal = .true.` sizes the k-grid about the centroid on the first update
   !>
   !> - Modes and spacings as `molecular_grid_set_kgrid` sizes them for the
   !>   constructor's `dr` and `kbuffer`
   !> - A translating update keeps the period and recenters `kref`
   !>
   !> @param[out] error  Test failure
   subroutine test_reciprocal_on(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid, manual
      type(structure_type) :: mol, moved
      type(mctc_error), allocatable :: merr
      integer :: modes(3)

      call make_water(mol)
      call becke_recipe(recipe, 12, 0.5_wp, 7, merr, rcut_upper=5.0_wp)
      if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, dr=0.6_wp, kbuffer=1.5_wp)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      call require_ok(error, merr, "reciprocal grid")
      if (allocated(error)) return
      call check(error, grid%has_kgrid .and. grid%npts_k == grid%nkx*grid%nky*grid%nkz &
         & .and. grid%npts_k > 0 .and. all(grid%kref == centroid(mol)), "auto-sized reciprocal grid")
      if (allocated(error)) return

      call new_molecular_grid(manual, merr, recipe=recipe, reciprocal=.false.)
      if (.not. allocated(merr)) call manual%update(mol, merr)
      if (.not. allocated(merr)) call molecular_grid_set_kgrid(manual, 0.6_wp, merr, buffer=1.5_wp, &
         & reference=centroid(mol))
      call require_ok(error, merr, "manual reciprocal grid")
      if (allocated(error)) return
      call check(error, grid%nkx == manual%nkx .and. grid%nky == manual%nky .and. grid%nkz == manual%nkz &
         & .and. grid%dkx == manual%dkx .and. grid%dky == manual%dky .and. grid%dkz == manual%dkz, &
         & "auto-sized modes differ from molecular_grid_set_kgrid")
      if (allocated(error)) return

      modes = [grid%nkx, grid%nky, grid%nkz]
      moved = mol
      moved%xyz = mol%xyz + spread([0.4_wp, 0.1_wp, -0.3_wp], 2, mol%nat)
      call grid%update(moved, merr)
      call require_ok(error, merr, "translated")
      if (allocated(error)) return
      call check(error, all([grid%nkx, grid%nky, grid%nkz] == modes) .and. all(grid%kref == centroid(moved)), &
         & "update must keep the period and recenter kref")
   end subroutine test_reciprocal_on

   !> `reciprocal = .false.` makes no k-grid until `molecular_grid_set_kgrid`, which stores its settings
   !>
   !> @param[out] error  Test failure
   subroutine test_reciprocal_off(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid, fresh
      type(structure_type) :: mol, moved
      type(mctc_error), allocatable :: merr
      integer :: modes(3)

      call make_water(mol)
      moved = mol
      moved%xyz = mol%xyz + spread([-0.2_wp, 0.5_wp, 0.3_wp], 2, mol%nat)
      call becke_recipe(recipe, 12, 0.5_wp, 7, merr, rcut_upper=5.0_wp)
      if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, reciprocal=.false.)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      call require_ok(error, merr, "integration-only grid")
      if (allocated(error)) return
      call check(error, .not. grid%has_kgrid .and. grid%npts_k == 0 .and. grid%nkx == 0, &
         & "reciprocal = .false. made a k-grid")
      if (allocated(error)) return
      call grid%update(moved, merr)
      if (.not. allocated(merr)) call grid%rebuild(merr)
      call require_ok(error, merr, "update and rebuild without k-grid")
      if (allocated(error)) return
      call check(error, .not. grid%has_kgrid .and. grid%npts_k == 0, "update or rebuild made a k-grid")
      if (allocated(error)) return

      call molecular_grid_set_kgrid(grid, 0.7_wp, merr, buffer=1.25_wp, reference=centroid(moved))
      call require_ok(error, merr, "set_kgrid")
      if (allocated(error)) return
      call check(error, grid%has_kgrid .and. grid%get_dr() == 0.7_wp .and. grid%get_kbuffer() == 1.25_wp, &
         & "set_kgrid must store dr and buffer")
      if (allocated(error)) return
      modes = [grid%nkx, grid%nky, grid%nkz]
      call grid%update(mol, merr)
      call require_ok(error, merr, "update with a set k-grid")
      if (allocated(error)) return
      call check(error, all([grid%nkx, grid%nky, grid%nkz] == modes) .and. all(grid%kref == centroid(mol)), &
         & "update must keep the set period about the centroid")
      if (allocated(error)) return
      ! Rebuild re-sizes with the stored values, as a fresh reciprocal grid would
      call grid%rebuild(merr)
      if (.not. allocated(merr)) call new_molecular_grid(fresh, merr, recipe=recipe, dr=0.7_wp, kbuffer=1.25_wp)
      if (.not. allocated(merr)) call fresh%update(mol, merr)
      call require_ok(error, merr, "rebuild and fresh grid")
      if (allocated(error)) return
      call check(error, grid%nkx == fresh%nkx .and. grid%nky == fresh%nky .and. grid%nkz == fresh%nkz &
         & .and. all(grid%kref == fresh%kref) .and. grid%dkx == fresh%dkx, &
         & "rebuild must re-size with the stored reciprocal settings")
   end subroutine test_reciprocal_off

   !> Per-element overrides select the recipe per atom; per-shell results follow it
   !>
   !> - First override listing an element wins
   !> - Requested shells reported before the cutoff, retained shells in the CSR arrays
   !>
   !> @param[out] error  Test failure
   subroutine test_element_overrides(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      type(moist_math_grid_atomic_type) :: local
      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      integer :: iat, idx, lo, hi, z

      call make_water(mol)
      call becke_recipe(recipe, 20, 0.5_wp, 11, merr)
      allocate (overrides(2))
      overrides(1)%elements = [2, 1]
      if (.not. allocated(merr)) call becke_recipe(overrides(1)%recipe, 16, 1.0_wp, 5, merr, rcut_upper=4.0_wp)
      overrides(2)%elements = [1, 8]
      if (.not. allocated(merr)) call becke_recipe(overrides(2)%recipe, 30, 0.5_wp, 29, merr)
      if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, overrides=overrides, &
         & reciprocal=.false.)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      call require_ok(error, merr, "overrides")
      if (allocated(error)) return
      call check(error, all(grid%nrad_per_atom == [30, 16, 16]) .and. all(grid%nang_per_atom == [302, 14, 14]), &
         & "per-atom sizes do not follow the first matching override")
      if (allocated(error)) return

      do iat = 1, mol%nat
         z = mol%num(mol%id(iat))
         idx = element_override_index(z, overrides)
         call new_atomic_grid(local, overrides(idx)%recipe, z, merr)
         call require_ok(error, merr, "atomic grid")
         if (allocated(error)) return
         lo = grid%atom_shell_offset(iat)
         hi = grid%atom_shell_offset(iat + 1) - 1
         call check(error, hi - lo + 1 == local%nshell .and. all(grid%shell_r(lo:hi) == local%shell_r) &
            & .and. all(grid%shell_npts(lo:hi) == local%shell_npts) &
            & .and. all(grid%shell_degree(lo:hi) == local%shell_degree) &
            & .and. all(grid%shell_spacing(lo:hi) == local%shell_spacing), "per-shell results")
         if (allocated(error)) return
         call check(error, grid%atom_offset(iat + 1) - grid%atom_offset(iat) <= local%npts &
            & .and. all(grid%owner(grid%atom_offset(iat):grid%atom_offset(iat + 1) - 1) == iat), &
            & "atom points")
         if (allocated(error)) return
      end do
      ! The hydrogen cutoff keeps fewer shells than requested
      call check(error, grid%atom_shell_offset(3) - grid%atom_shell_offset(2) < grid%nrad_per_atom(2) &
         & .and. maxval(grid%shell_r(grid%atom_shell_offset(2):grid%atom_shell_offset(3) - 1)) <= 4.0_wp, &
         & "rcut_upper of the hydrogen override")
   end subroutine test_element_overrides

   !* ================================================================================= *!
   !*                                  Geometry updates                                 *!
   !* ================================================================================= *!

   !> Updates translate fixed local grids, repartition, and prune, as an independent reconstruction
   !>
   !> Reference per atom: the element's atomic grid, points `r*u + R`,
   !> batched partition weights, weight `w_local*partition`, kept if
   !> `|w| >= 1e-14` and the partition weight reaches the threshold; at two
   !> geometries, the second crossing the threshold
   !>
   !> @param[out] error  Test failure
   subroutine test_fixed_local_quadrature(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: threshold = 0.1_wp
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_type) :: local(2)
      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: points(:, :), bw(:)
      real(wp) :: x(3), w
      integer :: step, iat, i, j, k, count_before, numbers(3)
      character(len=120) :: msg

      call make_water(mol)
      numbers = [8, 1, 1]
      call becke_recipe(recipe, 14, 0.5_wp, 11, merr)
      if (.not. allocated(merr)) call new_atomic_grid(local(1), recipe, 8, merr)
      if (.not. allocated(merr)) call new_atomic_grid(local(2), recipe, 1, merr)
      if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, becke_k=2, &
         & pruning_threshold=threshold, reciprocal=.false.)
      call require_ok(error, merr, "construction")
      if (allocated(error)) return

      count_before = 0
      do step = 1, 2
         call grid%update(mol, merr)
         call require_ok(error, merr, "update")
         if (allocated(error)) return
         j = 0
         do iat = 1, mol%nat
            associate (atom => local(merge(1, 2, iat == 1)))
               allocate (points(3, atom%npts), bw(atom%npts))
               do i = 1, atom%npts
                  do k = 1, 3
                     points(k, i) = atom%shell_r(atom%shell(i))*atom%u(k, i) + mol%xyz(k, iat)
                  end do
               end do
               call becke_partition_weights(iat, points, mol%xyz, numbers, bw, stiffness=2)
               do i = 1, atom%npts
                  w = atom%w(i)*bw(i)
                  if (abs(w) < wthr .or. bw(i) < threshold) cycle
                  j = j + 1
                  if (j > grid%ngrid) exit
                  x = points(:, i)
                  write (msg, "(a,i0,a,i0,a,i0)") "step ", step, ", atom ", iat, ", local point ", i
                  call check(error, grid%owner(j) == iat .and. &
                     & maxval(abs(grid%xyz(:, j) - x)) <= 4.0_wp*epsilon(1.0_wp)*(1.0_wp + maxval(abs(x))) &
                     & .and. abs(grid%w(j) - w) <= 1.0e-14_wp*abs(w), trim(msg))
                  if (allocated(error)) return
               end do
               deallocate (points, bw)
            end associate
         end do
         call check(error, j == grid%ngrid, "retained points differ from the reconstruction")
         if (allocated(error)) return
         if (step == 1) then
            count_before = grid%ngrid
            mol%xyz(:, 2) = mol%xyz(:, 2) + [0.0_wp, 0.35_wp, -0.2_wp]
         else
            call check(error, grid%ngrid /= count_before, "deformation must cross the pruning threshold")
            if (allocated(error)) return
         end if
      end do
   end subroutine test_fixed_local_quadrature

   !> An update to a molecule with a different atom count equals a fresh grid
   !>
   !> H2 to water to H2 on one grid, each against a grid constructed and
   !> updated once; the same code builds both, so the results must be equal
   !>
   !> @param[out] error  Test failure
   subroutine test_update_atom_count(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid, fresh
      type(structure_type) :: mols(3)
      type(mctc_error), allocatable :: merr
      character(len=32) :: label
      integer :: step

      call make_h2(mols(1))
      call make_water(mols(2))
      call make_h2(mols(3))
      call becke_recipe(recipe, 12, 0.5_wp, 7, merr, rcut_upper=5.0_wp)
      if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, reciprocal=.false.)
      call require_ok(error, merr, "construction")
      if (allocated(error)) return

      do step = 1, size(mols)
         write (label, "(a,i0)") "step ", step
         call grid%update(mols(step), merr)
         if (.not. allocated(merr)) call new_molecular_grid(fresh, merr, recipe=recipe, reciprocal=.false.)
         if (.not. allocated(merr)) call fresh%update(mols(step), merr)
         call require_ok(error, merr, trim(label))
         if (allocated(error)) return
         call check(error, grid%ngrid == fresh%ngrid .and. grid%natom == mols(step)%nat &
            & .and. size(grid%atom_offset) == mols(step)%nat + 1, trim(label)//": sizes")
         if (allocated(error)) return
         call check(error, all(grid%atom_offset == fresh%atom_offset) .and. all(grid%owner == fresh%owner) &
            & .and. all(grid%nrad_per_atom == fresh%nrad_per_atom) &
            & .and. all(grid%atom_shell_offset == fresh%atom_shell_offset), trim(label)//": ownership")
         if (allocated(error)) return
         call check(error, all(grid%xyz == fresh%xyz) .and. all(grid%w == fresh%w), &
            & trim(label)//": points and weights")
         if (allocated(error)) return
      end do
   end subroutine test_update_atom_count

   !> Translation keeps modes, weights, and owners; a large stretch needs rebuild
   !>
   !> @param[out] error  Test failure
   subroutine test_period_guard(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_shell_constant_type) :: shells
      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: xyz(:, :), w(:)
      integer, allocatable :: owner(:)
      real(wp) :: shift(3)
      integer :: modes(3)

      call make_h2(mol)
      call new_constant_shell_policy(shells, 17, merr)
      if (.not. allocated(merr)) call handymod_recipe(recipe, rule_midpoint, 12, 0.0_wp, 4.0_wp, 2.0_wp, &
         & shells, .true., merr)
      if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, dr=1.0_wp, kbuffer=1.0_wp)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      call require_ok(error, merr, "construction")
      if (allocated(error)) return
      xyz = grid%xyz
      w = grid%w
      owner = grid%owner
      modes = [grid%nkx, grid%nky, grid%nkz]

      shift = [0.125_wp, -0.25_wp, 0.5_wp]
      mol%xyz = mol%xyz + spread(shift, 2, mol%nat)
      call grid%update(mol, merr)
      call require_ok(error, merr, "translation")
      if (allocated(error)) return
      call check(error, grid%ngrid == size(w) .and. all(grid%owner == owner) &
         & .and. maxval(abs(grid%xyz - xyz - spread(shift, 2, grid%ngrid))) < 1.6e-14_wp &
         & .and. maxval(abs(grid%w - w)) < 1.0e-13_wp &
         & .and. all([grid%nkx, grid%nky, grid%nkz] == modes), "translation")
      if (allocated(error)) return

      xyz = grid%xyz
      mol%xyz(1, 2) = mol%xyz(1, 2) + 8.0_wp
      call grid%update(mol, merr)
      call require_error(error, merr, "rebuild", "large stretch")
      if (allocated(error)) return
      call check(error, all(grid%xyz == xyz), "failed update changed the committed points")
      if (allocated(error)) return
      call grid%rebuild(merr)
      call require_ok(error, merr, "rebuild")
      if (allocated(error)) return
      call check(error, grid%nkx > modes(1) .and. all(grid%kref == centroid(mol)), "rebuild must enlarge the period")
      if (allocated(error)) return
      call grid%update(mol, merr)
      call require_ok(error, merr, "update after rebuild")
   end subroutine test_period_guard

   !> Gaussian widths follow the weights at every geometry; point potentials publish none
   !>
   !> @param[out] error  Test failure
   subroutine test_gaussian_widths(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol, moved
      type(mctc_error), allocatable :: merr
      integer :: step

      call make_water(mol)
      call displace(mol, moved)
      call becke_recipe(recipe, 10, 0.5_wp, 7, merr, rcut_upper=5.0_wp)
      if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, gaussian=.true., &
         & xi0_factor=1.2_wp, ssf_a=0.64_wp, reciprocal=.false.)
      call require_ok(error, merr, "construction")
      if (allocated(error)) return
      do step = 1, 2
         if (step == 1) then
            call grid%update(mol, merr)
         else
            call grid%update(moved, merr)
         end if
         call require_ok(error, merr, "update")
         if (allocated(error)) return
         call check(error, allocated(grid%xi0) .and. grid%has_geometry_dependent_xi0(), "widths")
         if (allocated(error)) return
         call check(error, size(grid%xi0) == grid%ngrid .and. &
            & maxval(abs(grid%xi0**3*grid%w - 1.2_wp**3)) < 1.0e-13_wp, "width scale law")
         if (allocated(error)) return
      end do

      call new_molecular_grid(grid, merr, recipe=recipe)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      call require_ok(error, merr, "point potentials")
      if (allocated(error)) return
      call check(error, .not. allocated(grid%xi0) .and. .not. grid%has_geometry_dependent_xi0(), &
         & "point potentials must publish no widths")
   end subroutine test_gaussian_widths

   !> destroy frees the geometry and keeps the configuration; the next update equals a fresh grid
   !>
   !> @param[out] error  Test failure
   subroutine test_destroy(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid, fresh
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr

      call make_h2(mol)
      call becke_recipe(recipe, 8, 0.5_wp, 7, merr, rcut_upper=4.0_wp)
      if (.not. allocated(merr)) call new_molecular_grid(grid, merr, recipe=recipe, dr=1.0_wp)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      if (.not. allocated(merr)) call molecular_grid_set_kgrid(grid, grid%get_dr(), merr, &
         & reference=centroid(mol) + [0.5_wp, 0.2_wp, 0.0_wp])
      call require_ok(error, merr, "grid with an offset reference")
      if (allocated(error)) return
      call grid%destroy()
      call check(error, grid%ngrid == 0 .and. .not. allocated(grid%xyz) .and. .not. allocated(grid%atom_offset) &
         & .and. .not. allocated(grid%shell_r) .and. .not. grid%has_kgrid, "destroy left geometry behind")
      if (allocated(error)) return
      call grid%destroy()
      call grid%update(mol, merr)
      if (.not. allocated(merr)) call new_molecular_grid(fresh, merr, recipe=recipe, dr=1.0_wp)
      if (.not. allocated(merr)) call fresh%update(mol, merr)
      call require_ok(error, merr, "update after destroy")
      if (allocated(error)) return
      call check(error, grid%ngrid == fresh%ngrid .and. all(grid%xyz == fresh%xyz) .and. all(grid%w == fresh%w) &
         & .and. all(grid%kref == centroid(mol)) .and. grid%nkx == fresh%nkx, &
         & "update after destroy differs from a fresh grid")
   end subroutine test_destroy

end module test_math_grid_3d_molecular
