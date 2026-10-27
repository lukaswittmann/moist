!> Tests of constrained EC fitting, coupling requirements and fit adjoints
module test_moz_potential_ec
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type, fatal_error
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_channels_coupling, only: coupling_type, coupling_view_type, coupling_request_type, &
      & coupling_begin_registration, coupling_set_scope, coupling_snapshot, coupling_arm, coupling_invalidate, &
      & coupling_make_view, coupling_close_view, moist_phase_energy, moist_phase_response, moist_phase_gradient
   use moist_channels_response, only: response_type, response_item_type, potential_adjoint_response_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_model_moz_potential, only: moz_potential_type, ec_charges_type, custom_lj_type, new_custom_lj
   use test_helpers, only: check_moist_error, fd4_scalar, fd4_offsets
   use test_moz_potential_kernel, only: make_solute, make_solvent, make_grid, make_potential, &
      & sigma_u, eps_u, sigma_v, eps_v, q_v, weights_3d, check_close
   implicit none(type, external)
   private
   public :: collect_moz_potential_ec

contains

   !> Register EC tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_potential_ec(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("requirements_and_refresh", check_requirements), &
         & new_unittest("fit_shell", check_fit_shell), &
         & new_unittest("restraints", check_restraints), &
         & new_unittest("tail_dipole", check_tail_dipole), &
         & new_unittest("fit_adjoint_fd", check_adjoint), &
         & new_unittest("active_constraint_adjoint_fd", check_active_adjoint), &
         & new_unittest("bounds_and_release_fit", check_bounds_fit), &
         & new_unittest("bound_and_kink_adjoint_fd", check_bound_adjoint), &
         & new_unittest("bad_upper_bounds", check_bad_upper_bounds, should_fail=.true.), &
         & new_unittest("bad_inverted_bounds", check_bad_inverted_bounds, should_fail=.true.), &
         & new_unittest("guards_and_single_atom", check_guards), &
         & new_unittest("stale_fit_and_phase_guards", check_stale_guards), &
         & new_unittest("ng_tables_and_adjoint", check_tables)]
   end subroutine collect_moz_potential_ec

   !> Small asymmetric fit geometry with excluded and transition-region points
   !>
   !> @param[out] mol   solute structure
   !> @param[out] grid  volume samples
   subroutine fixture(mol, grid)
      !> Solute structure
      type(structure_type), intent(out) :: mol
      !> Volume samples
      type(moist_math_grid_3d_cartesian_type), intent(out) :: grid
      integer :: i
      call make_solute(mol)
      mol%charge = 0.5_wp
      grid%ngrid = 12
      grid%xyz = reshape([-1.0_wp, 0.0_wp, 0.0_wp, -2.4_wp, 0.2_wp, 0.3_wp, &
         & 2.2_wp, 1.1_wp, 0.1_wp, 0.1_wp, -2.3_wp, 0.3_wp, &
         & -3.0_wp, 1.0_wp, 1.0_wp, 3.0_wp, -1.0_wp, -0.5_wp, &
         & 1.0_wp, 3.0_wp, 0.7_wp, -1.0_wp, -3.0_wp, -0.8_wp, &
         & 6.0_wp, 0.5_wp, 2.0_wp, -6.0_wp, -0.5_wp, -1.0_wp, &
         & 0.2_wp, 6.0_wp, -2.0_wp, -0.3_wp, -6.0_wp, 1.0_wp], [3, 12])
      grid%w = [(0.3_wp + 0.1_wp*real(i, wp), i=1, 12)]
   end subroutine fixture

   !> Configure a fit shell of uniform radii
   !>
   !> @param[in,out] ec     term
   !> @param[in]     natom  number of atoms
   !> @param[in]     start  radius where each inner switch starts, bohr
   !> @param[in]     width  switch width, bohr
   !> @param[in]     outer  optional distance beyond the radii where the shell ends, bohr; far away by default
   subroutine set_shell(ec, natom, start, width, outer)
      !> Term
      type(ec_charges_type), intent(inout) :: ec
      !> Number of atoms
      integer, intent(in) :: natom
      !> Radius where each inner switch starts
      real(wp), intent(in) :: start
      !> Switch width
      real(wp), intent(in) :: width
      !> Optional distance beyond the radii where the shell ends
      real(wp), intent(in), optional :: outer
      ec%radii = spread(start, 1, natom)
      ec%inner_scale = 1.0_wp
      ec%switch_width = width
      ec%outer_distance = 100.0_wp
      if (present(outer)) ec%outer_distance = outer
   end subroutine set_shell

   !> Check default shell weights against an explicitly weighted fit and refuse invalid shells
   !>
   !> - Two hydrogens: default 1.2-scaled Bondi radius, inner and outer switch midpoints
   !> - Mask: zero inside the inner switch, one in the shell, zero beyond the outer switch
   !>
   !> @param[out] error Test failure
   subroutine check_fit_shell(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid, weighted
      type(ec_charges_type) :: ec, control
      type(moist_error), allocatable :: err
      real(wp), allocatable :: q(:), expected(:)
      real(wp) :: phi(5), settings(4, 6), radius
      integer :: i

      call new(mol, [1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      radius = 1.32_wp*aatoau
      grid%ngrid = 5
      allocate (grid%xyz(3, 5), source=0.0_wp)
      grid%xyz(1, :) = [1.0_wp, radius + 0.5_wp, 4.5_wp, radius + 4.5_wp, 9.0_wp]
      grid%w = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
      phi = [100.0_wp, 0.3_wp, -0.1_wp, 0.2_wp, 100.0_wp]
      ec%charge_penalty = 0.0_wp
      call check(error, ec%inner_scale == 1.0_wp .and. ec%outer_distance == 5.0_wp .and. ec%switch_width == 1.0_wp &
         & .and. ec%restraint_strength == 1.0e-4_wp, more="EC shell defaults changed")
      if (allocated(error)) return
      call fit_values(ec, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      weighted = grid
      weighted%w = [0.0_wp, 0.5_wp, 1.0_wp, 0.5_wp, 0.0_wp]
      control = ec
      call set_shell(control, mol%nat, 1.0e-3_wp, 1.0e-3_wp)
      call fit_values(control, mol, weighted, phi, expected, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_close(error, q, expected, 1.0e-12_wp, "EC default shell weights")
      if (allocated(error)) return

      ! Columns: inner_scale, outer_distance, switch_width, restraint_strength
      settings = reshape([0.0_wp, 5.0_wp, 1.0_wp, 1.0e-4_wp, 1.0_wp, -1.0_wp, 1.0_wp, 1.0e-4_wp, &
         & 1.0_wp, 5.0_wp, 0.0_wp, 1.0e-4_wp, 1.0_wp, 5.0_wp, 1.0_wp, 0.0_wp, &
         & ieee_value(1.0_wp, ieee_quiet_nan), 5.0_wp, 1.0_wp, 1.0e-4_wp, &
         & 1.0_wp, 5.0_wp, 1.0_wp, 1.0e-4_wp], [4, 6])
      do i = 1, size(settings, 2)
         ec%inner_scale = settings(1, i)
         ec%outer_distance = settings(2, i)
         ec%switch_width = settings(3, i)
         ec%restraint_strength = settings(4, i)
         ! Last case: valid shell, invalid radii
         if (i == size(settings, 2)) ec%radii = [1.0_wp, -1.0_wp]
         call ec%update(mol, err)
         call check(error, allocated(err), more="EC accepted an invalid shell setting")
         if (allocated(error)) return
      end do
      ec%radii = [1.0_wp, 1.0_wp, 1.0_wp]
      call ec%update(mol, err)
      call check(error, allocated(err), more="EC accepted radii of the wrong length")
   end subroutine check_fit_shell

   !> Analytic host potential with nonzero fit residuals
   !>
   !> @param[in]  grid  volume samples
   !> @param[out] phi   potential (ngrid)
   !> @param[out] dphi  point derivative (3, ngrid)
   subroutine host_field(grid, phi, dphi)
      !> Volume samples
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
      !> Potential
      real(wp), allocatable, intent(out) :: phi(:)
      !> Point derivative
      real(wp), allocatable, intent(out) :: dphi(:, :)
      real(wp) :: d(3), r
      integer :: g
      allocate (phi(grid%ngrid), dphi(3, grid%ngrid))
      do g = 1, grid%ngrid
         d = grid%xyz(:, g)
         r = sqrt(4.0_wp + sum(d**2))
         phi(g) = 0.1_wp*sin(0.4_wp*d(1)) + 0.01_wp*d(2)**2 - 0.01_wp*d(3) + 0.4_wp/r
         dphi(:, g) = -0.4_wp*d/r**3
         dphi(1, g) = dphi(1, g) + 0.04_wp*cos(0.4_wp*d(1))
         dphi(2, g) = dphi(2, g) + 0.02_wp*d(2)
         dphi(3, g) = dphi(3, g) - 0.01_wp
      end do
   end subroutine host_field

   !> Update with the grid, declare and answer an EC term without host charges
   !>
   !> @param[in,out] ec        term
   !> @param[in]     mol       solute structure
   !> @param[in]     grid      volume samples
   !> @param[in]     phi       host potential (ngrid)
   !> @param[in]     dphi      host point derivative (3, ngrid)
   !> @param[in,out] coupling  coupling
   !> @param[out]    err       setup failure
   !> @param[in]     phase     optional phase; gradient by default
   subroutine prepare(ec, mol, grid, phi, dphi, coupling, err, phase)
      !> Term
      type(ec_charges_type), intent(inout) :: ec
      !> Solute structure
      type(structure_type), intent(in) :: mol
      !> Volume samples
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
      !> Host potential
      real(wp), intent(in) :: phi(:)
      !> Host point derivative
      real(wp), intent(in) :: dphi(:, :)
      !> Coupling
      type(coupling_type), intent(inout) :: coupling
      !> Setup failure
      type(moist_error), allocatable, intent(out) :: err
      !> Optional phase
      integer, intent(in), optional :: phase
      integer :: selected
      call ec%update(mol, err)
      if (.not. allocated(err)) call ec%update_grid(grid, err)
      if (allocated(err)) return
      call coupling_begin_registration(coupling)
      call coupling_set_scope(coupling, 1)
      call ec%declare_pass(grid, coupling, err)
      if (allocated(err)) return
      call coupling_snapshot(coupling, ngrid=grid%ngrid, natom=mol%nat)
      selected = moist_phase_gradient
      if (present(phase)) selected = phase
      call answer(coupling, selected, phi, dphi, err)
   end subroutine prepare

   !> Answer only potential outputs in a staged coupling
   !>
   !> @param[in,out] coupling  coupling
   !> @param[in]     phase     phase to stage
   !> @param[in]     phi       potential (ngrid)
   !> @param[in]     dphi      point derivative (3, ngrid)
   !> @param[out]    err       staging, unexpected request or answer failure
   subroutine answer(coupling, phase, phi, dphi, err)
      !> Coupling
      type(coupling_type), intent(inout) :: coupling
      !> Phase to stage
      integer, intent(in) :: phase
      !> Potential
      real(wp), intent(in) :: phi(:)
      !> Point derivative
      real(wp), intent(in) :: dphi(:, :)
      !> Error handling
      type(moist_error), allocatable, intent(out) :: err
      class(coupling_request_type), allocatable :: item
      call coupling_arm(coupling, phase, err)
      if (allocated(err)) return
      do while (coupling%next())
         item = coupling%request()
         if (item%name() /= "point_potential") then
            call fatal_error(err, "EC requested host charges or an unexpected quantity")
            return
         end if
         if (item%is_missing("phi")) call coupling%answer("phi", phi, err)
         if (allocated(err)) return
         if (item%is_missing("dphi_dr")) call coupling%answer("dphi_dr", dphi, err)
         if (allocated(err)) return
      end do
   end subroutine answer

   !> Check phase requirements and refreshed fitted charges
   subroutine check_requirements(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(ec_charges_type) :: ec
      type(coupling_type), target :: coupling
      type(coupling_view_type) :: view
      class(coupling_request_type), allocatable :: item
      real(wp), allocatable :: phi(:), dphi(:, :), q(:), q_old(:), got(:)
      integer :: phase
      call fixture(mol, grid)
      call set_shell(ec, mol%nat, 1.0_wp, 1.0_wp)
      call host_field(grid, phi, dphi)
      call check(error,.not. allocated(ec%reference_charges), more="default EC requires no supplied charges")
      if (allocated(error)) return
      call prepare(ec, mol, grid, phi, dphi, coupling, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call coupling_begin_registration(coupling)
      call coupling_set_scope(coupling, 1)
      call ec%declare_pass(grid, coupling, err)
      call coupling_snapshot(coupling, ngrid=grid%ngrid, natom=mol%nat)
      call coupling_invalidate(coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do phase = moist_phase_energy, moist_phase_gradient
         call coupling_arm(coupling, phase, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, coupling%next(), more="phase request was absent")
         if (allocated(error)) return
         item = coupling%request()
         call check(error, item%name() == "point_potential" .and. item%is_missing("phi"), &
                    more="wrong potential request")
         if (allocated(error)) return
         call check(error, item%is_missing("dphi_dr") .eqv. (phase == moist_phase_gradient), &
                    more="wrong derivative phase")
         if (allocated(error)) return
         call check(error,.not. coupling%next(), more="EC registered an extra request")
         if (allocated(error)) return
      end do
      call answer(coupling, moist_phase_response, phi, dphi, err)
      call coupling_make_view(coupling, 1, view)
      if (.not. allocated(err)) call ec%tail_charges(view, q_old, err)
      if (.not. allocated(err)) call ec%field(grid, view, got, err)
      call coupling_close_view(coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(got - phi)) < 1.0e-14_wp)
      if (allocated(error)) return
      phi = phi + 0.02_wp*grid%xyz(1, :)
      call coupling_invalidate(coupling)
      call answer(coupling, moist_phase_response, phi, dphi, err)
      call coupling_make_view(coupling, 1, view)
      if (.not. allocated(err)) call ec%tail_charges(view, q, err)
      call coupling_close_view(coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, abs(sum(q) - mol%charge) < 1.0e-12_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(q - q_old)) > 1.0e-3_wp, more="EC cached host-dependent charges")
   end subroutine check_requirements

   !> Fit charges through the term interface
   !>
   !> @param[in]  config  configured term
   !> @param[in]  mol     structure
   !> @param[in]  grid    volume samples
   !> @param[in]  phi     host potential (ngrid)
   !> @param[out] q       fitted charges (natom)
   !> @param[out] err     fit failure
   subroutine fit_values(config, mol, grid, phi, q, err)
      !> Configured term
      type(ec_charges_type), intent(in) :: config
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Volume samples
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
      !> Host potential
      real(wp), intent(in) :: phi(:)
      !> Fitted charges
      real(wp), allocatable, intent(out) :: q(:)
      !> Fit failure
      type(moist_error), allocatable, intent(out) :: err
      type(ec_charges_type) :: ec
      type(coupling_type), target :: coupling
      type(coupling_view_type) :: view
      real(wp) :: dphi(3, grid%ngrid)
      ec = config
      dphi = 0.0_wp
      call prepare(ec, mol, grid, phi, dphi, coupling, err, moist_phase_response)
      if (allocated(err)) return
      call coupling_make_view(coupling, 1, view)
      call ec%tail_charges(view, q, err)
      call coupling_close_view(coupling)
   end subroutine fit_values

   !> Point-charge potential samples for controlled fit cases
   !>
   !> @param[in]  mol     charge centers
   !> @param[in]  grid    sample points
   !> @param[in]  q       generating charges (natom)
   !> @param[out] phi     sampled potential (ngrid)
   subroutine charge_field(mol, grid, q, phi)
      !> Charge centers
      type(structure_type), intent(in) :: mol
      !> Sample points
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
      !> Generating charges
      real(wp), intent(in) :: q(:)
      !> Sampled potential
      real(wp), allocatable, intent(out) :: phi(:)
      integer :: a, g
      allocate (phi(grid%ngrid), source=0.0_wp)
      do g = 1, grid%ngrid
         do a = 1, mol%nat
            phi(g) = phi(g) + q(a)/max(0.1_wp, norm2(grid%xyz(:, g) - mol%xyz(:, a)))
         end do
      end do
   end subroutine charge_field

   !> Check the reference restraint, its weight-scale invariance and element penalty thresholds
   subroutine check_restraints(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid, scaled
      type(ec_charges_type) :: ec
      real(wp), allocatable :: phi(:), q(:), baseline(:)
      integer :: species, z
      call new(mol, [6, 6], reshape([-1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      grid%ngrid = 4
      grid%xyz = reshape([-2.2_wp, 0.0_wp, 0.0_wp, 2.2_wp, 0.3_wp, 0.0_wp, &
         & -8.0_wp, 0.0_wp, 0.5_wp, 8.0_wp, 0.0_wp, 0.0_wp], [3, 4])
      grid%w = [1.0_wp, 1.0_wp, 512.0_wp, 512.0_wp]
      call set_shell(ec, mol%nat, 0.1_wp, 0.1_wp)
      ec%restraint_strength = 0.1_wp
      ec%charge_penalty = 0.0_wp
      ec%reference_charges = [0.1_wp, -0.1_wp]
      call charge_field(mol, grid, [0.3_wp, -0.3_wp], phi)
      phi(3:4) = phi(3:4)*3.0_wp
      call fit_values(ec, mol, grid, phi, baseline, err)
      ! Uniformly scaled quadrature weights leave the relative restraint and the charges unchanged
      scaled = grid
      scaled%w = 1.0e3_wp*grid%w
      if (.not. allocated(err)) call fit_values(ec, mol, scaled, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_close(error, q, baseline, 1.0e-12_wp, "EC fit under scaled quadrature weights")
      if (allocated(error)) return
      call check(error, maxval(abs(baseline - ec%reference_charges)) > 1.0e-2_wp, &
         & more="relative restraint pinned the charges to the reference")
      if (allocated(error)) return
      ec%restraint_strength = 1.0e8_wp
      call fit_values(ec, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(q - ec%reference_charges)) < 1.0e-7_wp)
      if (allocated(error)) return
      deallocate (ec%reference_charges)
      ec%restraint_strength = 1.0e-8_wp
      ec%charge_penalty = 100.0_wp
      do species = 1, 3
         select case (species)
         case (1)
            z = 1
         case (2)
            z = 3
         case default
            z = 8
         end select
         call new(mol, [z, z], reshape([-1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
         call charge_field(mol, grid, [8.0_wp, -8.0_wp], phi)
         call fit_values(ec, mol, grid, phi, q, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         if (z == 8) z = 4
         call check(error, abs(q(1) - real(z, wp)) < 1.0e-8_wp, more="incorrect element penalty threshold")
         if (allocated(error)) return
      end do
   end subroutine check_restraints

   !> Check that the default fit keeps the solute dipole in the tail charges
   !>
   !> - Host potential of point charges and axial point dipoles, so the charges cannot fit it exactly
   !> - Default shell and restraint on a Cartesian grid covering the shell
   !> - Tail dipole within 1 % of the host dipole and stable over a decade of restraint
   subroutine check_tail_dipole(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(ec_charges_type) :: ec
      real(wp), allocatable :: phi(:), q(:), q_weak(:)
      real(wp) :: source(3), dipole(3, 3), host(3), tail(3), weak(3), d(3)
      integer :: g, a
      ! HCN, bohr
      call new(mol, [1, 6, 7], reshape([0.0_wp, 0.0_wp, -2.015_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 2.185_wp], [3, 3]))
      call new_cartesian_point_grid(grid, err, nx=40, ny=40, nz=40, dr=0.6_wp, margin=0.0_wp)
      if (.not. allocated(err)) call grid%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      source = [0.25_wp, 0.05_wp, -0.30_wp]
      dipole = reshape([0.0_wp, 0.0_wp, 0.05_wp, 0.0_wp, 0.0_wp, -0.1_wp, 0.0_wp, 0.0_wp, -0.3_wp], [3, 3])
      allocate (phi(grid%ngrid), source=0.0_wp)
      do g = 1, grid%ngrid
         do a = 1, mol%nat
            d = grid%xyz(:, g) - mol%xyz(:, a)
            phi(g) = phi(g) + source(a)/max(0.1_wp, norm2(d)) + dot_product(dipole(:, a), d)/max(0.1_wp, norm2(d))**3
         end do
      end do
      host = matmul(mol%xyz, source) + sum(dipole, dim=2)
      ec%charge_penalty = 0.0_wp
      call fit_values(ec, mol, grid, phi, q, err)
      if (.not. allocated(err)) then
         ec%restraint_strength = 1.0e-5_wp
         call fit_values(ec, mol, grid, phi, q_weak, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) return
      tail = matmul(mol%xyz, q)
      weak = matmul(mol%xyz, q_weak)
      call check(error, norm2(tail - host) < 1.0e-2_wp*norm2(host), more="EC tail charges lost the host dipole")
      if (allocated(error)) return
      call check(error, norm2(tail - weak) < 1.0e-2_wp*norm2(host), &
         & more="EC tail dipole depends on the restraint strength")
   end subroutine check_tail_dipole

   !> Scalar functional of the full potential and fitted charges
   !>
   !> @param[in]  ec     configured term
   !> @param[in]  mol    structure
   !> @param[in]  grid   volume samples
   !> @param[in]  phi    host potential (ngrid)
   !> @param[in]  w_phi  field weights (ngrid)
   !> @param[in]  w_q    charge weights (natom)
   !> @param[out] value  scalar functional
   !> @param[out] err    fit failure
   subroutine functional(ec, mol, grid, phi, w_phi, w_q, value, err)
      !> Configured term
      type(ec_charges_type), intent(in) :: ec
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Volume samples
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
      !> Host potential
      real(wp), intent(in) :: phi(:)
      !> Field weights
      real(wp), intent(in) :: w_phi(:)
      !> Charge weights
      real(wp), intent(in) :: w_q(:)
      !> Scalar functional
      real(wp), intent(out) :: value
      !> Fit failure
      type(moist_error), allocatable, intent(out) :: err
      real(wp), allocatable :: q(:)
      call fit_values(ec, mol, grid, phi, q, err)
      if (.not. allocated(err)) value = dot_product(w_phi, phi) + dot_product(w_q, q)
   end subroutine functional

   !> Check the smooth-fit adjoints
   subroutine check_adjoint(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(ec_charges_type) :: ec
      call set_shell(ec, 3, 1.0_wp, 1.7_wp, 4.0_wp)
      call adjoint_case(error, ec, 0)
      if (allocated(error)) return
      ! Strong restraint off the fitted charges: the derivative of lambda through trace(H_data)
      ec%restraint_strength = 0.3_wp
      ec%reference_charges = [0.4_wp, -0.2_wp, 0.3_wp]
      call adjoint_case(error, ec, 0)
   end subroutine check_adjoint

   !> Check penalty-kink and hard-bound adjoints
   subroutine check_active_adjoint(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(ec_charges_type) :: ec
      ! Lower bounds, negative kinks, releases and the all-fixed branch: check_bound_adjoint, check_bounds_fit
      call set_shell(ec, 3, 1.0_wp, 1.0_wp)
      ec%charge_penalty = 0.0_wp
      ec%upper_bounds = [0.03_wp, 10.0_wp, 10.0_wp]
      call adjoint_case(error, ec, 1)
      if (allocated(error)) return
      deallocate (ec%upper_bounds)
      ec%charge_thresholds = [0.03_wp, 10.0_wp, 10.0_wp]
      ec%charge_penalty = 10.0_wp
      call adjoint_case(error, ec, 2)
   end subroutine check_active_adjoint

   !> Fit with active hard bounds, fixed charges and released starting constraints
   !>
   !> - Charges respect the bounds and sum exactly to the total charge
   !> - The generating charges violate the bound, so an unbounded fit leaves it
   subroutine check_bounds_fit(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(ec_charges_type) :: ec, bounded
      real(wp), allocatable :: phi(:), free(:), q(:)
      call set_shell(ec, 3, 1.0_wp, 1.0_wp)
      call fixture(mol, grid)
      call charge_field(mol, grid, [2.0_wp, -1.0_wp, -0.5_wp], phi)
      ec%charge_penalty = 0.0_wp
      call fit_values(ec, mol, grid, phi, free, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, free(2) < -0.25_wp .and. free(1) > 0.5_wp, more="control fit does not violate the bounds")
      if (allocated(error)) return
      ! Lower bound of the second charge alone
      bounded = ec
      bounded%lower_bounds = [-10.0_wp, -0.2_wp, -10.0_wp]
      call fit_values(bounded, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, q(2), -0.2_wp, more="lower bound not active", thr=1.0e-9_wp)
      if (allocated(error)) return
      call check(error, sum(q), mol%charge, more="total charge broken by a lower bound", thr=1.0e-12_wp)
      if (allocated(error)) return
      ! Lower and upper bounds pin two charges; the third follows from the total charge
      bounded%lower_bounds = [-10.0_wp, -0.05_wp, -10.0_wp]
      bounded%upper_bounds = [0.4_wp, 10.0_wp, 10.0_wp]
      call fit_values(bounded, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_close(error, q, [0.4_wp, -0.05_wp, mol%charge - 0.35_wp], 1.0e-9_wp, "two pinned charges")
      if (allocated(error)) return
      ! Bounds summing to the total charge leave a single feasible point
      bounded = ec
      bounded%lower_bounds = [0.375_wp, 0.0625_wp, 0.0625_wp]
      call fit_values(bounded, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_close(error, q, bounded%lower_bounds, 1.0e-12_wp, "single feasible point")
      if (allocated(error)) return
      ! Equal bounds fix one charge
      bounded = ec
      bounded%lower_bounds = [0.25_wp, -10.0_wp, -10.0_wp]
      bounded%upper_bounds = [0.25_wp, 10.0_wp, 10.0_wp]
      call fit_values(bounded, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, q(1), 0.25_wp, more="equal bounds did not fix the charge", thr=1.0e-12_wp)
      if (allocated(error)) return
      call check(error, sum(q), mol%charge, more="total charge broken by a fixed charge", thr=1.0e-12_wp)
      if (allocated(error)) return
      ! A reference below the bound starts the charge on it; the data release it again
      bounded = ec
      bounded%restraint_strength = 1.0e-4_wp
      bounded%reference_charges = [1.2_wp, -1.5_wp, 0.8_wp]
      bounded%lower_bounds = [-10.0_wp, -1.2_wp, -10.0_wp]
      call fit_values(bounded, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, q(2) > -1.1_wp, more="starting constraint was not released")
      if (allocated(error)) return
      call check(error, sum(q), mol%charge, more="total charge broken after a release", thr=1.0e-12_wp)
      if (allocated(error)) return
      ! References above the upper bounds start two charges on them; the data release them again
      bounded = ec
      bounded%restraint_strength = 1.0e-4_wp
      bounded%reference_charges = [3.5_wp, 0.5_wp, -3.5_wp]
      bounded%charge_thresholds = [1.0_wp, 0.3_wp, 10.0_wp]
      bounded%upper_bounds = [2.5_wp, -0.5_wp, 10.0_wp]
      call fit_values(bounded, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, q(1) < 2.4_wp .and. q(2) < -0.6_wp, more="starting upper constraints were not released")
      if (allocated(error)) return
      call check(error, sum(q), mol%charge, more="total charge broken after an upper release", thr=1.0e-12_wp)
      if (allocated(error)) return
      ! Every charge starts on a bound and the point is not optimal: all-fixed multiplier check releases two
      bounded = ec
      bounded%restraint_strength = 1.0e-4_wp
      bounded%reference_charges = [-3.0_wp, 3.0_wp, 3.0_wp]
      bounded%upper_bounds = [10.0_wp, 0.125_wp, 0.3125_wp]
      bounded%lower_bounds = [0.0625_wp, -10.0_wp, -10.0_wp]
      call fit_values(bounded, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(q <= bounded%upper_bounds + 1.0e-12_wp) .and. all(q >= bounded%lower_bounds - 1.0e-12_wp), &
         & more="bounds violated after releasing the fixed face")
      if (allocated(error)) return
      call check(error, q(1) > 0.1_wp, more="all-fixed start was not released")
      if (allocated(error)) return
      call check(error, sum(q), mol%charge, more="total charge broken after the all-fixed release", thr=1.0e-12_wp)
      if (allocated(error)) return
      ! Two excess-charge thresholds: one pinned at its kink, the other region released
      bounded = ec
      bounded%charge_penalty = 10.0_wp
      bounded%charge_thresholds = [0.03_wp, 0.1_wp, 10.0_wp]
      call fit_values(bounded, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, q(1), 0.03_wp, more="positive penalty kink not active", thr=1.0e-9_wp)
      if (allocated(error)) return
      call check(error, sum(q), mol%charge, more="total charge broken at a penalty kink", thr=1.0e-12_wp)
      if (allocated(error)) return
   end subroutine check_bounds_fit

   !> Fail an expected-failure test only on the targeted library error
   !>
   !> - Used by `should_fail=.true.` tests, where test-drive inverts the
   !>   verdict: a raised test failure passes, a clean return fails
   !> - Fails the test only when `err` names `expected`; no error or a
   !>   different one returns cleanly, so the test is reported as failed
   !>
   !> @param[out] error     test failure, set only on the expected error
   !> @param[in]  err       library error, possibly unallocated
   !> @param[in]  expected  substring the expected error message must contain
   subroutine expect_error(error, err, expected)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error, possibly unallocated
      type(moist_error), allocatable, intent(in) :: err
      !> Substring the expected error message must contain
      character(len=*), intent(in) :: expected

      if (.not. allocated(err)) return
      if (index(err%message, expected) > 0) call test_failed(error, err%message)
   end subroutine expect_error

   !> Upper bounds that sum below the total charge are refused when the term is updated
   subroutine check_bad_upper_bounds(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(ec_charges_type) :: ec
      call fixture(mol, grid)
      ec%charge_penalty = 0.0_wp
      ec%upper_bounds = [0.1_wp, 0.1_wp, 0.1_wp]
      call ec%update(mol, err)
      call expect_error(error, err, "incompatible with the total charge")
   end subroutine check_bad_upper_bounds

   !> A lower bound above its upper bound is refused when the term is updated
   subroutine check_bad_inverted_bounds(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(ec_charges_type) :: ec
      call fixture(mol, grid)
      ec%charge_penalty = 0.0_wp
      ec%lower_bounds = [0.2_wp, -10.0_wp, -10.0_wp]
      ec%upper_bounds = [0.1_wp, 10.0_wp, 10.0_wp]
      call ec%update(mol, err)
      call expect_error(error, err, "invalid thresholds or bounds")
   end subroutine check_bad_inverted_bounds

   !> Compare derivatives of fits with lower bounds, fixed charges, negative kinks and released constraints
   subroutine check_bound_adjoint(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(ec_charges_type) :: ec
      call set_shell(ec, 3, 1.0_wp, 1.0_wp)
      ec%charge_penalty = 0.0_wp
      ec%lower_bounds = [-10.0_wp, -0.2_wp, -10.0_wp]
      call adjoint_case(error, ec, 1, 2, -0.2_wp)
      if (allocated(error)) return
      deallocate (ec%lower_bounds)
      ec%lower_bounds = [0.25_wp, -10.0_wp, -10.0_wp]
      ec%upper_bounds = [0.25_wp, 10.0_wp, 10.0_wp]
      call adjoint_case(error, ec, 1, 1, 0.25_wp)
      if (allocated(error)) return
      deallocate (ec%lower_bounds, ec%upper_bounds)
      ec%lower_bounds = [0.375_wp, 0.0625_wp, 0.0625_wp]
      call adjoint_case(error, ec, 1, 1, 0.375_wp)
      if (allocated(error)) return
      deallocate (ec%lower_bounds)
      ec%restraint_strength = 1.0e-4_wp
      ec%reference_charges = [1.2_wp, -1.5_wp, 0.8_wp]
      ec%lower_bounds = [-10.0_wp, -1.2_wp, -10.0_wp]
      call adjoint_case(error, ec, 1, 0)
      if (allocated(error)) return
      deallocate (ec%reference_charges, ec%lower_bounds)
      ec%reference_charges = [3.5_wp, 0.5_wp, -3.5_wp]
      ec%charge_thresholds = [1.0_wp, 0.3_wp, 10.0_wp]
      ec%upper_bounds = [2.5_wp, -0.5_wp, 10.0_wp]
      call adjoint_case(error, ec, 1, 0)
      if (allocated(error)) return
      deallocate (ec%reference_charges, ec%charge_thresholds, ec%upper_bounds)
      ec%reference_charges = [-3.0_wp, 3.0_wp, 3.0_wp]
      ec%upper_bounds = [10.0_wp, 0.125_wp, 0.3125_wp]
      ec%lower_bounds = [0.0625_wp, -10.0_wp, -10.0_wp]
      call adjoint_case(error, ec, 1, 0)
      if (allocated(error)) return
      deallocate (ec%reference_charges, ec%upper_bounds, ec%lower_bounds)
      ec%restraint_strength = 0.01_wp
      ec%charge_penalty = 10.0_wp
      ec%charge_thresholds = [10.0_wp, 0.3_wp, 10.0_wp]
      call adjoint_case(error, ec, 2, 2, -0.3_wp)
      if (allocated(error)) return
      ! Two kinks: the positive pin must survive the release of the other penalty region
      deallocate (ec%charge_thresholds)
      ec%charge_thresholds = [0.03_wp, 0.1_wp, 10.0_wp]
      call adjoint_case(error, ec, 2, 1, 0.03_wp)
   end subroutine check_bound_adjoint

   !> Compare potential, nuclear, grid-position and volume-weight derivatives with finite differences
   !> @param[in]  config configured term
   !> @param[in]  mode   0 smooth, otherwise a constrained fit with the pinned charge checked
   !> @param[in]  pin    optional pinned atom; zero skips the pinned-charge check
   !> @param[in]  pinned optional pinned charge
   subroutine adjoint_case(error, config, mode, pin, pinned)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Configured term
      type(ec_charges_type), intent(in) :: config
      !> Fit mode
      integer, intent(in) :: mode
      !> Atom whose charge the active constraint pins; first atom by default
      integer, intent(in), optional :: pin
      !> Pinned charge value; 0.03 by default
      real(wp), intent(in), optional :: pinned
      type(moist_error), allocatable :: err
      type(structure_type) :: mol, displaced
      type(moist_math_grid_3d_cartesian_type) :: grid, moved
      type(ec_charges_type) :: ec
      type(coupling_type), target :: coupling
      type(coupling_view_type) :: view
      type(response_type) :: response
      type(volume_adjoint_type) :: acc
      class(response_item_type), allocatable :: item
      real(wp), allocatable :: phi(:), dphi(:, :), w_phi(:), returned(:), sampled(:), sampled_deriv(:, :), q(:)
      real(wp) :: w_q(3), gradient(3, 3), val(4), analytic, fd
      !> Stencil step; four-point central differences
      real(wp), parameter :: h = 1.0e-4_wp
      character(len=80) :: label
      integer :: i, j, kind, s, atom
      real(wp) :: expected
      atom = 1
      expected = 0.03_wp
      if (present(pin)) atom = pin
      if (present(pinned)) expected = pinned
      call fixture(mol, grid)
      call host_field(grid, phi, dphi)
      if (mode /= 0) call charge_field(mol, grid, [2.0_wp, -1.0_wp, -0.5_wp], phi)
      if (mode /= 0) dphi = 0.0_wp
      ec = config
      call prepare(ec, mol, grid, phi, dphi, coupling, err)
      call acc%init(grid%ngrid)
      w_phi = 0.02_wp*grid%xyz(2, :)
      w_q = [0.7_wp, -0.4_wp, 0.2_wp]
      gradient = 0.0_wp
      call coupling_make_view(coupling, 1, view)
      if (.not. allocated(err)) call ec%tail_charges(view, q, err)
      if (.not. allocated(err)) call ec%accumulate_adjoint(grid, view, w_phi, w_q, response, err, acc, gradient)
      call coupling_close_view(coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      if (mode /= 0 .and. atom > 0) then
         call check(error, abs(q(atom) - expected) < 1.0e-9_wp, more="active constraint was not exercised")
         if (allocated(error)) return
      end if
      call check(error, response%next())
      if (allocated(error)) return
      item = response%item()
      select type (item)
      type is (potential_adjoint_response_type)
         returned = item%w_phi
      class default
         call check(error, .false., more="EC returned a charge adjoint instead of a potential adjoint")
         return
      end select
      call check(error,.not. response%next())
      if (allocated(error)) return
      ! Routes: 1 potential, 2 nuclei, 3 grid points, 4 volume weights
      do kind = 1, 4
         do i = 1, 3
            do j = 1, 3
               if ((kind == 1 .or. kind == 4) .and. j > 1) cycle
               do s = 1, 4
                  displaced = mol
                  moved = grid
                  sampled = phi
                  select case (kind)
                  case (1)
                     sampled(i + 1) = sampled(i + 1) + fd4_offsets(s)*h
                  case (2)
                     displaced%xyz(j, i) = displaced%xyz(j, i) + fd4_offsets(s)*h
                  case (3)
                     moved%xyz(j, i + 1) = moved%xyz(j, i + 1) + fd4_offsets(s)*h
                     if (mode == 0) call host_field(moved, sampled, sampled_deriv)
                  case default
                     moved%w(i + 1) = moved%w(i + 1) + fd4_offsets(s)*h
                  end select
                  call functional(config, displaced, moved, sampled, w_phi, w_q, val(s), err)
                  call check_moist_error(error, err)
                  if (allocated(error)) return
               end do
               call fd4_scalar(val(1), val(2), val(3), val(4), h, fd, error)
               if (allocated(error)) return
               select case (kind)
               case (1)
                  analytic = returned(i + 1)
               case (2)
                  analytic = gradient(j, i)
               case (3)
                  analytic = acc%w_xyz(j, i + 1)
               case default
                  analytic = acc%w_w(i + 1)
               end select
               write (label, "(a, i0, a, i0, a, i0, a, i0)") "EC adjoint mode ", mode, " route ", kind, &
                  & " index ", i, " component ", j
               call check_close(error, [analytic], [fd], 1.0e-9_wp, trim(label))
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine adjoint_case

   !> Check invalid fits, radial refusal and the one-atom charge constraint
   subroutine check_guards(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      ! TODO: cover the iteration-limit refusal of solve_fit
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(moist_math_grid_radial_type) :: radial
      type(ec_charges_type) :: ec
      type(coupling_type) :: coupling
      real(wp), allocatable :: phi(:), q(:)
      call fixture(mol, grid)
      call set_shell(ec, mol%nat, 1.0_wp, 1.0_wp)
      call ec%update_grid(grid, err)
      call check(error, allocated(err), more="unupdated EC term built a fit")
      if (allocated(error)) return
      deallocate (err)
      call ec%declare_pass(grid, coupling, err)
      call check(error, allocated(err), more="EC term without a fit declared its request")
      if (allocated(error)) return
      ec%restraint_strength = 0.0_wp
      call ec%update(mol, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      ec%restraint_strength = 0.01_wp
      ec%lower_bounds = [1.0_wp, 1.0_wp, 1.0_wp]
      call ec%update(mol, err)
      call check(error, allocated(err), more="infeasible total charge accepted")
      if (allocated(error)) return
      deallocate (ec%lower_bounds)
      call ec%update(mol, err)
      if (.not. allocated(err)) call ec%update_grid(radial, err)
      call check(error, allocated(err), more="EC built a fit on a radial grid")
      if (allocated(error)) return
      call check(error, index(err%message, "volume grids only") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)
      call ec%declare_pass(radial, coupling, err)
      call check(error, allocated(err), more="EC declared on a radial grid")
      if (allocated(error)) return
      deallocate (err)
      grid%w = -1.0_wp
      call ec%update_grid(grid, err)
      call check(error, allocated(err), more="negative quadrature weights accepted")
      if (allocated(error)) return
      grid%w = 1.0_wp
      ec%inner_scale = 100.0_wp
      call ec%update(mol, err)
      if (.not. allocated(err)) call ec%update_grid(grid, err)
      call check(error, allocated(err), more="fully excluded fit accepted")
      if (allocated(error)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call set_shell(ec, mol%nat, 1.0_wp, 1.0_wp)
      mol%charge = 0.7_wp
      phi = grid%xyz(1, :)*0.1_wp
      call fit_values(ec, mol, grid, phi, q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, abs(q(1) - mol%charge) < 1.0e-12_wp)
   end subroutine check_guards

   !> Evaluate the full Ng-table functional through potential assembly
   !>
   !> @param[in]  mol    solute structure
   !> @param[in]  grid   volume grid
   !> @param[in]  phi    host potential (ngrid)
   !> @param[in]  solvent solvent potential
   !> @param[in]  w_sr   short-range table weights
   !> @param[in]  w_lr   real-space tail weights
   !> @param[in]  w_k    reciprocal-space tail weights
   !> @param[out] value  scalar functional
   !> @param[out] err    evaluation failure
   subroutine table_functional(mol, grid, phi, solvent, w_sr, w_lr, w_k, value, err)
      !> Solute structure
      type(structure_type), intent(in) :: mol
      !> Volume grid
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
      !> Host potential
      real(wp), intent(in) :: phi(:)
      !> Solvent potential
      type(moz_potential_type), intent(in) :: solvent
      !> Short-range weights
      real(wp), intent(in) :: w_sr(:, :)
      !> Real-space tail weights
      real(wp), intent(in) :: w_lr(:, :)
      !> Reciprocal-space tail weights
      complex(wp), intent(in) :: w_k(:, :)
      !> Scalar functional
      real(wp), intent(out) :: value
      !> Evaluation failure
      type(moist_error), allocatable, intent(out) :: err
      type(moz_potential_type) :: pot
      type(ec_charges_type) :: ec
      type(custom_lj_type) :: lj
      type(coupling_type), target :: coupling
      real(wp), allocatable :: sr(:, :), lr(:, :)
      complex(wp), allocatable :: uk(:, :)
      call set_shell(ec, mol%nat, 1.0_wp, 1.0_wp, 2.0_wp)
      call new_custom_lj(lj, sigma_u, eps_u, err)
      if (.not. allocated(err)) call pot%add(lj, err)
      if (.not. allocated(err)) call pot%add(ec, err)
      if (.not. allocated(err)) call pot%update(mol, grid, err)
      call coupling_begin_registration(coupling)
      if (.not. allocated(err)) call pot%declare(grid, coupling, err)
      call coupling_snapshot(coupling, ngrid=grid%ngrid, natom=mol%nat)
      if (.not. allocated(err)) call answer(coupling, moist_phase_response, phi, &
         & spread(phi*0.0_wp, dim=1, ncopies=3), err)
      if (.not. allocated(err)) call pot%compute(grid, sr, lr, uk, err, solvent, coupling)
      if (.not. allocated(err)) value = sum(w_sr*sr) + sum(w_lr*lr) + real(sum(conjg(w_k)*uk), wp)
   end subroutine table_functional

   !> Check Ng recombination and adjoint dispatch through the potential
   !>
   !> - Ng split on and off recombine to the same full table
   !> - Four-point central differences of the full table functional for every nuclear component,
   !>   grid-point positions and volume weights at probe points, and the host potential
   !> - Potential adjoint: one `potential_adjoint` item, kernel field adjoint plus fit part
   !> - Probe derivatives well above the tolerance, so the comparison discriminates
   subroutine check_tables(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      ! TODO: thread-count invariance of the fit on a grid spanning several blocks (1 vs several threads)
      type(moist_error), allocatable :: err
      type(structure_type) :: mol, solvent_mol, displaced
      type(moist_math_grid_3d_cartesian_type) :: grid, moved
      type(moz_potential_type) :: pot, solvent
      type(ec_charges_type) :: ec
      type(custom_lj_type) :: lj
      type(coupling_type), target :: coupling
      type(response_type) :: response
      type(volume_adjoint_type) :: acc
      class(response_item_type), allocatable :: item
      real(wp), allocatable :: phi(:), dphi(:, :), sampled(:), sampled_deriv(:, :), returned(:)
      real(wp), allocatable :: sr(:, :), lr(:, :), combined(:, :), w_sr(:, :), w_lr(:, :)
      complex(wp), allocatable :: uk(:, :), w_k(:, :)
      real(wp) :: gradient(3, 3), val(4), fd, analytic, signal
      !> Stencil step; four-point central differences
      real(wp), parameter :: h = 1.0e-4_wp
      character(len=80) :: label
      integer :: probe(3), a, c, s, ip, g, route
      call make_solute(mol)
      call set_shell(ec, mol%nat, 1.0_wp, 1.0_wp, 2.0_wp)
      mol%charge = 0.5_wp
      call make_solvent(solvent_mol)
      call make_grid(grid, mol, err)
      if (.not. allocated(err)) call make_potential(solvent, solvent_mol, err, sigma_v, eps_v, q_v)
      call host_field(grid, phi, dphi)
      if (.not. allocated(err)) call new_custom_lj(lj, sigma_u, eps_u, err)
      if (.not. allocated(err)) call pot%add(lj, err)
      if (.not. allocated(err)) call pot%add(ec, err)
      if (.not. allocated(err)) call pot%update(mol, grid, err)
      call coupling_begin_registration(coupling)
      if (.not. allocated(err)) call pot%declare(grid, coupling, err)
      call coupling_snapshot(coupling, ngrid=grid%ngrid, natom=mol%nat)
      if (.not. allocated(err)) call answer(coupling, moist_phase_gradient, phi, dphi, err)
      if (.not. allocated(err)) call pot%compute(grid, sr, lr, uk, err, solvent, coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      combined = sr + lr
      call pot%set_ng_split(.false., 1.0_wp, err)
      if (.not. allocated(err)) call pot%compute(grid, sr, lr, uk, err, solvent, coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(combined - sr)) < 1.0e-12_wp)
      if (allocated(error)) return
      call pot%set_ng_split(.true., 1.0_wp, err)
      call weights_3d(grid%ngrid, grid%npts_k, solvent_mol%nat, w_sr, w_lr, w_k)
      gradient = 0.0_wp
      call acc%init(grid%ngrid)
      if (.not. allocated(err)) call pot%compute_adjoint(grid, solvent, coupling, w_sr, w_lr, w_k, response, err, &
         & gradient=gradient, grid_adjoint=acc)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, response%next(), more="EC pushed no potential adjoint through the potential")
      if (allocated(error)) return
      item = response%item()
      select type (item)
      type is (potential_adjoint_response_type)
         returned = item%w_phi
      class default
         call check(error, .false., more="EC pushed an item other than the potential adjoint")
         return
      end select
      call check(error,.not. response%next(), more="the potential pushed more than the EC potential adjoint")
      if (allocated(error)) return

      ! Nuclear gradient, every component
      signal = 0.0_wp
      do a = 1, mol%nat
         do c = 1, 3
            do s = 1, 4
               displaced = mol
               displaced%xyz(c, a) = displaced%xyz(c, a) + fd4_offsets(s)*h
               call table_functional(displaced, grid, phi, solvent, w_sr, w_lr, w_k, val(s), err)
               call check_moist_error(error, err)
               if (allocated(error)) return
            end do
            call fd4_scalar(val(1), val(2), val(3), val(4), h, fd, error)
            if (allocated(error)) return
            write (label, "(a, i0, a, i0)") "EC tables nuclear gradient atom ", a, " component ", c
            call check_close(error, [gradient(c, a)], [fd], 1.0e-9_wp, trim(label))
            if (allocated(error)) return
            signal = max(signal, abs(fd))
         end do
      end do
      call check(error, signal > 1.0e-4_wp, more="EC nuclear gradient probe is below the discriminating scale")
      if (allocated(error)) return

      ! Host potential, grid-point positions (host field moved along) and volume weights at probe points
      probe = [grid%ngrid/2 + 4, grid%ngrid/3, grid%ngrid]
      do route = 1, 3
         signal = 0.0_wp
         do ip = 1, size(probe)
            g = probe(ip)
            do c = 1, 3
               if (route /= 2 .and. c > 1) cycle
               do s = 1, 4
                  moved = grid
                  sampled = phi
                  select case (route)
                  case (1)
                     sampled(g) = sampled(g) + fd4_offsets(s)*h
                  case (2)
                     moved%xyz(c, g) = moved%xyz(c, g) + fd4_offsets(s)*h
                     call host_field(moved, sampled, sampled_deriv)
                  case default
                     moved%w(g) = moved%w(g) + fd4_offsets(s)*h
                  end select
                  call table_functional(mol, moved, sampled, solvent, w_sr, w_lr, w_k, val(s), err)
                  call check_moist_error(error, err)
                  if (allocated(error)) return
               end do
               call fd4_scalar(val(1), val(2), val(3), val(4), h, fd, error)
               if (allocated(error)) return
               select case (route)
               case (1)
                  analytic = returned(g)
               case (2)
                  analytic = acc%w_xyz(c, g)
               case default
                  analytic = acc%w_w(g)
               end select
               write (label, "(a, i0, a, i0, a, i0)") "EC tables route ", route, " point ", g, " component ", c
               call check_close(error, [analytic], [fd], 1.0e-9_wp, trim(label))
               if (allocated(error)) return
               signal = max(signal, abs(fd))
            end do
         end do
         write (label, "(a, i0, a)") "EC tables route ", route, " probe is below the discriminating scale"
         call check(error, signal > 1.0e-4_wp, more=trim(label))
         if (allocated(error)) return
      end do
   end subroutine check_tables

   !> Check refusals of stale fits, other grids, missing phase outputs and nonfinite settings
   !>
   !> - Charges after a new `update` without `update_grid`: refusal
   !> - Declaration, field and adjoint on a grid other than the grid of the update: refusal, no rebuild
   !> - Grid adjoint in the response phase: `dphi_dr` missing, refusal instead of zeros
   !> - Nonfinite fit setting: refusal at `update`
   subroutine check_stale_guards(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid, moved
      type(ec_charges_type) :: ec
      !> Gradient-phase coupling; answers outlive phases, so the response-phase case gets its own
      type(coupling_type), target :: coupling, response_only
      type(coupling_view_type) :: view
      type(response_type) :: response
      type(volume_adjoint_type) :: acc
      real(wp), allocatable :: phi(:), dphi(:, :), q(:), w_phi(:)
      real(wp) :: w_q(3)
      call fixture(mol, grid)
      call set_shell(ec, mol%nat, 1.0_wp, 1.0_wp)
      call host_field(grid, phi, dphi)
      w_phi = 0.02_wp*grid%xyz(2, :)
      w_q = [0.7_wp, -0.4_wp, 0.2_wp]

      call prepare(ec, mol, grid, phi, dphi, coupling, err)
      if (.not. allocated(err)) call ec%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call coupling_make_view(coupling, 1, view)
      call ec%tail_charges(view, q, err)
      call coupling_close_view(coupling)
      call check(error, allocated(err), more="EC fitted charges after an update without update_grid")
      if (allocated(error)) return
      call check(error, index(err%message, "fit is not built") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      call prepare(ec, mol, grid, phi, dphi, coupling, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      moved = grid
      moved%xyz(1, 5) = moved%xyz(1, 5) + 0.1_wp
      call coupling_begin_registration(coupling)
      call coupling_set_scope(coupling, 1)
      call ec%declare_pass(moved, coupling, err)
      call check(error, allocated(err), more="EC declared on a grid other than the grid of its update")
      if (allocated(error)) return
      call check(error, err%message, "EC grid differs from the grid of its update")
      if (allocated(error)) return
      deallocate (err)
      call prepare(ec, mol, grid, phi, dphi, coupling, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call coupling_make_view(coupling, 1, view)
      call ec%field(moved, view, w_phi, err)
      call coupling_close_view(coupling)
      call check(error, allocated(err), more="EC field read on a grid other than the grid of its update")
      if (allocated(error)) return
      call check(error, err%message, "EC grid differs from the grid of its update")
      if (allocated(error)) return
      deallocate (err)
      w_phi = 0.02_wp*grid%xyz(2, :)
      call acc%init(grid%ngrid)
      call coupling_make_view(coupling, 1, view)
      call ec%accumulate_adjoint(moved, view, w_phi, w_q, response, err, acc)
      call coupling_close_view(coupling)
      call check(error, allocated(err), more="EC adjoint accepted a grid other than the grid of its update")
      if (allocated(error)) return
      call check(error, err%message, "EC grid differs from the grid of its update")
      if (allocated(error)) return
      deallocate (err)

      call prepare(ec, mol, grid, phi, dphi, response_only, err, moist_phase_response)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call coupling_make_view(response_only, 1, view)
      call ec%accumulate_adjoint(grid, view, w_phi, w_q, response, err, acc)
      call coupling_close_view(response_only)
      call check(error, allocated(err), more="EC grid adjoint ran without the gradient-phase dphi_dr")
      if (allocated(error)) return
      call check(error, index(err%message, "dphi_dr") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      ec%outer_distance = ieee_value(1.0_wp, ieee_quiet_nan)
      call ec%update(mol, err)
      call check(error, allocated(err), more="EC accepted a nonfinite shell distance")
      if (allocated(error)) return
      call check(error, err%message, "EC fit settings must be finite")
   end subroutine check_stale_guards

end module test_moz_potential_ec
