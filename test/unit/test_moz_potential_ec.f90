!> Tests of constrained EC fitting, coupling requirements and fit adjoints
module test_moz_potential_ec
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type, fatal_error
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_channels_coupling, only: coupling_type, coupling_view_type, coupling_request_type, &
      & coupling_begin_registration, coupling_set_scope, coupling_snapshot, coupling_arm, coupling_invalidate, &
      & coupling_make_view, coupling_close_view, moist_phase_energy, moist_phase_response, moist_phase_gradient
   use moist_channels_response, only: response_type, response_item_type, potential_adjoint_response_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_model_moz_potential, only: moz_potential_type, ec_charges_type, custom_lj_type, new_custom_lj
   use test_helpers, only: check_moist_error, fd4_scalar, fd4_offsets
   use test_moz_potential_kernel, only: make_solute, make_solvent, make_grid, make_potential, &
      & sigma_u, eps_u, sigma_v, eps_v, q_v, weights_3d, check_close
   implicit none(type, external)
   private
   public :: collect_moz_potential_ec

   ! TODO: this should be covered better by tests.

contains

   !> Register EC tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_potential_ec(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("requirements_and_refresh", check_requirements), &
         & new_unittest("long_range_and_restraints", check_restraints), &
         & new_unittest("fit_adjoint_fd", check_adjoint), &
         & new_unittest("active_constraint_adjoint_fd", check_active_adjoint), &
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
   !>
   !> @param[out] error  test error
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

   !> Check distant-sample priority, reference restraint and element penalty thresholds
   !>
   !> @param[out] error  test error
   subroutine check_restraints(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(ec_charges_type) :: ec
      real(wp), allocatable :: phi(:), q(:), baseline(:)
      integer :: species, z
      call new(mol, [6, 6], reshape([-1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      grid%ngrid = 4
      grid%xyz = reshape([-2.2_wp, 0.0_wp, 0.0_wp, 2.2_wp, 0.0_wp, 0.0_wp, &
         & -8.0_wp, 0.0_wp, 0.0_wp, 8.0_wp, 0.0_wp, 0.0_wp], [3, 4])
      grid%w = [1.0_wp, 1.0_wp, 512.0_wp, 512.0_wp]
      ec%exclusion_radius = 0.1_wp
      ec%restraint_strength = 1.0e-8_wp
      ec%charge_penalty = 0.0_wp
      ec%long_range_bias = 0.0_wp
      call charge_field(mol, grid, [0.3_wp, -0.3_wp], phi)
      phi(3:4) = phi(3:4)*3.0_wp
      call fit_values(ec, mol, grid, phi, baseline, err)
      if (.not. allocated(err)) then
         ec%long_range_bias = 10.0_wp
         call fit_values(ec, mol, grid, phi, q, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, q(1) > baseline(1) + 0.01_wp, more="distant samples did not receive greater weight")
      if (allocated(error)) return
      ec%reference_charges = [0.2_wp, -0.2_wp]
      ec%restraint_strength = 1.0e6_wp
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
   !>
   !> @param[out] error  test error
   subroutine check_adjoint(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(ec_charges_type) :: ec
      call adjoint_case(error, ec, 0)
   end subroutine check_adjoint

   !> Check penalty-kink and hard-bound adjoints
   !>
   !> @param[out] error  test error
   subroutine check_active_adjoint(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(ec_charges_type) :: ec
      ! TODO: cover lower-bound hits, negative penalty kinks, constraint release and the all-fixed
      ! multiplier branch of solve_fit; mutants of each survive the suite
      ec%charge_penalty = 0.0_wp
      ec%upper_bounds = [0.03_wp, 10.0_wp, 10.0_wp]
      call adjoint_case(error, ec, 1)
      if (allocated(error)) return
      deallocate (ec%upper_bounds)
      ec%charge_thresholds = [0.03_wp, 10.0_wp, 10.0_wp]
      ec%charge_penalty = 10.0_wp
      call adjoint_case(error, ec, 2)
   end subroutine check_active_adjoint

   !> Compare potential, nuclear, grid-position and volume-weight derivatives with finite differences
   !>
   !> @param[out] error  test error
   !> @param[in]  config configured term
   !> @param[in]  mode   0 smooth, 1 hard bound, 2 penalty kink
   subroutine adjoint_case(error, config, mode)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Configured term
      type(ec_charges_type), intent(in) :: config
      !> Fit mode
      integer, intent(in) :: mode
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
      integer :: i, j, kind, s
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
      if (mode /= 0) then
         call check(error, abs(q(1) - 0.03_wp) < 1.0e-9_wp, more="active constraint was not exercised")
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
   !>
   !> @param[out] error  test error
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
      ec%exclusion_radius = 100.0_wp
      call ec%update(mol, err)
      if (.not. allocated(err)) call ec%update_grid(grid, err)
      call check(error, allocated(err), more="fully excluded fit accepted")
      if (allocated(error)) return
      ec%exclusion_radius = 1.0_wp
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
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
   !>
   !> @param[out] error  test error
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
      call check(error, .not. response%next(), more="the potential pushed more than the EC potential adjoint")
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
   !>
   !> @param[out] error  test error
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

      ec%exclusion_radius = ieee_value(1.0_wp, ieee_quiet_nan)
      call ec%update(mol, err)
      call check(error, allocated(err), more="EC accepted a nonfinite exclusion radius")
      if (allocated(error)) return
      call check(error, err%message, "EC fit settings must be finite")
   end subroutine check_stale_guards

end module test_moz_potential_ec
