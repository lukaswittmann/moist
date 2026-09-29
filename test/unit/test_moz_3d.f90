!> 3D MOZ model and grid ownership tests
module test_moz_3d
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_grid_3d
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type
   use moist_channels_fields, only: field_query_type
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type
   implicit none(type, external)
   private
   public :: collect_moz_3d
contains
   !> Register 3D MOZ tests
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_moz_3d(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("contract", check_moz_3d), &
         & new_unittest("owned_grid", test_grid_model)]
   end subroutine collect_moz_3d

   !> Drive a 3D MOZ model through the coupling protocol
   !>
   !> The model owns a copy of the configured grid, drives the
   !> gaussian_potential request through the usual coupling protocol, then
   !> reports pending theory once the mandatory output is answered; a
   !> wrong-shape answer is rejected first
   !>
   !> @param[out] error Test error
   subroutine check_moz_3d(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: template
      class(moist_math_grid_3d_type), pointer :: grid
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      real(wp), allocatable :: phi(:)
      real(wp) :: energy
      integer :: i
      call new_context(ctx, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      template%nx = 2
      template%ny = 2
      template%nz = 2
      template%dr = 0.5_wp
      call new_moz_3d_model(model, template, ctx, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      grid => model%grid
      call check(error, associated(grid) .and. model%is_updated())
      if (allocated(error)) return
      call check(error, grid%ngrid == 8 .and. model%atom_count() == 1)
      if (allocated(error)) return
      template%nx = 3
      call check(error, grid%ngrid == 8, more="model must own an independent grid copy")
      if (allocated(error)) return

      call model%new_coupling(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, coupling%next(), more="gaussian_potential must be pending")
      if (allocated(error)) return

      ! A wrong-shape answer is rejected and leaves the output missing
      call coupling%answer("phi", [1.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return

      ! The correctly shaped answer (one value per grid point) completes the walk
      allocate (phi(grid%ngrid))
      phi = [(0.1_wp*real(i, wp), i=1, grid%ngrid)]
      call coupling%answer("phi", phi, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, .not. coupling%next(), more="gaussian_potential is the only request")
      if (allocated(error)) return

      energy = 9.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err) .and. energy == 9.0_wp)
      if (allocated(error)) return
      call check(error, index(err%message, "3D MOZ energy") > 0)
      if (allocated(error)) return

      call model%prepare_response(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call model%get_response(coupling, response, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "3D MOZ response") > 0)
      if (allocated(error)) return

      call model%release_coupling(coupling)
   end subroutine check_moz_3d

   !> Concrete grids work through the 3D MOZ model and preserve copy ownership
   !>
   !> @param[out] error Test error
   subroutine test_grid_model(error)
      !> Borrowed typed model domain
      class(moist_math_grid_3d_type), pointer :: model_domain
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(moist_math_grid_3d_molecular_type) :: molecular
      type(model_moz_3d_type), target :: model
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(moist_math_grid_3d_type), allocatable :: domain
      type(field_query_type) :: query
      real(wp), allocatable :: original(:, :)
      integer :: kind

      call new_context(ctx, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_cartesian_grid_3d(cart, 4, 6, 8, 0.5_wp, error=err)
      call require_success(error, err)
      if (allocated(error)) return
      molecular%nrad = 8
      molecular%lebedev_degree = 5
      molecular%rmax = 5.0_wp
      call molecular%validate(err)
      call require_success(error, err)
      if (allocated(error)) return
      do kind = 1, 2
         if (kind == 1) then
            allocate (domain, source=cart)
         else
            allocate (domain, source=molecular)
         end if
         call new_moz_3d_model(model, domain, ctx, err)
         call require_success(error, err)
         if (allocated(error)) return
         call model%update(mol, err)
         call require_success(error, err)
         if (allocated(error)) return
         model_domain => model%grid
         call check(error, model_domain%natom, 1)
         if (allocated(error)) return
         call query%fetch("w")
         call model%list_fields(query)
         call check(error, query%found)
         if (allocated(error)) return
         model_domain => model%grid
         call check(error, size(query%rvals), model_domain%ngrid)
         if (allocated(error)) return
         model_domain => model%grid
         original = model_domain%xyz
         mol%xyz(1, 1) = mol%xyz(1, 1) + 0.25_wp
         call domain%update(mol, err)
         call require_success(error, err)
         if (allocated(error)) return
         model_domain => model%grid
         call check(error, all(model_domain%xyz == original), "model must own an independent grid copy")
         if (allocated(error)) return
         model_domain => model%grid
         select type (g => model_domain)
         type is (moist_math_grid_3d_cartesian_type)
            call check(error, g%ngrid, 4*6*8)
         type is (moist_math_grid_3d_molecular_type)
            call check(error, g%has_kgrid)
         class default
            call test_failed(error, "model did not retain the concrete grid domain")
         end select
         if (allocated(error)) return
         deallocate (domain)
      end do
   end subroutine test_grid_model

   !> Forward a library error into the test framework
   !>
   !> @param[out] error Test error
   !> @param[in] err Library error
   subroutine require_success(error, err)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable, intent(in) :: err
      if (allocated(err)) call test_failed(error, err%message)
   end subroutine require_success
end module test_moz_3d
