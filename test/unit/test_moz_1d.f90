!> 1D MOZ model contract tests
module test_moz_1d
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type
   implicit none(type, external)
   private
   public :: collect_moz_1d
contains
   !> Register 1D MOZ tests
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_moz_1d(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("contract", check_moz_1d), &
         new_unittest("atom_count_follows_update", check_atom_count_follows_update)]
   end subroutine collect_moz_1d

   !> Drive a 1D MOZ model through the coupling protocol
   !>
   !> Mints a coupling, drives the atomic_multipoles request through the
   !> usual coupling protocol, then reports pending theory once every
   !> mandatory output is answered; a wrong-shape answer is rejected first
   !>
   !> @param[out] error Test error
   subroutine check_moz_1d(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      real(wp) :: energy, gradient(3, 2)

      call new_context(ctx, verbosity=0)
      call new_moz_1d_model(model, ctx, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call new(mol, [1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, model%atom_count(), 2)
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
      call check(error, coupling%next(), more="atomic_multipoles must be pending")
      if (allocated(error)) return

      ! A wrong-shape answer is rejected and leaves the output missing
      call coupling%answer("q", [1.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return

      ! The correctly shaped answer (one charge per atom) completes the walk
      call coupling%answer("q", [0.3_wp, -0.3_wp], err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, .not. coupling%next(), more="atomic_multipoles is the only request")
      if (allocated(error)) return

      energy = 7.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err) .and. energy == 7.0_wp)
      if (allocated(error)) return
      call check(error, index(err%message, "1D MOZ energy") > 0)
      if (allocated(error)) return

      call model%prepare_gradient(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      gradient = 3.0_wp
      call model%get_gradient(coupling, response, gradient, err)
      call check(error, allocated(err) .and. all(gradient == 3.0_wp))
      if (allocated(error)) return
      call check(error, index(err%message, "1D MOZ gradient") > 0)
      if (allocated(error)) return

      call model%release_coupling(coupling)
   end subroutine check_moz_1d

   !> An existing coupling follows the model to a new atom count
   !>
   !> The coupling is minted for two atoms; after an update to three the
   !> next staging takes three charges and refuses two
   !>
   !> @param[out] error Test error
   subroutine check_atom_count_follows_update(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling

      call new_context(ctx, verbosity=0)
      call new_moz_1d_model(model, ctx, err)
      if (.not. allocated(err)) then
         call new(mol, [1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
         call model%update(mol, err)
      end if
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, coupling%next(), more="atomic_multipoles must be pending")
      if (allocated(error)) return
      call coupling%answer("q", [0.3_wp, -0.3_wp], err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if

      ! Same model and coupling, one more atom
      call new(mol, [8, 1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp, 0.0_wp, 0.0_wp, &
         & -0.5_wp, 1.7_wp, 0.0_wp], [3, 3]))
      call model%update(mol, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, coupling%next(), more="the update kept the two-atom answer")
      if (allocated(error)) return
      call coupling%answer("q", [0.3_wp, -0.3_wp], err)
      call check(error, allocated(err), more="two charges were accepted for three atoms")
      if (allocated(error)) return
      deallocate (err)
      call coupling%answer("q", [-0.6_wp, 0.3_wp, 0.3_wp], err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, .not. coupling%next(), more="three charges complete the walk")
      if (allocated(error)) return

      call model%release_coupling(coupling)
   end subroutine check_atom_count_follows_update

end module test_moz_1d
