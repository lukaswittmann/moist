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
         new_unittest("atom_count_follows_update", check_atom_count_follows_update), &
         new_unittest("qat_charges_in_every_phase", check_qat_phases), &
         new_unittest("unsupported_coupling_modes", check_mode_failures)]
   end subroutine collect_moz_1d

   !> Drive a 1D MOZ model through the coupling protocol
   !>
   !> Mints a coupling, drives the atomic_charges request through the
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
      call check(error, coupling%next(), more="atomic_charges must be pending")
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
      call check(error, .not. coupling%next(), more="atomic_charges is the only request")
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
      call check(error, coupling%next(), more="atomic_charges must be pending")
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

   !> The default "qat" source declares the solute's partial charges, with
   !> `q` pending in the energy, response and gradient phases
   !>
   !> @param[out] error Test error
   subroutine check_qat_phases(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(coupling_type), pointer :: coupling
      !> Phase labels for diagnostics
      character(len=8), parameter :: phases(3) = [character(len=8) :: "energy", "response", "gradient"]
      integer :: phase

      call new_context(ctx, verbosity=0)
      call new_updated_model(ctx, model, err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, model%coupling_mode == "qat", more="the 1D default source is not 'qat'")
      if (allocated(error)) return
      ! Nothing is answered, so every phase finds q missing
      do phase = 1, size(phases)
         select case (phase)
         case (1)
            call model%prepare_energy(coupling, err)
         case (2)
            call model%prepare_response(coupling, err)
         case default
            call model%prepare_gradient(coupling, err)
         end select
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         call check(error, coupling%next(), more="no request pending in the "//trim(phases(phase))//" phase")
         if (allocated(error)) return
         associate (item => coupling%request())
            call check(error, item%name() == "atomic_charges" .and. item%is_missing("q"), &
               & more="q of atomic_charges is not pending in the "//trim(phases(phase))//" phase")
         end associate
         if (allocated(error)) return
         call check(error, .not. coupling%next(), more="atomic_charges is not the only request")
         if (allocated(error)) return
      end do

      call model%release_coupling(coupling)
   end subroutine check_qat_phases

   !> "ec" waits for the radial grid, "multipoles" is 3D only and an unknown
   !> source is refused; each failure is named and leaves no coupling or
   !> staging behind
   !>
   !> @param[out] error Test error
   subroutine check_mode_failures(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(coupling_type), pointer :: coupling
      !> Refused coupling sources
      character(len=16), parameter :: modes(3) = [character(len=16) :: "ec", "multipoles", "bogus"]
      !> Part of the refusal each source must name
      character(len=64), parameter :: reasons(3) = [character(len=64) :: &
         & "EC coupling requires the radial grid (pending grid refactor)", &
         & "coupling_mode 'multipoles' is not supported", &
         & "coupling_mode 'bogus' is not supported"]
      integer :: i

      call new_context(ctx, verbosity=0)
      call new_updated_model(ctx, model, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      do i = 1, size(modes)
         model%coupling_mode = modes(i)
         call model%new_coupling(coupling, err)
         call check(error, allocated(err) .and. .not. associated(coupling), &
            & more="coupling_mode '"//trim(modes(i))//"' was accepted")
         if (allocated(error)) return
         call check(error, index(err%message, "1D MOZ") > 0 .and. index(err%message, trim(reasons(i))) > 0, &
            & more=err%message)
         if (allocated(error)) return
         deallocate (err)
      end do

      ! A coupling staged under "qat" refuses a later unknown source and keeps no staging
      model%coupling_mode = "qat"
      call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      model%coupling_mode = "bogus"
      call model%prepare_response(coupling, err)
      call check(error, allocated(err), more="an unknown source was staged")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="the refused declaration left a staging behind")
      if (allocated(error)) return

      call model%release_coupling(coupling)
   end subroutine check_mode_failures

   !> Construct a 1D MOZ model and update it to two hydrogen atoms
   !>
   !> @param[in] ctx Run context, outlives the model
   !> @param[out] model Updated model
   !> @param[out] err Construction or update error
   subroutine new_updated_model(ctx, model, err)
      !> Run context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Updated model
      type(model_moz_1d_type), intent(out) :: model
      !> Construction or update error
      type(moist_error), allocatable, intent(out) :: err
      type(structure_type) :: mol
      call new_moz_1d_model(model, ctx, err)
      if (allocated(err)) return
      call new(mol, [1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call model%update(mol, err)
   end subroutine new_updated_model

end module test_moz_1d
