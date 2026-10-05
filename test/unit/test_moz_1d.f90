!> 1D MOZ model contract tests
module test_moz_1d
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type, atomic_charge_adjoint_response_type, response_accumulate
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
         new_unittest("getter_guards", check_getter_guards), &
         new_unittest("update_guards", check_update_guards), &
         new_unittest("atom_count_follows_update", check_atom_count_follows_update), &
         new_unittest("qat_charges_in_every_phase", check_qat_phases), &
         new_unittest("unsupported_coupling_modes", check_mode_failures)]
   end subroutine collect_moz_1d

   !> Drive a 1D MOZ model through the coupling protocol
   !>
   !> Mints a coupling, drives the atomic_charges request through the
   !> usual coupling protocol, then reports pending theory once every
   !> mandatory output is answered; a wrong-shape answer is rejected first
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

      associate (item => coupling%request())
         call check(error, item%is_missing("q"), more="wrong shape satisfied q")
      end associate
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
      call check(error, .not. associated(coupling), more="release left a coupling pointer")
   end subroutine check_moz_1d

   !> An existing coupling follows the model to a new atom count
   !>
   !> The coupling is minted for two atoms; after an update to three the
   !> next staging takes three charges and refuses two
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
      call check(error, .not. associated(coupling), more="release left a coupling pointer")
   end subroutine check_atom_count_follows_update

   !> The default "qat" source declares partial charges
   !>
   !> `q` pending in the energy, response and gradient phases
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
      call check(error, .not. associated(coupling), more="release left a coupling pointer")
   end subroutine check_qat_phases

   !> Unsupported electrostatic sources are refused
   !>
   !> "ec" waits for the radial grid, "multipoles" is 3D only and an unknown
   !> source is refused; each failure is named and leaves no coupling or
   !> staging behind
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
      call check(error, .not. associated(coupling), more="release left a coupling pointer")
   end subroutine check_mode_failures

   !> Getter guards precede pending theory and preserve caller outputs
   !>
   !> Every getter meets a foreign coupling, a model that is not updated, an
   !> unstaged coupling, missing charges, the wrong staging and, past every
   !> guard, the pending theory; each refusal is named and leaves the energy,
   !> the gradient and the seeded response untouched
   subroutine check_getter_guards(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model, foreign
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling, other
      type(response_type) :: response
      type(atomic_charge_adjoint_response_type) :: seed
      real(wp) :: energy, gradient(3, 2)
      integer :: phase, scenario
      !> Getter and scenario of the current case for diagnostics
      character(len=64) :: label
      character(len=8), parameter :: phases(3) = [character(len=8) :: "energy", "response", "gradient"]
      !> Refusal each scenario must name
      character(len=40), parameter :: reasons(6) = [character(len=40) :: &
         & "different model", "updated first", "not staged", "missing required outputs", &
         & "staged for the", "is not implemented"]

      call new_context(ctx, verbosity=0)
      call new_updated_model(ctx, model, err)
      if (.not. allocated(err)) call new_updated_model(ctx, foreign, err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call foreign%new_coupling(other, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call new(mol, [1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      seed%dg_dq = [4.0_wp, -2.0_wp]
      do phase = 1, 3
         do scenario = 1, 6
            label = trim(phases(phase))//" getter, case '"//trim(reasons(scenario))//"'"
            ! Every case starts from an updated model; the update ends every walk
            call model%update(mol, err)
            if (allocated(err)) then
               call test_failed(error, err%message)
               return
            end if
            select case (scenario)
            case (2)
               call model%invalidate()
            case (4, 6)
               select case (phase)
               case (1)
                  call model%prepare_energy(coupling, err)
               case (2)
                  call model%prepare_response(coupling, err)
               case default
                  call model%prepare_gradient(coupling, err)
               end select
            case (5)
               if (phase == 1) then
                  call model%prepare_response(coupling, err)
               else
                  call model%prepare_energy(coupling, err)
               end if
            case default
               ! Foreign and unstaged scenarios require no preparation
            end select
            if (allocated(err)) then
               call test_failed(error, err%message)
               return
            end if
            if (scenario == 6) then
               call check(error, coupling%next(), more=trim(label)//": charges not pending")
               if (allocated(error)) return
               call coupling%answer("q", [0.3_wp, -0.3_wp], err)
               if (allocated(err)) then
                  call test_failed(error, err%message)
                  return
               end if
            end if
            call response_accumulate(response, seed, err)
            if (allocated(err)) then
               call test_failed(error, err%message)
               return
            end if
            energy = 7.0_wp
            gradient = 3.0_wp
            if (scenario == 1) then
               call invoke_getter(model, other, phase, response, energy, gradient, err)
            else
               call invoke_getter(model, coupling, phase, response, energy, gradient, err)
            end if
            call check(error, allocated(err), more=trim(label)//": accepted")
            if (allocated(error)) return
            call check(error, index(err%message, trim(reasons(scenario))) > 0, &
               & more=trim(label)//": "//err%message)
            if (allocated(error)) return
            if (scenario == 6) then
               call check(error, index(err%message, "1D MOZ "//trim(phases(phase))) > 0, &
                  & more=trim(label)//": "//err%message)
               if (allocated(error)) return
            end if
            call check(error, energy == 7.0_wp .and. all(gradient == 3.0_wp), &
               & more=trim(label)//": the refusal wrote the energy or gradient")
            if (allocated(error)) return
            ! A rejected request, including pending theory, keeps the seeded response
            call check(error, response%next(), &
               & more=trim(label)//": the refusal cleared the response")
            if (allocated(error)) return
            select type (item => response%item())
            type is (atomic_charge_adjoint_response_type)
               call check(error, all(item%dg_dq == real(scenario, wp)*seed%dg_dq), &
                  & more=trim(label)//": the refusal changed the response")
            class default
               call test_failed(error, trim(label)//": the refusal replaced the response type")
            end select
            if (allocated(error)) return
            call check(error, .not. response%next(), &
               & more=trim(label)//": the refusal appended a response item")
            if (allocated(error)) return
            deallocate (err)
         end do
         ! Restart the seed accumulation for the next getter channel
         block
            type(response_type) :: empty
            response = empty
         end block
      end do
      call model%release_coupling(coupling)
      call foreign%release_coupling(other)
      call check(error, .not. associated(coupling) .and. .not. associated(other))
   end subroutine check_getter_guards

   !> Invoke one model getter for a common guard table
   !>
   !> @param[in,out] model Model under test
   !> @param[in,out] coupling Owned or foreign coupling
   !> @param[in] phase Getter index
   !> @param[in,out] response Seeded response
   !> @param[in,out] energy Seeded energy
   !> @param[in,out] gradient Seeded gradient
   !> @param[out] err Getter error
   subroutine invoke_getter(model, coupling, phase, response, energy, gradient, err)
      !> Model under test
      type(model_moz_1d_type), intent(inout) :: model
      !> Owned or foreign coupling
      type(coupling_type), intent(inout), target :: coupling
      !> Getter index
      integer, intent(in) :: phase
      !> Seeded response
      type(response_type), intent(inout) :: response
      !> Seeded energy
      real(wp), intent(inout) :: energy
      !> Seeded gradient
      real(wp), intent(inout) :: gradient(:, :)
      !> Getter error
      type(moist_error), allocatable, intent(out) :: err
      select case (phase)
      case (1)
         call model%get_energy(coupling, energy, err)
      case (2)
         call model%get_response(coupling, response, err)
      case default
         call model%get_gradient(coupling, response, gradient, err)
      end select
   end subroutine invoke_getter

   !> Failed updates reset count and status and end the staged walk
   subroutine check_update_guards(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling

      call check(error, model%atom_count(), 0)
      if (allocated(error)) return
      call check(error, .not. model%is_updated())
      if (allocated(error)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call model%update(mol, err)
      call check(error, allocated(err), more="unconstructed update accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Construct the 1D MOZ model") > 0)
      if (allocated(error)) return
      call new_context(ctx, verbosity=0)
      call new_updated_model(ctx, model, err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      ! q stays unanswered; the full pass rewinds, so a staged walk offers it again
      call check(error, coupling%next(), more="atomic_charges must be pending")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="atomic_charges is the only request")
      if (allocated(error)) return
      call new(mol, [integer ::], reshape([real(wp) ::], [3, 0]))
      call model%update(mol, err)
      call check(error, allocated(err), more="empty structure accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "at least one solute atom") > 0, more=err%message)
      if (allocated(error)) return
      call check(error, model%atom_count(), 0)
      if (allocated(error)) return
      call check(error, .not. model%is_updated(), more="failed update kept valid status")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="failed update kept the staged walk")
      if (allocated(error)) return
      call model%release_coupling(coupling)
      call check(error, .not. associated(coupling))
   end subroutine check_update_guards

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
