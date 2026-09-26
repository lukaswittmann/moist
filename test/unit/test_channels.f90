!> Raw primitive protocol, host walk and response accumulation regressions
module test_channels
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & coupling_request_type, coupling_registry_type, point_potential_request_type, &
      & gaussian_potential_request_type, gaussian_moment_request_type, moist_phase_energy, &
      & moist_phase_response, moist_phase_gradient, request_name_len, coupling_register, &
      & request_require, coupling_begin_registration, coupling_set_scope, &
      & coupling_snapshot, coupling_arm, coupling_invalidate, coupling_check_mandatory, &
      & coupling_make_view, coupling_close_view
   use moist_channels_response, only: response_type, response_item_type, &
      & potential_adjoint_response_type, density_response_type, &
      & gostshyp_amplitude_response_type, current_response_item, response_accumulate, &
      & response_clear
   use test_helpers, only: check_moist_error
   implicit none(type, external)
   private
   public :: collect_channels

   !> Diagnostic of an answer or query outside a `next()` window
   character(len=*), parameter :: no_current = "No current coupling request - call next() first"

contains

   !> Register protocol regressions
   subroutine collect_channels(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("registry_release", test_registry_release), &
         new_unittest("registry_release_tail_and_foreign", test_registry_release_tail), &
         new_unittest("shared_requests_and_compaction", test_shared_requests), &
         new_unittest("direct_potential_read", test_direct_potential_read), &
         new_unittest("partial_outputs_and_phase_reuse", test_partial), &
         new_unittest("one_output_per_answer", test_partial_copies), &
         new_unittest("snapshot_keeps_answers_unless_resized", test_snapshot_change), &
         new_unittest("staging_supersedes_the_walk", test_restaging), &
         new_unittest("rejected_answer_and_retry", test_rejected), &
         new_unittest("unknown_output_and_rank", test_unknown_output), &
         new_unittest("failed_preparation_leaves_no_staging", test_failed_preparation), &
         new_unittest("declaration_errors_are_returned", test_declaration_errors), &
         new_unittest("register_ignores_callers_answers", test_register_neutralises_copy), &
         new_unittest("scoped_registration_and_expired_view", test_scope), &
         new_unittest("point_and_gaussian_are_distinct", test_distinct_kinds), &
         new_unittest("energy_spatial_derivative", test_energy_derivative), &
         new_unittest("unrequired_registration_survives", test_empty_registration), &
         new_unittest("walk_unstaged_is_empty", test_walk_unstaged), &
         new_unittest("walk_one_visit_per_pass_in_order", test_walk_order), &
         new_unittest("walk_rewind_retries_rejected_output", test_walk_rewind), &
         new_unittest("walk_answered_pass_visits_nothing", test_walk_complete), &
         new_unittest("walk_early_exit_resumes", test_walk_resume), &
         new_unittest("walk_reset_by_staging_and_invalidation", test_walk_reset), &
         new_unittest("answer_without_current_request", test_answer_outside_window), &
         new_unittest("request_outside_window_is_placeholder", test_placeholder), &
         new_unittest("unknown_output_is_not_missing", test_unknown_not_missing), &
         new_unittest("staging_messages", test_staging_messages), &
         new_unittest("missing_outputs_named", test_missing_named), &
         new_unittest("register_without_begin", test_register_unprepared), &
         new_unittest("placeholder_request_not_registered", test_register_placeholder), &
         new_unittest("view_errors_named", test_view_errors), &
         new_unittest("response_walk_empty", test_response_walk_empty), &
         new_unittest("response_accumulate", test_response_accumulate), &
         new_unittest("response_walk_rewinds_and_resumes", test_response_walk_passes), &
         new_unittest("response_walk_reset_by_accumulate_and_clear", test_response_walk_reset), &
         new_unittest("response_item_is_a_copy", test_response_item_copy), &
         new_unittest("response_add_rejects_shape", test_response_add_rejects_shape), &
         new_unittest("response_density_every_rank", test_response_density_ranks), &
         new_unittest("response_add_rejects_vector_shape", test_response_vector_shape), &
         new_unittest("response_placeholder_cannot_accumulate", test_response_placeholder_add)]
   end subroutine collect_channels

   !* ================================================================================= *!
   !*                                      Fixtures                                     *!
   !* ================================================================================= *!

   !> A point request needing phi for energy and phi plus dphi_dr for the
   !> gradient, on two grid points, staged for the energy phase
   subroutine fixture(coupling, err)
      type(coupling_type), intent(out) :: coupling
      type(moist_error_type), allocatable, intent(out) :: err
      type(point_potential_request_type) :: point
      call coupling_begin_registration(coupling)
      call request_require(point, moist_phase_energy, "phi", err)
      if (allocated(err)) return
      call request_require(point, moist_phase_gradient, "phi", err)
      if (allocated(err)) return
      call request_require(point, moist_phase_gradient, "dphi_dr", err)
      if (allocated(err)) return
      call coupling_register(coupling, "potential", point, err)
      if (allocated(err)) return
      call coupling_snapshot(coupling, 2)
      call coupling_arm(coupling, moist_phase_energy, err)
   end subroutine fixture

   !> Three requests on two grid points, declared in the order point,
   !> gaussian, moments, each by its own scope; registered but not staged
   !>
   !> Every request needs one output in the energy phase; the Gaussian one
   !> also needs `dphi_dr` for the gradient
   subroutine three_requests(coupling, err)
      type(coupling_type), intent(out) :: coupling
      type(moist_error_type), allocatable, intent(out) :: err
      type(point_potential_request_type) :: point
      type(gaussian_potential_request_type) :: gaussian
      type(gaussian_moment_request_type) :: moments
      call coupling_begin_registration(coupling)
      call request_require(point, moist_phase_energy, "phi", err)
      if (allocated(err)) return
      call coupling_set_scope(coupling, 1)
      call coupling_register(coupling, "point", point, err)
      if (allocated(err)) return
      call request_require(gaussian, moist_phase_energy, "phi", err)
      if (allocated(err)) return
      call request_require(gaussian, moist_phase_gradient, "dphi_dr", err)
      if (allocated(err)) return
      call coupling_set_scope(coupling, 2)
      call coupling_register(coupling, "gaussian", gaussian, err)
      if (allocated(err)) return
      moments%width = [1.0_wp, 2.0_wp]
      call request_require(moments, moist_phase_energy, "gt", err)
      if (allocated(err)) return
      call coupling_set_scope(coupling, 3)
      call coupling_register(coupling, "moments", moments, err)
      if (allocated(err)) return
      call coupling_snapshot(coupling, 2)
   end subroutine three_requests

   !> Scientific name of the current request, "no_current_request" outside a `next()` window
   function current_name(coupling) result(name)
      type(coupling_type), intent(in) :: coupling
      character(len=request_name_len) :: name
      class(coupling_request_type), allocatable :: item
      item = coupling%request()
      name = item%name()
   end function current_name

   !> Answer the energy output of whichever request is current
   subroutine answer_energy_output(coupling, err)
      type(coupling_type), intent(inout) :: coupling
      type(moist_error_type), allocatable, intent(out) :: err
      select case (current_name(coupling))
      case ("gaussian_moments")
         call coupling%answer("gt", [1.0_wp, 2.0_wp], err)
      case default
         call coupling%answer("phi", [1.0_wp, 2.0_wp], err)
      end select
   end subroutine answer_energy_output

   !> Number of requests one full pass of the walk visits, answering none
   function count_visits(coupling) result(visits)
      type(coupling_type), intent(inout) :: coupling
      integer :: visits
      visits = 0
      do while (coupling%next())
         visits = visits + 1
      end do
   end function count_visits

   !* ================================================================================= *!
   !*                                 Declaration and data                              *!
   !* ================================================================================= *!

   subroutine test_partial(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call fixture(coupling, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      call check(error, coupling%next(), more="energy walk visits the potential")
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("phi") .and. .not. item%is_missing("dphi_dr"))
      end associate
      if (allocated(error)) return
      call coupling%answer("phi", [2.0_wp, 3.0_wp], err)
      call check_moist_error(error, err, "phi answer")
      if (allocated(error)) return
      ! The gradient phase keeps phi and asks for the derivative only
      call coupling_arm(coupling, moist_phase_gradient, err)
      call check(error, coupling%next(), more="gradient walk visits the potential")
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("dphi_dr") .and. .not. item%is_missing("phi"))
      end associate
      if (allocated(error)) return
      call coupling%answer("dphi_dr", spread([1.0_wp, 2.0_wp], 1, 3), err)
      call check_moist_error(error, err, "dphi_dr answer")
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_gradient, err)
      call check_moist_error(error, err, "partial fulfilment")
   end subroutine test_partial

   !> Outputs arrive one call at a time, in any order; a request copy is the
   !> state at the moment it was taken and never sees later answers
   subroutine test_partial_copies(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      class(coupling_request_type), allocatable :: first, second
      call fixture(coupling, err)
      call coupling_arm(coupling, moist_phase_gradient, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      first = coupling%request()
      call coupling%answer("dphi_dr", spread([3.0_wp, 4.0_wp], 1, 3), err)
      call check_moist_error(error, err, "dphi_dr first")
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_gradient, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "missing required outputs: phi") > 0)
      if (allocated(error)) return
      call coupling%answer("phi", [1.0_wp, 2.0_wp], err)
      call coupling_check_mandatory(coupling, moist_phase_gradient, err)
      call check_moist_error(error, err, "both outputs")
      if (allocated(error)) return
      second = coupling%request()
      call check(error, first%is_missing("phi") .and. first%is_missing("dphi_dr"), &
         & more="the earlier copy saw a later answer")
      if (allocated(error)) return
      call check(error, .not. second%is_missing("phi") .and. .not. second%is_missing("dphi_dr"))
   end subroutine test_partial_copies

   !> Only a new grid size drops answers; the coupling holds no grid values
   subroutine test_snapshot_change(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call fixture(coupling, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      call coupling%answer("phi", [1.0_wp, 2.0_wp], err)
      call check_moist_error(error, err, "phi answer")
      if (allocated(error)) return
      call coupling_snapshot(coupling, 2)
      call coupling_arm(coupling, moist_phase_gradient, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, .not. item%is_missing("phi") .and. item%is_missing("dphi_dr"))
      end associate
      if (allocated(error)) return
      call coupling_snapshot(coupling, 3)
      call coupling_arm(coupling, moist_phase_gradient, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("phi"), more="a resized grid kept an answer")
      end associate
   end subroutine test_snapshot_change

   !> A new staging leaves no current request, so an answer meant for the
   !> previous staging is refused and stores nothing
   subroutine test_restaging(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call fixture(coupling, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      call coupling_arm(coupling, moist_phase_energy, err)
      call coupling%answer("phi", [9.0_wp, 9.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, no_current)
      if (allocated(error)) return
      call check(error, coupling%next())
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("phi"), more="a refused answer was stored")
      end associate
      if (allocated(error)) return
      call coupling_invalidate(coupling)
      call coupling%answer("phi", [9.0_wp, 9.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, no_current)
   end subroutine test_restaging

   !> A rejected replacement discards the earlier answer until a valid retry
   subroutine test_rejected(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call fixture(coupling, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      call coupling%answer("phi", [1.0_wp, 2.0_wp], err)
      call check_moist_error(error, err, "phi answer")
      if (allocated(error)) return
      call coupling%answer("phi", [1.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "point_potential: phi shape mismatch") > 0)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "missing required outputs: phi") > 0)
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("phi"))
      end associate
      if (allocated(error)) return
      call coupling%answer("phi", [3.0_wp, 4.0_wp], err)
      call coupling_check_mandatory(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "retry")
      if (allocated(error)) return
      call coupling%answer("phi", [ieee_nan(), 0.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "phi must be finite") > 0)
   end subroutine test_rejected

   !> A wrong name or rank is refused by name and leaves other outputs intact
   subroutine test_unknown_output(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type), target :: coupling
      type(coupling_view_type) :: view
      real(wp), allocatable :: phi(:), dphi_dr(:, :)
      call fixture(coupling, err)
      call coupling_arm(coupling, moist_phase_gradient, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      call coupling%answer("phi", [1.0_wp, 2.0_wp], err)
      call check_moist_error(error, err, "phi answer")
      if (allocated(error)) return
      call coupling%answer("gt", [1.0_wp, 2.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "point_potential has no output 'gt'") > 0)
      if (allocated(error)) return
      call coupling%answer("dphi_dr", [1.0_wp, 2.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "dphi_dr rank mismatch") > 0)
      if (allocated(error)) return
      call coupling%answer("dphi_dr", reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp], [2, 2]), err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "dphi_dr shape mismatch") > 0)
      if (allocated(error)) return
      ! The valid phi survived all three refusals
      associate (item => coupling%request())
         call check(error, item%is_missing("dphi_dr") .and. .not. item%is_missing("phi"))
      end associate
      if (allocated(error)) return
      call coupling%answer("dphi_dr", spread([1.0_wp, 2.0_wp], 1, 3), err)
      call check_moist_error(error, err, "dphi_dr answer")
      if (allocated(error)) return
      call coupling_make_view(coupling, 1, view)
      call view%read("potential", "phi", dphi_dr, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "rank mismatch") > 0)
      if (allocated(error)) return
      call view%read("potential", "phi", phi, err)
      call check_moist_error(error, err, "phi read")
      if (allocated(error)) return
      call check(error, maxval(abs(phi - [1.0_wp, 2.0_wp])), 0.0_wp, thr=epsilon(1.0_wp))
      call coupling_close_view(coupling)
   end subroutine test_unknown_output

   !> A preparation that fails after registration leaves no staging behind but
   !> keeps the answers for the next successful one
   subroutine test_failed_preparation(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      type(point_potential_request_type) :: point
      call fixture(coupling, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      call coupling%answer("phi", [1.0_wp, 2.0_wp], err)
      call check_moist_error(error, err, "phi answer")
      if (allocated(error)) return
      ! The new pass fails on a registration the collection refuses
      call coupling_begin_registration(coupling)
      call request_require(point, moist_phase_energy, "phi", err)
      call coupling_register(coupling, "", point, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="a failed preparation left a walk")
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "is not staged") > 0)
      if (allocated(error)) return
      call coupling%answer("phi", [1.0_wp, 2.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, no_current)
      if (allocated(error)) return
      ! A correct pass finds phi still answered
      call coupling_begin_registration(coupling)
      call coupling_register(coupling, "potential", point, err)
      call coupling_snapshot(coupling, 2)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "re-preparation")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="phi was lost by the failed preparation")
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "retained phi")
   end subroutine test_failed_preparation

   !> Declaration mistakes come back as errors
   subroutine test_declaration_errors(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      type(point_potential_request_type) :: point
      call request_require(point, moist_phase_energy, "chi", err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "point_potential has no output 'chi'") > 0)
      if (allocated(error)) return
      call request_require(point, 0, "phi", err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "invalid coupling phase") > 0)
      if (allocated(error)) return
      ! A requirement under a false condition registers nothing to answer
      call request_require(point, moist_phase_energy, "phi", .false., err)
      call check_moist_error(error, err, "conditional requirement")
      if (allocated(error)) return
      call coupling_begin_registration(coupling)
      call coupling_register(coupling, "potential", point, err)
      call check_moist_error(error, err, "registration")
      if (allocated(error)) return
      call coupling_register(coupling, repeat("x", request_name_len + 1), point, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "Invalid local request name") > 0)
      if (allocated(error)) return
      call coupling_snapshot(coupling, 1)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check(error, .not. coupling%next(), more="a false condition still required phi")
      if (allocated(error)) return
      call coupling_arm(coupling, 4, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "Invalid coupling phase") > 0)
   end subroutine test_declaration_errors

   !> A registered copy contributes its declaration only, never its answers
   subroutine test_register_neutralises_copy(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: first
      type(coupling_type), target :: second
      type(coupling_view_type) :: view
      class(coupling_request_type), allocatable :: item
      real(wp), allocatable :: phi(:)
      call fixture(first, err)
      call check(error, first%next())
      if (allocated(error)) return
      call first%answer("phi", [1.0_wp, 2.0_wp], err)
      call check_moist_error(error, err, "phi answer")
      if (allocated(error)) return
      item = first%request()
      call check(error, .not. item%is_missing("phi"), more="the copy should carry the answer state")
      if (allocated(error)) return
      call coupling_begin_registration(second)
      call coupling_register(second, "potential", item, err)
      call coupling_snapshot(second, 2)
      call coupling_arm(second, moist_phase_energy, err)
      call check_moist_error(error, err, "registration of a copy")
      if (allocated(error)) return
      call check(error, count_visits(second), 1)
      if (allocated(error)) return
      call check(error, second%next())
      if (allocated(error)) return
      associate (copy => second%request())
         call check(error, copy%is_missing("phi"), more="the registration carried an answer over")
      end associate
      if (allocated(error)) return
      call coupling_make_view(second, 1, view)
      call view%read("potential", "phi", phi, err)
      call coupling_close_view(second)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "missing") > 0)
   end subroutine test_register_neutralises_copy

   !> Moment requests with different widths stay apart; each scope reads its own
   !> answer through a view that expires with the component call
   subroutine test_scope(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type), target :: coupling
      type(coupling_view_type) :: view, captured
      type(gaussian_moment_request_type) :: moments
      real(wp), allocatable :: gt(:)
      integer :: i, visits
      call coupling_begin_registration(coupling)
      call request_require(moments, moist_phase_energy, "gt", err)
      do i = 1, 2
         moments%width = [real(i, wp)]
         call coupling_set_scope(coupling, i)
         call coupling_register(coupling, "moments", moments, err)
      end do
      call coupling_snapshot(coupling, 1)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "staging")
      if (allocated(error)) return
      visits = 0
      do while (coupling%next())
         visits = visits + 1
         select type (item => coupling%request())
         type is (gaussian_moment_request_type)
            call coupling%answer("gt", item%width, err)
         end select
         call check_moist_error(error, err, "moment answer")
         if (allocated(error)) return
      end do
      call check(error, visits, 2, more="different widths must not share one request")
      if (allocated(error)) return
      call coupling_make_view(coupling, 2, view)
      captured = view
      call view%read("moments", "gt", gt, err)
      call check_moist_error(error, err, "scoped read")
      if (allocated(error)) return
      call check(error, gt(1), 2.0_wp)
      if (allocated(error)) return
      gt = -9.0_wp
      call view%read("moments", "gt", gt, err)
      call check_moist_error(error, err, "independent moment output")
      if (allocated(error)) return
      call check(error, gt(1), 2.0_wp)
      if (allocated(error)) return
      call coupling_close_view(coupling)
      call captured%read("moments", "gt", gt, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "Expired") > 0)
   end subroutine test_scope

   !> Point and Gaussian potentials are different calculations, never shared
   subroutine test_distinct_kinds(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      type(point_potential_request_type) :: point
      type(gaussian_potential_request_type) :: gaussian
      call coupling_begin_registration(coupling)
      call coupling_snapshot(coupling, 0)
      call request_require(point, moist_phase_energy, "phi", err)
      call request_require(gaussian, moist_phase_energy, "phi", err)
      call coupling_register(coupling, "point", point, err)
      call coupling_register(coupling, "gaussian", gaussian, err)
      call check_moist_error(error, err, "registration")
      if (allocated(error)) return
      call coupling_snapshot(coupling, 1)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check(error, count_visits(coupling), 2)
   end subroutine test_distinct_kinds

   !> An energy-phase requirement on a derivative output is honoured
   subroutine test_energy_derivative(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      type(point_potential_request_type) :: point
      call coupling_begin_registration(coupling)
      call request_require(point, moist_phase_energy, "dphi_dr", err)
      call coupling_register(coupling, "potential", point, err)
      call coupling_snapshot(coupling, 1)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("dphi_dr") .and. .not. item%is_missing("phi"))
      end associate
   end subroutine test_energy_derivative

   !> A registration without requirements survives: nothing to answer, but a
   !> component still resolves it and learns that the output is missing
   subroutine test_empty_registration(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type), target :: coupling
      type(coupling_view_type) :: view
      type(point_potential_request_type) :: point
      real(wp), allocatable :: phi(:)
      call coupling_begin_registration(coupling)
      call coupling_register(coupling, "unused", point, err)
      call coupling_snapshot(coupling, 0)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "staging")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="an unrequired output was walked")
      if (allocated(error)) return
      call coupling_make_view(coupling, 1, view)
      call view%check_mandatory(view%phase, err)
      call check_moist_error(error, err, "optional registration survives")
      if (allocated(error)) return
      call view%read("unused", "phi", phi, err)
      call coupling_close_view(coupling)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "missing required output phi") > 0)
   end subroutine test_empty_registration

   !> A coupling that never began a registration pass still accepts one, and
   !> an empty local name is refused like an overlong one
   !>
   !> @param[out] error Test failure
   subroutine test_register_unprepared(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Collection never passed to coupling_begin_registration
      type(coupling_type), target :: coupling
      !> Scoped read interface
      type(coupling_view_type) :: view
      !> Calculation registered on it
      type(point_potential_request_type) :: point
      !> Answer read back through the view
      real(wp), allocatable :: phi(:)
      call request_require(point, moist_phase_energy, "phi", err)
      call check_moist_error(error, err, "requirement")
      if (allocated(error)) return
      call coupling_register(coupling, "", point, err)
      call check(error, allocated(err), more="an empty local name was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Invalid local request name") > 0)
      if (allocated(error)) return
      call coupling_register(coupling, "potential", point, err)
      call check_moist_error(error, err, "registration without a begun pass")
      if (allocated(error)) return
      call coupling_snapshot(coupling, 2)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "staging")
      if (allocated(error)) return
      call check(error, count_visits(coupling), 1, more="the refused name left a request behind")
      if (allocated(error)) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call coupling%answer("phi", [5.0_wp, 6.0_wp], err)
      call check_moist_error(error, err, "phi answer")
      if (allocated(error)) return
      call coupling_make_view(coupling, 1, view)
      call view%read("potential", "phi", phi, err)
      call coupling_close_view(coupling)
      call check_moist_error(error, err, "phi read")
      if (allocated(error)) return
      call check(error, all(phi == [5.0_wp, 6.0_wp]), more="the read returned a different answer")
   end subroutine test_register_unprepared

   !> The placeholder returned outside a `next()` window cannot be registered,
   !> and the refusal leaves neither a request nor a local name behind
   !>
   !> @param[out] error Test failure
   subroutine test_register_placeholder(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Collection receiving the placeholder
      type(coupling_type), target :: coupling
      !> Scoped read interface
      type(coupling_view_type) :: view
      !> Placeholder request
      class(coupling_request_type), allocatable :: placeholder
      !> Read buffer
      real(wp), allocatable :: phi(:)
      call coupling_begin_registration(coupling)
      placeholder = coupling%request()
      call check(error, placeholder%name(), "no_current_request")
      if (allocated(error)) return
      call coupling_register(coupling, "ghost", placeholder, err)
      call check(error, allocated(err), more="the placeholder was registered")
      if (allocated(error)) return
      call check(error, index(err%message, "placeholder request 'no_current_request'") > 0 &
         & .and. index(err%message, "cannot be declared") > 0)
      if (allocated(error)) return
      call coupling_snapshot(coupling, 1)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "staging")
      if (allocated(error)) return
      call check(error, count_visits(coupling), 0, more="the refused placeholder left a request")
      if (allocated(error)) return
      call coupling_make_view(coupling, 1, view)
      call view%read("ghost", "phi", phi, err)
      call coupling_close_view(coupling)
      call check(error, allocated(err), more="the refused placeholder left a local name")
      if (allocated(error)) return
      call check(error, index(err%message, "did not register request 'ghost'") > 0)
   end subroutine test_register_placeholder

   !> Component reads name what went wrong: an unbound view, an invalid phase,
   !> a name another scope registered and an output the request does not declare
   !>
   !> @param[out] error Test failure
   subroutine test_view_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Collection with one request per scope
      type(coupling_type), target :: coupling
      !> View bound to scope 1, and one never bound
      type(coupling_view_type) :: view, unbound
      !> Read buffer
      real(wp), allocatable :: phi(:)
      call unbound%check_mandatory(moist_phase_energy, err)
      call check(error, allocated(err), more="an unbound view passed its check")
      if (allocated(error)) return
      call check(error, index(err%message, "Unbound component coupling view") > 0)
      if (allocated(error)) return
      call unbound%read("point", "phi", phi, err)
      call check(error, allocated(err), more="an unbound view was read")
      if (allocated(error)) return
      call check(error, index(err%message, "Unbound component coupling view") > 0)
      if (allocated(error)) return
      call three_requests(coupling, err)
      call check_moist_error(error, err, "three requests")
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, 0, err)
      call check(error, allocated(err), more="phase 0 passed the completeness check")
      if (allocated(error)) return
      call check(error, index(err%message, "Invalid coupling phase") > 0)
      if (allocated(error)) return
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "staging")
      if (allocated(error)) return
      do while (coupling%next())
         call answer_energy_output(coupling, err)
         call check_moist_error(error, err, "energy answer")
         if (allocated(error)) return
      end do
      call coupling_make_view(coupling, 1, view)
      call view%check_mandatory(4, err)
      call check(error, allocated(err), more="phase 4 passed the view check")
      if (allocated(error)) return
      call check(error, index(err%message, "Invalid coupling phase") > 0)
      if (allocated(error)) return
      call view%check_mandatory(moist_phase_energy, err)
      call check_moist_error(error, err, "energy check through the view")
      if (allocated(error)) return
      ! "gaussian" is registered, but by scope 2
      call view%read("gaussian", "phi", phi, err)
      call check(error, allocated(err), more="scope 1 read another scope's request")
      if (allocated(error)) return
      call check(error, index(err%message, "Component did not register request 'gaussian'") > 0)
      if (allocated(error)) return
      call view%read("point", "bogus", phi, err)
      call check(error, allocated(err), more="an undeclared output was read")
      if (allocated(error)) return
      call check(error, index(err%message, "point_potential: unknown output 'bogus'") > 0)
      if (allocated(error)) return
      call view%read("point", "phi", phi, err)
      call coupling_close_view(coupling)
      call check_moist_error(error, err, "phi read")
      if (allocated(error)) return
      call check(error, all(phi == [1.0_wp, 2.0_wp]), more="the read returned a different answer")
   end subroutine test_view_errors

   !* ================================================================================= *!
   !*                                     Host walk                                     *!
   !* ================================================================================= *!

   !> A coupling never staged has nothing to visit, and that is no error
   subroutine test_walk_unstaged(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: empty, registered
      call check(error, .not. empty%next(), more="a fresh coupling was walked")
      if (allocated(error)) return
      call three_requests(registered, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      call check(error, .not. registered%next(), more="an unstaged coupling was walked")
      if (allocated(error)) return
      call check(error, current_name(registered), "no_current_request")
   end subroutine test_walk_unstaged

   !> One pass visits every request with something missing once, in
   !> declaration order, then ends
   subroutine test_walk_order(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      character(len=request_name_len) :: seen(4)
      integer :: visits
      call three_requests(coupling, err)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      visits = 0
      seen = ""
      do while (coupling%next())
         visits = visits + 1
         if (visits > size(seen)) exit
         seen(visits) = current_name(coupling)
      end do
      call check(error, visits, 3)
      if (allocated(error)) return
      call check(error, seen(1), "point_potential")
      if (allocated(error)) return
      call check(error, seen(2), "gaussian_potential")
      if (allocated(error)) return
      call check(error, seen(3), "gaussian_moments")
      if (allocated(error)) return
      ! The gradient phase skips the requests with nothing new to supply
      call coupling_arm(coupling, moist_phase_gradient, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "gaussian_potential")
      if (allocated(error)) return
      call check(error, .not. coupling%next())
   end subroutine test_walk_order

   !> The end of a pass rewinds: a second loop retries exactly what is still
   !> missing, here an answer that was rejected
   subroutine test_walk_rewind(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      integer :: visits
      call three_requests(coupling, err)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      do while (coupling%next())
         if (current_name(coupling) == "gaussian_potential") then
            call coupling%answer("phi", [1.0_wp, 2.0_wp, 3.0_wp], err)
            call check(error, allocated(err), more="a mis-sized answer was accepted")
            if (allocated(error)) return
            deallocate (err)
         else
            call answer_energy_output(coupling, err)
            call check_moist_error(error, err, "energy answer")
            if (allocated(error)) return
         end if
      end do
      call coupling_check_mandatory(coupling, moist_phase_energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "gaussian_potential: missing required outputs: phi") > 0)
      if (allocated(error)) return
      visits = 0
      do while (coupling%next())
         visits = visits + 1
         call check(error, current_name(coupling), "gaussian_potential")
         if (allocated(error)) return
         call answer_energy_output(coupling, err)
         call check_moist_error(error, err, "retry")
         if (allocated(error)) return
      end do
      call check(error, visits, 1, more="the retry pass must visit only the rejected request")
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "complete after the retry")
   end subroutine test_walk_rewind

   !> Once everything is answered, every further pass is empty
   subroutine test_walk_complete(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call three_requests(coupling, err)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      do while (coupling%next())
         call answer_energy_output(coupling, err)
         call check_moist_error(error, err, "energy answer")
         if (allocated(error)) return
      end do
      call check(error, count_visits(coupling), 0)
      if (allocated(error)) return
      call check(error, count_visits(coupling), 0)
   end subroutine test_walk_complete

   !> Leaving a pass early resumes it at the next call; its end still rewinds
   subroutine test_walk_resume(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call three_requests(coupling, err)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      do while (coupling%next())
         exit
      end do
      call check(error, current_name(coupling), "point_potential")
      if (allocated(error)) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "gaussian_potential", more="the pass did not resume")
      if (allocated(error)) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "gaussian_moments")
      if (allocated(error)) return
      call check(error, .not. coupling%next())
      if (allocated(error)) return
      ! Nothing was answered, so the next pass starts over from the first request
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "point_potential")
   end subroutine test_walk_resume

   !> Staging restarts a half-walked pass; invalidation ends it
   subroutine test_walk_reset(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call three_requests(coupling, err)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      ! One step per statement: `next` moves the cursor, so it never shares an expression
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, coupling%next(), more="two steps into the pass")
      if (allocated(error)) return
      call check(error, current_name(coupling), "gaussian_potential")
      if (allocated(error)) return
      call coupling_arm(coupling, moist_phase_energy, err)
      call check(error, current_name(coupling), "no_current_request", more="staging kept a current request")
      if (allocated(error)) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "point_potential", more="staging did not restart the walk")
      if (allocated(error)) return
      ! A new declaration pass, as every `prepare_*` runs, also restarts it
      call check(error, coupling%next())
      if (allocated(error)) return
      call coupling_begin_registration(coupling)
      call check(error, current_name(coupling), "no_current_request")
      if (allocated(error)) return
      call three_requests(coupling, err)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "point_potential")
      if (allocated(error)) return
      ! Invalidation, as every model update pushes it, ends the walk for good
      call coupling_invalidate(coupling)
      call check(error, current_name(coupling), "no_current_request")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="an invalidated coupling was walked")
   end subroutine test_walk_reset

   !> An answer outside a `next()` window is refused with the documented text
   subroutine test_answer_outside_window(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call fixture(coupling, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      ! Staged, but no step taken yet
      call coupling%answer("phi", [1.0_wp, 2.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, no_current)
      if (allocated(error)) return
      ! After the pass has ended
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, .not. coupling%next())
      if (allocated(error)) return
      call coupling%answer("phi", spread([1.0_wp, 2.0_wp], 1, 1), err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, no_current)
      if (allocated(error)) return
      call coupling%answer("phi", spread(spread([1.0_wp, 2.0_wp], 1, 1), 1, 1), err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, no_current)
      if (allocated(error)) return
      ! The refusals stored nothing
      call check(error, coupling%next())
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("phi"))
      end associate
   end subroutine test_answer_outside_window

   !> Outside a `next()` window the request is a placeholder with no outputs
   subroutine test_placeholder(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      logical :: matched
      call fixture(coupling, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      matched = .false.
      select type (request => coupling%request())
      type is (point_potential_request_type)
         matched = .true.
      class default
         call check(error, request%name(), "no_current_request")
         if (allocated(error)) return
         call check(error, .not. request%is_missing("phi"))
         if (allocated(error)) return
      end select
      call check(error, .not. matched, more="the placeholder matched a request kind")
      if (allocated(error)) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "point_potential")
      if (allocated(error)) return
      call check(error, .not. coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "no_current_request", more="the end of a pass kept a request")
   end subroutine test_placeholder

   !> An output the request does not declare is never missing, but answering
   !> it is an error
   subroutine test_unknown_not_missing(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call fixture(coupling, err)
      call check(error, coupling%next())
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, .not. item%is_missing("bogus") .and. item%is_missing("phi"))
      end associate
      if (allocated(error)) return
      call coupling%answer("bogus", [1.0_wp, 2.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "point_potential has no output 'bogus'")
   end subroutine test_unknown_not_missing

   !> The staging check names the phase it wants and the one it found
   subroutine test_staging_messages(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call coupling_check_mandatory(coupling, moist_phase_energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "get_energy requires a coupling staged by prepare_energy; "// &
         & "this coupling is not staged")
      if (allocated(error)) return
      call fixture(coupling, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_gradient, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "get_gradient requires a coupling staged by prepare_gradient; "// &
         & "this coupling is staged for the energy phase")
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_response, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "get_response requires a coupling staged by prepare_response; "// &
         & "this coupling is staged for the energy phase")
   end subroutine test_staging_messages

   !> The completeness check lists every missing output of the request by name
   subroutine test_missing_named(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(coupling_type) :: coupling
      call fixture(coupling, err)
      call coupling_check_mandatory(coupling, moist_phase_energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "point_potential: missing required outputs: phi")
      if (allocated(error)) return
      call coupling_arm(coupling, moist_phase_gradient, err)
      call coupling_check_mandatory(coupling, moist_phase_gradient, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "point_potential: missing required outputs: phi, dphi_dr")
   end subroutine test_missing_named

   !* ================================================================================= *!
   !*                                  Response items                                   *!
   !* ================================================================================= *!

   !> Names of the items one full pass of the response walk visits, comma separated
   function walk_names(response) result(names)
      type(response_type), intent(inout) :: response
      character(len=:), allocatable :: names
      class(response_item_type), allocatable :: item
      names = ""
      do while (response%next())
         item = response%item()
         names = names//trim(item%name())//","
      end do
   end function walk_names

   !> Name of what `item()` returns now, "no_current_item" outside a pass
   function item_name(response) result(name)
      type(response_type), intent(in) :: response
      character(len=:), allocatable :: name
      class(response_item_type), allocatable :: item
      item = response%item()
      name = trim(item%name())
   end function item_name

   !> A potential adjoint accumulated before a density item
   subroutine two_items(response, err)
      type(response_type), intent(inout) :: response
      type(moist_error_type), allocatable, intent(out) :: err
      type(potential_adjoint_response_type) :: charge
      type(density_response_type) :: density
      allocate (charge%w_phi, source=[1.0_wp, 2.0_wp])
      call response_accumulate(response, charge, err)
      if (allocated(err)) return
      density%w_rho = [3.0_wp, 4.0_wp]
      call response_accumulate(response, density, err)
   end subroutine two_items

   !> A never-filled and a cleared response have nothing to visit; outside a
   !> pass `item()` is the "no_current_item" placeholder and the checked form says so
   subroutine test_response_walk_empty(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(response_type) :: response
      class(response_item_type), allocatable :: item
      call check(error, .not. response%next(), more="a never-filled response has an item")
      if (allocated(error)) return
      call check(error, item_name(response), "no_current_item")
      if (allocated(error)) return
      call current_response_item(response, item, err)
      call check(error, allocated(err), more="no item is current, yet no error")
      if (allocated(error)) return
      call check(error, err%message, "No current response item - call next() first")
      if (allocated(error)) return
      call response_clear(response)
      call check(error, .not. response%next(), more="a cleared response has an item")
   end subroutine test_response_walk_empty

   !> Accumulation keeps one item per kind in first-seen order: a repeated kind
   !> sums, and a partially filled item contributes only what it carries
   subroutine test_response_accumulate(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(response_type) :: response
      type(potential_adjoint_response_type) :: charge
      type(gostshyp_amplitude_response_type) :: amplitude
      type(density_response_type) :: density

      allocate (charge%w_phi, source=[1.0_wp, 2.0_wp, 3.0_wp])
      call response_accumulate(response, charge, err)
      call check_moist_error(error, err, "first accumulate")
      if (allocated(error)) return
      call check(error, walk_names(response), "potential_adjoint,")
      if (allocated(error)) return

      charge%w_phi = [1.0_wp, 1.0_wp, 1.0_wp]
      call response_accumulate(response, charge, err)
      call check_moist_error(error, err, "second accumulate")
      if (allocated(error)) return
      call check(error, walk_names(response), "potential_adjoint,", more="summing added an item")
      if (allocated(error)) return

      amplitude%w_overlap = [0.5_wp, 0.5_wp]
      call response_accumulate(response, amplitude, err)
      call check_moist_error(error, err, "amplitude accumulate")
      if (allocated(error)) return
      amplitude%w_normal_deriv = [7.0_wp, 8.0_wp]
      call response_accumulate(response, amplitude, err)
      call check_moist_error(error, err, "amplitude accumulate again")
      if (allocated(error)) return
      density%w_rho = [9.0_wp]
      call response_accumulate(response, density, err)
      call check_moist_error(error, err, "density accumulate")
      if (allocated(error)) return
      call check(error, walk_names(response), "potential_adjoint,gostshyp_amplitude,density,")
      if (allocated(error)) return

      ! One pass contracting each kind, as a host does
      do while (response%next())
         select type (item => response%item())
         type is (potential_adjoint_response_type)
            call check(error, all(item%w_phi == [2.0_wp, 3.0_wp, 4.0_wp]), "two accumulations sum")
         type is (gostshyp_amplitude_response_type)
            call check(error, all(item%w_overlap == 1.0_wp) &
               & .and. all(item%w_normal_deriv == [7.0_wp, 8.0_wp]), &
               & "overlap summed twice, normal derivative copied once")
         type is (density_response_type)
            call check(error, allocated(item%w_rho) .and. .not. allocated(item%w_grad_rho), &
               & "density carries only what was accumulated")
         class default
            call test_failed(error, "unexpected response item "//item%name())
         end select
         if (allocated(error)) return
      end do
   end subroutine test_response_accumulate

   !> A pass visits every item once and rewinds; a pass left early resumes at
   !> the next call
   subroutine test_response_walk_passes(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(response_type) :: response
      call two_items(response, err)
      call check_moist_error(error, err, "two items")
      if (allocated(error)) return
      call check(error, walk_names(response), "potential_adjoint,density,")
      if (allocated(error)) return
      call check(error, walk_names(response), "potential_adjoint,density,", &
         & more="the finished pass did not rewind")
      if (allocated(error)) return
      call check(error, response%next())
      if (allocated(error)) return
      call check(error, item_name(response), "potential_adjoint")
      if (allocated(error)) return
      call check(error, walk_names(response), "density,", more="a pass left early did not resume")
      if (allocated(error)) return
      call check(error, walk_names(response), "potential_adjoint,density,")
   end subroutine test_response_walk_passes

   !> `accumulate` and `clear` reset the cursor, so a model call always starts
   !> a fresh walk
   subroutine test_response_walk_reset(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(response_type) :: response
      type(gostshyp_amplitude_response_type) :: amplitude
      call two_items(response, err)
      call check_moist_error(error, err, "two items")
      if (allocated(error)) return
      call check(error, response%next())
      if (allocated(error)) return
      amplitude%w_overlap = [5.0_wp, 6.0_wp]
      amplitude%w_normal_deriv = [7.0_wp, 8.0_wp]
      call response_accumulate(response, amplitude, err)
      call check_moist_error(error, err, "amplitude accumulate")
      if (allocated(error)) return
      call check(error, item_name(response), "no_current_item", more="accumulate kept the current item")
      if (allocated(error)) return
      call check(error, walk_names(response), "potential_adjoint,density,gostshyp_amplitude,")
      if (allocated(error)) return
      call check(error, response%next())
      if (allocated(error)) return
      call check(error, response%next())
      if (allocated(error)) return
      call response_clear(response)
      call check(error, item_name(response), "no_current_item", more="clear kept the current item")
      if (allocated(error)) return
      call check(error, .not. response%next(), more="a cleared response has an item")
   end subroutine test_response_walk_reset

   !> `item()` is a copy: editing it leaves the response unchanged, and past the
   !> last item it is the placeholder again
   subroutine test_response_item_copy(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(response_type) :: response
      type(potential_adjoint_response_type) :: charge
      class(response_item_type), allocatable :: copy
      allocate (charge%w_phi, source=[1.0_wp, 2.0_wp])
      call response_accumulate(response, charge, err)
      call check_moist_error(error, err, "accumulate")
      if (allocated(error)) return
      call check(error, response%next())
      if (allocated(error)) return
      copy = response%item()
      select type (copy)
      type is (potential_adjoint_response_type)
         copy%w_phi(1) = 99.0_wp
      class default
         call test_failed(error, "the current item is not the potential adjoint")
      end select
      if (allocated(error)) return
      select type (item => response%item())
      type is (potential_adjoint_response_type)
         call check(error, all(item%w_phi == [1.0_wp, 2.0_wp]), "editing the copy changed the response")
      end select
      if (allocated(error)) return
      call current_response_item(response, copy, err)
      call check_moist_error(error, err, "current_response_item inside the pass")
      if (allocated(error)) return
      call check(error, trim(copy%name()), "potential_adjoint")
      if (allocated(error)) return
      call check(error, .not. response%next(), more="the pass did not end after the only item")
      if (allocated(error)) return
      call check(error, item_name(response), "no_current_item", more="an item is current after the pass")
   end subroutine test_response_item_copy

   !> Accumulating a differently shaped array is an error naming item and array
   subroutine test_response_add_rejects_shape(error)
      type(error_type), allocatable, intent(out) :: error
      type(response_type) :: response
      type(moist_error_type), allocatable :: err
      type(density_response_type) :: density

      allocate (density%w_grad_rho(3, 2))
      density%w_grad_rho = reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp], [3, 2])
      call response_accumulate(response, density, err)
      call check_moist_error(error, err, "first accumulate")
      if (allocated(error)) return
      deallocate (density%w_grad_rho)
      allocate (density%w_grad_rho(3, 1))
      density%w_grad_rho = reshape([1.0_wp, 2.0_wp, 3.0_wp], [3, 1])
      call response_accumulate(response, density, err)
      call check(error, allocated(err), "shape mismatch is an error")
      if (allocated(error)) return
      call check(error, index(err%message, "'density'") > 0 &
         & .and. index(err%message, "'w_grad_rho'") > 0, "error names item and array")
      if (allocated(error)) return
      call check(error, response%next())
      if (allocated(error)) return
      select type (item => response%item())
      type is (density_response_type)
         call check(error, all(shape(item%w_grad_rho) == [3, 2]), "stored item unchanged")
      class default
         call test_failed(error, "the stored item is not the density")
      end select

   end subroutine test_response_add_rejects_shape

   !> Density weights of every rank accumulate independently: an array absent
   !> on the stored item is copied, one absent on the new item contributes
   !> nothing, and a rank-3 shape mismatch is refused without touching the sums
   !>
   !> @param[out] error Test failure
   subroutine test_response_density_ranks(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Response under test
      type(response_type) :: response
      !> Density items, one per accumulation, so each carries only its own arrays
      type(density_response_type) :: rho_only, full, hess_only, grad_only, bad_hess
      !> Hessian weights of the full item
      real(wp) :: hess(3, 3, 2)
      !> Element index
      integer :: i

      allocate (rho_only%w_rho, source=[1.0_wp, 2.0_wp])
      call response_accumulate(response, rho_only, err)
      call check_moist_error(error, err, "density value weights")
      if (allocated(error)) return

      hess = reshape([(real(i, wp), i = 1, 18)], [3, 3, 2])
      full%w_rho = [10.0_wp, 20.0_wp]
      full%w_grad_rho = reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp], [3, 2])
      full%w_hess_rho = hess
      call response_accumulate(response, full, err)
      call check_moist_error(error, err, "full density item")
      if (allocated(error)) return

      allocate (hess_only%w_hess_rho(3, 3, 2), source=1.0_wp)
      call response_accumulate(response, hess_only, err)
      call check_moist_error(error, err, "Hessian weights only")
      if (allocated(error)) return

      allocate (grad_only%w_grad_rho(3, 2), source=100.0_wp)
      call response_accumulate(response, grad_only, err)
      call check_moist_error(error, err, "gradient weights only")
      if (allocated(error)) return

      allocate (bad_hess%w_hess_rho(3, 3, 1), source=0.0_wp)
      call response_accumulate(response, bad_hess, err)
      call check(error, allocated(err), more="a differently shaped Hessian was accumulated")
      if (allocated(error)) return
      call check(error, index(err%message, "'density'") > 0 &
         & .and. index(err%message, "'w_hess_rho'") > 0 &
         & .and. index(err%message, "different shape than the stored one") > 0, &
         & "error names item and array")
      if (allocated(error)) return

      call check(error, walk_names(response), "density,", more="density items were not merged")
      if (allocated(error)) return
      call check(error, response%next())
      if (allocated(error)) return
      select type (item => response%item())
      type is (density_response_type)
         call check(error, all(item%w_rho == [11.0_wp, 22.0_wp]), "value weights sum")
         if (allocated(error)) return
         call check(error, all(item%w_grad_rho == full%w_grad_rho + 100.0_wp), &
            & "gradient weights copied, then summed")
         if (allocated(error)) return
         call check(error, all(shape(item%w_hess_rho) == [3, 3, 2]), "stored Hessian shape unchanged")
         if (allocated(error)) return
         call check(error, all(item%w_hess_rho == hess + 1.0_wp), "Hessian weights copied, then summed")
      class default
         call test_failed(error, "the stored item is not the density")
      end select
   end subroutine test_response_density_ranks

   !> A differently sized rank-1 array is refused for the potential adjoint and
   !> for the first GOSTSHYP amplitude, and the stored item is left unchanged
   !>
   !> @param[out] error Test failure
   subroutine test_response_vector_shape(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Response under test
      type(response_type) :: response
      !> Potential adjoint item
      type(potential_adjoint_response_type) :: charge
      !> GOSTSHYP amplitude item
      type(gostshyp_amplitude_response_type) :: amplitude

      allocate (charge%w_phi, source=[1.0_wp, 2.0_wp, 3.0_wp])
      call response_accumulate(response, charge, err)
      call check_moist_error(error, err, "first potential adjoint")
      if (allocated(error)) return
      charge%w_phi = [1.0_wp, 2.0_wp]
      call response_accumulate(response, charge, err)
      call check(error, allocated(err), more="a shorter w_phi was accumulated")
      if (allocated(error)) return
      call check(error, index(err%message, "'potential_adjoint'") > 0 &
         & .and. index(err%message, "'w_phi'") > 0 &
         & .and. index(err%message, "different shape than the stored one") > 0, &
         & "error names item and array")
      if (allocated(error)) return

      amplitude%w_overlap = [1.0_wp, 2.0_wp]
      amplitude%w_normal_deriv = [3.0_wp, 4.0_wp]
      call response_accumulate(response, amplitude, err)
      call check_moist_error(error, err, "first amplitude")
      if (allocated(error)) return
      amplitude%w_overlap = [1.0_wp]
      call response_accumulate(response, amplitude, err)
      call check(error, allocated(err), more="a shorter w_overlap was accumulated")
      if (allocated(error)) return
      call check(error, index(err%message, "'gostshyp_amplitude'") > 0 &
         & .and. index(err%message, "'w_overlap'") > 0, "error names item and array")
      if (allocated(error)) return

      do while (response%next())
         select type (item => response%item())
         type is (potential_adjoint_response_type)
            call check(error, all(item%w_phi == [1.0_wp, 2.0_wp, 3.0_wp]), "stored w_phi unchanged")
         type is (gostshyp_amplitude_response_type)
            ! The refusal stops before the normal derivative is summed
            call check(error, all(item%w_overlap == [1.0_wp, 2.0_wp]) &
               & .and. all(item%w_normal_deriv == [3.0_wp, 4.0_wp]), "stored amplitudes unchanged")
         class default
            call test_failed(error, "unexpected response item "//item%name())
         end select
         if (allocated(error)) return
      end do
   end subroutine test_response_vector_shape

   !> The "no_current_item" placeholder stands for no item, so nothing
   !> accumulates into it
   !>
   !> @param[out] error Test failure
   subroutine test_response_placeholder_add(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Response under test
      type(response_type) :: response
      !> Placeholder returned outside a pass
      class(response_item_type), allocatable :: placeholder

      placeholder = response%item()
      call check(error, trim(placeholder%name()), "no_current_item")
      if (allocated(error)) return
      ! The first accumulation stores a copy; the second has to add into it
      call response_accumulate(response, placeholder, err)
      call response_accumulate(response, placeholder, err)
      call check(error, allocated(err), more="the placeholder accumulated an item")
      if (allocated(error)) return
      call check(error, index(err%message, "placeholder response item 'no_current_item'") > 0 &
         & .and. index(err%message, "cannot accumulate 'no_current_item'") > 0)
   end subroutine test_response_placeholder_add

   !* ================================================================================= *!
   !*                                   Model registry                                  *!
   !* ================================================================================= *!

   !> Releasing middle, head and tail nodes leaves the remaining collections usable
   !>
   !> @param[out] error Test failure
   subroutine test_registry_release(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Registry under test
      type(coupling_registry_type), target :: registry
      !> Borrowed collections
      type(coupling_type), pointer :: first, middle, last
      !> Repeated allocation index
      integer :: i
      call registry%mint(first, err)
      call registry%mint(middle, err)
      call registry%mint(last, err)
      call check_moist_error(error, err, "mint")
      if (allocated(error)) return
      call registry%release(middle)
      call check(error, .not. associated(middle))
      if (allocated(error)) return
      call check(error, registry%owns(first) .and. registry%owns(last))
      if (allocated(error)) return
      call registry%release(last)
      call check(error, .not. associated(last))
      if (allocated(error)) return
      call registry%release(first)
      call check(error, .not. associated(first))
      if (allocated(error)) return
      do i = 1, 100
         call registry%mint(first, err)
         call fixture(first, err)
         call check_moist_error(error, err, "reuse after release")
         if (allocated(error)) return
         call registry%release(first)
         call check(error, .not. associated(first))
         if (allocated(error)) return
      end do
      call registry%clear()
   end subroutine test_registry_release

   !> Releasing the oldest of three collections walks past the newer ones; a
   !> collection the registry never minted is left alone
   !>
   !> @param[out] error Test failure
   subroutine test_registry_release_tail(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Registry under test
      type(coupling_registry_type), target :: registry
      !> Borrowed collections, oldest first
      type(coupling_type), pointer :: oldest, middle, newest
      !> Collection minted elsewhere
      type(coupling_type), target :: foreign
      !> Pointer to the foreign collection
      type(coupling_type), pointer :: stray
      call registry%mint(oldest, err)
      call check_moist_error(error, err, "mint oldest")
      if (allocated(error)) return
      call registry%mint(middle, err)
      call check_moist_error(error, err, "mint middle")
      if (allocated(error)) return
      call registry%mint(newest, err)
      call check_moist_error(error, err, "mint newest")
      if (allocated(error)) return
      ! The newest node is the head, so the oldest is two links down
      call registry%release(oldest)
      call check(error, .not. associated(oldest), more="the released tail kept its pointer")
      if (allocated(error)) return
      call check(error, registry%owns(middle) .and. registry%owns(newest), &
         & more="releasing the tail dropped another collection")
      if (allocated(error)) return
      stray => foreign
      call registry%release(stray)
      call check(error, associated(stray, foreign), more="releasing a foreign collection nullified it")
      if (allocated(error)) return
      call check(error, .not. registry%owns(foreign), more="the registry claims a foreign collection")
      if (allocated(error)) return
      call check(error, registry%owns(middle) .and. registry%owns(newest), &
         & more="releasing a foreign collection dropped an owned one")
      if (allocated(error)) return
      call fixture(middle, err)
      call check_moist_error(error, err, "survivor still usable")
      if (allocated(error)) return
      call registry%clear()
   end subroutine test_registry_release_tail

   !> Requests with matching inputs are shared across scopes, and a request no
   !> scope declares any more is dropped at the next grid snapshot
   !>
   !> @param[out] error Test failure
   subroutine test_shared_requests(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Collection and the calculations registered on it
      type(coupling_type) :: coupling
      type(point_potential_request_type) :: point
      type(gaussian_moment_request_type) :: moments
      call coupling_begin_registration(coupling)
      call request_require(point, moist_phase_energy, "phi", err)
      call coupling_register(coupling, "potential", point, err)
      call coupling_set_scope(coupling, 2)
      call coupling_register(coupling, "potential", point, err)
      call check_moist_error(error, err, "shared registration")
      if (allocated(error)) return
      call coupling_snapshot(coupling, 1)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check(error, count_visits(coupling), 1, more="matching inputs were not shared")
      if (allocated(error)) return
      moments%width = [1.0_wp]
      call request_require(moments, moist_phase_energy, "gt", err)
      call coupling_begin_registration(coupling)
      call coupling_register(coupling, "potential", point, err)
      call coupling_set_scope(coupling, 3)
      call coupling_register(coupling, "moments", moments, err)
      call coupling_snapshot(coupling, 1)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "moments added")
      if (allocated(error)) return
      call check(error, count_visits(coupling), 2)
      if (allocated(error)) return
      ! The moments are no longer declared; the next snapshot drops them
      call coupling_begin_registration(coupling)
      call coupling_set_scope(coupling, 1)
      call coupling_register(coupling, "potential", point, err)
      call coupling_snapshot(coupling, 1)
      call coupling_arm(coupling, moist_phase_energy, err)
      call check_moist_error(error, err, "moments dropped")
      if (allocated(error)) return
      call check(error, count_visits(coupling), 1, more="an obsolete request was kept")
      if (allocated(error)) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "point_potential")
   end subroutine test_shared_requests

   !> Direct potential reads return independent outputs and respect view lifetime
   !>
   !> @param[out] error Test failure
   subroutine test_direct_potential_read(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Collection and scoped view
      type(coupling_type), target :: coupling
      type(coupling_view_type) :: view
      !> Requested output copy
      real(wp), allocatable :: phi(:)
      call fixture(coupling, err)
      call check_moist_error(error, err, "fixture")
      if (allocated(error)) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call coupling%answer("phi", [2.0_wp, 3.0_wp], err)
      call check_moist_error(error, err, "phi answer")
      if (allocated(error)) return
      call coupling_make_view(coupling, 1, view)
      call view%read("potential", "phi", phi, err)
      call check_moist_error(error, err, "direct potential read")
      if (allocated(error)) return
      phi = -9.0_wp
      call view%read("potential", "phi", phi, err)
      call check_moist_error(error, err, "independent potential output")
      if (allocated(error)) return
      call check(error, maxval(abs(phi - [2.0_wp, 3.0_wp])), 0.0_wp, thr=epsilon(1.0_wp))
      if (allocated(error)) return
      call coupling_close_view(coupling)
      call view%read("potential", "phi", phi, err)
      call check(error, allocated(err))
   end subroutine test_direct_potential_read

   !> A quiet NaN, built without arithmetic that would trap
   function ieee_nan() result(nan)
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      real(wp) :: nan
      nan = ieee_value(1.0_wp, ieee_quiet_nan)
   end function ieee_nan

end module test_channels
