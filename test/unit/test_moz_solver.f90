!> Numerical and lifecycle checks for the shared MOZ solvers
module test_moz_solver
   !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
   use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_associated, c_null_ptr
   use mctc_env, only: wp, mctc_error => error_type, fatal_error
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_math_blas, only: gemv_rows, gemm_rows_tn, gemm_rows_nn
   use moist_math_orthogonalization, only: moist_orthogonalize_mgs
   use moist_model_moz_solver, only: moist_model_moz_solver_type, &
      & moist_model_moz_solver_statistics_type, moist_model_moz_solver_evaluator_type, &
      & moist_model_moz_solver_preconditioner_type, moist_model_moz_solver_printer_type, &
      & moist_model_moz_solver_parameters_type, moist_model_moz_solver_iteration_type, &
      & moist_model_moz_solver_deflation_type, new_solver, run_iterations, solve_adjoint, &
      & adjoint_parameters => moist_model_moz_solver_jfnk_adjoint_parameters_type
   use moist_model_moz_solver_jfnk, only: moist_model_moz_solver_jfnk_type, new_jfnk_solver, &
      & jfnk_parameters => moist_model_moz_solver_jfnk_parameters_type
   use moist_model_moz_solver_picard, only: moist_model_moz_solver_picard_type, new_picard_solver, &
      & picard_parameters => moist_model_moz_solver_picard_parameters_type
   use moist_model_moz_solver_diis, only: moist_model_moz_solver_diis_type, new_diis_solver, &
      & diis_parameters => moist_model_moz_solver_diis_parameters_type
   use moist_model_moz_solver_core_preconditioner, only: check_preconditioner_name
   use moist_model_moz_solver_jfnk_workspace, only: moist_model_moz_solver_workspace_type
   implicit none(type, external)
   private
   public :: collect_moz_solver

   type, extends(moist_model_moz_solver_evaluator_type) :: diagonal_evaluator_type
   contains
      procedure :: residual => diagonal_residual
   end type diagonal_evaluator_type

   type, extends(moist_model_moz_solver_evaluator_type) :: identity_evaluator_type
   contains
      procedure :: residual => identity_residual
   end type identity_evaluator_type

   !> Small nonsymmetric problem with controlled callback failures
   type, extends(moist_model_moz_solver_evaluator_type) :: linear_evaluator_type
      !> Nonsymmetric dense Jacobian
      real(wp) :: matrix(2, 2) = reshape([2.0_wp, 1.0_wp, -0.5_wp, 3.0_wp], [2, 2])
      !> Exact root
      real(wp) :: target(2) = [1.0_wp, -2.0_wp]
      !> Largest allowed trial component
      real(wp) :: trial_limit = huge(1.0_wp)
      !> Residual: 1 = error, 2 = NaN; JVP: 1 = error, 2 = NaN, 3 = zero
      integer :: residual_mode = 0, jvp_mode = 0
      !> Observe attempts to change the frozen primal linearization.
      integer, pointer :: linearizations => null()
   contains
      procedure :: residual => linear_residual
      procedure :: has_jvp => linear_has_jvp
      procedure :: linearize => linear_linearize
      procedure :: jvp => linear_jvp
      procedure :: has_adjoint => linear_has_jvp
      procedure :: adjoint => linear_adjoint
      procedure :: deflation_setup => linear_deflation_setup
   end type linear_evaluator_type

   type, extends(moist_model_moz_solver_evaluator_type) :: infeasible_evaluator_type
   contains
      procedure :: residual => infeasible_residual
   end type infeasible_evaluator_type

   type :: soft_linearization_type
      real(wp), allocatable :: point(:)
   end type soft_linearization_type

   type, extends(moist_model_moz_solver_evaluator_type) :: soft_evaluator_type
      type(soft_linearization_type), pointer :: state => null()
   contains
      procedure :: residual => soft_residual
      procedure :: has_jvp => soft_has_jvp
      procedure :: linearize => soft_linearize
      procedure :: jvp => soft_jvp
   end type soft_evaluator_type

   type, extends(moist_model_moz_solver_evaluator_type) :: fixed_evaluator_type
      real(wp) :: value(1, 2) = 0.0_wp
   contains
      procedure :: residual => fixed_residual
   end type fixed_evaluator_type

   type, extends(moist_model_moz_solver_evaluator_type) :: coupled_evaluator_type
   contains
      procedure :: residual => coupled_residual
   end type coupled_evaluator_type

   type :: quadratic_linearization_type
      real(wp) :: point = 0.0_wp
   end type quadratic_linearization_type

   type, extends(moist_model_moz_solver_evaluator_type) :: quadratic_evaluator_type
      type(quadratic_linearization_type), pointer :: state => null()
      real(wp) :: root = 1.0_wp
      real(wp) :: max_argument = huge(1.0_wp)
      logical :: accept_linearization = .true.
   contains
      procedure :: residual => quadratic_residual
      procedure :: has_jvp => quadratic_has_jvp
      procedure :: linearize => quadratic_linearize
      procedure :: jvp => quadratic_jvp
   end type quadratic_evaluator_type

   type, extends(moist_model_moz_solver_preconditioner_type) :: diagonal_preconditioner_type
      integer :: napplied = 0
   contains
      procedure :: ready => diagonal_ready
      procedure :: apply => diagonal_apply
      procedure :: describe => diagonal_describe
   end type diagonal_preconditioner_type

   !> Non-diagonal left operator, with scaling and unusable-output cases
   type, extends(moist_model_moz_solver_preconditioner_type) :: matrix_preconditioner_type
      !> Fixed inverse operator; it does not commute with the test Jacobian
      real(wp) :: matrix(2, 2) = reshape([1.0_wp, -0.5_wp, 0.25_wp, 2.0_wp], [2, 2])
      !> Arbitrary scalar multiplier of the operator
      real(wp) :: scale = 1.0_wp
      !> Whether ready allows applications
      logical :: enabled = .true.
      !> Application count and bad-output mode: 1 = NaN, 2 = zero
      integer :: napplied = 0, output_mode = 0
      !> Optional observation of the forcing term actually used by GMRES
      real(wp), pointer :: forcing => null()
      !> Last observed forcing term
      real(wp) :: applied_eta = -1.0_wp
   contains
      procedure :: ready => matrix_ready
      procedure :: apply => matrix_apply
      procedure :: describe => matrix_describe
   end type matrix_preconditioner_type

   type, extends(moist_model_moz_solver_iteration_type) :: iteration_context_type
      integer :: nsteps = 0, nobserved = 0, ncheckpoints = 0, fail_at = 0
      logical :: fail_checkpoint = .false.
      !> Smallest observed advance time, s; huge before any observation
      real(wp) :: min_elapsed = huge(1.0_wp)
   contains
      procedure :: advance => context_advance
      procedure :: observe => context_observe
      procedure :: checkpoint => context_checkpoint
   end type iteration_context_type

   real(wp), parameter :: soft_cubic = 0.2_wp

contains

   subroutine collect_moz_solver(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("factory", test_factory), &
         & new_unittest("parameters_file", test_parameters_file), &
         & new_unittest("picard_update", test_picard_update), &
         & new_unittest("diis_extrapolation", test_diis_extrapolation), &
         & new_unittest("diis_convergence", test_diis_convergence), &
         & new_unittest("diis_stages", test_diis_stages), &
         & new_unittest("diis_divergence_restart", test_diis_divergence_restart), &
         & new_unittest("adjoint_solve", test_adjoint_solve), &
         & new_unittest("adjoint_failures", test_adjoint_failures), &
         & new_unittest("adjoint_shared_workspace", test_adjoint_shared_workspace), &
         & new_unittest("adjoint_workspace_recovery", test_adjoint_workspace_recovery), &
         & new_unittest("jfnk_nonlinear", test_jfnk_nonlinear), &
         & new_unittest("jfnk_preconditioner", test_jfnk_preconditioner), &
         & new_unittest("jfnk_line_search", test_jfnk_line_search), &
         & new_unittest("jfnk_forcing", test_jfnk_forcing), &
         & new_unittest("jfnk_deflation", test_jfnk_deflation), &
         & new_unittest("jfnk_deflation_warnings", test_jfnk_deflation_warnings), &
         & new_unittest("jfnk_recycle_warnings", test_jfnk_recycle_warnings), &
         & new_unittest("jfnk_stage_summary", test_jfnk_stage_summary), &
         & new_unittest("jfnk_exact_newton", test_jfnk_exact_newton), &
         & new_unittest("jfnk_exact_required", test_jfnk_exact_required), &
         & new_unittest("jfnk_validate", test_jfnk_validate), &
         & new_unittest("hybrid_switch", test_hybrid_switch), &
         & new_unittest("hybrid_convergence", test_hybrid_convergence), &
         & new_unittest("hybrid_lifecycle", test_hybrid_lifecycle), &
         & new_unittest("solver_failures", test_solver_failures), &
         & new_unittest("jfnk_callback_errors", test_jfnk_callback_errors), &
         & new_unittest("preconditioner_contract", test_preconditioner_contract), &
         & new_unittest("workspace_reconfigure", test_workspace_reconfigure), &
         & new_unittest("diis_residual_precision", test_diis_residual_precision), &
         & new_unittest("iteration_controller", test_iteration_controller), &
         & new_unittest("iteration_timer", test_iteration_timer), &
         & new_unittest("deflation_algebra", test_deflation_algebra), &
         & new_unittest("concurrent_solves", test_concurrent_solves), &
         & new_unittest("mgs_streaming", test_mgs_streaming), &
         & new_unittest("blas_row_blocks", test_blas_row_blocks), &
         & new_unittest("solver_lifecycle", test_solver_lifecycle), &
         & new_unittest("reduction_thread_independence", test_reduction_threads), &
         & new_unittest("workspace_and_parallel_step", test_parallel_step), &
         & new_unittest("parallel_restart", test_parallel_restart), &
         & new_unittest("parallel_cgs2", test_parallel_cgs2), &
         & new_unittest("happy_breakdown", test_happy_breakdown), &
         & new_unittest("infeasible_jvp", test_infeasible_jvp), &
         & new_unittest("recycle_small", test_recycle_small), &
         & new_unittest("recycle_reconfigure", test_recycle_reconfigure), &
         & new_unittest("recycle_parallel", test_recycle_parallel)]
   end subroutine collect_moz_solver

   !> Only the requested solver families are constructed, with validated controls
   subroutine test_factory(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      class(moist_model_moz_solver_type), allocatable :: solver
      type(moist_model_moz_solver_parameters_type) :: p
      type(moist_context_type), target :: ctx, pinned
      character(len=8), parameter :: names(5) = [character(len=8) :: "Picard", "DIIS", " mdiis ", "jFnK", " hYbRiD "]
      character(len=8), parameter :: rejected(3) = [character(len=8) :: "LBFGS", "GMRES", "unknown"]
      integer :: i

      call new_context(ctx, nthreads=0, verbosity=0)
      do i = 1, size(names)
         p%name = names(i)
         call new_solver(solver, ctx, merr, p)
         call check(error, .not. allocated(merr) .and. allocated(solver))
         if (allocated(error)) return
         select case (i)
         case (1)
            call check(error, solver%name() == "Picard" .and. solver%allows_postprocessing())
         case (2, 3)
            call check(error, solver%name() == "DIIS" .and. solver%allows_postprocessing())
         case (4)
            call check(error, solver%name() == "JFNK" .and. .not. solver%allows_postprocessing())
         case (5)
            call check(error, solver%name() == "HYBRID" .and. solver%allows_postprocessing())
         case default
            call test_failed(error, "Unexpected solver test case")
         end select
         if (allocated(error)) return
      end do
      do i = 1, size(rejected)
         p%name = rejected(i)
         call new_solver(solver, ctx, merr, p)
         call check(error, allocated(merr) .and. .not. allocated(solver))
         if (allocated(error)) return
      end do

      ! The factory hands the context thread count to the solver
      call new_context(pinned, nthreads=3, verbosity=0)
      p%name = "HYBRID"
      call new_solver(solver, pinned, merr, p)
      call check(error, .not. allocated(merr) .and. solver%team_size() == 3, "Context thread count")
      if (allocated(error)) return

      p%picard%damping = ieee_value(0.0_wp, ieee_quiet_nan)
      p%name = "PICARD"
      call new_solver(solver, ctx, merr, p)
      call check(error, allocated(merr) .and. .not. allocated(solver))
      if (allocated(error)) return
      call p%init_defaults()
      p%diis%depth = 0
      p%name = "DIIS"
      call new_solver(solver, ctx, merr, p)
      call check(error, allocated(merr) .and. .not. allocated(solver))
      if (allocated(error)) return
      p%diis%depth = 3
      p%diis%step_limit_tau = 1.0_wp
      p%name = "DIIS"
      call new_solver(solver, ctx, merr, p)
      call check(error, allocated(merr) .and. .not. allocated(solver))
      if (allocated(error)) return
      call p%init_defaults()
      p%jfnk%restart = 0
      p%name = "JFNK"
      call new_solver(solver, ctx, merr, p)
      call check(error, allocated(merr) .and. .not. allocated(solver))
      if (allocated(error)) return
      p%jfnk%restart = 3
      p%jfnk%eta_min = 0.2_wp
      p%name = "JFNK"
      call new_solver(solver, ctx, merr, p)
      call check(error, allocated(merr) .and. .not. allocated(solver))
      if (allocated(error)) return
      call p%init_defaults()
      p%picard%damping = 0.4_wp
      p%diis%depth = 3
      p%jfnk%restart = 4
      p%jfnk%max_krylov = 8
      p%jfnk%use_exact_jvp = .false.
      p%jfnk%recycle = 2
      p%name = "JFNK"
      call new_solver(solver, ctx, merr, p)
      call check(error, .not. allocated(merr))
      if (allocated(error)) return
      select type (solver)
      type is (moist_model_moz_solver_jfnk_type)
         call check(error, solver%params%restart == 4 .and. solver%params%max_krylov == 8 .and. &
            & .not. solver%params%use_exact_jvp .and. solver%params%recycle == 2)
      class default
         call test_failed(error, "JFNK factory returned the wrong solver type")
      end select
      if (allocated(error)) return

      ! Invalid line-search and forcing controls are rejected
      call p%init_defaults()
      p%jfnk%line_search_tau = 1.0_wp
      p%name = "JFNK"
      call new_solver(solver, ctx, merr, p)
      call check(error, allocated(merr) .and. .not. allocated(solver), "Backtracking factor must be below one")
      if (allocated(error)) return
      call p%init_defaults()
      p%jfnk%forcing_power = 3.0_wp
      p%name = "JFNK"
      call new_solver(solver, ctx, merr, p)
      call check(error, allocated(merr) .and. .not. allocated(solver), "Forcing exponent must lie in (1, 2]")
      if (allocated(error)) return

      ! Hybrid construction validates both component solvers and its threshold.
      do i = 1, 5
         call p%init_defaults()
         select case (i)
         case (1)
            p%hybrid%switch_rmse = 0.0_wp
         case (2)
            p%hybrid%switch_rmse = -1.0_wp
         case (3)
            p%hybrid%switch_rmse = ieee_value(0.0_wp, ieee_quiet_nan)
         case (4)
            p%diis%depth = 0
         case (5)
            p%jfnk%restart = 0
         case default
            call test_failed(error, "Unexpected hybrid parameter test case")
            return
         end select
         p%name = "HYBRID"
         call new_solver(solver, ctx, merr, p)
         call check(error, allocated(merr) .and. .not. allocated(solver))
         if (allocated(error)) return
      end do
      call pinned%delete()
      call ctx%delete()
   end subroutine test_factory

   !> Read nested JSON keys, round-trip JSON and TOML, and validate on read
   !>
   !> @param[out] error Test failure information
   subroutine test_parameters_file(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Parameters read from files and their defaults
      type(moist_model_moz_solver_parameters_type) :: p, defaults
      !> Library error
      type(mctc_error), allocatable :: merr
      !> File unit and I/O statuses
      integer :: unit, stat, close_stat, i
      character(len=*), parameter :: path = "moist-test-moz-solver-parameters.json"
      character(len=*), parameter :: stem = "moist-test-moz-solver-parameters-roundtrip"
      character(len=5), parameter :: extensions(2) = [".json", ".toml"]

      ! Omitted keys keep their defaults
      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)", iostat=stat) '{"diis": {"depth": 5}, '// &
         & '"jfnk": {"restart": 7, "eta0": 0.25, "use_exact_jvp": false, "recycle": 3, "line_search_max_iter": 4, '// &
         & '"adjoint": {"restart": 4, "max_krylov": 13}}}'
      close (unit, iostat=close_stat)
      call check(error, stat == 0 .and. close_stat == 0)
      if (allocated(error)) return
      call p%read_file(path, merr)
      call check(error, .not. allocated(merr))
      if (allocated(error)) return
      call check(error, p%diis%depth == 5 .and. p%jfnk%restart == 7 .and. &
         & p%jfnk%eta0 == 0.25_wp .and. .not. p%jfnk%use_exact_jvp .and. p%jfnk%line_search_max_iter == 4 .and. &
         & p%jfnk%recycle == 3 .and. &
         & p%jfnk%adjoint%restart == 4 .and. p%jfnk%adjoint%max_krylov == 13 .and. &
         & p%jfnk%max_krylov == defaults%jfnk%max_krylov .and. p%picard%damping == defaults%picard%damping, &
         & "Nested keys must set their fields and omitted keys keep defaults")
      if (allocated(error)) return

      ! Reading validates: eta_min above eta_max is rejected
      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)", iostat=stat) '{"jfnk": {"eta_min": 0.5, "eta_max": 0.1}}'
      close (unit, iostat=close_stat)
      call check(error, stat == 0 .and. close_stat == 0)
      if (allocated(error)) return
      call p%read_file(path, merr)
      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      call check(error, allocated(merr), "Invalid forcing bounds must be rejected on read")
      if (allocated(error)) return

      do i = 1, size(extensions)
         call p%init_defaults()
         p%diis%step_limit = .false.
         p%jfnk%deflation_rmse = 2.5e-4_wp
         p%jfnk%adjoint%restart = 6
         p%jfnk%adjoint%max_krylov = 17
         call p%write_file(stem//extensions(i), merr)
         call check(error, .not. allocated(merr))
         if (allocated(error)) return
         call p%init_defaults()
         call p%read_file(stem//extensions(i), merr)
         open (newunit=unit, file=stem//extensions(i), status="old", iostat=stat)
         if (stat == 0) close (unit, status="delete", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call check(error, .not. allocated(merr) .and. .not. p%diis%step_limit .and. &
            & p%jfnk%deflation_rmse == 2.5e-4_wp .and. &
            & p%jfnk%adjoint%restart == 6 .and. p%jfnk%adjoint%max_krylov == 17, "Round trip "//extensions(i))
         if (allocated(error)) return
      end do

      open (newunit=unit, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      call p%print_parameters(merr, unit)
      close (unit)
      call check(error, .not. allocated(merr), "Printing the parameters")
   end subroutine test_parameters_file

   !> Use input RMS for an immediate one-way switch and restart each stage in MDIIS
   !>
   !> @param[out] error Test failure information
   subroutine test_hybrid_switch(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      class(moist_model_moz_solver_type), allocatable :: solver
      type(moist_model_moz_solver_parameters_type) :: p
      type(moist_model_moz_solver_statistics_type) :: stats
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(1, 2), updated(1, 2), rmse
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      p%picard%damping = 0.5_wp
      p%jfnk%fd_eps = 1.0_wp
      p%jfnk%use_exact_jvp = .false.
      p%name = "HYBRID"
      call new_solver(solver, ctx, merr, p)
      call check(error, .not. allocated(merr))
      if (allocated(error)) return
      gamma = 2.0e-3_wp
      call solver%outer_step(identity_evaluator_type(), gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. solver%allows_postprocessing())
      if (allocated(error)) return
      stats = solver%get_statistics()
      call check(error, all(updated == 1.0e-3_wp) .and. stats%ndiis == 1 .and. stats%nresidual == 1)
      if (allocated(error)) return

      ! The maximum component exceeds the threshold, but the RMS is below it.
      gamma(1, :) = [1.2e-3_wp, 0.0_wp]
      call solver%outer_step(identity_evaluator_type(), gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. .not. solver%allows_postprocessing())
      if (allocated(error)) return
      stats = solver%get_statistics()
      call check(error, maxval(abs(updated)) < 1.0e-12_wp .and. stats%ndiis == 0 .and. stats%ngmres > 0)
      if (allocated(error)) return
      call check(error, abs(solver%prev_step_maxabs - maxval(abs(updated - gamma))) < 1.0e-14_wp)
      if (allocated(error)) return

      ! A later large residual stays in JFNK.
      gamma = 0.01_wp
      call solver%outer_step(identity_evaluator_type(), gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. .not. solver%allows_postprocessing())
      if (allocated(error)) return
      call check(error, maxval(abs(updated)) < 1.0e-12_wp)
      if (allocated(error)) return
      call solver%begin_stage()
      stats = solver%get_statistics()
      call check(error, solver%allows_postprocessing() .and. stats%ndiis == 0 .and. stats%ngmres == 0)
      if (allocated(error)) return
      gamma = 2.0e-3_wp
      call solver%outer_step(identity_evaluator_type(), gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. all(updated == 1.0e-3_wp))
      if (allocated(error)) return
      call solver%begin_stage()
      stats = solver%get_statistics()
      call check(error, stats%ndiis == 0)
      if (allocated(error)) return

      ! An exactly representable boundary distinguishes <= from <.
      p%hybrid%switch_rmse = 0.5_wp
      p%name = "HYBRID"
      call new_solver(solver, ctx, merr, p)
      gamma = 0.5_wp
      call solver%outer_step(identity_evaluator_type(), gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. rmse == 0.5_wp &
         & .and. .not. solver%allows_postprocessing() .and. maxval(abs(updated)) < 1.0e-12_wp)
      call ctx%delete()
   end subroutine test_hybrid_switch

   !> Converge a nonlinear problem through both hybrid phases
   !>
   !> @param[out] error Test failure information
   subroutine test_hybrid_convergence(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Evaluator failure information
      type(mctc_error), allocatable :: merr
      !> Solver used through the common interface
      class(moist_model_moz_solver_type), allocatable :: solver
      !> Default handoff threshold, short startup history
      type(moist_model_moz_solver_parameters_type) :: p
      !> Nonlinear residual and its borrowed linearization storage
      type(quadratic_evaluator_type) :: evaluator
      type(quadratic_linearization_type), target :: state
      !> Quiet run context borrowed by the solver
      type(moist_context_type), target :: ctx
      !> Input, update and input residual norm
      real(wp) :: gamma(1, 1), updated(1, 1), rmse
      !> Attempt count and observed phases
      integer :: it
      logical :: feasible, saw_mdiis, saw_jfnk

      call new_context(ctx, nthreads=0, verbosity=0)
      evaluator%state => state
      p%diis%depth = 3
      p%name = "HYBRID"
      call new_solver(solver, ctx, merr, p)
      call check(error, .not. allocated(merr))
      if (allocated(error)) return
      gamma = 0.1_wp
      saw_mdiis = .false.
      saw_jfnk = .false.
      do it = 1, 128
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         saw_mdiis = saw_mdiis .or. solver%allows_postprocessing()
         saw_jfnk = saw_jfnk .or. .not. solver%allows_postprocessing()
         if (rmse < 1.0e-10_wp) exit
         gamma = updated
      end do
      call check(error, saw_mdiis .and. saw_jfnk .and. it <= 128 .and. abs(gamma(1, 1) - 1.0_wp) < 1.0e-9_wp, &
         & "Hybrid must converge through MDIIS startup and Newton refinement")
      call ctx%delete()
   end subroutine test_hybrid_convergence

   !> Forward borrowed operators and diagnostics, and retain exact products at handoff
   !>
   !> @param[out] error Test failure information
   subroutine test_hybrid_lifecycle(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      class(moist_model_moz_solver_type), allocatable :: solver
      type(moist_model_moz_solver_parameters_type) :: p
      type(linear_evaluator_type) :: evaluator
      type(matrix_preconditioner_type), target :: precond, replacement
      type(moist_model_moz_solver_printer_type), target :: printer
      type(moist_model_moz_solver_statistics_type) :: stats
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(2, 1), updated(2, 1), rmse
      integer :: unit, calls
      logical :: feasible

      p%hybrid%switch_rmse = 10.0_wp
      p%jfnk%eta0 = 1.0e-12_wp
      p%jfnk%restart = 2
      p%jfnk%max_krylov = 2
      ! Verbosity 4 runs the switch note; the scratch unit keeps it off the console
      open (newunit=unit, status="scratch", action="readwrite")
      call new_context(ctx, nthreads=0, verbosity=4, unit=unit)
      p%name = "HYBRID"
      call new_solver(solver, ctx, merr, p)
      call check(error, .not. allocated(merr) .and. solver%uses_preconditioner())
      if (allocated(error)) then
         close (unit)
         return
      end if
      call printer%setup(unit=unit)
      call solver%set_iteration_printer(printer)
      call check(error, printer%diis_column .and. .not. printer%gmres_block, "An attached hybrid shows the DIIS columns")
      if (allocated(error)) then
         close (unit)
         return
      end if
      call solver%set_preconditioner(precond)
      gamma = 0.0_wp
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. precond%napplied > 0)
      if (allocated(error)) then
         close (unit)
         return
      end if
      stats = solver%get_statistics()
      call check(error, stats%njvp > 0 .and. maxval(abs(updated(:, 1) - evaluator%target)) < 1.0e-11_wp)
      call solver%set_iteration_printer()
      close (unit)
      if (allocated(error)) return
      call check(error, printer%gmres_block .and. .not. printer%diis_column, "The hybrid switches to the JFNK columns")
      if (allocated(error)) return

      calls = precond%napplied
      call solver%set_preconditioner(replacement)
      call solver%begin_stage()
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. replacement%napplied > 0 &
         & .and. precond%napplied == calls)
      if (allocated(error)) return
      calls = replacement%napplied
      call solver%set_preconditioner()
      call solver%begin_stage()
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. replacement%napplied == calls)
      if (allocated(error)) return

      ! Rejected MDIIS inputs do not advance the phase or mixing history.
      call solver%begin_stage()
      evaluator%residual_mode = 1
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      stats = solver%get_statistics()
      call check(error, allocated(merr) .and. solver%allows_postprocessing() .and. stats%ndiis == 0)
      call ctx%delete()
   end subroutine test_hybrid_lifecycle

   !> Shape errors and unusable residuals fail before an update can be accepted
   !>
   !> @param[out] error Test failure information
   subroutine test_solver_failures(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      class(moist_model_moz_solver_type), allocatable :: solver
      type(moist_model_moz_solver_parameters_type) :: p
      type(linear_evaluator_type) :: evaluator
      type(moist_context_type), target :: ctx
      character(len=6), parameter :: names(*) = ["PICARD", "DIIS  ", "JFNK  ", "HYBRID"]
      real(wp) :: gamma(2, 1), updated(2, 1), wrong(1, 1), empty(0, 1), empty_out(0, 1), rmse
      integer :: i, mode
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      gamma = 0.0_wp
      do i = 1, size(names)
         p%name = names(i)
         call new_solver(solver, ctx, merr, p)
         call solver%outer_step(evaluator, gamma, wrong, rmse, feasible, merr)
         call check(error, allocated(merr) .and. .not. feasible, "Reject mismatched solver fields")
         if (allocated(error)) return
         call solver%outer_step(identity_evaluator_type(), empty, empty_out, rmse, feasible, merr)
         call check(error, allocated(merr) .and. .not. feasible, "Reject empty solver fields")
         if (allocated(error)) return
         do mode = 1, 2
            evaluator%residual_mode = mode
            call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
            call check(error, allocated(merr) .and. .not. feasible, "Reject errors and non-finite residuals")
            if (allocated(error)) return
         end do
         evaluator%residual_mode = 0
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr), "Recover after a rejected input")
         if (allocated(error)) return
      end do
      call ctx%delete()
   end subroutine test_solver_failures

   !> Callback errors override a true feasibility flag, including recoverable trial errors
   !>
   !> @param[out] error Test failure information
   subroutine test_jfnk_callback_errors(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_jfnk_type) :: solver
      type(linear_evaluator_type) :: evaluator
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(2, 1), updated(2, 1), residual(2, 1), rmse
      real(wp) :: previous_step, previous_alpha
      integer :: mode
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      gamma = 0.0_wp
      do mode = 1, 3
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2))
         evaluator%jvp_mode = mode
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, .not. feasible .and. allocated(merr), "Reject unusable Jacobian products")
         if (allocated(error)) return
      end do
      evaluator%jvp_mode = 0
      evaluator%trial_limit = 0.0_wp
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, use_exact_jvp=.false., fd_eps=1.0_wp))
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, .not. feasible .and. allocated(merr), "A failed finite-difference probe must propagate its error")
      if (allocated(error)) return

      evaluator%trial_limit = huge(1.0_wp)
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, eta0=1.0e-12_wp))
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr))
      if (allocated(error)) return
      previous_step = solver%prev_step_maxabs
      previous_alpha = solver%prev_alpha
      evaluator%trial_limit = 0.0_wp
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, .not. feasible .and. allocated(merr), "Report failure when every trial is outside the domain")
      if (allocated(error)) return
      call check(error, solver%prev_step_maxabs == previous_step .and. solver%prev_alpha == previous_alpha, &
         & "A failed update must preserve the last accepted-step diagnostics")
      if (allocated(error)) return
      evaluator%trial_limit = 0.6_wp
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, eta0=1.0e-12_wp))
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. maxval(abs(updated)) <= evaluator%trial_limit)
      if (allocated(error)) return
      call evaluator%residual(updated, residual, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. sqrt(sum(residual**2)/2.0_wp) < rmse)
      call ctx%delete()
   end subroutine test_jfnk_callback_errors

   !> Non-diagonal left operators work independently of scaling and readiness
   !>
   !> @param[out] error Test failure information
   subroutine test_preconditioner_contract(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_jfnk_type) :: solver
      type(linear_evaluator_type) :: evaluator
      type(matrix_preconditioner_type), target :: precond
      real(wp), parameter :: scales(*) = [1.0_wp, 1.0e-20_wp, 1.0e20_wp, -2.0_wp]
      character(len=5), parameter :: allowed(*) = [character(len=5) :: "none", "block"]
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(2, 1), updated(2, 1), rmse
      integer :: i
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      gamma = 0.0_wp
      do i = 1, size(scales)
         precond%scale = scales(i)
         precond%napplied = 0
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, eta0=1.0e-12_wp))
         call solver%set_preconditioner(precond)
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr) .and. precond%napplied > 0)
         if (allocated(error)) return
         call check(error, maxval(abs(updated(:, 1) - evaluator%target)) < 1.0e-11_wp, &
            & "A scalar rescaling of the left operator must preserve the Newton solution")
         if (allocated(error)) return
      end do
      precond%enabled = .false.
      precond%napplied = 0
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, eta0=1.0e-12_wp))
      call solver%set_preconditioner(precond)
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. precond%napplied == 0 &
         & .and. maxval(abs(updated(:, 1) - evaluator%target)) < 1.0e-11_wp)
      if (allocated(error)) return
      precond%enabled = .true.
      precond%scale = 1.0_wp
      do i = 1, 2
         precond%output_mode = i
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2))
         call solver%set_preconditioner(precond)
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, .not. feasible .and. allocated(merr), "Reject unusable left operators")
         if (allocated(error)) return
      end do
      call check_preconditioner_name(" BLOCK ", allowed, "test model", merr)
      call check(error, .not. allocated(merr))
      if (allocated(error)) return
      call check_preconditioner_name("unknown", allowed, "test model", merr)
      call check(error, allocated(merr))
      if (allocated(error)) return
      call check_preconditioner_name("none", allowed(:0), "test model", merr)
      call check(error, allocated(merr), "An empty preconditioner list must report an error")
      call ctx%delete()
   end subroutine test_preconditioner_contract

   !> Reconfigure restart/history sizes without changing the field shape
   !>
   !> @param[out] error Test failure information
   subroutine test_workspace_reconfigure(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_jfnk_type) :: jfnk
      type(moist_model_moz_solver_diis_type) :: diis
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(2, 1), updated(2, 1), rmse
      integer, parameter :: sizes(*) = [2, 5, 1]
      integer :: i
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      gamma = 1.0_wp
      jfnk = new_jfnk_solver(ctx, jfnk_parameters(restart=3))
      call jfnk%ensure_workspace(gamma)
      do i = 1, size(sizes)
         jfnk%params%restart = sizes(i)
         call jfnk%ensure_workspace(gamma)
         call check(error, size(jfnk%workspace%ws_v, 2) == sizes(i) + 1 .and. size(jfnk%workspace%ws_h, 2) == sizes(i))
         if (allocated(error)) return
      end do
      jfnk%params%restart = 5
      do i = 2, 3
         jfnk%params%recycle = i
         call jfnk%ensure_workspace(gamma)
         call check(error, allocated(jfnk%workspace%ws_u))
         if (allocated(error)) return
         call check(error, size(jfnk%workspace%ws_u, 2) == i .and. size(jfnk%workspace%ws_bproj, 1) == i)
         if (allocated(error)) return
      end do
      jfnk%params%recycle = 0
      call jfnk%ensure_workspace(gamma)
      call check(error, .not. allocated(jfnk%workspace%ws_u))
      if (allocated(error)) return
      diis = new_diis_solver(ctx, diis_parameters(depth=2))
      do i = 1, size(sizes)
         diis%params%depth = sizes(i)
         call diis%outer_step(identity_evaluator_type(), gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr) .and. size(diis%fr, 2) == sizes(i))
         if (allocated(error)) return
      end do
      call ctx%delete()
   end subroutine test_workspace_reconfigure

   !> Preserve a sub-ulp residual when the resulting damped update is representable
   !>
   !> @param[out] error Test failure information
   subroutine test_diis_residual_precision(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_diis_type) :: solver
      type(fixed_evaluator_type) :: evaluator
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(1, 2), updated(1, 2), rmse
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      gamma = 1.0_wp
      evaluator%value = epsilon(1.0_wp)/4.0_wp
      solver = new_diis_solver(ctx, diis_parameters(depth=2), picard=picard_parameters(damping=4.0_wp))
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr))
      if (allocated(error)) return
      call check(error, all(updated == 1.0_wp - epsilon(1.0_wp)), &
         & "DIIS must use the evaluator residual without reconstructing it from rounded iterates")
      call ctx%delete()
   end subroutine test_diis_residual_precision

   !> Picard mixes the fixed-point map with the prescribed damping
   subroutine test_picard_update(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_picard_type) :: solver
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(3, 2), updated(3, 2), rmse
      logical :: feasible
      integer :: i

      call new_context(ctx, nthreads=0, verbosity=0)
      solver = new_picard_solver(ctx, picard_parameters(damping=0.25_wp))
      gamma = 2.0_wp
      call solver%outer_step(identity_evaluator_type(), gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr))
      if (allocated(error)) return
      call check(error, maxval(abs(updated - 1.5_wp)) < 1.0e-14_wp .and. &
         & abs(rmse - 2.0_wp) < 1.0e-14_wp .and. solver%prev_alpha == 0.25_wp)
      if (allocated(error)) return
      do i = 1, 100
         gamma = updated
         call solver%outer_step(identity_evaluator_type(), gamma, updated, rmse, feasible, merr)
      end do
      call check(error, maxval(abs(updated)) < 1.0e-12_wp)
      call ctx%delete()
   end subroutine test_picard_update

   !> Compare DIIS against analytic coefficients for orthogonal and correlated residuals
   subroutine test_diis_extrapolation(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_diis_type) :: solver
      type(fixed_evaluator_type) :: evaluator
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(1, 2), updated(1, 2), expected(2), rmse, damping
      logical :: feasible
      integer :: mode

      call new_context(ctx, nthreads=0, verbosity=0)
      do mode = 1, 3
         damping = 0.1_wp
         if (mode == 2) damping = 0.8_wp
         if (mode == 3) damping = 0.5_wp
         solver = new_diis_solver(ctx, diis_parameters(depth=2, damping=damping, step_limit=.false.), &
            & picard=picard_parameters(damping=0.5_wp))
         gamma = 0.0_wp
         evaluator%value = 0.0_wp
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         gamma(1, :) = [10.0_wp, 20.0_wp]
         evaluator%value(1, :) = [-1.0_wp, 0.0_wp]
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         gamma(1, :) = [-5.0_wp, 4.0_wp]
         evaluator%value(1, :) = [0.0_wp, -2.0_wp]
         expected = [7.8_wp, 17.2_wp]
         if (mode == 3) then
            evaluator%value(1, :) = [-0.5_wp, -2.0_wp]
            expected = (15.0_wp*[11.0_wp, 20.0_wp] + 2.0_wp*[-4.5_wp, 6.0_wp])/17.0_wp
         end if
         expected = gamma(1, :) + damping*(expected - gamma(1, :))
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         call check(error, maxval(abs(updated(1, :) - expected)) < 1.0e-12_wp, &
            & "DIIS must minimize the full residual Gram matrix")
         if (allocated(error)) return
      end do
      call ctx%delete()
   end subroutine test_diis_extrapolation

   !> DIIS converges on a coupled map and uses Picard mixing for singular histories
   subroutine test_diis_convergence(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_diis_type) :: solver
      type(fixed_evaluator_type) :: fixed
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(2, 1), updated(2, 1), flat(1, 2), mixed(1, 2), rmse
      logical :: feasible
      integer :: i, mode

      call new_context(ctx, nthreads=0, verbosity=0)
      do mode = 1, 2
         solver = new_diis_solver(ctx, diis_parameters(depth=2, damping=0.8_wp, step_limit=mode == 1))
         gamma = 0.0_wp
         do i = 1, 200
            call solver%outer_step(coupled_evaluator_type(), gamma, updated, rmse, feasible, merr)
            call check(error, feasible .and. .not. allocated(merr))
            if (allocated(error)) return
            if (rmse < 1.0e-10_wp) exit
            gamma = updated
         end do
         call check(error, maxval(abs(gamma(:, 1) - [1.0_wp, -2.0_wp])) < 1.0e-8_wp)
         if (allocated(error)) return
         call solver%reset()
         call check(error, solver%nstored == 0 .and. .not. allocated(solver%fr))
         if (allocated(error)) return
      end do

      ! A constant residual makes the DIIS system singular once the history is full;
      ! each such step falls back to Picard mixing
      solver = new_diis_solver(ctx, diis_parameters(depth=2), picard=picard_parameters(damping=0.5_wp))
      fixed%value(1, :) = [1.0_wp, -2.0_wp]
      flat = 0.0_wp
      do i = 1, 4
         call solver%outer_step(fixed, flat, mixed, rmse, feasible, merr)
         call check(error, maxval(abs(mixed - (flat - 0.5_wp*fixed%value))) < 1.0e-14_wp)
         if (allocated(error)) exit
         flat = mixed
      end do
      call ctx%delete()
   end subroutine test_diis_convergence

   !> A new equation clears DIIS residuals while preserving controls and storage
   subroutine test_diis_stages(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_diis_type), target :: solver
      class(moist_model_moz_solver_type), pointer :: driver
      type(fixed_evaluator_type) :: evaluator
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(1, 2), updated(1, 2), rmse
      logical :: feasible
      integer :: i

      call new_context(ctx, nthreads=0, verbosity=0)
      solver = new_diis_solver(ctx, diis_parameters(depth=2), picard=picard_parameters(damping=0.25_wp))
      driver => solver
      call driver%begin_stage()
      gamma = 0.0_wp
      evaluator%value(1, :) = [1.0_wp, -2.0_wp]
      do i = 1, 2
         call driver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         gamma = updated
      end do
      call check(error, solver%nstored == 2)
      if (allocated(error)) return
      call driver%begin_stage()
      call driver%begin_stage()
      call check(error, solver%nstored == 0 .and. allocated(solver%fr) .and. &
         & solver%prev_alpha == 0.0_wp .and. solver%prev_step_maxabs == 0.0_wp)
      if (allocated(error)) return
      evaluator%value(1, :) = [-0.3_wp, 0.8_wp]
      call driver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr))
      if (allocated(error)) return
      call check(error, maxval(abs(updated - (gamma - 0.25_wp*evaluator%value))) < 1.0e-14_wp .and. &
         & solver%nstored == 1, "new-stage DIIS update used residuals from the old equation")
      call ctx%delete()
   end subroutine test_diis_stages

   !> Residual growth clears a full DIIS history and returns a Picard update.
   !> @param[out] error Test failure.
   subroutine test_diis_divergence_restart(error)
      !> Test and library failures.
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_diis_type) :: solver
      type(fixed_evaluator_type) :: evaluator
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(1, 2), updated(1, 2), rmse
      logical :: feasible
      integer :: i

      call new_context(ctx, nthreads=0, verbosity=0)
      solver = new_diis_solver(ctx, diis_parameters(depth=2), picard=picard_parameters(damping=0.25_wp))
      gamma = 0.0_wp
      evaluator%value(1, :) = [1.0_wp, -2.0_wp]
      do i = 1, 2
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         gamma = updated
      end do
      evaluator%value = 100.0_wp*evaluator%value
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr))
      if (allocated(error)) return
      call check(error, solver%nstored == 0 .and. .not. allocated(solver%fr), &
         & "Divergent DIIS history must be released")
      if (allocated(error)) return
      call check(error, maxval(abs(updated - gamma + 0.25_wp*evaluator%value)) < 1.0e-13_wp, &
         & "Divergent DIIS must fall back to the configured Picard update")
      call ctx%delete()
   end subroutine test_diis_divergence_restart

   !> Frozen nonsymmetric transpose systems, including complex eigenvalues.
   !> @param[out] error Test failure.
   subroutine test_adjoint_solve(error)
      !> Test and library diagnostics.
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(linear_evaluator_type) :: evaluator
      type(matrix_preconditioner_type), target :: preconditioner
      !> A known solution provides a reference independent of GMRES.
      real(wp), parameter :: expected(2) = [0.75_wp, -1.25_wp]
      real(wp) :: rhs(2, 1)
      real(wp), allocatable :: solution(:, :)
      integer, target :: linearizations
      integer :: matrix_case, mode
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      linearizations = 0
      evaluator%linearizations => linearizations
      do matrix_case = 1, 2
         if (matrix_case == 2) evaluator%matrix = reshape([1.0_wp, 2.0_wp, -2.0_wp, 1.0_wp], [2, 2])
         rhs(:, 1) = matmul(transpose(evaluator%matrix), expected)
         do mode = 1, 2
            if (mode == 1) then
               call solve_adjoint(ctx, evaluator, rhs, solution, merr, tolerance=1.0e-11_wp)
            else
               call solve_adjoint(ctx, evaluator, rhs, solution, merr, tolerance=1.0e-11_wp, &
                  & preconditioner=preconditioner)
            end if
            call check(error, .not. allocated(merr))
            if (allocated(error)) return
            call check(error, allocated(solution))
            if (allocated(error)) return
            call check(error, maxval(abs(solution(:, 1) - expected)) < 1.0e-11_wp, &
               & "Adjoint must solve the transpose of a nonsymmetric matrix")
            if (allocated(error)) return
            call check(error, maxval(abs(matmul(transpose(evaluator%matrix), solution(:, 1)) - rhs(:, 1))) &
               & < 1.0e-11_wp)
            if (allocated(error)) return
         end do
      end do
      call check(error, linearizations == 0 .and. preconditioner%napplied > 0, &
         & "Adjoint must preserve the primal linearization and use its borrowed preconditioner")
      if (allocated(error)) return
      ! A converged guess is returned as it is; a nearby guess converges to the same solution
      call solve_adjoint(ctx, evaluator, rhs, solution, merr, tolerance=1.0e-11_wp, guess=reshape(expected, [2, 1]))
      call check(error, .not. allocated(merr) .and. allocated(solution))
      if (allocated(error)) return
      call check(error, all(solution(:, 1) == expected), "A converged adjoint guess must be kept")
      if (allocated(error)) return
      call solve_adjoint(ctx, evaluator, rhs, solution, merr, tolerance=1.0e-11_wp, &
         & guess=reshape(expected + [0.1_wp, -0.2_wp], [2, 1]))
      call check(error, .not. allocated(merr) .and. allocated(solution))
      if (allocated(error)) return
      call check(error, maxval(abs(solution(:, 1) - expected)) < 1.0e-11_wp, "Adjoint from a nearby guess")
      if (allocated(error)) return
      call solve_adjoint(ctx, evaluator, rhs, solution, merr, guess=reshape([1.0_wp, 2.0_wp, 3.0_wp], [3, 1]))
      call check(error, allocated(merr) .and. .not. allocated(solution), "Adjoint guess of the wrong shape")
      if (allocated(error)) return
      rhs = 0.0_wp
      call solve_adjoint(ctx, evaluator, rhs, solution, merr)
      call check(error, .not. allocated(merr))
      if (allocated(error)) return
      call check(error, allocated(solution))
      if (allocated(error)) return
      call check(error, all(solution == 0.0_wp))
      call ctx%delete()
   end subroutine test_adjoint_solve

   !> Unsupported operators, invalid controls and callback failures are explicit.
   !> @param[out] error Test failure.
   subroutine test_adjoint_failures(error)
      !> Test and library diagnostics.
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(linear_evaluator_type) :: evaluator
      type(identity_evaluator_type) :: unsupported
      type(moist_context_type), target :: ctx
      real(wp) :: rhs(2, 1)
      real(wp), allocatable :: solution(:, :)
      integer :: mode

      call new_context(ctx, nthreads=0, verbosity=0)
      rhs(:, 1) = [1.0_wp, -2.0_wp]
      call solve_adjoint(ctx, unsupported, rhs, solution, merr)
      call check(error, allocated(merr) .and. .not. allocated(solution))
      if (allocated(error)) return
      call solve_adjoint(ctx, evaluator, rhs, solution, merr, tolerance=0.0_wp)
      call check(error, allocated(merr) .and. .not. allocated(solution))
      if (allocated(error)) return
      call solve_adjoint(ctx, evaluator, rhs, solution, merr, max_iter=0)
      call check(error, allocated(merr) .and. .not. allocated(solution))
      if (allocated(error)) return
      call solve_adjoint(ctx, evaluator, rhs, solution, merr, parameters=adjoint_parameters(restart=0))
      call check(error, allocated(merr) .and. .not. allocated(solution), "Invalid adjoint restart")
      if (allocated(error)) return
      call solve_adjoint(ctx, evaluator, rhs, solution, merr, parameters=adjoint_parameters(max_krylov=0))
      call check(error, allocated(merr) .and. .not. allocated(solution), "Invalid adjoint Krylov budget")
      if (allocated(error)) return
      rhs(1, 1) = ieee_value(0.0_wp, ieee_quiet_nan)
      call solve_adjoint(ctx, evaluator, rhs, solution, merr)
      call check(error, allocated(merr) .and. .not. allocated(solution))
      if (allocated(error)) return
      rhs(:, 1) = [1.0_wp, -2.0_wp]
      do mode = 1, 3
         evaluator%jvp_mode = mode
         call solve_adjoint(ctx, evaluator, rhs, solution, merr, max_iter=1)
         call check(error, allocated(merr) .and. .not. allocated(solution), &
            & "An unusable transpose callback must not return a solution")
         if (allocated(error)) return
      end do
      call ctx%delete()
   end subroutine test_adjoint_failures

   !> Shared scratch preserves storage, frozen operators and the next primal update.
   !> @param[out] error Test failure
   subroutine test_adjoint_shared_workspace(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      class(moist_model_moz_solver_type), allocatable, target :: owner, fresh
      type(moist_model_moz_solver_parameters_type) :: parameters
      type(moist_model_moz_solver_statistics_type) :: before, after
      type(linear_evaluator_type) :: evaluator
      type(c_ptr) :: first, address
      real(wp) :: gamma(2, 1), updated(2, 1), reference(2, 1), rhs(2, 1), expected(2), rms, ref_rms
      real(wp), allocatable :: solution(:, :)
      character(len=6), parameter :: keys(*) = [character(len=6) :: "JFNK", "HYBRID", "Picard"]
      type(moist_context_type), target :: ctx
      integer, target :: linearizations
      integer :: mode, repeat, frozen, capacity
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      parameters%hybrid%switch_rmse = 100.0_wp
      parameters%jfnk%recycle = 2
      linearizations = 0
      evaluator%linearizations => linearizations
      do mode = 1, size(keys)
         parameters%name = keys(mode)
         call new_solver(owner, ctx, merr, parameters)
         call check(error, .not. allocated(merr))
         if (allocated(error)) return
         gamma = 0.0_wp
         call owner%outer_step(evaluator, gamma, updated, rms, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         first = c_null_ptr
         if (mode /= 3) call krylov_address(owner, first, capacity)
         do repeat = 1, 2
            expected = [0.75_wp, -1.25_wp] + real(repeat - 1, wp)*[0.1_wp, 0.2_wp]
            rhs(:, 1) = matmul(transpose(evaluator%matrix), expected)
            frozen = linearizations
            before = owner%get_statistics()
            call solve_adjoint(ctx, evaluator, rhs, solution, merr, tolerance=1.0e-11_wp, workspace_owner=owner)
            call check(error, .not. allocated(merr))
            if (allocated(error)) return
            call check(error, maxval(abs(solution(:, 1) - expected)) < 1.0e-11_wp)
            if (allocated(error)) return
            after = owner%get_statistics()
            call check(error, frozen == linearizations .and. after%nkrylov == before%nkrylov &
               & .and. after%njvp == before%njvp .and. after%nresidual == before%nresidual, &
               & "Transpose solve must retain primal settings, counters and linearization")
            if (allocated(error)) return
            call krylov_address(owner, address, capacity)
            call check(error, c_associated(address) .and. capacity == 51)
            if (allocated(error)) return
            if (c_associated(first)) then
               call check(error, c_associated(first, address), "Response must return the same Krylov allocation")
               if (allocated(error)) return
            end if
            first = address
         end do
         parameters%name = keys(mode)
         call new_solver(fresh, ctx, merr, parameters)
         call check(error, .not. allocated(merr))
         if (allocated(error)) return
         gamma(:, 1) = [3.0_wp, -4.0_wp]
         call owner%begin_stage()
         call owner%outer_step(evaluator, gamma, updated, rms, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         call fresh%outer_step(evaluator, gamma, reference, ref_rms, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         call check(error, maxval(abs(updated - reference)) < 1.0e-12_wp .and. abs(rms - ref_rms) < 1.0e-12_wp, &
            & "Primal update after response must match a fresh solver")
         if (allocated(error)) return
         call krylov_address(owner, address, capacity)
         call check(error, c_associated(first, address), "Next primal stage must reuse the allocation")
         if (allocated(error)) return
      end do
      call ctx%delete()
   end subroutine test_adjoint_shared_workspace

   !> Callback failures and differing restarts return usable, independent storage.
   !> @param[out] error Test failure
   subroutine test_adjoint_workspace_recovery(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_jfnk_type), target :: owner, independent
      type(linear_evaluator_type) :: evaluator
      type(c_ptr) :: first, address, other
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(2, 1), updated(2, 1), rhs(2, 1), rms, small(3, 1)
      real(wp), allocatable :: solution(:, :)
      integer :: mode, capacity
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      owner = new_jfnk_solver(ctx, jfnk_parameters(restart=5, max_krylov=20, recycle=2, eta0=1.0e-10_wp))
      gamma = 0.0_wp
      call owner%outer_step(evaluator, gamma, updated, rms, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr))
      if (allocated(error)) return
      rhs(:, 1) = matmul(transpose(evaluator%matrix), evaluator%target)
      call solve_adjoint(ctx, evaluator, rhs, solution, merr, tolerance=1.0e-11_wp, workspace_owner=owner)
      call check(error, .not. allocated(merr))
      if (allocated(error)) return
      call krylov_address(owner, first, capacity)
      call check(error, capacity == 51 .and. owner%params%restart == 5 .and. .not. owner%params%cgs2, &
         & "Shared capacity must preserve the distinct numerical settings")
      if (allocated(error)) return
      independent = new_jfnk_solver(ctx, jfnk_parameters(restart=5))
      call independent%ensure_workspace(gamma)
      call krylov_address(independent, other, capacity)
      call check(error, .not. c_associated(first, other), "Independent models must own distinct scratch")
      if (allocated(error)) return
      call solve_adjoint(ctx, evaluator, rhs, solution, merr, tolerance=1.0e-11_wp, workspace_owner=independent, &
         & parameters=adjoint_parameters(restart=3, max_krylov=8))
      call check(error, .not. allocated(merr) .and. allocated(solution))
      if (allocated(error)) return
      call check(error, maxval(abs(matmul(transpose(evaluator%matrix), solution(:, 1)) - rhs(:, 1))) < 1.0e-11_wp)
      if (allocated(error)) return
      call krylov_address(independent, address, capacity)
      call check(error, capacity == 6 .and. c_associated(address, other), &
         & "Configured response restart must not grow a sufficient primal basis")
      if (allocated(error)) return
      evaluator%jvp_mode = 1
      call solve_adjoint(ctx, evaluator, rhs, solution, merr, workspace_owner=independent, &
         & parameters=adjoint_parameters(restart=3, max_krylov=8))
      call check(error, allocated(merr) .and. .not. allocated(solution))
      if (allocated(error)) return
      call krylov_address(independent, address, capacity)
      call check(error, capacity == 6 .and. c_associated(address, other), &
         & "Failed configured adjoint must return its basis")
      if (allocated(error)) return
      do mode = 1, 3
         evaluator%jvp_mode = mode
         call solve_adjoint(ctx, evaluator, rhs, solution, merr, max_iter=1, workspace_owner=owner)
         call check(error, allocated(merr) .and. .not. allocated(solution))
         if (allocated(error)) return
         call krylov_address(owner, address, capacity)
         call check(error, c_associated(first, address), "Failed response must return the allocation")
         if (allocated(error)) return
      end do
      evaluator%jvp_mode = 0
      do mode = 1, 2
         owner%params%restart = 2
         call owner%begin_stage()
         call owner%outer_step(evaluator, gamma, updated, rms, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         call check(error, maxval(abs(updated(:, 1) - evaluator%target)) < 1.0e-10_wp)
         if (allocated(error)) return
         call solve_adjoint(ctx, evaluator, rhs, solution, merr, tolerance=1.0e-11_wp, workspace_owner=owner)
         call check(error, .not. allocated(merr) .and. owner%params%restart == 2 .and. owner%nrecycled == 0)
         if (allocated(error)) return
         call krylov_address(owner, address, capacity)
         call check(error, c_associated(first, address) .and. capacity == 51, &
            & "Differing restarts and recycling must not repeatedly resize the basis")
         if (allocated(error)) return
      end do
      small = 0.0_wp
      call owner%ensure_workspace(small)
      call check(error, owner%workspace%ws_nflat == 3 .and. size(owner%workspace%ws_v, 2) == 3, &
         & "Changed field shape must replace incompatible scratch")
      call ctx%delete()
   end subroutine test_adjoint_workspace_recovery

   !> Observe the allocation through ownership transfer, including Hybrid's child.
   !> @param[inout] solver Solver between completed stages
   !> @param[out] address Opaque address, never dereferenced
   !> @param[out] capacity Allocated basis columns
   subroutine krylov_address(solver, address, capacity)
      class(moist_model_moz_solver_type), intent(inout) :: solver
      type(c_ptr), intent(out) :: address
      integer, intent(out) :: capacity
      type(moist_model_moz_solver_workspace_type), allocatable, target :: scratch
      call solver%exchange_workspace(scratch)
      address = c_null_ptr
      capacity = 0
      if (allocated(scratch)) then
         if (allocated(scratch%ws_v)) then
            address = c_loc(scratch%ws_v)
            capacity = size(scratch%ws_v, 2)
         end if
      end if
      call solver%exchange_workspace(scratch)
   end subroutine krylov_address

   !> Finite differences and exact products reach the same root
   subroutine test_jfnk_nonlinear(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_jfnk_type) :: solver
      type(quadratic_evaluator_type) :: evaluator
      type(quadratic_linearization_type), target :: state
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(1, 1), updated(1, 1), rmse
      logical :: feasible
      integer :: mode, i

      call new_context(ctx, nthreads=0, verbosity=0)
      evaluator%state => state
      do mode = 1, 2
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=4, use_exact_jvp=mode == 2))
         gamma = 2.0_wp
         do i = 1, 30
            call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
            call check(error, feasible .and. .not. allocated(merr))
            if (allocated(error)) return
            if (rmse < 1.0e-10_wp) exit
            gamma = updated
         end do
         call check(error, abs(gamma(1, 1) - 1.0_wp) < 1.0e-9_wp)
         if (allocated(error)) return
         call check(error, (solver%njvp_total > 0) .eqv. (mode == 2))
         if (allocated(error)) return
         call check(error, allocated(solver%workspace%ws_f_pert) .eqv. (mode == 1), &
            & "FD scratch must be lazy")
         if (allocated(error)) return
      end do
      call ctx%delete()
   end subroutine test_jfnk_nonlinear

   !> JFNK with exact products and a tight fixed forcing term is exact Newton: quadratic convergence
   !>
   !> - F(x) = x**2 - 1 gives F(x_new) = F(x)**2/(4 x**2) <= F(x)**2 for x >= 1
   !>
   !> @param[out] error Test failure information
   subroutine test_jfnk_exact_newton(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(mctc_error), allocatable :: merr
      !> Solver under test
      type(moist_model_moz_solver_jfnk_type) :: solver
      !> Evaluator with exact products and its linearization
      type(quadratic_evaluator_type) :: evaluator
      type(quadratic_linearization_type), target :: state
      !> Quiet run context borrowed by the solver
      type(moist_context_type), target :: ctx
      !> Iterate, update and residuals of this and the previous step
      real(wp) :: gamma(1, 1), updated(1, 1), rmse, rmse_prev
      !> Step usability
      logical :: feasible
      !> Newton step counter
      integer :: i

      call new_context(ctx, nthreads=0, verbosity=0)
      evaluator%state => state
      solver = new_jfnk_solver(ctx, jfnk_parameters(use_exact_jvp=.true., eta0=1.0e-8_wp, eta_min=1.0e-8_wp, eta_max=1.0e-8_wp, &
         & alpha_max=1.0_wp))
      gamma = 2.0_wp
      rmse_prev = huge(1.0_wp)
      do i = 1, 10
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         if (i > 1) then
            call check(error, rmse <= rmse_prev**2*(1.0_wp + 1.0e-6_wp), "Newton steps must converge quadratically")
            if (allocated(error)) return
         end if
         if (rmse < 1.0e-12_wp) exit
         call check(error, solver%prev_alpha == 1.0_wp, "Newton steps must be full steps")
         if (allocated(error)) return
         rmse_prev = rmse
         gamma = updated
      end do
      call check(error, i <= 6 .and. abs(gamma(1, 1) - 1.0_wp) < 1.0e-12_wp)
      if (allocated(error)) return
      call check(error, solver%njvp_total > 0 .and. .not. allocated(solver%workspace%ws_f_pert), &
         & "Newton must use exact products only")
      call ctx%delete()
   end subroutine test_jfnk_exact_newton

   !> With jvp = exact, steps fail instead of falling back to finite-difference products
   !>
   !> @param[out] error Test failure information
   subroutine test_jfnk_exact_required(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(mctc_error), allocatable :: merr
      !> Solver under test
      type(moist_model_moz_solver_jfnk_type) :: solver
      !> Evaluator without products
      type(diagonal_evaluator_type) :: plain
      !> Evaluator whose linearization is refused
      type(quadratic_evaluator_type) :: refusing
      type(quadratic_linearization_type), target :: state
      !> Quiet run context borrowed by the solver
      type(moist_context_type), target :: ctx
      !> Iterate and update
      real(wp) :: gamma(1, 1), updated(1, 1), rmse
      !> Step usability
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      solver = new_jfnk_solver(ctx, jfnk_parameters(use_exact_jvp=.true.))
      gamma = 2.0_wp
      call solver%outer_step(plain, gamma, updated, rmse, feasible, merr)
      call check(error, .not. feasible .and. allocated(merr), "Evaluator without exact products")
      if (allocated(error)) return

      refusing%state => state
      refusing%accept_linearization = .false.
      call solver%outer_step(refusing, gamma, updated, rmse, feasible, merr)
      call check(error, .not. feasible .and. allocated(merr), "Refused exact linearization")
      if (allocated(error)) return
      call check(error, index(merr%message, "Deliberately unavailable") > 0, "Linearization error must be reported")
      if (allocated(error)) return
      call check(error, solver%nresidual_total == 2 .and. .not. allocated(solver%workspace%ws_f_pert), &
         & "No finite-difference probes may run")
      call ctx%delete()
   end subroutine test_jfnk_exact_required

   !> Validation rejects a negative recycling dimension and an expansion cap below one
   !>
   !> @param[out] error Test failure information
   subroutine test_jfnk_validate(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Controls under test
      type(jfnk_parameters) :: p
      !> Library error
      type(mctc_error), allocatable :: merr

      call p%validate(merr)
      call check(error, .not. allocated(merr), "Default controls must be valid")
      if (allocated(error)) return

      p%recycle = -1
      call p%validate(merr)
      call check(error, allocated(merr), "Negative recycling dimension must be rejected")
      if (allocated(error)) return

      call p%init_defaults()
      p%alpha_max = 0.5_wp
      call p%validate(merr)
      call check(error, allocated(merr), "Expansion cap below one must be rejected")
   end subroutine test_jfnk_validate

   !> Borrowed preconditioners survive workspace allocation, replacement and release
   subroutine test_jfnk_preconditioner(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_jfnk_type) :: solver
      type(diagonal_preconditioner_type), target :: precond, replacement
      type(moist_context_type), target :: ctx
      real(wp) :: gamma(7, 2), updated(7, 2), rmse
      logical :: feasible
      integer :: i, j

      call new_context(ctx, nthreads=0, verbosity=0)
      gamma = 0.0_wp
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=4, max_krylov=8, use_exact_jvp=.false., fd_eps=1.0_wp, eta0=1.0e-9_wp))
      call solver%set_preconditioner(precond)
      call solver%outer_step(diagonal_evaluator_type(), gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. precond%napplied > 0)
      if (allocated(error)) return
      do j = 1, 2
         do i = 1, 7
            call check(error, abs(updated(i, j) - target_value(i, j)) < 1.0e-10_wp)
            if (allocated(error)) return
         end do
      end do
      call solver%release_workspace()
      call check(error, associated(solver%precond, precond))
      if (allocated(error)) return
      call solver%set_preconditioner(replacement)
      call check(error, associated(solver%precond, replacement))
      if (allocated(error)) return
      call solver%set_fallback_filter(precond)
      call check(error, associated(solver%smoother, precond))
      if (allocated(error)) return
      call solver%set_fallback_filter()
      call solver%set_preconditioner()
      call check(error, .not. associated(solver%precond) .and. .not. associated(solver%smoother))
      call ctx%delete()
   end subroutine test_jfnk_preconditioner

   !> Backtracking handles both increasing merit and infeasible trial fields
   subroutine test_jfnk_line_search(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_jfnk_type) :: solver
      type(quadratic_evaluator_type) :: evaluator
      type(quadratic_linearization_type), target :: state
      real(wp) :: gamma(1, 1), updated(1, 1), rmse
      type(linear_evaluator_type) :: linear
      type(matrix_preconditioner_type), target :: precond, filter
      real(wp) :: field(2, 1), next(2, 1)
      type(moist_context_type), target :: ctx
      logical :: feasible
      integer :: mode

      call new_context(ctx, nthreads=0, verbosity=0)
      evaluator%state => state
      do mode = 1, 2
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=4))
         evaluator%max_argument = huge(1.0_wp)
         if (mode == 2) evaluator%max_argument = 2.0_wp
         gamma = 0.1_wp
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         call check(error, solver%prev_alpha > 0.0_wp .and. solver%prev_alpha < 1.0_wp .and. &
            & abs(updated(1, 1)**2 - 1.0_wp) < abs(gamma(1, 1)**2 - 1.0_wp))
         if (allocated(error)) return
      end do

      ! This invertible left operator gives a one-vector GMRES direction
      ! that increases the physical merit. As in pysofft the residual replaces
      ! it and is searched first; for this identity Jacobian it is the exact
      ! Newton step, so the Picard fallback and its filter are not needed.
      linear%matrix = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [2, 2])
      linear%target = [-1.0_wp, 0.0_wp]
      precond%matrix = reshape([1.0_wp, 1.0_wp, 0.0_wp, -3.0_wp], [2, 2])
      filter%matrix = linear%matrix
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=1, max_krylov=1))
      call solver%set_preconditioner(precond)
      call solver%set_fallback_filter(filter)
      field = 0.0_wp
      call solver%outer_step(linear, field, next, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. filter%napplied == 0)
      if (allocated(error)) return
      call check(error, solver%prev_alpha == 1.0_wp .and. maxval(abs(next(:, 1) - linear%target)) < 1.0e-14_wp, &
         & "An ascent direction must be replaced by a residual line search")
      if (allocated(error)) return

      ! The full Newton step leaves the allowed domain and the single Newton
      ! trial fails, so the filtered residual (Picard) step is taken at 1/4
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=1, max_krylov=1))
      solver%params%line_search_max_iter = 1
      call solver%set_fallback_filter(filter)
      linear%trial_limit = 0.5_wp
      call solver%outer_step(linear, field, next, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. filter%napplied == 1 .and. &
         & solver%prev_alpha == solver%params%fallback_alpha .and. &
         & maxval(abs(next(:, 1) - 0.25_wp*linear%target)) < 1.0e-14_wp, &
         & "A failed Newton search must use the filtered Picard fallback")
      call ctx%delete()
   end subroutine test_jfnk_line_search

   !> Lag the inner tolerance by one update, tighten it after stalls and detect a root
   !>
   !> @param[out] error Test failure information
   subroutine test_jfnk_forcing(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Model failure information
      type(mctc_error), allocatable :: merr
      !> One Krylov step leaves a nonzero residual for the next Newton update
      type(moist_model_moz_solver_jfnk_type), target :: solver
      !> Dense linear problem and borrowed operators
      type(linear_evaluator_type) :: evaluator
      type(matrix_preconditioner_type), target :: precond, filter
      !> Fields and successive input residual norms
      real(wp) :: gamma(2, 1), updated(2, 1), rmse, previous, expected
      !> Application count before a converged input
      integer :: calls
      !> Stalled outer steps, and zero steps returned at eta_min
      integer :: nstall, nmin
      logical :: feasible
      !> Silent run context borrowed by the solver
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=1, eta_max=0.9_wp))
      precond%forcing => solver%eta_current
      call solver%set_preconditioner(precond)
      gamma = 0.0_wp
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. precond%applied_eta == solver%params%eta0)
      if (allocated(error)) return
      ! As in pysofft the forcing term lags one update: the second inner solve
      ! still uses eta0, the third uses the ratio of the first two residuals
      previous = rmse
      gamma = updated
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. rmse < previous .and. &
         & precond%applied_eta == solver%params%eta0, "The second inner solve keeps eta0")
      if (allocated(error)) return
      expected = solver%params%forcing_beta*(rmse/previous)**solver%params%forcing_power
      expected = max(solver%params%eta_min, min(solver%params%eta_max, expected))
      gamma = updated
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. abs(precond%applied_eta - expected) < 1.0e-14_wp, &
         & "Later inner solves use the ratio of the previous two input residuals")
      if (allocated(error)) return

      calls = precond%napplied
      gamma(:, 1) = evaluator%target
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. rmse == 0.0_wp &
         & .and. all(updated == gamma) .and. solver%gmres_iters_last == 0 .and. precond%napplied == calls)
      if (allocated(error)) return

      ! Every Newton trial leaves the allowed domain. A zero fallback filter
      ! produces feasible but unchanged fields, exercising the stall control.
      call solver%begin_stage()
      solver%params%line_search_max_iter = 1
      evaluator%trial_limit = 0.1_wp
      filter%output_mode = 2
      call solver%set_fallback_filter(filter)
      gamma = 0.0_wp
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. solver%prev_alpha == 0.0_wp)
      if (allocated(error)) return
      expected = max(solver%params%eta_min, 0.1_wp*precond%applied_eta)
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr) .and. abs(precond%applied_eta - expected) < 1.0e-14_wp, &
         & "A stall must tighten the next inner solve")
      if (allocated(error)) return

      ! Further stalls tighten down to eta_min. The first zero step there is
      ! returned, so a driver can still accept a converged input; the second
      ! cannot tighten and ends the stage as an infeasible step without an error
      nmin = 0
      do nstall = 1, 8
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         if (allocated(merr) .or. .not. feasible) exit
         if (precond%applied_eta == solver%params%eta_min) nmin = nmin + 1
      end do
      call check(error, .not. allocated(merr) .and. .not. feasible .and. nstall <= 8 .and. nmin == 1 &
         & .and. solver%eta_current == solver%params%eta_min .and. precond%applied_eta == solver%params%eta_min, &
         & "One zero step at the lowest forcing term is returned, a second ends the stage unconverged")
      call ctx%delete()
   end subroutine test_jfnk_forcing

   !> Compare DEF1 Newton updates with the dense root, including a coarse-only RHS
   !>
   !> @param[out] error Test failure information
   subroutine test_jfnk_deflation(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Model failure information
      type(mctc_error), allocatable :: merr
      !> Solver, dense evaluator and noncommuting left preconditioner
      type(moist_model_moz_solver_jfnk_type) :: solver
      type(linear_evaluator_type) :: evaluator
      type(matrix_preconditioner_type), target :: precond
      !> Driver-provided basis and weighted basis
      real(wp) :: profile(2, 1), weighted(2, 1)
      !> Fields and residual norm
      real(wp) :: gamma(2, 1), updated(2, 1), rmse
      !> Quiet run context borrowed by the solver
      type(moist_context_type), target :: ctx
      integer :: mode
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      profile(:, 1) = [1.0_wp, 0.0_wp]
      weighted(:, 1) = [2.0_wp, 0.0_wp]
      gamma = 0.0_wp
      do mode = 1, 3
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, eta0=1.0e-12_wp, &
            & deflation=.true., deflation_rmse=10.0_wp))
         call solver%deflation%init_basis(profile, weighted, 1)
         if (mode > 1) call solver%set_preconditioner(precond)
         if (mode == 3) evaluator%target = profile(:, 1)
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr) .and. solver%deflation%active)
         if (allocated(error)) return
         call check(error, maxval(abs(updated(:, 1) - evaluator%target)) < 1.0e-11_wp, &
            & "DEF1 must reconstruct the full Newton direction after the projected solve")
         if (allocated(error)) return
         if (mode == 3) then
            call check(error, solver%gmres_iters_last == 0, "A coarse-only RHS needs no Arnoldi iterations")
            if (allocated(error)) return
         end if
      end do
      call ctx%delete()
   end subroutine test_jfnk_deflation

   !> A configured deflation that cannot take effect is rejected or counted for the stage summary
   !>
   !> - Deflation with finite-difference products fails validation
   !> - Nothing is counted while the residual stays above deflation_rmse
   !> - An evaluator without a coarse space, a failed or a singular coarse matrix
   !>   leave the step undeflated and are counted until the next stage
   !>
   !> @param[out] error Test failure information
   subroutine test_jfnk_deflation_warnings(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Model failure information
      type(mctc_error), allocatable :: merr
      !> Solver and the dense evaluator with a coarse space
      type(moist_model_moz_solver_jfnk_type) :: solver
      type(linear_evaluator_type) :: evaluator
      !> Evaluator with exact products but no coarse space
      type(quadratic_evaluator_type) :: quadratic
      type(quadratic_linearization_type), target :: state
      !> Driver basis and weights; zero weights make the coarse matrix singular
      real(wp) :: profile(2, 1), weighted(2, 1)
      !> Fields and residual norm
      real(wp) :: gamma(2, 1), updated(2, 1), point(1, 1), next(1, 1), rmse
      logical :: feasible
      !> Silent run context
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      call run_checks()
      call ctx%delete()

   contains

      !> Checks in order; the first failure ends them
      subroutine run_checks()
         !> Deflation controls with finite-difference products
         type(jfnk_parameters) :: params
         !> Validation failure
         type(mctc_error), allocatable :: invalid

         params = jfnk_parameters(use_exact_jvp=.false., deflation=.true.)
         call params%validate(invalid)
         call check(error, allocated(invalid), "Deflation with finite-difference products is rejected")
         if (allocated(error)) return

         gamma = 0.0_wp
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, deflation=.true., &
            & deflation_rmse=1.0e-30_wp))
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr) .and. solver%nundeflated_stage == 0, &
            & "Above deflation_rmse the coarse space is not due, so nothing is counted")
         if (allocated(error)) return

         quadratic%state => state
         point = 2.0_wp
         solver = new_jfnk_solver(ctx, jfnk_parameters(deflation=.true., deflation_rmse=10.0_wp))
         call solver%outer_step(quadratic, point, next, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr) .and. .not. solver%deflation%active .and. &
            & solver%nundeflated_stage == 1, "An evaluator without deflation_setup leaves the step undeflated")
         if (allocated(error)) return
         call solver%begin_stage()
         call check(error, solver%nundeflated_stage == 0, "A new stage clears the count")
         if (allocated(error)) return

         profile(:, 1) = [1.0_wp, 0.0_wp]
         weighted = 0.0_wp
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, eta0=1.0e-12_wp, deflation=.true., &
            & deflation_rmse=10.0_wp))
         call solver%deflation%init_basis(profile, weighted, 1)
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr) .and. .not. solver%deflation%active .and. &
            & solver%nundeflated_stage == 1, "A singular coarse matrix leaves the step undeflated")
         if (allocated(error)) return

         ! Without a driver basis the test evaluator reports a coarse-matrix failure
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, eta0=1.0e-12_wp, deflation=.true., &
            & deflation_rmse=10.0_wp))
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr) .and. .not. solver%deflation%active .and. &
            & solver%nundeflated_stage == 1, "A coarse-matrix failure leaves the step undeflated")
      end subroutine run_checks

   end subroutine test_jfnk_deflation_warnings

   !> Recycling that cannot be honoured is capped at construction, and deflated
   !> solves that skip it are counted for the stage summary
   !>
   !> @param[out] error Test failure information
   subroutine test_jfnk_recycle_warnings(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Model failure information
      type(mctc_error), allocatable :: merr
      !> Solver and the dense evaluator with a coarse space
      type(moist_model_moz_solver_jfnk_type) :: solver
      type(linear_evaluator_type) :: evaluator
      !> Driver basis and weights
      real(wp) :: profile(2, 1), weighted(2, 1)
      !> Fields and residual norm
      real(wp) :: gamma(2, 1), updated(2, 1), rmse
      logical :: feasible
      !> Silent run context
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      call run_checks()
      call ctx%delete()

   contains

      !> Checks in order; the first failure ends them
      subroutine run_checks()
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=4, max_krylov=8, recycle=10))
         call check(error, solver%recycle_dim(), 3, "restart caps the recycle dimension")
         if (allocated(error)) return
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=1, max_krylov=8, recycle=2))
         call check(error, solver%recycle_dim(), 0, "a single-step cycle leaves no room for recycling")
         if (allocated(error)) return
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=4, max_krylov=8, recycle=3))
         call check(error, solver%recycle_dim(), 3, "a recycle dimension within the cycle is kept")
         if (allocated(error)) return

         profile(:, 1) = [1.0_wp, 0.0_wp]
         weighted(:, 1) = [2.0_wp, 0.0_wp]
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, recycle=1, eta0=1.0e-12_wp, &
            & deflation=.true., deflation_rmse=10.0_wp))
         call solver%deflation%init_basis(profile, weighted, 1)
         gamma = 0.0_wp
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr) .and. solver%deflation%active .and. &
            & solver%nrecycled == 0 .and. solver%nunrecycled_stage > 0, "A deflated stage skips recycling")
         if (allocated(error)) return
         call solver%begin_stage()
         call check(error, solver%nunrecycled_stage == 0, "A new stage clears the count")
         if (allocated(error)) return
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, solver%nunrecycled_stage > 0, "A new deflated stage counts again")
      end subroutine run_checks

   end subroutine test_jfnk_recycle_warnings

   !> The stage summary counts fallbacks and is delegated by HYBRID
   !>
   !> - An inner solve capped at max_krylov short of its tolerance
   !> - A Picard fallback that ends in a zero step
   !> - Empty after begin_stage
   !>
   !> @param[out] error Test failure information
   subroutine test_jfnk_stage_summary(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Model failure information
      type(mctc_error), allocatable :: merr
      !> JFNK and HYBRID solvers
      type(moist_model_moz_solver_jfnk_type) :: solver
      class(moist_model_moz_solver_type), allocatable :: hybrid
      type(moist_model_moz_solver_parameters_type) :: p
      !> Soft test problem, dense evaluator and a zero fallback filter
      type(soft_evaluator_type) :: soft
      type(soft_linearization_type), target :: state
      type(linear_evaluator_type) :: evaluator
      type(matrix_preconditioner_type), target :: filter
      !> Fields and residual norm
      real(wp) :: field(200, 1), field_new(200, 1), gamma(2, 1), updated(2, 1), rmse
      !> Summary text
      character(len=:), allocatable :: text
      logical :: feasible
      !> Silent run context
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      call run_checks()
      call ctx%delete()

   contains

      !> Checks in order; the first failure ends them
      subroutine run_checks()
         ! A tight forcing term makes the first inner solve use its whole budget
         soft%state => state
         allocate (state%point(size(field)))
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=12, max_krylov=8, eta0=1.0e-8_wp))
         field = 0.0_wp
         call solver%outer_step(soft, field, field_new, rmse, feasible, merr)
         text = solver%stage_summary()
         call check(error, feasible .and. .not. allocated(merr) .and. text == "JFNK: 1 inner solve hit max_krylov", &
            & "A capped inner solve is counted: "//text)
         if (allocated(error)) return
         call solver%begin_stage()
         call check(error, len(solver%stage_summary()) == 0, "A new stage starts with an empty summary")
         if (allocated(error)) return

         ! Newton trials leave the domain; the zero filter makes the Picard step a zero step
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=2, line_search_max_iter=1))
         evaluator%trial_limit = 0.1_wp
         filter%output_mode = 2
         call solver%set_fallback_filter(filter)
         gamma = 0.0_wp
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         text = solver%stage_summary()
         call check(error, feasible .and. .not. allocated(merr) .and. solver%prev_alpha == 0.0_wp .and. &
            & text == "JFNK: 1 of 1 steps used the Picard fallback, 1 with a zero step", &
            & "A Picard fallback and its zero step are counted: "//text)
         if (allocated(error)) return

         ! HYBRID reports its JFNK phase only
         evaluator%trial_limit = huge(1.0_wp)
         p%hybrid%switch_rmse = 10.0_wp
         p%jfnk%eta0 = 1.0e-12_wp
         p%jfnk%restart = 1
         p%jfnk%max_krylov = 1
         p%name = "HYBRID"
         call new_solver(hybrid, ctx, merr, p)
         call check(error, .not. allocated(merr), "HYBRID construction")
         if (allocated(error)) return
         call check(error, len(hybrid%stage_summary()) == 0, "No JFNK phase, no summary")
         if (allocated(error)) return
         call hybrid%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         text = hybrid%stage_summary()
         call check(error, feasible .and. .not. allocated(merr) .and. index(text, "JFNK: ") == 1 .and. &
            & index(text, "1 inner solve hit max_krylov") > 0, "HYBRID delegates to its JFNK phase: "//text)
         if (allocated(error)) return
         call hybrid%begin_stage()
         call check(error, len(hybrid%stage_summary()) == 0, "A new HYBRID stage starts in MDIIS without a summary")
      end subroutine run_checks

   end subroutine test_jfnk_stage_summary

   !> Preserve convergence, checkpoint ordering and failure semantics of the shared loop
   subroutine test_iteration_controller(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(iteration_context_type) :: context
      !> Quiet run context whose timer records the advance time
      type(moist_context_type), target :: ctx
      real(wp) :: rmse
      integer :: niter
      logical :: converged

      call new_context(ctx, nthreads=0, verbosity=0)
      rmse = 7.0_wp
      call run_iterations(context, ctx, 0, 0.3_wp, 0, rmse, niter, converged, merr)
      call check(error, niter == 0 .and. .not. converged .and. rmse == 7.0_wp)
      if (allocated(error)) return
      call run_iterations(context, ctx, 10, 0.3_wp, 2, rmse, niter, converged, merr)
      call check(error, converged .and. niter == 4 .and. context%nobserved == 4 .and. &
         & context%ncheckpoints == 1 .and. .not. allocated(merr))
      if (allocated(error)) return
      context = iteration_context_type()
      call run_iterations(context, ctx, 10, 0.3_wp, 2, rmse, niter, converged, merr, &
         & checkpoint_before_convergence=.true.)
      call check(error, converged .and. context%ncheckpoints == 2)
      if (allocated(error)) return
      context = iteration_context_type(fail_at=2)
      call run_iterations(context, ctx, 10, 0.3_wp, 1, rmse, niter, converged, merr)
      call check(error, .not. converged .and. allocated(merr) .and. niter == 2 .and. &
         & context%nobserved == 1 .and. context%ncheckpoints == 1)
      if (allocated(error)) return
      context = iteration_context_type(fail_checkpoint=.true.)
      call run_iterations(context, ctx, 10, 2.0_wp, 1, rmse, niter, converged, merr, &
         & checkpoint_before_convergence=.true.)
      call check(error, .not. converged .and. allocated(merr) .and. niter == 1)
      call ctx%delete()
   end subroutine test_iteration_controller

   !> Advance timing under the open node, balanced after a failed advance
   subroutine test_iteration_timer(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(iteration_context_type) :: context
      !> Quiet run context whose timer records the advance time
      type(moist_context_type), target :: ctx
      real(wp) :: rmse
      !> Attempted steps, open node and depth before the stage, node count and advance node
      integer :: niter, outer, depth, nnodes, id
      logical :: converged

      call new_context(ctx, nthreads=0, verbosity=0)
      call ctx%timer%start("test solve")
      outer = ctx%timer%current()
      depth = ctx%timer%current_depth()
      context = iteration_context_type(fail_at=3)
      call run_iterations(context, ctx, 10, 0.1_wp, 0, rmse, niter, converged, merr)
      call check(error, allocated(merr) .and. .not. converged .and. niter == 3 .and. context%nobserved == 2, &
         & "The third advance must fail after two observed steps")
      if (allocated(error)) return
      call check(error, context%min_elapsed >= 0.0_wp, "Observed advance times must be nonnegative")
      if (allocated(error)) return
      call check(error, ctx%timer%current_depth() == depth .and. ctx%timer%current() == outer, &
         & "A failed advance must leave the timer stack unchanged")
      if (allocated(error)) return
      nnodes = ctx%timer%num_nodes()
      id = ctx%timer%resolve("Iterations", outer)
      call check(error, ctx%timer%num_nodes() == nnodes, "The advance node must sit under the open node")
      if (allocated(error)) return
      call check(error, .not. ieee_is_nan(ctx%timer%node_time(id)) .and. ctx%timer%last(id) >= 0.0_wp, &
         & "A failed advance must not leave its node running")
      if (allocated(error)) return
      call ctx%timer%stop("test solve")
      call ctx%delete()
   end subroutine test_iteration_timer

   !> Weighted coarse projection and reconstruction agree with an independent dense solve
   subroutine test_deflation_algebra(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_solver_deflation_type) :: deflation
      real(wp) :: basis(3, 2), weighted(3, 2), rhs(3), coordinates(2, 1), expanded(3), gram(2, 2)
      real(wp) :: determinant, expected(2)

      basis(:, 1) = [1.0_wp, 2.0_wp, 0.0_wp]
      basis(:, 2) = [0.0_wp, 1.0_wp, 0.0_wp]
      weighted = basis*spread([2.0_wp, 3.0_wp, 4.0_wp], 2, 2)
      rhs = [3.0_wp, -1.0_wp, 9.0_wp]
      gram = matmul(transpose(weighted), basis)
      determinant = gram(1, 1)*gram(2, 2) - gram(1, 2)*gram(2, 1)
      expected = matmul(transpose(weighted), rhs)
      expected = [gram(2, 2)*expected(1) - gram(1, 2)*expected(2), &
         & gram(1, 1)*expected(2) - gram(2, 1)*expected(1)]/determinant
      call deflation%init_basis(basis, weighted, 1)
      deflation%coarse = gram
      call deflation%factorize()
      call check(error, deflation%active .and. deflation%ncore == 2)
      if (allocated(error)) return
      call deflation%coordinates(rhs, coordinates)
      call check(error, maxval(abs(coordinates(:, 1) - expected)) < 1.0e-14_wp)
      if (allocated(error)) return
      expanded = 0.0_wp
      call deflation%expand(coordinates, 0.5_wp, expanded)
      call check(error, maxval(abs(expanded - 0.5_wp*matmul(basis, expected))) < 1.0e-14_wp)
      if (allocated(error)) return
      deflation%coarse = 0.0_wp
      call deflation%factorize()
      call check(error, .not. deflation%active)
      call deflation%release()
   end subroutine test_deflation_algebra

   !> Concurrent solves, each with its own run context, keep their state independent
   subroutine test_concurrent_solves(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=6), parameter :: names(3) = [character(len=6) :: "PICARD", "DIIS", "JFNK"]
      logical :: success(3, 2)
      integer :: family, problem

      !$omp parallel do collapse(2) num_threads(2) schedule(static)
      do problem = 1, 2
         do family = 1, 3
            call solve_one(names(family), 0.6_wp + 0.4_wp*problem, success(family, problem))
         end do
      end do
      !$omp end parallel do
      call check(error, all(success), "Independent solver contexts must converge to their own roots")

   contains

      recursive subroutine solve_one(name, root, success)
         character(len=*), intent(in) :: name
         real(wp), intent(in) :: root
         logical, intent(out) :: success
         class(moist_model_moz_solver_type), allocatable :: solver
         type(moist_model_moz_solver_parameters_type) :: p
         type(quadratic_evaluator_type) :: evaluator
         type(quadratic_linearization_type), target :: state
         type(mctc_error), allocatable :: merr
         real(wp) :: gamma(1, 1), updated(1, 1), rmse
         logical :: feasible
         integer :: iter
         !> Run context of this solve; contexts are not shared between threads
         type(moist_context_type), target :: ctx

         success = .false.
         call new_context(ctx, nthreads=1, verbosity=0)
         solve: block
            p%diis%depth = 2
            p%name = name
            call new_solver(solver, ctx, merr, p)
            if (allocated(merr)) exit solve
            evaluator%state => state
            evaluator%root = root
            gamma = 0.5_wp
            call solver%begin_stage()
            do iter = 1, 400
               call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
               if (.not. feasible .or. allocated(merr)) exit solve
               if (rmse < 1.0e-11_wp) exit
               gamma = updated
            end do
            success = abs(gamma(1, 1) - root) < 1.0e-10_wp
         end block solve
         call ctx%delete()
      end subroutine solve_one

   end subroutine test_concurrent_solves

   !> Linear callbacks can report errors without also changing their feasibility flag
   !>
   !> @param[in] self Dense problem and controlled failure modes
   !> @param[in] gamma Current field (2, 1)
   !> @param[out] residual Dense residual (2, 1)
   !> @param[out] feasible Deliberately true even on callback errors
   !> @param[out] error Deliberate evaluation failure
   subroutine linear_residual(self, gamma, residual, feasible, error)
      !> Dense problem and failure modes
      class(linear_evaluator_type), intent(in) :: self
      !> Current field
      real(wp), intent(in) :: gamma(:, :)
      !> Dense residual
      real(wp), intent(out) :: residual(:, :)
      !> Deliberately true on the error path
      logical, intent(out) :: feasible
      !> Deliberate evaluation failure
      type(mctc_error), allocatable, intent(out) :: error

      feasible = .true.
      if (self%residual_mode == 1 .or. maxval(abs(gamma)) > self%trial_limit) then
         call fatal_error(error, "Deliberate residual callback error")
         return
      end if
      residual(:, 1) = matmul(self%matrix, gamma(:, 1) - self%target)
      if (self%residual_mode == 2) residual = ieee_value(0.0_wp, ieee_quiet_nan)
   end subroutine linear_residual

   !> The dense evaluator supplies exact products
   !>
   !> @param[in] self Dense evaluator
   logical function linear_has_jvp(self) result(available)
      !> Dense evaluator
      class(linear_evaluator_type), intent(in) :: self
      available = .true.
   end function linear_has_jvp

   !> Accept the constant dense Jacobian without storing a point
   !>
   !> @param[in] self Dense evaluator
   !> @param[in] gamma Linearization point (2, 1)
   !> @param[out] feasible Always true
   !> @param[out] error Unallocated
   subroutine linear_linearize(self, gamma, feasible, error)
      !> Dense evaluator
      class(linear_evaluator_type), intent(in) :: self
      !> Linearization point
      real(wp), intent(in) :: gamma(:, :)
      !> Constant Jacobian is available
      logical, intent(out) :: feasible
      !> Unallocated failure information
      type(mctc_error), allocatable, intent(out) :: error
      feasible = .true.
      if (associated(self%linearizations)) self%linearizations = self%linearizations + 1
   end subroutine linear_linearize

   !> Apply the dense Jacobian or a controlled bad product
   !>
   !> @param[in] self Dense evaluator and failure mode
   !> @param[in] direction Flat vector (2)
   !> @param[out] jv Dense product (2)
   !> @param[out] feasible Deliberately true on callback errors
   !> @param[out] error Deliberate product failure
   subroutine linear_jvp(self, direction, jv, feasible, error)
      !> Dense evaluator and failure mode
      class(linear_evaluator_type), intent(in) :: self
      !> Flat direction
      real(wp), intent(in) :: direction(:)
      !> Dense product
      real(wp), intent(out) :: jv(:)
      !> Deliberately true on the error path
      logical, intent(out) :: feasible
      !> Deliberate product failure
      type(mctc_error), allocatable, intent(out) :: error

      feasible = .true.
      if (self%jvp_mode == 1) then
         call fatal_error(error, "Deliberate Jacobian callback error")
         return
      end if
      jv = matmul(self%matrix, direction)
      if (self%jvp_mode == 2) jv = ieee_value(0.0_wp, ieee_quiet_nan)
      if (self%jvp_mode == 3) jv = 0.0_wp
   end subroutine linear_jvp

   !> Apply the exact transpose or the same controlled failures as the primal.
   !> @param[in] self Frozen matrix and failure mode.
   !> @param[in] direction Euclidean direction.
   !> @param[out] jtv Transpose product.
   !> @param[out] feasible Product usability.
   !> @param[out] error Deliberate callback failure.
   subroutine linear_adjoint(self, direction, jtv, feasible, error)
      !> Frozen matrix and failure mode.
      class(linear_evaluator_type), intent(in) :: self
      !> Flat direction and transpose product.
      real(wp), intent(in) :: direction(:)
      real(wp), intent(out) :: jtv(:)
      !> Deliberately true even when returning an error.
      logical, intent(out) :: feasible
      !> Controlled callback error.
      type(mctc_error), allocatable, intent(out) :: error

      feasible = .true.
      if (self%jvp_mode == 1) then
         call fatal_error(error, "Deliberate transpose callback error")
         return
      end if
      jtv = matmul(transpose(self%matrix), direction)
      if (self%jvp_mode == 2) jtv = ieee_value(0.0_wp, ieee_quiet_nan)
      if (self%jvp_mode == 3) jtv = 0.0_wp
   end subroutine linear_adjoint

   !> Fill the coarse matrix using the driver's profiles and the dense Jacobian
   !>
   !> @param[in] self Dense test evaluator
   !> @param[inout] deflation Driver basis and evaluator-owned coarse factorization
   !> @param[out] error Missing driver basis
   subroutine linear_deflation_setup(self, deflation, error)
      !> Dense test evaluator
      class(linear_evaluator_type), intent(in) :: self
      !> Coarse space with driver-supplied profiles
      type(moist_model_moz_solver_deflation_type), intent(inout) :: deflation
      !> Missing driver basis
      type(mctc_error), allocatable, intent(out) :: error

      if (.not. deflation%has_basis(2, 1)) then
         call fatal_error(error, "Deliberate missing coarse basis")
         return
      end if
      deflation%coarse = matmul(transpose(deflation%wprofile), matmul(self%matrix, deflation%profile))
      call deflation%factorize()
   end subroutine linear_deflation_setup

   !> Enable or disable the borrowed test operator
   !>
   !> @param[in] self Operator settings
   !> @param[in] rmse Previous residual RMS (unused)
   logical function matrix_ready(self, rmse) result(usable)
      !> Operator settings
      class(matrix_preconditioner_type), intent(in) :: self
      !> Previous residual RMS
      real(wp), intent(in) :: rmse
      usable = self%enabled
   end function matrix_ready

   !> Apply a fixed dense inverse operator and record its use
   !>
   !> @param[inout] self Operator, counters and optional forcing-term observation
   !> @param[inout] vflat Flat vector (2)
   subroutine matrix_apply(self, vflat)
      !> Operator and counters
      class(matrix_preconditioner_type), intent(inout) :: self
      !> Flat vector to replace
      real(wp), contiguous, target, intent(inout) :: vflat(:)

      self%napplied = self%napplied + 1
      if (associated(self%forcing)) self%applied_eta = self%forcing
      select case (self%output_mode)
      case (1)
         vflat = ieee_value(0.0_wp, ieee_quiet_nan)
      case (2)
         vflat = 0.0_wp
      case default
         vflat = self%scale*matmul(self%matrix, vflat)
      end select
   end subroutine matrix_apply

   !> Short status for the test printer
   !>
   !> @param[in] self Dense operator (unused)
   !> @param[in] rmse Previous residual RMS (unused)
   function matrix_describe(self, rmse) result(text)
      !> Dense operator
      class(matrix_preconditioner_type), intent(in) :: self
      !> Previous residual RMS
      real(wp), intent(in) :: rmse
      !> Status text
      character(len=:), allocatable :: text
      text = "test matrix"
   end function matrix_describe

   subroutine fixed_residual(self, gamma, residual, feasible, error)
      class(fixed_evaluator_type), intent(in) :: self
      real(wp), intent(in) :: gamma(:, :)
      real(wp), intent(out) :: residual(:, :)
      logical, intent(out) :: feasible
      type(mctc_error), allocatable, intent(out) :: error
      residual = self%value
      feasible = .true.
   end subroutine fixed_residual

   subroutine coupled_residual(self, gamma, residual, feasible, error)
      class(coupled_evaluator_type), intent(in) :: self
      real(wp), intent(in) :: gamma(:, :)
      real(wp), intent(out) :: residual(:, :)
      logical, intent(out) :: feasible
      type(mctc_error), allocatable, intent(out) :: error
      real(wp), parameter :: matrix(2, 2) = reshape([0.3_wp, 0.1_wp, 0.05_wp, 0.6_wp], [2, 2])
      residual(:, 1) = matmul(matrix, gamma(:, 1) - [1.0_wp, -2.0_wp])
      feasible = .true.
   end subroutine coupled_residual

   subroutine quadratic_residual(self, gamma, residual, feasible, error)
      class(quadratic_evaluator_type), intent(in) :: self
      real(wp), intent(in) :: gamma(:, :)
      real(wp), intent(out) :: residual(:, :)
      logical, intent(out) :: feasible
      type(mctc_error), allocatable, intent(out) :: error
      feasible = maxval(abs(gamma)) <= self%max_argument
      residual = gamma**2 - self%root**2
   end subroutine quadratic_residual

   logical function quadratic_has_jvp(self) result(available)
      class(quadratic_evaluator_type), intent(in) :: self
      available = .true.
   end function quadratic_has_jvp

   subroutine quadratic_linearize(self, gamma, feasible, error)
      class(quadratic_evaluator_type), intent(in) :: self
      real(wp), intent(in) :: gamma(:, :)
      logical, intent(out) :: feasible
      type(mctc_error), allocatable, intent(out) :: error
      feasible = self%accept_linearization
      if (.not. feasible) then
         call fatal_error(error, "Deliberately unavailable exact linearization")
         return
      end if
      self%state%point = gamma(1, 1)
   end subroutine quadratic_linearize

   subroutine quadratic_jvp(self, direction, jv, feasible, error)
      class(quadratic_evaluator_type), intent(in) :: self
      real(wp), intent(in) :: direction(:)
      real(wp), intent(out) :: jv(:)
      logical, intent(out) :: feasible
      type(mctc_error), allocatable, intent(out) :: error
      jv = 2.0_wp*self%state%point*direction
      feasible = .true.
   end subroutine quadratic_jvp

   logical function diagonal_ready(self, rmse) result(usable)
      class(diagonal_preconditioner_type), intent(in) :: self
      real(wp), intent(in) :: rmse
      usable = .true.
   end function diagonal_ready

   subroutine diagonal_apply(self, vflat)
      class(diagonal_preconditioner_type), intent(inout) :: self
      real(wp), contiguous, target, intent(inout) :: vflat(:)
      self%napplied = self%napplied + 1
      vflat = 0.5_wp*vflat
   end subroutine diagonal_apply

   function diagonal_describe(self, rmse) result(text)
      class(diagonal_preconditioner_type), intent(in) :: self
      real(wp), intent(in) :: rmse
      character(len=:), allocatable :: text
      text = "test diagonal"
   end function diagonal_describe

   subroutine context_advance(self, iter, rmse, feasible, error)
      class(iteration_context_type), intent(inout) :: self
      integer, intent(in) :: iter
      real(wp), intent(out) :: rmse
      logical, intent(out) :: feasible
      type(mctc_error), allocatable, intent(out) :: error
      self%nsteps = self%nsteps + 1
      rmse = 1.0_wp/real(iter, wp)
      feasible = iter /= self%fail_at
      if (.not. feasible) call fatal_error(error, "Deliberate iteration failure")
   end subroutine context_advance

   subroutine context_observe(self, iter, rmse, elapsed)
      class(iteration_context_type), intent(inout) :: self
      integer, intent(in) :: iter
      real(wp), intent(in) :: rmse, elapsed
      self%nobserved = self%nobserved + 1
      self%min_elapsed = min(self%min_elapsed, elapsed)
   end subroutine context_observe

   subroutine context_checkpoint(self, error)
      class(iteration_context_type), intent(inout) :: self
      type(mctc_error), allocatable, intent(out) :: error
      self%ncheckpoints = self%ncheckpoints + 1
      if (self%fail_checkpoint) call fatal_error(error, "Deliberate checkpoint failure")
   end subroutine context_checkpoint

   subroutine test_mgs_streaming(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      !> Uneven block tail above the threading threshold.
      integer, parameter :: n = 65539
      !> Orthonormal basis, input vector, expected coefficients, and remainder.
      real(wp), allocatable :: v(:, :), w(:), original(:)
      real(wp) :: h(4), coefficients(4)
      !> Kernel selector, thread count, and original thread budget.
      integer :: mode, threads, saved, i, j, k
      saved = 1
!$    saved = omp_get_max_threads()
      allocate(v(n, 5), w(n), original(n))
      do j = 1, 5
         do i = 1, n
            v(i, j) = sin(0.013_wp*i*j)
         end do
         do k = 1, j - 1
            v(:, j) = v(:, j) - dot_product(v(:, j), v(:, k))*v(:, k)
         end do
         v(:, j) = v(:, j)/sqrt(dot_product(v(:, j), v(:, j)))
      end do
      coefficients = [0.2_wp, -0.7_wp, 0.5_wp, 0.1_wp]
      original = matmul(v(:, 1:4), coefficients) + 0.3_wp*v(:, 5)
      do mode = 0, 1
         do threads = 1, 4, 3
!$          call omp_set_num_threads(threads)
            w = original
            call moist_orthogonalize_mgs(v(:, 1:4), w, h, fused=mode == 1)
            call check(error, maxval(abs(h - coefficients)) < 1.0e-12_wp)
            if (allocated(error)) exit
            call check(error, maxval(abs(w - 0.3_wp*v(:, 5))) < 1.0e-13_wp)
            if (allocated(error)) exit
         end do
         if (allocated(error)) exit
      end do
      if (.not. allocated(error)) then
         w = original
         call moist_orthogonalize_mgs(v(:, 1:0), w, h(1:0), fused=.true.)
         call check(error, all(w == original))
      end if
!$    call omp_set_num_threads(saved)
   end subroutine test_mgs_streaming

   subroutine test_blas_row_blocks(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      !> Full matrix, source vectors, and independent expected products.
      real(wp) :: a(17, 5), x(17), h(5), y(17), projected(5), expected(5)
      real(wp) :: b(17, 3), small(5, 3), gram(5, 3), product(8, 3)
      !> Row and column indices, and block bounds.
      integer :: i, j, first, count
      do j = 1, 5
         do i = 1, 17
            a(i, j) = sin(real(i*j, wp))
         end do
         h(j) = 0.1_wp*j
      end do
      x = a(:, 2)
      do j = 1, 3
         b(:, j) = a(:, j + 1)
         small(:, j) = real(j, wp)*h
      end do
      do first = 1, 17, 8
         count = min(8, 18 - first)
         y = 0.3_wp
         call gemv_rows(a, first, count, h, y(first:first + count - 1), -0.7_wp, 0.5_wp, "N")
         call check(error, maxval(abs(y(first:first + count - 1) &
            & - (0.15_wp - 0.7_wp*matmul(a(first:first + count - 1, :), h)))) < 1.0e-14_wp)
         if (allocated(error)) return
         projected = h
         expected = 0.5_wp*h + 0.7_wp*matmul(transpose(a(first:first + count - 1, :)), &
            & x(first:first + count - 1))
         call gemv_rows(a, first, count, x(first:first + count - 1), projected, 0.7_wp, 0.5_wp, "T")
         call check(error, maxval(abs(projected - expected)) < 1.0e-14_wp)
         if (allocated(error)) return
         gram = 0.2_wp
         call gemm_rows_tn(a, b, first, count, gram, 0.7_wp, 0.5_wp)
         call check(error, maxval(abs(gram - (0.1_wp + 0.7_wp*matmul( &
            & transpose(a(first:first + count - 1, :)), b(first:first + count - 1, :))))) < 1.0e-14_wp)
         if (allocated(error)) return
         product = 0.2_wp
         call gemm_rows_nn(a, first, count, small, product, -0.7_wp, 0.5_wp)
         call check(error, maxval(abs(product(1:count, :) - (0.1_wp - 0.7_wp*matmul( &
            & a(first:first + count - 1, :), small)))) < 1.0e-14_wp)
         if (allocated(error)) return
      end do

      ! An empty inner dimension still applies C := beta*C, as dgemm does;
      ! beta = 0 overwrites C without reading it
      gram = ieee_value(0.0_wp, ieee_quiet_nan)
      call gemm_rows_tn(a, b, 18, 0, gram, 0.7_wp, 0.0_wp)
      call check(error, all(gram == 0.0_wp), "An empty row block with beta = 0 must clear the Gram result")
      if (allocated(error)) return
      gram = 0.2_wp
      call gemm_rows_tn(a, b, 1, 0, gram, 0.7_wp, 0.5_wp)
      call check(error, maxval(abs(gram - 0.1_wp)) < 1.0e-15_wp, "An empty row block must scale the Gram result by beta")
      if (allocated(error)) return
      product = 0.2_wp
      call gemm_rows_nn(a(:, 1:0), 1, 8, small(1:0, :), product, -0.7_wp, 0.5_wp)
      call check(error, maxval(abs(product - 0.1_wp)) < 1.0e-15_wp, "An empty inner dimension must scale the product by beta")
   end subroutine test_blas_row_blocks

   subroutine test_solver_lifecycle(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      !> Concrete solvers inspected through the abstract contract.
      type(moist_model_moz_solver_jfnk_type), target :: gmres
      type(moist_model_moz_solver_picard_type), target :: picard
      type(moist_model_moz_solver_diis_type), target :: mdiis
      !> Borrowed polymorphic solver and printer targets.
      class(moist_model_moz_solver_type), pointer :: solver
      type(moist_model_moz_solver_printer_type), target :: printer, replacement
      !> Allocation-free diagnostic snapshot.
      type(moist_model_moz_solver_statistics_type) :: stats
      !> Run context borrowed by the three solvers
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      gmres = new_jfnk_solver(ctx)
      mdiis = new_diis_solver(ctx)
      picard = new_picard_solver(ctx)

      solver => gmres
      gmres%eta_current = 0.3_wp
      gmres%gmres_iters_last = 7
      call solver%set_iteration_printer(printer)
      call check(error, associated(gmres%iter_printer, printer))
      if (allocated(error)) return
      call solver%set_iteration_printer(replacement)
      call check(error, associated(gmres%iter_printer, replacement))
      if (allocated(error)) return
      call solver%begin_stage()
      call check(error, gmres%eta_current == -1.0_wp)
      if (allocated(error)) return
      stats = solver%get_statistics()
      call check(error, stats%ngmres == 7)
      if (allocated(error)) return
      call solver%set_iteration_printer()
      call check(error, .not. associated(gmres%iter_printer))
      if (allocated(error)) return

      solver => mdiis
      mdiis%nstored = 3
      call solver%set_iteration_printer(printer)
      call check(error, associated(mdiis%iter_printer, printer))
      if (allocated(error)) return
      call solver%set_iteration_printer()
      call check(error, .not. associated(mdiis%iter_printer))
      if (allocated(error)) return
      call solver%begin_stage()
      stats = solver%get_statistics()
      call check(error, stats%ndiis == 0 .and. stats%ngmres == 0)
      if (allocated(error)) return

      solver => picard
      call solver%begin_stage()
      call solver%set_iteration_printer(printer)
      call solver%set_iteration_printer()
      stats = solver%get_statistics()
      call check(error, stats%ndiis == 0 .and. stats%ngmres == 0)
      if (allocated(error)) return
      call ctx%delete()
   end subroutine test_solver_lifecycle

   subroutine test_reduction_threads(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      !> Independent solver instances with identical initial state.
      type(moist_model_moz_solver_jfnk_type) :: serial, parallel, pinned
      !> Run contexts of the serial, parallel and pinned solvers, one solver each
      type(moist_context_type), target :: ctx1, ctx4, ctx3
      !> Controls shared by all three solvers
      type(jfnk_parameters) :: controls
      !> Residual callback errors.
      type(mctc_error), allocatable :: err1, err2
      !> Initial and solved fields with incomplete reduction blocks.
      real(wp), allocatable :: gamma(:, :), one(:, :), four(:, :), three(:, :)
      !> Residual norms.
      real(wp) :: rmse1, rmse2
      !> Feasibility indicators.
      logical :: ok1, ok2
      !> Original OpenMP team size.
      integer :: saved_threads

      saved_threads = 1
      !$ saved_threads = omp_get_max_threads()
      allocate (gamma(32769, 3), one(32769, 3), four(32769, 3), three(32769, 3))
      gamma = 0.0_wp
      controls = jfnk_parameters(restart=4, max_krylov=8, eta0=1.0e-8_wp, use_exact_jvp=.false., fd_eps=1.0e-3_wp)
      call new_context(ctx1, nthreads=1, verbosity=0)
      call new_context(ctx4, nthreads=4, verbosity=0)
      call new_context(ctx3, nthreads=3, verbosity=0)
      serial = new_jfnk_solver(ctx1, controls)
      parallel = new_jfnk_solver(ctx4, controls)
      ! The context count sizes the teams regardless of the host OpenMP setting
      pinned = new_jfnk_solver(ctx3, controls)
      !$ call omp_set_num_threads(1)
      call serial%outer_step(diagonal_evaluator_type(), gamma, one, rmse1, ok1, err1)
      call pinned%outer_step(diagonal_evaluator_type(), gamma, three, rmse2, ok2, err2)
      !$ call omp_set_num_threads(saved_threads)
      call check(error, pinned%team_size() == 3 .and. .not. allocated(err2) .and. ok2 .and. all(three == one), &
         & "An explicit thread count must reproduce the serial result")
      call parallel%outer_step(diagonal_evaluator_type(), gamma, four, rmse2, ok2, err2)
      if (allocated(error)) return
      call check(error, parallel%team_size() == 4 .and. .not. allocated(err1) .and. .not. allocated(err2) &
         & .and. ok1 .and. ok2)
      if (allocated(error)) return
      call check(error, all(one == four) .and. rmse1 == rmse2)
      if (allocated(error)) return
      call check(error, serial%gmres_iters_last == parallel%gmres_iters_last)
      call ctx3%delete()
      call ctx4%delete()
      call ctx1%delete()
   end subroutine test_reduction_threads

   subroutine test_parallel_step(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      !> Solver with a small Krylov basis and no preconditioner.
      type(moist_model_moz_solver_jfnk_type) :: solver
      !> Residual callback failure information.
      type(mctc_error), allocatable :: merr
      !> Input and result fields; size exceeds the parallel-vector threshold.
      real(wp), allocatable :: gamma(:, :), result(:, :)
      !> Small field used to exercise workspace shape changes.
      real(wp) :: small(7, 2)
      !> Residual norm and analytic error
      real(wp) :: rmse, max_error
      !> Residual feasibility flag.
      logical :: feasible
      !> Spatial and pair indices.
      integer :: i, j
      !> Quiet run context borrowed by the solver
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=4, max_krylov=8, eta0=1.0e-8_wp, use_exact_jvp=.false., &
         & fd_eps=1.0e-3_wp))
      allocate (gamma(32769, 3), result(32769, 3))
      gamma = 0.0_wp
      small = 0.0_wp
      call solver%ensure_workspace(gamma)
      call check(error, all(solver%workspace%ws_v == 0.0_wp) .and. .not. allocated(solver%workspace%ws_f_pert))
      if (allocated(error)) return
      call check(error, solver%workspace%ws_nflat == size(gamma))
      if (allocated(error)) return
      solver%workspace%ws_v(1, 1) = 7.0_wp
      call solver%ensure_workspace(gamma)
      call check(error, solver%workspace%ws_v(1, 1) == 7.0_wp)
      if (allocated(error)) return
      call solver%ensure_workspace(small)
      call check(error, solver%workspace%ws_npts == 7 .and. solver%workspace%ws_npair == 2)
      if (allocated(error)) return
      call check(error, all(solver%workspace%ws_v == 0.0_wp))
      if (allocated(error)) return

      call solver%outer_step(diagonal_evaluator_type(), gamma, result, rmse, feasible, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, feasible)
      if (allocated(error)) return
      max_error = 0.0_wp
      do j = 1, size(result, 2)
         do i = 1, size(result, 1)
            max_error = max(max_error, abs(result(i, j) - target_value(i, j)))
         end do
      end do
      call check(error, max_error < 1.0e-7_wp)
      if (allocated(error)) return
      call check(error, solver%gmres_iters_last >= 2 .and. solver%gmres_iters_last <= 4)
      call solver%release_workspace()
      call ctx%delete()
   end subroutine test_parallel_step

   subroutine test_parallel_restart(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      !> Solver whose basis is smaller than the number of distinct eigenvalues.
      type(moist_model_moz_solver_jfnk_type) :: solver
      !> Residual callback failure information.
      type(mctc_error), allocatable :: merr
      !> Input and result fields above the parallel-vector threshold.
      real(wp), allocatable :: gamma(:, :), result(:, :)
      !> Initial residual norm and maximum analytic error.
      real(wp) :: rmse, max_error
      !> Residual feasibility indicator.
      logical :: feasible
      !> Spatial and pair indices.
      integer :: i, j
      !> Quiet run context borrowed by the solver
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=20, eta0=1.0e-8_wp, use_exact_jvp=.false., &
         & fd_eps=1.0e-3_wp))
      allocate (gamma(32769, 3), result(32769, 3))
      gamma = 0.0_wp
      call solver%outer_step(diagonal_evaluator_type(), gamma, result, rmse, feasible, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, feasible .and. solver%gmres_iters_last > solver%params%restart)
      if (allocated(error)) return
      max_error = 0.0_wp
      do j = 1, size(result, 2)
         do i = 1, size(result, 1)
            max_error = max(max_error, abs(result(i, j) - target_value(i, j)))
         end do
      end do
      call check(error, max_error < 1.0e-7_wp)
      call ctx%delete()
   end subroutine test_parallel_restart

   subroutine test_parallel_cgs2(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      !> Solver whose basis is smaller than the number of distinct eigenvalues.
      type(moist_model_moz_solver_jfnk_type) :: solver
      !> Residual callback failure information.
      type(mctc_error), allocatable :: merr
      !> Input and result fields above the parallel-vector threshold.
      real(wp), allocatable :: gamma(:, :), result(:, :)
      !> Initial residual norm and maximum analytic error.
      real(wp) :: rmse, max_error
      !> Residual feasibility indicator.
      logical :: feasible
      !> Spatial and pair indices.
      integer :: i, j
      !> Quiet run context borrowed by the solver
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=20, eta0=1.0e-8_wp, use_exact_jvp=.false., &
         & fd_eps=1.0e-3_wp))
      solver%params%cgs2 = .true.
      allocate (gamma(32769, 3), result(32769, 3))
      gamma = 0.0_wp
      call solver%outer_step(diagonal_evaluator_type(), gamma, result, rmse, feasible, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, feasible .and. solver%gmres_iters_last > solver%params%restart)
      if (allocated(error)) return
      max_error = 0.0_wp
      do j = 1, size(result, 2)
         do i = 1, size(result, 1)
            max_error = max(max_error, abs(result(i, j) - target_value(i, j)))
         end do
      end do
      call check(error, max_error < 1.0e-7_wp)
      call ctx%delete()
   end subroutine test_parallel_cgs2

   subroutine test_happy_breakdown(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      !> Solver using exact finite differences for the identity operator.
      type(moist_model_moz_solver_jfnk_type) :: solver
      !> Residual callback failure information.
      type(mctc_error), allocatable :: merr
      !> Zero state and residual, unit right-hand side, and computed solution.
      real(wp), allocatable :: gamma(:, :), residual(:, :), rhs(:), delta(:)
      !> Feasibility and preconditioner activation flags.
      logical :: feasible, precond_valid
      !> Number of inner iterations used.
      integer :: iterations
      !> Quiet run context borrowed by the solver
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=4, use_exact_jvp=.false., fd_eps=1.0_wp))
      allocate (gamma(32769, 2), residual(32769, 2), rhs(65538), delta(65538))
      gamma = 0.0_wp
      residual = 0.0_wp
      rhs = 0.0_wp
      rhs(1) = 1.0_wp
      call solver%ensure_workspace(gamma)
      solver%workspace%ws_v(:, 2) = 3.0_wp
      precond_valid = .false.
      call solver%linear_solve(identity_evaluator_type(), gamma, residual, rhs, delta, &
                               1.0e-8_wp, precond_valid, feasible, merr, iterations)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, feasible .and. iterations == 1)
      if (allocated(error)) return
      call check(error, all(delta == rhs))
      if (allocated(error)) return
      call check(error, all(solver%workspace%ws_v(:, 2) == 0.0_wp))
      call ctx%delete()
   end subroutine test_happy_breakdown

   subroutine test_infeasible_jvp(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      !> Solver for a deliberately failing residual callback.
      type(moist_model_moz_solver_jfnk_type) :: solver
      !> Expected callback error.
      type(mctc_error), allocatable :: merr
      !> State, current residual, right-hand side, and output direction.
      real(wp), allocatable :: gamma(:, :), residual(:, :), rhs(:), delta(:)
      !> Feasibility and preconditioner activation flags.
      logical :: feasible, precond_valid
      !> Number of attempted inner iterations.
      integer :: iterations
      !> Quiet run context borrowed by the solver
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=2, max_krylov=4, use_exact_jvp=.false., fd_eps=1.0_wp))
      allocate (gamma(32769, 2), residual(32769, 2), rhs(65538), delta(65538))
      gamma = 0.0_wp
      residual = 0.0_wp
      rhs = 0.0_wp
      rhs(1) = 1.0_wp
      delta = 7.0_wp
      call solver%ensure_workspace(gamma)
      precond_valid = .false.
      call solver%linear_solve(infeasible_evaluator_type(), gamma, residual, rhs, delta, &
                               1.0e-8_wp, precond_valid, feasible, merr, iterations)
      call check(error, .not. feasible .and. allocated(merr))
      if (allocated(error)) return
      call check(error, merr%message == "deliberately infeasible JVP")
      if (allocated(error)) return
      call check(error, iterations == 1 .and. all(delta == 0.0_wp))
      call ctx%delete()
   end subroutine test_infeasible_jvp

   subroutine identity_residual(self, gamma, residual, feasible, error)
      !> Test evaluator context.
      class(identity_evaluator_type), intent(in) :: self
      !> Current field.
      real(wp), intent(in) :: gamma(:, :)
      !> Identity residual.
      real(wp), intent(out) :: residual(:, :)
      !> Feasibility indicator.
      logical, intent(out) :: feasible
      !> Error handle (unused for the identity residual).
      type(mctc_error), allocatable, intent(out) :: error
      residual = gamma
      feasible = .true.
   end subroutine identity_residual

   subroutine infeasible_residual(self, gamma, residual, feasible, error)
      !> Test evaluator context.
      class(infeasible_evaluator_type), intent(in) :: self
      !> Current field (unused).
      real(wp), intent(in) :: gamma(:, :)
      !> Invalid output residual.
      real(wp), intent(out) :: residual(:, :)
      !> Feasibility indicator.
      logical, intent(out) :: feasible
      !> Expected callback failure information.
      type(mctc_error), allocatable, intent(out) :: error
      residual = 0.0_wp
      feasible = .false.
      call fatal_error(error, "deliberately infeasible JVP")
   end subroutine infeasible_residual

   pure function target_value(i, j) result(value)
      !> Spatial and pair indices.
      integer, intent(in) :: i, j
      !> Root component, distinct across both field dimensions.
      real(wp) :: value
      value = 0.125_wp*real(j, wp) + 0.01_wp*real(mod(i, 13), wp)
   end function target_value

   subroutine diagonal_residual(self, gamma, residual, feasible, error)
      !> Test evaluator context.
      class(diagonal_evaluator_type), intent(in) :: self
      !> Current field.
      real(wp), intent(in) :: gamma(:, :)
      !> Residual with the same shape.
      real(wp), intent(out) :: residual(:, :)
      !> Feasibility indicator.
      logical, intent(out) :: feasible
      !> Error handle (unused for this analytic residual).
      type(mctc_error), allocatable, intent(out) :: error
      !> Spatial and pair indices.
      integer :: i, j
      do j = 1, size(gamma, 2)
         do i = 1, size(gamma, 1)
            residual(i, j) = (1.0_wp + 0.25_wp*real(j, wp))*(gamma(i, j) - target_value(i, j))
         end do
      end do
      feasible = .true.
   end subroutine diagonal_residual

   subroutine test_recycle_small(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      call check_recycle(error, 200, 1)
      if (allocated(error)) return
      call check_recycle_budget(error)
   end subroutine test_recycle_small

   !> A first solve capped below the restart length still harvests a recycled subspace
   !>
   !> @param[out] error Test failure information
   subroutine check_recycle_budget(error)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Budget-limited recycling solver and the soft test problem
      type(moist_model_moz_solver_jfnk_type) :: solver
      type(soft_evaluator_type) :: evaluator
      type(soft_linearization_type), target :: state
      !> Residual callback failure information
      type(mctc_error), allocatable :: merr
      !> Iterate, update and residual RMS
      real(wp) :: gamma(200, 1), gamma_new(200, 1), rmse
      logical :: feasible
      !> Quiet run context borrowed by the solver
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      evaluator%state => state
      allocate (state%point(size(gamma)))
      ! A tight forcing term makes the first solve use its whole 8-iteration
      ! budget, which is shorter than the 12-vector restart cycle
      solver = new_jfnk_solver(ctx, jfnk_parameters(restart=12, max_krylov=8, recycle=4, eta0=1.0e-8_wp))
      gamma = 0.0_wp
      call solver%outer_step(evaluator, gamma, gamma_new, rmse, feasible, merr)
      call check(error, feasible .and. .not. allocated(merr))
      if (allocated(error)) return
      call check(error, solver%gmres_iters_last == 8 .and. solver%nrecycled > 0, &
                 more="a budget-limited first solve must harvest a recycled subspace")
      call ctx%delete()
   end subroutine check_recycle_budget

   !> Reducing active recycling keeps shared storage and drops the oversized basis.
   !> @param[out] error Test failure information
   subroutine test_recycle_reconfigure(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(moist_model_moz_solver_jfnk_type), target :: solver
      type(soft_evaluator_type) :: evaluator
      type(soft_linearization_type), target :: state
      type(c_ptr) :: basis, recycled
      real(wp) :: gamma(200, 1), updated(200, 1), residual(200, 1), rmse, after
      type(moist_context_type), target :: ctx
      integer :: limit
      logical :: feasible

      call new_context(ctx, nthreads=0, verbosity=0)
      evaluator%state => state
      allocate (state%point(size(gamma)))
      do limit = 0, 1
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=12, max_krylov=8, recycle=4, eta0=1.0e-8_wp))
         gamma = 0.0_wp
         call solver%ensure_workspace(gamma)
         solver%workspace%shared = .true.
         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         call check(error, solver%nrecycled > 1, "The test must retain an active recycled basis")
         if (allocated(error)) return
         basis = c_loc(solver%workspace%ws_v)
         recycled = c_loc(solver%workspace%ws_u)
         gamma = updated

         ! Reconfigure within the active stage, without clearing it through begin_stage.
         solver%params%recycle = limit
         solver%params%max_krylov = 24
         call solver%ensure_workspace(gamma)
         call check(error, solver%nrecycled == 0, "An oversized recycled basis must be discarded")
         if (allocated(error)) return
         call check(error, c_associated(basis, c_loc(solver%workspace%ws_v)) .and. &
            & c_associated(recycled, c_loc(solver%workspace%ws_u)), &
            & "Reducing recycling must preserve shared allocations")
         if (allocated(error)) return

         call solver%outer_step(evaluator, gamma, updated, rmse, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         call evaluator%residual(updated, residual, feasible, merr)
         call check(error, feasible .and. .not. allocated(merr))
         if (allocated(error)) return
         after = sqrt(sum(residual*residual)/real(size(residual), wp))
         call check(error, after < rmse .and. solver%nrecycled <= limit, &
            & "The next update must reduce the residual within the new recycling limit")
         if (allocated(error)) return
         call check(error, c_associated(basis, c_loc(solver%workspace%ws_v)) .and. &
            & c_associated(recycled, c_loc(solver%workspace%ws_u)), &
            & "The next update must reuse shared allocations")
         if (allocated(error)) return
      end do
      call ctx%delete()
   end subroutine test_recycle_reconfigure

   subroutine test_recycle_parallel(error)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      call check_recycle(error, 35000, 2)
   end subroutine test_recycle_parallel

   subroutine check_recycle(error, npts, npair)
      !> Test failure information.
      type(error_type), allocatable, intent(out) :: error
      !> Field shape
      integer, intent(in) :: npts, npair
      !> Solvers with and without recycling
      type(moist_model_moz_solver_jfnk_type) :: solver
      type(soft_evaluator_type) :: evaluator
      type(soft_linearization_type), target :: state
      !> Residual callback failure information.
      type(mctc_error), allocatable :: merr
      !> Iterate and update
      real(wp), allocatable :: gamma(:, :), gamma_new(:, :)
      !> Residual RMS and largest error
      real(wp) :: rmse, max_error
      !> Krylov iterations and outer steps of both runs
      integer :: nkrylov(2), nouter(2), irun, it, i, j
      logical :: feasible
      !> Quiet run context borrowed by both solvers
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      evaluator%state => state
      allocate (state%point(npts*npair))
      do irun = 1, 2
         solver = new_jfnk_solver(ctx, jfnk_parameters(restart=10, max_krylov=10, recycle=merge(0, 4, irun == 1)))
         allocate (gamma(npts, npair), gamma_new(npts, npair))
         gamma = 0.0_wp
         nouter(irun) = 0
         do it = 1, 80
            call solver%outer_step(evaluator, gamma, gamma_new, rmse, feasible, merr)
            if (allocated(merr)) then
               call test_failed(error, merr%message)
               return
            end if
            if (rmse < 1.0e-11_wp) exit
            gamma = gamma_new
            nouter(irun) = it
         end do
         nkrylov(irun) = solver%nkrylov_total
         max_error = 0.0_wp
         do j = 1, npair
            do i = 1, npts
               max_error = max(max_error, abs(gamma(i, j) - soft_target(i + (j - 1)*npts)))
            end do
         end do
         ! The error along the smallest eigenvalue (1e-3) is bounded by
         ! ||R|| / lambda_min = rmse*sqrt(n)/1e-3, ~2e-6 for n = 70000.
         call check(error, rmse < 1.0e-11_wp .and. max_error < 1.0e-5_wp, &
                    more="Newton-GMRES did not reach the root")
         deallocate (gamma, gamma_new)
         if (allocated(error)) return
      end do
      call check(error, nkrylov(2) < nkrylov(1), &
                 more="recycling did not reduce the Krylov iterations")
      call ctx%delete()
   end subroutine check_recycle

   pure function soft_diag(i) result(d)
      !> Flat index
      integer, intent(in) :: i
      !> Diagonal entry
      real(wp) :: d
      real(wp), parameter :: small(4) = [1.0e-3_wp, 3.0e-3_wp, 1.0e-2_wp, 3.0e-2_wp]
      if (i <= 4) then
         d = small(i)
      else
         d = 1.0_wp + real(mod(i, 97), wp)/97.0_wp
      end if
   end function soft_diag

   pure function soft_target(i) result(t)
      !> Flat index
      integer, intent(in) :: i
      !> Root component
      real(wp) :: t
      t = 0.5_wp*sin(0.1_wp*real(i, wp))
   end function soft_target

   subroutine soft_residual(self, gamma, residual, feasible, error)
      !> Test evaluator context.
      class(soft_evaluator_type), intent(in) :: self
      !> Current field.
      real(wp), intent(in) :: gamma(:, :)
      !> Residual.
      real(wp), intent(out) :: residual(:, :)
      !> Feasibility indicator.
      logical, intent(out) :: feasible
      !> Error handle (unused).
      type(mctc_error), allocatable, intent(out) :: error
      !> Indices and offset from the root
      integer :: i, j, l
      real(wp) :: e
      do j = 1, size(gamma, 2)
         do i = 1, size(gamma, 1)
            l = i + (j - 1)*size(gamma, 1)
            e = gamma(i, j) - soft_target(l)
            residual(i, j) = soft_diag(l)*e + soft_cubic*e**3
         end do
      end do
      feasible = .true.
   end subroutine soft_residual

   logical function soft_has_jvp(self) result(available)
      !> Test evaluator context.
      class(soft_evaluator_type), intent(in) :: self
      available = .true.
   end function soft_has_jvp

   subroutine soft_linearize(self, gamma, feasible, error)
      !> Test evaluator context.
      class(soft_evaluator_type), intent(in) :: self
      !> Linearization point.
      real(wp), intent(in) :: gamma(:, :)
      !> Feasibility indicator.
      logical, intent(out) :: feasible
      !> Error handle (unused).
      type(mctc_error), allocatable, intent(out) :: error
      self%state%point(:) = reshape(gamma, [size(gamma)])
      feasible = .true.
   end subroutine soft_linearize

   subroutine soft_jvp(self, direction, jv, feasible, error)
      !> Test evaluator context.
      class(soft_evaluator_type), intent(in) :: self
      !> Direction.
      real(wp), intent(in) :: direction(:)
      !> Product.
      real(wp), intent(out) :: jv(:)
      !> Feasibility indicator.
      logical, intent(out) :: feasible
      !> Error handle (unused).
      type(mctc_error), allocatable, intent(out) :: error
      !> Flat index
      integer :: l
      do l = 1, size(direction)
         jv(l) = (soft_diag(l) + 3.0_wp*soft_cubic*(self%state%point(l) - soft_target(l))**2) &
                & *direction(l)
      end do
      feasible = .true.
   end subroutine soft_jvp

end module test_moz_solver
