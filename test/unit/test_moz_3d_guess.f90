!> 3D UV solves seeded from supplied guesses
!>
!> - A converged guess is the solution, without iterations
!> - A guess for part of the points starts from the selected cold field elsewhere
!> - Unusable guesses start cold; a failed warm attempt retries cold
!> - A solve that stops on a closure switch converges only on the requested closure
module test_moz_3d_guess
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_context, only: moist_context_type, new_context
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_model_moz_3d_driver, only: uv3d_problem_type, uv3d_evaluator_type, new_uv3d_problem
   use test_helpers, only: check_moist_error
   use test_moz_fixtures, only: make_uv3d_grid, make_uv3d_data
   implicit none(type, external)
   private
   public :: collect_moz_3d_guess

   !> Solvent site densities of the synthetic problems, 1/bohr^3
   real(wp), parameter :: rho(2) = [0.01_wp, 0.014_wp]
   !> Inverse thermal energy of the synthetic problems, 1/Hartree
   real(wp), parameter :: beta = 1.7_wp
   !> Cold-guess modes
   character(len=3), parameter :: modes(2) = ["cs0", "ts0"]

contains

   !> Register the guess tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_3d_guess(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("guess_full", check_full), &
         & new_unittest("guess_partial", check_partial), &
         & new_unittest("guess_rejected", check_rejected), &
         & new_unittest("guess_failed_warm", check_failed_warm), &
         & new_unittest("ramp_switch_end", ramp_switch_end)]
   end subroutine collect_moz_3d_guess

   !> A converged field as the guess is the solution, without iterations
   !>
   !> @param[out] error  test failure
   subroutine check_full(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type) :: cold, warm
      type(model_moz_3d_parameters_type) :: params
      real(wp), allocatable :: sr(:, :), lr(:, :)
      complex(wp), allocatable :: uk(:, :)
      real(wp) :: rmse
      integer :: niter
      logical :: converged

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_uv3d_grid(1, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_uv3d_data(grid, response, sr, lr, uk)
      call new_uv3d_problem(cold, ctx, grid, response, sr, lr, uk, rho, beta, params, err)
      if (.not. allocated(err)) call cold%solve(1.0e-12_wp, 128, converged, rmse, niter, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter > 0, more="cold reference")
      if (allocated(error)) return
      call new_uv3d_problem(warm, ctx, grid, response, sr, lr, uk, rho, beta, params, err)
      if (.not. allocated(err)) call warm%solve(1.0e-12_wp, 128, converged, rmse, niter, err, guess=cold%tr, &
         & matched=spread(.true., 1, grid%ngrid))
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == 0, more="converged guess was not accepted")
      if (allocated(error)) return
      call check(error, all(warm%tr == cold%tr), more="accepted guess changed")
      call ctx%delete()
   end subroutine check_full

   !> Unmatched points start from the selected cold field, in both cold modes
   !>
   !> - A masked guess solves exactly like the full field built by hand from the cold
   !>   field and the guess
   !> - A guess matching no point is a cold start
   !>
   !> @param[out] error  test failure
   subroutine check_partial(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type) :: cold, masked, mixed, none
      type(model_moz_3d_parameters_type) :: params
      real(wp), allocatable :: sr(:, :), lr(:, :), reference(:, :), start(:, :), zero(:, :)
      complex(wp), allocatable :: uk(:, :)
      logical, allocatable :: matched(:)
      real(wp) :: rmse
      integer :: mode, i, col, cold_steps, n_masked, n_mixed, n_none
      logical :: converged

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_uv3d_grid(1, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_uv3d_data(grid, response, sr, lr, uk)
      allocate (matched(grid%ngrid))
      do i = 1, grid%ngrid
         matched(i) = mod(i, 3) /= 0
      end do
      params%solver%jfnk%restart = 8
      do mode = 1, size(modes)
         params%cold_guess = modes(mode)
         call new_uv3d_problem(cold, ctx, grid, response, sr, lr, uk, rho, beta, params, err)
         if (.not. allocated(err)) call cold%solve(1.0e-12_wp, 128, converged, rmse, cold_steps, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged, more="cold reference "//modes(mode))
         if (allocated(error)) return
         ! Perturbed so that the guess still needs iterations
         reference = cold%tr*(1.0_wp + 0.05_wp*sin(spread(sr(:, 1), 2, 2)))

         call new_uv3d_problem(mixed, ctx, grid, response, sr, lr, uk, rho, beta, params, err)
         if (allocated(err)) exit
         allocate (start, mold=reference)
         if (modes(mode) == "cs0") then
            allocate (zero, mold=reference)
            zero = 0.0_wp
            call mixed%equation%map(zero, start, err)
            deallocate (zero)
            if (allocated(err)) exit
         else
            start = 0.0_wp
         end if
         do col = 1, size(start, 2)
            where (matched) start(:, col) = reference(:, col)
         end do
         call mixed%solve(1.0e-12_wp, 128, converged, rmse, n_mixed, err, guess=start)
         deallocate (start)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged, more="hand-built start "//modes(mode))
         if (allocated(error)) return

         call new_uv3d_problem(masked, ctx, grid, response, sr, lr, uk, rho, beta, params, err)
         if (.not. allocated(err)) call masked%solve(1.0e-12_wp, 128, converged, rmse, n_masked, err, &
            & guess=reference, matched=matched)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged .and. n_masked == n_mixed, &
            & more="masked guess did not start from the cold field "//modes(mode))
         if (allocated(error)) return
         call check(error, maxval(abs(masked%tr - mixed%tr)), 0.0_wp, thr=1.0e-12_wp, &
            & more="masked and hand-built starts differ "//modes(mode))
         if (allocated(error)) return

         call new_uv3d_problem(none, ctx, grid, response, sr, lr, uk, rho, beta, params, err)
         if (.not. allocated(err)) call none%solve(1.0e-12_wp, 128, converged, rmse, n_none, err, &
            & guess=reference, matched=spread(.false., 1, grid%ngrid))
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged .and. n_none == cold_steps, &
            & more="a guess matching no point was not a cold start "//modes(mode))
         if (allocated(error)) return
         call check(error, maxval(abs(none%tr - cold%tr)), 0.0_wp, thr=1.0e-12_wp, &
            & more="unmatched start differs from the cold solve "//modes(mode))
         if (allocated(error)) return
      end do
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine check_partial

   !> Guesses of the wrong shape or non-finite at a matched point start cold
   !>
   !> A non-finite value at an unmatched point is not used and does not reject the guess
   !>
   !> @param[out] error  test failure
   subroutine check_rejected(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type) :: problem
      type(model_moz_3d_parameters_type) :: params
      real(wp), allocatable :: sr(:, :), lr(:, :), reference(:, :), bad(:, :)
      complex(wp), allocatable :: uk(:, :)
      logical, allocatable :: matched(:)
      real(wp) :: rmse
      integer :: cold_steps, niter
      logical :: converged

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_uv3d_grid(1, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_uv3d_data(grid, response, sr, lr, uk)
      call new_uv3d_problem(problem, ctx, grid, response, sr, lr, uk, rho, beta, params, err)
      if (.not. allocated(err)) call problem%solve(1.0e-12_wp, 128, converged, rmse, cold_steps, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. cold_steps > 1, more="cold reference")
      if (allocated(error)) return
      reference = problem%tr

      call problem%solve(1.0e-12_wp, 128, converged, rmse, niter, err, guess=reference(2:, :))
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == cold_steps, more="wrong guess shape was used")
      if (allocated(error)) return

      call problem%solve(1.0e-12_wp, 128, converged, rmse, niter, err, guess=reference, &
         & matched=spread(.true., 1, grid%ngrid - 1))
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == cold_steps, more="wrong mask size was used")
      if (allocated(error)) return

      bad = reference
      bad(1, 2) = ieee_value(0.0_wp, ieee_quiet_nan)
      call problem%solve(1.0e-12_wp, 128, converged, rmse, niter, err, guess=bad)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == cold_steps, more="non-finite guess was used")
      if (allocated(error)) return

      allocate (matched(grid%ngrid), source=.true.)
      matched(1) = .false.
      call problem%solve(1.0e-12_wp, 128, converged, rmse, niter, err, guess=bad, matched=matched)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter < cold_steps, &
         & more="non-finite value at an unmatched point rejected the guess")
      if (allocated(error)) return
      call check(error, maxval(abs(problem%tr - reference)), 0.0_wp, thr=1.0e-10_wp, more="masked guess solution")
      call ctx%delete()
   end subroutine check_rejected

   !> A failed warm attempt from a masked guess retries entirely cold in the selected mode
   !>
   !> Identity response: the equation map ignores c, so damped Picard steps from a far
   !> guess run out of iterations while the cold start converges
   !>
   !> @param[out] error  test failure
   subroutine check_failed_warm(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type) :: problem
      type(model_moz_3d_parameters_type) :: params
      real(wp), allocatable :: sr(:, :), lr(:, :), reference(:, :)
      complex(wp), allocatable :: uk(:, :)
      logical, allocatable :: matched(:)
      real(wp) :: rmse
      integer :: i, mode, cold_steps, niter
      logical :: converged

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_uv3d_grid(1, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_uv3d_data(grid, response, sr, lr, uk)
      response%values = 0.0_wp
      ! The far but feasible seed exercises an exhausted warm attempt, not HNC's
      ! immediate rejection of an overflowing exponential.
      params%closure = "KH"
      do i = 1, size(response%values, 3)
         response%values(1, 1, i) = 1.0_wp
         response%values(2, 2, i) = 1.0_wp
      end do
      allocate (matched(grid%ngrid))
      do i = 1, grid%ngrid
         matched(i) = mod(i, 2) == 0
      end do
      params%solver%name = "PICARD"
      params%solver%picard%damping = 0.5_wp
      do mode = 1, size(modes)
         params%cold_guess = modes(mode)
         call new_uv3d_problem(problem, ctx, grid, response, sr, lr, uk, rho, beta, params, err)
         if (.not. allocated(err)) call problem%solve(1.0e-8_wp, 32, converged, rmse, cold_steps, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged, more="cold reference "//modes(mode))
         if (allocated(error)) return
         reference = problem%tr
         call problem%solve(1.0e-8_wp, 32, converged, rmse, niter, err, guess=reference + 1.0e8_wp, matched=matched)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged .and. niter == 32 + cold_steps, &
            & more="failed warm attempt did not retry cold "//modes(mode))
         if (allocated(error)) return
         call check(error, maxval(abs(problem%tr - reference)), 0.0_wp, thr=1.0e-13_wp, &
            & more="cold retry changed the cold solution "//modes(mode))
         if (allocated(error)) return
      end do
      call ctx%delete()
   end subroutine check_failed_warm

   !> A solve that stops on a closure switch converges only on the requested closure
   !>
   !> - A tolerance far above every residual makes the HNC ramp switch after each step:
   !>   KH, PSE2, PSE3, HNC
   !> - One step ends on KH -> PSE2: the residual meets the tolerance but the solve is
   !>   not converged
   !> - Three steps end on PSE3 -> HNC: converged, with the residual of the kept field
   !>   measured under HNC
   !>
   !> @param[out] error  test failure
   subroutine ramp_switch_end(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type), target :: problem
      type(uv3d_evaluator_type) :: evaluator
      type(model_moz_3d_parameters_type) :: params
      real(wp), allocatable :: sr(:, :), lr(:, :), residual(:, :)
      complex(wp), allocatable :: uk(:, :)
      real(wp) :: rmse, reference
      integer :: niter
      logical :: converged, feasible
      real(wp), parameter :: loose = 1.0e10_wp

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_uv3d_grid(1, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_uv3d_data(grid, response, sr, lr, uk)
      params%closure = "HNC"
      params%closure_ramp = .true.
      params%solver%name = "PICARD"
      params%solver%picard%damping = 0.5_wp
      call new_uv3d_problem(problem, ctx, grid, response, sr, lr, uk, rho, beta, params, err)
      if (.not. allocated(err)) call problem%solve(loose, 1, converged, rmse, niter, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. converged .and. niter == 1 .and. rmse < loose, &
         & more="a solve stopped on an intermediate closure converged")
      if (allocated(error)) return
      call new_uv3d_problem(problem, ctx, grid, response, sr, lr, uk, rho, beta, params, err)
      if (.not. allocated(err)) call problem%solve(loose, 3, converged, rmse, niter, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == 3, more="a solve stopped on the switch to HNC did not converge")
      if (allocated(error)) return
      evaluator%problem => problem
      allocate (residual, mold=problem%tr)
      call evaluator%residual(problem%tr, residual, feasible, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      reference = sqrt(sum(residual**2)/real(size(residual), wp))
      call check(error, feasible .and. abs(rmse - reference) <= 1.0e-12_wp*max(1.0_wp, reference), &
         & more="the residual after the switch was not measured under HNC")
      call ctx%delete()
   end subroutine ramp_switch_end

end module test_moz_3d_guess
