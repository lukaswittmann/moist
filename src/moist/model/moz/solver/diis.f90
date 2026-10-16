!> DIIS iteration using modified direct inversion in the iterative subspace
module moist_model_moz_solver_diis
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_model_moz_solver_base, only: moist_model_moz_solver_type, validate_solver_arrays
   use moist_model_moz_solver_base, only: moist_model_moz_solver_statistics_type
   use moist_math_lapack_gesv, only: dgesv
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_blas, only: gemv, syrk
   use moist_model_moz_solver_evaluator, only: moist_model_moz_solver_evaluator_type
   use moist_model_moz_solver_printer, only: moist_model_moz_solver_printer_type
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_diis_type, new_diis_solver

   !> Residual-minimizing mixing with a bounded circular history
   type, extends(moist_model_moz_solver_type) :: moist_model_moz_solver_diis_type
      !> Maximum stored iterates and residuals
      integer :: depth = 12
      !> Extrapolation damping when backtracking is disabled
      real(wp) :: damping = 0.2_wp
      !> Damping while filling history and on fallback
      real(wp) :: picard_damping = 0.1_wp
      !> Bound the extrapolated update relative to the residual
      logical :: backtracking = .true.
      !> Step reduction factor during backtracking
      real(wp) :: backtrack_tau = 0.5_wp
      !> Minimum accepted step fraction before fallback and history reset
      real(wp) :: min_alpha = 1.0e-3_wp
      !> Upper bound for step RMSE relative to fixed-point residual RMSE
      real(wp) :: max_step_ratio = 200.0_wp
      !> Stored iterates (nflat, depth)
      real(wp), allocatable :: fr(:, :)
      !> Stored fixed-point updates, equal to minus the evaluator residual
      real(wp), allocatable :: res(:, :)
      !> RMS residuals of stored iterates (depth)
      real(wp), allocatable :: rms(:)
      !> Number of occupied history slots
      integer :: nstored = 0
      !> Number of elements in each stored field
      integer :: nflat = 0
      !> Next slot to replace after the history fills
      integer :: history_head = 1
      !> Effective step fraction from the last update (for logging)
      real(wp) :: last_alpha = 0.0_wp
      !> Borrowed iteration printer, detached by omitting the setter argument
      type(moist_model_moz_solver_printer_type), pointer :: iter_printer => null()
   contains
      procedure, private :: mix => diis_mix
      procedure :: get_statistics => diis_get_statistics
      procedure :: outer_step => diis_outer_step
      procedure :: name => diis_name
      procedure :: begin_stage => diis_begin_stage
      procedure :: reset => diis_reset
      procedure :: set_iteration_printer => diis_set_iteration_printer
   end type moist_model_moz_solver_diis_type

contains

   !> Display name of the solver
   !>
   !> @param[in] self Solver to name
   pure function diis_name(self) result(text)
      !> Solver to name
      class(moist_model_moz_solver_diis_type), intent(in) :: self
      !> Display name
      character(len=:), allocatable :: text
      text = "DIIS"
   end function diis_name

   !> Return an allocation-free snapshot of solver statistics
   !>
   !> @param[in] self Solver to inspect
   !> @return stats Counters and cumulative timings, or unavailable sentinels
   function diis_get_statistics(self) result(stats)
      !> Solver to inspect
      class(moist_model_moz_solver_diis_type), intent(in) :: self
      !> Statistics snapshot
      type(moist_model_moz_solver_statistics_type) :: stats
      stats = moist_model_moz_solver_statistics_type()
      stats%ndiis = self%nstored
   end function diis_get_statistics

   !> Construct a DIIS solver with optional history and mixing controls
   !>
   !> @param[in] depth Maximum stored iterates and residuals
   !> @param[in] damping Extrapolation damping without backtracking
   !> @param[in] picard_damping Warm-up and fallback damping
   !> @param[in] backtracking Whether to bound extrapolated updates
   !> @param[in] backtrack_tau Backtracking reduction factor
   !> @param[in] min_alpha Minimum accepted step fraction
   !> @param[in] max_step_ratio Maximum update RMS relative to residual RMS
   !> @param[in] verbosity Diagnostic verbosity
   function new_diis_solver(depth, damping, picard_damping, backtracking, &
                             backtrack_tau, min_alpha, max_step_ratio, verbosity) result(self)
      !> Maximum stored iterates and residuals
      integer, intent(in), optional :: depth
      !> Extrapolation damping without backtracking
      real(wp), intent(in), optional :: damping
      !> Warm-up and fallback damping
      real(wp), intent(in), optional :: picard_damping
      !> Whether to bound extrapolated updates
      logical, intent(in), optional :: backtracking
      !> Backtracking reduction factor
      real(wp), intent(in), optional :: backtrack_tau
      !> Minimum accepted step fraction
      real(wp), intent(in), optional :: min_alpha
      !> Maximum update RMS relative to residual RMS
      real(wp), intent(in), optional :: max_step_ratio
      !> Diagnostic verbosity
      integer, intent(in), optional :: verbosity
      !> Constructed solver
      type(moist_model_moz_solver_diis_type) :: self

      if (present(depth)) self%depth = depth
      if (present(damping)) self%damping = damping
      if (present(picard_damping)) self%picard_damping = picard_damping
      if (present(backtracking)) self%backtracking = backtracking
      if (present(backtrack_tau)) self%backtrack_tau = backtrack_tau
      if (present(min_alpha)) self%min_alpha = min_alpha
      if (present(max_step_ratio)) self%max_step_ratio = max_step_ratio
      if (present(verbosity)) self%verbosity = verbosity
   end function new_diis_solver

   !> Borrow a printer whose target outlives the attachment; omit to detach
   !>
   !> @param[inout] self Solver to configure
   !> @param[in] printer Borrowed printer
   subroutine diis_set_iteration_printer(self, printer)
      !> Solver to configure
      class(moist_model_moz_solver_diis_type), intent(inout) :: self
      !> Borrowed printer
      type(moist_model_moz_solver_printer_type), target, optional, intent(in) :: printer
      nullify (self%iter_printer)
      if (present(printer)) self%iter_printer => printer
   end subroutine diis_set_iteration_printer

   !> Discard residuals from the previous equation while retaining allocated storage
   !>
   !> @param[inout] self Solver starting a new solve or coupling stage
   subroutine diis_begin_stage(self)
      !> DIIS history and update diagnostics
      class(moist_model_moz_solver_diis_type), intent(inout) :: self

      self%nstored = 0
      self%history_head = 1
      self%last_alpha = 0.0_wp
      self%prev_alpha = 0.0_wp
      self%prev_step_maxabs = 0.0_wp
   end subroutine diis_begin_stage

   !> Clear mixing history and release its storage
   !>
   !> @param[inout] self Solver history to clear
   subroutine diis_reset(self)
      !> Solver history to clear
      class(moist_model_moz_solver_diis_type), intent(inout) :: self

      call self%begin_stage()
      if (allocated(self%fr)) deallocate (self%fr)
      if (allocated(self%res)) deallocate (self%res)
      if (allocated(self%rms)) deallocate (self%rms)
      self%nflat = 0
   end subroutine diis_reset

   !> Store the evaluator residual and extrapolate from the current history
   !>
   !> @param[inout] self Solver and mixing history
   !> @param[in] residual Evaluator residual, c_prev - c_map (npts, npair)
   !> @param[in] c_prev Current iterate (npts, npair)
   !> @param[out] c_next Mixed iterate (npts, npair)
   subroutine diis_mix(self, residual, c_prev, c_next)
      !> Solver and mixing history
      class(moist_model_moz_solver_diis_type), intent(inout) :: self
      !> Evaluator residual (npts, npair)
      real(wp), intent(in) :: residual(:, :)
      !> Current iterate (npts, npair)
      real(wp), intent(in) :: c_prev(:, :)
      !> Mixed iterate (npts, npair)
      real(wp), intent(out) :: c_next(:, :)

      integer :: nflat, m, i, j, ip, idx, npts, npair, slot
      integer(lapack_ik) :: info
      logical :: filling_history
      integer(lapack_ik), allocatable :: ipiv(:)
      real(wp), allocatable :: A(:, :), A_gram(:, :), b(:)
      real(wp), allocatable :: diff(:), c_flat(:), c_comb(:), res_comb(:)
      real(wp) :: rmse, step_rmse, step_rmse_base, max_step_rmse, alpha_bt, sumsq

      nflat = size(c_prev)
      m = self%depth
      npts = size(c_prev, 1)
      npair = size(c_prev, 2)

      if (allocated(self%fr)) then
         if (self%nflat /= nflat .or. size(self%fr, 2) /= m) call self%reset()
      end if
      if (.not. allocated(self%fr)) then
         allocate (self%fr(nflat, m))
         allocate (self%res(nflat, m))
         allocate (self%rms(m))
         self%nflat = nflat
      end if

      allocate (diff(nflat), c_flat(nflat))
      idx = 0
      sumsq = 0.0_wp
      do ip = 1, npair
         do i = 1, npts
            idx = idx + 1
            ! Keep the supplied residual even when c_prev - residual rounds to c_prev.
            diff(idx) = -residual(i, ip)
            c_flat(idx) = c_prev(i, ip)
            sumsq = sumsq + diff(idx)*diff(idx)
         end do
      end do
      rmse = sqrt(sumsq/real(nflat, wp))

      filling_history = self%nstored < m
      if (filling_history) then
         slot = self%nstored + 1
         self%nstored = slot
      else
         slot = self%history_head
         self%history_head = modulo(self%history_head, m) + 1
      end if

      self%fr(:, slot) = c_flat
      self%res(:, slot) = diff
      self%rms(slot) = rmse

      if (filling_history) then
         idx = 0
         do ip = 1, npair
            do i = 1, npts
               idx = idx + 1
               c_next(i, ip) = c_prev(i, ip) + self%picard_damping*diff(idx)
            end do
         end do
         self%last_alpha = self%picard_damping
         return
      end if

      if (rmse > 10.0_wp*minval(self%rms(1:self%nstored))) then
         if (self%verbosity >= 2 .and. associated(self%iter_printer)) then
            call self%iter_printer%note("Resetting DIIS")
         end if
         call self%reset()
         idx = 0
         do ip = 1, npair
            do i = 1, npts
               idx = idx + 1
               c_next(i, ip) = c_prev(i, ip) + self%picard_damping*diff(idx)
            end do
         end do
         self%last_alpha = self%picard_damping
         return
      end if

      allocate (A_gram(m, m), A(m + 1, m + 1), b(m + 1), ipiv(m + 1), c_comb(nflat), res_comb(nflat))
      A = 0.0_wp
      A_gram = 0.0_wp
      b = 0.0_wp
      b(m + 1) = -1.0_wp

      call syrk(self%res(:, 1:m), A_gram, uplo="l", trans="t", alpha=1.0_wp, beta=0.0_wp)
      A(1:m, 1:m) = A_gram
      ! Complete the upper triangle from the lower triangle computed by SYRK
      do j = 1, m
         do i = j + 1, m
            A(j, i) = A(i, j)
         end do
      end do

      do j = 1, m + 1
         A(j, m + 1) = -1.0_wp
         A(m + 1, j) = -1.0_wp
      end do
      A(m + 1, m + 1) = 0.0_wp

      call dgesv(int(m + 1, lapack_ik), 1_lapack_ik, A, int(m + 1, lapack_ik), &
         ipiv, b, int(m + 1, lapack_ik), info)

      if (info /= 0) then
         idx = 0
         do ip = 1, npair
            do i = 1, npts
               idx = idx + 1
               c_next(i, ip) = c_prev(i, ip) + self%picard_damping*diff(idx)
            end do
         end do
         self%last_alpha = self%picard_damping
         return
      end if

      call gemv(self%fr(:, 1:m), b(1:m), c_comb, alpha=1.0_wp, beta=0.0_wp, trans="n")
      call gemv(self%res(:, 1:m), b(1:m), res_comb, alpha=1.0_wp, beta=0.0_wp, trans="n")

      c_comb = c_comb + res_comb

      if (self%backtracking) then
         sumsq = 0.0_wp
         do idx = 1, nflat
            res_comb(idx) = c_comb(idx) - c_flat(idx)
            sumsq = sumsq + res_comb(idx)*res_comb(idx)
         end do
         step_rmse_base = sqrt(sumsq/real(nflat, wp))
         max_step_rmse = self%max_step_ratio*max(rmse, tiny(1.0_wp))

         alpha_bt = 1.0_wp
         do while (alpha_bt >= self%min_alpha)
            step_rmse = alpha_bt*step_rmse_base
            if (step_rmse <= max_step_rmse) exit
            alpha_bt = alpha_bt*self%backtrack_tau
         end do

         if (alpha_bt < self%min_alpha) then
            call self%reset()
            idx = 0
            do ip = 1, npair
               do i = 1, npts
                  idx = idx + 1
                  c_next(i, ip) = c_prev(i, ip) + self%picard_damping*diff(idx)
               end do
            end do
            self%last_alpha = self%picard_damping
            return
         end if

         do idx = 1, nflat
            c_comb(idx) = c_flat(idx) + alpha_bt*res_comb(idx)
         end do
         self%last_alpha = alpha_bt
      else
         do idx = 1, nflat
            c_comb(idx) = c_flat(idx) + self%damping*(c_comb(idx) - c_flat(idx))
         end do
         self%last_alpha = self%damping
      end if

      idx = 0
      do ip = 1, npair
         do i = 1, npts
            idx = idx + 1
            c_next(i, ip) = c_comb(idx)
         end do
      end do
   end subroutine diis_mix

   !> Evaluate the residual and mix the resulting fixed-point estimate
   !>
   !> @param[inout] self Solver and mixing history
   !> @param[in] residual_fn Model residual evaluator
   !> @param[in] gamma Current field (npts, npair)
   !> @param[out] gamma_new Mixed field (npts, npair)
   !> @param[out] rmse RMS residual of the input field
   !> @param[out] feasible Whether the residual is usable
   !> @param[out] error Evaluation failure
   subroutine diis_outer_step(self, residual_fn, gamma, gamma_new, rmse, feasible, error)
      !> Solver and mixing history
      class(moist_model_moz_solver_diis_type), intent(inout) :: self
      !> Model residual evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_fn
      !> Current field (npts, npair)
      real(wp), intent(in) :: gamma(:, :)
      !> Mixed field (npts, npair)
      real(wp), intent(out) :: gamma_new(:, :)
      !> RMS residual of the input field
      real(wp), intent(out) :: rmse
      !> Whether the residual is usable
      logical, intent(out) :: feasible
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: F(:, :)

      feasible = .false.
      call validate_solver_arrays(gamma, gamma_new, error)
      if (allocated(error)) return
      allocate(F, mold=gamma)

      call residual_fn%residual(gamma, F, feasible, error)
      if (allocated(error)) feasible = .false.
      if (.not. feasible) return

      rmse = sqrt(sum(F**2) / real(size(F), wp))
      if (.not. ieee_is_finite(rmse)) then
         call fatal_error(error, "Solver residual norm must be finite")
         feasible = .false.
         return
      end if
      call self%mix(F, gamma, gamma_new)
      self%prev_step_maxabs = maxval(abs(gamma_new - gamma))
      self%prev_alpha = self%last_alpha
   end subroutine diis_outer_step

end module moist_model_moz_solver_diis
