!> Numerical controls shared by the MOZ solver factory
!>
!> JSON/TOML keys: verbosity, nthreads, picard.damping, diis.*, hybrid.*, jfnk.*, adjoint.*
module moist_model_moz_solver_parameters
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_model_parameters, only: moist_model_parameters_type
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_parameters_type

   !> Controls for Picard, DIIS, JFNK and HYBRID
   type, extends(moist_model_parameters_type) :: moist_model_moz_solver_parameters_type
      !> Diagnostic verbosity; zero is quiet
      integer :: verbosity = 0
      !> OpenMP team size; models should pass ctx%get_num_threads(), zero takes omp_get_max_threads
      integer :: nthreads = 0
      !> Picard damping, also used while DIIS fills its history
      real(wp) :: picard_damping = 0.1_wp
      !> DIIS history size
      integer :: diis_depth = 12
      !> DIIS damping when backtracking is disabled
      real(wp) :: diis_damping = 0.2_wp
      !> Bound the DIIS update relative to the fixed-point residual
      logical :: diis_backtracking = .true.
      !> DIIS backtracking factor, in (0, 1)
      real(wp) :: diis_backtrack_tau = 0.5_wp
      !> Smallest DIIS step fraction before the Picard update is used
      real(wp) :: diis_min_alpha = 1.0e-3_wp
      !> Largest DIIS step relative to the fixed-point residual
      real(wp) :: diis_max_step_ratio = 200.0_wp
      !> Hybrid DIIS-to-JFNK switch at or below this RMS residual
      real(wp) :: hybrid_switch_rmse = 1.0e-3_wp
      !> GMRES restart length
      integer :: jfnk_restart = 50
      !> Inner GMRES iterations per Newton update
      integer :: jfnk_max_krylov = 50
      !> Transpose GMRES restart length, independent of the nonlinear solve
      integer :: adjoint_restart = 50
      !> Inner GMRES iterations per transpose refinement
      integer :: adjoint_max_krylov = 50
      !> Forcing term of the first two Newton updates of a stage
      real(wp) :: jfnk_eta0 = 0.5_wp
      !> Lower bound of later forcing terms
      real(wp) :: jfnk_eta_min = 1.0e-3_wp
      !> Upper bound of later forcing terms
      real(wp) :: jfnk_eta_max = 0.1_wp
      !> Finite-difference displacement; zero selects an automatic value
      real(wp) :: jfnk_fd_eps = 0.0_wp
      !> Prefer exact matrix-free products when the evaluator provides them
      logical :: jfnk_exact_jvp = .true.
      !> Classical Gram-Schmidt with one reorthogonalization instead of MGS
      logical :: jfnk_cgs2 = .false.
      !> Recycled Krylov vectors (GCRO-DR); zero disables recycling
      integer :: jfnk_recycle = 0
      !> Coarse-space (DEF1) deflation of the inner solves
      logical :: jfnk_deflation = .false.
      !> Deflation is used below this RMS residual
      real(wp) :: jfnk_deflation_rmse = 1.0e-3_wp
      !> Relative decrease below which a restart cycle ends the inner solve; zero disables it
      real(wp) :: jfnk_stagnation_tol = 0.0_wp
      !> Maximum line-search fraction; one disables step extension
      real(wp) :: jfnk_alpha_max = 2.0_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Declare fields for JSON/TOML input, output and printing
      procedure :: register_entries => register_parameter_entries
      !> Check the ranges of all controls
      procedure :: validate => validate_solver_parameters
   end type moist_model_moz_solver_parameters_type

contains

   !> Restore compiled defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_parameter_defaults(self)
      !> Parameter values
      class(moist_model_moz_solver_parameters_type), intent(inout) :: self
      type(moist_model_moz_solver_parameters_type) :: defaults

      self%verbosity = defaults%verbosity
      self%nthreads = defaults%nthreads
      self%picard_damping = defaults%picard_damping
      self%diis_depth = defaults%diis_depth
      self%diis_damping = defaults%diis_damping
      self%diis_backtracking = defaults%diis_backtracking
      self%diis_backtrack_tau = defaults%diis_backtrack_tau
      self%diis_min_alpha = defaults%diis_min_alpha
      self%diis_max_step_ratio = defaults%diis_max_step_ratio
      self%hybrid_switch_rmse = defaults%hybrid_switch_rmse
      self%jfnk_restart = defaults%jfnk_restart
      self%jfnk_max_krylov = defaults%jfnk_max_krylov
      self%adjoint_restart = defaults%adjoint_restart
      self%adjoint_max_krylov = defaults%adjoint_max_krylov
      self%jfnk_eta0 = defaults%jfnk_eta0
      self%jfnk_eta_min = defaults%jfnk_eta_min
      self%jfnk_eta_max = defaults%jfnk_eta_max
      self%jfnk_fd_eps = defaults%jfnk_fd_eps
      self%jfnk_exact_jvp = defaults%jfnk_exact_jvp
      self%jfnk_cgs2 = defaults%jfnk_cgs2
      self%jfnk_recycle = defaults%jfnk_recycle
      self%jfnk_deflation = defaults%jfnk_deflation
      self%jfnk_deflation_rmse = defaults%jfnk_deflation_rmse
      self%jfnk_stagnation_tol = defaults%jfnk_stagnation_tol
      self%jfnk_alpha_max = defaults%jfnk_alpha_max
   end subroutine init_parameter_defaults

   !> Declare parameter fields for JSON/TOML input, output and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_parameter_entries(self)
      !> Parameter values
      class(moist_model_moz_solver_parameters_type), intent(inout), target :: self

      call self%register_int_scalar("verbosity", self%verbosity)
      call self%register_int_scalar("nthreads", self%nthreads)
      call self%register_real_scalar("picard.damping", self%picard_damping)
      call self%register_int_scalar("diis.depth", self%diis_depth)
      call self%register_real_scalar("diis.damping", self%diis_damping)
      call self%register_logical("diis.backtracking", self%diis_backtracking)
      call self%register_real_scalar("diis.backtrack_tau", self%diis_backtrack_tau)
      call self%register_real_scalar("diis.min_alpha", self%diis_min_alpha)
      call self%register_real_scalar("diis.max_step_ratio", self%diis_max_step_ratio)
      call self%register_real_scalar("hybrid.switch_rmse", self%hybrid_switch_rmse)
      call self%register_int_scalar("jfnk.restart", self%jfnk_restart)
      call self%register_int_scalar("jfnk.max_krylov", self%jfnk_max_krylov)
      call self%register_int_scalar("adjoint.restart", self%adjoint_restart)
      call self%register_int_scalar("adjoint.max_krylov", self%adjoint_max_krylov)
      call self%register_real_scalar("jfnk.eta0", self%jfnk_eta0)
      call self%register_real_scalar("jfnk.eta_min", self%jfnk_eta_min)
      call self%register_real_scalar("jfnk.eta_max", self%jfnk_eta_max)
      call self%register_real_scalar("jfnk.fd_eps", self%jfnk_fd_eps)
      call self%register_logical("jfnk.exact_jvp", self%jfnk_exact_jvp)
      call self%register_logical("jfnk.cgs2", self%jfnk_cgs2)
      call self%register_int_scalar("jfnk.recycle", self%jfnk_recycle)
      call self%register_logical("jfnk.deflation", self%jfnk_deflation)
      call self%register_real_scalar("jfnk.deflation_rmse", self%jfnk_deflation_rmse)
      call self%register_real_scalar("jfnk.stagnation_tol", self%jfnk_stagnation_tol)
      call self%register_real_scalar("jfnk.alpha_max", self%jfnk_alpha_max)
   end subroutine register_parameter_entries

   !> Check the ranges of all controls
   !>
   !> Backtracking and deflation controls are checked only when enabled
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error First invalid control
   subroutine validate_solver_parameters(self, error)
      !> Parameter values
      class(moist_model_moz_solver_parameters_type), intent(inout) :: self
      !> First invalid control
      type(error_type), allocatable, intent(out) :: error

      if (self%verbosity < 0) then
         call fatal_error(error, "Solver verbosity must be nonnegative")
      else if (self%nthreads < 0) then
         call fatal_error(error, "Solver thread count must be nonnegative")
      else if (.not. positive(self%picard_damping)) then
         call fatal_error(error, "Picard damping must be finite and positive")
      else if (self%diis_depth < 1 .or. .not. positive(self%diis_damping)) then
         call fatal_error(error, "DIIS needs a positive history size and finite positive damping")
      else if (self%diis_backtracking .and. (.not. unit_fraction(self%diis_backtrack_tau) .or. &
         & .not. positive(self%diis_min_alpha) .or. self%diis_min_alpha > 1.0_wp .or. &
         & .not. positive(self%diis_max_step_ratio))) then
         call fatal_error(error, "Invalid DIIS backtracking controls")
      else if (.not. positive(self%hybrid_switch_rmse)) then
         call fatal_error(error, "Hybrid switch RMS residual must be finite and positive")
      else if (self%jfnk_restart < 1 .or. self%jfnk_max_krylov < 1 .or. self%jfnk_recycle < 0) then
         call fatal_error(error, "JFNK needs positive Krylov limits and a nonnegative recycling dimension")
      else if (self%adjoint_restart < 1 .or. self%adjoint_max_krylov < 1) then
         call fatal_error(error, "Adjoint needs positive Krylov limits")
      else if (.not. unit_fraction(self%jfnk_eta0) .or. .not. unit_fraction(self%jfnk_eta_min) .or. &
         & .not. unit_fraction(self%jfnk_eta_max) .or. self%jfnk_eta_min > self%jfnk_eta_max) then
         call fatal_error(error, "JFNK forcing terms must lie between zero and one, with eta_min <= eta_max")
      else if (.not. ieee_is_finite(self%jfnk_fd_eps) .or. self%jfnk_fd_eps < 0.0_wp .or. &
         & .not. ieee_is_finite(self%jfnk_stagnation_tol) .or. self%jfnk_stagnation_tol < 0.0_wp .or. &
         & self%jfnk_stagnation_tol >= 1.0_wp .or. .not. ieee_is_finite(self%jfnk_alpha_max) .or. &
         & self%jfnk_alpha_max < 1.0_wp) then
         call fatal_error(error, "Invalid JFNK finite-difference, stagnation or line-search controls")
      else if (self%jfnk_deflation .and. .not. positive(self%jfnk_deflation_rmse)) then
         call fatal_error(error, "JFNK deflation threshold must be finite and positive")
      end if
   end subroutine validate_solver_parameters

   !> Finite positive value
   !>
   !> @param[in] value Scalar to check
   pure logical function positive(value)
      !> Scalar to check
      real(wp), intent(in) :: value
      positive = ieee_is_finite(value) .and. value > 0.0_wp
   end function positive

   !> Fraction strictly between zero and one
   !>
   !> @param[in] value Fraction to check
   pure logical function unit_fraction(value)
      !> Fraction to check
      real(wp), intent(in) :: value
      unit_fraction = positive(value) .and. value < 1.0_wp
   end function unit_fraction

end module moist_model_moz_solver_parameters
