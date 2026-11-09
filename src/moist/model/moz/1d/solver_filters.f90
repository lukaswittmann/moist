!> Small-k filters of the 1D UV solves, owned by the model and lent to its solver
!>
!> - spectral: JFNK left preconditioner when selected, and always JFNK's Picard fallback filter
!> - damping: filter of the updates PICARD and DIIS (also HYBRID's DIIS phase) allow to post-process
!> - Both apply w(k) = w0 + (1 - w0) k^2/(k^2 + k0^2) with their own k0 and weight
!> - Built from settings values, not the model parameters, which already use this module's names
!> - check refuses filters that were not set up or belong to another radial pair
module moist_model_moz_1d_solver_filters
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_model_moz_solver, only: moist_model_moz_solver_type, moist_model_moz_solver_preconditioner_type, &
      & check_preconditioner_name
   use moist_model_moz_1d_preconditioner, only: moist_model_moz_1d_preconditioner_spectral_type, &
      & select_1d_preconditioner, preconditioners_1d_uv
   use moist_model_moz_1d_operators, only: same_radial_grid
   implicit none(type, external)
   private

   public :: moist_model_moz_1d_solver_filters_type

   !> Preconditioner and damping filters of the 1D UV solves
   !>
   !> - attach lends them to a solver for one solve; detach keeps them for the response and later solves
   !> - Targets of the lent pointers: the object must not move while attached
   type :: moist_model_moz_1d_solver_filters_type
      !> Lowercase preconditioner name: none or spectral
      character(len=8) :: precond = "none"
      !> Spectral filter: preconditioner when selected and the Picard fallback filter
      type(moist_model_moz_1d_preconditioner_spectral_type) :: spectral
      !> Whether post-processed updates are damped
      logical :: damp_small_k = .false.
      !> Small-k damping filter
      type(moist_model_moz_1d_preconditioner_spectral_type) :: damping
      !> Flat update of the damping, reused across iterations and solves
      real(wp), allocatable, private :: delta(:)
      !> Radial pair of the last successful setup
      type(moist_math_grid_radial_type), private :: rgrid, kgrid
      !> Whether setup succeeded
      logical, private :: built = .false.
   contains
      !> Build both filters on a radial pair
      procedure :: setup => filters_setup
      !> Refuse filters that are not set up or belong to another radial pair
      procedure :: check => filters_check
      !> Lend the preconditioner and the fallback filter to a solver
      procedure :: attach => filters_attach
      !> Take both operators back from a solver
      procedure :: detach => filters_detach
      !> Left preconditioner of a response solve, null for none
      procedure :: response_preconditioner => filters_response_preconditioner
      !> Damp an update the solver allows to post-process
      procedure :: damp => filters_damp
   end type moist_model_moz_1d_solver_filters_type

contains

   !> Build both filters on a radial pair
   !>
   !> @param[out] self             filters; damping is built only when enabled
   !> @param[in]  rgrid            real-space uniform DST-IV grid
   !> @param[in]  kgrid            reciprocal grid of the pair
   !> @param[in]  nthreads         OpenMP team size of the filters
   !> @param[in]  precond          preconditioner name: none or spectral
   !> @param[in]  precond_k0       spectral-filter wave number, 1/bohr
   !> @param[in]  precond_w0       spectral-filter weight at k = 0
   !> @param[in]  damp_small_k     whether post-processed updates are damped
   !> @param[in]  damp_k0          damping wave number, 1/bohr
   !> @param[in]  damp_min_weight  damping weight at k = 0
   !> @param[out] error            unknown preconditioner name or filter setup failure
   subroutine filters_setup(self, rgrid, kgrid, nthreads, precond, precond_k0, precond_w0, damp_small_k, damp_k0, &
      & damp_min_weight, error)
      !> Filters
      class(moist_model_moz_1d_solver_filters_type), intent(out) :: self
      !> Radial pair
      type(moist_math_grid_radial_type), intent(in) :: rgrid, kgrid
      !> OpenMP team size
      integer, intent(in) :: nthreads
      !> Preconditioner name
      character(len=*), intent(in) :: precond
      !> Spectral filter
      real(wp), intent(in) :: precond_k0, precond_w0
      !> Damping switch
      logical, intent(in) :: damp_small_k
      !> Damping filter
      real(wp), intent(in) :: damp_k0, damp_min_weight
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_preconditioner_name(precond, preconditioners_1d_uv, "1D UV", error)
      if (allocated(error)) return
      self%precond = to_lower(trim(adjustl(precond)))
      self%spectral%nthreads = nthreads
      call self%spectral%setup(rgrid, kgrid, precond_k0, precond_w0, error)
      if (allocated(error)) return
      self%damp_small_k = damp_small_k
      if (self%damp_small_k) then
         self%damping%nthreads = nthreads
         call self%damping%setup(rgrid, kgrid, damp_k0, damp_min_weight, error)
         if (allocated(error)) return
      end if
      self%rgrid = rgrid
      self%kgrid = kgrid
      self%built = .true.
   end subroutine filters_setup

   !> Refuse filters that are not set up or belong to another radial pair
   !>
   !> An unbuilt filter would leave vectors unchanged, and one from another pair would
   !> apply other weights; both silently change the solve
   !>
   !> @param[in]  self   filters
   !> @param[in]  rgrid  real-space grid of the problem
   !> @param[in]  kgrid  reciprocal grid of the problem
   !> @param[out] error  filters not set up or built on another pair
   subroutine filters_check(self, rgrid, kgrid, error)
      !> Filters
      class(moist_model_moz_1d_solver_filters_type), intent(in) :: self
      !> Radial pair of the problem
      type(moist_math_grid_radial_type), intent(in) :: rgrid, kgrid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%built) then
         call fatal_error(error, "1D UV filters are not set up")
      else if (.not. (same_radial_grid(self%rgrid, rgrid) .and. same_radial_grid(self%kgrid, kgrid))) then
         call fatal_error(error, "1D UV filters belong to another radial grid pair")
      end if
   end subroutine filters_check

   !> Lend the preconditioner and the fallback filter to a solver
   !>
   !> - Only solvers that use a preconditioner take them
   !> - The fallback filter is lent also when no preconditioner is selected
   !>
   !> @param[inout] self    filters; must stay in place until detach
   !> @param[inout] solver  solver of the next solve
   subroutine filters_attach(self, solver)
      !> Filters
      class(moist_model_moz_1d_solver_filters_type), intent(inout), target :: self
      !> Solver
      class(moist_model_moz_solver_type), intent(inout) :: solver
      class(moist_model_moz_solver_preconditioner_type), pointer :: precond

      if (.not. solver%uses_preconditioner()) return
      call select_1d_preconditioner(self%precond, self%spectral, precond)
      if (associated(precond)) call solver%set_preconditioner(precond)
      call solver%set_fallback_filter(self%spectral)
   end subroutine filters_attach

   !> Take both operators back from a solver; the filters stay built
   !>
   !> @param[in]    self    filters
   !> @param[inout] solver  solver of the finished solve
   subroutine filters_detach(self, solver)
      !> Filters
      class(moist_model_moz_1d_solver_filters_type), intent(in) :: self
      !> Solver
      class(moist_model_moz_solver_type), intent(inout) :: solver

      associate (unused => self); end associate
      call solver%set_preconditioner()
      call solver%set_fallback_filter()
   end subroutine filters_detach

   !> Left preconditioner of a response solve, whatever solver ran the primal solve
   !>
   !> @param[inout] self  filters; must stay in place while the result is used
   function filters_response_preconditioner(self) result(precond)
      !> Filters
      class(moist_model_moz_1d_solver_filters_type), intent(inout), target :: self
      !> Selected operator, null for none
      class(moist_model_moz_solver_preconditioner_type), pointer :: precond

      call select_1d_preconditioner(self%precond, self%spectral, precond)
   end function filters_response_preconditioner

   !> Damp the update gamma_new - gamma at small k when enabled and the solver allows it
   !>
   !> @param[inout] self       filters; keeps the update scratch
   !> @param[in]    solver     solver that proposed the update
   !> @param[in]    gamma      input iterate, (npts, npair)
   !> @param[inout] gamma_new  proposed iterate; damped on return
   subroutine filters_damp(self, solver, gamma, gamma_new)
      !> Filters
      class(moist_model_moz_1d_solver_filters_type), intent(inout) :: self
      !> Solver
      class(moist_model_moz_solver_type), intent(in) :: solver
      !> Input iterate
      real(wp), intent(in) :: gamma(:, :)
      !> Proposed iterate
      real(wp), intent(inout) :: gamma_new(:, :)

      if (.not. (self%damp_small_k .and. solver%allows_postprocessing())) return
      if (allocated(self%delta)) then
         if (size(self%delta) /= size(gamma)) deallocate (self%delta)
      end if
      if (.not. allocated(self%delta)) allocate (self%delta(size(gamma)))
      self%delta = reshape(gamma_new - gamma, [size(self%delta)])
      call self%damping%apply(self%delta)
      gamma_new = gamma + reshape(self%delta, shape(gamma))
   end subroutine filters_damp

end module moist_model_moz_1d_solver_filters
