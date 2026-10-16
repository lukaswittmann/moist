!> 1D VV solve: residual evaluator, outer iterations and coupling stages
!>
!> - The iterate is the short-range t_s(r); the closure gives c_s(r) and the
!>   equation maps it back, R(t_s) = t_s - T(C(t_s))
!> - Stages lambda = i/n scale the short-range potential in the closure; each
!>   starts from t_s = T(c_s) of the previous stage, the first from c_s = 0
!> - Exact Jacobian and transpose products, J v = v - T'(c' v), and the
!>   Gaussian coarse space of the deflated Newton-Krylov steps
!> - Small-k damping of the updates the solver allows to post-process
!>
!> Adapted from moist_dev/pysofft (44874e7e), model/rism/1d/vv/driver.f90 and
!> model/rism/1d/deflation.f90
module moist_model_moz_solvent_vv_driver
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: resolve_num_threads
   use mctc_io_utils, only: to_lower
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_packed_sym, only: npair_from_ns
   use moist_model_moz_closure, only: moist_model_moz_closure_type, new_closure
   use moist_model_moz_solver, only: moist_model_moz_solver_type, moist_model_moz_solver_evaluator_type, &
      & moist_model_moz_solver_iteration_type, moist_model_moz_solver_printer_type, &
      & moist_model_moz_solver_preconditioner_type, moist_model_moz_solver_deflation_type, &
      & moist_model_moz_solver_parameters_type, moist_model_moz_solver_statistics_type, new_solver, &
      & run_iterations
   use moist_model_moz_1d_preconditioner, only: moist_model_moz_1d_preconditioner_spectral_type, &
      & moist_model_moz_1d_preconditioner_kspace_type, select_1d_preconditioner
   use moist_model_moz_solvent_vv_equation, only: vv_equation_type
   use moist_model_moz_solvent_vv_parameters, only: solvent_vv_parameters_type
   use moist_model_moz_fields, only: state_feasible, scale_flat, subtract_flat, unflatten, residual_difference
   implicit none(type, external)
   private

   public :: vv_problem_type, vv_evaluator_type, new_vv_problem

   !> Deflation profiles are cut where all of them fall below this fraction of their peak
   real(wp), parameter :: support_cutoff = 1.0e-16_wp

   !> Fixed VV data, operators and the iterate of one solve
   type :: vv_problem_type
      !> Equation with the per-k tables and transforms
      type(vv_equation_type) :: equation
      !> Closure and nonlinear solver
      class(moist_model_moz_closure_type), allocatable :: closure
      class(moist_model_moz_solver_type), allocatable :: solver
      !> Uniform DST-IV radial pair
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      !> Short-range potential, Hartree, w(k) and long-range potential, Hartree bohr^3, packed (npts, npair)
      real(wp), allocatable :: u_sr(:, :), wk(:, :), uk_lr(:, :)
      !> Site number densities, 1/bohr^3 (ns)
      real(wp), allocatable :: rho(:)
      !> Inverse thermal energy, 1/Hartree, and coupling of the current stage
      real(wp) :: beta = 0.0_wp, lambda = 1.0_wp
      !> Accepted iterate t_s(r) and its c_s(r), packed (npts, npair)
      real(wp), allocatable :: tr(:, :), cr(:, :)
      !> Closure derivative c'(t_s) of the last linearization and scratch, packed (npts, npair)
      real(wp), allocatable :: dcdt(:, :), c_work(:, :), dc(:, :), dt(:, :)
      !> Whether dcdt and the equation hold a linearization
      logical :: linearized = .false.
      !> Lowercase preconditioner name: none, spectral or kspace
      character(len=8) :: precond = "none"
      !> Small-k damping of post-processed updates
      logical :: damp_small_k = .false.
      !> Gaussian width, bohr, and profiles per site pair of the deflation space
      real(wp) :: deflation_sigma = 2.35_wp
      integer :: deflation_nfun = 2
      !> Spectral filter: preconditioner when selected and the Picard fallback filter
      type(moist_model_moz_1d_preconditioner_spectral_type) :: spectral
      !> Block k-space inverse, rebuilt by the outer iteration when selected
      type(moist_model_moz_1d_preconditioner_kspace_type) :: kspace
      !> Small-k damping filter
      type(moist_model_moz_1d_preconditioner_spectral_type) :: damping
   contains
      !> Run the coupling stages
      procedure :: solve => problem_solve
   end type vv_problem_type

   !> Residual and Jacobian of a borrowed problem; the target must outlive evaluator calls
   type, extends(moist_model_moz_solver_evaluator_type) :: vv_evaluator_type
      !> Problem whose workspaces the evaluator uses
      type(vv_problem_type), pointer :: problem => null()
   contains
      procedure :: residual => evaluator_residual
      procedure :: has_jvp => evaluator_has_jvp
      procedure :: linearize => evaluator_linearize
      procedure :: jvp => evaluator_jvp
      procedure :: has_adjoint => evaluator_has_adjoint
      procedure :: adjoint => evaluator_adjoint
      procedure :: deflation_setup => evaluator_deflation_setup
   end type vv_evaluator_type

   !> Solve-local outer iteration state
   type, extends(moist_model_moz_solver_iteration_type) :: vv_iteration_type
      !> Evaluator borrowing the problem
      type(vv_evaluator_type) :: evaluator
      !> Proposed iterate and its residual, packed (npts, npair), and damping scratch
      real(wp), allocatable :: gamma_new(:, :), residual(:, :), delta(:)
      !> Borrowed iteration table; null when quiet
      type(moist_model_moz_solver_printer_type), pointer :: printer => null()
      !> Stage tolerance; on the converging step the update is checked against its input
      real(wp) :: tolerance = 0.0_wp
   contains
      procedure :: advance => iteration_advance
      procedure :: observe => iteration_observe
      procedure :: checkpoint => iteration_checkpoint
   end type vv_iteration_type

contains

   !> Build a VV problem from the fixed tables and the settings
   !>
   !> @param[out] self      problem with c_s = 0
   !> @param[in]  rgrid     real-space uniform DST-IV grid
   !> @param[in]  kgrid     reciprocal grid of the pair
   !> @param[in]  wk        effective bulk w(k), packed (npts, npair)
   !> @param[in]  u_sr      short-range potential, Hartree, packed (npts, npair)
   !> @param[in]  uk_lr     long-range reciprocal potential, Hartree bohr^3, packed (npts, npair)
   !> @param[in]  rho       site number densities, 1/bohr^3 (ns)
   !> @param[in]  beta      inverse thermal energy, 1/Hartree
   !> @param[in]  params    validated settings
   !> @param[in]  nthreads  OpenMP team size of the closure, equation, solver and kspace inverse
   !> @param[in]  verbosity solver verbosity
   !> @param[out] error     closure, solver, equation or filter setup failure
   !> @param[in]  h_offset  optional fixed additive h(k) correction
   subroutine new_vv_problem(self, rgrid, kgrid, wk, u_sr, uk_lr, rho, beta, params, nthreads, verbosity, &
      & error, h_offset)
      !> Problem
      type(vv_problem_type), intent(out) :: self
      !> Radial pair
      type(moist_math_grid_radial_type), intent(in) :: rgrid, kgrid
      !> Packed w(k), short-range and long-range potentials
      real(wp), intent(in) :: wk(:, :), u_sr(:, :), uk_lr(:, :)
      !> Site number densities
      real(wp), intent(in) :: rho(:)
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Settings
      type(solvent_vv_parameters_type), intent(in) :: params
      !> Team size and solver verbosity
      integer, intent(in) :: nthreads, verbosity
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: h_offset(:, :)
      type(moist_model_moz_solver_parameters_type) :: controls

      if (any(shape(u_sr) /= shape(wk)) .or. size(wk, 2) /= npair_from_ns(size(rho))) then
         call fatal_error(error, "1D VV potential tables must have shape (npts, ns*(ns+1)/2)")
         return
      end if
      call new_closure(self%closure, params%closure, error, nthreads)
      if (allocated(error)) return
      controls = params%solver
      controls%nthreads = nthreads
      controls%verbosity = verbosity
      call new_solver(self%solver, params%solver_name, error, controls)
      if (allocated(error)) return
      self%equation%nthreads = nthreads
      call self%equation%setup(rgrid, kgrid, wk, rho, beta, uk_lr, error, h_offset)
      if (allocated(error)) return
      self%spectral%nthreads = nthreads
      call self%spectral%setup(rgrid, kgrid, params%precond_k0, params%precond_w0, error)
      if (allocated(error)) return
      self%damp_small_k = params%damp_small_k
      if (self%damp_small_k) then
         self%damping%nthreads = nthreads
         call self%damping%setup(rgrid, kgrid, params%damp_small_k_k0, params%damp_small_k_min_weight, error)
         if (allocated(error)) return
      end if
      self%kspace%nthreads = nthreads
      self%precond = to_lower(trim(adjustl(params%precond)))
      self%deflation_sigma = params%deflation_sigma
      self%deflation_nfun = params%deflation_nfun
      self%rgrid = rgrid
      self%kgrid = kgrid
      self%wk = wk
      self%u_sr = u_sr
      self%uk_lr = uk_lr
      self%rho = rho
      self%beta = beta
      allocate (self%tr, self%dcdt, self%c_work, self%dc, self%dt, mold=wk)
      allocate (self%cr(size(wk, 1), size(wk, 2)), source=0.0_wp)
   end subroutine new_vv_problem

   !> Run the coupling stages lambda = 1/n, ..., 1 from the current c_s
   !>
   !> - Each stage starts from t_s = T(c_s) and a fresh solver stage
   !> - A converged stage ends on a field whose residual was measured below the
   !>   tolerance, with cr = C(tr); rmse is the residual of that field
   !> - Stops at the first stage that does not converge; tr and cr hold its last iterate
   !>
   !> @param[inout] self       problem; returns the last accepted t_s and c_s
   !> @param[in]    tolerance  RMS residual of a converged stage
   !> @param[in]    maxiter    outer iterations per stage
   !> @param[in]    nlambda    number of stages
   !> @param[out]   converged  whether every stage converged
   !> @param[out]   rmse       last RMS residual
   !> @param[out]   niter      outer iterations of all stages
   !> @param[out]   error      map, residual, closure or preconditioner failure
   !> @param[in]    unit       output unit of the iteration table; absent for a quiet solve
   subroutine problem_solve(self, tolerance, maxiter, nlambda, converged, rmse, niter, error, unit)
      !> Problem
      class(vv_problem_type), intent(inout), target :: self
      !> Stage tolerance
      real(wp), intent(in) :: tolerance
      !> Iteration limit and stage count
      integer, intent(in) :: maxiter, nlambda
      !> Convergence flag
      logical, intent(out) :: converged
      !> Last RMS residual
      real(wp), intent(out) :: rmse
      !> Outer iterations
      integer, intent(out) :: niter
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Output unit of the iteration table
      integer, intent(in), optional :: unit
      type(vv_iteration_type) :: context
      type(moist_model_moz_solver_printer_type), target :: printer
      type(moist_model_moz_solver_statistics_type) :: work0, work1
      class(moist_model_moz_solver_preconditioner_type), pointer :: precond
      character(len=32) :: buf
      integer :: ilam, stage_iter

      converged = .false.
      rmse = huge(1.0_wp)
      niter = 0
      context%evaluator%problem => self
      context%tolerance = tolerance
      allocate (context%gamma_new, context%residual, mold=self%tr)
      if (self%damp_small_k) allocate (context%delta(size(self%tr)))
      if (present(unit)) then
         call printer%setup(unit=unit)
         call printer%header("1D VV")
         call printer%blank()
         context%printer => printer
         call self%solver%set_iteration_printer(printer)
      end if
      if (self%solver%uses_preconditioner()) then
         call select_1d_preconditioner(self%precond, self%spectral, precond, self%kspace)
         if (associated(precond)) call self%solver%set_preconditioner(precond)
         call self%solver%set_fallback_filter(self%spectral)
      end if
      work0 = self%solver%get_statistics()

      do ilam = 1, nlambda
         self%lambda = real(ilam, wp)/real(nlambda, wp)
         converged = .false.
         call self%equation%map(self%cr, self%tr, error)
         if (allocated(error)) exit
         self%linearized = .false.
         call self%solver%begin_stage()
         if (associated(context%printer)) then
            if (ilam > 1) then
               call printer%separator()
               call printer%blank()
            end if
            write (buf, "(f0.3)") self%lambda
            call printer%note("Solver: "//self%solver%name()//", lambda: "//trim(buf))
            call printer%blank()
            call printer%print_header()
            call printer%separator()
         end if
         call run_iterations(context, maxiter, tolerance, 0, rmse, stage_iter, converged, error)
         niter = niter + stage_iter
         if (allocated(error) .or. .not. converged) exit
      end do

      call self%solver%set_iteration_printer()
      call self%solver%set_preconditioner()
      call self%solver%set_fallback_filter()
      if (associated(context%printer)) then
         call printer%separator()
         call printer%blank()
         work1 = self%solver%get_statistics()
         if (work1%nresidual >= 0) then
            write (unit, "(2x,a,'...',i10)") "Krylov iterations         ", work1%nkrylov - work0%nkrylov
            write (unit, "(2x,a,'...',i10)") "Residual evaluations      ", work1%nresidual - work0%nresidual
            write (unit, "(2x,a,'...',i10)") "Exact Jacobian products   ", work1%njvp - work0%njvp
            write (unit, "(a)") ""
         end if
      end if
      converged = converged .and. .not. allocated(error)
   end subroutine problem_solve

   !> Take one solver update, damp it if allowed, and refresh c_s and the kspace inverse
   !>
   !> When the input iterate meets the stage tolerance the stage ends on a
   !> checked field: see finish_stage
   !>
   !> @param[inout] self     solve-local context
   !> @param[in]    iter     stage-local iteration
   !> @param[out]   rmse     RMS residual of the input iterate
   !> @param[out]   feasible whether an update was accepted
   !> @param[out]   error    solver, closure or preconditioner failure
   subroutine iteration_advance(self, iter, rmse, feasible, error)
      !> Solve-local context
      class(vv_iteration_type), intent(inout) :: self
      !> Stage-local iteration
      integer, intent(in) :: iter
      !> RMS residual
      real(wp), intent(out) :: rmse
      !> Whether the update was accepted
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      associate (p => self%evaluator%problem)
         call p%solver%outer_step(self%evaluator, p%tr, self%gamma_new, rmse, feasible, error)
         if (allocated(error) .or. .not. feasible) return
         if (p%damp_small_k .and. p%solver%allows_postprocessing()) then
            self%delta = reshape(self%gamma_new - p%tr, [size(self%delta)])
            call p%damping%apply(self%delta)
            self%gamma_new = p%tr + reshape(self%delta, shape(p%tr))
         end if
         if (rmse < self%tolerance) then
            call finish_stage(self, rmse, feasible, error)
            return
         end if
         p%tr = self%gamma_new
         call p%closure%apply(p%u_sr, p%tr, p%beta, p%lambda, p%cr, error)
         if (allocated(error)) then
            feasible = .false.
            return
         end if
         if (p%precond == "kspace" .and. p%solver%uses_preconditioner()) then
            if (p%kspace%wants_update(rmse, iter)) then
               call p%closure%derivative(p%u_sr, p%tr, p%beta, p%lambda, p%dc, error)
               if (.not. allocated(error)) call p%kspace%update(p%wk, p%cr, p%rho, p%beta, p%uk_lr, p%dc, error)
               if (allocated(error)) feasible = .false.
            end if
         end if
      end associate
   end subroutine iteration_advance

   !> End a converged stage on the better of the input iterate and the update
   !>
   !> The solver measured the input; the update gets a fresh residual and is
   !> kept only when it is not worse. rmse becomes the residual of the kept
   !> field and cr = C(tr), so the stage returns the field it reports
   !>
   !> @param[inout] self     solve-local context with the update in gamma_new
   !> @param[inout] rmse     residual of the input iterate; of the kept field on return
   !> @param[out]   feasible whether c_s of the kept field is usable
   !> @param[out]   error    closure failure on the kept field
   subroutine finish_stage(self, rmse, feasible, error)
      !> Solve-local context
      class(vv_iteration_type), intent(inout) :: self
      !> RMS residual
      real(wp), intent(inout) :: rmse
      !> Whether the kept field is usable
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(error_type), allocatable :: trial_error
      real(wp) :: fresh
      logical :: usable

      associate (p => self%evaluator%problem)
         call self%evaluator%residual(self%gamma_new, self%residual, usable, trial_error)
         if (usable .and. .not. allocated(trial_error)) then
            fresh = sqrt(sum(self%residual**2)/real(size(self%residual), wp))
            if (fresh <= rmse) then
               p%tr = self%gamma_new
               rmse = fresh
            end if
         end if
         call p%closure%apply(p%u_sr, p%tr, p%beta, p%lambda, p%cr, error)
         feasible = .not. allocated(error)
      end associate
   end subroutine finish_stage

   !> Print one row of the iteration table
   !>
   !> @param[inout] self    solve-local context
   !> @param[in]    iter    stage-local iteration
   !> @param[in]    rmse    RMS residual
   !> @param[in]    elapsed advance time, s
   subroutine iteration_observe(self, iter, rmse, elapsed)
      !> Solve-local context
      class(vv_iteration_type), intent(inout) :: self
      !> Stage-local iteration
      integer, intent(in) :: iter
      !> RMS residual and advance time
      real(wp), intent(in) :: rmse, elapsed
      type(moist_model_moz_solver_statistics_type) :: stats

      if (.not. associated(self%printer)) return
      associate (p => self%evaluator%problem)
         stats = p%solver%get_statistics()
         call self%printer%print_outer_row(iter, elapsed, rmse, p%solver%prev_step_maxabs, &
            & hmax=maxval(p%cr + p%tr), lambda=p%lambda, alpha=p%solver%prev_alpha, ndiis=stats%ndiis, &
            & gmres_total=stats%ngmres)
      end associate
   end subroutine iteration_observe

   !> No file checkpoints; converged bulk fields are retained in memory
   !>
   !> @param[inout] self  solve-local context
   !> @param[out]   error always set
   subroutine iteration_checkpoint(self, error)
      !> Solve-local context
      class(vv_iteration_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      associate (s => self); end associate
      call fatal_error(error, "1D VV solve keeps no checkpoints")
   end subroutine iteration_checkpoint

   !> Residual R(t_s) = t_s - T(C(t_s))
   !>
   !> @param[in]  self     evaluator bound to a problem
   !> @param[in]  gamma    iterate t_s(r), packed (npts, npair)
   !> @param[out] residual residual, packed (npts, npair)
   !> @param[out] feasible whether the residual is usable
   !> @param[out] error    closure failure, 1 + h below zero, map failure or non-finite residual
   subroutine evaluator_residual(self, gamma, residual, feasible, error)
      !> Evaluator
      class(vv_evaluator_type), intent(in) :: self
      !> Iterate
      real(wp), intent(in) :: gamma(:, :)
      !> Residual
      real(wp), intent(out) :: residual(:, :)
      !> Whether the residual is usable
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      feasible = .false.
      if (.not. associated(self%problem)) then
         call fatal_error(error, "1D VV evaluator is not bound to a problem")
         return
      end if
      associate (p => self%problem)
         call p%closure%apply(p%u_sr, gamma, p%beta, p%lambda, p%c_work, error)
         if (allocated(error)) return
         if (.not. state_feasible(p%c_work, gamma, resolve_num_threads(p%equation%nthreads))) then
            call fatal_error(error, "1D VV closure gave 1 + h < 0")
            return
         end if
         call p%equation%map(p%c_work, residual, error)
         if (allocated(error)) return
         call residual_difference(gamma, residual, resolve_num_threads(p%equation%nthreads), feasible)
         if (.not. feasible) call fatal_error(error, "1D VV residual is not finite")
      end associate
   end subroutine evaluator_residual

   !> Exact Jacobian-vector products are available once bound
   !>
   !> @param[in] self Evaluator
   logical function evaluator_has_jvp(self) result(available)
      !> Evaluator
      class(vv_evaluator_type), intent(in) :: self
      available = associated(self%problem)
   end function evaluator_has_jvp

   !> Cache c'(t_s) and the per-k linearization of the map at t_s
   !>
   !> @param[in]  self     evaluator bound to a problem
   !> @param[in]  gamma    linearization point t_s(r), packed (npts, npair)
   !> @param[out] feasible whether a linearization was stored
   !> @param[out] error    closure failure, 1 + h below zero or singular point solve
   subroutine evaluator_linearize(self, gamma, feasible, error)
      !> Evaluator
      class(vv_evaluator_type), intent(in) :: self
      !> Linearization point
      real(wp), intent(in) :: gamma(:, :)
      !> Whether the linearization is usable
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      feasible = .false.
      if (.not. associated(self%problem)) then
         call fatal_error(error, "1D VV evaluator is not bound to a problem")
         return
      end if
      associate (p => self%problem)
         p%linearized = .false.
         call p%closure%apply_and_derivative(p%u_sr, gamma, p%beta, p%lambda, p%c_work, p%dcdt, error)
         if (allocated(error)) return
         if (.not. state_feasible(p%c_work, gamma, resolve_num_threads(p%equation%nthreads))) then
            call fatal_error(error, "1D VV closure gave 1 + h < 0 at the linearization point")
            return
         end if
         call p%equation%linearize(p%c_work, error)
         feasible = .not. allocated(error)
         p%linearized = feasible
      end associate
   end subroutine evaluator_linearize

   !> Exact product J v = v - T'(c' v) at the last linearization
   !>
   !> @param[in]  self      evaluator with a linearization
   !> @param[in]  direction v, flat (npts*npair)
   !> @param[out] jv        J v, flat (npts*npair)
   !> @param[out] feasible  whether the product is usable
   !> @param[out] error     no linearization, transform failure or non-finite product
   subroutine evaluator_jvp(self, direction, jv, feasible, error)
      !> Evaluator
      class(vv_evaluator_type), intent(in) :: self
      !> Direction
      real(wp), intent(in) :: direction(:)
      !> Product
      real(wp), intent(out) :: jv(:)
      !> Whether the product is usable
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      feasible = .false.
      if (.not. linearized(self, size(direction), size(jv), error)) then
         jv = 0.0_wp
         return
      end if
      associate (p => self%problem)
         call scale_flat(p%dcdt, direction, p%dc, resolve_num_threads(p%equation%nthreads))
         call p%equation%jvp(p%dc, p%dt, error)
         if (allocated(error)) then
            jv = 0.0_wp
            return
         end if
         call subtract_flat(direction, p%dt, jv, resolve_num_threads(p%equation%nthreads), feasible)
         if (.not. feasible) call fatal_error(error, "1D VV Jacobian product is not finite")
      end associate
   end subroutine evaluator_jvp

   !> Transpose products are available once bound
   !>
   !> @param[in] self Evaluator
   logical function evaluator_has_adjoint(self) result(available)
      !> Evaluator
      class(vv_evaluator_type), intent(in) :: self
      available = associated(self%problem)
   end function evaluator_has_adjoint

   !> Exact transpose product J^T u = u - c' T'^T(u) at the last linearization
   !>
   !> @param[in]  self      evaluator with a linearization
   !> @param[in]  direction u, flat (npts*npair)
   !> @param[out] jtv       J^T u, flat (npts*npair)
   !> @param[out] feasible  whether the product is usable
   !> @param[out] error     no linearization, transform failure or non-finite product
   subroutine evaluator_adjoint(self, direction, jtv, feasible, error)
      !> Evaluator
      class(vv_evaluator_type), intent(in) :: self
      !> Direction
      real(wp), intent(in) :: direction(:)
      !> Transpose product
      real(wp), intent(out) :: jtv(:)
      !> Whether the product is usable
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      feasible = .false.
      if (.not. linearized(self, size(direction), size(jtv), error)) then
         jtv = 0.0_wp
         return
      end if
      associate (p => self%problem)
         call unflatten(direction, p%dt, resolve_num_threads(p%equation%nthreads))
         call p%equation%adjoint(p%dt, p%dc, error)
         if (allocated(error)) then
            jtv = 0.0_wp
            return
         end if
         call subtract_flat(direction, p%dc, jtv, resolve_num_threads(p%equation%nthreads), feasible, p%dcdt)
         if (.not. feasible) call fatal_error(error, "1D VV transpose product is not finite")
      end associate
   end subroutine evaluator_adjoint

   !> Fill and factorize the coarse matrix at the last linearization
   !>
   !> Builds the Gaussian profiles on first use; without a linearization, or on
   !> a coarse-matrix failure, the deflation stays inactive
   !>
   !> @param[in]    self      evaluator
   !> @param[inout] deflation basis, coarse matrix and its factors
   subroutine evaluator_deflation_setup(self, deflation)
      !> Evaluator
      class(vv_evaluator_type), intent(in) :: self
      !> Deflation object
      type(moist_model_moz_solver_deflation_type), intent(inout) :: deflation
      type(error_type), allocatable :: error

      call deflation%deactivate()
      if (.not. associated(self%problem)) return
      associate (p => self%problem)
         if (.not. p%linearized) return
         if (.not. deflation%has_basis(size(p%tr, 1), size(p%tr, 2))) then
            call build_profiles(p%rgrid%r, p%deflation_sigma, p%deflation_nfun, size(p%tr, 2), deflation)
         end if
         call p%equation%coarse_matrix(deflation%profile, deflation%wprofile, p%dcdt, deflation%coarse, error)
         if (allocated(error)) return
         call deflation%factorize()
      end associate
   end subroutine evaluator_deflation_setup

   !> Whether the evaluator holds a linearization for flat vectors of this size
   !>
   !> @param[in]  self   evaluator
   !> @param[in]  nin    size of the direction
   !> @param[in]  nout   size of the product
   !> @param[out] error  unbound, not linearized or wrong size
   logical function linearized(self, nin, nout, error)
      !> Evaluator
      class(vv_evaluator_type), intent(in) :: self
      !> Vector sizes
      integer, intent(in) :: nin, nout
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      linearized = .false.
      if (.not. associated(self%problem)) then
         call fatal_error(error, "1D VV evaluator is not bound to a problem")
      else if (.not. self%problem%linearized) then
         call fatal_error(error, "1D VV Jacobian product requested without a linearization")
      else if (nin /= size(self%problem%tr) .or. nout /= nin) then
         call fatal_error(error, "1D VV Jacobian vectors must have npts*npair entries")
      else
         linearized = .true.
      end if
   end function linearized

   !> Gaussian radial profiles of the coarse space, one set per site pair
   !>
   !> - g_1 = exp(-r^2/(2 sigma^2)), g_j = (r/sigma)^(2(j-1)) g_1, normalized
   !> - Weighted with Omega = diag(r^2) and cut where all profiles have decayed
   !>
   !> @param[in]    r          radial nodes, bohr (npts)
   !> @param[in]    sigma      Gaussian width, bohr
   !> @param[in]    nfun       profiles per site pair
   !> @param[in]    npair      site pairs
   !> @param[inout] deflation  receives the basis
   subroutine build_profiles(r, sigma, nfun, npair, deflation)
      !> Radial nodes
      real(wp), intent(in) :: r(:)
      !> Gaussian width
      real(wp), intent(in) :: sigma
      !> Profiles per pair and site pairs
      integer, intent(in) :: nfun, npair
      !> Deflation object
      type(moist_model_moz_solver_deflation_type), intent(inout) :: deflation
      real(wp), allocatable :: profile(:, :), wprofile(:, :)
      real(wp) :: profile_peak(nfun)
      integer :: i, j, ncore

      allocate (profile(size(r), nfun), wprofile(size(r), nfun))
      do j = 1, nfun
         profile(:, j) = exp(-r**2/(2.0_wp*sigma**2))*(r/sigma)**(2*(j - 1))
         profile(:, j) = profile(:, j)/norm2(profile(:, j))
      end do
      ! Profiles remain fixed throughout the backward support scan.
      profile_peak = maxval(abs(profile), dim=1)
      ncore = 0
      do i = size(r), 1, -1
         if (any(abs(profile(i, :)) > support_cutoff*profile_peak)) then
            ncore = i
            exit
         end if
      end do
      profile(ncore + 1:, :) = 0.0_wp
      do j = 1, nfun
         wprofile(:, j) = r**2*profile(:, j)
      end do
      call deflation%init_basis(profile, wprofile, npair)
   end subroutine build_profiles

end module moist_model_moz_solvent_vv_driver
