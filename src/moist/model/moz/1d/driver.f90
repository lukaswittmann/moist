!> 1D UV solve: residual evaluator, outer iterations and solvation free energy
!>
!> - Iterate is the short-range t_s(r)
!> - Closure gives c_s(r), equation maps it back, R(t_s) = t_s - T(C(t_s))
!> - One stage at full coupling, from a supplied t_s guess or c_s = 0
!> - Exact Jacobian and transpose products, J v = v - T'(c' v), T affine
!> - Small-k damping of updates the solver lets us post-process
!> - Solvation free energy dG = (1/beta) sum_ia rho_a int Phi(h_ia, c_ia) d^3r
!>   with physical h and c
!>   - Coulomb tails of c cancel in the sum over the sites of a neutral solvent
!> - Direct-correlation moment D = sum_ia rho_a int c_ia d^3r of the physical c,
!>   for the partial molar volume
!>   - Same Coulomb tails cancel in it
!> - Response: derivatives of dG with respect to u_sr, u_lr(r) and u_lr(k)
!>   at a converged iterate, from one transposed solve J^T lambda = dE/dt_s
module moist_model_moz_1d_driver
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use mctc_io_convert, only: autokcal
   use mctc_io_utils, only: to_lower
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_packed_sym, only: npair_from_ns
   use moist_model_moz_closure, only: moist_model_moz_closure_type, new_closure
   use moist_model_moz_solver, only: moist_model_moz_solver_type, moist_model_moz_solver_evaluator_type, &
      & moist_model_moz_solver_jfnk_adjoint_parameters_type, &
      & moist_model_moz_solver_iteration_type, moist_model_moz_solver_printer_type, &
      & moist_model_moz_solver_preconditioner_type, &
      & moist_model_moz_solver_statistics_type, new_solver, run_iterations, solve_adjoint, &
      & check_preconditioner_name
   use moist_model_moz_1d_preconditioner, only: moist_model_moz_1d_preconditioner_spectral_type, &
      & select_1d_preconditioner, preconditioners_1d_uv
   use moist_model_moz_1d_equation, only: uv_equation_type
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_response_kernel, only: uv_response_rhs, uv_response_sr, uv_response_lr
   use moist_model_moz_thermodynamics_kernel, only: scale_energy_derivative
   use moist_model_moz_fields, only: state_feasible, scale_flat, subtract_flat, unflatten, residual_difference
   implicit none(type, external)
   private

   public :: uv_problem_type, uv_evaluator_type, new_uv_problem


   !> Fixed UV data, operators and the iterate of one solve
   type :: uv_problem_type
      !> Run context, borrowed for the lifetime of the problem
      type(moist_context_type), pointer :: ctx => null()
      !> Equation with the per-k data and transforms
      type(uv_equation_type) :: equation
      !> Closure
      class(moist_model_moz_closure_type), allocatable :: closure
      !> Nonlinear solver
      class(moist_model_moz_solver_type), allocatable :: solver
      !> Transpose Krylov controls, independent of the nonlinear solver
      type(moist_model_moz_solver_jfnk_adjoint_parameters_type) :: adjoint
      !> Uniform DST-IV radial pair
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      !> Short- and long-range real-space potential, Hartree, (npts, nu nv)
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :)
      !> Solvent site number densities, 1/bohr^3 (nv)
      real(wp), allocatable :: rho(:)
      !> Inverse thermal energy, 1/Hartree
      real(wp) :: beta = 0.0_wp
      !> Accepted iterate t_s(r) and its c_s(r), (npts, nu nv)
      real(wp), allocatable :: tr(:, :), cr(:, :)
      !> Closure derivative c'(t_s) of the last linearization and scratch, (npts, nu nv)
      real(wp), allocatable :: dcdt(:, :), c_work(:, :), dc(:, :), dt(:, :)
      !> Whether dcdt holds a linearization
      logical :: linearized = .false.
      !> Whether the current t_s can seed another solve on the same radial pair
      logical, private :: has_guess = .false.
      !> Lowercase preconditioner name: none or spectral
      character(len=8) :: precond = "none"
      !> Small-k damping of post-processed updates
      logical :: damp_small_k = .false.
      !> Spectral filter: preconditioner when selected and the Picard fallback filter
      type(moist_model_moz_1d_preconditioner_spectral_type) :: spectral
      !> Small-k damping filter
      type(moist_model_moz_1d_preconditioner_spectral_type) :: damping
   contains
      !> Replace fixed-grid potential tables while retaining operators and solver storage
      procedure :: update_potential => problem_update_potential
      !> Solve the UV equations from an optional t_s guess, otherwise c_s = 0
      procedure :: solve => problem_solve
      !> Solvation free energy of the current iterate
      procedure :: free_energy => problem_free_energy
      !> Direct-correlation moment of the current iterate
      procedure :: direct_moment => problem_direct_moment
      !> Table adjoints of the free energy at a converged iterate
      procedure :: response => problem_response
   end type uv_problem_type

   !> Residual and Jacobian of a borrowed problem; the target must outlive evaluator calls
   type, extends(moist_model_moz_solver_evaluator_type) :: uv_evaluator_type
      !> Problem whose workspaces the evaluator uses
      type(uv_problem_type), pointer :: problem => null()
   contains
      procedure :: residual => evaluator_residual
      procedure :: has_jvp => evaluator_has_jvp
      procedure :: linearize => evaluator_linearize
      procedure :: jvp => evaluator_jvp
      procedure :: has_adjoint => evaluator_has_adjoint
      procedure :: adjoint => evaluator_adjoint
   end type uv_evaluator_type

   !> Solve-local outer iteration state
   type, extends(moist_model_moz_solver_iteration_type) :: uv_iteration_type
      !> Evaluator borrowing the problem
      type(uv_evaluator_type) :: evaluator
      !> Proposed iterate and its residual, (npts, nu nv), and damping scratch
      real(wp), allocatable :: gamma_new(:, :), residual(:, :), delta(:)
      !> Borrowed iteration table; null when quiet
      type(moist_model_moz_solver_printer_type), pointer :: printer => null()
      !> Tolerance; on the converging step the update is checked against its input
      real(wp) :: tolerance = 0.0_wp
   contains
      procedure :: advance => iteration_advance
      procedure :: observe => iteration_observe
      procedure :: checkpoint => iteration_checkpoint
   end type uv_iteration_type

contains

   !> Build a UV problem from the fixed tables and the settings
   !>
   !> @param[out] self       prepared problem; solve initializes the cold or retained guess
   !> @param[in]  ctx        run context; verbosity, output unit and OpenMP team size
   !> @param[in]  rgrid      real-space uniform DST-IV grid
   !> @param[in]  kgrid      reciprocal grid of the pair
   !> @param[in]  wu         solute W(k), packed (npts, nu(nu+1)/2)
   !> @param[in]  response   solvent response P(k) on kgrid, (npts, nv, nv)
   !> @param[in]  u_sr       short-range potential, Hartree, (npts, nu nv)
   !> @param[in]  ur_lr      long-range real-space potential, Hartree, (npts, nu nv)
   !> @param[in]  uk_lr      long-range reciprocal potential, Hartree bohr^3, (npts, nu nv)
   !> @param[in]  rho        solvent site number densities, 1/bohr^3 (nv)
   !> @param[in]  beta       inverse thermal energy, 1/Hartree
   !> @param[in]  params     validated settings
   !> @param[out] error      closure, solver, equation or filter setup failure
   subroutine new_uv_problem(self, ctx, rgrid, kgrid, wu, response, u_sr, ur_lr, uk_lr, rho, beta, params, error)
      !> Problem
      type(uv_problem_type), intent(out) :: self
      !> Run context, borrowed for the lifetime of the problem
      type(moist_context_type), intent(in), target :: ctx
      !> Radial pair
      type(moist_math_grid_radial_type), intent(in) :: rgrid, kgrid
      !> Solute W(k) and solvent response
      real(wp), intent(in) :: wu(:, :), response(:, :, :)
      !> Potential tables
      real(wp), intent(in) :: u_sr(:, :), ur_lr(:, :), uk_lr(:, :)
      !> Solvent site densities
      real(wp), intent(in) :: rho(:)
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Settings
      type(model_moz_1d_parameters_type), intent(in) :: params
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      self%ctx => ctx
      if (size(rho) < 1 .or. size(rho) /= size(response, 2)) then
         call fatal_error(error, "1D UV needs one density per solvent site in the response")
         return
      end if
      if (any(shape(u_sr) /= shape(uk_lr)) .or. any(shape(ur_lr) /= shape(uk_lr)) .or. &
         & mod(size(uk_lr, 2), size(rho)) /= 0) then
         call fatal_error(error, "1D UV potential tables must have shape (npts, nu nv)")
         return
      end if
      call new_closure(self%closure, params%closure, error, self%ctx)
      if (allocated(error)) return
      call check_preconditioner_name(params%precond, preconditioners_1d_uv, "1D UV", error)
      if (allocated(error)) return
      self%adjoint = params%solver%jfnk%adjoint
      call new_solver(self%solver, ctx, error, params%solver)
      if (allocated(error)) return
      self%equation%nthreads = ctx%get_num_threads()
      call self%equation%setup(rgrid, kgrid, wu, response, beta, uk_lr, error)
      if (allocated(error)) return
      self%spectral%nthreads = ctx%get_num_threads()
      call self%spectral%setup(rgrid, kgrid, params%precond_k0, params%precond_w0, error)
      if (allocated(error)) return
      self%damp_small_k = params%damp_small_k
      if (self%damp_small_k) then
         self%damping%nthreads = ctx%get_num_threads()
         call self%damping%setup(rgrid, kgrid, params%damp_small_k_k0, params%damp_small_k_min_weight, error)
         if (allocated(error)) return
      end if
      self%precond = to_lower(trim(adjustl(params%precond)))
      self%rgrid = rgrid
      self%kgrid = kgrid
      self%u_sr = u_sr
      self%ur_lr = ur_lr
      self%rho = rho
      self%beta = beta
      allocate (self%tr, self%dcdt, self%c_work, self%dc, self%dt, mold=u_sr)
      allocate (self%cr(size(u_sr, 1), size(u_sr, 2)), source=0.0_wp)
   end subroutine new_uv_problem

   !> Update host-dependent tables at fixed geometry, bulk, beta and solver settings
   !>
   !> - Preserves physical t: t_s(new) = t_s(old) + beta*(u_lr(old) - u_lr(new))
   !> - Next solve checks the retained field against the new equation and closure
   !>
   !> @param[inout] self    constructed problem; resources retained
   !> @param[in]    u_sr    short-range potential, Hartree, (npts, nu nv)
   !> @param[in]    ur_lr   real-space long-range potential, Hartree, (npts, nu nv)
   !> @param[in]    uk_lr   reciprocal long-range potential, Hartree bohr^3, (npts, nu nv)
   !> @param[out]   error   unconstructed problem or wrong shape; old tables retained on refusal
   subroutine problem_update_potential(self, u_sr, ur_lr, uk_lr, error)
      !> Problem
      class(uv_problem_type), intent(inout) :: self
      !> Short-range, real-space long-range (Hartree) and reciprocal long-range (Hartree bohr^3) potentials, (npts, nu nv)
      real(wp), intent(in) :: u_sr(:, :), ur_lr(:, :), uk_lr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(self%u_sr)) then
         call fatal_error(error, "Construct the 1D UV problem before replacing its potential")
         return
      end if
      if (any(shape(u_sr) /= shape(self%u_sr)) .or. any(shape(ur_lr) /= shape(self%ur_lr))) then
         call fatal_error(error, "1D UV potential update must keep the field shape")
         return
      end if
      call self%equation%update_potential(uk_lr, error)
      if (allocated(error)) return
      if (self%has_guess) self%tr = self%tr + self%beta*(self%ur_lr - ur_lr)
      self%u_sr(:, :) = u_sr
      self%ur_lr(:, :) = ur_lr
      self%linearized = .false.
   end subroutine problem_update_potential

   !> Solve the UV equations at full coupling
   !>
   !> - Starts from a usable t_s guess or t_s = T(0), with a fresh solver stage
   !> - A guess is checked against the current closure and equation; an incompatible,
   !>   non-finite or infeasible guess falls back to c_s = 0 with a warning; a failed
   !>   feasible warm attempt retries that cold path once, with fresh solver history,
   !>   and warns with the dropped error, if any
   !> - A converged solve ends on a field whose residual was measured below the
   !>   tolerance, with cr = C(tr); rmse is the residual of that field
   !> - Iteration table on the context unit at verbosity 4, with the solver's
   !>   fallback summary after each attempt
   !>
   !> @param[inout] self       problem; returns the last accepted t_s and c_s
   !> @param[in]    tolerance  RMS residual of a converged solve
   !> @param[in]    maxiter    outer iterations
   !> @param[out]   converged  whether the solve converged
   !> @param[out]   rmse       last RMS residual
   !> @param[out]   niter      outer iterations
   !> @param[out]   error      map, residual or closure failure
   !> @param[in]    guess      initial short-range t_s(r), (npts, nu nv); absent for a fresh start
   !> @param[in]    reuse_initial_guess Try the retained field when guess is absent; default false
   subroutine problem_solve(self, tolerance, maxiter, converged, rmse, niter, error, guess, reuse_initial_guess)
      !> Problem
      class(uv_problem_type), intent(inout), target :: self
      !> Tolerance
      real(wp), intent(in) :: tolerance
      !> Iteration limit
      integer, intent(in) :: maxiter
      !> Convergence flag
      logical, intent(out) :: converged
      !> Last RMS residual
      real(wp), intent(out) :: rmse
      !> Outer iterations
      integer, intent(out) :: niter
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Initial short-range indirect correlation
      real(wp), intent(in), optional :: guess(:, :)
      !> Try the previous converged field on the same radial pair
      logical, intent(in), optional :: reuse_initial_guess
      type(uv_iteration_type) :: context
      type(moist_model_moz_solver_printer_type), target :: printer
      type(moist_model_moz_solver_statistics_type) :: work0
      class(moist_model_moz_solver_preconditioner_type), pointer :: precond
      type(error_type), allocatable :: guess_error
      logical :: warm
      integer :: attempt, attempted_steps

      converged = .false.
      rmse = huge(1.0_wp)
      niter = 0
      context%evaluator%problem => self
      context%tolerance = tolerance
      allocate (context%gamma_new, context%residual, mold=self%tr)
      if (self%damp_small_k) allocate (context%delta(size(self%tr)))
      if (self%ctx%verbosity(4)) then
         call printer%setup(unit=self%ctx%unit, energy_label="DG kcal", denergy_label="dDG kcal", &
            & inner_rows=self%ctx%verbosity(5))
         context%printer => printer
         call self%solver%set_iteration_printer(printer)
         call printer%header("1D UV")
         call printer%blank()
      end if
      if (self%solver%uses_preconditioner()) then
         call select_1d_preconditioner(self%precond, self%spectral, precond)
         if (associated(precond)) call self%solver%set_preconditioner(precond)
         call self%solver%set_fallback_filter(self%spectral)
      end if
      work0 = self%solver%get_statistics()

      self%linearized = .false.
      warm = .false.
      if (present(guess)) then
         ! A rejected caller guess falls back to c_s = 0 and warns why
         if (any(shape(guess) /= shape(self%tr))) then
            call self%ctx%warning("1D UV: guess shape does not match the grid, ignoring guess")
         else if (.not. all(ieee_is_finite(guess))) then
            call self%ctx%warning("1D UV: guess is not finite, ignoring guess")
         else
            call context%evaluator%residual(guess, context%residual, warm, guess_error)
            warm = warm .and. .not. allocated(guess_error)
            if (warm) then
               self%tr = guess
               self%cr = self%c_work
            else if (allocated(guess_error)) then
               call self%ctx%warning("1D UV: guess is infeasible, ignoring guess: "//guess_error%message)
            else
               call self%ctx%warning("1D UV: guess is infeasible, ignoring guess")
            end if
         end if
      else if (present(reuse_initial_guess)) then
         if (reuse_initial_guess .and. self%has_guess) then
            call context%evaluator%residual(self%tr, context%residual, warm, guess_error)
            warm = warm .and. .not. allocated(guess_error)
            if (warm) self%cr = self%c_work
         end if
      end if
      self%has_guess = .false.
      do attempt = 1, 2
         if (.not. warm) then
            self%cr = 0.0_wp
            call self%equation%map(self%cr, self%tr, error)
            if (allocated(error)) exit
         end if
         call self%solver%begin_stage()
         if (associated(context%printer)) then
            if (attempt > 1) then
               call printer%separator()
               call printer%blank()
               call printer%note("Retry from c_s = 0")
            end if
            call printer%note("Solver: "//self%solver%name())
            call printer%blank()
            call printer%begin_stage()
         end if
         call run_iterations(context, self%ctx, maxiter, tolerance, 0, rmse, attempted_steps, &
            & converged, error)
         niter = niter + attempted_steps
         if (associated(context%printer)) call printer%note(self%solver%stage_summary())
         if (converged .and. .not. allocated(error)) exit
         if (.not. warm) exit
         if (allocated(error)) then
            call self%ctx%warning("1D UV: warm start failed: "//error%message//"; retrying cold")
            deallocate (error)
         end if
         warm = .false.
         self%linearized = .false.
      end do

      call self%solver%set_iteration_printer()
      call self%solver%set_preconditioner()
      call self%solver%set_fallback_filter()
      if (associated(context%printer)) then
         call printer%separator()
         call printer%blank()
         call self%solver%print_statistics(printer, since=work0)
      end if
      converged = converged .and. .not. allocated(error)
      self%has_guess = converged
   end subroutine problem_solve

   !> Solvation free energy of the current iterate, Hartree
   !>
   !> dG = (1/beta) sum_ia rho_a int Phi(h_ia, c_ia) d^3r with the physical
   !> c = c_s - beta u_lr and h = c_s + t_s
   !>
   !> @param[in]  self    problem with an accepted iterate
   !> @param[out] energy  solvation free energy, Hartree
   !> @param[out] error   invalid correlations
   subroutine problem_free_energy(self, energy, error)
      !> Problem
      class(uv_problem_type), intent(in) :: self
      !> Solvation free energy
      real(wp), intent(out) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: phi(:, :)
      real(wp) :: integral
      integer :: col, nv

      energy = 0.0_wp
      nv = size(self%rho)
      allocate (phi, mold=self%cr)
      call self%closure%free_energy_integrands(self%cr + self%tr, self%cr - self%beta*self%ur_lr, phi, error)
      if (allocated(error)) return
      do col = 1, size(phi, 2)
         call self%rgrid%integrate_field(phi(:, col), integral, error)
         if (allocated(error)) return
         energy = energy + self%rho(mod(col - 1, nv) + 1)*integral
      end do
      energy = energy/self%beta
   end subroutine problem_free_energy

   !> Direct-correlation moment of the current iterate
   !>
   !> - D = sum_ia rho_a int c_ia d^3r
   !> - Physical c = c_s - beta u_lr
   !>
   !> @param[in]  self    problem with an accepted iterate
   !> @param[out] moment  D, dimensionless
   !> @param[out] error   integration failure
   subroutine problem_direct_moment(self, moment, error)
      !> Problem
      class(uv_problem_type), intent(in) :: self
      !> Moment
      real(wp), intent(out) :: moment
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: integral
      integer :: col, nv

      moment = 0.0_wp
      nv = size(self%rho)
      do col = 1, size(self%cr, 2)
         call self%rgrid%integrate_field(self%cr(:, col) - self%beta*self%ur_lr(:, col), integral, error)
         if (allocated(error)) return
         moment = moment + self%rho(mod(col - 1, nv) + 1)*integral
      end do
   end subroutine problem_direct_moment

   !> Derivatives of the weighted solvation free energy with respect to the potential tables
   !>
   !> - Valid at a converged iterate; the closure derivative is taken there
   !> - Right-hand side g = dE/dt_s from energy-weighted Phi derivatives at physical h and c
   !> - Transposed solve J^T lambda = g with the problem's preconditioner
   !>   - Starts from lambda on entry when it has the field shape, else from zero
   !> - With dh = energy_seed*dPhi/dh and dc = energy_seed*dPhi/dc:
   !>   w_sr = dc_s/du_sr [a (dh + dc) + mu], mu = T'^T lambda;
   !>   w_lr = -beta a dc; w_k = -beta W^T k2r^T(lambda) P^T,
   !>   with a = rho_a 4 pi r^2 dr / beta
   !> - Adjoint convention: dE = sum w_sr du_sr + sum w_lr du_lr(r) + sum w_k du_lr(k)
   !>
   !> @param[inout] self       problem with a converged iterate
   !> @param[in]    tolerance  RMS of the transpose residual relative to the RMS of g
   !> @param[out]   w_sr       adjoint of u_sr, (npts, nu nv)
   !> @param[out]   w_lr       adjoint of u_lr(r), (npts, nu nv)
   !> @param[out]   w_k        adjoint of u_lr(k), (npts, nu nv)
   !> @param[inout] lambda     transposed solution; a field-shaped value on entry is the start
   !> @param[out]   error      infeasible iterate, closure failure or unconverged transposed solve
   !> @param[in]    free_energy_weight  optional free-energy sensitivity; default one
   !> @param[out]   w_wu       optional adjoint of packed solute W(k), Hartree
   subroutine problem_response(self, tolerance, w_sr, w_lr, w_k, lambda, error, free_energy_weight, w_wu)
      !> Problem
      class(uv_problem_type), intent(inout), target :: self
      !> Relative tolerance
      real(wp), intent(in) :: tolerance
      !> Table adjoints
      real(wp), allocatable, intent(out) :: w_sr(:, :), w_lr(:, :), w_k(:, :)
      !> Transposed solution
      real(wp), allocatable, intent(inout) :: lambda(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Free-energy sensitivity; defaults to one
      real(wp), optional, intent(in) :: free_energy_weight
      !> Optional intramolecular sensitivities from the same transposed solution
      real(wp), allocatable, intent(out), optional :: w_wu(:, :)
      type(uv_evaluator_type), target :: evaluator
      class(moist_model_moz_solver_preconditioner_type), pointer :: precond
      real(wp), allocatable :: weight(:), dh(:, :), dc(:, :), dcdu(:, :), rhs(:, :), mu(:, :), solution(:, :)
      real(wp) :: tol, energy_seed
      !> Per-point measure, private to each OpenMP worker
      real(wp) :: a
      integer :: i, col, nv, nt
      logical :: feasible, warm

      if (.not. ieee_is_finite(tolerance) .or. tolerance <= 0.0_wp) then
         call fatal_error(error, "1D UV response needs a positive finite tolerance")
         return
      end if
      energy_seed = 1.0_wp
      if (present(free_energy_weight)) energy_seed = free_energy_weight
      if (.not. ieee_is_finite(energy_seed)) then
         call fatal_error(error, "1D UV response needs a finite free-energy sensitivity")
         return
      end if
      nt = self%ctx%get_num_threads()
      nv = size(self%rho)
      evaluator%problem => self
      call evaluator%linearize(self%tr, feasible, error)
      if (allocated(error)) return
      allocate (dh, dc, dcdu, rhs, mold=self%tr)
      call self%closure%free_energy_derivative_fields(self%c_work + self%tr, self%c_work - self%beta*self%ur_lr, &
         & dh, dc, error)
      if (allocated(error)) return
      call self%closure%potential_derivative(self%u_sr, self%tr, self%beta, 1.0_wp, dcdu, error)
      if (allocated(error)) return
      weight = 4.0_wp*acos(-1.0_wp)*self%rgrid%r**2*self%rgrid%w/self%beta
      !$omp parallel do collapse(2) schedule(static) num_threads(nt)
      do col = 1, size(rhs, 2)
         do i = 1, size(rhs, 1)
            if (energy_seed /= 1.0_wp) then
               dh(i, col) = scale_energy_derivative(energy_seed, dh(i, col))
               dc(i, col) = scale_energy_derivative(energy_seed, dc(i, col))
            end if
            rhs(i, col) = uv_response_rhs(self%rho(mod(col - 1, nv) + 1)*weight(i), dh(i, col), dc(i, col), &
               & self%dcdt(i, col))
         end do
      end do
      !$omp end parallel do

      tol = max(tolerance*sqrt(sum(rhs**2)/real(size(rhs), wp)), tiny(1.0_wp))
      warm = .false.
      if (allocated(lambda)) warm = all(shape(lambda) == shape(rhs))
      ! A zero objective has an exact zero solution; discard a nonzero warm start
      warm = warm .and. any(rhs /= 0.0_wp)
      call select_1d_preconditioner(self%precond, self%spectral, precond)
      if (associated(precond) .and. warm) then
         call solve_adjoint(self%ctx, evaluator, rhs, solution, error, tolerance=tol, preconditioner=precond, &
            & guess=lambda, workspace_owner=self%solver, parameters=self%adjoint)
      else if (associated(precond)) then
         call solve_adjoint(self%ctx, evaluator, rhs, solution, error, tolerance=tol, preconditioner=precond, &
            & workspace_owner=self%solver, parameters=self%adjoint)
      else if (warm) then
         call solve_adjoint(self%ctx, evaluator, rhs, solution, error, tolerance=tol, guess=lambda, &
            & workspace_owner=self%solver, parameters=self%adjoint)
      else
         call solve_adjoint(self%ctx, evaluator, rhs, solution, error, tolerance=tol, workspace_owner=self%solver, &
            & parameters=self%adjoint)
      end if
      if (allocated(error)) return
      call move_alloc(solution, lambda)

      allocate (mu, w_sr, w_lr, w_k, mold=self%tr)
      if (present(w_wu)) then
         allocate (w_wu(size(self%tr, 1), npair_from_ns(size(self%tr, 2)/nv)))
         call self%equation%adjoint(lambda, mu, error, wk=w_k, cr=self%c_work, ww=w_wu)
      else
         call self%equation%adjoint(lambda, mu, error, wk=w_k)
      end if
      if (allocated(error)) return
      !$omp parallel do collapse(2) schedule(static) num_threads(nt) private(a)
      do col = 1, size(mu, 2)
         do i = 1, size(mu, 1)
            a = self%rho(mod(col - 1, nv) + 1)*weight(i)
            w_sr(i, col) = uv_response_sr(a, dh(i, col), dc(i, col), mu(i, col), dcdu(i, col))
            w_lr(i, col) = uv_response_lr(a, dc(i, col), self%beta)
            w_k(i, col) = -self%beta*w_k(i, col)
         end do
      end do
      !$omp end parallel do
   end subroutine problem_response

   !> Take one solver update, damp it if allowed, and refresh c_s
   !>
   !> When the input iterate meets the tolerance the solve ends on a checked
   !> field: see finish_solve
   !>
   !> @param[inout] self     solve-local context
   !> @param[in]    iter     outer iteration
   !> @param[out]   rmse     RMS residual of the input iterate
   !> @param[out]   feasible whether an update was accepted
   !> @param[out]   error    solver or closure failure
   subroutine iteration_advance(self, iter, rmse, feasible, error)
      !> Solve-local context
      class(uv_iteration_type), intent(inout) :: self
      !> Outer iteration
      integer, intent(in) :: iter
      !> RMS residual
      real(wp), intent(out) :: rmse
      !> Whether the update was accepted
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      associate (it => iter); end associate
      associate (p => self%evaluator%problem)
         call p%solver%outer_step(self%evaluator, p%tr, self%gamma_new, rmse, feasible, error)
         if (allocated(error) .or. .not. feasible) return
         if (p%damp_small_k .and. p%solver%allows_postprocessing()) then
            self%delta = reshape(self%gamma_new - p%tr, [size(self%delta)])
            call p%damping%apply(self%delta)
            self%gamma_new = p%tr + reshape(self%delta, shape(p%tr))
         end if
         if (rmse < self%tolerance) then
            call finish_solve(self, rmse, feasible, error)
            return
         end if
         p%tr = self%gamma_new
         call p%closure%apply(p%u_sr, p%tr, p%beta, 1.0_wp, p%cr, error)
         if (allocated(error)) feasible = .false.
      end associate
   end subroutine iteration_advance

   !> End a converged solve on the better of the input iterate and the update
   !>
   !> The solver measured the input; the update gets a fresh residual and is
   !> kept only when it is not worse. rmse becomes the residual of the kept
   !> field and cr = C(tr)
   !>
   !> @param[inout] self     solve-local context with the update in gamma_new
   !> @param[inout] rmse     residual of the input iterate; of the kept field on return
   !> @param[out]   feasible whether c_s of the kept field is usable
   !> @param[out]   error    closure failure on the kept field
   subroutine finish_solve(self, rmse, feasible, error)
      !> Solve-local context
      class(uv_iteration_type), intent(inout) :: self
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
         call p%closure%apply(p%u_sr, p%tr, p%beta, 1.0_wp, p%cr, error)
         feasible = .not. allocated(error)
      end associate
   end subroutine finish_solve

   !> Print one row of the iteration table with the running free energy, kcal/mol
   !>
   !> @param[inout] self    solve-local context
   !> @param[in]    iter    outer iteration
   !> @param[in]    rmse    RMS residual
   !> @param[in]    elapsed advance time, s
   subroutine iteration_observe(self, iter, rmse, elapsed)
      !> Solve-local context
      class(uv_iteration_type), intent(inout) :: self
      !> Outer iteration
      integer, intent(in) :: iter
      !> RMS residual and advance time
      real(wp), intent(in) :: rmse, elapsed
      type(error_type), allocatable :: energy_error
      real(wp) :: energy

      if (.not. associated(self%printer)) return
      associate (p => self%evaluator%problem)
         call p%free_energy(energy, energy_error)
         if (allocated(energy_error)) then
            call p%solver%print_row(self%printer, iter, elapsed, rmse, hmax=maxval(p%cr + p%tr))
         else
            call p%solver%print_row(self%printer, iter, elapsed, rmse, hmax=maxval(p%cr + p%tr), &
               & energy=energy*autokcal)
         end if
      end associate
   end subroutine iteration_observe

   !> No file checkpoints: the owning model retains converged fields in memory
   !>
   !> @param[inout] self  solve-local context
   !> @param[out]   error always set
   subroutine iteration_checkpoint(self, error)
      !> Solve-local context
      class(uv_iteration_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      associate (s => self); end associate
      call fatal_error(error, "1D UV solve keeps no checkpoints")
   end subroutine iteration_checkpoint

   !> Residual R(t_s) = t_s - T(C(t_s))
   !>
   !> @param[in]  self     evaluator bound to a problem
   !> @param[in]  gamma    iterate t_s(r), (npts, nu nv)
   !> @param[out] residual residual, (npts, nu nv)
   !> @param[out] feasible whether the residual is usable
   !> @param[out] error    closure failure, 1 + h below zero, map failure or non-finite residual
   subroutine evaluator_residual(self, gamma, residual, feasible, error)
      !> Evaluator
      class(uv_evaluator_type), intent(in) :: self
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
         call fatal_error(error, "1D UV evaluator is not bound to a problem")
         return
      end if
      associate (p => self%problem)
         call p%closure%apply(p%u_sr, gamma, p%beta, 1.0_wp, p%c_work, error)
         if (allocated(error)) return
         if (.not. state_feasible(p%c_work, gamma, p%ctx%get_num_threads())) then
            call fatal_error(error, "1D UV closure gave 1 + h < 0")
            return
         end if
         call p%equation%map(p%c_work, residual, error)
         if (allocated(error)) return
         call residual_difference(gamma, residual, p%ctx%get_num_threads(), feasible)
         if (.not. feasible) call fatal_error(error, "1D UV residual is not finite")
      end associate
   end subroutine evaluator_residual

   !> Exact Jacobian-vector products are available once bound
   !>
   !> @param[in] self Evaluator
   logical function evaluator_has_jvp(self) result(available)
      !> Evaluator
      class(uv_evaluator_type), intent(in) :: self
      available = associated(self%problem)
   end function evaluator_has_jvp

   !> Cache the closure derivative c'(t_s); the map itself is affine
   !>
   !> @param[in]  self     evaluator bound to a problem
   !> @param[in]  gamma    linearization point t_s(r), (npts, nu nv)
   !> @param[out] feasible whether a linearization was stored
   !> @param[out] error    closure failure or 1 + h below zero
   subroutine evaluator_linearize(self, gamma, feasible, error)
      !> Evaluator
      class(uv_evaluator_type), intent(in) :: self
      !> Linearization point
      real(wp), intent(in) :: gamma(:, :)
      !> Whether the linearization is usable
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      feasible = .false.
      if (.not. associated(self%problem)) then
         call fatal_error(error, "1D UV evaluator is not bound to a problem")
         return
      end if
      associate (p => self%problem)
         p%linearized = .false.
         call p%closure%apply_and_derivative(p%u_sr, gamma, p%beta, 1.0_wp, p%c_work, p%dcdt, error)
         if (allocated(error)) return
         if (.not. state_feasible(p%c_work, gamma, p%ctx%get_num_threads())) then
            call fatal_error(error, "1D UV closure gave 1 + h < 0 at the linearization point")
            return
         end if
         feasible = .true.
         p%linearized = .true.
      end associate
   end subroutine evaluator_linearize

   !> Exact product J v = v - T'(c' v) at the last linearization
   !>
   !> @param[in]  self      evaluator with a linearization
   !> @param[in]  direction v, flat (npts nu nv)
   !> @param[out] jv        J v, flat (npts nu nv)
   !> @param[out] feasible  whether the product is usable
   !> @param[out] error     no linearization, transform failure or non-finite product
   subroutine evaluator_jvp(self, direction, jv, feasible, error)
      !> Evaluator
      class(uv_evaluator_type), intent(in) :: self
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
         call scale_flat(p%dcdt, direction, p%dc, p%ctx%get_num_threads())
         call p%equation%jvp(p%dc, p%dt, error)
         if (allocated(error)) then
            jv = 0.0_wp
            return
         end if
         call subtract_flat(direction, p%dt, jv, p%ctx%get_num_threads(), feasible)
         if (.not. feasible) call fatal_error(error, "1D UV Jacobian product is not finite")
      end associate
   end subroutine evaluator_jvp

   !> Transpose products are available once bound
   !>
   !> @param[in] self Evaluator
   logical function evaluator_has_adjoint(self) result(available)
      !> Evaluator
      class(uv_evaluator_type), intent(in) :: self
      available = associated(self%problem)
   end function evaluator_has_adjoint

   !> Exact transpose product J^T u = u - c' T'^T(u) at the last linearization
   !>
   !> @param[in]  self      evaluator with a linearization
   !> @param[in]  direction u, flat (npts nu nv)
   !> @param[out] jtv       J^T u, flat (npts nu nv)
   !> @param[out] feasible  whether the product is usable
   !> @param[out] error     no linearization, transform failure or non-finite product
   subroutine evaluator_adjoint(self, direction, jtv, feasible, error)
      !> Evaluator
      class(uv_evaluator_type), intent(in) :: self
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
         call unflatten(direction, p%dt, p%ctx%get_num_threads())
         call p%equation%adjoint(p%dt, p%dc, error)
         if (allocated(error)) then
            jtv = 0.0_wp
            return
         end if
         call subtract_flat(direction, p%dc, jtv, p%ctx%get_num_threads(), feasible, p%dcdt)
         if (.not. feasible) call fatal_error(error, "1D UV transpose product is not finite")
      end associate
   end subroutine evaluator_adjoint

   !> Whether the evaluator holds a linearization for flat vectors of this size
   !>
   !> @param[in]  self   evaluator
   !> @param[in]  nin    size of the direction
   !> @param[in]  nout   size of the product
   !> @param[out] error  unbound, not linearized or wrong size
   logical function linearized(self, nin, nout, error)
      !> Evaluator
      class(uv_evaluator_type), intent(in) :: self
      !> Vector sizes
      integer, intent(in) :: nin, nout
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      linearized = .false.
      if (.not. associated(self%problem)) then
         call fatal_error(error, "1D UV evaluator is not bound to a problem")
      else if (.not. self%problem%linearized) then
         call fatal_error(error, "1D UV Jacobian product requested without a linearization")
      else if (nin /= size(self%problem%tr) .or. nout /= nin) then
         call fatal_error(error, "1D UV Jacobian vectors must have npts nu nv entries")
      else
         linearized = .true.
      end if
   end function linearized

end module moist_model_moz_1d_driver
