!> 3D UV solve: residual evaluator, outer iterations and solvation free energy
!>
!> - Iterate is the short-range t_s(r) on the grid
!> - Closure gives c_s(r), equation maps it back, R(t_s) = t_s - T(C(t_s)), at full coupling
!> - Cold starts use T(c_s = 0) or t_s = 0
!>   - Supplied guess or, when reuse is requested, the checked previous solution takes precedence
!>   - Starting field meeting the tolerance is the solution, without iterations
!>   - Failed warm attempt retries cold once
!> - Optional checkpoint store receives intermediate and final fields
!> - Exact Jacobian products J v = v - T'(c' v), T affine
!> - Solvation free energy dG = (1/beta) sum_a rho_a int Phi(h_a, c_a) dV
!>   with physical h and c, by the grid's quadrature
!>   - Coulomb tails of c cancel in the sum over the sites of a neutral solvent
!> - Direct-correlation moment D = sum_a rho_a int c_a dV of the physical c,
!>   for the partial molar volume
module moist_model_moz_3d_driver
!$ use omp_lib, only: omp_get_thread_num, omp_in_parallel
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_convert, only: autokcal
   use mctc_io_utils, only: to_lower
   use moist_context, only: moist_context_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_grid_3d_integration, only: integration_block_size, integration_parallel_threshold
   use moist_math_grid_3d_kernel_integration, only: merge_compensated_sums
   use moist_model_moz_closure, only: moist_model_moz_closure_type, new_closure
   use moist_model_moz_solver, only: moist_model_moz_solver_type, moist_model_moz_solver_evaluator_type, &
      & moist_model_moz_solver_iteration_type, moist_model_moz_solver_printer_type, &
      & moist_model_moz_solver_statistics_type, new_solver, run_iterations, solve_adjoint
   use moist_model_moz_fields, only: state_feasible, scale_flat, subtract_flat, residual_difference, unflatten
   use moist_model_moz_response_kernel, only: uv_response_rhs, uv_response_sr, uv_response_lr, uv_response_measure
   use moist_model_moz_thermodynamics_kernel, only: scale_energy_derivative
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type
   use moist_model_moz_3d_equation, only: uv3d_equation_type
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type, cold_guess_modes
   use moist_model_moz_checkpoint, only: moz_checkpoint_type
   implicit none(type, external)
   private

   public :: uv3d_problem_type, uv3d_evaluator_type, new_uv3d_problem

   !> Fixed 3D UV data, operators and the iterate of one solve
   !>
   !> - Do not copy a prepared problem with allocate(source=...)
   !> - Bypasses the equation's assignment, can share transform plans and grid bindings
   type :: uv3d_problem_type
      !> Run context, borrowed for the lifetime of the problem
      type(moist_context_type), pointer :: ctx => null()
      !> Equation with the owned grid, the response table and the transforms
      type(uv3d_equation_type) :: equation
      !> Closure
      class(moist_model_moz_closure_type), allocatable :: closure
      !> Nonlinear solver
      class(moist_model_moz_solver_type), allocatable :: solver
      !> Short- and long-range real-space potential, Hartree, (ngrid, nv)
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :)
      !> Solvent site number densities, 1/bohr^3 (nv)
      real(wp), allocatable :: rho(:)
      !> Inverse thermal energy, 1/Hartree
      real(wp) :: beta = 0.0_wp
      !> Accepted iterate t_s(r) and its c_s(r), (ngrid, nv)
      real(wp), allocatable :: tr(:, :), cr(:, :)
      !> Closure derivative c'(t_s) of the last linearization and scratch, (ngrid, nv)
      real(wp), allocatable :: dcdt(:, :), c_work(:, :), dc(:, :), dt(:, :)
      !> Whether dcdt holds a linearization
      logical :: linearized = .false.
      !> Whether tr can seed the next solve on the same grid
      logical, private :: has_guess = .false.
      !> Cold path used for first solves and failed or infeasible retained guesses
      character(len=8), private :: cold_guess = "cs0"
   contains
      !> Replace potential tables while retaining the fixed-grid equation and solver storage
      procedure :: update_potential => problem_update_potential
      !> Solve with the selected cold guess or the previous converged fields
      procedure :: solve => problem_solve
      !> Solvation free energy of the current iterate
      procedure :: free_energy => problem_free_energy
      !> Direct-correlation moment of the current iterate
      procedure :: direct_moment => problem_direct_moment
      !> Energy response to potential tables, optionally also to the grid geometry
      procedure :: response => problem_response
   end type uv3d_problem_type

   !> Per-block error ownership for deterministic reporting after parallel energy evaluation
   type :: uv3d_energy_error_type
      !> Error of one energy block, unallocated on success
      type(error_type), allocatable :: error
   end type uv3d_energy_error_type

   !> Residual and Jacobian of a borrowed problem; the target must outlive evaluator calls
   type, extends(moist_model_moz_solver_evaluator_type) :: uv3d_evaluator_type
      !> Problem whose workspaces the evaluator uses
      type(uv3d_problem_type), pointer :: problem => null()
   contains
      procedure :: residual => evaluator_residual
      procedure :: has_jvp => evaluator_has_jvp
      procedure :: linearize => evaluator_linearize
      procedure :: jvp => evaluator_jvp
      procedure :: has_adjoint => evaluator_has_adjoint
      procedure :: adjoint => evaluator_adjoint
   end type uv3d_evaluator_type

   !> Solve-local outer iteration state
   type, extends(moist_model_moz_solver_iteration_type) :: uv3d_iteration_type
      !> Evaluator borrowing the problem
      type(uv3d_evaluator_type) :: evaluator
      !> Proposed iterate and its residual, (ngrid, nv)
      real(wp), allocatable :: gamma_new(:, :), residual(:, :)
      !> Borrowed iteration table; null when quiet
      type(moist_model_moz_solver_printer_type), pointer :: printer => null()
      !> Tolerance; on the converging step the update is checked against its input
      real(wp) :: tolerance = 0.0_wp
      !> Borrowed store of intermediate saves; null without checkpoints
      class(moz_checkpoint_type), pointer :: store => null()
      !> Iterations of the earlier attempts; makes the saved iteration counts cumulative
      integer :: offset = 0
   contains
      procedure :: advance => iteration_advance
      procedure :: observe => iteration_observe
      procedure :: checkpoint => iteration_checkpoint
   end type uv3d_iteration_type

contains

   !> Build a 3D UV problem from the fixed tables and the settings
   !>
   !> @param[out] self       prepared problem; solve initializes the cold or retained guess
   !> @param[in]  ctx        run context; verbosity, output unit and OpenMP team size
   !> @param[in]  grid       updated grid, copied into the equation
   !> @param[in]  response   solvent response on the grid's k-points
   !> @param[in]  u_sr       short-range potential, Hartree, (ngrid, nv)
   !> @param[in]  ur_lr      long-range real-space potential, Hartree, (ngrid, nv)
   !> @param[in]  uk_lr      long-range reciprocal potential, Hartree bohr^3, (npts_k, nv)
   !> @param[in]  rho        solvent site number densities, 1/bohr^3 (nv)
   !> @param[in]  beta       inverse thermal energy, 1/Hartree
   !> @param[in]  params     validated settings
   !> @param[out] error      shape mismatch, closure, solver or equation setup failure
   subroutine new_uv3d_problem(self, ctx, grid, response, u_sr, ur_lr, uk_lr, rho, beta, params, error)
      !> Problem
      type(uv3d_problem_type), intent(out) :: self
      !> Run context, borrowed for the lifetime of the problem
      type(moist_context_type), intent(in), target :: ctx
      !> Updated grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Solvent response
      type(uv3d_susceptibility_type), intent(in) :: response
      !> Real-space potential tables
      real(wp), intent(in) :: u_sr(:, :), ur_lr(:, :)
      !> Long-range reciprocal potential
      complex(wp), intent(in) :: uk_lr(:, :)
      !> Solvent site densities
      real(wp), intent(in) :: rho(:)
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Settings
      type(model_moz_3d_parameters_type), intent(in) :: params
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      self%ctx => ctx
      self%cold_guess = to_lower(trim(adjustl(params%cold_guess)))
      if (.not. any(self%cold_guess == cold_guess_modes)) then
         call fatal_error(error, "Unknown 3D MOZ cold guess: "//trim(params%cold_guess)//" (choose cs0 or ts0)")
         return
      end if
      if (size(rho) < 1 .or. size(rho) /= response%nv) then
         call fatal_error(error, "3D UV needs one density per solvent site in the response")
         return
      end if
      if (any(shape(u_sr) /= [grid%ngrid, size(rho)]) .or. any(shape(ur_lr) /= shape(u_sr))) then
         call fatal_error(error, "3D UV potential tables must have shape (ngrid, nv)")
         return
      end if
      call new_closure(self%closure, params%closure, error, self%ctx)
      if (allocated(error)) return
      call new_solver(self%solver, ctx, error, params%solver)
      if (allocated(error)) return
      self%equation%nthreads = ctx%get_num_threads()
      call self%equation%setup(grid, response, uk_lr, beta, error)
      if (allocated(error)) return
      self%u_sr = u_sr
      self%ur_lr = ur_lr
      self%rho = rho
      self%beta = beta
      allocate (self%tr, self%dcdt, self%c_work, self%dc, self%dt, mold=u_sr)
      allocate (self%cr(size(u_sr, 1), size(u_sr, 2)), source=0.0_wp)
   end subroutine new_uv3d_problem

   !> Replace host-dependent tables at fixed geometry, bulk, beta and solver settings
   !>
   !> - Keeps physical t = h - c when changing the Ng split
   !>   - t_s(new) = t_s(old) + beta*(u_lr(old) - u_lr(new))
   !> - Next solve checks this guess against the new closure and equation
   !>
   !> @param[inout] self   constructed problem; equation, plans and solver scratch retained
   !> @param[in]    u_sr   short-range potential, Hartree, (ngrid, nv)
   !> @param[in]    ur_lr  long-range real-space potential, Hartree, (ngrid, nv)
   !> @param[in]    uk_lr  long-range reciprocal potential, Hartree bohr^3, (npts_k, nv)
   !> @param[out]   error  unset problem or shape mismatch; old tables retained on refusal
   subroutine problem_update_potential(self, u_sr, ur_lr, uk_lr, error)
      !> Problem
      class(uv3d_problem_type), intent(inout) :: self
      !> Short-range and long-range real-space potentials, Hartree, (ngrid, nv)
      real(wp), intent(in) :: u_sr(:, :), ur_lr(:, :)
      !> Long-range reciprocal potential, Hartree bohr^3, (npts_k, nv)
      complex(wp), intent(in) :: uk_lr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(self%u_sr)) then
         call fatal_error(error, "Construct the 3D UV problem before replacing its potential")
         return
      end if
      if (any(shape(u_sr) /= shape(self%u_sr)) .or. any(shape(ur_lr) /= shape(self%ur_lr))) then
         call fatal_error(error, "3D UV potential update must keep the field shape")
         return
      end if
      call self%equation%update_potential(uk_lr, error)
      if (allocated(error)) return
      if (self%has_guess) self%tr = self%tr + self%beta*(self%ur_lr - ur_lr)
      self%u_sr(:, :) = u_sr
      self%ur_lr(:, :) = ur_lr
      self%linearized = .false.
   end subroutine problem_update_potential

   !> Solve the 3D UV equations at full coupling
   !>
   !> - Cold starts use t_s = T(0) (cs0, default) or t_s = 0 (ts0)
   !> - A supplied guess, else a checked retained guess when requested, takes precedence
   !> - A usable starting field meeting the tolerance is the solution, without iterations
   !> - Infeasible guesses and failed warm solves use the selected cold path
   !> - A rejected supplied guess, or a warm solve that failed with an error, warns
   !> - Solver history is cleared for each attempt; plans and allocated workspace persist
   !> - A converged solve ends on a field whose residual was measured below the
   !>   tolerance, with cr = C(tr); rmse is the residual of that field
   !> - With a store: saves every `interval` iterations, and the final field when
   !>   the solve iterated without error
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
   !> @param[in]    reuse_initial_guess  try the retained converged field; default false for standalone calls
   !> @param[in]    guess      initial short-range t_s(r), (ngrid, nv); takes precedence over reuse
   !> @param[inout] checkpoint store of the saved states, with its header set
   subroutine problem_solve(self, tolerance, maxiter, converged, rmse, niter, error, reuse_initial_guess, guess, &
      & checkpoint)
      !> Problem
      class(uv3d_problem_type), intent(inout), target :: self
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
      !> Use the previously converged fields on this geometry
      logical, intent(in), optional :: reuse_initial_guess
      !> Initial short-range indirect correlation
      real(wp), intent(in), optional :: guess(:, :)
      !> Store of the saved states
      class(moz_checkpoint_type), allocatable, intent(inout), optional, target :: checkpoint
      type(uv3d_iteration_type) :: context
      type(moist_model_moz_solver_printer_type), target :: printer
      type(moist_model_moz_solver_statistics_type) :: work0
      type(error_type), allocatable :: guess_error
      integer :: attempt, attempted_steps, interval
      logical :: warm

      converged = .false.
      rmse = huge(1.0_wp)
      niter = 0
      context%evaluator%problem => self
      context%tolerance = tolerance
      allocate (context%gamma_new, context%residual, mold=self%tr)
      interval = 0
      if (present(checkpoint)) then
         if (allocated(checkpoint)) then
            context%store => checkpoint
            interval = checkpoint%interval
         end if
      end if
      if (self%ctx%verbosity(4)) then
         call printer%setup(unit=self%ctx%unit, energy_label="DG kcal", denergy_label="dDG kcal", &
            & inner_rows=self%ctx%verbosity(5))
         context%printer => printer
         call self%solver%set_iteration_printer(printer)
         call printer%header("3D UV")
         call printer%blank()
         if (associated(context%store)) then
            if (len(context%store%note()) > 0) call printer%note(context%store%note())
         end if
      end if
      work0 = self%solver%get_statistics()

      self%linearized = .false.
      warm = .false.
      if (present(guess)) then
         ! A rejected caller guess falls back to the cold path and warns why
         if (any(shape(guess) /= shape(self%tr))) then
            call self%ctx%warning("3D UV: guess shape does not match the grid, ignoring guess")
         else
            call context%evaluator%residual(guess, context%residual, warm, guess_error)
            warm = warm .and. .not. allocated(guess_error)
            if (warm) then
               self%tr = guess
               self%cr = self%c_work
            else if (allocated(guess_error)) then
               call self%ctx%warning("3D UV: guess is infeasible, ignoring guess: "//guess_error%message)
            else
               call self%ctx%warning("3D UV: guess is infeasible, ignoring guess")
            end if
         end if
      else if (present(reuse_initial_guess)) then
         if (reuse_initial_guess .and. self%has_guess) then
            call context%evaluator%residual(self%tr, context%residual, warm, guess_error)
            warm = warm .and. .not. allocated(guess_error)
            if (warm) self%cr = self%c_work
         end if
      end if
      if (warm) call accept_guess(context, tolerance, rmse, converged)
      self%has_guess = .false.
      do attempt = 1, 2
         if (converged) exit
         if (.not. warm) then
            self%cr = 0.0_wp
            if (self%cold_guess == "ts0") then
               self%tr = 0.0_wp
            else
               call self%equation%map(self%cr, self%tr, error)
               if (allocated(error)) exit
            end if
         end if
         call self%solver%begin_stage()
         if (associated(context%printer)) then
            if (attempt > 1) then
               call printer%separator()
               call printer%blank()
               call printer%note("Retry with cold guess "//trim(self%cold_guess))
            end if
            call printer%note("Solver: "//self%solver%name())
            if (.not. warm) call printer%note("Cold guess: "//trim(self%cold_guess))
            call printer%blank()
            call printer%begin_stage()
         end if
         context%offset = niter
         call run_iterations(context, self%ctx, maxiter, tolerance, interval, rmse, attempted_steps, &
            & converged, error)
         niter = niter + attempted_steps
         if (associated(context%printer)) call printer%note(self%solver%stage_summary())
         if (converged .and. .not. allocated(error)) exit
         if (.not. warm) exit
         ! A poor previous field must not prevent the selected cold path from succeeding
         if (allocated(error)) then
            call self%ctx%warning("3D UV: warm start failed: "//error%message//"; retrying cold")
            deallocate (error)
         end if
         warm = .false.
         self%linearized = .false.
      end do
      if (associated(context%store) .and. niter > 0 .and. .not. allocated(error)) then
         call save_state(self, context%store, rmse, niter, converged)
      end if

      call self%solver%set_iteration_printer()
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
   !> dG = (1/beta) sum_a rho_a int Phi(h_a, c_a) dV with the physical
   !> c = c_s - beta u_lr and h = c_s + t_s
   !>
   !> @param[in]  self    problem with an accepted iterate
   !> @param[out] energy  solvation free energy, Hartree
   !> @param[out] error   invalid correlations or quadrature failure
   subroutine problem_free_energy(self, energy, error)
      !> Problem
      class(uv3d_problem_type), intent(in) :: self
      !> Solvation free energy
      real(wp), intent(out) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: phi(:, :)
      real(wp) :: integral
      integer :: col

      energy = 0.0_wp
      if (size(self%cr, 1) >= integration_parallel_threshold) then
         call problem_blocked_free_energy(self, energy, error)
         return
      end if
      allocate (phi, mold=self%cr)
      call self%closure%free_energy_integrands(self%cr + self%tr, self%cr - self%beta*self%ur_lr, phi, error)
      if (allocated(error)) return
      do col = 1, size(phi, 2)
         call self%equation%integrate(phi(:, col), integral, error)
         if (allocated(error)) return
         energy = energy + self%rho(col)*integral
      end do
      energy = energy/self%beta
   end subroutine problem_free_energy

   !> Evaluate physical correlations and closure energy in bounded worker scratch
   !>
   !> - Grid owns the measure and each compensated block
   !> - Block partition and final merge order independent of the worker count
   !> - Checked closure binding retains custom integrands and invalid-state diagnostics
   !>
   !> @param[in]  self    problem with accepted c_s and t_s
   !> @param[out] energy  solvation free energy, Hartree, defined on success
   !> @param[out] error   first invalid block in grid order
   subroutine problem_blocked_free_energy(self, energy, error)
      !> Problem
      class(uv3d_problem_type), intent(in) :: self
      !> Solvation free energy, Hartree
      real(wp), intent(out) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: h(:, :, :), c(:, :, :), phi(:, :, :), partial(:, :, :)
      type(uv3d_energy_error_type), allocatable :: errors(:)
      integer :: npts, nv, nb, nt, ib, first, n, it, col, i

      energy = 0.0_wp
      npts = size(self%cr, 1)
      nv = size(self%cr, 2)
      nb = 1 + (npts - 1)/integration_block_size
      nt = min(self%ctx%get_num_threads(), nb)
!$    if (omp_in_parallel()) nt = 1
      allocate (h(integration_block_size, nv, nt), c(integration_block_size, nv, nt), &
         & phi(integration_block_size, nv, nt), partial(2, nb, nv), errors(nb))
      !$omp parallel do schedule(static) num_threads(nt) if(nt > 1) private(first, n, it, col, i)
      do ib = 1, nb
         first = 1 + (ib - 1)*integration_block_size
         n = min(integration_block_size, npts - first + 1)
         it = 1
!$       it = omp_get_thread_num() + 1
         do col = 1, nv
            do i = 1, n
               h(i, col, it) = self%cr(first + i - 1, col) + self%tr(first + i - 1, col)
               c(i, col, it) = self%cr(first + i - 1, col) - self%beta*self%ur_lr(first + i - 1, col)
            end do
         end do
         call self%closure%free_energy_integrands(h(:n, :, it), c(:n, :, it), phi(:n, :, it), errors(ib)%error)
         if (allocated(errors(ib)%error)) cycle
         do col = 1, nv
            call self%equation%integrate_block(phi(:n, col, it), first, partial(:, ib, col), errors(ib)%error)
            if (allocated(errors(ib)%error)) exit
         end do
      end do
      !$omp end parallel do
      do ib = 1, nb
         if (allocated(errors(ib)%error)) then
            call move_alloc(errors(ib)%error, error)
            return
         end if
      end do
      do col = 1, nv
         energy = energy + self%rho(col)*merge_compensated_sums(partial(:, :, col))
      end do
      energy = energy/self%beta
   end subroutine problem_blocked_free_energy

   !> Direct-correlation moment of the current iterate
   !>
   !> - D = sum_a rho_a int c_a dV
   !> - Physical c = c_s - beta u_lr
   !>
   !> @param[in]  self    problem with an accepted iterate
   !> @param[out] moment  D, dimensionless
   !> @param[out] error   quadrature failure
   subroutine problem_direct_moment(self, moment, error)
      !> Problem
      class(uv3d_problem_type), intent(in) :: self
      !> Moment
      real(wp), intent(out) :: moment
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: integral
      integer :: col

      moment = 0.0_wp
      do col = 1, size(self%cr, 2)
         call self%equation%integrate(self%cr(:, col) - self%beta*self%ur_lr(:, col), integral, error)
         if (allocated(error)) return
         moment = moment + self%rho(col)*integral
      end do
   end subroutine problem_direct_moment

   !> Pull the energy response back to the potential tables on the fixed grid and solvent
   !>
   !> - Uses the accumulated energy and moment seeds, one solve J^T lambda = dE/dt_s
   !> - shift = beta*dE/dD, dimensionless; free_energy_weight = dE/dGraw
   !> - Quadrature weights bohr^3, rho bohr^-3, beta Hartree^-1
   !> - w_k = Re(sum(conjg(w_k)*du_lr(k))); its -beta factor is in the equation
   !> - tolerance relative to the initial RHS RMS; convergence checks the true residual
   !> - Field-shaped lambda on entry is a warm start, else zero
   !> - Optional grid accumulator also receives explicit quadrature and implicit transform sensitivities
   !>   - R = t - T(C) gives -lambda^T R_R = +lambda^T T_R
   !>
   !> @param[inout] self                problem with a converged iterate
   !> @param[in]    shift               beta*dE/dD, dimensionless
   !> @param[in]    tolerance           RMS of the transpose residual relative to the RMS of the right-hand side
   !> @param[out]   w_sr                adjoint of u_sr, (ngrid, nv)
   !> @param[out]   w_lr                adjoint of u_lr(r), (ngrid, nv)
   !> @param[out]   w_k                 adjoint of u_lr(k), (npts_k, nv)
   !> @param[inout] lambda              transposed solution, (ngrid, nv); a field-shaped value on entry is the start
   !> @param[out]   error               infeasible iterate, bad input, closure failure or unconverged transposed solve
   !> @param[in]    free_energy_weight  optional free-energy sensitivity dE/dGraw; default one
   !> @param[inout] grid_adjoint        optional grid-geometry accumulator, ngrid entries
   subroutine problem_response(self, shift, tolerance, w_sr, w_lr, w_k, lambda, error, free_energy_weight, grid_adjoint)
      !> Problem
      class(uv3d_problem_type), intent(inout), target :: self
      !> Moment sensitivity beta*dE/dD, and relative transpose tolerance
      real(wp), intent(in) :: shift, tolerance
      !> Adjoints of u_sr and u_lr(r), (ngrid, nv)
      real(wp), allocatable, intent(out) :: w_sr(:, :), w_lr(:, :)
      !> Adjoint of u_lr(k), (npts_k, nv)
      complex(wp), allocatable, intent(out) :: w_k(:, :)
      !> Transposed solution, (ngrid, nv)
      real(wp), allocatable, intent(inout) :: lambda(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Free-energy sensitivity dE/dGraw; default one
      real(wp), optional, intent(in) :: free_energy_weight
      !> Grid-geometry adjoint accumulator
      type(volume_adjoint_type), optional, intent(inout) :: grid_adjoint
      type(uv3d_evaluator_type), target :: evaluator
      real(wp), allocatable :: weight(:), dphi_dh(:, :), dphi_dc(:, :), rhs(:, :), solution(:, :)
      real(wp) :: tol, energy_seed, a
      integer :: i, col, nt
      logical :: feasible, warm

      if (.not. (ieee_is_finite(shift) .and. ieee_is_finite(tolerance)) .or. tolerance <= 0.0_wp) then
         call fatal_error(error, "3D UV response needs a finite shift and a positive tolerance")
         return
      end if
      energy_seed = 1.0_wp
      if (present(free_energy_weight)) energy_seed = free_energy_weight
      if (.not. ieee_is_finite(energy_seed)) then
         call fatal_error(error, "3D UV response needs a finite free-energy sensitivity")
         return
      end if
      if (.not. allocated(self%tr) .or. .not. allocated(self%closure)) then
         call fatal_error(error, "3D UV response needs a constructed problem")
         return
      end if
      if (present(grid_adjoint)) then
         if (.not. grid_adjoint%is_initialized() .or. grid_adjoint%size() /= size(self%tr, 1)) then
            call fatal_error(error, "3D UV response needs a grid accumulator with ngrid entries")
            return
         end if
      end if
      nt = self%ctx%get_num_threads()
      evaluator%problem => self
      call evaluator%linearize(self%tr, feasible, error)
      if (allocated(error)) return
      if (.not. feasible) then
         call fatal_error(error, "3D UV response needs a feasible iterate")
         return
      end if
      allocate (dphi_dh, dphi_dc, w_sr, rhs, mold=self%tr)
      ! Phase-local names expose the physical fields without extra whole-field storage
      associate (h_phys => self%dc, c_phys => self%dt)
         !$omp parallel do collapse(2) schedule(static) num_threads(nt)
         do col = 1, size(self%tr, 2)
            do i = 1, size(self%tr, 1)
               h_phys(i, col) = self%c_work(i, col) + self%tr(i, col)
               c_phys(i, col) = self%c_work(i, col) - self%beta*self%ur_lr(i, col)
            end do
         end do
         !$omp end parallel do
         if (present(grid_adjoint)) then
            ! The RHS buffer is available for Phi until the adjoint RHS is assembled
            associate (phi => rhs)
               call self%closure%free_energy_integrands(h_phys, c_phys, phi, error)
               if (allocated(error)) return
               !$omp parallel do schedule(static) num_threads(nt)
               do i = 1, size(phi, 1)
                  grid_adjoint%w_w(i) = grid_adjoint%w_w(i) + &
                     & sum(uv_response_measure(self%rho, self%beta, energy_seed, phi(i, :), shift, c_phys(i, :)))
               end do
               !$omp end parallel do
            end associate
         end if
         call self%closure%free_energy_derivative_fields(h_phys, c_phys, dphi_dh, dphi_dc, error)
         if (allocated(error)) return
      end associate
      associate (dcdu => w_sr)
         call self%closure%potential_derivative(self%u_sr, self%tr, self%beta, 1.0_wp, dcdu, error)
         if (allocated(error)) return
      end associate
      call self%equation%integration_weights(weight, error)
      if (allocated(error)) return
      weight = weight/self%beta
      ! Seed the integrand derivatives, then add the thermodynamic moment sensitivity
      !$omp parallel do collapse(2) schedule(static) num_threads(nt) private(a)
      do col = 1, size(rhs, 2)
         do i = 1, size(rhs, 1)
            if (energy_seed /= 1.0_wp) then
               dphi_dh(i, col) = scale_energy_derivative(energy_seed, dphi_dh(i, col))
               dphi_dc(i, col) = scale_energy_derivative(energy_seed, dphi_dc(i, col))
            end if
            dphi_dc(i, col) = dphi_dc(i, col) + shift
            a = self%rho(col)*weight(i)
            rhs(i, col) = uv_response_rhs(a, dphi_dh(i, col), dphi_dc(i, col), self%dcdt(i, col))
         end do
      end do
      !$omp end parallel do

      tol = max(tolerance*sqrt(sum(rhs**2)/real(size(rhs), wp)), tiny(1.0_wp))
      warm = .false.
      if (allocated(lambda)) warm = all(shape(lambda) == shape(rhs))
      ! A zero objective has an exact zero solution; discard a nonzero warm start
      warm = warm .and. any(rhs /= 0.0_wp)
      if (warm) then
         call solve_adjoint(self%ctx, evaluator, rhs, solution, error, tolerance=tol, guess=lambda, &
            & workspace_owner=self%solver)
      else
         call solve_adjoint(self%ctx, evaluator, rhs, solution, error, tolerance=tol, workspace_owner=self%solver)
      end if
      if (allocated(error)) return
      call move_alloc(solution, lambda)
      deallocate (rhs)

      allocate (w_lr, mold=self%tr)
      ! mu = T'^T lambda; consume dc_s/du_sr pointwise as its buffer becomes w_sr
      associate (mu => self%dc, dcdu => w_sr)
         if (present(grid_adjoint)) then
            call self%equation%adjoint(lambda, mu, error, wk=w_k, cr=self%c_work, grid_adjoint=grid_adjoint)
         else
            call self%equation%adjoint(lambda, mu, error, wk=w_k)
         end if
         if (allocated(error)) return
         !$omp parallel do collapse(2) schedule(static) num_threads(nt) private(a)
         do col = 1, size(self%tr, 2)
            do i = 1, size(self%tr, 1)
               a = self%rho(col)*weight(i)
               w_sr(i, col) = uv_response_sr(a, dphi_dh(i, col), dphi_dc(i, col), mu(i, col), dcdu(i, col))
               w_lr(i, col) = uv_response_lr(a, dphi_dc(i, col), self%beta)
            end do
         end do
         !$omp end parallel do
      end associate
   end subroutine problem_response

   !> Take one solver update and refresh c_s
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
      class(uv3d_iteration_type), intent(inout) :: self
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
      class(uv3d_iteration_type), intent(inout) :: self
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
      class(uv3d_iteration_type), intent(inout) :: self
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

   !> Save the accepted iterate as an unconverged state
   !>
   !> - A failed save warns; the solve continues
   !>
   !> @param[inout] self  solve-local context
   !> @param[in]    iter  attempt-local iteration of the accepted iterate
   !> @param[in]    rmse  residual measured by its advance
   !> @param[out]   error never set
   subroutine iteration_checkpoint(self, iter, rmse, error)
      !> Solve-local context
      class(uv3d_iteration_type), intent(inout) :: self
      !> Outer iteration
      integer, intent(in) :: iter
      !> RMS residual
      real(wp), intent(in) :: rmse
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. associated(self%store)) return
      call save_state(self%evaluator%problem, self%store, rmse, self%offset + iter, .false.)
   end subroutine iteration_checkpoint

   !> Accept a feasible starting field that already meets the tolerance
   !>
   !> - Uses the residual the guess check left in the context
   !>
   !> @param[inout] context    solve-local context with the residual of the starting field
   !> @param[in]    tolerance  RMS residual of a converged solve
   !> @param[out]   rmse       RMS residual of the starting field
   !> @param[out]   converged  whether it meets the tolerance
   subroutine accept_guess(context, tolerance, rmse, converged)
      !> Solve-local context
      type(uv3d_iteration_type), intent(inout) :: context
      !> Tolerance
      real(wp), intent(in) :: tolerance
      !> RMS residual of the starting field
      real(wp), intent(out) :: rmse
      !> Whether it meets the tolerance
      logical, intent(out) :: converged
      character(len=16) :: buf

      rmse = sqrt(sum(context%residual**2)/real(size(context%residual), wp))
      converged = rmse < tolerance
      if (.not. associated(context%printer)) return
      write (buf, "(es9.2)") rmse
      if (converged) then
         call context%printer%note("Initial guess: RMS residual "//trim(adjustl(buf))//", accepted")
      else
         call context%printer%note("Initial guess: RMS residual "//trim(adjustl(buf)))
      end if
   end subroutine accept_guess

   !> Save the current iterate as a state of the store; a failed save warns
   !>
   !> @param[in]    problem    problem with an accepted iterate
   !> @param[inout] store      store with its header set
   !> @param[in]    rmse       last measured RMS residual
   !> @param[in]    niter      outer iterations so far
   !> @param[in]    converged  whether the field is converged
   subroutine save_state(problem, store, rmse, niter, converged)
      !> Problem
      type(uv3d_problem_type), intent(in) :: problem
      !> Store
      class(moz_checkpoint_type), intent(inout) :: store
      !> Last measured RMS residual
      real(wp), intent(in) :: rmse
      !> Outer iterations so far
      integer, intent(in) :: niter
      !> Whether the field is converged
      logical, intent(in) :: converged
      type(error_type), allocatable :: error

      call store%record(problem%tr, problem%beta, problem%ur_lr, rmse, 1.0_wp, niter, converged, error)
      if (allocated(error)) call problem%ctx%warning("3D UV checkpoint not saved: "//error%message)
   end subroutine save_state

   !> Residual R(t_s) = t_s - T(C(t_s))
   !>
   !> @param[in]  self     evaluator bound to a problem
   !> @param[in]  gamma    iterate t_s(r), (ngrid, nv)
   !> @param[out] residual residual, (ngrid, nv)
   !> @param[out] feasible whether the residual is usable
   !> @param[out] error    closure failure, 1 + h below zero, map failure or non-finite residual
   subroutine evaluator_residual(self, gamma, residual, feasible, error)
      !> Evaluator
      class(uv3d_evaluator_type), intent(in) :: self
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
         call fatal_error(error, "3D UV evaluator is not bound to a problem")
         return
      end if
      associate (p => self%problem)
         call p%closure%apply(p%u_sr, gamma, p%beta, 1.0_wp, p%c_work, error)
         if (allocated(error)) return
         if (.not. state_feasible(p%c_work, gamma, p%ctx%get_num_threads())) then
            call fatal_error(error, "3D UV closure gave 1 + h < 0")
            return
         end if
         call p%equation%map(p%c_work, residual, error)
         if (allocated(error)) return
         call residual_difference(gamma, residual, p%ctx%get_num_threads(), feasible)
         if (.not. feasible) call fatal_error(error, "3D UV residual is not finite")
      end associate
   end subroutine evaluator_residual

   !> Exact Jacobian-vector products are available once bound
   !>
   !> @param[in] self Evaluator
   logical function evaluator_has_jvp(self) result(available)
      !> Evaluator
      class(uv3d_evaluator_type), intent(in) :: self
      available = associated(self%problem)
   end function evaluator_has_jvp

   !> Cache the closure derivative c'(t_s); the map itself is affine
   !>
   !> @param[in]  self     evaluator bound to a problem
   !> @param[in]  gamma    linearization point t_s(r), (ngrid, nv)
   !> @param[out] feasible whether a linearization was stored
   !> @param[out] error    closure failure or 1 + h below zero
   subroutine evaluator_linearize(self, gamma, feasible, error)
      !> Evaluator
      class(uv3d_evaluator_type), intent(in) :: self
      !> Linearization point
      real(wp), intent(in) :: gamma(:, :)
      !> Whether the linearization is usable
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      feasible = .false.
      if (.not. associated(self%problem)) then
         call fatal_error(error, "3D UV evaluator is not bound to a problem")
         return
      end if
      associate (p => self%problem)
         p%linearized = .false.
         call p%closure%apply_and_derivative(p%u_sr, gamma, p%beta, 1.0_wp, p%c_work, p%dcdt, error)
         if (allocated(error)) return
         if (.not. state_feasible(p%c_work, gamma, p%ctx%get_num_threads())) then
            call fatal_error(error, "3D UV closure gave 1 + h < 0 at the linearization point")
            return
         end if
         feasible = .true.
         p%linearized = .true.
      end associate
   end subroutine evaluator_linearize

   !> Exact product J v = v - T'(c' v) at the last linearization
   !>
   !> @param[in]  self      evaluator with a linearization
   !> @param[in]  direction v, flat (ngrid nv)
   !> @param[out] jv        J v, flat (ngrid nv)
   !> @param[out] feasible  whether the product is usable
   !> @param[out] error     no linearization, transform failure or non-finite product
   subroutine evaluator_jvp(self, direction, jv, feasible, error)
      !> Evaluator
      class(uv3d_evaluator_type), intent(in) :: self
      !> Direction
      real(wp), intent(in) :: direction(:)
      !> Product
      real(wp), intent(out) :: jv(:)
      !> Whether the product is usable
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      feasible = .false.
      if (.not. associated(self%problem)) then
         jv = 0.0_wp
         call fatal_error(error, "3D UV evaluator is not bound to a problem")
         return
      end if
      associate (p => self%problem)
         if (.not. p%linearized) then
            jv = 0.0_wp
            call fatal_error(error, "3D UV Jacobian product requested without a linearization")
            return
         end if
         if (size(direction) /= size(p%tr) .or. size(jv) /= size(direction)) then
            jv = 0.0_wp
            call fatal_error(error, "3D UV Jacobian vectors must have ngrid nv entries")
            return
         end if
         call scale_flat(p%dcdt, direction, p%dc, p%ctx%get_num_threads())
         call p%equation%jvp(p%dc, p%dt, error)
         if (allocated(error)) then
            jv = 0.0_wp
            return
         end if
         call subtract_flat(direction, p%dt, jv, p%ctx%get_num_threads(), feasible)
         if (.not. feasible) call fatal_error(error, "3D UV Jacobian product is not finite")
      end associate
   end subroutine evaluator_jvp

   !> Transpose products are available once bound
   !>
   !> @param[in] self  evaluator
   pure logical function evaluator_has_adjoint(self) result(has)
      !> Evaluator
      class(uv3d_evaluator_type), intent(in) :: self
      has = associated(self%problem)
   end function evaluator_has_adjoint

   !> Exact transpose J^T v = v - c' T'^T v at the last linearization
   !>
   !> @param[in]  self       evaluator with a linearization
   !> @param[in]  direction  v, flat (ngrid nv)
   !> @param[out] jtv        J^T v, flat (ngrid nv)
   !> @param[out] feasible   whether the product is usable
   !> @param[out] error      no linearization, transform failure or non-finite product
   subroutine evaluator_adjoint(self, direction, jtv, feasible, error)
      !> Evaluator
      class(uv3d_evaluator_type), intent(in) :: self
      !> Direction
      real(wp), intent(in) :: direction(:)
      !> Transpose product
      real(wp), intent(out) :: jtv(:)
      !> Whether the product is usable
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      feasible = .false.
      if (.not. associated(self%problem)) then
         jtv = 0.0_wp
         call fatal_error(error, "3D UV evaluator is not bound to a problem")
         return
      end if
      associate (p => self%problem)
         if (.not. p%linearized) then
            jtv = 0.0_wp
            call fatal_error(error, "3D UV transpose requested without a linearization")
            return
         end if
         if (size(direction) /= size(p%tr) .or. size(jtv) /= size(direction)) then
            jtv = 0.0_wp
            call fatal_error(error, "3D UV transpose vectors must have ngrid nv entries")
            return
         end if
         call unflatten(direction, p%dt, p%ctx%get_num_threads())
         call p%equation%adjoint(p%dt, p%dc, error)
         if (allocated(error)) then
            jtv = 0.0_wp
            return
         end if
         call subtract_flat(direction, p%dc, jtv, p%ctx%get_num_threads(), feasible, p%dcdt)
         if (.not. feasible) call fatal_error(error, "3D UV transpose product is not finite")
      end associate
   end subroutine evaluator_adjoint

end module moist_model_moz_3d_driver
