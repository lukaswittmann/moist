!> UV transpose, host response and geometry-gradient finite differences
module test_moz_3d_response
   use, intrinsic :: iso_fortran_env, only: output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid, molecular_grid_set_kgrid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   use moist_model_moz_3d_equation, only: uv3d_equation_type
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_model_moz_3d_driver, only: uv3d_problem_type, uv3d_evaluator_type, new_uv3d_problem
   use moist_model_moz_solver, only: moist_model_moz_solver_statistics_type
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_model_moz_solvent_vv, only: solvent_vv_type
   use moist_model_moz_potential, only: lj_cspce, coulomb_spce, monopole_type, ec_charges_type, potential_sites_type, &
      & custom_lj_type, new_custom_lj, fixed_charges_type, new_fixed_charges
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type, response_item_type, atomic_charge_adjoint_response_type, &
      & potential_adjoint_response_type
   use test_moz_fixtures, only: new_cspce_water_solvent
   use test_helpers, only: check_moist_error, get_qc_handymod_recipe
   implicit none(type, external)
   private
   public :: collect_moz_3d_response

contains

   subroutine collect_moz_3d_response(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("equation_cartesian", equation_cartesian), &
         & new_unittest("equation_molecular", equation_molecular), &
         & new_unittest("tables_cartesian", tables_cartesian), &
         & new_unittest("tables_molecular", tables_molecular), &
         & new_unittest("tables_weighted_cartesian", tables_weighted_cartesian), &
         & new_unittest("tables_weighted_molecular", tables_weighted_molecular), &
         & new_unittest("host_charges", host_charges), &
         & new_unittest("host_field", host_field), &
         & new_unittest("fixed_terms", fixed_terms), &
         & new_unittest("scf_tables_cartesian", scf_tables_cartesian), &
         & new_unittest("scf_tables_molecular", scf_tables_molecular), &
         & new_unittest("scf_fallback", scf_fallback), &
         & new_unittest("cold_guesses", cold_guesses), &
         & new_unittest("scf_host", scf_host), &
         & new_unittest("geometry_cartesian", geometry_cartesian), &
         & new_unittest("geometry_molecular", geometry_molecular), &
         & new_unittest("gradient_cartesian", gradient_cartesian), &
         & new_unittest("gradient_charges", gradient_charges), &
         & new_unittest("gradient_molecular", gradient_molecular), &
         & new_unittest("gradient_ec", gradient_ec), &
         & new_unittest("gradient_drism_cartesian", gradient_drism_cartesian), &
         & new_unittest("gradient_drism_charges", gradient_drism_charges), &
         & new_unittest("gradient_drism_molecular", gradient_drism_molecular), &
         & new_unittest("gradient_drism_ec", gradient_drism_ec)]
   end subroutine collect_moz_3d_response

   subroutine make_grid(kind, grid, err)
      integer, intent(in) :: kind
      class(moist_math_grid_3d_type), allocatable, intent(out) :: grid
      type(moist_error), allocatable, intent(out) :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(moist_math_grid_3d_molecular_type) :: mole
      type(moist_math_grid_atomic_recipe_type) :: recipe
      integer :: i

      call new(mol, [1, 8], reshape([-0.4_wp, 0.2_wp, 0.1_wp, 0.9_wp, -0.1_wp, 0.3_wp], [3, 2]))
      if (kind == 1 .or. kind == 4) then
         call new_cartesian_point_grid(cart, err, 9, 8, 7, 0.9_wp, margin=0.0_wp)
         if (.not. allocated(err)) call cart%update(mol, err)
         if (.not. allocated(err)) allocate (grid, source=cart)
      else
         call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=4.0_wp)
         if (.not. allocated(err)) call new_molecular_point_grid(mole, err, recipe=recipe, &
            & dr=1.1_wp, nufft_tol=1.0e-12_wp)
         if (.not. allocated(err)) call mole%update(mol, err)
         if (.not. allocated(err)) allocate (grid, source=mole)
      end if
      if (allocated(err)) return
      grid%nthreads = 2
      if (kind >= 3) then
         do i = 1, grid%ngrid
            grid%w(i) = grid%w(i)*(1.0_wp + 0.125_wp*real(mod(i, 5), wp))
         end do
      end if
   end subroutine make_grid

   !> Small nonsymmetric, nonconstant solvent response, including a stable tail.
   !> This isolates reverse algebra from bulk convergence and molecular resolution.
   subroutine make_data(grid, response, sr, lr, uk)
      class(moist_math_grid_3d_type), intent(in) :: grid
      type(uv3d_susceptibility_type), intent(out) :: response
      real(wp), allocatable, intent(out) :: sr(:, :), lr(:, :)
      complex(wp), allocatable, intent(out) :: uk(:, :)
      real(wp) :: r(3), k(3), k2
      integer :: i, col

      response%nv = 2
      allocate (response%values(2, 2, grid%npts_k), response%index(grid%npts_k), response%radius(grid%npts_k))
      allocate (response%charge_response(2, grid%npts_k))
      response%charge = [-0.7_wp, 0.7_wp]
      allocate (sr(grid%ngrid, 2), lr(grid%ngrid, 2), uk(grid%npts_k, 2))
      do i = 1, grid%npts_k
         k = grid%kpoint(i)
         k2 = sum(k**2)
         response%index(i) = i
         response%radius(i) = sqrt(k2)
         response%values(:, :, i) = reshape([0.04_wp, -0.02_wp, 0.03_wp, 0.06_wp], [2, 2])*exp(-k2)
         response%values(1, 1, i) = response%values(1, 1, i) + 1.0_wp
         response%values(2, 2, i) = response%values(2, 2, i) + 1.0_wp
         response%charge_response(:, i) = [0.13_wp, -0.09_wp]/(1.0_wp + k2)
         do col = 1, 2
            uk(i, col) = cmplx(0.03_wp*col, 0.02_wp*k(1), wp)*exp(-k2)
         end do
      end do
      do i = 1, grid%ngrid
         r = grid%point(i)
         do col = 1, 2
            sr(i, col) = 0.08_wp*cos(r(1) + col)*exp(-sum(r**2)/8.0_wp)
            lr(i, col) = 0.03_wp*sin(r(2) + col)/(1.0_wp + sum(r**2))
         end do
      end do
   end subroutine make_data

   subroutine equation_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      call equation_case(1, error)
   end subroutine equation_cartesian

   subroutine equation_molecular(error)
      type(error_type), allocatable, intent(out) :: error
      call equation_case(2, error)
   end subroutine equation_molecular

   !> Equation dot identities, complex tail finite differences and J^T pairing.
   subroutine equation_case(kind, error)
      integer, intent(in) :: kind
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response, regular
      type(uv3d_equation_type) :: eq, plus, minus
      type(uv3d_problem_type), target :: problem
      type(uv3d_evaluator_type) :: evaluator
      type(model_moz_3d_parameters_type) :: params
      real(wp), allocatable :: sr(:, :), lr(:, :), v(:, :), seed(:, :), jv(:, :), wc(:, :), tp(:, :), tm(:, :)
      real(wp), allocatable :: flat(:), other(:), out(:), adj(:)
      complex(wp), allocatable :: uk(:, :), du(:, :), wk(:, :)
      real(wp) :: lhs, rhs
      real(wp), parameter :: eps = 1.0e-4_wp, beta = 1.7_wp
      logical :: feasible
      integer :: route, i, col
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_grid(kind, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_data(grid, response, sr, lr, uk)
      regular = response
      deallocate (regular%charge_response)
      allocate (v, seed, jv, wc, tp, tm, mold=sr)
      allocate (du, mold=uk)
      do col = 1, 2
         do i = 1, grid%ngrid
            v(i, col) = sin(0.17_wp*i + col)
            seed(i, col) = cos(0.23_wp*i*col)
         end do
         do i = 1, grid%npts_k
            du(i, col) = cmplx(cos(0.31_wp*i*col), sin(0.19_wp*i + col), wp)
         end do
      end do
      do route = 1, 2
         if (route == 1) then
            call eq%setup(grid, regular, uk, beta, err)
            if (.not. allocated(err)) call plus%setup(grid, regular, uk + eps*du, beta, err)
            if (.not. allocated(err)) call minus%setup(grid, regular, uk - eps*du, beta, err)
         else
            call eq%setup(grid, response, uk, beta, err)
            if (.not. allocated(err)) call plus%setup(grid, response, uk + eps*du, beta, err)
            if (.not. allocated(err)) call minus%setup(grid, response, uk - eps*du, beta, err)
         end if
         if (.not. allocated(err)) call eq%jvp(v, jv, err)
         if (.not. allocated(err)) call eq%adjoint(seed, wc, err, wk=wk)
         if (.not. allocated(err)) call plus%map(sr, tp, err)
         if (.not. allocated(err)) call minus%map(sr, tm, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         lhs = sum(seed*jv)
         rhs = sum(wc*v)
         call check(error, lhs, rhs, thr=2.0e-9_wp*max(1.0_wp, abs(lhs)), more="equation transpose")
         if (allocated(error)) return
         lhs = sum(seed*(tp - tm))/(2.0_wp*eps)
         rhs = real(sum(conjg(wk)*du), wp)
         call check(error, lhs, rhs, thr=2.0e-8_wp*max(1.0_wp, abs(lhs)), more="complex tail transpose")
         if (allocated(error)) return
      end do
      call new_uv3d_problem(problem, ctx, grid, response, sr, lr, uk, [0.01_wp, 0.014_wp], beta, params, err)
      if (.not. allocated(err)) call problem%equation%map(problem%cr, problem%tr, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      flat = reshape(v, [size(v)])
      other = reshape(seed, [size(seed)])
      allocate (out, adj, mold=flat)
      out = 1.0_wp
      call evaluator%jvp(flat, out, feasible, err)
      call check(error, allocated(err) .and. .not. feasible .and. all(out == 0.0_wp), &
         & more="an unbound Jacobian product must return a zero field")
      if (allocated(error)) return
      deallocate (err)
      evaluator%problem => problem
      adj = 1.0_wp
      call evaluator%adjoint(other, adj, feasible, err)
      call check(error, allocated(err) .and. .not. feasible .and. all(adj == 0.0_wp), &
         & more="a transpose without a linearization must return a zero field")
      if (allocated(error)) return
      deallocate (err)
      if (.not. allocated(err)) then
         call evaluator%linearize(problem%tr, feasible, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, feasible, more="initial UV linearization is infeasible")
      if (allocated(error)) return
      out = 1.0_wp
      call evaluator%jvp(flat(:size(flat) - 1), out, feasible, err)
      call check(error, allocated(err) .and. .not. feasible .and. all(out == 0.0_wp), &
         & more="a Jacobian product with a wrong size must return a zero field")
      if (allocated(error)) return
      deallocate (err)
      call evaluator%jvp(flat, out, feasible, err)
      if (.not. allocated(err)) call evaluator%adjoint(other, adj, feasible, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      lhs = sum(other*out)
      rhs = sum(flat*adj)
      call check(error, lhs, rhs, thr=2.0e-9_wp*max(1.0_wp, abs(lhs)), more="closure and equation J transpose")
      call ctx%delete()
   end subroutine equation_case

   subroutine tables_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      call tables_case(1, error)
   end subroutine tables_cartesian

   subroutine tables_molecular(error)
      type(error_type), allocatable, intent(out) :: error
      call tables_case(2, error)
   end subroutine tables_molecular

   !> Check response consistency with nonuniform stored Cartesian weights
   !>
   !> @param[out] error Test failure
   subroutine tables_weighted_cartesian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call tables_case(4, error)
   end subroutine tables_weighted_cartesian

   !> Check response consistency with modified stored molecular weights
   !>
   !> @param[out] error Test failure
   subroutine tables_weighted_molecular(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call tables_case(3, error)
   end subroutine tables_weighted_molecular

   !> All closures, all tables (real and imaginary k seeds), varied energy and
   !> moment seeds, warm starts and a zero objective, against converged energies.
   subroutine tables_case(kind, error)
      integer, intent(in) :: kind
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type) :: problem
      type(model_moz_3d_parameters_type) :: params
      character(len=4), parameter :: closures(5) = ["HNC ", "KH  ", "PSE2", "PSE3", "PSE4"]
      real(wp), parameter :: shifts(4) = [0.0_wp, 0.7_wp, 0.4_wp, 0.7_wp]
      real(wp), parameter :: weights(4) = [1.0_wp, 1.0_wp, 2.5_wp, 0.0_wp]
      real(wp), parameter :: beta = 1.7_wp
      real(wp), allocatable :: sr(:, :), lr(:, :), wr(:, :), wl(:, :), lambda(:, :), direction(:, :)
      real(wp), allocatable :: vs(:, :), vl(:, :), old_wr(:, :), old_wl(:, :)
      complex(wp), allocatable :: uk(:, :), wk(:, :), dk(:, :), vk(:, :), old_wk(:, :)
      real(wp) :: gp, gm, dp, dm, fd, dot, max_error, eps
      integer :: ic, is, itab, i, col
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_grid(kind, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_data(grid, response, sr, lr, uk)
      ! Keep the weighted Cartesian KH stencil local to its closure branch
      eps = 1.0e-4_wp
      if (kind == 4) eps = 1.0e-5_wp
      allocate (direction, mold=sr)
      allocate (dk, mold=uk)
      do col = 1, 2
         do i = 1, grid%ngrid
            direction(i, col) = sin(0.17_wp*i + col)
         end do
         do i = 1, grid%npts_k
            dk(i, col) = cmplx(cos(0.23_wp*i + col), sin(0.31_wp*i*col), wp)
         end do
      end do
      params%tolerance = 1.0e-12_wp
      max_error = 0.0_wp
      params%solver%jfnk%restart = 12
      params%solver%jfnk%adjoint%restart = 9
      params%solver%jfnk%adjoint%max_krylov = 50
      do ic = 1, size(closures)
         params%closure = trim(closures(ic))
         do is = 1, size(shifts)
            call solve_terms(sr, lr, uk, gp, dp)
            if (allocated(error)) return
            call problem%response(shifts(is), 1.0e-10_wp, wr, wl, wk, lambda, err, free_energy_weight=weights(is))
            if (allocated(err)) then
               write (*, "(a,a,a,i0)") "  response failed for ", trim(closures(ic)), ", seed=", is
            end if
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, size(problem%solver%workspace%ws_v, 2) == 13, &
               & more="3D response ignored its configured restart")
            if (allocated(error)) return
            do itab = 1, 4
               vs = sr
               vl = lr
               vk = uk
               select case (itab)
               case (1)
                  vs = sr + eps*direction
                  dot = sum(wr*direction)
               case (2)
                  vl = lr + eps*direction
                  dot = sum(wl*direction)
               case (3)
                  vk = uk + eps*real(dk, wp)
                  dot = real(sum(conjg(wk)*real(dk, wp)), wp)
               case default
                  vk = uk + cmplx(0.0_wp, eps, wp)*aimag(dk)
                  dot = real(sum(conjg(wk)*cmplx(0.0_wp, 1.0_wp, wp)*aimag(dk)), wp)
               end select
               call solve_terms(vs, vl, vk, gp, dp)
               if (allocated(error)) return
               select case (itab)
               case (1)
                  vs = sr - eps*direction
               case (2)
                  vl = lr - eps*direction
               case (3)
                  vk = uk - eps*real(dk, wp)
               case default
                  vk = uk - cmplx(0.0_wp, eps, wp)*aimag(dk)
               end select
               call solve_terms(vs, vl, vk, gm, dm)
               if (allocated(error)) return
               fd = (weights(is)*(gp - gm) + shifts(is)*(dp - dm)/beta)/(2.0_wp*eps)
               max_error = max(max_error, abs(fd - dot)/max(1.0e-3_wp, abs(fd)))
               call check(error, dot, fd, thr=2.0e-6_wp*max(1.0e-3_wp, abs(fd)), &
                  & more=trim(closures(ic))//" potential-table energy response")
               if (allocated(error)) return
            end do
         end do
      end do
      call solve_terms(sr, lr, uk, gp, dp)
      if (allocated(error)) return
      old_wr = wr
      old_wl = wl
      old_wk = wk
      call problem%response(shifts(4), 1.0e-10_wp, wr, wl, wk, lambda, err, free_energy_weight=weights(4))
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(wr - old_wr)) + maxval(abs(wl - old_wl)) + maxval(abs(wk - old_wk)) &
         & <= 1.0e-8_wp*max(1.0e-3_wp, maxval(abs(old_wr))), more="warm-started response differs")
      if (allocated(error)) return
      call problem%response(0.0_wp, 1.0e-10_wp, wr, wl, wk, lambda, err, free_energy_weight=0.0_wp)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(wr)) + maxval(abs(wl)) + maxval(abs(wk)) + maxval(abs(lambda)), 0.0_wp, &
         & thr=1.0e-25_wp, more="zero objective has nonzero adjoints")
      write (*, "(a,i0,a,es12.4)") "  3D table FD grid=", kind, " max scaled error=", max_error
      call ctx%delete()

   contains

      subroutine solve_terms(t_sr, t_lr, t_k, g, d)
         real(wp), intent(in) :: t_sr(:, :), t_lr(:, :)
         complex(wp), intent(in) :: t_k(:, :)
         real(wp), intent(out) :: g, d
         logical :: converged
         real(wp) :: rmse
         integer :: niter
         call new_uv3d_problem(problem, ctx, grid, response, t_sr, t_lr, t_k, [0.01_wp, 0.014_wp], &
            & beta, params, err)
         if (.not. allocated(err)) call problem%solve(params%tolerance, params%maxiter, converged, rmse, niter, err)
         if (.not. allocated(err)) then
            if (.not. converged) then
               call test_failed(error, "3D UV solve of the finite-difference test did not converge")
               return
            end if
            call problem%free_energy(g, err)
         end if
         if (.not. allocated(err)) call problem%direct_moment(d, err)
         call check_moist_error(error, err)
      end subroutine solve_terms

   end subroutine tables_case

   subroutine scf_tables_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      call scf_tables_case(1, error)
   end subroutine scf_tables_cartesian

   subroutine scf_tables_molecular(error)
      type(error_type), allocatable, intent(out) :: error
      call scf_tables_case(2, error)
   end subroutine scf_tables_molecular

   !> Changing both parts of the tail preserves physical t and converges to cold references.
   !> All closures and solvers, fixed Cartesian/molecular grids and retained JFNK storage are exercised.
   subroutine scf_tables_case(kind, error)
      integer, intent(in) :: kind
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type) :: warm, cold
      type(model_moz_3d_parameters_type) :: params
      type(moist_model_moz_solver_statistics_type) :: before, after
      character(len=4), parameter :: closures(9) = ["HNC ", "KH  ", "PSE2", "PSE3", "PSE4", "KH  ", "KH  ", "KH  ", &
         & "KH  "]
      character(len=6), parameter :: solvers(9) = ["JFNK  ", "JFNK  ", "JFNK  ", "JFNK  ", "JFNK  ", &
         & "PICARD", "DIIS  ", "HYBRID", "JFNK  "]
      real(wp), allocatable :: sr(:, :), lr(:, :), changed_sr(:, :), changed_lr(:, :)
      real(wp), allocatable :: previous_t(:, :), previous_lr(:, :)
      complex(wp), allocatable :: uk(:, :), changed_k(:, :)
      real(wp), parameter :: beta = 1.7_wp
      real(wp) :: rmse, ew, ec
      logical :: converged, kept
      integer :: ic, step, nw, nc, warm_steps, cold_steps
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_grid(kind, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_data(grid, response, sr, lr, uk)
      params%tolerance = 1.0e-12_wp
      params%solver%jfnk%restart = 8
      params%solver%picard%damping = 0.5_wp
      warm_steps = 0
      cold_steps = 0
      do ic = 1, size(closures)
         params%closure = closures(ic)
         params%solver%name = solvers(ic)
         ! The last case is exact Newton: exact products and a tight fixed forcing term
         if (ic == size(closures)) then
            params%solver%jfnk%use_exact_jvp = .true.
            params%solver%jfnk%eta0 = 1.0e-8_wp
            params%solver%jfnk%eta_min = 1.0e-8_wp
            params%solver%jfnk%eta_max = 1.0e-8_wp
            params%solver%jfnk%alpha_max = 1.0_wp
         end if
         call new_uv3d_problem(warm, ctx, grid, response, sr, lr, uk, [0.01_wp, 0.014_wp], beta, params, err)
         if (.not. allocated(err)) call warm%solve(params%tolerance, params%maxiter, converged, rmse, nw, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged, more="initial SCF table solve did not converge")
         if (allocated(error)) return
         previous_lr = lr
         do step = 1, 3
            changed_sr = sr*(1.0_wp + 0.01_wp*step)
            changed_lr = lr*(1.0_wp + 0.02_wp*step)
            changed_k = uk*(1.0_wp + 0.02_wp*step)
            previous_t = warm%tr
            before = warm%solver%get_statistics()
            if (params%solver%name == "JFNK") then
               ! Sentinel in deflation scratch: the undeflated 3D solve never touches it and
               ! a workspace reallocation drops it; a deflated solve would resize it and fail
               call check(error, allocated(warm%solver%workspace), more="JFNK storage released after solve")
               if (allocated(error)) return
               if (.not. allocated(warm%solver%workspace%ws_z)) allocate (warm%solver%workspace%ws_z(1))
               warm%solver%workspace%ws_z(1) = 7.0_wp
            end if
            call warm%update_potential(changed_sr, changed_lr, changed_k, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, maxval(abs(warm%tr - previous_t - beta*(previous_lr - changed_lr))), &
               & 0.0_wp, thr=1.0e-14_wp, more="Ng update did not preserve physical t")
            if (allocated(error)) return
            call warm%solve(params%tolerance, params%maxiter, converged, rmse, nw, err, reuse_initial_guess=.true.)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, converged .and. rmse < params%tolerance, more="warm SCF table solve did not converge")
            if (allocated(error)) return
            after = warm%solver%get_statistics()
            if (params%solver%name == "JFNK") then
               kept = .false.
               if (allocated(warm%solver%workspace)) then
                  if (allocated(warm%solver%workspace%ws_z)) kept = warm%solver%workspace%ws_z(1) == 7.0_wp
               end if
               call check(error, after%nresidual > before%nresidual .and. kept, &
                  & more="SCF update rebuilt or reset JFNK storage")
               if (allocated(error)) return
            end if
            call new_uv3d_problem(cold, ctx, grid, response, changed_sr, changed_lr, changed_k, &
               & [0.01_wp, 0.014_wp], beta, params, err)
            if (.not. allocated(err)) call cold%solve(params%tolerance, params%maxiter, converged, rmse, nc, err)
            if (.not. allocated(err)) call warm%free_energy(ew, err)
            if (.not. allocated(err)) call cold%free_energy(ec, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, converged, more="cold SCF reference did not converge")
            if (allocated(error)) return
            call check(error, ew, ec, thr=5.0e-12_wp, more="warm and cold SCF energies differ")
            if (allocated(error)) return
            call check(error, maxval(abs(warm%tr - cold%tr)), 0.0_wp, thr=5.0e-10_wp, more="SCF t fields differ")
            if (allocated(error)) return
            call check(error, maxval(abs(warm%cr - cold%cr)), 0.0_wp, thr=5.0e-10_wp, more="SCF c fields differ")
            if (allocated(error)) return
            warm_steps = warm_steps + nw
            cold_steps = cold_steps + nc
            previous_lr = changed_lr
         end do
      end do
      write (*, "(a,2i6)") "  SCF table outer iterations warm/cold=", warm_steps, cold_steps
      call ctx%delete()
   end subroutine scf_tables_case

   !> Invalid guesses, a failed feasible warm attempt, rejected updates and recovery.
   subroutine scf_fallback(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type) :: problem
      type(model_moz_3d_parameters_type) :: params
      real(wp), allocatable :: sr(:, :), lr(:, :), old_t(:, :), bad_sr(:, :)
      complex(wp), allocatable :: uk(:, :)
      real(wp) :: rmse
      logical :: converged
      integer :: i, niter
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_grid(1, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_data(grid, response, sr, lr, uk)
      ! Identity chi makes T independent of c. The cold field has exactly zero residual.
      response%values = 0.0_wp
      do i = 1, size(response%values, 3)
         response%values(1, 1, i) = 1.0_wp
         response%values(2, 2, i) = 1.0_wp
      end do
      params%solver%name = "PICARD"
      call new_uv3d_problem(problem, ctx, grid, response, sr, lr, uk, [0.01_wp, 0.014_wp], 1.7_wp, params, err)
      if (.not. allocated(err)) call problem%solve(1.0e-12_wp, 2, converged, rmse, niter, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged, more="fallback reference did not converge")
      if (allocated(error)) return
      old_t = problem%tr
      call problem%update_potential(sr, lr, uk(:, :1), err)
      call check(error, allocated(err) .and. all(problem%tr == old_t), more="invalid update damaged the kept fields")
      if (allocated(error)) return
      deallocate (err)
      problem%tr(1, 1) = ieee_value(0.0_wp, ieee_quiet_nan)
      call problem%solve(1.0e-12_wp, 2, converged, rmse, niter, err, reuse_initial_guess=.true.)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == 1 .and. all(problem%tr == old_t), more="nonfinite guess did not fall back")
      if (allocated(error)) return
      problem%tr = old_t + 1000.0_wp
      call problem%solve(1.0e-12_wp, 2, converged, rmse, niter, err, reuse_initial_guess=.true.)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == 3 .and. all(problem%tr == old_t), more="failed warm solve did not retry cold")
      if (allocated(error)) return
      bad_sr = sr
      bad_sr(1, 1) = ieee_value(0.0_wp, ieee_quiet_nan)
      call problem%update_potential(bad_sr, lr, uk, err)
      if (.not. allocated(err)) call problem%solve(1.0e-12_wp, 2, converged, rmse, niter, err, reuse_initial_guess=.true.)
      call check(error, allocated(err) .and. .not. converged, more="invalid cold and warm potentials were accepted")
      if (allocated(error)) return
      deallocate (err)
      call problem%update_potential(sr, lr, uk, err)
      if (.not. allocated(err)) call problem%solve(1.0e-12_wp, 2, converged, rmse, niter, err, reuse_initial_guess=.true.)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == 1 .and. all(problem%tr == old_t), more="failed solve poisoned recovery")
      call ctx%delete()
   end subroutine scf_fallback

   !> Both cold starts solve the same equations and retain warm-start semantics
   subroutine cold_guesses(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type) :: cs0, ts0
      type(model_moz_3d_parameters_type) :: params
      character(len=4), parameter :: closures(4) = ["KH  ", "HNC ", "PSE3", "KH  "]
      character(len=6), parameter :: solvers(4) = ["JFNK  ", "JFNK  ", "DIIS  ", "PICARD"]
      real(wp), allocatable :: sr(:, :), lr(:, :), wcs(:, :), lcs(:, :), acs(:, :)
      real(wp), allocatable :: wts(:, :), lts(:, :), ats(:, :)
      complex(wp), allocatable :: uk(:, :), kcs(:, :), kts(:, :)
      real(wp) :: rmse, ecs, ets, dcs, dts
      integer :: kind, ic, niter
      logical :: converged
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=2, verbosity=0)
      params%solver%picard%damping = 0.5_wp
      params%solver%jfnk%restart = 8
      do kind = 1, 2
         call make_grid(kind, grid, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call make_data(grid, response, sr, lr, uk)
         do ic = 1, size(closures)
            params%closure = closures(ic)
            params%solver%name = solvers(ic)
            params%cold_guess = "cs0"
            call new_uv3d_problem(cs0, ctx, grid, response, sr, lr, uk, [0.01_wp, 0.014_wp], 1.7_wp, params, err)
            if (.not. allocated(err)) call cs0%solve(1.0e-12_wp, 128, converged, rmse, niter, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, converged, more="cs0 reference did not converge")
            if (allocated(error)) return
            params%cold_guess = "TS0"
            call new_uv3d_problem(ts0, ctx, grid, response, sr, lr, uk, [0.01_wp, 0.014_wp], 1.7_wp, params, err)
            if (.not. allocated(err)) call ts0%solve(1.0e-12_wp, 128, converged, rmse, niter, err)
            if (.not. allocated(err)) call cs0%free_energy(ecs, err)
            if (.not. allocated(err)) call ts0%free_energy(ets, err)
            if (.not. allocated(err)) call cs0%direct_moment(dcs, err)
            if (.not. allocated(err)) call ts0%direct_moment(dts, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, converged, more="ts0 solve did not converge")
            if (allocated(error)) return
            call check(error, ecs, ets, thr=5.0e-12_wp, more="cold-guess energies differ")
            if (allocated(error)) return
            call check(error, dcs, dts, thr=5.0e-12_wp, more="cold-guess direct moments differ")
            if (allocated(error)) return
            call check(error, maxval(abs(cs0%tr - ts0%tr)), 0.0_wp, thr=5.0e-10_wp, more="cold-guess t fields differ")
            if (allocated(error)) return
            call check(error, maxval(abs(cs0%cr - ts0%cr)), 0.0_wp, thr=5.0e-10_wp, more="cold-guess c fields differ")
            if (allocated(error)) return
            if (ic == 1) then
               if (allocated(acs)) deallocate (acs, ats)
               call cs0%response(0.04_wp, 1.0e-11_wp, wcs, lcs, kcs, acs, err)
               if (.not. allocated(err)) call ts0%response(0.04_wp, 1.0e-11_wp, wts, lts, kts, ats, err)
               call check_moist_error(error, err)
               if (allocated(error)) return
               call check(error, maxval(abs(wcs - wts)), 0.0_wp, thr=5.0e-10_wp, more="cold-guess SR responses differ")
               if (allocated(error)) return
               call check(error, maxval(abs(lcs - lts)), 0.0_wp, thr=5.0e-10_wp, more="cold-guess LR responses differ")
               if (allocated(error)) return
               call check(error, maxval(abs(kcs - kts)), 0.0_wp, thr=5.0e-10_wp, more="cold-guess k responses differ")
               if (allocated(error)) return
            end if
            ! With only one step, a valid retained solution must take precedence over ts0.
            call ts0%solve(1.0e-12_wp, 1, converged, rmse, niter, err, reuse_initial_guess=.true.)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, converged .and. niter == 1, more="ts0 displaced a valid warm guess")
            if (allocated(error)) return
         end do
      end do
      params%cold_guess = "unknown"
      call new_uv3d_problem(ts0, ctx, grid, response, sr, lr, uk, [0.01_wp, 0.014_wp], 1.7_wp, params, err)
      call check(error, allocated(err), more="standalone constructor accepted an unknown cold guess")
      if (allocated(error)) return
      call ctx%delete()
      call ts0_fallback(error)
   end subroutine cold_guesses

   !> A failed warm attempt and an infeasible guess both recover through ts0
   subroutine ts0_fallback(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(moist_math_grid_3d_type), allocatable :: grid
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type) :: problem
      type(model_moz_3d_parameters_type) :: params
      real(wp), allocatable :: sr(:, :), lr(:, :), reference(:, :)
      complex(wp), allocatable :: uk(:, :)
      real(wp) :: rmse
      integer :: i, niter, cold_steps
      logical :: converged
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_grid(1, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_data(grid, response, sr, lr, uk)
      response%values = 0.0_wp
      do i = 1, size(response%values, 3)
         response%values(1, 1, i) = 1.0_wp
         response%values(2, 2, i) = 1.0_wp
      end do
      params%cold_guess = "ts0"
      params%solver%name = "PICARD"
      params%solver%picard%damping = 0.5_wp
      call new_uv3d_problem(problem, ctx, grid, response, sr, lr, uk, [0.01_wp, 0.014_wp], 1.7_wp, params, err)
      if (.not. allocated(err)) call problem%solve(1.0e-8_wp, 32, converged, rmse, cold_steps, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. cold_steps > 1, more="ts0 fallback reference did not converge")
      if (allocated(error)) return
      reference = problem%tr
      problem%tr = reference + 1.0e8_wp
      call problem%solve(1.0e-8_wp, 32, converged, rmse, niter, err, reuse_initial_guess=.true.)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == 32 + cold_steps, more="failed warm solve did not retry the selected ts0 path")
      if (allocated(error)) return
      call check(error, maxval(abs(problem%tr - reference)), 0.0_wp, thr=1.0e-13_wp, more="ts0 retry changed the cold solution")
      if (allocated(error)) return
      problem%tr(1, 1) = ieee_value(0.0_wp, ieee_quiet_nan)
      call problem%solve(1.0e-8_wp, 32, converged, rmse, niter, err, reuse_initial_guess=.true.)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == cold_steps, more="infeasible warm guess did not start from ts0")
      call ctx%delete()
   end subroutine ts0_fallback

   !> Real EC fields at fixed geometry: warm/cold energies, fields and response, then an update.
   subroutine scf_host(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_3d_type), target :: warm, cold
      type(potential_sites_type) :: sites
      type(structure_type) :: mol
      real(wp), allocatable :: phi(:), field(:), direction(:), dw(:), dc(:)
      real(wp), allocatable :: cw(:, :), tw(:, :), hw(:, :), cc(:, :), tc(:, :), hc(:, :)
      real(wp) :: r(3), ew, ec, repeated
      integer :: step, i, a

      call new_context(ctx, nthreads=2, verbosity=0)
      call host_model(ctx, "PC+", 2, warm, solvent, err, reuse_initial_guess=.true.)
      if (.not. allocated(err)) call host_model(ctx, "PC+", 2, cold, solvent, err, reuse_initial_guess=.false.)
      if (.not. allocated(err)) call solvent%potential%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = solvent%structure()
      allocate (phi(warm%grid%ngrid), direction(warm%grid%ngrid))
      do i = 1, size(phi)
         r = warm%grid%point(i)
         phi(i) = 0.0_wp
         do a = 1, mol%nat
            phi(i) = phi(i) + sites%q(a)/max(0.1_wp, norm2(r - mol%xyz(:, a)))
         end do
         direction(i) = 0.1_wp*sin(0.4_wp*r(1))*exp(-sum(r**2)/40.0_wp)
      end do
      do step = 0, 3
         field = phi + 0.001_wp*step*direction
         call host_energy(warm, field, .true., ew, err)
         if (.not. allocated(err)) call host_energy(cold, field, .true., ec, err)
         if (.not. allocated(err)) call warm%correlations(cw, tw, hw, err)
         if (.not. allocated(err)) call cold%correlations(cc, tc, hc, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, ew, ec, thr=1.0e-9_wp, more="EC SCF warm/cold energies differ")
         if (allocated(error)) return
         call check(error, maxval(abs(hw - hc)), 0.0_wp, thr=1.0e-8_wp, more="EC SCF warm/cold correlations differ")
         if (allocated(error)) return
      end do
      call host_response(warm, field, .true., dw, err)
      if (.not. allocated(err)) call host_response(cold, field, .true., dc, err)
      if (.not. allocated(err)) call host_energy(warm, field, .true., repeated, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, repeated == ew, more="EC response discarded the cached UV energy")
      if (allocated(error)) return
      call check(error, maxval(abs(dw - dc)), 0.0_wp, thr=1.0e-8_wp, more="EC SCF warm/cold responses differ")
      if (allocated(error)) return
      mol%xyz(1, 2) = mol%xyz(1, 2) + 0.01_wp
      call warm%update(mol, err)
      if (.not. allocated(err)) call cold%update(mol, err)
      if (.not. allocated(err)) call warm%correlations(cw, tw, hw, err)
      call check(error, allocated(err), more="geometry update kept a current UV solution")
      if (allocated(error)) return
      deallocate (err)
      call host_energy(warm, field, .true., ew, err)
      if (.not. allocated(err)) call host_energy(cold, field, .true., ec, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, ew == ec, more="geometry update reused an old initial guess")
   end subroutine scf_host

   !> Real cSPCE-water solution, with fixed bulk shared across the FD perturbations.
   subroutine host_model(ctx, correction, mode, model, solvent, err, reuse_initial_guess, drism, reuse_bulk)
      type(moist_context_type), intent(in), target :: ctx
      character(len=*), intent(in) :: correction
      integer, intent(in) :: mode
      type(model_moz_3d_type), intent(out) :: model
      type(solvent_vv_type), intent(inout) :: solvent
      type(moist_error), allocatable, intent(out) :: err
      logical, intent(in), optional :: reuse_initial_guess
      logical, intent(in), optional :: drism
      !> Reuse a fixed bulk solution with matching settings across test cases
      logical, intent(in), optional :: reuse_bulk
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(model_moz_3d_parameters_type) :: params
      type(monopole_type) :: charges
      type(ec_charges_type) :: ec
      type(custom_lj_type) :: zero_lj
      type(fixed_charges_type) :: weak_charges
      type(potential_sites_type) :: sites
      logical :: keep_bulk

      params%tolerance = 1.0e-11_wp
      params%solvent%npts = 1024
      params%solvent%dr = 0.06_wp
      params%solvent%tolerance = 1.0e-11_wp
      if (present(drism)) then
         if (drism) params%solvent%theory = "DRISM"
      end if
      params%correction = correction
      if (present(reuse_initial_guess)) params%reuse_initial_guess = reuse_initial_guess
      keep_bulk = .false.
      if (present(reuse_bulk)) keep_bulk = reuse_bulk
      if (.not. keep_bulk .or. .not. solvent%is_solved(params%solvent)) then
         call new_cspce_water_solvent(ctx, solvent, err)
      end if
      if (keep_bulk .and. .not. allocated(err)) then
         if (.not. solvent%is_solved(params%solvent)) call solvent%solve(err, params%solvent)
      end if
      if (.not. allocated(err)) call new_cartesian_point_grid(grid, err, 32, 32, 32, 0.6_wp, margin=8.0_wp)
      if (.not. allocated(err)) call new_moz_3d_model(model, ctx, grid, solvent, err, params)
      if (.not. allocated(err)) then
         if (mode == 3) then
            call solvent%potential%sites(sites, err)
            if (.not. allocated(err)) call new_custom_lj(zero_lj, sites%pair%sigma, 0.0_wp*sites%pair%epsilon, err)
            if (.not. allocated(err)) call model%potential%add(zero_lj, err)
         else
            call model%potential%add(lj_cspce, err)
         end if
      end if
      if (.not. allocated(err)) then
         select case (mode)
         case (1)
            call model%potential%add(charges, err)
         case (2)
            ec%exclusion_radius = 1.5_wp
            ec%charge_penalty = 0.0_wp
            call model%potential%add(ec, err)
         case (3)
            call new_fixed_charges(weak_charges, 0.002_wp*sites%q, err)
            if (.not. allocated(err)) call model%potential%add(weak_charges, err)
         case default
            call model%potential%add(coulomb_spce, err)
         end select
      end if
      if (.not. allocated(err)) call model%update(solvent%structure(), err)
   end subroutine host_model

   !> Answer either atomic charges or a point potential in the active phase.
   subroutine host_energy(model, data, field, energy, err)
      type(model_moz_3d_type), intent(inout), target :: model
      real(wp), intent(in) :: data(:)
      logical, intent(in) :: field
      real(wp), intent(out) :: energy
      type(moist_error), allocatable, intent(out) :: err
      type(coupling_type), pointer :: coupling

      energy = 0.0_wp
      call model%new_coupling(coupling, err)
      if (allocated(err)) return
      call model%prepare_energy(coupling, err)
      if (.not. allocated(err)) then
         if (coupling%next()) then
            if (field) then
               call coupling%answer("phi", data, err)
            else
               call coupling%answer("q", data, err)
            end if
         end if
      end if
      if (.not. allocated(err)) call model%get_energy(coupling, energy, err)
      call model%release_coupling(coupling)
   end subroutine host_energy

   subroutine host_response(model, data, field, derivative, err)
      type(model_moz_3d_type), intent(inout), target :: model
      real(wp), intent(in) :: data(:)
      logical, intent(in) :: field
      real(wp), allocatable, intent(out) :: derivative(:)
      type(moist_error), allocatable, intent(out) :: err
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      class(response_item_type), allocatable :: item

      call model%new_coupling(coupling, err)
      if (allocated(err)) return
      call model%prepare_response(coupling, err)
      if (.not. allocated(err)) then
         if (coupling%next()) then
            if (field) then
               call coupling%answer("phi", data, err)
            else
               call coupling%answer("q", data, err)
            end if
         end if
      end if
      if (.not. allocated(err)) call model%get_response(coupling, response, err)
      call model%release_coupling(coupling)
      if (allocated(err)) return
      do while (response%next())
         item = response%item()
         select type (item)
         type is (atomic_charge_adjoint_response_type)
            derivative = item%dg_dq
         type is (potential_adjoint_response_type)
            derivative = item%w_phi
         end select
      end do
   end subroutine host_response

   !> Host atomic-charge response for none/PC/PC+, including net-charge directions.
   !> Requesting the response before any energy must solve and cache the primal.
   subroutine host_charges(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_3d_type), target :: model
      type(potential_sites_type) :: sites
      character(len=4), parameter :: corrections(3) = ["none", "PC  ", "PC+ "]
      real(wp), parameter :: eps = 1.0e-4_wp
      real(wp), allocatable :: q(:), derivative(:), repeated(:)
      real(wp) :: ep, em, e0, e1, fd, max_error
      integer :: ic, i

      call new_context(ctx, nthreads=2, verbosity=0)
      max_error = 0.0_wp
      do ic = 1, size(corrections)
         call host_model(ctx, trim(corrections(ic)), 1, model, solvent, err)
         if (.not. allocated(err)) call solvent%potential%sites(sites, err)
         if (.not. allocated(err)) call host_response(model, sites%q, .false., derivative, err)
         if (.not. allocated(err)) call host_energy(model, sites%q, .false., e0, err)
         if (.not. allocated(err)) call host_energy(model, sites%q, .false., e1, err)
         if (.not. allocated(err)) call host_response(model, sites%q, .false., repeated, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, allocated(derivative) .and. allocated(repeated), more="missing host charge response")
         if (allocated(error)) return
         call check(error, e0 == e1, more="repeated energy changed after response")
         if (allocated(error)) return
         call check(error, maxval(abs(derivative - repeated)), 0.0_wp, thr=1.0e-9_wp, &
            & more="repeated response differs")
         if (allocated(error)) return
         do i = 1, size(sites%q)
            q = sites%q
            q(i) = q(i) + eps
            call host_energy(model, q, .false., ep, err)
            q(i) = q(i) - 2.0_wp*eps
            if (.not. allocated(err)) call host_energy(model, q, .false., em, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            fd = (ep - em)/(2.0_wp*eps)
            max_error = max(max_error, abs(fd - derivative(i))/max(1.0e-4_wp, abs(fd)))
            call check(error, derivative(i), fd, thr=2.0e-5_wp*max(1.0e-4_wp, abs(fd)), &
               & more=trim(corrections(ic))//" host charge response")
            if (allocated(error)) return
         end do
      end do
      write (*, "(a,es12.4)") "  3D host-charge FD max scaled error=", max_error
   end subroutine host_charges

   !> EC field response includes both direct field and fitted-tail charge paths.
   subroutine host_field(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_3d_type), target :: model
      type(potential_sites_type) :: sites
      type(structure_type) :: mol
      real(wp), parameter :: eps = 1.0e-5_wp
      real(wp), allocatable :: phi(:), direction(:), derivative(:)
      real(wp) :: r(3), ep, em, fd, dot
      integer :: i, a

      call new_context(ctx, nthreads=2, verbosity=0)
      call host_model(ctx, "PC+", 2, model, solvent, err)
      if (.not. allocated(err)) call solvent%potential%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = solvent%structure()
      allocate (phi(model%grid%ngrid), direction(model%grid%ngrid))
      phi = 0.0_wp
      do i = 1, size(phi)
         r = model%grid%point(i)
         do a = 1, mol%nat
            phi(i) = phi(i) + sites%q(a)/max(0.1_wp, norm2(r - mol%xyz(:, a)))
         end do
         direction(i) = 0.1_wp*sin(0.4_wp*r(1))*exp(-sum(r**2)/40.0_wp)
      end do
      call host_response(model, phi, .true., derivative, err)
      if (.not. allocated(err)) call host_energy(model, phi + eps*direction, .true., ep, err)
      if (.not. allocated(err)) call host_energy(model, phi - eps*direction, .true., em, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, allocated(derivative), more="missing host potential adjoint")
      if (allocated(error)) return
      dot = sum(derivative*direction)
      fd = (ep - em)/(2.0_wp*eps)
      call check(error, dot, fd, thr=2.0e-5_wp*max(1.0e-4_wp, abs(fd)), more="EC field and fitted-tail response")
      write (*, "(a,es12.4)") "  3D host-field FD scaled error=", abs(dot - fd)/max(1.0e-4_wp, abs(fd))
   end subroutine host_field

   !> Fixed potentials yield an empty host response despite their RISM adjoint.
   subroutine fixed_terms(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_3d_type), target :: model
      real(wp), allocatable :: derivative(:)

      call new_context(ctx, nthreads=2, verbosity=0)
      call host_model(ctx, "none", 0, model, solvent, err)
      if (.not. allocated(err)) call host_response(model, [0.0_wp], .false., derivative, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. allocated(derivative), more="fixed potential returned a host response")
   end subroutine fixed_terms

   subroutine geometry_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      call geometry_case(1, error)
   end subroutine geometry_cartesian

   subroutine geometry_molecular(error)
      type(error_type), allocatable, intent(out) :: error
      call geometry_case(2, error)
   end subroutine geometry_molecular

   !> Synthetic, converged RISM energies isolate geometry from potential tables.
   !> All closures, both energy/moment seeds and a zero objective are exercised.
   subroutine geometry_case(kind, error)
      integer, intent(in) :: kind
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(moist_math_grid_3d_type), allocatable :: grid, moved
      type(uv3d_susceptibility_type) :: response
      type(uv3d_problem_type) :: problem
      type(model_moz_3d_parameters_type) :: params
      type(volume_adjoint_type) :: acc
      character(len=4), parameter :: closures(5) = ["HNC ", "KH  ", "PSE2", "PSE3", "PSE4"]
      real(wp), parameter :: beta = 1.7_wp, eps = 2.0e-5_wp
      real(wp), allocatable :: sr(:, :), lr(:, :), wr(:, :), wl(:, :), lambda(:, :), dx(:, :), dw(:)
      real(wp), allocatable :: old_wr(:, :), old_wl(:, :)
      complex(wp), allocatable :: uk(:, :), wk(:, :), old_wk(:, :)
      real(wp) :: ep, em, g, d, dot, seed, shift, max_error
      integer :: i, ic, is
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=2, verbosity=0)
      call make_grid(kind, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call make_data(grid, response, sr, lr, uk)
      allocate (dx(3, grid%ngrid), dw(grid%ngrid))
      do i = 1, grid%ngrid
         dx(:, i) = [sin(0.13_wp*i), cos(0.17_wp*i), sin(0.23_wp*i)]
         dw(i) = 0.1_wp*cos(0.11_wp*i)
      end do
      if (kind == 1) then
         ! Cartesian quadrature is fixed by dr, so it has no independent
         ! geometry direction. The model test below moves its centroid.
         dx = 0.0_wp
         dw = 0.0_wp
      end if
      params%tolerance = 1.0e-12_wp
      max_error = 0.0_wp
      do ic = 1, size(closures)
         params%closure = trim(closures(ic))
         do is = 1, 3
            seed = 1.0_wp
            shift = 0.0_wp
            if (is == 2) then
               seed = 2.5_wp
               shift = 0.7_wp
            else if (is == 3) then
               seed = 0.0_wp
            end if
            call solve_geometry(grid, g, d)
            if (allocated(error)) return
            call problem%response(shift, 1.0e-11_wp, old_wr, old_wl, old_wk, lambda, err, free_energy_weight=seed)
            call acc%init(grid%ngrid)
            if (.not. allocated(err)) call problem%response(shift, 1.0e-11_wp, wr, wl, wk, lambda, err, &
               & free_energy_weight=seed, grid_adjoint=acc)
            call check_moist_error(error, err)
            if (allocated(error)) return
            if (is == 3) then
               call check(error, maxval(abs(acc%w_xyz)) + maxval(abs(acc%w_w)) + maxval(abs(acc%w_kref)), &
                  & 0.0_wp, thr=0.0_wp, more="zero objective has geometry sensitivities")
               if (allocated(error)) return
            end if
            call check(error, maxval(abs(wr - old_wr)) + maxval(abs(wl - old_wl)) + maxval(abs(wk - old_wk)), &
               & 0.0_wp, thr=1.0e-10_wp, more="geometry request changed table response")
            if (allocated(error)) return
            dot = sum(acc%w_xyz*dx) + sum(acc%w_w*dw)
            allocate (moved, source=grid)
            moved%xyz = grid%xyz + eps*dx
            moved%w = grid%w + eps*dw
            call solve_geometry(moved, g, d)
            ep = seed*g + shift*d/beta
            moved%xyz = grid%xyz - eps*dx
            moved%w = grid%w - eps*dw
            call solve_geometry(moved, g, d)
            em = seed*g + shift*d/beta
            deallocate (moved)
            if (allocated(error)) return
            d = (ep - em)/(2.0_wp*eps)
            max_error = max(max_error, abs(dot - d)/max(1.0e-5_wp, abs(d)))
            call check(error, dot, d, thr=2.0e-6_wp*max(1.0e-5_wp, abs(d)), &
               & more=trim(closures(ic))//" geometry energy/moment derivative")
            if (allocated(error)) return
         end do
      end do
      write (*, "(a,i0,a,es12.4)") "  3D geometry FD grid=", kind, " max scaled error=", max_error
      call ctx%delete()

   contains

      subroutine solve_geometry(active, energy, moment)
         class(moist_math_grid_3d_type), intent(in) :: active
         real(wp), intent(out) :: energy, moment
         logical :: converged
         real(wp) :: rmse
         integer :: niter
         call new_uv3d_problem(problem, ctx, active, response, sr, lr, uk, [0.01_wp, 0.014_wp], beta, params, err)
         if (.not. allocated(err)) call problem%solve(params%tolerance, params%maxiter, converged, rmse, niter, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged, more="geometry FD primal did not converge")
         if (allocated(error)) return
         call problem%free_energy(energy, err)
         if (.not. allocated(err)) call problem%direct_moment(moment, err)
         call check_moist_error(error, err)
      end subroutine solve_geometry

   end subroutine geometry_case

   subroutine gradient_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      call gradient_case(1, 0, error)
   end subroutine gradient_cartesian

   subroutine gradient_charges(error)
      type(error_type), allocatable, intent(out) :: error
      call gradient_case(1, 1, error)
   end subroutine gradient_charges

   subroutine gradient_molecular(error)
      type(error_type), allocatable, intent(out) :: error
      call gradient_case(2, 3, error)
   end subroutine gradient_molecular

   subroutine gradient_ec(error)
      type(error_type), allocatable, intent(out) :: error
      call gradient_case(1, 2, error)
   end subroutine gradient_ec

   subroutine gradient_drism_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      call gradient_case(1, 3, error, drism=.true.)
   end subroutine gradient_drism_cartesian

   subroutine gradient_drism_charges(error)
      type(error_type), allocatable, intent(out) :: error
      call gradient_case(1, 1, error, drism=.true.)
   end subroutine gradient_drism_charges

   subroutine gradient_drism_molecular(error)
      type(error_type), allocatable, intent(out) :: error
      call gradient_case(2, 3, error, drism=.true.)
   end subroutine gradient_drism_molecular

   subroutine gradient_drism_ec(error)
      type(error_type), allocatable, intent(out) :: error
      call gradient_case(1, 2, error, drism=.true.)
   end subroutine gradient_drism_ec

   !> Physical fixed-bulk water, with a weak solute on the coarse molecular grid;
   !> molecular point membership and k lattice stay
   !> fixed, with both reference owners tested. The EC case adds the host's
   !> fixed-point contraction of a geometry-dependent analytic potential.
   !> Molecular cases use a smaller grid, one deformation and one FD step
   subroutine gradient_case(kind, mode, error, drism)
      integer, intent(in) :: kind, mode
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in), optional :: drism
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(potential_sites_type) :: sites
      type(structure_type) :: mol, moved
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      class(response_item_type), allocatable :: item
      character(len=4), parameter :: corrections(3) = ["none", "PC  ", "PC+ "]
      real(wp), parameter :: steps(2) = [2.0e-4_wp, 1.0e-4_wp]
      real(wp), parameter :: translation(3) = [0.17_wp, -0.13_wp, 0.11_wp]
      real(wp), allocatable :: grad(:, :), delta(:, :), phi(:), spatial(:, :), host_delta(:), wphi(:)
      real(wp), allocatable :: q(:), dq(:), wq(:)
      real(wp) :: e0, e1, ep, em, dot, fd, eps, max_error
      integer :: ic, ref, nref, direction, istep, nstep, ndirection, fd_evaluations, nrad, degree
      logical :: compact

      call new_context(ctx, nthreads=2, verbosity=0)
      max_error = 0.0_wp
      fd_evaluations = 0
      nref = 1
      if (kind == 2) nref = 2
      compact = kind == 2 .and. mode == 3
      nstep = merge(1, size(steps), compact)
      ndirection = merge(1, 2, compact)
      nrad = merge(18, 36, compact)
      degree = merge(7, 11, compact)
      do ref = 1, nref
         do ic = 1, size(corrections)
            if (kind == 2) then
               write (output_unit, "(a,i0,a,a)") "  Molecular gradient reference=", ref, &
                  & " correction=", trim(corrections(ic))
               flush (output_unit)
            end if
            call host_model(ctx, trim(corrections(ic)), mode, model, solvent, err, drism=drism, &
               & reuse_bulk=kind == 2 .and. mode == 3)
            if (.not. allocated(err)) call solvent%potential%sites(sites, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            mol = solvent%structure()
            if (kind == 2) then
               call get_qc_handymod_recipe(recipe, err, nrad=nrad, degree=degree, rmax=6.0_wp)
               if (.not. allocated(err)) call new_molecular_point_grid(grid, err, recipe=recipe, dr=0.8_wp, &
                  & nufft_tol=1.0e-12_wp, weight_threshold=0.0_wp, pruning_threshold=0.0_wp)
               if (.not. allocated(err)) call grid%update(mol, err)
               if (ref == 2 .and. .not. allocated(err)) call molecular_grid_set_kgrid(grid, 0.8_wp, err, buffer=2.0_wp, &
                  & reference=sum(mol%xyz, dim=2)/real(mol%nat, wp))
               if (allocated(model%grid)) deallocate (model%grid)
               allocate (model%grid, source=grid)
               model%grid%nthreads = 2
               if (.not. allocated(err)) call model%update(mol, err)
               call check_moist_error(error, err)
               if (allocated(error)) return
            end if
            allocate (grad(3, mol%nat), delta(3, mol%nat))
            do direction = 1, ndirection
               delta = 0.0_wp
               if (direction == 1) then
                  delta(:, 1) = [0.13_wp, -0.21_wp, 0.17_wp]
                  delta(:, 2) = [-0.19_wp, 0.11_wp, -0.07_wp]
                  delta(:, 3) = [0.23_wp, 0.09_wp, -0.13_wp]
               else
                  delta = spread(translation, 2, mol%nat)
               end if
               if (direction == 1 .and. .not. compact) then
                  call model%update(mol, err)
               else
                  call geometry_energy(mol, e0)
                  if (allocated(error)) return
               end if
               if (.not. allocated(err)) call model%new_coupling(coupling, err)
               if (.not. allocated(err)) call model%prepare_gradient(coupling, err)
               if (.not. allocated(err) .and. mode == 2) then
                  call analytic_host(model%grid, mol, sites%q, delta, phi, spatial, host_delta)
                  if (coupling%next()) then
                     call coupling%answer("phi", phi, err)
                     if (.not. allocated(err)) call coupling%answer("dphi_dr", spatial, err)
                  end if
               else if (.not. allocated(err) .and. mode == 1) then
                  call analytic_charges(mol, sites%q, delta, q, dq)
                  if (coupling%next()) call coupling%answer("q", q, err)
               end if
               grad = 0.37_wp
               if (.not. allocated(err)) call model%get_gradient(coupling, response, grad, err)
               call model%release_coupling(coupling)
               call check_moist_error(error, err)
               if (allocated(error)) return
               grad = grad - 0.37_wp
               dot = sum(grad*delta)
               if (mode == 2) then
                  if (allocated(wphi)) deallocate (wphi)
                  do while (response%next())
                     item = response%item()
                     select type (item)
                     type is (potential_adjoint_response_type)
                        wphi = item%w_phi
                     end select
                  end do
                  call check(error, allocated(wphi), more="gradient omitted the EC host response")
                  if (allocated(error)) return
                  dot = dot + sum(wphi*host_delta)
               else if (mode == 1) then
                  if (allocated(wq)) deallocate (wq)
                  do while (response%next())
                     item = response%item()
                     select type (item)
                     type is (atomic_charge_adjoint_response_type)
                        wq = item%dg_dq
                     end select
                  end do
                  call check(error, allocated(wq), more="gradient omitted the host charge response")
                  if (allocated(error)) return
                  dot = dot + sum(wq*dq)
               else
                  call check(error, .not. response%next(), more="fixed terms returned a host response")
                  if (allocated(error)) return
               end if
               if (direction == 1 .and. .not. compact) call current_energy(mol, e0)
               call current_energy(mol, e1)
               if (allocated(error)) return
               call check(error, e0 == e1, more="gradient changed cached energy")
               if (allocated(error)) return
               if (compact) then
                  call check(error, sum(sum(grad, dim=2)*translation), 0.0_wp, thr=1.0e-9_wp, &
                     & more="whole-system translation changes energy")
                  if (allocated(error)) return
               end if
               do istep = 1, nstep
                  eps = steps(istep)
                  moved = mol
                  moved%xyz = mol%xyz + eps*delta
                  call geometry_energy(moved, ep)
                  moved%xyz = mol%xyz - eps*delta
                  call geometry_energy(moved, em)
                  if (allocated(error)) return
                  fd_evaluations = fd_evaluations + 2
                  fd = (ep - em)/(2.0_wp*eps)
                  max_error = max(max_error, abs(dot - fd)/max(1.0e-5_wp, abs(fd)))
                  call check(error, dot, fd, thr=5.0e-5_wp*max(1.0e-5_wp, abs(fd)), &
                     & more=trim(corrections(ic))//" model nuclear/host gradient")
                  if (allocated(error)) return
               end do
               if (direction == 2) then
                  call check(error, dot, 0.0_wp, thr=1.0e-9_wp, more="whole-system translation changes energy")
                  if (allocated(error)) return
               end if
            end do
            deallocate (grad, delta)
         end do
      end do
      write (*, "(a,i0,a,i0,a,es12.4)") "  3D model gradient FD grid=", kind, " host=", mode, " max scaled error=", max_error
      write (*, "(a,i0)") "  Displaced energy evaluations: ", fd_evaluations

   contains

      subroutine geometry_energy(active, energy)
         type(structure_type), intent(in) :: active
         real(wp), intent(out) :: energy
         call model%update(active, err)
         if (.not. allocated(err)) call current_energy(active, energy)
         call check_moist_error(error, err)
      end subroutine geometry_energy

      subroutine current_energy(active, energy)
         type(structure_type), intent(in) :: active
         real(wp), intent(out) :: energy
         if (mode == 2) call analytic_host(model%grid, active, sites%q, delta, phi, spatial, host_delta)
         if (mode == 1) call analytic_charges(active, sites%q, delta, q, dq)
         if (.not. allocated(err)) then
            if (mode == 2) then
               call host_energy(model, phi, .true., energy, err)
            else if (mode == 1) then
               call host_energy(model, q, .false., energy, err)
            else
               call host_energy(model, [0.0_wp], .false., energy, err)
            end if
         end if
         call check_moist_error(error, err)
      end subroutine current_energy

   end subroutine gradient_case

   !> Neutral charge transfer depending on the solute geometry, and its direction.
   subroutine analytic_charges(mol, charges, delta, q, dq)
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: charges(:), delta(:, :)
      real(wp), allocatable, intent(out) :: q(:), dq(:)
      q = charges
      q(1) = q(1) + 0.01_wp*(mol%xyz(1, 2) - mol%xyz(1, 3))
      q(2) = q(2) - 0.01_wp*(mol%xyz(1, 2) - mol%xyz(1, 3))
      allocate (dq(size(charges)), source=0.0_wp)
      dq(1) = 0.01_wp*(delta(1, 2) - delta(1, 3))
      dq(2) = -dq(1)
   end subroutine analytic_charges

   !> Smooth analytic host with moving charge centers and geometry-dependent
   !> charges. spatial is dphi/dr; host_delta holds dphi/dR at fixed r.
   subroutine analytic_host(grid, mol, charges, delta, phi, spatial, host_delta)
      class(moist_math_grid_3d_type), intent(in) :: grid
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: charges(:), delta(:, :)
      real(wp), allocatable, intent(out) :: phi(:), spatial(:, :), host_delta(:)
      real(wp), allocatable :: q(:), dq(:)
      real(wp) :: d(3), radius
      integer :: i, a
      call analytic_charges(mol, charges, delta, q, dq)
      allocate (phi(grid%ngrid), spatial(3, grid%ngrid), host_delta(grid%ngrid))
      phi = 0.0_wp
      spatial = 0.0_wp
      host_delta = 0.0_wp
      do i = 1, grid%ngrid
         do a = 1, mol%nat
            d = grid%point(i) - mol%xyz(:, a)
            radius = sqrt(sum(d**2) + 0.04_wp)
            phi(i) = phi(i) + q(a)/radius
            spatial(:, i) = spatial(:, i) - q(a)*d/radius**3
            host_delta(i) = host_delta(i) + q(a)*sum(d*delta(:, a))/radius**3 + dq(a)/radius
         end do
      end do
   end subroutine analytic_host

end module test_moz_3d_response
