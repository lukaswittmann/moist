!> Restarted flexible right-preconditioned GMRES inner solve
submodule(moist_model_moz_solver_jfnk) moist_model_moz_solver_jfnk_fgmres
   use moist_math_orthogonalization, only: orthogonalize_mgs => moist_orthogonalize_mgs
   use moist_model_moz_solver_jfnk_kernels, only: blocked_dot, copy_vector, divide_vector, project_basis, update_basis, &
      & generate_givens, apply_givens, upper_backsolve
   implicit none(type, external)

contains

   module subroutine fgmres_linear_solve(self, residual_cb, gamma_curr, f_curr, rhs_flat, delta_flat, &
                                         eta, feasible, error, gmres_iters)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      real(wp), intent(in) :: gamma_curr(:, :)
      real(wp), intent(in) :: f_curr(:, :)
      real(wp), intent(in) :: rhs_flat(:)
      real(wp), intent(out) :: delta_flat(:)
      real(wp), intent(in) :: eta
      logical, intent(out) :: feasible
      type(error_type), allocatable, intent(out) :: error
      integer, intent(out), optional :: gmres_iters

      !> Vector length, cycle length, cycle step, steps kept by the least-squares solve
      integer :: nflat, m, k, kk, i, j
      !> Arnoldi steps so far
      integer :: gmres_iter
      !> Cycle outcome
      logical :: converged, breakdown, cycle_end, usable
      !> Residual norms and tolerance
      real(wp) :: beta, rho, h_jj, tol, rho_prev_restart, znorm

      feasible = .true.
      nflat = size(rhs_flat)
      m = max(1, self%params%restart)
      gmres_iter = 0
      if (present(gmres_iters)) gmres_iters = 0
      self%nrecycled = 0
      call self%workspace%ensure_flexible(nflat, m)

      !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size))
      do i = 1, nflat
         delta_flat(i) = 0.0_wp
         self%workspace%ws_r(i) = rhs_flat(i)
      end do
      !$omp end parallel do
      beta = sqrt(blocked_dot(self%workspace%ws_r, self%workspace%ws_r, self%team_size()))
      if (.not. ieee_is_finite(beta)) then
         feasible = .false.
         call fatal_error(error, "Residual norm must be finite")
         return
      end if
      if (beta == 0.0_wp) return
      tol = max(tiny(1.0_wp), eta*beta)
      rho_prev_restart = beta

      call start_cycle()
      do while (gmres_iter < self%params%max_krylov)
         gmres_iter = gmres_iter + 1
         k = k + 1

         ! z_k = M^{-1} v_k, kept for the update; an unusable direction becomes v_k
         call copy_vector(self%workspace%ws_v(:, k), self%workspace%ws_zbasis(:, k), self%team_size())
         call self%apply_precond(self%workspace%ws_zbasis(:, k))
         znorm = sqrt(blocked_dot(self%workspace%ws_zbasis(:, k), self%workspace%ws_zbasis(:, k), self%team_size()))
         usable = ieee_is_finite(znorm) .and. znorm > 0.0_wp
         if (.not. usable) then
            self%nreplaced_stage = self%nreplaced_stage + 1
            call copy_vector(self%workspace%ws_v(:, k), self%workspace%ws_zbasis(:, k), self%team_size())
         end if

         ! w = J z_k
         call self%apply_jacobian(residual_cb, gamma_curr, f_curr, self%workspace%ws_zbasis(:, k), &
                                  self%workspace%ws_w, feasible, error)
         if (.not. feasible) then
            if (present(gmres_iters)) gmres_iters = gmres_iter
            return
         end if

         ! Arnoldi orthogonalization of w against v_1..v_k
         if (self%params%cgs2) then
            call project_basis(nflat, k, self%workspace%ws_v(:, 1:k), self%workspace%ws_w, &
               & self%workspace%ws_proj(1:k), self%team_size())
            call update_basis(nflat, k, self%workspace%ws_v(:, 1:k), self%workspace%ws_proj(1:k), -1.0_wp, &
               & self%workspace%ws_w, self%team_size())
            call project_basis(nflat, k, self%workspace%ws_v(:, 1:k), self%workspace%ws_w, &
               & self%workspace%ws_proj2(1:k), self%team_size())
            call update_basis(nflat, k, self%workspace%ws_v(:, 1:k), self%workspace%ws_proj2(1:k), -1.0_wp, &
               & self%workspace%ws_w, self%team_size())
            self%workspace%ws_h(1:k, k) = self%workspace%ws_proj(1:k) + self%workspace%ws_proj2(1:k)
         else
            call orthogonalize_mgs(self%workspace%ws_v(:, 1:k), self%workspace%ws_w, self%workspace%ws_h(1:k, k), &
               & nthreads=self%team_size())
         end if

         h_jj = sqrt(blocked_dot(self%workspace%ws_w, self%workspace%ws_w, self%team_size()))
         if (.not. ieee_is_finite(h_jj) .or. .not. all(ieee_is_finite(self%workspace%ws_h(1:k, k)))) then
            feasible = .false.
            if (present(gmres_iters)) gmres_iters = gmres_iter
            call fatal_error(error, "GMRES Jacobian product and preconditioner must be finite")
            return
         end if
         self%workspace%ws_h(k + 1, k) = h_jj
         breakdown = h_jj <= 10.0_wp*epsilon(1.0_wp)*maxval(abs(self%workspace%ws_h(1:k, k)))
         if (.not. breakdown) then
            call divide_vector(self%workspace%ws_w, h_jj, self%workspace%ws_v(:, k + 1), self%team_size())
         end if

         do j = 1, k - 1
            call apply_givens(self%workspace%ws_cs(j), self%workspace%ws_sn(j), self%workspace%ws_h(j, k), &
               & self%workspace%ws_h(j + 1, k))
         end do
         call generate_givens(self%workspace%ws_h(k, k), self%workspace%ws_h(k + 1, k), &
            & self%workspace%ws_cs(k), self%workspace%ws_sn(k))
         call apply_givens(self%workspace%ws_cs(k), self%workspace%ws_sn(k), self%workspace%ws_h(k, k), &
            & self%workspace%ws_h(k + 1, k))
         call apply_givens(self%workspace%ws_cs(k), self%workspace%ws_sn(k), self%workspace%ws_g(k), &
            & self%workspace%ws_g(k + 1))

         rho = abs(self%workspace%ws_g(k + 1))
         cycle_end = rho < tol .or. breakdown .or. k == m .or. gmres_iter == self%params%max_krylov

         if (associated(self%iter_printer)) then
            if (self%iter_printer%inner_rows .and. &
               & (gmres_iter <= 4 .or. mod(gmres_iter, 5) == 0 .or. cycle_end)) then
               call self%iter_printer%print_inner_gmres_row(gmres_iter, k, rho, tol)
            end if
         end if
         if (.not. cycle_end) cycle

         ! Least-squares solution; a projected matrix singular up to rounding drops its last
         ! columns, on the scale of the breakdown test
         kk = k
         do while (kk > 0)
            call upper_backsolve(self%workspace%ws_h(1:kk, 1:kk), self%workspace%ws_g(1:kk), &
               & self%workspace%ws_y(1:kk), usable, rank_tol=10.0_wp*epsilon(1.0_wp))
            if (usable) exit
            kk = kk - 1
         end do
         if (kk > 0) then
            call update_basis(nflat, kk, self%workspace%ws_zbasis(:, 1:kk), self%workspace%ws_y(1:kk), 1.0_wp, &
               & delta_flat, self%team_size())
         end if
         converged = rho < tol .and. kk == k
         if (converged) exit
         if (gmres_iter == self%params%max_krylov) then
            self%ncapped_stage = self%ncapped_stage + 1
            exit
         end if
         ! No usable column: another cycle from the same residual would repeat it
         if (kk == 0) exit
         if (rho_prev_restart > 0.0_wp .and. self%params%stagnation_tol > 0.0_wp .and. kk == k) then
            if ((rho_prev_restart - rho)/rho_prev_restart < self%params%stagnation_tol) then
               self%nstagnated_stage = self%nstagnated_stage + 1
               exit
            end if
         end if

         ! Restart from the unpreconditioned residual b - J delta
         call self%apply_jacobian(residual_cb, gamma_curr, f_curr, delta_flat, self%workspace%ws_w, feasible, error)
         if (.not. feasible) then
            if (present(gmres_iters)) gmres_iters = gmres_iter
            return
         end if
         !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size))
         do i = 1, nflat
            self%workspace%ws_r(i) = rhs_flat(i) - self%workspace%ws_w(i)
         end do
         !$omp end parallel do
         beta = sqrt(blocked_dot(self%workspace%ws_r, self%workspace%ws_r, self%team_size()))
         if (.not. ieee_is_finite(beta)) then
            feasible = .false.
            if (present(gmres_iters)) gmres_iters = gmres_iter
            call fatal_error(error, "Restarted residual norm must be finite")
            return
         end if
         ! The next cycle is measured from this residual: after dropped columns the estimate
         ! rho belongs to a step that was not taken
         rho_prev_restart = beta
         if (beta < tol) exit
         call start_cycle()
      end do

      if (present(gmres_iters)) gmres_iters = gmres_iter

   contains

      !> First basis vector from ws_r and its norm beta; zero Hessenberg matrix and
      !> rotated right-hand side beta e1
      subroutine start_cycle()
         call divide_vector(self%workspace%ws_r, beta, self%workspace%ws_v(:, 1), self%team_size())
         self%workspace%ws_g = 0.0_wp
         self%workspace%ws_g(1) = beta
         self%workspace%ws_h = 0.0_wp
         k = 0
      end subroutine start_cycle

   end subroutine fgmres_linear_solve

end submodule moist_model_moz_solver_jfnk_fgmres
