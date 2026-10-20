!> Restarted GMRES inner solve with optional DEF1 coarse-space deflation
submodule(moist_model_moz_solver_jfnk) moist_model_moz_solver_jfnk_gmres
   use moist_math_orthogonalization, only: orthogonalize_mgs => moist_orthogonalize_mgs
   use moist_model_moz_solver_jfnk_kernels, only: blocked_dot, copy_vector, divide_vector, project_basis, update_basis
   implicit none(type, external)

contains

   module subroutine gmres_linear_solve(self, residual_cb, gamma_curr, f_curr, rhs_flat, delta_flat, &
                                        eta, precond_valid, feasible, error, gmres_iters)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      real(wp), intent(in) :: gamma_curr(:, :)
      real(wp), intent(in) :: f_curr(:, :)
      real(wp), intent(in) :: rhs_flat(:)
      real(wp), intent(out) :: delta_flat(:)
      real(wp), intent(in) :: eta
      logical, intent(inout) :: precond_valid
      logical, intent(out) :: feasible
      type(error_type), allocatable, intent(out) :: error
      integer, intent(out), optional :: gmres_iters

      integer :: nflat, m, k, j, i, gmres_iter, nrestart, last_iter
      !> Recycled vectors in use, cycle length (Arnoldi steps), maximum recycle dimension,
      !> Jacobian products of the recycle refresh
      integer :: kr, mk, kmax, njv
      !> Cycle outcome
      logical :: converged, last_cycle, happy_breakdown
      !> Whether this Arnoldi step ends the cycle
      logical :: cycle_end
      real(wp) :: beta, rho, h_jj, tol, rho_prev_restart
      !> B y, the recycled-space part of the cycle correction
      real(wp), allocatable :: by(:)
      !> Whether this solve is deflated (DEF1)
      logical :: defl
      !> |P b| of the undeflated right-hand side
      real(wp) :: beta_rhs
      !> Coarse coordinates E^{-1} Z^T Omega w, (kz, 1)
      real(wp), allocatable :: ycoarse(:, :)

      feasible = .true.
      nflat = size(rhs_flat)
      !> Krylov dimension of one cycle; matches `jfnk_ensure_workspace`
      m = max(1, self%params%restart)
      kmax = self%recycle_dim()
      last_iter = 0
      gmres_iter = 0
      kr = 0
      if (present(gmres_iters)) gmres_iters = 0

      ! A deflated solve skips recycling and drops the recycled subspace
      defl = self%deflation%active
      if (defl) then
         if (kmax > 0) self%nunrecycled_stage = self%nunrecycled_stage + 1
         kmax = 0
         self%nrecycled = 0
         allocate (ycoarse(self%deflation%kz, 1))
         if (allocated(self%workspace%ws_z)) then
            if (size(self%workspace%ws_z) /= nflat) deallocate (self%workspace%ws_z)
         end if
         if (.not. allocated(self%workspace%ws_z)) allocate (self%workspace%ws_z(nflat))
      end if

      !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size))
      do i = 1, nflat
         delta_flat(i) = 0.0_wp
         self%workspace%ws_r(i) = rhs_flat(i)
      end do
      !$omp end parallel do
      if (defl) then
         ! Newton residual of the final direction equals the deflated GMRES
         ! residual, P (b - J delta) = P (I - J Q)(b - J x),
         ! so the inexact-Newton tolerance stays eta*|P b| of the
         ! undeflated right-hand side
         ! The oblique projector can inflate the right-hand side, so keep
         ! the tolerance relative to the undeflated preconditioned residual
         call copy_vector(rhs_flat, self%workspace%ws_w, self%team_size())
         if (precond_valid) call self%apply_precond(self%workspace%ws_w)
         beta_rhs = sqrt(blocked_dot(self%workspace%ws_w, self%workspace%ws_w, self%team_size()))
         if (.not. ieee_is_finite(beta_rhs)) then
            feasible = .false.
            call fatal_error(error, "Preconditioned residual norm must be finite")
            return
         end if
         if (beta_rhs == 0.0_wp .and. any(rhs_flat /= 0.0_wp)) then
            feasible = .false.
            call fatal_error(error, "Preconditioner annihilated a nonzero residual")
            return
         end if
         ! r := (I - J Q) b
         call self%deflate(residual_cb, gamma_curr, f_curr, self%workspace%ws_r, ycoarse, feasible, error)
         if (.not. feasible) return
      end if
      if (precond_valid) call self%apply_precond(self%workspace%ws_r)
      beta =sqrt(blocked_dot(self%workspace%ws_r, self%workspace%ws_r, self%team_size()))
      if (.not. ieee_is_finite(beta)) then
         feasible = .false.
         call fatal_error(error, "Preconditioned residual norm must be finite")
         return
      end if
      if (beta == 0.0_wp) then
         ! An exact coarse-space correction can legitimately remove the RHS
         ! An invertible preconditioner alone cannot remove a nonzero RHS
         if (precond_valid .and. .not. defl .and. any(rhs_flat /= 0.0_wp)) then
            feasible = .false.
            call fatal_error(error, "Preconditioner annihilated a nonzero residual")
            return
         end if
         if (defl) call self%deflation_finish(residual_cb, gamma_curr, f_curr, rhs_flat, &
                                              delta_flat, .true., ycoarse, feasible, error)
         return
      end if

      tol = eta*beta
      if (defl) tol = eta*beta_rhs
      ! Preserve relative accuracy when the left operator is rescaled
      tol = max(tiny(1.0_wp), tol)

      ! Recycled subspace: rebuild C = A U for this Jacobian, then take
      ! the minimal-residual correction in span U,
      ! delta = U C^T r, r := r - C C^T r
      ! A subspace built for the other preconditioning state is dropped
      if (kmax > 0 .and. self%nrecycled > 0) then
         if (precond_valid .eqv. self%recycle_precond) then
            call self%recycle_refresh(residual_cb, gamma_curr, f_curr, precond_valid, njv, &
                                      feasible, error)
            gmres_iter = njv
            last_iter = njv
            if (.not. feasible) then
               if (present(gmres_iters)) gmres_iters = last_iter
               return
            end if
            kr = self%nrecycled
            if (kr > 0) then
               call project_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_r, &
                  & self%workspace%ws_proj(1:kr), self%team_size())
               call update_basis(nflat, kr, self%workspace%ws_u(:, 1:kr), self%workspace%ws_proj(1:kr), 1.0_wp, &
                                 delta_flat, self%team_size())
               call update_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_proj(1:kr), -1.0_wp, &
                                 self%workspace%ws_r, self%team_size())
               beta = sqrt(blocked_dot(self%workspace%ws_r, self%workspace%ws_r, self%team_size()))
            end if
         else
            self%nrecycled = 0
         end if
      end if
      if (kmax > 0) allocate (by(kmax))

      if (.not. ieee_is_finite(beta)) then
         feasible = .false.
         call fatal_error(error, "Recycled residual norm must be finite")
         return
      end if
      if (beta < tol) then
         if (present(gmres_iters)) gmres_iters = last_iter
         if (defl) call self%deflation_finish(residual_cb, gamma_curr, f_curr, rhs_flat, &
                                              delta_flat, .true., ycoarse, feasible, error)
         return
      end if

      call divide_vector(self%workspace%ws_r, beta, self%workspace%ws_v(:, 1), self%team_size())
      self%workspace%ws_g = 0.0_wp
      self%workspace%ws_g(1) = beta
      self%workspace%ws_h = 0.0_wp
      if (kmax > 0) then
         self%workspace%ws_hbar = 0.0_wp
         self%workspace%ws_bproj = 0.0_wp
      end if
      nrestart = 0
      rho_prev_restart = beta
      !> Arnoldi steps per cycle: the recycled vectors take part of the m budget
      mk = max(1, m - kr)
      k = 0

      do while (gmres_iter < self%params%max_krylov)
         gmres_iter = gmres_iter + 1
         last_iter = gmres_iter
         k = k + 1

         ! Arnoldi step: compute w = M^{-1} J v_k using preallocated workspace
         call self%apply_jacobian(residual_cb, gamma_curr, f_curr, self%workspace%ws_v(:, k), self%workspace%ws_w, &
                                  feasible, error)
         if (.not. feasible) then
            if (present(gmres_iters)) gmres_iters = last_iter
            return
         end if

         ! DEF1: w := (I - J Q) J v,
         ! computed as J (v - Q J v), one more product
         if (defl) then
            call self%deflation%coordinates(self%workspace%ws_w, ycoarse, nthreads=self%team_size())
            call copy_vector(self%workspace%ws_v(:, k), self%workspace%ws_z, self%team_size())
            call self%deflation%expand(ycoarse, -1.0_wp, self%workspace%ws_z, nthreads=self%team_size())
            call self%apply_jacobian(residual_cb, gamma_curr, f_curr, self%workspace%ws_z, self%workspace%ws_w, &
                                     feasible, error)
            if (.not. feasible) then
               if (present(gmres_iters)) gmres_iters = last_iter
               return
            end if
         end if

         ! Apply left preconditioner
         if (precond_valid) call self%apply_precond(self%workspace%ws_w)

         ! GCRO deflation: w := (I - C C^T) w,
         ! keeping B(:, k) = C^T A v_k
         if (kr > 0) then
            call project_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_w, &
               & self%workspace%ws_bproj(1:kr, k), self%team_size())
            call update_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_bproj(1:kr, k), -1.0_wp, &
                              self%workspace%ws_w, self%team_size())
         end if

         ! Arnoldi orthogonalization of w against v_1..v_k
         if (self%params%cgs2) then
            ! Classical Gram-Schmidt with one reorthogonalization pass (CGS2)
            !
            ! - each pass streams the large w vector only twice: a BLAS-2
            !   projection s = V^T w and a back-update w := w - V s, instead
            !   of the 2*k separate sweeps modified Gram-Schmidt makes
            ! - the second pass reduces loss of orthogonality
            ! - fewer sweeps do not guarantee a faster solve
            ! - rounding can also change the Newton path
            ! - CGS2 is opt-in and MGS is the default
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
            if (present(gmres_iters)) gmres_iters = last_iter
            call fatal_error(error, "GMRES Jacobian product and preconditioner must be finite")
            return
         end if
         self%workspace%ws_h(k + 1, k) = h_jj

         ! Compare against this operator column, rather than an absolute
         ! cutoff that changes meaning when the preconditioner is rescaled
         happy_breakdown = h_jj <= 10.0_wp*epsilon(1.0_wp)*maxval(abs(self%workspace%ws_h(1:k, k)))
         if (.not. happy_breakdown) then
            call divide_vector(self%workspace%ws_w, h_jj, self%workspace%ws_v(:, k + 1), self%team_size())
         else
            ! Happy breakdown: the Krylov space is invariant
            ! Continuing would otherwise reuse whatever stale content
            ! ws_v(:,k+1) holds as an orthonormal basis vector
            !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size))
            do i = 1, nflat
               self%workspace%ws_v(i, k + 1) = 0.0_wp
            end do
            !$omp end parallel do
         end if
         ! Unrotated Hessenberg column for the harmonic Ritz harvest
         if (kmax > 0) self%workspace%ws_hbar(1:k + 1, k) = self%workspace%ws_h(1:k + 1, k)

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
         ! Reaching the iteration cap also ends the cycle, so a capped solve
         ! (max_krylov < restart) still returns its Newton direction
         cycle_end = rho < tol .or. happy_breakdown .or. k == mk .or. gmres_iter == self%params%max_krylov

         ! Inner rows: the first four, every fifth and the last of each cycle
         if (associated(self%iter_printer)) then
            if (self%iter_printer%inner_rows .and. &
               & (gmres_iter <= 4 .or. mod(gmres_iter, 5) == 0 .or. cycle_end)) then
               call self%iter_printer%print_inner_gmres_row(gmres_iter, k, rho, tol)
            end if
         end if

         ! Least-squares solution is only formed here
         if (cycle_end) then
            call upper_backsolve(self%workspace%ws_h(1:k, 1:k), self%workspace%ws_g(1:k), &
               & self%workspace%ws_y(1:k), feasible)
            if (.not. feasible) then
               if (present(gmres_iters)) gmres_iters = last_iter
               call fatal_error(error, "GMRES encountered a singular or non-finite Krylov system")
               return
            end if
            call update_basis(nflat, k, self%workspace%ws_v(:, 1:k), self%workspace%ws_y(1:k), 1.0_wp, &
               & delta_flat, self%team_size())
            ! GCRO: the recycled coordinates of the correction are -B y
            if (kr > 0) then
               by(1:kr) = matmul(self%workspace%ws_bproj(1:kr, 1:k), self%workspace%ws_y(1:k))
               call update_basis(nflat, kr, self%workspace%ws_u(:, 1:kr), by(1:kr), -1.0_wp, delta_flat, &
                  & self%team_size())
            end if
            converged = rho < tol
            last_cycle = converged .or. happy_breakdown .or. gmres_iter == self%params%max_krylov
            ! Iteration cap reached short of the tolerance, for the stage summary
            if (.not. converged .and. .not. happy_breakdown .and. gmres_iter == self%params%max_krylov) then
               self%ncapped_stage = self%ncapped_stage + 1
            end if

            ! Harvest the recycled subspace from this cycle: always once one
            ! exists, otherwise from the first cycle that filled the restart
            ! length or the iteration budget, so a solve capped below the
            ! restart length (max_krylov < restart) still starts recycling
            ! C is only needed for another cycle of this solve; the next
            ! Newton step rebuilds it for its own Jacobian
            if (kmax > 0 .and. (kr > 0 .or. k == mk .or. gmres_iter == self%params%max_krylov) &
               & .and. kr + k > kmax) then
               call self%recycle_harvest(kr, k, .not. last_cycle)
               self%recycle_precond = precond_valid
            end if

            if (last_cycle) exit

            ! Check for restart stagnation: bail out if this cycle barely reduced residual
            if (rho_prev_restart > 0.0_wp .and. self%params%stagnation_tol > 0.0_wp) then
               if ((rho_prev_restart - rho)/rho_prev_restart < self%params%stagnation_tol) then
                  self%nstagnated_stage = self%nstagnated_stage + 1
                  exit
               end if
            end if
            rho_prev_restart = rho

            ! Restart: recompute preconditioned residual r = M^{-1}(b - J*delta)
            nrestart = nrestart + 1
            call copy_vector(rhs_flat, self%workspace%ws_r, self%team_size())
            call self%apply_jacobian(residual_cb, gamma_curr, f_curr, delta_flat, self%workspace%ws_w, &
                                     feasible, error)
            if (.not. feasible) then
               if (present(gmres_iters)) gmres_iters = last_iter
               return
            end if
            !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size))
            do i = 1, nflat
               self%workspace%ws_r(i) = self%workspace%ws_r(i) - self%workspace%ws_w(i)
            end do
            !$omp end parallel do
            if (defl) then
               call self%deflate(residual_cb, gamma_curr, f_curr, self%workspace%ws_r, ycoarse, feasible, &
                                 error)
               if (.not. feasible) then
                  if (present(gmres_iters)) gmres_iters = last_iter
                  return
               end if
            end if
            ! Apply left preconditioner to restarted residual
            if (precond_valid) call self%apply_precond(self%workspace%ws_r)
            ! Deflate the restarted residual against the harvested C
            kr = self%nrecycled
            if (kr > 0) then
               call project_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_r, &
                  & self%workspace%ws_proj(1:kr), self%team_size())
               call update_basis(nflat, kr, self%workspace%ws_u(:, 1:kr), self%workspace%ws_proj(1:kr), 1.0_wp, &
                                 delta_flat, self%team_size())
               call update_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_proj(1:kr), -1.0_wp, &
                                 self%workspace%ws_r, self%team_size())
            end if
            beta = sqrt(blocked_dot(self%workspace%ws_r, self%workspace%ws_r, self%team_size()))
            if (.not. ieee_is_finite(beta)) then
               feasible = .false.
               if (present(gmres_iters)) gmres_iters = last_iter
               call fatal_error(error, "Restarted residual norm must be finite")
               return
            end if
            if (beta < tol) exit

            call divide_vector(self%workspace%ws_r, beta, self%workspace%ws_v(:, 1), self%team_size())
            self%workspace%ws_g = 0.0_wp
            self%workspace%ws_g(1) = beta
            self%workspace%ws_h = 0.0_wp
            if (kmax > 0) then
               self%workspace%ws_hbar = 0.0_wp
               self%workspace%ws_bproj = 0.0_wp
            end if
            mk = max(1, m - kr)
            k = 0
         end if
      end do

      if (present(gmres_iters)) gmres_iters = last_iter
      if (defl) call self%deflation_finish(residual_cb, gamma_curr, f_curr, rhs_flat, delta_flat, &
                                           .false., ycoarse, feasible, error)

   end subroutine gmres_linear_solve

   module subroutine gmres_deflate(self, residual_cb, gamma_curr, f_curr, vec, y, feasible, error)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      real(wp), intent(in) :: gamma_curr(:, :)
      real(wp), intent(in) :: f_curr(:, :)
      real(wp), intent(inout) :: vec(:)
      real(wp), intent(out) :: y(:, :)
      logical, intent(out) :: feasible
      type(error_type), allocatable, intent(out) :: error

      !> Vector index
      integer :: i

      call self%deflation%coordinates(vec, y, nthreads=self%team_size())
      !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, size(vec) >= parallel_vector_size))
      do i = 1, size(vec)
         self%workspace%ws_z(i) = 0.0_wp
      end do
      !$omp end parallel do
      call self%deflation%expand(y, 1.0_wp, self%workspace%ws_z, nthreads=self%team_size())
      call self%apply_jacobian(residual_cb, gamma_curr, f_curr, self%workspace%ws_z, self%workspace%ws_w, feasible, &
                               error)
      if (.not. feasible) return
      !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, size(vec) >= parallel_vector_size))
      do i = 1, size(vec)
         vec(i) = vec(i) - self%workspace%ws_w(i)
      end do
      !$omp end parallel do
   end subroutine gmres_deflate

   module subroutine gmres_deflation_finish(self, residual_cb, gamma_curr, f_curr, rhs_flat, delta_flat, &
                                            x_zero, y, feasible, error)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      real(wp), intent(in) :: gamma_curr(:, :)
      real(wp), intent(in) :: f_curr(:, :)
      real(wp), intent(in) :: rhs_flat(:)
      real(wp), intent(inout) :: delta_flat(:)
      logical, intent(in) :: x_zero
      real(wp), intent(out) :: y(:, :)
      logical, intent(out) :: feasible
      type(error_type), allocatable, intent(out) :: error

      !> Vector index
      integer :: i

      feasible = .true.
      if (x_zero) then
         call copy_vector(rhs_flat, self%workspace%ws_w, self%team_size())
      else
         call self%apply_jacobian(residual_cb, gamma_curr, f_curr, delta_flat, self%workspace%ws_w, &
                                  feasible, error)
         if (.not. feasible) return
         !$omp parallel do schedule(static) &
         !$omp & num_threads(merge(self%team_size(), 1, size(rhs_flat) >= parallel_vector_size))
         do i = 1, size(rhs_flat)
            self%workspace%ws_w(i) = rhs_flat(i) - self%workspace%ws_w(i)
         end do
         !$omp end parallel do
      end if
      call self%deflation%coordinates(self%workspace%ws_w, y, nthreads=self%team_size())
      call self%deflation%expand(y, 1.0_wp, delta_flat, nthreads=self%team_size())
   end subroutine gmres_deflation_finish

   !> Form a Givens rotation that removes b without an absolute scale cutoff
   !>
   !> @param[in] a First entry
   !> @param[in] b Entry to eliminate
   !> @param[out] c Rotation cosine
   !> @param[out] s Rotation sine
   subroutine generate_givens(a, b, c, s)
      !> First entry
      real(wp), intent(in) :: a
      !> Entry to eliminate
      real(wp), intent(in) :: b
      !> Rotation cosine
      real(wp), intent(out) :: c
      !> Rotation sine
      real(wp), intent(out) :: s
      real(wp) :: r

      if (b == 0.0_wp) then
         c = 1.0_wp
         s = 0.0_wp
      else if (abs(b) > abs(a)) then
         r = -a/b
         s = 1.0_wp/sqrt(1.0_wp + r**2)
         c = s*r
      else
         r = -b/a
         c = 1.0_wp/sqrt(1.0_wp + r**2)
         s = c*r
      end if
   end subroutine generate_givens

   !> Apply a Givens rotation to the pair (a, b)
   !>
   !> - a := c*a - s*b
   !> - b := s*a + c*b, with the old a
   !>
   !> @param[in]    c  rotation cosine
   !> @param[in]    s  rotation sine
   !> @param[inout] a  first component, rotated in place
   !> @param[inout] b  second component, rotated in place
   subroutine apply_givens(c, s, a, b)
      !> Rotation cosine
      real(wp), intent(in) :: c
      !> Rotation sine
      real(wp), intent(in) :: s
      !> First component, rotated in place
      real(wp), intent(inout) :: a
      !> Second component, rotated in place
      real(wp), intent(inout) :: b
      real(wp) :: tmp

      tmp = c*a - s*b
      b = s*a + c*b
      a = tmp
   end subroutine apply_givens

   !> Solve the small triangular Krylov system without dividing by a zero pivot
   !>
   !> @param[in] r Upper triangular matrix
   !> @param[in] rhs Rotated least-squares right-hand side
   !> @param[out] x Solution; only usable when feasible
   !> @param[out] feasible Nonsingular finite system and finite solution
   subroutine upper_backsolve(r, rhs, x, feasible)
      !> Upper triangular matrix
      real(wp), intent(in) :: r(:, :)
      !> Rotated right-hand side
      real(wp), intent(in) :: rhs(:)
      !> Solution coefficients
      real(wp), intent(out) :: x(:)
      !> Whether the solution is usable
      logical, intent(out) :: feasible

      integer :: n, i, j
      real(wp) :: sum_val

      n = size(rhs)
      feasible = .false.
      if (.not. all(ieee_is_finite(r)) .or. .not. all(ieee_is_finite(rhs))) return
      do i = n, 1, -1
         if (r(i, i) == 0.0_wp) return
         sum_val = rhs(i)
         do j = i + 1, n
            sum_val = sum_val - r(i, j)*x(j)
         end do
         x(i) = sum_val/r(i, i)
         if (.not. ieee_is_finite(x(i))) return
      end do
      feasible = .true.
   end subroutine upper_backsolve

end submodule moist_model_moz_solver_jfnk_gmres
