!> GCRO-DR recycled Krylov subspace: refresh of its images and harvest of harmonic Ritz vectors
submodule(moist_model_moz_solver_jfnk) moist_model_moz_solver_jfnk_recycle
   use moist_math_blas, only: gemm_rows_tn, gemm_rows_nn
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_lapack_potrf, only: wrap_potrf
   use moist_math_lapack_ggev, only: dggev
   use moist_model_moz_solver_jfnk_kernels, only: blocked_dot, reduction_block_size
   implicit none(type, external)

contains

   pure module function gmres_recycle_dim(self) result(kmax)
      class(moist_model_moz_solver_jfnk_type), intent(in) :: self
      integer :: kmax
      kmax = max(0, min(self%params%recycle, self%params%restart - 1, self%params%max_krylov - 1))
   end function gmres_recycle_dim

   module subroutine gmres_recycle_refresh(self, residual_cb, gamma_curr, f_curr, precond_valid, njv, &
                                           feasible, error)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      real(wp), intent(in) :: gamma_curr(:, :)
      real(wp), intent(in) :: f_curr(:, :)
      logical, intent(in) :: precond_valid
      integer, intent(out) :: njv
      logical, intent(out) :: feasible
      type(error_type), allocatable, intent(out) :: error

      !> Recycled count, column and pass indices, Cholesky status
      integer :: kr, i, pass, info
      !> Gram matrix / Cholesky factor and its inverse
      real(wp), allocatable :: g(:, :), rinv(:, :)
      !> Smallest and largest diagonal entry of the Cholesky factor
      real(wp) :: dmin, dmax

      feasible = .true.
      njv = 0
      kr = self%nrecycled
      if (kr == 0) return
      do i = 1, kr
         call self%apply_jacobian(residual_cb, gamma_curr, f_curr, self%workspace%ws_u(:, i), &
                                  self%workspace%ws_c(:, i), feasible, error)
         njv = njv + 1
         if (.not. feasible) then
            self%nrecycled = 0
            return
         end if
         if (precond_valid) then
            call self%apply_precond(self%workspace%ws_c(:, i))
         end if
      end do

      allocate (g(kr, kr), rinv(kr, kr))
      do pass = 1, 2
         call gram_basis(self%workspace%ws_c, kr, self%workspace%ws_c, kr, g, self%team_size())
         call wrap_potrf(g, info, "U")
         if (info == 0) then
            do i = 1, kr
               g(i + 1:kr, i) = 0.0_wp
            end do
            dmin = minval([(g(i, i), i=1, kr)])
            dmax = maxval([(g(i, i), i=1, kr)])
            if (dmin < sqrt(epsilon(1.0_wp))*dmax) info = 1
         end if
         if (info /= 0) then
            self%nrecycled = 0
            return
         end if
         call upper_inverse(g, rinv)
         call combine_basis(self%workspace%ws_c, kr, rinv, kr, self%team_size())
         call combine_basis(self%workspace%ws_u, kr, rinv, kr, self%team_size())
      end do
   end subroutine gmres_recycle_refresh

   module subroutine gmres_recycle_harvest(self, kr, j, need_c)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      integer, intent(in) :: kr
      integer, intent(in) :: j
      logical, intent(in) :: need_c

      !> Search-space and image dimensions, targets, counters, workspace size
      integer :: p, q, kmax, nsel, nkeep, i, ii, js, lwork
      !> LAPACK status
      integer(lapack_ik) :: lapack_info
      !> Eigenvalue ordering, kept QR columns
      integer, allocatable :: order(:), keep(:)
      logical, allocatable :: taken(:)
      !> Column scaling D, harmonic Ritz magnitudes, generalized eigenvalues
      real(wp), allocatable :: dsc(:), theta(:), alphar(:), alphai(:), betav(:), work(:)
      !> Small matrices of the harvest
      real(wp), allocatable :: gbar(:, :), wtv(:, :), aeig(:, :), beig(:, :), vr(:, :), vl(:, :)
      real(wp), allocatable :: psel(:, :), gp(:, :), qmat(:, :), rmat(:, :), rinv(:, :), pr(:, :)
      real(wp), allocatable :: t1(:, :), t2(:, :), pu(:, :)
      real(wp) :: tmp

      kmax = self%recycle_dim()
      p = kr + j
      q = p + 1
      if (p <= kmax .or. kmax == 0) return

      allocate (dsc(kr), gbar(q, p), wtv(q, p), source=0.0_wp)
      do i = 1, kr
         dsc(i) = 1.0_wp/sqrt(blocked_dot(self%workspace%ws_u(:, i), self%workspace%ws_u(:, i), self%team_size()))
         gbar(i, i) = dsc(i)
      end do
      if (kr > 0) gbar(1:kr, kr + 1:p) = self%workspace%ws_bproj(1:kr, 1:j)
      gbar(kr + 1:q, kr + 1:p) = self%workspace%ws_hbar(1:j + 1, 1:j)
      if (kr > 0) then
         allocate (t1(kr, kr), t2(j + 1, kr))
         call gram_basis(self%workspace%ws_c, kr, self%workspace%ws_u, kr, t1, self%team_size())
         call gram_basis(self%workspace%ws_v, j + 1, self%workspace%ws_u, kr, t2, self%team_size())
         do i = 1, kr
            wtv(1:kr, i) = t1(:, i)*dsc(i)
            wtv(kr + 1:q, i) = t2(:, i)*dsc(i)
         end do
      end if
      do i = 1, j
         wtv(kr + i, kr + i) = 1.0_wp
      end do

      ! Harmonic Ritz pairs: G^T G z = theta G^T (W^T Vh) z
      aeig = matmul(transpose(gbar), gbar)
      beig = matmul(transpose(gbar), wtv)
      lwork = 16*p + 16
      allocate (alphar(p), alphai(p), betav(p), vr(p, p), vl(1, 1), work(lwork), source=0.0_wp)
      call dggev("N", "V", int(p, lapack_ik), aeig, int(p, lapack_ik), beig, int(p, lapack_ik), &
                 alphar, alphai, betav, vl, 1_lapack_ik, vr, int(p, lapack_ik), work, &
                 int(lwork, lapack_ik), lapack_info)
      if (lapack_info /= 0) return

      ! Order by |theta|; infinite eigenvalues (beta = 0) last
      allocate (theta(p), order(p), taken(p))
      do i = 1, p
         tmp = hypot(alphar(i), alphai(i))
         if (abs(betav(i)) > tiny(1.0_wp)*max(1.0_wp, tmp)) then
            theta(i) = tmp/abs(betav(i))
         else
            theta(i) = huge(1.0_wp)
         end if
         order(i) = i
      end do
      do i = 2, p
         ii = order(i)
         js = i - 1
         do while (js >= 1)
            if (theta(order(js)) <= theta(ii)) exit
            order(js + 1) = order(js)
            js = js - 1
         end do
         order(js + 1) = ii
      end do

      ! Select up to k vectors; a conjugate pair enters as its real and imaginary parts
      allocate (psel(p, kmax), source=0.0_wp)
      taken = .false.
      nsel = 0
      do ii = 1, p
         i = order(ii)
         if (taken(i)) cycle
         if (theta(i) >= huge(1.0_wp)) exit
         if (alphai(i) == 0.0_wp) then
            if (nsel + 1 > kmax) exit
            nsel = nsel + 1
            psel(:, nsel) = vr(:, i)
            taken(i) = .true.
         else
            js = i
            if (alphai(i) < 0.0_wp) js = i - 1
            if (nsel + 2 > kmax) exit
            psel(:, nsel + 1) = vr(:, js)
            psel(:, nsel + 2) = vr(:, js + 1)
            nsel = nsel + 2
            taken(js) = .true.
            taken(js + 1) = .true.
         end if
      end do
      if (nsel == 0) return

      ! G P = Q R, dropping directions that G P does not resolve
      gp = matmul(gbar, psel(:, 1:nsel))
      allocate (qmat(q, nsel), rmat(nsel, nsel), keep(nsel))
      call small_qr(gp, qmat, rmat, keep, nkeep)
      if (nkeep == 0) return
      allocate (rinv(nkeep, nkeep))
      call upper_inverse(rmat(1:nkeep, 1:nkeep), rinv)
      pr = matmul(psel(:, keep(1:nkeep)), rinv)

      ! U := U~ P_u R^{-1} + V_j P_v R^{-1},
      ! C := C Q_c + V_{j+1} Q_v (in place)
      allocate (pu(kr, nkeep))
      do i = 1, kr
         pu(i, :) = dsc(i)*pr(i, :)
      end do
      call combine_basis(self%workspace%ws_u, kr, pu, nkeep, self%team_size(), self%workspace%ws_v, j, pr(kr + 1:p, :))
      if (need_c) then
         call combine_basis(self%workspace%ws_c, kr, qmat(1:kr, 1:nkeep), nkeep, self%team_size(), &
            self%workspace%ws_v, j + 1, qmat(kr + 1:q, 1:nkeep))
      end if
      self%nrecycled = nkeep
   end subroutine gmres_recycle_harvest

   !> Gram block G = X(:, 1:kx)^T Y(:, 1:ky) of two tall bases
   !>
   !> - Reduced over fixed row blocks in a fixed order, independent of the thread count
   !>
   !> @param[in]  x   tall basis (n, :)
   !> @param[in]  kx  columns of x used
   !> @param[in]  y   tall basis with the same row count
   !> @param[in]  ky  columns of y used
   !> @param[out] g   Gram block (kx, ky)
   !> @param[in]  nt  OpenMP team size
   subroutine gram_basis(x, kx, y, ky, g, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Tall basis (n, :)
      real(wp), intent(in), contiguous :: x(:, :)
      !> Tall basis with the same row count
      real(wp), intent(in), contiguous :: y(:, :)
      !> Columns of x used
      integer, intent(in) :: kx
      !> Columns of y used
      integer, intent(in) :: ky
      !> Gram block
      real(wp), intent(out) :: g(:, :)
      !> Per-block partial products
      real(wp), allocatable :: partial(:, :, :)
      !> Row count, block count, block index and bounds
      integer :: n, nb, b, first, count

      g = 0.0_wp
      n = size(x, 1)
      if (kx == 0 .or. ky == 0 .or. n == 0) return
      nb = (n - 1)/reduction_block_size + 1
      allocate (partial(kx, ky, nb))
      !$omp parallel do schedule(static) private(first, count) num_threads(merge(nt, 1, n >= parallel_vector_size))
      do b = 1, nb
         first = (b - 1)*reduction_block_size + 1
         count = min(reduction_block_size, n - first + 1)
         call gemm_rows_tn(x(:, 1:kx), y(:, 1:ky), first, count, partial(:, :, b), &
            & 1.0_wp, 0.0_wp)
      end do
      !$omp end parallel do
      do b = 1, nb
         g(1:kx, 1:ky) = g(1:kx, 1:ky) + partial(:, :, b)
      end do
   end subroutine gram_basis

   !> In-place basis combination X(:, 1:kout) := X(:, 1:kx) mx [+ Y(:, 1:ky) my]
   !>
   !> - Row block by row block through a small per-thread buffer
   !>
   !> @param[inout] x     tall basis, overwritten in its first kout columns
   !> @param[in]    kx    columns of x used as input (may be 0)
   !> @param[in]    mx    coefficients (kx, kout)
   !> @param[in]    kout  output columns
   !> @param[in]    nt    OpenMP team size
   !> @param[in]    y     optional second tall basis
   !> @param[in]    ky    columns of y used
   !> @param[in]    my    coefficients (ky, kout)
   subroutine combine_basis(x, kx, mx, kout, nt, y, ky, my)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Tall basis
      real(wp), intent(inout), contiguous :: x(:, :)
      !> Columns of x used as input, may be 0
      integer, intent(in) :: kx
      !> Output columns
      integer, intent(in) :: kout
      !> Coefficients of x
      real(wp), intent(in) :: mx(:, :)
      !> Second basis
      real(wp), intent(in), contiguous, optional :: y(:, :)
      !> Columns of y used
      integer, intent(in), optional :: ky
      !> Coefficients of y
      real(wp), intent(in), optional :: my(:, :)
      !> Rows per block of the buffer
      integer, parameter :: blk = 2048
      !> Per-thread result block
      real(wp), allocatable :: buf(:, :)
      !> Contiguous copies of the coefficients
      real(wp), allocatable :: cx(:, :), cy(:, :)
      !> Row count, columns of y, block index and bounds
      integer :: n, nyc, b, first, count
      !> BLAS accumulation factor: retain the x term if present, otherwise overwrite
      real(wp) :: beta_y

      n = size(x, 1)
      if (kout == 0 .or. n == 0) return
      nyc = 0
      if (present(y) .and. present(ky) .and. present(my)) nyc = ky
      cx = mx(1:kx, 1:kout)
      if (nyc > 0) cy = my(1:nyc, 1:kout)
      beta_y = merge(1.0_wp, 0.0_wp, kx > 0)
      !$omp parallel private(buf, b, first, count) num_threads(merge(nt, 1, n >= parallel_vector_size))
      allocate (buf(blk, kout))
      !$omp do schedule(static)
      do b = 1, (n - 1)/blk + 1
         first = (b - 1)*blk + 1
         count = min(blk, n - first + 1)
         if (kx > 0) call gemm_rows_nn(x(:, 1:kx), first, count, cx, buf, 1.0_wp, 0.0_wp)
         if (nyc > 0) call gemm_rows_nn(y(:, 1:nyc), first, count, cy, buf, 1.0_wp, beta_y)
         x(first:first + count - 1, 1:kout) = buf(1:count, 1:kout)
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel
   end subroutine combine_basis

   !> Thin QR of a small matrix by modified Gram-Schmidt with one reorthogonalization pass
   !>
   !> - Column whose remaining norm is below 1e-10 of the largest column norm is dropped
   !> - Result satisfies a(:, keep(1:nkeep)) = Q R
   !>
   !> @param[in]  a     Matrix (nrow, ncol)
   !> @param[out] q     Orthonormal columns (nrow, ncol), first nkeep valid
   !> @param[out] r     Upper triangular factor (ncol, ncol), leading nkeep block valid
   !> @param[out] keep  Indices of the kept columns
   !> @param[out] nkeep Number of kept columns
   pure subroutine small_qr(a, q, r, keep, nkeep)
      !> Matrix to factor
      real(wp), intent(in) :: a(:, :)
      !> Orthonormal factor
      real(wp), intent(out) :: q(:, :)
      !> Triangular factor
      real(wp), intent(out) :: r(:, :)
      !> Kept column indices
      integer, intent(out) :: keep(:)
      !> Number of kept columns
      integer, intent(out) :: nkeep
      !> Working column
      real(wp) :: w(size(a, 1))
      !> Largest column norm, projection, remaining norm
      real(wp) :: amax, s, nrm
      !> Column, pass and basis indices
      integer :: c, pass, l

      q = 0.0_wp
      r = 0.0_wp
      keep = 0
      nkeep = 0
      amax = 0.0_wp
      do c = 1, size(a, 2)
         amax = max(amax, norm2(a(:, c)))
      end do
      if (amax <= 0.0_wp) return
      do c = 1, size(a, 2)
         w = a(:, c)
         do pass = 1, 2
            do l = 1, nkeep
               s = dot_product(q(:, l), w)
               r(l, nkeep + 1) = r(l, nkeep + 1) + s
               w = w - s*q(:, l)
            end do
         end do
         nrm = norm2(w)
         if (nrm <= 1.0e-10_wp*amax) then
            r(:, nkeep + 1) = 0.0_wp
            cycle
         end if
         nkeep = nkeep + 1
         r(nkeep, nkeep) = nrm
         q(:, nkeep) = w/nrm
         keep(nkeep) = c
      end do
   end subroutine small_qr

   !> Inverse of a nonsingular upper triangular matrix by back substitution
   !>
   !> @param[in]  r    Upper triangular matrix (n, n)
   !> @param[out] rinv Its inverse (n, n), upper triangular
   pure subroutine upper_inverse(r, rinv)
      !> Upper triangular matrix
      real(wp), intent(in) :: r(:, :)
      !> Inverse
      real(wp), intent(out) :: rinv(:, :)
      !> Row, column and summation indices
      integer :: i, c, l
      !> Accumulated right-hand side
      real(wp) :: s

      rinv = 0.0_wp
      do c = 1, size(r, 2)
         do i = c, 1, -1
            s = merge(1.0_wp, 0.0_wp, i == c)
            do l = i + 1, c
               s = s - r(i, l)*rinv(l, c)
            end do
            rinv(i, c) = s/r(i, i)
         end do
      end do
   end subroutine upper_inverse

end submodule moist_model_moz_solver_jfnk_recycle
