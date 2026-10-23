!> Block-diagonal k-space inverse of the 1D MOZ site-site Jacobian
!>
!> - Factorizes `M(k) = I - w*(alpha*c*rho)` per k-point (LU)
!> - Strong preconditioner for the ill-conditioned many-site XRISM Jacobian
!> - Once active, small systems (npair <= 100) use the borrowed spectral filter
!> - Before activation `ready` is false and JFNK runs unpreconditioned
!> - Transforms reuse the spectral filter's uniform radial pair, run outside the OpenMP regions
!> - Adapted from moist_dev/pysofft (44874e7e)
module moist_model_moz_1d_preconditioner_kspace
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: resolve_num_threads
   use moist_math_packed_sym, only: packed_index_lut
   use moist_math_lapack, only: getrf, getrs
   use moist_model_moz_solver_core_preconditioner, only: moist_model_moz_solver_preconditioner_type
   use moist_model_moz_1d_preconditioner_spectral, only: moist_model_moz_1d_preconditioner_spectral_type
   implicit none(type, external)
   private
   public :: moist_model_moz_1d_preconditioner_kspace_type

   !> Site pairs above which the block solve replaces the spectral filter
   integer, parameter :: npair_block_min = 100

   !> Per-k LU factors of the analytic 1D site-site Jacobian
   type, extends(moist_model_moz_solver_preconditioner_type) :: moist_model_moz_1d_preconditioner_kspace_type
      !> Spectral filter providing the radial pair, the r2k filter and the fallback
      type(moist_model_moz_1d_preconditioner_spectral_type), pointer :: spectral => null()
      !> Residual RMS below which the factors are rebuilt and used
      real(wp) :: rmse_activate = 1.0e-3_wp
      !> Rebuild every `refresh` outer iterations once active
      integer :: refresh = 1
      !> OpenMP team size; models should set ctx%get_num_threads(), zero takes omp_get_max_threads
      integer :: nthreads = 0
      !> Number of sites (block dimension)
      integer :: ns = 0
      !> Whether the LU factors have been built
      logical :: valid = .false.
      !> LU factors of M(k) at each k-point (ns, ns, npts)
      real(wp), allocatable :: lu(:, :, :)
      !> LU pivot indices (ns, npts)
      integer, allocatable :: ipiv(:, :)
      !> Packed (is, js) -> pair index lookup for the blocks
      type(packed_index_lut) :: lut
   contains
      !> Borrow the spectral filter
      procedure :: bind => kspace_bind
      !> Whether the driver should rebuild the factors
      procedure :: wants_update => kspace_wants_update
      !> Build and factorize M(k)
      procedure :: update => kspace_update
      !> Whether the factors are built and active at this residual
      procedure :: ready => kspace_ready
      !> Whether `apply` uses the block solve for a number of site pairs
      procedure :: uses_block_solve => kspace_uses_block_solve
      !> Apply M(k)^-1, or the spectral filter for few site pairs
      procedure :: apply => kspace_apply
      !> Status for the solver notes
      procedure :: describe => kspace_describe
      !> Release the factors
      procedure :: release => kspace_release
   end type moist_model_moz_1d_preconditioner_kspace_type

contains

   !> Borrow the spectral filter used on the r2k leg and as the fallback
   !>
   !> @param[inout] self     Preconditioner
   !> @param[in]    spectral Spectral filter; must outlive the preconditioner
   subroutine kspace_bind(self, spectral)
      !> Preconditioner
      class(moist_model_moz_1d_preconditioner_kspace_type), intent(inout) :: self
      !> Spectral filter
      type(moist_model_moz_1d_preconditioner_spectral_type), target, intent(in) :: spectral
      self%spectral => spectral
   end subroutine kspace_bind

   !> Whether the driver should rebuild the factors at this outer iteration
   !>
   !> - True below `rmse_activate` on every `refresh`-th iteration
   !>
   !> @param[in] self  preconditioner
   !> @param[in] rmse  residual RMS of the current iterate
   !> @param[in] iter  stage-local outer iteration
   logical function kspace_wants_update(self, rmse, iter) result(due)
      !> Preconditioner
      class(moist_model_moz_1d_preconditioner_kspace_type), intent(in) :: self
      !> Current residual RMS
      real(wp), intent(in) :: rmse
      !> Outer iteration
      integer, intent(in) :: iter
      due = rmse < self%rmse_activate .and. mod(iter, max(1, self%refresh)) == 0
   end function kspace_wants_update

   !> Whether the linear solve is preconditioned
   !>
   !> - Usable once factorized on the filter's current grid
   !> - Previous residual below `rmse_activate`
   !>
   !> @param[in] self  preconditioner
   !> @param[in] rmse  residual RMS of the previous Newton step
   logical function kspace_ready(self, rmse) result(usable)
      !> Preconditioner
      class(moist_model_moz_1d_preconditioner_kspace_type), intent(in) :: self
      !> Residual RMS of the previous Newton step
      real(wp), intent(in) :: rmse
      usable = .false.
      if (.not. self%valid .or. .not. associated(self%spectral)) return
      if (.not. self%spectral%ready(0.0_wp)) return
      if (size(self%spectral%k) /= size(self%lu, 3)) return
      usable = rmse < self%rmse_activate
   end function kspace_ready

   !> Status for the solver notes: the operator in use, or the activation threshold
   !>
   !> @param[in] self Preconditioner
   !> @param[in] rmse Residual RMS of the previous Newton step
   function kspace_describe(self, rmse) result(text)
      !> Preconditioner
      class(moist_model_moz_1d_preconditioner_kspace_type), intent(in) :: self
      !> Residual RMS of the previous Newton step
      real(wp), intent(in) :: rmse
      !> Status text
      character(len=:), allocatable :: text
      !> Formatted numbers
      character(len=16) :: buf, bmin
      !> Site pairs, as the flat vectors in `apply` see them
      integer :: npair

      if (.not. self%ready(rmse)) then
         write (buf, "(es8.1)") self%rmse_activate
         text = "kspace, off until RMSR < "//trim(adjustl(buf))
         return
      end if
      npair = self%ns*(self%ns + 1)/2
      write (buf, "(i0)") npair
      write (bmin, "(i0)") npair_block_min
      if (self%uses_block_solve(npair)) then
         text = "kspace, block solve on "//trim(buf)//" site pairs"
      else
         text = "kspace, spectral filter for "//trim(buf)//" site pairs (block solve above "// &
            & trim(bmin)//")"
      end if
   end function kspace_describe

   !> Whether `apply` uses the block solve for this many site pairs
   !>
   !> - True above `npair_block_min` site pairs; fewer take the spectral filter once active
   !> - Known before the first `update`, so a driver can check its configuration
   !>
   !> @param[in] self   preconditioner
   !> @param[in] npair  site pairs of the iterate, ns*(ns + 1)/2
   pure logical function kspace_uses_block_solve(self, npair) result(uses)
      !> Preconditioner
      class(moist_model_moz_1d_preconditioner_kspace_type), intent(in) :: self
      !> Site pairs of the iterate
      integer, intent(in) :: npair

      associate (s => self); end associate
      uses = npair > npair_block_min
   end function kspace_uses_block_solve

   !> Apply M(k)^-1 for many-pair systems on the filter's grid, else the spectral filter
   !>
   !> - Spectral filter is a no-op while unbuilt
   !>
   !> @param[inout] self   preconditioner
   !> @param[inout] vflat  flat (npts*npair) vector, modified in place
   subroutine kspace_apply(self, vflat)
      !> Preconditioner
      class(moist_model_moz_1d_preconditioner_kspace_type), intent(inout) :: self
      !> Flat vector, modified in place
      real(wp), contiguous, target, intent(inout) :: vflat(:)

      !> Site pairs of the factors and grid size
      integer :: npair, npts

      if (.not. associated(self%spectral)) return
      if (self%ready(0.0_wp)) then
         npts = size(self%lu, 3)
         npair = self%ns*(self%ns + 1)/2
         if (self%uses_block_solve(npair) .and. mod(size(vflat), npts) == 0) then
            if (size(vflat)/npts == npair) then
               call apply_blocks(self, vflat)
               return
            end if
         end if
      end if
      call self%spectral%apply(vflat)
   end subroutine kspace_apply

   !> Block solve M(k) z = v at every k-point
   !>
   !> - Filtered r2k per pair, solve, symmetrize the packed result, k2r per pair
   !> - A forward transform failure leaves the vector unchanged
   !> - One field-sized k-space buffer; the back transform writes into the vector
   !>
   !> @param[inout] self   preconditioner with built factors
   !> @param[inout] vflat  flat (npts*npair) vector, modified in place
   subroutine apply_blocks(self, vflat)
      !> Preconditioner
      class(moist_model_moz_1d_preconditioner_kspace_type), intent(inout) :: self
      !> Flat vector, modified in place
      real(wp), contiguous, target, intent(inout) :: vflat(:)
      !> Pair-major view of vflat (npts, npair)
      real(wp), pointer :: v(:, :)
      !> k-space field (npts, npair)
      real(wp), allocatable :: vk(:, :)
      !> Per-thread block right-hand side (ns, ns)
      real(wp), allocatable :: rhs(:, :)
      type(error_type), allocatable :: error
      integer :: npts, npair, ns, ik, is, js, ip, info, nt

      npts = size(self%lu, 3)
      ns = self%ns
      npair = ns*(ns + 1)/2
      v(1:npts, 1:npair) => vflat
      allocate (vk(npts, npair))
      call self%spectral%forward(v, vk, error)
      if (allocated(error)) return

      ! Every k-point is an independent ns-by-ns solve
      nt = resolve_num_threads(self%nthreads)
      !$omp parallel num_threads(nt) default(none) &
      !$omp shared(npts, ns, vk, self) private(ik, is, js, ip, info, rhs)
      allocate (rhs(ns, ns))
      !$omp do schedule(static)
      do ik = 1, npts
         do js = 1, ns
            do is = 1, ns
               rhs(is, js) = vk(ik, self%lut%idx(is, js))
            end do
         end do
         call getrs(self%lu(:, :, ik), rhs, self%ipiv(:, ik), info)
         do js = 1, ns
            do is = 1, js
               ip = self%lut%idx(is, js)
               if (is == js) then
                  vk(ik, ip) = rhs(is, js)
               else
                  vk(ik, ip) = 0.5_wp*(rhs(is, js) + rhs(js, is))
               end if
            end do
         end do
      end do
      !$omp end do
      deallocate (rhs)
      !$omp end parallel

      call self%spectral%backward(vk, v, error)
   end subroutine apply_blocks

   !> Build and LU-factorize M(k) = I - w(k)*diag(alpha)*c(k)*diag(rho) per k-point
   !>
   !> - c(k) = FBT(c(r)) - beta*u_lr(k)
   !> - alpha(k) is the volume-weighted running average of dc/dgamma per pair up
   !>   to node nint(r_char/dr + 0.5), r_char = min(r_max, pi/k), clamped to [-1, 0]
   !> - A singular block is replaced by the identity; sets `valid` on success
   !>
   !> @param[inout] self  Preconditioner bound to a built spectral filter
   !> @param[in]    wk    Intramolecular correlation in k-space (npts, npair)
   !> @param[in]    cr    Direct correlation in r-space (npts, npair)
   !> @param[in]    rho   Site number densities (ns), 1/bohr^3
   !> @param[in]    beta  Inverse thermal energy, reciprocal units of uk_lr
   !> @param[in]    uk_lr Long-range potential in k-space (npts, npair)
   !> @param[in]    dcdr  Closure derivative dc/dgamma in r-space (npts, npair)
   !> @param[out]   error Unbound filter, inconsistent shapes or transform failure
   subroutine kspace_update(self, wk, cr, rho, beta, uk_lr, dcdr, error)
      !> Preconditioner
      class(moist_model_moz_1d_preconditioner_kspace_type), intent(inout) :: self
      !> Intramolecular w(k), direct c(r), long-range u_lr(k), closure derivative dc/dgamma; (npts, npair)
      real(wp), intent(in) :: wk(:, :), cr(:, :), uk_lr(:, :), dcdr(:, :)
      !> Site number densities (ns)
      real(wp), intent(in) :: rho(:)
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Packed k-dependent closure scaling and c(k) (npts, npair)
      real(wp), allocatable :: alpha_k(:, :), ck(:, :)
      !> Running volume integrals of dc/dgamma and of the volume
      real(wp), allocatable :: cum_num(:), cum_den(:)
      !> Node at the characteristic length pi/k at each k-point
      integer, allocatable :: j_char(:)
      !> Per-thread unpacked blocks (ns, ns)
      real(wp), allocatable :: w_pt(:, :), ck_rho(:, :)
      real(wp) :: pi_val
      integer :: npts, npair, ns, ik, is, js, ip, j, info, nt

      self%valid = .false.
      if (.not. associated(self%spectral)) then
         call fatal_error(error, "kspace preconditioner needs a bound spectral filter")
         return
      end if
      if (.not. self%spectral%ready(0.0_wp)) then
         call fatal_error(error, "kspace preconditioner needs a built spectral filter")
         return
      end if
      npts = size(self%spectral%k)
      ns = size(rho)
      npair = ns*(ns + 1)/2
      if (ns == 0 .or. any(shape(cr) /= [npts, npair]) .or. any(shape(wk) /= shape(cr)) .or. &
         & any(shape(uk_lr) /= shape(cr)) .or. any(shape(dcdr) /= shape(cr))) then
         call fatal_error(error, "kspace preconditioner fields must have shape (npts, ns*(ns+1)/2)")
         return
      end if

      if (self%ns /= ns .or. .not. allocated(self%lu)) then
         if (allocated(self%lu)) deallocate (self%lu)
         if (allocated(self%ipiv)) deallocate (self%ipiv)
         allocate (self%lu(ns, ns, npts), self%ipiv(ns, npts))
         self%ns = ns
         call self%lut%init(ns)
      else if (size(self%lu, 3) /= npts) then
         deallocate (self%lu, self%ipiv)
         allocate (self%lu(ns, ns, npts), self%ipiv(ns, npts))
      end if

      ! Sample the running average of dc/dgamma at the length scale pi/k,
      ! approximating the diagonal of FBT . diag(dc/dgamma) . FBT^-1; the node
      ! rule is pysofft's, including its rounding of ties to the upper node
      pi_val = acos(-1.0_wp)
      allocate (j_char(npts))
      do ik = 1, npts
         j_char(ik) = min(npts, max(1, nint(min(self%spectral%r(npts), pi_val/self%spectral%k(ik)) &
            & /self%spectral%dr + 0.5_wp)))
      end do
      allocate (alpha_k(npts, npair), cum_num(npts), cum_den(npts))
      do ip = 1, npair
         cum_num(1) = dcdr(1, ip)*self%spectral%r(1)**2*self%spectral%w(1)
         cum_den(1) = self%spectral%r(1)**2*self%spectral%w(1)
         do j = 2, npts
            cum_num(j) = cum_num(j - 1) + dcdr(j, ip)*self%spectral%r(j)**2*self%spectral%w(j)
            cum_den(j) = cum_den(j - 1) + self%spectral%r(j)**2*self%spectral%w(j)
         end do
         do ik = 1, npts
            alpha_k(ik, ip) = max(-1.0_wp, min(0.0_wp, cum_num(j_char(ik))/cum_den(j_char(ik))))
         end do
      end do

      allocate (ck(npts, npair))
      call self%spectral%trafo%fbt_r2k_all(cr, ck, error)
      if (allocated(error)) return
      ck = ck - beta*uk_lr

      nt = resolve_num_threads(self%nthreads)
      !$omp parallel num_threads(nt) default(none) &
      !$omp shared(npts, ns, alpha_k, wk, ck, rho, self) &
      !$omp private(ik, is, js, ip, info, w_pt, ck_rho)
      allocate (w_pt(ns, ns), ck_rho(ns, ns))
      !$omp do schedule(static)
      do ik = 1, npts
         do js = 1, ns
            do is = 1, ns
               ip = self%lut%idx(is, js)
               w_pt(is, js) = wk(ik, ip)
               ck_rho(is, js) = alpha_k(ik, ip)*ck(ik, ip)*rho(js)
            end do
         end do
         self%lu(:, :, ik) = -matmul(w_pt, ck_rho)
         do is = 1, ns
            self%lu(is, is, ik) = self%lu(is, is, ik) + 1.0_wp
         end do
         call getrf(self%lu(:, :, ik), self%ipiv(:, ik), info)
         if (info /= 0) then
            self%lu(:, :, ik) = 0.0_wp
            do is = 1, ns
               self%lu(is, is, ik) = 1.0_wp
               self%ipiv(is, ik) = is
            end do
         end if
      end do
      !$omp end do
      deallocate (w_pt, ck_rho)
      !$omp end parallel

      self%valid = .true.
   end subroutine kspace_update

   !> Release the factors, safe on repeated calls; the spectral filter stays bound
   !>
   !> @param[inout] self Preconditioner
   subroutine kspace_release(self)
      !> Preconditioner
      class(moist_model_moz_1d_preconditioner_kspace_type), intent(inout) :: self
      if (allocated(self%lu)) deallocate (self%lu)
      if (allocated(self%ipiv)) deallocate (self%ipiv)
      self%valid = .false.
      self%ns = 0
   end subroutine kspace_release

end module moist_model_moz_1d_preconditioner_kspace
