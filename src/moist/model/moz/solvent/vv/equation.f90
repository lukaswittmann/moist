!> 1D VV equation on a uniform radial pair: the map c_s(r) -> t_s(r) and its linearization
!>
!> - Per k: h = (I - w c R)^-1 w c w + h0, R = diag(rho)
!>   w = omega and h0 = 0 for RISM; w = omega + rho zeta and h0 = zeta for DRISM
!> - Ng split: the map takes short-range c_s(r); c(k) = c_s(k) - beta u_lr(k)
!>   and t_s(k) = h(k) - c(k) - beta u_lr(k)
!>   The nonlinear driver iterates on t_s and supplies c_s through the closure
!> - Linearization dt(k) = L dc M - dc with L = A^-1 w, M = w + R (h - h0), A = I - w c R
!> - Coarse matrix of the residual Jacobian I - T' c' on Gaussian profiles
!>
!> Each parallel region uses one transform clone per thread and num_threads(nt).
!> Adapted from moist_dev/pysofft (44874e7e), model/rism/1d/equation.f90
module moist_model_moz_solvent_vv_equation
   !$ use omp_lib, only: omp_get_thread_num
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: resolve_num_threads
   use moist_math_packed_sym, only: packed_index_lut, npair_from_ns
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, transform_dst4
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo
   use moist_math_lapack_getrf, only: wrap_getrf
   use moist_math_lapack_getrs, only: wrap_getrs
   use moist_model_moz_solvent_vv_kernel, only: vv_kblock, vv_build_blocks, vv_solve_block, vv_solve_point
   use moist_model_moz_fields, only: field_num_threads
   implicit none(type, external)
   private

   public :: vv_equation_type

   !> One radial transform owned by one thread
   type :: trafo_slot_type
      !> Transform with its own scratch
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
   end type trafo_slot_type

   !> Fixed VV data, the per-thread transforms and the cached linearization
   type :: vv_equation_type
      !> OpenMP team size; zero takes omp_get_max_threads, small fields run serially
      integer :: nthreads = 0
      !> Sites, grid points, packed pairs and k blocks; 0 before setup
      integer, private :: ns = 0, npts = 0, npair = 0, nblock = 0
      !> Inverse thermal energy, 1/Hartree
      real(wp), private :: beta = 0.0_wp
      !> Site number densities, 1/bohr^3 (ns)
      real(wp), allocatable, private :: rho(:)
      !> Long-range reciprocal potential, Hartree bohr^3, packed (npts, npair)
      real(wp), allocatable, private :: uk_lr(:, :)
      !> Fixed additive total-correlation correction h0(k), packed (npts, npair)
      real(wp), allocatable, private :: h_offset(:, :)
      !> w(k), (vv_kblock, ns, ns, nblock)
      real(wp), allocatable, private :: w_blk(:, :, :, :)
      !> L(k) and M(k) of the latest linearization, same layout
      real(wp), allocatable, private :: left_blk(:, :, :, :), right_blk(:, :, :, :)
      !> Packed-index lookup for ns sites
      type(packed_index_lut), private :: lut
      !> Transform of the radial pair and its per-thread clones
      type(trafo_slot_type), allocatable, private :: slots(:)
      !> k-space scratch, packed (npts, npair)
      real(wp), allocatable, private :: ck_work(:, :), hk_work(:, :), tk_work(:, :)
   contains
      !> Bind the radial pair and the fixed solve data
      procedure :: setup => equation_setup
      !> Whether setup succeeded
      procedure :: ready => equation_ready
      !> Number of sites
      procedure :: nsite => equation_nsite
      !> Map c_s(r) -> t_s(r), optionally h(k)
      procedure :: map => equation_map
      !> Cache L and M at c_s(r)
      procedure :: linearize => equation_linearize
      !> Whether a linearization is cached
      procedure :: linearized => equation_linearized
      !> Apply the cached linearization to dc(r)
      procedure :: jvp => equation_jvp
      !> Apply the transpose of the cached linearization
      procedure :: adjoint => equation_adjoint
      !> Coarse matrix of the residual Jacobian on radial profiles
      procedure :: coarse_matrix => equation_coarse_matrix
      !> Grow the per-thread transform clones to nt
      procedure, private :: ensure_slots => equation_ensure_slots
   end type vv_equation_type

contains

   !> Bind a uniform DST-IV radial pair and the fixed solve data
   !>
   !> @param[inout] self   equation; a previous setup and linearization are dropped
   !> @param[in]    rgrid  real-space radial grid
   !> @param[in]    kgrid  reciprocal radial grid of the pair
   !> @param[in]    wk     effective bulk w(k), packed (npts, npair)
   !> @param[in]    rho    site number densities, 1/bohr^3 (ns)
   !> @param[in]    beta   inverse thermal energy, 1/Hartree
   !> @param[in]    uk_lr  long-range reciprocal potential, Hartree bohr^3, packed (npts, npair)
   !> @param[out]   error  not a uniform DST-IV pair, inconsistent shapes or transform failure
   !> @param[in]    h_offset optional fixed additive h(k) correction; zero by default
   subroutine equation_setup(self, rgrid, kgrid, wk, rho, beta, uk_lr, error, h_offset)
      !> Equation
      class(vv_equation_type), intent(inout) :: self
      !> Real-space and reciprocal radial grids
      type(moist_math_grid_radial_type), intent(in) :: rgrid, kgrid
      !> Packed w(k) and long-range potential
      real(wp), intent(in) :: wk(:, :), uk_lr(:, :)
      !> Site number densities
      real(wp), intent(in) :: rho(:)
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: h_offset(:, :)
      integer :: npts, npair

      self%ns = 0
      if (allocated(self%slots)) deallocate (self%slots)
      if (allocated(self%left_blk)) deallocate (self%left_blk, self%right_blk)
      if (rgrid%transform /= transform_dst4 .or. kgrid%transform /= transform_dst4 .or. &
         & .not. allocated(rgrid%r) .or. .not. allocated(kgrid%r)) then
         call fatal_error(error, "1D VV equation needs the uniform DST-IV radial pair")
         return
      end if
      npts = size(rgrid%r)
      npair = npair_from_ns(size(rho))
      if (size(rho) < 1 .or. size(kgrid%r) /= npts .or. any(shape(wk) /= [npts, npair]) .or. &
         & any(shape(uk_lr) /= [npts, npair])) then
         call fatal_error(error, "1D VV equation data must have shape (npts, ns*(ns+1)/2)")
         return
      end if
      if (present(h_offset)) then
         if (any(shape(h_offset) /= shape(wk))) then
            call fatal_error(error, "1D VV h offset must match the packed intramolecular table")
            return
         end if
      end if
      allocate (self%slots(1))
      call new_radial_trafo(self%slots(1)%trafo, rgrid, kgrid, error)
      if (allocated(error)) then
         deallocate (self%slots)
         return
      end if
      call self%lut%init(size(rho))
      call vv_build_blocks(wk, self%lut, self%w_blk)
      self%rho = rho
      self%beta = beta
      self%uk_lr = uk_lr
      if (allocated(self%h_offset)) deallocate (self%h_offset)
      allocate (self%h_offset(npts, npair), source=0.0_wp)
      if (present(h_offset)) self%h_offset = h_offset
      self%npts = npts
      self%npair = npair
      self%nblock = size(self%w_blk, 4)
      if (allocated(self%ck_work)) deallocate (self%ck_work, self%hk_work, self%tk_work)
      allocate (self%ck_work(npts, npair), self%hk_work(npts, npair), self%tk_work(npts, npair))
      self%ns = size(rho)
   end subroutine equation_setup

   !> Whether setup succeeded
   !>
   !> @param[in] self Equation
   pure logical function equation_ready(self) result(ready)
      !> Equation
      class(vv_equation_type), intent(in) :: self
      ready = self%ns > 0
   end function equation_ready

   !> Number of sites; 0 before setup
   !>
   !> @param[in] self Equation
   pure integer function equation_nsite(self) result(ns)
      !> Equation
      class(vv_equation_type), intent(in) :: self
      ns = self%ns
   end function equation_nsite

   !> Whether a linearization is cached
   !>
   !> @param[in] self Equation
   pure logical function equation_linearized(self) result(cached)
      !> Equation
      class(vv_equation_type), intent(in) :: self
      cached = self%ns > 0 .and. allocated(self%left_blk)
   end function equation_linearized

   !> Grow the per-thread transform clones to at least nt
   !>
   !> @param[inout] self Equation after setup
   !> @param[in]    nt   team size
   subroutine equation_ensure_slots(self, nt)
      !> Equation
      class(vv_equation_type), intent(inout) :: self
      !> Team size
      integer, intent(in) :: nt
      type(trafo_slot_type), allocatable :: slots(:)
      integer :: it, n

      n = size(self%slots)
      if (n >= nt) return
      allocate (slots(nt))
      do it = 1, n
         call move_alloc(self%slots(it)%trafo, slots(it)%trafo)
      end do
      do it = n + 1, nt
         allocate (slots(it)%trafo, source=slots(1)%trafo)
      end do
      call move_alloc(slots, self%slots)
   end subroutine equation_ensure_slots

   !> Map the short-range direct correlation c_s(r) to t_s(r)
   !>
   !> - Blocked solve; the rare points that need a row interchange are re-solved
   !>   with LAPACK, chosen per point so the map stays smooth for difference quotients
   !> - One parallel region: forward transforms, per-k solve, inverse transforms
   !>
   !> @param[inout] self   equation after setup
   !> @param[in]    cr     c_s(r), packed (npts, npair)
   !> @param[out]   tr     t_s(r), packed (npts, npair)
   !> @param[out]   error  singular point solve or transform failure
   !> @param[out]   hk     optional physical h(k), packed (npts, npair)
   subroutine equation_map(self, cr, tr, error, hk)
      !> Equation
      class(vv_equation_type), intent(inout), target :: self
      !> Short-range direct correlation
      real(wp), intent(in) :: cr(:, :)
      !> Short-range indirect correlation
      real(wp), intent(out) :: tr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional total correlation in k-space
      real(wp), intent(out), optional :: hk(:, :)
      real(wp), pointer :: ck(:, :), hkw(:, :), tk(:, :)
      real(wp), allocatable :: ck_blk(:, :, :), hk_blk(:, :, :)
      real(wp), allocatable :: w_pt(:, :), ck_pt(:, :), hk_pt(:, :)
      logical :: need_pivot(vv_kblock)
      integer :: nt, it, ip, ib, ik, ik0, nk, kk, is, js, info, error_flag, trafo_flag
      type(error_type), allocatable :: terr

      if (.not. self%ready()) then
         call fatal_error(error, "1D VV equation is not set up")
         return
      end if
      if (any(shape(cr) /= [self%npts, self%npair]) .or. any(shape(tr) /= shape(cr))) then
         call fatal_error(error, "1D VV map fields must have shape (npts, ns*(ns+1)/2)")
         return
      end if
      nt = field_num_threads(self%npts*self%npair, resolve_num_threads(self%nthreads))
      call self%ensure_slots(nt)
      ck => self%ck_work
      hkw => self%hk_work
      tk => self%tk_work
      error_flag = 0
      trafo_flag = 0

      !$omp parallel num_threads(nt) default(none) &
      !$omp shared(self, cr, tr, ck, hkw, tk, error_flag, trafo_flag) &
      !$omp private(it, ip, ib, ik, ik0, nk, kk, is, js, info, need_pivot, ck_blk, hk_blk, &
      !$omp         w_pt, ck_pt, hk_pt, terr)
      it = 1
      !$ it = omp_get_thread_num() + 1
      allocate (ck_blk(vv_kblock, self%ns, self%ns), hk_blk(vv_kblock, self%ns, self%ns))
      allocate (w_pt(self%ns, self%ns), ck_pt(self%ns, self%ns), hk_pt(self%ns, self%ns))
      ! Lanes past the end of the grid stay zero: w is zero there, A = I
      ck_blk = 0.0_wp

      !$omp do schedule(static) reduction(max:trafo_flag)
      do ip = 1, self%npair
         call self%slots(it)%trafo%fbt_r2k(cr(:, ip), ck(:, ip), terr)
         if (allocated(terr)) then
            trafo_flag = 1
            deallocate (terr)
         end if
      end do
      !$omp end do

      !$omp do schedule(static) reduction(max:error_flag)
      do ib = 1, self%nblock
         ik0 = (ib - 1)*vv_kblock
         nk = min(vv_kblock, self%npts - ik0)
         do js = 1, self%ns
            do is = 1, js
               ip = self%lut%idx(is, js)
               ! OZ uses physical c(k); retain c_s(k) in ck for the Ng relation.
               ck_blk(1:nk, is, js) = ck(ik0 + 1:ik0 + nk, ip) - &
                  & self%beta*self%uk_lr(ik0 + 1:ik0 + nk, ip)
               ck_blk(1:nk, js, is) = ck_blk(1:nk, is, js)
            end do
         end do
         call vv_solve_block(self%ns, self%w_blk(:, :, :, ib), ck_blk, self%rho, hk_blk, need_pivot)
         do kk = 1, nk
            if (.not. need_pivot(kk)) cycle
            w_pt = self%w_blk(kk, :, :, ib)
            ck_pt = ck_blk(kk, :, :)
            call vv_solve_point(self%ns, w_pt, ck_pt, self%rho, hk_pt, info)
            if (info /= 0) error_flag = max(error_flag, ik0 + kk)
            hk_blk(kk, :, :) = hk_pt
         end do
         do js = 1, self%ns
            do is = 1, js
               ip = self%lut%idx(is, js)
               do kk = 1, nk
                  ik = ik0 + kk
                  hkw(ik, ip) = hk_blk(kk, is, js) + self%h_offset(ik, ip)
                  tk(ik, ip) = hkw(ik, ip) - ck(ik, ip)
               end do
            end do
         end do
      end do
      !$omp end do

      !$omp do schedule(static) reduction(max:trafo_flag)
      do ip = 1, self%npair
         call self%slots(it)%trafo%fbt_k2r(tk(:, ip), tr(:, ip), terr)
         if (allocated(terr)) then
            trafo_flag = 1
            deallocate (terr)
         end if
      end do
      !$omp end do
      !$omp end parallel

      if (trafo_flag /= 0) then
         call fatal_error(error, "Radial transform failed in the 1D VV equation")
      else if (error_flag /= 0) then
         call fatal_error(error, "Matrix factorization or solve failed in the 1D VV equation")
      else if (present(hk)) then
         hk = hkw
      end if
   end subroutine equation_map

   !> Cache L = A^-1 w and M = w + R (h - h0) at c_s(r)
   !>
   !> - From A (h - h0) = w c w and dA = -w dc R: dh = L dc M
   !> - h0 is fixed (zero for RISM, zeta for DRISM); site densities may differ
   !> - Point solves with LAPACK partial pivoting, once per Newton step
   !>
   !> @param[inout] self   equation after setup
   !> @param[in]    cr     c_s(r) at the linearization point, packed (npts, npair)
   !> @param[out]   error  singular or non-finite point solve, or transform failure
   subroutine equation_linearize(self, cr, error)
      !> Equation
      class(vv_equation_type), intent(inout), target :: self
      !> Short-range direct correlation
      real(wp), intent(in) :: cr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), pointer :: ck(:, :), left(:, :, :, :), right(:, :, :, :)
      real(wp), allocatable :: w_pt(:, :), c_pt(:, :), a_pt(:, :), rhs(:, :)
      integer, allocatable :: ipiv(:)
      integer :: nt, it, ik, ip, ib, kk, is, js, info, error_flag, trafo_flag
      type(error_type), allocatable :: terr

      if (.not. self%ready()) then
         call fatal_error(error, "1D VV equation is not set up")
         return
      end if
      if (any(shape(cr) /= [self%npts, self%npair])) then
         call fatal_error(error, "1D VV linearization field must have shape (npts, ns*(ns+1)/2)")
         return
      end if
      if (.not. allocated(self%left_blk)) then
         allocate (self%left_blk, mold=self%w_blk)
         allocate (self%right_blk, mold=self%w_blk)
      end if
      nt = field_num_threads(self%npts*self%npair, resolve_num_threads(self%nthreads))
      call self%ensure_slots(nt)
      ck => self%ck_work
      left => self%left_blk
      right => self%right_blk
      error_flag = 0
      trafo_flag = 0

      !$omp parallel num_threads(nt) default(none) &
      !$omp shared(self, cr, ck, left, right, error_flag, trafo_flag) &
      !$omp private(it, ik, ip, ib, kk, is, js, info, ipiv, w_pt, c_pt, a_pt, rhs, terr)
      it = 1
      !$ it = omp_get_thread_num() + 1
      allocate (ipiv(self%ns), w_pt(self%ns, self%ns), c_pt(self%ns, self%ns), a_pt(self%ns, self%ns), &
         & rhs(self%ns, 2*self%ns))

      !$omp do schedule(static)
      do ib = 1, self%nblock
         left(:, :, :, ib) = 0.0_wp
         right(:, :, :, ib) = 0.0_wp
      end do
      !$omp end do nowait

      !$omp do schedule(static) reduction(max:trafo_flag)
      do ip = 1, self%npair
         call self%slots(it)%trafo%fbt_r2k(cr(:, ip), ck(:, ip), terr)
         if (allocated(terr)) then
            trafo_flag = 1
            deallocate (terr)
         end if
         ck(:, ip) = ck(:, ip) - self%beta*self%uk_lr(:, ip)
      end do
      !$omp end do

      !$omp do schedule(static) reduction(max:error_flag)
      do ik = 1, self%npts
         kk = mod(ik - 1, vv_kblock) + 1
         ib = (ik - 1)/vv_kblock + 1
         w_pt = self%w_blk(kk, :, :, ib)
         do js = 1, self%ns
            do is = 1, js
               ip = self%lut%idx(is, js)
               c_pt(is, js) = ck(ik, ip)
               c_pt(js, is) = ck(ik, ip)
            end do
         end do
         ! A = I - w (c R) and the right-hand sides [w c w | w]
         a_pt = matmul(w_pt, c_pt)
         do js = 1, self%ns
            a_pt(:, js) = -a_pt(:, js)*self%rho(js)
            a_pt(js, js) = a_pt(js, js) + 1.0_wp
         end do
         rhs(:, 1:self%ns) = matmul(w_pt, matmul(c_pt, w_pt))
         rhs(:, self%ns + 1:2*self%ns) = w_pt
         call wrap_getrf(a_pt, ipiv, info)
         if (info == 0) call wrap_getrs(a_pt, rhs, ipiv, info)
         if (info == 0 .and. .not. all(abs(rhs) <= huge(1.0_wp))) info = 1
         if (info /= 0) then
            error_flag = max(error_flag, ik)
            cycle
         end if
         left(kk, :, :, ib) = rhs(:, self%ns + 1:2*self%ns)
         do js = 1, self%ns
            right(kk, :, js, ib) = w_pt(:, js) + self%rho(:)*rhs(:, js)
         end do
      end do
      !$omp end do
      !$omp end parallel

      if (trafo_flag /= 0) then
         call fatal_error(error, "Radial transform failed in the 1D VV linearization")
      else if (error_flag /= 0) then
         call fatal_error(error, "Singular or non-finite k-point solve in the 1D VV linearization")
      end if
      if (allocated(error)) deallocate (self%left_blk, self%right_blk)
   end subroutine equation_linearize

   !> Apply the cached linearization: dt(r) = k2r(L dc(k) M - dc(k))
   !>
   !> @param[inout] self   equation with a cached linearization
   !> @param[in]    dcr    direction dc(r), packed (npts, npair)
   !> @param[out]   dtr    response dt(r), packed (npts, npair)
   !> @param[out]   error  no linearization or transform failure
   subroutine equation_jvp(self, dcr, dtr, error)
      !> Equation
      class(vv_equation_type), intent(inout), target :: self
      !> Direction
      real(wp), intent(in) :: dcr(:, :)
      !> Response
      real(wp), intent(out) :: dtr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), pointer :: dck(:, :), dtk(:, :)
      real(wp), allocatable :: dc_blk(:, :, :), y_blk(:, :, :)
      real(wp) :: acc(vv_kblock)
      integer :: nt, it, ip, ib, ik0, nk, is, js, l, trafo_flag
      type(error_type), allocatable :: terr

      if (.not. self%linearized()) then
         call fatal_error(error, "1D VV equation is not linearized")
         return
      end if
      if (any(shape(dcr) /= [self%npts, self%npair]) .or. any(shape(dtr) /= shape(dcr))) then
         call fatal_error(error, "1D VV Jacobian fields must have shape (npts, ns*(ns+1)/2)")
         return
      end if
      nt = field_num_threads(self%npts*self%npair, resolve_num_threads(self%nthreads))
      call self%ensure_slots(nt)
      dck => self%ck_work
      dtk => self%tk_work
      trafo_flag = 0

      !$omp parallel num_threads(nt) default(none) &
      !$omp shared(self, dcr, dtr, dck, dtk, trafo_flag) &
      !$omp private(it, ip, ib, ik0, nk, is, js, l, dc_blk, y_blk, acc, terr)
      it = 1
      !$ it = omp_get_thread_num() + 1
      allocate (dc_blk(vv_kblock, self%ns, self%ns), y_blk(vv_kblock, self%ns, self%ns))
      dc_blk = 0.0_wp

      !$omp do schedule(static) reduction(max:trafo_flag)
      do ip = 1, self%npair
         call self%slots(it)%trafo%fbt_r2k(dcr(:, ip), dck(:, ip), terr)
         if (allocated(terr)) then
            trafo_flag = 1
            deallocate (terr)
         end if
      end do
      !$omp end do

      !$omp do schedule(static)
      do ib = 1, self%nblock
         ik0 = (ib - 1)*vv_kblock
         nk = min(vv_kblock, self%npts - ik0)
         do js = 1, self%ns
            do is = 1, js
               ip = self%lut%idx(is, js)
               dc_blk(1:nk, is, js) = dck(ik0 + 1:ik0 + nk, ip)
               dc_blk(1:nk, js, is) = dck(ik0 + 1:ik0 + nk, ip)
            end do
         end do
         ! Y = dc M
         y_blk = 0.0_wp
         do js = 1, self%ns
            do l = 1, self%ns
               do is = 1, self%ns
                  y_blk(:, is, js) = y_blk(:, is, js) + dc_blk(:, is, l)*self%right_blk(:, l, js, ib)
               end do
            end do
         end do
         ! dt(k) = (L Y)(is, js) - dc(k), upper triangle
         do js = 1, self%ns
            do is = 1, js
               ip = self%lut%idx(is, js)
               acc = 0.0_wp
               do l = 1, self%ns
                  acc = acc + self%left_blk(:, is, l, ib)*y_blk(:, l, js)
               end do
               dtk(ik0 + 1:ik0 + nk, ip) = acc(1:nk) - dck(ik0 + 1:ik0 + nk, ip)
            end do
         end do
      end do
      !$omp end do

      !$omp do schedule(static) reduction(max:trafo_flag)
      do ip = 1, self%npair
         call self%slots(it)%trafo%fbt_k2r(dtk(:, ip), dtr(:, ip), terr)
         if (allocated(terr)) then
            trafo_flag = 1
            deallocate (terr)
         end if
      end do
      !$omp end do
      !$omp end parallel

      if (trafo_flag /= 0) call fatal_error(error, "Radial transform failed in the 1D VV Jacobian product")
   end subroutine equation_jvp

   !> Apply the transpose of the cached linearization in the packed Euclidean product
   !>
   !> - The forward product unpacks dc(k) to a symmetric matrix and keeps the
   !>   upper triangle of L dc M, so per k the transpose places u(k) in the upper
   !>   triangle U and folds Z = L^T U M^T: z_p = Z_ij + Z_ji off the diagonal
   !> - dc(r) = r2k^T(z - u(k)) with u(k) = k2r^T u(r)
   !>
   !> @param[inout] self   equation with a cached linearization
   !> @param[in]    dtr    direction u(r), packed (npts, npair)
   !> @param[out]   dcr    transpose product, packed (npts, npair)
   !> @param[out]   error  no linearization or transform failure
   subroutine equation_adjoint(self, dtr, dcr, error)
      !> Equation
      class(vv_equation_type), intent(inout), target :: self
      !> Direction
      real(wp), intent(in) :: dtr(:, :)
      !> Transpose product
      real(wp), intent(out) :: dcr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), pointer :: uk(:, :), zk(:, :)
      real(wp), allocatable :: u_blk(:, :, :), y_blk(:, :, :)
      real(wp) :: acc(vv_kblock)
      integer :: nt, it, ip, ib, ik0, nk, is, js, l, trafo_flag
      type(error_type), allocatable :: terr

      if (.not. self%linearized()) then
         call fatal_error(error, "1D VV equation is not linearized")
         return
      end if
      if (any(shape(dtr) /= [self%npts, self%npair]) .or. any(shape(dcr) /= shape(dtr))) then
         call fatal_error(error, "1D VV Jacobian fields must have shape (npts, ns*(ns+1)/2)")
         return
      end if
      nt = field_num_threads(self%npts*self%npair, resolve_num_threads(self%nthreads))
      call self%ensure_slots(nt)
      uk => self%tk_work
      zk => self%ck_work
      trafo_flag = 0

      !$omp parallel num_threads(nt) default(none) &
      !$omp shared(self, dtr, dcr, uk, zk, trafo_flag) &
      !$omp private(it, ip, ib, ik0, nk, is, js, l, u_blk, y_blk, acc, terr)
      it = 1
      !$ it = omp_get_thread_num() + 1
      allocate (u_blk(vv_kblock, self%ns, self%ns), y_blk(vv_kblock, self%ns, self%ns))
      u_blk = 0.0_wp

      !$omp do schedule(static) reduction(max:trafo_flag)
      do ip = 1, self%npair
         call self%slots(it)%trafo%fbt_k2r_adj(dtr(:, ip), uk(:, ip), terr)
         if (allocated(terr)) then
            trafo_flag = 1
            deallocate (terr)
         end if
      end do
      !$omp end do

      !$omp do schedule(static)
      do ib = 1, self%nblock
         ik0 = (ib - 1)*vv_kblock
         nk = min(vv_kblock, self%npts - ik0)
         ! U = u(k) on the upper triangle, zero below
         do js = 1, self%ns
            do is = 1, js
               u_blk(1:nk, is, js) = uk(ik0 + 1:ik0 + nk, self%lut%idx(is, js))
            end do
         end do
         ! Y = U M^T
         y_blk = 0.0_wp
         do js = 1, self%ns
            do l = 1, self%ns
               do is = 1, l
                  y_blk(:, is, js) = y_blk(:, is, js) + u_blk(:, is, l)*self%right_blk(:, js, l, ib)
               end do
            end do
         end do
         ! z_p = (L^T Y)_ij + (L^T Y)_ji, minus u(k)
         do js = 1, self%ns
            do is = 1, js
               ip = self%lut%idx(is, js)
               acc = 0.0_wp
               do l = 1, self%ns
                  acc = acc + self%left_blk(:, l, is, ib)*y_blk(:, l, js)
               end do
               if (is /= js) then
                  do l = 1, self%ns
                     acc = acc + self%left_blk(:, l, js, ib)*y_blk(:, l, is)
                  end do
               end if
               zk(ik0 + 1:ik0 + nk, ip) = acc(1:nk) - uk(ik0 + 1:ik0 + nk, ip)
            end do
         end do
      end do
      !$omp end do

      !$omp do schedule(static) reduction(max:trafo_flag)
      do ip = 1, self%npair
         call self%slots(it)%trafo%fbt_r2k_adj(zk(:, ip), dcr(:, ip), terr)
         if (allocated(terr)) then
            trafo_flag = 1
            deallocate (terr)
         end if
      end do
      !$omp end do
      !$omp end parallel

      if (trafo_flag /= 0) call fatal_error(error, "Radial transform failed in the 1D VV transpose product")
   end subroutine equation_adjoint

   !> Coarse matrix E = Z^T Omega J Z of the residual Jacobian J = I - T' c'
   !>
   !> - Coarse space z_jp = g_j e_p: one radial profile on one site pair
   !> - With d_jp = r2k(c'_p g_j) and the pair response K_qp = L_sa M_bt + L_sb M_at
   !>   (L_sa M_at for a = b) of pair p = (a, b) at q = (s, t), and phi_i = k2r^T(Omega g_i),
   !>   E_(iq,jp) = delta_qp (<g_i, Omega g_j> + sum_k phi_i d_jp) - sum_k phi_i d_jp K_qp
   !> - phi_i is Gaussian in k, so the k sums stop where every phi_i has decayed
   !> - Each column pair is summed by one thread in a fixed order: E does not
   !>   depend on the team size
   !>
   !> @param[inout] self      equation with a cached linearization
   !> @param[in]    profile   radial profiles g_j, (npts, nfun)
   !> @param[in]    wprofile  weighted profiles Omega g_j, (npts, nfun)
   !> @param[in]    dcdt      closure derivative c' at the linearization point, packed (npts, npair)
   !> @param[out]   coarse    E, (nfun*npair, nfun*npair); (j, p) has index (p - 1) nfun + j
   !> @param[out]   error     no linearization, inconsistent shapes or transform failure
   subroutine equation_coarse_matrix(self, profile, wprofile, dcdt, coarse, error)
      !> Equation
      class(vv_equation_type), intent(inout), target :: self
      !> Radial profiles and weighted profiles
      real(wp), intent(in) :: profile(:, :), wprofile(:, :)
      !> Closure derivative
      real(wp), intent(in) :: dcdt(:, :)
      !> Coarse matrix
      real(wp), intent(out) :: coarse(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> phi_i below this fraction of its peak is dropped from the k sums
      real(wp), parameter :: phi_cutoff = 1.0e-16_wp
      integer, allocatable :: site1(:), site2(:)
      real(wp), allocatable :: phi(:, :), overlap(:, :)
      real(wp), allocatable :: dk(:, :), work(:), acc(:, :, :), wv(:, :, :), wsum(:, :)
      real(wp) :: kresp(vv_kblock), phi_peak(size(profile, 2))
      integer :: nt, it, nf, p, q, a, b, s, t, i, j, ib, ik0, nk, kcut, nblk, trafo_flag
      type(error_type), allocatable :: terr

      if (.not. self%linearized()) then
         call fatal_error(error, "1D VV equation is not linearized")
         return
      end if
      nf = size(profile, 2)
      if (size(profile, 1) /= self%npts .or. any(shape(wprofile) /= shape(profile)) .or. &
         & any(shape(dcdt) /= [self%npts, self%npair]) .or. &
         & any(shape(coarse) /= [nf*self%npair, nf*self%npair])) then
         call fatal_error(error, "1D VV coarse matrix arguments have inconsistent shapes")
         return
      end if
      nt = field_num_threads(self%npts*self%npair, resolve_num_threads(self%nthreads))
      call self%ensure_slots(nt)

      allocate (site1(self%npair), site2(self%npair))
      do t = 1, self%ns
         do s = 1, t
            site1(self%lut%idx(s, t)) = s
            site2(self%lut%idx(s, t)) = t
         end do
      end do
      allocate (phi(self%npts, nf), overlap(nf, nf))
      do i = 1, nf
         call self%slots(1)%trafo%fbt_k2r_adj(wprofile(:, i), phi(:, i), error)
         if (allocated(error)) return
         do j = 1, nf
            overlap(i, j) = dot_product(wprofile(:, i), profile(:, j))
         end do
      end do
      ! All transformed profiles are complete and fixed during the scan.
      phi_peak = maxval(abs(phi), dim=1)
      kcut = 1
      do ik0 = self%npts, 1, -1
         if (any(abs(phi(ik0, :)) > phi_cutoff*phi_peak)) then
            kcut = ik0
            exit
         end if
      end do
      nblk = (kcut - 1)/vv_kblock + 1
      trafo_flag = 0

      !$omp parallel num_threads(nt) default(none) &
      !$omp shared(self, profile, dcdt, coarse, nf, site1, site2, phi, overlap, kcut, nblk, trafo_flag) &
      !$omp private(it, p, q, a, b, s, t, i, j, ib, ik0, nk, dk, work, acc, wv, wsum, kresp, terr)
      it = 1
      !$ it = omp_get_thread_num() + 1
      allocate (dk(self%npts, nf), work(self%npts), acc(nf, nf, self%npair), wv(vv_kblock, nf, nf), &
         & wsum(nf, nf))

      !$omp do schedule(dynamic) reduction(max:trafo_flag)
      do p = 1, self%npair
         a = site1(p)
         b = site2(p)
         do j = 1, nf
            work = dcdt(:, p)*profile(:, j)
            call self%slots(it)%trafo%fbt_r2k(work, dk(:, j), terr)
            if (allocated(terr)) then
               trafo_flag = 1
               deallocate (terr)
            end if
         end do
         acc = 0.0_wp
         wsum = 0.0_wp
         do ib = 1, nblk
            ik0 = (ib - 1)*vv_kblock
            nk = min(vv_kblock, kcut - ik0)
            wv = 0.0_wp
            do j = 1, nf
               do i = 1, nf
                  wv(1:nk, i, j) = phi(ik0 + 1:ik0 + nk, i)*dk(ik0 + 1:ik0 + nk, j)
                  wsum(i, j) = wsum(i, j) + sum(wv(:, i, j))
               end do
            end do
            do q = 1, self%npair
               s = site1(q)
               t = site2(q)
               if (a /= b) then
                  kresp = self%left_blk(:, s, a, ib)*self%right_blk(:, b, t, ib) + &
                     & self%left_blk(:, s, b, ib)*self%right_blk(:, a, t, ib)
               else
                  kresp = self%left_blk(:, s, a, ib)*self%right_blk(:, a, t, ib)
               end if
               do j = 1, nf
                  do i = 1, nf
                     acc(i, j, q) = acc(i, j, q) + sum(wv(:, i, j)*kresp)
                  end do
               end do
            end do
         end do
         do q = 1, self%npair
            do j = 1, nf
               do i = 1, nf
                  coarse((q - 1)*nf + i, (p - 1)*nf + j) = -acc(i, j, q)
               end do
            end do
         end do
         do j = 1, nf
            do i = 1, nf
               coarse((p - 1)*nf + i, (p - 1)*nf + j) = coarse((p - 1)*nf + i, (p - 1)*nf + j) + &
                  & overlap(i, j) + wsum(i, j)
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel

      if (trafo_flag /= 0) call fatal_error(error, "Radial transform failed in the 1D VV coarse matrix")
   end subroutine equation_coarse_matrix

end module moist_model_moz_solvent_vv_equation
