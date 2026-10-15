!> 1D UV equation on a uniform radial pair: the map c_s(r) -> t_s(r) and its derivative
!>
!> - Per k: h = W c P with the solute W(k) and the fixed solvent response
!>   P_ba = w_ba + rho_b h_ba of the VV solution
!> - Ng split: the map takes short-range c_s(r); c(k) = c_s(k) - beta u_lr(k)
!>   and t_s(k) = h(k) - c(k) - beta u_lr(k)
!>   The nonlinear driver iterates on t_s and supplies c_s through the closure
!> - The map is affine in c: dt = k2r(W dc(k) P - dc(k)) needs no linearization,
!>   and its transpose is u -> r2k^T(W^T k2r^T(u) P^T - k2r^T(u))
!> - The map depends on u_lr(k) through -beta k2r(W u_lr(k) P); the transpose
!>   optionally returns W^T k2r^T(u) P^T for that derivative
!> - UV columns (i - 1) nv + a hold solute site i and solvent site a
!>
!> Each parallel region uses one transform clone per thread and num_threads(nt)
module moist_model_moz_1d_equation
   !$ use omp_lib, only: omp_get_thread_num
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: resolve_num_threads
   use moist_math_packed_sym, only: packed_index_lut, npair_from_ns
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, transform_dst4
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo
   use moist_model_moz_1d_kernel, only: uv_kblock, uv_block_product, uv_block_product_transpose, &
      & uv_block_intramolecular_adjoint
   use moist_model_moz_fields, only: field_num_threads
   implicit none(type, external)
   private

   public :: uv_equation_type

   !> Transform directions of transform_columns
   integer, parameter :: to_k = 1, to_r = 2, to_k_adjoint = 3, to_r_adjoint = 4

   !> One radial transform owned by one thread
   type :: trafo_slot_type
      !> Transform with its own scratch
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
   end type trafo_slot_type

   !> Fixed UV data and the per-thread transforms
   type :: uv_equation_type
      !> OpenMP team size; zero takes omp_get_max_threads, small fields run serially
      integer :: nthreads = 0
      !> Solute sites, solvent sites, UV columns and grid points; 0 before setup
      integer, private :: nu = 0, nv = 0, nuv = 0, npts = 0
      !> Inverse thermal energy, 1/Hartree
      real(wp), private :: beta = 0.0_wp
      !> Solute W(k), packed (npts, nu(nu+1)/2)
      real(wp), allocatable, private :: wu(:, :)
      !> Solvent response P(k), (npts, nv, nv)
      real(wp), allocatable, private :: response(:, :, :)
      !> Long-range reciprocal potential, Hartree bohr^3, (npts, nu nv)
      real(wp), allocatable, private :: uk_lr(:, :)
      !> Packed-index lookup for the solute
      type(packed_index_lut), private :: lut
      !> Transform of the radial pair and its per-thread clones
      type(trafo_slot_type), allocatable, private :: slots(:)
      !> k-space scratch, (npts, nu nv)
      real(wp), allocatable, private :: ck_work(:, :), hk_work(:, :), tk_work(:, :)
   contains
      !> Bind the radial pair and the fixed data
      procedure :: setup => equation_setup
      !> Whether setup succeeded
      procedure :: ready => equation_ready
      !> Replace the reciprocal tail on the fixed radial pair, retaining transforms and scratch
      procedure :: update_potential => equation_update_potential
      !> Map c_s(r) -> t_s(r), optionally h(k)
      procedure :: map => equation_map
      !> Derivative dc(r) -> dt(r)
      procedure :: jvp => equation_jvp
      !> Transpose of the derivative
      procedure :: adjoint => equation_adjoint
      !> Grow the per-thread transform clones to nt
      procedure, private :: ensure_slots => equation_ensure_slots
      !> Transform all columns in one direction
      procedure, private :: transform_columns => equation_transform_columns
      !> Block products h = W c P or their transpose
      procedure, private :: products => equation_products
      !> Pull back h = W c P into the packed solute W
      procedure, private :: intramolecular_adjoint => equation_intramolecular_adjoint
      !> Check the shapes of a pair of fields
      procedure, private :: check_fields => equation_check_fields
   end type uv_equation_type

contains

   !> Bind a uniform DST-IV radial pair and the fixed UV data
   !>
   !> @param[inout] self      equation; a previous setup is dropped
   !> @param[in]    rgrid     real-space radial grid
   !> @param[in]    kgrid     reciprocal radial grid of the pair
   !> @param[in]    wu        solute W(k), packed (npts, nu(nu+1)/2)
   !> @param[in]    response  solvent response P(k) on kgrid, (npts, nv, nv)
   !> @param[in]    beta      inverse thermal energy, 1/Hartree
   !> @param[in]    uk_lr     long-range reciprocal potential, Hartree bohr^3, (npts, nu nv)
   !> @param[out]   error     not a uniform DST-IV pair, inconsistent shapes or transform failure
   subroutine equation_setup(self, rgrid, kgrid, wu, response, beta, uk_lr, error)
      !> Equation
      class(uv_equation_type), intent(inout) :: self
      !> Real-space and reciprocal radial grids
      type(moist_math_grid_radial_type), intent(in) :: rgrid, kgrid
      !> Solute W(k), solvent response and long-range potential
      real(wp), intent(in) :: wu(:, :), response(:, :, :), uk_lr(:, :)
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer :: npts, nu, nv

      self%nu = 0
      if (allocated(self%slots)) deallocate (self%slots)
      if (rgrid%transform /= transform_dst4 .or. kgrid%transform /= transform_dst4 .or. &
         & .not. allocated(rgrid%r) .or. .not. allocated(kgrid%r)) then
         call fatal_error(error, "1D UV equation needs the uniform DST-IV radial pair")
         return
      end if
      npts = size(rgrid%r)
      nv = size(response, 2)
      nu = 0
      do while (npair_from_ns(nu + 1) <= size(wu, 2))
         nu = nu + 1
      end do
      if (nu < 1 .or. nv < 1 .or. size(kgrid%r) /= npts .or. npair_from_ns(nu) /= size(wu, 2) .or. &
         & size(wu, 1) /= npts .or. any(shape(response) /= [npts, nv, nv]) .or. &
         & any(shape(uk_lr) /= [npts, nu*nv])) then
         call fatal_error(error, "1D UV equation data must have shapes (npts, nu(nu+1)/2), (npts, nv, nv) "// &
            & "and (npts, nu nv)")
         return
      end if
      allocate (self%slots(1))
      call new_radial_trafo(self%slots(1)%trafo, rgrid, kgrid, error)
      if (allocated(error)) then
         deallocate (self%slots)
         return
      end if
      call self%lut%init(nu)
      self%wu = wu
      self%response = response
      self%uk_lr = uk_lr
      self%beta = beta
      self%npts = npts
      self%nv = nv
      self%nuv = nu*nv
      if (allocated(self%ck_work)) deallocate (self%ck_work, self%hk_work, self%tk_work)
      allocate (self%ck_work(npts, nu*nv), self%hk_work(npts, nu*nv), self%tk_work(npts, nu*nv))
      self%nu = nu
   end subroutine equation_setup

   !> Whether setup succeeded
   !>
   !> @param[in] self Equation
   pure logical function equation_ready(self) result(ready)
      !> Equation
      class(uv_equation_type), intent(in) :: self
      ready = self%nu > 0
   end function equation_ready

   !> Replace the tail while retaining geometry, solvent response, beta and transforms
   !>
   !> @param[inout] self Prepared equation
   !> @param[in] uk_lr Long-range reciprocal potential, Hartree bohr^3, (npts, nu nv)
   !> @param[out] error Unprepared equation or wrong shape; old table retained on refusal
   subroutine equation_update_potential(self, uk_lr, error)
      class(uv_equation_type), intent(inout) :: self
      real(wp), intent(in) :: uk_lr(:, :)
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%ready()) then
         call fatal_error(error, "Update a set-up 1D UV equation's potential")
      else if (any(shape(uk_lr) /= [self%npts, self%nuv])) then
         call fatal_error(error, "1D UV long-range potential must have shape (npts, nu nv)")
      else
         self%uk_lr(:, :) = uk_lr
      end if
   end subroutine equation_update_potential

   !> Map the short-range direct correlation c_s(r) to t_s(r)
   !>
   !> @param[inout] self   equation after setup
   !> @param[in]    cr     c_s(r), (npts, nu nv)
   !> @param[out]   tr     t_s(r), (npts, nu nv)
   !> @param[out]   error  not set up, shape mismatch or transform failure
   !> @param[out]   hk     optional physical h(k), (npts, nu nv)
   subroutine equation_map(self, cr, tr, error, hk)
      !> Equation
      class(uv_equation_type), intent(inout) :: self
      !> Short-range direct correlation
      real(wp), intent(in) :: cr(:, :)
      !> Short-range indirect correlation
      real(wp), intent(out) :: tr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional total correlation in k-space
      real(wp), intent(out), optional :: hk(:, :)

      call self%check_fields(cr, tr, error)
      if (allocated(error)) return
      call self%transform_columns(to_k, cr, self%ck_work, error)
      if (allocated(error)) return
      ! Preserve c_s(k): t_s(k) = h(k) - c_s(k), without cancelling the Coulomb tail.
      self%tk_work = self%ck_work
      call shift(field_num_threads(self%npts*self%nuv, resolve_num_threads(self%nthreads)), &
         & self%ck_work, -self%beta, self%uk_lr)
      call self%products(self%ck_work, self%hk_work, .false.)
      self%tk_work = self%hk_work - self%tk_work
      call self%transform_columns(to_r, self%tk_work, tr, error)
      if (allocated(error)) return
      if (present(hk)) hk = self%hk_work
   end subroutine equation_map

   !> Derivative of the map: dt(r) = k2r(W dc(k) P - dc(k)), dc(k) = r2k(dc(r))
   !>
   !> @param[inout] self   equation after setup
   !> @param[in]    dcr    direction dc(r), (npts, nu nv)
   !> @param[out]   dtr    response dt(r), (npts, nu nv)
   !> @param[out]   error  not set up, shape mismatch or transform failure
   subroutine equation_jvp(self, dcr, dtr, error)
      !> Equation
      class(uv_equation_type), intent(inout) :: self
      !> Direction
      real(wp), intent(in) :: dcr(:, :)
      !> Response
      real(wp), intent(out) :: dtr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%check_fields(dcr, dtr, error)
      if (allocated(error)) return
      call self%transform_columns(to_k, dcr, self%ck_work, error)
      if (allocated(error)) return
      call self%products(self%ck_work, self%hk_work, .false.)
      call difference(field_num_threads(self%npts*self%nuv, resolve_num_threads(self%nthreads)), &
         & self%hk_work, self%ck_work, self%tk_work)
      call self%transform_columns(to_r, self%tk_work, dtr, error)
   end subroutine equation_jvp

   !> Transpose of the derivative in the Euclidean product
   !>
   !> dc(r) = r2k^T(W^T u(k) P^T - u(k)), u(k) = k2r^T(u(r)); the long-range
   !> part W^T u(k) P^T is the transpose of the map's u_lr(k) dependence up to -beta
   !>
   !> @param[inout] self   equation after setup
   !> @param[in]    dtr    direction u(r), (npts, nu nv)
   !> @param[out]   dcr    transpose product, (npts, nu nv)
   !> @param[out]   error  not set up, shape mismatch or transform failure
   !> @param[out]   wk     optional W^T u(k) P^T, (npts, nu nv)
   !> @param[in]    cr     optional primal c_s(r), required with ww
   !> @param[out]   ww     optional packed W adjoint, (npts, nu(nu+1)/2), required with cr
   subroutine equation_adjoint(self, dtr, dcr, error, wk, cr, ww)
      !> Equation
      class(uv_equation_type), intent(inout) :: self
      !> Direction
      real(wp), intent(in) :: dtr(:, :)
      !> Transpose product
      real(wp), intent(out) :: dcr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional long-range transpose
      real(wp), intent(out), optional :: wk(:, :)
      !> Primal short-range correlation and packed intramolecular adjoint
      real(wp), intent(in), optional :: cr(:, :)
      real(wp), intent(out), optional :: ww(:, :)

      call self%check_fields(dtr, dcr, error)
      if (allocated(error)) return
      if (present(cr) .neqv. present(ww)) then
         call fatal_error(error, "1D UV intramolecular adjoint needs both primal c_s and packed W adjoint")
         return
      end if
      if (present(ww)) then
         call self%check_fields(cr, dcr, error)
         if (allocated(error)) return
         if (any(shape(ww) /= shape(self%wu))) then
            call fatal_error(error, "1D UV intramolecular adjoint shape mismatch")
            return
         end if
      end if
      if (present(wk)) then
         call self%check_fields(dtr, wk, error)
         if (allocated(error)) return
      end if
      call self%transform_columns(to_r_adjoint, dtr, self%tk_work, error)
      if (allocated(error)) return
      if (present(ww)) then
         call self%transform_columns(to_k, cr, self%ck_work, error)
         if (allocated(error)) return
         call shift(field_num_threads(self%npts*self%nuv, resolve_num_threads(self%nthreads)), &
            & self%ck_work, -self%beta, self%uk_lr)
         call self%intramolecular_adjoint(self%ck_work, self%tk_work, ww)
      end if
      call self%products(self%tk_work, self%hk_work, .true.)
      if (present(wk)) wk = self%hk_work
      call difference(field_num_threads(self%npts*self%nuv, resolve_num_threads(self%nthreads)), &
         & self%hk_work, self%tk_work, self%ck_work)
      call self%transform_columns(to_k_adjoint, self%ck_work, dcr, error)
   end subroutine equation_adjoint

   !> Check that the equation is set up and both fields are (npts, nu nv)
   !>
   !> @param[in]  self   equation
   !> @param[in]  fin    input field
   !> @param[in]  fout   output field
   !> @param[out] error  not set up or shape mismatch
   subroutine equation_check_fields(self, fin, fout, error)
      !> Equation
      class(uv_equation_type), intent(in) :: self
      !> Input and output fields
      real(wp), intent(in) :: fin(:, :), fout(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%ready()) then
         call fatal_error(error, "1D UV equation is not set up")
      else if (any(shape(fin) /= [self%npts, self%nuv]) .or. any(shape(fout) /= shape(fin))) then
         call fatal_error(error, "1D UV fields must have shape (npts, nu nv)")
      end if
   end subroutine equation_check_fields

   !> Threaded in-place shift a = a + s u
   !>
   !> @param[in]    nt  OpenMP team size
   !> @param[inout] a   field (npts, nu nv)
   !> @param[in]    s   scale
   !> @param[in]    u   field of the same shape
   subroutine shift(nt, a, s, u)
      !> Team size
      integer, intent(in) :: nt
      !> Field
      real(wp), intent(inout) :: a(:, :)
      !> Scale
      real(wp), intent(in) :: s
      !> Shift field
      real(wp), intent(in) :: u(:, :)
      integer :: i, col
      !$omp parallel do collapse(2) schedule(static) num_threads(nt)
      do col = 1, size(a, 2)
         do i = 1, size(a, 1)
            a(i, col) = a(i, col) + s*u(i, col)
         end do
      end do
      !$omp end parallel do
   end subroutine shift

   !> Threaded difference out = a - b, plus s u when given
   !>
   !> @param[in]  nt   OpenMP team size
   !> @param[in]  a    field (npts, nu nv)
   !> @param[in]  b    field of the same shape
   !> @param[out] out  difference
   !> @param[in]  s    optional scale of u
   !> @param[in]  u    optional field added with scale s
   subroutine difference(nt, a, b, out, s, u)
      !> Team size
      integer, intent(in) :: nt
      !> Fields
      real(wp), intent(in) :: a(:, :), b(:, :)
      !> Difference
      real(wp), intent(out) :: out(:, :)
      !> Optional shift
      real(wp), intent(in), optional :: s
      real(wp), intent(in), optional :: u(:, :)
      integer :: i, col
      if (present(s) .and. present(u)) then
         !$omp parallel do collapse(2) schedule(static) num_threads(nt)
         do col = 1, size(a, 2)
            do i = 1, size(a, 1)
               out(i, col) = (a(i, col) - b(i, col)) + s*u(i, col)
            end do
         end do
         !$omp end parallel do
      else
         !$omp parallel do collapse(2) schedule(static) num_threads(nt)
         do col = 1, size(a, 2)
            do i = 1, size(a, 1)
               out(i, col) = a(i, col) - b(i, col)
            end do
         end do
         !$omp end parallel do
      end if
   end subroutine difference

   !> Grow the per-thread transform clones to at least nt
   !>
   !> @param[inout] self Equation after setup
   !> @param[in]    nt   team size
   subroutine equation_ensure_slots(self, nt)
      !> Equation
      class(uv_equation_type), intent(inout) :: self
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

   !> Transform all columns in one direction, one transform clone per thread
   !>
   !> @param[inout] self       equation after setup
   !> @param[in]    direction  to_k, to_r, to_k_adjoint or to_r_adjoint
   !> @param[in]    fin        input columns (npts, nu nv)
   !> @param[out]   fout       output columns (npts, nu nv)
   !> @param[out]   error      transform failure
   subroutine equation_transform_columns(self, direction, fin, fout, error)
      !> Equation
      class(uv_equation_type), intent(inout) :: self
      !> Transform direction
      integer, intent(in) :: direction
      !> Input columns
      real(wp), intent(in) :: fin(:, :)
      !> Output columns
      real(wp), intent(out) :: fout(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(error_type), allocatable :: terr
      integer :: nt, it, ip, trafo_flag

      nt = field_num_threads(self%npts*self%nuv, resolve_num_threads(self%nthreads))
      call self%ensure_slots(nt)
      trafo_flag = 0
      !$omp parallel num_threads(nt) default(none) shared(self, direction, fin, fout) &
      !$omp private(it, ip, terr) reduction(max:trafo_flag)
      it = 1
      !$ it = omp_get_thread_num() + 1
      !$omp do schedule(static)
      do ip = 1, self%nuv
         select case (direction)
         case (to_k)
            call self%slots(it)%trafo%fbt_r2k(fin(:, ip), fout(:, ip), terr)
         case (to_r)
            call self%slots(it)%trafo%fbt_k2r(fin(:, ip), fout(:, ip), terr)
         case (to_k_adjoint)
            call self%slots(it)%trafo%fbt_r2k_adj(fin(:, ip), fout(:, ip), terr)
         case default
            call self%slots(it)%trafo%fbt_k2r_adj(fin(:, ip), fout(:, ip), terr)
         end select
         if (allocated(terr)) then
            trafo_flag = 1
            deallocate (terr)
         end if
      end do
      !$omp end do
      !$omp end parallel
      if (trafo_flag /= 0) call fatal_error(error, "Radial transform failed in the 1D UV equation")
   end subroutine equation_transform_columns

   !> Block products g = W f P, or g = W^T f P^T when transposed, over k lanes
   !>
   !> @param[in]  self       equation after setup
   !> @param[in]  fk         k-space columns (npts, nu nv)
   !> @param[out] gk         k-space columns (npts, nu nv)
   !> @param[in]  transposed  apply the transpose
   subroutine equation_products(self, fk, gk, transposed)
      !> Equation
      class(uv_equation_type), intent(in) :: self
      !> Input and output columns
      real(wp), intent(in) :: fk(:, :)
      real(wp), intent(out) :: gk(:, :)
      !> Apply the transpose
      logical, intent(in) :: transposed
      real(wp), allocatable :: w_blk(:, :), f_blk(:, :), p_blk(:, :, :), g_blk(:, :)
      integer :: nt, ib, ik0, nk, nblock

      nt = field_num_threads(self%npts*self%nuv, resolve_num_threads(self%nthreads))
      nblock = (self%npts - 1)/uv_kblock + 1
      !$omp parallel num_threads(nt) default(none) shared(self, fk, gk, transposed, nblock) &
      !$omp private(ib, ik0, nk, w_blk, f_blk, p_blk, g_blk)
      allocate (w_blk(uv_kblock, size(self%wu, 2)), f_blk(uv_kblock, self%nuv), &
         & p_blk(uv_kblock, self%nv, self%nv), g_blk(uv_kblock, self%nuv))
      ! Lanes past the end of the grid stay zero
      w_blk = 0.0_wp
      f_blk = 0.0_wp
      p_blk = 0.0_wp
      !$omp do schedule(static)
      do ib = 1, nblock
         ik0 = (ib - 1)*uv_kblock
         nk = min(uv_kblock, self%npts - ik0)
         w_blk(1:nk, :) = self%wu(ik0 + 1:ik0 + nk, :)
         f_blk(1:nk, :) = fk(ik0 + 1:ik0 + nk, :)
         p_blk(1:nk, :, :) = self%response(ik0 + 1:ik0 + nk, :, :)
         if (transposed) then
            call uv_block_product_transpose(self%lut, self%nv, w_blk, f_blk, p_blk, g_blk)
         else
            call uv_block_product(self%lut, self%nv, w_blk, f_blk, p_blk, g_blk)
         end if
         gk(ik0 + 1:ik0 + nk, :) = g_blk(1:nk, :)
      end do
      !$omp end do
      !$omp end parallel
   end subroutine equation_products

   !> Pull back the physical h(k) into packed W(k), independently over k blocks
   !>
   !> @param[in]  self  equation after setup
   !> @param[in]  ck    physical direct correlation, (npts, nu nv)
   !> @param[in]  seed  adjoint of h(k), (npts, nu nv)
   !> @param[out] ww    adjoint of packed W(k), (npts, nu(nu+1)/2)
   subroutine equation_intramolecular_adjoint(self, ck, seed, ww)
      class(uv_equation_type), intent(in) :: self
      real(wp), intent(in) :: ck(:, :), seed(:, :)
      real(wp), intent(out) :: ww(:, :)
      real(wp), allocatable :: c_blk(:, :), p_blk(:, :, :), u_blk(:, :), ww_blk(:, :)
      integer :: nt, ib, ik0, nk, nblock

      nt = field_num_threads(self%npts*self%nuv, resolve_num_threads(self%nthreads))
      nblock = (self%npts - 1)/uv_kblock + 1
      !$omp parallel num_threads(nt) default(none) shared(self, ck, seed, ww, nblock) &
      !$omp private(ib, ik0, nk, c_blk, p_blk, u_blk, ww_blk)
      allocate (c_blk(uv_kblock, self%nuv), p_blk(uv_kblock, self%nv, self%nv), &
         & u_blk(uv_kblock, self%nuv), ww_blk(uv_kblock, size(self%wu, 2)))
      c_blk = 0.0_wp
      p_blk = 0.0_wp
      u_blk = 0.0_wp
      !$omp do schedule(static)
      do ib = 1, nblock
         ik0 = (ib - 1)*uv_kblock
         nk = min(uv_kblock, self%npts - ik0)
         c_blk(1:nk, :) = ck(ik0 + 1:ik0 + nk, :)
         p_blk(1:nk, :, :) = self%response(ik0 + 1:ik0 + nk, :, :)
         u_blk(1:nk, :) = seed(ik0 + 1:ik0 + nk, :)
         call uv_block_intramolecular_adjoint(self%lut, self%nv, c_blk, p_blk, u_blk, ww_blk)
         ww(ik0 + 1:ik0 + nk, :) = ww_blk(1:nk, :)
      end do
      !$omp end do
      !$omp end parallel
   end subroutine equation_intramolecular_adjoint

end module moist_model_moz_1d_equation
