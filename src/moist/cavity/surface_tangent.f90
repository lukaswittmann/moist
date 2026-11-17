!> Forward-mode dual of [[cavity_surface_adjoint_type]] plus callback-free
!> second-order host/stage-1 exchange data for cavity HVP and Hessians
!>
!> - One slot per observable channel the adjoint accumulator owns, holding the
!>   directional derivatives along a batch of `ndir` nuclear directions
!> - Curvature channels (`d_k1`, `d_k2`) allocated only when `want_curvature`
!>   is `.true.`, gated by `have_curvature`
!> - [[cavity_hessian_exchange_type]] carries pre-staged Stage-1 data
!>   (`host_jets`, optional pre-computed `tangent` + `d_wbranch`) into the
!>   cavity Hessian traversal and collects surface-motion / LSF-adjoint
!>   tangent outputs without callbacks
module moist_cavity_surface_tangent
   use mctc_env_accuracy, only: wp

   implicit none(type, external)
   private

   public :: cavity_surface_tangent_type
   public :: cavity_hessian_exchange_type

   !> Directional derivatives of the surface observables along a direction batch
   type :: cavity_surface_tangent_type
      !> Tangent of the point positions `dr_i/dv` (3, ngrid, ndir)
      real(wp), allocatable :: d_xyz(:, :, :)
      !> Tangent of the Gaussian widths `dxi_i/dv` (ngrid, ndir)
      real(wp), allocatable :: d_xi(:, :)
      !> Tangent of the outward unit normals `dn_i/dv` (3, ngrid, ndir)
      real(wp), allocatable :: d_n(:, :, :)
      !> Tangent of the surface areas `da_i/dv` (ngrid, ndir)
      real(wp), allocatable :: d_a(:, :)
      !> Tangent of the integration weights `dw_i/dv` (ngrid, ndir)
      real(wp), allocatable :: d_w(:, :)
      !> Tangent of the switching factors `df_i/dv` (ngrid, ndir)
      real(wp), allocatable :: d_f(:, :)
      !> Tangent of the first principal curvature `dk1_i/dv` (ngrid, ndir);
      !> allocated only when `have_curvature` is `.true.`
      real(wp), allocatable :: d_k1(:, :)
      !> Tangent of the second principal curvature `dk2_i/dv` (ngrid, ndir);
      !> allocated only when `have_curvature` is `.true.`
      real(wp), allocatable :: d_k2(:, :)
      !> Whether `d_k1` and `d_k2` are allocated and carry meaningful values
      logical :: have_curvature = .false.
   contains
      !> Allocate the channels for `(ngrid, ndir)` and initialize them to zero
      procedure :: init => init_surface_tangent
      !> Reset every allocated channel to zero
      procedure :: zero => zero_surface_tangent
      !> Deallocate every channel
      procedure :: destroy => destroy_surface_tangent
      !> Report whether every channel is allocated with consistent shapes
      procedure :: is_initialized => surface_tangent_is_initialized
      !> Number of directions of the batch, zero when uninitialized
      procedure :: ndir => surface_tangent_ndir
      !> Copy a contiguous direction slice `ilo:ihi` into `sub`
      procedure :: slice => slice_surface_tangent
   end type cavity_surface_tangent_type

   !> Callback-free second-order exchange container for cavity HVP / Hessians
   !>
   !> Carries Stage-1 host partial LSF jet tangents at fixed `xyz` (`host_jets`),
   !> optional pre-computed Stage-1 surface tangents (`tangent`, `d_wbranch`) so
   !> Stage-2 HVP avoids recomputing Pass 1, and optional output buffers for
   !> surface point motion (`d_xyz`) and gradient-path LSF adjoint weights +
   !> tangents (`w_*`, `dw_*`).
   type :: cavity_hessian_exchange_type
      !> Whether `host_jets` holds partial LSF jet tangents `(40, ngrid, ndir)`
      logical :: have_host_jets = .false.
      !> Partial LSF jet tangents at fixed `xyz`: `(40, ngrid, ndir)`
      real(wp), allocatable :: host_jets(:, :, :)

      !> Whether `tangent` and `d_wbranch` hold pre-computed Stage-1 tangents
      logical :: have_tangent = .false.
      !> Pre-computed Stage-1 surface tangent for `(ngrid, ndir)`
      type(cavity_surface_tangent_type) :: tangent
      !> Pre-computed Stage-1 softmax branch-weight tangent `(ngrid, ndir)`
      real(wp), allocatable :: d_wbranch(:, :)

      !> Whether output channels (`d_xyz`, plus `w_*`/`dw_*` when `want_lsf`) are active
      logical :: want_outputs = .false.
      !> Whether gradient-path LSF adjoint weights and their tangents are requested
      logical :: want_lsf = .false.
      !> Surface point motion `dr_i/dv` `(3, ngrid, ndir)`
      real(wp), allocatable :: d_xyz(:, :, :)
      !> Base gradient-path LSF adjoint weights: `(ngrid)`, `(3, ngrid)`, `(3, 3, ngrid)`
      real(wp), allocatable :: w_value(:), w_gradient(:, :), w_hessian(:, :, :)
      !> Tangents of LSF adjoint weights: `(ngrid, ndir)`, `(3, ngrid, ndir)`, `(3, 3, ngrid, ndir)`
      real(wp), allocatable :: dw_value(:, :), dw_gradient(:, :, :), dw_hessian(:, :, :, :)
   contains
      !> Allocate requested output channels to zero for `(ngrid, ndir)`
      procedure :: init_outputs => init_exchange_outputs
      !> Release output channels only
      procedure :: clear_outputs => clear_exchange_outputs
      !> Release all buffers (inputs and outputs)
      procedure :: destroy => destroy_hessian_exchange
      !> Whether the active output channels match `(ngrid, ndir)`
      procedure :: outputs_initialized => exchange_outputs_initialized
      !> Whether this exchange requires the per-direction fixed Hessian channel
      procedure :: requires_per_dir => exchange_requires_per_dir
   end type cavity_hessian_exchange_type

contains

   !> Allocate the surface-tangent channels and initialize them to zero
   !>
   !> @param[inout] self           Surface-tangent container
   !> @param[in]    ngrid          Number of surface grid points
   !> @param[in]    ndir           Number of nuclear directions
   !> @param[in]    want_curvature Allocate `d_k1` and `d_k2` as well (default `.false.`)
   subroutine init_surface_tangent(self, ngrid, ndir, want_curvature)
      !> Surface-tangent container
      class(cavity_surface_tangent_type), intent(inout) :: self
      !> Number of surface grid points
      integer, intent(in) :: ngrid
      !> Number of nuclear directions
      integer, intent(in) :: ndir
      !> Whether to allocate the curvature channels
      logical, intent(in), optional :: want_curvature

      !> Resolved curvature flag
      logical :: curv

      curv = .false.
      if (present(want_curvature)) curv = want_curvature

      call self%destroy()

      allocate (self%d_xyz(3, ngrid, ndir), source=0.0_wp)
      allocate (self%d_xi(ngrid, ndir), source=0.0_wp)
      allocate (self%d_n(3, ngrid, ndir), source=0.0_wp)
      allocate (self%d_a(ngrid, ndir), source=0.0_wp)
      allocate (self%d_w(ngrid, ndir), source=0.0_wp)
      allocate (self%d_f(ngrid, ndir), source=0.0_wp)
      allocate (self%d_k1(ngrid, ndir), source=0.0_wp)
      allocate (self%d_k2(ngrid, ndir), source=0.0_wp)
      self%have_curvature = curv

   end subroutine init_surface_tangent

   !> Reset every allocated channel to zero
   !>
   !> @param[inout] self  Surface-tangent container
   subroutine zero_surface_tangent(self)
      !> Surface-tangent container
      class(cavity_surface_tangent_type), intent(inout) :: self

      if (allocated(self%d_xyz)) self%d_xyz = 0.0_wp
      if (allocated(self%d_xi)) self%d_xi = 0.0_wp
      if (allocated(self%d_n)) self%d_n = 0.0_wp
      if (allocated(self%d_a)) self%d_a = 0.0_wp
      if (allocated(self%d_w)) self%d_w = 0.0_wp
      if (allocated(self%d_f)) self%d_f = 0.0_wp
      if (allocated(self%d_k1)) self%d_k1 = 0.0_wp
      if (allocated(self%d_k2)) self%d_k2 = 0.0_wp

   end subroutine zero_surface_tangent

   !> Deallocate every channel and clear `have_curvature`
   !>
   !> @param[inout] self  Surface-tangent container
   subroutine destroy_surface_tangent(self)
      !> Surface-tangent container
      class(cavity_surface_tangent_type), intent(inout) :: self

      if (allocated(self%d_xyz)) deallocate (self%d_xyz)
      if (allocated(self%d_xi)) deallocate (self%d_xi)
      if (allocated(self%d_n)) deallocate (self%d_n)
      if (allocated(self%d_a)) deallocate (self%d_a)
      if (allocated(self%d_w)) deallocate (self%d_w)
      if (allocated(self%d_f)) deallocate (self%d_f)
      if (allocated(self%d_k1)) deallocate (self%d_k1)
      if (allocated(self%d_k2)) deallocate (self%d_k2)
      self%have_curvature = .false.

   end subroutine destroy_surface_tangent

   !> Check whether the surface-tangent channels are allocated with consistent shapes
   !>
   !> @param[in] self  Surface-tangent container
   pure logical function surface_tangent_is_initialized(self) result(initialized)
      !> Surface-tangent container
      class(cavity_surface_tangent_type), intent(in) :: self

      !> Grid and direction extents read from `d_xi`
      integer :: ngrid, ndir

      initialized = allocated(self%d_xyz) .and. allocated(self%d_xi) .and. &
                    allocated(self%d_n) .and. allocated(self%d_a) .and. &
                    allocated(self%d_w) .and. allocated(self%d_f) .and. &
                    allocated(self%d_k1) .and. allocated(self%d_k2)
      if (.not. initialized) return

      ngrid = size(self%d_xi, 1)
      ndir = size(self%d_xi, 2)
      initialized = all(shape(self%d_xyz) == [3, ngrid, ndir]) .and. &
                    all(shape(self%d_n) == [3, ngrid, ndir]) .and. &
                    all(shape(self%d_a) == [ngrid, ndir]) .and. &
                    all(shape(self%d_w) == [ngrid, ndir]) .and. &
                    all(shape(self%d_f) == [ngrid, ndir]) .and. &
                    all(shape(self%d_k1) == [ngrid, ndir]) .and. &
                    all(shape(self%d_k2) == [ngrid, ndir])

   end function surface_tangent_is_initialized

   !> Number of directions of the batch, zero when uninitialized
   pure integer function surface_tangent_ndir(self) result(ndir)
      class(cavity_surface_tangent_type), intent(in) :: self

      ndir = 0
      if (allocated(self%d_xi)) ndir = size(self%d_xi, 2)

   end function surface_tangent_ndir

   !> Copy direction slice `ilo:ihi` into `sub`
   !>
   !> @param[in]  self            Full-batch surface tangent
   !> @param[in]  ilo             First direction index (1-based)
   !> @param[in]  ihi             Last direction index (inclusive)
   !> @param[in]  want_curvature  Whether `sub` should carry curvature channels
   !> @param[out] sub             Sliced surface tangent for `nblk = ihi - ilo + 1`
   subroutine slice_surface_tangent(self, ilo, ihi, want_curvature, sub)
      !> Full-batch surface tangent
      class(cavity_surface_tangent_type), intent(in) :: self
      !> Direction slice bounds
      integer, intent(in) :: ilo, ihi
      !> Curvature channel request
      logical, intent(in) :: want_curvature
      !> Sliced surface tangent
      type(cavity_surface_tangent_type), intent(out) :: sub

      integer :: ngrid, nblk

      ngrid = size(self%d_xi, 1)
      nblk = ihi - ilo + 1
      call sub%init(ngrid, nblk, want_curvature)
      if (ngrid <= 0 .or. nblk <= 0) return

      sub%d_xyz = self%d_xyz(:, :, ilo:ihi)
      sub%d_xi = self%d_xi(:, ilo:ihi)
      sub%d_n = self%d_n(:, :, ilo:ihi)
      sub%d_a = self%d_a(:, ilo:ihi)
      sub%d_w = self%d_w(:, ilo:ihi)
      sub%d_f = self%d_f(:, ilo:ihi)
      if (want_curvature .and. self%have_curvature) then
         sub%d_k1 = self%d_k1(:, ilo:ihi)
         sub%d_k2 = self%d_k2(:, ilo:ihi)
      end if
   end subroutine slice_surface_tangent

   !> Allocate requested output channels of `cavity_hessian_exchange_type` to zero
   !>
   !> @param[inout] self      Exchange container
   !> @param[in]    ngrid     Number of surface grid points
   !> @param[in]    ndir      Number of directions in the batch
   !> @param[in]    want_lsf  Whether gradient-path LSF adjoint weights & tangents are requested
   subroutine init_exchange_outputs(self, ngrid, ndir, want_lsf)
      class(cavity_hessian_exchange_type), intent(inout) :: self
      integer, intent(in) :: ngrid, ndir
      logical, intent(in) :: want_lsf

      call self%clear_outputs()
      self%want_outputs = .true.
      self%want_lsf = want_lsf
      allocate (self%d_xyz(3, ngrid, ndir), source=0.0_wp)
      if (want_lsf) then
         allocate (self%w_value(ngrid), source=0.0_wp)
         allocate (self%w_gradient(3, ngrid), source=0.0_wp)
         allocate (self%w_hessian(3, 3, ngrid), source=0.0_wp)
         allocate (self%dw_value(ngrid, ndir), source=0.0_wp)
         allocate (self%dw_gradient(3, ngrid, ndir), source=0.0_wp)
         allocate (self%dw_hessian(3, 3, ngrid, ndir), source=0.0_wp)
      end if
   end subroutine init_exchange_outputs

   !> Release output channels of `cavity_hessian_exchange_type`
   subroutine clear_exchange_outputs(self)
      class(cavity_hessian_exchange_type), intent(inout) :: self

      self%want_outputs = .false.
      self%want_lsf = .false.
      if (allocated(self%d_xyz)) deallocate (self%d_xyz)
      if (allocated(self%w_value)) deallocate (self%w_value)
      if (allocated(self%w_gradient)) deallocate (self%w_gradient)
      if (allocated(self%w_hessian)) deallocate (self%w_hessian)
      if (allocated(self%dw_value)) deallocate (self%dw_value)
      if (allocated(self%dw_gradient)) deallocate (self%dw_gradient)
      if (allocated(self%dw_hessian)) deallocate (self%dw_hessian)
   end subroutine clear_exchange_outputs

   !> Release all inputs and outputs of `cavity_hessian_exchange_type`
   subroutine destroy_hessian_exchange(self)
      class(cavity_hessian_exchange_type), intent(inout) :: self

      self%have_host_jets = .false.
      if (allocated(self%host_jets)) deallocate (self%host_jets)
      self%have_tangent = .false.
      call self%tangent%destroy()
      if (allocated(self%d_wbranch)) deallocate (self%d_wbranch)
      call self%clear_outputs()
   end subroutine destroy_hessian_exchange

   !> Check whether active output channels match `(ngrid, ndir)`
   pure logical function exchange_outputs_initialized(self, ngrid, ndir) result(ok)
      class(cavity_hessian_exchange_type), intent(in) :: self
      integer, intent(in) :: ngrid, ndir

      if (.not. self%want_outputs) then
         ok = .true.
         return
      end if
      ok = allocated(self%d_xyz)
      if (.not. ok) return
      ok = all(shape(self%d_xyz) == [3, ngrid, ndir])
      if (.not. ok) return
      if (self%want_lsf) then
         ok = allocated(self%w_value) .and. allocated(self%w_gradient) .and. &
              allocated(self%w_hessian) .and. allocated(self%dw_value) .and. &
              allocated(self%dw_gradient) .and. allocated(self%dw_hessian)
         if (.not. ok) return
         ok = size(self%w_value) == ngrid .and. &
              all(shape(self%w_gradient) == [3, ngrid]) .and. &
              all(shape(self%w_hessian) == [3, 3, ngrid]) .and. &
              all(shape(self%dw_value) == [ngrid, ndir]) .and. &
              all(shape(self%dw_gradient) == [3, ngrid, ndir]) .and. &
              all(shape(self%dw_hessian) == [3, 3, ngrid, ndir])
      end if
   end function exchange_outputs_initialized

   !> Whether this exchange forces the per-direction fixed Hessian channel
   pure logical function exchange_requires_per_dir(self) result(per_dir)
      class(cavity_hessian_exchange_type), intent(in) :: self

      per_dir = self%have_host_jets .or. self%want_lsf
   end function exchange_requires_per_dir

end module moist_cavity_surface_tangent
