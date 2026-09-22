!> Raw host calculations, per-output requirements, and scoped component reads
module moist_channels_request
   use, intrinsic :: iso_c_binding, only: c_int64_t
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   implicit none(type, external)
   private
   public :: coupling_type, coupling_view_type, coupling_registry_type
   public :: coupling_request_type, point_potential_request_type, gaussian_potential_request_type
   public :: gaussian_moment_request_type
   public :: moist_phase_none, moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, moist_n_phases
   public :: output_phi, output_dphi_dr, output_dphi_dxi, output_gt, output_pt, output_mt, output_rt
   public :: grid_xyz, grid_normal, grid_xi, grid_f, grid_area, request_name_len


   !> No phase is staged
   integer, parameter :: moist_phase_none = 0
   !> Energy phase
   integer, parameter :: moist_phase_energy = 1
   !> Electronic response phase
   integer, parameter :: moist_phase_response = 2
   !> Nuclear gradient phase
   integer, parameter :: moist_phase_gradient = 3
   !> Supported phases
   integer, parameter :: moist_n_phases = 3
   !> Diagnostic name capacity
   integer, parameter :: request_name_len = 32
   !> Output bit for phi
   integer, parameter :: output_phi = 1
   !> Output bit for dphi_dr
   integer, parameter :: output_dphi_dr = 2
   !> Output bit for dphi_dxi
   integer, parameter :: output_dphi_dxi = 4
   !> Output bit for gt
   integer, parameter :: output_gt = 8
   !> Output bit for pt
   integer, parameter :: output_pt = 16
   !> Output bit for mt
   integer, parameter :: output_mt = 32
   !> Output bit for rt
   integer, parameter :: output_rt = 64
   !> Grid field bit for xyz
   integer, parameter :: grid_xyz = 1
   !> Grid field bit for normal
   integer, parameter :: grid_normal = 2
   !> Grid field bit for xi
   integer, parameter :: grid_xi = 4
   !> Grid field bit for f
   integer, parameter :: grid_f = 8
   !> Grid field bit for area
   integer, parameter :: grid_area = 16

   !> Scientific request base; all answers are independent of MOIST charges
   type, abstract :: coupling_request_type
      !> Opaque host token, replaced on every staging or invalidation
      integer(c_int64_t) :: handle = 0_c_int64_t
      !> Required output bits for each phase
      integer :: required(moist_n_phases) = 0
      !> Valid output bits
      integer :: available = 0
      !> Output bits touched by the latest answer operation
      integer :: touched = 0
      !> Shared grid fields used by this calculation
      integer :: fields = grid_xyz
      !> Expected grid size
      integer :: ngrid = -1
      !> Missing output bits for the staged phase
      integer :: missing = 0
      !> Last rejected answer, cleared by a successful retry
      type(error_type), allocatable :: failure
   contains
      procedure(request_name), deferred :: name
      procedure :: invalidate => request_invalidate
      procedure :: same_inputs => request_same_inputs
      procedure(request_merge), deferred :: merge
   end type coupling_request_type

   !> Polymorphic scientific request operations
   abstract interface

      !> Return a scientific calculation name
      !>
      !> @param[in] self Request to describe
      function request_name(self) result(name)
         import :: coupling_request_type, request_name_len
         implicit none(type, external)
         !> Request to describe
         class(coupling_request_type), intent(in) :: self
         !> Diagnostic name
         character(len=request_name_len) :: name
      end function request_name

      !> Merge only the outputs touched by this answer
      !>
      !> @param[in,out] self Live scientific request
      !> @param[in] answer Independent host answer value
      subroutine request_merge(self, answer)
         import :: coupling_request_type
         implicit none(type, external)
         !> Live scientific request
         class(coupling_request_type), intent(inout) :: self
         !> Independent host answer value
         class(coupling_request_type), intent(in) :: answer
      end subroutine request_merge
   end interface

   !> Raw point_potential calculation and its independently valid outputs
   type, extends(coupling_request_type) :: point_potential_request_type
      !> Raw phi answer
      real(wp), allocatable, private :: phi(:)
      !> Raw dphi_dr answer
      real(wp), allocatable, private :: dphi_dr(:, :)
   contains
      procedure :: merge => point_potential_merge
      procedure :: name => point_potential_name
      procedure :: require => point_potential_require
      procedure :: set => point_potential_set
      procedure :: get => point_potential_get
   end type point_potential_request_type

   !> Raw gaussian_potential calculation and its independently valid outputs
   type, extends(coupling_request_type) :: gaussian_potential_request_type
      !> Raw phi answer
      real(wp), allocatable, private :: phi(:)
      !> Raw dphi_dr answer
      real(wp), allocatable, private :: dphi_dr(:, :)
      !> Raw dphi_dxi answer
      real(wp), allocatable, private :: dphi_dxi(:)
   contains
      procedure :: merge => gaussian_potential_merge
      procedure :: name => gaussian_potential_name
      procedure :: require => gaussian_potential_require
      procedure :: set => gaussian_potential_set
      procedure :: get => gaussian_potential_get
   end type gaussian_potential_request_type

   !> Raw gaussian_moments calculation and its independently valid outputs
   type, extends(coupling_request_type) :: gaussian_moment_request_type
      !> Moment exponents in bohr**(-2), supplied by the component
      real(wp), allocatable :: width(:)
      !> Raw gt answer
      real(wp), allocatable, private :: gt(:)
      !> Raw pt answer
      real(wp), allocatable, private :: pt(:, :)
      !> Raw mt answer
      real(wp), allocatable, private :: mt(:, :, :)
      !> Raw rt answer
      real(wp), allocatable, private :: rt(:, :)
   contains
      procedure :: merge => gaussian_moment_merge
      procedure :: same_inputs => moment_same_inputs
      procedure :: name => gaussian_moment_name
      procedure :: require => gaussian_moment_require
      procedure :: set => gaussian_moment_set
      procedure :: get => gaussian_moment_get
   end type gaussian_moment_request_type

   !> Owning slot for a heterogeneous calculation
   type :: request_slot
      !> Concrete calculation
      class(coupling_request_type), allocatable :: item
   end type request_slot

   !> Registration local to a cavity or component
   type :: registration_type
      !> Component number; zero denotes the cavity
      integer :: scope = 0
      !> Local scientific name chosen during registration
      character(len=request_name_len) :: key = ""
      !> Registered request slot
      integer :: slot = 0
   end type registration_type
   !> Individually optional grid inputs
   type :: grid_type
      !> Evaluation points in bohr
      real(wp), allocatable :: xyz(:, :)
      !> Outward normals
      real(wp), allocatable :: normal(:, :)
      !> Gaussian inverse lengths in bohr**(-1)
      real(wp), allocatable :: xi(:)
      !> Switching factors
      real(wp), allocatable :: f(:)
      !> Surface areas in bohr**2
      real(wp), allocatable :: area(:)
   end type grid_type

   !> Model-owned collection exposed to hosts through borrowed handles
   type :: coupling_type
      private
      !> Scientific calculations
      type(request_slot), allocatable :: requests(:)
      !> Component-local registrations
      type(registration_type), allocatable :: registrations(:)
      !> Next never-reused handle
      integer(c_int64_t) :: serial = 0_c_int64_t
      !> Epoch of scoped views
      integer(c_int64_t) :: epoch = 0_c_int64_t
      !> Scope currently registering requirements
      integer :: scope = 1
      !> Frozen pending slot indices
      integer, allocatable :: pending_index(:)
      !> Grid inputs
      type(grid_type) :: grid
      !> Whether the declared grid snapshot is current
      logical :: grid_valid = .false.
      !> Current phase
      integer, public :: phase = moist_phase_none
      !> Number of evaluation points
      integer, public :: ngrid = 0
   contains
      procedure :: register => coupling_register
      procedure :: begin_registration
      procedure :: set_scope
      procedure :: snapshot => coupling_snapshot
      procedure :: arm => coupling_arm
      procedure :: invalidate => coupling_invalidate
      procedure :: n_requests
      procedure :: n_pending
      procedure :: declared
      procedure :: pending
      procedure :: lookup
      procedure :: answer
      procedure :: check_mandatory
      procedure :: make_view
      procedure :: close_view
      procedure, private :: compact
      procedure :: read_grid
      procedure :: grid_fields
      procedure :: has_snapshot
   end type coupling_type

   !> Short-lived read-only component view; never exposes pointers to requests
   type :: coupling_view_type
      private
      !> Borrowed collection, valid only during the component call
      type(coupling_type), pointer :: owner => null()
      !> Component-local registration scope
      integer :: scope = 0
      !> Epoch at construction, checked on every read
      integer(c_int64_t) :: epoch = -1_c_int64_t
      !> Current phase
      integer, public :: phase = moist_phase_none
   contains
      procedure :: read => view_read
      procedure :: check_mandatory => view_check
      procedure :: read_potential
      procedure :: read_moments
      procedure, private :: find_slot => view_find_slot
   end type coupling_view_type

   !> Model-owned node; hosts only borrow its collection
   type :: coupling_node
      !> Stable collection allocation
      type(coupling_type), allocatable :: item
      !> Next owned node
      type(coupling_node), allocatable :: next
   end type coupling_node

   !> Registry of every coupling minted by one model
   type :: coupling_registry_type
      private
      !> First owned node
      type(coupling_node), allocatable :: first
   contains
      procedure :: mint
      procedure :: release => registry_release
      procedure :: invalidate => registry_invalidate
      procedure :: owns
      procedure :: clear => registry_clear
   end type coupling_registry_type
contains

   !> Invalidate all outputs without changing their allocated storage
   !>
   !> @param[in,out] self Request to invalidate
   subroutine request_invalidate(self)
      !> Request to invalidate
      class(coupling_request_type), intent(inout) :: self
      self%available = 0
      self%touched = 0
      if (allocated(self%failure)) deallocate (self%failure)
   end subroutine request_invalidate

   !> Name of the point_potential calculation
   !>
   !> @param[in] self Request to describe
   function point_potential_name(self) result(name)
      !> Request to describe
      class(point_potential_request_type), intent(in) :: self
      !> Diagnostic name
      character(len=request_name_len) :: name
      name = "point_potential"
   end function point_potential_name

   !> Declare individual scientific outputs for one phase
   !>
   !> @param[in,out] self Request requirements
   !> @param[in] phase Phase index
   !> @param[in] phi Whether this output is required
   !> @param[in] dphi_dr Whether this output is required
   subroutine point_potential_require(self, phase, phi, dphi_dr)
      !> Request requirements
      class(point_potential_request_type), intent(inout) :: self
      !> Phase index
      integer, intent(in) :: phase
      !> Require phi
      logical, intent(in), optional :: phi
      !> Require dphi_dr
      logical, intent(in), optional :: dphi_dr
      if (phase < 1 .or. phase > moist_n_phases) return
      if (present(phi)) then
         if (phi) self%required(phase) = ior(self%required(phase), output_phi)
      end if
      if (present(dphi_dr)) then
         if (dphi_dr) self%required(phase) = ior(self%required(phase), output_dphi_dr)
      end if
   end subroutine point_potential_require

   !> Store selected raw outputs with shape and validity checks
   !>
   !> @param[in,out] self Scientific request
   !> @param[in] phi Raw phi array on the evaluation grid
   !> @param[in] dphi_dr Raw dphi_dr array on the evaluation grid
   subroutine point_potential_set(self, phi, dphi_dr)
      !> Scientific request
      class(point_potential_request_type), intent(inout) :: self
      !> Raw phi array
      real(wp), intent(in), optional :: phi(:)
      !> Raw dphi_dr array
      real(wp), intent(in), optional :: dphi_dr(:, :)
      if (present(phi)) then
         self%touched = ior(self%touched, output_phi)
         self%available = iand(self%available, not(output_phi))
         if (any(shape(phi) /= [self%ngrid])) then
            call fatal_error(self%failure, "point_potential: phi shape mismatch")
         else if (.not. all(ieee_is_finite(phi))) then
            call fatal_error(self%failure, "point_potential: phi must be finite")
         else
            self%phi = phi
            self%available = ior(self%available, output_phi)
         end if
      end if
      if (present(dphi_dr)) then
         self%touched = ior(self%touched, output_dphi_dr)
         self%available = iand(self%available, not(output_dphi_dr))
         if (any(shape(dphi_dr) /= [3, self%ngrid])) then
            call fatal_error(self%failure, "point_potential: dphi_dr shape mismatch")
         else if (.not. all(ieee_is_finite(dphi_dr))) then
            call fatal_error(self%failure, "point_potential: dphi_dr must be finite")
         else
            self%dphi_dr = dphi_dr
            self%available = ior(self%available, output_dphi_dr)
         end if
      end if
      if (iand(self%touched, not(self%available)) == 0) then
         if (allocated(self%failure)) deallocate (self%failure)
      end if
   end subroutine point_potential_set
   
   !> Read selected raw outputs with shape and validity checks
   !>
   !> @param[in] self Scientific request
   !> @param[out] phi Raw phi array on the evaluation grid
   !> @param[out] dphi_dr Raw dphi_dr array on the evaluation grid
   !> @param[out] error Missing or invalid output
   subroutine point_potential_get(self, phi, dphi_dr, error)
      !> Scientific request
      class(point_potential_request_type), intent(in) :: self
      !> Raw phi array
      real(wp), allocatable, intent(out), optional :: phi(:)
      !> Raw dphi_dr array
      real(wp), allocatable, intent(out), optional :: dphi_dr(:, :)
      !> Missing or invalid output
      type(error_type), allocatable, intent(out) :: error
      if (present(phi)) then
         if (iand(self%available, output_phi) == 0) then
            call fatal_error(error, "point_potential: missing required output phi")
            return
         end if
         phi = self%phi
      end if
      if (present(dphi_dr)) then
         if (iand(self%available, output_dphi_dr) == 0) then
            call fatal_error(error, "point_potential: missing required output dphi_dr")
            return
         end if
         dphi_dr = self%dphi_dr
      end if
   end subroutine point_potential_get

   !> Name of the gaussian_potential calculation
   !>
   !> @param[in] self Request to describe
   function gaussian_potential_name(self) result(name)
      !> Request to describe
      class(gaussian_potential_request_type), intent(in) :: self
      !> Diagnostic name
      character(len=request_name_len) :: name
      name = "gaussian_potential"
   end function gaussian_potential_name

   !> Declare individual scientific outputs for one phase
   !>
   !> @param[in,out] self Request requirements
   !> @param[in] phase Phase index
   !> @param[in] phi Whether this output is required
   !> @param[in] dphi_dr Whether this output is required
   !> @param[in] dphi_dxi Whether this output is required
   subroutine gaussian_potential_require(self, phase, phi, dphi_dr, dphi_dxi)
      !> Request requirements
      class(gaussian_potential_request_type), intent(inout) :: self
      !> Phase index
      integer, intent(in) :: phase
      !> Require phi
      logical, intent(in), optional :: phi
      !> Require dphi_dr
      logical, intent(in), optional :: dphi_dr
      !> Require dphi_dxi
      logical, intent(in), optional :: dphi_dxi
      if (phase < 1 .or. phase > moist_n_phases) return
      self%fields = ior(grid_xyz, grid_xi)
      if (present(phi)) then
         if (phi) self%required(phase) = ior(self%required(phase), output_phi)
      end if
      if (present(dphi_dr)) then
         if (dphi_dr) self%required(phase) = ior(self%required(phase), output_dphi_dr)
      end if
      if (present(dphi_dxi)) then
         if (dphi_dxi) self%required(phase) = ior(self%required(phase), output_dphi_dxi)
      end if
   end subroutine gaussian_potential_require

   !> Store selected raw outputs with shape and validity checks
   !>
   !> @param[in,out] self Scientific request
   !> @param[in] phi Raw phi array on the evaluation grid
   !> @param[in] dphi_dr Raw dphi_dr array on the evaluation grid
   !> @param[in] dphi_dxi Raw dphi_dxi array on the evaluation grid
   subroutine gaussian_potential_set(self, phi, dphi_dr, dphi_dxi)
      !> Scientific request
      class(gaussian_potential_request_type), intent(inout) :: self
      !> Raw phi array
      real(wp), intent(in), optional :: phi(:)
      !> Raw dphi_dr array
      real(wp), intent(in), optional :: dphi_dr(:, :)
      !> Raw dphi_dxi array
      real(wp), intent(in), optional :: dphi_dxi(:)
      if (present(phi)) then
         self%touched = ior(self%touched, output_phi)
         self%available = iand(self%available, not(output_phi))
         if (any(shape(phi) /= [self%ngrid])) then
            call fatal_error(self%failure, "gaussian_potential: phi shape mismatch")
         else if (.not. all(ieee_is_finite(phi))) then
            call fatal_error(self%failure, "gaussian_potential: phi must be finite")
         else
            self%phi = phi
            self%available = ior(self%available, output_phi)
         end if
      end if
      if (present(dphi_dr)) then
         self%touched = ior(self%touched, output_dphi_dr)
         self%available = iand(self%available, not(output_dphi_dr))
         if (any(shape(dphi_dr) /= [3, self%ngrid])) then
            call fatal_error(self%failure, "gaussian_potential: dphi_dr shape mismatch")
         else if (.not. all(ieee_is_finite(dphi_dr))) then
            call fatal_error(self%failure, "gaussian_potential: dphi_dr must be finite")
         else
            self%dphi_dr = dphi_dr
            self%available = ior(self%available, output_dphi_dr)
         end if
      end if
      if (present(dphi_dxi)) then
         self%touched = ior(self%touched, output_dphi_dxi)
         self%available = iand(self%available, not(output_dphi_dxi))
         if (any(shape(dphi_dxi) /= [self%ngrid])) then
            call fatal_error(self%failure, "gaussian_potential: dphi_dxi shape mismatch")
         else if (.not. all(ieee_is_finite(dphi_dxi))) then
            call fatal_error(self%failure, "gaussian_potential: dphi_dxi must be finite")
         else
            self%dphi_dxi = dphi_dxi
            self%available = ior(self%available, output_dphi_dxi)
         end if
      end if
      if (iand(self%touched, not(self%available)) == 0) then
         if (allocated(self%failure)) deallocate (self%failure)
      end if
   end subroutine gaussian_potential_set

   !> Read selected raw outputs with shape and validity checks
   !>
   !> @param[in] self Scientific request
   !> @param[out] phi Raw phi array on the evaluation grid
   !> @param[out] dphi_dr Raw dphi_dr array on the evaluation grid
   !> @param[out] dphi_dxi Raw dphi_dxi array on the evaluation grid
   !> @param[out] error Missing or invalid output
   subroutine gaussian_potential_get(self, phi, dphi_dr, dphi_dxi, error)
      !> Scientific request
      class(gaussian_potential_request_type), intent(in) :: self
      !> Raw phi array
      real(wp), allocatable, intent(out), optional :: phi(:)
      !> Raw dphi_dr array
      real(wp), allocatable, intent(out), optional :: dphi_dr(:, :)
      !> Raw dphi_dxi array
      real(wp), allocatable, intent(out), optional :: dphi_dxi(:)
      !> Missing or invalid output
      type(error_type), allocatable, intent(out) :: error
      if (present(phi)) then
         if (iand(self%available, output_phi) == 0) then
            call fatal_error(error, "gaussian_potential: missing required output phi")
            return
         end if
         phi = self%phi
      end if
      if (present(dphi_dr)) then
         if (iand(self%available, output_dphi_dr) == 0) then
            call fatal_error(error, "gaussian_potential: missing required output dphi_dr")
            return
         end if
         dphi_dr = self%dphi_dr
      end if
      if (present(dphi_dxi)) then
         if (iand(self%available, output_dphi_dxi) == 0) then
            call fatal_error(error, "gaussian_potential: missing required output dphi_dxi")
            return
         end if
         dphi_dxi = self%dphi_dxi
      end if
   end subroutine gaussian_potential_get

   !> Name of the gaussian_moments calculation
   !>
   !> @param[in] self Request to describe
   function gaussian_moment_name(self) result(name)
      !> Request to describe
      class(gaussian_moment_request_type), intent(in) :: self
      !> Diagnostic name
      character(len=request_name_len) :: name
      name = "gaussian_moments"
   end function gaussian_moment_name

   !> Declare individual scientific outputs for one phase
   !>
   !> @param[in,out] self Request requirements
   !> @param[in] phase Phase index
   !> @param[in] gt Whether this output is required
   !> @param[in] pt Whether this output is required
   !> @param[in] mt Whether this output is required
   !> @param[in] rt Whether this output is required
   subroutine gaussian_moment_require(self, phase, gt, pt, mt, rt)
      !> Request requirements
      class(gaussian_moment_request_type), intent(inout) :: self
      !> Phase index
      integer, intent(in) :: phase
      !> Require gt
      logical, intent(in), optional :: gt
      !> Require pt
      logical, intent(in), optional :: pt
      !> Require mt
      logical, intent(in), optional :: mt
      !> Require rt
      logical, intent(in), optional :: rt
      if (phase < 1 .or. phase > moist_n_phases) return
      if (present(gt)) then
         if (gt) self%required(phase) = ior(self%required(phase), output_gt)
      end if
      if (present(pt)) then
         if (pt) self%required(phase) = ior(self%required(phase), output_pt)
      end if
      if (present(mt)) then
         if (mt) self%required(phase) = ior(self%required(phase), output_mt)
      end if
      if (present(rt)) then
         if (rt) self%required(phase) = ior(self%required(phase), output_rt)
      end if
   end subroutine gaussian_moment_require

   !> Store selected raw outputs with shape and validity checks
   !>
   !> @param[in,out] self Scientific request
   !> @param[in] gt Raw gt array on the evaluation grid
   !> @param[in] pt Raw pt array on the evaluation grid
   !> @param[in] mt Raw mt array on the evaluation grid
   !> @param[in] rt Raw rt array on the evaluation grid
   subroutine gaussian_moment_set(self, gt, pt, mt, rt)
      !> Scientific request
      class(gaussian_moment_request_type), intent(inout) :: self
      !> Raw gt array
      real(wp), intent(in), optional :: gt(:)
      !> Raw pt array
      real(wp), intent(in), optional :: pt(:, :)
      !> Raw mt array
      real(wp), intent(in), optional :: mt(:, :, :)
      !> Raw rt array
      real(wp), intent(in), optional :: rt(:, :)
      if (present(gt)) then
         self%touched = ior(self%touched, output_gt)
         self%available = iand(self%available, not(output_gt))
         if (any(shape(gt) /= [self%ngrid])) then
            call fatal_error(self%failure, "gaussian_moments: gt shape mismatch")
         else if (.not. all(ieee_is_finite(gt))) then
            call fatal_error(self%failure, "gaussian_moments: gt must be finite")
         else
            self%gt = gt
            self%available = ior(self%available, output_gt)
         end if
      end if
      if (present(pt)) then
         self%touched = ior(self%touched, output_pt)
         self%available = iand(self%available, not(output_pt))
         if (any(shape(pt) /= [3, self%ngrid])) then
            call fatal_error(self%failure, "gaussian_moments: pt shape mismatch")
         else if (.not. all(ieee_is_finite(pt))) then
            call fatal_error(self%failure, "gaussian_moments: pt must be finite")
         else
            self%pt = pt
            self%available = ior(self%available, output_pt)
         end if
      end if
      if (present(mt)) then
         self%touched = ior(self%touched, output_mt)
         self%available = iand(self%available, not(output_mt))
         if (any(shape(mt) /= [3, 3, self%ngrid])) then
            call fatal_error(self%failure, "gaussian_moments: mt shape mismatch")
         else if (.not. all(ieee_is_finite(mt))) then
            call fatal_error(self%failure, "gaussian_moments: mt must be finite")
         else
            self%mt = mt
            self%available = ior(self%available, output_mt)
         end if
      end if
      if (present(rt)) then
         self%touched = ior(self%touched, output_rt)
         self%available = iand(self%available, not(output_rt))
         if (any(shape(rt) /= [3, self%ngrid])) then
            call fatal_error(self%failure, "gaussian_moments: rt shape mismatch")
         else if (.not. all(ieee_is_finite(rt))) then
            call fatal_error(self%failure, "gaussian_moments: rt must be finite")
         else
            self%rt = rt
            self%available = ior(self%available, output_rt)
         end if
      end if
      if (iand(self%touched, not(self%available)) == 0) then
         if (allocated(self%failure)) deallocate (self%failure)
      end if
   end subroutine gaussian_moment_set

   !> Read selected raw outputs with shape and validity checks
   !>
   !> @param[in] self Scientific request
   !> @param[out] gt Raw gt array on the evaluation grid
   !> @param[out] pt Raw pt array on the evaluation grid
   !> @param[out] mt Raw mt array on the evaluation grid
   !> @param[out] rt Raw rt array on the evaluation grid
   !> @param[out] error Missing or invalid output
   subroutine gaussian_moment_get(self, gt, pt, mt, rt, error)
      !> Scientific request
      class(gaussian_moment_request_type), intent(in) :: self
      !> Raw gt array
      real(wp), allocatable, intent(out), optional :: gt(:)
      !> Raw pt array
      real(wp), allocatable, intent(out), optional :: pt(:, :)
      !> Raw mt array
      real(wp), allocatable, intent(out), optional :: mt(:, :, :)
      !> Raw rt array
      real(wp), allocatable, intent(out), optional :: rt(:, :)
      !> Missing or invalid output
      type(error_type), allocatable, intent(out) :: error
      if (present(gt)) then
         if (iand(self%available, output_gt) == 0) then
            call fatal_error(error, "gaussian_moments: missing required output gt")
            return
         end if
         gt = self%gt
      end if
      if (present(pt)) then
         if (iand(self%available, output_pt) == 0) then
            call fatal_error(error, "gaussian_moments: missing required output pt")
            return
         end if
         pt = self%pt
      end if
      if (present(mt)) then
         if (iand(self%available, output_mt) == 0) then
            call fatal_error(error, "gaussian_moments: missing required output mt")
            return
         end if
         mt = self%mt
      end if
      if (present(rt)) then
         if (iand(self%available, output_rt) == 0) then
            call fatal_error(error, "gaussian_moments: missing required output rt")
            return
         end if
         rt = self%rt
      end if
   end subroutine gaussian_moment_get

   !> Start a fresh declaration pass while retaining charge-independent answers
   !>
   !> @param[in,out] self Collection to declare
   subroutine begin_registration(self)
      !> Collection to declare
      class(coupling_type), intent(inout) :: self
      integer :: i
      if (allocated(self%registrations)) deallocate (self%registrations)
      allocate (self%registrations(0))
      if (.not. allocated(self%requests)) allocate (self%requests(0))
      do i = 1, size(self%requests)
         self%requests(i)%item%required = 0
      end do
   end subroutine begin_registration

   !> Select the component whose declarations follow
   !>
   !> @param[in,out] self Collection to declare
   !> @param[in] scope Component index, zero for the cavity
   subroutine set_scope(self, scope)
      !> Collection to declare
      class(coupling_type), intent(inout) :: self
      !> Component index
      integer, intent(in) :: scope
      self%scope = scope
   end subroutine set_scope

   !> Register by local name, sharing calculations only when their inputs match
   !>
   !> @param[in,out] self Collection receiving the declaration
   !> @param[in] key Component-local scientific name
   !> @param[in] item Calculation and output requirements
   !> @param[out] error Registration error
   subroutine coupling_register(self, key, item, error)
      !> Collection receiving the declaration
      class(coupling_type), intent(inout) :: self
      !> Component-local scientific name
      character(len=*), intent(in) :: key
      !> Calculation and requirements
      class(coupling_request_type), intent(in) :: item
      !> Registration error
      type(error_type), allocatable, intent(out) :: error
      type(request_slot), allocatable :: grown(:)
      integer :: i, slot, stat
      logical :: equal
      if (.not. allocated(self%requests)) allocate (self%requests(0))
      if (.not. allocated(self%registrations)) allocate (self%registrations(0))
      if (len_trim(key) > request_name_len .or. len_trim(key) == 0) then
         call fatal_error(error, "Invalid local request name")
         return
      end if
      slot = 0
      do i = 1, size(self%requests)
         equal = self%requests(i)%item%same_inputs(item)
         if (equal) then
            slot = i
            exit
         end if
      end do
      if (slot == 0) then
         slot = size(self%requests) + 1
         allocate (grown(slot), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot allocate request registry")
            return
         end if
         allocate (grown(slot)%item, source=item, stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot allocate scientific request")
            return
         end if
         do i = 1, slot - 1
            call move_alloc(self%requests(i)%item, grown(i)%item)
         end do
         call move_alloc(grown, self%requests)
      else
         self%requests(slot)%item%required = ior(self%requests(slot)%item%required, item%required)
         if (any(self%registrations%slot == slot)) then
            self%requests(slot)%item%fields = ior(self%requests(slot)%item%fields, item%fields)
         else
            self%requests(slot)%item%fields = item%fields
         end if
      end if
      self%registrations = [self%registrations, registration_type(self%scope, key, slot)]
   end subroutine coupling_register
   
   !> Union of grid fields used by active declarations
   !>
   !> @param[in] self Collection to inspect
   function grid_fields(self) result(fields)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Grid field bit mask
      integer :: fields
      integer :: i
      fields = 0
      if (.not. allocated(self%requests)) return
      do i = 1, size(self%requests)
         fields = ior(fields, self%requests(i)%item%fields)
      end do
   end function grid_fields

   !> Copy only declared grid fields; absent unused fields are valid
   !>
   !> @param[in,out] self Collection to refresh
   !> @param[in] ngrid Evaluation point count
   !> @param[out] error Missing declared field or incorrect shape
   !> @param[in] xyz Evaluation points
   !> @param[in] normal Outward normals
   !> @param[in] xi Grid Gaussian inverse lengths
   !> @param[in] f Switching factors
   !> @param[in] area Surface areas
   subroutine coupling_snapshot(self, ngrid, error, xyz, normal, xi, f, area)
      !> Collection to refresh
      class(coupling_type), intent(inout) :: self
      !> Evaluation point count
      integer, intent(in) :: ngrid
      !> Missing declared field or incorrect shape
      type(error_type), allocatable, intent(out) :: error
      !> Optional xyz grid input
      real(wp), intent(in), optional :: xyz(:, :)
      !> Optional normal grid input
      real(wp), intent(in), optional :: normal(:, :)
      !> Optional xi grid input
      real(wp), intent(in), optional :: xi(:)
      !> Optional f grid input
      real(wp), intent(in), optional :: f(:)
      !> Optional area grid input
      real(wp), intent(in), optional :: area(:)
      integer :: fields, i
      self%grid_valid = .false.
      if (ngrid < 0) then
         call fatal_error(error, "Grid size must be nonnegative")
         return
      end if
      call self%compact()
      fields = self%grid_fields()
      if (self%ngrid /= ngrid) call self%invalidate()
      if (iand(fields, grid_xyz) /= 0) then
         if (.not. same_grid_2(self%grid%xyz, xyz)) call self%invalidate()
      end if
      if (iand(fields, grid_normal) /= 0) then
         if (.not. same_grid_2(self%grid%normal, normal)) call self%invalidate()
      end if
      if (iand(fields, grid_xi) /= 0) then
         if (.not. same_grid_1(self%grid%xi, xi)) call self%invalidate()
      end if
      if (iand(fields, grid_f) /= 0) then
         if (.not. same_grid_1(self%grid%f, f)) call self%invalidate()
      end if
      if (iand(fields, grid_area) /= 0) then
         if (.not. same_grid_1(self%grid%area, area)) call self%invalidate()
      end if

      self%ngrid = ngrid
      if (allocated(self%grid%xyz)) deallocate (self%grid%xyz)
      if (iand(fields, grid_xyz) /= 0) then
         if (.not. present(xyz)) then
            call fatal_error(error, "Declared grid field xyz is unavailable")
            return
         end if
         if (any(shape(xyz) /= [3, ngrid])) then
            call fatal_error(error, "Grid field xyz has incorrect shape")
            return
         end if
         self%grid%xyz = xyz
      end if
      if (allocated(self%grid%normal)) deallocate (self%grid%normal)
      if (iand(fields, grid_normal) /= 0) then
         if (.not. present(normal)) then
            call fatal_error(error, "Declared grid field normal is unavailable")
            return
         end if
         if (any(shape(normal) /= [3, ngrid])) then
            call fatal_error(error, "Grid field normal has incorrect shape")
            return
         end if
         self%grid%normal = normal
      end if
      if (allocated(self%grid%xi)) deallocate (self%grid%xi)
      if (iand(fields, grid_xi) /= 0) then
         if (.not. present(xi)) then
            call fatal_error(error, "Declared grid field xi is unavailable")
            return
         end if
         if (any(shape(xi) /= [ngrid])) then
            call fatal_error(error, "Grid field xi has incorrect shape")
            return
         end if
         self%grid%xi = xi
      end if
      if (allocated(self%grid%f)) deallocate (self%grid%f)
      if (iand(fields, grid_f) /= 0) then
         if (.not. present(f)) then
            call fatal_error(error, "Declared grid field f is unavailable")
            return
         end if
         if (any(shape(f) /= [ngrid])) then
            call fatal_error(error, "Grid field f has incorrect shape")
            return
         end if
         self%grid%f = f
      end if
      if (allocated(self%grid%area)) deallocate (self%grid%area)
      if (iand(fields, grid_area) /= 0) then
         if (.not. present(area)) then
            call fatal_error(error, "Declared grid field area is unavailable")
            return
         end if
         if (any(shape(area) /= [ngrid])) then
            call fatal_error(error, "Grid field area has incorrect shape")
            return
         end if
         self%grid%area = area
      end if
      do i = 1, self%n_requests()
         self%requests(i)%item%ngrid = ngrid
      end do
      self%grid_valid = .true.
   end subroutine coupling_snapshot

   !> Whether a current geometry snapshot is available, including an empty grid
   !>
   !> @param[in] self Collection to inspect
   function has_snapshot(self) result(valid)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Snapshot validity
      logical :: valid
      valid = self%grid_valid
   end function has_snapshot

   !> Invalidate every answer and every outstanding host token or scoped view
   !>
   !> @param[in,out] self Collection to invalidate
   subroutine coupling_invalidate(self)
      !> Collection to invalidate
      class(coupling_type), intent(inout) :: self
      integer :: i
      self%epoch = self%epoch + 1_c_int64_t
      self%phase = moist_phase_none
      self%grid_valid = .false.
      do i = 1, self%n_requests()
         call self%requests(i)%item%invalidate()
         self%requests(i)%item%handle = 0_c_int64_t
      end do
      if (allocated(self%pending_index)) deallocate (self%pending_index)
   end subroutine coupling_invalidate

   !> Stage missing outputs and mint fresh, never-reused answer tokens
   !>
   !> @param[in,out] self Collection to stage
   !> @param[in] phase Requested phase
   !> @param[out] error Invalid phase or exhausted token space
   subroutine coupling_arm(self, phase, error)
      !> Collection to stage
      class(coupling_type), intent(inout) :: self
      !> Requested phase
      integer, intent(in) :: phase
      !> Invalid phase or exhausted token space
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      if (phase < 1 .or. phase > moist_n_phases) then
         call fatal_error(error, "Invalid coupling phase")
         return
      end if
      if (self%serial > huge(self%serial) - int(self%n_requests(), c_int64_t)) then
         call fatal_error(error, "Coupling handle space exhausted")
         return
      end if
      self%phase = phase
      self%epoch = self%epoch + 1_c_int64_t
      if (allocated(self%pending_index)) deallocate (self%pending_index)
      allocate (self%pending_index(0))
      do i = 1, self%n_requests()
         self%serial = self%serial + 1_c_int64_t
         associate (item => self%requests(i)%item)
            item%handle = self%serial
            item%ngrid = self%ngrid
            item%missing = iand(item%required(phase), not(item%available))
            if (item%missing /= 0) self%pending_index = [self%pending_index, i]
         end associate
      end do
   end subroutine coupling_arm

   !> Number of registered calculations
   !>
   !> @param[in] self Collection to inspect
   function n_requests(self) result(n)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Calculation count
      integer :: n
      n = 0
      if (allocated(self%requests)) n = size(self%requests)
   end function n_requests

   !> Size of the frozen pending work list
   !>
   !> @param[in] self Collection to inspect
   function n_pending(self) result(n)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Pending count
      integer :: n
      n = 0
      if (allocated(self%pending_index)) n = size(self%pending_index)
   end function n_pending

   !> Copy one declaration for host inspection
   !>
   !> @param[in] self Collection to inspect
   !> @param[in] i One-based declaration index
   function declared(self, i) result(item)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> One-based declaration index
      integer, intent(in) :: i
      !> Independent request value, absent for an invalid index
      class(coupling_request_type), allocatable :: item
      if (i < 1 .or. i > self%n_requests()) return
      allocate (item, source=self%requests(i)%item)
      item%touched = 0
   end function declared

   !> Copy one pending calculation; submit its answer with answer()
   !>
   !> @param[in] self Collection to inspect
   !> @param[in] i One-based pending index
   function pending(self, i) result(item)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> One-based pending index
      integer, intent(in) :: i
      !> Independent request value
      class(coupling_request_type), allocatable :: item
      if (i < 1 .or. i > self%n_pending()) return
      allocate (item, source=self%requests(self%pending_index(i))%item)
      item%touched = 0
   end function pending

   !> Resolve a current host token, refusing every superseded staging
   !>
   !> @param[in] self Collection to inspect
   !> @param[in] handle Opaque token from pending()
   !> @param[out] item Independent scientific request
   !> @param[out] error Superseded or unknown token
   subroutine lookup(self, handle, item, error)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Opaque token
      integer(c_int64_t), intent(in) :: handle
      !> Scientific request copy
      class(coupling_request_type), allocatable, intent(out) :: item
      !> Invalid token error
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      if (handle /= 0_c_int64_t) then
         do i = 1, self%n_requests()
            if (self%requests(i)%item%handle /= handle) cycle
            allocate (item, source=self%requests(i)%item)
            item%touched = 0
            return
         end do
      end if
      call fatal_error(error, "Superseded or unknown coupling request handle")
   end subroutine lookup

   !> Submit a raw answer only into the staging that issued its handle
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] item Answered request value
   !> @param[out] error Rejected handle, type, or scientific output
   subroutine answer(self, item, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Answered request value
      class(coupling_request_type), intent(in) :: item
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      integer :: i, required(moist_n_phases), fields
      do i = 1, self%n_requests()
         if (item%handle == 0_c_int64_t) exit
         if (item%handle /= self%requests(i)%item%handle) cycle
         if (.not. same_type_as(item, self%requests(i)%item)) then
            self%requests(i)%item%available = 0
            self%requests(i)%item%missing = self%requests(i)%item%required(self%phase)
            call fatal_error(error, "Request handle and answer type disagree")
            return
         end if
         required = self%requests(i)%item%required
         fields = self%requests(i)%item%fields
         if (allocated(item%failure)) then
            self%requests(i)%item%available = iand(self%requests(i)%item%available, &
               & not(iand(item%touched, not(item%available))))
            call self%requests(i)%item%merge(item)
            self%requests(i)%item%missing = iand(required(self%phase), &
               & not(self%requests(i)%item%available))
            self%requests(i)%item%failure = item%failure
            call fatal_error(error, item%failure%message)
            return
         end if
         call self%requests(i)%item%merge(item)
         if (allocated(self%requests(i)%item%failure)) deallocate (self%requests(i)%item%failure)
         self%requests(i)%item%required = required
         self%requests(i)%item%fields = fields
         self%requests(i)%item%missing = iand(required(self%phase), &
            & not(self%requests(i)%item%available))
         return
      end do
      call fatal_error(error, "Superseded or unknown coupling request handle")
   end subroutine answer

   !> Check completeness against the active phase, never against optional stale data
   !>
   !> @param[in] self Collection to validate
   !> @param[in] phase Phase to consume
   !> @param[out] error Missing scientific output or unstaged collection
   subroutine check_mandatory(self, phase, error)
      !> Collection to validate
      class(coupling_type), intent(in) :: self
      !> Phase to consume
      integer, intent(in) :: phase
      !> Missing scientific output
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      if (phase < 1 .or. phase > moist_n_phases .or. phase /= self%phase) then
         call fatal_error(error, "Coupling is not staged for this phase")
         return
      end if
      do i = 1, self%n_requests()
         associate (item => self%requests(i)%item)
            if (iand(item%required(phase), not(item%available)) == 0) cycle
            if (allocated(item%failure)) then
               call fatal_error(error, trim(item%name())//": missing required outputs: " &
                  & //item%failure%message)
            else
               call fatal_error(error, trim(item%name())//": missing required outputs")
            end if
            return
         end associate
      end do
   end subroutine check_mandatory

   !> Create a read-only component scope for one call
   !>
   !> @param[in,out] self Collection owning the answers
   !> @param[in] scope Component registration index
   !> @param[out] view Scoped read interface
   subroutine make_view(self, scope, view)
      !> Collection owning the answers
      class(coupling_type), intent(inout), target :: self
      !> Component registration index
      integer, intent(in) :: scope
      !> Scoped read interface
      type(coupling_view_type), intent(out) :: view
      self%epoch = self%epoch + 1_c_int64_t
      view%owner => self
      view%scope = scope
      view%epoch = self%epoch
      view%phase = self%phase
   end subroutine make_view

   !> End the component call and invalidate any captured value copy of its view
   !>
   !> @param[in,out] self Collection owning the answers
   subroutine close_view(self)
      !> Collection owning the answers
      class(coupling_type), intent(inout) :: self
      self%epoch = self%epoch + 1_c_int64_t
   end subroutine close_view

   !> Validate a scoped read
   !>
   !> @param[in] self Component view
   !> @param[in] phase Phase to consume
   !> @param[out] error Expired view or missing outputs
   subroutine view_check(self, phase, error)
      !> Component view
      class(coupling_view_type), intent(in) :: self
      !> Phase to consume
      integer, intent(in) :: phase
      !> Expired view or missing outputs
      type(error_type), allocatable, intent(out) :: error
      if (.not. associated(self%owner)) then
         call fatal_error(error, "Unbound component coupling view")
         return
      end if
      if (self%epoch /= self%owner%epoch) then
         call fatal_error(error, "Expired component coupling view")
         return
      end if
      call self%owner%check_mandatory(phase, error)
   end subroutine view_check

   !> Copy a complete request for inspection within its registration scope
   !>
   !> Use read_potential or read_moments to copy only selected output arrays
   !>
   !> @param[in] self Component view
   !> @param[in] key Local name used in register()
   !> @param[out] item Independent request value
   !> @param[out] error Expired view or undeclared name
   subroutine view_read(self, key, item, error)
      !> Component view
      class(coupling_view_type), intent(in) :: self
      !> Local registration name
      character(len=*), intent(in) :: key
      !> Independent request value
      class(coupling_request_type), allocatable, intent(out) :: item
      !> Read error
      type(error_type), allocatable, intent(out) :: error
      !> Resolved request slot
      integer :: slot
      slot = self%find_slot(key, error)
      if (allocated(error)) return
      allocate (item, source=self%owner%requests(slot)%item)
      item%touched = 0
      item%handle = 0_c_int64_t
   end subroutine view_read

   !> Resolve a local registration without copying its answer arrays
   !>
   !> @param[in] self Component view
   !> @param[in] key Local registration name
   !> @param[out] error Expired view or undeclared name
   function view_find_slot(self, key, error) result(slot)
      !> Component view
      class(coupling_view_type), intent(in) :: self
      !> Local registration name
      character(len=*), intent(in) :: key
      !> Read error
      type(error_type), allocatable, intent(out) :: error
      !> Resolved slot, zero on failure
      integer :: slot
      !> Registration index
      integer :: i
      slot = 0
      call self%check_mandatory(self%phase, error)
      if (allocated(error)) return
      if (.not. allocated(self%owner%registrations)) then
         call fatal_error(error, "Component has no registered requests")
         return
      end if
      do i = 1, size(self%owner%registrations)
         associate (reg => self%owner%registrations(i))
            if (reg%scope /= self%scope .or. reg%key /= key) cycle
            if (reg%slot < 1 .or. reg%slot > self%owner%n_requests()) exit
            slot = reg%slot
            return
         end associate
      end do
      call fatal_error(error, "Component did not register request '"//key//"'")
   end function view_find_slot

   !> Copy selected potential outputs directly from the registered request
   !>
   !> @param[in] self Component view
   !> @param[out] error Missing or incompatible calculation
   !> @param[out] phi Potential on the grid
   !> @param[out] dphi_dr Spatial derivative of the potential
   !> @param[out] dphi_dxi Gaussian inverse-length derivative
   subroutine read_potential(self, error, phi, dphi_dr, dphi_dxi)
      !> Component view
      class(coupling_view_type), intent(in) :: self
      !> Read error
      type(error_type), allocatable, intent(out) :: error
      !> Potential on the grid
      real(wp), allocatable, optional, intent(out) :: phi(:)
      !> Spatial derivative
      real(wp), allocatable, optional, intent(out) :: dphi_dr(:, :)
      !> Inverse-length derivative
      real(wp), allocatable, optional, intent(out) :: dphi_dxi(:)
      !> Resolved request slot
      integer :: slot
      slot = self%find_slot("potential", error)
      if (allocated(error)) return
      select type (item => self%owner%requests(slot)%item)
      type is (point_potential_request_type)
         if (present(dphi_dxi)) then
            call fatal_error(error, "Point potential has no Gaussian width derivative")
            return
         end if
         call item%get(phi, dphi_dr, error)
      type is (gaussian_potential_request_type)
         call item%get(phi, dphi_dr, dphi_dxi, error)
      class default
         call fatal_error(error, "Registered potential has an incompatible calculation type")
      end select
   end subroutine read_potential

   !> Copy selected Gaussian moments directly from the registered request
   !>
   !> @param[in] self Component view
   !> @param[out] error Missing or incompatible calculation
   !> @param[out] gt Scalar moments
   !> @param[out] pt Vector moments
   !> @param[out] mt Second moments
   !> @param[out] rt Third contracted moments
   subroutine read_moments(self, error, gt, pt, mt, rt)
      !> Component view
      class(coupling_view_type), intent(in) :: self
      !> Read error
      type(error_type), allocatable, intent(out) :: error
      !> Scalar moments
      real(wp), allocatable, optional, intent(out) :: gt(:)
      !> Vector moments
      real(wp), allocatable, optional, intent(out) :: pt(:, :)
      !> Second moments
      real(wp), allocatable, optional, intent(out) :: mt(:, :, :)
      !> Third contracted moments
      real(wp), allocatable, optional, intent(out) :: rt(:, :)
      !> Resolved request slot
      integer :: slot
      slot = self%find_slot("moments", error)
      if (allocated(error)) return
      select type (item => self%owner%requests(slot)%item)
      type is (gaussian_moment_request_type)
         call item%get(gt, pt, mt, rt, error)
      class default
         call fatal_error(error, "Registered moments have an incompatible calculation type")
      end select
   end subroutine read_moments

   !> Allocate a model-owned coupling and return a borrowed pointer
   !>
   !> @param[in,out] self Model registry
   !> @param[out] coupling Borrowed collection
   !> @param[out] error Allocation failure
   subroutine mint(self, coupling, error)
      !> Model registry
      class(coupling_registry_type), intent(inout), target :: self
      !> Borrowed collection
      type(coupling_type), pointer, intent(out) :: coupling
      !> Allocation failure
      type(error_type), allocatable, intent(out) :: error
      type(coupling_node), allocatable :: node
      integer :: stat
      nullify (coupling)
      allocate (node, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate coupling registry node")
         return
      end if
      allocate (node%item, stat=stat)
      if (stat /= 0) then
         deallocate (node)
         call fatal_error(error, "Cannot allocate coupling")
         return
      end if
      call move_alloc(self%first, node%next)
      call move_alloc(node, self%first)
      coupling => self%first%item
   end subroutine mint

   !> Release one collection and nullify its borrowed pointer
   !>
   !> @param[in,out] self Model registry
   !> @param[in,out] coupling Collection to release; other aliases become invalid
   subroutine registry_release(self, coupling)
      !> Model registry
      class(coupling_registry_type), intent(inout), target :: self
      !> Collection to release
      type(coupling_type), pointer, intent(inout) :: coupling
      !> Node before the one being removed
      type(coupling_node), pointer :: previous
      !> Detached node; its collection is freed on return
      type(coupling_node), allocatable :: removed
      if (.not. associated(coupling)) return
      if (.not. allocated(self%first)) return
      if (associated(coupling, self%first%item)) then
         call move_alloc(self%first, removed)
         call move_alloc(removed%next, self%first)
         nullify (coupling)
         return
      end if
      previous => self%first
      do while (allocated(previous%next))
         if (associated(coupling, previous%next%item)) then
            call move_alloc(previous%next, removed)
            call move_alloc(removed%next, previous%next)
            nullify (coupling)
            return
         end if
         previous => previous%next
      end do
   end subroutine registry_release

   !> Push invalidation into every coupling issued by this model
   !>
   !> @param[in,out] self Model registry
   subroutine registry_invalidate(self)
      !> Model registry
      class(coupling_registry_type), intent(inout), target :: self
      type(coupling_node), pointer :: node
      if (.not. allocated(self%first)) return
      node => self%first
      do

         call node%item%invalidate()
         if (.not. allocated(node%next)) exit
         node => node%next
      end do
   end subroutine registry_invalidate

   !> Check model ownership before staging or consuming a borrowed collection
   !>
   !> @param[in] self Model registry
   !> @param[in] coupling Borrowed collection
   function owns(self, coupling) result(found)
      !> Model registry
      class(coupling_registry_type), intent(in), target :: self
      !> Borrowed collection
      type(coupling_type), intent(in), target :: coupling
      !> Whether this model minted the collection
      logical :: found
      type(coupling_node), pointer :: node
      type(coupling_type), pointer :: candidate
      found = .false.
      if (.not. allocated(self%first)) return
      node => self%first
      do

         candidate => node%item
         if (associated(candidate, coupling)) then
            found = .true.
            return
         end if
         if (.not. allocated(node%next)) exit
         node => node%next
      end do
   end function owns

   !> Release all owned collections when their model is destroyed
   !>
   !> @param[in,out] self Model registry
   subroutine registry_clear(self)
      !> Model registry
      class(coupling_registry_type), intent(inout), target :: self
      type(coupling_node), allocatable :: node
      do while (allocated(self%first))
         call move_alloc(self%first, node)
         call move_alloc(node%next, self%first)
         deallocate (node)
      end do
   end subroutine registry_clear

   !> Default input identity for calculations with no private scientific inputs
   !>
   !> @param[in] self Registered calculation
   !> @param[in] other Candidate calculation
   function request_same_inputs(self, other) result(equal)
      !> Registered calculation
      class(coupling_request_type), intent(in) :: self
      !> Candidate calculation
      class(coupling_request_type), intent(in) :: other
      !> Whether both calculations have identical inputs
      logical :: equal
      equal = same_type_as(self, other)
   end function request_same_inputs

   !> Moment calculations share answers only when their exponents match exactly
   !>
   !> @param[in] self Registered moments
   !> @param[in] other Candidate calculation
   function moment_same_inputs(self, other) result(equal)
      !> Registered moments
      class(gaussian_moment_request_type), intent(in) :: self
      !> Candidate calculation
      class(coupling_request_type), intent(in) :: other
      !> Whether inputs match
      logical :: equal
      equal = .false.
      select type (other)
      type is (gaussian_moment_request_type)
         if (.not. allocated(self%width) .or. .not. allocated(other%width)) return
         if (size(self%width) /= size(other%width)) return
         equal = all(self%width == other%width)
      end select
   end function moment_same_inputs

   !> Remove calculations no longer used and remap local registrations
   !>
   !> @param[in,out] self Collection after component declaration
   subroutine compact(self)
      !> Collection after component declaration
      class(coupling_type), intent(inout) :: self
      type(request_slot), allocatable :: active(:)
      integer, allocatable :: mapping(:)
      integer :: i, n
      if (.not. allocated(self%requests)) return
      if (.not. allocated(self%registrations)) return
      allocate (mapping(size(self%requests)), source=0)
      n = 0
      do i = 1, size(self%requests)
         if (.not. any(self%registrations%slot == i)) cycle
         n = n + 1
         mapping(i) = n
      end do
      if (n == size(self%requests)) return
      allocate (active(n))
      do i = 1, size(self%requests)
         if (mapping(i) == 0) cycle
         call move_alloc(self%requests(i)%item, active(mapping(i))%item)
      end do
      call move_alloc(active, self%requests)
      do i = 1, size(self%registrations)
         self%registrations(i)%slot = mapping(self%registrations(i)%slot)
      end do
      self%registrations = pack(self%registrations, self%registrations%slot > 0)
   end subroutine compact

   !> Merge partial outputs without overwriting independently supplied answers
   !>
   !> @param[in,out] self Live scientific request
   !> @param[in] answer Host answer value
   subroutine point_potential_merge(self, answer)
      !> Live scientific request
      class(point_potential_request_type), intent(inout) :: self
      !> Host answer value
      class(coupling_request_type), intent(in) :: answer
      select type (answer)
      type is (point_potential_request_type)
         if (iand(iand(answer%touched, answer%available), output_phi) /= 0) then
            self%phi = answer%phi
            self%available = ior(self%available, output_phi)
         end if
         if (iand(iand(answer%touched, answer%available), output_dphi_dr) /= 0) then
            self%dphi_dr = answer%dphi_dr
            self%available = ior(self%available, output_dphi_dr)
         end if
      end select
   end subroutine point_potential_merge

   !> Merge partial outputs without overwriting independently supplied answers
   !>
   !> @param[in,out] self Live scientific request
   !> @param[in] answer Host answer value
   subroutine gaussian_potential_merge(self, answer)
      !> Live scientific request
      class(gaussian_potential_request_type), intent(inout) :: self
      !> Host answer value
      class(coupling_request_type), intent(in) :: answer
      select type (answer)
      type is (gaussian_potential_request_type)
         if (iand(iand(answer%touched, answer%available), output_phi) /= 0) then
            self%phi = answer%phi
            self%available = ior(self%available, output_phi)
         end if
         if (iand(iand(answer%touched, answer%available), output_dphi_dr) /= 0) then
            self%dphi_dr = answer%dphi_dr
            self%available = ior(self%available, output_dphi_dr)
         end if
         if (iand(iand(answer%touched, answer%available), output_dphi_dxi) /= 0) then
            self%dphi_dxi = answer%dphi_dxi
            self%available = ior(self%available, output_dphi_dxi)
         end if
      end select
   end subroutine gaussian_potential_merge

   !> Merge partial outputs without overwriting independently supplied answers
   !>
   !> @param[in,out] self Live scientific request
   !> @param[in] answer Host answer value
   subroutine gaussian_moment_merge(self, answer)
      !> Live scientific request
      class(gaussian_moment_request_type), intent(inout) :: self
      !> Host answer value
      class(coupling_request_type), intent(in) :: answer
      select type (answer)
      type is (gaussian_moment_request_type)
         if (iand(iand(answer%touched, answer%available), output_gt) /= 0) then
            self%gt = answer%gt
            self%available = ior(self%available, output_gt)
         end if
         if (iand(iand(answer%touched, answer%available), output_pt) /= 0) then
            self%pt = answer%pt
            self%available = ior(self%available, output_pt)
         end if
         if (iand(iand(answer%touched, answer%available), output_mt) /= 0) then
            self%mt = answer%mt
            self%available = ior(self%available, output_mt)
         end if
         if (iand(iand(answer%touched, answer%available), output_rt) /= 0) then
            self%rt = answer%rt
            self%available = ior(self%available, output_rt)
         end if
      end select
   end subroutine gaussian_moment_merge

   !> Compare an optional grid input with its preceding snapshot
   !>
   !> @param[in] previous Owned grid snapshot
   !> @param[in] current Optional input for the next staging
   function same_grid_1(previous, current) result(equal)
      !> Owned grid snapshot
      real(wp), allocatable, intent(in) :: previous(:)
      !> Optional current grid input
      real(wp), optional, intent(in) :: current(:)
      !> Whether values and extents agree exactly
      logical :: equal
      equal = .false.
      if (.not. allocated(previous) .or. .not. present(current)) return
      if (any(shape(previous) /= shape(current))) return
      equal = all(previous == current)
   end function same_grid_1

   !> Compare an optional grid input with its preceding snapshot
   !>
   !> @param[in] previous Owned grid snapshot
   !> @param[in] current Optional input for the next staging
   function same_grid_2(previous, current) result(equal)
      !> Owned grid snapshot
      real(wp), allocatable, intent(in) :: previous(:, :)
      !> Optional current grid input
      real(wp), optional, intent(in) :: current(:, :)
      !> Whether values and extents agree exactly
      logical :: equal
      equal = .false.
      if (.not. allocated(previous) .or. .not. present(current)) return
      if (any(shape(previous) /= shape(current))) return
      equal = all(previous == current)
   end function same_grid_2

   !> Copy declared grid inputs with availability and lifetime checks
   !>
   !> @param[in] self Host collection
   !> @param[out] error Unavailable or undeclared grid field
   !> @param[out] xyz Optional grid field to copy
   !> @param[out] normal Optional grid field to copy
   !> @param[out] xi Optional grid field to copy
   !> @param[out] f Optional grid field to copy
   !> @param[out] area Optional grid field to copy
   subroutine read_grid(self, error, xyz, normal, xi, f, area)
      !> Host collection
      class(coupling_type), intent(in) :: self
      !> Missing grid field
      type(error_type), allocatable, intent(out) :: error
      !> Requested xyz values
      real(wp), allocatable, optional, intent(out) :: xyz(:, :)
      !> Requested normal values
      real(wp), allocatable, optional, intent(out) :: normal(:, :)
      !> Requested xi values
      real(wp), allocatable, optional, intent(out) :: xi(:)
      !> Requested f values
      real(wp), allocatable, optional, intent(out) :: f(:)
      !> Requested area values
      real(wp), allocatable, optional, intent(out) :: area(:)
      if (.not. self%grid_valid) then
         call fatal_error(error, "Grid snapshot is stale or unavailable")
         return
      end if
      if (present(xyz)) then
         if (.not. allocated(self%grid%xyz)) then
            call fatal_error(error, "Grid field xyz is not declared or unavailable")
            return
         end if
         xyz = self%grid%xyz
      end if
      if (present(normal)) then
         if (.not. allocated(self%grid%normal)) then
            call fatal_error(error, "Grid field normal is not declared or unavailable")
            return
         end if
         normal = self%grid%normal
      end if
      if (present(xi)) then
         if (.not. allocated(self%grid%xi)) then
            call fatal_error(error, "Grid field xi is not declared or unavailable")
            return
         end if
         xi = self%grid%xi
      end if
      if (present(f)) then
         if (.not. allocated(self%grid%f)) then
            call fatal_error(error, "Grid field f is not declared or unavailable")
            return
         end if
         f = self%grid%f
      end if
      if (present(area)) then
         if (.not. allocated(self%grid%area)) then
            call fatal_error(error, "Grid field area is not declared or unavailable")
            return
         end if
         area = self%grid%area
      end if
   end subroutine read_grid

end module moist_channels_request
