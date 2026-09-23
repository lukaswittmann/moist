!> Host -> MOIST coupling requests and their answers
!>
!> A coupling is the model-owned list of calculations the host performs on the
!> cavity grid, one request per distinct set of inputs
!>
!> - a request is a declaration: named outputs such as `phi` or `dphi_dr(3, ngrid)`
!>   and the cavity fields it reads; it holds no data, so copies are cheap
!> - each phase stages the outputs still missing and mints one handle per
!>   request; the host answers by handle and output name,
!>   `coupling%answer(request%handle, "phi", phi, error)`, and a superseded
!>   handle is refused
!> - the coupling owns the answer arrays; they survive across phases until the
!>   owner invalidates them, and a rejected answer leaves its output missing
!> - components never see the coupling; each reads its own registrations through
!>   a short-lived `coupling_view_type` after one completeness check per call
module moist_channels_coupling
   use, intrinsic :: iso_c_binding, only: c_int64_t
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   implicit none(type, external)
   private
   public :: coupling_type, coupling_view_type, coupling_registry_type
   public :: coupling_request_type, output_slot_type
   public :: point_potential_request_type, gaussian_potential_request_type
   public :: gaussian_moment_request_type
   public :: moist_phase_none, moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, moist_n_phases
   public :: request_name_len, output_name_len, grid_name_len

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
   !> Output name capacity
   integer, parameter :: output_name_len = 16
   !> Grid field name capacity; grid fields are cavity fields, so this matches
   !> MOIST_FIELD_NAME_MAX in moist.h
   integer, parameter :: grid_name_len = 64
   !> Leading extents an output may declare before ngrid
   integer, parameter :: max_lead = 2

   !> One declared output of a request; the coupling owns its answer array
   !>
   !> Answers are stored as (rows, ngrid) with rows = product(lead); the
   !> declared leading extents restore the caller's rank on read
   type :: output_slot_type
      !> Scientific name, e.g. "phi" or "dphi_dr"
      character(len=output_name_len) :: name = ""
      !> Leading extents before ngrid; unused entries stay 1
      integer :: lead(max_lead) = 1
      !> Number of leading extents: 0 for (ngrid), 1 for (lead(1), ngrid), ...
      integer :: nlead = 0
      !> Whether each phase requires this output
      logical :: required(moist_n_phases) = .false.
      !> Whether the coupling holds a valid answer
      logical :: available = .false.
   end type output_slot_type

   !> Scientific request: what a host computes, never the data it returns
   !>
   !> Kinds declare their outputs once in `declare`; every accessor works by
   !> output name, so a new kind adds no per-output code. Hosts answer through
   !> `coupling%answer(handle, name, values, error)`, so copies of a request are
   !> cheap declarations
   type, abstract :: coupling_request_type
      !> Opaque host token, replaced on every staging or invalidation
      integer(c_int64_t) :: handle = 0_c_int64_t
      !> Named outputs, filled by `declare` on first use
      type(output_slot_type), allocatable :: outputs(:)
      !> Cavity grid fields this calculation reads, by the cavity's field names
      character(len=grid_name_len), allocatable :: fields(:)
      !> Expected grid size
      integer :: ngrid = -1
      !> Phase staged by the latest arm
      integer :: phase = moist_phase_none
   contains
      procedure(request_name), deferred :: name
      procedure(request_declare), deferred :: declare
      procedure :: ensure_outputs => request_ensure_outputs
      procedure, private :: request_add_vector, request_add_array
      generic :: add => request_add_vector, request_add_array
      procedure :: uses => request_uses_grid
      procedure :: n_outputs => request_n_outputs
      procedure :: find => request_find
      procedure, private :: request_require_always, request_require_when
      generic :: require => request_require_always, request_require_when
      procedure :: clear_requirements => request_clear_requirements
      procedure :: is_required => request_is_required
      procedure :: is_available => request_is_available
      procedure :: is_missing => request_is_missing
      procedure :: n_missing => request_n_missing
      procedure :: rows => request_rows
      procedure :: invalidate => request_invalidate
      procedure :: same_inputs => request_same_inputs
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

      !> Declare the kind's outputs and grid fields
      !>
      !> @param[in,out] self Request to declare
      !> @param[out] error Invalid declaration
      subroutine request_declare(self, error)
         import :: coupling_request_type, error_type
         implicit none(type, external)
         !> Request to declare
         class(coupling_request_type), intent(inout) :: self
         !> Invalid declaration
         type(error_type), allocatable, intent(out) :: error
      end subroutine request_declare
   end interface

   !> Raw point_potential calculation: phi(ngrid), dphi_dr(3, ngrid)
   type, extends(coupling_request_type) :: point_potential_request_type
   contains
      procedure :: name => point_potential_name
      procedure :: declare => point_potential_declare
   end type point_potential_request_type

   !> Raw gaussian_potential calculation: phi, dphi_dr(3, ngrid), dphi_dxi
   type, extends(coupling_request_type) :: gaussian_potential_request_type
   contains
      procedure :: name => gaussian_potential_name
      procedure :: declare => gaussian_potential_declare
   end type gaussian_potential_request_type

   !> Raw gaussian_moments calculation: gt, pt(3, ngrid), mt(3, 3, ngrid), rt(3, ngrid)
   type, extends(coupling_request_type) :: gaussian_moment_request_type
      !> Moment exponents in bohr**(-2), supplied by the component
      real(wp), allocatable :: width(:)
   contains
      procedure :: name => gaussian_moment_name
      procedure :: declare => gaussian_moment_declare
      procedure :: same_inputs => moment_same_inputs
   end type gaussian_moment_request_type

   !> Stored answer of one output
   type :: output_values_type
      !> Flattened answer (rows, ngrid)
      real(wp), allocatable :: values(:, :)
   end type output_values_type

   !> Owning slot for a heterogeneous calculation and its answers
   type :: request_slot_type
      !> Concrete calculation
      class(coupling_request_type), allocatable :: item
      !> One answer array per declared output
      type(output_values_type), allocatable :: answers(:)
   end type request_slot_type

   !> Registration local to a cavity or component
   type :: registration_type
      !> Component number; zero denotes the cavity
      integer :: scope = 0
      !> Local scientific name chosen during registration
      character(len=request_name_len) :: key = ""
      !> Registered request slot
      integer :: slot = 0
   end type registration_type
   !> Model-owned collection exposed to hosts through borrowed handles
   type :: coupling_type
      private
      !> Scientific calculations
      type(request_slot_type), allocatable :: requests(:)
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
      !> Whether the declared grid fields are current on their owner
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
      procedure, private :: coupling_answer_1, coupling_answer_2, coupling_answer_3
      generic :: answer => coupling_answer_1, coupling_answer_2, coupling_answer_3
      procedure :: answer_flat => coupling_answer_flat
      procedure :: reject => coupling_reject
      procedure, private :: index_of
      procedure, private :: resolve
      procedure, private :: store
      procedure :: check_mandatory
      procedure :: make_view
      procedure :: close_view
      procedure, private :: compact
      procedure :: grid_fields
      procedure :: uses_grid
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
      procedure :: request => view_request
      procedure :: check_mandatory => view_check
      procedure, private :: validate => view_validate
      procedure, private :: view_read_1, view_read_2, view_read_3
      generic :: read => view_read_1, view_read_2, view_read_3
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

   !> Mark every output unanswered; the coupling keeps the storage
   !>
   !> @param[in,out] self Request to invalidate
   subroutine request_invalidate(self)
      !> Request to invalidate
      class(coupling_request_type), intent(inout) :: self
      integer :: i
      do i = 1, self%n_outputs()
         self%outputs(i)%available = .false.
      end do
   end subroutine request_invalidate

   !> Declare the outputs once; a failed declaration leaves the request undeclared
   !>
   !> @param[in,out] self Request to declare
   !> @param[out] error Invalid declaration
   subroutine request_ensure_outputs(self, error)
      !> Request to declare
      class(coupling_request_type), intent(inout) :: self
      !> Invalid declaration
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(self%outputs)) then
         call self%declare(error)
         if (allocated(error)) then
            if (allocated(self%outputs)) deallocate (self%outputs)
            return
         end if
      end if
      if (.not. allocated(self%fields)) allocate (self%fields(0))
   end subroutine request_ensure_outputs

   !> Declare that this calculation reads one cavity grid field
   !>
   !> @param[in,out] self Request being declared
   !> @param[in] name Cavity field name, e.g. "xyz" or "xi0"
   !> @param[out] error Name longer than a cavity field name can be
   subroutine request_uses_grid(self, name, error)
      !> Request being declared
      class(coupling_request_type), intent(inout) :: self
      !> Cavity field name
      character(len=*), intent(in) :: name
      !> Invalid name
      type(error_type), allocatable, intent(out) :: error
      !> Name padded to the stored length; gfortran ignores the constructor
      !> type-spec for an assumed-length dummy, so `[character(len=n) :: name]`
      !> would keep the actual length
      character(len=grid_name_len) :: field
      if (len_trim(name) > grid_name_len) then
         call fatal_error(error, trim(self%name())//": grid field name '"//name//"' is too long")
         return
      end if
      field = name
      if (.not. allocated(self%fields)) allocate (self%fields(0))
      self%fields = union_names(self%fields, [field])
   end subroutine request_uses_grid

   !> Names of b not yet in a, appended to a in first-seen order
   !>
   !> @param[in] a Existing names
   !> @param[in] b Names to merge
   pure function union_names(a, b) result(union)
      !> Existing names
      character(len=grid_name_len), intent(in) :: a(:)
      !> Names to merge
      character(len=grid_name_len), intent(in) :: b(:)
      !> Merged names
      character(len=grid_name_len), allocatable :: union(:)
      integer :: i
      union = a
      do i = 1, size(b)
         if (.not. any(union == b(i))) union = [union, b(i)]
      end do
   end function union_names

   !> Sort names so that hosts see a stable order across preparations
   !>
   !> @param[in,out] names Names to sort in place
   pure subroutine sort_names(names)
      !> Names to sort
      character(len=grid_name_len), intent(inout) :: names(:)
      character(len=grid_name_len) :: key
      integer :: i, j
      do i = 2, size(names)
         key = names(i)
         j = i - 1
         do while (j >= 1)
            if (names(j) <= key) exit
            names(j + 1) = names(j)
            j = j - 1
         end do
         names(j + 1) = key
      end do
   end subroutine sort_names

   !> Declare one vector output (ngrid); called from a kind's `declare`
   !>
   !> @param[in,out] self Request being declared
   !> @param[in] name Output name
   !> @param[out] error Invalid declaration
   subroutine request_add_vector(self, name, error)
      !> Request being declared
      class(coupling_request_type), intent(inout) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Invalid declaration
      type(error_type), allocatable, intent(out) :: error
      !> No leading extents
      integer :: lead(0)
      call request_add_array(self, name, lead, error)
   end subroutine request_add_vector

   !> Declare one output (lead, ngrid); called from a kind's `declare`
   !>
   !> @param[in,out] self Request being declared
   !> @param[in] name Output name
   !> @param[in] lead Leading extents before ngrid, at most `max_lead`
   !> @param[out] error Invalid declaration
   subroutine request_add_array(self, name, lead, error)
      !> Request being declared
      class(coupling_request_type), intent(inout) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Leading extents before ngrid
      integer, intent(in) :: lead(:)
      !> Invalid declaration
      type(error_type), allocatable, intent(out) :: error
      type(output_slot_type), allocatable :: grown(:)
      integer :: n
      if (len_trim(name) > output_name_len) then
         call fatal_error(error, trim(self%name())//": output name '"//name//"' is too long")
      else if (self%find(name) /= 0) then
         call fatal_error(error, trim(self%name())//": output '"//name//"' declared twice")
      else if (size(lead) > max_lead) then
         call fatal_error(error, &
            & trim(self%name())//": output '"//name//"' has too many leading extents")
      else if (any(lead < 1)) then
         call fatal_error(error, &
            & trim(self%name())//": output '"//name//"' has an empty leading extent")
      end if
      if (allocated(error)) return
      n = self%n_outputs()
      allocate (grown(n + 1))
      if (n > 0) grown(:n) = self%outputs
      grown(n + 1)%name = name
      grown(n + 1)%nlead = size(lead)
      grown(n + 1)%lead(:size(lead)) = lead
      call move_alloc(grown, self%outputs)
   end subroutine request_add_array

   !> Number of declared outputs; zero before `declare`
   !>
   !> @param[in] self Request to inspect
   function request_n_outputs(self) result(n)
      !> Request to inspect
      class(coupling_request_type), intent(in) :: self
      !> Output count
      integer :: n
      n = 0
      if (allocated(self%outputs)) n = size(self%outputs)
   end function request_n_outputs

   !> Index of a named output, zero when unknown
   !>
   !> @param[in] self Request to inspect
   !> @param[in] name Output name
   function request_find(self, name) result(i)
      !> Request to inspect
      class(coupling_request_type), intent(in) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Output index
      integer :: i
      do i = 1, self%n_outputs()
         if (self%outputs(i)%name == name) return
      end do
      i = 0
   end function request_find


   !> Flattened row count of one output
   !>
   !> @param[in] self Request to inspect
   !> @param[in] i Output index
   function request_rows(self, i) result(rows)
      !> Request to inspect
      class(coupling_request_type), intent(in) :: self
      !> Output index
      integer, intent(in) :: i
      !> Rows per grid point
      integer :: rows
      rows = product(self%outputs(i)%lead)
   end function request_rows

   !> Declare one output required in one phase
   !>
   !> @param[in,out] self Request requirements
   !> @param[in] phase Phase index
   !> @param[in] name Output name
   !> @param[out] error Invalid phase, unknown output or failed declaration
   subroutine request_require_always(self, phase, name, error)
      !> Request requirements
      class(coupling_request_type), intent(inout) :: self
      !> Phase index
      integer, intent(in) :: phase
      !> Output name
      character(len=*), intent(in) :: name
      !> Invalid requirement
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      if (phase < 1 .or. phase > moist_n_phases) then
         call fatal_error(error, trim(self%name())//": invalid coupling phase")
         return
      end if
      call self%ensure_outputs(error)
      if (allocated(error)) return
      i = self%find(name)
      if (i == 0) then
         call fatal_error(error, trim(self%name())//" has no output '"//name//"'")
         return
      end if
      self%outputs(i)%required(phase) = .true.
   end subroutine request_require_always

   !> Declare one output required in one phase under a condition
   !>
   !> @param[in,out] self Request requirements
   !> @param[in] phase Phase index
   !> @param[in] name Output name
   !> @param[in] when Whether the requirement applies
   !> @param[out] error Invalid phase, unknown output or failed declaration
   subroutine request_require_when(self, phase, name, when, error)
      !> Request requirements
      class(coupling_request_type), intent(inout) :: self
      !> Phase index
      integer, intent(in) :: phase
      !> Output name
      character(len=*), intent(in) :: name
      !> Whether the requirement applies
      logical, intent(in) :: when
      !> Invalid requirement
      type(error_type), allocatable, intent(out) :: error
      if (when) call request_require_always(self, phase, name, error)
   end subroutine request_require_when

   !> Drop every phase requirement before a new declaration pass
   !>
   !> @param[in,out] self Request requirements
   subroutine request_clear_requirements(self)
      !> Request requirements
      class(coupling_request_type), intent(inout) :: self
      integer :: i
      do i = 1, self%n_outputs()
         self%outputs(i)%required = .false.
      end do
   end subroutine request_clear_requirements

   !> Whether an output is required in a phase, the staged one by default
   !>
   !> @param[in] self Request to inspect
   !> @param[in] name Output name
   !> @param[in] phase Phase index, staged phase when absent
   function request_is_required(self, name, phase) result(required)
      !> Request to inspect
      class(coupling_request_type), intent(in) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Phase index
      integer, intent(in), optional :: phase
      !> Requirement
      logical :: required
      integer :: i, p
      required = .false.
      p = self%phase
      if (present(phase)) p = phase
      if (p < 1 .or. p > moist_n_phases) return
      i = self%find(name)
      if (i > 0) required = self%outputs(i)%required(p)
   end function request_is_required

   !> Whether an output holds a valid answer
   !>
   !> @param[in] self Request to inspect
   !> @param[in] name Output name
   function request_is_available(self, name) result(available)
      !> Request to inspect
      class(coupling_request_type), intent(in) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Availability
      logical :: available
      integer :: i
      available = .false.
      i = self%find(name)
      if (i > 0) available = self%outputs(i)%available
   end function request_is_available

   !> Whether an output is required by the staged phase and still unanswered
   !>
   !> @param[in] self Request to inspect
   !> @param[in] name Output name
   function request_is_missing(self, name) result(missing)
      !> Request to inspect
      class(coupling_request_type), intent(in) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Whether the output is missing
      logical :: missing
      missing = self%is_required(name) .and. .not. self%is_available(name)
   end function request_is_missing

   !> Number of required outputs without valid answers
   !>
   !> @param[in] self Request to inspect
   !> @param[in] phase Phase index, staged phase when absent
   function request_n_missing(self, phase) result(n)
      !> Request to inspect
      class(coupling_request_type), intent(in) :: self
      !> Phase index
      integer, intent(in), optional :: phase
      !> Missing count
      integer :: n
      integer :: i, p
      n = 0
      p = self%phase
      if (present(phase)) p = phase
      if (p < 1 .or. p > moist_n_phases) return
      do i = 1, self%n_outputs()
         if (self%outputs(i)%required(p) .and. .not. self%outputs(i)%available) n = n + 1
      end do
   end function request_n_missing

   !> Resolve an available output of the caller's rank
   !>
   !> @param[in] self Scientific request
   !> @param[in] name Output name
   !> @param[in] nlead Leading extents of the caller's array
   !> @param[out] error Unknown, mismatched or missing output
   function request_available_index(self, name, nlead, error) result(i)
      !> Scientific request
      class(coupling_request_type), intent(in) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Leading extents of the caller's array
      integer, intent(in) :: nlead
      !> Read error
      type(error_type), allocatable, intent(out) :: error
      !> Output index, zero on error
      integer :: i
      i = self%find(name)
      if (i == 0) then
         call fatal_error(error, trim(self%name())//": unknown output '"//name//"'")
      else if (self%outputs(i)%nlead /= nlead) then
         call fatal_error(error, trim(self%name())//": "//trim(name)//" rank mismatch")
         i = 0
      else if (.not. self%outputs(i)%available) then
         call fatal_error(error, trim(self%name())//": missing required output "//trim(name))
         i = 0
      end if
   end function request_available_index

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

   !> Outputs of the point_potential calculation
   !>
   !> @param[in,out] self Request to declare
   !> @param[out] error Invalid declaration
   subroutine point_potential_declare(self, error)
      !> Request to declare
      class(point_potential_request_type), intent(inout) :: self
      !> Invalid declaration
      type(error_type), allocatable, intent(out) :: error
      call self%add("phi", error)
      if (allocated(error)) return
      call self%add("dphi_dr", [3], error)
      if (allocated(error)) return
      call self%uses("xyz", error)
   end subroutine point_potential_declare

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

   !> Outputs of the gaussian_potential calculation; needs the grid widths
   !>
   !> @param[in,out] self Request to declare
   subroutine gaussian_potential_declare(self, error)
      !> Request to declare
      class(gaussian_potential_request_type), intent(inout) :: self
      !> Invalid declaration
      type(error_type), allocatable, intent(out) :: error
      call self%add("phi", error)
      if (allocated(error)) return
      call self%add("dphi_dr", [3], error)
      if (allocated(error)) return
      call self%add("dphi_dxi", error)
      if (allocated(error)) return
      call self%uses("xyz", error)
      if (allocated(error)) return
      call self%uses("xi0", error)
   end subroutine gaussian_potential_declare

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

   !> Outputs of the gaussian_moments calculation
   !>
   !> @param[in,out] self Request to declare
   !> @param[out] error Invalid declaration
   subroutine gaussian_moment_declare(self, error)
      !> Request to declare
      class(gaussian_moment_request_type), intent(inout) :: self
      !> Invalid declaration
      type(error_type), allocatable, intent(out) :: error
      call self%add("gt", error)
      if (allocated(error)) return
      call self%add("pt", [3], error)
      if (allocated(error)) return
      call self%add("mt", [3, 3], error)
      if (allocated(error)) return
      call self%add("rt", [3], error)
      if (allocated(error)) return
      call self%uses("xyz", error)
   end subroutine gaussian_moment_declare

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
      ! The previous staging ends here, not at arm: a pass that fails in between
      ! must leave no pending list, live handle or matching phase behind
      self%phase = moist_phase_none
      if (allocated(self%pending_index)) deallocate (self%pending_index)
      do i = 1, size(self%requests)
         call self%requests(i)%item%clear_requirements()
         self%requests(i)%item%handle = 0_c_int64_t
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
      type(request_slot_type), allocatable :: grown(:)
      !> Declared copy of the registered value
      class(coupling_request_type), allocatable :: declared
      integer :: i, slot, stat
      logical :: equal
      if (.not. allocated(self%requests)) allocate (self%requests(0))
      if (.not. allocated(self%registrations)) allocate (self%registrations(0))
      if (len_trim(key) > request_name_len .or. len_trim(key) == 0) then
         call fatal_error(error, "Invalid local request name")
         return
      end if
      allocate (declared, source=item, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate scientific request")
         return
      end if
      call declared%ensure_outputs(error)
      if (allocated(error)) return
      ! A registration is a declaration: availability and tokens are the coupling's
      call declared%invalidate()
      declared%handle = 0_c_int64_t
      slot = 0
      do i = 1, size(self%requests)
         equal = self%requests(i)%item%same_inputs(declared)
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
         allocate (grown(slot)%answers(declared%n_outputs()))
         call move_alloc(declared, grown(slot)%item)
         do i = 1, slot - 1
            call move_alloc(self%requests(i)%item, grown(i)%item)
            call move_alloc(self%requests(i)%answers, grown(i)%answers)
         end do
         call move_alloc(grown, self%requests)
      else
         associate (live => self%requests(slot)%item)
            do i = 1, min(live%n_outputs(), declared%n_outputs())
               live%outputs(i)%required = live%outputs(i)%required .or. declared%outputs(i)%required
            end do
            if (any(self%registrations%slot == slot)) then
               live%fields = union_names(live%fields, declared%fields)
            else
               live%fields = declared%fields
            end if
         end associate
      end if
      self%registrations = [self%registrations, registration_type(self%scope, key, slot)]
   end subroutine coupling_register

   !> Sorted union of the cavity grid fields used by active declarations
   !>
   !> @param[in] self Collection to inspect
   function grid_fields(self) result(fields)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Field names
      character(len=grid_name_len), allocatable :: fields(:)
      integer :: i
      allocate (fields(0))
      if (.not. allocated(self%registrations)) return
      do i = 1, size(self%registrations)
         fields = union_names(fields, self%requests(self%registrations(i)%slot)%item%fields)
      end do
      call sort_names(fields)
   end function grid_fields

   !> Whether any active declaration reads one cavity grid field
   !>
   !> @param[in] self Collection to inspect
   !> @param[in] name Cavity field name
   function uses_grid(self, name) result(used)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Cavity field name
      character(len=*), intent(in) :: name
      !> Whether the field is declared
      logical :: used
      used = any(self%grid_fields() == name)
   end function uses_grid

   !> Mark the owner's declared grid fields current; no grid values are copied
   !>
   !> Staleness is the owner's duty: it must invalidate whenever the grid changes
   !>
   !> @param[in,out] self Collection to refresh
   !> @param[in] ngrid Evaluation point count
   !> @param[in] available Grid field names the owner can supply
   !> @param[out] error Missing declared field
   subroutine coupling_snapshot(self, ngrid, available, error)
      !> Collection to refresh
      class(coupling_type), intent(inout) :: self
      !> Evaluation point count
      integer, intent(in) :: ngrid
      !> Grid field names the owner can supply
      character(len=*), intent(in) :: available(:)
      !> Missing declared field
      type(error_type), allocatable, intent(out) :: error
      !> Declared field names
      character(len=grid_name_len), allocatable :: declared(:)
      integer :: i
      self%grid_valid = .false.
      if (ngrid < 0) then
         call fatal_error(error, "Grid size must be nonnegative")
         return
      end if
      call self%compact()
      if (self%ngrid /= ngrid) call self%invalidate()
      self%ngrid = ngrid
      declared = self%grid_fields()
      do i = 1, size(declared)
         if (any(available == declared(i))) cycle
         call fatal_error(error, "Declared grid field "//trim(declared(i))//" is unavailable")
         return
      end do
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
            item%phase = phase
            if (item%n_missing() > 0) self%pending_index = [self%pending_index, i]
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
      i = self%index_of(handle)
      if (i == 0) then
         call fatal_error(error, "Superseded or unknown coupling request handle")
         return
      end if
      allocate (item, source=self%requests(i)%item)
   end subroutine lookup

   !> Slot of a current host token, zero for a superseded or unknown one
   !>
   !> @param[in] self Collection to inspect
   !> @param[in] handle Opaque token from pending()
   function index_of(self, handle) result(i)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Opaque token
      integer(c_int64_t), intent(in) :: handle
      !> Slot index
      integer :: i
      if (handle /= 0_c_int64_t) then
         do i = 1, self%n_requests()
            if (self%requests(i)%item%handle == handle) return
         end do
      end if
      i = 0
   end function index_of

   !> Validate and store one flattened answer (rows, ngrid) into the live request
   !>
   !> A rejected answer leaves the output unanswered, so a wrong retry never
   !> keeps an older value alive
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] handle Opaque token from pending()
   !> @param[in] name Output name
   !> @param[in] lead Caller's leading extents, absent for a pre-flattened array
   !> @param[in] values Flattened answer
   !> @param[out] error Rejected handle, name or answer
   subroutine store(self, handle, name, lead, values, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Opaque token
      integer(c_int64_t), intent(in) :: handle
      !> Output name
      character(len=*), intent(in) :: name
      !> Caller's leading extents
      integer, intent(in), optional :: lead(:)
      !> Flattened answer
      real(wp), intent(in) :: values(:, :)
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      !> Why the answer is rejected, empty when it is stored
      character(len=:), allocatable :: reason
      integer :: i, j
      call self%resolve(handle, name, i, j, error)
      if (allocated(error)) return
      associate (item => self%requests(i)%item)
         associate (out => item%outputs(j))
            out%available = .false.
            reason = ""
            if (present(lead)) then
               if (size(lead) /= out%nlead) then
                  reason = "rank mismatch"
               else if (any(lead /= out%lead(:size(lead)))) then
                  reason = "shape mismatch"
               end if
            end if
            if (len(reason) == 0) then
               if (any(shape(values) /= [item%rows(j), item%ngrid])) then
                  reason = "shape mismatch"
               else if (.not. all(ieee_is_finite(values))) then
                  reason = "must be finite"
               end if
            end if
            if (len(reason) > 0) then
               call fatal_error(error, trim(item%name())//": "//trim(name)//" "//reason)
               return
            end if
            self%requests(i)%answers(j)%values = values
            out%available = .true.
         end associate
      end associate
   end subroutine store

   !> Answer a vector output (ngrid) of a current request
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] handle Opaque token from pending()
   !> @param[in] name Output name
   !> @param[in] values Answer on the grid
   !> @param[out] error Rejected handle, name or answer
   subroutine coupling_answer_1(self, handle, name, values, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Opaque token
      integer(c_int64_t), intent(in) :: handle
      !> Output name
      character(len=*), intent(in) :: name
      !> Answer on the grid
      real(wp), intent(in) :: values(:)
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      !> No leading extents
      integer :: lead(0)
      call self%store(handle, name, lead, reshape(values, [1, size(values)]), error)
   end subroutine coupling_answer_1

   !> Answer a (lead(1), ngrid) output of a current request
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] handle Opaque token from pending()
   !> @param[in] name Output name
   !> @param[in] values Answer on the grid
   !> @param[out] error Rejected handle, name or answer
   subroutine coupling_answer_2(self, handle, name, values, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Opaque token
      integer(c_int64_t), intent(in) :: handle
      !> Output name
      character(len=*), intent(in) :: name
      !> Answer on the grid
      real(wp), intent(in) :: values(:, :)
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      call self%store(handle, name, [size(values, 1)], values, error)
   end subroutine coupling_answer_2

   !> Answer a (lead(1), lead(2), ngrid) output of a current request
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] handle Opaque token from pending()
   !> @param[in] name Output name
   !> @param[in] values Answer on the grid
   !> @param[out] error Rejected handle, name or answer
   subroutine coupling_answer_3(self, handle, name, values, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Opaque token
      integer(c_int64_t), intent(in) :: handle
      !> Output name
      character(len=*), intent(in) :: name
      !> Answer on the grid
      real(wp), intent(in) :: values(:, :, :)
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      call self%store(handle, name, [size(values, 1), size(values, 2)], &
         & reshape(values, [size(values, 1)*size(values, 2), size(values, 3)]), error)
   end subroutine coupling_answer_3

   !> Answer one output already flattened to (rows, ngrid), as C hosts supply it
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] handle Opaque token from pending()
   !> @param[in] name Output name
   !> @param[in] values Flattened answer
   !> @param[out] error Rejected handle, name or answer
   subroutine coupling_answer_flat(self, handle, name, values, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Opaque token
      integer(c_int64_t), intent(in) :: handle
      !> Output name
      character(len=*), intent(in) :: name
      !> Flattened answer
      real(wp), intent(in) :: values(:, :)
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      call self%store(handle, name, values=values, error=error)
   end subroutine coupling_answer_flat

   !> Mark one output unanswered without reading anything, reporting why
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] handle Opaque token from pending()
   !> @param[in] name Output name
   !> @param[in] message Reason, e.g. "grid size must match exactly"
   !> @param[out] error Rejected handle, name or the reason
   subroutine coupling_reject(self, handle, name, message, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Opaque token
      integer(c_int64_t), intent(in) :: handle
      !> Output name
      character(len=*), intent(in) :: name
      !> Reason
      character(len=*), intent(in) :: message
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      integer :: i, j
      call self%resolve(handle, name, i, j, error)
      if (allocated(error)) return
      associate (item => self%requests(i)%item)
         item%outputs(j)%available = .false.
         call fatal_error(error, trim(item%name())//": "//trim(name)//" "//message)
      end associate
   end subroutine coupling_reject

   !> Resolve a current token and an output name to slot and output indices
   !>
   !> @param[in] self Collection to inspect
   !> @param[in] handle Opaque token from pending()
   !> @param[in] name Output name
   !> @param[out] i Slot index
   !> @param[out] j Output index
   !> @param[out] error Superseded handle or unknown name
   subroutine resolve(self, handle, name, i, j, error)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Opaque token
      integer(c_int64_t), intent(in) :: handle
      !> Output name
      character(len=*), intent(in) :: name
      !> Slot index
      integer, intent(out) :: i
      !> Output index
      integer, intent(out) :: j
      !> Resolution error
      type(error_type), allocatable, intent(out) :: error
      j = 0
      i = self%index_of(handle)
      if (i == 0) then
         call fatal_error(error, "Superseded or unknown coupling request handle")
         return
      end if
      j = self%requests(i)%item%find(name)
      if (j == 0) call fatal_error(error, &
         & trim(self%requests(i)%item%name())//" has no output '"//name//"'")
   end subroutine resolve

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
      !> Missing output names, comma separated
      character(len=:), allocatable :: names
      integer :: i, j
      if (phase < 1 .or. phase > moist_n_phases .or. phase /= self%phase) then
         call fatal_error(error, "Coupling is not staged for this phase")
         return
      end if
      do i = 1, self%n_requests()
         associate (item => self%requests(i)%item)
            if (item%n_missing(phase) == 0) cycle
            names = ""
            do j = 1, item%n_outputs()
               if (.not. item%outputs(j)%required(phase) .or. item%outputs(j)%available) cycle
               if (len(names) > 0) names = names//", "
               names = names//trim(item%outputs(j)%name)
            end do
            call fatal_error(error, trim(item%name())//": missing required outputs: "//names)
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
      call self%validate(error)
      if (allocated(error)) return
      call self%owner%check_mandatory(phase, error)
   end subroutine view_check

   !> Refuse an unbound or expired view
   !>
   !> @param[in] self Component view
   !> @param[out] error Unbound or expired view
   subroutine view_validate(self, error)
      !> Component view
      class(coupling_view_type), intent(in) :: self
      !> Unbound or expired view
      type(error_type), allocatable, intent(out) :: error
      if (.not. associated(self%owner)) then
         call fatal_error(error, "Unbound component coupling view")
      else if (self%epoch /= self%owner%epoch) then
         call fatal_error(error, "Expired component coupling view")
      end if
   end subroutine view_validate

   !> Copy a complete request for inspection within its registration scope
   !>
   !> Use `read` to copy only one output array
   !>
   !> @param[in] self Component view
   !> @param[in] key Local name used in register()
   !> @param[out] item Independent request value
   !> @param[out] error Expired view or undeclared name
   subroutine view_request(self, key, item, error)
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
      item%handle = 0_c_int64_t
   end subroutine view_request

   !> Resolve a local registration without copying its answer arrays
   !>
   !> Completeness is the component's entry check, not a per-read rule
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
      call self%validate(error)
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

   !> Copy one vector output (ngrid) of a registered request
   !>
   !> @param[in] self Component view
   !> @param[in] key Local registration name
   !> @param[in] output Output name
   !> @param[out] values Independent copy
   !> @param[out] error Missing output, expired view or undeclared name
   subroutine view_read_1(self, key, output, values, error)
      !> Component view
      class(coupling_view_type), intent(in) :: self
      !> Local registration name
      character(len=*), intent(in) :: key
      !> Output name
      character(len=*), intent(in) :: output
      !> Independent copy
      real(wp), allocatable, intent(out) :: values(:)
      !> Read error
      type(error_type), allocatable, intent(out) :: error
      integer :: slot, j
      slot = self%find_slot(key, error)
      if (allocated(error)) return
      j = request_available_index(self%owner%requests(slot)%item, output, 0, error)
      if (j == 0) return
      values = self%owner%requests(slot)%answers(j)%values(1, :)
   end subroutine view_read_1

   !> Copy one (lead, ngrid) output of a registered request
   !>
   !> @param[in] self Component view
   !> @param[in] key Local registration name
   !> @param[in] output Output name
   !> @param[out] values Independent copy
   !> @param[out] error Missing output, expired view or undeclared name
   subroutine view_read_2(self, key, output, values, error)
      !> Component view
      class(coupling_view_type), intent(in) :: self
      !> Local registration name
      character(len=*), intent(in) :: key
      !> Output name
      character(len=*), intent(in) :: output
      !> Independent copy
      real(wp), allocatable, intent(out) :: values(:, :)
      !> Read error
      type(error_type), allocatable, intent(out) :: error
      integer :: slot, j
      slot = self%find_slot(key, error)
      if (allocated(error)) return
      j = request_available_index(self%owner%requests(slot)%item, output, 1, error)
      if (j == 0) return
      values = self%owner%requests(slot)%answers(j)%values
   end subroutine view_read_2

   !> Copy one (lead(1), lead(2), ngrid) output of a registered request
   !>
   !> @param[in] self Component view
   !> @param[in] key Local registration name
   !> @param[in] output Output name
   !> @param[out] values Independent copy
   !> @param[out] error Missing output, expired view or undeclared name
   subroutine view_read_3(self, key, output, values, error)
      !> Component view
      class(coupling_view_type), intent(in) :: self
      !> Local registration name
      character(len=*), intent(in) :: key
      !> Output name
      character(len=*), intent(in) :: output
      !> Independent copy
      real(wp), allocatable, intent(out) :: values(:, :, :)
      !> Read error
      type(error_type), allocatable, intent(out) :: error
      integer :: slot, j
      slot = self%find_slot(key, error)
      if (allocated(error)) return
      associate (live => self%owner%requests(slot))
         j = request_available_index(live%item, output, 2, error)
         if (j == 0) return
         associate (lead => live%item%outputs(j)%lead)
            values = reshape(live%answers(j)%values, &
               & [lead(1), lead(2), size(live%answers(j)%values, 2)])
         end associate
      end associate
   end subroutine view_read_3

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
      type(request_slot_type), allocatable :: active(:)
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
         call move_alloc(self%requests(i)%answers, active(mapping(i))%answers)
      end do
      call move_alloc(active, self%requests)
      do i = 1, size(self%registrations)
         self%registrations(i)%slot = mapping(self%registrations(i)%slot)
      end do
      self%registrations = pack(self%registrations, self%registrations%slot > 0)
   end subroutine compact

end module moist_channels_coupling
