!> Host -> MOIST coupling requests and their answers
!>
!> A coupling is the model-owned list of calculations the host performs on the
!> cavity grid, one request per distinct set of inputs
!>
!> - a request is a declaration: named outputs such as `phi` or `dphi_dr(3, ngrid)`;
!>   it holds no data, so copies are cheap
!> - each phase stages its requirements; the host walks them with a cursor,
!>   `do while (coupling%next())`, inspects `coupling%request()` and answers
!>   by output name, `coupling%answer("phi", phi, error)`
!> - the coupling owns the answer arrays; they survive across phases until the
!>   owner invalidates them, and a rejected answer leaves its output missing
!> - components never see the coupling; each reads its own registrations through
!>   a short-lived `coupling_view_type` after one completeness check per call
!> - request kinds are a closed set: `declare` is private, so new kinds are
!>   defined in this module only
module moist_channels_coupling
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   implicit none(type, external)
   private
   public :: coupling_type, coupling_view_type, coupling_registry_type
   public :: coupling_request_type
   public :: point_potential_request_type, gaussian_potential_request_type
   public :: gaussian_moment_request_type
   public :: current_request, answer_flat
   public :: coupling_register, request_require
   public :: coupling_begin_registration, coupling_set_scope, coupling_snapshot
   public :: coupling_arm, coupling_invalidate, coupling_check_mandatory
   public :: coupling_make_view, coupling_close_view
   public :: moist_phase_energy, moist_phase_response, moist_phase_gradient
   public :: request_name_len, output_name_len

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
   !> Leading extents an output may declare before ngrid
   integer, parameter :: max_lead = 2
   !> Diagnostic for an answer or query outside a `next()` window
   character(len=*), parameter :: no_current_request = &
      & "No current coupling request - call next() first"

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
   !> output name, so a new kind adds no per-output code. A host sees a copy
   !> through `coupling%request()`: its kind, `name()` and `is_missing(name)`
   type, abstract :: coupling_request_type
      private
      !> Named outputs, filled by `declare` on first use
      type(output_slot_type), allocatable :: outputs(:)
      !> Phase staged by the latest arm
      integer :: phase = moist_phase_none
   contains
      !> Scientific name, e.g. "gaussian_potential"
      procedure(request_name), deferred :: name
      !> Whether an output is required by the staged phase and still unanswered
      procedure :: is_missing => request_is_missing
      procedure(request_declare), deferred, private :: declare
      procedure, private :: ensure_outputs => request_ensure_outputs
      procedure, private :: request_add_vector, request_add_array
      generic, private :: add => request_add_vector, request_add_array
      procedure, private :: n_outputs => request_n_outputs
      procedure, private :: find => request_find
      procedure, private :: clear_requirements => request_clear_requirements
      procedure, private :: is_required => request_is_required
      procedure, private :: is_available => request_is_available
      procedure, private :: n_missing => request_n_missing
      procedure, private :: rows => request_rows
      procedure, private :: invalidate => request_invalidate
      procedure, private :: same_inputs => request_same_inputs
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

      !> Declare the kind's outputs
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
      procedure, private :: declare => point_potential_declare
   end type point_potential_request_type

   !> Raw gaussian_potential calculation: phi, dphi_dr(3, ngrid), dphi_dxi
   type, extends(coupling_request_type) :: gaussian_potential_request_type
   contains
      procedure :: name => gaussian_potential_name
      procedure, private :: declare => gaussian_potential_declare
   end type gaussian_potential_request_type

   !> Raw gaussian_moments calculation: gt, pt(3, ngrid), mt(3, 3, ngrid), rt(3, ngrid)
   type, extends(coupling_request_type) :: gaussian_moment_request_type
      !> Moment exponents in bohr**(-2), supplied by the component
      real(wp), allocatable :: width(:)
   contains
      procedure :: name => gaussian_moment_name
      procedure, private :: declare => gaussian_moment_declare
      procedure, private :: same_inputs => moment_same_inputs
   end type gaussian_moment_request_type

   !> Placeholder `coupling%request()` returns outside a `next()` window; no outputs
   type, extends(coupling_request_type) :: no_request_type
   contains
      procedure :: name => no_request_name
      procedure, private :: declare => no_request_declare
   end type no_request_type

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

   !> Model-owned collection exposed to hosts through borrowed pointers
   !>
   !> The type carries only what a host calls: `next`, `request` and `answer`.
   !> Declaring, staging and reading belong to the model and its components and
   !> are module procedures instead (`coupling_arm`, `coupling_register`, ...),
   !> which the `moist` umbrella does not re-export -- Fortran has no friend
   !> access, so privilege is expressed by what a host can reach, not by a
   !> binding it is asked not to call
   type :: coupling_type
      private
      !> Scientific calculations
      type(request_slot_type), allocatable :: requests(:)
      !> Component-local registrations
      type(registration_type), allocatable :: registrations(:)
      !> Epoch of scoped views
      integer(int64) :: epoch = 0_int64
      !> Scope currently registering requirements
      integer :: scope = 1
      !> Staged phase; none before staging and after invalidation
      integer :: phase = moist_phase_none
      !> Number of evaluation points
      integer :: ngrid = 0
      !> Request of the host walk; zero before the first and after the last
      integer :: cursor = 0
   contains
      !> Advance to the next request with a missing output
      procedure :: next => coupling_next
      !> Copy of the current request
      procedure :: request => coupling_request
      procedure, private :: coupling_answer_1, coupling_answer_2, coupling_answer_3
      !> Answer one output of the current request
      generic :: answer => coupling_answer_1, coupling_answer_2, coupling_answer_3
      procedure, private :: n_requests
      procedure, private :: resolve
      procedure, private :: store
      procedure, private :: reject
      procedure, private :: compact
   end type coupling_type

   !> Declare one output required in one phase, used before `coupling_register`
   interface request_require
      module procedure :: request_require_always
      module procedure :: request_require_when
   end interface request_require

   !> Short-lived read-only component view; never exposes pointers to requests
   type :: coupling_view_type
      private
      !> Borrowed collection, valid only during the component call
      type(coupling_type), pointer :: owner => null()
      !> Component-local registration scope
      integer :: scope = 0
      !> Epoch at construction, checked on every read
      integer(int64) :: epoch = -1_int64
      !> Current phase
      integer, public :: phase = moist_phase_none
   contains
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

   !* ================================================================================= *!
   !*                                 Request declarations                              *!
   !* ================================================================================= *!

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
      if (allocated(self%outputs)) return
      call self%declare(error)
      if (allocated(error)) then
         if (allocated(self%outputs)) deallocate (self%outputs)
      end if
   end subroutine request_ensure_outputs

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
   !> @param[in,out] request Request requirements
   !> @param[in] phase Phase index
   !> @param[in] name Output name
   !> @param[out] error Invalid phase, unknown output or failed declaration
   subroutine request_require_always(request, phase, name, error)
      !> Request requirements
      class(coupling_request_type), intent(inout) :: request
      !> Phase index
      integer, intent(in) :: phase
      !> Output name
      character(len=*), intent(in) :: name
      !> Invalid requirement
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      if (phase < 1 .or. phase > moist_n_phases) then
         call fatal_error(error, trim(request%name())//": invalid coupling phase")
         return
      end if
      call request%ensure_outputs(error)
      if (allocated(error)) return
      i = request%find(name)
      if (i == 0) then
         call fatal_error(error, trim(request%name())//" has no output '"//name//"'")
         return
      end if
      request%outputs(i)%required(phase) = .true.
   end subroutine request_require_always

   !> Declare one output required in one phase under a condition
   !>
   !> @param[in,out] request Request requirements
   !> @param[in] phase Phase index
   !> @param[in] name Output name
   !> @param[in] when Whether the requirement applies
   !> @param[out] error Invalid phase, unknown output or failed declaration
   subroutine request_require_when(request, phase, name, when, error)
      !> Request requirements
      class(coupling_request_type), intent(inout) :: request
      !> Phase index
      integer, intent(in) :: phase
      !> Output name
      character(len=*), intent(in) :: name
      !> Whether the requirement applies
      logical, intent(in) :: when
      !> Invalid requirement
      type(error_type), allocatable, intent(out) :: error
      if (when) call request_require_always(request, phase, name, error)
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
   !> State of the copy, taken when `request()` was called; an unknown name
   !> is not missing
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

   !* ================================================================================= *!
   !*                                   Request kinds                                   *!
   !* ================================================================================= *!

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

   !> Outputs of the gaussian_potential calculation
   !>
   !> @param[in,out] self Request to declare
   !> @param[out] error Invalid declaration
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
   end subroutine gaussian_moment_declare

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

   !> Name of the placeholder returned outside a `next()` window
   !>
   !> @param[in] self Placeholder to describe
   function no_request_name(self) result(name)
      !> Placeholder to describe
      class(no_request_type), intent(in) :: self
      !> Diagnostic name
      character(len=request_name_len) :: name
      name = "no_current_request"
   end function no_request_name

   !> The placeholder stands for no calculation, so it cannot be declared
   !>
   !> @param[in,out] self Placeholder to declare
   !> @param[out] error Always set
   subroutine no_request_declare(self, error)
      !> Placeholder to declare
      class(no_request_type), intent(inout) :: self
      !> Always set
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "The placeholder request '"//trim(self%name())//"' cannot be declared")
   end subroutine no_request_declare

   !* ================================================================================= *!
   !*                              Declaration and staging                              *!
   !* ================================================================================= *!

   !> Start a fresh declaration pass while retaining the stored answers
   !>
   !> @param[in,out] coupling Collection to declare
   subroutine coupling_begin_registration(coupling)
      !> Collection to declare
      class(coupling_type), intent(inout) :: coupling
      integer :: i
      if (allocated(coupling%registrations)) deallocate (coupling%registrations)
      allocate (coupling%registrations(0))
      if (.not. allocated(coupling%requests)) allocate (coupling%requests(0))
      ! The previous staging ends here, not at arm: a pass that fails in between
      ! must leave no current request or matching phase behind
      coupling%phase = moist_phase_none
      coupling%cursor = 0
      do i = 1, size(coupling%requests)
         call coupling%requests(i)%item%clear_requirements()
      end do
   end subroutine coupling_begin_registration

   !> Select the component whose declarations follow
   !>
   !> @param[in,out] coupling Collection to declare
   !> @param[in] scope Component index, zero for the cavity
   subroutine coupling_set_scope(coupling, scope)
      !> Collection to declare
      class(coupling_type), intent(inout) :: coupling
      !> Component index
      integer, intent(in) :: scope
      coupling%scope = scope
   end subroutine coupling_set_scope

   !> Register by local name, sharing calculations only when their inputs match
   !>
   !> @param[in,out] coupling Collection receiving the declaration
   !> @param[in] key Component-local scientific name
   !> @param[in] item Calculation and output requirements
   !> @param[out] error Registration error
   subroutine coupling_register(coupling, key, item, error)
      !> Collection receiving the declaration
      class(coupling_type), intent(inout) :: coupling
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
      if (.not. allocated(coupling%requests)) allocate (coupling%requests(0))
      if (.not. allocated(coupling%registrations)) allocate (coupling%registrations(0))
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
      ! A registration is a declaration: availability is the coupling's
      call declared%invalidate()
      slot = 0
      do i = 1, size(coupling%requests)
         if (coupling%requests(i)%item%same_inputs(declared)) then
            slot = i
            exit
         end if
      end do
      if (slot == 0) then
         slot = size(coupling%requests) + 1
         allocate (grown(slot), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot allocate request registry")
            return
         end if
         allocate (grown(slot)%answers(declared%n_outputs()))
         call move_alloc(declared, grown(slot)%item)
         do i = 1, slot - 1
            call move_alloc(coupling%requests(i)%item, grown(i)%item)
            call move_alloc(coupling%requests(i)%answers, grown(i)%answers)
         end do
         call move_alloc(grown, coupling%requests)
      else
         associate (live => coupling%requests(slot)%item)
            do i = 1, min(live%n_outputs(), declared%n_outputs())
               live%outputs(i)%required = live%outputs(i)%required .or. declared%outputs(i)%required
            end do
         end associate
      end if
      coupling%registrations = [coupling%registrations, registration_type(coupling%scope, key, slot)]
   end subroutine coupling_register

   !> Record the grid size after a declaration pass; a new size drops every answer
   !>
   !> Staleness at an unchanged size is the owner's duty: a model update
   !> invalidates every coupling it minted
   !>
   !> @param[in,out] coupling Collection to refresh
   !> @param[in] ngrid Evaluation point count
   subroutine coupling_snapshot(coupling, ngrid)
      !> Collection to refresh
      class(coupling_type), intent(inout) :: coupling
      !> Evaluation point count
      integer, intent(in) :: ngrid
      call coupling%compact()
      if (coupling%ngrid /= ngrid) call coupling_invalidate(coupling)
      coupling%ngrid = ngrid
   end subroutine coupling_snapshot

   !> Invalidate every answer, the staging and any outstanding scoped view
   !>
   !> @param[in,out] coupling Collection to invalidate
   subroutine coupling_invalidate(coupling)
      !> Collection to invalidate
      class(coupling_type), intent(inout) :: coupling
      integer :: i
      coupling%epoch = coupling%epoch + 1_int64
      coupling%phase = moist_phase_none
      coupling%cursor = 0
      do i = 1, coupling%n_requests()
         call coupling%requests(i)%item%invalidate()
      end do
   end subroutine coupling_invalidate

   !> Stage the phase requirements and start a new host walk
   !>
   !> @param[in,out] coupling Collection to stage
   !> @param[in] phase Requested phase
   !> @param[out] error Invalid phase
   subroutine coupling_arm(coupling, phase, error)
      !> Collection to stage
      class(coupling_type), intent(inout) :: coupling
      !> Requested phase
      integer, intent(in) :: phase
      !> Invalid phase
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      if (phase < 1 .or. phase > moist_n_phases) then
         call fatal_error(error, "Invalid coupling phase")
         return
      end if
      coupling%phase = phase
      coupling%cursor = 0
      coupling%epoch = coupling%epoch + 1_int64
      do i = 1, coupling%n_requests()
         coupling%requests(i)%item%phase = phase
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

   !* ================================================================================= *!
   !*                                     Host walk                                     *!
   !* ================================================================================= *!

   !> Advance to the next request with a missing output of the staged phase
   !>
   !> - each request is visited at most once per pass
   !> - false ends the pass and rewinds, so the next call starts a new one
   !> - a pass left early resumes here; staging resets the cursor
   !> - moves the cursor, so it must stand alone in a loop condition
   !>
   !> @param[in,out] self Collection to walk
   function coupling_next(self) result(more)
      !> Collection to walk
      class(coupling_type), intent(inout) :: self
      !> Whether a request is now current
      logical :: more
      integer :: i
      more = .false.
      if (self%phase /= moist_phase_none) then
         do i = self%cursor + 1, self%n_requests()
            if (self%requests(i)%item%n_missing(self%phase) == 0) cycle
            self%cursor = i
            more = .true.
            return
         end do
      end if
      self%cursor = 0
   end function coupling_next

   !> Copy of the current request; answer it with `answer`
   !>
   !> The copy carries the requirement state at this call, so `is_missing`
   !> tells a host what to compute; outside a `next()` window it is a
   !> placeholder named "no_current_request" with no outputs
   !>
   !> @param[in] self Collection to inspect
   function coupling_request(self) result(item)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Independent request value
      class(coupling_request_type), allocatable :: item
      if (self%cursor < 1 .or. self%cursor > self%n_requests()) then
         allocate (no_request_type :: item)
         return
      end if
      allocate (item, source=self%requests(self%cursor)%item)
   end function coupling_request

   !> Copy the current request, or report that none is current
   !>
   !> The checked form of `coupling%request()` for the C layer
   !>
   !> @param[in] coupling Collection to inspect
   !> @param[out] item Independent request value
   !> @param[out] error No current request
   subroutine current_request(coupling, item, error)
      !> Collection to inspect
      class(coupling_type), intent(in) :: coupling
      !> Independent request value
      class(coupling_request_type), allocatable, intent(out) :: item
      !> No current request
      type(error_type), allocatable, intent(out) :: error
      if (coupling%cursor < 1 .or. coupling%cursor > coupling%n_requests()) then
         call fatal_error(error, no_current_request)
         return
      end if
      allocate (item, source=coupling%requests(coupling%cursor)%item)
   end subroutine current_request

   !> Resolve an output name of the current request to slot and output indices
   !>
   !> @param[in] self Collection to inspect
   !> @param[in] name Output name
   !> @param[out] i Slot index
   !> @param[out] j Output index
   !> @param[out] error No current request or unknown name
   subroutine resolve(self, name, i, j, error)
      !> Collection to inspect
      class(coupling_type), intent(in) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Slot index
      integer, intent(out) :: i
      !> Output index
      integer, intent(out) :: j
      !> Resolution error
      type(error_type), allocatable, intent(out) :: error
      j = 0
      i = self%cursor
      if (i < 1 .or. i > self%n_requests()) then
         i = 0
         call fatal_error(error, no_current_request)
         return
      end if
      j = self%requests(i)%item%find(name)
      if (j == 0) call fatal_error(error, &
         & trim(self%requests(i)%item%name())//" has no output '"//name//"'")
   end subroutine resolve

   !> Mark one output unanswered and report why its answer was refused
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] i Slot index
   !> @param[in] j Output index
   !> @param[in] reason Why the answer is refused, e.g. "shape mismatch"
   !> @param[out] error The refusal
   subroutine reject(self, i, j, reason, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Slot index
      integer, intent(in) :: i
      !> Output index
      integer, intent(in) :: j
      !> Why the answer is refused
      character(len=*), intent(in) :: reason
      !> The refusal
      type(error_type), allocatable, intent(out) :: error
      associate (item => self%requests(i)%item)
         item%outputs(j)%available = .false.
         call fatal_error(error, trim(item%name())//": "//trim(item%outputs(j)%name)//" "//reason)
      end associate
   end subroutine reject

   !> Validate and store one flattened answer (rows, ngrid)
   !>
   !> A rejected answer leaves the output unanswered, so a wrong retry never
   !> keeps an older value alive
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] i Slot index
   !> @param[in] j Output index
   !> @param[in] values Flattened answer
   !> @param[out] error Rejected answer
   !> @param[in] lead Caller's leading extents, absent for a pre-flattened array
   subroutine store(self, i, j, values, error, lead)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Slot index
      integer, intent(in) :: i
      !> Output index
      integer, intent(in) :: j
      !> Flattened answer
      real(wp), intent(in) :: values(:, :)
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      !> Caller's leading extents
      integer, intent(in), optional :: lead(:)
      !> Why the answer is rejected, empty when it is stored
      character(len=:), allocatable :: reason
      reason = ""
      associate (item => self%requests(i)%item)
         if (present(lead)) then
            if (size(lead) /= item%outputs(j)%nlead) then
               reason = "rank mismatch"
            else if (any(lead /= item%outputs(j)%lead(:size(lead)))) then
               reason = "shape mismatch"
            end if
         end if
         if (len(reason) == 0) then
            if (any(shape(values) /= [item%rows(j), self%ngrid])) then
               reason = "shape mismatch"
            else if (.not. all(ieee_is_finite(values))) then
               reason = "must be finite"
            end if
         end if
      end associate
      if (len(reason) > 0) then
         call self%reject(i, j, reason, error)
         return
      end if
      self%requests(i)%answers(j)%values = values
      self%requests(i)%item%outputs(j)%available = .true.
   end subroutine store

   !> Answer a vector output (ngrid) of the current request
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] name Output name
   !> @param[in] values Answer on the grid
   !> @param[out] error No current request or rejected answer
   subroutine coupling_answer_1(self, name, values, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Answer on the grid
      real(wp), intent(in) :: values(:)
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      !> No leading extents
      integer :: lead(0)
      integer :: i, j
      call self%resolve(name, i, j, error)
      if (allocated(error)) return
      call self%store(i, j, reshape(values, [1, size(values)]), error, lead)
   end subroutine coupling_answer_1

   !> Answer a (lead(1), ngrid) output of the current request
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] name Output name
   !> @param[in] values Answer on the grid
   !> @param[out] error No current request or rejected answer
   subroutine coupling_answer_2(self, name, values, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Answer on the grid
      real(wp), intent(in) :: values(:, :)
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      integer :: i, j
      call self%resolve(name, i, j, error)
      if (allocated(error)) return
      call self%store(i, j, values, error, [size(values, 1)])
   end subroutine coupling_answer_2

   !> Answer a (lead(1), lead(2), ngrid) output of the current request
   !>
   !> @param[in,out] self Collection to answer
   !> @param[in] name Output name
   !> @param[in] values Answer on the grid
   !> @param[out] error No current request or rejected answer
   subroutine coupling_answer_3(self, name, values, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: self
      !> Output name
      character(len=*), intent(in) :: name
      !> Answer on the grid
      real(wp), intent(in) :: values(:, :, :)
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      integer :: i, j
      call self%resolve(name, i, j, error)
      if (allocated(error)) return
      call self%store(i, j, reshape(values, [size(values, 1)*size(values, 2), size(values, 3)]), &
         & error, [size(values, 1), size(values, 2)])
   end subroutine coupling_answer_3

   !> Answer one output of the current request from a flat C-order buffer
   !>
   !> - `values` holds rows * ngrid numbers in the declared (rows, ngrid) order,
   !>   with ngrid the grid size of the coupling; the caller's buffer is trusted
   !>   to be that long, as a pointer cannot say otherwise
   !> - assumed size on purpose: the rows of an output are private here, so
   !>   the C layer cannot build the (rows, ngrid) view itself
   !>
   !> @param[in,out] coupling Collection to answer
   !> @param[in] name Output name
   !> @param[in] values Flattened answer, rows * ngrid long
   !> @param[out] error No current request or rejected answer
   subroutine answer_flat(coupling, name, values, error)
      !> Collection to answer
      class(coupling_type), intent(inout) :: coupling
      !> Output name
      character(len=*), intent(in) :: name
      !> Flattened answer
      real(wp), intent(in) :: values(*)
      !> Rejected answer
      type(error_type), allocatable, intent(out) :: error
      integer :: i, j, rows
      call coupling%resolve(name, i, j, error)
      if (allocated(error)) return
      rows = coupling%requests(i)%item%rows(j)
      call coupling%store(i, j, reshape(values(:rows*coupling%ngrid), [rows, coupling%ngrid]), &
         & error)
   end subroutine answer_flat

   !> Check the staging and the required outputs before a phase is consumed
   !>
   !> The staging comes first and is named: a neighbouring phase can carry
   !> every answer by accident, e.g. the response answers are still fresh when
   !> the gradient phase is read
   !>
   !> @param[in] coupling Collection to validate
   !> @param[in] phase Phase to consume
   !> @param[out] error Wrong staging or missing scientific output
   subroutine coupling_check_mandatory(coupling, phase, error)
      !> Collection to validate
      class(coupling_type), intent(in) :: coupling
      !> Phase to consume
      integer, intent(in) :: phase
      !> Wrong staging or missing scientific output
      type(error_type), allocatable, intent(out) :: error
      !> Missing output names, comma separated
      character(len=:), allocatable :: names
      integer :: i, j
      if (phase < 1 .or. phase > moist_n_phases) then
         call fatal_error(error, "Invalid coupling phase")
         return
      end if
      if (coupling%phase /= phase) then
         if (coupling%phase == moist_phase_none) then
            call fatal_error(error, "get_"//trim(phase_name(phase))// &
               & " requires a coupling staged by prepare_"//trim(phase_name(phase))// &
               & "; this coupling is not staged")
         else
            call fatal_error(error, "get_"//trim(phase_name(phase))// &
               & " requires a coupling staged by prepare_"//trim(phase_name(phase))// &
               & "; this coupling is staged for the "//trim(phase_name(coupling%phase))//" phase")
         end if
         return
      end if
      do i = 1, coupling%n_requests()
         associate (item => coupling%requests(i)%item)
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
   end subroutine coupling_check_mandatory

   !> Name of a phase for diagnostics, e.g. "gradient"
   !>
   !> Fixed length, so that no deferred-length temporary is created for the
   !> result; callers `trim` it
   !>
   !> @param[in] phase Phase index
   pure function phase_name(phase) result(name)
      !> Phase index
      integer, intent(in) :: phase
      !> Phase name, blank padded
      character(len=8) :: name
      select case (phase)
      case (moist_phase_energy)
         name = "energy"
      case (moist_phase_response)
         name = "response"
      case (moist_phase_gradient)
         name = "gradient"
      case default
         name = "none"
      end select
   end function phase_name

   !* ================================================================================= *!
   !*                                  Component views                                  *!
   !* ================================================================================= *!

   !> Create a read-only component scope for one call
   !>
   !> @param[in,out] coupling Collection owning the answers
   !> @param[in] scope Component registration index
   !> @param[out] view Scoped read interface
   subroutine coupling_make_view(coupling, scope, view)
      !> Collection owning the answers
      class(coupling_type), intent(inout), target :: coupling
      !> Component registration index
      integer, intent(in) :: scope
      !> Scoped read interface
      type(coupling_view_type), intent(out) :: view
      coupling%epoch = coupling%epoch + 1_int64
      view%owner => coupling
      view%scope = scope
      view%epoch = coupling%epoch
      view%phase = coupling%phase
   end subroutine coupling_make_view

   !> End the component call and invalidate any captured value copy of its view
   !>
   !> @param[in,out] coupling Collection owning the answers
   subroutine coupling_close_view(coupling)
      !> Collection owning the answers
      class(coupling_type), intent(inout) :: coupling
      coupling%epoch = coupling%epoch + 1_int64
   end subroutine coupling_close_view

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
      call coupling_check_mandatory(self%owner, phase, error)
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

   !* ================================================================================= *!
   !*                                   Model registry                                  *!
   !* ================================================================================= *!

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
         call coupling_invalidate(node%item)
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

end module moist_channels_coupling
