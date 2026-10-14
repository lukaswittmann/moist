!> Abstract cavity: the surface grid every discretization shares
!>
!> Update and build, the same for every cavity type
!>
!> - `update(mol, error)` takes the structure and does everything that needs
!>   no host state; a geometric cavity then holds its surface, a cavity whose
!>   level set requires build-phase answers (isodensity) defers it
!> - `prepare_build(coupling, error)` clears the coupling, declares the
!>   build-phase requests and stages them; the host walks and answers them
!> - `build(coupling, error)` reads those answers, completes a deferred
!>   surface and invalidates every other coupling the cavity minted; on a
!>   geometric cavity the walk is empty and nothing is rebuilt, a deferring
!>   cavity re-reads its answers and re-projects on every build
!> - every `update` calls `update_cavity_base` first, every `build` calls
!>   `build_cavity_base` first; each sets `updated` and `built` once its own
!>   work succeeded
!> - `updated` and `built` record the latest successful update and build;
!>   a deferring cavity is not built from `prepare_build` until a `build`
!>   succeeds, so a skipped or failed build leaves nothing readable
module moist_cavity_type
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use mctc_io, only: structure_type
   use moist_radius_type, only: radius_type
   use moist_context, only: moist_context_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_channels_response, only: response_type
   use moist_channels_coupling, only: coupling_type, coupling_registry_type, &
      & moist_phase_build, coupling_begin_registration, coupling_set_scope, coupling_snapshot, &
      & coupling_arm, coupling_invalidate, coupling_check_mandatory
   use moist_channels_fields, only: field_query_type
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter

   implicit none(type, external)
   private

   public :: cavity_type, list_cavity_fields_base, add_volume_origin_gradient
   public :: update_cavity_base, build_cavity_base, cavity_declare_pass
   public :: write_cavity_xyz_debug, write_cavity_csv_debug, write_cavity_pqr_debug

   !> Abstract base type containing minimal cavity/surface information
   !>
   !> Cavities within moist are per default discretized using Gaussians
   type, abstract :: cavity_type
      !> Borrowed run context (verbosity/debug/timer); owned by the caller,
      !> never allocated or freed by the cavity
      !>
      !> - optional at construction; a model points a cavity copy without one
      !>   at its own context
      type(moist_context_type), pointer :: ctx => null()

      !> Short name of the discretization, e.g. `DROP`; set by the constructor
      character(len=:), allocatable :: label

      !> Sphere radii, bohr (nat)
      real(wp), allocatable :: radii(:)
      !> Radii model used to update cached radii
      class(radius_type), allocatable :: radius_model

      !> Number of atomic spheres
      integer :: nsph = 0

      !> Area per sphere (nsph)
      real(wp), allocatable :: asph(:)

      !> Number of cavity points
      integer :: ngrid = 0

      !> Cartesian coordinates of atomic sphere centers (3, nsph)
      real(wp), allocatable :: sphxyz(:, :)

      !> Owner of each grid point (ngrid)
      integer, allocatable :: owner(:)

      !> Cartesian coords of points (3,ngrid)
      real(wp), allocatable :: xyz(:, :)
      !> Point area, bohr**2 (ngrid)
      real(wp), allocatable :: a(:)
      !> Gaussian switching factor of each surface point (ngrid)
      real(wp), allocatable :: f(:)
      !> Gaussian width of each surface point (ngrid)
      real(wp), allocatable :: xi0(:)
      !> Outward unit normal of each surface point (3, ngrid)
      real(wp), allocatable :: normal0(:, :)
      !> Point volume element (divergence theorem), bohr**3 (ngrid)
      real(wp), allocatable :: v(:)

      !> Total surface area, bohr**2
      real(wp), allocatable :: total_area
      !> Total cavity volume, bohr**3
      real(wp), allocatable :: total_volume

      !> Nuclear derivatives of surface Gaussian widths (3, nsph, ngrid)
      real(wp), allocatable :: xi1_rA(:, :, :)
      !> Nuclear derivatives of Gaussian switching factors (3, nsph, ngrid)
      real(wp), allocatable :: f1_rA(:, :, :)
      !> Nuclear derivatives of surface positions (3, 3, nsph, ngrid)
      real(wp), allocatable :: xyz1_rA(:, :, :, :)
      !> Nuclear derivatives of surface point volumes (3, nsph, ngrid)
      !>
      !> - contracting over the grid gives the total-volume nuclear gradient
      !> - include the motion of the volume origin, so every atom is reached
      real(wp), allocatable :: v1_rA(:, :, :)

      !> Molecular structure of the latest successful `update`
      type(structure_type), allocatable :: mol
      !> Whether the latest update completed successfully
      logical :: updated = .false.
      !> Whether the surface is complete: built by the latest update, or by the
      !> latest `build` when the update defers it
      logical :: built = .false.
      !> Registry of the couplings minted by this cavity
      type(coupling_registry_type), private :: couplings
   contains

      !* ------------------------------ Update and build ------------------------------ *!

      !> Take a new structure and do everything that needs no host state;
      !> calls `update_cavity_base` first
      procedure(update_cavity), deferred :: update
      !> Whether the surface waits for build-phase answers after `update`
      procedure :: defers_build => cavity_defers_build_default
      !> Stage the build phase: clear and declare the coupling
      procedure, non_overridable :: prepare_build => cavity_prepare_build
      !> Read the build-phase answers of a staged coupling; a deferring cavity
      !> overrides it, calling `build_cavity_base` first
      procedure :: build => build_cavity_base

      !* ------------------------------- Surface state -------------------------------- *!

      !> Drop the surface: `ngrid` zero, every surface array of the base
      !> deallocated, not built
      procedure :: clear_surface => clear_cavity_surface_base
      !> Whether the latest update completed
      procedure, non_overridable :: is_updated => cavity_is_updated
      !> Whether the surface is built
      procedure, non_overridable :: is_built => cavity_is_built
      !> Require a successfully updated cavity
      procedure, non_overridable :: require_updated => cavity_require_updated
      !> Require a built surface
      procedure, non_overridable :: require_built => cavity_require_built

      !* ------------------------------- Host couplings ------------------------------- *!

      !> Declare the cavity's own host requests (none by default)
      procedure :: declare_coupling => declare_cavity_coupling_default
      !> Extents a coupling of this cavity records: grid and atoms
      procedure, non_overridable :: snapshot_counts => cavity_snapshot_counts
      !> Build a cavity-owned coupling and declare the cavity's requests
      procedure, non_overridable :: new_coupling => cavity_new_coupling
      !> Release a coupling minted by this cavity
      procedure, non_overridable :: release_coupling => cavity_release_coupling
      !> Release every coupling minted by this cavity
      procedure, non_overridable :: clear_couplings => cavity_clear_couplings
      !> Whether this cavity minted a coupling
      procedure, non_overridable :: owns_coupling => cavity_owns_coupling

      !* -------------------------- Derivatives and adjoints -------------------------- *!

      procedure(get_cavity_gradient), deferred :: get_gradient
      !> Map accumulated surface-observable adjoints to host response channels
      procedure :: get_surface_response => get_cavity_surface_response_default
      !> Contract accumulated surface-observable adjoints into the nuclear gradient:
      !> the cavity's surface channels, then the volume-origin channel
      procedure, non_overridable :: get_surface_gradient => get_cavity_surface_gradient
      !> Contract the per-point surface channels into the nuclear gradient
      procedure :: contract_surface_gradient => contract_cavity_surface_gradient_default
      !> Whether the surface geometry responds to the host field (DROP), so that
      !> host surface weights are consumed in the response phase
      procedure :: has_field_dependent_geometry => cavity_field_dependent_default
      !> Origin of the divergence-theorem volume, the centroid of the sphere centers
      procedure, non_overridable :: get_volume_origin => get_cavity_volume_origin

      !* ----------------------- Readable results and printing ------------------------ *!

      !> Declare the readable results this cavity currently holds
      procedure :: list_fields => list_cavity_fields_base
      !> Print basic cavity information
      procedure :: print => print_cavity_info
      !> Solvent-independent settings; null for a cavity without any
      procedure :: parameters => cavity_parameters_default
      !> Print the cavity kind and its registered settings
      procedure :: print_parameters => print_cavity_parameters
      !> Print the radius model of the current structure; called by every `update`
      procedure, non_overridable :: print_radii => print_cavity_radii
      !> Write grid to XYZ file for visualization
      procedure :: write_xyz_debug => write_cavity_xyz_debug
      !> Write grid to CSV file for visualization
      procedure :: write_csv_debug => write_cavity_csv_debug
      !> Write grid to PQR file for visualization
      procedure :: write_pqr_debug => write_cavity_pqr_debug

      !* ---------------------------------- Context ----------------------------------- *!

      !> Fail without a run context; called first by every `update`
      procedure, non_overridable :: require_context => require_cavity_context

   end type cavity_type

   ! Abstract interfaces for deferred procedures
   abstract interface

      !> Take a new structure and do everything that needs no host state
      !>
      !> - calls `update_cavity_base` first
      !> - the whole surface unless `defers_build`; sets `updated` and `built`
      !>   on success
      !>
      !> @param[in,out] self Cavity instance
      !> @param[in] mol Molecular structure
      !> @param[out] error Error handling
      subroutine update_cavity(self, mol, error)
         import :: cavity_type, structure_type, error_type
         implicit none(type, external)
         !> Cavity instance
         class(cavity_type), intent(inout) :: self
         !> Molecular structure
         type(structure_type), intent(in) :: mol
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_cavity

      subroutine get_cavity_gradient(self, error)
         import :: cavity_type, error_type
         implicit none(type, external)
         class(cavity_type), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_cavity_gradient

   end interface

contains

   !* ================================================================================= *!
   !*                                 Update and build                                  *!
   !* ================================================================================= *!

   !> Shared first step of every `update`
   !>
   !> - invalidates every coupling the cavity minted; neither updated nor
   !>   built until the caller's `update` succeeds
   !> - copies `mol`, which may be `self%mol` itself
   !>
   !> @param[in,out] self  Cavity instance
   !> @param[in]     mol   Molecular structure
   !> @param[out]    error Error handling
   subroutine update_cavity_base(self, mol, error)
      !> Cavity instance
      class(cavity_type), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Copy of `mol`, which may be `self%mol` itself
      type(structure_type) :: tmp

      call self%couplings%invalidate()
      self%updated = .false.
      self%built = .false.
      tmp = mol
      self%mol = tmp

   end subroutine update_cavity_base

   !> Default: the surface needs no build-phase answer
   !>
   !> @param[in] self Cavity instance
   function cavity_defers_build_default(self) result(defers)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Whether `update` leaves the surface to `build`
      logical :: defers

      defers = .false.

   end function cavity_defers_build_default

   !> Stage the build phase of a coupling
   !>
   !> - clears every answer of the coupling; the atom count comes from the
   !>   structure of the latest `update`
   !> - a deferring cavity is not built until the next `build` succeeds
   !>
   !> @param[in,out] self     Updated cavity
   !> @param[in,out] coupling Coupling to stage
   !> @param[out]    error    Error handling
   subroutine cavity_prepare_build(self, coupling, error)
      !> Cavity instance
      class(cavity_type), intent(inout) :: self
      !> Coupling to stage
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%require_updated(error)
      if (allocated(error)) return
      if (self%defers_build()) self%built = .false.
      call coupling_invalidate(coupling)
      call cavity_declare_pass(self, coupling, error)
      if (allocated(error)) return
      call coupling_arm(coupling, moist_phase_build, error)

   end subroutine cavity_prepare_build

   !> Read the build-phase answers of `coupling`; shared first step of every
   !> `build`
   !>
   !> - the coupling must be staged by `prepare_build` and fully answered
   !> - its build-phase answers stay; every other coupling the cavity minted
   !>   is invalidated
   !> - not built from entry, so a failed build leaves the cavity unbuilt
   !> - a cavity that `update` built completely is built again and the
   !>   coupling records its extents; one that `defers_build` overrides
   !>   `build`, calls this first and then completes its surface: built, and
   !>   the coupling re-declared and so unstaged, only at the end of the
   !>   override
   !>
   !> @param[in,out] self     Updated cavity
   !> @param[in,out] coupling Coupling staged by `prepare_build`
   !> @param[out]    error    Error handling
   subroutine build_cavity_base(self, coupling, error)
      !> Cavity instance
      class(cavity_type), intent(inout) :: self
      !> Coupling staged by `prepare_build`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%require_updated(error)
      if (allocated(error)) return
      self%built = .false.
      call coupling_check_mandatory(coupling, moist_phase_build, error)
      if (allocated(error)) return
      call self%couplings%invalidate(except=coupling)
      if (self%defers_build()) return
      self%built = .true.
      call cavity_declare_pass(self, coupling, error)

   end subroutine build_cavity_base

   !* ================================================================================= *!
   !*                                   Surface state                                   *!
   !* ================================================================================= *!

   !> Drop the surface: `ngrid` zero, every surface array of the base
   !> deallocated, not built
   !>
   !> Arrays of a concrete cavity may stay; nothing reads them while the
   !> cavity is not built
   !>
   !> @param[in,out] self Cavity instance
   subroutine clear_cavity_surface_base(self)
      !> Cavity instance
      class(cavity_type), intent(inout) :: self

      self%built = .false.
      self%ngrid = 0
      if (allocated(self%asph)) deallocate (self%asph)
      if (allocated(self%owner)) deallocate (self%owner)
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%a)) deallocate (self%a)
      if (allocated(self%f)) deallocate (self%f)
      if (allocated(self%xi0)) deallocate (self%xi0)
      if (allocated(self%normal0)) deallocate (self%normal0)
      if (allocated(self%v)) deallocate (self%v)
      if (allocated(self%total_area)) deallocate (self%total_area)
      if (allocated(self%total_volume)) deallocate (self%total_volume)
      if (allocated(self%xi1_rA)) deallocate (self%xi1_rA)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      if (allocated(self%xyz1_rA)) deallocate (self%xyz1_rA)
      if (allocated(self%v1_rA)) deallocate (self%v1_rA)

   end subroutine clear_cavity_surface_base

   !> Whether the latest update completed
   !>
   !> @param[in] self Cavity instance
   function cavity_is_updated(self) result(updated)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Update status
      logical :: updated

      updated = self%updated

   end function cavity_is_updated

   !> Whether the surface is built
   !>
   !> @param[in] self Cavity instance
   function cavity_is_built(self) result(built)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Build status
      logical :: built

      built = self%built

   end function cavity_is_built

   !> Require a successfully updated cavity
   !>
   !> @param[in]  self  Cavity instance
   !> @param[out] error Error handling
   subroutine cavity_require_updated(self, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%updated) call fatal_error(error, "Cavity must be updated first")

   end subroutine cavity_require_updated

   !> Require a built surface
   !>
   !> @param[in]  self  Cavity instance
   !> @param[out] error Error handling
   subroutine cavity_require_built(self, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%built) call fatal_error(error, "Cavity must be built first")

   end subroutine cavity_require_built

   !* ================================================================================= *!
   !*                                  Host couplings                                   *!
   !* ================================================================================= *!

   !> Default cavity-side declaration: a cavity requests nothing of its own
   !>
   !> @param[in]    self     Cavity instance
   !> @param[inout] coupling Coupling being declared
   !> @param[out]   error    Error handling
   subroutine declare_cavity_coupling_default(self, coupling, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Coupling being declared
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine declare_cavity_coupling_default

   !> Declare the cavity's requests at scope 0, then record the extents
   !>
   !> - `natom` once updated, from the spheres of the structure
   !> - `ngrid` once built, so grid outputs are refused before
   !>
   !> @param[in]    self     Cavity instance
   !> @param[inout] coupling Coupling to declare
   !> @param[out]   error    Error handling
   subroutine cavity_declare_pass(self, coupling, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Coupling to declare
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Counts the snapshot records; unallocated ones are absent
      integer, allocatable :: ngrid, natom

      call coupling_begin_registration(coupling)
      call coupling_set_scope(coupling, 0)
      call self%declare_coupling(coupling, error)
      if (allocated(error)) return
      call self%snapshot_counts(ngrid, natom)
      call coupling_snapshot(coupling, ngrid=ngrid, natom=natom)

   end subroutine cavity_declare_pass

   !> Extents a coupling of this cavity records; an unallocated count is absent
   !>
   !> - `ngrid` once built, so grid outputs are refused before
   !> - `natom` once updated
   !>
   !> @param[in]  self  Cavity instance
   !> @param[out] ngrid Surface point count
   !> @param[out] natom Sphere count
   subroutine cavity_snapshot_counts(self, ngrid, natom)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Surface point count
      integer, allocatable, intent(out) :: ngrid
      !> Sphere count
      integer, allocatable, intent(out) :: natom

      if (self%built) ngrid = self%ngrid
      if (self%updated) natom = self%nsph

   end subroutine cavity_snapshot_counts

   !> Build a cavity-owned coupling, declared but not staged
   !>
   !> The snapshot records the atom count once updated and the point count
   !> once built; outputs over a count it lacks are refused
   !>
   !> @param[in,out] self     Cavity instance
   !> @param[out]    coupling Borrowed coupling, released by `release_coupling`
   !> @param[out]    error    Error handling
   subroutine cavity_new_coupling(self, coupling, error)
      !> Cavity instance
      class(cavity_type), intent(inout), target :: self
      !> Borrowed coupling
      type(coupling_type), pointer, intent(out) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      nullify (coupling)
      call self%couplings%mint(coupling, error)
      if (allocated(error)) return
      call cavity_declare_pass(self, coupling, error)
      if (allocated(error)) call self%couplings%release(coupling)

   end subroutine cavity_new_coupling

   !> Release a coupling minted by this cavity
   !>
   !> @param[in,out] self     Cavity instance
   !> @param[in,out] coupling Coupling to release, nullified
   subroutine cavity_release_coupling(self, coupling)
      !> Cavity instance
      class(cavity_type), intent(inout), target :: self
      !> Coupling to release
      type(coupling_type), pointer, intent(inout) :: coupling

      call self%couplings%release(coupling)

   end subroutine cavity_release_coupling

   !> Release every coupling minted by this cavity, e.g. on a copy
   !>
   !> @param[in,out] self Cavity instance
   subroutine cavity_clear_couplings(self)
      !> Cavity instance
      class(cavity_type), intent(inout) :: self

      call self%couplings%clear()

   end subroutine cavity_clear_couplings

   !> Whether this cavity minted a coupling
   !>
   !> @param[in] self     Cavity instance
   !> @param[in] coupling Coupling to check
   function cavity_owns_coupling(self, coupling) result(owned)
      !> Cavity instance
      class(cavity_type), intent(in), target :: self
      !> Coupling to check
      type(coupling_type), intent(in), target :: coupling
      !> Whether the cavity minted it
      logical :: owned

      owned = self%couplings%owns(coupling)

   end function cavity_owns_coupling

   !* ================================================================================= *!
   !*                             Derivatives and adjoints                              *!
   !* ================================================================================= *!

   !> Default surface-response hook for cavities without field-dependent geometry
   !>
   !> @param[inout] self      Cavity instance, unchanged
   !> @param[in]    acc       Surface-observable adjoints, unused
   !> @param[inout] response  Response accumulator, unchanged
   !> @param[out]   error     Error handling
   subroutine get_cavity_surface_response_default(self, acc, response, error)
      !> Cavity instance
      class(cavity_type), intent(inout) :: self
      !> Surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Response accumulator
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_cavity_surface_response_default

   !> Contract accumulated surface-observable adjoints into the nuclear gradient
   !>
   !> - the cavity's `contract_surface_gradient` maps the per-point channels
   !> - the volume-origin channel is mapped here, the same for every cavity
   !>
   !> @param[in]    self     Cavity instance
   !> @param[in]    acc      Surface-observable adjoints
   !> @param[inout] gradient Nuclear-gradient accumulator (3, nsph)
   !> @param[out]   error    Error handling
   subroutine get_cavity_surface_gradient(self, acc, gradient, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%contract_surface_gradient(acc, gradient, error)
      if (allocated(error)) return
      call add_volume_origin_gradient(acc%w_c, gradient)

   end subroutine get_cavity_surface_gradient

   !> Default reverse-mode surface contraction
   !>
   !> Cavities that do not implement the surface-adjoint contraction must be
   !> reached through the forward path instead; returning silently here would
   !> hand back a zero gradient, so it errors
   !>
   !> @param[in]    self     Cavity instance
   !> @param[in]    acc      Surface-observable adjoints, unused
   !> @param[inout] gradient Nuclear-gradient accumulator, unchanged
   !> @param[out]   error    Error handling
   subroutine contract_cavity_surface_gradient_default(self, acc, gradient, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "This cavity does not provide a reverse-mode surface gradient")

   end subroutine contract_cavity_surface_gradient_default

   !> Default geometry query: the surface does not respond to the host field
   !>
   !> A cavity that overrides `get_surface_response` (DROP) overrides this to
   !> `.true.`; it decides whether host surface weights are declared for the
   !> response phase or only for the gradient
   !>
   !> @param[in] self Cavity instance
   function cavity_field_dependent_default(self) result(field_dependent)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Whether the surface geometry depends on the host field
      logical :: field_dependent

      field_dependent = .false.

   end function cavity_field_dependent_default

   !> Origin `c` of the divergence-theorem volume `V = sum_i a_i ((r_i - c) . n_i)/3`
   !>
   !> - the discrete surface leaves `sum_i a_i n_i` nonzero, so `V` depends on `c`
   !>   and a fixed origin makes it change linearly under rigid translation
   !> - the centroid of the sphere centers moves with the molecule; each center
   !>   enters with weight `1/nsph`, which the volume derivatives must carry
   !>
   !> @param[in] self Cavity instance holding the sphere centers
   !> @return Volume origin, bohr
   pure function get_cavity_volume_origin(self) result(origin)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Volume origin, bohr
      real(wp) :: origin(3)

      origin = sum(self%sphxyz(:, :self%nsph), dim=2)/real(max(1, self%nsph), wp)

   end function get_cavity_volume_origin

   !> Add a volume-origin adjoint to the nuclear gradient
   !>
   !> - the origin of `get_volume_origin` moves with every atom by `dc/dR_A = I/nsph`,
   !>   so `w_c = dE/dc` reaches each atom as `w_c/nsph`
   !> - serves the reverse contraction and the forward `*_rA` arrays alike; the
   !>   map is constant, so a second-order pass reuses it unchanged
   !>
   !> @param[in]    w_c      Adjoint of the volume origin (3)
   !> @param[inout] gradient Nuclear-gradient accumulator (3, nsph)
   pure subroutine add_volume_origin_gradient(w_c, gradient)
      !> Adjoint of the volume origin
      real(wp), intent(in) :: w_c(3)
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)

      !> Atom index
      integer :: iat

      do iat = 1, size(gradient, 2)
         gradient(:, iat) = gradient(:, iat) + w_c/real(size(gradient, 2), wp)
      end do

   end subroutine add_volume_origin_gradient

   !* ================================================================================= *!
   !*                           Readable results and printing                           *!
   !* ================================================================================= *!

   !> Declare the results every discretized cavity carries
   !>
   !> Before the build only the spheres: `nsph`, `radii` and `sphxyz`
   !>
   !> @param[in]    self   Cavity instance
   !> @param[inout] query  Walker collecting or fetching the declarations
   subroutine list_cavity_fields_base(self, query)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Walker collecting or fetching the declarations
      type(field_query_type), intent(inout) :: query

      call query%add_int_value("nsph", "Number of atomic spheres", self%nsph)
      call query%add_real("radii", "Sphere radii, bohr (nsph)", self%radii)
      call query%add_real2("sphxyz", "Sphere center coordinates, bohr (3, nsph)", self%sphxyz)
      if (.not. self%built) return

      call query%add_int_value("ngrid", "Number of surface grid points", self%ngrid)
      call query%add_real_scalar("area", "Total surface area, bohr**2", self%total_area)
      call query%add_real_scalar("volume", "Total enclosed volume, bohr**3", self%total_volume)

      call query%add_real2("xyz", "Grid point coordinates, bohr (3, ngrid)", self%xyz)
      call query%add_real("a", "Grid point area, bohr**2 (ngrid)", self%a)
      call query%add_real("f", "Gaussian switching factor per grid point (ngrid)", self%f)
      call query%add_real("xi0", "Gaussian width per grid point (ngrid)", self%xi0)
      call query%add_real2("normal0", "Outward unit normal per grid point (3, ngrid)", self%normal0)
      call query%add_real("v", "Grid point volume element, bohr**3 (ngrid)", self%v)
      call query%add_int("owner", "Sphere each grid point belongs to, 0-based (ngrid)", &
         & self%owner, zero_based=.true.)
      call query%add_real("asph", "Surface area per sphere, bohr**2 (nsph)", self%asph)

   end subroutine list_cavity_fields_base

   !> Print basic cavity information (grid points, total area, total volume)
   subroutine print_cavity_info(self, unit)
      class(cavity_type), intent(in) :: self
      integer, intent(in), optional :: unit
      integer :: iunit
      type(prettyprinter) :: pp

      iunit = cavity_unit(self)
      if (present(unit)) iunit = unit

      if (.not. allocated(self%total_area) .or. .not. allocated(self%total_volume)) then
         if (associated(self%ctx)) then
            call self%ctx%warning("Cavity not fully initialized")
         else
            write (iunit, "(a)") "[Warning] Cavity not fully initialized"
         end if
         return
      end if

      pp = new_prettyprinter(unit=iunit, fmt_len=20)

      call pp%blank()
      call pp%push("Results:")
      call pp%kv("Cavity points", self%ngrid)
      call pp%kv("Total area", self%total_area, "bohr^2")
      call pp%kv("Total volume", self%total_volume, "bohr^3")
      call pp%pop()
      call pp%blank()

   end subroutine print_cavity_info

   !> Default settings hook: the cavity has no parameter object
   !>
   !> @param[in] self Cavity instance
   function cavity_parameters_default(self) result(param)
      !> Cavity instance
      class(cavity_type), intent(in), target :: self
      !> Settings of the cavity, never associated here
      class(moist_model_parameters_type), pointer :: param

      param => null()

   end function cavity_parameters_default

   !> Print the cavity section: its kind, then the registered settings of
   !> `parameters()`
   !>
   !> @param[in] self Cavity instance
   !> @param[in] unit Output unit; defaults to the run context's unit
   !> @param[in] level Verbosity level; settings above it are left out, nothing
   !>                  prints at 0 and 1 prints the normal settings (default: the context's
   !>                  verbosity, every setting without a context)
   subroutine print_cavity_parameters(self, unit, level)
      !> Cavity instance
      class(cavity_type), intent(in), target :: self
      !> Output unit
      integer, intent(in), optional :: unit
      !> Verbosity level of the printout
      integer, intent(in), optional :: level

      !> Section printer
      type(prettyprinter) :: pp
      !> Registered settings, when the cavity has any
      class(moist_model_parameters_type), pointer :: param
      !> Effective output unit and verbosity level
      integer :: iu, lvl

      iu = cavity_unit(self)
      if (present(unit)) iu = unit
      lvl = huge(1)
      if (associated(self%ctx)) lvl = self%ctx%verbosity()
      if (present(level)) lvl = level
      ! Silent at 0; from 1 an explicit request shows at least the normal settings
      if (lvl < 1) return
      lvl = max(2, lvl)

      pp = new_prettyprinter(unit=iu)
      if (allocated(self%label)) then
         call pp%push("Cavity ("//self%label//"):")
      else
         call pp%push("Cavity:")
      end if
      param => self%parameters()
      if (associated(param)) call param%print_table(pp, lvl)
      call pp%pop()
      call pp%blank()

   end subroutine print_cavity_parameters

   !> Print the radius model at verbosity 2 and above
   !>
   !> Called by each cavity's `update` right after its radius model was
   !> updated: the printout lists the radii of the current structure
   !>
   !> @param[in] self Cavity instance
   subroutine print_cavity_radii(self)
      !> Cavity instance
      class(cavity_type), intent(in) :: self

      if (.not. associated(self%ctx) .or. .not. allocated(self%radius_model)) return
      if (self%ctx%verbosity(2)) call self%radius_model%print(unit=self%ctx%unit)

   end subroutine print_cavity_radii

   !> Write grid points to an XYZ file as helium atoms (debug visualization)
   subroutine write_cavity_xyz_debug(self, filename, error)
      use mctc_io_convert, only: autoaa
      class(cavity_type), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, i

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, "write_xyz_debug: cavity grid not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "write_xyz_debug: no grid points to write")
         return
      end if

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open XYZ file for writing: "//trim(filename))
         return
      end if

      write (unit, "(i0)") self%ngrid
      write (unit, "(a)") "drop cavity grid points as He (Angstrom)"
      do i = 1, self%ngrid
         write (unit, "(a2,1x,3f16.8)") "He", &
            self%xyz(1, i)*autoaa, &
            self%xyz(2, i)*autoaa, &
            self%xyz(3, i)*autoaa
      end do
      close (unit)

      write (cavity_unit(self), "(a,1x,a)") "[Info] Wrote cavity grid to", trim(filename)

   end subroutine write_cavity_xyz_debug

   !> Write grid points to a CSV file (debug visualization)
   subroutine write_cavity_csv_debug(self, filename, error)
      class(cavity_type), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: stat, i, unit

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, "write_csv_debug: cavity grid not allocated")
         return
      end if
      if (.not. allocated(self%a)) then
         call fatal_error(error, "write_csv_debug: point areas not allocated")
         return
      end if
      if (.not. allocated(self%owner)) then
         call fatal_error(error, "write_csv_debug: point owners not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "write_csv_debug: no grid points to write")
         return
      end if

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open CSV file for writing: "//trim(filename))
         return
      end if

      write (unit, "(a)") "ngrid,x,y,z,owner,area"

      do i = 1, self%ngrid
         write (unit, "(i0,7(',',g0))") i, &
            self%xyz(1, i), self%xyz(2, i), self%xyz(3, i), &
            self%owner(i), self%a(i)
      end do
      close (unit)

      write (cavity_unit(self), "(a,1x,a)") "[Info] Wrote cavity grid to", trim(filename)

   end subroutine write_cavity_csv_debug

   !> Write grid points to a PQR file (debug visualization)
   !>
   !> Grid points are written as HETATM records with:
   !> - positions converted from bohr to Angstrom
   !> - charge set to 0.0
   !> - radius set to the final adapted integration weight `a(i)` (area element
   !>   with switching function applied), also converted to Angstrom
   !>
   !> @param[in]  self      Cavity instance
   !> @param[in]  filename  Output PQR file path
   subroutine write_cavity_pqr_debug(self, filename, error)
      use mctc_io_convert, only: autoaa
      class(cavity_type), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, i

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, "write_pqr_debug: cavity grid not allocated")
         return
      end if
      if (.not. allocated(self%a)) then
         call fatal_error(error, "write_pqr_debug: point areas not allocated")
         return
      end if
      if (.not. allocated(self%owner)) then
         call fatal_error(error, "write_pqr_debug: point owners not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "write_pqr_debug: no grid points to write")
         return
      end if

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open PQR file for writing: "//trim(filename))
         return
      end if

      do i = 1, self%ngrid
         write (unit, "(a6,i5,1x,a4,a1,a3,1x,a1,i4,4x,3f8.3,f8.4,f7.4)") &
            "HETATM", i, "GP  ", " ", "GRD", "A", self%owner(i), &
            self%xyz(1, i)*autoaa, &
            self%xyz(2, i)*autoaa, &
            self%xyz(3, i)*autoaa, &
            0.0_wp, &
            (sqrt(self%a(i)/(2.0_wp*pi))*autoaa + 0.0001_wp)
      end do
      write (unit, "(a)") "END"
      close (unit)

      write (cavity_unit(self), "(a,1x,a)") "[Info] Wrote cavity PQR to", trim(filename)

   end subroutine write_cavity_pqr_debug

   !* ================================================================================= *!
   !*                                      Context                                      *!
   !* ================================================================================= *!

   !> Output unit for a cavity: the borrowed run context's unit when one is
   !> attached, otherwise the standard output unit
   !>
   !> - the base procedures below are reachable on a cavity that never went
   !>   through a constructor, so the association is checked, not assumed
   pure function cavity_unit(self) result(iunit)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Unit to write to
      integer :: iunit

      iunit = output_unit
      if (associated(self%ctx)) iunit = self%ctx%unit

   end function cavity_unit

   !> Fail unless a run context is attached
   !>
   !> @param[in]  self  Cavity instance
   !> @param[out] error Missing-context error
   subroutine require_cavity_context(self, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. associated(self%ctx)) call fatal_error(error, &
         & "Cavity has no context: pass one at construction or use it through a model")

   end subroutine require_cavity_context


end module moist_cavity_type
