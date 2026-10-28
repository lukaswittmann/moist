!> Checkpoints of the MOZ solves: field states and the abstract store
!>
!> - A state holds the physical t(r) of one solve, (npoint, ncol), with the
!>   discretization it lives on, its sites and provenance; c(r) and h(r) follow
!>   from t(r) by the closure and are not stored
!> - Physical t(r) stays valid across Ng-split changes; the short-range
!>   iterate is t_s = t - beta u_lr with the current long-range potential
!> - Only the discretization gates loading: kind, grid shape and spacing, site lists
!> - The solve decides convergence from the residual of the loaded field under
!>   the current equations; closure, theory, temperature and residual are provenance
!> - Models and solvents hold an optional store; constructors stay file-free
!>
!>    call solvent%set_checkpoint(store)
!>    call solvent%solve(error, parameters)
module moist_model_moz_checkpoint
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   implicit none(type, external)
   private

   public :: moz_checkpoint_state_type, moz_checkpoint_type
   public :: checkpoint_kind_vv, checkpoint_kind_uv1d, checkpoint_kind_uv3d

   !> Equation families of the states
   character(len=*), parameter :: checkpoint_kind_vv = "1D VV", checkpoint_kind_uv1d = "1D UV", &
      & checkpoint_kind_uv3d = "3D UV"

   !> Physical correlations of one solve with their discretization and provenance
   type :: moz_checkpoint_state_type
      !> Equation family: "1D VV", "1D UV" or "3D UV"
      character(len=:), allocatable :: kind
      !> Grid dimensions: [npts] radial, [nx, ny, nz] Cartesian, [ngrid] other point grids
      integer, allocatable :: grid_shape(:)
      !> Grid spacing, bohr; 0 for point grids without one
      real(wp) :: spacing = 0.0_wp
      !> Box corner of a Cartesian grid, bohr (3); provenance
      real(wp), allocatable :: origin(:)
      !> Atomic numbers of the solvent sites
      integer, allocatable :: solvent_sites(:)
      !> Atomic numbers of the solute sites; unallocated for VV
      integer, allocatable :: solute_sites(:)
      !> Solute coordinates, bohr (3, nsolute); provenance
      real(wp), allocatable :: solute_xyz(:, :)
      !> Closure and bulk theory; provenance
      character(len=:), allocatable :: closure, theory
      !> Temperature, K; provenance
      real(wp) :: temperature = 0.0_wp
      !> Last measured RMS residual and coupling of the stored field
      real(wp) :: residual = huge(1.0_wp), lambda = 1.0_wp
      !> Outer iterations of the solve that wrote the state
      integer :: iterations = 0
      !> Whether the solve had converged
      logical :: converged = .false.
      !> Shape of the fields, (npoint, ncol): expected by a header, stored by a loaded state
      integer, allocatable :: field_shape(:)
      !> Physical indirect correlation, (npoint, ncol)
      real(wp), allocatable :: tr(:, :)
   contains
      !> Whether a loaded state lives on the discretization of a header
      procedure :: check_compatible => state_check_compatible
      !> One-line provenance for messages
      procedure :: summary => state_summary
   end type moz_checkpoint_state_type

   !> Store of the states of one solve family
   !>
   !> - The owner sets `header` before each solve; drivers save copies with fields
   !> - Save and load failures are reported, never fatal to a solve
   type, abstract :: moz_checkpoint_type
      !> Iterations between intermediate saves; 0 saves only the final field
      integer :: interval = 0
      !> Discretization, sites and provenance of the states of the current solve
      type(moz_checkpoint_state_type) :: header
   contains
      !> Save one state
      procedure(save_interface), deferred :: save
      !> Load a state fitting the header
      procedure(load_interface), deferred :: load
      !> Location of saved states, for messages; empty when nothing is saved
      procedure(describe_interface), deferred :: describe
      !> Load an initial guess for the header, as the short-range iterate of the current solve
      procedure :: load_guess => checkpoint_load_guess
      !> Save the current short-range iterate of a solve with the header
      procedure :: record => checkpoint_record
      !> Saved location and interval for the solver notes; empty when nothing is saved
      procedure :: note => checkpoint_note
   end type moz_checkpoint_type

   !> Abstract interfaces of the deferred store bindings
   abstract interface
      !> Save one state
      !>
      !> @param[inout] self   store
      !> @param[in]    state  state with fields
      !> @param[out]   error  failed save
      subroutine save_interface(self, state, error)
         import :: moz_checkpoint_type, moz_checkpoint_state_type, error_type
         implicit none(type, external)
         !> Store
         class(moz_checkpoint_type), intent(inout) :: self
         !> State with fields
         type(moz_checkpoint_state_type), intent(in) :: state
         !> Failed save
         type(error_type), allocatable, intent(out) :: error
      end subroutine save_interface

      !> Load the first available state fitting the header
      !>
      !> @param[inout] self     store with a header
      !> @param[out]   state    loaded state; defined when found
      !> @param[out]   found    whether a fitting state was loaded
      !> @param[out]   source   location of the loaded state, for messages
      !> @param[out]   warning  unreadable or unfitting states that were skipped
      subroutine load_interface(self, state, found, source, warning)
         import :: moz_checkpoint_type, moz_checkpoint_state_type, error_type
         implicit none(type, external)
         !> Store
         class(moz_checkpoint_type), intent(inout) :: self
         !> Loaded state
         type(moz_checkpoint_state_type), intent(out) :: state
         !> Whether a state was loaded
         logical, intent(out) :: found
         !> Location of the loaded state
         character(len=:), allocatable, intent(out) :: source
         !> Skipped states
         type(error_type), allocatable, intent(out) :: warning
      end subroutine load_interface

      !> Location of saved states
      !>
      !> @param[in] self  store
      function describe_interface(self) result(location)
         import :: moz_checkpoint_type
         implicit none(type, external)
         !> Store
         class(moz_checkpoint_type), intent(in) :: self
         !> Location; empty when nothing is saved
         character(len=:), allocatable :: location
      end function describe_interface
   end interface

contains

   !> Check that a loaded state lives on the discretization of a header
   !>
   !> - Kind, grid shape, spacing and both site lists must agree
   !> - The state must store fields, of the header's field shape when it has one
   !> - Needs no fields in memory, so stores check before reading them
   !>
   !> @param[in]  self    loaded state
   !> @param[in]  header  discretization of the current solve
   !> @param[out] error   first mismatch
   subroutine state_check_compatible(self, header, error)
      !> Loaded state
      class(moz_checkpoint_state_type), intent(in) :: self
      !> Discretization of the current solve
      type(moz_checkpoint_state_type), intent(in) :: header
      !> First mismatch
      type(error_type), allocatable, intent(out) :: error

      if (.not. same_string(self%kind, header%kind)) then
         call fatal_error(error, "equation kind differs")
      else if (.not. same_integers(self%grid_shape, header%grid_shape)) then
         call fatal_error(error, "grid shape differs")
      else if (abs(self%spacing - header%spacing) > 1.0e-12_wp*max(1.0_wp, abs(header%spacing))) then
         call fatal_error(error, "grid spacing differs")
      else if (.not. same_integers(self%solvent_sites, header%solvent_sites)) then
         call fatal_error(error, "solvent sites differ")
      else if (.not. same_integers(self%solute_sites, header%solute_sites)) then
         call fatal_error(error, "solute sites differ")
      else if (.not. allocated(self%field_shape)) then
         call fatal_error(error, "no indirect correlation")
      else if (allocated(header%field_shape)) then
         if (.not. same_integers(self%field_shape, header%field_shape)) call fatal_error(error, "field shape differs")
      end if
   end subroutine state_check_compatible

   !> One-line provenance: closure, theory, temperature, residual and convergence
   !>
   !> @param[in] self  state
   function state_summary(self) result(text)
      !> State
      class(moz_checkpoint_state_type), intent(in) :: self
      !> Provenance
      character(len=:), allocatable :: text
      character(len=32) :: buf

      text = ""
      if (allocated(self%closure)) text = trim(self%closure)//", "
      if (allocated(self%theory)) text = text//trim(self%theory)//", "
      write (buf, "(f0.2,' K')") self%temperature
      text = text//trim(buf)
      write (buf, "(es9.2)") self%residual
      text = text//", residual "//trim(adjustl(buf))
      if (self%converged) then
         text = text//", converged"
      else
         text = text//", not converged"
      end if
   end function state_summary

   !> Load an initial guess for the header, as the short-range iterate of the current solve
   !>
   !> - Skipped states and loaded fields of another shape warn; the solve then starts cold
   !> - Source and provenance of a loaded guess are written at `level`
   !> - Messages carry the header kind as label
   !>
   !> @param[inout] self   store with its header set
   !> @param[in]    ctx    run context for warnings and messages
   !> @param[in]    beta   inverse thermal energy, 1/Hartree
   !> @param[in]    ur_lr  current long-range real-space potential, Hartree, (npoint, ncol)
   !> @param[in]    level  verbosity of the source message
   !> @param[out]   guess  short-range t_s(r), (npoint, ncol); unallocated without a guess
   subroutine checkpoint_load_guess(self, ctx, beta, ur_lr, level, guess)
      !> Store
      class(moz_checkpoint_type), intent(inout) :: self
      !> Run context
      type(moist_context_type), intent(in) :: ctx
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Current long-range real-space potential
      real(wp), intent(in) :: ur_lr(:, :)
      !> Verbosity of the source message
      integer, intent(in) :: level
      !> Short-range guess
      real(wp), allocatable, intent(out) :: guess(:, :)
      type(moz_checkpoint_state_type) :: state
      type(error_type), allocatable :: warning
      character(len=:), allocatable :: source, label
      logical :: found

      label = "MOZ"
      if (allocated(self%header%kind)) label = self%header%kind
      self%header%field_shape = shape(ur_lr)
      call self%load(state, found, source, warning)
      if (allocated(warning)) call ctx%warning(label//" checkpoint: "//warning%message)
      if (.not. found) return
      if (any(shape(state%tr) /= shape(ur_lr))) then
         call ctx%warning(label//" checkpoint: "//source//" has fields of another shape; ignored")
         return
      end if
      guess = state%tr - beta*ur_lr
      call ctx%message(label//" initial guess: "//source//" ("//state%summary()//")", level)
   end subroutine checkpoint_load_guess

   !> Save the current short-range iterate of a solve with the header
   !>
   !> - Stores the physical t = t_s + beta u_lr
   !>
   !> @param[inout] self       store with its header set
   !> @param[in]    ts         short-range indirect correlation, (npoint, ncol)
   !> @param[in]    beta       inverse thermal energy, 1/Hartree
   !> @param[in]    ur_lr      long-range real-space potential, Hartree, (npoint, ncol)
   !> @param[in]    rmse       last measured RMS residual
   !> @param[in]    lambda     coupling of the field
   !> @param[in]    niter      outer iterations so far
   !> @param[in]    converged  whether the field is converged
   !> @param[out]   error      failed save
   subroutine checkpoint_record(self, ts, beta, ur_lr, rmse, lambda, niter, converged, error)
      !> Store
      class(moz_checkpoint_type), intent(inout) :: self
      !> Short-range indirect correlation
      real(wp), intent(in) :: ts(:, :)
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Long-range real-space potential
      real(wp), intent(in) :: ur_lr(:, :)
      !> Last measured RMS residual
      real(wp), intent(in) :: rmse
      !> Coupling of the field
      real(wp), intent(in) :: lambda
      !> Outer iterations so far
      integer, intent(in) :: niter
      !> Whether the field is converged
      logical, intent(in) :: converged
      !> Failed save
      type(error_type), allocatable, intent(out) :: error
      type(moz_checkpoint_state_type) :: state

      state = self%header
      state%field_shape = shape(ts)
      state%tr = ts + beta*ur_lr
      state%residual = rmse
      state%lambda = lambda
      state%iterations = niter
      state%converged = converged
      call self%save(state, error)
   end subroutine checkpoint_record

   !> Saved location and interval for the solver notes
   !>
   !> @param[in] self  store
   function checkpoint_note(self) result(text)
      !> Store
      class(moz_checkpoint_type), intent(in) :: self
      !> Note; empty when nothing is saved
      character(len=:), allocatable :: text
      character(len=16) :: buf

      text = self%describe()
      if (len(text) == 0) return
      text = "Checkpoint: "//text
      if (self%interval > 0) then
         write (buf, "(i0)") self%interval
         text = text//", every "//trim(buf)//" iterations"
      end if
   end function checkpoint_note

   !> Whether two optional strings are both absent or equal
   !>
   !> @param[in] a  first string
   !> @param[in] b  second string
   pure function same_string(a, b) result(same)
      !> First string
      character(len=:), allocatable, intent(in) :: a
      !> Second string
      character(len=:), allocatable, intent(in) :: b
      !> Equality
      logical :: same
      same = allocated(a) .eqv. allocated(b)
      if (same .and. allocated(a)) same = a == b
   end function same_string

   !> Whether two optional integer lists are both absent or equal
   !>
   !> @param[in] a  first list
   !> @param[in] b  second list
   pure function same_integers(a, b) result(same)
      !> First list
      integer, allocatable, intent(in) :: a(:)
      !> Second list
      integer, allocatable, intent(in) :: b(:)
      !> Equality
      logical :: same
      same = allocated(a) .eqv. allocated(b)
      if (.not. same .or. .not. allocated(a)) return
      same = size(a) == size(b)
      if (same) same = all(a == b)
   end function same_integers

end module moist_model_moz_checkpoint
