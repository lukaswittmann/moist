!> Checkpoints of the MOZ solves: field states and the abstract store
!>
!> - A state holds the physical t(r) of one solve, (npoint, ncol), with the
!>   discretization it lives on (its grid layout), its sites and provenance; c(r) and
!>   h(r) follow from t(r) by the closure and are not stored
!> - Physical t(r) stays valid across Ng-split changes; the short-range
!>   iterate is t_s = t - beta u_lr with the current long-range potential
!> - A state seeds a solve of the same kind with the same site lists, on any grid:
!>   - `identity`: atom-centered grids with the same atomic grids on built-in Lebedev
!>     angular grids, both with the position of each point among the unpruned points,
!>     which moves with its atom; the points retained at both structures are seeded, and
!>     no other method is tried for such a state
!>   - `copy`: current points that coincide with a stored point after the rigid fit of the
!>     solute take its value, for every kind of grid
!>   - from Cartesian grids `fourier` onto Cartesian grids and `trilinear` onto others,
!>     `points` (moving least squares) from molecular and other point grids: the stored
!>     solute is fitted rigidly onto the current one, the short-range part t - tau of the
!>     stored field is evaluated at the current points in the stored frame, and the stored
!>     tail tau is added back
!>   - `spline`: radial grids, spline of t through the stored nodes
!>   - points the stored field does not cover start cold: c_s = 0 in 1D, the selected cold
!>     guess in 3D
!> - The solve decides convergence from the residual of the loaded field under
!>   the current equations; closure, theory and residual are provenance
!> - Models and solvents hold an optional store; constructors stay file-free
!>
!>    call solvent%set_checkpoint(store)
!>    call solvent%solve(error, parameters)
module moist_model_moz_checkpoint
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use mctc_io_codata2018, only: Boltzmann_constant, atomic_unit_of_energy
   use moist_math_grid_layout, only: moist_math_grid_layout_type, layout_cartesian, layout_radial, layout_points
   use moist_math_grid_transfer, only: rigid_fit_type, fit_rigid, grid_transfer_type, grid_transfer_options_type, &
      & grid_transfer_report_type
   use moist_model_moz_potential_kernel_ng, only: ng_real_tail
   implicit none(type, external)
   private

   public :: moz_checkpoint_state_type, moz_checkpoint_type
   public :: checkpoint_kind_vv, checkpoint_kind_uv1d, checkpoint_kind_uv3d
   public :: tail_known, tail_zero, tail_unknown
   public :: moz_checkpoint_guess

   !> Equation families of the states
   character(len=*), parameter :: checkpoint_kind_vv = "1D VV", checkpoint_kind_uv1d = "1D UV", &
      & checkpoint_kind_uv3d = "3D UV"

   !> Long-range tail of a 3D UV state: reproduced from its data, absent, or not reproducible
   character(len=*), parameter :: tail_known = "known", tail_zero = "zero", tail_unknown = "unknown"

   !> Physical correlations of one solve with their discretization and provenance
   type :: moz_checkpoint_state_type
      !> Equation family: "1D VV", "1D UV" or "3D UV"
      character(len=:), allocatable :: kind
      !> Discretization of the fields: dimensions, spacing, and the points of point grids
      type(moist_math_grid_layout_type) :: grid
      !> Atomic numbers of the solvent sites
      integer, allocatable :: solvent_sites(:)
      !> Atomic numbers of the solute sites; unallocated for VV
      integer, allocatable :: solute_sites(:)
      !> Solute coordinates, bohr (3, nsolute)
      real(wp), allocatable :: solute_xyz(:, :)
      !> 3D UV long-range tail: "known", "zero" or "unknown"; unallocated for 1D states
      character(len=:), allocatable :: tail
      !> Ng screening parameter of a known tail, 1/bohr
      real(wp) :: alpha = 0.0_wp
      !> Tail charges of the solute atoms (nsolute) and of the solvent sites (ncol) of a known tail, e
      real(wp), allocatable :: solute_charges(:), solvent_charges(:)
      !> Closure and bulk theory; provenance
      character(len=:), allocatable :: closure, theory
      !> Temperature, K; provenance, and needed to rebuild a known tail
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
      !> Whether a state with fields can be saved and loaded again
      procedure :: check_complete => state_check_complete
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
      !> Record the long-range tail of the current 3D solve in the header
      procedure :: set_tail => checkpoint_set_tail
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

   !> Check that a loaded state can seed the solve of a header
   !>
   !> - Kind and both site lists must agree, and the column count of the fields
   !> - The stored layout must be consistent: one field row per stored point, point
   !>   coordinates, weights, owners, nodes and Fourier basis of matching sizes,
   !>   increasing point positions; the layout of the header must be valid too
   !> - With the field shape of the header set, a transfer method must apply and the state
   !>   must cover at least one point of the header grid; on 3D grids every method except
   !>   `identity` needs the solute coordinates of both sides and their rigid fit
   !> - Needs no fields in memory, so stores check before reading them, and skip states
   !>   that cannot seed the solve
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
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      type(rigid_fit_type) :: fit
      character(len=:), allocatable :: reason, method

      call check_metadata(self, header, error)
      if (allocated(error) .or. .not. allocated(header%field_shape)) return
      ! Whether any point is covered; point grids fit stencils only up to the first accepted one
      call plan_transfer(self, header, header%field_shape(1), method, reason, transfer, report, fit, &
         & check_only=.true.)
      if (len(method) == 0) call fatal_error(error, reason)
   end subroutine state_check_compatible

   !> Whether a state with fields can be saved and loaded again
   !>
   !> - Kind, grid shape, solvent sites and t(r); a stored field shape equals that of t(r)
   !> - The layout passes the check a load applies to the stored grid, with one point per
   !>   row of t(r)
   !>
   !> @param[in]  self   state to save
   !> @param[out] error  first missing or inconsistent part
   subroutine state_check_complete(self, error)
      !> State to save
      class(moz_checkpoint_state_type), intent(in) :: self
      !> First missing or inconsistent part
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: reason

      if (.not. allocated(self%kind) .or. .not. allocated(self%grid%shape) .or. &
         & .not. allocated(self%solvent_sites) .or. .not. allocated(self%tr)) then
         call fatal_error(error, "Checkpoint state needs a kind, grid shape, solvent sites and t(r)")
         return
      end if
      if (allocated(self%field_shape)) then
         if (size(self%field_shape) /= 2) then
            call fatal_error(error, "Checkpoint state: field shape differs from t(r)")
            return
         else if (any(self%field_shape /= shape(self%tr))) then
            call fatal_error(error, "Checkpoint state: field shape differs from t(r)")
            return
         end if
      end if
      if (.not. valid_layout(self%grid, size(self%tr, 1), reason)) then
         call fatal_error(error, "Checkpoint state: "//reason)
      end if
   end subroutine state_check_complete

   !> The checks of `state_check_compatible` that need no transfer: kind, sites, field
   !> shapes and the consistency of both layouts
   !>
   !> @param[in]  self    loaded state
   !> @param[in]  header  discretization of the current solve
   !> @param[out] error   first mismatch
   subroutine check_metadata(self, header, error)
      !> Loaded state
      class(moz_checkpoint_state_type), intent(in) :: self
      !> Discretization of the current solve
      type(moz_checkpoint_state_type), intent(in) :: header
      !> First mismatch
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: fshape(:)
      character(len=:), allocatable :: reason

      call state_field_shape(self, fshape)
      if (.not. same_string(self%kind, header%kind)) then
         call fatal_error(error, "equation kind differs")
      else if (.not. same_integers(self%solvent_sites, header%solvent_sites)) then
         call fatal_error(error, "solvent sites differ")
      else if (.not. same_integers(self%solute_sites, header%solute_sites)) then
         call fatal_error(error, "solute sites differ")
      else if (.not. allocated(fshape)) then
         call fatal_error(error, "no indirect correlation")
      else if (size(fshape) /= 2) then
         call fatal_error(error, "field shape differs")
      else if (allocated(header%field_shape)) then
         if (size(header%field_shape) /= 2) then
            call fatal_error(error, "field shape differs")
         else if (fshape(2) /= header%field_shape(2)) then
            call fatal_error(error, "field shape differs")
         end if
      end if
      if (allocated(error)) return
      if (.not. valid_layout(self%grid, fshape(1), reason)) then
         call fatal_error(error, "stored grid: "//reason)
         return
      end if
      if (.not. allocated(header%field_shape)) return
      if (.not. valid_layout(header%grid, header%field_shape(1), reason)) then
         call fatal_error(error, "current grid: "//reason)
         return
      end if
   end subroutine check_metadata

   !> Whether a layout is consistent with `nrow` field rows
   !>
   !> - Cartesian, radial and point layouts pass the layout's own `validate` and count
   !>   `nrow` points; Cartesian and radial layouts give their dimensions
   !> - Positions among the unpruned points, when given, one per row and increasing
   !> - Layouts of no kind: dimensions, when given, multiply to the rows
   !>
   !> @param[in]  grid    layout
   !> @param[in]  nrow    field rows
   !> @param[out] reason  first inconsistency
   function valid_layout(grid, nrow, reason) result(valid)
      !> Layout
      type(moist_math_grid_layout_type), intent(in) :: grid
      !> Field rows
      integer, intent(in) :: nrow
      !> First inconsistency
      character(len=:), allocatable, intent(out) :: reason
      !> Consistency
      logical :: valid
      integer :: n

      valid = .false.
      reason = "fields do not match the grid"
      if (allocated(grid%point_index)) then
         n = size(grid%point_index)
         if (n /= nrow) return
         if (n < 1) then
            reason = "no points"
            return
         end if
         if (any(grid%point_index(2:) <= grid%point_index(:n - 1))) then
            reason = "point positions do not increase"
            return
         end if
      end if
      select case (layout_kind(grid))
      case (layout_cartesian, layout_radial, layout_points)
         if (layout_kind(grid) /= layout_points .and. .not. allocated(grid%shape)) return
         call grid%validate(n, reason)
         if (len(reason) > 0) return
         reason = "fields do not match the grid"
         if (n /= nrow) return
      case default
         if (allocated(grid%shape) .and. .not. allocated(grid%point_index)) then
            if (dimension_product(grid%shape) /= int(nrow, int64)) return
         end if
      end select
      valid = .true.
      reason = ""
   end function valid_layout

   !> Product of grid dimensions in 64 bits; -1 for a dimension below one or an overflow
   !>
   !> @param[in] dims  dimensions
   pure function dimension_product(dims) result(n)
      !> Dimensions
      integer, intent(in) :: dims(:)
      !> Product
      integer(int64) :: n
      integer :: i

      n = 1_int64
      do i = 1, size(dims)
         if (dims(i) < 1) then
            n = -1_int64
            return
         end if
         if (n > huge(n)/int(dims(i), int64)) then
            n = -1_int64
            return
         end if
         n = n*int(dims(i), int64)
      end do
   end function dimension_product

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
   !> - Stores offer the first state that can seed the solve; skipped states warn
   !> - `moz_checkpoint_guess` moves it onto the current grid; `matched` marks the points with a
   !>   value, and without `matched` only a state covering every point is used
   !> - Source, provenance, method and coverage of a loaded guess are written at `level`
   !> - Messages carry the header kind as label
   !>
   !> @param[inout] self     store with its header set; its field shape becomes that of
   !>                        `ur_lr`
   !> @param[in]    ctx      run context for warnings, messages and threads
   !> @param[in]    beta     inverse thermal energy, 1/Hartree
   !> @param[in]    ur_lr    current long-range real-space potential, Hartree, (npoint, ncol)
   !> @param[in]    level    verbosity of the source message
   !> @param[out]   guess    short-range t_s(r), (npoint, ncol); unallocated without a guess
   !> @param[out]   matched  points with a stored value, (npoint); unallocated without a guess
   subroutine checkpoint_load_guess(self, ctx, beta, ur_lr, level, guess, matched)
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
      !> Points with a stored value
      logical, allocatable, intent(out), optional :: matched(:)
      type(moz_checkpoint_state_type) :: state
      type(error_type), allocatable :: warning
      character(len=:), allocatable :: source, label, method, reason
      logical, allocatable :: covered(:)
      character(len=64) :: buf
      integer :: ncovered
      logical :: found

      label = "MOZ"
      if (allocated(self%header%kind)) label = self%header%kind
      self%header%field_shape = shape(ur_lr)
      call self%load(state, found, source, warning)
      if (allocated(warning)) call ctx%warning(label//" checkpoint: "//warning%message)
      if (.not. found) return
      call moz_checkpoint_guess(state, self%header, ctx, beta, ur_lr, guess, covered, method, reason)
      if (.not. allocated(guess)) then
         call ctx%warning(label//" checkpoint: "//source//" cannot seed this grid ("//reason//"); ignored")
         return
      end if
      ncovered = count(covered)
      if (ncovered < size(covered) .and. .not. present(matched)) then
         deallocate (guess)
         call ctx%warning(label//" checkpoint: "//source//" covers only part of this grid; ignored")
         return
      end if
      if (present(matched)) call move_alloc(covered, matched)
      if (ncovered < size(ur_lr, 1)) then
         write (buf, "(i0,' of ',i0,' points')") ncovered, size(ur_lr, 1)
         method = method//", "//trim(buf)//", the others cold"
      end if
      call ctx%message(label//" initial guess: "//source//" ("//state%summary()//", "//method//")", level)
   end subroutine checkpoint_load_guess

   !> Initial guess for the header from a stored state
   !>
   !> - The method is `identity` when it applies, else the grid transfer (`copy`,
   !>   `fourier`, `trilinear`, `points`, `spline`; see the module description); without
   !>   one, or without a covered point, `guess` stays unallocated and `reason` says why
   !> - `guess` is the short-range t_s = t - beta ur_lr on the header points that the
   !>   state covers and zero elsewhere; `matched` marks the covered points
   !> - Needs only the state and the header: host codes may call it on states they keep
   !>
   !> @param[in]  state    stored state with fields
   !> @param[in]  header   discretization of the current solve: kind, sites, grid layout and
   !>                      solute coordinates
   !> @param[in]  ctx      run context for the thread count
   !> @param[in]  beta     inverse thermal energy of the current solve, 1/Hartree
   !> @param[in]  ur_lr    current long-range real-space potential, Hartree, (npoint, ncol)
   !> @param[out] guess    short-range t_s on the header points, (npoint, ncol)
   !> @param[out] matched  covered points (npoint)
   !> @param[out] method   method used; empty without one
   !> @param[out] reason   why the state cannot seed the header; empty with a method
   subroutine moz_checkpoint_guess(state, header, ctx, beta, ur_lr, guess, matched, method, reason)
      !> Stored state
      type(moz_checkpoint_state_type), intent(in) :: state
      !> Current discretization
      type(moz_checkpoint_state_type), intent(in) :: header
      !> Run context
      type(moist_context_type), intent(in) :: ctx
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Current long-range real-space potential
      real(wp), intent(in) :: ur_lr(:, :)
      !> Short-range guess
      real(wp), allocatable, intent(out) :: guess(:, :)
      !> Covered points
      logical, allocatable, intent(out) :: matched(:)
      !> Method used
      character(len=:), allocatable, intent(out) :: method
      !> Why no method applies
      character(len=:), allocatable, intent(out) :: reason
      type(moz_checkpoint_state_type) :: target
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      type(rigid_fit_type) :: fit
      type(error_type), allocatable :: error
      real(wp), allocatable :: tr(:, :), xs(:, :), ys(:, :), tau(:, :)
      logical, allocatable :: covered(:)
      integer, allocatable :: idx(:), fshape(:)
      integer :: npoint, i

      method = ""
      npoint = size(ur_lr, 1)
      target = header
      target%field_shape = shape(ur_lr)
      ! The transfer below decides the method; one prepare serves the check and the move
      call check_metadata(state, target, error)
      if (allocated(error)) then
         reason = error%message
         return
      end if
      if (.not. allocated(state%tr)) then
         reason = "no indirect correlation"
         return
      end if
      call state_field_shape(state, fshape)
      if (any(shape(state%tr) /= fshape)) then
         reason = "fields do not match the stored grid"
         return
      end if
      call plan_transfer(state, target, npoint, method, reason, transfer, report, fit, ctx%get_num_threads())
      if (len(method) == 0) return
      if (method == "identity") then
         call match_points(state%grid%point_index, target%grid%point_index, state%tr, tr, covered)
      else
         allocate (tr(npoint, size(state%tr, 2)), source=0.0_wp)
         if (report%ninterpolated == 0 .or. layout_kind(state%grid) == layout_radial) then
            ! Copies and radial splines move the physical t
            call transfer%apply(state%tr, tr, error, covered=covered)
         else
            ! The short-range part t - tau moves; the stored tail is added back at the
            ! current points in the stored frame
            call state%grid%points(xs)
            call stored_tail(state, xs, tau)
            call transfer%apply(state%tr - tau, tr, error, covered=covered)
            if (.not. allocated(error)) then
               call target%grid%points(xs)
               call fit%map(xs, ys)
               idx = pack([(i, i=1, npoint)], covered)
               call stored_tail(state, ys(:, idx), tau)
               tr(idx, :) = tr(idx, :) + tau
            end if
         end if
      end if
      if (allocated(error)) then
         method = ""
         reason = error%message
         return
      end if
      allocate (guess(npoint, size(ur_lr, 2)), source=0.0_wp)
      where (spread(covered, 2, size(ur_lr, 2))) guess = tr - beta*ur_lr
      call move_alloc(covered, matched)
      reason = ""
   end subroutine moz_checkpoint_guess

   !> Field shape of a state: the stored one, else that of its t (a state built in memory);
   !> left unallocated without either
   !>
   !> @param[in]  state   state
   !> @param[out] fshape  field shape
   pure subroutine state_field_shape(state, fshape)
      !> State
      class(moz_checkpoint_state_type), intent(in) :: state
      !> Field shape
      integer, allocatable, intent(out) :: fshape(:)
      if (allocated(state%field_shape)) then
         fshape = state%field_shape
      else if (allocated(state%tr)) then
         fshape = shape(state%tr)
      end if
   end subroutine state_field_shape

   !> Method that moves a state onto the points of a header, from metadata alone
   !>
   !> - `identity` for atom-centered grids with the same atomic grids on built-in Lebedev
   !>   angular grids, both with point positions among the unpruned points; it is final,
   !>   also when it covers no point; otherwise the grid transfer from the stored layout
   !>   onto the header layout, after the rigid fit of the solute for 3D grids
   !> - Interpolating a 3D field needs its long-range tail; without one only exact copies
   !>   apply
   !> - A method that covers no point is no method
   !>
   !> @param[in]  state     loaded state; its fields are not needed
   !> @param[in]  header    discretization of the current solve
   !> @param[in]  npoint    points of the current solve
   !> @param[out] method    "identity" or the methods of the grid transfer; empty without one
   !> @param[out] reason    why no method applies
   !> @param[out] transfer  prepared grid transfer; not prepared for `identity`
   !> @param[out] report    coverage and methods of the grid transfer
   !> @param[out] fit       rigid map from the current frame into the stored frame (3D)
   !> @param[in]  nthreads  threads of the interpolation; 0 or absent for the default
   !> @param[in]  check_only  decide only whether a point is covered; a transfer from point
   !>                         layouts that leaves targets to the fit then cannot be applied
   subroutine plan_transfer(state, header, npoint, method, reason, transfer, report, fit, nthreads, check_only)
      !> Loaded state
      class(moz_checkpoint_state_type), intent(in) :: state
      !> Current discretization
      type(moz_checkpoint_state_type), intent(in) :: header
      !> Points of the current solve
      integer, intent(in) :: npoint
      !> Method
      character(len=:), allocatable, intent(out) :: method
      !> Why no method applies
      character(len=:), allocatable, intent(out) :: reason
      !> Prepared grid transfer
      type(grid_transfer_type), intent(out) :: transfer
      !> Coverage and methods of the grid transfer
      type(grid_transfer_report_type), intent(out) :: report
      !> Rigid map
      type(rigid_fit_type), intent(out) :: fit
      !> Threads of the interpolation
      integer, intent(in), optional :: nthreads
      !> Only whether a point is covered
      logical, intent(in), optional :: check_only
      type(grid_transfer_options_type) :: options
      type(error_type), allocatable :: error
      character(len=:), allocatable :: tail_reason, layout_reason
      logical, allocatable :: covered(:)
      integer :: n

      method = ""
      reason = "no method for these grids"
      if (present(nthreads)) options%nthreads = nthreads
      if (present(check_only)) options%check_only = check_only
      ! Atom-centered grids with the same atomic grids: points follow their atoms
      if (allocated(state%grid%point_index) .and. allocated(header%grid%point_index)) then
         if (same_shells(state, header) .and. known_angular_grid(state) .and. known_angular_grid(header)) then
            if (size(header%grid%point_index) == npoint) then
               allocate (covered(npoint))
               call match_cover(state%grid%point_index, header%grid%point_index, covered)
               if (any(covered)) then
                  method = "identity"
                  reason = ""
               else
                  reason = "covers no point of this grid"
               end if
               return
            end if
         end if
      end if
      select case (layout_kind(state%grid))
      case (layout_radial)
         if (layout_kind(header%grid) /= layout_radial .or. .not. allocated(state%grid%r) &
            & .or. .not. allocated(header%grid%r)) then
            reason = "radial nodes missing"
            return
         end if
         if (size(header%grid%r) /= npoint) return
         call transfer%prepare(state%grid, header%grid, error, options=options)
      case (layout_cartesian, layout_points)
         if (.not. allocated(state%solute_xyz) .or. .not. allocated(header%solute_xyz)) then
            reason = "solute coordinates missing"
            return
         end if
         if (any(shape(state%solute_xyz) /= shape(header%solute_xyz))) then
            reason = "solute coordinates differ in shape"
            return
         end if
         if (.not. all(ieee_is_finite(state%solute_xyz)) .or. .not. all(ieee_is_finite(header%solute_xyz))) then
            reason = "solute coordinates not finite"
            return
         end if
         call header%grid%validate(n, layout_reason, coordinates=.true.)
         if (n /= npoint) then
            reason = "the current grid has no points"
            return
         end if
         call fit_rigid(state%solute_xyz, header%solute_xyz, fit, error)
         if (allocated(error)) then
            reason = error%message
            return
         end if
         ! Without a tail only exact copies move the physical t
         options%copy_only = .not. tail_usable(state, tail_reason)
         call transfer%prepare(state%grid, header%grid, error, options=options, map=fit)
      case default
         reason = "the stored grid is not described"
         return
      end select
      if (allocated(error)) then
         reason = error%message
         return
      end if
      call transfer%get_report(report)
      if (report%ncopied + report%ninterpolated > 0) then
         method = report%method
         reason = ""
      else if (report%nexcluded > 0) then
         reason = tail_reason
      else
         reason = report%reason
      end if
   end subroutine plan_transfer

   !> Whether the long-range tail of a stored 3D state can be rebuilt
   !>
   !> - "zero" needs nothing; "known" needs alpha > 0, a temperature, one solute charge
   !>   per stored solute atom and one solvent charge per field column
   !>
   !> @param[in]    state   loaded state
   !> @param[inout] reason  why not; unchanged when usable as "zero"
   logical function tail_usable(state, reason)
      !> Loaded state
      class(moz_checkpoint_state_type), intent(in) :: state
      !> Why not
      character(len=:), allocatable, intent(inout) :: reason
      integer, allocatable :: fshape(:)

      tail_usable = .false.
      if (.not. allocated(state%tail)) then
         reason = "long-range tail of the stored state unknown"
         return
      end if
      if (state%tail == tail_zero) then
         tail_usable = .true.
         return
      end if
      if (state%tail /= tail_known) then
         reason = "long-range tail of the stored state unknown"
         return
      end if
      reason = "long-range data of the stored state incomplete"
      if (.not. (state%alpha > 0.0_wp) .or. .not. (state%temperature > 0.0_wp)) return
      if (.not. ieee_is_finite(state%alpha) .or. .not. ieee_is_finite(state%temperature)) return
      if (.not. allocated(state%solute_charges) .or. .not. allocated(state%solvent_charges)) return
      call state_field_shape(state, fshape)
      if (.not. allocated(state%solute_xyz) .or. .not. allocated(fshape)) return
      if (size(state%solute_charges) /= size(state%solute_xyz, 2)) return
      if (size(fshape) /= 2) return
      if (size(state%solvent_charges) /= fshape(2)) return
      if (.not. all(ieee_is_finite(state%solute_charges)) .or. .not. all(ieee_is_finite(state%solvent_charges))) return
      tail_usable = .true.
      reason = ""
   end function tail_usable

   !> Long-range tail beta_s q_v sum_a q_a erf(alpha r_a)/r_a of a stored state at points
   !> in its frame; zero for a "zero" tail; the tail must pass `tail_usable`
   !>
   !> @param[in]  state  loaded state with a known or zero tail
   !> @param[in]  x      points in the stored frame, bohr (3, n)
   !> @param[out] tau    tail, (n, ncol)
   subroutine stored_tail(state, x, tau)
      !> Loaded state
      type(moz_checkpoint_state_type), intent(in) :: state
      !> Points
      real(wp), intent(in) :: x(:, :)
      !> Tail
      real(wp), allocatable, intent(out) :: tau(:, :)
      real(wp), allocatable :: phi(:), f(:)
      real(wp) :: beta_s
      integer :: a, v

      ! Called after tail_usable: a zero tail, or complete data of a known one
      allocate (tau(size(x, 2), size(state%tr, 2)), source=0.0_wp)
      if (state%tail /= tail_known) return
      beta_s = 1.0_wp/(Boltzmann_constant/atomic_unit_of_energy*state%temperature)
      allocate (phi(size(x, 2)), source=0.0_wp)
      allocate (f(size(x, 2)))
      do a = 1, size(state%solute_charges)
         call ng_real_tail(norm2(x - spread(state%solute_xyz(:, a), 2, size(x, 2)), dim=1), state%alpha, f)
         phi = phi + state%solute_charges(a)*f
      end do
      do v = 1, size(tau, 2)
         tau(:, v) = beta_s*state%solvent_charges(v)*phi
      end do
   end subroutine stored_tail

   !> Current points whose position is among the stored point positions
   !>
   !> @param[in]  stored   increasing positions of the stored points
   !> @param[in]  target   increasing positions of the current points
   !> @param[out] covered  current points with a stored position
   pure subroutine match_cover(stored, target, covered)
      !> Positions of the stored points
      integer, intent(in) :: stored(:)
      !> Positions of the current points
      integer, intent(in) :: target(:)
      !> Covered current points
      logical, intent(out) :: covered(:)
      integer :: i, j

      covered = .false.
      i = 1
      do j = 1, size(target)
         do while (i <= size(stored))
            if (stored(i) >= target(j)) exit
            i = i + 1
         end do
         if (i > size(stored)) exit
         covered(j) = stored(i) == target(j)
      end do
   end subroutine match_cover

   !> Kind of a layout, empty when undescribed
   !>
   !> @param[in] grid  layout
   pure function layout_kind(grid) result(kind)
      !> Layout
      type(moist_math_grid_layout_type), intent(in) :: grid
      !> Kind
      character(len=:), allocatable :: kind
      kind = ""
      if (allocated(grid%kind)) kind = grid%kind
   end function layout_kind

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

   !> Record the long-range tail of the current 3D solve in the header
   !>
   !> - The tail `q_v sum_a q_a erf(alpha r_a)/r_a`, rebuilt on the header points with the
   !>   header solute coordinates, must reproduce `ur_lr`; then it is "known"
   !> - Without tail charges or with `alpha = 0`, a vanishing `ur_lr` is "zero"
   !> - Otherwise the tail is "unknown", also with tail charges on a header without points
   !>   or solute coordinates
   !>
   !> @param[inout] self             store with its header set
   !> @param[in]    ur_lr            long-range real-space potential of the solve, Hartree, (npoint, ncol)
   !> @param[in]    alpha            Ng screening parameter, 1/bohr; 0 without a tail
   !> @param[in]    solute_charges   tail charges of the solute atoms, e (nsolute)
   !> @param[in]    solvent_charges  charges of the solvent sites, e (ncol)
   subroutine checkpoint_set_tail(self, ur_lr, alpha, solute_charges, solvent_charges)
      !> Store
      class(moz_checkpoint_type), intent(inout) :: self
      !> Long-range real-space potential
      real(wp), intent(in) :: ur_lr(:, :)
      !> Ng screening parameter
      real(wp), intent(in) :: alpha
      !> Tail charges of the solute atoms
      real(wp), intent(in), optional :: solute_charges(:)
      !> Charges of the solvent sites
      real(wp), intent(in), optional :: solvent_charges(:)
      real(wp), allocatable :: xyz(:, :), phi(:), f(:)
      real(wp) :: scale
      integer :: a, v

      associate (header => self%header)
         header%alpha = 0.0_wp
         if (allocated(header%solute_charges)) deallocate (header%solute_charges)
         if (allocated(header%solvent_charges)) deallocate (header%solvent_charges)
         header%tail = tail_unknown
         scale = 1.0_wp
         if (size(ur_lr) > 0) scale = max(1.0_wp, maxval(abs(ur_lr)))
         if (alpha <= 0.0_wp .or. .not. present(solute_charges) .or. .not. present(solvent_charges)) then
            if (all(ur_lr == 0.0_wp)) header%tail = tail_zero
            return
         end if
         if (.not. allocated(header%solute_xyz)) return
         if (size(solute_charges) /= size(header%solute_xyz, 2) .or. size(solvent_charges) /= size(ur_lr, 2)) return
         call header%grid%points(xyz)
         if (size(xyz, 2) == 0 .or. size(xyz, 2) /= size(ur_lr, 1)) return
         allocate (phi(size(xyz, 2)), source=0.0_wp)
         allocate (f(size(xyz, 2)))
         do a = 1, size(solute_charges)
            call ng_real_tail(norm2(xyz - spread(header%solute_xyz(:, a), 2, size(xyz, 2)), dim=1), alpha, f)
            phi = phi + solute_charges(a)*f
         end do
         do v = 1, size(ur_lr, 2)
            if (any(abs(solvent_charges(v)*phi - ur_lr(:, v)) > 1.0e-10_wp*scale)) return
         end do
         header%tail = tail_known
         header%alpha = alpha
         header%solute_charges = solute_charges
         header%solvent_charges = solvent_charges
      end associate
   end subroutine checkpoint_set_tail

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

   !> Whether shell point counts name the built-in angular directions
   !>
   !> @param[in] state  molecular state or header
   pure logical function known_angular_grid(state) result(known)
      !> Molecular state or header
      class(moz_checkpoint_state_type), intent(in) :: state
      known = .false.
      if (allocated(state%grid%angular_grid)) known = state%grid%angular_grid == "lebedev"
   end function known_angular_grid

   !> Stored values at the points with the same position
   !>
   !> @param[in]  stored   increasing positions of the stored points
   !> @param[in]  target   increasing positions of the points of the current grid
   !> @param[in]  values   stored field, (nstored, ncol)
   !> @param[out] tr       field on the current points, (ntarget, ncol); zero where not covered
   !> @param[out] covered  current points with a stored value, (ntarget)
   pure subroutine match_points(stored, target, values, tr, covered)
      !> Positions of the stored points
      integer, intent(in) :: stored(:)
      !> Positions of the current points
      integer, intent(in) :: target(:)
      !> Stored field
      real(wp), intent(in) :: values(:, :)
      !> Field on the current points
      real(wp), allocatable, intent(out) :: tr(:, :)
      !> Covered current points
      logical, allocatable, intent(out) :: covered(:)
      integer :: i, j

      allocate (tr(size(target), size(values, 2)), source=0.0_wp)
      allocate (covered(size(target)), source=.false.)
      ! Both lists increase: one merge pass
      i = 1
      do j = 1, size(target)
         do while (i <= size(stored))
            if (stored(i) >= target(j)) exit
            i = i + 1
         end do
         if (i > size(stored)) exit
         if (stored(i) == target(j)) then
            tr(j, :) = values(i, :)
            covered(j) = .true.
         end if
      end do
   end subroutine match_points

   !> Whether two states use the same atomic grids, or neither has any
   !>
   !> @param[in] a  first state
   !> @param[in] b  second state
   pure function same_shells(a, b) result(same)
      !> First state
      class(moz_checkpoint_state_type), intent(in) :: a
      !> Second state
      type(moz_checkpoint_state_type), intent(in) :: b
      !> Equality
      logical :: same
      associate (ga => a%grid, gb => b%grid)
         same = same_integers(ga%shell_offset, gb%shell_offset) .and. same_integers(ga%shell_points, gb%shell_points)
         if (.not. same) return
         same = allocated(ga%shell_radius) .eqv. allocated(gb%shell_radius)
         if (.not. same .or. .not. allocated(ga%shell_radius)) return
         same = size(ga%shell_radius) == size(gb%shell_radius)
         if (same) same = all(abs(ga%shell_radius - gb%shell_radius) <= 1.0e-12_wp*max(1.0_wp, abs(gb%shell_radius)))
      end associate
   end function same_shells

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
