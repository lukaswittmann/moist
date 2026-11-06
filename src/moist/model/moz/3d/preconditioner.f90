!> Cartesian auxiliary-space preconditioner of the 3D UV Jacobian
!>
!> - Approximates the inverse of J v = v - T'(c' v) with periodic Cartesian lattices; the
!>   molecular residual, Jacobian and convergence test stay as they are
!> - Lattice: a cube of side L = max_d 2 pi/dk_d of the grid's Fourier modes (n_d dr on a
!>   Cartesian grid), centred on the points, with n cells per axis; n is capped at the
!>   finest lattice the grid's modes resolve, L/dr rounded down to even
!> - Transfers (`moist_math_grid_deposit`): P trilinear from the cell centres, R the
!>   weighted deposition w_i P_ia / m_a; the same for every solvent site
!> - Lattice Jacobian J_n x = x - k2r[(P^T - I) r2k(D_n x)] with the bulk response tabulated
!>   at the lattice's own k-points and D_n = R c', frozen per Newton step (`refresh`)
!> - Nested level (coarse_n = n/2): D averaged over 2x2x2 cells, the same Jacobian on the
!>   coarse lattice; B(q) = q + P_2(y - R_2 q) with y from bounded GMRES on J_coarse,
!>   R_2 the block average and P_2 the injection
!> - Fine level: bounded flexible GMRES on J_n preconditioned by B (B = I without nesting)
!> - Molecular action: M^{-1} v = v + P(z - R v) with z ~ J_n^{-1} R v; v when the lattice
!>   solve is not accepted. Molecular modes R does not see pass unchanged
!> - The inner solves make M^{-1} depend on v: use it as a flexible right preconditioner
!> - Storage stays inside: `setup` places the lattices, has a response source tabulate the
!>   bulk response at their k-points and builds the operators; callers query sizes and work
module moist_model_moz_3d_preconditioner
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_utils_timer, only: cat_preconditioner
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_box_grid
   use moist_math_grid_deposit, only: grid_deposit_type
   use moist_math_grid_layout, only: moist_math_grid_layout_type, layout_cartesian
   use moist_model_moz_3d_equation, only: uv3d_equation_type
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type, new_uv3d_susceptibility
   use moist_model_moz_solvent_vv, only: solvent_vv_type
   use moist_model_moz_3d_preconditioner_krylov, only: bounded_operator_type, bounded_gmres_type, &
      & bounded_gmres_result_type
   use moist_model_moz_solver_core_preconditioner, only: moist_model_moz_solver_preconditioner_type
   implicit none(type, external)
   private

   public :: uv3d_aux_preconditioner_type, uv3d_aux_statistics_type, uv3d_aux_level_type
   public :: uv3d_aux_response_source_type, uv3d_aux_solvent_source_type, new_uv3d_aux_level

   !> Bulk response tables at the k-points of a lattice
   type, abstract :: uv3d_aux_response_source_type
   contains
      !> Response at the k-points of `grid`
      procedure(tabulate_interface), deferred :: tabulate
   end type uv3d_aux_response_source_type

   abstract interface
      !> Response at the k-points of a lattice
      !>
      !> @param[in]  self      source
      !> @param[in]  grid      lattice with its k-points
      !> @param[out] response  response table
      !> @param[out] error     tabulation failure
      subroutine tabulate_interface(self, grid, response, error)
         import :: uv3d_aux_response_source_type, moist_math_grid_3d_type, uv3d_susceptibility_type, error_type
         implicit none(type, external)
         !> Source
         class(uv3d_aux_response_source_type), intent(in) :: self
         !> Lattice with its k-points
         class(moist_math_grid_3d_type), intent(in) :: grid
         !> Response table
         type(uv3d_susceptibility_type), intent(out) :: response
         !> Tabulation failure
         type(error_type), allocatable, intent(out) :: error
      end subroutine tabulate_interface
   end interface

   !> The solvent's own response (`new_uv3d_susceptibility`), as the molecular grid uses
   type, extends(uv3d_aux_response_source_type) :: uv3d_aux_solvent_source_type
      !> Solved bulk solvent, borrowed for the setup
      type(solvent_vv_type), pointer :: solvent => null()
      !> Susceptibility mode: oz, direct or linear
      character(len=:), allocatable :: mode
   contains
      !> Response of the solved bulk solvent at the k-points of a lattice
      procedure :: tabulate => solvent_tabulate
   end type uv3d_aux_solvent_source_type

   !> Work of the preconditioner; index 1 the fine lattice, 2 the coarse one
   type :: uv3d_aux_statistics_type
      !> Molecular applications
      integer(int64) :: napply = 0
      !> Bounded solves per lattice
      integer(int64) :: nsolve(2) = 0
      !> Arnoldi steps per lattice
      integer(int64) :: niteration(2) = 0
      !> Lattice Jacobian products per lattice
      integer(int64) :: nproduct(2) = 0
      !> Solves that returned the identity action, per lattice
      integer(int64) :: nfallback(2) = 0
   end type uv3d_aux_statistics_type

   !> One lattice: J x = x - T'(D x)
   type, extends(bounded_operator_type) :: uv3d_aux_level_type
      private
      !> Cells per axis
      integer :: n = 0
      !> Cells and solvent sites
      integer :: ncell = 0, nv = 0
      !> OpenMP team size
      integer :: nthreads = 1
      !> Lattice as a Cartesian grid
      type(moist_math_grid_3d_cartesian_type) :: grid
      !> Map derivative on the lattice with its own response table
      type(uv3d_equation_type) :: equation
      !> Frozen closure derivative D (ncell, nv)
      real(wp), allocatable :: d(:, :)
      !> Scratch (ncell, nv)
      real(wp), allocatable :: work(:, :)
   contains
      !> J x = x - T'(D x) on this lattice
      procedure :: matvec => level_matvec
      !> Freeze a closure derivative on the cells (ncell, nv)
      procedure :: set_derivative => level_set_derivative
   end type uv3d_aux_level_type

   !> Fine lattice, preconditioned by the coarse correction when nested
   type, extends(uv3d_aux_level_type) :: uv3d_aux_fine_level_type
      private
      !> Whether the coarse lattice is used
      logical :: nested = .false.
      !> Residual reduction of the coarse solves
      real(wp) :: rtol = 0.1_wp
      !> Coarse lattice
      type(uv3d_aux_level_type) :: coarse
      !> Work space of the coarse solves
      type(bounded_gmres_type) :: coarse_solver
      !> Coarse right-hand side and solution (ncell_coarse*nv)
      real(wp), allocatable :: bc(:), yc(:)
      !> Coarse solves, steps, products and fallbacks
      integer(int64) :: nsolve = 0, niteration = 0, nproduct = 0, nfallback = 0
   contains
      !> Coarse correction of a fine-lattice vector; the vector itself without nesting
      procedure :: precondition => fine_precondition
   end type uv3d_aux_fine_level_type

   !> Cartesian auxiliary-space preconditioner
   type, extends(moist_model_moz_solver_preconditioner_type) :: uv3d_aux_preconditioner_type
      private
      !> Run context for timers, borrowed; may stay null
      type(moist_context_type), pointer :: ctx => null()
      !> Requested fine cells per axis
      integer :: aux_n = 64
      !> Residual reduction of the lattice solves
      real(wp) :: rtol = 0.1_wp
      !> Molecular points and solvent sites
      integer :: ngrid = 0, nv = 0
      !> OpenMP team size
      integer :: nthreads = 1
      !> Lattice side, bohr
      real(wp) :: lattice_side = 0.0_wp
      !> Transfers between the molecular points and the fine lattice
      type(grid_deposit_type) :: deposit
      !> Fine lattice and, inside it, the coarse one
      type(uv3d_aux_fine_level_type) :: fine
      !> Work space of the fine solves
      type(bounded_gmres_type) :: solver
      !> Fine right-hand side and solution (ncell*nv)
      real(wp), allocatable :: b(:), z(:)
      !> Molecular correction (ngrid, nv)
      real(wp), allocatable :: correction(:, :)
      !> Lattices placed and operators built
      logical :: placed = .false., has_operators = .false.
      !> Lent to a solve, and refreshed since
      logical :: lent = .false., refreshed = .false.
      !> Fine-lattice work
      type(uv3d_aux_statistics_type) :: fine_stats
   contains
      !> Place the lattices on a grid, tabulate their responses and build the operators
      procedure :: setup => aux_setup
      !> Whether the operators are built
      procedure :: built => aux_built
      !> Whether it was built for a solve with these points and sites
      procedure :: check => aux_check
      !> Whether the preconditioner is lent to a solve
      procedure :: attached => aux_attached
      !> Whether the coarse lattice is used
      procedure :: nested => aux_nested
      !> Cells per axis of the fine (1) or coarse (2) lattice; 0 without it
      procedure :: lattice_size => aux_lattice_size
      !> Lattice side, bohr
      procedure :: side => aux_side
      !> Freeze the closure derivative of the current Newton step
      procedure :: refresh => aux_refresh
      !> Lend the preconditioner to a solve
      procedure :: attach => aux_attach
      !> End the solve the preconditioner was lent to
      procedure :: detach => aux_detach
      !> Work so far
      procedure :: statistics => aux_statistics
      !> Whether the operators are built and hold this Newton step's derivative
      procedure :: ready => aux_ready
      !> M^{-1} v = v + P(z - R v), z ~ J^{-1} R v
      procedure :: apply => aux_apply
      !> Name and lattice sizes for the solver notes
      procedure :: describe => aux_describe
      !> The inner solves make the action depend on the vector
      procedure :: flexible => aux_flexible
   end type uv3d_aux_preconditioner_type

contains

   !> Place the lattices on an updated grid, tabulate their responses and build the operators
   !>
   !> @param[inout] self      preconditioner; a previous setup is dropped
   !> @param[in]    grid      updated grid with Fourier modes
   !> @param[in]    source    bulk response at the lattices' k-points
   !> @param[in]    beta      inverse thermal energy, 1/Hartree
   !> @param[in]    aux_n     requested fine cells per axis, even, at least 4
   !> @param[in]    coarse_n  aux_n, or aux_n/2 for a nested coarse lattice
   !> @param[in]    max_iter  most Arnoldi steps per lattice solve, at least 1
   !> @param[in]    rtol      residual reduction of the lattice solves, in (0, 1)
   !> @param[in]    nthreads  OpenMP team size
   !> @param[out]   error     invalid settings or weights, grid not updated or without Fourier
   !>                          modes, fewer than 4 usable cells per axis, points the lattice
   !>                          box does not enclose, a response source without solvent,
   !>                          tables, or lattice equation setup
   !> @param[in]    ctx       run context for timers, borrowed
   subroutine aux_setup(self, grid, source, beta, aux_n, coarse_n, max_iter, rtol, nthreads, error, ctx)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(inout) :: self
      !> Updated grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Response source
      class(uv3d_aux_response_source_type), intent(in) :: source
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Requested fine and coarse cells per axis
      integer, intent(in) :: aux_n, coarse_n
      !> Most Arnoldi steps
      integer, intent(in) :: max_iter
      !> Residual reduction
      real(wp), intent(in) :: rtol
      !> OpenMP team size
      integer, intent(in) :: nthreads
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Run context
      type(moist_context_type), intent(in), target, optional :: ctx
      type(uv3d_susceptibility_type) :: fine, coarse

      self%placed = .false.
      self%has_operators = .false.
      self%refreshed = .false.
      if (max_iter < 1 .or. .not. (rtol > 0.0_wp .and. rtol < 1.0_wp)) then
         call fatal_error(error, "cartesian-aux: inner steps must be at least 1 and the tolerance in (0, 1)")
         return
      end if
      call place_lattices(self, grid, aux_n, coarse_n, nthreads, error)
      if (allocated(error)) return
      call source%tabulate(self%fine%grid, fine, error)
      if (.not. allocated(error) .and. self%fine%nested) call source%tabulate(self%fine%coarse%grid, coarse, error)
      if (allocated(error)) return
      if (self%fine%nested) then
         call build_operators(self, fine, coarse, beta, max_iter, rtol, error, ctx)
      else
         call build_operators(self, fine, fine, beta, max_iter, rtol, error, ctx)
      end if
   end subroutine aux_setup

   !> Place the lattices on an updated grid and build the transfers
   !>
   !> @param[inout] self      preconditioner
   !> @param[in]    grid      updated grid with Fourier modes
   !> @param[in]    aux_n     requested fine cells per axis, even, at least 4
   !> @param[in]    coarse_n  aux_n, or aux_n/2 for a nested coarse lattice
   !> @param[in]    nthreads  OpenMP team size
   !> @param[out]   error     grid not updated or without Fourier modes, invalid sizes or weights
   subroutine place_lattices(self, grid, aux_n, coarse_n, nthreads, error)
      !> Preconditioner
      type(uv3d_aux_preconditioner_type), intent(inout) :: self
      !> Updated grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Requested fine cells per axis
      integer, intent(in) :: aux_n
      !> Requested coarse cells per axis
      integer, intent(in) :: coarse_n
      !> OpenMP team size
      integer, intent(in) :: nthreads
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: two_pi = 8.0_wp*atan(1.0_wp)
      type(moist_math_grid_layout_type) :: layout
      real(wp) :: spacing
      integer :: nmax, n

      if (aux_n < 4 .or. mod(aux_n, 2) /= 0) then
         call fatal_error(error, "cartesian-aux: auxiliary cells per axis must be even and at least 4")
         return
      end if
      if (coarse_n /= aux_n .and. 2*coarse_n /= aux_n) then
         call fatal_error(error, "cartesian-aux: coarse cells per axis must equal the auxiliary ones or half of them")
         return
      end if
      if (grid%ngrid < 1 .or. .not. allocated(grid%xyz) .or. .not. allocated(grid%w)) then
         call fatal_error(error, "cartesian-aux needs an updated grid")
         return
      end if
      call grid%get_layout(layout)
      if (layout%kind == layout_cartesian) then
         spacing = layout%spacing
         self%lattice_side = real(maxval(layout%shape), wp)*spacing
      else if (allocated(layout%kspacing) .and. allocated(layout%kmodes)) then
         self%lattice_side = two_pi/minval(layout%kspacing)
         spacing = two_pi/(real(layout%kmodes(1), wp)*layout%kspacing(1))
      else
         call fatal_error(error, "cartesian-aux needs a grid with Fourier modes")
         return
      end if
      nmax = 2*int(self%lattice_side/spacing/2.0_wp + 1.0e-9_wp)
      n = min(aux_n, nmax)
      if (n < 4) then
         call fatal_error(error, "cartesian-aux: the grid's Fourier modes resolve fewer than 4 cells per axis")
         return
      end if
      self%aux_n = aux_n
      self%nthreads = max(1, nthreads)
      self%ngrid = grid%ngrid
      call self%deposit%build(grid%xyz, grid%w, n, self%lattice_side, error, self%nthreads)
      if (allocated(error)) then
         error%message = "cartesian-aux: "//error%message
         return
      end if
      self%fine%nested = coarse_n < aux_n
      call place_level(self%fine%uv3d_aux_level_type, n, self%lattice_side/real(n, wp), self%deposit%origin, &
         & self%nthreads, error)
      if (allocated(error)) return
      if (self%fine%nested) then
         call place_level(self%fine%coarse, n/2, 2.0_wp*self%lattice_side/real(n, wp), self%deposit%origin, &
            & self%nthreads, error)
         if (allocated(error)) return
      end if
      self%placed = .true.
   end subroutine place_lattices

   !> A lattice as a Cartesian box grid
   !>
   !> @param[inout] level     lattice
   !> @param[in]    n         cells per axis
   !> @param[in]    spacing   cell size, bohr
   !> @param[in]    origin    box corner, bohr
   !> @param[in]    nthreads  OpenMP team size
   !> @param[out]   error     invalid box
   subroutine place_level(level, n, spacing, origin, nthreads, error)
      !> Lattice
      type(uv3d_aux_level_type), intent(inout) :: level
      !> Cells per axis
      integer, intent(in) :: n
      !> Cell size
      real(wp), intent(in) :: spacing
      !> Box corner
      real(wp), intent(in) :: origin(3)
      !> OpenMP team size
      integer, intent(in) :: nthreads
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      level%n = n
      level%ncell = n**3
      level%nthreads = max(1, nthreads)
      call new_cartesian_box_grid(level%grid, n, spacing, origin, error)
   end subroutine place_level

   !> One lattice operator J x = x - T'(D x), for use outside the preconditioner
   !>
   !> - D starts at zero; `set_derivative` freezes another
   !>
   !> @param[out] level     lattice operator
   !> @param[in]  n         cells per axis
   !> @param[in]  spacing   cell size, bohr
   !> @param[in]  origin    box corner, bohr
   !> @param[in]  source    bulk response at the lattice's k-points
   !> @param[in]  beta      inverse thermal energy, 1/Hartree
   !> @param[in]  nthreads  OpenMP team size
   !> @param[out] error     invalid box, table or transform failure
   subroutine new_uv3d_aux_level(level, n, spacing, origin, source, beta, nthreads, error)
      !> Lattice operator
      type(uv3d_aux_level_type), intent(out) :: level
      !> Cells per axis
      integer, intent(in) :: n
      !> Cell size
      real(wp), intent(in) :: spacing
      !> Box corner
      real(wp), intent(in) :: origin(3)
      !> Response source
      class(uv3d_aux_response_source_type), intent(in) :: source
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> OpenMP team size
      integer, intent(in) :: nthreads
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(uv3d_susceptibility_type) :: response

      call place_level(level, n, spacing, origin, nthreads, error)
      if (.not. allocated(error)) call source%tabulate(level%grid, response, error)
      if (.not. allocated(error)) call build_level(level, response, beta, level%nthreads, error)
   end subroutine new_uv3d_aux_level

   !> Freeze a closure derivative on the cells of a lattice
   !>
   !> @param[inout] self   lattice operator
   !> @param[in]    d      derivative (ncell, nv)
   !> @param[out]   error  operator not built or wrong shape
   subroutine level_set_derivative(self, d, error)
      !> Lattice operator
      class(uv3d_aux_level_type), intent(inout) :: self
      !> Derivative
      real(wp), intent(in) :: d(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(self%d)) then
         call fatal_error(error, "cartesian-aux: the lattice operator is not built")
      else if (size(d, 1) /= self%ncell .or. size(d, 2) /= self%nv) then
         call fatal_error(error, "cartesian-aux: the derivative does not match the lattice")
      else
         self%d = d
      end if
   end subroutine level_set_derivative

   !> Response of the solved bulk solvent at the k-points of a lattice
   !>
   !> @param[in]  self      source with a borrowed solvent
   !> @param[in]  grid      lattice
   !> @param[out] response  response table
   !> @param[out] error     no solvent, or the tabulation failed
   subroutine solvent_tabulate(self, grid, response, error)
      !> Source
      class(uv3d_aux_solvent_source_type), intent(in) :: self
      !> Lattice
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Response table
      type(uv3d_susceptibility_type), intent(out) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. associated(self%solvent) .or. .not. allocated(self%mode)) then
         call fatal_error(error, "cartesian-aux: the response source has no solvent")
         return
      end if
      call new_uv3d_susceptibility(response, grid, self%solvent, self%mode, error)
   end subroutine solvent_tabulate

   !> Bind the response tables of the lattices and build their operators
   !>
   !> @param[inout] self      placed preconditioner
   !> @param[in]    fine      response at the fine lattice's k-points
   !> @param[in]    coarse    response at the coarse lattice's k-points (ignored unless nested)
   !> @param[in]    beta      inverse thermal energy, 1/Hartree
   !> @param[in]    max_iter  most Arnoldi steps per lattice solve, at least 1
   !> @param[in]    rtol      residual reduction of the lattice solves, in (0, 1)
   !> @param[out]   error     lattice equation setup, or responses of different solvent sites
   !> @param[in]    ctx       run context for timers, borrowed
   subroutine build_operators(self, fine, coarse, beta, max_iter, rtol, error, ctx)
      !> Preconditioner
      type(uv3d_aux_preconditioner_type), intent(inout) :: self
      !> Fine response
      type(uv3d_susceptibility_type), intent(in) :: fine
      !> Coarse response
      type(uv3d_susceptibility_type), intent(in) :: coarse
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Most Arnoldi steps
      integer, intent(in) :: max_iter
      !> Residual reduction
      real(wp), intent(in) :: rtol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Run context
      type(moist_context_type), intent(in), target, optional :: ctx

      self%ctx => null()
      if (present(ctx)) self%ctx => ctx
      if (associated(self%ctx)) call self%ctx%timer%set_category_summary(cat_preconditioner, .true.)
      self%rtol = rtol
      self%nv = fine%nv
      self%fine_stats = uv3d_aux_statistics_type()
      self%fine%nsolve = 0
      self%fine%niteration = 0
      self%fine%nproduct = 0
      self%fine%nfallback = 0
      call build_level(self%fine%uv3d_aux_level_type, fine, beta, self%nthreads, error)
      if (allocated(error)) return
      call self%solver%setup(self%fine%ncell*self%nv, max_iter, self%nthreads)
      if (self%fine%nested) then
         if (coarse%nv /= self%nv) then
            call fatal_error(error, "cartesian-aux: the lattice responses have different solvent sites")
            return
         end if
         call build_level(self%fine%coarse, coarse, beta, self%nthreads, error)
         if (allocated(error)) return
         call self%fine%coarse_solver%setup(self%fine%coarse%ncell*self%nv, max_iter, self%nthreads)
         self%fine%rtol = rtol
         if (allocated(self%fine%bc)) deallocate (self%fine%bc)
         if (allocated(self%fine%yc)) deallocate (self%fine%yc)
         allocate (self%fine%bc(self%fine%coarse%ncell*self%nv), self%fine%yc(self%fine%coarse%ncell*self%nv))
      end if
      if (allocated(self%b)) deallocate (self%b)
      if (allocated(self%z)) deallocate (self%z)
      if (allocated(self%correction)) deallocate (self%correction)
      allocate (self%b(self%fine%ncell*self%nv), self%z(self%fine%ncell*self%nv))
      allocate (self%correction(self%ngrid, self%nv))
      self%has_operators = .true.
   end subroutine build_operators

   !> Set up the equation of one lattice with a zero long-range table
   !>
   !> @param[inout] level     lattice with its grid
   !> @param[in]    response  response at its k-points
   !> @param[in]    beta      inverse thermal energy
   !> @param[in]    nthreads  OpenMP team size
   !> @param[out]   error     invalid table or transform failure
   subroutine build_level(level, response, beta, nthreads, error)
      !> Lattice
      type(uv3d_aux_level_type), intent(inout) :: level
      !> Response
      type(uv3d_susceptibility_type), intent(in) :: response
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> OpenMP team size
      integer, intent(in) :: nthreads
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      complex(wp), allocatable :: uk(:, :)

      level%nv = response%nv
      allocate (uk(level%grid%npts_k, max(1, response%nv)), source=(0.0_wp, 0.0_wp))
      level%equation%nthreads = nthreads
      call level%equation%setup(level%grid, response, uk, beta, error)
      if (allocated(error)) then
         error%message = "cartesian-aux: "//error%message
         return
      end if
      if (allocated(level%d)) deallocate (level%d)
      if (allocated(level%work)) deallocate (level%work)
      allocate (level%d(level%ncell, level%nv), source=0.0_wp)
      allocate (level%work(level%ncell, level%nv))
   end subroutine build_level

   !> Freeze the closure derivative of the current Newton step: D = R c', and its
   !> 2x2x2 block average on the coarse lattice
   !>
   !> @param[inout] self  built preconditioner
   !> @param[in]    dcdt  closure derivative on the molecular points (ngrid, nv)
   subroutine aux_refresh(self, dcdt)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(inout) :: self
      !> Closure derivative
      real(wp), intent(in) :: dcdt(:, :)

      if (.not. self%has_operators) return
      if (associated(self%ctx)) call self%ctx%timer%start("Preconditioner", category=cat_preconditioner)
      call self%deposit%restrict(dcdt, self%fine%d, self%nthreads)
      if (self%fine%nested) call block_average(self%fine%n, self%nv, self%fine%d, self%fine%coarse%d, self%nthreads)
      if (associated(self%ctx)) call self%ctx%timer%stop("Preconditioner")
      self%refreshed = .true.
   end subroutine aux_refresh

   !> Lend the preconditioner to a solve; it is ready after the first refresh
   !>
   !> @param[inout] self  preconditioner
   subroutine aux_attach(self)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(inout) :: self
      self%lent = .true.
      self%refreshed = .false.
   end subroutine aux_attach

   !> End the solve the preconditioner was lent to
   !>
   !> @param[inout] self  preconditioner
   subroutine aux_detach(self)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(inout) :: self
      self%lent = .false.
      self%refreshed = .false.
   end subroutine aux_detach

   !> Whether the operators are built
   !>
   !> @param[in] self  preconditioner
   pure logical function aux_built(self)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(in) :: self
      aux_built = self%has_operators
   end function aux_built

   !> Whether the preconditioner was built for a solve with these points and sites
   !>
   !> - A guard against lending it to another problem; a stale lattice would only make the
   !>   preconditioner weaker, the molecular solve stays the same
   !>
   !> @param[in]  self   preconditioner
   !> @param[in]  ngrid  points of the solve
   !> @param[in]  nv     solvent sites of the solve
   !> @param[out] error  not built, or built for other points or sites
   subroutine aux_check(self, ngrid, nv, error)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(in) :: self
      !> Points
      integer, intent(in) :: ngrid
      !> Solvent sites
      integer, intent(in) :: nv
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%has_operators) then
         call fatal_error(error, "cartesian-aux: the preconditioner is not built")
      else if (self%ngrid /= ngrid .or. self%nv /= nv) then
         call fatal_error(error, "cartesian-aux: the preconditioner was built for other points or solvent sites")
      end if
   end subroutine aux_check

   !> Whether the preconditioner is lent to a solve
   !>
   !> @param[in] self  preconditioner
   pure logical function aux_attached(self)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(in) :: self
      aux_attached = self%lent
   end function aux_attached

   !> Whether the coarse lattice is used
   !>
   !> @param[in] self  preconditioner
   pure logical function aux_nested(self)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(in) :: self
      aux_nested = self%placed .and. self%fine%nested
   end function aux_nested

   !> Cells per axis of the fine (1) or coarse (2) lattice; 0 when it is not placed
   !>
   !> @param[in] self   preconditioner
   !> @param[in] level  1 fine, 2 coarse
   pure integer function aux_lattice_size(self, level) result(n)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(in) :: self
      !> Lattice
      integer, intent(in) :: level
      n = 0
      if (.not. self%placed) return
      if (level == 1) n = self%fine%n
      if (level == 2 .and. self%fine%nested) n = self%fine%coarse%n
   end function aux_lattice_size

   !> Lattice side, bohr; 0 before the lattices are placed
   !>
   !> @param[in] self  preconditioner
   pure real(wp) function aux_side(self)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(in) :: self
      aux_side = 0.0_wp
      if (self%placed) aux_side = self%lattice_side
   end function aux_side

   !> Work since the build
   !>
   !> @param[in] self  preconditioner
   function aux_statistics(self) result(stats)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(in) :: self
      !> Work
      type(uv3d_aux_statistics_type) :: stats
      stats = self%fine_stats
      stats%nsolve(2) = self%fine%nsolve
      stats%niteration(2) = self%fine%niteration
      stats%nproduct(2) = self%fine%nproduct
      stats%nfallback(2) = self%fine%nfallback
   end function aux_statistics

   !> Whether the operators are built and hold this Newton step's derivative
   !>
   !> @param[in] self  preconditioner
   !> @param[in] rmse  previous residual RMS, unused
   logical function aux_ready(self, rmse) result(usable)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(in) :: self
      !> Previous residual RMS
      real(wp), intent(in) :: rmse
      usable = self%has_operators .and. self%refreshed
   end function aux_ready

   !> The action depends on the vector through the bounded inner solves
   !>
   !> @param[in] self  preconditioner
   pure logical function aux_flexible(self) result(flexible)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(in) :: self
      flexible = .true.
   end function aux_flexible

   !> Name and lattice sizes for the solver notes
   !>
   !> @param[in] self  preconditioner
   !> @param[in] rmse  previous residual RMS, unused
   function aux_describe(self, rmse) result(text)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(in) :: self
      !> Previous residual RMS
      real(wp), intent(in) :: rmse
      !> Status text
      character(len=:), allocatable :: text
      character(len=32) :: sizes

      if (.not. self%has_operators) then
         text = "cartesian-aux, not built"
         return
      end if
      if (self%fine%nested) then
         write (sizes, "(i0,'^3/',i0,'^3')") self%fine%n, self%fine%coarse%n
      else
         write (sizes, "(i0,'^3')") self%fine%n
      end if
      text = "cartesian-aux "//trim(sizes)
      if (self%fine%n < self%aux_n) text = text//", capped at the grid's Fourier modes"
      if (.not. self%refreshed) text = text//", waiting for a linearization"
   end function aux_describe

   !> M^{-1} v = v + P(z - R v), z ~ J^{-1} R v; v unchanged when the lattice solve fails
   !>
   !> @param[inout] self   preconditioner and its scratch
   !> @param[inout] vflat  molecular Krylov vector, i + (site - 1)*ngrid
   subroutine aux_apply(self, vflat)
      !> Preconditioner
      class(uv3d_aux_preconditioner_type), intent(inout) :: self
      !> Flat Krylov vector
      real(wp), contiguous, target, intent(inout) :: vflat(:)
      real(wp), pointer :: v(:, :)
      type(bounded_gmres_result_type) :: result

      if (.not. (self%has_operators .and. self%refreshed)) return
      if (size(vflat) /= self%ngrid*self%nv) return
      if (associated(self%ctx)) call self%ctx%timer%start("Preconditioner", category=cat_preconditioner)
      v(1:self%ngrid, 1:self%nv) => vflat
      call restrict_into(self%deposit, v, self%b, self%fine%ncell, self%nv, self%nthreads)
      call self%solver%solve(self%fine, self%b, self%z, self%rtol, result)
      self%fine_stats%napply = self%fine_stats%napply + 1
      self%fine_stats%nsolve(1) = self%fine_stats%nsolve(1) + 1
      self%fine_stats%niteration(1) = self%fine_stats%niteration(1) + result%iterations
      self%fine_stats%nproduct(1) = self%fine_stats%nproduct(1) + result%nmatvec
      if (result%accepted) then
         self%z = self%z - self%b
         call interpolate_from(self%deposit, self%z, self%correction, self%fine%ncell, self%nv, self%nthreads)
         v = v + self%correction
      else
         self%fine_stats%nfallback(1) = self%fine_stats%nfallback(1) + 1
      end if
      if (associated(self%ctx)) call self%ctx%timer%stop("Preconditioner")
   end subroutine aux_apply

   !> R of molecular columns into a flat lattice vector
   !>
   !> @param[in]  deposit   transfers
   !> @param[in]  v         molecular columns (ngrid, nv)
   !> @param[out] b         lattice columns (ncell, nv)
   !> @param[in]  ncell     cells per lattice
   !> @param[in]  nv        solvent sites
   !> @param[in]  nthreads  OpenMP team size
   subroutine restrict_into(deposit, v, b, ncell, nv, nthreads)
      !> Transfers
      type(grid_deposit_type), intent(in) :: deposit
      !> Molecular columns (ngrid, nv)
      real(wp), intent(in) :: v(:, :)
      !> Cells per lattice and solvent sites
      integer, intent(in) :: ncell, nv
      !> Lattice columns
      real(wp), intent(out) :: b(ncell, nv)
      !> OpenMP team size
      integer, intent(in) :: nthreads
      call deposit%restrict(v, b, nthreads)
   end subroutine restrict_into

   !> P of a flat lattice vector into molecular columns
   !>
   !> @param[in]  deposit   transfers
   !> @param[in]  z         lattice columns (ncell, nv)
   !> @param[out] v         molecular columns (ngrid, nv)
   !> @param[in]  ncell     cells per lattice
   !> @param[in]  nv        solvent sites
   !> @param[in]  nthreads  OpenMP team size
   subroutine interpolate_from(deposit, z, v, ncell, nv, nthreads)
      !> Transfers
      type(grid_deposit_type), intent(in) :: deposit
      !> Cells per lattice and solvent sites
      integer, intent(in) :: ncell, nv
      !> Lattice columns
      real(wp), intent(in) :: z(ncell, nv)
      !> Molecular columns (ngrid, nv)
      real(wp), intent(out) :: v(:, :)
      !> OpenMP team size
      integer, intent(in) :: nthreads
      call deposit%interpolate(z, v, nthreads)
   end subroutine interpolate_from

   !> J x = x - T'(D x) on one lattice
   !>
   !> @param[inout] self  lattice and its scratch
   !> @param[in]    x     flat lattice vector, cell + (site - 1)*ncell
   !> @param[out]   y     product
   !> @param[out]   ok    false when the transform failed
   subroutine level_matvec(self, x, y, ok)
      !> Lattice
      class(uv3d_aux_level_type), intent(inout) :: self
      !> Vector
      real(wp), contiguous, intent(in) :: x(:)
      !> Product
      real(wp), contiguous, intent(out) :: y(:)
      !> Whether the product was formed
      logical, intent(out) :: ok
      call lattice_jacobian(self, x, y, self%ncell, self%nv, ok)
   end subroutine level_matvec

   !> J x = x - T'(D x) on lattice columns
   !>
   !> @param[inout] level  lattice and its scratch
   !> @param[in]    x      vector (ncell, nv)
   !> @param[out]   y      product (ncell, nv)
   !> @param[in]    ncell  cells
   !> @param[in]    nv     columns
   !> @param[out]   ok     false when the transform failed
   subroutine lattice_jacobian(level, x, y, ncell, nv, ok)
      !> Lattice
      type(uv3d_aux_level_type), intent(inout) :: level
      !> Cells and columns
      integer, intent(in) :: ncell, nv
      !> Vector
      real(wp), intent(in) :: x(ncell, nv)
      !> Product
      real(wp), intent(out), target :: y(ncell, nv)
      !> Whether the product was formed
      logical, intent(out) :: ok
      type(error_type), allocatable :: error
      integer :: i, s

      do s = 1, nv
         !$omp parallel do schedule(static) num_threads(level%nthreads)
         do i = 1, ncell
            level%work(i, s) = level%d(i, s)*x(i, s)
         end do
         !$omp end parallel do
      end do
      call level%equation%jvp(level%work, y, error)
      ok = .not. allocated(error)
      if (.not. ok) return
      do s = 1, nv
         !$omp parallel do schedule(static) num_threads(level%nthreads)
         do i = 1, ncell
            y(i, s) = x(i, s) - y(i, s)
         end do
         !$omp end parallel do
      end do
   end subroutine lattice_jacobian

   !> Coarse correction B(q) = q + P_2(y - R_2 q), y ~ J_coarse^{-1} R_2 q; q without nesting
   !>
   !> @param[inout] self  fine lattice with the coarse one
   !> @param[in]    q     flat fine vector
   !> @param[out]   z     preconditioned vector
   !> @param[out]   ok    always true; a failed coarse solve gives z = q
   subroutine fine_precondition(self, q, z, ok)
      !> Fine lattice
      class(uv3d_aux_fine_level_type), intent(inout) :: self
      !> Vector
      real(wp), contiguous, intent(in) :: q(:)
      !> Preconditioned vector
      real(wp), contiguous, intent(out) :: z(:)
      !> Whether z was formed
      logical, intent(out) :: ok
      type(bounded_gmres_result_type) :: result

      ok = .true.
      z = q
      if (.not. self%nested) return
      call block_average(self%n, self%nv, q, self%bc, self%nthreads)
      call self%coarse_solver%solve(self%coarse, self%bc, self%yc, self%rtol, result)
      self%nsolve = self%nsolve + 1
      self%niteration = self%niteration + result%iterations
      self%nproduct = self%nproduct + result%nmatvec
      if (.not. result%accepted) then
         self%nfallback = self%nfallback + 1
         return
      end if
      self%yc = self%yc - self%bc
      call inject_add(self%n, self%nv, self%yc, z, self%nthreads)
   end subroutine fine_precondition

   !> Average over 2x2x2 blocks: coarse(I) = mean of fine cells (2I-1 .. 2I) per axis
   !>
   !> @param[in]  n         fine cells per axis, even
   !> @param[in]  nv        columns
   !> @param[in]  fine      fine values, x fastest, column after column ((n^3)*nv)
   !> @param[out] coarse    coarse values (((n/2)^3)*nv)
   !> @param[in]  nthreads  OpenMP team size
   subroutine block_average(n, nv, fine, coarse, nthreads)
      !> Fine cells per axis
      integer, intent(in) :: n
      !> Columns
      integer, intent(in) :: nv
      !> Fine values
      real(wp), intent(in) :: fine(n, n, n, nv)
      !> Coarse values
      real(wp), intent(out) :: coarse(n/2, n/2, n/2, nv)
      !> OpenMP team size
      integer, intent(in) :: nthreads
      integer :: i, j, k, s

      do s = 1, nv
         !$omp parallel do schedule(static) num_threads(nthreads) private(i, j)
         do k = 1, n/2
            do j = 1, n/2
               do i = 1, n/2
                  coarse(i, j, k, s) = 0.125_wp*sum(fine(2*i - 1:2*i, 2*j - 1:2*j, 2*k - 1:2*k, s))
               end do
            end do
         end do
         !$omp end parallel do
      end do
   end subroutine block_average

   !> Add the coarse values to every fine cell of their block (injection)
   !>
   !> @param[in]    n         fine cells per axis, even
   !> @param[in]    nv        columns
   !> @param[in]    coarse    coarse values (((n/2)^3)*nv)
   !> @param[inout] fine      fine values ((n^3)*nv)
   !> @param[in]    nthreads  OpenMP team size
   subroutine inject_add(n, nv, coarse, fine, nthreads)
      !> Fine cells per axis
      integer, intent(in) :: n
      !> Columns
      integer, intent(in) :: nv
      !> Coarse values
      real(wp), intent(in) :: coarse(n/2, n/2, n/2, nv)
      !> Fine values
      real(wp), intent(inout) :: fine(n, n, n, nv)
      !> OpenMP team size
      integer, intent(in) :: nthreads
      integer :: i, j, k, s

      do s = 1, nv
         !$omp parallel do schedule(static) num_threads(nthreads) private(i, j)
         do k = 1, n
            do j = 1, n
               do i = 1, n
                  fine(i, j, k, s) = fine(i, j, k, s) + coarse((i + 1)/2, (j + 1)/2, (k + 1)/2, s)
               end do
            end do
         end do
         !$omp end parallel do
      end do
   end subroutine inject_add

end module moist_model_moz_3d_preconditioner
