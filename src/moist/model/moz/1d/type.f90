!> One-dimensional MOZ model state
!>
!> - Owned radial grid pair `rgrid`, `kgrid`, independent of the solute structure
!>   and of the grid of the solvent's VV solve
!> - Solute potential `potential` (public) and a private copy of the updated bulk solvent
!> - Settings of the UV solve and, nested as `solvent`, of the VV solve of the
!>   solvent copy; `prepare_solvent` solves it unless its state and residual match
!> - Energy: 1D UV solve on the model's grid with the solvent susceptibility
!>   of the VV solution there; the tables come from
!>   `potential%compute(rgrid, kgrid, ..., solvent=..., coupling=...)`
!> - Solvation free energy and partial molar volume of every energy
!> - Energy and response share one UV solve; tables equal to those of the kept
!>   solution at this structure reuse it; changed tables on the same grids start
!>   from its converged short-range iterate
!> - Fixed solvent response and solute W are cached by bulk identity, geometry
!>   and reciprocal grid, independently of the current UV solve and host tables
!> - Response: derivative of the reported energy with respect to the host data of
!>   the solute terms, from one transposed solve and the potential adjoint; the
!>   transposed solution is kept as the start of the next response
!> - Nuclear gradient: intramolecular W(k) pullback from the same transposed solve,
!>   at fixed radial grids, site parameters and bulk solvent; host derivatives returned
!> - Coupling scope of solute term i: i
!> - Monopole term owns its host-charge requirements
!> - The model does not watch its potential: update it again after a change
module moist_model_moz_1d_type
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_convert, only: autokcal, autoaa
   use moist_math_packed_sym, only: packed_index_lut
   use moist_context, only: moist_context_type
   use moist_model_type, only: solvation_model_type, invalidate_model
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_channels_fields, only: field_query_type
   use moist_channels_coupling, only: coupling_type, coupling_begin_registration, moist_phase_energy, &
      & moist_phase_response, moist_phase_gradient, coupling_snapshot, coupling_check_mandatory
   use moist_channels_response, only: response_type, response_clear
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_potential_kernel_table_1d, only: require_radial_grids
   use moist_model_moz_solvent_vv, only: solvent_vv_type, solvent_vv_thermodynamics_type
   use moist_model_moz_1d_thermodynamics, only: model_moz_1d_thermodynamics_type, new_moz_1d_thermodynamics
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_1d_operators, only: uv_operators_type, same_grid => same_radial_grid
   use moist_model_moz_1d_driver, only: uv_problem_type, new_uv_problem
   use moist_model_moz_solvent_vv_kernel, only: vv_intramolecular_adjoint
   implicit none(type, external)
   private
   public :: model_moz_1d_type, new_moz_1d_model

   !> Direct 1D MOZ model; energy, host response and nuclear gradient from the UV solve
   type, extends(solvation_model_type) :: model_moz_1d_type
      !> Solute potential; add terms before `update`
      type(moz_potential_type) :: potential
      !> Owned real-space radial grid, bohr
      type(moist_math_grid_radial_type) :: rgrid
      !> Owned reciprocal radial grid, 1/bohr
      type(moist_math_grid_radial_type) :: kgrid
      !> Copy of the updated bulk solvent
      type(solvent_vv_type), allocatable, private :: solvent
      !> Settings of the UV solve and of the VV solve of the solvent copy
      type(model_moz_1d_parameters_type), private :: params
      !> Solute coordinates of the latest update, bohr (3, nu)
      real(wp), allocatable, private :: xyz(:, :)
      !> Whether the correlations belong to the latest energy at this structure
      logical, private :: solved = .false.
      !> Bulk identity of the kept UV solution
      integer(int64), private :: bulk_generation = -1_int64
      !> Physical c(r), t(r) = h - c and h(r) of the latest UV solve, (npts, nu nv)
      real(wp), allocatable, private :: cr(:, :), tr(:, :), hr(:, :)
      !> Free energy and partial molar volume of the latest UV solve
      type(model_moz_1d_thermodynamics_type), private :: thermo
      !> Latest UV problem, kept for the response and for reuse
      type(uv_problem_type), allocatable, private :: problem
      !> Fixed operators survive host updates and UV problem reconstruction
      type(uv_operators_type), private :: operators
      !> Potential tables of the kept solution, (npts, nu nv)
      real(wp), allocatable, private :: u_sr(:, :), ur_lr(:, :), uk_lr(:, :)
      !> Transposed solution of the latest response, start of the next one
      real(wp), allocatable, private :: lambda(:, :)
   contains
      final :: destroy_moz_1d_model
      !> Invalidate the UV solution, response guess and host couplings
      procedure :: invalidate => moz_1d_invalidate
      procedure :: update => moz_1d_update
      procedure :: get_energy => moz_1d_energy
      procedure :: get_response => moz_1d_response
      procedure :: get_gradient => moz_1d_gradient
      procedure :: atom_count => moz_1d_atom_count
      !> Declare the coupling requests of every solute term in its own scope
      procedure :: declare_pass => moz_1d_declare_pass
      !> Publish the radial grid as named fields
      procedure :: list_fields => moz_1d_list_fields
      !> Whether the model holds a solvent
      procedure :: has_solvent => moz_1d_has_solvent
      !> Number of solvent sites; 0 before construction
      procedure :: solvent_nsite => moz_1d_solvent_nsite
      !> Solve the solvent copy unless its state and residual satisfy the requested VV settings
      procedure :: prepare_solvent => moz_1d_prepare_solvent
      !> Whether the solvent copy holds a VV solution satisfying the requested settings
      procedure :: solvent_solved => moz_1d_solvent_solved
      !> Copies of the physical correlations of the latest UV solve
      procedure :: correlations => moz_1d_correlations
      !> Free energy and partial molar volume of the latest UV solve
      procedure :: thermodynamics => moz_1d_thermodynamics
      !> UV solution for the tables of the coupling, reused when they are unchanged
      procedure, private :: solve_uv => moz_1d_solve_uv
   end type model_moz_1d_type

contains

   !> Copy the radial grid pair, an updated bulk solvent and the settings
   !>
   !> `prepare_solvent` reuses a solved copy only when the bulk state matches
   !> and its achieved residual meets the requested tolerance; otherwise it solves again.
   !> Without explicit parameters, a solved solvent supplies its accepted VV settings.
   !>
   !> @param[out] self        model
   !> @param[in]  ctx         run context, borrowed for the lifetime of the model
   !> @param[in]  rgrid       real-space radial grid, nodes r > 0; e.g. from `new_uniform_radial_pair`
   !> @param[in]  kgrid       reciprocal radial grid of the pair
   !> @param[in]  solvent     bulk solvent with an updated potential, copied into the model
   !> @param[out] error       unbuilt grid, node at r <= 0, solvent not updated, invalid settings or copy failure
   !> @param[in]  parameters  UV/VV settings; absent parameters retain accepted bulk settings or use defaults
   subroutine new_moz_1d_model(self, ctx, rgrid, kgrid, solvent, error, parameters)
      !> Model
      type(model_moz_1d_type), intent(out) :: self
      !> Run context, borrowed for the lifetime of the model
      type(moist_context_type), intent(in), target :: ctx
      !> Real-space radial grid, bohr
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Reciprocal radial grid, 1/bohr
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Bulk solvent with an updated potential
      type(solvent_vv_type), intent(in) :: solvent
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Settings of the UV and VV solves
      type(model_moz_1d_parameters_type), intent(in), optional :: parameters
      integer :: stat
      call require_radial_grids(rgrid, kgrid, error)
      if (allocated(error)) return
      if (.not. solvent%potential%is_updated()) then
         call fatal_error(error, "Update the 1D VV solvent before constructing the 1D MOZ model")
         return
      end if
      if (present(parameters)) then
         self%params = parameters
      else if (solvent%is_solved()) then
         self%params%solvent = solvent%settings()
      end if
      call self%params%validate(error)
      if (allocated(error)) return
      ! TODO: share the solvent check and copy with the 3D model
      allocate (self%solvent, source=solvent, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to copy the 1D VV solvent")
         return
      end if
      self%rgrid = rgrid
      self%kgrid = kgrid
      self%ctx => ctx
   end subroutine new_moz_1d_model

   !> Invalidate the UV solution, response guess and host couplings
   !>
   !> @param[in,out] self Model
   subroutine moz_1d_invalidate(self)
      !> Model
      class(model_moz_1d_type), intent(inout) :: self
      call invalidate_model(self)
      self%solved = .false.
      if (allocated(self%lambda)) deallocate (self%lambda)
      if (allocated(self%problem)) deallocate (self%problem)
   end subroutine moz_1d_invalidate

   !> Update the solute potential from the structure on the radial grid and print it at verbosity 2
   !>
   !> @param[in,out] self   model
   !> @param[in]     mol    solute structure
   !> @param[out]    error  unconstructed model, empty structure, no terms, term failure or a term refusing the radial grid
   subroutine moz_1d_update(self, mol, error)
      !> Model
      class(model_moz_1d_type), intent(inout) :: self
      !> Solute structure
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%invalidate()
      if (.not. associated(self%ctx) .or. .not. allocated(self%solvent)) then
         call fatal_error(error, "Construct the 1D MOZ model before update")
         return
      end if
      call self%potential%update(mol, self%rgrid, error)
      if (allocated(error)) return
      self%xyz = mol%xyz
      if (self%ctx%verbosity(2)) then
         call self%ctx%message("Solute potential parameters", 2)
         call self%potential%print_table(self%ctx%unit, error)
         if (allocated(error)) return
         call self%ctx%message("", 2)
      end if
      self%updated = .true.
      self%built = .true.
   end subroutine moz_1d_update

   !> Site count of the latest valid structure, from the solute potential
   !>
   !> @param[in] self  model
   function moz_1d_atom_count(self) result(natom)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Site count
      integer :: natom
      natom = self%potential%natom()
   end function moz_1d_atom_count

   !> Declare the coupling requests of every solute term in its own scope
   !>
   !> - Coupling scope of term i: i
   !> - No requests from parameter terms
   !> - Monopole term declares host charges for tables
   !>
   !> @param[in,out] self      updated model
   !> @param[in,out] coupling  coupling to declare
   !> @param[out]    error     failed declaration
   subroutine moz_1d_declare_pass(self, coupling, error)
      !> Updated model
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling to declare
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call coupling_begin_registration(coupling)
      call self%potential%declare(self%rgrid, coupling, error)
      if (allocated(error)) return
      call coupling_snapshot(coupling, ngrid=self%rgrid%npts, natom=self%potential%natom())
   end subroutine moz_1d_declare_pass

   !> Publish the radial grid as named fields, the model's evaluation domain
   !>
   !> - Every solute site carries the same radii
   !>
   !> @param[in]     self   model
   !> @param[in,out] query  field walker
   subroutine moz_1d_list_fields(self, query)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query
      if (.not. allocated(self%rgrid%r)) return
      call query%add_int_value("ngrid", "Number of radial grid points", self%rgrid%npts)
      call query%add_int_value("natom", "Number of solute atoms", self%potential%natom())
      call query%add_real("r", "Radii from each solute site, bohr (ngrid)", self%rgrid%r)
   end subroutine moz_1d_list_fields

   !> Whether the model holds a solvent
   !>
   !> @param[in] self  model
   pure function moz_1d_has_solvent(self) result(has)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Presence
      logical :: has
      has = allocated(self%solvent)
   end function moz_1d_has_solvent

   !> Number of solvent sites; 0 before construction
   !>
   !> @param[in] self  model
   pure function moz_1d_solvent_nsite(self) result(nsite)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Site count
      integer :: nsite
      nsite = 0
      if (allocated(self%solvent)) nsite = self%solvent%nsite()
   end function moz_1d_solvent_nsite

   !> Solve the VV equations of the solvent copy unless it holds a current solution
   !>
   !> Uses the `solvent` settings and their grid; the solvent of the caller is
   !> not touched. Called by the energy; callers may call it early to fail fast
   !>
   !> @param[in,out] self   constructed model
   !> @param[out]    error  unconstructed model or VV solve failure
   subroutine moz_1d_prepare_solvent(self, error)
      !> Model
      class(model_moz_1d_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(self%solvent)) then
         call fatal_error(error, "Construct the 1D MOZ model before preparing its solvent")
         return
      end if
      if (self%solvent%is_solved(self%params%solvent)) return
      call self%solvent%solve(error, self%params%solvent)
   end subroutine moz_1d_prepare_solvent

   !> Whether the solvent copy holds a VV solution satisfying the requested settings
   !>
   !> @param[in] self  model
   pure function moz_1d_solvent_solved(self) result(solved)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Solution flag
      logical :: solved
      solved = .false.
      if (allocated(self%solvent)) solved = self%solvent%is_solved(self%params%solvent)
   end function moz_1d_solvent_solved

   !> Add the solvation free energy of the UV solve
   !>
   !> - UV solution for the tables of the coupling (solve_uv)
   !> - Adds (1/beta) sum_ia rho_a int Phi(h_ia, c_ia) d^3r in Hartree
   !>
   !> @param[in,out] self      updated model
   !> @param[in,out] coupling  coupling staged for the energy, with every mandatory answer
   !> @param[in,out] energy    energy, Hartree; unchanged on error
   !> @param[out]    error     foreign or incomplete coupling, model not updated, VV or UV failure
   subroutine moz_1d_energy(self, coupling, energy, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Energy
      real(wp), intent(inout) :: energy
      !> Error
      type(error_type), allocatable, intent(out) :: error

      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_energy, error)
      if (allocated(error)) return
      call self%solve_uv(coupling, error)
      if (allocated(error)) return
      energy = energy + self%thermo%free_energy
   end subroutine moz_1d_energy

   !> UV solution for the potential tables of the coupling
   !>
   !> - Solves the solvent copy first when needed (prepare_solvent)
   !> - Tables from the solute potential with the host answers on the coupling;
   !>   tables equal to those of the kept solution at this structure reuse it
   !> - Solvent response P_ba = w_ba + rho_b h_ba on the model's k grid
   !>   (susceptibility_at) with the charge-moment correction of its small-k modes
   !> - XRISM UV solve with the model settings, from the previous converged t_s
   !>   at this structure and grid pair when available, otherwise c_s = 0; keeps
   !>   the problem, the tables, the physical correlations and the thermodynamics
   !>
   !> @param[in,out] self      updated model
   !> @param[in,out] coupling  staged coupling with every mandatory answer
   !> @param[out]    error     VV or UV failure; a failed solve leaves no solution
   subroutine moz_1d_solve_uv(self, coupling, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Error
      type(error_type), allocatable, intent(out) :: error
      type(solvent_vv_thermodynamics_type) :: bulk
      type(model_moz_1d_thermodynamics_type) :: thermo
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :), uk_lr(:, :)
      character(len=96) :: buf
      real(wp) :: rmse, free_energy, moment
      integer :: nv, niter, nthreads, verbosity
      logical :: converged, reuse_guess

      call self%prepare_solvent(error)
      if (allocated(error)) return
      if (self%bulk_generation /= self%solvent%solution_generation()) then
         self%solved = .false.
         if (allocated(self%problem)) deallocate (self%problem)
         if (allocated(self%lambda)) deallocate (self%lambda)
      end if
      call self%potential%compute(self%rgrid, self%kgrid, u_sr, ur_lr, uk_lr, error, &
         & solvent=self%solvent%potential, coupling=coupling)
      if (allocated(error)) return
      if (allocated(self%problem)) then
         if (same_grid(self%rgrid, self%problem%rgrid) .and. same_grid(self%kgrid, self%problem%kgrid)) then
            if (self%solved) then
               if (same_table(u_sr, self%u_sr) .and. same_table(ur_lr, self%ur_lr) .and. &
                  & same_table(uk_lr, self%uk_lr)) return
            end if
         else
            self%solved = .false.
            deallocate (self%problem)
            if (allocated(self%lambda)) deallocate (self%lambda)
         end if
      end if
      reuse_guess = self%solved .and. self%params%reuse_initial_guess
      self%solved = .false.

      call self%operators%prepare(self%kgrid, self%solvent, self%xyz, error)
      if (allocated(error)) return
      nv = self%solvent%nsite()

      nthreads = self%ctx%get_num_threads()
      verbosity = self%ctx%verbosity()
      if (allocated(self%problem)) then
         if (.not. self%problem%equation%ready()) deallocate (self%problem)
      end if
      if (allocated(self%problem)) then
         call self%problem%update_potential(u_sr, ur_lr, uk_lr, error)
      else
         allocate (self%problem)
         call new_uv_problem(self%problem, self%rgrid, self%kgrid, self%operators%wu, self%operators%response, &
            & u_sr, ur_lr, uk_lr, &
            & spread(self%solvent%density(), 1, nv), self%solvent%beta(), self%params, nthreads, verbosity, error)
      end if
      if (allocated(error)) return
      if (verbosity >= 2) then
         call self%problem%solve(self%params%tolerance, self%params%maxiter, converged, rmse, niter, error, &
            & unit=self%ctx%unit, reuse_initial_guess=reuse_guess)
      else
         call self%problem%solve(self%params%tolerance, self%params%maxiter, converged, rmse, niter, error, &
            & reuse_initial_guess=reuse_guess)
      end if
      if (allocated(error)) return
      write (buf, "(i0,' iterations, RMS residual ',es9.2)") niter, rmse
      if (.not. converged) then
         call fatal_error(error, "1D UV solve did not converge: "//trim(buf))
         return
      end if
      call self%problem%free_energy(free_energy, error)
      if (allocated(error)) return
      call self%problem%direct_moment(moment, error)
      if (allocated(error)) return
      call self%solvent%thermodynamics(bulk, error)
      if (allocated(error)) return
      call new_moz_1d_thermodynamics(thermo, free_energy, moment, bulk%compressibility, self%problem%beta, error)
      if (allocated(error)) return

      self%cr = self%problem%cr - self%problem%beta*ur_lr
      self%tr = self%problem%tr + self%problem%beta*ur_lr
      self%hr = self%problem%cr + self%problem%tr
      call move_alloc(u_sr, self%u_sr)
      call move_alloc(ur_lr, self%ur_lr)
      call move_alloc(uk_lr, self%uk_lr)
      self%thermo = thermo
      self%solved = .true.
      self%bulk_generation = self%solvent%solution_generation()
      if (verbosity >= 1) then
         call self%ctx%message("1D UV converged: "//trim(buf), 1)
         write (buf, "(f0.4)") free_energy*autokcal
         call self%ctx%message("  solvation free energy       "//trim(buf)//" kcal/mol", 1)
         write (buf, "(f0.4)") thermo%partial_molar_volume*autoaa**3
         call self%ctx%message("  partial molar volume        "//trim(buf)//" A^3", 1)
      end if
   end subroutine moz_1d_solve_uv

   !> Whether two tables are equal in shape and value
   !>
   !> @param[in] a  new table
   !> @param[in] b  kept table; unallocated compares unequal
   pure logical function same_table(a, b)
      !> Tables
      real(wp), intent(in) :: a(:, :)
      real(wp), allocatable, intent(in) :: b(:, :)
      same_table = .false.
      if (.not. allocated(b)) return
      if (any(shape(a) /= shape(b))) return
      same_table = all(a == b)
   end function same_table

   !> Copies of the physical correlations of the latest UV solve, (npts, nu nv)
   !>
   !> UV columns (i - 1) nv + a hold solute site i and solvent site a
   !>
   !> @param[in]  self   model after an energy at its current structure
   !> @param[out] cr     direct correlation c(r)
   !> @param[out] tr     indirect correlation t(r) = h(r) - c(r)
   !> @param[out] hr     total correlation h(r)
   !> @param[out] error  no UV solution at the current structure
   subroutine moz_1d_correlations(self, cr, tr, hr, error)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Physical correlations
      real(wp), allocatable, intent(out) :: cr(:, :), tr(:, :), hr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%solved) then
         call fatal_error(error, "1D MOZ model has no UV solution at its current structure")
         return
      end if
      cr = self%cr
      tr = self%tr
      hr = self%hr
   end subroutine moz_1d_correlations

   !> Free energy and partial molar volume of the latest UV solve
   !>
   !> @param[in]  self    model after an energy at its current structure
   !> @param[out] thermo  thermodynamics of the reported energy
   !> @param[out] error   no UV solution at the current structure
   subroutine moz_1d_thermodynamics(self, thermo, error)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Thermodynamics
      type(model_moz_1d_thermodynamics_type), intent(out) :: thermo
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%solved) then
         call fatal_error(error, "1D MOZ model has no UV solution at its current structure")
         return
      end if
      thermo = self%thermo
   end subroutine moz_1d_thermodynamics

   !> Derivative of the reported energy with respect to the host data of the solute terms
   !>
   !> - UV solution for the tables of the coupling (solve_uv), reused after an
   !>   energy with the same host answers
   !> - Table adjoints from one transposed solve at the `tolerance` setting, relative
   !>   to its right-hand side
   !> - The potential hands them to its terms: monopole returns atomic_charge_adjoint;
   !>   terms without host data return nothing
   !> - `response` is cleared once the adjoints are known; a refusal leaves it untouched
   !>
   !> @param[in,out] self      updated model
   !> @param[in,out] coupling  coupling staged for the response, with every mandatory answer
   !> @param[in,out] response  host response items
   !> @param[out]    error     foreign or incomplete coupling, model not updated, VV, UV or adjoint failure
   subroutine moz_1d_response(self, coupling, response, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Response
      type(response_type), intent(inout) :: response
      !> Error
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), w_k(:, :)

      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_response, error)
      if (allocated(error)) return
      call self%solve_uv(coupling, error)
      if (allocated(error)) return
      call self%problem%response(self%params%tolerance, w_sr, w_lr, w_k, self%lambda, error)
      if (allocated(error)) return
      call response_clear(response)
      call self%potential%compute_adjoint(self%rgrid, self%kgrid, self%solvent%potential, coupling, &
         & w_sr, w_lr, w_k, response, error)
   end subroutine moz_1d_response

   !> Nuclear gradient of the reported UV energy at fixed bulk and radial grids
   !>
   !> Site potentials are centered on each radial origin; at fixed parameters their
   !> explicit geometry derivative is zero. The solute W(k) carries pair-distance
   !> dependence. Host charge derivatives are returned for contraction by the host.
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] response Response
   !> @param[in,out] gradient Gradient
   !> @param[out] error Error
   subroutine moz_1d_gradient(self, coupling, response, gradient, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Response
      type(response_type), intent(inout) :: response
      !> Gradient
      real(wp), intent(inout) :: gradient(:, :)
      !> Error
      type(error_type), allocatable, intent(out) :: error
      type(packed_index_lut) :: lut
      real(wp), allocatable :: local(:, :), w_sr(:, :), w_lr(:, :), w_k(:, :), w_wu(:, :)
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_gradient, error)
      if (allocated(error)) return
      if (any(shape(gradient) /= [3, self%atom_count()])) then
         call fatal_error(error, "1D MOZ gradient shape mismatch")
         return
      end if
      call self%solve_uv(coupling, error)
      if (allocated(error)) return
      call self%problem%response(self%params%tolerance, w_sr, w_lr, w_k, self%lambda, error, w_wu=w_wu)
      if (allocated(error)) return
      allocate (local(3, self%atom_count()), source=0.0_wp)
      call lut%init(self%atom_count())
      call vv_intramolecular_adjoint(self%kgrid%r, self%xyz, lut, w_wu, local)
      call response_clear(response)
      call self%potential%compute_adjoint(self%rgrid, self%kgrid, self%solvent%potential, coupling, &
         & w_sr, w_lr, w_k, response, error)
      if (allocated(error)) return
      gradient = gradient + local
   end subroutine moz_1d_gradient

   !> Release the coupling registry
   !>
   !> @param[in,out] self Model being destroyed
   subroutine destroy_moz_1d_model(self)
      !> Model being destroyed
      type(model_moz_1d_type), intent(inout) :: self
      call self%clear_couplings()
   end subroutine destroy_moz_1d_model

end module moist_model_moz_1d_type
