!> Three-dimensional MOZ model state
!>
!> - Owned spatial grid copy; the theory uses only its points, weights, k-points
!>   and transform engine, so any point grid serves
!> - Solute potential `potential` (public) and a private copy of the updated bulk solvent
!> - Settings of the UV solve and, nested as `solvent`, of the VV solve of the
!>   solvent copy; `prepare_solvent` solves it unless its state and residual match
!> - Energy: 3D UV solve on the grid with the solvent response at its k-points; the
!>   tables come from `potential%compute(grid, ..., solvent=..., coupling=...)`
!> - Partial molar volume and PC/PC+ terms of every energy; the `correction`
!>   setting selects the term added to the energy
!> - Optional leading net-charge finite-size quantity; fixed charge, bulk and lattice
!> - Tables equal to those of the kept solution at this structure reuse it;
!>   changed host answers keep solver storage/plans and can reuse the converged fields
!> - Fixed-grid, fixed-bulk response through one shared RISM adjoint solve
!> - Nuclear gradient through potential and transform pullbacks and one grid contraction
!> - Coupling scope of solute term i: i
!> - Point grids only, by design: Gaussian widths give the tables nothing and cost more
!> - Optional checkpoint store: initial guess of a solve without a reusable
!>   field, saves during and after each UV solve
!> - The model does not watch its potential: update it again after a change
module moist_model_moz_3d_type
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_convert, only: autokcal, autoaa
   use mctc_io_utils, only: to_lower
   use moist_context, only: moist_context_type
   use moist_output_citations, only: citation
   use moist_utils_timer, only: cat_setup, cat_solve, cat_energy, cat_properties, cat_gradient
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_model_type, only: solvation_model_type, invalidate_model
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_channels_fields, only: field_query_type
   use moist_channels_coupling, only: coupling_type, coupling_begin_registration, moist_phase_energy, &
      & moist_phase_response, moist_phase_gradient, coupling_snapshot, coupling_check_mandatory
   use moist_channels_response, only: response_type, response_clear
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_solvent_vv, only: solvent_vv_type, solvent_vv_thermodynamics_type
   use moist_model_moz_thermodynamics, only: model_moz_3d_thermodynamics_type, new_moz_3d_thermodynamics, &
      & thermo_list_type, accumulate_solvation_adjoint, huenenberger_correction_type, new_huenenberger_correction
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type, new_uv3d_susceptibility
   use moist_model_moz_3d_driver, only: uv3d_problem_type, new_uv3d_problem
   use moist_model_moz_checkpoint, only: moz_checkpoint_type, checkpoint_kind_uv3d
   implicit none(type, external)
   private
   public :: model_moz_3d_type, new_moz_3d_model

   !> Construct from a grid and an updated solvent; the input overload is in moist_model_moz_3d_input
   interface new_moz_3d_model
      module procedure new_moz_3d_model_grid
   end interface new_moz_3d_model

   !> Refusal of a grid with Gaussian widths
   character(len=*), parameter :: gaussian_grid_refused = &
      & "MOZ 3D model supports point grids only; Gaussian-width grids are not supported yet"

   !> Direct 3D MOZ model; UV energy, host response and nuclear gradient
   type, extends(solvation_model_type) :: model_moz_3d_type
      !> Solute potential; add terms before `update`
      type(moz_potential_type) :: potential
      !> Owned spatial grid
      class(moist_math_grid_3d_type), allocatable :: grid
      !> Copy of the updated bulk solvent
      type(solvent_vv_type), allocatable, private :: solvent
      !> Settings of the UV solve and of the VV solve of the solvent copy
      type(model_moz_3d_parameters_type), private :: params
      !> Solvent response on the k-points of the grid of the latest update
      type(uv3d_susceptibility_type), private :: response
      !> Whether the response belongs to the grid of the latest update
      logical, private :: response_ready = .false.
      !> Whether the correlations belong to the latest energy at this structure
      logical, private :: solved = .false.
      !> Physical c, t = h - c, h of last UV solve, (ngrid, nv)
      real(wp), allocatable, private :: cr(:, :), tr(:, :), hr(:, :)
      !> Free energy, partial molar volume and correction terms of the latest UV solve
      type(model_moz_3d_thermodynamics_type), private :: thermo
      !> Cached thermodynamic dependency graph of the same UV solution
      type(thermo_list_type), private :: thermo_items
      !> Solute of the latest update; its total charge is independent of host density
      type(structure_type), private :: mol
      !> Cartesian box lengths captured on successful grid update, bohr; zero otherwise
      real(wp), private :: box_lengths(3) = 0.0_wp
      !> Finite-size quantity cached for the current structure, lattice and bulk solution
      type(huenenberger_correction_type), allocatable, private :: finite_size
      !> Latest UV problem, kept for reuse
      type(uv3d_problem_type), allocatable, private :: problem
      !> Potential tables of the kept solution
      real(wp), allocatable, private :: u_sr(:, :), ur_lr(:, :)
      !> Long-range reciprocal potential of the kept solution, (npts_k, nv)
      complex(wp), allocatable, private :: uk_lr(:, :)
      !> Last transposed solution, reused as a starting point on the same grid
      real(wp), allocatable, private :: lambda(:, :)
      !> Optional store of UV states
      class(moz_checkpoint_type), allocatable, private :: checkpoint
   contains
      final :: destroy_moz_3d_model
      procedure :: invalidate => moz_3d_invalidate
      !> Attach a UV checkpoint store, or detach it
      procedure :: set_checkpoint => moz_3d_set_checkpoint
      procedure :: update => moz_3d_update
      procedure :: get_energy => moz_3d_energy
      procedure :: get_response => moz_3d_response
      procedure :: get_gradient => moz_3d_gradient
      procedure :: atom_count => moz_3d_atom_count
      !> Declare the coupling requests of every solute term in its own scope
      procedure :: declare_pass => moz_3d_declare_pass
      !> Publish the grid geometry as named fields
      procedure :: list_fields => moz_3d_list_fields
      !> Whether the model holds a solvent
      procedure :: has_solvent => moz_3d_has_solvent
      !> Number of solvent sites; 0 before construction
      procedure :: solvent_nsite => moz_3d_solvent_nsite
      !> Solve the solvent copy unless its state and residual satisfy the requested VV settings
      procedure :: prepare_solvent => moz_3d_prepare_solvent
      !> Whether the solvent copy holds a VV solution satisfying the requested settings
      procedure :: solvent_solved => moz_3d_solvent_solved
      !> Copies of the physical correlations of the latest UV solve
      procedure :: correlations => moz_3d_correlations
      !> Free energy, partial molar volume and correction terms of the latest UV solve
      procedure :: thermodynamics => moz_3d_thermodynamics
      !> UV solution for the tables of the coupling, reused when they are unchanged
      procedure, private :: solve_uv => moz_3d_solve_uv
      procedure, private :: prepare_finite_size => moz_3d_prepare_finite_size
   end type model_moz_3d_type

contains

   !> Copy the configured spatial grid, an updated bulk solvent and the settings
   !>
   !> - Gaussian widths already present: refusal at construction
   !> - Geometry-dependent Gaussian widths: refusal at `update`
   !> - `prepare_solvent` reuses a solved copy only when the bulk state matches
   !>   and its achieved residual meets the requested tolerance; otherwise it solves again
   !> - Without explicit parameters, a solved solvent supplies its accepted VV settings
   !>
   !> @param[out] self        model
   !> @param[in]  ctx         run context, borrowed for the lifetime of the model
   !> @param[in]  grid        spatial point grid template
   !> @param[in]  solvent     bulk solvent with an updated potential, copied into the model
   !> @param[out] error       gaussian-width grid, solvent not updated, invalid settings or allocation failure
   !> @param[in]  parameters  UV/VV settings; absent parameters retain accepted bulk settings or use defaults
   subroutine new_moz_3d_model_grid(self, ctx, grid, solvent, error, parameters)
      !> Model
      type(model_moz_3d_type), intent(out) :: self
      !> Run context, borrowed for the lifetime of the model
      type(moist_context_type), intent(in), target :: ctx
      !> Spatial point grid template, copied into the model
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Bulk solvent with an updated potential
      type(solvent_vv_type), intent(in) :: solvent
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Settings of the UV and VV solves
      type(model_moz_3d_parameters_type), intent(in), optional :: parameters
      !> Allocation status of the copies
      integer :: stat
      if (allocated(grid%xi0)) then
         call fatal_error(error, gaussian_grid_refused)
         return
      end if
      if (.not. solvent%potential%is_updated()) then
         call fatal_error(error, "Update the 1D VV solvent before constructing the 3D MOZ model")
         return
      end if
      if (present(parameters)) then
         self%params = parameters
      else if (solvent%is_solved()) then
         self%params%solvent = solvent%settings()
      end if
      call self%params%validate(error)
      if (allocated(error)) return
      ! TODO: share the solvent check and copy with the 1D model
      allocate (self%solvent, source=solvent, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to copy the 1D VV solvent")
         return
      end if
      self%ctx => ctx
      allocate (self%grid, source=grid, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to copy 3D MOZ grid")
         return
      end if
      self%grid%nthreads = ctx%get_num_threads()
   end subroutine new_moz_3d_model_grid

   !> Drop UV results and the response cache when the model is invalidated
   !>
   !> @param[inout] self  model
   subroutine moz_3d_invalidate(self)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      call invalidate_model(self)
      self%solved = .false.
      self%response_ready = .false.
      if (allocated(self%lambda)) deallocate (self%lambda)
      if (allocated(self%problem)) deallocate (self%problem)
      if (allocated(self%finite_size)) deallocate (self%finite_size)
      self%box_lengths = 0.0_wp
   end subroutine moz_3d_invalidate

   !> Update the owned grid geometry, then the solute potential on it
   !>
   !> - Potential table printed at verbosity 2
   !> - Grid-dependent term data (EC fit) built from the updated grid
   !> - Solute charges must sum to the structure's total charge; host-data charges are fitted later
   !>
   !> @param[in,out] self   model
   !> @param[in]     mol    solute structure
   !> @param[out]    error  grid error, Gaussian-width grid, empty structure, no terms, term failure
   !>                       or charges not summing to the total charge
   subroutine moz_3d_update(self, mol, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Solute structure
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%invalidate()
      if (.not. allocated(self%grid) .or. .not. allocated(self%solvent)) then
         call fatal_error(error, "Construct the 3D MOZ model before update")
         return
      end if
      call self%grid%update(mol, error)
      if (allocated(error)) return
      if (allocated(self%grid%xi0)) then
         call fatal_error(error, gaussian_grid_refused)
         return
      end if
      select type (grid => self%grid)
      class is (moist_math_grid_3d_cartesian_type)
         self%box_lengths = grid%dr*real([grid%nx, grid%ny, grid%nz], wp)
      end select
      if (to_lower(trim(adjustl(self%params%finite_size))) /= "none") then
         if (any(self%box_lengths <= 0.0_wp)) then
            call fatal_error(error, "huenenberger correction requires a Cartesian grid")
            return
         end if
      end if
      call self%potential%update(mol, self%grid, error)
      if (.not. allocated(error)) call self%potential%check_total_charge(error)
      if (allocated(error)) return
      self%mol = mol
      call self%ctx%cite(citation%rism3d)
      if (self%potential%ng_split()) call self%ctx%cite(citation%ng)
      if (to_lower(trim(adjustl(self%params%correction))) /= "none") call self%ctx%cite(citation%pressure_correction)
      if (to_lower(trim(adjustl(self%params%finite_size))) /= "none") call self%ctx%cite(citation%finite_size)
      if (self%ctx%verbosity(2)) then
         call self%ctx%message("Solute potential parameters", 2)
         call self%potential%print_table(self%ctx%unit, error)
         if (allocated(error)) return
         call self%ctx%message("", 2)
      end if
      self%updated = .true.
      self%built = .true.
   end subroutine moz_3d_update

   !> Site count of the latest valid structure, from the solute potential
   !>
   !> @param[in] self  model
   function moz_3d_atom_count(self) result(natom)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Site count
      integer :: natom
      natom = self%potential%natom()
   end function moz_3d_atom_count

   !> Declare the coupling requests of every solute term in its own scope
   !>
   !> - Coupling scope of term i: i
   !> - No requests from parameter terms
   !> - EC term declares point-potential probes at grid points
   !>
   !> @param[in,out] self      updated model
   !> @param[in,out] coupling  coupling to declare
   !> @param[out]    error     failed declaration
   subroutine moz_3d_declare_pass(self, coupling, error)
      !> Updated model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling to declare
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call coupling_begin_registration(coupling)
      call self%potential%declare(self%grid, coupling, error)
      if (allocated(error)) return
      call coupling_snapshot(coupling, ngrid=self%grid%ngrid, natom=self%grid%natom)
   end subroutine moz_3d_declare_pass

   !> Whether the model holds a solvent
   !>
   !> @param[in] self  model
   pure function moz_3d_has_solvent(self) result(has)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Presence
      logical :: has
      has = allocated(self%solvent)
   end function moz_3d_has_solvent

   !> Number of solvent sites; 0 before construction
   !>
   !> @param[in] self  model
   pure function moz_3d_solvent_nsite(self) result(nsite)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Site count
      integer :: nsite
      nsite = 0
      if (allocated(self%solvent)) nsite = self%solvent%nsite()
   end function moz_3d_solvent_nsite

   !> Publish the grid geometry as named fields, the model's evaluation domain
   !>
   !> @param[in]     self   model
   !> @param[in,out] query  field walker
   subroutine moz_3d_list_fields(self, query)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query
      if (.not. allocated(self%grid)) return
      call query%add_int_value("ngrid", "Number of volume grid points", self%grid%ngrid)
      call query%add_int_value("natom", "Number of solute atoms", self%grid%natom)
      call query%add_real2("xyz", "Grid coordinates, bohr (3, ngrid)", self%grid%xyz)
      call query%add_real("w", "Integration weights, bohr**3 (ngrid)", self%grid%w)
      call query%add_int("owner", "Atom owner, zero based; -1 for unowned points", &
         & self%grid%owner, zero_based=.true.)
   end subroutine moz_3d_list_fields

   !> Add the solvation free energy of the 3D UV solve
   !>
   !> - UV solution for the tables of the coupling (solve_uv)
   !> - Adds (1/beta) sum_a rho_a int Phi(h_a, c_a) dV in Hartree plus the term the
   !>   `correction` setting selects
   !>
   !> @param[in,out] self      updated model
   !> @param[in,out] coupling  coupling staged for the energy, with every mandatory answer
   !> @param[in,out] energy    energy, Hartree; unchanged on error
   !> @param[out]    error     foreign or incomplete coupling, model not updated, VV or UV failure
   subroutine moz_3d_energy(self, coupling, energy, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling owned by the model
      class(coupling_type), intent(inout), target :: coupling
      !> Energy
      real(wp), intent(inout) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: evaluated_energy
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_energy, error)
      if (allocated(error)) return
      call self%solve_uv(coupling, error)
      if (allocated(error)) return
      call self%thermo_items%energy(evaluated_energy, error)
      if (allocated(error)) return
      energy = energy + evaluated_energy
   end subroutine moz_3d_energy

   !> Solve the VV equations of the solvent copy unless it holds a current solution
   !>
   !> - Uses the `solvent` settings and their grid; the solvent of the caller is not touched
   !> - Called by the energy; callers may call it early to fail fast
   !>
   !> @param[in,out] self   constructed model
   !> @param[out]    error  unconstructed model or VV solve failure
   subroutine moz_3d_prepare_solvent(self, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(self%solvent)) then
         call fatal_error(error, "Construct the 3D MOZ model before preparing its solvent")
         return
      end if
      if (self%solvent%is_solved(self%params%solvent)) return
      self%response_ready = .false.
      call self%solvent%solve(error, self%params%solvent)
   end subroutine moz_3d_prepare_solvent

   !> Whether the solvent copy holds a VV solution satisfying the requested settings
   !>
   !> @param[in] self  model
   pure function moz_3d_solvent_solved(self) result(solved)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Solution flag
      logical :: solved
      solved = .false.
      if (allocated(self%solvent)) solved = self%solvent%is_solved(self%params%solvent)
   end function moz_3d_solvent_solved

   !> Set up the constant finite-size correction once
   !>
   !> - Host charge partitions must conserve Q
   !>
   !> @param[inout] self    model
   !> @param[in]    charge  net charge Q of the potential
   !> @param[out]   error   missing analytical Ng Coulomb split or correction failure
   subroutine moz_3d_prepare_finite_size(self, charge, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Total solute charge Q
      real(wp), intent(in) :: charge
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(huenenberger_correction_type) :: quantity
      real(wp) :: screening

      if (to_lower(trim(adjustl(self%params%finite_size))) == "none") return
      if (.not. self%potential%ng_split()) then
         call fatal_error(error, "huenenberger correction requires the analytical Ng Coulomb split")
         return
      end if
      if (abs(charge - self%mol%charge) > 1.0e-6_wp) then
         call fatal_error(error, "huenenberger correction requires solute tail charges summing to mol%charge")
         return
      end if
      if (allocated(self%finite_size)) return
      call self%solvent%charge_screening(screening, error)
      if (allocated(error)) return
      call new_huenenberger_correction(quantity, self%mol%charge, screening, self%box_lengths, error)
      if (allocated(error)) return
      allocate (self%finite_size, source=quantity)
   end subroutine moz_3d_prepare_finite_size

   !> UV solution for the potential tables of the coupling
   !>
   !> - Solves the solvent copy first when needed (prepare_solvent)
   !> - Tables from the solute potential with the host answers on the coupling;
   !>   tables equal to those of the kept solution at this structure reuse it
   !> - Solvent response on the grid's k-points, once per update (`susceptibility` setting)
   !> - Changed host tables retain the equation and solver resources on the same grid;
   !>   reuse_initial_guess starts from the previous accepted fields with a cold fallback
   !> - Timed on the context timer as `3D UV RISM`, after the solvent solve, with
   !>   its setup, iteration, energy and property phases
   !>
   !> @param[in,out] self      updated model
   !> @param[in,out] coupling  staged coupling with every mandatory answer
   !> @param[out]    error     VV, response or UV failure; a failed solve leaves no solution
   subroutine moz_3d_solve_uv(self, coupling, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Timer depth on entry, restored on every return
      integer :: d0

      call self%prepare_solvent(error)
      if (allocated(error)) return
      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("3D UV RISM", category=cat_solve)
      call solve_uv_tables(self, coupling, error)
      call self%ctx%timer%unwind(d0)
   end subroutine moz_3d_solve_uv

   !> Body of solve_uv after the solvent solve; an early return leaves its phase
   !> timer for the caller to unwind
   !>
   !> @param[in,out] self      updated model with a solved solvent
   !> @param[in,out] coupling  staged coupling with every mandatory answer
   !> @param[out]    error     response or UV failure; a failed solve leaves no solution
   subroutine solve_uv_tables(self, coupling, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(solvent_vv_thermodynamics_type) :: bulk
      type(model_moz_3d_thermodynamics_type) :: thermo
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :), guess(:, :)
      complex(wp), allocatable :: uk_lr(:, :)
      type(prettyprinter) :: pp
      character(len=96) :: buf
      character(len=:), allocatable :: correction_label
      real(wp) :: rmse, free_energy, moment, reported_energy, tail_charge
      integer :: nv, niter
      logical :: converged, reuse_guess

      call self%ctx%timer%start("Setup", category=cat_setup)
      if (self%response_ready) then
         if (.not. self%response%matches(self%grid, self%solvent, self%params%susceptibility)) then
            self%response_ready = .false.
            self%solved = .false.
            if (allocated(self%problem)) deallocate (self%problem)
            if (allocated(self%lambda)) deallocate (self%lambda)
            if (allocated(self%finite_size)) deallocate (self%finite_size)
         end if
      end if
      call self%potential%compute(self%grid, u_sr, ur_lr, uk_lr, error, &
         & solvent=self%solvent%potential, coupling=coupling, net_charge=tail_charge)
      if (allocated(error)) return
      call self%prepare_finite_size(tail_charge, error)
      if (allocated(error)) then
         self%solved = .false.
         return
      end if
      if (self%solved .and. allocated(self%problem)) then
         if (same_real(u_sr, self%u_sr) .and. same_real(ur_lr, self%ur_lr) .and. same_complex(uk_lr, self%uk_lr)) then
            return
         end if
      end if
      reuse_guess = self%solved .and. self%params%reuse_initial_guess
      self%solved = .false.
      if (.not. self%response_ready) then
         call new_uv3d_susceptibility(self%response, self%grid, self%solvent, self%params%susceptibility, error)
         if (allocated(error)) return
         self%response_ready = .true.
      end if

      nv = self%solvent%nsite()
      if (allocated(self%problem)) then
         if (.not. self%problem%equation%ready()) deallocate (self%problem)
      end if
      if (allocated(self%problem)) then
         call self%problem%update_potential(u_sr, ur_lr, uk_lr, error)
      else
         allocate (self%problem)
         call new_uv3d_problem(self%problem, self%ctx, self%grid, self%response, u_sr, ur_lr, uk_lr, &
            & spread(self%solvent%density(), 1, nv), self%solvent%beta(), self%params, error)
      end if
      if (allocated(error)) return
      if (allocated(self%checkpoint)) then
         call set_checkpoint_header(self)
         if (.not. reuse_guess) call self%checkpoint%load_guess(self%ctx, self%problem%beta, ur_lr, 4, guess)
      end if
      call self%ctx%timer%stop("Setup")
      call self%problem%solve(self%params%tolerance, self%params%maxiter, converged, rmse, niter, error, &
         & reuse_initial_guess=reuse_guess, guess=guess, checkpoint=self%checkpoint)
      if (allocated(error)) return
      write (buf, "(i0,' iterations, RMS residual ',es9.2)") niter, rmse
      if (.not. converged) then
         call fatal_error(error, "3D UV solve did not converge: "//trim(buf))
         return
      end if
      call self%ctx%timer%start("Energy", category=cat_energy)
      call self%problem%free_energy(free_energy, error)
      if (allocated(error)) return
      call self%ctx%timer%stop("Energy")
      call self%ctx%timer%start("Properties", category=cat_properties)
      call self%problem%direct_moment(moment, error)
      if (allocated(error)) return
      call self%solvent%thermodynamics(bulk, error)
      if (allocated(error)) return
      ! Branch instead of passing an unallocated actual: ifx treats it as present for a class dummy
      if (allocated(self%finite_size)) then
         call new_moz_3d_thermodynamics(thermo, free_energy, moment, bulk, self%problem%beta, error, &
            & items=self%thermo_items, correction=self%params%correction, finite_size=self%finite_size)
      else
         call new_moz_3d_thermodynamics(thermo, free_energy, moment, bulk, self%problem%beta, error, &
            & items=self%thermo_items, correction=self%params%correction)
      end if
      if (allocated(error)) return
      call self%thermo_items%energy(reported_energy, error)
      if (allocated(error)) return
      call self%ctx%timer%stop("Properties")

      self%cr = self%problem%cr - self%problem%beta*ur_lr
      self%tr = self%problem%tr + self%problem%beta*ur_lr
      self%hr = self%problem%cr + self%problem%tr
      call move_alloc(u_sr, self%u_sr)
      call move_alloc(ur_lr, self%ur_lr)
      call move_alloc(uk_lr, self%uk_lr)
      self%thermo = thermo
      self%solved = .true.
      if (self%ctx%verbosity(3)) then
         pp = new_prettyprinter(unit=self%ctx%unit)
         call pp%push("3D UV solution:")
         call pp%kv("Iterations", niter)
         call pp%kv("RMS residual", rmse, use_exp=.true.)
         call pp%kv("Solvation free energy", free_energy*autokcal, "kcal/mol")
         call pp%kv("Partial molar volume", thermo%partial_molar_volume*autoaa**3, "A^3")
         call pp%kv("PC correction", thermo%pc*autokcal, "kcal/mol")
         call pp%kv("PC+ correction", thermo%pc_plus*autokcal, "kcal/mol")
         if (allocated(self%finite_size)) then
            call pp%kv("huenenberger correction", thermo%finite_size*autokcal, "kcal/mol")
         end if
         correction_label = trim(self%params%correction)
         if (allocated(self%finite_size)) then
            if (to_lower(correction_label) == "none") then
               correction_label = trim(self%params%finite_size)
            else
               correction_label = correction_label//", "//trim(self%params%finite_size)
            end if
         end if
         call pp%kv("Reported free energy", reported_energy*autokcal, "kcal/mol ("//correction_label//")")
         call pp%pop()
         call pp%blank()
      end if
   end subroutine solve_uv_tables

   !> Copies of the physical correlations of the latest UV solve, (ngrid, nv)
   !>
   !> Column a holds solvent site a on the model's grid
   !>
   !> @param[in]  self   model after an energy at its current structure
   !> @param[out] cr     direct correlation c(r)
   !> @param[out] tr     indirect correlation t(r) = h(r) - c(r)
   !> @param[out] hr     total correlation h(r)
   !> @param[out] error  no UV solution at the current structure
   subroutine moz_3d_correlations(self, cr, tr, hr, error)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Physical correlations
      real(wp), allocatable, intent(out) :: cr(:, :), tr(:, :), hr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%solved) then
         call fatal_error(error, "3D MOZ model has no UV solution at its current structure")
         return
      end if
      cr = self%cr
      tr = self%tr
      hr = self%hr
   end subroutine moz_3d_correlations

   !> Free energy, partial molar volume and correction terms of the latest UV solve
   !>
   !> @param[in]  self    model after an energy at its current structure
   !> @param[out] thermo  thermodynamics; the free energy is uncorrected
   !> @param[out] error   no UV solution at the current structure
   subroutine moz_3d_thermodynamics(self, thermo, error)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Thermodynamics
      type(model_moz_3d_thermodynamics_type), intent(out) :: thermo
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%solved) then
         call fatal_error(error, "3D MOZ model has no UV solution at its current structure")
         return
      end if
      thermo = self%thermo
   end subroutine moz_3d_thermodynamics

   !> Attach a copy of a UV checkpoint store, or detach the store when absent
   !>
   !> - The store supplies the initial guess of a UV solve without a reusable
   !>   field and receives its states
   !> - The VV store travels with the solvent passed to the constructor
   !>
   !> @param[inout] self        model
   !> @param[in]    checkpoint  store to copy
   subroutine moz_3d_set_checkpoint(self, checkpoint)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Store to copy
      class(moz_checkpoint_type), intent(in), optional :: checkpoint
      if (allocated(self%checkpoint)) deallocate (self%checkpoint)
      if (present(checkpoint)) allocate (self%checkpoint, source=checkpoint)
   end subroutine moz_3d_set_checkpoint

   !> Describe the states of the current UV solve in the store header
   !>
   !> - Cartesian grids by their dimensions, spacing and origin; other point grids
   !>   by their point count
   !> - Solvent and solute sites, solute coordinates, closure and temperature
   !>
   !> @param[inout] self  updated model with a store
   subroutine set_checkpoint_header(self)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      type(structure_type) :: solvent

      solvent = self%solvent%structure()
      associate (header => self%checkpoint%header)
         header%kind = checkpoint_kind_uv3d
         if (allocated(header%origin)) deallocate (header%origin)
         select type (grid => self%grid)
         class is (moist_math_grid_3d_cartesian_type)
            header%grid_shape = [grid%nx, grid%ny, grid%nz]
            header%spacing = grid%dr
            header%origin = grid%get_origin()
         class default
            header%grid_shape = [self%grid%ngrid]
            header%spacing = 0.0_wp
         end select
         header%solvent_sites = solvent%num(solvent%id)
         header%solute_sites = self%mol%num(self%mol%id)
         header%solute_xyz = self%mol%xyz
         header%closure = trim(adjustl(self%params%closure))
         header%theory = trim(self%solvent%theory_name())
         header%temperature = self%solvent%temperature()
      end associate
   end subroutine set_checkpoint_header

   !> Whether two real tables are equal in shape and value
   !>
   !> @param[in] a  new table
   !> @param[in] b  kept table; unallocated compares unequal
   pure logical function same_real(a, b)
      !> New table
      real(wp), intent(in) :: a(:, :)
      !> Kept table
      real(wp), allocatable, intent(in) :: b(:, :)
      same_real = .false.
      if (.not. allocated(b)) return
      if (any(shape(a) /= shape(b))) return
      same_real = all(a == b)
   end function same_real

   !> Whether two complex tables are equal in shape and value
   !>
   !> @param[in] a  new table
   !> @param[in] b  kept table; unallocated compares unequal
   pure logical function same_complex(a, b)
      !> New table
      complex(wp), intent(in) :: a(:, :)
      !> Kept table
      complex(wp), allocatable, intent(in) :: b(:, :)
      same_complex = .false.
      if (.not. allocated(b)) return
      if (any(shape(a) /= shape(b))) return
      same_complex = all(a == b)
   end function same_complex

   !> Host response of the selected energy at fixed geometry and bulk solvent
   !>
   !> @param[in,out] self      model
   !> @param[in,out] coupling  coupling owned by the model
   !> @param[in,out] response  host field and charge sensitivities
   !> @param[out]    error     coupling, UV solve, adjoint solve or potential pullback failure
   subroutine moz_3d_response(self, coupling, response, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling owned by the model
      class(coupling_type), intent(inout), target :: coupling
      !> Host response
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: energy_seed, shift
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :)
      complex(wp), allocatable :: w_k(:, :)
      !> Timer depth before the response phases
      integer :: d0
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_response, error)
      if (allocated(error)) return
      call self%solve_uv(coupling, error)
      if (allocated(error)) return
      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("3D UV response", category=cat_gradient)
      phases: block
         call self%ctx%timer%start("Adjoint solve")
         call accumulate_solvation_adjoint(self%thermo_items, energy_seed, shift, error)
         if (allocated(error)) exit phases
         call self%problem%response(shift, self%params%tolerance, w_sr, w_lr, w_k, self%lambda, error, &
            & free_energy_weight=energy_seed)
         if (allocated(error)) exit phases
         call self%ctx%timer%stop("Adjoint solve")
         call self%ctx%timer%start("Contraction")
         call response_clear(response)
         call self%potential%compute_adjoint(self%grid, self%solvent%potential, coupling, &
            & w_sr, w_lr, w_k, response, error)
      end block phases
      call self%ctx%timer%unwind(d0)
   end subroutine moz_3d_response

   !> Nuclear contribution at fixed bulk and reciprocal lattice, plus host response
   !>
   !> @param[in,out] self      model
   !> @param[in,out] coupling  coupling owned by the model
   !> @param[in,out] response  host sensitivities, cleared after primal/adjoint success
   !> @param[in,out] gradient  nuclear gradient accumulator, Hartree/bohr
   !> @param[out]    error     invalid request, primal/adjoint solve or geometry pullback failure
   subroutine moz_3d_gradient(self, coupling, response, gradient, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling owned by the model
      class(coupling_type), intent(inout), target :: coupling
      !> Host response
      type(response_type), intent(inout) :: response
      !> Nuclear gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: energy_seed, shift
      real(wp), allocatable :: local(:, :), w_sr(:, :), w_lr(:, :)
      complex(wp), allocatable :: w_k(:, :)
      type(volume_adjoint_type) :: acc
      !> Timer depth before the gradient phases
      integer :: d0
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_gradient, error)
      if (allocated(error)) return
      if (any(shape(gradient) /= [3, self%grid%natom])) then
         call fatal_error(error, "3D MOZ gradient shape mismatch")
         return
      end if
      call self%solve_uv(coupling, error)
      if (allocated(error)) return
      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("3D UV gradient", category=cat_gradient)
      phases: block
         call self%ctx%timer%start("Adjoint solve")
         call accumulate_solvation_adjoint(self%thermo_items, energy_seed, shift, error)
         if (allocated(error)) exit phases
         call acc%init(self%grid%ngrid)
         call self%problem%response(shift, self%params%tolerance, w_sr, w_lr, w_k, self%lambda, error, &
            & free_energy_weight=energy_seed, grid_adjoint=acc)
         if (allocated(error)) exit phases
         call self%ctx%timer%stop("Adjoint solve")
         call self%ctx%timer%start("Contraction")
         allocate (local(3, self%grid%natom), source=0.0_wp)
         call response_clear(response)
         call self%potential%compute_adjoint(self%grid, self%solvent%potential, coupling, &
            & w_sr, w_lr, w_k, response, error, gradient=local, grid_adjoint=acc)
         if (allocated(error)) exit phases
         call self%grid%get_volume_gradient(acc, local, error)
         if (allocated(error)) exit phases
         gradient = gradient + local
      end block phases
      call self%ctx%timer%unwind(d0)
   end subroutine moz_3d_gradient

   !> Release the coupling registry
   !>
   !> @param[in,out] self  model being destroyed
   subroutine destroy_moz_3d_model(self)
      !> Model being destroyed
      type(model_moz_3d_type), intent(inout) :: self
      call self%clear_couplings()
   end subroutine destroy_moz_3d_model

end module moist_model_moz_3d_type
