!> Bulk solvent of the MOZ models
!>
!> - Not a solvation model: holds the solvent system (solvent id, temperature,
!>   number density, structure) and the solvent potential
!> - Potential refuses host-data terms and carries the mixing rule and Ng split
!>   of the VV tables
!> - One site per atom
!> - Custom solvent: hand-filled `solvent_system_type` (solvent id 0,
!>   temperature, number density and structure; DRISM also needs epsilon >= 1)
!> - DRISM takes its dielectric constant from the solvent system
!> - Owns the 1D VV solve, the solvent susceptibility and the bulk thermodynamics
!> - The UV models copy it
!>
!>    call solvent%potential%add(lj_cspce, error)
!>    call solvent%potential%add(coulomb_cspce, error)
!>    call solvent%update(error)
!>    call solvent%potential%compute(rgrid, kgrid, u_sr, ur_lr, uk_lr, error)
!>    call solvent%solve(error, parameters)
!>
!> - The solvent does not watch its potential: update it again after a change
!> - A change also makes a VV solution stale
!> - An attached checkpoint store supplies the initial guess of each solve and
!>   receives its states; copies of the solvent carry the store
module moist_model_moz_solvent_vv
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_codata2018, only: Boltzmann_constant, atomic_unit_of_energy
   use mctc_io_convert, only: autokcal
   use moist_context, only: moist_context_type
   use moist_utils_timer, only: cat_setup, cat_solve, cat_properties
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_solvent_system, only: solvent_system_type
   use moist_output_citations, only: citation
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo
   use moist_math_packed_sym, only: packed_index_lut, npair_from_ns
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_potential_sites, only: potential_sites_type
   use moist_model_moz_potential_kernel_ng, only: ng_real_tail
   use moist_model_moz_solvent_vv_kernel, only: vv_intramolecular, vv_transform_at
   use moist_model_moz_solvent_vv_response, only: vv_response_type, new_vv_response
   use moist_model_moz_solvent_vv_parameters, only: solvent_vv_parameters_type
   use moist_model_moz_solvent_theory, only: bulk_theory_type, new_bulk_theory
   use moist_model_moz_solvent_vv_driver, only: vv_problem_type, new_vv_problem
   use moist_model_moz_solvent_vv_thermodynamics, only: solvent_vv_thermodynamics_type, new_vv_thermodynamics
   use moist_model_moz_checkpoint, only: moz_checkpoint_type, checkpoint_kind_vv
   implicit none(type, external)
   private

   public :: solvent_vv_type, new_vv_solvent, solvent_vv_parameters_type, solvent_vv_thermodynamics_type

   !> Unique identifiers for successful bulk solves, including solves with unchanged potentials
   integer(int64) :: last_solution_stamp = 0_int64

   !> Bulk solvent: system and potential
   type :: solvent_vv_type
      !> Run context, borrowed for the lifetime of the solvent
      type(moist_context_type), pointer :: ctx => null()
      !> Solvent system: id, temperature, number density and structure (bohr)
      type(solvent_system_type), private :: system
      !> Solvent potential; refuses host-data terms. Second side of the UV tables
      type(moz_potential_type) :: potential
      !> Whether the solvent was constructed
      logical, private :: constructed = .false.
      !> Whether the VV solution is available
      logical, private :: solved = .false.
      !> Potential change stamp the VV solution was computed for
      integer(int64), private :: solved_generation = -1_int64
      !> Identity of the accepted bulk solution
      integer(int64), private :: solution_stamp = -1_int64
      !> Achieved RMS residual of the accepted bulk solution
      real(wp), private :: achieved_residual = huge(1.0_wp)
      !> Radial pair of the VV solution
      type(moist_math_grid_radial_type), private :: rgrid, kgrid
      !> Physical c(r), t(r) = h - c, h(r) of the solution, packed (npts, npair)
      real(wp), allocatable, private :: cr(:, :), tr(:, :), hr(:, :)
      !> Solvent susceptibility chi(k) = w(k) + rho h(k), packed (npts, npair)
      real(wp), allocatable, private :: chi(:, :)
      !> Native short-range direct correlation transform and its cubic interpolant
      type(vv_response_type), private :: oz_response
      !> Fixed theory of the accepted bulk solution
      class(bulk_theory_type), allocatable, private :: theory
      !> Settings of the accepted bulk solution
      type(solvent_vv_parameters_type), private :: solved_parameters
      !> Bulk thermodynamics of the VV solution
      type(solvent_vv_thermodynamics_type), private :: thermo
      !> Optional store of VV states: initial guess of each solve, saves during and after it
      class(moz_checkpoint_type), allocatable, private :: checkpoint
   contains
      !> Update the potential from the solvent structure
      procedure :: update => vv_update
      !> Solve the 1D VV equations
      procedure :: solve => vv_solve
      !> Attach a checkpoint store, or detach it
      procedure :: set_checkpoint => vv_set_checkpoint
      !> Whether the VV solution is available and matches the potential
      procedure :: is_solved => vv_is_solved
      !> Copy of the radial pair of the VV solution
      procedure :: radial_grid => vv_radial_grid
      !> Copies of the physical correlations of the VV solution
      procedure :: correlations => vv_correlations
      !> Copy of the solvent susceptibility of the VV solution
      procedure :: susceptibility => vv_susceptibility
      !> Solvent susceptibility at arbitrary wave numbers
      procedure :: susceptibility_at => vv_susceptibility_at
      !> Susceptibility reconstructed from cubic-interpolated short-range c(k)
      procedure :: susceptibility_oz_at => vv_susceptibility_oz_at
      !> Long-wavelength screening 1 - 1/epsilon of the solved bulk OZ response
      procedure :: charge_screening => vv_charge_screening
      !> Canonical bulk theory name, RISM before the first solve
      procedure :: theory_name => vv_theory_name
      !> Settings of the accepted bulk solve; defaults before a current solve
      procedure :: settings => vv_settings
      !> Achieved bulk RMS residual; huge when no current solution is available
      procedure :: residual => vv_residual
      !> Bulk solution identity for response caches; -1 when stale
      procedure :: solution_generation => vv_solution_generation
      !> Bulk thermodynamics of the VV solution
      procedure :: thermodynamics => vv_thermodynamics
      !> Temperature, K
      procedure :: temperature => vv_temperature
      !> Inverse thermal energy 1/(kB T), 1/Hartree
      procedure :: beta => vv_beta
      !> Molecular number density, equal to every site density, 1/bohr^3
      procedure :: density => vv_density
      !> Number of solvent sites
      procedure :: nsite => vv_nsite
      !> Copy of the solvent structure
      procedure :: structure => vv_structure
      !> Solvent id; 0 for a custom solvent
      procedure :: solvent_id => vv_solvent_id
   end type solvent_vv_type

contains

   !> Construct a bulk solvent from a solvent system
   !>
   !> Reads the solvent id, temperature, number density
   !> (`solvent%number_density_au`) and structure (`solv_mol`, bohr) of the
   !> system; a table entry comes from `new_solvent_system`
   !>
   !> @param[out] self Solvent
   !> @param[in] ctx Run context, borrowed for the lifetime of the solvent
   !> @param[in] system Solvent system
   !> @param[out] error Missing structure or invalid state point
   subroutine new_vv_solvent(self, ctx, system, error)
      !> Solvent
      type(solvent_vv_type), intent(out) :: self
      !> Run context, borrowed for the lifetime of the solvent
      type(moist_context_type), intent(in), target :: ctx
      !> Solvent system
      type(solvent_system_type), intent(in) :: system
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(system%solv_mol)) then
         call fatal_error(error, "1D VV solvent system has no structure (solv_mol)")
         return
      end if
      if (system%solv_mol%nat < 1) then
         call fatal_error(error, "1D VV solvent requires at least one atom")
         return
      end if
      if (.not. ieee_is_finite(system%temperature)) then
         call fatal_error(error, "1D VV solvent temperature must be finite")
         return
      end if
      if (system%temperature <= 0.0_wp) then
         call fatal_error(error, "1D VV solvent temperature must be positive")
         return
      end if
      if (.not. ieee_is_finite(system%solvent%number_density_au)) then
         call fatal_error(error, "1D VV solvent number density must be finite")
         return
      end if
      if (system%solvent%number_density_au <= 0.0_wp) then
         call fatal_error(error, "1D VV solvent number density must be positive")
         return
      end if
      self%system = system
      call self%potential%refuse_host_fed_terms()
      self%ctx => ctx
      self%constructed = .true.
   end subroutine new_vv_solvent

   !> Update the potential from the solvent structure
   !>
   !> Terms keyed by solvent receive the solvent id; none for a custom solvent
   !> (id 0). Drops any VV solution
   !>
   !> @param[in,out] self Constructed solvent
   !> @param[out] error Unconstructed solvent, no terms or term failure
   subroutine vv_update(self, error)
      !> Solvent
      class(solvent_vv_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      self%solved = .false.
      if (.not. self%constructed) then
         call fatal_error(error, "Construct the 1D VV solvent before update")
         return
      end if
      if (self%system%solvent%id > 0) then
         call self%potential%update(self%system%solv_mol, error, self%system%solvent%id)
      else
         call self%potential%update(self%system%solv_mol, error)
      end if
      if (allocated(error)) return
      if (.not. associated(self%ctx)) return
      if (self%ctx%verbosity(2)) then
         call self%ctx%message("Solvent potential parameters", 2)
         call self%potential%print_table(self%ctx%unit, error)
         if (allocated(error)) return
         call self%ctx%message("", 2)
      end if
   end subroutine vv_update

   !> Solve the 1D VV equations on the uniform radial pair of the settings
   !>
   !> - w denotes the physical intramolecular matrix, not the DRISM effective matrix
   !> - Ng-split tables from the updated potential, w(k) from the rigid structure
   !>   and one site density per atom
   !> - Solver verbosity and thread count come from the run context
   !> - Iteration table and result block written at verbosity 2
   !> - Keeps the physical c, t and h, the susceptibility chi(k) = w(k) + rho h(k)
   !>   and the bulk thermodynamics on success
   !> - A stage that does not converge is an error and leaves the solvent unsolved
   !> - Timed on the context timer as `1D VV RISM` with its setup, iteration and
   !>   property phases
   !>
   !> @param[in,out] self        Updated solvent
   !> @param[out]    error       Unupdated solvent, invalid settings, setup or solve failure
   !> @param[in]     parameters  Settings; defaults when absent
   subroutine vv_solve(self, error, parameters)
      !> Solvent
      class(solvent_vv_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Settings
      type(solvent_vv_parameters_type), intent(in), optional :: parameters
      !> Timer depth on entry, restored on every return
      integer :: d0

      self%solved = .false.
      if (.not. self%constructed) then
         call fatal_error(error, "Construct the 1D VV solvent before solve")
         return
      end if
      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("1D VV RISM", category=cat_solve)
      call solve_bulk(self, error, parameters)
      call self%ctx%timer%unwind(d0)
   end subroutine vv_solve

   !> Body of vv_solve; an early return leaves its phase timer for the caller to unwind
   !>
   !> @param[in,out] self        Constructed solvent
   !> @param[out]    error       Unupdated solvent, invalid settings, setup or solve failure
   !> @param[in]     parameters  Settings; defaults when absent
   subroutine solve_bulk(self, error, parameters)
      !> Solvent
      class(solvent_vv_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Settings
      type(solvent_vv_parameters_type), intent(in), optional :: parameters
      type(solvent_vv_parameters_type) :: params
      type(vv_problem_type), target :: problem
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(packed_index_lut) :: lut
      type(potential_sites_type) :: sites
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :), uk_lr(:, :), wk(:, :), rho(:), hk(:, :)
      real(wp), allocatable :: cr_short(:, :), q(:), tail(:), wbar(:, :), zeta(:, :), zmat(:, :), zscaled(:, :)
      real(wp), allocatable :: guess(:, :)
      type(prettyprinter) :: pp
      character(len=64) :: buf
      real(wp) :: rmse, alpha
      integer :: ns, niter, i, j, a
      logical :: converged

      if (.not. self%potential%is_updated()) then
         call fatal_error(error, "Update the 1D VV solvent before solve")
         return
      end if
      if (present(parameters)) params = parameters
      call params%validate(error)
      if (allocated(error)) return
      call self%ctx%timer%start("Setup", category=cat_setup)
      call new_uniform_radial_pair(rgrid, kgrid, params%npts, params%dr, error)
      if (allocated(error)) return
      call self%potential%compute(rgrid, kgrid, u_sr, ur_lr, uk_lr, error)
      if (allocated(error)) return
      ns = self%nsite()
      call lut%init(ns)
      allocate (wk(params%npts, npair_from_ns(ns)))
      call vv_intramolecular(kgrid%r, self%system%solv_mol%xyz, lut, wk)
      rho = spread(self%density(), 1, ns)
      call self%potential%sites(sites, error)
      if (allocated(error)) return
      allocate (q(ns), source=0.0_wp)
      if (allocated(sites%q)) q = sites%q
      call new_bulk_theory(self%theory, params%theory, self%system%solv_mol%xyz, q, self%density(), &
         & self%beta(), self%system%solvent%epsilon, params%drism_smear, error)
      if (allocated(error)) return
      call self%ctx%cite(citation%rism)
      if (self%theory%name() == "DRISM") call self%ctx%cite(citation%drism)
      if (self%potential%ng_split()) call self%ctx%cite(citation%ng)
      allocate (zeta, mold=wk)
      allocate (zmat(ns, ns), zscaled(ns, ns))
      do i = 1, size(wk, 1)
         call self%theory%correction_at(kgrid%r(i), zmat, zscaled)
         do j = 1, ns
            do a = 1, j
               zeta(i, lut%idx(a, j)) = zmat(a, j)
            end do
         end do
      end do
      wbar = wk + self%density()*zeta
      call new_vv_problem(problem, self%ctx, rgrid, kgrid, wbar, u_sr, ur_lr, uk_lr, rho, self%beta(), params, &
         & error, h_offset=zeta)
      if (allocated(error)) return
      if (allocated(self%checkpoint)) call load_guess(self, params, ur_lr, guess)
      call self%ctx%timer%stop("Setup")
      call problem%solve(params%tolerance, params%maxiter, params%nlambda, converged, rmse, niter, error, &
         & guess=guess, checkpoint=self%checkpoint)
      if (allocated(error)) return
      write (buf, "(i0,' iterations, RMS residual ',es9.2)") niter, rmse
      if (.not. converged .or. .not. ieee_is_finite(rmse) .or. rmse >= params%tolerance) then
         call fatal_error(error, "1D VV solve did not converge: "//trim(buf))
         return
      end if
      call self%ctx%timer%start("Properties", category=cat_properties)
      self%rgrid = rgrid
      self%kgrid = kgrid
      self%cr = problem%cr - self%beta()*ur_lr
      self%tr = problem%tr + self%beta()*ur_lr
      self%hr = problem%cr + problem%tr
      allocate (hk, mold=wk)
      call new_radial_trafo(trafo, rgrid, kgrid, error)
      if (allocated(error)) return
      call trafo%fbt_r2k_all(self%hr, hk, error)
      if (allocated(error)) return
      self%chi = wk + self%density()*hk
      cr_short = problem%cr
      alpha = 1.0_wp
      if (self%potential%ng_split()) then
         alpha = self%potential%alpha()
      else if (any(q /= 0.0_wp)) then
         allocate (tail(rgrid%npts))
         call ng_real_tail(rgrid%r, alpha, tail)
         do j = 1, ns
            do i = 1, j
               cr_short(:, lut%idx(i, j)) = cr_short(:, lut%idx(i, j)) + self%beta()*q(i)*q(j)*tail
            end do
         end do
      end if
      call new_vv_response(self%oz_response, rgrid, kgrid, cr_short, self%system%solv_mol%xyz, q, &
         & self%density(), self%beta(), alpha, self%ctx%get_num_threads(), error, self%theory)
      if (allocated(error)) return
      call new_vv_thermodynamics(self%thermo, problem%closure, rgrid, lut, self%hr, self%cr, self%density(), &
         & self%beta(), error)
      if (allocated(error)) return
      call self%ctx%timer%stop("Properties")
      self%solved = .true.
      self%solved_generation = self%potential%generation()
      self%solved_parameters = params
      self%achieved_residual = rmse
      !$omp atomic capture
      last_solution_stamp = last_solution_stamp + 1_int64
      self%solution_stamp = last_solution_stamp
      !$omp end atomic
      if (self%ctx%verbosity(2)) then
         pp = new_prettyprinter(unit=self%ctx%unit)
         call pp%push("1D VV solution:")
         call pp%kv("Iterations", niter)
         call pp%kv("RMS residual", rmse, use_exp=.true.)
         call pp%kv("Excess chemical potential", self%thermo%chemical_potential*autokcal, "kcal/mol")
         call pp%kv("Pressure", self%thermo%pressure, "Eh/bohr^3", fmt="ES16.5")
         call pp%kv("Isothermal compressibility", self%thermo%compressibility, "bohr^3/Eh", fmt="ES16.5")
         call pp%pop()
         call pp%blank()
      end if
   end subroutine solve_bulk

   !> Set the store header of this solve and load an initial guess from the store
   !>
   !> - Header: radial grid, solvent sites, closure, theory and temperature
   !> - Loaded guess reported at verbosity 2, with the VV iteration table
   !>
   !> @param[inout] self    solvent with a store and a constructed theory
   !> @param[in]    params  validated settings of the solve
   !> @param[in]    ur_lr   long-range real-space potential, Hartree, packed (npts, npair)
   !> @param[out]   guess   short-range t_s(r), packed (npts, npair); unallocated without a guess
   subroutine load_guess(self, params, ur_lr, guess)
      !> Solvent
      class(solvent_vv_type), intent(inout) :: self
      !> Settings of the solve
      type(solvent_vv_parameters_type), intent(in) :: params
      !> Long-range real-space potential
      real(wp), intent(in) :: ur_lr(:, :)
      !> Short-range guess
      real(wp), allocatable, intent(out) :: guess(:, :)

      associate (header => self%checkpoint%header)
         header%kind = checkpoint_kind_vv
         header%grid_shape = [params%npts]
         header%spacing = params%dr
         header%solvent_sites = self%system%solv_mol%num(self%system%solv_mol%id)
         header%closure = trim(adjustl(params%closure))
         header%theory = trim(self%theory%name())
         header%temperature = self%temperature()
      end associate
      call self%checkpoint%load_guess(self%ctx, self%beta(), ur_lr, 2, guess)
   end subroutine load_guess

   !> Attach a copy of a checkpoint store, or detach the store when absent
   !>
   !> - The store supplies the initial guess of each solve and receives its states
   !> - Detached: solves read and write nothing
   !>
   !> @param[inout] self        solvent
   !> @param[in]    checkpoint  store to copy
   subroutine vv_set_checkpoint(self, checkpoint)
      !> Solvent
      class(solvent_vv_type), intent(inout) :: self
      !> Store to copy
      class(moz_checkpoint_type), intent(in), optional :: checkpoint
      if (allocated(self%checkpoint)) deallocate (self%checkpoint)
      if (present(checkpoint)) allocate (self%checkpoint, source=checkpoint)
   end subroutine vv_set_checkpoint

   !> Whether the VV solution is available and matches the potential
   !>
   !> - A potential changed after the solve (new term, update, Ng split or mixing rule)
   !>   makes the solution stale
   !> - So does a replacement by a potential with another change history
   !> - A copy of the solving potential keeps it
   !>
   !> @param[in] self        Solvent
   !> @param[in] parameters  Optional requested bulk settings and residual tolerance
   pure function vv_is_solved(self, parameters) result(solved)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Requested bulk settings, compared with the accepted solve
      type(solvent_vv_parameters_type), intent(in), optional :: parameters
      !> Solution flag
      logical :: solved
      solved = self%solved .and. self%potential%is_updated() .and. &
         & self%potential%generation() == self%solved_generation
      if (present(parameters)) then
         solved = solved .and. self%solved_parameters%same_state(parameters) .and. &
            & self%achieved_residual < parameters%tolerance
      end if
   end function vv_is_solved

   !> Canonical name of the bulk theory; RISM before the first solve
   !>
   !> @param[in] self  bulk solvent
   pure function vv_theory_name(self) result(name)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Canonical theory name
      character(len=8) :: name
      name = "RISM"
      if (allocated(self%theory)) name = self%theory%name()
   end function vv_theory_name

   !> Copy the accepted bulk settings for a model borrowing a solved solvent
   !>
   !> @param[in] self  bulk solvent
   function vv_settings(self) result(parameters)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Accepted settings; defaults before a current solve
      type(solvent_vv_parameters_type) :: parameters
      if (self%is_solved()) parameters = self%solved_parameters
   end function vv_settings

   !> Achieved residual of the accepted solution, with staleness accounted for
   !>
   !> @param[in] self  bulk solvent
   pure function vv_residual(self) result(rmse)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Achieved RMS residual; huge when no current solution is available
      real(wp) :: rmse
      rmse = huge(1.0_wp)
      if (self%is_solved()) rmse = self%achieved_residual
   end function vv_residual

   !> Identity of the bulk solution; copies keep it, each solve takes a new one
   !>
   !> @param[in] self  bulk solvent
   pure function vv_solution_generation(self) result(generation)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Solution stamp; -1 when stale
      integer(int64) :: generation
      generation = -1_int64
      if (self%is_solved()) generation = self%solution_stamp
   end function vv_solution_generation

   !> Reconstruct OZ response from short-range cbar on arbitrary wave numbers
   !>
   !> - Optional charge response holds chi q / k^2, including its k -> 0 limit
   !>
   !> @param[in]  self             solved solvent
   !> @param[in]  k                wave numbers, 1/bohr (nk)
   !> @param[out] chi              susceptibility, packed (nk, npair)
   !> @param[out] error            no current VV solution or response failure
   !> @param[out] charge_response  chi q / k^2 (nk, ns), optional
   subroutine vv_susceptibility_oz_at(self, k, chi, error, charge_response)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Wave numbers, 1/bohr
      real(wp), intent(in) :: k(:)
      !> Susceptibility
      real(wp), allocatable, intent(out) :: chi(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Charge response chi q / k^2
      real(wp), allocatable, intent(out), optional :: charge_response(:, :)
      if (.not. self%is_solved()) then
         call fatal_error(error, "Solve the 1D VV solvent before requesting the OZ response")
         return
      end if
      call self%oz_response%at(k, chi, error, charge_response)
   end subroutine vv_susceptibility_oz_at

   !> Dielectric screening from the same analytical limit as the OZ susceptibility
   !>
   !> @param[in]  self       solved solvent
   !> @param[out] screening  long-wavelength screening, dimensionless
   !> @param[out] error      no current VV solution or response failure
   subroutine vv_charge_screening(self, screening, error)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Screening 1 - 1/epsilon
      real(wp), intent(out) :: screening
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      screening = 0.0_wp
      if (.not. self%is_solved()) then
         call fatal_error(error, "Solve the 1D VV solvent before requesting dielectric screening")
         return
      end if
      call self%oz_response%charge_screening(screening, error)
   end subroutine vv_charge_screening

   !> Copy of the radial pair of the VV solution
   !>
   !> @param[in]  self   Solved solvent
   !> @param[out] rgrid  Real-space grid
   !> @param[out] kgrid  Reciprocal grid
   !> @param[out] error  No current VV solution
   subroutine vv_radial_grid(self, rgrid, kgrid, error)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Radial pair
      type(moist_math_grid_radial_type), intent(out) :: rgrid, kgrid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%is_solved()) then
         call fatal_error(error, "1D VV solvent has no current solution")
         return
      end if
      rgrid = self%rgrid
      kgrid = self%kgrid
   end subroutine vv_radial_grid

   !> Copies of the physical correlations of the VV solution, packed (npts, npair)
   !>
   !> @param[in]  self   Solved solvent
   !> @param[out] cr     Direct correlation c(r)
   !> @param[out] tr     Indirect correlation t(r) = h(r) - c(r)
   !> @param[out] hr     Total correlation h(r)
   !> @param[out] error  No current VV solution
   subroutine vv_correlations(self, cr, tr, hr, error)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Physical correlations
      real(wp), allocatable, intent(out) :: cr(:, :), tr(:, :), hr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%is_solved()) then
         call fatal_error(error, "1D VV solvent has no current solution")
         return
      end if
      cr = self%cr
      tr = self%tr
      hr = self%hr
   end subroutine vv_correlations

   !> Copy of the solvent susceptibility chi(k) = w(k) + rho h(k), packed (npts, npair)
   !>
   !> The k-space response the UV models consume, on the reciprocal grid of radial_grid
   !>
   !> @param[in]  self   Solved solvent
   !> @param[out] chi    Susceptibility
   !> @param[out] error  No current VV solution
   subroutine vv_susceptibility(self, chi, error)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Susceptibility
      real(wp), allocatable, intent(out) :: chi(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%is_solved()) then
         call fatal_error(error, "1D VV solvent has no current solution")
         return
      end if
      chi = self%chi
   end subroutine vv_susceptibility

   !> Solvent susceptibility chi(k) = w(k) + rho h(k) at arbitrary wave numbers, packed (nk, npair)
   !>
   !> - w(k) exactly from the rigid structure; h(k) by the radial transform of
   !>   the stored h(r) at each k (vv_transform_at), the quadrature the solution's
   !>   DST-IV pair evaluates at its nodes, with h(r) = 0 beyond the VV grid
   !> - The complete native k-grid in its stored order reuses chi without quadrature;
   !>   other requests cost one pass over the r grid per requested wave number
   !> - For 0 <= k <= pi/dr of the VV grid; direct quadrature retains radial
   !>   discretization, finite-cutoff and bulk-solve errors
   !>
   !> @param[in]  self   Solved solvent
   !> @param[in]  k      wave numbers, 1/bohr (nk)
   !> @param[out] chi    susceptibility, packed (nk, npair)
   !> @param[out] error  no current VV solution, or a nonfinite/out-of-band wave number
   subroutine vv_susceptibility_at(self, k, chi, error)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Wave numbers
      real(wp), intent(in) :: k(:)
      !> Susceptibility
      real(wp), allocatable, intent(out) :: chi(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(packed_index_lut) :: lut
      real(wp), allocatable :: wk(:, :)
      real(wp) :: band_edge
      integer :: nthreads

      if (.not. self%is_solved()) then
         call fatal_error(error, "1D VV solvent has no current solution")
         return
      end if
      if (.not. all(ieee_is_finite(k))) then
         call fatal_error(error, "1D VV susceptibility needs finite wave numbers")
         return
      end if
      if (size(k) == size(self%kgrid%r)) then
         if (all(k == self%kgrid%r)) then
            chi = self%chi
            return
         end if
      end if
      band_edge = real(size(self%kgrid%r), wp)*self%kgrid%spacing
      if (size(k) > 0) then
         if (.not. (minval(k) >= 0.0_wp .and. maxval(k) <= band_edge)) then
            call fatal_error(error, "1D VV susceptibility is defined for wave numbers in [0, pi/dr] of the VV grid")
            return
         end if
      end if
      nthreads = 0
      if (associated(self%ctx)) nthreads = self%ctx%get_num_threads()
      allocate (chi(size(k), size(self%chi, 2)))
      call vv_transform_at(self%rgrid%r, self%rgrid%w, self%hr, k, chi, nthreads)
      chi = self%density()*chi
      call lut%init(self%nsite())
      allocate (wk, mold=chi)
      call vv_intramolecular(k, self%system%solv_mol%xyz, lut, wk)
      chi = chi + wk
   end subroutine vv_susceptibility_at

   !> Bulk thermodynamics of the VV solution
   !>
   !> @param[in]  self    Solved solvent
   !> @param[out] thermo  Excess chemical potential, pressure, compressibility
   !> @param[out] error   No current VV solution
   subroutine vv_thermodynamics(self, thermo, error)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Thermodynamics
      type(solvent_vv_thermodynamics_type), intent(out) :: thermo
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%is_solved()) then
         call fatal_error(error, "1D VV solvent has no current solution")
         return
      end if
      thermo = self%thermo
   end subroutine vv_thermodynamics

   !> Temperature, K; 0 before construction
   !>
   !> @param[in] self Solvent
   pure function vv_temperature(self) result(temperature)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Temperature, K
      real(wp) :: temperature
      temperature = 0.0_wp
      if (self%constructed) temperature = self%system%temperature
   end function vv_temperature

   !> Inverse thermal energy 1/(kB T), 1/Hartree; 0 before construction
   !>
   !> @param[in] self Solvent
   pure function vv_beta(self) result(beta)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Inverse thermal energy, 1/Hartree
      real(wp) :: beta
      beta = 0.0_wp
      if (self%constructed) then
         beta = 1.0_wp/(Boltzmann_constant/atomic_unit_of_energy*self%system%temperature)
      end if
   end function vv_beta

   !> Molecular number density, equal to every site density, 1/bohr^3; 0 before construction
   !>
   !> @param[in] self Solvent
   pure function vv_density(self) result(density)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Number density, 1/bohr^3
      real(wp) :: density
      density = 0.0_wp
      if (self%constructed) density = self%system%solvent%number_density_au
   end function vv_density

   !> Number of solvent sites, one per atom; 0 before construction
   !>
   !> @param[in] self Solvent
   pure function vv_nsite(self) result(nsite)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Site count
      integer :: nsite
      nsite = 0
      if (self%constructed) nsite = self%system%solv_mol%nat
   end function vv_nsite

   !> Copy of the solvent structure, bohr
   !>
   !> @param[in] self Constructed solvent
   function vv_structure(self) result(mol)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Solvent structure; empty before construction
      type(structure_type) :: mol
      if (self%constructed) mol = self%system%solv_mol
   end function vv_structure

   !> Solvent id; 0 for a custom solvent or before construction
   !>
   !> @param[in] self Solvent
   pure function vv_solvent_id(self) result(id)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Solvent id
      integer :: id
      id = 0
      if (self%constructed) id = self%system%solvent%id
   end function vv_solvent_id

end module moist_model_moz_solvent_vv
