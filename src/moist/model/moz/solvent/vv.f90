!> Bulk solvent of the MOZ models
!>
!> Not a solvation model: holds the solvent system (solvent id, temperature,
!> number density, structure), the Lennard-Jones mixing rule of every
!> interaction with this solvent and the solvent potential set, which refuses
!> host-fed terms. One site per atom. A custom solvent is a hand-filled
!> `solvation_system_type` (solvent id 0, temperature, number density and
!> structure). Will own the 1D VV solve and the solvent susceptibility; the UV
!> models copy it
module moist_model_moz_solvent_vv
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_data_solvents, only: solvation_system_type
   use moist_math_packed_sym, only: packed_index, npair_from_ns
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_model_moz_potential_base, only: potential_type
   use moist_model_moz_potential_set, only: potential_set_type
   use moist_model_moz_potential_lj_base, only: lj_mixing_lorentz_berthelot, check_lj_mixing
   use moist_model_moz_potential_kernel_category, only: require_split_alpha
   use moist_model_moz_potential_kernel_table_1d, only: assemble_table_1d
   implicit none(type, external)
   private

   public :: solvent_vv_type, new_vv_solvent

   !> Bulk solvent: system, mixing rule and potential set
   type :: solvent_vv_type
      !> Run context, borrowed for the lifetime of the solvent
      type(moist_context_type), pointer :: ctx => null()
      !> Solvent system: id, temperature, number density and structure (bohr)
      type(solvation_system_type), private :: system
      !> Solvent potential set; refuses host-fed terms. Read by the UV models
      !> and the kernels; add terms through `add_potential` only
      type(potential_set_type) :: set
      !> Lennard-Jones mixing rule of every interaction with this solvent
      integer, private :: mixing_rule = lj_mixing_lorentz_berthelot
      !> Whether the solvent was constructed
      logical, private :: constructed = .false.
      !> Whether the VV solution is available
      logical, private :: solved = .false.
   contains
      !> Append a copy of a potential term
      procedure :: add_potential => vv_add_potential
      !> Build the potential set from the solvent structure
      procedure :: build => vv_build
      !> Whether the potential set is built
      procedure :: is_built => vv_is_built
      !> Packed VV site pairs of the upper triangle
      procedure :: site_pairs => vv_site_pairs
      !> Ng-split VV potential tables
      procedure :: potential_tables => vv_potential_tables
      !> Solve the 1D VV equations; pending
      procedure :: solve => vv_solve
      !> Whether the VV solution is available
      procedure :: is_solved => vv_is_solved
      !> Number of solvent sites
      procedure :: nsite => vv_nsite
      !> Copy of the solvent structure
      procedure :: structure => vv_structure
      !> Solvent id; 0 for a custom solvent
      procedure :: solvent_id => vv_solvent_id
      !> Temperature, K
      procedure :: temperature => vv_temperature
      !> Number density, 1/bohr^3
      procedure :: density => vv_density
      !> Lennard-Jones mixing rule
      procedure :: mixing => vv_mixing
   end type solvent_vv_type

contains

   !> Construct a bulk solvent from a solvation system
   !>
   !> Reads the solvent id, temperature, number density
   !> (`solvent_number_density_au`) and structure (`solv_mol`, bohr) of the
   !> system; a table entry comes from `new_solvation_system`
   !>
   !> @param[out] self Solvent
   !> @param[in] ctx Run context, borrowed for the lifetime of the solvent
   !> @param[in] system Solvent system
   !> @param[out] error Missing structure, invalid state point or unknown mixing rule
   !> @param[in] mixing Optional Lennard-Jones mixing rule; Lorentz-Berthelot when absent
   subroutine new_vv_solvent(self, ctx, system, error, mixing)
      !> Solvent
      type(solvent_vv_type), intent(out) :: self
      !> Run context, borrowed for the lifetime of the solvent
      type(moist_context_type), intent(in), target :: ctx
      !> Solvent system
      type(solvation_system_type), intent(in) :: system
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional Lennard-Jones mixing rule
      integer, intent(in), optional :: mixing
      if (present(mixing)) then
         call check_lj_mixing(mixing, error)
         if (allocated(error)) return
         self%mixing_rule = mixing
      end if
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
      if (.not. ieee_is_finite(system%solvent_number_density_au)) then
         call fatal_error(error, "1D VV solvent number density must be finite")
         return
      end if
      if (system%solvent_number_density_au <= 0.0_wp) then
         call fatal_error(error, "1D VV solvent number density must be positive")
         return
      end if
      self%system = system
      self%set%refuse_host_fed = .true.
      self%ctx => ctx
      self%constructed = .true.
   end subroutine new_vv_solvent

   !> Append a copy of a potential term; host-fed terms are refused
   !>
   !> Drops any VV solution; the set needs a new build
   !>
   !> @param[in,out] self Solvent
   !> @param[in] term Term to copy
   !> @param[out] error Host-fed term or copy failure
   subroutine vv_add_potential(self, term, error)
      !> Solvent
      class(solvent_vv_type), intent(inout) :: self
      !> Term to copy
      class(potential_type), intent(in) :: term
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      self%set%refuse_host_fed = .true.
      call self%set%add(term, error)
      if (allocated(error)) return
      self%solved = .false.
   end subroutine vv_add_potential

   !> Build the potential set from the solvent structure
   !>
   !> Terms keyed by solvent receive the solvent id; none for a custom solvent (id 0)
   !>
   !> @param[in,out] self Constructed solvent
   !> @param[out] error Unconstructed solvent, empty set or term failure
   subroutine vv_build(self, error)
      !> Solvent
      class(solvent_vv_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      self%solved = .false.
      if (.not. self%constructed) then
         call fatal_error(error, "Construct the 1D VV solvent before build")
         return
      end if
      if (self%set%n_terms() < 1) then
         call fatal_error(error, "1D VV solvent has no potential terms; add one before build")
         return
      end if
      if (self%system%solvent_id > 0) then
         call self%set%build(self%system%solv_mol, error, self%system%solvent_id)
      else
         call self%set%build(self%system%solv_mol, error)
      end if
      if (allocated(error)) return
      if (.not. associated(self%ctx)) return
      if (self%ctx%writes(2)) then
         call self%ctx%message("Solvent potential parameters", 2)
         call self%set%print_table(self%system%solv_mol, self%ctx%unit, error)
         if (allocated(error)) return
         call self%ctx%message("", 2)
      end if
   end subroutine vv_build

   !> Whether the potential set is built
   !>
   !> @param[in] self Solvent
   pure function vv_is_built(self) result(built)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Build flag
      logical :: built
      built = self%set%natom > 0
   end function vv_is_built

   !> Packed VV site pairs of the upper triangle
   !>
   !> Column `packed_index(i, j)` holds `[i, j]` with i <= j
   !>
   !> @param[in] self Solvent
   !> @param[out] pairs Site pairs (2, nsite (nsite + 1)/2)
   subroutine vv_site_pairs(self, pairs)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Site pairs
      integer, allocatable, intent(out) :: pairs(:, :)
      integer :: i, j, ns
      ns = self%nsite()
      allocate (pairs(2, npair_from_ns(ns)))
      do j = 1, ns
         do i = 1, j
            pairs(:, packed_index(i, j)) = [i, j]
         end do
      end do
   end subroutine vv_site_pairs

   !> Ng-split VV potential tables of the packed site pairs
   !>
   !> @param[in] self Built solvent
   !> @param[in] rgrid Real-space radial grid
   !> @param[in] kgrid Reciprocal radial grid
   !> @param[in] ng_split Whether the Coulomb tail is split off
   !> @param[in] alpha Ng split parameter, 1/bohr; positive and finite when the split is on
   !> @param[out] u_sr Short-range potential, Hartree (npts, npair)
   !> @param[out] ur_lr Long-range real-space potential, Hartree (npts, npair)
   !> @param[out] uk_lr Long-range reciprocal potential, Hartree bohr^3 (nk, npair)
   !> @param[out] error Unbuilt set, invalid split parameter or kernel failure
   subroutine vv_potential_tables(self, rgrid, kgrid, ng_split, alpha, u_sr, ur_lr, uk_lr, error)
      !> Built solvent
      class(solvent_vv_type), intent(in) :: self
      !> Real-space radial grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Ng split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> Short-range potential (npts, npair)
      real(wp), allocatable, intent(out) :: u_sr(:, :)
      !> Long-range real-space potential (npts, npair)
      real(wp), allocatable, intent(out) :: ur_lr(:, :)
      !> Long-range reciprocal potential (nk, npair)
      real(wp), allocatable, intent(out) :: uk_lr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: pairs(:, :)
      call self%set%require_built(error)
      if (allocated(error)) return
      call require_split_alpha(ng_split, alpha, error)
      if (allocated(error)) return
      call self%site_pairs(pairs)
      call assemble_table_1d(rgrid, kgrid, self%set, self%set, pairs, self%mixing_rule, ng_split, alpha, &
         & u_sr, ur_lr, uk_lr, error)
   end subroutine vv_potential_tables

   !> Solve the 1D VV equations; pending
   !>
   !> @param[in,out] self Built solvent
   !> @param[out] error Always set
   subroutine vv_solve(self, error)
      !> Solvent
      class(solvent_vv_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      self%solved = .false.
      call fatal_error(error, "1D VV solve is pending")
   end subroutine vv_solve

   !> Whether the VV solution is available
   !>
   !> @param[in] self Solvent
   pure function vv_is_solved(self) result(solved)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Solution flag
      logical :: solved
      solved = self%solved
   end function vv_is_solved

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
      if (self%constructed) id = self%system%solvent_id
   end function vv_solvent_id

   !> Temperature, K
   !>
   !> @param[in] self Constructed solvent
   function vv_temperature(self) result(temperature)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Temperature, K
      real(wp) :: temperature
      if (.not. self%constructed) error stop "1D VV solvent temperature read before construction"
      temperature = self%system%temperature
   end function vv_temperature

   !> Number density, 1/bohr^3
   !>
   !> @param[in] self Constructed solvent
   function vv_density(self) result(density)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Number density, 1/bohr^3
      real(wp) :: density
      if (.not. self%constructed) error stop "1D VV solvent density read before construction"
      density = self%system%solvent_number_density_au
   end function vv_density

   !> Lennard-Jones mixing rule of every interaction with this solvent
   !>
   !> @param[in] self Solvent
   pure function vv_mixing(self) result(mixing)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Mixing rule, `lj_mixing_lorentz_berthelot` or `lj_mixing_geometric`
      integer :: mixing
      mixing = self%mixing_rule
   end function vv_mixing

end module moist_model_moz_solvent_vv
