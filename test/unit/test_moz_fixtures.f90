!> Shared fixtures of the MOZ tests (not a suite)
!>
!> - Table water solvent with custom LJ and fixed charges
!> - Pair-only term with switchable declaration failure
!> - Reference 12-6 potential with Lorentz-Berthelot mixing
!> - Answers to pending atomic_charges and radial_potential requests
module test_moz_fixtures
   use mctc_env, only: wp, moist_error => error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_channels_coupling, only: coupling_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_data_solvents, only: solvation_system_type, new_solvation_system
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent
   use moist_model_moz_potential_term, only: potential_term_type, potential_name_len
   use moist_model_moz_potential_lj_base, only: new_lj_12_6
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use moist_model_moz_potential_electrostatic_fixed, only: fixed_charges_type, new_fixed_charges
   implicit none(type, external)
   private

   public :: water_id, water_sigma, water_epsilon, water_charges
   public :: new_water_system, new_water_solvent, refusing_term_type, ref_lj, answer_charges, answer_radial_potential

   !> Water solvent id of the table
   integer, parameter :: water_id = 175
   !> Distinct per-site sigma, bohr; swapped pair columns detectable
   real(wp), parameter :: water_sigma(3) = [6.0_wp, 1.5_wp, 2.0_wp]
   !> Epsilon per site, Hartree
   real(wp), parameter :: water_epsilon(3) = [2.4e-4_wp, 1.0e-5_wp, 2.0e-5_wp]
   !> Charges per site, e
   real(wp), parameter :: water_charges(3) = [-0.8_wp, 0.3_wp, 0.5_wp]

   !> Pair-only test term whose declaration can be made to fail
   !>
   !> - Pointer-shared switch; copied term follows test flag inside model
   type, extends(potential_term_type) :: refusing_term_type
      !> Whether `declare_pass` fails; never when unassociated
      logical, pointer :: refuse => null()
   contains
      procedure :: name => refusing_name
      procedure :: update => refusing_update
      procedure :: declare_pass_1d => refusing_declare_pass_1d
      procedure :: declare_pass_3d => refusing_declare_pass_3d
   end type refusing_term_type

contains

   !> Construct table water system
   !>
   !> @param[out] system  water solvation system, structure in bohr
   !> @param[out] err     table failure
   subroutine new_water_system(system, err)
      !> Water solvation system
      type(solvation_system_type), intent(out) :: system
      !> Table failure
      type(moist_error), allocatable, intent(out) :: err
      call new_solvation_system(system, water_id, error=err)
   end subroutine new_water_system

   !> Updated table water solvent with custom LJ and fixed charges
   !>
   !> @param[in]  ctx      run context, outlives the solvent
   !> @param[out] solvent  updated solvent
   !> @param[out] err      construction or update error
   !> @param[in]  mixing   optional Lennard-Jones mixing rule; Lorentz-Berthelot when absent
   subroutine new_water_solvent(ctx, solvent, err, mixing)
      !> Run context, outlives the solvent
      type(moist_context_type), intent(in), target :: ctx
      !> Updated solvent
      type(solvent_vv_type), intent(out) :: solvent
      !> Construction or update error
      type(moist_error), allocatable, intent(out) :: err
      !> Optional Lennard-Jones mixing rule
      integer, intent(in), optional :: mixing
      type(solvation_system_type) :: system
      type(custom_lj_type) :: lj
      type(fixed_charges_type) :: charges
      call new_water_system(system, err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      if (present(mixing) .and. .not. allocated(err)) call solvent%potential%set_mixing(mixing, err)
      if (.not. allocated(err)) call new_custom_lj(lj, water_sigma, water_epsilon, err)
      if (.not. allocated(err)) call new_fixed_charges(charges, water_charges, err)
      if (.not. allocated(err)) call solvent%potential%add(lj, err)
      if (.not. allocated(err)) call solvent%potential%add(charges, err)
      if (.not. allocated(err)) call solvent%update(err)
   end subroutine new_water_solvent

   !> Reference 12-6 potential with Lorentz-Berthelot mixing
   !>
   !> @param[in] r   distance, bohr
   !> @param[in] e1  epsilon of site 1, Hartree
   !> @param[in] s1  sigma of site 1, bohr
   !> @param[in] e2  epsilon of site 2, Hartree
   !> @param[in] s2  sigma of site 2, bohr
   elemental function ref_lj(r, e1, s1, e2, s2) result(u)
      !> Distance, bohr
      real(wp), intent(in) :: r
      !> Epsilon of site 1
      real(wp), intent(in) :: e1
      !> Sigma of site 1
      real(wp), intent(in) :: s1
      !> Epsilon of site 2
      real(wp), intent(in) :: e2
      !> Sigma of site 2
      real(wp), intent(in) :: s2
      !> Potential, Hartree
      real(wp) :: u
      real(wp) :: sr6
      sr6 = (0.5_wp*(s1 + s2)/r)**6
      u = 4.0_wp*sqrt(e1*e2)*(sr6*sr6 - sr6)
   end function ref_lj

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function refusing_name(self) result(name)
      !> Term
      class(refusing_term_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "refusing"
   end function refusing_name

   !> Fill Lennard-Jones data covering every atom
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         structure of this side
   !> @param[out]    error       parameter failure
   !> @param[in]     solvent_id  unused
   subroutine refusing_update(self, mol, error, solvent_id)
      !> Term
      class(refusing_term_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(moist_error), allocatable, intent(out) :: error
      !> Unused
      integer, intent(in), optional :: solvent_id
      real(wp), allocatable :: sigma(:), epsilon(:)
      allocate (sigma(mol%nat), source=3.0_wp)
      allocate (epsilon(mol%nat), source=1.0e-4_wp)
      if (.not. allocated(self%pair)) allocate (self%pair)
      call new_lj_12_6(self%pair, sigma, epsilon, error)
   end subroutine refusing_update

   !> Declare nothing for a radial grid, or fail when switched to refuse
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      radial grid
   !> @param[in,out] coupling  coupling, scope set by the model
   !> @param[out]    error     refusal
   subroutine refusing_declare_pass_1d(self, grid, coupling, error)
      !> Term
      class(refusing_term_type), intent(inout) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Coupling, scope set by the model
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(moist_error), allocatable, intent(out) :: error
      call refuse_if_switched(self, error)
   end subroutine refusing_declare_pass_1d

   !> Declare nothing for a volume grid, or fail when switched to refuse
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      volume grid
   !> @param[in,out] coupling  coupling, scope set by the model
   !> @param[out]    error     refusal
   subroutine refusing_declare_pass_3d(self, grid, coupling, error)
      !> Term
      class(refusing_term_type), intent(inout) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Coupling, scope set by the model
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(moist_error), allocatable, intent(out) :: error
      call refuse_if_switched(self, error)
   end subroutine refusing_declare_pass_3d

   !> Fail when the term is switched to refuse
   !>
   !> @param[in]  self   term
   !> @param[out] error  refusal
   subroutine refuse_if_switched(self, error)
      !> Term
      class(refusing_term_type), intent(in) :: self
      !> Error handling
      type(moist_error), allocatable, intent(out) :: error
      if (.not. associated(self%refuse)) return
      if (self%refuse) call fatal_error(error, "Test term: refused declaration")
   end subroutine refuse_if_switched

   !> Answer the pending atomic charges of a staged walk
   !>
   !> @param[in,out] coupling  staged coupling
   !> @param[in]     q         charges per atom
   !> @param[out]    err       missing request or rejected answer
   subroutine answer_charges(coupling, q, err)
      !> Staged coupling
      type(coupling_type), intent(inout) :: coupling
      !> Charges per atom
      real(wp), intent(in) :: q(:)
      !> Missing request or rejected answer
      type(moist_error), allocatable, intent(out) :: err
      do while (coupling%next())
         associate (item => coupling%request())
            if (item%name() /= "atomic_charges") cycle
         end associate
         call coupling%answer("q", q, err)
         return
      end do
      call fatal_error(err, "atomic_charges is not pending")
   end subroutine answer_charges

   !> Answer the pending radial potential of a staged walk
   !>
   !> @param[in,out] coupling  staged coupling
   !> @param[in]     phi       radial potential per atom (ngrid, natom)
   !> @param[out]    err       missing request or rejected answer
   subroutine answer_radial_potential(coupling, phi, err)
      !> Staged coupling
      type(coupling_type), intent(inout) :: coupling
      !> Radial potential per atom (ngrid, natom)
      real(wp), intent(in) :: phi(:, :)
      !> Missing request or rejected answer
      type(moist_error), allocatable, intent(out) :: err
      do while (coupling%next())
         associate (item => coupling%request())
            if (item%name() /= "radial_potential") cycle
         end associate
         call coupling%answer("phi", phi, err)
         return
      end do
      call fatal_error(err, "radial_potential is not pending")
   end subroutine answer_radial_potential

end module test_moz_fixtures
