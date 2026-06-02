module moist_model_component_cds
   use mctc_env, only: wp, error_type
   use moist_type, only: solvation_model_component, cavity_type
   use moist_type, only: potential_type
   use mctc_io_convert, only: autokcal
   use mctc_io, only: structure_type
   use moist_model_component_cds_sigma, only: cds_surfacetensions
   use moist_model_component_cds_calculator, only: calc_cds
   implicit none
!    private
!    public :: moist_model_component_cds, new_moist_model_component_cds

!    type, extends(solvation_model_component) :: moist_model_component_cds
!       class(structure_type), allocatable :: mol
!       class(cavity_type), allocatable :: cavity
!       type(cds_surfacetensions), allocatable :: surft
!    contains
!       procedure :: update => moist_model_component_cds_update
!       procedure :: get_energy => moist_model_component_cds_get_energy
!       procedure :: get_potential => moist_model_component_cds_get_potential
!       procedure :: get_gradient => moist_model_component_cds_get_gradient
!       procedure :: set_surft => moist_model_component_cds_set_surft
!    end type moist_model_component_cds

!    interface new_moist_model_component_cds
!       module procedure new_moist_model_component_cds
!    end interface

! contains

!    subroutine new_moist_model_component_cds(error, self)
!       type(error_type), allocatable, intent(out) :: error
!       type(moist_model_component_cds), intent(out) :: self
!       allocate (character(len=1) :: self%name)
!       self%name = 'cds'
!    end subroutine new_moist_model_component_cds

!    subroutine moist_model_component_cds_set_surft(self, surft_in)
!       class(moist_model_component_cds), intent(inout) :: self
!       type(cds_surfacetensions), intent(in) :: surft_in
!       if (allocated(self%surft)) deallocate (self%surft)
!       allocate (self%surft)
!       self%surft = surft_in
!    end subroutine moist_model_component_cds_set_surft

!    subroutine moist_model_component_cds_update(self, mol, cavity)
!       class(moist_model_component_cds), intent(inout) :: self
!       type(structure_type), intent(in) :: mol
!       class(cavity_type), intent(inout) :: cavity
!       self%mol = mol
!       self%cavity = cavity
!    end subroutine moist_model_component_cds_update

!    subroutine moist_model_component_cds_get_energy(self, energy)
!       class(moist_model_component_cds), intent(inout) :: self
!       real(wp), intent(inout) :: energy(:)
!       real(wp), allocatable :: cds(:)
!       real(wp) :: cds_sm, cds_total_cal, cds_total_kcal

!       if (.not. allocated(self%surft) .or. .not. allocated(self%cavity%a)) then
!          energy(:) = 0.0_wp
!          return
!       end if

!       call calc_cds(self%surft, self%cavity%a, cds, cds_sm)
!       cds_total_cal = sum(cds) + cds_sm
!       cds_total_kcal = cds_total_cal/1000.0_wp
!       energy(1) = cds_total_kcal/autokcal
!       if (allocated(cds)) deallocate (cds)
!    end subroutine moist_model_component_cds_get_energy

!    subroutine moist_model_component_cds_get_potential(self, potential)
!       class(moist_model_component_cds), intent(inout) :: self
!       type(potential_type), intent(inout) :: potential
!       ! not implemented
!    end subroutine moist_model_component_cds_get_potential

!    subroutine moist_model_component_cds_get_gradient(self, gradient)
!       class(moist_model_component_cds), intent(inout) :: self
!       real(wp), intent(inout) :: gradient(:, :)
!       gradient = 0.0_wp
!    end subroutine moist_model_component_cds_get_gradient

end module moist_model_component_cds
