!> CDS energy contribution for SMD
!>
!> Computes the cavity/dispersion/solvent (CDS) contribution as in SMD.
!> Input surface areas are the per-atom SASA obtained from the iSwiG
!> cavity generator. Energies are accumulated from atom-resolved
!> surface tensions (sigma_k) and a solvent macroscopic term.
module moist_model_component_cds_calculator
   use mctc_env, only: wp
   use mctc_io_convert, only: autoaa
   use moist_model_component_cds_sigma, only: cds_surfacetensions
   implicit none
   private
   public :: calc_cds

contains

   !> Compute CDS terms (returns cal/mol)
   subroutine calc_cds(surft, surface, cds, cds_sm)
      !> Per-atom and solvent surface tensions
      type(cds_surfacetensions), intent(in) :: surft
      !> Per-atom SASA (iSwiG) in bohr^2
      real(wp), intent(in) :: surface(:)
      !> CDS contribution per atom (cal/mol)
      real(wp), allocatable, intent(out) :: cds(:)
      !> Solvent CDS contribution (cal/mol)
      real(wp), intent(out) :: cds_sm
      integer :: i, nat

      nat = size(surface)
      allocate (cds(nat))
      cds = 0.0_wp
      ! Convert area from bohr^2 to Å^2 (autoaa^2) to match sigma_k units
      cds_sm = surft%sm*sum(surface)*autoaa**2
      do i = 1, nat
         cds(i) = surft%sk(i)*surface(i)*autoaa**2
      end do
   end subroutine calc_cds

end module moist_model_component_cds_calculator
