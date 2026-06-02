!> SMD CDS parameter initialization for moist
!>
!> Initializes the subset of SMD parameters required for the CDS
!> (cavity/dispersion/solvent) term using available solvent system
!> properties. For water, a dedicated parameter set is used.
module moist_model_component_cds_parameters
   use mctc_env, only: wp
   use iso_fortran_env, only: output_unit
   use mctc_env_error, only: error_type
   use moist_data_solvents, only: solvation_system_parameters
   implicit none
   private

   public :: smd_param, new_smd_model_parameters

   integer, parameter :: max_elem = 118

   !> SMD parameters required for sigma_k and solvent term
   type :: smd_param
      !> Element-dependent sigma term
      real(wp) :: zk(max_elem)
      !> Pair-dependent sigma term
      real(wp) :: zkk(max_elem, max_elem)
      !> Switching-function pair radii (Å)
      real(wp) :: rzkk(max_elem, max_elem)
      !> Switching-function pair widths (Å)
      real(wp) :: drzkk(max_elem, max_elem)
      !> Additional nonclassical NC3 parameter and switching
      real(wp) :: nc3
      real(wp) :: rnc3
      real(wp) :: drnc3
      !> Macroscopic solvent term coefficient
      real(wp) :: s_m
      !> Solvent alpha (used for oxygen radius scaling in some variants)
      real(wp) :: alpha
   end type smd_param

contains

   !> Initialize SMD parameters using solvent system properties
   subroutine new_smd_model_parameters(self, system)
      type(smd_param), intent(out) :: self
      type(solvation_system_parameters), intent(in) :: system

      ! local tables mirroring numsa defaults (trimmed to what sigma.f90 uses)
      real(wp) :: ref_zk_h2o(max_elem)
      real(wp) :: ref_zkk_h2o(max_elem, max_elem)
      real(wp) :: ref_zk(3, max_elem)
      real(wp) :: ref_zkk(3, max_elem, max_elem)
      real(wp) :: ref_rzkk(max_elem, max_elem)
      real(wp) :: ref_drzkk(max_elem, max_elem)
      real(wp) :: ref_nc3, ref_rnc3, ref_drnc3
      real(wp) :: ref_sg, ref_sr2, ref_sp2, ref_sb2

      integer :: Z, Z2
      logical :: is_water
      real(wp) :: n, alpha, beta, msurft_mN, arom, fclbr

      call init_default_tables(.false., ref_zk_h2o, ref_zkk_h2o, ref_zk, ref_zkk, &
                               ref_rzkk, ref_drzkk, ref_nc3, ref_rnc3, ref_drnc3, ref_sg, ref_sr2, ref_sp2, ref_sb2)

      ! Simple water detection by name
      is_water = index(system%solvent_name, 'water') > 0 .or. index(system%solvent_name, 'h2o') > 0

      if (is_water) then
         ! H2O special: s_m = 0.0 and dedicated zk/zkk set
         self%alpha = 0.82_wp
         self%zk = 0.0_wp; self%zkk = 0.0_wp
         do Z = 1, max_elem
            self%zk(Z) = ref_zk_h2o(Z)
         end do
         do Z = 1, max_elem
            do Z2 = 1, max_elem
               self%zkk(Z, Z2) = ref_zkk_h2o(Z, Z2)
            end do
         end do
         self%rzkk = ref_rzkk
         self%drzkk = ref_drzkk
         self%nc3 = ref_nc3
         self%rnc3 = ref_rnc3
         self%drnc3 = ref_drnc3
         self%s_m = 0.0_wp
         return
      end if

      ! Other solvents: linear combinations in (n, alpha, beta)
      n = system%solvent_refractive_index
      alpha = system%solvent_alpha
      beta = system%solvent_beta
      ! Convert surface tension from N/m to mN/m (as used in SMD parametrization)
      msurft_mN = system%solvent_surface_tension_si*1000.0_wp
      ! unavailable from our system data -> default to 0
      arom = 0.0_wp
      fclbr = 0.0_wp

      self%alpha = alpha
      self%zk = 0.0_wp; self%zkk = 0.0_wp
      do Z = 1, max_elem
         self%zk(Z) = ref_zk(1, Z)*n + ref_zk(2, Z)*alpha + ref_zk(3, Z)*beta
         do Z2 = 1, max_elem
            self%zkk(Z, Z2) = ref_zkk(1, Z, Z2)*n + ref_zkk(2, Z, Z2)*alpha + ref_zkk(3, Z, Z2)*beta
         end do
      end do

      self%rzkk = ref_rzkk
      self%drzkk = ref_drzkk
      self%nc3 = ref_nc3
      self%rnc3 = ref_rnc3
      self%drnc3 = ref_drnc3
      self%s_m = ref_sg*msurft_mN + ref_sr2*(arom**2) + ref_sp2*(fclbr**2) + ref_sb2*(beta**2)
   end subroutine new_smd_model_parameters

   !> Fill default/reference tables (subset of numsa/src/smd/init.f90)
   subroutine init_default_tables(h2o, ref_zk_h2o, ref_zkk_h2o, ref_zk, ref_zkk, &
                                  ref_rzkk, ref_drzkk, ref_nc3, ref_rnc3, ref_drnc3, ref_sg, ref_sr2, ref_sp2, ref_sb2)
      logical, intent(in) :: h2o
      real(wp), intent(out) :: ref_zk_h2o(max_elem)
      real(wp), intent(out) :: ref_zkk_h2o(max_elem, max_elem)
      real(wp), intent(out) :: ref_zk(3, max_elem)
      real(wp), intent(out) :: ref_zkk(3, max_elem, max_elem)
      real(wp), intent(out) :: ref_rzkk(max_elem, max_elem)
      real(wp), intent(out) :: ref_drzkk(max_elem, max_elem)
      real(wp), intent(out) :: ref_nc3, ref_rnc3, ref_drnc3
      real(wp), intent(out) :: ref_sg, ref_sr2, ref_sp2, ref_sb2

      integer :: i
      ref_zk_h2o = 0.0_wp
      ref_zkk_h2o = 0.0_wp
      ref_zk = 0.0_wp
      ref_zkk = 0.0_wp
      ref_rzkk = 0.0_wp
      ref_drzkk = 0.0_wp

      if (h2o) then
         ! not used from here, we always init with h2o=.false. and fill both sets
      end if

      ! H2O special terms
      ref_zk_h2o(1) = 48.69_wp
      ref_zk_h2o(6) = 129.74_wp
      ref_zk_h2o(9) = 38.18_wp
      ref_zk_h2o(17) = 9.82_wp
      ref_zk_h2o(35) = -8.72_wp
      ref_zk_h2o(16) = -9.10_wp
      ref_zkk_h2o(1, 6) = -60.77_wp
      ref_zkk_h2o(6, 6) = -72.95_wp
      ref_zkk_h2o(8, 6) = 68.69_wp
      ref_zkk_h2o(7, 6) = -48.22_wp
      ref_zkk_h2o(8, 7) = 121.98_wp
      ref_zkk_h2o(8, 15) = 68.85_wp
      ref_nc3 = 84.10_wp

      ! Non-water linear combination coefficients
      ! n-dependence
      ref_zk(1, 6) = 58.10_wp
      ref_zk(1, 8) = -17.56_wp
      ref_zk(1, 7) = 32.62_wp
      ref_zk(1, 17) = -24.31_wp
      ref_zk(1, 35) = -35.42_wp
      ref_zk(1, 16) = -33.17_wp
      ref_zk(1, 14) = -18.04_wp
      ref_zkk(1, 1, 6) = -36.37_wp
      ref_zkk(1, 6, 6) = -62.05_wp
      ref_zkk(1, 1, 8) = -19.39_wp
      ref_zkk(1, 8, 6) = -15.70_wp
      ref_zkk(1, 6, 7) = -99.76_wp
      ! alpha-dependence
      ref_zk(2, 6) = 48.10_wp
      ref_zk(2, 8) = 193.06_wp
      ref_zkk(2, 8, 6) = 95.99_wp
      ref_zkk(2, 6, 7) = 152.20_wp
      ref_zkk(2, 7, 6) = -41.00_wp
      ! beta-dependence
      ref_zk(3, 6) = 32.87_wp
      ref_zk(3, 8) = -43.79_wp
      ref_zkk(3, 8, 8) = -128.16_wp
      ref_zkk(3, 8, 7) = 79.13_wp

      ! Switching radii and widths (subset used in sigma.f90 T)
      ref_rzkk(1, 6) = 1.55_wp
      ref_rzkk(1, 8) = 1.55_wp
      ref_rzkk(6, 1) = 1.55_wp
      ref_rzkk(6, 6) = 1.84_wp
      ref_rzkk(6, 7) = 1.84_wp
      ref_rzkk(6, 8) = 1.84_wp
      ref_rzkk(6, 9) = 1.84_wp
      ref_rzkk(6, 15) = 2.2_wp
      ref_rzkk(6, 16) = 2.2_wp
      ref_rzkk(6, 17) = 2.1_wp
      ref_rzkk(6, 35) = 2.3_wp
      ref_rzkk(6, 53) = 2.6_wp
      ref_rzkk(7, 6) = 1.84_wp
      ref_rzkk(8, 6) = 1.33_wp
      ref_rzkk(8, 7) = 1.5_wp
      ref_rzkk(8, 8) = 1.8_wp
      ref_rzkk(8, 15) = 2.1_wp
      ref_rnc3 = 1.225_wp

      ref_drzkk(1, 6) = 0.3_wp
      ref_drzkk(1, 8) = 0.3_wp
      ref_drzkk(6, 1) = 0.3_wp
      ref_drzkk(6, 6) = 0.3_wp
      ref_drzkk(6, 7) = 0.3_wp
      ref_drzkk(6, 8) = 0.3_wp
      ref_drzkk(6, 9) = 0.3_wp
      ref_drzkk(6, 15) = 0.3_wp
      ref_drzkk(6, 16) = 0.3_wp
      ref_drzkk(6, 17) = 0.3_wp
      ref_drzkk(6, 35) = 0.3_wp
      ref_drzkk(6, 53) = 0.3_wp
      ref_drzkk(7, 6) = 0.3_wp
      ref_drzkk(8, 6) = 0.1_wp
      ref_drzkk(8, 7) = 0.3_wp
      ref_drzkk(8, 8) = 0.3_wp
      ref_drzkk(8, 15) = 0.3_wp
      ref_drnc3 = 0.065_wp

      ref_sg = 0.35_wp
      ref_sr2 = -4.19_wp
      ref_sp2 = -6.68_wp
      ref_sb2 = 0.00_wp
   end subroutine init_default_tables

end module moist_model_component_cds_parameters
