!> Free energy and partial molar volume of a 1D UV solution, in atomic units
module moist_model_moz_1d_thermodynamics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_model_moz_thermodynamics_kernel, only: partial_molar_volume
   implicit none(type, external)
   private

   public :: model_moz_1d_thermodynamics_type, new_moz_1d_thermodynamics

   type :: model_moz_1d_thermodynamics_type
      !> Solvation free energy, Hartree
      real(wp) :: free_energy = 0.0_wp
      !> Partial molar volume, bohr^3
      real(wp) :: partial_molar_volume = 0.0_wp
   end type model_moz_1d_thermodynamics_type

contains

   !> Construct the diagnostics from a converged UV solution and bulk compressibility
   !>
   !> @param[out] self         Thermodynamics
   !> @param[in]  free_energy  Solvation free energy, Hartree
   !> @param[in]  moment       D = sum_ia rho_a int c_ia d^3r, dimensionless
   !> @param[in]  kappa        Bulk isothermal compressibility, bohr^3/Hartree
   !> @param[in]  beta         Inverse thermal energy, Hartree^-1
   !> @param[out] error        Non-finite input or volume, or non-positive beta
   subroutine new_moz_1d_thermodynamics(self, free_energy, moment, kappa, beta, error)
      type(model_moz_1d_thermodynamics_type), intent(out) :: self
      real(wp), intent(in) :: free_energy, moment, kappa, beta
      type(error_type), allocatable, intent(out) :: error

      if (.not. all(ieee_is_finite([free_energy, moment, kappa, beta]))) then
         call fatal_error(error, "1D UV thermodynamics needs finite free energy, moment and bulk compressibility")
         return
      end if
      if (beta <= 0.0_wp) then
         call fatal_error(error, "1D UV thermodynamics needs a positive inverse thermal energy")
         return
      end if
      self%free_energy = free_energy
      self%partial_molar_volume = partial_molar_volume(moment, kappa, beta)
      if (.not. ieee_is_finite(self%partial_molar_volume)) then
         call fatal_error(error, "1D UV partial molar volume is non-finite")
      end if
   end subroutine new_moz_1d_thermodynamics

end module moist_model_moz_1d_thermodynamics
