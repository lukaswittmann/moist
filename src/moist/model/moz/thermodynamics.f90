!> Thermodynamic quantities and corrections of a 3D UV solvation free energy
!>
!> - Partial molar volume Vbar = kT kappa_T (1 - D)
!>   from the k = 0 direct-correlation moment of the UV solution
!> - D = sum_a rho_a int c_a d^3r over physical c
!> - Additive correction terms, Hartree, from the bulk pressure P, the ideal-gas
!>   pressure P_ideal = rho kT and kappa_T of the solvent:
!>   PC = -P Vbar and PC+ = -(P - P_ideal) Vbar
!> - Reference: Sergiievskyi, Jeanmairet, Levesque, Borgis, JCP 143, 184116 (2015)
!> - Selectors none, PC and PC+, in any case
!> - Response: both corrections are linear in D
!> - dE/dD = (P or P - P_ideal) kappa_T/beta
!> - The response adds beta dE/dD to dPhi/dc
!> - A selected leading finite-size quantity has fixed charge, box and bulk inputs
!> - Its energy is included in the same list and gives no host-response seed
module moist_model_moz_thermodynamics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_moz_solvent_vv_thermodynamics, only: solvent_vv_thermodynamics_type
   use moist_model_moz_thermodynamics_kernel, only: pc_moment_shift, pc_plus_moment_shift, moment_response_shift
   use moist_model_moz_thermodynamics_type, only: thermodynamic_quantity_type, thermo_input_type
   use moist_model_moz_thermodynamics_list, only: thermo_list_type
   use moist_model_moz_thermodynamics_quantities, only: new_uv_thermo_list, new_bulk_thermo_list, &
      & huenenberger_correction_type, new_huenenberger_correction
   implicit none(type, external)
   private

   public :: model_moz_3d_thermodynamics_type, new_moz_3d_thermodynamics, check_correction_name
   public :: correction_moment_shift
   public :: thermodynamic_quantity_type, thermo_input_type, thermo_list_type, new_uv_thermo_list, new_bulk_thermo_list
   public :: accumulate_solvation_adjoint
   public :: huenenberger_correction_type, new_huenenberger_correction, check_finite_size_name

   !> Solvation free energy of a 3D UV solution and its correction terms, atomic units
   type :: model_moz_3d_thermodynamics_type
      !> Uncorrected solvation free energy, Hartree
      real(wp) :: free_energy = 0.0_wp
      !> Partial molar volume, bohr^3
      real(wp) :: partial_molar_volume = 0.0_wp
      !> Correction terms -P Vbar (PC) and -(P - P_ideal) Vbar (PC+), Hartree
      real(wp) :: pc = 0.0_wp, pc_plus = 0.0_wp
      !> Selected electrostatic finite-size correction, Hartree; fixed in the host response
      real(wp) :: finite_size = 0.0_wp
   contains
      !> Free energy with the selected correction term added
      procedure :: corrected_energy
   end type model_moz_3d_thermodynamics_type

contains

   !> Partial molar volume and correction terms of a 3D UV solution
   !>
   !> @param[out] self         thermodynamics
   !> @param[in]  free_energy  uncorrected solvation free energy, Hartree
   !> @param[in]  moment       D = sum_a rho_a int c_a d^3r of the 3D UV solution
   !> @param[in]  solvent      bulk thermodynamics of the solvent's VV solution
   !> @param[in]  beta         inverse thermal energy, 1/Hartree
   !> @param[out] error        non-finite input or non-positive beta
   !> @param[out] items        optional owned forward cache of the evaluated list
   !> @param[in]  correction   optional single selected correction, none, PC or PC+
   !> @param[in]  finite_size  optional finite-size quantity, registered and selected with the pressure term
   subroutine new_moz_3d_thermodynamics(self, free_energy, moment, solvent, beta, error, items, correction, finite_size)
      !> Thermodynamics
      type(model_moz_3d_thermodynamics_type), intent(out) :: self
      !> Uncorrected free energy and direct-correlation moment
      real(wp), intent(in) :: free_energy, moment
      !> Bulk thermodynamics of the solvent
      type(solvent_vv_thermodynamics_type), intent(in) :: solvent
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional owned forward cache of the evaluated list
      type(thermo_list_type), optional, intent(out), target :: items
      !> Optional single selected correction: none, PC or PC+
      character(len=*), optional, intent(in) :: correction
      !> Optional finite-size quantity, registered and selected alongside the pressure term
      class(thermodynamic_quantity_type), optional, intent(in) :: finite_size
      !> Owned cache when no output list is requested
      type(thermo_list_type), target :: local_items
      !> Build in the output directly; ifort loses cache state in the optional-output copy
      type(thermo_list_type), pointer :: evaluated

      if (present(items)) then
         evaluated => items
      else
         evaluated => local_items
      end if

      if (.not. (ieee_is_finite(free_energy) .and. ieee_is_finite(moment) .and. &
         & ieee_is_finite(solvent%pressure) .and. ieee_is_finite(solvent%pressure_ideal) .and. &
         & ieee_is_finite(solvent%compressibility) .and. ieee_is_finite(beta))) then
         call fatal_error(error, "3D UV thermodynamics need finite free energy, moment and solvent thermodynamics")
         return
      end if
      if (beta <= 0.0_wp) then
         call fatal_error(error, "3D UV thermodynamics need a positive inverse thermal energy")
         return
      end if
      self%free_energy = free_energy
      if (present(correction)) then
         call new_uv_thermo_list(evaluated, free_energy, moment, solvent%pressure, solvent%pressure_ideal, &
            & solvent%compressibility, beta, error, corrections=[correction])
      else
         call new_uv_thermo_list(evaluated, free_energy, moment, solvent%pressure, solvent%pressure_ideal, &
            & solvent%compressibility, beta, error)
      end if
      if (allocated(error)) return
      if (present(finite_size)) then
         call evaluated%add(finite_size, error, selected=.true.)
         if (allocated(error)) return
      end if
      call evaluated%evaluate(error)
      if (.not. allocated(error)) call evaluated%get_value("partial_molar_volume", self%partial_molar_volume, error)
      if (.not. allocated(error)) call evaluated%get_value("pc", self%pc, error)
      if (.not. allocated(error)) call evaluated%get_value("pc+", self%pc_plus, error)
      if (allocated(error)) return
      if (present(finite_size)) then
         call evaluated%get_value(finite_size%name, self%finite_size, error)
         if (allocated(error)) return
      end if
   end subroutine new_moz_3d_thermodynamics

   !> Accumulate ordinary energy derivatives, then convert dE/dD once for the UV driver
   !>
   !> - Bulk pressure/compressibility leaves are fixed in the current host response
   !>
   !> @param[inout] items        evaluated list from new_moz_3d_thermodynamics
   !> @param[out]   energy_seed  adjoint of the free energy
   !> @param[out]   shift        UV response shift beta dE/dD, dimensionless
   !> @param[out]   error        missing or non-positive beta, or failed accumulation
   subroutine accumulate_solvation_adjoint(items, energy_seed, shift, error)
      !> Evaluated list, accumulated in place
      type(thermo_list_type), intent(inout) :: items
      !> Adjoint of the free energy; UV response shift beta dE/dD, dimensionless
      real(wp), intent(out) :: energy_seed, shift
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: derivative, beta
      energy_seed = 0.0_wp
      shift = 0.0_wp
      call items%get_value("beta", beta, error)
      if (allocated(error)) return
      if (.not. ieee_is_finite(beta) .or. beta <= 0.0_wp) then
         call fatal_error(error, "Thermodynamic response needs a positive inverse thermal energy")
         return
      end if
      call items%accumulate_adjoint(error)
      if (.not. allocated(error)) call items%get_adjoint("free_energy", energy_seed, error)
      if (.not. allocated(error)) call items%get_adjoint("uv_moment", derivative, error)
      if (allocated(error)) return
      shift = moment_response_shift(beta, derivative)
   end subroutine accumulate_solvation_adjoint

   !> Free energy with the selected pressure term and configured finite-size term, Hartree
   !>
   !> @param[in]    self        thermodynamics
   !> @param[in]    correction  none, PC or PC+, case-insensitive
   !> @param[inout] energy      corrected free energy; unchanged on error
   !> @param[out]   error       unknown correction selector
   subroutine corrected_energy(self, correction, energy, error)
      !> Thermodynamics
      class(model_moz_3d_thermodynamics_type), intent(in) :: self
      !> Selector
      character(len=*), intent(in) :: correction
      !> Free energy
      real(wp), intent(inout) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_correction_name(correction, error)
      if (allocated(error)) return

      energy = self%free_energy + self%finite_size
      select case (to_lower(trim(adjustl(correction))))
      case ("pc")
         energy = energy + self%pc
      case ("pc+")
         energy = energy + self%pc_plus
      case default
         ! check_correction_name above leaves only "none" here
         continue
      end select
   end subroutine corrected_energy

   !> Moment term of a correction for the response: beta dE/dD, dimensionless
   !>
   !> P kappa_T for PC, (P - P_ideal) kappa_T for PC+ and zero for none
   !>
   !> @param[in] correction  none, PC or PC+; any other name gives zero
   !> @param[in] solvent     bulk thermodynamics of the solvent's VV solution
   pure function correction_moment_shift(correction, solvent) result(shift)
      !> Selector
      character(len=*), intent(in) :: correction
      !> Bulk thermodynamics of the solvent
      type(solvent_vv_thermodynamics_type), intent(in) :: solvent
      !> Moment term
      real(wp) :: shift

      shift = 0.0_wp
      select case (to_lower(trim(adjustl(correction))))
      case ("pc")
         shift = pc_moment_shift(solvent%pressure, solvent%compressibility)
      case ("pc+")
         shift = pc_plus_moment_shift(solvent%pressure, solvent%pressure_ideal, solvent%compressibility)
      case default
         continue
      end select
   end function correction_moment_shift

   !> Refuse an unknown correction selector
   !>
   !> @param[in]  correction  selector
   !> @param[out] error       name other than none, PC or PC+
   subroutine check_correction_name(correction, error)
      !> Selector
      character(len=*), intent(in) :: correction
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select case (to_lower(trim(adjustl(correction))))
      case ("none", "pc", "pc+")
         continue
      case default
         call fatal_error(error, "Unknown solvation free-energy correction: "//trim(correction)// &
            & " (choose none, PC or PC+)")
      end select
   end subroutine check_correction_name

   !> Refuse an unknown finite-size selector
   !>
   !> - The finite-size term has a separate selector and may accompany PC or PC+
   !>
   !> @param[in]  name   selector, case-insensitive
   !> @param[out] error  name other than none or huenenberger
   subroutine check_finite_size_name(name, error)
      !> Selector
      character(len=*), intent(in) :: name
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      select case (to_lower(trim(adjustl(name))))
      case ("none", "huenenberger")
         continue
      case default
         call fatal_error(error, "Unknown electrostatic finite-size correction: "//trim(name)// &
            & " (choose none or huenenberger)")
      end select
   end subroutine check_finite_size_name

end module moist_model_moz_thermodynamics
