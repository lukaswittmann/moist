!> Factory for the mutually exclusive RISM and DRISM bulk theories
module moist_model_moz_solvent_theory
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_moz_solvent_theory_type, only: bulk_theory_type, check_bulk_charge_neutrality
   use moist_model_moz_solvent_theory_rism, only: rism_theory_type
   use moist_model_moz_solvent_theory_drism, only: drism_theory_type, new_drism_theory
   implicit none(type, external)
   private
   public :: bulk_theory_type, rism_theory_type, new_bulk_theory, canonical_bulk_theory
contains
   !> Canonical name; XRISM is an alias for ordinary RISM
   !> @param[in] name Requested name
   pure function canonical_bulk_theory(name) result(canonical)
      character(len=*), intent(in) :: name
      character(len=8) :: canonical
      select case (to_lower(trim(adjustl(name))))
      case ("rism", "xrism")
         canonical = "RISM"
      case ("drism")
         canonical = "DRISM"
      case default
         canonical = ""
      end select
   end function canonical_bulk_theory

   !> Construct one fixed bulk theory from the solvent state and site charges
   !> @param[out] self Selected theory
   !> @param[in] name RISM (default) or DRISM
   !> @param[in] xyz Site coordinates, bohr (3, ns)
   !> @param[in] q Site charges, e (ns)
   !> @param[in] rho Molecular number density, 1/bohr^3
   !> @param[in] beta Inverse thermal energy, 1/Hartree
   !> @param[in] dielectric DRISM dielectric constant
   !> @param[in] smear DRISM Gaussian width, bohr
   !> @param[out] error Unknown theory, invalid bulk charges or invalid DRISM state
   subroutine new_bulk_theory(self, name, xyz, q, rho, beta, dielectric, smear, error)
      class(bulk_theory_type), allocatable, intent(out) :: self
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: xyz(:, :), q(:), rho, beta, dielectric, smear
      type(error_type), allocatable, intent(out) :: error
      select case (canonical_bulk_theory(name))
      case ("RISM")
         call check_bulk_charge_neutrality(q, error)
         if (.not. allocated(error)) allocate (rism_theory_type :: self)
      case ("DRISM")
         allocate (drism_theory_type :: self)
         select type (self)
         type is (drism_theory_type)
            call new_drism_theory(self, xyz, q, rho, beta, dielectric, smear, error)
         end select
      case default
         call fatal_error(error, "Unknown bulk theory: "//trim(name)//" (choose RISM or DRISM)")
      end select
      if (allocated(error)) then
         if (allocated(self)) deallocate (self)
      end if
   end subroutine new_bulk_theory
end module moist_model_moz_solvent_theory
