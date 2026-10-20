!> Ordinary RISM: the bulk intramolecular matrix is unchanged
module moist_model_moz_solvent_theory_rism
   use mctc_env, only: wp
   use moist_model_moz_solvent_theory_type, only: bulk_theory_type
   implicit none(type, external)
   private
   public :: rism_theory_type

   !> Default bulk theory with no dielectric correction
   type, extends(bulk_theory_type) :: rism_theory_type
   contains
      procedure :: name => rism_name
      procedure :: correction_at => rism_correction_at
   end type rism_theory_type
contains
   !> Canonical name
   !>
   !> @param[in] self Fixed theory
   pure function rism_name(self) result(name)
      !> Fixed theory
      class(rism_theory_type), intent(in) :: self
      !> Canonical name, "RISM"
      character(len=8) :: name
      name = "RISM"
   end function rism_name

   !> Ordinary RISM has zero dielectric correction at every wave number
   !>
   !> @param[in] self Fixed theory
   !> @param[in] k Wave number, 1/bohr
   !> @param[out] zeta Zero correction (ns, ns)
   !> @param[out] zeta_over_k2 Zero regularized correction (ns, ns)
   pure subroutine rism_correction_at(self, k, zeta, zeta_over_k2)
      !> Fixed theory
      class(rism_theory_type), intent(in) :: self
      !> Wave number, 1/bohr
      real(wp), intent(in) :: k
      !> Zero correction and zero regularized correction, both (ns, ns)
      real(wp), intent(out) :: zeta(:, :), zeta_over_k2(:, :)
      zeta = 0.0_wp
      zeta_over_k2 = 0.0_wp
   end subroutine rism_correction_at
end module moist_model_moz_solvent_theory_rism
