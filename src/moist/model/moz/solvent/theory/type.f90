!> Fixed bulk-theory correction shared by the VV solve and solvent response
module moist_model_moz_solvent_theory_type
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   implicit none(type, external)
   private
   public :: bulk_theory_type
   public :: bulk_charges_neutral, check_bulk_charge_neutrality

   !> One theory for one rigid bulk solvent at a fixed state point
   type, abstract :: bulk_theory_type
   contains
      !> Canonical name of the theory
      procedure(theory_name_interface), deferred :: name
      !> Dielectric correction and its analytical quotient by k squared
      procedure(correction_interface), deferred :: correction_at
   end type bulk_theory_type

   !> Deferred interfaces of a bulk theory
   abstract interface
      !> Name of the selected theory
      !>
      !> @param[in] self Fixed theory
      pure function theory_name_interface(self) result(name)
         import :: bulk_theory_type
         implicit none(type, external)
         !> Fixed theory
         class(bulk_theory_type), intent(in) :: self
         !> Canonical name
         character(len=8) :: name
      end function theory_name_interface

      !> Evaluate zeta(k) and zeta(k)/k^2, including the finite zero limit
      !>
      !> @param[in] self Fixed theory
      !> @param[in] k Nonnegative wave number, 1/bohr
      !> @param[out] zeta Symmetric correction, bohr^3 (ns, ns)
      !> @param[out] zeta_over_k2 Regularized correction, bohr^5 (ns, ns)
      pure subroutine correction_interface(self, k, zeta, zeta_over_k2)
         import :: bulk_theory_type, wp
         implicit none(type, external)
         !> Fixed theory
         class(bulk_theory_type), intent(in) :: self
         !> Nonnegative wave number, 1/bohr
         real(wp), intent(in) :: k
         !> Symmetric correction (bohr^3) and regularized correction (bohr^5), both (ns, ns)
         real(wp), intent(out) :: zeta(:, :), zeta_over_k2(:, :)
      end subroutine correction_interface
   end interface

contains

   !> Neutral molecular bulk charges within charge-table precision and summation error
   !>
   !> @param[in] q  site charges, e (ns)
   pure function bulk_charges_neutral(q) result(neutral)
      !> Site charges, e (ns)
      real(wp), intent(in) :: q(:)
      !> Neutral within tolerance
      logical :: neutral
      real(wp), parameter :: charge_floor = 1.0e-6_wp
      real(wp) :: total, magnitude, roundoff

      neutral = .false.
      if (size(q) < 1) return
      if (.not. all(ieee_is_finite(q))) return
      total = sum(q)
      magnitude = sum(abs(q))
      if (.not. ieee_is_finite(total) .or. .not. ieee_is_finite(magnitude)) return
      roundoff = 64.0_wp*epsilon(1.0_wp)*max(1.0_wp, magnitude)
      ! Summation error also applies to a table residual at the decimal boundary
      neutral = abs(total) <= charge_floor + roundoff
   end function bulk_charges_neutral

   !> Reject unsupported bulk charges without changing the supplied charge model
   !>
   !> @param[in]  q      site charges, e (ns)
   !> @param[out] error  unsupported charges: empty, nonfinite or net charge
   subroutine check_bulk_charge_neutrality(q, error)
      !> Site charges, e (ns)
      real(wp), intent(in) :: q(:)
      !> Unsupported charge failure, unallocated if neutral
      type(error_type), allocatable, intent(out) :: error
      character(len=24) :: charge

      if (bulk_charges_neutral(q)) return
      if (size(q) < 1) then
         call fatal_error(error, "Bulk solvent requires at least one site charge")
      else if (.not. all(ieee_is_finite(q))) then
         call fatal_error(error, "Bulk solvent site charges must be finite")
      else if (.not. ieee_is_finite(sum(q)) .or. .not. ieee_is_finite(sum(abs(q)))) then
         call fatal_error(error, "Bulk solvent charge totals must be finite")
      else
         write (charge, "(es24.16)") sum(q)
         call fatal_error(error, "Bulk solvent must be electrically neutral; net charge = "// &
            & trim(adjustl(charge))//" e (tolerance: 1e-6 e plus floating-point summation error)")
      end if
   end subroutine check_bulk_charge_neutrality
end module moist_model_moz_solvent_theory_type
