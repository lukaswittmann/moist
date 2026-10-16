!> One scalar thermodynamic quantity or energy correction, in atomic units
module moist_model_moz_thermodynamics_type
   use mctc_env, only: wp, error_type
   implicit none(type, external)
   private
   public :: thermodynamic_quantity_type, thermo_input_type

   !> One base for diagnostics and corrections; dependencies are named scalar items
   type, abstract :: thermodynamic_quantity_type
      !> Name unique within a list, and the atomic-unit label of the scalar value
      character(len=:), allocatable :: name, unit
      !> Named inputs, in the order expected by evaluate and accumulate_adjoint
      character(len=:), allocatable :: dependencies(:)
      !> Only corrections can be selected as additive energy contributions
      logical :: correction = .false.
      !> Nonempty groups permit at most one selected member
      character(len=:), allocatable :: exclusive_group
   contains
      procedure(evaluate_item), deferred :: evaluate
      procedure(accumulate_item), deferred :: accumulate_adjoint
   end type thermodynamic_quantity_type

   !> A scalar already supplied by a grid reduction or the bulk solvent
   type, extends(thermodynamic_quantity_type) :: thermo_input_type
      real(wp) :: scalar = 0.0_wp
   contains
      procedure :: evaluate => input_evaluate
      procedure :: accumulate_adjoint => input_accumulate
   end type thermo_input_type

   abstract interface
      !> Dependencies are in the item's declared order; value has the item's unit
      subroutine evaluate_item(self, values, value, error)
         import thermodynamic_quantity_type, wp, error_type
         implicit none(type, external)
         class(thermodynamic_quantity_type), intent(inout) :: self
         real(wp), intent(in) :: values(:)
         real(wp), intent(out) :: value
         type(error_type), allocatable, intent(out) :: error
      end subroutine evaluate_item

      !> Add seed * d(value)/d(dependency) to every dependency sensitivity
      !> seed = dE/dvalue, with E in Hartree: units are Hartree per item unit
      !> Dependency adjoints are Hartree per dependency unit; the UV driver applies beta scaling
      pure subroutine accumulate_item(self, values, value, seed, adjoints)
         import thermodynamic_quantity_type, wp
         implicit none(type, external)
         class(thermodynamic_quantity_type), intent(in) :: self
         real(wp), intent(in) :: values(:), value, seed
         real(wp), intent(inout) :: adjoints(:)
      end subroutine accumulate_item
   end interface

contains

   subroutine input_evaluate(self, values, value, error)
      class(thermo_input_type), intent(inout) :: self
      real(wp), intent(in) :: values(:)
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out) :: error
      value = self%scalar
   end subroutine input_evaluate

   pure subroutine input_accumulate(self, values, value, seed, adjoints)
      class(thermo_input_type), intent(in) :: self
      real(wp), intent(in) :: values(:), value, seed
      real(wp), intent(inout) :: adjoints(:)
      ! The list retains the accumulated leaf sensitivity; there are no dependencies.
   end subroutine input_accumulate

end module moist_model_moz_thermodynamics_type
