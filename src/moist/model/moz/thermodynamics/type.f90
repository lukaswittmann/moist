!> One scalar thermodynamic quantity or energy correction, in atomic units
module moist_model_moz_thermodynamics_type
   use mctc_env, only: wp, error_type
   implicit none(type, external)
   private
   public :: thermodynamic_quantity_type, thermo_input_type

   !> Abstract scalar item with named scalar dependencies
   !>
   !> - Base of diagnostics, inputs and energy corrections
   type, abstract :: thermodynamic_quantity_type
      !> Item name, unique within a list, and the atomic-unit label of its scalar value
      character(len=:), allocatable :: name, unit
      !> Dependency item names (n_dep), in the order expected by evaluate and accumulate_adjoint
      character(len=:), allocatable :: dependencies(:)
      !> Energy correction, selectable as an additive energy term; default false
      logical :: correction = .false.
      !> Group label; a nonempty group permits at most one selected member
      character(len=:), allocatable :: exclusive_group
   contains
      procedure(evaluate_item), deferred :: evaluate
      procedure(accumulate_item), deferred :: accumulate_adjoint
   end type thermodynamic_quantity_type

   !> Scalar supplied by a grid reduction or the bulk solvent
   type, extends(thermodynamic_quantity_type) :: thermo_input_type
      !> Input value, in the declared unit; default 0
      real(wp) :: scalar = 0.0_wp
   contains
      procedure :: evaluate => input_evaluate
      procedure :: accumulate_adjoint => input_accumulate
   end type thermo_input_type

   !> Deferred item methods
   abstract interface
      !> Evaluate the item from its dependency values
      !>
      !> - Dependency values arrive in the item's declared order
      !> - value has the item's unit
      !>
      !> @param[inout] self    item to evaluate
      !> @param[in]    values  dependency values (n_dep), declared order
      !> @param[out]   value   item value, in the item unit
      !> @param[out]   error   error handler, allocated on failure
      subroutine evaluate_item(self, values, value, error)
         import thermodynamic_quantity_type, wp, error_type
         implicit none(type, external)
         !> Item to evaluate
         class(thermodynamic_quantity_type), intent(inout) :: self
         !> Dependency values (n_dep), declared order
         real(wp), intent(in) :: values(:)
         !> Item value, in the item unit
         real(wp), intent(out) :: value
         !> Error handler, allocated on failure
         type(error_type), allocatable, intent(out) :: error
      end subroutine evaluate_item

      !> Add reverse sensitivities of the item to its dependencies
      !>
      !> - Adds seed * d(value)/d(dependency) to every dependency sensitivity
      !> - seed = dE/dvalue with E in Hartree, so Hartree per item unit
      !> - Dependency adjoints are Hartree per dependency unit; the UV driver applies beta scaling
      !>
      !> @param[in]    self      item
      !> @param[in]    values    dependency values (n_dep), declared order
      !> @param[in]    value     forward value of the item, in the item unit
      !> @param[in]    seed      dE/dvalue, Hartree per item unit
      !> @param[inout] adjoints  dependency sensitivities (n_dep), accumulated in place, Hartree per dependency unit
      pure subroutine accumulate_item(self, values, value, seed, adjoints)
         import thermodynamic_quantity_type, wp
         implicit none(type, external)
         !> Item
         class(thermodynamic_quantity_type), intent(in) :: self
         !> Dependency values (n_dep), item value and seed dE/dvalue [Hartree per item unit]
         real(wp), intent(in) :: values(:), value, seed
         !> Dependency sensitivities (n_dep), accumulated in place, Hartree per dependency unit
         real(wp), intent(inout) :: adjoints(:)
      end subroutine accumulate_item
   end interface

contains

   !> Return the stored scalar
   !>
   !> @param[inout] self    scalar input
   !> @param[in]    values  dependency values, empty
   !> @param[out]   value   stored scalar, in the item unit
   !> @param[out]   error   error handler, never set
   subroutine input_evaluate(self, values, value, error)
      !> Scalar input
      class(thermo_input_type), intent(inout) :: self
      !> Dependency values, empty
      real(wp), intent(in) :: values(:)
      !> Stored scalar, in the item unit
      real(wp), intent(out) :: value
      !> Error handler, never set
      type(error_type), allocatable, intent(out) :: error
      value = self%scalar
   end subroutine input_evaluate

   !> No-op; the list keeps the leaf sensitivity
   !>
   !> @param[in]    self      scalar input
   !> @param[in]    values    dependency values, empty, unused
   !> @param[in]    value     forward value, unused
   !> @param[in]    seed      dE/dvalue, unused
   !> @param[inout] adjoints  dependency sensitivities, empty, unused
   pure subroutine input_accumulate(self, values, value, seed, adjoints)
      !> Scalar input
      class(thermo_input_type), intent(in) :: self
      !> Dependency values, forward value and seed, all unused
      real(wp), intent(in) :: values(:), value, seed
      !> Dependency sensitivities, empty, unused
      real(wp), intent(inout) :: adjoints(:)
      ! The list retains the accumulated leaf sensitivity; there are no dependencies
   end subroutine input_accumulate

end module moist_model_moz_thermodynamics_type
