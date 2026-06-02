module moist_radius_type
   use mctc_env, only: wp
   use mctc_env, only: error_type
   use mctc_io, only: structure_type
   implicit none
   private

   public :: radius_type

   ! TODO: Add error propagation

   !> Abstract base class for cached atomic radii models.
   type, abstract :: radius_type
      !> Number of atoms in the cached state.
      integer :: nat = 0
      !> Verbosity level controlling model diagnostics.
      integer :: verbosity = 0
      !> Cached radii vector in atomic units.
      real(wp), allocatable :: f0(:)
      !> Cached derivatives dR_i/dR_A with shape (3,nat,nat).
      real(wp), allocatable :: f1_rA(:, :, :)
   contains
      !> Update internal radii state for a molecular structure.
      procedure(radii_update_interface), deferred :: update
      !> Print model-specific radii information.
      procedure(radii_print_interface), deferred :: print
   end type radius_type

   abstract interface

      !> Update model state and cache radii-dependent data.
      !> @param[inout] self   radii model instance
      !> @param[in]    mol    molecular structure
      !> @param[out]   error  error handle on failure
      subroutine radii_update_interface(self, mol, error)
         import :: radius_type, structure_type, error_type
         class(radius_type), intent(inout) :: self
         type(structure_type), intent(in) :: mol
         type(error_type), allocatable, intent(out) :: error
      end subroutine radii_update_interface

      !> Print model-specific radii information.
      !> @param[in] self  radii model instance
      !> @param[in] unit  optional Fortran output unit
      subroutine radii_print_interface(self, unit)
         import :: radius_type
         class(radius_type), intent(in) :: self
         integer, intent(in), optional :: unit
      end subroutine radii_print_interface

   end interface

end module moist_radius_type
