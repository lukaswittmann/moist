module moist_radii_draco
   use mctc_env, only: error_type, fatal_error
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_io, only: structure_type
   use moist_radius_type, only: radius_type
   implicit none
   private

   public :: draco_radius_type
   public :: new_draco_radii

   !> Placeholder radii model for DRACO.
   type, extends(radius_type) :: draco_radius_type
   contains
      !> Placeholder update routine.
      procedure :: update => update_draco_radii
      !> Print DRACO model status.
      procedure :: print => print_draco_radii
   end type draco_radius_type

contains

   !> Constructor for the DRACO placeholder model.
   !> @param[out] self       DRACO radii model
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_draco_radii(self, verbosity)
      !> DRACO radii model
      type(draco_radius_type), intent(out) :: self
      !> Optional print level
      integer, intent(in), optional :: verbosity

      self%nat = 0
      if (present(verbosity)) then
         self%verbosity = verbosity
      else
         self%verbosity = 0
      end if
   end subroutine new_draco_radii

   !> Print DRACO radii model status.
   !> @param[in] self  DRACO radii model
   !> @param[in] unit  optional output unit
   subroutine print_draco_radii(self, unit)
      !> DRACO radii model
      class(draco_radius_type), intent(in) :: self
      !> Optional output unit
      integer, intent(in), optional :: unit

      integer :: iu

      iu = output_unit
      if (present(unit)) iu = unit

      write (iu, '(a)') "DRACO radii model: not implemented yet."
      write (iu, '(a,i0)') "Cached nat: ", self%nat
   end subroutine print_draco_radii

   !> Placeholder update implementation.
   !> @param[inout] self   DRACO radii model
   !> @param[in]    mol    molecular structure
   !> @param[out]   error  not-implemented error
   subroutine update_draco_radii(self, mol, error)
      !> DRACO radii model
      class(draco_radius_type), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      self%nat = mol%nat
      if (allocated(self%f0)) deallocate (self%f0)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      call fatal_error(error, "DRACO radii are not implemented yet")
   end subroutine update_draco_radii

end module moist_radii_draco
