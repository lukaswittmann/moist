!> Property request flags for the DROP cavity
!>
!> Controls which optional properties are computed and stored during cavity construction
!> Core quantities (area, volume, Gaussian widths, ...) are always computed
module moist_cavity_drop_request
   use, intrinsic :: iso_fortran_env, only: output_unit
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter

   implicit none

   public :: drop_property_request
   public :: drop_request_default
   public :: drop_request_diagnostics
   public :: drop_request_fine
   private

   !> Property request flags for DROP cavity computation
   !>
   !> Flags control optional quantities beyond the core pipeline
   !> (area, volume, Gaussian widths, gradient). All flags default
   !> to `.false.` for minimal overhead.
   type :: drop_property_request

      !> Compute local grid-point density (diagnostic)
      logical :: grid_point_density = .false.

      !> Compute principal, mean, and Gaussian curvatures (diagnostic)
      logical :: curvature = .false.

      !> Store surface normal vectors at grid points
      logical :: normal = .false.

      !> Store sphere-center to grid-point distances
      logical :: r_iI = .false.

      !> Store anchor-to-projected-point displacement distances
      logical :: rho = .false.

      !> Compute CPCM solvation energy
      logical :: cpcm = .false.

      !> Compute marching-cubes reference area and volume
      logical :: mc = .false.

   contains
      !> Print active flags to output
      procedure :: print => print_property_request

   end type drop_property_request

contains

   !> Print active property flags to standard output
   !>
   !> @param[in] self  Request to display
   subroutine print_property_request(self)
      class(drop_property_request), intent(in) :: self
      type(prettyprinter) :: pp

      pp = new_prettyprinter(unit=output_unit)

      call pp%blank()
      call pp%push('Property Request:')
      call pp%kv('Grid point density', self%grid_point_density)
      call pp%kv('Curvature', self%curvature)
      call pp%kv('CPCM', self%cpcm)
      call pp%kv('Surface normals', self%normal)
      call pp%kv('r_iI distances', self%r_iI)
      call pp%kv('rho displacements', self%rho)
      call pp%kv('Marching cubes', self%mc)
      call pp%pop()
      call pp%blank()

   end subroutine print_property_request

   !> Default request: no optional properties
   !>
   !> @return  Request with all flags off
   type(drop_property_request) function drop_request_default() result(req)
   end function drop_request_default

   !> Diagnostic request: curvature, grid density, normals, sphere-center
   !> distances, and anchor displacements
   !>
   !> @return  Request with diagnostic properties enabled
   type(drop_property_request) function drop_request_diagnostics() result(req)
      req%grid_point_density = .true.
      req%curvature = .true.
      req%normal = .true.
      req%r_iI = .true.
      req%rho = .true.
   end function drop_request_diagnostics

   !> Full request: all geometric diagnostic properties enabled
   !> (CPCM energy and marching-cubes reference remain off; opt in separately)
   !>
   !> @return  Request with all geometric diagnostic properties enabled
   type(drop_property_request) function drop_request_fine() result(req)
      req%grid_point_density = .true.
      req%curvature = .true.
      req%normal = .true.
      req%r_iI = .true.
      req%rho = .true.
   end function drop_request_fine

end module moist_cavity_drop_request
