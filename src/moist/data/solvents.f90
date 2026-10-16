
!> Solvent tables: descriptors, geometries, atomic charges and multipoles
!>
!> The data lives in the includes under data/solvents; the solvent objects
!> built from it are in moist_solvent
module moist_data_solvents
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use mctc_io_utils, only: to_lower
   implicit none(type, external)

   !> Number of entries in the solvent descriptor table
   integer, parameter, public :: max_solvents = 180

   !> Longest canonical solvent name, MOIST_NAME_MAX in the C API
   integer, parameter, public :: max_solvent_name_len = 32

   public :: get_solvent_id, get_solvent_table, get_solvent_geometry
   public :: solvent_multipole_data_type, get_solvent_charges, get_solvent_multipoles

   private

   !> Atom-centered Cartesian multipoles in e*bohr^rank
   type :: solvent_multipole_data_type
      !> Charges, (nat)
      real(wp), allocatable :: monopole(:)
      !> Dipoles X Y Z, (3, nat)
      real(wp), allocatable :: dipole(:, :)
      !> Quadrupoles XX YY ZZ XY XZ YZ, (6, nat)
      real(wp), allocatable :: quadrupole(:, :)
      !> Octupoles XXX YYY ZZZ XXY XXZ XYY XYZ XZZ YYZ YZZ, (10, nat)
      real(wp), allocatable :: octupole(:, :)
   end type solvent_multipole_data_type

contains

   !> Get the solvent ID from an alias
   !>
   !> @param[in]  alias      solvent name or alias; case and surrounding blanks ignored
   !> @param[out] solvent_id table ID, 0 on error
   !> @param[out] error      blank or unknown alias
   subroutine get_solvent_id(alias, solvent_id, error)

      !> Solvent alias (case insensitive)
      character(len=*), intent(in) :: alias

      !> Solvent ID
      integer, intent(out) :: solvent_id

      !> Error output for unknown aliases
      type(error_type), allocatable, intent(out) :: error

      !> Iterables
      integer :: i, j

      character(len=64) :: name_list(max_solvents)
      character(len=64) :: alias_list(10, max_solvents)
      character(:), allocatable :: query

      integer, dimension(max_solvents) :: id_list
      real(wp), dimension(max_solvents) :: eps, refr, A, B, g, rho

      !> Get basic solvent information
      call get_solvent_table(id_list, name_list, alias_list, eps, refr, A, B, g, rho)

      solvent_id = 0

      ! Reject blank queries, which would otherwise match the blank alias padding
      if (len_trim(alias) == 0) then
         call fatal_error(error, message="Empty solvent alias", stat=1)
         return
      end if

      ! Normalise
      query = trim(adjustl(to_lower(alias)))

      ! Search for the alias in the alias list, normalised like the query
      do i = 1, max_solvents
         do j = 1, 10
            if (len_trim(alias_list(j, i)) == 0) cycle
            if (query == trim(adjustl(to_lower(alias_list(j, i))))) then
               solvent_id = id_list(i)
               return
            end if
         end do
      end do

      ! If we reach here, the alias was not found
      call fatal_error(error, message="Unknown solvent: "//trim(alias), stat=1)

   end subroutine get_solvent_id

   !> Solvent descriptor table; entry i belongs to table ID ids(i)
   !>
   !> @param[out] ids              table IDs
   !> @param[out] names            canonical names
   !> @param[out] aliases          aliases per entry, blank padded
   !> @param[out] permittivity     static relative permittivity
   !> @param[out] refractive_index refractive index at the sodium D line
   !> @param[out] abraham_alpha    Abraham hydrogen-bond acidity
   !> @param[out] abraham_beta     Abraham hydrogen-bond basicity
   !> @param[out] surface_tension  surface tension in mN/m
   !> @param[out] mass_density     mass density in kg/m^3
   include "solvents/descriptors.inc"

   !> Solvent geometry in Angstrom
   include "solvents/geometries.inc"

   !> Atomic charges in e
   include "solvents/charges.inc"

   !> Signed, atom-centered Cartesian moments in e*bohr^rank
   include "solvents/multipoles.inc"

end module moist_data_solvents
