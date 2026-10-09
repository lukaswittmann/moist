!> Resolved per-atom data of one side of an interaction
!>
!> - Resolved plain data for kernels: merged LJ parameters, tail charges and field-presence flag
!> - Per-atom data from owning terms
!> - Field passed separately from sites (3D only)
!> - Whole-side field term owns every charge; tail charges for Ng tails only
!> - No point-charge field beside a supplied field
module moist_model_moz_potential_sites
   use mctc_env, only: wp
   use moist_model_moz_potential_lj_base, only: lj_12_6_type
   implicit none(type, external)
   private

   public :: potential_sites_type

   !> Resolved per-atom data of one side
   type :: potential_sites_type
      !> Number of atoms of the side
      integer :: natom = 0
      !> Merged Lennard-Jones parameters; absent when no term has any
      type(lj_12_6_type), allocatable :: pair
      !> Merged tail charges, e (natom); absent when no term has charges
      real(wp), allocatable :: q(:)
      !> Whether a term of the side supplies a grid field
      logical :: has_field = .false.
   end type potential_sites_type

end module moist_model_moz_potential_sites
