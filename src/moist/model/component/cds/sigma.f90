!> SMD per-atom surface tension coefficients (sigma_k)
!>
!> This module computes the atom-resolved SMD surface tension values
!> based on local environment (neighbors) using the switching function
!> and solvent-dependent parameters as in the SMD model. These sigma_k
!> values, together with the molecular SASA (from the iSwiG surface),
!> form the CDS contribution to the solvation energy.
module moist_model_component_cds_sigma
   use mctc_env, only: wp
   use mctc_io_symbols, only: to_number
   use mctc_io_convert, only: autoaa
   use moist_model_component_cds_parameters, only: smd_param
   implicit none
   private
   public :: cds_surfacetensions, calc_surfacetensions

   type :: cds_info
      !> Number of atoms in the solute
      integer :: nat
      !> Atomic numbers per atom (Z)
      integer, allocatable :: Z(:)
   end type cds_info

   type :: cds_surfacetensions
      !> Per-atom surface tension coefficient sigma_k (energy per area)
      real(wp), allocatable :: sk(:)
      !> Solvent-dependent macroscopic surface tension contribution
      real(wp) :: sm
   end type cds_surfacetensions

contains

   !> Compute SMD per-atom surface tension values
   subroutine calc_surfacetensions(coords, species, ident, param, surft)
      !> Atomic coordinates in bohr [3,nat]
      real(wp), intent(in) :: coords(:, :)
      integer, intent(in) :: species(:)
      !> Element symbols table (lookup via identifiers of `species`)
      character(len=*), intent(in) :: ident(:)
      type(smd_param), intent(in) :: param
      type(cds_surfacetensions), intent(out) :: surft

      type(cds_info) :: info
      integer :: i, j, k
      real(wp) :: s_temp1, s_temp2, s_temp3, s_temp4
      real(wp) :: nc_temp, add_temp

      call init_info(species, ident, info)
      allocate (surft%sk(info%nat))

      s_temp1 = 0.0_wp; s_temp2 = 0.0_wp; s_temp3 = 0.0_wp; s_temp4 = 0.0_wp
      nc_temp = 0.0_wp; add_temp = 0.0_wp
      surft%sm = param%s_m

      do i = 1, info%nat
         select case (info%Z(i))
         case (1) ! H
            do j = 1, info%nat
               select case (info%Z(j))
               case (6)  ! H,C
                  s_temp1 = s_temp1 + switch_T(coords(:, i), coords(:, j), param%rzkk(1, 6), param%drzkk(1, 6))
               case (8)  ! H,O
                  s_temp2 = s_temp2 + switch_T(coords(:, i), coords(:, j), param%rzkk(1, 8), param%drzkk(1, 8))
               case default
               end select
            end do
            surft%sk(i) = param%zk(1) &
                          + param%zkk(1, 6)*s_temp1 &
                          + param%zkk(1, 8)*s_temp2
            s_temp1 = 0.0_wp; s_temp2 = 0.0_wp
         case (6) ! C
            do j = 1, info%nat
               select case (info%Z(j))
               case (6)  ! C,C
                  if (i /= j) s_temp1 = s_temp1 + switch_T(coords(:, i), coords(:, j), param%rzkk(6, 6), param%drzkk(6, 6))
               case (7)  ! C,N
                  s_temp2 = s_temp2 + switch_T(coords(:, i), coords(:, j), param%rzkk(6, 7), param%drzkk(6, 7))
               case default
               end select
            end do
            surft%sk(i) = param%zk(6) &
                          + param%zkk(6, 6)*s_temp1 &
                          + param%zkk(6, 7)*(s_temp2**2)
            s_temp1 = 0.0_wp; s_temp2 = 0.0_wp
         case (7) ! N
            do j = 1, info%nat
               select case (info%Z(j))
               case (6)  ! N,C
                  add_temp = add_temp + switch_T(coords(:, i), coords(:, j), param%rzkk(6, 6), param%drzkk(6, 6))
                  do k = 1, info%nat
                     if ((k /= j) .and. (k /= i)) then
                       nc_temp = nc_temp + switch_T(coords(:, j), coords(:, k), param%rzkk(6, info%Z(k)), param%drzkk(6, info%Z(k)))
                     end if
                  end do
                  add_temp = add_temp*(nc_temp**2)
                  s_temp1 = s_temp1 + add_temp
                  add_temp = 0.0_wp; nc_temp = 0.0_wp
                  s_temp2 = s_temp2 + switch_T(coords(:, i), coords(:, j), param%rnc3, param%drnc3)
               case default
               end select
            end do
            surft%sk(i) = param%zk(7) &
                          + param%zkk(7, 6)*(s_temp1**1.3_wp) &
                          + param%nc3*s_temp2
            s_temp1 = 0.0_wp; s_temp2 = 0.0_wp
         case (8) ! O
            do j = 1, info%nat
               select case (info%Z(j))
               case (6)  ! O,C
                  s_temp1 = s_temp1 + switch_T(coords(:, i), coords(:, j), param%rzkk(8, 6), param%drzkk(8, 6))
               case (7)  ! O,N
                  s_temp2 = s_temp2 + switch_T(coords(:, i), coords(:, j), param%rzkk(8, 7), param%drzkk(8, 7))
               case (8)  ! O,O
                  if (i /= j) s_temp3 = s_temp3 + switch_T(coords(:, i), coords(:, j), param%rzkk(8, 8), param%drzkk(8, 8))
               case (15) ! O,P
                  s_temp4 = s_temp4 + switch_T(coords(:, i), coords(:, j), param%rzkk(8, 15), param%drzkk(8, 15))
               case default
               end select
            end do
            surft%sk(i) = param%zk(8) &
                          + param%zkk(8, 6)*s_temp1 &
                          + param%zkk(8, 7)*s_temp2 &
                          + param%zkk(8, 8)*s_temp3 &
                          + param%zkk(8, 15)*s_temp4
            s_temp1 = 0.0_wp; s_temp2 = 0.0_wp; s_temp3 = 0.0_wp; s_temp4 = 0.0_wp
         case default
            surft%sk(i) = param%zk(info%Z(i))
         end select
      end do
   end subroutine calc_surfacetensions

   !> Switching function used in SMD sigma_k definitions
   function switch_T(coord1, coord2, rzkk, drzkk) result(Tval)
      !> Atom coordinates [3]
      real(wp), intent(in) :: coord1(3), coord2(3)
      real(wp), intent(in) :: rzkk, drzkk
      real(wp) :: Tval
      real(wp) :: R, denom
      ! distances are converted to Å for the empirical switching definition
      R = sqrt((coord1(1) - coord2(1))**2 + (coord1(2) - coord2(2))**2 + (coord1(3) - coord2(3))**2)*autoaa
      if (R < (rzkk + drzkk)) then
         denom = R - drzkk - rzkk
         Tval = exp(drzkk/denom)
      else
         Tval = 0.0_wp
      end if
   end function switch_T

   subroutine init_info(species, ident, self)
      integer, intent(in) :: species(:)
      character(len=*), intent(in) :: ident(:)
      type(cds_info), intent(out) :: self
      integer :: elem
      self%nat = size(species)
      allocate (self%Z(self%nat))
      do elem = 1, self%nat
         self%Z(elem) = to_number(ident(species(elem)))
      end do
   end subroutine init_info

end module moist_model_component_cds_sigma
