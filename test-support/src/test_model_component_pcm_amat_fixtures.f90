!> Shared fixtures for the Gaussian PCM interaction-matrix test suites
!>
!> Split out of `test_model_component_pcm_amat` so the assembly and adjoint
!> suites, which sit one directory below it, can reach them: fpm only lets a
!> test source use modules from its own directory or below.
module test_model_component_pcm_amat_fixtures
   use mctc_env, only: wp
   use moist_model_component_pcm_amat_kernel, only: pcm_amat_x_far
   implicit none (type, external)
   private

   public :: count_branches
   public :: nmol, nleb_survey

   !> Number of mstore structures drawn per test
   integer, parameter :: nmol = 10

   !> Lebedev grid used for the survey tests
   integer, parameter :: nleb_survey = 50

contains

   !> Count the pairs on either side of the erf-saturation threshold
   !>
   !> @param[in]  xi    Gaussian widths (ngrid)
   !> @param[in]  xyz   Surface positions (3, ngrid)
   !> @param[out] nfar  Number of saturated pairs
   !> @param[out] nnear Number of unsaturated pairs
   subroutine count_branches(xi, xyz, nfar, nnear)
      !> Gaussian widths and surface positions
      real(wp), intent(in) :: xi(:), xyz(:, :)
      !> Saturated and unsaturated pair counts
      integer, intent(out) :: nfar, nnear

      !> Pair indices and grid size
      integer :: i, j, ngrid
      !> Squared separation
      real(wp) :: r2
      !> Per-point saturation bounds
      real(wp), allocatable :: bound(:)

      ngrid = size(xi)
      allocate (bound(ngrid))
      bound = pcm_amat_x_far/(xi*xi)

      nfar = 0
      nnear = 0
      do i = 1, ngrid
         do j = 1, i - 1
            r2 = sum((xyz(:, i) - xyz(:, j))**2)
            if (r2 >= bound(i) + bound(j)) then
               nfar = nfar + 1
            else
               nnear = nnear + 1
            end if
         end do
      end do

   end subroutine count_branches

end module test_model_component_pcm_amat_fixtures
