!> Interface-level unit tests and shared harness for the Gaussian PCM interaction matrix
module test_model_component_pcm_amat
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, &
                                             pcm_amat_surface_weights
   use moist_model_component_pcm_amat_kernel, only: pcm_amat_x_far
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none (type, external)
   private

   public :: collect_model_component_pcm_amat
   public :: count_branches
   public :: nmol, nleb_survey

   !> Number of mstore structures drawn per test
   integer, parameter :: nmol = 10

   !> Lebedev grid used for the survey tests
   integer, parameter :: nleb_survey = 50

contains

   !> Collect the interface-level test suite
   subroutine collect_model_component_pcm_amat(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rejects_invalid_surface", test_rejects_invalid_surface) &
                  ]

   end subroutine collect_model_component_pcm_amat

   !* --------------------------------- Shared helpers -------------------------------- *!

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

   !* ------------------------------------- Tests ------------------------------------- *!

   !> Malformed surfaces are rejected and leave a defined output behind
   subroutine test_rejects_invalid_surface(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Synthetic surface
      real(wp) :: xi(3), f(3), xyz(3, 3), amat(3, 3), amat_small(2, 2)
      !> Contraction vectors and adjoints
      real(wp) :: q1(3), q2(3), w_xi(3), w_f(3), w_xyz(3, 3)

      xi = [1.5_wp, 2.0_wp, 0.9_wp]
      f = [0.8_wp, 0.95_wp, 0.6_wp]
      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.7_wp, -0.4_wp, 0.9_wp]
      xyz(:, 3) = [-3.1_wp, 2.2_wp, -1.4_wp]
      q1 = [0.3_wp, -0.2_wp, 0.1_wp]
      q2 = [0.1_wp, 0.4_wp, -0.3_wp]

      ! A non-positive width has no pair width, so the kernel is undefined.
      amat = 1.0_wp
      call assemble_pcm_amat([1.5_wp, 0.0_wp, 0.9_wp], f, xyz, amat, err)
      call check(error, allocated(err), more="a zero Gaussian width was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(amat)), 0.0_wp, thr=0.0_wp, &
                 more="rejected assembly left the matrix untouched")
      if (allocated(error)) return
      deallocate (err)

      ! A non-positive switching factor divides the diagonal by zero.
      amat = 1.0_wp
      call assemble_pcm_amat(xi, [0.8_wp, -0.1_wp, 0.6_wp], xyz, amat, err)
      call check(error, allocated(err), more="a negative switching factor was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(amat)), 0.0_wp, thr=0.0_wp, &
                 more="rejected assembly left the matrix untouched")
      if (allocated(error)) return
      deallocate (err)

      ! The output matrix must match the surface it is assembled on.
      amat_small = 1.0_wp
      call assemble_pcm_amat(xi, f, xyz, amat_small, err)
      call check(error, allocated(err), more="a mis-shaped output matrix was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(amat_small)), 0.0_wp, thr=0.0_wp, &
                 more="rejected assembly left the matrix untouched")
      if (allocated(error)) return
      deallocate (err)

      ! The adjoint route validates the same surface plus its own vectors.
      w_xi = 1.0_wp
      w_f = 1.0_wp
      w_xyz = 1.0_wp
      call pcm_amat_surface_weights(xi, f, xyz, q1(1:2), q2, w_xi, w_f, w_xyz, err)
      call check(error, allocated(err), &
                 more="a mis-shaped contraction vector was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(w_xi)) + maxval(abs(w_f)) + maxval(abs(w_xyz)), &
                 0.0_wp, thr=0.0_wp, more="rejected contraction left the weights untouched")

   end subroutine test_rejects_invalid_surface

end module test_model_component_pcm_amat
