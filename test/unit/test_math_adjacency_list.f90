!> Test suite for adjacency-list utilities in moist_math_adjacency_list
module test_math_adjacency_list
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_adjacency_list, only: adjacency_list_type
   implicit none (type, external)
   private

   public :: collect_math_adjacency_list

   real(wp), parameter :: thr = 10.0_wp * epsilon(1.0_wp)

contains

   !> Collect all adjacency-list tests.
   subroutine collect_math_adjacency_list(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("pair_membership_bruteforce", test_pair_membership_bruteforce), &
         new_unittest("distances_disabled_neighbour_content", test_neighbour_content) &
      ]
   end subroutine collect_math_adjacency_list

   !> Compare adjacency membership against brute-force distance checks.
   subroutine test_pair_membership_bruteforce(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 6)
      real(wp), parameter :: cutoff = 1.5_wp
      real(wp) :: cutoff2, d2, tol
      integer, allocatable :: ids(:)
      integer :: i, j, expected_count
      logical :: in_expected, in_list

      ! Mix of inside, outside, and exact-cutoff pairs.
      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.5_wp, 0.0_wp, 0.0_wp]     ! exactly at cutoff from 1
      xyz(:, 3) = [1.5002_wp, 0.0_wp, 0.0_wp]  ! just outside cutoff from 1
      xyz(:, 4) = [0.0_wp, 1.0_wp, 0.0_wp]
      xyz(:, 5) = [0.0_wp, 0.0_wp, 2.0_wp]
      xyz(:, 6) = [0.0_wp, 1.5_wp, 0.0_wp]     ! exactly at cutoff from 1

      call nlist%init(cutoff=cutoff)
      call nlist%update(xyz)

      cutoff2 = cutoff * cutoff
      tol = 10.0_wp * epsilon(cutoff2)

      do i = 1, size(xyz, 2)
         ids = nlist%get_neighbours(i)

         expected_count = 0
         do j = 1, size(xyz, 2)
            if (j == i) cycle
            d2 = sum((xyz(:, i) - xyz(:, j))**2)
            if (d2 <= cutoff2 + tol) expected_count = expected_count + 1
         end do

         call check(error, size(ids) == expected_count, &
            "Neighbour count mismatch against brute-force reference")
         if (allocated(error)) return

         do j = 1, size(xyz, 2)
            in_list = any(ids == j)

            if (j == i) then
               call check(error, .not. in_list, "Self index must never appear in neighbour list")
            else
               d2 = sum((xyz(:, i) - xyz(:, j))**2)
               in_expected = (d2 <= cutoff2 + tol)
               call check(error, in_list .eqv. in_expected, &
                  "Pair membership mismatch against brute-force reference")
            end if
            if (allocated(error)) return
         end do
      end do

      call nlist%destroy()
   end subroutine test_pair_membership_bruteforce


   !> Neighbour list content should stay correct when distance storage is disabled.
   subroutine test_neighbour_content(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 4)
      integer, allocatable :: ids(:)

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [3.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 4) = [0.0_wp, 4.0_wp, 0.0_wp]

      call nlist%init(cutoff=2.1_wp)
      call nlist%update(xyz)

      ids = nlist%get_neighbours(1)
      call check(error, size(ids) == 1, "Point 1 should have one neighbour")
      if (allocated(error)) return
      call check(error, any(ids == 2), "Point 1 should be connected to point 2")
      if (allocated(error)) return

      ids = nlist%get_neighbours(2)
      call check(error, size(ids) == 2, "Point 2 should have two neighbours")
      if (allocated(error)) return
      call check(error, any(ids == 1) .and. any(ids == 3), &
         "Point 2 should be connected to points 1 and 3")
      if (allocated(error)) return

      ids = nlist%get_neighbours(3)
      call check(error, size(ids) == 1, "Point 3 should have one neighbour")
      if (allocated(error)) return
      call check(error, any(ids == 2), "Point 3 should be connected to point 2")
      if (allocated(error)) return

      ids = nlist%get_neighbours(4)
      call check(error, size(ids) == 0, "Point 4 should have no neighbours")
      if (allocated(error)) return

      call nlist%destroy()
   end subroutine test_neighbour_content

end module test_math_adjacency_list
