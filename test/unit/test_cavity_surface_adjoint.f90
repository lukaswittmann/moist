!> Test suite for the cavity-surface adjoint accumulator
!>
!> - Typed storage: init sizes and zeroes every channel, re-init resizes,
!>   zero clears
!> - add_surface_weights sums repeated contributions channel by channel and
!>   leaves omitted channels untouched
!> - is_initialized rejects storage with a missing channel or a channel whose
!>   size disagrees with the others
!> - add_surface_weights rejects a mis-shaped input on every channel, each with
!>   its own error message
!> - Error paths (`should_fail`): an uninitialized accumulator and mis-shaped
!>   vector or scalar weights each fail through their own named error
module test_cavity_surface_adjoint
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed, to_string
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   implicit none(type, external)
   private

   public :: collect_cavity_surface_adjoint

contains

   !> Collect all cavity_surface_adjoint tests
   !>
   !> @param[out] testsuite  collected unit tests
   subroutine collect_cavity_surface_adjoint(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("typed_storage", test_storage), &
                  new_unittest("add_weights_accumulates", test_accumulate), &
                  new_unittest("initialized_shapes", test_initialized_shapes), &
                  new_unittest("all_weight_shapes", test_all_weight_shapes), &
                  new_unittest("bad_uninitialized", test_uninit_fails, should_fail=.true.), &
                  new_unittest("bad_xyz_shape", test_xyz_fails, should_fail=.true.), &
                  new_unittest("bad_normal_shape", test_normal_fails, should_fail=.true.), &
                  new_unittest("bad_area_size", test_area_fails, should_fail=.true.) &
                  ]
   end subroutine collect_cavity_surface_adjoint

   !> Fail an expected-failure test only on the targeted library error
   !>
   !> - Used by `should_fail=.true.` tests, where test-drive inverts the
   !>   verdict: a raised test failure passes, a clean return fails
   !> - Fails the test only when `err` names `expected`; no error or a
   !>   different one returns cleanly, so the test is reported as failed
   !>
   !> @param[out] error     test failure, set only on the expected error
   !> @param[in]  err       library error, possibly unallocated
   !> @param[in]  expected  substring the expected error message must contain
   subroutine expect_error(error, err, expected)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error, possibly unallocated
      type(mctc_error), allocatable, intent(in) :: err
      !> Substring the expected error message must contain
      character(len=*), intent(in) :: expected

      if (.not. allocated(err)) return
      if (index(err%message, expected) > 0) call test_failed(error, err%message)
   end subroutine expect_error

   !> init sizes and zeroes every channel, re-init resizes, zero clears
   subroutine test_storage(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc

      call check(error, .not. acc%is_initialized(), "default surface adjoint claims to be initialized")
      if (allocated(error)) return
      call check(error, acc%size() == 0, "default accumulator size must be zero")
      if (allocated(error)) return
      call acc%init(3)
      call check(error, acc%is_initialized() .and. acc%size() == 3)
      if (allocated(error)) return
      call check(error, all(acc%w_xyz == 0.0_wp) .and. all(acc%w_n == 0.0_wp) &
         & .and. all(acc%w_w == 0.0_wp) .and. all(acc%w_xi == 0.0_wp) .and. all(acc%w_f == 0.0_wp) &
         & .and. all(acc%w_a == 0.0_wp) .and. all(acc%w_k1 == 0.0_wp) .and. all(acc%w_k2 == 0.0_wp))
      if (allocated(error)) return

      acc%w_a = 1.0_wp
      call acc%init(5)
      call check(error, acc%is_initialized() .and. acc%size() == 5 .and. all(acc%w_a == 0.0_wp), &
         & "re-init must resize and zero every channel")
      if (allocated(error)) return

      acc%w_xyz = 1.0_wp
      acc%w_n = 2.0_wp
      acc%w_w = 3.0_wp
      acc%w_xi = 4.0_wp
      acc%w_f = 5.0_wp
      acc%w_a = 6.0_wp
      acc%w_k1 = 7.0_wp
      acc%w_k2 = 8.0_wp
      call acc%zero()
      call check(error, all(acc%w_xyz == 0.0_wp) .and. all(acc%w_n == 0.0_wp) &
         & .and. all(acc%w_w == 0.0_wp) .and. all(acc%w_xi == 0.0_wp) .and. all(acc%w_f == 0.0_wp) &
         & .and. all(acc%w_a == 0.0_wp) .and. all(acc%w_k1 == 0.0_wp) .and. all(acc%w_k2 == 0.0_wp))
   end subroutine test_storage

   !> add_surface_weights sums repeated contributions channel by channel
   !>
   !> - Every channel supplied twice, the switching-factor channel three times
   !> - Omitted channels keep their previous value
   subroutine test_accumulate(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp) :: v3(3, 2), s(2)

      call acc%init(2)
      v3 = reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp], [3, 2])
      s = [0.5_wp, -1.5_wp]
      call acc%add_surface_weights(err, w_xi=s, w_f=s, w_xyz=v3, w_n=-v3, w_k1=2.0_wp*s, &
         & w_k2=3.0_wp*s, w_a=4.0_wp*s, w_w=5.0_wp*s)
      if (.not. allocated(err)) call acc%add_surface_weights(err, w_xi=s, w_f=s, w_xyz=v3, w_n=-v3, &
         & w_k1=2.0_wp*s, w_k2=3.0_wp*s, w_a=4.0_wp*s, w_w=5.0_wp*s)
      if (.not. allocated(err)) call acc%add_surface_weights(err, w_f=s)
      if (.not. allocated(err)) call acc%add_surface_weights(err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, all(acc%w_xyz == 2.0_wp*v3) .and. all(acc%w_n == -2.0_wp*v3) .and. all(acc%w_xi == 2.0_wp*s) &
         & .and. all(acc%w_f == 3.0_wp*s) .and. all(acc%w_k1 == 4.0_wp*s) .and. all(acc%w_k2 == 6.0_wp*s) &
         & .and. all(acc%w_a == 8.0_wp*s) .and. all(acc%w_w == 10.0_wp*s), &
         & "add_surface_weights did not accumulate channel by channel")
   end subroutine test_accumulate

   !> Missing channels and inconsistent dimensions invalidate storage
   subroutine test_initialized_shapes(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc
      integer :: channel

      do channel = 1, 8
         call acc%init(3)
         select case (channel)
         case (1)
            deallocate (acc%w_xyz)
         case (2)
            deallocate (acc%w_n)
         case (3)
            deallocate (acc%w_w)
         case (4)
            deallocate (acc%w_xi)
         case (5)
            deallocate (acc%w_f)
         case (6)
            deallocate (acc%w_a)
         case (7)
            deallocate (acc%w_k1)
         case (8)
            deallocate (acc%w_k2)
         case default
            call test_failed(error, "invalid surface channel selector")
            return
         end select
         call check(error, .not. acc%is_initialized(), &
                    "missing surface channel "//to_string(channel)//" accepted")
         if (allocated(error)) return
      end do
      do channel = 1, 10
         call acc%init(3)
         select case (channel)
         case (1)
            deallocate (acc%w_xyz)
            allocate (acc%w_xyz(3, 4), source=0.0_wp)
         case (2)
            deallocate (acc%w_n)
            allocate (acc%w_n(3, 4), source=0.0_wp)
         case (3)
            deallocate (acc%w_w)
            allocate (acc%w_w(4), source=0.0_wp)
         case (4)
            deallocate (acc%w_xi)
            allocate (acc%w_xi(4), source=0.0_wp)
         case (5)
            deallocate (acc%w_f)
            allocate (acc%w_f(4), source=0.0_wp)
         case (6)
            deallocate (acc%w_a)
            allocate (acc%w_a(4), source=0.0_wp)
         case (7)
            deallocate (acc%w_k1)
            allocate (acc%w_k1(4), source=0.0_wp)
         case (8)
            deallocate (acc%w_k2)
            allocate (acc%w_k2(4), source=0.0_wp)
         case (9)
            deallocate (acc%w_xyz)
            allocate (acc%w_xyz(2, 3), source=0.0_wp)
         case (10)
            deallocate (acc%w_n)
            allocate (acc%w_n(2, 3), source=0.0_wp)
         case default
            call test_failed(error, "invalid surface channel selector")
            return
         end select
         call check(error, .not. acc%is_initialized(), &
                    "inconsistent surface channel case "//to_string(channel)//" accepted")
         if (allocated(error)) return
      end do
   end subroutine test_initialized_shapes

   !> Every optional channel rejects mismatched input dimensions
   subroutine test_all_weight_shapes(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      integer :: channel
      real(wp) :: scalar(2), short_vector(2, 3), long_vector(3, 4)
      character(len=:), allocatable :: expected

      scalar = 1.0_wp
      short_vector = 1.0_wp
      long_vector = 1.0_wp
      do channel = 1, 10
         call acc%init(3)
         select case (channel)
         case (1)
            call acc%add_surface_weights(err, w_xi=scalar)
            expected = "xi weight size mismatch"
         case (2)
            call acc%add_surface_weights(err, w_f=scalar)
            expected = "f weight size mismatch"
         case (3)
            call acc%add_surface_weights(err, w_k1=scalar)
            expected = "k1 weight size mismatch"
         case (4)
            call acc%add_surface_weights(err, w_k2=scalar)
            expected = "k2 weight size mismatch"
         case (5)
            call acc%add_surface_weights(err, w_a=scalar)
            expected = "area weight size mismatch"
         case (6)
            call acc%add_surface_weights(err, w_w=scalar)
            expected = "integration weight size mismatch"
         case (7)
            call acc%add_surface_weights(err, w_xyz=short_vector)
            expected = "xyz weight shape mismatch"
         case (8)
            call acc%add_surface_weights(err, w_xyz=long_vector)
            expected = "xyz weight shape mismatch"
         case (9)
            call acc%add_surface_weights(err, w_n=short_vector)
            expected = "normal weight shape mismatch"
         case (10)
            call acc%add_surface_weights(err, w_n=long_vector)
            expected = "normal weight shape mismatch"
         case default
            call test_failed(error, "invalid surface channel selector")
            return
         end select
         call check(error, allocated(err), &
                    "mismatched surface input accepted in case "//to_string(channel))
         if (allocated(error)) return
         call check(error, index(err%message, expected) > 0, &
                    "wrong surface shape error in case "//to_string(channel)//": "//err%message)
         if (allocated(error)) return
      end do
   end subroutine test_all_weight_shapes

   !> add_surface_weights on an uninitialized accumulator fails
   subroutine test_uninit_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: err

      call acc%add_surface_weights(err, w_a=[1.0_wp])
      call expect_error(error, err, "not initialized")
   end subroutine test_uninit_fails

   !> Position weights of the wrong shape fail
   subroutine test_xyz_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp) :: w_xyz(3, 4)

      call acc%init(3)
      w_xyz = 1.0_wp
      call acc%add_surface_weights(err, w_xyz=w_xyz)
      call expect_error(error, err, "xyz weight shape mismatch")
   end subroutine test_xyz_fails

   !> Normal weights of the wrong shape fail
   subroutine test_normal_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp) :: w_n(2, 3)

      call acc%init(3)
      w_n = 1.0_wp
      call acc%add_surface_weights(err, w_n=w_n)
      call expect_error(error, err, "normal weight shape mismatch")
   end subroutine test_normal_fails

   !> Area weights of the wrong length fail
   subroutine test_area_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: err

      call acc%init(3)
      call acc%add_surface_weights(err, w_a=[1.0_wp, 2.0_wp])
      call expect_error(error, err, "area weight size mismatch")
   end subroutine test_area_fails

end module test_cavity_surface_adjoint
