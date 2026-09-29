!> Test suite for the cavity-surface adjoint accumulator
!>
!> - Typed storage: init sizes and zeroes every channel, re-init resizes,
!>   zero clears
!> - add_surface_weights sums repeated contributions channel by channel and
!>   leaves omitted channels untouched
!> - Error paths (`should_fail`): an uninitialized accumulator and mis-shaped
!>   vector or scalar weights each fail through their own named error
module test_cavity_surface_adjoint
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
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
   !>
   !> @param[out] error  test failure
   subroutine test_storage(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc

      call check(error, .not. acc%is_initialized(), "default surface adjoint claims to be initialized")
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
      acc%w_k2 = 2.0_wp
      call acc%zero()
      call check(error, all(acc%w_xyz == 0.0_wp) .and. all(acc%w_k2 == 0.0_wp))
   end subroutine test_storage

   !> add_surface_weights sums repeated contributions channel by channel
   !>
   !> - Every channel supplied once, the switching-factor channel twice
   !> - Omitted channels keep their previous value
   !>
   !> @param[out] error  test failure
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
      if (.not. allocated(err)) call acc%add_surface_weights(err, w_f=s)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, all(acc%w_xyz == v3) .and. all(acc%w_n == -v3) .and. all(acc%w_xi == s) &
         & .and. all(acc%w_f == 2.0_wp*s) .and. all(acc%w_k1 == 2.0_wp*s) .and. all(acc%w_k2 == 3.0_wp*s) &
         & .and. all(acc%w_a == 4.0_wp*s) .and. all(acc%w_w == 5.0_wp*s), &
         & "add_surface_weights did not accumulate channel by channel")
   end subroutine test_accumulate

   !> add_surface_weights on an uninitialized accumulator fails
   !>
   !> @param[out] error  test failure, set on the expected error
   subroutine test_uninit_fails(error)
      !> Test failure, set on the expected error
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: err

      call acc%add_surface_weights(err, w_a=[1.0_wp])
      call expect_error(error, err, "not initialized")
   end subroutine test_uninit_fails

   !> Position weights of the wrong shape fail
   !>
   !> @param[out] error  test failure, set on the expected error
   subroutine test_xyz_fails(error)
      !> Test failure, set on the expected error
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
   !>
   !> @param[out] error  test failure, set on the expected error
   subroutine test_normal_fails(error)
      !> Test failure, set on the expected error
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
   !>
   !> @param[out] error  test failure, set on the expected error
   subroutine test_area_fails(error)
      !> Test failure, set on the expected error
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: err

      call acc%init(3)
      call acc%add_surface_weights(err, w_a=[1.0_wp, 2.0_wp])
      call expect_error(error, err, "area weight size mismatch")
   end subroutine test_area_fails

end module test_cavity_surface_adjoint
