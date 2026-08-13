!> Tests for the DROP grid filter and branch-weight computation
module test_cavity_drop_filter
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use moist_context, only: moist_context_type
   use moist_cavity_drop, only: cavity_type_drop
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none(type, external)
   private

   public :: collect_cavity_drop_filter

   !> Array identities used by the tag helpers
   integer, parameter :: SLOT_XYZ = 1
   integer, parameter :: SLOT_ANCHORXYZ = 2
   integer, parameter :: SLOT_NORMAL0 = 3
   integer, parameter :: SLOT_ANCHOR_WLEB0 = 4
   integer, parameter :: SLOT_LAMBDA0 = 5
   integer, parameter :: SLOT_ISWIG_F0 = 6
   integer, parameter :: SLOT_ANCHOR_XI0 = 7
   integer, parameter :: SLOT_RHO = 8
   integer, parameter :: SLOT_R_II0 = 9
   integer, parameter :: SLOT_WBRANCH = 10
   integer, parameter :: SLOT_PHI0 = 11
   integer, parameter :: SLOT_CPJAC = 12
   integer, parameter :: SLOT_W_F0 = 13
   integer, parameter :: SLOT_NUMBERING = 14
   integer, parameter :: SLOT_OWNER = 15
   integer, parameter :: SLOT_BRANCH = 16
   integer, parameter :: SLOT_ANCHOR_ID = 17
   integer, parameter :: SLOT_BRANCH_COUNT = 18

   !> Softmax weights theshold
   real(wp), parameter :: SOFTMAX_THR = 1.0e-14_wp

contains

   !> Collect all DROP filter tests
   subroutine collect_cavity_drop_filter(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("filter_permutation_consistency", test_filter_permutation), &
                  new_unittest("filter_compacts_every_array", test_filter_size_audit), &
                  new_unittest("filter_cutoff_boundary", test_filter_cutoff_boundary), &
                  new_unittest("filter_keeps_all", test_filter_keeps_all), &
                  new_unittest("filter_keeps_none", test_filter_keeps_none), &
                  new_unittest("filter_skips_unallocated", test_filter_skips_unallocated), &
                  new_unittest("filter_reports_no_error", test_filter_reports_no_error), &
                  new_unittest("branch_partition_of_unity", test_branch_partition_of_unity), &
                  new_unittest("branch_folds_into_wleb", test_branch_folds_into_wleb), &
                  new_unittest("branch_singletons_untouched", test_branch_singletons_untouched), &
                  new_unittest("branch_uniform_phi", test_branch_uniform_phi), &
                  new_unittest("branch_dominant_sibling", test_branch_dominant_sibling), &
                  new_unittest("branch_all_below_cut_fallback", test_branch_fallback), &
                  new_unittest("branch_softmax_width_limits", test_branch_width_limits) &
                  ]

   end subroutine collect_cavity_drop_filter

   !* ================================================================================= *!
   !*                                   filter_arrays                                   *!
   !* ================================================================================= *!

   !> A scattered keep mask must compact every array with the same permutation:
   !> each surviving point still carries all of its own field values
   subroutine test_filter_permutation(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 9
      !> Points 1, 4, 5 and 9 fall below the cutoff: leading, trailing and a run
      integer, parameter :: survivor(5) = [2, 3, 6, 7, 8]
      real(wp) :: wleb(n), f(n)
      integer :: i, k

      do i = 1, n
         wleb(i) = real(i, wp)
         f(i) = 1.0_wp + 0.01_wp*real(i, wp)
      end do
      !> Kept points sit at wleb*f >= 2.02, dropped ones at <= 9.0e-6
      f([1, 4, 5, 9]) = 1.0e-6_wp

      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 0.5_wp

      call cav%filter_arrays("Test", failed)
      call check(error, .not. allocated(failed), "filtering succeeds")
      if (allocated(error)) return

      call check(error, cav%ngrid, size(survivor), "surviving point count")
      if (allocated(error)) return

      do k = 1, size(survivor)
         i = survivor(k)
         call check_point(error, cav, k, i, wleb(i), f(i))
         if (allocated(error)) return
      end do

   end subroutine test_filter_permutation

   !> Every array listed in `compact_grid_arrays` comes back at the new length
   subroutine test_filter_size_audit(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 6
      integer, parameter :: nkept = 3
      real(wp) :: wleb(n), f(n)

      wleb = 1.0_wp
      f = 1.0_wp
      !> Keep points 2, 3 and 5
      f([1, 4, 6]) = 0.0_wp

      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 0.5_wp

      call cav%filter_arrays("Test", failed)
      call check(error, .not. allocated(failed), "filtering succeeds")
      if (allocated(error)) return

      call check(error, cav%ngrid, nkept, "ngrid follows the keep count")
      if (allocated(error)) return

      call check_grid_sizes(error, cav, nkept)

   end subroutine test_filter_size_audit

   !> The keep criterion is `wleb*f > wleb_cut`, strictly and on the product.
   subroutine test_filter_cutoff_boundary(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 4
      real(wp) :: wleb(n), f(n)

      !> 1: product exactly at the cutoff, 2: just above, 3: just below,
      !> 4: wleb far above the cutoff but killed by a tiny switching value
      wleb = [0.5_wp, 0.5_wp, 0.5_wp, 100.0_wp]
      f = [0.5_wp, 0.5000001_wp, 0.4999999_wp, 1.0e-6_wp]

      call build_cavity(cav, ctx, n, wleb, f)
      !> 0.5*0.5 is exact in binary, so point 1 sits precisely on the cutoff
      cav%param%wleb_cut = 0.25_wp

      call cav%filter_arrays("Test", failed)
      call check(error, .not. allocated(failed), "filtering succeeds")
      if (allocated(error)) return

      call check(error, cav%ngrid, 1, "only the point strictly above survives")
      if (allocated(error)) return

      call check_point(error, cav, 1, 2, wleb(2), f(2))

   end subroutine test_filter_cutoff_boundary

   !> Keeping everything leaves the grid untouched
   subroutine test_filter_keeps_all(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 5
      real(wp) :: wleb(n), f(n)
      integer :: i

      wleb = 1.0_wp
      f = 1.0_wp

      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 1.0e-12_wp

      call cav%filter_arrays("Test", failed)
      call check(error, .not. allocated(failed), "filtering succeeds")
      if (allocated(error)) return

      call check(error, cav%ngrid, n, "no point removed")
      if (allocated(error)) return

      call check_grid_sizes(error, cav, n)
      if (allocated(error)) return

      do i = 1, n
         call check_point(error, cav, i, i, wleb(i), f(i))
         if (allocated(error)) return
      end do

   end subroutine test_filter_keeps_all

   !> An empty result is a valid one: zero-length arrays, ngrid zero, no crash,
   !> and a second pass over the empty grid still behaves
   subroutine test_filter_keeps_none(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 4
      real(wp) :: wleb(n), f(n)

      wleb = 1.0_wp
      f = 1.0e-30_wp

      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 1.0e-12_wp

      call cav%filter_arrays("Test", failed)
      call check(error, .not. allocated(failed), "filtering succeeds")
      if (allocated(error)) return

      call check(error, cav%ngrid, 0, "every point removed")
      if (allocated(error)) return

      call check_grid_sizes(error, cav, 0)
      if (allocated(error)) return

      call cav%filter_arrays("Test", failed)
      call check(error, .not. allocated(failed), "filtering an empty grid succeeds")
      if (allocated(error)) return

      call check(error, cav%ngrid, 0, "empty grid stays empty")

   end subroutine test_filter_keeps_none

   !> Arrays that do not exist yet at the pre-filter stage are skipped rather
   !> than faulting, and stay unallocated
   subroutine test_filter_skips_unallocated(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 4
      integer :: i

      ctx%verbosity = 0
      cav%ctx => ctx
      cav%ngrid = n
      cav%param%wleb_cut = 0.5_wp

      !> Only the subset that exists before projection
      allocate (cav%xyz(3, n), cav%wleb(n), cav%f(n), cav%owner(n))
      do i = 1, n
         cav%xyz(:, i) = [tag2(SLOT_XYZ, i, 1), tag2(SLOT_XYZ, i, 2), tag2(SLOT_XYZ, i, 3)]
         cav%owner(i) = itag(SLOT_OWNER, i)
      end do
      cav%wleb = 1.0_wp
      cav%f = [1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]

      call cav%filter_arrays("Test", failed)
      call check(error, .not. allocated(failed), "filtering a partial cavity succeeds")
      if (allocated(error)) return

      call check(error, cav%ngrid, 2, "surviving point count")
      if (allocated(error)) return

      call check(error, size(cav%xyz, 2), 2, "xyz compacted")
      if (allocated(error)) return
      call check(error, cav%xyz(1, 2), tag2(SLOT_XYZ, 3, 1), "second survivor is point 3")
      if (allocated(error)) return
      call check(error, cav%owner(2), itag(SLOT_OWNER, 3), "owner follows the same permutation")
      if (allocated(error)) return

      call check(error, .not. allocated(cav%numbering), "numbering stays unallocated")
      if (allocated(error)) return
      call check(error, .not. allocated(cav%phi0), "phi0 stays unallocated")
      if (allocated(error)) return
      call check(error, .not. allocated(cav%converged), "converged stays unallocated")

   end subroutine test_filter_skips_unallocated

   !> Neither entry point has a failure mode; both leave `error` unallocated
   subroutine test_filter_reports_no_error(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 3
      real(wp) :: wleb(n), f(n)

      wleb = 1.0_wp
      f = 1.0_wp

      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 1.0e-12_wp
      call cav%branch_weight%init(1.0_wp)
      cav%anchor_id = [1, 2, 2]
      cav%branch_count = [1, 2, 2]
      cav%phi0 = [0.0_wp, 0.0_wp, 0.1_wp]

      call cav%compute_branch_weights(failed)
      call check(error, .not. allocated(failed), "branch weighting reports no error")
      if (allocated(error)) return

      call cav%filter_arrays("Test", failed)
      call check(error, .not. allocated(failed), "filtering reports no error")

   end subroutine test_filter_reports_no_error

   !* ================================================================================= *!
   !*                              compute_branch_weights                               *!
   !* ================================================================================= *!

   !> Branch weights of the kept siblings sum to one, both for a group that
   !> survives intact and for one that loses a sibling and is re-normalized
   subroutine test_branch_partition_of_unity(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 8
      real(wp), parameter :: sigma = 1.0_wp
      real(wp) :: wleb(n), f(n), wref(3)

      wleb = 1.0_wp
      f = 1.0_wp

      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 1.0e-3_wp
      call cav%branch_weight%init(sigma)

      !> Singleton, a three-fold group, singleton, a two-fold group whose
      !> second sibling is far too weak to survive, singleton
      cav%anchor_id = [10, 11, 11, 11, 12, 13, 13, 14]
      cav%branch_count = [1, 3, 3, 3, 1, 2, 2, 1]
      cav%phi0 = [0.0_wp, 0.0_wp, 0.2_wp, 0.5_wp, 0.0_wp, 0.0_wp, 12.0_wp, 0.0_wp]

      call cav%compute_branch_weights(failed)
      call check(error, .not. allocated(failed), "branch weighting succeeds")
      if (allocated(error)) return

      !> Intact group: weights match the plain softmax and sum to one
      wref = ref_softmax(cav%phi0(2:4), sigma)
      call check(error, cav%wbranch(2), wref(1), thr_abs=SOFTMAX_THR, thr_rel=SOFTMAX_THR)
      if (allocated(error)) return
      call check(error, cav%wbranch(3), wref(2), thr_abs=SOFTMAX_THR, thr_rel=SOFTMAX_THR)
      if (allocated(error)) return
      call check(error, cav%wbranch(4), wref(3), thr_abs=SOFTMAX_THR, thr_rel=SOFTMAX_THR)
      if (allocated(error)) return
      call check(error, sum(cav%wbranch(2:4)), 1.0_wp, thr_abs=SOFTMAX_THR, thr_rel=SOFTMAX_THR)
      if (allocated(error)) return
      call check(error, cav%branch_count(2), 3, "intact group keeps its branch count")
      if (allocated(error)) return

      !> Pruned group: the survivor carries the whole anchor weight
      call check(error, cav%wbranch(6), 1.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, cav%wbranch(7), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, cav%wleb(7), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, cav%branch_count(6), 1, "collapsed group reports one branch")
      if (allocated(error)) return
      call check(error, cav%branch_count(7), 1, "dropped sibling reports the kept count")

   end subroutine test_branch_partition_of_unity

   !> The branch weight is folded into wleb by an exact multiplication
   subroutine test_branch_folds_into_wleb(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 5
      real(wp) :: wleb(n), f(n), before(n)
      integer :: i

      do i = 1, n
         wleb(i) = 0.1_wp*real(i, wp)
      end do
      f = 1.0_wp
      before = wleb

      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 1.0e-3_wp
      call cav%branch_weight%init(1.0_wp)

      cav%anchor_id = [20, 20, 20, 21, 21]
      cav%branch_count = [3, 3, 3, 2, 2]
      cav%phi0 = [0.0_wp, 0.3_wp, 0.7_wp, 0.0_wp, 20.0_wp]

      call cav%compute_branch_weights(failed)
      call check(error, .not. allocated(failed), "branch weighting succeeds")
      if (allocated(error)) return

      do i = 1, n
         call check(error, cav%wleb(i), before(i)*cav%wbranch(i), &
                    thr_abs=0.0_wp, thr_rel=0.0_wp)
         if (allocated(error)) return
      end do

      !> The pruned sibling is zeroed, which is what makes the follow-up
      !> filter_arrays call drop it
      call check(error, cav%wleb(5), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)

   end subroutine test_branch_folds_into_wleb

   !> Unbranched points are skipped outright: no weight, no scaling, no count
   !> change, whatever their objective values look like
   subroutine test_branch_singletons_untouched(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 4
      real(wp) :: wleb(n), f(n), before(n)
      integer :: i

      do i = 1, n
         wleb(i) = 0.25_wp*real(i, wp)
      end do
      f = 1.0_wp
      before = wleb

      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 1.0e-3_wp
      call cav%branch_weight%init(1.0_wp)

      cav%anchor_id = [30, 31, 32, 33]
      cav%branch_count = [1, 1, 1, 1]
      cav%phi0 = [0.0_wp, 5.0_wp, -2.0_wp, 100.0_wp]

      call cav%compute_branch_weights(failed)
      call check(error, .not. allocated(failed), "branch weighting succeeds")
      if (allocated(error)) return

      do i = 1, n
         call check(error, cav%wleb(i), before(i), thr_abs=0.0_wp, thr_rel=0.0_wp)
         if (allocated(error)) return
         call check(error, cav%wbranch(i), tag1(SLOT_WBRANCH, i), &
                    thr_abs=0.0_wp, thr_rel=0.0_wp)
         if (allocated(error)) return
         call check(error, cav%branch_count(i), 1, "singleton branch count unchanged")
         if (allocated(error)) return
      end do

   end subroutine test_branch_singletons_untouched

   !> Degenerate siblings share the anchor weight equally
   subroutine test_branch_uniform_phi(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 4
      real(wp) :: wleb(n), f(n)
      integer :: i

      wleb = 1.0_wp
      f = 1.0_wp

      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 1.0e-3_wp
      call cav%branch_weight%init(1.0_wp)

      cav%anchor_id = [40, 40, 40, 40]
      cav%branch_count = [4, 4, 4, 4]
      cav%phi0 = 0.75_wp

      call cav%compute_branch_weights(failed)
      call check(error, .not. allocated(failed), "branch weighting succeeds")
      if (allocated(error)) return

      do i = 1, n
         call check(error, cav%wbranch(i), 0.25_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
         if (allocated(error)) return
         call check(error, cav%branch_count(i), 4, "all siblings retained")
         if (allocated(error)) return
      end do

   end subroutine test_branch_uniform_phi

   !> One clearly better branch takes the whole anchor weight and the group
   !> collapses to a single point
   subroutine test_branch_dominant_sibling(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 3
      real(wp) :: wleb(n), f(n)
      integer :: i

      wleb = 1.0_wp
      f = 1.0_wp

      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 1.0e-3_wp
      call cav%branch_weight%init(1.0_wp)

      cav%anchor_id = [50, 50, 50]
      cav%branch_count = [3, 3, 3]
      cav%phi0 = [0.0_wp, 30.0_wp, 30.0_wp]

      call cav%compute_branch_weights(failed)
      call check(error, .not. allocated(failed), "branch weighting succeeds")
      if (allocated(error)) return

      call check(error, cav%wbranch(1), 1.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, cav%wleb(1), 1.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return

      do i = 2, 3
         call check(error, cav%wbranch(i), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
         if (allocated(error)) return
         call check(error, cav%wleb(i), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
         if (allocated(error)) return
      end do

      do i = 1, n
         call check(error, cav%branch_count(i), 1, "group collapsed to one branch")
         if (allocated(error)) return
      end do

   end subroutine test_branch_dominant_sibling

   !> When the cutoff would throw away every sibling, the strongest one is kept
   !> instead of discarding the anchor entirely
   subroutine test_branch_fallback(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 4
      real(wp) :: wleb(n), f(n)
      integer :: i

      wleb = 1.0_wp
      f = 1.0_wp

      call build_cavity(cav, ctx, n, wleb, f)
      !> Softmax weights of the group below are 0.087/0.644/0.237/0.032, so a
      !> cutoff of 0.7 rejects all four and forces the safety fallback
      cav%param%wleb_cut = 0.7_wp
      call cav%branch_weight%init(1.0_wp)

      cav%anchor_id = [60, 60, 60, 60]
      cav%branch_count = [4, 4, 4, 4]
      cav%phi0 = [3.0_wp, 1.0_wp, 2.0_wp, 4.0_wp]

      call cav%compute_branch_weights(failed)
      call check(error, .not. allocated(failed), "branch weighting succeeds")
      if (allocated(error)) return

      !> The lowest objective value wins the maxloc pick
      call check(error, cav%wbranch(2), 1.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, cav%wleb(2), 1.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return

      do i = 1, n
         if (i == 2) cycle
         call check(error, cav%wbranch(i), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
         if (allocated(error)) return
      end do

      call check(error, sum(cav%wbranch(1:4)), 1.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return

      do i = 1, n
         call check(error, cav%branch_count(i), 1, "fallback keeps a single branch")
         if (allocated(error)) return
      end do

   end subroutine test_branch_fallback

   !> The softmax width interpolates between argmax and uniform weighting
   subroutine test_branch_width_limits(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cav
      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: failed
      integer, parameter :: n = 2
      real(wp) :: wleb(n), f(n)

      wleb = 1.0_wp
      f = 1.0_wp

      !> A zero cutoff keeps both siblings alive so the raw weights are visible
      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 0.0_wp
      call cav%branch_weight%init(0.01_wp)

      cav%anchor_id = [70, 70]
      cav%branch_count = [2, 2]
      cav%phi0 = [0.0_wp, 1.0_wp]

      call cav%compute_branch_weights(failed)
      call check(error, .not. allocated(failed), "narrow softmax succeeds")
      if (allocated(error)) return

      call check(error, cav%wbranch(1), 1.0_wp, thr_abs=1.0e-30_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, cav%wbranch(2) < 1.0e-40_wp, "narrow softmax suppresses the loser")
      if (allocated(error)) return
      call check(error, cav%wbranch(2) > 0.0_wp, "the loser keeps a positive weight")
      if (allocated(error)) return

      !> A wide softmax washes the objective difference out
      call build_cavity(cav, ctx, n, wleb, f)
      cav%param%wleb_cut = 0.0_wp
      call cav%branch_weight%init(1.0e6_wp)

      cav%anchor_id = [70, 70]
      cav%branch_count = [2, 2]
      cav%phi0 = [0.0_wp, 1.0_wp]

      call cav%compute_branch_weights(failed)
      call check(error, .not. allocated(failed), "wide softmax succeeds")
      if (allocated(error)) return

      call check(error, cav%wbranch(1), 0.5_wp, thr_abs=1.0e-6_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, cav%wbranch(2), 0.5_wp, thr_abs=1.0e-6_wp, thr_rel=0.0_wp)

   end subroutine test_branch_width_limits

   !* ================================================================================= *!
   !*                                      Fixture                                      *!
   !* ================================================================================= *!

   !> Build a cavity carrying every array `compact_grid_arrays` filters
   !>
   !> All arrays except `wleb` and `f` are tagged with `tag*`, so a surviving
   !> point can be traced back to its original grid index field by field.
   !>
   !> @param[out]   cav  Cavity to populate
   !> @param[inout] ctx  Run context the cavity borrows (silenced here)
   !> @param[in]    n    Number of grid points
   !> @param[in]    wleb Quadrature weights driving the keep mask
   !> @param[in]    f    Switching values driving the keep mask
   subroutine build_cavity(cav, ctx, n, wleb, f)
      !> Cavity to populate
      type(cavity_type_drop), intent(out) :: cav
      !> Borrowed run context
      type(moist_context_type), intent(inout), target :: ctx
      !> Number of grid points
      integer, intent(in) :: n
      !> Quadrature weights and switching values
      real(wp), intent(in) :: wleb(:), f(:)

      integer :: i, j

      ctx%verbosity = 0
      cav%ctx => ctx
      cav%ngrid = n

      allocate (cav%xyz(3, n), cav%anchorxyz(3, n), cav%normal0(3, n))
      allocate (cav%wleb(n), cav%anchor_wleb0(n), cav%lambda0(n), cav%iswig_f0(n))
      allocate (cav%f(n), cav%anchor_xi0(n), cav%rho(n), cav%r_iI0(n))
      allocate (cav%wbranch(n), cav%phi0(n), cav%cpjac_scal0(n), cav%w_f0(n))
      allocate (cav%numbering(n), cav%owner(n), cav%branch(n), cav%anchor_id(n))
      allocate (cav%branch_count(n), cav%converged(n))

      do i = 1, n
         do j = 1, 3
            cav%xyz(j, i) = tag2(SLOT_XYZ, i, j)
            cav%anchorxyz(j, i) = tag2(SLOT_ANCHORXYZ, i, j)
            cav%normal0(j, i) = tag2(SLOT_NORMAL0, i, j)
         end do

         cav%wleb(i) = wleb(i)
         cav%f(i) = f(i)

         cav%anchor_wleb0(i) = tag1(SLOT_ANCHOR_WLEB0, i)
         cav%lambda0(i) = tag1(SLOT_LAMBDA0, i)
         cav%iswig_f0(i) = tag1(SLOT_ISWIG_F0, i)
         cav%anchor_xi0(i) = tag1(SLOT_ANCHOR_XI0, i)
         cav%rho(i) = tag1(SLOT_RHO, i)
         cav%r_iI0(i) = tag1(SLOT_R_II0, i)
         cav%wbranch(i) = tag1(SLOT_WBRANCH, i)
         cav%phi0(i) = tag1(SLOT_PHI0, i)
         cav%cpjac_scal0(i) = tag1(SLOT_CPJAC, i)
         cav%w_f0(i) = tag1(SLOT_W_F0, i)

         cav%numbering(i) = itag(SLOT_NUMBERING, i)
         cav%owner(i) = itag(SLOT_OWNER, i)
         cav%branch(i) = itag(SLOT_BRANCH, i)
         cav%anchor_id(i) = itag(SLOT_ANCHOR_ID, i)
         cav%branch_count(i) = itag(SLOT_BRANCH_COUNT, i)

         cav%converged(i) = mod(i, 3) == 0
      end do

   end subroutine build_cavity

   !> Assert that grid point `inew` carries every field of original point `iold`
   !>
   !> @param[out] error    Error handle
   !> @param[in]  cav      Filtered cavity
   !> @param[in]  inew     Index in the compacted grid
   !> @param[in]  iold     Index the point had before filtering
   !> @param[in]  wleb_ref Quadrature weight the point started with
   !> @param[in]  f_ref    Switching value the point started with
   subroutine check_point(error, cav, inew, iold, wleb_ref, f_ref)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Filtered cavity
      type(cavity_type_drop), intent(in) :: cav
      !> Compacted and original grid index
      integer, intent(in) :: inew, iold
      !> Untagged values the point started with
      real(wp), intent(in) :: wleb_ref, f_ref

      integer :: j

      do j = 1, 3
         call check(error, cav%xyz(j, inew), tag2(SLOT_XYZ, iold, j), "xyz")
         if (allocated(error)) return
         call check(error, cav%anchorxyz(j, inew), tag2(SLOT_ANCHORXYZ, iold, j), "anchorxyz")
         if (allocated(error)) return
         call check(error, cav%normal0(j, inew), tag2(SLOT_NORMAL0, iold, j), "normal0")
         if (allocated(error)) return
      end do

      call check(error, cav%wleb(inew), wleb_ref, "wleb")
      if (allocated(error)) return
      call check(error, cav%f(inew), f_ref, "f")
      if (allocated(error)) return

      call check(error, cav%anchor_wleb0(inew), tag1(SLOT_ANCHOR_WLEB0, iold), "anchor_wleb0")
      if (allocated(error)) return
      call check(error, cav%lambda0(inew), tag1(SLOT_LAMBDA0, iold), "lambda0")
      if (allocated(error)) return
      call check(error, cav%iswig_f0(inew), tag1(SLOT_ISWIG_F0, iold), "iswig_f0")
      if (allocated(error)) return
      call check(error, cav%anchor_xi0(inew), tag1(SLOT_ANCHOR_XI0, iold), "anchor_xi0")
      if (allocated(error)) return
      call check(error, cav%rho(inew), tag1(SLOT_RHO, iold), "rho")
      if (allocated(error)) return
      call check(error, cav%r_iI0(inew), tag1(SLOT_R_II0, iold), "r_iI0")
      if (allocated(error)) return
      call check(error, cav%wbranch(inew), tag1(SLOT_WBRANCH, iold), "wbranch")
      if (allocated(error)) return
      call check(error, cav%phi0(inew), tag1(SLOT_PHI0, iold), "phi0")
      if (allocated(error)) return
      call check(error, cav%cpjac_scal0(inew), tag1(SLOT_CPJAC, iold), "cpjac_scal0")
      if (allocated(error)) return
      call check(error, cav%w_f0(inew), tag1(SLOT_W_F0, iold), "w_f0")
      if (allocated(error)) return

      call check(error, cav%numbering(inew), itag(SLOT_NUMBERING, iold), "numbering")
      if (allocated(error)) return
      call check(error, cav%owner(inew), itag(SLOT_OWNER, iold), "owner")
      if (allocated(error)) return
      call check(error, cav%branch(inew), itag(SLOT_BRANCH, iold), "branch")
      if (allocated(error)) return
      call check(error, cav%anchor_id(inew), itag(SLOT_ANCHOR_ID, iold), "anchor_id")
      if (allocated(error)) return
      call check(error, cav%branch_count(inew), itag(SLOT_BRANCH_COUNT, iold), "branch_count")
      if (allocated(error)) return

      call check(error, cav%converged(inew) .eqv. (mod(iold, 3) == 0), "converged")

   end subroutine check_point

   !> Assert that every filtered array has exactly `nvalid` grid entries
   !>
   !> @param[out] error  Error handle
   !> @param[in]  cav    Filtered cavity
   !> @param[in]  nvalid Expected number of grid points
   subroutine check_grid_sizes(error, cav, nvalid)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Filtered cavity
      type(cavity_type_drop), intent(in) :: cav
      !> Expected number of grid points
      integer, intent(in) :: nvalid

      call check(error, size(cav%xyz, 1), 3, "xyz keeps its leading extent")
      if (allocated(error)) return
      call check(error, size(cav%xyz, 2), nvalid, "xyz")
      if (allocated(error)) return
      call check(error, size(cav%anchorxyz, 1), 3, "anchorxyz keeps its leading extent")
      if (allocated(error)) return
      call check(error, size(cav%anchorxyz, 2), nvalid, "anchorxyz")
      if (allocated(error)) return
      call check(error, size(cav%normal0, 1), 3, "normal0 keeps its leading extent")
      if (allocated(error)) return
      call check(error, size(cav%normal0, 2), nvalid, "normal0")
      if (allocated(error)) return

      call check(error, size(cav%wleb), nvalid, "wleb")
      if (allocated(error)) return
      call check(error, size(cav%anchor_wleb0), nvalid, "anchor_wleb0")
      if (allocated(error)) return
      call check(error, size(cav%lambda0), nvalid, "lambda0")
      if (allocated(error)) return
      call check(error, size(cav%iswig_f0), nvalid, "iswig_f0")
      if (allocated(error)) return
      call check(error, size(cav%f), nvalid, "f")
      if (allocated(error)) return
      call check(error, size(cav%anchor_xi0), nvalid, "anchor_xi0")
      if (allocated(error)) return
      call check(error, size(cav%rho), nvalid, "rho")
      if (allocated(error)) return
      call check(error, size(cav%r_iI0), nvalid, "r_iI0")
      if (allocated(error)) return
      call check(error, size(cav%wbranch), nvalid, "wbranch")
      if (allocated(error)) return
      call check(error, size(cav%phi0), nvalid, "phi0")
      if (allocated(error)) return
      call check(error, size(cav%cpjac_scal0), nvalid, "cpjac_scal0")
      if (allocated(error)) return
      call check(error, size(cav%w_f0), nvalid, "w_f0")
      if (allocated(error)) return

      call check(error, size(cav%numbering), nvalid, "numbering")
      if (allocated(error)) return
      call check(error, size(cav%owner), nvalid, "owner")
      if (allocated(error)) return
      call check(error, size(cav%branch), nvalid, "branch")
      if (allocated(error)) return
      call check(error, size(cav%anchor_id), nvalid, "anchor_id")
      if (allocated(error)) return
      call check(error, size(cav%branch_count), nvalid, "branch_count")
      if (allocated(error)) return
      call check(error, size(cav%converged), nvalid, "converged")

   end subroutine check_grid_sizes

   !> Tag carried by 1D real array `slot` at grid index `i`
   pure function tag1(slot, i) result(val)
      !> Array identity and grid index
      integer, intent(in) :: slot, i
      !> Tag value
      real(wp) :: val

      val = real(1000*slot + 10*i, wp)

   end function tag1

   !> Tag carried by component `j` of 2D real array `slot` at grid index `i`
   pure function tag2(slot, i, j) result(val)
      !> Array identity, grid index and component
      integer, intent(in) :: slot, i, j
      !> Tag value
      real(wp) :: val

      val = real(1000*slot + 10*i + j, wp)

   end function tag2

   !> Tag carried by integer array `slot` at grid index `i`
   pure function itag(slot, i) result(val)
      !> Array identity and grid index
      integer, intent(in) :: slot, i
      !> Tag value
      integer :: val

      val = 1000*slot + 10*i

   end function itag

   !> Softmax written straight from the definition, without the shift by the
   !> minimum that the implementation uses
   pure function ref_softmax(phi, sigma) result(weights)
      !> Objective values of one branch group
      real(wp), intent(in) :: phi(:)
      !> Softmax width
      real(wp), intent(in) :: sigma
      !> Reference weights
      real(wp) :: weights(size(phi))

      weights = exp(-phi/sigma)
      weights = weights/sum(weights)

   end function ref_softmax

end module test_cavity_drop_filter
