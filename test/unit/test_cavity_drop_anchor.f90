!> Tests for the DROP anchor grid, before any surface row exists
!>
!> - `fill`: Lebedev layout on every sphere
!> - `compute_switching`: iSwiG limits of isolated and buried anchors
!> - `filter`: the pre-filter compacts every anchor array with one permutation
module test_cavity_drop_anchor
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use moist_cavity_bornzeta, only: get_born_zeta
   use moist_cavity_drop_anchor, only: drop_anchor_type
   use moist_cavity_drop_gaussian, only: moist_cavity_drop_iswig, new_iswig
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_type, new_lebedev_grid
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none(type, external)
   private

   public :: collect_cavity_drop_anchor

   !> Array identities used by the tag helpers
   integer, parameter :: SLOT_XYZ = 1
   integer, parameter :: SLOT_XI0 = 7
   integer, parameter :: SLOT_OWNER = 15
   integer, parameter :: SLOT_ID = 17

   !> Lebedev points per sphere of the geometric fixtures
   integer, parameter :: NLEB = 302

   !> Radius of both spheres of the overlapping pair, bohr
   real(wp), parameter :: PAIR_RADIUS = 2.0_wp
   !> Distance of an anchor from the other sphere's surface beyond which it
   !> counts as clearly buried or clearly exposed, bohr
   real(wp), parameter :: PAIR_MARGIN = 0.5_wp

contains

   !> Collect all DROP anchor tests
   subroutine collect_cavity_drop_anchor(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fill_sphere_layout", test_fill_sphere_layout), &
                  new_unittest("fill_drops_stale_switching", test_fill_drops_stale_switching), &
                  new_unittest("switching_isolated_spheres", test_switching_isolated), &
                  new_unittest("switching_buried_and_exposed", test_switching_buried_exposed), &
                  new_unittest("prefilter_drops_buried_anchors", test_prefilter_drops_buried), &
                  new_unittest("filter_permutation_consistency", test_filter_permutation), &
                  new_unittest("filter_compacts_every_array", test_filter_size_audit), &
                  new_unittest("filter_cutoff_boundary", test_filter_cutoff_boundary), &
                  new_unittest("filter_keeps_all", test_filter_keeps_all), &
                  new_unittest("filter_keeps_none", test_filter_keeps_none), &
                  new_unittest("filter_skips_unallocated", test_filter_skips_unallocated), &
                  new_unittest("filter_reports_removed_count", test_filter_removed_count) &
                  ]

   end subroutine collect_cavity_drop_anchor

   !* ================================================================================= *!
   !*                               drop_anchor_type%fill                               *!
   !* ================================================================================= *!

   !> Anchor `(i - 1)*num_leb + j` is Lebedev node `j` on sphere `i`
   !>
   !> - lies on sphere `i` in the direction of node `j`
   !> - carries the node's weight, owner `i` and its own index as id
   !> - the weights of every sphere sum to 4 pi
   subroutine test_fill_sphere_layout(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Two spheres of different radius away from the origin
      real(wp), parameter :: centers(3, 2) = reshape([0.3_wp, -0.2_wp, 0.1_wp, &
                                                      1.0_wp, -2.0_wp, 3.5_wp], [3, 2])
      real(wp), parameter :: radii(2) = [2.5_wp, 1.75_wp]
      real(wp), parameter :: thr = 1.0e-12_wp

      type(drop_anchor_type) :: anchors
      type(moist_math_grid_angular_type) :: leb
      real(wp) :: offset(3)
      integer :: i, jj, ii

      call lebedev(error, leb)
      if (allocated(error)) return

      call anchors%fill(centers, radii, leb%points, leb%weights)

      call check(error, anchors%ngrid, 2*NLEB, "one anchor per node and sphere")
      if (allocated(error)) return
      call check_anchor_sizes(error, anchors, 2*NLEB, with_switching=.false.)
      if (allocated(error)) return

      do i = 1, 2
         do jj = 1, NLEB
            ii = (i - 1)*NLEB + jj
            call check(error, anchors%id(ii), ii, "id is the anchor index")
            if (allocated(error)) return
            call check(error, anchors%owner(ii), i, "owner is the sphere")
            if (allocated(error)) return
            call check(error, anchors%wleb0(ii), leb%weights(jj), "weight of the node")
            if (allocated(error)) return

            offset = anchors%xyz(:, ii) - centers(:, i)
            call check(error, norm2(offset), radii(i), thr_abs=thr*radii(i), thr_rel=10.0_wp*thr, &
                       message="anchor lies on its sphere")
            if (allocated(error)) return
            call check(error, dot_product(offset, leb%points(:, jj))/radii(i), 1.0_wp, thr_abs=thr, thr_rel=10.0_wp*thr, &
                       message="anchor lies in the direction of its node")
            if (allocated(error)) return
         end do

         call check(error, sum(anchors%wleb0((i - 1)*NLEB + 1:i*NLEB)), 4.0_wp*pi, thr_abs=thr, thr_rel=10.0_wp*thr, &
                    message="solid-angle weights of one sphere")
         if (allocated(error)) return
      end do

   end subroutine test_fill_sphere_layout

   !> A refill drops widths and switching values of the previous grid, so a
   !> stale `f0` of another length never reaches the pre-filter
   subroutine test_fill_drops_stale_switching(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: centers(3, 2) = reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                                      0.0_wp, 0.0_wp, 3.0_wp], [3, 2])
      real(wp), parameter :: radii(2) = [2.0_wp, 2.0_wp]

      type(drop_anchor_type) :: anchors
      type(moist_math_grid_angular_type) :: leb

      call lebedev(error, leb)
      if (allocated(error)) return

      call anchors%fill(centers, radii, leb%points, leb%weights)
      allocate (anchors%xi0(anchors%ngrid), source=1.0_wp)
      allocate (anchors%f0(anchors%ngrid), source=1.0_wp)

      call anchors%fill(centers(:, 1:1), radii(1:1), leb%points, leb%weights)

      call check(error, anchors%ngrid, NLEB, "one sphere after the refill")
      if (allocated(error)) return
      call check(error, .not. allocated(anchors%xi0), "stale widths dropped")
      if (allocated(error)) return
      call check(error, .not. allocated(anchors%f0), "stale switching values dropped")
      if (allocated(error)) return
      call check_anchor_sizes(error, anchors, NLEB, with_switching=.false.)

   end subroutine test_fill_drops_stale_switching

   !* ================================================================================= *!
   !*                         drop_anchor_type%compute_switching                        *!
   !* ================================================================================= *!

   !> Spheres far beyond each other's reach switch nothing off: every anchor
   !> keeps `f0 = 1`, also showing that no sphere switches its own anchors
   subroutine test_switching_isolated(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_anchor_type) :: anchors
      type(structure_type) :: mol

      call new(mol, [6, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 100.0_wp], [3, 2]))
      call switched_anchors(error, anchors, mol, [2.5_wp, 1.75_wp])
      if (allocated(error)) return

      call check_anchor_sizes(error, anchors, 2*NLEB, with_switching=.true.)
      if (allocated(error)) return
      call check(error, all(anchors%f0 == 1.0_wp), "isolated anchors are not switched")
      if (allocated(error)) return
      call check(error, all(anchors%xi0 > 0.0_wp .and. anchors%xi0 < huge(1.0_wp)), &
                 "finite positive widths")

   end subroutine test_switching_isolated

   !> In an overlapping pair, an anchor well inside the other sphere is switched
   !> off and one well outside it is untouched; every value lies in [0, 1]
   subroutine test_switching_buried_exposed(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_anchor_type) :: anchors
      type(structure_type) :: mol
      real(wp) :: r_other
      integer :: ia, n_buried, n_exposed

      call overlapping_pair(mol)
      call switched_anchors(error, anchors, mol, [PAIR_RADIUS, PAIR_RADIUS])
      if (allocated(error)) return

      call check(error, all(anchors%f0 >= 0.0_wp .and. anchors%f0 <= 1.0_wp), "f0 within [0, 1]")
      if (allocated(error)) return

      n_buried = 0
      n_exposed = 0
      do ia = 1, anchors%ngrid
         r_other = norm2(anchors%xyz(:, ia) - mol%xyz(:, 3 - anchors%owner(ia)))
         if (r_other < PAIR_RADIUS - PAIR_MARGIN) then
            n_buried = n_buried + 1
            call check(error, anchors%f0(ia) < 1.0e-12_wp, "buried anchor switched off")
            if (allocated(error)) return
         else if (r_other > PAIR_RADIUS + PAIR_MARGIN) then
            n_exposed = n_exposed + 1
            call check(error, anchors%f0(ia) > 1.0_wp - 1.0e-12_wp, "exposed anchor untouched")
            if (allocated(error)) return
         end if
      end do

      ! The fixture must probe both limits
      call check(error, n_buried > 0, "fixture has buried anchors")
      if (allocated(error)) return
      call check(error, n_exposed > 0, "fixture has exposed anchors")

   end subroutine test_switching_buried_exposed

   !> Fill, switching and pre-filter in production order: buried anchors are
   !> gone, every exposed anchor survives, and the survivors keep their order
   subroutine test_prefilter_drops_buried(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_anchor_type) :: anchors
      type(structure_type) :: mol
      real(wp) :: r_other
      integer :: ia, nremoved, n_exposed_before, n_exposed_after
      logical, allocatable :: exposed(:)

      call overlapping_pair(mol)
      call switched_anchors(error, anchors, mol, [PAIR_RADIUS, PAIR_RADIUS])
      if (allocated(error)) return

      allocate (exposed(anchors%ngrid))
      do ia = 1, anchors%ngrid
         r_other = norm2(anchors%xyz(:, ia) - mol%xyz(:, 3 - anchors%owner(ia)))
         exposed(ia) = r_other > PAIR_RADIUS + PAIR_MARGIN
      end do
      n_exposed_before = count(exposed)

      call anchors%filter(1.0e-10_wp, nremoved)

      call check(error, nremoved > 0, "buried anchors removed")
      if (allocated(error)) return
      call check(error, all(anchors%id(2:anchors%ngrid) > anchors%id(1:anchors%ngrid - 1)), &
                 "survivors keep their order")
      if (allocated(error)) return

      n_exposed_after = 0
      do ia = 1, anchors%ngrid
         r_other = norm2(anchors%xyz(:, ia) - mol%xyz(:, 3 - anchors%owner(ia)))
         call check(error, r_other >= PAIR_RADIUS - PAIR_MARGIN, "no buried anchor survives")
         if (allocated(error)) return
         if (exposed(anchors%id(ia))) n_exposed_after = n_exposed_after + 1
      end do
      call check(error, n_exposed_after, n_exposed_before, "every exposed anchor survives")

   end subroutine test_prefilter_drops_buried

   !* ================================================================================= *!
   !*                              drop_anchor_type%filter                              *!
   !* ================================================================================= *!

   !> A scattered keep mask must compact every array with the same permutation:
   !> each surviving anchor still carries all of its own field values
   subroutine test_filter_permutation(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_anchor_type) :: anchors
      integer, parameter :: n = 9
      !> Anchors 1, 4, 5 and 9 fall below the cutoff: leading, trailing and a run
      integer, parameter :: survivor(5) = [2, 3, 6, 7, 8]
      real(wp) :: wleb0(n), f0(n)
      integer :: i, k, nremoved

      do i = 1, n
         wleb0(i) = real(i, wp)
         f0(i) = 1.0_wp + 0.01_wp*real(i, wp)
      end do
      !> Kept anchors sit at wleb0*f0 >= 2.02, dropped ones at <= 9.0e-6
      f0([1, 4, 5, 9]) = 1.0e-6_wp

      call build_anchors(anchors, n, wleb0, f0)
      call anchors%filter(0.5_wp, nremoved)

      call check(error, anchors%ngrid, size(survivor), "surviving anchor count")
      if (allocated(error)) return
      call check(error, nremoved, n - size(survivor), "removed anchor count")
      if (allocated(error)) return

      do k = 1, size(survivor)
         i = survivor(k)
         call check_anchor(error, anchors, k, i, wleb0(i), f0(i))
         if (allocated(error)) return
      end do

   end subroutine test_filter_permutation

   !> Every anchor array comes back at the new length
   subroutine test_filter_size_audit(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_anchor_type) :: anchors
      integer, parameter :: n = 6
      integer, parameter :: nkept = 3
      real(wp) :: wleb0(n), f0(n)
      integer :: nremoved

      wleb0 = 1.0_wp
      f0 = 1.0_wp
      !> Keep anchors 2, 3 and 5
      f0([1, 4, 6]) = 0.0_wp

      call build_anchors(anchors, n, wleb0, f0)
      call anchors%filter(0.5_wp, nremoved)

      call check(error, anchors%ngrid, nkept, "ngrid follows the keep count")
      if (allocated(error)) return

      call check_anchor_sizes(error, anchors, nkept)

   end subroutine test_filter_size_audit

   !> The keep criterion is `wleb0*f0 > wleb_cut`, strictly and on the product
   subroutine test_filter_cutoff_boundary(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_anchor_type) :: anchors
      integer, parameter :: n = 4
      real(wp) :: wleb0(n), f0(n)
      integer :: nremoved

      !> 1: product exactly at the cutoff, 2: just above, 3: just below,
      !> 4: wleb0 far above the cutoff but killed by a tiny switching value
      wleb0 = [0.5_wp, 0.5_wp, 0.5_wp, 100.0_wp]
      f0 = [0.5_wp, 0.5000001_wp, 0.4999999_wp, 1.0e-6_wp]

      call build_anchors(anchors, n, wleb0, f0)
      !> 0.5*0.5 is exact in binary, so anchor 1 sits precisely on the cutoff
      call anchors%filter(0.25_wp, nremoved)

      call check(error, anchors%ngrid, 1, "only the anchor strictly above survives")
      if (allocated(error)) return

      call check_anchor(error, anchors, 1, 2, wleb0(2), f0(2))

   end subroutine test_filter_cutoff_boundary

   !> Keeping everything leaves the anchors untouched
   subroutine test_filter_keeps_all(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_anchor_type) :: anchors
      integer, parameter :: n = 5
      real(wp) :: wleb0(n), f0(n)
      integer :: i, nremoved

      wleb0 = 1.0_wp
      f0 = 1.0_wp

      call build_anchors(anchors, n, wleb0, f0)
      call anchors%filter(1.0e-12_wp, nremoved)

      call check(error, anchors%ngrid, n, "no anchor removed")
      if (allocated(error)) return
      call check(error, nremoved, 0, "nothing reported removed")
      if (allocated(error)) return

      call check_anchor_sizes(error, anchors, n)
      if (allocated(error)) return

      do i = 1, n
         call check_anchor(error, anchors, i, i, wleb0(i), f0(i))
         if (allocated(error)) return
      end do

   end subroutine test_filter_keeps_all

   !> An empty result is a valid one: zero-length arrays, ngrid zero, no crash,
   !> and a second pass over the empty set still behaves
   subroutine test_filter_keeps_none(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_anchor_type) :: anchors
      integer, parameter :: n = 4
      real(wp) :: wleb0(n), f0(n)
      integer :: nremoved

      wleb0 = 1.0_wp
      f0 = 1.0e-30_wp

      call build_anchors(anchors, n, wleb0, f0)
      call anchors%filter(1.0e-12_wp, nremoved)

      call check(error, anchors%ngrid, 0, "every anchor removed")
      if (allocated(error)) return
      call check(error, nremoved, n, "every anchor reported removed")
      if (allocated(error)) return

      call check_anchor_sizes(error, anchors, 0)
      if (allocated(error)) return

      call anchors%filter(1.0e-12_wp, nremoved)

      call check(error, anchors%ngrid, 0, "empty set stays empty")
      if (allocated(error)) return
      call check(error, nremoved, 0, "nothing left to remove")

   end subroutine test_filter_keeps_none

   !> Optional arrays that are not allocated are skipped rather than faulting,
   !> and stay unallocated
   subroutine test_filter_skips_unallocated(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_anchor_type) :: anchors
      integer, parameter :: n = 4
      integer :: i, nremoved

      anchors%ngrid = n

      !> The keep mask inputs plus positions and owners; no widths, no ids
      allocate (anchors%xyz(3, n), anchors%wleb0(n), anchors%f0(n), anchors%owner(n))
      do i = 1, n
         anchors%xyz(:, i) = [tag2(SLOT_XYZ, i, 1), tag2(SLOT_XYZ, i, 2), tag2(SLOT_XYZ, i, 3)]
         anchors%owner(i) = itag(SLOT_OWNER, i)
      end do
      anchors%wleb0 = 1.0_wp
      anchors%f0 = [1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]

      call anchors%filter(0.5_wp, nremoved)

      call check(error, anchors%ngrid, 2, "surviving anchor count")
      if (allocated(error)) return

      call check(error, size(anchors%xyz, 2), 2, "xyz compacted")
      if (allocated(error)) return
      call check(error, anchors%xyz(1, 2), tag2(SLOT_XYZ, 3, 1), "second survivor is anchor 3")
      if (allocated(error)) return
      call check(error, anchors%owner(2), itag(SLOT_OWNER, 3), "owner follows the same permutation")
      if (allocated(error)) return

      call check(error, .not. allocated(anchors%xi0), "xi0 stays unallocated")
      if (allocated(error)) return
      call check(error, .not. allocated(anchors%id), "id stays unallocated")

   end subroutine test_filter_skips_unallocated

   !> The removed count and the surviving count add up to the anchors before the filter
   subroutine test_filter_removed_count(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_anchor_type) :: anchors
      integer, parameter :: n = 7
      real(wp) :: wleb0(n), f0(n)
      integer :: nremoved

      wleb0 = 1.0_wp
      f0 = [1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 1.0_wp]

      call build_anchors(anchors, n, wleb0, f0)
      call anchors%filter(0.5_wp, nremoved)

      call check(error, nremoved, 3, "removed anchor count")
      if (allocated(error)) return
      call check(error, nremoved + anchors%ngrid, n, "removed and kept add up")

   end subroutine test_filter_removed_count

   !* ================================================================================= *!
   !*                                      Fixture                                      *!
   !* ================================================================================= *!

   !> Lebedev grid of the geometric fixtures
   !>
   !> @param[out] error Error handle
   !> @param[out] leb   Unit-sphere grid with `NLEB` nodes
   subroutine lebedev(error, leb)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Unit-sphere grid
      type(moist_math_grid_angular_type), intent(out) :: leb

      type(moist_error_type), allocatable :: grid_error

      call new_lebedev_grid(leb, grid_error, npts=NLEB, positive_weights_only=.true.)
      if (allocated(grid_error)) call test_failed(error, grid_error%message)

   end subroutine lebedev

   !> Two equal spheres whose centers are one radius apart
   !>
   !> @param[out] mol Carbon pair along z
   subroutine overlapping_pair(mol)
      !> Carbon pair along z
      type(structure_type), intent(out) :: mol

      call new(mol, [6, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, PAIR_RADIUS], [3, 2]))

   end subroutine overlapping_pair

   !> Anchors of `mol` with widths and switching values, as the cavity sets them up
   !>
   !> @param[out] error   Error handle
   !> @param[out] anchors Filled and switched anchor grid
   !> @param[in]  mol     Sphere centers
   !> @param[in]  radii   Sphere radii, bohr
   subroutine switched_anchors(error, anchors, mol, radii)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Filled and switched anchor grid
      type(drop_anchor_type), intent(out) :: anchors
      !> Sphere centers
      type(structure_type), intent(in) :: mol
      !> Sphere radii, bohr
      real(wp), intent(in) :: radii(:)

      type(moist_math_grid_angular_type) :: leb
      type(moist_cavity_drop_iswig) :: iswig
      real(wp) :: zeta
      logical :: found

      call lebedev(error, leb)
      if (allocated(error)) return
      call get_born_zeta(NLEB, zeta, found)
      if (.not. found) then
         call test_failed(error, "no fitted Born zeta for the fixture grid")
         return
      end if

      call anchors%fill(mol%xyz, radii, leb%points, leb%weights)
      call new_iswig(iswig, zeta)
      call iswig%update(mol, radii, wleb_max=maxval(anchors%wleb0), nthreads=1)
      call anchors%compute_switching(iswig, 1)

   end subroutine switched_anchors

   !> Build an anchor set carrying every array the pre-filter compacts
   !>
   !> All arrays except `wleb0` and `f0` are tagged with `tag*`, so a surviving
   !> anchor can be traced back to its original index field by field
   !>
   !> @param[out] anchors Anchor set to populate
   !> @param[in]  n       Number of anchors
   !> @param[in]  wleb0   Quadrature weights driving the keep mask
   !> @param[in]  f0      Switching values driving the keep mask
   subroutine build_anchors(anchors, n, wleb0, f0)
      !> Anchor set to populate
      type(drop_anchor_type), intent(out) :: anchors
      !> Number of anchors
      integer, intent(in) :: n
      !> Quadrature weights and switching values
      real(wp), intent(in) :: wleb0(:), f0(:)

      integer :: i, j

      anchors%ngrid = n

      allocate (anchors%xyz(3, n), anchors%wleb0(n), anchors%f0(n))
      allocate (anchors%xi0(n), anchors%owner(n), anchors%id(n))

      do i = 1, n
         do j = 1, 3
            anchors%xyz(j, i) = tag2(SLOT_XYZ, i, j)
         end do
         anchors%wleb0(i) = wleb0(i)
         anchors%f0(i) = f0(i)
         anchors%xi0(i) = tag1(SLOT_XI0, i)
         anchors%owner(i) = itag(SLOT_OWNER, i)
         anchors%id(i) = itag(SLOT_ID, i)
      end do

   end subroutine build_anchors

   !> Assert that anchor `inew` carries every field of original anchor `iold`
   !>
   !> @param[out] error     Error handle
   !> @param[in]  anchors   Filtered anchor set
   !> @param[in]  inew      Index after filtering
   !> @param[in]  iold      Index before filtering
   !> @param[in]  wleb0_ref Quadrature weight the anchor started with
   !> @param[in]  f0_ref    Switching value the anchor started with
   subroutine check_anchor(error, anchors, inew, iold, wleb0_ref, f0_ref)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Filtered anchor set
      type(drop_anchor_type), intent(in) :: anchors
      !> Compacted and original index
      integer, intent(in) :: inew, iold
      !> Untagged values the anchor started with
      real(wp), intent(in) :: wleb0_ref, f0_ref

      integer :: j

      do j = 1, 3
         call check(error, anchors%xyz(j, inew), tag2(SLOT_XYZ, iold, j), "xyz")
         if (allocated(error)) return
      end do
      call check(error, anchors%wleb0(inew), wleb0_ref, "wleb0")
      if (allocated(error)) return
      call check(error, anchors%f0(inew), f0_ref, "f0")
      if (allocated(error)) return
      call check(error, anchors%xi0(inew), tag1(SLOT_XI0, iold), "xi0")
      if (allocated(error)) return
      call check(error, anchors%owner(inew), itag(SLOT_OWNER, iold), "owner")
      if (allocated(error)) return
      call check(error, anchors%id(inew), itag(SLOT_ID, iold), "id")

   end subroutine check_anchor

   !> Assert that every anchor array has exactly `nvalid` entries
   !>
   !> @param[out] error          Error handle
   !> @param[in]  anchors        Anchor set
   !> @param[in]  nvalid         Expected number of anchors
   !> @param[in]  with_switching Whether `xi0` and `f0` are checked, default true
   subroutine check_anchor_sizes(error, anchors, nvalid, with_switching)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Anchor set
      type(drop_anchor_type), intent(in) :: anchors
      !> Expected number of anchors
      integer, intent(in) :: nvalid
      !> Whether `xi0` and `f0` are checked
      logical, intent(in), optional :: with_switching

      logical :: switching

      switching = .true.
      if (present(with_switching)) switching = with_switching

      call check(error, size(anchors%xyz, 1), 3, "xyz keeps its leading extent")
      if (allocated(error)) return
      call check(error, size(anchors%xyz, 2), nvalid, "xyz")
      if (allocated(error)) return
      call check(error, size(anchors%wleb0), nvalid, "wleb0")
      if (allocated(error)) return
      call check(error, size(anchors%owner), nvalid, "owner")
      if (allocated(error)) return
      call check(error, size(anchors%id), nvalid, "id")
      if (allocated(error) .or. .not. switching) return
      call check(error, size(anchors%f0), nvalid, "f0")
      if (allocated(error)) return
      call check(error, size(anchors%xi0), nvalid, "xi0")

   end subroutine check_anchor_sizes

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

end module test_cavity_drop_anchor
