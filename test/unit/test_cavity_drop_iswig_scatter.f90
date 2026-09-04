!> Tests for the sparse iSwiG row and block scatter
!>
!> The scatters are pure index arithmetic over small arrays, so every check
!> here is exact rather than a finite difference: each one compares the scatter
!> against a dense reference built the obvious way -- fill a full `(3, nsph)` or
!> `(3, nsph, 3, nsph)` array with one contribution, then reduce -- and the
!> reference reduces in the *same* slot order the scatter writes in, so the two
!> agree to the bit even where a repeated atom id makes a column accumulate.
!>
!> The dense reference is what the production path may not do, not what it may
!> not be compared against: it exists only in this file.
!>
!> Two defects a fixture of distinct ids over a zeroed accumulator would miss
!> are covered deliberately:
!>
!> * an influence set may name the same atom twice -- from a repeated id, or
!>   because the owner also appears among the neighbours -- and the two
!>   contributions must sum, not overwrite;
!> * a pass reduces many grid points into one buffer, so a scatter that assigns
!>   would silently discard everything already there.
!>
!> Every synthetic fixture is guarded by an anti-vacuity floor, because a test
!> that runs a term without seeing it is worse than no test at all.
module test_cavity_drop_iswig_scatter
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use test_helpers, only: get_test_structures, get_test_radii
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   use moist_cavity_drop_gaussian, only: moist_cavity_drop_iswig, new_iswig, &
                                         iswig_workspace_type
   use moist_cavity_drop_derivatives_iswig_scatter, only: scatter_iswig_rows, &
                                                          scatter_iswig_block

   implicit none(type, external)
   private

   public :: collect_cavity_drop_iswig_scatter

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Atoms in the synthetic fixtures; large enough that `1` and `nsph` are
   !> genuinely distinct boundary ids and that most columns stay untouched
   integer, parameter :: nsph_fix = 7

   !> Neighbours in the synthetic sparse fixtures
   integer, parameter :: n_nb_fix = 4

   !> Influence-set size of the synthetic block fixtures, `n_nb + 1`
   integer, parameter :: n_blk_fix = 5

   !> Scalar the synthetic fixtures are weighted by
   !>
   !> Neither zero nor one, so a dropped or doubled weight cannot hide, and
   !> `-3/4` exactly, which is what keeps the comparison below exact; see
   !> [[fill_value]]
   real(wp), parameter :: weight_fix = -0.75_wp

   !> Smallest magnitude a fixture, or the change it makes, may have and still
   !> count as having been seen
   real(wp), parameter :: vacuity_floor = 1.0e-6_wp

   !> Agreement demanded of the two routes through the real second derivative:
   !> the sparse rows collapse the rank-one term into one scalar per direction
   !> while the block multiplies the stored outer product, so the two factor the
   !> same algebra differently and agree tightly rather than bitwise
   real(wp), parameter :: real_data_tol = 1.0e-12_wp

   !> Born width of the nleb = 110 Lebedev grid, as `parameters.f90:24` sets it
   real(wp), parameter :: iswig_swx = 4.900490_wp

   !> Representative largest raw Lebedev weight, as `setup.f90:35` passes it
   real(wp), parameter :: iswig_wleb_max = 1.25e-2_wp

   !> Gaussian width the real-data probes are evaluated at
   real(wp), parameter :: iswig_xi = 1.6_wp

   !> Structures drawn for the real-data probe; must be a multiple of five
   integer, parameter :: n_real_structures = 5

   !> Switching-function window a real-data probe must land in to be usable: a
   !> saturated point has zero rows and would scatter nothing
   real(wp), parameter :: f_window(2) = [0.02_wp, 0.98_wp]

   !> Usable probes the real-data test must find before it counts as having run
   integer, parameter :: min_real_probes = 20

contains

   !> Collect the suite
   subroutine collect_cavity_drop_iswig_scatter(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rows_scatter_matches_dense", test_rows_matches_dense), &
                  new_unittest("rows_scatter_accumulates_aliases", test_rows_aliased), &
                  new_unittest("block_scatter_matches_dense", test_block_matches_dense), &
                  new_unittest("block_scatter_accumulates_aliases", test_block_aliased), &
                  new_unittest("scatter_real_iswig_data", test_real_iswig_data) &
                  ]
   end subroutine collect_cavity_drop_iswig_scatter

   !* ================================================================================= *!
   !*                          Synthetic fixtures                                       *!
   !* ================================================================================= *!

   !> Deterministic fill whose whole arithmetic is exact in double precision
   !>
   !> Every value is an odd multiple of `2**-11` bounded by `512`, so it is never
   !> zero, and `weight_fix` is `-3/4`. Products are therefore multiples of
   !> `3 * 2**-13` bounded by `384`, and the longest sum any fixture forms -- one
   !> pre-fill plus four block contributions -- is a multiple of `2**-13`
   !> bounded by `2**11`. That needs 24 mantissa bits, so *every* intermediate
   !> on both sides of the comparison is exactly representable.
   !>
   !> This is not decoration. The scatter's `acc + weight*row` is contracted into
   !> a single FMA at the project's optimisation level while the dense reference
   !> rounds the product before adding it, so the two differ by an ulp on
   !> arbitrary data and the exactness the spec asks for is unreachable through
   !> a merely matching summation order. On dyadic data the FMA and the rounded
   !> product are the same number, and equality holds to the bit whatever the
   !> compiler contracts.
   !>
   !> A period of `1048573` over at most 225 slots keeps the values distinct, so
   !> a swapped index cannot cancel out of a comparison.
   !>
   !> @param[in]  seed  Fixture seed
   !> @param[in]  islot Linear slot index within the fixture
   !> @return     val   Fill value; dyadic, never zero, of magnitude up to 512
   pure function fill_value(seed, islot) result(val)
      !> Fixture seed
      integer, intent(in) :: seed
      !> Linear slot index
      integer, intent(in) :: islot
      !> Fill value
      real(wp) :: val

      !> Prime period of the generator, just under `2**20`
      integer(int64), parameter :: period = 1048573_int64
      !> Half the period, offset by a half so that no value lands on zero
      real(wp), parameter :: centre = 524286.5_wp
      !> Scale to a multiple of `2**-11` after the half-integer offset
      real(wp), parameter :: scale = 0.0009765625_wp

      integer(int64) :: state

      state = modulo(int(seed, int64)*7919_int64 + int(islot, int64)*104729_int64, period)
      val = (real(state, wp) - centre)*scale

   end function fill_value

   !> Fill a compact row fixture and the accumulator it is scattered into
   !>
   !> @param[in]  seed      Fixture seed
   !> @param[out] rows      Neighbour rows (3, n_nb_fix)
   !> @param[out] owner_row Owner row (3)
   !> @param[out] acc0      Nonzero pre-fill of the accumulator (3, nsph_fix)
   pure subroutine fill_rows_fixture(seed, rows, owner_row, acc0)
      !> Fixture seed
      integer, intent(in) :: seed
      !> Neighbour rows
      real(wp), intent(out) :: rows(ndim, n_nb_fix)
      !> Owner row
      real(wp), intent(out) :: owner_row(ndim)
      !> Nonzero pre-fill of the accumulator
      real(wp), intent(out) :: acc0(ndim, nsph_fix)

      integer :: jj, iaxis

      do jj = 1, n_nb_fix
         do iaxis = 1, ndim
            rows(iaxis, jj) = fill_value(seed, iaxis + ndim*jj)
         end do
      end do
      do iaxis = 1, ndim
         owner_row(iaxis) = fill_value(seed + 5, iaxis)
      end do
      do jj = 1, nsph_fix
         do iaxis = 1, ndim
            acc0(iaxis, jj) = fill_value(seed + 11, iaxis + ndim*jj)
         end do
      end do

   end subroutine fill_rows_fixture

   !> Fill a local-block fixture and the accumulator it is scattered into
   !>
   !> @param[in]  seed Fixture seed
   !> @param[out] blk  Local second-derivative block (3, n, 3, n)
   !> @param[out] acc0 Nonzero pre-fill of the accumulator (3, nsph, 3, nsph)
   pure subroutine fill_block_fixture(seed, blk, acc0)
      !> Fixture seed
      integer, intent(in) :: seed
      !> Local second-derivative block
      real(wp), intent(out) :: blk(ndim, n_blk_fix, ndim, n_blk_fix)
      !> Nonzero pre-fill of the accumulator
      real(wp), intent(out) :: acc0(ndim, nsph_fix, ndim, nsph_fix)

      integer :: ia, ib, iaxis, jaxis, islot

      islot = 0
      do ib = 1, n_blk_fix
         do jaxis = 1, ndim
            do ia = 1, n_blk_fix
               do iaxis = 1, ndim
                  islot = islot + 1
                  blk(iaxis, ia, jaxis, ib) = fill_value(seed, islot)
               end do
            end do
         end do
      end do

      islot = 0
      do ib = 1, nsph_fix
         do jaxis = 1, ndim
            do ia = 1, nsph_fix
               do iaxis = 1, ndim
                  islot = islot + 1
                  acc0(iaxis, ia, jaxis, ib) = fill_value(seed + 17, islot)
               end do
            end do
         end do
      end do

   end subroutine fill_block_fixture

   !* ================================================================================= *!
   !*                          Dense references                                         *!
   !* ================================================================================= *!

   !> Dense reference for the compact row scatter
   !>
   !> One full `(3, nsph)` array per compact slot, reduced in slot order. This
   !> is the route the production scatter exists to avoid, and reducing in the
   !> same order the scatter writes in is what makes the comparison exact even
   !> when a column receives several contributions.
   !>
   !> @param[in]  n_nb      Number of neighbours
   !> @param[in]  nb_idx    Neighbour atom ids
   !> @param[in]  rows      Neighbour rows
   !> @param[in]  owner     Owner atom id
   !> @param[in]  owner_row Owner row
   !> @param[in]  weight    Scalar the point is weighted by
   !> @param[in]  acc0      Pre-fill of the accumulator
   !> @param[out] ref       Reference accumulator (3, nsph)
   pure subroutine rows_reference(n_nb, nb_idx, rows, owner, owner_row, weight, acc0, ref)
      !> Number of neighbours
      integer, intent(in) :: n_nb
      !> Neighbour atom ids
      integer, intent(in) :: nb_idx(:)
      !> Neighbour rows
      real(wp), intent(in) :: rows(:, :)
      !> Owner atom id
      integer, intent(in) :: owner
      !> Owner row
      real(wp), intent(in) :: owner_row(ndim)
      !> Scalar the point is weighted by
      real(wp), intent(in) :: weight
      !> Pre-fill of the accumulator
      real(wp), intent(in) :: acc0(:, :)
      !> Reference accumulator
      real(wp), intent(out) :: ref(:, :)

      real(wp) :: dense(ndim, size(acc0, 2))
      integer :: jj

      ref(:, :) = acc0(:, :)
      do jj = 1, n_nb
         dense(:, :) = 0.0_wp
         dense(:, nb_idx(jj)) = weight*rows(:, jj)
         ref(:, :) = ref(:, :) + dense(:, :)
      end do

      dense(:, :) = 0.0_wp
      dense(:, owner) = weight*owner_row(:)
      ref(:, :) = ref(:, :) + dense(:, :)

   end subroutine rows_reference

   !> Dense reference for the influence-set block scatter
   !>
   !> One full `(3, nsph, 3, nsph)` array per `(ia, ib)` pair of the influence
   !> set, reduced with the column outermost so the per-element order matches
   !> the scatter's.
   !>
   !> @param[in]  n      Influence-set size
   !> @param[in]  idx    Atom ids of the influence set, owner first
   !> @param[in]  blk    Local second-derivative block
   !> @param[in]  weight Scalar the point is weighted by
   !> @param[in]  acc0   Pre-fill of the accumulator
   !> @param[out] ref    Reference accumulator (3, nsph, 3, nsph)
   pure subroutine block_reference(n, idx, blk, weight, acc0, ref)
      !> Influence-set size
      integer, intent(in) :: n
      !> Atom ids of the influence set, owner first
      integer, intent(in) :: idx(:)
      !> Local second-derivative block
      real(wp), intent(in) :: blk(:, :, :, :)
      !> Scalar the point is weighted by
      real(wp), intent(in) :: weight
      !> Pre-fill of the accumulator
      real(wp), intent(in) :: acc0(:, :, :, :)
      !> Reference accumulator
      real(wp), intent(out) :: ref(:, :, :, :)

      real(wp) :: dense(ndim, size(acc0, 2), ndim, size(acc0, 4))
      integer :: ia, ib

      ref(:, :, :, :) = acc0(:, :, :, :)
      do ib = 1, n
         do ia = 1, n
            dense(:, :, :, :) = 0.0_wp
            dense(:, idx(ia), :, idx(ib)) = weight*blk(:, ia, :, ib)
            ref(:, :, :, :) = ref(:, :, :, :) + dense(:, :, :, :)
         end do
      end do

   end subroutine block_reference

   !* ================================================================================= *!
   !*                          Compact row scatter                                      *!
   !* ================================================================================= *!

   !> Distinct ids, boundary ids and a pre-filled accumulator
   !>
   !> `nb_idx` spans atom `1` and atom `nsph` so an off-by-one in the scatter
   !> would run off an end, and the accumulator is pre-filled so that a scatter
   !> which assigned would lose the pre-fill.
   subroutine test_rows_matches_dense(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Distinct ids, spanning both boundaries; owner outside the neighbour set
      integer, parameter :: nb_idx(n_nb_fix) = [1, 3, 5, nsph_fix]
      integer, parameter :: owner = 2

      real(wp) :: rows(ndim, n_nb_fix), owner_row(ndim)
      real(wp) :: acc0(ndim, nsph_fix), acc(ndim, nsph_fix), ref(ndim, nsph_fix)

      call fill_rows_fixture(3, rows, owner_row, acc0)

      call check_alive(error, maxval(abs(rows)), "neighbour rows")
      if (allocated(error)) return
      call check_alive(error, maxval(abs(owner_row)), "owner row")
      if (allocated(error)) return
      call check_alive(error, maxval(abs(acc0)), "accumulator pre-fill")
      if (allocated(error)) return

      call rows_reference(n_nb_fix, nb_idx, rows, owner, owner_row, weight_fix, acc0, ref)
      call check_alive(error, maxval(abs(ref - acc0)), "row scatter contribution")
      if (allocated(error)) return

      acc(:, :) = acc0(:, :)
      call scatter_iswig_rows(n_nb_fix, nb_idx, rows, owner, owner_row, weight_fix, acc)

      call check_exact(error, maxval(abs(acc - ref)), "row scatter against dense reference")
      if (allocated(error)) return

      ! The pre-fill must still be in there: an assigning scatter would have
      ! replaced every touched column instead of adding to it.
      call check_exact(error, &
                       maxval(abs(acc(:, owner) - (acc0(:, owner) + weight_fix*owner_row))), &
                       "owner column adds to its pre-fill")
      if (allocated(error)) return
      call check_alive(error, maxval(abs(acc0(:, owner))), "owner column pre-fill")

   end subroutine test_rows_matches_dense

   !> Repeated ids, including the owner among the neighbours
   !>
   !> `nb_idx` names atom `6` twice and atom `6` is also the owner, so that
   !> column takes three separate contributions. A scatter written as one
   !> vector-subscripted section, or one that assigns, keeps only the last of
   !> them; the explicit three-term check below is what distinguishes the two.
   subroutine test_rows_aliased(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Atom 6 appears twice and is also the owner; both boundaries present
      integer, parameter :: nb_idx(n_nb_fix) = [6, 1, 6, nsph_fix]
      integer, parameter :: owner = 6

      real(wp) :: rows(ndim, n_nb_fix), owner_row(ndim)
      real(wp) :: acc0(ndim, nsph_fix), acc(ndim, nsph_fix), ref(ndim, nsph_fix)
      real(wp) :: expect(ndim)

      call fill_rows_fixture(41, rows, owner_row, acc0)

      call check_alive(error, maxval(abs(rows)), "neighbour rows")
      if (allocated(error)) return

      call rows_reference(n_nb_fix, nb_idx, rows, owner, owner_row, weight_fix, acc0, ref)
      call check_alive(error, maxval(abs(ref - acc0)), "row scatter contribution")
      if (allocated(error)) return

      acc(:, :) = acc0(:, :)
      call scatter_iswig_rows(n_nb_fix, nb_idx, rows, owner, owner_row, weight_fix, acc)

      call check_exact(error, maxval(abs(acc - ref)), "aliased row scatter against dense reference")
      if (allocated(error)) return

      ! Spelled out rather than left to the reference: the aliased column is the
      ! sum of all three contributions on top of its pre-fill, in slot order.
      expect(:) = acc0(:, owner)
      expect(:) = expect(:) + weight_fix*rows(:, 1)
      expect(:) = expect(:) + weight_fix*rows(:, 3)
      expect(:) = expect(:) + weight_fix*owner_row(:)
      call check_exact(error, maxval(abs(acc(:, owner) - expect)), "aliased column accumulates")
      if (allocated(error)) return

      ! ... and the sum is not any one of them, so overwriting cannot pass.
      call check_alive(error, maxval(abs(weight_fix*rows(:, 1))), "first aliased contribution")
      if (allocated(error)) return
      call check_alive(error, maxval(abs(weight_fix*rows(:, 3))), "second aliased contribution")
      if (allocated(error)) return
      call check_alive(error, maxval(abs(weight_fix*owner_row)), "owner contribution")

   end subroutine test_rows_aliased

   !* ================================================================================= *!
   !*                          Influence-set block scatter                              *!
   !* ================================================================================= *!

   !> Distinct ids, boundary ids and a pre-filled Hessian accumulator
   subroutine test_block_matches_dense(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Owner first, then distinct neighbours spanning both boundaries
      integer, parameter :: idx(n_blk_fix) = [4, 1, 3, nsph_fix, 5]

      real(wp) :: blk(ndim, n_blk_fix, ndim, n_blk_fix)
      real(wp) :: acc0(ndim, nsph_fix, ndim, nsph_fix)
      real(wp) :: acc(ndim, nsph_fix, ndim, nsph_fix)
      real(wp) :: ref(ndim, nsph_fix, ndim, nsph_fix)

      call fill_block_fixture(7, blk, acc0)

      call check_alive(error, maxval(abs(blk)), "local block")
      if (allocated(error)) return
      call check_alive(error, maxval(abs(acc0)), "accumulator pre-fill")
      if (allocated(error)) return

      call block_reference(n_blk_fix, idx, blk, weight_fix, acc0, ref)
      call check_alive(error, maxval(abs(ref - acc0)), "block scatter contribution")
      if (allocated(error)) return

      acc(:, :, :, :) = acc0(:, :, :, :)
      call scatter_iswig_block(n_blk_fix, idx, blk, weight_fix, acc)

      call check_exact(error, maxval(abs(acc - ref)), "block scatter against dense reference")
      if (allocated(error)) return

      ! Untouched atoms must be untouched to the bit: atoms 2 and 6 are in no
      ! slot of `idx`, so an overrun of the scatter would show here.
      call check_exact(error, maxval(abs(acc(:, 2, :, :) - acc0(:, 2, :, :))), &
                       "row of an atom outside the influence set")
      if (allocated(error)) return
      call check_exact(error, maxval(abs(acc(:, :, :, 6) - acc0(:, :, :, 6))), &
                       "column of an atom outside the influence set")

   end subroutine test_block_matches_dense

   !> Repeated ids in the influence set
   !>
   !> Atom `2` is the owner and appears again in slot 3, and atom `nsph` fills
   !> slots 2 and 5. The `(2, 2)` diagonal block therefore takes four separate
   !> `(ia, ib)` contributions and the `(nsph, nsph)` block another four; a
   !> scatter that overwrote would keep one of each.
   subroutine test_block_aliased(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Owner first; atom 2 and atom nsph each appear twice, atom 1 once
      integer, parameter :: idx(n_blk_fix) = [2, nsph_fix, 2, 1, nsph_fix]

      real(wp) :: blk(ndim, n_blk_fix, ndim, n_blk_fix)
      real(wp) :: acc0(ndim, nsph_fix, ndim, nsph_fix)
      real(wp) :: acc(ndim, nsph_fix, ndim, nsph_fix)
      real(wp) :: ref(ndim, nsph_fix, ndim, nsph_fix)
      real(wp) :: expect(ndim, ndim)

      call fill_block_fixture(23, blk, acc0)

      call check_alive(error, maxval(abs(blk)), "local block")
      if (allocated(error)) return

      call block_reference(n_blk_fix, idx, blk, weight_fix, acc0, ref)
      call check_alive(error, maxval(abs(ref - acc0)), "block scatter contribution")
      if (allocated(error)) return

      acc(:, :, :, :) = acc0(:, :, :, :)
      call scatter_iswig_block(n_blk_fix, idx, blk, weight_fix, acc)

      call check_exact(error, maxval(abs(acc - ref)), &
                       "aliased block scatter against dense reference")
      if (allocated(error)) return

      ! The aliased diagonal block, spelled out in the scatter's own order:
      ! column `ib` outermost, row `ia` innermost.
      expect(:, :) = acc0(:, 2, :, 2)
      expect(:, :) = expect(:, :) + weight_fix*blk(:, 1, :, 1)
      expect(:, :) = expect(:, :) + weight_fix*blk(:, 3, :, 1)
      expect(:, :) = expect(:, :) + weight_fix*blk(:, 1, :, 3)
      expect(:, :) = expect(:, :) + weight_fix*blk(:, 3, :, 3)
      call check_exact(error, maxval(abs(acc(:, 2, :, 2) - expect)), &
                       "aliased diagonal block accumulates all four pairs")
      if (allocated(error)) return

      ! The off-diagonal strip between the two aliased ids collects both of
      ! atom 2's slots against both of atom nsph's.
      expect(:, :) = acc0(:, 2, :, nsph_fix)
      expect(:, :) = expect(:, :) + weight_fix*blk(:, 1, :, 2)
      expect(:, :) = expect(:, :) + weight_fix*blk(:, 3, :, 2)
      expect(:, :) = expect(:, :) + weight_fix*blk(:, 1, :, 5)
      expect(:, :) = expect(:, :) + weight_fix*blk(:, 3, :, 5)
      call check_exact(error, maxval(abs(acc(:, 2, :, nsph_fix) - expect)), &
                       "aliased off-diagonal strip accumulates all four pairs")
      if (allocated(error)) return

      call check_alive(error, maxval(abs(weight_fix*blk(:, 1, :, 1))), "first aliased pair")

   end subroutine test_block_aliased

   !* ================================================================================= *!
   !*                          Genuine iSwiG output                                     *!
   !* ================================================================================= *!

   !> Scatter real second-derivative output of both producers and cross-check
   !>
   !> The synthetic fixtures test the index arithmetic; this one tests that the
   !> arithmetic is wired to the conventions the producers actually use. The
   !> block scatter fills a dense `(3, nsph, 3, nsph)` Hessian, which is then
   !> contracted with a nuclear direction and compared to the row scatter of
   !> [[iswig_swi_f2_rArB_sparse]] for the same direction. The two disagree the
   !> moment either scatter puts a slot on the wrong atom, because
   !> `swi2_rArB_block` is owner-first over `n_nb + 1` slots while
   !> `swi2_rArB_sparse` is neighbour-only over `n_nb` with the owner apart.
   subroutine test_real_iswig_data(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Weight applied to the point; neither zero nor one
      real(wp), parameter :: weight = 0.371_wp

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      type(iswig_workspace_type) :: work
      real(wp), allocatable :: radii(:), v(:, :)
      real(wp), allocatable :: rows2(:, :), blk(:, :, :, :), mix(:, :)
      real(wp), allocatable :: grad(:, :), hess(:, :, :, :), col(:, :)
      integer, allocatable :: idx(:)
      real(wp) :: pos(ndim), owner_row2(ndim), dxi2, d2xi, f_val, scale
      integer :: icase, owner, iatom, jatom, iaxis, n, nat, nprobe

      nprobe = 0
      call get_test_structures(mols, n_real_structures)

      do icase = 1, size(mols)
         mol = mols(icase)
         nat = mol%nat
         if (nat < 3) cycle
         call get_test_radii(mol, radii)

         call new_iswig(iswig, iswig_swx)
         call iswig%update(mol, radii, wleb_max=iswig_wleb_max)
         call work%init(iswig)

         if (allocated(rows2)) deallocate (rows2, blk, mix, idx, grad, hess, col, v)
         allocate (rows2(ndim, nat), blk(ndim, nat + 1, ndim, nat + 1))
         allocate (mix(ndim, nat + 1), idx(nat + 1))
         allocate (grad(ndim, nat), hess(ndim, nat, ndim, nat), col(ndim, nat))
         allocate (v(ndim, nat))

         do iatom = 1, nat
            do iaxis = 1, ndim
               v(iaxis, iatom) = sin(0.7_wp*real(11*icase + 13*iatom + 29*iaxis, wp))
            end do
         end do

         do owner = 1, nat
            call unsaturated_probe(iswig, mol, radii, owner, work, pos, f_val)
            if (f_val < f_window(1) .or. f_val > f_window(2)) cycle
            if (work%n_nb < 1) cycle
            nprobe = nprobe + 1

            ! Re-filled deliberately rather than relying on the scan having
            ! left the cache at `pos`, so the producers below are driven by a
            ! cache this loop can see being built.
            call iswig%swi_collect(pos, owner, iswig_xi, f_val, work)
            call iswig%swi2_rArB_block(work, n, idx, blk, mix, d2xi)
            hess(:, :, :, :) = 0.0_wp
            call scatter_iswig_block(n, idx, blk, weight, hess)

            call iswig%swi2_rArB_sparse(work, v, rows2, owner_row2, dxi2)
            grad(:, :) = 0.0_wp
            call scatter_iswig_rows(work%n_nb, work%idx, rows2, work%owner, &
                                    owner_row2, weight, grad)

            do iatom = 1, nat
               col(:, iatom) = 0.0_wp
               do jatom = 1, nat
                  col(:, iatom) = col(:, iatom) &
                                  + matmul(hess(:, iatom, :, jatom), v(:, jatom))
               end do
            end do

            scale = maxval(abs(grad))
            call check_alive(error, scale, "scattered real second-derivative rows")
            if (allocated(error)) return
            if (maxval(abs(col - grad)) > real_data_tol*scale) then
               call test_failed(error, "block scatter and row scatter disagree on real data")
               return
            end if
         end do

         call work%destroy()
      end do

      ! The draw is deterministic and yields 117 usable probes today; the floor
      ! is well under that but far above zero, so a fixture change that stopped
      ! the scan from landing in the switching window would be caught here
      ! rather than passing as a test that ran nothing.
      if (nprobe < min_real_probes) then
         call test_failed(error, "too few unsaturated iSwiG probes; the test saw almost nothing")
      end if

   end subroutine test_real_iswig_data

   !> Place a probe on the owner's sphere where the switching function is partial
   !>
   !> A point deep inside a neighbour has `f = 0` and every derivative row zero,
   !> and a point far outside every neighbour has `f = 1` with no neighbours
   !> cached; either would scatter nothing. A handful of fixed directions on the
   !> owner's sphere is scanned and the first partial one returned.
   !>
   !> @param[in]    iswig Switching function
   !> @param[in]    mol   Structure
   !> @param[in]    radii Atomic radii
   !> @param[in]    owner Owner atom index
   !> @param[inout] work  Caller-owned neighbour cache
   !> @param[out]   pos   Probe position
   !> @param[out]   f_val Switching value there; outside the window if none found
   subroutine unsaturated_probe(iswig, mol, radii, owner, work, pos, f_val)
      !> Switching function
      type(moist_cavity_drop_iswig), intent(in) :: iswig
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Atomic radii
      real(wp), intent(in) :: radii(:)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Caller-owned neighbour cache
      type(iswig_workspace_type), intent(inout) :: work
      !> Probe position
      real(wp), intent(out) :: pos(ndim)
      !> Switching value at the probe
      real(wp), intent(out) :: f_val

      !> Directions scanned on the owner's sphere, in spherical angles
      integer, parameter :: n_scan = 12

      integer :: iscan
      real(wp) :: theta, phi, dir(ndim)

      f_val = -1.0_wp
      pos(:) = mol%xyz(:, owner)

      do iscan = 1, n_scan
         theta = 0.31_wp + 0.41_wp*real(iscan, wp)
         phi = 0.17_wp + 1.13_wp*real(iscan, wp)
         dir(1) = sin(theta)*cos(phi)
         dir(2) = sin(theta)*sin(phi)
         dir(3) = cos(theta)
         pos(:) = mol%xyz(:, owner) + radii(owner)*dir(:)
         call iswig%swi_collect(pos, owner, iswig_xi, f_val, work)
         if (f_val > f_window(1) .and. f_val < f_window(2) .and. work%n_nb >= 1) return
      end do

   end subroutine unsaturated_probe

   !* ================================================================================= *!
   !*                          Assertion helpers                                        *!
   !* ================================================================================= *!

   !> Fail unless a deviation is exactly zero
   !>
   !> Every synthetic comparison here is a plain copy, a single multiply and an
   !> ordered sum on both sides, so equality holds to the bit and a tolerance
   !> would only hide a defect.
   !>
   !> @param[out] error Error handle
   !> @param[in]  dev   Observed deviation
   !> @param[in]  what  What was compared
   subroutine check_exact(error, dev, what)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Observed deviation
      real(wp), intent(in) :: dev
      !> What was compared
      character(len=*), intent(in) :: what

      ! Fixed length rather than a deferred-length allocatable: the suite runs
      ! its tests inside an `!$omp parallel do`, where gfortran has been seen to
      ! share the length temporary of a deferred-length result across threads.
      character(len=24) :: shown

      if (.not. dev <= 0.0_wp) then
         write (shown, '(es16.8)') dev
         call test_failed(error, "not exact: "//what, "deviation "//trim(adjustl(shown)))
      end if

   end subroutine check_exact

   !> Fail unless a magnitude is above the anti-vacuity floor
   !>
   !> @param[out] error Error handle
   !> @param[in]  val   Observed magnitude
   !> @param[in]  what  What was measured
   subroutine check_alive(error, val, what)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Observed magnitude
      real(wp), intent(in) :: val
      !> What was measured
      character(len=*), intent(in) :: what

      if (.not. val > vacuity_floor) call test_failed(error, "vacuous fixture: "//what)

   end subroutine check_alive

end module test_cavity_drop_iswig_scatter
