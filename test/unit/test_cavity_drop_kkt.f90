!> Test suite for the bordered KKT sensitivity solve of the DROP derivatives
!>
!> [[drop_kkt_factor_type]] factors
!>
!>   K = [ H_L  -g ]
!>       [ g^T   0 ]
!>
!> once per grid point and then solves batches of right-hand sides against the
!> stored factors. `solve` takes the primal batch, `solve_tangent` the tangent
!> one: differentiating `K x = b` gives `K dx = db - dK x`, with `dK` in the
!> same bordered layout, so the two `dlsf1_r` terms enter with opposite signs.
!>
!> The object is pure linear algebra, so these tests need no cavity, no level
!> set and no grid; `H0`, `g0` and the tangents below are plain numbers chosen
!> to make `K` well conditioned (`test_primal_well_conditioned` pins that).
!>
!> Three checks carry the weight:
!> * `test_tangent_explicit_dk` assembles `dK` by hand and forms `db - dK x`
!>   without the routine, which pins the relative sign of the two `dlsf1_r`
!>   terms algebraically -- the zero-tangent case cannot, since `dK` vanishes
!>   there and only the sign of `db` itself survives.
!> * `test_batched_matches_per_direction` re-solves one direction at a time and
!>   compares column for column, which localises a transposed or off-by-`nseed`
!>   batch column index that finite differences would only blur.
!> * the two finite-difference tests are the ground truth for the whole formula.
module test_cavity_drop_kkt
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, test_failed
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_factor_type
   use test_helpers, only: fd4_scalar, fd4_offsets
   implicit none(type, external)
   private

   public :: collect_cavity_drop_kkt

   !> Seeds and directions of the fixture; both larger than one, so the
   !> `(idir-1)*nseed + iseed` column convention is exercised rather than
   !> trivially satisfied
   integer, parameter :: nseed = 3, ndir = 3

   !> Tolerance for the exact algebraic comparisons (same arithmetic, same order)
   real(wp), parameter :: exact_tol = 1.0e-13_wp
   !> Tolerance for the finite-difference comparisons, absolute on components of
   !> order one; the project-wide `1e-10` target, held since 2026-09-03.
   !>
   !> The references use the 4-point `fd4_scalar` stencil. They were 2-point
   !> central differences until 2026-09-07, which could not honestly reach this
   !> bound: every stencil point refactorizes and re-solves, so the solver's own
   !> error is amplified by `1/2h` and *rises* as the step falls, while
   !> truncation falls as `h^2`. Measured worst deviation, 2-point:
   !>
   !>     h    1e-3    1e-4    3e-5    1e-5     1e-6
   !>     dev  1.5e-6  1.5e-8  1.4e-9  1.6e-10  9.1e-11
   !>
   !> Clean `h^2` down to `3e-5`, then the fall stops: `1e-6` gives `9.1e-11`
   !> where truncation alone predicts `1.5e-12`. That is the roundoff floor, so
   !> `tangent_fd_h1em6` was not passing on accuracy -- it sat 10% under the
   !> bound with nowhere better to go, and `tangent_fd_h1em5` failed at
   !> `1.62e-10`. The same fixture with the 4-point stencil:
   !>
   !>     h    1e-2    3e-3     1e-3     3e-4     1e-4     3e-5     1e-6
   !>     dev  2.8e-8  2.3e-10  2.9e-12  7.7e-13  1.8e-12  7.2e-12  1.2e-10
   !>
   !> `h^4` truncation (123x across a 3.33x step) reaches the floor four decades
   !> coarser, where the `1/h` amplification is four decades smaller. The
   !> minimum is `7.7e-13`, 130x under the bound.
   real(wp), parameter :: fd_tol = 1.0e-10_wp

contains

   !> Register the bordered KKT tests
   subroutine collect_cavity_drop_kkt(testsuite)
      !> Registered tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("primal_well_conditioned", test_primal_well_conditioned), &
                  new_unittest("zero_tangent_matches_solve", test_zero_tangent), &
                  new_unittest("tangent_explicit_dk", test_tangent_explicit_dk), &
                  new_unittest("batched_matches_per_direction", test_batched_per_direction), &
                  new_unittest("tangent_fd_h1em3", test_tangent_fd_h1em3), &
                  new_unittest("tangent_fd_h1em4", test_tangent_fd_h1em4), &
                  new_unittest("shape_mismatch_reports_error", test_shape_mismatch) &
                  ]
   end subroutine collect_cavity_drop_kkt

   !> Build the fixture: a symmetric Lagrangian Hessian, a level-set gradient,
   !> their tangents, the primal right-hand sides and the tangent ones
   !>
   !> The third `dH` direction is deliberately **not** symmetric. A physical
   !> `d(H_L)` is symmetric, so a row/column transposition inside the routine
   !> would be invisible with symmetric tangents alone; the bordered matrix is
   !> factored with a general `getrf`, so a non-symmetric direction costs the
   !> reference nothing and pins the orientation.
   !>
   !> `db` varies with the seed *and* the direction, so a column index that mixes
   !> the two up cannot pass by accident.
   !>
   !> @param[out] H0 Lagrangian Hessian
   !> @param[out] g0 Level-set gradient
   !> @param[out] dH Lagrangian Hessian tangents, one per direction
   !> @param[out] dg Level-set gradient tangents, one per direction
   !> @param[out] b0 Primal right-hand sides, one per seed
   !> @param[out] db Tangent right-hand sides, column `(idir-1)*nseed + iseed`
   subroutine kkt_fixture(H0, g0, dH, dg, b0, db)
      !> Lagrangian Hessian and its tangents
      real(wp), intent(out) :: H0(3, 3), dH(3, 3, ndir)
      !> Level-set gradient and its tangents
      real(wp), intent(out) :: g0(3), dg(3, ndir)
      !> Primal and tangent right-hand sides
      real(wp), intent(out) :: b0(4, nseed), db(4, nseed*ndir)

      !> Seed, direction and column indices
      integer :: iseed, idir, icol

      H0 = reshape([2.3_wp, 0.4_wp, -0.2_wp, &
                    0.4_wp, 1.7_wp, 0.35_wp, &
                    -0.2_wp, 0.35_wp, 3.1_wp], [3, 3])
      g0 = [0.6_wp, -1.3_wp, 0.9_wp]

      ! Symmetric, as a physical d(H_L) is
      dH(:, :, 1) = reshape([0.9_wp, -0.25_wp, 0.4_wp, &
                             -0.25_wp, 1.3_wp, 0.15_wp, &
                             0.4_wp, 0.15_wp, -0.7_wp], [3, 3])
      dH(:, :, 2) = reshape([-0.35_wp, 0.6_wp, 0.1_wp, &
                             0.6_wp, 0.2_wp, -0.45_wp, &
                             0.1_wp, -0.45_wp, 0.8_wp], [3, 3])
      ! Non-symmetric on purpose: pins the row-versus-column contraction
      dH(:, :, 3) = reshape([0.2_wp, 0.05_wp, 0.9_wp, &
                             1.1_wp, -0.6_wp, 0.25_wp, &
                             -0.3_wp, 0.7_wp, 0.45_wp], [3, 3])

      dg(:, 1) = [0.7_wp, 0.2_wp, -0.6_wp]
      dg(:, 2) = [-0.4_wp, 0.95_wp, 0.3_wp]
      dg(:, 3) = [0.15_wp, -0.55_wp, 1.05_wp]

      ! Three structurally different seeds, mirroring the shapes nuclear.f90
      ! builds: a pure constraint row, a pure position row, and a dense one
      b0(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp]
      b0(:, 2) = [0.85_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      b0(:, 3) = [0.3_wp, -0.7_wp, 1.1_wp, 0.25_wp]

      do idir = 1, ndir
         do iseed = 1, nseed
            icol = (idir - 1)*nseed + iseed
            db(1, icol) = 0.1_wp + 0.3_wp*(iseed - 1) - 0.2_wp*(idir - 1)
            db(2, icol) = -0.45_wp + 0.15_wp*(idir - 1) + 0.05_wp*(iseed - 1)
            db(3, icol) = 0.8_wp - 0.25_wp*(iseed - 1) + 0.1_wp*(idir - 1)
            db(4, icol) = 0.35_wp + 0.2_wp*(idir - 1) - 0.4_wp*(iseed - 1)
         end do
      end do
   end subroutine kkt_fixture

   !> Assemble the bordered matrix the way [[drop_kkt_factor]] does
   !>
   !> @param[in]  H Hessian block
   !> @param[in]  g Gradient block
   !> @returns    K Bordered matrix
   function bordered(H, g) result(K)
      !> Hessian block
      real(wp), intent(in) :: H(3, 3)
      !> Gradient block
      real(wp), intent(in) :: g(3)
      !> Bordered matrix
      real(wp) :: K(4, 4)

      K = 0.0_wp
      K(1:3, 1:3) = H
      K(1:3, 4) = -g
      K(4, 1:3) = g
   end function bordered

   !> Promote a failed solver call into a test failure
   !>
   !> @param[inout] error   Test error
   !> @param[in]    merr    Solver error, if any
   !> @param[in]    context Calling site, used to prefix the diagnostic
   !> @param[out]   failed  Whether the caller should return immediately
   subroutine promote(error, merr, context, failed)
      !> Test error
      type(error_type), allocatable, intent(inout) :: error
      !> Solver error
      type(mctc_error), allocatable, intent(in) :: merr
      !> Calling site
      character(len=*), intent(in) :: context
      !> Whether to return
      logical, intent(out) :: failed

      failed = allocated(merr)
      if (failed) call test_failed(error, context//": "//merr%message)
   end subroutine promote

   !> Compare two `(4, ncol)` batches component by component
   !>
   !> @param[inout] error    Test error
   !> @param[in]    actual   Batch under test
   !> @param[in]    expected Reference batch
   !> @param[in]    thr      Absolute tolerance
   !> @param[in]    label    Prefix identifying the comparison
   subroutine check_batch(error, actual, expected, thr, label)
      !> Test error
      type(error_type), allocatable, intent(inout) :: error
      !> Batch under test and its reference
      real(wp), intent(in) :: actual(:, :), expected(:, :)
      !> Absolute tolerance
      real(wp), intent(in) :: thr
      !> Comparison label
      character(len=*), intent(in) :: label

      !> Column and component indices
      integer :: icol, icomp

      do icol = 1, size(actual, 2)
         do icomp = 1, 4
            call check(error, actual(icomp, icol), expected(icomp, icol), thr=thr, &
                       message=label//": column "//to_string(icol)//" component "// &
                       to_string(icomp)//" is "//to_string(actual(icomp, icol))// &
                       ", expected "//to_string(expected(icomp, icol)))
            if (allocated(error)) return
         end do
      end do
   end subroutine check_batch

   !> The fixture must stay a test of a well-conditioned system
   !>
   !> Solving against the four unit vectors materialises `K^-1`, so both its
   !> largest entry and the residual of `K K^-1 = I` are pinned here. A future
   !> edit to `H0` or `g0` that drifts towards a degenerate projected point
   !> fails here rather than silently loosening every other test.
   !>
   !> @param[out] error Test error
   subroutine test_primal_well_conditioned(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      type(drop_kkt_factor_type) :: fac
      type(mctc_error), allocatable :: merr
      real(wp) :: H0(3, 3), dH(3, 3, ndir), g0(3), dg(3, ndir)
      real(wp) :: b0(4, nseed), db(4, nseed*ndir)
      real(wp) :: K(4, 4), Kinv(4, 4), prod(4, 4), rhs(4, nseed)
      logical :: failed
      integer :: i

      call kkt_fixture(H0, g0, dH, dg, b0, db)
      K = bordered(H0, g0)

      call fac%factor(H0, g0, "test_primal_well_conditioned", merr)
      call promote(error, merr, "factor", failed)
      if (failed) return

      Kinv = 0.0_wp
      do i = 1, 4
         Kinv(i, i) = 1.0_wp
      end do
      call fac%solve(Kinv, "test_primal_well_conditioned", merr)
      call promote(error, merr, "solve", failed)
      if (failed) return

      prod = matmul(K, Kinv)
      do i = 1, 4
         prod(i, i) = prod(i, i) - 1.0_wp
      end do
      call check(error, maxval(abs(prod)) < 1.0e-13_wp, &
                 "K*K^-1 deviates from the identity by "//to_string(maxval(abs(prod))))
      if (allocated(error)) return

      ! max|K| is 3.1 and max|K^-1| is ~0.52 here, so the system is benign
      call check(error, maxval(abs(Kinv)) < 2.0_wp, &
                 "Bordered inverse is larger than expected, max entry "// &
                 to_string(maxval(abs(Kinv)))//"; the fixture is no longer well conditioned")
      if (allocated(error)) return

      rhs = b0
      call fac%solve(rhs, "test_primal_well_conditioned", merr)
      call promote(error, merr, "solve", failed)
      if (failed) return

      call check_batch(error, matmul(K, rhs), b0, 1.0e-13_wp, "primal residual")
   end subroutine test_primal_well_conditioned

   !> A zero tangent must reduce to a plain solve of `db`
   !>
   !> With `dH = 0` and `dg = 0` the whole `dK x` correction drops out, so this
   !> pins that `db` enters the batch unnegated and that the factors are reused
   !> rather than rebuilt. It says nothing about the relative sign of the two
   !> `dlsf1_r` terms -- `test_tangent_explicit_dk` is what covers that.
   !>
   !> @param[out] error Test error
   subroutine test_zero_tangent(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      type(drop_kkt_factor_type) :: fac
      type(mctc_error), allocatable :: merr
      real(wp) :: H0(3, 3), dH(3, 3, ndir), g0(3), dg(3, ndir)
      real(wp) :: b0(4, nseed), db(4, nseed*ndir)
      real(wp) :: x(4, nseed), batch(4, nseed*ndir), reference(4, nseed*ndir)
      real(wp) :: dH_zero(3, 3, ndir), dg_zero(3, ndir)
      logical :: failed

      call kkt_fixture(H0, g0, dH, dg, b0, db)
      dH_zero = 0.0_wp
      dg_zero = 0.0_wp

      call fac%factor(H0, g0, "test_zero_tangent", merr)
      call promote(error, merr, "factor", failed)
      if (failed) return

      x = b0
      call fac%solve(x, "test_zero_tangent", merr)
      call promote(error, merr, "solve", failed)
      if (failed) return

      reference = db
      call fac%solve(reference, "test_zero_tangent", merr)
      call promote(error, merr, "solve", failed)
      if (failed) return

      batch = db
      call fac%solve_tangent(dH_zero, dg_zero, x, batch, "test_zero_tangent", merr)
      call promote(error, merr, "solve_tangent", failed)
      if (failed) return

      call check_batch(error, batch, reference, exact_tol, "zero tangent")
   end subroutine test_zero_tangent

   !> The tangent must match `K^-1 (db - dK x)` with `dK` assembled by hand
   !>
   !> The reference builds the bordered `dK` with the same layout
   !> [[drop_kkt_factor]] uses for `K`, contracts it with `x` as an ordinary
   !> 4x4 matrix-vector product and solves the result with the primal `solve`.
   !> Nothing of `solve_tangent`'s own arithmetic is reused, so a flipped sign
   !> on either `dlsf1_r` term, or a transposed `dH` contraction, shows up here
   !> exactly rather than as a finite-difference discrepancy.
   !>
   !> This test also pins that `x` survives the call unmodified.
   !>
   !> @param[out] error Test error
   subroutine test_tangent_explicit_dk(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      type(drop_kkt_factor_type) :: fac
      type(mctc_error), allocatable :: merr
      real(wp) :: H0(3, 3), dH(3, 3, ndir), g0(3), dg(3, ndir)
      real(wp) :: b0(4, nseed), db(4, nseed*ndir)
      real(wp) :: x(4, nseed), x_before(4, nseed)
      real(wp) :: batch(4, nseed*ndir), reference(4, nseed*ndir)
      real(wp) :: dK(4, 4)
      logical :: failed
      integer :: idir, iseed, icol

      call kkt_fixture(H0, g0, dH, dg, b0, db)

      call fac%factor(H0, g0, "test_tangent_explicit_dk", merr)
      call promote(error, merr, "factor", failed)
      if (failed) return

      x = b0
      call fac%solve(x, "test_tangent_explicit_dk", merr)
      call promote(error, merr, "solve", failed)
      if (failed) return
      x_before = x

      do idir = 1, ndir
         dK = bordered(dH(:, :, idir), dg(:, idir))
         do iseed = 1, nseed
            icol = (idir - 1)*nseed + iseed
            reference(:, icol) = db(:, icol) - matmul(dK, x(:, iseed))
         end do
      end do
      call fac%solve(reference, "test_tangent_explicit_dk", merr)
      call promote(error, merr, "solve", failed)
      if (failed) return

      batch = db
      call fac%solve_tangent(dH, dg, x, batch, "test_tangent_explicit_dk", merr)
      call promote(error, merr, "solve_tangent", failed)
      if (failed) return

      call check_batch(error, batch, reference, exact_tol, "explicit dK")
      if (allocated(error)) return

      call check_batch(error, x, x_before, 0.0_wp, "primal solution modified")
   end subroutine test_tangent_explicit_dk

   !> The batch must agree with a per-direction loop, column for column
   !>
   !> One `getrs` for `nseed*ndir` columns is the whole point of the batched
   !> signature, and a transposed or off-by-`nseed` column index is its
   !> characteristic failure. The reference calls the same routine with a single
   !> direction at a time -- where both a correct `(idir-1)*nseed + iseed` and a
   !> transposed `(iseed-1)*ndir + idir` collapse to `iseed` -- so a mixed-up
   !> batch index appears as a column permutation against it.
   !>
   !> @param[out] error Test error
   subroutine test_batched_per_direction(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      type(drop_kkt_factor_type) :: fac
      type(mctc_error), allocatable :: merr
      real(wp) :: H0(3, 3), dH(3, 3, ndir), g0(3), dg(3, ndir)
      real(wp) :: b0(4, nseed), db(4, nseed*ndir)
      real(wp) :: x(4, nseed)
      real(wp) :: batch(4, nseed*ndir), reference(4, nseed*ndir)
      logical :: failed
      integer :: idir, lo, hi

      call kkt_fixture(H0, g0, dH, dg, b0, db)

      call fac%factor(H0, g0, "test_batched_per_direction", merr)
      call promote(error, merr, "factor", failed)
      if (failed) return

      x = b0
      call fac%solve(x, "test_batched_per_direction", merr)
      call promote(error, merr, "solve", failed)
      if (failed) return

      batch = db
      call fac%solve_tangent(dH, dg, x, batch, "test_batched_per_direction", merr)
      call promote(error, merr, "solve_tangent", failed)
      if (failed) return

      reference = db
      do idir = 1, ndir
         lo = (idir - 1)*nseed + 1
         hi = idir*nseed
         call fac%solve_tangent(dH(:, :, idir:idir), dg(:, idir:idir), x, &
                                reference(:, lo:hi), "test_batched_per_direction", merr)
         call promote(error, merr, "solve_tangent", failed)
         if (failed) return
      end do

      call check_batch(error, batch, reference, exact_tol, "batched versus per-direction")
   end subroutine test_batched_per_direction

   !> Finite-difference check at `h = 1e-3`
   !> @param[out] error Test error
   subroutine test_tangent_fd_h1em3(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      call run_tangent_fd(error, 1.0e-3_wp)
   end subroutine test_tangent_fd_h1em3

   !> Finite-difference check at `h = 1e-4`
   !> @param[out] error Test error
   subroutine test_tangent_fd_h1em4(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      call run_tangent_fd(error, 1.0e-4_wp)
   end subroutine test_tangent_fd_h1em4

   !> Compare the tangent batch against central differences of the primal solve
   !>
   !> `K(t) = bordered(H0 + t dH, g0 + t dg)` and `b(t) = b0 + t db` are factored
   !> and solved afresh at each displacement, so the reference shares nothing
   !> with `solve_tangent` but the primal path.
   !>
   !> Both step sizes straddle the `eps^(1/3)` optimum of a central first
   !> difference, so neither is expected to be the better of the two; each is
   !> asserted against the same absolute tolerance and no monotone improvement
   !> is claimed.
   !>
   !> @param[inout] error Test error
   !> @param[in]    h     Displacement
   subroutine run_tangent_fd(error, h)
      !> Test error
      type(error_type), allocatable, intent(inout) :: error
      !> Displacement
      real(wp), intent(in) :: h

      type(drop_kkt_factor_type) :: fac, fac_step
      type(mctc_error), allocatable :: merr
      real(wp) :: H0(3, 3), dH(3, 3, ndir), g0(3), dg(3, ndir)
      real(wp) :: b0(4, nseed), db(4, nseed*ndir)
      real(wp) :: x(4, nseed)
      real(wp) :: batch(4, nseed*ndir), reference(4, nseed*ndir)
      real(wp) :: stencil(4, nseed*ndir, 4)
      logical :: failed
      integer :: idir, iseed, icol, ioff

      call kkt_fixture(H0, g0, dH, dg, b0, db)

      call fac%factor(H0, g0, "run_tangent_fd", merr)
      call promote(error, merr, "factor", failed)
      if (failed) return

      x = b0
      call fac%solve(x, "run_tangent_fd", merr)
      call promote(error, merr, "solve", failed)
      if (failed) return

      batch = db
      call fac%solve_tangent(dH, dg, x, batch, "run_tangent_fd", merr)
      call promote(error, merr, "solve_tangent", failed)
      if (failed) return

      ! The right-hand side, the Hessian and the gradient all move together at
      ! each stencil point, and the factorisation is redone there, so `1/h`
      ! amplifies the solver's own error -- see `fd_tol` for why that makes the
      ! 4-point stencil the only one that reaches the bound on this fixture
      do ioff = 1, size(fd4_offsets)
         associate (s_h => fd4_offsets(ioff)*h)
            do idir = 1, ndir
               do iseed = 1, nseed
                  icol = (idir - 1)*nseed + iseed
                  stencil(:, icol, ioff) = b0(:, iseed) + s_h*db(:, icol)
               end do

               call fac_step%factor(H0 + s_h*dH(:, :, idir), g0 + s_h*dg(:, idir), &
                                    "run_tangent_fd", merr)
               call promote(error, merr, "displaced factor", failed)
               if (failed) return
               call fac_step%solve( &
                  stencil(:, (idir - 1)*nseed + 1:idir*nseed, ioff), &
                  "run_tangent_fd", merr)
               call promote(error, merr, "displaced solve", failed)
               if (failed) return
            end do
         end associate
      end do

      reference = fd4_scalar(stencil(:, :, 1), stencil(:, :, 2), stencil(:, :, 3), &
                             stencil(:, :, 4), h)

      call check_batch(error, batch, reference, fd_tol, &
                       "4-point central difference at h = "//to_string(h))
   end subroutine run_tangent_fd

   !> An inconsistent batch shape must be reported rather than read out of bounds
   !>
   !> The `getrs` status of a factored system is unreachable from a correct
   !> caller, but the shape contract is not: `rhs` is the caller's buffer, sized
   !> once per thread, and `(4, nseed*ndir)` is the one thing `solve_tangent`
   !> cannot infer on its own.
   !>
   !> @param[out] error Test error
   subroutine test_shape_mismatch(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      type(drop_kkt_factor_type) :: fac
      type(mctc_error), allocatable :: merr
      real(wp) :: H0(3, 3), dH(3, 3, ndir), g0(3), dg(3, ndir)
      real(wp) :: b0(4, nseed), db(4, nseed*ndir)
      real(wp) :: x(4, nseed), short(4, nseed*ndir - 1)
      logical :: failed

      call kkt_fixture(H0, g0, dH, dg, b0, db)

      call fac%factor(H0, g0, "test_shape_mismatch", merr)
      call promote(error, merr, "factor", failed)
      if (failed) return

      x = b0
      call fac%solve(x, "test_shape_mismatch", merr)
      call promote(error, merr, "solve", failed)
      if (failed) return

      short = db(:, 1:nseed*ndir - 1)
      call fac%solve_tangent(dH, dg, x, short, "test_shape_mismatch", merr)
      call check(error, allocated(merr), &
                 "A batch of "//to_string(nseed*ndir - 1)//" columns for "// &
                 to_string(nseed)//" seeds and "//to_string(ndir)// &
                 " directions was accepted")
   end subroutine test_shape_mismatch

end module test_cavity_drop_kkt
