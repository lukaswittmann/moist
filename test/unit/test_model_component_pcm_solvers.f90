!> Unit tests for the PCM operator and solver types
!>
!> - pcm_dense_solver: factorization chosen per operator (Cholesky, Bunch-Kaufman,
!>   LU), normal and transpose solves, block solves, factor reuse and failures
!> - pcm_gmres_solver: real and complex operators against the dense solver in
!>   both orientations, warm start, zero right-hand side, non-convergence,
!>   restarts, statistics, rebuilt preconditioner on setup, invalid calls
!> - pcm_block_jacobi: block size 1 is Jacobi, clusters of a two-field block
!>   operator by coordinate bisection and by contiguous ranges
!> - sytrf/sytrs wrappers of moist_math_lapack on a symmetric indefinite system
!> - Complex getrf/getrs/sytrf/sytrs wrappers (plain transpose, complex symmetric)
!> - setup of an id that holds factors (repeated medium phase) factors anew
!> - solve_complex: real operators with complex right-hand sides (dense, GMRES
!>   and the pcm_solver default), complex operators through zsytrf/zgetrf and
!>   complex GMRES, the default of a non-overriding extension
!> - pcm_nonlinear_problem: Newton with the Armijo line search on the cubic
!>   x^3 - 2 x + 2, converging from x = -2 and stalling from x = 0 at the
!>   local minimum of |F| (singular Jacobian, root out of reach), iteration limit
module test_model_component_pcm_solvers
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, sp
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_env_error, only: moist_error_type => error_type
   use test_helpers, only: check_moist_error
   use moist_math_lapack, only: getrf, getrs, sytrf, sytrs
   use moist_model_continuum_component_pcm_solvers, only: pcm_operator, pcm_complex_operator, &
      & pcm_solver, pcm_dense_solver, pcm_gmres_solver, pcm_block_jacobi, &
      & pcm_factor_cholesky, pcm_factor_symmetric, pcm_factor_lu, &
      & pcm_factor_complex_symmetric, pcm_factor_complex_lu, pcm_storage_dense, &
      & pcm_storage_matrix_free, pcm_nonlinear_problem, solve_pcm_lu
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_solvers

   !> Dimension of the small test operators
   integer, parameter :: n = 6
   !> Dimension of the random test operators of the iterative solver
   integer, parameter :: nr = 60
   !> Relative residual threshold after a direct solve of a well-conditioned system
   real(wp), parameter :: thr = 1.0e-13_wp
   !> Relative residual threshold for single precision
   real(sp), parameter :: thr_sp = 1.0e-5_sp

   !> Operator backed by a dense matrix
   type, extends(pcm_operator) :: dense_test_operator
      !> Operator matrix, (npt, npt)
      real(wp), allocatable :: mat(:, :)
   contains
      !> Apply the matrix or its transpose
      procedure :: apply => test_operator_apply
      !> Copy the matrix
      procedure :: assemble => test_operator_assemble
      !> Principal submatrix of the points idx
      procedure :: assemble_block => test_operator_assemble_block
   end type dense_test_operator

   !> Complex operator backed by a dense matrix
   type, extends(pcm_complex_operator) :: complex_test_operator
      !> Operator matrix, (npt, npt)
      complex(wp), allocatable :: mat(:, :)
   contains
      !> Apply the matrix or its plain transpose
      procedure :: apply => complex_test_operator_apply
      !> Copy the matrix
      procedure :: assemble => complex_test_operator_assemble
      !> Principal submatrix of the points idx
      procedure :: assemble_block => complex_test_operator_assemble_block
   end type complex_test_operator

   !> Solver extension that only implements solve_block (x = 2 rhs, 3 rhs transposed)
   type, extends(pcm_solver) :: scaling_test_solver
   contains
      !> Scale the right-hand sides
      procedure :: solve_block => scaling_test_solver_solve_block
   end type scaling_test_solver

   !> Decoupled cubics F_i(x) = x_i^3 - 2 x_i + offset
   !>
   !> - Root -1.7692923542386314; local minimum of |F| at sqrt(2/3), where
   !>   F = 0.9113 and the Jacobian is singular
   !> - Newton from 0 without a line search cycles between 0 and 1
   type, extends(pcm_nonlinear_problem) :: cubic_test_problem
      !> Constant term, zero or non-finite for the forcing tests
      real(wp) :: offset = 2.0_wp
   contains
      !> Residual
      procedure :: residual => cubic_test_problem_residual
      !> Newton step
      procedure :: step => cubic_test_problem_step
   end type cubic_test_problem

   !> Root of x^3 - 2 x + 2
   real(wp), parameter :: cubic_root = -1.7692923542386314_wp

contains

   !> Collect all exported unit tests
   !>
   !> @param[out] testsuite  registered tests
   subroutine collect_model_component_pcm_solvers(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("lapack_sytrf_sytrs", test_lapack_sytrf_sytrs), &
         & new_unittest("lapack_sytrf_sytrs_sp", test_lapack_sytrf_sytrs_sp), &
         & new_unittest("operator_apply", test_operator_apply_matches_matrix), &
         & new_unittest("dense_cholesky", test_dense_cholesky), &
         & new_unittest("dense_symmetric", test_dense_symmetric), &
         & new_unittest("dense_lu", test_dense_lu), &
         & new_unittest("dense_block_equals_single", test_dense_block_equals_single), &
         & new_unittest("dense_factor_reuse", test_dense_factor_reuse), &
         & new_unittest("dense_several_operators", test_dense_several_operators), &
         & new_unittest("dense_cholesky_failure", test_dense_cholesky_failure), &
         & new_unittest("dense_invalid_calls", test_dense_invalid_calls), &
         & new_unittest("dense_nonfinite_not_counted", test_dense_nonfinite_not_counted), &
         & new_unittest("dense_rejects_nan_rhs", test_dense_rejects_nan_rhs, should_fail=.true.), &
         & new_unittest("dense_rejects_inf_rhs", test_dense_rejects_inf_rhs, should_fail=.true.), &
         & new_unittest("dense_rejects_nan_operator", test_dense_rejects_nan_operator, should_fail=.true.), &
         & new_unittest("dense_complex_rejects_nan_rhs", test_dense_complex_rejects_nan_rhs, should_fail=.true.), &
         & new_unittest("lu_rejects_nan_rhs", test_lu_rejects_nan_rhs, should_fail=.true.), &
         & new_unittest("lu_rejects_nan_matrix", test_lu_rejects_nan_matrix, should_fail=.true.), &
         & new_unittest("lapack_complex", test_lapack_complex), &
         & new_unittest("complex_operator_apply", test_complex_operator_apply_matches_matrix), &
         & new_unittest("dense_complex_rhs", test_dense_complex_rhs), &
         & new_unittest("dense_complex_symmetric", test_dense_complex_symmetric), &
         & new_unittest("dense_complex_lu", test_dense_complex_lu), &
         & new_unittest("dense_complex_invalid_calls", test_dense_complex_invalid_calls), &
         & new_unittest("solver_complex_defaults", test_solver_complex_defaults), &
         & new_unittest("dense_setup_replaces", test_dense_setup_replaces), &
         & new_unittest("dense_setup_complex_replaces", test_dense_setup_complex_replaces), &
         & new_unittest("inexact", test_inexact), &
         & new_unittest("matrix_free", test_matrix_free), &
         & new_unittest("gmres_defaults", test_gmres_defaults), &
         & new_unittest("gmres_real_vs_dense", test_gmres_real_vs_dense), &
         & new_unittest("gmres_complex_vs_dense", test_gmres_complex_vs_dense), &
         & new_unittest("gmres_complex_rhs", test_gmres_complex_rhs), &
         & new_unittest("gmres_warm_start", test_gmres_warm_start), &
         & new_unittest("gmres_zero_rhs", test_gmres_zero_rhs), &
         & new_unittest("gmres_no_convergence", test_gmres_no_convergence), &
         & new_unittest("gmres_short_restart", test_gmres_short_restart), &
         & new_unittest("gmres_setup_rebuilds", test_gmres_setup_rebuilds), &
         & new_unittest("gmres_invalid_calls", test_gmres_invalid_calls), &
         & new_unittest("block_jacobi_size_one", test_block_jacobi_size_one), &
         & new_unittest("block_jacobi_clusters", test_block_jacobi_clusters), &
         & new_unittest("block_jacobi_complex", test_block_jacobi_complex), &
         & new_unittest("nonlinear_converges", test_nonlinear_converges), &
         & new_unittest("nonlinear_stalls", test_nonlinear_stalls), &
         & new_unittest("nonlinear_max_iter", test_nonlinear_max_iter), &
         & new_unittest("nonlinear_rejects_nan_forcing", test_nonlinear_rejects_nan_forcing, should_fail=.true.), &
         & new_unittest("nonlinear_rejects_inf_forcing", test_nonlinear_rejects_inf_forcing, should_fail=.true.), &
         & new_unittest("nonlinear_zero_forcing", test_nonlinear_zero_forcing) &
         & ]
   end subroutine collect_model_component_pcm_solvers

   !> Apply the matrix of a test operator or its transpose
   !>
   !> @param[in]  self   operator
   !> @param[in]  trans  apply the transpose if true
   !> @param[in]  x      vector of length npt
   !> @param[out] y      vector of length npt
   subroutine test_operator_apply(self, trans, x, y)
      !> Operator
      class(dense_test_operator), intent(in) :: self
      !> Apply the transpose
      logical, intent(in) :: trans
      !> Input vector
      real(wp), contiguous, intent(in) :: x(:)
      !> Result vector
      real(wp), contiguous, intent(out) :: y(:)

      if (trans) then
         y = matmul(x, self%mat)
      else
         y = matmul(self%mat, x)
      end if
   end subroutine test_operator_apply

   !> Copy the matrix of a test operator
   !>
   !> @param[in]  self   operator
   !> @param[out] a      matrix of shape (npt, npt)
   !> @param[out] error  never allocated
   subroutine test_operator_assemble(self, a, error)
      !> Operator
      class(dense_test_operator), intent(in) :: self
      !> Dense operator matrix
      real(wp), contiguous, intent(out) :: a(:, :)
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      a = self%mat
      ! intent(out) already deallocated it; this marks the dummy as defined
      if (allocated(error)) deallocate (error)
   end subroutine test_operator_assemble

   !> Principal submatrix of a test operator on the points idx
   !>
   !> - Block operators (npt = m npoint): unknowns (a - 1) npoint + idx, a = 1..m
   !>
   !> @param[in]  self   operator
   !> @param[in]  idx    point indices
   !> @param[out] a      block, (m nb, m nb)
   !> @param[out] error  never allocated
   subroutine test_operator_assemble_block(self, idx, a, error)
      !> Operator
      class(dense_test_operator), intent(in) :: self
      !> Point indices
      integer, contiguous, intent(in) :: idx(:)
      !> Block
      real(wp), contiguous, intent(out) :: a(:, :)
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      integer, allocatable :: u(:)

      call block_unknowns(self%npt, self%npoint, idx, u)
      a = self%mat(u, u)
      ! intent(out) already deallocated it; this marks the dummy as defined
      if (allocated(error)) deallocate (error)
   end subroutine test_operator_assemble_block

   !> Unknowns of the points idx of a (block) operator
   !>
   !> @param[in]  npt     operator dimension
   !> @param[in]  npoint  number of points, 0 for npt
   !> @param[in]  idx     point indices
   !> @param[out] u       unknowns (a - 1) npoint + idx, a = 1..npt / npoint
   pure subroutine block_unknowns(npt, npoint, idx, u)
      !> Operator dimension
      integer, intent(in) :: npt
      !> Number of points
      integer, intent(in) :: npoint
      !> Point indices
      integer, intent(in) :: idx(:)
      !> Unknowns
      integer, allocatable, intent(out) :: u(:)

      integer :: np, a, nb

      np = npoint
      if (np == 0) np = npt
      nb = size(idx)
      allocate (u(npt/np*nb))
      do a = 1, npt/np
         u((a - 1)*nb + 1:a*nb) = (a - 1)*np + idx
      end do
   end subroutine block_unknowns

   !> Build a test operator
   !>
   !> @param[out] op                 operator
   !> @param[in]  mat                operator matrix, (npt, npt)
   !> @param[in]  id                 operator id
   !> @param[in]  symmetric          declared symmetric
   !> @param[in]  positive_definite  declared positive definite
   subroutine new_test_operator(op, mat, id, symmetric, positive_definite)
      !> Operator
      type(dense_test_operator), intent(out) :: op
      !> Operator matrix
      real(wp), intent(in) :: mat(:, :)
      !> Operator id
      integer, intent(in) :: id
      !> Declared symmetric
      logical, intent(in) :: symmetric
      !> Declared positive definite
      logical, intent(in) :: positive_definite

      op%mat = mat
      op%npt = size(mat, 1)
      op%id = id
      op%symmetric = symmetric
      op%positive_definite = positive_definite
   end subroutine new_test_operator

   !> Symmetric positive definite matrix: Hilbert matrix plus identity
   !>
   !> @param[out] a  matrix, (n, n)
   pure subroutine spd_matrix(a)
      !> Matrix
      real(wp), intent(out) :: a(n, n)

      integer :: i, j

      do j = 1, n
         do i = 1, n
            a(i, j) = 1.0_wp/real(i + j - 1, wp)
         end do
         a(j, j) = a(j, j) + 1.0_wp
      end do
   end subroutine spd_matrix

   !> Symmetric indefinite matrix Q diag(d) Q**T with a Householder reflector Q
   !>
   !> - Eigenvalues 3, -2, 1, -4, 2.5, -1.5 (condition number 4)
   !>
   !> @param[out] a  matrix, (n, n)
   pure subroutine indefinite_matrix(a)
      !> Matrix
      real(wp), intent(out) :: a(n, n)

      real(wp), parameter :: d(n) = [3.0_wp, -2.0_wp, 1.0_wp, -4.0_wp, 2.5_wp, -1.5_wp]
      real(wp), parameter :: v(n) = [1.0_wp, -2.0_wp, 0.5_wp, 3.0_wp, -1.0_wp, 2.0_wp]
      real(wp) :: q(n, n)
      integer :: i, j

      do j = 1, n
         do i = 1, n
            q(i, j) = -2.0_wp*v(i)*v(j)/sum(v**2)
         end do
         q(j, j) = q(j, j) + 1.0_wp
      end do
      do j = 1, n
         do i = 1, n
            a(i, j) = sum(q(i, :)*d*q(j, :))
         end do
      end do
   end subroutine indefinite_matrix

   !> General non-symmetric matrix with a dominant diagonal
   !>
   !> @param[out] a  matrix, (n, n)
   pure subroutine general_matrix(a)
      !> Matrix
      real(wp), intent(out) :: a(n, n)

      integer :: i, j

      do j = 1, n
         do i = 1, n
            a(i, j) = 1.0_wp/real(i + 2*j - 2, wp) + 0.1_wp*real(i - j, wp)
         end do
         a(j, j) = a(j, j) + 2.0_wp
      end do
   end subroutine general_matrix

   !> Right-hand sides with distinct columns
   !>
   !> @param[out] b  right-hand sides, (n, nrhs)
   pure subroutine rhs_block(b)
      !> Right-hand sides
      real(wp), intent(out) :: b(:, :)

      integer :: i, k

      do k = 1, size(b, 2)
         do i = 1, size(b, 1)
            b(i, k) = sin(real(i*k, wp)) + 0.25_wp*real(k, wp)
         end do
      end do
   end subroutine rhs_block

   !> Relative residual norm(a x - b) / norm(b)
   !>
   !> @param[in] a  matrix, (n, n)
   !> @param[in] x  solution, (n)
   !> @param[in] b  right-hand side, (n)
   pure function rel_residual(a, x, b) result(r)
      !> Matrix
      real(wp), intent(in) :: a(:, :)
      !> Solution
      real(wp), intent(in) :: x(:)
      !> Right-hand side
      real(wp), intent(in) :: b(:)
      !> Relative residual
      real(wp) :: r

      r = norm2(matmul(a, x) - b)/norm2(b)
   end function rel_residual

   !> Solve op x = b and op**T x = b and check both residuals
   !>
   !> @param[out]    error   test failure
   !> @param[in,out] solver  solver under test
   !> @param[in]     op      operator
   !> @param[in]     tol     relative residual threshold
   subroutine check_both_orientations(error, solver, op, tol)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Solver under test
      class(pcm_solver), intent(inout) :: solver
      !> Operator
      type(dense_test_operator), intent(in) :: op
      !> Relative residual threshold
      real(wp), intent(in) :: tol

      type(moist_error_type), allocatable :: err
      real(wp) :: b(n, 1), x(n), xt(n)

      call rhs_block(b)
      x = 0.0_wp
      call solver%solve(op, .false., b(:, 1), x, err)
      call check_moist_error(error, err, "normal solve")
      if (allocated(error)) return
      call check(error, rel_residual(op%mat, x, b(:, 1)) < tol, "normal solve residual")
      if (allocated(error)) return

      xt = 0.0_wp
      call solver%solve(op, .true., b(:, 1), xt, err)
      call check_moist_error(error, err, "transpose solve")
      if (allocated(error)) return
      call check(error, rel_residual(transpose(op%mat), xt, b(:, 1)) < tol, &
         & "transpose solve residual")
      if (allocated(error)) return
   end subroutine check_both_orientations

   !> sytrf/sytrs wrappers on a zero-diagonal symmetric system, both triangles
   subroutine test_lapack_sytrf_sytrs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      ! Zero diagonal forces a 2x2 pivot; determinant 12
      real(wp), parameter :: a(3, 3) = reshape([0.0_wp, 1.0_wp, 2.0_wp, 1.0_wp, 0.0_wp, &
                                                3.0_wp, 2.0_wp, 3.0_wp, 0.0_wp], [3, 3])
      real(wp), parameter :: x_ref(3, 2) = reshape([1.0_wp, -2.0_wp, 0.5_wp, &
                                                    -0.25_wp, 3.0_wp, 1.5_wp], [3, 2])
      character(len=1), parameter :: triangles(2) = ["l", "u"]
      real(wp) :: f(3, 3), x(3, 2)
      integer :: ipiv(3), info, it

      do it = 1, size(triangles)
         f = a
         call sytrf(f, ipiv, info, uplo=triangles(it))
         call check(error, info, 0, "sytrf info, uplo = "//triangles(it))
         if (allocated(error)) return
         x = matmul(a, x_ref)
         call sytrs(f, x, ipiv, info, uplo=triangles(it))
         call check(error, info, 0, "sytrs info, uplo = "//triangles(it))
         if (allocated(error)) return
         call check(error, maxval(abs(x - x_ref)) < thr, "sytrs solution, uplo = "//triangles(it))
         if (allocated(error)) return
      end do

      ! Default triangle is upper
      f = a
      call sytrf(f, ipiv, info)
      call check(error, info, 0, "sytrf info, default uplo")
      if (allocated(error)) return
      x = matmul(a, x_ref)
      call sytrs(f, x, ipiv, info)
      call check(error, info, 0, "sytrs info, default uplo")
      if (allocated(error)) return
      call check(error, maxval(abs(x - x_ref)) < thr, "sytrs solution, default uplo")
      if (allocated(error)) return

      ! Exactly singular matrix: positive info
      f = 0.0_wp
      call sytrf(f, ipiv, info, uplo="l")
      call check(error, info > 0, "sytrf reports a singular matrix")
   end subroutine test_lapack_sytrf_sytrs

   !> Single precision sytrf/sytrs wrappers
   subroutine test_lapack_sytrf_sytrs_sp(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(sp), parameter :: a(3, 3) = reshape([0.0_sp, 1.0_sp, 2.0_sp, 1.0_sp, 0.0_sp, &
                                                3.0_sp, 2.0_sp, 3.0_sp, 0.0_sp], [3, 3])
      real(sp), parameter :: x_ref(3, 1) = reshape([1.0_sp, -2.0_sp, 0.5_sp], [3, 1])
      real(sp) :: f(3, 3), x(3, 1)
      integer :: ipiv(3), info

      f = a
      call sytrf(f, ipiv, info, uplo="l")
      call check(error, info, 0, "ssytrf info")
      if (allocated(error)) return
      x = matmul(a, x_ref)
      call sytrs(f, x, ipiv, info, uplo="l")
      call check(error, info, 0, "ssytrs info")
      if (allocated(error)) return
      call check(error, maxval(abs(x - x_ref)) < thr_sp, "ssytrs solution")
   end subroutine test_lapack_sytrf_sytrs_sp

   !> The test operator applies its matrix and transpose consistently with assemble
   subroutine test_operator_apply_matches_matrix(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(moist_error_type), allocatable :: err
      real(wp) :: g(n, n), a(n, n), b(n, 1), y(n)

      call general_matrix(g)
      call new_test_operator(op, g, 1, .false., .false.)
      call op%assemble(a, err)
      call check_moist_error(error, err, "assemble")
      if (allocated(error)) return
      call rhs_block(b)
      call op%apply(.false., b(:, 1), y)
      call check(error, maxval(abs(y - matmul(a, b(:, 1)))) < thr, "apply")
      if (allocated(error)) return
      call op%apply(.true., b(:, 1), y)
      call check(error, maxval(abs(y - matmul(transpose(a), b(:, 1)))) < thr, "apply transpose")
   end subroutine test_operator_apply_matches_matrix

   !> Positive definite operator: Cholesky, both orientations
   subroutine test_dense_cholesky(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n)

      call spd_matrix(a)
      call new_test_operator(op, a, 1, .true., .true.)
      call solver%reset(1)
      call solver%setup(op, err)
      call check_moist_error(error, err, "setup")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_cholesky, "factorization kind")
      if (allocated(error)) return
      call check_both_orientations(error, solver, op, thr)
   end subroutine test_dense_cholesky

   !> Symmetric indefinite operator: Bunch-Kaufman, both orientations
   subroutine test_dense_symmetric(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n)

      call indefinite_matrix(a)
      call new_test_operator(op, a, 1, .true., .false.)
      call solver%reset(1)
      call solver%setup(op, err)
      call check_moist_error(error, err, "setup")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_symmetric, "factorization kind")
      if (allocated(error)) return
      call check_both_orientations(error, solver, op, thr)
   end subroutine test_dense_symmetric

   !> General operator: LU; the transpose solve differs from the normal one
   subroutine test_dense_lu(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 1), x(n), xt(n)

      call general_matrix(a)
      call new_test_operator(op, a, 1, .false., .false.)
      call solver%reset(1)
      call solver%setup(op, err)
      call check_moist_error(error, err, "setup")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_lu, "factorization kind")
      if (allocated(error)) return
      call check_both_orientations(error, solver, op, thr)
      if (allocated(error)) return

      call rhs_block(b)
      call solver%solve(op, .false., b(:, 1), x, err)
      call check_moist_error(error, err, "normal solve")
      if (allocated(error)) return
      call solver%solve(op, .true., b(:, 1), xt, err)
      call check_moist_error(error, err, "transpose solve")
      if (allocated(error)) return
      call check(error, maxval(abs(x - xt)) > 1.0e-3_wp, &
         & "transpose solve differs from normal solve")
   end subroutine test_dense_lu

   !> solve_block with several right-hand sides equals repeated solve, all kinds
   subroutine test_dense_block_equals_single(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nrhs = 3
      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, nrhs), x(n, nrhs), xk(n)
      integer :: kind, k, it
      logical :: trans

      do kind = 1, 3
         select case (kind)
         case (pcm_factor_cholesky)
            call spd_matrix(a)
            call new_test_operator(op, a, 1, .true., .true.)
         case (pcm_factor_symmetric)
            call indefinite_matrix(a)
            call new_test_operator(op, a, 1, .true., .false.)
         case default
            call general_matrix(a)
            call new_test_operator(op, a, 1, .false., .false.)
         end select
         call solver%reset(1)
         call rhs_block(b)
         do it = 1, 2
            trans = it == 2
            call solver%solve_block(op, trans, b, x, err)
            call check_moist_error(error, err, "solve_block")
            if (allocated(error)) return
            call check(error, solver%factorization(1), kind, "factorization kind")
            if (allocated(error)) return
            do k = 1, nrhs
               call solver%solve(op, trans, b(:, k), xk, err)
               call check_moist_error(error, err, "solve")
               if (allocated(error)) return
               call check(error, norm2(x(:, k) - xk) <= thr*norm2(xk), &
                  & "block column equals single solve")
               if (allocated(error)) return
               if (trans) then
                  call check(error, rel_residual(transpose(a), x(:, k), b(:, k)) < thr, &
                     & "block transpose residual")
               else
                  call check(error, rel_residual(a, x(:, k), b(:, k)) < thr, "block residual")
               end if
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_dense_block_equals_single

   !> Each setup factors anew; solves reuse the factors and factor on first use
   subroutine test_dense_factor_reuse(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 2), x(n, 2), x1(n)

      call general_matrix(a)
      call new_test_operator(op, a, 1, .false., .false.)
      call rhs_block(b)
      call solver%reset(1)
      call check(error, solver%factorization(1), 0, "no factors after reset")
      if (allocated(error)) return

      call solver%setup(op, err)
      call check_moist_error(error, err, "setup")
      if (allocated(error)) return
      call solver%setup(op, err)
      call check_moist_error(error, err, "second setup")
      if (allocated(error)) return
      call solver%solve(op, .false., b(:, 1), x1, err)
      call check_moist_error(error, err, "solve")
      if (allocated(error)) return
      call solver%solve(op, .true., b(:, 1), x1, err)
      call check_moist_error(error, err, "transpose solve")
      if (allocated(error)) return
      call solver%solve_block(op, .false., b, x, err)
      call check_moist_error(error, err, "solve_block")
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 2, "one factorization per setup")
      if (allocated(error)) return
      call check(error, solver%n_solves, 4, "solves count right-hand sides")
      if (allocated(error)) return

      ! After reset the factors are gone; the next solve factors again
      call solver%reset(1)
      call check(error, solver%factorization(1), 0, "factors discarded by reset")
      if (allocated(error)) return
      call solver%solve(op, .false., b(:, 1), x1, err)
      call check_moist_error(error, err, "solve after reset")
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 3, "refactored after reset")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_lu, "factors after lazy setup")
      if (allocated(error)) return
      call check(error, rel_residual(a, x1, b(:, 1)) < thr, "residual after reset")
   end subroutine test_dense_factor_reuse

   !> Factors of several operators are kept apart by id
   subroutine test_dense_several_operators(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op_spd, op_sym, op_gen
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n)

      call spd_matrix(a)
      call new_test_operator(op_spd, a, 1, .true., .true.)
      call indefinite_matrix(a)
      call new_test_operator(op_sym, a, 2, .true., .false.)
      call general_matrix(a)
      call new_test_operator(op_gen, a, 3, .false., .false.)
      call solver%reset(3)
      call solver%setup(op_spd, err)
      if (.not. allocated(err)) call solver%setup(op_sym, err)
      if (.not. allocated(err)) call solver%setup(op_gen, err)
      call check_moist_error(error, err, "setup")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_cholesky, "operator 1 kind")
      if (allocated(error)) return
      call check(error, solver%factorization(2), pcm_factor_symmetric, "operator 2 kind")
      if (allocated(error)) return
      call check(error, solver%factorization(3), pcm_factor_lu, "operator 3 kind")
      if (allocated(error)) return
      call check(error, solver%factorization(4), 0, "id outside the reset range")
      if (allocated(error)) return

      call check_both_orientations(error, solver, op_spd, thr)
      if (allocated(error)) return
      call check_both_orientations(error, solver, op_sym, thr)
      if (allocated(error)) return
      call check_both_orientations(error, solver, op_gen, thr)
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 3, "one factorization per operator")
   end subroutine test_dense_several_operators

   !> A non-SPD operator declared positive definite fails in Cholesky, no LU fallback
   subroutine test_dense_cholesky_failure(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 1), x(n)

      call indefinite_matrix(a)
      call new_test_operator(op, a, 1, .true., .true.)
      call solver%reset(1)
      call solver%setup(op, err)
      call check(error, allocated(err), "setup reports a failed Cholesky factorization")
      if (allocated(error)) return
      call check(error, index(err%message, "[pcm_dense_setup]") == 1, "error names the routine")
      if (allocated(error)) return
      call check(error, solver%factorization(1), 0, "operator stays unfactored")
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 0, "no factorization counted")
      if (allocated(error)) return

      call rhs_block(b)
      call solver%solve(op, .false., b(:, 1), x, err)
      call check(error, allocated(err), "solve propagates the setup error")
      if (allocated(error)) return
      call check(error, solver%n_solves, 0, "no solve counted")
   end subroutine test_dense_cholesky_failure

   !> Missing reset, ids outside the reset range and shape mismatches are errors
   subroutine test_dense_invalid_calls(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 2), x(n, 1), x1(n - 1)

      call general_matrix(a)
      call new_test_operator(op, a, 2, .false., .false.)
      call rhs_block(b)

      call solver%setup(op, err)
      call check(error, allocated(err), "setup before reset")
      if (allocated(error)) return
      call solver%reset(1)
      call solver%setup(op, err)
      call check(error, allocated(err), "id outside the reset range")
      if (allocated(error)) return

      call solver%reset(2)
      call solver%solve_block(op, .false., b, x, err)
      call check(error, allocated(err), "solution block of the wrong shape")
      if (allocated(error)) return
      call solver%solve(op, .false., b(:, 1), x1, err)
      call check(error, allocated(err), "solution of the wrong length")
      if (allocated(error)) return
      call solver%solve(op, .false., b(:n - 1, 1), x1, err)
      call check(error, allocated(err), "right-hand side of the wrong length")
      if (allocated(error)) return
      call check(error, solver%n_solves, 0, "no solve counted")
   end subroutine test_dense_invalid_calls

   !> Dense solve rejects a NaN right-hand side
   subroutine test_dense_rejects_nan_rhs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 1), x(n)

      call general_matrix(a)
      call rhs_block(b)
      b(2, 1) = ieee_value(1.0_wp, ieee_quiet_nan)
      call new_test_operator(op, a, 1, .false., .false.)
      call solver%reset(1)
      call solver%solve(op, .false., b(:, 1), x, err)
      call check_moist_error(error, err)
   end subroutine test_dense_rejects_nan_rhs

   !> Dense transposed solve rejects an infinite right-hand side
   subroutine test_dense_rejects_inf_rhs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 1), x(n)

      call general_matrix(a)
      call rhs_block(b)
      b(2, 1) = ieee_value(1.0_wp, ieee_positive_inf)
      call new_test_operator(op, a, 1, .false., .false.)
      call solver%reset(1)
      call solver%solve(op, .true., b(:, 1), x, err)
      call check_moist_error(error, err)
   end subroutine test_dense_rejects_inf_rhs

   !> Dense solve rejects a NaN operator with a finite forcing
   subroutine test_dense_rejects_nan_operator(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 1), x(n)

      call general_matrix(a)
      call rhs_block(b)
      a(1, 2) = ieee_value(1.0_wp, ieee_quiet_nan)
      call new_test_operator(op, a, 1, .false., .false.)
      call solver%reset(1)
      call solver%solve(op, .false., b(:, 1), x, err)
      call check_moist_error(error, err)
   end subroutine test_dense_rejects_nan_operator

   !> Dense complex solve rejects a NaN right-hand side
   subroutine test_dense_complex_rejects_nan_rhs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(complex_test_operator) :: zop
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      complex(wp) :: za(nr, nr), zb(nr), zx(nr)
      real(wp) :: re(nr, 2)

      call random_complex_matrix(11, .false., za)
      call new_complex_test_operator(zop, za, 1, .false.)
      call random_rhs(13, re)
      zb = cmplx(re(:, 1), re(:, 2), kind=wp)
      zb(3) = cmplx(1.0_wp, ieee_value(1.0_wp, ieee_quiet_nan), kind=wp)
      call solver%reset(1)
      call solver%solve_complex(zop, .false., zb, zx, err)
      call check_moist_error(error, err)
   end subroutine test_dense_complex_rejects_nan_rhs

   !> Legacy LU wrapper rejects a NaN right-hand side
   subroutine test_lu_rejects_nan_rhs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 1), q(n)

      call general_matrix(a)
      call rhs_block(b)
      b(4, 1) = ieee_value(1.0_wp, ieee_quiet_nan)
      call solve_pcm_lu(a, b(:, 1), q, err)
      call check_moist_error(error, err)
   end subroutine test_lu_rejects_nan_rhs

   !> Legacy LU wrapper rejects a NaN matrix
   subroutine test_lu_rejects_nan_matrix(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 1), q(n)

      call general_matrix(a)
      call rhs_block(b)
      a(2, 2) = ieee_value(1.0_wp, ieee_quiet_nan)
      call solve_pcm_lu(a, b(:, 1), q, err)
      call check_moist_error(error, err)
   end subroutine test_lu_rejects_nan_matrix

   !> A rejected non-finite right-hand side is not counted as a solve
   subroutine test_dense_nonfinite_not_counted(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 1), x(n)

      call general_matrix(a)
      call rhs_block(b)
      b(2, 1) = ieee_value(1.0_wp, ieee_positive_inf)
      call new_test_operator(op, a, 1, .false., .false.)
      call solver%reset(1)
      call solver%solve(op, .true., b(:, 1), x, err)
      call check(error, solver%n_solves, 0, "no failed solve counted")
   end subroutine test_dense_nonfinite_not_counted

   !> Apply the matrix of a complex test operator or its plain transpose
   !>
   !> @param[in]  self   operator
   !> @param[in]  trans  apply the transpose (not conjugated) if true
   !> @param[in]  x      vector of length npt
   !> @param[out] y      vector of length npt
   subroutine complex_test_operator_apply(self, trans, x, y)
      !> Operator
      class(complex_test_operator), intent(in) :: self
      !> Apply the transpose
      logical, intent(in) :: trans
      !> Input vector
      complex(wp), contiguous, intent(in) :: x(:)
      !> Result vector
      complex(wp), contiguous, intent(out) :: y(:)

      if (trans) then
         y = matmul(x, self%mat)
      else
         y = matmul(self%mat, x)
      end if
   end subroutine complex_test_operator_apply

   !> Copy the matrix of a complex test operator
   !>
   !> @param[in]  self   operator
   !> @param[out] a      matrix of shape (npt, npt)
   !> @param[out] error  never allocated
   subroutine complex_test_operator_assemble(self, a, error)
      !> Operator
      class(complex_test_operator), intent(in) :: self
      !> Dense operator matrix
      complex(wp), contiguous, intent(out) :: a(:, :)
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      a = self%mat
      ! intent(out) already deallocated it; this marks the dummy as defined
      if (allocated(error)) deallocate (error)
   end subroutine complex_test_operator_assemble

   !> Principal submatrix of a complex test operator on the points idx
   !>
   !> @param[in]  self   operator
   !> @param[in]  idx    point indices
   !> @param[out] a      block, (m nb, m nb)
   !> @param[out] error  never allocated
   subroutine complex_test_operator_assemble_block(self, idx, a, error)
      !> Operator
      class(complex_test_operator), intent(in) :: self
      !> Point indices
      integer, contiguous, intent(in) :: idx(:)
      !> Block
      complex(wp), contiguous, intent(out) :: a(:, :)
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      integer, allocatable :: u(:)

      call block_unknowns(self%npt, self%npoint, idx, u)
      a = self%mat(u, u)
      ! intent(out) already deallocated it; this marks the dummy as defined
      if (allocated(error)) deallocate (error)
   end subroutine complex_test_operator_assemble_block

   !> Build a complex test operator
   !>
   !> @param[out] op         operator
   !> @param[in]  mat        operator matrix, (npt, npt)
   !> @param[in]  id         operator id
   !> @param[in]  symmetric  declared symmetric (mat = mat**T)
   subroutine new_complex_test_operator(op, mat, id, symmetric)
      !> Operator
      type(complex_test_operator), intent(out) :: op
      !> Operator matrix
      complex(wp), intent(in) :: mat(:, :)
      !> Operator id
      integer, intent(in) :: id
      !> Declared symmetric
      logical, intent(in) :: symmetric

      op%mat = mat
      op%npt = size(mat, 1)
      op%id = id
      op%symmetric = symmetric
   end subroutine new_complex_test_operator

   !> Scale the right-hand sides: x = 2 rhs, or x = 3 rhs for the transpose
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     operator, unused
   !> @param[in]     trans  solve with the transpose if true
   !> @param[in]     rhs    right-hand sides, (npt, nrhs)
   !> @param[in,out] x      ignored on entry, solutions on exit, (npt, nrhs)
   !> @param[out]    error  never allocated
   subroutine scaling_test_solver_solve_block(self, op, trans, rhs, x, error)
      !> Solver state
      class(scaling_test_solver), intent(inout) :: self
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Solve with the transpose
      logical, intent(in) :: trans
      !> Right-hand sides
      real(wp), contiguous, intent(in) :: rhs(:, :)
      !> Ignored on entry, solutions on exit
      real(wp), contiguous, intent(inout) :: x(:, :)
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      associate (unused_self => self, unused_op => op)
      end associate
      x = merge(3.0_wp, 2.0_wp, trans)*rhs
      ! intent(out) already deallocated it; this marks the dummy as defined
      if (allocated(error)) deallocate (error)
   end subroutine scaling_test_solver_solve_block

   !> Residual of the decoupled cubics
   !>
   !> @param[in]  self   problem
   !> @param[in]  x      unknowns, length npt
   !> @param[out] f      residual, length npt
   !> @param[out] error  never allocated
   subroutine cubic_test_problem_residual(self, x, f, error)
      !> Problem
      class(cubic_test_problem), intent(in) :: self
      !> Unknowns
      real(wp), contiguous, intent(in) :: x(:)
      !> Residual
      real(wp), contiguous, intent(out) :: f(:)
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      f = x**3 - 2.0_wp*x + self%offset
      ! intent(out) already deallocated it; this marks the dummy as defined
      if (allocated(error)) deallocate (error)
   end subroutine cubic_test_problem_residual

   !> Newton step of the decoupled cubics, p = -f/(3 x^2 - 2)
   !>
   !> @param[in]  self   problem
   !> @param[in]  x      unknowns, length npt
   !> @param[in]  f      residual at x, length npt
   !> @param[out] p      step, length npt
   !> @param[out] error  never allocated
   subroutine cubic_test_problem_step(self, x, f, p, error)
      !> Problem
      class(cubic_test_problem), intent(in) :: self
      !> Unknowns
      real(wp), contiguous, intent(in) :: x(:)
      !> Residual at x
      real(wp), contiguous, intent(in) :: f(:)
      !> Newton step
      real(wp), contiguous, intent(out) :: p(:)
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      associate (unused_self => self)
      end associate
      p = -f/(3.0_wp*x**2 - 2.0_wp)
      ! intent(out) already deallocated it; this marks the dummy as defined
      if (allocated(error)) deallocate (error)
   end subroutine cubic_test_problem_step

   !> Complex symmetric, non-Hermitian matrix: indefinite + i spd
   !>
   !> - Exactly symmetric; the imaginary part is positive definite with
   !>   eigenvalues >= 1, so the smallest singular value is at least 1
   !>
   !> @param[out] a  matrix, (n, n)
   pure subroutine complex_symmetric_matrix(a)
      !> Matrix
      complex(wp), intent(out) :: a(n, n)

      real(wp) :: re(n, n), im(n, n)

      call indefinite_matrix(re)
      call spd_matrix(im)
      ! indefinite_matrix is symmetric only up to rounding
      re = 0.5_wp*(re + transpose(re))
      a = cmplx(re, im, kind=wp)
   end subroutine complex_symmetric_matrix

   !> Complex general matrix: general + i spd / 2
   !>
   !> @param[out] a  matrix, (n, n)
   pure subroutine complex_general_matrix(a)
      !> Matrix
      complex(wp), intent(out) :: a(n, n)

      real(wp) :: re(n, n), im(n, n)

      call general_matrix(re)
      call spd_matrix(im)
      a = cmplx(re, 0.5_wp*im, kind=wp)
   end subroutine complex_general_matrix

   !> Complex right-hand side with nonzero real and imaginary parts
   !>
   !> @param[out] b  right-hand side, (n)
   pure subroutine complex_rhs(b)
      !> Right-hand side
      complex(wp), intent(out) :: b(:)

      integer :: i

      do i = 1, size(b)
         b(i) = cmplx(sin(real(i, wp)), cos(real(2*i, wp)) + 0.5_wp, kind=wp)
      end do
   end subroutine complex_rhs

   !> Euclidean norm of a complex vector
   !>
   !> @param[in] v  vector
   pure function znorm(v) result(r)
      !> Vector
      complex(wp), intent(in) :: v(:)
      !> Norm
      real(wp) :: r

      r = norm2(abs(v))
   end function znorm

   !> Real operator, complex right-hand side: both orientations and parts
   !>
   !> - The residual of the complex solve is below tol
   !> - Real and imaginary parts equal the real solves of the parts
   !>
   !> @param[out]    error   test failure
   !> @param[in,out] solver  solver under test
   !> @param[in]     op      real operator
   !> @param[in]     tol     relative residual threshold
   subroutine check_complex_rhs(error, solver, op, tol)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Solver under test
      class(pcm_solver), intent(inout) :: solver
      !> Operator
      type(dense_test_operator), intent(in) :: op
      !> Relative residual threshold
      real(wp), intent(in) :: tol

      type(moist_error_type), allocatable :: err
      complex(wp) :: b(n), x(n)
      real(wp) :: a(n, n), x_re(n), x_im(n)
      integer :: it
      logical :: trans

      call complex_rhs(b)
      do it = 1, 2
         trans = it == 2
         if (trans) then
            a = transpose(op%mat)
         else
            a = op%mat
         end if
         x = (0.0_wp, 0.0_wp)
         call solver%solve_complex(op, trans, b, x, err)
         call check_moist_error(error, err, "complex solve, trans = "//merge("t", "f", trans))
         if (allocated(error)) return
         call check(error, znorm(matmul(a, x) - b) < tol*znorm(b), &
            & "complex residual, trans = "//merge("t", "f", trans))
         if (allocated(error)) return

         x_re = 0.0_wp
         call solver%solve(op, trans, real(b, wp), x_re, err)
         call check_moist_error(error, err, "real part solve")
         if (allocated(error)) return
         x_im = 0.0_wp
         call solver%solve(op, trans, aimag(b), x_im, err)
         call check_moist_error(error, err, "imaginary part solve")
         if (allocated(error)) return
         call check(error, norm2(real(x, wp) - x_re) <= tol*norm2(x_re), &
            & "real part equals real solve, trans = "//merge("t", "f", trans))
         if (allocated(error)) return
         call check(error, norm2(aimag(x) - x_im) <= tol*norm2(x_im), &
            & "imaginary part equals real solve, trans = "//merge("t", "f", trans))
         if (allocated(error)) return
      end do
   end subroutine check_complex_rhs

   !> Complex operator: solve op x = b and op**T x = b, check with op%apply
   !>
   !> @param[out]    error   test failure
   !> @param[in,out] solver  solver under test
   !> @param[in]     op      complex operator
   !> @param[in]     tol     relative residual threshold
   !> @param[out]    x       normal solution, (n)
   !> @param[out]    xt      transposed solution, (n)
   subroutine check_complex_orientations(error, solver, op, tol, x, xt)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Solver under test
      class(pcm_solver), intent(inout) :: solver
      !> Operator
      type(complex_test_operator), intent(in) :: op
      !> Relative residual threshold
      real(wp), intent(in) :: tol
      !> Normal solution
      complex(wp), intent(out) :: x(n)
      !> Transposed solution
      complex(wp), intent(out) :: xt(n)

      type(moist_error_type), allocatable :: err
      complex(wp) :: b(n), y(n)

      call complex_rhs(b)
      x = (0.0_wp, 0.0_wp)
      call solver%solve_complex(op, .false., b, x, err)
      call check_moist_error(error, err, "normal complex solve")
      if (allocated(error)) return
      call op%apply(.false., x, y)
      call check(error, znorm(y - b) < tol*znorm(b), "normal complex residual")
      if (allocated(error)) return

      xt = (0.0_wp, 0.0_wp)
      call solver%solve_complex(op, .true., b, xt, err)
      call check_moist_error(error, err, "transpose complex solve")
      if (allocated(error)) return
      call op%apply(.true., xt, y)
      call check(error, znorm(y - b) < tol*znorm(b), "transpose complex residual")
   end subroutine check_complex_orientations

   !> Complex getrf/getrs and sytrf/sytrs wrappers
   !>
   !> - zsytrf/zsytrs on a complex symmetric (non-Hermitian) system with a
   !>   zero diagonal, both triangles and the default
   !> - zgetrs with trans = "t" solves the plain transpose
   subroutine test_lapack_complex(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      ! Zero diagonal forces a 2x2 pivot; symmetric, not Hermitian
      complex(wp), parameter :: s(3, 3) = reshape([ &
         & (0.0_wp, 0.0_wp), (1.0_wp, 2.0_wp), (2.0_wp, -1.0_wp), &
         & (1.0_wp, 2.0_wp), (0.0_wp, 0.0_wp), (3.0_wp, 0.5_wp), &
         & (2.0_wp, -1.0_wp), (3.0_wp, 0.5_wp), (0.0_wp, 0.0_wp)], [3, 3])
      complex(wp), parameter :: g(3, 3) = reshape([ &
         & (4.0_wp, 1.0_wp), (1.0_wp, -2.0_wp), (0.5_wp, 0.0_wp), &
         & (-1.0_wp, 0.5_wp), (3.0_wp, -1.0_wp), (2.0_wp, 1.0_wp), &
         & (0.0_wp, 2.0_wp), (1.0_wp, 1.0_wp), (5.0_wp, 0.5_wp)], [3, 3])
      complex(wp), parameter :: x_ref(3, 2) = reshape([ &
         & (1.0_wp, 0.5_wp), (-2.0_wp, 1.0_wp), (0.5_wp, -3.0_wp), &
         & (-0.25_wp, 0.0_wp), (3.0_wp, 2.0_wp), (1.5_wp, -0.5_wp)], [3, 2])
      character(len=1), parameter :: triangles(2) = ["l", "u"]
      complex(wp) :: f(3, 3), x(3, 2)
      integer :: ipiv(3), info, it

      do it = 1, size(triangles)
         f = s
         call sytrf(f, ipiv, info, uplo=triangles(it))
         call check(error, info, 0, "zsytrf info, uplo = "//triangles(it))
         if (allocated(error)) return
         x = matmul(s, x_ref)
         call sytrs(f, x, ipiv, info, uplo=triangles(it))
         call check(error, info, 0, "zsytrs info, uplo = "//triangles(it))
         if (allocated(error)) return
         call check(error, maxval(abs(x - x_ref)) < thr, "zsytrs solution, uplo = "//triangles(it))
         if (allocated(error)) return
      end do

      ! Default triangle is upper
      f = s
      call sytrf(f, ipiv, info)
      call check(error, info, 0, "zsytrf info, default uplo")
      if (allocated(error)) return
      x = matmul(s, x_ref)
      call sytrs(f, x, ipiv, info)
      call check(error, info, 0, "zsytrs info, default uplo")
      if (allocated(error)) return
      call check(error, maxval(abs(x - x_ref)) < thr, "zsytrs solution, default uplo")
      if (allocated(error)) return

      ! Exactly singular matrix: positive info
      f = (0.0_wp, 0.0_wp)
      call sytrf(f, ipiv, info, uplo="l")
      call check(error, info > 0, "zsytrf reports a singular matrix")
      if (allocated(error)) return

      ! LU: normal, plain transpose and default orientation
      f = g
      call getrf(f, ipiv, info)
      call check(error, info, 0, "zgetrf info")
      if (allocated(error)) return
      x = matmul(g, x_ref)
      call getrs(f, x, ipiv, info, trans="n")
      call check(error, info, 0, "zgetrs info, trans = n")
      if (allocated(error)) return
      call check(error, maxval(abs(x - x_ref)) < thr, "zgetrs solution, trans = n")
      if (allocated(error)) return
      x = matmul(transpose(g), x_ref)
      call getrs(f, x, ipiv, info, trans="t")
      call check(error, info, 0, "zgetrs info, trans = t")
      if (allocated(error)) return
      call check(error, maxval(abs(x - x_ref)) < thr, "zgetrs solution, trans = t")
      if (allocated(error)) return
      x = matmul(g, x_ref)
      call getrs(f, x, ipiv, info)
      call check(error, info, 0, "zgetrs info, default trans")
      if (allocated(error)) return
      call check(error, maxval(abs(x - x_ref)) < thr, "zgetrs solution, default trans")
      if (allocated(error)) return

      f = (0.0_wp, 0.0_wp)
      call getrf(f, ipiv, info)
      call check(error, info > 0, "zgetrf reports a singular matrix")
   end subroutine test_lapack_complex

   !> The complex test operator applies its matrix and plain transpose
   subroutine test_complex_operator_apply_matches_matrix(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(complex_test_operator) :: op
      type(moist_error_type), allocatable :: err
      complex(wp) :: g(n, n), a(n, n), b(n), y(n)

      call complex_general_matrix(g)
      call new_complex_test_operator(op, g, 1, .false.)
      call op%assemble(a, err)
      call check_moist_error(error, err, "assemble")
      if (allocated(error)) return
      call complex_rhs(b)
      call op%apply(.false., b, y)
      call check(error, maxval(abs(y - matmul(a, b))) < thr, "apply")
      if (allocated(error)) return
      call op%apply(.true., b, y)
      call check(error, maxval(abs(y - matmul(transpose(a), b))) < thr, "apply transpose")
   end subroutine test_complex_operator_apply_matches_matrix

   !> Dense solver, real operator, complex right-hand side: all real kinds
   !> - The real factors are reused; each complex solve counts two solves
   subroutine test_dense_complex_rhs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n)
      complex(wp) :: b(n), x(n)
      integer :: kind, n_solves

      do kind = 1, 3
         select case (kind)
         case (pcm_factor_cholesky)
            call spd_matrix(a)
            call new_test_operator(op, a, 1, .true., .true.)
         case (pcm_factor_symmetric)
            call indefinite_matrix(a)
            call new_test_operator(op, a, 1, .true., .false.)
         case default
            call general_matrix(a)
            call new_test_operator(op, a, 1, .false., .false.)
         end select
         call solver%reset(1)
         call check_complex_rhs(error, solver, op, thr)
         if (allocated(error)) return
         call check(error, solver%factorization(1), kind, "factorization kind")
         if (allocated(error)) return
         call check(error, solver%n_factorizations, kind, "one factorization per operator")
         if (allocated(error)) return
      end do

      call complex_rhs(b)
      n_solves = solver%n_solves
      call solver%solve_complex(op, .false., b, x, err)
      call check_moist_error(error, err, "complex solve")
      if (allocated(error)) return
      call check(error, solver%n_solves, n_solves + 2, "complex solve counts two solves")
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 3, "factors reused")
      if (allocated(error)) return

      call solver%solve_complex(op, .false., b(:n - 1), x, err)
      call check(error, allocated(err), "right-hand side of the wrong length")
      if (allocated(error)) return
      call solver%solve_complex(op, .false., b, x(:n - 1), err)
      call check(error, allocated(err), "solution of the wrong length")
   end subroutine test_dense_complex_rhs

   !> Complex symmetric operator: zsytrf, both orientations, reuse and reset
   subroutine test_dense_complex_symmetric(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(complex_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      complex(wp) :: a(n, n), x(n), xt(n), b(n)

      call complex_symmetric_matrix(a)
      call check(error, maxval(abs(a - transpose(a))) <= 0.0_wp, "matrix is symmetric")
      if (allocated(error)) return
      call check(error, maxval(abs(a - conjg(transpose(a)))) > 1.0_wp, "matrix is not Hermitian")
      if (allocated(error)) return
      call new_complex_test_operator(op, a, 1, .true.)
      call solver%reset(1)
      call check(error, solver%factorization(1), 0, "no factors after reset")
      if (allocated(error)) return
      call solver%setup_complex(op, err)
      call check_moist_error(error, err, "setup_complex")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_complex_symmetric, "factorization kind")
      if (allocated(error)) return
      call solver%setup_complex(op, err)
      call check_moist_error(error, err, "second setup_complex")
      if (allocated(error)) return

      call check_complex_orientations(error, solver, op, thr, x, xt)
      if (allocated(error)) return
      call check(error, znorm(x - xt) <= thr*znorm(x), "symmetric: transpose solve equals normal")
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 2, "one factorization per setup_complex")
      if (allocated(error)) return
      call check(error, solver%n_solves, 2, "one solve per complex right-hand side")
      if (allocated(error)) return

      ! Reset discards the complex factors; the next solve factors again
      call solver%reset(1)
      call check(error, solver%factorization(1), 0, "factors discarded by reset")
      if (allocated(error)) return
      call complex_rhs(b)
      call solver%solve_complex(op, .false., b, x, err)
      call check_moist_error(error, err, "solve after reset")
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 3, "refactored after reset")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_complex_symmetric, &
         & "factors after lazy setup_complex")
   end subroutine test_dense_complex_symmetric

   !> Complex general operator: zgetrf; the transpose is plain, not conjugated
   subroutine test_dense_complex_lu(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(complex_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      complex(wp) :: a(n, n), x(n), xt(n)

      call complex_general_matrix(a)
      call new_complex_test_operator(op, a, 1, .false.)
      call solver%reset(1)
      call solver%setup_complex(op, err)
      call check_moist_error(error, err, "setup_complex")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_complex_lu, "factorization kind")
      if (allocated(error)) return

      call check_complex_orientations(error, solver, op, thr, x, xt)
      if (allocated(error)) return
      call check(error, znorm(x - xt) > 1.0e-3_wp, "transpose solve differs from normal solve")
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 1, "one factorization")
      if (allocated(error)) return
      call check(error, solver%n_solves, 2, "one solve per complex right-hand side")
   end subroutine test_dense_complex_lu

   !> Complex operators: invalid ids, shapes, singular matrices, family clashes
   subroutine test_dense_complex_invalid_calls(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(complex_test_operator) :: zop
      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      complex(wp) :: za(n, n), b(n), x(n)
      real(wp) :: a(n, n), rb(n), rx(n)

      call complex_general_matrix(za)
      call new_complex_test_operator(zop, za, 2, .false.)
      call complex_rhs(b)

      call solver%setup_complex(zop, err)
      call check(error, allocated(err), "setup_complex before reset")
      if (allocated(error)) return
      call solver%reset(1)
      call solver%setup_complex(zop, err)
      call check(error, allocated(err), "id outside the reset range")
      if (allocated(error)) return
      call check(error, index(err%message, "[pcm_dense_setup_complex]") == 1, &
         & "error names the routine")
      if (allocated(error)) return
      call solver%solve_complex(zop, .false., b, x, err)
      call check(error, allocated(err), "solve with an id outside the reset range")
      if (allocated(error)) return

      call solver%reset(2)
      call solver%solve_complex(zop, .false., b(:n - 1), x, err)
      call check(error, allocated(err), "right-hand side of the wrong length")
      if (allocated(error)) return
      call solver%solve_complex(zop, .false., b, x(:n - 1), err)
      call check(error, allocated(err), "solution of the wrong length")
      if (allocated(error)) return
      call check(error, solver%n_solves, 0, "no solve counted")
      if (allocated(error)) return

      ! Singular matrices, symmetric and general
      za = (0.0_wp, 0.0_wp)
      call new_complex_test_operator(zop, za, 1, .true.)
      call solver%setup_complex(zop, err)
      call check(error, allocated(err), "singular complex symmetric operator")
      if (allocated(error)) return
      call check(error, index(err%message, "zsytrf") > 0, "error names zsytrf")
      if (allocated(error)) return
      call new_complex_test_operator(zop, za, 1, .false.)
      call solver%solve_complex(zop, .false., b, x, err)
      call check(error, allocated(err), "singular complex general operator")
      if (allocated(error)) return
      call check(error, index(err%message, "zgetrf") > 0, "error names zgetrf")
      if (allocated(error)) return
      call check(error, solver%factorization(1), 0, "operator stays unfactored")
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 0, "no factorization counted")
      if (allocated(error)) return

      ! One id for a real and a complex operator in the same phase
      call general_matrix(a)
      call new_test_operator(op, a, 1, .false., .false.)
      call complex_general_matrix(za)
      call new_complex_test_operator(zop, za, 1, .false.)
      call solver%setup(op, err)
      call check_moist_error(error, err, "real setup")
      if (allocated(error)) return
      call solver%solve_complex(zop, .false., b, x, err)
      call check(error, allocated(err), "complex solve on an id holding real factors")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_lu, "real factors kept")
      if (allocated(error)) return

      ! An explicit setup_complex replaces the real factors of the id
      call solver%setup_complex(zop, err)
      call check_moist_error(error, err, "complex setup")
      if (allocated(error)) return
      rb = real(b, wp)
      call solver%solve(op, .false., rb, rx, err)
      call check(error, allocated(err), "real solve on an id holding complex factors")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_complex_lu, "complex factors kept")
   end subroutine test_dense_complex_invalid_calls

   !> Defaults of an extension overriding only solve_block
   !> - Real operator, complex right-hand side: parts go through solve_block
   !> - setup_complex does nothing; complex operators are an error
   subroutine test_solver_complex_defaults(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(complex_test_operator) :: zop
      type(scaling_test_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n)
      complex(wp) :: za(n, n), b(n), x(n)

      call general_matrix(a)
      call new_test_operator(op, a, 1, .false., .false.)
      call complex_rhs(b)
      call solver%solve_complex(op, .false., b, x, err)
      call check_moist_error(error, err, "complex right-hand side")
      if (allocated(error)) return
      call check(error, maxval(abs(x - 2.0_wp*b)) <= 0.0_wp, "normal solve through solve_block")
      if (allocated(error)) return
      call solver%solve_complex(op, .true., b, x, err)
      call check_moist_error(error, err, "complex right-hand side, transpose")
      if (allocated(error)) return
      call check(error, maxval(abs(x - 3.0_wp*b)) <= 0.0_wp, "transpose passed to solve_block")
      if (allocated(error)) return

      call complex_general_matrix(za)
      call new_complex_test_operator(zop, za, 1, .false.)
      call solver%setup_complex(zop, err)
      call check_moist_error(error, err, "default setup_complex")
      if (allocated(error)) return
      call solver%solve_complex(zop, .false., b, x, err)
      call check(error, allocated(err), "default refuses complex operators")
      if (allocated(error)) return
      call check(error, index(err%message, "[pcm_solver]") == 1, "error names pcm_solver")
   end subroutine test_solver_complex_defaults

   !> A second setup of the same id (repeated medium phase) factors the new matrix
   !> - Solves after the second setup use the new matrix
   !> - A failed setup discards the previous factors
   subroutine test_dense_setup_replaces(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a1(n, n), a2(n, n), b(n, 1), x1(n), x2(n)

      call spd_matrix(a1)
      call general_matrix(a2)
      call rhs_block(b)
      call solver%reset(1)

      call new_test_operator(op, a1, 1, .true., .true.)
      call solver%setup(op, err)
      call check_moist_error(error, err, "first setup")
      if (allocated(error)) return
      call solver%solve(op, .false., b(:, 1), x1, err)
      call check_moist_error(error, err, "first solve")
      if (allocated(error)) return
      call check(error, rel_residual(a1, x1, b(:, 1)) < thr, "first residual")
      if (allocated(error)) return

      call new_test_operator(op, a2, 1, .false., .false.)
      call solver%setup(op, err)
      call check_moist_error(error, err, "second setup")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_lu, "factorization of the new operator")
      if (allocated(error)) return
      call solver%solve(op, .false., b(:, 1), x2, err)
      call check_moist_error(error, err, "second solve")
      if (allocated(error)) return
      call check(error, rel_residual(a2, x2, b(:, 1)) < thr, "second solve uses the new matrix")
      if (allocated(error)) return
      call check(error, norm2(x2 - x1) > 1.0e-3_wp, "second solution differs from the first")
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 2, "both setups factor")
      if (allocated(error)) return
      call check(error, solver%n_solves, 2, "solves counted")
      if (allocated(error)) return

      ! Declared positive definite but indefinite: the old factors are gone
      call indefinite_matrix(a2)
      call new_test_operator(op, a2, 1, .true., .true.)
      call solver%setup(op, err)
      call check(error, allocated(err), "failed setup reports an error")
      if (allocated(error)) return
      call check(error, solver%factorization(1), 0, "failed setup discards the old factors")
   end subroutine test_dense_setup_replaces

   !> A second setup_complex of the same id factors the new complex matrix
   subroutine test_dense_setup_complex_replaces(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(complex_test_operator) :: op
      type(pcm_dense_solver) :: solver
      type(moist_error_type), allocatable :: err
      complex(wp) :: a1(n, n), a2(n, n), b(n), x1(n), x2(n), y(n)

      call complex_symmetric_matrix(a1)
      call complex_general_matrix(a2)
      call complex_rhs(b)
      call solver%reset(1)

      call new_complex_test_operator(op, a1, 1, .true.)
      call solver%setup_complex(op, err)
      call check_moist_error(error, err, "first setup_complex")
      if (allocated(error)) return
      call solver%solve_complex(op, .true., b, x1, err)
      call check_moist_error(error, err, "first solve")
      if (allocated(error)) return
      call op%apply(.true., x1, y)
      call check(error, znorm(y - b) < thr*znorm(b), "first residual")
      if (allocated(error)) return

      call new_complex_test_operator(op, a2, 1, .false.)
      call solver%setup_complex(op, err)
      call check_moist_error(error, err, "second setup_complex")
      if (allocated(error)) return
      call check(error, solver%factorization(1), pcm_factor_complex_lu, &
         & "factorization of the new operator")
      if (allocated(error)) return
      call solver%solve_complex(op, .true., b, x2, err)
      call check_moist_error(error, err, "second solve")
      if (allocated(error)) return
      call op%apply(.true., x2, y)
      call check(error, znorm(y - b) < thr*znorm(b), "second solve uses the new matrix")
      if (allocated(error)) return
      call check(error, znorm(x2 - x1) > 1.0e-3_wp, "second solution differs from the first")
      if (allocated(error)) return
      call check(error, solver%n_factorizations, 2, "both setups factor")
      if (allocated(error)) return

      ! Singular replacement: the old factors are gone
      a2 = (0.0_wp, 0.0_wp)
      call new_complex_test_operator(op, a2, 1, .false.)
      call solver%setup_complex(op, err)
      call check(error, allocated(err), "failed setup_complex reports an error")
      if (allocated(error)) return
      call check(error, solver%factorization(1), 0, "failed setup_complex discards the old factors")
   end subroutine test_dense_setup_complex_replaces

   !> Pseudo-random numbers in [-1, 1) from the MINSTD generator
   !>
   !> @param[in,out] seed  generator state in 1..2147483646
   !> @param[out]    v     random numbers
   pure subroutine minstd_fill(seed, v)
      !> Generator state
      integer(int64), intent(inout) :: seed
      !> Random numbers
      real(wp), intent(out) :: v(:)

      integer(int64), parameter :: modulus = 2147483647_int64
      integer :: i

      do i = 1, size(v)
         seed = mod(48271_int64*seed, modulus)
         v(i) = 2.0_wp*real(seed, wp)/real(modulus, wp) - 1.0_wp
      end do
   end subroutine minstd_fill

   !> Well-conditioned random matrix: diagonal 2..3 plus a random part of norm ~0.3
   !>
   !> - The symmetric part is positive definite, so restarted GMRES converges
   !>
   !> @param[in]  seed       generator seed
   !> @param[in]  symmetric  symmetrize the random part
   !> @param[out] a          matrix, (nr, nr)
   pure subroutine random_matrix(seed, symmetric, a)
      !> Generator seed
      integer, intent(in) :: seed
      !> Symmetrize the random part
      logical, intent(in) :: symmetric
      !> Matrix
      real(wp), intent(out) :: a(nr, nr)

      integer(int64) :: state
      integer :: i

      state = int(seed, int64)
      call minstd_fill(state, a(:, 1))
      do i = 2, nr
         call minstd_fill(state, a(:, i))
      end do
      a = 0.5_wp/sqrt(real(nr, wp))*a
      if (symmetric) a = 0.5_wp*(a + transpose(a))
      do i = 1, nr
         a(i, i) = a(i, i) + 2.0_wp + real(i, wp)/real(nr, wp)
      end do
   end subroutine random_matrix

   !> Well-conditioned random complex matrix, symmetric (not Hermitian) on request
   !>
   !> @param[in]  seed       generator seed
   !> @param[in]  symmetric  make the matrix equal its plain transpose
   !> @param[out] a          matrix, (nr, nr)
   pure subroutine random_complex_matrix(seed, symmetric, a)
      !> Generator seed
      integer, intent(in) :: seed
      !> Make the matrix equal its plain transpose
      logical, intent(in) :: symmetric
      !> Matrix
      complex(wp), intent(out) :: a(nr, nr)

      real(wp) :: re(nr, nr), im(nr, nr)

      call random_matrix(seed, symmetric, re)
      call random_matrix(seed + 7919, symmetric, im)
      a = cmplx(re, 0.5_wp*im, kind=wp)
   end subroutine random_complex_matrix

   !> Right-hand sides with distinct columns, (nr, nrhs)
   !>
   !> @param[in]  seed  generator seed
   !> @param[out] b     right-hand sides
   pure subroutine random_rhs(seed, b)
      !> Generator seed
      integer, intent(in) :: seed
      !> Right-hand sides
      real(wp), intent(out) :: b(:, :)

      integer(int64) :: state
      integer :: k

      state = int(seed, int64)
      do k = 1, size(b, 2)
         call minstd_fill(state, b(:, k))
      end do
   end subroutine random_rhs

   !> GMRES solver with a block-Jacobi preconditioner of a given block size
   !>
   !> @param[out] solver      solver, reset for n_operators
   !> @param[in]  block_size  points per cluster
   !> @param[in]  tol         relative residual threshold
   !> @param[in]  n_operators number of operator slots
   subroutine new_gmres(solver, block_size, tol, n_operators)
      !> Solver
      type(pcm_gmres_solver), intent(out) :: solver
      !> Points per cluster
      integer, intent(in) :: block_size
      !> Relative residual threshold
      real(wp), intent(in) :: tol
      !> Number of operator slots
      integer, intent(in) :: n_operators

      type(pcm_block_jacobi) :: jacobi

      jacobi%block_size = block_size
      allocate (solver%preconditioner, source=jacobi)
      solver%tol = tol
      call solver%reset(n_operators)
   end subroutine new_gmres

   !> Two-field block operator on 4 separated groups of 8 points
   !>
   !> - Points of group g lie near x = 100 g; point p belongs to group
   !>   mod(p - 1, 4), so contiguous index ranges mix the groups
   !> - Points couple only within their group (exp(-distance)), the fields
   !>   through a non-symmetric 2 x 2 weight: the matrix is block diagonal
   !>   over the groups
   !>
   !> @param[out] xyz  positions, (3, 32)
   !> @param[out] a    matrix, (64, 64)
   pure subroutine grouped_block_matrix(xyz, a)
      !> Positions
      real(wp), intent(out) :: xyz(3, 32)
      !> Matrix
      real(wp), intent(out) :: a(64, 64)

      real(wp), parameter :: weight(2, 2) = reshape([1.0_wp, -0.3_wp, 0.4_wp, 0.8_wp], [2, 2])
      integer(int64) :: state
      real(wp) :: coupling
      integer :: i, j, fa, fb

      state = 12345_int64
      call minstd_fill(state, xyz(:, 1))
      do i = 2, 32
         call minstd_fill(state, xyz(:, i))
      end do
      do i = 1, 32
         xyz(1, i) = xyz(1, i) + 100.0_wp*real(mod(i - 1, 4), wp)
      end do
      a = 0.0_wp
      do j = 1, 32
         do i = 1, 32
            if (mod(i - 1, 4) /= mod(j - 1, 4)) cycle
            coupling = 0.3_wp*exp(-norm2(xyz(:, i) - xyz(:, j)))
            do fb = 1, 2
               do fa = 1, 2
                  a((fa - 1)*32 + i, (fb - 1)*32 + j) = weight(fa, fb)*coupling
               end do
            end do
         end do
      end do
      do i = 1, 64
         a(i, i) = a(i, i) + 3.0_wp + real(i, wp)/64.0_wp
      end do
   end subroutine grouped_block_matrix

   !> inexact is false for direct solvers and true for GMRES
   subroutine test_inexact(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(pcm_dense_solver) :: dense
      type(pcm_gmres_solver) :: gmres
      type(scaling_test_solver) :: scaling

      call check(error, dense%inexact(), .false., "dense solver is exact")
      if (allocated(error)) return
      call check(error, gmres%inexact(), "GMRES is inexact")
      if (allocated(error)) return
      call check(error, scaling%inexact(), .false., "pcm_solver default is exact")
   end subroutine test_inexact

   !> Matrix-free kernels: never for direct solvers, by explicit storage for GMRES
   subroutine test_matrix_free(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(pcm_dense_solver) :: dense
      type(pcm_gmres_solver) :: gmres

      call check(error, dense%matrix_free(100000, 4), .false., "dense solver stores matrices")
      if (allocated(error)) return
      call check(error, gmres%storage, pcm_storage_dense, "default dense storage")
      if (allocated(error)) return
      call check(error, gmres%matrix_free(600000000, 4), .false., "default stores matrices at any size")
      if (allocated(error)) return
      gmres%storage = pcm_storage_dense
      call check(error, gmres%matrix_free(100000, 4), .false., "dense storage")
      if (allocated(error)) return
      gmres%storage = pcm_storage_matrix_free
      call check(error, gmres%matrix_free(10, 1), "matrix-free storage")
   end subroutine test_matrix_free

   !> Default settings; the default preconditioner is created by the first setup
   subroutine test_gmres_defaults(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(pcm_gmres_solver) :: solver
      type(pcm_block_jacobi) :: jacobi
      type(dense_test_operator) :: op
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n)

      call check(error, solver%tol, 1.0e-8_wp, "default tolerance")
      if (allocated(error)) return
      call check(error, solver%restart, 100, "default restart")
      if (allocated(error)) return
      call check(error, solver%max_iter, 1000, "default iteration limit")
      if (allocated(error)) return
      call check(error, jacobi%block_size, 64, "default block size")
      if (allocated(error)) return
      call check(error, allocated(solver%preconditioner), .false., "no preconditioner yet")
      if (allocated(error)) return

      call general_matrix(a)
      call new_test_operator(op, a, 1, .false., .false.)
      call solver%reset(1)
      call solver%setup(op, err)
      call check_moist_error(error, err, "setup")
      if (allocated(error)) return
      call check(error, allocated(solver%preconditioner), "preconditioner created by setup")
      if (allocated(error)) return
      select type (p => solver%preconditioner)
      type is (pcm_block_jacobi)
         call check(error, p%n_blocks(1), 1, "one block for 6 points")
      class default
         call test_failed(error, "default preconditioner is not pcm_block_jacobi")
      end select
   end subroutine test_gmres_defaults

   !> Real operators: GMRES equals the dense solution in both orientations
   !> - Non-symmetric and symmetric operators, three right-hand sides each
   !> - Block size 2, so GMRES iterates
   subroutine test_gmres_real_vs_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nrhs = 3
      real(wp), parameter :: tol = 1.0e-12_wp
      type(dense_test_operator) :: op
      type(pcm_gmres_solver) :: solver
      type(pcm_dense_solver) :: direct
      type(moist_error_type), allocatable :: err
      real(wp) :: a(nr, nr), b(nr, nrhs), x(nr, nrhs), x_ref(nr, nrhs), at(nr, nr)
      integer :: isym, it, k, n_iterations
      logical :: trans

      call random_rhs(11, b)
      do isym = 1, 2
         call random_matrix(3, isym == 2, a)
         call new_test_operator(op, a, 1, isym == 2, .false.)
         call new_gmres(solver, 2, tol, 1)
         call direct%reset(1)
         do it = 1, 2
            trans = it == 2
            x = 0.0_wp
            n_iterations = solver%n_iterations
            call solver%solve_block(op, trans, b, x, err)
            call check_moist_error(error, err, "GMRES solve_block")
            if (allocated(error)) return
            call check(error, solver%last_iterations > 1, "GMRES iterates")
            if (allocated(error)) return
            call check(error, solver%n_iterations - n_iterations >= nrhs*2, &
               & "iterations are counted")
            if (allocated(error)) return
            x_ref = 0.0_wp
            call direct%solve_block(op, trans, b, x_ref, err)
            call check_moist_error(error, err, "dense solve_block")
            if (allocated(error)) return
            at = a
            if (trans) at = transpose(a)
            do k = 1, nrhs
               call check(error, rel_residual(at, x(:, k), b(:, k)) <= 2.0_wp*tol, &
                  & "GMRES residual")
               if (allocated(error)) return
               call check(error, norm2(x(:, k) - x_ref(:, k)) <= 1.0e-10_wp*norm2(x_ref(:, k)), &
                  & "GMRES equals dense, symmetric = "//merge("t", "f", isym == 2) &
                  & //", trans = "//merge("t", "f", trans))
               if (allocated(error)) return
            end do
         end do
         call check(error, solver%n_solves, 2*nrhs, "right-hand sides are counted")
         if (allocated(error)) return
         call check(error, solver%last_residual <= tol, "last residual")
         if (allocated(error)) return
      end do
   end subroutine test_gmres_real_vs_dense

   !> Complex operators: complex GMRES equals the dense solution, plain transpose
   subroutine test_gmres_complex_vs_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: tol = 1.0e-12_wp
      type(complex_test_operator) :: op
      type(pcm_gmres_solver) :: solver
      type(pcm_dense_solver) :: direct
      type(moist_error_type), allocatable :: err
      complex(wp) :: a(nr, nr), b(nr), x(nr), x_ref(nr), y(nr)
      real(wp) :: re(nr, 2)
      integer :: isym, it
      logical :: trans

      call random_rhs(17, re)
      b = cmplx(re(:, 1), re(:, 2), kind=wp)
      do isym = 1, 2
         call random_complex_matrix(5, isym == 2, a)
         call new_complex_test_operator(op, a, 1, isym == 2)
         call new_gmres(solver, 2, tol, 1)
         call direct%reset(1)
         do it = 1, 2
            trans = it == 2
            x = (0.0_wp, 0.0_wp)
            call solver%solve_complex(op, trans, b, x, err)
            call check_moist_error(error, err, "complex GMRES solve")
            if (allocated(error)) return
            call check(error, solver%last_iterations > 1, "complex GMRES iterates")
            if (allocated(error)) return
            call op%apply(trans, x, y)
            call check(error, znorm(y - b) <= 2.0_wp*tol*znorm(b), "complex GMRES residual")
            if (allocated(error)) return
            x_ref = (0.0_wp, 0.0_wp)
            call direct%solve_complex(op, trans, b, x_ref, err)
            call check_moist_error(error, err, "dense complex solve")
            if (allocated(error)) return
            call check(error, znorm(x - x_ref) <= 1.0e-10_wp*znorm(x_ref), &
               & "complex GMRES equals dense, symmetric = "//merge("t", "f", isym == 2) &
               & //", trans = "//merge("t", "f", trans))
            if (allocated(error)) return
         end do
         call check(error, solver%n_solves, 2, "complex right-hand sides are counted")
         if (allocated(error)) return
      end do
   end subroutine test_gmres_complex_vs_dense

   !> Real operator, complex right-hand side through the default solve_complex
   subroutine test_gmres_complex_rhs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_gmres_solver) :: solver
      real(wp) :: a(n, n)

      call general_matrix(a)
      call new_test_operator(op, a, 1, .false., .false.)
      call new_gmres(solver, 1, 1.0e-13_wp, 1)
      call check_complex_rhs(error, solver, op, 1.0e-10_wp)
   end subroutine test_gmres_complex_rhs

   !> Warm start: the exact solution needs no iteration, a nearby guess fewer
   subroutine test_gmres_warm_start(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(complex_test_operator) :: zop
      type(pcm_gmres_solver) :: solver
      type(pcm_dense_solver) :: direct
      type(moist_error_type), allocatable :: err
      real(wp) :: a(nr, nr), b(nr, 1), x(nr), x_ref(nr), dx(nr, 1)
      complex(wp) :: za(nr, nr), zb(nr), zx(nr), zx_ref(nr)
      integer :: its_zero, its_near

      call random_matrix(23, .false., a)
      call new_test_operator(op, a, 1, .false., .false.)
      call random_rhs(29, b)
      call new_gmres(solver, 1, 1.0e-10_wp, 1)
      call direct%reset(1)
      x_ref = 0.0_wp
      call direct%solve(op, .true., b(:, 1), x_ref, err)
      call check_moist_error(error, err, "dense solve")
      if (allocated(error)) return

      x = x_ref
      call solver%solve(op, .true., b(:, 1), x, err)
      call check_moist_error(error, err, "solve from the solution")
      if (allocated(error)) return
      call check(error, solver%last_iterations, 0, "no iteration from the solution")
      if (allocated(error)) return
      call check(error, norm2(x - x_ref) <= 0.0_wp, "solution unchanged")
      if (allocated(error)) return

      x = 0.0_wp
      call solver%solve(op, .true., b(:, 1), x, err)
      call check_moist_error(error, err, "solve from zero")
      if (allocated(error)) return
      its_zero = solver%last_iterations

      call random_rhs(31, dx)
      x = x_ref + 1.0e-6_wp*norm2(x_ref)*dx(:, 1)/norm2(dx(:, 1))
      call solver%solve(op, .true., b(:, 1), x, err)
      call check_moist_error(error, err, "solve from a nearby guess")
      if (allocated(error)) return
      its_near = solver%last_iterations
      call check(error, its_near > 0 .and. its_near < its_zero, &
         & "nearby guess needs fewer iterations than zero")
      if (allocated(error)) return
      call check(error, rel_residual(transpose(a), x, b(:, 1)) <= 2.0_wp*solver%tol, &
         & "residual from a nearby guess")
      if (allocated(error)) return

      ! Complex operator from its solution
      call random_complex_matrix(37, .false., za)
      call new_complex_test_operator(zop, za, 1, .false.)
      zb = cmplx(b(:, 1), dx(:, 1), kind=wp)
      zx_ref = (0.0_wp, 0.0_wp)
      call direct%reset(1)
      call direct%solve_complex(zop, .false., zb, zx_ref, err)
      call check_moist_error(error, err, "dense complex solve")
      if (allocated(error)) return
      zx = zx_ref
      call solver%solve_complex(zop, .false., zb, zx, err)
      call check_moist_error(error, err, "complex solve from the solution")
      if (allocated(error)) return
      call check(error, solver%last_iterations, 0, "no complex iteration from the solution")
   end subroutine test_gmres_warm_start

   !> A zero right-hand side gives x = 0 without iterating
   subroutine test_gmres_zero_rhs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(complex_test_operator) :: zop
      type(pcm_gmres_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(nr, nr), b(nr), x(nr)
      complex(wp) :: za(nr, nr), zb(nr), zx(nr)

      call random_matrix(41, .false., a)
      call new_test_operator(op, a, 1, .false., .false.)
      call new_gmres(solver, 4, 1.0e-10_wp, 1)
      b = 0.0_wp
      x = 1.0_wp
      call solver%solve(op, .false., b, x, err)
      call check_moist_error(error, err, "zero right-hand side")
      if (allocated(error)) return
      call check(error, maxval(abs(x)) <= 0.0_wp, "zero solution")
      if (allocated(error)) return
      call check(error, solver%last_iterations, 0, "no iteration")
      if (allocated(error)) return

      call random_complex_matrix(43, .true., za)
      call new_complex_test_operator(zop, za, 1, .true.)
      call solver%reset(1)
      zb = (0.0_wp, 0.0_wp)
      zx = (1.0_wp, -1.0_wp)
      call solver%solve_complex(zop, .true., zb, zx, err)
      call check_moist_error(error, err, "zero complex right-hand side")
      if (allocated(error)) return
      call check(error, maxval(abs(zx)) <= 0.0_wp, "zero complex solution")
   end subroutine test_gmres_zero_rhs

   !> Non-convergence within max_iter is an error that names the iterations
   subroutine test_gmres_no_convergence(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_gmres_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(nr, nr), b(nr, 1), x(nr)

      call random_matrix(47, .false., a)
      call new_test_operator(op, a, 3, .false., .false.)
      call random_rhs(53, b)
      call new_gmres(solver, 1, 1.0e-12_wp, 3)
      solver%max_iter = 3
      x = 0.0_wp
      call solver%solve(op, .false., b(:, 1), x, err)
      call check(error, allocated(err), "three iterations do not converge")
      if (allocated(error)) return
      call check(error, index(err%message, "[pcm_gmres_solve_block] operator 3") == 1, &
         & "error names the routine and operator")
      if (allocated(error)) return
      call check(error, index(err%message, "did not converge in 3 iterations") > 0, &
         & "error reports the iterations")
      if (allocated(error)) return
      call check(error, solver%last_iterations, 3, "iterations recorded")
      if (allocated(error)) return
      call check(error, solver%last_residual > solver%tol, "residual recorded")
      if (allocated(error)) return

      ! max_iter 0 only accepts a converged initial guess
      solver%max_iter = 0
      x = 0.0_wp
      call solver%solve(op, .false., b(:, 1), x, err)
      call check(error, allocated(err), "max_iter 0 from zero")
   end subroutine test_gmres_no_convergence

   !> A restart shorter than the iterations needed still converges
   subroutine test_gmres_short_restart(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_gmres_solver) :: solver, reference
      type(moist_error_type), allocatable :: err
      real(wp) :: a(nr, nr), b(nr, 1), x(nr), x_ref(nr)

      call random_matrix(59, .true., a)
      call new_test_operator(op, a, 1, .true., .false.)
      call random_rhs(61, b)
      call new_gmres(solver, 1, 1.0e-11_wp, 1)
      solver%restart = 2
      call new_gmres(reference, 1, 1.0e-11_wp, 1)
      x = 0.0_wp
      call solver%solve(op, .false., b(:, 1), x, err)
      call check_moist_error(error, err, "GMRES(2)")
      if (allocated(error)) return
      x_ref = 0.0_wp
      call reference%solve(op, .false., b(:, 1), x_ref, err)
      call check_moist_error(error, err, "GMRES(100)")
      if (allocated(error)) return
      call check(error, solver%last_iterations > 2*solver%restart, "several cycles")
      if (allocated(error)) return
      call check(error, reference%last_iterations < solver%last_iterations, &
         & "restarts cost iterations")
      if (allocated(error)) return
      call check(error, rel_residual(a, x, b(:, 1)) <= 2.0_wp*solver%tol, "GMRES(2) residual")
      if (allocated(error)) return
      call check(error, norm2(x - x_ref) <= 1.0e-9_wp*norm2(x_ref), "GMRES(2) equals GMRES(100)")
   end subroutine test_gmres_short_restart

   !> setup rebuilds the preconditioner when the operator of an id changes
   !>
   !> - Block size >= npt: the preconditioner is the exact inverse, so a
   !>   current preconditioner converges in one iteration and a stale one not
   subroutine test_gmres_setup_rebuilds(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op1, op2
      type(complex_test_operator) :: zop1, zop2
      type(pcm_gmres_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(nr, nr), b(nr, 1), x(nr)
      complex(wp) :: za(nr, nr), zb(nr), zx(nr)

      call random_rhs(67, b)
      call random_matrix(71, .false., a)
      call new_test_operator(op1, a, 1, .false., .false.)
      call random_matrix(73, .false., a)
      call new_test_operator(op2, a, 1, .false., .false.)
      call new_gmres(solver, nr, 1.0e-10_wp, 1)

      call solver%setup(op1, err)
      call check_moist_error(error, err, "first setup")
      if (allocated(error)) return
      x = 0.0_wp
      call solver%solve(op1, .false., b(:, 1), x, err)
      call check_moist_error(error, err, "first solve")
      if (allocated(error)) return
      call check(error, solver%last_iterations, 1, "exact preconditioner: one iteration")
      if (allocated(error)) return

      ! Stale preconditioner of op1 for op2: converges, but not in one iteration
      x = 0.0_wp
      call solver%solve(op2, .false., b(:, 1), x, err)
      call check_moist_error(error, err, "solve with a stale preconditioner")
      if (allocated(error)) return
      call check(error, solver%last_iterations > 1, "stale preconditioner iterates")
      if (allocated(error)) return

      call solver%setup(op2, err)
      call check_moist_error(error, err, "second setup")
      if (allocated(error)) return
      x = 0.0_wp
      call solver%solve(op2, .true., b(:, 1), x, err)
      call check_moist_error(error, err, "second solve")
      if (allocated(error)) return
      call check(error, solver%last_iterations, 1, "setup rebuilt the preconditioner")
      if (allocated(error)) return

      ! After reset a solve sets the preconditioner up again
      call solver%reset(1)
      x = 0.0_wp
      call solver%solve(op1, .false., b(:, 1), x, err)
      call check_moist_error(error, err, "solve after reset")
      if (allocated(error)) return
      call check(error, solver%last_iterations, 1, "lazy setup after reset")
      if (allocated(error)) return

      ! Complex operators: setup_complex rebuilds as well
      call random_complex_matrix(79, .false., za)
      call new_complex_test_operator(zop1, za, 1, .false.)
      call random_complex_matrix(83, .false., za)
      call new_complex_test_operator(zop2, za, 1, .false.)
      zb = cmplx(b(:, 1), 1.0_wp, kind=wp)
      call solver%setup_complex(zop1, err)
      call check_moist_error(error, err, "first setup_complex")
      if (allocated(error)) return
      call solver%setup_complex(zop2, err)
      call check_moist_error(error, err, "second setup_complex")
      if (allocated(error)) return
      zx = (0.0_wp, 0.0_wp)
      call solver%solve_complex(zop2, .true., zb, zx, err)
      call check_moist_error(error, err, "complex solve")
      if (allocated(error)) return
      call check(error, solver%last_iterations, 1, "setup_complex rebuilt the preconditioner")
      if (allocated(error)) return

      ! The id now holds a complex preconditioner; a real solve sets up anew
      x = 0.0_wp
      call solver%solve(op1, .false., b(:, 1), x, err)
      call check_moist_error(error, err, "real solve on a complex slot")
      if (allocated(error)) return
      call check(error, solver%last_iterations, 1, "real preconditioner set up again")
   end subroutine test_gmres_setup_rebuilds

   !> Missing reset, invalid settings and shape mismatches are errors
   subroutine test_gmres_invalid_calls(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(complex_test_operator) :: zop
      type(pcm_gmres_solver) :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: a(n, n), b(n, 2), x(n, 1), x1(n)
      complex(wp) :: za(n, n), zb(n), zx(n - 1)

      call general_matrix(a)
      call new_test_operator(op, a, 2, .false., .false.)
      call rhs_block(b)
      x1 = 0.0_wp
      call solver%solve(op, .false., b(:, 1), x1, err)
      call check(error, allocated(err), "solve before reset")
      if (allocated(error)) return
      call check(error, index(err%message, "reset the solver first") > 0, "error asks for reset")
      if (allocated(error)) return
      call solver%reset(1)
      call solver%setup(op, err)
      call check(error, allocated(err), "id outside the reset range")
      if (allocated(error)) return

      call solver%reset(2)
      x = 0.0_wp
      call solver%solve_block(op, .false., b, x, err)
      call check(error, allocated(err), "solution block of the wrong shape")
      if (allocated(error)) return
      solver%tol = 0.0_wp
      call solver%solve(op, .false., b(:, 1), x1, err)
      call check(error, allocated(err), "zero tolerance")
      if (allocated(error)) return
      call check(error, index(err%message, "invalid settings") > 0, "error names the settings")
      if (allocated(error)) return
      solver%tol = 1.0e-8_wp
      solver%restart = 0
      call solver%solve(op, .false., b(:, 1), x1, err)
      call check(error, allocated(err), "zero restart")
      if (allocated(error)) return
      solver%restart = 10

      call complex_general_matrix(za)
      call new_complex_test_operator(zop, za, 1, .false.)
      call complex_rhs(zb)
      zx = (0.0_wp, 0.0_wp)
      call solver%solve_complex(zop, .false., zb, zx, err)
      call check(error, allocated(err), "complex solution of the wrong length")
      if (allocated(error)) return
      call check(error, solver%n_solves, 0, "no solve counted")
   end subroutine test_gmres_invalid_calls

   !> Block size 1 is Jacobi: y = x / diag(op), both orientations
   subroutine test_block_jacobi_size_one(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_block_jacobi) :: jacobi
      type(moist_error_type), allocatable :: err
      real(wp) :: a(nr, nr), x(nr, 1), y(nr), d(nr)
      integer :: i, it

      call random_matrix(89, .false., a)
      call new_test_operator(op, a, 1, .false., .false.)
      call random_rhs(97, x)
      do i = 1, nr
         d(i) = a(i, i)
      end do
      jacobi%block_size = 1
      call jacobi%reset(1)
      call jacobi%setup(op, err)
      call check_moist_error(error, err, "setup")
      if (allocated(error)) return
      call check(error, jacobi%n_blocks(1), nr, "one block per point")
      if (allocated(error)) return
      do it = 1, 2
         call jacobi%apply(op, it == 2, x(:, 1), y, err)
         call check_moist_error(error, err, "apply")
         if (allocated(error)) return
         call check(error, maxval(abs(y - x(:, 1)/d)) <= 1.0e-15_wp*maxval(abs(y)), &
            & "Jacobi, trans = "//merge("t", "f", it == 2))
         if (allocated(error)) return
      end do

      ! Applying an id that is not set up is an error
      call jacobi%reset(1)
      call jacobi%apply(op, .false., x(:, 1), y, err)
      call check(error, allocated(err), "apply before setup")
   end subroutine test_block_jacobi_size_one

   !> Clusters of a two-field block operator by bisection and by contiguous ranges
   !>
   !> - With points, bisection recovers the 4 separated groups: the block
   !>   preconditioner is the exact inverse (M^-1 op x = x) and GMRES takes one
   !>   iteration, which also checks the unknowns (a - 1) npoint + idx
   !> - Without points, contiguous ranges mix the groups; GMRES still converges
   !> - Points of the wrong shape and npt not a multiple of npoint are errors
   subroutine test_block_jacobi_clusters(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(dense_test_operator) :: op
      type(pcm_block_jacobi) :: jacobi
      type(pcm_gmres_solver) :: solver
      type(pcm_dense_solver) :: direct
      type(moist_error_type), allocatable :: err
      real(wp), target :: xyz(3, 32)
      real(wp) :: a(64, 64), x(64, 2), y(64), ax(64), x_ref(64)
      integer :: it

      call grouped_block_matrix(xyz, a)
      call new_test_operator(op, a, 1, .false., .false.)
      op%npoint = 32
      op%points => xyz
      call random_rhs(101, x)

      jacobi%block_size = 8
      call jacobi%reset(1)
      call jacobi%setup(op, err)
      call check_moist_error(error, err, "setup with points")
      if (allocated(error)) return
      call check(error, jacobi%n_blocks(1), 4, "four clusters")
      if (allocated(error)) return
      do it = 1, 2
         if (it == 2) then
            ax = matmul(x(:, 1), a)
         else
            ax = matmul(a, x(:, 1))
         end if
         call jacobi%apply(op, it == 2, ax, y, err)
         call check_moist_error(error, err, "apply")
         if (allocated(error)) return
         call check(error, maxval(abs(y - x(:, 1))) <= 1.0e-13_wp*maxval(abs(x(:, 1))), &
            & "clusters are the groups, trans = "//merge("t", "f", it == 2))
         if (allocated(error)) return
      end do

      call new_gmres(solver, 8, 1.0e-10_wp, 1)
      call direct%reset(1)
      x_ref = 0.0_wp
      call direct%solve(op, .true., x(:, 2), x_ref, err)
      call check_moist_error(error, err, "dense solve")
      if (allocated(error)) return
      y = 0.0_wp
      call solver%solve(op, .true., x(:, 2), y, err)
      call check_moist_error(error, err, "GMRES with clusters")
      if (allocated(error)) return
      call check(error, solver%last_iterations, 1, "exact block preconditioner")
      if (allocated(error)) return
      call check(error, norm2(y - x_ref) <= 1.0e-9_wp*norm2(x_ref), "clusters: solution")
      if (allocated(error)) return

      ! Contiguous index ranges mix the groups
      op%points => null()
      call jacobi%setup(op, err)
      call check_moist_error(error, err, "setup without points")
      if (allocated(error)) return
      call check(error, jacobi%n_blocks(1), 4, "four ranges")
      if (allocated(error)) return
      call solver%setup(op, err)
      call check_moist_error(error, err, "GMRES setup without points")
      if (allocated(error)) return
      y = 0.0_wp
      call solver%solve(op, .true., x(:, 2), y, err)
      call check_moist_error(error, err, "GMRES with ranges")
      if (allocated(error)) return
      call check(error, solver%last_iterations > 1, "ranges are not the groups")
      if (allocated(error)) return
      call check(error, norm2(y - x_ref) <= 1.0e-9_wp*norm2(x_ref), "ranges: solution")
      if (allocated(error)) return

      ! Invalid layouts
      op%points => xyz(:, 1:31)
      call jacobi%setup(op, err)
      call check(error, allocated(err), "points of the wrong shape")
      if (allocated(error)) return
      op%points => null()
      op%npoint = 30
      call jacobi%setup(op, err)
      call check(error, allocated(err), "npt not a multiple of npoint")
      if (allocated(error)) return
      call check(error, jacobi%n_blocks(1), 0, "failed setup leaves the id empty")
   end subroutine test_block_jacobi_clusters

   !> Complex block-Jacobi: one block is exact, small blocks converge to the dense solution
   subroutine test_block_jacobi_complex(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(complex_test_operator) :: op
      type(pcm_gmres_solver) :: solver
      type(pcm_dense_solver) :: direct
      type(moist_error_type), allocatable :: err
      complex(wp) :: a(nr, nr), b(nr), x(nr), x_ref(nr)
      real(wp) :: re(nr, 2)
      integer :: it, ib, block_size
      logical :: trans

      call random_complex_matrix(103, .false., a)
      call new_complex_test_operator(op, a, 1, .false.)
      call random_rhs(107, re)
      b = cmplx(re(:, 1), re(:, 2), kind=wp)
      call direct%reset(1)
      do ib = 1, 2
         block_size = merge(nr, 7, ib == 1)
         call new_gmres(solver, block_size, 1.0e-11_wp, 1)
         do it = 1, 2
            trans = it == 2
            x_ref = (0.0_wp, 0.0_wp)
            call direct%solve_complex(op, trans, b, x_ref, err)
            call check_moist_error(error, err, "dense complex solve")
            if (allocated(error)) return
            x = (0.0_wp, 0.0_wp)
            call solver%solve_complex(op, trans, b, x, err)
            call check_moist_error(error, err, "complex GMRES")
            if (allocated(error)) return
            if (ib == 1) then
               call check(error, solver%last_iterations, 1, "one complex block is exact")
               if (allocated(error)) return
            end if
            call check(error, znorm(x - x_ref) <= 1.0e-9_wp*znorm(x_ref), &
               & "complex block-Jacobi GMRES, block size "//merge("nr", "7 ", ib == 1) &
               & //", trans = "//merge("t", "f", trans))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_block_jacobi_complex

   !> Newton converges to the root from x = -2 in a few full steps
   subroutine test_nonlinear_converges(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cubic_test_problem) :: problem
      type(pcm_dense_solver), target :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: x(2)

      problem%npt = 2
      problem%solver => solver
      solver%nonlinear_tol = 1.0e-13_wp
      x = [-2.0_wp, -1.5_wp]
      call problem%solve(x, err)
      call check_moist_error(error, err, "Newton from x = -2")
      if (allocated(error)) return
      call check(error, all(abs(x - cubic_root) <= 1.0e-12_wp), "root of the cubic")
      if (allocated(error)) return
      call check(error, solver%nonlinear_iterations <= 6, "a few Newton steps")
      if (allocated(error)) return
      call check(error, solver%nonlinear_residual <= solver%nonlinear_tol, "residual recorded")
   end subroutine test_nonlinear_converges

   !> Newton from x = 0 stalls at the local minimum of |F| and reports it at once
   subroutine test_nonlinear_stalls(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cubic_test_problem) :: problem
      type(pcm_dense_solver), target :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: x(1)

      problem%npt = 1
      problem%solver => solver
      x = 0.0_wp
      call problem%solve(x, err)
      call check(error, allocated(err), "no root reachable from x = 0")
      if (allocated(error)) return
      call check(error, index(err%message, "[pcm_nonlinear_solve] Newton stalled") == 1, &
         & "error reports the stall")
      if (allocated(error)) return
      call check(error, abs(x(1) - sqrt(2.0_wp/3.0_wp)) <= 1.0e-3_wp, &
         & "iterate at the local minimum of |F|")
      if (allocated(error)) return
      call check(error, solver%nonlinear_iterations < solver%nonlinear_max_iter, &
         & "stops before the iteration limit")
      if (allocated(error)) return
      call check(error, solver%nonlinear_residual > 0.4_wp, "residual recorded")
   end subroutine test_nonlinear_stalls

   !> Too few Newton steps for the tolerance are an error
   subroutine test_nonlinear_max_iter(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cubic_test_problem) :: problem
      type(pcm_dense_solver), target :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: x(1)

      problem%npt = 1
      problem%solver => solver
      solver%nonlinear_max_iter = 2
      x = -2.0_wp
      call problem%solve(x, err)
      call check(error, allocated(err), "two steps do not converge")
      if (allocated(error)) return
      call check(error, index(err%message, "did not converge") > 0 &
         & .and. index(err%message, "after 2 steps") > 0, "error reports the steps")
      if (allocated(error)) return
      call check(error, solver%nonlinear_iterations, 2, "iterations recorded")
   end subroutine test_nonlinear_max_iter

   !> Newton rejects a NaN forcing
   subroutine test_nonlinear_rejects_nan_forcing(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cubic_test_problem) :: problem
      type(pcm_dense_solver), target :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: x(2)

      problem%npt = 2
      problem%solver => solver
      x = 1.0_wp
      problem%offset = ieee_value(1.0_wp, ieee_quiet_nan)
      call problem%solve(x, err)
      call check_moist_error(error, err)
   end subroutine test_nonlinear_rejects_nan_forcing

   !> Newton rejects an infinite forcing
   subroutine test_nonlinear_rejects_inf_forcing(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cubic_test_problem) :: problem
      type(pcm_dense_solver), target :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: x(2)

      problem%npt = 2
      problem%solver => solver
      x = 1.0_wp
      problem%offset = ieee_value(1.0_wp, ieee_positive_inf)
      call problem%solve(x, err)
      call check_moist_error(error, err)
   end subroutine test_nonlinear_rejects_inf_forcing

   !> Exact zero forcing returns x = 0 without error
   subroutine test_nonlinear_zero_forcing(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cubic_test_problem) :: problem
      type(pcm_dense_solver), target :: solver
      type(moist_error_type), allocatable :: err
      real(wp) :: x(2)

      problem%npt = 2
      problem%solver => solver
      x = 1.0_wp
      problem%offset = 0.0_wp
      call problem%solve(x, err)
      call check_moist_error(error, err, "zero forcing")
      if (allocated(error)) return
      call check(error, all(x == 0.0_wp), "zero forcing gives x = 0")
   end subroutine test_nonlinear_zero_forcing
end module test_model_component_pcm_solvers
