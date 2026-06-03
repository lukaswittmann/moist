!> Test suite for linalg utilities in moist_math_linalg
module test_math_linalg
   use mctc_env, only: wp
   use mctc_io_utils, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_linalg, only: mat3x3_inv, setup_tangent_frame, sym3_21, &
                                outer3, outer3_linear, outer_matrix, logaddexp, eig_2x2_symmetric, &
                                outer4, sym4_31, sym4_22, sym4_211
   use moist_math_lapack, only: getrf, getri
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_lapack_syev, only: dsyev
   ! Raw vendored linalg solver APIs (upstream test ports folded into this suite).
   use mctc_env_accuracy, only: ip => i4
   use moist_math_linalg_lusol_ez, only: solve
   use moist_math_linalg_lsmr, only: lsmr
   use moist_math_linalg_lsqr, only: lsqr_solver_ez
   implicit none(type, external)
   private

   public :: collect_math_linalg

   !> Shared tolerance for tensor tests.
   real(wp), parameter :: tensor_tol = 1.0e-12_wp

   ! Raw-kernel test constants (vendored APIs ported from upstream test suites)
   !> Residual tolerance for the COO solve ports (upstream uses 1e-12; relaxed)
   real(wp), parameter :: lusol_thr = 1.0e-10_wp
   real(wp), parameter :: lsqr_thr = 1.0e-10_wp
   !> LSMR Paige-Saunders problem size (nbar reduced from upstream 1000) and
   !> singular-value duplication count
   integer(ip), parameter :: lsmr_nbar = 100, lsmr_nduplc = 40
   !> LSMR relative solution-error tolerance (upstream etol)
   real(wp), parameter :: lsmr_etol = 1.0e-3_wp

contains

   !> Register linalg tests.
   subroutine collect_math_linalg(testsuite)
      !> Registered tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mat3x3_inv_identity", test_mat3x3_inv_identity), &
                  new_unittest("mat3x3_inv_diagonal", test_mat3x3_inv_diagonal), &
                  new_unittest("mat3x3_inv_full", test_mat3x3_inv_full), &
                  new_unittest("mat3x3_inv_singular", test_mat3x3_inv_singular), &
                  new_unittest("mat3x3_inv_lapack", test_mat3x3_inv_vs_lapack), &
                  new_unittest("tangent_frame_orthonormal", test_tangent_frame_orthonormal), &
                  new_unittest("tangent_frame_axes", test_tangent_frame_axes), &
                  new_unittest("sym3_21", test_sym3_21), &
                  new_unittest("outer3", test_outer3), &
                  new_unittest("outer3_linear", test_outer3_linear), &
                  new_unittest("outer_matrix", test_outer_matrix), &
                  new_unittest("logaddexp_stability", test_logaddexp_stability), &
                  new_unittest("logaddexp_values", test_logaddexp_values), &
                  new_unittest("eig_2x2_symmetric", test_eig_2x2_symmetric), &
                  new_unittest("outer4_brute_force", test_outer4_brute_force), &
                  new_unittest("outer4_full_symmetry", test_outer4_full_symmetry), &
                  new_unittest("outer4_outer_matrix", test_outer4_outer_matrix), &
                  new_unittest("outer4_zero_vector", test_outer4_zero_vector), &
                  new_unittest("sym4_31_brute_force", test_sym4_31_brute_force), &
                  new_unittest("sym4_31_symmetric_input", test_sym4_31_symmetric_input), &
                  new_unittest("sym4_31_linearity", test_sym4_31_linearity), &
                  new_unittest("sym4_22_brute_force", test_sym4_22_brute_force), &
                  new_unittest("sym4_22_left_biased", test_sym4_22_left_biased), &
                  new_unittest("sym4_22_jkl_symmetry", test_sym4_22_jkl_symmetry), &
                  new_unittest("sym4_22_full_symmetric", test_sym4_22_full_symmetric), &
                  new_unittest("sym4_22_contraction", test_sym4_22_contraction), &
                  new_unittest("sym4_211_brute_force", test_sym4_211_brute_force), &
                  new_unittest("sym4_211_full_symmetry", test_sym4_211_full_symmetry), &
                  new_unittest("sym4_211_scaling", test_sym4_211_scaling), &
                  ! Raw-kernel ports (vendored sparse linear solvers):
                  new_unittest("lusol-dense-3x3", test_lusol_dense_3x3), &
                  new_unittest("lusol-rectangular-3x4", test_lusol_rectangular_3x4), &
                  new_unittest("lsqr-dense-3x3", test_lsqr_dense_3x3), &
                  new_unittest("lsqr-rectangular-3x4", test_lsqr_rectangular_3x4), &
                  new_unittest("lsmr-over-determined", test_lsmr_over_determined), &
                  new_unittest("lsmr-square", test_lsmr_square), &
                  new_unittest("lsmr-under-determined", test_lsmr_under_determined) &
                  ]
   end subroutine collect_math_linalg

   !> Invert the identity matrix.
   subroutine test_mat3x3_inv_identity(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(3, 3), Ainv(3, 3), det
      logical :: success
      integer :: i

      ! Identity matrix.
      A = 0.0_wp
      do i = 1, 3
         A(i, i) = 1.0_wp
      end do

      success = mat3x3_inv(A, Ainv, det)

      call check(error, abs(det - 1.0_wp) < 1.0e-12_wp, "Determinant should be 1")
      if (allocated(error)) return

      do i = 1, 3
         call check(error, abs(Ainv(i, i) - 1.0_wp) < 1.0e-12_wp, "Diagonal should be 1")
         if (allocated(error)) return
      end do
   end subroutine test_mat3x3_inv_identity

   !> Invert a diagonal matrix.
   subroutine test_mat3x3_inv_diagonal(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(3, 3), Ainv(3, 3), det
      logical :: success

      A = 0.0_wp
      A(1, 1) = 2.0_wp
      A(2, 2) = 3.0_wp
      A(3, 3) = 4.0_wp

      success = mat3x3_inv(A, Ainv, det)

      call check(error, abs(det - 24.0_wp) < 1.0e-12_wp, "Determinant should be 24")
      if (allocated(error)) return

      call check(error, abs(Ainv(1, 1) - 0.5_wp) < 1.0e-12_wp, "Ainv(1,1) should be 0.5")
      if (allocated(error)) return
      call check(error, abs(Ainv(2, 2) - 1.0_wp/3.0_wp) < 1.0e-12_wp, "Ainv(2,2) should be 1/3")
      if (allocated(error)) return
      call check(error, abs(Ainv(3, 3) - 0.25_wp) < 1.0e-12_wp, "Ainv(3,3) should be 0.25")
   end subroutine test_mat3x3_inv_diagonal

   !> Invert a dense matrix.
   subroutine test_mat3x3_inv_full(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(3, 3), Ainv(3, 3), det, prod(3, 3)
      logical :: success
      integer :: i, j

      ! Dense nonsingular matrix.
      A = reshape([ &
                  1.0_wp, 2.0_wp, 3.0_wp, &
                  0.0_wp, 4.0_wp, 5.0_wp, &
                  1.0_wp, 0.0_wp, 6.0_wp &
                  ], [3, 3])

      success = mat3x3_inv(A, Ainv, det)

      ! Check A * Ainv = I.
      prod = matmul(A, Ainv)
      do i = 1, 3
         do j = 1, 3
            if (i == j) then
               call check(error, abs(prod(i, j) - 1.0_wp) < 1.0e-10_wp, &
                          "Product diagonal should be 1")
            else
               call check(error, abs(prod(i, j)) < 1.0e-10_wp, &
                          "Product off-diagonal should be 0")
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_mat3x3_inv_full

   !> Detect a singular matrix.
   subroutine test_mat3x3_inv_singular(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(3, 3), Ainv(3, 3), det
      logical :: success

      ! Two identical rows.
      A = reshape([ &
                  1.0_wp, 2.0_wp, 3.0_wp, &
                  1.0_wp, 2.0_wp, 3.0_wp, &
                  4.0_wp, 5.0_wp, 6.0_wp &
                  ], [3, 3])

      success = mat3x3_inv(A, Ainv, det)

      call check(error, abs(det) < 1.0e-10_wp, "Singular matrix should have zero determinant")
      if (allocated(error)) return

      ! Singular matrices return a zero inverse.
      call check(error, maxval(abs(Ainv)) < 1.0e-10_wp, "Inverse of singular matrix should be zero")
   end subroutine test_mat3x3_inv_singular

   !> Compare inversion with LAPACK.
   subroutine test_mat3x3_inv_vs_lapack(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(3, 3), Ainv_analytic(3, 3), Ainv_lapack(3, 3), det
      logical :: success
      integer :: ipiv(3), info

      ! Symmetric positive definite matrix.
      A = reshape([ &
                  4.0_wp, 1.0_wp, 2.0_wp, &
                  1.0_wp, 5.0_wp, 3.0_wp, &
                  2.0_wp, 3.0_wp, 6.0_wp &
                  ], [3, 3])

      ! Analytic inversion.
      success = mat3x3_inv(A, Ainv_analytic, det)

      ! LAPACK inversion.
      Ainv_lapack = A
      call getrf(Ainv_lapack, ipiv, info)
      if (info /= 0) then
         call check(error, .false., "LAPACK getrf failed")
         return
      end if
      call getri(Ainv_lapack, ipiv, info)
      if (info /= 0) then
         call check(error, .false., "LAPACK getri failed")
         return
      end if

      ! Compare inverses.
      call check(error, maxval(abs(Ainv_analytic - Ainv_lapack)) < 1.0e-10_wp, &
                 "Analytic and LAPACK inversions should match")
   end subroutine test_mat3x3_inv_vs_lapack

   !> Check tangent-frame orthonormality.
   subroutine test_tangent_frame_orthonormal(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: normal(3), t1(3), t2(3)
      real(wp) :: dot_n_t1, dot_n_t2, dot_t1_t2, norm_t1, norm_t2

      ! Generic normal.
      normal = [1.0_wp, 2.0_wp, 3.0_wp]
      normal = normal/sqrt(sum(normal**2))

      call setup_tangent_frame(normal, t1, t2)

      ! Orthogonality.
      dot_n_t1 = dot_product(normal, t1)
      dot_n_t2 = dot_product(normal, t2)
      dot_t1_t2 = dot_product(t1, t2)

      call check(error, abs(dot_n_t1) < 1.0e-12_wp, "t1 should be orthogonal to normal")
      if (allocated(error)) return
      call check(error, abs(dot_n_t2) < 1.0e-12_wp, "t2 should be orthogonal to normal")
      if (allocated(error)) return
      call check(error, abs(dot_t1_t2) < 1.0e-12_wp, "t1 and t2 should be orthogonal")
      if (allocated(error)) return

      ! Normalization.
      norm_t1 = sqrt(sum(t1**2))
      norm_t2 = sqrt(sum(t2**2))
      call check(error, abs(norm_t1 - 1.0_wp) < 1.0e-12_wp, "t1 should be normalized")
      if (allocated(error)) return
      call check(error, abs(norm_t2 - 1.0_wp) < 1.0e-12_wp, "t2 should be normalized")
   end subroutine test_tangent_frame_orthonormal

   !> Check tangent frame for an axis normal.
   subroutine test_tangent_frame_axes(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: normal(3), t1(3), t2(3)

      ! Normal along z.
      normal = [0.0_wp, 0.0_wp, 1.0_wp]
      call setup_tangent_frame(normal, t1, t2)

      ! t1 should be x.
      call check(error, abs(t1(1) - 1.0_wp) < 1.0e-12_wp, "t1(1) should be 1")
      if (allocated(error)) return
      call check(error, abs(t1(2)) < 1.0e-12_wp, "t1(2) should be 0")
      if (allocated(error)) return
      call check(error, abs(t1(3)) < 1.0e-12_wp, "t1(3) should be 0")
   end subroutine test_tangent_frame_axes

   !> Check selected sym3_21 entries.
   subroutine test_sym3_21(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: hess(3, 3), grad(3), tensor(3, 3, 3)
      real(wp) :: expected

      ! Diagonal Hessian.
      hess = reshape([1.0_wp, 0.0_wp, 0.0_wp, &
                      0.0_wp, 2.0_wp, 0.0_wp, &
                      0.0_wp, 0.0_wp, 3.0_wp], [3, 3])
      grad = [1.0_wp, 1.0_wp, 1.0_wp]

      tensor = sym3_21(hess, grad)

      ! Diagonal entry.
      expected = 3.0_wp
      call check(error, abs(tensor(1, 1, 1) - expected) < 1.0e-12_wp, "T_111 should be 3")
      if (allocated(error)) return

      ! Off-diagonal entry.
      call check(error, abs(tensor(1, 2, 3)) < 1.0e-12_wp, "T_123 should be 0")
   end subroutine test_sym3_21

   !> Check selected outer3 entries.
   subroutine test_outer3(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: vec(3), tensor(3, 3, 3)

      vec = [1.0_wp, 2.0_wp, 3.0_wp]
      tensor = outer3(vec)

      ! First diagonal entry.
      call check(error, abs(tensor(1, 1, 1) - 1.0_wp) < 1.0e-12_wp, "T_111 should be 1")
      if (allocated(error)) return

      ! Mixed entry.
      call check(error, abs(tensor(1, 2, 3) - 6.0_wp) < 1.0e-12_wp, "T_123 should be 6")
      if (allocated(error)) return

      ! Second diagonal entry.
      call check(error, abs(tensor(2, 2, 2) - 8.0_wp) < 1.0e-12_wp, "T_222 should be 8")
   end subroutine test_outer3

   !> Check selected outer3_linear entries.
   subroutine test_outer3_linear(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: vec(3), dvec(3), dtensor(3, 3, 3)

      vec = [1.0_wp, 2.0_wp, 3.0_wp]
      dvec = [0.1_wp, 0.2_wp, 0.3_wp]

      dtensor = outer3_linear(vec, dvec)

      ! First diagonal derivative.
      call check(error, abs(dtensor(1, 1, 1) - 0.3_wp) < 1.0e-12_wp, "dT_111 should be 0.3")
      if (allocated(error)) return

      ! Mixed derivative.
      call check(error, abs(dtensor(1, 2, 3) - 1.8_wp) < 1.0e-12_wp, "dT_123 should be 1.8")
   end subroutine test_outer3_linear

   !> Check selected outer_matrix entries.
   subroutine test_outer_matrix(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: left(3), right(3), mat(3, 3)

      left = [1.0_wp, 2.0_wp, 3.0_wp]
      right = [4.0_wp, 5.0_wp, 6.0_wp]

      mat = outer_matrix(left, right)

      ! First entry.
      call check(error, abs(mat(1, 1) - 4.0_wp) < 1.0e-12_wp, "M_11 should be 4")
      if (allocated(error)) return

      ! Mixed entry.
      call check(error, abs(mat(2, 3) - 12.0_wp) < 1.0e-12_wp, "M_23 should be 12")
      if (allocated(error)) return

      ! Last diagonal entry.
      call check(error, abs(mat(3, 3) - 18.0_wp) < 1.0e-12_wp, "M_33 should be 18")
   end subroutine test_outer_matrix

   !> Check logaddexp at large magnitudes.
   subroutine test_logaddexp_stability(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: x, result

      ! Large positive x.
      x = 100.0_wp
      result = logaddexp(x)
      call check(error, abs(result - x) < 1.0e-10_wp, "For large x, log(1+e^x) ≈ x")
      if (allocated(error)) return

      ! Large negative x.
      x = -100.0_wp
      result = logaddexp(x)
      call check(error, result >= 0.0_wp, "Result should be non-negative")
      if (allocated(error)) return
      call check(error, result < 1.0e-40_wp, "For large negative x, result should be near 0")
   end subroutine test_logaddexp_stability

   !> Check known logaddexp values.
   subroutine test_logaddexp_values(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: x, result, expected

      ! x = 0.
      x = 0.0_wp
      result = logaddexp(x)
      expected = log(2.0_wp)
      call check(error, abs(result - expected) < 1.0e-12_wp, "logaddexp(0) should be log(2)")
      if (allocated(error)) return

      ! x = 1.
      x = 1.0_wp
      result = logaddexp(x)
      expected = log(1.0_wp + exp(1.0_wp))
      call check(error, abs(result - expected) < 1.0e-12_wp, &
                 "logaddexp(1) should match direct computation")
   end subroutine test_logaddexp_values

   !> Compare eig_2x2_symmetric with LAPACK.
   subroutine test_eig_2x2_symmetric(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: R(2, 2), R_lapack(2, 2)
      real(wp) :: lambda_min_analytic, lambda_max_analytic
      real(wp) :: v_min_analytic(2), v_max_analytic(2)
      real(wp) :: evals_lapack(2), evecs_lapack(2, 2)
      real(wp) :: work_lapack(8)
      integer(lapack_ik) :: info
      real(wp) :: dot1, dot2

      ! Non-diagonal symmetric matrix.
      ! R = [4  1]
      !     [1  3]
      R = reshape([4.0_wp, 1.0_wp, 1.0_wp, 3.0_wp], [2, 2])

      ! Analytic solution.
      call eig_2x2_symmetric(R(1, 1), R(1, 2), R(2, 2), &
                             lambda_min_analytic, lambda_max_analytic, v_min_analytic, v_max_analytic)

      ! LAPACK solution.
      R_lapack = R
      call dsyev('V', 'U', int(2, lapack_ik), R_lapack, int(2, lapack_ik), evals_lapack, &
                 work_lapack, int(8, lapack_ik), info)
      if (info /= 0) then
         call check(error, .false., "LAPACK dsyev failed for test case 1")
         return
      end if
      evecs_lapack = R_lapack  ! dsyev returns eigenvectors in R_lapack.

      ! Compare eigenvalues.
      call check(error, abs(lambda_min_analytic - evals_lapack(1)) < 1.0e-12_wp, &
                 "Case 1: Smallest eigenvalue should match LAPACK")
      if (allocated(error)) return

      call check(error, abs(lambda_max_analytic - evals_lapack(2)) < 1.0e-12_wp, &
                 "Case 1: Largest eigenvalue should match LAPACK")
      if (allocated(error)) return

      ! Compare eigenvectors up to sign.
      dot1 = abs(dot_product(v_min_analytic, evecs_lapack(:, 1)))
      call check(error, abs(dot1 - 1.0_wp) < 1.0e-10_wp, &
                 "Case 1: Smallest eigenvector should match LAPACK (up to sign)")
      if (allocated(error)) return

      dot2 = abs(dot_product(v_max_analytic, evecs_lapack(:, 2)))
      call check(error, abs(dot2 - 1.0_wp) < 1.0e-10_wp, &
                 "Case 1: Largest eigenvector should match LAPACK (up to sign)")
      if (allocated(error)) return

      ! Diagonal matrix.
      ! R = [5  0]
      !     [0  2]
      R = reshape([5.0_wp, 0.0_wp, 0.0_wp, 2.0_wp], [2, 2])

      ! Analytic solution.
      call eig_2x2_symmetric(R(1, 1), R(1, 2), R(2, 2), &
                             lambda_min_analytic, lambda_max_analytic, v_min_analytic, v_max_analytic)

      ! LAPACK solution.
      R_lapack = R
      call dsyev('V', 'U', int(2, lapack_ik), R_lapack, int(2, lapack_ik), evals_lapack, &
                 work_lapack, int(8, lapack_ik), info)
      if (info /= 0) then
         call check(error, .false., "LAPACK dsyev failed for test case 2")
         return
      end if
      evecs_lapack = R_lapack

      ! Compare eigenvalues.
      call check(error, abs(lambda_min_analytic - evals_lapack(1)) < 1.0e-12_wp, &
                 "Case 2: Smallest eigenvalue should match LAPACK")
      if (allocated(error)) return

      call check(error, abs(lambda_max_analytic - evals_lapack(2)) < 1.0e-12_wp, &
                 "Case 2: Largest eigenvalue should match LAPACK")
      if (allocated(error)) return

      ! Compare eigenvectors up to sign.
      dot1 = abs(dot_product(v_min_analytic, evecs_lapack(:, 1)))
      call check(error, abs(dot1 - 1.0_wp) < 1.0e-10_wp, &
                 "Case 2: Smallest eigenvector should match LAPACK (up to sign)")
      if (allocated(error)) return

      dot2 = abs(dot_product(v_max_analytic, evecs_lapack(:, 2)))
      call check(error, abs(dot2 - 1.0_wp) < 1.0e-10_wp, &
                 "Case 2: Largest eigenvector should match LAPACK (up to sign)")
   end subroutine test_eig_2x2_symmetric

   ! outer4(v) = v_i v_j v_k v_l.

   !> Compare outer4 with an explicit reference.
   subroutine test_outer4_brute_force(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: v(3), t(3, 3, 3, 3), ref(3, 3, 3, 3)
      integer :: i, j, k, l

      v = [1.5_wp, -2.0_wp, 0.75_wp]

      t = outer4(v)
      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            ref(i, j, k, l) = v(i)*v(j)*v(k)*v(l)
         end do; end do; end do; end do

      call check(error, maxval(abs(t - ref)) < tensor_tol, &
                 "outer4: element-wise mismatch vs brute-force reference")
   end subroutine test_outer4_brute_force

   !> Check full outer4 symmetry.
   subroutine test_outer4_full_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: v(3), t(3, 3, 3, 3)
      integer :: i, j, k, l

      v = [0.7_wp, 1.3_wp, -0.4_wp]
      t = outer4(v)

      ! Adjacent swaps generate S_4.
      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            call check(error, abs(t(i, j, k, l) - t(j, i, k, l)) < tensor_tol, "outer4: i<->j")
            if (allocated(error)) return
            call check(error, abs(t(i, j, k, l) - t(i, k, j, l)) < tensor_tol, "outer4: j<->k")
            if (allocated(error)) return
            call check(error, abs(t(i, j, k, l) - t(i, j, l, k)) < tensor_tol, "outer4: k<->l")
            if (allocated(error)) return
         end do; end do; end do; end do
   end subroutine test_outer4_full_symmetry

   !> Cross-check outer4 with outer_matrix.
   subroutine test_outer4_outer_matrix(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: v(3), vvT(3, 3), t(3, 3, 3, 3), ref(3, 3, 3, 3)
      integer :: i, j, k, l

      v = [-1.1_wp, 2.2_wp, 0.3_wp]
      vvT = outer_matrix(v, v)
      t = outer4(v)

      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            ref(i, j, k, l) = vvT(i, j)*vvT(k, l)
         end do; end do; end do; end do

      call check(error, maxval(abs(t - ref)) < tensor_tol, &
                 "outer4(v) should factor as outer_matrix(v,v) on (i,j) and (k,l)")
   end subroutine test_outer4_outer_matrix

   !> Check outer4 for the zero vector.
   subroutine test_outer4_zero_vector(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: v(3), t(3, 3, 3, 3)

      v = 0.0_wp
      t = outer4(v)
      call check(error, maxval(abs(t)) < tensor_tol, "outer4(0) should be the zero tensor")
   end subroutine test_outer4_zero_vector

   ! sym4_31(g, h3): four-term symmetrisation.

   !> Compare sym4_31 with an explicit reference.
   subroutine test_sym4_31_brute_force(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: g(3), h3(3, 3, 3), t(3, 3, 3, 3), ref(3, 3, 3, 3)
      integer :: i, j, k, l, a, b, c

      g = [0.4_wp, -0.6_wp, 1.2_wp]
      ! Asymmetric h3 separates the four terms.
      do c = 1, 3; do b = 1, 3; do a = 1, 3
            h3(a, b, c) = real(a, wp) + 0.5_wp*real(b, wp) - 0.25_wp*real(c, wp)
         end do; end do; end do

      t = sym4_31(g, h3)
      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            ref(i, j, k, l) = g(i)*h3(j, k, l) + g(j)*h3(i, k, l) &
                              + g(k)*h3(i, j, l) + g(l)*h3(i, j, k)
         end do; end do; end do; end do

      call check(error, maxval(abs(t - ref)) < tensor_tol, &
                 "sym4_31: element-wise mismatch vs brute-force reference")
   end subroutine test_sym4_31_brute_force

   !> Symmetric h3 should give full output symmetry.
   subroutine test_sym4_31_symmetric_input(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: g(3), v(3), h3(3, 3, 3), t(3, 3, 3, 3)
      integer :: i, j, k, l

      g = [1.1_wp, -0.3_wp, 0.7_wp]
      v = [0.5_wp, 2.0_wp, -1.5_wp]
      h3 = outer3(v)
      t = sym4_31(g, h3)

      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            call check(error, abs(t(i, j, k, l) - t(j, i, k, l)) < tensor_tol, &
                       "sym4_31 with symmetric h3: must be symmetric in i<->j")
            if (allocated(error)) return
            call check(error, abs(t(i, j, k, l) - t(i, k, j, l)) < tensor_tol, &
                       "sym4_31 with symmetric h3: must be symmetric in j<->k")
            if (allocated(error)) return
            call check(error, abs(t(i, j, k, l) - t(i, j, l, k)) < tensor_tol, &
                       "sym4_31 with symmetric h3: must be symmetric in k<->l")
            if (allocated(error)) return
         end do; end do; end do; end do
   end subroutine test_sym4_31_symmetric_input

   !> Check sym4_31 scaling in each argument.
   subroutine test_sym4_31_linearity(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: g(3), h3(3, 3, 3), t(3, 3, 3, 3), t_scaled_g(3, 3, 3, 3), t_scaled_h(3, 3, 3, 3)
      integer :: a, b, c

      g = [0.2_wp, -1.0_wp, 0.5_wp]
      do c = 1, 3; do b = 1, 3; do a = 1, 3
            h3(a, b, c) = 0.3_wp*a - 0.1_wp*b + 0.7_wp*c
         end do; end do; end do

      t = sym4_31(g, h3)
      t_scaled_g = sym4_31(2.5_wp*g, h3)
      t_scaled_h = sym4_31(g, -3.0_wp*h3)

      call check(error, maxval(abs(t_scaled_g - 2.5_wp*t)) < tensor_tol, &
                 "sym4_31 must be linear in g")
      if (allocated(error)) return
      call check(error, maxval(abs(t_scaled_h - (-3.0_wp)*t)) < tensor_tol, &
                 "sym4_31 must be linear in h3")
   end subroutine test_sym4_31_linearity

   ! sym4_22(A, B): three 2+2 pair partitions.

   !> Compare sym4_22 with its 3-term reference.
   subroutine test_sym4_22_brute_force(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(3, 3), B(3, 3), t(3, 3, 3, 3), ref(3, 3, 3, 3)
      integer :: i, j, k, l

      call build_symmetric_pair(A, B)

      t = sym4_22(A, B)
      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            ref(i, j, k, l) = A(i, j)*B(k, l) + A(i, k)*B(j, l) + A(i, l)*B(j, k)
         end do; end do; end do; end do

      call check(error, maxval(abs(t - ref)) < tensor_tol, &
                 "sym4_22: element-wise mismatch vs documented 3-term reference")
   end subroutine test_sym4_22_brute_force

   !> Pin the partial, left-biased contract.
   subroutine test_sym4_22_left_biased(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(3, 3), B(3, 3), t_partial(3, 3, 3, 3), t_full(3, 3, 3, 3)

      call build_symmetric_pair(A, B)
      ! Avoid accidental equality with the full symmetrisation.
      B = B + reshape([0.0_wp, 0.4_wp, 0.0_wp, &
                       0.4_wp, 0.0_wp, 0.0_wp, &
                       0.0_wp, 0.0_wp, 0.0_wp], [3, 3])

      t_partial = sym4_22(A, B)
      t_full = sym4_22_full_ref(A, B)

      call check(error, maxval(abs(t_partial - t_full)) > 1.0e-3_wp, &
                 "sym4_22 must be left-biased (partial != fully symmetric) for A /= B")
   end subroutine test_sym4_22_left_biased

   !> Check j/k/l symmetry for symmetric inputs.
   subroutine test_sym4_22_jkl_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(3, 3), B(3, 3), t(3, 3, 3, 3)
      integer :: i, j, k, l

      call build_symmetric_pair(A, B)
      t = sym4_22(A, B)

      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            call check(error, abs(t(i, j, k, l) - t(i, k, j, l)) < tensor_tol, &
                       "sym4_22: must be symmetric in j<->k")
            if (allocated(error)) return
            call check(error, abs(t(i, j, k, l) - t(i, j, l, k)) < tensor_tol, &
                       "sym4_22: must be symmetric in k<->l")
            if (allocated(error)) return
            call check(error, abs(t(i, j, k, l) - t(i, l, k, j)) < tensor_tol, &
                       "sym4_22: must be symmetric in j<->l")
            if (allocated(error)) return
         end do; end do; end do; end do
   end subroutine test_sym4_22_jkl_symmetry

   !> Swapped calls should produce the full 6-term tensor.
   subroutine test_sym4_22_full_symmetric(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(3, 3), B(3, 3), t_sum(3, 3, 3, 3), t_full(3, 3, 3, 3)
      integer :: i, j, k, l

      call build_symmetric_pair(A, B)

      t_sum = sym4_22(A, B) + sym4_22(B, A)
      t_full = sym4_22_full_ref(A, B)

      call check(error, maxval(abs(t_sum - t_full)) < tensor_tol, &
                 "sym4_22(A,B) + sym4_22(B,A) must equal the full 6-term symmetric tensor")
      if (allocated(error)) return

      ! The sum must be fully symmetric.
      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            call check(error, abs(t_sum(i, j, k, l) - t_sum(j, i, k, l)) < tensor_tol, &
                       "sym4_22 + swapped: must be symmetric in i<->j")
            if (allocated(error)) return
         end do; end do; end do; end do
   end subroutine test_sym4_22_full_symmetric

   !> Check the equivalent four-vector contraction.
   subroutine test_sym4_22_contraction(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(3, 3), B(3, 3), t(3, 3, 3, 3)
      real(wp) :: u(3), v(3), w(3), x(3)
      real(wp) :: contracted, expected
      integer :: i, j, k, l

      call build_symmetric_pair(A, B)
      u = [1.0_wp, 0.5_wp, -0.3_wp]
      v = [0.2_wp, -1.1_wp, 0.4_wp]
      w = [-0.7_wp, 0.6_wp, 1.0_wp]
      x = [0.9_wp, 0.0_wp, -0.5_wp]

      t = sym4_22(A, B)
      contracted = 0.0_wp
      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            contracted = contracted + t(i, j, k, l)*u(i)*v(j)*w(k)*x(l)
         end do; end do; end do; end do

      expected = dot_product(u, matmul(A, v))*dot_product(w, matmul(B, x)) &
                 + dot_product(u, matmul(A, w))*dot_product(v, matmul(B, x)) &
                 + dot_product(u, matmul(A, x))*dot_product(v, matmul(B, w))

      call check(error, abs(contracted - expected) < tensor_tol, &
                 "sym4_22 contraction with four vectors must match bilinear-form sum")
   end subroutine test_sym4_22_contraction

   ! sym4_211(g, H): six-term symmetrisation.

   !> Compare sym4_211 with its 6-term reference.
   subroutine test_sym4_211_brute_force(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: g(3), H(3, 3), t(3, 3, 3, 3), ref(3, 3, 3, 3)
      integer :: i, j, k, l

      g = [0.6_wp, -1.4_wp, 0.9_wp]
      H = reshape([2.0_wp, 0.3_wp, -0.5_wp, &
                   0.3_wp, 1.1_wp, 0.8_wp, &
                   -0.5_wp, 0.8_wp, 3.0_wp], [3, 3])

      t = sym4_211(g, H)
      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            ref(i, j, k, l) = g(i)*g(j)*H(k, l) + g(i)*g(k)*H(j, l) + g(i)*g(l)*H(j, k) &
                              + g(j)*g(k)*H(i, l) + g(j)*g(l)*H(i, k) + g(k)*g(l)*H(i, j)
         end do; end do; end do; end do

      call check(error, maxval(abs(t - ref)) < tensor_tol, &
                 "sym4_211: element-wise mismatch vs brute-force reference")
   end subroutine test_sym4_211_brute_force

   !> Check full sym4_211 symmetry.
   subroutine test_sym4_211_full_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: g(3), H(3, 3), t(3, 3, 3, 3)
      integer :: i, j, k, l

      g = [0.6_wp, -1.4_wp, 0.9_wp]
      H = reshape([2.0_wp, 0.3_wp, -0.5_wp, &
                   0.3_wp, 1.1_wp, 0.8_wp, &
                   -0.5_wp, 0.8_wp, 3.0_wp], [3, 3])
      t = sym4_211(g, H)

      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            call check(error, abs(t(i, j, k, l) - t(j, i, k, l)) < tensor_tol, &
                       "sym4_211: must be symmetric in i<->j")
            if (allocated(error)) return
            call check(error, abs(t(i, j, k, l) - t(i, k, j, l)) < tensor_tol, &
                       "sym4_211: must be symmetric in j<->k")
            if (allocated(error)) return
            call check(error, abs(t(i, j, k, l) - t(i, j, l, k)) < tensor_tol, &
                       "sym4_211: must be symmetric in k<->l")
            if (allocated(error)) return
         end do; end do; end do; end do
   end subroutine test_sym4_211_full_symmetry

   !> Check sym4_211 scaling in each argument.
   subroutine test_sym4_211_scaling(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: g(3), H(3, 3), t(3, 3, 3, 3), t_g(3, 3, 3, 3), t_H(3, 3, 3, 3)

      g = [0.6_wp, -1.4_wp, 0.9_wp]
      H = reshape([2.0_wp, 0.3_wp, -0.5_wp, &
                   0.3_wp, 1.1_wp, 0.8_wp, &
                   -0.5_wp, 0.8_wp, 3.0_wp], [3, 3])

      t = sym4_211(g, H)
      t_g = sym4_211(2.0_wp*g, H)
      t_H = sym4_211(g, -1.5_wp*H)

      call check(error, maxval(abs(t_g - 4.0_wp*t)) < tensor_tol, &
                 "sym4_211 must scale as g^2 (quadratic in g)")
      if (allocated(error)) return
      call check(error, maxval(abs(t_H - (-1.5_wp)*t)) < tensor_tol, &
                 "sym4_211 must scale linearly in H")
   end subroutine test_sym4_211_scaling

   ! Tensor-test helpers.

   !> Build distinct symmetric matrices.
   pure subroutine build_symmetric_pair(A, B)
      real(wp), intent(out) :: A(3, 3), B(3, 3)

      A = reshape([1.0_wp, 0.5_wp, -0.2_wp, &
                   0.5_wp, 2.0_wp, 0.7_wp, &
                   -0.2_wp, 0.7_wp, 1.5_wp], [3, 3])

      B = reshape([3.0_wp, -0.4_wp, 0.6_wp, &
                   -0.4_wp, 1.2_wp, 0.1_wp, &
                   0.6_wp, 0.1_wp, 2.5_wp], [3, 3])
   end subroutine build_symmetric_pair

   !> Full 6-term 2+2 symmetric reference.
   pure function sym4_22_full_ref(A, B) result(t)
      real(wp), intent(in) :: A(3, 3), B(3, 3)
      real(wp) :: t(3, 3, 3, 3)
      integer :: i, j, k, l

      do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
            t(i, j, k, l) = A(i, j)*B(k, l) + A(i, k)*B(j, l) + A(i, l)*B(j, k) &
                            + B(i, j)*A(k, l) + B(i, k)*A(j, l) + B(i, l)*A(j, k)
         end do; end do; end do; end do
   end function sym4_22_full_ref

   !===========================================================================
   ! Raw-kernel ports: jacobwilliams/lusol (moist_math_linalg_lusol_ez%solve)
   !===========================================================================
   ! Solve sparse systems supplied in coordinate (COO) form and verify a
   ! near-zero residual ||A*x - b||. Reference problems are from lusol_test.f90.

   !> lusol_test test_1: a 3x3 system A*x = b in COO form.
   subroutine test_lusol_dense_3x3(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: m = 3, n = 3
      real(wp), parameter :: b(m) = real([1, 2, 3], wp)
      integer, parameter :: icol(m*n) = [1, 1, 1, 2, 2, 2, 3, 3, 3]
      integer, parameter :: irow(m*n) = [1, 2, 3, 1, 2, 3, 1, 2, 3]
      real(wp), parameter :: a(m*n) = real([1, 4, 7, 2, 5, 88, 3, 66, 9], wp)
      real(wp), parameter :: a_mat(m, n) = reshape(a, [m, n])
      real(wp) :: x(n)
      integer :: istat

      call solve(n, m, m*n, irow, icol, a, b, x, istat)
      call check(error, maxval(abs(matmul(a_mat, x) - b)) < lusol_thr, &
                 "LUSOL 3x3 residual too large")
   end subroutine test_lusol_dense_3x3

   !> lusol_test test_2: a rectangular 3x4 system (n > m) in COO form.
   subroutine test_lusol_rectangular_3x4(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: m = 3, n = 4
      real(wp), parameter :: b(m) = real([1, 2, 3], wp)
      integer, parameter :: icol(m*n) = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
      integer, parameter :: irow(m*n) = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]
      real(wp), parameter :: a(m*n) = [4.1_wp, 1.1_wp, 11.1_wp, &
                                       5.1_wp, -3.1_wp, 3.1_wp, &
                                       66.1_wp, 8.1_wp, -87.1_wp, &
                                       0.1_wp, -9.1_wp, 2.1_wp]
      real(wp), parameter :: a_mat(m, n) = reshape(a, [m, n])
      real(wp) :: x(n)
      integer :: istat

      call solve(n, m, m*n, irow, icol, a, b, x, istat)
      call check(error, maxval(abs(matmul(a_mat, x) - b)) < lusol_thr, &
                 "LUSOL 3x4 residual too large")
   end subroutine test_lusol_rectangular_3x4

   !===========================================================================
   ! Raw-kernel ports: jacobwilliams/LSQR (lsqr_solver_ez)
   !===========================================================================
   ! Same COO systems as the LUSOL ports, driven through the object-oriented
   ! lsqr_solver_ez. The Paige-Saunders generator path is covered by the LSMR
   ! ports below.

   !> lsqrtest_ez test_1: a 3x3 system A*x = b in COO form.
   subroutine test_lsqr_dense_3x3(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: m = 3, n = 3
      real(wp), parameter :: b(m) = real([1, 2, 3], wp)
      integer, parameter :: icol(m*n) = [1, 1, 1, 2, 2, 2, 3, 3, 3]
      integer, parameter :: irow(m*n) = [1, 2, 3, 1, 2, 3, 1, 2, 3]
      real(wp), parameter :: a(m*n) = real([1, 4, 7, 2, 5, 88, 3, 66, 9], wp)
      real(wp), parameter :: a_mat(m, n) = reshape(a, [m, n])
      type(lsqr_solver_ez) :: solver
      real(wp) :: x(n)
      integer :: istop

      call solver%initialize(m, n, a, irow, icol, itnlim=100)
      call solver%solve(b, 0.0_wp, x, istop)
      call check(error, maxval(abs(matmul(a_mat, x) - b)) < lsqr_thr, &
                 "LSQR 3x3 residual too large")
   end subroutine test_lsqr_dense_3x3

   !> lsqrtest_ez test_2: a rectangular 3x4 system (n > m) in COO form.
   subroutine test_lsqr_rectangular_3x4(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: m = 3, n = 4
      real(wp), parameter :: b(m) = real([1, 2, 3], wp)
      integer, parameter :: icol(m*n) = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
      integer, parameter :: irow(m*n) = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]
      real(wp), parameter :: a(m*n) = [4.1_wp, 1.1_wp, 11.1_wp, &
                                       5.1_wp, -3.1_wp, 3.1_wp, &
                                       66.1_wp, 8.1_wp, -87.1_wp, &
                                       0.1_wp, -9.1_wp, 2.1_wp]
      real(wp), parameter :: a_mat(m, n) = reshape(a, [m, n])
      type(lsqr_solver_ez) :: solver
      real(wp) :: x(n)
      integer :: istop

      call solver%initialize(m, n, a, irow, icol, itnlim=100)
      call solver%solve(b, 0.0_wp, x, istop)
      call check(error, maxval(abs(matmul(a_mat, x) - b)) < lsqr_thr, &
                 "LSQR 3x4 residual too large")
   end subroutine test_lsqr_rectangular_3x4

   !===========================================================================
   ! Raw-kernel ports: jacobwilliams/LSMR (matrix-free Paige-Saunders problem)
   !===========================================================================
   ! A = Y*D*Z applied matrix-free via the Aprod callbacks; lstp builds b from a
   ! known xtrue and LSMR must recover it. Three groups sweep damping values.

   !> m = 2*nbar, n = nbar.
   subroutine test_lsmr_over_determined(error)
      type(error_type), allocatable, intent(out) :: error

      call lsmr_run_group(2*lsmr_nbar, lsmr_nbar, "over", error)
   end subroutine test_lsmr_over_determined

   !> m = n = nbar.
   subroutine test_lsmr_square(error)
      type(error_type), allocatable, intent(out) :: error

      call lsmr_run_group(lsmr_nbar, lsmr_nbar, "square", error)
   end subroutine test_lsmr_square

   !> m = nbar, n = 2*nbar.
   subroutine test_lsmr_under_determined(error)
      type(error_type), allocatable, intent(out) :: error

      call lsmr_run_group(lsmr_nbar, 2*lsmr_nbar, "under", error)
   end subroutine test_lsmr_under_determined

   !> Generate and solve the test problem for a sweep of damping values,
   !> asserting that LSMR recovers the true solution. The factored operator and
   !> scratch vectors are local automatic arrays reached by the nested Aprod
   !> callbacks via host association: LSMR drives the solve through those
   !> callbacks, whose fixed argument list cannot carry the operator.
   subroutine lsmr_run_group(m, n, label, error)
      integer(ip), intent(in) :: m
      integer(ip), intent(in) :: n
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error

      !> Factored operator A = Y*D*Z: Householder vectors Y, Z and singular values D.
      real(wp) :: lsmr_d(min(m, n)), lsmr_hy(m), lsmr_hz(n)
      !> Reverse-communication scratch vectors (lengths m and n).
      real(wp) :: lsmr_wm(m), lsmr_wn(n)
      real(wp) :: b(m), x(n), xtrue(n)
      real(wp) :: damp, atol, btol, conlim, normA, condA, normr, normAr, normx
      real(wp) :: norme
      integer(ip) :: istop, itn, itnlim, j, ndamp, npower

      do ndamp = 3, 8
         npower = ndamp
         damp = 10.0_wp**(-ndamp)
         if (ndamp == 8) damp = 0.0_wp

         do j = 1, n
            xtrue(j) = 0.1_wp*j
         end do

         call lsmr_lstp(m, n, lsmr_nduplc, npower, damp, xtrue, b, condA, normr)

         atol = 1.0e-12_wp
         btol = atol
         conlim = 1000.0_wp*condA
         itnlim = 4*(m + n + 50)

         call lsmr(m, n, lsmr_aprod1, lsmr_aprod2, b, damp, atol, btol, conlim, itnlim, &
                   0_ip, 0_ip, x, istop, itn, normA, condA, normr, normAr, normx)

         call check(error, istop >= 0 .and. istop <= 7, &
                    label//": LSMR error istop, ndamp="//to_string(ndamp))
         if (allocated(error)) return

         norme = sqrt(dot_product(x - xtrue, x - xtrue))/(1.0_wp + sqrt(dot_product(xtrue, xtrue)))
         call check(error, norme < lsmr_etol, &
                    label//": solution error too large, ndamp="//to_string(ndamp))
         if (allocated(error)) return
      end do

   contains

      !> y := y + A*x with A = Y*D*Z.
      subroutine lsmr_aprod1(m, n, x, y)
         integer(ip), intent(in) :: m
         integer(ip), intent(in) :: n
         real(wp), intent(in) :: x(n)
         real(wp), intent(inout) :: y(m)

         integer(ip) :: minmn

         minmn = min(m, n)
         lsmr_wn = x
         call lsmr_hprod(n, lsmr_hz, lsmr_wn)
         lsmr_wm(1:minmn) = lsmr_d(1:minmn)*lsmr_wn(1:minmn)
         lsmr_wm(n + 1:m) = 0.0_wp
         call lsmr_hprod(m, lsmr_hy, lsmr_wm)
         y = y + lsmr_wm
      end subroutine lsmr_aprod1

      !> x := x + A'*y with A = Y*D*Z.
      subroutine lsmr_aprod2(m, n, x, y)
         integer(ip), intent(in) :: m
         integer(ip), intent(in) :: n
         real(wp), intent(inout) :: x(n)
         real(wp), intent(in) :: y(m)

         integer(ip) :: minmn

         minmn = min(m, n)
         lsmr_wm = y
         call lsmr_hprod(m, lsmr_hy, lsmr_wm)
         lsmr_wn(1:minmn) = lsmr_d(1:minmn)*lsmr_wm(1:minmn)
         lsmr_wn(m + 1:n) = 0.0_wp
         call lsmr_hprod(n, lsmr_hz, lsmr_wn)
         x = x + lsmr_wn
      end subroutine lsmr_aprod2

      !> Generate the least-squares test problem (A in factored form, the
      !> right-hand side b, the condition number and the residual norm).
      subroutine lsmr_lstp(m, n, nduplc, npower, damp, x, b, condA, normr)
         integer(ip), intent(in) :: m
         integer(ip), intent(in) :: n
         integer(ip), intent(in) :: nduplc
         integer(ip), intent(in) :: npower
         real(wp), intent(in) :: damp
         real(wp), intent(inout) :: x(n)
         real(wp), intent(out) :: b(m)
         real(wp), intent(out) :: condA
         real(wp), intent(out) :: normr

         integer(ip) :: i, j, minmn
         real(wp) :: alfa, beta, dampsq, fourpi, t

         minmn = min(m, n)
         dampsq = damp**2
         fourpi = 4.0_wp*3.141592_wp
         alfa = fourpi/m
         beta = fourpi/n

         do i = 1, m
            lsmr_hy(i) = sin(alfa*i)
         end do
         do i = 1, n
            lsmr_hz(i) = cos(beta*i)
         end do

         alfa = sqrt(dot_product(lsmr_hy, lsmr_hy))
         beta = sqrt(dot_product(lsmr_hz, lsmr_hz))
         lsmr_hy = (-1.0_wp/alfa)*lsmr_hy
         lsmr_hz = (-1.0_wp/beta)*lsmr_hz

         ! Singular values of A.
         do i = 1, minmn
            j = (i - 1 + nduplc)/nduplc
            t = real(j*nduplc, wp)/real(minmn, wp)
            lsmr_d(i) = t**npower
         end do

         condA = sqrt((lsmr_d(minmn)**2 + dampsq)/(lsmr_d(1)**2 + dampsq))

         lsmr_wn = x
         call lsmr_hprod(n, lsmr_hz, lsmr_wn)
         if (m < n) then
            lsmr_wn(m + 1:n) = 0.0_wp
            call lsmr_hprod(n, lsmr_hz, lsmr_wn)
            x = lsmr_wn
         end if

         lsmr_wm(1:minmn) = dampsq*lsmr_wn(1:minmn)/lsmr_d(1:minmn)
         lsmr_wm(minmn + 1:m) = 1.0_wp
         call lsmr_hprod(m, lsmr_hy, lsmr_wm)

         normr = sqrt(dot_product(lsmr_wm, lsmr_wm))
         b = lsmr_wm
         call lsmr_aprod1(m, n, x, b)
      end subroutine lsmr_lstp

   end subroutine lsmr_run_group

   !> Apply a Householder transformation: x := (I - 2*z*z')*x.
   subroutine lsmr_hprod(n, z, x)
      integer(ip), intent(in) :: n
      real(wp), intent(in) :: z(n)
      real(wp), intent(inout) :: x(n)

      real(wp) :: s

      s = 2.0_wp*dot_product(z, x)
      x = x - s*z
   end subroutine lsmr_hprod

end module test_math_linalg
