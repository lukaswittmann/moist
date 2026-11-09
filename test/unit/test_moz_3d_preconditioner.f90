!> Cartesian auxiliary-space preconditioner of the 3D UV Jacobian: moz_3d_preconditioner
!>
!> - Bounded flexible GMRES: a dense system solved to roundoff, also under a preconditioner
!>   that changes every step; a bounded solve accepted with its step and product counts;
!>   the identity accepted at once; breakdown in an invariant subspace; an inconsistent
!>   and a non-finite operator rejected with the identity action
!> - Lattice Jacobian against a dense J = I - F^-1 (P^T - I) F D assembled from explicit
!>   DFT matrices and solved by LAPACK, two coupled sites with a nonsymmetric response, on
!>   the direct and the nested path; the molecular action against v + P(J^-1 R v - R v),
!>   with the lattice operator and the transfers built by the test from the same inputs
!> - A zero closure derivative and vectors R does not see leave v unchanged on both paths
!> - A failed lattice solve returns v and counts a fallback
!> - The lattice is capped at the grid's Fourier modes; invalid sizes and inner settings,
!>   grids without Fourier modes or never updated and negative weights are refused, each in a
!>   test of its own
module test_moz_3d_preconditioner
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_lapack, only: getrf, getrs, lapack_gels
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type
   use moist_math_grid_deposit, only: grid_deposit_type
   use moist_model_moz_3d_preconditioner, only: uv3d_aux_preconditioner_type, uv3d_aux_statistics_type, &
      & uv3d_aux_response_source_type, uv3d_aux_level_type, new_uv3d_aux_level
   use moist_model_moz_3d_preconditioner_krylov, only: bounded_operator_type, bounded_gmres_type, &
      & bounded_gmres_result_type
   use moist_model_moz_solver_jfnk_kernels, only: upper_backsolve
   use test_helpers, only: get_qc_handymod_recipe, expect_moist_error, get_test_water
   implicit none(type, external)
   private

   public :: collect_moz_3d_preconditioner

   !> Dense test operator with optional faults and a preconditioner that changes every step
   type, extends(bounded_operator_type) :: dense_operator_type
      !> Matrix
      real(wp), allocatable :: a(:, :)
      !> Products so far
      integer :: nproduct = 0
      !> Products after this one are wrong
      integer :: corrupt_after = huge(0)
      !> Products contain a NaN
      logical :: nonfinite = .false.
      !> Preconditioner alternates between two diagonal scalings
      logical :: varying = .false.
      !> Preconditioner applications so far
      integer :: nprec = 0
   contains
      procedure :: matvec => dense_matvec
      procedure :: precondition => dense_precondition
   end type dense_operator_type

   !> Two-site response of `coupling` at the k-points of a lattice
   type, extends(uv3d_aux_response_source_type) :: synthetic_source_type
   contains
      procedure :: tabulate => synthetic_tabulate
   end type synthetic_source_type

contains

   !> Collect the preconditioner tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_3d_preconditioner(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("krylov_dense", test_krylov_dense), &
         & new_unittest("krylov_flexible", test_krylov_flexible), &
         & new_unittest("krylov_bounded", test_krylov_bounded), &
         & new_unittest("krylov_identity", test_krylov_identity), &
         & new_unittest("krylov_breakdown", test_krylov_breakdown), &
         & new_unittest("krylov_rejected", test_krylov_rejected), &
         & new_unittest("dense_direct", test_dense_direct), &
         & new_unittest("dense_nested", test_dense_nested), &
         & new_unittest("identity_derivative", test_identity_derivative), &
         & new_unittest("unseen_modes", test_unseen_modes), &
         & new_unittest("fallback", test_fallback), &
         & new_unittest("size_cap", test_size_cap), &
         & new_unittest("refused_odd_cells", test_refused_odd_cells, should_fail=.true.), &
         & new_unittest("refused_coarse_size", test_refused_coarse_size, should_fail=.true.), &
         & new_unittest("refused_inner_steps", test_refused_inner_steps, should_fail=.true.), &
         & new_unittest("refused_inner_rtol", test_refused_inner_rtol, should_fail=.true.), &
         & new_unittest("refused_negative_weight", test_refused_negative_weight, should_fail=.true.), &
         & new_unittest("refused_no_modes", test_refused_no_modes, should_fail=.true.), &
         & new_unittest("refused_unbuilt_grid", test_refused_unbuilt_grid, should_fail=.true.) &
         & ]
   end subroutine collect_moz_3d_preconditioner

   !> y = A x, wrong after `corrupt_after` products, with a NaN when `nonfinite`
   subroutine dense_matvec(self, x, y, ok)
      class(dense_operator_type), intent(inout) :: self
      real(wp), contiguous, intent(in) :: x(:)
      real(wp), contiguous, intent(out) :: y(:)
      logical, intent(out) :: ok
      self%nproduct = self%nproduct + 1
      y = matmul(self%a, x)
      if (self%nproduct > self%corrupt_after) y = -3.0_wp*y
      if (self%nonfinite) y(1) = ieee_value(1.0_wp, ieee_quiet_nan)
      ok = .true.
   end subroutine dense_matvec

   !> z = q, or alternately two different diagonal scalings of q
   subroutine dense_precondition(self, q, z, ok)
      class(dense_operator_type), intent(inout) :: self
      real(wp), contiguous, intent(in) :: q(:)
      real(wp), contiguous, intent(out) :: z(:)
      logical, intent(out) :: ok
      integer :: i
      ok = .true.
      if (.not. self%varying) then
         z = q
         return
      end if
      self%nprec = self%nprec + 1
      do i = 1, size(q)
         if (mod(self%nprec, 2) == 1) then
            z(i) = q(i)/(1.0_wp + 0.3_wp*real(i, wp))
         else
            z(i) = q(i)*(1.0_wp + 0.1_wp*real(mod(i, 3), wp))
         end if
      end do
   end subroutine dense_precondition

   !> Nonsymmetric, well-conditioned 6x6 matrix
   !>
   !> @param[out] a  matrix (6, 6)
   subroutine dense_matrix(a)
      !> Matrix (6, 6)
      real(wp), allocatable, intent(out) :: a(:, :)
      integer :: i, j
      allocate (a(6, 6))
      do j = 1, 6
         do i = 1, 6
            a(i, j) = 0.3_wp*sin(real(i + 2*j, wp))
         end do
         a(j, j) = a(j, j) + 2.0_wp + 0.5_wp*real(j, wp)
      end do
   end subroutine dense_matrix

   !> Solution of a dense system by LAPACK
   !>
   !> @param[in]  a      matrix (n, n)
   !> @param[in]  b      right-hand side (n)
   !> @param[out] x      solution (n)
   !> @param[out] error  reference-solve failure
   subroutine dense_solve(a, b, x, error)
      !> Reference-solve failure
      type(error_type), allocatable, intent(out) :: error
      !> Matrix (n, n) and right-hand side (n)
      real(wp), intent(in) :: a(:, :), b(:)
      !> Solution (n)
      real(wp), allocatable, intent(out) :: x(:)
      real(wp), allocatable :: lu(:, :), rhs(:, :)
      integer, allocatable :: ipiv(:)
      integer :: info
      lu = a
      allocate (ipiv(size(b)))
      call getrf(lu, ipiv, info)
      call check(error, info == 0, "LAPACK dense reference solve failed")
      if (allocated(error)) return
      rhs = reshape(b, [size(b), 1])
      call getrs(lu, rhs, ipiv, info)
      call check(error, info == 0, "LAPACK dense reference solve failed")
      if (allocated(error)) return
      x = rhs(:, 1)
   end subroutine dense_solve

   !> Six steps solve a 6x6 system to roundoff
   subroutine test_krylov_dense(error)
      type(error_type), allocatable, intent(out) :: error
      type(dense_operator_type) :: op
      type(bounded_gmres_type) :: gmres
      type(bounded_gmres_result_type) :: result
      real(wp), allocatable :: b(:), x(:), ref(:)

      call dense_matrix(op%a)
      b = [1.0_wp, -2.0_wp, 0.5_wp, 3.0_wp, -1.0_wp, 0.25_wp]
      allocate (x(6))
      call gmres%setup(6, 6, 1)
      call gmres%solve(op, b, x, 1.0e-14_wp, result)
      call dense_solve(op%a, b, ref, error)
      if (allocated(error)) return
      call check(error, result%accepted, "dense solve not accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(x - ref)) <= 1.0e-12_wp*maxval(abs(ref)), "dense solve off the LAPACK solution")
      if (allocated(error)) return
      call check(error, result%nmatvec == result%iterations + 2, "products: start, steps and the check")
      if (allocated(error)) return

      ! Four of six steps: the least-squares solve reads the 4x4 block, below the band too
      call gmres%setup(6, 4, 1)
      call gmres%solve(op, b, x, 1.0e-14_wp, result)
      call check(error, result%accepted .and. result%iterations == 4 .and. result%nmatvec == 6 .and. &
         & norm2(b - matmul(op%a, x)) < norm2(b - matmul(op%a, b)), "four bounded steps not accepted")
      if (allocated(error)) return

      ! Independent least-squares minimiser over span{r0, A r0, A**2 r0, A**3 r0}.
      block
         real(wp) :: q(6, 4), aq(6, 4), rhs(6, 1), work(64), r0(6), expected(6)
         integer(lapack_ik) :: info
         integer :: k

         r0 = b - matmul(op%a, b)
         q(:, 1) = r0/norm2(r0)
         do k = 2, 4
            q(:, k) = matmul(op%a, q(:, k - 1))
            q(:, k) = q(:, k)/norm2(q(:, k))
         end do
         aq = matmul(op%a, q)
         rhs(:, 1) = r0
         call lapack_gels("N", 6_lapack_ik, 4_lapack_ik, 1_lapack_ik, aq, 6_lapack_ik, &
            & rhs, 6_lapack_ik, work, 64_lapack_ik, info)
         call check(error, info == 0, "LAPACK Krylov least-squares reference failed")
         if (allocated(error)) return
         expected = b + matmul(q, rhs(1:4, 1))
         call check(error, norm2(b - matmul(op%a, x)), norm2(b - matmul(op%a, expected)), &
            & thr=1.0e-12_wp, more="four-step residual differs from the least-squares minimum")
         if (allocated(error)) return
         call check(error, norm2(x - expected) <= 1.0e-10_wp*norm2(expected), &
            & "four-step solution differs from the independent Krylov minimiser")
         if (allocated(error)) return
      end block

      ! The back-solve reads the upper triangle only
      block
         real(wp) :: r(3, 3), y(3), nan
         logical :: feasible
         nan = ieee_value(1.0_wp, ieee_quiet_nan)
         r = reshape([2.0_wp, nan, nan, 1.0_wp, 4.0_wp, nan, 0.5_wp, 1.0_wp, 5.0_wp], [3, 3])
         call upper_backsolve(r, [3.5_wp, 5.0_wp, 5.0_wp], y, feasible)
         call check(error, feasible .and. maxval(abs(y - [1.0_wp, 1.0_wp, 1.0_wp])) < 1.0e-15_wp, &
            & "entries below the diagonal spoiled the back-solve")
      end block
   end subroutine test_krylov_dense

   !> A preconditioner that changes every step does not spoil the flexible solve
   subroutine test_krylov_flexible(error)
      type(error_type), allocatable, intent(out) :: error
      type(dense_operator_type) :: op
      type(bounded_gmres_type) :: gmres
      type(bounded_gmres_result_type) :: result
      real(wp), allocatable :: b(:), x(:), ref(:)

      call dense_matrix(op%a)
      op%varying = .true.
      b = [0.5_wp, 1.0_wp, -1.5_wp, 2.0_wp, 0.0_wp, -0.75_wp]
      allocate (x(6))
      call gmres%setup(6, 6, 1)
      call gmres%solve(op, b, x, 1.0e-14_wp, result)
      call dense_solve(op%a, b, ref, error)
      if (allocated(error)) return
      call check(error, result%accepted .and. op%nprec >= 2, "flexible solve not accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(x - ref)) <= 1.0e-12_wp*maxval(abs(ref)), "flexible solve off the LAPACK solution")
   end subroutine test_krylov_flexible

   !> Two steps reduce the residual and are accepted with four products
   subroutine test_krylov_bounded(error)
      type(error_type), allocatable, intent(out) :: error
      type(dense_operator_type) :: op
      type(bounded_gmres_type) :: gmres
      type(bounded_gmres_result_type) :: result
      real(wp), allocatable :: b(:), x(:)

      call dense_matrix(op%a)
      b = [1.0_wp, 1.0_wp, -1.0_wp, 2.0_wp, 0.5_wp, -0.5_wp]
      allocate (x(6))
      call gmres%setup(6, 2, 1)
      call gmres%solve(op, b, x, 1.0e-3_wp, result)
      call check(error, result%accepted .and. result%iterations == 2 .and. result%nmatvec == 4, &
         & "bounded solve: accepted after two steps and four products")
      if (allocated(error)) return
      call check(error, norm2(b - matmul(op%a, x)) < norm2(b - matmul(op%a, b)), "residual not reduced")
   end subroutine test_krylov_bounded

   !> The identity operator: accepted at once with x = b
   subroutine test_krylov_identity(error)
      type(error_type), allocatable, intent(out) :: error
      type(dense_operator_type) :: op
      type(bounded_gmres_type) :: gmres
      type(bounded_gmres_result_type) :: result
      real(wp), allocatable :: b(:), x(:)
      integer :: i

      allocate (op%a(5, 5), source=0.0_wp)
      do i = 1, 5
         op%a(i, i) = 1.0_wp
      end do
      b = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
      allocate (x(5))
      call gmres%setup(5, 2, 1)
      call gmres%solve(op, b, x, 0.1_wp, result)
      call check(error, result%accepted .and. result%iterations == 0 .and. result%nmatvec == 1 .and. all(x == b), &
         & "identity: x = b after one product")
   end subroutine test_krylov_identity

   !> A residual in a two-dimensional invariant subspace: breakdown after two steps, exact
   subroutine test_krylov_breakdown(error)
      type(error_type), allocatable, intent(out) :: error
      type(dense_operator_type) :: op
      type(bounded_gmres_type) :: gmres
      type(bounded_gmres_result_type) :: result
      real(wp), allocatable :: b(:), x(:)
      integer :: i

      allocate (op%a(6, 6), source=0.0_wp)
      do i = 1, 6
         op%a(i, i) = real(i, wp)
      end do
      b = [1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      allocate (x(6))
      call gmres%setup(6, 6, 1)
      call gmres%solve(op, b, x, 0.0_wp, result)
      call check(error, result%accepted .and. result%iterations == 2 .and. result%nmatvec == 4, "breakdown after two steps")
      if (allocated(error)) return
      call check(error, maxval(abs(x - [1.0_wp, 0.5_wp, 1.0_wp/3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp])) <= 1.0e-13_wp, &
         & "breakdown solution not exact")
   end subroutine test_krylov_breakdown

   !> A product that changes before the check, or a NaN, gives the identity action
   subroutine test_krylov_rejected(error)
      type(error_type), allocatable, intent(out) :: error
      type(dense_operator_type) :: op
      type(bounded_gmres_type) :: gmres
      type(bounded_gmres_result_type) :: result
      real(wp), allocatable :: b(:), x(:)

      call dense_matrix(op%a)
      op%corrupt_after = 3
      b = [1.0_wp, -1.0_wp, 2.0_wp, 0.0_wp, 1.0_wp, 0.5_wp]
      allocate (x(6))
      call gmres%setup(6, 2, 1)
      call gmres%solve(op, b, x, 0.1_wp, result)
      call check(error, .not. result%accepted .and. all(x == b) .and. result%nmatvec == 4 .and. result%iterations == 2, &
         & "a worse recomputed residual was accepted, or its steps not reported")
      if (allocated(error)) return
      op%corrupt_after = huge(0)
      op%nonfinite = .true.
      call gmres%solve(op, b, x, 0.1_wp, result)
      call check(error, .not. result%accepted .and. all(x == b), "a non-finite product was accepted")
   end subroutine test_krylov_rejected

   !> Small molecular grid of water with its Fourier modes
   !>
   !> @param[out] grid        molecular grid
   !> @param[out] err         error handling
   !> @param[in]  reciprocal  configure the reciprocal geometry; passed to the grid constructor
   subroutine water_grid(grid, err, reciprocal)
      !> Molecular grid
      type(moist_math_grid_3d_molecular_type), intent(out) :: grid
      !> Error handling
      type(moist_error), allocatable, intent(out) :: err
      !> Configure the reciprocal geometry
      logical, intent(in), optional :: reciprocal
      type(structure_type) :: mol
      type(moist_math_grid_atomic_recipe_type) :: recipe

      call get_test_water(mol)
      call get_qc_handymod_recipe(recipe, err, nrad=10, degree=7, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(grid, err, recipe=recipe, pruning_threshold=0.1_wp, &
         & dr=0.8_wp, reciprocal=reciprocal)
      if (.not. allocated(err)) call grid%update(mol, err)
   end subroutine water_grid

   !> Coupling P(|k|) of two sites, nonsymmetric
   !>
   !> @param[in] k2  squared wave number, 1/bohr^2
   pure function coupling(k2) result(p)
      !> Squared wave number, 1/bohr^2
      real(wp), intent(in) :: k2
      real(wp) :: p(2, 2)
      p = reshape([0.4_wp, -0.2_wp, 0.3_wp, 0.6_wp], [2, 2])*exp(-0.5_wp*k2)
      p(1, 1) = p(1, 1) + 1.0_wp
      p(2, 2) = p(2, 2) + 1.0_wp
   end function coupling

   !> Two-site response table of `coupling` at the k-points of a grid
   !>
   !> @param[in]  grid      grid with its k-points
   !> @param[out] response  response table
   subroutine synthetic_response(grid, response)
      !> Grid with its k-points
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Response table
      type(uv3d_susceptibility_type), intent(out) :: response
      real(wp) :: k(3)
      integer :: i

      response%nv = 2
      allocate (response%values(2, 2, grid%npts_k), response%index(grid%npts_k), response%radius(grid%npts_k))
      do i = 1, grid%npts_k
         k = grid%kpoint(i)
         response%index(i) = i
         response%radius(i) = norm2(k)
         response%values(:, :, i) = coupling(sum(k**2))
      end do
   end subroutine synthetic_response

   !> Response of `coupling` for the lattices of the preconditioner
   subroutine synthetic_tabulate(self, grid, response, error)
      class(synthetic_source_type), intent(in) :: self
      class(moist_math_grid_3d_type), intent(in) :: grid
      type(uv3d_susceptibility_type), intent(out) :: response
      type(moist_error), allocatable, intent(out) :: error
      call synthetic_response(grid, response)
   end subroutine synthetic_tabulate

   !> Closure derivative of a KH-like range on the grid points
   !>
   !> @param[in]  grid  grid
   !> @param[out] dcdt  closure derivative (ngrid, 2)
   subroutine derivative(grid, dcdt)
      !> Grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Closure derivative (ngrid, 2)
      real(wp), allocatable, intent(out) :: dcdt(:, :)
      integer :: i
      allocate (dcdt(grid%ngrid, 2))
      do i = 1, grid%ngrid
         dcdt(i, 1) = -0.5_wp - 0.4_wp*sin(0.7_wp*grid%xyz(1, i) + 0.3_wp*grid%xyz(3, i))
         dcdt(i, 2) = -0.3_wp + 0.25_wp*cos(0.5_wp*grid%xyz(2, i))
      end do
   end subroutine derivative

   !> A molecular vector of two columns
   !>
   !> @param[in]  grid  grid
   !> @param[out] v     flat vector (2 ngrid)
   subroutine molecular_vector(grid, v)
      !> Grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Flat vector (2 ngrid)
      real(wp), allocatable, intent(out) :: v(:)
      integer :: i
      allocate (v(2*grid%ngrid))
      do i = 1, grid%ngrid
         v(i) = sin(real(i, wp)) + 0.3_wp*grid%xyz(1, i)
         v(grid%ngrid + i) = cos(1.3_wp*real(i, wp)) - 0.2_wp*grid%xyz(2, i)
      end do
   end subroutine molecular_vector

   !> Preconditioner on the water grid with synthetic responses, refreshed
   !>
   !> @param[out] grid             water grid
   !> @param[out] aux              refreshed preconditioner
   !> @param[in]  n                fine cells per axis
   !> @param[in]  coarse_n         coarse cells per axis
   !> @param[in]  max_iter         most inner steps
   !> @param[in]  rtol             inner residual reduction
   !> @param[out] err              error handling
   !> @param[in]  zero_derivative  refresh with a zero closure derivative; default false
   subroutine aux_fixture(grid, aux, n, coarse_n, max_iter, rtol, err, zero_derivative)
      !> Water grid
      type(moist_math_grid_3d_molecular_type), intent(out) :: grid
      !> Refreshed preconditioner
      type(uv3d_aux_preconditioner_type), intent(out) :: aux
      !> Fine and coarse cells per axis, most inner steps
      integer, intent(in) :: n, coarse_n, max_iter
      !> Inner residual reduction
      real(wp), intent(in) :: rtol
      !> Error handling
      type(moist_error), allocatable, intent(out) :: err
      !> Refresh with a zero closure derivative
      logical, intent(in), optional :: zero_derivative
      type(synthetic_source_type) :: source
      real(wp), allocatable :: dcdt(:, :)

      call water_grid(grid, err)
      if (allocated(err)) return
      call aux%setup(grid, source, 1.7_wp, n, coarse_n, max_iter, rtol, 1, err)
      if (allocated(err)) return
      call derivative(grid, dcdt)
      if (present(zero_derivative)) then
         if (zero_derivative) dcdt = 0.0_wp
      end if
      call aux%attach()
      call aux%refresh(dcdt)
   end subroutine aux_fixture

   !> Dense lattice Jacobian J = I - K D with K from explicit DFT matrices:
   !> K_ga(r, r') = (1/N) sum_k exp(i k.(r - r')) (P_ga(|k|) - delta_ga)
   !>
   !> @param[in]  n     cells per axis
   !> @param[in]  side  lattice side, bohr
   !> @param[in]  d     closure derivative on the cells (n^3, 2)
   !> @param[out] jac   dense Jacobian (2 n^3, 2 n^3)
   subroutine dense_lattice_jacobian(n, side, d, jac)
      !> Cells per axis
      integer, intent(in) :: n
      !> Lattice side, bohr
      real(wp), intent(in) :: side
      !> Closure derivative on the cells (n^3, 2)
      real(wp), intent(in) :: d(:, :)
      !> Dense Jacobian (2 n^3, 2 n^3)
      real(wp), allocatable, intent(out) :: jac(:, :)
      real(wp), parameter :: two_pi = 8.0_wp*atan(1.0_wp)
      real(wp), allocatable :: kmat(:, :, :, :)
      integer :: ncell, a, b, mx, my, mz, g, s, ia(3), ib(3)
      real(wp) :: k(3), p(2, 2), phase
      integer :: m(3)

      ncell = n**3
      allocate (kmat(ncell, ncell, 2, 2), source=0.0_wp)
      do mz = 0, n - 1
         do my = 0, n - 1
            do mx = 0, n - 1
               m = [mx, my, mz]
               where (m > n/2) m = m - n
               k = two_pi*real(m, wp)/side
               p = coupling(sum(k**2))
               p(1, 1) = p(1, 1) - 1.0_wp
               p(2, 2) = p(2, 2) - 1.0_wp
               do b = 1, ncell
                  ib = [mod(b - 1, n), mod((b - 1)/n, n), (b - 1)/n**2]
                  do a = 1, ncell
                     ia = [mod(a - 1, n), mod((a - 1)/n, n), (a - 1)/n**2]
                     phase = two_pi*real(dot_product([mx, my, mz], ia - ib), wp)/real(n, wp)
                     do s = 1, 2
                        do g = 1, 2
                           kmat(a, b, s, g) = kmat(a, b, s, g) + cos(phase)*p(g, s)/real(ncell, wp)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      ! Row (a, s), column (b, g): delta - K_gs(a, b) D_g(b)
      allocate (jac(2*ncell, 2*ncell), source=0.0_wp)
      do g = 1, 2
         do b = 1, ncell
            do s = 1, 2
               do a = 1, ncell
                  jac(a + ncell*(s - 1), b + ncell*(g - 1)) = -kmat(a, b, s, g)*d(b, g)
               end do
            end do
            jac(b + ncell*(g - 1), b + ncell*(g - 1)) = jac(b + ncell*(g - 1), b + ncell*(g - 1)) + 1.0_wp
         end do
      end do
   end subroutine dense_lattice_jacobian

   !> Lattice solve and molecular action against the dense inverse
   !>
   !> @param[in]  coarse_n  coarse cells of the 4-cell lattice: 2 nests, 4 does not
   !> @param[out] error     test failure
   subroutine dense_case(coarse_n, error)
      !> Coarse cells of the 4-cell lattice
      integer, intent(in) :: coarse_n
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: grid
      type(uv3d_aux_preconditioner_type) :: aux
      type(uv3d_aux_statistics_type) :: stats
      type(grid_deposit_type) :: deposit
      type(uv3d_aux_level_type) :: level
      type(synthetic_source_type) :: source
      type(bounded_gmres_type) :: gmres
      type(bounded_gmres_result_type) :: result
      real(wp), allocatable :: jac(:, :), b(:), z(:), ref(:), v(:), out(:), rv(:), pz(:), expect(:), dcdt(:, :), d(:, :)
      integer :: i, ncell

      call aux_fixture(grid, aux, 4, coarse_n, 128, 1.0e-14_wp, err)
      call check(error, .not. allocated(err), "fixture failed")
      if (allocated(error)) return
      call check(error, aux%lattice_size(1) == 4 .and. (aux%nested() .eqv. (coarse_n == 2)), "lattice sizes")
      if (allocated(error)) return
      ncell = 4**3
      ! The transfers and the frozen derivative of the preconditioner, from the same inputs
      call deposit%build(grid%xyz, grid%w, 4, aux%side(), err, 1)
      call check(error, .not. allocated(err), "deposit failed")
      if (allocated(error)) return
      call derivative(grid, dcdt)
      allocate (d(ncell, 2))
      call deposit%restrict(dcdt, d, 1)
      call dense_lattice_jacobian(4, aux%side(), d, jac)
      ! The lattice operator alone
      call new_uv3d_aux_level(level, 4, aux%side()/4.0_wp, deposit%origin, source, 1.7_wp, 1, err)
      if (.not. allocated(err)) call level%set_derivative(d, err)
      call check(error, .not. allocated(err), "lattice operator failed")
      if (allocated(error)) return
      allocate (b(2*ncell), z(2*ncell))
      do i = 1, 2*ncell
         b(i) = sin(0.37_wp*real(i, wp)) + 0.1_wp
      end do
      call gmres%setup(2*ncell, 128, 1)
      call gmres%solve(level, b, z, 1.0e-14_wp, result)
      call dense_solve(jac, b, ref, error)
      if (allocated(error)) return
      call check(error, result%accepted, "lattice solve not accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(z - ref)) <= 1.0e-9_wp*maxval(abs(ref)), "lattice solve off the dense inverse")
      if (allocated(error)) return
      ! Molecular action, through the nested solve on that path
      call molecular_vector(grid, v)
      out = v
      call aux%apply(out)
      allocate (rv(2*ncell), pz(2*grid%ngrid))
      call restrict_columns(deposit, grid%ngrid, ncell, v, rv)
      call dense_solve(jac, rv, ref, error)
      if (allocated(error)) return
      call interpolate_columns(deposit, grid%ngrid, ncell, ref - rv, pz)
      expect = v + pz
      call check(error, maxval(abs(out - expect)) <= 1.0e-9_wp*maxval(abs(expect)), "molecular action off")
      if (allocated(error)) return
      stats = aux%statistics()
      if (coarse_n == 2) then
         call check(error, stats%nsolve(2) == stats%niteration(1) .and. stats%nsolve(2) > 0, &
            & "each fine Krylov step must solve on the coarse lattice")
      else
         call check(error, stats%nsolve(2) == 0, "a direct lattice solve used a coarse lattice")
      end if
   end subroutine dense_case

   !> R of two molecular columns
   !>
   !> @param[in]  deposit  transfers
   !> @param[in]  ngrid    points
   !> @param[in]  ncell    cells
   !> @param[in]  v        molecular columns (ngrid, 2)
   !> @param[out] rv       lattice columns (ncell, 2)
   subroutine restrict_columns(deposit, ngrid, ncell, v, rv)
      !> Transfers
      type(grid_deposit_type), intent(in) :: deposit
      !> Points and cells
      integer, intent(in) :: ngrid, ncell
      !> Molecular columns (ngrid, 2)
      real(wp), intent(in) :: v(ngrid, 2)
      !> Lattice columns (ncell, 2)
      real(wp), intent(out) :: rv(ncell, 2)
      call deposit%restrict(v, rv, 1)
   end subroutine restrict_columns

   !> P of two lattice columns
   !>
   !> @param[in]  deposit  transfers
   !> @param[in]  ngrid    points
   !> @param[in]  ncell    cells
   !> @param[in]  z        lattice columns (ncell, 2)
   !> @param[out] pz       molecular columns (ngrid, 2)
   subroutine interpolate_columns(deposit, ngrid, ncell, z, pz)
      !> Transfers
      type(grid_deposit_type), intent(in) :: deposit
      !> Points and cells
      integer, intent(in) :: ngrid, ncell
      !> Lattice columns (ncell, 2)
      real(wp), intent(in) :: z(ncell, 2)
      !> Molecular columns (ngrid, 2)
      real(wp), intent(out) :: pz(ngrid, 2)
      call deposit%interpolate(z, pz, 1)
   end subroutine interpolate_columns

   !> Direct path against the dense inverse
   subroutine test_dense_direct(error)
      type(error_type), allocatable, intent(out) :: error
      call dense_case(4, error)
   end subroutine test_dense_direct

   !> Nested path against the dense inverse
   subroutine test_dense_nested(error)
      type(error_type), allocatable, intent(out) :: error
      call dense_case(2, error)
   end subroutine test_dense_nested

   !> A zero closure derivative gives the identity Jacobian and leaves v unchanged
   subroutine test_identity_derivative(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: grid
      type(uv3d_aux_preconditioner_type) :: aux
      type(uv3d_aux_statistics_type) :: stats
      real(wp), allocatable :: v(:), out(:)
      integer :: coarse(2), ic
      character(len=120) :: text

      coarse = [8, 4]
      do ic = 1, 2
         call aux_fixture(grid, aux, 8, coarse(ic), 2, 0.1_wp, err, zero_derivative=.true.)
         call check(error, .not. allocated(err), "fixture failed")
         if (allocated(error)) return
         call molecular_vector(grid, v)
         out = v
         call aux%apply(out)
         stats = aux%statistics()
         call check(error, all(out == v), "a zero derivative changed v")
         if (allocated(error)) return
         write (text, "(a,5(1x,i0))") "identity: one accepted solve without steps; got", stats%napply, &
            & stats%nsolve(1), stats%nfallback(1), stats%niteration(1), stats%nproduct(1)
         call check(error, stats%napply == 1 .and. stats%nsolve(1) == 1 .and. stats%nfallback(1) == 0 .and. &
            & stats%niteration(1) == 0, trim(text))
         if (allocated(error)) return
      end do
   end subroutine test_identity_derivative

   !> A molecular vector with R v = 0 passes unchanged, on both paths
   subroutine test_unseen_modes(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: grid
      type(uv3d_aux_preconditioner_type) :: aux
      type(grid_deposit_type) :: deposit
      real(wp), allocatable :: u(:), v(:), out(:), ru(:), y(:), py(:), rp(:, :), unit(:), col(:), rv(:)
      integer, allocatable :: occupied(:)
      integer :: coarse(2), ic, ncell, nocc, j, ngrid

      coarse = [4, 2]
      do ic = 1, 2
         call aux_fixture(grid, aux, 4, coarse(ic), 2, 0.1_wp, err)
         call check(error, .not. allocated(err), "fixture failed")
         if (allocated(error)) return
         ncell = aux%lattice_size(1)**3
         ngrid = grid%ngrid
         ! The preconditioner's transfers, from the same inputs
         call deposit%build(grid%xyz, grid%w, aux%lattice_size(1), aux%side(), err, 1)
         call check(error, .not. allocated(err), "deposit failed")
         if (allocated(error)) return
         occupied = pack([(j, j=1, ncell)], deposit%mass > 0.0_wp)
         nocc = size(occupied)
         ! R P on the occupied cells, column by column
         allocate (rp(nocc, nocc), unit(2*ncell), col(2*ngrid), rv(2*ncell))
         do j = 1, nocc
            unit = 0.0_wp
            unit(occupied(j)) = 1.0_wp
            call interpolate_columns(deposit, ngrid, ncell, unit, col)
            call restrict_columns(deposit, ngrid, ncell, col, rv)
            rp(:, j) = rv(occupied)
         end do
         ! v = u - P (R P)^{-1} R u on the first column, 0 in the second
         call molecular_vector(grid, u)
         u(ngrid + 1:) = 0.0_wp
         allocate (ru(2*ncell))
         call restrict_columns(deposit, ngrid, ncell, u, ru)
         call dense_solve(rp, ru(occupied), y, error)
         if (allocated(error)) return
         unit = 0.0_wp
         unit(occupied) = y
         allocate (py(2*ngrid))
         call interpolate_columns(deposit, ngrid, ncell, unit, py)
         v = u - py
         call restrict_columns(deposit, ngrid, ncell, v, rv)
         call check(error, maxval(abs(rv)) <= 1.0e-12_wp*maxval(abs(u)) .and. maxval(abs(v)) > 0.1_wp, &
            & "test vector not in the kernel of R")
         if (allocated(error)) return
         out = v
         call aux%apply(out)
         call check(error, maxval(abs(out - v)) <= 1.0e-10_wp*maxval(abs(v)), "a vector R does not see changed")
         if (allocated(error)) return
         deallocate (rp, unit, col, rv, ru, py)
      end do
   end subroutine test_unseen_modes

   !> A NaN in the closure derivative fails the lattice solve: v unchanged, one fallback
   subroutine test_fallback(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: grid
      type(uv3d_aux_preconditioner_type) :: aux
      type(uv3d_aux_statistics_type) :: stats
      real(wp), allocatable :: v(:), out(:), dcdt(:, :)

      call aux_fixture(grid, aux, 8, 4, 2, 0.1_wp, err)
      call check(error, .not. allocated(err), "fixture failed")
      if (allocated(error)) return
      call derivative(grid, dcdt)
      dcdt(1, 1) = ieee_value(1.0_wp, ieee_quiet_nan)
      call aux%refresh(dcdt)
      call molecular_vector(grid, v)
      out = v
      call aux%apply(out)
      stats = aux%statistics()
      call check(error, all(out == v) .and. stats%nfallback(1) == 1 .and. all(ieee_is_finite(out)), &
         & "a failed lattice solve did not return v")
   end subroutine test_fallback

   !> More cells than the grid's Fourier modes resolve: capped and reported
   subroutine test_size_cap(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: grid
      type(uv3d_aux_preconditioner_type) :: aux
      type(synthetic_source_type) :: source
      integer :: nmax

      call water_grid(grid, err)
      call check(error, .not. allocated(err), "grid failed")
      if (allocated(error)) return
      nmax = maxval(grid%get_kmodes())
      call aux%setup(grid, source, 1.0_wp, 64, 32, 2, 0.1_wp, 1, err)
      call check(error, .not. allocated(err), "setup failed")
      if (allocated(error)) return
      call check(error, nmax < 64 .and. aux%lattice_size(1) == nmax .and. aux%lattice_size(2) == nmax/2, &
         & "lattice not capped")
      if (allocated(error)) return
      call check(error, index(aux%describe(1.0_wp), "capped") > 0, aux%describe(1.0_wp))
   end subroutine test_size_cap

   !> Setup from refused sizes, inner settings or grids
   !>
   !> - A failing fixture, a refusal for another reason than the expected one, or a
   !>   preconditioner left built return no test error, so they show as an unexpected pass
   !>
   !> @param[out] error    the refusal
   !> @param[in]  variant  input: 1 odd cells, 2 coarse lattice not half, 3 zero inner steps,
   !>                      4 inner tolerance one, 5 negative weight, 6 no Fourier modes,
   !>                      7 a grid never updated
   subroutine setup_refused(error, variant)
      !> The refusal
      type(error_type), allocatable, intent(out) :: error
      !> Input
      integer, intent(in) :: variant
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_3d_cartesian_type) :: box
      type(uv3d_aux_preconditioner_type) :: aux
      type(synthetic_source_type) :: source
      character(len=:), allocatable :: expected

      call water_grid(grid, err, reciprocal=variant /= 6)
      if (allocated(err)) return
      select case (variant)
      case (1)
         call aux%setup(grid, source, 1.0_wp, 7, 7, 2, 0.1_wp, 1, err)
         expected = "even and at least 4"
      case (2)
         call aux%setup(grid, source, 1.0_wp, 16, 4, 2, 0.1_wp, 1, err)
         expected = "coarse cells per axis must equal"
      case (3)
         call aux%setup(grid, source, 1.0_wp, 8, 4, 0, 0.1_wp, 1, err)
         expected = "inner steps must be at least 1"
      case (4)
         call aux%setup(grid, source, 1.0_wp, 8, 4, 2, 1.0_wp, 1, err)
         expected = "inner steps must be at least 1"
      case (5)
         grid%w(5) = -abs(grid%w(5))
         call aux%setup(grid, source, 1.0_wp, 8, 4, 2, 0.1_wp, 1, err)
         expected = "quadrature weights must be positive"
      case (7)
         call new_cartesian_point_grid(box, err, 8, 8, 8, 0.5_wp)
         if (allocated(err)) return
         call aux%setup(box, source, 1.0_wp, 8, 4, 2, 0.1_wp, 1, err)
         expected = "needs an updated grid"
      case default
         call aux%setup(grid, source, 1.0_wp, 8, 4, 2, 0.1_wp, 1, err)
         expected = "needs a grid with Fourier modes"
      end select
      if (aux%built()) return
      call expect_moist_error(error, err, expected)
   end subroutine setup_refused

   !> Refused: setup with odd cells per axis
   subroutine test_refused_odd_cells(error)
      type(error_type), allocatable, intent(out) :: error
      call setup_refused(error, 1)
   end subroutine test_refused_odd_cells

   !> Refused: setup with a coarse lattice other than half the fine one
   subroutine test_refused_coarse_size(error)
      type(error_type), allocatable, intent(out) :: error
      call setup_refused(error, 2)
   end subroutine test_refused_coarse_size

   !> Refused: setup with zero inner steps
   subroutine test_refused_inner_steps(error)
      type(error_type), allocatable, intent(out) :: error
      call setup_refused(error, 3)
   end subroutine test_refused_inner_steps

   !> Refused: setup with an inner tolerance of one
   subroutine test_refused_inner_rtol(error)
      type(error_type), allocatable, intent(out) :: error
      call setup_refused(error, 4)
   end subroutine test_refused_inner_rtol

   !> Refused: setup with a negative quadrature weight
   subroutine test_refused_negative_weight(error)
      type(error_type), allocatable, intent(out) :: error
      call setup_refused(error, 5)
   end subroutine test_refused_negative_weight

   !> Refused: setup with a grid without Fourier modes
   subroutine test_refused_no_modes(error)
      type(error_type), allocatable, intent(out) :: error
      call setup_refused(error, 6)
   end subroutine test_refused_no_modes

   !> Refused: setup on a Cartesian grid that was never updated
   subroutine test_refused_unbuilt_grid(error)
      type(error_type), allocatable, intent(out) :: error
      call setup_refused(error, 7)
   end subroutine test_refused_unbuilt_grid

end module test_moz_3d_preconditioner
