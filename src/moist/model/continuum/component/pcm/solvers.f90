!> PCM linear solvers
!>
!> - Operator and solver types for pcmgen-generated kernels: pcm_operator,
!>   pcm_complex_operator, pcm_solver, pcm_dense_solver, pcm_gmres_solver and
!>   the preconditioners pcm_preconditioner and pcm_block_jacobi (see pcmgen
!>   docs/moist-solver-interface.md)
!> - pcm_dense_solver factors each operator once (Cholesky, Bunch-Kaufman or
!>   LU) and reuses the factors for both orientations and many right-hand sides
!> - pcm_gmres_solver runs restarted GMRES with right preconditioning on the
!>   operator applications; pcm_block_jacobi factors the diagonal blocks of
!>   point clusters found by recursive coordinate bisection
!> - Solutions are updated in place: x holds the initial guess on entry (used
!>   by iterative solvers, ignored by direct ones) and the solution on exit
!> - Complex arithmetic is holomorphic: complex operators are applied and
!>   solved with the plain transpose, never the conjugate transpose, and a
!>   symmetric complex operator equals its transpose (not Hermitian); only the
!>   Krylov inner products of GMRES are Hermitian
!> - solve_complex with a real operator and a complex right-hand side solves
!>   the real and imaginary parts as two real right-hand sides (every solver)
!> - solve_complex with a complex operator: pcm_dense_solver factors with
!>   zsytrf (symmetric) or zgetrf, pcm_gmres_solver iterates in complex
!>   arithmetic
!> - Matrix-based routines used by the native PCM code: solve_pcm_lu,
!>   solve_pcm_cholesky, solve_pcm_iterative and solve_pcm_inversion
!> - LAPACK and BLAS are called only through moist_math_lapack and
!>   moist_math_blas, so ILP64 builds keep working
!> - pcm_nonlinear_problem: nonlinear systems of the generated kernels, solved
!>   by Newton's method with an Armijo backtracking line search; the kernel
!>   computes each step with the solver it points to
module moist_model_continuum_component_pcm_solvers
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int64, output_unit
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use moist_math_lapack, only: getrf, getrs, getri, potrf, potrs, sytrf, sytrs
   use moist_math_blas, only: dot, gemv, nrm2
   implicit none(type, external)
   private

   public :: solve_pcm_lu
   public :: solve_pcm_cholesky
   public :: solve_pcm_iterative
   public :: solve_pcm_inversion
   public :: pcm_operator, pcm_complex_operator
   public :: pcm_solver, pcm_dense_solver, pcm_gmres_solver
   public :: pcm_nonlinear_problem
   public :: pcm_preconditioner, pcm_block_jacobi
   public :: pcm_factor_cholesky, pcm_factor_symmetric, pcm_factor_lu
   public :: pcm_factor_complex_symmetric, pcm_factor_complex_lu
   public :: pcm_storage_dense, pcm_storage_matrix_free

   !> Factorization of a symmetric positive definite operator (potrf, lower)
   integer, parameter :: pcm_factor_cholesky = 1
   !> Factorization of a symmetric indefinite operator (sytrf, lower)
   integer, parameter :: pcm_factor_symmetric = 2
   !> Factorization of a general operator (getrf)
   integer, parameter :: pcm_factor_lu = 3
   !> Factorization of a complex symmetric operator (zsytrf, lower)
   integer, parameter :: pcm_factor_complex_symmetric = 4
   !> Factorization of a complex general operator (zgetrf)
   integer, parameter :: pcm_factor_complex_lu = 5

   !> Operator storage of the generated kernels: dense matrices
   integer, parameter :: pcm_storage_dense = 1
   !> Operator storage of the generated kernels: matrix-free applications
   integer, parameter :: pcm_storage_matrix_free = 2

   !> Sufficient-decrease constant c of the Newton line search on |F|^2
   real(wp), parameter :: newton_armijo_c = 1.0e-4_wp
   !> Shortest Newton step length tried before the iteration counts as stalled
   real(wp), parameter :: newton_step_min = 2.0_wp**(-20)

   !> Per-id work of an iterative solver or preconditioner: none cached
   integer, parameter :: slot_empty = 0
   !> Per-id work of an iterative solver or preconditioner: real operator
   integer, parameter :: slot_real = 1
   !> Per-id work of an iterative solver or preconditioner: complex operator
   integer, parameter :: slot_complex = 2

   !> Square linear operator of a generated PCM model
   !>
   !> - Valid only during the evaluate call that created it
   !> - Extensions hold pointers to the caller's inputs and intermediates
   !> - A block operator stacks m fields of npoint surface points each
   !>   (npt = m npoint); unknown (a - 1) npoint + i is field a at point i
   type, abstract :: pcm_operator
      !> Dimension of the operator (m npoint)
      integer :: npt = 0
      !> Number of surface points; 0 means npt (a single field)
      integer :: npoint = 0
      !> Operator index within the generated module, 1-based
      integer :: id = 0
      !> The operator equals its transpose (shown by the generator)
      logical :: symmetric = .false.
      !> The operator is symmetric positive definite (declared in the model)
      logical :: positive_definite = .false.
      !> Positions of the surface points, (3, npoint); null if the model has none
      real(wp), pointer, contiguous :: points(:, :) => null()
   contains
      !> Apply the operator or its transpose to a vector
      procedure(apply_iface), deferred :: apply
      !> Form the dense operator matrix
      procedure(assemble_iface), deferred :: assemble
      !> Form the operator of a sub-surface (near-field block)
      procedure(assemble_block_iface), deferred :: assemble_block
   end type pcm_operator

   !> Dense matrix adapter for the legacy iterative entry point
   type, extends(pcm_operator) :: pcm_matrix_operator
      real(wp), allocatable :: a(:, :)
   contains
      procedure :: apply => pcm_matrix_apply
      procedure :: assemble => pcm_matrix_assemble
      procedure :: assemble_block => pcm_matrix_assemble_block
   end type pcm_matrix_operator

   !> Deferred procedures of pcm_operator
   abstract interface
      !> Apply an operator or its transpose
      !>
      !> @param[in]  self   operator
      !> @param[in]  trans  apply the transpose if true
      !> @param[in]  x      vector of length npt
      !> @param[out] y      vector of length npt
      subroutine apply_iface(self, trans, x, y)
         import :: pcm_operator, wp
         implicit none(type, external)
         !> Operator
         class(pcm_operator), intent(in) :: self
         !> Apply the transpose
         logical, intent(in) :: trans
         !> Input vector
         real(wp), contiguous, intent(in) :: x(:)
         !> Result vector
         real(wp), contiguous, intent(out) :: y(:)
      end subroutine apply_iface

      !> Form the dense operator matrix
      !>
      !> @param[in]  self   operator
      !> @param[out] a      matrix of shape (npt, npt)
      !> @param[out] error  allocated if workspace allocation failed
      subroutine assemble_iface(self, a, error)
         import :: pcm_operator, wp, error_type
         implicit none(type, external)
         !> Operator
         class(pcm_operator), intent(in) :: self
         !> Dense operator matrix
         real(wp), contiguous, intent(out) :: a(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine assemble_iface

      !> Form the operator of the sub-surface of the points idx
      !>
      !> - Every kernel, diagonal and point vector of the operator expression is
      !>   restricted to idx; this is the principal submatrix only for operators
      !>   without matrix products, else the near-field approximation of it
      !> - Block (a, b) of a block operator occupies rows (a - 1) nb + 1 .. a nb
      !>   and the matching columns, points in the order of idx
      !>
      !> @param[in]  self   operator
      !> @param[in]  idx    nb distinct point indices in 1..npoint
      !> @param[out] a      matrix of shape (m nb, m nb)
      !> @param[out] error  allocated if workspace allocation failed
      subroutine assemble_block_iface(self, idx, a, error)
         import :: pcm_operator, wp, error_type
         implicit none(type, external)
         !> Operator
         class(pcm_operator), intent(in) :: self
         !> Point indices of the sub-surface
         integer, contiguous, intent(in) :: idx(:)
         !> Operator matrix of the sub-surface
         real(wp), contiguous, intent(out) :: a(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine assemble_block_iface
   end interface

   !> Square complex linear operator of a generated PCM model
   !>
   !> - Valid only during the evaluate call that created it
   !> - Extensions hold pointers to the caller's inputs and intermediates
   !> - Holomorphic: trans applies the plain transpose, never the conjugate one
   type, abstract :: pcm_complex_operator
      !> Dimension of the operator (m npoint)
      integer :: npt = 0
      !> Number of surface points; 0 means npt (a single field)
      integer :: npoint = 0
      !> Operator index within the generated module, 1-based
      integer :: id = 0
      !> The operator equals its plain transpose (shown by the generator)
      logical :: symmetric = .false.
      !> Positions of the surface points, (3, npoint); null if the model has none
      real(wp), pointer, contiguous :: points(:, :) => null()
   contains
      !> Apply the operator or its transpose to a vector
      procedure(complex_apply_iface), deferred :: apply
      !> Form the dense operator matrix
      procedure(complex_assemble_iface), deferred :: assemble
      !> Form the operator of a sub-surface (near-field block)
      procedure(complex_assemble_block_iface), deferred :: assemble_block
   end type pcm_complex_operator

   !> Deferred procedures of pcm_complex_operator
   abstract interface
      !> Apply a complex operator or its plain transpose
      !>
      !> @param[in]  self   operator
      !> @param[in]  trans  apply the transpose (not conjugated) if true
      !> @param[in]  x      vector of length npt
      !> @param[out] y      vector of length npt
      subroutine complex_apply_iface(self, trans, x, y)
         import :: pcm_complex_operator, wp
         implicit none(type, external)
         !> Operator
         class(pcm_complex_operator), intent(in) :: self
         !> Apply the transpose
         logical, intent(in) :: trans
         !> Input vector
         complex(wp), contiguous, intent(in) :: x(:)
         !> Result vector
         complex(wp), contiguous, intent(out) :: y(:)
      end subroutine complex_apply_iface

      !> Form the dense complex operator matrix
      !>
      !> @param[in]  self   operator
      !> @param[out] a      matrix of shape (npt, npt)
      !> @param[out] error  allocated if workspace allocation failed
      subroutine complex_assemble_iface(self, a, error)
         import :: pcm_complex_operator, wp, error_type
         implicit none(type, external)
         !> Operator
         class(pcm_complex_operator), intent(in) :: self
         !> Dense operator matrix
         complex(wp), contiguous, intent(out) :: a(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine complex_assemble_iface

      !> Form the complex operator of the sub-surface of the points idx
      !>
      !> - Same restriction and block layout as assemble_block of pcm_operator
      !>
      !> @param[in]  self   operator
      !> @param[in]  idx    nb distinct point indices in 1..npoint
      !> @param[out] a      matrix of shape (m nb, m nb)
      !> @param[out] error  allocated if workspace allocation failed
      subroutine complex_assemble_block_iface(self, idx, a, error)
         import :: pcm_complex_operator, wp, error_type
         implicit none(type, external)
         !> Operator
         class(pcm_complex_operator), intent(in) :: self
         !> Point indices of the sub-surface
         integer, contiguous, intent(in) :: idx(:)
         !> Operator matrix of the sub-surface
         complex(wp), contiguous, intent(out) :: a(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine complex_assemble_block_iface
   end interface

   !> Linear solver used by generated PCM code
   !>
   !> - reset is called at the start of every assemble phase and discards the
   !>   work cached per operator id
   !> - setup (setup_complex for complex operators) is called for an operator
   !>   before any solve with it (factorization, preconditioner) and again
   !>   whenever the operator of that id changed, e.g. by a repeated medium
   !>   phase with new solvent parameters; it always discards the work cached
   !>   for that id and prepares the operator anew
   !> - Solves prepare an operator on first use if nothing is cached for its id
   !> - x holds the initial guess on entry and the solution on exit; direct
   !>   solvers ignore the guess
   !> - solve_complex with a real operator solves the real and imaginary parts
   !>   through solve_block; with a complex operator it is an error unless an
   !>   extension overrides solve_complex_complex
   !> - inexact is true for solvers whose solutions satisfy the system only to a
   !>   tolerance (the generated code then adds an adjoint energy correction)
   !> - matrix_free asks whether the generated kernels should apply operators
   !>   without storing their matrices; direct solvers need the matrices
   !> - Outputs are undefined when error is allocated
   type, abstract :: pcm_solver
      !> Relative residual of nonlinear systems: |F(x)| <= nonlinear_tol |F(0)|
      real(wp) :: nonlinear_tol = 1.0e-10_wp
      !> Maximum number of Newton steps per nonlinear system
      integer :: nonlinear_max_iter = 50
      !> Newton steps of the last nonlinear system
      integer :: nonlinear_iterations = 0
      !> Relative residual |F(x)|/|F(0)| of the last nonlinear system
      real(wp) :: nonlinear_residual = 0.0_wp
   contains
      !> Discard state cached for the operators of a previous phase
      procedure :: reset => pcm_solver_reset
      !> Prepare an operator for solves
      procedure :: setup => pcm_solver_setup
      !> Prepare a complex operator for solves
      procedure :: setup_complex => pcm_solver_setup_complex
      !> Solve op x = rhs or op**T x = rhs for one right-hand side
      procedure :: solve => pcm_solver_solve
      !> Solve op x = rhs or op**T x = rhs for several right-hand sides
      procedure(solve_block_iface), deferred :: solve_block
      !> Solve with a real operator and a complex right-hand side
      procedure :: solve_complex_real => pcm_solver_solve_complex_real
      !> Solve with a complex operator and a complex right-hand side
      procedure :: solve_complex_complex => pcm_solver_solve_complex_complex
      !> Solve op x = rhs or op**T x = rhs for one complex right-hand side
      generic :: solve_complex => solve_complex_real, solve_complex_complex
      !> Solutions satisfy the system only to a tolerance
      procedure :: inexact => pcm_solver_inexact
      !> The generated kernels should not store the operator matrices
      procedure :: matrix_free => pcm_solver_matrix_free
   end type pcm_solver

   !> Nonlinear system F(x) = 0 of a generated kernel
   !>
   !> - The kernel evaluates the residual and computes Newton steps
   !>   p = -J(x)^-1 F(x), setting up its Jacobian operator with the solver it
   !>   points to; solve runs Newton's method with an Armijo line search
   !> - Converged when |F(x)| <= solver%nonlinear_tol |F(0)|, in at most
   !>   solver%nonlinear_max_iter steps; the statistics go to the solver
   !> - Components point to the caller's data; valid only during that call
   type, abstract :: pcm_nonlinear_problem
      !> Number of unknowns
      integer :: npt = 0
      !> Linear solver of the Newton steps; its controls and statistics
      class(pcm_solver), pointer :: solver => null()
   contains
      !> Residual F(x)
      procedure(nonlinear_residual_iface), deferred :: residual
      !> Newton step p = -J(x)^-1 F(x)
      procedure(nonlinear_step_iface), deferred :: step
      !> Solve F(x) = 0 from the initial guess in x
      procedure :: solve => pcm_nonlinear_solve
   end type pcm_nonlinear_problem

   !> Deferred procedures of pcm_nonlinear_problem
   abstract interface
      !> Residual of a nonlinear system
      !>
      !> @param[in]  self   problem
      !> @param[in]  x      unknowns, length npt
      !> @param[out] f      residual, length npt
      !> @param[out] error  allocated if a workspace allocation failed
      subroutine nonlinear_residual_iface(self, x, f, error)
         import :: pcm_nonlinear_problem, wp, error_type
         implicit none(type, external)
         !> Problem
         class(pcm_nonlinear_problem), intent(in) :: self
         !> Unknowns
         real(wp), contiguous, intent(in) :: x(:)
         !> Residual
         real(wp), contiguous, intent(out) :: f(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine nonlinear_residual_iface

      !> Newton step of a nonlinear system
      !>
      !> @param[in]  self   problem
      !> @param[in]  x      unknowns, length npt
      !> @param[in]  f      residual at x, length npt
      !> @param[out] p      step -J(x)^-1 f, length npt
      !> @param[out] error  allocated if the linear solve or an allocation failed
      subroutine nonlinear_step_iface(self, x, f, p, error)
         import :: pcm_nonlinear_problem, wp, error_type
         implicit none(type, external)
         !> Problem
         class(pcm_nonlinear_problem), intent(in) :: self
         !> Unknowns
         real(wp), contiguous, intent(in) :: x(:)
         !> Residual at x
         real(wp), contiguous, intent(in) :: f(:)
         !> Newton step
         real(wp), contiguous, intent(out) :: p(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine nonlinear_step_iface
   end interface

   !> Deferred procedure of pcm_solver
   abstract interface
      !> Solve a linear system with several right-hand sides
      !>
      !> @param[in,out] self   solver state
      !> @param[in]     op     operator
      !> @param[in]     trans  solve with the transpose if true
      !> @param[in]     rhs    right-hand sides, (npt, nrhs)
      !> @param[in,out] x      initial guesses on entry, solutions on exit, (npt, nrhs)
      !> @param[out]    error  allocated on failure
      subroutine solve_block_iface(self, op, trans, rhs, x, error)
         import :: pcm_solver, pcm_operator, wp, error_type
         implicit none(type, external)
         !> Solver state
         class(pcm_solver), intent(inout) :: self
         !> Operator
         class(pcm_operator), intent(in) :: op
         !> Solve with the transpose
         logical, intent(in) :: trans
         !> Right-hand sides
         real(wp), contiguous, intent(in) :: rhs(:, :)
         !> Initial guesses on entry, solutions on exit
         real(wp), contiguous, intent(inout) :: x(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_block_iface
   end interface

   !> Dense factors of one operator
   type :: pcm_dense_factors
      !> Factors are available
      logical :: ready = .false.
      !> Factorization: one of the pcm_factor_* kinds
      integer :: kind = 0
      !> Real factors from potrf (lower), sytrf (lower) or getrf, (npt, npt)
      real(wp), allocatable :: a(:, :)
      !> Complex factors from zsytrf (lower) or zgetrf, (npt, npt)
      complex(wp), allocatable :: z(:, :)
      !> Pivot indices from sytrf, getrf, zsytrf or zgetrf, (npt)
      integer, allocatable :: ipiv(:)
   end type pcm_dense_factors

   !> Dense direct solver: assemble and factor once per operator
   !>
   !> - Cholesky for positive definite operators, Bunch-Kaufman for other
   !>   symmetric ones, else LU
   !> - Complex operators: complex Bunch-Kaufman (zsytrf) if symmetric, else
   !>   complex LU (zgetrf); transposed solves use the plain transpose
   !> - A failed Cholesky factorization is an error, never a fallback to LU
   !> - The factors serve both orientations and any number of right-hand sides
   !> - Real and complex operators share the slots indexed by id; setup and
   !>   setup_complex replace whatever is cached, a solve on a slot holding
   !>   factors of the other kind is an error
   type, extends(pcm_solver) :: pcm_dense_solver
      !> Factors indexed by operator id
      type(pcm_dense_factors), allocatable, private :: factors(:)
      !> Number of factorizations since construction
      integer :: n_factorizations = 0
      !> Number of right-hand sides solved since construction
      integer :: n_solves = 0
   contains
      !> Discard factors of a previous phase
      procedure :: reset => pcm_dense_reset
      !> Assemble and factor an operator
      procedure :: setup => pcm_dense_setup
      !> Solve with the cached factors
      procedure :: solve_block => pcm_dense_solve_block
      !> Assemble and factor a complex operator
      procedure :: setup_complex => pcm_dense_setup_complex
      !> Solve with the cached complex factors
      procedure :: solve_complex_complex => pcm_dense_solve_complex_complex
      !> Factorization used for an operator id, 0 if not factored
      procedure :: factorization => pcm_dense_factorization
   end type pcm_dense_solver

   !> Preconditioner M of an iterative solver, applied as y = M^-1 x
   !>
   !> - Work is cached per operator id in slots created by reset, like the
   !>   factors of pcm_dense_solver
   !> - setup (setup_complex) always discards what is cached for the id of the
   !>   operator and builds the preconditioner anew
   !> - apply (apply_complex) with trans applies M^-T, the preconditioner of
   !>   op**T (plain transpose for complex operators)
   type, abstract :: pcm_preconditioner
   contains
      !> Discard the work of a previous phase and create n_operators slots
      procedure :: reset => pcm_preconditioner_reset
      !> Build the preconditioner of an operator
      procedure(preconditioner_setup_iface), deferred :: setup
      !> Build the preconditioner of a complex operator
      procedure(preconditioner_setup_complex_iface), deferred :: setup_complex
      !> Apply M^-1 or M^-T
      procedure(preconditioner_apply_iface), deferred :: apply
      !> Apply M^-1 or M^-T of a complex operator
      procedure(preconditioner_apply_complex_iface), deferred :: apply_complex
   end type pcm_preconditioner

   !> Deferred procedures of pcm_preconditioner
   abstract interface
      !> Build the preconditioner of an operator, replacing the one of its id
      !>
      !> @param[in,out] self   preconditioner state
      !> @param[in]     op     operator
      !> @param[out]    error  allocated on failure
      subroutine preconditioner_setup_iface(self, op, error)
         import :: pcm_preconditioner, pcm_operator, error_type
         implicit none(type, external)
         !> Preconditioner state
         class(pcm_preconditioner), intent(inout) :: self
         !> Operator
         class(pcm_operator), intent(in) :: op
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine preconditioner_setup_iface

      !> Build the preconditioner of a complex operator, replacing the one of its id
      !>
      !> @param[in,out] self   preconditioner state
      !> @param[in]     op     complex operator
      !> @param[out]    error  allocated on failure
      subroutine preconditioner_setup_complex_iface(self, op, error)
         import :: pcm_preconditioner, pcm_complex_operator, error_type
         implicit none(type, external)
         !> Preconditioner state
         class(pcm_preconditioner), intent(inout) :: self
         !> Complex operator
         class(pcm_complex_operator), intent(in) :: op
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine preconditioner_setup_complex_iface

      !> Apply the preconditioner of an operator: y = M^-1 x or M^-T x
      !>
      !> @param[in]  self   preconditioner state
      !> @param[in]  op     operator the preconditioner was set up for
      !> @param[in]  trans  apply M^-T if true
      !> @param[in]  x      vector of length npt
      !> @param[out] y      vector of length npt
      !> @param[out] error  allocated if the operator is not set up
      subroutine preconditioner_apply_iface(self, op, trans, x, y, error)
         import :: pcm_preconditioner, pcm_operator, wp, error_type
         implicit none(type, external)
         !> Preconditioner state
         class(pcm_preconditioner), intent(in) :: self
         !> Operator
         class(pcm_operator), intent(in) :: op
         !> Apply the transpose
         logical, intent(in) :: trans
         !> Input vector
         real(wp), contiguous, intent(in) :: x(:)
         !> Result vector
         real(wp), contiguous, intent(out) :: y(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine preconditioner_apply_iface

      !> Apply the preconditioner of a complex operator: y = M^-1 x or M^-T x
      !>
      !> @param[in]  self   preconditioner state
      !> @param[in]  op     complex operator the preconditioner was set up for
      !> @param[in]  trans  apply M^-T (not conjugated) if true
      !> @param[in]  x      vector of length npt
      !> @param[out] y      vector of length npt
      !> @param[out] error  allocated if the operator is not set up
      subroutine preconditioner_apply_complex_iface(self, op, trans, x, y, error)
         import :: pcm_preconditioner, pcm_complex_operator, wp, error_type
         implicit none(type, external)
         !> Preconditioner state
         class(pcm_preconditioner), intent(in) :: self
         !> Complex operator
         class(pcm_complex_operator), intent(in) :: op
         !> Apply the transpose
         logical, intent(in) :: trans
         !> Input vector
         complex(wp), contiguous, intent(in) :: x(:)
         !> Result vector
         complex(wp), contiguous, intent(out) :: y(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine preconditioner_apply_complex_iface
   end interface

   !> LU factors of one diagonal block of pcm_block_jacobi
   type :: pcm_jacobi_block
      !> Unknowns of the block in the order of its rows, (m nb)
      integer, allocatable :: unknowns(:)
      !> Real LU factors from getrf, (m nb, m nb)
      real(wp), allocatable :: a(:, :)
      !> Complex LU factors from zgetrf, (m nb, m nb)
      complex(wp), allocatable :: z(:, :)
      !> Pivot indices, (m nb)
      integer, allocatable :: ipiv(:)
   end type pcm_jacobi_block

   !> Block-Jacobi preconditioner of one operator
   type :: pcm_jacobi_slot
      !> Kind of the factors: slot_empty, slot_real or slot_complex
      integer :: kind = slot_empty
      !> Operator dimension the factors were built for
      integer :: npt = 0
      !> Diagonal blocks, one per point cluster
      type(pcm_jacobi_block), allocatable :: blocks(:)
   end type pcm_jacobi_slot

   !> Block-Jacobi preconditioner on clusters of surface points
   !>
   !> - Clusters of at most block_size points by recursive coordinate
   !>   bisection of op%points (split across the longest extent); contiguous
   !>   index ranges if op%points is null
   !> - Each diagonal block is formed with op%assemble_block and LU-factored
   !>   (getrf, zgetrf); M^-T uses getrs with trans = "t"
   !> - Block operators (npt = m npoint): a cluster couples all m fields of its
   !>   points, unknowns (a - 1) npoint + idx
   !> - block_size 1 is Jacobi; setup costs O(npt block_size^2 m), memory
   !>   O(npt block_size m)
   type, extends(pcm_preconditioner) :: pcm_block_jacobi
      !> Maximum number of points per cluster
      integer :: block_size = 64
      !> Factors indexed by operator id
      type(pcm_jacobi_slot), allocatable, private :: slots(:)
   contains
      !> Discard the factors of a previous phase
      procedure :: reset => pcm_block_jacobi_reset
      !> Cluster the points and factor the diagonal blocks
      procedure :: setup => pcm_block_jacobi_setup
      !> Cluster the points and factor the complex diagonal blocks
      procedure :: setup_complex => pcm_block_jacobi_setup_complex
      !> Solve with the diagonal blocks
      procedure :: apply => pcm_block_jacobi_apply
      !> Solve with the complex diagonal blocks
      procedure :: apply_complex => pcm_block_jacobi_apply_complex
      !> Number of diagonal blocks of an operator id, 0 if not set up
      procedure :: n_blocks => pcm_block_jacobi_n_blocks
   end type pcm_block_jacobi

   !> Restarted GMRES with right preconditioning
   !>
   !> - Applies the operator with op%apply only; never assembles or factors it
   !> - GMRES(restart): modified Gram-Schmidt Arnoldi, Givens rotations;
   !>   right preconditioning makes the Arnoldi residual the residual of the
   !>   unpreconditioned system, which is checked against the true residual
   !>   norm(rhs - op x) <= tol norm(rhs) at the end of every cycle
   !> - Starts from x (warm start); a zero right-hand side gives x = 0
   !> - max_iter bounds the iterations per right-hand side; non-convergence is
   !>   an error
   !> - Complex operators: the same algorithm in complex arithmetic; the
   !>   operator is applied with the plain transpose, the Krylov inner products
   !>   are Hermitian
   !> - preconditioner: pcm_block_jacobi unless set before reset; setup and
   !>   setup_complex always rebuild it for the id, a solve on an id without
   !>   preconditioner sets it up first
   !> - storage chooses whether the generated kernels store the operator
   !>   matrices: pcm_storage_dense or pcm_storage_matrix_free
   type, extends(pcm_solver) :: pcm_gmres_solver
      !> Relative residual threshold norm(rhs - op x) / norm(rhs)
      real(wp) :: tol = 1.0e-8_wp
      !> Krylov dimension per cycle
      integer :: restart = 100
      !> Maximum number of iterations per right-hand side
      integer :: max_iter = 1000
      !> Preconditioner; pcm_block_jacobi if not allocated at the first setup
      class(pcm_preconditioner), allocatable :: preconditioner
      !> Operator storage of the generated kernels, one of the pcm_storage_* kinds
      integer :: storage = pcm_storage_dense
      !> Number of right-hand sides solved since construction
      integer :: n_solves = 0
      !> Number of iterations since construction
      integer :: n_iterations = 0
      !> Iterations of the last right-hand side
      integer :: last_iterations = 0
      !> Relative true residual of the last right-hand side
      real(wp) :: last_residual = 0.0_wp
      !> Number of operator slots created by the last reset
      integer, private :: n_operators = -1
      !> Kind of the preconditioner per operator id: slot_empty, slot_real or slot_complex
      integer, allocatable, private :: prepared(:)
   contains
      !> Discard the preconditioners of a previous phase
      procedure :: reset => pcm_gmres_reset
      !> Build the preconditioner of an operator
      procedure :: setup => pcm_gmres_setup
      !> Build the preconditioner of a complex operator
      procedure :: setup_complex => pcm_gmres_setup_complex
      !> Solve each right-hand side by GMRES
      procedure :: solve_block => pcm_gmres_solve_block
      !> Solve a complex operator by complex GMRES
      procedure :: solve_complex_complex => pcm_gmres_solve_complex_complex
      !> Solutions satisfy the system to tol only
      procedure :: inexact => pcm_gmres_inexact
      !> Matrix-free kernels by storage
      procedure :: matrix_free => pcm_gmres_matrix_free
   end type pcm_gmres_solver

contains

   !> Solve a nonlinear system F(x) = 0 by Newton's method with an Armijo line search
   !>
   !> - The problem computes each step p = -J(x)^-1 F(x)
   !> - Step lengths t = 1, 1/2, 1/4, ... until
   !>   |F(x + t p)|^2 <= (1 - 2 c t) |F(x)|^2, c = newton_armijo_c
   !> - Converged when |F(x)| <= nonlinear_tol |F(0)| (x = 0 if F(0) = 0)
   !> - Stalled when no t >= newton_step_min decreases |F|: the Newton step is
   !>   far too long, so J(x) is near singular at a point that is not a root
   !>   (typically past a fold of the solution branch); error at once
   !> - Not converged in nonlinear_max_iter steps, or a failed residual or step: error
   !> - Statistics go to the solver: nonlinear_iterations, nonlinear_residual
   !>
   !> @param[in]     self   problem
   !> @param[in,out] x      initial guess on entry, solution on exit
   !> @param[out]    error  allocated if Newton failed
   subroutine pcm_nonlinear_solve(self, x, error)
      !> Problem
      class(pcm_nonlinear_problem), intent(in) :: self
      !> Initial guess on entry, solution on exit
      real(wp), contiguous, intent(inout) :: x(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_nonlinear_solve]"
      real(wp), allocatable :: f(:), p(:), x_trial(:), f_trial(:)
      real(wp) :: scale, norm_f, norm_trial, t
      integer :: stat, steps

      if (.not. associated(self%solver)) then
         call fatal_error(error, here//" the problem has no solver")
         return
      end if
      allocate (f(size(x)), stat=stat)
      if (stat == 0) allocate (p(size(x)), stat=stat)
      if (stat == 0) allocate (x_trial(size(x)), stat=stat)
      if (stat == 0) allocate (f_trial(size(x)), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, here//" allocation failed")
         return
      end if
      x_trial = 0.0_wp
      call self%residual(x_trial, f, error)
      if (allocated(error)) return
      scale = nrm2(f)
      self%solver%nonlinear_iterations = 0
      self%solver%nonlinear_residual = 0.0_wp
      if (.not. scale > 0.0_wp) then
         x = 0.0_wp
         return
      end if

      call self%residual(x, f, error)
      if (allocated(error)) return
      norm_f = nrm2(f)
      steps = 0
      do while (.not. norm_f <= self%solver%nonlinear_tol*scale)
         if (steps >= self%solver%nonlinear_max_iter) then
            call newton_statistics(self%solver, steps, norm_f/scale)
            call fatal_error(error, here//" Newton did not converge: |F|/|F(0)| = " &
               & //real_str(norm_f/scale)//" after "//int_str(steps)//" steps")
            return
         end if
         call self%step(x, f, p, error)
         if (allocated(error)) return
         steps = steps + 1
         t = 1.0_wp
         do
            x_trial = x + t*p
            call self%residual(x_trial, f_trial, error)
            if (allocated(error)) return
            norm_trial = nrm2(f_trial)
            ! Fails for a non-finite trial residual too
            if (norm_trial <= sqrt(1.0_wp - 2.0_wp*newton_armijo_c*t)*norm_f) exit
            t = 0.5_wp*t
            if (t < newton_step_min) then
               call newton_statistics(self%solver, steps, norm_f/scale)
               call fatal_error(error, here//" Newton stalled at |F|/|F(0)| = " &
                  & //real_str(norm_f/scale)//" after "//int_str(steps)//" steps: no step" &
                  & //" length down to "//real_str(newton_step_min)//" decreases |F|; the" &
                  & //" Jacobian is near singular away from a root (fold of the solution branch)")
               return
            end if
         end do
         x = x_trial
         f = f_trial
         norm_f = norm_trial
      end do
      call newton_statistics(self%solver, steps, norm_f/scale)
   end subroutine pcm_nonlinear_solve

   !> Record the statistics of a nonlinear solve
   !>
   !> @param[in,out] solver    solver of the Newton steps
   !> @param[in]     steps     Newton steps taken
   !> @param[in]     residual  final |F|/|F(0)|
   subroutine newton_statistics(solver, steps, residual)
      !> Solver of the Newton steps
      class(pcm_solver), intent(inout) :: solver
      !> Newton steps taken
      integer, intent(in) :: steps
      !> Final relative residual
      real(wp), intent(in) :: residual

      solver%nonlinear_iterations = steps
      solver%nonlinear_residual = residual
   end subroutine newton_statistics

   !> Solve PCM system using LU factorization
   !> Solves A*q = rhs via LAPACK's LU decomposition (DGETRF + DGETRS).
   subroutine solve_pcm_lu(amat, rhs, q, error, unit)
      !> System matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)
      !> Right-hand side vector (ngrid)
      real(wp), intent(in) :: rhs(:)
      !> Solution vector - surface charges (ngrid)
      real(wp), intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Output unit for factorization diagnostics
      integer, intent(in), optional :: unit

      integer :: n, info, iunit
      integer, allocatable :: ipiv(:)
      real(wp), allocatable :: amat_copy(:, :), q_mat(:, :)

      iunit = output_unit
      if (present(unit)) iunit = unit
      n = size(amat, 1)
      ! Copy matrix (LAPACK overwrites input)
      allocate (amat_copy(n, n))
      amat_copy = amat

      ! Reshape RHS into 2D matrix for getrs (n, 1)
      allocate (q_mat(n, 1))
      q_mat(:, 1) = rhs

      ! Allocate pivot indices
      allocate (ipiv(n))

      ! LU factorization
      call getrf(amat_copy, ipiv, info)
      if (info /= 0) then
         write (iunit, "(A,I0)") "[solve_pcm_lu] LAPACK getrf failed with info = ", info
         write (iunit, "(A,I0)") "[solve_pcm_lu] Matrix size n = ", n
         call fatal_error(error, "[solve_pcm_lu] LAPACK getrf failed")
         return
      end if

      ! Solve using factorization (getrs expects 2D matrix)
      call getrs(amat_copy, q_mat, ipiv, info)
      if (info /= 0) then
         call fatal_error(error, "[solve_pcm_lu] LAPACK getrs failed")
         return
      end if

      ! Extract solution from 2D matrix
      q = q_mat(:, 1)

   end subroutine solve_pcm_lu

   !> Solve PCM system using Cholesky factorization
   !> Solves A*q = rhs via LAPACK's Cholesky decomposition (DPOTRF + DPOTRS).
   !> Assumes A is symmetric positive definite - faster than LU for such matrices.
   subroutine solve_pcm_cholesky(amat, rhs, q, error)
      !> System matrix (ngrid, ngrid) - must be symmetric positive definite
      real(wp), intent(in) :: amat(:, :)
      !> Right-hand side vector (ngrid)
      real(wp), intent(in) :: rhs(:)
      !> Solution vector - surface charges (ngrid)
      real(wp), intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n, info
      real(wp), allocatable :: amat_copy(:, :), q_mat(:, :)

      n = size(amat, 1)

      ! Copy matrix (LAPACK overwrites input)
      allocate (amat_copy(n, n))
      amat_copy = amat

      ! Reshape RHS into 2D matrix for potrs (n, 1)
      allocate (q_mat(n, 1))
      q_mat(:, 1) = rhs

      ! Cholesky factorization: A = L*L^T (lower triangular)
      call potrf(amat_copy, info, uplo="l")
      if (info /= 0) then
         if (info > 0) then
            call fatal_error(error, "[solve_pcm_cholesky] Matrix not positive definite")
         else
            call fatal_error(error, "[solve_pcm_cholesky] LAPACK potrf failed")
         end if
         return
      end if

      ! Solve using factorization (potrs expects 2D matrix)
      call potrs(amat_copy, q_mat, info, uplo="l")
      if (info /= 0) then
         call fatal_error(error, "[solve_pcm_cholesky] LAPACK potrs failed")
         return
      end if

      ! Extract solution from 2D matrix
      q = q_mat(:, 1)

   end subroutine solve_pcm_cholesky

   !> Solve PCM system using matrix inversion
   !> Computes A^(-1) and then q = A^(-1)*rhs. Efficient if matrix is reused many times.
   subroutine solve_pcm_inversion(amat, rhs, q, error)
      !> System matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)
      !> Right-hand side vector (ngrid)
      real(wp), intent(in) :: rhs(:)
      !> Solution vector - surface charges (ngrid)
      real(wp), intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n, info
      integer, allocatable :: ipiv(:)
      real(wp), allocatable :: amat_inv(:, :)

      n = size(amat, 1)

      ! Copy matrix for inversion
      allocate (amat_inv(n, n))
      amat_inv = amat

      ! Allocate pivot indices
      allocate (ipiv(n))

      ! LU factorization
      call getrf(amat_inv, ipiv, info)
      if (info /= 0) then
         call fatal_error(error, "[solve_pcm_inversion] LAPACK getrf failed")
         return
      end if

      ! Compute inverse
      call getri(amat_inv, ipiv, info=info)
      if (info /= 0) then
         call fatal_error(error, "[solve_pcm_inversion] LAPACK getri failed")
         return
      end if

      ! Multiply: q = A^(-1) * rhs
      call gemv(amat_inv, rhs, q)

   end subroutine solve_pcm_inversion

   !> Solve the legacy matrix PCM system with restarted GMRES.
   !>
   !> Uses block-Jacobi preconditioning through pcm_gmres_solver.
   subroutine solve_pcm_iterative(amat, rhs, q, tol, maxiter, error)
      !> Dense system matrix
      real(wp), intent(in) :: amat(:, :)
      !> Right-hand side
      real(wp), intent(in) :: rhs(:)
      !> Solution
      real(wp), intent(out) :: q(:)
      !> Relative tolerance
      real(wp), intent(in) :: tol
      !> Maximum Krylov iterations
      integer, intent(in) :: maxiter
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(pcm_matrix_operator) :: op
      type(pcm_gmres_solver) :: solver
      real(wp), allocatable :: b(:, :), x(:, :)
      real(wp) :: rhs_norm
      integer :: n

      n = size(rhs)
      if (size(amat, 1) /= n .or. size(amat, 2) /= n .or. size(q) /= n) then
         call fatal_error(error, "[solve_pcm_iterative] inconsistent matrix and vector dimensions")
         return
      end if
      if (tol <= 0.0_wp .or. maxiter < 1) then
         call fatal_error(error, "[solve_pcm_iterative] tolerance and maxiter must be positive")
         return
      end if
      op%npt = n
      op%npoint = n
      op%id = 1
      op%a = amat
      solver%tol = tol
      rhs_norm = nrm2(rhs)
      if (rhs_norm > 0.0_wp) solver%tol = tol*max(1.0_wp, rhs_norm)/rhs_norm
      solver%max_iter = maxiter
      solver%restart = max(1, min(100, n))
      call solver%reset(1)
      allocate (b(n, 1), x(n, 1))
      b(:, 1) = rhs
      x = 0.0_wp
      call solver%solve_block(op, .false., b, x, error)
      if (allocated(error)) return
      q = x(:, 1)
   end subroutine solve_pcm_iterative

   !> Apply a dense matrix or its transpose.
   subroutine pcm_matrix_apply(self, trans, x, y)
      class(pcm_matrix_operator), intent(in) :: self
      logical, intent(in) :: trans
      real(wp), contiguous, intent(in) :: x(:)
      real(wp), contiguous, intent(out) :: y(:)

      if (trans) then
         y = matmul(transpose(self%a), x)
      else
         y = matmul(self%a, x)
      end if
   end subroutine pcm_matrix_apply

   !> Copy the dense matrix.
   subroutine pcm_matrix_assemble(self, a, error)
      class(pcm_matrix_operator), intent(in) :: self
      real(wp), contiguous, intent(out) :: a(:, :)
      type(error_type), allocatable, intent(out) :: error

      a = self%a
   end subroutine pcm_matrix_assemble

   !> Extract one principal submatrix for block-Jacobi.
   subroutine pcm_matrix_assemble_block(self, idx, a, error)
      class(pcm_matrix_operator), intent(in) :: self
      integer, contiguous, intent(in) :: idx(:)
      real(wp), contiguous, intent(out) :: a(:, :)
      type(error_type), allocatable, intent(out) :: error

      a = self%a(idx, idx)
   end subroutine pcm_matrix_assemble_block

   !> Default reset: the solver keeps no per-operator state
   !>
   !> @param[in,out] self         solver state
   !> @param[in]     n_operators  number of operators of the generated module
   subroutine pcm_solver_reset(self, n_operators)
      !> Solver state
      class(pcm_solver), intent(inout) :: self
      !> Number of operators
      integer, intent(in) :: n_operators

      associate (unused_self => self, unused_n => n_operators)
      end associate
   end subroutine pcm_solver_reset

   !> Default setup: nothing to prepare
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     operator
   !> @param[out]    error  never allocated
   subroutine pcm_solver_setup(self, op, error)
      !> Solver state
      class(pcm_solver), intent(inout) :: self
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      associate (unused_self => self, unused_op => op)
      end associate
      ! intent(out) already deallocated it; this marks the dummy as defined
      if (allocated(error)) deallocate (error)
   end subroutine pcm_solver_setup

   !> Solve for one right-hand side through solve_block
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     operator
   !> @param[in]     trans  solve with the transpose if true
   !> @param[in]     rhs    right-hand side, (npt)
   !> @param[in,out] x      initial guess on entry, solution on exit, (npt)
   !> @param[out]    error  allocated on failure
   subroutine pcm_solver_solve(self, op, trans, rhs, x, error)
      !> Solver state
      class(pcm_solver), intent(inout) :: self
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Solve with the transpose
      logical, intent(in) :: trans
      !> Right-hand side
      real(wp), contiguous, intent(in) :: rhs(:)
      !> Initial guess on entry, solution on exit
      real(wp), contiguous, intent(inout) :: x(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_solver_solve]"
      real(wp), allocatable :: rhs_block(:, :), x_block(:, :)
      integer :: stat

      if (size(x) /= size(rhs)) then
         call fatal_error(error, here//" solution of length "//int_str(size(x)) &
            & //" does not match right-hand side of length "//int_str(size(rhs)))
         return
      end if
      allocate (rhs_block(size(rhs), 1), stat=stat)
      if (stat == 0) allocate (x_block(size(rhs), 1), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, here//" failed to allocate the right-hand side block")
         return
      end if
      rhs_block(:, 1) = rhs
      x_block(:, 1) = x
      call self%solve_block(op, trans, rhs_block, x_block, error)
      if (allocated(error)) return
      x = x_block(:, 1)
   end subroutine pcm_solver_solve

   !> Default complex setup: nothing to prepare
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     complex operator
   !> @param[out]    error  never allocated
   subroutine pcm_solver_setup_complex(self, op, error)
      !> Solver state
      class(pcm_solver), intent(inout) :: self
      !> Complex operator
      class(pcm_complex_operator), intent(in) :: op
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      associate (unused_self => self, unused_op => op)
      end associate
      ! intent(out) already deallocated it; this marks the dummy as defined
      if (allocated(error)) deallocate (error)
   end subroutine pcm_solver_setup_complex

   !> Solve a real operator with a complex right-hand side through solve_block
   !>
   !> - Real and imaginary parts form a block of two real right-hand sides, so
   !>   the real factors or preconditioner of op are reused
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     real operator
   !> @param[in]     trans  solve with the transpose if true
   !> @param[in]     rhs    complex right-hand side, (npt)
   !> @param[in,out] x      complex initial guess on entry, solution on exit, (npt)
   !> @param[out]    error  allocated on failure
   subroutine pcm_solver_solve_complex_real(self, op, trans, rhs, x, error)
      !> Solver state
      class(pcm_solver), intent(inout) :: self
      !> Real operator
      class(pcm_operator), intent(in) :: op
      !> Solve with the transpose
      logical, intent(in) :: trans
      !> Complex right-hand side
      complex(wp), contiguous, intent(in) :: rhs(:)
      !> Complex initial guess on entry, solution on exit
      complex(wp), contiguous, intent(inout) :: x(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_solver_solve_complex]"
      real(wp), allocatable :: rhs_block(:, :), x_block(:, :)
      integer :: stat

      if (size(x) /= size(rhs)) then
         call fatal_error(error, here//" solution of length "//int_str(size(x)) &
            & //" does not match right-hand side of length "//int_str(size(rhs)))
         return
      end if
      allocate (rhs_block(size(rhs), 2), stat=stat)
      if (stat == 0) allocate (x_block(size(rhs), 2), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, here//" failed to allocate the right-hand side block")
         return
      end if
      rhs_block(:, 1) = real(rhs, wp)
      rhs_block(:, 2) = aimag(rhs)
      x_block(:, 1) = real(x, wp)
      x_block(:, 2) = aimag(x)
      call self%solve_block(op, trans, rhs_block, x_block, error)
      if (allocated(error)) return
      x = cmplx(x_block(:, 1), x_block(:, 2), kind=wp)
   end subroutine pcm_solver_solve_complex_real

   !> Default complex solve with a complex operator: not supported
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     complex operator
   !> @param[in]     trans  solve with the transpose if true
   !> @param[in]     rhs    complex right-hand side, (npt)
   !> @param[in,out] x      complex initial guess, unchanged
   !> @param[out]    error  always allocated
   subroutine pcm_solver_solve_complex_complex(self, op, trans, rhs, x, error)
      !> Solver state
      class(pcm_solver), intent(inout) :: self
      !> Complex operator
      class(pcm_complex_operator), intent(in) :: op
      !> Solve with the transpose
      logical, intent(in) :: trans
      !> Complex right-hand side
      complex(wp), contiguous, intent(in) :: rhs(:)
      !> Complex initial guess, unchanged
      complex(wp), contiguous, intent(inout) :: x(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      associate (unused_self => self, unused_op => op, unused_trans => trans, &
         & unused_rhs => rhs, unused_x => x)
      end associate
      call fatal_error(error, "[pcm_solver] complex operators are not supported by this solver")
   end subroutine pcm_solver_solve_complex_complex

   !> Solutions of the default solver are exact up to rounding (direct solvers)
   !>
   !> @param[in] self  solver state
   pure function pcm_solver_inexact(self) result(inexact)
      !> Solver state
      class(pcm_solver), intent(in) :: self
      !> Solutions satisfy the system only to a tolerance
      logical :: inexact

      associate (unused_self => self)
      end associate
      inexact = .false.
   end function pcm_solver_inexact

   !> By default the generated kernels store the operator matrices
   !>
   !> @param[in] self        solver state
   !> @param[in] npt         operator dimension
   !> @param[in] n_matrices  number of (npt, npt) matrices the kernels would store
   pure function pcm_solver_matrix_free(self, npt, n_matrices) result(matrix_free)
      !> Solver state
      class(pcm_solver), intent(in) :: self
      !> Operator dimension
      integer, intent(in) :: npt
      !> Number of matrices the kernels would store
      integer, intent(in) :: n_matrices
      !> Apply the operators without storing their matrices
      logical :: matrix_free

      associate (unused_self => self, unused_npt => npt, unused_n => n_matrices)
      end associate
      matrix_free = .false.
   end function pcm_solver_matrix_free

   !> Discard the factors of a previous phase
   !>
   !> @param[in,out] self         solver state
   !> @param[in]     n_operators  number of operators of the generated module
   subroutine pcm_dense_reset(self, n_operators)
      !> Solver state
      class(pcm_dense_solver), intent(inout) :: self
      !> Number of operators
      integer, intent(in) :: n_operators

      if (allocated(self%factors)) deallocate (self%factors)
      allocate (self%factors(max(0, n_operators)))
   end subroutine pcm_dense_reset

   !> Assemble and factor an operator, replacing any factors cached for its id
   !>
   !> - positive_definite: Cholesky; a failure is an error (no LU fallback)
   !> - symmetric: Bunch-Kaufman; otherwise LU
   !> - On failure the operator stays unfactored
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     operator
   !> @param[out]    error  allocated on failure
   subroutine pcm_dense_setup(self, op, error)
      !> Solver state
      class(pcm_dense_solver), intent(inout) :: self
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_dense_setup]"
      real(wp), allocatable :: a(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n_slots, kind, info, stat

      n_slots = 0
      if (allocated(self%factors)) n_slots = size(self%factors)
      call check_operator(op%id, op%npt, n_slots, here, error)
      if (allocated(error)) return
      self%factors(op%id) = pcm_dense_factors()

      allocate (a(op%npt, op%npt), stat=stat)
      if (stat == 0 .and. .not. op%positive_definite) allocate (ipiv(op%npt), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, here//" failed to allocate the matrix of operator " &
            & //int_str(op%id))
         return
      end if
      call op%assemble(a, error)
      if (allocated(error)) return

      if (op%positive_definite) then
         kind = pcm_factor_cholesky
         call potrf(a, info, uplo="l")
         if (info > 0) then
            call fatal_error(error, here//" operator "//int_str(op%id) &
               & //" is declared positive definite but is not (potrf info = " &
               & //int_str(info)//")")
            return
         end if
      else if (op%symmetric) then
         kind = pcm_factor_symmetric
         call sytrf(a, ipiv, info, uplo="l")
      else
         kind = pcm_factor_lu
         call getrf(a, ipiv, info)
      end if
      call check_factor_info(info, kind, op%id, here, error)
      if (allocated(error)) return

      associate (f => self%factors(op%id))
         call move_alloc(a, f%a)
         if (allocated(ipiv)) call move_alloc(ipiv, f%ipiv)
         f%kind = kind
         f%ready = .true.
      end associate
      self%n_factorizations = self%n_factorizations + 1
   end subroutine pcm_dense_setup

   !> Solve with the factors of an operator, factoring it on first use
   !>
   !> - Symmetric factors serve both orientations; LU uses getrs with trans = "t"
   !> - The initial guesses in x are ignored
   !> - n_solves grows by the number of right-hand sides
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     operator
   !> @param[in]     trans  solve with the transpose if true
   !> @param[in]     rhs    right-hand sides, (npt, nrhs)
   !> @param[in,out] x      ignored on entry, solutions on exit, (npt, nrhs)
   !> @param[out]    error  allocated on failure
   subroutine pcm_dense_solve_block(self, op, trans, rhs, x, error)
      !> Solver state
      class(pcm_dense_solver), intent(inout) :: self
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Solve with the transpose
      logical, intent(in) :: trans
      !> Right-hand sides
      real(wp), contiguous, intent(in) :: rhs(:, :)
      !> Ignored on entry, solutions on exit
      real(wp), contiguous, intent(inout) :: x(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_dense_solve_block]"
      integer :: n_slots, info

      call check_block_shapes(op, rhs, x, here, error)
      if (allocated(error)) return
      n_slots = 0
      if (allocated(self%factors)) n_slots = size(self%factors)
      call check_operator(op%id, op%npt, n_slots, here, error)
      if (allocated(error)) return
      if (.not. self%factors(op%id)%ready) then
         call self%setup(op, error)
         if (allocated(error)) return
      end if

      x = rhs
      associate (f => self%factors(op%id))
         select case (f%kind)
         case (pcm_factor_cholesky)
            call potrs(f%a, x, info, uplo="l")
         case (pcm_factor_symmetric)
            call sytrs(f%a, x, f%ipiv, info, uplo="l")
         case (pcm_factor_lu)
            call getrs(f%a, x, f%ipiv, info, trans=merge("t", "n", trans))
         case default
            call fatal_error(error, here//" operator id "//int_str(op%id) &
               & //" holds complex factors; set up the real operator first")
            return
         end select
         if (info /= 0) then
            call fatal_error(error, here//" "//factor_name(f%kind) &
               & //" solve failed for operator "//int_str(op%id)//" (info = " &
               & //int_str(info)//")")
            return
         end if
      end associate
      self%n_solves = self%n_solves + size(rhs, 2)
   end subroutine pcm_dense_solve_block

   !> Assemble and factor a complex operator, replacing any factors cached for its id
   !>
   !> - symmetric (A = A**T): complex Bunch-Kaufman (zsytrf, lower); otherwise
   !>   complex LU (zgetrf)
   !> - A singular or rejected factorization is an error and the operator stays
   !>   unfactored
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     complex operator
   !> @param[out]    error  allocated on failure
   subroutine pcm_dense_setup_complex(self, op, error)
      !> Solver state
      class(pcm_dense_solver), intent(inout) :: self
      !> Complex operator
      class(pcm_complex_operator), intent(in) :: op
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_dense_setup_complex]"
      complex(wp), allocatable :: a(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n_slots, kind, info, stat

      n_slots = 0
      if (allocated(self%factors)) n_slots = size(self%factors)
      call check_operator(op%id, op%npt, n_slots, here, error)
      if (allocated(error)) return
      self%factors(op%id) = pcm_dense_factors()

      allocate (a(op%npt, op%npt), stat=stat)
      if (stat == 0) allocate (ipiv(op%npt), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, here//" failed to allocate the matrix of operator " &
            & //int_str(op%id))
         return
      end if
      call op%assemble(a, error)
      if (allocated(error)) return

      if (op%symmetric) then
         kind = pcm_factor_complex_symmetric
         call sytrf(a, ipiv, info, uplo="l")
      else
         kind = pcm_factor_complex_lu
         call getrf(a, ipiv, info)
      end if
      call check_factor_info(info, kind, op%id, here, error)
      if (allocated(error)) return

      associate (f => self%factors(op%id))
         call move_alloc(a, f%z)
         call move_alloc(ipiv, f%ipiv)
         f%kind = kind
         f%ready = .true.
      end associate
      self%n_factorizations = self%n_factorizations + 1
   end subroutine pcm_dense_setup_complex

   !> Solve with the complex factors of an operator, factoring it on first use
   !>
   !> - Symmetric factors serve both orientations; LU uses zgetrs with trans = "t"
   !>   (plain transpose, not conjugated)
   !> - The initial guess in x is ignored
   !> - n_solves grows by one
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     complex operator
   !> @param[in]     trans  solve with the transpose if true
   !> @param[in]     rhs    complex right-hand side, (npt)
   !> @param[in,out] x      ignored on entry, complex solution on exit, (npt)
   !> @param[out]    error  allocated on failure
   subroutine pcm_dense_solve_complex_complex(self, op, trans, rhs, x, error)
      !> Solver state
      class(pcm_dense_solver), intent(inout) :: self
      !> Complex operator
      class(pcm_complex_operator), intent(in) :: op
      !> Solve with the transpose
      logical, intent(in) :: trans
      !> Complex right-hand side
      complex(wp), contiguous, intent(in) :: rhs(:)
      !> Ignored on entry, complex solution on exit
      complex(wp), contiguous, intent(inout) :: x(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_dense_solve_complex]"
      complex(wp), allocatable :: x_block(:, :)
      integer :: n_slots, info, stat

      if (size(rhs) /= op%npt .or. size(x) /= size(rhs)) then
         call fatal_error(error, here//" right-hand side of length "//int_str(size(rhs)) &
            & //" and solution of length "//int_str(size(x))//" do not match operator " &
            & //int_str(op%id)//" of dimension "//int_str(op%npt))
         return
      end if
      n_slots = 0
      if (allocated(self%factors)) n_slots = size(self%factors)
      call check_operator(op%id, op%npt, n_slots, here, error)
      if (allocated(error)) return
      if (.not. self%factors(op%id)%ready) then
         call self%setup_complex(op, error)
         if (allocated(error)) return
      end if
      allocate (x_block(op%npt, 1), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, here//" failed to allocate the right-hand side block")
         return
      end if

      x_block(:, 1) = rhs
      associate (f => self%factors(op%id))
         select case (f%kind)
         case (pcm_factor_complex_symmetric)
            call sytrs(f%z, x_block, f%ipiv, info, uplo="l")
         case (pcm_factor_complex_lu)
            call getrs(f%z, x_block, f%ipiv, info, trans=merge("t", "n", trans))
         case default
            call fatal_error(error, here//" operator id "//int_str(op%id) &
               & //" holds real factors; set up the complex operator first")
            return
         end select
         if (info /= 0) then
            call fatal_error(error, here//" "//factor_name(f%kind) &
               & //" solve failed for operator "//int_str(op%id)//" (info = " &
               & //int_str(info)//")")
            return
         end if
      end associate
      x = x_block(:, 1)
      self%n_solves = self%n_solves + 1
   end subroutine pcm_dense_solve_complex_complex

   !> Factorization used for an operator
   !>
   !> @param[in] self  solver state
   !> @param[in] id    operator id
   function pcm_dense_factorization(self, id) result(kind)
      !> Solver state
      class(pcm_dense_solver), intent(in) :: self
      !> Operator id
      integer, intent(in) :: id
      !> One of the pcm_factor_* kinds; 0 if not factored
      integer :: kind

      kind = 0
      if (.not. allocated(self%factors)) return
      if (id < 1 .or. id > size(self%factors)) return
      if (self%factors(id)%ready) kind = self%factors(id)%kind
   end function pcm_dense_factorization

   !> Default reset: the preconditioner keeps no per-operator state
   !>
   !> @param[in,out] self         preconditioner state
   !> @param[in]     n_operators  number of operators of the generated module
   subroutine pcm_preconditioner_reset(self, n_operators)
      !> Preconditioner state
      class(pcm_preconditioner), intent(inout) :: self
      !> Number of operators
      integer, intent(in) :: n_operators

      associate (unused_self => self, unused_n => n_operators)
      end associate
   end subroutine pcm_preconditioner_reset

   !> Discard the block factors of a previous phase
   !>
   !> @param[in,out] self         preconditioner state
   !> @param[in]     n_operators  number of operators of the generated module
   subroutine pcm_block_jacobi_reset(self, n_operators)
      !> Preconditioner state
      class(pcm_block_jacobi), intent(inout) :: self
      !> Number of operators
      integer, intent(in) :: n_operators

      if (allocated(self%slots)) deallocate (self%slots)
      allocate (self%slots(max(0, n_operators)))
   end subroutine pcm_block_jacobi_reset

   !> Cluster the points of an operator and factor its diagonal blocks
   !>
   !> - Replaces the blocks cached for the operator id; on failure the id
   !>   stays empty
   !>
   !> @param[in,out] self   preconditioner state
   !> @param[in]     op     operator
   !> @param[out]    error  allocated on failure
   subroutine pcm_block_jacobi_setup(self, op, error)
      !> Preconditioner state
      class(pcm_block_jacobi), intent(inout) :: self
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_block_jacobi_setup]"
      type(pcm_jacobi_block), allocatable :: blocks(:)
      integer :: n_slots, m, ib, nb, info, stat

      n_slots = 0
      if (allocated(self%slots)) n_slots = size(self%slots)
      call check_operator(op%id, op%npt, n_slots, here, error)
      if (allocated(error)) return
      self%slots(op%id) = pcm_jacobi_slot()
      call jacobi_layout(self%block_size, op%id, op%npt, op%npoint, op%points, here, blocks, m, &
         & error)
      if (allocated(error)) return

      do ib = 1, size(blocks)
         associate (b => blocks(ib))
            nb = size(b%unknowns)/m
            allocate (b%a(m*nb, m*nb), stat=stat)
            if (stat /= 0) then
               call fatal_error(error, here//" failed to allocate block "//int_str(ib) &
                  & //" of operator "//int_str(op%id))
               return
            end if
            call op%assemble_block(b%unknowns(1:nb), b%a, error)
            if (allocated(error)) return
            call getrf(b%a, b%ipiv, info)
            call check_block_info(info, ib, op%id, "getrf", here, error)
            if (allocated(error)) return
         end associate
      end do

      associate (slot => self%slots(op%id))
         call move_alloc(blocks, slot%blocks)
         slot%npt = op%npt
         slot%kind = slot_real
      end associate
   end subroutine pcm_block_jacobi_setup

   !> Cluster the points of a complex operator and factor its diagonal blocks
   !>
   !> - Replaces the blocks cached for the operator id; on failure the id
   !>   stays empty
   !>
   !> @param[in,out] self   preconditioner state
   !> @param[in]     op     complex operator
   !> @param[out]    error  allocated on failure
   subroutine pcm_block_jacobi_setup_complex(self, op, error)
      !> Preconditioner state
      class(pcm_block_jacobi), intent(inout) :: self
      !> Complex operator
      class(pcm_complex_operator), intent(in) :: op
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_block_jacobi_setup_complex]"
      type(pcm_jacobi_block), allocatable :: blocks(:)
      integer :: n_slots, m, ib, nb, info, stat

      n_slots = 0
      if (allocated(self%slots)) n_slots = size(self%slots)
      call check_operator(op%id, op%npt, n_slots, here, error)
      if (allocated(error)) return
      self%slots(op%id) = pcm_jacobi_slot()
      call jacobi_layout(self%block_size, op%id, op%npt, op%npoint, op%points, here, blocks, m, &
         & error)
      if (allocated(error)) return

      do ib = 1, size(blocks)
         associate (b => blocks(ib))
            nb = size(b%unknowns)/m
            allocate (b%z(m*nb, m*nb), stat=stat)
            if (stat /= 0) then
               call fatal_error(error, here//" failed to allocate block "//int_str(ib) &
                  & //" of operator "//int_str(op%id))
               return
            end if
            call op%assemble_block(b%unknowns(1:nb), b%z, error)
            if (allocated(error)) return
            call getrf(b%z, b%ipiv, info)
            call check_block_info(info, ib, op%id, "zgetrf", here, error)
            if (allocated(error)) return
         end associate
      end do

      associate (slot => self%slots(op%id))
         call move_alloc(blocks, slot%blocks)
         slot%npt = op%npt
         slot%kind = slot_complex
      end associate
   end subroutine pcm_block_jacobi_setup_complex

   !> Apply the block-Jacobi preconditioner: y = M^-1 x or M^-T x
   !>
   !> @param[in]  self   preconditioner state
   !> @param[in]  op     operator the preconditioner was set up for
   !> @param[in]  trans  apply M^-T if true
   !> @param[in]  x      vector of length npt
   !> @param[out] y      vector of length npt
   !> @param[out] error  allocated if the operator is not set up
   subroutine pcm_block_jacobi_apply(self, op, trans, x, y, error)
      !> Preconditioner state
      class(pcm_block_jacobi), intent(in) :: self
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Apply the transpose
      logical, intent(in) :: trans
      !> Input vector
      real(wp), contiguous, intent(in) :: x(:)
      !> Result vector
      real(wp), contiguous, intent(out) :: y(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_block_jacobi_apply]"
      real(wp), allocatable :: w(:, :)
      integer :: ib, nu, info, stat

      call check_jacobi_slot(self, op%id, op%npt, slot_real, size(x), size(y), here, error)
      if (allocated(error)) return
      associate (blocks => self%slots(op%id)%blocks)
         nu = 0
         do ib = 1, size(blocks)
            nu = max(nu, size(blocks(ib)%unknowns))
         end do
         allocate (w(nu, 1), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, here//" failed to allocate the block work array")
            return
         end if
         do ib = 1, size(blocks)
            associate (b => blocks(ib))
               nu = size(b%unknowns)
               w(1:nu, 1) = x(b%unknowns)
               call getrs(b%a, w(1:nu, 1:1), b%ipiv, info, trans=merge("t", "n", trans))
               if (info /= 0) then
                  call fatal_error(error, here//" getrs failed for block "//int_str(ib) &
                     & //" of operator "//int_str(op%id)//" (info = "//int_str(info)//")")
                  return
               end if
               y(b%unknowns) = w(1:nu, 1)
            end associate
         end do
      end associate
   end subroutine pcm_block_jacobi_apply

   !> Apply the complex block-Jacobi preconditioner: y = M^-1 x or M^-T x
   !>
   !> @param[in]  self   preconditioner state
   !> @param[in]  op     complex operator the preconditioner was set up for
   !> @param[in]  trans  apply M^-T (not conjugated) if true
   !> @param[in]  x      vector of length npt
   !> @param[out] y      vector of length npt
   !> @param[out] error  allocated if the operator is not set up
   subroutine pcm_block_jacobi_apply_complex(self, op, trans, x, y, error)
      !> Preconditioner state
      class(pcm_block_jacobi), intent(in) :: self
      !> Complex operator
      class(pcm_complex_operator), intent(in) :: op
      !> Apply the transpose
      logical, intent(in) :: trans
      !> Input vector
      complex(wp), contiguous, intent(in) :: x(:)
      !> Result vector
      complex(wp), contiguous, intent(out) :: y(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_block_jacobi_apply_complex]"
      complex(wp), allocatable :: w(:, :)
      integer :: ib, nu, info, stat

      call check_jacobi_slot(self, op%id, op%npt, slot_complex, size(x), size(y), here, error)
      if (allocated(error)) return
      associate (blocks => self%slots(op%id)%blocks)
         nu = 0
         do ib = 1, size(blocks)
            nu = max(nu, size(blocks(ib)%unknowns))
         end do
         allocate (w(nu, 1), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, here//" failed to allocate the block work array")
            return
         end if
         do ib = 1, size(blocks)
            associate (b => blocks(ib))
               nu = size(b%unknowns)
               w(1:nu, 1) = x(b%unknowns)
               call getrs(b%z, w(1:nu, 1:1), b%ipiv, info, trans=merge("t", "n", trans))
               if (info /= 0) then
                  call fatal_error(error, here//" zgetrs failed for block "//int_str(ib) &
                     & //" of operator "//int_str(op%id)//" (info = "//int_str(info)//")")
                  return
               end if
               y(b%unknowns) = w(1:nu, 1)
            end associate
         end do
      end associate
   end subroutine pcm_block_jacobi_apply_complex

   !> Number of diagonal blocks of an operator id
   !>
   !> @param[in] self  preconditioner state
   !> @param[in] id    operator id
   function pcm_block_jacobi_n_blocks(self, id) result(n)
      !> Preconditioner state
      class(pcm_block_jacobi), intent(in) :: self
      !> Operator id
      integer, intent(in) :: id
      !> Number of blocks; 0 if the id is not set up
      integer :: n

      n = 0
      if (.not. allocated(self%slots)) return
      if (id < 1 .or. id > size(self%slots)) return
      if (self%slots(id)%kind == slot_empty) return
      n = size(self%slots(id)%blocks)
   end function pcm_block_jacobi_n_blocks

   !> Check that an operator id holds block factors of the expected kind
   !>
   !> @param[in]  self     preconditioner state
   !> @param[in]  id       operator id
   !> @param[in]  npt      operator dimension
   !> @param[in]  kind     expected kind, slot_real or slot_complex
   !> @param[in]  nx       length of the input vector
   !> @param[in]  ny       length of the result vector
   !> @param[in]  routine  "[routine]" prefix of error messages
   !> @param[out] error    allocated if the id is not set up for this operator
   subroutine check_jacobi_slot(self, id, npt, kind, nx, ny, routine, error)
      !> Preconditioner state
      class(pcm_block_jacobi), intent(in) :: self
      !> Operator id
      integer, intent(in) :: id
      !> Operator dimension
      integer, intent(in) :: npt
      !> Expected kind of the factors
      integer, intent(in) :: kind
      !> Length of the input vector
      integer, intent(in) :: nx
      !> Length of the result vector
      integer, intent(in) :: ny
      !> Prefix of error messages
      character(len=*), intent(in) :: routine
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n_slots

      n_slots = 0
      if (allocated(self%slots)) n_slots = size(self%slots)
      call check_operator(id, npt, n_slots, routine, error)
      if (allocated(error)) return
      if (self%slots(id)%kind /= kind .or. self%slots(id)%npt /= npt) then
         call fatal_error(error, routine//" operator "//int_str(id) &
            & //" is not set up for this preconditioner")
         return
      end if
      if (nx /= npt .or. ny /= npt) then
         call fatal_error(error, routine//" vectors of length "//int_str(nx)//" and " &
            & //int_str(ny)//" do not match operator "//int_str(id)//" of dimension " &
            & //int_str(npt))
         return
      end if
   end subroutine check_jacobi_slot

   !> Turn the info of a block factorization into an error
   !>
   !> @param[in]  info     info returned by getrf or zgetrf
   !> @param[in]  ib       block index
   !> @param[in]  id       operator id
   !> @param[in]  name     name of the factorization routine
   !> @param[in]  routine  "[routine]" prefix of error messages
   !> @param[out] error    allocated if info is nonzero
   subroutine check_block_info(info, ib, id, name, routine, error)
      !> Factorization status
      integer, intent(in) :: info
      !> Block index
      integer, intent(in) :: ib
      !> Operator id
      integer, intent(in) :: id
      !> Name of the factorization routine
      character(len=*), intent(in) :: name
      !> Prefix of error messages
      character(len=*), intent(in) :: routine
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (info > 0) then
         call fatal_error(error, routine//" diagonal block "//int_str(ib)//" of operator " &
            & //int_str(id)//" is singular ("//name//" info = "//int_str(info)//")")
         return
      else if (info < 0) then
         call fatal_error(error, routine//" "//name//" rejected argument "//int_str(-info) &
            & //" for block "//int_str(ib)//" of operator "//int_str(id))
         return
      end if
   end subroutine check_block_info

   !> Point clusters and unknowns of the diagonal blocks of an operator
   !>
   !> - npoint 0 means npt (a single field); npt must be a multiple of npoint
   !> - Clusters of at most block_size points: recursive coordinate bisection
   !>   of points if associated, else contiguous index ranges
   !> - Block unknowns (a - 1) npoint + idx(k) in row order a = 1..m, k = 1..nb;
   !>   the first nb of them are the point indices idx
   !>
   !> @param[in]  block_size  maximum number of points per cluster
   !> @param[in]  id          operator id
   !> @param[in]  npt         operator dimension
   !> @param[in]  npoint      number of surface points of the operator, 0 for npt
   !> @param[in]  points      positions (3, npoint), may be null
   !> @param[in]  routine     "[routine]" prefix of error messages
   !> @param[out] blocks      blocks with unknowns and pivot arrays allocated
   !> @param[out] m           number of stacked fields, npt / npoint
   !> @param[out] error       allocated on invalid input
   subroutine jacobi_layout(block_size, id, npt, npoint, points, routine, blocks, m, error)
      !> Maximum number of points per cluster
      integer, intent(in) :: block_size
      !> Operator id
      integer, intent(in) :: id
      !> Operator dimension
      integer, intent(in) :: npt
      !> Number of surface points of the operator
      integer, intent(in) :: npoint
      !> Positions of the surface points
      real(wp), pointer, contiguous, intent(in) :: points(:, :)
      !> Prefix of error messages
      character(len=*), intent(in) :: routine
      !> Blocks with unknowns and pivot arrays
      type(pcm_jacobi_block), allocatable, intent(out) :: blocks(:)
      !> Number of stacked fields
      integer, intent(out) :: m
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: perm(:), first(:)
      integer :: np, n_clusters, ib, ic, a, nb, stat

      m = 1
      if (block_size < 1) then
         call fatal_error(error, routine//" block_size must be positive, got " &
            & //int_str(block_size))
         return
      end if
      np = npoint
      if (np == 0) np = npt
      if (npt > 0 .and. (np < 1 .or. mod(npt, max(np, 1)) /= 0)) then
         call fatal_error(error, routine//" operator "//int_str(id)//" of dimension " &
            & //int_str(npt)//" is not a multiple of its "//int_str(npoint)//" points")
         return
      end if
      if (npt > 0) m = npt/np
      if (associated(points)) then
         if (size(points, 1) /= 3 .or. size(points, 2) /= np) then
            call fatal_error(error, routine//" points of operator "//int_str(id) &
               & //" have shape ("//int_str(size(points, 1))//", " &
               & //int_str(size(points, 2))//"), expected (3, "//int_str(np)//")")
            return
         end if
      end if

      n_clusters = 0
      if (npt > 0) n_clusters = (np + block_size - 1)/block_size
      allocate (perm(max(np, 0)), stat=stat)
      if (stat == 0) allocate (first(n_clusters + 1), stat=stat)
      if (stat == 0) allocate (blocks(n_clusters), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, routine//" failed to allocate the clusters of operator " &
            & //int_str(id))
         return
      end if
      perm = [(ib, ib=1, size(perm))]
      if (n_clusters > 0) then
         if (associated(points)) then
            ic = 0
            call bisect_clusters(points, perm, 1, np, n_clusters, first, ic)
         else
            do ic = 1, n_clusters
               first(ic) = 1 + int((int(ic - 1, int64)*np)/n_clusters)
            end do
         end if
      end if
      first(n_clusters + 1) = np + 1

      do ib = 1, n_clusters
         nb = first(ib + 1) - first(ib)
         allocate (blocks(ib)%unknowns(m*nb), stat=stat)
         if (stat == 0) allocate (blocks(ib)%ipiv(m*nb), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, routine//" failed to allocate the clusters of operator " &
               & //int_str(id))
            return
         end if
         associate (idx => perm(first(ib):first(ib + 1) - 1))
            do a = 1, m
               blocks(ib)%unknowns((a - 1)*nb + 1:a*nb) = (a - 1)*np + idx
            end do
         end associate
      end do
   end subroutine jacobi_layout

   !> Recursive coordinate bisection of perm(lo:hi) into n_clusters clusters
   !>
   !> - Splits across the longest extent of the bounding box into parts
   !>   proportional to the cluster counts floor(n_clusters/2) and the rest, so
   !>   every cluster holds at most ceiling(size / n_clusters) points
   !> - Clusters are appended in depth-first order: first(ic) is the start of
   !>   cluster ic in perm
   !>
   !> @param[in]     points      positions, (3, npoint)
   !> @param[in,out] perm        point indices, reordered within lo..hi
   !> @param[in]     lo          first position of the range in perm
   !> @param[in]     hi          last position of the range in perm
   !> @param[in]     n_clusters  number of clusters of the range, >= 1
   !> @param[in,out] first       start positions of the clusters
   !> @param[in,out] ic          number of clusters found so far
   recursive subroutine bisect_clusters(points, perm, lo, hi, n_clusters, first, ic)
      !> Positions of the surface points
      real(wp), intent(in) :: points(:, :)
      !> Point indices
      integer, intent(inout) :: perm(:)
      !> First position of the range
      integer, intent(in) :: lo
      !> Last position of the range
      integer, intent(in) :: hi
      !> Number of clusters of the range
      integer, intent(in) :: n_clusters
      !> Start positions of the clusters
      integer, intent(inout) :: first(:)
      !> Number of clusters found so far
      integer, intent(inout) :: ic

      real(wp) :: extent(3)
      integer :: n_left, n_left_clusters, k, axis

      if (n_clusters <= 1) then
         ic = ic + 1
         first(ic) = lo
         return
      end if
      n_left_clusters = n_clusters/2
      n_left = int((int(hi - lo + 1, int64)*n_left_clusters)/n_clusters)
      do k = 1, 3
         extent(k) = maxval(points(k, perm(lo:hi))) - minval(points(k, perm(lo:hi)))
      end do
      axis = 1
      if (extent(2) > extent(axis)) axis = 2
      if (extent(3) > extent(axis)) axis = 3
      call select_smallest(points(axis, :), perm(lo:hi), n_left)
      call bisect_clusters(points, perm, lo, lo + n_left - 1, n_left_clusters, first, ic)
      call bisect_clusters(points, perm, lo + n_left, hi, n_clusters - n_left_clusters, &
         & first, ic)
   end subroutine bisect_clusters

   !> Reorder perm so that its first k entries have the k smallest keys
   !>
   !> - Selection by partitioning (Wirth); perm(k) has the k-th smallest key,
   !>   keys of perm(:k) are <= those of perm(k+1:)
   !>
   !> @param[in]     key   key per point index
   !> @param[in,out] perm  point indices
   !> @param[in]     k     number of smallest keys to gather, 1..size(perm)
   pure subroutine select_smallest(key, perm, k)
      !> Key per point index
      real(wp), intent(in) :: key(:)
      !> Point indices
      integer, intent(inout) :: perm(:)
      !> Number of smallest keys to gather
      integer, intent(in) :: k

      real(wp) :: pivot
      integer :: lo, hi, i, j, tmp

      lo = 1
      hi = size(perm)
      do while (lo < hi)
         pivot = key(perm(k))
         i = lo
         j = hi
         do while (i <= j)
            do while (key(perm(i)) < pivot)
               i = i + 1
            end do
            do while (pivot < key(perm(j)))
               j = j - 1
            end do
            if (i <= j) then
               tmp = perm(i)
               perm(i) = perm(j)
               perm(j) = tmp
               i = i + 1
               j = j - 1
            end if
         end do
         if (j < k) lo = i
         if (k < i) hi = j
      end do
   end subroutine select_smallest

   !> Discard the preconditioners of a previous phase
   !>
   !> @param[in,out] self         solver state
   !> @param[in]     n_operators  number of operators of the generated module
   subroutine pcm_gmres_reset(self, n_operators)
      !> Solver state
      class(pcm_gmres_solver), intent(inout) :: self
      !> Number of operators
      integer, intent(in) :: n_operators

      self%n_operators = max(0, n_operators)
      if (allocated(self%prepared)) deallocate (self%prepared)
      allocate (self%prepared(self%n_operators), source=slot_empty)
      if (allocated(self%preconditioner)) call self%preconditioner%reset(self%n_operators)
   end subroutine pcm_gmres_reset

   !> Build the preconditioner of an operator, replacing the one of its id
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     operator
   !> @param[out]    error  allocated on failure
   subroutine pcm_gmres_setup(self, op, error)
      !> Solver state
      class(pcm_gmres_solver), intent(inout) :: self
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call gmres_prepare(self, op%id, op%npt, "[pcm_gmres_setup]", error)
      if (allocated(error)) return
      self%prepared(op%id) = slot_empty
      call self%preconditioner%setup(op, error)
      if (allocated(error)) return
      self%prepared(op%id) = slot_real
   end subroutine pcm_gmres_setup

   !> Build the preconditioner of a complex operator, replacing the one of its id
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     complex operator
   !> @param[out]    error  allocated on failure
   subroutine pcm_gmres_setup_complex(self, op, error)
      !> Solver state
      class(pcm_gmres_solver), intent(inout) :: self
      !> Complex operator
      class(pcm_complex_operator), intent(in) :: op
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call gmres_prepare(self, op%id, op%npt, "[pcm_gmres_setup_complex]", error)
      if (allocated(error)) return
      self%prepared(op%id) = slot_empty
      call self%preconditioner%setup_complex(op, error)
      if (allocated(error)) return
      self%prepared(op%id) = slot_complex
   end subroutine pcm_gmres_setup_complex

   !> Solve each right-hand side by preconditioned GMRES, starting from x
   !>
   !> - Sets up the preconditioner first if the id has none
   !> - n_solves grows by the number of right-hand sides
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     operator
   !> @param[in]     trans  solve with the transpose if true
   !> @param[in]     rhs    right-hand sides, (npt, nrhs)
   !> @param[in,out] x      initial guesses on entry, solutions on exit, (npt, nrhs)
   !> @param[out]    error  allocated on failure
   subroutine pcm_gmres_solve_block(self, op, trans, rhs, x, error)
      !> Solver state
      class(pcm_gmres_solver), intent(inout) :: self
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Solve with the transpose
      logical, intent(in) :: trans
      !> Right-hand sides
      real(wp), contiguous, intent(in) :: rhs(:, :)
      !> Initial guesses on entry, solutions on exit
      real(wp), contiguous, intent(inout) :: x(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_gmres_solve_block]"
      integer :: k

      call check_block_shapes(op, rhs, x, here, error)
      if (allocated(error)) return
      call gmres_prepare(self, op%id, op%npt, here, error)
      if (allocated(error)) return
      if (self%prepared(op%id) /= slot_real) then
         call self%setup(op, error)
         if (allocated(error)) return
      end if
      do k = 1, size(rhs, 2)
         call gmres_real(self, op, trans, rhs(:, k), x(:, k), &
            & here//" operator "//int_str(op%id)//", right-hand side "//int_str(k)//":", error)
         if (allocated(error)) return
      end do
   end subroutine pcm_gmres_solve_block

   !> Solve a complex operator by preconditioned complex GMRES, starting from x
   !>
   !> - Sets up the preconditioner first if the id has none
   !> - n_solves grows by one
   !>
   !> @param[in,out] self   solver state
   !> @param[in]     op     complex operator
   !> @param[in]     trans  solve with the transpose (not conjugated) if true
   !> @param[in]     rhs    complex right-hand side, (npt)
   !> @param[in,out] x      complex initial guess on entry, solution on exit, (npt)
   !> @param[out]    error  allocated on failure
   subroutine pcm_gmres_solve_complex_complex(self, op, trans, rhs, x, error)
      !> Solver state
      class(pcm_gmres_solver), intent(inout) :: self
      !> Complex operator
      class(pcm_complex_operator), intent(in) :: op
      !> Solve with the transpose
      logical, intent(in) :: trans
      !> Complex right-hand side
      complex(wp), contiguous, intent(in) :: rhs(:)
      !> Complex initial guess on entry, solution on exit
      complex(wp), contiguous, intent(inout) :: x(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: here = "[pcm_gmres_solve_complex]"

      if (size(rhs) /= op%npt .or. size(x) /= size(rhs)) then
         call fatal_error(error, here//" right-hand side of length "//int_str(size(rhs)) &
            & //" and solution of length "//int_str(size(x))//" do not match operator " &
            & //int_str(op%id)//" of dimension "//int_str(op%npt))
         return
      end if
      call gmres_prepare(self, op%id, op%npt, here, error)
      if (allocated(error)) return
      if (self%prepared(op%id) /= slot_complex) then
         call self%setup_complex(op, error)
         if (allocated(error)) return
      end if
      call gmres_complex(self, op, trans, rhs, x, here//" operator "//int_str(op%id)//":", &
         & error)
   end subroutine pcm_gmres_solve_complex_complex

   !> GMRES solutions satisfy the system to the tolerance only
   !>
   !> @param[in] self  solver state
   pure function pcm_gmres_inexact(self) result(inexact)
      !> Solver state
      class(pcm_gmres_solver), intent(in) :: self
      !> Solutions satisfy the system only to a tolerance
      logical :: inexact

      associate (unused_self => self)
      end associate
      inexact = .true.
   end function pcm_gmres_inexact

   !> Matrix-free kernels by the storage setting
   !>
   !> @param[in] self        solver state
   !> @param[in] npt         operator dimension
   !> @param[in] n_matrices  number of (npt, npt) matrices the kernels would store
   pure function pcm_gmres_matrix_free(self, npt, n_matrices) result(matrix_free)
      !> Solver state
      class(pcm_gmres_solver), intent(in) :: self
      !> Operator dimension
      integer, intent(in) :: npt
      !> Number of matrices the kernels would store
      integer, intent(in) :: n_matrices
      !> Apply the operators without storing their matrices
      logical :: matrix_free

      associate (unused_npt => npt, unused_n_matrices => n_matrices)
      end associate
      matrix_free = self%storage == pcm_storage_matrix_free
   end function pcm_gmres_matrix_free

   !> Check the settings and the operator id; create the default preconditioner
   !>
   !> @param[in,out] self     solver state
   !> @param[in]     id       operator id
   !> @param[in]     npt      operator dimension
   !> @param[in]     routine  "[routine]" prefix of error messages
   !> @param[out]    error    allocated on invalid settings or id
   subroutine gmres_prepare(self, id, npt, routine, error)
      !> Solver state
      class(pcm_gmres_solver), intent(inout) :: self
      !> Operator id
      integer, intent(in) :: id
      !> Operator dimension
      integer, intent(in) :: npt
      !> Prefix of error messages
      character(len=*), intent(in) :: routine
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n_slots

      if (.not. (self%tol > 0.0_wp) .or. self%restart < 1 .or. self%max_iter < 0) then
         call fatal_error(error, routine//" invalid settings: tol = "//real_str(self%tol) &
            & //", restart = "//int_str(self%restart)//", max_iter = " &
            & //int_str(self%max_iter)//" (need tol > 0, restart >= 1, max_iter >= 0)")
         return
      end if
      n_slots = 0
      if (allocated(self%prepared)) n_slots = size(self%prepared)
      call check_operator(id, npt, n_slots, routine, error)
      if (allocated(error)) return
      if (.not. allocated(self%preconditioner)) then
         allocate (pcm_block_jacobi :: self%preconditioner)
         call self%preconditioner%reset(self%n_operators)
      end if
   end subroutine gmres_prepare

   !> Restarted right-preconditioned GMRES for one real right-hand side
   !>
   !> - Converged when norm(b - op x) <= tol norm(b), checked with the true
   !>   residual at the start and after every cycle
   !> - Updates the statistics of self
   !>
   !> @param[in,out] self     solver state with a set-up preconditioner
   !> @param[in]     op       operator
   !> @param[in]     trans    solve with the transpose if true
   !> @param[in]     b        right-hand side, (npt)
   !> @param[in,out] x        initial guess on entry, solution on exit, (npt)
   !> @param[in]     context  prefix of error messages
   !> @param[out]    error    allocated on non-convergence or breakdown
   subroutine gmres_real(self, op, trans, b, x, context, error)
      !> Solver state
      class(pcm_gmres_solver), intent(inout) :: self
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Solve with the transpose
      logical, intent(in) :: trans
      !> Right-hand side
      real(wp), contiguous, intent(in) :: b(:)
      !> Initial guess on entry, solution on exit
      real(wp), contiguous, intent(inout) :: x(:)
      !> Prefix of error messages
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: v(:, :), h(:, :), g(:), cs(:), sn(:), y(:), r(:), w(:), z(:)
      real(wp) :: b_norm, beta, threshold, h_next, tmp
      integer :: n, kdim, its, i, j, k, stat

      n = size(b)
      self%n_solves = self%n_solves + 1
      self%last_iterations = 0
      self%last_residual = 0.0_wp
      b_norm = nrm2(b)
      if (.not. ieee_is_finite(b_norm)) then
         call fatal_error(error, context//" right-hand side is not finite")
         return
      end if
      if (b_norm <= 0.0_wp) then
         x = 0.0_wp
         return
      end if

      kdim = max(1, min(self%restart, n))
      allocate (v(n, kdim + 1), stat=stat)
      if (stat == 0) allocate (h(kdim + 1, kdim), stat=stat)
      if (stat == 0) allocate (g(kdim + 1), stat=stat)
      if (stat == 0) allocate (cs(kdim), stat=stat)
      if (stat == 0) allocate (sn(kdim), stat=stat)
      if (stat == 0) allocate (y(kdim), stat=stat)
      if (stat == 0) allocate (r(n), stat=stat)
      if (stat == 0) allocate (w(n), stat=stat)
      if (stat == 0) allocate (z(n), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, context//" failed to allocate GMRES work arrays")
         return
      end if
      threshold = self%tol*b_norm
      call op%apply(trans, x, w)
      r = b - w
      beta = nrm2(r)
      if (.not. ieee_is_finite(beta)) then
         call fatal_error(error, context//" residual of the initial guess is not finite")
         return
      end if
      its = 0
      do
         if (beta <= threshold) exit
         if (its >= self%max_iter .or. .not. ieee_is_finite(beta)) then
            call gmres_statistics(self, its, beta/b_norm)
            call fatal_error(error, context//" GMRES did not converge in "//int_str(its) &
               & //" iterations (relative residual "//real_str(beta/b_norm) &
               & //", tolerance "//real_str(self%tol)//")")
            return
         end if

         ! One cycle of GMRES(kdim) from the current x
         v(:, 1) = r/beta
         g = 0.0_wp
         g(1) = beta
         k = 0
         do j = 1, kdim
            its = its + 1
            k = j
            call self%preconditioner%apply(op, trans, v(:, j), z, error)
            if (allocated(error)) return
            call op%apply(trans, z, w)
            do i = 1, j
               h(i, j) = dot(v(:, i), w)
               w = w - h(i, j)*v(:, i)
            end do
            h_next = nrm2(w)
            do i = 1, j - 1
               tmp = cs(i)*h(i, j) + sn(i)*h(i + 1, j)
               h(i + 1, j) = -sn(i)*h(i, j) + cs(i)*h(i + 1, j)
               h(i, j) = tmp
            end do
            call givens_real(h(j, j), h_next, cs(j), sn(j))
            g(j + 1) = -sn(j)*g(j)
            g(j) = cs(j)*g(j)
            if (abs(g(j + 1)) <= threshold .or. .not. (h_next > 0.0_wp) &
               & .or. its >= self%max_iter) exit
            v(:, j + 1) = w/h_next
         end do

         ! x += M^-1 V y with the least-squares solution y of the cycle
         do i = k, 1, -1
            if (.not. (abs(h(i, i)) > 0.0_wp)) then
               call gmres_statistics(self, its, beta/b_norm)
               call fatal_error(error, context//" GMRES breakdown after "//int_str(its) &
                  & //" iterations; the preconditioned operator is singular")
               return
            end if
            y(i) = (g(i) - dot_product(h(i, i + 1:k), y(i + 1:k)))/h(i, i)
         end do
         call gemv(v(:, 1:k), y(1:k), w)
         call self%preconditioner%apply(op, trans, w, z, error)
         if (allocated(error)) return
         x = x + z
         call op%apply(trans, x, w)
         r = b - w
         beta = nrm2(r)
      end do
      call gmres_statistics(self, its, beta/b_norm)
   end subroutine gmres_real

   !> Restarted right-preconditioned GMRES for one complex operator
   !>
   !> - The operator is applied with the plain transpose; the Krylov basis is
   !>   orthonormal in the Hermitian inner product
   !> - Converged when norm(b - op x) <= tol norm(b), checked with the true
   !>   residual at the start and after every cycle
   !> - Updates the statistics of self
   !>
   !> @param[in,out] self     solver state with a set-up preconditioner
   !> @param[in]     op       complex operator
   !> @param[in]     trans    solve with the transpose if true
   !> @param[in]     b        right-hand side, (npt)
   !> @param[in,out] x        initial guess on entry, solution on exit, (npt)
   !> @param[in]     context  prefix of error messages
   !> @param[out]    error    allocated on non-convergence or breakdown
   subroutine gmres_complex(self, op, trans, b, x, context, error)
      !> Solver state
      class(pcm_gmres_solver), intent(inout) :: self
      !> Complex operator
      class(pcm_complex_operator), intent(in) :: op
      !> Solve with the transpose
      logical, intent(in) :: trans
      !> Right-hand side
      complex(wp), contiguous, intent(in) :: b(:)
      !> Initial guess on entry, solution on exit
      complex(wp), contiguous, intent(inout) :: x(:)
      !> Prefix of error messages
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      complex(wp), allocatable :: v(:, :), h(:, :), g(:), sn(:), y(:), r(:), w(:), z(:)
      real(wp), allocatable :: cs(:)
      complex(wp) :: tmp
      real(wp) :: b_norm, beta, threshold, h_next
      integer :: n, kdim, its, i, j, k, stat

      n = size(b)
      self%n_solves = self%n_solves + 1
      self%last_iterations = 0
      self%last_residual = 0.0_wp
      b_norm = znorm2(b)
      if (.not. ieee_is_finite(b_norm)) then
         call fatal_error(error, context//" right-hand side is not finite")
         return
      end if
      if (b_norm <= 0.0_wp) then
         x = (0.0_wp, 0.0_wp)
         return
      end if

      kdim = max(1, min(self%restart, n))
      allocate (v(n, kdim + 1), stat=stat)
      if (stat == 0) allocate (h(kdim + 1, kdim), stat=stat)
      if (stat == 0) allocate (g(kdim + 1), stat=stat)
      if (stat == 0) allocate (cs(kdim), stat=stat)
      if (stat == 0) allocate (sn(kdim), stat=stat)
      if (stat == 0) allocate (y(kdim), stat=stat)
      if (stat == 0) allocate (r(n), stat=stat)
      if (stat == 0) allocate (w(n), stat=stat)
      if (stat == 0) allocate (z(n), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, context//" failed to allocate GMRES work arrays")
         return
      end if
      threshold = self%tol*b_norm
      call op%apply(trans, x, w)
      r = b - w
      beta = znorm2(r)
      if (.not. ieee_is_finite(beta)) then
         call fatal_error(error, context//" residual of the initial guess is not finite")
         return
      end if
      its = 0
      do
         if (beta <= threshold) exit
         if (its >= self%max_iter .or. .not. ieee_is_finite(beta)) then
            call gmres_statistics(self, its, beta/b_norm)
            call fatal_error(error, context//" GMRES did not converge in "//int_str(its) &
               & //" iterations (relative residual "//real_str(beta/b_norm) &
               & //", tolerance "//real_str(self%tol)//")")
            return
         end if

         ! One cycle of GMRES(kdim) from the current x
         v(:, 1) = r/beta
         g = (0.0_wp, 0.0_wp)
         g(1) = beta
         k = 0
         do j = 1, kdim
            its = its + 1
            k = j
            call self%preconditioner%apply_complex(op, trans, v(:, j), z, error)
            if (allocated(error)) return
            call op%apply(trans, z, w)
            do i = 1, j
               ! Hermitian inner product: dot_product conjugates v
               h(i, j) = dot_product(v(:, i), w)
               w = w - h(i, j)*v(:, i)
            end do
            h_next = znorm2(w)
            do i = 1, j - 1
               tmp = cs(i)*h(i, j) + sn(i)*h(i + 1, j)
               h(i + 1, j) = -conjg(sn(i))*h(i, j) + cs(i)*h(i + 1, j)
               h(i, j) = tmp
            end do
            call givens_complex(h(j, j), h_next, cs(j), sn(j))
            g(j + 1) = -conjg(sn(j))*g(j)
            g(j) = cs(j)*g(j)
            if (abs(g(j + 1)) <= threshold .or. .not. (h_next > 0.0_wp) &
               & .or. its >= self%max_iter) exit
            v(:, j + 1) = w/h_next
         end do

         ! x += M^-1 V y with the least-squares solution y of the cycle
         do i = k, 1, -1
            if (.not. (abs(h(i, i)) > 0.0_wp)) then
               call gmres_statistics(self, its, beta/b_norm)
               call fatal_error(error, context//" GMRES breakdown after "//int_str(its) &
                  & //" iterations; the preconditioned operator is singular")
               return
            end if
            ! plain (unconjugated) sum: back substitution, not an inner product
            y(i) = (g(i) - sum(h(i, i + 1:k)*y(i + 1:k)))/h(i, i)
         end do
         w = matmul(v(:, 1:k), y(1:k))
         call self%preconditioner%apply_complex(op, trans, w, z, error)
         if (allocated(error)) return
         x = x + z
         call op%apply(trans, x, w)
         r = b - w
         beta = znorm2(r)
      end do
      call gmres_statistics(self, its, beta/b_norm)
   end subroutine gmres_complex

   !> Record the iterations and residual of a right-hand side
   !>
   !> @param[in,out] self      solver state
   !> @param[in]     its       iterations of the right-hand side
   !> @param[in]     residual  relative true residual
   subroutine gmres_statistics(self, its, residual)
      !> Solver state
      class(pcm_gmres_solver), intent(inout) :: self
      !> Iterations of the right-hand side
      integer, intent(in) :: its
      !> Relative true residual
      real(wp), intent(in) :: residual

      self%n_iterations = self%n_iterations + its
      self%last_iterations = its
      self%last_residual = residual
   end subroutine gmres_statistics

   !> Real Givens rotation zeroing b against a
   !>
   !> - [c s; -s c] [a; b] = [r; 0]; a is overwritten by r
   !>
   !> @param[in,out] a  first entry, r on exit
   !> @param[in]     b  entry to eliminate
   !> @param[out]    c  cosine
   !> @param[out]    s  sine
   pure subroutine givens_real(a, b, c, s)
      !> First entry, r on exit
      real(wp), intent(inout) :: a
      !> Entry to eliminate
      real(wp), intent(in) :: b
      !> Cosine
      real(wp), intent(out) :: c
      !> Sine
      real(wp), intent(out) :: s

      real(wp) :: r

      if (.not. (abs(b) > 0.0_wp)) then
         c = 1.0_wp
         s = 0.0_wp
         return
      end if
      r = hypot(a, b)
      c = a/r
      s = b/r
      a = r
   end subroutine givens_real

   !> Complex Givens rotation zeroing a real b >= 0 against a
   !>
   !> - [c s; -conjg(s) c] [a; b] = [r; 0] with real c; a is overwritten by r
   !>
   !> @param[in,out] a  first entry, r on exit
   !> @param[in]     b  entry to eliminate, a norm
   !> @param[out]    c  cosine
   !> @param[out]    s  complex sine
   pure subroutine givens_complex(a, b, c, s)
      !> First entry, r on exit
      complex(wp), intent(inout) :: a
      !> Entry to eliminate
      real(wp), intent(in) :: b
      !> Cosine
      real(wp), intent(out) :: c
      !> Complex sine
      complex(wp), intent(out) :: s

      complex(wp) :: phase
      real(wp) :: t

      if (.not. (abs(b) > 0.0_wp)) then
         c = 1.0_wp
         s = (0.0_wp, 0.0_wp)
         return
      end if
      if (.not. (abs(a) > 0.0_wp)) then
         c = 0.0_wp
         s = (1.0_wp, 0.0_wp)
         a = b
         return
      end if
      t = hypot(abs(a), b)
      phase = a/abs(a)
      c = abs(a)/t
      s = phase*b/t
      a = phase*t
   end subroutine givens_complex

   !> Euclidean norm of a complex vector
   !>
   !> @param[in] v  vector
   pure function znorm2(v) result(r)
      !> Vector
      complex(wp), intent(in) :: v(:)
      !> Norm
      real(wp) :: r

      r = sqrt(real(dot_product(v, v), wp))
   end function znorm2

   !> Check that an operator id addresses one of the slots created by reset
   !>
   !> @param[in]  id       operator id
   !> @param[in]  npt      operator dimension
   !> @param[in]  n_slots  number of slots created by the last reset
   !> @param[in]  routine  "[routine]" prefix of error messages
   !> @param[out] error    allocated if the id or dimension is invalid
   subroutine check_operator(id, npt, n_slots, routine, error)
      !> Operator id
      integer, intent(in) :: id
      !> Operator dimension
      integer, intent(in) :: npt
      !> Number of slots
      integer, intent(in) :: n_slots
      !> Prefix of error messages
      character(len=*), intent(in) :: routine
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (id < 1 .or. id > n_slots) then
         call fatal_error(error, routine//" operator id "//int_str(id) &
            & //" is outside 1.."//int_str(n_slots)//"; reset the solver first")
         return
      end if
      if (npt < 0) then
         call fatal_error(error, routine//" operator "//int_str(id) &
            & //" has negative dimension "//int_str(npt))
         return
      end if
   end subroutine check_operator

   !> Turn the info of a factorization into an error
   !>
   !> @param[in]  info     info returned by the factorization
   !> @param[in]  kind     factorization kind, one of the pcm_factor_* kinds
   !> @param[in]  id       operator id
   !> @param[in]  routine  "[routine]" prefix of error messages
   !> @param[out] error    allocated if info is nonzero
   subroutine check_factor_info(info, kind, id, routine, error)
      !> Factorization status
      integer, intent(in) :: info
      !> Factorization kind
      integer, intent(in) :: kind
      !> Operator id
      integer, intent(in) :: id
      !> Prefix of error messages
      character(len=*), intent(in) :: routine
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (info > 0) then
         call fatal_error(error, routine//" operator "//int_str(id) &
            & //" is singular ("//factor_name(kind)//" info = "//int_str(info)//")")
         return
      else if (info < 0) then
         call fatal_error(error, routine//" "//factor_name(kind) &
            & //" rejected argument "//int_str(-info)//" for operator "//int_str(id))
         return
      end if
   end subroutine check_factor_info

   !> Check the shapes of right-hand sides and solutions against an operator
   !>
   !> @param[in]  op       operator
   !> @param[in]  rhs      right-hand sides, (npt, nrhs)
   !> @param[in]  x        solutions, (npt, nrhs)
   !> @param[in]  routine  "[routine]" prefix of error messages
   !> @param[out] error    allocated on a shape mismatch
   subroutine check_block_shapes(op, rhs, x, routine, error)
      !> Operator
      class(pcm_operator), intent(in) :: op
      !> Right-hand sides
      real(wp), intent(in) :: rhs(:, :)
      !> Solutions
      real(wp), intent(in) :: x(:, :)
      !> Prefix of error messages
      character(len=*), intent(in) :: routine
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (size(rhs, 1) /= op%npt .or. any(shape(x) /= shape(rhs))) then
         call fatal_error(error, routine//" right-hand sides ("//int_str(size(rhs, 1)) &
            & //", "//int_str(size(rhs, 2))//") and solutions ("//int_str(size(x, 1)) &
            & //", "//int_str(size(x, 2))//") do not match operator " &
            & //int_str(op%id)//" of dimension "//int_str(op%npt))
         return
      end if
   end subroutine check_block_shapes

   !> LAPACK name of a factorization kind
   !>
   !> @param[in] kind  one of the pcm_factor_* kinds
   pure function factor_name(kind) result(name)
      !> Factorization kind
      integer, intent(in) :: kind
      !> Name of the factorization routine
      character(len=:), allocatable :: name

      select case (kind)
      case (pcm_factor_cholesky)
         name = "potrf"
      case (pcm_factor_symmetric)
         name = "sytrf"
      case (pcm_factor_complex_symmetric)
         name = "zsytrf"
      case (pcm_factor_complex_lu)
         name = "zgetrf"
      case default
         name = "getrf"
      end select
   end function factor_name

   !> Integer as text
   !>
   !> @param[in] i  value
   pure function int_str(i) result(s)
      !> Value
      integer, intent(in) :: i
      !> Decimal representation without blanks
      character(len=:), allocatable :: s

      character(len=24) :: buffer

      write (buffer, "(i0)") i
      s = trim(buffer)
   end function int_str

   !> Real as text in scientific notation
   !>
   !> @param[in] r  value
   pure function real_str(r) result(s)
      !> Value
      real(wp), intent(in) :: r
      !> Representation without leading blanks
      character(len=:), allocatable :: s

      character(len=32) :: buffer

      write (buffer, "(es12.4)") r
      s = trim(adjustl(buffer))
   end function real_str

end module moist_model_continuum_component_pcm_solvers
