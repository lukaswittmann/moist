!> Access to generated conductor matrix values and derivatives for regression tests
!>
!> The public solver callback receives each assembled operator. For derivative
!> probes we choose x and supply v = A x, so its solution is known exactly.
!> At feps = 1 the generated geometry gradient is (1/2) x^T (dA) x.
!> With x = q1, a potential tangent dv = A q2 gives q1^T (dA) q2 directly.
!> This avoids subtracting large diagonal weights to recover pair derivatives
!> and introduces no linear-solve error into the existing matrix comparisons.
module test_pcm_conductor_probe
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_boys, only: boys_x_far
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver, pcm_operator
   use pcmgen_gauss_conductor, only: conductor_state, &
      & conductor_assemble, conductor_medium, conductor_solve, conductor_weights, conductor_tangent
   implicit none(type, external)
   private

   public :: count_branches, nmol, nleb_survey
   public :: probe_conductor_matrix, probe_conductor_matrix_gradient
   public :: probe_conductor_bilinear_weights, probe_conductor_nuclear_gradient
   public :: probe_conductor_fixed_seed_tangent

   !> Number of mstore structures drawn per test
   integer, parameter :: nmol = 10

   !> Lebedev grid used for the survey tests
   integer, parameter :: nleb_survey = 50

   !> Capture an operator and return the known solution of a manufactured RHS
   type, extends(pcm_solver) :: conductor_probe_solver
      !> Assembled operator
      real(wp), allocatable :: matrix(:, :)
      !> Prescribed solution, only used for v = A seed
      real(wp), allocatable :: seed(:)
      !> Residual allowance for cancellation in a manufactured zero tangent
      real(wp) :: rhs_atol = 0.0_wp
   contains
      procedure :: setup => probe_setup
      procedure :: solve_block => probe_solve_block
   end type conductor_probe_solver

contains

   !> Count the pairs on either side of the erf-saturation threshold
   !>
   !> @param[in]  xi    Gaussian widths (ngrid)
   !> @param[in]  xyz   Surface positions (3, ngrid)
   !> @param[out] nfar  Number of saturated pairs
   !> @param[out] nnear Number of unsaturated pairs
   subroutine count_branches(xi, xyz, nfar, nnear)
      !> Gaussian inverse lengths
      real(wp), intent(in) :: xi(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Saturated pair count
      integer, intent(out) :: nfar
      !> Unsaturated pair count
      integer, intent(out) :: nnear

      !> Pair indices and grid size
      integer :: i, j, ngrid
      !> Squared separation
      real(wp) :: r2
      !> Per-point saturation bounds
      real(wp), allocatable :: bound(:)

      ngrid = size(xi)
      allocate (bound(ngrid))
      bound = boys_x_far/(xi*xi)

      nfar = 0
      nnear = 0
      do i = 1, ngrid
         do j = 1, i - 1
            r2 = sum((xyz(:, i) - xyz(:, j))**2)
            if (r2 >= bound(i) + bound(j)) then
               nfar = nfar + 1
            else
               nnear = nnear + 1
            end if
         end do
      end do

   end subroutine count_branches

   !> Copy the operator supplied by the generated assembly phase
   !>
   !> @param[in,out] self Operator capture and prescribed solution
   !> @param[in] op Symmetric operator of the current surface
   !> @param[out] error Allocated on failure
   subroutine probe_setup(self, op, error)
      !> Probe solver
      class(conductor_probe_solver), intent(inout) :: self
      !> Generated conductor operator
      class(pcm_operator), intent(in) :: op
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer :: stat

      if (allocated(self%matrix)) deallocate (self%matrix)
      allocate (self%matrix(op%npt, op%npt), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Conductor probe matrix allocation failed")
         return
      end if
      call op%assemble(self%matrix, error)
   end subroutine probe_setup

   !> Return the prescribed solution after checking the manufactured RHS
   !>
   !> @param[in,out] self Operator capture and prescribed solution
   !> @param[in] op Symmetric operator of the current surface
   !> @param[in] trans Must be false for the manufactured solve
   !> @param[in] rhs A times the prescribed solution, (ngrid, 1)
   !> @param[in,out] x Prescribed solution on return, (ngrid, 1)
   !> @param[out] error Allocated on failure
   subroutine probe_solve_block(self, op, trans, rhs, x, error)
      !> Probe solver
      class(conductor_probe_solver), intent(inout) :: self
      !> Generated conductor operator
      class(pcm_operator), intent(in) :: op
      !> Transpose flag
      logical, intent(in) :: trans
      !> Manufactured right-hand side
      real(wp), contiguous, intent(in) :: rhs(:, :)
      !> Solution block
      real(wp), contiguous, intent(inout) :: x(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(self%seed) .or. .not. allocated(self%matrix)) then
         call fatal_error(error, "Conductor probe solve requires a prescribed solution")
         return
      end if
      if (size(rhs, 2) /= 1 .or. size(rhs, 1) /= size(self%seed) &
          .or. any(shape(x) /= shape(rhs)) .or. .not. op%symmetric .or. trans) then
         call fatal_error(error, "Unexpected conductor probe solve")
         return
      end if
      if (.not. all(abs(rhs(:, 1) - matmul(self%matrix, self%seed)) <= self%rhs_atol)) then
         call fatal_error(error, "Conductor probe RHS is not A times the prescribed solution")
         return
      end if
      x(:, 1) = self%seed
   end subroutine probe_solve_block

   !> Read the generated operator without factoring it
   !>
   !> @param[in] xi Positive widths (ngrid), bohr^-1
   !> @param[in] f Positive dimensionless factors (ngrid)
   !> @param[in] xyz Cartesian coordinates (3, ngrid), bohr
   !> @param[out] matrix Output matrix (ngrid, ngrid)
   !> @param[out] error Allocated on failure
   subroutine probe_conductor_matrix(xi, f, xyz, matrix, error)
      !> Gaussian inverse lengths
      real(wp), intent(in) :: xi(:)
      !> Switching factors
      real(wp), intent(in) :: f(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Dense conductor operator
      real(wp), intent(out) :: matrix(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(conductor_state) :: state
      type(conductor_probe_solver) :: solver

      call conductor_assemble(state, xyz, xi, f, error, solver)
      if (allocated(error)) return
      select type (captured => state%solver)
      type is (conductor_probe_solver)
         if (any(shape(matrix) /= shape(captured%matrix))) then
            call fatal_error(error, "Conductor probe matrix output has the wrong shape")
            return
         end if
         matrix = captured%matrix
      end select
   end subroutine probe_conductor_matrix

   !> Geometry gradient for a known solution of the generated conductor system
   !>
   !> @param[in,out] state Assembled operator with unit dielectric scaling
   !> @param[in] seed Manufactured solution vector (ngrid)
   !> @param[out] w_xi Surface weights (ngrid)
   !> @param[out] w_f Surface weights (ngrid)
   !> @param[out] w_xyz Cartesian surface weights (3, ngrid)
   !> @param[out] error Allocated on failure
   subroutine seed_weights(state, seed, w_xi, w_f, w_xyz, error)
      !> Generated conductor state
      type(conductor_state), intent(inout) :: state
      !> Prescribed solution
      real(wp), intent(in) :: seed(:)
      !> Width derivatives
      real(wp), intent(out) :: w_xi(:)
      !> Switching-factor derivatives
      real(wp), intent(out) :: w_f(:)
      !> Position derivatives
      real(wp), intent(out) :: w_xyz(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: v(:), q(:), grad_v(:)
      real(wp) :: energy, grad_feps
      integer :: stat

      allocate (q(size(seed)), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Conductor probe solve allocation failed")
         return
      end if
      allocate (grad_v(size(seed)), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Conductor probe solve allocation failed")
         return
      end if
      select type (solver => state%solver)
      type is (conductor_probe_solver)
         solver%seed = seed
         v = matmul(solver%matrix, seed)
      class default
         call fatal_error(error, "Conductor probe state has the wrong solver")
         return
      end select
      call conductor_solve(state, v, energy, q, grad_v, error)
      if (allocated(error)) return
      call conductor_weights(state, w_xyz, w_xi, w_f, grad_feps, error)
   end subroutine seed_weights

   !> Generated surface derivatives of q1^T A q2 at fixed q1 and q2
   !>
   !> @param[in] xi Positive widths (ngrid), bohr^-1
   !> @param[in] f Positive dimensionless factors (ngrid)
   !> @param[in] xyz Cartesian coordinates (3, ngrid), bohr
   !> @param[in] q1 Fixed vector (ngrid)
   !> @param[in] q2 Fixed vector (ngrid)
   !> @param[out] w_xi Surface weights (ngrid)
   !> @param[out] w_f Surface weights (ngrid)
   !> @param[out] w_xyz Cartesian surface weights (3, ngrid)
   !> @param[out] error Allocated on failure
   subroutine probe_conductor_bilinear_weights(xi, f, xyz, q1, q2, w_xi, w_f, w_xyz, error)
      !> Gaussian inverse lengths
      real(wp), intent(in) :: xi(:)
      !> Switching factors
      real(wp), intent(in) :: f(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Left contraction vector
      real(wp), intent(in) :: q1(:)
      !> Right contraction vector
      real(wp), intent(in) :: q2(:)
      !> Width derivatives
      real(wp), intent(out) :: w_xi(:)
      !> Switching-factor derivatives
      real(wp), intent(out) :: w_f(:)
      !> Position derivatives
      real(wp), intent(out) :: w_xyz(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(conductor_state) :: state
      type(conductor_probe_solver) :: solver
      real(wp) :: zero_r(3, size(xi), 1), zero_s(size(xi), 1)
      real(wp) :: dv(size(xi), 1), hvp_r(3, size(xi), 1)
      real(wp) :: hvp_xi(size(xi), 1), hvp_f(size(xi), 1), hvp_v(size(xi), 1)

      call conductor_assemble(state, xyz, xi, f, error, solver)
      if (allocated(error)) return
      call conductor_medium(state, 1.0_wp, error)
      if (allocated(error)) return
      call seed_weights(state, q1, w_xi, w_f, w_xyz, error)
      if (allocated(error)) return
      select type (captured => state%solver)
      type is (conductor_probe_solver)
         captured%seed = q2
         dv(:, 1) = matmul(captured%matrix, q2)
      class default
         call fatal_error(error, "Conductor probe state has the wrong solver")
         return
      end select
      zero_r = 0.0_wp
      zero_s = 0.0_wp
      call conductor_tangent(state, zero_r, zero_s, zero_s, dv, hvp_r, hvp_xi, hvp_f, hvp_v, error)
      if (allocated(error)) return
      w_xi = hvp_xi(:, 1)
      w_f = hvp_f(:, 1)
      w_xyz = hvp_r(:, :, 1)
   end subroutine probe_conductor_bilinear_weights

   !> Geometry Hessian of (1/2) seed^T A seed with a fixed solution
   !>
   !> Supply dv = (dA) seed so the solution tangent vanishes. The callback
   !> checks its zero RHS to roundoff before returning the prescribed zero.
   !> Matrix derivatives come from the public bilinear probe, not a formula.
   !>
   !> @param[in] xi Positive widths (ngrid), bohr^-1
   !> @param[in] f Positive dimensionless factors (ngrid)
   !> @param[in] xyz Cartesian coordinates (3, ngrid), bohr
   !> @param[in] seed Fixed vector (ngrid)
   !> @param[in] dr Position directions (3, ngrid, nblk)
   !> @param[in] dxi Width directions (ngrid, nblk)
   !> @param[in] df Switching directions (ngrid, nblk)
   !> @param[out] h_r Position-gradient tangents (3, ngrid, nblk)
   !> @param[out] h_xi Width-gradient tangents (ngrid, nblk)
   !> @param[out] h_f Switching-gradient tangents (ngrid, nblk)
   !> @param[out] error Allocated on failure
   subroutine probe_conductor_fixed_seed_tangent(xi, f, xyz, seed, dr, dxi, df, h_r, h_xi, h_f, error)
      !> Gaussian inverse lengths
      real(wp), intent(in) :: xi(:)
      !> Switching factors
      real(wp), intent(in) :: f(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Prescribed solution
      real(wp), intent(in) :: seed(:)
      !> Position directions
      real(wp), intent(in) :: dr(:, :, :)
      !> Width directions
      real(wp), intent(in) :: dxi(:, :)
      !> Switching directions
      real(wp), intent(in) :: df(:, :)
      !> Position-gradient tangents
      real(wp), intent(out) :: h_r(:, :, :)
      !> Width-gradient tangents
      real(wp), intent(out) :: h_xi(:, :)
      !> Switching-gradient tangents
      real(wp), intent(out) :: h_f(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(conductor_state) :: state
      type(conductor_probe_solver) :: solver
      real(wp) :: basis(size(xi)), w_xi(size(xi)), w_f(size(xi)), w_r(3, size(xi))
      real(wp) :: dv(size(xi), size(dxi, 2)), h_v(size(xi), size(dxi, 2))
      real(wp) :: rhs_scale
      integer :: point, direction

      call conductor_assemble(state, xyz, xi, f, error, solver)
      if (allocated(error)) return
      call conductor_medium(state, 1.0_wp, error)
      if (allocated(error)) return
      call seed_weights(state, seed, w_xi, w_f, w_r, error)
      if (allocated(error)) return
      do point = 1, size(xi)
         basis = 0.0_wp
         basis(point) = 1.0_wp
         call probe_conductor_bilinear_weights(xi, f, xyz, basis, seed, w_xi, w_f, w_r, error)
         if (allocated(error)) return
         do direction = 1, size(dxi, 2)
            dv(point, direction) = dot_product(w_xi, dxi(:, direction)) &
               & + dot_product(w_f, df(:, direction)) + sum(w_r*dr(:, :, direction))
         end do
      end do
      select type (captured => state%solver)
      type is (conductor_probe_solver)
         captured%seed = 0.0_wp
      end select
      do direction = 1, size(dxi, 2)
         rhs_scale = max(tiny(1.0_wp), maxval(abs(dv(:, direction))))
         select type (captured => state%solver)
         type is (conductor_probe_solver)
            captured%rhs_atol = 32.0_wp*epsilon(1.0_wp)*rhs_scale
         end select
         call conductor_tangent(state, dr(:, :, direction:direction), dxi(:, direction:direction), &
            & df(:, direction:direction), dv(:, direction:direction), h_r(:, :, direction:direction), &
            & h_xi(:, direction:direction), h_f(:, direction:direction), h_v(:, direction:direction), error)
         if (allocated(error)) return
      end do
   end subroutine probe_conductor_fixed_seed_tangent

   !> Contract surface derivatives into nuclear coordinates
   !>
   !> @param[in] xi1_rA Forward derivatives (3, nat, ngrid)
   !> @param[in] f1_rA Forward derivatives (3, nat, ngrid)
   !> @param[in] xyz1_rA Forward derivatives (3, 3, nat, ngrid)
   !> @param[in] w_xi Surface weights (ngrid)
   !> @param[in] w_f Surface weights (ngrid)
   !> @param[in] w_xyz Cartesian surface weights (3, ngrid)
   !> @param[out] grad_rA Output gradient (3, nat)
   !> @param[out] error Allocated on failure
   subroutine probe_conductor_nuclear_gradient(xi1_rA, f1_rA, xyz1_rA, w_xi, w_f, w_xyz, grad_rA, error)
      !> Nuclear width derivatives
      real(wp), intent(in) :: xi1_rA(:, :, :)
      !> Nuclear switching-factor derivatives
      real(wp), intent(in) :: f1_rA(:, :, :)
      !> Nuclear position derivatives
      real(wp), intent(in) :: xyz1_rA(:, :, :, :)
      !> Width derivatives
      real(wp), intent(in) :: w_xi(:)
      !> Switching-factor derivatives
      real(wp), intent(in) :: w_f(:)
      !> Position derivatives
      real(wp), intent(in) :: w_xyz(:, :)
      !> Contracted nuclear gradient
      real(wp), intent(out) :: grad_rA(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer :: atom, axis, point

      if (any(shape(grad_rA) /= [3, size(xi1_rA, 2)])) then
         call fatal_error(error, "Conductor probe nuclear gradient has the wrong shape")
         return
      end if
      grad_rA = 0.0_wp
      do point = 1, size(w_xi)
         do atom = 1, size(grad_rA, 2)
            do axis = 1, 3
               grad_rA(axis, atom) = grad_rA(axis, atom) &
                  & + w_xi(point)*xi1_rA(axis, atom, point) &
                  & + w_f(point)*f1_rA(axis, atom, point) &
                  & + dot_product(xyz1_rA(:, axis, atom, point), w_xyz(:, point))
            end do
         end do
      end do
   end subroutine probe_conductor_nuclear_gradient

   !> Read matrix values and recover each entry's generated nuclear derivative
   !>
   !> Each entry depends on one point or one pair. Two-point probes retain
   !> quadratic cost without adding a derivative interface to production code.
   !>
   !> @param[in] xi Positive widths (ngrid), bohr^-1
   !> @param[in] f Positive dimensionless factors (ngrid)
   !> @param[in] xyz Cartesian coordinates (3, ngrid), bohr
   !> @param[in] xi1_rA Forward derivatives (3, nat, ngrid)
   !> @param[in] f1_rA Forward derivatives (3, nat, ngrid)
   !> @param[in] xyz1_rA Forward derivatives (3, 3, nat, ngrid)
   !> @param[out] matrix Output matrix (ngrid, ngrid)
   !> @param[out] derivative Output derivatives (3, nat, ngrid, ngrid)
   !> @param[out] error Allocated on failure
   subroutine probe_conductor_matrix_gradient(xi, f, xyz, xi1_rA, f1_rA, xyz1_rA, matrix, derivative, error)
      !> Gaussian inverse lengths
      real(wp), intent(in) :: xi(:)
      !> Switching factors
      real(wp), intent(in) :: f(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Nuclear width derivatives
      real(wp), intent(in) :: xi1_rA(:, :, :)
      !> Nuclear switching-factor derivatives
      real(wp), intent(in) :: f1_rA(:, :, :)
      !> Nuclear position derivatives
      real(wp), intent(in) :: xyz1_rA(:, :, :, :)
      !> Dense conductor operator
      real(wp), intent(out) :: matrix(:, :)
      !> Nuclear matrix derivatives
      real(wp), intent(out) :: derivative(:, :, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: w_xi(2), w_f(2), w_xyz(3, 2)
      real(wp) :: grad(3, size(xi1_rA, 2))
      integer :: i, j, pair(2)

      call probe_conductor_matrix(xi, f, xyz, matrix, error)
      if (allocated(error)) return
      do j = 1, size(xi)
         call probe_conductor_bilinear_weights(xi(j:j), f(j:j), xyz(:, j:j), &
            & [1.0_wp], [1.0_wp], w_xi(:1), w_f(:1), w_xyz(:, :1), error)
         if (allocated(error)) return
         call probe_conductor_nuclear_gradient(xi1_rA(:, :, j:j), f1_rA(:, :, j:j), &
            & xyz1_rA(:, :, :, j:j), w_xi(:1), w_f(:1), w_xyz(:, :1), grad, error)
         if (allocated(error)) return
         derivative(:, :, j, j) = grad
         do i = 1, j - 1
            pair = [i, j]
            call probe_conductor_bilinear_weights(xi(pair), f(pair), xyz(:, pair), &
               & [1.0_wp, 0.0_wp], [0.0_wp, 1.0_wp], w_xi, w_f, w_xyz, error)
            if (allocated(error)) return
            call probe_conductor_nuclear_gradient(xi1_rA(:, :, pair), f1_rA(:, :, pair), &
               & xyz1_rA(:, :, :, pair), w_xi, w_f, w_xyz, grad, error)
            if (allocated(error)) return
            derivative(:, :, i, j) = grad
            derivative(:, :, j, i) = grad
         end do
      end do
   end subroutine probe_conductor_matrix_gradient

end module test_pcm_conductor_probe
