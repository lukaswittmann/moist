!> Test suite for the softmax branch weights and their derivatives
!>
!> - Direct checks of [[softmax_weights]], [[softmax_weights_grad]] and
!>   [[softmax_weights_hess]], no cavity involved
!> - Objective quadratic in a parameter vector `q`, per branch `m`:
!>
!>     phi_m(q)  = phi0_m + dphi0(:,m) . q + 0.5 q . ddphi0(:,:,m) q
!>     dphi_m(q) = dphi0(:,m) + ddphi0(:,:,m) q
!>     ddphi_m   = ddphi0(:,:,m)
!>
!> - Every finite difference moves both `phi` and `dphi` with `q`; holding
!>   `dphi` fixed drops the `ddphi` terms from the reference and hides an
!>   error in them
!> - Four branches, balanced weights, width other than one, evaluation at a
!>   non-zero `q`; `fixture_is_balanced` guards these properties
!> - Second derivatives checked twice: against differences of the analytic
!>   gradient and against differences of the weights alone
module test_cavity_drop_branching
   use mctc_env_accuracy, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string
   use moist_cavity_drop_branching, only: branch_weight_type, softmax_weights, &
      & softmax_weights_grad, softmax_weights_hess
   implicit none(type, external)
   private

   public :: collect_cavity_drop_branching

   !> Number of branches of the fixture group
   integer, parameter :: nbranch = 4
   !> Number of parameters of the fixture objective
   integer, parameter :: nparam = 3

   !> Softmax width of the fixture; not one, so a missing `1/sigma` shows
   real(wp), parameter :: sigma_phi = 0.7_wp

   !> Relative tolerance of the sum rules and the symmetry
   !>
   !> - Scaled by the largest entry
   !> - Measured 3.5e-17 (`dweights`), 9.0e-17 (`ddweights`) and 4.2e-17
   !>   (asymmetry) against entries up to 0.19, about 20x below the bound
   real(wp), parameter :: rule_tol = 1.0e-14_wp

   !> Tolerance of comparisons between identical arithmetic
   !>
   !> - Same call on the same inputs, bitwise equality expected; small slack
   !>   for compiler contraction
   real(wp), parameter :: exact_tol = 1.0e-15_wp

   !> Coefficient of the `fd_coef*h^2` bound of every finite-difference check
   !>
   !> - Measured truncation coefficients 0.023 (`dweights`), 0.039 (`ddweights`
   !>   from the gradient) and 0.040 (`ddweights` from the weights alone)
   !> - Headroom between 5x and 9x
   real(wp), parameter :: fd_coef = 0.2_wp

   !> Smallest accepted error ratio between the two steps of one check
   !>
   !> - Steps differ by 10x and every stencil is second order, measured ratio
   !>   99.5 to 100.6
   !> - A wrong analytic term leaves a step-independent error, ratio near one
   real(wp), parameter :: fd_ratio_min = 50.0_wp

   !> Steps of the central differences of the weights and of the gradient
   !>
   !> - Measured worst error 2.3e-8 and 2.3e-10 (`dweights`), 3.9e-8 and
   !>   3.9e-10 (`ddweights` from the gradient); clean `h^2` down to `1e-4`,
   !>   roundoff floor near 1e-11 at `1e-5`
   real(wp), parameter :: first_steps(2) = [1.0e-3_wp, 1.0e-4_wp]

   !> Steps of the second differences of the weights alone
   !>
   !> - Measured worst error 4.0e-6 and 4.0e-8; `1e-4` already roundoff
   !>   dominated at 1.2e-8, against 4.0e-10 expected from truncation
   real(wp), parameter :: second_steps(2) = [1.0e-2_wp, 1.0e-3_wp]

contains

   !> Register the softmax branch-weight tests
   !>
   !> @param[out] testsuite collected tests
   subroutine collect_cavity_drop_branching(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fixture_is_balanced", test_fixture_balanced), &
                  new_unittest("sum_rules", test_sum_rules), &
                  new_unittest("hessian_symmetry", test_symmetry), &
                  new_unittest("hess_first_derivative_matches_grad", test_hess_matches_grad), &
                  new_unittest("dweights_fd_of_weights", test_dweights_fd), &
                  new_unittest("ddweights_fd_of_gradient", test_ddweights_fd_gradient), &
                  new_unittest("ddweights_fd_of_weights", test_ddweights_fd_weights), &
                  new_unittest("type_bound_matches_free", test_type_bound), &
                  new_unittest("zero_width", test_zero_width), &
                  new_unittest("single_branch", test_single_branch) &
                  ]
   end subroutine collect_cavity_drop_branching

   !* ================================================================================= *!
   !*                                      Fixture                                      *!
   !* ================================================================================= *!

   !> Build the coefficients of the quadratic objective
   !>
   !> - Values of order one, no two equal, so a transposed index shows
   !> - `ddphi0` filled from its upper triangle, symmetric by construction
   !>
   !> @param[out] phi0   objective values at `q = 0`
   !> @param[out] dphi0  first derivatives at `q = 0`
   !> @param[out] ddphi0 constant second derivatives
   !> @param[out] q0     evaluation point
   subroutine objective_fixture(phi0, dphi0, ddphi0, q0)
      !> Objective values at the origin
      real(wp), intent(out) :: phi0(nbranch)
      !> First derivatives at the origin
      real(wp), intent(out) :: dphi0(nparam, nbranch)
      !> Constant second derivatives
      real(wp), intent(out) :: ddphi0(nparam, nparam, nbranch)
      !> Evaluation point
      real(wp), intent(out) :: q0(nparam)

      !> Upper triangle per branch, order 11, 22, 33, 12, 13, 23
      real(wp) :: upper(6, nbranch)
      !> Branch index
      integer :: ibranch

      phi0 = [0.42_wp, 0.18_wp, 0.71_wp, 0.55_wp]
      dphi0 = reshape([0.63_wp, -0.27_wp, 0.41_wp, &
                       -0.38_wp, 0.52_wp, 0.19_wp, &
                       0.24_wp, 0.36_wp, -0.58_wp, &
                       -0.47_wp, -0.21_wp, 0.33_wp], [nparam, nbranch])
      upper = reshape([0.83_wp, 0.47_wp, 0.62_wp, 0.21_wp, -0.34_wp, 0.15_wp, &
                       0.39_wp, 0.91_wp, 0.28_wp, -0.26_wp, 0.17_wp, 0.44_wp, &
                       0.68_wp, 0.33_wp, 0.75_wp, 0.31_wp, 0.12_wp, -0.29_wp, &
                       0.52_wp, 0.64_wp, 0.41_wp, -0.18_wp, 0.37_wp, 0.23_wp], [6, nbranch])
      do ibranch = 1, nbranch
         ddphi0(1, 1, ibranch) = upper(1, ibranch)
         ddphi0(2, 2, ibranch) = upper(2, ibranch)
         ddphi0(3, 3, ibranch) = upper(3, ibranch)
         ddphi0(1, 2, ibranch) = upper(4, ibranch)
         ddphi0(2, 1, ibranch) = upper(4, ibranch)
         ddphi0(1, 3, ibranch) = upper(5, ibranch)
         ddphi0(3, 1, ibranch) = upper(5, ibranch)
         ddphi0(2, 3, ibranch) = upper(6, ibranch)
         ddphi0(3, 2, ibranch) = upper(6, ibranch)
      end do
      q0 = [0.35_wp, -0.28_wp, 0.46_wp]
   end subroutine objective_fixture

   !> Evaluate the quadratic objective and its gradient at `q`
   !>
   !> @param[in]  q     parameter vector
   !> @param[out] phi   objective values
   !> @param[out] dphi  first derivatives
   !> @param[out] ddphi second derivatives
   subroutine objective_at(q, phi, dphi, ddphi)
      !> Parameter vector
      real(wp), intent(in) :: q(nparam)
      !> Objective values
      real(wp), intent(out) :: phi(nbranch)
      !> First derivatives
      real(wp), intent(out) :: dphi(nparam, nbranch)
      !> Second derivatives
      real(wp), intent(out) :: ddphi(nparam, nparam, nbranch)

      !> Coefficients at the origin
      real(wp) :: phi0(nbranch), dphi0(nparam, nbranch)
      !> Unused evaluation point of the fixture
      real(wp) :: q0(nparam)
      !> Branch index
      integer :: ibranch

      call objective_fixture(phi0, dphi0, ddphi, q0)
      do ibranch = 1, nbranch
         dphi(:, ibranch) = dphi0(:, ibranch) + matmul(ddphi(:, :, ibranch), q)
         phi(ibranch) = phi0(ibranch) + dot_product(dphi0(:, ibranch), q) &
                        + 0.5_wp*dot_product(q, matmul(ddphi(:, :, ibranch), q))
      end do
   end subroutine objective_at

   !> Evaluate the analytic weights and derivatives at the fixture point
   !>
   !> @param[out] weights   softmax weights
   !> @param[out] dweights  first derivatives
   !> @param[out] ddweights second derivatives
   subroutine reference_hess(weights, dweights, ddweights)
      !> Softmax weights
      real(wp), intent(out) :: weights(nbranch)
      !> First derivatives
      real(wp), intent(out) :: dweights(nparam, nbranch)
      !> Second derivatives
      real(wp), intent(out) :: ddweights(nparam, nparam, nbranch)

      !> Objective at the fixture point
      real(wp) :: phi(nbranch), dphi(nparam, nbranch), ddphi(nparam, nparam, nbranch)
      !> Fixture coefficients, only `q0` used
      real(wp) :: phi0(nbranch), dphi0(nparam, nbranch), ddphi0(nparam, nparam, nbranch)
      !> Evaluation point
      real(wp) :: q0(nparam)

      call objective_fixture(phi0, dphi0, ddphi0, q0)
      call objective_at(q0, phi, dphi, ddphi)
      call softmax_weights_hess(phi, dphi, ddphi, sigma_phi, weights, dweights, ddweights)
   end subroutine reference_hess

   !> Return the fixture evaluation point
   !>
   !> @param[out] q0 evaluation point
   subroutine fixture_point(q0)
      !> Evaluation point
      real(wp), intent(out) :: q0(nparam)

      !> Unused fixture coefficients
      real(wp) :: phi0(nbranch), dphi0(nparam, nbranch), ddphi0(nparam, nparam, nbranch)

      call objective_fixture(phi0, dphi0, ddphi0, q0)
   end subroutine fixture_point

   !> Evaluate the softmax weights alone at `q`
   !>
   !> @param[in] q parameter vector
   function weights_at(q) result(weights)
      !> Parameter vector
      real(wp), intent(in) :: q(nparam)
      !> Softmax weights
      real(wp) :: weights(nbranch)

      !> Objective at `q`
      real(wp) :: phi(nbranch), dphi(nparam, nbranch), ddphi(nparam, nparam, nbranch)

      call objective_at(q, phi, dphi, ddphi)
      call softmax_weights(phi, sigma_phi, weights)
   end function weights_at

   !> Evaluate the analytic first derivatives of [[softmax_weights_grad]] at `q`
   !>
   !> @param[in] q parameter vector
   function gradient_at(q) result(dweights)
      !> Parameter vector
      real(wp), intent(in) :: q(nparam)
      !> First derivatives of the weights
      real(wp) :: dweights(nparam, nbranch)

      !> Objective at `q`, with `dphi` moving along
      real(wp) :: phi(nbranch), dphi(nparam, nbranch), ddphi(nparam, nparam, nbranch)
      !> Softmax weights, unused
      real(wp) :: weights(nbranch)
      !> Constant width
      real(wp) :: dsigma_phi(nparam)

      call objective_at(q, phi, dphi, ddphi)
      dsigma_phi = 0.0_wp
      call softmax_weights_grad(phi, dphi, sigma_phi, dsigma_phi, weights, dweights)
   end function gradient_at

   !> Assert the absolute bound at both steps and second-order error decay
   !>
   !> @param[in,out] error test error
   !> @param[in]     err   worst deviation per step
   !> @param[in]     steps finite-difference steps, coarse first
   !> @param[in]     label prefix identifying the comparison
   subroutine check_fd_pair(error, err, steps, label)
      !> Test error
      type(error_type), allocatable, intent(inout) :: error
      !> Worst deviation per step
      real(wp), intent(in) :: err(2)
      !> Finite-difference steps
      real(wp), intent(in) :: steps(2)
      !> Comparison label
      character(len=*), intent(in) :: label

      !> Step index
      integer :: istep

      do istep = 1, 2
         call check(error, err(istep) <= fd_coef*steps(istep)*steps(istep), &
                    label//": deviation "//to_string(err(istep))//" at h = "// &
                    to_string(steps(istep))//" exceeds "// &
                    to_string(fd_coef*steps(istep)*steps(istep)))
         if (allocated(error)) return
      end do

      call check(error, err(1) >= fd_ratio_min*err(2), &
                 label//": deviation fell from "//to_string(err(1))//" to "// &
                 to_string(err(2))//", not second order")
   end subroutine check_fd_pair

   !* ================================================================================= *!
   !*                                       Tests                                       *!
   !* ================================================================================= *!

   !> Check that the fixture keeps every term of the Hessian live
   !>
   !> - Weights within `[0.1, 0.6]`, none saturated
   !> - `ddphi` contribution to `ddweights` far above every tolerance
   !> - No `ddweights` entry close to zero
   !>
   !> @param[out] error test error
   subroutine test_fixture_balanced(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      !> Analytic weights and derivatives
      real(wp) :: weights(nbranch), dweights(nparam, nbranch), ddweights(nparam, nparam, nbranch)
      !> Same without the `ddphi` term
      real(wp) :: w_flat(nbranch), dw_flat(nparam, nbranch), ddw_flat(nparam, nparam, nbranch)
      !> Objective at the fixture point
      real(wp) :: phi(nbranch), dphi(nparam, nbranch), ddphi(nparam, nparam, nbranch)
      !> Zero second derivatives
      real(wp) :: ddphi_zero(nparam, nparam, nbranch)
      !> Evaluation point
      real(wp) :: q0(nparam)

      call fixture_point(q0)
      call objective_at(q0, phi, dphi, ddphi)
      call softmax_weights_hess(phi, dphi, ddphi, sigma_phi, weights, dweights, ddweights)
      ddphi_zero = 0.0_wp
      call softmax_weights_hess(phi, dphi, ddphi_zero, sigma_phi, w_flat, dw_flat, ddw_flat)

      ! Measured weights 0.134 to 0.464
      call check(error, minval(weights) > 0.1_wp .and. maxval(weights) < 0.6_wp, &
                 "Fixture weights not balanced, range "//to_string(minval(weights))// &
                 " to "//to_string(maxval(weights)))
      if (allocated(error)) return

      ! Measured 0.159
      call check(error, maxval(abs(ddweights - ddw_flat)) > 0.05_wp, &
                 "ddphi contribution to ddweights only "// &
                 to_string(maxval(abs(ddweights - ddw_flat))))
      if (allocated(error)) return

      ! Measured smallest entry 4.5e-3
      call check(error, minval(abs(ddweights)) > 1.0e-3_wp, &
                 "ddweights entry close to zero, "//to_string(minval(abs(ddweights))))
   end subroutine test_fixture_balanced

   !> Check that the weights sum to one and every derivative sums to zero
   !>
   !> @param[out] error test error
   subroutine test_sum_rules(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      !> Analytic weights and derivatives
      real(wp) :: weights(nbranch), dweights(nparam, nbranch), ddweights(nparam, nparam, nbranch)
      !> Branch sums
      real(wp) :: dsum(nparam), ddsum(nparam, nparam)
      !> Parameter indices
      integer :: kparam, lparam

      call reference_hess(weights, dweights, ddweights)
      dsum = sum(dweights, dim=2)
      ddsum = sum(ddweights, dim=3)

      call check(error, sum(weights), 1.0_wp, thr=rule_tol, &
                 message="Weights sum to "//to_string(sum(weights)))
      if (allocated(error)) return

      do kparam = 1, nparam
         call check(error, abs(dsum(kparam)) <= rule_tol*maxval(abs(dweights)), &
                    "sum_m dweights("//to_string(kparam)//", m) = "//to_string(dsum(kparam)))
         if (allocated(error)) return
      end do

      do lparam = 1, nparam
         do kparam = 1, nparam
            call check(error, abs(ddsum(kparam, lparam)) <= rule_tol*maxval(abs(ddweights)), &
                       "sum_m ddweights("//to_string(kparam)//", "//to_string(lparam)// &
                       ", m) = "//to_string(ddsum(kparam, lparam)))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_sum_rules

   !> Check that `ddweights` is symmetric in its parameter indices
   !>
   !> @param[out] error test error
   subroutine test_symmetry(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      !> Analytic weights and derivatives
      real(wp) :: weights(nbranch), dweights(nparam, nbranch), ddweights(nparam, nparam, nbranch)
      !> Parameter and branch indices
      integer :: kparam, lparam, ibranch

      call reference_hess(weights, dweights, ddweights)

      do ibranch = 1, nbranch
         do lparam = 1, nparam
            do kparam = lparam + 1, nparam
               call check(error, ddweights(kparam, lparam, ibranch), &
                          ddweights(lparam, kparam, ibranch), &
                          thr=rule_tol*maxval(abs(ddweights)), &
                          message="ddweights asymmetric in ("//to_string(kparam)//", "// &
                          to_string(lparam)//") for branch "//to_string(ibranch))
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_symmetry

   !> Check the Hessian routine against [[softmax_weights_grad]]
   !>
   !> - Same weights and first derivatives expected
   !>
   !> @param[out] error test error
   subroutine test_hess_matches_grad(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      !> Analytic weights and derivatives from the Hessian routine
      real(wp) :: weights(nbranch), dweights(nparam, nbranch), ddweights(nparam, nparam, nbranch)
      !> Weights from the plain routine
      real(wp) :: w_plain(nbranch)
      !> Weights and first derivatives from the gradient routine
      real(wp) :: w_grad(nbranch), dw_grad(nparam, nbranch)
      !> Objective at the fixture point
      real(wp) :: phi(nbranch), dphi(nparam, nbranch), ddphi(nparam, nparam, nbranch)
      !> Constant width
      real(wp) :: dsigma_phi(nparam)
      !> Evaluation point
      real(wp) :: q0(nparam)

      call fixture_point(q0)
      call objective_at(q0, phi, dphi, ddphi)
      call softmax_weights_hess(phi, dphi, ddphi, sigma_phi, weights, dweights, ddweights)
      call softmax_weights(phi, sigma_phi, w_plain)
      dsigma_phi = 0.0_wp
      call softmax_weights_grad(phi, dphi, sigma_phi, dsigma_phi, w_grad, dw_grad)

      call check(error, maxval(abs(weights - w_plain)) <= exact_tol, &
                 "Weights differ from softmax_weights by "// &
                 to_string(maxval(abs(weights - w_plain))))
      if (allocated(error)) return

      call check(error, maxval(abs(weights - w_grad)) <= exact_tol, &
                 "Weights differ from softmax_weights_grad by "// &
                 to_string(maxval(abs(weights - w_grad))))
      if (allocated(error)) return

      call check(error, maxval(abs(dweights - dw_grad)) <= exact_tol, &
                 "dweights differ from softmax_weights_grad by "// &
                 to_string(maxval(abs(dweights - dw_grad))))
   end subroutine test_hess_matches_grad

   !> Check `dweights` against central differences of the weights
   !>
   !> @param[out] error test error
   subroutine test_dweights_fd(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      !> Analytic weights and derivatives
      real(wp) :: weights(nbranch), dweights(nparam, nbranch), ddweights(nparam, nparam, nbranch)
      !> Evaluation point and displaced point
      real(wp) :: q0(nparam), q(nparam)
      !> Central difference of the weights
      real(wp) :: fd(nbranch)
      !> Worst deviation per step
      real(wp) :: err(2)
      !> Current step
      real(wp) :: h
      !> Step and parameter indices
      integer :: istep, kparam

      call fixture_point(q0)
      call reference_hess(weights, dweights, ddweights)

      err = 0.0_wp
      do istep = 1, 2
         h = first_steps(istep)
         do kparam = 1, nparam
            q = q0
            q(kparam) = q0(kparam) + h
            fd = weights_at(q)
            q(kparam) = q0(kparam) - h
            fd = (fd - weights_at(q))/(2.0_wp*h)
            err(istep) = max(err(istep), maxval(abs(fd - dweights(kparam, :))))
         end do
      end do

      call check_fd_pair(error, err, first_steps, "dweights vs FD of weights")
   end subroutine test_dweights_fd

   !> Check `ddweights` against central differences of the analytic gradient
   !>
   !> - `ddweights(k, l, :)` against the difference in `q_l` of `dweights(k, :)`
   !> - Both `phi` and `dphi` displaced, see [[gradient_at]]
   !>
   !> @param[out] error test error
   subroutine test_ddweights_fd_gradient(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      !> Analytic weights and derivatives
      real(wp) :: weights(nbranch), dweights(nparam, nbranch), ddweights(nparam, nparam, nbranch)
      !> Evaluation point and displaced point
      real(wp) :: q0(nparam), q(nparam)
      !> Central difference of the gradient in one parameter
      real(wp) :: fd(nparam, nbranch)
      !> Worst deviation per step
      real(wp) :: err(2)
      !> Current step
      real(wp) :: h
      !> Step and parameter indices
      integer :: istep, kparam, lparam

      call fixture_point(q0)
      call reference_hess(weights, dweights, ddweights)

      err = 0.0_wp
      do istep = 1, 2
         h = first_steps(istep)
         do lparam = 1, nparam
            q = q0
            q(lparam) = q0(lparam) + h
            fd = gradient_at(q)
            q(lparam) = q0(lparam) - h
            fd = (fd - gradient_at(q))/(2.0_wp*h)
            do kparam = 1, nparam
               err(istep) = max(err(istep), &
                                maxval(abs(fd(kparam, :) - ddweights(kparam, lparam, :))))
            end do
         end do
      end do

      call check_fd_pair(error, err, first_steps, "ddweights vs FD of gradient")
   end subroutine test_ddweights_fd_gradient

   !> Check `ddweights` against second differences of the weights alone
   !>
   !> - Four-point mixed stencil off the diagonal, three-point on it
   !> - Independent of [[softmax_weights_grad]]
   !>
   !> @param[out] error test error
   subroutine test_ddweights_fd_weights(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      !> Analytic weights and derivatives
      real(wp) :: weights(nbranch), dweights(nparam, nbranch), ddweights(nparam, nparam, nbranch)
      !> Evaluation point and displaced point
      real(wp) :: q0(nparam), q(nparam)
      !> Second difference of the weights
      real(wp) :: fd(nbranch)
      !> Worst deviation per step
      real(wp) :: err(2)
      !> Current step
      real(wp) :: h
      !> Step and parameter indices
      integer :: istep, kparam, lparam

      call fixture_point(q0)
      call reference_hess(weights, dweights, ddweights)

      err = 0.0_wp
      do istep = 1, 2
         h = second_steps(istep)
         do lparam = 1, nparam
            do kparam = 1, nparam
               q = q0
               if (kparam == lparam) then
                  q(kparam) = q0(kparam) + h
                  fd = weights_at(q)
                  q(kparam) = q0(kparam) - h
                  fd = (fd + weights_at(q) - 2.0_wp*weights_at(q0))/(h*h)
               else
                  q(kparam) = q0(kparam) + h
                  q(lparam) = q0(lparam) + h
                  fd = weights_at(q)
                  q(lparam) = q0(lparam) - h
                  fd = fd - weights_at(q)
                  q(kparam) = q0(kparam) - h
                  fd = fd + weights_at(q)
                  q(lparam) = q0(lparam) + h
                  fd = (fd - weights_at(q))/(4.0_wp*h*h)
               end if
               err(istep) = max(err(istep), maxval(abs(fd - ddweights(kparam, lparam, :))))
            end do
         end do
      end do

      call check_fd_pair(error, err, second_steps, "ddweights vs FD of weights")
   end subroutine test_ddweights_fd_weights

   !> Check the type-bound `weights_hess` against the free routine
   !>
   !> @param[out] error test error
   subroutine test_type_bound(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      !> Branch-weight model under test
      type(branch_weight_type) :: model
      !> Free-routine weights and derivatives
      real(wp) :: weights(nbranch), dweights(nparam, nbranch), ddweights(nparam, nparam, nbranch)
      !> Type-bound weights and derivatives
      real(wp) :: w_tb(nbranch), dw_tb(nparam, nbranch), ddw_tb(nparam, nparam, nbranch)
      !> Objective at the fixture point
      real(wp) :: phi(nbranch), dphi(nparam, nbranch), ddphi(nparam, nparam, nbranch)
      !> Evaluation point
      real(wp) :: q0(nparam)

      call fixture_point(q0)
      call objective_at(q0, phi, dphi, ddphi)
      call softmax_weights_hess(phi, dphi, ddphi, sigma_phi, weights, dweights, ddweights)

      call model%init(sigma_phi)
      call model%weights_hess(phi, dphi, ddphi, w_tb, dw_tb, ddw_tb)

      call check(error, maxval(abs(w_tb - weights)) <= exact_tol, &
                 "Type-bound weights differ by "//to_string(maxval(abs(w_tb - weights))))
      if (allocated(error)) return

      call check(error, maxval(abs(dw_tb - dweights)) <= exact_tol, &
                 "Type-bound dweights differ by "//to_string(maxval(abs(dw_tb - dweights))))
      if (allocated(error)) return

      call check(error, maxval(abs(ddw_tb - ddweights)) <= exact_tol, &
                 "Type-bound ddweights differ by "// &
                 to_string(maxval(abs(ddw_tb - ddweights))))
   end subroutine test_type_bound

   !> Check that a vanishing width gives zero first and second derivatives
   !>
   !> - Weights not asserted, `0/0` in [[softmax_weights]] for the lowest branch
   !>
   !> @param[out] error test error
   subroutine test_zero_width(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      !> Weights and derivatives at zero width
      real(wp) :: weights(nbranch), dweights(nparam, nbranch), ddweights(nparam, nparam, nbranch)
      !> Objective at the fixture point
      real(wp) :: phi(nbranch), dphi(nparam, nbranch), ddphi(nparam, nparam, nbranch)
      !> Evaluation point
      real(wp) :: q0(nparam)

      call fixture_point(q0)
      call objective_at(q0, phi, dphi, ddphi)
      call softmax_weights_hess(phi, dphi, ddphi, 0.0_wp, weights, dweights, ddweights)

      call check(error, all(dweights == 0.0_wp), &
                 "dweights not zero at zero width, max "//to_string(maxval(abs(dweights))))
      if (allocated(error)) return

      call check(error, all(ddweights == 0.0_wp), &
                 "ddweights not zero at zero width, max "//to_string(maxval(abs(ddweights))))
   end subroutine test_zero_width

   !> Check that a single branch has weight one and zero derivatives
   !>
   !> - Non-zero `dphi` and `ddphi`, so the zero is a cancellation
   !>
   !> @param[out] error test error
   subroutine test_single_branch(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      !> Objective of the single branch
      real(wp) :: phi(1), dphi(nparam, 1), ddphi(nparam, nparam, 1)
      !> Weights and derivatives of the single branch
      real(wp) :: weights(1), dweights(nparam, 1), ddweights(nparam, nparam, 1)

      phi = [0.42_wp]
      dphi(:, 1) = [0.63_wp, -0.27_wp, 0.41_wp]
      ddphi(:, :, 1) = reshape([0.83_wp, 0.21_wp, -0.34_wp, &
                                0.21_wp, 0.47_wp, 0.15_wp, &
                                -0.34_wp, 0.15_wp, 0.62_wp], [nparam, nparam])

      call softmax_weights_hess(phi, dphi, ddphi, sigma_phi, weights, dweights, ddweights)

      call check(error, weights(1), 1.0_wp, thr=exact_tol, &
                 message="Single-branch weight is "//to_string(weights(1)))
      if (allocated(error)) return

      call check(error, maxval(abs(dweights)) <= exact_tol, &
                 "Single-branch dweights not zero, max "//to_string(maxval(abs(dweights))))
      if (allocated(error)) return

      call check(error, maxval(abs(ddweights)) <= exact_tol, &
                 "Single-branch ddweights not zero, max "//to_string(maxval(abs(ddweights))))
   end subroutine test_single_branch

end module test_cavity_drop_branching
