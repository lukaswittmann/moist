module moist_cavity_drop_branching
   use mctc_env_accuracy, only: wp
   implicit none
   private

   public :: branch_weight_type
   public :: softmax_weights
   public :: softmax_weights_grad

   !> Branch-weight model for degenerate closest-point branches
   !>
   !> Uses the softmax model
   !> p_m = exp(-Phi_m/sigma_phi) / sum_n exp(-Phi_n/sigma_phi)
   !> with softmax width sigma_phi = s.
   type :: branch_weight_type
      !> Softmax scale parameter; also serves as the softmax width sigma_phi.
      real(wp) :: s = 1.0_wp
   contains
      !> Initialize softmax scale parameter.
      procedure :: init => branch_weight_init
      !> Compute branch weights for one branch group.
      procedure :: weights => branch_weight_weights
      !> Compute Shannon entropy of branch weights for one branch group.
      procedure :: branch_weights_entropy
      !> Compute branch weights and derivatives for one branch group.
      procedure :: weights_grad => branch_weight_weights_grad
   end type branch_weight_type

contains

   !> Initialize branch-weight model
   !> @param[inout] self Branch-weight instance
   !> @param[in]    s    Softmax scale parameter
   subroutine branch_weight_init(self, s)
      class(branch_weight_type), intent(inout) :: self
      real(wp), intent(in) :: s

      self%s = s
   end subroutine branch_weight_init

   !> Compute softmax branch weights from objective values
   !> @param[in]  phi       Objective values Phi_m for one branch group
   !> @param[in]  sigma_phi Softmax width
   !> @param[out] weights   Softmax weights p_m
   pure subroutine softmax_weights(phi, sigma_phi, weights)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: sigma_phi
      real(wp), intent(out) :: weights(:)

      real(wp) :: phi_min, z
      real(wp) :: exponents(size(phi))

      phi_min = minval(phi)
      exponents = exp(-(phi - phi_min)/sigma_phi)
      z = sum(exponents)

      if (z <= tiny(1.0_wp)) then
         weights = 1.0_wp/real(size(phi), wp)
      else
         weights = exponents/z
      end if
   end subroutine softmax_weights

   !> Compute softmax weights and their derivatives.
   !>
   !> dweights(k,m) corresponds to derivative of p_m w.r.t. k
   !> @param[in]  phi        Objective values Phi_m (nbranch)
   !> @param[in]  dphi       Derivatives dPhi/dq_k (nparam,nbranch)
   !> @param[in]  sigma_phi  Softmax width
   !> @param[in]  dsigma_phi Derivatives dsigma_phi/dq_k (nparam)
   !> @param[out] weights    Softmax weights p_m (nbranch)
   !> @param[out] dweights   Derivatives dp_m/dq_k (nparam,nbranch)
   pure subroutine softmax_weights_grad(phi, dphi, sigma_phi, dsigma_phi, weights, dweights)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: dphi(:, :)
      real(wp), intent(in) :: sigma_phi
      real(wp), intent(in) :: dsigma_phi(:)
      real(wp), intent(out) :: weights(:)
      real(wp), intent(out) :: dweights(:, :)

      integer :: nbranch, nparam
      integer :: iparam, ibranch
      real(wp) :: inv_sigma, inv_sigma2, mean_q
      real(wp) :: qk(size(phi))

      nbranch = size(phi)
      nparam = size(dphi, dim=1)

      call softmax_weights(phi, sigma_phi, weights)

      if (sigma_phi <= tiny(1.0_wp)) then
         dweights = 0.0_wp
         return
      end if

      inv_sigma = 1.0_wp/sigma_phi
      inv_sigma2 = inv_sigma*inv_sigma

      do iparam = 1, nparam
         do ibranch = 1, nbranch
            qk(ibranch) = -dphi(iparam, ibranch)*inv_sigma &
                          + phi(ibranch)*dsigma_phi(iparam)*inv_sigma2
         end do

         mean_q = dot_product(weights, qk)
         do ibranch = 1, nbranch
            dweights(iparam, ibranch) = weights(ibranch)*(qk(ibranch) - mean_q)
         end do
      end do
   end subroutine softmax_weights_grad

   !> Type-bound wrapper for softmax weights
   !> @param[in]  self    Branch-weight instance
   !> @param[in]  phi     Objective values Phi_m
   !> @param[out] weights Softmax weights p_m
   pure subroutine branch_weight_weights(self, phi, weights)
      class(branch_weight_type), intent(in) :: self
      real(wp), intent(in) :: phi(:)
      real(wp), intent(out) :: weights(:)

      call softmax_weights(phi, self%s, weights)
   end subroutine branch_weight_weights

   !> Type-bound Shannon entropy diagnostic for branch weights
   !>
   !> Computes H = -sum_m p_m log(p_m) with p_m from the branch softmax model.
   !> @param[in]  self    Branch-weight instance
   !> @param[in]  phi     Objective values Phi_m
   !> @param[out] entropy Branch entropy H
   pure subroutine branch_weights_entropy(self, phi, entropy)
      class(branch_weight_type), intent(in) :: self
      real(wp), intent(in) :: phi(:)
      real(wp), intent(out) :: entropy

      real(wp) :: weights(size(phi))
      integer :: ibranch

      entropy = 0.0_wp
      if (size(phi) <= 0) return

      call softmax_weights(phi, self%s, weights)

      do ibranch = 1, size(weights)
         if (weights(ibranch) > tiny(1.0_wp)) then
            entropy = entropy - weights(ibranch)*log(weights(ibranch))
         end if
      end do
   end subroutine branch_weights_entropy

   !> Type-bound wrapper for softmax weights and derivatives
   !>
   !> The softmax width sigma_phi = s is independent of the nuclear
   !> coordinates, so dsigma_phi = 0 and only the dPhi term survives.
   !> @param[in]  self     Branch-weight instance
   !> @param[in]  phi      Objective values Phi_m
   !> @param[in]  dphi     Derivatives dPhi/dq_k
   !> @param[out] weights  Softmax weights p_m
   !> @param[out] dweights Derivatives dp_m/dq_k
   pure subroutine branch_weight_weights_grad(self, phi, dphi, weights, dweights)
      class(branch_weight_type), intent(in) :: self
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: dphi(:, :)
      real(wp), intent(out) :: weights(:)
      real(wp), intent(out) :: dweights(:, :)

      real(wp) :: dsigma_phi(size(dphi, dim=1))

      dsigma_phi = 0.0_wp
      call softmax_weights_grad(phi, dphi, self%s, dsigma_phi, weights, dweights)
   end subroutine branch_weight_weights_grad

end module moist_cavity_drop_branching
