!> Branch weights of the DROP projection
!>
!> An anchor whose closest-point problem has several solutions (branches) shares
!> its quadrature weight among them. Every branch carries a distance coordinate
!> `x_i` and a gate `g_i` in `[0, 1]`, the product of its own switches; a closer
!> branch takes weight from a farther one in proportion to its own gate:
!>
!>     x_i    = sqrt(pdist_i^2 + eps^2),  pdist_i^2 = 2 Phi_i/alpha,  eps = x_max/2
!>
!>     h(d)   = 1                                  d <= 0
!>            = eta(1 - t)/(eta(t) + eta(1 - t))   0 < t = d/x_max < 1
!>            = 0                                  d >= x_max
!>     eta(t) = exp(-1/t)
!>
!>     m_ij   = g_j (1 - h(x_i - x_j))
!>     a_i    = prod_{j /= i} (1 - m_ij)
!>     Z      = sum_k g_k a_k + prod_k (1 - g_k)
!>     w_i    = g_i a_i / Z
!>
!> - `w_i` is the effective weight of the branch, its own gate included
!> - a branch with `g_j = 0` takes nothing: `m_ij = 0` and its term in `Z`
!>   vanishes, so the other weights are those of the group without it; the
!>   gates being flat at zero, their first derivatives are as well
!> - a branch `x_max` or more behind a fully open one (`g_j = 1`) has `a_i = 0`:
!>   weight zero, and no other weight depends on it (`h = 1` for every closer
!>   branch); a search may stop at that distance
!> - behind a partly open branch the far one keeps the factor `1 - g_j`; with
!>   every pair separated by `x_max` or more `Z = 1` and the weights are those
!>   of front-to-back compositing, `w_i = g_i prod_{j closer} (1 - g_j)`
!> - `eps` keeps the coordinate smooth where an anchor crosses the surface
!>   (`pdist = 0`, where the distance itself has a cusp); the window in `pdist`
!>   then widens, to `sqrt(2) x_max` at `pdist = 0`
!> - `Z >= 1 >= a_i`, so `w_i <= g_i`
!> - a single branch has `w = g`; open branches (`g = 1`) share as
!>   `m = 1 - h`, equally at one distance
module moist_cavity_drop_branching
   use mctc_env_accuracy, only: wp
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   implicit none(type, external)
   private

   public :: branch_weight_type, branch_pdist2_max, branch_soft_fraction

   !> Softening length of the distance coordinate, in units of the window
   real(wp), parameter :: branch_soft_fraction = 0.5_wp

   !> Branch-weight model for degenerate closest-point branches
   type :: branch_weight_type
      !> Distance excess at which a branch loses all weight (Bohr)
      real(wp) :: x_max = 0.0_wp
      !> Softening length of the distance coordinate (Bohr)
      real(wp) :: soft = 0.0_wp
      !> Objective weight, `Phi = (alpha/2) pdist^2`
      real(wp) :: alpha = 1.0_wp
      !> Separation below which two points are one branch (Bohr)
      real(wp) :: sep_cut = 1.0e-8_wp
      !> Pair factor `h`, falling from 1 to 0 on `[0, x_max]`
      type(moist_math_smoothstep_sigmoid_type), private :: h
   contains
      !> Set the window
      procedure :: init => branch_weight_init
      !> Branch weights of one group
      procedure :: weights => branch_weight_weights
      !> Branch weights and their Jacobians in the objectives and gates
      procedure :: weights_jac => branch_weight_weights_jac
      !> Branch weights and their derivatives along given input derivatives
      procedure :: weights_grad => branch_weight_weights_grad
      !> Largest distance with weight behind a fully open branch
      procedure :: pdist_max => branch_weight_pdist_max
   end type branch_weight_type

contains

   !> Initialize the branch-weight model
   !>
   !> A non-positive window lets a closer branch take its full gate from every
   !> farther one; exact ties share
   !>
   !> @param[inout] self  Branch-weight instance
   !> @param[in]    x_max Distance excess at which a branch loses all weight (Bohr)
   !> @param[in]    alpha Objective weight, `Phi = (alpha/2) pdist^2`
   !> @param[in]    sep_cut Separation below which two points are one branch (Bohr);
   !>                       kept when absent
   subroutine branch_weight_init(self, x_max, alpha, sep_cut)
      class(branch_weight_type), intent(inout) :: self
      real(wp), intent(in) :: x_max
      real(wp), intent(in) :: alpha
      real(wp), intent(in), optional :: sep_cut

      self%x_max = max(x_max, 0.0_wp)
      self%soft = branch_soft_fraction*self%x_max
      self%alpha = alpha
      if (present(sep_cut)) self%sep_cut = sep_cut
      call new_smoothstep_sigmoid(self%h, x_off=self%x_max, x_on=0.0_wp, &
                                  p_hi=1.0_wp, a_hi=1.0_wp, p_lo=1.0_wp, a_lo=1.0_wp)
   end subroutine branch_weight_init

   !> Distance coordinate of an objective value
   !>
   !> @param[in] self Branch-weight instance
   !> @param[in] phi  Objective
   !> @return    x    `sqrt(2 phi/alpha + eps^2)`
   elemental function coordinate(self, phi) result(x)
      type(branch_weight_type), intent(in) :: self
      real(wp), intent(in) :: phi
      real(wp) :: x

      x = sqrt(2.0_wp*max(phi, 0.0_wp)/self%alpha + self%soft**2)
   end function coordinate

   !> Largest distance with weight behind a fully open branch at `pdist_min`
   !>
   !> A branch has an exactly vanishing weight once its coordinate exceeds the
   !> one of a fully open branch by `x_max`, whatever its other weights:
   !>
   !>     pdist_max = sqrt( (sqrt(pdist_min^2 + soft^2) + x_max)^2 - soft^2 )
   !>
   !> A shell `x_max` thick far from the surface, `sqrt(2)` times that at
   !> `pdist_min = 0`
   !>
   !> @param[in] self      Branch-weight instance
   !> @param[in] pdist_min Projection distance of the closest fully open branch (Bohr)
   !> @return    pdist_max Largest distance inside its window (Bohr)
   elemental function branch_weight_pdist_max(self, pdist_min) result(pdist_max)
      class(branch_weight_type), intent(in) :: self
      real(wp), intent(in) :: pdist_min
      real(wp) :: pdist_max

      pdist_max = sqrt(branch_pdist2_max(pdist_min*pdist_min, self%x_max**2, self%x_max, self%soft))
   end function branch_weight_pdist_max

   !> Largest squared distance inside the branch window of a closest branch
   !>
   !> - `pdist_slack` is the distance window `x_max` of the branch weight and
   !>   `pdist_soft` its softening; a branch beyond the returned value has no
   !>   weight behind the branch at `pdist2_min`
   !> - `pdist2_slack` is a constant squared-distance part; the weight's window
   !>   has `pdist2_slack = pdist_slack^2`
   !>
   !> @param[in] pdist2_min   Squared distance of the closest branch (Bohr^2)
   !> @param[in] pdist2_slack Constant squared-distance window (Bohr^2)
   !> @param[in] pdist_slack  Distance window (Bohr)
   !> @param[in] pdist_soft   Softening length (Bohr)
   !> @return    pdist2_max   Squared distance bounding the window (Bohr^2)
   elemental function branch_pdist2_max(pdist2_min, pdist2_slack, pdist_slack, pdist_soft) result(pdist2_max)
      real(wp), intent(in) :: pdist2_min
      real(wp), intent(in) :: pdist2_slack
      real(wp), intent(in) :: pdist_slack
      real(wp), intent(in) :: pdist_soft
      real(wp) :: pdist2_max

      pdist2_max = pdist2_min + pdist2_slack + 2.0_wp*pdist_slack*sqrt(pdist2_min + pdist_soft**2)
   end function branch_pdist2_max

   !> Branch weights of one group
   !>
   !> @param[in]  self    Branch-weight instance
   !> @param[in]  phi     Objective per branch
   !> @param[in]  gate    Gate per branch, in `[0, 1]`
   !> @param[out] weights Branch weights, the gate included
   pure subroutine branch_weight_weights(self, phi, gate, weights)
      class(branch_weight_type), intent(in) :: self
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: gate(:)
      real(wp), intent(out) :: weights(:)

      integer :: n, i, j
      real(wp) :: hij, x(size(phi)), a(size(phi))

      n = size(phi)
      if (n <= 0) return
      x = coordinate(self, phi)
      do i = 1, n
         a(i) = 1.0_wp
         do j = 1, n
            if (j == i) cycle
            call pair_factor(self, x(i) - x(j), hij)
            a(i) = a(i)*(1.0_wp - gate(j)*(1.0_wp - hij))
         end do
      end do
      weights(1:n) = gate(1:n)*a/(sum(gate(1:n)*a) + product(1.0_wp - gate(1:n)))
   end subroutine branch_weight_weights

   !> Pair factor `h` of a coordinate difference and its derivative
   !>
   !> @param[in]  self Branch-weight instance
   !> @param[in]  d    Coordinate difference `x_i - x_j`
   !> @param[out] h    Pair factor
   !> @param[out] hp   Its derivative in `d`
   pure subroutine pair_factor(self, d, h, hp)
      type(branch_weight_type), intent(in) :: self
      real(wp), intent(in) :: d
      real(wp), intent(out) :: h
      real(wp), intent(out), optional :: hp

      if (self%x_max > 0.0_wp) then
         call self%h%eval(d, h, hp)
      else
         ! Zero window: a step, the closer branch takes everything
         h = merge(1.0_wp, 0.0_wp, d <= 0.0_wp)
         if (present(hp)) hp = 0.0_wp
      end if
   end subroutine pair_factor

   !> Branch weights and their Jacobians in the objectives and gates
   !>
   !> Products are formed directly and never divided by a factor, so a vanishing
   !> factor gives an exactly vanishing weight and derivative
   !>
   !> @param[in]  self     Branch-weight instance
   !> @param[in]  phi      Objective per branch (n)
   !> @param[in]  gate     Gate per branch (n)
   !> @param[out] weights  Branch weights (n)
   !> @param[out] jac_phi  `jac_phi(i, l) = d weights(i) / d phi(l)` (n, n)
   !> @param[out] jac_gate `jac_gate(i, l) = d weights(i) / d gate(l)` (n, n)
   pure subroutine branch_weight_weights_jac(self, phi, gate, weights, jac_phi, jac_gate)
      class(branch_weight_type), intent(in) :: self
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: gate(:)
      real(wp), intent(out) :: weights(:)
      real(wp), intent(out) :: jac_phi(:, :)
      real(wp), intent(out) :: jac_gate(:, :)

      integer :: n, i, k, l
      ! `1 - m(i, j)` and the partials of `m(i, j)` in `g_j` and in `x_j`, with `d = x_i - x_j`
      real(wp) :: om(size(phi), size(phi)), dm_g(size(phi), size(phi)), dm_x(size(phi), size(phi))
      ! `a_i` and its partials in the coordinates and the gates
      real(wp) :: a(size(phi)), da_x(size(phi), size(phi)), da_g(size(phi), size(phi))
      ! `u_k = g_k a_k` and `prod_{q /= k} (1 - g_q)`
      real(wp) :: u(size(phi)), others(size(phi))
      real(wp) :: rest, total, dtotal, h, hp
      ! Distance coordinates and their derivatives in the objectives
      real(wp) :: x(size(phi)), dx(size(phi))

      n = size(phi)
      if (n <= 0) return

      x = coordinate(self, phi)
      ! Zero only without a window, where no weight depends on the coordinates
      dx = 0.0_wp
      where (x > 0.0_wp) dx = 1.0_wp/(self%alpha*x)

      do l = 1, n
         do i = 1, n
            if (i == l) then
               om(i, l) = 1.0_wp
               dm_g(i, l) = 0.0_wp
               dm_x(i, l) = 0.0_wp
            else
               call pair_factor(self, x(i) - x(l), h, hp)
               om(i, l) = 1.0_wp - gate(l)*(1.0_wp - h)
               dm_g(i, l) = 1.0_wp - h
               ! dm/dx_l = -g_l h' d(x_i - x_l)/dx_l
               dm_x(i, l) = gate(l)*hp
            end if
         end do
      end do

      do i = 1, n
         a(i) = product(om(i, 1:n))
         da_x(i, i) = 0.0_wp
         da_g(i, i) = 0.0_wp
         do l = 1, n
            if (l == i) cycle
            ! Product over the factors of a_i other than the one with l
            rest = 1.0_wp
            do k = 1, n
               if (k == i .or. k == l) cycle
               rest = rest*om(i, k)
            end do
            da_g(i, l) = -dm_g(i, l)*rest
            da_x(i, l) = -dm_x(i, l)*rest
            ! x_i enters every pair with the opposite sign
            da_x(i, i) = da_x(i, i) + dm_x(i, l)*rest
         end do
      end do

      u = gate(1:n)*a
      do k = 1, n
         others(k) = 1.0_wp
         do l = 1, n
            if (l /= k) others(k) = others(k)*(1.0_wp - gate(l))
         end do
      end do
      total = sum(u) + product(1.0_wp - gate(1:n))
      weights(1:n) = u/total

      ! dZ = sum_k (du_k - others_k dg_k), du_k = g_k da_k + a_k dg_k
      do l = 1, n
         dtotal = 0.0_wp
         do k = 1, n
            dtotal = dtotal + gate(k)*da_x(k, l)
         end do
         do i = 1, n
            jac_phi(i, l) = (gate(i)*da_x(i, l) - weights(i)*dtotal)/total*dx(l)
         end do
         dtotal = a(l) - others(l)
         do k = 1, n
            dtotal = dtotal + gate(k)*da_g(k, l)
         end do
         do i = 1, n
            jac_gate(i, l) = (gate(i)*da_g(i, l) - weights(i)*dtotal)/total
         end do
         jac_gate(l, l) = jac_gate(l, l) + a(l)/total
      end do
   end subroutine branch_weight_weights_jac

   !> Branch weights and their derivatives along given input derivatives
   !>
   !> @param[in]  self     Branch-weight instance
   !> @param[in]  phi      Objective per branch (nbranch)
   !> @param[in]  gate     Gate per branch (nbranch)
   !> @param[in]  dphi     Objective derivatives (nparam, nbranch)
   !> @param[in]  dgate    Gate derivatives (nparam, nbranch)
   !> @param[out] weights  Branch weights (nbranch)
   !> @param[out] dweights Weight derivatives (nparam, nbranch)
   pure subroutine branch_weight_weights_grad(self, phi, gate, dphi, dgate, weights, dweights)
      class(branch_weight_type), intent(in) :: self
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: gate(:)
      real(wp), intent(in) :: dphi(:, :)
      real(wp), intent(in) :: dgate(:, :)
      real(wp), intent(out) :: weights(:)
      real(wp), intent(out) :: dweights(:, :)

      integer :: n, i, l
      real(wp) :: jac_phi(size(phi), size(phi)), jac_gate(size(phi), size(phi))

      n = size(phi)
      call self%weights_jac(phi, gate, weights, jac_phi, jac_gate)
      do i = 1, n
         dweights(:, i) = 0.0_wp
         do l = 1, n
            if (jac_phi(i, l) /= 0.0_wp) dweights(:, i) = dweights(:, i) + jac_phi(i, l)*dphi(:, l)
            if (jac_gate(i, l) /= 0.0_wp) dweights(:, i) = dweights(:, i) + jac_gate(i, l)*dgate(:, l)
         end do
      end do
   end subroutine branch_weight_weights_grad

end module moist_cavity_drop_branching
