!> Switching functions for DROP
!>
!> All switching functions in this module share an abstract base type [[moist_cavity_drop_swif]]
!> that provides chain-rule helpers for propagating derivatives through the switching function
!>
!> Implemented functions:
!>
!> * [[moist_cavity_drop_sigmoid_bump_swif]] - partition-of-unity bump
!>
!>     Xlo = min(from, to); Xhi = max(from, to); W = Xhi - Xlo
!>     u(X) = 1                                        if X <= Xlo
!>          = 0                                        if X >= Xhi
!>          = (Xhi - X) / W                            otherwise
!>     A(u) = exp(-aHi / u^pHi)                        for 0 < u < 1
!>     B(u) = exp(-aLo / (1-u)^pLo)                    for 0 < u < 1
!>     S_up(X) = B / (A + B)                           for 0 < u < 1
!>             = 1 at u = 0 (X = Xhi), 0 at u = 1 (X = Xlo)
!>
!>   Numerically S_up is evaluated via the equivalent stable sigmoid form
!>     S_up = 1 / (1 + exp(g)),  g = -aHi*u^(-pHi) + aLo*(1-u)^(-pLo)
!>
!> * [[moist_cavity_drop_smooth_step_swif]] - erf-based bump (not used)
!>
!>     lo = min(from, to); hi = max(from, to)
!>     B(x) = 0                                        if x <= lo
!>          = 1                                        if x >= hi
!>          = 0.5 * (1 + erf(arg(x)))                 otherwise
!>     arg(x) = k * (hi - lo) * (2*x - lo - hi) / ((x - lo) * (hi - x))
!>
!> The actual value follows the `from`/`to` orientation:
!>   from < to -> rising 0 -> 1; from > to -> falling 1 -> 0.
module moist_cavity_drop_switching
   use mctc_env_accuracy, only: wp
   use mctc_io_constants, only: pi

   implicit none
   private

   public :: moist_cavity_drop_swif
   public :: moist_cavity_drop_smooth_step_swif, new_smooth_step_swif
   public :: moist_cavity_drop_sigmoid_bump_swif, new_sigmoid_bump_swif

   integer, parameter :: ndim = 3
   real(wp), parameter :: inv_sqrtpi = 1.0_wp/sqrt(pi)
   !> Saturate erf(arg) and exp(-arg^2) once the derivatives are numerically zero
   real(wp), parameter :: arg_cutoff = 26.0_wp
   !> Saturate the sigmoid bump once |g| exceeds this
   real(wp), parameter :: g_cutoff = 50.0_wp

   !> Abstract base type for scalar switching functions `S: R -> [0, 1]`
   type, abstract :: moist_cavity_drop_swif
   contains
      procedure(swif_eval_iface), deferred :: eval
      procedure :: f0 => swif_f0
      procedure :: f1_rA => swif_f1_rA
      procedure :: f2_rArB => swif_f2_rArB
   end type moist_cavity_drop_swif

   abstract interface
      !> Evaluate a scalar switching function and its scalar derivatives.
      !> @param[in]  self  Switching-function instance
      !> @param[in]  x0    Input value
      !> @param[out] f0    Switching value
      !> @param[out] f1    First derivative with respect to x0
      !> @param[out] f2    Second derivative with respect to x0
      pure subroutine swif_eval_iface(self, x0, f0, f1, f2)
         import :: wp, moist_cavity_drop_swif
         class(moist_cavity_drop_swif), intent(in) :: self
         !> Input value
         real(wp), intent(in) :: x0
         !> Switching value
         real(wp), intent(out) :: f0
         !> First derivative with respect to x0
         real(wp), intent(out), optional :: f1
         !> Second derivative with respect to x0
         real(wp), intent(out), optional :: f2
      end subroutine swif_eval_iface
   end interface

   !> Erf-based bump-function switching function
   !>
   !> Transitions from 0 to 1 on the interval [from, to] using
   !> a compact-support bump construction
   type, extends(moist_cavity_drop_swif) :: moist_cavity_drop_smooth_step_swif
      !> Start of transition region (0 before)
      real(wp) :: from = 0.0_wp
      !> End of transition region (1 after)
      real(wp) :: to = 1.0_wp
      !> Dimensionless bump steepness parameter
      real(wp) :: k = 1.0_wp
   contains
      procedure :: eval => smooth_step_eval
   end type moist_cavity_drop_smooth_step_swif

   !> Partition-of-unity bump switching function
   !>
   !> Transitions from 0 to 1 on the interval [from, to] using
   !> mutually-vanishing exponentials `A(u) = exp(-aHi/u^pHi)` and
   !> `B(u) = exp(-aLo/(1-u)^pLo)` glued as `B/(A+B)`
   type, extends(moist_cavity_drop_swif) :: moist_cavity_drop_sigmoid_bump_swif
      !> Start of transition region (0 before)
      real(wp) :: from = 0.0_wp
      !> End of transition region (1 after)
      real(wp) :: to = 1.0_wp
      !> Exponent controlling sharpness of the upper plateau entry
      real(wp) :: p_hi = 2.0_wp
      !> Prefactor controlling sharpness of the upper plateau entry
      real(wp) :: a_hi = 1.0_wp
      !> Exponent controlling sharpness of the lower plateau exit
      real(wp) :: p_lo = 2.0_wp
      !> Prefactor controlling sharpness of the lower plateau exit
      real(wp) :: a_lo = 1.0_wp
   contains
      procedure :: eval => sigmoid_bump_eval
   end type moist_cavity_drop_sigmoid_bump_swif

contains

   !> Construct a new smooth step (erf-bump) switching function
   !>
   !> If from < to: transitions 0 -> 1 on [from, to]
   !> If from > to: transitions 1 -> 0 on [to, from]
   !>
   !> @param[inout] self      Smooth step instance to initialise
   !> @param[in]    from_val  Start value of the transition region
   !> @param[in]    to_val    End value of the transition region
   !> @param[in]    k_val     Optional bump steepness parameter
   subroutine new_smooth_step_swif(self, from_val, to_val, k_val)
      class(moist_cavity_drop_smooth_step_swif), intent(inout) :: self
      !> Start of transition region
      real(wp), intent(in) :: from_val
      !> End of transition region
      real(wp), intent(in) :: to_val
      !> Optional bump steepness parameter
      real(wp), intent(in), optional :: k_val

      self%from = from_val
      self%to = to_val
      self%k = 1.0_wp
      if (present(k_val)) self%k = k_val

   end subroutine new_smooth_step_swif

   !> Construct a new sigmoid bump switching function
   !>
   !> If from < to: transitions 0 -> 1 on [from, to]
   !> If from > to: transitions 1 -> 0 on [to, from]
   !>
   !> @param[inout] self      Sigmoid bump instance to initialise
   !> @param[in]    from_val  Start value of the transition region
   !> @param[in]    to_val    End value of the transition region
   !> @param[in]    p_hi      Optional upper plateau exponent
   !> @param[in]    a_hi      Optional upper plateau prefactor
   !> @param[in]    p_lo      Optional lower plateau exponent
   !> @param[in]    a_lo      Optional lower plateau prefactor
   subroutine new_sigmoid_bump_swif(self, from_val, to_val, p_hi, a_hi, p_lo, a_lo)
      class(moist_cavity_drop_sigmoid_bump_swif), intent(inout) :: self
      !> Start of transition region
      real(wp), intent(in) :: from_val
      !> End of transition region
      real(wp), intent(in) :: to_val
      !> Optional upper plateau exponent (default 2)
      real(wp), intent(in), optional :: p_hi
      !> Optional upper plateau prefactor (default 1)
      real(wp), intent(in), optional :: a_hi
      !> Optional lower plateau exponent (default 2)
      real(wp), intent(in), optional :: p_lo
      !> Optional lower plateau prefactor (default 1)
      real(wp), intent(in), optional :: a_lo

      self%from = from_val
      self%to = to_val
      self%p_hi = 2.0_wp
      self%a_hi = 1.0_wp
      self%p_lo = 2.0_wp
      self%a_lo = 1.0_wp
      if (present(p_hi)) self%p_hi = p_hi
      if (present(a_hi)) self%a_hi = a_hi
      if (present(p_lo)) self%p_lo = p_lo
      if (present(a_lo)) self%a_lo = a_lo

   end subroutine new_sigmoid_bump_swif

   !> Evaluate the erf-bump switching value and its first two scalar
   !> derivatives
   !>
   !> @param[in]  self  Smooth step instance
   !> @param[in]  x0    Input value
   !> @param[out] f0    Switching value
   !> @param[out] f1    First derivative with respect to x0
   !> @param[out] f2    Second derivative with respect to x0
   pure subroutine smooth_step_eval(self, x0, f0, f1, f2)
      class(moist_cavity_drop_smooth_step_swif), intent(in) :: self
      !> Input value
      real(wp), intent(in) :: x0
      !> Switching value
      real(wp), intent(out) :: f0
      !> First derivative with respect to x0
      real(wp), intent(out), optional :: f1
      !> Second derivative with respect to x0
      real(wp), intent(out), optional :: f2
      real(wp) :: lo, hi, width, p, q, arg, darg, d2arg, bump, dbump, d2bump
      real(wp) :: direction

      lo = min(self%from, self%to)
      hi = max(self%from, self%to)
      width = hi - lo
      if (width <= epsilon(1.0_wp)) then
         if (x0 <= lo) then
            bump = 0.0_wp
         else
            bump = 1.0_wp
         end if
         dbump = 0.0_wp
         d2bump = 0.0_wp
      else if (x0 <= lo) then
         bump = 0.0_wp
         dbump = 0.0_wp
         d2bump = 0.0_wp
      else if (x0 >= hi) then
         bump = 1.0_wp
         dbump = 0.0_wp
         d2bump = 0.0_wp
      else
         p = x0 - lo
         q = hi - x0
         arg = self%k*width*(2.0_wp*x0 - lo - hi)/(p*q)

         if (arg <= -arg_cutoff) then
            bump = 0.0_wp
            dbump = 0.0_wp
            d2bump = 0.0_wp
         else if (arg >= arg_cutoff) then
            bump = 1.0_wp
            dbump = 0.0_wp
            d2bump = 0.0_wp
         else
            darg = self%k*width*(1.0_wp/(p*p) + 1.0_wp/(q*q))
            d2arg = 2.0_wp*self%k*width*(1.0_wp/(q*q*q) - 1.0_wp/(p*p*p))

            bump = 0.5_wp*(1.0_wp + erf(arg))
            dbump = inv_sqrtpi*exp(-(arg*arg))*darg
            d2bump = inv_sqrtpi*exp(-(arg*arg))*(d2arg - 2.0_wp*arg*darg*darg)
         end if
      end if

      direction = 1.0_wp
      if (self%from > self%to) direction = -1.0_wp

      f0 = merge(bump, 1.0_wp - bump, direction > 0.0_wp)

      if (present(f1)) f1 = direction*dbump
      if (present(f2)) f2 = direction*d2bump

   end subroutine smooth_step_eval

   !> Evaluate the sigmoid bump switching value and its first two scalar
   !> derivatives
   !>
   !> Uses the stable logistic form `S_up = 1/(1 + exp(g))` with
   !> `g(u) = -aHi*u^(-pHi) + aLo*(1-u)^(-pLo)` and `u = (hi - x0)/width`
   !>
   !> @param[in]  self  Sigmoid bump instance
   !> @param[in]  x0    Input value
   !> @param[out] f0    Switching value
   !> @param[out] f1    First derivative with respect to x0
   !> @param[out] f2    Second derivative with respect to x0
   pure subroutine sigmoid_bump_eval(self, x0, f0, f1, f2)
      class(moist_cavity_drop_sigmoid_bump_swif), intent(in) :: self
      !> Input value
      real(wp), intent(in) :: x0
      !> Switching value
      real(wp), intent(out) :: f0
      !> First derivative with respect to x0
      real(wp), intent(out), optional :: f1
      !> Second derivative with respect to x0
      real(wp), intent(out), optional :: f2
      real(wp) :: lo, hi, width, u, one_minus_u
      real(wp) :: p_hi, a_hi, p_lo, a_lo
      real(wp) :: g, gp, gpp, s, s1m, exp_neg_g
      real(wp) :: dS_du, d2S_du, du_dx, d2u_dx2
      real(wp) :: bump, dbump, d2bump
      real(wp) :: direction

      lo = min(self%from, self%to)
      hi = max(self%from, self%to)
      width = hi - lo
      p_hi = self%p_hi
      a_hi = self%a_hi
      p_lo = self%p_lo
      a_lo = self%a_lo

      if (width <= epsilon(1.0_wp)) then
         if (x0 <= lo) then
            bump = 0.0_wp
         else
            bump = 1.0_wp
         end if
         dbump = 0.0_wp
         d2bump = 0.0_wp
      else if (x0 <= lo) then
         bump = 0.0_wp
         dbump = 0.0_wp
         d2bump = 0.0_wp
      else if (x0 >= hi) then
         bump = 1.0_wp
         dbump = 0.0_wp
         d2bump = 0.0_wp
      else
         u = (hi - x0)/width
         one_minus_u = 1.0_wp - u
         g = -a_hi*u**(-p_hi) + a_lo*one_minus_u**(-p_lo)

         if (g >= g_cutoff) then
            bump = 0.0_wp
            dbump = 0.0_wp
            d2bump = 0.0_wp
         else if (g <= -g_cutoff) then
            bump = 1.0_wp
            dbump = 0.0_wp
            d2bump = 0.0_wp
         else
            !> Stable sigmoid evaluation for any sign of g.
            if (g >= 0.0_wp) then
               exp_neg_g = exp(-g)
               s = exp_neg_g/(1.0_wp + exp_neg_g)
            else
               s = 1.0_wp/(1.0_wp + exp(g))
            end if
            s1m = s*(1.0_wp - s)

            gp = p_hi*a_hi*u**(-p_hi - 1.0_wp) &
                 + p_lo*a_lo*one_minus_u**(-p_lo - 1.0_wp)
            gpp = -p_hi*(p_hi + 1.0_wp)*a_hi*u**(-p_hi - 2.0_wp) &
                  + p_lo*(p_lo + 1.0_wp)*a_lo*one_minus_u**(-p_lo - 2.0_wp)

            dS_du = -s1m*gp
            d2S_du = s1m*((1.0_wp - 2.0_wp*s)*gp*gp - gpp)

            !> u(x) = (hi - x)/width  =>  du/dx = -1/width,  d2u/dx2 = 0
            du_dx = -1.0_wp/width
            d2u_dx2 = 0.0_wp

            bump = s
            dbump = dS_du*du_dx
            d2bump = d2S_du*du_dx*du_dx + dS_du*d2u_dx2
         end if
      end if

      direction = 1.0_wp
      if (self%from > self%to) direction = -1.0_wp

      f0 = merge(bump, 1.0_wp - bump, direction > 0.0_wp)

      if (present(f1)) f1 = direction*dbump
      if (present(f2)) f2 = direction*d2bump

   end subroutine sigmoid_bump_eval

   !> Evaluate the switching function value via the abstract `eval` hook
   !>
   !> When from < to the value rises smoothly from 0 to 1.
   !> When from > to the value falls smoothly from 1 to 0.
   !> @param[in] self  Switching-function instance
   !> @param[in] x0    Input value
   pure function swif_f0(self, x0) result(val)
      class(moist_cavity_drop_swif), intent(in) :: self
      !> Input value
      real(wp), intent(in) :: x0
      real(wp) :: val

      call self%eval(x0, val)

   end function swif_f0

   !> First derivative of the switching function with respect to nuclear
   !> coordinates, via the chain rule through x0
   !>
   !> dS/dR_A = S'(x0) * dx0/dR_A
   !> @param[in] self    Switching-function instance
   !> @param[in] x0      Input scalar value
   !> @param[in] x1      Gradient of x0 w.r.t. nuclear coordinates (3, ncenters)
   !> @param[in] active  Optional list of active atom indices
   pure function swif_f1_rA(self, x0, x1, active) result(grad)
      class(moist_cavity_drop_swif), intent(in) :: self
      !> Input scalar value
      real(wp), intent(in) :: x0
      !> Gradient of x0 w.r.t. nuclear positions (3, ncenters)
      real(wp), intent(in) :: x1(:, :)
      !> Optional active atom indices for screening
      integer, intent(in), optional :: active(:)
      real(wp) :: grad(size(x1, 1), size(x1, 2))
      real(wp) :: f0_dummy, dS
      integer :: ii, iatom

      call self%eval(x0, f0_dummy, dS)

      if (present(active)) then
         grad = 0.0_wp
         do ii = 1, size(active)
            iatom = active(ii)
            grad(:, iatom) = dS*x1(:, iatom)
         end do
      else
         grad = dS*x1
      end if

   end function swif_f1_rA

   !> Second derivative (Hessian) of the switching function with respect
   !> to nuclear coordinates
   !>
   !> d^2S/(dR_A dR_B) = S''(x0) * dx0/dR_A * dx0/dR_B
   !>                   + S'(x0) * d^2x0/(dR_A dR_B)
   !> @param[in] self  Switching-function instance
   !> @param[in] x0    Input scalar value
   !> @param[in] x1    Gradient of x0 (3, ncenters)
   !> @param[in] x2    Hessian of x0 (3, 3, ncenters, ncenters)
   pure function swif_f2_rArB(self, x0, x1, x2) result(hess)
      class(moist_cavity_drop_swif), intent(in) :: self
      !> Input scalar value
      real(wp), intent(in) :: x0
      !> Gradient of x0 w.r.t. nuclear positions
      real(wp), intent(in) :: x1(:, :)
      !> Hessian of x0 w.r.t. nuclear positions
      real(wp), intent(in) :: x2(:, :, :, :)
      real(wp) :: hess(size(x2, 1), size(x2, 2), size(x2, 3), size(x2, 4))
      real(wp) :: f0_dummy, dS, d2S
      integer :: ncenters, ai, bi, i, j

      ncenters = size(x1, 2)

      hess = 0.0_wp
      call self%eval(x0, f0_dummy, dS, d2S)

      do j = 1, ncenters
         do i = 1, ncenters
            do bi = 1, ndim
               do ai = 1, ndim
                  hess(ai, bi, i, j) = d2S*x1(ai, i)*x1(bi, j) &
                     & + dS*x2(ai, bi, i, j)
               end do
            end do
         end do
      end do

   end function swif_f2_rArB

end module moist_cavity_drop_switching
