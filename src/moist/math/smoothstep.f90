!> Smoothstep (switching) functions `S: R -> [0, 1]` on a band `[from, to]`
!>
!> All smoothsteps share the abstract base type [[moist_math_smoothstep_type]], which
!> provides chain-rule helpers for propagating derivatives through the smoothstep
!>
!> Orientation follows `from`/`to`:
!>   from < to -> rising 0 -> 1 on [from, to]; from > to -> falling 1 -> 0 on [to, from]
!>   exactly 0/1 outside the band, all derivatives exactly 0 there
!>
!> Implemented functions:
!>
!> * [[moist_math_smoothstep_polynomial_type]] - polynomial smoothstep, `C^n` at both ends
!>
!>     lo = min(from, to); hi = max(from, to); t = clamp((x - lo)/(hi - lo), 0, 1)
!>     S_n(t) = t^(n+1) sum_{k=0}^{n} C(n+k, k) C(2n+1, n-k) (-t)^k,  degree 2n+1
!>     coefficients of t^(n+1), ..., t^(2n+1):
!>       S_1: 3, -2
!>       S_2: 10, -15, 6
!>       S_3: 35, -84, 70, -20
!>       S_4: 126, -420, 540, -315, 70
!>       S_5: 462, -1980, 3465, -3080, 1386, -252
!>
!> * [[moist_math_smoothstep_sigmoid_type]] - partition-of-unity bump
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
!> * [[moist_math_smoothstep_erf_type]] - erf-based bump
!>
!>     lo = min(from, to); hi = max(from, to)
!>     B(x) = 0                                        if x <= lo
!>          = 1                                        if x >= hi
!>          = 0.5 * (1 + erf(arg(x)))                 otherwise
!>     arg(x) = k * (hi - lo) * (2*x - lo - hi) / ((x - lo) * (hi - x))
!>
!> Derivatives: polynomial through order 4; sigmoid and erf through order 2
module moist_math_smoothstep
   use mctc_env_accuracy, only: wp
   use mctc_io_constants, only: pi
   use, intrinsic :: iso_fortran_env, only: int64

   implicit none(type, external)
   private

   public :: moist_math_smoothstep_type
   public :: moist_math_smoothstep_erf_type, new_smoothstep_erf
   public :: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   public :: moist_math_smoothstep_polynomial_type, new_smoothstep_polynomial

   integer, parameter :: ndim = 3
   real(wp), parameter :: inv_sqrtpi = 1.0_wp/sqrt(pi)
   !> Saturate erf(arg) and exp(-arg^2) once the derivatives are numerically zero
   real(wp), parameter :: arg_cutoff = 26.0_wp
   !> Saturate the sigmoid bump once |g| exceeds this
   real(wp), parameter :: g_cutoff = 50.0_wp
   !> Highest polynomial smoothstep order n
   integer, parameter :: max_order = 5
   !> Highest polynomial degree, 2*max_order + 1
   integer, parameter :: max_degree = 2*max_order + 1
   !> Highest derivative order of the polynomial smoothstep
   integer, parameter :: max_deriv = 4

   !> Abstract base type for scalar smoothsteps `S: R -> [0, 1]`
   type, abstract :: moist_math_smoothstep_type
   contains
      procedure(smoothstep_eval_iface), deferred :: eval
      procedure :: f0 => smoothstep_f0
      procedure :: f1_rA => smoothstep_f1_rA
      procedure :: f2_rArB => smoothstep_f2_rArB
   end type moist_math_smoothstep_type

   abstract interface
      !> Evaluate a scalar smoothstep and its scalar derivatives
      !>
      !> - Types without a third/fourth derivative `error stop` when `f3` or `f4` is present
      !>
      !> @param[in]  self  Smoothstep instance
      !> @param[in]  x0    Input value
      !> @param[out] f0    Smoothstep value
      !> @param[out] f1    First derivative with respect to x0
      !> @param[out] f2    Second derivative with respect to x0
      !> @param[out] f3    Third derivative with respect to x0
      !> @param[out] f4    Fourth derivative with respect to x0
      pure subroutine smoothstep_eval_iface(self, x0, f0, f1, f2, f3, f4)
         import :: wp, moist_math_smoothstep_type
         implicit none(type, external)
         class(moist_math_smoothstep_type), intent(in) :: self
         !> Input value
         real(wp), intent(in) :: x0
         !> Smoothstep value
         real(wp), intent(out) :: f0
         !> First derivative with respect to x0
         real(wp), intent(out), optional :: f1
         !> Second derivative with respect to x0
         real(wp), intent(out), optional :: f2
         !> Third derivative with respect to x0
         real(wp), intent(out), optional :: f3
         !> Fourth derivative with respect to x0
         real(wp), intent(out), optional :: f4
      end subroutine smoothstep_eval_iface
   end interface

   !> Erf-based bump smoothstep
   !>
   !> Transitions from 0 to 1 on the interval [from, to] using
   !> a compact-support bump construction
   !>
   !> - No production consumer, exercised by the unit tests only
   type, extends(moist_math_smoothstep_type) :: moist_math_smoothstep_erf_type
      !> Start of transition region (0 before)
      real(wp) :: from = 0.0_wp
      !> End of transition region (1 after)
      real(wp) :: to = 1.0_wp
      !> Dimensionless bump steepness parameter
      real(wp) :: k = 1.0_wp
   contains
      procedure :: eval => smoothstep_erf_eval
   end type moist_math_smoothstep_erf_type

   !> Partition-of-unity bump smoothstep
   !>
   !> Transitions from 0 to 1 on the interval [from, to] using
   !> mutually-vanishing exponentials `A(u) = exp(-aHi/u^pHi)` and
   !> `B(u) = exp(-aLo/(1-u)^pLo)` glued as `B/(A+B)`
   type, extends(moist_math_smoothstep_type) :: moist_math_smoothstep_sigmoid_type
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
      procedure :: eval => smoothstep_sigmoid_eval
   end type moist_math_smoothstep_sigmoid_type

   !> Polynomial smoothstep `S_n` of degree `2n+1`, `C^n` at both band edges
   !>
   !> - Integer coefficients; plateau values exactly 0/1 and derivatives exactly 0
   !> - Derivatives through order 4
   type, extends(moist_math_smoothstep_type) :: moist_math_smoothstep_polynomial_type
      !> Start of transition region (0 before)
      real(wp) :: from = 0.0_wp
      !> End of transition region (1 after)
      real(wp) :: to = 1.0_wp
      !> Continuity order n, 1..5; 0 until constructed
      integer :: order = 0
      !> Polynomial degree 2n+1; -1 until constructed
      integer :: degree = -1
      !> `coef(j, m)` multiplies `t^j` in `d^m S_n/dt^m`
      real(wp) :: coef(0:max_degree, 0:max_deriv) = 0.0_wp
   contains
      procedure :: eval => smoothstep_polynomial_eval
   end type moist_math_smoothstep_polynomial_type

contains

   !* ================================================================================= *!
   !*                                   Constructors                                    *!
   !* ================================================================================= *!

   !> Construct a new erf-bump smoothstep
   !>
   !> If from < to: transitions 0 -> 1 on [from, to]
   !> If from > to: transitions 1 -> 0 on [to, from]
   !>
   !> @param[inout] self      Smoothstep instance to initialise
   !> @param[in]    from_val  Start value of the transition region
   !> @param[in]    to_val    End value of the transition region
   !> @param[in]    k_val     Optional bump steepness parameter
   subroutine new_smoothstep_erf(self, from_val, to_val, k_val)
      class(moist_math_smoothstep_erf_type), intent(inout) :: self
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

   end subroutine new_smoothstep_erf

   !> Construct a new sigmoid bump smoothstep
   !>
   !> If from < to: transitions 0 -> 1 on [from, to]
   !> If from > to: transitions 1 -> 0 on [to, from]
   !>
   !> @param[inout] self      Smoothstep instance to initialise
   !> @param[in]    from_val  Start value of the transition region
   !> @param[in]    to_val    End value of the transition region
   !> @param[in]    p_hi      Optional upper plateau exponent
   !> @param[in]    a_hi      Optional upper plateau prefactor
   !> @param[in]    p_lo      Optional lower plateau exponent
   !> @param[in]    a_lo      Optional lower plateau prefactor
   subroutine new_smoothstep_sigmoid(self, from_val, to_val, p_hi, a_hi, p_lo, a_lo)
      class(moist_math_smoothstep_sigmoid_type), intent(inout) :: self
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

   end subroutine new_smoothstep_sigmoid

   !> Construct a new polynomial smoothstep `S_n`
   !>
   !> If from < to: transitions 0 -> 1 on [from, to]
   !> If from > to: transitions 1 -> 0 on [to, from]
   !>
   !> - Coefficients `(-1)^k C(n+k, k) C(2n+1, n-k)` of `t^(n+1+k)`, k = 0..n, in integers
   !> - t-derivative tables by exact integer differentiation
   !> - `order` outside 1..5 is an `error stop`
   !>
   !> @param[inout] self      Smoothstep instance to initialise
   !> @param[in]    from_val  Start value of the transition region
   !> @param[in]    to_val    End value of the transition region
   !> @param[in]    order     Optional continuity order n, 1..5 (default 4)
   subroutine new_smoothstep_polynomial(self, from_val, to_val, order)
      class(moist_math_smoothstep_polynomial_type), intent(inout) :: self
      !> Start of transition region
      real(wp), intent(in) :: from_val
      !> End of transition region
      real(wp), intent(in) :: to_val
      !> Optional continuity order n (default 4)
      integer, intent(in), optional :: order
      integer(int64) :: binom(0:max_degree, 0:max_degree)
      integer(int64) :: icoef(0:max_degree)
      integer :: n, degree, j, k, m

      n = 4
      if (present(order)) n = order
      if (n < 1 .or. n > max_order) error stop "moist smoothstep: polynomial order must be in 1..5"
      degree = 2*n + 1

      ! Pascal triangle, binom(j, k) = C(j, k)
      binom = 0_int64
      binom(0, 0) = 1_int64
      do j = 1, max_degree
         binom(j, 0) = 1_int64
         do k = 1, j
            binom(j, k) = binom(j - 1, k - 1) + binom(j - 1, k)
         end do
      end do

      icoef = 0_int64
      do k = 0, n
         icoef(n + 1 + k) = merge(1_int64, -1_int64, mod(k, 2) == 0)*binom(n + k, k)*binom(degree, n - k)
      end do

      self%from = from_val
      self%to = to_val
      self%order = n
      self%degree = degree
      self%coef = 0.0_wp
      self%coef(:, 0) = real(icoef, wp)
      do m = 1, max_deriv
         ! Ascending j reads icoef(j + 1) before it is overwritten
         do j = 0, max_degree - 1
            icoef(j) = icoef(j + 1)*int(j + 1, int64)
         end do
         icoef(max_degree) = 0_int64
         self%coef(:, m) = real(icoef, wp)
      end do

   end subroutine new_smoothstep_polynomial

   !* ================================================================================= *!
   !*                                    Evaluators                                     *!
   !* ================================================================================= *!

   !> Evaluate the erf-bump smoothstep value and its first two scalar
   !> derivatives
   !>
   !> @param[in]  self  Smoothstep instance
   !> @param[in]  x0    Input value
   !> @param[out] f0    Smoothstep value
   !> @param[out] f1    First derivative with respect to x0
   !> @param[out] f2    Second derivative with respect to x0
   !> @param[out] f3    Not available, `error stop` when present
   !> @param[out] f4    Not available, `error stop` when present
   pure subroutine smoothstep_erf_eval(self, x0, f0, f1, f2, f3, f4)
      class(moist_math_smoothstep_erf_type), intent(in) :: self
      !> Input value
      real(wp), intent(in) :: x0
      !> Smoothstep value
      real(wp), intent(out) :: f0
      !> First derivative with respect to x0
      real(wp), intent(out), optional :: f1
      !> Second derivative with respect to x0
      real(wp), intent(out), optional :: f2
      !> Third derivative with respect to x0, not available
      real(wp), intent(out), optional :: f3
      !> Fourth derivative with respect to x0, not available
      real(wp), intent(out), optional :: f4
      real(wp) :: lo, hi, width, p, q, arg, darg, d2arg, bump, dbump, d2bump
      real(wp) :: direction

      if (present(f3) .or. present(f4)) error stop "moist smoothstep: erf type has no third/fourth derivative"

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

            bump = 0.5_wp*erfc(-arg)
            dbump = inv_sqrtpi*exp(-(arg*arg))*darg
            d2bump = inv_sqrtpi*exp(-(arg*arg))*(d2arg - 2.0_wp*arg*darg*darg)
         end if
      end if

      direction = 1.0_wp
      if (self%from > self%to) direction = -1.0_wp

      f0 = merge(bump, 1.0_wp - bump, direction > 0.0_wp)

      if (present(f1)) f1 = direction*dbump
      if (present(f2)) f2 = direction*d2bump

   end subroutine smoothstep_erf_eval

   !> Evaluate the sigmoid bump smoothstep value and its first two scalar
   !> derivatives
   !>
   !> Uses the stable logistic form `S_up = 1/(1 + exp(g))` with
   !> `g(u) = -aHi*u^(-pHi) + aLo*(1-u)^(-pLo)` and `u = (hi - x0)/width`
   !>
   !> @param[in]  self  Smoothstep instance
   !> @param[in]  x0    Input value
   !> @param[out] f0    Smoothstep value
   !> @param[out] f1    First derivative with respect to x0
   !> @param[out] f2    Second derivative with respect to x0
   !> @param[out] f3    Not available, `error stop` when present
   !> @param[out] f4    Not available, `error stop` when present
   pure subroutine smoothstep_sigmoid_eval(self, x0, f0, f1, f2, f3, f4)
      class(moist_math_smoothstep_sigmoid_type), intent(in) :: self
      !> Input value
      real(wp), intent(in) :: x0
      !> Smoothstep value
      real(wp), intent(out) :: f0
      !> First derivative with respect to x0
      real(wp), intent(out), optional :: f1
      !> Second derivative with respect to x0
      real(wp), intent(out), optional :: f2
      !> Third derivative with respect to x0, not available
      real(wp), intent(out), optional :: f3
      !> Fourth derivative with respect to x0, not available
      real(wp), intent(out), optional :: f4
      real(wp) :: lo, hi, width, u, one_minus_u
      real(wp) :: p_hi, a_hi, p_lo, a_lo
      real(wp) :: g, gp, gpp, s, s1m, exp_neg_g
      real(wp) :: dS_du, d2S_du, du_dx, d2u_dx2
      real(wp) :: bump, dbump, d2bump
      real(wp) :: direction

      if (present(f3) .or. present(f4)) error stop "moist smoothstep: sigmoid type has no third/fourth derivative"

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
         one_minus_u = (x0 - lo)/width
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
            ! Stable sigmoid evaluation for any sign of g
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

            ! u(x) = (hi - x)/width  =>  du/dx = -1/width,  d2u/dx2 = 0
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

   end subroutine smoothstep_sigmoid_eval

   !> Evaluate the polynomial smoothstep value and its first four scalar derivatives
   !>
   !> - `t` clamped to [0, 1]; value exactly 0/1 on the plateaus
   !> - Upper half by point symmetry, `S(t) = 1 - S(s)`, `S^(k)(t) = (-1)^(k+1) S^(k)(s)`,
   !>   `s = (hi - x0)/width`; plain Horner near t = 1 cancels to ~1e-10 and can exit [0, 1],
   !>   mirroring keeps the Horner argument at s <= 1/2
   !> - Derivatives masked to exactly 0 outside the open band, where the clamped
   !>   polynomial would leave orders above n nonzero
   !> - Band scaling 1/width per derivative order
   !>
   !> @param[in]  self  Polynomial smoothstep instance
   !> @param[in]  x0    Input value
   !> @param[out] f0    Smoothstep value
   !> @param[out] f1    First derivative with respect to x0
   !> @param[out] f2    Second derivative with respect to x0
   !> @param[out] f3    Third derivative with respect to x0
   !> @param[out] f4    Fourth derivative with respect to x0
   pure subroutine smoothstep_polynomial_eval(self, x0, f0, f1, f2, f3, f4)
      class(moist_math_smoothstep_polynomial_type), intent(in) :: self
      !> Input value
      real(wp), intent(in) :: x0
      !> Smoothstep value
      real(wp), intent(out) :: f0
      !> First derivative with respect to x0
      real(wp), intent(out), optional :: f1
      !> Second derivative with respect to x0
      real(wp), intent(out), optional :: f2
      !> Third derivative with respect to x0
      real(wp), intent(out), optional :: f3
      !> Fourth derivative with respect to x0
      real(wp), intent(out), optional :: f4
      real(wp) :: lo, hi, width, inv_w, inv_w2, t, s, inband, direction, bump, scale, mirror, p0
      logical :: upper
      integer :: degree

      degree = self%degree
      if (degree < 0) error stop "moist smoothstep: polynomial type used before new_smoothstep_polynomial"

      lo = min(self%from, self%to)
      hi = max(self%from, self%to)
      width = hi - lo
      direction = merge(-1.0_wp, 1.0_wp, self%from > self%to)

      if (width <= epsilon(1.0_wp)) then
         bump = merge(0.0_wp, 1.0_wp, x0 <= lo)
         f0 = merge(bump, 1.0_wp - bump, direction > 0.0_wp)
         if (present(f1)) f1 = 0.0_wp
         if (present(f2)) f2 = 0.0_wp
         if (present(f3)) f3 = 0.0_wp
         if (present(f4)) f4 = 0.0_wp
         return
      end if

      t = min(max((x0 - lo)/width, 0.0_wp), 1.0_wp)
      inband = merge(1.0_wp, 0.0_wp, x0 > lo .and. x0 < hi)
      inv_w = 1.0_wp/width
      inv_w2 = inv_w*inv_w

      ! Evaluate at the distance to the nearer edge; even orders flip sign when mirrored
      upper = t > 0.5_wp
      if (upper) then
         s = min(max((hi - x0)/width, 0.0_wp), 1.0_wp)
      else
         s = t
      end if
      mirror = merge(-1.0_wp, 1.0_wp, upper)

      p0 = horner(self%coef(:, 0), degree, s)
      bump = merge(1.0_wp - p0, p0, upper)
      f0 = merge(bump, 1.0_wp - bump, direction > 0.0_wp)

      scale = direction*inband
      if (present(f1)) f1 = scale*inv_w*horner(self%coef(:, 1), degree - 1, s)
      if (present(f2)) f2 = scale*mirror*inv_w2*horner(self%coef(:, 2), degree - 2, s)
      if (present(f3)) f3 = scale*(inv_w2*inv_w)*horner(self%coef(:, 3), degree - 3, s)
      if (present(f4)) f4 = scale*mirror*(inv_w2*inv_w2)*horner(self%coef(:, 4), degree - 4, s)

   end subroutine smoothstep_polynomial_eval

   !> Horner evaluation of `sum_{j=0}^{degree} coef(j) t^j`
   !>
   !> @param[in] coef    Coefficients, lower bound 0
   !> @param[in] degree  Polynomial degree, < 0 gives 0
   !> @param[in] t       Argument
   pure function horner(coef, degree, t) result(p)
      !> Coefficients, lower bound 0
      real(wp), intent(in) :: coef(0:)
      !> Polynomial degree
      integer, intent(in) :: degree
      !> Argument
      real(wp), intent(in) :: t
      !> Polynomial value
      real(wp) :: p
      integer :: j

      p = 0.0_wp
      do j = degree, 0, -1
         p = p*t + coef(j)
      end do

   end function horner

   !* ================================================================================= *!
   !*                      Base-type helpers: value and chain rule                      *!
   !* ================================================================================= *!

   !> Evaluate the smoothstep value via the abstract `eval` hook
   !>
   !> When from < to the value rises smoothly from 0 to 1
   !> When from > to the value falls smoothly from 1 to 0
   !>
   !> @param[in] self  Smoothstep instance
   !> @param[in] x0    Input value
   pure function smoothstep_f0(self, x0) result(val)
      class(moist_math_smoothstep_type), intent(in) :: self
      !> Input value
      real(wp), intent(in) :: x0
      real(wp) :: val

      call self%eval(x0, val)

   end function smoothstep_f0

   !> First derivative of the smoothstep with respect to centre
   !> coordinates, via the chain rule through x0
   !>
   !> dS/dR_A = S'(x0) * dx0/dR_A
   !>
   !> - Layout (3, ncenters), one Cartesian triple per centre
   !> - With `active`, only those centre columns are filled, the rest are 0
   !>
   !> @param[in] self    Smoothstep instance
   !> @param[in] x0      Input scalar value
   !> @param[in] x1      Gradient of x0 w.r.t. centre coordinates (3, ncenters)
   !> @param[in] active  Optional selected centre indices
   pure function smoothstep_f1_rA(self, x0, x1, active) result(grad)
      class(moist_math_smoothstep_type), intent(in) :: self
      !> Input scalar value
      real(wp), intent(in) :: x0
      !> Gradient of x0 w.r.t. centre coordinates (3, ncenters)
      real(wp), intent(in) :: x1(:, :)
      !> Optional selected centre indices
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

   end function smoothstep_f1_rA

   !> Second derivative (Hessian) of the smoothstep with respect
   !> to centre coordinates
   !>
   !> d^2S/(dR_A dR_B) = S''(x0) * dx0/dR_A * dx0/dR_B
   !>                   + S'(x0) * d^2x0/(dR_A dR_B)
   !>
   !> - Layout (3, 3, ncenters, ncenters), Cartesian components first
   !>
   !> @param[in] self  Smoothstep instance
   !> @param[in] x0    Input scalar value
   !> @param[in] x1    Gradient of x0 w.r.t. centre coordinates (3, ncenters)
   !> @param[in] x2    Hessian of x0 w.r.t. centre coordinates (3, 3, ncenters, ncenters)
   pure function smoothstep_f2_rArB(self, x0, x1, x2) result(hess)
      class(moist_math_smoothstep_type), intent(in) :: self
      !> Input scalar value
      real(wp), intent(in) :: x0
      !> Gradient of x0 w.r.t. centre coordinates (3, ncenters)
      real(wp), intent(in) :: x1(:, :)
      !> Hessian of x0 w.r.t. centre coordinates (3, 3, ncenters, ncenters)
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

   end function smoothstep_f2_rArB

end module moist_math_smoothstep
