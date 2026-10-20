!> Numerical kernels for HNC, KH and PSE2-4 MOZ closures
!>
!> - Adapted from generated kernels in moist_dev/pysofft (44874e7e)
!> - Closure objects, input checks and array dispatch live in moist_model_moz_closure
!> - Values and closure derivatives take matching full (u, t) or Ng (u_sr, t_s) inputs
!>   - return c or c_s in that representation
!>   - caller owns the split
!> - Local energy integrands and derivatives require physical h and c
!> - Direct correlations are evaluated by sweeps over one-dimensional blocks
!> - Loop keeps exponential arithmetic together for compiler vectorization
!> - Input/output arrays must have equal sizes and separate output storage
!> - Derivatives and local energy integrands are pointwise; integration belongs to the caller
!> - Exponential range reduction requires IEEE binary64, round-to-nearest and no fast-math
module moist_model_moz_closure_kernel
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: closure_hnc_c, closure_hnc_dc, closure_hnc_d2c, closure_hnc_phi
   public :: closure_kh_c, closure_kh_dc, closure_kh_d2c, closure_kh_phi
   public :: closure_pse2_c, closure_pse2_dc, closure_pse2_d2c, closure_pse2_phi
   public :: closure_pse3_c, closure_pse3_dc, closure_pse3_d2c, closure_pse3_phi
   public :: closure_pse4_c, closure_pse4_dc, closure_pse4_d2c, closure_pse4_phi


   !> Reciprocal of ln2, log2(e)
   real(wp), parameter :: exp_log2e = 1.4426950408889634_wp
   !> High part of ln2 for range reduction
   real(wp), parameter :: exp_ln2_hi = 6.93147180369123816490e-01_wp
   !> Low part of ln2 for range reduction
   real(wp), parameter :: exp_ln2_lo = 1.90821492927058770002e-10_wp
   !> Rounding shifter 1.5*2**52 for integer rounding of the exponent
   real(wp), parameter :: exp_shifter = 6755399441055744.0_wp
   !> IEEE bits of exp_shifter
   integer(int64), parameter :: exp_shifter_bits = int(z'4338000000000000', int64)
   !> 2**52, moves the exponent field into place
   integer(int64), parameter :: exp_two52 = 4503599627370496_int64

   public :: closure_du, closure_phi_dc, closure_hnc_phi_dh, closure_kh_phi_dh
   public :: closure_pse2_phi_dh, closure_pse3_phi_dh, closure_pse4_phi_dh
   public :: closure_dc_sweep, closure_d2c_sweep, closure_c_dc_sweep

contains

   !> HNC derivative with respect to t
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   pure elemental function closure_hnc_dc(u, t, beta, lambda) result(dc)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Derivative dc/dt
      real(wp) :: dc
      !> Closure exponent t_star = -beta*lambda*u + t
      real(wp) :: t_star

      t_star = -beta*lambda*u + t

      dc = fast_exp(t_star) - 1.0_wp
   end function closure_hnc_dc

   !> HNC second derivative with respect to t
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   pure elemental function closure_hnc_d2c(u, t, beta, lambda) result(d2c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Second derivative d2c/dt2
      real(wp) :: d2c
      !> Closure exponent t_star = -beta*lambda*u + t
      real(wp) :: t_star

      t_star = -beta*lambda*u + t

      d2c = fast_exp(t_star)
   end function closure_hnc_d2c

   !> HNC local energy integrand
   !>
   !> - Density factors and spatial integration are supplied by the caller
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function closure_hnc_phi(h, c) result(phi)
      !> Physical total correlation, dimensionless
      real(wp), intent(in) :: h
      !> Physical direct correlation, dimensionless
      real(wp), intent(in) :: c
      !> Dimensionless local energy integrand
      real(wp) :: phi

      phi = -1.0_wp/2.0_wp*c*h + (1.0_wp/2.0_wp)*h**2 - c
   end function closure_hnc_phi

   !> KH derivative with respect to t
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   pure elemental function closure_kh_dc(u, t, beta, lambda) result(dc)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Derivative dc/dt
      real(wp) :: dc
      !> Closure exponent t_star = -beta*lambda*u + t
      real(wp) :: t_star

      t_star = -beta*lambda*u + t

      if (t_star <= 0.0_wp) then
         dc = fast_exp(t_star) - 1.0_wp
      else
         dc = 0.0_wp
      end if
   end function closure_kh_dc

   !> KH second derivative with respect to t
   !>
   !> - Jumps at t_star = 0
   !> - c and dc/dt are continuous there
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   pure elemental function closure_kh_d2c(u, t, beta, lambda) result(d2c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Second derivative d2c/dt2
      real(wp) :: d2c
      !> Closure exponent t_star = -beta*lambda*u + t
      real(wp) :: t_star

      t_star = -beta*lambda*u + t

      if (t_star <= 0.0_wp) then
         d2c = fast_exp(t_star)
      else
         d2c = 0.0_wp
      end if
   end function closure_kh_d2c

   !> KH local energy integrand
   !>
   !> - Density factors and spatial integration are supplied by the caller
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function closure_kh_phi(h, c) result(phi)
      !> Physical total correlation, dimensionless
      real(wp), intent(in) :: h
      !> Physical direct correlation, dimensionless
      real(wp), intent(in) :: c
      !> Dimensionless local energy integrand
      real(wp) :: phi

      if (h <= 0.0_wp) then
         phi = -1.0_wp/2.0_wp*c*h + (1.0_wp/2.0_wp)*h**2 - c
      else
         phi = -1.0_wp/2.0_wp*c*h - c
      end if
   end function closure_kh_phi

   !> PSE2 derivative with respect to t
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   pure elemental function closure_pse2_dc(u, t, beta, lambda) result(dc)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Derivative dc/dt
      real(wp) :: dc
      !> Closure exponent t_star = -beta*lambda*u + t
      real(wp) :: t_star

      t_star = -beta*lambda*u + t

      if (t_star <= 0.0_wp) then
         dc = fast_exp(t_star) - 1.0_wp
      else
         dc = t_star
      end if
   end function closure_pse2_dc

   !> PSE2 second derivative with respect to t
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   pure elemental function closure_pse2_d2c(u, t, beta, lambda) result(d2c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Second derivative d2c/dt2
      real(wp) :: d2c
      !> Closure exponent t_star = -beta*lambda*u + t
      real(wp) :: t_star

      t_star = -beta*lambda*u + t

      if (t_star <= 0.0_wp) then
         d2c = fast_exp(t_star)
      else
         d2c = 1.0_wp
      end if
   end function closure_pse2_d2c

   !> PSE2 local energy integrand
   !>
   !> - Density factors and spatial integration are supplied by the caller
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function closure_pse2_phi(h, c) result(phi)
      !> Physical total correlation, dimensionless
      real(wp), intent(in) :: h
      !> Physical direct correlation, dimensionless
      real(wp), intent(in) :: c
      !> Dimensionless local energy integrand
      real(wp) :: phi

      if (h <= 0.0_wp) then
         phi = -1.0_wp/2.0_wp*c*h + (1.0_wp/2.0_wp)*h**2 - c
      else
         phi = -1.0_wp/2.0_wp*c*h + (1.0_wp/2.0_wp)*h**2 - 1.0_wp/6.0_wp*pse_invert(h, &
            2)**3 - c
      end if
   end function closure_pse2_phi

   !> PSE3 derivative with respect to t
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   pure elemental function closure_pse3_dc(u, t, beta, lambda) result(dc)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Derivative dc/dt
      real(wp) :: dc
      !> Closure exponent t_star = -beta*lambda*u + t
      real(wp) :: t_star

      t_star = -beta*lambda*u + t

      if (t_star <= 0.0_wp) then
         dc = fast_exp(t_star) - 1.0_wp
      else
         dc = t_star*((1.0_wp/6.0_wp)*t_star + 1.0_wp/2.0_wp) + t_star*((1.0_wp/3.0_wp) &
            *t_star + 1.0_wp/2.0_wp)
      end if
   end function closure_pse3_dc

   !> PSE3 second derivative with respect to t
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   pure elemental function closure_pse3_d2c(u, t, beta, lambda) result(d2c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Second derivative d2c/dt2
      real(wp) :: d2c
      !> Closure exponent t_star = -beta*lambda*u + t
      real(wp) :: t_star

      t_star = -beta*lambda*u + t

      if (t_star <= 0.0_wp) then
         d2c = fast_exp(t_star)
      else
         d2c = 1.0_wp + t_star
      end if
   end function closure_pse3_d2c

   !> PSE3 local energy integrand
   !>
   !> - Density factors and spatial integration are supplied by the caller
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function closure_pse3_phi(h, c) result(phi)
      !> Physical total correlation, dimensionless
      real(wp), intent(in) :: h
      !> Physical direct correlation, dimensionless
      real(wp), intent(in) :: c
      !> Dimensionless local energy integrand
      real(wp) :: phi

      if (h <= 0.0_wp) then
         phi = -1.0_wp/2.0_wp*c*h + (1.0_wp/2.0_wp)*h**2 - c
      else
         phi = -1.0_wp/2.0_wp*c*h + (1.0_wp/2.0_wp)*h**2 - 1.0_wp/24.0_wp*pse_invert(h, &
            3)**4 - c
      end if
   end function closure_pse3_phi

   !> PSE4 derivative with respect to t
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   pure elemental function closure_pse4_dc(u, t, beta, lambda) result(dc)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Derivative dc/dt
      real(wp) :: dc
      !> Closure exponent t_star = -beta*lambda*u + t
      real(wp) :: t_star
      real(wp) :: x0

      t_star = -beta*lambda*u + t

      if (t_star <= 0.0_wp) then
         dc = fast_exp(t_star) - 1.0_wp
      else
         x0 = t_star*((1.0_wp/24.0_wp)*t_star + 1.0_wp/6.0_wp) + 1.0_wp/2.0_wp
         dc = t_star*x0 + t_star*(t_star*((1.0_wp/12.0_wp)*t_star + 1.0_wp/6.0_wp) + &
            x0)
      end if
   end function closure_pse4_dc

   !> PSE4 second derivative with respect to t
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   pure elemental function closure_pse4_d2c(u, t, beta, lambda) result(d2c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Second derivative d2c/dt2
      real(wp) :: d2c
      !> Closure exponent t_star = -beta*lambda*u + t
      real(wp) :: t_star
      real(wp) :: x0

      t_star = -beta*lambda*u + t

      if (t_star <= 0.0_wp) then
         d2c = fast_exp(t_star)
      else
         x0 = 2.0_wp*t_star
         d2c = t_star*((1.0_wp/4.0_wp)*t_star + 1.0_wp/3.0_wp) + x0*((1.0_wp/24.0_wp)* &
            t_star + 1.0_wp/6.0_wp) + x0*((1.0_wp/12.0_wp)*t_star + 1.0_wp/6.0_wp) &
            + 1.0_wp
      end if
   end function closure_pse4_d2c

   !> PSE4 local energy integrand
   !>
   !> - Density factors and spatial integration are supplied by the caller
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function closure_pse4_phi(h, c) result(phi)
      !> Physical total correlation, dimensionless
      real(wp), intent(in) :: h
      !> Physical direct correlation, dimensionless
      real(wp), intent(in) :: c
      !> Dimensionless local energy integrand
      real(wp) :: phi

      if (h <= 0.0_wp) then
         phi = -1.0_wp/2.0_wp*c*h + (1.0_wp/2.0_wp)*h**2 - c
      else
         phi = -1.0_wp/2.0_wp*c*h + (1.0_wp/2.0_wp)*h**2 - 1.0_wp/120.0_wp*pse_invert(h &
            , 4)**5 - c
      end if
   end function closure_pse4_phi

   !> HNC direct correlation over a block of points
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   !> @param[out] c Full c or Ng short-range c_s, dimensionless
   pure subroutine closure_hnc_c(u, t, beta, lambda, c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u(:)
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t(:)
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Full c or Ng short-range c_s, dimensionless
      real(wp), intent(out) :: c(:)
      !> Point index
      integer :: i
      !> Closure exponent t_star = -beta*lambda*u + t, and exp of the branch argument
      real(wp) :: t_star, e
      !> Clamped argument, multiple of ln2, reduced argument, polynomial
      real(wp) :: ex_c, ex_k, ex_r, ex_p
      !> IEEE bits of the shifted multiple of ln2
      integer(int64) :: ex_n

      do i = 1, size(t)
         t_star = -beta*lambda*u(i) + t(i)
         ex_c = min(max(t_star, -708.0_wp), 710.0_wp)
         ex_k = ex_c*exp_log2e + exp_shifter
         ex_n = transfer(ex_k, ex_n)
         ex_k = ex_k - exp_shifter
         ex_r = (ex_c - ex_k*exp_ln2_hi) - ex_k*exp_ln2_lo
         ex_p = 1.60590438368216133e-10_wp
         ex_p = 2.08767569878681002e-9_wp + ex_r*ex_p
         ex_p = 2.50521083854417202e-8_wp + ex_r*ex_p
         ex_p = 2.75573192239858883e-7_wp + ex_r*ex_p
         ex_p = 2.75573192239858925e-6_wp + ex_r*ex_p
         ex_p = 2.48015873015873016e-5_wp + ex_r*ex_p
         ex_p = 1.98412698412698413e-4_wp + ex_r*ex_p
         ex_p = 1.38888888888888894e-3_wp + ex_r*ex_p
         ex_p = 8.33333333333333322e-3_wp + ex_r*ex_p
         ex_p = 4.16666666666666644e-2_wp + ex_r*ex_p
         ex_p = 1.66666666666666657e-1_wp + ex_r*ex_p
         ex_p = 5.00000000000000000e-1_wp + ex_r*ex_p
         ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
         ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
         e = (2.0_wp*ex_p)*transfer((ex_n - exp_shifter_bits + 1022_int64)*exp_two52, 1.0_wp)
         c(i) = (e - 1.0_wp) - t(i)
      end do
   end subroutine closure_hnc_c

   !> KH direct correlation over a block of points
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   !> @param[out] c Full c or Ng short-range c_s, dimensionless
   pure subroutine closure_kh_c(u, t, beta, lambda, c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u(:)
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t(:)
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Full c or Ng short-range c_s, dimensionless
      real(wp), intent(out) :: c(:)
      !> Point index
      integer :: i
      !> Closure exponent t_star = -beta*lambda*u + t, and exp of the branch argument
      real(wp) :: t_star, e
      !> Clamped argument, multiple of ln2, reduced argument, polynomial
      real(wp) :: ex_c, ex_k, ex_r, ex_p
      !> IEEE bits of the shifted multiple of ln2
      integer(int64) :: ex_n

      do i = 1, size(t)
         t_star = -beta*lambda*u(i) + t(i)
         ex_c = min(max(min(t_star, 0.0_wp), -708.0_wp), 710.0_wp)
         ex_k = ex_c*exp_log2e + exp_shifter
         ex_n = transfer(ex_k, ex_n)
         ex_k = ex_k - exp_shifter
         ex_r = (ex_c - ex_k*exp_ln2_hi) - ex_k*exp_ln2_lo
         ex_p = 1.60590438368216133e-10_wp
         ex_p = 2.08767569878681002e-9_wp + ex_r*ex_p
         ex_p = 2.50521083854417202e-8_wp + ex_r*ex_p
         ex_p = 2.75573192239858883e-7_wp + ex_r*ex_p
         ex_p = 2.75573192239858925e-6_wp + ex_r*ex_p
         ex_p = 2.48015873015873016e-5_wp + ex_r*ex_p
         ex_p = 1.98412698412698413e-4_wp + ex_r*ex_p
         ex_p = 1.38888888888888894e-3_wp + ex_r*ex_p
         ex_p = 8.33333333333333322e-3_wp + ex_r*ex_p
         ex_p = 4.16666666666666644e-2_wp + ex_r*ex_p
         ex_p = 1.66666666666666657e-1_wp + ex_r*ex_p
         ex_p = 5.00000000000000000e-1_wp + ex_r*ex_p
         ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
         ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
         e = (2.0_wp*ex_p)*transfer((ex_n - exp_shifter_bits + 1022_int64)*exp_two52, 1.0_wp)
         c(i) = ((e - 1.0_wp) - t(i)) + max(t_star, 0.0_wp)
      end do
   end subroutine closure_kh_c

   !> PSE2 direct correlation over a block of points
   !>
   !> - Combine exp(min(x, 0)) with the polynomial at max(x, 0) in one loop
   !> - Evaluates both branches without a per-point conditional
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   !> @param[out] c Full c or Ng short-range c_s, dimensionless
   pure subroutine closure_pse2_c(u, t, beta, lambda, c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u(:)
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t(:)
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Full c or Ng short-range c_s, dimensionless
      real(wp), intent(out) :: c(:)
      !> Point index
      integer :: i
      !> Closure exponent t_star = -beta*lambda*u + t, and exp of the branch argument
      real(wp) :: t_star, e, positive
      !> Clamped argument, multiple of ln2, reduced argument, polynomial
      real(wp) :: ex_c, ex_k, ex_r, ex_p
      !> IEEE bits of the shifted multiple of ln2
      integer(int64) :: ex_n

      do i = 1, size(t)
         t_star = -beta*lambda*u(i) + t(i)
         ex_c = min(max(min(t_star, 0.0_wp), -708.0_wp), 710.0_wp)
         ex_k = ex_c*exp_log2e + exp_shifter
         ex_n = transfer(ex_k, ex_n)
         ex_k = ex_k - exp_shifter
         ex_r = (ex_c - ex_k*exp_ln2_hi) - ex_k*exp_ln2_lo
         ex_p = 1.60590438368216133e-10_wp
         ex_p = 2.08767569878681002e-9_wp + ex_r*ex_p
         ex_p = 2.50521083854417202e-8_wp + ex_r*ex_p
         ex_p = 2.75573192239858883e-7_wp + ex_r*ex_p
         ex_p = 2.75573192239858925e-6_wp + ex_r*ex_p
         ex_p = 2.48015873015873016e-5_wp + ex_r*ex_p
         ex_p = 1.98412698412698413e-4_wp + ex_r*ex_p
         ex_p = 1.38888888888888894e-3_wp + ex_r*ex_p
         ex_p = 8.33333333333333322e-3_wp + ex_r*ex_p
         ex_p = 4.16666666666666644e-2_wp + ex_r*ex_p
         ex_p = 1.66666666666666657e-1_wp + ex_r*ex_p
         ex_p = 5.00000000000000000e-1_wp + ex_r*ex_p
         ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
         ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
         e = (2.0_wp*ex_p)*transfer((ex_n - exp_shifter_bits + 1022_int64)*exp_two52, 1.0_wp)
         positive = max(t_star, 0.0_wp)
         c(i) = ((e - 1.0_wp) - t(i)) + positive*((1.0_wp/2.0_wp)*positive + 1.0_wp)
      end do
   end subroutine closure_pse2_c

   !> PSE3 direct correlation over a block of points
   !>
   !> - Combine exp(min(x, 0)) with the polynomial at max(x, 0) in one loop
   !> - Evaluates both branches without a per-point conditional
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   !> @param[out] c Full c or Ng short-range c_s, dimensionless
   pure subroutine closure_pse3_c(u, t, beta, lambda, c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u(:)
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t(:)
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Full c or Ng short-range c_s, dimensionless
      real(wp), intent(out) :: c(:)
      !> Point index
      integer :: i
      !> Closure exponent t_star = -beta*lambda*u + t, and exp of the branch argument
      real(wp) :: t_star, e, positive
      !> Clamped argument, multiple of ln2, reduced argument, polynomial
      real(wp) :: ex_c, ex_k, ex_r, ex_p
      !> IEEE bits of the shifted multiple of ln2
      integer(int64) :: ex_n

      do i = 1, size(t)
         t_star = -beta*lambda*u(i) + t(i)
         ex_c = min(max(min(t_star, 0.0_wp), -708.0_wp), 710.0_wp)
         ex_k = ex_c*exp_log2e + exp_shifter
         ex_n = transfer(ex_k, ex_n)
         ex_k = ex_k - exp_shifter
         ex_r = (ex_c - ex_k*exp_ln2_hi) - ex_k*exp_ln2_lo
         ex_p = 1.60590438368216133e-10_wp
         ex_p = 2.08767569878681002e-9_wp + ex_r*ex_p
         ex_p = 2.50521083854417202e-8_wp + ex_r*ex_p
         ex_p = 2.75573192239858883e-7_wp + ex_r*ex_p
         ex_p = 2.75573192239858925e-6_wp + ex_r*ex_p
         ex_p = 2.48015873015873016e-5_wp + ex_r*ex_p
         ex_p = 1.98412698412698413e-4_wp + ex_r*ex_p
         ex_p = 1.38888888888888894e-3_wp + ex_r*ex_p
         ex_p = 8.33333333333333322e-3_wp + ex_r*ex_p
         ex_p = 4.16666666666666644e-2_wp + ex_r*ex_p
         ex_p = 1.66666666666666657e-1_wp + ex_r*ex_p
         ex_p = 5.00000000000000000e-1_wp + ex_r*ex_p
         ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
         ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
         e = (2.0_wp*ex_p)*transfer((ex_n - exp_shifter_bits + 1022_int64)*exp_two52, 1.0_wp)
         positive = max(t_star, 0.0_wp)
         c(i) = ((e - 1.0_wp) - t(i)) + positive*(positive*((1.0_wp/6.0_wp)*positive + 1.0_wp/2.0_wp) + 1.0_wp)
      end do
   end subroutine closure_pse3_c

   !> PSE4 direct correlation over a block of points
   !>
   !> - Combine exp(min(x, 0)) with the polynomial at max(x, 0) in one loop
   !> - Evaluates both branches without a per-point conditional
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   !> @param[out] c Full c or Ng short-range c_s, dimensionless
   pure subroutine closure_pse4_c(u, t, beta, lambda, c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u(:)
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t(:)
      !> Inverse thermal energy, reciprocal units of the potential
      real(wp), intent(in) :: beta
      !> Potential scaling parameter in [0, 1]
      real(wp), intent(in) :: lambda
      !> Full c or Ng short-range c_s, dimensionless
      real(wp), intent(out) :: c(:)
      !> Point index
      integer :: i
      !> Closure exponent t_star = -beta*lambda*u + t, and exp of the branch argument
      real(wp) :: t_star, e, positive
      !> Clamped argument, multiple of ln2, reduced argument, polynomial
      real(wp) :: ex_c, ex_k, ex_r, ex_p
      !> IEEE bits of the shifted multiple of ln2
      integer(int64) :: ex_n

      do i = 1, size(t)
         t_star = -beta*lambda*u(i) + t(i)
         ex_c = min(max(min(t_star, 0.0_wp), -708.0_wp), 710.0_wp)
         ex_k = ex_c*exp_log2e + exp_shifter
         ex_n = transfer(ex_k, ex_n)
         ex_k = ex_k - exp_shifter
         ex_r = (ex_c - ex_k*exp_ln2_hi) - ex_k*exp_ln2_lo
         ex_p = 1.60590438368216133e-10_wp
         ex_p = 2.08767569878681002e-9_wp + ex_r*ex_p
         ex_p = 2.50521083854417202e-8_wp + ex_r*ex_p
         ex_p = 2.75573192239858883e-7_wp + ex_r*ex_p
         ex_p = 2.75573192239858925e-6_wp + ex_r*ex_p
         ex_p = 2.48015873015873016e-5_wp + ex_r*ex_p
         ex_p = 1.98412698412698413e-4_wp + ex_r*ex_p
         ex_p = 1.38888888888888894e-3_wp + ex_r*ex_p
         ex_p = 8.33333333333333322e-3_wp + ex_r*ex_p
         ex_p = 4.16666666666666644e-2_wp + ex_r*ex_p
         ex_p = 1.66666666666666657e-1_wp + ex_r*ex_p
         ex_p = 5.00000000000000000e-1_wp + ex_r*ex_p
         ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
         ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
         e = (2.0_wp*ex_p)*transfer((ex_n - exp_shifter_bits + 1022_int64)*exp_two52, 1.0_wp)
         positive = max(t_star, 0.0_wp)
         c(i) = ((e - 1.0_wp) - t(i)) + &
            & positive*(positive*(positive*((1.0_wp/24.0_wp)*positive + 1.0_wp/6.0_wp) + 1.0_wp/2.0_wp) + 1.0_wp)
      end do
   end subroutine closure_pse4_c

   !> Closure derivative over a block of points
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   !> @param[in] kind Closure kind: 0 HNC, 1 KH, 2-4 PSE order; other values give NaN
   !> @param[out] dc Derivative dc/dt
   pure subroutine closure_dc_sweep(u, t, beta, lambda, kind, dc)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u(:)
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t(:)
      !> Inverse thermal energy and potential scaling
      real(wp), intent(in) :: beta, lambda
      !> Closure kind: 0 HNC, 1 KH, 2-4 PSE order
      integer, intent(in) :: kind
      !> Derivative dc/dt
      real(wp), intent(out) :: dc(:)

      select case (kind)
      case (0)
         dc = closure_hnc_dc(u, t, beta, lambda)
      case (1)
         dc = closure_kh_dc(u, t, beta, lambda)
      case (2)
         dc = closure_pse2_dc(u, t, beta, lambda)
      case (3)
         dc = closure_pse3_dc(u, t, beta, lambda)
      case (4)
         dc = closure_pse4_dc(u, t, beta, lambda)
      case default
         dc = ieee_value(0.0_wp, ieee_quiet_nan)
      end select
   end subroutine closure_dc_sweep

   !> Closure second derivative over a block of points
   !>
   !> - KH jumps at t_star = 0; PSE-n is continuous there
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   !> @param[in] kind Closure kind: 0 HNC, 1 KH, 2-4 PSE order; other values give NaN
   !> @param[out] d2c Second derivative d2c/dt2
   pure subroutine closure_d2c_sweep(u, t, beta, lambda, kind, d2c)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u(:)
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t(:)
      !> Inverse thermal energy and potential scaling
      real(wp), intent(in) :: beta, lambda
      !> Closure kind: 0 HNC, 1 KH, 2-4 PSE order
      integer, intent(in) :: kind
      !> Second derivative d2c/dt2
      real(wp), intent(out) :: d2c(:)

      select case (kind)
      case (0)
         d2c = closure_hnc_d2c(u, t, beta, lambda)
      case (1)
         d2c = closure_kh_d2c(u, t, beta, lambda)
      case (2)
         d2c = closure_pse2_d2c(u, t, beta, lambda)
      case (3)
         d2c = closure_pse3_d2c(u, t, beta, lambda)
      case (4)
         d2c = closure_pse4_d2c(u, t, beta, lambda)
      case default
         d2c = ieee_value(0.0_wp, ieee_quiet_nan)
      end select
   end subroutine closure_d2c_sweep

   !> Direct correlation and its derivative over a block, sharing one exponential
   !>
   !> - Values equal the closure_*_c sweeps; derivatives equal closure_*_dc up to rounding
   !>
   !> @param[in] u Full u or Ng short-range u_sr, energy units
   !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of the potential
   !> @param[in] lambda Potential scaling parameter in [0, 1]
   !> @param[in] kind Closure kind: 0 HNC, 1 KH, 2-4 PSE order; other values give NaN
   !> @param[out] c Full c or Ng short-range c_s, dimensionless
   !> @param[out] dc Derivative dc/dt
   pure subroutine closure_c_dc_sweep(u, t, beta, lambda, kind, c, dc)
      !> Full u or Ng short-range u_sr, energy units
      real(wp), intent(in) :: u(:)
      !> Matching full t = h - c or Ng t_s, dimensionless
      real(wp), intent(in) :: t(:)
      !> Inverse thermal energy and potential scaling
      real(wp), intent(in) :: beta, lambda
      !> Closure kind: 0 HNC, 1 KH, 2-4 PSE order
      integer, intent(in) :: kind
      !> Full c or Ng short-range c_s, dimensionless
      real(wp), intent(out) :: c(:)
      !> Derivative dc/dt
      real(wp), intent(out) :: dc(:)
      !> Point index
      integer :: i
      !> Closure exponent, exp of the branch argument, positive part
      real(wp) :: t_star, e, positive

      select case (kind)
      case (0)
         do i = 1, size(t)
            t_star = -beta*lambda*u(i) + t(i)
            e = fast_exp(t_star)
            c(i) = (e - 1.0_wp) - t(i)
            dc(i) = e - 1.0_wp
         end do
      case (1)
         ! exp(min(t_star, 0)) is exactly one on the linear branch
         do i = 1, size(t)
            t_star = -beta*lambda*u(i) + t(i)
            e = fast_exp(min(t_star, 0.0_wp))
            c(i) = ((e - 1.0_wp) - t(i)) + max(t_star, 0.0_wp)
            dc(i) = e - 1.0_wp
         end do
      case (2)
         do i = 1, size(t)
            t_star = -beta*lambda*u(i) + t(i)
            e = fast_exp(min(t_star, 0.0_wp))
            positive = max(t_star, 0.0_wp)
            c(i) = ((e - 1.0_wp) - t(i)) + positive*((1.0_wp/2.0_wp)*positive + 1.0_wp)
            dc(i) = (e - 1.0_wp) + positive
         end do
      case (3)
         do i = 1, size(t)
            t_star = -beta*lambda*u(i) + t(i)
            e = fast_exp(min(t_star, 0.0_wp))
            positive = max(t_star, 0.0_wp)
            c(i) = ((e - 1.0_wp) - t(i)) + positive*(positive*((1.0_wp/6.0_wp)*positive + 1.0_wp/2.0_wp) + 1.0_wp)
            dc(i) = (e - 1.0_wp) + positive*((1.0_wp/2.0_wp)*positive + 1.0_wp)
         end do
      case (4)
         do i = 1, size(t)
            t_star = -beta*lambda*u(i) + t(i)
            e = fast_exp(min(t_star, 0.0_wp))
            positive = max(t_star, 0.0_wp)
            c(i) = ((e - 1.0_wp) - t(i)) + &
               & positive*(positive*(positive*((1.0_wp/24.0_wp)*positive + 1.0_wp/6.0_wp) + 1.0_wp/2.0_wp) + 1.0_wp)
            dc(i) = (e - 1.0_wp) + positive*(positive*((1.0_wp/6.0_wp)*positive + 1.0_wp/2.0_wp) + 1.0_wp)
         end do
      case default
         c = ieee_value(0.0_wp, ieee_quiet_nan)
         dc = c
      end select
   end subroutine closure_c_dc_sweep

   !> Exponential approximation shared by the closure kernels
   !>
   !> - Arguments below -708 are clamped
   !> - Result scaled as (2p)*2**(n-1), exponent field stays below the Inf code 2047
   !> - Finite up to ln(huge) ~ 709.78
   !> - Above that, up to the 710 clamp, the final product overflows to +Inf
   !>
   !> @param[in] x Dimensionless exponent
   pure elemental function fast_exp(x) result(y)
      !> Dimensionless exponent
      real(wp), intent(in) :: x
      !> Approximate exponential
      real(wp) :: y
      !> Clamped argument, multiple of ln2, reduced argument, polynomial
      real(wp) :: ex_c, ex_k, ex_r, ex_p
      !> IEEE bits of the shifted multiple of ln2
      integer(int64) :: ex_n

      ex_c = min(max(x, -708.0_wp), 710.0_wp)
      ex_k = ex_c*exp_log2e + exp_shifter
      ex_n = transfer(ex_k, ex_n)
      ex_k = ex_k - exp_shifter
      ex_r = (ex_c - ex_k*exp_ln2_hi) - ex_k*exp_ln2_lo
      ex_p = 1.60590438368216133e-10_wp
      ex_p = 2.08767569878681002e-9_wp + ex_r*ex_p
      ex_p = 2.50521083854417202e-8_wp + ex_r*ex_p
      ex_p = 2.75573192239858883e-7_wp + ex_r*ex_p
      ex_p = 2.75573192239858925e-6_wp + ex_r*ex_p
      ex_p = 2.48015873015873016e-5_wp + ex_r*ex_p
      ex_p = 1.98412698412698413e-4_wp + ex_r*ex_p
      ex_p = 1.38888888888888894e-3_wp + ex_r*ex_p
      ex_p = 8.33333333333333322e-3_wp + ex_r*ex_p
      ex_p = 4.16666666666666644e-2_wp + ex_r*ex_p
      ex_p = 1.66666666666666657e-1_wp + ex_r*ex_p
      ex_p = 5.00000000000000000e-1_wp + ex_r*ex_p
      ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
      ex_p = 1.00000000000000000e0_wp + ex_r*ex_p
      y = (2.0_wp*ex_p)*transfer((ex_n - exp_shifter_bits + 1022_int64)*exp_two52, 1.0_wp)
   end function fast_exp

   !> Truncated exponential series evaluated in Horner form
   !>
   !> @param[in] x Dimensionless series argument
   !> @param[in] n Truncation order
   pure function pse_series(x, n) result(s)
      !> Dimensionless series argument
      real(wp), intent(in) :: x
      !> Truncation order
      integer, intent(in) :: n
      !> Truncated series value
      real(wp) :: s
      integer :: i

      s = 1.0_wp
      do i = n, 1, -1
         s = 1.0_wp + x/real(i, wp)*s
      end do
   end function pse_series

   !> Recover the positive PSE exponent from the total correlation
   !>
   !> Newton iteration solves s_n(t_star) = 1 + h for h > 0
   !>
   !> @param[in] h Physical total correlation, dimensionless; h > 0
   !> @param[in] n PSE expansion order
   pure function pse_invert(h, n) result(tstar)
      !> Physical total correlation, dimensionless; h > 0
      real(wp), intent(in) :: h
      !> PSE expansion order
      integer, intent(in) :: n
      !> Recovered closure exponent
      real(wp) :: tstar
      real(wp) :: f, df, dt
      integer :: iter
      real(wp), parameter :: tol = 1.0e-13_wp

      tstar = h
      do iter = 1, 100
         df = pse_series(tstar, n - 1)
         if (df <= 0.0_wp) exit
         f = pse_series(tstar, n) - (1.0_wp + h)
         dt = f/df
         tstar = tstar - dt
         if (tstar < 0.0_wp) tstar = 0.0_wp
         if (abs(dt) <= tol*(abs(tstar) + 1.0_wp)) exit
      end do
   end function pse_invert

   !> Potential derivative for closures depending on x = t - beta*lambda*u
   !>
   !> @param[in] dc_dt Derivative dc/dt
   !> @param[in] beta Inverse thermal energy, reciprocal units of u
   !> @param[in] lambda Potential scaling parameter, between zero and one
   pure elemental function closure_du(dc_dt, beta, lambda) result(value)
      !> Derivative with respect to t, inverse temperature and potential scale
      real(wp), intent(in) :: dc_dt, beta, lambda
      !> Derivative with respect to the supplied potential at fixed t and split
      real(wp) :: value
      value = -beta*lambda*(dc_dt + 1.0_wp)
   end function closure_du

   !> Partial derivative of the HNC, KH and PSE integrands with respect to c
   !>
   !> @param[in] h Physical total correlation, dimensionless
   pure elemental function closure_phi_dc(h) result(value)
      !> Physical total correlation
      real(wp), intent(in) :: h
      !> Dimensionless local derivative
      real(wp) :: value
      value = -1.0_wp - 0.5_wp*h
   end function closure_phi_dc

   !> HNC integrand partial derivative with respect to h at fixed c
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function closure_hnc_phi_dh(h, c) result(value)
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Dimensionless local derivative
      real(wp) :: value
      value = h - 0.5_wp*c
   end function closure_hnc_phi_dh

   !> KH integrand partial derivative with respect to h at fixed c
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function closure_kh_phi_dh(h, c) result(value)
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Dimensionless local derivative
      real(wp) :: value
      value = min(h, 0.0_wp) - 0.5_wp*c
   end function closure_kh_phi_dh

   !> PSE integrand derivative, differentiating the positive-branch inverse
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   !> @param[in] order PSE expansion order, 2 to 4
   pure elemental function pse_phi_dh(h, c, order) result(value)
      !> Physical correlations and polynomial order
      real(wp), intent(in) :: h, c
      !> PSE expansion order, 2 to 4
      integer, intent(in) :: order
      !> Dimensionless local derivative
      real(wp) :: value
      real(wp) :: x, term, poly
      integer :: j
      value = h - 0.5_wp*c
      if (h <= 0.0_wp) return
      x = pse_invert(h, order)
      term = 1.0_wp
      poly = term
      do j = 1, order - 1
         term = term*x/real(j, wp)
         poly = poly + term
      end do
      term = term*x/real(order, wp)
      value = value - term/poly
   end function pse_phi_dh

   !> PSE2 integrand partial derivative with respect to h at fixed c
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function closure_pse2_phi_dh(h, c) result(value)
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Dimensionless local derivative
      real(wp) :: value
      value = pse_phi_dh(h, c, 2)
   end function closure_pse2_phi_dh

   !> PSE3 integrand partial derivative with respect to h at fixed c
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function closure_pse3_phi_dh(h, c) result(value)
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Dimensionless local derivative
      real(wp) :: value
      value = pse_phi_dh(h, c, 3)
   end function closure_pse3_phi_dh

   !> PSE4 integrand partial derivative with respect to h at fixed c
   !>
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function closure_pse4_phi_dh(h, c) result(value)
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Dimensionless local derivative
      real(wp) :: value
      value = pse_phi_dh(h, c, 4)
   end function closure_pse4_phi_dh

end module moist_model_moz_closure_kernel
