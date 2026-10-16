!> Pointwise UV response kernels, independent of grid dimension and layout
!>
!> h = c_s + t_s, c = c_s - beta u_lr and c_s = C(u_sr,t_s).
!> The closure supplies dh = dPhi/dh (fixed c) and dc = dPhi/dc (fixed h),
!> where Phi is the dimensionless free-energy integrand. Callers apply energy
!> seeds and may add the moment correction seed beta*dE/dD to dc.
!> Neither derivative includes density or quadrature weights:
!> a = rho_site * quadrature_weight / beta is the point's energy weight.
module moist_model_moz_response_kernel
   use mctc_env, only: wp
   implicit none(type, external)
   private
   public :: uv_response_rhs, uv_response_sr, uv_response_lr, uv_response_measure

contains

   !> Explicit quadrature sensitivity of the selected energy and direct moment.
   !> phi is the dimensionless closure integrand; shift = beta*dE/dD.
   !> rho/beta has energy/volume units.
   pure elemental function uv_response_measure(rho, beta, energy_seed, phi, shift, c) result(w)
      real(wp), intent(in) :: rho, beta, energy_seed, phi, shift, c
      real(wp) :: w
      w = (rho/beta)*(energy_seed*phi + shift*c)
   end function uv_response_measure

   !> dE/dt_s at fixed potential tables; dcdt = dc_s/dt_s at fixed u_sr.
   pure elemental function uv_response_rhs(a, dh, dc, dcdt) result(g)
      real(wp), intent(in) :: a, dh, dc, dcdt
      real(wp) :: g
      g = a*(dh*(1.0_wp + dcdt) + dc*dcdt)
   end function uv_response_rhs

   !> Adjoint of u_sr; dcdu = dc_s/du_sr at fixed t_s.
   !> mu = T'^T lambda is the contribution from the equation adjoint.
   pure elemental function uv_response_sr(a, dh, dc, mu, dcdu) result(w)
      real(wp), intent(in) :: a, dh, dc, mu, dcdu
      real(wp) :: w
      w = dcdu*(a*(dh + dc) + mu)
   end function uv_response_sr

   !> Adjoint of u_lr(r), which enters the energy through physical c
   pure elemental function uv_response_lr(a, dc, beta) result(w)
      real(wp), intent(in) :: a, dc, beta
      real(wp) :: w
      w = -beta*a*dc
   end function uv_response_lr

end module moist_model_moz_response_kernel
