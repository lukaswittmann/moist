!> Pointwise UV response kernels, independent of grid dimension and layout
!>
!> - h = c_s + t_s, c = c_s - beta u_lr, c_s = C(u_sr, t_s)
!> - Closure supplies dh = dPhi/dh (fixed c) and dc = dPhi/dc (fixed h)
!> - Phi is the dimensionless free-energy integrand
!> - Callers apply energy seeds, may add the moment correction seed beta*dE/dD to dc
!> - Derivatives exclude density and quadrature weights
!> - a = rho_site*quadrature_weight/beta is the point's energy weight
module moist_model_moz_response_kernel
   use mctc_env, only: wp
   implicit none(type, external)
   private
   public :: uv_response_rhs, uv_response_sr, uv_response_lr, uv_response_measure

contains

   !> Explicit quadrature sensitivity of the selected energy and direct moment
   !>
   !> - phi is the dimensionless closure integrand
   !> - shift = beta*dE/dD
   !> - rho/beta has energy/volume units
   !>
   !> @param[in] rho Site number density
   !> @param[in] beta Inverse thermal energy
   !> @param[in] energy_seed Seed weight of the energy term
   !> @param[in] phi Closure integrand, dimensionless
   !> @param[in] shift Moment correction seed beta*dE/dD
   !> @param[in] c Physical direct correlation, dimensionless
   pure elemental function uv_response_measure(rho, beta, energy_seed, phi, shift, c) result(w)
      !> Site number density, inverse thermal energy, energy seed, closure integrand, moment shift and direct correlation
      real(wp), intent(in) :: rho, beta, energy_seed, phi, shift, c
      !> Quadrature sensitivity (rho/beta)*(energy_seed*phi + shift*c)
      real(wp) :: w
      w = (rho/beta)*(energy_seed*phi + shift*c)
   end function uv_response_measure

   !> dE/dt_s at fixed potential tables
   !>
   !> @param[in] a Energy weight rho_site*quadrature_weight/beta of the point
   !> @param[in] dh Closure derivative dPhi/dh at fixed c
   !> @param[in] dc Closure derivative dPhi/dc at fixed h
   !> @param[in] dcdt dc_s/dt_s at fixed u_sr
   pure elemental function uv_response_rhs(a, dh, dc, dcdt) result(g)
      !> Energy weight, closure derivatives dh and dc, and dc_s/dt_s
      real(wp), intent(in) :: a, dh, dc, dcdt
      !> Energy derivative dE/dt_s
      real(wp) :: g
      g = a*(dh*(1.0_wp + dcdt) + dc*dcdt)
   end function uv_response_rhs

   !> Adjoint of u_sr
   !>
   !> - mu = T'^T lambda is the contribution from the equation adjoint
   !>
   !> @param[in] a Energy weight rho_site*quadrature_weight/beta of the point
   !> @param[in] dh Closure derivative dPhi/dh at fixed c
   !> @param[in] dc Closure derivative dPhi/dc at fixed h
   !> @param[in] mu Equation-adjoint contribution T'^T lambda
   !> @param[in] dcdu dc_s/du_sr at fixed t_s
   pure elemental function uv_response_sr(a, dh, dc, mu, dcdu) result(w)
      !> Energy weight, closure derivatives dh and dc, equation adjoint mu and dc_s/du_sr
      real(wp), intent(in) :: a, dh, dc, mu, dcdu
      !> Adjoint of u_sr
      real(wp) :: w
      w = dcdu*(a*(dh + dc) + mu)
   end function uv_response_sr

   !> Adjoint of u_lr(r), which enters the energy through physical c
   !>
   !> @param[in] a Energy weight rho_site*quadrature_weight/beta of the point
   !> @param[in] dc Closure derivative dPhi/dc at fixed h
   !> @param[in] beta Inverse thermal energy
   pure elemental function uv_response_lr(a, dc, beta) result(w)
      !> Energy weight, closure derivative dc and inverse thermal energy
      real(wp), intent(in) :: a, dc, beta
      !> Adjoint of u_lr
      real(wp) :: w
      w = -beta*a*dc
   end function uv_response_lr

end module moist_model_moz_response_kernel
