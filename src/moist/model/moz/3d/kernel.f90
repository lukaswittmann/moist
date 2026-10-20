!> Per-k kernels of the 3D UV equation
!>
!> - One k-point of the solvent-site columns, h_a(k) = sum_g c_g(k) P_ga(|k|)
!> - P(g, a): solvent response w_ga + rho_g h_ga of the VV solution at |k|
!> - Ng split: c(k) = c_s(k) - beta u_lr(k)
!>   - t_s(k) = h(k) - c(k) - beta u_lr(k)
!>
!> - Products follow the 3D UV equation of moist_dev/pysofft (44874e7e),
!>   model/rism/3d/uv/equation.F90
module moist_model_moz_3d_kernel
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: uv3d_point_map, uv3d_point_jvp, uv3d_point_adjoint, uv3d_tail_adjoint

contains

   !> t_s(k) from c_s(k) at one k-point
   !>
   !> - Optional OZ route projects uk onto charge and uses its stable chi*q/k^2
   !> - Physical Ng tables have uk = amplitude*charge; custom tables must satisfy this too
   !>
   !> @param[in]  p                response P(g, a), (nv, nv)
   !> @param[in]  cs               c_s(k) of each solvent site, (nv)
   !> @param[in]  uk               long-range potential u_lr(k) of each site, Hartree bohr^3, (nv)
   !> @param[in]  beta             inverse thermal energy, 1/Hartree
   !> @param[out] ts               t_s(k), (nv)
   !> @param[in]  charge           optional solvent site charges of the OZ route, e (nv)
   !> @param[in]  charge_response  optional stable chi*q/k^2 of the OZ route, (nv)
   !> @param[in]  k2               optional |k|^2, 1/bohr^2
   pure subroutine uv3d_point_map(p, cs, uk, beta, ts, charge, charge_response, k2)
      !> Response
      real(wp), intent(in) :: p(:, :)
      !> Short-range direct correlation and long-range potential
      complex(wp), intent(in) :: cs(:), uk(:)
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Short-range indirect correlation
      complex(wp), intent(out) :: ts(:)
      !> Optional OZ route: site charges (nv), stable chi q / k^2 (nv), |k|^2
      real(wp), intent(in), optional :: charge(:), charge_response(:), k2
      complex(wp) :: amplitude
      real(wp) :: q2
      integer :: a

      ! t_s = c_s chi - c_s - beta u_lr chi
      ! Keep short-range terms separate from u_lr, whose entries grow as 1/k^2
      call uv3d_point_jvp(p, cs, ts)
      if (present(charge) .and. present(charge_response) .and. present(k2)) then
         q2 = sum(charge**2)
         if (q2 > 0.0_wp) then
            amplitude = sum(charge*uk)/q2
            ts = ts - beta*(amplitude*k2)*charge_response
            return
         end if
      end if
      do a = 1, size(cs)
         ts(a) = ts(a) - beta*sum(uk*p(:, a))
      end do
   end subroutine uv3d_point_map

   !> Derivative of the map at one k-point: dt(k) = sum_g dc_g P_ga - dc_a
   !>
   !> @param[in]  p    response P(g, a), (nv, nv)
   !> @param[in]  dc   direction dc_s(k), (nv)
   !> @param[out] dt   dt_s(k), (nv)
   pure subroutine uv3d_point_jvp(p, dc, dt)
      !> Response
      real(wp), intent(in) :: p(:, :)
      !> Direction
      complex(wp), intent(in) :: dc(:)
      !> Response of the map
      complex(wp), intent(out) :: dt(:)
      integer :: a

      do a = 1, size(dc)
         dt(a) = sum(dc*p(:, a)) - dc(a)
      end do
   end subroutine uv3d_point_jvp

   !> Transpose of dc_s -> dt_s under Re(sum(conjg(seed)*direction))
   !>
   !> @param[in]  p     response P(g, a), (nv, nv)
   !> @param[in]  seed  adjoint of dt_s(k), (nv)
   !> @param[out] wc    adjoint of dc_s(k), (nv)
   pure subroutine uv3d_point_adjoint(p, seed, wc)
      !> Response
      real(wp), intent(in) :: p(:, :)
      !> Adjoint of the map output
      complex(wp), intent(in) :: seed(:)
      !> Adjoint of the direction
      complex(wp), intent(out) :: wc(:)
      integer :: g

      do g = 1, size(seed)
         wc(g) = sum(p(g, :)*seed) - seed(g)
      end do
   end subroutine uv3d_point_adjoint

   !> Transpose of the long-range potential contribution, including -beta
   !>
   !> - Stable OZ route projects uk onto the solvent charge vector
   !> - Its transpose must differentiate that projection, including the zero mode
   !>
   !> @param[in]  p                response P(g, a), (nv, nv)
   !> @param[in]  seed             adjoint of t_s(k), (nv)
   !> @param[in]  beta             inverse thermal energy, 1/Hartree
   !> @param[out] wk               adjoint of u_lr(k), (nv)
   !> @param[in]  charge           optional solvent site charges of the OZ route, e (nv)
   !> @param[in]  charge_response  optional stable chi*q/k^2 of the OZ route, (nv)
   !> @param[in]  k2               optional |k|^2, 1/bohr^2
   pure subroutine uv3d_tail_adjoint(p, seed, beta, wk, charge, charge_response, k2)
      !> Response, and inverse thermal energy
      real(wp), intent(in) :: p(:, :), beta
      !> Adjoint of the map output
      complex(wp), intent(in) :: seed(:)
      !> Adjoint of the long-range potential
      complex(wp), intent(out) :: wk(:)
      !> Optional OZ route: site charges (nv), stable chi q / k^2 (nv), |k|^2
      real(wp), intent(in), optional :: charge(:), charge_response(:), k2
      real(wp) :: q2
      integer :: g

      if (present(charge) .and. present(charge_response) .and. present(k2)) then
         q2 = sum(charge**2)
         if (q2 > 0.0_wp) then
            wk = (-beta*k2/q2)*charge*sum(charge_response*seed)
            return
         end if
      end if
      do g = 1, size(seed)
         wk(g) = -beta*sum(p(g, :)*seed)
      end do
   end subroutine uv3d_tail_adjoint

end module moist_model_moz_3d_kernel
