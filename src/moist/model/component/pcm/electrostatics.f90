!> Nuclear gradient of the PCM electrostatic surface coupling
!>
!> The electrostatic coupling energy `E = sum_i q_i phi(r_i)` depends on the
!> nuclei twice: through the potential's own sources at fixed surface, and
!> through the surface positions `r_i(R)`. The host owns the potential, so the
!> second route is the host total `w_xyz_i = dE_host/dr_i = q_i grad phi(r_i)`
!> (nuclear plus electronic); MOIST forms only the first, the nuclear charges
!> moving under the fixed surface charges.
module moist_model_component_pcm_electrostatics
   use mctc_env, only: wp, error_type, fatal_error
   implicit none(type, external)
   private

   public :: pcm_electrostatic_nuclear_gradient
   public :: pcm_electrostatic_direct_gradient

contains

   !> Nuclear gradient of the electrostatic surface coupling at fixed surface
   !>
   !> The solute nuclei move under the fixed surface charges:
   !>
   !> `grad_rA = sum_i q_i Z_A (r_i - R_A)/|r_i - R_A|^3`. This is the only
   !> nuclear term MOIST forms itself; the surface-position term is the host
   !> total `w_xyz_i = dE_host/dr_i`, which the cavity contracts in reverse
   !> mode or [[pcm_electrostatic_nuclear_gradient]] contracts in forward mode.
   !>
   !> @param[in]  xyz        Surface positions (3, ngrid)
   !> @param[in]  sphxyz     Atomic sphere centers (3, nsph)
   !> @param[in]  surface_q  Surface charges (ngrid)
   !> @param[in]  za         Nuclear charges (nsph)
   !> @param[out] grad_rA    Direct nuclear gradient at fixed surface (3, nsph)
   !> @param[in] xi         Optional Gaussian inverse lengths (ngrid)
   !> @param[out] error      Error handling
   subroutine pcm_electrostatic_direct_gradient(xyz, sphxyz, surface_q, za, grad_rA, error, xi)
      !> Surface positions and sphere centers
      real(wp), intent(in) :: xyz(:, :)
      !> Sphxyz.
      real(wp), intent(in) :: sphxyz(:, :)
      !> Surface charges and nuclear charges
      real(wp), intent(in) :: surface_q(:)
      !> Za.
      real(wp), intent(in) :: za(:)
      !> Direct nuclear gradient
      real(wp), intent(out) :: grad_rA(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Gaussian inverse lengths; absent selects the point operator.
      real(wp), intent(in), optional :: xi(:)

      !> Surface, source-atom, and extent indices
      integer :: i, katom, ngrid, nsph
      !> Surface charge and displacement data
      real(wp) :: qi, rvec(3), r2, inv_r3, x
      !> Squared-distance threshold for coincident sources
      real(wp), parameter :: r2tol = 1.0e-30_wp

      grad_rA = 0.0_wp
      ngrid = size(surface_q)
      nsph = size(za)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(grad_rA, 1) /= 3 .or. size(grad_rA, 2) /= nsph) then
         call fatal_error(error, "pcm_electrostatic_direct_gradient: array shape mismatch")
         return
      end if

      if (present(xi)) then
         if (size(xi) /= ngrid) then
            call fatal_error(error, "pcm_electrostatic_direct_gradient: width shape mismatch")
            return
         end if
      end if

      ! Serial: the work is O(ngrid*nsph) with a handful of flops per pair and
      ! the reduction target is tiny, so an OpenMP array reduction here would
      ! cost more than it saves
      do i = 1, ngrid
         qi = surface_q(i)
         do katom = 1, nsph
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            if (r2 <= r2tol) cycle
            inv_r3 = 1.0_wp/(sqrt(r2)*r2)
            if (present(xi)) then
               x = xi(i)*sqrt(r2)
               if (abs(x) < 1.0e-3_wp) then
                  inv_r3 = 4.0_wp*xi(i)*xi(i)*xi(i)/(3.0_wp*sqrt(acos(-1.0_wp))) &
                          & *(1.0_wp - 0.6_wp*x*x + 3.0_wp*x*x*x*x/14.0_wp)
               else
                  inv_r3 = inv_r3*(erf(x) - 2.0_wp*x*exp(-x*x)/sqrt(acos(-1.0_wp)))
               end if
            end if
            grad_rA(:, katom) = grad_rA(:, katom) + qi*za(katom)*inv_r3*rvec
         end do
      end do

   end subroutine pcm_electrostatic_direct_gradient

   !> Forward-mode nuclear gradient of the electrostatic surface coupling
   !>
   !> The direct term of [[pcm_electrostatic_direct_gradient]] plus the host
   !> total surface-position weight contracted with the surface's response to
   !> nuclear motion: `grad_rA += sum_i xyz1_rA(:, :, A, i)^T w_xyz(:, i)`.
   !>
   !> Used by cavities that expose `xyz1_rA`; the reverse-mode path hands
   !> `w_xyz` to the cavity instead.
   !>
   !> @param[in]  xyz         Surface positions (3, ngrid)
   !> @param[in]  sphxyz      Atomic sphere centers (3, nsph)
   !> @param[in]  xyz1_rA     Surface-position derivatives (3, 3, nsph, ngrid)
   !> @param[in]  surface_q   Surface charges (ngrid)
   !> @param[in]  w_xyz       Host total surface-position weight (3, ngrid)
   !> @param[in]  za          Nuclear charges (nsph)
   !> @param[out] grad_rA     Nuclear gradient (3, nsph)
   !> @param[in] xi          Optional Gaussian inverse lengths (ngrid)
   !> @param[out] error       Error handling
   subroutine pcm_electrostatic_nuclear_gradient(xyz, sphxyz, xyz1_rA, &
                                                 surface_q, w_xyz, za, grad_rA, error, xi)
      !> Surface positions, sphere centers, and surface-position derivatives
      real(wp), intent(in) :: xyz(:, :)
      !> Sphxyz.
      real(wp), intent(in) :: sphxyz(:, :)
      !> Xyz1 ra.
      real(wp), intent(in) :: xyz1_rA(:, :, :, :)
      !> Surface charges, host surface-position weights, and nuclear charges
      real(wp), intent(in) :: surface_q(:)
      !> W xyz.
      real(wp), intent(in) :: w_xyz(:, :)
      !> Za.
      real(wp), intent(in) :: za(:)
      !> Contracted nuclear gradient
      real(wp), intent(out) :: grad_rA(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Gaussian inverse lengths; absent selects the point operator.
      real(wp), intent(in), optional :: xi(:)

      !> Surface, moving-atom, source-atom, and extent indices
      integer :: i, iatom, ngrid, nsph
      !> Tmp reduction target; explicit shape
      real(wp) :: acc(3, size(za))

      grad_rA = 0.0_wp
      ngrid = size(surface_q)
      nsph = size(za)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(xyz1_rA, 1) /= 3 .or. size(xyz1_rA, 2) /= 3 .or. &
          size(xyz1_rA, 3) /= nsph .or. size(xyz1_rA, 4) /= ngrid .or. &
          size(w_xyz, 1) /= 3 .or. size(w_xyz, 2) /= ngrid .or. &
          size(grad_rA, 1) /= 3 .or. size(grad_rA, 2) /= nsph) then
         call fatal_error(error, "pcm_electrostatic_nuclear_gradient: array shape mismatch")
         return
      end if

      call pcm_electrostatic_direct_gradient(xyz, sphxyz, surface_q, za, grad_rA, error, xi)
      if (allocated(error)) return
      acc = 0.0_wp
      !$omp parallel do default(none) &
      !$omp shared(xyz1_rA, w_xyz, ngrid, nsph) &
      !$omp private(i, iatom) &
      !$omp reduction(+:acc) schedule(static)
      do i = 1, ngrid
         do iatom = 1, nsph
            acc(1, iatom) = acc(1, iatom) &
                            + dot_product(xyz1_rA(:, 1, iatom, i), w_xyz(:, i))
            acc(2, iatom) = acc(2, iatom) &
                            + dot_product(xyz1_rA(:, 2, iatom, i), w_xyz(:, i))
            acc(3, iatom) = acc(3, iatom) &
                            + dot_product(xyz1_rA(:, 3, iatom, i), w_xyz(:, i))
         end do
      end do
      !$omp end parallel do

      grad_rA(:, :) = grad_rA + acc
   end subroutine pcm_electrostatic_nuclear_gradient

end module moist_model_component_pcm_electrostatics
