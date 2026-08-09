!> Nuclear gradient of the PCM electrostatic surface coupling
module moist_model_component_pcm_electrostatics
   use mctc_env, only: wp, error_type, fatal_error
   !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
   implicit none (type, external)
   private

   public :: pcm_electrostatic_nuclear_gradient

contains

   !> Contract direct nuclear and electronic surface-field contributions
   !>
   !> @param[in]  xyz         Surface positions
   !> @param[in]  sphxyz      Atomic sphere centers
   !> @param[in]  xyz1_rA     Surface-position derivatives
   !> @param[in]  surface_q   Surface charges
   !> @param[in]  qefield     Charge-weighted electronic field
   !> @param[in]  za          Nuclear charges
   !> @param[out] grad_rA     Nuclear gradient
   !> @param[out] error       Error handling
   subroutine pcm_electrostatic_nuclear_gradient(xyz, sphxyz, xyz1_rA, &
      surface_q, qefield, za, grad_rA, error)
      !> Surface positions, sphere centers, and surface-position derivatives
      real(wp), intent(in) :: xyz(:, :), sphxyz(:, :), xyz1_rA(:, :, :, :)
      !> Surface charges, electronic field weights, and nuclear charges
      real(wp), intent(in) :: surface_q(:), qefield(:, :), za(:)
      !> Contracted nuclear gradient
      real(wp), intent(out) :: grad_rA(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface, moving-atom, source-atom, and extent indices
      integer :: i, iatom, katom, ngrid, nsph
      !> Surface charge, displacement data, nuclear field, and chain-rule field
      real(wp) :: qi, rvec(3), r2, inv_r3, enuc(3), chain(3)
      !> Squared-distance threshold for coincident sources
      real(wp), parameter :: r2tol = 1.0e-30_wp

      grad_rA = 0.0_wp
      ngrid = size(surface_q)
      nsph = size(za)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(xyz1_rA, 1) /= 3 .or. size(xyz1_rA, 2) /= 3 .or. &
          size(xyz1_rA, 3) /= nsph .or. size(xyz1_rA, 4) /= ngrid .or. &
          size(qefield, 1) /= 3 .or. size(qefield, 2) /= ngrid .or. &
          size(grad_rA, 1) /= 3 .or. size(grad_rA, 2) /= nsph) then
         call fatal_error(error, "pcm_electrostatic_nuclear_gradient: array shape mismatch")
         return
      end if

      ! TODO: Reverse-mode (z-vector) migration
      !
      ! Using xyz1_rA(3, 3, nsph, ngrid) scales quadratically and thus this does not parallelize
      ! well >4 threads
      !
      ! The chain vector below is the surface adjoint w_xyz, and the direct nuclear term is
      ! independent of that; we can thus use w_xyz = qefield - q_i*E_nuc via
      ! cavity_surface_adjoint_type (as pcm_base_get_surface_weights already does) and let the
      ! cavity contract it via contract_surface_lsf_weights which solves the same per-point
      !bordered KKT system with 4 RHS instead of 3*n_active (i.e. not forming xyz1_rA explicitly)

      !$omp parallel do default(none) reduction(+:grad_rA) &
      !$omp shared(xyz, sphxyz, xyz1_rA, surface_q, qefield, za, ngrid, nsph) &
      !$omp private(i, iatom, katom, qi, rvec, r2, inv_r3, enuc, chain) schedule(static)
      do i = 1, ngrid
         qi = surface_q(i)
         enuc = 0.0_wp
         do katom = 1, nsph
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            if (r2 <= r2tol) cycle
            inv_r3 = 1.0_wp/(sqrt(r2)*r2)
            enuc = enuc + za(katom)*inv_r3*rvec
            grad_rA(:, katom) = grad_rA(:, katom) + qi*za(katom)*inv_r3*rvec
         end do
         chain = qefield(:, i) - qi*enuc
         do iatom = 1, nsph
            grad_rA(1, iatom) = grad_rA(1, iatom) &
                                + dot_product(xyz1_rA(:, 1, iatom, i), chain)
            grad_rA(2, iatom) = grad_rA(2, iatom) &
                                + dot_product(xyz1_rA(:, 2, iatom, i), chain)
            grad_rA(3, iatom) = grad_rA(3, iatom) &
                                + dot_product(xyz1_rA(:, 3, iatom, i), chain)
         end do
      end do
      !$omp end parallel do
   end subroutine pcm_electrostatic_nuclear_gradient

end module moist_model_component_pcm_electrostatics
