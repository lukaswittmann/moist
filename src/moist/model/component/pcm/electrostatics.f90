!> Nuclear gradient of the PCM electrostatic surface coupling
module moist_model_component_pcm_electrostatics
   use mctc_env, only: wp, error_type, fatal_error
   implicit none (type, external)
   private

   public :: pcm_electrostatic_nuclear_gradient
   public :: pcm_electrostatic_surface_weights
   public :: pcm_electrostatic_potential_tangent
   public :: pcm_electrostatic_surface_weights_response

contains

   !> Split the electrostatic surface coupling into adjoint and direct parts
   !>
   !> The reverse-mode counterpart of [[pcm_electrostatic_nuclear_gradient]].
   !> That routine contracts the surface-position adjoint `qefield - q_i*E_nuc`
   !> with `xyz1_rA` on the spot; here the same vector is handed back as a
   !> surface weight for the cavity to contract, and only the term that does
   !> *not* flow through the surface -- the nuclei moving under the fixed
   !> surface charges -- is returned as a gradient.
   !>
   !> Deliberately kept as a separate routine rather than shared with the
   !> forward version: folding them together would reorder the floating-point
   !> accumulation and perturb the legacy path.
   !>
   !> @param[in]  xyz        Surface positions (3, ngrid)
   !> @param[in]  sphxyz     Atomic sphere centers (3, nsph)
   !> @param[in]  surface_q  Surface charges (ngrid)
   !> @param[in]  qefield    Charge-weighted electronic field (3, ngrid)
   !> @param[in]  za         Nuclear charges (nsph)
   !> @param[out] w_xyz      Surface-position adjoint (3, ngrid)
   !> @param[out] grad_rA    Direct nuclear gradient at fixed surface (3, nsph)
   !> @param[out] error      Error handling
   subroutine pcm_electrostatic_surface_weights(xyz, sphxyz, surface_q, qefield, za, &
                                                w_xyz, grad_rA, error)
      !> Surface positions and sphere centers
      real(wp), intent(in) :: xyz(:, :), sphxyz(:, :)
      !> Surface charges, electronic field weights, and nuclear charges
      real(wp), intent(in) :: surface_q(:), qefield(:, :), za(:)
      !> Surface-position adjoint
      real(wp), intent(out) :: w_xyz(:, :)
      !> Direct nuclear gradient
      real(wp), intent(out) :: grad_rA(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface, source-atom, and extent indices
      integer :: i, katom, ngrid, nsph
      !> Surface charge, displacement data, and nuclear field
      real(wp) :: qi, rvec(3), r2, inv_r3, enuc(3)
      !> Squared-distance threshold for coincident sources
      real(wp), parameter :: r2tol = 1.0e-30_wp

      grad_rA = 0.0_wp
      w_xyz = 0.0_wp
      ngrid = size(surface_q)
      nsph = size(za)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(qefield, 1) /= 3 .or. size(qefield, 2) /= ngrid .or. &
          size(w_xyz, 1) /= 3 .or. size(w_xyz, 2) /= ngrid .or. &
          size(grad_rA, 1) /= 3 .or. size(grad_rA, 2) /= nsph) then
         call fatal_error(error, "pcm_electrostatic_surface_weights: array shape mismatch")
         return
      end if

      ! Serial: the work is O(ngrid*nsph) with a handful of flops per pair and
      ! the reduction target is tiny, so an OpenMP array reduction here would
      ! cost more than it saves
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
         w_xyz(:, i) = qefield(:, i) - qi*enuc
      end do

   end subroutine pcm_electrostatic_surface_weights

   !> Directional derivative of the point-charge potential at the surface
   !>
   !> For `phi_i = sum_k s_k / |r_i - R_k|` the potential moves along a
   !> direction on which the surface points move by `d_i` and the nuclei by
   !> `v_k` as
   !>
   !>     dphi_i = sum_k s_k g_ik . (v_k - d_i),   g_ik = (r_i - R_k)/|r_i - R_k|^3
   !>
   !> The coincidence threshold is the one of the gradient, so a source that
   !> the gradient skips contributes no response either.
   !>
   !> @param[in]  xyz    Surface positions (3, ngrid)
   !> @param[in]  sphxyz Atomic sphere centers (3, nsph)
   !> @param[in]  za     Source charges (nsph)
   !> @param[in]  d_xyz  Surface-position tangents (3, ngrid, ndir)
   !> @param[in]  dirs   Nuclear directions (3, nsph, ndir)
   !> @param[out] dphi   Potential response per direction (ngrid, ndir)
   !> @param[out] error  Error handling
   subroutine pcm_electrostatic_potential_tangent(xyz, sphxyz, za, d_xyz, dirs, dphi, error)
      !> Surface positions and sphere centers
      real(wp), intent(in) :: xyz(:, :), sphxyz(:, :)
      !> Source charges
      real(wp), intent(in) :: za(:)
      !> Surface-position tangents
      real(wp), intent(in) :: d_xyz(:, :, :)
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Potential response per direction
      real(wp), intent(out) :: dphi(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface, source-atom, direction and extent indices
      integer :: i, katom, idir, ngrid, nsph, ndir
      !> Displacement data and the scaled field of one source
      real(wp) :: rvec(3), r2, inv_r3, g(3)
      !> Squared-distance threshold for coincident sources
      real(wp), parameter :: r2tol = 1.0e-30_wp

      dphi = 0.0_wp
      ngrid = size(xyz, 2)
      nsph = size(za)
      ndir = size(dirs, 3)
      if (size(xyz, 1) /= 3 .or. size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(d_xyz, 1) /= 3 .or. size(d_xyz, 2) /= ngrid .or. size(d_xyz, 3) /= ndir .or. &
          size(dirs, 1) /= 3 .or. size(dirs, 2) /= nsph .or. &
          size(dphi, 1) /= ngrid .or. size(dphi, 2) /= ndir) then
         call fatal_error(error, "pcm_electrostatic_potential_tangent: array shape mismatch")
         return
      end if

      !$omp parallel do default(none) shared(xyz, sphxyz, za, d_xyz, dirs, dphi, ngrid, nsph, ndir) &
      !$omp private(i, katom, idir, rvec, r2, inv_r3, g) schedule(static)
      do i = 1, ngrid
         do katom = 1, nsph
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            if (r2 <= r2tol) cycle
            inv_r3 = 1.0_wp/(sqrt(r2)*r2)
            g = za(katom)*inv_r3*rvec
            do idir = 1, ndir
               dphi(i, idir) = dphi(i, idir) &
                               + dot_product(g, dirs(:, katom, idir) - d_xyz(:, i, idir))
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine pcm_electrostatic_potential_tangent

   !> Response of the electrostatic surface adjoint and of the direct gradient
   !>
   !> The second-order counterpart of [[pcm_electrostatic_surface_weights]] for
   !> the point-charge sources: along a direction on which the surface points
   !> move by `d_i`, the nuclei by `v_k` and the surface charges by `dq_i`, the
   !> surface-position adjoint `-q_i E_i` and the direct gradient
   !> `sum_i q_i s_k g_ik` respond by
   !>
   !>     dw_i    = -dq_i E_i - q_i sum_k s_k T_ik (d_i - v_k)
   !>     dg_k    = sum_i [ dq_i s_k g_ik + q_i s_k T_ik (d_i - v_k) ]
   !>
   !> with `T_ik = dg_ik/dr_i = I/|r|^3 - 3 r r^T/|r|^5`. The host's electronic
   !> field is not part of this: it is host data, held fixed here.
   !>
   !> @param[in]  xyz        Surface positions (3, ngrid)
   !> @param[in]  sphxyz     Atomic sphere centers (3, nsph)
   !> @param[in]  surface_q  Surface charges (ngrid)
   !> @param[in]  dq         Surface-charge response per direction (ngrid, ndir)
   !> @param[in]  za         Source charges (nsph)
   !> @param[in]  d_xyz      Surface-position tangents (3, ngrid, ndir)
   !> @param[in]  dirs       Nuclear directions (3, nsph, ndir)
   !> @param[out] dw_xyz     Surface-position adjoint response (3, ngrid, ndir)
   !> @param[out] dgrad_rA   Direct nuclear-gradient response (3, nsph, ndir)
   !> @param[out] error      Error handling
   subroutine pcm_electrostatic_surface_weights_response(xyz, sphxyz, surface_q, dq, za, &
                                                         d_xyz, dirs, dw_xyz, dgrad_rA, error)
      !> Surface positions and sphere centers
      real(wp), intent(in) :: xyz(:, :), sphxyz(:, :)
      !> Surface charges and their response per direction
      real(wp), intent(in) :: surface_q(:), dq(:, :)
      !> Source charges
      real(wp), intent(in) :: za(:)
      !> Surface-position tangents
      real(wp), intent(in) :: d_xyz(:, :, :)
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Surface-position adjoint response
      real(wp), intent(out) :: dw_xyz(:, :, :)
      !> Direct nuclear-gradient response
      real(wp), intent(out) :: dgrad_rA(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface, source-atom, direction and extent indices
      integer :: i, katom, idir, ngrid, nsph, ndir
      !> Surface charge, displacement data, and the field of one source
      real(wp) :: qi, rvec(3), r2, inv_r3, inv_r5, g(3), enuc(3)
      !> Relative motion of the pair and the field response `T (d - v)`
      real(wp) :: rel(3), tv(3)
      !> Squared-distance threshold for coincident sources
      real(wp), parameter :: r2tol = 1.0e-30_wp

      dw_xyz = 0.0_wp
      dgrad_rA = 0.0_wp
      ngrid = size(surface_q)
      nsph = size(za)
      ndir = size(dirs, 3)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(dq, 1) /= ngrid .or. size(dq, 2) /= ndir .or. &
          size(d_xyz, 1) /= 3 .or. size(d_xyz, 2) /= ngrid .or. size(d_xyz, 3) /= ndir .or. &
          size(dirs, 1) /= 3 .or. size(dirs, 2) /= nsph .or. &
          size(dw_xyz, 1) /= 3 .or. size(dw_xyz, 2) /= ngrid .or. size(dw_xyz, 3) /= ndir .or. &
          size(dgrad_rA, 1) /= 3 .or. size(dgrad_rA, 2) /= nsph .or. &
          size(dgrad_rA, 3) /= ndir) then
         call fatal_error(error, "pcm_electrostatic_surface_weights_response: array shape mismatch")
         return
      end if

      ! Surface-position adjoint: every point owns its row, no reduction
      !$omp parallel do default(none) &
      !$omp shared(xyz, sphxyz, surface_q, dq, za, d_xyz, dirs, dw_xyz, ngrid, nsph, ndir) &
      !$omp private(i, katom, idir, qi, rvec, r2, inv_r3, inv_r5, g, enuc, rel, tv) schedule(static)
      do i = 1, ngrid
         qi = surface_q(i)
         enuc = 0.0_wp
         do katom = 1, nsph
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            if (r2 <= r2tol) cycle
            inv_r3 = 1.0_wp/(sqrt(r2)*r2)
            inv_r5 = inv_r3/r2
            enuc = enuc + za(katom)*inv_r3*rvec
            do idir = 1, ndir
               rel = d_xyz(:, i, idir) - dirs(:, katom, idir)
               tv = inv_r3*rel - 3.0_wp*inv_r5*dot_product(rvec, rel)*rvec
               dw_xyz(:, i, idir) = dw_xyz(:, i, idir) - qi*za(katom)*tv
            end do
         end do
         do idir = 1, ndir
            dw_xyz(:, i, idir) = dw_xyz(:, i, idir) - dq(i, idir)*enuc
         end do
      end do
      !$omp end parallel do

      ! Direct term: every source owns its column, no reduction
      !$omp parallel do default(none) &
      !$omp shared(xyz, sphxyz, surface_q, dq, za, d_xyz, dirs, dgrad_rA, ngrid, nsph, ndir) &
      !$omp private(i, katom, idir, qi, rvec, r2, inv_r3, inv_r5, g, rel, tv) schedule(static)
      do katom = 1, nsph
         do i = 1, ngrid
            qi = surface_q(i)
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            if (r2 <= r2tol) cycle
            inv_r3 = 1.0_wp/(sqrt(r2)*r2)
            inv_r5 = inv_r3/r2
            g = za(katom)*inv_r3*rvec
            do idir = 1, ndir
               rel = d_xyz(:, i, idir) - dirs(:, katom, idir)
               tv = inv_r3*rel - 3.0_wp*inv_r5*dot_product(rvec, rel)*rvec
               dgrad_rA(:, katom, idir) = dgrad_rA(:, katom, idir) &
                                          + dq(i, idir)*g + qi*za(katom)*tv
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine pcm_electrostatic_surface_weights_response

   !> Contract direct nuclear and electronic surface-field contributions
   !>
   !> Should not be used; is the legacy forward path implementation
   !> of [[pcm_electrostatic_nuclear_gradient]]
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
      !> Tmp reduction target; explicit shape
      real(wp) :: acc(3, size(za))

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

      acc = 0.0_wp
      !$omp parallel do default(none) &
      !$omp shared(xyz, sphxyz, xyz1_rA, surface_q, qefield, za, ngrid, nsph) &
      !$omp private(i, iatom, katom, qi, rvec, r2, inv_r3, enuc, chain) &
      !$omp reduction(+:acc) schedule(static)
      do i = 1, ngrid
         qi = surface_q(i)
         enuc = 0.0_wp
         do katom = 1, nsph
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            if (r2 <= r2tol) cycle
            inv_r3 = 1.0_wp/(sqrt(r2)*r2)
            enuc = enuc + za(katom)*inv_r3*rvec
            acc(:, katom) = acc(:, katom) + qi*za(katom)*inv_r3*rvec
         end do
         chain = qefield(:, i) - qi*enuc
         do iatom = 1, nsph
            acc(1, iatom) = acc(1, iatom) &
                            + dot_product(xyz1_rA(:, 1, iatom, i), chain)
            acc(2, iatom) = acc(2, iatom) &
                            + dot_product(xyz1_rA(:, 2, iatom, i), chain)
            acc(3, iatom) = acc(3, iatom) &
                            + dot_product(xyz1_rA(:, 3, iatom, i), chain)
         end do
      end do
      !$omp end parallel do

      grad_rA(:, :) = acc
   end subroutine pcm_electrostatic_nuclear_gradient

end module moist_model_component_pcm_electrostatics
