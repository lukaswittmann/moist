!> Nuclear gradient of the PCM electrostatic surface coupling
module moist_model_component_pcm_electrostatics
   use mctc_env, only: wp, error_type, fatal_error
   !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
   implicit none (type, external)
   private

   public :: pcm_electrostatic_nuclear_gradient
   public :: pcm_electrostatic_surface_weights

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
      !> Thread identity and count for the manual reduction
      integer :: ithread, nthreads
      !> Per-thread gradient accumulator and the buffer collecting them
      real(wp), allocatable :: acc(:, :), grad_buf(:, :, :)

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

      ! gfortran miscompiles reduction(+:) on an assumed-shape array dummy, so
      ! the reduction is done by hand: each thread sums into a private
      ! accumulator, parks it in its own slice, and the slices are added up
      ! serially. That also makes the result independent of the thread count.
      nthreads = 1
      !$ nthreads = omp_get_max_threads()
      allocate (grad_buf(3, nsph, nthreads), source=0.0_wp)

      !$omp parallel default(none) &
      !$omp shared(xyz, sphxyz, xyz1_rA, surface_q, qefield, za, ngrid, nsph, grad_buf) &
      !$omp private(i, iatom, katom, qi, rvec, r2, inv_r3, enuc, chain, ithread, acc)
      ithread = 1
      !$ ithread = omp_get_thread_num() + 1
      allocate (acc(3, nsph), source=0.0_wp)
      !$omp do schedule(static)
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
      !$omp end do
      grad_buf(:, :, ithread) = acc
      !$omp end parallel

      do ithread = 1, nthreads
         grad_rA(:, :) = grad_rA(:, :) + grad_buf(:, :, ithread)
      end do
   end subroutine pcm_electrostatic_nuclear_gradient

end module moist_model_component_pcm_electrostatics
