!> Adjoint contractions for the Gaussian PCM interaction matrix
submodule(moist_model_component_pcm_amat) moist_model_component_pcm_amat_adjoint
   use mctc_env, only: fatal_error
   use moist_model_component_pcm_amat_kernel, only: pcm_amat_near_grad, pcm_amat_diag_grad
   !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
   implicit none (type, external)

   !> Threshold below which a charge product contributes no useful precision
   real(wp), parameter :: qtol = 1.0e-30_wp

contains

   !> Contract a Gaussian PCM matrix derivative to surface-variable weights
   !>
   !> @param[in]  xi     Gaussian widths
   !> @param[in]  f      Gaussian switching factors
   !> @param[in]  xyz    Surface positions
   !> @param[in]  q1     Left contraction vector
   !> @param[in]  q2     Right contraction vector
   !> @param[out] w_xi   Width weights
   !> @param[out] w_f    Switching-factor weights
   !> @param[out] w_xyz  Position weights
   !> @param[out] error  Error handling
   module subroutine pcm_amat_surface_weights(xi, f, xyz, q1, q2, w_xi, w_f, &
                                              w_xyz, error)
      !> Gaussian widths
      real(wp), intent(in) :: xi(:)
      !> Gaussian switching factors
      real(wp), intent(in) :: f(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Left and right contraction vectors
      real(wp), intent(in) :: q1(:), q2(:)
      !> Gaussian-width and switching-factor weights
      real(wp), intent(out) :: w_xi(:), w_f(:)
      !> Surface-position weights
      real(wp), intent(out) :: w_xyz(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface-point indices and total number of points
      integer :: i, j, ngrid
      !> Row-point properties and contraction-vector entries
      real(wp) :: xi_i, xyz_i(3), bound_i, q1i, q2i
      !> Charge product, masked separation, and position factor
      real(wp) :: qsym, r2s, scale
      !> Kernel value and its three derivative channels
      real(wp) :: a, a_xi_i, a_xi_j, a_r2
      !> Diagonal kernel value and its two derivative channels
      real(wp) :: a_diag, a_diag_xi, a_diag_f
      !> Accumulated position weight for the row point
      real(wp) :: acc(3)
      !> Saturation flag for the current pair
      logical :: is_far
      !> Per-point saturation bounds
      real(wp), allocatable :: bound(:)
      !> Per-thread squared-separation scratch for one matrix row
      real(wp), allocatable :: r2(:)

      w_xi = 0.0_wp
      w_f = 0.0_wp
      w_xyz = 0.0_wp
      call validate_pcm_surface(xi, f, xyz, error)
      if (allocated(error)) return
      ngrid = size(xi)
      if (size(q1) /= ngrid .or. size(q2) /= ngrid .or. &
          size(w_xi) /= ngrid .or. size(w_f) /= ngrid .or. &
          size(w_xyz, 1) /= 3 .or. size(w_xyz, 2) /= ngrid) then
         call fatal_error(error, "pcm_amat_surface_weights: array shape mismatch")
         return
      end if

      allocate (bound(ngrid))
      call saturation_bounds(xi, bound)

      ! The scratch row is a private allocatable: each thread allocates its own
      ! copy inside the region. A block-local declaration would be cleaner, but
      ! ifx rejects it under default(none).
      !$omp parallel default(none) &
      !$omp shared(xi, f, xyz, q1, q2, w_xi, w_f, w_xyz, bound, ngrid) &
      !$omp private(i, j, xi_i, xyz_i, bound_i, q1i, q2i, qsym, r2s, scale, r2, &
      !$omp         a, a_xi_i, a_xi_j, a_r2, a_diag, a_diag_xi, a_diag_f, acc, is_far)
      allocate (r2(ngrid))
      !$omp do schedule(dynamic, 8)
      do i = 1, ngrid
         xi_i = xi(i)
         xyz_i = xyz(:, i)
         bound_i = bound(i)
         q1i = q1(i)
         q2i = q2(i)

         call pcm_amat_diag_grad(xi_i, f(i), a_diag, a_diag_xi, a_diag_f)
         w_xi(i) = q1i*q2i*a_diag_xi
         w_f(i) = q1i*q2i*a_diag_f

         do j = 1, ngrid
            r2(j) = max((xyz_i(1) - xyz(1, j))**2 + (xyz_i(2) - xyz(2, j))**2 &
                        + (xyz_i(3) - xyz(3, j))**2, r2_floor)
         end do

         ! The saturated pass has no width channel. Near pairs and the self
         ! term are masked and handled below.
         acc = 0.0_wp
         do j = 1, ngrid
            is_far = r2(j) >= bound_i + bound(j)
            qsym = q1i*q2(j) + q1(j)*q2i
            r2s = merge(r2(j), 1.0_wp, is_far)
            scale = merge(-qsym/(r2s*sqrt(r2s)), 0.0_wp, is_far)
            acc(1) = acc(1) + scale*(xyz_i(1) - xyz(1, j))
            acc(2) = acc(2) + scale*(xyz_i(2) - xyz(2, j))
            acc(3) = acc(3) + scale*(xyz_i(3) - xyz(3, j))
         end do

         do j = 1, ngrid
            if (j == i) cycle
            if (r2(j) >= bound_i + bound(j)) cycle
            qsym = q1i*q2(j) + q1(j)*q2i
            if (abs(qsym) <= qtol) cycle
            call pcm_amat_near_grad(xi_i, xi(j), r2(j), a, a_xi_i, a_xi_j, a_r2)
            w_xi(i) = w_xi(i) + qsym*a_xi_i
            scale = 2.0_wp*qsym*a_r2
            acc(1) = acc(1) + scale*(xyz_i(1) - xyz(1, j))
            acc(2) = acc(2) + scale*(xyz_i(2) - xyz(2, j))
            acc(3) = acc(3) + scale*(xyz_i(3) - xyz(3, j))
         end do

         w_xyz(:, i) = acc
      end do
      !$omp end do
      !$omp end parallel
   end subroutine pcm_amat_surface_weights

   !> Contract surface-variable weights with nuclear derivative arrays
   !>
   !> Should not be used; is the legacy forward path implementation
   !> of [[pcm_amat_nuclear_gradient]]
   !>
   !> @param[in]  xi1_rA      Width derivatives
   !> @param[in]  f1_rA       Switching-factor derivatives
   !> @param[in]  xyz1_rA     Surface-position derivatives
   !> @param[in]  w_xi        Width weights
   !> @param[in]  w_f         Switching-factor weights
   !> @param[in]  w_xyz       Position weights
   !> @param[out] grad_rA     Nuclear gradient
   !> @param[out] error       Error handling
   module subroutine pcm_amat_nuclear_gradient(xi1_rA, f1_rA, xyz1_rA, w_xi, &
                                                w_f, w_xyz, grad_rA, error)
      !> Nuclear derivatives of widths, switching factors, and positions
      real(wp), intent(in) :: xi1_rA(:, :, :), f1_rA(:, :, :), xyz1_rA(:, :, :, :)
      !> Surface-variable adjoint weights
      real(wp), intent(in) :: w_xi(:), w_f(:), w_xyz(:, :)
      !> Contracted nuclear gradient
      real(wp), intent(out) :: grad_rA(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface, atom, Cartesian, and extent indices
      integer :: i, iatom, iaxis, ngrid, nsph
      !> Thread identity and count for the manual reduction
      integer :: ithread, nthreads
      !> Per-thread gradient accumulator and the buffer collecting them
      real(wp), allocatable :: acc(:, :), grad_buf(:, :, :)

      grad_rA = 0.0_wp
      ngrid = size(w_xi)
      nsph = size(xi1_rA, 2)
      if (size(xi1_rA, 1) /= 3 .or. size(f1_rA, 1) /= 3 .or. &
          size(xyz1_rA, 1) /= 3 .or. size(xyz1_rA, 2) /= 3 .or. &
          size(f1_rA, 2) /= nsph .or. size(xyz1_rA, 3) /= nsph .or. &
          size(xi1_rA, 3) /= ngrid .or. size(f1_rA, 3) /= ngrid .or. &
          size(xyz1_rA, 4) /= ngrid .or. size(w_f) /= ngrid .or. &
          size(w_xyz, 1) /= 3 .or. size(w_xyz, 2) /= ngrid .or. &
          size(grad_rA, 1) /= 3 .or. size(grad_rA, 2) /= nsph) then
         call fatal_error(error, "pcm_amat_nuclear_gradient: array shape mismatch")
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
      !$omp shared(xi1_rA, f1_rA, xyz1_rA, w_xi, w_f, w_xyz, ngrid, nsph, grad_buf) &
      !$omp private(i, iatom, iaxis, ithread, acc)
      ithread = 1
      !$ ithread = omp_get_thread_num() + 1
      allocate (acc(3, nsph), source=0.0_wp)
      !$omp do schedule(static)
      do i = 1, ngrid
         do iatom = 1, nsph
            do iaxis = 1, 3
               acc(iaxis, iatom) = acc(iaxis, iatom) &
                                   + w_xi(i)*xi1_rA(iaxis, iatom, i) &
                                   + w_f(i)*f1_rA(iaxis, iatom, i) &
                                   + dot_product( &
                                      xyz1_rA(:, iaxis, iatom, i), w_xyz(:, i))
            end do
         end do
      end do
      !$omp end do
      grad_buf(:, :, ithread) = acc
      !$omp end parallel

      do ithread = 1, nthreads
         grad_rA(:, :) = grad_rA(:, :) + grad_buf(:, :, ithread)
      end do
   end subroutine pcm_amat_nuclear_gradient

end submodule moist_model_component_pcm_amat_adjoint
