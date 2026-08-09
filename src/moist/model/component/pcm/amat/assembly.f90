!> Assembly of the Gaussian PCM interaction matrix and dense derivatives.
submodule(moist_model_component_pcm_amat) moist_model_component_pcm_amat_assembly
   use mctc_env, only: fatal_error
   use moist_model_component_pcm_amat_kernel, only: pcm_amat_x_far, &
      pcm_amat_far_value_row, pcm_amat_near_value, pcm_amat_near_grad, &
      pcm_amat_diag_value, pcm_amat_diag_grad
   implicit none(type, external)

   !> Cache-blocking tile used when mirroring the assembled triangle
   integer, parameter :: sym_tile = 64

contains

   !> Validate the common Gaussian surface arrays
   !>
   !> @param[in]  xi     Gaussian widths
   !> @param[in]  f      Gaussian switching factors
   !> @param[in]  xyz    Surface positions
   !> @param[out] error  Error handling
   module subroutine validate_pcm_surface(xi, f, xyz, error)
      !> Gaussian widths
      real(wp), intent(in) :: xi(:)
      !> Gaussian switching factors
      real(wp), intent(in) :: f(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (size(xi) <= 0) then
         call fatal_error(error, "Gaussian PCM surface is empty")
      else if (size(f) /= size(xi) .or. size(xyz, 1) /= 3 .or. &
               size(xyz, 2) /= size(xi)) then
         call fatal_error(error, "Gaussian PCM surface arrays have inconsistent shapes")
      else if (any(xi <= 0.0_wp)) then
         call fatal_error(error, "Gaussian PCM widths must be positive")
      else if (any(f <= 0.0_wp)) then
         call fatal_error(error, "Gaussian PCM switching factors must be positive")
      end if
   end subroutine validate_pcm_surface

   !> Compute per-point saturation bounds for the near/far test
   !>
   !> @param[in]  xi     Gaussian widths
   !> @param[out] bound  Per-point saturation bounds
   module pure subroutine saturation_bounds(xi, bound)
      !> Gaussian widths
      real(wp), intent(in) :: xi(:)
      !> Per-point saturation bounds
      real(wp), intent(out) :: bound(:)

      bound = pcm_amat_x_far/(xi*xi)
   end subroutine saturation_bounds

   !> Mirror the upper triangle of a square matrix onto the lower one.
   !> @param[inout] amat Matrix whose upper triangle is filled on entry
   subroutine mirror_upper_triangle(amat)
      !> Matrix whose upper triangle is filled on entry
      real(wp), intent(inout) :: amat(:, :)

      !> Block indices, point ranges, and block count
      integer :: ib, jb, i0, i1, j0, j1, nblocks
      !> Point indices and total number of points
      integer :: i, j, ngrid

      ngrid = size(amat, 1)
      nblocks = (ngrid + sym_tile - 1)/sym_tile

      !$omp parallel do default(none) shared(amat, ngrid, nblocks) &
      !$omp private(ib, jb, i0, i1, j0, j1, i, j) schedule(dynamic)
      do ib = 1, nblocks
         i0 = (ib - 1)*sym_tile + 1
         i1 = min(ib*sym_tile, ngrid)
         do jb = 1, ib
            j0 = (jb - 1)*sym_tile + 1
            j1 = min(jb*sym_tile, ngrid)
            do j = j0, min(j1, i1 - 1)
               do i = max(i0, j + 1), i1
                  amat(i, j) = amat(j, i)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine mirror_upper_triangle

   !> Assemble the Gaussian PCM interaction matrix
   !>
   !> @param[in]  xi     Gaussian widths
   !> @param[in]  f      Gaussian switching factors
   !> @param[in]  xyz    Surface positions
   !> @param[out] amat   Interaction matrix
   !> @param[out] error  Error handling
   module subroutine assemble_pcm_amat(xi, f, xyz, amat, error)
      !> Gaussian widths
      real(wp), intent(in) :: xi(:)
      !> Gaussian switching factors
      real(wp), intent(in) :: f(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Interaction matrix
      real(wp), intent(out) :: amat(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface-point indices and total number of points
      integer :: i, j, ngrid
      !> Row-point width, position, and saturation bound
      real(wp) :: xi_i, xyz_i(3), bound_i
      !> Per-point saturation bounds
      real(wp), allocatable :: bound(:)

      ! Every matrix element is written on the success path. Avoiding an
      ! initial zeroing pass preserves first-touch performance.
      call validate_pcm_surface(xi, f, xyz, error)
      if (allocated(error)) then
         amat = 0.0_wp
         return
      end if

      ngrid = size(xi)
      if (size(amat, 1) /= ngrid .or. size(amat, 2) /= ngrid) then
         amat = 0.0_wp
         call fatal_error(error, "assemble_pcm_amat: matrix shape mismatch")
         return
      end if

      allocate (bound(ngrid))
      call saturation_bounds(xi, bound)

      !$omp parallel default(none) shared(xi, f, xyz, amat, bound, ngrid) &
      !$omp private(i, j, xi_i, xyz_i, bound_i)
      block
         !> Per-thread squared-separation scratch for one matrix row
         real(wp), allocatable :: r2(:)

         allocate (r2(ngrid))
         !$omp do schedule(dynamic, 8)
         do i = 1, ngrid
            call pcm_amat_diag_value(xi(i), f(i), amat(i, i))

            xi_i = xi(i)
            xyz_i = xyz(:, i)
            bound_i = bound(i)

            do j = 1, i - 1
               r2(j) = max((xyz_i(1) - xyz(1, j))**2 + (xyz_i(2) - xyz(2, j))**2 &
                           + (xyz_i(3) - xyz(3, j))**2, r2_floor)
            end do

            call pcm_amat_far_value_row(i - 1, r2, amat(1:i - 1, i))

            do j = 1, i - 1
               if (r2(j) >= bound_i + bound(j)) cycle
               call pcm_amat_near_value(xi_i, xi(j), r2(j), amat(j, i))
            end do
         end do
         !$omp end do
      end block
      !$omp end parallel

      call mirror_upper_triangle(amat)
   end subroutine assemble_pcm_amat

   !> Assemble the Gaussian PCM matrix and its nuclear derivative tensor
   !>
   !> @param[in]  xi       Gaussian widths
   !> @param[in]  f        Gaussian switching factors
   !> @param[in]  xyz      Surface positions
   !> @param[in]  xi1_rA   Nuclear derivatives of widths
   !> @param[in]  f1_rA    Nuclear derivatives of switching factors
   !> @param[in]  xyz1_rA  Nuclear derivatives of surface positions
   !> @param[out] amat     Interaction matrix
   !> @param[out] amat1_rA Nuclear derivative tensor
   !> @param[out] error    Error handling
   module subroutine assemble_pcm_amat_with_gradient(xi, f, xyz, xi1_rA, f1_rA, &
                                                      xyz1_rA, amat, amat1_rA, error)
      !> Gaussian widths
      real(wp), intent(in) :: xi(:)
      !> Gaussian switching factors
      real(wp), intent(in) :: f(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Nuclear derivatives of Gaussian widths
      real(wp), intent(in) :: xi1_rA(:, :, :)
      !> Nuclear derivatives of switching factors
      real(wp), intent(in) :: f1_rA(:, :, :)
      !> Nuclear derivatives of surface positions
      real(wp), intent(in) :: xyz1_rA(:, :, :, :)
      !> Interaction matrix
      real(wp), intent(out) :: amat(:, :)
      !> Nuclear derivative tensor
      real(wp), intent(out) :: amat1_rA(:, :, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface, atom, Cartesian, and extent indices
      integer :: i, j, iatom, iaxis, ngrid, nsph
      !> Kernel value and its three derivative channels
      real(wp) :: a, a_xi_i, a_xi_j, a_r2
      !> Diagonal kernel value and its two derivative channels
      real(wp) :: a_diag, a_diag_xi, a_diag_f
      !> Pair displacement, squared separation, and nuclear response
      real(wp) :: rvec(3), r2, dr2
      !> Per-point saturation bounds
      real(wp), allocatable :: bound(:)

      amat = 0.0_wp
      amat1_rA = 0.0_wp
      call validate_pcm_surface(xi, f, xyz, error)
      if (allocated(error)) return
      ngrid = size(xi)
      if (size(xi1_rA, 1) /= 3 .or. size(f1_rA, 1) /= 3 .or. &
          size(xyz1_rA, 1) /= 3 .or. size(xyz1_rA, 2) /= 3 .or. &
          size(xi1_rA, 2) /= size(f1_rA, 2) .or. &
          size(xi1_rA, 2) /= size(xyz1_rA, 3) .or. &
          size(xi1_rA, 3) /= ngrid .or. size(f1_rA, 3) /= ngrid .or. &
          size(xyz1_rA, 4) /= ngrid) then
         call fatal_error(error, "assemble_pcm_amat_with_gradient: derivative shape mismatch")
         return
      end if
      nsph = size(xi1_rA, 2)
      if (size(amat, 1) /= ngrid .or. size(amat, 2) /= ngrid) then
         call fatal_error(error, "assemble_pcm_amat_with_gradient: matrix shape mismatch")
         return
      end if
      if (size(amat1_rA, 1) /= 3 .or. size(amat1_rA, 2) /= nsph .or. &
          size(amat1_rA, 3) /= ngrid .or. size(amat1_rA, 4) /= ngrid) then
         call fatal_error(error, "assemble_pcm_amat_with_gradient: output shape mismatch")
         return
      end if

      allocate (bound(ngrid))
      call saturation_bounds(xi, bound)

      !$omp parallel do default(none) &
      !$omp shared(xi, f, xyz, xi1_rA, f1_rA, xyz1_rA, amat, amat1_rA, bound, ngrid, nsph) &
      !$omp private(i, j, iatom, iaxis, a, a_xi_i, a_xi_j, a_r2, a_diag, a_diag_xi, &
      !$omp         a_diag_f, rvec, r2, dr2) schedule(dynamic, 8)
      do i = 1, ngrid
         call pcm_amat_diag_grad(xi(i), f(i), a_diag, a_diag_xi, a_diag_f)
         amat(i, i) = a_diag
         do iatom = 1, nsph
            do iaxis = 1, 3
               amat1_rA(iaxis, iatom, i, i) = a_diag_xi*xi1_rA(iaxis, iatom, i) &
                                              + a_diag_f*f1_rA(iaxis, iatom, i)
            end do
         end do

         do j = 1, i - 1
            rvec = xyz(:, i) - xyz(:, j)
            r2 = max(sum(rvec*rvec), r2_floor)
            if (r2 >= bound(i) + bound(j)) then
               a = 1.0_wp/sqrt(r2)
               a_xi_i = 0.0_wp
               a_xi_j = 0.0_wp
               a_r2 = -0.5_wp/(r2*sqrt(r2))
            else
               call pcm_amat_near_grad(xi(i), xi(j), r2, a, a_xi_i, a_xi_j, a_r2)
            end if
            amat(i, j) = a
            amat(j, i) = a

            do iatom = 1, nsph
               do iaxis = 1, 3
                  dr2 = 2.0_wp*dot_product(rvec, xyz1_rA(:, iaxis, iatom, i) &
                                           - xyz1_rA(:, iaxis, iatom, j))
                  amat1_rA(iaxis, iatom, i, j) = a_xi_i*xi1_rA(iaxis, iatom, i) &
                                                 + a_xi_j*xi1_rA(iaxis, iatom, j) &
                                                 + a_r2*dr2
                  amat1_rA(iaxis, iatom, j, i) = amat1_rA(iaxis, iatom, i, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine assemble_pcm_amat_with_gradient

end submodule moist_model_component_pcm_amat_assembly
