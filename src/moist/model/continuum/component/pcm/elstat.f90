!> Nuclear gradient of the PCM electrostatic surface coupling
!>
!> The electrostatic coupling energy depends on the nuclei through the host
!> potential twice: through the potential's own sources at fixed surface, and
!> through the surface positions `r_i(R)`. Both routes are weighted by the
!> potential adjoint `w_phi_i = dE/dphi_i`, which is the surface charge `q_i`
!> only for stationary potential-only models and the symmetrized adjoint of a
!> non-symmetric one; normal-trace models add the `dE/dvn` channel. The
!> host owns the potential, so the second route is the host total
!> `w_xyz_i = dE_host/dr_i = w_phi_i grad phi(r_i)` (nuclear plus electronic);
!> MOIST forms only the first, the nuclear charges moving under the fixed
!> surface
module moist_model_continuum_component_pcm_elstat
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: resolve_num_threads
   implicit none(type, external)
   private

   public :: pcm_elstat_nuclear_gradient
   public :: pcm_elstat_direct_gradient
   public :: pcm_elstat_normal_direct_gradient
   public :: pcm_elstat_direct_hessian
   public :: pcm_elstat_normal_direct_hessian
   public :: pcm_elstat_potential_derivatives
   public :: pcm_elstat_potential_tangent

   !> Packed symmetric `(3, 3)` index pairs `[xx, xy, xz, yy, yz, zz]`
   integer, parameter :: sym2_idx(2, 6) = reshape([ &
      & 1, 1, 1, 2, 1, 3, 2, 2, 2, 3, 3, 3], [2, 6])

contains

   !> Fixed-surface nuclear derivative of a potential spatial-gradient adjoint
   !>
   !> The Gaussian branch differentiates Z*erf(xi*r)/r, matching the potential
   !> used by the host. Omitting xi selects the point-charge potential.
   subroutine pcm_elstat_normal_direct_gradient(xyz, sphxyz, w_grad, za, grad_rA, error, xi)
      real(wp), intent(in) :: xyz(:, :)
      real(wp), intent(in) :: sphxyz(:, :)
      real(wp), intent(in) :: w_grad(:, :)
      real(wp), intent(in) :: za(:)
      real(wp), intent(out) :: grad_rA(:, :)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: xi(:)
      integer :: i, katom, ngrid, nsph
      real(wp) :: rvec(3), r2, r, x, aterm, c0, c1, projection, sqrt_pi

      grad_rA = 0.0_wp
      ngrid = size(w_grad, 2)
      nsph = size(za)
      if (size(w_grad, 1) /= 3 .or. size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(grad_rA, 1) /= 3 .or. size(grad_rA, 2) /= nsph) then
         call fatal_error(error, "pcm_elstat_normal_direct_gradient: array shape mismatch")
         return
      end if
      if (present(xi)) then
         if (size(xi) /= ngrid) then
            call fatal_error(error, "pcm_elstat_normal_direct_gradient: width shape mismatch")
            return
         end if
      end if

      sqrt_pi = sqrt(acos(-1.0_wp))
      do i = 1, ngrid
         do katom = 1, nsph
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            r = sqrt(r2)
            if (present(xi)) then
               x = xi(i)*r
               if (abs(x) < 1.0e-3_wp) then
                  c0 = 4.0_wp*xi(i)**3/(3.0_wp*sqrt_pi) &
                     & *(1.0_wp - 0.6_wp*x*x + 3.0_wp*x**4/14.0_wp)
                  c1 = -8.0_wp*xi(i)**5/(5.0_wp*sqrt_pi) &
                     & + 8.0_wp*xi(i)**7*r2/(7.0_wp*sqrt_pi)
               else
                  aterm = erf(x) - 2.0_wp*x*exp(-x*x)/sqrt_pi
                  c0 = aterm/(r2*r)
                  c1 = 4.0_wp*xi(i)**3*exp(-x*x)/(sqrt_pi*r2) &
                     & - 3.0_wp*aterm/(r2*r2*r)
               end if
            else
               if (r2 <= 1.0e-30_wp) cycle
               c0 = 1.0_wp/(r2*r)
               c1 = -3.0_wp/(r2*r2*r)
            end if
            projection = dot_product(rvec, w_grad(:, i))
            grad_rA(:, katom) = grad_rA(:, katom) + za(katom) &
               & *(c0*w_grad(:, i) + c1*projection*rvec)
         end do
      end do
   end subroutine pcm_elstat_normal_direct_gradient

   !> Nuclear gradient of the electrostatic surface coupling at fixed surface
   !>
   !> The solute nuclei move under the fixed surface:
   !>
   !> `grad_rA = sum_i w_phi_i Z_A (r_i - R_A)/|r_i - R_A|^3`. This is the only
   !> nuclear term MOIST forms itself; the surface-position term is the host
   !> total `w_xyz_i = dE_host/dr_i`, which the cavity contracts in reverse
   !> mode or [[pcm_elstat_nuclear_gradient]] contracts in forward mode
   !>
   !> @param[in]  xyz        Surface positions (3, ngrid)
   !> @param[in]  sphxyz     Atomic sphere centers (3, nsph)
   !> @param[in]  w_phi      Potential adjoint (ngrid)
   !> @param[in]  za         Nuclear charges (nsph)
   !> @param[out] grad_rA    Direct nuclear gradient at fixed surface (3, nsph)
   !> @param[in] xi         Optional Gaussian inverse lengths (ngrid)
   !> @param[out] error      Error handling
   subroutine pcm_elstat_direct_gradient(xyz, sphxyz, w_phi, za, grad_rA, error, xi)
      !> Surface positions and sphere centers
      real(wp), intent(in) :: xyz(:, :)
      !> Sphxyz
      real(wp), intent(in) :: sphxyz(:, :)
      !> Potential adjoint and nuclear charges
      real(wp), intent(in) :: w_phi(:)
      !> Za
      real(wp), intent(in) :: za(:)
      !> Direct nuclear gradient
      real(wp), intent(out) :: grad_rA(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Gaussian inverse lengths; absent selects the point operator
      real(wp), intent(in), optional :: xi(:)

      !> Surface, source-atom, and extent indices
      integer :: i, katom, ngrid, nsph
      !> Adjoint weight and displacement data
      real(wp) :: qi, rvec(3), r2, inv_r3, x
      !> Squared-distance threshold for coincident sources
      real(wp), parameter :: r2tol = 1.0e-30_wp

      grad_rA = 0.0_wp
      ngrid = size(w_phi)
      nsph = size(za)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(grad_rA, 1) /= 3 .or. size(grad_rA, 2) /= nsph) then
         call fatal_error(error, "pcm_elstat_direct_gradient: array shape mismatch")
         return
      end if

      if (present(xi)) then
         if (size(xi) /= ngrid) then
            call fatal_error(error, "pcm_elstat_direct_gradient: width shape mismatch")
            return
         end if
      end if

      ! Serial: the work is O(ngrid*nsph) with a handful of flops per pair and
      ! the reduction target is tiny, so an OpenMP array reduction here would
      ! cost more than it saves
      do i = 1, ngrid
         qi = w_phi(i)
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

   end subroutine pcm_elstat_direct_gradient

   !> Forward-mode nuclear gradient of the electrostatic surface coupling
   !>
   !> The direct term of [[pcm_elstat_direct_gradient]] plus the host
   !> total surface-position weight contracted with the surface's response to
   !> nuclear motion: `grad_rA += sum_i xyz1_rA(:, :, A, i)^T w_xyz(:, i)`
   !>
   !> Used by cavities that expose `xyz1_rA`; the reverse-mode path hands
   !> `w_xyz` to the cavity instead
   !>
   !> @param[in]  xyz         Surface positions (3, ngrid)
   !> @param[in]  sphxyz      Atomic sphere centers (3, nsph)
   !> @param[in]  xyz1_rA     Surface-position derivatives (3, 3, nsph, ngrid)
   !> @param[in]  w_phi       Potential adjoint (ngrid)
   !> @param[in]  w_xyz       Host total surface-position weight (3, ngrid)
   !> @param[in]  za          Nuclear charges (nsph)
   !> @param[out] grad_rA     Nuclear gradient (3, nsph)
   !> @param[in] xi          Optional Gaussian inverse lengths (ngrid)
   !> @param[out] error       Error handling
   !> @param[in]  nthreads    OpenMP team size; absent takes omp_get_max_threads
   subroutine pcm_elstat_nuclear_gradient(xyz, sphxyz, xyz1_rA, &
                                          w_phi, w_xyz, za, grad_rA, error, xi, nthreads)
      !> Surface positions, sphere centers, and surface-position derivatives
      real(wp), intent(in) :: xyz(:, :)
      !> Sphxyz
      real(wp), intent(in) :: sphxyz(:, :)
      !> Xyz1 ra
      real(wp), intent(in) :: xyz1_rA(:, :, :, :)
      !> Potential adjoint, host surface-position weights, and nuclear charges
      real(wp), intent(in) :: w_phi(:)
      !> W xyz
      real(wp), intent(in) :: w_xyz(:, :)
      !> Za
      real(wp), intent(in) :: za(:)
      !> Contracted nuclear gradient
      real(wp), intent(out) :: grad_rA(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Gaussian inverse lengths; absent selects the point operator
      real(wp), intent(in), optional :: xi(:)
      !> OpenMP team size
      integer, intent(in), optional :: nthreads

      !> Surface, moving-atom, source-atom, and extent indices, and team size
      integer :: i, iatom, ngrid, nsph, nt
      !> Tmp reduction target; explicit shape
      real(wp) :: acc(3, size(za))

      grad_rA = 0.0_wp
      ngrid = size(w_phi)
      nsph = size(za)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(xyz1_rA, 1) /= 3 .or. size(xyz1_rA, 2) /= 3 .or. &
          size(xyz1_rA, 3) /= nsph .or. size(xyz1_rA, 4) /= ngrid .or. &
          size(w_xyz, 1) /= 3 .or. size(w_xyz, 2) /= ngrid .or. &
          size(grad_rA, 1) /= 3 .or. size(grad_rA, 2) /= nsph) then
         call fatal_error(error, "pcm_elstat_nuclear_gradient: array shape mismatch")
         return
      end if

      call pcm_elstat_direct_gradient(xyz, sphxyz, w_phi, za, grad_rA, error, xi)
      if (allocated(error)) return
      acc = 0.0_wp
      nt = resolve_num_threads(nthreads)
      !$omp parallel do num_threads(nt) default(none) &
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
   end subroutine pcm_elstat_nuclear_gradient

   !> Radial kernel coefficients `c0`, `c1`, `c2` of `psi(r, xi) = erf(xi*r)/r` or `1/r`
   pure subroutine elstat_radial_kernels(r2, has_xi, xi_i, c0, c1, c2, exp_x2)
      real(wp), intent(in) :: r2
      logical, intent(in) :: has_xi
      real(wp), intent(in) :: xi_i
      real(wp), intent(out) :: c0, c1, c2, exp_x2

      real(wp) :: r, x, x2, aterm, sqrt_pi, xi2, xi3, xi5, xi7

      sqrt_pi = sqrt(acos(-1.0_wp))
      c0 = 0.0_wp
      c1 = 0.0_wp
      c2 = 0.0_wp
      exp_x2 = 0.0_wp
      if (has_xi) then
         r = sqrt(r2)
         x = xi_i*r
         x2 = x*x
         exp_x2 = exp(-x2)
         xi2 = xi_i*xi_i
         xi3 = xi2*xi_i
         xi5 = xi3*xi2
         xi7 = xi5*xi2
         if (abs(x) < 0.15_wp) then
            c0 = (4.0_wp*xi3/sqrt_pi)*(1.0_wp/3.0_wp - x2/5.0_wp + x2*x2/14.0_wp &
               & - x2**3/54.0_wp + x2**4/264.0_wp - x2**5/1560.0_wp + x2**6/10800.0_wp)
            c1 = -(8.0_wp*xi5/sqrt_pi)*(1.0_wp/5.0_wp - x2/7.0_wp + x2*x2/18.0_wp &
               & - x2**3/66.0_wp + x2**4/312.0_wp - x2**5/1800.0_wp + x2**6/12240.0_wp)
            c2 = (16.0_wp*xi7/sqrt_pi)*(1.0_wp/7.0_wp - x2/9.0_wp + x2*x2/22.0_wp &
               & - x2**3/78.0_wp + x2**4/360.0_wp - x2**5/2040.0_wp + x2**6/13680.0_wp)
         else
            aterm = erf(x) - 2.0_wp*x*exp_x2/sqrt_pi
            c0 = aterm/(r2*r)
            c1 = 4.0_wp*xi3*exp_x2/(sqrt_pi*r2) - 3.0_wp*aterm/(r2*r2*r)
            c2 = 15.0_wp*aterm/(r2*r2*r2*r) - 4.0_wp*xi3*exp_x2*(5.0_wp + 2.0_wp*x2)/(sqrt_pi*r2*r2)
         end if
      else
         if (r2 <= 1.0e-30_wp) return
         r = sqrt(r2)
         c0 = 1.0_wp/(r2*r)
         c1 = -3.0_wp/(r2*r2*r)
         c2 = 15.0_wp/(r2*r2*r2*r)
      end if

   end subroutine elstat_radial_kernels

   !> Directional derivative of [[pcm_elstat_direct_gradient]] along moving surface & nuclei
   !>
   !> @param[in]     xyz       surface positions `(3, ngrid)`
   !> @param[in]     sphxyz    atomic sphere centers `(3, nsph)`
   !> @param[in]     w_phi     potential adjoint `(ngrid)`
   !> @param[in]     dw_phi    potential-adjoint tangents `(ngrid, ndir)`
   !> @param[in]     za        nuclear charges `(nsph)`
   !> @param[in]     d_xyz     surface position tangents `(3, ngrid, ndir)`
   !> @param[in]     dirs      nuclear directions `(3, nsph, ndir)`
   !> @param[in,out] hvp_rA    nuclear HVP accumulator `(3, nsph, ndir)`
   !> @param[out]    error     error handling
   !> @param[in]     xi        Gaussian inverse lengths `(ngrid)`, optional
   !> @param[in]     d_xi      Gaussian width tangents `(ngrid, ndir)`, optional
   subroutine pcm_elstat_direct_hessian(xyz, sphxyz, w_phi, dw_phi, za, d_xyz, dirs, &
                                        hvp_rA, error, xi, d_xi)
      real(wp), intent(in) :: xyz(:, :), sphxyz(:, :)
      real(wp), intent(in) :: w_phi(:), dw_phi(:, :), za(:)
      real(wp), intent(in) :: d_xyz(:, :, :), dirs(:, :, :)
      real(wp), intent(inout) :: hvp_rA(:, :, :)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: xi(:), d_xi(:, :)

      integer :: i, katom, idir, ngrid, nsph, ndir
      logical :: has_xi
      real(wp) :: rvec(3), uvec(3), r2, c0, c1, c2, exp_x2, xi_i, dxi_i, sqrt_pi, dc0_dxi

      ngrid = size(w_phi)
      nsph = size(za)
      ndir = size(dirs, 3)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(dw_phi, 1) /= ngrid .or. size(dw_phi, 2) /= ndir .or. &
          size(d_xyz, 1) /= 3 .or. size(d_xyz, 2) /= ngrid .or. size(d_xyz, 3) /= ndir .or. &
          size(dirs, 1) /= 3 .or. size(dirs, 2) /= nsph .or. &
          size(hvp_rA, 1) /= 3 .or. size(hvp_rA, 2) /= nsph .or. size(hvp_rA, 3) /= ndir) then
         call fatal_error(error, "pcm_elstat_direct_hessian: array shape mismatch")
         return
      end if
      has_xi = present(xi)
      if (has_xi) then
         if (size(xi) /= ngrid) then
            call fatal_error(error, "pcm_elstat_direct_hessian: width shape mismatch")
            return
         end if
         if (present(d_xi)) then
            if (size(d_xi, 1) /= ngrid .or. size(d_xi, 2) /= ndir) then
               call fatal_error(error, "pcm_elstat_direct_hessian: width-tangent shape mismatch")
               return
            end if
         end if
      end if

      sqrt_pi = sqrt(acos(-1.0_wp))
      do katom = 1, nsph
         do i = 1, ngrid
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            if (.not. has_xi .and. r2 <= 1.0e-30_wp) cycle
            xi_i = 0.0_wp
            if (has_xi) xi_i = xi(i)
            call elstat_radial_kernels(r2, has_xi, xi_i, c0, c1, c2, exp_x2)
            dc0_dxi = 0.0_wp
            if (has_xi) dc0_dxi = 4.0_wp*xi_i*xi_i*exp_x2/sqrt_pi
            do idir = 1, ndir
               uvec = d_xyz(:, i, idir) - dirs(:, katom, idir)
               dxi_i = 0.0_wp
               if (has_xi .and. present(d_xi)) dxi_i = d_xi(i, idir)
               hvp_rA(:, katom, idir) = hvp_rA(:, katom, idir) + za(katom)*( &
                  & dw_phi(i, idir)*c0*rvec + w_phi(i)*( &
                  & c0*uvec + c1*dot_product(rvec, uvec)*rvec + dxi_i*dc0_dxi*rvec))
            end do
         end do
      end do

   end subroutine pcm_elstat_direct_hessian

   !> Directional derivative of [[pcm_elstat_normal_direct_gradient]] along moving surface & nuclei
   !>
   !> @param[in]     xyz       surface positions `(3, ngrid)`
   !> @param[in]     sphxyz    atomic sphere centers `(3, nsph)`
   !> @param[in]     w_grad    normal-trace spatial weights `w_vn * n` `(3, ngrid)`
   !> @param[in]     dw_grad   tangents `dw_vn * n + w_vn * d_n` `(3, ngrid, ndir)`
   !> @param[in]     za        nuclear charges `(nsph)`
   !> @param[in]     d_xyz     surface position tangents `(3, ngrid, ndir)`
   !> @param[in]     dirs      nuclear directions `(3, nsph, ndir)`
   !> @param[in,out] hvp_rA    nuclear HVP accumulator `(3, nsph, ndir)`
   !> @param[out]    error     error handling
   !> @param[in]     xi        Gaussian inverse lengths `(ngrid)`, optional
   !> @param[in]     d_xi      Gaussian width tangents `(ngrid, ndir)`, optional
   subroutine pcm_elstat_normal_direct_hessian(xyz, sphxyz, w_grad, dw_grad, za, d_xyz, dirs, &
                                               hvp_rA, error, xi, d_xi)
      real(wp), intent(in) :: xyz(:, :), sphxyz(:, :)
      real(wp), intent(in) :: w_grad(:, :), dw_grad(:, :, :), za(:)
      real(wp), intent(in) :: d_xyz(:, :, :), dirs(:, :, :)
      real(wp), intent(inout) :: hvp_rA(:, :, :)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: xi(:), d_xi(:, :)

      integer :: i, katom, idir, ngrid, nsph, ndir
      logical :: has_xi
      real(wp) :: rvec(3), uvec(3), wvec(3), dwvec(3), r2, c0, c1, c2, exp_x2
      real(wp) :: xi_i, dxi_i, sqrt_pi, dc0_dxi, rw, ru, wu, rdw

      ngrid = size(w_grad, 2)
      nsph = size(za)
      ndir = size(dirs, 3)
      if (size(w_grad, 1) /= 3 .or. size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(dw_grad, 1) /= 3 .or. size(dw_grad, 2) /= ngrid .or. size(dw_grad, 3) /= ndir .or. &
          size(d_xyz, 1) /= 3 .or. size(d_xyz, 2) /= ngrid .or. size(d_xyz, 3) /= ndir .or. &
          size(dirs, 1) /= 3 .or. size(dirs, 2) /= nsph .or. &
          size(hvp_rA, 1) /= 3 .or. size(hvp_rA, 2) /= nsph .or. size(hvp_rA, 3) /= ndir) then
         call fatal_error(error, "pcm_elstat_normal_direct_hessian: array shape mismatch")
         return
      end if
      has_xi = present(xi)
      if (has_xi) then
         if (size(xi) /= ngrid) then
            call fatal_error(error, "pcm_elstat_normal_direct_hessian: width shape mismatch")
            return
         end if
         if (present(d_xi)) then
            if (size(d_xi, 1) /= ngrid .or. size(d_xi, 2) /= ndir) then
               call fatal_error(error, "pcm_elstat_normal_direct_hessian: width-tangent shape mismatch")
               return
            end if
         end if
      end if

      sqrt_pi = sqrt(acos(-1.0_wp))
      do katom = 1, nsph
         do i = 1, ngrid
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            if (.not. has_xi .and. r2 <= 1.0e-30_wp) cycle
            xi_i = 0.0_wp
            if (has_xi) xi_i = xi(i)
            call elstat_radial_kernels(r2, has_xi, xi_i, c0, c1, c2, exp_x2)
            dc0_dxi = 0.0_wp
            if (has_xi) dc0_dxi = 4.0_wp*xi_i*xi_i*exp_x2/sqrt_pi
            wvec = w_grad(:, i)
            rw = dot_product(rvec, wvec)
            do idir = 1, ndir
               uvec = d_xyz(:, i, idir) - dirs(:, katom, idir)
               dwvec = dw_grad(:, i, idir)
               ru = dot_product(rvec, uvec)
               wu = dot_product(wvec, uvec)
               rdw = dot_product(rvec, dwvec)
               dxi_i = 0.0_wp
               if (has_xi .and. present(d_xi)) dxi_i = d_xi(i, idir)
               hvp_rA(:, katom, idir) = hvp_rA(:, katom, idir) + za(katom)*( &
                  & c0*dwvec + c1*rdw*rvec + c1*(ru*wvec + rw*uvec + wu*rvec) &
                  & + c2*rw*ru*rvec + dxi_i*dc0_dxi*(wvec - 2.0_wp*xi_i*xi_i*rw*rvec))
            end do
         end do
      end do

   end subroutine pcm_elstat_normal_direct_hessian

   !> Point- or Gaussian-charge electrostatic potential and its surface derivatives
   subroutine pcm_elstat_potential_derivatives(xyz, sphxyz, za, phi, dphi_dr, error, &
                                               xi, dphi_dxi, d2phi_dr2, d2phi_dr_dxi)
      real(wp), intent(in) :: xyz(:, :), sphxyz(:, :), za(:)
      real(wp), intent(out) :: phi(:), dphi_dr(:, :)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: xi(:)
      real(wp), intent(out), optional :: dphi_dxi(:), d2phi_dr2(:, :, :), d2phi_dr_dxi(:, :)

      integer :: i, katom, a, b, ngrid, nsph
      logical :: has_xi
      real(wp) :: rvec(3), r2, r, x, x2, psi, c0, c1, c2, exp_x2, xi_i, sqrt_pi, dc0_dxi

      ngrid = size(phi)
      nsph = size(za)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(dphi_dr, 1) /= 3 .or. size(dphi_dr, 2) /= ngrid) then
         call fatal_error(error, "pcm_elstat_potential_derivatives: array shape mismatch")
         return
      end if
      has_xi = present(xi)
      if (has_xi) then
         if (size(xi) /= ngrid) then
            call fatal_error(error, "pcm_elstat_potential_derivatives: width shape mismatch")
            return
         end if
      end if

      phi = 0.0_wp
      dphi_dr = 0.0_wp
      if (present(dphi_dxi)) dphi_dxi = 0.0_wp
      if (present(d2phi_dr2)) d2phi_dr2 = 0.0_wp
      if (present(d2phi_dr_dxi)) d2phi_dr_dxi = 0.0_wp

      sqrt_pi = sqrt(acos(-1.0_wp))
      do i = 1, ngrid
         xi_i = 0.0_wp
         if (has_xi) xi_i = xi(i)
         do katom = 1, nsph
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            if (.not. has_xi .and. r2 <= 1.0e-30_wp) cycle
            call elstat_radial_kernels(r2, has_xi, xi_i, c0, c1, c2, exp_x2)
            r = sqrt(r2)
            if (has_xi) then
               x = xi_i*r
               x2 = x*x
               if (abs(x) < 0.05_wp) then
                  psi = (2.0_wp*xi_i/sqrt_pi)*(1.0_wp - x2/3.0_wp + x2*x2/10.0_wp - x2**3/42.0_wp)
               else
                  psi = erf(x)/r
               end if
            else
               psi = 1.0_wp/r
            end if
            phi(i) = phi(i) + za(katom)*psi
            dphi_dr(:, i) = dphi_dr(:, i) - za(katom)*c0*rvec
            if (has_xi .and. present(dphi_dxi)) then
               dphi_dxi(i) = dphi_dxi(i) + za(katom)*(2.0_wp/sqrt_pi)*exp_x2
            end if
            if (present(d2phi_dr2)) then
               do b = 1, 3
                  do a = 1, 3
                     d2phi_dr2(a, b, i) = d2phi_dr2(a, b, i) - za(katom)*c1*rvec(a)*rvec(b)
                  end do
                  d2phi_dr2(b, b, i) = d2phi_dr2(b, b, i) - za(katom)*c0
               end do
            end if
            if (has_xi .and. present(d2phi_dr_dxi)) then
               dc0_dxi = 4.0_wp*xi_i*xi_i*exp_x2/sqrt_pi
               d2phi_dr_dxi(:, i) = d2phi_dr_dxi(:, i) - za(katom)*dc0_dxi*rvec
            end if
         end do
      end do

   end subroutine pcm_elstat_potential_derivatives

   !> Total directional tangents of [[pcm_elstat_potential_derivatives]] along moving surface & nuclei
   subroutine pcm_elstat_potential_tangent(xyz, sphxyz, za, d_xyz, dirs, d_phi, d_dphi_dr, &
                                           error, xi, d_xi, d_dphi_dxi, d_d2phi_dr2, d_d2phi_dr_dxi)
      real(wp), intent(in) :: xyz(:, :), sphxyz(:, :), za(:)
      real(wp), intent(in) :: d_xyz(:, :, :), dirs(:, :, :)
      real(wp), intent(out) :: d_phi(:, :), d_dphi_dr(:, :, :)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: xi(:), d_xi(:, :)
      real(wp), intent(out), optional :: d_dphi_dxi(:, :), d_d2phi_dr2(:, :, :), d_d2phi_dr_dxi(:, :, :)

      integer :: i, katom, idir, k, a, b, ngrid, nsph, ndir
      logical :: has_xi
      real(wp) :: rvec(3), uvec(3), r2, c0, c1, c2, exp_x2, xi_i, dxi_i, sqrt_pi
      real(wp) :: ru, dc0_dxi, dpsi_dxi, entry, diag_delta

      ngrid = size(d_phi, 1)
      ndir = size(d_phi, 2)
      nsph = size(za)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid .or. &
          size(sphxyz, 1) /= 3 .or. size(sphxyz, 2) /= nsph .or. &
          size(d_xyz, 1) /= 3 .or. size(d_xyz, 2) /= ngrid .or. size(d_xyz, 3) /= ndir .or. &
          size(dirs, 1) /= 3 .or. size(dirs, 2) /= nsph .or. size(dirs, 3) /= ndir .or. &
          size(d_dphi_dr, 1) /= 3 .or. size(d_dphi_dr, 2) /= ngrid .or. size(d_dphi_dr, 3) /= ndir) then
         call fatal_error(error, "pcm_elstat_potential_tangent: array shape mismatch")
         return
      end if
      has_xi = present(xi)
      if (has_xi) then
         if (size(xi) /= ngrid) then
            call fatal_error(error, "pcm_elstat_potential_tangent: width shape mismatch")
            return
         end if
      end if

      d_phi = 0.0_wp
      d_dphi_dr = 0.0_wp
      if (present(d_dphi_dxi)) d_dphi_dxi = 0.0_wp
      if (present(d_d2phi_dr2)) d_d2phi_dr2 = 0.0_wp
      if (present(d_d2phi_dr_dxi)) d_d2phi_dr_dxi = 0.0_wp

      sqrt_pi = sqrt(acos(-1.0_wp))
      do i = 1, ngrid
         xi_i = 0.0_wp
         if (has_xi) xi_i = xi(i)
         do katom = 1, nsph
            rvec = xyz(:, i) - sphxyz(:, katom)
            r2 = sum(rvec*rvec)
            if (.not. has_xi .and. r2 <= 1.0e-30_wp) cycle
            call elstat_radial_kernels(r2, has_xi, xi_i, c0, c1, c2, exp_x2)
            dc0_dxi = 0.0_wp
            dpsi_dxi = 0.0_wp
            if (has_xi) then
               dc0_dxi = 4.0_wp*xi_i*xi_i*exp_x2/sqrt_pi
               dpsi_dxi = 2.0_wp*exp_x2/sqrt_pi
            end if
            do idir = 1, ndir
               uvec = d_xyz(:, i, idir) - dirs(:, katom, idir)
               ru = dot_product(rvec, uvec)
               dxi_i = 0.0_wp
               if (has_xi .and. present(d_xi)) dxi_i = d_xi(i, idir)

               d_phi(i, idir) = d_phi(i, idir) + za(katom)*(-c0*ru + dxi_i*dpsi_dxi)
               d_dphi_dr(:, i, idir) = d_dphi_dr(:, i, idir) - za(katom)*( &
                  & c0*uvec + c1*ru*rvec + dxi_i*dc0_dxi*rvec)

               if (has_xi .and. present(d_dphi_dxi)) then
                  d_dphi_dxi(i, idir) = d_dphi_dxi(i, idir) - za(katom)*( &
                     & dc0_dxi*ru + dxi_i*(4.0_wp*xi_i/sqrt_pi)*r2*exp_x2)
               end if
               if (present(d_d2phi_dr2)) then
                  do k = 1, 6
                     a = sym2_idx(1, k)
                     b = sym2_idx(2, k)
                     diag_delta = 0.0_wp
                     if (a == b) diag_delta = 1.0_wp
                     entry = c1*ru*diag_delta + c2*ru*rvec(a)*rvec(b) &
                        & + c1*(uvec(a)*rvec(b) + rvec(a)*uvec(b)) &
                        & + dxi_i*dc0_dxi*(diag_delta - 2.0_wp*xi_i*xi_i*rvec(a)*rvec(b))
                     d_d2phi_dr2(k, i, idir) = d_d2phi_dr2(k, i, idir) - za(katom)*entry
                  end do
               end if
               if (has_xi .and. present(d_d2phi_dr_dxi)) then
                  d_d2phi_dr_dxi(:, i, idir) = d_d2phi_dr_dxi(:, i, idir) - za(katom)*( &
                     & dc0_dxi*(uvec - 2.0_wp*xi_i*xi_i*ru*rvec) &
                     & + dxi_i*(8.0_wp*xi_i/sqrt_pi)*exp_x2*(1.0_wp - xi_i*xi_i*r2)*rvec)
               end if
            end do
         end do
      end do

   end subroutine pcm_elstat_potential_tangent

end module moist_model_continuum_component_pcm_elstat
