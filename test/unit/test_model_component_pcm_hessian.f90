!> Second-order unit tests for all nine PCM continuum components and electrostatic kernels
module test_model_component_pcm_hessian
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: check_moist_error, submit
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_surface_tangent, only: cavity_surface_tangent_type
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & coupling_make_view, coupling_close_view, gaussian_potential_request_type, &
      & point_potential_request_type, gaussian_potential_tangent_request_type, &
      & point_potential_tangent_request_type, nuclear_charge_request_type
   use moist_channels_response, only: response_type, response_item_type, response_clear, &
      & potential_adjoint_response_type, potential_adjoint_tangent_response_type
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, solver_type
   use moist_model_continuum_component_pcm_solvers, only: pcm_storage_dense, &
      & pcm_storage_matrix_free
   use moist_model_continuum_component_pcm_elstat, only: pcm_elstat_potential_derivatives, &
      & pcm_elstat_potential_tangent
   use test_model_component_pcm_fixtures, only: npcm, new_pcm, fd_cavity, seed_drop, &
      & get_seed, mol_triatomic, get_solute, host_charges, pcm_cpcm, pcm_iefpcm, pcm_dpcm

   implicit none(type, external)
   private

   public :: collect_model_component_pcm_hessian

   integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
   real(wp), parameter :: COEFF(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp
   integer, parameter :: sym2_idx(2, 6) = reshape([1, 1, 1, 2, 1, 3, 2, 2, 2, 3, 3, 3], [2, 6])

contains

   subroutine collect_model_component_pcm_hessian(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("elstat_potential_tangent_fd", test_elstat_potential_tangent_fd), &
         & new_unittest("all_nine_pcms_gaussian_fd", test_all_pcms_gaussian_fd), &
         & new_unittest("all_nine_pcms_point_potential_fd", test_all_pcms_point_potential_fd), &
         & new_unittest("all_nine_pcms_point_surface_fd", test_all_pcms_point_surface_fd), &
         & new_unittest("iterative_matrix_free_hvp_matches_dense", &
         &              test_iterative_mf_hvp_matches_dense) &
         & ]
   end subroutine collect_model_component_pcm_hessian

   !> Verify `pcm_elstat_potential_tangent` (Gaussian short/long range + point) against 4-pt FD
   subroutine test_elstat_potential_tangent_fd(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      integer, parameter :: ngrid = 4, nsph = 2, ndir = 2
      real(wp), parameter :: step = 2.0e-4_wp, tol = 1.0e-9_wp
      real(wp) :: xyz(3, ngrid), sphxyz(3, nsph), za(nsph), xi(ngrid)
      real(wp) :: d_xyz(3, ngrid, ndir), dirs(3, nsph, ndir), d_xi(ngrid, ndir)
      real(wp) :: d_phi(ngrid, ndir), d_dphi_dr(3, ngrid, ndir), d_dphi_dxi(ngrid, ndir)
      real(wp) :: d_d2phi_dr2(6, ngrid, ndir), d_d2phi_dr_dxi(3, ngrid, ndir)
      real(wp) :: phi_t(ngrid), dphi_dr_t(3, ngrid), dphi_dxi_t(ngrid)
      real(wp) :: d2phi_dr2_t(3, 3, ngrid), d2phi_dr_dxi_t(3, ngrid)
      real(wp) :: fd_phi(ngrid), fd_dr(3, ngrid), fd_dxi(ngrid)
      real(wp) :: fd_dr2(6, ngrid), fd_dr_dxi(3, ngrid)
      real(wp) :: xyz_t(3, ngrid), sph_t(3, nsph), xi_t(ngrid), t
      integer :: igrid, idir, ioff, iaxis, iatom, k, a, b, mode
      logical :: use_xi

      sphxyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.4_wp, -0.7_wp, 0.5_wp], [3, nsph])
      ! Point 1 is inside the |x| < 0.15 Taylor regime for atom 1 (Gaussian mode),
      ! while well-separated for point mode when shifted slightly
      xyz = reshape([ &
         & 0.03_wp, -0.02_wp, 0.04_wp, &
         & 0.65_wp, 0.45_wp, -0.35_wp, &
         & -1.2_wp, 0.9_wp, 0.7_wp, &
         & 2.3_wp, -1.5_wp, 1.1_wp], [3, ngrid])
      za = [2.2_wp, -1.3_wp]
      xi = [0.85_wp, 1.15_wp, 0.65_wp, 1.45_wp]

      do idir = 1, ndir
         do igrid = 1, ngrid
            d_xi(igrid, idir) = 0.25_wp*cos(0.7_wp*igrid + 1.1_wp*idir)
            do iaxis = 1, 3
               d_xyz(iaxis, igrid, idir) = 0.3_wp*sin(0.5_wp*igrid + 0.9_wp*iaxis - idir)
            end do
         end do
         do iatom = 1, nsph
            do iaxis = 1, 3
               dirs(iaxis, iatom, idir) = 0.35_wp*cos(0.8_wp*iatom + 0.6_wp*iaxis + idir)
            end do
         end do
      end do

      do mode = 1, 2
         use_xi = (mode == 1)
         if (.not. use_xi) xyz(:, 1) = [0.75_wp, -0.65_wp, 0.85_wp]
         if (use_xi) then
            call pcm_elstat_potential_tangent(xyz, sphxyz, za, d_xyz, dirs, d_phi, d_dphi_dr, &
               & err, xi=xi, d_xi=d_xi, d_dphi_dxi=d_dphi_dxi, d_d2phi_dr2=d_d2phi_dr2, &
               & d_d2phi_dr_dxi=d_d2phi_dr_dxi)
         else
            call pcm_elstat_potential_tangent(xyz, sphxyz, za, d_xyz, dirs, d_phi, d_dphi_dr, &
               & err, d_d2phi_dr2=d_d2phi_dr2)
         end if
         call check_moist_error(error, err, "pcm_elstat_potential_tangent")
         if (allocated(error)) return

         do idir = 1, ndir
            fd_phi = 0.0_wp
            fd_dr = 0.0_wp
            fd_dxi = 0.0_wp
            fd_dr2 = 0.0_wp
            fd_dr_dxi = 0.0_wp
            do ioff = 1, 4
               t = real(OFFSET(ioff), wp)*step
               xyz_t = xyz + t*d_xyz(:, :, idir)
               sph_t = sphxyz + t*dirs(:, :, idir)
               xi_t = xi + t*d_xi(:, idir)
               if (use_xi) then
                  call pcm_elstat_potential_derivatives(xyz_t, sph_t, za, phi_t, dphi_dr_t, err, &
                     & xi=xi_t, dphi_dxi=dphi_dxi_t, d2phi_dr2=d2phi_dr2_t, &
                     & d2phi_dr_dxi=d2phi_dr_dxi_t)
               else
                  call pcm_elstat_potential_derivatives(xyz_t, sph_t, za, phi_t, dphi_dr_t, err, &
                     & d2phi_dr2=d2phi_dr2_t)
               end if
               call check_moist_error(error, err, "pcm_elstat_potential_derivatives")
               if (allocated(error)) return

               fd_phi = fd_phi + COEFF(ioff)*phi_t/step
               fd_dr = fd_dr + COEFF(ioff)*dphi_dr_t/step
               if (use_xi) then
                  fd_dxi = fd_dxi + COEFF(ioff)*dphi_dxi_t/step
                  fd_dr_dxi = fd_dr_dxi + COEFF(ioff)*d2phi_dr_dxi_t/step
               end if
               do k = 1, 6
                  a = sym2_idx(1, k)
                  b = sym2_idx(2, k)
                  fd_dr2(k, :) = fd_dr2(k, :) + COEFF(ioff)*d2phi_dr2_t(a, b, :)/step
               end do
            end do

            call check(error, maxval(abs(d_phi(:, idir) - fd_phi)), 0.0_wp, thr=tol, &
               & more="d_phi FD mismatch")
            if (allocated(error)) return
            call check(error, maxval(abs(d_dphi_dr(:, :, idir) - fd_dr)), 0.0_wp, thr=tol, &
               & more="d_dphi_dr FD mismatch")
            if (allocated(error)) return
            call check(error, maxval(abs(d_d2phi_dr2(:, :, idir) - fd_dr2)), 0.0_wp, thr=tol, &
               & more="d_d2phi_dr2 FD mismatch")
            if (allocated(error)) return
            if (use_xi) then
               call check(error, maxval(abs(d_dphi_dxi(:, idir) - fd_dxi)), 0.0_wp, thr=tol, &
                  & more="d_dphi_dxi FD mismatch")
               if (allocated(error)) return
               call check(error, maxval(abs(d_d2phi_dr_dxi(:, :, idir) - fd_dr_dxi)), 0.0_wp, &
                  & thr=tol, more="d_d2phi_dr_dxi FD mismatch")
               if (allocated(error)) return
            end if
         end do
      end do
   end subroutine test_elstat_potential_tangent_fd

   subroutine test_all_pcms_gaussian_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: kind

      do kind = 1, npcm
         call verify_pcm_component_hessian(kind, smeared=.true., point_pot=.false., error=error)
         if (allocated(error)) return
      end do
   end subroutine test_all_pcms_gaussian_fd

   subroutine test_all_pcms_point_potential_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: kind

      do kind = 1, npcm
         call verify_pcm_component_hessian(kind, smeared=.true., point_pot=.true., error=error)
         if (allocated(error)) return
      end do
   end subroutine test_all_pcms_point_potential_fd

   subroutine test_all_pcms_point_surface_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: kind

      do kind = 1, npcm
         call verify_pcm_component_hessian(kind, smeared=.false., point_pot=.false., error=error)
         if (allocated(error)) return
      end do
   end subroutine test_all_pcms_point_surface_fd

   !> Verify iterative GMRES (dense and matrix-free storage) matches dense HVP hooks
   subroutine test_iterative_mf_hvp_matches_dense(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: kinds(3) = [pcm_cpcm, pcm_iefpcm, pcm_dpcm]
      integer, parameter :: ndir = 2
      real(wp), parameter :: tol = 1.0e-9_wp

      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: err
      type(structure_type) :: mol
      type(fd_cavity) :: cavity
      class(model_continuum_component_pcm), allocatable :: pcm_ref, pcm_it
      type(coupling_type), target :: cpl_ref, cpl_it
      type(coupling_view_type) :: view
      type(cavity_surface_tangent_type) :: tangent
      type(cavity_surface_adjoint_type), allocatable :: dacc_ref(:), dacc_it(:)
      type(response_type) :: resp_ref, resp_it
      class(response_item_type), allocatable :: it_ref, it_chk
      real(wp), allocatable :: qsrc(:), z(:), dirs(:, :, :), hvp_ref(:, :, :), hvp_it(:, :, :)
      real(wp) :: dn_raw(3)
      integer :: ikind, imode, ngrid, nat, idir, igrid, iaxis
      character(len=:), allocatable :: tag

      call new_context(ctx, verbosity=0)
      call get_solute(mol_triatomic, mol)
      qsrc = host_charges(mol)
      z = real(mol%num(mol%id), wp)
      call get_seed(seed_drop, .true., mol, cavity, error)
      if (allocated(error)) return
      ngrid = cavity%ngrid
      nat = mol%nat

      allocate (dirs(3, nat, ndir))
      do idir = 1, ndir
         do igrid = 1, nat
            do iaxis = 1, 3
               dirs(iaxis, igrid, idir) = 0.25_wp*cos(0.5_wp*iaxis + 0.9_wp*igrid - 1.1_wp*idir)
            end do
         end do
      end do
      call tangent%init(ngrid, ndir, want_curvature=.true.)
      do idir = 1, ndir
         do igrid = 1, ngrid
            tangent%d_a(igrid, idir) = 0.08_wp*cavity%a(igrid)*sin(0.31_wp*igrid + idir)
            tangent%d_w(igrid, idir) = 0.05_wp*cos(0.27_wp*igrid - idir)
            tangent%d_k1(igrid, idir) = 0.07_wp*sin(0.41_wp*igrid + 0.7_wp*idir)
            tangent%d_k2(igrid, idir) = 0.06_wp*cos(0.53_wp*igrid - 0.5_wp*idir)
            tangent%d_xi(igrid, idir) = 0.07_wp*cavity%xi0(igrid)*cos(0.37_wp*igrid + idir)
            tangent%d_f(igrid, idir) = 0.04_wp*cavity%f(igrid)*sin(0.29_wp*igrid - idir)
            do iaxis = 1, 3
               tangent%d_xyz(iaxis, igrid, idir) = 0.15_wp*sin(0.23_wp*igrid + 0.7_wp*iaxis + idir)
               dn_raw(iaxis) = 0.18_wp*cos(0.19_wp*igrid + 0.8_wp*iaxis - idir)
            end do
            tangent%d_n(:, igrid, idir) = dn_raw &
               & - dot_product(cavity%normal0(:, igrid), dn_raw)*cavity%normal0(:, igrid)
         end do
      end do

      allocate (dacc_ref(ndir), dacc_it(ndir), hvp_ref(3, nat, ndir), hvp_it(3, nat, ndir))
      do ikind = 1, size(kinds)
         call new_pcm(kinds(ikind), ctx, pcm_ref, error)
         if (allocated(error)) return
         call pcm_ref%update(mol, cavity, err)
         if (.not. allocated(err)) call pcm_ref%new_coupling(cavity, cpl_ref, err)
         if (.not. allocated(err)) call pcm_ref%prepare_hvp(cavity, cpl_ref, dirs, tangent, err)
         call check_moist_error(error, err, "dense reference setup")
         if (allocated(error)) return
         call answer_hvp_coupling(cavity, cpl_ref, mol, qsrc, z, error)
         if (allocated(error)) return

         do idir = 1, ndir
            call dacc_ref(idir)%init(ngrid)
         end do
         hvp_ref = 0.0_wp
         call response_clear(resp_ref)
         call coupling_make_view(cpl_ref, 1, view)
         call pcm_ref%get_hessian_surface_weights(view, cavity, dirs, tangent, dacc_ref, err)
         if (.not. allocated(err)) call pcm_ref%get_direct_hessian(view, cavity, dirs, tangent, &
            & hvp_ref, err)
         if (.not. allocated(err)) call pcm_ref%get_trace_response_tangent(view, cavity, tangent, &
            & resp_ref, err)
         call coupling_close_view(cpl_ref)
         call check_moist_error(error, err, "dense reference eval")
         if (allocated(error)) return
         call check(error, resp_ref%next(), more="dense missing potential_adjoint_tangent")
         if (allocated(error)) return
         it_ref = resp_ref%item()

         do imode = 1, 2
            call new_pcm(kinds(ikind), ctx, pcm_it, error)
            if (allocated(error)) return
            pcm_it%param%solver = solver_type%iterative
            pcm_it%param%solver_tol = 1.0e-12_wp
            pcm_it%param%solver_maxiter = 500
            if (imode == 1) then
               pcm_it%param%storage = pcm_storage_dense
               tag = trim(pcm_it%name)//" iterative_dense"
            else
               pcm_it%param%storage = pcm_storage_matrix_free
               tag = trim(pcm_it%name)//" iterative_mf"
            end if
            call pcm_it%update(mol, cavity, err)
            if (.not. allocated(err)) call pcm_it%new_coupling(cavity, cpl_it, err)
            if (.not. allocated(err)) call pcm_it%prepare_hvp(cavity, cpl_it, dirs, tangent, err)
            call check_moist_error(error, err, tag//" setup")
            if (allocated(error)) return
            call answer_hvp_coupling(cavity, cpl_it, mol, qsrc, z, error)
            if (allocated(error)) return

            do idir = 1, ndir
               call dacc_it(idir)%init(ngrid)
            end do
            hvp_it = 0.0_wp
            call response_clear(resp_it)
            call coupling_make_view(cpl_it, 1, view)
            call pcm_it%get_hessian_surface_weights(view, cavity, dirs, tangent, dacc_it, err)
            if (.not. allocated(err)) call pcm_it%get_direct_hessian(view, cavity, dirs, tangent, &
               & hvp_it, err)
            if (.not. allocated(err)) call pcm_it%get_trace_response_tangent(view, cavity, tangent, &
               & resp_it, err)
            call coupling_close_view(cpl_it)
            call check_moist_error(error, err, tag//" eval")
            if (allocated(error)) return

            do idir = 1, ndir
               call check_close(error, reshape(dacc_it(idir)%w_xyz, [3*ngrid]), &
                  & reshape(dacc_ref(idir)%w_xyz, [3*ngrid]), tol, tol, tag//" w_xyz")
               if (allocated(error)) return
               call check_close(error, dacc_it(idir)%w_xi, dacc_ref(idir)%w_xi, tol, tol, tag//" w_xi")
               if (allocated(error)) return
               call check_close(error, dacc_it(idir)%w_f, dacc_ref(idir)%w_f, tol, tol, tag//" w_f")
               if (allocated(error)) return
               call check_close(error, dacc_it(idir)%w_a, dacc_ref(idir)%w_a, tol, tol, tag//" w_a")
               if (allocated(error)) return
               if (pcm_ref%uses_normal_trace()) then
                  call check_close(error, reshape(dacc_it(idir)%w_n, [3*ngrid]), &
                     & reshape(dacc_ref(idir)%w_n, [3*ngrid]), tol, tol, tag//" w_n")
                  if (allocated(error)) return
                  call check_close(error, dacc_it(idir)%w_k1, dacc_ref(idir)%w_k1, tol, tol, tag//" w_k1")
                  if (allocated(error)) return
                  call check_close(error, dacc_it(idir)%w_k2, dacc_ref(idir)%w_k2, tol, tol, tag//" w_k2")
                  if (allocated(error)) return
               end if
            end do
            call check_close(error, reshape(hvp_it, [3*nat*ndir]), &
               & reshape(hvp_ref, [3*nat*ndir]), tol, tol, tag//" hvp")
            if (allocated(error)) return

            call check(error, resp_it%next(), more=tag//": missing potential_adjoint_tangent")
            if (allocated(error)) return
            it_chk = resp_it%item()
            select type (it_ref)
            type is (potential_adjoint_tangent_response_type)
               select type (it_chk)
               type is (potential_adjoint_tangent_response_type)
                  call check_close(error, reshape(it_chk%dw_phi, [ngrid*ndir]), &
                     & reshape(it_ref%dw_phi, [ngrid*ndir]), tol, tol, tag//" dw_phi")
                  if (allocated(error)) return
                  if (pcm_ref%uses_normal_trace()) then
                     call check_close(error, reshape(it_chk%dw_dphi_dr, [3*ngrid*ndir]), &
                        & reshape(it_ref%dw_dphi_dr, [3*ngrid*ndir]), tol, tol, tag//" dw_dphi_dr")
                     if (allocated(error)) return
                  end if
               class default
                  call test_failed(error, tag//": unexpected response item")
                  return
               end select
            end select
         end do
      end do
   end subroutine test_iterative_mf_hvp_matches_dense

   !> Check `get_hessian_surface_weights`, `get_direct_hessian` and `get_trace_response_tangent`
   !> against 4-point central differences for one PCM family and surface/potential configuration
   subroutine verify_pcm_component_hessian(kind, smeared, point_pot, error)
      integer, intent(in) :: kind
      logical, intent(in) :: smeared, point_pot
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: step = 1.5e-4_wp, atol = 5.0e-8_wp, rtol = 5.0e-7_wp
      integer, parameter :: ndir = 2

      type(moist_context_type), target :: ctx
      type(moist_error_type), allocatable :: err
      type(structure_type) :: mol, mol_t
      type(fd_cavity) :: cavity, cav_t
      class(model_continuum_component_pcm), allocatable :: pcm, pcm_t
      type(coupling_type), target :: coupling, cpl_t
      type(coupling_view_type) :: view
      type(cavity_surface_tangent_type) :: tangent
      type(cavity_surface_adjoint_type), allocatable :: dacc(:)
      type(cavity_surface_adjoint_type) :: w_t, fd_w
      type(response_type) :: resp_tan, resp_t
      class(response_item_type), allocatable :: item
      real(wp), allocatable :: qsrc(:), z(:), dirs(:, :, :), hvp(:, :, :)
      real(wp), allocatable :: grad_t(:, :), fd_hvp(:, :)
      real(wp), allocatable :: fd_dw_phi(:), fd_dw_dphi_dr(:, :)
      real(wp), allocatable :: an_dw_phi(:, :), an_dw_dphi_dr(:, :, :)
      real(wp) :: dn_raw(3), n_pert(3), t
      integer :: ngrid, nat, igrid, idir, ioff, iaxis, iatom
      character(len=:), allocatable :: tag

      call new_context(ctx, verbosity=0)
      call get_solute(mol_triatomic, mol)
      qsrc = host_charges(mol)
      z = real(mol%num(mol%id), wp)
      call get_seed(seed_drop, smeared, mol, cavity, error)
      if (allocated(error)) return
      ngrid = cavity%ngrid
      nat = mol%nat

      call new_pcm(kind, ctx, pcm, error)
      if (allocated(error)) return
      pcm%param%point_potential = point_pot
      call pcm%update(mol, cavity, err)
      call check_moist_error(error, err, "PCM update")
      if (allocated(error)) return
      tag = trim(pcm%name)

      ! Build non-trivial nuclear directions and surface tangents (with d_n tangent to unit sphere)
      allocate (dirs(3, nat, ndir))
      do idir = 1, ndir
         do iatom = 1, nat
            do iaxis = 1, 3
               dirs(iaxis, iatom, idir) = 0.25_wp*cos(0.5_wp*iaxis + 0.9_wp*iatom - 1.1_wp*idir)
            end do
         end do
      end do

      call tangent%init(ngrid, ndir, want_curvature=.true.)
      do idir = 1, ndir
         do igrid = 1, ngrid
            tangent%d_a(igrid, idir) = 0.08_wp*cavity%a(igrid)*sin(0.31_wp*igrid + idir)
            tangent%d_w(igrid, idir) = 0.05_wp*cos(0.27_wp*igrid - idir)
            tangent%d_k1(igrid, idir) = 0.07_wp*sin(0.41_wp*igrid + 0.7_wp*idir)
            tangent%d_k2(igrid, idir) = 0.06_wp*cos(0.53_wp*igrid - 0.5_wp*idir)
            if (smeared) then
               tangent%d_xi(igrid, idir) = 0.07_wp*cavity%xi0(igrid)*cos(0.37_wp*igrid + idir)
               tangent%d_f(igrid, idir) = 0.04_wp*cavity%f(igrid)*sin(0.29_wp*igrid - idir)
            end if
            do iaxis = 1, 3
               tangent%d_xyz(iaxis, igrid, idir) = 0.15_wp*sin(0.23_wp*igrid + 0.7_wp*iaxis + idir)
               dn_raw(iaxis) = 0.18_wp*cos(0.19_wp*igrid + 0.8_wp*iaxis - idir)
            end do
            tangent%d_n(:, igrid, idir) = dn_raw &
               & - dot_product(cavity%normal0(:, igrid), dn_raw)*cavity%normal0(:, igrid)
         end do
      end do

      ! Stage and answer HVP on the bare component
      call pcm%new_coupling(cavity, coupling, err)
      if (.not. allocated(err)) call pcm%prepare_hvp(cavity, coupling, dirs, tangent, err)
      call check_moist_error(error, err, tag//" prepare_hvp")
      if (allocated(error)) return
      call answer_hvp_coupling(cavity, coupling, mol, qsrc, z, error)
      if (allocated(error)) return

      allocate (dacc(ndir))
      allocate (hvp(3, nat, ndir), source=0.0_wp)
      do idir = 1, ndir
         call dacc(idir)%init(ngrid)
      end do

      call coupling_make_view(coupling, 1, view)
      call pcm%get_hessian_surface_weights(view, cavity, dirs, tangent, dacc, err)
      if (.not. allocated(err)) then
         call pcm%get_direct_hessian(view, cavity, dirs, tangent, hvp, err)
      end if
      if (.not. allocated(err)) then
         call pcm%get_trace_response_tangent(view, cavity, tangent, resp_tan, err)
      end if
      call coupling_close_view(coupling)
      call check_moist_error(error, err, tag//" analytic HVP hooks")
      if (allocated(error)) return

      call check(error, resp_tan%next(), more=tag//": missing potential_adjoint_tangent")
      if (allocated(error)) return
      item = resp_tan%item()
      select type (item)
      type is (potential_adjoint_tangent_response_type)
         an_dw_phi = item%dw_phi
         if (pcm%uses_normal_trace()) then
            call check(error, allocated(item%dw_dphi_dr), more=tag//": missing dw_dphi_dr")
            if (allocated(error)) return
            an_dw_dphi_dr = item%dw_dphi_dr
         end if
      class default
         call test_failed(error, tag//": unexpected response item")
         return
      end select

      ! Finite-difference check over all perturbed observables + nuclei
      allocate (grad_t(3, nat), fd_hvp(3, nat), fd_dw_phi(ngrid), fd_dw_dphi_dr(3, ngrid))
      do idir = 1, ndir
         call fd_w%init(ngrid)
         fd_hvp = 0.0_wp
         fd_dw_phi = 0.0_wp
         fd_dw_dphi_dr = 0.0_wp

         do ioff = 1, 4
            t = real(OFFSET(ioff), wp)*step
            mol_t = mol
            mol_t%xyz = mol%xyz + t*dirs(:, :, idir)
            cav_t = cavity
            cav_t%xyz = cavity%xyz + t*tangent%d_xyz(:, :, idir)
            cav_t%a = cavity%a + t*tangent%d_a(:, idir)
            cav_t%k1 = cavity%k1 + t*tangent%d_k1(:, idir)
            cav_t%k2 = cavity%k2 + t*tangent%d_k2(:, idir)
            if (smeared) then
               cav_t%xi0 = cavity%xi0 + t*tangent%d_xi(:, idir)
               cav_t%f = cavity%f + t*tangent%d_f(:, idir)
            end if
            do igrid = 1, ngrid
               n_pert = cavity%normal0(:, igrid) + t*tangent%d_n(:, igrid, idir)
               cav_t%normal0(:, igrid) = n_pert/norm2(n_pert)
            end do

            call new_pcm(kind, ctx, pcm_t, error)
            if (allocated(error)) return
            pcm_t%param%point_potential = point_pot
            call pcm_t%update(mol_t, cav_t, err)
            if (.not. allocated(err)) call pcm_t%new_coupling(cav_t, cpl_t, err)
            if (.not. allocated(err)) call pcm_t%prepare_gradient(cav_t, cpl_t, err)
            call check_moist_error(error, err, tag//" displaced gradient setup")
            if (allocated(error)) return
            call answer_hvp_coupling(cav_t, cpl_t, mol_t, qsrc, z, error)
            if (allocated(error)) return

            call w_t%init(ngrid)
            call response_clear(resp_t)
            grad_t = 0.0_wp
            call coupling_make_view(cpl_t, 1, view)
            call pcm_t%get_gradient_surface_weights(view, cav_t, w_t, err)
            if (.not. allocated(err)) call pcm_t%get_direct_gradient(view, cav_t, grad_t, err)
            if (.not. allocated(err)) call pcm_t%get_trace_response(view, cav_t, resp_t, err)
            call coupling_close_view(cpl_t)
            call check_moist_error(error, err, tag//" displaced gradient eval")
            if (allocated(error)) return

            fd_w%w_xyz = fd_w%w_xyz + COEFF(ioff)*w_t%w_xyz/step
            fd_w%w_n = fd_w%w_n + COEFF(ioff)*w_t%w_n/step
            fd_w%w_a = fd_w%w_a + COEFF(ioff)*w_t%w_a/step
            fd_w%w_xi = fd_w%w_xi + COEFF(ioff)*w_t%w_xi/step
            fd_w%w_f = fd_w%w_f + COEFF(ioff)*w_t%w_f/step
            fd_w%w_k1 = fd_w%w_k1 + COEFF(ioff)*w_t%w_k1/step
            fd_w%w_k2 = fd_w%w_k2 + COEFF(ioff)*w_t%w_k2/step
            fd_hvp = fd_hvp + COEFF(ioff)*grad_t/step

            call check(error, resp_t%next(), more=tag//": missing displaced potential_adjoint")
            if (allocated(error)) return
            item = resp_t%item()
            select type (item)
            type is (potential_adjoint_response_type)
               fd_dw_phi = fd_dw_phi + COEFF(ioff)*item%w_phi/step
               if (pcm%uses_normal_trace()) then
                  fd_dw_dphi_dr = fd_dw_dphi_dr + COEFF(ioff)*item%w_dphi_dr/step
               end if
            end select
         end do

         call check_close(error, reshape(dacc(idir)%w_xyz, [3*ngrid]), reshape(fd_w%w_xyz, [3*ngrid]), &
            & atol, rtol, tag//" w_xyz")
         if (allocated(error)) return
         call check_close(error, reshape(dacc(idir)%w_n, [3*ngrid]), reshape(fd_w%w_n, [3*ngrid]), &
            & atol, rtol, tag//" w_n")
         if (allocated(error)) return
         call check_close(error, dacc(idir)%w_a, fd_w%w_a, atol, rtol, tag//" w_a")
         if (allocated(error)) return
         call check_close(error, dacc(idir)%w_xi, fd_w%w_xi, atol, rtol, tag//" w_xi")
         if (allocated(error)) return
         call check_close(error, dacc(idir)%w_f, fd_w%w_f, atol, rtol, tag//" w_f")
         if (allocated(error)) return
         call check_close(error, dacc(idir)%w_k1, fd_w%w_k1, atol, rtol, tag//" w_k1")
         if (allocated(error)) return
         call check_close(error, dacc(idir)%w_k2, fd_w%w_k2, atol, rtol, tag//" w_k2")
         if (allocated(error)) return
         call check_close(error, reshape(hvp(:, :, idir), [3*nat]), reshape(fd_hvp, [3*nat]), &
            & atol, rtol, tag//" direct_hessian")
         if (allocated(error)) return
         call check_close(error, an_dw_phi(:, idir), fd_dw_phi, atol, rtol, tag//" dw_phi")
         if (allocated(error)) return
         if (pcm%uses_normal_trace()) then
            call check_close(error, reshape(an_dw_dphi_dr(:, :, idir), [3*ngrid]), &
               & reshape(fd_dw_dphi_dr, [3*ngrid]), atol, rtol, tag//" dw_dphi_dr")
            if (allocated(error)) return
         end if
      end do
   end subroutine verify_pcm_component_hessian

   !> Answer primal and tangent potential/nuclear requests from point/Gaussian charges `qsrc`
   subroutine answer_hvp_coupling(cavity, coupling, src, qsrc, z, error)
      type(fd_cavity), intent(in) :: cavity
      type(coupling_type), intent(inout), target :: coupling
      type(structure_type), intent(in) :: src
      real(wp), intent(in) :: qsrc(:), z(:)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: phi(:), dphi_dr(:, :), dphi_dxi(:), d2phi_dr2(:, :, :), d2phi_dr_dxi(:, :)
      real(wp), allocatable :: d_phi(:, :), d_dphi_dr(:, :, :), d_dphi_dxi(:, :)
      real(wp), allocatable :: d_d2phi_dr2(:, :, :), d_d2phi_dr_dxi(:, :, :)
      integer :: ngrid, ndir

      ngrid = cavity%ngrid
      do while (coupling%next())
         select type (req => coupling%request())
         type is (gaussian_potential_request_type)
            allocate (phi(ngrid), dphi_dr(3, ngrid), dphi_dxi(ngrid), &
               & d2phi_dr2(3, 3, ngrid), d2phi_dr_dxi(3, ngrid))
            call pcm_elstat_potential_derivatives(cavity%xyz, src%xyz, qsrc, phi, dphi_dr, err, &
               & xi=cavity%xi0, dphi_dxi=dphi_dxi, d2phi_dr2=d2phi_dr2, d2phi_dr_dxi=d2phi_dr_dxi)
            call check_moist_error(error, err, "gaussian_potential derivatives")
            if (allocated(error)) return
            call submit(coupling, "phi", phi)
            call submit(coupling, "dphi_dr", dphi_dr)
            call submit(coupling, "dphi_dxi", dphi_dxi)
            call submit(coupling, "d2phi_dr2", d2phi_dr2)
            call submit(coupling, "d2phi_dr_dxi", d2phi_dr_dxi)
            deallocate (phi, dphi_dr, dphi_dxi, d2phi_dr2, d2phi_dr_dxi)
         type is (point_potential_request_type)
            allocate (phi(ngrid), dphi_dr(3, ngrid), d2phi_dr2(3, 3, ngrid))
            call pcm_elstat_potential_derivatives(cavity%xyz, src%xyz, qsrc, phi, dphi_dr, err, &
               & d2phi_dr2=d2phi_dr2)
            call check_moist_error(error, err, "point_potential derivatives")
            if (allocated(error)) return
            call submit(coupling, "phi", phi)
            call submit(coupling, "dphi_dr", dphi_dr)
            call submit(coupling, "d2phi_dr2", d2phi_dr2)
            deallocate (phi, dphi_dr, d2phi_dr2)
         type is (gaussian_potential_tangent_request_type)
            ndir = size(req%dirs, 3)
            allocate (d_phi(ngrid, ndir), d_dphi_dr(3, ngrid, ndir), d_dphi_dxi(ngrid, ndir), &
               & d_d2phi_dr2(6, ngrid, ndir), d_d2phi_dr_dxi(3, ngrid, ndir))
            call pcm_elstat_potential_tangent(req%xyz, src%xyz, qsrc, req%d_xyz, req%dirs, &
               & d_phi, d_dphi_dr, err, xi=req%xi, d_xi=req%d_xi, d_dphi_dxi=d_dphi_dxi, &
               & d_d2phi_dr2=d_d2phi_dr2, d_d2phi_dr_dxi=d_d2phi_dr_dxi)
            call check_moist_error(error, err, "gaussian_potential_tangent")
            if (allocated(error)) return
            call submit(coupling, "d_phi", d_phi)
            call submit(coupling, "d_dphi_dr", d_dphi_dr)
            call submit(coupling, "d_dphi_dxi", d_dphi_dxi)
            if (req%is_missing("d_d2phi_dr2")) call submit(coupling, "d_d2phi_dr2", d_d2phi_dr2)
            if (req%is_missing("d_d2phi_dr_dxi")) call submit(coupling, "d_d2phi_dr_dxi", d_d2phi_dr_dxi)
            deallocate (d_phi, d_dphi_dr, d_dphi_dxi, d_d2phi_dr2, d_d2phi_dr_dxi)
         type is (point_potential_tangent_request_type)
            ndir = size(req%dirs, 3)
            allocate (d_phi(ngrid, ndir), d_dphi_dr(3, ngrid, ndir), d_d2phi_dr2(6, ngrid, ndir))
            call pcm_elstat_potential_tangent(req%xyz, src%xyz, qsrc, req%d_xyz, req%dirs, &
               & d_phi, d_dphi_dr, err, d_d2phi_dr2=d_d2phi_dr2)
            call check_moist_error(error, err, "point_potential_tangent")
            if (allocated(error)) return
            call submit(coupling, "d_phi", d_phi)
            call submit(coupling, "d_dphi_dr", d_dphi_dr)
            if (req%is_missing("d_d2phi_dr2")) call submit(coupling, "d_d2phi_dr2", d_d2phi_dr2)
            deallocate (d_phi, d_dphi_dr, d_d2phi_dr2)
         type is (nuclear_charge_request_type)
            call submit(coupling, "z", z)
         end select
      end do
   end subroutine answer_hvp_coupling

   !> Combined absolute + relative elementwise check
   subroutine check_close(error, val, ref, atol, rtol, label)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: val(:), ref(:), atol, rtol
      character(len=*), intent(in) :: label
      integer :: i
      real(wp) :: diff, bound

      do i = 1, size(val)
         diff = abs(val(i) - ref(i))
         bound = max(atol, rtol*abs(ref(i)))
         if (diff > bound) then
            call check(error, val(i), ref(i), thr=bound, more=label)
            return
         end if
      end do
   end subroutine check_close

end module test_model_component_pcm_hessian
