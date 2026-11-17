!> End-to-end nuclear Hessian and two-stage HVP tests of `model_continuum_type`
module test_model_continuum_hessian
   use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_associated, &
      & c_f_pointer, c_funloc, c_loc
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: check_moist_error, submit, drop_fixture_geometry, &
      & build_drop_test_cavity, LSF_SVDW, LSF_CFC, FIX_PLAIN, FIX_CROSS
   use test_model_component_pcm_fixtures, only: pcm_cpcm, pcm_iefpcm, pcm_dpcm, &
      & pcm_slic, pcm_nonlocal, new_pcm
   use moist_context, only: moist_context_type, new_context
   use moist_radii, only: default_cpcm_radii
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, &
      & moist_cavity_iswig_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop, drop_hvp_per_dir_max
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_isodensity_param, only: moist_cavity_drop_lsf_isodensity_param_type
   use moist_cavity_drop_lsf_isodensity_internal, only: &
      & moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_drop_lsf_isodensity_callback, only: &
      & moist_cavity_drop_lsf_isodensity_callback_type
   use moist_channels_coupling, only: coupling_type, density_matrix_request_type, &
      & nuclear_charge_request_type, gaussian_potential_request_type, &
      & point_potential_request_type, gaussian_potential_tangent_request_type, &
      & point_potential_tangent_request_type, density_tangent_request_type
   use moist_channels_response, only: response_type, response_item_type, &
      & potential_adjoint_response_type, density_response_type, &
      & potential_adjoint_tangent_response_type, density_tangent_response_type
   use moist_model_continuum_type, only: model_continuum_type, new_continuum_model
   use moist_model_continuum_component_pv, only: model_continuum_component_pv, new_component_pv
   use moist_model_continuum_component_gostshyp, only: model_continuum_component_gostshyp, &
      & new_component_gostshyp
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm
   use moist_model_continuum_component_pcm_elstat, only: pcm_elstat_potential_derivatives, &
      & pcm_elstat_potential_tangent

   implicit none(type, external)
   private

   public :: collect_model_continuum_hessian

   integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
   real(wp), parameter :: COEFF(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp
   integer, parameter :: sym2_idx(2, 6) = reshape([1, 1, 1, 2, 1, 3, 2, 2, 2, 3, 3, 3], [2, 6])
   integer, parameter :: sym3_idx(3, 10) = reshape([ &
      & 1, 1, 1, 1, 1, 2, 1, 1, 3, 1, 2, 2, 1, 2, 3, &
      & 1, 3, 3, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3], [3, 10])

   !> Callback context for atom-centered s-Gaussian density rho(r) = sum_A d_A exp(-2*alpha_A*|r - R_A|^2)
   type :: s_density_ctx_type
      real(wp), allocatable :: xyz(:, :)
      real(wp), allocatable :: alpha(:)
      real(wp), allocatable :: diag_d(:)
   end type s_density_ctx_type

contains

   subroutine collect_model_continuum_hessian(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("pv_hessian_fd", test_pv_hessian_fd), &
         & new_unittest("pv_cross_hessian_fd", test_pv_cross_hessian_fd), &
         & new_unittest("pv_hvp_matches_dense", test_pv_hvp_matches_dense), &
         & new_unittest("pcm_svdw_and_cfc_hvp_fd", test_pcm_geometric_hvp_fd), &
         & new_unittest("pcm_cross_branch_hvp_fd", test_pcm_cross_hvp_fd), &
         & new_unittest("multi_component_hvp_fd", test_multi_component_hvp_fd), &
         & new_unittest("isodensity_internal_and_callback_hvp_fd", test_isodensity_hvp_fd), &
         & new_unittest("pure_density_and_charge_hvp_fd", test_pure_density_and_charge_hvp_fd), &
         & new_unittest("pcm_multiblock_chunk_hvp", test_pcm_multiblock_chunk_hvp), &
         & new_unittest("hessian_guards", test_hessian_guards) &
         & ]
   end subroutine collect_model_continuum_hessian

   !* ================================================================================= *!
   !*                                   PV Hessians                                     *!
   !* ================================================================================= *!

   subroutine test_pv_hessian_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_pv_hessian_fd(FIX_PLAIN, LSF_SVDW, 1.0e-9_wp, error)
      if (allocated(error)) return
      call run_pv_hessian_fd(FIX_PLAIN, LSF_CFC, 1.0e-9_wp, error)
   end subroutine test_pv_hessian_fd

   subroutine test_pv_cross_hessian_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_pv_hessian_fd(FIX_CROSS, LSF_SVDW, 1.0e-8_wp, error)
   end subroutine test_pv_cross_hessian_fd

   subroutine run_pv_hessian_fd(fix_kind, lsf_kind, tol, error)
      integer, intent(in) :: fix_kind, lsf_kind
      real(wp), intent(in) :: tol
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: pressure = 2.5_wp, step = 2.5e-4_wp, cross_s = 2.0_wp
      type(structure_type) :: mol, mol_t
      type(moist_context_type), target :: ctx, ctx_t
      type(cavity_type_drop), allocatable :: cav, cav_t
      type(model_continuum_type) :: model, model_t
      type(model_continuum_component_pv) :: pv
      type(coupling_type), pointer :: cpl, cpl_t
      type(response_type) :: resp
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: hess(:, :, :, :), fd(:, :, :, :), grad_t(:, :)
      integer :: nat, iatom, iaxis, ioff
      real(wp) :: asym

      nullify (cpl, cpl_t)

      call drop_fixture_geometry(fix_kind, mol)
      nat = mol%nat
      call build_drop_test_cavity(cav, ctx, mol, fix_kind, lsf_kind, error, &
         & cross_branch_s=cross_s, want_fine=.true.)
      if (allocated(error)) return

      call new_continuum_model(model, ctx, cav, err)
      call new_component_pv(pv, pressure)
      if (.not. allocated(err)) call model%add_component(pv, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (.not. allocated(err)) call model%new_coupling(cpl, err)
      if (.not. allocated(err)) call model%prepare_gradient(cpl, err)
      call check_moist_error(error, err, "PV model setup")
      if (allocated(error)) return

      allocate (hess(3, nat, 3, nat), fd(3, nat, 3, nat), grad_t(3, nat), source=0.0_wp)
      call model%get_hessian(cpl, hess, err)
      call check_moist_error(error, err, "PV get_hessian")
      if (allocated(error)) return
      call model%release_coupling(cpl)

      asym = maxval(abs(reshape(hess, [3*nat, 3*nat]) - transpose(reshape(hess, [3*nat, 3*nat]))))
      call check(error, asym, 0.0_wp, thr=1.0e-11_wp*maxval(abs(hess)), more="PV Hessian symmetry")
      if (allocated(error)) return

      do iatom = 1, nat
         do iaxis = 1, 3
            do ioff = 1, 4
               mol_t = mol
               mol_t%xyz(iaxis, iatom) = mol%xyz(iaxis, iatom) + real(OFFSET(ioff), wp)*step
               call build_drop_test_cavity(cav_t, ctx_t, mol_t, fix_kind, lsf_kind, error, &
                  & cross_branch_s=cross_s, want_fine=.true.)
               if (allocated(error)) return
               call new_continuum_model(model_t, ctx_t, cav_t, err)
               if (.not. allocated(err)) call model_t%add_component(pv, err)
               if (.not. allocated(err)) call model_t%update(mol_t, err)
               if (.not. allocated(err)) call model_t%new_coupling(cpl_t, err)
               if (.not. allocated(err)) call model_t%prepare_gradient(cpl_t, err)
               grad_t = 0.0_wp
               if (.not. allocated(err)) call model_t%get_gradient(cpl_t, resp, grad_t, err)
               call check_moist_error(error, err, "PV displaced gradient")
               if (allocated(error)) return
               call model_t%release_coupling(cpl_t)
               fd(:, :, iaxis, iatom) = fd(:, :, iaxis, iatom) + COEFF(ioff)*grad_t/step
            end do
         end do
      end do

      call check_close(error, reshape(hess, [size(hess)]), reshape(fd, [size(fd)]), &
         & tol, tol, "PV dense Hessian vs FD")
   end subroutine run_pv_hessian_fd

   subroutine test_pv_hvp_matches_dense(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(structure_type) :: mol
      type(moist_context_type), target :: ctx
      type(cavity_type_drop), allocatable :: cav
      type(model_continuum_type) :: model
      type(model_continuum_component_pv) :: pv
      type(coupling_type), pointer :: cpl
      type(response_type) :: resp
      real(wp), allocatable :: hess(:, :, :, :), dirs(:, :, :), hvp(:, :, :), ref(:, :, :)
      integer :: nat, ndir, idir, iatom, iaxis

      nullify (cpl)

      call drop_fixture_geometry(FIX_PLAIN, mol)
      nat = mol%nat
      call build_drop_test_cavity(cav, ctx, mol, FIX_PLAIN, LSF_SVDW, error, want_fine=.true.)
      if (allocated(error)) return
      call new_continuum_model(model, ctx, cav, err)
      call new_component_pv(pv, 2.5_wp)
      if (.not. allocated(err)) call model%add_component(pv, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (.not. allocated(err)) call model%new_coupling(cpl, err)
      call check_moist_error(error, err, "PV setup")
      if (allocated(error)) return

      allocate (hess(3, nat, 3, nat), source=0.0_wp)
      call model%get_hessian(cpl, hess, err)
      call check_moist_error(error, err, "PV get_hessian")
      if (allocated(error)) return

      ! General directions (per-direction fixed channel)
      ndir = min(drop_hvp_per_dir_max, 3*nat - 1)
      allocate (dirs(3, nat, ndir), hvp(3, nat, ndir), ref(3, nat, ndir), source=0.0_wp)
      do idir = 1, ndir
         do iatom = 1, nat
            do iaxis = 1, 3
               dirs(iaxis, iatom, idir) = 0.4_wp*sin(0.5_wp*iaxis + 0.8_wp*iatom + idir)
               ref(:, :, idir) = ref(:, :, idir) + hess(:, :, iaxis, iatom)*dirs(iaxis, iatom, idir)
            end do
         end do
      end do

      call model%prepare_hvp_geometry(cpl, dirs, err)
      if (.not. allocated(err)) call check(error, .not. cpl%next(), more="PV stage 1 must be empty")
      if (allocated(error)) return
      if (.not. allocated(err)) call model%prepare_hvp(cpl, err)
      if (.not. allocated(err)) call check(error, .not. cpl%next(), more="PV stage 2 must be empty")
      if (allocated(error)) return
      if (.not. allocated(err)) call model%get_hvp(cpl, resp, hvp, err)
      call check_moist_error(error, err, "PV two-stage HVP")
      if (allocated(error)) return
      call check(error, .not. resp%next(), more="PV HVP must emit no response items")
      if (allocated(error)) return
      call check_close(error, reshape(hvp, [size(hvp)]), reshape(ref, [size(ref)]), &
         & 1.0e-12_wp, 1.0e-12_wp, "PV two-stage HVP vs dense")
      call model%release_coupling(cpl)
   end subroutine test_pv_hvp_matches_dense

   !* ================================================================================= *!
   !*                            PCM & Multi-component HVP                              *!
   !* ================================================================================= *!

   subroutine test_pcm_geometric_hvp_fd(error)
      type(error_type), allocatable, intent(out) :: error

      call run_pcm_model_hvp_fd(FIX_PLAIN, LSF_SVDW, pcm_cpcm, .false., 1.0e-9_wp, error)
      if (allocated(error)) return
      call run_pcm_model_hvp_fd(FIX_PLAIN, LSF_SVDW, pcm_cpcm, .true., 1.0e-9_wp, error)
      if (allocated(error)) return
      call run_pcm_model_hvp_fd(FIX_PLAIN, LSF_CFC, pcm_iefpcm, .false., 1.0e-9_wp, error)
      if (allocated(error)) return
      call run_pcm_model_hvp_fd(FIX_PLAIN, LSF_SVDW, pcm_dpcm, .false., 1.0e-8_wp, error)
      if (allocated(error)) return
      call run_pcm_model_hvp_fd(FIX_PLAIN, LSF_SVDW, pcm_slic, .false., 1.0e-8_wp, error)
      if (allocated(error)) return
      call run_pcm_model_hvp_fd(FIX_PLAIN, LSF_SVDW, pcm_nonlocal, .false., 1.0e-8_wp, error)
   end subroutine test_pcm_geometric_hvp_fd

   subroutine test_pcm_cross_hvp_fd(error)
      type(error_type), allocatable, intent(out) :: error

      call run_pcm_model_hvp_fd(FIX_CROSS, LSF_SVDW, pcm_cpcm, .false., 1.0e-8_wp, error)
      if (allocated(error)) return
      call run_pcm_model_hvp_fd(FIX_CROSS, LSF_SVDW, pcm_dpcm, .false., 2.0e-8_wp, error)
   end subroutine test_pcm_cross_hvp_fd

   !> Two-stage HVP of a single-PCM continuum model with nuclear-only charges `z` vs 4-pt FD
   subroutine run_pcm_model_hvp_fd(fix_kind, lsf_kind, pcm_kind, point_pot, tol, error)
      integer, intent(in) :: fix_kind, lsf_kind, pcm_kind
      logical, intent(in) :: point_pot
      real(wp), intent(in) :: tol
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: step = 2.0e-4_wp, cross_s = 2.0_wp
      integer, parameter :: ndir = 2
      type(structure_type) :: mol, mol_t
      type(moist_context_type), target :: ctx, ctx_t
      type(cavity_type_drop), allocatable :: cav, cav_t
      class(model_continuum_component_pcm), allocatable :: pcm
      type(model_continuum_type) :: model, model_t
      type(coupling_type), pointer :: cpl, cpl_t
      type(response_type) :: resp, resp_t
      class(response_item_type), allocatable :: item
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: z(:), dirs(:, :, :), hvp(:, :, :), fd_hvp(:, :, :), grad_t(:, :)
      real(wp), allocatable :: an_dw_phi(:, :), fd_dw_phi(:), w_phi_t(:)
      real(wp), allocatable :: an_dw_dphi_dr(:, :, :), fd_dw_dphi_dr(:, :), w_dphi_dr_t(:, :)
      integer, allocatable :: num_base(:), num_t(:)
      integer :: nat, ngrid, idir, ioff, iatom, iaxis, i, idx
      logical :: has_normal
      character(len=:), allocatable :: tag

      nullify (cpl, cpl_t)

      call drop_fixture_geometry(fix_kind, mol)
      nat = mol%nat
      allocate (z(nat), dirs(3, nat, ndir), hvp(3, nat, ndir), fd_hvp(3, nat, ndir), &
         & grad_t(3, nat), source=0.0_wp)
      do iatom = 1, nat
         z(iatom) = 0.8_wp*sin(1.3_wp*iatom + 0.4_wp)
      end do
      z = z - sum(z)/real(nat, wp)
      do idir = 1, ndir
         do iatom = 1, nat
            do iaxis = 1, 3
               dirs(iaxis, iatom, idir) = 0.35_wp*cos(0.45_wp*iaxis + 0.85_wp*iatom - idir)
            end do
         end do
      end do

      call build_drop_test_cavity(cav, ctx, mol, fix_kind, lsf_kind, error, &
         & cross_branch_s=cross_s, want_fine=.true.)
      if (allocated(error)) return
      call new_pcm(pcm_kind, ctx, pcm, error)
      if (allocated(error)) return
      pcm%param%point_potential = point_pot
      tag = trim(pcm%name)

      call new_continuum_model(model, ctx, cav, err)
      if (.not. allocated(err)) call model%add_component(pcm, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (.not. allocated(err)) call model%new_coupling(cpl, err)
      if (.not. allocated(err)) call model%prepare_hvp_geometry(cpl, dirs, err)
      call check_moist_error(error, err, tag//" stage 1")
      if (allocated(error)) return
      call check(error, .not. cpl%next(), more=tag//": geometric stage 1 must be empty")
      if (allocated(error)) return

      call model%prepare_hvp(cpl, err)
      call check_moist_error(error, err, tag//" stage 2")
      if (allocated(error)) return
      call answer_charge_model_coupling(model, cpl, mol, z, error)
      if (allocated(error)) return

      call model%get_hvp(cpl, resp, hvp, err)
      call check_moist_error(error, err, tag//" get_hvp")
      if (allocated(error)) return
      call check(error, resp%next(), more=tag//": missing potential_adjoint_tangent")
      if (allocated(error)) return
      item = resp%item()
      has_normal = .false.
      select type (item)
      type is (potential_adjoint_tangent_response_type)
         an_dw_phi = item%dw_phi
         has_normal = allocated(item%dw_dphi_dr)
         if (has_normal) an_dw_dphi_dr = item%dw_dphi_dr
      class default
         call test_failed(error, tag//": wrong response item")
         return
      end select
      select type (mcav => model%cavity)
      type is (cavity_type_drop)
         ngrid = mcav%ngrid
         num_base = mcav%numbering(1:ngrid)
      end select
      call model%release_coupling(cpl)

      allocate (fd_dw_phi(ngrid))
      if (has_normal) allocate (fd_dw_dphi_dr(3, ngrid))
      do idir = 1, ndir
         fd_dw_phi = 0.0_wp
         if (has_normal) fd_dw_dphi_dr = 0.0_wp
         do ioff = 1, 4
            mol_t = mol
            mol_t%xyz = mol%xyz + real(OFFSET(ioff), wp)*step*dirs(:, :, idir)
            call build_drop_test_cavity(cav_t, ctx_t, mol_t, fix_kind, lsf_kind, error, &
               & cross_branch_s=cross_s, want_fine=.true.)
            if (allocated(error)) return
            call new_continuum_model(model_t, ctx_t, cav_t, err)
            if (.not. allocated(err)) call model_t%add_component(pcm, err)
            if (.not. allocated(err)) call model_t%update(mol_t, err)
            if (.not. allocated(err)) call model_t%new_coupling(cpl_t, err)
            if (.not. allocated(err)) call model_t%prepare_gradient(cpl_t, err)
            call check_moist_error(error, err, tag//" displaced setup")
            if (allocated(error)) return
            call answer_charge_model_coupling(model_t, cpl_t, mol_t, z, error)
            if (allocated(error)) return

            grad_t = 0.0_wp
            call model_t%get_gradient(cpl_t, resp_t, grad_t, err)
            call check_moist_error(error, err, tag//" displaced gradient")
            if (allocated(error)) return
            fd_hvp(:, :, idir) = fd_hvp(:, :, idir) + COEFF(ioff)*grad_t/step

            call check(error, resp_t%next(), more=tag//": missing displaced potential_adjoint")
            if (allocated(error)) return
            item = resp_t%item()
            select type (item)
            type is (potential_adjoint_response_type)
               w_phi_t = item%w_phi
               if (has_normal) w_dphi_dr_t = item%w_dphi_dr
            end select
            select type (tcav => model_t%cavity)
            type is (cavity_type_drop)
               num_t = tcav%numbering(1:tcav%ngrid)
            end select
            do i = 1, ngrid
               idx = findloc(num_t, num_base(i), dim=1)
               if (idx > 0) then
                  fd_dw_phi(i) = fd_dw_phi(i) + COEFF(ioff)*w_phi_t(idx)/step
                  if (has_normal) then
                     fd_dw_dphi_dr(:, i) = fd_dw_dphi_dr(:, i) &
                        & + COEFF(ioff)*w_dphi_dr_t(:, idx)/step
                  end if
               end if
            end do
            call model_t%release_coupling(cpl_t)
         end do
         call check_close(error, an_dw_phi(:, idir), fd_dw_phi, 5.0e-7_wp, 5.0e-6_wp, &
            & tag//" model dw_phi vs FD")
         if (allocated(error)) return
         if (has_normal) then
            call check_close(error, reshape(an_dw_dphi_dr(:, :, idir), [3*ngrid]), &
               & reshape(fd_dw_dphi_dr, [3*ngrid]), 5.0e-7_wp, 5.0e-6_wp, &
               & tag//" model dw_dphi_dr vs FD")
            if (allocated(error)) return
         end if
      end do

      call check_close(error, reshape(hvp, [size(hvp)]), reshape(fd_hvp, [size(fd_hvp)]), &
         & tol, tol*10.0_wp, tag//" model HVP vs FD")
   end subroutine run_pcm_model_hvp_fd

   !> Multi-component `CPCM + DPCM + PV` continuum model HVP vs 4-point FD
   subroutine test_multi_component_hvp_fd(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      real(wp), parameter :: step = 2.0e-4_wp, pressure = 1.2_wp
      integer, parameter :: ndir = 2
      type(structure_type) :: mol, mol_t
      type(moist_context_type), target :: ctx, ctx_t
      type(cavity_type_drop), allocatable :: cav, cav_t
      class(model_continuum_component_pcm), allocatable :: cpcm, dpcm
      type(model_continuum_component_pv) :: pv
      type(model_continuum_type) :: model, model_t
      type(coupling_type), pointer :: cpl, cpl_t
      type(response_type) :: resp, resp_t
      real(wp), allocatable :: z(:), dirs(:, :, :), hvp(:, :, :), fd_hvp(:, :, :), grad_t(:, :)
      integer :: nat, idir, ioff, iatom, iaxis

      nullify (cpl, cpl_t)

      call drop_fixture_geometry(FIX_PLAIN, mol)
      nat = mol%nat
      z = [1.1_wp, -1.7_wp, 0.6_wp]
      allocate (dirs(3, nat, ndir), hvp(3, nat, ndir), fd_hvp(3, nat, ndir), &
         & grad_t(3, nat), source=0.0_wp)
      do idir = 1, ndir
         do iatom = 1, nat
            do iaxis = 1, 3
               dirs(iaxis, iatom, idir) = 0.3_wp*sin(0.6_wp*iaxis + 0.9_wp*iatom + idir)
            end do
         end do
      end do

      call build_drop_test_cavity(cav, ctx, mol, FIX_PLAIN, LSF_SVDW, error, want_fine=.true.)
      if (allocated(error)) return
      call new_pcm(pcm_cpcm, ctx, cpcm, error)
      if (allocated(error)) return
      call new_pcm(pcm_dpcm, ctx, dpcm, error)
      if (allocated(error)) return
      call new_component_pv(pv, pressure)

      call new_continuum_model(model, ctx, cav, err)
      if (.not. allocated(err)) call model%add_component(cpcm, err)
      if (.not. allocated(err)) call model%add_component(dpcm, err)
      if (.not. allocated(err)) call model%add_component(pv, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (.not. allocated(err)) call model%new_coupling(cpl, err)
      if (.not. allocated(err)) call model%prepare_hvp_geometry(cpl, dirs, err)
      if (.not. allocated(err)) call model%prepare_hvp(cpl, err)
      call check_moist_error(error, err, "multi-component HVP setup")
      if (allocated(error)) return
      call answer_charge_model_coupling(model, cpl, mol, z, error)
      if (allocated(error)) return
      call model%get_hvp(cpl, resp, hvp, err)
      call check_moist_error(error, err, "multi-component get_hvp")
      if (allocated(error)) return
      call model%release_coupling(cpl)

      do idir = 1, ndir
         do ioff = 1, 4
            mol_t = mol
            mol_t%xyz = mol%xyz + real(OFFSET(ioff), wp)*step*dirs(:, :, idir)
            call build_drop_test_cavity(cav_t, ctx_t, mol_t, FIX_PLAIN, LSF_SVDW, error, want_fine=.true.)
            if (allocated(error)) return
            call new_continuum_model(model_t, ctx_t, cav_t, err)
            if (.not. allocated(err)) call model_t%add_component(cpcm, err)
            if (.not. allocated(err)) call model_t%add_component(dpcm, err)
            if (.not. allocated(err)) call model_t%add_component(pv, err)
            if (.not. allocated(err)) call model_t%update(mol_t, err)
            if (.not. allocated(err)) call model_t%new_coupling(cpl_t, err)
            if (.not. allocated(err)) call model_t%prepare_gradient(cpl_t, err)
            call check_moist_error(error, err, "multi displaced setup")
            if (allocated(error)) return
            call answer_charge_model_coupling(model_t, cpl_t, mol_t, z, error)
            if (allocated(error)) return
            grad_t = 0.0_wp
            call model_t%get_gradient(cpl_t, resp_t, grad_t, err)
            call check_moist_error(error, err, "multi displaced gradient")
            if (allocated(error)) return
            fd_hvp(:, :, idir) = fd_hvp(:, :, idir) + COEFF(ioff)*grad_t/step
            call model_t%release_coupling(cpl_t)
         end do
      end do

      call check_close(error, reshape(hvp, [size(hvp)]), reshape(fd_hvp, [size(fd_hvp)]), &
         & 1.0e-8_wp, 1.0e-7_wp, "multi-component CPCM+DPCM+PV HVP vs FD")
   end subroutine test_multi_component_hvp_fd

   !* ================================================================================= *!
   !*                     Isodensity (Internal & Callback) + CPCM HVP                   *!
   !* ================================================================================= *!

   subroutine test_isodensity_hvp_fd(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: hvp_int(:, :, :), hvp_cb(:, :, :), fd_hvp(:, :, :)

      call eval_isodensity_hvp(.false., hvp_int, fd_hvp, .true., error)
      if (allocated(error)) return
      call check_close(error, reshape(hvp_int, [size(hvp_int)]), reshape(fd_hvp, [size(fd_hvp)]), &
         & 2.0e-8_wp, 2.0e-7_wp, "isodensity_internal CPCM HVP vs FD")
      if (allocated(error)) return

      call eval_isodensity_hvp(.true., hvp_cb, fd_hvp, .false., error)
      if (allocated(error)) return
      call check_close(error, reshape(hvp_cb, [size(hvp_cb)]), reshape(hvp_int, [size(hvp_int)]), &
         & 1.0e-12_wp, 1.0e-12_wp, "isodensity_callback matches isodensity_internal HVP")
   end subroutine test_isodensity_hvp_fd

   subroutine eval_isodensity_hvp(use_callback, total_hvp, fd_hvp, compute_fd, error)
      logical, intent(in) :: use_callback, compute_fd
      real(wp), allocatable, intent(out) :: total_hvp(:, :, :), fd_hvp(:, :, :)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: step = 2.0e-4_wp
      integer, parameter :: nat = 2, ndir = 2
      real(wp), parameter :: coords(3, nat) = reshape([0.0_wp, 0.0_wp, 0.0_wp, 2.4_wp, 0.3_wp, -0.2_wp], &
         & [3, nat])
      real(wp), parameter :: alpha(nat) = [0.28_wp, 0.32_wp]
      real(wp), parameter :: diag_d(nat) = [0.35_wp, 0.30_wp]
      real(wp), parameter :: z(nat) = [1.1_wp, -1.1_wp]

      type(structure_type) :: mol, mol_t
      type(s_density_ctx_type), target :: dctx, dctx_t
      type(moist_context_type), target :: ctx, ctx_t
      type(model_continuum_type) :: model, model_t
      type(coupling_type), pointer :: cpl, cpl_t
      type(response_type) :: resp, resp_r, resp_g
      class(response_item_type), allocatable :: item
      type(moist_error_type), allocatable :: err
      real(wp) :: dirs(3, nat, ndir), grad_t(3, nat)
      integer :: idir, ioff, iatom, iaxis

      nullify (cpl, cpl_t)

      call new(mol, [8, 6], coords)
      dctx%xyz = coords
      dctx%alpha = alpha
      dctx%diag_d = diag_d

      do idir = 1, ndir
         do iatom = 1, nat
            do iaxis = 1, 3
               dirs(iaxis, iatom, idir) = 0.3_wp*cos(0.7_wp*iaxis + 1.1_wp*iatom - idir)
            end do
         end do
      end do

      call build_iso_cpcm_model(use_callback, mol, dctx, ctx, model, cpl, error)
      if (allocated(error)) return

      ! Stage 1: density_tangent
      call model%prepare_hvp_geometry(cpl, dirs, err)
      call check_moist_error(error, err, "iso prepare_hvp_geometry")
      if (allocated(error)) return
      call answer_iso_hvp_coupling(model, cpl, mol, dctx, z, error)
      if (allocated(error)) return

      ! Stage 2: potential + potential_tangent
      call model%prepare_hvp(cpl, err)
      call check_moist_error(error, err, "iso prepare_hvp")
      if (allocated(error)) return
      call answer_iso_hvp_coupling(model, cpl, mol, dctx, z, error)
      if (allocated(error)) return

      allocate (total_hvp(3, nat, ndir), fd_hvp(3, nat, ndir), source=0.0_wp)
      call model%get_hvp(cpl, resp, total_hvp, err)
      call check_moist_error(error, err, "iso get_hvp")
      if (allocated(error)) return

      ! Add host density HVP contribution from `density_tangent` response item
      do while (resp%next())
         item = resp%item()
         select type (item)
         type is (density_tangent_response_type)
            call contract_s_density_hvp(model%cavity%xyz(:, 1:model%cavity%ngrid), dctx, dirs, &
               & item, total_hvp)
         end select
      end do
      call model%release_coupling(cpl)

      if (.not. compute_fd) return

      do idir = 1, ndir
         do ioff = 1, 4
            mol_t = mol
            mol_t%xyz = mol%xyz + real(OFFSET(ioff), wp)*step*dirs(:, :, idir)
            dctx_t = dctx
            dctx_t%xyz = mol_t%xyz
            call build_iso_cpcm_model(use_callback, mol_t, dctx_t, ctx_t, model_t, cpl_t, error)
            if (allocated(error)) return

            ! Response phase for density_response_type, then gradient phase for nuclear gradient
            call model_t%prepare_response(cpl_t, err)
            call check_moist_error(error, err, "iso displaced prepare_response")
            if (allocated(error)) return
            call answer_iso_hvp_coupling(model_t, cpl_t, mol_t, dctx_t, z, error)
            if (allocated(error)) return
            call model_t%get_response(cpl_t, resp_r, err)
            call check_moist_error(error, err, "iso displaced get_response")
            if (allocated(error)) return

            call model_t%prepare_gradient(cpl_t, err)
            call check_moist_error(error, err, "iso displaced prepare_gradient")
            if (allocated(error)) return
            call answer_iso_hvp_coupling(model_t, cpl_t, mol_t, dctx_t, z, error)
            if (allocated(error)) return
            grad_t = 0.0_wp
            call model_t%get_gradient(cpl_t, resp_g, grad_t, err)
            call check_moist_error(error, err, "iso displaced get_gradient")
            if (allocated(error)) return

            do while (resp_r%next())
               item = resp_r%item()
               select type (item)
               type is (density_response_type)
                  call contract_s_density_gradient(model_t%cavity%xyz(:, 1:model_t%cavity%ngrid), &
                     & dctx_t, item, grad_t)
               end select
            end do
            fd_hvp(:, :, idir) = fd_hvp(:, :, idir) + COEFF(ioff)*grad_t/step
            call model_t%release_coupling(cpl_t)
         end do
      end do
   end subroutine eval_isodensity_hvp

   subroutine build_iso_cpcm_model(use_callback, mol, dctx, ctx, model, cpl, error)
      logical, intent(in) :: use_callback
      type(structure_type), intent(in) :: mol
      type(s_density_ctx_type), target, intent(in) :: dctx
      type(moist_context_type), target, intent(out) :: ctx
      type(model_continuum_type), intent(out) :: model
      type(coupling_type), pointer, intent(out) :: cpl
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      type(cavity_type_drop) :: cav
      class(model_continuum_component_pcm), allocatable :: cpcm
      type(moist_cavity_drop_lsf_isodensity_param_type) :: iso_param
      type(moist_cavity_drop_parameters_type) :: drop_param
      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf_int
      type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf_cb
      real(wp) :: dcart(2, 2)

      call new_context(ctx, verbosity=0)
      iso_param%rho_iso = 1.0e-3_wp
      iso_param%scale = 1.0e3_wp
      iso_param%core_fill = .false.
      drop_param = moist_cavity_drop_parameters_type(num_leb=26, do_fine=.true., &
         & tolerance=1.0e-13_wp, proj_maxiter=150, proj_level=2)

      if (use_callback) then
         call lsf_cb%new(c_funloc(s_density_callback), c_loc(dctx), param=iso_param)
         call new_cavity_drop(cav, radius_model=default_cpcm_radii(), lsf_model=lsf_cb, &
            & error=err, param=drop_param, ctx=ctx)
      else
         call lsf_int%new([1, 2], [0, 0], [1, 1], dctx%alpha, [1.0_wp, 1.0_wp], err, param=iso_param)
         if (.not. allocated(err)) then
            call new_cavity_drop(cav, radius_model=default_cpcm_radii(), lsf_model=lsf_int, &
               & error=err, param=drop_param, ctx=ctx)
         end if
      end if
      call check_moist_error(error, err, "iso cavity init")
      if (allocated(error)) return

      call new_pcm(pcm_cpcm, ctx, cpcm, error)
      if (allocated(error)) return
      call new_continuum_model(model, ctx, cav, err)
      if (.not. allocated(err)) call model%add_component(cpcm, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (.not. allocated(err)) call model%new_coupling(cpl, err)
      if (.not. allocated(err)) call model%prepare_build(cpl, err)
      call check_moist_error(error, err, "iso model prepare_build")
      if (allocated(error)) return

      dcart = 0.0_wp
      dcart(1, 1) = dctx%diag_d(1)
      dcart(2, 2) = dctx%diag_d(2)
      do while (cpl%next())
         select type (req => cpl%request())
         type is (density_matrix_request_type)
            call submit(cpl, "dcart", dcart)
         type is (nuclear_charge_request_type)
            call submit(cpl, "z", real(mol%num(mol%id), wp))
         end select
      end do
      call model%build(cpl, err)
      call check_moist_error(error, err, "iso model build")
   end subroutine build_iso_cpcm_model

   !> Pure density and charge perturbations (`dirs = 0`, matching CPHF iterations) vs 4-pt FD
   subroutine test_pure_density_and_charge_hvp_fd(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: step = 2.0e-4_wp
      integer, parameter :: nat = 2, ndir = 2
      real(wp), parameter :: coords(3, nat) = reshape([0.0_wp, 0.0_wp, 0.0_wp, &
         & 2.4_wp, 0.3_wp, -0.2_wp], [3, nat])
      real(wp), parameter :: alpha(nat) = [0.28_wp, 0.32_wp]
      real(wp), parameter :: diag_d(nat) = [0.35_wp, 0.30_wp]
      real(wp), parameter :: z(nat) = [1.1_wp, -1.1_wp]
      real(wp), parameter :: d_diag_d(nat, ndir) = reshape([0.08_wp, -0.05_wp, &
         & 0.04_wp, 0.06_wp], [nat, ndir])

      type(structure_type) :: mol
      type(s_density_ctx_type), target :: dctx, dctx_t, dctx_dir
      type(moist_context_type), target :: ctx, ctx_t
      type(model_continuum_type) :: model, model_t
      type(coupling_type), pointer :: cpl, cpl_t
      type(response_type) :: resp, resp_r, resp_g
      class(response_item_type), allocatable :: item
      type(moist_error_type), allocatable :: err
      type(density_response_type) :: dens_item
      real(wp) :: dirs(3, nat, ndir), total_hvp(3, nat, ndir), fd_hvp(3, nat, ndir)
      real(wp) :: grad_t(3, nat), t
      real(wp), allocatable :: an_dw_phi(:, :), fd_dw_phi(:), w_phi_t(:)
      real(wp), allocatable :: d_rho(:, :), d_grad_rho(:, :, :), d_hess_rho(:, :, :), d_d3_rho(:, :, :)
      real(wp), allocatable :: d_phi(:, :), d_dphi_dr(:, :, :), d_dphi_dxi(:, :)
      real(wp), allocatable :: phi(:), dphi_dr(:, :), dphi_dxi(:)
      integer, allocatable :: num_base(:), num_t(:)
      integer :: idir, ioff, ngrid, i, idx

      nullify (cpl, cpl_t)
      call new(mol, [8, 6], coords)
      dctx%xyz = coords
      dctx%alpha = alpha
      dctx%diag_d = diag_d
      dirs = 0.0_wp

      call build_iso_cpcm_model(.false., mol, dctx, ctx, model, cpl, error)
      if (allocated(error)) return

      ! Stage 1: density_tangent from d_diag_d at fixed surface (dirs = 0)
      call model%prepare_hvp_geometry(cpl, dirs, err)
      call check_moist_error(error, err, "pure density prepare_hvp_geometry")
      if (allocated(error)) return
      do while (cpl%next())
         select type (req => cpl%request())
         type is (density_tangent_request_type)
            ngrid = size(req%xyz, 2)
            allocate (d_rho(ngrid, ndir), d_grad_rho(3, ngrid, ndir), &
               & d_hess_rho(6, ngrid, ndir), d_d3_rho(10, ngrid, ndir))
            call eval_s_density_param_tangent(req%xyz, dctx, d_diag_d, d_rho, d_grad_rho, &
               & d_hess_rho, d_d3_rho)
            call submit(cpl, "d_rho", d_rho)
            call submit(cpl, "d_grad_rho", d_grad_rho)
            call submit(cpl, "d_hess_rho", d_hess_rho)
            call submit(cpl, "d_d3_rho", d_d3_rho)
            deallocate (d_rho, d_grad_rho, d_hess_rho, d_d3_rho)
         end select
      end do

      ! Stage 2: primal potential + potential_tangent (surface motion + dz)
      call model%prepare_hvp(cpl, err)
      call check_moist_error(error, err, "pure density prepare_hvp")
      if (allocated(error)) return
      ngrid = model%cavity%ngrid
      do while (cpl%next())
         select type (req => cpl%request())
         type is (gaussian_potential_request_type)
            allocate (phi(ngrid), dphi_dr(3, ngrid), dphi_dxi(ngrid))
            call pcm_elstat_potential_derivatives(model%cavity%xyz(:, 1:ngrid), mol%xyz, z, &
               & phi, dphi_dr, err, xi=model%cavity%xi0(1:ngrid), dphi_dxi=dphi_dxi)
            call check_moist_error(error, err, "pure density gaussian_potential")
            if (allocated(error)) return
            call submit(cpl, "phi", phi)
            call submit(cpl, "dphi_dr", dphi_dr)
            call submit(cpl, "dphi_dxi", dphi_dxi)
            deallocate (phi, dphi_dr, dphi_dxi)
         type is (gaussian_potential_tangent_request_type)
            allocate (d_phi(ngrid, ndir), d_dphi_dr(3, ngrid, ndir), d_dphi_dxi(ngrid, ndir))
            call pcm_elstat_potential_tangent(req%xyz, mol%xyz, z, req%d_xyz, req%dirs, &
               & d_phi, d_dphi_dr, err, xi=req%xi, d_xi=req%d_xi, d_dphi_dxi=d_dphi_dxi)
            call check_moist_error(error, err, "pure density potential_tangent")
            if (allocated(error)) return
            call submit(cpl, "d_phi", d_phi)
            call submit(cpl, "d_dphi_dr", d_dphi_dr)
            call submit(cpl, "d_dphi_dxi", d_dphi_dxi)
            deallocate (d_phi, d_dphi_dr, d_dphi_dxi)
         type is (nuclear_charge_request_type)
            call submit(cpl, "z", z)
         end select
      end do

      total_hvp = 0.0_wp
      call model%get_hvp(cpl, resp, total_hvp, err)
      call check_moist_error(error, err, "pure density get_hvp")
      if (allocated(error)) return

      do while (resp%next())
         item = resp%item()
         select type (item)
         type is (potential_adjoint_tangent_response_type)
            an_dw_phi = item%dw_phi
         type is (density_tangent_response_type)
            call contract_s_density_hvp(model%cavity%xyz(:, 1:ngrid), dctx, dirs, item, total_hvp)
            ! Add host density derivative w.r.t. d_diag_d at fixed geometry
            dens_item%w_rho = item%w_rho
            dens_item%w_grad_rho = item%w_grad_rho
            dens_item%w_hess_rho = item%w_hess_rho
            dctx_dir = dctx
            do idir = 1, ndir
               dctx_dir%diag_d = d_diag_d(:, idir)
               call contract_s_density_gradient(model%cavity%xyz(:, 1:ngrid), dctx_dir, &
                  & dens_item, total_hvp(:, :, idir))
            end do
         end select
      end do
      select type (mcav => model%cavity)
      type is (cavity_type_drop)
         num_base = mcav%numbering(1:ngrid)
      end select
      call model%release_coupling(cpl)

      fd_hvp = 0.0_wp
      allocate (fd_dw_phi(ngrid))
      do idir = 1, ndir
         fd_dw_phi = 0.0_wp
         do ioff = 1, 4
            t = real(OFFSET(ioff), wp)*step
            dctx_t = dctx
            dctx_t%diag_d = diag_d + t*d_diag_d(:, idir)
            call build_iso_cpcm_model(.false., mol, dctx_t, ctx_t, model_t, cpl_t, error)
            if (allocated(error)) return

            call model_t%prepare_response(cpl_t, err)
            if (.not. allocated(err)) call answer_iso_hvp_coupling(model_t, cpl_t, mol, dctx_t, z, error)
            if (allocated(error)) return
            if (.not. allocated(err)) call model_t%get_response(cpl_t, resp_r, err)
            call check_moist_error(error, err, "pure density displaced get_response")
            if (allocated(error)) return

            call model_t%prepare_gradient(cpl_t, err)
            if (.not. allocated(err)) call answer_iso_hvp_coupling(model_t, cpl_t, mol, dctx_t, z, error)
            if (allocated(error)) return
            grad_t = 0.0_wp
            if (.not. allocated(err)) call model_t%get_gradient(cpl_t, resp_g, grad_t, err)
            call check_moist_error(error, err, "pure density displaced get_gradient")
            if (allocated(error)) return

            do while (resp_r%next())
               item = resp_r%item()
               select type (item)
               type is (density_response_type)
                  call contract_s_density_gradient(model_t%cavity%xyz(:, 1:model_t%cavity%ngrid), &
                     & dctx_t, item, grad_t)
               type is (potential_adjoint_response_type)
                  w_phi_t = item%w_phi
               end select
            end do
            fd_hvp(:, :, idir) = fd_hvp(:, :, idir) + COEFF(ioff)*grad_t/step

            select type (tcav => model_t%cavity)
            type is (cavity_type_drop)
               num_t = tcav%numbering(1:tcav%ngrid)
            end select
            do i = 1, ngrid
               idx = findloc(num_t, num_base(i), dim=1)
               if (idx > 0) fd_dw_phi(i) = fd_dw_phi(i) + COEFF(ioff)*w_phi_t(idx)/step
            end do
            call model_t%release_coupling(cpl_t)
         end do
         call check_close(error, an_dw_phi(:, idir), fd_dw_phi, 5.0e-7_wp, 5.0e-6_wp, &
            & "pure density/charge dw_phi vs FD")
         if (allocated(error)) return
      end do
      call check_close(error, reshape(total_hvp, [size(total_hvp)]), &
         & reshape(fd_hvp, [size(fd_hvp)]), 2.0e-8_wp, 2.0e-7_wp, "pure density/charge HVP vs FD")
   end subroutine test_pure_density_and_charge_hvp_fd

   !> Multi-block chunking (`ndir = 194 > drop_hvp_chunk_dirs = 192`) matches single-dir HVPs
   subroutine test_pcm_multiblock_chunk_hvp(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nat = 2, ndir = 194
      integer, parameter :: probe(4) = [1, 65, 193, 194]
      real(wp), parameter :: coords(3, nat) = reshape([0.0_wp, 0.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 3.2_wp], [3, nat])
      real(wp), parameter :: z(nat) = [1.2_wp, -1.2_wp]

      type(structure_type) :: mol
      type(moist_context_type), target :: ctx
      type(cavity_type_drop), allocatable :: cav
      class(model_continuum_component_pcm), allocatable :: cpcm
      type(model_continuum_type) :: model
      type(coupling_type), pointer :: cpl
      type(response_type) :: resp_blk, resp_one
      class(response_item_type), allocatable :: item
      type(moist_error_type), allocatable :: err
      real(wp) :: dirs(3, nat, ndir), hvp_blk(3, nat, ndir), hvp_one(3, nat, 1)
      real(wp), allocatable :: dw_blk(:, :), dw_one(:, :)
      integer :: idir, iatom, iaxis, k

      nullify (cpl)
      call new(mol, [8, 1], coords)
      do idir = 1, ndir
         do iatom = 1, nat
            do iaxis = 1, 3
               dirs(iaxis, iatom, idir) = 0.25_wp*sin(0.37_wp*idir + 0.9_wp*iatom + 0.5_wp*iaxis)
            end do
         end do
      end do

      call build_drop_test_cavity(cav, ctx, mol, FIX_PLAIN, LSF_SVDW, error, want_fine=.true.)
      if (allocated(error)) return
      call new_pcm(pcm_cpcm, ctx, cpcm, error)
      if (allocated(error)) return
      call new_continuum_model(model, ctx, cav, err)
      if (.not. allocated(err)) call model%add_component(cpcm, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (.not. allocated(err)) call model%new_coupling(cpl, err)
      call check_moist_error(error, err, "multiblock setup")
      if (allocated(error)) return

      call model%prepare_hvp_geometry(cpl, dirs, err)
      if (.not. allocated(err)) call model%prepare_hvp(cpl, err)
      call check_moist_error(error, err, "multiblock prepare_hvp")
      if (allocated(error)) return
      call answer_charge_model_coupling(model, cpl, mol, z, error)
      if (allocated(error)) return
      hvp_blk = 0.0_wp
      call model%get_hvp(cpl, resp_blk, hvp_blk, err)
      call check_moist_error(error, err, "multiblock get_hvp")
      if (allocated(error)) return
      call check(error, resp_blk%next(), more="multiblock missing potential_adjoint_tangent")
      if (allocated(error)) return
      item = resp_blk%item()
      select type (item)
      type is (potential_adjoint_tangent_response_type)
         dw_blk = item%dw_phi
      class default
         call test_failed(error, "multiblock wrong response item")
         return
      end select

      do k = 1, size(probe)
         idir = probe(k)
         call model%prepare_hvp_geometry(cpl, dirs(:, :, idir:idir), err)
         if (.not. allocated(err)) call model%prepare_hvp(cpl, err)
         call check_moist_error(error, err, "single probe prepare_hvp")
         if (allocated(error)) return
         call answer_charge_model_coupling(model, cpl, mol, z, error)
         if (allocated(error)) return
         hvp_one = 0.0_wp
         call model%get_hvp(cpl, resp_one, hvp_one, err)
         call check_moist_error(error, err, "single probe get_hvp")
         if (allocated(error)) return
         call check(error, resp_one%next(), more="single probe missing potential_adjoint_tangent")
         if (allocated(error)) return
         item = resp_one%item()
         select type (item)
         type is (potential_adjoint_tangent_response_type)
            dw_one = item%dw_phi
         end select
         call check_close(error, reshape(hvp_blk(:, :, idir), [3*nat]), &
            & reshape(hvp_one(:, :, 1), [3*nat]), 1.0e-12_wp, 1.0e-12_wp, "multiblock hvp parity")
         if (allocated(error)) return
         call check_close(error, dw_blk(:, idir), dw_one(:, 1), &
            & 1.0e-12_wp, 1.0e-12_wp, "multiblock dw_phi parity")
         if (allocated(error)) return
      end do
      call model%release_coupling(cpl)
   end subroutine test_pcm_multiblock_chunk_hvp

   !* ================================================================================= *!
   !*                                    Guards                                         *!
   !* ================================================================================= *!

   subroutine test_hessian_guards(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(moist_context_type), target :: ctx
      type(cavity_type_drop), allocatable :: cav
      type(cavity_type_iswig) :: iswig
      class(model_continuum_component_pcm), allocatable :: cpcm
      type(model_continuum_component_pv) :: pv
      type(model_continuum_component_gostshyp) :: gost
      type(model_continuum_type) :: m_cpcm, m_gost, m_iswig, m_fwd
      type(coupling_type), pointer :: cpl
      type(response_type) :: resp
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :), hess(:, :, :, :)
      integer :: nat

      nullify (cpl)

      call drop_fixture_geometry(FIX_PLAIN, mol)
      nat = mol%nat
      allocate (dirs(3, nat, 1), hvp(3, nat, 1), hess(3, nat, 3, nat))
      dirs = 1.0_wp
      hvp = 7.0_wp
      hess = 9.0_wp

      call build_drop_test_cavity(cav, ctx, mol, FIX_PLAIN, LSF_SVDW, error, want_fine=.true.)
      if (allocated(error)) return
      call new_pcm(pcm_cpcm, ctx, cpcm, error)
      if (allocated(error)) return

      ! Dense get_hessian refused on PCM; ordering guards on two-stage HVP
      call new_continuum_model(m_cpcm, ctx, cav, err)
      if (.not. allocated(err)) call m_cpcm%add_component(cpcm, err)
      if (.not. allocated(err)) call m_cpcm%update(mol, err)
      if (.not. allocated(err)) call m_cpcm%new_coupling(cpl, err)
      call check_moist_error(error, err, "guard setup")
      if (allocated(error)) return

      call m_cpcm%get_hessian(cpl, hess, err)
      call check(error, allocated(err) .and. all(hess == 9.0_wp), &
         & more="PCM get_hessian must refuse cleanly without touching hessian")
      if (allocated(error)) return
      deallocate (err)

      call m_cpcm%prepare_hvp(cpl, err)
      call check(error, allocated(err), more="prepare_hvp before prepare_hvp_geometry must fail")
      if (allocated(error)) return
      deallocate (err)

      call m_cpcm%get_hvp(cpl, resp, hvp, err)
      call check(error, allocated(err) .and. all(hvp == 7.0_wp), &
         & more="get_hvp before prepare_hvp must fail cleanly")
      if (allocated(error)) return
      deallocate (err)
      call m_cpcm%release_coupling(cpl)

      ! GOSTSHYP component refused
      call new_component_gostshyp(gost, 1.0e-5_wp)
      call new_continuum_model(m_gost, ctx, cav, err)
      if (.not. allocated(err)) call m_gost%add_component(gost, err)
      if (.not. allocated(err)) call m_gost%update(mol, err)
      if (.not. allocated(err)) call m_gost%new_coupling(cpl, err)
      call check_moist_error(error, err, "gostshyp setup")
      if (allocated(error)) return
      call m_gost%prepare_hvp_geometry(cpl, dirs, err)
      call check(error, allocated(err), more="GOSTSHYP HVP must refuse cleanly")
      if (allocated(error)) return
      deallocate (err)
      call m_gost%release_coupling(cpl)

      ! iSwiG cavity refused
      call new_component_pv(pv, 1.0_wp)
      call new_cavity_iswig(iswig, radius_model=default_cpcm_radii(), error=err, &
         & param=moist_cavity_iswig_parameters_type(num_leb=26), ctx=ctx)
      if (.not. allocated(err)) call new_continuum_model(m_iswig, ctx, iswig, err)
      if (.not. allocated(err)) call m_iswig%add_component(pv, err)
      if (.not. allocated(err)) call m_iswig%update(mol, err)
      if (.not. allocated(err)) call m_iswig%new_coupling(cpl, err)
      call check_moist_error(error, err, "iSwiG model setup")
      if (allocated(error)) return
      call m_iswig%get_hessian(cpl, hess, err)
      call check(error, allocated(err) .and. all(hess == 9.0_wp), more="iSwiG Hessian must refuse")
      if (allocated(error)) return
      deallocate (err)
      call m_iswig%release_coupling(cpl)

      ! force_forward_gradient refused
      call new_continuum_model(m_fwd, ctx, cav, err)
      if (.not. allocated(err)) call m_fwd%add_component(pv, err)
      if (.not. allocated(err)) call m_fwd%update(mol, err)
      call m_fwd%use_forward_gradient(.true.)
      if (.not. allocated(err)) call m_fwd%new_coupling(cpl, err)
      call check_moist_error(error, err, "forward model setup")
      if (allocated(error)) return
      call m_fwd%get_hessian(cpl, hess, err)
      call check(error, allocated(err) .and. all(hess == 9.0_wp), &
         & more="forward-gradient model Hessian must refuse")
      call m_fwd%release_coupling(cpl)
   end subroutine test_hessian_guards

   !* ================================================================================= *!
   !*                            Analytic Host Helpers                                  *!
   !* ================================================================================= *!

   subroutine answer_charge_model_coupling(model, coupling, mol, z, error)
      type(model_continuum_type), intent(in) :: model
      type(coupling_type), intent(inout), target :: coupling
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: z(:)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: phi(:), dphi_dr(:, :), dphi_dxi(:), d2phi_dr2(:, :, :), d2phi_dr_dxi(:, :)
      real(wp), allocatable :: d_phi(:, :), d_dphi_dr(:, :, :), d_dphi_dxi(:, :)
      real(wp), allocatable :: d_d2phi_dr2(:, :, :), d_d2phi_dr_dxi(:, :, :)
      integer :: ngrid, ndir

      ngrid = model%cavity%ngrid
      do while (coupling%next())
         select type (req => coupling%request())
         type is (gaussian_potential_request_type)
            allocate (phi(ngrid), dphi_dr(3, ngrid), dphi_dxi(ngrid), &
               & d2phi_dr2(3, 3, ngrid), d2phi_dr_dxi(3, ngrid))
            call pcm_elstat_potential_derivatives(model%cavity%xyz(:, 1:ngrid), mol%xyz, z, &
               & phi, dphi_dr, err, xi=model%cavity%xi0(1:ngrid), dphi_dxi=dphi_dxi, &
               & d2phi_dr2=d2phi_dr2, d2phi_dr_dxi=d2phi_dr_dxi)
            call check_moist_error(error, err, "model gaussian_potential")
            if (allocated(error)) return
            call submit(coupling, "phi", phi)
            call submit(coupling, "dphi_dr", dphi_dr)
            call submit(coupling, "dphi_dxi", dphi_dxi)
            if (req%is_missing("d2phi_dr2")) call submit(coupling, "d2phi_dr2", d2phi_dr2)
            if (req%is_missing("d2phi_dr_dxi")) call submit(coupling, "d2phi_dr_dxi", d2phi_dr_dxi)
            deallocate (phi, dphi_dr, dphi_dxi, d2phi_dr2, d2phi_dr_dxi)
         type is (point_potential_request_type)
            allocate (phi(ngrid), dphi_dr(3, ngrid), d2phi_dr2(3, 3, ngrid))
            call pcm_elstat_potential_derivatives(model%cavity%xyz(:, 1:ngrid), mol%xyz, z, &
               & phi, dphi_dr, err, d2phi_dr2=d2phi_dr2)
            call check_moist_error(error, err, "model point_potential")
            if (allocated(error)) return
            call submit(coupling, "phi", phi)
            call submit(coupling, "dphi_dr", dphi_dr)
            if (req%is_missing("d2phi_dr2")) call submit(coupling, "d2phi_dr2", d2phi_dr2)
            deallocate (phi, dphi_dr, d2phi_dr2)
         type is (gaussian_potential_tangent_request_type)
            ndir = size(req%dirs, 3)
            allocate (d_phi(ngrid, ndir), d_dphi_dr(3, ngrid, ndir), d_dphi_dxi(ngrid, ndir), &
               & d_d2phi_dr2(6, ngrid, ndir), d_d2phi_dr_dxi(3, ngrid, ndir))
            call pcm_elstat_potential_tangent(req%xyz, mol%xyz, z, req%d_xyz, req%dirs, &
               & d_phi, d_dphi_dr, err, xi=req%xi, d_xi=req%d_xi, d_dphi_dxi=d_dphi_dxi, &
               & d_d2phi_dr2=d_d2phi_dr2, d_d2phi_dr_dxi=d_d2phi_dr_dxi)
            call check_moist_error(error, err, "model gaussian_potential_tangent")
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
            call pcm_elstat_potential_tangent(req%xyz, mol%xyz, z, req%d_xyz, req%dirs, &
               & d_phi, d_dphi_dr, err, d_d2phi_dr2=d_d2phi_dr2)
            call check_moist_error(error, err, "model point_potential_tangent")
            if (allocated(error)) return
            call submit(coupling, "d_phi", d_phi)
            call submit(coupling, "d_dphi_dr", d_dphi_dr)
            if (req%is_missing("d_d2phi_dr2")) call submit(coupling, "d_d2phi_dr2", d_d2phi_dr2)
            deallocate (d_phi, d_dphi_dr, d_d2phi_dr2)
         type is (nuclear_charge_request_type)
            call submit(coupling, "z", z)
         end select
      end do
   end subroutine answer_charge_model_coupling

   subroutine answer_iso_hvp_coupling(model, coupling, mol, dctx, z, error)
      type(model_continuum_type), intent(in) :: model
      type(coupling_type), intent(inout), target :: coupling
      type(structure_type), intent(in) :: mol
      type(s_density_ctx_type), intent(in) :: dctx
      real(wp), intent(in) :: z(:)
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: d_rho(:, :), d_grad_rho(:, :, :), d_hess_rho(:, :, :), d_d3_rho(:, :, :)
      integer :: ngrid, ndir

      do while (coupling%next())
         select type (req => coupling%request())
         type is (density_tangent_request_type)
            ngrid = size(req%xyz, 2)
            ndir = size(req%dirs, 3)
            allocate (d_rho(ngrid, ndir), d_grad_rho(3, ngrid, ndir), &
               & d_hess_rho(6, ngrid, ndir), d_d3_rho(10, ngrid, ndir))
            call eval_s_density_tangent(req%xyz, dctx, req%dirs, d_rho, d_grad_rho, &
               & d_hess_rho, d_d3_rho)
            call submit(coupling, "d_rho", d_rho)
            call submit(coupling, "d_grad_rho", d_grad_rho)
            call submit(coupling, "d_hess_rho", d_hess_rho)
            call submit(coupling, "d_d3_rho", d_d3_rho)
            deallocate (d_rho, d_grad_rho, d_hess_rho, d_d3_rho)
         end select
      end do
      call answer_charge_model_coupling(model, coupling, mol, z, error)
   end subroutine answer_iso_hvp_coupling

   !> 1D Hermite-Gaussian factors for exp(-beta * u^2) up to order 4
   pure subroutine gaussian_1d_polys(beta, s, poly)
      real(wp), intent(in) :: beta, s(3)
      real(wp), intent(out) :: poly(0:4, 3)
      integer :: ax
      real(wp) :: u, u2, b2

      b2 = beta*beta
      do ax = 1, 3
         u = s(ax)
         u2 = u*u
         poly(0, ax) = 1.0_wp
         poly(1, ax) = -2.0_wp*beta*u
         poly(2, ax) = 4.0_wp*b2*u2 - 2.0_wp*beta
         poly(3, ax) = -8.0_wp*b2*beta*u2*u + 12.0_wp*b2*u
         poly(4, ax) = 16.0_wp*b2*b2*u2*u2 - 48.0_wp*b2*beta*u2 + 12.0_wp*b2
      end do
   end subroutine gaussian_1d_polys

   pure real(wp) function lobe_deriv(pref, poly, axes) result(val)
      real(wp), intent(in) :: pref, poly(0:4, 3)
      integer, intent(in) :: axes(:)
      integer :: cnt(3), j

      cnt = 0
      do j = 1, size(axes)
         cnt(axes(j)) = cnt(axes(j)) + 1
      end do
      val = pref*poly(cnt(1), 1)*poly(cnt(2), 2)*poly(cnt(3), 3)
   end function lobe_deriv

   function s_density_callback(ctx_ptr, point, rho, drho, d2rho_ptr, d3rho_ptr, d4rho_ptr) &
         result(status) bind(C)
      type(c_ptr), value :: ctx_ptr, d2rho_ptr, d3rho_ptr, d4rho_ptr
      real(c_double), intent(in) :: point(3)
      real(c_double), intent(out) :: rho, drho(3)
      integer(c_int) :: status

      type(s_density_ctx_type), pointer :: dctx
      real(c_double), pointer :: d2(:, :), d3(:, :, :), d4(:, :, :, :)
      real(wp) :: poly(0:4, 3), s(3), beta, pref
      integer :: iat, a, b, c, d

      status = 0_c_int
      call c_f_pointer(ctx_ptr, dctx)
      rho = 0.0_c_double
      drho = 0.0_c_double
      if (c_associated(d2rho_ptr)) then
         call c_f_pointer(d2rho_ptr, d2, [3, 3])
         d2 = 0.0_c_double
      end if
      if (c_associated(d3rho_ptr)) then
         call c_f_pointer(d3rho_ptr, d3, [3, 3, 3])
         d3 = 0.0_c_double
      end if
      if (c_associated(d4rho_ptr)) then
         call c_f_pointer(d4rho_ptr, d4, [3, 3, 3, 3])
         d4 = 0.0_c_double
      end if

      do iat = 1, size(dctx%alpha)
         beta = 2.0_wp*dctx%alpha(iat)
         s = real(point, wp) - dctx%xyz(:, iat)
         pref = dctx%diag_d(iat)*exp(-beta*sum(s*s))
         call gaussian_1d_polys(beta, s, poly)
         rho = rho + pref
         do a = 1, 3
            drho(a) = drho(a) + lobe_deriv(pref, poly, [a])
         end do
         if (c_associated(d2rho_ptr)) then
            do b = 1, 3
               do a = 1, 3
                  d2(a, b) = d2(a, b) + lobe_deriv(pref, poly, [a, b])
               end do
            end do
         end if
         if (c_associated(d3rho_ptr)) then
            do c = 1, 3
               do b = 1, 3
                  do a = 1, 3
                     d3(a, b, c) = d3(a, b, c) + lobe_deriv(pref, poly, [a, b, c])
                  end do
               end do
            end do
         end if
         if (c_associated(d4rho_ptr)) then
            do d = 1, 3
               do c = 1, 3
                  do b = 1, 3
                     do a = 1, 3
                        d4(a, b, c, d) = d4(a, b, c, d) + lobe_deriv(pref, poly, [a, b, c, d])
                     end do
                  end do
               end do
            end do
         end if
      end do
   end function s_density_callback

   pure subroutine eval_s_density_tangent(xyz, dctx, dirs, d_rho, d_grad, d_hess, d_d3)
      real(wp), intent(in) :: xyz(:, :), dirs(:, :, :)
      type(s_density_ctx_type), intent(in) :: dctx
      real(wp), intent(out) :: d_rho(:, :), d_grad(:, :, :), d_hess(:, :, :), d_d3(:, :, :)

      real(wp) :: poly(0:4, 3), s(3), beta, pref, vk
      integer :: igrid, idir, iat, k, a, idx, p, q, r

      d_rho = 0.0_wp
      d_grad = 0.0_wp
      d_hess = 0.0_wp
      d_d3 = 0.0_wp
      do igrid = 1, size(xyz, 2)
         do iat = 1, size(dctx%alpha)
            beta = 2.0_wp*dctx%alpha(iat)
            s = xyz(:, igrid) - dctx%xyz(:, iat)
            pref = dctx%diag_d(iat)*exp(-beta*sum(s*s))
            call gaussian_1d_polys(beta, s, poly)
            do idir = 1, size(dirs, 3)
               do k = 1, 3
                  vk = dirs(k, iat, idir)
                  d_rho(igrid, idir) = d_rho(igrid, idir) - vk*lobe_deriv(pref, poly, [k])
                  do a = 1, 3
                     d_grad(a, igrid, idir) = d_grad(a, igrid, idir) &
                        & - vk*lobe_deriv(pref, poly, [a, k])
                  end do
                  do idx = 1, 6
                     p = sym2_idx(1, idx)
                     q = sym2_idx(2, idx)
                     d_hess(idx, igrid, idir) = d_hess(idx, igrid, idir) &
                        & - vk*lobe_deriv(pref, poly, [p, q, k])
                  end do
                  do idx = 1, 10
                     p = sym3_idx(1, idx)
                     q = sym3_idx(2, idx)
                     r = sym3_idx(3, idx)
                     d_d3(idx, igrid, idir) = d_d3(idx, igrid, idir) &
                        & - vk*lobe_deriv(pref, poly, [p, q, r, k])
                  end do
               end do
            end do
         end do
      end do
   end subroutine eval_s_density_tangent

   pure subroutine eval_s_density_param_tangent(xyz, dctx, d_diag_d, d_rho, d_grad, d_hess, d_d3)
      real(wp), intent(in) :: xyz(:, :), d_diag_d(:, :)
      type(s_density_ctx_type), intent(in) :: dctx
      real(wp), intent(out) :: d_rho(:, :), d_grad(:, :, :), d_hess(:, :, :), d_d3(:, :, :)

      real(wp) :: poly(0:4, 3), s(3), beta, exp_r2, dv
      integer :: igrid, idir, iat, a, idx, p, q, r

      d_rho = 0.0_wp
      d_grad = 0.0_wp
      d_hess = 0.0_wp
      d_d3 = 0.0_wp
      do igrid = 1, size(xyz, 2)
         do iat = 1, size(dctx%alpha)
            beta = 2.0_wp*dctx%alpha(iat)
            s = xyz(:, igrid) - dctx%xyz(:, iat)
            exp_r2 = exp(-beta*sum(s*s))
            call gaussian_1d_polys(beta, s, poly)
            do idir = 1, size(d_diag_d, 2)
               dv = d_diag_d(iat, idir)*exp_r2
               d_rho(igrid, idir) = d_rho(igrid, idir) + dv
               do a = 1, 3
                  d_grad(a, igrid, idir) = d_grad(a, igrid, idir) + lobe_deriv(dv, poly, [a])
               end do
               do idx = 1, 6
                  p = sym2_idx(1, idx)
                  q = sym2_idx(2, idx)
                  d_hess(idx, igrid, idir) = d_hess(idx, igrid, idir) &
                     & + lobe_deriv(dv, poly, [p, q])
               end do
               do idx = 1, 10
                  p = sym3_idx(1, idx)
                  q = sym3_idx(2, idx)
                  r = sym3_idx(3, idx)
                  d_d3(idx, igrid, idir) = d_d3(idx, igrid, idir) &
                     & + lobe_deriv(dv, poly, [p, q, r])
               end do
            end do
         end do
      end do
   end subroutine eval_s_density_param_tangent

   pure subroutine contract_s_density_gradient(xyz, dctx, item, grad)
      real(wp), intent(in) :: xyz(:, :)
      type(s_density_ctx_type), intent(in) :: dctx
      type(density_response_type), intent(in) :: item
      real(wp), intent(inout) :: grad(:, :)

      real(wp) :: poly(0:4, 3), s(3), beta, pref, term
      integer :: igrid, iat, k, a, b

      do igrid = 1, size(xyz, 2)
         do iat = 1, size(dctx%alpha)
            beta = 2.0_wp*dctx%alpha(iat)
            s = xyz(:, igrid) - dctx%xyz(:, iat)
            pref = dctx%diag_d(iat)*exp(-beta*sum(s*s))
            call gaussian_1d_polys(beta, s, poly)
            do k = 1, 3
               term = item%w_rho(igrid)*lobe_deriv(pref, poly, [k])
               do a = 1, 3
                  term = term + item%w_grad_rho(a, igrid)*lobe_deriv(pref, poly, [a, k])
                  do b = 1, 3
                     term = term + item%w_hess_rho(a, b, igrid)*lobe_deriv(pref, poly, [a, b, k])
                  end do
               end do
               grad(k, iat) = grad(k, iat) - term
            end do
         end do
      end do
   end subroutine contract_s_density_gradient

   pure subroutine contract_s_density_hvp(xyz, dctx, dirs, item, hvp)
      real(wp), intent(in) :: xyz(:, :), dirs(:, :, :)
      type(s_density_ctx_type), intent(in) :: dctx
      type(density_tangent_response_type), intent(in) :: item
      real(wp), intent(inout) :: hvp(:, :, :)

      real(wp) :: poly(0:4, 3), s(3), u(3), dw2(3, 3), beta, pref, term1, term2, val
      integer :: igrid, idir, iat, k, a, b, m, idx

      do igrid = 1, size(xyz, 2)
         do iat = 1, size(dctx%alpha)
            beta = 2.0_wp*dctx%alpha(iat)
            s = xyz(:, igrid) - dctx%xyz(:, iat)
            pref = dctx%diag_d(iat)*exp(-beta*sum(s*s))
            call gaussian_1d_polys(beta, s, poly)
            do idir = 1, size(dirs, 3)
               do idx = 1, 6
                  a = sym2_idx(1, idx)
                  b = sym2_idx(2, idx)
                  val = item%dw_hess_rho(idx, igrid, idir)
                  dw2(a, b) = val
                  dw2(b, a) = val
               end do
               u = item%d_xyz(:, igrid, idir) - dirs(:, iat, idir)
               do k = 1, 3
                  term1 = item%dw_rho(igrid, idir)*lobe_deriv(pref, poly, [k])
                  do a = 1, 3
                     term1 = term1 + item%dw_grad_rho(a, igrid, idir)*lobe_deriv(pref, poly, [a, k])
                     do b = 1, 3
                        term1 = term1 + dw2(a, b)*lobe_deriv(pref, poly, [a, b, k])
                     end do
                  end do
                  term2 = 0.0_wp
                  do m = 1, 3
                     term2 = term2 + u(m)*(item%w_rho(igrid)*lobe_deriv(pref, poly, [m, k]))
                     do a = 1, 3
                        term2 = term2 + u(m)*(item%w_grad_rho(a, igrid)*lobe_deriv(pref, poly, [m, a, k]))
                        do b = 1, 3
                           term2 = term2 + u(m)*(item%w_hess_rho(a, b, igrid) &
                              & *lobe_deriv(pref, poly, [m, a, b, k]))
                        end do
                     end do
                  end do
                  hvp(k, iat, idir) = hvp(k, iat, idir) - (term1 + term2)
               end do
            end do
         end do
      end do
   end subroutine contract_s_density_hvp

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

end module test_model_continuum_hessian
