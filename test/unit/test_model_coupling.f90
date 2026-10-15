!> Unit tests for coupling declaration and staging
!>
!> Asserts the request set every model configuration declares through
!> `new_coupling` (dedupe, phase masks, mandatory flags), what the
!> `prepare_*` calls leave missing, and how the host walk behaves at model
!> level: staging restarts it, a model update ends it, and every `get_*`
!> names a wrong staging or a missing output
module test_model_coupling
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_model_continuum_component_pcm_type, only: moist_pcm_parameters_type
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_channels_coupling, only: coupling_type, coupling_request_type, &
      gaussian_potential_request_type, gaussian_moment_request_type, request_name_len, &
      nuclear_charge_request_type, density_matrix_request_type, coupling_view_type, &
      coupling_make_view, coupling_close_view
   use moist_cavity_drop_lsf_isodensity_internal, only: moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_drop_lsf_isodensity_param, only: moist_cavity_drop_lsf_isodensity_param_type
   use moist_channels_response, only: response_type, potential_adjoint_response_type, &
      & response_accumulate
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_type, only: cavity_type
   use moist_model_continuum_component_pcm_type, only: solver_type
   use moist_model_continuum_component_pcm_cpcm, only: model_continuum_component_cpcm, new_component_cpcm
   use moist_model_continuum_component, only: model_continuum_component_pv, new_component_pv, &
      & model_continuum_component_gostshyp, new_component_gostshyp
   use moist_model_continuum, only: model_continuum_type, new_continuum_model
   use moist_model_type, only: solvation_model_type
   use moist_cavity_iswig, only: cavity_type_iswig
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: radius_type_static, default_cpcm_radii, new_radii_custom_atoms, radius_type
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: build_test_cavity, fill_point_charge_potential, fill_point_charge_field, &
      & copy_potential_adjoint
   use test_model_component_helper, only: sw_ngrid => fixture_ngrid_param, &
      & sw_areas => fixture_areas_param, sw_xis => fixture_xis_param, &
      & sw_fs => fixture_fs_param, sw_xyz => fixture_xyz_param, &
      & sw_normals => fixture_normals_param

   implicit none(type, external)
   private

   public :: collect_model_coupling

   !> Tolerance for values that must agree to roundoff
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   !> Dielectric constant of the CPCM fixtures
   real(wp), parameter :: eps_model = 32.0_wp
   !> Lebedev order of the iSwiG fixtures
   integer, parameter :: nleb_iswig = 14
   !> Pressure of the PV and GOSTSHYP fixtures
   real(wp), parameter :: test_pressure = 1.67_wp
   !> Point charges filling the potential trace of the MB16-43 fixture
   real(wp), parameter :: qat_vals(*) = [ &
      &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
      &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]

   !> DROP (SvdW) parameters of the O-C-H fixture, as in `cavity_drop_nuclear_adjoint`
   integer, parameter :: nleb_drop = 50
   real(wp), parameter :: drop_blend_k = 2.5_wp
   real(wp), parameter :: drop_blend_3b = 1.0_wp
   real(wp), parameter :: drop_proj_tol = 1.0e-14_wp
   integer, parameter :: drop_proj_maxiter = 1000
   integer, parameter :: drop_proj_level = 2
   integer, parameter :: drop_prune_level = 4

contains
   subroutine collect_model_coupling(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("rejected_calls_preserve_response", test_preserve_response), &
         new_unittest("uninitialized_drop_response", test_uninitialized_drop), &
         new_unittest("multiple_couplings_push_invalidation", test_invalidation), &
         new_unittest("two_distinct_pcm_components_order", test_order), &
         new_unittest("moments_energy_requirements_and_disabled", test_moments), &
         new_unittest("request_free_pv", test_pv), &
         new_unittest("model_copy_owns_independent_couplings", test_model_copy), &
         new_unittest("model_staging_restarts_walk", test_model_restart), &
         new_unittest("phase_answer_reuse_and_energy_refresh", test_phase_reuse), &
         new_unittest("model_update_ends_walk", test_model_update_walk), &
         new_unittest("model_accessors_name_the_problem", test_model_messages), &
         new_unittest("build_answers_survive_energy_cycles", test_build_answers_kept), &
         new_unittest("gradient_charges_go_with_energy_staging", test_gradient_charges_reasked), &
         new_unittest("failed_build_leaves_the_model_unbuilt", test_model_failed_build), &
         new_unittest("component_failure_drops_the_surface", test_component_failure_drops_surface), &
         new_unittest("update_drops_build_answers", test_update_drops_build_answers)]
   end subroutine collect_model_coupling

   subroutine test_invalidation(error)
            type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: model
      type(coupling_type), pointer :: first, second
      !> Borrowed model cavity
      class(cavity_type), pointer :: cav
      real(wp) :: energy
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(model, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      cav => model%cavity
      call model%new_coupling(first, err)
      call model%new_coupling(second, err)
      call model%prepare_energy(first, err)
      call model%prepare_energy(second, err)
      call fill_point_charge_potential(cav, first, qat_vals, mol)
      call fill_point_charge_potential(cav, second, qat_vals, mol)
      energy = 0.0_wp
      call model%get_energy(first, energy, err)
      if (failed(error, err, "energy")) return
      mol%xyz(1, 1) = mol%xyz(1, 1) + 1.0e-4_wp
      call model%update(mol, err)
      if (failed(error, err, "geometry update")) return
      call model%get_energy(first, energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call model%get_energy(second, energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, .not. first%next(), more="an invalidated coupling was walked")
      if (allocated(error)) return
      call model%prepare_response(second, err)
      if (failed(error, err, "response directly after update")) return
      call check(error, second%next(), more="the update kept the potential answer")
      if (allocated(error)) return
      associate (item => second%request())
         call check(error, item%is_missing("phi"))
      end associate
      if (allocated(error)) return
      call check(error, .not. second%next(), more="the potential is the only request")
   end subroutine test_invalidation

   subroutine test_model_copy(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: original
      type(model_continuum_type), allocatable, target :: copied
      class(solvation_model_type), allocatable, target :: polymorphic_copy
      class(solvation_model_type), pointer :: polymorphic_source
      type(coupling_type), pointer :: coupling
      !> Borrowed model cavity
      class(cavity_type), pointer :: cav
      real(wp) :: energy
      integer :: mode
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(original, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      cav => original%cavity
      call original%new_coupling(coupling, err)
      call original%prepare_energy(coupling, err)
      call fill_point_charge_potential(cav, coupling, qat_vals, mol)
      polymorphic_source => original
      do mode = 1, 3
         if (mode == 3) then
            allocate (model_continuum_type :: polymorphic_copy)
            polymorphic_copy = polymorphic_source
            energy = 0.0_wp
            call polymorphic_copy%get_energy(coupling, energy, err)
            call check(error, allocated(err))
            if (allocated(error)) return
            deallocate (polymorphic_copy)
         else
            if (mode == 1) then
               allocate (copied, source=original)
            else
               allocate (copied)
               copied = original
            end if
            energy = 0.0_wp
            call copied%get_energy(coupling, energy, err)
            call check(error, allocated(err))
            if (allocated(error)) return
            deallocate (copied)
         end if
         energy = 0.0_wp
         call original%get_energy(coupling, energy, err)
         if (failed(error, err, "original survives destroying its copy")) return
      end do
      call original%update(mol, err)
      if (failed(error, err, "original registry remains owned")) return
   end subroutine test_model_copy

   subroutine test_order(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: model
      type(model_continuum_component_cpcm) :: pcm1, pcm2
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      !> Borrowed model cavity
      class(cavity_type), pointer :: cav
      real(wp) :: energies(4), repeated_energy
      !> Independent accumulator for repeated calls
      real(wp), allocatable :: repeated_gradient(:, :)
      real(wp), allocatable :: gradients(:, :, :)
      integer :: k
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call new_component_cpcm(pcm1, epsilon=2.0_wp, error=err, ctx=ctx)
      call new_component_cpcm(pcm2, epsilon=30.0_wp, error=err, ctx=ctx)
      allocate(gradients(3, mol%nat, 4), source=0.0_wp)
      allocate(repeated_gradient(3, mol%nat), source=0.0_wp)
      energies = 0.0_wp
      do k = 1, 4
         call new_continuum_model(model, ctx, cavity, err)
         if (k == 1 .or. k == 3) call model%add_component(pcm1, err)
         if (k /= 3) call model%add_component(pcm2, err)
         if (k == 2) call model%add_component(pcm1, err)
         call model%update(mol, err)
         if (failed(error, err, "two PCM update")) return
         cav => model%cavity
         call model%new_coupling(coupling, err)
         call model%prepare_energy(coupling, err)
         call check(error, count_visits(coupling), 1, more="two PCMs must share one potential request")
         if (allocated(error)) return
         call fill_point_charge_potential(cav, coupling, qat_vals, mol)
         call model%get_energy(coupling, energies(k), err)
         if (failed(error, err, "two PCM energy")) return
         repeated_energy = energies(k)
         call model%get_energy(coupling, repeated_energy, err)
         if (failed(error, err, "repeated energy")) return
         call check(error, repeated_energy, 2.0_wp*energies(k), thr=thr)
         if (allocated(error)) return
         call model%prepare_gradient(coupling, err)
         call fill_point_charge_field(cav, coupling, qat_vals, mol)
         call model%get_gradient(coupling, response, gradients(:,:,k), err)
         if (failed(error, err, "two PCM gradient")) return
         repeated_gradient = gradients(:, :, k)
         call model%get_gradient(coupling, response, repeated_gradient, err)
         if (failed(error, err, "repeated gradient")) return
         call check(error, maxval(abs(repeated_gradient - 2.0_wp*gradients(:, :, k))), &
            & 0.0_wp, thr=1.0e-12_wp)
         if (allocated(error)) return
      end do
      call check(error, energies(1), energies(2), thr=thr)
      if (allocated(error)) return
      call check(error, energies(1), energies(3)+energies(4), thr=thr)
      if (allocated(error)) return
      call check(error, maxval(abs(gradients(:,:,1)-gradients(:,:,2))) < 1.0e-12_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(gradients(:,:,1)-gradients(:,:,3)-gradients(:,:,4))) < 1.0e-12_wp)
   end subroutine test_order

   subroutine test_moments(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: model
      type(model_continuum_component_gostshyp) :: component
      type(coupling_type), pointer :: coupling
      !> Borrowed model cavity
      class(cavity_type), pointer :: cav
      real(wp) :: energy
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call new_continuum_model(model, ctx, cavity, err)
      call new_component_gostshyp(component, test_pressure)
      call model%add_component(component, err)
      call model%update(mol, err)
      cav => model%cavity
      call model%new_coupling(coupling, err)
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "moment staging")) return
      call check(error, coupling%next(), more="the energy walk skipped the moments")
      if (allocated(error)) return
      select type (item => coupling%request())
      type is (gaussian_moment_request_type)
         call check(error, item%is_missing("gt") .and. item%is_missing("pt"))
         if (allocated(error)) return
         call check(error, .not. item%is_missing("mt") .and. .not. item%is_missing("rt"), &
            & more="a fixed cavity needs no higher moments for the energy")
         if (allocated(error)) return
         call check(error, size(item%width), cav%ngrid)
         if (allocated(error)) return
         call check(error, all(item%width >= 0.0_wp))
         if (allocated(error)) return
         !> Narrow points are switched off at zero width; they are the smallest ones
         call check(error, any(item%width > 0.0_wp))
         if (allocated(error)) return
         call check(error, maxval(abs(item%width*cav%a - pi*log(2.0_wp)), &
            & mask=item%width > 0.0_wp), 0.0_wp, thr=thr)
         if (allocated(error)) return
         call check(error, maxval(cav%a, mask=item%width == 0.0_wp) &
            & < minval(cav%a, mask=item%width > 0.0_wp))
      class default
         call test_failed(error, "unexpected request "//item%name())
      end select
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="GOSTSHYP declares one request")
      if (allocated(error)) return
      call model%prepare_response(coupling, err)
      if (failed(error, err, "moment response staging")) return
      call check(error, coupling%next())
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("gt") .and. item%is_missing("pt"))
         if (allocated(error)) return
         call check(error, .not. item%is_missing("mt") .and. .not. item%is_missing("rt"))
      end associate
      if (allocated(error)) return
      call check(error, .not. coupling%next())
      if (allocated(error)) return
      call model%prepare_gradient(coupling, err)
      if (failed(error, err, "moment gradient staging")) return
      call check(error, coupling%next())
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("gt") .and. item%is_missing("pt") .and. &
            & item%is_missing("mt") .and. item%is_missing("rt"))
      end associate
      if (allocated(error)) return
      call check(error, .not. coupling%next())
      if (allocated(error)) return
      call new_continuum_model(model, ctx, cavity, err)
      call new_component_gostshyp(component, 0.0_wp)
      call model%add_component(component, err)
      call model%update(mol, err)
      call model%new_coupling(coupling, err)
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "disabled moment staging")) return
      call check(error, count_visits(coupling), 0, more="zero pressure still asked for moments")
      if (allocated(error)) return
      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      if (failed(error, err, "disabled GOSTSHYP energy")) return
      call check(error, energy, 0.0_wp, thr=0.0_wp)
      if (allocated(error)) return
      call new_continuum_model(model, ctx, cavity, err)
      call new_component_gostshyp(component, test_pressure)
      component%scale = 0.0_wp
      call model%add_component(component, err)
      call model%update(mol, err)
      call model%new_coupling(coupling, err)
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "zero scale moment staging")) return
      call check(error, count_visits(coupling), 0, more="zero scale still asked for moments")
      if (allocated(error)) return
      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      if (failed(error, err, "zero scale GOSTSHYP energy")) return
      call check(error, energy, 0.0_wp, thr=0.0_wp)
   end subroutine test_moments

   subroutine test_pv(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: model
      type(model_continuum_component_pv) :: component
      type(coupling_type), pointer :: coupling
      real(wp) :: energy
      !> Pressure times cavity volume, the expected PV energy
      real(wp) :: pv
      !> Absolute tolerance on p*V, with a 10x looser relative fallback
      real(wp), parameter :: pv_abs_thr = 1.0e-10_wp, pv_rel_thr = 1.0e-9_wp
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call new_continuum_model(model, ctx, cavity, err)
      call new_component_pv(component, test_pressure)
      call model%add_component(component, err)
      call model%update(mol, err)
      call model%new_coupling(coupling, err)
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "PV staging")) return
      call check(error, .not. coupling%next(), more="PV asked the host for something")
      if (allocated(error)) return
      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      if (failed(error, err, "PV energy without any answer")) return
      call check(error, energy > 0.0_wp, more="PV energy must be positive at positive pressure")
      if (allocated(error)) return
      ! Measured |E - p*V| <= 4e-14 at E ~ 2.0e3 (below one ulp), macOS arm64
      ! gfortran 14.3, 1-4 threads
      pv = test_pressure*model%cavity%total_volume
      call check(error, abs(energy - pv) <= pv_abs_thr .or. &
         & abs(energy - pv) <= pv_rel_thr*abs(pv), more="PV energy differs from p*V")
   end subroutine test_pv

   logical function failed(error, err, context)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> moist error
      type(moist_error_type), allocatable, intent(inout) :: err
      !> What was attempted
      character(len=*), intent(in) :: context

      failed = allocated(err)
      if (failed) then
         call test_failed(error, context//" failed: "//err%message)
         err%stat = 0
         deallocate (err)
      end if

   end function failed

   !> MB16-43 "01" on the iSwiG fixture cavity of `test_model_continuum`
   !>
   !> @param[out] ctx          Run context owned by the caller
   !> @param[out] mol          Structure
   !> @param[out] radius_model Radius model storage
   !> @param[out] cavity       Updated cavity
   !> @param[out] error        testdrive failure
   subroutine iswig_fixture(ctx, mol, radius_model, cavity, error)
      !> Run context owned by the caller
      type(moist_context_type), intent(out), target :: ctx
      !> Structure
      type(structure_type), intent(out) :: mol
      !> Radius model storage
      type(radius_type_static), intent(out) :: radius_model
      !> Updated cavity
      type(cavity_type_iswig), intent(out) :: cavity
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err

      call new_context(ctx, nthreads=0)
      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, nleb_iswig, ctx, radius_model, cavity, err)
      if (failed(error, err, "iSwiG cavity setup")) return

   end subroutine iswig_fixture

   !> DROP (SvdW) cavity on the O-C-H fixture of `test_cavity_drop_nuclear_adjoint`
   !>
   !> @param[out] ctx    Run context owned by the caller
   !> @param[out] mol    Structure
   !> @param[out] cavity Configured cavity (updated by the model)
   !> @param[out] error  testdrive failure
   subroutine drop_fixture(ctx, mol, cavity, error)
      !> Run context owned by the caller
      type(moist_context_type), intent(out), target :: ctx
      !> Structure
      type(structure_type), intent(out) :: mol
      !> Configured cavity
      type(cavity_type_drop), intent(out) :: cavity
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err

      call new(mol, [8, 6, 1], reshape([ &
                                       0.00_wp, 0.00_wp, 0.00_wp, &
                                       0.00_wp, 0.00_wp, 4.60_wp, &
                                       2.60_wp, 0.40_wp, -1.10_wp], [3, 3]))
      call new_context(ctx, nthreads=0, verbosity=0)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=drop_blend_k, &
            blend_3b=drop_blend_3b))
         call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
            error=err, param=moist_cavity_drop_parameters_type(num_leb=nleb_drop, tolerance=drop_proj_tol, &
            proj_maxiter=drop_proj_maxiter, proj_level=drop_proj_level, wleb_prune_level=drop_prune_level, &
            do_fine=.true.), ctx=ctx)
      end block
      if (failed(error, err, "DROP setup")) return

   end subroutine drop_fixture

   !> Updated continuum model with a CPCM component and optionally a PV component
   !>
   !> @param[out] model          Model to build
   !> @param[in]  cavity         Cavity template copied into the model
   !> @param[in]  ctx            Run context owned by the caller
   !> @param[in]  mol            Structure
   !> @param[in]  with_pv        Whether to append a PV component
   !> @param[out] error          testdrive failure
   !> @param[in]  configure_only Leave the model un-updated, so that more
   !>                            components can still be added
   subroutine cpcm_model(model, cavity, ctx, mol, with_pv, error, configure_only)
      !> Model to build
      type(model_continuum_type), intent(out) :: model
      !> Cavity template copied into the model
      class(cavity_type), intent(in) :: cavity
      !> Run context owned by the caller
      type(moist_context_type), intent(in), target :: ctx
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Whether to append a PV component
      logical, intent(in) :: with_pv
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> Leave the model un-updated
      logical, intent(in), optional :: configure_only

      type(moist_error_type), allocatable :: err
      type(model_continuum_component_cpcm) :: cpcm
      type(model_continuum_component_pv) :: pv

      call new_continuum_model(model, ctx, cavity, err)
      if (failed(error, err, "model setup")) return
      if (with_pv) then
         call new_component_pv(pv, test_pressure)
         call model%add_component(pv, err)
         if (failed(error, err, "adding PV")) return
      end if
      call new_component_cpcm(cpcm, epsilon=eps_model, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%cholesky), ctx=ctx)
      if (failed(error, err, "CPCM setup")) return
      call model%add_component(cpcm, err)
      if (failed(error, err, "adding CPCM")) return
      if (present(configure_only)) then
         if (configure_only) return
      end if
      call model%update(mol, err)
      if (failed(error, err, "model update")) return

   end subroutine cpcm_model

   !> Invalid model calls leave the caller's previous response intact
   !>
   !> The foreign coupling is staged and answered by its own model, so only
   !> the ownership guard refuses it; its owner accepts it as a control
   subroutine test_preserve_response(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Shared fixture state
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: model, other
      type(coupling_type), pointer :: coupling, foreign
      !> Cavity of the foreign coupling's owner
      class(cavity_type), pointer :: other_cav
      !> Previous response and sentinel charge
      type(response_type) :: response
      type(potential_adjoint_response_type) :: charge
      type(potential_adjoint_response_type), allocatable :: saved
      !> Gradient accumulator, unchanged on failure
      real(wp), allocatable :: gradient(:, :)
      !> Response and gradient of the owner's control call
      type(response_type) :: control_response
      real(wp), allocatable :: control_gradient(:, :)
      !> Rejection case
      integer :: scenario
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(model, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      call cpcm_model(other, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      call model%new_coupling(coupling, err)
      call other%new_coupling(foreign, err)
      if (failed(error, err, "coupling setup")) return
      other_cav => other%cavity
      charge%w_phi = [7.0_wp]
      call response_accumulate(response, charge, err)
      allocate (gradient(3, mol%nat), source=9.0_wp)
      allocate (control_gradient(3, mol%nat), source=0.0_wp)
      do scenario = 1, 5
         select case (scenario)
         case (1)
            call other%prepare_response(foreign, err)
            if (failed(error, err, "foreign response staging")) return
            call fill_point_charge_potential(other_cav, foreign, qat_vals, mol)
            call model%get_response(foreign, response, err)
         case (2)
            call model%get_gradient(coupling, response, gradient, err)
         case (3)
            call model%prepare_response(coupling, err)
            call model%get_response(coupling, response, err)
         case (4)
            call model%update(mol, err)
            call model%get_gradient(coupling, response, gradient, err)
         case (5)
            call other%prepare_gradient(foreign, err)
            if (failed(error, err, "foreign gradient staging")) return
            call fill_point_charge_potential(other_cav, foreign, qat_vals, mol)
            call fill_point_charge_field(other_cav, foreign, qat_vals, mol)
            call model%get_gradient(foreign, response, gradient, err)
         case default
            error stop "test_model_coupling: unhandled scenario"
         end select
         call check(error, allocated(err))
         if (allocated(error)) return
         if (scenario == 1 .or. scenario == 5) then
            call check(error, index(err%message, "different model") > 0, more=err%message)
            if (allocated(error)) return
         end if
         call copy_potential_adjoint(response, saved)
         call check(error, allocated(saved))
         if (allocated(error)) return
         call check(error, saved%w_phi(1), 7.0_wp, thr=0.0_wp)
         if (allocated(error)) return
         call check(error, maxval(abs(gradient - 9.0_wp)), 0.0_wp, thr=0.0_wp)
         if (allocated(error)) return
         ! The owner accepts the same staging, so ownership alone refused it
         if (scenario == 1) then
            call other%get_response(foreign, control_response, err)
            if (failed(error, err, "owner response on its staged coupling")) return
         else if (scenario == 5) then
            call other%get_gradient(foreign, control_response, control_gradient, err)
            if (failed(error, err, "owner gradient on its staged coupling")) return
         end if
      end do
      call model%release_coupling(coupling)
      call other%release_coupling(foreign)
   end subroutine test_preserve_response

   !> A bare DROP cavity reports its missing level set without dereferencing it
   subroutine test_uninitialized_drop(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Uninitialized cavity and empty accumulators
      type(cavity_type_drop) :: cavity
      type(cavity_surface_adjoint_type) :: weights
      type(response_type) :: response
      call cavity%get_surface_response(weights, response, err)
      call check(error, allocated(err))
   end subroutine test_uninitialized_drop

   !> Number of requests one full pass of the walk visits, answering none
   function count_visits(coupling) result(visits)
      !> Staged coupling
      type(coupling_type), intent(inout) :: coupling
      !> Visited requests
      integer :: visits
      visits = 0
      do while (coupling%next())
         visits = visits + 1
      end do
   end function count_visits

   !> Canonical name of the current request, "no_current_request" outside a `next()` window
   function current_name(coupling) result(name)
      !> Coupling to inspect
      type(coupling_type), intent(in) :: coupling
      !> Name of the current request
      character(len=request_name_len) :: name
      class(coupling_request_type), allocatable :: item
      item = coupling%request()
      name = item%name()
   end function current_name

   !> A new `prepare_*` restarts a half-walked pass at the first request
   subroutine test_model_restart(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Fixture: a CPCM and a GOSTSHYP component, two requests
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: model
      type(model_continuum_component_gostshyp) :: gostshyp
      type(coupling_type), pointer :: coupling
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(model, cavity, ctx, mol, .false., error, configure_only=.true.)
      if (allocated(error)) return
      call new_component_gostshyp(gostshyp, test_pressure)
      call model%add_component(gostshyp, err)
      if (failed(error, err, "adding GOSTSHYP")) return
      call model%update(mol, err)
      if (failed(error, err, "model update")) return
      call model%new_coupling(coupling, err)
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "energy staging")) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "gaussian_potential")
      if (allocated(error)) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "gaussian_moments")
      if (allocated(error)) return
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "energy restaging")) return
      call check(error, current_name(coupling), "no_current_request", more="staging kept a current request")
      if (allocated(error)) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "gaussian_potential", &
         & more="prepare_energy did not restart the walk")
      if (allocated(error)) return
      ! The same for the later phases, which keep their answers
      call model%prepare_gradient(coupling, err)
      if (failed(error, err, "gradient staging")) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call check(error, current_name(coupling), "gaussian_potential")
      if (allocated(error)) return
      call model%release_coupling(coupling)
   end subroutine test_model_restart

   !> Later phases reuse phi; a new energy evaluation discards it
   subroutine test_phase_reuse(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Fixture
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: model
      type(coupling_type), pointer :: coupling
      class(cavity_type), pointer :: cav
      !> Host response
      type(response_type) :: response
      type(potential_adjoint_response_type), allocatable :: charge
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(model, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      cav => model%cavity
      call model%new_coupling(coupling, err)
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "energy staging")) return
      call fill_point_charge_potential(cav, coupling, qat_vals, mol)
      call model%prepare_response(coupling, err)
      if (failed(error, err, "response staging")) return
      call check(error, count_visits(coupling), 0, more="response discarded valid phi")
      if (allocated(error)) return
      call model%get_response(coupling, response, err)
      if (failed(error, err, "response with reused phi")) return
      call copy_potential_adjoint(response, charge)
      call check(error, allocated(charge))
      if (allocated(error)) return
      call check(error, maxval(abs(charge%w_phi)) > 0.0_wp)
      if (allocated(error)) return
      call model%prepare_gradient(coupling, err)
      if (failed(error, err, "gradient staging")) return
      call check(error, coupling%next())
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, .not. item%is_missing("phi"))
         if (allocated(error)) return
         call check(error, item%is_missing("dphi_dr") .and. item%is_missing("dphi_dxi"))
      end associate
      if (allocated(error)) return
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "new energy staging")) return
      call check(error, coupling%next(), more="new energy retained the previous phi")
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("phi"))
      end associate
      if (allocated(error)) return
      call model%prepare_gradient(coupling, err)
      if (failed(error, err, "gradient without energy answers")) return
      call check(error, coupling%next())
      if (allocated(error)) return
      associate (item => coupling%request())
         call check(error, item%is_missing("phi") .and. item%is_missing("dphi_dr") .and. &
            & item%is_missing("dphi_dxi"))
      end associate
      if (allocated(error)) return
      call model%release_coupling(coupling)
      call check(error, .not. associated(coupling))
   end subroutine test_phase_reuse

   !> A model update ends every walk; answering is refused until the next staging
   subroutine test_model_update_walk(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Fixture
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: model
      type(coupling_type), pointer :: coupling
      !> Borrowed model cavity
      class(cavity_type), pointer :: cav
      !> Potential answer on the current grid
      real(wp), allocatable :: phi(:)
      real(wp) :: energy
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(model, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      cav => model%cavity
      call model%new_coupling(coupling, err)
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "energy staging")) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call model%update(mol, err)
      if (failed(error, err, "model update")) return
      call check(error, current_name(coupling), "no_current_request", more="the update kept a current request")
      if (allocated(error)) return
      allocate (phi(cav%ngrid), source=0.0_wp)
      call coupling%answer("phi", phi, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "No current coupling request - call next() first")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="an invalidated coupling was walked")
      if (allocated(error)) return
      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "this coupling is not staged") > 0)
      if (allocated(error)) return
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "energy restaging")) return
      call check(error, count_visits(coupling), 1)
      if (allocated(error)) return
      call model%release_coupling(coupling)
   end subroutine test_model_update_walk

   !> Every accessor names a wrong staging and the outputs still missing
   subroutine test_model_messages(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      !> Fixture
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: model
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      !> Borrowed model cavity
      class(cavity_type), pointer :: cav
      real(wp), allocatable :: gradient(:, :)
      real(wp) :: energy
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(model, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      cav => model%cavity
      call model%new_coupling(coupling, err)
      if (failed(error, err, "coupling setup")) return
      allocate (gradient(3, mol%nat), source=0.0_wp)
      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "get_energy requires a coupling staged by prepare_energy; "// &
         & "this coupling is not staged")
      if (allocated(error)) return
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "energy staging")) return
      call model%get_gradient(coupling, response, gradient, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "get_gradient requires a coupling staged by prepare_gradient; "// &
         & "this coupling is staged for the energy phase")
      if (allocated(error)) return
      call model%get_response(coupling, response, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "get_response requires a coupling staged by prepare_response; "// &
         & "this coupling is staged for the energy phase")
      if (allocated(error)) return
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "gaussian_potential: missing required outputs: phi")
      if (allocated(error)) return
      call fill_point_charge_potential(cav, coupling, qat_vals, mol)
      call model%get_energy(coupling, energy, err)
      if (failed(error, err, "energy")) return
      call model%prepare_gradient(coupling, err)
      if (failed(error, err, "gradient staging")) return
      call model%get_gradient(coupling, response, gradient, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, err%message, "gaussian_potential: missing required outputs: dphi_dr, dphi_dxi")
      if (allocated(error)) return
      call check(error, maxval(abs(gradient)), 0.0_wp, thr=0.0_wp, &
         & more="a refused gradient touched the accumulator")
      if (allocated(error)) return
      call model%release_coupling(coupling)
   end subroutine test_model_messages

   !> Count the nuclear_charges requests a pass visits, answering nothing
   !>
   !> @param[inout] coupling Staged coupling; the pass ends and rewinds
   function count_charge_visits(coupling) result(visits)
      !> Staged coupling
      type(coupling_type), intent(inout) :: coupling
      !> Visited nuclear_charges requests
      integer :: visits
      visits = 0
      do while (coupling%next())
         select type (request => coupling%request())
         type is (nuclear_charge_request_type)
            visits = visits + 1
         end select
      end do
   end function count_charge_visits

   !> The build phase answers the density matrix and the nuclear charges once:
   !> two energy cycles later the gradient reads `z` without asking for it
   !>
   !> Internal isodensity cavity around a pseudo atom (one d shell) with the
   !> core fill on, so `z` is required in the build phase; `z = Z` marks no
   !> ECP nucleus, the surface stays that of the bare density. Before the
   !> build, the energy staging names the build phase and the missing output
   subroutine test_build_answers_kept(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(cavity_type_drop) :: cavity
      class(radius_type), allocatable :: radii
      type(model_continuum_type), target :: model
      type(coupling_type), pointer :: coupling
      class(cavity_type), pointer :: cav
      type(response_type) :: response
      real(wp) :: dcart(6, 6), energy
      real(wp), allocatable :: gradient(:, :)
      integer :: i, icycle, nvisit

      call new_context(ctx, nthreads=0, verbosity=0)
      call new(mol, [54], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      dcart = 0.0_wp
      do i = 1, 6
         dcart(i, i) = 1.0_wp
      end do
      call lsf%new(sh_atom=[1], sh_l=[2], sh_nprim=[1], exps=[0.5_wp], coeffs=[1.0_wp], error=err, &
                   param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=0.4_wp, scale=2.5_wp))
      if (failed(error, err, "LSF setup")) return
      call new_radii_custom_atoms([4.0_wp], radii, err)
      if (failed(error, err, "radii setup")) return
      call new_cavity_drop(cavity, radius_model=radii, lsf_model=lsf, error=err, &
                           param=moist_cavity_drop_parameters_type(num_leb=6, proj_level=1), ctx=ctx)
      if (failed(error, err, "cavity setup")) return
      call cpcm_model(model, cavity, ctx, mol, .false., error, configure_only=.true.)
      if (allocated(error)) return

      ! The update leaves the cavity without a surface until the build
      call model%update(mol, err)
      if (failed(error, err, "model update")) return
      cav => model%cavity
      call check(error, model%is_updated() .and. .not. model%is_built(), &
         & more="the model built without its density matrix")
      if (allocated(error)) return
      call model%new_coupling(coupling, err)
      if (failed(error, err, "coupling of the updated model")) return
      call check_model_energy_refused(model, coupling, error)
      if (allocated(error)) return

      call model%prepare_build(coupling, err)
      if (failed(error, err, "build staging")) return
      nvisit = 0
      do while (coupling%next())
         nvisit = nvisit + 1
         select type (request => coupling%request())
         type is (density_matrix_request_type)
            call coupling%answer("dcart", dcart, err)
         type is (nuclear_charge_request_type)
            call coupling%answer("z", [54.0_wp], err)
         class default
            call test_failed(error, "the build phase asked for "//trim(request%name()))
            return
         end select
         if (failed(error, err, "build answer")) return
      end do
      call check(error, nvisit, 2, more="the build phase asks for the density matrix and the charges")
      if (allocated(error)) return
      call model%build(coupling, err)
      if (failed(error, err, "model build")) return
      call check(error, model%is_built() .and. cav%ngrid > 0, more="the pseudo atom has no surface")
      if (allocated(error)) return

      do icycle = 1, 2
         call model%prepare_energy(coupling, err)
         if (failed(error, err, "energy staging")) return
         call check(error, count_charge_visits(coupling), 0, more="the energy phase asked for z")
         if (allocated(error)) return
         call fill_point_charge_potential(cav, coupling, [0.3_wp], mol)
         energy = 0.0_wp
         call model%get_energy(coupling, energy, err)
         if (failed(error, err, "energy")) return
      end do

      call model%prepare_gradient(coupling, err)
      if (failed(error, err, "gradient staging")) return
      call check(error, count_charge_visits(coupling), 0, &
         & more="the gradient asked again for z answered in the build phase")
      if (allocated(error)) return
      call fill_point_charge_field(cav, coupling, [0.3_wp], mol)
      allocate (gradient(3, mol%nat), source=0.0_wp)
      call model%get_gradient(coupling, response, gradient, err)
      if (failed(error, err, "gradient with build-phase charges")) return
      call model%release_coupling(coupling)
   end subroutine test_build_answers_kept

   !> Negative control: on a cavity without build-phase requests, `z` comes
   !> in the gradient phase and goes with the next energy staging
   subroutine test_gradient_charges_reasked(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(model_continuum_type), target :: model
      type(coupling_type), pointer :: coupling
      class(cavity_type), pointer :: cav
      type(response_type) :: response
      real(wp), allocatable :: gradient(:, :)
      real(wp) :: energy
      integer :: icycle

      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(model, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      cav => model%cavity
      call model%new_coupling(coupling, err)
      if (failed(error, err, "coupling")) return
      allocate (gradient(3, mol%nat), source=0.0_wp)
      do icycle = 1, 2
         call model%prepare_energy(coupling, err)
         if (failed(error, err, "energy staging")) return
         call fill_point_charge_potential(cav, coupling, qat_vals, mol)
         energy = 0.0_wp
         call model%get_energy(coupling, energy, err)
         if (failed(error, err, "energy")) return
         call model%prepare_gradient(coupling, err)
         if (failed(error, err, "gradient staging")) return
         call check(error, count_charge_visits(coupling), 1, more="the gradient phase did not ask for z")
         if (allocated(error)) return
         call fill_point_charge_field(cav, coupling, qat_vals, mol)
         call model%get_gradient(coupling, response, gradient, err)
         if (failed(error, err, "gradient")) return
      end do
      call model%release_coupling(coupling)
   end subroutine test_gradient_charges_reasked

   !> Stage the energy phase of an updated model that is not built
   !>
   !> The staging is refused by name
   !>
   !> @param[inout] model    Updated model, not built
   !> @param[inout] coupling Coupling of the model
   !> @param[out]   error    testdrive failure
   subroutine check_model_energy_refused(model, coupling, error)
      !> Updated model, not built
      type(model_continuum_type), intent(inout) :: model
      !> Coupling of the model
      type(coupling_type), intent(inout), target :: coupling
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err

      call check(error, model%is_updated() .and. .not. model%is_built(), &
         & more="the model is not updated or already built")
      if (allocated(error)) return
      call model%prepare_energy(coupling, err)
      call check(error, allocated(err), more="an unbuilt model staged its energy")
      if (allocated(error)) return
      call check(error, err%message, "Solvation model must be built first")
      if (allocated(error)) return
      err%stat = 0
      deallocate (err)

   end subroutine check_model_energy_refused

   !> CPCM model over an internal isodensity pseudo atom, not updated yet
   !>
   !> One d shell with the core fill on, so the build phase asks for the
   !> density matrix and the nuclear charges
   !>
   !> @param[out] model Model to build
   !> @param[in]  ctx   Run context owned by the caller
   !> @param[out] mol   Single atom at the origin
   !> @param[out] error testdrive failure
   subroutine pseudo_atom_model(model, ctx, mol, error)
      !> Model to build
      type(model_continuum_type), intent(out) :: model
      !> Run context owned by the caller
      type(moist_context_type), intent(in), target :: ctx
      !> Single atom at the origin
      type(structure_type), intent(out) :: mol
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(cavity_type_drop) :: cavity
      class(radius_type), allocatable :: radii

      call new(mol, [54], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call lsf%new(sh_atom=[1], sh_l=[2], sh_nprim=[1], exps=[0.5_wp], coeffs=[1.0_wp], error=err, &
                   param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=0.4_wp, scale=2.5_wp))
      if (failed(error, err, "LSF setup")) return
      call new_radii_custom_atoms([4.0_wp], radii, err)
      if (failed(error, err, "radii setup")) return
      call new_cavity_drop(cavity, radius_model=radii, lsf_model=lsf, error=err, &
                           param=moist_cavity_drop_parameters_type(num_leb=6, proj_level=1), ctx=ctx)
      if (failed(error, err, "cavity setup")) return
      call cpcm_model(model, cavity, ctx, mol, .false., error, configure_only=.true.)

   end subroutine pseudo_atom_model

   !> Answer the build phase of the pseudo atom: every cartesian d component
   !> once, `z = Z`
   !>
   !> @param[inout] coupling Coupling staged for the build phase
   !> @param[out]   err      Refused answer
   subroutine answer_pseudo_atom_build(coupling, err)
      !> Coupling staged for the build phase
      type(coupling_type), intent(inout) :: coupling
      !> Refused answer
      type(moist_error_type), allocatable, intent(out) :: err

      real(wp) :: dcart(6, 6)
      integer :: i

      dcart = 0.0_wp
      do i = 1, 6
         dcart(i, i) = 1.0_wp
      end do
      do while (coupling%next())
         select type (request => coupling%request())
         type is (density_matrix_request_type)
            call coupling%answer("dcart", dcart, err)
         type is (nuclear_charge_request_type)
            call coupling%answer("z", [54.0_wp], err)
         end select
         if (allocated(err)) return
      end do

   end subroutine answer_pseudo_atom_build

   !> A model whose cavity defers its surface is not built from `prepare_build`
   !> until a `build` succeeds
   !>
   !> - staging the build phase unbuilds it, so skipping the build leaves the
   !>   energy refused
   !> - a build refused for a missing output leaves it unbuilt
   !> - so does a build of a coupling staged for another phase
   !> - control: a geometric model stays built through `prepare_build`, and
   !>   after an answered build the energy is evaluated
   subroutine test_model_failed_build(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx, ctx_iswig
      type(structure_type) :: mol, mol_iswig
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: iswig
      type(model_continuum_type), target :: model, geometric
      type(coupling_type), pointer :: coupling
      class(cavity_type), pointer :: cav
      real(wp) :: energy

      ! Control: the geometric model keeps its surface
      call iswig_fixture(ctx_iswig, mol_iswig, radii, iswig, error)
      if (allocated(error)) return
      call cpcm_model(geometric, iswig, ctx_iswig, mol_iswig, .false., error)
      if (allocated(error)) return
      call geometric%new_coupling(coupling, err)
      if (failed(error, err, "geometric coupling")) return
      call geometric%prepare_build(coupling, err)
      if (failed(error, err, "geometric build staging")) return
      cav => geometric%cavity
      call check(error, geometric%is_built() .and. cav%is_built(), more="prepare_build dropped a geometric surface")
      if (allocated(error)) return
      call geometric%prepare_energy(coupling, err)
      if (failed(error, err, "geometric energy staging without a build")) return
      call fill_point_charge_potential(cav, coupling, qat_vals, mol_iswig)
      energy = 0.0_wp
      call geometric%get_energy(coupling, energy, err)
      if (failed(error, err, "geometric energy without a build")) return
      call geometric%release_coupling(coupling)

      call new_context(ctx, nthreads=0, verbosity=0)
      call pseudo_atom_model(model, ctx, mol, error)
      if (allocated(error)) return
      call model%update(mol, err)
      if (failed(error, err, "model update")) return
      cav => model%cavity
      call model%new_coupling(coupling, err)
      if (failed(error, err, "coupling")) return
      call build_and_evaluate(error)
      if (allocated(error)) return

      ! Staged but not built
      call model%prepare_build(coupling, err)
      if (failed(error, err, "build staging")) return
      call check(error, .not. model%is_built(), more="prepare_build kept the previous surface")
      if (allocated(error)) return
      call check_model_energy_refused(model, coupling, error)
      if (allocated(error)) return

      ! Built without the answers
      call model%prepare_build(coupling, err)
      if (.not. allocated(err)) call model%build(coupling, err)
      call check(error, allocated(err), more="a build without its answers succeeded")
      if (allocated(error)) return
      call check(error, index(err%message, "density_matrix: missing required outputs: dcart") > 0, &
         & more=err%message)
      if (allocated(error)) return
      err%stat = 0
      deallocate (err)
      call check(error, .not. model%is_built(), more="a failed build left a surface")
      if (allocated(error)) return
      call check_model_energy_refused(model, coupling, error)
      if (allocated(error)) return

      ! Built again, then a build of the coupling staged for the energy phase
      call build_and_evaluate(error)
      if (allocated(error)) return
      call model%build(coupling, err)
      call check(error, allocated(err), more="a build of an energy staging succeeded")
      if (allocated(error)) return
      call check(error, index(err%message, "staged for the energy phase") > 0, more=err%message)
      if (allocated(error)) return
      err%stat = 0
      deallocate (err)
      call check(error, .not. model%is_built(), more="a refused build left a surface")
      if (allocated(error)) return
      call check_model_energy_refused(model, coupling, error)
      if (allocated(error)) return
      call model%release_coupling(coupling)

   contains

      !> Run the answered build phase, then evaluate the energy
      subroutine build_and_evaluate(error)
         !> Test failure
         type(error_type), allocatable, intent(out) :: error

         call model%prepare_build(coupling, err)
         if (.not. allocated(err)) call answer_pseudo_atom_build(coupling, err)
         if (.not. allocated(err)) call model%build(coupling, err)
         if (.not. allocated(err)) call model%prepare_energy(coupling, err)
         if (failed(error, err, "control: build and energy staging")) return
         call check(error, model%is_built(), more="control: the build left no surface")
         if (allocated(error)) return
         call fill_point_charge_potential(cav, coupling, [0.3_wp], mol)
         energy = 0.0_wp
         call model%get_energy(coupling, energy, err)
         if (failed(error, err, "control: energy after the build")) return
      end subroutine build_and_evaluate

   end subroutine test_model_failed_build

   !> A component that fails to update in the build leaves the model unbuilt
   !>
   !> The CPCM component carries a 1x1 external matrix, the pseudo atom's
   !> surface more points, so the component update in `build` fails after the
   !> cavity projected; the model stays unbuilt and the energy is refused.
   !> Negative control: the same model without an external matrix builds and
   !> evaluates
   subroutine test_component_failure_drops_surface(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(cavity_type_drop) :: cavity
      class(radius_type), allocatable :: radii
      type(model_continuum_component_cpcm) :: cpcm
      type(model_continuum_type), target :: model, control
      type(coupling_type), pointer :: coupling
      class(cavity_type), pointer :: cav
      real(wp) :: energy

      call new_context(ctx, nthreads=0, verbosity=0)

      ! Control: no external matrix
      call pseudo_atom_model(control, ctx, mol, error)
      if (allocated(error)) return
      call control%update(mol, err)
      if (.not. allocated(err)) call control%new_coupling(coupling, err)
      if (.not. allocated(err)) call control%prepare_build(coupling, err)
      if (.not. allocated(err)) call answer_pseudo_atom_build(coupling, err)
      if (.not. allocated(err)) call control%build(coupling, err)
      if (.not. allocated(err)) call control%prepare_energy(coupling, err)
      if (failed(error, err, "control: build and energy staging")) return
      cav => control%cavity
      call fill_point_charge_potential(cav, coupling, [0.3_wp], mol)
      energy = 0.0_wp
      call control%get_energy(coupling, energy, err)
      if (failed(error, err, "control: energy")) return
      call control%release_coupling(coupling)

      call new(mol, [54], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call lsf%new(sh_atom=[1], sh_l=[2], sh_nprim=[1], exps=[0.5_wp], coeffs=[1.0_wp], error=err, &
                   param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=0.4_wp, scale=2.5_wp))
      if (failed(error, err, "LSF setup")) return
      call new_radii_custom_atoms([4.0_wp], radii, err)
      if (failed(error, err, "radii setup")) return
      call new_cavity_drop(cavity, radius_model=radii, lsf_model=lsf, error=err, &
                           param=moist_cavity_drop_parameters_type(num_leb=6, proj_level=1), ctx=ctx)
      if (failed(error, err, "cavity setup")) return
      call new_continuum_model(model, ctx, cavity, err)
      if (failed(error, err, "model setup")) return
      call new_component_cpcm(cpcm, epsilon=eps_model, external_matrix=reshape([1.0_wp], [1, 1]), &
         error=err, param=moist_pcm_parameters_type(solver=solver_type%cholesky), ctx=ctx)
      if (.not. allocated(err)) call model%add_component(cpcm, err)
      if (failed(error, err, "CPCM with an external matrix")) return

      call model%update(mol, err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_build(coupling, err)
      if (.not. allocated(err)) call answer_pseudo_atom_build(coupling, err)
      if (failed(error, err, "build staging")) return
      call model%build(coupling, err)
      call check(error, allocated(err), more="a build with a mis-sized external matrix succeeded")
      if (allocated(error)) return
      call check(error, index(err%message, "External PCM matrix dimension mismatch") > 0, more=err%message)
      if (allocated(error)) return
      err%stat = 0
      deallocate (err)
      cav => model%cavity
      call check(error, .not. model%is_built(), more="a failed component update left the model built")
      if (allocated(error)) return
      call check_model_energy_refused(model, coupling, error)
      if (allocated(error)) return
      call model%release_coupling(coupling)
   end subroutine test_component_failure_drops_surface

   !> A model update drops the answers of the build phase
   !>
   !> The density matrix the cavity scope reads after the build is gone after
   !> the next update, the walk ends, the energy is refused and the build
   !> phase asks for both outputs again
   subroutine test_update_drops_build_answers(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(model_continuum_type), target :: model
      type(coupling_type), pointer :: coupling
      type(coupling_view_type) :: view
      real(wp), allocatable :: dcart(:, :)
      integer :: nvisit

      call new_context(ctx, nthreads=0, verbosity=0)
      call pseudo_atom_model(model, ctx, mol, error)
      if (allocated(error)) return
      call model%update(mol, err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (failed(error, err, "update and coupling")) return
      call model%prepare_build(coupling, err)
      if (.not. allocated(err)) call answer_pseudo_atom_build(coupling, err)
      if (.not. allocated(err)) call model%build(coupling, err)
      if (failed(error, err, "build")) return
      call coupling_make_view(coupling, 0, view)
      call view%read("density_matrix", "dcart", dcart, err)
      call coupling_close_view(coupling)
      if (failed(error, err, "control: density read after the build")) return

      call model%update(mol, err)
      if (failed(error, err, "second update")) return
      call check(error, .not. coupling%next(), more="the update left the walk running")
      if (allocated(error)) return
      call coupling_make_view(coupling, 0, view)
      call view%read("density_matrix", "dcart", dcart, err)
      call coupling_close_view(coupling)
      call check(error, allocated(err), more="the update kept the density matrix of the build")
      if (allocated(error)) return
      err%stat = 0
      deallocate (err)
      call check_model_energy_refused(model, coupling, error)
      if (allocated(error)) return
      call model%prepare_build(coupling, err)
      if (failed(error, err, "second build staging")) return
      nvisit = 0
      do while (coupling%next())
         nvisit = nvisit + 1
      end do
      call check(error, nvisit, 2, more="the build phase did not ask for both outputs again")
      if (allocated(error)) return
      call model%release_coupling(coupling)
   end subroutine test_update_drops_build_answers

end module test_model_coupling
