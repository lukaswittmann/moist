!> Construction, lifecycle and solver controls of the PCM components
!>
!> - Constructors, context inheritance, cavity requirements and rejected
!>   material parameters of every model
!> - Guards: missing or stale surface data, energies before `update`, missing
!>   host input, wrong grids, invalid solvers, tangents without current charges
!> - State: reuse across grid changes, re-solves on a fresh potential, the
!>   shared timer tree
!> - Solver choice and potential input: every CPCM solver on neutral and
!>   charged systems, a potential set directly on the request
module test_model_component_pcm_interface
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, &
      & ieee_is_finite
   use mctc_env, only: wp
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: check_moist_error, component_view, build_test_cavity, &
      & fill_point_charge_potential, fill_point_charge_field, set_host_potential, &
      & fill_missing_with_zeros, stage_point_charge_energy
   use test_model_component_pcm_fixtures, only: npcm, new_pcm, make_component, seed_iswig, seed_drop, &
      & fd_cavity, get_seed, pcm_cpcm, pcm_cosmo, pcm_iefpcm, pcm_ssvpe, pcm_lpb, pcm_nonlocal, &
      & pcm_slic, mol_triatomic, get_solute, host_energy, host_charges, stage_energy, &
      & answer_sources
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, &
      & moist_cavity_iswig_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: default_cpcm_radii, radius_type_static, new_cosmo_radii
   use moist_channels_coupling, only: coupling_type, gaussian_potential_request_type
   use moist_channels_response, only: response_type, parameter_adjoint_response_type, &
      & response_clear
   use moist_model_continuum_type, only: model_continuum_type, new_continuum_model
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, &
      & moist_pcm_parameters_type, solver_type
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, &
      & new_component_cpcm, model_continuum_component_cosmo, new_component_cosmo
   use moist_model_continuum_component_pcm_iefpcm, only: model_continuum_component_iefpcm, &
      & new_component_iefpcm
   use moist_model_continuum_component_pcm_ssvpe, only: model_continuum_component_ssvpe, &
      & new_component_ssvpe
   use moist_model_continuum_component_pcm_dpcm, only: model_continuum_component_dpcm, &
      & new_component_dpcm
   use moist_model_continuum_component_pcm_dpcm_original, only: &
      & model_continuum_component_dpcm_original, new_component_dpcm_original
   use moist_model_continuum_component_pcm_lpb, only: model_continuum_component_lpb, &
      & new_component_lpb, debye_kappa
   use moist_model_continuum_component_pcm_slic, only: model_continuum_component_slic, &
      & new_component_slic
   use moist_model_continuum_component_pcm_nonlocal, only: model_continuum_component_nonlocal, &
      & new_component_nonlocal
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_interface

   !> Tolerance of identities that hold to roundoff
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   !> Absolute tolerance for charge and energy identities
   real(wp), parameter :: thr2 = 1.0e-12_wp

contains

   !> Collect the construction, lifecycle and solver-control tests
   !>
   !> @param[out] testsuite Collection of tests
   subroutine collect_model_component_pcm_interface(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("constructors", test_constructors), &
         & new_unittest("optional_context", test_optional_context), &
         & new_unittest("cavity_requirements", test_cavity_requirements), &
         & new_unittest("automatic_cavity_geometry", test_automatic_cavity_geometry), &
         & new_unittest("pcm_rejects_dielectric_below_one", test_pcm_invalid_epsilon), &
         & new_unittest("debye_kappa_reference", test_debye_kappa), &
         & new_unittest("slic_default_mu", test_slic_default_mu), &
         & new_unittest("layer_lifecycle_guards", test_layer_lifecycle_guards), &
         & new_unittest("parameter_derivatives_opt_in", test_parameter_derivatives_opt_in), &
         & new_unittest("conductor_point_cavity_energy", test_conductor_point_cavity), &
         & new_unittest("cpcm_requires_update", test_cpcm_requires_update, should_fail=.true.), &
         & new_unittest("cpcm_external_potential_requires_input", &
         & test_cpcm_external_potential_requires_input, should_fail=.true.), &
         & new_unittest("cpcm_invalid_solver", test_cpcm_invalid_solver, should_fail=.true.), &
         & new_unittest("cpcm_grid_mismatch_guards", test_cpcm_grid_mismatch_guards), &
         & new_unittest("tangent_requires_solution", test_tangent_requires_solution), &
         & new_unittest("gradient_failure_preserves_accumulator", &
         & test_gradient_failure_preserves_accumulator), &
         & new_unittest("cpcm_reallocates_on_grid_change", test_cpcm_reallocate_on_ngrid_change), &
         & new_unittest("cpcm_stale_charge_regression", test_cpcm_stale_charge_regression), &
         & new_unittest("cpcm_records_shared_timer_tree", test_cpcm_timer_tree), &
         & new_unittest("cpcm_lu_energy", test_cpcm_energy_lu), &
         & new_unittest("cpcm_charged_system", test_cpcm_charged), &
         & new_unittest("cpcm_solver_comparison", test_cpcm_solver_comparison), &
         & new_unittest("cpcm_external_potential", test_cpcm_external_potential), &
         & new_unittest("cpcm_accepts_infinite_eps", test_cpcm_accepts_infinite_eps), &
         & new_unittest("lpb_accepts_zero_kappa", test_lpb_accepts_zero_kappa), &
         & new_unittest("slic_flux_regular", test_slic_flux_regular), &
         & new_unittest("slic_recovers_after_bad_trace", test_slic_recovers_after_bad_trace), &
         & new_unittest("lpb_nonfinite_operator", test_lpb_nonfinite_operator, should_fail=.true.), &
         & new_unittest("iefpcm_rejects_eps_below_one", test_iefpcm_rejects_eps_below_one, should_fail=.true.), &
         & new_unittest("lpb_rejects_negative_kappa", test_lpb_rejects_negative_kappa, should_fail=.true.), &
         & new_unittest("nonlocal_rejects_zero_lambda", test_nonlocal_rejects_zero_lambda, should_fail=.true.), &
         & new_unittest("cpcm_rejects_eps_nan", test_cpcm_rejects_eps_nan, should_fail=.true.), &
         & new_unittest("cosmo_rejects_eps_nan", test_cosmo_rejects_eps_nan, should_fail=.true.), &
         & new_unittest("iefpcm_rejects_eps_nan", test_iefpcm_rejects_eps_nan, should_fail=.true.), &
         & new_unittest("ssvpe_rejects_eps_nan", test_ssvpe_rejects_eps_nan, should_fail=.true.), &
         & new_unittest("dpcm_rejects_eps_nan", test_dpcm_rejects_eps_nan, should_fail=.true.), &
         & new_unittest("dpcm_original_rejects_eps_nan", &
         & test_dpcm_original_rejects_eps_nan, should_fail=.true.), &
         & new_unittest("lpb_rejects_eps_nan", test_lpb_rejects_eps_nan, should_fail=.true.), &
         & new_unittest("slic_rejects_eps_nan", test_slic_rejects_eps_nan, should_fail=.true.), &
         & new_unittest("lpb_rejects_kappa_nan", test_lpb_rejects_kappa_nan, should_fail=.true.), &
         & new_unittest("slic_rejects_alpha_nan", test_slic_rejects_alpha_nan, should_fail=.true.), &
         & new_unittest("slic_rejects_beta_nan", test_slic_rejects_beta_nan, should_fail=.true.), &
         & new_unittest("slic_rejects_gamma_nan", test_slic_rejects_gamma_nan, should_fail=.true.), &
         & new_unittest("slic_rejects_mu_nan", test_slic_rejects_mu_nan, should_fail=.true.), &
         & new_unittest("slic_rejects_phi_static_nan", test_slic_rejects_phi_static_nan, should_fail=.true.), &
         & new_unittest("nonlocal_rejects_eps_solvent_nan", &
         & test_nonlocal_rejects_eps_solvent_nan, should_fail=.true.), &
         & new_unittest("nonlocal_rejects_eps_inf_nan", test_nonlocal_rejects_eps_inf_nan, should_fail=.true.), &
         & new_unittest("nonlocal_rejects_lambda_nan", test_nonlocal_rejects_lambda_nan, should_fail=.true.), &
         & new_unittest("nonlocal_rejects_eps_solute_nan", &
         & test_nonlocal_rejects_eps_solute_nan, should_fail=.true.), &
         & new_unittest("lpb_rejects_kappa_inf", test_lpb_rejects_kappa_inf, should_fail=.true.), &
         & new_unittest("slic_rejects_alpha_inf", test_slic_rejects_alpha_inf, should_fail=.true.), &
         & new_unittest("slic_rejects_beta_inf", test_slic_rejects_beta_inf, should_fail=.true.), &
         & new_unittest("slic_rejects_gamma_inf", test_slic_rejects_gamma_inf, should_fail=.true.), &
         & new_unittest("slic_rejects_mu_inf", test_slic_rejects_mu_inf, should_fail=.true.), &
         & new_unittest("slic_rejects_phi_static_inf", test_slic_rejects_phi_static_inf, should_fail=.true.), &
         & new_unittest("nonlocal_rejects_eps_solvent_inf", &
         & test_nonlocal_rejects_eps_solvent_inf, should_fail=.true.), &
         & new_unittest("nonlocal_rejects_eps_inf_inf", test_nonlocal_rejects_eps_inf_inf, should_fail=.true.), &
         & new_unittest("nonlocal_rejects_lambda_inf", test_nonlocal_rejects_lambda_inf, should_fail=.true.), &
         & new_unittest("nonlocal_rejects_eps_solute_inf", &
         & test_nonlocal_rejects_eps_solute_inf, should_fail=.true.), &
         & new_unittest("slic_rejects_singular_flux_boundary", &
         & test_slic_rejects_singular_flux_boundary, should_fail=.true.), &
         & new_unittest("slic_rejects_singular_flux_negative_alpha", &
         & test_slic_rejects_singular_flux_negative_alpha, should_fail=.true.), &
         & new_unittest("slic_rejects_singular_flux_large_alpha", &
         & test_slic_rejects_singular_flux_large_alpha, should_fail=.true.), &
         & new_unittest("slic_rejects_singular_flux_conductor", &
         & test_slic_rejects_singular_flux_conductor, should_fail=.true.), &
         & new_unittest("cpcm_rejects_gaussian_without_switching_update", &
         & test_cpcm_rejects_gaussian_without_switching_update, should_fail=.true.), &
         & new_unittest("cosmo_rejects_gaussian_without_switching_update", &
         & test_cosmo_rejects_gaussian_without_switching_update, should_fail=.true.), &
         & new_unittest("iefpcm_rejects_gaussian_without_switching_update", &
         & test_iefpcm_rejects_gaussian_without_switching_update, should_fail=.true.), &
         & new_unittest("ssvpe_rejects_gaussian_without_switching_update", &
         & test_ssvpe_rejects_gaussian_without_switching_update, should_fail=.true.), &
         & new_unittest("lpb_rejects_gaussian_without_switching_update", &
         & test_lpb_rejects_gaussian_without_switching_update, should_fail=.true.), &
         & new_unittest("nonlocal_rejects_gaussian_without_switching_update", &
         & test_nonlocal_rejects_gaussian_without_switching_update, should_fail=.true.), &
         & new_unittest("cpcm_rejects_gaussian_without_switching_coupling", &
         & test_cpcm_rejects_gaussian_without_switching_coupling, should_fail=.true.), &
         & new_unittest("cosmo_rejects_gaussian_without_switching_coupling", &
         & test_cosmo_rejects_gaussian_without_switching_coupling, should_fail=.true.), &
         & new_unittest("slic_rejects_nan_normal_trace_point", &
         & test_slic_rejects_nan_normal_trace_point, should_fail=.true.), &
         & new_unittest("slic_rejects_inf_normal_trace_point", &
         & test_slic_rejects_inf_normal_trace_point, should_fail=.true.), &
         & new_unittest("slic_rejects_nan_normal_trace_gaussian", &
         & test_slic_rejects_nan_normal_trace_gaussian, should_fail=.true.), &
         & new_unittest("slic_rejects_inf_normal_trace_gaussian", &
         & test_slic_rejects_inf_normal_trace_gaussian, should_fail=.true.) &
         & ]
   end subroutine collect_model_component_pcm_interface

   !* ================================================================================= *!
   !*                                   Construction                                    *!
   !* ================================================================================= *!

   !> Verify each concrete type copies the shared solver controls
   subroutine test_constructors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moist_pcm_parameters_type) :: settings
      type(model_continuum_component_iefpcm) :: ief
      type(model_continuum_component_ssvpe) :: ssv
      type(model_continuum_component_dpcm) :: dpcm
      type(model_continuum_component_dpcm_original) :: dpcm0
      type(model_continuum_component_lpb) :: lpb
      type(model_continuum_component_slic) :: slic
      type(model_continuum_component_nonlocal) :: nonlocal

      call new_context(ctx, verbosity=0)
      settings%solver = solver_type%iterative
      settings%solver_tol = 2.0e-9_wp
      settings%solver_maxiter = 37
      settings%solver_restart = 11
      settings%newton_tol = 3.0e-11_wp
      settings%newton_maxiter = 17
      call new_component_iefpcm(ief, 80.0_wp, err, settings, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, ief%param%solver_restart, 11)
      if (allocated(error)) return
      call new_component_ssvpe(ssv, 80.0_wp, err, settings, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, ssv%param%solver_maxiter, 37)
      if (allocated(error)) return
      call new_component_dpcm(dpcm, 80.0_wp, err, settings, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, dpcm%uses_normal_trace())
      if (allocated(error)) return
      call new_component_dpcm_original(dpcm0, 80.0_wp, err, settings, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, dpcm0%uses_normal_trace())
      if (allocated(error)) return
      call new_component_lpb(lpb, 1.0_wp, 0.1_wp, err, settings, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call new_component_slic(slic, 80.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, err, settings, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, slic%param%newton_maxiter, 17)
      if (allocated(error)) return
      call new_component_nonlocal(nonlocal, 80.0_wp, 2.0_wp, 1.0_wp, err, settings, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, nonlocal%uses_normal_trace())
      if (allocated(error)) return
      call check(error, nonlocal%needs_curvature)
      if (allocated(error)) return
      call ctx%delete()
   end subroutine test_constructors

   !> All PCM models inherit missing contexts and preserve explicit ones
   subroutine test_optional_context(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library diagnostic
      type(moist_error_type), allocatable :: err
      !> Contexts owned by this test
      type(moist_context_type), target :: model_ctx, own_ctx
      !> Disassociated actual argument, as passed by the C binding for NULL
      type(moist_context_type), pointer :: absent_ctx
      !> Model and minimal cavity used only for context inheritance
      type(model_continuum_type) :: model
      type(cavity_type_iswig) :: cavity
      !> Empty structure; the missing-context guard must precede geometry access
      type(structure_type) :: mol
      !> Concrete PCM component under test
      class(model_continuum_component_pcm), allocatable :: component
      !> Constructor and context mode
      integer :: kind, mode

      nullify (absent_ctx)
      call new_context(model_ctx, verbosity=0)
      call new_context(own_ctx, verbosity=0)
      do mode = 1, 3
         do kind = 1, 9
            select case (mode)
            case (1)
               call make_component(kind, component, err)
            case (2)
               call make_component(kind, component, err, ctx=own_ctx)
            case (3)
               call make_component(kind, component, err, ctx=absent_ctx)
            case default
               call test_failed(error, "test_optional_context: invalid mode")
               return
            end select
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            if (mode == 2) then
               call check(error, associated(component%ctx, own_ctx))
            else
               call check(error,.not. associated(component%ctx))
               if (allocated(error)) return
               call component%update(mol, cavity, err)
               call check(error, allocated(err))
               if (allocated(error)) return
               call check(error, index(err%message, "Component has no context") > 0)
            end if
            if (allocated(error)) return
            call new_continuum_model(model, model_ctx, cavity, err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            call model%add_component(component, err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            if (mode == 2) then
               call check(error, associated(model%components(1)%item%ctx, own_ctx))
            else
               call check(error, associated(model%components(1)%item%ctx, model_ctx))
               if (allocated(error)) return
               call check(error,.not. associated(component%ctx))
            end if
            if (allocated(error)) return
            call model%components(1)%item%require_context(err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            deallocate (component)
         end do
      end do

   end subroutine test_optional_context

   !> Required geometry is enabled on model copies without clearing requested extras
   subroutine test_cavity_requirements(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library diagnostic
      type(moist_error_type), allocatable :: err
      !> Context owned by this test
      type(moist_context_type), target :: ctx
      !> Original cavity configurations
      type(cavity_type_drop) :: drop
      type(cavity_type_iswig) :: iswig
      !> Model owns copies of each cavity and component
      type(model_continuum_type) :: model
      !> PCM family under test
      class(model_continuum_component_pcm), allocatable :: component
      !> Constructor and pre-enabled geometry selection
      integer :: kind, preset
      !> Geometry enabled on the original cavity
      logical :: had_normal, had_curvature

      call new_context(ctx, verbosity=0)
      do preset = 0, 2
         had_normal = preset == 1
         had_curvature = preset == 2
         do kind = 1, 9
            call make_component(kind, component, err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            drop%param%do_normal = had_normal
            drop%param%do_curvature = had_curvature
            drop%param%do_rho = .true.
            call new_continuum_model(model, ctx, drop, err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            call model%add_component(component, err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            call model%components(1)%item%configure_cavity(model%cavity, err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            select type (cavity => model%cavity)
            type is (cavity_type_drop)
               call check(error, cavity%param%do_normal .eqv. (had_normal .or. kind >= 3))
               if (allocated(error)) return
               call check(error, cavity%param%do_curvature .eqv. (had_curvature .or. kind >= 3))
               if (allocated(error)) return
               call check(error, cavity%param%do_rho .and. .not. cavity%param%do_fine &
                  & .and. .not. cavity%param%do_grid_density .and. .not. cavity%param%do_r_iI)
               if (allocated(error)) return
               call cavity%configure(normal=.false., curvature=.false., error=err)
               call check(error,.not. allocated(err))
               if (allocated(error)) return
               call check(error, cavity%param%do_normal .eqv. (had_normal .or. kind >= 3))
               if (allocated(error)) return
               call check(error, cavity%param%do_curvature .eqv. (had_curvature .or. kind >= 3))
            end select
            if (allocated(error)) return
            call check(error, (drop%param%do_normal .eqv. had_normal) &
               & .and. (drop%param%do_curvature .eqv. had_curvature))
            if (allocated(error)) return

            iswig%need_curvature = had_curvature
            call new_continuum_model(model, ctx, iswig, err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            call model%add_component(component, err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            call model%components(1)%item%configure_cavity(model%cavity, err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            select type (cavity => model%cavity)
            type is (cavity_type_iswig)
               call check(error, cavity%need_curvature .eqv. (had_curvature .or. kind >= 3))
               if (allocated(error)) return
               call cavity%configure(curvature=.false., error=err)
               call check(error,.not. allocated(err))
               if (allocated(error)) return
               call check(error, cavity%need_curvature .eqv. (had_curvature .or. kind >= 3))
            end select
            if (allocated(error)) return
            call check(error, iswig%need_curvature .eqv. had_curvature)
            if (allocated(error)) return
            deallocate (component)
         end do
      end do
   end subroutine test_cavity_requirements

   !> Model update computes required normals and curvatures without do_fine
   subroutine test_automatic_cavity_geometry(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library diagnostic
      type(moist_error_type), allocatable :: err
      !> Context owned by this test
      type(moist_context_type), target :: ctx
      !> One-sphere geometry
      type(structure_type) :: mol
      !> Original cavities, created without optional geometry
      type(cavity_type_drop) :: drop
      type(cavity_type_iswig) :: iswig
      !> Smooth vdW surface
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Layer PCM requires both normals and curvature
      type(model_continuum_component_iefpcm) :: component
      !> Model owns its cavity copy
      type(model_continuum_type) :: model
      !> Cavity and repeated update indices
      integer :: kind, repeat

      call new_context(ctx, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call lsf%new()
      call new_cavity_drop(drop, radius_model=default_cpcm_radii(), lsf_model=lsf, error=err, &
         & param=moist_cavity_drop_parameters_type(num_leb=14, do_rho=.true.))
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call new_cavity_iswig(iswig, radius_model=default_cpcm_radii(), error=err, &
         & param=moist_cavity_iswig_parameters_type(num_leb=14))
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call new_component_iefpcm(component, 80.0_wp, error=err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      do kind = 1, 2
         if (kind == 1) then
            call new_continuum_model(model, ctx, drop, err)
         else
            call new_continuum_model(model, ctx, iswig, err)
         end if
         call check(error,.not. allocated(err))
         if (allocated(error)) return
         call model%add_component(component, err)
         call check(error,.not. allocated(err))
         if (allocated(error)) return
         do repeat = 1, 2
            call model%update(mol, err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
            associate (cavity => model%cavity)
               call check(error, allocated(cavity%normal0) .and. allocated(cavity%k1) .and. allocated(cavity%k2))
               if (allocated(error)) return
               call check(error, size(cavity%normal0, 2) == cavity%ngrid .and. size(cavity%k1) == cavity%ngrid)
               if (allocated(error)) return
               call check(error, all(abs(sum(cavity%normal0**2, dim=1) - 1.0_wp) < 1.0e-10_wp))
               if (allocated(error)) return
               call check(error, all(abs(cavity%k1 - 1.0_wp/cavity%radii(cavity%owner)) < 1.0e-10_wp))
               if (allocated(error)) return
               call check(error, all(abs(cavity%k2 - cavity%k1) < 1.0e-10_wp))
               if (allocated(error)) return
            end associate
            select type (cavity => model%cavity)
            type is (cavity_type_drop)
               call check(error, cavity%param%do_rho .and. allocated(cavity%rho))
               if (allocated(error)) return
               call check(error,.not. cavity%param%do_fine)
               if (allocated(error)) return
            end select
         end do
      end do
      call check(error,.not. drop%param%do_normal .and. .not. drop%param%do_curvature)
      if (allocated(error)) return
      call check(error,.not. allocated(drop%normal0) .and. .not. allocated(drop%k1))
      if (allocated(error)) return
      call check(error,.not. iswig%need_curvature .and. .not. allocated(iswig%k1))
   end subroutine test_automatic_cavity_geometry

   !> IEF-PCM rejects eps < 1
   subroutine test_iefpcm_rejects_eps_below_one(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_iefpcm) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_iefpcm(pcm, 0.5_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_iefpcm_rejects_eps_below_one

   !> LPB rejects a negative inverse Debye length
   subroutine test_lpb_rejects_negative_kappa(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_lpb) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_lpb(pcm, 80.0_wp, -1.0_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_lpb_rejects_negative_kappa

   !> Nonlocal PCM rejects a zero correlation length
   subroutine test_nonlocal_rejects_zero_lambda(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_nonlocal) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_nonlocal(pcm, 80.0_wp, 2.0_wp, 0.0_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_nonlocal_rejects_zero_lambda

   !> CPCM rejects eps = NaN
   subroutine test_cpcm_rejects_eps_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_cpcm) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_cpcm(pcm, ieee_value(1.0_wp, ieee_quiet_nan), error=err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_cpcm_rejects_eps_nan

   !> COSMO rejects eps = NaN
   subroutine test_cosmo_rejects_eps_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_cosmo) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_cosmo(pcm, ieee_value(1.0_wp, ieee_quiet_nan), error=err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_cosmo_rejects_eps_nan

   !> IEF-PCM rejects eps = NaN
   subroutine test_iefpcm_rejects_eps_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_iefpcm) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_iefpcm(pcm, ieee_value(1.0_wp, ieee_quiet_nan), err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_iefpcm_rejects_eps_nan

   !> SS(V)PE rejects eps = NaN
   subroutine test_ssvpe_rejects_eps_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_ssvpe) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_ssvpe(pcm, ieee_value(1.0_wp, ieee_quiet_nan), err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_ssvpe_rejects_eps_nan

   !> DPCM rejects eps = NaN
   subroutine test_dpcm_rejects_eps_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_dpcm) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_dpcm(pcm, ieee_value(1.0_wp, ieee_quiet_nan), err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_dpcm_rejects_eps_nan

   !> DPCM original rejects eps = NaN
   subroutine test_dpcm_original_rejects_eps_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_dpcm_original) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_dpcm_original(pcm, ieee_value(1.0_wp, ieee_quiet_nan), err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_dpcm_original_rejects_eps_nan

   !> LPB rejects eps = NaN
   subroutine test_lpb_rejects_eps_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_lpb) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_lpb(pcm, ieee_value(1.0_wp, ieee_quiet_nan), 0.1_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_lpb_rejects_eps_nan

   !> SLIC rejects eps = NaN
   subroutine test_slic_rejects_eps_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, ieee_value(1.0_wp, ieee_quiet_nan), 0.1_wp, 0.2_wp, 0.3_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_eps_nan

   !> LPB rejects kappa = NaN
   subroutine test_lpb_rejects_kappa_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_lpb) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_lpb(pcm, 80.0_wp, ieee_value(1.0_wp, ieee_quiet_nan), err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_lpb_rejects_kappa_nan

   !> SLIC rejects alpha = NaN
   subroutine test_slic_rejects_alpha_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, ieee_value(1.0_wp, ieee_quiet_nan), 0.2_wp, 0.3_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_alpha_nan

   !> SLIC rejects beta = NaN
   subroutine test_slic_rejects_beta_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, 0.1_wp, ieee_value(1.0_wp, ieee_quiet_nan), 0.3_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_beta_nan

   !> SLIC rejects gamma = NaN
   subroutine test_slic_rejects_gamma_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, 0.1_wp, 0.2_wp, ieee_value(1.0_wp, ieee_quiet_nan), err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_gamma_nan

   !> SLIC rejects mu = NaN
   subroutine test_slic_rejects_mu_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, err, mu=ieee_value(1.0_wp, &
         & ieee_quiet_nan), ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_mu_nan

   !> SLIC rejects phi_static = NaN
   subroutine test_slic_rejects_phi_static_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, err, &
         & phi_static=ieee_value(1.0_wp, ieee_quiet_nan), ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_phi_static_nan

   !> Nonlocal PCM rejects eps_solvent = NaN
   subroutine test_nonlocal_rejects_eps_solvent_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_nonlocal) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_nonlocal(pcm, ieee_value(1.0_wp, ieee_quiet_nan), 2.0_wp, 1.0_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_nonlocal_rejects_eps_solvent_nan

   !> Nonlocal PCM rejects eps_inf = NaN
   subroutine test_nonlocal_rejects_eps_inf_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_nonlocal) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_nonlocal(pcm, 80.0_wp, ieee_value(1.0_wp, ieee_quiet_nan), 1.0_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_nonlocal_rejects_eps_inf_nan

   !> Nonlocal PCM rejects lambda = NaN
   subroutine test_nonlocal_rejects_lambda_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_nonlocal) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_nonlocal(pcm, 80.0_wp, 2.0_wp, ieee_value(1.0_wp, ieee_quiet_nan), err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_nonlocal_rejects_lambda_nan

   !> Nonlocal PCM rejects eps_solute = NaN
   subroutine test_nonlocal_rejects_eps_solute_nan(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_nonlocal) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_nonlocal(pcm, 80.0_wp, 2.0_wp, 1.0_wp, err, epsilon_solute=ieee_value(1.0_wp, &
         & ieee_quiet_nan), ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_nonlocal_rejects_eps_solute_nan

   !> LPB rejects kappa = +Inf
   subroutine test_lpb_rejects_kappa_inf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_lpb) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_lpb(pcm, 80.0_wp, ieee_value(1.0_wp, ieee_positive_inf), err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_lpb_rejects_kappa_inf

   !> SLIC rejects alpha = +Inf
   subroutine test_slic_rejects_alpha_inf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, ieee_value(1.0_wp, ieee_positive_inf), 0.2_wp, 0.3_wp, &
         & err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_alpha_inf

   !> SLIC rejects beta = +Inf
   subroutine test_slic_rejects_beta_inf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, 0.1_wp, ieee_value(1.0_wp, ieee_positive_inf), 0.3_wp, &
         & err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_beta_inf

   !> SLIC rejects gamma = +Inf
   subroutine test_slic_rejects_gamma_inf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, 0.1_wp, 0.2_wp, ieee_value(1.0_wp, ieee_positive_inf), &
         & err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_gamma_inf

   !> SLIC rejects mu = +Inf
   subroutine test_slic_rejects_mu_inf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, err, mu=ieee_value(1.0_wp, &
         & ieee_positive_inf), ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_mu_inf

   !> SLIC rejects phi_static = +Inf
   subroutine test_slic_rejects_phi_static_inf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, err, &
         & phi_static=ieee_value(1.0_wp, ieee_positive_inf), ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_phi_static_inf

   !> Nonlocal PCM rejects eps_solvent = +Inf
   subroutine test_nonlocal_rejects_eps_solvent_inf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_nonlocal) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_nonlocal(pcm, ieee_value(1.0_wp, ieee_positive_inf), 2.0_wp, 1.0_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_nonlocal_rejects_eps_solvent_inf

   !> Nonlocal PCM rejects eps_inf = +Inf
   subroutine test_nonlocal_rejects_eps_inf_inf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_nonlocal) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_nonlocal(pcm, 80.0_wp, ieee_value(1.0_wp, ieee_positive_inf), 1.0_wp, err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_nonlocal_rejects_eps_inf_inf

   !> Nonlocal PCM rejects lambda = +Inf
   subroutine test_nonlocal_rejects_lambda_inf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_nonlocal) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_nonlocal(pcm, 80.0_wp, 2.0_wp, ieee_value(1.0_wp, ieee_positive_inf), err, ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_nonlocal_rejects_lambda_inf

   !> Nonlocal PCM rejects eps_solute = +Inf
   subroutine test_nonlocal_rejects_eps_solute_inf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_nonlocal) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_nonlocal(pcm, 80.0_wp, 2.0_wp, 1.0_wp, err, epsilon_solute=ieee_value(1.0_wp, &
         & ieee_positive_inf), ctx=ctx)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_nonlocal_rejects_eps_solute_inf

   !> CPCM accepts the conductor limit eps = +Inf
   subroutine test_cpcm_accepts_infinite_eps(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_cpcm) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_cpcm(pcm, ieee_value(1.0_wp, ieee_positive_inf), error=err, ctx=ctx)
      call check_moist_error(error, err, "CPCM eps = +Inf")
      call ctx%delete()
   end subroutine test_cpcm_accepts_infinite_eps

   !> LPB accepts kappa = 0, the electrolyte-free limit
   subroutine test_lpb_accepts_zero_kappa(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_lpb) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_lpb(pcm, 80.0_wp, 0.0_wp, err, ctx=ctx)
      call check_moist_error(error, err, "LPB kappa = 0")
      call ctx%delete()
   end subroutine test_lpb_accepts_zero_kappa

   !> Debye inverse screening length against the SI formula and its scalings
   !>
   !> Reference from CODATA 2018 SI constants, kappa = sqrt(2 e^2 NA I/(eps eps0 kB T)) in bohr^-1
   subroutine test_debye_kappa(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> 78.4, 0.1 mol/L, 298.15 K
      real(wp), parameter :: ref = 0.0550440690847816_wp
      real(wp), parameter :: tol = 1.0e-9_wp
      real(wp) :: k0

      k0 = debye_kappa(78.4_wp, 0.1_wp, 298.15_wp)
      call check(error, abs(k0 - ref)/ref < tol, message="Debye kappa deviates from the SI value")
      if (allocated(error)) return
      call check(error, abs(debye_kappa(78.4_wp, 0.4_wp, 298.15_wp)/k0 - 2.0_wp) < 1.0e-12_wp, &
         & message="Debye kappa is not proportional to sqrt(I)")
      if (allocated(error)) return
      call check(error, abs(debye_kappa(4*78.4_wp, 0.1_wp, 298.15_wp)/k0 - 0.5_wp) < 1.0e-12_wp, &
         & message="Debye kappa is not proportional to 1/sqrt(eps)")
      if (allocated(error)) return
      call check(error, abs(debye_kappa(78.4_wp, 0.1_wp, 4*298.15_wp)/k0 - 0.5_wp) < 1.0e-12_wp, &
         & message="Debye kappa is not proportional to 1/sqrt(T)")
      if (allocated(error)) return
      call check(error, debye_kappa(78.4_wp, 0.0_wp, 298.15_wp) == 0.0_wp, &
         & message="Debye kappa does not vanish without electrolyte")
   end subroutine test_debye_kappa


   !> Omitted SLIC mu makes the perturbation vanish at zero field, h(0) = 0
   subroutine test_slic_default_mu(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moist_pcm_parameters_type) :: param
      type(model_continuum_component_slic) :: omitted, explicit, plain
      type(structure_type) :: mol
      type(fd_cavity) :: cavity
      real(wp), parameter :: alpha = 0.3_wp, beta = -2.0_wp, gamma = -0.4_wp
      real(wp) :: e_omitted, e_explicit, e_plain

      call new_context(ctx, verbosity=0)
      param%solver = solver_type%dense
      call new_component_slic(omitted, 80.0_wp, alpha, beta, gamma, err, param, ctx=ctx)
      call check_moist_error(error, err, "slic omitted mu")
      if (allocated(error)) return
      ! h(E=0) = alpha*tanh(-gamma) + mu
      call check(error, alpha*tanh(-gamma) + omitted%mu, 0.0_wp, thr=1.0e-15_wp, &
         & message="SLIC perturbation does not vanish at zero field")
      if (allocated(error)) return
      call check(error, omitted%mu /= 0.0_wp, message="SLIC default mu is zero")
      if (allocated(error)) return
      call new_component_slic(explicit, 80.0_wp, alpha, beta, gamma, err, param, &
         & mu=alpha*tanh(gamma), ctx=ctx)
      call check_moist_error(error, err, "slic explicit mu")
      if (allocated(error)) return
      call new_component_slic(plain, 80.0_wp, alpha, beta, gamma, err, param, mu=0.0_wp, ctx=ctx)
      call check_moist_error(error, err, "slic zero mu")
      if (allocated(error)) return

      call get_solute(mol_triatomic, mol)
      call get_seed(seed_drop, .true., mol, cavity, error)
      if (allocated(error)) return
      call host_energy(omitted, mol, cavity, mol, host_charges(mol), e_omitted, error)
      if (allocated(error)) return
      call host_energy(explicit, mol, cavity, mol, host_charges(mol), e_explicit, error)
      if (allocated(error)) return
      call host_energy(plain, mol, cavity, mol, host_charges(mol), e_plain, error)
      if (allocated(error)) return
      call check(error, e_omitted, e_explicit, thr=1.0e-10_wp, &
         & message="SLIC omitted mu differs from alpha*tanh(gamma)")
      if (allocated(error)) return
      call check(error, abs(e_omitted - e_plain) > 1.0e-6_wp, &
         & message="SLIC omitted mu is indistinguishable from mu = 0")
      call ctx%delete()
   end subroutine test_slic_default_mu

   !> SLIC rejects the singular flux boundary 2*feps*(mu + |alpha|) = 1
   subroutine test_slic_rejects_singular_flux_boundary(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 3.0_wp, 0.5_wp, 0.2_wp, 0.3_wp, err, mu=0.5_wp)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_singular_flux_boundary

   !> SLIC rejects a singular flux with alpha < 0
   subroutine test_slic_rejects_singular_flux_negative_alpha(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 3.0_wp, -0.5_wp, 0.2_wp, 0.3_wp, err, mu=0.51_wp)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_singular_flux_negative_alpha

   !> SLIC rejects alpha = 1 at eps = 80 without a shift, 2*feps > 1
   subroutine test_slic_rejects_singular_flux_large_alpha(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, 80.0_wp, 1.0_wp, 0.2_wp, 0.3_wp, err, mu=0.0_wp)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_singular_flux_large_alpha

   !> SLIC rejects the singular conductor flux, boundary mu + |alpha| = 1/2
   subroutine test_slic_rejects_singular_flux_conductor(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: pcm

      call new_context(ctx, verbosity=0)
      call new_component_slic(pcm, ieee_value(1.0_wp, ieee_positive_inf), 0.2_wp, 0.2_wp, 0.3_wp, &
         & err, mu=0.3_wp)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine test_slic_rejects_singular_flux_conductor

   !> SLIC flux guard 2*feps*(mu + |alpha|) < 1 accepts regular parameters
   subroutine test_slic_flux_regular(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_slic) :: slic

      call new_context(ctx, verbosity=0)
      ! eps = 3: feps = 1/2, boundary at mu + |alpha| = 1
      call new_component_slic(slic, 3.0_wp, 0.5_wp, 0.2_wp, 0.3_wp, err, mu=0.49_wp, ctx=ctx)
      call check_moist_error(error, err, "slic below the flux boundary")
      if (allocated(error)) return
      call new_component_slic(slic, 3.0_wp, -0.5_wp, 0.2_wp, 0.3_wp, err, mu=0.49_wp, ctx=ctx)
      call check_moist_error(error, err, "slic below the flux boundary, alpha < 0")
      if (allocated(error)) return
      ! a negative shift lowers the flux: 2*(79/81)*(1 - 0.5) < 1
      call new_component_slic(slic, 80.0_wp, 1.0_wp, 0.2_wp, 0.3_wp, err, mu=-0.5_wp, ctx=ctx)
      call check_moist_error(error, err, "slic with a negative shift")
      if (allocated(error)) return
      ! conductor limit: feps = 1, boundary at mu + |alpha| = 1/2
      call new_component_slic(slic, ieee_value(1.0_wp, ieee_positive_inf), 0.2_wp, 0.2_wp, 0.3_wp, &
         & err, mu=0.29_wp, ctx=ctx)
      call check_moist_error(error, err, "slic conductor below the flux boundary")
      call ctx%delete()
   end subroutine test_slic_flux_regular

!> Every PCM flavor must reject a dielectric below one and accept eps = 1
!> Below eps = 1 the scaling factor f(eps) = (eps-1)/eps turns negative (and
!> diverges at eps = 0), so the model is undefined there; eps = 1 is the
!> vacuum limit, a valid model with f(eps) = 0
   subroutine test_pcm_invalid_epsilon(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(model_continuum_component_cpcm) :: pcm_model
      type(model_continuum_component_cosmo) :: cosmo_model
      !> Dielectric constants that must be rejected
      real(wp), parameter :: bad_epsilon(*) = [0.0_wp, 0.5_wp, -1.0_wp]
      !> Index over the rejected dielectric constants
      integer :: ieps

      !> Run context owned here and borrowed by the components
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      do ieps = 1, size(bad_epsilon)
         call new_component_cpcm(pcm_model, bad_epsilon(ieps), error=err, ctx=ctx)
         call check(error, allocated(err), "CPCM accepted a dielectric below one")
         if (allocated(error)) return
         call check(error, index(err%message, "must be >= 1") > 0, &
            & "CPCM rejected the dielectric with an unexpected message: "//err%message)
         if (allocated(error)) return
         deallocate (err)

         call new_component_cosmo(cosmo_model, bad_epsilon(ieps), error=err, ctx=ctx)
         call check(error, allocated(err), "COSMO accepted a dielectric below one")
         if (allocated(error)) return
         call check(error, index(err%message, "must be >= 1") > 0, &
            & "COSMO rejected the dielectric with an unexpected message: "//err%message)
         if (allocated(error)) return
         deallocate (err)
      end do

      ! The vacuum limit is a valid model, not an error
      call new_component_cpcm(pcm_model, 1.0_wp, error=err, ctx=ctx)
      call check(error,.not. allocated(err), "CPCM rejected the vacuum limit eps = 1")
      if (allocated(error)) return
      call check(error, pcm_model%feps, 0.0_wp, thr=thr, &
         & message="CPCM f(eps) is not zero at eps = 1")
      if (allocated(error)) return

      call new_component_cosmo(cosmo_model, 1.0_wp, error=err, ctx=ctx)
      call check(error,.not. allocated(err), "COSMO rejected the vacuum limit eps = 1")
      if (allocated(error)) return
      call check(error, cosmo_model%feps, 0.0_wp, thr=thr, &
         & message="COSMO f(eps) is not zero at eps = 1")

   end subroutine test_pcm_invalid_epsilon

   !* ================================================================================= *!
   !*                                      Guards                                       *!
   !* ================================================================================= *!

   !> Refuse missing or stale surface data and unsolved adjoints
   subroutine test_layer_lifecycle_guards(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(model_continuum_component_iefpcm) :: component
      type(cavity_type_drop) :: cavity
      type(coupling_type), target :: coupling
      type(cavity_surface_adjoint_type) :: weights
      real(wp), allocatable :: adjoint(:)
      real(wp) :: energy

      call new_context(ctx, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_component_iefpcm(component, 80.0_wp, err, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call component%potential_adjoint(adjoint, err)
      call expect_error("adjoint is unavailable")
      if (allocated(error)) return
      call component%normal_trace_adjoint(adjoint, err)
      call expect_error("normal adjoint is unavailable")
      if (allocated(error)) return
      energy = 7.0_wp
      call component%get_energy(component_view(coupling), cavity, energy, err)
      call expect_error("updated before solving")
      if (allocated(error)) return
      call check(error, energy, 7.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return

      cavity%ngrid = 1
      cavity%nsph = 1
      cavity%updated = .true.
      cavity%built = .true.
      call component%update(mol, cavity, err)
      call expect_error("positions, areas, and normals")
      if (allocated(error)) return
      cavity%xyz = reshape([2.0_wp, 0.0_wp, 0.0_wp], [3, 1])
      cavity%normal0 = reshape([1.0_wp, 0.0_wp, 0.0_wp], [3, 1])
      cavity%a = [1.0_wp]
      call component%update(mol, cavity, err)
      call expect_error("principal curvatures")
      if (allocated(error)) return
      cavity%k1 = [0.4_wp]
      cavity%k2 = [0.5_wp]
      cavity%ngrid = 2
      call component%update(mol, cavity, err)
      call expect_error("geometry is stale")
      if (allocated(error)) return
      cavity%xyz = spread(cavity%xyz(:, 1), 2, 2)
      cavity%normal0 = spread(cavity%normal0(:, 1), 2, 2)
      cavity%a = [1.0_wp, 1.0_wp]
      call component%update(mol, cavity, err)
      call expect_error("curvature is stale")
      if (allocated(error)) return
      cavity%ngrid = 1
      cavity%xi0 = [0.7_wp]
      call component%update(mol, cavity, err)
      call expect_error("requires switching factors")
      if (allocated(error)) return
      cavity%f = [1.0_wp]
      call component%update(mol, cavity, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call component%new_coupling(cavity, coupling, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call component%prepare_energy(cavity, coupling, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call set_host_potential(coupling, [0.1_wp])
      cavity%ngrid = 2
      call component%get_energy(component_view(coupling), cavity, energy, err)
      call expect_error("potential grid mismatch")
      if (allocated(error)) return
      call check(error, energy, 7.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return

      cavity%ngrid = 1
      call component%prepare_gradient(cavity, coupling, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            call coupling%answer("dphi_dr", reshape([0.1_wp, 0.2_wp, 0.3_wp], [3, 1]), err)
            call check(error,.not. allocated(err))
            if (allocated(error)) return
         end select
      end do
      call weights%init(2)
      call component%get_host_surface_weights(component_view(coupling), weights, 2, err)
      call expect_error("spatial derivative grid mismatch")
      if (allocated(error)) return
      call weights%init(1)
      call component%get_host_surface_weights(component_view(coupling), weights, 1, err)
      call expect_error("dphi_dxi")

   contains

      !> Match the guard diagnostic before the next rejected operation
      !>
      !> @param[in] message Expected diagnostic substring
      subroutine expect_error(message)
         !> Expected diagnostic substring
         character(len=*), intent(in) :: message
         call check(error, allocated(err), message="Expected layer guard: "//message)
         if (allocated(error)) return
         call check(error, index(err%message, message) > 0, message=err%message)
      end subroutine expect_error

   end subroutine test_layer_lifecycle_guards

   !> CPCM refuses Gaussian widths without switching factors at `update`
   subroutine test_cpcm_rejects_gaussian_without_switching_update(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call reject_incomplete_update(error, pcm_cpcm)
   end subroutine test_cpcm_rejects_gaussian_without_switching_update

   !> COSMO refuses Gaussian widths without switching factors at `update`
   subroutine test_cosmo_rejects_gaussian_without_switching_update(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call reject_incomplete_update(error, pcm_cosmo)
   end subroutine test_cosmo_rejects_gaussian_without_switching_update

   !> IEF-PCM refuses Gaussian widths without switching factors at `update`
   subroutine test_iefpcm_rejects_gaussian_without_switching_update(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call reject_incomplete_update(error, pcm_iefpcm)
   end subroutine test_iefpcm_rejects_gaussian_without_switching_update

   !> SS(V)PE refuses Gaussian widths without switching factors at `update`
   subroutine test_ssvpe_rejects_gaussian_without_switching_update(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call reject_incomplete_update(error, pcm_ssvpe)
   end subroutine test_ssvpe_rejects_gaussian_without_switching_update

   !> LPB refuses Gaussian widths without switching factors at `update`
   subroutine test_lpb_rejects_gaussian_without_switching_update(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call reject_incomplete_update(error, pcm_lpb)
   end subroutine test_lpb_rejects_gaussian_without_switching_update

   !> Nonlocal PCM refuses Gaussian widths without switching factors at `update`
   subroutine test_nonlocal_rejects_gaussian_without_switching_update(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call reject_incomplete_update(error, pcm_nonlocal)
   end subroutine test_nonlocal_rejects_gaussian_without_switching_update

   !> CPCM refuses Gaussian widths without switching factors at `declare_coupling`
   subroutine test_cpcm_rejects_gaussian_without_switching_coupling(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call reject_incomplete_coupling(error, pcm_cpcm)
   end subroutine test_cpcm_rejects_gaussian_without_switching_coupling

   !> COSMO refuses Gaussian widths without switching factors at `declare_coupling`
   subroutine test_cosmo_rejects_gaussian_without_switching_coupling(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call reject_incomplete_coupling(error, pcm_cosmo)
   end subroutine test_cosmo_rejects_gaussian_without_switching_coupling

   !> Gaussian cavity with widths but without switching factors
   !>
   !> @param[out] mol    solute
   !> @param[out] cavity incomplete cavity
   subroutine incomplete_gaussian_cavity(mol, cavity)
      !> Solute
      type(structure_type), intent(out) :: mol
      !> Incomplete cavity
      type(cavity_type_drop), intent(out) :: cavity

      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      cavity%ngrid = 1
      cavity%nsph = 1
      cavity%updated = .true.
      cavity%built = .true.
      cavity%xyz = reshape([2.0_wp, 0.0_wp, 0.0_wp], [3, 1])
      cavity%normal0 = reshape([1.0_wp, 0.0_wp, 0.0_wp], [3, 1])
      cavity%a = [1.0_wp]
      cavity%k1 = [0.4_wp]
      cavity%k2 = [0.5_wp]
      cavity%xi0 = [0.7_wp]
   end subroutine incomplete_gaussian_cavity

   !> Update of one model on the incomplete cavity, the single call under test
   !>
   !> A failing setup returns without error and so reports an unexpected pass
   !>
   !> @param[out] error test failure
   !> @param[in]  kind  model identifier
   subroutine reject_incomplete_update(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind

      type(moist_error_type), allocatable :: err, setup
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(model_continuum_component_pcm), allocatable :: component
      type(cavity_type_drop) :: cavity

      call new_context(ctx, verbosity=0)
      call incomplete_gaussian_cavity(mol, cavity)
      call make_component(kind, component, setup, ctx)
      if (allocated(setup)) return
      call component%update(mol, cavity, err)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine reject_incomplete_update

   !> Coupling declaration of one model on the incomplete cavity, the single call under test
   !>
   !> A failing setup returns without error and so reports an unexpected pass
   !>
   !> @param[out] error test failure
   !> @param[in]  kind  model identifier
   subroutine reject_incomplete_coupling(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind

      type(moist_error_type), allocatable :: err, setup
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(model_continuum_component_pcm), allocatable :: component
      type(cavity_type_drop) :: cavity
      type(coupling_type), target :: coupling

      call new_context(ctx, verbosity=0)
      call incomplete_gaussian_cavity(mol, cavity)
      call make_component(kind, component, setup, ctx)
      if (allocated(setup)) return
      call component%new_coupling(cavity, coupling, err)
      call check_moist_error(error, err)
      call ctx%delete()
   end subroutine reject_incomplete_coupling

   !> Conductor models on a point cavity without widths still update and solve
   subroutine test_conductor_point_cavity(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(model_continuum_component_pcm), allocatable :: component
      type(cavity_type_drop) :: cavity
      type(coupling_type), target :: coupling
      real(wp) :: energy
      integer :: i

      call new_context(ctx, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      cavity%ngrid = 1
      cavity%nsph = 1
      cavity%updated = .true.
      cavity%built = .true.
      cavity%xyz = reshape([2.0_wp, 0.0_wp, 0.0_wp], [3, 1])
      cavity%normal0 = reshape([1.0_wp, 0.0_wp, 0.0_wp], [3, 1])
      cavity%a = [1.0_wp]
      cavity%k1 = [0.4_wp]
      cavity%k2 = [0.5_wp]

      do i = 1, 2
         call make_component(i, component, err, ctx)
         call check(error,.not. allocated(err), "conductor construction")
         if (allocated(error)) return
         call component%update(mol, cavity, err)
         call check(error,.not. allocated(err), "point cavity rejected by update")
         if (allocated(error)) return
         call component%new_coupling(cavity, coupling, err)
         call check(error,.not. allocated(err), "point cavity rejected by new_coupling")
         if (allocated(error)) return
         call component%prepare_energy(cavity, coupling, err)
         call check(error,.not. allocated(err), "point cavity rejected by prepare_energy")
         if (allocated(error)) return
         call set_host_potential(coupling, [0.1_wp])
         energy = 0.0_wp
         call component%get_energy(component_view(coupling), cavity, energy, err)
         call check(error,.not. allocated(err), "point cavity energy failed")
         if (allocated(error)) return
         call check(error, energy == energy .and. abs(energy) < huge(1.0_wp) .and. energy /= 0.0_wp, &
            & "point cavity energy is not finite and nonzero")
         if (allocated(error)) return
      end do
   end subroutine test_conductor_point_cavity

!> Test that get_energy() reports an error when update() has not been called
   subroutine test_cpcm_requires_update(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(model_continuum_component_cpcm) :: pcm_model
      type(coupling_type) :: coupling
      !> Never updated; only its (empty) grid data reach the component
      type(cavity_type_iswig) :: cavity
      real(wp) :: energy_array

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      ! No coupling is built either: the missing matrix must be reported first
      call new_component_cpcm(pcm_model, 78.4_wp, error=err, ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if

      energy_array = 0.0_wp
      call pcm_model%get_energy(component_view(coupling), cavity, energy_array, err)
      if (.not. allocated(err)) return

      call test_failed(error, "Energy before update was correctly rejected: "//err%message)

   end subroutine test_cpcm_requires_update

!> Test that get_energy() is rejected when no molecular potential was supplied:
!> neither on a coupling that never declared the request, nor on one staged
!> for the energy phase whose potential request was left unanswered
   subroutine test_cpcm_external_potential_requires_input(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(model_continuum_component_cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      !> Never built: declares no potential request at all
      type(coupling_type) :: coupling
      !> Built and staged, but its potential request is never answered
      type(coupling_type) :: staged
      type(radius_type_static) :: radius_model
      real(wp) :: energy_array

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      call new_component_cpcm(pcm_model, epsilon=78.4_wp, error=err, &
                              param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_model%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if

      energy_array = 0.0_wp
      call pcm_model%get_energy(component_view(coupling), cavity, energy_array, err)
      if (.not. allocated(err)) return
      deallocate (err)

      call pcm_model%new_coupling(cavity, staged, err)
      if (allocated(err)) then
         call test_failed(error, "Coupling setup failed: "//err%message)
         return
      end if
      call pcm_model%prepare_energy(cavity, staged, err)
      if (allocated(err)) then
         call test_failed(error, "Energy staging failed: "//err%message)
         return
      end if
      energy_array = 0.0_wp
      call pcm_model%get_energy(component_view(staged), cavity, energy_array, err)
      if (.not. allocated(err)) return

      call test_failed(error, "Missing molecular potential was correctly rejected: "//err%message)

   end subroutine test_cpcm_external_potential_requires_input

!> Test that invalid solver identifiers are rejected during the solve step
   subroutine test_cpcm_invalid_solver(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(model_continuum_component_cpcm) :: pcm_model
      type(moist_context_type), target :: ctx

      call new_context(ctx)
      call new_component_cpcm(pcm_model, epsilon=78.4_wp, error=err, &
                              param=moist_pcm_parameters_type(solver=-1), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "Invalid solver correctly rejected: "//err%message)
      end if
   end subroutine test_cpcm_invalid_solver

!> Test that charges, host potentials and host derivatives from the wrong grid
!> are refused by name rather than read
   subroutine test_cpcm_grid_mismatch_guards(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err

      !> Molecular structure
      type(structure_type) :: mol
      !> Component under test
      type(model_continuum_component_cpcm) :: pcm_model
      !> Two cavities of the same molecule on different Lebedev grids
      type(cavity_type_iswig) :: cavity_small, cavity_large
      !> One potential trace per cavity grid
      type(coupling_type) :: coupling_small, coupling_large
      !> Radius model storage of each cavity
      type(radius_type_static) :: radius_small, radius_large
      !> Surface adjoints offered to the refused host contraction
      type(cavity_surface_adjoint_type) :: acc
      !> Potential adjoint requested before any solve
      real(wp), allocatable :: w_phi(:)
      !> Atomic point charges
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      !> Energy accumulator
      real(wp) :: energy

      !> Run context owned here and borrowed by the cavities and component
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 14, ctx, radius_small, cavity_small, err)
      if (allocated(err)) then
         call test_failed(error, "Small cavity setup failed: "//err%message)
         return
      end if
      call build_test_cavity(mol, 26, ctx, radius_large, cavity_large, err)
      if (allocated(err)) then
         call test_failed(error, "Large cavity setup failed: "//err%message)
         return
      end if
      call check(error, cavity_small%ngrid /= cavity_large%ngrid, &
         & more="both cavities have the same grid size, the test is vacuous")
      if (allocated(error)) return

      call new_component_cpcm(pcm_model, epsilon=78.4_wp, error=err, &
                              param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if

      ! No charges have been solved yet
      call pcm_model%potential_adjoint(w_phi, err)
      call check_refusal("surface charges are unavailable")
      if (allocated(error)) return

      ! Matrix on the large grid, host potential answered on the small one
      call pcm_model%update(mol, cavity_large, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if
      call stage_point_charge_energy(error, pcm_model, cavity_small, qat_vals, mol, coupling_small)
      if (allocated(error)) return
      energy = 0.0_wp
      call pcm_model%get_energy(component_view(coupling_small), cavity_large, energy, err)
      call check_refusal("potential size mismatch")
      if (allocated(error)) return
      call check(error, energy, 0.0_wp, thr=0.0_wp, &
         & more="refused energy was accumulated")
      if (allocated(error)) return

      ! Valid charges on the large grid, host derivatives offered for another size
      call stage_point_charge_energy(error, pcm_model, cavity_large, qat_vals, mol, coupling_large)
      if (allocated(error)) return
      call pcm_model%get_energy(component_view(coupling_large), cavity_large, energy, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM energy failed: "//err%message)
         return
      end if
      call pcm_model%prepare_gradient(cavity_large, coupling_large, err)
      if (allocated(err)) then
         call test_failed(error, "Gradient staging failed: "//err%message)
         return
      end if
      call fill_missing_with_zeros(cavity_large, coupling_large)
      call fill_point_charge_field(cavity_large, coupling_large, qat_vals, mol)
      call pcm_model%get_host_surface_weights(component_view(coupling_large), acc, &
         & cavity_large%ngrid + 1, err)
      call check_refusal("PCM spatial derivative grid mismatch")

   contains

      !> Require a library error naming the expected refusal
      !>
      !> @param[in] expected Distinctive substring of the expected message
      subroutine check_refusal(expected)
         !> Distinctive substring of the expected message
         character(len=*), intent(in) :: expected

         if (.not. allocated(err)) then
            call test_failed(error, "expected refusal was not raised: "//expected)
            return
         end if
         call check(error, index(err%message, expected) > 0, &
            & more="unexpected error message: "//err%message)
         deallocate (err)
      end subroutine check_refusal

   end subroutine test_cpcm_grid_mismatch_guards

   !> The tangent needs current charges and nothing else, in every model
   !>
   !> - After `update` the operators are new and no charges exist yet: refused
   !> - Directly after the energy, with no surface weights requested: accepted
   !> - After a second `update` the charges belong to the previous operators: refused
   subroutine test_tangent_requires_solution(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the components
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Solute and synthetic cavity
      type(structure_type) :: mol
      type(fd_cavity) :: cavity
      !> Zero direction and tangent accumulator
      type(cavity_surface_adjoint_type) :: direction, tangent
      !> Potential direction, adjoint tangent and energy
      real(wp), allocatable :: dphi(:), dw_phi(:)
      real(wp) :: energy
      !> Model index
      integer :: kind

      call new_context(ctx, verbosity=0)
      call get_solute(mol_triatomic, mol)
      call get_seed(seed_iswig, .true., mol, cavity, error)
      if (allocated(error)) return
      allocate (dphi(cavity%ngrid), source=0.0_wp)
      do kind = 1, npcm
         call new_pcm(kind, ctx, pcm, error)
         if (allocated(error)) return
         call pcm%update(mol, cavity, err)
         call check_moist_error(error, err, "update")
         if (allocated(error)) return
         call tangent%init(cavity%ngrid)
         call pcm%surface_weight_tangent(cavity, direction, dphi, tangent, dw_phi, err)
         call check(error, allocated(err), &
            & message=trim(pcm%name)//": tangent accepted operators without charges")
         if (allocated(error)) return
         call host_energy(pcm, mol, cavity, mol, host_charges(mol), energy, error)
         if (allocated(error)) return
         call pcm%surface_weight_tangent(cavity, direction, dphi, tangent, dw_phi, err)
         call check_moist_error(error, err, trim(pcm%name)//": tangent directly after the energy")
         if (allocated(error)) return
         call pcm%update(mol, cavity, err)
         call check_moist_error(error, err, "update")
         if (allocated(error)) return
         call pcm%surface_weight_tangent(cavity, direction, dphi, tangent, dw_phi, err)
         call check(error, allocated(err), &
            & message=trim(pcm%name)//": tangent accepted charges of stale operators")
         if (allocated(error)) return
         deallocate (err)
      end do
   end subroutine test_tangent_requires_solution

   !> The parameter derivatives are in the response only on request
   !>
   !> Default off: the response of every model holds the potential adjoint and
   !> no parameter item; with `parameter_derivatives` the item appears next to it
   !> in `get_response`, never in a read of the gradient phase
   subroutine test_parameter_derivatives_opt_in(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the components
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Solute and synthetic cavity
      type(structure_type) :: mol
      type(fd_cavity) :: cavity
      !> Host coupling and response
      type(coupling_type), target :: coupling
      type(response_type) :: response
      real(wp) :: energy
      integer :: kind, mode, nitem
      logical :: found

      call new_context(ctx, verbosity=0)
      call get_solute(mol_triatomic, mol)
      call get_seed(seed_iswig, .true., mol, cavity, error)
      if (allocated(error)) return
      do kind = 1, npcm
         do mode = 0, 1
            call new_pcm(kind, ctx, pcm, error)
            if (allocated(error)) return
            call check(error, .not. pcm%param%parameter_derivatives, &
               & message=trim(pcm%name)//": parameter derivatives off by default")
            if (allocated(error)) return
            pcm%param%parameter_derivatives = mode == 1
            call stage_energy(pcm, mol, cavity, coupling, error)
            if (allocated(error)) return
            call answer_sources(cavity, coupling, mol, host_charges(mol), real(mol%num(mol%id), wp))
            energy = 0.0_wp
            call pcm%get_energy(component_view(coupling), cavity, energy, err)
            call check_moist_error(error, err, trim(pcm%name)//": energy")
            if (allocated(error)) return
            call pcm%prepare_response(cavity, coupling, err)
            call check_moist_error(error, err, trim(pcm%name)//": response staging")
            if (allocated(error)) return
            call answer_sources(cavity, coupling, mol, host_charges(mol), real(mol%num(mol%id), wp))
            call response_clear(response)
            call pcm%get_response(component_view(coupling), cavity, response, err)
            call check_moist_error(error, err, trim(pcm%name)//": response")
            if (allocated(error)) return
            found = .false.
            nitem = 0
            do while (response%next())
               nitem = nitem + 1
               select type (item => response%item())
               type is (parameter_adjoint_response_type)
                  found = .true.
               end select
            end do
            call check(error, found .eqv. (mode == 1), &
               & message=trim(pcm%name)//": parameter item only on request "//merge("on ", "off", mode == 1))
            if (allocated(error)) return
            call check(error, nitem, 1 + mode, &
               & message=trim(pcm%name)//": response items")
            if (allocated(error)) return
            ! gradient phase: the reverse path reads the host part through get_response
            call pcm%prepare_gradient(cavity, coupling, err)
            call check_moist_error(error, err, trim(pcm%name)//": gradient staging")
            if (allocated(error)) return
            call answer_sources(cavity, coupling, mol, host_charges(mol), real(mol%num(mol%id), wp))
            call response_clear(response)
            call pcm%get_response(component_view(coupling), cavity, response, err)
            call check_moist_error(error, err, trim(pcm%name)//": gradient-phase response")
            if (allocated(error)) return
            nitem = 0
            do while (response%next())
               nitem = nitem + 1
               select type (item => response%item())
               type is (parameter_adjoint_response_type)
                  call check(error, .false., message=trim(pcm%name)//": parameter item in a gradient read")
                  if (allocated(error)) return
               end select
            end do
            call check(error, nitem, 1, message=trim(pcm%name)//": gradient-phase response items")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_parameter_derivatives_opt_in

   !> A failed gradient leaves the caller's accumulator untouched, in every model
   !>
   !> The synthetic cavity has no reverse-mode surface gradient, so the contraction
   !> fails after the direct nuclear term is known
   subroutine test_gradient_failure_preserves_accumulator(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the components
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Solute and synthetic cavity
      type(structure_type) :: mol
      type(fd_cavity) :: cavity
      !> Host coupling and the host part of the gradient phase
      type(coupling_type), target :: coupling
      type(response_type) :: response
      !> Accumulator and its value on entry
      real(wp), allocatable :: gradient(:, :), sentinel(:, :)
      !> Energy
      real(wp) :: energy
      !> Model, atom and axis indices
      integer :: kind, iat, ix

      call new_context(ctx, verbosity=0)
      call get_solute(mol_triatomic, mol)
      call get_seed(seed_iswig, .true., mol, cavity, error)
      if (allocated(error)) return
      allocate (sentinel(3, mol%nat))
      do iat = 1, mol%nat
         do ix = 1, 3
            sentinel(ix, iat) = real(10*iat + ix, wp)
         end do
      end do
      do kind = 1, npcm
         call new_pcm(kind, ctx, pcm, error)
         if (allocated(error)) return
         call stage_energy(pcm, mol, cavity, coupling, error)
         if (allocated(error)) return
         call answer_sources(cavity, coupling, mol, host_charges(mol), real(mol%num(mol%id), wp))
         energy = 0.0_wp
         call pcm%get_energy(component_view(coupling), cavity, energy, err)
         call check_moist_error(error, err, trim(pcm%name)//": energy")
         if (allocated(error)) return
         call pcm%prepare_gradient(cavity, coupling, err)
         call check_moist_error(error, err, trim(pcm%name)//": gradient staging")
         if (allocated(error)) return
         call answer_sources(cavity, coupling, mol, host_charges(mol), real(mol%num(mol%id), wp))
         gradient = sentinel
         call pcm%get_gradient(component_view(coupling), cavity, response, gradient, err)
         call check(error, allocated(err), &
            & message=trim(pcm%name)//": gradient on a cavity without surface gradient")
         if (allocated(error)) return
         deallocate (err)
         call check(error, all(gradient == sentinel), &
            & message=trim(pcm%name)//": failed gradient changed the accumulator")
         if (allocated(error)) return
      end do
   end subroutine test_gradient_failure_preserves_accumulator

   !* ================================================================================= *!
   !*                                       State                                       *!
   !* ================================================================================= *!

!> Test that a PCM instance can be reused safely when the cavity grid size changes
   subroutine test_cpcm_reallocate_on_ngrid_change(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(model_continuum_component_cpcm) :: pcm_reused, pcm_fresh
      type(cavity_type_iswig) :: cavity_small, cavity_large
      !> One potential trace per cavity grid
      type(coupling_type) :: coupling_small, coupling_large
      type(radius_type_static) :: radius_small, radius_large
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_reused, energy_fresh, charge_rms

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 14, ctx, radius_small, cavity_small, err)
      if (allocated(err)) then
         call test_failed(error, "Small cavity setup failed: "//err%message)
         return
      end if
      call build_test_cavity(mol, 50, ctx, radius_large, cavity_large, err)
      if (allocated(err)) then
         call test_failed(error, "Large cavity setup failed: "//err%message)
         return
      end if

      call new_component_cpcm(pcm_reused, epsilon=78.4_wp, error=err, &
                              param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "Reused CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_reused%update(mol, cavity_small, err)
      if (allocated(err)) then
         call test_failed(error, "Update with small cavity failed: "//err%message)
         return
      end if
      call stage_point_charge_energy(error, pcm_reused, cavity_small, qat_vals, mol, coupling_small)
      if (allocated(error)) return
      energy_reused = 0.0_wp
      call pcm_reused%get_energy(component_view(coupling_small), cavity_small, energy_reused, err)
      if (allocated(err)) then
         call test_failed(error, "Energy with small cavity failed: "//err%message)
         return
      end if

      call pcm_reused%update(mol, cavity_large, err)
      if (allocated(err)) then
         call test_failed(error, "Update with large cavity failed after reuse: "//err%message)
         return
      end if
      call stage_point_charge_energy(error, pcm_reused, cavity_large, qat_vals, mol, coupling_large)
      if (allocated(error)) return
      energy_reused = 0.0_wp
      call pcm_reused%get_energy(component_view(coupling_large), cavity_large, energy_reused, err)
      if (allocated(err)) then
         call test_failed(error, "Energy with reused large cavity failed: "//err%message)
         return
      end if

      if (size(pcm_reused%q) /= cavity_large%ngrid) then
         call test_failed(error, "Surface charge array was not resized to the new grid")
         return
      end if

      call new_component_cpcm(pcm_fresh, epsilon=78.4_wp, error=err, &
                              param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "Fresh CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_fresh%update(mol, cavity_large, err)
      if (allocated(err)) then
         call test_failed(error, "Fresh update with large cavity failed: "//err%message)
         return
      end if
      energy_fresh = 0.0_wp
      call pcm_fresh%get_energy(component_view(coupling_large), cavity_large, energy_fresh, err)
      if (allocated(err)) then
         call test_failed(error, "Fresh energy with large cavity failed: "//err%message)
         return
      end if

      call check(error, energy_reused, energy_fresh, thr=thr, &
         & message="Reused CPCM model changed energy after grid-size change")
      if (allocated(error)) return

      charge_rms = sqrt(sum((pcm_reused%q - pcm_fresh%q)**2)/real(size(pcm_fresh%q), wp))
      call check(error, charge_rms, 0.0_wp, thr=thr2, &
         & message="Reused CPCM model changed surface charges after grid-size change")
      if (allocated(error)) return

   end subroutine test_cpcm_reallocate_on_ngrid_change

!> Regression: changing the molecular potential without an intervening update()
!> must re-solve for the surface charges
!>
!> This is the ordinary SCF pattern - the geometry is fixed, so the host calls
!> update() once and then answers the potential request afresh on every
!> cycle. [[pcm_ensure_charges]] used to take the fresh phi but leave
!> `charges_valid` untouched, so the second cycle paired a fresh phi with the
!> stale q and the energy came out linear in phi instead of quadratic. The
!> potential here is the point-charge trace of `qat_vals`, so scaling the
!> charges scales phi uniformly
!>
!> The caching half is tested too: an unchanged potential answer
!> must still short-circuit before the linear solve, otherwise every accessor
!> call - and the general model issues two per get_response - pays for another
!> solve
   subroutine test_cpcm_stale_charge_regression(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(model_continuum_component_cpcm) :: pcm_model, pcm_fresh
      type(cavity_type_iswig) :: cavity
      type(radius_type_static) :: radius_model
      type(coupling_type) :: coupling
      real(wp), parameter :: epsilon = 78.4_wp
      !> Uniform factor applied to the charges between the two "SCF cycles"
      real(wp), parameter :: scale = 3.0_wp
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_first, energy_repeat, energy_cached, energy_scaled
      real(wp) :: energy_reference

      !> Run context owned here and borrowed by the cavity and components
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      call new_component_cpcm(pcm_model, epsilon=epsilon, error=err, &
                              param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if
      call pcm_model%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if

      ! Cycle 1
      call stage_point_charge_energy(error, pcm_model, cavity, qat_vals, mol, coupling)
      if (allocated(error)) return
      energy_first = 0.0_wp
      call pcm_model%get_energy(component_view(coupling), cavity, energy_first, err)
      if (allocated(err)) then
         call test_failed(error, "First cycle energy failed: "//err%message)
         return
      end if

      ! Repeating the accessor with unchanged charges must return the very same
      ! number, bit for bit
      energy_repeat = 0.0_wp
      call pcm_model%get_energy(component_view(coupling), cavity, energy_repeat, err)
      if (allocated(err)) then
         call test_failed(error, "Repeated cycle energy failed: "//err%message)
         return
      end if
      call check(error, energy_repeat, energy_first, thr=0.0_wp, &
         & message="repeating an accessor with unchanged charges changed the energy")
      if (allocated(error)) return

      ! ... and it must do so without re-solving. Poisoning the cached charges
      ! and solve energy makes the solve count observable: a short circuit keeps
      ! the poisoned values (energy doubles), a fresh solve would overwrite them
      ! and hand back the unpoisoned energy
      pcm_model%q(:) = 2.0_wp*pcm_model%q(:)
      pcm_model%kernel_energy = 2.0_wp*pcm_model%kernel_energy
      energy_cached = 0.0_wp
      call pcm_model%get_energy(component_view(coupling), cavity, energy_cached, err)
      if (allocated(err)) then
         call test_failed(error, "Cached cycle energy failed: "//err%message)
         return
      end if
      call check(error, energy_cached, 2.0_wp*energy_first, thr=thr, &
         & message="unchanged molecular potential triggered a fresh PCM solve")
      if (allocated(error)) return
      pcm_model%q(:) = 0.5_wp*pcm_model%q(:)
      pcm_model%kernel_energy = 0.5_wp*pcm_model%kernel_energy

      ! Cycle 2: new potential, same geometry, deliberately no update() in between
      call pcm_model%prepare_energy(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Second cycle staging failed: "//err%message)
         return
      end if
      call fill_point_charge_potential(cavity, coupling, scale*qat_vals, mol)
      energy_scaled = 0.0_wp
      call pcm_model%get_energy(component_view(coupling), cavity, energy_scaled, err)
      if (allocated(err)) then
         call test_failed(error, "Second cycle energy failed: "//err%message)
         return
      end if

      ! E = -0.5 f(eps) phi^T A^-1 phi is quadratic in the potential at fixed geometry
      call check(error, energy_scaled, scale**2*energy_first, thr=thr2, &
         & message="PCM energy is not quadratic in the molecular potential - "// &
         & "surface charges were not re-solved after the potential changed")
      if (allocated(error)) return

      ! Same statement against an independently constructed model
      call new_component_cpcm(pcm_fresh, epsilon=epsilon, error=err, &
                              param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "Fresh CPCM initialization failed: "//err%message)
         return
      end if
      call pcm_fresh%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Fresh CPCM update failed: "//err%message)
         return
      end if
      energy_reference = 0.0_wp
      call pcm_fresh%get_energy(component_view(coupling), cavity, energy_reference, err)
      if (allocated(err)) then
         call test_failed(error, "Fresh CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy_scaled, energy_reference, thr=thr, &
         & message="reused CPCM model does not match a freshly built one "// &
         & "after the molecular potential changed")
      if (allocated(error)) return

      call check(error, maxval(abs(pcm_model%q - pcm_fresh%q)), 0.0_wp, thr=thr, &
         & message="reused CPCM model kept stale surface charges")
      if (allocated(error)) return

   end subroutine test_cpcm_stale_charge_regression

!> Test that the component records its work on the run context's shared timer,
!> nested under whatever the caller had open. This is the observable proof that
!> the component borrows the same context as the cavity: without it the PCM
!> nodes would either be absent or sit at top level
   subroutine test_cpcm_timer_tree(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(model_continuum_component_cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(radius_type_static) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_array
      integer :: i
      logical :: found_setup, found_energy, found_solve
      character(len=:), allocatable :: name

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=2)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      call new_component_cpcm(pcm_model, epsilon=78.4_wp, error=err, &
                              param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if

      ! Open an enclosing node so the PCM nodes have to nest below it
      call ctx%timer%start("host")

      call pcm_model%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if
      call stage_point_charge_energy(error, pcm_model, cavity, qat_vals, mol, coupling)
      if (allocated(error)) return

      energy_array = 0.0_wp
      call pcm_model%get_energy(component_view(coupling), cavity, energy_array, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM energy failed: "//err%message)
         return
      end if

      call ctx%timer%stop("host")

      found_setup = .false.
      found_energy = .false.
      found_solve = .false.
      do i = 1, ctx%timer%num_nodes()
         name = ctx%timer%node_name(i)
         if (name == "PCM setup") then
            found_setup = .true.
            call check(error, ctx%timer%node_depth(i), 1, &
               & message="PCM setup timer is not nested under the enclosing node")
            if (allocated(error)) return
         else if (name == "PCM energy") then
            found_energy = .true.
            call check(error, ctx%timer%node_depth(i), 1, &
               & message="PCM energy timer is not nested under the enclosing node")
            if (allocated(error)) return
         else if (name == "PCM solve") then
            found_solve = .true.
            call check(error, ctx%timer%node_depth(i), 2, &
               & message="PCM solve timer is not nested under PCM energy")
            if (allocated(error)) return
         end if
      end do

      call check(error, found_setup, message="PCM setup timer node missing")
      if (allocated(error)) return
      call check(error, found_energy, message="PCM energy timer node missing")
      if (allocated(error)) return
      call check(error, found_solve, message="PCM solve timer node missing")
      if (allocated(error)) return

   end subroutine test_cpcm_timer_tree

   !* ================================================================================= *!
   !*                              Solvers and potential input                          *!
   !* ================================================================================= *!

!> Test CPCM energy calculation with all solvers on neutral system
   subroutine test_cpcm_energy_lu(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), parameter :: qat(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp), parameter :: epsilon = 78.4_wp
      real(wp), parameter :: ref_energy = -1.1251382967601465E-2_wp

      call get_structure(mol, "MB16-43", "01")
      call test_all_solvers(error, mol, qat, epsilon, ref_energy, "neutral system")

   end subroutine test_cpcm_energy_lu

!> Test CPCM with charged system (non-zero total charge) using all solvers
   subroutine test_cpcm_charged(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), parameter :: qat(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -1.1_wp, 0.1_wp, -0.1_wp]
      real(wp), parameter :: epsilon = 78.4_wp
      real(wp), parameter :: ref_energy = -8.9947197526553738E-2_wp

      call get_structure(mol, "MB16-43", "01")
      call test_all_solvers(error, mol, qat, epsilon, ref_energy, "charged system")

   end subroutine test_cpcm_charged

!> Test all solvers and compare to reference energy
   subroutine test_all_solvers(error, mol, qat, epsilon, ref_energy, system_name)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Molecular structure
      type(structure_type), intent(in) :: mol

      !> Atomic charges
      real(wp), intent(in) :: qat(:)

      !> Dielectric constant
      real(wp), intent(in) :: epsilon

      !> Reference energy for comparison
      real(wp), intent(in) :: ref_energy

      !> System name for error messages
      character(len=*), intent(in) :: system_name

      type(model_continuum_component_cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      real(wp) :: energy_array
      type(radius_type_static) :: radius_model
      integer :: solvers(2)
      character(len=20) :: solver_names(2)
      real(wp) :: energies(2)
      integer :: i
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      ! Setup solver types and names
      solvers = [solver_type%dense, solver_type%iterative]
      solver_names = ["dense     ", "iterative "]

      ! Build cavity and point-charge potential trace
      call new_cosmo_radii(radius_model)
      call new_cavity_iswig(cavity, radius_model=radius_model, error=err, &
                            param=moist_cavity_iswig_parameters_type(num_leb=50), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "Cavity initialization failed: "//err%message)
         return
      end if
      call cavity%update(mol, error=err)
      if (allocated(err)) then
         call test_failed(error, "Cavity update failed: "//err%message)
         return
      end if

      ! Test both generated solver choices
      do i = 1, 2
         call new_component_cpcm(pcm_model, epsilon=epsilon, error=err, &
                                 param=moist_pcm_parameters_type(solver=solvers(i)), ctx=ctx)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed for "// &
                             trim(solver_names(i))//" solver")
            return
         end if

         ! Set CG tolerance
         pcm_model%param%solver_tol = 1.0e-14_wp
         pcm_model%param%solver_maxiter = 10000

         call pcm_model%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed for "//trim(solver_names(i))//": "//err%message)
            return
         end if
         call stage_point_charge_energy(error, pcm_model, cavity, qat, mol, coupling)
         if (allocated(error)) return
         energy_array = 0.0_wp
         call pcm_model%get_energy(component_view(coupling), cavity, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed for "//trim(solver_names(i))//": "//err%message)
            return
         end if
         energies(i) = energy_array

         ! Basic sanity checks
         if (abs(energies(i)) < 1.0e-12_wp) then
            call test_failed(error, trim(solver_names(i))// &
                             " solver: Energy unexpectedly small ("//system_name//")")
            return
         end if

         if (energies(i) > 0.0_wp) then
            call test_failed(error, trim(solver_names(i))// &
                             " solver: Solvation energy should be negative ("//system_name//")")
            return
         end if

         ! Compare to reference energy
         call check(error, energies(i), ref_energy, thr=thr*10.0_wp, &
                    message=trim(solver_names(i))//" solver energy deviates from reference ("// &
                    system_name//")")
         if (allocated(error)) return
      end do

   end subroutine test_all_solvers

!> Test that all solvers give identical results
   subroutine test_cpcm_solver_comparison(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(model_continuum_component_cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(radius_type_static) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp), parameter :: epsilon = 78.4_wp
      integer :: solvers(2)
      character(len=20) :: solver_names(2)
      real(wp) :: energies(2)
      real(wp), allocatable :: charges(:, :)
      real(wp) :: energy_array
      real(wp) :: charge_rms
      integer :: i
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      ! Setup solver types and names
      solvers = [solver_type%dense, solver_type%iterative]
      solver_names = ["dense     ", "iterative "]

      ! Get test molecule
      call get_structure(mol, "MB16-43", "01")

      ! Build cavity
      call new_cosmo_radii(radius_model)
      call new_cavity_iswig(cavity, radius_model=radius_model, error=err, &
                            param=moist_cavity_iswig_parameters_type(num_leb=50), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "Cavity initialization failed: "//err%message)
         return
      end if
      call cavity%update(mol, error=err)
      if (allocated(err)) then
         call test_failed(error, "Cavity update failed: "//err%message)
         return
      end if

      ! Allocate storage for all solver results
      allocate (charges(cavity%ngrid, 2))

      ! Compare the dense and iterative generated solvers
      do i = 1, 2
         call new_component_cpcm(pcm_model, epsilon=epsilon, error=err, &
                                 param=moist_pcm_parameters_type(solver=solvers(i)), ctx=ctx)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed ("//trim(solver_names(i))//")")
            return
         end if

         ! Tighten CG tolerance
         pcm_model%param%solver_tol = 1.0e-14_wp
         pcm_model%param%solver_maxiter = 10000

         call pcm_model%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed ("//trim(solver_names(i))//"): "//err%message)
            return
         end if
         ! One coupling per component, staged with the point-charge trace
         call stage_point_charge_energy(error, pcm_model, cavity, qat_vals, mol, coupling)
         if (allocated(error)) return
         energy_array = 0.0_wp
         call pcm_model%get_energy(component_view(coupling), cavity, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed ("//trim(solver_names(i))//"): "//err%message)
            return
         end if
         energies(i) = energy_array
         charges(:, i) = pcm_model%q
      end do

      ! Compare with the dense solver
      do i = 2, 2
         call check(error, energies(i), energies(1), thr=thr, &
                    message=trim(solver_names(i))//" energy differs from "// &
                    trim(solver_names(1)))
         if (allocated(error)) return

         ! Compare charges to the dense solver
         charge_rms = sqrt(sum((charges(:, i) - charges(:, 1))**2)/real(cavity%ngrid, wp))
         call check(error, charge_rms, 0.0_wp, thr=thr2, &
                    message=trim(solver_names(i))//" charges differ from "// &
                    trim(solver_names(1)))
         if (allocated(error)) return
      end do

      deallocate (charges)

   end subroutine test_cpcm_solver_comparison

!> Test that a potential array set directly on the potential request reproduces
!> the result obtained from the point-charge fixture
   subroutine test_cpcm_external_potential(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(model_continuum_component_cpcm) :: pcm_internal, pcm_external
      type(cavity_type_iswig) :: cavity
      !> Coupling carrying the point-charge trace, and one whose potential
      !> request is answered directly with the array read back from it
      type(coupling_type) :: coupling, host_coupling
      type(radius_type_static) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_internal, energy_external, charge_rms
      real(wp), allocatable :: phi_ref(:), q_ref(:)

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      call new_component_cpcm(pcm_internal, epsilon=78.4_wp, error=err, &
                              param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_internal%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM update failed: "//err%message)
         return
      end if
      call stage_point_charge_energy(error, pcm_internal, cavity, qat_vals, mol, coupling)
      if (allocated(error)) return

      energy_internal = 0.0_wp
      call pcm_internal%get_energy(component_view(coupling), cavity, energy_internal, err)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM energy failed: "//err%message)
         return
      end if

      allocate (phi_ref, source=pcm_internal%phi)
      allocate (q_ref, source=pcm_internal%q)

      call new_component_cpcm(pcm_external, epsilon=78.4_wp, error=err, &
                              param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "External-potential CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_external%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "External-potential CPCM update failed: "//err%message)
         return
      end if

      call pcm_external%new_coupling(cavity, host_coupling, err)
      if (allocated(err)) then
         call test_failed(error, "External-potential coupling setup failed: "//err%message)
         return
      end if
      call pcm_external%prepare_energy(cavity, host_coupling, err)
      if (allocated(err)) then
         call test_failed(error, "External-potential energy staging failed: "//err%message)
         return
      end if
      call set_host_potential(host_coupling, phi_ref)
      energy_external = 0.0_wp
      call pcm_external%get_energy(component_view(host_coupling), cavity, energy_external, err)
      if (allocated(err)) then
         call test_failed(error, "External-potential CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy_external, energy_internal, thr=thr, &
         & message="directly set potential energy differs from the point-charge fixture")
      if (allocated(error)) return

      charge_rms = sqrt(sum((pcm_external%q - q_ref)**2)/real(size(q_ref), wp))
      call check(error, charge_rms, 0.0_wp, thr=thr2, &
         & message="External potential surface charges differ from internal reference")
      if (allocated(error)) return

   end subroutine test_cpcm_external_potential


   !> SLIC reports a NaN normal trace on a point surface instead of a bad energy
   subroutine test_slic_rejects_nan_normal_trace_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(error_type), allocatable :: setup
      type(moist_context_type), target :: ctx
      class(model_continuum_component_pcm), allocatable :: pcm
      type(structure_type) :: mol
      type(fd_cavity), target :: cavity
      type(coupling_type), target :: coupling
      real(wp) :: normal(3), energy

      call new_context(ctx, verbosity=0)
      call stage_poisoned_trace(ctx, .false., ieee_value(1.0_wp, ieee_quiet_nan), mol, cavity, pcm, &
         & coupling, normal, setup)
      if (allocated(setup)) return
      energy = 0.0_wp
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err)
   end subroutine test_slic_rejects_nan_normal_trace_point

   !> SLIC reports a +Inf normal trace on a point surface instead of a bad energy
   subroutine test_slic_rejects_inf_normal_trace_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(error_type), allocatable :: setup
      type(moist_context_type), target :: ctx
      class(model_continuum_component_pcm), allocatable :: pcm
      type(structure_type) :: mol
      type(fd_cavity), target :: cavity
      type(coupling_type), target :: coupling
      real(wp) :: normal(3), energy

      call new_context(ctx, verbosity=0)
      call stage_poisoned_trace(ctx, .false., ieee_value(1.0_wp, ieee_positive_inf), mol, cavity, &
         & pcm, coupling, normal, setup)
      if (allocated(setup)) return
      energy = 0.0_wp
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err)
   end subroutine test_slic_rejects_inf_normal_trace_point

   !> SLIC reports a NaN normal trace on a gaussian surface instead of a bad energy
   subroutine test_slic_rejects_nan_normal_trace_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(error_type), allocatable :: setup
      type(moist_context_type), target :: ctx
      class(model_continuum_component_pcm), allocatable :: pcm
      type(structure_type) :: mol
      type(fd_cavity), target :: cavity
      type(coupling_type), target :: coupling
      real(wp) :: normal(3), energy

      call new_context(ctx, verbosity=0)
      call stage_poisoned_trace(ctx, .true., ieee_value(1.0_wp, ieee_quiet_nan), mol, cavity, pcm, &
         & coupling, normal, setup)
      if (allocated(setup)) return
      energy = 0.0_wp
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err)
   end subroutine test_slic_rejects_nan_normal_trace_gaussian

   !> SLIC reports a +Inf normal trace on a gaussian surface instead of a bad energy
   subroutine test_slic_rejects_inf_normal_trace_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(error_type), allocatable :: setup
      type(moist_context_type), target :: ctx
      class(model_continuum_component_pcm), allocatable :: pcm
      type(structure_type) :: mol
      type(fd_cavity), target :: cavity
      type(coupling_type), target :: coupling
      real(wp) :: normal(3), energy

      call new_context(ctx, verbosity=0)
      call stage_poisoned_trace(ctx, .true., ieee_value(1.0_wp, ieee_positive_inf), mol, cavity, &
         & pcm, coupling, normal, setup)
      if (allocated(setup)) return
      energy = 0.0_wp
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err)
   end subroutine test_slic_rejects_inf_normal_trace_gaussian

   !> Stage SLIC and poison the normal of the first point
   !>
   !> The coupling refuses a non-finite host answer, so the trace is poisoned through the
   !> cavity normal; Newton then sees a non-finite forcing. Setup failures come back in
   !> `setup`, never in the error of the call under test
   !>
   !> @param[in]  ctx      run context borrowed by the component
   !> @param[in]  gaussian Gaussian instead of point surface
   !> @param[in]  bad      value written into the normal
   !> @param[out] mol      solute
   !> @param[out] cavity   cavity carrying the poisoned normal
   !> @param[out] pcm      staged component
   !> @param[out] coupling staged host coupling
   !> @param[out] normal   original normal of the first point
   !> @param[out] setup    allocated if the setup failed
   subroutine stage_poisoned_trace(ctx, gaussian, bad, mol, cavity, pcm, coupling, normal, setup)
      !> Run context
      type(moist_context_type), target, intent(in) :: ctx
      !> Surface kind
      logical, intent(in) :: gaussian
      !> Bad value
      real(wp), intent(in) :: bad
      !> Solute
      type(structure_type), intent(out) :: mol
      !> Cavity
      type(fd_cavity), target, intent(out) :: cavity
      !> Component
      class(model_continuum_component_pcm), allocatable, intent(out) :: pcm
      !> Coupling
      type(coupling_type), target, intent(out) :: coupling
      !> Original normal
      real(wp), intent(out) :: normal(3)
      !> Setup failure
      type(error_type), allocatable, intent(out) :: setup

      call get_solute(mol_triatomic, mol)
      call get_seed(seed_drop, gaussian, mol, cavity, setup)
      if (allocated(setup)) return
      call new_pcm(pcm_slic, ctx, pcm, setup)
      if (allocated(setup)) return
      call stage_energy(pcm, mol, cavity, coupling, setup)
      if (allocated(setup)) return
      call answer_sources(cavity, coupling, mol, host_charges(mol), real(mol%num(mol%id), wp))
      normal = cavity%normal0(:, 1)
      cavity%normal0(:, 1) = bad
   end subroutine stage_poisoned_trace

   !> SLIC recovers from a bad normal trace and keeps no bad energy
   !>
   !> The failed energy call leaves the output untouched, and the component solves
   !> again once the normal is restored; point and Gaussian surfaces
   subroutine test_slic_recovers_after_bad_trace(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      class(model_continuum_component_pcm), allocatable :: pcm
      type(structure_type) :: mol
      type(fd_cavity), target :: cavity
      type(coupling_type), target :: coupling
      real(wp) :: normal(3), energy
      integer :: ismear

      call new_context(ctx, verbosity=0)
      do ismear = 1, 2
         call stage_poisoned_trace(ctx, ismear == 2, ieee_value(1.0_wp, ieee_quiet_nan), mol, &
            & cavity, pcm, coupling, normal, error)
         if (allocated(error)) return
         energy = 0.0_wp
         call pcm%get_energy(component_view(coupling), cavity, energy, err)
         call check(error, ieee_is_finite(energy), "SLIC: energy set from a bad trace")
         if (allocated(error)) return
         cavity%normal0(:, 1) = normal
         call host_energy(pcm, mol, cavity, mol, host_charges(mol), energy, error)
         if (allocated(error)) return
         call check(error, ieee_is_finite(energy), "SLIC: no recovery")
         if (allocated(error)) return
      end do
      call ctx%delete()
   end subroutine test_slic_recovers_after_bad_trace

   !> A dense LPB solve of a non-finite operator is an error, never a NaN energy
   !>
   !> A non-finite curvature reaches only the operator self term, not the host traces
   subroutine test_lpb_nonfinite_operator(error)
      !> Test failure, the library error of the dense solve
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(error_type), allocatable :: setup
      type(moist_context_type), target :: ctx
      class(model_continuum_component_pcm), allocatable :: pcm
      type(structure_type) :: mol
      type(fd_cavity) :: cavity
      type(coupling_type), target :: coupling
      real(wp) :: energy

      call new_context(ctx, verbosity=0)
      call get_solute(mol_triatomic, mol)
      call get_seed(seed_drop, .true., mol, cavity, setup)
      if (allocated(setup)) return
      cavity%k1(1) = ieee_value(1.0_wp, ieee_quiet_nan)
      call new_pcm(pcm_lpb, ctx, pcm, setup)
      if (allocated(setup)) return
      call stage_energy(pcm, mol, cavity, coupling, setup)
      if (allocated(setup)) return
      call answer_sources(cavity, coupling, mol, host_charges(mol), real(mol%num(mol%id), wp))
      energy = 0.0_wp
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err)
   end subroutine test_lpb_nonfinite_operator
end module test_model_component_pcm_interface
