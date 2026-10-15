!> Constructor and parameter tests for generated layer PCM components
module test_model_component_pcm_layers
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, fatal_error
   use mctc_io, only: structure_type, new
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, moist_cavity_iswig_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_isodensity_internal, only: moist_cavity_drop_lsf_isodensity_internal_type
   use moist_radii, only: default_cpcm_radii
   use moist_model_continuum_type, only: model_continuum_type, new_continuum_model
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, new_component_cpcm, &
      & model_continuum_component_cosmo, new_component_cosmo
   use mctc_env_error, only: moist_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_model_continuum_component_pcm_type, only: moist_pcm_parameters_type, solver_type
   use moist_model_continuum_component_pcm_iefpcm, only: model_continuum_component_iefpcm, new_component_iefpcm
   use moist_model_continuum_component_pcm_ssvpe, only: model_continuum_component_ssvpe, new_component_ssvpe
   use moist_model_continuum_component_pcm_dpcm, only: model_continuum_component_dpcm, new_component_dpcm
   use moist_model_continuum_component_pcm_dpcm_original, only: model_continuum_component_dpcm_original, &
      & new_component_dpcm_original
   use moist_model_continuum_component_pcm_lpb, only: model_continuum_component_lpb, new_component_lpb, debye_kappa
   use moist_model_continuum_component_pcm_slic, only: model_continuum_component_slic, new_component_slic
   use moist_model_continuum_component_pcm_nonlocal, only: model_continuum_component_nonlocal, new_component_nonlocal
   use moist_channels_coupling, only: coupling_type, gaussian_potential_request_type, &
      & point_potential_request_type, nuclear_charge_request_type
   use moist_channels_response, only: response_type, potential_adjoint_response_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use test_helpers, only: component_view, fd4_scalar, fd4_offsets, set_host_potential
   use test_model_component_helper, only: nfixture => fixture_ngrid_param, &
      & fixture_xyz => fixture_xyz_param, fixture_normals => fixture_normals_param, &
      & fixture_areas => fixture_areas_param, fixture_xis => fixture_xis_param, fixture_fs => fixture_fs_param
   use moist_model_continuum_component_pcm_solvers, only: pcm_dense_solver, pcm_storage_matrix_free
   use pcmgen_gauss_ssvpe, only: ssvpe_state, ssvpe_assemble, &
      & ssvpe_medium, ssvpe_solve
   use pcmgen_gauss_iefpcm, only: iefpcm_state, iefpcm_assemble, &
      & iefpcm_medium, iefpcm_solve
   use test_model_component_pcm_layer_golden, only: golden_eps, hf_n, hf_xyz, hf_normal, &
      & hf_curv, hf_xi, hf_f, hf_area, hf_v, hf_q, hf_qsym, hf_energy, &
      & hf_ief_q, hf_ief_qsym, hf_ief_energy
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_layers

contains

   !> Register focused model tests
   subroutine collect_model_component_pcm_layers(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("constructors", test_constructors), &
         & new_unittest("optional-context", test_optional_context), &
         & new_unittest("cavity-requirements", test_cavity_requirements), &
         & new_unittest("automatic-cavity-geometry", test_automatic_cavity_geometry), &
         & new_unittest("invalid-physical-parameters", test_invalid_parameters), &
         & new_unittest("PySCF-layer-goldens", test_layer_goldens), &
         & new_unittest("layer-lifecycle-guards", test_layer_lifecycle_guards), &
         & new_unittest("Gaussian-component-paths", test_gaussian_component_paths), &
         & new_unittest("point-component-paths", test_point_component_paths)]
   end subroutine collect_model_component_pcm_layers

   !> All PCM models inherit missing contexts and preserve explicit ones
   !>
   !> @param[out] error Test failure
   subroutine test_optional_context(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library diagnostic
      type(moist_error), allocatable :: err
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
   !>
   !> @param[out] error Test failure
   subroutine test_cavity_requirements(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library diagnostic
      type(moist_error), allocatable :: err
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
               call check(error, cavity%param%do_curvature .eqv. (had_curvature .or. (kind >= 3 .and. kind <= 8)))
               if (allocated(error)) return
               call check(error, cavity%param%do_rho .and. .not. cavity%param%do_fine &
                  & .and. .not. cavity%param%do_grid_density .and. .not. cavity%param%do_r_iI)
               if (allocated(error)) return
               call cavity%configure(normal=.false., curvature=.false., error=err)
               call check(error,.not. allocated(err))
               if (allocated(error)) return
               call check(error, cavity%param%do_normal .eqv. (had_normal .or. kind >= 3))
               if (allocated(error)) return
               call check(error, cavity%param%do_curvature .eqv. (had_curvature .or. (kind >= 3 .and. kind <= 8)))
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
               call check(error, cavity%need_curvature .eqv. (had_curvature .or. (kind >= 3 .and. kind <= 8)))
               if (allocated(error)) return
               call cavity%configure(curvature=.false., error=err)
               call check(error,.not. allocated(err))
               if (allocated(error)) return
               call check(error, cavity%need_curvature .eqv. (had_curvature .or. (kind >= 3 .and. kind <= 8)))
            end select
            if (allocated(error)) return
            call check(error, iswig%need_curvature .eqv. had_curvature)
            if (allocated(error)) return
            deallocate (component)
         end do
      end do
   end subroutine test_cavity_requirements

   !> Model update computes required normals and curvatures without do_fine
   !>
   !> @param[out] error Test failure
   subroutine test_automatic_cavity_geometry(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library diagnostic
      type(moist_error), allocatable :: err
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

   !> Construct each PCM family with optional context forwarding
   !>
   !> @param[in] kind Constructor index
   !> @param[out] component Constructed component
   !> @param[out] err Library diagnostic
   !> @param[in] ctx Optional borrowed context
   subroutine make_component(kind, component, err, ctx)
      !> Constructor index
      integer, intent(in) :: kind
      !> Constructed component
      class(model_continuum_component_pcm), allocatable, intent(out) :: component
      !> Library diagnostic
      type(moist_error), allocatable, intent(out) :: err
      !> Optional borrowed context
      type(moist_context_type), target, optional, intent(in) :: ctx

      select case (kind)
      case (1)
         allocate (model_continuum_component_cpcm :: component)
         select type (component)
         type is (model_continuum_component_cpcm)
            call new_component_cpcm(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (2)
         allocate (model_continuum_component_cosmo :: component)
         select type (component)
         type is (model_continuum_component_cosmo)
            call new_component_cosmo(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (3)
         allocate (model_continuum_component_iefpcm :: component)
         select type (component)
         type is (model_continuum_component_iefpcm)
            call new_component_iefpcm(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (4)
         allocate (model_continuum_component_ssvpe :: component)
         select type (component)
         type is (model_continuum_component_ssvpe)
            call new_component_ssvpe(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (5)
         allocate (model_continuum_component_dpcm :: component)
         select type (component)
         type is (model_continuum_component_dpcm)
            call new_component_dpcm(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (6)
         allocate (model_continuum_component_dpcm_original :: component)
         select type (component)
         type is (model_continuum_component_dpcm_original)
            call new_component_dpcm_original(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (7)
         allocate (model_continuum_component_lpb :: component)
         select type (component)
         type is (model_continuum_component_lpb)
            call new_component_lpb(component, 80.0_wp, 0.1_wp, error=err, ctx=ctx)
         end select
      case (8)
         allocate (model_continuum_component_slic :: component)
         select type (component)
         type is (model_continuum_component_slic)
            call new_component_slic(component, 80.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, error=err, ctx=ctx)
         end select
      case (9)
         allocate (model_continuum_component_nonlocal :: component)
         select type (component)
         type is (model_continuum_component_nonlocal)
            call new_component_nonlocal(component, 80.0_wp, 2.0_wp, 1.0_wp, error=err, ctx=ctx)
         end select
      case default
         call fatal_error(err, "Unknown PCM component fixture")
      end select
   end subroutine make_component

   !> Verify each concrete type copies the shared solver controls
   subroutine test_constructors(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
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
      call check(error,.not. lpb%is_inert())
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
      call check(error,.not. nonlocal%needs_curvature)
      if (allocated(error)) return
      call check(error, debye_kappa(78.4_wp, 0.1_wp, 298.15_wp), 0.055_wp, thr=0.01_wp)
      call ctx%delete()
   end subroutine test_constructors

   !> Reject nonphysical material values before assembling a cavity
   subroutine test_invalid_parameters(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_continuum_component_iefpcm) :: ief
      type(model_continuum_component_lpb) :: lpb
      type(model_continuum_component_nonlocal) :: nonlocal

      call new_context(ctx, verbosity=0)
      call new_component_iefpcm(ief, 0.5_wp, err, ctx=ctx)
      call check(error, allocated(err))
      if (allocated(error)) return
      call new_component_lpb(lpb, 80.0_wp, -1.0_wp, err, ctx=ctx)
      call check(error, allocated(err))
      if (allocated(error)) return
      call new_component_nonlocal(nonlocal, 80.0_wp, 2.0_wp, 0.0_wp, err, ctx=ctx)
      call check(error, allocated(err))
      call ctx%delete()
   end subroutine test_invalid_parameters

   !> Generated Gaussian IEF-PCM and SS(V)PE kernels reproduce PySCF charges
   subroutine test_layer_goldens(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(ssvpe_state) :: ssv
      type(iefpcm_state) :: ief
      type(pcm_dense_solver) :: solver
      real(wp) :: q(hf_n), qsym(hf_n), energy

      call ssvpe_assemble(ssv, hf_xyz, hf_normal, hf_curv, hf_curv, hf_xi, hf_f, &
         & hf_area, err, solver)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call ssvpe_medium(ssv, golden_eps, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call ssvpe_solve(ssv, hf_v, energy, q, qsym, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, maxval(abs(q - hf_q)), 0.0_wp, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(qsym - hf_qsym)), 0.0_wp, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, energy, hf_energy, thr=1.0e-10_wp)
      if (allocated(error)) return

      call iefpcm_assemble(ief, hf_xyz, hf_normal, hf_curv, hf_curv, hf_xi, hf_f, &
         & hf_area, err, solver)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call iefpcm_medium(ief, golden_eps, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call iefpcm_solve(ief, hf_v, energy, q, qsym, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, maxval(abs(q - hf_ief_q)), 0.0_wp, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(qsym - hf_ief_qsym)), 0.0_wp, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, energy, hf_ief_energy, thr=1.0e-10_wp)
   end subroutine test_layer_goldens

   !> Refuse missing or stale surface data and unsolved adjoints
   subroutine test_layer_lifecycle_guards(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
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
      call expect_error("position derivative shape mismatch")
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

   !> Gaussian component solves, host responses and geometry adjoints
   subroutine test_gaussian_component_paths(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: kind

      do kind = 1, 9
         call check_component_path(kind, .true., error)
         if (allocated(error)) return
      end do
   end subroutine test_gaussian_component_paths

   !> Point component solves, host responses and geometry adjoints
   subroutine test_point_component_paths(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: kind

      do kind = 1, 9
         call check_component_path(kind, .false., error)
         if (allocated(error)) return
      end do
   end subroutine test_point_component_paths

   !> Check a component through its public coupling and surface interfaces
   !>
   !> @param[in] kind Component constructor index
   !> @param[in] gaussian Whether the surface has Gaussian widths
   !> @param[out] error Test failure
   subroutine check_component_path(kind, gaussian, error)
      !> Component constructor index
      integer, intent(in) :: kind
      !> Whether the surface has Gaussian widths
      logical, intent(in) :: gaussian
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(model_continuum_component_pcm), allocatable :: component
      type(cavity_type_drop) :: cavity, displaced
      type(cavity_surface_adjoint_type) :: weights, response_weights
      type(coupling_type), target :: coupling
      type(response_type) :: response
      type(potential_adjoint_response_type) :: trace
      real(wp), parameter :: h = 1.0e-4_wp
      real(wp), parameter :: direction(nfixture) = &
         & [0.02_wp, -0.03_wp, 0.04_wp, -0.01_wp, 0.05_wp, -0.02_wp, 0.03_wp]
      real(wp) :: energy, cached, iterative, samples(4), expected, numerical, shift
      real(wp) :: dxyz(3, nfixture), dn(3, nfixture), dg(3, nfixture)
      real(wp), allocatable :: dense_phi(:), dense_normal(:), iterative_normal(:)
      integer :: channel, offset, i
      character(len=:), allocatable :: label

      call new_context(ctx, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call make_component(kind, component, err, ctx)
      call require_success("constructor")
      if (allocated(error)) return
      label = trim(component%name)//merge(" Gaussian", " point   ", gaussian)
      component%param%solver_tol = 1.0e-12_wp
      component%param%newton_tol = 1.0e-13_wp
      cavity%ngrid = nfixture
      cavity%nsph = 1
      cavity%updated = .true.
      cavity%built = .true.
      cavity%xyz = fixture_xyz
      cavity%normal0 = fixture_normals
      cavity%a = fixture_areas
      allocate (cavity%k1(nfixture), source=0.4_wp)
      allocate (cavity%k2(nfixture), source=0.5_wp)
      if (gaussian) then
         cavity%xi0 = fixture_xis
         cavity%f = fixture_fs
      end if
      do i = 1, nfixture
         dxyz(:, i) = direction(i)*[1.0_wp, -0.7_wp, 0.3_wp]
         dn(:, i) = direction(i)*[-0.4_wp, 0.2_wp, 0.8_wp]
         dg(:, i) = direction(i)*[0.3_wp, 0.6_wp, -0.2_wp]
      end do

      call component%update(mol, cavity, err)
      call require_success("update")
      if (allocated(error)) return
      call sample_energy(cavity, 0.0_wp, energy)
      if (allocated(error)) return
      call check(error, ieee_is_finite(energy) .and. abs(energy) > 1.0e-10_wp, &
         & message=label//" nontrivial finite energy")
      if (allocated(error)) return
      call sample_energy(cavity, 0.0_wp, cached)
      if (allocated(error)) return
      call check(error, cached, energy, thr=1.0e-14_wp, message=label//" cached energy")
      if (allocated(error)) return

      call component%prepare_response(cavity, coupling, err)
      call require_success("prepare response")
      if (allocated(error)) return
      call answer_host(cavity, 0.0_wp)
      if (allocated(error)) return
      call component%get_response(component_view(coupling), cavity, response, err)
      call require_success("response")
      if (allocated(error)) return
      call check(error, response%next(), message=label//" potential response exists")
      if (allocated(error)) return
      select type (item => response%item())
      type is (potential_adjoint_response_type)
         trace = item
      class default
         call test_failed(error, label//" unexpected response item")
         return
      end select
      call check(error, allocated(trace%w_phi), message=label//" potential adjoint exists")
      if (allocated(error)) return
      call check(error, allocated(trace%w_dphi_dr) .eqv. component%uses_normal_trace(), &
         & message=label//" normal-trace response requirement")
      if (allocated(error)) return
      dense_phi = trace%w_phi
      call component%normal_trace_adjoint(dense_normal, err)
      call require_success("normal adjoint")
      if (allocated(error)) return
      expected = dot_product(trace%w_phi, direction)
      if (allocated(trace%w_dphi_dr)) expected = expected + sum(trace%w_dphi_dr*dg)
      do offset = 1, 4
         call sample_energy(cavity, h*fd4_offsets(offset), samples(offset))
         if (allocated(error)) return
      end do
      call fd4_scalar(samples(1), samples(2), samples(3), samples(4), h, numerical, error)
      if (allocated(error)) return
      call check(error, numerical, expected, thr=2.0e-8_wp, message=label//" host-trace adjoint")
      if (allocated(error)) return
      expected = 0.0_wp
      if (allocated(trace%w_dphi_dr)) expected = sum(trace%w_dphi_dr*dg)
      do offset = 1, 4
         call sample_energy(cavity, h*fd4_offsets(offset), samples(offset), normal_only=.true.)
         if (allocated(error)) return
      end do
      call fd4_scalar(samples(1), samples(2), samples(3), samples(4), h, numerical, error)
      if (allocated(error)) return
      call check(error, numerical, expected, thr=2.0e-8_wp, message=label//" changed normal trace")
      if (allocated(error)) return

      call sample_energy(cavity, 0.0_wp, cached)
      if (allocated(error)) return
      call component%prepare_gradient(cavity, coupling, err)
      call require_success("prepare gradient")
      if (allocated(error)) return
      call answer_host(cavity, 0.0_wp)
      if (allocated(error)) return
      call weights%init(nfixture)
      call component%get_gradient_surface_weights(component_view(coupling), cavity, weights, err)
      call require_success("surface adjoints")
      if (allocated(error)) return
      if (gaussian) then
         allocate (moist_cavity_drop_lsf_isodensity_internal_type :: cavity%lsf_model)
         call component%prepare_response(cavity, coupling, err)
         call require_success("prepare field-dependent response")
         if (allocated(error)) return
         call answer_host(cavity, 0.0_wp)
         if (allocated(error)) return
         call response_weights%init(nfixture)
         call component%get_surface_weights(component_view(coupling), cavity, response_weights, err)
         call require_success("field-dependent surface adjoints")
         if (allocated(error)) return
         call check(error, maxval(abs(response_weights%w_xyz - weights%w_xyz)), &
            & 0.0_wp, thr=1.0e-12_wp, message=label//" field-dependent position adjoints")
         if (allocated(error)) return
         call check(error, maxval(abs(response_weights%w_xi - weights%w_xi)), &
            & 0.0_wp, thr=1.0e-12_wp, message=label//" field-dependent width adjoints")
         if (allocated(error)) return
         call check(error, maxval(abs(response_weights%w_n - weights%w_n)), &
            & 0.0_wp, thr=1.0e-12_wp, message=label//" field-dependent normal adjoints")
         if (allocated(error)) return
      end if
      do channel = 1, merge(7, 5, gaussian)
         select case (channel)
         case (1)
            expected = sum(weights%w_xyz*dxyz)
         case (2)
            expected = sum(weights%w_n*dn)
         case (3)
            expected = dot_product(weights%w_a, direction)
         case (4)
            expected = dot_product(weights%w_k1, direction)
         case (5)
            expected = dot_product(weights%w_k2, direction)
         case (6)
            expected = dot_product(weights%w_xi, direction)
         case (7)
            expected = dot_product(weights%w_f, direction)
         case default
            call test_failed(error, label//" unknown geometry channel")
            return
         end select
         do offset = 1, 4
            displaced = cavity
            shift = h*fd4_offsets(offset)
            select case (channel)
            case (1)
               displaced%xyz = cavity%xyz + shift*dxyz
            case (2)
               displaced%normal0 = cavity%normal0 + shift*dn
            case (3)
               displaced%a = cavity%a + shift*direction
            case (4)
               displaced%k1 = cavity%k1 + shift*direction
            case (5)
               displaced%k2 = cavity%k2 + shift*direction
            case (6)
               displaced%xi0 = cavity%xi0 + shift*direction
            case (7)
               displaced%f = cavity%f + shift*direction
            case default
               call test_failed(error, label//" unknown geometry channel")
               return
            end select
            call component%update(mol, displaced, err)
            call require_success("displaced update")
            if (allocated(error)) return
            call sample_energy(displaced, 0.0_wp, samples(offset))
            if (allocated(error)) return
         end do
         call fd4_scalar(samples(1), samples(2), samples(3), samples(4), h, numerical, error)
         if (allocated(error)) return
         call check(error, numerical, expected, thr=2.0e-8_wp, &
            & message=label//" geometry channel "//achar(iachar("0") + channel))
         if (allocated(error)) return
      end do

      component%param%solver = solver_type%iterative
      component%param%storage = pcm_storage_matrix_free
      call component%update(mol, cavity, err)
      call require_success("iterative update")
      if (allocated(error)) return
      call sample_energy(cavity, 0.0_wp, iterative)
      if (allocated(error)) return
      call check(error, iterative, energy, thr=2.0e-9_wp, message=label//" dense/iterative energy")
      if (allocated(error)) return
      call component%potential_adjoint(trace%w_phi, err)
      call require_success("iterative adjoint")
      if (allocated(error)) return
      call check(error, maxval(abs(trace%w_phi - dense_phi)), 0.0_wp, thr=2.0e-9_wp, &
         & message=label//" dense/iterative potential adjoint")
      if (allocated(error)) return
      call component%normal_trace_adjoint(iterative_normal, err)
      call require_success("iterative normal adjoint")
      if (allocated(error)) return
      call check(error, maxval(abs(iterative_normal - dense_normal)), 0.0_wp, thr=2.0e-9_wp, &
         & message=label//" dense/iterative normal adjoint")
      if (allocated(error)) return
      call component%prepare_gradient(cavity, coupling, err)
      call require_success("prepare iterative gradient")
      if (allocated(error)) return
      call answer_host(cavity, 0.0_wp)
      if (allocated(error)) return
      call response_weights%init(nfixture)
      call component%get_gradient_surface_weights(component_view(coupling), cavity, response_weights, err)
      call require_success("iterative surface adjoints")
      if (allocated(error)) return
      call check(error, max(maxval(abs(response_weights%w_xyz - weights%w_xyz)), &
         & maxval(abs(response_weights%w_n - weights%w_n)), &
         & maxval(abs(response_weights%w_a - weights%w_a)), &
         & maxval(abs(response_weights%w_k1 - weights%w_k1)), &
         & maxval(abs(response_weights%w_k2 - weights%w_k2)), &
         & maxval(abs(response_weights%w_xi - weights%w_xi)), &
         & maxval(abs(response_weights%w_f - weights%w_f))), &
         & 0.0_wp, thr=2.0e-9_wp, message=label//" dense/iterative surface adjoints")

   contains

      !> Turn a library diagnostic into a model-labelled test failure
      !>
      !> @param[in] phase Operation that produced the diagnostic
      subroutine require_success(phase)
         !> Operation that produced the diagnostic
         character(len=*), intent(in) :: phase
         if (.not. allocated(err)) return
         if (allocated(label)) then
            call test_failed(error, label//" "//phase//": "//err%message)
         else
            call test_failed(error, "PCM component "//phase//": "//err%message)
         end if
      end subroutine require_success

      !> Solve a staged potential, preserving the current assembled operators
      !>
      !> @param[inout] surface Synthetic cavity surface
      !> @param[in] trace_shift Independent host-potential perturbation
      !> @param[out] value Component energy
      !> @param[in] normal_only Change the normal trace while keeping phi fixed
      subroutine sample_energy(surface, trace_shift, value, normal_only)
         !> Synthetic cavity surface
         type(cavity_type_drop), intent(inout) :: surface
         !> Independent host-potential perturbation
         real(wp), intent(in) :: trace_shift
         !> Component energy
         real(wp), intent(out) :: value
         !> Change the normal trace while keeping phi fixed
         logical, intent(in), optional :: normal_only

         call component%new_coupling(surface, coupling, err)
         call require_success("new coupling")
         if (allocated(error)) return
         call component%prepare_energy(surface, coupling, err)
         call require_success("prepare energy")
         if (allocated(error)) return
         call answer_host(surface, trace_shift, normal_only)
         if (allocated(error)) return
         value = 0.0_wp
         call component%get_energy(component_view(coupling), surface, value, err)
         call require_success("energy")
      end subroutine sample_energy

      !> Answer a quadratic host potential and its analytic raw derivatives
      !>
      !> @param[in] surface Synthetic cavity surface
      !> @param[in] trace_shift Independent host-potential perturbation
      !> @param[in] normal_only Change the normal trace while keeping phi fixed
      subroutine answer_host(surface, trace_shift, normal_only)
         !> Synthetic cavity surface
         type(cavity_type_drop), intent(in) :: surface
         !> Independent host-potential perturbation
         real(wp), intent(in) :: trace_shift
         !> Change the normal trace while keeping phi fixed
         logical, intent(in), optional :: normal_only
         real(wp) :: phi(nfixture), grad(3, nfixture), width(nfixture)
         real(wp) :: hessian(3, 3, nfixture), mixed(3, nfixture), dx(3), dxi, base_grad(3), phi_shift
         integer :: point, axis

         hessian = 0.0_wp
         mixed = 0.003_wp
         phi_shift = trace_shift
         if (present(normal_only)) then
            if (normal_only) phi_shift = 0.0_wp
         end if
         do point = 1, nfixture
            dx = surface%xyz(:, point) - fixture_xyz(:, point)
            dxi = 0.0_wp
            if (gaussian) dxi = surface%xi0(point) - fixture_xis(point)
            base_grad = 0.01_wp*real(point, wp)*[1.0_wp, -0.5_wp, 0.7_wp]
            do axis = 1, 3
               hessian(axis, axis, point) = 0.002_wp*real(axis + point, wp)
            end do
            width(point) = 0.013_wp*real(point, wp) + dot_product(mixed(:, point), dx)
            phi(point) = (-1.0_wp)**point*0.05_wp*real(point, wp) + &
               & dot_product(base_grad, dx) + 0.5_wp*dot_product(dx, matmul(hessian(:, :, point), dx)) + &
               & width(point)*dxi + phi_shift*direction(point)
            grad(:, point) = base_grad + matmul(hessian(:, :, point), dx) + &
               & mixed(:, point)*dxi + trace_shift*dg(:, point)
         end do
         do while (coupling%next())
            select type (request => coupling%request())
            type is (gaussian_potential_request_type)
               call coupling%answer("phi", phi, err)
               call require_success("submit phi")
               if (allocated(error)) return
               call coupling%answer("dphi_dr", grad, err)
               call require_success("submit dphi_dr")
               if (allocated(error)) return
               call coupling%answer("dphi_dxi", width, err)
               call require_success("submit dphi_dxi")
               if (allocated(error)) return
               call coupling%answer("d2phi_dr2", hessian, err)
               call require_success("submit d2phi_dr2")
               if (allocated(error)) return
               call coupling%answer("d2phi_dr_dxi", mixed, err)
               call require_success("submit d2phi_dr_dxi")
               if (allocated(error)) return
            type is (point_potential_request_type)
               call coupling%answer("phi", phi, err)
               call require_success("submit phi")
               if (allocated(error)) return
               call coupling%answer("dphi_dr", grad, err)
               call require_success("submit dphi_dr")
               if (allocated(error)) return
               call coupling%answer("d2phi_dr2", hessian, err)
               call require_success("submit d2phi_dr2")
               if (allocated(error)) return
            type is (nuclear_charge_request_type)
               call coupling%answer("z", [1.0_wp], err)
               call require_success("submit nuclear charge")
               if (allocated(error)) return
            end select
         end do
      end subroutine answer_host

   end subroutine check_component_path

end module test_model_component_pcm_layers
