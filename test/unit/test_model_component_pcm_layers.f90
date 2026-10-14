!> Constructor and parameter tests for generated layer PCM components
module test_model_component_pcm_layers
   use mctc_env, only: wp
   use mctc_io, only: structure_type, new
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, moist_cavity_iswig_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: default_cpcm_radii
   use moist_model_continuum_type, only: model_continuum_type, new_continuum_model
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, new_component_cpcm, &
      & model_continuum_component_cosmo, new_component_cosmo
   use mctc_env_error, only: moist_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check
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
   use moist_model_continuum_component_pcm_solvers, only: pcm_dense_solver
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
         & new_unittest("PySCF-layer-goldens", test_layer_goldens)]
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
      type(moist_context_type), pointer :: absent_ctx => null()
      !> Model and minimal cavity used only for context inheritance
      type(model_continuum_type) :: model
      type(cavity_type_iswig) :: cavity
      !> Empty structure; the missing-context guard must precede geometry access
      type(structure_type) :: mol
      !> Concrete PCM component under test
      class(model_continuum_component_pcm), allocatable :: component
      !> Constructor and context mode
      integer :: kind, mode

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
               call fatal_error(error, "test_optional_context: invalid mode")
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

end module test_model_component_pcm_layers
