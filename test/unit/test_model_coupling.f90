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
   use moist_model_component_pcm_type, only: moist_pcm_parameters_type
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_channels_coupling, only: coupling_type, coupling_request_type, &
      gaussian_potential_request_type, gaussian_moment_request_type, request_name_len
   use moist_channels_response, only: response_type, potential_adjoint_response_type, &
      & response_accumulate
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_type, only: cavity_type
   use moist_model_component_pcm_type, only: solver_type
   use moist_model_component_pcm_cpcm, only: solvation_model_component_cpcm, new_component_cpcm
   use moist_model_components, only: solvation_model_component_pv, new_component_pv, &
      & solvation_model_component_gostshyp, new_component_gostshyp
   use moist_model_general, only: solvation_model_general, new_model_general
   use moist_cavity_iswig, only: cavity_type_iswig
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: radius_type_static, default_cpcm_radii
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: build_test_cavity, fill_point_charge_potential, fill_point_charge_field, &
      & copy_potential_adjoint
   use test_model_component_helper, only: sw_ngrid => fixture_ngrid_param, &
      & sw_areas => fixture_areas_param, sw_xis => fixture_xis_param, &
      & sw_fs => fixture_fs_param, sw_xyz => fixture_xyz_param, &
      & sw_normals => fixture_normals_param

   implicit none
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

   !> DROP (SvdW) parameters of the O-C-H fixture, as in the golden suite
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
         new_unittest("model_update_ends_walk", test_model_update_walk), &
         new_unittest("model_accessors_name_the_problem", test_model_messages)]
   end subroutine collect_model_coupling

   subroutine test_invalidation(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(solvation_model_general), target :: model
      type(coupling_type), pointer :: first, second
      real(wp) :: energy
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(model, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      call model%new_coupling(first, err)
      call model%new_coupling(second, err)
      call model%prepare_energy(first, err)
      call model%prepare_energy(second, err)
      call fill_point_charge_potential(model%cavity, first, qat_vals, mol)
      call fill_point_charge_potential(model%cavity, second, qat_vals, mol)
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
      type(solvation_model_general), target :: original
      type(solvation_model_general), allocatable, target :: copied
      type(coupling_type), pointer :: coupling
      real(wp) :: energy
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(original, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      call original%new_coupling(coupling, err)
      call original%prepare_energy(coupling, err)
      call fill_point_charge_potential(original%cavity, coupling, qat_vals, mol)
      allocate(copied, source=original)
      energy = 0.0_wp
      call copied%get_energy(coupling, energy, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      deallocate(copied)
      call original%get_energy(coupling, energy, err)
      if (failed(error, err, "original survives destroying its copy")) return
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
      type(solvation_model_general), target :: model
      type(solvation_model_component_cpcm) :: pcm1, pcm2
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      real(wp) :: energies(4)
      real(wp), allocatable :: gradients(:, :, :)
      integer :: k
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call new_component_cpcm(pcm1, ctx, epsilon=2.0_wp, error=err)
      call new_component_cpcm(pcm2, ctx, epsilon=30.0_wp, error=err)
      allocate(gradients(3, mol%nat, 4), source=0.0_wp)
      energies = 0.0_wp
      do k = 1, 4
         call new_model_general(model, cavity, ctx, err)
         if (k == 1 .or. k == 3) call model%add_component(pcm1, err)
         if (k /= 3) call model%add_component(pcm2, err)
         if (k == 2) call model%add_component(pcm1, err)
         call model%update(mol, err)
         if (failed(error, err, "two PCM update")) return
         call model%new_coupling(coupling, err)
         call model%prepare_energy(coupling, err)
         call check(error, count_visits(coupling), 1, more="two PCMs must share one potential request")
         if (allocated(error)) return
         call fill_point_charge_potential(model%cavity, coupling, qat_vals, mol)
         call model%get_energy(coupling, energies(k), err)
         if (failed(error, err, "two PCM energy")) return
         call model%prepare_gradient(coupling, err)
         call fill_point_charge_field(model%cavity, coupling, qat_vals, mol)
         call model%get_gradient(coupling, response, gradients(:,:,k), err)
         if (failed(error, err, "two PCM gradient")) return
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
      type(solvation_model_general), target :: model
      type(solvation_model_component_gostshyp) :: component
      type(coupling_type), pointer :: coupling
      real(wp) :: energy
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call new_model_general(model, cavity, ctx, err)
      call new_component_gostshyp(component, test_pressure)
      call model%add_component(component, err)
      call model%update(mol, err)
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
         call check(error, size(item%width), model%cavity%ngrid)
         if (allocated(error)) return
         call check(error, all(item%width >= 0.0_wp))
      class default
         call test_failed(error, "unexpected request "//item%name())
      end select
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="GOSTSHYP declares one request")
      if (allocated(error)) return
      call new_model_general(model, cavity, ctx, err)
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
   end subroutine test_moments

   subroutine test_pv(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(radius_type_static) :: radii
      type(cavity_type_iswig) :: cavity
      type(solvation_model_general), target :: model
      type(solvation_model_component_pv) :: component
      type(coupling_type), pointer :: coupling
      real(wp) :: energy
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call new_model_general(model, cavity, ctx, err)
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

   !> MB16-43 "01" on the iSwiG fixture cavity of `test_model_general`
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

      call new_context(ctx)
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
      call new_context(ctx, verbosity=0)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=drop_blend_k, &
            blend_3b=drop_blend_3b))
         call new_cavity_drop(cavity, ctx, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
            error=err, param=moist_cavity_drop_parameters_type(num_leb=nleb_drop, tolerance=drop_proj_tol, &
            proj_maxiter=drop_proj_maxiter, proj_level=drop_proj_level, wleb_prune_level=drop_prune_level))
      end block
      if (failed(error, err, "DROP setup")) return
      call cavity%properties(do_fine=.true.)

   end subroutine drop_fixture

   !> Updated general model with a CPCM component and optionally a PV component
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
      type(solvation_model_general), intent(out) :: model
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
      type(solvation_model_component_cpcm) :: cpcm
      type(solvation_model_component_pv) :: pv

      call new_model_general(model, cavity, ctx, err)
      if (failed(error, err, "model setup")) return
      if (with_pv) then
         call new_component_pv(pv, test_pressure)
         call model%add_component(pv, err)
         if (failed(error, err, "adding PV")) return
      end if
      call new_component_cpcm(cpcm, ctx, epsilon=eps_model, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%cholesky))
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
   !> @param[out] error Test failure
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
      type(solvation_model_general), target :: model, other
      type(coupling_type), pointer :: coupling, foreign
      !> Previous response and sentinel charge
      type(response_type) :: response
      type(potential_adjoint_response_type) :: charge
      type(potential_adjoint_response_type), allocatable :: saved
      !> Gradient accumulator, unchanged on failure
      real(wp), allocatable :: gradient(:, :)
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
      charge%w_phi = [7.0_wp]
      call response_accumulate(response, charge, err)
      allocate (gradient(3, mol%nat), source=9.0_wp)
      do scenario = 1, 4
         select case (scenario)
         case (1)
            call model%get_response(foreign, response, err)
         case (2)
            call model%get_gradient(coupling, response, gradient, err)
         case (3)
            call model%prepare_response(coupling, err)
            call model%get_response(coupling, response, err)
         case (4)
            call model%update(mol, err)
            call model%get_gradient(coupling, response, gradient, err)
         end select
         call check(error, allocated(err))
         if (allocated(error)) return
         call copy_potential_adjoint(response, saved)
         call check(error, allocated(saved))
         if (allocated(error)) return
         call check(error, saved%w_phi(1), 7.0_wp, thr=0.0_wp)
         if (allocated(error)) return
         call check(error, maxval(abs(gradient - 9.0_wp)), 0.0_wp, thr=0.0_wp)
         if (allocated(error)) return
      end do
      call model%release_coupling(coupling)
      call other%release_coupling(foreign)
   end subroutine test_preserve_response

   !> A bare DROP cavity reports its missing level set without dereferencing it
   !>
   !> @param[out] error Test failure
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

   !> Scientific name of the current request, "no_current_request" outside a `next()` window
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
   !>
   !> @param[out] error Test failure
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
      type(solvation_model_general), target :: model
      type(solvation_model_component_gostshyp) :: gostshyp
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

   !> A model update ends every walk; answering is refused until the next staging
   !>
   !> @param[out] error Test failure
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
      type(solvation_model_general), target :: model
      type(coupling_type), pointer :: coupling
      !> Potential answer on the current grid
      real(wp), allocatable :: phi(:)
      real(wp) :: energy
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(model, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      call model%new_coupling(coupling, err)
      call model%prepare_energy(coupling, err)
      if (failed(error, err, "energy staging")) return
      call check(error, coupling%next())
      if (allocated(error)) return
      call model%update(mol, err)
      if (failed(error, err, "model update")) return
      call check(error, current_name(coupling), "no_current_request", more="the update kept a current request")
      if (allocated(error)) return
      allocate (phi(model%cavity%ngrid), source=0.0_wp)
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
   !>
   !> @param[out] error Test failure
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
      type(solvation_model_general), target :: model
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      real(wp), allocatable :: gradient(:, :)
      real(wp) :: energy
      call iswig_fixture(ctx, mol, radii, cavity, error)
      if (allocated(error)) return
      call cpcm_model(model, cavity, ctx, mol, .false., error)
      if (allocated(error)) return
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
      call fill_point_charge_potential(model%cavity, coupling, qat_vals, mol)
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

end module test_model_coupling
