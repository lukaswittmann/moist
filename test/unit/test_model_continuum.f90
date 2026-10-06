!> Unit tests for the continuum list-based solvation model
!>
!> Covers the model container itself rather than any single component
!>
!> - a component driven through the model reproduces the procedural result
!> - several components sum, and the lifecycle guards fire
!> - per-component numerics live in the `test_model_component_*` suites
module test_model_continuum
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_c_binding, only: c_null_funptr, c_null_ptr
   use moist_model_continuum_component_pcm_type, only: moist_pcm_parameters_type
   use test_helpers, only: component_view
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_channels_fields, only: field_query_type, field_real
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type, potential_adjoint_response_type
   use moist_model_continuum_component_pcm_type, only: solver_type
   use moist_model_continuum_component_pcm_cpcm, only: model_continuum_component_cpcm, new_component_cpcm
   use moist_model_continuum_component, only: model_continuum_component_pv, new_component_pv, &
      & model_continuum_component_gostshyp, new_component_gostshyp
   use moist_model_continuum_component_type, only: autogpa
   use mctc_io_codata2018, only: Hartree_energy, Bohr_radius
   use moist_model_continuum, only: model_continuum_type, new_continuum_model
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, moist_cavity_iswig_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_isodensity_callback, only: moist_cavity_drop_lsf_isodensity_callback_type
   use moist_cavity_drop_lsf_isodensity_internal, only: moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_type, only: cavity_type
   use moist_radii, only: radius_type_static, new_cosmo_radii
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: build_test_cavity, stage_model_point_charge_energy, &
      & fill_missing_with_zeros, fill_point_charge_field, copy_potential_adjoint, &
      & fill_point_charge_potential, read_printout, printed_entry

   implicit none(type, external)
   private

   public :: collect_model_continuum

   !> Tolerance for values that must agree to roundoff
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   !> Absolute tolerance for equivalent energy
   real(wp), parameter :: thr2 = 1.0e-12_wp

contains

!> Collect the continuum-model test suite
!>
!> @param[out] testsuite Collection of tests
   subroutine collect_model_continuum(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("continuum_model_cpcm", test_continuum_model_smoke), &
         & new_unittest("continuum_model_cpcm_pv", test_continuum_model_pv_smoke), &
         & new_unittest("continuum_model_guards", test_continuum_model_guards), &
         & new_unittest("continuum_model_component_energies", test_continuum_model_component_energies), &
         & new_unittest("continuum_model_parameter_printout", test_continuum_model_parameter_printout) &
         & ]

   end subroutine collect_model_continuum

!> Single-component coverage for the continuum list-based solvation model
   subroutine test_continuum_model_smoke(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(model_continuum_type), target :: model
      type(model_continuum_component_cpcm) :: pcm_component, pcm_reference
      type(cavity_type_iswig) :: cavity
      type(radius_type_static) :: radius_model
      type(coupling_type), pointer :: coupling
      !> Coupling this model never minted, for the pre-update guard
      type(coupling_type), target :: foreign
      type(response_type) :: response
      type(potential_adjoint_response_type), allocatable :: charge
      real(wp) :: energy, reference_energy
      !> Independent forward component gradient and model result
      real(wp), allocatable :: gradient(:, :), reference_gradient(:, :)

      real(wp), parameter :: epsilon = 32.0_wp
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]

      !> Run context owned here and borrowed by the cavity, model and components
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)
      call get_structure(mol, "MB16-43", "01")

      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call new_continuum_model(model, ctx, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Continuum-model construction failed: "//err%message)
         return
      end if

      call check(error, .not. model%is_updated(), more="new model must start invalid")
      if (allocated(error)) return

      ! An accessor must refuse to run before the first update
      ! (no coupling can be built yet, so a foreign one stands in)
      energy = 0.0_wp
      call model%get_energy(foreign, energy, err)
      call check(error, allocated(err), &
         & more="continuum-model energy was available before the first update")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)

      call new_component_cpcm(pcm_component, epsilon=epsilon, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%cholesky), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "Continuum-model CPCM construction failed: "//err%message)
         return
      end if
      call model%add_component(pcm_component, err)
      if (allocated(err)) then
         call test_failed(error, "Adding CPCM component failed: "//err%message)
         return
      end if

      ! Move the solute after building the cavity template to test model update
      mol%xyz(1, :) = mol%xyz(1, :) + 0.25_wp
      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "Continuum-model update failed: "//err%message)
         return
      end if
      call check(error, maxval(abs(model%cavity%sphxyz - mol%xyz)), 0.0_wp, &
         & thr=thr, message="model update did not move the cavity centers")
      if (allocated(error)) return
      call stage_model_point_charge_energy(error, model, qat_vals, mol, coupling)
      if (allocated(error)) return

      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      if (allocated(err)) then
         call test_failed(error, "Continuum-model energy failed: "//err%message)
         return
      end if

      ! Procedural reference on an independently updated cavity
      call new_component_cpcm(pcm_reference, epsilon=epsilon, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%cholesky), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "Reference CPCM construction failed: "//err%message)
         return
      end if
      call cavity%update(mol, err)
      if (.not. allocated(err)) call pcm_reference%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Reference CPCM update failed: "//err%message)
         return
      end if
      reference_energy = 0.0_wp
      call pcm_reference%get_energy(component_view(coupling), cavity, reference_energy, err)
      if (allocated(err)) then
         call test_failed(error, "Reference CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy, reference_energy, thr=thr2, &
         & message="continuum model did not reproduce the procedural CPCM energy")
      if (allocated(error)) return

      energy = 3.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call model%get_energy(coupling, energy, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, energy, 3.0_wp + 2.0_wp*reference_energy, thr=thr2, &
                 message="Repeated energy getters must accumulate")
      if (allocated(error)) return

      ! The response must carry the potential adjoint (the CPCM charges)
      call model%prepare_response(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Continuum-model response staging failed: "//err%message)
         return
      end if
      call model%get_response(coupling, response, err)
      if (allocated(err)) then
         call test_failed(error, "Continuum-model response failed: "//err%message)
         return
      end if
      call copy_potential_adjoint(response, charge)
      call check(error, allocated(charge), &
         & more="continuum model did not expose the CPCM potential adjoint")
      if (allocated(error)) return
      call check(error, maxval(abs(charge%w_phi - pcm_reference%q)), 0.0_wp, &
         & thr=thr2, &
         & message="continuum-model CPCM charges differ from the procedural reference")
      if (allocated(error)) return

      ! Repeated response calls replace the host adjoint rather than doubling it
      call model%get_response(coupling, response, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call copy_potential_adjoint(response, charge)
      call check(error, maxval(abs(charge%w_phi - pcm_reference%q)), 0.0_wp, thr=thr2, &
         & message="Repeated response getter accumulated stale charges")
      if (allocated(error)) return

      call model%prepare_gradient(coupling, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call fill_point_charge_field(model%cavity, coupling, qat_vals, mol)
      allocate (gradient(3, mol%nat), reference_gradient(3, mol%nat), source=0.0_wp)
      call pcm_reference%get_gradient(component_view(coupling), cavity, response, reference_gradient, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call model%get_gradient(coupling, response, gradient, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, maxval(abs(gradient - reference_gradient)), 0.0_wp, thr=thr2, &
         & message="model gradient differs from independent forward component")
      if (allocated(error)) return
      call check(error, .not. allocated(model%cavity%xyz1_rA), &
         & more="reverse model gradient built a surface nuclear Jacobian")
      if (allocated(error)) return
      call copy_potential_adjoint(response, charge)
      call check(error, allocated(charge), more="model gradient omitted host charges")
      if (allocated(error)) return
      call check(error, maxval(abs(charge%w_phi - pcm_reference%q)), 0.0_wp, thr=thr2, &
         & message="model gradient retained stale response charges")
      if (allocated(error)) return

      call model%use_forward_gradient(.true.)
      gradient = 0.0_wp
      call model%get_gradient(coupling, response, gradient, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      ! The reverse path matches to 6e-17 as well; only the Jacobian tells the paths apart
      call check(error, allocated(model%cavity%xyz1_rA), &
         & more="forward selection did not take the forward path")
      if (allocated(error)) return
      call check(error, maxval(abs(gradient - reference_gradient)), 0.0_wp, thr=thr2, &
         & message="forward model gradient differs from independent component")
      if (allocated(error)) return
      call model%use_forward_gradient(.false.)

      ! Components are frozen once the model has been updated
      call model%add_component(pcm_component, err)
      call check(error, allocated(err), &
         & more="adding a component after the first update was not rejected")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)

      ! A failed update must invalidate a model that was previously usable
      deallocate (model%cavity)
      call model%update(mol, err)
      call check(error, allocated(err), &
         & more="continuum-model update without a cavity was not rejected")
      if (allocated(error)) return
      call check(error, .not. model%is_updated(), &
         & more="failed continuum-model update left the model marked usable")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)
      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err), &
         & more="continuum-model energy remained available after a failed update")

   end subroutine test_continuum_model_smoke

!> Smoke test for a two-component (CPCM + PV) continuum solvation model
   subroutine test_continuum_model_pv_smoke(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(model_continuum_type), target :: model_pcm, model_pv, model_zero
      type(model_continuum_component_cpcm) :: pcm_component
      type(cavity_type_iswig) :: cavity
      type(radius_type_static) :: radius_model
      type(coupling_type), pointer :: coupling, coupling_pcm, coupling_zero
      type(response_type) :: response
      type(potential_adjoint_response_type), allocatable :: charge
      !> Borrowed surface views of the three model domains
      class(cavity_type), pointer :: cav_pv, cav_pcm, cav_zero
      real(wp) :: energy_pcm, energy_pv, energy_zero, volume
      real(wp), allocatable :: gradient_pcm(:, :), gradient_pv(:, :)
      real(wp), allocatable :: gradient_zero(:, :), volume_gradient(:, :)

      real(wp), parameter :: epsilon = 32.0_wp
      real(wp), parameter :: pressure = 0.75_wp
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]

      !> Run context owned here and borrowed by the cavity, models and components
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)
      call get_structure(mol, "MB16-43", "01")

      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call new_component_cpcm(pcm_component, epsilon=epsilon, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%cholesky), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "CPCM construction failed: "//err%message)
         return
      end if

      ! Reference model: CPCM only
      call build_pv_model(model_pcm, cavity, ctx, pcm_component, .false., 0.0_wp, mol, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM-only model setup failed: "//err%message)
         return
      end if

      ! Two-component model: CPCM + PV at finite pressure
      call build_pv_model(model_pv, cavity, ctx, pcm_component, .true., pressure, mol, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV model setup failed: "//err%message)
         return
      end if

      ! Two-component model with a vanishing pressure
      call build_pv_model(model_zero, cavity, ctx, pcm_component, .true., 0.0_wp, mol, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV(0) model setup failed: "//err%message)
         return
      end if

      cav_pcm => model_pcm%cavity
      cav_pv => model_pv%cavity
      cav_zero => model_zero%cavity

      ! One coupling for all three models
      call stage_model_point_charge_energy(error, model_pv, qat_vals, mol, coupling)
      if (allocated(error)) return

      volume = cav_pv%total_volume
      call check(error, cav_pcm%ngrid == cav_pv%ngrid .and. cav_zero%ngrid == cav_pv%ngrid, &
         & "the three models were not built on the same cavity surface")
      if (allocated(error)) return

      energy_pcm = 0.0_wp
      call stage_model_point_charge_energy(error, model_pcm, qat_vals, mol, coupling_pcm)
      if (allocated(error)) return
      call model_pcm%get_energy(coupling_pcm, energy_pcm, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM-only energy failed: "//err%message)
         return
      end if
      energy_pv = 0.0_wp
      call model_pv%get_energy(coupling, energy_pv, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV energy failed: "//err%message)
         return
      end if
      energy_zero = 0.0_wp
      call stage_model_point_charge_energy(error, model_zero, qat_vals, mol, coupling_zero)
      if (allocated(error)) return
      call model_zero%get_energy(coupling_zero, energy_zero, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV(0) energy failed: "//err%message)
         return
      end if

      ! The component energies must simply add up
      call check(error, energy_pv, energy_pcm + pressure*volume, thr=thr2, &
         & message="CPCM+PV model energy is not the sum of its components")
      if (allocated(error)) return

      ! A vanishing pressure must leave the CPCM energy untouched
      call check(error, energy_zero, energy_pcm, thr=thr, &
         & message="PV at zero pressure changed the model energy")
      if (allocated(error)) return

      ! The shared surface-adjoint accumulator must survive two components
      call model_pv%prepare_response(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV response staging failed: "//err%message)
         return
      end if
      call model_pv%get_response(coupling, response, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV potential failed: "//err%message)
         return
      end if
      call copy_potential_adjoint(response, charge)
      call check(error, allocated(charge), &
         & more="CPCM+PV model produced no potential adjoint item")
      if (allocated(error)) return

      ! Asking a second time must still carry the CPCM potential adjoint when a
      ! second, non-electrostatic component shares the accumulator
      call model_pv%get_response(coupling, response, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV second response failed: "//err%message)
         return
      end if
      call copy_potential_adjoint(response, charge)
      call check(error, allocated(charge), &
         & more="CPCM+PV model did not expose the CPCM potential adjoint")
      if (allocated(error)) return

      ! Gradient phase
      call model_pv%prepare_gradient(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Gradient staging failed: "//err%message)
         return
      end if
      call fill_missing_with_zeros(cav_pv, coupling)
      call fill_point_charge_field(cav_pv, coupling, qat_vals, mol)

      allocate (gradient_pcm(3, mol%nat), source=0.0_wp)
      allocate (gradient_pv(3, mol%nat), source=0.0_wp)
      allocate (gradient_zero(3, mol%nat), source=0.0_wp)
      allocate (volume_gradient(3, mol%nat), source=0.0_wp)

      call model_pcm%prepare_gradient(coupling_pcm, err)
      call fill_point_charge_field(cav_pcm, coupling_pcm, qat_vals, mol)
      call model_pcm%get_gradient(coupling_pcm, response, gradient_pcm, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM-only gradient failed: "//err%message)
         return
      end if
      call model_pv%get_gradient(coupling, response, gradient_pv, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV gradient failed: "//err%message)
         return
      end if
      call model_zero%prepare_gradient(coupling_zero, err)
      call fill_point_charge_field(cav_zero, coupling_zero, qat_vals, mol)
      call model_zero%get_gradient(coupling_zero, response, gradient_zero, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV(0) gradient failed: "//err%message)
         return
      end if
      ! The reverse-mode gradient contracts the surface adjoints directly
      call cav_pv%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, "Cavity forward gradient failed: "//err%message)
         return
      end if
      if (.not. allocated(cav_pv%v1_rA)) then
         call test_failed(error, "Cavity produced no per-point volume derivatives")
         return
      end if
      volume_gradient = sum(cav_pv%v1_rA, dim=3)

      ! The PV gradient is the pressure-scaled cavity-volume gradient
      call check(error, maxval(abs(gradient_pv - gradient_pcm - pressure*volume_gradient)), &
         & 0.0_wp, thr=thr2, &
         & message="CPCM+PV gradient is not the CPCM gradient plus p*dV/dR")
      if (allocated(error)) return

      ! A vanishing pressure must leave the CPCM gradient untouched
      call check(error, maxval(abs(gradient_zero - gradient_pcm)), 0.0_wp, thr=thr, &
         & message="PV at zero pressure changed the model gradient")

      if (allocated(error)) return
      gradient_zero = 3.0_wp
      call model_pcm%get_gradient(coupling_pcm, response, gradient_zero, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call model_pcm%get_gradient(coupling_pcm, response, gradient_zero, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, maxval(abs(gradient_zero - 3.0_wp - 2.0_wp*gradient_pcm)), &
                 0.0_wp, thr=thr2, message="Repeated gradient getters must accumulate")

   end subroutine test_continuum_model_pv_smoke

!> Foreign-coupling and isodensity-density guard checks
!>
!> A coupling minted by one model is refused by every other model, and a
!> model without an internal isodensity cavity refuses a density
   subroutine test_continuum_model_guards(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err

      !> Molecular structure
      type(structure_type) :: mol
      !> Owning model and a second, independent model
      type(model_continuum_type), target :: model_a, model_b, unconstructed
      !> Only component of both models
      type(model_continuum_component_pv) :: pv_component
      !> Cavity template copied into both models
      type(cavity_type_iswig) :: cavity
      !> Radius model storage
      type(radius_type_static) :: radius_model
      !> Coupling minted by `model_a`
      type(coupling_type), pointer :: coupling
      !> Response list handed to the refused gradient
      type(response_type) :: response
      !> Nuclear-gradient accumulator, must stay untouched
      real(wp), allocatable :: gradient(:, :)
      !> Density matrix offered to a model that cannot take one
      real(wp) :: density(1, 1)
      !> Energy accumulator, must stay untouched
      real(wp) :: energy
      !> Named fields fetched through the public model interface
      type(field_query_type) :: fields

      !> Run context owned here and borrowed by the cavity and models
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)
      call get_structure(mol, "MB16-43", "01")

      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call new_component_pv(pv_component, 0.5_wp)
      call unconstructed%add_component(pv_component, err)
      call expect_error("Construct the model", "unconstructed add")
      if (allocated(error)) return
      call build_model(model_a)
      if (allocated(error)) return
      call build_model(model_b)
      if (allocated(error)) return

      call model_a%new_coupling(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Coupling setup failed: "//err%message)
         return
      end if

      call check(error, model_a%atom_count() == mol%nat, more="wrong model atom count")
      if (allocated(error)) return
      call fields%fetch("sphxyz")
      call model_a%list_fields(fields)
      call check(error, fields%found, more="model did not publish cavity fields")
      if (allocated(error)) return
      call check(error, maxval(abs(fields%rvals - reshape(mol%xyz, [3*mol%nat]))), &
         & 0.0_wp, thr=thr, message="model published wrong sphere centers")
      if (allocated(error)) return
      call model_a%use_forward_gradient(.true.)
      call check(error, model_a%force_forward_gradient)
      if (allocated(error)) return
      call model_a%use_forward_gradient(.false.)
      call check(error, .not. model_a%force_forward_gradient)
      if (allocated(error)) return

      allocate (gradient(3, mol%nat), source=0.0_wp)
      energy = 2.0_wp
      call model_a%get_energy(coupling, energy, err)
      call expect_error("requires a coupling staged", "unstaged energy")
      if (allocated(error)) return
      call model_a%get_response(coupling, response, err)
      call expect_error("requires a coupling staged", "unstaged response")
      if (allocated(error)) return
      call model_a%get_gradient(coupling, response, gradient, err)
      call expect_error("requires a coupling staged", "unstaged gradient")
      if (allocated(error)) return

      call model_b%get_energy(coupling, energy, err)
      call check_foreign("energy")
      if (allocated(error)) return
      call model_b%get_response(coupling, response, err)
      call check_foreign("response")
      if (allocated(error)) return
      call check(error, energy, 2.0_wp, thr=0.0_wp, &
         & more="refused energy wrote into the accumulator")
      if (allocated(error)) return

      ! Both models are updated, so the ownership check is what refuses
      call model_b%prepare_energy(coupling, err)
      call check_foreign("staging")
      if (allocated(error)) return

      call model_b%get_gradient(coupling, response, gradient, err)
      call check_foreign("gradient")
      if (allocated(error)) return
      call check(error, maxval(abs(gradient)), 0.0_wp, thr=0.0_wp, &
         & more="refused gradient wrote into the accumulator")
      if (allocated(error)) return

      ! The owner itself stages the same coupling
      call model_a%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Owner staging failed: "//err%message)
         return
      end if
      call model_a%prepare_gradient(coupling, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call model_a%get_gradient(coupling, response, gradient(:, :mol%nat-1), err)
      call expect_error("shape mismatch", "gradient shape")
      if (allocated(error)) return
      call model_a%invalidate()
      call model_a%get_energy(coupling, energy, err)
      call expect_error("must be updated first", "invalidated energy")
      if (allocated(error)) return
      call model_a%get_response(coupling, response, err)
      call expect_error("must be updated first", "invalidated response")
      if (allocated(error)) return
      call model_a%get_gradient(coupling, response, gradient, err)
      call expect_error("must be updated first", "invalidated gradient")
      if (allocated(error)) return
      call check(error, energy, 2.0_wp, thr=0.0_wp, &
         & more="invalidated energy wrote into the accumulator")
      if (allocated(error)) return
      call model_a%release_coupling(coupling)

      density = 0.0_wp
      call model_b%set_isodensity_density(density, err)
      if (.not. allocated(err)) then
         call test_failed(error, "iSwiG model accepted an isodensity density")
         return
      end if
      call check(error, index(err%message, "requires an internal isodensity cavity") > 0, &
         & more="unexpected error message: "//err%message)
      if (allocated(error)) return
      call check(error, .not. model_b%is_updated(), &
         & more="a refused density left the model marked usable")

   contains

      !> Require a diagnostic and release the library error
      !>
      !> @param[in] text Required diagnostic fragment
      !> @param[in] label Operation name
      subroutine expect_error(text, label)
         !> Required diagnostic fragment
         character(len=*), intent(in) :: text
         !> Operation name
         character(len=*), intent(in) :: label

         if (.not. allocated(err)) then
            call test_failed(error, label//" was accepted")
            return
         end if
         call check(error, index(err%message, text) > 0, more=label//": "//err%message)
         deallocate (err)
      end subroutine expect_error

      !> Assemble and update a PV-only continuum model
      !>
      !> @param[out] model Model to build
      subroutine build_model(model)
         !> Model to build
         type(model_continuum_type), intent(out) :: model

         call new_continuum_model(model, ctx, cavity, err)
         if (.not. allocated(err)) call model%add_component(pv_component, err)
         if (.not. allocated(err)) call model%update(mol, err)
         if (allocated(err)) call test_failed(error, "Model setup failed: "//err%message)
      end subroutine build_model

      !> Require the foreign-coupling refusal for one accessor
      !>
      !> @param[in] label Accessor name used in failure messages
      subroutine check_foreign(label)
         !> Accessor name used in failure messages
         character(len=*), intent(in) :: label

         if (.not. allocated(err)) then
            call test_failed(error, "foreign coupling accepted by "//label)
            return
         end if
         call check(error, index(err%message, "Coupling belongs to a different model") > 0, &
            & more=label//": "//err%message)
         deallocate (err)
      end subroutine check_foreign

   end subroutine test_continuum_model_guards

!> Per-component energies published by the continuum model
!>
!> CPCM(32), PV(p), CPCM(4), PV(0): repeated names told apart by index, a
!> zero contribution, and a second component that fails once the cavity
!> volume is taken away
   subroutine test_continuum_model_component_energies(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err

      !> Molecular structure
      type(structure_type) :: mol
      !> Model under test
      type(model_continuum_type), target :: model
      !> CPCM templates at two dielectric constants
      type(model_continuum_component_cpcm) :: pcm_high, pcm_low
      !> PV templates at a finite and a vanishing pressure
      type(model_continuum_component_pv) :: pv_finite, pv_zero
      !> Cavity template copied into the model
      type(cavity_type_iswig) :: cavity
      !> Radius model storage
      type(radius_type_static) :: radius_model
      !> Coupling minted by the model
      type(coupling_type), pointer :: coupling
      !> Component name and description
      character(len=:), allocatable :: name, description
      !> Energy accumulator
      real(wp) :: energy
      !> Component energies of the first and the doubled-charge evaluation
      real(wp) :: first(4), second(4)
      !> Energy of one component
      real(wp) :: value
      !> Whether a component lists its energy
      logical :: found
      !> Component index
      integer :: i

      !> Nonzero accumulator seed; a failed evaluation must leave it bit for bit
      real(wp), parameter :: seed = 1.25_wp
      real(wp), parameter :: pressure = 0.75_wp
      character(len=4), parameter :: names(4) = [character(len=4) :: "CPCM", "PV", "CPCM", "PV"]
      character(len=*), parameter :: cpcm_about = "Conductor-like polarizable continuum, f(eps) = (eps - 1)/eps"
      character(len=*), parameter :: pv_about = "Pressure-volume work, pressure times cavity volume"
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]

      !> Run context owned here and borrowed by the cavity, model and components
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)
      call get_structure(mol, "MB16-43", "01")

      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call new_component_cpcm(pcm_high, epsilon=32.0_wp, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%cholesky), ctx=ctx)
      if (.not. allocated(err)) call new_component_cpcm(pcm_low, epsilon=4.0_wp, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%cholesky), ctx=ctx)
      call new_component_pv(pv_finite, pressure)
      call new_component_pv(pv_zero, 0.0_wp)
      if (.not. allocated(err)) call new_continuum_model(model, ctx, cavity, err)
      if (.not. allocated(err)) call model%add_component(pcm_high, err)
      if (.not. allocated(err)) call model%add_component(pv_finite, err)
      if (.not. allocated(err)) call model%add_component(pcm_low, err)
      if (.not. allocated(err)) call model%add_component(pv_zero, err)
      if (allocated(err)) then
         call test_failed(error, "Model setup failed: "//err%message)
         return
      end if

      ! Count and names before the first update; nothing is evaluated yet
      call check(error, model%component_count() == 4, more="wrong component count")
      if (allocated(error)) return
      do i = 1, 4
         call model%component_name(i, name)
         call check(error, allocated(name), more="component name missing before update")
         if (allocated(error)) return
         call check(error, name == trim(names(i)) .and. len(name) == len_trim(names(i)), &
            & more="wrong component name: "//name)
         if (allocated(error)) return
         call model%component_description(i, description)
         call check(error, allocated(description), more="component description missing before update")
         if (allocated(error)) return
         if (mod(i, 2) == 1) then
            call check(error, description == cpcm_about, more="wrong CPCM description: "//description)
         else
            call check(error, description == pv_about, more="wrong PV description: "//description)
         end if
         if (allocated(error)) return
      end do
      call model%component_name(0, name)
      call check(error, .not. allocated(name), more="named a component below the list")
      if (allocated(error)) return
      call model%component_name(5, name)
      call check(error, .not. allocated(name), more="named a component past the list")
      if (allocated(error)) return
      call model%component_description(5, description)
      call check(error, .not. allocated(description), more="described a component past the list")
      if (allocated(error)) return
      call check_unavailable("before update")
      if (allocated(error)) return

      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "Model update failed: "//err%message)
         return
      end if
      call check_unavailable("before evaluation")
      if (allocated(error)) return

      ! Every component publishes its own energy; their ordered sum is the total
      call stage_model_point_charge_energy(error, model, qat_vals, mol, coupling)
      if (allocated(error)) return
      energy = seed
      call model%get_energy(coupling, energy, err)
      if (allocated(err)) then
         call test_failed(error, "Model energy failed: "//err%message)
         return
      end if
      call read_all("first evaluation", first)
      if (allocated(error)) return
      call check(error, ieee_is_finite(energy), more="model energy is not finite")
      if (allocated(error)) return
      call check(error, energy - seed, ((first(1) + first(2)) + first(3)) + first(4), thr=thr2, &
         & message="component energies do not add up to the model energy")
      if (allocated(error)) return
      call check(error, first(2), pressure*model%cavity%total_volume, thr=thr, &
         & message="PV energy is not pressure times volume")
      if (allocated(error)) return
      call check(error, first(4) == 0.0_wp, more="PV at zero pressure published a nonzero energy")
      if (allocated(error)) return
      ! The two CPCM components differ only in f(epsilon) = (epsilon - 1)/epsilon
      call check(error, first(1) /= 0.0_wp .and. first(3) /= 0.0_wp, more="CPCM energies vanish")
      if (allocated(error)) return
      call check(error, first(3)/first(1), (3.0_wp/4.0_wp)/(31.0_wp/32.0_wp), thr=thr2, &
         & message="repeated CPCM names were not told apart by index")
      if (allocated(error)) return

      ! An evaluation refused at the phase check keeps the stored energies
      call model%prepare_response(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Response staging failed: "//err%message)
         return
      end if
      energy = seed
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err), more="energy of a response-staged coupling was accepted")
      if (allocated(error)) return
      deallocate (err)
      call check(error, energy == seed, more="refused evaluation wrote into the accumulator")
      if (allocated(error)) return
      call read_all("refused evaluation", second)
      if (allocated(error)) return
      call check(error, all(second == first), more="refused evaluation changed the stored energies")
      if (allocated(error)) return

      ! A second evaluation replaces the first: doubled charges, four times the CPCM energy
      call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Energy restaging failed: "//err%message)
         return
      end if
      call fill_point_charge_potential(model%cavity, coupling, 2.0_wp*qat_vals, mol)
      energy = seed
      call model%get_energy(coupling, energy, err)
      if (allocated(err)) then
         call test_failed(error, "Second model energy failed: "//err%message)
         return
      end if
      call read_all("second evaluation", second)
      if (allocated(error)) return
      call check(error, second(1), 4.0_wp*first(1), thr=thr2, &
         & message="second evaluation did not replace the first CPCM energy")
      if (allocated(error)) return
      call check(error, second(3), 4.0_wp*first(3), thr=thr2, &
         & message="second evaluation did not replace the second CPCM energy")
      if (allocated(error)) return
      call check(error, second(2) == first(2) .and. second(4) == 0.0_wp, &
         & more="PV energies changed with the host potential")
      if (allocated(error)) return

      ! Fail fast: without a cavity volume the first PV component fails
      deallocate (model%cavity%total_volume)
      energy = seed
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err), more="energy without a cavity volume was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Cavity volume is unavailable") > 0, &
         & more="unexpected failure: "//err%message)
      if (allocated(error)) return
      deallocate (err)
      call check(error, energy == seed, more="failed evaluation wrote into the accumulator")
      if (allocated(error)) return
      call read_energy(1, found, value)
      call check(error, found, more="component before the failure lost its energy")
      if (allocated(error)) return
      call check(error, value, second(1), thr=thr2, &
         & message="component before the failure has a wrong energy")
      if (allocated(error)) return
      do i = 2, 4
         call read_energy(i, found, value)
         call check(error, .not. found, more="failing or later component still lists an energy")
         if (allocated(error)) return
      end do

      ! An update clears every energy and restores the volume
      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "Model re-update failed: "//err%message)
         return
      end if
      call check_unavailable("after update")
      if (allocated(error)) return

      ! So does an update that fails before reaching the components
      call evaluate("evaluation before a failed update")
      if (allocated(error)) return
      deallocate (model%cavity)
      call model%update(mol, err)
      call check(error, allocated(err), more="update without a cavity was accepted")
      if (allocated(error)) return
      deallocate (err)
      call check_unavailable("after a failed update")
      if (allocated(error)) return
      call model%release_coupling(coupling)

   contains

      !> Stage the point-charge energy phase and evaluate it
      !>
      !> @param[in] label Situation named in failure messages
      subroutine evaluate(label)
         !> Situation named in failure messages
         character(len=*), intent(in) :: label
         !> Component energies, read only to require them listed
         real(wp) :: values(4)

         call model%prepare_energy(coupling, err)
         if (.not. allocated(err)) then
            call fill_point_charge_potential(model%cavity, coupling, qat_vals, mol)
            energy = seed
            call model%get_energy(coupling, energy, err)
         end if
         if (allocated(err)) then
            call test_failed(error, label//" failed: "//err%message)
            return
         end if
         call read_all(label, values)
      end subroutine evaluate

      !> Fetch the stored energy of one component by name
      !>
      !> @param[in] index 1-based component index
      !> @param[out] listed Whether the component lists its energy
      !> @param[out] stored Energy, zero when not listed
      subroutine read_energy(index, listed, stored)
         !> 1-based component index
         integer, intent(in) :: index
         !> Whether the component lists its energy
         logical, intent(out) :: listed
         !> Energy, zero when not listed
         real(wp), intent(out) :: stored
         !> Field walker
         type(field_query_type) :: query

         stored = 0.0_wp
         call query%fetch("energy")
         call model%list_component_fields(index, query)
         listed = query%found
         if (listed) stored = query%rvals(1)
      end subroutine read_energy

      !> Require every component to list nothing, by enumeration and by name
      !>
      !> @param[in] label Situation named in the failure message
      subroutine check_unavailable(label)
         !> Situation named in the failure message
         character(len=*), intent(in) :: label
         !> Field walker
         type(field_query_type) :: query
         !> Component index
         integer :: j
         !> Whether the energy is listed
         logical :: listed
         !> Fetched energy
         real(wp) :: stored

         do j = 1, model%component_count()
            call query%enumerate()
            call model%list_component_fields(j, query)
            call read_energy(j, listed, stored)
            call check(error, query%nfield == 0 .and. .not. listed, &
               & more=label//": a component still lists an energy")
            if (allocated(error)) return
         end do
      end subroutine check_unavailable

      !> Read every component energy, each a listed, finite real scalar
      !>
      !> @param[in] label Situation named in the failure message
      !> @param[out] values Component energies in list order
      subroutine read_all(label, values)
         !> Situation named in the failure message
         character(len=*), intent(in) :: label
         !> Component energies in list order
         real(wp), intent(out) :: values(:)
         !> Field walker
         type(field_query_type) :: query
         !> Component index
         integer :: j
         !> Whether the energy is listed
         logical :: listed

         values = 0.0_wp
         do j = 1, size(values)
            call query%enumerate()
            call model%list_component_fields(j, query)
            call check(error, query%nfield == 1, more=label//": a component does not list one result")
            if (allocated(error)) return
            call check(error, query%info(1)%name == "energy" .and. query%info(1)%rank == 0 &
               & .and. query%info(1)%dtype == field_real, more=label//": wrong energy descriptor")
            if (allocated(error)) return
            call read_energy(j, listed, values(j))
            call check(error, listed .and. ieee_is_finite(values(j)), &
               & more=label//": a component energy is missing or not finite")
            if (allocated(error)) return
         end do
      end subroutine read_all

   end subroutine test_continuum_model_component_energies

!> Settings printout of a model, before any update
!>
!> iSwiG + CPCM + PV + GOSTSHYP: the cavity section, then one section per
!> component headed by its 1-based index, with its description, `scale`
!> only when it is not 1, the solvent inputs and the registered settings of
!> the copy the model holds; a DROP cavity is labelled by its level set and
!> prints its grouped settings
   subroutine test_continuum_model_parameter_printout(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err

      !> Model under test
      type(model_continuum_type), target :: model
      !> Component templates
      type(model_continuum_component_cpcm) :: pcm
      type(model_continuum_component_pv) :: pv
      type(model_continuum_component_gostshyp) :: gostshyp
      !> Cavity templates
      type(cavity_type_iswig) :: cavity
      type(cavity_type_drop) :: drop
      !> Level sets of the DROP cavities
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_cfc_type) :: cfc
      type(moist_cavity_drop_lsf_isodensity_callback_type) :: rho_callback
      type(moist_cavity_drop_lsf_isodensity_internal_type) :: rho_internal
      !> Radius model storage
      type(radius_type_static) :: radius_model
      !> Lines read back from the scratch unit
      character(len=160) :: lines(120)
      !> Line count and scratch unit
      integer :: nline, unit
      !> Expected GPa values, in the printer's fixed format
      character(len=16) :: pv_gpa, gostshyp_gpa
      !> Run context owned here and borrowed by the cavities, model and components
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cosmo_radii(radius_model)
      call new_cavity_iswig(cavity, radius_model, err, moist_cavity_iswig_parameters_type(num_leb=50), ctx=ctx)
      if (.not. allocated(err)) call new_component_cpcm(pcm, epsilon=4.0_wp, error=err, &
         & param=moist_pcm_parameters_type(solver=solver_type%iterative, solver_tol=1.0e-9_wp, &
         & solver_maxiter=77), ctx=ctx)
      call new_component_pv(pv, 2.5e-4_wp)
      call new_component_gostshyp(gostshyp, 3.0e-5_wp)
      gostshyp%scale = 0.5_wp
      if (.not. allocated(err)) call new_continuum_model(model, ctx, cavity, err)
      if (.not. allocated(err)) call model%add_component(pcm, err)
      if (.not. allocated(err)) call model%add_component(pv, err)
      if (.not. allocated(err)) call model%add_component(gostshyp, err)
      if (allocated(err)) then
         call test_failed(error, "Model setup failed: "//err%message)
         return
      end if

      ! The model's copy keeps the solver settings in its parameter object
      select type (item => model%components(1)%item)
      type is (model_continuum_component_cpcm)
         call check(error, item%param%solver == solver_type%iterative .and. item%param%solver_maxiter == 77, &
            & more="model copy lost the PCM solver settings")
         if (allocated(error)) return
         call check(error, item%param%solver_tol, 1.0e-9_wp, thr=1.0e-15_wp, &
            & message="model copy lost the PCM solver tolerance")
      class default
         call test_failed(error, "first component is not CPCM")
      end select
      if (allocated(error)) return

      open (newunit=unit, status="scratch", action="readwrite")
      call model%print_parameters(unit)
      call read_printout(unit, lines, nline)
      close (unit)

      call require("Cavity (vdW iSwiG):")
      call require("Component 1 (CPCM):")
      call require("Conductor-like polarizable continuum, f(eps) = (eps - 1)/eps")
      call require("Component 2 (PV):")
      call require("Pressure-volume work, pressure times cavity volume")
      call require("Component 3 (GOSTSHYP):")
      call require("Gaussians on surface tesserae to simulate hydrostatic pressure")
      if (allocated(error)) return
      call check(error, printed_entry(lines(:nline), "num_leb", "50") .and. &
         & printed_entry(lines(:nline), "cut_f", "1.00E-10"), more="iSwiG settings missing")
      if (allocated(error)) return
      call check(error, printed_entry(lines(:nline), "Epsilon", "4.000000") .and. &
         & printed_entry(lines(:nline), "f(eps)", "0.750000"), more="PCM inputs missing")
      if (allocated(error)) return
      call check(error, printed_entry(lines(:nline), "solver", "4") .and. &
         & printed_entry(lines(:nline), "solver_tol", "1.00E-09") .and. &
         & printed_entry(lines(:nline), "solver_maxiter", "77"), more="PCM settings missing")
      if (allocated(error)) return
      ! Pressures in atomic units, then in GPa through the one conversion factor;
      ! it agrees with Eh/a_0^3 to the rounding of the CODATA values (4e-12)
      call check(error, autogpa, Hartree_energy/Bohr_radius**3*1.0e-9_wp, rel=.true., thr=1.0e-11_wp, &
         & message="Eh/bohr^3 -> GPa factor")
      if (allocated(error)) return
      write (pv_gpa, "(f16.6)") 2.5e-4_wp*autogpa
      write (gostshyp_gpa, "(f16.6)") 3.0e-5_wp*autogpa
      call check(error, printed_entry(lines(:nline), "Pressure", "2.50E-04 Eh/bohr^3") .and. &
         & printed_entry(lines(:nline), "Pressure", "3.00E-05 Eh/bohr^3"), more="pressures missing")
      if (allocated(error)) return
      call check(error, printed_entry(lines(:nline), "Pressure", trim(adjustl(pv_gpa))//" GPa") .and. &
         & printed_entry(lines(:nline), "Pressure", trim(adjustl(gostshyp_gpa))//" GPa"), &
         & more="pressures in GPa missing: "//trim(adjustl(pv_gpa))//", "//trim(adjustl(gostshyp_gpa)))
      if (allocated(error)) return
      ! Only GOSTSHYP was rescaled
      call check(error, printed_entry(lines(:nline), "scale", "0.500000") .and. &
         & count(index(adjustl(lines(:nline)), "scale ") == 1) == 1, more="scale printed for an unscaled component")
      if (allocated(error)) return
      ! Sections follow the list order
      call check(error, first_line("Cavity (") < first_line("Component 1 ") .and. &
         & first_line("Component 1 ") < first_line("Component 2 ") .and. &
         & first_line("Component 2 ") < first_line("Component 3 "), more="sections out of order")
      if (allocated(error)) return

      ! Every level-set model is named; the two isodensity models share a name
      call lsf%new()
      call cfc%new()
      call rho_callback%new(c_null_funptr, c_null_ptr)
      call rho_internal%new([1], [0], [1], [1.0_wp], [1.0_wp], err)
      if (allocated(err)) then
         call test_failed(error, "Isodensity setup failed: "//err%message)
         return
      end if
      call check(error, lsf%name == "SvdW" .and. cfc%name == "CFC" .and. rho_callback%name == "Isodensity" &
         & .and. rho_internal%name == "Isodensity", more="unnamed or misnamed level-set model")
      if (allocated(error)) return

      ! DROP is labelled by its level set and prints its grouped settings, without an update
      call new_cavity_drop(drop, radius_model, cfc, err, ctx=ctx)
      if (.not. allocated(err)) then
         call check(error, drop%label == "CFC-DROP", more="CFC DROP label: "//drop%label)
         if (allocated(error)) return
         call new_cavity_drop(drop, radius_model, lsf, err, ctx=ctx)
      end if
      if (allocated(err)) then
         call test_failed(error, "DROP setup failed: "//err%message)
         return
      end if
      call check(error, drop%label == "SvdW-DROP", more="SvdW DROP label: "//drop%label)
      if (allocated(error)) return
      open (newunit=unit, status="scratch", action="readwrite")
      call drop%print_parameters(unit)
      call read_printout(unit, lines, nline)
      close (unit)
      call require("Cavity (SvdW-DROP):")
      call require("DROP Parameters:")
      call require("Implicit surface (SvdW):")
      if (allocated(error)) return
      call check(error, printed_entry(lines(:nline), "Level set", "SvdW") .and. &
         & printed_entry(lines(:nline), "Number of Leb. points", "194"), more="DROP settings missing")
      if (allocated(error)) return

      ! The density source tells the two isodensity models apart
      call new_cavity_drop(drop, radius_model, rho_callback, err, ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "Isodensity DROP setup failed: "//err%message)
         return
      end if
      open (newunit=unit, status="scratch", action="readwrite")
      call drop%print_parameters(unit)
      call read_printout(unit, lines, nline)
      close (unit)
      call require("Cavity (Isodensity-DROP):")
      call require("Implicit surface (isodensity):")
      if (allocated(error)) return
      call check(error, printed_entry(lines(:nline), "Level set", "Isodensity (callback)"), &
         & more="isodensity density source missing")

   contains

      !> Require one printed line to hold `text`
      !>
      !> @param[in] text Expected text
      subroutine require(text)
         !> Expected text
         character(len=*), intent(in) :: text

         if (allocated(error)) return
         call check(error, first_line(text) > 0, more="printout lacks '"//text//"'")
      end subroutine require

      !> Index of the first printed line holding `text`, 0 when none does
      !>
      !> @param[in] text Text to find
      integer function first_line(text) result(iline)
         !> Text to find
         character(len=*), intent(in) :: text
         !> Line index
         integer :: i

         iline = 0
         do i = 1, nline
            if (index(lines(i), text) > 0) then
               iline = i
               return
            end if
         end do
      end function first_line

   end subroutine test_continuum_model_parameter_printout

!> Assemble an updated continuum model from CPCM and an optional PV component
!>
!> PV is appended at the requested pressure when `with_pv` is set
!>
!> @param[out] model Model to build
!> @param[in] cavity Cavity template copied into the model
!> @param[in] ctx Run context owned by the caller
!> @param[in] pcm_component CPCM component template
!> @param[in] with_pv Whether to append a PV component
!> @param[in] pressure Pressure of the PV component
!> @param[in] mol Molecular structure
!> @param[out] error Error handling
   subroutine build_pv_model(model, cavity, ctx, pcm_component, with_pv, pressure, mol, error)

      !> Model to build
      type(model_continuum_type), intent(out) :: model

      !> Cavity template copied into the model
      type(cavity_type_iswig), intent(in) :: cavity

      !> Run context owned by the caller
      type(moist_context_type), intent(in), target :: ctx

      !> CPCM component template
      type(model_continuum_component_cpcm), intent(in) :: pcm_component

      !> Whether to append a PV component
      logical, intent(in) :: with_pv

      !> Pressure of the PV component
      real(wp), intent(in) :: pressure

      !> Molecular structure
      type(structure_type), intent(in) :: mol

      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      type(model_continuum_component_pv) :: pv_component

      call new_continuum_model(model, ctx, cavity, error)
      if (allocated(error)) return
      call model%add_component(pcm_component, error)
      if (allocated(error)) return
      if (with_pv) then
         call new_component_pv(pv_component, pressure)
         call model%add_component(pv_component, error)
         if (allocated(error)) return
      end if
      call model%update(mol, error)

   end subroutine build_pv_model

end module test_model_continuum
