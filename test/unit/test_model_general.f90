!> Unit tests for the general list-based solvation model
!>
!> Covers the model container itself rather than any single component: that a
!> component driven through the model reproduces the procedural result, that
!> several components sum, and that the lifecycle guards fire. Per-component
!> numerics live in the `test_model_component_*` suites
module test_model_general
   use moist_model_component_pcm_type, only: moist_pcm_parameters_type
   use test_helpers, only: component_view
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type, potential_adjoint_response_type
   use moist_model_component_pcm_type, only: solver_type
   use moist_model_component_pcm_cpcm, only: solvation_model_component_cpcm, new_component_cpcm
   use moist_model_components, only: solvation_model_component_pv, new_component_pv
   use moist_model_general, only: solvation_model_general, new_model_general
   use moist_cavity_iswig, only: cavity_type_iswig
   use moist_radii, only: radius_type_static
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: build_test_cavity, stage_model_point_charge_energy, &
      & fill_missing_with_zeros, fill_point_charge_field, copy_potential_adjoint

   implicit none(type, external)
   private

   public :: collect_model_general

   !> Tolerance for values that must agree to roundoff
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   !> Loose tolerance for values that pass through a linear solve
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

!> Collect the general-model test suite
   subroutine collect_model_general(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("general_model_cpcm", test_general_model_smoke), &
         & new_unittest("general_model_cpcm_pv", test_general_model_pv_smoke), &
         & new_unittest("general_model_guards", test_general_model_guards) &
         & ]

   end subroutine collect_model_general

!> Single-component coverage for the general list-based solvation model
   subroutine test_general_model_smoke(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(solvation_model_general), target :: model
      type(solvation_model_component_cpcm) :: pcm_component, pcm_reference
      type(cavity_type_iswig) :: cavity
      type(radius_type_static) :: radius_model
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      type(potential_adjoint_response_type), allocatable :: charge
      real(wp) :: energy, reference_energy

      real(wp), parameter :: epsilon = 32.0_wp
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]

      !> Run context owned here and borrowed by the cavity, model and components
      type(moist_context_type), target :: ctx

      call new_context(ctx)
      call get_structure(mol, "MB16-43", "01")

      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call new_model_general(model, cavity, ctx, err)
      if (allocated(err)) then
         call test_failed(error, "General-model construction failed: "//err%message)
         return
      end if

      ! An accessor must refuse to run before the first update
      ! (and before a coupling can be built)
      nullify (coupling)
      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err), &
         & more="general-model energy was available before the first update")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)

      call new_component_cpcm(pcm_component, ctx, epsilon=epsilon, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%cholesky))
      if (allocated(err)) then
         call test_failed(error, "General-model CPCM construction failed: "//err%message)
         return
      end if
      call model%add_component(pcm_component, err)
      if (allocated(err)) then
         call test_failed(error, "Adding CPCM component failed: "//err%message)
         return
      end if

      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "General-model update failed: "//err%message)
         return
      end if
      call stage_model_point_charge_energy(error, model, qat_vals, mol, coupling)
      if (allocated(error)) return

      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      if (allocated(err)) then
         call test_failed(error, "General-model energy failed: "//err%message)
         return
      end if

      ! Procedural reference on an independently updated cavity
      call new_component_cpcm(pcm_reference, ctx, epsilon=epsilon, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%cholesky))
      if (allocated(err)) then
         call test_failed(error, "Reference CPCM construction failed: "//err%message)
         return
      end if
      call pcm_reference%update(mol, cavity, err)
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
         & message="general model did not reproduce the procedural CPCM energy")
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
         call test_failed(error, "General-model response staging failed: "//err%message)
         return
      end if
      call model%get_response(coupling, response, err)
      if (allocated(err)) then
         call test_failed(error, "General-model response failed: "//err%message)
         return
      end if
      call copy_potential_adjoint(response, charge)
      call check(error, allocated(charge), &
         & more="general model did not expose the CPCM potential adjoint")
      if (allocated(error)) return
      call check(error, maxval(abs(charge%w_phi - pcm_reference%q)), 0.0_wp, &
         & thr=thr2, &
         & message="general-model CPCM charges differ from the procedural reference")
      if (allocated(error)) return

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
         & more="general-model update without a cavity was not rejected")
      if (allocated(error)) return
      call check(error, .not. model%updated, &
         & more="failed general-model update left the model marked usable")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)
      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err), &
         & more="general-model energy remained available after a failed update")

   end subroutine test_general_model_smoke

!> Smoke test for a two-component (CPCM + PV) general solvation model
   subroutine test_general_model_pv_smoke(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(solvation_model_general), target :: model_pcm, model_pv, model_zero
      type(solvation_model_component_cpcm) :: pcm_component
      type(cavity_type_iswig) :: cavity
      type(radius_type_static) :: radius_model
      type(coupling_type), pointer :: coupling, coupling_pcm, coupling_zero
      type(response_type) :: response
      type(potential_adjoint_response_type), allocatable :: charge
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

      call new_context(ctx)
      call get_structure(mol, "MB16-43", "01")

      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call new_component_cpcm(pcm_component, ctx, epsilon=epsilon, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%cholesky))
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

      ! One coupling for all three models
      call stage_model_point_charge_energy(error, model_pv, qat_vals, mol, coupling)
      if (allocated(error)) return

      volume = model_pv%cavity%total_volume

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
      call fill_missing_with_zeros(model_pv%cavity, coupling)
      call fill_point_charge_field(model_pv%cavity, coupling, qat_vals, mol)

      allocate (gradient_pcm(3, mol%nat), source=0.0_wp)
      allocate (gradient_pv(3, mol%nat), source=0.0_wp)
      allocate (gradient_zero(3, mol%nat), source=0.0_wp)
      allocate (volume_gradient(3, mol%nat), source=0.0_wp)

      call model_pcm%prepare_gradient(coupling_pcm, err)
      call fill_point_charge_field(model_pcm%cavity, coupling_pcm, qat_vals, mol)
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
      call fill_point_charge_field(model_zero%cavity, coupling_zero, qat_vals, mol)
      call model_zero%get_gradient(coupling_zero, response, gradient_zero, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV(0) gradient failed: "//err%message)
         return
      end if
      ! The reverse-mode gradient contracts the surface adjoints directly
      call model_pv%cavity%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, "Cavity forward gradient failed: "//err%message)
         return
      end if
      if (.not. allocated(model_pv%cavity%v1_rA)) then
         call test_failed(error, "Cavity produced no per-point volume derivatives")
         return
      end if
      volume_gradient = sum(model_pv%cavity%v1_rA, dim=3)

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

   end subroutine test_general_model_pv_smoke

!> A coupling minted by one model is refused by every other model, and a model
!> without an internal isodensity cavity refuses a density
   subroutine test_general_model_guards(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err

      !> Molecular structure
      type(structure_type) :: mol
      !> Owning model and a second, independent model
      type(solvation_model_general), target :: model_a, model_b
      !> Only component of both models
      type(solvation_model_component_pv) :: pv_component
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

      !> Run context owned here and borrowed by the cavity and models
      type(moist_context_type), target :: ctx

      call new_context(ctx)
      call get_structure(mol, "MB16-43", "01")

      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call new_component_pv(pv_component, 0.5_wp)
      call build_model(model_a)
      if (allocated(error)) return
      call build_model(model_b)
      if (allocated(error)) return

      call model_a%new_coupling(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Coupling setup failed: "//err%message)
         return
      end if

      ! Both models are updated, so the ownership check is what refuses
      call model_b%prepare_energy(coupling, err)
      call check_foreign("staging")
      if (allocated(error)) return

      allocate (gradient(3, mol%nat), source=0.0_wp)
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
      call check(error, .not. model_b%updated, &
         & more="a refused density left the model marked usable")

   contains

      !> Assemble and update a PV-only general model
      !>
      !> @param[out] model Model to build
      subroutine build_model(model)
         !> Model to build
         type(solvation_model_general), intent(out) :: model

         call new_model_general(model, cavity, ctx, err)
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

   end subroutine test_general_model_guards

!> Assemble an updated general model from a CPCM component and an optional
!> PV component at the requested pressure
   subroutine build_pv_model(model, cavity, ctx, pcm_component, with_pv, pressure, mol, error)

      !> Model to build
      type(solvation_model_general), intent(out) :: model

      !> Cavity template copied into the model
      type(cavity_type_iswig), intent(in) :: cavity

      !> Run context owned by the caller
      type(moist_context_type), intent(in), target :: ctx

      !> CPCM component template
      type(solvation_model_component_cpcm), intent(in) :: pcm_component

      !> Whether to append a PV component
      logical, intent(in) :: with_pv

      !> Pressure of the PV component
      real(wp), intent(in) :: pressure

      !> Molecular structure
      type(structure_type), intent(in) :: mol

      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      type(solvation_model_component_pv) :: pv_component

      call new_model_general(model, cavity, ctx, error)
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

end module test_model_general
