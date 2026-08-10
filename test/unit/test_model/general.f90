!> Unit tests for the general list-based solvation model.
!>
!> Covers the model container itself rather than any single component: that a
!> component driven through the model reproduces the procedural result, that
!> several components sum, and that the lifecycle guards fire. Per-component
!> numerics live in the `test_model_component_*` suites.
module test_model_general
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_type, only: coupling_type, potential_type
   use moist_model_component_pcm_type, only: solver_type
   use moist_model_component_pcm_cpcm, only: cpcm, new_cpcm
   use moist_model_components, only: pv, new_pv
   use moist_model_general, only: general_solvation_model, new_general_model
   use moist_cavity_iswig, only: cavity_type_iswig
   use moist_radii, only: static_radius_type
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: build_test_cavity, make_charge_coupling

   implicit none
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
         & new_unittest("general_model_cpcm_pv", test_general_model_pv_smoke) &
         & ]

   end subroutine collect_model_general

!> Single-component coverage for the general list-based solvation model
   subroutine test_general_model_smoke(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(general_solvation_model) :: model
      type(cpcm) :: pcm_component, pcm_reference
      type(cavity_type_iswig) :: cavity
      type(static_radius_type) :: radius_model
      type(coupling_type) :: coupling
      type(potential_type) :: potential
      real(wp) :: energy, reference_energy

      real(wp), parameter :: epsilon = 32.0_wp
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]

      !> Run context owned here and borrowed by the cavity, model and components
      type(moist_context_type), target :: ctx

      call new_context(ctx)
      call get_structure(mol, "MB16-43", "01")
      call make_charge_coupling(qat_vals, coupling)

      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      call new_general_model(model, cavity, ctx, err)
      if (allocated(err)) then
         call test_failed(error, "General-model construction failed: "//err%message)
         return
      end if

      ! An accessor must refuse to run before the first update.
      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err), &
         & more="general-model energy was available before the first update")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)

      call new_cpcm(pcm_component, ctx, epsilon, solver=solver_type%cholesky, error=err)
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

      energy = 0.0_wp
      call model%get_energy(coupling, energy, err)
      if (allocated(err)) then
         call test_failed(error, "General-model energy failed: "//err%message)
         return
      end if

      ! Procedural reference on an independently updated cavity.
      call new_cpcm(pcm_reference, ctx, epsilon, solver=solver_type%cholesky, error=err)
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
      call pcm_reference%get_energy(coupling, cavity, reference_energy, err)
      if (allocated(err)) then
         call test_failed(error, "Reference CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy, reference_energy, thr=thr2, &
         & message="general model did not reproduce the procedural CPCM energy")
      if (allocated(error)) return

      ! The direct trace channel must carry the component surface charges.
      call model%get_trace_potential(coupling, potential, err)
      if (allocated(err)) then
         call test_failed(error, "General-model trace potential failed: "//err%message)
         return
      end if
      call check(error, allocated(potential%w_elstat_umol), &
         & more="general model did not expose CPCM surface charges")
      if (allocated(error)) return
      call check(error, maxval(abs(potential%w_elstat_umol - pcm_reference%q)), 0.0_wp, &
         & thr=thr2, &
         & message="general-model CPCM charges differ from the procedural reference")
      if (allocated(error)) return

      ! Components are frozen once the model has been updated.
      call model%add_component(pcm_component, err)
      call check(error, allocated(err), &
         & more="adding a component after the first update was not rejected")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)

      ! A failed update must invalidate a model that was previously usable.
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
      type(general_solvation_model) :: model_pcm, model_pv, model_zero
      type(cpcm) :: pcm_component
      type(cavity_type_iswig) :: cavity
      type(static_radius_type) :: radius_model
      type(coupling_type) :: coupling
      type(potential_type) :: potential
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
      call make_charge_coupling(qat_vals, coupling)

      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      call new_cpcm(pcm_component, ctx, epsilon, solver=solver_type%cholesky, error=err)
      if (allocated(err)) then
         call test_failed(error, "CPCM construction failed: "//err%message)
         return
      end if

      ! Reference model: CPCM only.
      call build_pv_model(model_pcm, cavity, ctx, pcm_component, .false., 0.0_wp, mol, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM-only model setup failed: "//err%message)
         return
      end if

      ! Two-component model: CPCM + PV at finite pressure.
      call build_pv_model(model_pv, cavity, ctx, pcm_component, .true., pressure, mol, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV model setup failed: "//err%message)
         return
      end if

      ! Two-component model with a vanishing pressure.
      call build_pv_model(model_zero, cavity, ctx, pcm_component, .true., 0.0_wp, mol, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV(0) model setup failed: "//err%message)
         return
      end if

      volume = model_pv%cavity%total_volume

      energy_pcm = 0.0_wp
      call model_pcm%get_energy(coupling, energy_pcm, err)
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
      call model_zero%get_energy(coupling, energy_zero, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV(0) energy failed: "//err%message)
         return
      end if

      ! The component energies must simply add up.
      call check(error, energy_pv, energy_pcm + pressure*volume, thr=thr2, &
         & message="CPCM+PV model energy is not the sum of its components")
      if (allocated(error)) return

      ! A vanishing pressure must leave the CPCM energy untouched.
      call check(error, energy_zero, energy_pcm, thr=thr, &
         & message="PV at zero pressure changed the model energy")
      if (allocated(error)) return

      ! The shared surface-adjoint accumulator must survive two components.
      call model_pv%get_potential(coupling, potential, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV potential failed: "//err%message)
         return
      end if
      call check(error, allocated(potential%w_elstat_umol), &
         & more="CPCM+PV model produced no electrostatic potential channel")
      if (allocated(error)) return

      ! The direct trace channel must still carry the CPCM surface charges when
      ! a second, non-electrostatic component shares the accumulator.
      call model_pv%get_trace_potential(coupling, potential, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV trace potential failed: "//err%message)
         return
      end if
      call check(error, allocated(potential%w_elstat_umol), &
         & more="CPCM+PV model did not expose CPCM surface charges")
      if (allocated(error)) return

      allocate (gradient_pcm(3, mol%nat), source=0.0_wp)
      allocate (gradient_pv(3, mol%nat), source=0.0_wp)
      allocate (gradient_zero(3, mol%nat), source=0.0_wp)
      allocate (volume_gradient(3, mol%nat), source=0.0_wp)

      call model_pcm%get_gradient(coupling, gradient_pcm, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM-only gradient failed: "//err%message)
         return
      end if
      call model_pv%get_gradient(coupling, gradient_pv, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV gradient failed: "//err%message)
         return
      end if
      call model_zero%get_gradient(coupling, gradient_zero, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM+PV(0) gradient failed: "//err%message)
         return
      end if
      if (.not. allocated(model_pv%cavity%v1_rA)) then
         call test_failed(error, "Cavity produced no per-point volume derivatives")
         return
      end if
      volume_gradient = sum(model_pv%cavity%v1_rA, dim=3)

      ! The PV gradient is the pressure-scaled cavity-volume gradient.
      call check(error, maxval(abs(gradient_pv - gradient_pcm - pressure*volume_gradient)), &
         & 0.0_wp, thr=thr2, &
         & message="CPCM+PV gradient is not the CPCM gradient plus p*dV/dR")
      if (allocated(error)) return

      ! A vanishing pressure must leave the CPCM gradient untouched.
      call check(error, maxval(abs(gradient_zero - gradient_pcm)), 0.0_wp, thr=thr, &
         & message="PV at zero pressure changed the model gradient")

   end subroutine test_general_model_pv_smoke

!> Assemble an updated general model from a CPCM component and an optional
!> PV component at the requested pressure.
   subroutine build_pv_model(model, cavity, ctx, pcm_component, with_pv, pressure, mol, error)

      !> Model to build
      type(general_solvation_model), intent(out) :: model

      !> Cavity template copied into the model
      type(cavity_type_iswig), intent(in) :: cavity

      !> Run context owned by the caller
      type(moist_context_type), intent(in), target :: ctx

      !> CPCM component template
      type(cpcm), intent(in) :: pcm_component

      !> Whether to append a PV component
      logical, intent(in) :: with_pv

      !> Pressure of the PV component
      real(wp), intent(in) :: pressure

      !> Molecular structure
      type(structure_type), intent(in) :: mol

      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      type(pv) :: pv_component

      call new_general_model(model, cavity, ctx, error)
      if (allocated(error)) return
      call model%add_component(pcm_component, error)
      if (allocated(error)) return
      if (with_pv) then
         call new_pv(pv_component, pressure)
         call model%add_component(pv_component, error)
         if (allocated(error)) return
      end if
      call model%update(mol, error)

   end subroutine build_pv_model

end module test_model_general
