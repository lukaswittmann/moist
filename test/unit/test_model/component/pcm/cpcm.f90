!> Unit tests for the CPCM solvation model component
module test_model_component_pcm_cpcm
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use moist_type, only: coupling_type
   use moist_model_component_pcm_type, only: solver_type, potential_source
   use moist_model_component_pcm_cpcm, only: cpcm, new_cpcm
   use moist_model_component_pcm_cosmo, only: cosmo, new_cosmo
   use moist_model_component_pcm_solvers, only: solve_pcm_lu
   use moist_model_component_pcm_amat, only: assemble_pcm_amat
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_cavity_drop, only: cavity_type_drop
   use moist_radii, only: static_radius_type, new_cosmo_radii, &
      & new_radii_custom_atoms, radius_type
   use moist_context, only: moist_context_type, new_context
   use test_model_component_helper, only: surface_fixture, &
      & new_surface_fixture, check_surface_weights, &
      & ngrid_sw => fixture_ngrid_param, sw_areas => fixture_areas_param, &
      & sw_xis => fixture_xis_param, sw_fs => fixture_fs_param, &
      & sw_xyz => fixture_xyz_param, sw_normals => fixture_normals_param
   use test_helpers, only: get_test_structures, center_at_origin, &
      & fd4_scalar, fd4_offsets, get_test_cavity_iswig, build_test_cavity, &
      & make_charge_coupling
   implicit none (type, external)
   private

   public :: collect_model_component_pcm_cpcm

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   !> Host surface potential driving the surface-weight test
   real(wp), parameter :: sw_phi(ngrid_sw) = [-0.31_wp, 0.22_wp, -0.17_wp, 0.41_wp, &
                                              -0.28_wp, 0.13_wp, -0.36_wp]

contains

!> Collect all exported unit tests
   subroutine collect_model_component_pcm_cpcm(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("cpcm_born_ion", test_cpcm_born), &
         & new_unittest("cpcm_lu_energy", test_cpcm_energy_lu), &
         & new_unittest("cpcm_charged_system", test_cpcm_charged), &
         & new_unittest("cpcm_solver_comparison", test_cpcm_solver_comparison), &
         & new_unittest("cpcm_vacuum_limit", test_cpcm_vacuum_limit), &
         & new_unittest("cpcm_external_potential", test_cpcm_external_potential), &
         & new_unittest("cpcm_external_potential_requires_input", test_cpcm_external_potential_requires_input, &
            should_fail=.true.), &
         & new_unittest("cpcm_external_matrix", test_cpcm_external_matrix), &
         & new_unittest("cpcm_spin_resolved_charges", test_cpcm_spin_resolved_charges), &
         & new_unittest("cpcm_requires_update", test_cpcm_requires_update, should_fail=.true.), &
         & new_unittest("cpcm_invalid_solver", test_cpcm_invalid_solver, should_fail=.true.), &
         & new_unittest("pcm_rejects_dielectric_below_one", test_pcm_invalid_epsilon), &
         & new_unittest("cpcm_iterative_rejects_non_spd_matrix", test_cpcm_iterative_not_spd, should_fail=.true.), &
         & new_unittest("cpcm_reallocates_on_grid_change", test_cpcm_reallocate_on_ngrid_change), &
         & new_unittest("cpcm_records_shared_timer_tree", test_cpcm_timer_tree), &
         & new_unittest("cpcm_stale_charge_regression", test_cpcm_stale_charge_regression), &
         & new_unittest("cpcm_surface_weights", test_cpcm_surface_weights), &
         & new_unittest("cpcm_molecular_surface_weights", test_cpcm_molecular_surface_weights), &
         & new_unittest("cpcm_nuclear_gradient", test_cpcm_nuclear_gradient), &
         & new_unittest("cpcm_external_nuclear_gradient", test_cpcm_external_nuclear_gradient), &
         & new_unittest("cpcm_translation_invariance", test_cpcm_translation_invariance), &
         & new_unittest("cpcm_dielectric_scaling", test_cpcm_dielectric_scaling), &
         & new_unittest("cpcm_coincident_points", test_cpcm_coincident_points) &
         ! & new_unittest("cpcm_solver_timing", test_cpcm_timing) &
         & ]

   end subroutine collect_model_component_pcm_cpcm

   !> Test CPCM against the analytic Born ion
   subroutine test_cpcm_born(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      class(radius_type), allocatable :: radius_model

      !> Born sphere radius (Bohr)
      real(wp), parameter :: rad = 4.0_wp
      !> Dielectric constant
      real(wp), parameter :: epsilon = 78.4_wp
      !> Lebedev grid sizes (coarse and fine, both iSwiG-supported)
      integer, parameter :: nlebs(2) = [302, 1202]

      real(wp) :: xyz(3, 1), energy_array, e_ref, errs(2)
      integer :: ileb
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      xyz(:, 1) = 0.0_wp
      call new(mol, [1], xyz)
      call make_charge_coupling([1.0_wp], coupling)

      call new_radii_custom_atoms([rad], radius_model, err)
      if (allocated(err)) then
         call test_failed(error, "Radius model setup failed: "//err%message)
         return
      end if

      ! Analytic CPCM Born energy: E = -1/2 * (eps-1)/eps * q^2 / R
      e_ref = -0.5_wp*(epsilon - 1.0_wp)/epsilon/rad

      do ileb = 1, size(nlebs)
         call new_cavity_iswig(cavity, ctx, nleb=nlebs(ileb), &
            & radius_model=radius_model, error=err)
         if (allocated(err)) then
            call test_failed(error, "Cavity initialization failed: "//err%message)
            return
         end if
         call cavity%update(mol, error=err)
         if (allocated(err)) then
            call test_failed(error, "Cavity update failed: "//err%message)
            return
         end if

         call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%cholesky, error=err)
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
         call pcm_model%get_energy(coupling, cavity, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed: "//err%message)
            return
         end if

         errs(ileb) = abs(energy_array - e_ref)
      end do

      ! Both grids must reproduce the analytic Born energy
      do ileb = 1, size(nlebs)
         call check(error, errs(ileb)/abs(e_ref), 0.0_wp, thr=1.0E-6_wp, &
            & message="Born ion energy deviates from analytic CPCM result")
         if (allocated(error)) return
      end do

   end subroutine test_cpcm_born

!> Test CPCM energy calculation with all solvers on neutral system
   subroutine test_cpcm_energy_lu(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), parameter :: qat(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp), parameter :: epsilon = 78.4_wp
      real(wp), parameter :: ref_energy = -1.125110188578178E-2_wp

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
      real(wp), parameter :: ref_energy = -8.994600196814041E-2_wp

      call get_structure(mol, "MB16-43", "01")
      call test_all_solvers(error, mol, qat, epsilon, ref_energy, "charged system")

   end subroutine test_cpcm_charged

!> Test that all solvers give identical results
   subroutine test_cpcm_solver_comparison(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp), parameter :: epsilon = 78.4_wp
      integer :: solvers(4)
      character(len=20) :: solver_names(4)
      real(wp) :: energies(4)
      real(wp), allocatable :: charges(:, :)
      real(wp) :: energy_array
      real(wp) :: energy_diff, charge_rms
      integer :: i, j
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      ! Setup solver types and names
      solvers = [solver_type%inversion, solver_type%lu, &
                 solver_type%cholesky, solver_type%iterative]
      solver_names = ["inversion ", "lu        ", &
                      "cholesky  ", "iterative "]

      ! Get test molecule
      call get_structure(mol, "MB16-43", "01")

      ! Build cavity
      call new_cosmo_radii(radius_model)
      call new_cavity_iswig(cavity, ctx, nleb=50, radius_model=radius_model, error=err)
      if (allocated(err)) then
         call test_failed(error, "Cavity initialization failed: "//err%message)
         return
      end if
      call cavity%update(mol, error=err)
      if (allocated(err)) then
         call test_failed(error, "Cavity update failed: "//err%message)
         return
      end if

      ! Setup wavefunction
      coupling%qat = reshape(qat_vals, [size(qat_vals), 1])

      ! Allocate storage for all solver results
      allocate (charges(cavity%ngrid, 4))

      ! Test all 4 solvers
      do i = 1, 4
         call new_cpcm(pcm_model, ctx, epsilon, solver=solvers(i), error=err)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed ("//trim(solver_names(i))//")")
            return
         end if

         ! Tighten CG tolerance
         pcm_model%solver_tol = 1.0e-14_wp
         pcm_model%solver_maxiter = 10000

         call pcm_model%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed ("//trim(solver_names(i))//"): "//err%message)
            return
         end if
         energy_array = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed ("//trim(solver_names(i))//"): "//err%message)
            return
         end if
         energies(i) = energy_array
         charges(:, i) = pcm_model%q
      end do

      ! Compare all solver pairs (using first solver as reference)
      do i = 2, 4
         ! Compare energies to first solver (inversion)
         call check(error, energies(i), energies(1), thr=thr, &
                    message=trim(solver_names(i))//" energy differs from "// &
                    trim(solver_names(1)))
         if (allocated(error)) return

         ! Compare charges to first solver (inversion)
         charge_rms = sqrt(sum((charges(:, i) - charges(:, 1))**2)/real(cavity%ngrid, wp))
         call check(error, charge_rms, 0.0_wp, thr=thr2, &
                    message=trim(solver_names(i))//" charges differ from "// &
                    trim(solver_names(1)))
         if (allocated(error)) return
      end do

      deallocate (charges)

   end subroutine test_cpcm_solver_comparison

!> Test CPCM vacuum limit, where f(epsilon)=0 should give zero charges and energy.
   subroutine test_cpcm_vacuum_limit(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      integer :: solvers(4)
      character(len=20) :: solver_names(4)
      real(wp) :: energy_array
      integer :: i

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      solvers = [solver_type%inversion, solver_type%lu, &
                 solver_type%cholesky, solver_type%iterative]
      solver_names = ["inversion ", "lu        ", &
                      "cholesky  ", "iterative "]

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_coupling(qat_vals, coupling)

      do i = 1, 4
         call new_cpcm(pcm_model, ctx, 1.0_wp, solver=solvers(i), error=err)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed ("//trim(solver_names(i))//")")
            return
         end if

         call pcm_model%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed ("//trim(solver_names(i))//"): "//err%message)
            return
         end if

         energy_array = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed ("//trim(solver_names(i))//"): "//err%message)
            return
         end if

         call check(error, energy_array, 0.0_wp, thr=thr, &
            & message=trim(solver_names(i))//" vacuum energy should vanish")
         if (allocated(error)) return

         call check(error, maxval(abs(pcm_model%q)), 0.0_wp, thr=thr, &
            & message=trim(solver_names(i))//" vacuum surface charges should vanish")
         if (allocated(error)) return
      end do

   end subroutine test_cpcm_vacuum_limit

!> Test that externally supplied potentials reproduce the internally computed result.
   subroutine test_cpcm_external_potential(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_internal, pcm_external
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_internal, energy_external, charge_rms
      real(wp), allocatable :: phi_ref(:), q_ref(:)

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_coupling(qat_vals, coupling)

      call new_cpcm(pcm_internal, ctx, 78.4_wp, solver=solver_type%cholesky, error=err)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_internal%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM update failed: "//err%message)
         return
      end if

      energy_internal = 0.0_wp
      call pcm_internal%get_energy(coupling, cavity, energy_internal, err)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM energy failed: "//err%message)
         return
      end if

      allocate (phi_ref, source=pcm_internal%phi)
      allocate (q_ref, source=pcm_internal%q)

      call new_cpcm(pcm_external, ctx, 78.4_wp, solver=solver_type%cholesky, &
         & phi_source=potential_source%external, error=err)
      if (allocated(err)) then
         call test_failed(error, "External-potential CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_external%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "External-potential CPCM update failed: "//err%message)
         return
      end if

      call pcm_external%input_potential(phi_ref)
      energy_external = 0.0_wp
      call pcm_external%get_energy(coupling, cavity, energy_external, err)
      if (allocated(err)) then
         call test_failed(error, "External-potential CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy_external, energy_internal, thr=thr, &
         & message="External potential energy differs from internally computed potential")
      if (allocated(error)) return

      charge_rms = sqrt(sum((pcm_external%q - q_ref)**2)/real(size(q_ref), wp))
      call check(error, charge_rms, 0.0_wp, thr=thr2, &
         & message="External potential surface charges differ from internal reference")
      if (allocated(error)) return

   end subroutine test_cpcm_external_potential

!> Test that external potential mode requires input_potential() before get_energy().
   subroutine test_cpcm_external_potential_requires_input(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_array

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_coupling(qat_vals, coupling)

      call new_cpcm(pcm_model, ctx, 78.4_wp, solver=solver_type%cholesky, &
         & phi_source=potential_source%external, error=err)
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
      call pcm_model%get_energy(coupling, cavity, energy_array, err)
      if (.not. allocated(err)) return

      call test_failed(error, "External potential without input was correctly rejected: "//err%message)

   end subroutine test_cpcm_external_potential_requires_input

!> Test that supplying an external matrix reproduces the cavity-assembled result.
   subroutine test_cpcm_external_matrix(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_internal, pcm_external
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_internal, energy_external, charge_rms
      real(wp), allocatable :: amat_ref(:, :), q_ref(:)

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_coupling(qat_vals, coupling)

      call new_cpcm(pcm_internal, ctx, 78.4_wp, solver=solver_type%lu, error=err)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_internal%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM update failed: "//err%message)
         return
      end if

      energy_internal = 0.0_wp
      call pcm_internal%get_energy(coupling, cavity, energy_internal, err)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM energy failed: "//err%message)
         return
      end if

      allocate (amat_ref, source=pcm_internal%amat)
      allocate (q_ref, source=pcm_internal%q)

      call new_cpcm(pcm_external, ctx, 78.4_wp, solver=solver_type%lu, &
         & external_matrix=amat_ref, error=err)
      if (allocated(err)) then
         call test_failed(error, "External-matrix CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_external%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "External-matrix CPCM update failed: "//err%message)
         return
      end if

      energy_external = 0.0_wp
      call pcm_external%get_energy(coupling, cavity, energy_external, err)
      if (allocated(err)) then
         call test_failed(error, "External-matrix CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy_external, energy_internal, thr=thr, &
         & message="External matrix energy differs from cavity-assembled matrix")
      if (allocated(error)) return

      charge_rms = sqrt(sum((pcm_external%q - q_ref)**2)/real(size(q_ref), wp))
      call check(error, charge_rms, 0.0_wp, thr=thr2, &
         & message="External matrix surface charges differ from cavity-assembled matrix")
      if (allocated(error)) return

   end subroutine test_cpcm_external_matrix

!> Test that spin-resolved charges are summed consistently when building the potential.
   subroutine test_cpcm_spin_resolved_charges(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_scalar, pcm_spin
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling_scalar, coupling_spin
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_scalar, energy_spin, charge_rms

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      call make_charge_coupling(qat_vals, coupling_scalar)
      allocate (coupling_spin%qat(size(qat_vals), 2))
      coupling_spin%qat(:, 1) = 0.25_wp*qat_vals
      coupling_spin%qat(:, 2) = 0.75_wp*qat_vals

      call new_cpcm(pcm_scalar, ctx, 78.4_wp, solver=solver_type%cholesky, error=err)
      if (allocated(err)) then
         call test_failed(error, "Scalar-charge CPCM initialization failed: "//err%message)
         return
      end if
      call pcm_scalar%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Scalar-charge CPCM update failed: "//err%message)
         return
      end if
      energy_scalar = 0.0_wp
      call pcm_scalar%get_energy(coupling_scalar, cavity, energy_scalar, err)
      if (allocated(err)) then
         call test_failed(error, "Scalar-charge CPCM energy failed: "//err%message)
         return
      end if

      call new_cpcm(pcm_spin, ctx, 78.4_wp, solver=solver_type%cholesky, error=err)
      if (allocated(err)) then
         call test_failed(error, "Spin-charge CPCM initialization failed: "//err%message)
         return
      end if
      call pcm_spin%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Spin-charge CPCM update failed: "//err%message)
         return
      end if
      energy_spin = 0.0_wp
      call pcm_spin%get_energy(coupling_spin, cavity, energy_spin, err)
      if (allocated(err)) then
         call test_failed(error, "Spin-charge CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy_spin, energy_scalar, thr=thr, &
         & message="Spin-resolved charges should reproduce scalar-charge energy")
      if (allocated(error)) return

      charge_rms = sqrt(sum((pcm_spin%q - pcm_scalar%q)**2)/real(size(pcm_scalar%q), wp))
      call check(error, charge_rms, 0.0_wp, thr=thr2, &
         & message="Spin-resolved charges should reproduce scalar-charge surface charges")
      if (allocated(error)) return

   end subroutine test_cpcm_spin_resolved_charges

!> Test that get_energy() reports an error when update() has not been called.
   subroutine test_cpcm_requires_update(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(cpcm) :: pcm_model
      type(coupling_type) :: coupling
      !> Never updated; only its (empty) grid data reach the component
      type(cavity_type_iswig) :: cavity
      real(wp), parameter :: qat_vals(1) = [0.0_wp]
      real(wp) :: energy_array

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call make_charge_coupling(qat_vals, coupling)
      call new_cpcm(pcm_model, ctx, 78.4_wp, error=err)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if

      energy_array = 0.0_wp
      call pcm_model%get_energy(coupling, cavity, energy_array, err)
      if (.not. allocated(err)) return

      call test_failed(error, "Energy before update was correctly rejected: "//err%message)

   end subroutine test_cpcm_requires_update

!> Test that invalid solver identifiers are rejected during the solve step.
   subroutine test_cpcm_invalid_solver(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_array

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_coupling(qat_vals, coupling)

      call new_cpcm(pcm_model, ctx, 78.4_wp, solver=-1, error=err)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization unexpectedly failed: "//err%message)
         return
      end if

      call pcm_model%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if

      energy_array = 0.0_wp
      call pcm_model%get_energy(coupling, cavity, energy_array, err)
      if (.not. allocated(err)) return

      if (index(err%message, "Unknown solver type") == 0) return

      call test_failed(error, "Invalid solver correctly rejected: "//trim(err%message))

   end subroutine test_cpcm_invalid_solver

   !> Test that the CG solver reports an error for a non-SPD system matrix
   subroutine test_cpcm_iterative_not_spd(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp), allocatable :: bad_amat(:, :)
      real(wp) :: energy_array
      integer :: i

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_coupling(qat_vals, coupling)

      ! Negative-definite external matrix: -identity
      allocate (bad_amat(cavity%ngrid, cavity%ngrid), source=0.0_wp)
      do i = 1, cavity%ngrid
         bad_amat(i, i) = -1.0_wp
      end do

      call new_cpcm(pcm_model, ctx, 78.4_wp, solver=solver_type%iterative, &
         & external_matrix=bad_amat, error=err)
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
      call pcm_model%get_energy(coupling, cavity, energy_array, err)
      if (.not. allocated(err)) return

      call test_failed(error, "Non-SPD matrix was correctly rejected: "//err%message)

   end subroutine test_cpcm_iterative_not_spd

!> Test that a PCM instance can be reused safely when the cavity grid size changes.
   subroutine test_cpcm_reallocate_on_ngrid_change(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_reused, pcm_fresh
      type(cavity_type_iswig) :: cavity_small, cavity_large
      type(coupling_type) :: coupling
      type(static_radius_type) :: radius_small, radius_large
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_reused, energy_fresh, charge_rms

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

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
      call make_charge_coupling(qat_vals, coupling)

      call new_cpcm(pcm_reused, ctx, 78.4_wp, solver=solver_type%lu, error=err)
      if (allocated(err)) then
         call test_failed(error, "Reused CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_reused%update(mol, cavity_small, err)
      if (allocated(err)) then
         call test_failed(error, "Update with small cavity failed: "//err%message)
         return
      end if
      energy_reused = 0.0_wp
      call pcm_reused%get_energy(coupling, cavity_small, energy_reused, err)
      if (allocated(err)) then
         call test_failed(error, "Energy with small cavity failed: "//err%message)
         return
      end if

      call pcm_reused%update(mol, cavity_large, err)
      if (allocated(err)) then
         call test_failed(error, "Update with large cavity failed after reuse: "//err%message)
         return
      end if
      energy_reused = 0.0_wp
      call pcm_reused%get_energy(coupling, cavity_large, energy_reused, err)
      if (allocated(err)) then
         call test_failed(error, "Energy with reused large cavity failed: "//err%message)
         return
      end if

      if (size(pcm_reused%q) /= cavity_large%ngrid) then
         call test_failed(error, "Surface charge array was not resized to the new grid")
         return
      end if

      if (size(pcm_reused%amat, 1) /= cavity_large%ngrid .or. &
          & size(pcm_reused%amat, 2) /= cavity_large%ngrid) then
         call test_failed(error, "PCM matrix was not resized to the new grid")
         return
      end if

      call new_cpcm(pcm_fresh, ctx, 78.4_wp, solver=solver_type%lu, error=err)
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
      call pcm_fresh%get_energy(coupling, cavity_large, energy_fresh, err)
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

!> Test that the component records its work on the run context's shared timer,
!> nested under whatever the caller had open. This is the observable proof that
!> the component borrows the same context as the cavity: without it the PCM
!> nodes would either be absent or sit at top level.
   subroutine test_cpcm_timer_tree(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_array
      integer :: i
      logical :: found_setup, found_energy, found_solve
      character(len=:), allocatable :: name

      !> Run context owned here and borrowed by the cavity and component
      type(moist_context_type), target :: ctx

      call new_context(ctx, verbosity=2)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_coupling(qat_vals, coupling)

      call new_cpcm(pcm_model, ctx, 78.4_wp, solver=solver_type%cholesky, error=err)
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

      energy_array = 0.0_wp
      call pcm_model%get_energy(coupling, cavity, energy_array, err)
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

!> Regression: changing the atomic charges without an intervening update() must
!> re-solve for the surface charges.
!>
!> This is the ordinary SCF pattern - the geometry is fixed, so the host calls
!> update() once and then feeds new `coupling%qat` on every cycle.
!> [[pcm_ensure_charges]] used to recompute phi from those charges but leave
!> `charges_valid` untouched, so the second cycle paired a fresh phi with the
!> stale q and the energy came out linear in q instead of quadratic.
!>
!> The caching half is tested too: an unchanged `coupling%qat` must still
!> short-circuit before the linear solve, otherwise every accessor call - and
!> the general model issues two per get_potential - pays for another solve.
   subroutine test_cpcm_stale_charge_regression(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_model, pcm_fresh
      type(cavity_type_iswig) :: cavity
      type(static_radius_type) :: radius_model
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

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 14, ctx, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%cholesky, error=err)
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
      call make_charge_coupling(qat_vals, coupling)
      energy_first = 0.0_wp
      call pcm_model%get_energy(coupling, cavity, energy_first, err)
      if (allocated(err)) then
         call test_failed(error, "First cycle energy failed: "//err%message)
         return
      end if

      ! Repeating the accessor with unchanged charges must return the very same
      ! number, bit for bit.
      energy_repeat = 0.0_wp
      call pcm_model%get_energy(coupling, cavity, energy_repeat, err)
      if (allocated(err)) then
         call test_failed(error, "Repeated cycle energy failed: "//err%message)
         return
      end if
      call check(error, energy_repeat, energy_first, thr=0.0_wp, &
         & message="repeating an accessor with unchanged charges changed the energy")
      if (allocated(error)) return

      ! ... and it must do so without re-solving. Poisoning the cached charges
      ! makes the solve count observable: a short circuit keeps the poisoned
      ! charges (energy doubles with them), a fresh solve would overwrite them
      ! and hand back the unpoisoned energy.
      pcm_model%q(:) = 2.0_wp*pcm_model%q(:)
      energy_cached = 0.0_wp
      call pcm_model%get_energy(coupling, cavity, energy_cached, err)
      if (allocated(err)) then
         call test_failed(error, "Cached cycle energy failed: "//err%message)
         return
      end if
      call check(error, energy_cached, 2.0_wp*energy_first, thr=thr, &
         & message="unchanged atomic charges triggered a fresh PCM solve")
      if (allocated(error)) return
      pcm_model%q(:) = 0.5_wp*pcm_model%q(:)

      ! Cycle 2: new charges, same geometry, deliberately no update() in between
      call make_charge_coupling(scale*qat_vals, coupling)
      energy_scaled = 0.0_wp
      call pcm_model%get_energy(coupling, cavity, energy_scaled, err)
      if (allocated(err)) then
         call test_failed(error, "Second cycle energy failed: "//err%message)
         return
      end if

      ! E = 0.5 * q^T A^-1 q is quadratic in the atomic charges at fixed geometry
      call check(error, energy_scaled, scale**2*energy_first, thr=thr2, &
         & message="PCM energy is not quadratic in the atomic charges - "// &
         & "surface charges were not re-solved after the charges changed")
      if (allocated(error)) return

      ! Same statement against an independently constructed model
      call new_cpcm(pcm_fresh, ctx, epsilon, solver=solver_type%cholesky, error=err)
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
      call pcm_fresh%get_energy(coupling, cavity, energy_reference, err)
      if (allocated(err)) then
         call test_failed(error, "Fresh CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy_scaled, energy_reference, thr=thr, &
         & message="reused CPCM model does not match a freshly built one "// &
         & "after the atomic charges changed")
      if (allocated(error)) return

      call check(error, maxval(abs(pcm_model%q - pcm_fresh%q)), 0.0_wp, thr=thr, &
         & message="reused CPCM model kept stale surface charges")
      if (allocated(error)) return

   end subroutine test_cpcm_stale_charge_regression

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

      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      real(wp) :: energy_array
      real(wp) :: energy_lu, energy_cholesky, energy_iterative, energy_inversion
      type(static_radius_type) :: radius_model
      integer :: solvers(4)
      character(len=20) :: solver_names(4)
      real(wp) :: energies(4)
      real(wp) :: energy_diff
      integer :: i
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      ! Setup solver types and names
      solvers = [solver_type%inversion, solver_type%lu, &
                 solver_type%cholesky, solver_type%iterative]
      solver_names = ["inversion ", "lu        ", &
                      "cholesky  ", "iterative "]

      ! Build cavity and wavefunction
      call new_cosmo_radii(radius_model)
      call new_cavity_iswig(cavity, ctx, nleb=50, radius_model=radius_model, error=err)
      if (allocated(err)) then
         call test_failed(error, "Cavity initialization failed: "//err%message)
         return
      end if
      call cavity%update(mol, error=err)
      if (allocated(err)) then
         call test_failed(error, "Cavity update failed: "//err%message)
         return
      end if
      coupling%qat = reshape(qat, [size(qat), 1])

      ! Test all 4 solvers
      do i = 1, 4
         call new_cpcm(pcm_model, ctx, epsilon, solver=solvers(i), error=err)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed for "// &
                             trim(solver_names(i))//" solver")
            return
         end if

         ! Set CG tolerance
         pcm_model%solver_tol = 1.0e-14_wp
         pcm_model%solver_maxiter = 10000

         call pcm_model%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed for "//trim(solver_names(i))//": "//err%message)
            return
         end if
         energy_array = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy_array, err)
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

!> Benchmark timing for all PCM solvers on polyalanine structures
   subroutine test_cpcm_timing(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      real(wp) :: energy_array
      real(wp), allocatable :: qat(:)
      type(static_radius_type) :: radius_model
      real(wp), parameter :: epsilon = 78.4_wp
      integer :: n_ala, ii, i
      integer(8) :: t1, t2, rate
      real(wp) :: time_lu, time_cholesky, time_iterative, time_inversion
      real(wp) :: energy_lu, energy_cholesky, energy_iterative, energy_inversion
      real(wp) :: dE_lu, dE_cholesky, dE_iterv
      real(wp) :: rms_lu, rms_cholesky, rms_iterv
      real(wp), allocatable :: charges_lu(:), charges_cholesky(:)
      real(wp), allocatable :: charges_iterative(:), charges_inversion(:)
      character(len=10) :: n_str
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      ! Print header
      print '(a)', ""
      print '(a)', "CPCM Solver Timing Benchmark (Polyalanine)"
      print '(a)', "==========================================="
      print '(a7, a10, 4(a11), 3(a11), 3(a11))', "N_Ala", "N_Grid", &
         "t_inv (s)", "t_lu (s)", "t_chol (s)", "t_iter (s)", &
         "dE_lu", "dE_chol", "dE_iter", &
         "rms_lu", "rms_chol", "rms_iter"
      print '(a7, a10, 10(a11))', "------", "---------", &
         "----------", "----------", "----------", "----------", &
         "----------", "----------", "----------", &
         "----------", "----------", "----------"

      call system_clock(count_rate=rate)

      ! Loop over polyalanine structures: 4, 8, 12, ..., 100
      ! Use steps that give us structures around 50-100
      do n_ala = 4, 100, 4

         ! Get structure from mstore
         write (n_str, '(a, i2.2)') 'polyala_', n_ala
         call get_structure(mol, "POLYALANINE", trim(n_str))

         ! Prepare charges (simple uniform distribution)
         allocate (qat(mol%nat))
         qat = 0.0_wp
         do ii = 1, mol%nat
            qat(ii) = 0.3_wp*sin(real(ii, wp))
         end do
         coupling%qat = reshape(qat, [size(qat), 1])

         ! Build cavity
         call new_cosmo_radii(radius_model)
         call new_cavity_iswig(cavity, ctx, nleb=50, radius_model=radius_model, error=err)
         if (allocated(err)) then
            call test_failed(error, "Cavity initialization failed: "//err%message)
            return
         end if
         call cavity%update(mol, error=err)
         if (allocated(err)) then
            call test_failed(error, "Cavity update failed: "//err%message)
            return
         end if

         ! ===== Time LU solver (reference) =====
         call system_clock(t1)
         call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%lu, error=err)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed (lu)")
            return
         end if
         call pcm_model%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed (lu): "//err%message)
            return
         end if
         energy_array = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed (lu): "//err%message)
            return
         end if
         energy_lu = energy_array
         allocate (charges_lu(cavity%ngrid))
         charges_lu = pcm_model%q
         call system_clock(t2)
         time_lu = real(t2 - t1, wp)/real(rate, wp)

         ! ===== Time Cholesky solver =====
         call system_clock(t1)
         call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%cholesky, error=err)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed (cholesky)")
            return
         end if
         call pcm_model%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed (cholesky): "//err%message)
            return
         end if
         energy_array = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed (cholesky): "//err%message)
            return
         end if
         energy_cholesky = energy_array
         allocate (charges_cholesky(cavity%ngrid))
         charges_cholesky = pcm_model%q
         call system_clock(t2)
         time_cholesky = real(t2 - t1, wp)/real(rate, wp)

         ! ===== Time iterative solver =====
         call system_clock(t1)
         call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%iterative, error=err)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed (iterative)")
            return
         end if
         call pcm_model%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed (iterative): "//err%message)
            return
         end if
         energy_array = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed (iterative): "//err%message)
            return
         end if
         energy_iterative = energy_array
         allocate (charges_iterative(cavity%ngrid))
         charges_iterative = pcm_model%q
         call system_clock(t2)
         time_iterative = real(t2 - t1, wp)/real(rate, wp)

         ! ===== Time inversion solver =====
         call system_clock(t1)
         call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%inversion, error=err)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed (inversion)")
            return
         end if
         call pcm_model%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed (inversion): "//err%message)
            return
         end if
         energy_array = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed (inversion): "//err%message)
            return
         end if
         energy_inversion = energy_array
         allocate (charges_inversion(cavity%ngrid))
         charges_inversion = pcm_model%q
         call system_clock(t2)
         time_inversion = real(t2 - t1, wp)/real(rate, wp)

         ! Compute comparison metrics (all relative to inversion solver as reference)
         dE_lu = abs(energy_lu - energy_inversion)
         dE_cholesky = abs(energy_cholesky - energy_inversion)
         dE_iterv = abs(energy_iterative - energy_inversion)

         rms_lu = sqrt(sum((charges_lu - charges_inversion)**2)/real(cavity%ngrid, wp))
         rms_cholesky = sqrt(sum((charges_cholesky - charges_inversion)**2)/real(cavity%ngrid, wp))
         rms_iterv = sqrt(sum((charges_iterative - charges_inversion)**2)/real(cavity%ngrid, wp))

         ! Print results: N_Ala, N_Grid, times (invers, lu, cholesky, iterv),
         !                dE (lu, cholesky, iterv), rms (lu, cholesky, iterv)
         print '(i7, i10, 4(f11.3), 3(e11.2), 3(e11.2))', n_ala, cavity%ngrid, &
            time_inversion, time_lu, time_cholesky, time_iterative, &
            dE_lu, dE_cholesky, dE_iterv, &
            rms_lu, rms_cholesky, rms_iterv

         ! Cleanup for next iteration
         deallocate (qat, charges_lu, charges_cholesky, &
                     charges_iterative, charges_inversion)

      end do

      print '(a)', "==========================================="
      print '(a)', ""

   end subroutine test_cpcm_timing


!* ================================================================================= *!
!*                    Dielectric validation, adjoints and gradients                  *!
!* ================================================================================= *!

!> Every PCM flavor must reject a dielectric below one and accept eps = 1.
!> Below eps = 1 the scaling factor f(eps) = (eps-1)/eps turns negative (and
!> diverges at eps = 0), so the model is undefined there; eps = 1 is the
!> vacuum limit, a valid model with f(eps) = 0.
   subroutine test_pcm_invalid_epsilon(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(cpcm) :: pcm_model
      type(cosmo) :: cosmo_model
      !> Dielectric constants that must be rejected
      real(wp), parameter :: bad_epsilon(*) = [0.0_wp, 0.5_wp, -1.0_wp]
      !> Index over the rejected dielectric constants
      integer :: ieps

      !> Run context owned here and borrowed by the components
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      do ieps = 1, size(bad_epsilon)
         call new_cpcm(pcm_model, ctx, bad_epsilon(ieps), error=err)
         call check(error, allocated(err), "CPCM accepted a dielectric below one")
         if (allocated(error)) return
         call check(error, index(err%message, "must be >= 1") > 0, &
            & "CPCM rejected the dielectric with an unexpected message: "//err%message)
         if (allocated(error)) return
         deallocate (err)

         call new_cosmo(cosmo_model, ctx, bad_epsilon(ieps), error=err)
         call check(error, allocated(err), "COSMO accepted a dielectric below one")
         if (allocated(error)) return
         call check(error, index(err%message, "must be >= 1") > 0, &
            & "COSMO rejected the dielectric with an unexpected message: "//err%message)
         if (allocated(error)) return
         deallocate (err)
      end do

      ! The vacuum limit is a valid model, not an error.
      call new_cpcm(pcm_model, ctx, 1.0_wp, error=err)
      call check(error, .not. allocated(err), "CPCM rejected the vacuum limit eps = 1")
      if (allocated(error)) return
      call check(error, pcm_model%feps, 0.0_wp, thr=thr, &
         & message="CPCM f(eps) is not zero at eps = 1")
      if (allocated(error)) return

      call new_cosmo(cosmo_model, ctx, 1.0_wp, error=err)
      call check(error, .not. allocated(err), "COSMO rejected the vacuum limit eps = 1")
      if (allocated(error)) return
      call check(error, cosmo_model%feps, 0.0_wp, thr=thr, &
         & message="COSMO f(eps) is not zero at eps = 1")

   end subroutine test_pcm_invalid_epsilon

!> Component surface adjoints on the shared synthetic seven-point surface.
!>
!> Drives `get_surface_weights` on a hand-built DROP surface and differences it
!> with the shared 4-point harness against an independently coded CPCM energy.
!> The type-bound `amat_surface_weights` / `amat_nuclear_gradient` pair is then
!> cross-checked against the same weights contracted by hand.
   subroutine test_cpcm_surface_weights(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_model
      type(cavity_type_drop) :: cavity
      type(cavity_surface_adjoint_type) :: weights
      type(surface_fixture) :: surface
      type(coupling_type) :: coupling

      !> Dielectric constant
      real(wp), parameter :: epsilon = 78.4_wp
      !> Dielectric scaling shared with the reference energy
      real(wp), parameter :: feps_ref = (epsilon - 1.0_wp)/epsilon

      !> Dummy molecular geometry
      real(wp) :: xyz_mol(3, 1)
      !> Raw A-matrix surface weights from the type-bound interface
      real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :)
      !> Type-bound and independently contracted A-matrix nuclear gradients
      real(wp) :: grad_rA(3, 1), grad_ref(3, 1)
      !> Grid point, nuclear-displacement and Cartesian-component indices
      integer :: igrid, iaxis, icoord
      !> Largest deviation between the two nuclear-gradient contractions
      real(wp) :: gradient_deviation

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      ! A single dummy centre: the external potential source means the molecular
      ! geometry never enters the energy, but update() stores it.
      xyz_mol(:, 1) = 0.0_wp
      call new(mol, [1], xyz_mol)

      ! Synthetic surface carrying only the fields the CPCM matrix reads.
      cavity%ngrid = ngrid_sw
      allocate (cavity%a, source=sw_areas)
      allocate (cavity%xi0, source=sw_xis)
      allocate (cavity%f, source=sw_fs)
      allocate (cavity%xyz, source=sw_xyz)
      allocate (cavity%normal0, source=sw_normals)

      call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%lu, &
         & phi_source=potential_source%external, error=err)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_model%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if

      call pcm_model%input_potential(sw_phi)

      call weights%init(ngrid_sw)
      call pcm_model%get_surface_weights(coupling, cavity, weights, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM surface-weight assembly failed: "//err%message)
         return
      end if

      call new_surface_fixture(surface, sw_areas, sw_xis, sw_fs, sw_xyz, sw_normals)
      call check_surface_weights(error, surface, weights, cpcm_surface_energy, &
         & "cpcm", step=1.0e-4_wp, thr_abs=2.0e-11_wp, thr_rel=2.0e-9_wp, &
         & check_normal=.false.)
      if (allocated(error)) return

      allocate (w_xi(ngrid_sw), w_f(ngrid_sw), w_xyz(3, ngrid_sw))
      call pcm_model%amat_surface_weights(cavity, w_xi, w_f, w_xyz, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM A-matrix surface weights failed: "//err%message)
         return
      end if
      call check(error, maxval(abs(weights%w_xi - 0.5_wp*w_xi/feps_ref)), &
         & 0.0_wp, thr=thr, message="type-bound xi weights disagree with PCM weights")
      if (allocated(error)) return
      call check(error, maxval(abs(weights%w_f - 0.5_wp*w_f/feps_ref)), &
         & 0.0_wp, thr=thr, message="type-bound f weights disagree with PCM weights")
      if (allocated(error)) return
      call check(error, maxval(abs(weights%w_xyz - 0.5_wp*w_xyz/feps_ref)), &
         & 0.0_wp, thr=thr, message="type-bound xyz weights disagree with PCM weights")
      if (allocated(error)) return

      cavity%nsph = 1
      allocate (cavity%xi1_rA(3, 1, ngrid_sw))
      allocate (cavity%f1_rA(3, 1, ngrid_sw))
      allocate (cavity%xyz1_rA(3, 3, 1, ngrid_sw))
      do igrid = 1, ngrid_sw
         do iaxis = 1, 3
            cavity%xi1_rA(iaxis, 1, igrid) = 0.01_wp*real(iaxis + igrid, wp)
            cavity%f1_rA(iaxis, 1, igrid) = -0.02_wp*real(2*iaxis + igrid, wp)
            do icoord = 1, 3
               cavity%xyz1_rA(icoord, iaxis, 1, igrid) = &
                  & 0.005_wp*real(icoord - iaxis + igrid, wp)
            end do
         end do
      end do

      call pcm_model%amat_nuclear_gradient(cavity, grad_rA, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM A-matrix nuclear gradient failed: "//err%message)
         return
      end if

      grad_ref = 0.0_wp
      do igrid = 1, ngrid_sw
         do iaxis = 1, 3
            grad_ref(iaxis, 1) = grad_ref(iaxis, 1) &
               & + w_xi(igrid)*cavity%xi1_rA(iaxis, 1, igrid) &
               & + w_f(igrid)*cavity%f1_rA(iaxis, 1, igrid)
            do icoord = 1, 3
               grad_ref(iaxis, 1) = grad_ref(iaxis, 1) &
                  & + w_xyz(icoord, igrid)*cavity%xyz1_rA(icoord, iaxis, 1, igrid)
            end do
         end do
      end do
      gradient_deviation = maxval(abs(grad_rA - grad_ref))
      call check(error, gradient_deviation, 0.0_wp, thr=thr, &
         & message="type-bound A-matrix nuclear gradient disagrees with surface contraction")

   contains

      !> Rebuild the CPCM energy independently on a perturbed surface.
      function cpcm_surface_energy(trial) result(energy)
         !> Perturbed surface fixture
         type(surface_fixture), intent(in) :: trial
         !> CPCM polarization energy
         real(wp) :: energy

         real(wp) :: amat_trial(ngrid_sw, ngrid_sw)
         real(wp) :: q_trial(ngrid_sw), rhs(ngrid_sw)
         real(wp) :: xi_i, xi_j, xi_ij, r_ij
         integer :: i, j
         type(moist_error_type), allocatable :: solve_error

         amat_trial = 0.0_wp
         do i = 1, ngrid_sw
            amat_trial(i, i) = trial%xi(i)*sqrt(2.0_wp/pi)/trial%f(i)
            do j = 1, i - 1
               xi_i = trial%xi(i)
               xi_j = trial%xi(j)
               xi_ij = xi_i*xi_j/sqrt(xi_i*xi_i + xi_j*xi_j)
               r_ij = norm2(trial%xyz(:, i) - trial%xyz(:, j))
               amat_trial(i, j) = erf(xi_ij*r_ij)/r_ij
               amat_trial(j, i) = amat_trial(i, j)
            end do
         end do

         rhs(:) = -feps_ref*sw_phi(:)
         call solve_pcm_lu(amat_trial, rhs, q_trial, solve_error)
         if (allocated(solve_error)) error stop "cpcm_surface_energy: CPCM solve failed"

         energy = 0.5_wp*dot_product(q_trial, sw_phi)
      end function cpcm_surface_energy

   end subroutine test_cpcm_surface_weights

!> Component surface adjoints against finite differences on a molecular cavity.
!>
!> `test_cpcm_surface_weights` above drives the same code path on a seven-point
!> synthetic surface; this one runs it on an iSwiG cavity built from an mstore
!> structure, where a grid point has hundreds of partners spread across both
!> branches of the pair kernel and the switching factors span several orders of
!> magnitude. The potential is supplied externally and held fixed, so the
!> energy depends on the surface only through the CPCM matrix and the weights
!> `get_surface_weights` accumulates are exactly its gradient:
!>
!>    E = -f(eps)/2 phi^T A^{-1} phi  =>  dE/dp = 1/(2 f(eps)) q^T (dA/dp) q
!>
!> Each channel is differenced with the shared 4-point central formula by
!> perturbing the cavity in place and re-running the component end to end,
!> matrix assembly and solve included.
   subroutine test_cpcm_molecular_surface_weights(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type), allocatable :: mols(:)
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(cavity_surface_adjoint_type) :: weights
      type(coupling_type) :: coupling
      real(wp), allocatable :: phi(:)

      !> Dielectric constant
      real(wp), parameter :: epsilon = 78.4_wp
      !> Lebedev grid size. Every stencil point costs one O(ngrid**3) solve.
      integer, parameter :: nleb = 26
      !> Grid pointe differenced, strided across the exposed surface
      integer, parameter :: n_fd_points = 4
      !> Smallest switching factor a differenced grid point may have.
      !>
      !> An iSwiG surface keeps nearly buried  grid points with f down to ~1e-10, whose
      !> diagonal entry sqrt(2/pi)*xi/f then reaches ~1e10. Their analytic weights
      !> are as valid as any other -- the component treats every grid point alike --
      !> but a relative step on such an f is ~1e-13 wide, so the difference
      !> quotient is pure rounding noise. The sample is therefore drawn from the
      !> exposed part of the surface; the buried diagonal is covered analytically
      !> by the diagonal tests in `amat/kernel.f90`.
      real(wp), parameter :: fd_min_f = 0.1_wp
      !> Smallest exposed surface the sampling accepts, so that a cavity change
      !> which empties the filter fails loudly instead of thinning the test.
      integer, parameter :: fd_min_points = 20
      !> Relative step for the width and switch channels, absolute for positions
      real(wp), parameter :: rel_step = 1.0e-3_wp
      real(wp), parameter :: xyz_step = 1.0e-3_wp
      !> 4-point central FD tolerances
      real(wp), parameter :: fd_atol = 1.0e-12_wp
      real(wp), parameter :: fd_rtol = 1.0e-9_wp

      integer :: ngrid, ip, ig, iax, k, stride, nexposed
      integer, allocatable :: exposed(:)
      real(wp) :: vals(4), step, fd, saved
      character(len=96) :: context

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_test_structures(mols, 5)
      call center_at_origin(mols(1))
      call get_test_cavity_iswig(mols(1), cavity, err, nleb=nleb)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      ngrid = cavity%ngrid
      allocate (phi(ngrid))
      do ig = 1, ngrid
         phi(ig) = 0.05_wp*sin(0.83_wp*real(ig, wp)) - 0.01_wp
      end do

      call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%lu, &
         & phi_source=potential_source%external, error=err)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if

      call weights%init(ngrid)
      call surface_energy(vals(1))
      if (allocated(error)) return
      call pcm_model%get_surface_weights(coupling, cavity, weights, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM surface-weight assembly failed: "//err%message)
         return
      end if

      nexposed = count(cavity%f >= fd_min_f)
      call check(error, nexposed >= fd_min_points, &
         & more="exposed surface too small to difference meaningfully")
      if (allocated(error)) return
      allocate (exposed(nexposed))
      nexposed = 0
      do ig = 1, ngrid
         if (cavity%f(ig) < fd_min_f) cycle
         nexposed = nexposed + 1
         exposed(nexposed) = ig
      end do

      stride = max(1, nexposed/n_fd_points)
      do ip = 1, n_fd_points
         if (1 + (ip - 1)*stride > nexposed) exit
         ig = exposed(1 + (ip - 1)*stride)

         step = rel_step*cavity%xi0(ig)
         saved = cavity%xi0(ig)
         do k = 1, 4
            cavity%xi0(ig) = saved + fd4_offsets(k)*step
            call surface_energy(vals(k))
            if (allocated(error)) return
         end do
         cavity%xi0(ig) = saved
         fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), step)
         write (context, "(a,i0)") "dE/dxi at grid point ", ig
         call check(error, weights%w_xi(ig), fd, thr=fd_atol + fd_rtol*abs(fd), &
            & more=trim(context))
         if (allocated(error)) return

         step = rel_step*cavity%f(ig)
         saved = cavity%f(ig)
         do k = 1, 4
            cavity%f(ig) = saved + fd4_offsets(k)*step
            call surface_energy(vals(k))
            if (allocated(error)) return
         end do
         cavity%f(ig) = saved
         fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), step)
         write (context, "(a,i0)") "dE/df at grid point ", ig
         call check(error, weights%w_f(ig), fd, thr=fd_atol + fd_rtol*abs(fd), &
            & more=trim(context))
         if (allocated(error)) return

         do iax = 1, 3
            saved = cavity%xyz(iax, ig)
            do k = 1, 4
               cavity%xyz(iax, ig) = saved + fd4_offsets(k)*xyz_step
               call surface_energy(vals(k))
               if (allocated(error)) return
            end do
            cavity%xyz(iax, ig) = saved
            fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), xyz_step)
            write (context, "(a,i0,a,i0)") "dE/dxyz axis ", iax, " at grid point ", ig
            call check(error, weights%w_xyz(iax, ig), fd, thr=fd_atol + fd_rtol*abs(fd), &
               & more=trim(context))
            if (allocated(error)) return
         end do
      end do

   contains

      !> Re-run the component on the current cavity and return its energy.
      !> Reports through the enclosing `error`, so callers only have to test it.
      subroutine surface_energy(energy)

         !> CPCM polarization energy on the current cavity
         real(wp), intent(out) :: energy

         type(moist_error_type), allocatable :: local_err

         energy = 0.0_wp
         call pcm_model%update(mols(1), cavity, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "CPCM update failed: "//local_err%message)
            return
         end if
         call pcm_model%input_potential(phi)
         call pcm_model%get_energy(coupling, cavity, energy, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "CPCM energy failed: "//local_err%message)
            return
         end if

      end subroutine surface_energy

   end subroutine test_cpcm_molecular_surface_weights

!> Check the complete CPCM nuclear gradient against finite differences.
!>
!> This exercises the component entry point end to end with the internal
!> point-charge potential. Both terms in
!>
!>    dE/dR_A = q^T dphi/dR_A
!>              + 1/(2 f(eps)) q^T (dA/dR_A) q
!>
!> are active because rebuilding the iSwiG cavity moves its  grid points and changes
!> its switching factors while the source charges move with the nuclei.
!>
!> @param[out] error Test failure
   subroutine test_cpcm_nuclear_gradient(error)

      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Reference and displaced molecular structures
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: trial
      !> PCM component and molecular cavity
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      !> Fixed point-charge coupling data
      type(coupling_type) :: coupling
      !> Deterministic neutral point charges
      real(wp), allocatable :: qat(:)
      !> Analytic PCM nuclear gradient
      real(wp), allocatable :: gradient(:, :)
      !> Four-point finite-difference energies
      real(wp) :: values(4)
      !> Finite-difference derivative and saved coordinate
      real(wp) :: fd, saved
      !> Atom, Cartesian-axis, and stencil indices
      integer :: iatom, iaxis, k
      !> Failure context
      character(len=96) :: context

      !> Dielectric constant
      real(wp), parameter :: epsilon = 78.4_wp
      !> Lebedev grid size
      integer, parameter :: nleb = 26
      !> Nuclear-displacement step in bohr
      real(wp), parameter :: step = 1.0e-4_wp
      !> Tight tolerances accounting for cavity and linear-solve roundoff
      real(wp), parameter :: fd_atol = 5.0e-9_wp
      real(wp), parameter :: fd_rtol = 5.0e-8_wp

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_test_structures(mols, 5)
      call center_at_origin(mols(1))
      call make_neutral_charges(mols(1)%nat, qat)
      coupling%qat = reshape(qat, [size(qat), 1])

      call get_test_cavity_iswig(mols(1), cavity, err, nleb=nleb)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%cholesky, error=err)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if
      call pcm_model%update(mols(1), cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if

      allocate (gradient(3, mols(1)%nat), source=0.0_wp)
      call pcm_model%get_gradient(coupling, cavity, gradient, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM nuclear gradient failed: "//err%message)
         return
      end if
      call check(error, maxval(abs(gradient)) > 1.0e-8_wp, &
         & more="CPCM nuclear gradient is zero, the test is vacuous")
      if (allocated(error)) return

      do iatom = 1, min(2, mols(1)%nat)
         do iaxis = 1, 3
            saved = mols(1)%xyz(iaxis, iatom)
            do k = 1, 4
               trial = mols(1)
               trial%xyz(iaxis, iatom) = saved + fd4_offsets(k)*step
               call displaced_energy(trial, values(k))
               if (allocated(error)) return
            end do
            fd = fd4_scalar(values(1), values(2), values(3), values(4), step)
            write (context, "(a,i0,a,i0)") "CPCM gradient atom ", iatom, &
               & ", axis ", iaxis
            call check(error, gradient(iaxis, iatom), fd, &
               & thr_abs=fd_atol, thr_rel=fd_rtol, more=trim(context))
            if (allocated(error)) return
         end do
      end do

   contains

      !> Evaluate the CPCM energy for one displaced molecular structure.
      !> @param[in]  displaced_mol Displaced molecular structure
      !> @param[out] energy        CPCM polarization energy
      subroutine displaced_energy(displaced_mol, energy)
         !> Displaced molecular structure
         type(structure_type), intent(in) :: displaced_mol
         !> CPCM polarization energy
         real(wp), intent(out) :: energy

         !> Local library error handling
         type(moist_error_type), allocatable :: local_err

         call get_test_cavity_iswig(displaced_mol, cavity, local_err, nleb=nleb)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced cavity setup failed: "//local_err%message)
            return
         end if
         call pcm_model%update(displaced_mol, cavity, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced CPCM update failed: "//local_err%message)
            return
         end if
         energy = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced CPCM energy failed: "//local_err%message)
         end if
      end subroutine displaced_energy

   end subroutine test_cpcm_nuclear_gradient

!> Check the external-potential CPCM nuclear gradient against finite differences.
!>
!> The supplied potential is the sum of the solute nuclear potential and a
!> linear electronic potential. The latter has a constant spatial derivative,
!> so `coupling%elstat_qefield(:,i) = q_i*dphi_elec/dr_i` is known exactly. This
!> exercises the external branch's nuclear-charge construction and electronic
!> field contribution through the complete component entry point.
!>
!> @param[out] error Test failure
   subroutine test_cpcm_external_nuclear_gradient(error)

      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Reference and displaced molecular structures
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: trial
      !> PCM component and molecular cavity
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      !> Valid and deliberately malformed external coupling data
      type(coupling_type) :: coupling, bad_coupling
      !> External potential on the current cavity
      real(wp), allocatable :: phi(:)
      !> Analytic gradient and rejected-call sentinel
      real(wp), allocatable :: gradient(:, :), bad_gradient(:, :)
      !> Four-point finite-difference energies
      real(wp) :: values(4)
      !> Finite-difference derivative and saved coordinate
      real(wp) :: fd, saved
      !> Surface, Cartesian-axis, and stencil indices
      integer :: igrid, iaxis, k
      !> Failure context
      character(len=96) :: context

      !> Dielectric constant
      real(wp), parameter :: epsilon = 78.4_wp
      !> Lebedev grid size
      integer, parameter :: nleb = 26
      !> Constant gradient of the synthetic electronic potential
      real(wp), parameter :: electronic_field(3) = [0.013_wp, -0.017_wp, 0.009_wp]
      !> Nuclear-displacement step in bohr
      real(wp), parameter :: step = 1.0e-4_wp
      !> Tight tolerances accounting for cavity and linear-solve roundoff
      real(wp), parameter :: fd_atol = 5.0e-9_wp
      real(wp), parameter :: fd_rtol = 5.0e-8_wp

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_test_structures(mols, 5)
      call center_at_origin(mols(1))
      call get_test_cavity_iswig(mols(1), cavity, err, nleb=nleb)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%cholesky, &
         & phi_source=potential_source%external, error=err)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if
      call pcm_model%update(mols(1), cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if

      call build_external_potential(mols(1), cavity, phi)
      call pcm_model%input_potential(phi)
      call pcm_model%ensure_charges(coupling, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM charge solution failed: "//err%message)
         return
      end if

      allocate (coupling%elstat_qefield(3, cavity%ngrid))
      do igrid = 1, cavity%ngrid
         coupling%elstat_qefield(:, igrid) = pcm_model%q(igrid)*electronic_field
      end do
      allocate (gradient(3, mols(1)%nat), source=0.0_wp)
      call pcm_model%get_gradient(coupling, cavity, gradient, err)
      if (allocated(err)) then
         call test_failed(error, "External CPCM nuclear gradient failed: "//err%message)
         return
      end if
      call check(error, maxval(abs(coupling%elstat_qefield)) > 1.0e-8_wp, &
         & more="external electronic-field contribution is zero, the test is vacuous")
      if (allocated(error)) return

      bad_coupling = coupling
      deallocate (bad_coupling%elstat_qefield)
      allocate (bad_coupling%elstat_qefield(2, cavity%ngrid), source=0.0_wp)
      allocate (bad_gradient(3, mols(1)%nat), source=1.0_wp)
      call pcm_model%get_gradient(bad_coupling, cavity, bad_gradient, err)
      call check(error, allocated(err), &
         & more="external CPCM gradient accepted a malformed electronic field")
      if (allocated(error)) return
      call check(error, maxval(abs(bad_gradient - 1.0_wp)), 0.0_wp, thr=0.0_wp, &
         & more="rejected external CPCM gradient modified the accumulator")
      if (allocated(error)) return
      deallocate (err)

      do iaxis = 1, 3
         saved = mols(1)%xyz(iaxis, 1)
         do k = 1, 4
            trial = mols(1)
            trial%xyz(iaxis, 1) = saved + fd4_offsets(k)*step
            call displaced_external_energy(trial, values(k))
            if (allocated(error)) return
         end do
         fd = fd4_scalar(values(1), values(2), values(3), values(4), step)
         write (context, "(a,i0)") "external CPCM gradient atom 1, axis ", iaxis
         call check(error, gradient(iaxis, 1), fd, &
            & thr_abs=fd_atol, thr_rel=fd_rtol, more=trim(context))
         if (allocated(error)) return
      end do

   contains

      !> Build the synthetic nuclear-plus-electronic external potential.
      !> @param[in]  potential_mol    Molecular structure providing nuclear sources
      !> @param[in]  potential_cavity Cavity carrying evaluation points
      !> @param[out] potential        External molecular potential
      subroutine build_external_potential(potential_mol, potential_cavity, potential)
         !> Molecular structure providing nuclear sources
         type(structure_type), intent(in) :: potential_mol
         !> Cavity carrying evaluation points
         type(cavity_type_iswig), intent(in) :: potential_cavity
         !> External molecular potential
         real(wp), allocatable, intent(out) :: potential(:)

         !> Surface and atom indices
         integer :: i, iatom
         !> Nuclear charge
         real(wp) :: za

         allocate (potential(potential_cavity%ngrid), source=0.0_wp)
         do i = 1, potential_cavity%ngrid
            do iatom = 1, potential_mol%nat
               za = real(potential_mol%num(potential_mol%id(iatom)), wp)
               potential(i) = potential(i) + za/norm2( &
                  & potential_cavity%xyz(:, i) - potential_mol%xyz(:, iatom))
            end do
            potential(i) = potential(i) + &
               & dot_product(electronic_field, potential_cavity%xyz(:, i))
         end do
      end subroutine build_external_potential

      !> Evaluate the external-potential CPCM energy for one displaced structure.
      !> @param[in]  displaced_mol Displaced molecular structure
      !> @param[out] energy        CPCM polarization energy
      subroutine displaced_external_energy(displaced_mol, energy)
         !> Displaced molecular structure
         type(structure_type), intent(in) :: displaced_mol
         !> CPCM polarization energy
         real(wp), intent(out) :: energy

         !> Local library error handling
         type(moist_error_type), allocatable :: local_err
         !> External potential on the displaced cavity
         real(wp), allocatable :: displaced_phi(:)

         call get_test_cavity_iswig(displaced_mol, cavity, local_err, nleb=nleb)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced cavity setup failed: "//local_err%message)
            return
         end if
         call pcm_model%update(displaced_mol, cavity, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced CPCM update failed: "//local_err%message)
            return
         end if
         call build_external_potential(displaced_mol, cavity, displaced_phi)
         call pcm_model%input_potential(displaced_phi)
         energy = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced CPCM energy failed: "//local_err%message)
         end if
      end subroutine displaced_external_energy

   end subroutine test_cpcm_external_nuclear_gradient

!> The CPCM energy is invariant under rigid translation of the solute.
!>
!> Every iSwiG grid point is placed relative to an atomic centre, so shifting the
!> whole structure shifts the surface with it and leaves every distance in both
!> the matrix and the potential trace unchanged. Unlike rotation -- which
!> reorients the molecule against a lab-fixed Lebedev grid and therefore only
!> holds up to discretization error -- this is an exact identity, and it fails
!> the moment any part of the assembly picks up an absolute coordinate.
   subroutine test_cpcm_translation_invariance(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: shifted
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      real(wp), allocatable :: qat(:)

      !> Dielectric constant
      real(wp), parameter :: epsilon = 78.4_wp
      !> Lebedev grid size
      integer, parameter :: nleb = 26
      !> Rigid shift applied to the whole structure, bohr
      real(wp), parameter :: shift(3) = [1.7_wp, -2.3_wp, 0.9_wp]

      integer :: imol, iat
      real(wp) :: energy, energy_shifted
      character(len=64) :: context

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_test_structures(mols, 5)

      do imol = 1, size(mols)
         call center_at_origin(mols(imol))
         call make_neutral_charges(mols(imol)%nat, qat)
         coupling%qat = reshape(qat, [size(qat), 1])

         call get_test_cavity_iswig(mols(imol), cavity, err, nleb=nleb)
         if (allocated(err)) then
            call test_failed(error, "Cavity setup failed: "//err%message)
            return
         end if
         call new_cpcm(pcm_model, ctx, epsilon, solver=solver_type%cholesky, error=err)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed: "//err%message)
            return
         end if
         call pcm_model%update(mols(imol), cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed: "//err%message)
            return
         end if
         energy = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed: "//err%message)
            return
         end if

         shifted = mols(imol)
         do iat = 1, shifted%nat
            shifted%xyz(:, iat) = shifted%xyz(:, iat) + shift
         end do

         call get_test_cavity_iswig(shifted, cavity, err, nleb=nleb)
         if (allocated(err)) then
            call test_failed(error, "Shifted cavity setup failed: "//err%message)
            return
         end if
         call pcm_model%update(shifted, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "Shifted CPCM update failed: "//err%message)
            return
         end if
         energy_shifted = 0.0_wp
         call pcm_model%get_energy(coupling, cavity, energy_shifted, err)
         if (allocated(err)) then
            call test_failed(error, "Shifted CPCM energy failed: "//err%message)
            return
         end if

         ! Guard against a vacuous comparison: a structure whose polarization
         ! energy happened to vanish would pass any invariance test.
         write (context, "(a,i0)") "structure ", imol
         call check(error, abs(energy) > 1.0e-6_wp, &
            & more=trim(context)//": CPCM energy is zero, invariance is vacuous")
         if (allocated(error)) return
         call check(error, energy_shifted, energy, thr=1.0e-13_wp*abs(energy), &
            & more=trim(context)//": CPCM energy is not translation invariant")
         if (allocated(error)) return

         deallocate (qat)
      end do

   end subroutine test_cpcm_translation_invariance

!> The CPCM energy is exactly linear in the dielectric factor f(eps).
!>
!> CPCM solves A q = -f(eps) phi and reports E = 1/2 q.phi, so at fixed geometry
!> and fixed potential the whole dielectric dependence is the single scalar
!> f(eps) = (eps-1)/eps. Energies computed at different dielectrics must
!> therefore rescale exactly, with no discretization error of their own: the
!> matrix, its factorization and the trace are identical in all of them. This
!> pins down the one place the dielectric enters and would catch, for instance,
!> f(eps) being applied twice or to the wrong side of the solve.
   subroutine test_cpcm_dielectric_scaling(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type), allocatable :: mols(:)
      type(cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      real(wp), allocatable :: qat(:)

      !> Dielectrics spanning the weakly polar to the conductor limit
      real(wp), parameter :: epsilons(4) = [2.0_wp, 8.93_wp, 78.4_wp, 1.0e6_wp]
      !> Lebedev grid size
      integer, parameter :: nleb = 26

      integer :: imol, ieps
      real(wp) :: energy, e_ref, feps, feps_ref, predicted
      character(len=96) :: context

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_test_structures(mols, 5)

      do imol = 1, size(mols)
         call center_at_origin(mols(imol))
         call make_neutral_charges(mols(imol)%nat, qat)
         coupling%qat = reshape(qat, [size(qat), 1])

         call get_test_cavity_iswig(mols(imol), cavity, err, nleb=nleb)
         if (allocated(err)) then
            call test_failed(error, "Cavity setup failed: "//err%message)
            return
         end if

         e_ref = 0.0_wp
         feps_ref = 0.0_wp
         do ieps = 1, size(epsilons)
            feps = (epsilons(ieps) - 1.0_wp)/epsilons(ieps)
            call new_cpcm(pcm_model, ctx, epsilons(ieps), solver=solver_type%cholesky, error=err)
            if (allocated(err)) then
               call test_failed(error, "CPCM initialization failed: "//err%message)
               return
            end if
            call pcm_model%update(mols(imol), cavity, err)
            if (allocated(err)) then
               call test_failed(error, "CPCM update failed: "//err%message)
               return
            end if
            energy = 0.0_wp
            call pcm_model%get_energy(coupling, cavity, energy, err)
            if (allocated(err)) then
               call test_failed(error, "CPCM energy failed: "//err%message)
               return
            end if

            if (ieps == 1) then
               e_ref = energy
               feps_ref = feps
               write (context, "(a,i0)") "structure ", imol
               call check(error, abs(e_ref) > 1.0e-6_wp, &
                  & more=trim(context)//": CPCM energy is zero, scaling is vacuous")
               if (allocated(error)) return
               cycle
            end if

            predicted = e_ref*feps/feps_ref
            write (context, "(a,i0,a,f0.2)") "structure ", imol, ", eps = ", epsilons(ieps)
            call check(error, energy, predicted, thr=1.0e-14_wp*abs(predicted), &
               & more=trim(context)//": CPCM energy does not scale with f(eps)")
            if (allocated(error)) return
         end do

         deallocate (qat)
      end do

   end subroutine test_cpcm_dielectric_scaling

!> Bit-coincident  grid points must not poison the CPCM matrix.
!>
!> In symmetric molecules two Lebedev points owned by different spheres can land
!> at bit-identical coordinates (the Td axis of NH4+/PH4+, a linear anion at
!> certain radii scales). A literal erf(xi_ij r_ij)/r_ij kernel evaluates 0/0
!> there and the NaN propagates into every downstream solve. moist evaluates the
!> off-diagonal through the Boys function instead, A_ij = 2/sqrt(pi) xi_ij F0(T)
!> with T = (xi_ij r_ij)^2, which is regular at T = 0 and returns exactly the
!> zero-separation limit 2 xi_ij/sqrt(pi). This pins that property.
   subroutine test_cpcm_coincident_points(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(cavity_type_drop) :: cavity

      !> Three  grid points, the first two bit-identical with distinct widths
      integer, parameter :: ngrid_c = 3
      real(wp), parameter :: xis(ngrid_c) = [0.65_wp, 0.83_wp, 0.71_wp]
      real(wp), parameter :: fs(ngrid_c) = [0.91_wp, 0.88_wp, 0.79_wp]
      real(wp), parameter :: xyz(3, ngrid_c) = reshape([ &
                                               1.5_wp, 0.0_wp, 0.0_wp, &
                                               1.5_wp, 0.0_wp, 0.0_wp, &
                                               -1.1_wp, 0.7_wp, 0.3_wp], [3, ngrid_c])

      real(wp), allocatable :: amat(:, :)
      real(wp) :: xi_ij, expected

      cavity%ngrid = ngrid_c
      allocate (cavity%xi0, source=xis)
      allocate (cavity%f, source=fs)
      allocate (cavity%xyz, source=xyz)

      allocate (amat(ngrid_c, ngrid_c))
      call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, amat, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM matrix assembly failed: "//err%message)
         return
      end if

      call check(error, all(ieee_is_finite_matrix(amat)), &
         & "coincident  grid points must not produce non-finite CPCM matrix entries")
      if (allocated(error)) return

      ! The coincident pair must carry the analytic zero-separation limit.
      xi_ij = xis(1)*xis(2)/sqrt(xis(1)**2 + xis(2)**2)
      expected = 2.0_wp*xi_ij/sqrt(pi)
      call check(error, amat(1, 2), expected, thr=thr, &
         & message="coincident pair must equal the analytic r -> 0 limit")
      if (allocated(error)) return

      call check(error, amat(2, 1), expected, thr=thr, &
         & message="the CPCM matrix must stay symmetric at coincident  grid points")

   contains

      !> Elementwise finite test for a rank-2 array.
      pure function ieee_is_finite_matrix(a) result(finite)
         !> Matrix to test
         real(wp), intent(in) :: a(:, :)
         !> Elementwise finiteness
         logical :: finite(size(a, 1), size(a, 2))

         finite = (a == a) .and. (abs(a) <= huge(1.0_wp))
      end function ieee_is_finite_matrix

   end subroutine test_cpcm_coincident_points

!> Deterministic, charge-neutral atomic charges for a structure.
!> @param[in]  nat  Number of atoms
!> @param[out] qat  Atomic charges summing to zero
   subroutine make_neutral_charges(nat, qat)

      !> Number of atoms
      integer, intent(in) :: nat

      !> Atomic charges summing to zero
      real(wp), allocatable, intent(out) :: qat(:)

      integer :: iat

      allocate (qat(nat))
      do iat = 1, nat
         qat(iat) = 0.2_wp*sin(1.3_wp*real(iat, wp))
      end do
      qat = qat - sum(qat)/real(nat, wp)

   end subroutine make_neutral_charges

end module test_model_component_pcm_cpcm
