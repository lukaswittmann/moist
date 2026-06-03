!> Unit tests for the CPCM solvation model component.
!! Tests energy evaluation, solver consistency, and PCM API behavior.
module test_component_pcm_cpcm
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use moist_type, only: wavefunction_type, potential_type
   use moist_model_component_pcm_type, only: solver_type, potential_source
   use moist_model_component_pcm_cpcm, only: cpcm, new_cpcm
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_radii, only: static_radius_type, new_cosmo_radii
   implicit none
   private

   public :: collect_component_pcm_cpcm

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

!> Convert integer to string
   pure function to_string(i) result(str)
      integer, intent(in) :: i
      character(len=20) :: str
      write (str, '(i0)') i
   end function to_string

!> Collect all exported unit tests
   subroutine collect_component_pcm_cpcm(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("CPCM LU Energy", test_cpcm_energy_lu), &
         & new_unittest("CPCM Charged System", test_cpcm_charged), &
         & new_unittest("CPCM Solver Comparison", test_cpcm_solver_comparison), &
         & new_unittest("CPCM Vacuum Limit", test_cpcm_vacuum_limit), &
         & new_unittest("CPCM External Potential", test_cpcm_external_potential), &
         & new_unittest("CPCM External Potential Requires Input", test_cpcm_external_potential_requires_input, &
            should_fail=.true.), &
         & new_unittest("CPCM External Matrix", test_cpcm_external_matrix), &
         & new_unittest("CPCM Spin-Resolved Charges", test_cpcm_spin_resolved_charges), &
         & new_unittest("CPCM Requires Update", test_cpcm_requires_update, should_fail=.true.), &
         & new_unittest("CPCM Invalid Solver", test_cpcm_invalid_solver, should_fail=.true.), &
         & new_unittest("CPCM Reallocates On Grid Change", test_cpcm_reallocate_on_ngrid_change) &
         ! & new_unittest("CPCM Solver Timing", test_cpcm_timing) &
         & ]

   end subroutine collect_component_pcm_cpcm

!> Build a CPCM test cavity for a given molecule and Lebedev grid size.
   subroutine build_test_cavity(mol, nleb, radius_model, cavity, error)

      !> Molecular structure
      type(structure_type), intent(in) :: mol

      !> Lebedev grid size
      integer, intent(in) :: nleb

      !> Radius model storage
      type(static_radius_type), intent(out) :: radius_model

      !> Constructed cavity
      type(cavity_type_iswig), intent(out) :: cavity

      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      call new_cosmo_radii(radius_model)
      call new_cavity_iswig(cavity, nleb=nleb, radius_model=radius_model, error=error)
      if (allocated(error)) return

      call cavity%update(mol, error=error)

   end subroutine build_test_cavity

!> Build a single-column wavefunction charge array.
   subroutine make_charge_wfn(qat, wfn)

      !> Atomic charges
      real(wp), intent(in) :: qat(:)

      !> Wavefunction to populate
      type(wavefunction_type), intent(out) :: wfn

      wfn%qat = reshape(qat, [size(qat), 1])

   end subroutine make_charge_wfn

!> Test CPCM energy calculation with all solvers on neutral system
   subroutine test_cpcm_energy_lu(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), parameter :: qat(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp), parameter :: epsilon = 78.4_wp
      real(wp), parameter :: ref_energy = -5.1913485531103667E-3_wp

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
      real(wp), parameter :: ref_energy = -0.12166653890956987_wp

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
      type(wavefunction_type) :: wfn
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp), parameter :: epsilon = 78.4_wp
      integer :: solvers(4)
      character(len=20) :: solver_names(4)
      real(wp) :: energies(4)
      real(wp), allocatable :: charges(:, :)
      real(wp) :: energy_array(1)
      real(wp) :: energy_diff, charge_rms
      integer :: i, j

      ! Setup solver types and names
      solvers = [solver_type%inversion, solver_type%lu, &
                 solver_type%cholesky, solver_type%iterative]
      solver_names = ["inversion ", "lu        ", &
                      "cholesky  ", "iterative "]

      ! Get test molecule
      call get_structure(mol, "MB16-43", "01")

      ! Build cavity
      call new_cosmo_radii(radius_model)
      call new_cavity_iswig(cavity, nleb=50, radius_model=radius_model, error=err)
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
      wfn%qat = reshape(qat_vals, [size(qat_vals), 1])

      ! Allocate storage for all solver results
      allocate (charges(cavity%ngrid, 4))

      ! Test all 4 solvers
      do i = 1, 4
         call new_cpcm(pcm_model, epsilon, solver=solvers(i), error=err)
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
         call pcm_model%get_energy(wfn, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed ("//trim(solver_names(i))//"): "//err%message)
            return
         end if
         energies(i) = energy_array(1)
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
      type(wavefunction_type) :: wfn
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      integer :: solvers(4)
      character(len=20) :: solver_names(4)
      real(wp) :: energy_array(1)
      integer :: i

      solvers = [solver_type%inversion, solver_type%lu, &
                 solver_type%cholesky, solver_type%iterative]
      solver_names = ["inversion ", "lu        ", &
                      "cholesky  ", "iterative "]

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_wfn(qat_vals, wfn)

      do i = 1, 4
         call new_cpcm(pcm_model, 1.0_wp, solver=solvers(i), error=err)
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
         call pcm_model%get_energy(wfn, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed ("//trim(solver_names(i))//"): "//err%message)
            return
         end if

         call check(error, energy_array(1), 0.0_wp, thr=thr, &
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
      type(wavefunction_type) :: wfn
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_internal(1), energy_external(1), charge_rms
      real(wp), allocatable :: phi_ref(:), q_ref(:)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_wfn(qat_vals, wfn)

      call new_cpcm(pcm_internal, 78.4_wp, solver=solver_type%cholesky, error=err)
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
      call pcm_internal%get_energy(wfn, energy_internal, err)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM energy failed: "//err%message)
         return
      end if

      allocate (phi_ref, source=pcm_internal%phi)
      allocate (q_ref, source=pcm_internal%q)

      call new_cpcm(pcm_external, 78.4_wp, solver=solver_type%cholesky, &
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
      call pcm_external%get_energy(wfn, energy_external, err)
      if (allocated(err)) then
         call test_failed(error, "External-potential CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy_external(1), energy_internal(1), thr=thr, &
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
      type(wavefunction_type) :: wfn
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_array(1)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_wfn(qat_vals, wfn)

      call new_cpcm(pcm_model, 78.4_wp, solver=solver_type%cholesky, &
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
      call pcm_model%get_energy(wfn, energy_array, err)
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
      type(wavefunction_type) :: wfn
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_internal(1), energy_external(1), charge_rms
      real(wp), allocatable :: amat_ref(:, :), q_ref(:)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_wfn(qat_vals, wfn)

      call new_cpcm(pcm_internal, 78.4_wp, solver=solver_type%lu, error=err)
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
      call pcm_internal%get_energy(wfn, energy_internal, err)
      if (allocated(err)) then
         call test_failed(error, "Internal CPCM energy failed: "//err%message)
         return
      end if

      allocate (amat_ref, source=pcm_internal%amat)
      allocate (q_ref, source=pcm_internal%q)

      call new_cpcm(pcm_external, 78.4_wp, solver=solver_type%lu, &
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
      call pcm_external%get_energy(wfn, energy_external, err)
      if (allocated(err)) then
         call test_failed(error, "External-matrix CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy_external(1), energy_internal(1), thr=thr, &
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
      type(wavefunction_type) :: wfn_scalar, wfn_spin
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_scalar(1), energy_spin(1), charge_rms

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if

      call make_charge_wfn(qat_vals, wfn_scalar)
      allocate (wfn_spin%qat(size(qat_vals), 2))
      wfn_spin%qat(:, 1) = 0.25_wp*qat_vals
      wfn_spin%qat(:, 2) = 0.75_wp*qat_vals

      call new_cpcm(pcm_scalar, 78.4_wp, solver=solver_type%cholesky, error=err)
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
      call pcm_scalar%get_energy(wfn_scalar, energy_scalar, err)
      if (allocated(err)) then
         call test_failed(error, "Scalar-charge CPCM energy failed: "//err%message)
         return
      end if

      call new_cpcm(pcm_spin, 78.4_wp, solver=solver_type%cholesky, error=err)
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
      call pcm_spin%get_energy(wfn_spin, energy_spin, err)
      if (allocated(err)) then
         call test_failed(error, "Spin-charge CPCM energy failed: "//err%message)
         return
      end if

      call check(error, energy_spin(1), energy_scalar(1), thr=thr, &
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
      type(wavefunction_type) :: wfn
      real(wp), parameter :: qat_vals(1) = [0.0_wp]
      real(wp) :: energy_array(1)

      call make_charge_wfn(qat_vals, wfn)
      call new_cpcm(pcm_model, 78.4_wp, error=err)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if

      energy_array = 0.0_wp
      call pcm_model%get_energy(wfn, energy_array, err)
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
      type(wavefunction_type) :: wfn
      type(static_radius_type) :: radius_model
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_array(1)

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 50, radius_model, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Cavity setup failed: "//err%message)
         return
      end if
      call make_charge_wfn(qat_vals, wfn)

      call new_cpcm(pcm_model, 78.4_wp, solver=-1, error=err)
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
      call pcm_model%get_energy(wfn, energy_array, err)
      if (.not. allocated(err)) return

      if (index(err%message, "Unknown solver type") == 0) return

      call test_failed(error, "Invalid solver correctly rejected: "//trim(err%message))

   end subroutine test_cpcm_invalid_solver

!> Test that a PCM instance can be reused safely when the cavity grid size changes.
   subroutine test_cpcm_reallocate_on_ngrid_change(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(cpcm) :: pcm_reused, pcm_fresh
      type(cavity_type_iswig) :: cavity_small, cavity_large
      type(wavefunction_type) :: wfn
      type(static_radius_type) :: radius_small, radius_large
      real(wp), parameter :: qat_vals(*) = [&
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, &
         &  0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp, 0.1_wp, -0.1_wp]
      real(wp) :: energy_reused(1), energy_fresh(1), charge_rms

      call get_structure(mol, "MB16-43", "01")
      call build_test_cavity(mol, 14, radius_small, cavity_small, err)
      if (allocated(err)) then
         call test_failed(error, "Small cavity setup failed: "//err%message)
         return
      end if
      call build_test_cavity(mol, 50, radius_large, cavity_large, err)
      if (allocated(err)) then
         call test_failed(error, "Large cavity setup failed: "//err%message)
         return
      end if
      call make_charge_wfn(qat_vals, wfn)

      call new_cpcm(pcm_reused, 78.4_wp, solver=solver_type%lu, error=err)
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
      call pcm_reused%get_energy(wfn, energy_reused, err)
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
      call pcm_reused%get_energy(wfn, energy_reused, err)
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

      call new_cpcm(pcm_fresh, 78.4_wp, solver=solver_type%lu, error=err)
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
      call pcm_fresh%get_energy(wfn, energy_fresh, err)
      if (allocated(err)) then
         call test_failed(error, "Fresh energy with large cavity failed: "//err%message)
         return
      end if

      call check(error, energy_reused(1), energy_fresh(1), thr=thr, &
         & message="Reused CPCM model changed energy after grid-size change")
      if (allocated(error)) return

      charge_rms = sqrt(sum((pcm_reused%q - pcm_fresh%q)**2)/real(size(pcm_fresh%q), wp))
      call check(error, charge_rms, 0.0_wp, thr=thr2, &
         & message="Reused CPCM model changed surface charges after grid-size change")
      if (allocated(error)) return

   end subroutine test_cpcm_reallocate_on_ngrid_change

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
      type(wavefunction_type) :: wfn
      real(wp) :: energy_array(1)
      real(wp) :: energy_lu, energy_cholesky, energy_iterative, energy_inversion
      type(static_radius_type) :: radius_model
      integer :: solvers(4)
      character(len=20) :: solver_names(4)
      real(wp) :: energies(4)
      real(wp) :: energy_diff
      integer :: i

      ! Setup solver types and names
      solvers = [solver_type%inversion, solver_type%lu, &
                 solver_type%cholesky, solver_type%iterative]
      solver_names = ["inversion ", "lu        ", &
                      "cholesky  ", "iterative "]

      ! Build cavity and wavefunction
      call new_cosmo_radii(radius_model)
      call new_cavity_iswig(cavity, nleb=50, radius_model=radius_model, error=err)
      if (allocated(err)) then
         call test_failed(error, "Cavity initialization failed: "//err%message)
         return
      end if
      call cavity%update(mol, error=err)
      if (allocated(err)) then
         call test_failed(error, "Cavity update failed: "//err%message)
         return
      end if
      wfn%qat = reshape(qat, [size(qat), 1])

      ! Test all 4 solvers
      do i = 1, 4
         call new_cpcm(pcm_model, epsilon, solver=solvers(i), error=err)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed for "// &
                             trim(solver_names(i))//" solver")
            return
         end if

         call pcm_model%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed for "//trim(solver_names(i))//": "//err%message)
            return
         end if
         energy_array = 0.0_wp
         call pcm_model%get_energy(wfn, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed for "//trim(solver_names(i))//": "//err%message)
            return
         end if
         energies(i) = energy_array(1)

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
      type(wavefunction_type) :: wfn
      real(wp) :: energy_array(1)
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
         wfn%qat = reshape(qat, [size(qat), 1])

         ! Build cavity
         call new_cosmo_radii(radius_model)
         call new_cavity_iswig(cavity, nleb=50, radius_model=radius_model, error=err)
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
         call new_cpcm(pcm_model, epsilon, solver=solver_type%lu, error=err)
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
         call pcm_model%get_energy(wfn, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed (lu): "//err%message)
            return
         end if
         energy_lu = energy_array(1)
         allocate (charges_lu(cavity%ngrid))
         charges_lu = pcm_model%q
         call system_clock(t2)
         time_lu = real(t2 - t1, wp)/real(rate, wp)

         ! ===== Time Cholesky solver =====
         call system_clock(t1)
         call new_cpcm(pcm_model, epsilon, solver=solver_type%cholesky, error=err)
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
         call pcm_model%get_energy(wfn, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed (cholesky): "//err%message)
            return
         end if
         energy_cholesky = energy_array(1)
         allocate (charges_cholesky(cavity%ngrid))
         charges_cholesky = pcm_model%q
         call system_clock(t2)
         time_cholesky = real(t2 - t1, wp)/real(rate, wp)

         ! ===== Time iterative solver =====
         call system_clock(t1)
         call new_cpcm(pcm_model, epsilon, solver=solver_type%iterative, error=err)
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
         call pcm_model%get_energy(wfn, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed (iterative): "//err%message)
            return
         end if
         energy_iterative = energy_array(1)
         allocate (charges_iterative(cavity%ngrid))
         charges_iterative = pcm_model%q
         call system_clock(t2)
         time_iterative = real(t2 - t1, wp)/real(rate, wp)

         ! ===== Time inversion solver =====
         call system_clock(t1)
         call new_cpcm(pcm_model, epsilon, solver=solver_type%inversion, error=err)
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
         call pcm_model%get_energy(wfn, energy_array, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed (inversion): "//err%message)
            return
         end if
         energy_inversion = energy_array(1)
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

end module test_component_pcm_cpcm
