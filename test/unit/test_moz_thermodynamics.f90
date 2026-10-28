!> Tests of MOZ thermodynamic kernels, units and pressure-correction responses
module test_moz_thermodynamics
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type, fatal_error
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_model_moz_thermodynamics_kernel, only: bulk_compressibility_denominator, bulk_compressibility, &
      & bulk_pressure, bulk_ideal_pressure, partial_molar_volume, pc_correction, pc_plus_correction, &
      & pc_moment_shift, pc_plus_moment_shift, coulomb_box_mean, electrostatic_finite_size_correction
   use moist_model_moz, only: model_moz_3d_thermodynamics_type
   use moist_model_moz_thermodynamics, only: new_moz_3d_thermodynamics, &
      & correction_moment_shift, thermodynamic_quantity_type, thermo_list_type, new_uv_thermo_list, new_bulk_thermo_list, &
      & accumulate_solvation_adjoint, huenenberger_correction_type, new_huenenberger_correction
   use moist_model_moz_solvent_vv_thermodynamics, only: solvent_vv_thermodynamics_type
   use test_helpers, only: check_moist_error
   implicit none(type, external)
   private
   public :: collect_moz_thermodynamics

   !> Test-only linear quantity/correction, with an optional guard against reevaluation
   type, extends(thermodynamic_quantity_type) :: linear_item_type
      real(wp), allocatable :: coefficients(:)
      integer :: evaluations = 0
      logical :: once = .false.
   contains
      procedure :: evaluate => linear_evaluate
      procedure :: accumulate_adjoint => linear_accumulate
   end type linear_item_type

contains

   subroutine collect_moz_thermodynamics(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("reference_values", check_reference_values), &
         & new_unittest("corrected_energy", check_corrected_energy), &
         & new_unittest("finite_size", check_finite_size), &
         & new_unittest("finite_size_response", check_finite_size_response), &
         & new_unittest("beta_scaling", check_beta_scaling), &
         & new_unittest("moment_derivatives", check_moment_derivatives), &
         & new_unittest("shared_dependencies", check_shared_dependencies), &
         & new_unittest("selection", check_selection), &
         & new_unittest("graph_derivatives", check_graph_derivatives), &
         & new_unittest("list_errors", check_list_errors), &
         & new_unittest("parallel_lists", check_parallel_lists)]
   end subroutine collect_moz_thermodynamics

   !> Independent cubic Coulomb integral, charge parity and rectangular-box symmetries
   subroutine check_finite_size(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(huenenberger_correction_type) :: quantity
      real(wp) :: mean, energy, empty(0)

      mean = coulomb_box_mean([40.0_wp, 40.0_wp, 40.0_wp])
      call check(error, mean, 2.380077363979553_wp/40.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, coulomb_box_mean([20.0_wp, 30.0_wp, 50.0_wp]), &
         & coulomb_box_mean([50.0_wp, 20.0_wp, 30.0_wp]), thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, coulomb_box_mean([40.0_wp, 60.0_wp, 100.0_wp]), &
         & 0.5_wp*coulomb_box_mean([20.0_wp, 30.0_wp, 50.0_wp]), thr=1.0e-14_wp)
      if (allocated(error)) return
      call new_huenenberger_correction(quantity, -1.0_wp, 0.95_wp, [40.0_wp, 40.0_wp, 40.0_wp], err)
      if (.not. allocated(err)) call quantity%evaluate(empty, energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, energy, -0.5_wp*0.95_wp*mean, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, energy, electrostatic_finite_size_correction(1.0_wp, 0.95_wp, mean), thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, electrostatic_finite_size_correction(0.0_wp, 0.95_wp, mean), 0.0_wp)
      if (allocated(error)) return
      call check(error, electrostatic_finite_size_correction(1.0_wp, 0.0_wp, mean), 0.0_wp)
      if (allocated(error)) return
      call new_huenenberger_correction(quantity, 1.0_wp, 0.95_wp, [40.0_wp, 0.0_wp, 40.0_wp], err)
      call check(error, allocated(err), more="a finite-size correction without a periodic box was accepted")
   end subroutine check_finite_size

   !> A fixed charge correction changes energy but neither UV adjoint seed
   subroutine check_finite_size_response(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(solvent_vv_thermodynamics_type) :: bulk
      type(model_moz_3d_thermodynamics_type) :: thermo
      type(thermo_list_type) :: items
      type(huenenberger_correction_type) :: quantity
      real(wp) :: energy, record_energy, seed, shift, correction
      integer :: i
      character(len=4), parameter :: pressure_terms(3) = [character(len=4) :: "none", "PC", "PC+"]
      real(wp), parameter :: pressure_energy(3) = [0.0_wp, -0.5625_wp, -0.3125_wp]
      real(wp), parameter :: pressure_shift(3) = [0.0_wp, 1.5_wp, 0.8333333333333333_wp]

      bulk%pressure = 0.140625_wp
      bulk%pressure_ideal = 0.0625_wp
      bulk%compressibility = 10.666666666666666_wp
      correction = -0.5_wp*0.95_wp*2.380077363979553_wp/40.0_wp
      call new_huenenberger_correction(quantity, -1.0_wp, 0.95_wp, [40.0_wp, 40.0_wp, 40.0_wp], err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do i = 1, size(pressure_terms)
         call new_moz_3d_thermodynamics(thermo, 0.125_wp, 0.25_wp, bulk, 2.0_wp, err, items=items, &
            & correction=pressure_terms(i), finite_size=quantity)
         if (.not. allocated(err)) call items%energy(energy, err)
         if (.not. allocated(err)) call accumulate_solvation_adjoint(items, seed, shift, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, energy, 0.125_wp + pressure_energy(i) + correction, thr=1.0e-14_wp)
         if (allocated(error)) return
         call thermo%corrected_energy(pressure_terms(i), record_energy, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, record_energy, energy, thr=1.0e-14_wp)
         if (allocated(error)) return
         call check(error, seed, 1.0_wp)
         if (allocated(error)) return
         call check(error, shift, pressure_shift(i), thr=1.0e-14_wp)
         if (allocated(error)) return
      end do
      call items%select("huenenberger", .false., err)
      if (.not. allocated(err)) call items%energy(energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, energy, 0.125_wp - 0.3125_wp, thr=1.0e-14_wp)
   end subroutine check_finite_size_response

   !> Public selectors are normalized and invalid requests preserve the caller's energy
   subroutine check_corrected_energy(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_3d_thermodynamics_type) :: thermo
      character(len=8), parameter :: selectors(3) = [character(len=8) :: " NoNe ", " pC ", " pC+ "]
      character(len=8), parameter :: invalid(3) = [character(len=8) :: "PC++", "", "pressure"]
      real(wp), parameter :: expected(3) = [1.125_wp, 0.625_wp, 0.875_wp]
      real(wp) :: energy
      integer :: i

      thermo%free_energy = 1.0_wp
      thermo%pc = -0.5_wp
      thermo%pc_plus = -0.25_wp
      thermo%finite_size = 0.125_wp
      do i = 1, size(selectors)
         energy = 7.0_wp
         call thermo%corrected_energy(selectors(i), energy, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, energy, expected(i), thr=1.0e-14_wp, more="incorrect public corrected energy")
         if (allocated(error)) return
      end do
      do i = 1, size(invalid)
         energy = 7.0_wp
         call thermo%corrected_energy(invalid(i), energy, err)
         call check(error, allocated(err), more="an unknown public correction selector was accepted")
         if (allocated(error)) return
         call check(error, energy, 7.0_wp, thr=0.0_wp, more="a rejected correction changed the caller's energy")
         if (allocated(error)) return
      end do
   end subroutine check_corrected_energy

   !> Fixed scalar reference and the one-site ideal-gas limit
   subroutine check_reference_values(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: denominator, kappa, pressure, ideal, volume

      denominator = bulk_compressibility_denominator(0.125_wp, -4.0_wp)
      kappa = bulk_compressibility(0.125_wp, 2.0_wp, denominator)
      pressure = bulk_pressure(3, 0.125_wp, 2.0_wp, -4.0_wp)
      ideal = bulk_ideal_pressure(0.125_wp, 2.0_wp)
      volume = partial_molar_volume(0.25_wp, kappa, 2.0_wp)
      call check(error, denominator, 1.5_wp)
      if (allocated(error)) return
      call check(error, kappa, 10.666666666666666_wp)
      if (allocated(error)) return
      call check(error, pressure, 0.140625_wp)
      if (allocated(error)) return
      call check(error, ideal, 0.0625_wp)
      if (allocated(error)) return
      call check(error, volume, 4.0_wp)
      if (allocated(error)) return
      call check(error, pc_correction(pressure, volume), -0.5625_wp)
      if (allocated(error)) return
      call check(error, pc_plus_correction(pressure, ideal, volume), -0.3125_wp)
      if (allocated(error)) return
      call check(error, pc_moment_shift(pressure, kappa), 1.5_wp)
      if (allocated(error)) return
      call check(error, pc_plus_moment_shift(pressure, ideal, kappa), 0.8333333333333333_wp)
      if (allocated(error)) return

      denominator = bulk_compressibility_denominator(0.125_wp, 0.0_wp)
      kappa = bulk_compressibility(0.125_wp, 2.0_wp, denominator)
      pressure = bulk_pressure(1, 0.125_wp, 2.0_wp, 0.0_wp)
      volume = partial_molar_volume(0.0_wp, kappa, 2.0_wp)
      call check(error, pressure, ideal, more="one-site ideal-gas pressure")
      if (allocated(error)) return
      call check(error, volume, 8.0_wp, more="ideal-gas volume per molecule")
      if (allocated(error)) return
      call check(error, pc_plus_correction(pressure, ideal, volume), 0.0_wp)
   end subroutine check_reference_values

   !> A beta rescaling at fixed C and D preserves volumes and dimensionless shifts
   !> Pressure and correction energies carry exactly one inverse beta
   subroutine check_beta_scaling(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: rho = 0.007_wp, beta = 1.7_wp, factor = 11.0_wp
      real(wp) :: denominator, k0, k1, p0, p1, i0, i1, v0, v1

      denominator = bulk_compressibility_denominator(rho, -15.0_wp)
      k0 = bulk_compressibility(rho, beta, denominator)
      k1 = bulk_compressibility(rho, factor*beta, denominator)
      p0 = bulk_pressure(3, rho, beta, -15.0_wp)
      p1 = bulk_pressure(3, rho, factor*beta, -15.0_wp)
      i0 = bulk_ideal_pressure(rho, beta)
      i1 = bulk_ideal_pressure(rho, factor*beta)
      v0 = partial_molar_volume(-0.4_wp, k0, beta)
      v1 = partial_molar_volume(-0.4_wp, k1, factor*beta)
      call check(error, k1, factor*k0, rel=.true., thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, p1, p0/factor, rel=.true., thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, i1, i0/factor, rel=.true., thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, v1, v0, rel=.true., thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, pc_correction(p1, v1), pc_correction(p0, v0)/factor, rel=.true., thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, pc_plus_correction(p1, i1, v1), pc_plus_correction(p0, i0, v0)/factor, &
         & rel=.true., thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, pc_moment_shift(p1, k1), pc_moment_shift(p0, k0), rel=.true., thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, pc_plus_moment_shift(p1, i1, k1), pc_plus_moment_shift(p0, i0, k0), &
         & rel=.true., thr=1.0e-14_wp)
   end subroutine check_beta_scaling

   !> Public corrected energies differentiated with respect to the UV moment D
   !> Check beta dE/dD against the response shift over independent beta and D values
   subroutine check_moment_derivatives(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(solvent_vv_thermodynamics_type) :: bulk
      type(model_moz_3d_thermodynamics_type) :: plus, minus
      character(len=4), parameter :: corrections(3) = ["none", "PC  ", "PC+ "]
      real(wp), parameter :: betas(4) = [0.7_wp, 1.7_wp, 17.0_wp, 170.0_wp]
      real(wp), parameter :: moments(4) = [-0.6_wp, 0.1_wp, 0.9_wp, 1.3_wp], step = 1.0e-5_wp
      real(wp) :: fd, shift, energy_plus, energy_minus
      integer :: ib, im, ic

      bulk%pressure = 0.31_wp
      bulk%pressure_ideal = 0.043_wp
      bulk%compressibility = 23.7_wp
      do ib = 1, size(betas)
         do im = 1, size(moments)
            call new_moz_3d_thermodynamics(plus, 0.125_wp, moments(im) + step, bulk, betas(ib), err)
            if (.not. allocated(err)) then
               call new_moz_3d_thermodynamics(minus, 0.125_wp, moments(im) - step, bulk, betas(ib), err)
            end if
            call check_moist_error(error, err)
            if (allocated(error)) return
            do ic = 1, size(corrections)
               call plus%corrected_energy(corrections(ic), energy_plus, err)
               if (.not. allocated(err)) call minus%corrected_energy(corrections(ic), energy_minus, err)
               call check_moist_error(error, err)
               if (allocated(error)) return
               fd = betas(ib)*(energy_plus - energy_minus)/(2.0_wp*step)
               shift = correction_moment_shift(corrections(ic), bulk)
               call check(error, shift, fd, thr=2.0e-9_wp, &
                  & more="beta dE/dD differs for "//trim(corrections(ic)))
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine check_moment_derivatives

   !> Independent list construction and reverse passes remain safe under concurrent calls
   subroutine check_parallel_lists(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: success, passed
      integer :: i

      success = .true.
      !$omp parallel do num_threads(4) schedule(static) private(passed) reduction(.and.:success)
      do i = 1, 128
         call evaluate_list(i, passed)
         success = success .and. passed
      end do
      !$omp end parallel do
      call check(error, success, more="concurrent thermodynamic lists lost their values or sensitivities")

   contains

      subroutine evaluate_list(index, passed)
         integer, intent(in) :: index
         logical, intent(out) :: passed
         type(moist_error), allocatable :: err
         type(solvent_vv_thermodynamics_type) :: bulk
         type(model_moz_3d_thermodynamics_type) :: thermo
         type(thermo_list_type) :: items
         real(wp) :: energy, seed, shift, base, expected, expected_shift

         passed = .false.
         base = real(index, wp)/256.0_wp
         bulk%pressure = 0.140625_wp
         bulk%pressure_ideal = 0.0625_wp
         bulk%compressibility = 10.666666666666666_wp
         expected = base
         expected_shift = 0.0_wp
         if (mod(index, 2) == 0) then
            call new_moz_3d_thermodynamics(thermo, base, 0.25_wp, bulk, 2.0_wp, err, &
               & items=items, correction="PC+")
            expected = base - 0.3125_wp
            expected_shift = 0.8333333333333333_wp
         else
            call new_moz_3d_thermodynamics(thermo, base, 0.25_wp, bulk, 2.0_wp, err, &
               & items=items, correction="none")
         end if
         if (.not. allocated(err)) call items%energy(energy, err)
         if (.not. allocated(err)) call accumulate_solvation_adjoint(items, seed, shift, err)
         if (allocated(err)) return
         passed = abs(energy - expected) < 1.0e-13_wp .and. &
            & abs(seed - 1.0_wp) < 1.0e-13_wp .and. abs(shift - expected_shift) < 1.0e-13_wp
      end subroutine evaluate_list

   end subroutine check_parallel_lists

   !> Shared nodes are evaluated once, and all selected branches accumulate before reverse mode
   subroutine check_shared_dependencies(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(thermo_list_type) :: items
      type(linear_item_type) :: item
      real(wp) :: energy, derivative

      ! Deliberately register dependents before their common dependency.
      item%name = "a"
      item%unit = "Hartree"
      item%dependencies = ["shared"]
      item%coefficients = [2.0_wp]
      item%correction = .true.
      call items%add(item, err, selected=.true.)
      item%name = "b"
      item%coefficients = [-0.5_wp]
      if (.not. allocated(err)) call items%add(item, err, selected=.true.)
      item%name = "diagnostic"
      item%coefficients = [100.0_wp]
      if (.not. allocated(err)) call items%add(item, err)
      item%name = "shared"
      item%dependencies = ["x"]
      item%coefficients = [3.0_wp]
      item%correction = .false.
      item%once = .true.
      if (.not. allocated(err)) call items%add(item, err)
      if (.not. allocated(err)) call items%add_input("x", 2.0_wp, err)
      if (.not. allocated(err)) call items%add_input("base", 0.125_wp, err, unit="Hartree")
      if (.not. allocated(err)) call items%set_base_energy("base", err)
      if (.not. allocated(err)) call items%evaluate(err)
      ! A second evaluation must reuse the common node rather than call it again.
      if (.not. allocated(err)) call items%evaluate(err)
      if (.not. allocated(err)) call items%energy(energy, err)
      if (.not. allocated(err)) call items%accumulate_adjoint(err)
      if (.not. allocated(err)) call items%get_adjoint("x", derivative, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, energy, 9.125_wp, more="diagnostic contributed to selected energy")
      if (allocated(error)) return
      call check(error, derivative, 4.5_wp, more="shared dependency lost an accumulated contribution")
      if (allocated(error)) return
      call items%get_adjoint("diagnostic", derivative, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, derivative, 0.0_wp, more="diagnostic seeded the reverse pass")
      if (allocated(error)) return
      call items%accumulate_adjoint(err)
      if (.not. allocated(err)) call items%get_adjoint("x", derivative, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, derivative, 4.5_wp, more="a repeated reverse pass retained old sensitivities")
   end subroutine check_shared_dependencies

   !> PC/PC+ selection is exclusive, diagnostics are inert and copies own their state
   subroutine check_selection(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(thermo_list_type) :: items, copy
      type(linear_item_type) :: correction
      real(wp) :: energy, original, seed, shift, derivative

      call new_uv_thermo_list(items, 0.125_wp, 0.25_wp, 0.140625_wp, 0.0625_wp, &
         & 10.666666666666666_wp, 2.0_wp, err)
      if (.not. allocated(err)) call items%evaluate(err)
      if (.not. allocated(err)) call items%energy(energy, err)
      if (.not. allocated(err)) call accumulate_solvation_adjoint(items, seed, shift, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, energy, 0.125_wp)
      if (allocated(error)) return
      call check(error, seed, 1.0_wp)
      if (allocated(error)) return
      call check(error, shift, 0.0_wp, more="unselected corrections changed the response")
      if (allocated(error)) return

      call items%select(" PC ", .true., err)
      if (.not. allocated(err)) call items%energy(original, err)
      if (.not. allocated(err)) call accumulate_solvation_adjoint(items, seed, shift, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, original, -0.4375_wp)
      if (allocated(error)) return
      call check(error, shift, 1.5_wp)
      if (allocated(error)) return
      call items%select("pc+", .true., err)
      call check(error, allocated(err), more="simultaneous PC and PC+ were accepted")
      if (allocated(error)) return
      deallocate (err)
      call items%energy(energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, energy, original, more="a rejected selection mutated the active corrections")
      if (allocated(error)) return
      copy = items
      call copy%select("pc", .false., err)
      if (.not. allocated(err)) call copy%select("pc+", .true., err)
      if (.not. allocated(err)) call accumulate_solvation_adjoint(copy, seed, shift, err)
      if (.not. allocated(err)) call copy%energy(energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, energy, -0.1875_wp)
      if (allocated(error)) return
      call check(error, shift, 0.8333333333333333_wp)
      if (allocated(error)) return
      call items%energy(energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, energy, original, more="a copied list changed its source")
      if (allocated(error)) return

      ! A second selected branch can depend on the raw energy; its seed must accumulate.
      correction%name = "extra"
      correction%unit = "Hartree"
      correction%dependencies = ["free_energy"]
      correction%coefficients = [2.0_wp]
      correction%correction = .true.
      if (.not. allocated(error)) call copy%add(correction, err, selected=.true.)
      if (.not. allocated(err)) call copy%evaluate(err)
      if (.not. allocated(err)) call accumulate_solvation_adjoint(copy, seed, shift, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, seed, 3.0_wp, more="raw energy sensitivity was not accumulated")
      if (allocated(error)) return
      call copy%set_input("uv_moment", 0.5_wp, err)
      if (.not. allocated(err)) call copy%get_adjoint("uv_moment", derivative, err)
      call check(error, allocated(err), more="a changed input retained stale sensitivities")
      if (allocated(error)) return
      deallocate (err)
      call copy%evaluate(err)
      if (.not. allocated(err)) call copy%get_value("partial_molar_volume", derivative, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, derivative, 2.6666666666666665_wp)
   end subroutine check_selection

   !> All primitive sensitivities of the bulk and UV graphs agree with finite differences
   subroutine check_graph_derivatives(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(thermo_list_type) :: items, plus, minus
      type(linear_item_type) :: observable
      character(len=16), parameter :: bulk_names(3) = [character(len=16) :: "rho", "beta", "bulk_moment"]
      real(wp), parameter :: bulk_values(3) = [0.125_wp, 2.0_wp, -4.0_wp]
      character(len=16), parameter :: uv_names(6) = [character(len=16) :: &
         & "free_energy", "uv_moment", "pressure", "pressure_ideal", "compressibility", "beta"]
      real(wp), parameter :: uv_values(6) = [0.125_wp, 0.25_wp, 0.140625_wp, 0.0625_wp, 10.666666666666666_wp, 2.0_wp]
      character(len=4), parameter :: corrections(3) = ["none", "PC  ", "PC+ "]
      integer :: i, ic

      call new_bulk_thermo_list(items, 3, bulk_values(1), bulk_values(2), bulk_values(3), err)
      observable%name = "observable"
      observable%unit = "Hartree"
      observable%dependencies = [character(len=32) :: "pressure", "pressure_ideal", "compressibility"]
      observable%coefficients = [0.3_wp, 0.2_wp, 0.05_wp]
      observable%correction = .true.
      if (.not. allocated(err)) call items%add(observable, err, selected=.true.)
      if (.not. allocated(err)) call items%evaluate(err)
      if (.not. allocated(err)) call items%accumulate_adjoint(err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do i = 1, size(bulk_names)
         call compare_derivative(bulk_names(i), bulk_values(i))
         if (allocated(error)) return
      end do
      do ic = 1, size(corrections)
         call new_uv_thermo_list(items, uv_values(1), uv_values(2), uv_values(3), uv_values(4), &
            & uv_values(5), uv_values(6), err, corrections=[corrections(ic)])
         if (.not. allocated(err)) call items%evaluate(err)
         if (.not. allocated(err)) call items%accumulate_adjoint(err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         do i = 1, size(uv_names)
            call compare_derivative(uv_names(i), uv_values(i))
            if (allocated(error)) return
         end do
      end do

   contains

      subroutine compare_derivative(name, value)
         character(len=*), intent(in) :: name
         real(wp), intent(in) :: value
         real(wp) :: step, ep, em, derivative, fd
         step = 1.0e-6_wp*max(abs(value), 1.0_wp)
         plus = items
         minus = items
         call plus%set_input(name, value + step, err)
         if (.not. allocated(err)) call minus%set_input(name, value - step, err)
         if (.not. allocated(err)) call plus%evaluate(err)
         if (.not. allocated(err)) call minus%evaluate(err)
         if (.not. allocated(err)) call plus%energy(ep, err)
         if (.not. allocated(err)) call minus%energy(em, err)
         if (.not. allocated(err)) call items%get_adjoint(name, derivative, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         fd = (ep - em)/(2.0_wp*step)
         call check(error, derivative, fd, thr=1.0e-7_wp*max(abs(fd), 1.0_wp), &
            & more="reverse sensitivity differs for "//trim(name))
      end subroutine compare_derivative

   end subroutine check_graph_derivatives

   subroutine check_list_errors(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(thermo_list_type) :: items, empty
      type(linear_item_type) :: item
      real(wp) :: value

      call empty%evaluate(err)
      if (.not. allocated(err)) call empty%accumulate_adjoint(err)
      if (.not. allocated(err)) call empty%energy(value, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, value, 0.0_wp)
      if (allocated(error)) return
      call new_uv_thermo_list(items, 0.0_wp, 0.0_wp, 0.1_wp, 0.02_wp, 5.0_wp, 2.0_wp, &
         & err, corrections=[character(len=4) :: "PC", "PC+"])
      call expect_error("PC and PC+ factory selection was accepted")
      if (allocated(error)) return
      call items%get_value("pc", value, err)
      call expect_error("unevaluated values were available")
      if (allocated(error)) return
      call items%select("pressure", .true., err)
      call expect_error("diagnostic quantity was selected as an energy")
      if (allocated(error)) return
      call items%set_base_energy("pressure", err)
      call expect_error("a pressure was accepted as a Hartree base energy")
      if (allocated(error)) return
      call items%add_input(" PC ", 0.0_wp, err)
      call expect_error("a duplicate name was accepted")
      if (allocated(error)) return
      call items%select("missing", .true., err)
      call expect_error("an unknown correction was accepted")
      if (allocated(error)) return
      call items%set_input("pc", 0.0_wp, err)
      call expect_error("a computed item was changed as an input")
      if (allocated(error)) return

      item%name = "missing_dependency"
      item%dependencies = ["absent"]
      item%coefficients = [1.0_wp]
      item%correction = .true.
      item%unit = "bohr^3"
      call empty%add(item, err, selected=.true.)
      call expect_error("a volume was accepted as an energy correction")
      if (allocated(error)) return
      item%correction = .false.
      call empty%add(item, err)
      if (.not. allocated(err)) call empty%evaluate(err)
      call expect_error("a missing dependency was accepted")
      if (allocated(error)) return
      call empty%add_input("absent", 1.0_wp, err)
      if (.not. allocated(err)) call empty%evaluate(err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call empty%set_input("absent", ieee_value(0.0_wp, ieee_quiet_nan), err)
      if (.not. allocated(err)) call empty%evaluate(err)
      call expect_error("a non-finite input was accepted")
      if (allocated(error)) return

      call new_bulk_thermo_list(items, 3, 0.125_wp, 2.0_wp, 8.0_wp, err)
      if (.not. allocated(err)) call items%evaluate(err)
      call expect_error("zero compressibility denominator was accepted")
      if (allocated(error)) return
      item%name = "cycle_a"
      item%dependencies = ["cycle_b"]
      call items%add(item, err)
      item%name = "cycle_b"
      item%dependencies = ["cycle_a"]
      if (.not. allocated(err)) call items%add(item, err)
      if (.not. allocated(err)) call items%evaluate(err)
      call expect_error("a dependency cycle was accepted")

   contains

      subroutine expect_error(message)
         character(len=*), intent(in) :: message
         call check(error, allocated(err), more=message)
         if (allocated(err)) deallocate (err)
      end subroutine expect_error

   end subroutine check_list_errors

   subroutine linear_evaluate(self, values, value, error)
      class(linear_item_type), intent(inout) :: self
      real(wp), intent(in) :: values(:)
      real(wp), intent(out) :: value
      type(moist_error), allocatable, intent(out) :: error
      self%evaluations = self%evaluations + 1
      value = dot_product(self%coefficients, values)
      if (self%once .and. self%evaluations > 1) call fatal_error(error, "shared dependency evaluated twice")
   end subroutine linear_evaluate

   pure subroutine linear_accumulate(self, values, value, seed, adjoints)
      class(linear_item_type), intent(in) :: self
      real(wp), intent(in) :: values(:), value, seed
      real(wp), intent(inout) :: adjoints(:)
      adjoints = adjoints + seed*self%coefficients
   end subroutine linear_accumulate

end module test_moz_thermodynamics
