!> Unit tests for the pressure-volume solvation-model component
!>
!> We are testing
!>   * the energy is `p * V`, with `V` the enclosed cavity volume
!>   * the nuclear gradient is `p * dV/dR`
!>   * the surface adjoints are `p *` the cavity volume adjoints
module test_model_component_pv
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_type, only: coupling_type, potential_type
   use moist_model_components, only: pv, new_pv
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_cavity_drop, only: cavity_type_drop
   use moist_cavity_numsa, only: cavity_type_numsa, new_cavity_numsa
   use moist_radii, only: radius_type, static_radius_type, new_cosmo_radii, &
      & new_radii_custom_atoms
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: get_test_structures, center_at_origin, &
      & get_test_cavity_iswig, fd4_scalar, fd4_offsets
   use test_model_component_helper, only: surface_fixture, &
      & new_surface_fixture, check_surface_weights, fixture_radial_normals, &
      & ngrid_sw => fixture_ngrid_param, sw_areas => fixture_areas_param, &
      & sw_xis => fixture_xis_param, sw_fs => fixture_fs_param, &
      & sw_xyz => fixture_xyz_param
   implicit none (type, external)
   private

   public :: collect_model_component_pv

   !> Tolerance for values that must agree to roundoff
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)

contains

!> Collect the PV component test suite.
!>
!> @param[out] testsuite Collected unit tests
   subroutine collect_model_component_pv(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("pv_sphere_volume", test_pv_sphere_volume), &
         & new_unittest("pv_nuclear_gradient", test_pv_nuclear_gradient), &
         & new_unittest("pv_surface_weights", test_pv_surface_weights), &
         & new_unittest("pv_zero_pressure_short_circuit", test_pv_short_circuit), &
         & new_unittest("pv_lifecycle_guards", test_pv_guards) &
         & ]

   end subroutine collect_model_component_pv

!> The PV energy of a spherical cavity against the analytic `p * 4/3 pi R^3`
!>
!> A one-atom iSwiG surface carries no switching contributions, so its
!> divergence-theorem volume `1/3 sum_i a_i (n_i . r_i)` collapses to the exact
!> sphere volume at any Lebedev order and any sphere position -- the sphere is
!> deliberately placed off the origin so a dropped center offset in that
!> contraction would show up. The pressures span ambient to extreme; one atomic
!> unit of pressure is about 29.4 TPa, so the top of the sweep is far above any
!> physical solvation pressure and makes the PV term dominate outright.
!>
!> @param[out] error Error handling
   subroutine test_pv_sphere_volume(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Single-sphere and two-sphere test structures
      type(structure_type) :: mol, mol_pair
      !> Spherical test cavity
      type(cavity_type_iswig) :: cavity
      !> Radius model pinning the sphere radius exactly
      class(radius_type), allocatable :: radius_model
      !> Component under test
      type(pv) :: pv_component
      !> Host coupling data, never read by PV
      type(coupling_type) :: coupling
      !> Energy accumulator and the analytic reference volume
      real(wp) :: energy, volume_ref
      !> Failure context
      character(len=128) :: context
      !> Radius, Lebedev-order and pressure indices
      integer :: irad, ileb, ipres

      !> Sphere radii to sweep, bohr
      real(wp), parameter :: test_radii(*) = [1.5_wp, 3.0_wp, 5.5_wp]
      !> Lebedev orders to sweep
      integer, parameter :: test_nlebs(*) = [26, 110]
      !> Pressures to sweep, atomic units: zero, 1 atm, ~29 GPa, extreme
      real(wp), parameter :: test_pressures(*) = &
         & [0.0_wp, 3.4e-9_wp, 1.0e-3_wp, 1.0_wp]
      !> Preloaded accumulator value, so an assignment cannot pass for an add
      real(wp), parameter :: sentinel = 0.125_wp
      !> Deliberately off-origin sphere center, bohr
      real(wp), parameter :: center(3) = [0.7_wp, -1.3_wp, 2.1_wp]
      !> Separation of the two spheres in the additivity check, bohr
      real(wp), parameter :: separation = 60.0_wp
      !> Absolute and relative tolerances of the analytic comparison
      real(wp), parameter :: vol_atol = 1.0e-12_wp
      real(wp), parameter :: vol_rtol = 1.0e-11_wp

      call new (mol, [6], reshape(center, [3, 1]))

      do irad = 1, size(test_radii)
         volume_ref = 4.0_wp/3.0_wp*pi*test_radii(irad)**3

         call new_radii_custom_atoms([test_radii(irad)], radius_model, err)
         if (allocated(err)) then
            call test_failed(error, "Radius model setup failed: "//err%message)
            return
         end if

         do ileb = 1, size(test_nlebs)
            call get_test_cavity_iswig(mol, cavity, err, nleb=test_nlebs(ileb), &
               & radius_model=radius_model)
            if (allocated(err)) then
               call test_failed(error, "Spherical cavity setup failed: "//err%message)
               return
            end if

            do ipres = 1, size(test_pressures)
               call new_pv(pv_component, test_pressures(ipres))
               call pv_component%update(mol, cavity, err)
               if (allocated(err)) then
                  call test_failed(error, "PV update failed: "//err%message)
                  return
               end if

               energy = sentinel
               call pv_component%get_energy(coupling, cavity, energy, err)
               if (allocated(err)) then
                  call test_failed(error, "PV energy failed: "//err%message)
                  return
               end if

               write (context, '(a,f0.2,a,i0,a,es9.2)') &
                  & "PV sphere R = ", test_radii(irad), ", nleb = ", &
                  & test_nlebs(ileb), ", p = ", test_pressures(ipres)
               call check(error, energy - sentinel, &
                  & test_pressures(ipres)*volume_ref, &
                  & thr_abs=vol_atol, thr_rel=vol_rtol, more=trim(context))
               if (allocated(error)) return
            end do
         end do
      end do

      ! The sweep is only meaningful if the top pressure moved the accumulator.
      call check(error, abs(energy - sentinel) > 1.0_wp, &
         & more="PV sphere energy is negligible, the sweep is vacuous")
      if (allocated(error)) return

      ! Two spheres far enough apart not to switch each other off must give
      ! twice the single-sphere volume, so PV really reads the cavity total.
      call new (mol_pair, [6, 6], reshape([ &
                                          center, center + [separation, 0.0_wp, 0.0_wp]], [3, 2]))
      call new_radii_custom_atoms([test_radii(2), test_radii(2)], radius_model, err)
      if (allocated(err)) then
         call test_failed(error, "Pair radius model setup failed: "//err%message)
         return
      end if
      call get_test_cavity_iswig(mol_pair, cavity, err, nleb=test_nlebs(1), &
         & radius_model=radius_model)
      if (allocated(err)) then
         call test_failed(error, "Two-sphere cavity setup failed: "//err%message)
         return
      end if

      call new_pv(pv_component, test_pressures(size(test_pressures)))
      call pv_component%update(mol_pair, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "PV update on the two-sphere cavity failed: "//err%message)
         return
      end if
      energy = 0.0_wp
      call pv_component%get_energy(coupling, cavity, energy, err)
      if (allocated(err)) then
         call test_failed(error, "PV energy on the two-sphere cavity failed: "//err%message)
         return
      end if
      call check(error, energy, &
         & test_pressures(size(test_pressures))*2.0_wp*4.0_wp/3.0_wp*pi*test_radii(2)**3, &
         & thr_abs=vol_atol, thr_rel=vol_rtol, &
         & more="PV energy of two separated spheres is not twice the single-sphere term")

   end subroutine test_pv_sphere_volume

!> The PV nuclear gradient against fourth-order central differences.
!>
!> The cavity is rebuilt from scratch at every stencil point, so this differences
!> the same quantity the component reports rather than a frozen surface. The
!> gradient accumulator carries a sentinel throughout, which turns an assignment
!> where an accumulation was meant into a failure; the exact pressure scaling is
!> then checked separately, where it can be asserted to roundoff instead of to
!> finite-difference accuracy.
!>
!> @param[out] error Error handling
   subroutine test_pv_nuclear_gradient(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Sampled test structures and the displaced copy driven through the cavity
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: trial
      !> Cavity rebuilt at the reference and at every displaced geometry
      type(cavity_type_iswig) :: cavity
      !> Components at unit and at scaled pressure
      type(pv) :: pv_component, pv_scaled
      !> Host coupling data, never read by PV
      type(coupling_type) :: coupling
      !> Gradient accumulators at unit and at scaled pressure
      real(wp), allocatable :: gradient(:, :), gradient_scaled(:, :)
      !> Stencil samples, in `fd4_offsets` order, and the resulting derivative
      real(wp) :: values(4), fd
      !> Saved reference coordinate
      real(wp) :: saved
      !> Failure context
      character(len=128) :: context
      !> Atom, axis and stencil indices
      integer :: iatom, iaxis, k

      !> Lebedev grid size
      integer, parameter :: nleb = 26
      !> Nuclear-displacement step in bohr
      real(wp), parameter :: step = 1.0e-3_wp
      !> Tolerances accounting for the cavity rebuild roundoff
      real(wp), parameter :: fd_atol = 5.0e-9_wp
      real(wp), parameter :: fd_rtol = 5.0e-8_wp
      !> Reference pressure of the finite-difference comparison
      real(wp), parameter :: unit_pressure = 1.0_wp
      !> Second pressure, for the scaling check
      real(wp), parameter :: scaled_pressure = 3.25_wp
      !> Preloaded accumulator value, so an assignment cannot pass for an add
      real(wp), parameter :: sentinel = 0.5_wp

      call get_test_structures(mols, 5)
      call center_at_origin(mols(1))

      call get_test_cavity_iswig(mols(1), cavity, err, nleb=nleb)
      if (allocated(err)) then
         call test_failed(error, "Reference cavity setup failed: "//err%message)
         return
      end if

      call new_pv(pv_component, unit_pressure)
      call pv_component%update(mols(1), cavity, err)
      if (allocated(err)) then
         call test_failed(error, "PV update failed: "//err%message)
         return
      end if

      allocate (gradient(3, mols(1)%nat), source=sentinel)
      call pv_component%get_gradient(coupling, cavity, gradient, err)
      if (allocated(err)) then
         call test_failed(error, "PV gradient failed: "//err%message)
         return
      end if
      call check(error, maxval(abs(gradient - sentinel)) > 1.0e-8_wp, &
         & more="PV nuclear gradient is zero, the test is vacuous")
      if (allocated(error)) return

      ! The pressure enters as a pure prefactor, so this holds to roundoff.
      call new_pv(pv_scaled, scaled_pressure)
      call pv_scaled%update(mols(1), cavity, err)
      if (allocated(err)) then
         call test_failed(error, "Scaled PV update failed: "//err%message)
         return
      end if
      allocate (gradient_scaled(3, mols(1)%nat), source=sentinel)
      call pv_scaled%get_gradient(coupling, cavity, gradient_scaled, err)
      if (allocated(err)) then
         call test_failed(error, "Scaled PV gradient failed: "//err%message)
         return
      end if
      call check(error, maxval(abs((gradient_scaled - sentinel) &
         & - scaled_pressure*(gradient - sentinel)/unit_pressure)), 0.0_wp, thr=thr, &
         & message="PV nuclear gradient does not scale linearly with the pressure")
      if (allocated(error)) return

      do iatom = 1, min(2, mols(1)%nat)
         do iaxis = 1, 3
            saved = mols(1)%xyz(iaxis, iatom)
            do k = 1, size(fd4_offsets)
               trial = mols(1)
               trial%xyz(iaxis, iatom) = saved + fd4_offsets(k)*step
               call displaced_energy(trial, values(k))
               if (allocated(error)) return
            end do
            fd = fd4_scalar(values(1), values(2), values(3), values(4), step)
            write (context, '(a,i0,a,i0)') "PV gradient atom ", iatom, ", axis ", iaxis
            call check(error, gradient(iaxis, iatom) - sentinel, fd, &
               & thr_abs=fd_atol, thr_rel=fd_rtol, more=trim(context))
            if (allocated(error)) return
         end do
      end do

   contains

      !> Evaluate the PV energy for one displaced molecular structure.
      !>
      !> @param[in]  displaced_mol Displaced structure
      !> @param[out] energy        PV energy on the rebuilt cavity
      subroutine displaced_energy(displaced_mol, energy)

         !> Displaced structure
         type(structure_type), intent(in) :: displaced_mol

         !> PV energy on the rebuilt cavity
         real(wp), intent(out) :: energy

         type(moist_error_type), allocatable :: local_err

         energy = 0.0_wp
         call get_test_cavity_iswig(displaced_mol, cavity, local_err, nleb=nleb)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced cavity setup failed: "//local_err%message)
            return
         end if
         call pv_component%update(displaced_mol, cavity, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced PV update failed: "//local_err%message)
            return
         end if
         call pv_component%get_energy(coupling, cavity, energy, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced PV energy failed: "//local_err%message)
            return
         end if

      end subroutine displaced_energy

   end subroutine test_pv_nuclear_gradient

!> The PV surface adjoints against fourth-order central differences.
!>
!> Driven on the shared synthetic seven-point DROP surface, carrying the radial
!> normal field: the volume is the divergence-theorem integral
!> `1/3 sum_i a_i (r_i . n_i)`, which only measures an enclosed volume on a
!> closed radial surface, so the harness's tilted field would change what is
!> being tested. Unlike CPCM, the volume depends on the normals, so all four
!> channels stay enabled. The pressure is deliberately not one, so a dropped
!> prefactor cannot hide.
!>
!> @param[out] error Error handling
   subroutine test_pv_surface_weights(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Dummy structure; PV only stores it
      type(structure_type) :: mol
      !> Synthetic DROP surface carrying the volume adjoint fields
      type(cavity_type_drop) :: cavity
      !> Components at finite and at zero pressure
      type(pv) :: pv_component, pv_zero
      !> Host coupling data, never read by PV
      type(coupling_type) :: coupling
      !> Analytic surface weights, and a prefilled accumulator PV must not touch
      type(cavity_surface_adjoint_type) :: weights, prefilled
      !> Fixture mirroring the synthetic surface for the harness
      type(surface_fixture) :: surface
      !> Radial normal field of the fixture
      real(wp) :: normals(3, ngrid_sw)
      !> Dummy molecular geometry
      real(wp) :: xyz_mol(3, 1)
      !> Grid point index
      integer :: igrid

      !> Pressure of the component under test
      real(wp), parameter :: pressure = 2.5_wp
      !> Finite-difference displacement, near the roundoff/truncation optimum
      real(wp), parameter :: step = 3.0e-4_wp
      !> Tolerances of the finite-difference comparison
      real(wp), parameter :: fd_atol = 2.0e-10_wp
      real(wp), parameter :: fd_rtol = 2.0e-9_wp
      !> Value prefilled into every channel of the untouched accumulator
      real(wp), parameter :: prefill = 0.75_wp

      normals = fixture_radial_normals()

      xyz_mol(:, 1) = 0.0_wp
      call new (mol, [1], xyz_mol)

      ! Synthetic surface carrying only the fields the volume adjoint reads,
      ! plus the total volume PV's own lifecycle guards require.
      cavity%ngrid = ngrid_sw
      cavity%nsph = 1
      allocate (cavity%a, source=sw_areas)
      allocate (cavity%xi0, source=sw_xis)
      allocate (cavity%f, source=sw_fs)
      allocate (cavity%xyz, source=sw_xyz)
      allocate (cavity%normal0, source=normals)
      allocate (cavity%total_volume)
      cavity%total_volume = 0.0_wp
      do igrid = 1, ngrid_sw
         cavity%total_volume = cavity%total_volume &
            & + sw_areas(igrid)*dot_product(sw_xyz(:, igrid), normals(:, igrid))/3.0_wp
      end do

      call new_pv(pv_component, pressure)
      call pv_component%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "PV update failed: "//err%message)
         return
      end if

      call weights%init(ngrid_sw)
      call pv_component%get_surface_weights(coupling, cavity, weights, err)
      if (allocated(err)) then
         call test_failed(error, "PV surface-weight assembly failed: "//err%message)
         return
      end if
      call check(error, maxval(abs(weights%w_n)) > 0.0_wp, &
         & more="PV wrote no normal weights, the cavity hook is missing")
      if (allocated(error)) return

      call new_surface_fixture(surface, sw_areas, sw_xis, sw_fs, sw_xyz, normals)
      call check_surface_weights(error, surface, weights, pv_surface_energy, &
         & "pv", step=step, thr_abs=fd_atol, thr_rel=fd_rtol)
      if (allocated(error)) return

      ! A zero pressure must return before the accumulator is written at all.
      call prefilled%init(ngrid_sw)
      prefilled%w_xi = prefill
      prefilled%w_f = prefill
      prefilled%w_a = prefill
      prefilled%w_w = prefill
      prefilled%w_k1 = prefill
      prefilled%w_k2 = prefill
      prefilled%w_xyz = prefill
      prefilled%w_n = prefill

      call new_pv(pv_zero, 0.0_wp)
      call pv_zero%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "PV(0) update failed: "//err%message)
         return
      end if
      call pv_zero%get_surface_weights(coupling, cavity, prefilled, err)
      if (allocated(err)) then
         call test_failed(error, "PV(0) surface-weight assembly failed: "//err%message)
         return
      end if
      call check(error, max( &
         & maxval(abs(prefilled%w_xi - prefill)), maxval(abs(prefilled%w_f - prefill)), &
         & maxval(abs(prefilled%w_a - prefill)), maxval(abs(prefilled%w_w - prefill)), &
         & maxval(abs(prefilled%w_k1 - prefill)), maxval(abs(prefilled%w_k2 - prefill)), &
         & maxval(abs(prefilled%w_xyz - prefill)), maxval(abs(prefilled%w_n - prefill))), &
         & 0.0_wp, thr=0.0_wp, &
         & message="PV at zero pressure wrote to the surface-adjoint accumulator")

   contains

      !> Rebuild the PV energy independently on a perturbed surface.
      !>
      !> @param[in] trial Perturbed surface fixture
      !> @return PV energy of the perturbed surface
      function pv_surface_energy(trial) result(energy)

         !> Perturbed surface fixture
         type(surface_fixture), intent(in) :: trial

         !> PV energy of the perturbed surface
         real(wp) :: energy

         !> Quadrature areas reconstructed from the perturbed fixture
         real(wp), allocatable :: areas(:)
         !> Grid point index
         integer :: i

         areas = trial%areas()
         energy = 0.0_wp
         do i = 1, trial%ngrid()
            energy = energy &
               & + areas(i)*dot_product(trial%xyz(:, i), trial%normal(:, i))
         end do
         energy = pressure*energy/3.0_wp

      end function pv_surface_energy

   end subroutine test_pv_surface_weights

!> A zero pressure must short-circuit before the cavity is asked for anything.
!> Driven on a NUMSA cavity, which never fills the per-point volume derivatives,
!> so a missing short circuit is observable as the error PV raises without them.
!>
!> @param[out] error Error handling
   subroutine test_pv_short_circuit(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Molecular structure
      type(structure_type) :: mol
      !> Cavity without per-point volume derivatives
      type(cavity_type_numsa) :: cavity
      !> Radius model borrowed by the cavity
      type(static_radius_type) :: radius_model
      !> Host coupling data, never read by PV
      type(coupling_type) :: coupling
      !> Component under test
      type(pv) :: pv_component
      !> Gradient accumulator carrying a sentinel
      real(wp), allocatable :: gradient(:, :)

      !> Run context owned here and borrowed by the cavity
      type(moist_context_type), target :: ctx

      call new_context(ctx)
      call get_structure(mol, "MB16-43", "01")

      call new_cosmo_radii(radius_model)
      call new_cavity_numsa(cavity, ctx, nleb=110, radii=radius_model, error=err)
      if (allocated(err)) then
         call test_failed(error, "NUMSA cavity setup failed: "//err%message)
         return
      end if
      call cavity%update(mol, error=err)
      if (allocated(err)) then
         call test_failed(error, "NUMSA cavity update failed: "//err%message)
         return
      end if

      allocate (gradient(3, mol%nat), source=1.5_wp)

      ! Zero pressure: no cavity call at all, accumulator untouched.
      call new_pv(pv_component, 0.0_wp)
      call pv_component%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "PV(0) update failed: "//err%message)
         return
      end if
      call pv_component%get_gradient(coupling, cavity, gradient, err)
      call check(error, .not. allocated(err), &
         & more="PV at zero pressure queried the cavity volume gradient")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)
      call check(error, maxval(abs(gradient - 1.5_wp)), 0.0_wp, thr=thr, &
         & message="PV at zero pressure modified the gradient accumulator")
      if (allocated(error)) return

      ! Finite pressure: the missing cavity hook must surface as an error.
      call new_pv(pv_component, 0.75_wp)
      call pv_component%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "PV update failed: "//err%message)
         return
      end if
      call pv_component%get_gradient(coupling, cavity, gradient, err)
      call check(error, allocated(err), &
         & more="PV did not report the missing cavity volume-gradient hook")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)

   end subroutine test_pv_short_circuit

!> The lifecycle guards: a cavity that was never updated carries no volume, a
!> volume contribution produces no host-trace potential, and a mis-shaped
!> gradient accumulator is rejected without being written to.
!>
!> @param[out] error Error handling
   subroutine test_pv_guards(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Molecular structure
      type(structure_type) :: mol
      !> Cavity left unupdated, then updated for the shape guard
      type(cavity_type_iswig) :: cavity
      !> Radius model borrowed by the cavity
      type(static_radius_type) :: radius_model
      !> Host coupling data, never read by PV
      type(coupling_type) :: coupling
      !> Potential accumulator PV must leave alone
      type(potential_type) :: potential
      !> Component under test
      type(pv) :: pv_component
      !> Energy accumulator carrying a sentinel
      real(wp) :: energy
      !> Gradient accumulator of the wrong shape
      real(wp), allocatable :: bad_gradient(:, :)

      !> Pressure of the component under test
      real(wp), parameter :: pressure = 0.5_wp
      !> Preloaded accumulator value
      real(wp), parameter :: sentinel = 1.0_wp

      !> Run context owned here and borrowed by the cavity
      type(moist_context_type), target :: ctx

      call new_context(ctx)
      call get_structure(mol, "MB16-43", "01")
      call new_cosmo_radii(radius_model)

      call new_cavity_iswig(cavity, ctx, nleb=26, radius_model=radius_model, error=err)
      if (allocated(err)) then
         call test_failed(error, "Cavity construction failed: "//err%message)
         return
      end if

      ! Never updated: no total volume, so neither update nor get_energy may run.
      call new_pv(pv_component, pressure)
      call pv_component%update(mol, cavity, err)
      call check(error, allocated(err), &
         & more="PV accepted a cavity that was never updated")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)

      energy = sentinel
      call pv_component%get_energy(coupling, cavity, energy, err)
      call check(error, allocated(err), &
         & more="PV returned an energy for a cavity without a volume")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)
      call check(error, energy, sentinel, thr=thr, &
         & message="rejected PV energy still modified the accumulator")
      if (allocated(error)) return

      call cavity%update(mol, error=err)
      if (allocated(err)) then
         call test_failed(error, "Cavity update failed: "//err%message)
         return
      end if
      call pv_component%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "PV update failed: "//err%message)
         return
      end if

      ! A volume contribution carries no host-trace potential.
      call pv_component%get_potential(coupling, cavity, potential, err)
      if (allocated(err)) then
         call test_failed(error, "PV potential failed: "//err%message)
         return
      end if
      call check(error, .not. allocated(potential%vat) &
         & .and. .not. allocated(potential%w_lsf0) &
         & .and. .not. allocated(potential%w_elstat_umol), &
         & more="PV wrote to the host potential")
      if (allocated(error)) return

      ! A gradient sized for the wrong number of nuclei.
      allocate (bad_gradient(3, mol%nat + 1), source=sentinel)
      call pv_component%get_gradient(coupling, cavity, bad_gradient, err)
      call check(error, allocated(err), more="a mis-shaped gradient was accepted")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)
      call check(error, maxval(abs(bad_gradient - sentinel)), 0.0_wp, thr=0.0_wp, &
         & more="rejected PV gradient left the accumulator untouched")

   end subroutine test_pv_guards

end module test_model_component_pv
