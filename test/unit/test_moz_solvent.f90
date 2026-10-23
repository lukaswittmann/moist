!> Tests of the MOZ bulk solvent: moz_solvent
module test_moz_solvent
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use mctc_io_codata2018, only: Boltzmann_constant, atomic_unit_of_energy
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_solvent_system, only: solvent_system_type, new_solvent_system
   use moist_math_packed_sym, only: packed_index, npair_from_ns
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent
   use moist_model_moz_solvent_water, only: new_spce_water_geometry
   use moist_model_moz_potential, only: moz_potential_type, potential_sites_type
   use moist_model_moz_potential_lj_base, only: lj_mixing_lorentz_berthelot, lj_mixing_geometric
   use moist_model_moz_potential_electrostatic_monopole, only: monopole_type
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use test_helpers, only: check_moist_error
   use test_moz_fixtures, only: water_id, water_sigma, water_epsilon, water_charges, new_water_system, &
      & new_water_solvent, new_cspce_water_solvent
   implicit none(type, external)
   private
   public :: collect_moz_solvent

contains
   !> Register tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_solvent(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("water_geometry", check_water_geometry), &
         & new_unittest("construct_from_table", check_construct_from_table), &
         & new_unittest("construct_custom_system", check_construct_custom), &
         & new_unittest("density_derivation", check_density_derivation), &
         & new_unittest("host_fed_refused", check_host_fed_refused), &
         & new_unittest("empty_potential_update", check_empty_potential_update), &
         & new_unittest("packed_site_pairs", check_site_pairs), &
         & new_unittest("vv_potential_tables", check_potential_tables), &
         & new_unittest("vv_tables_geometric_mixing", check_geometric_tables), &
         & new_unittest("solve_requires_update", check_solve_requires_update), &
         & new_unittest("state_point", check_state_point), &
         & new_unittest("cspce_water", check_cspce_water), &
         & new_unittest("stale_solution", check_stale_solution)]
   end subroutine collect_moz_solvent

   !> Check canonical water bonds and angle, and preservation of supplied library geometry
   !>
   !> @param[out] error Test failure
   subroutine check_water_geometry(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, supplied
      type(solvent_system_type) :: system
      type(solvent_vv_type) :: solvent
      real(wp) :: bond1(3), bond2(3), angle

      call new_spce_water_geometry(mol)
      call check(error, mol%nat == 3 .and. all(mol%num(mol%id) == [8, 1, 1]) .and. mol%charge == 0.0_wp)
      if (allocated(error)) return
      bond1 = mol%xyz(:, 2) - mol%xyz(:, 1)
      bond2 = mol%xyz(:, 3) - mol%xyz(:, 1)
      call check(error, norm2(bond1)/aatoau, 1.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, norm2(bond2)/aatoau, 1.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      angle = acos(dot_product(bond1, bond2)/(norm2(bond1)*norm2(bond2)))*180.0_wp/acos(-1.0_wp)
      call check(error, angle, 109.47_wp, thr=1.0e-12_wp)
      if (allocated(error)) return

      call new_context(ctx, nthreads=1, verbosity=0)
      call new_water_system(system, err)
      supplied = mol
      supplied%xyz = 1.1_wp*supplied%xyz
      system%solv_mol = supplied
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = solvent%structure()
      call check(error, all(mol%xyz == supplied%xyz), "Library construction must preserve explicitly supplied geometry")
   end subroutine check_water_geometry

   !> Check table water construction and invalid state points
   !>
   !> - Sites O, H, H; system id and temperature
   !> - Lorentz-Berthelot by default; unsolved state
   !> - Invalid state points and unknown mixing rule: refusal
   subroutine check_construct_from_table(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(solvent_vv_type) :: solvent
      type(structure_type) :: mol

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_water_system(system, err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%nsite(), 3)
      if (allocated(error)) return
      call check(error, solvent%solvent_id(), water_id)
      if (allocated(error)) return
      call check(error, solvent%potential%mixing(), lj_mixing_lorentz_berthelot)
      if (allocated(error)) return
      call check(error, associated(solvent%ctx, ctx), more="the solvent must borrow its context")
      if (allocated(error)) return
      call check(error, .not. solvent%is_solved() .and. .not. solvent%potential%is_updated(), &
         & more="a fresh solvent is solved or built")
      if (allocated(error)) return
      mol = solvent%structure()
      call check(error, all(mol%num(mol%id) == [8, 1, 1]), more="water sites are not O, H, H")
      if (allocated(error)) return
      ! Table geometry in Angstrom, solvent geometry in bohr; O-H about 0.96 A
      call check(error, norm2(mol%xyz(:, 2) - mol%xyz(:, 1)), 0.96_wp*aatoau, thr=0.05_wp*aatoau, &
         & more="O-H distance is not in bohr")
      if (allocated(error)) return

      call new_solvent_system(system, water_id, temperature=350.0_wp, error=err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      if (.not. allocated(err)) call solvent%potential%set_mixing(lj_mixing_geometric, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%potential%mixing(), lj_mixing_geometric)
      if (allocated(error)) return

      call solvent%potential%set_mixing(0, err)
      call check(error, allocated(err), more="an unknown mixing rule was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "mixing rule") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      system%temperature = -1.0_wp
      call new_vv_solvent(solvent, ctx, system, err)
      call check(error, allocated(err), more="a negative temperature was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "temperature must be positive") > 0, more=err%message)
   end subroutine check_construct_from_table

   !> Check custom solvent construction
   !>
   !> - Hand-filled system: id 0, supplied number density and structure in bohr
   !> - Missing or empty structure, zero density: refusal
   subroutine check_construct_custom(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_system_type) :: system
      type(structure_type) :: mol, copy

      call new_context(ctx, nthreads=0, verbosity=0)
      call new(mol, [8, 1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp, 0.0_wp, 0.0_wp, &
         & -0.5_wp, 1.7_wp, 0.0_wp], [3, 3]))
      system%solvent%id = 0
      system%temperature = 298.15_wp
      system%solvent%number_density_au = 5.0e-3_wp
      call new_vv_solvent(solvent, ctx, system, err)
      call check(error, allocated(err), more="a system without structure was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "no structure") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      system%solv_mol = mol
      call new_vv_solvent(solvent, ctx, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%solvent_id(), 0)
      if (allocated(error)) return
      call check(error, solvent%nsite(), 3)
      if (allocated(error)) return
      copy = solvent%structure()
      call check(error, all(copy%xyz == mol%xyz), more="the structure is not kept as given (bohr)")
      if (allocated(error)) return

      system%solvent%number_density_au = 0.0_wp
      call new_vv_solvent(solvent, ctx, system, err)
      call check(error, allocated(err), more="a zero density was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "density must be positive") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)
      system%solvent%number_density_au = 5.0e-3_wp
      call new(system%solv_mol, [integer ::], reshape([real(wp) ::], [3, 0]))
      call new_vv_solvent(solvent, ctx, system, err)
      call check(error, allocated(err), more="an empty structure was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "at least one atom") > 0, more=err%message)
   end subroutine check_construct_custom

   !> Check water number density in the source system
   !>
   !> - About 3.33e-2 molecules per A^3, 4.94e-3 per bohr^3; tolerance 1%
   subroutine check_density_derivation(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(solvent_vv_type) :: solvent
      real(wp), parameter :: water_per_bohr3 = 3.33e-2_wp/aatoau**3

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_water_system(system, err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, abs(water_per_bohr3 - 4.94e-3_wp) < 0.01_wp*4.94e-3_wp, &
         & more="reference conversion of the water density is wrong")
      if (allocated(error)) return
      call check(error, system%solvent%number_density_au, water_per_bohr3, thr=0.01_wp, rel=.true., &
         & more="water number density is off by more than 1%")
   end subroutine check_density_derivation

   !> Check host-term refusal at solvent add, built or unbuilt; unchanged potential
   subroutine check_host_fed_refused(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(solvent_vv_type) :: solvent
      type(moz_potential_type) :: pot
      type(monopole_type) :: host
      type(custom_lj_type) :: lj
      integer :: pass

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_water_system(system, err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      if (.not. allocated(err)) call new_custom_lj(lj, water_sigma, water_epsilon, err)
      if (.not. allocated(err)) call solvent%potential%add(lj, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      do pass = 1, 2
         if (pass == 2) then
            call host%update(solvent%structure(), err)
            call check_moist_error(error, err)
            if (allocated(error)) return
         end if
         call solvent%potential%add(host, err)
         call check(error, allocated(err), more="a host-fed term was added to the solvent")
         if (allocated(error)) return
         call check(error, index(err%message, "cannot serve a solvent") > 0, more=err%message)
         if (allocated(error)) return
         deallocate (err)
         pot = solvent%potential
         call check(error, pot%n_terms(), 1)
         if (allocated(error)) return
      end do
   end subroutine check_host_fed_refused

   !> Check build refusal without solvent terms
   subroutine check_empty_potential_update(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(solvent_vv_type) :: solvent, unconstructed

      call unconstructed%update(err)
      call check(error, allocated(err), more="an unconstructed solvent was built")
      if (allocated(error)) return
      call check(error, index(err%message, "Construct the 1D VV solvent") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_water_system(system, err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call solvent%update(err)
      call check(error, allocated(err), more="an empty solvent potential was updated")
      if (allocated(error)) return
      call check(error, index(err%message, "has no terms") > 0, more=err%message)
   end subroutine check_empty_potential_update

   !> Check packed VV upper-triangle pairs
   !>
   !> - Column packed_index(i, j): [i, j], i <= j, every pair once
   subroutine check_site_pairs(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      integer, allocatable :: pairs(:, :)
      integer :: i, j, ip

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_water_solvent(ctx, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call solvent%potential%site_pairs(pairs, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(pairs, 1), 2)
      if (allocated(error)) return
      call check(error, size(pairs, 2), npair_from_ns(3))
      if (allocated(error)) return
      call check(error, size(pairs, 2), 6)
      if (allocated(error)) return
      do j = 1, 3
         do i = 1, j
            ip = packed_index(i, j)
            call check(error, pairs(1, ip) == i .and. pairs(2, ip) == j, more="packed column holds the wrong pair")
            if (allocated(error)) return
         end do
      end do
   end subroutine check_site_pairs

   !> Check LJ and fixed-charge water VV tables at alpha = 1
   !>
   !> - Per packed pair: u_sr + ur_lr = u_LJ + q_i q_j/r
   !> - Real tail: ur_lr = q_i q_j erf(r)/r
   !> - Reciprocal tail: uk_lr = q_i q_j 4 pi exp(-k^2/4)/k^2
   subroutine check_potential_tables(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :), uk_lr(:, :)
      real(wp), allocatable :: sr6(:), expected(:), k2(:)
      real(wp) :: sigma, epsilon, qq
      real(wp), parameter :: pi = acos(-1.0_wp)
      integer :: i, j, ip

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_uniform_radial_pair(rgrid, kgrid, 64, 0.1_wp, err)
      if (.not. allocated(err)) call new_water_solvent(ctx, solvent, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call solvent%potential%compute(rgrid, kgrid, u_sr, ur_lr, uk_lr, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, all(shape(u_sr) == [rgrid%npts, 6]) .and. all(shape(ur_lr) == [rgrid%npts, 6]) &
         & .and. all(shape(uk_lr) == [kgrid%npts, 6]), more="VV tables are not (npts, npair)")
      if (allocated(error)) return
      k2 = kgrid%r**2
      do j = 1, 3
         do i = 1, j
            ip = packed_index(i, j)
            sigma = 0.5_wp*(water_sigma(i) + water_sigma(j))
            epsilon = sqrt(water_epsilon(i)*water_epsilon(j))
            qq = water_charges(i)*water_charges(j)
            sr6 = (sigma/rgrid%r)**6
            expected = 4.0_wp*epsilon*(sr6*sr6 - sr6) + qq/rgrid%r
            call check(error, maxval(abs(u_sr(:, ip) + ur_lr(:, ip) - expected)/max(1.0_wp, abs(expected))) &
               & < 1.0e-12_wp, more="u_sr + ur_lr is not the full VV potential")
            if (allocated(error)) return
            call check(error, maxval(abs(ur_lr(:, ip) - qq*erf(rgrid%r)/rgrid%r)) < 1.0e-12_wp, &
               & more="ur_lr is not the Ng real-space tail")
            if (allocated(error)) return
            call check(error, maxval(abs(uk_lr(:, ip) - qq*4.0_wp*pi*exp(-0.25_wp*k2)/k2)) < 1.0e-10_wp, &
               & more="uk_lr is not the Ng reciprocal tail")
            if (allocated(error)) return
         end do
      end do
   end subroutine check_potential_tables

   !> Check geometric mixing and split validation
   !>
   !> - u_sr + ur_lr = u_LJ(sqrt(s_i s_j), sqrt(e_i e_j)) + q_i q_j/r
   !> - Invalid split parameter: refusal only with split enabled
   subroutine check_geometric_tables(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :), uk_lr(:, :), sr6(:), expected(:)
      real(wp) :: sigma, epsilon, qq
      integer :: i, j, ip

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_uniform_radial_pair(rgrid, kgrid, 64, 0.1_wp, err)
      if (.not. allocated(err)) call new_water_solvent(ctx, solvent, err, mixing=lj_mixing_geometric)
      if (.not. allocated(err)) call solvent%potential%compute(rgrid, kgrid, u_sr, ur_lr, uk_lr, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do j = 1, 3
         do i = 1, j
            ip = packed_index(i, j)
            sigma = sqrt(water_sigma(i)*water_sigma(j))
            epsilon = sqrt(water_epsilon(i)*water_epsilon(j))
            qq = water_charges(i)*water_charges(j)
            sr6 = (sigma/rgrid%r)**6
            expected = 4.0_wp*epsilon*(sr6*sr6 - sr6) + qq/rgrid%r
            call check(error, maxval(abs(u_sr(:, ip) + ur_lr(:, ip) - expected)/max(1.0_wp, abs(expected))) &
               & < 1.0e-12_wp, more="geometric mixing is not used in the VV potential")
            if (allocated(error)) return
         end do
      end do

      call solvent%potential%set_ng_split(.true., -1.0_wp, err)
      call check(error, allocated(err), more="a negative alpha was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "alpha must be positive") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)
      call solvent%potential%set_ng_split(.false., -1.0_wp, err)
      call check_moist_error(error, err)
   end subroutine check_geometric_tables

   !> Check that the VV solve needs a constructed, updated solvent
   !>
   !> - Unconstructed and constructed but not updated solvents are refused
   !> - A refused solve leaves the solvent unsolved
   subroutine check_solve_requires_update(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(solvent_vv_type) :: solvent

      call solvent%solve(err)
      call check(error, allocated(err), more="an unconstructed solvent was solved")
      if (allocated(error)) return
      call check(error, index(err%message, "Construct") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_water_system(system, err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call solvent%solve(err)
      call check(error, allocated(err), more="a solvent without an updated potential was solved")
      if (allocated(error)) return
      call check(error, index(err%message, "Update") > 0, more=err%message)
      if (allocated(error)) return
      call check(error, .not. solvent%is_solved(), more="a refused solve marked the solvent solved")
   end subroutine check_solve_requires_update

   !> Check the state-point accessors
   !>
   !> - Zero before construction
   !> - Temperature and density of the system
   !> - beta = 1/(kB T) in 1/Hartree
   subroutine check_state_point(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(solvent_vv_type) :: solvent

      call check(error, solvent%temperature() == 0.0_wp .and. solvent%beta() == 0.0_wp .and. &
         & solvent%density() == 0.0_wp, &
         & more="an unconstructed solvent reports a state point")
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_solvent_system(system, water_id, temperature=350.0_wp, error=err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%temperature(), 350.0_wp)
      if (allocated(error)) return
      call check(error, solvent%beta(), atomic_unit_of_energy/(Boltzmann_constant*350.0_wp), rel=.true., &
         & thr=1.0e-14_wp, more="beta is not 1/(kB T) in 1/Hartree")
      if (allocated(error)) return
      call check(error, solvent%density(), system%solvent%number_density_au)
   end subroutine check_state_point

   !> Check the SPC/E water fixture
   !>
   !> - Updated with three sites
   !> - SPC/E charges, LJ on every site
   subroutine check_cspce_water(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(potential_sites_type) :: sites

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call solvent%potential%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%potential%is_updated() .and. sites%natom == 3, &
         & more="cSPC/E water is not updated with three sites")
      if (allocated(error)) return
      call check(error, allocated(sites%q) .and. allocated(sites%pair), more="cSPC/E water lacks charges or LJ")
      if (allocated(error)) return
      call check(error, maxval(abs(sites%q - [-0.8476_wp, 0.4238_wp, 0.4238_wp])), 0.0_wp, &
         & more="SPC/E charges differ")
      if (allocated(error)) return
      call check(error, all(sites%pair%sigma > 0.0_wp) .and. all(sites%pair%epsilon > 0.0_wp), &
         & more="cSPC/E LJ is missing on a site")
   end subroutine check_cspce_water

   !> Check that changes to the potential make a VV solution stale
   !>
   !> - Generation increases on add, update, Ng split and mixing changes
   !> - A refused setting leaves it unchanged
   !> - The model copy carries the generation
   subroutine check_stale_solution(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_type), allocatable :: copy
      integer(int64) :: generation

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cspce_water_solvent(ctx, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      generation = solvent%potential%generation()
      call solvent%potential%set_mixing(lj_mixing_geometric, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%potential%generation() > generation, more="a mixing change kept the generation")
      if (allocated(error)) return
      generation = solvent%potential%generation()
      call solvent%potential%set_ng_split(.true., 0.5_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%potential%generation() > generation, more="an Ng split change kept the generation")
      if (allocated(error)) return
      generation = solvent%potential%generation()
      call solvent%potential%set_mixing(0, err)
      call check(error, allocated(err), more="an unknown mixing rule was accepted")
      if (allocated(error)) return
      deallocate (err)
      call check(error, solvent%potential%generation(), generation, more="a refused setting changed the generation")
      if (allocated(error)) return
      call solvent%update(err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%potential%generation() > generation, more="an update kept the generation")
      if (allocated(error)) return
      call check(error, .not. solvent%is_solved(), more="an unsolved solvent reports a solution")
      if (allocated(error)) return
      allocate (copy, source=solvent)
      call check(error, copy%potential%generation(), solvent%potential%generation(), &
         & more="the copy lost the generation")
   end subroutine check_stale_solution

end module test_moz_solvent
