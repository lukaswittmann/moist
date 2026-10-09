!> Tests of the MOZ bulk solvent: moz_solvent
module test_moz_solvent
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_data_solvents, only: solvation_system_type, new_solvation_system
   use moist_math_packed_sym, only: packed_index, npair_from_ns
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_potential_lj_base, only: lj_mixing_lorentz_berthelot, lj_mixing_geometric
   use moist_model_moz_potential_electrostatic_host_charges, only: host_charges_type
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use test_helpers, only: check_moist_error
   use test_moz_fixtures, only: water_id, water_sigma, water_epsilon, water_charges, new_water_system, &
      & new_water_solvent
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
         & new_unittest("construct_from_table", check_construct_from_table), &
         & new_unittest("construct_custom_system", check_construct_custom), &
         & new_unittest("density_derivation", check_density_derivation), &
         & new_unittest("host_fed_refused", check_host_fed_refused), &
         & new_unittest("empty_potential_update", check_empty_potential_update), &
         & new_unittest("packed_site_pairs", check_site_pairs), &
         & new_unittest("vv_potential_tables", check_potential_tables), &
         & new_unittest("vv_tables_geometric_mixing", check_geometric_tables), &
         & new_unittest("solve_pending", check_solve_pending)]
   end subroutine collect_moz_solvent

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
      type(solvation_system_type) :: system
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

      call new_solvation_system(system, water_id, temperature=350.0_wp, error=err)
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
      type(solvation_system_type) :: system
      type(structure_type) :: mol, copy

      call new_context(ctx, nthreads=0, verbosity=0)
      call new(mol, [8, 1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp, 0.0_wp, 0.0_wp, &
         & -0.5_wp, 1.7_wp, 0.0_wp], [3, 3]))
      system%solvent_id = 0
      system%temperature = 298.15_wp
      system%solvent_number_density_au = 5.0e-3_wp
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

      system%solvent_number_density_au = 0.0_wp
      call new_vv_solvent(solvent, ctx, system, err)
      call check(error, allocated(err), more="a zero density was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "density must be positive") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)
      system%solvent_number_density_au = 5.0e-3_wp
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
      type(solvation_system_type) :: system
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
      call check(error, system%solvent_number_density_au, water_per_bohr3, thr=0.01_wp, rel=.true., &
         & more="water number density is off by more than 1%")
   end subroutine check_density_derivation

   !> Check host-term refusal at solvent add, built or unbuilt; unchanged potential
   subroutine check_host_fed_refused(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvation_system_type) :: system
      type(solvent_vv_type) :: solvent
      type(moz_potential_type) :: pot
      type(host_charges_type) :: host
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
      type(solvation_system_type) :: system
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

   !> Check pending VV solve and unsolved state
   subroutine check_solve_pending(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_water_solvent(ctx, solvent, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call solvent%solve(err)
      call check(error, allocated(err), more="the pending VV solve succeeded")
      if (allocated(error)) return
      call check(error, index(err%message, "1D VV solve is pending") > 0, more=err%message)
      if (allocated(error)) return
      call check(error, .not. solvent%is_solved(), more="a failed solve marked the solvent solved")
   end subroutine check_solve_pending

end module test_moz_solvent
