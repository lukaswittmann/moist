!> Tests of the MOZ electrostatic terms: moz_potential_electrostatic
!>
!> - Fixed charges and host-fed monopoles
!> - Coupling declarations per phase and tail charges from scoped views
!> - Host point-charge field against fixed-charge tables, volume and radial
!> - Adjoint response items for host, with central-difference checks
module test_moz_potential_electrostatic
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_channels_coupling, only: coupling_type, coupling_view_type, coupling_request_type, &
      & coupling_begin_registration, coupling_set_scope, coupling_snapshot, coupling_arm, &
      & coupling_make_view, coupling_close_view, moist_phase_energy, moist_phase_response, moist_phase_gradient
   use moist_channels_response, only: response_type, response_item_type, potential_adjoint_response_type, &
      & atomic_charge_adjoint_response_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_model_moz_potential_term, only: potential_term_type
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_potential_sites, only: potential_sites_type
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use moist_model_moz_potential_electrostatic_fixed, only: fixed_charges_type, new_fixed_charges
   use moist_model_moz_potential_electrostatic_monopole, only: monopole_type
   use moist_model_moz_potential_electrostatic_ec, only: ec_charges_type
   use moist_model_moz_potential_electrostatic_multipole, only: multipole_type
   use moist_model_moz_potential_electrostatic_solvent, only: solvent_charges_type, solvent_multipoles_type, &
      & solvent_model_charges_type, coulomb_resp_solvent, coulomb_spce, multipoles_mbis_gas
   use moist_data_solvents, only: get_solvent_charges
   use test_moz_fixtures, only: water_id, new_water_system
   use moist_data_solvents, only: solvation_system_type, new_solvation_system, get_solvent_id
   use test_helpers, only: check_moist_error
   use test_moz_potential_kernel, only: make_potential, make_solute, make_solvent, make_grid, &
      & check_close, sigma_u, eps_u, q_u, sigma_v, eps_v, q_v, weights_3d, eval_functional, lb, weights_1d, &
      & tables_3d, adjoint_3d
   implicit none(type, external)
   private
   public :: collect_moz_potential_electrostatic

   !> Ng split parameter, 1/bohr
   real(wp), parameter :: alpha = 1.0_wp

   !> Declare a potential for a radial or a volume grid
   interface declare_potential
      module procedure declare_potential_1d
      module procedure declare_potential_3d
   end interface declare_potential

contains

   !> Register tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_potential_electrostatic(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("fixed_charges", check_fixed_charges), &
         & new_unittest("solvent_charges_from_table", check_solvent_charges), &
         & new_unittest("solvent_model_charges", check_solvent_model_charges), &
         & new_unittest("charge_fallback", check_charge_fallback), &
         & new_unittest("monopole_requests_by_phase", check_monopole_requests), &
         & new_unittest("monopole_tail_charges_from_view", check_monopole_tail_charges), &
         & new_unittest("pending_terms", check_pending_terms), &
         & new_unittest("monopole_tables_match_fixed", check_monopole_tables), &
         & new_unittest("monopole_adjoint_gradient_phase", check_monopole_adjoint), &
         & new_unittest("monopole_adjoint_response_phase", check_monopole_response_phase), &
         & new_unittest("monopole_1d", check_monopole_1d)]
   end subroutine collect_moz_potential_electrostatic

   !* ================================================================================= *!
   !*                                   Shared helpers                                  *!
   !* ================================================================================= *!

   !> Update a solute potential of test LJ plus one electrostatic term
   !>
   !> @param[out] pot   updated potential
   !> @param[in]  mol   solute structure
   !> @param[in]  term  electrostatic term
   !> @param[out] err   build failure
   subroutine make_monopole_potential(pot, mol, term, err)
      !> Updated potential
      type(moz_potential_type), intent(out) :: pot
      !> Solute structure
      type(structure_type), intent(in) :: mol
      !> Electrostatic term
      class(potential_term_type), intent(in) :: term
      !> Build failure
      type(moist_error), allocatable, intent(out) :: err
      type(custom_lj_type) :: lj
      call new_custom_lj(lj, sigma_u, eps_u, err)
      if (.not. allocated(err)) call pot%add(lj, err)
      if (.not. allocated(err)) call pot%add(term, err)
      if (.not. allocated(err)) call pot%update(mol, err)
   end subroutine make_monopole_potential

   !> Declare every term of a potential for a radial grid, term i in scope i, and record the extents
   !>
   !> @param[in,out] coupling  coupling
   !> @param[in,out] pot       updated potential
   !> @param[in]     grid      radial grid
   !> @param[out]    err       declaration failure
   subroutine declare_potential_1d(coupling, pot, grid, err)
      !> Coupling
      type(coupling_type), intent(inout) :: coupling
      !> Updated potential
      type(moz_potential_type), intent(inout) :: pot
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Declaration failure
      type(moist_error), allocatable, intent(out) :: err
      call coupling_begin_registration(coupling)
      call pot%declare(grid, coupling, err)
      if (allocated(err)) return
      call coupling_snapshot(coupling, ngrid=grid%npts, natom=pot%natom())
   end subroutine declare_potential_1d

   !> Declare every term of a potential for a volume grid, term i in scope i, and record the extents
   !>
   !> @param[in,out] coupling  coupling
   !> @param[in,out] pot       updated potential
   !> @param[in]     grid      volume grid
   !> @param[out]    err       declaration failure
   subroutine declare_potential_3d(coupling, pot, grid, err)
      !> Coupling
      type(coupling_type), intent(inout) :: coupling
      !> Updated potential
      type(moz_potential_type), intent(inout) :: pot
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Declaration failure
      type(moist_error), allocatable, intent(out) :: err
      call coupling_begin_registration(coupling)
      call pot%declare(grid, coupling, err)
      if (allocated(err)) return
      call coupling_snapshot(coupling, ngrid=grid%ngrid, natom=pot%natom())
   end subroutine declare_potential_3d

   !> Stage a phase and answer pending charge requests
   !>
   !> @param[in,out] coupling  declared coupling
   !> @param[in]     phase     phase to stage
   !> @param[in]     q         atomic charges, e (natom)
   !> @param[out]    err       staging or answer failure
   subroutine stage_charges(coupling, phase, q, err)
      !> Declared coupling
      type(coupling_type), intent(inout) :: coupling
      !> Phase to stage
      integer, intent(in) :: phase
      !> Atomic charges, e (natom)
      real(wp), intent(in) :: q(:)
      !> Staging or answer failure
      type(moist_error), allocatable, intent(out) :: err
      class(coupling_request_type), allocatable :: item
      call coupling_arm(coupling, phase, err)
      if (allocated(err)) return
      do while (coupling%next())
         item = coupling%request()
         if (item%is_missing("q")) call coupling%answer("q", q, err)
         if (allocated(err)) return
      end do
   end subroutine stage_charges

   !> Collect one unanswered coupling walk
   !>
   !> - Visited requests and declared outputs; missing output marked `*`
   !> - Subroutine output to avoid deferred-length function-result races between tests
   !>
   !> @param[in,out] coupling  staged coupling
   !> @param[out]    summary   visited requests and their outputs
   subroutine walk_summary(coupling, summary)
      !> Staged coupling
      type(coupling_type), intent(inout) :: coupling
      !> Visited requests and their outputs
      character(len=:), allocatable, intent(out) :: summary
      !> Outputs a host term may declare, in listing order
      character(len=8), parameter :: outputs(3) = [character(len=8) :: "phi", "dphi_dr", "q"]
      class(coupling_request_type), allocatable :: item
      character(len=:), allocatable :: listed
      integer :: i
      summary = ""
      do while (coupling%next())
         item = coupling%request()
         listed = ""
         do i = 1, size(outputs)
            if (item%output_extent(trim(outputs(i))) == 0) cycle
            if (len(listed) > 0) listed = listed//","
            listed = listed//trim(outputs(i))
            if (item%is_missing(trim(outputs(i)))) listed = listed//"*"
         end do
         summary = summary//trim(item%name())//"("//listed//");"
      end do
   end subroutine walk_summary

   !> Compare the walk of each phase of a freshly declared coupling
   !>
   !> @param[out]    error     test error
   !> @param[in,out] coupling  declared coupling
   !> @param[in]     walks     expected walks of the energy, response and gradient phases
   subroutine check_walks(error, coupling, walks)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Declared coupling
      type(coupling_type), intent(inout) :: coupling
      !> Expected walks of the energy, response and gradient phases
      character(len=*), intent(in) :: walks(:)
      type(moist_error), allocatable :: err
      character(len=:), allocatable :: summary
      integer, parameter :: phases(3) = [moist_phase_energy, moist_phase_response, moist_phase_gradient]
      integer :: ip
      do ip = 1, size(phases)
         call coupling_arm(coupling, phases(ip), err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call walk_summary(coupling, summary)
         call check(error, summary, trim(walks(ip)))
         if (allocated(error)) return
      end do
   end subroutine check_walks

   !> Extract the potential and charge adjoint items of a response
   !>
   !> @param[in,out] response  response
   !> @param[out]    w_phi     potential adjoint, unallocated when absent
   !> @param[out]    dg_dq     charge adjoint, unallocated when absent
   subroutine read_items(response, w_phi, dg_dq)
      !> Response
      type(response_type), intent(inout) :: response
      !> Potential adjoint
      real(wp), allocatable, intent(out) :: w_phi(:)
      !> Charge adjoint
      real(wp), allocatable, intent(out) :: dg_dq(:)
      class(response_item_type), allocatable :: item
      do while (response%next())
         item = response%item()
         select type (item)
         type is (potential_adjoint_response_type)
            w_phi = item%w_phi
         type is (atomic_charge_adjoint_response_type)
            dg_dq = item%dg_dq
         end select
      end do
   end subroutine read_items

   !* ================================================================================= *!
   !*                                       Tests                                       *!
   !* ================================================================================= *!

   !> Check fixed-charge validation and tail charges
   subroutine check_fixed_charges(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(fixed_charges_type) :: term
      type(structure_type) :: mol
      type(coupling_view_type) :: unbound
      real(wp), allocatable :: tail(:)
      real(wp) :: q(0)

      call make_solute(mol)
      call new_fixed_charges(term, q, err)
      call check(error, allocated(err), more="empty charges were accepted")
      if (allocated(error)) return
      call new_fixed_charges(term, [0.1_wp, ieee_value(1.0_wp, ieee_positive_inf)], err)
      call check(error, allocated(err), more="non-finite charges were accepted")
      if (allocated(error)) return
      call term%update(mol, err)
      call check(error, allocated(err), more="an unconstructed term was built")
      if (allocated(error)) return

      call new_fixed_charges(term, [0.1_wp, -0.1_wp], err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call term%update(mol, err)
      call check(error, allocated(err), more="charges of the wrong length were built")
      if (allocated(error)) return

      call new_fixed_charges(term, q_u, err)
      if (.not. allocated(err)) call term%update(mol, err)
      if (.not. allocated(err)) call term%tail_charges(unbound, tail, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. allocated(term%pair) .and. .not. term%is_host_fed() .and. .not. term%has_field(), &
         & more="fixed charges supply neither LJ data nor a field and read no host data")
      if (allocated(error)) return
      call check(error, allocated(tail), more="fixed charges supplied no tail charges")
      if (allocated(error)) return
      call check(error, all(tail == q_u), more="tail charges are the fixed charges")
   end subroutine check_fixed_charges

   !> Check table charges by solvent id, environment and model
   !>
   !> - No table solvent: no charges
   !> - Multipole stub: refusal
   subroutine check_solvent_charges(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(solvent_charges_type) :: term, unset
      type(solvent_multipoles_type) :: multipoles
      type(structure_type) :: mol, fragment
      type(solvation_system_type) :: system
      type(coupling_view_type) :: unbound
      real(wp), allocatable :: tail(:), reference(:)

      call new_water_system(system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = system%solv_mol
      call make_solvent(fragment)
      term = solvent_charges_type(environment="vacuum", model="resp")
      call term%update(mol, err, water_id)
      call check(error, allocated(err), more="an unknown environment was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Unknown solvent environment 'vacuum'") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)
      term = solvent_charges_type(environment="gas", model="mulliken")
      call term%update(mol, err, water_id)
      call check(error, allocated(err), more="an unknown charge model was accepted")
      if (allocated(error)) return
      deallocate (err)
      call unset%update(mol, err, water_id)
      call check(error, allocated(err), more="a term without a choice was built")
      if (allocated(error)) return
      deallocate (err)

      term = solvent_charges_type(environment="Solvent", model="RESP")
      ! Solute or custom solvent: no table charges, every atom uncovered
      call term%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. any(term%charges_covered(mol%nat)), more="a side without a table solvent was covered")
      if (allocated(error)) return
      call term%update(mol, err, 0)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. any(term%charges_covered(mol%nat)), more="a custom solvent was covered")
      if (allocated(error)) return
      call term%update(fragment, err, water_id)
      call check(error, allocated(err), more="a structure of another size than the table entry was built")
      if (allocated(error)) return
      call term%update(mol, err, water_id)
      if (.not. allocated(err)) call term%tail_charges(unbound, tail, err)
      if (.not. allocated(err)) call get_solvent_charges(water_id, "solvent", "resp", reference, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. term%is_host_fed() .and. .not. term%has_field(), more="solvent charges are parameter-built")
      if (allocated(error)) return
      call check(error, all(term%charges_covered(mol%nat)), more="water table charges must cover every atom")
      if (allocated(error)) return
      call check(error, size(tail), mol%nat, more="tail charge count")
      if (allocated(error)) return
      call check(error, all(tail == reference), more="tail charges differ from the table")
      if (allocated(error)) return
      call check(error, abs(sum(tail)) < 1.0e-5_wp, more="water charges do not sum to zero")
      if (allocated(error)) return

      multipoles = multipoles_mbis_gas
      call multipoles%update(mol, err, water_id)
      call check(error, allocated(err), more="the multipole stub built")
      if (allocated(error)) return
      call check(error, index(err%message, "not implemented") > 0, more=err%message)
   end subroutine check_solvent_charges

   !> Check neutral-water SPC/E coverage and RESP preset against table
   subroutine check_solvent_model_charges(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(solvent_model_charges_type) :: spce, unknown
      type(solvent_charges_type) :: resp
      type(solvation_system_type) :: system
      type(structure_type) :: mol, fragment
      type(coupling_view_type) :: unbound
      real(wp), allocatable :: tail(:), reference(:)

      call new_water_system(system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = system%solv_mol
      spce = coulomb_spce
      call spce%update(mol, err, water_id)
      if (.not. allocated(err)) call spce%tail_charges(unbound, tail, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(tail == merge(-0.8476_wp, 0.4238_wp, mol%num(mol%id) == 8)), &
         & more="SPC/E charges differ")
      if (allocated(error)) return
      call check(error, abs(sum(tail)) < 1.0e-12_wp, more="SPC/E water is not neutral")
      if (allocated(error)) return

      call check(error, all(spce%charges_covered(mol%nat)), more="SPC/E must cover water")
      if (allocated(error)) return

      ! Non-neutral or non-water structure: atoms left for later charge terms
      call make_solvent(fragment)
      call spce%update(fragment, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. any(spce%charges_covered(fragment%nat)), more="SPC/E covered a non-water structure")
      if (allocated(error)) return
      unknown = solvent_model_charges_type(set=7)
      call unknown%update(mol, err, water_id)
      call check(error, allocated(err), more="an unknown solvent charge model was built")
      if (allocated(error)) return
      deallocate (err)

      resp = coulomb_resp_solvent
      call resp%update(mol, err, water_id)
      if (.not. allocated(err)) call resp%tail_charges(unbound, tail, err)
      if (.not. allocated(err)) call get_solvent_charges(water_id, "solvent", "resp", reference, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(tail == reference), more="the RESP preset differs from the table")
   end subroutine check_solvent_model_charges

   !> Check ordered charge merging
   !>
   !> - First charge supplier per atom; later terms fill uncovered atoms
   !> - Atom left without charge: build failure
   subroutine check_charge_fallback(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(solvation_system_type) :: system
      type(moz_potential_type) :: pot, partial
      type(potential_sites_type) :: sites
      type(fixed_charges_type) :: first, rest
      type(structure_type) :: mol
      real(wp), allocatable :: q(:), reference(:)
      integer :: id

      ! Methanol: no SPC/E charges, full RESP table coverage
      call get_solvent_id("methanol", id, err)
      if (.not. allocated(err)) call new_solvation_system(system, id, error=err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = system%solv_mol
      call pot%add(coulomb_spce, err)
      if (.not. allocated(err)) call pot%add(coulomb_resp_solvent, err)
      if (.not. allocated(err)) call pot%update(mol, err, id)
      if (.not. allocated(err)) call pot%sites(sites, err)
      if (.not. allocated(err)) call get_solvent_charges(id, "solvent", "resp", reference, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      q = sites%q
      call check(error, all(q == reference), more="merged charges are not the RESP charges")
      if (allocated(error)) return

      ! First-term ownership retained; second term fills remaining atoms, no charge sum
      call new_fixed_charges(first, [0.5_wp, (0.0_wp, id=2, mol%nat)], err, [.true., (.false., id=2, mol%nat)])
      if (.not. allocated(err)) call new_fixed_charges(rest, [(0.1_wp, id=1, mol%nat)], err)
      if (.not. allocated(err)) call partial%add(first, err)
      if (.not. allocated(err)) call partial%add(rest, err)
      if (.not. allocated(err)) call partial%update(mol, err)
      if (.not. allocated(err)) call partial%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      q = sites%q
      call check(error, q(1) == 0.5_wp .and. all(q(2:) == 0.1_wp), more="charges were summed or not filled")
      if (allocated(error)) return

      ! Uncovered atoms without fallback: named build error
      partial = moz_potential_type()
      call partial%add(first, err)
      if (.not. allocated(err)) call partial%update(mol, err)
      call check(error, allocated(err), more="atoms without charges were accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "No charges for atoms: O2") > 0, more=err%message)
   end subroutine check_charge_fallback

   !> Check monopole declarations per phase and solvent-side refusal
   !>
   !> - Own charge request on radial and volume grids
   !>
   !> @param[out] error  test error
   subroutine check_monopole_requests(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moz_potential_type) :: pot, solvent_pot
      type(monopole_type) :: charges
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(coupling_type), target :: coupling_3d, coupling_1d
      character(len=96) :: walks(3)

      call make_solute(mol)
      call make_grid(grid, mol, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 8, 0.6_wp, err)
      if (.not. allocated(err)) call make_monopole_potential(pot, mol, charges, err)
      if (.not. allocated(err)) call declare_potential(coupling_3d, pot, grid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      walks = "atomic_charges(q*);"
      call check_walks(error, coupling_3d, walks)
      if (allocated(error)) return
      call declare_potential(coupling_1d, pot, rgrid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_walks(error, coupling_1d, walks)
      if (allocated(error)) return

      call solvent_pot%refuse_host_fed_terms()
      call solvent_pot%add(charges, err)
      call check(error, allocated(err), more="a solvent potential accepted host-fed monopoles")
      if (allocated(error)) return
      call check(error, index(err%message, "cannot serve a solvent") > 0, more=err%message)
   end subroutine check_monopole_requests

   !> Check host tail charges from scoped coupling view
   !>
   !> - Refusal before build or through unbound view
   !> - Host answer after staging
   subroutine check_monopole_tail_charges(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol_u
      type(moz_potential_type) :: pot_b
      type(potential_sites_type) :: sites
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(monopole_type) :: charges
      type(coupling_type), target :: coupling
      type(coupling_view_type) :: view, unbound
      real(wp), allocatable :: q(:)

      call make_solute(mol_u)
      call charges%tail_charges(unbound, q, err)
      call check(error, allocated(err) .and. .not. allocated(q), more="an unbuilt host term supplied charges")
      if (allocated(error)) return
      call check(error, index(err%message, "is not built") > 0, more=err%message)
      if (allocated(error)) return

      call make_grid(grid, mol_u, err)
      if (.not. allocated(err)) call make_monopole_potential(pot_b, mol_u, charges, err)
      if (.not. allocated(err)) call charges%update(mol_u, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call charges%tail_charges(unbound, q, err)
      call check(error, allocated(err) .and. .not. allocated(q), more="host charges read through an unbound view")
      if (allocated(error)) return
      deallocate (err)
      call pot_b%sites(sites, err)
      call check(error, allocated(err), more="the potential resolved host charges without the coupling")
      if (allocated(error)) return
      call check(error, index(err%message, "need the coupling") > 0, more=err%message)
      if (allocated(error)) return

      call declare_potential(coupling, pot_b, grid, err)
      if (.not. allocated(err)) call stage_charges(coupling, moist_phase_energy, q_u, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call coupling_make_view(coupling, 2, view)
      call charges%tail_charges(view, q, err)
      call coupling_close_view(coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_close(error, q, q_u, 0.0_wp, "host tail charges")
   end subroutine check_monopole_tail_charges

   !> Check EC construction and multipole stub refusal
   subroutine check_pending_terms(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(ec_charges_type) :: ec
      type(multipole_type) :: multipoles

      call make_solute(mol)
      call ec%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call multipoles%update(mol, err)
      call check(error, allocated(err), more="host multipoles built")
      if (allocated(error)) return
      call check(error, index(err%message, "not implemented in this build round") > 0, more=err%message)
   end subroutine check_pending_terms

   !> Compare host-fed and fixed-charge tables
   !>
   !> - Monopole tables from the same point charges
   subroutine check_monopole_tables(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_a, pot_b, pot_v
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(monopole_type) :: charges
      type(coupling_type), target :: coupling, unanswered, unused
      real(wp), allocatable :: sr_a(:, :), lr_a(:, :), sr_b(:, :), lr_b(:, :)
      complex(wp), allocatable :: uk_a(:, :), uk_b(:, :)
      logical :: split
      integer :: isplit

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_grid(grid, mol_u, err)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      if (.not. allocated(err)) call make_potential(pot_a, mol_u, err, sigma_u, eps_u, q_u)
      call check_moist_error(error, err)
      if (allocated(error)) return

      call make_monopole_potential(pot_b, mol_u, charges, err)
      if (.not. allocated(err)) call declare_potential(coupling, pot_b, grid, err)
      if (.not. allocated(err)) call stage_charges(coupling, moist_phase_energy, q_u, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      ! Missing host answers: refusal
      call tables_3d(grid, pot_b, pot_v, lb, .true., sr_b, lr_b, uk_b, err, coupling=unanswered)
      call check(error, allocated(err), more="a host-fed potential was updated without the host's answers")
      if (allocated(error)) return

      do isplit = 1, 2
         split = isplit == 1
         call tables_3d(grid, pot_a, pot_v, lb, split, sr_a, lr_a, uk_a, err, coupling=unused)
         if (.not. allocated(err)) call tables_3d(grid, pot_b, pot_v, lb, split, sr_b, lr_b, uk_b, err, coupling=coupling)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check_close(error, reshape(sr_b, [size(sr_b)]), reshape(sr_a, [size(sr_a)]), 1.0e-12_wp, "host u_sr")
         if (allocated(error)) return
         call check_close(error, reshape(lr_b, [size(lr_b)]), reshape(lr_a, [size(lr_a)]), 1.0e-12_wp, "host ur_lr")
         if (allocated(error)) return
         call check(error, maxval(abs(uk_b - uk_a)) < 1.0e-12_wp, more="host uk_lr")
         if (allocated(error)) return
      end do
   end subroutine check_monopole_tables

   !> Compare gradient-phase monopoles with fixed charges
   !>
   !> - Same nuclear and grid gradients as fixed charges
   !> - Charge adjoint against central differences
   !>
   !> @param[out] error  test error
   subroutine check_monopole_adjoint(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_a, pot_c, pot_v, pot_q
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(monopole_type) :: charges
      type(coupling_type), target :: coupling_a, coupling_c
      type(volume_adjoint_type) :: acc_a, acc_c
      type(response_type) :: resp_a, resp_c
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), g_a(:, :), g_c(:, :)
      real(wp), allocatable :: w_phi_a(:), dq_a(:), w_phi_c(:), dq_c(:), q(:)
      complex(wp), allocatable :: w_k(:, :)
      real(wp), parameter :: hq = 1.0e-3_wp
      real(wp) :: lp, lm, fd
      integer :: a

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_grid(grid, mol_u, err)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      if (.not. allocated(err)) call make_potential(pot_a, mol_u, err, sigma_u, eps_u, q_u)
      if (.not. allocated(err)) call make_monopole_potential(pot_c, mol_u, charges, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call declare_potential(coupling_c, pot_c, grid, err)
      if (.not. allocated(err)) call stage_charges(coupling_c, moist_phase_gradient, q_u, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      call weights_3d(grid%ngrid, grid%npts_k, 2, w_sr, w_lr, w_k)
      allocate (g_a(3, 3), g_c(3, 3), source=0.0_wp)
      call acc_a%init(grid%ngrid)
      call acc_c%init(grid%ngrid)
      call adjoint_3d(grid, pot_a, pot_v, lb, .true., coupling_a, w_sr, w_lr, w_k, resp_a, err, gradient=g_a, grid_adjoint=acc_a)
      if (.not. allocated(err)) call adjoint_3d(grid, pot_c, pot_v, lb, .true., coupling_c, w_sr, w_lr, w_k, resp_c, err, &
         & gradient=g_c, grid_adjoint=acc_c)
      call check_moist_error(error, err)
      if (allocated(error)) return

      call read_items(resp_a, w_phi_a, dq_a)
      call check(error, .not. allocated(w_phi_a) .and. .not. allocated(dq_a), &
         & more="fixed charges pushed response items")
      if (allocated(error)) return
      call read_items(resp_c, w_phi_c, dq_c)
      call check(error, .not. allocated(w_phi_c) .and. allocated(dq_c), &
         & more="monopoles must push only the charge adjoint")
      if (allocated(error)) return
      call check_close(error, reshape(g_c, [9]), reshape(g_a, [9]), 1.0e-12_wp, "monopole gradient")
      if (allocated(error)) return
      call check_close(error, reshape(acc_c%w_xyz, [size(acc_c%w_xyz)]), reshape(acc_a%w_xyz, [size(acc_a%w_xyz)]), &
         & 1.0e-12_wp, "monopole grid adjoint")
      if (allocated(error)) return

      do a = 1, 3
         q = q_u
         q(a) = q(a) + hq
         call make_potential(pot_q, mol_u, err, sigma_u, eps_u, q)
         if (.not. allocated(err)) call eval_functional(grid, mol_u, pot_q, pot_v, w_sr, w_lr, w_k, lp, err)
         q(a) = q(a) - 2.0_wp*hq
         if (.not. allocated(err)) call make_potential(pot_q, mol_u, err, sigma_u, eps_u, q)
         if (.not. allocated(err)) call eval_functional(grid, mol_u, pot_q, pot_v, w_sr, w_lr, w_k, lm, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         fd = (lp - lm)/(2.0_wp*hq)
         call check_close(error, [dq_c(a)], [fd], 1.0e-9_wp, "monopole charge adjoint")
         if (allocated(error)) return
      end do
   end subroutine check_monopole_adjoint

   !> Check response-phase adjoints without geometry derivatives
   !>
   !> - Charge adjoint without nuclear or grid gradients
   !> - Response items equal to gradient-phase items
   subroutine check_monopole_response_phase(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_b, pot_v
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(monopole_type) :: charges
      type(coupling_type), target :: coupling, coupling_g
      type(volume_adjoint_type) :: acc
      type(response_type) :: resp, resp_g
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), grad(:, :), w_phi(:), dq(:), w_phi_g(:), dq_g(:)
      complex(wp), allocatable :: w_k(:, :)

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_grid(grid, mol_u, err)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      if (.not. allocated(err)) call make_monopole_potential(pot_b, mol_u, charges, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call declare_potential(coupling, pot_b, grid, err)
      if (.not. allocated(err)) call stage_charges(coupling, moist_phase_response, q_u, err)
      if (.not. allocated(err)) call declare_potential(coupling_g, pot_b, grid, err)
      if (.not. allocated(err)) call stage_charges(coupling_g, moist_phase_gradient, q_u, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      call weights_3d(grid%ngrid, grid%npts_k, 2, w_sr, w_lr, w_k)
      call adjoint_3d(grid, pot_b, pot_v, lb, .true., coupling, w_sr, w_lr, w_k, resp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (grad(3, 3), source=0.0_wp)
      call acc%init(grid%ngrid)
      call adjoint_3d(grid, pot_b, pot_v, lb, .true., coupling_g, w_sr, w_lr, w_k, resp_g, err, gradient=grad, grid_adjoint=acc)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call read_items(resp, w_phi, dq)
      call read_items(resp_g, w_phi_g, dq_g)
      call check(error, .not. allocated(w_phi) .and. .not. allocated(w_phi_g) .and. allocated(dq), &
         & more="monopoles must push only the charge adjoint in every phase")
      if (allocated(error)) return
      call check_close(error, dq, dq_g, 1.0e-14_wp, "response-phase charge adjoint")
   end subroutine check_monopole_response_phase

   !> Check 1D host-charge tables and adjoints
   !>
   !> - Charges through coupling; tables match fixed charges, refusal without coupling
   !> - Charge adjoint against central differences of fixed-charge functional
   !> - Fixed-charge adjoint: no response items
   subroutine check_monopole_1d(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_a, pot_c, pot_v, pot_q
      type(monopole_type) :: charges
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(coupling_type), target :: coupling_a, coupling_c
      type(response_type) :: resp_a, resp_c
      real(wp), allocatable :: sr_a(:, :), lr_a(:, :), uk_a(:, :), sr_c(:, :), lr_c(:, :), uk_c(:, :)
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), w_k(:, :), q(:)
      real(wp), allocatable :: w_phi_a(:), dq_a(:), w_phi_c(:), dq_c(:)
      !> Site pairs, solute major
      integer, parameter :: pairs(2, 6) = reshape([1, 1, 1, 2, 2, 1, 2, 2, 3, 1, 3, 2], [2, 6])
      !> - L linear in charges; exact central differences for any step
      !> - Step large enough to exceed LJ-wall roundoff at first node
      real(wp), parameter :: hq = 0.1_wp
      real(wp) :: lp, lm, fd
      integer :: a

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      if (.not. allocated(err)) call make_potential(pot_a, mol_u, err, sigma_u, eps_u, q_u)
      if (.not. allocated(err)) call make_monopole_potential(pot_c, mol_u, charges, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 32, 0.6_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call declare_potential(coupling_c, pot_c, rgrid, err)
      if (.not. allocated(err)) call stage_charges(coupling_c, moist_phase_response, q_u, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      call pot_c%compute(rgrid, kgrid, sr_c, lr_c, uk_c, err, solvent=pot_v)
      call check(error, allocated(err), more="host charges were assembled in 1D without the coupling")
      if (allocated(error)) return
      call check(error, index(err%message, "need the coupling") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)
      call pot_a%compute(rgrid, kgrid, sr_a, lr_a, uk_a, err, solvent=pot_v)
      if (.not. allocated(err)) call pot_c%compute(rgrid, kgrid, sr_c, lr_c, uk_c, err, solvent=pot_v, coupling=coupling_c)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(sr_c == sr_a) .and. all(lr_c == lr_a) .and. all(uk_c == uk_a), &
         & more="1D tables of host charges differ from fixed charges")
      if (allocated(error)) return

      call weights_1d(rgrid%npts, kgrid%npts, size(pairs, 2), w_sr, w_lr, w_k)
      call pot_a%compute_adjoint(rgrid, kgrid, pot_v, coupling_a, w_sr, w_lr, w_k, resp_a, err)
      if (.not. allocated(err)) call pot_c%compute_adjoint(rgrid, kgrid, pot_v, coupling_c, w_sr, w_lr, w_k, resp_c, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call read_items(resp_a, w_phi_a, dq_a)
      call check(error, .not. allocated(w_phi_a) .and. .not. allocated(dq_a), &
         & more="fixed charges pushed response items in 1D")
      if (allocated(error)) return
      call read_items(resp_c, w_phi_c, dq_c)
      call check(error, .not. allocated(w_phi_c) .and. allocated(dq_c), &
         & more="1D host charges push only the charge adjoint")
      if (allocated(error)) return
      call check(error, size(dq_c) == 3, more="1D host charge adjoint is not per solute atom")
      if (allocated(error)) return

      ! Charge adjoint against central differences of fixed charges (L is linear in q)
      do a = 1, 3
         q = q_u
         q(a) = q(a) + hq
         call make_potential(pot_q, mol_u, err, sigma_u, eps_u, q)
         if (.not. allocated(err)) call pot_q%compute(rgrid, kgrid, sr_a, lr_a, uk_a, err, solvent=pot_v)
         if (allocated(err)) exit
         lp = sum(w_sr*sr_a) + sum(w_lr*lr_a) + sum(w_k*uk_a)
         q(a) = q(a) - 2.0_wp*hq
         call make_potential(pot_q, mol_u, err, sigma_u, eps_u, q)
         if (.not. allocated(err)) call pot_q%compute(rgrid, kgrid, sr_a, lr_a, uk_a, err, solvent=pot_v)
         if (allocated(err)) exit
         lm = sum(w_sr*sr_a) + sum(w_lr*lr_a) + sum(w_k*uk_a)
         fd = (lp - lm)/(2.0_wp*hq)
         call check_close(error, [dq_c(a)], [fd], 1.0e-9_wp, "1D host charge adjoint")
         if (allocated(error)) return
      end do
      call check_moist_error(error, err)
   end subroutine check_monopole_1d

end module test_moz_potential_electrostatic
