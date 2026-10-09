!> 1D MOZ model contract tests
module test_moz_1d
   use mctc_env, only: wp, moist_error => error_type, fatal_error
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_data_solvents, only: solvation_system_type
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use moist_model_moz_potential_electrostatic_fixed, only: fixed_charges_type, new_fixed_charges
   use moist_model_moz_potential_electrostatic_host_charges, only: host_charges_type
   use moist_model_moz_potential_electrostatic_host_potential, only: host_potential_type
   use moist_model_moz_potential_electrostatic_multipole, only: host_multipoles_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_channels_coupling, only: coupling_type, coupling_request_type
   use moist_channels_fields, only: field_query_type
   use moist_channels_response, only: response_type, atomic_charge_adjoint_response_type, response_accumulate
   use test_helpers, only: check_moist_error
   use test_moz_fixtures, only: water_sigma, water_epsilon, water_charges, new_water_system, &
      & new_water_solvent, refusing_term_type, answer_charges, answer_radial_potential
   implicit none(type, external)
   private
   public :: collect_moz_1d

contains
   !> Register 1D MOZ tests
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_moz_1d(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("contract", check_moz_1d), &
         new_unittest("getter_guards", check_getter_guards), &
         new_unittest("update_guards", check_update_guards), &
         new_unittest("atom_count_follows_update", check_atom_count_follows_update), &
         new_unittest("host_charges_in_every_phase", check_host_charges_phases), &
         new_unittest("host_potential_declared", check_host_potential), &
         new_unittest("term_failures", check_term_failures), &
         new_unittest("solvent_copy", check_solvent_copy), &
         new_unittest("uv_potential_tables", check_potential_tables)]
   end subroutine collect_moz_1d

   !> Drive a 1D MOZ model through the coupling protocol
   !>
   !> - Host charges: new coupling and atomic_charges request through coupling protocol
   !> - Theory-pending error after all mandatory answers
   !> - Wrong-shape answer rejected first
   subroutine check_moz_1d(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      real(wp) :: energy, gradient(3, 2)

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_host_charges_model(ctx, model, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call new(mol, [1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, model%atom_count(), 2)
      if (allocated(error)) return

      call model%new_coupling(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, coupling%next(), more="atomic_charges must be pending")
      if (allocated(error)) return

      ! A wrong-shape answer is rejected and leaves the output missing
      call coupling%answer("q", [1.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return

      associate (item => coupling%request())
         call check(error, item%is_missing("q"), more="wrong shape satisfied q")
      end associate
      if (allocated(error)) return

      ! The correctly shaped answer (one charge per atom) completes the walk
      call coupling%answer("q", [0.3_wp, -0.3_wp], err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, .not. coupling%next(), more="atomic_charges is the only request")
      if (allocated(error)) return

      energy = 7.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err) .and. energy == 7.0_wp)
      if (allocated(error)) return
      call check(error, index(err%message, "1D MOZ energy") > 0)
      if (allocated(error)) return

      call model%prepare_gradient(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      gradient = 3.0_wp
      call model%get_gradient(coupling, response, gradient, err)
      call check(error, allocated(err) .and. all(gradient == 3.0_wp))
      if (allocated(error)) return
      call check(error, index(err%message, "1D MOZ gradient") > 0)
      if (allocated(error)) return

      call model%release_coupling(coupling)
      call check(error, .not. associated(coupling), more="release left a coupling pointer")
   end subroutine check_moz_1d

   !> An existing coupling follows the model to a new atom count
   !>
   !> The coupling is minted for two atoms; after an update to three the
   !> next staging takes three charges and refuses two
   subroutine check_atom_count_follows_update(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_host_charges_model(ctx, model, err)
      if (.not. allocated(err)) then
         call new(mol, [1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
         call model%update(mol, err)
      end if
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, coupling%next(), more="atomic_charges must be pending")
      if (allocated(error)) return
      call coupling%answer("q", [0.3_wp, -0.3_wp], err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if

      ! Same model and coupling, one more atom
      call new(mol, [8, 1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp, 0.0_wp, 0.0_wp, &
         & -0.5_wp, 1.7_wp, 0.0_wp], [3, 3]))
      call model%update(mol, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, coupling%next(), more="the update kept the two-atom answer")
      if (allocated(error)) return
      call coupling%answer("q", [0.3_wp, -0.3_wp], err)
      call check(error, allocated(err), more="two charges were accepted for three atoms")
      if (allocated(error)) return
      deallocate (err)
      call coupling%answer("q", [-0.6_wp, 0.3_wp, 0.3_wp], err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, .not. coupling%next(), more="three charges complete the walk")
      if (allocated(error)) return

      call model%release_coupling(coupling)
      call check(error, .not. associated(coupling), more="release left a coupling pointer")
   end subroutine check_atom_count_follows_update

   !> Check host-charge declarations
   !>
   !> - `q` pending in energy, response and gradient phases
   subroutine check_host_charges_phases(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      character(len=96) :: walks(3)

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_updated_model(ctx, model, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      walks = "atomic_charges(q*);"
      call check_phase_walks(error, model, walks)
   end subroutine check_host_charges_phases

   !> Check host-potential declarations on the radial grid of the model
   !>
   !> - Radial potential pending in every phase, beside the host charges
   !> - Radial grid published as named fields
   !> - Answer over the wrong number of radii refused; correct shape accepted
   subroutine check_host_potential(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(host_potential_type) :: potential
      type(custom_lj_type) :: lj
      type(coupling_type), pointer :: coupling
      type(field_query_type) :: query
      real(wp), allocatable :: phi(:, :)
      character(len=96) :: walks(3)
      character(len=5), parameter :: names(3) = [character(len=5) :: "ngrid", "natom", "r"]
      integer :: i

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_custom_lj(lj, [3.0_wp, 3.0_wp], [1.0e-4_wp, 1.0e-4_wp], err)
      if (.not. allocated(err)) call new_model(ctx, model, err)
      if (.not. allocated(err)) call model%potential%add(potential, err)
      if (.not. allocated(err)) call model%potential%add(lj, err)
      if (.not. allocated(err)) call update_h2(model, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      walks = "atomic_charges(q*);radial_potential(phi*);"
      call check_phase_walks(error, model, walks)
      if (allocated(error)) return

      do i = 1, size(names)
         call query%fetch(trim(names(i)))
         call model%list_fields(query)
         call check(error, query%found, more="missing radial grid field "//trim(names(i)))
         if (allocated(error)) return
         select case (i)
         case (1)
            call check(error, query%ivals(1), model%rgrid%npts)
         case (2)
            call check(error, query%ivals(1), 2)
         case default
            call check(error, all(query%rvals == model%rgrid%r), more="published radii differ from the grid")
         end select
         if (allocated(error)) return
      end do

      call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      if (.not. allocated(err)) call answer_charges(coupling, [0.1_wp, -0.1_wp], err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, coupling%next(), more="radial_potential must be pending")
      if (allocated(error)) return
      allocate (phi(model%rgrid%npts + 1, 2), source=0.0_wp)
      call coupling%answer("phi", phi, err)
      call check(error, allocated(err), more="a radial potential over the wrong radii was accepted")
      if (allocated(error)) return
      deallocate (err, phi)
      allocate (phi(model%rgrid%npts, 2), source=0.0_wp)
      call coupling%answer("phi", phi, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call model%release_coupling(coupling)
   end subroutine check_host_potential

   !> Check named solute failures and cleared coupling state
   !>
   !> - Solute potential without terms: refusal at update
   !> - Host multipole stub: refusal at update (build)
   !> - Failed term declaration: coupling refused, previously staged walk ended
   subroutine check_term_failures(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(host_multipoles_type) :: multipoles
      type(host_charges_type) :: charges
      type(refusing_term_type) :: refusing
      !> Declaration switch shared with the refusing term
      logical, target :: refuse
      type(coupling_type), pointer :: coupling

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_model(ctx, model, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call update_h2(model, err)
      call check(error, allocated(err) .and. .not. model%is_updated(), more="an empty solute potential was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "has no terms") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      call model%potential%add(multipoles, err)
      if (.not. allocated(err)) call update_h2(model, err)
      call check(error, allocated(err) .and. .not. model%is_updated(), more="the multipoles stub was built")
      if (allocated(error)) return
      call check(error, index(err%message, "host_multipoles") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      ! Host charges first, then pair-only term with switchable declaration failure
      call new_model(ctx, model, err)
      if (.not. allocated(err)) call model%potential%add(charges, err)
      refuse = .false.
      refusing%refuse => refuse
      if (.not. allocated(err)) call model%potential%add(refusing, err)
      if (.not. allocated(err)) call update_h2(model, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      refuse = .true.
      call model%new_coupling(coupling, err)
      call check(error, allocated(err) .and. .not. associated(coupling), more="a failing declaration was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "refused declaration") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      ! Failed redeclaration: coupling refused, previous staging cleared
      refuse = .false.
      call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      refuse = .true.
      call model%prepare_response(coupling, err)
      call check(error, allocated(err), more="a failing declaration was staged")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="the refused declaration left a staging behind")
      if (allocated(error)) return

      call model%release_coupling(coupling)
      call check(error, .not. associated(coupling), more="release left a coupling pointer")
   end subroutine check_term_failures

   !> Check solvent ownership, model invalidation and split validation
   !>
   !> - Built solvent required; independent copy in model
   !> - Added term: model invalidation
   !> - Ng split validation at setter
   subroutine check_solvent_copy(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(solvation_system_type) :: system
      type(solvent_vv_type) :: solvent
      type(moz_potential_type) :: solvent_pot
      type(fixed_charges_type) :: extra
      type(host_charges_type) :: charges
      type(moist_math_grid_radial_type) :: rgrid, kgrid, empty

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_water_system(system, err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new_radial_pair(rgrid, kgrid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new_moz_1d_model(model, ctx, empty, empty, solvent, err)
      call check(error, allocated(err), more="an unbuilt radial grid was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "built real-space radial grid") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)
      call new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, err)
      call check(error, allocated(err), more="a solvent that is not updated was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Update the 1D VV solvent") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      call new_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, err)
      if (.not. allocated(err)) call new_fixed_charges(extra, [0.1_wp, 0.1_wp, -0.2_wp], err)
      if (.not. allocated(err)) call solvent%potential%add(extra, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model%has_solvent(), more="the model holds no solvent")
      if (allocated(error)) return
      solvent_pot = solvent%potential
      call check(error, solvent_pot%n_terms() == 3 .and. .not. solvent%potential%is_updated(), &
         & more="the extra term did not reach the original solvent")
      if (allocated(error)) return
      call check(error, model%solvent_nsite() == 3, more="the model solvent lost its sites")
      if (allocated(error)) return
      call check(error, model%potential%ng_split() .and. model%potential%alpha() == 1.0_wp, more="Ng split defaults changed")
      if (allocated(error)) return

      call model%potential%add(charges, err)
      if (.not. allocated(err)) call update_h2(model, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model%is_updated())
      if (allocated(error)) return
      call model%potential%set_ng_split(.true., -1.0_wp, err)
      call check(error, allocated(err) .and. model%is_updated(), more="a negative alpha was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "alpha must be positive") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)
      call check(error, model%potential%alpha() == 1.0_wp, more="a refused alpha replaced the split")
      if (allocated(error)) return
      call model%potential%set_ng_split(.false., 0.0_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. model%potential%ng_split() .and. model%potential%is_updated(), &
         & more="setting the split dropped the potential update")
      if (allocated(error)) return
      call update_h2(model, err)
      if (.not. allocated(err)) call model%potential%add(charges, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. model%potential%is_updated(), more="adding a term kept the potential updated")
      if (allocated(error)) return
      call check(error, model%potential%n_terms(), 2)
   end subroutine check_solvent_copy

   !> UV tables of an LJ + H2 solute against the water solvent, for fixed
   !> solute charges, for host charges and for the host potential
   !> phi(r, a) = q_a/r with host charges, both read from the coupling
   !>
   !> - Column (iu - 1) nv + iv holds [iu, iv]
   !> - At alpha = 1: u_sr + ur_lr = u_LJ + q_u q_v/r, ur_lr = q_u q_v erf(r)/r
   !> - Potential not updated, before the model update and after a later `add`: refusal
   subroutine check_potential_tables(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(solvent_vv_type) :: solvent
      type(custom_lj_type) :: lj
      type(fixed_charges_type) :: fixed
      type(host_charges_type) :: host
      type(host_potential_type) :: field
      type(coupling_type), pointer :: coupling
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :), uk_lr(:, :), sr6(:), expected(:), phi(:, :)
      integer, allocatable :: pairs(:, :)
      real(wp), parameter :: sigma_u(2) = [2.5_wp, 3.5_wp], epsilon_u(2) = [1.0e-4_wp, 3.0e-4_wp]
      real(wp), parameter :: q_u(2) = [0.25_wp, -0.25_wp]
      real(wp) :: sigma, epsilon, qq
      integer :: iu, iv, ip, source

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_uniform_radial_pair(rgrid, kgrid, 64, 0.1_wp, err)
      if (.not. allocated(err)) call new_custom_lj(lj, sigma_u, epsilon_u, err)
      if (.not. allocated(err)) call new_fixed_charges(fixed, q_u, err)
      if (.not. allocated(err)) call new_water_solvent(ctx, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      do source = 1, 3
         call new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, err)
         if (.not. allocated(err)) call model%potential%add(lj, err)
         if (.not. allocated(err)) then
            select case (source)
            case (1)
               call model%potential%add(fixed, err)
            case (2)
               call model%potential%add(host, err)
            case default
               call model%potential%add(field, err)
            end select
         end if
         call check_moist_error(error, err)
         if (allocated(error)) return

         call model%potential%compute(model%rgrid, model%kgrid, u_sr, ur_lr, uk_lr, err, solvent=solvent%potential)
         call check(error, allocated(err), more="tables of a potential that is not updated")
         if (allocated(error)) return
         call check(error, index(err%message, "not updated") > 0, more=err%message)
         if (allocated(error)) return
         deallocate (err)

         call update_h2(model, err)
         if (.not. allocated(err)) call model%new_coupling(coupling, err)
         if (.not. allocated(err)) call model%potential%site_pairs(pairs, err, solvent=solvent%potential)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, all(shape(pairs) == [2, 6]), more="UV pair list is not (2, nu nv)")
         if (allocated(error)) return
         do iu = 1, 2
            do iv = 1, 3
               call check(error, all(pairs(:, (iu - 1)*3 + iv) == [iu, iv]), more="UV pair list is not solute major")
               if (allocated(error)) return
            end do
         end do

         call model%prepare_energy(coupling, err)
         if (source >= 2) then
            if (.not. allocated(err)) call answer_charges(coupling, q_u, err)
         end if
         if (source == 3) then
            allocate (phi(model%rgrid%npts, 2))
            do iu = 1, 2
               phi(:, iu) = q_u(iu)/model%rgrid%r
            end do
            if (.not. allocated(err)) call answer_radial_potential(coupling, phi, err)
         end if
         if (.not. allocated(err)) call model%potential%compute(model%rgrid, model%kgrid, u_sr, ur_lr, uk_lr, err, &
            & solvent=solvent%potential, coupling=coupling)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, all(shape(u_sr) == [model%rgrid%npts, 6]) .and. all(shape(uk_lr) == [model%kgrid%npts, 6]), &
            & more="UV tables are not (npts, nu nv)")
         if (allocated(error)) return
         do iu = 1, 2
            do iv = 1, 3
               ip = (iu - 1)*3 + iv
               sigma = 0.5_wp*(sigma_u(iu) + water_sigma(iv))
               epsilon = sqrt(epsilon_u(iu)*water_epsilon(iv))
               qq = q_u(iu)*water_charges(iv)
               sr6 = (sigma/model%rgrid%r)**6
               expected = 4.0_wp*epsilon*(sr6*sr6 - sr6) + qq/model%rgrid%r
               call check(error, maxval(abs(u_sr(:, ip) + ur_lr(:, ip) - expected)/max(1.0_wp, abs(expected))) &
                  & < 1.0e-12_wp, more="u_sr + ur_lr is not the full UV potential")
               if (allocated(error)) return
               call check(error, maxval(abs(ur_lr(:, ip) - qq*erf(model%rgrid%r)/model%rgrid%r)) < 1.0e-12_wp, &
                  & more="ur_lr is not the Ng real-space tail")
               if (allocated(error)) return
            end do
         end do

         call model%potential%add(lj, err)
         if (.not. allocated(err)) call model%potential%compute(model%rgrid, model%kgrid, u_sr, ur_lr, uk_lr, err, &
            & solvent=solvent%potential, coupling=coupling)
         call check(error, allocated(err), more="tables after an add without a new update")
         if (allocated(error)) return
         call check(error, index(err%message, "not updated") > 0, more=err%message)
         if (allocated(error)) return
         deallocate (err)
         call model%release_coupling(coupling)
      end do
   end subroutine check_potential_tables

   !> Compare unanswered walks of a fresh coupling
   !>
   !> - Energy, response and gradient phases
   !>
   !> @param[in,out] model  updated model with its solute terms
   !> @param[in]     walks  expected walk summaries of the three phases
   subroutine check_phase_walks(error, model, walks)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Updated model with its solute terms
      type(model_moz_1d_type), intent(inout), target :: model
      !> Expected walk summaries of the three phases
      character(len=*), intent(in) :: walks(:)
      type(moist_error), allocatable :: err
      type(coupling_type), pointer :: coupling
      !> Phase labels for diagnostics
      character(len=8), parameter :: phases(3) = [character(len=8) :: "energy", "response", "gradient"]
      !> Walk of the current phase
      character(len=:), allocatable :: summary
      integer :: phase

      call model%new_coupling(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      do phase = 1, size(phases)
         select case (phase)
         case (1)
            call model%prepare_energy(coupling, err)
         case (2)
            call model%prepare_response(coupling, err)
         case default
            call model%prepare_gradient(coupling, err)
         end select
         if (allocated(err)) then
            call test_failed(error, err%message)
            exit
         end if
         call walk_summary(coupling, summary)
         call check(error, summary, trim(walks(phase)), more=trim(phases(phase))//" phase")
         if (allocated(error)) exit
      end do
      call model%release_coupling(coupling)
   end subroutine check_phase_walks

   !> Collect one unanswered coupling walk
   !>
   !> - Visited requests and declared outputs; missing output marked `*`, e.g. "atomic_charges(q*);"
   !> - Subroutine output to avoid gfortran's static deferred-length function-result temporary
   !> - Shared length temporary: race between tests
   !>
   !> @param[in,out] coupling  staged coupling
   !> @param[out]    summary   visited requests and their outputs
   subroutine walk_summary(coupling, summary)
      !> Staged coupling
      type(coupling_type), intent(inout) :: coupling
      !> Visited requests and their outputs
      character(len=:), allocatable, intent(out) :: summary
      !> Outputs a MOZ term may declare, in listing order
      character(len=8), parameter :: outputs(5) = [character(len=8) :: &
         & "phi", "dphi_dr", "q", "mu", "theta"]
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

   !> Getter guards precede pending theory and preserve caller outputs
   !>
   !> Every getter meets a foreign coupling, a model that is not updated, an
   !> unstaged coupling, missing charges, the wrong staging and, past every
   !> guard, the pending theory; each refusal is named and leaves the energy,
   !> the gradient and the seeded response untouched
   subroutine check_getter_guards(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model, foreign
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling, other
      type(response_type) :: response
      type(atomic_charge_adjoint_response_type) :: seed
      real(wp) :: energy, gradient(3, 2)
      integer :: phase, scenario
      !> Getter and scenario of the current case for diagnostics
      character(len=64) :: label
      character(len=8), parameter :: phases(3) = [character(len=8) :: "energy", "response", "gradient"]
      !> Refusal each scenario must name
      character(len=40), parameter :: reasons(6) = [character(len=40) :: &
         & "different model", "updated first", "not staged", "missing required outputs", &
         & "staged for the", "is not implemented"]

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_updated_model(ctx, model, err)
      if (.not. allocated(err)) call new_updated_model(ctx, foreign, err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call foreign%new_coupling(other, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call new(mol, [1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      seed%dg_dq = [4.0_wp, -2.0_wp]
      do phase = 1, 3
         do scenario = 1, 6
            label = trim(phases(phase))//" getter, case '"//trim(reasons(scenario))//"'"
            ! Every case starts from an updated model; the update ends every walk
            call model%update(mol, err)
            if (allocated(err)) then
               call test_failed(error, err%message)
               return
            end if
            select case (scenario)
            case (2)
               call model%invalidate()
            case (4, 6)
               select case (phase)
               case (1)
                  call model%prepare_energy(coupling, err)
               case (2)
                  call model%prepare_response(coupling, err)
               case default
                  call model%prepare_gradient(coupling, err)
               end select
            case (5)
               if (phase == 1) then
                  call model%prepare_response(coupling, err)
               else
                  call model%prepare_energy(coupling, err)
               end if
            case default
               ! Foreign and unstaged scenarios require no preparation
            end select
            if (allocated(err)) then
               call test_failed(error, err%message)
               return
            end if
            if (scenario == 6) then
               call check(error, coupling%next(), more=trim(label)//": charges not pending")
               if (allocated(error)) return
               call coupling%answer("q", [0.3_wp, -0.3_wp], err)
               if (allocated(err)) then
                  call test_failed(error, err%message)
                  return
               end if
            end if
            call response_accumulate(response, seed, err)
            if (allocated(err)) then
               call test_failed(error, err%message)
               return
            end if
            energy = 7.0_wp
            gradient = 3.0_wp
            if (scenario == 1) then
               call invoke_getter(model, other, phase, response, energy, gradient, err)
            else
               call invoke_getter(model, coupling, phase, response, energy, gradient, err)
            end if
            call check(error, allocated(err), more=trim(label)//": accepted")
            if (allocated(error)) return
            call check(error, index(err%message, trim(reasons(scenario))) > 0, &
               & more=trim(label)//": "//err%message)
            if (allocated(error)) return
            if (scenario == 6) then
               call check(error, index(err%message, "1D MOZ "//trim(phases(phase))) > 0, &
                  & more=trim(label)//": "//err%message)
               if (allocated(error)) return
            end if
            call check(error, energy == 7.0_wp .and. all(gradient == 3.0_wp), &
               & more=trim(label)//": the refusal wrote the energy or gradient")
            if (allocated(error)) return
            ! A rejected request, including pending theory, keeps the seeded response
            call check(error, response%next(), &
               & more=trim(label)//": the refusal cleared the response")
            if (allocated(error)) return
            select type (item => response%item())
            type is (atomic_charge_adjoint_response_type)
               call check(error, all(item%dg_dq == real(scenario, wp)*seed%dg_dq), &
                  & more=trim(label)//": the refusal changed the response")
            class default
               call test_failed(error, trim(label)//": the refusal replaced the response type")
            end select
            if (allocated(error)) return
            call check(error, .not. response%next(), &
               & more=trim(label)//": the refusal appended a response item")
            if (allocated(error)) return
            deallocate (err)
         end do
         ! Restart the seed accumulation for the next getter channel
         block
            type(response_type) :: empty
            response = empty
         end block
      end do
      call model%release_coupling(coupling)
      call foreign%release_coupling(other)
      call check(error, .not. associated(coupling) .and. .not. associated(other))
   end subroutine check_getter_guards

   !> Invoke one model getter for a common guard table
   !>
   !> @param[in,out] model Model under test
   !> @param[in,out] coupling Owned or foreign coupling
   !> @param[in] phase Getter index
   !> @param[in,out] response Seeded response
   !> @param[in,out] energy Seeded energy
   !> @param[in,out] gradient Seeded gradient
   !> @param[out] err Getter error
   subroutine invoke_getter(model, coupling, phase, response, energy, gradient, err)
      !> Model under test
      type(model_moz_1d_type), intent(inout) :: model
      !> Owned or foreign coupling
      type(coupling_type), intent(inout), target :: coupling
      !> Getter index
      integer, intent(in) :: phase
      !> Seeded response
      type(response_type), intent(inout) :: response
      !> Seeded energy
      real(wp), intent(inout) :: energy
      !> Seeded gradient
      real(wp), intent(inout) :: gradient(:, :)
      !> Getter error
      type(moist_error), allocatable, intent(out) :: err
      select case (phase)
      case (1)
         call model%get_energy(coupling, energy, err)
      case (2)
         call model%get_response(coupling, response, err)
      case default
         call model%get_gradient(coupling, response, gradient, err)
      end select
   end subroutine invoke_getter

   !> Failed updates reset count and status and end the staged walk
   subroutine check_update_guards(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_1d_type), target :: model
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling

      call check(error, model%atom_count(), 0)
      if (allocated(error)) return
      call check(error, .not. model%is_updated())
      if (allocated(error)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call model%update(mol, err)
      call check(error, allocated(err), more="unconstructed update accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Construct the 1D MOZ model") > 0)
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_updated_model(ctx, model, err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      ! q stays unanswered; the full pass rewinds, so a staged walk offers it again
      call check(error, coupling%next(), more="atomic_charges must be pending")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="atomic_charges is the only request")
      if (allocated(error)) return
      call new(mol, [integer ::], reshape([real(wp) ::], [3, 0]))
      call model%update(mol, err)
      call check(error, allocated(err), more="empty structure accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "at least one atom") > 0, more=err%message)
      if (allocated(error)) return
      call check(error, model%atom_count(), 0)
      if (allocated(error)) return
      call check(error, .not. model%is_updated(), more="failed update kept valid status")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="failed update kept the staged walk")
      if (allocated(error)) return
      call model%release_coupling(coupling)
      call check(error, .not. associated(coupling))
   end subroutine check_update_guards

   !> Construct a 1D MOZ model on the water solvent, without solute terms
   !>
   !> @param[in]  ctx    run context, outlives the model
   !> @param[out] model  constructed model
   !> @param[out] err    construction error
   subroutine new_model(ctx, model, err)
      !> Run context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Constructed model
      type(model_moz_1d_type), intent(out) :: model
      !> Construction error
      type(moist_error), allocatable, intent(out) :: err
      type(solvent_vv_type) :: solvent
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      call new_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call new_radial_pair(rgrid, kgrid, err)
      if (allocated(err)) return
      call new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, err)
   end subroutine new_model

   !> Construct a 1D MOZ model with a host charges term
   !>
   !> @param[in]  ctx    run context, outlives the model
   !> @param[out] model  constructed model
   !> @param[out] err    construction error
   subroutine new_host_charges_model(ctx, model, err)
      !> Run context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Constructed model
      type(model_moz_1d_type), intent(out) :: model
      !> Construction error
      type(moist_error), allocatable, intent(out) :: err
      type(host_charges_type) :: charges
      call new_model(ctx, model, err)
      if (allocated(err)) return
      call model%potential%add(charges, err)
   end subroutine new_host_charges_model

   !> Construct a 1D MOZ model with a host charges term, updated to two hydrogen atoms
   !>
   !> @param[in]  ctx    run context, outlives the model
   !> @param[out] model  updated model
   !> @param[out] err    construction or update error
   subroutine new_updated_model(ctx, model, err)
      !> Run context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Updated model
      type(model_moz_1d_type), intent(out) :: model
      !> Construction or update error
      type(moist_error), allocatable, intent(out) :: err
      call new_host_charges_model(ctx, model, err)
      if (allocated(err)) return
      call update_h2(model, err)
   end subroutine new_updated_model

   !> Update a model to two hydrogen atoms one bohr apart
   !>
   !> @param[in,out] model  constructed model
   !> @param[out]    err    update error
   subroutine update_h2(model, err)
      !> Constructed model
      type(model_moz_1d_type), intent(inout) :: model
      !> Update error
      type(moist_error), allocatable, intent(out) :: err
      type(structure_type) :: mol
      call new(mol, [1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call model%update(mol, err)
   end subroutine update_h2

   !> Uniform radial pair of 16 points, 0.2 bohr apart
   !>
   !> @param[out] rgrid  real-space radial grid
   !> @param[out] kgrid  reciprocal radial grid
   !> @param[out] err    grid failure
   subroutine new_radial_pair(rgrid, kgrid, err)
      !> Real-space radial grid
      type(moist_math_grid_radial_type), intent(out) :: rgrid
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(out) :: kgrid
      !> Grid failure
      type(moist_error), allocatable, intent(out) :: err
      call new_uniform_radial_pair(rgrid, kgrid, 16, 0.2_wp, err)
   end subroutine new_radial_pair

end module test_moz_1d
