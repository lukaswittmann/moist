!> 3D MOZ model and grid ownership tests
module test_moz_3d
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_data_solvents, only: solvation_system_type
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent
   use moist_model_moz_potential_term, only: potential_term_type
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use moist_model_moz_potential_electrostatic_fixed, only: fixed_charges_type, new_fixed_charges
   use moist_model_moz_potential_electrostatic_host_charges, only: host_charges_type
   use moist_model_moz_potential_electrostatic_host_potential, only: host_potential_type
   use moist_model_moz_potential_electrostatic_multipole, only: host_multipoles_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid, &
      & new_cartesian_gaussian_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid, &
      & new_molecular_gaussian_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   use test_helpers, only: get_qc_handymod_recipe, check_moist_error
   use test_moz_fixtures, only: water_sigma, water_epsilon, water_charges, new_water_system, answer_charges, &
      & new_water_solvent, refusing_term_type
   use moist_channels_fields, only: field_query_type
   use moist_channels_coupling, only: coupling_type, coupling_request_type
   use moist_channels_response, only: response_type, atomic_charge_adjoint_response_type, response_accumulate
   implicit none(type, external)
   private
   public :: collect_moz_3d

contains
   !> Register 3D MOZ tests
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_moz_3d(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("contract", check_moz_3d), &
         & new_unittest("getter_guards", check_getter_guards), &
         & new_unittest("update_guards", check_update_guards), &
         & new_unittest("owned_grid", test_grid_model), &
         & new_unittest("host_potential_point_requests", check_host_potential_requests), &
         & new_unittest("gaussian_grid_refused", check_gaussian_grid_refused), &
         & new_unittest("host_charges_and_multipoles", check_charge_sources), &
         & new_unittest("term_failures", check_term_failures), &
         & new_unittest("solvent_copy", check_solvent_copy), &
         & new_unittest("uv_potential_tables", check_potential_tables)]
   end subroutine collect_moz_3d

   !> Drive a 3D MOZ model through the coupling protocol
   !>
   !> - Owned copy of configured point grid
   !> - Host potential: atomic_charges and point_potential requests through coupling protocol
   !> - Theory-pending error after all mandatory answers
   !> - Wrong-shape answer rejected first
   subroutine check_moz_3d(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: template
      class(moist_math_grid_3d_type), pointer :: grid
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      type(solvent_vv_type) :: solvent
      type(host_potential_type) :: potential
      real(wp), allocatable :: phi(:)
      real(wp) :: energy
      integer :: i
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_water_solvent(ctx, solvent, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_cartesian_point_grid(template, err, nx=2, ny=2, nz=2, margin=0.0_wp)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      template%nx = 2
      template%ny = 2
      template%nz = 2
      template%dr = 0.5_wp
      call new_moz_3d_model(model, ctx, template, solvent, err)
      if (.not. allocated(err)) call model%potential%add(potential, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      grid => model%grid
      call check(error, associated(grid) .and. model%is_updated())
      if (allocated(error)) return
      call check(error, grid%ngrid == 8 .and. model%atom_count() == 1)
      if (allocated(error)) return
      template%nx = 3
      call check(error, grid%ngrid == 8, more="model must own an independent grid copy")
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
      call check(error, coupling%next(), more="the tail charges must be pending")
      if (allocated(error)) return
      call coupling%answer("q", [0.0_wp], err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, coupling%next(), more="point_potential must be pending")
      if (allocated(error)) return

      ! A wrong-shape answer is rejected and leaves the output missing
      call coupling%answer("phi", [1.0_wp], err)
      call check(error, allocated(err))
      if (allocated(error)) return

      ! The correctly shaped answer (one value per grid point) completes the walk
      allocate (phi(grid%ngrid))
      phi = [(0.1_wp*real(i, wp), i=1, grid%ngrid)]
      call coupling%answer("phi", phi, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, .not. coupling%next(), more="atomic_charges and point_potential are the only requests")
      if (allocated(error)) return

      energy = 9.0_wp
      call model%get_energy(coupling, energy, err)
      call check(error, allocated(err) .and. energy == 9.0_wp)
      if (allocated(error)) return
      call check(error, index(err%message, "3D MOZ energy") > 0)
      if (allocated(error)) return

      call model%prepare_response(coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call model%get_response(coupling, response, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      call check(error, index(err%message, "3D MOZ response") > 0)
      if (allocated(error)) return

      call model%release_coupling(coupling)
   end subroutine check_moz_3d

   !> Check point-grid operation and copy ownership
   !>
   !> - Cartesian and molecular point grids
   subroutine test_grid_model(error)
      !> Borrowed typed model grid
      class(moist_math_grid_3d_type), pointer :: model_grid
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_cartesian_type) :: cart
      !> Molecular point grid
      type(moist_math_grid_3d_molecular_type) :: bare
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(model_moz_3d_type), target :: model
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(moist_math_grid_3d_type), allocatable :: grid
      type(field_query_type) :: query
      real(wp), allocatable :: original(:, :)
      type(solvent_vv_type) :: solvent
      type(host_charges_type) :: charges
      integer :: kind, template_ngrid

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_water_solvent(ctx, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_cartesian_point_grid(cart, err, 4, 6, 8, 0.5_wp, margin=0.0_wp)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(bare, err, recipe=recipe)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do kind = 1, 2
         select case (kind)
         case (1)
            allocate (grid, source=cart)
         case default
            allocate (grid, source=bare)
         end select
         call new_moz_3d_model(model, ctx, grid, solvent, err)
         if (.not. allocated(err)) call model%potential%add(charges, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, associated(model%ctx, ctx), "model must retain its borrowed context")
         if (allocated(error)) return
         template_ngrid = grid%ngrid
         call model%update(mol, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, grid%ngrid, template_ngrid, "updating the model must leave the template untouched")
         if (allocated(error)) return
         model_grid => model%grid
         call check(error, model_grid%natom, 1)
         if (allocated(error)) return
         call check(error, .not. allocated(model_grid%xi0), more="a point grid carries Gaussian widths")
         if (allocated(error)) return
         call check_grid_fields(error, model)
         if (allocated(error)) return
         call query%fetch("w")
         call model%list_fields(query)
         call check(error, query%found)
         if (allocated(error)) return
         model_grid => model%grid
         call check(error, size(query%rvals), model_grid%ngrid)
         if (allocated(error)) return
         model_grid => model%grid
         original = model_grid%xyz
         mol%xyz(1, 1) = mol%xyz(1, 1) + 0.25_wp
         call grid%update(mol, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         model_grid => model%grid
         call check(error, all(model_grid%xyz == original), "model must own an independent grid copy")
         if (allocated(error)) return
         model_grid => model%grid
         select type (g => model_grid)
         type is (moist_math_grid_3d_cartesian_type)
            call check(error, g%ngrid, 4*6*8)
         type is (moist_math_grid_3d_molecular_type)
            call check(error, g%has_reciprocal())
         class default
            call test_failed(error, "model did not retain the concrete grid grid")
         end select
         if (allocated(error)) return
         deallocate (grid)
      end do
   end subroutine test_grid_model

   !> Check host-potential declarations on point grids
   !>
   !> - Tail charges, then point potential at grid points
   !> - phi in every phase; dphi_dr in gradient phase
   !> - Cartesian and molecular point grids
   subroutine check_host_potential_requests(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: cart
      !> Molecular point grid
      type(moist_math_grid_3d_molecular_type) :: points
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(host_potential_type) :: potential
      !> Expected walks of the energy, response and gradient phases
      character(len=96) :: walks(3)
      integer :: kind

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cartesian_point_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(points, err, recipe=recipe)
      call check_moist_error(error, err)
      if (allocated(error)) return
      walks(1) = "atomic_charges(q*);point_potential(phi*,dphi_dr);"
      walks(2) = walks(1)
      walks(3) = "atomic_charges(q*);point_potential(phi*,dphi_dr*);"
      do kind = 1, 2
         if (kind == 1) then
            call new_updated_model(ctx, cart, model, err, potential)
         else
            call new_updated_model(ctx, points, model, err, potential)
         end if
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check_phase_walks(error, model, walks)
         if (allocated(error)) return
      end do
   end subroutine check_host_potential_requests

   !> Check Gaussian-width grid refusal
   !>
   !> - Realized widths: refusal at construction
   !> - Geometry-dependent widths: refusal at update
   subroutine check_gaussian_grid_refused(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(moist_math_grid_3d_molecular_type) :: widths
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(solvent_vv_type) :: solvent
      type(host_charges_type) :: charges
      type(structure_type) :: mol
      character(len=*), parameter :: refusal = &
         & "MOZ 3D model supports point grids only; Gaussian-width grids are not supported yet"

      call new_context(ctx, nthreads=0, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call new_cartesian_gaussian_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_gaussian_grid(widths, err, recipe=recipe)
      call check_moist_error(error, err)
      if (allocated(error)) return

      ! Template without geometry: accepted at construction, refused at update
      call new_moz_3d_model(model, ctx, cart, solvent, err)
      if (.not. allocated(err)) call model%potential%add(charges, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call model%update(mol, err)
      call check(error, allocated(err) .and. .not. model%is_updated(), more="a Gaussian grid was updated")
      if (allocated(error)) return
      call check(error, err%message, refusal)
      if (allocated(error)) return
      deallocate (err)

      ! Realized widths: refused at construction
      call widths%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new_moz_3d_model(model, ctx, widths, solvent, err)
      call check(error, allocated(err), more="a grid with Gaussian widths was accepted")
      if (allocated(error)) return
      call check(error, err%message, refusal)
   end subroutine check_gaussian_grid_refused

   !> Check host term declarations and multipole refusal
   !>
   !> - Host charges: partial charges pending in every phase
   !> - Separate term scopes; shared charge request for identical inputs
   !> - Host multipole stub: refusal at update
   subroutine check_charge_sources(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(host_charges_type) :: charges
      type(host_potential_type) :: potential
      type(host_multipoles_type) :: multipoles
      type(structure_type) :: mol
      !> Expected walks of the energy, response and gradient phases
      character(len=96) :: walks(3)

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cartesian_point_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_updated_model(ctx, cart, model, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      walks = "atomic_charges(q*);"
      call check_phase_walks(error, model, walks)
      if (allocated(error)) return

      call new_updated_model(ctx, cart, model, err, potential)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      if (.not. allocated(err)) call model%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      ! Host potential requests: charges and point probe
      walks(1) = "atomic_charges(q*);point_potential(phi*,dphi_dr);"
      walks(2) = walks(1)
      walks(3) = "atomic_charges(q*);point_potential(phi*,dphi_dr*);"
      call check_phase_walks(error, model, walks)
      if (allocated(error)) return

      ! First charge supplier owns each atom
      ! Host potential after another charge term: field without owned charges
      call new_updated_model(ctx, cart, model, err, charges)
      if (.not. allocated(err)) call model%potential%add(potential, err)
      if (.not. allocated(err)) call model%update(mol, err)
      call check(error, allocated(err), more="a host potential behind another charge term was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "must supply every charge") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      call new_updated_model(ctx, cart, model, err, multipoles)
      call check(error, allocated(err) .and. .not. model%is_updated(), more="the multipoles stub was built")
      if (allocated(error)) return
      call check(error, index(err%message, "host_multipoles") > 0, more=err%message)
   end subroutine check_charge_sources

   !> Check named solute failures and cleared coupling state
   !>
   !> - Solute potential without terms: refusal at update
   !> - Failed term declaration: coupling refused, previously staged walk ended
   subroutine check_term_failures(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(solvent_vv_type) :: solvent
      type(refusing_term_type) :: refusing
      !> Declaration switch shared with the refusing term
      logical, target :: refuse
      type(coupling_type), pointer :: coupling
      type(structure_type) :: mol

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cartesian_point_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call new_moz_3d_model(model, ctx, cart, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call model%update(mol, err)
      call check(error, allocated(err) .and. .not. model%is_updated(), more="an empty solute potential was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "has no terms") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      ! Host charges first, then pair-only term with switchable declaration failure
      call new_updated_model(ctx, cart, model, err)
      refuse = .false.
      refusing%refuse => refuse
      if (.not. allocated(err)) call model%potential%add(refusing, err)
      if (.not. allocated(err)) call model%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      refuse = .true.
      call model%new_coupling(coupling, err)
      call check(error, allocated(err) .and. .not. associated(coupling), more="a failing declaration was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "refused declaration") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      ! Failed redeclaration: previous staging cleared
      refuse = .false.
      call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      refuse = .true.
      call model%prepare_response(coupling, err)
      call check(error, allocated(err), more="a failing declaration was staged")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="the refused declaration left a staging behind")
      if (allocated(error)) return
      call model%release_coupling(coupling)
      call check(error, .not. associated(coupling), more="release left a coupling pointer")
   end subroutine check_term_failures

   !> Check solvent ownership and model invalidation
   !>
   !> - Built solvent required; independent copy in model
   !> - Added term: model invalidation
   subroutine check_solvent_copy(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(solvation_system_type) :: system
      type(solvent_vv_type) :: solvent
      type(moz_potential_type) :: solvent_pot
      type(fixed_charges_type) :: extra
      type(host_charges_type) :: charges

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cartesian_point_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_water_system(system, err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new_moz_3d_model(model, ctx, cart, solvent, err)
      call check(error, allocated(err), more="a solvent that is not updated was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Update the 1D VV solvent") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      call new_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call new_updated_model(ctx, cart, model, err)
      if (.not. allocated(err)) call new_fixed_charges(extra, [0.1_wp, 0.1_wp, -0.2_wp], err)
      if (.not. allocated(err)) call solvent%potential%add(extra, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      solvent_pot = solvent%potential
      call check(error, solvent_pot%n_terms() == 3 .and. .not. solvent%potential%is_updated(), &
         & more="the extra term did not reach the original solvent")
      if (allocated(error)) return
      call check(error, model%solvent_nsite() == 3, more="the model solvent lost its sites")
      if (allocated(error)) return
      call check(error, model%potential%ng_split() .and. model%potential%alpha() == 1.0_wp, more="Ng split defaults changed")
      if (allocated(error)) return
      call check(error, model%is_updated())
      if (allocated(error)) return
      call model%potential%add(charges, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. model%potential%is_updated(), more="adding a term kept the potential updated")
      if (allocated(error)) return
      call check(error, model%potential%n_terms(), 2)
   end subroutine check_solvent_copy

   !> UV tables on a Cartesian grid against the water solvent at alpha = 1,
   !> for fixed solute charges and for host charges read from the coupling
   !>
   !> - Per solvent site v: u_sr + ur_lr = sum_a u_LJ(r_a) + q_v q_a/r_a
   !> - ur_lr = q_v sum_a q_a erf(r_a)/r_a
   !> - Potential not updated, before the model update: refusal
   subroutine check_potential_tables(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(custom_lj_type) :: lj
      type(fixed_charges_type) :: fixed
      type(host_charges_type) :: host
      type(coupling_type), pointer :: coupling
      type(structure_type) :: mol
      type(solvent_vv_type) :: solvent
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :)
      complex(wp), allocatable :: uk_lr(:, :)
      real(wp), parameter :: sigma_u(2) = [2.5_wp, 3.5_wp], epsilon_u(2) = [1.0e-4_wp, 3.0e-4_wp]
      real(wp), parameter :: q_u(2) = [0.25_wp, -0.25_wp]
      integer :: source

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cartesian_point_grid(cart, err, 4, 4, 4, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_custom_lj(lj, sigma_u, epsilon_u, err)
      if (.not. allocated(err)) call new_fixed_charges(fixed, q_u, err)
      if (.not. allocated(err)) call new_water_solvent(ctx, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      ! Off the grid points, which sit at half spacings around the solute
      call new(mol, [1, 1], reshape([0.1_wp, 0.05_wp, 0.0_wp, 1.3_wp, 0.0_wp, 0.1_wp], [3, 2]))

      do source = 1, 2
         call new_moz_3d_model(model, ctx, cart, solvent, err)
         if (.not. allocated(err)) call model%potential%add(lj, err)
         if (.not. allocated(err)) then
            if (source == 1) then
               call model%potential%add(fixed, err)
            else
               call model%potential%add(host, err)
            end if
         end if
         call check_moist_error(error, err)
         if (allocated(error)) return

         call model%potential%compute(cart, u_sr, ur_lr, uk_lr, err, solvent=solvent%potential)
         call check(error, allocated(err), more="tables of a potential that is not updated")
         if (allocated(error)) return
         call check(error, index(err%message, "not updated") > 0, more=err%message)
         if (allocated(error)) return
         deallocate (err)

         call model%update(mol, err)
         if (.not. allocated(err)) call model%new_coupling(coupling, err)
         if (.not. allocated(err)) call model%prepare_energy(coupling, err)
         if (source == 2) then
            if (.not. allocated(err)) call answer_charges(coupling, q_u, err)
         end if
         if (.not. allocated(err)) call model%potential%compute(model%grid, u_sr, ur_lr, uk_lr, err, &
            & solvent=solvent%potential, coupling=coupling)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check_tables(error, model, mol, sigma_u, epsilon_u, q_u, u_sr, ur_lr, uk_lr)
         if (allocated(error)) return
         call model%release_coupling(coupling)
      end do
   end subroutine check_potential_tables

   !> Compare 3D tables with the analytic LJ + point-charge potential
   !>
   !> @param[in] model      updated model
   !> @param[in] mol        solute structure
   !> @param[in] sigma_u    solute sigma per atom, bohr
   !> @param[in] epsilon_u  solute epsilon per atom, Hartree
   !> @param[in] q_u        solute charges per atom, e
   !> @param[in] u_sr       short-range table (ngrid, nv)
   !> @param[in] ur_lr      long-range real-space table (ngrid, nv)
   !> @param[in] uk_lr      long-range reciprocal table (npts_k, nv)
   subroutine check_tables(error, model, mol, sigma_u, epsilon_u, q_u, u_sr, ur_lr, uk_lr)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Updated model
      type(model_moz_3d_type), intent(in) :: model
      !> Solute structure
      type(structure_type), intent(in) :: mol
      !> Solute sigma per atom, bohr
      real(wp), intent(in) :: sigma_u(:)
      !> Solute epsilon per atom, Hartree
      real(wp), intent(in) :: epsilon_u(:)
      !> Solute charges per atom, e
      real(wp), intent(in) :: q_u(:)
      !> Short-range table (ngrid, nv)
      real(wp), intent(in) :: u_sr(:, :)
      !> Long-range real-space table (ngrid, nv)
      real(wp), intent(in) :: ur_lr(:, :)
      !> Long-range reciprocal table (npts_k, nv)
      complex(wp), intent(in) :: uk_lr(:, :)
      real(wp) :: r, sr6, sigma, epsilon, full, tail, worst_full, worst_tail
      integer :: ig, iv, ia
      call check(error, all(shape(u_sr) == [model%grid%ngrid, 3]) .and. all(shape(ur_lr) == shape(u_sr)) &
         & .and. all(shape(uk_lr) == [model%grid%npts_k, 3]), more="3D tables are not (ngrid, nv)")
      if (allocated(error)) return
      worst_full = 0.0_wp
      worst_tail = 0.0_wp
      do iv = 1, 3
         do ig = 1, model%grid%ngrid
            full = 0.0_wp
            tail = 0.0_wp
            do ia = 1, mol%nat
               r = norm2(model%grid%xyz(:, ig) - mol%xyz(:, ia))
               sigma = 0.5_wp*(sigma_u(ia) + water_sigma(iv))
               epsilon = sqrt(epsilon_u(ia)*water_epsilon(iv))
               sr6 = (sigma/r)**6
               full = full + 4.0_wp*epsilon*(sr6*sr6 - sr6) + water_charges(iv)*q_u(ia)/r
               tail = tail + water_charges(iv)*q_u(ia)*erf(r)/r
            end do
            worst_full = max(worst_full, abs(u_sr(ig, iv) + ur_lr(ig, iv) - full)/max(1.0_wp, abs(full)))
            worst_tail = max(worst_tail, abs(ur_lr(ig, iv) - tail))
         end do
      end do
      call check(error, worst_full < 1.0e-12_wp, more="u_sr + ur_lr is not the full UV potential")
      if (allocated(error)) return
      call check(error, worst_tail < 1.0e-12_wp, more="ur_lr is not the Ng real-space tail")
   end subroutine check_tables

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
      type(model_moz_3d_type), intent(inout), target :: model
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
      call check_moist_error(error, err)
      if (allocated(error)) return
      do phase = 1, size(phases)
         select case (phase)
         case (1)
            call model%prepare_energy(coupling, err)
         case (2)
            call model%prepare_response(coupling, err)
         case default
            call model%prepare_gradient(coupling, err)
         end select
         call check_moist_error(error, err)
         if (allocated(error)) exit
         call walk_summary(coupling, summary)
         call check(error, summary, trim(walks(phase)), &
            & more=trim(phases(phase))//" phase")
         if (allocated(error)) exit
      end do
      call model%release_coupling(coupling)
   end subroutine check_phase_walks

   !> One pass of the walk, answering nothing: each visited request with its
   !> declared outputs, a missing one marked `*`, e.g. "atomic_charges(q*);"
   !>
   !> A subroutine: gfortran returns a deferred-length function result through
   !> a static (thread-shared) length temporary, which races between tests
   !>
   !> @param[in,out] coupling Staged coupling
   !> @param[out] summary Visited requests and their outputs
   subroutine walk_summary(coupling, summary)
      !> Staged coupling
      type(coupling_type), intent(inout) :: coupling
      !> Visited requests and their outputs
      character(len=:), allocatable, intent(out) :: summary
      !> Outputs a MOZ source may declare, in listing order
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

   !> Construct a 3D MOZ model on copied grid and water solvent
   !>
   !> - One solute term, updated to one hydrogen atom
   !>
   !> @param[in]  ctx    run context, outlives the model
   !> @param[in]  grid   spatial grid template
   !> @param[out] model  updated model
   !> @param[out] err    construction or update error
   !> @param[in]  term   optional solute term; host charges when absent
   subroutine new_updated_model(ctx, grid, model, err, term)
      !> Run context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Spatial grid template
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Updated model
      type(model_moz_3d_type), intent(out) :: model
      !> Construction or update error
      type(moist_error), allocatable, intent(out) :: err
      !> Optional solute term; host charges when absent
      class(potential_term_type), intent(in), optional :: term
      type(solvent_vv_type) :: solvent
      type(host_charges_type) :: charges
      type(structure_type) :: mol
      call new_water_solvent(ctx, solvent, err)
      if (allocated(err)) return
      call new_moz_3d_model(model, ctx, grid, solvent, err)
      if (allocated(err)) return
      if (present(term)) then
         call model%potential%add(term, err)
      else
         call model%potential%add(charges, err)
      end if
      if (allocated(err)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call model%update(mol, err)
   end subroutine new_updated_model


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
      type(model_moz_3d_type), target :: model, foreign
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(structure_type) :: mol
      type(coupling_type), pointer :: coupling, other
      type(response_type) :: response
      type(atomic_charge_adjoint_response_type) :: seed
      real(wp) :: energy, gradient(3, 1)
      integer :: phase, scenario
      !> Getter and scenario of the current case for diagnostics
      character(len=64) :: label
      character(len=8), parameter :: phases(3) = [character(len=8) :: "energy", "response", "gradient"]
      !> Refusal each scenario must name
      character(len=40), parameter :: reasons(6) = [character(len=40) :: &
         & "different model", "updated first", "not staged", "missing required outputs", &
         & "staged for the", "is not implemented"]

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cartesian_point_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_updated_model(ctx, cart, model, err)
      if (.not. allocated(err)) call new_updated_model(ctx, cart, foreign, err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call foreign%new_coupling(other, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      seed%dg_dq = [4.0_wp]
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
               call coupling%answer("q", [0.3_wp], err)
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
               call check(error, index(err%message, "3D MOZ "//trim(phases(phase))) > 0, &
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
      type(model_moz_3d_type), intent(inout) :: model
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

   !> Failed grid updates end the staged walk and retain an invalid model
   subroutine check_update_guards(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: cart
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
      call check(error, index(err%message, "Construct the 3D MOZ model") > 0)
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cartesian_point_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_updated_model(ctx, cart, model, err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      ! q stays unanswered; the full pass rewinds, so a staged walk offers it again
      call check(error, coupling%next(), more="atomic_charges must be pending")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="atomic_charges is the only request")
      if (allocated(error)) return
      call new(mol, [integer ::], reshape([real(wp) ::], [3, 0]))
      call model%update(mol, err)
      call check(error, allocated(err), more="empty grid update accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "at least one solute atom") > 0, more=err%message)
      if (allocated(error)) return
      call check(error, .not. model%is_updated(), more="failed update kept valid status")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="failed update kept the staged walk")
      if (allocated(error)) return
      call model%release_coupling(coupling)
      call check(error, .not. associated(coupling))
   end subroutine check_update_guards

   !> Named model fields match the owned domain in both concrete grid types
   !>
   !> A grid without Gaussian widths declares no `xi0`
   !> @param[in] model Updated model
   subroutine check_grid_fields(error, model)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Updated model
      type(model_moz_3d_type), intent(in) :: model
      type(field_query_type) :: query
      character(len=5), parameter :: names(6) = [character(len=5) :: "ngrid", "natom", "xyz", "w", "xi0", "owner"]
      integer :: i
      do i = 1, size(names)
         call query%fetch(trim(names(i)))
         call model%list_fields(query)
         if (i == 5 .and. .not. allocated(model%grid%xi0)) then
            call check(error, .not. query%found, more="unallocated xi0 was declared")
            if (allocated(error)) return
            cycle
         end if
         call check(error, query%found, more="missing grid field "//trim(names(i)))
         if (allocated(error)) return
         select case (i)
         case (1)
            call check(error, query%ivals(1), model%grid%ngrid)
         case (2)
            call check(error, query%ivals(1), model%grid%natom)
         case (3)
            call check(error, all(query%rvals == reshape(model%grid%xyz, [size(model%grid%xyz)])))
         case (4)
            call check(error, all(query%rvals == model%grid%w))
         case (5)
            call check(error, all(query%rvals == model%grid%xi0))
         case default
            call check(error, all(query%ivals == model%grid%owner - 1))
         end select
         if (allocated(error)) return
      end do
   end subroutine check_grid_fields

end module test_moz_3d
