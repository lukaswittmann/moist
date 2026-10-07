!> 3D MOZ model and grid ownership tests
module test_moz_3d
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_gaussian_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid, new_molecular_gaussian_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   use test_helpers, only: get_qc_handymod_recipe
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
         & new_unittest("ec_requests_by_grid", check_ec_requests), &
         & new_unittest("qat_and_multipoles_requests", check_charge_sources), &
         & new_unittest("unknown_coupling_mode", check_unknown_mode)]
   end subroutine collect_moz_3d

   !> Drive a 3D MOZ model through the coupling protocol
   !>
   !> The model owns a copy of the configured grid, drives the
   !> gaussian_potential and atomic_charges requests of the default "ec"
   !> source through the usual coupling protocol, then reports pending theory
   !> once the mandatory outputs are answered; a wrong-shape answer is
   !> rejected first
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
      real(wp), allocatable :: phi(:)
      real(wp) :: energy
      integer :: i
      call new_context(ctx, nthreads=0, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_cartesian_gaussian_grid(template, err, nx=2, ny=2, nz=2, margin=0.0_wp)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      template%nx = 2
      template%ny = 2
      template%nz = 2
      template%dr = 0.5_wp
      call new_moz_3d_model(model, ctx, template, err)
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
      call check(error, coupling%next(), more="gaussian_potential must be pending")
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
      call check(error, coupling%next(), more="the tail charges must be pending")
      if (allocated(error)) return
      call coupling%answer("q", [0.0_wp], err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, .not. coupling%next(), more="gaussian_potential and atomic_charges are the only requests")
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

   !> Concrete grids work through the 3D MOZ model and preserve copy ownership
   !>
   !> Cartesian, molecular with Gaussian widths and bare molecular grids
   subroutine test_grid_model(error)
      !> Borrowed typed model grid
      class(moist_math_grid_3d_type), pointer :: model_grid
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_cartesian_type) :: cart
      !> Molecular grids with and without Gaussian widths
      type(moist_math_grid_3d_molecular_type) :: molecular, bare
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(model_moz_3d_type), target :: model
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(moist_math_grid_3d_type), allocatable :: grid
      type(field_query_type) :: query
      real(wp), allocatable :: original(:, :)
      integer :: kind, template_ngrid

      call new_context(ctx, nthreads=0, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_cartesian_gaussian_grid(cart, err, 4, 6, 8, 0.5_wp, margin=0.0_wp)
      call require_success(error, err)
      if (allocated(error)) return
      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_gaussian_grid(molecular, err, recipe=recipe)
      if (.not. allocated(err)) call new_molecular_point_grid(bare, err, recipe=recipe)
      call require_success(error, err)
      if (allocated(error)) return
      do kind = 1, 3
         select case (kind)
         case (1)
            allocate (grid, source=cart)
         case (2)
            allocate (grid, source=molecular)
         case default
            allocate (grid, source=bare)
         end select
         call new_moz_3d_model(model, ctx, grid, err)
         call require_success(error, err)
         if (allocated(error)) return
         call check(error, associated(model%ctx, ctx), "model must retain its borrowed context")
         if (allocated(error)) return
         template_ngrid = grid%ngrid
         call model%update(mol, err)
         call require_success(error, err)
         if (allocated(error)) return
         call check(error, grid%ngrid, template_ngrid, "updating the model must leave the template untouched")
         if (allocated(error)) return
         model_grid => model%grid
         call check(error, model_grid%natom, 1)
         if (allocated(error)) return
         call check(error, allocated(model_grid%xi0) .neqv. (kind == 3), &
            & more="only the bare molecular grid has no Gaussian widths")
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
         call require_success(error, err)
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

   !> "ec", the 3D default, declares the grid potential next to the tail
   !> charges; the grid picks the potential and whether its gradient needs
   !> `dphi_dxi`
   !>
   !> - Cartesian: Gaussian widths fixed by the spacing, no `dphi_dxi`
   !> - molecular with widths: they follow the weights, so `dphi_dxi` too
   !> - molecular without widths: a bare point potential
   subroutine check_ec_requests(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: cart
      !> Molecular grids with and without Gaussian widths
      type(moist_math_grid_3d_molecular_type) :: widths, points
      type(moist_math_grid_atomic_recipe_type) :: recipe
      !> Expected walks of the energy, response and gradient phases
      character(len=96) :: walks(3)
      integer :: kind

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cartesian_gaussian_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      call require_success(error, err)
      if (allocated(error)) return
      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_gaussian_grid(widths, err, recipe=recipe)
      if (.not. allocated(err)) call new_molecular_point_grid(points, err, recipe=recipe)
      call require_success(error, err)
      if (allocated(error)) return
      do kind = 1, 3
         select case (kind)
         case (1)
            call new_updated_model(ctx, cart, model, err)
            walks(1) = "gaussian_potential(phi*,dphi_dr,dphi_dxi);atomic_charges(q*);"
            walks(3) = "gaussian_potential(phi*,dphi_dr*,dphi_dxi);atomic_charges(q*);"
         case (2)
            call new_updated_model(ctx, widths, model, err)
            walks(1) = "gaussian_potential(phi*,dphi_dr,dphi_dxi);atomic_charges(q*);"
            walks(3) = "gaussian_potential(phi*,dphi_dr*,dphi_dxi*);atomic_charges(q*);"
         case default
            call new_updated_model(ctx, points, model, err)
            walks(1) = "point_potential(phi*,dphi_dr);atomic_charges(q*);"
            walks(3) = "point_potential(phi*,dphi_dr*);atomic_charges(q*);"
         end select
         walks(2) = walks(1)
         call require_success(error, err)
         if (allocated(error)) return
         call check(error, model%coupling_mode == "ec", more="the 3D default source is not 'ec'")
         if (allocated(error)) return
         call check_phase_walks(error, model, walks)
         if (allocated(error)) return
      end do
   end subroutine check_ec_requests

   !> "qat" declares only the partial charges and "multipoles" only the
   !> point multipoles, every output pending in every phase
   subroutine check_charge_sources(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: cart
      !> Expected walks of the energy, response and gradient phases
      character(len=96) :: walks(3)

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cartesian_gaussian_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_updated_model(ctx, cart, model, err)
      call require_success(error, err)
      if (allocated(error)) return
      model%coupling_mode = "qat"
      walks = "atomic_charges(q*);"
      call check_phase_walks(error, model, walks)
      if (allocated(error)) return
      model%coupling_mode = "multipoles"
      walks = "atomic_multipoles(q*,mu*,theta*);"
      call check_phase_walks(error, model, walks)
   end subroutine check_charge_sources

   !> An unknown coupling source is refused by name and leaves no coupling behind
   subroutine check_unknown_mode(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(model_moz_3d_type), target :: model
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(coupling_type), pointer :: coupling

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cartesian_gaussian_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_updated_model(ctx, cart, model, err)
      call require_success(error, err)
      if (allocated(error)) return
      model%coupling_mode = "bogus"
      call model%new_coupling(coupling, err)
      call check(error, allocated(err) .and. .not. associated(coupling), more="an unknown source was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "3D MOZ coupling_mode 'bogus' is not supported") > 0, more=err%message)
      if (allocated(error)) return

      ! A refused redeclaration must end the previously staged walk
      model%coupling_mode = "qat"
      call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      call require_success(error, err)
      if (allocated(error)) return
      model%coupling_mode = "bogus"
      call model%prepare_response(coupling, err)
      call check(error, allocated(err), more="an unknown source was staged")
      if (allocated(error)) return
      call check(error, .not. coupling%next(), more="the refused declaration left a staging behind")
      if (allocated(error)) return
      call model%release_coupling(coupling)
      call check(error, .not. associated(coupling), more="release left a coupling pointer")
   end subroutine check_unknown_mode

   !> Stage the energy, response and gradient phases of a fresh coupling and
   !> compare each walk, answering nothing, to the expected one
   !> @param[in,out] model Updated model with its coupling source set
   !> @param[in] walks Expected walk summaries of the three phases
   subroutine check_phase_walks(error, model, walks)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Updated model with its coupling source set
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
      call require_success(error, err)
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
         call require_success(error, err)
         if (allocated(error)) exit
         call walk_summary(coupling, summary)
         call check(error, summary, trim(walks(phase)), &
            & more=trim(model%coupling_mode)//" source, "//trim(phases(phase))//" phase")
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
      character(len=8), parameter :: outputs(6) = [character(len=8) :: &
         & "phi", "dphi_dr", "dphi_dxi", "q", "mu", "theta"]
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

   !> Construct a 3D MOZ model on a copy of `grid`, updated to one hydrogen atom
   !>
   !> @param[in] ctx Run context, outlives the model
   !> @param[in] grid Spatial grid template
   !> @param[out] model Updated model
   !> @param[out] err Construction or update error
   subroutine new_updated_model(ctx, grid, model, err)
      !> Run context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Spatial grid template
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Updated model
      type(model_moz_3d_type), intent(out) :: model
      !> Construction or update error
      type(moist_error), allocatable, intent(out) :: err
      type(structure_type) :: mol
      call new_moz_3d_model(model, ctx, grid, err)
      if (allocated(err)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call model%update(mol, err)
   end subroutine new_updated_model

   !> Forward a library error into the test framework
   !> @param[in] err Library error
   subroutine require_success(error, err)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable, intent(in) :: err
      if (allocated(err)) call test_failed(error, err%message)
   end subroutine require_success

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
      call new_cartesian_gaussian_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_updated_model(ctx, cart, model, err)
      if (.not. allocated(err)) call new_updated_model(ctx, cart, foreign, err)
      model%coupling_mode = "qat"
      foreign%coupling_mode = "qat"
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
      call new_cartesian_gaussian_grid(cart, err, 2, 2, 2, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_updated_model(ctx, cart, model, err)
      model%coupling_mode = "qat"
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_energy(coupling, err)
      call require_success(error, err)
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
