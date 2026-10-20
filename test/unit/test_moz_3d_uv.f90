!> Tests of the 3D MOZ (UV) settings and solvent response table: moz_3d_uv
!>
!> SPC/E water as the solvent on a small VV grid; Cartesian and atom-centred grids
!> around a water molecule; each test runs as a selected test
module test_moz_3d_uv
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_model_type, only: solvation_model_type
   use moist_math_packed_sym, only: packed_index_lut
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type, moist_math_grid_3d_trafo_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type, new_uv3d_susceptibility
   use moist_model_moz_3d_equation, only: uv3d_equation_type
   use moist_model_moz_3d_driver, only: uv3d_problem_type, new_uv3d_problem
   use moist_model_moz_closure, only: moist_model_moz_closure_hnc_type
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_model_moz, only: model_moz_3d_thermodynamics_type
   use moist_model_moz_potential, only: lj_cspce, coulomb_spce, custom_lj_type, new_custom_lj, &
      & fixed_charges_type, new_fixed_charges
   use moist_model_moz_thermodynamics_kernel, only: coulomb_box_mean, electrostatic_finite_size_correction
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type
   use mctc_io_convert, only: autokcal, aatoau
   use moist_model_moz_solvent_vv, only: solvent_vv_type, solvent_vv_parameters_type
   use test_helpers, only: check_moist_error, get_qc_handymod_recipe
   use test_moz_fixtures, only: new_cspce_water_solvent
   implicit none(type, external)
   private
   public :: collect_moz_3d_uv

   !> Parameter file of the round-trip test
   character(len=*), parameter :: parameter_file = "moist-test-moz-3d-parameters.toml"

   !> Custom local integrand must survive the block producer's polymorphic dispatch
   type, extends(moist_model_moz_closure_hnc_type) :: energy_probe_closure_type
   contains
      procedure :: local_integrand => energy_probe_integrand
   end type energy_probe_closure_type

contains

   !> Register tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_3d_uv(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("parameters", check_parameters), &
         & new_unittest("susceptibility_cartesian", check_susceptibility_cartesian), &
         & new_unittest("susceptibility_molecular", check_susceptibility_molecular), &
         & new_unittest("susceptibility_refusals", check_susceptibility_refusals), &
         & new_unittest("equation_map", check_equation_map), &
         & new_unittest("equation_copy", check_equation_copy), &
         & new_unittest("equation_setup_errors", check_equation_setup_errors), &
         & new_unittest("energy_cartesian", check_energy_cartesian), &
         & new_unittest("finite_size_box", check_finite_size_box), &
         & new_unittest("energy_threads", check_energy_threads), &
         & new_unittest("energy_blocks", check_energy_blocks), &
         & new_unittest("energy_refusals", check_energy_refusals)]
   end subroutine collect_moz_3d_uv

   !> Bounded energy production agrees with full-field quadrature for all closures
   subroutine check_energy_blocks(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(structure_type) :: mol
      type(uv3d_susceptibility_type) :: table
      type(uv3d_problem_type) :: problem
      type(model_moz_3d_parameters_type) :: params
      character(len=4), parameter :: closures(6) = ["HNC ", "KH  ", "PSE2", "PSE3", "PSE4", "HNC "]
      integer, parameter :: nv = 8
      real(wp), allocatable :: sr(:, :), lr(:, :), h(:, :), c(:, :), phi(:, :)
      complex(wp), allocatable :: uk(:, :)
      real(wp) :: reference, integral, actual, serial
      integer :: i, col, ic, nt
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=4, verbosity=0)
      ! More than 65536 points, with a final partial block and non-dyadic dV.
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_cartesian_point_grid(grid, err, 71, 37, 29, 0.37_wp, margin=0.0_wp)
      if (.not. allocated(err)) call grid%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      ! Eight sites also exercise closure sweeps above their parallel threshold.
      table%nv = nv
      allocate (table%values(nv, nv, 1), source=0.0_wp)
      allocate (table%index(grid%npts_k), source=1)
      allocate (table%radius(1), source=0.0_wp)
      allocate (sr(grid%ngrid, nv), lr(grid%ngrid, nv), h(grid%ngrid, nv), c(grid%ngrid, nv), phi(grid%ngrid, nv))
      allocate (uk(grid%npts_k, nv), source=cmplx(0.0_wp, 0.0_wp, wp))
      sr = 0.0_wp
      do col = 1, nv
         do i = 1, grid%ngrid
            h(i, col) = 0.2_wp*sin(0.013_wp*i + col)
            c(i, col) = 0.03_wp*cos(0.031_wp*i - col)
            lr(i, col) = 0.01_wp*sin(0.007_wp*i*col)
         end do
      end do
      do ic = 1, size(closures)
         params%closure = closures(ic)
         call new_uv3d_problem(problem, ctx, grid, table, sr, lr, uk, [(0.003_wp + 0.002_wp*col, col=1, nv)], &
            & 1.7_wp, params, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         if (ic == 6) then
            deallocate (problem%closure)
            allocate (energy_probe_closure_type :: problem%closure)
         end if
         problem%cr = c + problem%beta*lr
         problem%tr = h - problem%cr
         call problem%closure%free_energy_integrands(h, c, phi, err)
         reference = 0.0_wp
         do col = 1, nv
            if (.not. allocated(err)) call grid%integrate_field(phi(:, col), integral, err)
            if (allocated(err)) exit
            reference = reference + problem%rho(col)*integral
         end do
         reference = reference/problem%beta
         call check_moist_error(error, err)
         if (allocated(error)) return
         do nt = 1, 4, 3
            problem%equation%nthreads = nt
            problem%closure%nthreads = nt
            call problem%free_energy(actual, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, actual, reference, thr=2.0e-12_wp*max(1.0_wp, abs(reference)), &
               & more="blocked energy differs from full-field quadrature: "//closures(ic))
            if (allocated(error)) return
            if (nt == 1) serial = actual
            call check(error, actual == serial, more="energy block merge depends on worker count")
            if (allocated(error)) return
         end do
      end do
      problem%tr(grid%ngrid, 2) = -2.0_wp - problem%cr(grid%ngrid, 2)
      call problem%free_energy(actual, err)
      call check(error, allocated(err), more="last block accepted negative g")
      if (allocated(error)) return
      problem%tr = h - problem%cr
      problem%cr(1, 1) = ieee_value(0.0_wp, ieee_quiet_nan)
      call problem%free_energy(actual, err)
      call check(error, allocated(err), more="energy block accepted non-finite correlations")
      call ctx%delete()
   end subroutine check_energy_blocks

   pure function energy_probe_integrand(self, h, c) result(phi)
      class(energy_probe_closure_type), intent(in) :: self
      real(wp), intent(in) :: h, c
      real(wp) :: phi
      phi = h + 2.0_wp*c
   end function energy_probe_integrand

   !> Check the 3D settings
   !>
   !> - Defaults: KH, no correction, OZ response, JFNK without preconditioner
   !> - TOML round trip, including the nested solver and solvent groups
   !> - Unknown response modes, preconditioners and corrections are refused
   subroutine check_parameters(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_3d_parameters_type) :: params, copy
      integer :: unit, stat

      call check(error, params%closure == "KH" .and. params%correction == "none" .and. &
         & params%finite_size == "none" .and. &
         & params%susceptibility == "oz" .and. params%solver%name == "JFNK" .and. params%precond == "none" .and. &
         & params%solvent%tolerance == 1.0e-11_wp .and. params%reuse_initial_guess .and. params%cold_guess == "cs0", &
         & more="unexpected 3D defaults")
      if (allocated(error)) return
      call params%validate(err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      params%closure = "PSE3"
      params%correction = "PC+"
      params%finite_size = "hunenberger"
      params%tolerance = 1.0e-9_wp
      params%reuse_initial_guess = .false.
      params%cold_guess = "ts0"
      params%solver%jfnk%restart = 40
      params%solvent%npts = 2048
      call params%write_file(parameter_file, err)
      if (.not. allocated(err)) call copy%read_file(parameter_file, err)
      open (newunit=unit, file=parameter_file, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, copy%closure == "PSE3" .and. copy%correction == "PC+" .and. copy%susceptibility == "oz" &
         & .and. copy%finite_size == "hunenberger" &
         & .and. copy%tolerance == 1.0e-9_wp .and. .not. copy%reuse_initial_guess &
         & .and. copy%cold_guess == "ts0" &
         & .and. copy%solver%jfnk%restart == 40 .and. copy%solvent%npts == 2048, &
         & more="3D settings changed in the TOML round trip")
      if (allocated(error)) return

      copy%susceptibility = "linear"
      call copy%validate(err)
      call check(error, allocated(err), more="a finite-size correction with a non-OZ response was accepted")
      if (allocated(error)) return
      deallocate (err)
      copy%finite_size = "none"
      call copy%validate(err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      copy%finite_size = "unknown"
      call copy%validate(err)
      call check(error, allocated(err), more="an unknown finite-size correction was accepted")
      if (allocated(error)) return
      deallocate (err)
      call copy%init_defaults()
      call check(error, copy%finite_size == "none" .and. copy%cold_guess == "cs0", more="3D defaults were not restored")
      if (allocated(error)) return
      copy%cold_guess = "zero"
      call copy%validate(err)
      call check(error, allocated(err), more="an unknown cold guess was accepted")
      if (allocated(error)) return
      deallocate (err)
      copy%cold_guess = "TS0"
      call copy%validate(err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call copy%init_defaults()
      if (allocated(error)) return
      copy%susceptibility = "cubic"
      call copy%validate(err)
      call check(error, allocated(err), more="an unknown response mode was accepted")
      if (allocated(error)) return
      deallocate (err)
      call copy%init_defaults()
      copy%precond = "spectral"
      call copy%validate(err)
      call check(error, allocated(err), more="a preconditioner the 3D solve lacks was accepted")
      if (allocated(error)) return
      deallocate (err)
      call copy%init_defaults()
      copy%correction = "PC++"
      call copy%validate(err)
      call check(error, allocated(err), more="an unknown correction was accepted")
   end subroutine check_parameters

   !> Solved SPC/E water on a small VV grid
   !>
   !> @param[in]  ctx      run context
   !> @param[out] solvent  solved solvent
   !> @param[out] err      construction or solve failure
   subroutine solved_water(ctx, solvent, err)
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Solvent
      type(solvent_vv_type), intent(out) :: solvent
      !> Construction or solve failure
      type(moist_error), allocatable, intent(out) :: err
      type(solvent_vv_parameters_type) :: params
      params%npts = 1024
      params%dr = 0.06_wp
      params%tolerance = 1.0e-10_wp
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call solvent%solve(err, params)
   end subroutine solved_water

   !> Check a response table against the solvent response at each k-point
   !>
   !> - Every k-point's matrix is the response at its |k|, which its shell records
   !> - Matrices are symmetric for one solvent species
   !>
   !> @param[inout] error    test error
   !> @param[in]    grid     updated grid
   !> @param[in]    table    response table on the grid
   !> @param[in]    solvent  solved solvent
   !> @param[in]    label    grid label
   subroutine check_table(error, grid, table, solvent, label, oz)
      !> Test error
      type(error_type), allocatable, intent(inout) :: error
      !> Updated grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Response table
      type(uv3d_susceptibility_type), intent(in) :: table
      !> Solved solvent
      type(solvent_vv_type), intent(in) :: solvent
      !> Grid label
      character(len=*), intent(in) :: label
      logical, intent(in), optional :: oz
      type(moist_error), allocatable :: err
      type(packed_index_lut) :: lut
      real(wp), allocatable :: chi(:, :), k(:)
      integer :: j, g, a
      logical :: reconstruct

      call check(error, size(table%index) == grid%npts_k .and. table%nv == 3 .and. &
         & minval(table%index) >= 1 .and. maxval(table%index) == size(table%radius) .and. &
         & size(table%values, 3) == size(table%radius), more=label//": table shapes")
      if (allocated(error)) return
      allocate (k(grid%npts_k))
      do j = 1, grid%npts_k
         k(j) = norm2(grid%kpoint(j))
      end do
      call check(error, all(table%radius(table%index) == k), more=label//": shell radius differs from |k|")
      if (allocated(error)) return
      reconstruct = .false.
      if (present(oz)) reconstruct = oz
      if (reconstruct) then
         call solvent%susceptibility_oz_at(k, chi, err)
      else
         call solvent%susceptibility_at(k, chi, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) return
      call lut%init(3)
      do j = 1, grid%npts_k
         do a = 1, 3
            do g = 1, 3
               call check(error, table%values(g, a, table%index(j)), chi(j, lut%idx(g, a)), thr=1.0e-13_wp, &
                  & more=label//": response differs from susceptibility_at")
               if (allocated(error)) return
            end do
         end do
      end do
      call check(error, all(abs(table%values - reshape(table%values, shape(table%values), order=[2, 1, 3])) == 0.0_wp), &
         & more=label//": response matrices are not symmetric")
   end subroutine check_table

   !> Check the response table on a non-cubic Cartesian grid
   !>
   !> - direct: radial quadrature per distinct |k|, shared shells (nshell < npts_k)
   !> - linear: interpolation of the VV grid's response, close to direct but not equal
   subroutine check_susceptibility_cartesian(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(uv3d_susceptibility_type) :: direct, linear, oz
      real(wp) :: gap

      call new_context(ctx, nthreads=0, verbosity=0)
      call solved_water(ctx, solvent, err)
      if (.not. allocated(err)) call new_cartesian_point_grid(grid, err, 8, 10, 12, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call grid%update(solvent%structure(), err)
      if (.not. allocated(err)) call new_uv3d_susceptibility(direct, grid, solvent, "direct", err)
      if (.not. allocated(err)) call new_uv3d_susceptibility(linear, grid, solvent, "LINEAR", err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_table(error, grid, direct, solvent, "Cartesian")
      if (allocated(error)) return
      call check(error, size(direct%radius) < grid%npts_k, more="Cartesian k-points share no |k|")
      if (allocated(error)) return
      call check(error, all(linear%index == direct%index) .and. all(linear%radius == direct%radius), &
         & more="linear and direct tables differ in their shells")
      if (allocated(error)) return
      gap = maxval(abs(linear%values - direct%values))
      call check(error, gap > 0.0_wp .and. gap < 1.0e-2_wp, more="linear response not close to the direct one")
      if (allocated(error)) return
      call new_uv3d_susceptibility(oz, grid, solvent, "oz", err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_table(error, grid, oz, solvent, "Cartesian OZ", oz=.true.)
   end subroutine check_susceptibility_cartesian

   !> Check the response table on an atom-centred grid with its NUFFT k-points
   subroutine check_susceptibility_molecular(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid
      type(uv3d_susceptibility_type) :: table

      call new_context(ctx, nthreads=0, verbosity=0)
      call solved_water(ctx, solvent, err)
      if (.not. allocated(err)) call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(grid, err, recipe=recipe)
      if (.not. allocated(err)) call grid%update(solvent%structure(), err)
      if (.not. allocated(err)) call new_uv3d_susceptibility(table, grid, solvent, "direct", err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_table(error, grid, table, solvent, "molecular")
      if (allocated(error)) return
      call new_uv3d_susceptibility(table, grid, solvent, "oz", err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_table(error, grid, table, solvent, "molecular OZ", oz=.true.)
   end subroutine check_susceptibility_molecular

   !> Check refusals: no k-points, unknown mode, unsolved solvent, |k| beyond the VV band
   subroutine check_susceptibility_refusals(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent, unsolved
      type(moist_math_grid_3d_cartesian_type) :: grid, fine
      type(uv3d_susceptibility_type) :: table
      type(structure_type) :: mol, atom
      integer :: mode
      character(len=6), parameter :: modes(3) = [character(len=6) :: "oz", "direct", "linear"]

      call new_context(ctx, nthreads=0, verbosity=0)
      call solved_water(ctx, solvent, err)
      if (.not. allocated(err)) call new_cspce_water_solvent(ctx, unsolved, err)
      if (.not. allocated(err)) call new_cartesian_point_grid(grid, err, 8, 8, 8, 0.5_wp, margin=0.0_wp)
      ! pi/0.02 bohr per axis lies far beyond the band pi/0.06 bohr of the VV grid
      if (.not. allocated(err)) call new_cartesian_point_grid(fine, err, 4, 4, 4, 0.02_wp, margin=0.0_wp)
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = solvent%structure()
      ! A single atom fits the small fine box
      call new(atom, [10], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))

      call new_uv3d_susceptibility(table, grid, solvent, "direct", err)
      call check(error, allocated(err), more="a grid without k-points was accepted")
      if (allocated(error)) return
      deallocate (err)
      call grid%update(mol, err)
      if (.not. allocated(err)) call fine%update(atom, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new_uv3d_susceptibility(table, grid, solvent, "cubic", err)
      call check(error, allocated(err), more="an unknown mode was accepted")
      if (allocated(error)) return
      deallocate (err)
      do mode = 1, size(modes)
         call new_uv3d_susceptibility(table, grid, unsolved, modes(mode), err)
         call check(error, allocated(err), more=trim(modes(mode))//": an unsolved solvent was accepted")
         if (allocated(error)) return
         deallocate (err)
         call new_uv3d_susceptibility(table, fine, solvent, modes(mode), err)
         call check(error, allocated(err), more=trim(modes(mode))//": |k| beyond the VV band was accepted")
         if (allocated(error)) return
         deallocate (err)
      end do
   end subroutine check_susceptibility_refusals

   !> Synthetic 3D UV data on an updated grid
   !>
   !> - Response: the water table with nonsymmetric matrices, so a transposed
   !>   contraction shows
   !> - u_lr(k): smooth complex values decaying in |k|
   !> - c0, c1: smooth site fields around the first atom
   !>
   !> @param[in]  grid     updated grid
   !> @param[in]  solvent  solved solvent
   !> @param[out] table    response table
   !> @param[out] uk       long-range potential (npts_k, 3)
   !> @param[out] c0, c1   site fields (ngrid, 3)
   !> @param[out] err      table failure
   subroutine synthetic_data(grid, solvent, table, uk, c0, c1, err)
      !> Updated grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Solved solvent
      type(solvent_vv_type), intent(in) :: solvent
      !> Response table
      type(uv3d_susceptibility_type), intent(out) :: table
      !> Long-range potential
      complex(wp), allocatable, intent(out) :: uk(:, :)
      !> Site fields
      real(wp), allocatable, intent(out) :: c0(:, :), c1(:, :)
      !> Table failure
      type(moist_error), allocatable, intent(out) :: err
      real(wp) :: k(3), r(3), r0(3)
      integer :: j, g, a, s

      call new_uv3d_susceptibility(table, grid, solvent, "direct", err)
      if (allocated(err)) return
      do s = 1, size(table%values, 3)
         do a = 1, 3
            do g = 1, 3
               table%values(g, a, s) = table%values(g, a, s) + 0.05_wp*real(g - 2*a, wp)/real(s, wp)
            end do
         end do
      end do
      allocate (uk(grid%npts_k, 3), c0(grid%ngrid, 3), c1(grid%ngrid, 3))
      do j = 1, grid%npts_k
         k = grid%kpoint(j)
         do g = 1, 3
            uk(j, g) = cmplx(real(g, wp), 0.3_wp*k(1), wp)*exp(-0.5_wp*dot_product(k, k))
         end do
      end do
      r0 = [sum(grid%xyz(1, :)), sum(grid%xyz(2, :)), sum(grid%xyz(3, :))]/real(grid%ngrid, wp)
      do j = 1, grid%ngrid
         r = grid%xyz(:, j) - r0
         do g = 1, 3
            c0(j, g) = -real(g, wp)*exp(-0.5_wp*dot_product(r, r))
            c1(j, g) = c0(j, g) + 0.2_wp*sin(real(g, wp)*r(1))*exp(-0.3_wp*dot_product(r, r))
         end do
      end do
   end subroutine synthetic_data

   !> Check the 3D UV map and its derivative on Cartesian and atom-centred grids
   !>
   !> - The map equals a per-k evaluation written out here with its own engine
   !> - The map is affine: map(c1) - map(c0) = jvp(c1 - c0)
   !> - One and four threads agree
   subroutine check_equation_map(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_cartesian_type), target :: cart
      type(moist_math_grid_3d_molecular_type), target :: mole
      class(moist_math_grid_3d_type), pointer :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      type(uv3d_susceptibility_type) :: table
      type(uv3d_equation_type) :: eq1, eq4
      complex(wp), allocatable :: uk(:, :), ck(:, :), tk(:, :)
      real(wp), allocatable :: c0(:, :), c1(:, :), t0(:, :), t1(:, :), t4(:, :), dt(:, :), tref(:, :), work(:, :)
      real(wp), parameter :: beta = 1.3_wp
      complex(wp) :: c(3)
      real(wp) :: scale
      integer :: kind, j, a
      character(len=9) :: label

      call new_context(ctx, nthreads=0, verbosity=0)
      call solved_water(ctx, solvent, err)
      ! Exercise the threaded field copies as well as the smaller molecular case.
      if (.not. allocated(err)) call new_cartesian_point_grid(cart, err, 32, 32, 32, 0.4_wp, margin=0.0_wp)
      if (.not. allocated(err)) call cart%update(solvent%structure(), err)
      if (.not. allocated(err)) call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(mole, err, recipe=recipe)
      if (.not. allocated(err)) call mole%update(solvent%structure(), err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      do kind = 1, 2
         if (kind == 1) then
            grid => cart
            label = "Cartesian"
         else
            grid => mole
            label = "molecular"
         end if
         call synthetic_data(grid, solvent, table, uk, c0, c1, err)
         eq1%nthreads = 1
         eq4%nthreads = 4
         if (.not. allocated(err)) call eq1%setup(grid, table, uk, beta, err)
         if (.not. allocated(err)) call eq4%setup(grid, table, uk, beta, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         allocate (t0, t1, t4, dt, tref, work, mold=c0)
         work = c0
         call eq1%map(c0, t0, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, all(c0 == work), more=trim(label)//": map changed its input")
         if (allocated(error)) return
         work = c1
         if (.not. allocated(err)) call eq1%map(c1, t1, err)
         if (.not. allocated(err)) call eq4%map(c1, t4, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, all(c1 == work), more=trim(label)//": threaded map changed its input")
         if (allocated(error)) return
         work = c1 - c0
         call eq1%jvp(work, dt, err)
         if (.not. allocated(err)) call eq4%jvp(work, tref, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         scale = maxval(abs(t1))
         call check(error, all(work == c1 - c0), more=trim(label)//": Jacobian changed its input")
         if (allocated(error)) return
         call check(error, maxval(abs(tref - dt)) <= 1.0e-13_wp*scale, &
            & more=trim(label)//": Jacobian thread counts differ")
         if (allocated(error)) return
         call check(error, maxval(abs((t1 - t0) - dt)) <= 1.0e-12_wp*scale, &
            & more=trim(label)//": map differences differ from the derivative")
         if (allocated(error)) return
         call check(error, maxval(abs(t4 - t1)) <= 1.0e-13_wp*scale, more=trim(label)//": thread counts differ")
         if (allocated(error)) return

         ! Reference: t_s(k)_a = sum_g c_g P(g, a) - c_a - beta u_a with c = c_s(k) - beta u
         call grid%new_trafo(trafo, err)
         if (.not. allocated(err)) call trafo%prepare(3, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         allocate (ck(grid%npts_k, 3), tk(grid%npts_k, 3))
         work = c1
         call trafo%fft_r2k(work, ck, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         do j = 1, grid%npts_k
            c = ck(j, :) - beta*uk(j, :)
            do a = 1, 3
               tk(j, a) = c(1)*table%values(1, a, table%index(j)) + c(2)*table%values(2, a, table%index(j)) &
                  & + c(3)*table%values(3, a, table%index(j)) - c(a) - beta*uk(j, a)
            end do
         end do
         call trafo%fft_k2r(tk, tref, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, maxval(abs(t1 - tref)) <= 1.0e-12_wp*scale, &
            & more=trim(label)//": map differs from the per-k reference")
         if (allocated(error)) return
         work = c1
         call eq1%adjoint(c1, t0, err)
         if (.not. allocated(err)) call eq4%adjoint(c1, t4, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, all(c1 == work), more=trim(label)//": transpose changed its input")
         if (allocated(error)) return
         call check(error, maxval(abs(t4 - t0)) <= 1.0e-13_wp*maxval(abs(t0)), &
            & more=trim(label)//": transpose thread counts differ")
         if (allocated(error)) return
         call trafo%destroy()
         deallocate (trafo, ck, tk, t0, t1, t4, dt, tref, work)
      end do
   end subroutine check_equation_map

   !> Check that a copy rebuilds its own engine and outlives its source
   subroutine check_equation_copy(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid
      type(uv3d_susceptibility_type) :: table
      type(uv3d_equation_type), allocatable :: source
      type(uv3d_equation_type) :: copy
      complex(wp), allocatable :: uk(:, :)
      real(wp), allocatable :: c0(:, :), c1(:, :), t_source(:, :), t_copy(:, :)

      call new_context(ctx, nthreads=0, verbosity=0)
      call solved_water(ctx, solvent, err)
      if (.not. allocated(err)) call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(grid, err, recipe=recipe)
      if (.not. allocated(err)) call grid%update(solvent%structure(), err)
      if (.not. allocated(err)) call synthetic_data(grid, solvent, table, uk, c0, c1, err)
      allocate (source)
      if (.not. allocated(err)) call source%setup(grid, table, uk, 1.3_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (t_source, t_copy, mold=c1)
      call source%map(c1, t_source, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      copy = source
      deallocate (source)
      call check(error, copy%ready(), more="the copy is not set up")
      if (allocated(error)) return
      call copy%map(c1, t_copy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(t_copy == t_source), more="the copy maps differently from its source")
   end subroutine check_equation_copy

   !> Check refused setups and unset equations
   subroutine check_equation_setup_errors(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(moist_math_grid_3d_cartesian_type) :: grid, bare
      type(uv3d_susceptibility_type) :: table, empty
      type(uv3d_equation_type) :: eq
      complex(wp), allocatable :: uk(:, :)
      real(wp), allocatable :: c0(:, :), c1(:, :), t(:, :)

      call new_context(ctx, nthreads=0, verbosity=0)
      call solved_water(ctx, solvent, err)
      if (.not. allocated(err)) call new_cartesian_point_grid(grid, err, 8, 8, 8, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_cartesian_point_grid(bare, err, 8, 8, 8, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call grid%update(solvent%structure(), err)
      if (.not. allocated(err)) call synthetic_data(grid, solvent, table, uk, c0, c1, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (t, mold=c0)

      call eq%map(c0, t, err)
      call check(error, allocated(err), more="an unset equation mapped")
      if (allocated(error)) return
      deallocate (err)
      call eq%setup(bare, table, uk, 1.3_wp, err)
      call check(error, allocated(err) .and. .not. eq%ready(), more="a grid without k-points was accepted")
      if (allocated(error)) return
      deallocate (err)
      call eq%setup(grid, empty, uk, 1.3_wp, err)
      call check(error, allocated(err) .and. .not. eq%ready(), more="an empty response was accepted")
      if (allocated(error)) return
      deallocate (err)
      call eq%setup(grid, table, uk(:, :2), 1.3_wp, err)
      call check(error, allocated(err) .and. .not. eq%ready(), more="a mis-shaped u_lr(k) was accepted")
      if (allocated(error)) return
      deallocate (err)
      call eq%setup(grid, table, uk, 0.0_wp, err)
      call check(error, allocated(err) .and. .not. eq%ready(), more="a zero beta was accepted")
      if (allocated(error)) return
      deallocate (err)
      call eq%setup(grid, table, uk, 1.3_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call eq%map(c0(:, :2), t(:, :2), err)
      call check(error, allocated(err), more="mis-shaped fields were mapped")
   end subroutine check_equation_setup_errors

   !> Settings of the energy tests: a small VV grid and tight solves
   !>
   !> @param[out] params  settings
   subroutine energy_settings(params)
      !> Settings
      type(model_moz_3d_parameters_type), intent(out) :: params
      params%tolerance = 1.0e-8_wp
      params%solvent%npts = 1024
      params%solvent%dr = 0.06_wp
      params%solvent%tolerance = 1.0e-10_wp
   end subroutine energy_settings

   !> 3D model of a water solute in SPC/E water on the given grid, updated
   !>
   !> @param[in]  ctx      run context
   !> @param[in]  grid     grid template
   !> @param[in]  params   settings
   !> @param[out] model    updated model
   !> @param[out] err      construction or update failure
   subroutine new_water_model(ctx, grid, params, model, err)
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Grid template
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Settings
      type(model_moz_3d_parameters_type), intent(in) :: params
      !> Model
      type(model_moz_3d_type), intent(out) :: model
      !> Construction or update failure
      type(moist_error), allocatable, intent(out) :: err
      type(solvent_vv_type) :: solvent
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call new_moz_3d_model(model, ctx, grid, solvent, err, params)
      if (.not. allocated(err)) call model%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call model%potential%add(coulomb_spce, err)
      if (.not. allocated(err)) call model%update(solvent%structure(), err)
   end subroutine new_water_model

   !> Energy through a coupling of the model
   !>
   !> @param[inout] model   updated model
   !> @param[inout] energy  energy, Hartree
   !> @param[out]   err     coupling or energy failure
   subroutine model_energy(model, energy, err)
      !> Model
      class(solvation_model_type), intent(inout), target :: model
      !> Energy
      real(wp), intent(inout) :: energy
      !> Coupling or energy failure
      type(moist_error), allocatable, intent(out) :: err
      type(coupling_type), pointer :: coupling
      call model%new_coupling(coupling, err)
      if (allocated(err)) return
      call model%prepare_energy(coupling, err)
      if (.not. allocated(err)) call model%get_energy(coupling, energy, err)
      call model%release_coupling(coupling)
   end subroutine model_energy

   !> Cartesian grid of the energy tests: 32^3 points, 0.6 bohr, 8 bohr margin
   !>
   !> @param[out] grid  grid template
   !> @param[out] err   configuration failure
   subroutine energy_grid(grid, err)
      !> Grid template
      type(moist_math_grid_3d_cartesian_type), intent(out) :: grid
      !> Configuration failure
      type(moist_error), allocatable, intent(out) :: err
      call new_cartesian_point_grid(grid, err, 32, 32, 32, 0.6_wp, margin=8.0_wp)
   end subroutine energy_grid

   !> Check the energy of a water solute on a Cartesian grid
   !>
   !> - The solve converges; a second energy with the same tables reuses it
   !> - PC+ adds its term; the partial molar volume is positive
   !> - Correlations on the grid with 1 + h >= 0
   subroutine check_energy_cartesian(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(model_moz_3d_parameters_type) :: params
      type(model_moz_3d_type), target :: plain, corrected
      type(model_moz_3d_thermodynamics_type) :: thermo
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :)
      real(wp) :: e0, e1, e_plus

      call new_context(ctx, nthreads=0, verbosity=0)
      call energy_settings(params)
      call energy_grid(grid, err)
      if (.not. allocated(err)) call new_water_model(ctx, grid, params, plain, err)
      params%correction = "PC+"
      if (.not. allocated(err)) call new_water_model(ctx, grid, params, corrected, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      e0 = 0.0_wp
      e1 = 0.0_wp
      e_plus = 0.0_wp
      call model_energy(plain, e0, err)
      if (.not. allocated(err)) call model_energy(plain, e1, err)
      if (.not. allocated(err)) call model_energy(corrected, e_plus, err)
      if (.not. allocated(err)) call plain%thermodynamics(thermo, err)
      if (.not. allocated(err)) call plain%correlations(cr, tr, hr, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, e1 == e0, more="a repeated energy with the same tables changed")
      if (allocated(error)) return
      call check(error, e0, thermo%free_energy, thr=1.0e-14_wp, more="no correction changed the energy")
      if (allocated(error)) return
      call check(error, e_plus - e0, thermo%pc_plus, rel=.true., thr=1.0e-8_wp, more="PC+ energy lacks its term")
      if (allocated(error)) return
      call check(error, thermo%partial_molar_volume > 0.0_wp, more="non-positive partial molar volume")
      if (allocated(error)) return
      call check(error, all(shape(hr) == [grid%nx*grid%ny*grid%nz, 3]) .and. minval(hr) >= -1.0_wp - 1.0e-10_wp, &
         & more="correlations of the wrong shape or with 1 + h < 0")
   end subroutine check_energy_cartesian

   !> Finite-size energy uses the last successful Cartesian grid update
   !>
   !> A rectangular box tests all three lengths. Editing dr before evaluation
   !> keeps the built box; updating then doubles its sides and halves the term.
   !> Molecular grids remain unsupported by this box-only correction.
   subroutine check_finite_size_box(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(moist_math_grid_3d_molecular_type) :: molecular
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(model_moz_3d_parameters_type) :: params
      type(model_moz_3d_type), target :: model, refused
      type(model_moz_3d_thermodynamics_type) :: thermo
      type(solvent_vv_type) :: solvent
      type(custom_lj_type) :: lj
      type(fixed_charges_type) :: charges
      type(structure_type) :: ion
      real(wp) :: energy, screening, expected, first

      call new_context(ctx, nthreads=0, verbosity=0)
      call energy_settings(params)
      params%finite_size = "hunenberger"
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call solvent%solve(err, params%solvent)
      if (.not. allocated(err)) call solvent%charge_screening(screening, err)
      if (.not. allocated(err)) call new_cartesian_point_grid(grid, err, 24, 28, 32, 0.8_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_moz_3d_model(model, ctx, grid, solvent, err, params)
      if (.not. allocated(err)) call new_custom_lj(lj, [3.5_wp*aatoau], [0.1_wp/autokcal], err)
      if (.not. allocated(err)) call new_fixed_charges(charges, [-1.0_wp], err)
      if (.not. allocated(err)) call model%potential%add(lj, err)
      if (.not. allocated(err)) call model%potential%add(charges, err)
      call new(ion, [17], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      ion%charge = -1.0_wp
      if (.not. allocated(err)) call model%update(ion, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      select type (built => model%grid)
      type is (moist_math_grid_3d_cartesian_type)
         built%dr = 1.6_wp
      end select
      energy = 0.0_wp
      call model_energy(model, energy, err)
      if (.not. allocated(err)) call model%thermodynamics(thermo, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      expected = electrostatic_finite_size_correction(-1.0_wp, screening, &
         & coulomb_box_mean([19.2_wp, 22.4_wp, 25.6_wp]))
      call check(error, thermo%finite_size, expected, thr=1.0e-14_wp, &
         & more="finite-size correction used pending settings instead of the built rectangular box")
      if (allocated(error)) return
      call check(error, energy, thermo%free_energy + expected, thr=1.0e-12_wp)
      if (allocated(error)) return
      first = thermo%finite_size
      call model%update(ion, err)
      if (.not. allocated(err)) call model_energy(model, energy, err)
      if (.not. allocated(err)) call model%thermodynamics(thermo, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, thermo%finite_size, 0.5_wp*first, thr=1.0e-14_wp, &
         & more="updating doubled box lengths did not refresh the finite-size correction")
      if (allocated(error)) return

      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(molecular, err, recipe=recipe)
      if (.not. allocated(err)) call new_moz_3d_model(refused, ctx, molecular, solvent, err, params)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call refused%update(ion, err)
      call check(error, allocated(err), more="molecular quadrature was accepted for a Cartesian-box correction")
      if (allocated(error)) return
      call check(error, index(err%message, "requires a Cartesian grid") > 0, more=err%message)
   end subroutine check_finite_size_box

   !> Check that one and four threads give the same energy
   subroutine check_energy_threads(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx1, ctx4
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(model_moz_3d_parameters_type) :: params
      type(model_moz_3d_type), target :: one, four
      real(wp) :: e1, e4

      call new_context(ctx1, nthreads=1, verbosity=0)
      call new_context(ctx4, nthreads=4, verbosity=0)
      call energy_settings(params)
      call energy_grid(grid, err)
      if (.not. allocated(err)) call new_water_model(ctx1, grid, params, one, err)
      if (.not. allocated(err)) call new_water_model(ctx4, grid, params, four, err)
      e1 = 0.0_wp
      e4 = 0.0_wp
      if (.not. allocated(err)) call model_energy(one, e1, err)
      if (.not. allocated(err)) call model_energy(four, e4, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, e4, e1, rel=.true., thr=1.0e-10_wp, more="thread counts give different energies")
   end subroutine check_energy_threads

   !> Check refused energies and stale correlations
   !>
   !> - No correlations or thermodynamics before an energy
   !> - A solve that does not converge is an error and leaves the energy unchanged
   !> - An update drops the correlations
   subroutine check_energy_refusals(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(model_moz_3d_parameters_type) :: params
      type(model_moz_3d_type), target :: model, short
      type(model_moz_3d_thermodynamics_type) :: thermo
      type(solvent_vv_type) :: solvent
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :)
      real(wp) :: energy

      call new_context(ctx, nthreads=0, verbosity=0)
      call energy_settings(params)
      call energy_grid(grid, err)
      if (.not. allocated(err)) call new_water_model(ctx, grid, params, model, err)
      params%maxiter = 1
      if (.not. allocated(err)) call new_water_model(ctx, grid, params, short, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      call model%correlations(cr, tr, hr, err)
      call check(error, allocated(err), more="correlations before an energy")
      if (allocated(error)) return
      deallocate (err)
      call model%thermodynamics(thermo, err)
      call check(error, allocated(err), more="thermodynamics before an energy")
      if (allocated(error)) return
      deallocate (err)

      energy = 5.0_wp
      call model_energy(short, energy, err)
      call check(error, allocated(err) .and. energy == 5.0_wp, more="an unconverged solve changed the energy")
      if (allocated(error)) return
      call check(error, index(err%message, "did not converge") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      energy = 0.0_wp
      call model_energy(model, energy, err)
      if (.not. allocated(err)) call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call model%update(solvent%structure(), err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call model%correlations(cr, tr, hr, err)
      call check(error, allocated(err), more="an update kept the correlations")
   end subroutine check_energy_refusals

end module test_moz_3d_uv
