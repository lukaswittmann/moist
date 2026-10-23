!> Tests of the 1D MOZ (UV) settings, solvent preparation, equation and solve: moz_1d_uv
!>
!> SPC/E water as the solvent; the equation tests use a four-site solute with a
!> synthetic nonsymmetric solvent response; each test runs as a selected test
module test_moz_1d_uv
   use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_associated
   !$ use omp_lib, only: omp_in_parallel
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check, skip_test, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair, &
      & transform_quadrature
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo
   use moist_math_packed_sym, only: packed_index_lut, npair_from_ns
   use moist_model_moz_1d_kernel, only: uv_kblock, uv_block_product, uv_block_product_transpose, &
      & uv_charge_moment_correction
   use moist_model_moz_1d_equation, only: uv_equation_type
   use moist_model_moz_1d_driver, only: uv_problem_type, new_uv_problem
   use moist_model_moz_1d_operators, only: uv_operators_type
   use moist_model_moz_solvent_vv_kernel, only: vv_intramolecular
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_solvent_vv, only: solvent_vv_type, solvent_vv_thermodynamics_type, new_vv_solvent
   use moist_model_moz, only: model_moz_1d_thermodynamics_type
   use moist_model_moz_potential, only: lj_cspce, coulomb_spce, monopole_type, potential_sites_type, &
      & custom_lj_type, new_custom_lj, fixed_charges_type, new_fixed_charges, moz_potential_type
   use moist_solvent_system, only: solvent_system_type
   use mctc_io_convert, only: aatoau, autokcal
   use mctc_io_codata2018, only: Boltzmann_constant, atomic_unit_of_energy
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type, response_item_type, atomic_charge_adjoint_response_type
   use mctc_io, only: structure_type, new
   use test_helpers, only: check_moist_error
   use test_moz_fixtures, only: new_cspce_water_solvent
   implicit none(type, external)
   private
   public :: collect_moz_1d_uv

   !> Parameter file of the round-trip test
   character(len=*), parameter :: parameter_file = "moist-test-moz-1d-parameters.toml"
   !> Solute and solvent sites, grid points and spacing of the equation tests
   integer, parameter :: nu = 4, nv = 3, npts = 512
   real(wp), parameter :: dr = 0.05_wp
   !> Inverse thermal energy of the equation tests, 1/Hartree
   real(wp), parameter :: beta = 1.7_wp

   !> Synthetic UV data on the test grid
   type :: uv_case_type
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(packed_index_lut) :: lut
      !> Solute W(k), packed (npts, nu(nu+1)/2), and solvent response (npts, nv, nv)
      real(wp), allocatable :: wu(:, :), response(:, :, :)
      !> Long-range potential and model c_s(r), (npts, nu nv)
      real(wp), allocatable :: uk_lr(:, :), cr(:, :)
   end type uv_case_type

contains

   !> Register tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_1d_uv(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("parameters", check_parameters), &
         & new_unittest("prepare_solvent", check_prepare_solvent), &
         & new_unittest("block_products", check_block_products), &
         & new_unittest("map_consistency", check_map_consistency), &
         & new_unittest("jvp_affine", check_jvp_affine), &
         & new_unittest("thread_independence", check_thread_independence), &
         & new_unittest("setup_errors", check_setup_errors), &
         & new_unittest("problem_setup_errors", check_problem_setup_errors), &
         & new_unittest("warm_start", check_warm_start), &
         & new_unittest("scf_reuse", check_scf_reuse), &
         & new_unittest("operator_cache", check_operator_cache), &
         & new_unittest("charge_moment", check_charge_moment), &
         & new_unittest("water_in_water", check_water_in_water), &
         & new_unittest("response_tables", check_response_tables), &
         & new_unittest("response_charges", check_response_charges), &
         & new_unittest("response_fixed", check_response_fixed), &
         & new_unittest("host_charges", check_host_charges), &
         & new_unittest("energy_refusals", check_energy_refusals), &
         & new_unittest("pysofft_reference", check_pysofft_reference)]
   end subroutine collect_moz_1d_uv

   !> Check the 1D settings
   !>
   !> - Defaults: KH, JFNK with the spectral filter; VV defaults
   !>   in the solvent group
   !> - TOML round trip, including the nested solver and solvent groups
   !> - Names the UV solve does not implement and invalid nested settings are refused
   subroutine check_parameters(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_1d_parameters_type) :: params, copy
      integer :: unit, stat

      call check(error, params%closure == "KH" .and. params%solver%name == "JFNK" .and. &
         & params%precond == "spectral" .and. params%solvent%closure == "KH" .and. &
         & params%solvent%solver%jfnk%deflation, more="unexpected 1D defaults")
      if (allocated(error)) return
      call params%validate(err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      params%closure = "PSE3"
      params%tolerance = 1.0e-9_wp
      params%precond = "none"
      params%solver%jfnk%restart = 40
      params%solver%jfnk%adjoint%restart = 12
      params%solver%jfnk%adjoint%max_krylov = 36
      params%reuse_initial_guess = .false.
      params%solvent%closure = "HNC"
      params%solvent%npts = 2048
      params%solvent%solver%jfnk%restart = 30
      call params%write_file(parameter_file, err)
      if (.not. allocated(err)) call copy%read_file(parameter_file, err)
      open (newunit=unit, file=parameter_file, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, copy%closure == "PSE3" .and. &
         & copy%tolerance == 1.0e-9_wp .and. copy%precond == "none" .and. &
         & .not. copy%reuse_initial_guess .and. copy%solver%jfnk%adjoint%restart == 12 .and. &
         & copy%solver%jfnk%adjoint%max_krylov == 36 .and. &
         & copy%solver%jfnk%restart == 40 .and. copy%solvent%closure == "HNC" .and. copy%solvent%npts == 2048 .and. &
         & copy%solvent%solver%jfnk%restart == 30, more="1D settings changed in the TOML round trip")
      if (allocated(error)) return

      copy%precond = "kspace"
      call copy%validate(err)
      call check(error, allocated(err), more="the VV-only kspace preconditioner was accepted for UV")
      if (allocated(error)) return
      deallocate (err)
      call copy%init_defaults()
      copy%solver%jfnk%deflation = .true.
      call copy%validate(err)
      call check(error, allocated(err), more="the VV-only deflation was accepted for UV")
      if (allocated(error)) return
      deallocate (err)
      call copy%init_defaults()
      copy%solvent%nlambda = 0
      call copy%validate(err)
      call check(error, allocated(err), more="invalid nested solvent settings were accepted")
   end subroutine check_parameters

   !> Check the lazy VV solve of the solvent copy
   !>
   !> - Construction keeps an unsolved solvent unsolved
   !> - prepare_solvent solves the copy with the solvent settings; the caller's
   !>   solvent stays unsolved
   !> - A compatible solved solvent is reused, and invalid settings are refused
   subroutine check_prepare_solvent(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_1d_type) :: model, model2
      type(model_moz_1d_parameters_type) :: params
      type(moist_math_grid_radial_type) :: rgrid, kgrid

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 512, 0.05_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      params%solvent%npts = 1024
      params%solvent%dr = 0.06_wp

      call new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, err, params)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. model%solvent_solved(), more="construction solved the solvent")
      if (allocated(error)) return
      call model%prepare_solvent(err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model%solvent_solved() .and. .not. solvent%is_solved(), &
         & more="prepare_solvent did not solve the copy alone")
      if (allocated(error)) return

      call solvent%solve(err, params%solvent)
      if (.not. allocated(err)) call new_moz_1d_model(model2, ctx, rgrid, kgrid, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model2%solvent_solved(), more="a solved solvent lost its solution in the copy")
      if (allocated(error)) return
      call model2%prepare_solvent(err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      params%closure = "PY"
      call new_moz_1d_model(model2, ctx, rgrid, kgrid, solvent, err, params)
      call check(error, allocated(err), more="invalid settings were accepted")
   end subroutine check_prepare_solvent

   !> Build the synthetic UV case
   !>
   !> - Tetrahedral four-site solute, C-H 2.06 bohr
   !> - P_ba(k) = delta_ab + 0.2 (1 + 0.1 a + 0.05 b) exp(-k^2/(1 + 0.1 (a + b))), not symmetric
   !>
   !> @param[out] case  UV data
   !> @param[out] err   grid failure
   !> @param[in]  grid_npts optional radial grid size
   subroutine new_uv_case(case, err, grid_npts)
      !> UV data
      type(uv_case_type), intent(out) :: case
      !> Grid failure
      type(moist_error), allocatable, intent(out) :: err
      !> Optional radial grid size
      integer, intent(in), optional :: grid_npts
      real(wp), parameter :: d = 2.06_wp/sqrt(3.0_wp)
      real(wp), parameter :: xyz(3, nu) = reshape([0.0_wp, 0.0_wp, 0.0_wp, d, d, d, -d, -d, d, &
         & -d, d, -d, d, -d, -d], [3, nu])
      integer :: a, b, col, np

      np = npts
      if (present(grid_npts)) np = grid_npts
      call new_uniform_radial_pair(case%rgrid, case%kgrid, np, dr, err)
      if (allocated(err)) return
      call case%lut%init(nu)
      allocate (case%wu(np, npair_from_ns(nu)), case%response(np, nv, nv), case%uk_lr(np, nu*nv), &
         & case%cr(np, nu*nv))
      call vv_intramolecular(case%kgrid%r, xyz, case%lut, case%wu)
      do a = 1, nv
         do b = 1, nv
            case%response(:, b, a) = 0.2_wp*(1.0_wp + 0.1_wp*a + 0.05_wp*b)* &
               & exp(-case%kgrid%r**2/(1.0_wp + 0.1_wp*(a + b)))
            if (a == b) case%response(:, b, a) = case%response(:, b, a) + 1.0_wp
         end do
      end do
      do col = 1, nu*nv
         case%uk_lr(:, col) = 0.01_wp*exp(-case%kgrid%r**2)/real(col, wp)
         case%cr(:, col) = -(0.3_wp + 0.05_wp*col)*exp(-(case%rgrid%r/2.5_wp)**2)
      end do
   end subroutine new_uv_case

   !> Set up an equation from the case
   !>
   !> @param[in]  case      UV data
   !> @param[out] eq        equation
   !> @param[in]  nthreads  team size
   !> @param[out] err       setup failure
   subroutine setup_uv_equation(case, eq, nthreads, err)
      !> UV data
      type(uv_case_type), intent(in) :: case
      !> Equation
      type(uv_equation_type), intent(out) :: eq
      !> Team size
      integer, intent(in) :: nthreads
      !> Setup failure
      type(moist_error), allocatable, intent(out) :: err
      eq%nthreads = nthreads
      call eq%setup(case%rgrid, case%kgrid, case%wu, case%response, beta, case%uk_lr, err)
   end subroutine setup_uv_equation

   !> Deterministic smooth UV field, (npts, nu nv)
   !>
   !> @param[in]  r     radial nodes, bohr
   !> @param[in]  seed  pattern selector
   !> @param[out] v     field
   subroutine field(r, seed, v)
      !> Radial nodes
      real(wp), intent(in) :: r(:)
      !> Pattern selector
      integer, intent(in) :: seed
      !> Field
      real(wp), intent(out) :: v(:, :)
      integer :: col
      do col = 1, size(v, 2)
         v(:, col) = sin(0.6_wp*(col + seed)*r + 0.2_wp*seed)*exp(-r/3.0_wp)
      end do
   end subroutine field

   !> Full solute W of one lane from the packed block
   pure function full_w(lut, w) result(wm)
      type(packed_index_lut), intent(in) :: lut
      real(wp), intent(in) :: w(:)
      real(wp) :: wm(lut%ns, lut%ns)
      integer :: i, j
      do j = 1, lut%ns
         do i = 1, lut%ns
            wm(i, j) = w(lut%idx(i, j))
         end do
      end do
   end function full_w

   !> Check the block products against dense matrix products per lane
   !>
   !> - h = W c P with c and h as nu x nv matrices
   !> - <u, W c P> = <W^T u P^T, c>
   subroutine check_block_products(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(uv_case_type) :: case
      real(wp) :: c(uv_kblock, nu*nv), u(uv_kblock, nu*nv), h(uv_kblock, nu*nv), v(uv_kblock, nu*nv)
      real(wp) :: cm(nu, nv), hm(nu, nv)
      integer :: kk, i, a

      call new_uv_case(case, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do kk = 1, uv_kblock
         do i = 1, nu*nv
            c(kk, i) = cos(0.3_wp*kk*i)
            u(kk, i) = sin(0.7_wp*kk + 0.4_wp*i)
         end do
      end do
      call uv_block_product(case%lut, nv, case%wu(1:uv_kblock, :), c, case%response(1:uv_kblock, :, :), h)
      call uv_block_product_transpose(case%lut, nv, case%wu(1:uv_kblock, :), u, case%response(1:uv_kblock, :, :), v)
      do kk = 1, uv_kblock
         do i = 1, nu
            do a = 1, nv
               cm(i, a) = c(kk, (i - 1)*nv + a)
            end do
         end do
         hm = matmul(full_w(case%lut, case%wu(kk, :)), matmul(cm, case%response(kk, :, :)))
         do i = 1, nu
            do a = 1, nv
               call check(error, h(kk, (i - 1)*nv + a), hm(i, a), thr=1.0e-13_wp, more="block product differs")
               if (allocated(error)) return
            end do
         end do
      end do
      call check(error, abs(sum(u*h) - sum(v*c)) < 1.0e-12_wp*sqrt(sum(u**2)*sum(h**2)), &
         & more="transpose block product fails the dot-product test")
   end subroutine check_block_products

   !> Check the map against point products and the Ng relation
   !>
   !> - h(k) = W c(k) P at selected points, c(k) = r2k(c_s) - beta u_lr(k)
   !> - t_s(r) = k2r(h - c - beta u_lr)
   subroutine check_map_consistency(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(uv_case_type) :: case
      type(uv_equation_type) :: eq
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      real(wp), allocatable :: tr(:, :), hk(:, :), ck(:, :), tk(:, :), tref(:, :)
      real(wp) :: cm(nu, nv), hm(nu, nv)
      integer, parameter :: points(4) = [1, 2, 100, npts]
      integer :: col, ip, ik, i, a

      call new_uv_case(case, err)
      if (.not. allocated(err)) call setup_uv_equation(case, eq, 2, err)
      if (.not. allocated(err)) call new_radial_trafo(trafo, case%rgrid, case%kgrid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (tr, hk, ck, tk, tref, mold=case%cr)
      call eq%map(case%cr, tr, err, hk=hk)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do col = 1, nu*nv
         call trafo%fbt_r2k(case%cr(:, col), ck(:, col), err)
         if (allocated(err)) exit
         ck(:, col) = ck(:, col) - beta*case%uk_lr(:, col)
         tk(:, col) = hk(:, col) - ck(:, col) - beta*case%uk_lr(:, col)
         call trafo%fbt_k2r(tk(:, col), tref(:, col), err)
         if (allocated(err)) exit
      end do
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(tr - tref)), 0.0_wp, thr=1.0e-12_wp*maxval(abs(tref)), &
         & more="t_s(r) is not k2r(h - c - beta u_lr)")
      if (allocated(error)) return
      do ip = 1, size(points)
         ik = points(ip)
         do i = 1, nu
            do a = 1, nv
               cm(i, a) = ck(ik, (i - 1)*nv + a)
            end do
         end do
         hm = matmul(full_w(case%lut, case%wu(ik, :)), matmul(cm, case%response(ik, :, :)))
         do i = 1, nu
            do a = 1, nv
               call check(error, hk(ik, (i - 1)*nv + a), hm(i, a), thr=1.0e-12_wp*maxval(abs(hm)), &
                  & more="map h(k) differs from the point product")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine check_map_consistency

   !> Check the derivative and its transpose
   !>
   !> - The map is affine: map(c + v) - map(c) = jvp(v)
   !> - <u, jvp(v)> = <adjoint(u), v>
   subroutine check_jvp_affine(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(uv_case_type) :: case
      type(uv_equation_type) :: eq
      real(wp), allocatable :: v(:, :), u(:, :), t0(:, :), t1(:, :), jv(:, :), ju(:, :)

      call new_uv_case(case, err)
      if (.not. allocated(err)) call setup_uv_equation(case, eq, 2, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (v, u, t0, t1, jv, ju, mold=case%cr)
      call field(case%rgrid%r, 1, v)
      call field(case%rgrid%r, 5, u)
      call eq%map(case%cr, t0, err)
      if (.not. allocated(err)) call eq%map(case%cr + v, t1, err)
      if (.not. allocated(err)) call eq%jvp(v, jv, err)
      if (.not. allocated(err)) call eq%adjoint(u, ju, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(t1 - t0 - jv)), 0.0_wp, thr=1.0e-12_wp*maxval(abs(jv)), &
         & more="the derivative differs from the difference of the affine map")
      if (allocated(error)) return
      call check(error, abs(sum(u*jv) - sum(ju*v)) < 1.0e-12_wp*sqrt(sum(u**2)*sum(jv**2)), &
         & more="transpose product fails the dot-product test")
   end subroutine check_jvp_affine

   !> Check that map, derivative and transpose do not depend on the team size
   subroutine check_thread_independence(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(uv_case_type) :: case
      type(uv_equation_type) :: eq
      real(wp), allocatable :: v(:, :), t1(:, :), t4(:, :), j1(:, :), j4(:, :), a1(:, :), a4(:, :)
      logical :: in_parallel
      integer :: igrid

      in_parallel = .false.
      !$ in_parallel = omp_in_parallel()
      if (in_parallel) then
         call skip_test(error, "Thread test must run as a selected test outside the suite OpenMP team")
         return
      end if
      do igrid = 1, 2
         call new_uv_case(case, err, merge(npts, 8192, igrid == 1))
         call check_moist_error(error, err)
         if (allocated(error)) return
         allocate (v, t1, t4, j1, j4, a1, a4, mold=case%cr)
         call field(case%rgrid%r, 2, v)
         call setup_uv_equation(case, eq, 1, err)
         if (.not. allocated(err)) call eq%map(case%cr, t1, err)
         if (.not. allocated(err)) call eq%jvp(v, j1, err)
         if (.not. allocated(err)) call eq%adjoint(v, a1, err)
         if (.not. allocated(err)) call setup_uv_equation(case, eq, 4, err)
         if (.not. allocated(err)) call eq%map(case%cr, t4, err)
         if (.not. allocated(err)) call eq%jvp(v, j4, err)
         if (.not. allocated(err)) call eq%adjoint(v, a4, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, all(t1 == t4) .and. all(j1 == j4) .and. all(a1 == a4), &
            & more="results depend on the team size")
         if (allocated(error)) return
         deallocate (v, t1, t4, j1, j4, a1, a4)
      end do
   end subroutine check_thread_independence

   !> Check refusals: map before setup, wrong shapes, non-DST-IV pair
   subroutine check_setup_errors(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(uv_case_type) :: case
      type(uv_equation_type) :: eq
      real(wp), allocatable :: tr(:, :)

      call new_uv_case(case, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (tr, mold=case%cr)
      call eq%map(case%cr, tr, err)
      call check(error, allocated(err), more="map before setup succeeded")
      if (allocated(error)) return
      deallocate (err)
      call eq%setup(case%rgrid, case%kgrid, case%wu(:, 1:9), case%response, beta, case%uk_lr, err)
      call check(error, allocated(err) .and. .not. eq%ready(), more="a solute W of no triangle size was accepted")
      if (allocated(error)) return
      deallocate (err)
      call eq%setup(case%rgrid, case%kgrid, case%wu, case%response, beta, case%uk_lr(:, 1:nu*nv - 1), err)
      call check(error, allocated(err), more="a long-range table of the wrong width was accepted")
      if (allocated(error)) return
      deallocate (err)
      case%rgrid%transform = transform_quadrature
      call setup_uv_equation(case, eq, 1, err)
      call check(error, allocated(err), more="a quadrature grid was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "DST-IV") > 0, more=err%message)
   end subroutine check_setup_errors

   !> Refuse a density count different from the solvent response, even when it divides the field width
   !>
   !> @param[out] error Test failure
   subroutine check_problem_setup_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(uv_case_type) :: case
      type(uv_problem_type) :: problem
      type(model_moz_1d_parameters_type) :: params, defaults
      real(wp), allocatable :: rho(:), u_sr(:, :), ur_lr(:, :)
      integer, parameter :: wrong_counts(4) = [0, 2, 4, 6]
      !> VV-only and unknown names, refused by the UV problem itself
      character(len=8), parameter :: bad_precond(2) = [character(len=8) :: "kspace", "bogus"]
      !> Quiet four-thread run context borrowed by the problems
      type(moist_context_type), target :: ctx
      integer :: i

      call new_context(ctx, nthreads=4, verbosity=0)
      call new_uv_case(case, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      u_sr = 0.0_wp*case%cr
      ur_lr = u_sr
      do i = 1, size(wrong_counts)
         rho = spread(0.005_wp, 1, wrong_counts(i))
         call new_uv_problem(problem, ctx, case%rgrid, case%kgrid, case%wu, case%response, u_sr, ur_lr, &
            & case%uk_lr, rho, beta, params, err)
         call check(error, allocated(err), more="a density count different from the solvent response was accepted")
         if (allocated(error)) return
         call check(error, .not. problem%equation%ready() .and. .not. allocated(problem%rho), &
            & more="a refused density count left a usable problem")
         if (allocated(error)) return
      end do
      rho = [0.005_wp, 0.004_wp, 0.003_wp]
      do i = 1, size(bad_precond)
         params%precond = bad_precond(i)
         call new_uv_problem(problem, ctx, case%rgrid, case%kgrid, case%wu, case%response, u_sr, ur_lr, &
            & case%uk_lr, rho, beta, params, err)
         call check(error, allocated(err), more="1D UV accepted preconditioner '"//trim(bad_precond(i))//"'")
         if (allocated(error)) return
      end do
      params%precond = defaults%precond
      params%damp_small_k = .true.
      call new_uv_problem(problem, ctx, case%rgrid, case%kgrid, case%wu, case%response, u_sr, ur_lr, &
         & case%uk_lr, rho, beta, params, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, problem%equation%ready() .and. all(problem%rho == rho), &
         & more="a valid site density vector was refused or changed")
      if (allocated(error)) return
      call check(error, problem%spectral%nthreads == 4 .and. problem%damping%nthreads == 4, &
         & more="the filters did not inherit the model's thread count")
      call ctx%delete()
   end subroutine check_problem_setup_errors

   !> Warm and fresh solves agree after a potential change; unusable guesses fall back
   subroutine check_warm_start(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(uv_case_type) :: case
      type(uv_problem_type) :: cold, warm
      type(model_moz_1d_parameters_type) :: params
      !> Quiet run context borrowed by the problems
      type(moist_context_type), target :: ctx
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :), guess(:, :), invalid(:, :)
      real(wp) :: rmse, e_cold, e_warm
      integer :: niter, i
      logical :: converged

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_uv_case(case, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      params%tolerance = 1.0e-12_wp
      params%damp_small_k = .false.
      u_sr = -0.02_wp*case%cr
      ur_lr = 0.0_wp*case%cr
      call build_problem(cold)
      if (.not. allocated(err)) call cold%solve(params%tolerance, params%maxiter, converged, rmse, niter, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged, more="fresh UV solve did not converge")
      if (allocated(error)) return
      guess = cold%tr

      ! A converged field must be accepted as the initial iterate, with a fresh solver stage.
      call build_problem(warm)
      if (.not. allocated(err)) call warm%solve(params%tolerance, 1, converged, rmse, niter, err, guess=guess)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == 1, more="the converged starting guess was not used")
      if (allocated(error)) return

      ! A previous solution is only a guess: both solves must converge for the changed potential.
      u_sr = 1.01_wp*u_sr
      call build_problem(cold)
      if (.not. allocated(err)) call cold%solve(params%tolerance, params%maxiter, converged, rmse, niter, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged, more="fresh solve of the changed potential did not converge")
      if (allocated(error)) return
      call cold%free_energy(e_cold, err)
      if (.not. allocated(err)) call build_problem(warm)
      if (.not. allocated(err)) then
         call warm%solve(params%tolerance, params%maxiter, converged, rmse, niter, err, guess=guess)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged, more="warm solve of the changed potential did not converge")
      if (allocated(error)) return
      call warm%free_energy(e_warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(warm%tr - cold%tr)), 0.0_wp, thr=1.0e-10_wp, &
         & more="warm and fresh correlations differ")
      if (allocated(error)) return
      call check(error, e_warm, e_cold, thr=1.0e-10_wp, more="warm and fresh energies differ")
      if (allocated(error)) return

      do i = 1, 2
         if (i == 1) then
            invalid = guess(1:npts - 1, :)
         else
            invalid = guess
            invalid(1, 1) = ieee_value(0.0_wp, ieee_quiet_nan)
         end if
         call build_problem(warm)
         if (.not. allocated(err)) then
            call warm%solve(params%tolerance, params%maxiter, converged, rmse, niter, err, guess=invalid)
         end if
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged, more="an unusable guess prevented convergence")
         if (allocated(error)) return
         call check(error, maxval(abs(warm%tr - cold%tr)), 0.0_wp, thr=1.0e-10_wp, &
            & more="an unusable guess changed the fresh solve")
         if (allocated(error)) return
      end do
      call ctx%delete()

   contains

      !> Rebuild a problem with the current tables, without carrying solver history
      subroutine build_problem(problem)
         type(uv_problem_type), intent(out) :: problem
         call new_uv_problem(problem, ctx, case%rgrid, case%kgrid, case%wu, case%response, u_sr, ur_lr, &
            & case%uk_lr, spread(0.005_wp, 1, nv), beta, params, err)
      end subroutine build_problem

   end subroutine check_warm_start

   !> Changed host tables retain storage, preserve the split, and recover from bad warm solves
   subroutine check_scf_reuse(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(uv_case_type) :: case
      type(uv_problem_type), target :: warm, cold
      type(model_moz_1d_parameters_type) :: params
      type(c_ptr) :: basis
      real(wp), allocatable :: sr(:, :), lr(:, :), old_lr(:, :), old_t(:, :)
      real(wp), allocatable :: wr(:, :), wl(:, :), wk(:, :), lambda(:, :)
      real(wp) :: rmse, ew, ec
      logical :: converged
      integer :: niter, step, a
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=2, verbosity=0)
      call new_uv_case(case, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      params%tolerance = 1.0e-12_wp
      params%damp_small_k = .false.
      params%precond = "none"
      params%solver%jfnk%restart = 6
      params%solver%jfnk%adjoint%restart = 5
      sr = -0.02_wp*case%cr
      lr = 0.0_wp*case%cr
      call build(warm)
      if (.not. allocated(err)) call warm%solve(params%tolerance, params%maxiter, converged, rmse, niter, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged, more="1D reuse reference did not converge")
      if (allocated(error)) return
      basis = c_loc(warm%solver%workspace%ws_v)
      do step = 1, 2
         old_t = warm%tr
         old_lr = lr
         sr = 1.01_wp*sr
         lr = lr + 0.001_wp*case%cr
         case%uk_lr = 1.01_wp*case%uk_lr
         call warm%update_potential(sr, lr, case%uk_lr, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, maxval(abs(warm%tr - old_t - beta*(old_lr - lr))), 0.0_wp, thr=1.0e-14_wp, &
            & more="1D potential update did not preserve physical t")
         if (allocated(error)) return
         call warm%solve(params%tolerance, params%maxiter, converged, rmse, niter, err, reuse_initial_guess=.true.)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged .and. c_associated(basis, c_loc(warm%solver%workspace%ws_v)), &
            & more="1D host update lost convergence or replaced the Krylov allocation")
         if (allocated(error)) return
         call build(cold)
         if (.not. allocated(err)) call cold%solve(params%tolerance, params%maxiter, converged, rmse, niter, err)
         if (.not. allocated(err)) call warm%free_energy(ew, err)
         if (.not. allocated(err)) call cold%free_energy(ec, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, converged .and. maxval(abs(warm%tr - cold%tr)) < 1.0e-10_wp, &
            & more="1D reused and cold fields differ")
         if (allocated(error)) return
         call check(error, ew, ec, thr=1.0e-10_wp, more="1D reused and cold energies differ")
         if (allocated(error)) return
      end do
      call warm%response(1.0e-10_wp, wr, wl, wk, lambda, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(warm%solver%workspace%ws_v, 2) == 7 .and. &
         & c_associated(basis, c_loc(warm%solver%workspace%ws_v)), more="1D response ignored its restart setting")
      if (allocated(error)) return
      old_t = warm%tr
      call warm%update_potential(sr, lr, case%uk_lr(:, :1), err)
      call check(error, allocated(err) .and. all(warm%tr == old_t), more="rejected 1D update changed the state")
      if (allocated(error)) return
      deallocate (err)

      ! With W=P=I and zero tails, T is zero. A cold start converges in one step,
      ! while two damped Picard updates cannot repair a large feasible t_s guess.
      case%wu = 0.0_wp
      case%response = 0.0_wp
      do a = 1, nu
         case%wu(:, case%lut%idx(a, a)) = 1.0_wp
      end do
      do a = 1, nv
         case%response(:, a, a) = 1.0_wp
      end do
      sr = 0.0_wp
      lr = 0.0_wp
      case%uk_lr = 0.0_wp
      params%solver%name = "PICARD"
      call build(warm)
      if (.not. allocated(err)) call warm%solve(params%tolerance, 2, converged, rmse, niter, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == 1)
      if (allocated(error)) return
      warm%tr = 1000.0_wp
      call warm%solve(params%tolerance, 2, converged, rmse, niter, err, reuse_initial_guess=.true.)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == 3 .and. all(warm%tr == 0.0_wp), &
         & more="failed feasible 1D warm attempt did not retry cold")
      if (allocated(error)) return
      warm%tr(1, 1) = ieee_value(0.0_wp, ieee_quiet_nan)
      call warm%solve(params%tolerance, 2, converged, rmse, niter, err, reuse_initial_guess=.true.)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, converged .and. niter == 1 .and. all(warm%tr == 0.0_wp), &
         & more="infeasible retained 1D guess did not recover")
      call ctx%delete()

   contains

      subroutine build(problem)
         type(uv_problem_type), intent(out) :: problem
         call new_uv_problem(problem, ctx, case%rgrid, case%kgrid, case%wu, case%response, sr, lr, case%uk_lr, &
            & spread(0.005_wp, 1, nv), beta, params, err)
      end subroutine build

   end subroutine check_scf_reuse

   !> Fixed operators match fresh sampling across geometry, grid and bulk changes
   subroutine check_operator_cache(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(solvent_vv_type) :: solvent, bulk_copy
      type(model_moz_1d_parameters_type) :: params
      type(custom_lj_type) :: lj
      type(fixed_charges_type) :: charges
      type(uv_operators_type) :: operators, copied
      type(moist_math_grid_radial_type) :: rgrid, kgrid, original_grid, bad_grid
      real(wp), parameter :: bulk_xyz(3, 2) = reshape([0.0_wp, 0.0_wp, 0.0_wp, &
         & 1.7_wp, 0.0_wp, 0.0_wp], [3, 2])
      real(wp), parameter :: bulk_q(2) = [-0.02_wp, 0.02_wp]
      real(wp), parameter :: original_xyz(3, 3) = reshape([0.0_wp, 0.0_wp, 0.0_wp, &
         & 1.7_wp, 0.2_wp, 0.0_wp, -0.8_wp, 1.2_wp, 0.5_wp], [3, 3])
      real(wp), allocatable :: held_response(:, :, :), held_w(:, :)
      real(wp) :: xyz(3, 3)

      ! Weak polar bulk keeps this cache test inexpensive; it needs no UV solves.
      call new_context(ctx, nthreads=1, verbosity=0)
      system%solvent%id = 0
      system%temperature = 1.0e5_wp
      system%solvent%number_density_au = 0.005_wp
      ! Ideal-dipole dielectric of the weak polar bulk for DRISM
      system%solvent%epsilon = 1.0_wp + 4.0_wp*acos(-1.0_wp)*atomic_unit_of_energy &
         & /(Boltzmann_constant*system%temperature)*system%solvent%number_density_au &
         & *sum(matmul(bulk_xyz, bulk_q)**2)/3.0_wp
      allocate (system%solv_mol)
      call new(system%solv_mol, [1, 1], bulk_xyz)
      call new_vv_solvent(solvent, ctx, system, err)
      if (.not. allocated(err)) call new_custom_lj(lj, [1.0_wp, 1.0_wp], [0.0_wp, 0.0_wp], err)
      if (.not. allocated(err)) call new_fixed_charges(charges, bulk_q, err)
      if (.not. allocated(err)) call solvent%potential%add(lj, err)
      if (.not. allocated(err)) call solvent%potential%add(charges, err)
      if (.not. allocated(err)) call solvent%update(err)
      params%solvent%npts = 512
      params%solvent%dr = 0.25_wp
      params%solvent%nlambda = 1
      if (.not. allocated(err)) call solvent%solve(err, params%solvent)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 128, 1.0_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      xyz = original_xyz
      original_grid = kgrid
      call check(error, .not. operators%response_matches(kgrid, solvent) .and. &
         & .not. operators%intramolecular_matches(kgrid, xyz), more="empty operators matched their inputs")
      if (allocated(error)) return
      call operators%prepare(kgrid, solvent, xyz, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call verify_operators(xyz)
      if (allocated(error)) return
      bulk_copy = solvent
      copied = operators
      held_response = operators%response
      held_w = operators%wu
      call operators%prepare(kgrid, bulk_copy, xyz, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, operators%response_matches(kgrid, bulk_copy) .and. &
         & operators%intramolecular_matches(kgrid, xyz) .and. &
         & all(operators%response == held_response) .and. all(operators%wu == held_w), &
         & more="unchanged inputs or copying the bulk changed cached operators")
      if (allocated(error)) return

      xyz(1, 2) = xyz(1, 2) + 0.1_wp
      call check(error, operators%response_matches(kgrid, solvent) .and. &
         & .not. operators%intramolecular_matches(kgrid, xyz), more="solute motion invalidated the wrong operator")
      if (allocated(error)) return
      call operators%prepare(kgrid, solvent, xyz, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(operators%response == held_response) .and. any(operators%wu /= held_w), &
         & more="solute motion did not preserve the solvent response and refresh W")
      if (allocated(error)) return
      call verify_operators(xyz)
      if (allocated(error)) return
      call check(error, copied%intramolecular_matches(original_grid, original_xyz), &
         & more="updating one cache changed an independent copy")
      if (allocated(error)) return
      call check(error, .not. operators%intramolecular_matches(kgrid, xyz(:, :2)), &
         & more="changing the number of solute sites kept W current")
      if (allocated(error)) return
      call operators%prepare(kgrid, solvent, xyz(:, :2), err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call verify_operators(xyz(:, :2))
      if (allocated(error)) return

      ! Same point count with different nodes must invalidate both tables.
      call new_uniform_radial_pair(rgrid, kgrid, 128, 0.9_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. operators%response_matches(kgrid, solvent) .and. &
         & .not. operators%intramolecular_matches(kgrid, xyz), more="changed reciprocal grid kept operators current")
      if (allocated(error)) return
      call operators%prepare(kgrid, solvent, xyz, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call verify_operators(xyz)
      if (allocated(error)) return
      call check(error, copied%response_matches(original_grid, bulk_copy), &
         & more="changing a grid changed the key of an independent cache copy")
      if (allocated(error)) return

      call new_uniform_radial_pair(rgrid, kgrid, 64, 0.9_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. operators%response_matches(kgrid, solvent) .and. &
         & .not. operators%intramolecular_matches(kgrid, xyz), more="new grid size kept operators current")
      if (allocated(error)) return
      call operators%prepare(kgrid, solvent, xyz, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call verify_operators(xyz)
      if (allocated(error)) return
      held_response = operators%response
      held_w = operators%wu

      ! An out-of-band refresh fails without publishing new keys or data.
      call new_uniform_radial_pair(rgrid, bad_grid, 128, 0.1_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call operators%prepare(bad_grid, solvent, xyz, err)
      call check(error, allocated(err), more="out-of-band operator sampling was accepted")
      if (allocated(error)) return
      call check(error, operators%response_matches(kgrid, solvent) .and. &
         & operators%intramolecular_matches(kgrid, xyz) .and. &
         & all(operators%response == held_response) .and. all(operators%wu == held_w), &
         & more="failed preparation damaged the previously accepted cache")
      if (allocated(error)) return

      call solvent%potential%set_ng_split(.true., 0.8_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. operators%response_matches(kgrid, solvent) .and. &
         & operators%response_matches(kgrid, bulk_copy), more="stale bulk matched a cached response")
      if (allocated(error)) return
      call operators%prepare(kgrid, solvent, xyz, err)
      call check(error, allocated(err), more="operators accepted an unsolved bulk solvent")
      if (allocated(error)) return

      ! A fresh accepted bulk identity also changes RISM to the analytical DRISM path.
      params%solvent%theory = "DRISM"
      call solvent%update(err)
      if (.not. allocated(err)) call solvent%solve(err, params%solvent)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. operators%response_matches(kgrid, solvent) .and. &
         & operators%intramolecular_matches(kgrid, xyz), more="new bulk identity invalidated the wrong operator")
      if (allocated(error)) return
      call operators%prepare(kgrid, solvent, xyz, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call verify_operators(xyz)
      if (allocated(error)) return
      call check(error, all(operators%wu == held_w) .and. copied%response_matches(original_grid, bulk_copy), &
         & more="bulk refresh changed W or invalidated an independent cache")

   contains

      !> Compare cached values with the sampling and correction used before caching.
      subroutine verify_operators(coords)
         real(wp), intent(in) :: coords(:, :)
         type(packed_index_lut) :: lut
         type(potential_sites_type) :: sites
         real(wp), allocatable :: chi(:, :), z(:, :), expected(:, :, :), expected_w(:, :)
         integer :: a, b, ns
         logical :: applied

         if (solvent%theory_name() == "DRISM") then
            call solvent%susceptibility_oz_at(kgrid%r, chi, err, z)
         else
            call solvent%susceptibility_at(kgrid%r, chi, err)
         end if
         if (.not. allocated(err)) call solvent%potential%sites(sites, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         ns = solvent%nsite()
         call lut%init(ns)
         allocate (expected(kgrid%npts, ns, ns))
         do a = 1, ns
            do b = 1, ns
               expected(:, b, a) = chi(:, lut%idx(a, b))
            end do
         end do
         call uv_charge_moment_correction(kgrid%r, sites%q, expected, applied, z)
         call check(error, applied, more="operator cache fixture did not exercise charge-moment correction")
         if (allocated(error)) return
         call lut%init(size(coords, 2))
         allocate (expected_w(kgrid%npts, npair_from_ns(size(coords, 2))))
         call vv_intramolecular(kgrid%r, coords, lut, expected_w)
         call check(error, all(shape(operators%response) == shape(expected)) .and. &
            & all(shape(operators%wu) == shape(expected_w)), more="cached operators have incorrect shapes")
         if (allocated(error)) return
         call check(error, operators%response_matches(kgrid, solvent) .and. &
            & operators%intramolecular_matches(kgrid, coords), more="prepared operators do not match their inputs")
         if (allocated(error)) return
         call check(error, all(operators%response == expected) .and. all(operators%wu == expected_w), &
            & more="caching changed the existing sampling or intramolecular mathematics")
      end subroutine verify_operators

   end subroutine check_operator_cache

   !> Check the charge-moment correction of the solvent response
   !>
   !> - A symmetric positive-semidefinite P has a known quadratic S_a/k^2.
   !> - A symmetric low-k error is repaired while preserving symmetry and the
   !>   charge-neutral block; modes at and above 0.05/bohr are untouched.
   !> - Repeated correction and an uncharged solvent leave P unchanged.
   !> - Reversing all solvent charges gives the same corrected response.
   subroutine check_charge_moment(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: q(nv) = [-0.8476_wp, 0.4238_wp, 0.4238_wp]
      real(wp), parameter :: k(10) = [0.0_wp, 0.005_wp, 0.02_wp, 0.049_wp, 0.05_wp, &
         & 0.1_wp, 0.15_wp, 0.2_wp, 0.25_wp, 0.3_wp]
      real(wp), parameter :: aa(nv, nv) = reshape([1.4_wp, 0.2_wp, -0.1_wp, &
         & 0.2_wp, 1.1_wp, 0.15_wp, -0.1_wp, 0.15_wp, 0.9_wp], [nv, nv])
      real(wp), parameter :: neutral(nv, 2) = reshape([1.0_wp, 1.0_wp, 1.0_wp, &
         & 0.0_wp, 1.0_wp, -1.0_wp], [nv, 2])
      real(wp) :: p(size(k), nv, nv), corrupt(size(k), nv, nv), reversed(size(k), nv, nv)
      real(wp) :: change(nv, nv), s(nv), moment
      logical :: applied
      integer :: j, a, nlow

      ! A is symmetric positive definite; P = 11^T + k^2 f(k^2) A has P q = k^2 f A q.
      s = matmul(aa, q)
      do j = 1, size(k)
         p(j, :, :) = 1.0_wp + k(j)**2*(1.0_wp + 0.3_wp*k(j)**2 + 0.2_wp*k(j)**4)*aa
      end do
      nlow = count(k < 0.05_wp)
      corrupt = p
      do j = 1, nlow
         corrupt(j, 1, 1) = corrupt(j, 1, 1) + 0.02_wp
      end do
      p = corrupt
      call uv_charge_moment_correction(k, q, p, applied)
      call check(error, applied, more="the correction was not applied")
      if (allocated(error)) return
      do j = 1, nlow
         do a = 1, nv
            moment = s(a)*k(j)**2*(1.0_wp + 0.3_wp*k(j)**2 + 0.2_wp*k(j)**4)
            call check(error, dot_product(q, p(j, :, a)), moment, thr=1.0e-12_wp, more="charge moment not restored")
            if (allocated(error)) return
         end do
         call check(error, maxval(abs(p(j, :, :) - transpose(p(j, :, :)))) < 1.0e-14_wp, &
            & more="fitted charge-moment correction broke susceptibility symmetry")
         if (allocated(error)) return
         change = p(j, :, :) - corrupt(j, :, :)
         call check(error, maxval(abs(matmul(transpose(neutral), matmul(change, neutral)))) < 1.0e-14_wp, &
            & more="fitted charge-moment correction changed the charge-neutral block")
         if (allocated(error)) return
      end do
      call check(error, all(p(nlow + 1:, :, :) == corrupt(nlow + 1:, :, :)), &
         & more="the correction changed modes at or above its cutoff")
      if (allocated(error)) return
      reversed = corrupt
      call uv_charge_moment_correction(k, -q, reversed, applied)
      call check(error, applied .and. maxval(abs(reversed - p)) < 1.0e-14_wp, &
         & more="reversing charges changed the fitted response correction")
      if (allocated(error)) return
      corrupt = p
      call uv_charge_moment_correction(k, q, p, applied)
      call check(error, applied .and. maxval(abs(p - corrupt)) < 1.0e-14_wp, &
         & more="repeated fitted correction changed an already restored response")
      if (allocated(error)) return
      corrupt = p
      call uv_charge_moment_correction(k, [0.0_wp, 0.0_wp, 0.0_wp], p, applied)
      call check(error, .not. applied .and. all(p == corrupt), more="an uncharged solvent was corrected")
   end subroutine check_charge_moment

   !> Settings of the energy tests: the solvent and UV solves on the same coarse pair
   !>
   !> @param[out] params  settings
   subroutine energy_settings(params)
      !> Settings
      type(model_moz_1d_parameters_type), intent(out) :: params
      params%tolerance = 1.0e-10_wp
      params%solvent%npts = 1024
      params%solvent%dr = 0.06_wp
      params%solvent%tolerance = 1.0e-10_wp
   end subroutine energy_settings

   !> SPC/E water solvent and a 1D model on the solvent's pair with the given settings
   !>
   !> @param[in]  ctx      run context
   !> @param[in]  params   settings
   !> @param[out] solvent  updated, unsolved solvent
   !> @param[out] model    constructed model without solute terms
   !> @param[out] err      construction failure
   subroutine new_water_model(ctx, params, solvent, model, err)
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Settings
      type(model_moz_1d_parameters_type), intent(in) :: params
      !> Solvent
      type(solvent_vv_type), intent(out) :: solvent
      !> Model
      type(model_moz_1d_type), intent(out) :: model
      !> Construction failure
      type(moist_error), allocatable, intent(out) :: err
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, params%solvent%npts, params%solvent%dr, err)
      if (.not. allocated(err)) call new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, err, params)
   end subroutine new_water_model

   !> Energy through a coupling of the model, answering host charges when requested
   !>
   !> @param[inout] model   updated model
   !> @param[in]    charges host charges for an atomic_charges request; none when absent
   !> @param[inout] energy  energy, Hartree
   !> @param[out]   err     coupling or energy failure
   subroutine model_energy(model, energy, err, charges)
      !> Model
      type(model_moz_1d_type), intent(inout), target :: model
      !> Energy
      real(wp), intent(inout) :: energy
      !> Coupling or energy failure
      type(moist_error), allocatable, intent(out) :: err
      !> Host charges
      real(wp), intent(in), optional :: charges(:)
      type(coupling_type), pointer :: coupling

      call model%new_coupling(coupling, err)
      if (allocated(err)) return
      call model%prepare_energy(coupling, err)
      if (.not. allocated(err) .and. present(charges)) then
         if (coupling%next()) call coupling%answer("q", charges, err)
      end if
      if (.not. allocated(err)) call model%get_energy(coupling, energy, err)
      call model%release_coupling(coupling)
   end subroutine model_energy

   !> Check a water solute in water against the VV solution
   !>
   !> In ordinary RISM, a solute identical to the solvent obeys the VV equation, so at
   !> infinite dilution h_uv = h_vv and the solvation free energy equals the
   !> excess chemical potential of the solvent and the partial molar volume the
   !> molecular volume 1/rho
   subroutine check_water_in_water(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_1d_type), target :: model
      type(model_moz_1d_parameters_type) :: params
      type(solvent_vv_thermodynamics_type) :: thermo
      type(model_moz_1d_thermodynamics_type) :: uv
      type(packed_index_lut) :: lut
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :), cvv(:, :), tvv(:, :), hvv(:, :)
      real(wp) :: energy
      integer :: i, a

      call new_context(ctx, nthreads=0, verbosity=0)
      call energy_settings(params)
      call new_water_model(ctx, params, solvent, model, err)
      if (.not. allocated(err)) call model%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call model%potential%add(coulomb_spce, err)
      if (.not. allocated(err)) call model%update(solvent%structure(), err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call model%thermodynamics(uv, err)
      call check(error, allocated(err), more="thermodynamics before an energy were returned")
      if (allocated(error)) return
      deallocate (err)
      energy = 0.0_wp
      if (.not. allocated(err)) call model_energy(model, energy, err)
      if (.not. allocated(err)) call model%correlations(cr, tr, hr, err)
      if (.not. allocated(err)) call model%thermodynamics(uv, err)
      if (.not. allocated(err)) call solvent%solve(err, params%solvent)
      if (.not. allocated(err)) call solvent%thermodynamics(thermo, err)
      if (.not. allocated(err)) call solvent%correlations(cvv, tvv, hvv, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, energy, thermo%chemical_potential, rel=.true., thr=1.0e-7_wp, &
         & more="water in water differs from the solvent excess chemical potential")
      if (allocated(error)) return
      call check(error, energy, uv%free_energy, thr=1.0e-14_wp, more="thermodynamic free energy differs from the reported energy")
      if (allocated(error)) return
      call check(error, uv%partial_molar_volume*solvent%density(), 1.0_wp, rel=.true., thr=1.0e-6_wp, &
         & more="partial molar volume of water in water differs from 1/rho")
      if (allocated(error)) return
      call lut%init(3)
      do i = 1, 3
         do a = 1, 3
            call check(error, maxval(abs(hr(:, (i - 1)*3 + a) - hvv(:, lut%idx(i, a)))), 0.0_wp, thr=1.0e-8_wp, &
               & more="h_uv of water in water differs from h_vv")
            if (allocated(error)) return
         end do
      end do
      call model%invalidate()
      call model%thermodynamics(uv, err)
      call check(error, allocated(err), more="invalidation kept thermodynamics available")
   end subroutine check_water_in_water

   !> Check the table adjoints of the UV problem against central differences
   !>
   !> - Water solute in SPC/E water on the coarse pair of the energy tests
   !> - E = a dG for unit, nonunit and zero free-energy sensitivities; for a smooth perturbation
   !>   of each table, sum(w du) matches the central difference of E; the k-space
   !>   table moves E about 1e4 times less and takes a larger step
   !> - A warm start from the previous lambda gives the same adjoints
   subroutine check_response_tables(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_1d_parameters_type) :: params
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moz_potential_type) :: solute
      type(potential_sites_type) :: sites
      type(packed_index_lut) :: lut_u, lut_v
      type(uv_problem_type) :: problem
      type(structure_type) :: mol
      real(wp), parameter :: weights(3) = [1.0_wp, 0.0_wp, 2.5_wp]
      real(wp), parameter :: steps(3) = [1.0e-6_wp, 1.0e-6_wp, 1.0e-4_wp]
      real(wp) :: eps
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :), uk_lr(:, :), chi(:, :), response(:, :, :), wu(:, :)
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), w_k(:, :), lambda(:, :), du(:, :)
      real(wp), allocatable :: v_sr(:, :), v_lr(:, :), v_k(:, :)
      real(wp) :: g_plus, g_minus, fd, dot
      integer :: nu, nv, a, b, i, col, itab, is
      logical :: applied

      call new_context(ctx, nthreads=0, verbosity=0)
      call energy_settings(params)
      params%tolerance = 1.0e-12_wp
      params%solvent%tolerance = 1.0e-12_wp
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call solvent%solve(err, params%solvent)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, params%solvent%npts, params%solvent%dr, err)
      if (.not. allocated(err)) call solute%add(lj_cspce, err)
      if (.not. allocated(err)) call solute%add(coulomb_spce, err)
      if (.not. allocated(err)) call solute%update(solvent%structure(), rgrid, err)
      if (.not. allocated(err)) call solute%compute(rgrid, kgrid, u_sr, ur_lr, uk_lr, err, solvent=solvent%potential)
      if (.not. allocated(err)) call solvent%susceptibility_at(kgrid%r, chi, err)
      if (.not. allocated(err)) call solvent%potential%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      nv = solvent%nsite()
      nu = 3
      call lut_v%init(nv)
      allocate (response(size(chi, 1), nv, nv))
      do a = 1, nv
         do b = 1, nv
            response(:, b, a) = chi(:, lut_v%idx(a, b))
         end do
      end do
      call uv_charge_moment_correction(kgrid%r, sites%q, response, applied)
      call lut_u%init(nu)
      allocate (wu(size(kgrid%r), npair_from_ns(nu)))
      mol = solvent%structure()
      call vv_intramolecular(kgrid%r, mol%xyz, lut_u, wu)

      do is = 1, size(weights)
         call solve_terms(u_sr, ur_lr, uk_lr, g_plus)
         if (allocated(error)) return
         call problem%response(1.0e-10_wp, w_sr, w_lr, w_k, lambda, err, &
            & free_energy_weight=weights(is))
         call check_moist_error(error, err)
         if (allocated(error)) return
         do itab = 1, 3
            eps = steps(itab)
            allocate (du, mold=u_sr)
            do col = 1, size(du, 2)
               do i = 1, size(du, 1)
                  if (itab == 3) then
                     du(i, col) = (1.0_wp + 0.2_wp*col)*exp(-((kgrid%r(i) - 1.5_wp)/0.8_wp)**2)
                  else
                     du(i, col) = (1.0_wp + 0.2_wp*col)*exp(-((rgrid%r(i) - 6.0_wp)/1.5_wp)**2)
                  end if
               end do
            end do
            v_sr = u_sr
            v_lr = ur_lr
            v_k = uk_lr
            select case (itab)
            case (1)
               v_sr = u_sr + eps*du
               dot = sum(w_sr*du)
            case (2)
               v_lr = ur_lr + eps*du
               dot = sum(w_lr*du)
            case default
               v_k = uk_lr + eps*du
               dot = sum(w_k*du)
            end select
            call solve_terms(v_sr, v_lr, v_k, g_plus)
            if (allocated(error)) return
            select case (itab)
            case (1)
               v_sr = u_sr - eps*du
            case (2)
               v_lr = ur_lr - eps*du
            case default
               v_k = uk_lr - eps*du
            end select
            call solve_terms(v_sr, v_lr, v_k, g_minus)
            if (allocated(error)) return
            fd = weights(is)*(g_plus - g_minus)/(2.0_wp*eps)
            call check(error, dot, fd, rel=.true., thr=1.0e-5_wp, more="table adjoint differs from central differences")
            if (allocated(error)) return
            deallocate (du)
         end do
      end do

      ! Warm start from the last lambda at the last free-energy weight
      call solve_terms(u_sr, ur_lr, uk_lr, g_plus)
      if (allocated(error)) return
      v_sr = w_sr
      call problem%response(1.0e-10_wp, w_sr, w_lr, w_k, lambda, err, &
         & free_energy_weight=weights(size(weights)))
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(w_sr - v_sr)) <= 1.0e-8_wp*maxval(abs(v_sr)), &
         & more="warm-started response differs")
      if (allocated(error)) return

      ! A zero objective must also have a usable tolerance and zero table adjoints.
      call problem%response(1.0e-10_wp, w_sr, w_lr, w_k, lambda, err, free_energy_weight=0.0_wp)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(w_sr)) + maxval(abs(w_lr)) + maxval(abs(w_k)) + maxval(abs(lambda)), 0.0_wp, &
         & more="zero objective has nonzero table adjoints")

   contains

      !> Solve with the given tables and return the free energy
      subroutine solve_terms(t_sr, t_lr, t_k, g)
         real(wp), intent(in) :: t_sr(:, :), t_lr(:, :), t_k(:, :)
         real(wp), intent(out) :: g
         logical :: converged
         real(wp) :: rmse
         integer :: niter
         g = 0.0_wp
         call new_uv_problem(problem, ctx, rgrid, kgrid, wu, response, t_sr, t_lr, t_k, &
            & spread(solvent%density(), 1, nv), solvent%beta(), params, err)
         if (.not. allocated(err)) call problem%solve(params%tolerance, params%maxiter, converged, rmse, niter, err)
         if (.not. allocated(err) .and. .not. converged) then
            call test_failed(error, "UV solve of the response test did not converge")
            return
         end if
         if (.not. allocated(err)) call problem%free_energy(g, err)
         call check_moist_error(error, err)
      end subroutine solve_terms

   end subroutine check_response_tables

   !> Response of the model through a coupling, answering host charges in both phases
   !>
   !> @param[inout] model    updated model
   !> @param[in]    charges  host charges
   !> @param[out]   dg_dq    atomic charge adjoint; unallocated when the response has none
   !> @param[out]   nitem    number of response items
   !> @param[out]   err      coupling or response failure
   subroutine model_response(model, charges, dg_dq, nitem, err)
      !> Model
      type(model_moz_1d_type), intent(inout), target :: model
      !> Host charges
      real(wp), intent(in) :: charges(:)
      !> Atomic charge adjoint
      real(wp), allocatable, intent(out) :: dg_dq(:)
      !> Number of response items
      integer, intent(out) :: nitem
      !> Coupling or response failure
      type(moist_error), allocatable, intent(out) :: err
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      class(response_item_type), allocatable :: item

      nitem = 0
      call model%new_coupling(coupling, err)
      if (allocated(err)) return
      call model%prepare_response(coupling, err)
      if (.not. allocated(err)) then
         if (coupling%next()) call coupling%answer("q", charges, err)
      end if
      if (.not. allocated(err)) call model%get_response(coupling, response, err)
      call model%release_coupling(coupling)
      if (allocated(err)) return
      do while (response%next())
         nitem = nitem + 1
         item = response%item()
         select type (item)
         type is (atomic_charge_adjoint_response_type)
            dg_dq = item%dg_dq
         end select
      end do
   end subroutine model_response

   !> Check the charge response of a water solute with host charges
   !>
   !> - dE/dq_i of the reported energy against central differences; tight solves
   !>   keep the differences above their noise
   !> - A second energy with the same charges reuses the solution
   !> - A response without an energy before it equals the one after it
   subroutine check_response_charges(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_1d_type), target :: model, cold
      type(model_moz_1d_parameters_type) :: params
      type(monopole_type) :: charges
      type(potential_sites_type) :: sites
      real(wp), parameter :: eps = 1.0e-4_wp
      real(wp), allocatable :: q(:), dg_dq(:), dg_cold(:)
      real(wp) :: e0, e1, e_plus, e_minus, e_cold
      integer :: i, nitem

      call new_context(ctx, nthreads=0, verbosity=0)
      call energy_settings(params)
      params%tolerance = 1.0e-12_wp
      params%solvent%tolerance = 1.0e-12_wp
      call new_water_model(ctx, params, solvent, model, err)
      if (.not. allocated(err)) call new_water_model(ctx, params, solvent, cold, err)
      if (.not. allocated(err)) call model%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call model%potential%add(charges, err)
      if (.not. allocated(err)) call cold%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call cold%potential%add(charges, err)
      if (.not. allocated(err)) call model%update(solvent%structure(), err)
      if (.not. allocated(err)) call cold%update(solvent%structure(), err)
      if (.not. allocated(err)) call solvent%potential%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      q = sites%q

      e0 = 0.0_wp
      e1 = 0.0_wp
      call model_energy(model, e0, err, q)
      if (.not. allocated(err)) call model_energy(model, e1, err, q)
      if (.not. allocated(err)) call model_response(model, q, dg_dq, nitem, err)
      if (.not. allocated(err)) call model_response(cold, q, dg_cold, nitem, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, e1 == e0, more="a repeated energy with the same charges changed")
      if (allocated(error)) return
      call check(error, allocated(dg_dq) .and. allocated(dg_cold), more="no atomic charge adjoint")
      if (allocated(error)) return
      call check(error, size(dg_dq) == 3 .and. nitem == 1, more="unexpected response items")
      if (allocated(error)) return
      call check(error, maxval(abs(dg_cold - dg_dq)) <= 1.0e-8_wp*maxval(abs(dg_dq)), &
         & more="response without an energy before it differs")
      if (allocated(error)) return

      do i = 1, size(q)
         e_plus = 0.0_wp
         e_minus = 0.0_wp
         q(i) = sites%q(i) + eps
         call model_energy(model, e_plus, err, q)
         if (i == 1 .and. .not. allocated(err)) then
            call cold%update(solvent%structure(), err)
            e_cold = 0.0_wp
            if (.not. allocated(err)) call model_energy(cold, e_cold, err, q)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, e_plus, e_cold, thr=1.0e-10_wp, &
               & more="a host-charge change gives different warm and fresh energies")
            if (allocated(error)) return
         end if
         q(i) = sites%q(i) - eps
         if (.not. allocated(err)) call model_energy(model, e_minus, err, q)
         q(i) = sites%q(i)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, dg_dq(i), (e_plus - e_minus)/(2.0_wp*eps), rel=.true., thr=2.0e-6_wp, &
            & more="charge response differs from central differences")
         if (allocated(error)) return
      end do
   end subroutine check_response_charges

   !> Check that fixed charges give an empty response
   subroutine check_response_fixed(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_1d_type), target :: model
      type(model_moz_1d_parameters_type) :: params
      type(coupling_type), pointer :: coupling
      type(response_type) :: response

      call new_context(ctx, nthreads=0, verbosity=0)
      call energy_settings(params)
      call new_water_model(ctx, params, solvent, model, err)
      if (.not. allocated(err)) call model%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call model%potential%add(coulomb_spce, err)
      if (.not. allocated(err)) call model%update(solvent%structure(), err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_response(coupling, err)
      if (.not. allocated(err)) then
         call check(error, .not. coupling%next(), more="fixed charges requested host data")
         if (allocated(error)) return
         call model%get_response(coupling, response, err)
      end if
      if (associated(coupling)) call model%release_coupling(coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. response%next(), more="fixed charges gave a response item")
   end subroutine check_response_fixed

   !> Check that host charges through the coupling equal fixed table charges
   subroutine check_host_charges(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_1d_type), target :: fixed, host
      type(model_moz_1d_parameters_type) :: params
      type(monopole_type) :: charges
      type(potential_sites_type) :: sites
      real(wp) :: e_fixed, e_host

      call new_context(ctx, nthreads=0, verbosity=0)
      call energy_settings(params)
      call new_water_model(ctx, params, solvent, fixed, err)
      if (.not. allocated(err)) call new_water_model(ctx, params, solvent, host, err)
      if (.not. allocated(err)) call fixed%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call fixed%potential%add(coulomb_spce, err)
      if (.not. allocated(err)) call host%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call host%potential%add(charges, err)
      if (.not. allocated(err)) call fixed%update(solvent%structure(), err)
      if (.not. allocated(err)) call host%update(solvent%structure(), err)
      if (.not. allocated(err)) call solvent%potential%sites(sites, err)
      e_fixed = 0.0_wp
      e_host = 0.0_wp
      if (.not. allocated(err)) call model_energy(fixed, e_fixed, err)
      if (.not. allocated(err)) call model_energy(host, e_host, err, sites%q)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, e_host, e_fixed, rel=.true., thr=1.0e-12_wp, more="host charges changed the energy")
   end subroutine check_host_charges

   !> Check refused energies and stale correlations
   !>
   !> - No correlations before an energy
   !> - A UV solve that does not converge is an error and leaves the energy unchanged
   !> - An update drops the correlations
   subroutine check_energy_refusals(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_1d_type), target :: model
      type(model_moz_1d_parameters_type) :: params
      type(structure_type) :: mol
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :)
      real(wp) :: energy

      call new_context(ctx, nthreads=0, verbosity=0)
      call energy_settings(params)
      params%maxiter = 1
      call new_water_model(ctx, params, solvent, model, err)
      if (.not. allocated(err)) call model%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call model%potential%add(coulomb_spce, err)
      if (.not. allocated(err)) call model%update(solvent%structure(), err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call model%correlations(cr, tr, hr, err)
      call check(error, allocated(err), more="correlations before an energy were returned")
      if (allocated(error)) return
      deallocate (err)
      energy = 7.0_wp
      call model_energy(model, energy, err)
      call check(error, allocated(err) .and. energy == 7.0_wp, more="a one-iteration UV solve was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "did not converge") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      call energy_settings(params)
      call new_water_model(ctx, params, solvent, model, err)
      if (.not. allocated(err)) call model%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call model%potential%add(coulomb_spce, err)
      if (.not. allocated(err)) call model%update(solvent%structure(), err)
      if (.not. allocated(err)) call model_energy(model, energy, err)
      if (.not. allocated(err)) call model%correlations(cr, tr, hr, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = solvent%structure()
      mol%xyz(1, 2) = mol%xyz(1, 2) + 0.1_wp
      call model%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call model%correlations(cr, tr, hr, err)
      call check(error, allocated(err), more="an update kept the correlations")
   end subroutine check_energy_refusals

   !> Check methane and charged solutes in TIP3P water against pysofft, XRISM/KH
   !>
   !> Reference: moist_dev origin/pysofft (d3c3d076), new_rism_vv_problem and
   !> new_rism_uv_problem with the explicit force fields of its rism_smoke test,
   !> TIP3P at 298 K and 0.0334 1/A^3, 8192 x 0.025 A, tolerances 1e-10:
   !> 12.298556534797479 and -44.558275999089851 kcal/mol. The charged cases exercise
   !> the charge-moment correction of the solvent response.
   !> Fresh pysofft c15bd85a runs also include an anion made by reversing a_q
   !> (a synthetic fixture, not a physical alanine anion). Its legacy energy is
   !> -254.339837058849724 kcal/mol. Replacing pysofft's rank-one correction by an
   !> independent construction of the neutral block R P R, R = I - q q^T/(q^T q),
   !> and both prescribed charge blocks gives -44.558273534779538 (cation) and
   !> -254.339839980442065 (anion). These symmetric references retain the 1e-6
   !> kcal/mol check; the legacy comparisons allow the measured 3e-6 shift.
   subroutine check_pysofft_reference(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(solvent_vv_type) :: solvent
      type(model_moz_1d_parameters_type) :: params
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(custom_lj_type) :: lj
      type(fixed_charges_type) :: charges
      real(wp), parameter :: w_sig(3) = [3.15_wp, 0.40_wp, 0.40_wp], w_eps(3) = [0.15201_wp, 0.04601_wp, 0.04601_wp]
      real(wp), parameter :: w_q(3) = [-0.834_wp, 0.417_wp, 0.417_wp]
      real(wp), parameter :: w_xyz(3, 3) = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.7566950_wp, 0.5858800_wp, &
         & 0.0_wp, -0.7566950_wp, 0.5858800_wp], [3, 3])
      real(wp), parameter :: m_sig(5) = [3.5_wp, 2.5_wp, 2.5_wp, 2.5_wp, 2.5_wp]
      real(wp), parameter :: m_eps(5) = [0.066001_wp, 0.030000_wp, 0.030000_wp, 0.030000_wp, 0.030000_wp]
      real(wp), parameter :: m_q(5) = [-0.351617_wp, 0.087864_wp, 0.087821_wp, 0.087835_wp, 0.088096_wp]
      real(wp), parameter :: m_xyz(3, 5) = reshape([0.000000_wp, 0.000000_wp, 0.000000_wp, &
         & -0.245808_wp, 0.265985_wp, -1.028259_wp, -0.639868_wp, -0.820203_wp, 0.326070_wp, &
         & 1.043914_wp, -0.309230_wp, 0.055734_wp, -0.158239_wp, 0.863449_wp, 0.646454_wp], [3, 5])
      real(wp), parameter :: a_xyz(3, 13) = reshape([ &
         & -1.84407_wp, 2.05773_wp, 0.16140_wp, -0.48394_wp, 1.43397_wp, -0.17933_wp, &
         & -1.98200_wp, 2.11380_wp, 1.26276_wp, -2.66381_wp, 1.44534_wp, -0.27055_wp, &
         & -1.91853_wp, 3.08218_wp, -0.26201_wp, 0.63630_wp, 2.28786_wp, 0.37479_wp, &
         & 0.75185_wp, 3.44838_wp, 0.01369_wp, 1.52226_wp, 1.78208_wp, 1.25607_wp, &
         & -0.41366_wp, 0.06093_wp, 0.33742_wp, -0.53186_wp, 0.04945_wp, 1.37726_wp, &
         & 0.50757_wp, -0.36787_wp, 0.08832_wp, -0.37764_wp, 1.40201_wp, -1.28698_wp, &
         & -1.16776_wp, -0.51952_wp, -0.09704_wp], [3, 13])
      real(wp), parameter :: a_sig(13) = [3.39770953124363_wp, 3.39770953124363_wp, 2.60017699876439_wp, &
         & 2.60017699876439_wp, 2.60017699876439_wp, 3.31521230994383_wp, 3.04812087424536_wp, &
         & 3.04812087424536_wp, 2.85586493087067_wp, 1.10649620793030_wp, 1.10649620793030_wp, &
         & 2.60017699876439_wp, 1.10649620793030_wp]
      real(wp), parameter :: a_eps(13) = [0.1078_wp, 0.1078_wp, 0.0208_wp, 0.0208_wp, 0.0208_wp, 0.0988_wp, &
         & 0.1463_wp, 0.1463_wp, 0.7828_wp, 0.0100_wp, 0.0100_wp, 0.0208_wp, 0.0100_wp]
      real(wp), parameter :: a_q(13) = [-0.10_wp, -0.04_wp, 0.05_wp, 0.05_wp, 0.05_wp, 0.50_wp, -0.25_wp, &
         & -0.25_wp, -0.29_wp, 0.43_wp, 0.43_wp, -0.01_wp, 0.43_wp]
      real(wp) :: energy

      call new_context(ctx, nthreads=0, verbosity=0)
      system%solvent%id = 0
      system%temperature = 298.0_wp
      system%solvent%number_density_au = 0.0334_wp/aatoau**3
      allocate (system%solv_mol)
      call new(system%solv_mol, [8, 1, 1], w_xyz*aatoau)
      call new_vv_solvent(solvent, ctx, system, err)
      if (.not. allocated(err)) call new_custom_lj(lj, w_sig*aatoau, w_eps/autokcal, err)
      if (.not. allocated(err)) call solvent%potential%add(lj, err)
      if (.not. allocated(err)) call new_fixed_charges(charges, w_q, err)
      if (.not. allocated(err)) call solvent%potential%add(charges, err)
      if (.not. allocated(err)) call solvent%update(err)
      params%tolerance = 1.0e-10_wp
      params%maxiter = 500
      params%solvent%npts = 8192
      params%solvent%dr = 0.025_wp*aatoau
      params%solvent%tolerance = 1.0e-10_wp
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 8192, 0.025_wp*aatoau, err)
      if (.not. allocated(err)) call solvent%solve(err, params%solvent)
      call check_moist_error(error, err)
      if (allocated(error)) return

      call solute_energy([6, 1, 1, 1, 1], m_sig, m_eps, m_q, m_xyz, energy)
      if (allocated(error)) return
      call check(error, energy*autokcal, 12.298556534797479_wp, thr=1.0e-6_wp, more="methane differs from pysofft")
      if (allocated(error)) return
      call solute_energy([6, 6, 1, 1, 1, 6, 8, 8, 7, 1, 1, 1, 1], a_sig, a_eps, a_q, a_xyz, energy)
      if (allocated(error)) return
      call check(error, energy*autokcal, -44.558273534779538_wp, thr=1.0e-6_wp, &
         & more="alanine cation differs from the symmetric pysofft reference")
      if (allocated(error)) return
      call check(error, energy*autokcal, -44.558275999089851_wp, thr=5.0e-6_wp, &
         & more="alanine cation differs from pysofft")
      if (allocated(error)) return
      call solute_energy([6, 6, 1, 1, 1, 6, 8, 8, 7, 1, 1, 1, 1], a_sig, a_eps, -a_q, a_xyz, energy)
      if (allocated(error)) return
      call check(error, energy*autokcal, -254.339839980442065_wp, thr=1.0e-6_wp, &
         & more="charge-reversed anion differs from the symmetric pysofft reference")
      if (allocated(error)) return
      call check(error, energy*autokcal, -254.339837058849724_wp, thr=5.0e-6_wp, &
         & more="charge-reversed anion differs from pysofft")

   contains

      !> Solvation free energy of one solute with explicit parameters, Hartree
      subroutine solute_energy(num, sig, eps, q, xyz, dg)
         integer, intent(in) :: num(:)
         real(wp), intent(in) :: sig(:), eps(:), q(:), xyz(:, :)
         real(wp), intent(out) :: dg
         type(model_moz_1d_type), target :: model
         type(custom_lj_type) :: ulj
         type(fixed_charges_type) :: uq
         type(structure_type) :: mol
         dg = 0.0_wp
         call new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, err, params)
         if (.not. allocated(err)) call new_custom_lj(ulj, sig*aatoau, eps/autokcal, err)
         if (.not. allocated(err)) call model%potential%add(ulj, err)
         if (.not. allocated(err)) call new_fixed_charges(uq, q, err)
         if (.not. allocated(err)) call model%potential%add(uq, err)
         call new(mol, num, xyz*aatoau)
         mol%charge = sum(q)
         if (.not. allocated(err)) call model%update(mol, err)
         if (.not. allocated(err)) call model_energy(model, dg, err)
         call check_moist_error(error, err)
      end subroutine solute_energy

   end subroutine check_pysofft_reference

end module test_moz_1d_uv
