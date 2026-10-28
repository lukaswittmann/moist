!> Tests of the 1D VV equation, evaluator and solve: moz_solvent_vv
!>
!> SPC/E water on a uniform radial pair with Ng-split tables from the solvent
!> potential and a Gaussian model c_s(r), and full water solves on a coarse
!> pair; each test runs as a selected test
module test_moz_solvent_vv
   !$ use omp_lib, only: omp_in_parallel
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, &
      & ieee_positive_inf, ieee_negative_inf
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type
   use mctc_io_convert, only: aatoau
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed, skip_test
   use moist_context, only: moist_context_type, new_context
   use moist_math_packed_sym, only: packed_index_lut, npair_from_ns
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair, &
      & transform_quadrature
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo
   use moist_model_moz_potential_lj_base, only: lj_mixing_geometric
   use moist_model_moz_solver, only: moist_model_moz_solver_deflation_type
   use moist_model_moz_closure, only: moist_model_moz_closure_type, new_closure
   use moist_model_moz_potential, only: potential_sites_type
   use moist_model_moz_solvent_vv, only: solvent_vv_type, solvent_vv_parameters_type, solvent_vv_thermodynamics_type
   use moist_model_moz_solvent_vv_thermodynamics, only: new_vv_thermodynamics
   use moist_model_moz_solvent_vv_kernel, only: vv_kblock, vv_intramolecular, vv_solve_block, vv_solve_point
   use moist_model_moz_solvent_vv_equation, only: vv_equation_type
   use moist_model_moz_solvent_vv_driver, only: vv_problem_type, vv_evaluator_type, new_vv_problem
   use moist_model_moz_1d_preconditioner, only: moist_model_moz_1d_preconditioner_spectral_type
   use test_helpers, only: check_moist_error
   use test_moz_fixtures, only: new_cspce_water_solvent, new_water_solvent
   implicit none(type, external)
   private
   public :: collect_moz_solvent_vv

   !> Grid points and spacing of the test grid, bohr
   integer, parameter :: npts = 512
   real(wp), parameter :: dr = 0.05_wp
   !> Grid points and spacing of the solve tests, bohr
   integer, parameter :: solve_npts = 1024
   real(wp), parameter :: solve_dr = 0.06_wp
   !> Stage tolerance of the solve tests
   real(wp), parameter :: solve_tolerance = 1.0e-9_wp
   !> Parameter file of the round-trip test
   character(len=*), parameter :: parameter_file = "moist-test-solvent-vv-parameters.toml"

   !> SPC/E water VV data on the test grid
   type :: vv_case_type
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(packed_index_lut) :: lut
      !> Site coordinates, bohr (3, ns)
      real(wp), allocatable :: xyz(:, :)
      !> Site densities, 1/bohr^3 (ns)
      real(wp), allocatable :: rho(:)
      real(wp) :: beta = 0.0_wp
      !> Packed w(k), short- and long-range potentials and model c_s(r) (npts, npair)
      real(wp), allocatable :: wk(:, :), u_sr(:, :), ur_lr(:, :), uk_lr(:, :), cr(:, :)
   end type vv_case_type

contains

   !> Register tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_solvent_vv(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("intramolecular", check_intramolecular), &
         & new_unittest("block_solve", check_block_solve), &
         & new_unittest("map_consistency", check_map_consistency), &
         & new_unittest("jvp_differences", check_jvp_differences), &
         & new_unittest("coarse_matrix", check_coarse_matrix), &
         & new_unittest("thread_independence", check_thread_independence), &
         & new_unittest("setup_errors", check_setup_errors), &
         & new_unittest("parameters", check_parameters), &
         & new_unittest("evaluator_products", check_evaluator_products), &
         & new_unittest("damping_step", check_damping_step), &
         & new_unittest("solve_water", check_solve_water), &
         & new_unittest("solve_variants", check_solve_variants), &
         & new_unittest("closure_ramp", check_closure_ramp), &
         & new_unittest("solve_refusals", check_solve_refusals), &
         & new_unittest("thermodynamics_reduction", check_thermodynamics_reduction), &
         & new_unittest("solve_outputs", check_solve_outputs), &
         & new_unittest("reference_values", check_reference_values), &
         & new_unittest("converged_field", check_converged_field), &
         & new_unittest("solve_replacement", check_solve_replacement), &
         & new_unittest("susceptibility_at", check_susceptibility_at)]
   end subroutine collect_moz_solvent_vv

   !> Build the SPC/E water case
   !>
   !> Model c_s(r) = -a_p exp(-(r/3)^2) with a different amplitude per pair
   !>
   !> @param[out] case  VV data
   !> @param[out] err   construction failure
   !> @param[in]  grid_npts optional radial grid size
   subroutine new_vv_case(case, err, grid_npts)
      !> VV data
      type(vv_case_type), intent(out) :: case
      !> Construction failure
      type(moist_error), allocatable, intent(out) :: err
      !> Optional radial grid size
      integer, intent(in), optional :: grid_npts
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(structure_type) :: mol
      integer :: ns, ip, np

      np = npts
      if (present(grid_npts)) np = grid_npts
      call new_context(ctx, nthreads=1, verbosity=0)
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(case%rgrid, case%kgrid, np, dr, err)
      if (allocated(err)) return
      call solvent%potential%compute(case%rgrid, case%kgrid, case%u_sr, case%ur_lr, case%uk_lr, err)
      if (allocated(err)) return
      mol = solvent%structure()
      ns = mol%nat
      case%xyz = mol%xyz
      case%rho = spread(solvent%density(), 1, ns)
      case%beta = solvent%beta()
      call case%lut%init(ns)
      allocate (case%wk(np, npair_from_ns(ns)), case%cr(np, npair_from_ns(ns)))
      call vv_intramolecular(case%kgrid%r, case%xyz, case%lut, case%wk)
      do ip = 1, size(case%cr, 2)
         case%cr(:, ip) = -(0.4_wp + 0.1_wp*ip)*exp(-(case%rgrid%r/3.0_wp)**2)
      end do
   end subroutine new_vv_case

   !> Set up an equation from a case
   !>
   !> @param[in]  case      VV data
   !> @param[out] eq        equation
   !> @param[in]  nthreads  team size
   !> @param[out] err       setup failure
   subroutine setup_equation(case, eq, nthreads, err)
      !> VV data
      type(vv_case_type), intent(in) :: case
      !> Equation
      type(vv_equation_type), intent(out) :: eq
      !> Team size
      integer, intent(in) :: nthreads
      !> Setup failure
      type(moist_error), allocatable, intent(out) :: err
      eq%nthreads = nthreads
      call eq%setup(case%rgrid, case%kgrid, case%wk, case%rho, case%beta, case%uk_lr, err)
   end subroutine setup_equation

   !> Deterministic smooth direction, packed (npts, npair)
   !>
   !> @param[in]  r     radial nodes, bohr
   !> @param[in]  seed  pattern selector
   !> @param[out] v     direction
   subroutine direction(r, seed, v)
      !> Radial nodes, bohr
      real(wp), intent(in) :: r(:)
      !> Pattern selector
      integer, intent(in) :: seed
      !> Direction
      real(wp), intent(out) :: v(:, :)
      integer :: ip
      do ip = 1, size(v, 2)
         v(:, ip) = sin(0.7_wp*(ip + seed)*r + 0.3_wp*seed)*exp(-r/4.0_wp)
      end do
   end subroutine direction

   !> Check w(k): one on the diagonal, sin(kd)/(kd) between sites
   subroutine check_intramolecular(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(vv_case_type) :: case
      real(wp) :: d, x
      integer :: ik

      call new_vv_case(case, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(case%wk(:, case%lut%idx(2, 2)) == 1.0_wp), more="w_HH(k) is not one")
      if (allocated(error)) return
      d = norm2(case%xyz(:, 1) - case%xyz(:, 3))
      do ik = 1, npts, 97
         x = d*case%kgrid%r(ik)
         call check(error, case%wk(ik, case%lut%idx(1, 3)), sin(x)/x, thr=1.0e-15_wp)
         if (allocated(error)) return
         call check(error, case%wk(ik, case%lut%idx(3, 1)), case%wk(ik, case%lut%idx(1, 3)))
         if (allocated(error)) return
      end do
   end subroutine check_intramolecular

   !> Check the blocked OZ solve against the pivoted point solve
   !>
   !> - Lanes without interchanges agree with LAPACK
   !> - A lane whose first column needs a swap is flagged; the point solve handles it
   subroutine check_block_solve(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ns = 3, swap_lane = 5
      real(wp), parameter :: rho(ns) = 0.05_wp
      real(wp) :: w(vv_kblock, ns, ns), ck(vv_kblock, ns, ns)
      real(wp) :: hk(vv_kblock, ns, ns), hk_pt(ns, ns), a(ns, ns)
      logical :: need_pivot(vv_kblock)
      integer :: kk, i, j, info

      do kk = 1, vv_kblock
         do j = 1, ns
            do i = 1, ns
               w(kk, i, j) = merge(1.0_wp, 0.0_wp, i == j) + 0.1_wp*sin(real(kk + i + j, wp))
               ck(kk, i, j) = 2.0_wp*cos(real(kk*(i + j), wp))
            end do
         end do
         w(kk, :, :) = 0.5_wp*(w(kk, :, :) + transpose(w(kk, :, :)))
         ck(kk, :, :) = 0.5_wp*(ck(kk, :, :) + transpose(ck(kk, :, :)))
      end do
      ! w = I and c = (I - A)/rho give A = [[0.1, 1, 0], [1, 0.1, 0], [0, 0, 1]]
      a = reshape([0.1_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [ns, ns])
      w(swap_lane, :, :) = 0.0_wp
      do i = 1, ns
         w(swap_lane, i, i) = 1.0_wp
         ck(swap_lane, i, :) = -a(i, :)/rho(1)
         ck(swap_lane, i, i) = ck(swap_lane, i, i) + 1.0_wp/rho(1)
      end do

      call vv_solve_block(ns, w, ck, rho, hk, need_pivot)
      call check(error, need_pivot(swap_lane), more="a lane needing a row swap was not flagged")
      if (allocated(error)) return
      call check(error, count(need_pivot), 1, more="lanes without swaps were flagged")
      if (allocated(error)) return
      do kk = 1, vv_kblock
         call vv_solve_point(ns, w(kk, :, :), ck(kk, :, :), rho, hk_pt, info)
         call check(error, info, 0)
         if (allocated(error)) return
         if (need_pivot(kk)) then
            ! A h = w c w with w = I for the swap lane
            call check(error, maxval(abs(matmul(a, hk_pt) - ck(kk, :, :))), 0.0_wp, &
               & thr=1.0e-12_wp*maxval(abs(ck(kk, :, :))), more="pivoted point solve is wrong")
         else
            call check(error, maxval(abs(hk(kk, :, :) - hk_pt)), 0.0_wp, thr=1.0e-12_wp*maxval(abs(hk_pt)), &
               & more="blocked and pivoted solves differ")
         end if
         if (allocated(error)) return
      end do
   end subroutine check_block_solve

   !> Check the map against point solves and the Ng relation
   !>
   !> - h(k) at selected points equals the pivoted point solve of c(k) = c_s(k) - beta u_lr(k)
   !> - t_s(r) = k2r(h - c - beta u_lr)
   subroutine check_map_consistency(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(vv_case_type) :: case
      type(vv_equation_type) :: eq
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      real(wp), allocatable :: tr(:, :), hk(:, :), ck(:, :), tk(:, :), tref(:, :), w_pt(:, :), c_pt(:, :)
      real(wp), allocatable :: h_pt(:, :)
      integer, parameter :: points(5) = [1, 2, 3, 200, npts]
      integer :: ip, ik, i, j, n, info

      call new_vv_case(case, err)
      if (.not. allocated(err)) call setup_equation(case, eq, 2, err)
      if (.not. allocated(err)) call new_radial_trafo(trafo, case%rgrid, case%kgrid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      n = size(case%cr, 2)
      allocate (tr(npts, n), hk(npts, n), ck(npts, n), tk(npts, n), tref(npts, n))
      call eq%map(case%cr, tr, err, hk=hk)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do ip = 1, n
         call trafo%fbt_r2k(case%cr(:, ip), ck(:, ip), err)
         if (allocated(err)) exit
         ck(:, ip) = ck(:, ip) - case%beta*case%uk_lr(:, ip)
         tk(:, ip) = hk(:, ip) - ck(:, ip) - case%beta*case%uk_lr(:, ip)
         call trafo%fbt_k2r(tk(:, ip), tref(:, ip), err)
         if (allocated(err)) exit
      end do
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(tr - tref)), 0.0_wp, thr=1.0e-12_wp*maxval(abs(tref)), &
         & more="t_s(r) is not k2r(h - c - beta u_lr)")
      if (allocated(error)) return

      allocate (w_pt(3, 3), c_pt(3, 3), h_pt(3, 3))
      do i = 1, size(points)
         ik = points(i)
         do j = 1, 3
            do ip = 1, 3
               w_pt(ip, j) = case%wk(ik, case%lut%idx(ip, j))
               c_pt(ip, j) = ck(ik, case%lut%idx(ip, j))
            end do
         end do
         call vv_solve_point(3, w_pt, c_pt, case%rho, h_pt, info)
         call check(error, info, 0)
         if (allocated(error)) return
         do j = 1, 3
            do ip = 1, j
               call check(error, hk(ik, case%lut%idx(ip, j)), h_pt(ip, j), thr=1.0e-9_wp*maxval(abs(h_pt)), &
                  & more="map h(k) differs from the point solve")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine check_map_consistency

   !> Check the exact Jacobian product against central differences of the map
   !>
   !> Relative error below 1e-6
   subroutine check_jvp_differences(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(vv_case_type) :: case
      type(vv_equation_type) :: eq
      real(wp), allocatable :: v(:, :), jv(:, :), tp(:, :), tm(:, :), fd(:, :)
      real(wp) :: eps

      call new_vv_case(case, err)
      if (.not. allocated(err)) call setup_equation(case, eq, 2, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (v, jv, tp, tm, mold=case%cr)
      call direction(case%rgrid%r, 1, v)
      eps = 1.0e-5_wp*norm2(case%cr)/norm2(v)
      call eq%linearize(case%cr, err)
      if (.not. allocated(err)) call eq%jvp(v, jv, err)
      if (.not. allocated(err)) call eq%map(case%cr + eps*v, tp, err)
      if (.not. allocated(err)) call eq%map(case%cr - eps*v, tm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      fd = (tp - tm)/(2.0_wp*eps)
      call check(error, norm2(jv - fd)/norm2(fd) < 1.0e-6_wp, &
         & more="Jacobian product differs from central differences")
   end subroutine check_jvp_differences

   !> Check the analytic coarse matrix against Jacobian products
   !>
   !> E_(iq,jp) = sum_r Omega g_i (z - T'(c' z))_q with z = g_j e_p
   subroutine check_coarse_matrix(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(vv_case_type) :: case
      type(vv_equation_type) :: eq
      integer, parameter :: nf = 2
      real(wp), parameter :: sigma = 2.35_wp
      real(wp), allocatable :: profile(:, :), wprofile(:, :), dcdt(:, :), coarse(:, :), ref(:, :)
      real(wp), allocatable :: z(:, :), jz(:, :)
      integer :: n, i, j, p, q

      call new_vv_case(case, err)
      if (.not. allocated(err)) call setup_equation(case, eq, 2, err)
      if (.not. allocated(err)) call eq%linearize(case%cr, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      n = size(case%cr, 2)
      allocate (profile(npts, nf), wprofile(npts, nf), dcdt(npts, n), coarse(nf*n, nf*n), ref(nf*n, nf*n))
      do j = 1, nf
         profile(:, j) = exp(-case%rgrid%r**2/(2.0_wp*sigma**2))*(case%rgrid%r/sigma)**(2*(j - 1))
         profile(:, j) = profile(:, j)/norm2(profile(:, j))
         wprofile(:, j) = case%rgrid%r**2*profile(:, j)
      end do
      do p = 1, n
         dcdt(:, p) = -0.6_wp + 0.1_wp*p*exp(-case%rgrid%r/5.0_wp)
      end do
      call eq%coarse_matrix(profile, wprofile, dcdt, coarse, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      allocate (z(npts, n), jz(npts, n))
      do p = 1, n
         do j = 1, nf
            z = 0.0_wp
            z(:, p) = profile(:, j)
            call eq%jvp(dcdt*z, jz, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            jz = z - jz
            do q = 1, n
               do i = 1, nf
                  ref((q - 1)*nf + i, (p - 1)*nf + j) = dot_product(wprofile(:, i), jz(:, q))
               end do
            end do
         end do
      end do
      call check(error, maxval(abs(coarse - ref)), 0.0_wp, thr=1.0e-10_wp*maxval(abs(ref)), &
         & more="analytic coarse matrix differs from Jacobian products")
   end subroutine check_coarse_matrix

   !> Check that map, Jacobian product and coarse matrix do not depend on the team size
   subroutine check_thread_independence(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(vv_case_type) :: case
      type(vv_equation_type) :: eq
      real(wp), allocatable :: v(:, :)
      real(wp), allocatable :: tr1(:, :), tr4(:, :), jv1(:, :), jv4(:, :), profile(:, :), wprofile(:, :)
      real(wp), allocatable :: e1(:, :), e4(:, :)
      logical :: in_parallel
      integer :: n, it, nt, igrid

      in_parallel = .false.
      !$ in_parallel = omp_in_parallel()
      if (in_parallel) then
         call skip_test(error, "Thread test must run as a selected test outside the suite OpenMP team")
         return
      end if
      do igrid = 1, 2
         call new_vv_case(case, err, merge(npts, 16384, igrid == 1))
         call check_moist_error(error, err)
         if (allocated(error)) return
         n = size(case%cr, 2)
         allocate (v, mold=case%cr)
         call direction(case%rgrid%r, 3, v)
         allocate (profile(size(case%cr, 1), 1), wprofile(size(case%cr, 1), 1))
         profile(:, 1) = exp(-case%rgrid%r**2/(2.0_wp*2.35_wp**2))
         wprofile(:, 1) = case%rgrid%r**2*profile(:, 1)
         do it = 1, 2
            nt = merge(1, 4, it == 1)
            call setup_equation(case, eq, nt, err)
            if (it == 1) then
               allocate (tr1, jv1, mold=case%cr)
               allocate (e1(n, n))
               if (.not. allocated(err)) call eq%map(case%cr, tr1, err)
               if (.not. allocated(err)) call eq%linearize(case%cr, err)
               if (.not. allocated(err)) call eq%jvp(v, jv1, err)
               if (.not. allocated(err)) call eq%coarse_matrix(profile, wprofile, case%cr, e1, err)
            else
               allocate (tr4, jv4, mold=case%cr)
               allocate (e4(n, n))
               if (.not. allocated(err)) call eq%map(case%cr, tr4, err)
               if (.not. allocated(err)) call eq%linearize(case%cr, err)
               if (.not. allocated(err)) call eq%jvp(v, jv4, err)
               if (.not. allocated(err)) call eq%coarse_matrix(profile, wprofile, case%cr, e4, err)
            end if
            call check_moist_error(error, err)
            if (allocated(error)) return
         end do
         call check(error, all(tr1 == tr4) .and. all(jv1 == jv4) .and. all(e1 == e4), &
            & more="results depend on the team size")
         if (allocated(error)) return
         deallocate (v, tr1, tr4, jv1, jv4, profile, wprofile, e1, e4)
      end do
   end subroutine check_thread_independence

   !> Check refusals: non-DST-IV pair, wrong shapes, map before setup, product before linearization
   subroutine check_setup_errors(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(vv_case_type) :: case
      type(vv_equation_type) :: eq
      real(wp), allocatable :: tr(:, :)

      call new_vv_case(case, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (tr, mold=case%cr)
      call eq%map(case%cr, tr, err)
      call check(error, allocated(err), more="map before setup succeeded")
      if (allocated(error)) return
      deallocate (err)

      call eq%setup(case%rgrid, case%kgrid, case%wk(:, 1:5), case%rho, case%beta, case%uk_lr, err)
      call check(error, allocated(err), more="a w(k) of the wrong shape was accepted")
      if (allocated(error)) return
      call check(error, .not. eq%ready(), more="a failed setup left the equation ready")
      if (allocated(error)) return
      deallocate (err)

      case%rgrid%transform = transform_quadrature
      call setup_equation(case, eq, 1, err)
      call check(error, allocated(err), more="a quadrature grid was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "DST-IV") > 0, more=err%message)
      if (allocated(error)) return
      deallocate (err)

      call new_vv_case(case, err)
      if (.not. allocated(err)) call setup_equation(case, eq, 1, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call eq%jvp(case%cr, tr, err)
      call check(error, allocated(err), more="a product before linearization succeeded")
      if (allocated(error)) return
      call check(error, index(err%message, "not linearized") > 0, more=err%message)
   end subroutine check_setup_errors

   !> Check the VV settings
   !>
   !> - Defaults: KH, JFNK with the spectral filter and deflation
   !> - TOML round trip, including the nested solver group
   !> - Unknown names and invalid ranges are refused
   subroutine check_parameters(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(solvent_vv_parameters_type) :: params, copy
      integer :: unit, stat

      call check(error, params%closure == "KH" .and. params%solver%name == "JFNK" .and. &
         & params%precond == "spectral" .and. params%solver%jfnk%deflation, more="unexpected VV defaults")
      if (allocated(error)) return
      call params%validate(err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      params%closure = "HNC"
      params%npts = 4096
      params%closure_ramp_rmse = 1.0e-4_wp
      params%solver%name = "HYBRID"
      params%precond = "kspace"
      params%deflation_nfun = 3
      params%solver%jfnk%restart = 30
      params%solver%jfnk%deflation = .false.
      call params%write_file(parameter_file, err)
      if (.not. allocated(err)) call copy%read_file(parameter_file, err)
      open (newunit=unit, file=parameter_file, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, copy%closure == "HNC" .and. copy%npts == 4096 .and. copy%closure_ramp_rmse == 1.0e-4_wp .and. &
         & copy%solver%name == "HYBRID" .and. copy%precond == "kspace" .and. copy%deflation_nfun == 3 .and. &
         & copy%solver%jfnk%restart == 30 .and. .not. copy%solver%jfnk%deflation, &
         & more="VV settings changed in the TOML round trip")
      if (allocated(error)) return
      call copy%init_defaults()
      call check(error, copy%closure == "KH" .and. copy%solver%jfnk%restart == 50 .and. copy%solver%jfnk%deflation, &
         & more="init_defaults did not restore the VV and solver defaults")
      if (allocated(error)) return

      call refused("closure", "PY")
      if (allocated(error)) return
      call refused("solver", "LBFGS")
      if (allocated(error)) return
      call refused("precond", "jacobi")
      if (allocated(error)) return
      call refused("closure_ramp_rmse", "0")
      if (allocated(error)) return
      call refused("restart", "0")

   contains

      !> Check that one invalid setting is refused
      subroutine refused(field, value)
         !> Setting and its invalid value
         character(len=*), intent(in) :: field, value
         type(solvent_vv_parameters_type) :: bad
         select case (field)
         case ("closure")
            bad%closure = value
         case ("solver")
            bad%solver%name = value
         case ("precond")
            bad%precond = value
         case ("closure_ramp_rmse")
            bad%closure_ramp_rmse = 0.0_wp
         case default
            bad%solver%jfnk%restart = 0
         end select
         call bad%validate(err)
         call check(error, allocated(err), more="an invalid "//field//" was accepted: "//value)
      end subroutine refused

   end subroutine check_parameters

   !> Build a VV problem from the case
   !>
   !> @param[in]  ctx      run context, borrowed by the problem
   !> @param[in]  case     VV data
   !> @param[in]  closure  closure name
   !> @param[out] problem  problem
   !> @param[out] err      setup failure
   subroutine setup_problem(ctx, case, closure, problem, err)
      !> Run context, borrowed by the problem
      type(moist_context_type), intent(in), target :: ctx
      !> VV data
      type(vv_case_type), intent(in) :: case
      !> Closure name
      character(len=*), intent(in) :: closure
      !> Problem
      type(vv_problem_type), intent(out) :: problem
      !> Setup failure
      type(moist_error), allocatable, intent(out) :: err
      type(solvent_vv_parameters_type) :: params

      params%closure = closure
      call new_vv_problem(problem, ctx, case%rgrid, case%kgrid, case%wk, case%u_sr, case%ur_lr, case%uk_lr, &
         & case%rho, case%beta, params, err)
   end subroutine setup_problem

   !> Check the evaluator against differences of its residual
   !>
   !> - J v against central differences of R, relative error below 1e-6 (HNC)
   !> - <u, J v> = <J^T u, v> to roundoff
   !> - The coarse matrix of the Gaussian space is factorized after a linearization
   subroutine check_evaluator_products(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(vv_case_type) :: case
      type(vv_problem_type), target :: problem
      type(vv_evaluator_type) :: evaluator
      type(moist_model_moz_solver_deflation_type) :: deflation
      type(moist_context_type), target :: ctx
      real(wp), allocatable :: gamma(:, :), v(:, :), u(:, :), rp(:, :), rm(:, :), jv(:), jtu(:)
      real(wp) :: eps
      logical :: feasible

      call new_context(ctx, nthreads=2, verbosity=0)
      call new_vv_case(case, err)
      if (.not. allocated(err)) call setup_problem(ctx, case, "HNC", problem, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      evaluator%problem => problem
      allocate (gamma, v, u, rp, rm, mold=case%cr)
      allocate (jv(size(gamma)), jtu(size(gamma)))
      call problem%equation%map(case%cr, gamma, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call direction(case%rgrid%r, 1, v)
      call direction(case%rgrid%r, 4, u)
      eps = 1.0e-5_wp*norm2(gamma)/norm2(v)

      jv = 1.0_wp
      call evaluator%jvp(reshape(v, [size(v)]), jv, feasible, err)
      call check(error, allocated(err) .and. .not. feasible .and. all(jv == 0.0_wp), &
         & more="a rejected Jacobian product must return a zero field")
      if (allocated(error)) return
      deallocate (err)
      jtu = 1.0_wp
      call evaluator%adjoint(reshape(u, [size(u)]), jtu, feasible, err)
      call check(error, allocated(err) .and. .not. feasible .and. all(jtu == 0.0_wp), &
         & more="a rejected transpose product must return a zero field")
      if (allocated(error)) return
      deallocate (err)
      call evaluator%linearize(gamma, feasible, err)
      if (.not. allocated(err)) call evaluator%jvp(reshape(v, [size(v)]), jv, feasible, err)
      if (.not. allocated(err)) call evaluator%adjoint(reshape(u, [size(u)]), jtu, feasible, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call evaluator%residual(gamma + eps*v, rp, feasible, err)
      if (.not. allocated(err)) call evaluator%residual(gamma - eps*v, rm, feasible, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      rp = (rp - rm)/(2.0_wp*eps)
      call check(error, norm2(jv - reshape(rp, [size(rp)]))/norm2(rp) < 1.0e-6_wp, &
         & more="residual Jacobian product differs from central differences")
      if (allocated(error)) return
      call check(error, abs(dot_product(reshape(u, [size(u)]), jv) - dot_product(jtu, reshape(v, [size(v)]))) < &
         & 1.0e-12_wp*norm2(u)*norm2(jv), more="transpose product fails the dot-product test")
      if (allocated(error)) return

      call evaluator%deflation_setup(deflation, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, deflation%active .and. deflation%kz == 2*size(gamma, 2), &
         & more="the coarse matrix was not factorized")
      call ctx%delete()
   end subroutine check_evaluator_products

   !> Check one damped Picard step of the driver
   !>
   !> From c_s = 0: t_1 = t_0 + K2R(f(k) R2K(-a R(t_0))) with t_0 = T(0), the
   !> Picard fraction a and the small-k damping filter f
   subroutine check_damping_step(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(vv_case_type) :: case
      type(vv_problem_type), target :: problem
      type(vv_evaluator_type) :: evaluator
      type(solvent_vv_parameters_type) :: params
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_model_moz_1d_preconditioner_spectral_type) :: filter
      real(wp), allocatable :: t0(:, :), residual(:, :), step(:)
      real(wp) :: rmse
      integer :: niter
      logical :: converged, feasible
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=1, verbosity=0)
      call new_vv_case(case, err)
      params%solver%name = "PICARD"
      if (.not. allocated(err)) call new_vv_problem(problem, ctx, case%rgrid, case%kgrid, case%wk, case%u_sr, &
         & case%ur_lr, case%uk_lr, case%rho, case%beta, params, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call problem%solve(params%tolerance, 1, converged, rmse, niter, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. converged .and. niter == 1, more="one Picard step converged the water solve")
      if (allocated(error)) return

      allocate (t0, residual, mold=case%cr)
      call problem%equation%map(0.0_wp*case%cr, t0, err)
      evaluator%problem => problem
      if (.not. allocated(err)) call evaluator%residual(t0, residual, feasible, err)
      rgrid = case%rgrid
      kgrid = case%kgrid
      if (.not. allocated(err)) call filter%setup(rgrid, kgrid, params%damp_small_k_k0, &
         & params%damp_small_k_min_weight, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      step = reshape(-params%solver%picard%damping*residual, [size(residual)])
      call filter%apply(step)
      t0 = t0 + reshape(step, shape(t0))
      call check(error, maxval(abs(problem%tr - t0)), 0.0_wp, thr=1.0e-12_wp*maxval(abs(t0)), &
         & more="driver step is not the damped Picard step")
      call ctx%delete()
   end subroutine check_damping_step

   !> Solve SPC/E water on the coarse solve grid
   !>
   !> @param[in]  ctx      run context
   !> @param[in]  params   settings; the grid and tolerance are set here
   !> @param[out] solvent  solved solvent
   !> @param[out] err      setup or solve failure
   subroutine solve_water(ctx, params, solvent, err)
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Settings
      type(solvent_vv_parameters_type), intent(in) :: params
      !> Solvent
      type(solvent_vv_type), intent(out) :: solvent
      !> Setup or solve failure
      type(moist_error), allocatable, intent(out) :: err
      type(solvent_vv_parameters_type) :: p

      p = params
      p%npts = solve_npts
      p%dr = solve_dr
      p%tolerance = solve_tolerance
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call solvent%solve(err, p)
   end subroutine solve_water

   !> Check a default JFNK solve of SPC/E water
   !>
   !> - Solved, with h = c + t and 1 + h >= 0
   !> - O-O peak between 2.5 and 3.1 Angstrom, above 2
   !> - The stored physical fields give a residual below the tolerance after the Ng shift
   subroutine check_solve_water(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_parameters_type) :: params
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(vv_problem_type), target :: problem
      type(vv_evaluator_type) :: evaluator
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :), u_sr(:, :), ur_lr(:, :), uk_lr(:, :), wk(:, :)
      real(wp), allocatable :: residual(:, :)
      type(packed_index_lut) :: lut
      type(structure_type) :: mol
      real(wp) :: rpeak
      logical :: feasible
      integer :: oo

      call new_context(ctx, nthreads=0, verbosity=0)
      call solve_water(ctx, params, solvent, err)
      if (.not. allocated(err)) call solvent%radial_grid(rgrid, kgrid, err)
      if (.not. allocated(err)) call solvent%correlations(cr, tr, hr, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%is_solved(), more="a converged solve left the solvent unsolved")
      if (allocated(error)) return
      call check(error, maxval(abs(hr - cr - tr)) <= 1.0e-12_wp*maxval(abs(hr)), more="h is not c + t")
      if (allocated(error)) return
      call check(error, minval(hr) >= -1.0_wp - 1.0e-10_wp, more="1 + h is negative")
      if (allocated(error)) return
      call lut%init(3)
      oo = lut%idx(1, 1)
      rpeak = rgrid%r(maxloc(hr(:, oo), 1))/aatoau
      call check(error, rpeak > 2.5_wp .and. rpeak < 3.1_wp .and. maxval(hr(:, oo)) > 1.0_wp, &
         & more="unexpected O-O peak")
      if (allocated(error)) return

      call solvent%potential%compute(rgrid, kgrid, u_sr, ur_lr, uk_lr, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = solvent%structure()
      allocate (wk, residual, mold=cr)
      call vv_intramolecular(kgrid%r, mol%xyz, lut, wk)
      call new_vv_problem(problem, ctx, rgrid, kgrid, wk, u_sr, ur_lr, uk_lr, spread(solvent%density(), 1, 3), &
         & solvent%beta(), params, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      evaluator%problem => problem
      call evaluator%residual(tr - solvent%beta()*ur_lr, residual, feasible, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, norm2(residual)/sqrt(real(size(residual), wp)) < solve_tolerance, &
         & more="stored solution does not satisfy the VV equations")
   end subroutine check_solve_water

   !> Check that other solver paths reach the default solution
   !>
   !> HYBRID switching to JFNK at once (MDIIS diverges from the cold water
   !> start), JFNK with the kspace inverse, JFNK without preconditioner or
   !> deflation; h agrees to 1e-6
   !>
   !> - Only kspace warns: water has too few site pairs for its block solve
   subroutine check_solve_variants(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_parameters_type) :: params
      real(wp), allocatable :: cr(:, :), tr(:, :), href(:, :), hr(:, :)
      character(len=*), parameter :: label(3) = [character(len=16) :: "HYBRID", "kspace", "plain JFNK"]
      integer :: ivar

      call new_context(ctx, nthreads=0, verbosity=0)
      call solve_water(ctx, params, solvent, err)
      if (.not. allocated(err)) call solvent%correlations(cr, tr, href, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do ivar = 1, 3
         call params%init_defaults()
         select case (ivar)
         case (1)
            params%solver%name = "HYBRID"
            params%solver%hybrid%switch_rmse = 10.0_wp
         case (2)
            params%precond = "kspace"
         case default
            params%precond = "none"
            params%solver%jfnk%deflation = .false.
         end select
         call solve_water(ctx, params, solvent, err)
         if (.not. allocated(err)) call solvent%correlations(cr, tr, hr, err)
         if (allocated(err)) then
            call test_failed(error, trim(label(ivar))//": "//err%message)
            return
         end if
         call check(error, maxval(abs(hr - href)), 0.0_wp, thr=1.0e-6_wp, &
            & more=trim(label(ivar))//" reached a different solution")
         if (allocated(error)) return
      end do
      call ctx%delete()
   end subroutine check_solve_variants

   !> Check that the closure ramp reaches the solution of the requested closure
   !>
   !> - HNC from KH through PSE2-4 and HNC from c_s = 0 agree in h to 1e-6
   subroutine check_closure_ramp(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_parameters_type) :: params
      real(wp), allocatable :: cr(:, :), tr(:, :), href(:, :), hr(:, :)

      call new_context(ctx, nthreads=0, verbosity=0)
      params%closure = "HNC"
      params%closure_ramp = .false.
      call solve_water(ctx, params, solvent, err)
      if (.not. allocated(err)) call solvent%correlations(cr, tr, href, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      params%closure_ramp = .true.
      call solve_water(ctx, params, solvent, err)
      if (.not. allocated(err)) call solvent%correlations(cr, tr, hr, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(hr - href)), 0.0_wp, thr=1.0e-6_wp, &
         & more="the ramped HNC solve reached a different solution")
      call ctx%delete()
   end subroutine check_closure_ramp

   !> Check refused and stale solves
   !>
   !> - Invalid settings and a stage that does not converge leave the solvent unsolved
   !> - A potential change makes a solution stale; update and solve again
   subroutine check_solve_refusals(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_parameters_type) :: params
      type(solvent_vv_thermodynamics_type) :: thermo
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :)

      call new_context(ctx, nthreads=0, verbosity=0)
      params%precond = "jacobi"
      call solve_water(ctx, params, solvent, err)
      call check(error, allocated(err) .and. .not. solvent%is_solved(), more="invalid settings were accepted")
      if (allocated(error)) return
      deallocate (err)
      call params%init_defaults()
      params%maxiter = 1
      call solve_water(ctx, params, solvent, err)
      call check(error, allocated(err), more="one iteration converged the water solve")
      if (allocated(error)) return
      call check(error, index(err%message, "did not converge") > 0 .and. .not. solvent%is_solved(), &
         & more=err%message)
      if (allocated(error)) return
      deallocate (err)

      call params%init_defaults()
      call solve_water(ctx, params, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call solvent%potential%set_mixing(lj_mixing_geometric, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. solvent%is_solved(), more="a mixing change kept the solution")
      if (allocated(error)) return
      call solvent%correlations(cr, tr, hr, err)
      call check(error, allocated(err), more="a stale solution was returned")
      if (allocated(error)) return
      deallocate (err)
      call solvent%susceptibility(cr, err)
      call check(error, allocated(err), more="a stale susceptibility was returned")
      if (allocated(error)) return
      deallocate (err)
      call solvent%thermodynamics(thermo, err)
      call check(error, allocated(err), more="stale thermodynamics were returned")
      if (allocated(error)) return
      deallocate (err)
      call solvent%update(err)
      if (.not. allocated(err)) call solvent%solve(err, solve_settings())
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%is_solved(), more="a second solve after the update failed")

   contains

      !> Default settings on the solve grid
      function solve_settings() result(p)
         type(solvent_vv_parameters_type) :: p
         p%npts = solve_npts
         p%dr = solve_dr
         p%tolerance = solve_tolerance
      end function solve_settings

   end subroutine check_solve_refusals

   !> Check the thermodynamic reductions on Gaussian correlations
   !>
   !> HNC with h = 0 and c_p = -a_p exp(-r^2) per packed pair of two sites:
   !> Phi = -c, so mu = (rho/beta) pi^(3/2) A and C = -pi^(3/2) A with
   !> A = a_11 + a_22 + 2 a_12; kappa_T = beta/(rho (1 - rho C)) and
   !> P = (3/2 rho - rho^2 C/2)/beta
   subroutine check_thermodynamics_reduction(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(packed_index_lut) :: lut
      class(moist_model_moz_closure_type), allocatable :: closure
      type(solvent_vv_thermodynamics_type) :: thermo
      real(wp), parameter :: rho = 0.005_wp, beta = 1000.0_wp, amp(3) = [0.3_wp, 0.5_wp, 0.7_wp]
      real(wp), allocatable :: hr(:, :), cr(:, :)
      real(wp) :: pi32, asum, moment
      integer :: ip

      call new_uniform_radial_pair(rgrid, kgrid, 2048, 0.005_wp, err)
      if (.not. allocated(err)) call new_closure(closure, "HNC", err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call lut%init(2)
      allocate (hr(2048, 3), cr(2048, 3))
      hr = 0.0_wp
      do ip = 1, 3
         cr(:, ip) = -amp(ip)*exp(-rgrid%r**2)
      end do
      asum = amp(lut%idx(1, 1)) + amp(lut%idx(2, 2)) + 2.0_wp*amp(lut%idx(1, 2))
      pi32 = acos(-1.0_wp)**1.5_wp
      moment = -pi32*asum
      call new_vv_thermodynamics(thermo, closure, rgrid, lut, hr, cr, rho, beta, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, thermo%direct_moment, moment, rel=.true., thr=1.0e-10_wp, more="direct-correlation moment")
      if (allocated(error)) return
      call check(error, thermo%chemical_potential, rho*pi32*asum/beta, rel=.true., thr=1.0e-10_wp, &
         & more="excess chemical potential")
      if (allocated(error)) return
      call check(error, thermo%compressibility, beta/(rho*(1.0_wp - rho*moment)), rel=.true., thr=1.0e-10_wp, &
         & more="isothermal compressibility")
      if (allocated(error)) return
      call check(error, thermo%pressure, (1.5_wp*rho - 0.5_wp*rho**2*moment)/beta, rel=.true., thr=1.0e-10_wp, &
         & more="pressure")
      if (allocated(error)) return
      call new_vv_thermodynamics(thermo, closure, rgrid, lut, hr(:, 1:2), cr(:, 1:2), rho, beta, err)
      call check(error, allocated(err), more="unpaired correlations were accepted")
   end subroutine check_thermodynamics_reduction

   !> Check the susceptibility and thermodynamics of a water solve
   !>
   !> - chi(k) = w(k) + rho r2k(h)(k)
   !> - Charge structure factor sum_ab q_a q_b chi_ab(k) vanishes as k -> 0
   !> - Compressibility route matches the k = 0 OZ relation rho kT kappa_T = 1 + rho h(k = 0)
   !> - P - n rho kT/2 = 1/(2 kappa_T)
   subroutine check_solve_outputs(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_parameters_type) :: params
      type(solvent_vv_thermodynamics_type) :: thermo
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(potential_sites_type) :: sites
      type(packed_index_lut) :: lut
      type(structure_type) :: mol
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :), chi(:, :), wk(:, :), hk(:, :), szz(:)
      real(wp) :: rho, h0
      integer :: a, b

      call new_context(ctx, nthreads=0, verbosity=0)
      call solve_water(ctx, params, solvent, err)
      if (.not. allocated(err)) call solvent%radial_grid(rgrid, kgrid, err)
      if (.not. allocated(err)) call solvent%correlations(cr, tr, hr, err)
      if (.not. allocated(err)) call solvent%susceptibility(chi, err)
      if (.not. allocated(err)) call solvent%thermodynamics(thermo, err)
      if (.not. allocated(err)) call solvent%potential%sites(sites, err)
      if (.not. allocated(err)) call new_radial_trafo(trafo, rgrid, kgrid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      rho = solvent%density()
      mol = solvent%structure()
      call lut%init(3)
      allocate (wk, hk, mold=hr)
      call vv_intramolecular(kgrid%r, mol%xyz, lut, wk)
      call trafo%fbt_r2k_all(hr, hk, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(chi - wk - rho*hk)), 0.0_wp, thr=1.0e-12_wp*maxval(abs(chi)), &
         & more="chi is not w + rho h(k)")
      if (allocated(error)) return

      allocate (szz(size(chi, 1)))
      szz = 0.0_wp
      do b = 1, 3
         do a = 1, 3
            szz = szz + sites%q(a)*sites%q(b)*chi(:, lut%idx(a, b))
         end do
      end do
      call check(error, abs(szz(1)) < 1.0e-4_wp*maxval(abs(szz)), more="charge structure factor does not vanish at small k")
      if (allocated(error)) return

      call rgrid%integrate_field(hr(:, lut%idx(1, 1)), h0, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, rho*thermo%compressibility/solvent%beta(), 1.0_wp + rho*h0, rel=.true., thr=1.0e-4_wp, &
         & more="compressibility route disagrees with the k = 0 structure factor")
      if (allocated(error)) return
      call check(error, thermo%pressure - 1.5_wp*rho/solvent%beta(), 0.5_wp/thermo%compressibility, rel=.true., &
         & thr=1.0e-10_wp, more="pressure and compressibility are inconsistent")
      if (allocated(error)) return
      call check(error, ieee_is_finite(thermo%chemical_potential), more="chemical potential is not finite")
   end subroutine check_solve_outputs

   !> Check SPC/E water against reference values, XRISM/KH
   !>
   !> Reference: moist_dev, computed for the moist_dev to moist migration;
   !> solvent 175 with the default 8192-point, 0.015 bohr grid, GMRES
   !> without preconditioner, one stage and tolerance 1e-10:
   !> rho mu_ex = -2.5671601075586147e-5 Eh/bohr^3; h(r) at three radii.
   !> The default settings here (spectral filter, deflation) reach it to 2e-11
   !> in the energy and 4e-14 in h
   subroutine check_reference_values(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_parameters_type) :: params
      type(solvent_vv_thermodynamics_type) :: thermo
      type(packed_index_lut) :: lut
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :)
      real(wp), parameter :: ref_rho_mu = -2.5671601075586147e-5_wp
      integer, parameter :: points(3) = [374, 600, 1000]
      !> Reference h_OO, h_OH and h_HH at the points
      real(wp), parameter :: ref_h(3, 3) = reshape([ &
         & 1.3504911788524510_wp, -4.5288448942201853e-2_wp, -5.7089913235408574e-2_wp, &
         & -1.1860725557315277e-1_wp, -1.3351311106949737e-1_wp, -6.5987338816384522e-3_wp, &
         & 5.3985935339539992e-4_wp, -1.2532339583426477e-2_wp, -1.6827707636984712e-2_wp], [3, 3])
      integer :: i

      call new_context(ctx, nthreads=0, verbosity=0)
      params%tolerance = 1.0e-10_wp
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call solvent%solve(err, params)
      if (.not. allocated(err)) call solvent%correlations(cr, tr, hr, err)
      if (.not. allocated(err)) call solvent%thermodynamics(thermo, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%density()*thermo%chemical_potential, ref_rho_mu, rel=.true., thr=1.0e-9_wp, &
         & more="rho mu_ex differs from the reference")
      if (allocated(error)) return
      call lut%init(3)
      do i = 1, size(points)
         call check(error, hr(points(i), lut%idx(1, 1)), ref_h(1, i), thr=1.0e-10_wp, more="h_OO differs from the reference")
         if (allocated(error)) return
         call check(error, hr(points(i), lut%idx(1, 2)), ref_h(2, i), thr=1.0e-10_wp, more="h_OH differs from the reference")
         if (allocated(error)) return
         call check(error, hr(points(i), lut%idx(2, 2)), ref_h(3, i), thr=1.0e-10_wp, more="h_HH differs from the reference")
         if (allocated(error)) return
      end do
   end subroutine check_reference_values

   !> Check that a converged stage returns the field whose residual met the tolerance
   !>
   !> One soft site and a Picard fraction of 3 (overrelaxed): the update from
   !> an input below the tolerance has a larger residual than its input
   subroutine check_converged_field(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(vv_problem_type), target :: problem
      type(vv_evaluator_type) :: evaluator
      type(solvent_vv_parameters_type) :: params
      integer, parameter :: n = 65
      real(wp), allocatable :: wk(:, :), u_sr(:, :), ur_lr(:, :), uk_lr(:, :), residual(:, :)
      real(wp) :: initial, tolerance, rmse, final
      integer :: niter
      logical :: converged, feasible
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=1, verbosity=0)
      call new_uniform_radial_pair(rgrid, kgrid, n, 0.1_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (wk(n, 1), u_sr(n, 1), ur_lr(n, 1), uk_lr(n, 1), residual(n, 1))
      wk = 1.0_wp
      u_sr(:, 1) = 0.1_wp*exp(-rgrid%r**2)
      ur_lr = 0.0_wp
      uk_lr = 0.0_wp
      params%solver%name = "PICARD"
      params%solver%picard%damping = 3.0_wp
      params%damp_small_k = .false.
      call new_vv_problem(problem, ctx, rgrid, kgrid, wk, u_sr, ur_lr, uk_lr, [0.001_wp], 1.0_wp, params, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      evaluator%problem => problem
      call problem%equation%map(problem%cr, problem%tr, err)
      if (.not. allocated(err)) call evaluator%residual(problem%tr, residual, feasible, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      initial = norm2(residual)/sqrt(real(size(residual), wp))
      tolerance = 1.1_wp*initial

      call problem%solve(tolerance, 1, converged, rmse, niter, err)
      if (.not. allocated(err)) call evaluator%residual(problem%tr, residual, feasible, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      final = norm2(residual)/sqrt(real(size(residual), wp))
      call check(error, converged, more="the stage below the tolerance did not converge")
      if (allocated(error)) return
      call check(error, final < tolerance, more="the returned field is above the tolerance")
      if (allocated(error)) return
      call check(error, final, rmse, rel=.true., thr=1.0e-12_wp, more="reported residual is not that of the returned field")
      call ctx%delete()
   end subroutine check_converged_field

   !> Check that replacing the potential makes a solution stale
   !>
   !> - An independently built potential with as many changes as the solving one
   !>   does not match its stamp: the solvent is unsolved and refuses its fields
   !> - A copy of the solving potential keeps the solution, also in a copy of the solvent
   subroutine check_solve_replacement(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent, other
      type(solvent_vv_type), allocatable :: copy
      type(solvent_vv_parameters_type) :: params
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :)

      call new_context(ctx, nthreads=0, verbosity=0)
      call solve_water(ctx, params, solvent, err)
      if (.not. allocated(err)) call new_water_solvent(ctx, other, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, other%potential%generation() /= solvent%potential%generation(), &
         & more="independently built potentials share a change stamp")
      if (allocated(error)) return

      allocate (copy, source=solvent)
      copy%potential = solvent%potential
      call check(error, copy%is_solved(), more="a copy of the solving potential lost the solution")
      if (allocated(error)) return

      solvent%potential = other%potential
      call check(error, .not. solvent%is_solved(), more="a replaced potential kept the solution")
      if (allocated(error)) return
      call solvent%correlations(cr, tr, hr, err)
      call check(error, allocated(err), more="correlations of a replaced potential were returned")
      if (allocated(error)) return
      call check(error, copy%is_solved(), more="replacing the original's potential changed the copy")
   end subroutine check_solve_replacement

   !> Check the susceptibility at arbitrary wave numbers
   !>
   !> - The solution's own nodes return the stored chi exactly
   !> - A subset of the nodes, evaluated by the direct transform, agrees with the
   !>   FFT-based DST-IV values
   !> - Between the nodes, from k = 0 up to the band edge, chi is w(k) plus
   !>   rho h(k) with h(k) = 4 pi sum_i dr r_i^2 h_i j0(k r_i)
   !> - The charge moment sum_b q_b chi_ab vanishes at k = 0
   !> - Nonfinite, negative and out-of-band wave numbers are refused
   subroutine check_susceptibility_at(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_parameters_type) :: params
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(potential_sites_type) :: sites
      type(structure_type) :: mol
      type(packed_index_lut) :: lut
      integer, parameter :: nk = 600
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :), chi(:, :), node_chi(:, :), wk(:, :), ref(:, :), k(:)
      real(wp) :: dk, moment
      real(wp) :: bad_k(3), nonfinite(3)
      integer :: i, ik, ip, a, b, j

      call new_context(ctx, nthreads=0, verbosity=0)
      call solve_water(ctx, params, solvent, err)
      if (.not. allocated(err)) call solvent%radial_grid(rgrid, kgrid, err)
      if (.not. allocated(err)) call solvent%correlations(cr, tr, hr, err)
      if (.not. allocated(err)) call solvent%susceptibility(chi, err)
      if (.not. allocated(err)) call solvent%susceptibility_at(kgrid%r, node_chi, err)
      if (.not. allocated(err)) call solvent%potential%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(node_chi == chi), more="the solution's nodes did not return the stored chi")
      if (allocated(error)) return
      call solvent%susceptibility_at(kgrid%r(1:200), node_chi, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(node_chi - chi(1:200, :))), 0.0_wp, thr=1.0e-12_wp*maxval(abs(chi)), &
         & more="direct transform at the nodes differs from the DST-IV values")
      if (allocated(error)) return

      dk = kgrid%spacing
      allocate (k(nk + 1))
      k(1:nk) = [(real(i, wp)*dk, i = 0, nk - 1)]
      k(nk + 1) = real(size(kgrid%r), wp)*dk
      call solvent%susceptibility_at(k, chi, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = solvent%structure()
      call lut%init(3)
      allocate (wk(size(k), size(hr, 2)), ref(size(k), size(hr, 2)))
      call vv_intramolecular(k, mol%xyz, lut, wk)
      do ip = 1, size(hr, 2)
         do ik = 1, size(k)
            ref(ik, ip) = 4.0_wp*acos(-1.0_wp)*sum(rgrid%w*rgrid%r**2*hr(:, ip)*j0(k(ik)*rgrid%r))
         end do
      end do
      ref = wk + solvent%density()*ref
      call check(error, maxval(abs(chi - ref)), 0.0_wp, thr=1.0e-12_wp*maxval(abs(ref)), &
         & more="chi between the nodes differs from the radial transform")
      if (allocated(error)) return

      do a = 1, 3
         moment = 0.0_wp
         do b = 1, 3
            moment = moment + sites%q(b)*chi(1, lut%idx(a, b))
         end do
         call check(error, abs(moment) < 1.0e-3_wp*maxval(abs(chi(:, lut%idx(a, a)))), &
            & more="charge moment of chi does not vanish at k = 0")
         if (allocated(error)) return
      end do

      call solvent%susceptibility_at([-0.1_wp], chi, err)
      call check(error, allocated(err), more="a negative wave number was accepted")
      if (allocated(error)) return
      deallocate (err)
      call solvent%susceptibility_at([1.01_wp*k(nk + 1)], chi, err)
      call check(error, allocated(err), more="a wave number beyond pi/dr was accepted")
      if (allocated(error)) return

      nonfinite = [ieee_value(0.0_wp, ieee_quiet_nan), ieee_value(0.0_wp, ieee_positive_inf), &
         & ieee_value(0.0_wp, ieee_negative_inf)]
      do j = 1, size(nonfinite)
         do i = 1, size(bad_k)
            bad_k = [0.0_wp, 0.5_wp*dk, dk]
            bad_k(i) = nonfinite(j)
            call solvent%susceptibility_at(bad_k, chi, err)
            call check(error, allocated(err) .and. .not. allocated(chi), &
               & more="a nonfinite wave number was accepted")
            if (allocated(error)) return
         end do
      end do
      k = kgrid%r
      k(size(k)/2) = nonfinite(1)
      call solvent%susceptibility_at(k, chi, err)
      call check(error, allocated(err) .and. .not. allocated(chi), more="a grid containing a NaN was accepted")
      if (allocated(error)) return
      call solvent%susceptibility_at(k(:0), chi, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(chi, 1) == 0, more="an empty wave-number request failed")

   contains

      !> j0(x) = sin(x)/x, 1 at x = 0
      elemental real(wp) function j0(x)
         real(wp), intent(in) :: x
         if (abs(x) < 1.0e-6_wp) then
            j0 = 1.0_wp - x*x/6.0_wp
         else
            j0 = sin(x)/x
         end if
      end function j0

   end subroutine check_susceptibility_at

end module test_moz_solvent_vv
