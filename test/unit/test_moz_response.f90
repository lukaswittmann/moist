!> Short-range direct-correlation transfer and its bulk/grid cache
module test_moz_response
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, &
      & ieee_positive_inf, ieee_negative_inf
   use mctc_env, only: wp, moist_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_context, only: moist_context_type, new_context
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_packed_sym, only: packed_index_lut, npair_from_ns
   use moist_data_solvents, only: get_solvent_id, get_solvent_charges
   use moist_solvent_system, only: solvent_system_type, new_solvent_system
   use moist_model_moz_solvent_vv, only: solvent_vv_type, solvent_vv_parameters_type
   use moist_model_moz_solvent_theory, only: bulk_theory_type, new_bulk_theory
   use moist_model_moz_solvent_theory_drism, only: drism_theory_type, new_drism_theory
   use moist_model_moz_solvent_vv_response, only: vv_response_type, new_vv_response, vv_response_point
   use moist_model_moz_solvent_vv_kernel, only: vv_intramolecular, vv_solve_point
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type, new_uv3d_susceptibility
   use moist_model_moz_3d_kernel, only: uv3d_point_map
   use test_helpers, only: check_moist_error
   use test_moz_fixtures, only: new_cspce_water_solvent
   implicit none(type, external)
   private
   public :: collect_moz_response
contains

   subroutine collect_moz_response(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("cubic_accuracy", check_cubic_accuracy), &
         & new_unittest("oz_limits", check_oz_limits), &
         & new_unittest("charge_screening", check_charge_screening), &
         & new_unittest("bulk_neutrality", check_bulk_neutrality), &
         & new_unittest("bulk_response", check_bulk_response), &
         & new_unittest("cache_invalidation", check_cache_invalidation)]
   end subroutine collect_moz_response

   !> Bulk neutrality is required by both theory construction and regularized OZ
   subroutine check_bulk_neutrality(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(bulk_theory_type), allocatable :: theory
      type(drism_theory_type) :: drism
      type(vv_response_type) :: response
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: sampled(:, :), charge_response(:, :)
      real(wp), parameter :: rho = 0.025_wp, beta = 3.0_wp, alpha = 0.85_wp
      real(wp), parameter :: xyz(3, 2) = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.7_wp, 0.0_wp, 0.0_wp], [3, 2])
      real(wp), parameter :: d2(2, 2) = reshape([0.0_wp, 2.89_wp, 2.89_wp, 0.0_wp], [2, 2])
      real(wp), parameter :: valid(2, 6) = reshape([0.0_wp, 0.0_wp, -0.6_wp, 0.6_wp, &
         & -0.6_wp, 0.6_wp + epsilon(1.0_wp), -0.6_wp, 0.6_wp + 1.0e-12_wp, &
         & -0.6_wp, 0.6_wp + 1.0e-6_wp, -0.6_wp, 0.6_wp - 1.0e-6_wp], [2, 6])
      real(wp), parameter :: ks(3) = [0.0_wp, 1.0e-8_wp, 0.2_wp]
      character(len=5), parameter :: names(3) = [character(len=5) :: "RISM", "XRISM", "DRISM"]
      real(wp) :: invalid(2, 7), cr_short(32, 3), cbar(2, 2), chi(2, 2), z(2)
      integer :: i, j, ik, info

      call new_uniform_radial_pair(rgrid, kgrid, 32, 0.1_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      cr_short = 0.0_wp
      cbar = 0.0_wp
      do i = 1, size(valid, 2)
         do j = 1, size(names)
            ! DRISM requires a nonzero dipole in addition to neutrality.
            if (i == 1 .and. j == 3) cycle
            call new_bulk_theory(theory, names(j), xyz, valid(:, i), rho, beta, 10.0_wp, 0.75_wp, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, allocated(theory), more="neutral bulk theory was not constructed")
            if (allocated(error)) return
         end do
         call new_vv_response(response, rgrid, kgrid, cr_short, xyz, valid(:, i), rho, beta, alpha, 1, err)
         if (.not. allocated(err)) call response%at(ks, sampled, err, charge_response)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, all(ieee_is_finite(sampled)) .and. all(ieee_is_finite(charge_response)), &
            & more="neutral bulk reconstruction is not finite")
         if (allocated(error)) return
         if (i == 1) then
            call check(error, all(charge_response == 0.0_wp), more="uncharged bulk has a charge response")
            if (allocated(error)) return
         end if
         do ik = 1, size(ks)
            call vv_response_point(ks(ik), cbar, d2, valid(:, i), rho, beta, alpha, chi, z, info)
            call check(error, info == 0, more="bulk charges within the neutrality tolerance were rejected")
            if (allocated(error)) return
            call check(error, all(ieee_is_finite(chi)) .and. all(ieee_is_finite(z)), &
               & more="neutral bulk point reconstruction is not finite")
            if (allocated(error)) return
         end do
      end do

      call check_rounded_table("water")
      if (allocated(error)) return
      call check_rounded_table("1,1,1-trichloroethane")
      if (allocated(error)) return

      invalid(:, 1) = [-0.6_wp, 0.6_wp + 2.0e-6_wp]
      invalid(:, 2) = [-0.6_wp, 0.6_wp - 2.0e-6_wp]
      invalid(:, 3) = [-0.6_wp, 0.6_wp + 1.0e-6_wp + 1.0e-12_wp]
      invalid(:, 4) = [-0.6_wp, 0.6_wp - 1.0e-6_wp - 1.0e-12_wp]
      invalid(:, 5) = [ieee_value(0.0_wp, ieee_quiet_nan), 0.6_wp]
      invalid(:, 6) = [ieee_value(0.0_wp, ieee_positive_inf), 0.6_wp]
      invalid(:, 7) = [ieee_value(0.0_wp, ieee_negative_inf), 0.6_wp]
      do i = 1, size(invalid, 2)
         do j = 1, size(names)
            call new_bulk_theory(theory, names(j), xyz, invalid(:, i), rho, beta, 10.0_wp, 0.75_wp, err)
            call check(error, allocated(err), more="invalid bulk charges were accepted by "//trim(names(j)))
            if (allocated(error)) return
            call check(error, .not. allocated(theory), more="failed construction retained a bulk theory")
            if (allocated(error)) return
            if (i <= 4) then
               call check(error, index(err%message, "net charge") > 0, more="neutrality error omits the net charge")
               if (allocated(error)) return
               call check(error, index(err%message, "1e-6 e") > 0, more="neutrality error omits the tolerance")
               if (allocated(error)) return
            end if
         end do
         call new_drism_theory(drism, xyz, invalid(:, i), rho, beta, 10.0_wp, 0.75_wp, err)
         call check(error, allocated(err), more="direct DRISM construction accepted invalid bulk charges")
         if (allocated(error)) return
         call new_vv_response(response, rgrid, kgrid, cr_short, xyz, invalid(:, i), rho, beta, alpha, 1, err)
         call check(error, allocated(err), more="bulk response construction accepted invalid charges")
         if (allocated(error)) return
         call response%at(ks, sampled, err)
         call check(error, allocated(err) .and. .not. allocated(sampled), &
            & more="failed response construction retained a usable transform")
         if (allocated(error)) return
         do ik = 1, size(ks)
            chi = 1.0_wp
            z = 1.0_wp
            call vv_response_point(ks(ik), cbar, d2, invalid(:, i), rho, beta, alpha, chi, z, info)
            call check(error, info == -4 .and. all(chi == 0.0_wp) .and. all(z == 0.0_wp), &
               & more="direct OZ reconstruction did not reject invalid charges cleanly")
            if (allocated(error)) return
         end do
      end do
      call new_bulk_theory(theory, "RISM", xyz(:, :0), [real(wp) ::], rho, beta, 10.0_wp, 0.75_wp, err)
      call check(error, allocated(err) .and. .not. allocated(theory), more="empty bulk charge model was accepted")

   contains

      !> Rounded built-in RESP gas tables at both signs of the tolerance boundary
      subroutine check_rounded_table(name)
         character(len=*), intent(in) :: name
         type(solvent_system_type) :: system
         real(wp), allocatable :: q(:), cr(:, :), distances(:, :), direct(:, :), point(:, :), charge_point(:)
         integer :: id, ns, a, b, it, ik

         call get_solvent_id(name, id, err)
         if (.not. allocated(err)) call new_solvent_system(system, id, error=err)
         if (.not. allocated(err)) call get_solvent_charges(id, "gas", "resp", q, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, abs(sum(q)), 1.0e-6_wp, thr=1.0e-14_wp, more="rounded-table fixture changed")
         if (allocated(error)) return
         ns = size(q)
         allocate (cr(rgrid%npts, npair_from_ns(ns)), direct(ns, ns), source=0.0_wp)
         allocate (distances(ns, ns), point(ns, ns), charge_point(ns))
         do b = 1, ns
            do a = 1, ns
               distances(a, b) = sum((system%solv_mol%xyz(:, a) - system%solv_mol%xyz(:, b))**2)
            end do
         end do
         do it = 1, size(names)
            call new_bulk_theory(theory, names(it), system%solv_mol%xyz, q, rho, beta, 10.0_wp, 0.75_wp, err)
            if (.not. allocated(err)) call new_vv_response(response, rgrid, kgrid, cr, system%solv_mol%xyz, &
               & q, rho, beta, alpha, 1, err, theory=theory)
            if (.not. allocated(err)) call response%at(ks, sampled, err, charge_response)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, all(ieee_is_finite(sampled)) .and. all(ieee_is_finite(charge_response)), &
               & more=trim(names(it))//" rounded-table response is not finite: "//name)
            if (allocated(error)) return
            do ik = 1, size(ks)
               call vv_response_point(ks(ik), direct, distances, q, rho, beta, alpha, &
                  & point, charge_point, info, theory=theory)
               call check(error, info == 0 .and. all(ieee_is_finite(point)) .and. all(ieee_is_finite(charge_point)), &
                  & more=trim(names(it))//" rejected rounded-table charges: "//name)
               if (allocated(error)) return
            end do
         end do
      end subroutine check_rounded_table

   end subroutine check_bulk_neutrality

   !> A rigid dipole gives epsilon = 1 + 4 pi beta rho mu^2/3 in ordinary RISM
   subroutine check_charge_screening(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(vv_response_type) :: response
      real(wp), allocatable :: cbar(:, :)
      real(wp), parameter :: rho = 0.025_wp, beta = 3.0_wp
      real(wp), parameter :: q(2) = [-0.6_wp, 0.6_wp]
      real(wp), parameter :: xyz(3, 2) = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.7_wp, 0.0_wp, 0.0_wp], [3, 2])
      real(wp) :: screening, expected, dipole(3)

      call new_uniform_radial_pair(rgrid, kgrid, 32, 0.1_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (cbar(32, 3), source=0.0_wp)
      call new_vv_response(response, rgrid, kgrid, cbar, xyz, q, rho, beta, 0.85_wp, 1, err)
      if (.not. allocated(err)) call response%charge_screening(screening, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      dipole = matmul(xyz, q)
      expected = 4.0_wp*acos(-1.0_wp)*beta*rho*sum(dipole**2)/3.0_wp
      expected = expected/(1.0_wp + expected)
      call check(error, screening, expected, thr=1.0e-14_wp)
   end subroutine check_charge_screening

   !> An analytical Gaussian tests interpolation accuracy independently of OZ
   subroutine check_cubic_accuracy(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(vv_response_type) :: response
      real(wp), allocatable :: cr(:, :), chi(:, :), exact(:), z(:, :)
      real(wp), parameter :: k(7) = [0.0_wp, 0.037_wp, 0.155_wp, 0.271_wp, 0.923_wp, 1.717_wp, 2.01_wp]
      real(wp), parameter :: width = 0.3_wp, rho = 0.02_wp, amplitude = -0.45_wp
      real(wp) :: errors(2)
      integer :: i

      exact = 1.0_wp/(1.0_wp - rho*amplitude*(acos(-1.0_wp)/width)**1.5_wp*exp(-k*k/(4.0_wp*width)))
      do i = 1, 2
         call new_uniform_radial_pair(rgrid, kgrid, 128*4**(i - 1), 0.1_wp, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         cr = reshape(amplitude*exp(-width*rgrid%r**2), [rgrid%npts, 1])
         call new_vv_response(response, rgrid, kgrid, cr, reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), &
            & [0.0_wp], rho, 2.0_wp, 1.0_wp, 1, err)
         if (.not. allocated(err)) call response%at(k, chi, err, z)
         call check_moist_error(error, err)
         if (allocated(error)) return
         errors(i) = maxval(abs(chi(:, 1) - exact))
         call check(error, all(z == 0.0_wp), more="an uncharged solvent has a charge response")
         if (allocated(error)) return
         ! Native samples must be exact, independent of the interpolation accuracy.
         call response%at(kgrid%r, chi, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, maxval(abs(chi(:, 1) - 1.0_wp/(1.0_wp - &
            & rho*amplitude*(acos(-1.0_wp)/width)**1.5_wp*exp(-kgrid%r**2/(4.0_wp*width))))) < 1.0e-12_wp, &
            & more="native direct-correlation transform differs from the Gaussian transform")
         if (allocated(error)) return
      end do
      call check(error, errors(2) < errors(1)/32.0_wp .and. errors(2) < 2.0e-6_wp, &
         & more="short-range interpolation does not converge with cubic accuracy")
      if (allocated(error)) return
      call response%at([real(kgrid%npts, wp)*kgrid%w(1)], chi, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, abs(chi(1, 1) - 1.0_wp) < 1.0e-12_wp, more="band-edge sampling is not analytical")
      if (allocated(error)) return
      call response%at([real(wp) ::], chi, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(chi, 1) == 0, more="empty sampling did not return an empty table")
      if (allocated(error)) return
      call response%at([-0.1_wp], chi, err)
      call check(error, allocated(err), more="negative wave number was accepted")
      if (allocated(error)) return
      call response%at([ieee_value(0.0_wp, ieee_quiet_nan)], chi, err)
      call check(error, allocated(err), more="NaN wave number was accepted")
      if (allocated(error)) return
      call response%at([real(kgrid%npts, wp)*kgrid%w(1) + 0.1_wp], chi, err)
      call check(error, allocated(err), more="sampling beyond the radial band was accepted")
   end subroutine check_cubic_accuracy

   !> Compare ordinary OZ at finite k and an independent dipole limit at zero
   subroutine check_oz_limits(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: q(2) = [-0.6_wp, 0.6_wp], rho = 0.025_wp, beta = 3.0_wp, alpha = 0.85_wp
      real(wp), parameter :: d2(2, 2) = reshape([0.0_wp, 2.89_wp, 2.89_wp, 0.0_wp], [2, 2])
      real(wp), parameter :: cbar(2, 2) = reshape([-0.4_wp, 0.08_wp, 0.08_wp, -0.7_wp], [2, 2])
      real(wp), parameter :: ks(3) = [0.3_wp, 1.7_wp, 10.0_wp]
      real(wp) :: chi(2, 2), z(2), zero_chi(2, 2), zero_z(2), w(2, 2), c(2, 2), h(2, 2)
      real(wp) :: dq(2), expected_z(2), k, gamma
      complex(wp) :: cs(2), uk(2), ts(2), expected(2), amplitude
      integer :: ik, i, j, info

      do ik = 1, size(ks)
         k = ks(ik)
         do j = 1, 2
            do i = 1, 2
               w(i, j) = 1.0_wp
               if (d2(i, j) > 0.0_wp) w(i, j) = sin(k*sqrt(d2(i, j)))/(k*sqrt(d2(i, j)))
               c(i, j) = cbar(i, j) - beta*q(i)*q(j)*4.0_wp*acos(-1.0_wp)*exp(-k*k/(4*alpha**2))/k**2
            end do
         end do
         call vv_solve_point(2, w, c, [rho, rho], h, info)
         call check(error, info == 0, more="reference OZ solve failed")
         if (allocated(error)) return
         call vv_response_point(k, cbar, d2, q, rho, beta, alpha, chi, z, info)
         call check(error, info == 0, more="regularized OZ solve failed")
         if (allocated(error)) return
         call check(error, maxval(abs(chi - (w + rho*h))) < 2.0e-13_wp, &
            & more="analytical Coulomb restoration differs from ordinary OZ")
         if (allocated(error)) return
         call check(error, maxval(abs(matmul(chi, q)/k**2 - z)) < 2.0e-13_wp, &
            & more="stable charge response differs from the finite-k susceptibility")
         if (allocated(error)) return
      end do
      call vv_response_point(0.0_wp, cbar, d2, q, rho, beta, alpha, zero_chi, zero_z, info)
      call check(error, info == 0, more="zero-wavevector OZ limit failed")
      if (allocated(error)) return
      call check(error, all(zero_chi == 1.0_wp/(1.0_wp - rho*sum(cbar))), &
         & more="zero response does not match molecular compressibility")
      if (allocated(error)) return
      ! With cbar = 0, Sherman-Morrison gives Dq / (1 + rho gamma q^T Dq).
      dq = -matmul(d2, q)/6.0_wp
      gamma = 4.0_wp*acos(-1.0_wp)*beta
      expected_z = dq/(1.0_wp + rho*gamma*dot_product(q, dq))
      call vv_response_point(0.0_wp, 0.0_wp*cbar, d2, q, rho, beta, alpha, chi, z, info)
      call check(error, info == 0 .and. maxval(abs(z - expected_z)) < 1.0e-14_wp, &
         & more="zero charge response differs from the analytical dipole limit")
      if (allocated(error)) return
      k = 1.0e-12_wp
      call vv_response_point(k, cbar, d2, q, rho, beta, alpha, chi, z, info)
      call check(error, info == 0 .and. maxval(abs(z - zero_z)) < 1.0e-14_wp, &
         & more="small-k charge response loses its finite limit")
      if (allocated(error)) return
      cs = [cmplx(0.2_wp, -0.3_wp, wp), cmplx(0.7_wp, 0.1_wp, wp)]
      amplitude = cmplx(0.9_wp, -0.4_wp, wp)
      uk = amplitude*q/k**2
      call uv3d_point_map(chi, cs, uk, beta, ts, q, z, k*k)
      expected = matmul(cs, zero_chi) - cs - beta*amplitude*zero_z
      call check(error, maxval(abs(ts - expected)) < 1.0e-13_wp, &
         & more="charged-solute map amplifies low-k cancellation")
   end subroutine check_oz_limits

   !> Check the achieved default residual and OZ consistency of a physical bulk solve
   subroutine check_bulk_response(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_parameters_type) :: params
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(packed_index_lut) :: lut
      real(wp), allocatable :: cr(:, :), tr(:, :), hr(:, :), u(:, :), ur(:, :), uk(:, :), ck(:, :)
      real(wp), allocatable :: wk(:, :), chi(:, :), exported(:, :)
      real(wp) :: w(3, 3), c(3, 3), h(3, 3), tight_gap, loose_gap
      integer :: ik, i, j, info

      call new_context(ctx, nthreads=1, verbosity=0)
      params%npts = 1024
      params%dr = 0.06_wp
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call solvent%solve(err, params)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, params%tolerance == 1.0e-11_wp .and. solvent%residual() < params%tolerance, &
         & more="default bulk solve did not achieve the tight residual")
      if (allocated(error)) return
      call solvent%radial_grid(rgrid, kgrid, err)
      if (.not. allocated(err)) call solvent%correlations(cr, tr, hr, err)
      if (.not. allocated(err)) call solvent%potential%compute(rgrid, kgrid, u, ur, uk, err)
      if (.not. allocated(err)) call new_radial_trafo(trafo, rgrid, kgrid, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (ck, wk, mold=cr)
      call trafo%fbt_r2k_all(cr + solvent%beta()*ur, ck, err)
      if (.not. allocated(err)) call solvent%susceptibility_oz_at(kgrid%r, chi, err)
      if (.not. allocated(err)) call solvent%susceptibility(exported, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call lut%init(3)
      block
         use mctc_io, only: structure_type
         type(structure_type) :: mol
         mol = solvent%structure()
         call vv_intramolecular(kgrid%r, mol%xyz, lut, wk)
      end block
      do ik = 10, params%npts, 31
         do j = 1, 3
            do i = 1, 3
               w(i, j) = wk(ik, lut%idx(i, j))
               c(i, j) = ck(ik, lut%idx(i, j)) - solvent%beta()*uk(ik, lut%idx(i, j))
            end do
         end do
         call vv_solve_point(3, w, c, spread(solvent%density(), 1, 3), h, info)
         call check(error, info == 0, more="bulk reference OZ solve failed")
         if (allocated(error)) return
         do j = 1, 3
            do i = 1, j
               call check(error, chi(ik, lut%idx(i, j)), w(i, j) + solvent%density()*h(i, j), thr=2.0e-9_wp, &
                  & more="native reconstructed response is inconsistent with direct-correlation OZ")
               if (allocated(error)) return
            end do
         end do
      end do
      tight_gap = maxval(abs(chi - exported))
      call check(error, tight_gap < 1.0e-7_wp, more="tight bulk solution is inconsistent with reconstructed OZ")
      if (allocated(error)) return
      params%tolerance = 1.0e-4_wp
      call solvent%solve(err, params)
      if (.not. allocated(err)) call solvent%susceptibility_oz_at(kgrid%r, chi, err)
      if (.not. allocated(err)) call solvent%susceptibility(exported, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      loose_gap = maxval(abs(chi - exported))
      call check(error, solvent%residual() < params%tolerance .and. loose_gap > 100.0_wp*tight_gap, &
         & more="explicit loose bulk tolerance was not respected or did not expose exported-h inconsistency")
   end subroutine check_bulk_response

   !> Cache reuse and invalidation after bulk, transfer mode and target grid changes
   subroutine check_cache_invalidation(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent, copy
      type(solvent_vv_parameters_type) :: params
      type(moist_math_grid_3d_cartesian_type) :: grid, other
      type(uv3d_susceptibility_type) :: table

      call new_context(ctx, nthreads=1, verbosity=0)
      params%npts = 1024
      params%dr = 0.06_wp
      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call solvent%solve(err, params)
      if (.not. allocated(err)) call new_cartesian_point_grid(grid, err, 8, 10, 12, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call grid%update(solvent%structure(), err)
      if (.not. allocated(err)) call new_uv3d_susceptibility(table, grid, solvent, "oz", err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, table%matches(grid, solvent, "OZ"), more="current response cache was not reusable")
      if (allocated(error)) return
      copy = solvent
      call check(error, table%matches(grid, copy, "oz"), more="copying a bulk solution changed its identity")
      if (allocated(error)) return
      call check(error, .not. table%matches(grid, solvent, "direct"), more="transfer-mode change kept the cache")
      if (allocated(error)) return
      call new_cartesian_point_grid(other, err, 8, 10, 12, 0.7_wp, margin=0.0_wp)
      if (.not. allocated(err)) call other%update(solvent%structure(), err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. table%matches(other, solvent, "oz"), more="reciprocal-grid change kept the cache")
      if (allocated(error)) return
      params%tolerance = 1.0e-7_wp
      call solvent%solve(err, params)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. table%matches(grid, solvent, "oz"), more="bulk re-solve kept the old cache")
      if (allocated(error)) return
      call check(error, table%matches(grid, copy, "oz"), more="bulk re-solve invalidated an independent copy")
      if (allocated(error)) return
      call new_uv3d_susceptibility(table, grid, solvent, "oz", err)
      if (.not. allocated(err)) call solvent%potential%set_ng_split(.true., 0.9_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. table%matches(grid, solvent, "oz") .and. solvent%residual() == huge(1.0_wp), &
         & more="bulk potential change kept the response or its residual current")
   end subroutine check_cache_invalidation
end module test_moz_response
