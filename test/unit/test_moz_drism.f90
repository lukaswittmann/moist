!> Independent dielectric limits, native derivatives and bulk-state transitions
module test_moz_drism
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_context, only: moist_context_type, new_context
   use moist_math_packed_sym, only: packed_index_lut
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_model_moz_solvent_theory, only: bulk_theory_type, new_bulk_theory
   use moist_model_moz_solvent_vv, only: solvent_vv_type, solvent_vv_parameters_type
   use moist_model_moz_solvent_vv_kernel, only: vv_intramolecular, vv_solve_point
   use moist_model_moz_solvent_vv_equation, only: vv_equation_type
   use moist_model_moz_solvent_vv_response, only: vv_response_point
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type, new_uv3d_susceptibility
   use moist_model_moz_1d_kernel, only: uv_charge_moment_correction
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use test_helpers, only: check_moist_error
   use moist_solvent_system, only: solvent_system_type
   use test_moz_fixtures, only: new_cspce_water_solvent, new_water_system
   implicit none(type, external)
   private
   public :: collect_moz_drism

   real(wp), parameter :: xyz(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, 1.43_wp, 0.0_wp, 1.10_wp, -1.43_wp, 0.0_wp, 1.10_wp], [3, 3])
   real(wp), parameter :: q(3) = [-0.8476_wp, 0.4238_wp, 0.4238_wp]
   real(wp), parameter :: rho = 0.004938_wp, beta = 1059.0_wp, smear = 0.6682_wp
contains
   subroutine collect_moz_drism(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("parameters", parameters), &
         & new_unittest("correction_reference", correction_reference), &
         & new_unittest("rotation", rotation), &
         & new_unittest("dielectric_limits", dielectric_limits), &
         & new_unittest("native_products", native_products), &
         & new_unittest("charge_moment", charge_moment), &
         & new_unittest("bulk_state", bulk_state), &
         & new_unittest("model_selection", model_selection)]
   end subroutine collect_moz_drism

   !> Default RISM, aliases, explicit epsilon = 1 and TOML registration
   subroutine parameters(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(solvent_vv_parameters_type) :: params, copy
      class(bulk_theory_type), allocatable :: theory
      real(wp) :: zeta(3, 3), scaled(3, 3)
      integer :: unit, stat
      character(len=*), parameter :: path = "moist-test-drism-parameters.toml"

      call check(error, params%theory == "RISM" .and. params%tolerance == 1.0e-11_wp)
      if (allocated(error)) return
      call new_bulk_theory(theory, params%theory, xyz, q, rho, beta, 78.44_wp, smear, err)
      if (.not. allocated(err)) call theory%correction_at(0.1_wp, zeta, scaled)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(zeta == 0.0_wp) .and. all(scaled == 0.0_wp))
      if (allocated(error)) return
      copy%theory = "XRISM"
      call check(error, params%same_state(copy), more="XRISM alias changes the RISM state")
      if (allocated(error)) return
      params%theory = "DRISM"
      params%drism_smear = 0.75_wp
      call params%write_file(path, err)
      if (.not. allocated(err)) call copy%read_file(path, err, strict=.true.)
      open(newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close(unit, status="delete")
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, params%same_state(copy))
      if (allocated(error)) return
      copy%drism_smear = 0.8_wp
      call check(error, .not. params%same_state(copy))
      if (allocated(error)) return
      call new_bulk_theory(theory, "DRISM", xyz, q, rho, beta, 1.0_wp, smear, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call theory%correction_at(0.1_wp, zeta, scaled)
      call check(error, theory%name() == "DRISM" .and. any(zeta /= 0.0_wp), &
         & more="epsilon = 1 silently selected ordinary RISM")
      if (allocated(error)) return
      copy%theory = "bad"
      call copy%validate(err)
      call check(error, allocated(err), more="unknown theory was accepted")
      if (allocated(error)) return
      deallocate(err)
      call new_bulk_theory(theory, "DRISM", xyz, q, rho, beta, ieee_value(0.0_wp, ieee_quiet_nan), smear, err)
      call check(error, allocated(err) .and. .not. allocated(theory), more="NaN dielectric was accepted")
      if (allocated(error)) return
      deallocate(err)
      call new_bulk_theory(theory, "DRISM", xyz, q, rho, beta, 0.0_wp, smear, err)
      call check(error, allocated(err) .and. .not. allocated(theory), more="dielectric below 1 was accepted")
      if (allocated(error)) return
      deallocate(err)
      call new_bulk_theory(theory, "DRISM", xyz, [0.0_wp, 0.0_wp, 0.0_wp], rho, beta, 10.0_wp, smear, err)
      call check(error, allocated(err) .and. .not. allocated(theory), more="nonpolar DRISM was accepted")
   end subroutine parameters

   !> Frozen values evaluated independently with scipy.special.spherical_jn
   !> These coordinates already have RISMiCal's canonical dipole and transverse axes.
   subroutine correction_reference(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(bulk_theory_type), allocatable :: theory
      type(packed_index_lut) :: lut
      real(wp), parameter :: k(3) = [0.25_wp, 1.0_wp, 3.0_wp]
      real(wp), parameter :: reference(6, 3) = reshape([ &
         & 3.84528899818493031e0_wp, -3.76390216900350616e0_wp, 3.68423791931281208e0_wp, &
         & -3.76390216900350616e0_wp, 3.68423791931281208e0_wp, 3.68423791931281208e0_wp, &
         & 4.71419457164198121e1_wp, -3.26401786272691865e1_wp, 2.25994333629925173e1_wp, &
         & -3.26401786272691865e1_wp, 2.25994333629925173e1_wp, 2.25994333629925173e1_wp, &
         & 4.28460458504181076e1_wp, 9.10965065491509485e0_wp, 1.93683532301464845e0_wp, &
         & 9.10965065491509485e0_wp, 1.93683532301464845e0_wp, 1.93683532301464845e0_wp], [6, 3])
      real(wp) :: zeta(3, 3), scaled(3, 3), worst
      integer :: ik, i, j
      call new_bulk_theory(theory, "DRISM", xyz, q, rho, beta, 78.44_wp, smear, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call lut%init(3)
      worst = 0.0_wp
      do ik = 1, size(k)
         call theory%correction_at(k(ik), zeta, scaled)
         do j = 1, 3
            do i = 1, j
               worst = max(worst, abs(zeta(i, j) - reference(lut%idx(i, j), ik)))
            end do
         end do
      end do
      call check(error, worst < 2.0e-12_wp, more="DRISM differs from independent Bessel reference")
   end subroutine correction_reference

   !> Arbitrary rigid rotations, translations and water hydrogen exchange
   subroutine rotation(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(bulk_theory_type), allocatable :: original, moved
      real(wp) :: zeta(3, 3), scaled(3, 3), other(3, 3), other_scaled(3, 3), matrix(3, 3)
      real(wp) :: coords(3, 3), angle, axis(3), cross(3, 3), worst
      integer :: it, i
      call new_bulk_theory(original, "DRISM", xyz, q, rho, beta, 78.44_wp, smear, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      axis = [0.3_wp, -0.7_wp, 0.4_wp]
      axis = axis/norm2(axis)
      cross = reshape([0.0_wp, axis(3), -axis(2), -axis(3), 0.0_wp, axis(1), &
         & axis(2), -axis(1), 0.0_wp], [3, 3])
      worst = 0.0_wp
      do it = 0, 12
         angle = 0.37_wp*it
         matrix = (1.0_wp - cos(angle))*spread(axis, 2, 3)*spread(axis, 1, 3) + sin(angle)*cross
         do i = 1, 3
            matrix(i, i) = matrix(i, i) + cos(angle)
         end do
         coords = matmul(matrix, xyz) + spread([1.1_wp, -0.3_wp, 2.0_wp], 2, 3)
         call new_bulk_theory(moved, "DRISM", coords, q, rho, beta, 78.44_wp, smear, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call original%correction_at(1.1_wp, zeta, scaled)
         call moved%correction_at(1.1_wp, other, other_scaled)
         worst = max(worst, maxval(abs(zeta - other)))
      end do
      coords = xyz(:, [1, 3, 2])
      call new_bulk_theory(moved, "DRISM", coords, q, rho, beta, 78.44_wp, smear, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call moved%correction_at(1.1_wp, other, other_scaled)
      worst = max(worst, maxval(abs(zeta - other([1, 3, 2], [1, 3, 2]))))
      call check(error, worst < 2.0e-12_wp, more="DRISM changed under a rigid motion or hydrogen exchange")
   end subroutine rotation

   !> Dielectric screening independently fixes the low-k charge coefficient
   subroutine dielectric_limits(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(bulk_theory_type), allocatable :: theory
      type(packed_index_lut) :: lut
      real(wp), parameter :: eps(5) = [1.0_wp, 2.0_wp, 10.0_wp, 78.44_wp, 100.0_wp]
      real(wp), parameter :: k(6) = [0.0_wp, 1.0e-12_wp, 1.0e-8_wp, 1.0e-5_wp, 0.1_wp, 1.0_wp]
      real(wp) :: distance2(3, 3), cbar(3, 3), chi(3, 3), z(3), gamma, screening
      real(wp) :: zeta(3, 3), scaled(3, 3), wk(1, 6), w(3, 3), ck(3, 3), hk(3, 3), expected(3, 3)
      integer :: i, j, ie, ik, info
      do j = 1, 3
         do i = 1, 3
            distance2(i, j) = sum((xyz(:, i) - xyz(:, j))**2)
         end do
      end do
      cbar = reshape([-0.4_wp, 0.12_wp, -0.08_wp, 0.12_wp, -0.2_wp, 0.02_wp, &
         & -0.08_wp, 0.02_wp, -0.15_wp], [3, 3])
      call lut%init(3)
      do ie = 1, size(eps)
         call new_bulk_theory(theory, "DRISM", xyz, q, rho, beta, eps(ie), smear, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         do ik = 1, size(k)
            call vv_response_point(k(ik), cbar, distance2, q, rho, beta, 0.85_wp, chi, z, info, theory)
            call check(error, info, 0)
            if (allocated(error)) return
            if (ik <= 3) then
               screening = 4.0_wp*acos(-1.0_wp)*beta*rho*dot_product(q, z)
               call check(error, screening, 1.0_wp - 1.0_wp/eps(ie), thr=2.0e-13_wp)
            else if (ik >= 5) then
               call theory%correction_at(k(ik), zeta, scaled)
               call vv_intramolecular([k(ik)], xyz, lut, wk)
               do j = 1, 3
                  do i = 1, 3
                     w(i, j) = wk(1, lut%idx(i, j))
                  end do
               end do
               gamma = 4.0_wp*acos(-1.0_wp)*beta*exp(-k(ik)**2/(4.0_wp*0.85_wp**2))
               ck = cbar - gamma*spread(q, 2, 3)*spread(q, 1, 3)/k(ik)**2
               call vv_solve_point(3, w + rho*zeta, ck, spread(rho, 1, 3), hk, info)
               expected = w + rho*(hk + zeta)
               call check(error, maxval(abs(chi - expected)) < 5.0e-11_wp, &
                  & more="reconstructed susceptibility differs from physical omega + rho h")
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine dielectric_limits

   !> Physical h offset, finite-difference JVP and packed Euclidean transpose
   subroutine native_products(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      class(bulk_theory_type), allocatable :: theory
      type(packed_index_lut) :: lut
      type(vv_equation_type) :: equation
      integer, parameter :: npts = 128
      real(wp) :: wk(npts, 6), wbar(npts, 6), offset(npts, 6), lr(npts, 6), cr(npts, 6)
      real(wp) :: direction(npts, 6), seed(npts, 6), tr(npts, 6), plus(npts, 6), minus(npts, 6)
      real(wp) :: jv(npts, 6), adjoint(npts, 6), hk(npts, 6), short(npts, 6), ck(npts, 6)
      real(wp) :: zeta(3, 3), scaled(3, 3), w(3, 3), c(3, 3), h(3, 3), lhs, rhs
      real(wp), parameter :: step = 1.0e-5_wp, local_beta = 3.0_wp, local_rho = 0.025_wp
      integer :: i, j, ik, ip, info
      call new_uniform_radial_pair(rgrid, kgrid, npts, 0.2_wp, err)
      if (.not. allocated(err)) call new_radial_trafo(trafo, rgrid, kgrid, err)
      if (.not. allocated(err)) call new_bulk_theory(theory, "DRISM", xyz, q, local_rho, local_beta, &
         & 2.0_wp, smear, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call lut%init(3)
      call vv_intramolecular(kgrid%r, xyz, lut, wk)
      do ik = 1, npts
         call theory%correction_at(kgrid%r(ik), zeta, scaled)
         do j = 1, 3
            do i = 1, j
               ip = lut%idx(i, j)
               offset(ik, ip) = zeta(i, j)
               lr(ik, ip) = 4.0_wp*acos(-1.0_wp)*q(i)*q(j)* &
                  & exp(-kgrid%r(ik)**2/(4.0_wp*0.85_wp**2))/kgrid%r(ik)**2
            end do
         end do
      end do
      wbar = wk + local_rho*offset
      do ip = 1, 6
         cr(:, ip) = -(0.2_wp + 0.03_wp*ip)*exp(-(rgrid%r/2.0_wp)**2)
         direction(:, ip) = 0.03_wp*cos(0.3_wp*ip*rgrid%r)*exp(-(rgrid%r/3.0_wp)**2)
         seed(:, ip) = sin(0.2_wp*ip*rgrid%r)*exp(-(rgrid%r/4.0_wp)**2)
      end do
      call equation%setup(rgrid, kgrid, wbar, spread(local_rho, 1, 3), local_beta, lr, err, offset)
      if (.not. allocated(err)) call equation%map(cr, tr, err, hk)
      if (.not. allocated(err)) call trafo%fbt_r2k_all(cr, short, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      ck = short - local_beta*lr
      do ik = 1, npts
         do j = 1, 3
            do i = 1, 3
               w(i, j) = wbar(ik, lut%idx(i, j))
               c(i, j) = ck(ik, lut%idx(i, j))
            end do
         end do
         call vv_solve_point(3, w, c, spread(local_rho, 1, 3), h, info)
         call check(error, info, 0)
         if (allocated(error)) return
         do j = 1, 3
            do i = 1, j
               ip = lut%idx(i, j)
               call check(error, hk(ik, ip), h(i, j) + offset(ik, ip), thr=5.0e-12_wp)
               if (allocated(error)) return
            end do
         end do
      end do
      call equation%linearize(cr, err)
      if (.not. allocated(err)) call equation%jvp(direction, jv, err)
      if (.not. allocated(err)) call equation%adjoint(seed, adjoint, err)
      if (.not. allocated(err)) call equation%map(cr + step*direction, plus, err)
      if (.not. allocated(err)) call equation%map(cr - step*direction, minus, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(jv - (plus - minus)/(2.0_wp*step))) < 2.0e-8_wp, &
         & more="DRISM JVP disagrees with the residual map")
      if (allocated(error)) return
      lhs = sum(seed*jv)
      rhs = sum(direction*adjoint)
      call check(error, lhs, rhs, thr=1.0e-11_wp*max(1.0_wp, abs(lhs), abs(rhs)))
   end subroutine native_products

   !> The analytical moment preserves susceptibility symmetry without a fit
   subroutine charge_moment(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: neutral(3, 2) = reshape([1.0_wp, 1.0_wp, 1.0_wp, &
         & 0.0_wp, 1.0_wp, -1.0_wp], [3, 2])
      real(wp) :: p(2, 3, 3), z(2, 3), target(3), k(2), before(2, 3, 3)
      logical :: applied
      integer :: a
      p = 1.0_wp
      p(1, 1, 1) = p(1, 1, 1) + 1.0e-4_wp
      k = [1.0e-3_wp, 0.1_wp]
      z(1, :) = [-0.3_wp, 0.15_wp, 0.15_wp]
      z(2, :) = z(1, :)
      before = p
      call uv_charge_moment_correction(k, q, p, applied, z)
      call check(error, applied)
      if (allocated(error)) return
      do a = 1, 3
         target(a) = dot_product(q, p(1, :, a))
      end do
      call check(error, maxval(abs(target - k(1)**2*z(1, :))) < 3.0e-16_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(p(1, :, :) - transpose(p(1, :, :)))) < 3.0e-16_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(matmul(transpose(neutral), &
         & matmul(p(1, :, :) - before(1, :, :), neutral)))) < 1.0e-14_wp, &
         & more="analytical charge-moment correction changed the charge-neutral block")
      if (allocated(error)) return
      call check(error, all(p(2, :, :) == before(2, :, :)))
   end subroutine charge_moment

   !> Real water solves; theory, smear, tolerance and table identity affect reuse
   subroutine bulk_state(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_parameters_type) :: params, changed
      type(solvent_system_type) :: system
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(uv3d_susceptibility_type) :: response
      real(wp) :: screening
      call new_context(ctx, nthreads=2, verbosity=0)
      call new_cspce_water_solvent(ctx, solvent, err)
      params%npts = 1024
      params%dr = 0.06_wp
      params%maxiter = 256
      if (.not. allocated(err)) call solvent%solve(err, params)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%theory_name() == "RISM" .and. solvent%is_solved(params))
      if (allocated(error)) return
      call new_cartesian_point_grid(grid, err, 8, 8, 8, 0.6_wp, margin=0.0_wp)
      mol = solvent%structure()
      if (.not. allocated(err)) call grid%update(mol, err)
      if (.not. allocated(err)) call new_uv3d_susceptibility(response, grid, solvent, "oz", err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      params%theory = "DRISM"
      call check(error, .not. solvent%is_solved(params))
      if (allocated(error)) return
      call solvent%solve(err, params)
      if (.not. allocated(err)) call solvent%charge_screening(screening, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, solvent%is_solved(params) .and. solvent%theory_name() == "DRISM")
      if (allocated(error)) return
      call new_water_system(system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, screening, 1.0_wp - 1.0_wp/system%solvent%epsilon, thr=2.0e-13_wp, &
         & more="DRISM did not take the solvent-system dielectric")
      if (allocated(error)) return
      call check(error, .not. response%matches(grid, solvent, "oz"), more="RISM table survived a DRISM bulk solve")
      if (allocated(error)) return
      changed = params
      changed%drism_smear = 0.75_wp
      call check(error, .not. solvent%is_solved(changed), more="changed smear reused the old bulk solution")
      if (allocated(error)) return
      changed = params
      changed%tolerance = solvent%residual()*0.5_wp
      call check(error, .not. solvent%is_solved(changed), more="stricter tolerance reused an insufficient solution")
      if (allocated(error)) return
      params%tolerance = 1.0e-8_wp
      call check(error, solvent%is_solved(params), more="looser tolerance discarded a tight bulk solution")
   end subroutine bulk_state

   !> Requested DRISM overrides a supplied RISM solution in both UV models
   subroutine model_selection(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(solvent_vv_parameters_type) :: bulk
      type(model_moz_1d_type) :: model1d
      type(model_moz_3d_type) :: model3d
      type(model_moz_1d_parameters_type) :: params1d
      type(model_moz_3d_parameters_type) :: params3d
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_math_grid_3d_cartesian_type) :: grid

      call new_context(ctx, nthreads=2, verbosity=0)
      call new_cspce_water_solvent(ctx, solvent, err)
      bulk%npts = 1024
      bulk%dr = 0.06_wp
      if (.not. allocated(err)) call solvent%solve(err, bulk)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 64, 0.2_wp, err)
      if (.not. allocated(err)) call new_cartesian_point_grid(grid, err, 8, 8, 8, 0.6_wp)
      call check_moist_error(error, err)
      if (allocated(error)) return
      params1d%solvent = bulk
      params1d%solvent%theory = "DRISM"
      params3d%solvent = params1d%solvent
      call new_moz_1d_model(model1d, ctx, rgrid, kgrid, solvent, err, params1d)
      if (.not. allocated(err)) call new_moz_3d_model(model3d, ctx, grid, solvent, err, params3d)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. model1d%solvent_solved() .and. .not. model3d%solvent_solved(), &
         & more="RISM bulk was reported as ready for requested DRISM")
      if (allocated(error)) return
      call model1d%prepare_solvent(err)
      if (.not. allocated(err)) call model3d%prepare_solvent(err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model1d%solvent_solved() .and. model3d%solvent_solved())
      if (allocated(error)) return
      call check(error, solvent%is_solved(bulk) .and. solvent%theory_name() == "RISM", &
         & more="UV setup changed the caller's bulk solvent")
      if (allocated(error)) return
      call new_moz_1d_model(model1d, ctx, rgrid, kgrid, solvent, err)
      if (.not. allocated(err)) call new_moz_3d_model(model3d, ctx, grid, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model1d%solvent_solved() .and. model3d%solvent_solved(), &
         & more="a supplied solved solvent lost its settings without an override")
   end subroutine model_selection
end module test_moz_drism
