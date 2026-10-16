!> Discrete 1D UV nuclear gradients, from kernels through the host coupling
module test_moz_1d_gradient
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_math_packed_sym, only: packed_index_lut, npair_from_ns
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_model_moz_solvent_vv_kernel, only: vv_intramolecular, vv_intramolecular_adjoint
   use moist_model_moz_1d_equation, only: uv_equation_type
   use moist_model_moz_1d_driver, only: uv_problem_type, new_uv_problem
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_model_moz_solvent_vv, only: solvent_vv_type
   use moist_model_moz_potential, only: lj_cspce, coulomb_spce, monopole_type, potential_sites_type
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type, response_item_type, atomic_charge_adjoint_response_type
   use test_helpers, only: check_moist_error
   use test_moz_fixtures, only: new_cspce_water_solvent
   implicit none(type, external)
   private
   public :: collect_moz_1d_gradient

   real(wp), parameter :: coordinates(3, 3) = reshape([ &
      & 0.2_wp, -0.3_wp, 0.1_wp, 1.4_wp, 0.2_wp, -0.4_wp, -0.5_wp, 1.1_wp, 0.7_wp], [3, 3])
   real(wp), parameter :: direction(3, 3) = reshape([ &
      & 0.3_wp, -0.2_wp, 0.1_wp, -0.1_wp, 0.4_wp, -0.3_wp, 0.2_wp, 0.1_wp, -0.4_wp], [3, 3])

contains

   !> Register kernel, equation, all-closure and complete model checks
   subroutine collect_moz_1d_gradient(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("intramolecular", check_intramolecular), &
         new_unittest("equation_geometry", check_equation_geometry), &
         new_unittest("problem_geometry", check_problem_geometry), &
         new_unittest("gradient_fixed", check_gradient_fixed), &
         new_unittest("gradient_host", check_gradient_host), &
         new_unittest("gradient_drism_fixed", check_gradient_drism_fixed), &
         new_unittest("gradient_drism_host", check_gradient_drism_host)]
   end subroutine collect_moz_1d_gradient

   !> Zero modes, coincident sites, the small-x branch, force/torque and central differences
   subroutine check_intramolecular(error)
      type(error_type), allocatable, intent(out) :: error
      type(packed_index_lut) :: lut
      real(wp), parameter :: k(6) = [0.0_wp, 1.0e-6_wp, 2.0e-4_wp, 0.004_wp, 0.8_wp, 4.2_wp]
      real(wp) :: xyz(3, 3), ww(6, 6), plus(6, 6), minus(6, 6), grad(3, 3), torque(3)
      real(wp) :: eps, fd, dot, delta(3, 3)
      integer :: i, j, ik, ip, mode

      call lut%init(3)
      do ip = 1, 6
         do ik = 1, 6
            ww(ik, ip) = cos(real(ik + 2*ip, wp))/real(ik + ip, wp)
         end do
      end do
      do mode = 1, 2
         xyz = coordinates
         if (mode == 2) xyz(:, 2) = xyz(:, 1)
         grad = 0.0_wp
         call vv_intramolecular_adjoint(k, xyz, lut, ww, grad)
         call check(error, maxval(abs(sum(grad, dim=2))) < 1.0e-14_wp, more="nonzero net intramolecular force")
         if (allocated(error)) return
         torque = 0.0_wp
         do i = 1, 3
            torque = torque + cross(xyz(:, i), grad(:, i))
         end do
         call check(error, maxval(abs(torque)) < 1.0e-14_wp, more="nonzero intramolecular torque")
         if (allocated(error)) return
         do j = 1, 3
            do i = 1, 3
               eps = 1.0e-6_wp
               delta = 0.0_wp
               delta(i, j) = eps
               call vv_intramolecular(k, xyz + delta, lut, plus)
               call vv_intramolecular(k, xyz - delta, lut, minus)
               fd = sum(ww*(plus - minus))/(2.0_wp*eps)
               call check(error, abs(fd - grad(i, j)) < 2.0e-9_wp, more="sinc coordinate pullback differs")
               if (allocated(error)) return
            end do
         end do
         dot = sum(grad*direction)
         eps = 1.0e-6_wp
         call vv_intramolecular(k, xyz + eps*direction, lut, plus)
         call vv_intramolecular(k, xyz - eps*direction, lut, minus)
         fd = sum(ww*(plus - minus))/(2.0_wp*eps)
         call check(error, abs(fd - dot) < 2.0e-9_wp, more="combined sinc pullback differs")
         if (allocated(error)) return
      end do
      ww = 0.0_wp
      ww(1, :) = 3.0_wp
      grad = 0.0_wp
      call vv_intramolecular_adjoint(k, coordinates, lut, ww, grad)
      call check(error, all(grad == 0.0_wp), more="zero mode has a coordinate derivative")
   end subroutine check_intramolecular

   !> Packed W pullback includes the analytic Coulomb tail and a nonsymmetric P
   subroutine check_equation_geometry(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(packed_index_lut) :: lut
      type(uv_equation_type) :: eq
      integer, parameter :: n = 63, nv = 2
      real(wp) :: wu(n, 6), response(n, nv, nv), uk(n, 6), cr(n, 6), seed(n, 6)
      real(wp) :: mu(n, 6), ww(n, 6), ww_bad(n, 5), plus(n, 6), minus(n, 6), grad(3, 3), ref(3, 3)
      real(wp) :: eps, dot, fd
      integer :: col, a, b, nt, is

      call new_uniform_radial_pair(rgrid, kgrid, n, 0.12_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call lut%init(3)
      call vv_intramolecular(kgrid%r, coordinates, lut, wu)
      do a = 1, nv
         do b = 1, nv
            response(:, b, a) = (0.15_wp + 0.07_wp*a + 0.03_wp*b)*exp(-kgrid%r**2)
            if (a == b) response(:, b, a) = response(:, b, a) + 0.5_wp
         end do
      end do
      do col = 1, 6
         cr(:, col) = cos(real(col, wp))*exp(-rgrid%r**2)
         seed(:, col) = sin(real(col, wp))*rgrid%r*exp(-(rgrid%r/2.0_wp)**2)
         uk(:, col) = 0.04_wp*cos(real(col, wp))*exp(-kgrid%r**2)/kgrid%r**2
      end do
      do nt = 1, 4, 3
         eq%nthreads = nt
         call eq%setup(rgrid, kgrid, wu, response, 2.3_wp, uk, err)
         if (.not. allocated(err)) call eq%adjoint(seed, mu, err, cr=cr, ww=ww)
         call check_moist_error(error, err)
         if (allocated(error)) return
         grad = 0.0_wp
         call vv_intramolecular_adjoint(kgrid%r, coordinates, lut, ww, grad)
         if (nt == 1) ref = grad
         call check(error, maxval(abs(grad - ref)) < 1.0e-13_wp, more="W pullback changes with thread count")
         if (allocated(error)) return
         dot = sum(grad*direction)
         do is = 1, 2
            eps = 1.0e-4_wp/real(is, wp)
            call map_at(coordinates + eps*direction, plus)
            if (.not. allocated(err)) call map_at(coordinates - eps*direction, minus)
            call check_moist_error(error, err)
            if (allocated(error)) return
            fd = sum(seed*(plus - minus))/(2.0_wp*eps)
            call check(error, abs(fd - dot) < 1.0e-7_wp*max(1.0e-8_wp, abs(fd)), &
               & more="equation nuclear pullback differs from finite differences")
            if (allocated(error)) return
         end do
      end do
      call eq%adjoint(seed, mu, err, ww=ww)
      call check(error, allocated(err), more="missing primal correlation was accepted")
      if (allocated(error)) return
      call eq%adjoint(seed, mu, err, cr=cr)
      call check(error, allocated(err), more="missing W accumulator was accepted")
      if (allocated(error)) return
      call eq%adjoint(seed, mu, err, cr=cr, ww=ww_bad)
      call check(error, allocated(err), more="wrong W accumulator shape was accepted")
   contains
      subroutine map_at(xyz, tr)
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(out) :: tr(:, :)
         call vv_intramolecular(kgrid%r, xyz, lut, wu)
         call eq%setup(rgrid, kgrid, wu, response, 2.3_wp, uk, err)
         if (.not. allocated(err)) call eq%map(cr, tr, err)
      end subroutine map_at
   end subroutine check_equation_geometry

   !> One adjoint for every closure and signed/zero energy seeds, at two displacement steps
   subroutine check_problem_geometry(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(packed_index_lut) :: lut
      type(uv_problem_type) :: problem
      type(model_moz_1d_parameters_type) :: params
      integer, parameter :: n = 128, nv = 2
      character(len=4), parameter :: closures(5) = ["HNC ", "KH  ", "PSE2", "PSE3", "PSE4"]
      real(wp), parameter :: seeds(3) = [-0.7_wp, 0.0_wp, 1.8_wp]
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), w_k(:, :), ww(:, :), lambda(:, :)
      real(wp) :: wu(n, 6), response(n, nv, nv), uk(n, 6), sr(n, 6), lr(n, 6), grad(3, 3)
      real(wp) :: e0, after, plus, minus, dot, fd, eps, worst
      integer :: col, a, b, ic, is, ih

      call new_uniform_radial_pair(rgrid, kgrid, n, 0.08_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call lut%init(3)
      do a = 1, nv
         do b = 1, nv
            response(:, b, a) = (0.05_wp + 0.02_wp*a + 0.01_wp*b)*exp(-kgrid%r**2)
            if (a == b) response(:, b, a) = response(:, b, a) + 0.25_wp
         end do
      end do
      do col = 1, 6
         sr(:, col) = 0.06_wp*cos(real(col, wp))*exp(-rgrid%r**2)
         lr(:, col) = 0.005_wp*sin(real(col, wp))*exp(-(rgrid%r/2.0_wp)**2)
         uk(:, col) = 0.01_wp*sin(real(col, wp))*exp(-kgrid%r**2)
      end do
      params%precond = "none"
      params%damp_small_k = .false.
      params%tolerance = 1.0e-13_wp
      worst = 0.0_wp
      do ic = 1, size(closures)
         params%closure = closures(ic)
         do is = 1, size(seeds)
            call solve_at(coordinates, e0)
            if (allocated(error)) return
            call problem%response(1.0e-11_wp, w_sr, w_lr, w_k, lambda, err, &
               & free_energy_weight=seeds(is), w_wu=ww)
            if (.not. allocated(err)) call problem%free_energy(after, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, after == e0, more="adjoint changed the accepted energy")
            if (allocated(error)) return
            grad = 0.0_wp
            call vv_intramolecular_adjoint(kgrid%r, coordinates, lut, ww, grad)
            if (seeds(is) == 0.0_wp) then
               call check(error, all(grad == 0.0_wp) .and. all(ww == 0.0_wp), more="zero energy has a gradient")
               if (allocated(error)) return
               cycle
            end if
            dot = sum(grad*direction)
            do ih = 1, 2
               eps = 1.0e-4_wp/real(ih, wp)
               call solve_at(coordinates + eps*direction, plus)
               if (allocated(error)) return
               call solve_at(coordinates - eps*direction, minus)
               if (allocated(error)) return
               fd = seeds(is)*(plus - minus)/(2.0_wp*eps)
               worst = max(worst, abs(fd - dot)/max(1.0e-8_wp, abs(fd)))
               call check(error, abs(fd - dot) < 2.0e-6_wp*max(1.0e-8_wp, abs(fd)), &
                  & more=trim(closures(ic))//" nuclear derivative differs from finite differences")
               if (allocated(error)) return
            end do
         end do
      end do
      write (*, "(a,es12.4)") "All-closure 1D geometry maximum scaled FD discrepancy: ", worst
   contains
      subroutine solve_at(xyz, energy)
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(out) :: energy
         real(wp) :: rmse
         integer :: niter
         logical :: converged
         call vv_intramolecular(kgrid%r, xyz, lut, wu)
         call new_uv_problem(problem, rgrid, kgrid, wu, response, sr, lr, uk, [0.012_wp, 0.017_wp], &
            & 2.3_wp, params, 4, 0, err)
         if (.not. allocated(err)) call problem%solve(params%tolerance, params%maxiter, converged, rmse, niter, err)
         if (.not. allocated(err)) then
            if (.not. converged) then
               call test_failed(error, "Synthetic 1D gradient solve did not converge")
               return
            end if
            call problem%free_energy(energy, err)
         end if
         call check_moist_error(error, err)
      end subroutine solve_at
   end subroutine check_problem_geometry

   subroutine check_gradient_fixed(error)
      type(error_type), allocatable, intent(out) :: error
      call model_gradient_case(error, .false.)
   end subroutine check_gradient_fixed

   subroutine check_gradient_host(error)
      type(error_type), allocatable, intent(out) :: error
      call model_gradient_case(error, .true.)
   end subroutine check_gradient_host

   subroutine check_gradient_drism_fixed(error)
      type(error_type), allocatable, intent(out) :: error
      call model_gradient_case(error, .false., drism=.true.)
   end subroutine check_gradient_drism_fixed

   subroutine check_gradient_drism_host(error)
      type(error_type), allocatable, intent(out) :: error
      call model_gradient_case(error, .true., drism=.true.)
   end subroutine check_gradient_drism_host

   !> Complete KH water-solute gradients with fixed bulk and an optional bond-dependent host charge
   subroutine model_gradient_case(error, host, drism)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: host
      logical, intent(in), optional :: drism
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(solvent_vv_type) :: solvent
      type(model_moz_1d_type), target :: model
      type(model_moz_1d_parameters_type) :: params
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(structure_type) :: mol, moved
      type(monopole_type) :: charges
      type(potential_sites_type) :: sites
      type(coupling_type), pointer :: coupling
      type(response_type) :: response
      class(response_item_type), allocatable :: item
      real(wp) :: xyz(3, 3), grad(3, 3), repeated(3, 3), total(3, 3), q(3), dg_dq(3), torque(3), dq(3)
      real(wp) :: before, after, plus, minus, eps, fd, worst, worst_absolute, bond0, bond, angle, rotation(3, 3)
      real(wp), parameter :: transfer = 0.07_wp
      integer :: i, a, is, nitem

      call new_context(ctx, nthreads=4, verbosity=0)
      call new_cspce_water_solvent(ctx, solvent, err)
      params%tolerance = 1.0e-13_wp
      params%solvent%npts = 1024
      params%solvent%dr = 0.06_wp
      params%solvent%tolerance = 1.0e-12_wp
      if (present(drism)) then
         if (drism) params%solvent%theory = "DRISM"
      end if
      if (.not. allocated(err)) call solvent%solve(err, params%solvent)
      if (.not. allocated(err)) call solvent%potential%sites(sites, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 1024, 0.06_wp, err)
      if (.not. allocated(err)) call new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, err, params)
      if (.not. allocated(err)) call model%potential%add(lj_cspce, err)
      if (host) then
         if (.not. allocated(err)) call model%potential%add(charges, err)
      else
         if (.not. allocated(err)) call model%potential%add(coulomb_spce, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) return
      mol = solvent%structure()
      xyz = mol%xyz
      bond0 = norm2(xyz(:, 2) - xyz(:, 1))
      xyz(:, 2) = xyz(:, 2) + [0.13_wp, -0.07_wp, 0.10_wp]
      xyz(:, 3) = xyz(:, 3) + [-0.08_wp, 0.05_wp, -0.03_wp]
      call energy_at(xyz, before)
      if (allocated(error)) return
      call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_gradient(coupling, err)
      if (.not. allocated(err) .and. host) call answer_host(coupling, xyz)
      grad = 0.003_wp
      if (.not. allocated(err)) call model%get_gradient(coupling, response, grad, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      grad = grad - 0.003_wp
      total = grad
      dg_dq = 0.0_wp
      nitem = 0
      do while (response%next())
         item = response%item()
         select type (item)
         type is (atomic_charge_adjoint_response_type)
            dg_dq = item%dg_dq
            nitem = nitem + 1
         class default
            call test_failed(error, "Unexpected 1D gradient host response")
            return
         end select
      end do
      call check(error, nitem == merge(1, 0, host), more="wrong 1D gradient host response count")
      if (allocated(error)) return
      if (host) then
         bond = norm2(xyz(:, 2) - xyz(:, 1))
         dq = transfer*(dg_dq(2) - dg_dq(1))*(xyz(:, 2) - xyz(:, 1))/bond
         total(:, 2) = total(:, 2) + dq
         total(:, 1) = total(:, 1) - dq
      end if
      repeated = 0.0_wp
      call model%get_gradient(coupling, response, repeated, err)
      call model%release_coupling(coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(repeated - grad)) < 1.0e-10_wp, more="repeated gradient accumulated old state")
      if (allocated(error)) return
      call current_energy(after)
      if (allocated(error)) return
      call check(error, after == before, more="gradient changed the cached model energy")
      if (allocated(error)) return
      call check(error, maxval(abs(sum(total, dim=2))) < 1.0e-10_wp, more="nonzero total nuclear force")
      if (allocated(error)) return
      torque = 0.0_wp
      do i = 1, 3
         torque = torque + cross(xyz(:, i), total(:, i))
      end do
      call check(error, maxval(abs(torque)) < 1.0e-10_wp, more="nonzero total nuclear torque")
      if (allocated(error)) return
      worst = 0.0_wp
      worst_absolute = 0.0_wp
      do is = 1, 2
         eps = 2.0e-4_wp/real(is, wp)
         do i = 1, 3
            do a = 1, 3
               moved = mol
               moved%xyz = xyz
               moved%xyz(a, i) = moved%xyz(a, i) + eps
               call energy_at(moved%xyz, plus)
               if (allocated(error)) return
               moved%xyz(a, i) = xyz(a, i) - eps
               call energy_at(moved%xyz, minus)
               if (allocated(error)) return
               fd = (plus - minus)/(2.0_wp*eps)
               worst = max(worst, abs(fd - total(a, i))/max(1.0e-6_wp, abs(fd)))
               worst_absolute = max(worst_absolute, abs(fd - total(a, i)))
               if (abs(fd - total(a, i)) >= 2.0e-5_wp*max(1.0e-6_wp, abs(fd))) then
                  write (*, "(a,3i3,4es16.7)") "1D model FD (step,atom,axis; eps,analytic,FD,error): ", &
                     & is, i, a, eps, total(a, i), fd, abs(fd - total(a, i))
               end if
               call check(error, abs(fd - total(a, i)) < 2.0e-5_wp*max(1.0e-6_wp, abs(fd)), &
                  & more="complete 1D nuclear gradient differs from finite differences")
               if (allocated(error)) return
            end do
         end do
      end do
      call energy_at(xyz + spread([0.3_wp, -0.2_wp, 0.4_wp], 2, 3), plus)
      if (allocated(error)) return
      call check(error, abs(plus - before) < 1.0e-12_wp, more="translation changed the 1D energy")
      if (allocated(error)) return
      angle = 0.37_wp
      rotation = reshape([cos(angle), sin(angle), 0.0_wp, -sin(angle), cos(angle), 0.0_wp, &
         & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      call energy_at(matmul(rotation, xyz), plus)
      if (allocated(error)) return
      call check(error, abs(plus - before) < 1.0e-12_wp, more="rotation changed the 1D energy")
      if (allocated(error)) return
      call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_gradient(coupling, err)
      if (.not. allocated(err) .and. host) call answer_host(coupling, mol%xyz)
      repeated = 0.0_wp
      if (.not. allocated(err)) call model%get_gradient(coupling, response, repeated, err)
      if (associated(coupling)) call model%release_coupling(coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(repeated - matmul(rotation, grad))) < 1.0e-10_wp, &
         & more="rotation changed the gradient covariantly")
      if (allocated(error)) return
      mol%xyz = xyz
      call model%update(mol, err)
      if (.not. allocated(err)) call model%new_coupling(coupling, err)
      if (.not. allocated(err)) call model%prepare_gradient(coupling, err)
      if (.not. allocated(err) .and. host) call answer_host(coupling, xyz)
      repeated = 0.0_wp
      if (.not. allocated(err)) call model%get_gradient(coupling, response, repeated, err)
      if (associated(coupling)) call model%release_coupling(coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(repeated - grad)) < 1.0e-10_wp, more="gradient requires a preceding energy")
      if (allocated(error)) return
      write (*, "(a,l1,a,es12.4)") "Complete 1D model, host=", host, "; maximum scaled FD discrepancy: ", worst
      write (*, "(a,es12.4)") "Maximum absolute FD discrepancy, Hartree/bohr: ", worst_absolute
   contains
      subroutine answer_host(cpl, positions)
         type(coupling_type), intent(inout) :: cpl
         real(wp), intent(in) :: positions(:, :)
         q = sites%q
         bond = norm2(positions(:, 2) - positions(:, 1))
         q(1) = q(1) - transfer*(bond - bond0)
         q(2) = q(2) + transfer*(bond - bond0)
         if (cpl%next()) call cpl%answer("q", q, err)
      end subroutine answer_host

      subroutine current_energy(energy)
         real(wp), intent(out) :: energy
         call model%new_coupling(coupling, err)
         if (.not. allocated(err)) call model%prepare_energy(coupling, err)
         if (.not. allocated(err) .and. host) call answer_host(coupling, mol%xyz)
         energy = 0.0_wp
         if (.not. allocated(err)) call model%get_energy(coupling, energy, err)
         if (associated(coupling)) call model%release_coupling(coupling)
         call check_moist_error(error, err)
      end subroutine current_energy

      subroutine energy_at(positions, energy)
         real(wp), intent(in) :: positions(:, :)
         real(wp), intent(out) :: energy
         mol%xyz = positions
         call model%update(mol, err)
         if (.not. allocated(err)) call current_energy(energy)
         if (allocated(err)) call check_moist_error(error, err)
      end subroutine energy_at
   end subroutine model_gradient_case

   pure function cross(a, b) result(c)
      real(wp), intent(in) :: a(3), b(3)
      real(wp) :: c(3)
      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   end function cross

end module test_moz_1d_gradient
