!> Numerical checks of HNC, KH and PSE2-4 MOZ closures
module test_moz_closure
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads, omp_get_dynamic, &
!$    & omp_set_dynamic, omp_get_num_threads, omp_in_parallel, omp_get_thread_limit, &
!$    & omp_get_max_active_levels, omp_set_max_active_levels
   use mctc_env, only: wp, moist_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
   use moist_model_moz_closure, only: moist_model_moz_closure_type, &
      & moist_model_moz_closure_hnc_type, moist_model_moz_closure_kh_type, &
      & moist_model_moz_closure_pse_type, new_closure
   use moist_model_moz_closure_kernel, only: closure_hnc_c, closure_kh_c, &
      & closure_pse2_c, closure_pse3_c, closure_pse4_c
   implicit none(type, external)
   private
   public :: collect_moz_closure
   character(len=4), parameter :: names(*) = ["HNC ", "KH  ", "PSE2", "PSE3", "PSE4"]

   !> Only the four energy/solve operations: no response derivative bindings
   type, extends(moist_model_moz_closure_type) :: energy_only_closure_type
   contains
      procedure :: apply => custom_apply
      procedure :: derivative => custom_derivative
      procedure :: free_energy_integrand => custom_integrand
      procedure :: local_integrand => custom_local_integrand
   end type energy_only_closure_type

   !> A custom closure may supply just the checked scalar energy derivatives
   type, extends(energy_only_closure_type) :: checked_derivative_closure_type
   contains
      procedure :: free_energy_derivatives => custom_checked_derivatives
   end type checked_derivative_closure_type

   !> Or it may supply just the pure local energy derivatives
   type, extends(energy_only_closure_type) :: local_derivative_closure_type
   contains
      procedure :: local_energy_derivatives => custom_local_derivatives
   end type local_derivative_closure_type

contains

   !> Register construction, values, derivatives, energies and invalid inputs
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_moz_closure(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("construction", test_construction), &
         & new_unittest("reference_values", test_reference_values), &
         & new_unittest("first_derivatives", test_first_derivatives), &
         & new_unittest("second_derivatives", test_second_derivatives), &
         & new_unittest("branch_boundary", test_branch_boundary), &
         & new_unittest("energy_integrands", test_energy_integrands), &
         & new_unittest("physical_inputs", test_physical_inputs), &
         & new_unittest("invalid_inputs", test_invalid_inputs), &
         & new_unittest("noncontiguous_arrays", test_noncontiguous_arrays), &
         & new_unittest("block_sweeps", test_block_sweeps), &
         & new_unittest("nested_teams", test_nested_teams), &
         & new_unittest("fused_derivative", test_fused_derivative), &
         & new_unittest("array_integrands", test_array_integrands), &
         & new_unittest("energy_derivatives", test_energy_derivatives), &
         & new_unittest("custom_extension", test_custom_extension)]
   end subroutine collect_moz_closure

   !> Outer block workers own parallelism even when nested teams are enabled
   subroutine test_nested_teams(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: saved_levels
      logical :: saved_dynamic

      saved_levels = 1
      saved_dynamic = .false.
!$    saved_levels = omp_get_max_active_levels()
!$    saved_dynamic = omp_get_dynamic()
!$    call omp_set_max_active_levels(2)
!$    call omp_set_dynamic(.false.)
      call check_nested_teams(error)
!$    call omp_set_max_active_levels(saved_levels)
!$    call omp_set_dynamic(saved_dynamic)

   contains

      subroutine check_nested_teams(error)
         type(error_type), allocatable, intent(out) :: error
         type :: error_slot_type
            type(moist_error), allocatable :: error
         end type error_slot_type
         type(error_slot_type) :: errors(2)
         class(moist_model_moz_closure_type), allocatable :: closure
         type(moist_error), allocatable :: err
         ! Eight sites make one energy block reach the closure threading threshold.
         integer, parameter :: npts = 8192, npair = 8
         real(wp), allocatable :: h(:, :), c(:, :), reference(:, :), phi(:, :, :)
         integer :: kind, ib, i, j, teams(2), expected(2)
         logical :: active(2), expect_active

         allocate (h(npts, npair), c(npts, npair), reference(npts, npair), phi(npts, npair, 2))
         do j = 1, npair
            do i = 1, npts
               h(i, j) = 0.025_wp*real(mod(i + j, 41) - 20, wp)
               c(i, j) = 0.01_wp*real(mod(i - j, 17), wp)
            end do
         end do
         expect_active = .false.
!$       expect_active = omp_get_thread_limit() >= 2
         do kind = 1, size(names)
            call new_closure(closure, names(kind), err, nthreads=1)
            if (.not. allocated(err)) call closure%free_energy_integrands(h, c, reference, err)
            call check(error, .not. allocated(err), "Serial nested-team reference: "//trim(names(kind)))
            if (allocated(error)) return
            closure%nthreads = 2
            call check(error, closure%team_size(size(h)), 2, &
               & "Top-level closure parallelism must remain available")
            if (allocated(error)) return
            active = .false.
            expected = 2
            !$omp parallel do schedule(static) num_threads(2)
            do ib = 1, 2
!$             active(ib) = omp_in_parallel()
               if (active(ib)) expected(ib) = 1
               teams(ib) = closure%team_size(size(h))
               call closure%free_energy_integrands(h, c, phi(:, :, ib), errors(ib)%error)
            end do
            !$omp end parallel do
            call check(error, all(active .eqv. expect_active), "Outer energy-block team was not exercised")
            if (allocated(error)) return
            call check(error, all(teams == expected), "Closure requested a nested worker team: "//trim(names(kind)))
            if (allocated(error)) return
            do ib = 1, 2
               call check(error, .not. allocated(errors(ib)%error), "Nested energy sweep: "//trim(names(kind)))
               if (allocated(error)) return
               call check(error, all(phi(:, :, ib) == reference), "Nested energy differs: "//trim(names(kind)))
               if (allocated(error)) return
            end do
         end do
      end subroutine check_nested_teams

   end subroutine test_nested_teams

   !> Verify concrete types, orders, whitespace and factory failure cleanup
   !>
   !> @param[out] error Test failure
   subroutine test_construction(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      character(len=6), parameter :: identifiers(*) = &
         & [" HNC  ", " hnc  ", " KH   ", " kh   ", " PSE2 ", " pse2 ", &
         & " PSE3 ", " pse3 ", " PSE4 ", " pse4 "]
      character(len=7), parameter :: unsupported(*) = &
         & ["unknown", "HNC+MS ", "HNC+MV ", "MS     ", "MV     ", "PSE1   ", "PY     ", "VMB    "]
      integer :: i

      do i = 1, size(identifiers)
         call new_closure(closure, identifiers(i), err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call check(error, allocated(closure))
         if (allocated(error)) return
         select type (closure)
         type is (moist_model_moz_closure_hnc_type)
            call check(error, i <= 2)
         type is (moist_model_moz_closure_kh_type)
            call check(error, i >= 3 .and. i <= 4)
         type is (moist_model_moz_closure_pse_type)
            call check(error, closure%order, (i - 1)/2)
         class default
            call check(error, .false., "Unexpected closure type")
         end select
         if (allocated(error)) return
      end do
      do i = 1, size(unsupported)
         call new_closure(closure, "HNC", err)
         call new_closure(closure, unsupported(i), err)
         call check(error, allocated(err))
         if (allocated(error)) return
         call check(error, .not. allocated(closure))
         if (allocated(error)) return
         call check(error, index(err%message, "Unknown MOZ closure") > 0)
         if (allocated(error)) return
      end do
   end subroutine test_construction

   !> Check exponential and polynomial branches against independent values
   !>
   !> @param[out] error Test failure
   subroutine test_reference_values(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      real(wp) :: u(3, 2), t(3, 2), c(3, 2), expected(3, 2), x(3, 2), series, term, lambda
      real(wp), parameter :: beta = 2.0_wp
      integer :: k, i, j, n, order, coupling

      u = reshape([0.2_wp, 0.4_wp, 0.0_wp, -0.2_wp, 0.3_wp, -0.4_wp], shape(u))
      t = reshape([-0.1_wp, 0.6_wp, 0.0_wp, 0.3_wp, -0.2_wp, 0.4_wp], shape(t))
      do k = 1, size(names)
         call new_closure(closure, names(k), err)
         do coupling = 0, 2
            lambda = 0.5_wp*real(coupling, wp)
            x = t - beta*lambda*u
            call closure%apply(u, t, beta, lambda, c, err)
            call check(error, .not. allocated(err))
            if (allocated(error)) return
            order = k - 1
            do j = 1, size(t, 2)
               do i = 1, size(t, 1)
                  if (k == 1 .or. x(i, j) <= 0.0_wp) then
                     expected(i, j) = exp(x(i, j)) - 1.0_wp - t(i, j)
                  else
                     series = 1.0_wp
                     term = 1.0_wp
                     do n = 1, order
                        term = term*x(i, j)/real(n, wp)
                        series = series + term
                     end do
                     expected(i, j) = series - 1.0_wp - t(i, j)
                  end if
               end do
            end do
            call check(error, maxval(abs(c - expected)) < 2.0e-14_wp, trim(names(k)))
            if (allocated(error)) return
         end do
         call closure%apply(0.0_wp*u, 0.0_wp*t, beta, 1.0_wp, c, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call check(error, maxval(abs(c)) < 1.0e-15_wp)
         if (allocated(error)) return
      end do
   end subroutine test_reference_values

   !> Compare analytic first derivatives with centered finite differences
   !>
   !> @param[out] error Test failure
   subroutine test_first_derivatives(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      real(wp) :: u(3, 2), t(3, 2), dc(3, 2), plus(3, 2), minus(3, 2)
      real(wp), parameter :: step = 1.0e-5_wp
      integer :: k

      u = 0.3_wp
      t = reshape([-0.3_wp, 0.5_wp, 1.0_wp, -0.2_wp, 0.7_wp, 0.9_wp], shape(t))
      do k = 1, size(names)
         call new_closure(closure, names(k), err)
         call closure%derivative(u, t, 1.5_wp, 0.4_wp, dc, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call closure%apply(u, t + step, 1.5_wp, 0.4_wp, plus, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call closure%apply(u, t - step, 1.5_wp, 0.4_wp, minus, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call check(error, maxval(abs(dc - (plus - minus)/(2.0_wp*step))) < 2.0e-10_wp, trim(names(k)))
         if (allocated(error)) return
      end do
   end subroutine test_first_derivatives

   !> Fused c and dc/dt match separate apply and derivative, above the threading threshold
   !>
   !> @param[out] error Test failure
   subroutine test_fused_derivative(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      !> Above the parallel threshold, with both branches and a deep hard core
      integer, parameter :: npts = 40000, npair = 2
      real(wp), allocatable :: u(:, :), t(:, :), c(:, :), dc(:, :), fc(:, :), fdc(:, :)
      integer :: i, j, k

      allocate (u(npts, npair), t(npts, npair), c(npts, npair), dc(npts, npair), &
         & fc(npts, npair), fdc(npts, npair))
      do j = 1, npair
         do i = 1, npts
            u(i, j) = 0.05_wp*real(mod(i + 7*j, 41) - 20, wp)
            t(i, j) = 0.1_wp*real(mod(3*i + j, 23) - 11, wp)
         end do
      end do
      u(1, 1) = 60.0_wp
      do k = 1, size(names)
         call new_closure(closure, names(k), err)
         call closure%apply(u, t, 1.5_wp, 0.7_wp, c, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call closure%derivative(u, t, 1.5_wp, 0.7_wp, dc, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call closure%apply_and_derivative(u, t, 1.5_wp, 0.7_wp, fc, fdc, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call check(error, all(fc == c), "Fused values must equal apply: "//trim(names(k)))
         if (allocated(error)) return
         call check(error, maxval(abs(fdc - dc)/(1.0_wp + abs(dc))) < 1.0e-14_wp, &
            & "Fused derivatives must match derivative: "//trim(names(k)))
         if (allocated(error)) return
      end do
      call closure%apply_and_derivative(u, t, 1.5_wp, 0.7_wp, fc, fdc(:1, :), err)
      call check(error, allocated(err), "Derivative shape mismatch")
      if (allocated(error)) return

      call new_closure(closure, "PSE3", err, nthreads=-1)
      call check(error, allocated(err) .and. .not. allocated(closure), "Negative closure thread count")
      if (allocated(error)) return

      ! An explicit context thread count gives the same values as one thread
      call new_closure(closure, "PSE3", err, nthreads=1)
      call check(error, closure%team_size(size(t)) == 1 .and. closure%team_size(10) == 1)
      if (allocated(error)) return
      call closure%apply_and_derivative(u, t, 1.5_wp, 0.7_wp, c, dc, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call new_closure(closure, "PSE3", err, nthreads=3)
      call check(error, closure%nthreads == 3 .and. closure%team_size(size(t)) == 3 &
         & .and. closure%team_size(10) == 1, "Closure team size")
      if (allocated(error)) return
      call closure%apply_and_derivative(u, t, 1.5_wp, 0.7_wp, fc, fdc, err)
      call check(error, .not. allocated(err) .and. all(fc == c) .and. all(fdc == dc), &
         & "Thread count must not change closure values")
   end subroutine test_fused_derivative

   !> Array and pure integrands match the checked scalar integrand
   !>
   !> @param[out] error Test failure
   subroutine test_array_integrands(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      real(wp) :: h(3, 2), c(3, 2), phi(3, 2), scalar
      integer :: i, j, k

      ! Includes hard-core roundoff one ulp below h = -1 and both h branches
      h = reshape([nearest(-1.0_wp, -1.0_wp), -0.4_wp, 0.0_wp, 0.3_wp, 1.2_wp, 2.5_wp], shape(h))
      c = reshape([-2.2_wp, -0.1_wp, 0.2_wp, -0.5_wp, 0.4_wp, 1.1_wp], shape(c))
      do k = 1, size(names)
         call new_closure(closure, names(k), err)
         call closure%free_energy_integrands(h, c, phi, err)
         call check(error, .not. allocated(err), trim(names(k)))
         if (allocated(error)) return
         do j = 1, size(h, 2)
            do i = 1, size(h, 1)
               call closure%free_energy_integrand(h(i, j), c(i, j), scalar, err)
               call check(error, .not. allocated(err) .and. phi(i, j) == scalar &
                  & .and. closure%local_integrand(h(i, j), c(i, j)) == scalar, trim(names(k)))
               if (allocated(error)) return
            end do
         end do
         h(2, 1) = -1.1_wp
         call closure%free_energy_integrands(h, c, phi, err)
         call check(error, allocated(err), "Unphysical h in an array")
         if (allocated(error)) return
         h(2, 1) = -0.4_wp
         call closure%free_energy_integrands(h, c(:, :1), phi, err)
         call check(error, allocated(err), "Energy array shape mismatch")
         if (allocated(error)) return
      end do
   end subroutine test_array_integrands

   !> Check the integrand and potential derivatives used by response calculations
   !>
   !> - Array fields equal the scalar and the unchecked local derivatives
   !> - dPhi/dh and dPhi/dc match central differences of the integrand away from h = 0
   !> - dc/du at fixed t matches central differences of the closure
   !> - Unphysical h and mismatched shapes are refused
   !>
   !> @param[out] error Test failure
   subroutine test_energy_derivatives(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      real(wp), parameter :: step = 1.0e-6_wp
      real(wp) :: h(3, 2), c(3, 2), dh(3, 2), dc(3, 2), sh, sc, lh, lc
      real(wp) :: u(3, 2), t(3, 2), du(3, 2), plus(3, 2), minus(3, 2)
      integer :: i, j, k

      h = reshape([-0.6_wp, -0.2_wp, 0.3_wp, 0.8_wp, 1.5_wp, 2.2_wp], shape(h))
      c = reshape([-1.4_wp, -0.3_wp, 0.2_wp, -0.5_wp, 0.4_wp, 1.1_wp], shape(c))
      u = 0.3_wp
      t = reshape([-0.3_wp, 0.5_wp, 1.0_wp, -0.2_wp, 0.7_wp, 0.9_wp], shape(t))
      do k = 1, size(names)
         call new_closure(closure, names(k), err)
         call closure%free_energy_derivative_fields(h, c, dh, dc, err)
         call check(error, .not. allocated(err), trim(names(k)))
         if (allocated(error)) return
         do j = 1, size(h, 2)
            do i = 1, size(h, 1)
               call closure%free_energy_derivatives(h(i, j), c(i, j), sh, sc, err)
               call closure%local_energy_derivatives(h(i, j), c(i, j), lh, lc)
               call check(error, .not. allocated(err) .and. dh(i, j) == sh .and. dc(i, j) == sc &
                  & .and. lh == sh .and. lc == sc, trim(names(k))//" array and scalar derivatives")
               if (allocated(error)) return
               call check(error, dh(i, j), (closure%local_integrand(h(i, j) + step, c(i, j)) &
                  & - closure%local_integrand(h(i, j) - step, c(i, j)))/(2.0_wp*step), thr=1.0e-8_wp, &
                  & more=trim(names(k))//" dPhi/dh")
               if (allocated(error)) return
               call check(error, dc(i, j), (closure%local_integrand(h(i, j), c(i, j) + step) &
                  & - closure%local_integrand(h(i, j), c(i, j) - step))/(2.0_wp*step), thr=1.0e-8_wp, &
                  & more=trim(names(k))//" dPhi/dc")
               if (allocated(error)) return
            end do
         end do

         call closure%potential_derivative(u, t, 1.5_wp, 0.4_wp, du, err)
         call check(error, .not. allocated(err), trim(names(k)))
         if (allocated(error)) return
         call closure%apply(u + step, t, 1.5_wp, 0.4_wp, plus, err)
         if (.not. allocated(err)) call closure%apply(u - step, t, 1.5_wp, 0.4_wp, minus, err)
         call check(error, .not. allocated(err), trim(names(k)))
         if (allocated(error)) return
         call check(error, maxval(abs(du - (plus - minus)/(2.0_wp*step))) < 1.0e-8_wp, &
            & trim(names(k))//" dc/du at fixed t")
         if (allocated(error)) return

         h(2, 1) = -1.1_wp
         call closure%free_energy_derivative_fields(h, c, dh, dc, err)
         call check(error, allocated(err), "Unphysical h in derivative arrays")
         if (allocated(error)) return
         h(2, 1) = -0.2_wp
         call closure%free_energy_derivative_fields(h, c(:, :1), dh, dc, err)
         call check(error, allocated(err), "Energy derivative array shape mismatch")
         if (allocated(error)) return
      end do
   end subroutine test_energy_derivatives

   !> Energy-only extensions compile and work; response derivatives stay optional
   subroutine test_custom_extension(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(energy_only_closure_type) :: energy_only
      type(checked_derivative_closure_type) :: checked
      type(local_derivative_closure_type) :: local
      real(wp) :: u(2, 1), t(2, 1), values(2, 1), derivative(2, 1)
      real(wp) :: h(2, 1), c(2, 1), phi(2, 1), dh(2, 1), dc(2, 1), sh, sc, scalar_phi

      u(:, 1) = [0.2_wp, 0.4_wp]
      t(:, 1) = [-0.1_wp, 0.6_wp]
      call energy_only%apply_and_derivative(u, t, 2.0_wp, 0.5_wp, values, derivative, err)
      call check(error, .not. allocated(err), more="energy-only closure cannot evaluate")
      if (allocated(error)) return
      call check(error, maxval(abs(values(:, 1) - [-0.3_wp, 0.2_wp])) < 1.0e-15_wp &
         & .and. all(derivative == 1.0_wp), more="custom closure operations were not used")
      if (allocated(error)) return

      h(:, 1) = [0.3_wp, -0.4_wp]
      c(:, 1) = [-0.2_wp, 0.1_wp]
      call energy_only%free_energy_integrand(h(1, 1), c(1, 1), scalar_phi, err)
      call check(error, .not. allocated(err), more="energy-only scalar energy failed")
      if (allocated(error)) return
      call check(error, scalar_phi, 0.125_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call energy_only%free_energy_integrands(h, c, phi, err)
      call check(error, .not. allocated(err), more="energy-only array energy failed")
      if (allocated(error)) return
      call check(error, maxval(abs(phi - 0.125_wp)) < 1.0e-15_wp)
      if (allocated(error)) return

      call energy_only%free_energy_derivatives(h(1, 1), c(1, 1), sh, sc, err)
      call check(error, allocated(err), more="missing scalar response derivatives were accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "unsupported") > 0)
      if (allocated(error)) return
      call energy_only%free_energy_derivative_fields(h, c, dh, dc, err)
      call check(error, allocated(err), more="missing array response derivatives were accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "unsupported") > 0)
      if (allocated(error)) return
      call energy_only%local_energy_derivatives(h(1, 1), c(1, 1), sh, sc)
      call check(error, ieee_is_nan(sh) .and. ieee_is_nan(sc), more="missing pure derivatives invented values")
      if (allocated(error)) return
      call energy_only%potential_derivative(u, t, 2.0_wp, 0.5_wp, values, err)
      call check(error, allocated(err), more="missing potential derivative was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "does not provide a potential derivative") > 0)
      if (allocated(error)) return

      ! Neither response-capable extension has to implement both derivative hooks.
      call checked%free_energy_derivative_fields(h, c, dh, dc, err)
      call check(error, .not. allocated(err), more="checked-only custom derivatives failed")
      if (allocated(error)) return
      call check(error, maxval(abs(dh(:, 1) - [0.5_wp, -0.5_wp])) < 1.0e-15_wp &
         & .and. maxval(abs(dc(:, 1) - [-0.5_wp, 0.5_wp])) < 1.0e-15_wp)
      if (allocated(error)) return
      call local%free_energy_derivatives(h(1, 1), c(1, 1), sh, sc, err)
      call check(error, .not. allocated(err), more="pure-only checked scalar derivatives failed")
      if (allocated(error)) return
      call check(error, abs(sh - 0.5_wp) < 1.0e-15_wp .and. abs(sc + 0.5_wp) < 1.0e-15_wp)
      if (allocated(error)) return
      call local%free_energy_derivative_fields(h, c, dh, dc, err)
      call check(error, .not. allocated(err), more="pure-only custom derivative fields failed")
      if (allocated(error)) return
      call check(error, maxval(abs(dh(:, 1) - [0.5_wp, -0.5_wp])) < 1.0e-15_wp &
         & .and. maxval(abs(dc(:, 1) - [-0.5_wp, 0.5_wp])) < 1.0e-15_wp)
      if (allocated(error)) return

      call checked%free_energy_derivative_fields(h, c(:1, :), dh, dc, err)
      call check(error, allocated(err), more="custom derivative shape mismatch was accepted")
      if (allocated(error)) return
      h(1, 1) = -1.1_wp
      call checked%free_energy_derivative_fields(h, c, dh, dc, err)
      call check(error, allocated(err), more="checked-only derivatives accepted unphysical h")
      if (allocated(error)) return
      call local%free_energy_derivative_fields(h, c, dh, dc, err)
      call check(error, allocated(err), more="pure-only derivatives accepted unphysical h")
   end subroutine test_custom_extension

   !> Linear custom relation c = t - beta*lambda*u
   subroutine custom_apply(self, u, t, beta, lambda, values, error)
      class(energy_only_closure_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :), t(:, :), beta, lambda
      real(wp), intent(out) :: values(:, :)
      type(moist_error), allocatable, intent(out) :: error
      values = t - beta*lambda*u
   end subroutine custom_apply

   !> Its required solver derivative dc/dt is exactly one
   subroutine custom_derivative(self, u, t, beta, lambda, values, error)
      class(energy_only_closure_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :), t(:, :), beta, lambda
      real(wp), intent(out) :: values(:, :)
      type(moist_error), allocatable, intent(out) :: error
      values = 1.0_wp
   end subroutine custom_derivative

   !> Checked scalar interface to the custom quadratic energy
   subroutine custom_integrand(self, h, c, phi, error)
      class(energy_only_closure_type), intent(in) :: self
      real(wp), intent(in) :: h, c
      real(wp), intent(out) :: phi
      type(moist_error), allocatable, intent(out) :: error
      phi = self%local_integrand(h, c)
   end subroutine custom_integrand

   !> Quadratic energy distinct from the built-in closure integrands
   pure function custom_local_integrand(self, h, c) result(phi)
      class(energy_only_closure_type), intent(in) :: self
      real(wp), intent(in) :: h, c
      real(wp) :: phi
      phi = 0.5_wp*(h - c)**2
   end function custom_local_integrand

   !> Checked-only derivative provider; the array wrapper validates its inputs
   subroutine custom_checked_derivatives(self, h, c, dh, dc, error)
      class(checked_derivative_closure_type), intent(in) :: self
      real(wp), intent(in) :: h, c
      real(wp), intent(out) :: dh, dc
      type(moist_error), allocatable, intent(out) :: error
      dh = h - c
      dc = c - h
   end subroutine custom_checked_derivatives

   !> Pure-only derivative provider, with the checked interfaces inherited
   pure subroutine custom_local_derivatives(self, h, c, dh, dc)
      class(local_derivative_closure_type), intent(in) :: self
      real(wp), intent(in) :: h, c
      real(wp), intent(out) :: dh, dc
      dh = h - c
      dc = c - h
   end subroutine custom_local_derivatives

   !> Compare d2c/dt2 with finite differences of dc/dt away from the branch point
   !>
   !> @param[out] error Test failure
   subroutine test_second_derivatives(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      real(wp) :: u(3, 2), t(3, 2), d2c(3, 2), plus(3, 2), minus(3, 2), jump(2, 1)
      real(wp), parameter :: step = 1.0e-5_wp
      integer :: k

      ! t_star = t - 0.18 stays at least 0.3 away from the branch point
      u = 0.3_wp
      t = reshape([-0.3_wp, 0.5_wp, 1.0_wp, -0.2_wp, 0.7_wp, 0.9_wp], shape(t))
      do k = 1, size(names)
         call new_closure(closure, names(k), err)
         call closure%derivative2(u, t, 1.5_wp, 0.4_wp, d2c, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call closure%derivative(u, t + step, 1.5_wp, 0.4_wp, plus, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call closure%derivative(u, t - step, 1.5_wp, 0.4_wp, minus, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call check(error, maxval(abs(d2c - (plus - minus)/(2.0_wp*step))) < 2.0e-9_wp, trim(names(k)))
         if (allocated(error)) return
      end do

      ! KH: d2c/dt2 tends to one below the branch point and is zero above it
      call new_closure(closure, "KH", err)
      call closure%derivative2(u(:2, :1), reshape([-1.0e-12_wp, 1.0e-12_wp], [2, 1]), &
         & 1.0_wp, 0.0_wp, jump, err)
      call check(error, .not. allocated(err) .and. abs(jump(1, 1) - 1.0_wp) < 1.0e-11_wp &
         & .and. jump(2, 1) == 0.0_wp, "KH second derivative jump")
      if (allocated(error)) return
      call closure%derivative2(u, t, 1.5_wp, 1.5_wp, d2c, err)
      call check(error, allocated(err), "Invalid lambda")
   end subroutine test_second_derivatives

   !> Check values and first derivatives across x=0 for KH and PSE
   !>
   !> @param[out] error Test failure
   subroutine test_branch_boundary(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      real(wp) :: u(3, 1), t(3, 1), c(3, 1), dc(3, 1)
      integer :: k

      u = 0.5_wp
      t(:, 1) = [0.5_wp - 1.0e-8_wp, 0.5_wp, 0.5_wp + 1.0e-8_wp]
      do k = 2, 5
         call new_closure(closure, names(k), err)
         call closure%apply(u, t, 1.0_wp, 1.0_wp, c, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call closure%derivative(u, t, 1.0_wp, 1.0_wp, dc, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call check(error, maxval(abs(c + 0.5_wp)) < 1.0e-14_wp)
         if (allocated(error)) return
         call check(error, maxval(abs(dc)) < 1.1e-8_wp)
         if (allocated(error)) return
      end do
   end subroutine test_branch_boundary

   !> Check HNC, KH and rational PSE energy references
   !>
   !> @param[out] error Test failure
   subroutine test_energy_integrands(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      real(wp), parameter :: h(*) = [1.5_wp, 5.0_wp/3.0_wp, 41.0_wp/24.0_wp]
      real(wp), parameter :: expected(*) = [89.0_wp/60.0_wp, 683.0_wp/360.0_wp, 11561.0_wp/5760.0_wp]
      real(wp) :: phi
      integer :: k

      do k = 1, 5
         call new_closure(closure, names(k), err)
         call closure%free_energy_integrand(-0.4_wp, 0.2_wp, phi, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call check(error, abs(phi + 0.08_wp) < 1.0e-14_wp)
         if (allocated(error)) return
         if (k >= 3) then
            call closure%free_energy_integrand(h(k - 2), -0.3_wp, phi, err)
            call check(error, .not. allocated(err))
            if (allocated(error)) return
            call check(error, abs(phi - expected(k - 2)) < 1.0e-13_wp)
         else
            call closure%free_energy_integrand(0.4_wp, 0.2_wp, phi, err)
            call check(error, .not. allocated(err))
            if (allocated(error)) return
            if (k == 1) call check(error, abs(phi + 0.16_wp) < 1.0e-14_wp)
            if (k == 2) call check(error, abs(phi + 0.24_wp) < 1.0e-14_wp)
         end if
         if (allocated(error)) return
      end do

   end subroutine test_energy_integrands

   !> Check caller reconstruction of physical trial inputs and Ng renormalization
   !>
   !> @param[out] error Test failure
   subroutine test_physical_inputs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      real(wp) :: u(2, 1), t(2, 1), c(2, 1), shift(2, 1), lambda
      integer :: k, coupling

      ! The caller's short-range potential and iterate are both zero.
      u(:, 1) = [1.0_wp, -0.6_wp]
      do k = 1, size(names)
         call new_closure(closure, names(k), err)
         do coupling = 0, 2
            lambda = 0.5_wp*real(coupling, wp)
            shift = lambda*u
            t = shift
            call closure%apply(u, t, 1.0_wp, lambda, c, err)
            call check(error, .not. allocated(err))
            if (allocated(error)) return
            call check(error, maxval(abs(c + shift)) < 2.0e-14_wp, trim(names(k)))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_physical_inputs

   !> Reject invalid shapes, temperature, scaling and non-finite correlations
   !>
   !> @param[out] error Test failure
   subroutine test_invalid_inputs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      real(wp) :: u(2, 1), t(2, 1), c(2, 1), wrong(1, 1), phi, nan, inf
      integer :: k

      nan = ieee_value(0.0_wp, ieee_quiet_nan)
      inf = ieee_value(0.0_wp, ieee_positive_inf)
      do k = 1, size(names)
         u = 0.0_wp
         t = 0.0_wp
         call new_closure(closure, names(k), err)
         call closure%apply(u, t, 1.0_wp, 1.0_wp, wrong, err)
         call check(error, allocated(err), "Output shape mismatch")
         if (allocated(error)) return
         call closure%derivative(u, t(:1, :), 1.0_wp, 1.0_wp, c, err)
         call check(error, allocated(err), "Input shape mismatch")
         if (allocated(error)) return
         call closure%apply(u, t, 0.0_wp, 1.0_wp, c, err)
         call check(error, allocated(err), "Nonpositive beta")
         if (allocated(error)) return
         call closure%derivative(u, t, 1.0_wp, 1.1_wp, c, err)
         call check(error, allocated(err), "Invalid lambda")
         if (allocated(error)) return
         call closure%apply(u, t, nan, 1.0_wp, c, err)
         call check(error, allocated(err), "NaN beta")
         if (allocated(error)) return
         u(1, 1) = inf
         call closure%apply(u, t, 1.0_wp, 1.0_wp, c, err)
         call check(error, allocated(err), "Infinite potential")
         if (allocated(error)) return
         u = 0.0_wp
         t(1, 1) = nan
         call closure%derivative(u, t, 1.0_wp, 1.0_wp, c, err)
         call check(error, allocated(err), "NaN correlation")
         if (allocated(error)) return
         call closure%free_energy_integrand(-1.1_wp, 0.0_wp, phi, err)
         call check(error, allocated(err), "Unphysical h")
         if (allocated(error)) return
         ! h = c + t of a deep hard core can round to one ulp below -1
         call closure%free_energy_integrand(nearest(-1.0_wp, -1.0_wp), 0.0_wp, phi, err)
         call check(error, .not. allocated(err), "Hard-core roundoff below h = -1")
         if (allocated(error)) return
         call closure%free_energy_integrand(0.0_wp, nan, phi, err)
         call check(error, allocated(err), "NaN c")
         if (allocated(error)) return
      end do
      call new_closure(closure, "PSE2", err)
      select type (closure)
      type is (moist_model_moz_closure_pse_type)
         closure%order = 5
      end select
      t = 0.0_wp
      call closure%apply(u, t, 1.0_wp, 1.0_wp, c, err)
      call check(error, allocated(err), "Invalid PSE order")
      if (allocated(error)) return
      call new_closure(closure, "HNC", err)
      t = 1000.0_wp
      call closure%apply(u, t, 1.0_wp, 1.0_wp, c, err)
      call check(error, allocated(err), "Exponential overflow")
   end subroutine test_invalid_inputs

   !> Check strided sections and ensure untouched output rows remain intact
   !>
   !> @param[out] error Test failure
   subroutine test_noncontiguous_arrays(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure
      type(moist_error), allocatable :: err
      real(wp) :: u(4, 2), t(4, 2), c(4, 2), expected(4, 2)
      integer :: k

      u = reshape([0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp, -0.1_wp, -0.2_wp, -0.3_wp, -0.4_wp], shape(u))
      t = reshape([-0.2_wp, 0.0_wp, 0.4_wp, 0.7_wp, -0.3_wp, -0.1_wp, 0.3_wp, 0.6_wp], shape(t))
      do k = 1, size(names)
         call new_closure(closure, names(k), err)
         call closure%apply(u, t, 1.0_wp, 1.0_wp, expected, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         c = 42.0_wp
         call closure%apply(u(1::2, :), t(1::2, :), 1.0_wp, 1.0_wp, c(1::2, :), err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call check(error, maxval(abs(c(1::2, :) - expected(1::2, :))) < 1.0e-14_wp)
         if (allocated(error)) return
         call check(error, all(c(2::2, :) == 42.0_wp))
         if (allocated(error)) return
         call closure%derivative(u, t, 1.0_wp, 1.0_wp, expected, err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         c = 42.0_wp
         call closure%derivative(u(1::2, :), t(1::2, :), 1.0_wp, 1.0_wp, c(1::2, :), err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call check(error, maxval(abs(c(1::2, :) - expected(1::2, :))) < 1.0e-14_wp)
         if (allocated(error)) return
         call check(error, all(c(2::2, :) == 42.0_wp))
         if (allocated(error)) return
      end do
   end subroutine test_noncontiguous_arrays

   !> Compare block dispatch and strided batches with kernels using real OpenMP teams
   !>
   !> @param[out] error Test failure
   subroutine test_block_sweeps(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> OpenMP settings restored even when a comparison fails
      integer :: saved_threads
      logical :: saved_dynamic, has_openmp, in_parallel

      saved_threads = 1
      saved_dynamic = .false.
      has_openmp = .false.
      in_parallel = .false.
!$    saved_threads = omp_get_max_threads()
!$    saved_dynamic = omp_get_dynamic()
!$    has_openmp = .true.
!$    in_parallel = omp_in_parallel()
!$    call omp_set_dynamic(.false.)

      call check_block_sweeps(error, has_openmp, in_parallel)

!$    call omp_set_num_threads(saved_threads)
!$    call omp_set_dynamic(saved_dynamic)

   contains

      !> Failure returns to the caller, which restores OpenMP settings before returning
      !>
      !> @param[out] error Test failure
      !> @param[in] has_openmp Whether OpenMP is enabled
      !> @param[in] in_parallel Whether the caller is inside a parallel team
      subroutine check_block_sweeps(error, has_openmp, in_parallel)
         !> Test failure
         type(error_type), allocatable, intent(out) :: error
         !> Whether OpenMP is enabled and the caller is inside a parallel team
         logical, intent(in) :: has_openmp, in_parallel
         !> Above the threading threshold, with a three-point tail after sixteen blocks
         integer, parameter :: npts = 65539, npair = 3
         integer, parameter :: widths(2) = [1, npair], thread_counts(2) = [1, 4]
         real(wp), parameter :: beta = 2.0_wp, lambda = 0.5_wp, sentinel = 42.0_wp
         class(moist_model_moz_closure_type), allocatable :: closure
         type(moist_error), allocatable :: err
         real(wp), allocatable :: u(:, :), t(:, :), values(:, :), expected(:, :)
         real(wp), allocatable :: strided_u(:, :), strided_t(:, :), strided_values(:, :)
         real(wp) :: saved_value, nan, inf
         integer :: i, j, kind, run, batch, width, observed_threads, expected_threads

         call check(error, .not. in_parallel, &
            & "Block-sweep test must run as a selected test outside the suite OpenMP team")
         if (allocated(error)) return

         nan = ieee_value(0.0_wp, ieee_quiet_nan)
         inf = ieee_value(0.0_wp, ieee_positive_inf)
         allocate (u(npts, npair), t(npts, npair), values(npts, npair), expected(npts, npair))
         allocate (strided_u(2*npts, 2*npair - 1), strided_t(2*npts, 2*npair - 1), &
            & strided_values(2*npts, 2*npair - 1))
         do j = 1, npair
            do i = 1, npts
               u(i, j) = 0.03125_wp*real(mod(i + 3*j, 17) - 8, wp)
               ! Exact binary increments include both closure branches and their join
               t(i, j) = beta*lambda*u(i, j) + 0.0625_wp*real(mod(i + 5*j, 15) - 7, wp)
            end do
         end do
         strided_u = 0.0_wp
         strided_t = 0.0_wp
         strided_u(1::2, 1::2) = u
         strided_t(1::2, 1::2) = t

         do kind = 1, size(names)
            call new_closure(closure, names(kind), err)
            call check(error, .not. allocated(err))
            if (allocated(error)) return
            do j = 1, npair
               select case (kind)
               case (1)
                  call closure_hnc_c(u(:, j), t(:, j), beta, lambda, expected(:, j))
               case (2)
                  call closure_kh_c(u(:, j), t(:, j), beta, lambda, expected(:, j))
               case (3)
                  call closure_pse2_c(u(:, j), t(:, j), beta, lambda, expected(:, j))
               case (4)
                  call closure_pse3_c(u(:, j), t(:, j), beta, lambda, expected(:, j))
               case (5)
                  call closure_pse4_c(u(:, j), t(:, j), beta, lambda, expected(:, j))
               case default
                  call check(error, .false., "Unknown test closure")
                  return
               end select
            end do

            do run = 1, size(thread_counts)
!$             call omp_set_num_threads(thread_counts(run))
               expected_threads = 1
!$             expected_threads = min(thread_counts(run), omp_get_thread_limit())
               observed_threads = 1
               !$omp parallel reduction(max:observed_threads)
!$             observed_threads = omp_get_num_threads()
               !$omp end parallel
               if (has_openmp) then
                  call check(error, observed_threads == expected_threads, &
                     & "Block-sweep test must create the requested OpenMP team")
                  if (allocated(error)) return
               end if
               do batch = 1, size(widths)
                  width = widths(batch)
                  call closure%apply(u(:, 1:width), t(:, 1:width), beta, lambda, values(:, 1:width), err)
                  call check(error, .not. allocated(err))
                  if (allocated(error)) return
                  call check(error, all(values(:, 1:width) == expected(:, 1:width)), &
                     & "Contiguous block dispatch must match kernels exactly")
                  if (allocated(error)) return

                  strided_values = sentinel
                  call closure%apply(strided_u(1::2, 1:2*width - 1:2), &
                     & strided_t(1::2, 1:2*width - 1:2), beta, lambda, &
                     & strided_values(1::2, 1:2*width - 1:2), err)
                  call check(error, .not. allocated(err))
                  if (allocated(error)) return
                  call check(error, all(strided_values(1::2, 1:2*width - 1:2) == expected(:, 1:width)), &
                     & "Strided block dispatch must match kernels exactly")
                  if (allocated(error)) return
                  call check(error, all(strided_values(2::2, :) == sentinel) .and. &
                     & all(strided_values(:, 2::2) == sentinel), &
                     & "Strided block sweeps must preserve excluded rows and columns")
                  if (allocated(error)) return
                  if (width < npair) then
                     call check(error, all(strided_values(:, 2*width + 1:) == sentinel), &
                        & "Single-column block sweeps must preserve the other columns")
                     if (allocated(error)) return
                  end if
               end do

               ! Invalid values in later blocks and the tail must reach the caller.
               strided_u(2*npts - 1, 2*npair - 1) = inf
               call closure%apply(strided_u(1::2, 1::2), strided_t(1::2, 1::2), beta, lambda, &
                  & strided_values(1::2, 1::2), err)
               call check(error, allocated(err), "Infinite potential in a strided tail block")
               if (allocated(error)) return
               call check(error, index(err%message, "potential and correlation must be finite") > 0)
               if (allocated(error)) return
               strided_u(2*npts - 1, 2*npair - 1) = u(npts, npair)

               saved_value = t(4097, 1)
               t(4097, 1) = nan
               call closure%apply(u, t, beta, lambda, values, err)
               call check(error, allocated(err), "NaN correlation after the first block")
               if (allocated(error)) return
               call check(error, index(err%message, "potential and correlation must be finite") > 0)
               if (allocated(error)) return
               t(4097, 1) = saved_value

               saved_value = u(npts, npair)
               u(npts, npair) = -huge(1.0_wp)
               call closure%apply(u, t, beta, 1.0_wp, values, err)
               call check(error, allocated(err), "Non-finite output in the tail block")
               if (allocated(error)) return
               call check(error, index(err%message, "kernel returned a non-finite result") > 0)
               if (allocated(error)) return
               u(npts, npair) = saved_value

               call closure%apply(u, t, beta, lambda, values, err)
               call check(error, .not. allocated(err), "Valid call after a failed block")
               if (allocated(error)) return
               call check(error, all(values == expected), "Recovery must preserve the kernel values")
               if (allocated(error)) return
            end do
         end do
      end subroutine check_block_sweeps

   end subroutine test_block_sweeps

end module test_moz_closure
