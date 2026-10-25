!> Checks of the 1D MOZ spectral filter and block k-space preconditioner
module test_moz_1d_preconditioner
   use mctc_env, only: wp, moist_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair, &
      & transform_quadrature
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo
   use moist_context, only: moist_context_type, new_context
   use moist_model_moz_solver_core_preconditioner, only: moist_model_moz_solver_preconditioner_type
   use moist_model_moz_solver, only: moist_model_moz_solver_type, moist_model_moz_solver_evaluator_type, &
      & moist_model_moz_solver_parameters_type, new_solver
   use moist_model_moz_solver_jfnk, only: moist_model_moz_solver_jfnk_type, new_jfnk_solver
   use moist_model_moz_1d_preconditioner, only: moist_model_moz_1d_preconditioner_spectral_type, &
      & moist_model_moz_1d_preconditioner_kspace_type, select_1d_preconditioner, &
      & preconditioners_1d_uv, preconditioners_1d_vv
   use moist_model_moz_1d_solver_filters, only: moist_model_moz_1d_solver_filters_type
   implicit none(type, external)
   private
   public :: collect_moz_1d_preconditioner

   !> Residual R(gamma) = gamma with its root at zero
   type, extends(moist_model_moz_solver_evaluator_type) :: identity_residual_type
   contains
      procedure :: residual => identity_residual
   end type identity_residual_type

   !> Radial nodes and spacing of the test grids, bohr; 48 has an odd factor,
   !> so the characteristic-length index meets exact rounding ties
   integer, parameter :: npts = 48
   real(wp), parameter :: dr = 0.25_wp
   !> Filter parameters
   real(wp), parameter :: k0 = 1.3_wp, w0 = 0.2_wp

contains

   !> Register filter, block inverse, activation, threading and selection checks
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_moz_1d_preconditioner(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("spectral_identity", test_spectral_identity), &
         & new_unittest("spectral_filter", test_spectral_filter), &
         & new_unittest("spectral_thread_independence", test_spectral_threads), &
         & new_unittest("spectral_setup_errors", test_spectral_setup_errors), &
         & new_unittest("kspace_block_inverse", test_kspace_block_inverse), &
         & new_unittest("kspace_activation", test_kspace_activation), &
         & new_unittest("kspace_thread_independence", test_kspace_threads), &
         & new_unittest("kspace_update_errors", test_kspace_update_errors), &
         & new_unittest("kspace_reallocation", test_kspace_reallocation), &
         & new_unittest("kspace_singular_fallback", test_kspace_singular), &
         & new_unittest("kspace_release_guard", test_kspace_release_guard), &
         & new_unittest("kspace_block_query", test_kspace_block_query), &
         & new_unittest("selection", test_selection), &
         & new_unittest("filters_damping", test_filters_damping), &
         & new_unittest("filters_attach", test_filters_attach)]
   end subroutine collect_moz_1d_preconditioner

   !> Build the uniform test grid pair
   !>
   !> @param[out] rgrid r-space grid
   !> @param[out] kgrid k-space grid
   !> @param[out] error Test failure
   subroutine make_pair(rgrid, kgrid, error)
      !> r-space and k-space grids
      type(moist_math_grid_radial_type), intent(out) :: rgrid, kgrid
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      call new_uniform_radial_pair(rgrid, kgrid, npts, dr, err)
      call check(error, .not. allocated(err), "Uniform radial pair")
   end subroutine make_pair

   !> Smooth flat test vector, pair-major
   !>
   !> @param[in] npair Number of site pairs
   !> @return v Flat vector (npts*npair)
   function test_vector(npair) result(v)
      !> Number of site pairs
      integer, intent(in) :: npair
      !> Flat vector
      real(wp), allocatable :: v(:)
      integer :: i, ip
      allocate (v(npts*npair))
      do ip = 1, npair
         do i = 1, npts
            v(i + (ip - 1)*npts) = sin(0.3_wp*i + 0.7_wp*ip)*exp(-0.02_wp*i)
         end do
      end do
   end function test_vector

   !> A unit weight at k = 0 makes the filter the identity up to transform roundoff
   !>
   !> @param[out] error Test failure
   subroutine test_spectral_identity(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_1d_preconditioner_spectral_type) :: spectral
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_error), allocatable :: err
      real(wp), allocatable :: v(:), original(:)

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call check(error, .not. spectral%ready(0.0_wp) .and. spectral%describe(0.0_wp) == "spectral, not built")
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, 1.0_wp, err)
      call check(error, .not. allocated(err) .and. spectral%ready(huge(1.0_wp)) .and. &
         & spectral%describe(1.0_wp) == "spectral")
      if (allocated(error)) return
      original = test_vector(3)
      v = original
      call spectral%apply(v)
      call check(error, maxval(abs(v - original)) < 1.0e-12_wp, "w0 = 1 must leave the vector unchanged")
      if (allocated(error)) return
      ! A length that is not a multiple of npts is left untouched
      v = original(:npts + 5)
      call spectral%apply(v)
      call check(error, all(v == original(:npts + 5)))
      if (allocated(error)) return
      call spectral%release()
      call spectral%release()
      call check(error, .not. spectral%ready(0.0_wp), "Release is repeatable")
   end subroutine test_spectral_identity

   !> In k-space the filtered vector equals w(k) times the transform of the input
   !>
   !> @param[out] error Test failure
   subroutine test_spectral_filter(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_1d_preconditioner_spectral_type) :: spectral
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(moist_error), allocatable :: err
      real(wp), allocatable :: v(:), filtered(:), vk(:, :), fk(:, :), weight(:)
      integer, parameter :: npair = 3
      integer :: ip

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, w0, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call new_radial_trafo(trafo, rgrid, kgrid, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      v = test_vector(npair)
      filtered = v
      call spectral%apply(filtered)
      allocate (vk(npts, npair), fk(npts, npair))
      call trafo%fbt_r2k_all(reshape(v, [npts, npair]), vk, err)
      call trafo%fbt_r2k_all(reshape(filtered, [npts, npair]), fk, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      weight = w0 + (1.0_wp - w0)*kgrid%r**2/(kgrid%r**2 + k0**2)
      do ip = 1, npair
         call check(error, maxval(abs(fk(:, ip) - weight*vk(:, ip))) < 1.0e-11_wp*(1.0_wp + maxval(abs(vk))), &
            & "Filter must scale each k component by w(k)")
         if (allocated(error)) return
      end do
      call check(error, maxval(abs(filtered - v)) > 1.0e-3_wp, "The filter must change a generic vector")
   end subroutine test_spectral_filter

   !> Thread counts, repeated applications and changing field sizes agree with the batched transform
   !>
   !> @param[out] error Test failure
   subroutine test_spectral_threads(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_1d_preconditioner_spectral_type) :: spectral
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(moist_error), allocatable :: err
      real(wp), allocatable :: v(:), reference(:), vk(:, :), vr(:, :), weight(:)
      real(wp), allocatable :: strided_in(:, :), strided_k(:, :), strided_r(:, :)
      integer, parameter :: widths(3) = [5, 2, 7], threads(4) = [1, 4, 2, 0]
      integer :: iw, it, rep, ip, npair, i

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, w0, err)
      if (.not. allocated(err)) call new_radial_trafo(trafo, rgrid, kgrid, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      weight = w0 + (1.0_wp - w0)*kgrid%r**2/(kgrid%r**2 + k0**2)
      do iw = 1, size(widths)
         npair = widths(iw)
         allocate (vk(npts, npair), vr(npts, npair))
         do it = 1, size(threads)
            spectral%nthreads = threads(it)
            v = test_vector(npair)
            reference = v
            do rep = 1, 3
               call trafo%fbt_r2k_all(reshape(reference, [npts, npair]), vk, err)
               call check(error, .not. allocated(err))
               if (allocated(error)) return
               do ip = 1, npair
                  vk(:, ip) = vk(:, ip)*weight
               end do
               call trafo%fbt_k2r_all(vk, vr, err)
               call check(error, .not. allocated(err))
               if (allocated(error)) return
               reference = reshape(vr, [npts*npair])
               call spectral%apply(v)
               call check(error, maxval(abs(v - reference)) < 1.0e-12_wp, &
                  & "Parallel filtering or reused scratch changed the filtered field")
               if (allocated(error)) return
            end do
         end do
         deallocate (vk, vr)
      end do

      ! The public transforms also accept strided inputs and outputs.
      npair = 5
      allocate (strided_in(2*npts, 2*npair), strided_k(2*npts, 2*npair), strided_r(2*npts, 2*npair))
      strided_in = -999.0_wp
      strided_k = -999.0_wp
      strided_r = -999.0_wp
      strided_in(1::2, 1::2) = reshape(test_vector(npair), [npts, npair])
      spectral%nthreads = 4
      call spectral%forward(strided_in(1::2, 1::2), strided_k(1::2, 1::2), err)
      if (.not. allocated(err)) call spectral%backward(strided_k(1::2, 1::2), strided_r(1::2, 1::2), err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      v = test_vector(npair)
      call spectral%apply(v)
      call check(error, maxval(abs(reshape(strided_r(1::2, 1::2), [npts*npair]) - v)) < 1.0e-12_wp .and. &
         & all(strided_k(2::2, :) == -999.0_wp) .and. all(strided_k(:, 2::2) == -999.0_wp) .and. &
         & all(strided_r(2::2, :) == -999.0_wp) .and. all(strided_r(:, 2::2) == -999.0_wp), &
         & "Strided filtering changed the result or wrote outside the output view")
      if (allocated(error)) return
      call spectral%forward(strided_in(1::2, 1::2), strided_k(1::2, 1:2*npair - 2:2), err)
      call check(error, allocated(err), "Forward filtering accepted an output of the wrong width")
      if (allocated(error)) return
      call spectral%backward(strided_k(1::2, 1::2), strided_r(1:npts - 1, 1::2), err)
      call check(error, allocated(err), "Backward filtering accepted an output of the wrong height")
      if (allocated(error)) return

      ! Replacing the grid must discard scratch from the previous setup.
      call new_uniform_radial_pair(rgrid, kgrid, npts + 4, dr, err)
      if (.not. allocated(err)) call spectral%setup(rgrid, kgrid, k0, 1.0_wp, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      reference = [(sin(0.1_wp*i), i = 1, 3*(npts + 4))]
      v = reference
      call spectral%apply(v)
      call check(error, maxval(abs(v - reference)) < 1.0e-12_wp, "Replacing the grid left stale filter scratch")
      if (allocated(error)) return
      v = reference(:0)
      call spectral%apply(v)
      call check(error, size(v) == 0, "An empty vector failed")
   end subroutine test_spectral_threads

   !> Invalid filter parameters or grids are reported and leave the filter unbuilt
   !>
   !> @param[out] error Test failure
   subroutine test_spectral_setup_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_1d_preconditioner_spectral_type) :: spectral
      type(moist_math_grid_radial_type) :: rgrid, kgrid, other_r, other_k
      type(moist_error), allocatable :: err
      real(wp), allocatable :: v(:)

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, 0.0_wp, w0, err)
      call check(error, allocated(err) .and. .not. spectral%ready(0.0_wp), "k0 must be positive")
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, 1.5_wp, err)
      call check(error, allocated(err) .and. .not. spectral%ready(0.0_wp), "w0 must lie in [0, 1]")
      if (allocated(error)) return
      call new_uniform_radial_pair(other_r, other_k, npts + 4, dr, err)
      call spectral%setup(rgrid, other_k, k0, w0, err)
      call check(error, allocated(err) .and. .not. spectral%ready(0.0_wp), "Grids must be paired")
      if (allocated(error)) return
      ! Quadrature pairs (e.g. descending Chebyshev nodes) are not supported
      other_r = rgrid
      other_k = kgrid
      other_r%transform = transform_quadrature
      other_k%transform = transform_quadrature
      call spectral%setup(other_r, other_k, k0, w0, err)
      call check(error, allocated(err) .and. .not. spectral%ready(0.0_wp), "Only the uniform DST-IV pair")
      if (allocated(error)) return
      v = test_vector(2)
      call spectral%apply(v)
      call check(error, all(v == test_vector(2)), "An unbuilt filter is a no-op")
   end subroutine test_spectral_setup_errors

   !> Build packed model inputs for ns sites; dc/dgamma varies with r and pair
   !> and leaves [-1, 0], so the running average and its clamp are exercised
   !>
   !> @param[in]  ns    Number of sites
   !> @param[in]  r     r nodes (n); the fields get n rows
   !> @param[in]  k     k nodes (n)
   !> @param[out] wk    Intramolecular correlation (npts, npair)
   !> @param[out] cr    Direct correlation (npts, npair)
   !> @param[out] uk    Long-range potential (npts, npair)
   !> @param[out] dcdr  Closure derivative (npts, npair)
   !> @param[out] rho   Site densities (ns)
   subroutine model_inputs(ns, r, k, wk, cr, uk, dcdr, rho)
      !> Number of sites
      integer, intent(in) :: ns
      !> Radial and reciprocal nodes
      real(wp), intent(in) :: r(:), k(:)
      !> Packed pair fields (npts, npair)
      real(wp), allocatable, intent(out) :: wk(:, :), cr(:, :), uk(:, :), dcdr(:, :)
      !> Site densities
      real(wp), allocatable, intent(out) :: rho(:)
      integer :: is, js, ip

      allocate (wk(size(r), ns*(ns + 1)/2), cr(size(r), ns*(ns + 1)/2), uk(size(r), ns*(ns + 1)/2), &
         & dcdr(size(r), ns*(ns + 1)/2), rho(ns))
      do js = 1, ns
         rho(js) = 0.05_wp + 0.01_wp*js
         do is = 1, js
            ip = js*(js - 1)/2 + is
            if (is == js) then
               wk(:, ip) = 1.0_wp
            else
               wk(:, ip) = 0.2_wp*exp(-0.1_wp*k**2*real(is + js, wp)/real(ns, wp))
            end if
            cr(:, ip) = -0.3_wp*(1.0_wp + 0.1_wp*mod(ip, 7))*exp(-r**2)
            uk(:, ip) = 0.01_wp*exp(-k**2)/real(ip, wp)
            dcdr(:, ip) = -1.6_wp*exp(-r**2/(1.0_wp + 0.05_wp*ip)) + 0.4_wp*sin(0.7_wp*r)
         end do
      end do
   end subroutine model_inputs

   !> The block solve equals the symmetrized dense inverse of M(k) built here
   !>
   !> @param[out] error Test failure
   subroutine test_kspace_block_inverse(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ns = 14, npair = ns*(ns + 1)/2
      real(wp), parameter :: beta = 1.7_wp
      type(moist_model_moz_1d_preconditioner_spectral_type), target :: spectral
      type(moist_model_moz_1d_preconditioner_kspace_type) :: kspace
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(moist_error), allocatable :: err
      real(wp), allocatable :: wk(:, :), cr(:, :), uk(:, :), dcdr(:, :), rho(:)
      real(wp), allocatable :: v(:), result(:), vk(:, :), ck(:, :), zk(:, :), expected(:, :)
      real(wp), allocatable :: alpha(:, :), num(:), den(:)
      real(wp) :: m(ns, ns), w_pt(ns, ns), crho(ns, ns), rhs(ns, ns), weight
      integer :: ik, is, js, ip, j

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, w0, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call new_radial_trafo(trafo, rgrid, kgrid, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call model_inputs(ns, rgrid%r, kgrid%r, wk, cr, uk, dcdr, rho)
      call kspace%bind(spectral)
      call kspace%update(wk, cr, rho, beta, uk, dcdr, err)
      call check(error, .not. allocated(err) .and. kspace%valid .and. kspace%ready(1.0e-4_wp))
      if (allocated(error)) return

      v = test_vector(npair)
      result = v
      call kspace%apply(result)

      ! Independent reference: running average of dc/dgamma up to node
      ! nint(min(r_max, pi/k)/dr + 0.5), filtered transform and
      ! dense solve
      allocate (alpha(npts, npair), num(npts), den(npts))
      do ip = 1, npair
         num(1) = dcdr(1, ip)*rgrid%r(1)**2*dr
         den(1) = rgrid%r(1)**2*dr
         do j = 2, npts
            num(j) = num(j - 1) + dcdr(j, ip)*rgrid%r(j)**2*dr
            den(j) = den(j - 1) + rgrid%r(j)**2*dr
         end do
         do ik = 1, npts
            j = min(npts, max(1, nint(min(rgrid%r(npts), acos(-1.0_wp)/kgrid%r(ik))/dr + 0.5_wp)))
            alpha(ik, ip) = max(-1.0_wp, min(0.0_wp, num(j)/den(j)))
         end do
      end do
      call check(error, any(alpha == -1.0_wp) .and. any(alpha > -1.0_wp .and. alpha < 0.0_wp), &
         & "The inputs must exercise the clamp and the running average")
      if (allocated(error)) return
      allocate (vk(npts, npair), ck(npts, npair), zk(npts, npair), expected(npts, npair))
      call trafo%fbt_r2k_all(reshape(v, [npts, npair]), vk, err)
      call trafo%fbt_r2k_all(cr, ck, err)
      ck = ck - beta*uk
      do ik = 1, npts
         weight = w0 + (1.0_wp - w0)*kgrid%r(ik)**2/(kgrid%r(ik)**2 + k0**2)
         do js = 1, ns
            do is = 1, ns
               ip = max(is, js)*(max(is, js) - 1)/2 + min(is, js)
               w_pt(is, js) = wk(ik, ip)
               crho(is, js) = alpha(ik, ip)*ck(ik, ip)*rho(js)
               rhs(is, js) = weight*vk(ik, ip)
            end do
         end do
         m = -matmul(w_pt, crho)
         do is = 1, ns
            m(is, is) = m(is, is) + 1.0_wp
         end do
         call gauss_jordan(m, rhs)
         do js = 1, ns
            do is = 1, js
               zk(ik, js*(js - 1)/2 + is) = 0.5_wp*(rhs(is, js) + rhs(js, is))
            end do
         end do
      end do
      call trafo%fbt_k2r_all(zk, expected, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, maxval(abs(result - reshape(expected, [npts*npair]))) < &
         & 1.0e-10_wp*(1.0_wp + maxval(abs(expected))), "Block solve must invert M(k)")
      if (allocated(error)) return
      call check(error, maxval(abs(result - v)) > 1.0e-3_wp, "M(k) must differ from the identity")
   end subroutine test_kspace_block_inverse

   !> Activation threshold, refresh schedule, status text and the small-system path
   !>
   !> @param[out] error Test failure
   subroutine test_kspace_activation(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_1d_preconditioner_spectral_type), target :: spectral
      type(moist_model_moz_1d_preconditioner_kspace_type) :: kspace
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_error), allocatable :: err
      real(wp), allocatable :: wk(:, :), cr(:, :), uk(:, :), dcdr(:, :), rho(:)
      real(wp), allocatable :: v(:), filtered(:)

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, w0, err)
      call kspace%bind(spectral)
      call check(error, .not. kspace%ready(0.0_wp) .and. &
         & kspace%describe(0.0_wp) == "kspace, off until RMSR < 1.0E-03")
      if (allocated(error)) return
      kspace%refresh = 2
      call check(error, kspace%wants_update(1.0e-4_wp, 4) .and. .not. kspace%wants_update(1.0e-4_wp, 3) &
         & .and. .not. kspace%wants_update(1.0e-2_wp, 4), "Refresh schedule")
      if (allocated(error)) return

      ! Three sites (six pairs): the factors exist but the spectral filter is used
      call model_inputs(3, rgrid%r, kgrid%r, wk, cr, uk, dcdr, rho)
      call kspace%update(wk, cr, rho, 1.0_wp, uk, dcdr, err)
      call check(error, .not. allocated(err) .and. kspace%ready(1.0e-4_wp) .and. .not. kspace%ready(1.0e-2_wp))
      if (allocated(error)) return
      call check(error, kspace%describe(1.0e-4_wp) == &
         & "kspace, spectral filter for 6 site pairs (block solve above 100)")
      if (allocated(error)) return
      v = test_vector(6)
      filtered = v
      call spectral%apply(filtered)
      call kspace%apply(v)
      call check(error, all(v == filtered), "Few site pairs use the spectral filter")
      if (allocated(error)) return
      call kspace%release()
      call check(error, .not. kspace%valid .and. .not. kspace%ready(0.0_wp))
   end subroutine test_kspace_activation

   !> The per-k block solves give identical results for any team size
   !>
   !> @param[out] error Test failure
   subroutine test_kspace_threads(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ns = 14, npair = ns*(ns + 1)/2
      type(moist_model_moz_1d_preconditioner_spectral_type), target :: spectral
      type(moist_model_moz_1d_preconditioner_kspace_type) :: kspace
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_error), allocatable :: err
      real(wp), allocatable :: wk(:, :), cr(:, :), uk(:, :), dcdr(:, :), rho(:)
      real(wp), allocatable :: one(:), three(:), lu_one(:, :, :)

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, w0, err)
      call kspace%bind(spectral)
      call model_inputs(ns, rgrid%r, kgrid%r, wk, cr, uk, dcdr, rho)
      kspace%nthreads = 1
      call kspace%update(wk, cr, rho, 1.0_wp, uk, dcdr, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      lu_one = kspace%lu
      one = test_vector(npair)
      call kspace%apply(one)
      kspace%nthreads = 3
      call kspace%update(wk, cr, rho, 1.0_wp, uk, dcdr, err)
      three = test_vector(npair)
      call kspace%apply(three)
      call check(error, .not. allocated(err) .and. all(kspace%lu == lu_one) .and. all(one == three), &
         & "Team size must not change the factors or the block solve")
   end subroutine test_kspace_threads

   !> Missing filter and inconsistent shapes are reported
   !>
   !> @param[out] error Test failure
   subroutine test_kspace_update_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_1d_preconditioner_spectral_type), target :: spectral
      type(moist_model_moz_1d_preconditioner_kspace_type) :: kspace
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_error), allocatable :: err
      real(wp), allocatable :: wk(:, :), cr(:, :), uk(:, :), dcdr(:, :), rho(:)

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call model_inputs(3, rgrid%r, kgrid%r, wk, cr, uk, dcdr, rho)
      call kspace%update(wk, cr, rho, 1.0_wp, uk, dcdr, err)
      call check(error, allocated(err) .and. .not. kspace%valid, "Unbound filter")
      if (allocated(error)) return
      call kspace%bind(spectral)
      call kspace%update(wk, cr, rho, 1.0_wp, uk, dcdr, err)
      call check(error, allocated(err), "Unbuilt filter")
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, w0, err)
      call kspace%update(wk, cr, rho(:2), 1.0_wp, uk, dcdr, err)
      call check(error, allocated(err) .and. .not. kspace%valid, "Site count must match the pair count")
      if (allocated(error)) return
      call kspace%update(wk(:npts - 1, :), cr(:npts - 1, :), rho, 1.0_wp, uk(:npts - 1, :), &
         & dcdr(:npts - 1, :), err)
      call check(error, allocated(err) .and. .not. kspace%valid, "Fields must cover the radial grid")
   end subroutine test_kspace_update_errors

   !> Factors follow changes of the site count and of the radial grid
   !>
   !> @param[out] error Test failure
   subroutine test_kspace_reallocation(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_1d_preconditioner_spectral_type), target :: spectral
      type(moist_model_moz_1d_preconditioner_kspace_type) :: kspace
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_error), allocatable :: err
      real(wp), allocatable :: wk(:, :), cr(:, :), uk(:, :), dcdr(:, :), rho(:)

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, w0, err)
      call kspace%bind(spectral)
      call model_inputs(14, rgrid%r, kgrid%r, wk, cr, uk, dcdr, rho)
      call kspace%update(wk, cr, rho, 1.0_wp, uk, dcdr, err)
      call check(error, .not. allocated(err) .and. all(shape(kspace%lu) == [14, 14, npts]))
      if (allocated(error)) return
      call model_inputs(3, rgrid%r, kgrid%r, wk, cr, uk, dcdr, rho)
      call kspace%update(wk, cr, rho, 1.0_wp, uk, dcdr, err)
      call check(error, .not. allocated(err) .and. kspace%ns == 3 .and. all(shape(kspace%lu) == [3, 3, npts]) &
         & .and. all(shape(kspace%lut%idx) == [3, 3]), "A new site count reallocates the factors")
      if (allocated(error)) return

      ! A new radial grid: stale factors are not used until the next update
      call new_uniform_radial_pair(rgrid, kgrid, npts/2, dr, err)
      call spectral%setup(rgrid, kgrid, k0, w0, err)
      call check(error, .not. allocated(err) .and. .not. kspace%ready(0.0_wp), &
         & "Factors for another grid must not be used")
      if (allocated(error)) return
      call model_inputs(3, rgrid%r, kgrid%r, wk, cr, uk, dcdr, rho)
      call kspace%update(wk, cr, rho, 1.0_wp, uk, dcdr, err)
      call check(error, .not. allocated(err) .and. all(shape(kspace%lu) == [3, 3, npts/2]) &
         & .and. kspace%ready(0.0_wp), "A new grid size reallocates the factors")
   end subroutine test_kspace_reallocation

   !> A singular block becomes the identity factorization
   !>
   !> @param[out] error Test failure
   subroutine test_kspace_singular(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ns = 3, npair = ns*(ns + 1)/2
      type(moist_model_moz_1d_preconditioner_spectral_type), target :: spectral
      type(moist_model_moz_1d_preconditioner_kspace_type) :: kspace
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_error), allocatable :: err
      real(wp) :: wk(npts, npair), cr(npts, npair), uk(npts, npair), dcdr(npts, npair)
      real(wp) :: rho(ns)
      integer :: is, ip

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, w0, err)
      call kspace%bind(spectral)
      ! w = 1, c(r) = 0 and alpha = -1, so M(k) = I + diag(rho)*c(k)
      ! with c(k) = -beta*u(k); u(k) = 1/(beta*rho) on the diagonal at k_1 makes
      ! M(k_1) = 0 exactly (powers of two)
      rho = [0.125_wp, 0.25_wp, 0.5_wp]
      wk = 0.0_wp
      uk = 0.0_wp
      do is = 1, ns
         ip = is*(is + 1)/2
         wk(:, ip) = 1.0_wp
         uk(1, ip) = 1.0_wp/(2.0_wp*rho(is))
         uk(2:, ip) = 0.1_wp
      end do
      cr = 0.0_wp
      dcdr = -1.0_wp
      call kspace%update(wk, cr, rho, 2.0_wp, uk, dcdr, err)
      call check(error, .not. allocated(err) .and. kspace%valid)
      if (allocated(error)) return
      call check(error, all(kspace%ipiv(:, 1) == [1, 2, 3]) .and. kspace%lu(1, 1, 1) == 1.0_wp .and. &
         & kspace%lu(2, 2, 1) == 1.0_wp .and. kspace%lu(3, 3, 1) == 1.0_wp .and. &
         & kspace%lu(2, 1, 1) == 0.0_wp, "A singular block must fall back to the identity")
      if (allocated(error)) return
      call check(error, abs(kspace%lu(1, 1, 2) - (1.0_wp - 0.2_wp*rho(1))) < 1.0e-14_wp, &
         & "Regular blocks are factorized")
   end subroutine test_kspace_singular

   !> Releasing or failing to rebuild the filter disables the block solve
   !>
   !> @param[out] error Test failure
   subroutine test_kspace_release_guard(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ns = 14, npair = ns*(ns + 1)/2
      type(moist_model_moz_1d_preconditioner_spectral_type), target :: spectral
      type(moist_model_moz_1d_preconditioner_kspace_type) :: kspace
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_error), allocatable :: err
      real(wp), allocatable :: wk(:, :), cr(:, :), uk(:, :), dcdr(:, :), rho(:), v(:)

      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call spectral%setup(rgrid, kgrid, k0, w0, err)
      call kspace%bind(spectral)
      call model_inputs(ns, rgrid%r, kgrid%r, wk, cr, uk, dcdr, rho)
      call kspace%update(wk, cr, rho, 1.0_wp, uk, dcdr, err)
      call check(error, .not. allocated(err) .and. kspace%ready(0.0_wp) .and. &
         & kspace%describe(0.0_wp) == "kspace, block solve on 105 site pairs")
      if (allocated(error)) return
      ! A failed re-setup releases the filter first
      call spectral%setup(rgrid, kgrid, -1.0_wp, w0, err)
      call check(error, allocated(err) .and. .not. kspace%ready(0.0_wp) .and. &
         & kspace%describe(0.0_wp) == "kspace, off until RMSR < 1.0E-03")
      if (allocated(error)) return
      v = test_vector(npair)
      call kspace%apply(v)
      call check(error, all(v == test_vector(npair)), "An unbuilt filter leaves the vector unchanged")
      if (allocated(error)) return
      call spectral%release()
      call kspace%apply(v)
      call check(error, all(v == test_vector(npair)))
   end subroutine test_kspace_release_guard

   !> Name lists and selection of the lent operator
   !>
   !> @param[out] error Test failure
   subroutine test_selection(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_1d_preconditioner_spectral_type), target :: spectral
      type(moist_model_moz_1d_preconditioner_kspace_type), target :: kspace
      class(moist_model_moz_solver_preconditioner_type), pointer :: precond

      call check(error, any(preconditioners_1d_vv == "kspace") .and. .not. any(preconditioners_1d_uv == "kspace") &
         & .and. all(preconditioners_1d_uv(:2) == preconditioners_1d_vv(:2)))
      if (allocated(error)) return
      call select_1d_preconditioner("none", spectral, precond, kspace)
      call check(error, .not. associated(precond))
      if (allocated(error)) return
      call select_1d_preconditioner(" Spectral ", spectral, precond)
      call check(error, associated(precond, spectral))
      if (allocated(error)) return
      call select_1d_preconditioner("KSPACE", spectral, precond, kspace)
      call check(error, associated(precond, kspace) .and. associated(kspace%spectral, spectral))
      if (allocated(error)) return
      call select_1d_preconditioner("kspace", spectral, precond)
      call check(error, .not. associated(precond), "kspace needs the model-owned block inverse")
   end subroutine test_selection

   !> The block solve is known from the site-pair count before any update
   !>
   !> - Water (3 sites, 6 pairs) and the 100-pair boundary take the spectral filter
   !> - 14 sites (105 pairs) use the block solve
   !>
   !> @param[out] error Test failure
   subroutine test_kspace_block_query(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Unbuilt block inverse
      type(moist_model_moz_1d_preconditioner_kspace_type) :: kspace

      call check(error, .not. kspace%uses_block_solve(6) .and. .not. kspace%uses_block_solve(100) &
         & .and. kspace%uses_block_solve(101) .and. kspace%uses_block_solve(105), &
         & "The block solve needs more than 100 site pairs")
   end subroutine test_kspace_block_query

   !> Solve a X = b in place by Gauss-Jordan elimination with partial pivoting
   !>
   !> @param[inout] a Matrix (n, n), destroyed
   !> @param[inout] b Right-hand sides (n, m), replaced by the solution
   pure subroutine gauss_jordan(a, b)
      !> Matrix, destroyed
      real(wp), intent(inout) :: a(:, :)
      !> Right-hand sides, replaced by the solution
      real(wp), intent(inout) :: b(:, :)
      real(wp), allocatable :: row(:), brow(:)
      integer :: n, i, p

      n = size(a, 1)
      do i = 1, n
         p = i - 1 + maxloc(abs(a(i:n, i)), 1)
         if (p /= i) then
            row = a(i, :)
            a(i, :) = a(p, :)
            a(p, :) = row
            brow = b(i, :)
            b(i, :) = b(p, :)
            b(p, :) = brow
         end if
         b(i, :) = b(i, :)/a(i, i)
         a(i, :) = a(i, :)/a(i, i)
         do p = 1, n
            if (p == i) cycle
            b(p, :) = b(p, :) - a(p, i)*b(i, :)
            a(p, :) = a(p, :) - a(p, i)*a(i, :)
         end do
      end do
   end subroutine gauss_jordan

   !> Damp only updates the solver allows to post-process, with the damping filter's own weights
   !>
   !> - PICARD and DIIS: damped when enabled, unchanged when disabled
   !> - JFNK: never damped; HYBRID: damped in the MDIIS phase, not after the switch to JFNK
   !> - Damping and spectral filters differ, so applying the wrong one fails the check
   !> @param[out] error Test failure
   subroutine test_filters_damping(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_1d_solver_filters_type) :: filters, undamped
      type(moist_model_moz_1d_preconditioner_spectral_type) :: reference, spectral
      type(moist_model_moz_solver_parameters_type) :: parameters
      class(moist_model_moz_solver_type), allocatable :: solver
      type(identity_residual_type) :: evaluator
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_context_type), target :: ctx
      type(moist_error), allocatable :: err
      real(wp), allocatable :: gamma(:, :), proposed(:, :), damped(:, :), expected(:, :), wrong(:, :), delta(:)
      real(wp), allocatable :: updated(:, :)
      character(len=6), parameter :: postprocessed(2) = [character(len=6) :: "PICARD", "DIIS"]
      integer, parameter :: npair = 3
      real(wp), parameter :: damp_k0 = 2.0_wp*k0, damp_weight = 0.5_wp
      real(wp) :: rmse, scale
      logical :: feasible
      integer :: i

      call new_context(ctx, nthreads=0, verbosity=0)
      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      call filters%setup(rgrid, kgrid, 1, "spectral", k0, w0, .true., damp_k0, damp_weight, err)
      if (.not. allocated(err)) call undamped%setup(rgrid, kgrid, 1, "spectral", k0, w0, .false., damp_k0, &
         & damp_weight, err)
      if (.not. allocated(err)) call reference%setup(rgrid, kgrid, damp_k0, damp_weight, err)
      if (.not. allocated(err)) call spectral%setup(rgrid, kgrid, k0, w0, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      gamma = reshape(0.5_wp*test_vector(npair), [npts, npair])
      proposed = reshape(test_vector(npair), [npts, npair])
      delta = reshape(proposed - gamma, [npts*npair])
      call reference%apply(delta)
      expected = gamma + reshape(delta, [npts, npair])
      delta = reshape(proposed - gamma, [npts*npair])
      call spectral%apply(delta)
      wrong = gamma + reshape(delta, [npts, npair])
      scale = 1.0e-12_wp*(1.0_wp + maxval(abs(expected)))
      call check(error, maxval(abs(expected - wrong)) > 1.0e-6_wp, "The test needs distinct filters")
      if (allocated(error)) return

      do i = 1, size(postprocessed)
         parameters%name = postprocessed(i)
         call new_solver(solver, ctx, err, parameters)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         damped = proposed
         call filters%damp(solver, gamma, damped)
         call check(error, maxval(abs(damped - expected)) < scale, &
            & "Enabled damping must apply the damping filter: "//trim(postprocessed(i)))
         if (allocated(error)) return
         damped = proposed
         call undamped%damp(solver, gamma, damped)
         call check(error, all(damped == proposed), "Disabled damping must keep the update: "//trim(postprocessed(i)))
         if (allocated(error)) return
      end do

      parameters%name = "JFNK"
      call new_solver(solver, ctx, err, parameters)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      damped = proposed
      call filters%damp(solver, gamma, damped)
      call check(error, all(damped == proposed), "JFNK updates must not be damped")
      if (allocated(error)) return

      parameters%name = "HYBRID"
      parameters%hybrid%switch_rmse = 1.0e3_wp
      parameters%jfnk%use_exact_jvp = .false.
      call new_solver(solver, ctx, err, parameters)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      damped = proposed
      call filters%damp(solver, gamma, damped)
      call check(error, maxval(abs(damped - expected)) < scale, "HYBRID must damp in its MDIIS phase")
      if (allocated(error)) return
      allocate (updated, mold=gamma)
      call solver%outer_step(evaluator, gamma, updated, rmse, feasible, err)
      call check(error, feasible .and. .not. allocated(err) .and. .not. solver%allows_postprocessing(), &
         & "HYBRID must switch to JFNK below switch_rmse")
      if (allocated(error)) return
      damped = proposed
      call filters%damp(solver, gamma, damped)
      call check(error, all(damped == proposed), "HYBRID must not damp after switching to JFNK")
      call ctx%delete()
   end subroutine test_filters_damping

   !> Lend the selected preconditioner and always the fallback filter; detach keeps both built
   !> @param[out] error Test failure
   subroutine test_filters_attach(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_model_moz_1d_solver_filters_type), target :: filters
      type(moist_model_moz_solver_jfnk_type) :: jfnk
      class(moist_model_moz_solver_preconditioner_type), pointer :: precond
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_context_type), target :: ctx
      type(moist_error), allocatable :: err

      call new_context(ctx, nthreads=0, verbosity=0)
      call make_pair(rgrid, kgrid, error)
      if (allocated(error)) return
      jfnk = new_jfnk_solver(ctx)

      call filters%setup(rgrid, kgrid, 1, "none", k0, w0, .false., k0, w0, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call filters%attach(jfnk)
      call check(error, .not. associated(jfnk%precond) .and. associated(jfnk%smoother, filters%spectral), &
         & "Without a preconditioner the fallback filter must still be lent")
      if (allocated(error)) return
      precond => filters%response_preconditioner()
      call check(error, .not. associated(precond), "No response preconditioner for none")
      if (allocated(error)) return
      call filters%detach(jfnk)
      call check(error, .not. associated(jfnk%precond) .and. .not. associated(jfnk%smoother), &
         & "Detach must take both operators back")
      if (allocated(error)) return

      call filters%setup(rgrid, kgrid, 1, "Spectral", k0, w0, .false., k0, w0, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call filters%attach(jfnk)
      call check(error, associated(jfnk%precond, filters%spectral) .and. associated(jfnk%smoother, filters%spectral), &
         & "The spectral filter must precondition and filter the fallback")
      if (allocated(error)) return
      call filters%detach(jfnk)
      precond => filters%response_preconditioner()
      call check(error, associated(precond, filters%spectral) .and. filters%spectral%ready(1.0_wp) .and. &
         & .not. associated(jfnk%precond) .and. .not. associated(jfnk%smoother), &
         & "Detached filters must stay built for the response")
      if (allocated(error)) return

      call filters%setup(rgrid, kgrid, 1, "kspace", k0, w0, .false., k0, w0, err)
      call check(error, allocated(err), "1D UV filters must refuse the VV-only kspace")
      call ctx%delete()
   end subroutine test_filters_attach

   !> Residual R(gamma) = gamma
   !>
   !> @param[in]  self      evaluator
   !> @param[in]  gamma     iterate
   !> @param[out] residual  gamma
   !> @param[out] feasible  always true
   !> @param[out] error     never allocated
   subroutine identity_residual(self, gamma, residual, feasible, error)
      class(identity_residual_type), intent(in) :: self
      real(wp), intent(in) :: gamma(:, :)
      real(wp), intent(out) :: residual(:, :)
      logical, intent(out) :: feasible
      type(moist_error), allocatable, intent(out) :: error

      associate (unused => self); end associate
      residual = gamma
      feasible = .true.
   end subroutine identity_residual

end module test_moz_1d_preconditioner
