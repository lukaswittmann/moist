!> Electronic reverse path of the DROP cavity
!>
!> - [[contract_surface_lsf_weights]] folds a surface adjoint into level-set
!>   value, gradient and Hessian adjoints
!> - [[get_surface_response_drop]] scales them by the density chain-rule factor
!>   and accumulates a density item
!>
!> The weights are pinned against the forward nuclear Jacobian
!>
!> - An atom A that does not own grid point i enters the surface only through
!>   the level-set field
!> - The adjoint jet contracted with the nuclear derivatives of the level-set
!>   jet at the projected point must then reproduce, for every axis beta
!>
!>   w_xyz_i . xyz1_rA + w_xi_i * xi1_rA + w_a_i * a_i1_rA + w_w_i * wleb1_rA
!>   + w_n_i . normal1_rA + w_k1_i * k1_rA + w_k2_i * k2_rA
!>
!> - The owner atom is excluded, its anchor term is not part of the level-set
!>   adjoint
!> - `w_f` is excluded, switching is anchor-only
!> - The forward area derivative carries a switching part, `a = R^2 f wleb`
!>   with `f` moving with every neighbor of the anchor; that part belongs to
!>   `w_f`, so the reference removes `R^2 wleb f1_rA`
!> - The two sides are separate passes, the forward sensitivity pass against
!>   the level set's own nuclear jet, but both go through `drop_kkt_solve`, the
!>   packed contraction and the objective jet, so a defect shared by those
!>   (e.g. a swapped packed Hessian unpack) cancels here
!>
!> Error paths: the per-point failures of the tiled reverse loops, shared by
!> the response and the nuclear-gradient path, are checked on both; a host
!> callback density fails or flattens on demand at one grid point; an unbuilt
!> cavity is refused by every entry point, wrongly shaped weights by the
!> contraction
module test_cavity_drop_potential
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_funloc, c_loc, &
      & c_f_pointer, c_associated
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, &
      & test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_lsf_isodensity_internal, only: &
      & moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_drop_lsf_isodensity_callback, only: &
      & moist_cavity_drop_lsf_isodensity_callback_type
   use moist_cavity_drop_lsf_isodensity_param, only: &
      & moist_cavity_drop_lsf_isodensity_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_seed_result_type, build_seed_state, apply_seed, seed_state_ok
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   use moist_cavity_adjoint, only: cavity_adjoint_type
   use moist_channels_coupling, only: coupling_type, density_matrix_request_type, &
      & nuclear_charge_request_type
   use moist_channels_response, only: response_type, response_item_type, &
      & density_response_type, response_accumulate, current_response_item
   use moist_context, only: moist_context_type, new_context
   use moist_math_packed, only: packed_index
   use moist_radii, only: default_cpcm_radii, new_radii_custom_atoms, radius_type
   use test_helpers, only: fill_legacy_radii
   implicit none(type, external)
   private

   public :: collect_cavity_drop_potential

   !> SvdW blend sharpness of the fixture
   real(wp), parameter :: BLEND_K = 2.5_wp
   !> SvdW three-body blend of the fixture
   real(wp), parameter :: BLEND_3B = 1.0_wp
   !> Lebedev points per atom of the SvdW fixture
   integer, parameter :: NUM_LEB = 302
   !> Projection tolerance
   real(wp), parameter :: PROJ_TOL = 1.0E-14_wp
   !> Projection iteration cap
   integer, parameter :: PROJ_MAXITER = 1000
   !> Projection level
   integer, parameter :: PROJ_LEVEL = 1

   !> Absolute bound on the contraction against the forward sum
   !>
   !> - Measured worst deviation 2.2e-15 over 4152 entries, all channels and
   !>   each channel alone, against adjoint magnitudes up to ~2
   real(wp), parameter :: EQ_ABS = 1.0E-12_wp
   !> Relative bound on each forward reference entry
   real(wp), parameter :: EQ_REL = 1.0E-12_wp

   !> Density-item bound against the direct contraction
   !>
   !> - Relative to the larger of one and the entry size
   !> - Measured worst absolute deviation 6.8e-13 against entries up to ~9.1e3
   !>   (isodensity fixture at 302 points per atom), relative ~1e-16
   real(wp), parameter :: ITEM_TOL = 1.0E-12_wp

   !> Round-off bound for the switching channel not moving the adjoint
   !>
   !> - Measured deviation exactly zero
   real(wp), parameter :: SWITCH_TOL = 1.0E-12_wp

   !> Near-umbilic guard, principal curvatures are not differentiable at k1 == k2
   real(wp), parameter :: CURV_GAP_THR = 5.0E-2_wp
   !> Reference entries below this carry no weight in the vacuity guard
   real(wp), parameter :: VACUITY_THR = 1.0E-6_wp
   !> Minimum compared entries (measured 1842 to 4152)
   integer, parameter :: MIN_COMPARED = 1000
   !> Minimum live entries above the vacuity threshold (measured 1841 to 4066)
   integer, parameter :: MIN_LIVE = 500
   !> Violation (deviation over tolerance) the broken contraction must exceed
   !>
   !> - Measured 3e8 and 2e8 with the former tolerance of 5e-10
   real(wp), parameter :: CONTROL_MARGIN = 1.0E3_wp
   !> Relative bound on the linearity of the symmetric Hessian seed
   !>
   !> - Measured worst deviation 1.7e-16 against responses up to ~0.5
   !> - Control gap between the single-entry seeds measured 0.089 to 0.24
   real(wp), parameter :: SEED_LIN_TOL = 1.0E-13_wp

   !> Grid point the failure tests break: column 6 of the second surface tile
   integer, parameter :: BAD_POINT = 70
   !> First grid point of the tile holding [[BAD_POINT]]
   integer, parameter :: BAD_TILE_START = 65

   !> Host density behaves normally
   integer, parameter :: HOST_OK = 0
   !> Host callback reports failure at the broken point
   integer, parameter :: HOST_FAIL = 1
   !> Host density has vanishing derivatives at the broken point
   integer, parameter :: HOST_FLAT = 2

   !> Host side of the callback level set: two Gaussian centers
   type :: gauss_host_type
      !> Gaussian centers (3, 2)
      real(wp) :: center(3, 2) = 0.0_wp
      !> Gaussian exponents
      real(wp) :: expo(2) = 0.0_wp
      !> Behaviour at `bad_xyz`, one of the HOST_* modes
      integer :: mode = HOST_OK
      !> Point the failure modes act on
      real(wp) :: bad_xyz(3) = 0.0_wp
   end type gauss_host_type

   !> Channel identifier of w_xyz
   integer, parameter :: CH_XYZ = 1
   !> Channel identifier of w_xi
   integer, parameter :: CH_XI = 2
   !> Channel identifier of w_a
   integer, parameter :: CH_A = 3
   !> Channel identifier of w_w
   integer, parameter :: CH_W = 4
   !> Channel identifier of w_n
   integer, parameter :: CH_N = 5
   !> Channel identifier of w_k1
   integer, parameter :: CH_K1 = 6
   !> Channel identifier of w_k2
   integer, parameter :: CH_K2 = 7
   !> Channel identifier of w_f
   integer, parameter :: CH_F = 8

   !> Every channel that reaches the level set
   integer, parameter :: ALL_CHANNELS(7) = [CH_XYZ, CH_XI, CH_A, CH_W, CH_N, CH_K1, CH_K2]

   !> Outcome of one contraction-versus-forward comparison
   type :: compare_stats_type
      !> Entries compared
      integer :: ncompared = 0
      !> Entries whose forward reference exceeds the vacuity threshold
      integer :: nlive = 0
      !> Entries outside the tolerance, NaN included
      integer :: nfail = 0
      !> Largest absolute deviation
      real(wp) :: max_diff = 0.0_wp
      !> Largest deviation over its tolerance
      real(wp) :: max_viol = 0.0_wp
   end type compare_stats_type

contains

   !> Register the electronic reverse-path tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_cavity_drop_potential(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lsf_weights_vs_forward", test_lsf_weights_vs_forward), &
                  new_unittest("lsf_weights_single_channels", &
                               test_lsf_weights_single_channels), &
                  new_unittest("lsf_weights_offdiagonal_control", &
                               test_lsf_weights_offdiagonal_control), &
                  new_unittest("symmetric_hessian_seed", test_symmetric_hessian_seed), &
                  new_unittest("lsf_weights_ignore_switching", &
                               test_lsf_weights_ignore_switching), &
                  new_unittest("density_response_svdw_untouched", &
                               test_density_response_svdw_untouched), &
                  new_unittest("density_response_isodensity", &
                               test_density_response_isodensity), &
                  new_unittest("bad_accumulator_uninitialized", &
                               test_bad_accumulator_uninitialized, should_fail=.true.), &
                  new_unittest("bad_accumulator_size", test_bad_accumulator_size, &
                               should_fail=.true.), &
                  new_unittest("bad_response_unbuilt", test_bad_response_unbuilt, &
                               should_fail=.true.), &
                  new_unittest("bad_lsf_weights_unbuilt", test_bad_lsf_weights_unbuilt, &
                               should_fail=.true.), &
                  new_unittest("bad_gradient_unbuilt", test_bad_gradient_unbuilt, &
                               should_fail=.true.), &
                  new_unittest("bad_lsf_weights_shape", test_bad_lsf_weights_shape, &
                               should_fail=.true.), &
                  new_unittest("bad_response_unreached_point", test_bad_response_unreached, &
                               should_fail=.true.), &
                  new_unittest("bad_gradient_unreached_point", test_bad_gradient_unreached, &
                               should_fail=.true.), &
                  new_unittest("bad_response_host_failure", test_bad_response_host_failure, &
                               should_fail=.true.), &
                  new_unittest("bad_gradient_host_failure", test_bad_gradient_host_failure, &
                               should_fail=.true.), &
                  new_unittest("bad_response_degenerate_point", test_bad_response_degenerate, &
                               should_fail=.true.), &
                  new_unittest("bad_gradient_degenerate_point", test_bad_gradient_degenerate, &
                               should_fail=.true.) &
                  ]

   end subroutine collect_cavity_drop_potential

   !* ============================================================================== *!
   !*                              Contraction identity                               *!
   !* ============================================================================== *!

   !> Every contributing channel at once
   !>
   !> @param[out] error  test failure
   subroutine test_lsf_weights_vs_forward(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx

      call build_svdw_cavity(cavity, ctx, .true., error)
      if (allocated(error)) return
      call run_identity(cavity, ALL_CHANNELS, "all channels", error)

   end subroutine test_lsf_weights_vs_forward

   !> Each contributing channel on its own, so one cannot hide behind another
   !>
   !> @param[out] error  test failure
   subroutine test_lsf_weights_single_channels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      character(len=5), parameter :: labels(7) = [character(len=5) :: &
                                                 "w_xyz", "w_xi", "w_a", "w_w", "w_n", &
                                                 "w_k1", "w_k2"]
      integer :: ich

      call build_svdw_cavity(cavity, ctx, .true., error)
      if (allocated(error)) return
      do ich = 1, size(ALL_CHANNELS)
         call run_identity(cavity, [ALL_CHANNELS(ich)], trim(labels(ich)), error)
         if (allocated(error)) return
      end do

   end subroutine test_lsf_weights_single_channels

   !> Halving the off-diagonal Hessian adjoint must break the identity
   !>
   !> - A wrong packed-unpack convention (doubled or halved off-diagonal)
   !>   changes exactly this part of the result
   !> - The Hessian adjoint is exactly symmetric, one packed seed fills both
   !>   off-diagonal slots; it must also be non-vanishing
   !>
   !> @param[out] error  test failure
   subroutine test_lsf_weights_offdiagonal_control(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)
      type(compare_stats_type) :: stats
      integer, parameter :: sets(2, 2) = reshape([CH_K1, CH_K2, CH_N, CH_K1], [2, 2])
      real(wp) :: asym, size_w2
      integer :: iset, igrid, a, b

      call build_svdw_cavity(cavity, ctx, .true., error)
      if (allocated(error)) return

      do iset = 1, 2
         call fill_channels(acc, cavity%ngrid, sets(:, iset))
         call lsf_weights(cavity, acc, w0, w1, w2, error)
         if (allocated(error)) return
         call compare_to_forward(cavity, acc, w0, w1, w2, .true., stats, error)
         if (allocated(error)) return
         call require_healthy(stats, "intact Hessian adjoint, set "// &
                              to_string(iset), error)
         if (allocated(error)) return

         asym = 0.0_wp
         size_w2 = maxval(abs(w2))
         do igrid = 1, cavity%ngrid
            do a = 1, 3
               do b = 1, 3
                  asym = max(asym, abs(w2(a, b, igrid) - w2(b, a, igrid)))
               end do
            end do
         end do
         if (.not. (asym == 0.0_wp .and. size_w2 > VACUITY_THR)) then
            call test_failed(error, "Hessian adjoint is asymmetric or vanishing, set "// &
                             to_string(iset)//": max asymmetry "//to_string(asym)// &
                             ", max entry "//to_string(size_w2))
            return
         end if

         do igrid = 1, cavity%ngrid
            do a = 1, 3
               do b = 1, 3
                  if (a /= b) w2(a, b, igrid) = 0.5_wp*w2(a, b, igrid)
               end do
            end do
         end do
         call compare_to_forward(cavity, acc, w0, w1, w2, .true., stats, error)
         if (allocated(error)) return
         if (.not. (stats%max_viol > CONTROL_MARGIN)) then
            call test_failed(error, "halved off-diagonal Hessian adjoint was not "// &
                             "detected, set "//to_string(iset)//": violation "// &
                             to_string(stats%max_viol)//" is below the margin "// &
                             to_string(CONTROL_MARGIN))
            return
         end if
      end do

   end subroutine test_lsf_weights_offdiagonal_control

   !> Symmetric Hessian seed is the mean of the two single-entry seeds
   !>
   !> - `apply_seed` is linear in its seed, so for every off-diagonal pair the
   !>   response to `(E_ab + E_ba)/2` must be the mean of the `E_ab` and `E_ba`
   !>   responses in `dxi`, `dk1` and `dk2`, the fields the jet adjoint reads
   !> - Negative control: the `E_ab` and `E_ba` responses differ far beyond the
   !>   tolerance, so the symmetric seed changes the adjoint
   !> - Generic curved, non-umbilic point with live switching slopes
   !>
   !> @param[out] error  test failure
   subroutine test_symmetric_hessian_seed(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_seed_state_type) :: state
      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc, f_wleb
      type(drop_seed_result_type) :: r_ab, r_ba, r_sym
      integer, parameter :: pairs(2, 3) = reshape([1, 2, 1, 3, 2, 3], [2, 3])
      real(wp) :: dlsf1_r(3), dlsf2_rr(3, 3), dr(3)
      real(wp) :: v_ab(3), v_ba(3), v_sym(3), scale, dev, gap
      integer :: status, ipair, a, b

      state%lsf1_r = [0.3_wp, -0.5_wp, 0.8_wp]
      state%lsf2_rr = [0.7_wp, 0.2_wp, -0.4_wp, 0.15_wp, 0.35_wp, 0.5_wp]
      state%lambda_val = 0.6_wp
      state%alpha_coeff = 1.0_wp
      state%anchor = [0.4_wp, 0.1_wp, 1.2_wp]
      state%owner_xyz = [0.0_wp, 0.0_wp, 0.0_wp]
      state%anchor_wleb0 = 0.05_wp
      state%cpjac_scal0 = 1.3_wp
      state%w_f0 = 0.8_wp
      state%wbranch = 1.0_wp
      state%wleb = 0.04_wp
      state%xi0 = 2.0_wp
      state%want_curvature = .true.
      call new_smoothstep_sigmoid(f_crit, 0.5_wp, 1.5_wp)
      call new_smoothstep_sigmoid(f_foc, 0.0_wp, 2.0_wp)
      call new_smoothstep_sigmoid(f_wleb, 0.0_wp, 1.0_wp)

      call build_seed_state(state, f_crit, f_foc, f_wleb, .false., status)
      if (status /= seed_state_ok) then
         call test_failed(error, "seed state is degenerate, status "//to_string(status))
         return
      end if

      dlsf1_r = 0.0_wp
      dr = 0.0_wp
      do ipair = 1, size(pairs, 2)
         a = pairs(1, ipair)
         b = pairs(2, ipair)

         dlsf2_rr = 0.0_wp
         dlsf2_rr(a, b) = 1.0_wp
         call apply_seed(state, dlsf1_r, dlsf2_rr, dr, 0.0_wp, r_ab)
         dlsf2_rr = 0.0_wp
         dlsf2_rr(b, a) = 1.0_wp
         call apply_seed(state, dlsf1_r, dlsf2_rr, dr, 0.0_wp, r_ba)
         dlsf2_rr = 0.0_wp
         dlsf2_rr(a, b) = 0.5_wp
         dlsf2_rr(b, a) = 0.5_wp
         call apply_seed(state, dlsf1_r, dlsf2_rr, dr, 0.0_wp, r_sym)

         v_ab = [r_ab%dxi, r_ab%dk1, r_ab%dk2]
         v_ba = [r_ba%dxi, r_ba%dk1, r_ba%dk2]
         v_sym = [r_sym%dxi, r_sym%dk1, r_sym%dk2]
         scale = max(maxval(abs(v_ab)), maxval(abs(v_ba)))
         dev = maxval(abs(v_sym - 0.5_wp*(v_ab + v_ba)))
         gap = maxval(abs(v_ab - v_ba))

         if (.not. (dev <= SEED_LIN_TOL*scale)) then
            call test_failed(error, "symmetric seed is not the mean of the single-entry "// &
                             "seeds, pair "//to_string(a)//to_string(b)//": deviation "// &
                             to_string(dev)//", scale "//to_string(scale))
            return
         end if
         if (.not. (gap > CONTROL_MARGIN*SEED_LIN_TOL*scale)) then
            call test_failed(error, "control: single-entry seeds agree, pair "// &
                             to_string(a)//to_string(b)//": gap "//to_string(gap)// &
                             ", scale "//to_string(scale))
            return
         end if
      end do

   end subroutine test_symmetric_hessian_seed

   !> The switching channel never reaches the level set
   !>
   !> - Switching is anchor-only, so adding `w_f` changes nothing beyond round-off
   !> - Negative control: perturbing the width channel does change the result
   !>
   !> @param[out] error  test failure
   subroutine test_lsf_weights_ignore_switching(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      real(wp), allocatable :: a0(:), a1(:, :), a2(:, :, :)
      real(wp), allocatable :: b0(:), b1(:, :), b2(:, :, :)
      integer, parameter :: base(3) = [CH_XYZ, CH_XI, CH_N]
      real(wp) :: dev_f, dev_xi

      call build_svdw_cavity(cavity, ctx, .false., error)
      if (allocated(error)) return

      call fill_channels(acc, cavity%ngrid, base)
      call lsf_weights(cavity, acc, a0, a1, a2, error)
      if (allocated(error)) return

      call fill_channels(acc, cavity%ngrid, [base, CH_F])
      if (.not. any(abs(acc%w_f) > VACUITY_THR)) then
         call test_failed(error, "w_f channel was not populated")
         return
      end if
      call lsf_weights(cavity, acc, b0, b1, b2, error)
      if (allocated(error)) return
      dev_f = max(maxval(abs(a0 - b0)), maxval(abs(a1 - b1)), maxval(abs(a2 - b2)))

      call fill_channels(acc, cavity%ngrid, base)
      acc%w_xi = 1.01_wp*acc%w_xi
      call lsf_weights(cavity, acc, b0, b1, b2, error)
      if (allocated(error)) return
      dev_xi = max(maxval(abs(a0 - b0)), maxval(abs(a1 - b1)), maxval(abs(a2 - b2)))

      if (.not. (maxval(abs(a1)) > VACUITY_THR)) then
         call test_failed(error, "reference adjoint is vacuous")
         return
      end if
      if (.not. (dev_f <= SWITCH_TOL*max(1.0_wp, maxval(abs(a1))))) then
         call test_failed(error, "w_f changed the level-set adjoint by "// &
                          to_string(dev_f))
         return
      end if
      if (.not. (dev_xi > 1.0E-6_wp*maxval(abs(a1)))) then
         call test_failed(error, "control: a 1% change of w_xi moved the adjoint "// &
                          "by only "// &
                          to_string(dev_xi))
         return
      end if

   end subroutine test_lsf_weights_ignore_switching

   !> Contraction identity for one channel selection, with its guards
   !>
   !> @param[in]  cavity    Built cavity with the forward arrays filled
   !> @param[in]  channels  Channel identifiers to populate
   !> @param[in]  label     Channel description for failure messages
   !> @param[out] error     Test failure
   subroutine run_identity(cavity, channels, label, error)
      !> Built cavity with the forward arrays filled
      type(cavity_type_drop), intent(in) :: cavity
      !> Channel identifiers to populate
      integer, intent(in) :: channels(:)
      !> Channel description for failure messages
      character(len=*), intent(in) :: label
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cavity_adjoint_type) :: acc
      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)
      type(compare_stats_type) :: stats

      call fill_channels(acc, cavity%ngrid, channels)
      call lsf_weights(cavity, acc, w0, w1, w2, error)
      if (allocated(error)) return
      call compare_to_forward(cavity, acc, w0, w1, w2, &
                              any(channels == CH_K1 .or. channels == CH_K2), stats, error)
      if (allocated(error)) return
      call require_healthy(stats, label, error)

   end subroutine run_identity

   !> Fail on a deviation, or on a comparison too small to mean anything
   !>
   !> @param[in]  stats  Comparison outcome
   !> @param[in]  label  Channel description for failure messages
   !> @param[out] error  Test failure
   subroutine require_healthy(stats, label, error)
      !> Comparison outcome
      type(compare_stats_type), intent(in) :: stats
      !> Channel description for failure messages
      character(len=*), intent(in) :: label
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      if (stats%nfail > 0) then
         call test_failed(error, "contraction/forward mismatch for "//label//": "// &
                          to_string(stats%nfail)//" of "//to_string(stats%ncompared)// &
                          " entries outside tolerance, max deviation "// &
                          to_string(stats%max_diff))
         return
      end if
      if (stats%ncompared < MIN_COMPARED .or. stats%nlive < MIN_LIVE) then
         call test_failed(error, "comparison for "//label//" is vacuous: "// &
                          to_string(stats%ncompared)//" entries compared, "// &
                          to_string(stats%nlive)//" above the vacuity threshold")
      end if

   end subroutine require_healthy

   !> Compare the level-set adjoint contraction with the forward sum
   !>
   !> @param[in]  cavity       Built cavity with the forward arrays filled
   !> @param[in]  acc          Surface adjoints that produced the weights
   !> @param[in]  w0           LSF value adjoints (ngrid)
   !> @param[in]  w1           LSF gradient adjoints (3, ngrid)
   !> @param[in]  w2           LSF Hessian adjoints (3, 3, ngrid)
   !> @param[in]  skip_umbilic Skip near-umbilic points (curvature channels present)
   !> @param[out] stats        Comparison outcome
   !> @param[out] error        Test failure
   subroutine compare_to_forward(cavity, acc, w0, w1, w2, skip_umbilic, stats, error)
      !> Built cavity with the forward arrays filled
      type(cavity_type_drop), intent(in) :: cavity
      !> Surface adjoints that produced the weights
      type(cavity_adjoint_type), intent(in) :: acc
      !> LSF value adjoints
      real(wp), intent(in) :: w0(:)
      !> LSF gradient adjoints
      real(wp), intent(in) :: w1(:, :)
      !> LSF Hessian adjoints
      real(wp), intent(in) :: w2(:, :, :)
      !> Skip near-umbilic points
      logical, intent(in) :: skip_umbilic
      !> Comparison outcome
      type(compare_stats_type), intent(out) :: stats
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      type(mctc_error), allocatable :: lsf_err
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :), lsf3_rr_rA(:, :, :)
      real(wp) :: rev, fwd, diff, viol
      integer :: igrid, iatom, iaxis, islot, a, b, nsph

      nsph = cavity%nsph
      allocate (lsf, source=cavity%lsf_model)
      call lsf%set_max_deriv(3)
      ! Active-indexed outputs; the atom count bounds the active count at every point
      allocate (lsf1_rA(3, nsph), lsf2_r_rA(3, 3, nsph), lsf3_rr_rA(6, 3, nsph))

      do igrid = 1, cavity%ngrid
         if (skip_umbilic) then
            if (abs(cavity%k1(igrid) - cavity%k2(igrid)) < CURV_GAP_THR) cycle
         end if
         call lsf%prepare(cavity%xyz(:, igrid), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "LSF prepare failed: "//lsf_err%message)
            return
         end if
         call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)

         do iatom = 1, nsph
            if (iatom == cavity%owner(igrid)) cycle
            do iaxis = 1, 3
               ! Atoms dropped by screening keep a zero contraction
               rev = 0.0_wp
               do islot = 1, lsf%active_count()
                  if (lsf%active_atom(islot) /= iatom) cycle
                  rev = w0(igrid)*lsf1_rA(iaxis, islot)
                  do a = 1, 3
                     rev = rev + w1(a, igrid)*lsf2_r_rA(a, iaxis, islot)
                     do b = 1, 3
                        rev = rev + w2(a, b, igrid) &
                              *lsf3_rr_rA(packed_index(a, b), iaxis, islot)
                     end do
                  end do
               end do

               fwd = acc%w_xi(igrid)*cavity%xi1_rA(iaxis, iatom, igrid) &
                     + acc%w_a(igrid)*(cavity%a_i1_rA(iaxis, iatom, igrid) &
                                       - cavity%radii(cavity%owner(igrid))**2 &
                                       *cavity%wleb(igrid) &
                                       *cavity%f1_rA(iaxis, iatom, igrid)) &
                     + acc%w_w(igrid)*cavity%wleb1_rA(iaxis, iatom, igrid) &
                     + acc%w_k1(igrid)*cavity%k1_rA(iaxis, iatom, igrid) &
                     + acc%w_k2(igrid)*cavity%k2_rA(iaxis, iatom, igrid) &
                     + dot_product(acc%w_xyz(:, igrid), &
                                   cavity%xyz1_rA(:, iaxis, iatom, igrid)) &
                     + dot_product(acc%w_n(:, igrid), &
                                   cavity%normal1_rA(:, iatom, iaxis, igrid))

               diff = abs(rev - fwd)
               viol = diff/max(EQ_ABS, EQ_REL*abs(fwd))
               stats%ncompared = stats%ncompared + 1
               if (abs(fwd) > VACUITY_THR) stats%nlive = stats%nlive + 1
               ! Negated pass test, so a NaN on either side fails
               if (.not. (viol <= 1.0_wp)) stats%nfail = stats%nfail + 1
               if (diff > stats%max_diff) stats%max_diff = diff
               if (viol > stats%max_viol) stats%max_viol = viol
            end do
         end do
      end do

   end subroutine compare_to_forward

   !> Run the routine under test and allocate its outputs
   !>
   !> @param[in]  cavity  Built cavity
   !> @param[in]  acc     Surface adjoints
   !> @param[out] w0      LSF value adjoints (ngrid)
   !> @param[out] w1      LSF gradient adjoints (3, ngrid)
   !> @param[out] w2      LSF Hessian adjoints (3, 3, ngrid)
   !> @param[out] error   Test failure
   subroutine lsf_weights(cavity, acc, w0, w1, w2, error)
      !> Built cavity
      type(cavity_type_drop), intent(in) :: cavity
      !> Surface adjoints
      type(cavity_adjoint_type), intent(in) :: acc
      !> LSF value adjoints
      real(wp), allocatable, intent(out) :: w0(:)
      !> LSF gradient adjoints
      real(wp), allocatable, intent(out) :: w1(:, :)
      !> LSF Hessian adjoints
      real(wp), allocatable, intent(out) :: w2(:, :, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: cav_error

      allocate (w0(cavity%ngrid), w1(3, cavity%ngrid), w2(3, 3, cavity%ngrid))
      call cavity%contract_surface_lsf_weights(acc, w0, w1, w2, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "contract_surface_lsf_weights failed: "// &
                          cav_error%message)
      end if

   end subroutine lsf_weights

   !> Initialize the accumulator and fill the chosen channels
   !>
   !> Every channel and component gets its own phase, so no two are proportional
   !>
   !> @param[out] acc       Surface adjoints
   !> @param[in]  ngrid     Grid size
   !> @param[in]  channels  Channel identifiers to populate
   subroutine fill_channels(acc, ngrid, channels)
      !> Surface adjoints
      type(cavity_adjoint_type), intent(out) :: acc
      !> Grid size
      integer, intent(in) :: ngrid
      !> Channel identifiers to populate
      integer, intent(in) :: channels(:)

      integer :: igrid, ax, ich

      call acc%init(ngrid)
      do ich = 1, size(channels)
         do igrid = 1, ngrid
            select case (channels(ich))
            case (CH_XYZ)
               do ax = 1, 3
                  acc%w_xyz(ax, igrid) = weight(igrid, ax, CH_XYZ)
               end do
            case (CH_XI)
               acc%w_xi(igrid) = weight(igrid, 0, CH_XI)
            case (CH_A)
               acc%w_a(igrid) = weight(igrid, 0, CH_A)
            case (CH_W)
               acc%w_w(igrid) = weight(igrid, 0, CH_W)
            case (CH_N)
               do ax = 1, 3
                  acc%w_n(ax, igrid) = weight(igrid, ax, CH_N)
               end do
            case (CH_K1)
               acc%w_k1(igrid) = weight(igrid, 0, CH_K1)
            case (CH_K2)
               acc%w_k2(igrid) = weight(igrid, 0, CH_K2)
            case (CH_F)
               acc%w_f(igrid) = weight(igrid, 0, CH_F)
            case default
               error stop "fill_channels: unknown channel"
            end select
         end do
      end do

   end subroutine fill_channels

   !> Deterministic pseudo-random weight in roughly [-1.3, 1.3]
   !>
   !> @param[in] igrid  Grid index
   !> @param[in] comp   Component, zero for scalar channels
   !> @param[in] ichan  Channel identifier
   !> @returns          Weight value
   pure function weight(igrid, comp, ichan) result(w)
      !> Grid index
      integer, intent(in) :: igrid
      !> Component, zero for scalar channels
      integer, intent(in) :: comp
      !> Channel identifier
      integer, intent(in) :: ichan
      !> Weight value
      real(wp) :: w

      w = sin(0.7_wp*real(igrid, wp) + 1.3_wp*real(comp, wp) + 2.1_wp*real(ichan, wp)) &
          + 0.3_wp*cos(0.21_wp*real(igrid, wp)*real(comp + 1, wp) &
                       + 0.9_wp*real(ichan, wp))

   end function weight

   !* ============================================================================== *!
   !*                                Density response                                 *!
   !* ============================================================================== *!

   !> A level set that is not a density leaves the response untouched
   !>
   !> @param[out] error  test failure
   subroutine test_density_response_svdw_untouched(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(response_type) :: empty, seeded
      type(density_response_type) :: seed, item
      type(mctc_error), allocatable :: cav_error
      integer :: nitem, ngrid

      call build_svdw_cavity(cavity, ctx, .false., error)
      if (allocated(error)) return
      ngrid = cavity%ngrid

      call fill_channels(acc, ngrid, [CH_XYZ, CH_XI, CH_N, CH_K1, CH_K2])
      ! Guard: the contraction would produce a nonzero item if it ran
      if (.not. (maxval(abs(acc%w_n)) > VACUITY_THR)) then
         call test_failed(error, "adjoint is vacuous")
         return
      end if

      call cavity%get_surface_response(acc, empty, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "get_surface_response failed: "//cav_error%message)
         return
      end if
      call count_items(empty, nitem, item)
      call check(error, nitem, 0, more="an empty response gained an item")
      if (allocated(error)) return

      allocate (seed%w_rho(ngrid), source=1.5_wp)
      allocate (seed%w_grad_rho(3, ngrid), source=-0.25_wp)
      allocate (seed%w_hess_rho(3, 3, ngrid), source=0.75_wp)
      call response_accumulate(seeded, seed, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "seeding the response failed: "//cav_error%message)
         return
      end if
      call cavity%get_surface_response(acc, seeded, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "get_surface_response failed: "//cav_error%message)
         return
      end if
      call count_items(seeded, nitem, item)
      call check(error, nitem, 1, more="the response item count changed")
      if (allocated(error)) return
      if (.not. allocated(item%w_rho)) then
         call test_failed(error, "the seeded density item vanished")
         return
      end if
      if (any(item%w_rho /= 1.5_wp) .or. any(item%w_grad_rho /= -0.25_wp) &
          .or. any(item%w_hess_rho /= 0.75_wp)) then
         call test_failed(error, "the seeded density item was modified")
      end if

   end subroutine test_density_response_svdw_untouched

   !> The density item is the level-set adjoint scaled by -scale, accumulated
   !>
   !> @param[out] error  test failure
   subroutine test_density_response_isodensity(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(response_type) :: response
      type(density_response_type) :: item
      type(mctc_error), allocatable :: cav_error
      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)
      real(wp) :: scale
      integer :: nitem, pass

      call build_iso_cavity(cavity, ctx, error)
      if (allocated(error)) return

      select type (model => cavity%lsf_model)
      type is (moist_cavity_drop_lsf_isodensity_internal_type)
         scale = model%param%scale
      class default
         call test_failed(error, "fixture level set is not the internal isodensity one")
         return
      end select
      ! A unit scale would not tell -scale from -1
      if (.not. (abs(scale - 1.0_wp) > 0.1_wp)) then
         call test_failed(error, "fixture scale is too close to one: "//to_string(scale))
         return
      end if

      call fill_channels(acc, cavity%ngrid, ALL_CHANNELS)
      call lsf_weights(cavity, acc, w0, w1, w2, error)
      if (allocated(error)) return
      if (.not. (maxval(abs(w0)) > VACUITY_THR .and. maxval(abs(w1)) > VACUITY_THR &
                 .and. maxval(abs(w2)) > VACUITY_THR)) then
         call test_failed(error, "level-set adjoint is vacuous: "// &
                          to_string(maxval(abs(w0)))//" "//to_string(maxval(abs(w1)))// &
                          " "//to_string(maxval(abs(w2))))
         return
      end if

      do pass = 1, 2
         call cavity%get_surface_response(acc, response, cav_error)
         if (allocated(cav_error)) then
            call test_failed(error, "get_surface_response failed: "//cav_error%message)
            return
         end if
         call count_items(response, nitem, item)
         call check(error, nitem, 1, &
                    more="expected one density item after call "//to_string(pass))
         if (allocated(error)) return
         call check_item(item, -scale*real(pass, wp), w0, w1, w2, &
                         "call "//to_string(pass), error)
         if (allocated(error)) return
      end do

   end subroutine test_density_response_isodensity

   !> Compare a density item with a multiple of the level-set adjoint
   !>
   !> @param[in]  item   Density item
   !> @param[in]  factor Expected multiple
   !> @param[in]  w0     LSF value adjoints
   !> @param[in]  w1     LSF gradient adjoints
   !> @param[in]  w2     LSF Hessian adjoints
   !> @param[in]  label  Description for failure messages
   !> @param[out] error  Test failure
   subroutine check_item(item, factor, w0, w1, w2, label, error)
      !> Density item
      type(density_response_type), intent(in) :: item
      !> Expected multiple
      real(wp), intent(in) :: factor
      !> LSF value adjoints
      real(wp), intent(in) :: w0(:)
      !> LSF gradient adjoints
      real(wp), intent(in) :: w1(:, :)
      !> LSF Hessian adjoints
      real(wp), intent(in) :: w2(:, :, :)
      !> Description for failure messages
      character(len=*), intent(in) :: label
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: dev0, dev1, dev2, size0, size1, size2

      if (.not. (allocated(item%w_rho) .and. allocated(item%w_grad_rho) &
                 .and. allocated(item%w_hess_rho))) then
         call test_failed(error, "density item is incomplete, "//label)
         return
      end if
      if (size(item%w_rho) /= size(w0) .or. size(item%w_grad_rho, 2) /= size(w1, 2) &
          .or. size(item%w_hess_rho, 3) /= size(w2, 3)) then
         call test_failed(error, "density item has the wrong grid size, "//label)
         return
      end if
      dev0 = maxval(abs(item%w_rho - factor*w0))
      dev1 = maxval(abs(item%w_grad_rho - factor*w1))
      dev2 = maxval(abs(item%w_hess_rho - factor*w2))
      size0 = max(1.0_wp, maxval(abs(factor*w0)))
      size1 = max(1.0_wp, maxval(abs(factor*w1)))
      size2 = max(1.0_wp, maxval(abs(factor*w2)))
      if (.not. (dev0 <= ITEM_TOL*size0 .and. dev1 <= ITEM_TOL*size1 &
                 .and. dev2 <= ITEM_TOL*size2)) then
         call test_failed(error, "density item differs from "//to_string(factor)// &
                          " * w_lsf, "//label//": deviations "//to_string(dev0)//" "// &
                          to_string(dev1)//" "//to_string(dev2))
      end if

   end subroutine check_item

   !> Count the items of a response and copy out its density item
   !>
   !> @param[in,out] response  Response to walk
   !> @param[out]   nitem     Number of items
   !> @param[out]   item      Last density item found, unallocated arrays otherwise
   subroutine count_items(response, nitem, item)
      !> Response to walk
      type(response_type), intent(inout) :: response
      !> Number of items
      integer, intent(out) :: nitem
      !> Last density item found
      type(density_response_type), intent(out) :: item

      class(response_item_type), allocatable :: current
      type(mctc_error), allocatable :: err

      nitem = 0
      do while (response%next())
         nitem = nitem + 1
         call current_response_item(response, current, err)
         if (allocated(err)) cycle
         select type (current)
         type is (density_response_type)
            item = current
         end select
      end do

   end subroutine count_items

   !* ============================================================================== *!
   !*                                    Refusals                                     *!
   !* ============================================================================== *!

   !> Report the expected refusal as a test failure
   !>
   !> The should_fail registration turns that failure into a pass; an absent
   !> error or a different message leaves the test passing, which is reported
   !>
   !> @param[out] error     Test failure, set only on the expected error
   !> @param[in]  err       Library error, possibly unallocated
   !> @param[in]  expected  Substring the expected message must contain
   subroutine expect_error(error, err, expected)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error, possibly unallocated
      type(mctc_error), allocatable, intent(in) :: err
      !> Substring the expected message must contain
      character(len=*), intent(in) :: expected

      if (.not. allocated(err)) return
      if (index(err%message, expected) > 0) call test_failed(error, err%message)

   end subroutine expect_error

   !> An accumulator that was never initialized is refused
   !>
   !> @param[out] error  test failure
   subroutine test_bad_accumulator_uninitialized(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)

      ! A failed setup returns cleanly, so the test is reported as passing
      call build_svdw_cavity(cavity, ctx, .false., error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      allocate (w0(cavity%ngrid), w1(3, cavity%ngrid), w2(3, 3, cavity%ngrid))
      call cavity%contract_surface_lsf_weights(acc, w0, w1, w2, err)
      call expect_error(error, err, "accumulator is not initialized")

   end subroutine test_bad_accumulator_uninitialized

   !> An accumulator of the wrong grid size is refused
   !>
   !> @param[out] error  test failure
   subroutine test_bad_accumulator_size(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)

      call build_svdw_cavity(cavity, ctx, .false., error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call acc%init(cavity%ngrid + 1)
      allocate (w0(cavity%ngrid), w1(3, cavity%ngrid), w2(3, 3, cavity%ngrid))
      call cavity%contract_surface_lsf_weights(acc, w0, w1, w2, err)
      call expect_error(error, err, "accumulator grid size mismatch")

   end subroutine test_bad_accumulator_size

   !> A response request on a cavity that was never built is refused
   !>
   !> @param[out] error  test failure
   subroutine test_bad_response_unbuilt(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(response_type) :: response
      type(mctc_error), allocatable :: err

      allocate (cavity)
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_svdw_cavity(cavity, ctx, .false., err)
      if (allocated(err)) return
      call acc%init(1)
      call cavity%get_surface_response(acc, response, err)
      call expect_error(error, err, "Cavity must be built first")

   end subroutine test_bad_response_unbuilt

   !> An unbuilt cavity is refused by the weight contraction, naming it
   !>
   !> @param[out] error  test failure
   subroutine test_bad_lsf_weights_unbuilt(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)

      allocate (cavity)
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_svdw_cavity(cavity, ctx, .false., err)
      if (allocated(err)) return
      call acc%init(1)
      allocate (w0(1), w1(3, 1), w2(3, 3, 1))
      call cavity%contract_surface_lsf_weights(acc, w0, w1, w2, err)
      call expect_error(error, err, "contract_surface_lsf_weights: Cavity must be built first")

   end subroutine test_bad_lsf_weights_unbuilt

   !> An unbuilt cavity is refused by the nuclear gradient
   !>
   !> - the base-type wrapper refuses before the DROP contraction is reached
   !>
   !> @param[out] error  test failure
   subroutine test_bad_gradient_unbuilt(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: gradient(:, :)

      allocate (cavity)
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_svdw_cavity(cavity, ctx, .false., err)
      if (allocated(err)) return
      call acc%init(1)
      allocate (gradient(3, max(1, cavity%nsph)), source=0.0_wp)
      call cavity%get_surface_gradient(acc, gradient, err)
      call expect_error(error, err, "Cavity must be built first")

   end subroutine test_bad_gradient_unbuilt

   !> Weight buffers that do not match the grid are refused before any write
   !>
   !> @param[out] error  test failure
   subroutine test_bad_lsf_weights_shape(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)

      ! A failed setup returns cleanly, so the test is reported as passing
      call build_svdw_cavity(cavity, ctx, .false., error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call acc%init(cavity%ngrid)
      ! Hessian weights in the packed (6, ngrid) layout instead of (3, 3, ngrid)
      allocate (w0(cavity%ngrid), w1(3, cavity%ngrid), w2(6, 1, cavity%ngrid))
      call cavity%contract_surface_lsf_weights(acc, w0, w1, w2, err)
      call expect_error(error, err, "contract_surface_lsf_weights: weight shape mismatch")

   end subroutine test_bad_lsf_weights_shape

   !> A surface point no atom reaches is refused by the response path, naming it
   !>
   !> @param[out] error  test failure
   subroutine test_bad_response_unreached(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)

      call build_svdw_cavity(cavity, ctx, .false., error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      if (cavity%ngrid < BAD_POINT) return
      cavity%xyz(:, BAD_POINT) = 100.0_wp
      call acc%init(cavity%ngrid)
      allocate (w0(cavity%ngrid), w1(3, cavity%ngrid), w2(3, 3, cavity%ngrid))
      call cavity%contract_surface_lsf_weights(acc, w0, w1, w2, err)
      call expect_error(error, err, "contract_surface_lsf_weights: level-set evaluation at "// &
                        "surface grid point "//to_string(BAD_POINT)//" is unreached")

   end subroutine test_bad_response_unreached

   !> A surface point no atom reaches is refused by the gradient path, naming it
   !>
   !> @param[out] error  test failure
   subroutine test_bad_gradient_unreached(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: gradient(:, :)

      call build_svdw_cavity(cavity, ctx, .false., error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      if (cavity%ngrid < BAD_POINT) return
      cavity%xyz(:, BAD_POINT) = 100.0_wp
      call acc%init(cavity%ngrid)
      allocate (gradient(3, cavity%nsph), source=0.0_wp)
      call cavity%get_surface_gradient(acc, gradient, err)
      call expect_error(error, err, "get_surface_gradient_drop: level-set evaluation at "// &
                        "surface grid point "//to_string(BAD_POINT)//" is unreached")

   end subroutine test_bad_gradient_unreached

   !> A host density failure is refused by the response path, naming its tile and column
   !>
   !> @param[out] error  test failure
   subroutine test_bad_response_host_failure(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      call run_host_fault(HOST_FAIL, .false., "contract_surface_lsf_weights: level-set "// &
                          "evaluation failed on the surface tile from grid point "// &
                          to_string(BAD_TILE_START)//": evaluate_batch failed at point "// &
                          to_string(BAD_POINT - BAD_TILE_START + 1)//":", error)

   end subroutine test_bad_response_host_failure

   !> A host density failure is refused by the gradient path, naming its tile and column
   !>
   !> @param[out] error  test failure
   subroutine test_bad_gradient_host_failure(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      call run_host_fault(HOST_FAIL, .true., "get_surface_gradient_drop: level-set "// &
                          "evaluation failed on the surface tile from grid point "// &
                          to_string(BAD_TILE_START)//": evaluate_batch failed at point "// &
                          to_string(BAD_POINT - BAD_TILE_START + 1)//":", error)

   end subroutine test_bad_gradient_host_failure

   !> A point with a vanishing level-set gradient is refused by the response path
   !>
   !> Its weights would otherwise stay zero
   !>
   !> @param[out] error  test failure
   subroutine test_bad_response_degenerate(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      call run_host_fault(HOST_FLAT, .false., "contract_surface_lsf_weights: level set "// &
                          "gradient vanishes at grid point "//to_string(BAD_POINT), error)

   end subroutine test_bad_response_degenerate

   !> A point with a vanishing level-set gradient is refused by the gradient path
   !>
   !> @param[out] error  test failure
   subroutine test_bad_gradient_degenerate(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      call run_host_fault(HOST_FLAT, .true., "get_surface_gradient_drop: level set "// &
                          "gradient vanishes at grid point "//to_string(BAD_POINT), error)

   end subroutine test_bad_gradient_degenerate

   !> Break the callback density at [[BAD_POINT]] after the build, then run one reverse path
   !>
   !> @param[in]  mode      Host failure mode, HOST_FAIL or HOST_FLAT
   !> @param[in]  gradient  Run the nuclear gradient instead of the density response
   !> @param[in]  expected  Substring the expected message must contain
   !> @param[out] error     Test failure, set only on the expected error
   subroutine run_host_fault(mode, gradient, expected, error)
      !> Host failure mode
      integer, intent(in) :: mode
      !> Run the nuclear gradient instead of the density response
      logical, intent(in) :: gradient
      !> Substring the expected message must contain
      character(len=*), intent(in) :: expected
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error

      type(gauss_host_type), target :: host
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
      type(response_type) :: response
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: grad(:, :)

      ! A failed setup returns cleanly, so the test is reported as passing
      call build_callback_cavity(cavity, ctx, host, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      if (cavity%ngrid < BAD_POINT) return
      host%bad_xyz = cavity%xyz(:, BAD_POINT)
      host%mode = mode
      call acc%init(cavity%ngrid)
      if (gradient) then
         allocate (grad(3, cavity%nsph), source=0.0_wp)
         call cavity%get_surface_gradient(acc, grad, err)
      else
         call cavity%get_surface_response(acc, response, err)
      end if
      call expect_error(error, err, expected)

   end subroutine run_host_fault

   !* ============================================================================== *!
   !*                                    Fixtures                                     *!
   !* ============================================================================== *!

   !> Construct the SvdW DROP cavity without building it
   !>
   !> @param[in,out] cavity   Allocated cavity
   !> @param[in,out] ctx      Run context borrowed by the cavity
   !> @param[in]    forward  Enable every optional property and forward array
   !> @param[out]   err      Library error
   subroutine new_svdw_cavity(cavity, ctx, forward, err)
      !> Allocated cavity
      type(cavity_type_drop), intent(inout) :: cavity
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(inout) :: ctx
      !> Enable every optional property and forward array
      logical, intent(in) :: forward
      !> Library error
      type(mctc_error), allocatable, intent(out) :: err

      type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
      type(moist_cavity_drop_parameters_type) :: param

      call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type( &
                             blend_k=BLEND_K, blend_3b=BLEND_3B))
      param = moist_cavity_drop_parameters_type( &
              num_leb=NUM_LEB, tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
              proj_level=PROJ_LEVEL, wleb_prune_level=4, do_normal=.true., &
              do_curvature=.true., do_fine=forward)
      call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), &
                           lsf_model=svdw_template, error=err, param=param, ctx=ctx)

   end subroutine new_svdw_cavity

   !> Build the SvdW DROP cavity on a small asymmetric molecule
   !>
   !> The two heavy centers make an elongated cavity with well-separated
   !> principal curvatures and the off-axis hydrogen removes the residual
   !> rotational symmetry, so no channel is accidentally degenerate
   !>
   !> @param[out]   cavity   Built cavity, forward arrays filled when requested
   !> @param[in,out] ctx      Run context borrowed by the cavity
   !> @param[in]    forward  Enable every optional property and forward array
   !> @param[out]   error    Test failure
   subroutine build_svdw_cavity(cavity, ctx, forward, error)
      !> Built cavity
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(inout) :: ctx
      !> Enable every optional property and forward array
      logical, intent(in) :: forward
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(mctc_error), allocatable :: cav_error

      call new(mol, [8, 6, 1], reshape([ &
                                       0.00_wp, 0.00_wp, 0.00_wp, &
                                       0.00_wp, 0.00_wp, 4.60_wp, &
                                       2.60_wp, 0.40_wp, -1.10_wp], [3, 3]))
      allocate (cavity)
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_svdw_cavity(cavity, ctx, forward, cav_error)
      if (.not. allocated(cav_error)) call cavity%update(mol, error=cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "failed to build cavity: "//cav_error%message)
         return
      end if
      if (cavity%ngrid == 0) then
         call test_failed(error, "empty grid")
         return
      end if
      if (.not. forward) return

      ! Forward sensitivity pass fills every *_rA array the reference reads
      call cavity%get_gradient(cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "forward gradient failed: "//cav_error%message)
         return
      end if
      if (.not. (allocated(cavity%xyz1_rA) .and. allocated(cavity%normal1_rA) &
                 .and. allocated(cavity%k1_rA) .and. allocated(cavity%k2_rA))) then
         call test_failed(error, "forward derivative arrays are incomplete")
      end if

   end subroutine build_svdw_cavity

   !> Build a DROP cavity on a two-center isodensity level set
   !>
   !> Each center carries one d shell weighting every cartesian component once,
   !> with different exponents so the pair has no mirror symmetry; the isovalue
   !> sets a level-set scale well away from one
   !>
   !> @param[out]   cavity  Built cavity
   !> @param[in,out] ctx     Run context borrowed by the cavity
   !> @param[out]   error   Test failure
   subroutine build_iso_cavity(cavity, ctx, error)
      !> Built cavity
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(inout) :: ctx
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(moist_cavity_drop_lsf_isodensity_param_type) :: param
      type(moist_cavity_drop_parameters_type) :: drop_param
      type(structure_type) :: mol
      class(radius_type), allocatable :: radii
      type(coupling_type), target :: coupling
      type(mctc_error), allocatable :: err
      real(wp) :: dcart(12, 12)
      integer :: i

      call new(mol, [54, 54], &
               reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 5.0_wp], [3, 2]))
      param = moist_cavity_drop_lsf_isodensity_param_type(core_fill=.false.)
      param%rho_iso = 0.4_wp
      param%scale = 1.0_wp/param%rho_iso
      call lsf%new(sh_atom=[1, 2], sh_l=[2, 2], sh_nprim=[1, 1], exps=[0.5_wp, 0.4_wp], &
                   coeffs=[1.0_wp, 1.0_wp], error=err, param=param)
      if (.not. allocated(err)) then
         dcart = 0.0_wp
         do i = 1, 12
            dcart(i, i) = 1.0_wp
         end do
         call lsf%set_density(dcart, err)
      end if
      if (.not. allocated(err)) call lsf%set_nuclear_charges([54.0_wp, 54.0_wp])
      if (.not. allocated(err)) call lsf%update(mol, [4.0_wp, 4.0_wp])
      if (.not. allocated(err)) call new_radii_custom_atoms([4.0_wp, 4.0_wp], radii, err)
      allocate (cavity)
      call new_context(ctx, nthreads=0, verbosity=0)
      if (.not. allocated(err)) then
         drop_param = moist_cavity_drop_parameters_type( &
                      num_leb=NUM_LEB, proj_level=PROJ_LEVEL, do_normal=.true., &
                      do_curvature=.true.)
         call new_cavity_drop(cavity, radius_model=radii, lsf_model=lsf, error=err, &
                              param=drop_param, ctx=ctx)
      end if
      if (.not. allocated(err)) call cavity%update(mol, err)
      if (.not. allocated(err)) call cavity%prepare_build(coupling, err)
      if (.not. allocated(err)) then
         do while (coupling%next())
            select type (request => coupling%request())
            type is (density_matrix_request_type)
               call coupling%answer("dcart", dcart, err)
            type is (nuclear_charge_request_type)
               call coupling%answer("z", [54.0_wp, 54.0_wp], err)
            end select
            if (allocated(err)) exit
         end do
      end if
      if (.not. allocated(err)) call cavity%build(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "failed to build isodensity cavity: "//err%message)
         return
      end if
      if (cavity%ngrid == 0) call test_failed(error, "empty isodensity grid")

   end subroutine build_iso_cavity

   !> Build a DROP cavity on a host callback density of two Gaussian centers
   !>
   !> The host stays healthy through the build; the caller breaks it afterwards
   !>
   !> @param[out]   cavity  Built cavity
   !> @param[in,out] ctx     Run context borrowed by the cavity
   !> @param[in,out] host    Host density; the cavity keeps its address
   !> @param[out]   error   Test failure
   subroutine build_callback_cavity(cavity, ctx, host, error)
      !> Built cavity
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(inout) :: ctx
      !> Host density; must outlive the cavity
      type(gauss_host_type), target, intent(inout) :: host
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf
      type(moist_cavity_drop_parameters_type) :: drop_param
      type(structure_type) :: mol
      class(radius_type), allocatable :: radii
      type(coupling_type), target :: coupling
      type(mctc_error), allocatable :: err

      host%center = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 5.0_wp], [3, 2])
      host%expo = [0.5_wp, 0.4_wp]
      host%mode = HOST_OK
      call new(mol, [8, 8], host%center)
      call lsf%new(callback_ptr=c_funloc(gauss_host_callback), context=c_loc(host), &
                   param=moist_cavity_drop_lsf_isodensity_param_type(core_fill=.false.))
      call new_radii_custom_atoms([4.0_wp, 4.0_wp], radii, err)
      allocate (cavity)
      call new_context(ctx, nthreads=0, verbosity=0)
      if (.not. allocated(err)) then
         drop_param = moist_cavity_drop_parameters_type(num_leb=NUM_LEB, proj_level=PROJ_LEVEL)
         call new_cavity_drop(cavity, radius_model=radii, lsf_model=lsf, error=err, &
                              param=drop_param, ctx=ctx)
      end if
      if (.not. allocated(err)) call cavity%update(mol, err)
      if (.not. allocated(err)) call cavity%prepare_build(coupling, err)
      if (.not. allocated(err)) then
         ! The core fill is off, so the build phase asks for nothing
         do while (coupling%next())
         end do
         call cavity%build(coupling, err)
      end if
      if (allocated(err)) then
         call test_failed(error, "failed to build callback cavity: "//err%message)
         return
      end if
      if (cavity%ngrid == 0) call test_failed(error, "empty callback grid")

   end subroutine build_callback_cavity

   !> Host density callback: a sum of Gaussians, broken on demand at one point
   !>
   !> @param[in]  context Address of a [[gauss_host_type]]
   !> @param[in]  point   Evaluation point in Bohr
   !> @param[out] rho     Density
   !> @param[out] drho    Density gradient
   !> @param[out] d2rho   Density Hessian (Fortran (3,3), or NULL)
   !> @param[out] d3rho   Density third derivative (Fortran (3,3,3), or NULL)
   !> @returns            0, or 7 for HOST_FAIL at the broken point
   function gauss_host_callback(context, point, rho, drho, d2rho, d3rho) result(status) bind(C)
      !> Address of the host state
      type(c_ptr), value :: context
      !> Evaluation point in Bohr
      real(c_double), intent(in) :: point(3)
      !> Density
      real(c_double), intent(out) :: rho
      !> Density gradient
      real(c_double), intent(out) :: drho(3)
      !> Density Hessian, or NULL
      type(c_ptr), value :: d2rho
      !> Density third derivative, or NULL
      type(c_ptr), value :: d3rho
      !> Zero on success
      integer(c_int) :: status

      type(gauss_host_type), pointer :: host
      real(c_double), pointer :: hptr(:, :), tptr(:, :, :)
      real(wp) :: d(3), g, a, rho_w, drho_w(3), h(3, 3), t(3, 3, 3), delta(3, 3)
      logical :: broken
      integer :: ic, i, j, l

      status = 0_c_int
      call c_f_pointer(context, host)
      broken = host%mode /= HOST_OK .and. all(real(point, wp) == host%bad_xyz)
      if (broken .and. host%mode == HOST_FAIL) then
         rho = 0.0_c_double
         drho = 0.0_c_double
         status = 7_c_int
         return
      end if

      delta = 0.0_wp
      do i = 1, 3
         delta(i, i) = 1.0_wp
      end do
      rho_w = 0.0_wp
      drho_w = 0.0_wp
      h = 0.0_wp
      t = 0.0_wp
      do ic = 1, size(host%expo)
         a = host%expo(ic)
         d = real(point, wp) - host%center(:, ic)
         g = exp(-a*sum(d*d))
         rho_w = rho_w + g
         drho_w = drho_w - 2.0_wp*a*d*g
         do j = 1, 3
            do i = 1, 3
               h(i, j) = h(i, j) + (4.0_wp*a*a*d(i)*d(j) - 2.0_wp*a*delta(i, j))*g
               do l = 1, 3
                  t(i, j, l) = t(i, j, l) + (-8.0_wp*a*a*a*d(i)*d(j)*d(l) &
                                             + 4.0_wp*a*a*(delta(i, j)*d(l) + delta(i, l)*d(j) &
                                                           + delta(j, l)*d(i)))*g
               end do
            end do
         end do
      end do
      ! A flat density keeps its value but loses every derivative
      if (broken) then
         drho_w = 0.0_wp
         h = 0.0_wp
         t = 0.0_wp
      end if

      rho = real(rho_w, c_double)
      drho = real(drho_w, c_double)
      if (c_associated(d2rho)) then
         call c_f_pointer(d2rho, hptr, [3, 3])
         hptr = real(h, c_double)
      end if
      if (c_associated(d3rho)) then
         call c_f_pointer(d3rho, tptr, [3, 3, 3])
         tptr = real(t, c_double)
      end if
   end function gauss_host_callback

end module test_cavity_drop_potential
