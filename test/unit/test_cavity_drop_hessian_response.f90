!> Adjoint-response half of the DROP surface Hessian
!>
!> [[get_surface_hessian_response_drop]] is the `J^T (d omega/dv)` term: the
!> primal map is held fixed and the *folded* surface adjoints move, because
!> [[prepare_surface_weights]] builds them out of `a`, `wleb`, `xi0` and the
!> branch softmax. Driving `acc%w_a` and `acc%w_w` is what makes them move; an
!> accumulator without those two channels has a geometry-independent `eff` and
!> this half is then identically zero.
!>
!> ## Ground truth
!>
!> Write the shipped gradient as `G(R; acc) = Phi(R) . eff(R, acc)`, linear in
!> `eff`. Central-differencing `G` along a nuclear direction gives
!> `Phi' . eff + Phi . d(eff)`, the sum of *both* halves, so a reference for
!> this half alone needs the first term removed. It is removed by subtraction,
!> against a second accumulator whose folded weights are frozen:
!>
!>     acc_frozen:  w_xi := eff(acc)%w_xi |_base,  w_f := eff(acc)%w_f |_base,
!>                  w_a  := 0,                     w_w := 0
!>
!> and the remaining channels copied unchanged. [[prepare_surface_weights]]
!> then folds nothing, so `eff(R, acc_frozen)` reproduces `eff(acc)|_base`
!> exactly at every geometry -- except through `branch_phi_adj`, which it
!> re-derives from the moving grid. Hence the identity this suite asserts,
!> which holds on a branched grid as well as on an unbranched one:
!>
!>     response(acc) - response(acc_frozen)
!>         ==  d/dv [ G(R; acc) - G(R; acc_frozen) ]
!>
!> Both `Phi'` terms cancel exactly, because the two accumulators share their
!> base-geometry `eff`. On an unbranched grid `response(acc_frozen)` is
!> identically zero and the left-hand side collapses to this half on its own;
!> `frozen_response_is_zero` asserts that separately, so the subtraction cannot
!> be hiding a term.
!>
!> The subtraction is taken *inside* the stencil, one geometry at a time:
!> `G(R; acc) - G(R; acc_frozen)` vanishes at the base geometry and grows
!> linearly with the displacement, so everything the two gradients share is
!> gone before the division by `h` rather than after it.
!>
!> ## Which channels the subtraction tests drive, and why not all of them
!>
!> `w_xyz`, `w_n`, `w_k1` and `w_k2` are `source=`-copies of the raw adjoints
!> (`weights.f90`), so their tangent is identically zero and they **cannot
!> contribute to this half at all**. In the subtraction reference they are
!> therefore pure ballast: an identical term in both gradients, cancelling
!> analytically and costing precision numerically. The three subtraction cases
!> drive `w_xi`, `w_f`, `w_a` and `w_w` -- every channel that can move -- and
!> `svdw_plain_all_channels_fd` exists to *measure* what the four frozen ones
!> cost when they are added back. They cost three orders:
!>
!> | channels driven                | svdw/plain, h = 2.0e-4    |
!> |--------------------------------|---------------------------|
!> | w_xi, w_f, w_a, w_w            | 4.0e-11 abs / 0    rel    |
!> | + w_xyz, w_n                   | 3.1e-10 abs / 7.5e-11 rel |
!> | + w_k1, w_k2                   | 1.0e-08 abs / 8.3e-09 rel |
!>
!> The curvature row is the near-umbilic amplification of `hessian.md`
!> (`kernel.f90:489`) arriving by a different route than in
!> `test_cavity_drop_hessian_fixed`: `res%dk1` is not merely noisy on this
!> fixture, it is *large*, so `w_k1 res%dk1` dominates both gradients and the
!> live-frozen difference loses the digits it dominates. Nothing in this half
!> reads a curvature response -- `deff%have_wk` is `.false.` by construction --
!> so that row is a statement about the reference, not about the code under
!> test. It is left in the suite as a failing case with a measured number
!> rather than removed or given a tolerance of its own.
!>
!> ## The composite, and why it needs the surrogate accumulator
!>
!> `both_halves_svdw` / `both_halves_cfc` are the first assertions in the
!> project that run the two halves of the Hessian together:
!>
!>     d/dv [ G(R; acc) ]  ==  H_fixed . v  +  response(acc, v)
!>
!> `H_fixed` cannot be asked for with `acc` itself: [[check_frozen_weights]]
!> refuses any accumulator with a live `w_a` or `w_w`, which is exactly the
!> accumulator this half exists for. It is asked for with `acc_frozen`
!> instead, and that is not a workaround but an identity --
!> [[get_surface_hessian_fixed_drop]] reads its accumulator only through
!> [[prepare_surface_weights]], and the two produce the same `eff` at the base
!> geometry by construction. There is no subtraction in the composite's
!> reference, so `w_xyz` and `w_n` are driven here as well; the curvature
!> channels are not, because the fixed half is itself noise limited on them
!> (see its own suite's header).
!>
!> ## Fixtures
!>
!>   * `FIX_PLAIN`, the asymmetric OCH triple, never branches. `wbranch` is
!>     exactly one everywhere, `branch_phi_adj` and its tangent are identically
!>     zero, and the composite above is available. SvdW and CFC.
!>   * `FIX_CROSS`, the five-carbon cross at `proj_level = 7` with a softened
!>     softmax (`s = 0.5`), does branch. That is the only fixture on which
!>     `dbranch_phi_adj` is nonzero, and therefore the only one that can catch a
!>     driver which drops the branch channel of pass 2. The fixed half refuses
!>     this grid, so the composite is not available here.
!>
!> **The branched fixture does not reach the project bound, and the reason is
!> the fixture.** It plateaus at `1.5e-8` absolute across the whole window and
!> rises as `1/h` below it, which puts the round-off of the gradient difference
!> at `~7e-12` -- three orders above the same quantity on `FIX_PLAIN`. The
!> amplifier is the branches themselves: `get_test_cross` has no near-degenerate
!> minima, only 3.5-4.4 Bohr secondary solutions kept alive by the softened
!> softmax, and their `wleb` sits near the pruning floor with `xi0` three orders
!> above the grid's typical value. [[branch_point_adjoint]] forms
!> `-0.5 w_xi xi0 / wbranch` on exactly those points, so `branch_phi_adj` is
!> enormous there and cancels in the group mean. A better-conditioned branched
!> fixture would need branches of comparable weight, which this geometry does
!> not have; the case therefore ships as a failure with a measured number.
!>
!> ## The mutation that shaped the branched case
!>
!> `deff%branch_phi_adj` was zeroed after pass 2 in
!> `derivatives/hessian_response.f90` and the suite re-run:
!>
!>   * `svdw_cross_branching_fd` fails at `6.50e-3` absolute and `4.67e-3`
!>     relative -- worst component atom 2 axis 3, direction 1, analytic
!>     `12.5652314` against a reference of `12.5717282`. That is five and a half
!>     orders above the fixture's own `1.5e-8` floor;
!>   * the number is *identical in every digit at all of h = 2.5e-4, 2.0e-4 and
!>     1.5e-4*. A step-independent deviation is a missing term, not a stencil
!>     artefact, and that is the cleanest part of the signature;
!>   * every `FIX_PLAIN` case is unchanged, bit for bit, including the two
!>     composites and the all-channel case. `branch_phi_adj` is structurally
!>     zero there, so a driver that drops it is invisible to an unbranched
!>     fixture -- which is why one branched fixture is worth its cost even at
!>     `1.5e-8`.
!>
!> ## Step and grid guard
!>
!> `FD_STEPS` is measured; the sweep is in the comment on that parameter. The
!> grid is guarded at every stencil geometry ([[assert_grid_match]]) on
!> `numbering`, `owner`, `branch_count` and `anchor_id`, so a step large enough
!> to re-enumerate the grid fails loudly instead of putting a step into the
!> reference.
module test_cavity_drop_hessian_response
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_radii, only: default_cpcm_radii
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: fill_legacy_radii, get_test_cross

   implicit none(type, external)
   private

   public :: collect_cavity_drop_hessian_response

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Level-set model of a fixture
   integer, parameter :: LSF_SVDW = 1, LSF_CFC = 2

   !> Geometry of a fixture
   integer, parameter :: FIX_PLAIN = 1, FIX_CROSS = 2

   !> Directions pushed through in one call: one sparse, one dense
   integer, parameter :: NDIR = 2

   !> Surface-adjoint channels. `w_a` and `w_w` are the two whose fold moves
   !> with the geometry, and are what this half exists to differentiate
   integer, parameter :: CH_XI = 1, CH_F = 2, CH_XYZ = 3, CH_N = 4
   integer, parameter :: CH_K1 = 5, CH_K2 = 6, CH_A = 7, CH_W = 8
   integer, parameter :: NCHAN = 8

   !> Shared fixture settings
   real(wp), parameter :: BLEND_K = 2.5_wp
   real(wp), parameter :: BLEND_3B = 1.0_wp
   integer, parameter :: NUM_LEB = 50
   real(wp), parameter :: PROJ_TOL = 1.0E-14_wp
   integer, parameter :: PROJ_MAXITER = 1000
   integer, parameter :: PROJ_LEVEL = 2
   integer, parameter :: WLEB_PRUNE = 4

   !> Branching fixture: multistart level and softmax temperature
   real(wp), parameter :: CROSS_BLEND_K = 1.0_wp
   integer, parameter :: CROSS_PROJ_LEVEL = 7
   real(wp), parameter :: CROSS_BRANCH_S = 0.5_wp

   !> Production softmax temperature, used by the non-branching fixture
   real(wp), parameter :: PLAIN_BRANCH_S = 0.0025_wp

   !> Central-difference steps of the subtraction identity
   !>
   !> Measured over both directions, as `worst absolute / worst relative`, the
   !> relative one taken over exactly the components that already exceed
   !> `FD_ABS` -- that is, over the components that decide the test, so a `0`
   !> entry means no component is over the absolute bound at all:
   !>
   !> | h      | svdw/plain        | cfc/plain         | svdw/cross        |
   !> |--------|-------------------|-------------------|-------------------|
   !> | 1.0e-3 | 1.10e-08 / 2.5e-09| 2.34e-07 / 1.1e-08| 8.92e-09 / 6.1e-07|
   !> | 6.0e-4 | 1.43e-09 / 3.2e-10| 3.03e-08 / 1.0e-09| 1.09e-08 / 7.3e-08|
   !> | 5.0e-4 | 6.99e-10 / 1.6e-10| 1.46e-08 / 4.9e-10| 1.84e-08 / 4.1e-08|
   !> | 4.0e-4 | 3.01e-10 / 6.2e-11| 5.96e-09 / 2.0e-10| 1.73e-08 / 3.4e-08|
   !> | 3.5e-4 | 1.90e-10 / 4.1e-11| 3.52e-09 / 1.2e-10| 1.91e-08 / 3.0e-08|
   !> | 3.0e-4 | 8.69e-11 / 0      | 1.89e-09 / 6.4e-11| 1.66e-08 / 3.6e-08|
   !> | 2.5e-4 | 7.32e-11 / 0      | 8.75e-10 / 3.0e-11| 1.51e-08 / 6.8e-08|
   !> | 2.0e-4 | 4.04e-11 / 0      | 3.87e-10 / 1.3e-11| 1.51e-08 / 8.9e-08|
   !> | 1.5e-4 | 8.68e-11 / 0      | 1.31e-10 / 2.2e-10| 6.20e-08 / 8.2e-08|
   !> | 1.0e-4 | 1.62e-10 / 1.9e-09| 8.98e-11 / 0      | 3.84e-08 / 3.1e-08|
   !>
   !> A clean bowl on `FIX_PLAIN`: `O(h^4)` truncation on the way down -- the
   !> `cfc/plain` absolute column falls by a factor of 124 from `1e-3` to
   !> `3e-4`, which is `(10/3)^4` -- and round-off on the way up below `1.5e-4`.
   !> The two steps below are inside the flat bottom for both level-set models
   !> and are a factor of 1.25 apart. Both are required, so a value that agrees
   !> at one step only -- the signature of a step sitting on the round-off wall
   !> -- fails. `1.5e-4` is deliberately not one of them: `cfc/plain` has a
   !> round-off spike there (`2.2e-10` relative) that `2.0e-4` and `1.0e-4` on
   !> either side of it do not.
   !>
   !> The floor is two orders below what the sibling suites reach, and the
   !> subtraction is why: `G(acc) - G(acc_frozen)` is formed at each stencil
   !> geometry, so everything the two gradients share -- the whole `Phi' . eff`
   !> term, most of the projection's own round-off -- never reaches the
   !> difference quotient. `svdw/cross` is the case that does not benefit; see
   !> the module header for what its floor is made of.
   real(wp), parameter :: FD_STEPS(2) = [2.5E-4_wp, 2.0E-4_wp]

   !> Central-difference steps of the composite assertion
   !>
   !> The composite differences the shipped gradient itself rather than a
   !> difference of two gradients, so it inherits the fixed half's floor rather
   !> than this suite's. Measured worst deviation of
   !> `FD - (H_fixed . v + response)`, over both directions:
   !>
   !> | h      | svdw/plain        | cfc/plain         |
   !> |--------|-------------------|-------------------|
   !> | 3.0e-3 | 3.07e-05 / 2.4e-06| 2.16e-05 / 9.8e-07|
   !> | 1.0e-3 | 3.79e-07 / 3.0e-08| 2.66e-07 / 1.2e-08|
   !> | 6.0e-4 | 4.91e-08 / 3.9e-09| 3.45e-08 / 1.6e-09|
   !> | 4.0e-4 | 9.67e-09 / 7.8e-10| 6.84e-09 / 3.1e-10|
   !> | 3.0e-4 | 3.10e-09 / 2.6e-10| 2.22e-09 / 9.8e-11|
   !> | 2.5e-4 | 1.50e-09 / 1.3e-10| 1.02e-09 / 4.6e-11|
   !> | 2.2e-4 | 9.53e-10 / 2.0e-10| 6.45e-10 / 2.1e-11|
   !> | 1.9e-4 | 5.06e-10 / 4.2e-11| 2.75e-10 / 4.1e-11|
   !> | 1.7e-4 | 3.35e-10 / 1.3e-09| 2.49e-10 / 4.5e-10|
   !> | 1.5e-4 | 2.60e-10 / 6.9e-11| 2.28e-10 / 3.3e-11|
   !> | 1.4e-4 | 1.40e-10 / 1.4e-10| 1.29e-10 / 3.2e-10|
   !> | 1.2e-4 | 2.35e-10 / 2.5e-10| 2.55e-10 / 4.7e-10|
   !> | 1.0e-4 | 3.64e-10 / 2.6e-09| 2.16e-10 / 4.7e-10|
   !>
   !> Same bowl, one order shallower and with a much rougher bottom: the
   !> absolute column never goes under `1e-10`, so the composite is decided by
   !> the relative bound throughout, and the relative column below `2e-4` is
   !> round-off jitter rather than a trend -- `1.7e-4` misses by an order while
   !> `1.9e-4` and `1.5e-4` on either side of it pass. Those two are the pair,
   !> and they are the only pair in the sweep that passes for both level-set
   !> models. The margin is a factor of two, and that is a statement about the
   !> reference, not about either half.
   real(wp), parameter :: COMPOSITE_STEPS(2) = [1.9E-4_wp, 1.5E-4_wp]

   !> Finite-difference agreement bounds
   !>
   !> The project target, `1e-10` absolute and `1e-10` relative; a component
   !> fails only when it misses *both*. Met by `svdw_plain_fd`, `cfc_plain_fd`
   !> and both composites. Missed, deliberately and with the numbers in the
   !> module header, by `svdw_plain_all_channels_fd` (the curvature channels'
   !> known amplification, `1.0e-8`) and by `svdw_cross_branching_fd` (the
   !> branched fixture's own floor, `1.5e-8`). Neither is given a tolerance of
   !> its own.
   real(wp), parameter :: FD_ABS = 1.0E-10_wp
   real(wp), parameter :: FD_REL = 1.0E-10_wp

   !> Below this the differenced reference carries no information and a
   !> comparison against it would pass for the wrong reason
   real(wp), parameter :: VACUITY_THR = 1.0E-6_wp

contains

   !> Collect the suite
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_cavity_drop_hessian_response(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("svdw_plain_fd", test_svdw_plain), &
                  new_unittest("cfc_plain_fd", test_cfc_plain), &
                  new_unittest("svdw_plain_all_channels_fd", test_svdw_plain_all), &
                  new_unittest("svdw_cross_branching_fd", test_svdw_cross), &
                  new_unittest("both_halves_svdw", test_both_halves_svdw), &
                  new_unittest("both_halves_cfc", test_both_halves_cfc), &
                  new_unittest("frozen_response_is_zero", test_frozen_response_is_zero), &
                  new_unittest("shape_guards", test_shape_guards) &
                  ]
   end subroutine collect_cavity_drop_hessian_response

   !* ================================================================================= *!
   !*                       Subtraction identity: this half alone                       *!
   !* ================================================================================= *!

   !> SvdW on the non-branching fixture
   !>
   !> @param[out] error Error handle
   subroutine test_svdw_plain(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_response_fd(FIX_PLAIN, LSF_SVDW, moving_channels(), "svdw/plain", error)
   end subroutine test_svdw_plain

   !> CFC on the non-branching fixture
   !>
   !> @param[out] error Error handle
   subroutine test_cfc_plain(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_response_fd(FIX_PLAIN, LSF_CFC, moving_channels(), "cfc/plain", error)
   end subroutine test_cfc_plain

   !> Every channel at once, the four with a zero tangent included
   !>
   !> **Expected to fail at the project bound, at `1.0e-8`.** The four channels
   !> it adds cannot contribute to this half -- their tangent is identically
   !> zero -- so they only add a common term to both gradients of the reference,
   !> and the curvature pair's near-umbilic magnitude is what the subtraction
   !> then loses its digits to. Kept because the number is worth recording and
   !> would otherwise be invisible; see the module header for the
   !> channel-by-channel table.
   !>
   !> @param[out] error Error handle
   subroutine test_svdw_plain_all(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_response_fd(FIX_PLAIN, LSF_SVDW, all_channels(), "svdw/plain+all", error)
   end subroutine test_svdw_plain_all

   !> SvdW on the branching fixture
   !>
   !> The only case in which `dbranch_phi_adj` is nonzero, and **expected to
   !> fail at the project bound, at `1.5e-8`**: the fixture's far-branch points
   !> put a `7e-12` round-off floor under the gradient difference. The module
   !> header has the diagnosis, and the mutation that shows the branch channel
   !> is nevertheless present and right to five and a half orders better than a
   !> driver which drops it.
   !>
   !> @param[out] error Error handle
   subroutine test_svdw_cross(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_response_fd(FIX_CROSS, LSF_SVDW, moving_channels(), "svdw/cross", error)
   end subroutine test_svdw_cross

   !> Central-difference the frozen-weight subtraction against this half
   !>
   !> @param[in]  fix_kind Geometry of the fixture
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  channels Adjoint channels to populate
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine run_response_fd(fix_kind, lsf_kind, channels, label, error)
      !> Geometry of the fixture
      integer, intent(in) :: fix_kind
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Adjoint channels to populate
      integer, intent(in) :: channels(:)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol

      !> Base-geometry folded weights, which define the frozen accumulator
      real(wp), allocatable :: eff_xi(:), eff_f(:)
      !> Analytic response of both accumulators, and their difference
      real(wp), allocatable :: hvp(:, :, :), hvp_frozen(:, :, :), analytic(:, :, :)
      !> Differenced reference, one block per step
      real(wp), allocatable :: fd(:, :, :, :)
      !> Nuclear directions
      real(wp), allocatable :: dirs(:, :, :)
      !> Extents and loop indices
      integer :: nsph, idir, istep

      call fixture_geometry(fix_kind, mol)
      call build_cavity(cavity, ctx, mol, fix_kind, lsf_kind, error)
      if (allocated(error)) return
      call assert_branching(cavity, fix_kind, "base geometry", error)
      if (allocated(error)) return

      nsph = cavity%nsph
      call build_directions(nsph, dirs)
      call fold_effective(cavity, channels, eff_xi, eff_f)

      allocate (hvp(ndim, nsph, NDIR), source=0.0_wp)
      allocate (hvp_frozen(ndim, nsph, NDIR), source=0.0_wp)
      call response_of(cavity, channels, dirs, .false., eff_xi, eff_f, hvp, label, error)
      if (allocated(error)) return
      call response_of(cavity, channels, dirs, .true., eff_xi, eff_f, hvp_frozen, label, error)
      if (allocated(error)) return

      allocate (analytic(ndim, nsph, NDIR))
      analytic = hvp - hvp_frozen

      !* ------------------------- Central-difference reference ------------------------ *!
      allocate (fd(ndim, nsph, NDIR, size(FD_STEPS)), source=0.0_wp)
      do istep = 1, size(FD_STEPS)
         do idir = 1, NDIR
            call fd_frozen_difference(mol, cavity, fix_kind, lsf_kind, channels, &
                                      eff_xi, eff_f, dirs(:, :, idir), FD_STEPS(istep), &
                                      fd(:, :, idir, istep), label, error)
            if (allocated(error)) return
         end do
      end do

      call compare_to_reference(analytic, fd, FD_STEPS, label, error)

   end subroutine run_response_fd

   !* ================================================================================= *!
   !*                        Composite: both halves of the Hessian                      *!
   !* ================================================================================= *!

   !> SvdW: the fixed half plus this one reproduce the differenced gradient
   !>
   !> @param[out] error Error handle
   subroutine test_both_halves_svdw(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_composite_fd(LSF_SVDW, "svdw/plain", error)
   end subroutine test_both_halves_svdw

   !> CFC: the fixed half plus this one reproduce the differenced gradient
   !>
   !> @param[out] error Error handle
   subroutine test_both_halves_cfc(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_composite_fd(LSF_CFC, "cfc/plain", error)
   end subroutine test_both_halves_cfc

   !> Central-difference the shipped gradient against the sum of both halves
   !>
   !> Only on `FIX_PLAIN`: the fixed half refuses a branched grid, and on a
   !> branched grid it would also be missing the second-order branch term.
   !>
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine run_composite_fd(lsf_kind, label, error)
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_surface_adjoint_type) :: acc_frozen
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol

      !> Base-geometry folded weights, which define the frozen accumulator
      real(wp), allocatable :: eff_xi(:), eff_f(:)
      !> Fixed half, this half and their sum
      real(wp), allocatable :: hess(:, :, :, :), hvp(:, :, :), total(:, :, :)
      !> Differenced reference, one block per step
      real(wp), allocatable :: fd(:, :, :, :)
      !> Nuclear directions
      real(wp), allocatable :: dirs(:, :, :)
      !> Extents and loop indices
      integer :: nsph, iatom, iaxis, idir, istep
      !> Channels driven by the composite
      integer :: channels(6)

      channels = moving_channels()

      call fixture_geometry(FIX_PLAIN, mol)
      call build_cavity(cavity, ctx, mol, FIX_PLAIN, lsf_kind, error)
      if (allocated(error)) return
      call assert_branching(cavity, FIX_PLAIN, "base geometry", error)
      if (allocated(error)) return

      nsph = cavity%nsph
      call build_directions(nsph, dirs)
      call fold_effective(cavity, channels, eff_xi, eff_f)

      !* ------------------------------- The fixed half -------------------------------- *!
      ! Asked for with the frozen surrogate, because `check_frozen_weights`
      ! refuses a live `w_a` or `w_w`. The two accumulators fold to the same
      ! `eff` at this geometry, and the fixed half reads nothing else.
      call seed_adjoint(cavity, channels, .true., eff_xi, eff_f, acc_frozen, error)
      if (allocated(error)) return

      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      call cavity%get_surface_hessian_fixed(acc_frozen, hess, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "fixed-adjoint Hessian failed ("//label//"): "// &
                          cav_error%message)
         return
      end if

      !* ------------------------------ The response half ------------------------------ *!
      allocate (hvp(ndim, nsph, NDIR), source=0.0_wp)
      call response_of(cavity, channels, dirs, .false., eff_xi, eff_f, hvp, label, error)
      if (allocated(error)) return

      allocate (total(ndim, nsph, NDIR), source=0.0_wp)
      do idir = 1, NDIR
         total(:, :, idir) = hvp(:, :, idir)
         do iatom = 1, nsph
            do iaxis = 1, ndim
               total(:, :, idir) = total(:, :, idir) &
                                   + hess(:, :, iaxis, iatom)*dirs(iaxis, iatom, idir)
            end do
         end do
      end do

      ! Both halves have to carry something, or the sum could be right for the
      ! wrong reason
      if (maxval(abs(hvp)) <= VACUITY_THR) then
         call test_failed(error, "the response half is vacuous for "//label//" (max "// &
                          to_string(maxval(abs(hvp)))//")")
         return
      end if
      if (maxval(abs(hess)) <= VACUITY_THR) then
         call test_failed(error, "the fixed half is vacuous for "//label//" (max "// &
                          to_string(maxval(abs(hess)))//")")
         return
      end if

      !* ------------------------- Central-difference reference ------------------------ *!
      allocate (fd(ndim, nsph, NDIR, size(COMPOSITE_STEPS)), source=0.0_wp)
      do istep = 1, size(COMPOSITE_STEPS)
         do idir = 1, NDIR
            call fd_surface_gradient(mol, cavity, FIX_PLAIN, lsf_kind, channels, &
                                     eff_xi, eff_f, dirs(:, :, idir), &
                                     COMPOSITE_STEPS(istep), fd(:, :, idir, istep), &
                                     label, error)
            if (allocated(error)) return
         end do
      end do

      call compare_to_reference(total, fd, COMPOSITE_STEPS, "both halves, "//label, error)

   end subroutine run_composite_fd

   !* ================================================================================= *!
   !*                              Structural properties                                *!
   !* ================================================================================= *!

   !> A frozen accumulator has no response on an unbranched grid
   !>
   !> With `w_a` and `w_w` zero the folds do nothing, and without a branched
   !> anchor `branch_phi_adj` is identically zero, so every channel pass 2 emits
   !> is exactly zero and so is the contraction. This is an identity rather than
   !> a tolerance, and it is what licenses reading the subtraction identity of
   !> `svdw_plain_fd` as a statement about this half alone.
   !>
   !> @param[out] error Error handle
   subroutine test_frozen_response_is_zero(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol

      real(wp), allocatable :: eff_xi(:), eff_f(:)
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :)
      integer :: channels(NCHAN)

      channels = all_channels()

      call fixture_geometry(FIX_PLAIN, mol)
      call build_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error)
      if (allocated(error)) return
      call assert_branching(cavity, FIX_PLAIN, "single-branch fixture", error)
      if (allocated(error)) return

      call build_directions(cavity%nsph, dirs)
      call fold_effective(cavity, channels, eff_xi, eff_f)

      allocate (hvp(ndim, cavity%nsph, NDIR), source=0.0_wp)
      call response_of(cavity, channels, dirs, .true., eff_xi, eff_f, hvp, "frozen", error)
      if (allocated(error)) return

      if (maxval(abs(hvp)) /= 0.0_wp) then
         call test_failed(error, "a frozen accumulator produced a response on an"// &
                          " unbranched grid (max "//to_string(maxval(abs(hvp)))//")")
         return
      end if

      ! The live accumulator must not be zero as well, or the assertion above
      ! would be satisfied by a routine that returns nothing at all
      hvp = 0.0_wp
      call response_of(cavity, channels, dirs, .false., eff_xi, eff_f, hvp, "live", error)
      if (allocated(error)) return
      if (maxval(abs(hvp)) <= VACUITY_THR) then
         call test_failed(error, "the live accumulator has no response either (max "// &
                          to_string(maxval(abs(hvp)))//")")
         return
      end if

   end subroutine test_frozen_response_is_zero

   !> Mis-shaped arguments must be refused, and must not write anything
   !>
   !> @param[out] error Error handle
   subroutine test_shape_guards(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol

      real(wp), allocatable :: eff_xi(:), eff_f(:)
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :)
      integer :: nsph

      call fixture_geometry(FIX_PLAIN, mol)
      call build_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error)
      if (allocated(error)) return

      nsph = cavity%nsph
      call fold_effective(cavity, moving_channels(), eff_xi, eff_f)
      call seed_adjoint(cavity, moving_channels(), .false., eff_xi, eff_f, acc, error)
      if (allocated(error)) return

      ! Wrong number of spheres in `dirs`
      allocate (dirs(ndim, nsph + 1, NDIR), source=0.1_wp)
      allocate (hvp(ndim, nsph, NDIR), source=0.0_wp)
      call cavity%get_surface_hessian_response(acc, dirs, hvp, cav_error)
      if (.not. allocated(cav_error)) then
         call test_failed(error, "a mis-shaped direction array was accepted")
         return
      end if
      deallocate (cav_error, dirs)

      ! Accumulator with the wrong direction count
      allocate (dirs(ndim, nsph, NDIR), source=0.1_wp)
      deallocate (hvp)
      allocate (hvp(ndim, nsph, NDIR + 1), source=0.0_wp)
      call cavity%get_surface_hessian_response(acc, dirs, hvp, cav_error)
      if (.not. allocated(cav_error)) then
         call test_failed(error, "a mis-shaped accumulator was accepted")
         return
      end if
      deallocate (cav_error)

      ! The accumulator is `intent(inout)`: a refused call owes the caller the
      ! block it was handed, untouched
      if (maxval(abs(hvp)) /= 0.0_wp) then
         call test_failed(error, "a refused call wrote into the accumulator")
         return
      end if

   end subroutine test_shape_guards

   !* ================================================================================= *!
   !*                         Finite-difference references                              *!
   !* ================================================================================= *!

   !> Five-point central difference of `G(R; acc) - G(R; acc_frozen)`
   !>
   !> `O(h^4)`. A three-point stencil leaves a truncation error that, at a step
   !> small enough to keep the grid stable, still sits orders above the
   !> round-off floor of the subtraction. The subtraction itself is done at each
   !> stencil geometry rather than between two finished derivatives, so the two
   !> gradients' common content -- including the whole fixed half -- never
   !> reaches the difference quotient.
   !>
   !> @param[in]  mol      Base structure
   !> @param[in]  ref_cav  Base cavity, for the grid comparison
   !> @param[in]  fix_kind Geometry of the fixture
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  channels Adjoint channels to populate
   !> @param[in]  eff_xi   Base-geometry folded width adjoint
   !> @param[in]  eff_f    Base-geometry folded switching adjoint
   !> @param[in]  vdir     Nuclear direction `(3, nsph)`
   !> @param[in]  step     Central-difference step
   !> @param[out] deriv    Differenced difference `(3, nsph)`
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine fd_frozen_difference(mol, ref_cav, fix_kind, lsf_kind, channels, &
                                   eff_xi, eff_f, vdir, step, deriv, label, error)
      !> Base structure
      type(structure_type), intent(in) :: mol
      !> Base cavity
      type(cavity_type_drop), intent(in) :: ref_cav
      !> Geometry of the fixture
      integer, intent(in) :: fix_kind
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Adjoint channels to populate
      integer, intent(in) :: channels(:)
      !> Base-geometry folded weights
      real(wp), intent(in) :: eff_xi(:), eff_f(:)
      !> Nuclear direction
      real(wp), intent(in) :: vdir(:, :)
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Differenced difference
      real(wp), intent(out) :: deriv(:, :)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Five-point central stencil of the first derivative
      integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
      real(wp), parameter :: COEFF(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol_disp
      real(wp), allocatable :: grad_live(:, :), grad_frozen(:, :)
      integer :: iside
      character(len=32) :: side

      deriv = 0.0_wp
      allocate (grad_live(size(deriv, 1), size(deriv, 2)))
      allocate (grad_frozen(size(deriv, 1), size(deriv, 2)))

      do iside = 1, size(OFFSET)
         write (side, "(a, i0, a, es9.2, a)") "offset ", OFFSET(iside), " (h = ", step, ")"

         mol_disp = mol
         mol_disp%xyz = mol%xyz + real(OFFSET(iside), wp)*step*vdir

         call build_cavity(cavity, ctx, mol_disp, fix_kind, lsf_kind, error)
         if (allocated(error)) return
         call assert_grid_match(ref_cav, cavity, label//" "//trim(side), error)
         if (allocated(error)) return

         call surface_gradient(cavity, channels, .false., eff_xi, eff_f, grad_live, &
                               label, error)
         if (allocated(error)) return
         call surface_gradient(cavity, channels, .true., eff_xi, eff_f, grad_frozen, &
                               label, error)
         if (allocated(error)) return

         deriv = deriv + COEFF(iside)*(grad_live - grad_frozen)/step

         deallocate (cavity)
      end do

   end subroutine fd_frozen_difference

   !> Five-point central difference of the shipped surface gradient
   !>
   !> @param[in]  mol      Base structure
   !> @param[in]  ref_cav  Base cavity, for the grid comparison
   !> @param[in]  fix_kind Geometry of the fixture
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  channels Adjoint channels to populate
   !> @param[in]  eff_xi   Base-geometry folded width adjoint
   !> @param[in]  eff_f    Base-geometry folded switching adjoint
   !> @param[in]  vdir     Nuclear direction `(3, nsph)`
   !> @param[in]  step     Central-difference step
   !> @param[out] deriv    Differenced gradient `(3, nsph)`
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine fd_surface_gradient(mol, ref_cav, fix_kind, lsf_kind, channels, &
                                  eff_xi, eff_f, vdir, step, deriv, label, error)
      !> Base structure
      type(structure_type), intent(in) :: mol
      !> Base cavity
      type(cavity_type_drop), intent(in) :: ref_cav
      !> Geometry of the fixture
      integer, intent(in) :: fix_kind
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Adjoint channels to populate
      integer, intent(in) :: channels(:)
      !> Base-geometry folded weights
      real(wp), intent(in) :: eff_xi(:), eff_f(:)
      !> Nuclear direction
      real(wp), intent(in) :: vdir(:, :)
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Differenced gradient
      real(wp), intent(out) :: deriv(:, :)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Five-point central stencil of the first derivative
      integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
      real(wp), parameter :: COEFF(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol_disp
      real(wp), allocatable :: grad(:, :)
      integer :: iside
      character(len=32) :: side

      deriv = 0.0_wp
      allocate (grad(size(deriv, 1), size(deriv, 2)))

      do iside = 1, size(OFFSET)
         write (side, "(a, i0, a, es9.2, a)") "offset ", OFFSET(iside), " (h = ", step, ")"

         mol_disp = mol
         mol_disp%xyz = mol%xyz + real(OFFSET(iside), wp)*step*vdir

         call build_cavity(cavity, ctx, mol_disp, fix_kind, lsf_kind, error)
         if (allocated(error)) return
         call assert_grid_match(ref_cav, cavity, label//" "//trim(side), error)
         if (allocated(error)) return

         call surface_gradient(cavity, channels, .false., eff_xi, eff_f, grad, label, error)
         if (allocated(error)) return

         deriv = deriv + COEFF(iside)*grad/step

         deallocate (cavity)
      end do

   end subroutine fd_surface_gradient

   !> Compare an analytic block against the differenced reference of every step
   !>
   !> The whole block is scanned before anything is reported: the deviation this
   !> suite exists to expose is the worst one, and failing on the first
   !> component over the bound would name an arbitrary early one instead.
   !>
   !> @param[in]  analytic Analytic block `(3, nsph, NDIR)`
   !> @param[in]  fd       Differenced reference `(3, nsph, NDIR, nstep)`
   !> @param[in]  steps    Steps behind `fd`
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine compare_to_reference(analytic, fd, steps, label, error)
      !> Analytic block
      real(wp), intent(in) :: analytic(:, :, :)
      !> Differenced reference
      real(wp), intent(in) :: fd(:, :, :, :)
      !> Steps behind the reference
      real(wp), intent(in) :: steps(:)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer :: istep, idir, iatom, iaxis
      integer :: bad_step, bad_dir, bad_atom, bad_axis
      real(wp) :: ref, diff, worst, worst_rel, bad_ana, bad_ref

      ! Anti-vacuity, per direction: a reference at machine zero would be
      ! matched by anything, and a fixture in which one direction went quiet
      ! would still pass on the strength of the other
      do idir = 1, size(fd, 3)
         if (maxval(abs(fd(:, :, idir, 1))) <= VACUITY_THR) then
            call test_failed(error, "the differenced reference is vacuous for "//label// &
                             ", direction "//to_string(idir)//" (max "// &
                             to_string(maxval(abs(fd(:, :, idir, 1))))//")")
            return
         end if
      end do

      worst = 0.0_wp
      worst_rel = 0.0_wp
      bad_step = 0
      bad_dir = 0
      bad_atom = 0
      bad_axis = 0
      bad_ana = 0.0_wp
      bad_ref = 0.0_wp
      do istep = 1, size(steps)
         do idir = 1, size(fd, 3)
            do iatom = 1, size(fd, 2)
               do iaxis = 1, size(fd, 1)
                  ref = fd(iaxis, iatom, idir, istep)
                  diff = abs(analytic(iaxis, iatom, idir) - ref)
                  ! A component has to miss both bounds to be a failure: the
                  ! absolute one alone would condemn a large component, the
                  ! relative one alone a component that is numerically zero
                  if (diff > FD_ABS .and. diff > FD_REL*abs(ref)) then
                     if (diff > worst) then
                        worst = diff
                        worst_rel = diff/max(abs(ref), tiny(1.0_wp))
                        bad_step = istep
                        bad_dir = idir
                        bad_atom = iatom
                        bad_axis = iaxis
                        bad_ana = analytic(iaxis, iatom, idir)
                        bad_ref = ref
                     end if
                  end if
               end do
            end do
         end do
      end do

      if (bad_step > 0) then
         call test_failed(error, "adjoint-response Hessian mismatch for "//label// &
                          ": worst deviation "//to_string(worst)//" absolute, "// &
                          to_string(worst_rel)//" relative, at atom "//to_string(bad_atom)// &
                          " axis "//to_string(bad_axis)//", direction "//to_string(bad_dir)// &
                          " (h = "//to_string(steps(bad_step))//"): analytic "// &
                          to_string(bad_ana)//" finite difference "//to_string(bad_ref))
         return
      end if

   end subroutine compare_to_reference

   !* ================================================================================= *!
   !*                          Calls into the routines under test                       *!
   !* ================================================================================= *!

   !> Run the adjoint-response half for one accumulator
   !>
   !> @param[in]  cavity   Cavity to differentiate
   !> @param[in]  channels Adjoint channels to populate
   !> @param[in]  dirs     Nuclear directions `(3, nsph, NDIR)`
   !> @param[in]  frozen   Whether to freeze the folded weights
   !> @param[in]  eff_xi   Base-geometry folded width adjoint
   !> @param[in]  eff_f    Base-geometry folded switching adjoint
   !> @param[out] hvp      Response block `(3, nsph, NDIR)`
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine response_of(cavity, channels, dirs, frozen, eff_xi, eff_f, hvp, label, error)
      !> Cavity to differentiate
      type(cavity_type_drop), intent(in) :: cavity
      !> Adjoint channels to populate
      integer, intent(in) :: channels(:)
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Whether to freeze the folded weights
      logical, intent(in) :: frozen
      !> Base-geometry folded weights
      real(wp), intent(in) :: eff_xi(:), eff_f(:)
      !> Response block
      real(wp), intent(out) :: hvp(:, :, :)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error

      hvp = 0.0_wp
      call seed_adjoint(cavity, channels, frozen, eff_xi, eff_f, acc, error)
      if (allocated(error)) return

      call cavity%get_surface_hessian_response(acc, dirs, hvp, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "adjoint-response Hessian failed ("//label//"): "// &
                          cav_error%message)
         return
      end if

   end subroutine response_of

   !> Run the shipped surface gradient for one accumulator
   !>
   !> @param[in]  cavity   Cavity to contract
   !> @param[in]  channels Adjoint channels to populate
   !> @param[in]  frozen   Whether to freeze the folded weights
   !> @param[in]  eff_xi   Base-geometry folded width adjoint
   !> @param[in]  eff_f    Base-geometry folded switching adjoint
   !> @param[out] grad     Nuclear gradient `(3, nsph)`
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine surface_gradient(cavity, channels, frozen, eff_xi, eff_f, grad, label, error)
      !> Cavity to contract
      type(cavity_type_drop), intent(in) :: cavity
      !> Adjoint channels to populate
      integer, intent(in) :: channels(:)
      !> Whether to freeze the folded weights
      logical, intent(in) :: frozen
      !> Base-geometry folded weights
      real(wp), intent(in) :: eff_xi(:), eff_f(:)
      !> Nuclear gradient
      real(wp), intent(out) :: grad(:, :)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error

      grad = 0.0_wp
      call seed_adjoint(cavity, channels, frozen, eff_xi, eff_f, acc, error)
      if (allocated(error)) return

      call cavity%get_surface_gradient(acc, grad, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "surface gradient failed ("//label//"): "// &
                          cav_error%message)
         return
      end if

   end subroutine surface_gradient

   !* ================================================================================= *!
   !*                             Preconditions of the tests                            *!
   !* ================================================================================= *!

   !> Assert that the fixture branches exactly as much as it is meant to
   !>
   !> @param[in]  cavity   Cavity to inspect
   !> @param[in]  fix_kind Geometry of the fixture
   !> @param[in]  label    Human-readable geometry description
   !> @param[out] error    Error handle
   subroutine assert_branching(cavity, fix_kind, label, error)
      !> Cavity to inspect
      type(cavity_type_drop), intent(in) :: cavity
      !> Geometry of the fixture
      integer, intent(in) :: fix_kind
      !> Geometry description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      logical :: branched

      if (.not. allocated(cavity%branch_count)) then
         call test_failed(error, "no branch bookkeeping at "//label)
         return
      end if
      branched = any(cavity%branch_count(1:cavity%ngrid) > 1)

      select case (fix_kind)
      case (FIX_PLAIN)
         if (branched) then
            call test_failed(error, "the single-branch fixture branched at "//label// &
                             " (max branch_count "// &
                             to_string(maxval(cavity%branch_count(1:cavity%ngrid)))//")")
            return
         end if
         if (any(cavity%wbranch(1:cavity%ngrid) /= 1.0_wp)) then
            call test_failed(error, "wbranch is not exactly one on the single-branch"// &
                             " fixture at "//label)
            return
         end if
      case default
         if (.not. branched) then
            call test_failed(error, "the branching fixture did not branch at "//label// &
                             "; dbranch_phi_adj would then be structurally zero and the"// &
                             " branch channel of pass 2 untested")
            return
         end if
      end select

   end subroutine assert_branching

   !> Assert that two geometries carry the very same grid points
   !>
   !> The frozen weights are carried across the stencil slot by slot, so a point
   !> that appears, vanishes or changes branch identity between two geometries
   !> would silently put a step into the reference. `branch_count` and
   !> `anchor_id` are compared as well as `numbering` and `owner`, because a
   !> group that re-enumerates without changing the point set would corrupt the
   !> branch stage alone. This assertion is what bounds the step from above.
   !>
   !> @param[in]  ref   Reference cavity
   !> @param[in]  cav   Displaced cavity
   !> @param[in]  label Human-readable geometry description
   !> @param[out] error Error handle
   subroutine assert_grid_match(ref, cav, label, error)
      !> Reference cavity
      type(cavity_type_drop), intent(in) :: ref
      !> Displaced cavity
      type(cavity_type_drop), intent(in) :: cav
      !> Geometry description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer :: igrid

      if (cav%ngrid /= ref%ngrid) then
         call test_failed(error, "the grid changed size at "//label//" ("// &
                          to_string(ref%ngrid)//" -> "//to_string(cav%ngrid)// &
                          "); the step is above the branch-enumeration ceiling")
         return
      end if

      do igrid = 1, ref%ngrid
         if (cav%numbering(igrid) /= ref%numbering(igrid)) then
            call test_failed(error, "the grid was reordered or repopulated at "//label// &
                             ", slot "//to_string(igrid)//" ("// &
                             to_string(ref%numbering(igrid))//" -> "// &
                             to_string(cav%numbering(igrid))//")")
            return
         end if
         if (cav%owner(igrid) /= ref%owner(igrid)) then
            call test_failed(error, "the owner sphere changed at "//label//", slot "// &
                             to_string(igrid)//" ("//to_string(ref%owner(igrid))//" -> "// &
                             to_string(cav%owner(igrid))//")")
            return
         end if
         if (cav%branch_count(igrid) /= ref%branch_count(igrid)) then
            call test_failed(error, "the branch count changed at "//label//", slot "// &
                             to_string(igrid)//" ("//to_string(ref%branch_count(igrid))// &
                             " -> "//to_string(cav%branch_count(igrid))//")")
            return
         end if
         if (cav%anchor_id(igrid) /= ref%anchor_id(igrid)) then
            call test_failed(error, "the anchor group changed at "//label//", slot "// &
                             to_string(igrid)//" ("//to_string(ref%anchor_id(igrid))// &
                             " -> "//to_string(cav%anchor_id(igrid))//")")
            return
         end if
      end do

   end subroutine assert_grid_match

   !* ================================================================================= *!
   !*                                   Accumulators                                    *!
   !* ================================================================================= *!

   !> Channels the subtraction identity and the composite drive
   !>
   !> `w_a` and `w_w` are the two that make `eff` geometry dependent, so they
   !> are the point of the fixture rather than an addition to it, and `w_xi` and
   !> `w_f` are what they fold into. `w_xyz` and `w_n` have an identically zero
   !> tangent and ride along for realism: they cost the subtraction an order of
   !> precision and are still inside the bound, while the curvature pair costs
   !> two more and is not -- see the channel table in the module header.
   !>
   !> @return Channel identifiers
   pure function moving_channels() result(channels)
      !> Channel identifiers
      integer :: channels(6)

      channels = [CH_XI, CH_F, CH_XYZ, CH_N, CH_A, CH_W]
   end function moving_channels

   !> Every channel the accumulator carries
   !>
   !> @return Channel identifiers
   pure function all_channels() result(channels)
      !> Channel identifiers
      integer :: channels(NCHAN)

      channels = [CH_XI, CH_F, CH_XYZ, CH_N, CH_K1, CH_K2, CH_A, CH_W]
   end function all_channels

   !> Fold the base-geometry effective weights the frozen accumulator restores
   !>
   !> Mirrors [[prepare_surface_weights]] with `fold_switching = .true.`; the
   !> two folds are the only geometry-dependent part of it, and the guards it
   !> takes are unconditional here because [[point_weight]] is bounded well away
   !> from `seed_weight_tol = 1e-30`.
   !>
   !> This is a deliberate duplication of three lines of production code: `eff`
   !> is not reachable from a test, and reproducing it is what lets the frozen
   !> accumulator exist at all. It does not have to be bit-for-bit -- a one-ulp
   !> drift enters the identity as `Phi' . (E' - E)`, some `1e-14` of the
   !> reference.
   !>
   !> @param[in]  cavity   Cavity supplying the grid
   !> @param[in]  channels Adjoint channels the accumulator carries
   !> @param[out] eff_xi   Folded width adjoint (ngrid)
   !> @param[out] eff_f    Folded switching adjoint (ngrid)
   subroutine fold_effective(cavity, channels, eff_xi, eff_f)
      !> Cavity supplying the grid
      type(cavity_type_drop), intent(in) :: cavity
      !> Adjoint channels the accumulator carries
      integer, intent(in) :: channels(:)
      !> Folded weights
      real(wp), allocatable, intent(out) :: eff_xi(:), eff_f(:)

      real(wp) :: w_a, w_w, r_own
      integer :: ngrid, igrid

      ngrid = cavity%ngrid
      allocate (eff_xi(ngrid), eff_f(ngrid))

      do igrid = 1, ngrid
         eff_xi(igrid) = raw_weight(cavity%numbering(igrid), CH_XI, channels)
         eff_f(igrid) = raw_weight(cavity%numbering(igrid), CH_F, channels)

         w_a = raw_weight(cavity%numbering(igrid), CH_A, channels)
         w_w = raw_weight(cavity%numbering(igrid), CH_W, channels)
         r_own = cavity%radii(cavity%owner(igrid))

         eff_xi(igrid) = eff_xi(igrid) &
                         - 2.0_wp*cavity%a(igrid)*w_a/cavity%xi0(igrid) &
                         - 2.0_wp*cavity%wleb(igrid)*w_w/cavity%xi0(igrid)
         eff_f(igrid) = eff_f(igrid) + w_a*r_own*r_own*cavity%wleb(igrid)
      end do

   end subroutine fold_effective

   !> Populate a surface-adjoint accumulator, live or frozen
   !>
   !> The live accumulator is a pure function of the persistent point id
   !> `cavity%numbering`, so the same weights are reproduced on a displaced grid
   !> without an explicit mapping. The frozen one replaces the width and
   !> switching channels by the base-geometry *folded* weights and zeroes `w_a`
   !> and `w_w`, so [[prepare_surface_weights]] folds nothing and returns those
   !> same weights at every geometry.
   !>
   !> @param[in]  cavity   Cavity supplying the grid
   !> @param[in]  channels Channel identifiers to populate
   !> @param[in]  frozen   Whether to build the frozen surrogate
   !> @param[in]  eff_xi   Base-geometry folded width adjoint
   !> @param[in]  eff_f    Base-geometry folded switching adjoint
   !> @param[out] acc      Surface-adjoint accumulator
   !> @param[out] error    Error handle
   subroutine seed_adjoint(cavity, channels, frozen, eff_xi, eff_f, acc, error)
      !> Cavity supplying the grid
      type(cavity_type_drop), intent(in) :: cavity
      !> Channel identifiers to populate
      integer, intent(in) :: channels(:)
      !> Whether to build the frozen surrogate
      logical, intent(in) :: frozen
      !> Base-geometry folded weights
      real(wp), intent(in) :: eff_xi(:), eff_f(:)
      !> Surface-adjoint accumulator
      type(cavity_surface_adjoint_type), intent(out) :: acc
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: add_error
      real(wp), allocatable :: ws(:), wv(:, :)
      integer :: ngrid, igrid, iaxis, ichannel, ich

      ngrid = cavity%ngrid
      allocate (ws(ngrid), wv(ndim, ngrid))
      call acc%init(ngrid)

      do ichannel = 1, size(channels)
         ich = channels(ichannel)

         ! The frozen surrogate carries the folded width and switching channels
         ! and no area or weight channel at all
         if (frozen .and. (ich == CH_A .or. ich == CH_W)) cycle
         if (frozen .and. ich == CH_XI) then
            ws = eff_xi
         else if (frozen .and. ich == CH_F) then
            ws = eff_f
         else
            do igrid = 1, ngrid
               ws(igrid) = point_weight(cavity%numbering(igrid), ich)
            end do
         end if

         do igrid = 1, ngrid
            do iaxis = 1, ndim
               wv(iaxis, igrid) = point_weight(cavity%numbering(igrid), NCHAN*ich + iaxis)
            end do
         end do

         select case (ich)
         case (CH_XI)
            call acc%add_surface_weights(add_error, w_xi=ws)
         case (CH_F)
            call acc%add_surface_weights(add_error, w_f=ws)
         case (CH_XYZ)
            call acc%add_surface_weights(add_error, w_xyz=wv)
         case (CH_N)
            call acc%add_surface_weights(add_error, w_n=wv)
         case (CH_K1)
            call acc%add_surface_weights(add_error, w_k1=ws)
         case (CH_K2)
            call acc%add_surface_weights(add_error, w_k2=ws)
         case (CH_A)
            call acc%add_surface_weights(add_error, w_a=ws)
         case default
            call acc%add_surface_weights(add_error, w_w=ws)
         end select
         if (allocated(add_error)) then
            call test_failed(error, "failed to seed the surface adjoint: "// &
                             add_error%message)
            return
         end if
      end do

   end subroutine seed_adjoint

   !> Raw adjoint weight of one channel, zero when the fixture does not drive it
   !>
   !> @param[in] id       Persistent point id
   !> @param[in] channel  Channel selector
   !> @param[in] channels Channels the fixture drives
   !> @return             Raw adjoint weight
   pure function raw_weight(id, channel, channels) result(w)
      !> Persistent point id
      integer, intent(in) :: id
      !> Channel selector
      integer, intent(in) :: channel
      !> Channels the fixture drives
      integer, intent(in) :: channels(:)
      !> Raw adjoint weight
      real(wp) :: w

      w = 0.0_wp
      if (any(channels == channel)) w = point_weight(id, channel)
   end function raw_weight

   !> Reproducible adjoint weight of one persistent point and channel
   !>
   !> Varies across the grid so that a term which happens to cancel for uniform
   !> weights still shows up, and is bounded well away from zero so that no
   !> channel is accidentally switched off.
   !>
   !> @param[in] id      Persistent point id, `cavity%numbering`
   !> @param[in] channel Channel selector
   !> @return            Adjoint weight
   pure function point_weight(id, channel) result(w)
      !> Persistent point id
      integer, intent(in) :: id
      !> Channel selector
      integer, intent(in) :: channel
      !> Adjoint weight
      real(wp) :: w

      w = 0.60_wp + 0.35_wp*sin(0.7_wp*real(id, wp) + 1.3_wp*real(channel, wp)) &
          + 0.11_wp*cos(0.23_wp*real(id, wp)*real(channel + 2, wp))
   end function point_weight

   !* ================================================================================= *!
   !*                                     Fixture                                       *!
   !* ================================================================================= *!

   !> The two nuclear directions pushed through in one call
   !>
   !> Direction 1 moves a single atom along a single axis -- the sparsest
   !> column, and the one a wrong influence set would zero out. Direction 2
   !> moves every atom along its own vector and is not a rigid translation, so
   !> nothing about it cancels. Both go through in one call, so the per-point
   !> direction loop is exercised rather than a degenerate `ndir = 1` path.
   !>
   !> @param[in]  nsph Number of spheres
   !> @param[out] dirs Nuclear directions `(3, nsph, NDIR)`
   subroutine build_directions(nsph, dirs)
      !> Number of spheres
      integer, intent(in) :: nsph
      !> Nuclear directions
      real(wp), allocatable, intent(out) :: dirs(:, :, :)

      integer :: iatom, iaxis

      if (allocated(dirs)) deallocate (dirs)
      allocate (dirs(ndim, nsph, NDIR), source=0.0_wp)

      dirs(3, 2, 1) = 1.0_wp

      do iatom = 1, nsph
         do iaxis = 1, ndim
            dirs(iaxis, iatom, 2) = sin(1.1_wp*real(iatom, wp) + 0.6_wp*real(iaxis, wp))
         end do
      end do

   end subroutine build_directions

   !> Fixture geometry
   !>
   !> `FIX_PLAIN` is asymmetric on purpose: a symmetric geometry drives the
   !> multistart projection into sibling branches, which is the other fixture's
   !> job. `FIX_CROSS` is the shared five-carbon cross, whose concave seams give
   !> several minima per anchor.
   !>
   !> @param[in]  fix_kind Geometry selector
   !> @param[out] mol      Structure
   subroutine fixture_geometry(fix_kind, mol)
      !> Geometry selector
      integer, intent(in) :: fix_kind
      !> Structure
      type(structure_type), intent(out) :: mol

      select case (fix_kind)
      case (FIX_PLAIN)
         call new(mol, [8, 6, 1], reshape([ &
                                          0.00_wp, 0.00_wp, 0.00_wp, &
                                          0.00_wp, 0.00_wp, 4.60_wp, &
                                          2.60_wp, 0.40_wp, -1.10_wp], [3, 3]))
      case default
         call get_test_cross(mol)
      end select

   end subroutine fixture_geometry

   !> Build the DROP cavity for a fixture and a level-set model
   !>
   !> @param[out]   cavity   Constructed cavity
   !> @param[inout] ctx      Run context borrowed by the cavity; must outlive it
   !> @param[in]    mol      Structure to build on
   !> @param[in]    fix_kind Geometry of the fixture
   !> @param[in]    lsf_kind Level-set model
   !> @param[out]   error    Error handle
   subroutine build_cavity(cavity, ctx, mol, fix_kind, lsf_kind, error)
      !> Constructed cavity
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(inout) :: ctx
      !> Structure to build on
      type(structure_type), intent(in) :: mol
      !> Geometry of the fixture
      integer, intent(in) :: fix_kind
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: radii(:)
      type(mctc_error), allocatable :: cav_error
      !> Per-fixture settings. Named apart from the module parameters they are
      !> assigned from: Fortran is case insensitive, so a local `blend_k` would
      !> shadow `BLEND_K` and turn the assignment into a silent self-assignment.
      real(wp) :: blend_k_loc, branch_s_loc
      integer :: proj_level_loc

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      ! The softmax scale is only raised on the branching fixture; it also sets
      ! the admissible branch radius, and that is what keeps the cross's far
      ! siblings alive -- and what makes them the noisiest points on the grid.
      if (fix_kind == FIX_CROSS) then
         blend_k_loc = CROSS_BLEND_K
         proj_level_loc = CROSS_PROJ_LEVEL
         branch_s_loc = CROSS_BRANCH_S
      else
         blend_k_loc = BLEND_K
         proj_level_loc = PROJ_LEVEL
         branch_s_loc = PLAIN_BRANCH_S
      end if

      allocate (cavity)
      call new_context(ctx, verbosity=0)
      select case (lsf_kind)
      case (LSF_SVDW)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=blend_k_loc, blend_3b=BLEND_3B)
            call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                                 tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                                 proj_level=proj_level_loc, wleb_prune_level=WLEB_PRUNE, &
                                 branch_weight_s=branch_s_loc, &
                                 radius_model=default_cpcm_radii(), &
                                 lsf_model=svdw_template, error=cav_error)
         end block
      case default
         block
            type(moist_cavity_drop_lsf_cfc_type) :: cfc_template
            call cfc_template%new()
            call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                                 tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                                 proj_level=proj_level_loc, wleb_prune_level=WLEB_PRUNE, &
                                 branch_weight_s=branch_s_loc, &
                                 radius_model=default_cpcm_radii(), &
                                 lsf_model=cfc_template, error=cav_error)
         end block
      end select
      if (allocated(cav_error)) then
         call test_failed(error, "failed to initialize cavity: "//cav_error%message)
         return
      end if

      ! Curvature and normals are surface observables this suite drives, so the
      ! cavity has to be asked for them
      call cavity%properties(do_fine=.true.)

      call cavity%update(mol, error=cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "failed to build cavity: "//cav_error%message)
         return
      end if
      if (cavity%ngrid == 0) then
         call test_failed(error, "empty grid")
         return
      end if

   end subroutine build_cavity

end module test_cavity_drop_hessian_response
