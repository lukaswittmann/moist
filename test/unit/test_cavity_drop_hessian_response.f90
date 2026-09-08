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
!> The curvature row is the near-umbilic amplification of the `kernel.f90`
!> discriminant: the principal-curvature gap used to be formed as
!> `sqrt(max(KM^2 - KG, 0))`, a difference of two `KM^2`-sized quantities that
!> `apply_seed` then divides by (fixed 2026-09-07, so this table is the pre-fix
!> measurement). It arrives here by a different route than in
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
!>     softmax (`s = 2.0`), does branch. That is the only fixture on which
!>     `dbranch_phi_adj` is nonzero, and therefore the only one that can catch a
!>     driver which drops the branch channel of pass 2. The fixed half refuses
!>     this grid, so the composite is not available here.
!>
!> ## What the softmax temperature of the branched fixture buys
!>
!> `CROSS_BRANCH_S` was `0.5` until 2026-09-07, and at that value
!> `svdw_cross_branching_fd` could not be made to pass at any step: it bottomed
!> out at `1.5e-8`, a factor of 150 over the bound. The cause was not the
!> derivative and not the stencil. It was `max|G|`.
!>
!> **The floor of a differenced reference is `~4 eps max|G| / h`**, where `max|G|`
!> is the largest component of the gradient being differenced -- not of the
!> derivative being checked. Measured on this fixture across six softmax
!> temperatures, `floor * h / (eps max|G|)` stays in `3.5 .. 7.0` while `max|G|`
!> itself moves by a factor of 134, so the law fixes the whole design:
!>
!> | `s`  | branched pts | `max|branch_phi_adj|` | `max|G|` | `floor * h` |
!> |------|--------------|-----------------------|----------|-------------|
!> | 0.25 |      0       |  0                    |   35.2   |   4e-14     |
!> | 0.5  |     12       |  4.4e2                | **4721** |   5e-12     |
!> | 1.0  |     12       |  2.4e1                |    338   |   2.9e-13   |
!> | 2.0  |     12       |  4.0                  |   95.3   |   7.3e-14   |
!> | 4.0  |     12       |  1.1                  |   57.3   |   6.0e-14   |
!> | 8.0  |     12       |  3.8e-1               |   47.4   |   7.3e-14   |
!>
!> At `s = 0.5` a floor under `1e-10` needs `h > 4 eps * 4721 / 1e-10`, about
!> `5e-2` Bohr. No stencil reaches that -- the grid guard would fire first and
!> the truncation term of any order would be enormous there -- which is why two
!> rounds of step tuning on this case failed and had to.
!>
!> **What made `max|G|` 4721 rather than 35.** The gradient carries
!> `branch_phi_adj` and nothing else on this fixture is remotely that large; the
!> table's first and last rows bracket it, a grid with no branches at all and a
!> grid with the same twelve branch points but a negligible branch adjoint both
!> sitting near 40. The twelve points are six anchor groups of two, each a live
!> branch and one that is all but dead, and the dead one is the problem. At
!> `s = 0.5` it carried `wbranch = 1.5e-4` and `wleb = 8.6e-6` against the live
!> branch's `0.9998` and `0.268`. Since the Gaussian width goes as
!> `xi0 ~ wleb^(-1/2)`, that put `xi0 = 425` next to a live `2.50` -- a factor
!> of 170, which is `sqrt(3.1e4)` exactly. [[branch_point_adjoint]] forms
!> `-0.5 w_xi xi0 / wbranch` and the group reduction multiplies `wbranch` back
!> in, so what survives is `branch_phi_adj ~ w_xi xi0 / (2 s)`: the `wbranch`
!> cancels and the width does not. **A dying branch does not stop contributing
!> as its weight goes to zero -- its adjoint grows, because its Lebedev weight
!> is what is going to zero and `xi0` is that weight to the `-1/2`.**
!>
!> Raising `s` to `2.0` shares the group weight more evenly, which lifts the
!> dead branch's `wleb` and collapses its `xi0`. `max|G|` falls 50-fold, the
!> floor falls with it, and the case passes -- see `FD_STEPS` for the step this
!> then requires, which is larger, not smaller.
!>
!> Two things this is *not*. It is not a cancellation: `adj_m - mean_adj_branch`
!> inside [[compute_branch_phi_adj]] was checked at every branched point and the
!> two differ by four orders, never by less. And it is not a production concern:
!> at the production temperature the same branch is not kept alive at all -- the
!> `s = 0.25` row above already has no branches -- so a branch that survives a
!> sharp softmax is genuinely near-degenerate and has neither a tiny `wleb` nor
!> a huge `xi0`. The pathology belongs to a far solution held open by a softened
!> softmax, which is what this fixture is made of.
!>
!> ## The mutation that shaped the branched case
!>
!> `deff%branch_phi_adj` was zeroed after pass 2 in
!> `derivatives/hessian_response.f90` and the suite re-run:
!>
!>   * `svdw_cross_branching_fd` fails at `4.3069029e-1` absolute and `3.24e-2`
!>     relative -- worst component atom 2 axis 3, direction 1, analytic
!>     `12.8716481` against a reference of `13.3023384`. That is nine and a half
!>     orders above the fixture's floor of `~8e-11`;
!>   * the number is *identical in every printed digit at both steps*,
!>     `4.306903e-1` at `h = 1.1e-3` and at `h = 8.0e-4`. A step-independent
!>     deviation is a missing term, not a stencil artefact, and that is the
!>     cleanest part of the signature;
!>   * every `FIX_PLAIN` case is unchanged, bit for bit, including the two
!>     composites and the all-channel case. `branch_phi_adj` is structurally
!>     zero there, so a driver that drops it is invisible to an unbranched
!>     fixture -- which is why one branched fixture is worth its cost.
!>
!> ## Step and grid guard
!>
!> `FD_STEPS` and `COMPOSITE_STEPS` are measured; the sweeps are in the comments
!> on those parameters. The grid is guarded at every stencil geometry
!> ([[assert_grid_match]]) on `numbering`, `owner`, `branch_count` and
!> `anchor_id`, so a step large enough to re-enumerate the grid fails loudly
!> instead of putting a step into the reference.
!>
!> The guard matters more than it did: the shipped steps are four times the
!> fourth-order ones they replaced, and a seven-point stencil reaches `3h`
!> rather than `2h`. It was swept and does not fire anywhere the sweeps above
!> go -- the largest displacement any of them applies is `3 * 3.0e-3`, nine
!> times the `1e-3` reach of the shipped `FD_STEPS(1)`, on both fixtures and
!> both level-set models.
module test_cavity_drop_hessian_response
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use testdrive, only: new_unittest, unittest_type, error_type, to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_context, only: moist_context_type
   use test_helpers, only: drop_fixture_geometry, build_drop_test_cavity, &
                           LSF_SVDW, LSF_CFC, FIX_PLAIN, FIX_CROSS

   implicit none(type, external)
   private

   public :: collect_cavity_drop_hessian_response

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Directions pushed through in one call: one sparse, one dense
   integer, parameter :: NDIR = 2

   !> Surface-adjoint channels. `w_a` and `w_w` are the two whose fold moves
   !> with the geometry, and are what this half exists to differentiate
   integer, parameter :: CH_XI = 1, CH_F = 2, CH_XYZ = 3, CH_N = 4
   integer, parameter :: CH_K1 = 5, CH_K2 = 6, CH_A = 7, CH_W = 8
   integer, parameter :: NCHAN = 8

   !> Softmax temperature of the branching fixture
   real(wp), parameter :: CROSS_BRANCH_S = 2.0_wp

   !> Branch set the cross fixture is expected to find: six anchor groups of two
   !>
   !> Pinned rather than merely asserted nonzero, because `CROSS_BRANCH_S` also
   !> sets the admissible branch radius --
   !> `branch_rho_cut = log(1/wleb_cut) sqrt(branch_weight_s)` -- so raising it
   !> from `0.5` to `2.0` on 2026-09-07 doubled that radius. Branch *discovery*
   !> on this geometry is known to be platform dependent (the 4.4 Bohr secondary
   !> minimum sits at the edge of the multistart seed rings), and a platform that
   !> admitted a different set would change `max|G|` and with it the round-off
   !> floor the steps are chosen against. This count was stable at
   !> `s = 0.5, 1, 2, 4, 8` here, a sixteen-fold range; if it ever fires
   !> elsewhere, re-measure `max|G|` and the sweep in the comment on `FD_STEPS`
   !> rather than adjusting the number.
   integer, parameter :: CROSS_BRANCHED_PTS = 12
   integer, parameter :: CROSS_MAX_BRANCH = 2

   !> Central-difference steps of the subtraction identity
   !>
   !> Measured over both directions and all four subtraction cases, as the worst
   !> absolute deviation `|analytic - FD|`. A `*` marks a step at which some
   !> component misses the absolute bound but is saved by the relative one, so
   !> the case still passes; a value with no `*` is clean on both:
   !>
   !> | h      | svdw/plain | cfc/plain  | svdw/all   | svdw/cross |
   !> |--------|------------|------------|------------|------------|
   !> | 3.0e-3 | 3.36e-09 * | 4.18e-08 * | 3.35e-09 * | 9.68e-10 * |
   !> | 2.5e-3 | 1.13e-09 * | 1.40e-08 * | 1.13e-09 * | 3.21e-10 * |
   !> | 2.0e-3 | 2.90e-10 * | 3.68e-09 * | 2.85e-10 * | 1.06e-10   |
   !> | 1.5e-3 | 5.42e-11   | 6.63e-10   | 6.87e-11   | 8.62e-11   |
   !> | 1.3e-3 | 4.42e-11   | 2.77e-10   | 4.04e-11   | 1.11e-10   |
   !> | 1.2e-3 | 3.78e-11   | 1.63e-10   | 2.78e-11   | 7.61e-11   |
   !> | 1.1e-3 | 4.77e-11   | 1.02e-10   | 3.13e-11   | 8.20e-11   |
   !> | 1.0e-3 | 8.05e-11   | 6.52e-11   | 4.63e-11   | 1.17e-10   |
   !> | 9.0e-4 | 2.58e-11   | 4.53e-11   | 6.13e-11   | 9.83e-11   |
   !> | 8.0e-4 | 5.59e-11   | 5.21e-11   | 4.34e-11   | 7.87e-11   |
   !> | 7.0e-4 | 7.35e-11   | 4.90e-11   | 1.11e-10   | 7.22e-11   |
   !> | 6.0e-4 | 1.15e-10   | 3.77e-11   | 5.95e-11   | 1.65e-10 * |
   !> | 5.0e-4 | 1.18e-10   | 1.20e-10   | 1.55e-10   | 1.93e-10   |
   !> | 4.0e-4 | 1.00e-10   | 1.26e-10   | 1.40e-10   | 2.70e-10   |
   !> | 3.0e-4 | 2.40e-10   | 1.31e-10   | 1.06e-10   | 2.49e-10 * |
   !>
   !> **Read the table by columns, not by rows.** The two walls move
   !> independently. Truncation is `O(h^6)` and is set by the level-set model:
   !> `cfc/plain` is twelve times `svdw/plain` at `3e-3` and is what closes the
   !> window from above, at about `1.2e-3`. Round-off is `~4 eps max|G| / h` and
   !> is set by the fixture: `svdw/cross` carries the largest `max|G|` of the
   !> four and is what closes it from below, at about `7e-4`. Everything in
   !> between is one flat, jittery bottom near `5e-11`; the ordering inside it is
   !> noise and must not be read as a trend.
   !>
   !> The pair below sits in that window, a factor of 1.375 apart. Both are
   !> required, so a value that agrees at one step only -- the signature of a
   !> step sitting on a wall -- fails. Worst excess over the whole suite,
   !> `max(|d|/FD_ABS, |d|/(FD_REL |ref|))` minimised per component as the test
   !> does, is `0.57`, so both bounds can be tightened by `1.75` before anything
   !> breaks; it is bit-identical across repeated runs and thread counts.
   !>
   !> **Before moving these, read the module header's floor law.** The window is
   !> a property of `max|G|`, and `max|G|` is a property of the fixtures. If a
   !> fixture changes, `4 eps max|G| / 1e-10` is the smallest step that can still
   !> work, and no stencil order will rescue a step below it -- raising the order
   !> only moves the *upper* wall. That is the whole reason this reference is
   !> `O(h^6)`: the steps it needs are ten times larger than the fourth-order
   !> ones it replaced, which is the opposite of the usual instinct.
   !>
   !> The subtraction is why these columns bottom out two orders below what the
   !> sibling suites reach: `G(acc) - G(acc_frozen)` is formed at each stencil
   !> geometry, so everything the two gradients share -- the whole `Phi' . eff`
   !> term, most of the projection's own round-off -- never reaches the
   !> difference quotient. `svdw/cross` benefits like the rest but is still the
   !> column that closes the window from below, for the plain reason that its
   !> `max|G|` is the largest of the four: `95` against `40` on `FIX_PLAIN`.
   !>
   !> **Historical.** At fourth order and `CROSS_BRANCH_S = 0.5` the shipped pair
   !> was `[3.0e-4, 2.5e-4]` and no pair passed the suite; the failure counts of
   !> the last sweep taken there, all other parameters unchanged, were
   !> `[1.5e-4, 1.0e-4] 4`, `[2.0e-4, 1.5e-4] 4`, `[2.5e-4, 2.0e-4] 3`,
   !> `[3.0e-4, 2.5e-4] 1`, `[4.0e-4, 3.0e-4] 2`, `[6.0e-4, 4.0e-4] 4`,
   !> `[1.0e-3, 6.0e-4] 4`. The one survivor was `svdw_cross_branching_fd`, and
   !> the header explains why it could not have been tuned away.
   real(wp), parameter :: FD_STEPS(2) = [1.1E-3_wp, 8.0E-4_wp]

   !> Central-difference steps of the composite assertion
   !>
   !> The composite differences the shipped gradient itself rather than a
   !> difference of two gradients, so nothing cancels before the difference
   !> quotient and its floor is the bare `~4 eps max|G| / h`. Worst absolute
   !> deviation of `FD - (H_fixed . v + response)`, over both directions, at
   !> `O(h^6)`:
   !>
   !> | h      | svdw/plain | cfc/plain  |
   !> |--------|------------|------------|
   !> | 3.0e-3 | 5.23e-08   | 5.82e-08   |
   !> | 2.0e-3 | 4.59e-09   | 5.13e-09   |
   !> | 1.5e-3 | 8.19e-10   | 9.08e-10   |
   !> | 1.2e-3 | 2.20e-10   | 2.32e-10   |
   !> | 1.0e-3 | 1.02e-10   | 9.55e-11   |
   !> | 9.0e-4 | 5.18e-11   | 5.67e-11   |
   !> | 8.0e-4 | 5.92e-11   | 5.25e-11   |
   !> | 7.0e-4 | 2.50e-11   | 7.26e-11   |
   !> | 6.0e-4 | 5.12e-11   | 4.62e-11   |
   !> | 5.0e-4 | 4.57e-11   | 6.56e-11   |
   !> | 4.0e-4 | 7.24e-11   | 8.82e-11   |
   !>
   !> The two models track each other here, unlike in the subtraction sweep:
   !> both walls are properties of the shipped gradient, which is the same
   !> object in both. Truncation closes the window at about `1.1e-3` and the
   !> bottom is flat from there down past `4e-4` -- the round-off wall of
   !> `FIX_PLAIN` is far enough below that this sweep never reaches it. The pair
   !> below sits mid-window, a factor of 1.5 apart, worst excess `0.51`.
   !>
   !> **This sweep is the before-and-after of the stencil change.** At `O(h^4)`
   !> the bowl bottomed at `1.4e-4` and the absolute column never went under
   !> `1e-10` at any step, so the composite passed on the relative bound alone
   !> with a margin of two; the pair was `[1.9e-4, 1.5e-4]`, picked around a
   !> round-off spike at `1.7e-4`. Sixth order takes the bottom of the bowl from
   !> `1.29e-10` to `2.50e-11`, and the shipped pair from `5.06e-10`/`2.60e-10`
   !> to `5.18e-11`/`5.12e-11` -- five-fold at the bottom, five to ten at the
   !> pair -- and widens the window by an order, at four extra cavity builds per
   !> step.
   real(wp), parameter :: COMPOSITE_STEPS(2) = [9.0E-4_wp, 6.0E-4_wp]

   !> Finite-difference agreement bounds
   !>
   !> The project target, `1e-10` absolute and `1e-10` relative; a component
   !> fails only when it misses *both*. **Every case in this suite meets it**,
   !> with `1.75` in hand -- the worst excess over all eight tests, both steps
   !> and both directions is `0.57`, measured by the same `min` of the two
   !> ratios the comparison takes.
   !>
   !> The suite carried two documented exceptions and has neither any more.
   !> `svdw_plain_all_channels_fd` was a deliberate miss at `1.0e-8`, charged to
   !> "the curvature channels' known amplification"; that amplification was the
   !> `sqrt(KM^2 - KG)` cancellation in `kernel.f90`, fixed 2026-09-07.
   !> `svdw_cross_branching_fd` was a deliberate miss at `1.5e-8`, charged to
   !> the branched fixture's own round-off floor; that floor was `max|G|`, and
   !> the module header has what it was made of and what moved it.
   !>
   !> Neither was given a tolerance of its own, and neither should be. A
   !> per-case bound here would have hidden both diagnoses: the curvature one
   !> was a real defect in shipped code, and the branch one was a fixture
   !> parameter that also cost the case nine orders of discriminating power.
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
   !> The only case in which `dbranch_phi_adj` is nonzero, and the only reason
   !> `FIX_CROSS` is built at all. It failed at `1.5e-8` until 2026-09-07, when
   !> `CROSS_BRANCH_S` went from `0.5` to `2.0`; the module header has why that
   !> is a fifty-fold cut in the reference's round-off floor rather than a
   !> loosened fixture. It is now the *most* discriminating case in the suite as
   !> well as a passing one: zeroing the branch channel of pass 2 moves it by
   !> `0.43`, nine and a half orders over its floor and step-independent to
   !> every printed digit.
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

      call drop_fixture_geometry(fix_kind, mol)
      call build_drop_test_cavity(cavity, ctx, mol, fix_kind, lsf_kind, error, &
                                  cross_branch_s=CROSS_BRANCH_S)
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

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_drop_test_cavity(cavity, ctx, mol, FIX_PLAIN, lsf_kind, error)
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

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_drop_test_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error)
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

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_drop_test_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error)
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

      !> Seven-point central stencil of the first derivative, `O(h^6)`
      !>
      !> Sixth order rather than fourth because the round-off floor of this
      !> reference is set by `max|G|` and cannot be lowered by taking a smaller
      !> step; the only way under the project bound is to take a *larger* one,
      !> and that needs the truncation term of a higher-order stencil. See the
      !> comment on `FD_STEPS`.
      integer, parameter :: OFFSET(6) = [-3, -2, -1, 1, 2, 3]
      real(wp), parameter :: COEFF(6) = [-1.0_wp, 9.0_wp, -45.0_wp, &
                                         45.0_wp, -9.0_wp, 1.0_wp]/60.0_wp

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

         call build_drop_test_cavity(cavity, ctx, mol_disp, fix_kind, lsf_kind, error, &
                                     cross_branch_s=CROSS_BRANCH_S)
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

      !> Seven-point central stencil of the first derivative, `O(h^6)`
      !>
      !> The same order the subtraction reference uses, for the same reason and
      !> with a smaller payoff: this one differences the shipped gradient rather
      !> than a difference of two, so its floor was always `eps max|G| / h` and
      !> was already the binding wall at fourth order. The sweep in the comment
      !> on `COMPOSITE_STEPS` has the before and after.
      integer, parameter :: OFFSET(6) = [-3, -2, -1, 1, 2, 3]
      real(wp), parameter :: COEFF(6) = [-1.0_wp, 9.0_wp, -45.0_wp, &
                                         45.0_wp, -9.0_wp, 1.0_wp]/60.0_wp

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

         call build_drop_test_cavity(cavity, ctx, mol_disp, fix_kind, lsf_kind, error, &
                                     cross_branch_s=CROSS_BRANCH_S)
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
         if (count(cavity%branch_count(1:cavity%ngrid) > 1) /= CROSS_BRANCHED_PTS &
             .or. maxval(cavity%branch_count(1:cavity%ngrid)) /= CROSS_MAX_BRANCH) then
            call test_failed(error, "the branching fixture found a different branch set"// &
                             " at "//label//": "// &
                             to_string(count(cavity%branch_count(1:cavity%ngrid) > 1))// &
                             " branched points of max multiplicity "// &
                             to_string(maxval(cavity%branch_count(1:cavity%ngrid)))// &
                             ", expected "//to_string(CROSS_BRANCHED_PTS)//" of "// &
                             to_string(CROSS_MAX_BRANCH)//"; the step choice is tied to"// &
                             " this set through max|G| -- see CROSS_BRANCHED_PTS")
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

end module test_cavity_drop_hessian_response
