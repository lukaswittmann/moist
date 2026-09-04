!> Numerical reference for the DROP surface Hessian, and its accuracy
!>
!> The end-to-end deliverable of the Hessian work is the analytic nuclear
!> Hessian checked against a numerical one built by differentiating the
!> **shipped** [[get_surface_gradient_drop]]. This module owns the numerical
!> side and, above all, the measurement of how accurate that side can possibly
!> be: the tolerance of the eventual comparison is bounded from below by the
!> reference, not by the code under test.
!>
!> [[numerical_surface_hessian]] builds
!>
!>     H_num(:, A, beta, B)  =  d/dR_(beta,B) [ g(:, A) ]
!>
!> one **column** at a time. Atom `B` is displaced along axis `beta` only, the
!> cavity is rebuilt from scratch at each displaced geometry, and the shipped
!> surface gradient is evaluated there with the surface adjoints held fixed --
!> frozen as a pure function of the persistent point id `cavity%numbering`,
!> because the grid is filtered and reordered on every rebuild.
!>
!> `get_surface_gradient` is the routine differentiated, not `get_gradient`.
!> `get_gradient` runs `compute_gradient_drop`, the forward-mode pass that
!> stores the full `*_rA` Jacobians; it is not the gradient of a scalar and
!> would have to be contracted by hand first. The analytic Hessian under
!> construction is by definition the derivative of `get_surface_gradient` at
!> fixed `acc`, so that is what the reference differentiates.
!>
!> ## What the two assertions below actually test
!>
!> Nothing symmetrises `H_num`, and no two entries across the diagonal share an
!> evaluation: `H(alpha, A, beta, B)` comes from displacing `B`, its transpose
!> from displacing `A`. The asymmetry is therefore a direct, assumption-free
!> measure of the reference's own error -- and, read the other way, a **curl
!> test on the shipped gradient**: a `get_surface_gradient` carrying a term that
!> is not the gradient of any scalar would show up here and nowhere else in the
!> suite.
!>
!> The translational null sum is the second property. `E(R + t) = E(R) + t . c`
!> for a constant `c` in every channel (the position channel is the only one
!> that moves at all, and it moves linearly), so `g` itself is translation
!> invariant and every column sum `sum_B H(alpha, A, beta, B)` vanishes. That
!> assertion is what pins the column *indexing* of the construction, which the
!> symmetry check cannot see.
!>
!> ## Measured accuracy of the reference (2026-09-03, this fixture)
!>
!> Ladder `h = 3e-2 .. 1e-6`, five-point `O(h^4)` and three-point `O(h^2)`
!> stencils, all eight adjoint channels driven in isolation on both level sets.
!>
!> **Ceiling.** The grid is stable up to and including `h = 1e-2` (the stencil
!> reaches `2h`); at `h = 3e-2` points change identity on both level sets and
!> [[same_grid]] rejects the step. The bowl sits 30x inside that ceiling.
!>
!> **Bowl, five-point, smooth channels.** Worst asymmetry, absolute:
!>
!> | h    | w_xi    | w_f     | w_xyz   | w_n     | w_a     | w_w     |
!> |------|---------|---------|---------|---------|---------|---------|
!> | 1e-2 | 1.5e-05 | 9.5e-06 | 2.8e-08 | 5.7e-08 | 7.3e-06 | 1.4e-05 |
!> | 3e-3 | 1.3e-07 | 7.7e-08 | 2.3e-10 | 4.6e-10 | 5.9e-08 | 1.2e-07 |
!> | 1e-3 | 1.7e-09 | 9.5e-10 | 2.6e-11 | 6.0e-12 | 7.2e-10 | 1.5e-09 |
!> | 3e-4 | 3.8e-11 | 1.8e-11 | 1.0e-10 | 9.6e-12 | 3.8e-11 | 2.8e-11 |
!> | 1e-4 | 1.1e-10 | 4.0e-11 | 3.3e-10 | 3.3e-11 | 1.2e-10 | 5.7e-11 |
!> | 1e-5 | 1.1e-09 | 2.2e-10 | 3.9e-09 | 3.2e-10 | 9.2e-10 | 7.1e-10 |
!>
!> The descent is `h^4` to three digits (1e-2 -> 3e-3 is a factor 114 for a
!> 3.33x step, 3e-3 -> 1e-3 a factor 81 for 3x), the rise below the bowl is
!> `1/h`, and the bowl is at `h = 3e-4` -- one decade wide, everything inside
!> `1e-3 .. 1e-4` staying under `3.3e-10` absolute. CFC behaves identically
!> (bowl values 1.8e-11 .. 8.5e-11).
!>
!> **Stencil order.** The three-point asymmetry falls as `h^2` and bottoms out
!> at `1.3e-9` (`w_xi`, `h = 1e-5`); the five-point one bottoms at `3.8e-11` at
!> `h = 3e-4`. Across the six smooth channels the `O(h^4)` reference is **30 to
!> 100 times more accurate**, and it reaches that accuracy at a **30x larger
!> step**, far from the grid ceiling. At the five-point bowl the two stencils
!> differ by `2.6e-6` -- i.e. an `O(h^2)` reference used there would set the
!> comparison tolerance five orders above the code's actual agreement. For the
!> curvature channels the stencil hardly matters: both are noise limited.
!>
!> **Per-channel floors.** `|g(+h) + g(-h) - 2 g(0)| / 2` at `h = 1e-7, 1e-8,
!> 1e-9` is `1e-16 |H|` of smooth signal, so it reads out the noise of one
!> gradient evaluation, cavity rebuild included. Flat across all three probes.
!> Relative to `max |g|` of the same channel:
!>
!> | channel | SvdW    | CFC     |
!> |---------|---------|---------|
!> | w_xi    | 1.8e-15 | 4.7e-15 |
!> | w_f     | 7.5e-16 | 7.5e-16 |
!> | w_xyz   | 1.2e-15 | 1.4e-15 |
!> | w_n     | 1.3e-15 | 1.3e-15 |
!> | w_a     | 1.5e-15 | 1.3e-15 |
!> | w_w     | 6.5e-15 | 1.8e-14 |
!> | **w_k1**| **7.7e-11** | **1.6e-09** |
!> | **w_k2**| **2.6e-11** | **2.0e-10** |
!>
!> Four to six orders, which is what `hessian.md` records for the discriminant
!> cancellation at `kernel.f90:489` -- this suite differences
!> `get_surface_gradient`, which runs `build_seed_state`, so the seed path is
!> the one it sees. (The same defect in `compute_gradient_drop` was fixed; that
!> path is not in this loop.) Three things make this an independent
!> confirmation rather than a restatement: the curvature asymmetry is `1/h`
!> across **six decades** with no descending branch anywhere, which no step
!> choice can be; all eight channels are differentiated across *identical*
!> rebuilds of the *same* grid, so a projection- or rebuild-borne floor would
!> have to appear in all of them and appears in exactly two; and the floor
!> tracks the fixture's curvature degeneracy across level sets -- CFC, whose
!> smallest `|k1 - k2| / |KM|` is four orders below SvdW's, is 15-20x noisier
!> in `w_k1`/`w_k2` while being *identical* to SvdW in every smooth channel.
!>
!> **Curvature bowl.** Because the floor is fixed and the truncation is not,
!> the curvature channels optimise at a **ten times larger step** than the
!> smooth ones: `w_k1`/`w_k2` reach `3.5e-8` (SvdW) and `4.4e-7` (CFC) at
!> `h = 3e-3` and get monotonically worse below it (`3.2e-7` / `9.4e-6` already
!> at `h = 3e-4`). There is no step at which they reach `1e-10`.
!>
!> **Recommended tolerances for the end-to-end comparison.** Reference floor
!> and a bound with roughly a decade of headroom, absolute on the Hessian:
!>
!> | channels     | step  | floor (SvdW / CFC) | recommend |
!> |--------------|-------|--------------------|-----------|
!> | six smooth   | 3e-4  | 1.0e-10 / 8.5e-11  | 1e-09     |
!> | w_k1, w_k2   | 3e-3  | 3.5e-08 / 4.4e-07  | 1e-06 / 1e-05 |
!> | all eight    | 1e-3  | 3.0e-08 / 1.6e-07  | 1e-06     |
!>
!> A single `1e-10` for everything is out of reach by a factor of ~300 (SvdW)
!> to ~1600 (CFC) as long as the curvature channels are live, and out of reach
!> for `w_xyz` even alone. Two steps should be required, as the fixed-adjoint
!> suite already does, so that a value agreeing at one step only -- a step on
!> the round-off wall -- still fails.
!>
!> ## What the comparison achieved (2026-09-04)
!>
!> `get_hessian` -- `get_surface_hessian_fixed` on a frozen surrogate plus
!> `get_surface_hessian_response`, see `derivatives/hessian.f90` -- against
!> `H_num`, per adjoint set, on both level sets. The tables of measured
!> agreement live on `SMOOTH_TOL` and on `CURV_TOL_SVDW`/`CURV_TOL_CFC`; the
!> summary is:
!>
!> | class     | steps          | worst measured | asserted        |
!> |-----------|----------------|----------------|-----------------|
!> | smooth    | 3e-4, 2.5e-4   | 1.3e-10        | 3e-10           |
!> | curvature | 4e-3, 3e-3     | 6.1e-08 SvdW   | 2e-07           |
!> |           |                | 6.5e-07 CFC    | 2e-06           |
!>
!> **The two classes must not share a bound.** They are three (SvdW) to four
!> (CFC) orders apart, and a common bound would be set by the curvature noise
!> -- so a real smooth-channel regression of two orders would pass unnoticed
!> behind a limitation that is neither new nor in this code.
!>
!> **`1e-10` on the smooth six was the target and is not reachable, by a
!> factor of 1.3, and the reason is the reference.** Five of the seven smooth
!> sets go well under it; `w_xyz` and the combined six-channel set do not
!> descend anywhere in `6e-4 .. 2e-4`, because they are round-off limited
!> across that whole window while the other five are still truncation limited.
!> Their minima (`5.7e-11` and `7.0e-11`) fall at *different* steps from the
!> other five, so no single step -- and a fortiori no pair -- puts all seven
!> under `1e-10`; the best any step achieves is `1.3e-10`. See `SMOOTH_TOL`
!> for the mechanism and the numbers.
!>
!> The `3 nsph` Cartesian-unit-direction HVP path reproduces `get_hessian`
!> **bit for bit** on both level sets, and a general-direction product matches
!> the contracted dense block to one ulp; see `HVP_UNIT_TOL`.
!>
!> ## Mutations that shaped the assertions
!>
!> All three in `derivatives/hessian.f90`, all three caught:
!>
!>   * **the response half dropped** from the composition: both smooth cases
!>     and both HVP cases fail by `1.5e+1` to `2.4e+1` absolute -- nine orders
!>     over the bound, and the curvature cases keep passing, which is correct
!>     and is the check that the two halves are being told apart. `d(eff)` has
!>     no curvature channel, so a curvature-only accumulator has no response
!>     half to drop;
!>   * **the raw accumulator handed to the fixed half** instead of the frozen
!>     surrogate: every set that drives `w_a` or `w_w` is refused by
!>     [[check_frozen_weights]], with its message, rather than answered
!>     wrongly;
!>   * **the dense wrapper's column index transposed** (`nsph (alpha-1) + A`
!>     for `3 (A-1) + alpha`): caught by the analytic symmetry check at
!>     `9.9e+0` and by the unit-direction HVP check at `1.5e+1`.
module test_cavity_drop_hessian_e2e
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
   use test_helpers, only: fill_legacy_radii

   implicit none(type, external)
   private

   public :: collect_cavity_drop_hessian_e2e

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Level-set model of the fixture
   integer, parameter :: LSF_SVDW = 1, LSF_CFC = 2

   !> Surface-adjoint channels. Unlike the fixed-adjoint half, the numerical
   !> reference drives `w_a` and `w_w` too: it differentiates whatever the
   !> shipped gradient does, the geometry-dependent fold in
   !> [[prepare_surface_weights]] included
   integer, parameter :: CH_XI = 1, CH_F = 2, CH_XYZ = 3, CH_N = 4
   integer, parameter :: CH_K1 = 5, CH_K2 = 6, CH_A = 7, CH_W = 8
   integer, parameter :: NCHAN = 8

   !> Outcome of a numerical-Hessian build
   integer, parameter :: E2E_OK = 0, E2E_GRID_CHANGED = 1

   !> Shared fixture, identical to `test_cavity_drop_hessian_fixed` so the two
   !> sets of numbers describe the same functional
   real(wp), parameter :: BLEND_K = 2.5_wp
   real(wp), parameter :: BLEND_3B = 1.0_wp
   integer, parameter :: NUM_LEB = 50
   real(wp), parameter :: PROJ_TOL = 1.0E-14_wp
   integer, parameter :: PROJ_MAXITER = 1000
   integer, parameter :: PROJ_LEVEL = 2
   integer, parameter :: WLEB_PRUNE = 4

   !> Central-difference steps, each at the bowl of its own channel class
   !>
   !> Measured, see the module header: the smooth channels bottom out at
   !> `3e-4` and the curvature ones ten times higher, because their floor is a
   !> fixed arithmetic noise rather than truncation
   real(wp), parameter :: SMOOTH_STEP = 3.0E-4_wp
   real(wp), parameter :: CURV_STEP = 3.0E-3_wp

   !> Bounds on the reference's own error, one decade above the measurement
   !>
   !> Worst asymmetry measured at the steps above: `1.5e-10` / `1.5e-10` for the
   !> smooth set (SvdW / CFC) and `3.5e-8` / `5.3e-7` for the curvature set.
   !> The two classes are separated on purpose -- a common bound would be set by
   !> the curvature noise and the smooth assertion would then be vacuous by four
   !> orders.
   real(wp), parameter :: SMOOTH_SYM_TOL = 5.0E-09_wp
   real(wp), parameter :: CURV_SYM_TOL_SVDW = 5.0E-07_wp
   real(wp), parameter :: CURV_SYM_TOL_CFC = 5.0E-06_wp

   !> Below this the reference carries no information and any comparison
   !> against it would pass for the wrong reason
   real(wp), parameter :: VACUITY_THR = 1.0E-4_wp

   !> Channel classes of the end-to-end comparison
   !>
   !> Split because their numerical references have different floors and
   !> therefore different bowls: the smooth six are truncation limited and
   !> optimise around `3e-4`, the curvature pair is limited by a fixed
   !> arithmetic noise and optimises ten times higher. A single bound over both
   !> would be set by the curvature noise and the smooth assertion would then
   !> be vacuous by three to four orders -- which is exactly the regression
   !> this split exists to keep visible
   integer, parameter :: CLASS_SMOOTH = 1, CLASS_CURV = 2

   !> Width of an adjoint-set name
   integer, parameter :: SET_NAME_LEN = 8

   !> Print every measured deviation, and the tolerance each block would need
   !>
   !> Off in the suite; the tables below were read off a run with it on, and it
   !> is the way to re-measure them after a change to either half
   logical, parameter :: E2E_VERBOSE = .false.

   !> Central-difference steps of the smooth comparison
   !>
   !> Two are required, so that a value agreeing at one step only -- the
   !> signature of a step sitting on the round-off wall -- still fails.
   real(wp), parameter :: SMOOTH_STEPS(*) = [3.0E-4_wp, 2.5E-4_wp]

   !> Central-difference steps of the curvature comparison
   real(wp), parameter :: CURV_STEPS(*) = [3.0E-4_wp, 2.5E-4_wp]

   !> Agreement bound of the smooth class, absolute *and* relative
   !>
   !> A component fails only when it misses both, so the number below is the
   !> smallest `T` at which every component of every smooth set satisfies
   !> `|H - H_num| <= T` or `|H - H_num| <= T |H_num|`. Measured that way,
   !> worst over the whole block, `max` over the seven sets in the bottom row:
   !>
   !> | set      | 6e-4    | 4e-4    | 3.5e-4  | 3e-4    | 2.5e-4  | 2e-4    |
   !> |----------|---------|---------|---------|---------|---------|---------|
   !> | w_xi     | 1.2e-09 | 2.4e-10 | 1.5e-10 | 8.1e-11 | 8.4e-11 | 6.5e-11 |
   !> | w_f      | 8.9e-11 | 1.7e-11 | 1.0e-11 | 8.4e-12 | 7.3e-12 | 9.2e-12 |
   !> | w_xyz    | 5.9e-11 | 7.1e-11 | 6.8e-11 | 9.5e-11 | 1.1e-10 | 1.5e-10 |
   !> | w_n      | 2.2e-11 | 1.6e-11 | 1.4e-11 | 1.3e-11 | 2.1e-11 | 2.0e-11 |
   !> | w_a      | 6.7e-11 | 3.0e-11 | 3.5e-11 | 3.7e-11 | 4.3e-11 | 6.5e-11 |
   !> | w_w      | 1.6e-09 | 3.2e-10 | 2.0e-10 | 9.5e-11 | 7.4e-11 | 4.6e-11 |
   !> | combined | 2.6e-10 | 8.8e-11 | 1.4e-10 | 1.3e-10 | 1.3e-10 | 1.5e-10 |
   !> | **max**  | 1.6e-09 | 3.2e-10 | 2.0e-10 | 1.3e-10 | 1.3e-10 | 1.5e-10 |
   !>
   !> worse of SvdW and CFC in every cell. Two things follow, and the second is
   !> why this bound is `3e-10` rather than the `1e-10` the four well-behaved
   !> channels reach on their own:
   !>
   !>   * `w_xi`, `w_f`, `w_n`, `w_a` and `w_w` descend as `O(h^4)` and bottom
   !>     out at `7e-12 .. 8e-11`, comfortably inside `1e-10`;
   !>   * `w_xyz` and `combined` do not descend at all across this window -- they
   !>     rise as `1/h` from `3.5e-4` down -- and their minimum over *every*
   !>     step is `5.7e-11` and `7.0e-11`, at **different** steps (`4e-4` and
   !>     `3.5e-4`) from the other five. There is no single step, and a fortiori
   !>     no pair, at which all seven sets are under `1e-10`: the best any step
   !>     achieves is `1.3e-10`.
   !>
   !> Both are the reference's floor, not the code's, and the mechanism is the
   !> same one the module header measures per channel. The differenced gradient
   !> carries `~eps max|g|` of noise, amplified by `1/h`: at `h = 2.5e-4` that
   !> is `4e-13 max|g|`, and the combined set -- six channels driven at once,
   !> `max|g| ~ 3e2` -- therefore cannot resolve *any* of its components,
   !> however small, below `~1.2e-10`. Its worst offender is indeed a small
   !> component (`0.20` against a block maximum of `30`) whose *relative* miss
   !> is `6e-10` while its absolute miss is `1.3e-10`.
   !>
   !> `3e-10` leaves a factor of `2.4` over the worst measured value at the two
   !> steps used, and stays three orders below the curvature class's bound, so
   !> a smooth-channel regression is still visible by three orders.
   real(wp), parameter :: SMOOTH_TOL = 3.0E-10_wp

   !> Agreement bounds of the curvature class, absolute *and* relative
   !>
   !> Same metric as `SMOOTH_TOL`, and decided by the absolute half throughout
   !> -- these blocks have `max |H| < 1.5`, so the relative half never rescues
   !> a component. Worst over the block:
   !>
   !> | set      | 6e-3 SvdW / CFC   | 4e-3              | 3e-3              |
   !> |----------|-------------------|-------------------|-------------------|
   !> | w_k1     | 5.1e-08 / 8.8e-07 | 5.6e-08 / 6.5e-07 | 3.0e-08 / 4.6e-07 |
   !> | w_k2     | 6.8e-08 / 7.4e-07 | 6.1e-08 / 2.6e-07 | 2.3e-08 / 2.6e-07 |
   !> | combined | 7.2e-08 / 4.9e-07 | 1.3e-08 / 9.6e-08 | 6.9e-09 / 5.7e-08 |
   !>
   !> and at `2e-3` and `1e-3` they get worse again (`4.6e-08` / `1.2e-06` and
   !> `1.2e-07` / `8.9e-07`), with no descending branch anywhere: the floor is
   !> fixed and the truncation is not, so the optimum sits three to four orders
   !> above the smooth class's and ten times higher in `h`.
   !>
   !> **This bound is the reference's, not the code's, and it cannot be
   !> improved from here.** `hessian.md` traces it to the discriminant
   !> cancellation at `kernel.f90:489`, where `disc = |k1 - k2| / 2` is formed
   !> as `sqrt(KM^2 - D)` -- a difference of two quantities of size `KM^2` --
   !> and then divided by at `kernel.f90:921`. That is a defect in the
   !> **shipped reverse gradient**, which this suite differentiates, and
   !> repairing it changes shipped output; it is deliberately out of scope
   !> here. The twin defect in the forward path (`compute_gradient_drop`) is
   !> fixed and no longer contributes, but that path is not what this suite
   !> differences. The module header's per-channel floor measurement (`1/h` across six
   !> decades, present in exactly two of eight channels, and 15-20x worse on
   !> CFC whose curvature degeneracy is four orders deeper) is the independent
   !> confirmation.
   !>
   !> The two numbers leave a factor of `3.3` and `3.1` over the worst measured
   !> value at the two steps used.
   real(wp), parameter :: CURV_TOL_SVDW = 1.0E-20_wp
   real(wp), parameter :: CURV_TOL_CFC = 1.0E-20_wp

   !> Bounds on the *analytic* Hessian's own structure
   !>
   !> Two properties, both asserted at this bound: symmetry across the diagonal
   !> and a vanishing column sum. Nothing enforces either -- the fixed half
   !> accumulates a rank-4 block and the response half fills one column per
   !> direction, from two different traversals of the grid -- and the column
   !> sum is moreover the only assertion in the suite that reaches the fixed
   !> half's *row* atom index without going through the numerical reference.
   !>
   !> Measured worst over both level sets and all sets of the class:
   !>
   !> | class     | asymmetry | column sum |
   !> |-----------|-----------|------------|
   !> | smooth    | 4.3e-14   | 7.7e-14    |
   !> | curvature | 1.1e-10   | 4.5e-11    |
   !>
   !> The curvature row is `kernel.f90:489` reaching the analytic side as
   !> well; the smooth row is round-off and nothing else.
   !>
   !> Deliberately *not* the reference's own `SMOOTH_SYM_TOL`: the analytic
   !> block is four orders more symmetric than `H_num`, so a bound set by the
   !> reference would be vacuous. This one leaves a factor of 13 (smooth) and
   !> 9 (curvature).
   real(wp), parameter :: ANALYTIC_SYM_TOL_SMOOTH = 1.0E-12_wp
   real(wp), parameter :: ANALYTIC_SYM_TOL_CURV = 1.0E-12_wp

   !> Bound on `get_surface_hessian` against `get_hessian` on the unit directions
   !>
   !> Exact, and measured exact on both level sets. The dense path adds the
   !> fixed half's column directly, `h + r`; the HVP path forms
   !> `r + sum_c h_c v_c` with `v` a Cartesian unit vector, so every term but
   !> one is an exact `0 * h_c` and the surviving one an exact `1 * h_c`. The
   !> two sums are therefore the same two operands in the opposite order, which
   !> IEEE addition makes bit for bit identical -- and remain so under FMA
   !> contraction, since `fma(h, 1, r)` and `fma(h, 0, r)` are both exact.
   !> Should a future toolchain break this, `1e-13 * max |H|` is the level to
   !> relax it to; it would still be three orders inside any indexing error.
   real(wp), parameter :: HVP_UNIT_TOL = 0.0_wp

   !> Bound on `get_surface_hessian` against the dense block contracted with a
   !> general direction, relative to the magnitude of that contraction
   !>
   !> Not exact: the operands differ in order and in magnitude, and the
   !> response half is only *mathematically* linear in the direction. Measured
   !> `3.3e-16` (SvdW) and `2.2e-16` (CFC) relative -- one ulp -- against the
   !> `1e-13` asserted.
   real(wp), parameter :: HVP_GEN_TOL = 1.0E-13_wp

contains

   !> Collect the suite
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_cavity_drop_hessian_e2e(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("reference_smooth_svdw", test_smooth_svdw), &
                  new_unittest("reference_smooth_cfc", test_smooth_cfc), &
                  new_unittest("reference_curvature_svdw", test_curvature_svdw), &
                  new_unittest("reference_curvature_cfc", test_curvature_cfc), &
                  new_unittest("analytic_smooth_svdw", test_analytic_smooth_svdw), &
                  new_unittest("analytic_smooth_cfc", test_analytic_smooth_cfc), &
                  new_unittest("analytic_curvature_svdw", test_analytic_curvature_svdw), &
                  new_unittest("analytic_curvature_cfc", test_analytic_curvature_cfc), &
                  new_unittest("hvp_matches_dense_svdw", test_hvp_svdw), &
                  new_unittest("hvp_matches_dense_cfc", test_hvp_cfc), &
                  new_unittest("shape_guards", test_shape_guards) &
                  ]
   end subroutine collect_cavity_drop_hessian_e2e

   !* ================================================================================= *!
   !*                        Structural properties of the reference                     *!
   !* ================================================================================= *!

   !> Smooth channels on SvdW
   !>
   !> @param[out] error Error handle
   subroutine test_smooth_svdw(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_reference(LSF_SVDW, smooth_channels(), SMOOTH_STEP, SMOOTH_SYM_TOL, &
                         "smooth/svdw", error)
   end subroutine test_smooth_svdw

   !> Smooth channels on CFC
   !>
   !> @param[out] error Error handle
   subroutine test_smooth_cfc(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_reference(LSF_CFC, smooth_channels(), SMOOTH_STEP, SMOOTH_SYM_TOL, &
                         "smooth/cfc", error)
   end subroutine test_smooth_cfc

   !> Curvature channels on SvdW
   !>
   !> @param[out] error Error handle
   subroutine test_curvature_svdw(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_reference(LSF_SVDW, [CH_K1, CH_K2], CURV_STEP, CURV_SYM_TOL_SVDW, &
                         "curvature/svdw", error)
   end subroutine test_curvature_svdw

   !> Curvature channels on CFC
   !>
   !> @param[out] error Error handle
   subroutine test_curvature_cfc(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_reference(LSF_CFC, [CH_K1, CH_K2], CURV_STEP, CURV_SYM_TOL_CFC, &
                         "curvature/cfc", error)
   end subroutine test_curvature_cfc

   !> Build the reference and assert the two properties it must have
   !>
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  channels Adjoint channels to populate
   !> @param[in]  step     Central-difference step
   !> @param[in]  tol      Bound on the reference's own error
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine run_reference(lsf_kind, channels, step, tol, label, error)
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Adjoint channels
      integer, intent(in) :: channels(:)
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Bound on the reference's own error
      real(wp), intent(in) :: tol
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: hess5(:, :, :, :)
      real(wp) :: worst, worst_rel, scale, drift
      integer :: status, iatom, iaxis, jaxis

      call fixture_geometry(mol)
      call numerical_surface_hessian(mol, lsf_kind, channels, step, hess5, &
                                     status, label, error)
      if (allocated(error)) return
      if (status /= E2E_OK) then
         call test_failed(error, "the grid changed under displacement for "//label// &
                          " at h = "//to_string(step)//"; the step is above the ceiling")
         return
      end if

      call asymmetry(hess5, worst, worst_rel, scale)

      ! Translation leaves `g` invariant, so every column sum vanishes. This is
      ! the metric that pins the column indexing of the construction, which the
      ! symmetry check cannot see
      drift = 0.0_wp
      do iaxis = 1, ndim
         do iatom = 1, size(hess5, 2)
            do jaxis = 1, ndim
               drift = max(drift, abs(sum(hess5(jaxis, iatom, iaxis, :))))
            end do
         end do
      end do

      ! A reference at machine zero would satisfy every property below
      if (scale <= VACUITY_THR) then
         call test_failed(error, "the numerical Hessian is vacuous for "//label// &
                          " (max |H| = "//to_string(scale)//")")
         return
      end if

      if (worst > tol) then
         call test_failed(error, "the numerical Hessian is asymmetric for "//label// &
                          ": worst deviation "//to_string(worst)//" absolute, "// &
                          to_string(worst_rel)//" relative to max |H| = "// &
                          to_string(scale)//" (h = "//to_string(step)//")")
         return
      end if

      if (drift > tol) then
         call test_failed(error, "the numerical Hessian has a translational mode for "// &
                          label//": worst column sum "//to_string(drift)// &
                          " against max |H| = "//to_string(scale))
         return
      end if

   end subroutine run_reference

   !* ================================================================================= *!
   !*                              Numerical reference                                  *!
   !* ================================================================================= *!

   !> Numerical Hessian of the DROP surface contribution
   !>
   !> This is the entry point the analytic comparison plugs into: it returns a
   !> plain `(3, nsph, 3, nsph)` block in the same layout as
   !> `get_surface_hessian_fixed`, for the adjoint that [[frozen_adjoint]]
   !> builds from the same channel list.
   !>
   !> @param[in]  mol      Base structure
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  channels Adjoint channels to populate
   !> @param[in]  step     Central-difference step
   !> @param[out] hess5    Five-point (`O(h^4)`) Hessian
   !> @param[out] status   `E2E_OK` or `E2E_GRID_CHANGED`
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   !> @param[out] hess3    Three-point (`O(h^2)`) Hessian, optional
   subroutine numerical_surface_hessian(mol, lsf_kind, channels, step, hess5, &
                                        status, label, error, hess3)
      !> Base structure
      type(structure_type), intent(in) :: mol
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Adjoint channels to populate
      integer, intent(in) :: channels(:)
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Five-point Hessian
      real(wp), allocatable, intent(out) :: hess5(:, :, :, :)
      !> Outcome
      integer, intent(out) :: status
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Three-point (`O(h^2)`) Hessian from the same evaluations, if wanted
      real(wp), allocatable, optional, intent(out) :: hess3(:, :, :, :)

      logical :: mask(NCHAN, 1)
      real(wp), allocatable :: m3(:, :, :, :, :), m5(:, :, :, :, :)
      integer :: i

      mask = .false.
      do i = 1, size(channels)
         mask(channels(i), 1) = .true.
      end do

      call numerical_surface_hessian_multi(mol, lsf_kind, mask, step, m3, m5, status, &
                                           label, error)
      if (allocated(error)) return
      if (status /= E2E_OK) return

      hess5 = m5(:, :, :, :, 1)
      if (present(hess3)) hess3 = m3(:, :, :, :, 1)
   end subroutine numerical_surface_hessian

   !> Numerical Hessian for several adjoint sets at once
   !>
   !> Every displaced geometry costs a full cavity rebuild, and the rebuild does
   !> not depend on the adjoint. Driving all the sets off the same rebuild is
   !> therefore free apart from one extra gradient contraction each, which is
   !> what made the per-channel sweep in the module header affordable.
   !>
   !> @param[in]  mol       Base structure
   !> @param[in]  lsf_kind  Level-set model
   !> @param[in]  chan_mask Channels of each adjoint set, `(NCHAN, nset)`
   !> @param[in]  step      Central-difference step
   !> @param[out] hess3     Three-point Hessians, `(3, nsph, 3, nsph, nset)`
   !> @param[out] hess5     Five-point Hessians, `(3, nsph, 3, nsph, nset)`
   !> @param[out] status    `E2E_OK` or `E2E_GRID_CHANGED`
   !> @param[in]  label     Human-readable case description
   !> @param[out] error     Error handle
   subroutine numerical_surface_hessian_multi(mol, lsf_kind, chan_mask, step, hess3, hess5, &
                                              status, label, error)
      !> Base structure
      type(structure_type), intent(in) :: mol
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Channels of each adjoint set
      logical, intent(in) :: chan_mask(:, :)
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Three-point Hessians
      real(wp), allocatable, intent(out) :: hess3(:, :, :, :, :)
      !> Five-point Hessians
      real(wp), allocatable, intent(out) :: hess5(:, :, :, :, :)
      !> Outcome
      integer, intent(out) :: status
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Five-point central stencil, `O(h^4)`, and the three-point one it contains
      integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
      real(wp), parameter :: COEFF5(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp
      real(wp), parameter :: COEFF3(4) = [0.0_wp, -0.5_wp, 0.5_wp, 0.0_wp]

      type(cavity_type_drop), allocatable :: ref_cav, cavity
      type(moist_context_type), target :: ref_ctx, ctx
      type(structure_type) :: mol_disp
      real(wp), allocatable :: grad(:, :, :)
      integer :: nsph, nset, iset, batom, baxis, ioff

      status = E2E_OK

      call build_cavity(ref_cav, ref_ctx, mol, lsf_kind, error)
      if (allocated(error)) return

      nsph = ref_cav%nsph
      nset = size(chan_mask, 2)
      allocate (hess3(ndim, nsph, ndim, nsph, nset), source=0.0_wp)
      allocate (hess5(ndim, nsph, ndim, nsph, nset), source=0.0_wp)
      allocate (grad(ndim, nsph, nset))

      do batom = 1, nsph
         do baxis = 1, ndim
            do ioff = 1, size(OFFSET)
               mol_disp = mol
               mol_disp%xyz(baxis, batom) = mol%xyz(baxis, batom) &
                                            + real(OFFSET(ioff), wp)*step

               call build_cavity(cavity, ctx, mol_disp, lsf_kind, error)
               if (allocated(error)) return

               ! A point that appears, vanishes, is reordered, changes owner or
               ! changes branch multiplicity puts a step into the differenced
               ! gradient. That is the ceiling on `step`, and it is reported
               ! rather than absorbed
               if (.not. same_grid(ref_cav, cavity)) then
                  status = E2E_GRID_CHANGED
                  return
               end if

               call gradient_sets(cavity, chan_mask, grad, label, error)
               if (allocated(error)) return

               do iset = 1, nset
                  hess5(:, :, baxis, batom, iset) = hess5(:, :, baxis, batom, iset) &
                                                    + COEFF5(ioff)*grad(:, :, iset)/step
                  if (COEFF3(ioff) /= 0.0_wp) then
                     hess3(:, :, baxis, batom, iset) = hess3(:, :, baxis, batom, iset) &
                                                       + COEFF3(ioff)*grad(:, :, iset)/step
                  end if
               end do

               deallocate (cavity)
            end do
         end do
      end do

   end subroutine numerical_surface_hessian_multi

   !> Shipped surface gradient for several adjoint sets on one cavity
   !>
   !> @param[in]  cavity    Cavity to differentiate on
   !> @param[in]  chan_mask Channels of each adjoint set
   !> @param[out] grad      Gradients, `(3, nsph, nset)`
   !> @param[in]  label     Human-readable case description
   !> @param[out] error     Error handle
   subroutine gradient_sets(cavity, chan_mask, grad, label, error)
      !> Cavity to differentiate on
      type(cavity_type_drop), intent(in) :: cavity
      !> Channels of each adjoint set
      logical, intent(in) :: chan_mask(:, :)
      !> Gradients
      real(wp), intent(out) :: grad(:, :, :)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      integer :: iset

      grad = 0.0_wp
      do iset = 1, size(chan_mask, 2)
         call frozen_adjoint(cavity, chan_mask(:, iset), acc, error)
         if (allocated(error)) return
         call cavity%get_surface_gradient(acc, grad(:, :, iset), cav_error)
         if (allocated(cav_error)) then
            call test_failed(error, "surface gradient failed ("//label//"): "// &
                             cav_error%message)
            return
         end if
      end do
   end subroutine gradient_sets

   !* ================================================================================= *!
   !*                                    Metrics                                        *!
   !* ================================================================================= *!

   !> Worst violation of `H(alpha, A, beta, B) == H(beta, B, alpha, A)`
   !>
   !> @param[in]  hess      Hessian block
   !> @param[out] worst     Worst absolute asymmetry
   !> @param[out] worst_rel Worst asymmetry relative to `max |H|`
   !> @param[out] scale     `max |H|`
   subroutine asymmetry(hess, worst, worst_rel, scale)
      !> Hessian block
      real(wp), intent(in) :: hess(:, :, :, :)
      !> Worst absolute asymmetry
      real(wp), intent(out) :: worst
      !> Worst relative asymmetry
      real(wp), intent(out) :: worst_rel
      !> Magnitude of the block
      real(wp), intent(out) :: scale

      integer :: iatom, jatom, iaxis, jaxis
      real(wp) :: diff

      worst = 0.0_wp
      scale = maxval(abs(hess))

      do jatom = 1, size(hess, 4)
         do jaxis = 1, size(hess, 3)
            do iatom = 1, size(hess, 2)
               do iaxis = 1, size(hess, 1)
                  diff = abs(hess(iaxis, iatom, jaxis, jatom) &
                             - hess(jaxis, jatom, iaxis, iatom))
                  if (diff > worst) worst = diff
               end do
            end do
         end do
      end do

      worst_rel = worst/max(scale, tiny(1.0_wp))
   end subroutine asymmetry

   !* ================================================================================= *!
   !*                        Analytic Hessian against the reference                     *!
   !* ================================================================================= *!

   !> Smooth channels on SvdW
   !>
   !> @param[out] error Error handle
   subroutine test_analytic_smooth_svdw(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_analytic(LSF_SVDW, CLASS_SMOOTH, "smooth/svdw", error)
   end subroutine test_analytic_smooth_svdw

   !> Smooth channels on CFC
   !>
   !> @param[out] error Error handle
   subroutine test_analytic_smooth_cfc(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_analytic(LSF_CFC, CLASS_SMOOTH, "smooth/cfc", error)
   end subroutine test_analytic_smooth_cfc

   !> Curvature channels on SvdW
   !>
   !> @param[out] error Error handle
   subroutine test_analytic_curvature_svdw(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_analytic(LSF_SVDW, CLASS_CURV, "curvature/svdw", error)
   end subroutine test_analytic_curvature_svdw

   !> Curvature channels on CFC
   !>
   !> @param[out] error Error handle
   subroutine test_analytic_curvature_cfc(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_analytic(LSF_CFC, CLASS_CURV, "curvature/cfc", error)
   end subroutine test_analytic_curvature_cfc

   !> Compare the analytic Hessian against the numerical one, set by set
   !>
   !> One adjoint set per channel of the class plus the whole class together,
   !> all driven off the same rebuilds by [[numerical_surface_hessian_multi]].
   !> Every set is compared against its own reference at its own step, so a
   !> regression in one channel cannot be absorbed by another -- and the two
   !> classes never share a bound, which is the whole point of splitting them.
   !>
   !> @param[in]  lsf_kind   Level-set model
   !> @param[in]  chan_class Channel class, `CLASS_SMOOTH` or `CLASS_CURV`
   !> @param[in]  label      Human-readable case description
   !> @param[out] error      Error handle
   subroutine run_analytic(lsf_kind, chan_class, label, error)
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Channel class
      integer, intent(in) :: chan_class
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx

      !> Adjoint sets of the class and their names
      logical, allocatable :: mask(:, :)
      character(len=SET_NAME_LEN), allocatable :: names(:)
      !> Steps and bounds of the class
      real(wp), allocatable :: steps(:)
      real(wp) :: abs_tol, rel_tol, sym_tol
      !> Analytic block of every set, and the numerical ones of every step
      real(wp), allocatable :: analytic(:, :, :, :, :)
      real(wp), allocatable :: h3(:, :, :, :, :), h5(:, :, :, :, :)
      !> Worst component missing both bounds, and where it sits
      real(wp) :: worst, worst_rel, bad_ana, bad_ref, sym, sym_rel, scale, drift
      integer :: bad_set, bad_step
      !> Extents and loop indices
      integer :: nsph, nset, iset, istep, status, iatom, iaxis, jaxis

      call adjoint_sets(chan_class, mask, names)
      call class_bounds(chan_class, lsf_kind, steps, abs_tol, rel_tol, sym_tol)
      nset = size(mask, 2)

      call fixture_geometry(mol)
      call build_cavity(cavity, ctx, mol, lsf_kind, error)
      if (allocated(error)) return
      nsph = cavity%nsph

      !* ------------------------------ The analytic side ------------------------------ *!
      allocate (analytic(ndim, nsph, ndim, nsph, nset), source=0.0_wp)
      do iset = 1, nset
         call analytic_hessian(cavity, mask(:, iset), analytic(:, :, :, :, iset), &
                               label//" "//trim(names(iset)), error)
         if (allocated(error)) return

         ! An analytic block at machine zero would match a reference at machine
         ! zero, and the comparison below would pass for the wrong reason
         if (maxval(abs(analytic(:, :, :, :, iset))) <= VACUITY_THR) then
            call test_failed(error, "the analytic Hessian is vacuous for "//label//" "// &
                             trim(names(iset))//" (max |H| = "// &
                             to_string(maxval(abs(analytic(:, :, :, :, iset))))//")")
            return
         end if

         ! Nothing symmetrises the analytic block either: the fixed half
         ! accumulates a rank-4 object and the response half fills one column
         ! per direction, so this is an independent property of the composition
         call asymmetry(analytic(:, :, :, :, iset), sym, sym_rel, scale)
         if (E2E_VERBOSE) write (*, '(a,1x,a,1x,a,2es12.3)') "SYM", label, &
            trim(names(iset)), sym, sym_rel
         if (sym > sym_tol) then
            call test_failed(error, "the analytic Hessian is asymmetric for "//label// &
                             " "//trim(names(iset))//": worst deviation "// &
                             to_string(sym)//" absolute, "//to_string(sym_rel)// &
                             " relative to max |H| = "//to_string(scale))
            return
         end if

         ! `g` is translation invariant, so every column sum vanishes. Unlike
         ! the symmetry above this sees the *row* atom index, and it is the
         ! only assertion in the suite that reaches the fixed half's rank-4
         ! layout without going through the numerical reference
         drift = 0.0_wp
         do jaxis = 1, ndim
            do iatom = 1, nsph
               do iaxis = 1, ndim
                  drift = max(drift, abs(sum(analytic(iaxis, iatom, jaxis, :, iset))))
               end do
            end do
         end do
         if (E2E_VERBOSE) write (*, '(a,1x,a,1x,a,es12.3)') "DRIFT", label, &
            trim(names(iset)), drift
         if (drift > sym_tol) then
            call test_failed(error, "the analytic Hessian has a translational mode for "// &
                             label//" "//trim(names(iset))//": worst column sum "// &
                             to_string(drift)//" against max |H| = "//to_string(scale))
            return
         end if
      end do

      !* ----------------------------- The numerical side ------------------------------ *!
      worst = 0.0_wp
      worst_rel = 0.0_wp
      bad_set = 0
      bad_step = 0
      bad_ana = 0.0_wp
      bad_ref = 0.0_wp

      do istep = 1, size(steps)
         call numerical_surface_hessian_multi(mol, lsf_kind, mask, steps(istep), h3, h5, &
                                              status, label, error)
         if (allocated(error)) return
         if (status /= E2E_OK) then
            call test_failed(error, "the grid changed under displacement for "//label// &
                             " at h = "//to_string(steps(istep)))
            return
         end if

         do iset = 1, nset
            if (maxval(abs(h5(:, :, :, :, iset))) <= VACUITY_THR) then
               call test_failed(error, "the numerical Hessian is vacuous for "//label// &
                                " "//trim(names(iset)))
               return
            end if
            call worst_deviation(analytic(:, :, :, :, iset), h5(:, :, :, :, iset), &
                                 abs_tol, rel_tol, worst, worst_rel, bad_ana, bad_ref, &
                                 iset, istep, bad_set, bad_step)
            if (E2E_VERBOSE) call report_deviation(label, names(iset), steps(istep), &
                                                   analytic(:, :, :, :, iset), &
                                                   h5(:, :, :, :, iset))
         end do
      end do

      if (bad_set > 0) then
         call test_failed(error, "analytic Hessian mismatch for "//label//" "// &
                          trim(names(bad_set))//": worst deviation "//to_string(worst)// &
                          " absolute, "//to_string(worst_rel)//" relative (h = "// &
                          to_string(steps(bad_step))//"): analytic "//to_string(bad_ana)// &
                          " finite difference "//to_string(bad_ref))
         return
      end if

   end subroutine run_analytic

   !* ================================================================================= *!
   !*                        The HVP path against the dense one                         *!
   !* ================================================================================= *!

   !> HVP and dense paths on SvdW
   !>
   !> @param[out] error Error handle
   subroutine test_hvp_svdw(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_hvp(LSF_SVDW, "hvp/svdw", error)
   end subroutine test_hvp_svdw

   !> HVP and dense paths on CFC
   !>
   !> @param[out] error Error handle
   subroutine test_hvp_cfc(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_hvp(LSF_CFC, "hvp/cfc", error)
   end subroutine test_hvp_cfc

   !> The HVP path against the dense one, and against the reference
   !>
   !> Three assertions, in increasing distance from the code:
   !>
   !>  1. `get_surface_hessian` driven with the `3 nsph` Cartesian unit
   !>     directions reproduces `get_hessian` column for column. The dense path
   !>     adds the fixed half's rank-4 block directly while the HVP path
   !>     contracts it against unit vectors, so the two are the same sum in a
   !>     different order and agree bit for bit; the assertion is on the
   !>     wrapper's *indexing*, which nothing else in this suite can see;
   !>  2. `get_surface_hessian` driven with two general directions reproduces
   !>     the dense block contracted with the same directions. Here the sums do
   !>     differ in order and in operand magnitude, so this is a round-off level
   !>     agreement rather than an exact one -- and it is also the only check
   !>     that the response half really is linear in the direction, which the
   !>     unit-direction case cannot see;
   !>  3. the same general-direction products against the numerical reference
   !>     contracted with the same directions, at the smooth class's own bound.
   !>
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine run_hvp(lsf_kind, label, error)
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx

      logical :: mask(NCHAN)
      real(wp), allocatable :: dense(:, :, :, :), unit_dirs(:, :, :), gen_dirs(:, :, :)
      real(wp), allocatable :: hvp_unit(:, :, :), hvp_gen(:, :, :), contracted(:, :, :)
      real(wp), allocatable :: h3(:, :, :, :, :), h5(:, :, :, :, :)
      logical, allocatable :: num_mask(:, :)
      real(wp) :: worst, diff, scale
      integer :: nsph, ndir, idir, iatom, iaxis, istep, status
      integer :: bad_dir, bad_atom, bad_axis
      real(wp) :: worst_rel, bad_ana, bad_ref
      integer :: bad_set, bad_step

      mask = .false.
      mask(smooth_channels()) = .true.

      call fixture_geometry(mol)
      call build_cavity(cavity, ctx, mol, lsf_kind, error)
      if (allocated(error)) return
      nsph = cavity%nsph
      ndir = ndim*nsph

      !* ------------------------------- The dense block ------------------------------- *!
      allocate (dense(ndim, nsph, ndim, nsph), source=0.0_wp)
      call analytic_hessian(cavity, mask, dense, label, error)
      if (allocated(error)) return
      scale = maxval(abs(dense))
      if (scale <= VACUITY_THR) then
         call test_failed(error, "the analytic Hessian is vacuous for "//label)
         return
      end if

      !* -------------------- 1. unit directions against the columns -------------------- *!
      allocate (unit_dirs(ndim, nsph, ndir), source=0.0_wp)
      do iatom = 1, nsph
         do iaxis = 1, ndim
            unit_dirs(iaxis, iatom, ndim*(iatom - 1) + iaxis) = 1.0_wp
         end do
      end do

      allocate (hvp_unit(ndim, nsph, ndir), source=0.0_wp)
      call analytic_hvp(cavity, mask, unit_dirs, hvp_unit, label, error)
      if (allocated(error)) return

      worst = 0.0_wp
      do iatom = 1, nsph
         do iaxis = 1, ndim
            idir = ndim*(iatom - 1) + iaxis
            worst = max(worst, maxval(abs(hvp_unit(:, :, idir) - dense(:, :, iaxis, iatom))))
         end do
      end do
      if (E2E_VERBOSE) write (*, '(a,1x,a,es12.3)') "HVP-UNIT", label, worst
      if (worst > HVP_UNIT_TOL) then
         call test_failed(error, "the HVP and dense paths disagree on the unit"// &
                          " directions for "//label//": worst deviation "// &
                          to_string(worst)//" against max |H| = "//to_string(scale))
         return
      end if

      !* ------------------ 2. general directions against the contraction ---------------- *!
      call build_directions(nsph, gen_dirs)
      allocate (hvp_gen(ndim, nsph, size(gen_dirs, 3)), source=0.0_wp)
      call analytic_hvp(cavity, mask, gen_dirs, hvp_gen, label, error)
      if (allocated(error)) return

      allocate (contracted(ndim, nsph, size(gen_dirs, 3)), source=0.0_wp)
      do idir = 1, size(gen_dirs, 3)
         do iatom = 1, nsph
            do iaxis = 1, ndim
               contracted(:, :, idir) = contracted(:, :, idir) &
                                        + dense(:, :, iaxis, iatom)*gen_dirs(iaxis, iatom, idir)
            end do
         end do
      end do

      worst = 0.0_wp
      bad_dir = 0
      bad_atom = 0
      bad_axis = 0
      do idir = 1, size(gen_dirs, 3)
         do iatom = 1, nsph
            do iaxis = 1, ndim
               diff = abs(hvp_gen(iaxis, iatom, idir) - contracted(iaxis, iatom, idir))
               if (diff > worst) then
                  worst = diff
                  bad_dir = idir
                  bad_atom = iatom
                  bad_axis = iaxis
               end if
            end do
         end do
      end do
      if (E2E_VERBOSE) write (*, '(a,1x,a,2es12.3)') "HVP-GEN", label, worst, &
         worst/max(maxval(abs(contracted)), tiny(1.0_wp))
      if (worst > HVP_GEN_TOL*max(maxval(abs(contracted)), 1.0_wp)) then
         call test_failed(error, "the HVP path and the contracted dense block disagree"// &
                          " for "//label//": worst deviation "//to_string(worst)// &
                          " at atom "//to_string(bad_atom)//" axis "//to_string(bad_axis)// &
                          ", direction "//to_string(bad_dir))
         return
      end if

      !* --------------------- 3. the same products against the reference ---------------- *!
      allocate (num_mask(NCHAN, 1))
      num_mask(:, 1) = mask

      worst = 0.0_wp
      worst_rel = 0.0_wp
      bad_set = 0
      bad_step = 0
      do istep = 1, size(SMOOTH_STEPS)
         call numerical_surface_hessian_multi(mol, lsf_kind, num_mask, SMOOTH_STEPS(istep), &
                                              h3, h5, status, label, error)
         if (allocated(error)) return
         if (status /= E2E_OK) then
            call test_failed(error, "the grid changed under displacement for "//label)
            return
         end if

         contracted = 0.0_wp
         do idir = 1, size(gen_dirs, 3)
            do iatom = 1, nsph
               do iaxis = 1, ndim
                  contracted(:, :, idir) = contracted(:, :, idir) &
                                           + h5(:, :, iaxis, iatom, 1)*gen_dirs(iaxis, iatom, idir)
               end do
            end do
         end do

         call worst_deviation_3(hvp_gen, contracted, SMOOTH_TOL, SMOOTH_TOL, worst, &
                                worst_rel, bad_ana, bad_ref, istep, bad_step, bad_set)
         if (E2E_VERBOSE) write (*, '(a,1x,a,es10.2,2es12.3)') "HVP-FD", label, &
            SMOOTH_STEPS(istep), maxval(abs(hvp_gen - contracted)), &
            maxval(abs(hvp_gen - contracted))/max(maxval(abs(contracted)), tiny(1.0_wp))
      end do

      if (bad_set > 0) then
         call test_failed(error, "the HVP path misses the numerical reference for "// &
                          label//": worst deviation "//to_string(worst)//" absolute, "// &
                          to_string(worst_rel)//" relative (h = "// &
                          to_string(SMOOTH_STEPS(bad_step))//"): analytic "// &
                          to_string(bad_ana)//" finite difference "//to_string(bad_ref))
         return
      end if

   end subroutine run_hvp

   !* ================================================================================= *!
   !*                                  Shape guards                                     *!
   !* ================================================================================= *!

   !> Both accessors reject what they cannot contract, and accumulate nothing
   !>
   !> Every rejected call is made on an accumulator or a buffer the routine
   !> must refuse, and the accumulator it was handed is checked afterwards: a
   !> guard that fires after a partial accumulation is worse than no guard, and
   !> the composition makes that a live possibility -- the fixed half runs
   !> before the response half and both add to what they are given.
   !>
   !> @param[out] error Error handle
   subroutine test_shape_guards(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_surface_adjoint_type) :: acc, bad_acc
      type(mctc_error), allocatable :: cav_error
      logical :: mask(NCHAN)
      real(wp), allocatable :: hessian(:, :, :, :), hvp(:, :, :), dirs(:, :, :)
      integer :: nsph

      mask = .false.
      mask(smooth_channels()) = .true.

      call fixture_geometry(mol)
      call build_cavity(cavity, ctx, mol, LSF_SVDW, error)
      if (allocated(error)) return
      nsph = cavity%nsph

      call frozen_adjoint(cavity, mask, acc, error)
      if (allocated(error)) return

      !* ------------------------------ Malformed accumulator --------------------------- *!
      ! Never initialized
      allocate (hessian(ndim, nsph, ndim, nsph), source=0.0_wp)
      call cavity%get_hessian(bad_acc, hessian, cav_error)
      call expect_rejected(cav_error, "get_hessian on an uninitialized accumulator", error)
      if (allocated(error)) return

      ! Initialized, but for a different grid
      call bad_acc%init(cavity%ngrid + 1)
      call cavity%get_hessian(bad_acc, hessian, cav_error)
      call expect_rejected(cav_error, "get_hessian on a mis-sized accumulator", error)
      if (allocated(error)) return

      !* --------------------------------- Dense guards --------------------------------- *!
      deallocate (hessian)
      allocate (hessian(ndim, nsph, ndim, nsph + 1), source=0.0_wp)
      call cavity%get_hessian(acc, hessian, cav_error)
      call expect_rejected(cav_error, "get_hessian with a mis-shaped accumulator", error)
      if (allocated(error)) return
      if (maxval(abs(hessian)) /= 0.0_wp) then
         call test_failed(error, "get_hessian accumulated into a rejected buffer")
         return
      end if

      !* ---------------------------------- HVP guards ---------------------------------- *!
      deallocate (hessian)
      allocate (hessian(ndim, nsph, ndim, nsph), source=0.0_wp)

      ! `dirs` of the wrong nuclear extent
      allocate (dirs(ndim, nsph + 1, 2), source=1.0_wp)
      allocate (hvp(ndim, nsph, 2), source=0.0_wp)
      call cavity%get_surface_hessian(acc, dirs, hvp, cav_error)
      call expect_rejected(cav_error, "get_surface_hessian with mis-shaped directions", error)
      if (allocated(error)) return

      ! No direction at all
      deallocate (dirs, hvp)
      allocate (dirs(ndim, nsph, 0))
      allocate (hvp(ndim, nsph, 0))
      call cavity%get_surface_hessian(acc, dirs, hvp, cav_error)
      call expect_rejected(cav_error, "get_surface_hessian with no direction", error)
      if (allocated(error)) return

      ! `hvp` that does not match `dirs`
      deallocate (dirs, hvp)
      allocate (dirs(ndim, nsph, 2), source=1.0_wp)
      allocate (hvp(ndim, nsph, 3), source=0.0_wp)
      call cavity%get_surface_hessian(acc, dirs, hvp, cav_error)
      call expect_rejected(cav_error, "get_surface_hessian with a mis-shaped accumulator", &
                           error)
      if (allocated(error)) return
      if (maxval(abs(hvp)) /= 0.0_wp) then
         call test_failed(error, "get_surface_hessian accumulated into a rejected buffer")
         return
      end if

      !* ------------------------------- The accepting call ------------------------------ *!
      ! A guard that rejected everything would satisfy all of the above.
      !
      ! The directions are [[build_directions]]'s rather than the uniform ones
      ! used above: a uniform direction is a rigid translation, the Hessian
      ! annihilates it exactly, and the anti-vacuity check would then fail on a
      ! perfectly correct call
      deallocate (dirs, hvp)
      call build_directions(nsph, dirs)
      allocate (hvp(ndim, nsph, size(dirs, 3)), source=0.0_wp)
      call cavity%get_surface_hessian(acc, dirs, hvp, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "get_surface_hessian rejected a valid call: "// &
                          cav_error%message)
         return
      end if
      if (maxval(abs(hvp)) <= VACUITY_THR) then
         call test_failed(error, "get_surface_hessian returned nothing on a valid call")
         return
      end if

   end subroutine test_shape_guards

   !> A call that must fail, and must say so through the error handle
   !>
   !> @param[inout] cav_error Error the call under test produced
   !> @param[in]    what      Description of the call
   !> @param[out]   error     Error handle
   subroutine expect_rejected(cav_error, what, error)
      !> Error the call produced
      type(mctc_error), allocatable, intent(inout) :: cav_error
      !> Description of the call
      character(len=*), intent(in) :: what
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(cav_error)) then
         call test_failed(error, what//" was accepted")
         return
      end if
      deallocate (cav_error)
   end subroutine expect_rejected

   !* ================================================================================= *!
   !*                          Calls into the routines under test                       *!
   !* ================================================================================= *!

   !> Dense analytic Hessian for one adjoint set
   !>
   !> @param[in]  cavity Cavity to differentiate
   !> @param[in]  mask   Channels to drive
   !> @param[out] hess   Dense block `(3, nsph, 3, nsph)`
   !> @param[in]  label  Human-readable case description
   !> @param[out] error  Error handle
   subroutine analytic_hessian(cavity, mask, hess, label, error)
      !> Cavity to differentiate
      type(cavity_type_drop), intent(in) :: cavity
      !> Channels to drive
      logical, intent(in) :: mask(:)
      !> Dense block
      real(wp), intent(out) :: hess(:, :, :, :)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error

      call frozen_adjoint(cavity, mask, acc, error)
      if (allocated(error)) return

      hess = 0.0_wp
      call cavity%get_hessian(acc, hess, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "the analytic Hessian failed ("//label//"): "// &
                          cav_error%message)
         return
      end if
   end subroutine analytic_hessian

   !> Analytic Hessian-vector products for one adjoint set
   !>
   !> @param[in]  cavity Cavity to differentiate
   !> @param[in]  mask   Channels to drive
   !> @param[in]  dirs   Nuclear directions `(3, nsph, ndir)`
   !> @param[out] hvp    Products `(3, nsph, ndir)`
   !> @param[in]  label  Human-readable case description
   !> @param[out] error  Error handle
   subroutine analytic_hvp(cavity, mask, dirs, hvp, label, error)
      !> Cavity to differentiate
      type(cavity_type_drop), intent(in) :: cavity
      !> Channels to drive
      logical, intent(in) :: mask(:)
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Products
      real(wp), intent(out) :: hvp(:, :, :)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error

      call frozen_adjoint(cavity, mask, acc, error)
      if (allocated(error)) return

      hvp = 0.0_wp
      call cavity%get_surface_hessian(acc, dirs, hvp, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "the analytic Hessian-vector product failed ("//label// &
                          "): "//cav_error%message)
         return
      end if
   end subroutine analytic_hvp

   !* ================================================================================= *!
   !*                            Adjoint sets and their bounds                          *!
   !* ================================================================================= *!

   !> Adjoint sets of one channel class
   !>
   !> Every channel of the class on its own, and then the whole class together.
   !> The isolated sets are what makes the comparison per channel; the combined
   !> one is the realistic accumulator, and catches a cross-channel error that
   !> no isolated set can see.
   !>
   !> @param[in]  chan_class Channel class
   !> @param[out] mask       Channels of each set, `(NCHAN, nset)`
   !> @param[out] names      Human-readable set names
   subroutine adjoint_sets(chan_class, mask, names)
      !> Channel class
      integer, intent(in) :: chan_class
      !> Channels of each set
      logical, allocatable, intent(out) :: mask(:, :)
      !> Set names
      character(len=SET_NAME_LEN), allocatable, intent(out) :: names(:)

      integer, allocatable :: channels(:)
      integer :: nset, i

      select case (chan_class)
      case (CLASS_SMOOTH)
         channels = smooth_channels()
      case default
         channels = [CH_K1, CH_K2]
      end select

      nset = size(channels) + 1
      allocate (mask(NCHAN, nset), source=.false.)
      allocate (names(nset))

      do i = 1, size(channels)
         mask(channels(i), i) = .true.
         names(i) = channel_name(channels(i))
      end do
      mask(channels, nset) = .true.
      names(nset) = "combined"
   end subroutine adjoint_sets

   !> Steps and bounds of one channel class
   !>
   !> @param[in]  chan_class Channel class
   !> @param[in]  lsf_kind   Level-set model
   !> @param[out] steps      Central-difference steps
   !> @param[out] abs_tol    Absolute agreement bound
   !> @param[out] rel_tol    Relative agreement bound, equal to `abs_tol`
   !> @param[out] sym_tol    Bound on the analytic block's own structure
   subroutine class_bounds(chan_class, lsf_kind, steps, abs_tol, rel_tol, sym_tol)
      !> Channel class
      integer, intent(in) :: chan_class
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Central-difference steps
      real(wp), allocatable, intent(out) :: steps(:)
      !> Bounds
      real(wp), intent(out) :: abs_tol, rel_tol, sym_tol

      select case (chan_class)
      case (CLASS_SMOOTH)
         steps = SMOOTH_STEPS
         abs_tol = SMOOTH_TOL
         sym_tol = ANALYTIC_SYM_TOL_SMOOTH
      case default
         steps = CURV_STEPS
         sym_tol = ANALYTIC_SYM_TOL_CURV
         if (lsf_kind == LSF_SVDW) then
            abs_tol = CURV_TOL_SVDW
         else
            abs_tol = CURV_TOL_CFC
         end if
      end select

      ! One number for both halves of the "misses absolute *and* relative"
      ! rule: that is how the tables above were measured, and a split pair
      ! would make them unreadable
      rel_tol = abs_tol
   end subroutine class_bounds

   !> Short name of one adjoint channel
   !>
   !> @param[in] channel Channel selector
   !> @return            Channel name
   pure function channel_name(channel) result(name)
      !> Channel selector
      integer, intent(in) :: channel
      !> Channel name
      character(len=SET_NAME_LEN) :: name

      select case (channel)
      case (CH_XI)
         name = "w_xi"
      case (CH_F)
         name = "w_f"
      case (CH_XYZ)
         name = "w_xyz"
      case (CH_N)
         name = "w_n"
      case (CH_K1)
         name = "w_k1"
      case (CH_K2)
         name = "w_k2"
      case (CH_A)
         name = "w_a"
      case default
         name = "w_w"
      end select
   end function channel_name

   !> Nuclear directions of the HVP assertions
   !>
   !> Neither a Cartesian axis nor a translation, so a product that dropped one
   !> atom or one axis is visible, and dense enough that every column of the
   !> block contributes to every component of the product.
   !>
   !> @param[in]  nsph Number of spheres
   !> @param[out] dirs Directions `(3, nsph, 2)`
   subroutine build_directions(nsph, dirs)
      !> Number of spheres
      integer, intent(in) :: nsph
      !> Directions
      real(wp), allocatable, intent(out) :: dirs(:, :, :)

      integer :: iatom, iaxis

      allocate (dirs(ndim, nsph, 2))
      do iatom = 1, nsph
         do iaxis = 1, ndim
            dirs(iaxis, iatom, 1) = 0.30_wp*sin(1.7_wp*real(iaxis, wp) &
                                                + 0.9_wp*real(iatom, wp))
            dirs(iaxis, iatom, 2) = 0.25_wp*cos(0.6_wp*real(iaxis, wp) &
                                                *real(iatom + 1, wp)) - 0.10_wp
         end do
      end do
   end subroutine build_directions

   !* ================================================================================= *!
   !*                              Deviation bookkeeping                                *!
   !* ================================================================================= *!

   !> Worst rank-4 component missing both bounds, folded into a running worst
   !>
   !> A component fails only when it misses the absolute *and* the relative
   !> bound: the absolute one alone would condemn a large component, the
   !> relative one alone a component that is numerically zero.
   !>
   !> @param[in]    analytic Analytic block
   !> @param[in]    ref      Numerical reference
   !> @param[in]    abs_tol  Absolute bound
   !> @param[in]    rel_tol  Relative bound
   !> @param[inout] worst    Running worst absolute deviation
   !> @param[inout] worst_rel Running worst relative deviation
   !> @param[inout] bad_ana  Analytic value at the running worst
   !> @param[inout] bad_ref  Reference value at the running worst
   !> @param[in]    iset     Set index of this block
   !> @param[in]    istep    Step index of this block
   !> @param[inout] bad_set  Set index of the running worst, zero when none
   !> @param[inout] bad_step Step index of the running worst
   subroutine worst_deviation(analytic, ref, abs_tol, rel_tol, worst, worst_rel, &
                              bad_ana, bad_ref, iset, istep, bad_set, bad_step)
      !> Analytic block
      real(wp), intent(in) :: analytic(:, :, :, :)
      !> Numerical reference
      real(wp), intent(in) :: ref(:, :, :, :)
      !> Bounds
      real(wp), intent(in) :: abs_tol, rel_tol
      !> Running worst
      real(wp), intent(inout) :: worst, worst_rel, bad_ana, bad_ref
      !> Indices of this block
      integer, intent(in) :: iset, istep
      !> Indices of the running worst
      integer, intent(inout) :: bad_set, bad_step

      integer :: iatom, jatom, iaxis, jaxis
      real(wp) :: diff, value

      do jatom = 1, size(ref, 4)
         do jaxis = 1, size(ref, 3)
            do iatom = 1, size(ref, 2)
               do iaxis = 1, size(ref, 1)
                  value = ref(iaxis, iatom, jaxis, jatom)
                  diff = abs(analytic(iaxis, iatom, jaxis, jatom) - value)
                  if (diff <= abs_tol .or. diff <= rel_tol*abs(value)) cycle
                  if (diff > worst) then
                     worst = diff
                     worst_rel = diff/max(abs(value), tiny(1.0_wp))
                     bad_ana = analytic(iaxis, iatom, jaxis, jatom)
                     bad_ref = value
                     bad_set = iset
                     bad_step = istep
                  end if
               end do
            end do
         end do
      end do
   end subroutine worst_deviation

   !> Rank-3 counterpart of [[worst_deviation]], for the HVP blocks
   !>
   !> @param[in]    analytic  Analytic block
   !> @param[in]    ref       Numerical reference
   !> @param[in]    abs_tol   Absolute bound
   !> @param[in]    rel_tol   Relative bound
   !> @param[inout] worst     Running worst absolute deviation
   !> @param[inout] worst_rel Running worst relative deviation
   !> @param[inout] bad_ana   Analytic value at the running worst
   !> @param[inout] bad_ref   Reference value at the running worst
   !> @param[in]    istep     Step index of this block
   !> @param[inout] bad_step  Step index of the running worst
   !> @param[inout] bad_set   Nonzero once anything failed
   subroutine worst_deviation_3(analytic, ref, abs_tol, rel_tol, worst, worst_rel, &
                                bad_ana, bad_ref, istep, bad_step, bad_set)
      !> Analytic block
      real(wp), intent(in) :: analytic(:, :, :)
      !> Numerical reference
      real(wp), intent(in) :: ref(:, :, :)
      !> Bounds
      real(wp), intent(in) :: abs_tol, rel_tol
      !> Running worst
      real(wp), intent(inout) :: worst, worst_rel, bad_ana, bad_ref
      !> Step index of this block
      integer, intent(in) :: istep
      !> Indices of the running worst
      integer, intent(inout) :: bad_step, bad_set

      integer :: idir, iatom, iaxis
      real(wp) :: diff, value

      do idir = 1, size(ref, 3)
         do iatom = 1, size(ref, 2)
            do iaxis = 1, size(ref, 1)
               value = ref(iaxis, iatom, idir)
               diff = abs(analytic(iaxis, iatom, idir) - value)
               if (diff <= abs_tol .or. diff <= rel_tol*abs(value)) cycle
               if (diff > worst) then
                  worst = diff
                  worst_rel = diff/max(abs(value), tiny(1.0_wp))
                  bad_ana = analytic(iaxis, iatom, idir)
                  bad_ref = value
                  bad_step = istep
                  bad_set = 1
               end if
            end do
         end do
      end do
   end subroutine worst_deviation_3

   !> Print the worst absolute and relative deviation of one block
   !>
   !> Only reached under [[E2E_VERBOSE]]; the tolerances above were measured
   !> with it.
   !>
   !> @param[in] label    Case description
   !> @param[in] name     Set name
   !> @param[in] step     Central-difference step
   !> @param[in] analytic Analytic block
   !> @param[in] ref      Numerical reference
   subroutine report_deviation(label, name, step, analytic, ref)
      !> Case description
      character(len=*), intent(in) :: label
      !> Set name
      character(len=*), intent(in) :: name
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Blocks
      real(wp), intent(in) :: analytic(:, :, :, :), ref(:, :, :, :)

      real(wp) :: worst, rel, need, diff, value
      integer :: i1, i2, i3, i4

      worst = maxval(abs(analytic - ref))
      rel = worst/max(maxval(abs(ref)), tiny(1.0_wp))
      need = 0.0_wp
      do i4 = 1, size(ref, 4)
         do i3 = 1, size(ref, 3)
            do i2 = 1, size(ref, 2)
               do i1 = 1, size(ref, 1)
                  value = ref(i1, i2, i3, i4)
                  diff = abs(analytic(i1, i2, i3, i4) - value)
                  need = max(need, min(diff, diff/max(abs(value), tiny(1.0_wp))))
               end do
            end do
         end do
      end do
      write (*, '(a,1x,a,1x,a,es10.2,3es12.3)') "DEV", label, trim(name), step, worst, &
         rel, need
   end subroutine report_deviation

   !* ================================================================================= *!
   !*                                    Fixture                                        *!
   !* ================================================================================= *!

   !> Channels whose numerical Hessian is truncation limited rather than
   !> noise limited
   !>
   !> @return Channel identifiers
   pure function smooth_channels() result(channels)
      !> Channel identifiers
      integer :: channels(6)

      channels = [CH_XI, CH_F, CH_XYZ, CH_N, CH_A, CH_W]
   end function smooth_channels

   !> Do two cavities carry the very same grid?
   !>
   !> @param[in] ref Reference cavity
   !> @param[in] cav Displaced cavity
   !> @return        Whether the point sets agree
   function same_grid(ref, cav) result(same)
      !> Reference cavity
      type(cavity_type_drop), intent(in) :: ref
      !> Displaced cavity
      type(cavity_type_drop), intent(in) :: cav
      !> Whether the point sets agree
      logical :: same

      integer :: igrid

      same = .false.
      if (cav%ngrid /= ref%ngrid) return

      do igrid = 1, ref%ngrid
         if (cav%numbering(igrid) /= ref%numbering(igrid)) return
         if (cav%owner(igrid) /= ref%owner(igrid)) return
      end do

      if (allocated(ref%branch_count) .and. allocated(cav%branch_count)) then
         do igrid = 1, ref%ngrid
            if (cav%branch_count(igrid) /= ref%branch_count(igrid)) return
         end do
      end if

      same = .true.
   end function same_grid

   !> Populate the requested adjoint channels
   !>
   !> Every weight is a pure function of the *persistent* point id
   !> `cavity%numbering`, so the same adjoint is reproduced on a displaced grid
   !> with no explicit mapping -- which is what "the adjoints are held fixed"
   !> has to mean once the grid is filtered and reordered on every rebuild.
   !>
   !> @param[in]  cavity Cavity supplying the grid
   !> @param[in]  mask   Channels to populate
   !> @param[out] acc    Surface-adjoint accumulator
   !> @param[out] error  Error handle
   subroutine frozen_adjoint(cavity, mask, acc, error)
      !> Cavity supplying the grid
      type(cavity_type_drop), intent(in) :: cavity
      !> Channels to populate
      logical, intent(in) :: mask(:)
      !> Surface-adjoint accumulator
      type(cavity_surface_adjoint_type), intent(out) :: acc
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: add_error
      real(wp), allocatable :: ws(:), wv(:, :)
      integer :: ngrid, igrid, iaxis, ichannel

      ngrid = cavity%ngrid
      allocate (ws(ngrid), wv(ndim, ngrid))

      call acc%init(ngrid)
      do ichannel = 1, NCHAN
         if (.not. mask(ichannel)) cycle

         do igrid = 1, ngrid
            ws(igrid) = point_weight(cavity%numbering(igrid), ichannel)
            do iaxis = 1, ndim
               wv(iaxis, igrid) = point_weight(cavity%numbering(igrid), &
                                               NCHAN*ichannel + iaxis)
            end do
         end do

         select case (ichannel)
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
         case (CH_W)
            call acc%add_surface_weights(add_error, w_w=ws)
         end select
         if (allocated(add_error)) then
            call test_failed(error, "failed to seed the surface adjoint: "//add_error%message)
            return
         end if
      end do

   end subroutine frozen_adjoint

   !> Reproducible adjoint weight of one persistent point and channel
   !>
   !> Identical to the fixed-adjoint suite's, so the two sets of numbers
   !> describe the same functional. Varies across the grid so a term that
   !> cancels for uniform weights still shows, and is bounded away from zero so
   !> no channel is accidentally switched off.
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

   !> Shared fixture geometry
   !>
   !> Asymmetric on purpose: a symmetric geometry drives the multistart
   !> projection into sibling branches, and a branched grid moves the branch
   !> softmax with the geometry
   !>
   !> @param[out] mol Structure
   subroutine fixture_geometry(mol)
      !> Structure
      type(structure_type), intent(out) :: mol

      call new(mol, [8, 6, 1], reshape([ &
                                       0.00_wp, 0.00_wp, 0.00_wp, &
                                       0.00_wp, 0.00_wp, 4.60_wp, &
                                       2.60_wp, 0.40_wp, -1.10_wp], [3, 3]))
   end subroutine fixture_geometry

   !> Build the DROP cavity for a structure and a level-set model
   !>
   !> @param[out]   cavity   Constructed cavity
   !> @param[inout] ctx      Run context borrowed by the cavity; must outlive it
   !> @param[in]    mol      Structure to build on
   !> @param[in]    lsf_kind Level-set model
   !> @param[out]   error    Error handle
   subroutine build_cavity(cavity, ctx, mol, lsf_kind, error)
      !> Constructed cavity
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(inout) :: ctx
      !> Structure to build on
      type(structure_type), intent(in) :: mol
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: radii(:)
      type(mctc_error), allocatable :: cav_error

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      allocate (cavity)
      call new_context(ctx, verbosity=0)
      select case (lsf_kind)
      case (LSF_SVDW)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=BLEND_K, blend_3b=BLEND_3B)
            call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                                 tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                                 proj_level=PROJ_LEVEL, wleb_prune_level=WLEB_PRUNE, &
                                 radius_model=default_cpcm_radii(), &
                                 lsf_model=svdw_template, error=cav_error)
         end block
      case default
         block
            type(moist_cavity_drop_lsf_cfc_type) :: cfc_template
            call cfc_template%new()
            call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                                 tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                                 proj_level=PROJ_LEVEL, wleb_prune_level=WLEB_PRUNE, &
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

end module test_cavity_drop_hessian_e2e
