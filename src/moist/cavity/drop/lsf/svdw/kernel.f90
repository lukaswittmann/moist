!> Low-level Smooth van der Waals (SvdW) level-set kernel
!>
!> Model
!> -----
!> - atom A at R_A, radius R_a: d_A = ||r - R_A|| - R_a, u_A = exp(-(k/3) d_A), g_A = w(d_A) u_A,
!>   k = blend_k, fade w as below
!> - power sums p_m = sum_A g_A**m, q = sum_A g_A**(3/2)
!>
!>    Z = s_1 p3 + (s_2/2)(q**2 - p3) + (s_3/6)(p1**3 - 3 p1 p2 + 2 p3) ,   S = -(1/k) log Z
!>
!> - one-, two- and three-body terms in power-sum form (Newton's identities): O(n), never
!>   over pairs or triples
!> - kind slot: 1 -> m = 1 (p1), 2 -> m = 2 (p2), 3 -> m = 3 (p3), 4 -> m = 3/2 (q)
!>
!> Tensor families
!> ---------------
!> Each power sum is a plain sum of single-atom functions, P_k = sum_A g_A**m_k:
!>
!>    ps{n}(.., kind)  d^n/dr^n P_k          aggregate over atoms
!>    at{n}(.., kind)  d^n/dr^n g_A**m_k     one atom A
!>    bt{n}(.., kind)  d^n/dr^n g_B**m_k     one atom B
!>    ws{n}(.., kind)  sum_B v_B . d/dR_B of ps{n}
!>    aw{n}(.., kind)  v_A . d/dR_A of at{n}
!>    atd{n}, atdd{n}  d/dR_a and d^2/dR_a^2 of at{n}
!>    btd{n}           d/dR_b of bt{n}
!>    awd{n}           d/dR_a of aw{n}
!>
!> Exact rules (r field point, R_A nucleus A, R_a its radius):
!>
!>    d/dr_i    : F{n} -> F{n+1}, one more index          (all families)
!>    d/dR_A,b  : ps -> -at,  ws -> -aw,  atom-A families F -> -F{n+1},  B -> 0
!>    d/dR_B,t  : ps -> -bt,  bt -> -bt,  btd -> -btd,  atom-A families -> 0
!>    d/dR_a    : ps -> atd,  at -> atd,  atd -> atdd,  ws -> awd,  B -> 0
!>    d/dR_b    : ps -> btd,  bt -> btd,  atom-A families -> 0
!>    sum_B v_B . d/dR_B : ps -> ws,  at -> aw
!>
!> - radius rules bump no index: R_a enters only through d_A, so d/dR_a is the same spatial
!>   tensor of -dg/dd; hence `f4_rr_radrad` two ranks below `f4_rr_rArB`
!> - band at infinity: -dg/dd = (k m/3) g
!> - no cross-atom second derivative of a power sum: all A-B coupling in `f2_rArB`, `f2_radrad`
!>   and `f2_rA_rad` comes from Z and the log
!>
!> Routines and cost
!> -----------------
!> - `svdw_strip_powersums_o{L}`: power sums plus geometry cache (x, 1/x, u)
!> - `svdw_atom_*_eval`: per-atom tensors from that cache, O(1) per active atom
!> - `svdw_powersums`: fused libm reference at the band at infinity, for the reference drivers
!> - O(1) per point: `svdw_spatial_eval`, `svdw_tangent_eval`
!> - O(1) per active atom: `svdw_nuclear_eval`, `svdw_hvp_eval`, `svdw_normalized_eval`
!> - once per atom pair, result O(n_active**2): `svdw_pair_eval`, `svdw_pair_diag_eval`,
!>   `svdw_radpair_eval`, `svdw_nucrad_eval`; the `hvp_*` outputs avoid forming it
!> - `svdw_vjp_eval`: reverse mode of `svdw_tangent_eval`, spatial jet indices contracted with
!>   one point's adjoint weights before the nuclear derivative; 3 numbers per atom instead of
!>   3 + 9 + 27
!> - `svdw_hvp_vjp_eval`: same contraction of `svdw_hvp_eval`, direction left open as a linear
!>   form in `ws`; weighted mixed nuclear block for every direction from one call per atom
!>   plus a matrix product
!> - `svdw_vjp_diag_eval`: weighted own-atom nuclear Hessian block; with `svdw_hvp_vjp_eval`
!>   the weighted hvp row for any direction from per-atom quantities alone
!> - derivative stacks dispatch on `max_deriv`: out-of-range requests take the highest branch,
!>   results above `max_deriv` left untouched, never zeroed
!> - `svdw_vjp_eval`, `svdw_hvp_vjp_eval`, `svdw_vjp_diag_eval`: one order, no `max_deriv`
!>
!> Strip kernel
!> ------------
!> - `svdw_strip_powersums_o{L}` (stage 1) and `svdw_strip_lift_o{L}` (stage 2), L = 0..4:
!>   spatial jet of S through order L for `np <= svdw_strip_width(L)` points sharing one
!>   atom list; widths 32, 32, 32, 16, 8
!> - atom loop outer, point loop innermost, no branch inside a point loop: plateau, faded and
!>   skipped atoms and the `include_q` split are separate loops
!> - per-point path: a strip of one
!> - on-nucleus pairs counted (`n_zero`), not guarded; caller falls back to the guarded
!>   per-atom path for such a strip
!>
!> Layouts
!> - points `(3, np)`
!> - stage-1 power sums point-first: `ps0(np, nkind)`, `ps1(np, 3, nkind)`,
!>   `ps{n}(np, svdw_npack(n), nkind)`
!> - lift outputs point-last: `lsf0(np)`, `lsf1_r(3, np)`, `lsf2_rr(6, np)`,
!>   `lsf3_rrr(10, np)`, `lsf4_rrrr(15, np)`
!> - symmetric groups packed as sorted tuples a <= b <= c <= d, last index slowest, as
!>   `moist_math_packed`:
!>
!>    slot = a + b(b-1)/2 + (c-1)c(c+1)/6 + (d-1)d(d+1)(d+2)/24
!>
!> Fade
!> - `w = 1 - S_n(t) = S_n(s)`, n = 4, `t = clamp((d_A - fade_start)/fade_width, 0, 1)`,
!>   `s = 1 - t`
!> - evaluated as `s**(n+1) P_n`, `P_n = sum_k C(n+k, k) t**k`: positive coefficients, exact
!>   zero beyond the band
!> - band at infinity: `fade_start = huge(0.0_wp)`, finite `fade_width`, every atom on the
!>   plain body
!>
!> Inline exp
!> - `exp(x) = 2**k (1 + r + r**2 q(r))`, `k = anint(x/ln 2)`, two-constant Cody-Waite
!>   reduction to |r| <= ln2/2, q of degree 10, `2**k` from its exponent bits
!> - argument clamped below at -708; `blend_k*max(R_a)/3 <= 708` checked per strip
!> - <= 1 ulp, from the generator's `--check` sweep with and without FMA
!> - one exp serves the sums and the per-atom tensors (geometry cache); `svdw_powersums` and
!>   the reference drivers keep libm `exp`
!>
!> Batched nuclear vjp
!> -------------------
!> - `svdw_strip_powersums_o2` (with geometry cache) -> `svdw_strip_lift_vjp` ->
!>   `svdw_strip_vjp`: `res(:, a, p) = d/dR_a [w0 S + w1 . dS/dr + W : d2S/dr2]` for every
!>   atom of the shared list, `np <= 32` points
!> - seeds `w0(np)`, `w1(3, np)`, `w2(6, np)`: symmetric part of the (3, 3) seed packed xx, xy,
!>   yy, xz, yz, zz, off-diagonals not doubled (the kernel applies the factor 2)
!> - stage 2' emits the adjoint `lam` of `ps0`, `ps1`, `ps2` in their layout; stage 3 contracts
!>   it with every atom's radial chain against the stage-1 cache, never forming `at3`
!> - `res(3, nat, np)` point-last; exact zeros beyond the band
!>
!> Source layout
!> -------------
!> - this file: declarations, public names, interfaces
!> - `kernel/<group>.f90`: one submodule per group, compiled in parallel: atom, atan, spat,
!>   pair, nucr, tang, vjp, strip, orch (also the reference drivers)
module moist_cavity_drop_lsf_svdw_kernel
   use mctc_env_accuracy, only: wp
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none (type, external)
   private

   !> Cartesian dimension.
   integer, parameter :: ndim = 3

   !> Number of exponential kinds carried by every tensor family:
   !> 1 = u**1 (p1), 2 = u**2 (p2), 3 = u**3 (p3), 4 = u**(3/2) (q).
   integer, parameter, public :: nkind = 4

   !> Band at infinity of the reference drivers: fade onset, and any finite width
   real(wp), parameter :: ref_fade_start = huge(1.0_wp), ref_fade_width = 1.0_wp

   public :: svdw_atom_eval
   public :: svdw_atom_tangent_eval
   public :: svdw_atom_radius_tangent_eval
   public :: svdw_atom_rad_eval
   public :: svdw_atom_radrad_eval
   public :: svdw_atom_tangent_rad_eval
   public :: svdw_spatial_eval
   public :: svdw_nuclear_eval
   public :: svdw_radius_eval
   public :: svdw_pair_eval
   public :: svdw_pair_diag_eval
   public :: svdw_radpair_eval
   public :: svdw_radpair_diag_eval
   public :: svdw_nucrad_eval
   public :: svdw_nucrad_diag_eval
   public :: svdw_tangent_eval
   public :: svdw_hvp_eval
   public :: svdw_radius_hvp_eval
   public :: svdw_vjp_eval
   public :: svdw_radius_vjp_eval
   public :: svdw_hvp_vjp_eval
   public :: svdw_vjp_diag_eval
   public :: svdw_normalized_eval

   public :: svdw_powersums
   public :: svdw_atom_geometry
   public :: svdw_powersums_tangent
   public :: svdw_lsf_spatial
   public :: svdw_lsf_nuclear
   public :: svdw_lsf_radius
   public :: svdw_lsf_pair
   public :: svdw_lsf_tangent
   public :: svdw_lsf_hvp
   public :: svdw_lsf_radius_hvp
   public :: svdw_lsf_normalized

   !> Points per strip of the strip kernels, by derivative order 0..4
   integer, parameter, public :: svdw_strip_width(0:4) = [32, 32, 32, 16, 8]
   !> Packed components of a symmetric rank-n spatial tensor, n = 0..4
   integer, parameter, public :: svdw_npack(0:4) = [1, 3, 6, 10, 15]

   !> Inline exp: 1/ln 2 and the Cody-Waite split of ln 2, exact bit patterns
   real(wp), parameter :: exp_inv_ln2 = transfer(int(z'3FF71547652B82FE', int64), 1.0_wp)
   real(wp), parameter :: exp_ln2_hi = transfer(int(z'3FE62E42FEE00000', int64), 1.0_wp)
   real(wp), parameter :: exp_ln2_lo = transfer(int(z'3DEA39EF35793C76', int64), 1.0_wp)
   !> Inline exp: q(r) of exp(r) ~ 1 + r + r**2 q(r), ascending, exact bit patterns
   real(wp), parameter :: exp_q0 = transfer(int(z'3FE0000000000000', int64), 1.0_wp)
   real(wp), parameter :: exp_q1 = transfer(int(z'3FC5555555555557', int64), 1.0_wp)
   real(wp), parameter :: exp_q2 = transfer(int(z'3FA5555555555556', int64), 1.0_wp)
   real(wp), parameter :: exp_q3 = transfer(int(z'3F811111111100DF', int64), 1.0_wp)
   real(wp), parameter :: exp_q4 = transfer(int(z'3F56C16C16C162D6', int64), 1.0_wp)
   real(wp), parameter :: exp_q5 = transfer(int(z'3F2A01A01ABE62DD', int64), 1.0_wp)
   real(wp), parameter :: exp_q6 = transfer(int(z'3EFA01A01A6D7808', int64), 1.0_wp)
   real(wp), parameter :: exp_q7 = transfer(int(z'3EC71DE02375656C', int64), 1.0_wp)
   real(wp), parameter :: exp_q8 = transfer(int(z'3E927E4DB67B4303', int64), 1.0_wp)
   real(wp), parameter :: exp_q9 = transfer(int(z'3E5AF4DDD84882FE', int64), 1.0_wp)
   real(wp), parameter :: exp_q10 = transfer(int(z'3E21F72FC730B4FF', int64), 1.0_wp)
   !> Inline exp argument range; 2**k stays normal inside it
   real(wp), parameter :: exp_arg_min = -708.0_wp
   real(wp), parameter :: exp_arg_max = 708.0_wp

   public :: svdw_strip_powersums_o0
   public :: svdw_strip_lift_o0
   public :: svdw_strip_powersums_o1
   public :: svdw_strip_lift_o1
   public :: svdw_strip_powersums_o2
   public :: svdw_strip_lift_o2
   public :: svdw_strip_powersums_o3
   public :: svdw_strip_lift_o3
   public :: svdw_strip_powersums_o4
   public :: svdw_strip_lift_o4

   public :: svdw_strip_lift_vjp
   public :: svdw_strip_vjp

   !* ================================================================================= *!
   !*                                     Interfaces                                    *!
   !* ================================================================================= *!

   interface

   !> Accumulate the power-sum tensors of one evaluation point.
   !>
   !> Sets `ps{n}(i1..in, kind)` = sum over active atoms of
   !> `d^n/dr_i1..dr_in of u_A**m_kind`, for n = 0..max_deriv.  This is the
   !> single hottest loop in the DROP projection, so it is emitted *fused*:
   !>
   !>   * nothing per-atom is materialised.  The unfused spelling -- call
   !>     `svdw_atom_eval` into `at0..at4`, then add the buffers into `ps0..ps4`
   !>     -- writes and re-reads 52 doubles per atom at `max_deriv = 2` alone,
   !>     and `svdw_atom_eval` is far too large for gfortran to inline, so the
   !>     buffers cannot be optimised away.  Here every value goes straight from
   !>     a register into its accumulator.
   !>   * only the *independent* component of each symmetric block is
   !>     accumulated (6 of 9 at order 2, 10 of 27 at order 3, 15 of 81 at order
   !>     4).  The symmetric copies are broadcast once, after the loop, so the
   !>     caller still sees the full dense tensors and no consumer of `ps*` has
   !>     to know about the reduction.
   !>   * the four kinds share their geometric factors: the `n_i n_j` products
   !>     and the powers of `1/x` are common subexpressions of one atom, hoisted
   !>     by cse across all four exponents at once.
   !>   * the kind-4 (`q`) work sits behind `include_q`. Z reaches that kind only
   !>     through the two-body weight, and the library default sets that weight to
   !>     zero, so on the default path a quarter of the accumulators -- and the
   !>     only `sqrt` besides the distance -- drop out. The guard is loop-invariant
   !>     and perfectly predicted; the shared part is emitted once, so this costs
   !>     no duplicated code.
   !>
   !> `max_deriv` selects the branch, so the value, gradient and Hessian paths
   !> never touch third- or fourth-order work.  Orders above `max_deriv` are left
   !> untouched, never zeroed.
   !>
   !> Requires `||d(:, ia)|| > 0` for every active atom: the derivatives of the
   !> distance are singular on a nucleus.  A caller that cannot rule that out
   !> must fall back to a guarded loop over `svdw_atom_eval`.
   !>
   !> @param[in]    d            Displacements r - R_A of the active atoms [3, n_active]
   !> @param[in]    radii        Radii of the active atoms [n_active]
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    include_q    Accumulate the kind-4 (u**(3/2), `q`) power sums? Pass `blend_2b /=
   !>                            0`. When false the kind-4 slots are SET TO ZERO rather than
   !>                            computed, which is sound only because Z depends on that kind solely
   !>                            through the two-body weight, so every consumer multiplies them by
   !>                            that weight
   !> @param[in]    max_deriv    Highest spatial-derivative order to accumulate (0..4); a larger
   !>                            value is treated as 4
   !> @param[out]   ps0          Order-0 power-sum tensors
   !> @param[out]   ps1          Order-1 power-sum tensors
   !> @param[out]   ps2          Order-2 power-sum tensors
   !> @param[out]   ps3          Order-3 power-sum tensors
   !> @param[out]   ps4          Order-4 power-sum tensors
   pure module subroutine svdw_powersums(d, radii, blend_k, include_q, max_deriv, ps0, ps1, ps2, ps3, ps4)
      !> Displacements r - R_A of the active atoms [3, n_active]
      real(wp), contiguous, intent(in) :: d(:, :)
      !> Radii of the active atoms [n_active]
      real(wp), contiguous, intent(in) :: radii(:)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Accumulate the kind-4 (u**(3/2), `q`) power sums? Pass `blend_2b /= 0`. When false the
      !> kind-4 slots are SET TO ZERO rather than computed, which is sound only because Z depends on
      !> that kind solely through the two-body weight, so every consumer multiplies them by that
      !> weight
      logical, intent(in) :: include_q
      !> Highest spatial-derivative order to accumulate (0..4); a larger value is treated as 4
      integer, intent(in) :: max_deriv
      !> Order-0 power-sum tensors
      real(wp), intent(out) :: ps0(nkind)
      !> Order-1 power-sum tensors
      real(wp), intent(out) :: ps1(3, nkind)
      !> Order-2 power-sum tensors
      real(wp), intent(out) :: ps2(3, 3, nkind)
      !> Order-3 power-sum tensors
      real(wp), intent(out) :: ps3(3, 3, 3, nkind)
      !> Order-4 power-sum tensors
      real(wp), intent(out) :: ps4(3, 3, 3, 3, nkind)
   end subroutine svdw_powersums

   !> Per-atom kind tensors of one active atom.
   !>
   !> Sets `at{n}(i1..in, kind)` = d^n/dr_i1..dr_in of the atom function
   !> `g = (w(d_A) u_A)**m`, `d_A = ||d_a|| - R_a`, `u_A = exp(-(blend_k/3) d_A)`,
   !> where the kind slot selects the exponent: 1 -> m = 1, 2 -> m = 2,
   !> 3 -> m = 3, 4 -> m = 3/2.  The fade `w = 1 - S_n(t)`,
   !> `t = clamp((d_A - fade_start)/fade_width, 0, 1)`, is the strip kernel's.
   !>
   !> The distance, its reciprocal and `u_A` are read from the stage-1 geometry
   !> cache of `svdw_strip_powersums_o{L}` (`cache_x`, `cache_inv_x`, `cache_u` of
   !> this point and atom), so the per-atom tensors and the strip power sums see
   !> one exp per atom and point.  With `d_A <= fade_start` (always, for the band
   !> at infinity `fade_start = huge(0.0_wp)`) the plain atom function runs, else
   !> the faded one; beyond the band every result is an exact zero.
   !>
   !> The caller runs this once per active atom; the per-atom tensors are reused
   !> by `svdw_nuclear_eval`, `svdw_pair_eval` and `svdw_hvp_eval`.
   !>
   !> Nuclear derivatives never need a separate table: `g` depends on `r` and
   !> `R_A` only through `d_a = r - R_A`, hence `d/dR_A,b = -d/dr_b` and the
   !> nuclear index is just one more entry of the same tensor, with a sign.
   !>
   !> The plateau branch divides by `x_a`, the arithmetic of the unfaded kernel,
   !> so the band at infinity reproduces it bit for bit on the same `x` and `u`;
   !> only the faded branch reads `inv_x_a`, at orders >= 1.  Order 0 reads only
   !> `u_a`, so on a nucleus (`x_a = 0`, `inv_x_a` Inf in the strip cache) only
   !> order 0 is defined.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    d_a          Displacement r - R_A of the active atom
   !> @param[in]    R_a          Radius of the active atom
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    max_deriv    Highest spatial-derivative order to set (0..4); a larger value is
   !>                            treated as 4
   !> @param[out]   at0          Order-0 kind tensor of this atom
   !> @param[out]   at1          Order-1 kind tensor of this atom
   !> @param[out]   at2          Order-2 kind tensor of this atom
   !> @param[out]   at3          Order-3 kind tensor of this atom
   !> @param[out]   at4          Order-4 kind tensor of this atom
   !> @param[in]    x_a          Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and
   !>                            atom
   !> @param[in]    inv_x_a      Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded
   !>                            branch at orders >= 1 only (Inf on a nucleus)
   !> @param[in]    u_a          Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade,
   !>                            gate-masked: stage-1 `cache_u` (inline exp), 0 where the strip's
   !>                            hard screening dropped the atom
   !> @param[in]    fade_start   Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp)
   !>                            for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, as passed to the strip kernel; > 0
   pure module subroutine svdw_atom_eval(d_a, R_a, blend_k, max_deriv, at0, at1, at2, at3, at4, x_a, &
                                  inv_x_a, u_a, fade_start, fade_width)
      !> Displacement r - R_A of the active atom
      real(wp), intent(in) :: d_a(3)
      !> Radius of the active atom
      real(wp), intent(in) :: R_a
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Highest spatial-derivative order to set (0..4); a larger value is treated as 4
      integer, intent(in) :: max_deriv
      !> Order-0 kind tensor of this atom
      real(wp), intent(out) :: at0(nkind)
      !> Order-1 kind tensor of this atom
      real(wp), intent(out) :: at1(3, nkind)
      !> Order-2 kind tensor of this atom
      real(wp), intent(out) :: at2(3, 3, nkind)
      !> Order-3 kind tensor of this atom
      real(wp), intent(out) :: at3(3, 3, 3, nkind)
      !> Order-4 kind tensor of this atom
      real(wp), intent(out) :: at4(3, 3, 3, 3, nkind)
      !> Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and atom
      real(wp), intent(in) :: x_a
      !> Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded branch at orders >= 1
      !> only (Inf on a nucleus)
      real(wp), intent(in) :: inv_x_a
      !> Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade, gate-masked: stage-1
      !> `cache_u` (inline exp), 0 where the strip's hard screening dropped the atom
      real(wp), intent(in) :: u_a
      !> Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, as passed to the strip kernel; > 0
      real(wp), intent(in) :: fade_width
   end subroutine svdw_atom_eval

   !> Per-atom direction-contracted kind tensors of one active atom.
   !>
   !> Sets `aw{n}(i1..in, kind)` = `v_a . d/dR_A` of `at{n}(i1..in, kind)`
   !> = `-sum_alpha v_a(alpha) * at{n+1}(i1..in, alpha, kind)`, on the atom function
   !> and geometry cache of `svdw_atom_eval`.
   !>
   !> Contracting the nuclear index here, against the caller's per-atom
   !> displacement direction, is what keeps every `tangent_*` and `hvp_*` output
   !> O(n_active): the caller accumulates `ws* = ws* + aw*` in the same pass that
   !> builds `ps*`, and the full pairwise nuclear tensor is never formed.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    d_a          Displacement r - R_A of the active atom
   !> @param[in]    R_a          Radius of the active atom
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    v_a          Nuclear displacement direction of this atom
   !> @param[in]    max_deriv    Highest spatial-derivative order to set (0..3); a larger value is
   !>                            treated as 3
   !> @param[out]   aw0          Order-0 contracted kind tensor of this atom
   !> @param[out]   aw1          Order-1 contracted kind tensor of this atom
   !> @param[out]   aw2          Order-2 contracted kind tensor of this atom
   !> @param[out]   aw3          Order-3 contracted kind tensor of this atom
   !> @param[in]    x_a          Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and
   !>                            atom
   !> @param[in]    inv_x_a      Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded
   !>                            branch at orders >= 1 only (Inf on a nucleus)
   !> @param[in]    u_a          Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade,
   !>                            gate-masked: stage-1 `cache_u` (inline exp), 0 where the strip's
   !>                            hard screening dropped the atom
   !> @param[in]    fade_start   Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp)
   !>                            for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, as passed to the strip kernel; > 0
   pure module subroutine svdw_atom_tangent_eval(d_a, R_a, blend_k, v_a, max_deriv, aw0, aw1, aw2, aw3, &
                                          x_a, inv_x_a, u_a, fade_start, fade_width)
      !> Displacement r - R_A of the active atom
      real(wp), intent(in) :: d_a(3)
      !> Radius of the active atom
      real(wp), intent(in) :: R_a
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Nuclear displacement direction of this atom
      real(wp), intent(in) :: v_a(3)
      !> Highest spatial-derivative order to set (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> Order-0 contracted kind tensor of this atom
      real(wp), intent(out) :: aw0(nkind)
      !> Order-1 contracted kind tensor of this atom
      real(wp), intent(out) :: aw1(3, nkind)
      !> Order-2 contracted kind tensor of this atom
      real(wp), intent(out) :: aw2(3, 3, nkind)
      !> Order-3 contracted kind tensor of this atom
      real(wp), intent(out) :: aw3(3, 3, 3, nkind)
      !> Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and atom
      real(wp), intent(in) :: x_a
      !> Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded branch at orders >= 1
      !> only (Inf on a nucleus)
      real(wp), intent(in) :: inv_x_a
      !> Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade, gate-masked: stage-1
      !> `cache_u` (inline exp), 0 where the strip's hard screening dropped the atom
      real(wp), intent(in) :: u_a
      !> Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, as passed to the strip kernel; > 0
      real(wp), intent(in) :: fade_width
   end subroutine svdw_atom_tangent_eval

   !> Per-atom radius-contracted kind tensors of one active atom.
   !>
   !> Sets `awr{n}(i1..in, kind)` = `vr_a * d/dR_a` of `at{n}(i1..in, kind)`
   !> = `vr_a * atd{n}(i1..in, kind)`, where `R_a` is the *radius* of atom A (see
   !> `svdw_atom_rad_eval`); on the plateau that is
   !> `(blend_k * m_kind / 3) * vr_a * at{n}(i1..in, kind)`.
   !>
   !> This is the radius half of the joint direction `(v_B, vr_B)` over nuclear
   !> positions and radii.  It lands in the same family as
   !> `svdw_atom_tangent_eval`, so a caller that needs the radius channel simply
   !> accumulates `ws* = ws* + aw* + awr*` in the same pass and every downstream
   !> `tangent_*` and `hvp_*` output picks the radius term up unchanged -- no
   !> aggregate routine in this module knows the channel exists.
   !>
   !> It is a separate routine rather than an extra term inside
   !> `svdw_atom_tangent_eval` so that a geometry-*independent* radius model,
   !> which simply does not call it, keeps the nuclear contraction's arithmetic:
   !> folding the two halves into one cse pool perturbs the heavily cancelling
   !> nuclear contraction by ~1e-12 relative even when `vr_a` is exactly zero.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    d_a          Displacement r - R_A of the active atom
   !> @param[in]    R_a          Radius of the active atom
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    vr_a         Radius direction of this atom, vr_A = sum_B (dR_A/dX_B) . v_B
   !> @param[in]    max_deriv    Highest spatial-derivative order to set (0..3); a larger value is
   !>                            treated as 3
   !> @param[out]   awr0         Order-0 radius-contracted kind tensor of this atom
   !> @param[out]   awr1         Order-1 radius-contracted kind tensor of this atom
   !> @param[out]   awr2         Order-2 radius-contracted kind tensor of this atom
   !> @param[out]   awr3         Order-3 radius-contracted kind tensor of this atom
   !> @param[in]    x_a          Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and
   !>                            atom
   !> @param[in]    inv_x_a      Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded
   !>                            branch at orders >= 1 only (Inf on a nucleus)
   !> @param[in]    u_a          Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade,
   !>                            gate-masked: stage-1 `cache_u` (inline exp), 0 where the strip's
   !>                            hard screening dropped the atom
   !> @param[in]    fade_start   Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp)
   !>                            for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, as passed to the strip kernel; > 0
   pure module subroutine svdw_atom_radius_tangent_eval(d_a, R_a, blend_k, vr_a, max_deriv, awr0, awr1, &
                                                 awr2, awr3, x_a, inv_x_a, u_a, fade_start, &
                                                 fade_width)
      !> Displacement r - R_A of the active atom
      real(wp), intent(in) :: d_a(3)
      !> Radius of the active atom
      real(wp), intent(in) :: R_a
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Radius direction of this atom, vr_A = sum_B (dR_A/dX_B) . v_B
      real(wp), intent(in) :: vr_a
      !> Highest spatial-derivative order to set (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> Order-0 radius-contracted kind tensor of this atom
      real(wp), intent(out) :: awr0(nkind)
      !> Order-1 radius-contracted kind tensor of this atom
      real(wp), intent(out) :: awr1(3, nkind)
      !> Order-2 radius-contracted kind tensor of this atom
      real(wp), intent(out) :: awr2(3, 3, nkind)
      !> Order-3 radius-contracted kind tensor of this atom
      real(wp), intent(out) :: awr3(3, 3, 3, nkind)
      !> Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and atom
      real(wp), intent(in) :: x_a
      !> Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded branch at orders >= 1
      !> only (Inf on a nucleus)
      real(wp), intent(in) :: inv_x_a
      !> Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade, gate-masked: stage-1
      !> `cache_u` (inline exp), 0 where the strip's hard screening dropped the atom
      real(wp), intent(in) :: u_a
      !> Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, as passed to the strip kernel; > 0
      real(wp), intent(in) :: fade_width
   end subroutine svdw_atom_radius_tangent_eval

   !> Per-atom radius tensors of one active atom.
   !>
   !> Sets `atd{n}(i1..in, kind)` = `d/dR_a` of `at{n}(i1..in, kind)`, where `R_a`
   !> is the *radius* of atom A.  The atom function depends on `R_a` only through
   !> `d_A = ||d_a|| - R_a`, so this is the same spatial derivative of `-dg/dd`.
   !> On the plateau `-dg/dd = (blend_k m_kind / 3) g` and `atd` is that multiple
   !> of `at`; inside the band the fade breaks the proportionality.
   !>
   !> Same geometry cache, fade and plateau/faded branch as `svdw_atom_eval`.
   !> Read by `svdw_radius_eval`, `svdw_radius_vjp_eval`, `svdw_radpair*_eval`
   !> and `svdw_nucrad_diag_eval` (as `atd`, and as `btd` for the second atom).
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    d_a          Displacement r - R_A of the active atom
   !> @param[in]    R_a          Radius of the active atom
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    max_deriv    Highest spatial-derivative order to set (0..3); a larger value is
   !>                            treated as 3
   !> @param[out]   atd0         Order-0 radius kind tensor of this atom
   !> @param[out]   atd1         Order-1 radius kind tensor of this atom
   !> @param[out]   atd2         Order-2 radius kind tensor of this atom
   !> @param[out]   atd3         Order-3 radius kind tensor of this atom
   !> @param[in]    x_a          Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and
   !>                            atom
   !> @param[in]    inv_x_a      Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded
   !>                            branch at orders >= 1 only (Inf on a nucleus)
   !> @param[in]    u_a          Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade,
   !>                            gate-masked: stage-1 `cache_u` (inline exp), 0 where the strip's
   !>                            hard screening dropped the atom
   !> @param[in]    fade_start   Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp)
   !>                            for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, as passed to the strip kernel; > 0
   pure module subroutine svdw_atom_rad_eval(d_a, R_a, blend_k, max_deriv, atd0, atd1, atd2, atd3, x_a, &
                                      inv_x_a, u_a, fade_start, fade_width)
      !> Displacement r - R_A of the active atom
      real(wp), intent(in) :: d_a(3)
      !> Radius of the active atom
      real(wp), intent(in) :: R_a
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Highest spatial-derivative order to set (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> Order-0 radius kind tensor of this atom
      real(wp), intent(out) :: atd0(nkind)
      !> Order-1 radius kind tensor of this atom
      real(wp), intent(out) :: atd1(3, nkind)
      !> Order-2 radius kind tensor of this atom
      real(wp), intent(out) :: atd2(3, 3, nkind)
      !> Order-3 radius kind tensor of this atom
      real(wp), intent(out) :: atd3(3, 3, 3, nkind)
      !> Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and atom
      real(wp), intent(in) :: x_a
      !> Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded branch at orders >= 1
      !> only (Inf on a nucleus)
      real(wp), intent(in) :: inv_x_a
      !> Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade, gate-masked: stage-1
      !> `cache_u` (inline exp), 0 where the strip's hard screening dropped the atom
      real(wp), intent(in) :: u_a
      !> Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, as passed to the strip kernel; > 0
      real(wp), intent(in) :: fade_width
   end subroutine svdw_atom_rad_eval

   !> Per-atom second radius tensors of one active atom.
   !>
   !> Sets `atdd{n}(i1..in, kind)` = `d^2/dR_a^2` of `at{n}(i1..in, kind)`, the same
   !> spatial derivative of `d^2g/dd^2`; `(blend_k m_kind / 3)**2 at{n}` on the
   !> plateau.  Same geometry cache, fade and branch as `svdw_atom_eval`.  Read by
   !> `svdw_radpair_diag_eval`.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    d_a          Displacement r - R_A of the active atom
   !> @param[in]    R_a          Radius of the active atom
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    max_deriv    Highest spatial-derivative order to set (0..2); a larger value is
   !>                            treated as 2
   !> @param[out]   atdd0        Order-0 second radius kind tensor of this atom
   !> @param[out]   atdd1        Order-1 second radius kind tensor of this atom
   !> @param[out]   atdd2        Order-2 second radius kind tensor of this atom
   !> @param[in]    x_a          Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and
   !>                            atom
   !> @param[in]    inv_x_a      Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded
   !>                            branch at orders >= 1 only (Inf on a nucleus)
   !> @param[in]    u_a          Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade,
   !>                            gate-masked: stage-1 `cache_u` (inline exp), 0 where the strip's
   !>                            hard screening dropped the atom
   !> @param[in]    fade_start   Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp)
   !>                            for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, as passed to the strip kernel; > 0
   pure module subroutine svdw_atom_radrad_eval(d_a, R_a, blend_k, max_deriv, atdd0, atdd1, atdd2, x_a, &
                                         inv_x_a, u_a, fade_start, fade_width)
      !> Displacement r - R_A of the active atom
      real(wp), intent(in) :: d_a(3)
      !> Radius of the active atom
      real(wp), intent(in) :: R_a
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Highest spatial-derivative order to set (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> Order-0 second radius kind tensor of this atom
      real(wp), intent(out) :: atdd0(nkind)
      !> Order-1 second radius kind tensor of this atom
      real(wp), intent(out) :: atdd1(3, nkind)
      !> Order-2 second radius kind tensor of this atom
      real(wp), intent(out) :: atdd2(3, 3, nkind)
      !> Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and atom
      real(wp), intent(in) :: x_a
      !> Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded branch at orders >= 1
      !> only (Inf on a nucleus)
      real(wp), intent(in) :: inv_x_a
      !> Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade, gate-masked: stage-1
      !> `cache_u` (inline exp), 0 where the strip's hard screening dropped the atom
      real(wp), intent(in) :: u_a
      !> Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, as passed to the strip kernel; > 0
      real(wp), intent(in) :: fade_width
   end subroutine svdw_atom_radrad_eval

   !> Radius derivative of the joint-direction contracted tensors of one atom.
   !>
   !> Sets `awd{n}(i1..in, kind)` = `d/dR_a` of atom A's share `aw{n} + awr{n}` of
   !> the joint-direction `ws{n}`:
   !>
   !>    awd{n} = -sum_alpha v_a(alpha) atd{n+1}(.., alpha, kind) + vr_a atdd{n}
   !>
   !> what `svdw_radius_hvp_eval` reads for the radius derivative of `ws`.  Pass
   !> `vr_a = 0` for a pure nuclear direction.  Same geometry cache, fade and
   !> branch as `svdw_atom_eval`.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    d_a          Displacement r - R_A of the active atom
   !> @param[in]    R_a          Radius of the active atom
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    v_a          Nuclear displacement direction of this atom
   !> @param[in]    vr_a         Radius direction of this atom, vr_A = sum_B (dR_A/dX_B) . v_B
   !> @param[in]    max_deriv    Highest spatial-derivative order to set (0..2); a larger value is
   !>                            treated as 2
   !> @param[out]   awd0         Order-0 radius-differentiated joint contracted kind tensor of this
   !>                            atom
   !> @param[out]   awd1         Order-1 radius-differentiated joint contracted kind tensor of this
   !>                            atom
   !> @param[out]   awd2         Order-2 radius-differentiated joint contracted kind tensor of this
   !>                            atom
   !> @param[in]    x_a          Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and
   !>                            atom
   !> @param[in]    inv_x_a      Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded
   !>                            branch at orders >= 1 only (Inf on a nucleus)
   !> @param[in]    u_a          Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade,
   !>                            gate-masked: stage-1 `cache_u` (inline exp), 0 where the strip's
   !>                            hard screening dropped the atom
   !> @param[in]    fade_start   Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp)
   !>                            for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, as passed to the strip kernel; > 0
   pure module subroutine svdw_atom_tangent_rad_eval(d_a, R_a, blend_k, v_a, vr_a, max_deriv, awd0, &
                                              awd1, awd2, x_a, inv_x_a, u_a, fade_start, &
                                              fade_width)
      !> Displacement r - R_A of the active atom
      real(wp), intent(in) :: d_a(3)
      !> Radius of the active atom
      real(wp), intent(in) :: R_a
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Nuclear displacement direction of this atom
      real(wp), intent(in) :: v_a(3)
      !> Radius direction of this atom, vr_A = sum_B (dR_A/dX_B) . v_B
      real(wp), intent(in) :: vr_a
      !> Highest spatial-derivative order to set (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> Order-0 radius-differentiated joint contracted kind tensor of this atom
      real(wp), intent(out) :: awd0(nkind)
      !> Order-1 radius-differentiated joint contracted kind tensor of this atom
      real(wp), intent(out) :: awd1(3, nkind)
      !> Order-2 radius-differentiated joint contracted kind tensor of this atom
      real(wp), intent(out) :: awd2(3, 3, nkind)
      !> Distance ||d_a|| of this atom, stage-1 `cache_x` of this point and atom
      real(wp), intent(in) :: x_a
      !> Reciprocal distance 1/x_a, stage-1 `cache_inv_x`; read by the faded branch at orders >= 1
      !> only (Inf on a nucleus)
      real(wp), intent(in) :: inv_x_a
      !> Screening factor exp(-(blend_k/3)(x_a - R_a)) before the fade, gate-masked: stage-1
      !> `cache_u` (inline exp), 0 where the strip's hard screening dropped the atom
      real(wp), intent(in) :: u_a
      !> Fade onset X_on in d_A, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, as passed to the strip kernel; > 0
      real(wp), intent(in) :: fade_width
   end subroutine svdw_atom_tangent_rad_eval

   !> Level-set value and its pure spatial derivatives.
   !>
   !> Assembles `Z = s_1 p3 + (s_2/2)(q^2 - p3) + (s_3/6)(p1^3 - 3 p1 p2 + 2 p3)`
   !> from the accumulated power sums and returns `S = -log(Z)/blend_k` together
   !> with `d^n S / dr^n`, n = 1..max_deriv.  O(1): the O(n_active) work is the
   !> accumulation of `ps*` that the caller did beforehand.
   !>
   !> `f0` is the single shared expression `-log(Z)/blend_k`; there is no second
   !> spelling of it anywhere in this module, so value and gradient stay
   !> bit-for-bit consistent under `-ffp-contract=fast`.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps3          Accumulated power-sum tensor d^3 p_kind / dr^3, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps4          Accumulated power-sum tensor d^4 p_kind / dr^4, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    max_deriv    Highest spatial-derivative order to compute (0..4); a larger value
   !>                            is treated as 4
   !> @param[out]   f0           Level-set value S = -log(Z)/blend_k
   !> @param[out]   f1_r         Spatial gradient dS/dr
   !> @param[out]   f2_rr        Spatial Hessian d^2S/dr^2
   !> @param[out]   f3_rrr       Third spatial derivative d^3S/dr^3
   !> @param[out]   f4_rrrr      Fourth spatial derivative d^4S/dr^4
   pure module subroutine svdw_spatial_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, ps3, ps4, max_deriv, &
                                     f0, f1_r, f2_rr, f3_rrr, f4_rrrr)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Accumulated power-sum tensor d^3 p_kind / dr^3, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps3(3, 3, 3, nkind)
      !> Accumulated power-sum tensor d^4 p_kind / dr^4, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps4(3, 3, 3, 3, nkind)
      !> Highest spatial-derivative order to compute (0..4); a larger value is treated as 4
      integer, intent(in) :: max_deriv
      !> Level-set value S = -log(Z)/blend_k
      real(wp), intent(out) :: f0
      !> Spatial gradient dS/dr
      real(wp), intent(out) :: f1_r(3)
      !> Spatial Hessian d^2S/dr^2
      real(wp), intent(out) :: f2_rr(3, 3)
      !> Third spatial derivative d^3S/dr^3
      real(wp), intent(out) :: f3_rrr(3, 3, 3)
      !> Fourth spatial derivative d^4S/dr^4
      real(wp), intent(out) :: f4_rrrr(3, 3, 3, 3)
   end subroutine svdw_spatial_eval

   !> Nuclear derivative of the spatial derivatives, for one active atom A.
   !>
   !> Returns `d/dR_A,s  d^n S / dr^n` for n = 0..max_deriv, i.e. `f1_rA`,
   !> `f2_r_rA`, `f3_rr_rA` and `f4_rrr_rA`.  The trailing index of every result
   !> is the nuclear component s; the leading ones are spatial and fully
   !> symmetric among themselves.
   !>
   !> Because each power sum is a sum over single-atom terms, `d/dR_A` of an
   !> aggregate tensor collapses onto atom A's own tensor with one extra index
   !> and a minus sign.  The routine is therefore O(1) per atom and O(n_active)
   !> over the atom loop, with no pair or triple expansion anywhere.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps3          Accumulated power-sum tensor d^3 p_kind / dr^3, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    at0          Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at1          Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at2          Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at3          Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at4          Atom-A kind tensor d^4 u_A^m / dr^4 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    max_deriv    Highest *spatial* derivative order to compute (0..3); a larger
   !>                            value is treated as 3
   !> @param[out]   f1_rA        dS/dR_A
   !> @param[out]   f2_r_rA      d^2S/(dr dR_A); index order (spatial, nuclear)
   !> @param[out]   f3_rr_rA     d^3S/(dr^2 dR_A); index order (spatial, spatial, nuclear)
   !> @param[out]   f4_rrr_rA    d^4S/(dr^3 dR_A); index order (spatial x 3, nuclear)
   pure module subroutine svdw_nuclear_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, ps3, at0, at1, at2, &
                                     at3, at4, max_deriv, f1_rA, f2_r_rA, f3_rr_rA, f4_rrr_rA)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Accumulated power-sum tensor d^3 p_kind / dr^3, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps3(3, 3, 3, nkind)
      !> Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at0(nkind)
      !> Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at1(3, nkind)
      !> Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at2(3, 3, nkind)
      !> Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at3(3, 3, 3, nkind)
      !> Atom-A kind tensor d^4 u_A^m / dr^4 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at4(3, 3, 3, 3, nkind)
      !> Highest *spatial* derivative order to compute (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> dS/dR_A
      real(wp), intent(out) :: f1_rA(3)
      !> d^2S/(dr dR_A); index order (spatial, nuclear)
      real(wp), intent(out) :: f2_r_rA(3, 3)
      !> d^3S/(dr^2 dR_A); index order (spatial, spatial, nuclear)
      real(wp), intent(out) :: f3_rr_rA(3, 3, 3)
      !> d^4S/(dr^3 dR_A); index order (spatial x 3, nuclear)
      real(wp), intent(out) :: f4_rrr_rA(3, 3, 3, 3)
   end subroutine svdw_nuclear_eval

   !> Radius derivative of the spatial derivatives, for one active atom A.
   !>
   !> Returns `d/dR_a  d^n S / dr^n` for n = 0..max_deriv, i.e. `f1_rad`,
   !> `f2_r_rad`, `f3_rr_rad` and `f4_rrr_rad`, where `R_a` is the *radius* of
   !> atom A.  A radius is a scalar, so unlike `svdw_nuclear_eval` these results
   !> carry no trailing derivative index: the rank equals the spatial order.
   !>
   !> The radius enters only through `d_A = ||r - R_A|| - R_a`, so `d/dR_a` of a
   !> power sum collapses onto atom A's radius tensors `atd` (the spatial
   !> derivatives of `-dg/dd`, from `svdw_atom_rad_eval`).  No derivative order is
   !> consumed: this routine reaches the same spatial order as
   !> `svdw_nuclear_eval`, with `atd` in place of `at`.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps3          Accumulated power-sum tensor d^3 p_kind / dr^3, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    atd0         Atom-A radius tensor d/dR_a of the order-0 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd1         Atom-A radius tensor d/dR_a of the order-1 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd2         Atom-A radius tensor d/dR_a of the order-2 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd3         Atom-A radius tensor d/dR_a of the order-3 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    max_deriv    Highest *spatial* derivative order to compute (0..3); a larger
   !>                            value is treated as 3
   !> @param[out]   f1_rad       dS/dR_a
   !> @param[out]   f2_r_rad     d^2S/(dr dR_a)
   !> @param[out]   f3_rr_rad    d^3S/(dr^2 dR_a)
   !> @param[out]   f4_rrr_rad   d^4S/(dr^3 dR_a)
   pure module subroutine svdw_radius_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, ps3, atd0, atd1, atd2, &
                                    atd3, max_deriv, f1_rad, f2_r_rad, f3_rr_rad, f4_rrr_rad)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Accumulated power-sum tensor d^3 p_kind / dr^3, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps3(3, 3, 3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-0 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd0(nkind)
      !> Atom-A radius tensor d/dR_a of the order-1 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd1(3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-2 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd2(3, 3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-3 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd3(3, 3, 3, nkind)
      !> Highest *spatial* derivative order to compute (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> dS/dR_a
      real(wp), intent(out) :: f1_rad
      !> d^2S/(dr dR_a)
      real(wp), intent(out) :: f2_r_rad(3)
      !> d^3S/(dr^2 dR_a)
      real(wp), intent(out) :: f3_rr_rad(3, 3)
      !> d^4S/(dr^3 dR_a)
      real(wp), intent(out) :: f4_rrr_rad(3, 3, 3)
   end subroutine svdw_radius_eval

   !> Second nuclear derivative of the spatial derivatives, for A /= B.
   !>
   !> Returns `d^2/(dR_A,s dR_B,t) d^n S / dr^n` for n = 0..max_deriv, i.e.
   !> `f2_rArB`, `f3_r_rArB` and `f4_rr_rArB`.  Index order is
   !> (spatial..., s, t); the spatial indices are fully symmetric.
   !>
   !> A power sum has *no* cross-atom second nuclear derivative, so every term
   !> here is bilinear in the two atoms' own tensors: all A-B coupling comes from
   !> the nonlinearity of Z in the power sums and from the log.  Call
   !> `svdw_pair_diag_eval` instead when A == B, which adds the extra
   !> single-tensor term that `d/dR_A at = -at` produces.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    at0          Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at1          Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at2          Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at3          Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    bt0          Atom-B kind tensor d^0 u_B^m / dr^0 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    bt1          Atom-B kind tensor d^1 u_B^m / dr^1 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    bt2          Atom-B kind tensor d^2 u_B^m / dr^2 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    bt3          Atom-B kind tensor d^3 u_B^m / dr^3 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    max_deriv    Highest *spatial* derivative order to compute (0..2); a larger
   !>                            value is treated as 2
   !> @param[out]   f2_rArB      d^2S/(dR_A dR_B); index order (s, t)
   !> @param[out]   f3_r_rArB    d^3S/(dr dR_A dR_B); index order (spatial, s, t)
   !> @param[out]   f4_rr_rArB   d^4S/(dr^2 dR_A dR_B); index order (spatial, spatial, s, t)
   pure module subroutine svdw_pair_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, at0, at1, at2, at3, bt0, &
                                  bt1, bt2, bt3, max_deriv, f2_rArB, f3_r_rArB, f4_rr_rArB)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at0(nkind)
      !> Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at1(3, nkind)
      !> Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at2(3, 3, nkind)
      !> Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at3(3, 3, 3, nkind)
      !> Atom-B kind tensor d^0 u_B^m / dr^0 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: bt0(nkind)
      !> Atom-B kind tensor d^1 u_B^m / dr^1 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: bt1(3, nkind)
      !> Atom-B kind tensor d^2 u_B^m / dr^2 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: bt2(3, 3, nkind)
      !> Atom-B kind tensor d^3 u_B^m / dr^3 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: bt3(3, 3, 3, nkind)
      !> Highest *spatial* derivative order to compute (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d^2S/(dR_A dR_B); index order (s, t)
      real(wp), intent(out) :: f2_rArB(3, 3)
      !> d^3S/(dr dR_A dR_B); index order (spatial, s, t)
      real(wp), intent(out) :: f3_r_rArB(3, 3, 3)
      !> d^4S/(dr^2 dR_A dR_B); index order (spatial, spatial, s, t)
      real(wp), intent(out) :: f4_rr_rArB(3, 3, 3, 3)
   end subroutine svdw_pair_eval

   !> Second nuclear derivative of the spatial derivatives, for A == B.
   !>
   !> Same outputs and index order as `svdw_pair_eval`, but with both nuclear
   !> derivatives taken with respect to the *same* atom, which contributes the
   !> extra term `d/dR_A at_k[J] = -at_k[J + s]` that the off-diagonal routine
   !> does not have.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    at0          Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at1          Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at2          Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at3          Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at4          Atom-A kind tensor d^4 u_A^m / dr^4 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    max_deriv    Highest *spatial* derivative order to compute (0..2); a larger
   !>                            value is treated as 2
   !> @param[out]   f2_rArB      d^2S/(dR_A dR_B); index order (s, t)
   !> @param[out]   f3_r_rArB    d^3S/(dr dR_A dR_B); index order (spatial, s, t)
   !> @param[out]   f4_rr_rArB   d^4S/(dr^2 dR_A dR_B); index order (spatial, spatial, s, t)
   pure module subroutine svdw_pair_diag_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, at0, at1, at2, at3, &
                                       at4, max_deriv, f2_rArB, f3_r_rArB, f4_rr_rArB)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at0(nkind)
      !> Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at1(3, nkind)
      !> Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at2(3, 3, nkind)
      !> Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at3(3, 3, 3, nkind)
      !> Atom-A kind tensor d^4 u_A^m / dr^4 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at4(3, 3, 3, 3, nkind)
      !> Highest *spatial* derivative order to compute (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d^2S/(dR_A dR_B); index order (s, t)
      real(wp), intent(out) :: f2_rArB(3, 3)
      !> d^3S/(dr dR_A dR_B); index order (spatial, s, t)
      real(wp), intent(out) :: f3_r_rArB(3, 3, 3)
      !> d^4S/(dr^2 dR_A dR_B); index order (spatial, spatial, s, t)
      real(wp), intent(out) :: f4_rr_rArB(3, 3, 3, 3)
   end subroutine svdw_pair_diag_eval

   !> Second *radius* derivative of the spatial derivatives, for A /= B.
   !>
   !> Returns `d^2/(dR_a dR_b) d^n S / dr^n` for n = 0..max_deriv, i.e.
   !> `f2_radrad`, `f3_r_radrad` and `f4_rr_radrad`, where `R_a` and `R_b` are the
   !> *radii* of atoms A and B.  Radii are scalars, so neither derivative appends
   !> an index: the rank of every result is the spatial order, two less than the
   !> corresponding `svdw_pair_eval` output.
   !>
   !> Reads the radius tensors `atd` and `btd` of the two atoms
   !> (`svdw_atom_rad_eval`).  The atom-A families carry no `R_b` dependence, so the
   !> entire cross-radius coupling comes from the nonlinearity of Z in the power
   !> sums and from the log -- the same statement as for `f2_rArB`.  Call
   !> `svdw_radpair_diag_eval` instead when A == B, which adds the extra
   !> single-tensor term.
   !>
   !> This is the uncontracted form of `svdw_radius_hvp_eval`'s radius-radius
   !> half; prefer that routine unless the O(n_active**2) tensor itself is wanted.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    atd0         Atom-A radius tensor d/dR_a of the order-0 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd1         Atom-A radius tensor d/dR_a of the order-1 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd2         Atom-A radius tensor d/dR_a of the order-2 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    btd0         Atom-B radius tensor d/dR_b of the order-0 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    btd1         Atom-B radius tensor d/dR_b of the order-1 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    btd2         Atom-B radius tensor d/dR_b of the order-2 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    max_deriv    Highest *spatial* derivative order to compute (0..2); a larger
   !>                            value is treated as 2
   !> @param[out]   f2_radrad    d^2S/(dR_a dR_b)
   !> @param[out]   f3_r_radrad  d^3S/(dr dR_a dR_b)
   !> @param[out]   f4_rr_radrad d^4S/(dr^2 dR_a dR_b)
   pure module subroutine svdw_radpair_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, atd0, atd1, atd2, &
                                     btd0, btd1, btd2, max_deriv, f2_radrad, f3_r_radrad, &
                                     f4_rr_radrad)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-0 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd0(nkind)
      !> Atom-A radius tensor d/dR_a of the order-1 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd1(3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-2 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd2(3, 3, nkind)
      !> Atom-B radius tensor d/dR_b of the order-0 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: btd0(nkind)
      !> Atom-B radius tensor d/dR_b of the order-1 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: btd1(3, nkind)
      !> Atom-B radius tensor d/dR_b of the order-2 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: btd2(3, 3, nkind)
      !> Highest *spatial* derivative order to compute (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d^2S/(dR_a dR_b)
      real(wp), intent(out) :: f2_radrad
      !> d^3S/(dr dR_a dR_b)
      real(wp), intent(out) :: f3_r_radrad(3)
      !> d^4S/(dr^2 dR_a dR_b)
      real(wp), intent(out) :: f4_rr_radrad(3, 3)
   end subroutine svdw_radpair_eval

   !> Second *radius* derivative of the spatial derivatives, for A == B.
   !>
   !> Same outputs and ranks as `svdw_radpair_eval`, but with both derivatives
   !> taken with respect to the same radius, which contributes the extra term
   !> `d/dR_a atd_k[J] = atdd_k[J]` (`svdw_atom_radrad_eval`) that the
   !> off-diagonal routine does not have.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    atd0         Atom-A radius tensor d/dR_a of the order-0 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd1         Atom-A radius tensor d/dR_a of the order-1 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd2         Atom-A radius tensor d/dR_a of the order-2 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atdd0        Atom-A second radius tensor d^2/dR_a^2 of the order-0 kind tensor
   !>                            from svdw_atom_radrad_eval; last index selects the kind
   !> @param[in]    atdd1        Atom-A second radius tensor d^2/dR_a^2 of the order-1 kind tensor
   !>                            from svdw_atom_radrad_eval; last index selects the kind
   !> @param[in]    atdd2        Atom-A second radius tensor d^2/dR_a^2 of the order-2 kind tensor
   !>                            from svdw_atom_radrad_eval; last index selects the kind
   !> @param[in]    max_deriv    Highest *spatial* derivative order to compute (0..2); a larger
   !>                            value is treated as 2
   !> @param[out]   f2_radrad    d^2S/(dR_a dR_b)
   !> @param[out]   f3_r_radrad  d^3S/(dr dR_a dR_b)
   !> @param[out]   f4_rr_radrad d^4S/(dr^2 dR_a dR_b)
   pure module subroutine svdw_radpair_diag_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, atd0, atd1, &
                                          atd2, atdd0, atdd1, atdd2, max_deriv, f2_radrad, &
                                          f3_r_radrad, f4_rr_radrad)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-0 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd0(nkind)
      !> Atom-A radius tensor d/dR_a of the order-1 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd1(3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-2 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd2(3, 3, nkind)
      !> Atom-A second radius tensor d^2/dR_a^2 of the order-0 kind tensor from
      !> svdw_atom_radrad_eval; last index selects the kind
      real(wp), intent(in) :: atdd0(nkind)
      !> Atom-A second radius tensor d^2/dR_a^2 of the order-1 kind tensor from
      !> svdw_atom_radrad_eval; last index selects the kind
      real(wp), intent(in) :: atdd1(3, nkind)
      !> Atom-A second radius tensor d^2/dR_a^2 of the order-2 kind tensor from
      !> svdw_atom_radrad_eval; last index selects the kind
      real(wp), intent(in) :: atdd2(3, 3, nkind)
      !> Highest *spatial* derivative order to compute (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d^2S/(dR_a dR_b)
      real(wp), intent(out) :: f2_radrad
      !> d^3S/(dr dR_a dR_b)
      real(wp), intent(out) :: f3_r_radrad(3)
      !> d^4S/(dr^2 dR_a dR_b)
      real(wp), intent(out) :: f4_rr_radrad(3, 3)
   end subroutine svdw_radpair_diag_eval

   !> Mixed nuclear-radius derivative of the spatial derivatives, for A /= B.
   !>
   !> Returns `d^2/(dR_A,s dR_b) d^n S / dr^n` for n = 0..max_deriv, i.e.
   !> `f2_rA_rad`, `f3_r_rA_rad` and `f4_rr_rA_rad`: a *position* derivative on
   !> atom A and a *radius* derivative on atom B.  Index order is (spatial..., s);
   !> the radius contributes no index, so the rank is one above the spatial order.
   !>
   !> The two slots are not interchangeable.  `f{n}_rA_rad(.., s)` for the atom
   !> pair (A, B) is *not* the transpose of the same array for (B, A) -- the first
   !> asks for A's position and B's radius, the second for B's position and A's
   !> radius.  A caller that needs the full nuclear-radius block must call this routine for
   !> both orderings, and `svdw_nucrad_diag_eval` on the diagonal.
   !>
   !> Reads atom A's `at` and atom B's radius tensors `btd` (`svdw_atom_rad_eval`).
   !> This is the uncontracted form of the nuclear-radius block that
   !> `svdw_hvp_eval` and `svdw_radius_hvp_eval` already supply contracted.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    at0          Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at1          Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at2          Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at3          Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    btd0         Atom-B radius tensor d/dR_b of the order-0 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    btd1         Atom-B radius tensor d/dR_b of the order-1 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    btd2         Atom-B radius tensor d/dR_b of the order-2 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    max_deriv    Highest *spatial* derivative order to compute (0..2); a larger
   !>                            value is treated as 2
   !> @param[out]   f2_rA_rad    d^2S/(dR_A dR_b); trailing index is A's component
   !> @param[out]   f3_r_rA_rad  d^3S/(dr dR_A dR_b); index order (spatial, s)
   !> @param[out]   f4_rr_rA_rad d^4S/(dr^2 dR_A dR_b); index order (spatial, spatial, s)
   pure module subroutine svdw_nucrad_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, at0, at1, at2, at3, &
                                    btd0, btd1, btd2, max_deriv, f2_rA_rad, f3_r_rA_rad, &
                                    f4_rr_rA_rad)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at0(nkind)
      !> Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at1(3, nkind)
      !> Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at2(3, 3, nkind)
      !> Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at3(3, 3, 3, nkind)
      !> Atom-B radius tensor d/dR_b of the order-0 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: btd0(nkind)
      !> Atom-B radius tensor d/dR_b of the order-1 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: btd1(3, nkind)
      !> Atom-B radius tensor d/dR_b of the order-2 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: btd2(3, 3, nkind)
      !> Highest *spatial* derivative order to compute (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d^2S/(dR_A dR_b); trailing index is A's component
      real(wp), intent(out) :: f2_rA_rad(3)
      !> d^3S/(dr dR_A dR_b); index order (spatial, s)
      real(wp), intent(out) :: f3_r_rA_rad(3, 3)
      !> d^4S/(dr^2 dR_A dR_b); index order (spatial, spatial, s)
      real(wp), intent(out) :: f4_rr_rA_rad(3, 3, 3)
   end subroutine svdw_nucrad_eval

   !> Mixed nuclear-radius derivative of the spatial derivatives, for A == B.
   !>
   !> Same outputs and index order as `svdw_nucrad_eval`, but with the position and
   !> the radius belonging to the *same* atom, so `d/dR_a` still sees atom A's own
   !> tensors instead of annihilating them: `at -> atd` (`svdw_atom_rad_eval`).
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    at0          Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at1          Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at2          Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at3          Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    atd0         Atom-A radius tensor d/dR_a of the order-0 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd1         Atom-A radius tensor d/dR_a of the order-1 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd2         Atom-A radius tensor d/dR_a of the order-2 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd3         Atom-A radius tensor d/dR_a of the order-3 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    max_deriv    Highest *spatial* derivative order to compute (0..2); a larger
   !>                            value is treated as 2
   !> @param[out]   f2_rA_rad    d^2S/(dR_A dR_b); trailing index is A's component
   !> @param[out]   f3_r_rA_rad  d^3S/(dr dR_A dR_b); index order (spatial, s)
   !> @param[out]   f4_rr_rA_rad d^4S/(dr^2 dR_A dR_b); index order (spatial, spatial, s)
   pure module subroutine svdw_nucrad_diag_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, at0, at1, at2, &
                                         at3, atd0, atd1, atd2, atd3, max_deriv, f2_rA_rad, &
                                         f3_r_rA_rad, f4_rr_rA_rad)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at0(nkind)
      !> Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at1(3, nkind)
      !> Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at2(3, 3, nkind)
      !> Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at3(3, 3, 3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-0 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd0(nkind)
      !> Atom-A radius tensor d/dR_a of the order-1 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd1(3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-2 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd2(3, 3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-3 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd3(3, 3, 3, nkind)
      !> Highest *spatial* derivative order to compute (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d^2S/(dR_A dR_b); trailing index is A's component
      real(wp), intent(out) :: f2_rA_rad(3)
      !> d^3S/(dr dR_A dR_b); index order (spatial, s)
      real(wp), intent(out) :: f3_r_rA_rad(3, 3)
      !> d^4S/(dr^2 dR_A dR_b); index order (spatial, spatial, s)
      real(wp), intent(out) :: f4_rr_rA_rad(3, 3, 3)
   end subroutine svdw_nucrad_diag_eval

   !> Directional nuclear derivative of the spatial derivatives (contracted).
   !>
   !> Given one displacement direction v_B per atom, returns
   !>
   !>    tangent_f{n}(i1..in) = sum_B sum_alpha v_B,alpha d/dR_B,alpha
   !>                           d^n S / dr_i1..dr_in
   !>
   !> for n = 0..max_deriv: the single nuclear index is contracted away and the
   !> result is a purely spatial tensor.  These are derivatives of `S` itself,
   !> i.e. the `-log(Z)/blend_k` lift is already applied; nothing is left for the
   !> caller to lift.
   !>
   !> The contraction is done at the `u_A` level by `svdw_atom_tangent_eval`, so
   !> the whole chain is O(n_active) and no pairwise tensor is ever built.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps3          Accumulated power-sum tensor d^3 p_kind / dr^3, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ws0          Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the
   !>                            order-0 power-sum tensor; last index selects the kind
   !> @param[in]    ws1          Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the
   !>                            order-1 power-sum tensor; last index selects the kind
   !> @param[in]    ws2          Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the
   !>                            order-2 power-sum tensor; last index selects the kind
   !> @param[in]    ws3          Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the
   !>                            order-3 power-sum tensor; last index selects the kind
   !> @param[in]    max_deriv    Highest spatial-derivative order to compute (0..3); a larger value
   !>                            is treated as 3
   !> @param[out]   tangent_f0   sum_B v_B . dS/dR_B
   !> @param[out]   tangent_f1_r sum_B v_B . d/dR_B of dS/dr
   !> @param[out]   tangent_f2_rr sum_B v_B . d/dR_B of d^2S/dr^2
   !> @param[out]   tangent_f3_rrr sum_B v_B . d/dR_B of d^3S/dr^3
   pure module subroutine svdw_tangent_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, ps3, ws0, ws1, ws2, &
                                     ws3, max_deriv, tangent_f0, tangent_f1_r, tangent_f2_rr, &
                                     tangent_f3_rrr)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Accumulated power-sum tensor d^3 p_kind / dr^3, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps3(3, 3, 3, nkind)
      !> Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the order-0 power-sum tensor;
      !> last index selects the kind
      real(wp), intent(in) :: ws0(nkind)
      !> Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the order-1 power-sum tensor;
      !> last index selects the kind
      real(wp), intent(in) :: ws1(3, nkind)
      !> Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the order-2 power-sum tensor;
      !> last index selects the kind
      real(wp), intent(in) :: ws2(3, 3, nkind)
      !> Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the order-3 power-sum tensor;
      !> last index selects the kind
      real(wp), intent(in) :: ws3(3, 3, 3, nkind)
      !> Highest spatial-derivative order to compute (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_B v_B . dS/dR_B
      real(wp), intent(out) :: tangent_f0
      !> sum_B v_B . d/dR_B of dS/dr
      real(wp), intent(out) :: tangent_f1_r(3)
      !> sum_B v_B . d/dR_B of d^2S/dr^2
      real(wp), intent(out) :: tangent_f2_rr(3, 3)
      !> sum_B v_B . d/dR_B of d^3S/dr^3
      real(wp), intent(out) :: tangent_f3_rrr(3, 3, 3)
   end subroutine svdw_tangent_eval

   !> Nuclear Hessian-vector product of the spatial derivatives, for one atom A.
   !>
   !> Given one displacement direction v_B per atom, returns
   !>
   !>    hvp_f{n+1}(i1..in, s) = sum_B sum_alpha v_B,alpha
   !>                            d^2/(dR_A,s dR_B,alpha) d^n S / dr_i1..dr_in
   !>
   !> for n = 0..max_deriv: exactly one of the two nuclear indices is contracted,
   !> the other is retained as the trailing index.  These are derivatives of `S`
   !> itself, with the log lift already applied.
   !>
   !> One call per active atom turns what would be an O(n_active^2) tensor into
   !> O(n_active) storage and O(n_active) work.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ws0          Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the
   !>                            order-0 power-sum tensor; last index selects the kind
   !> @param[in]    ws1          Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the
   !>                            order-1 power-sum tensor; last index selects the kind
   !> @param[in]    ws2          Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the
   !>                            order-2 power-sum tensor; last index selects the kind
   !> @param[in]    at0          Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at1          Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at2          Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at3          Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    aw0          Atom-A direction-contracted kind tensor from
   !>                            svdw_atom_tangent_eval; last index selects the kind
   !> @param[in]    aw1          Atom-A direction-contracted kind tensor from
   !>                            svdw_atom_tangent_eval; last index selects the kind
   !> @param[in]    aw2          Atom-A direction-contracted kind tensor from
   !>                            svdw_atom_tangent_eval; last index selects the kind
   !> @param[in]    aw3          Atom-A direction-contracted kind tensor from
   !>                            svdw_atom_tangent_eval; last index selects the kind
   !> @param[in]    max_deriv    Highest spatial-derivative order to compute (0..2); a larger value
   !>                            is treated as 2
   !> @param[out]   hvp_f1_rA    sum_B v_B . d^2S/(dR_A dR_B); trailing index is A's component
   !> @param[out]   hvp_f2_r_rA  sum_B v_B . d^3S/(dr dR_A dR_B); trailing index is A's component
   !> @param[out]   hvp_f3_rr_rA sum_B v_B . d^4S/(dr^2 dR_A dR_B); trailing index is A's component
   pure module subroutine svdw_hvp_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, ws0, ws1, ws2, at0, at1, &
                                 at2, at3, aw0, aw1, aw2, aw3, max_deriv, hvp_f1_rA, hvp_f2_r_rA, &
                                 hvp_f3_rr_rA)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the order-0 power-sum tensor;
      !> last index selects the kind
      real(wp), intent(in) :: ws0(nkind)
      !> Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the order-1 power-sum tensor;
      !> last index selects the kind
      real(wp), intent(in) :: ws1(3, nkind)
      !> Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the order-2 power-sum tensor;
      !> last index selects the kind
      real(wp), intent(in) :: ws2(3, 3, nkind)
      !> Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at0(nkind)
      !> Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at1(3, nkind)
      !> Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at2(3, 3, nkind)
      !> Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at3(3, 3, 3, nkind)
      !> Atom-A direction-contracted kind tensor from svdw_atom_tangent_eval; last index selects the
      !> kind
      real(wp), intent(in) :: aw0(nkind)
      !> Atom-A direction-contracted kind tensor from svdw_atom_tangent_eval; last index selects the
      !> kind
      real(wp), intent(in) :: aw1(3, nkind)
      !> Atom-A direction-contracted kind tensor from svdw_atom_tangent_eval; last index selects the
      !> kind
      real(wp), intent(in) :: aw2(3, 3, nkind)
      !> Atom-A direction-contracted kind tensor from svdw_atom_tangent_eval; last index selects the
      !> kind
      real(wp), intent(in) :: aw3(3, 3, 3, nkind)
      !> Highest spatial-derivative order to compute (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> sum_B v_B . d^2S/(dR_A dR_B); trailing index is A's component
      real(wp), intent(out) :: hvp_f1_rA(3)
      !> sum_B v_B . d^3S/(dr dR_A dR_B); trailing index is A's component
      real(wp), intent(out) :: hvp_f2_r_rA(3, 3)
      !> sum_B v_B . d^4S/(dr^2 dR_A dR_B); trailing index is A's component
      real(wp), intent(out) :: hvp_f3_rr_rA(3, 3, 3)
   end subroutine svdw_hvp_eval

   !> Joint-direction Hessian-vector product with a retained *radius* derivative.
   !>
   !> Returns `sum_B (v_B . d/dR_B + vr_B d/dR_b) d/dR_a d^n S / dr^n` for
   !> n = 0..max_deriv, i.e. `hvp_f1_rad`, `hvp_f2_r_rad` and `hvp_f3_rr_rad`,
   !> where `R_a` is the radius of atom A.  The results carry no trailing index:
   !> the rank is the spatial order.
   !>
   !> This is the radius row of the cavity Hessian.  With `R = R(X)` the chain
   !> rule needs four blocks -- nuclear-nuclear, nuclear-radius, radius-nuclear
   !> and radius-radius -- and this routine supplies the two that retain a radius
   !> index, because the contraction direction is the joint pair `(v_B, vr_B)`.
   !> `svdw_hvp_eval` supplies the other two from the same `ws` inputs.
   !>
   !> `d/dR_a` takes `ps` to atom A's `atd` (`svdw_atom_rad_eval`) and the
   !> joint-direction `ws` to atom A's `awd` (`svdw_atom_tangent_rad_eval`).
   !>
   !> The radius-radius coupling between *different* atoms is nonzero and is not
   !> recoverable from per-atom data: the atom-A families carry no `R_b`
   !> dependence, so all of it comes from the nonlinearity of Z in the power sums
   !> and from the log, exactly as for `f2_rArB`.
   !>
   !> Results above `max_deriv` are left untouched.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ws0          Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the
   !>                            order-0 power-sum tensor; last index selects the kind
   !> @param[in]    ws1          Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the
   !>                            order-1 power-sum tensor; last index selects the kind
   !> @param[in]    ws2          Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the
   !>                            order-2 power-sum tensor; last index selects the kind
   !> @param[in]    atd0         Atom-A radius tensor d/dR_a of the order-0 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd1         Atom-A radius tensor d/dR_a of the order-1 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd2         Atom-A radius tensor d/dR_a of the order-2 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    awd0         Atom-A radius derivative d/dR_a of the order-0 joint contracted
   !>                            kind tensor (aw + awr) from svdw_atom_tangent_rad_eval; last index
   !>                            selects the kind
   !> @param[in]    awd1         Atom-A radius derivative d/dR_a of the order-1 joint contracted
   !>                            kind tensor (aw + awr) from svdw_atom_tangent_rad_eval; last index
   !>                            selects the kind
   !> @param[in]    awd2         Atom-A radius derivative d/dR_a of the order-2 joint contracted
   !>                            kind tensor (aw + awr) from svdw_atom_tangent_rad_eval; last index
   !>                            selects the kind
   !> @param[in]    max_deriv    Highest spatial-derivative order to compute (0..2); a larger value
   !>                            is treated as 2
   !> @param[out]   hvp_f1_rad   sum_B v_B . d^2S/(dR_a dR_B) + sum_B vr_B d^2S/(dR_a dR_b)
   !> @param[out]   hvp_f2_r_rad the same for dS/dr
   !> @param[out]   hvp_f3_rr_rad the same for d^2S/dr^2
   pure module subroutine svdw_radius_hvp_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, ws0, ws1, ws2, &
                                        atd0, atd1, atd2, awd0, awd1, awd2, max_deriv, &
                                        hvp_f1_rad, hvp_f2_r_rad, hvp_f3_rr_rad)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the order-0 power-sum tensor;
      !> last index selects the kind
      real(wp), intent(in) :: ws0(nkind)
      !> Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the order-1 power-sum tensor;
      !> last index selects the kind
      real(wp), intent(in) :: ws1(3, nkind)
      !> Accumulated direction-contracted tensor sum_B v_B . d/dR_B of the order-2 power-sum tensor;
      !> last index selects the kind
      real(wp), intent(in) :: ws2(3, 3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-0 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd0(nkind)
      !> Atom-A radius tensor d/dR_a of the order-1 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd1(3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-2 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd2(3, 3, nkind)
      !> Atom-A radius derivative d/dR_a of the order-0 joint contracted kind tensor (aw + awr) from
      !> svdw_atom_tangent_rad_eval; last index selects the kind
      real(wp), intent(in) :: awd0(nkind)
      !> Atom-A radius derivative d/dR_a of the order-1 joint contracted kind tensor (aw + awr) from
      !> svdw_atom_tangent_rad_eval; last index selects the kind
      real(wp), intent(in) :: awd1(3, nkind)
      !> Atom-A radius derivative d/dR_a of the order-2 joint contracted kind tensor (aw + awr) from
      !> svdw_atom_tangent_rad_eval; last index selects the kind
      real(wp), intent(in) :: awd2(3, 3, nkind)
      !> Highest spatial-derivative order to compute (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> sum_B v_B . d^2S/(dR_a dR_B) + sum_B vr_B d^2S/(dR_a dR_b)
      real(wp), intent(out) :: hvp_f1_rad
      !> the same for dS/dr
      real(wp), intent(out) :: hvp_f2_r_rad(3)
      !> the same for d^2S/dr^2
      real(wp), intent(out) :: hvp_f3_rr_rad(3, 3)
   end subroutine svdw_radius_hvp_eval

   !> Jet-contracted nuclear gradient of one active atom (reverse mode).
   !>
   !> Given the adjoint weights (`w0`, `w1`, `w2`) of one evaluation point,
   !> returns the single vector
   !>
   !>    vjp_f1_rA(s) = w0 * f1_rA(s) + sum_a w1(a) * f2_r_rA(a, s)
   !>                   + sum_a sum_b w2(a, b) * f3_rr_rA(a, b, s)
   !>
   !> i.e. the `f1_rA` rung of `svdw_nuclear_eval` with the spatial jet indices
   !> contracted away and the nuclear index retained.  `w2` is a *general* 3x3:
   !> all nine entries are contracted, no symmetry is assumed and no factor of
   !> two is folded in.
   !>
   !> This is the reverse-mode mirror of `svdw_tangent_eval`, which contracts the
   !> nuclear index and keeps the spatial ones.  The contraction is applied to the
   !> *spatial jet* before the nuclear derivative is taken: `d/dR_A` of the single
   !> scalar `w0 S + w1 . dS/dr + w2 : d^2S/dr^2`.  Since `d/dR_A` is a
   !> derivation, that is identical to differentiating all 1 + 3 + 9 jet
   !> components and contracting afterwards, but the emitted code is the partials
   !> of one collapsed scalar rather than thirty tensor components -- and those
   !> partials do not depend on the nuclear component, so all three components of
   !> the result share them.
   !>
   !> One call per active atom therefore replaces the (3 + 9 + 27) doubles per
   !> atom that a caller would otherwise fetch from `svdw_nuclear_eval` and
   !> contract itself.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    at0          Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at1          Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at2          Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at3          Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    w0           Adjoint weight of the level-set value S
   !> @param[in]    w1           Adjoint weights of the spatial gradient dS/dr
   !> @param[in]    w2           Adjoint weights of the spatial Hessian d^2S/dr^2; a general 3x3,
   !>                            contracted over all nine entries
   !> @param[out]   vjp_f1_rA    w0 f1_rA + w1 . f2_r_rA + w2 : f3_rr_rA of this atom; the index is
   !>                            the nuclear component s
   pure module subroutine svdw_vjp_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, at0, at1, at2, at3, w0, &
                                 w1, w2, vjp_f1_rA)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at0(nkind)
      !> Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at1(3, nkind)
      !> Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at2(3, 3, nkind)
      !> Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at3(3, 3, 3, nkind)
      !> Adjoint weight of the level-set value S
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient dS/dr
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the spatial Hessian d^2S/dr^2; a general 3x3, contracted over all nine
      !> entries
      real(wp), intent(in) :: w2(3, 3)
      !> w0 f1_rA + w1 . f2_r_rA + w2 : f3_rr_rA of this atom; the index is the nuclear component s
      real(wp), intent(out) :: vjp_f1_rA(3)
   end subroutine svdw_vjp_eval

   !> Jet-contracted radius derivative of one active atom (reverse mode).
   !>
   !> Given the adjoint weights (`w0`, `w1`, `w2`) of one evaluation point,
   !> returns the single scalar
   !>
   !>    vjp_f1_rad = w0 * f1_rad + sum_a w1(a) * f2_r_rad(a)
   !>                 + sum_a sum_b w2(a, b) * f3_rr_rad(a, b)
   !>
   !> i.e. the `f1_rad` rung of `svdw_radius_eval` with the spatial jet indices
   !> contracted away.  `w2` is a *general* 3x3: all nine entries are contracted,
   !> no symmetry is assumed and no factor of two is folded in.
   !>
   !> The radius twin of `svdw_vjp_eval`, and it inherits that routine's whole
   !> argument: the contraction is applied to the *spatial jet* before the radius
   !> derivative is taken, which is legitimate because `d/dR_a` is a derivation.
   !> What differs is the retained index -- a radius is a scalar, so there is
   !> none, and one call per active atom returns one number where the nuclear
   !> version returns three.  Against the uncontracted `svdw_radius_eval` that is
   !> 1 double per atom instead of 1 + 3 + 9.
   !>
   !> Because the radius consumes no derivative order, this routine reads the
   !> power sums and the radius tensors `atd` (`svdw_atom_rad_eval`) only up to
   !> order 2 -- one order less than `svdw_vjp_eval`, exactly as
   !> `svdw_radius_eval` needs one less than `svdw_nuclear_eval`.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    atd0         Atom-A radius tensor d/dR_a of the order-0 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd1         Atom-A radius tensor d/dR_a of the order-1 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    atd2         Atom-A radius tensor d/dR_a of the order-2 kind tensor from
   !>                            svdw_atom_rad_eval; last index selects the kind
   !> @param[in]    w0           Adjoint weight of the level-set value S
   !> @param[in]    w1           Adjoint weights of the spatial gradient dS/dr
   !> @param[in]    w2           Adjoint weights of the spatial Hessian d^2S/dr^2; a general 3x3,
   !>                            contracted over all nine entries
   !> @param[out]   vjp_f1_rad   w0 f1_rad + w1 . f2_r_rad + w2 : f3_rr_rad of this atom; a scalar,
   !>                            since a radius carries no derivative index
   pure module subroutine svdw_radius_vjp_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, atd0, atd1, atd2, &
                                        w0, w1, w2, vjp_f1_rad)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-0 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd0(nkind)
      !> Atom-A radius tensor d/dR_a of the order-1 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd1(3, nkind)
      !> Atom-A radius tensor d/dR_a of the order-2 kind tensor from svdw_atom_rad_eval; last index
      !> selects the kind
      real(wp), intent(in) :: atd2(3, 3, nkind)
      !> Adjoint weight of the level-set value S
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient dS/dr
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the spatial Hessian d^2S/dr^2; a general 3x3, contracted over all nine
      !> entries
      real(wp), intent(in) :: w2(3, 3)
      !> w0 f1_rad + w1 . f2_r_rad + w2 : f3_rr_rad of this atom; a scalar, since a radius carries
      !> no derivative index
      real(wp), intent(out) :: vjp_f1_rad
   end subroutine svdw_radius_vjp_eval

   !> The weighted hvp row of one active atom as a linear form in `ws`.
   !>
   !> `svdw_hvp_eval` contracted against the adjoint weights (`w0`, `w1`, `w2`)
   !> of one evaluation point is the row
   !>
   !>    row(s) = w0 hvp_f1_rA(s) + sum_a w1(a) hvp_f2_r_rA(a, s)
   !>             + sum_a sum_b w2(a, b) hvp_f3_rr_rA(a, b, s) ,
   !>
   !> and it is *linear* in its direction inputs, the aggregate `ws{n}` and the
   !> atom's own `aw{n}`, because a directional derivative is linear in the
   !> direction.  This routine returns the coefficients of the `ws` half,
   !>
   !>    row(s) = sum_kind [ hvp_vjp_ws0(s, kind) ws0(kind)
   !>                      + sum_a hvp_vjp_ws1(a, s, kind) ws1(a, kind)
   !>                      + sum_a sum_b hvp_vjp_ws2(a, b, s, kind) ws2(a, b, kind) ]
   !>             + row(s) evaluated by svdw_hvp_eval with ws = 0 ,
   !>
   !> valid for any symmetric `ws2`: the two orderings of a mixed index share one
   !> coefficient, split evenly between them (an exact operation), so contracting
   !> all nine entries counts the pair once.
   !>
   !> Why it exists.  A caller that needs the row for *many* directions at one
   !> point -- the direction-free nuclear Hessian block, one column per Cartesian
   !> unit direction of every active atom -- would otherwise call `svdw_hvp_eval`
   !> once per direction per atom, O(n_active**2) kernel evaluations per point.
   !> With `ws(v) = sum_B aw_B(v)` and `aw_B(e_t) = -bt{n+1}(.., t)` under the
   !> contraction rule, the whole `ws` half of the block is one call per atom
   !> followed by a rank-52 matrix product of per-atom quantities; the `aw`
   !> half, live only where the moving atom is the atom itself, is the own-atom
   !> block of `svdw_vjp_diag_eval` less this form applied to `aw_A` (see there).
   !>
   !> Reads the power sums and the atom's tensors to the same orders as
   !> `svdw_hvp_eval` at `max_deriv = 2`, and takes no `max_deriv` of its own.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    at0          Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at1          Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at2          Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at3          Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    w0           Adjoint weight of the level-set value S
   !> @param[in]    w1           Adjoint weights of the spatial gradient dS/dr
   !> @param[in]    w2           Adjoint weights of the spatial Hessian d^2S/dr^2; a general 3x3,
   !>                            contracted over all nine entries
   !> @param[out]   hvp_vjp_ws0  d/d ws0 of the weighted hvp row; index order (nuclear s, kind)
   !> @param[out]   hvp_vjp_ws1  d/d ws1 of the weighted hvp row; index order (spatial, nuclear s,
   !>                            kind)
   !> @param[out]   hvp_vjp_ws2  d/d ws2 of the weighted hvp row; index order (spatial, spatial,
   !>                            nuclear s, kind), the symmetric pair's coefficient split evenly
   !>                            over both orderings
   pure module subroutine svdw_hvp_vjp_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, at0, at1, at2, at3, &
                                     w0, w1, w2, hvp_vjp_ws0, hvp_vjp_ws1, hvp_vjp_ws2)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at0(nkind)
      !> Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at1(3, nkind)
      !> Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at2(3, 3, nkind)
      !> Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at3(3, 3, 3, nkind)
      !> Adjoint weight of the level-set value S
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient dS/dr
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the spatial Hessian d^2S/dr^2; a general 3x3, contracted over all nine
      !> entries
      real(wp), intent(in) :: w2(3, 3)
      !> d/d ws0 of the weighted hvp row; index order (nuclear s, kind)
      real(wp), intent(out) :: hvp_vjp_ws0(3, nkind)
      !> d/d ws1 of the weighted hvp row; index order (spatial, nuclear s, kind)
      real(wp), intent(out) :: hvp_vjp_ws1(3, 3, nkind)
      !> d/d ws2 of the weighted hvp row; index order (spatial, spatial, nuclear s, kind), the
      !> symmetric pair's coefficient split evenly over both orderings
      real(wp), intent(out) :: hvp_vjp_ws2(3, 3, 3, nkind)
   end subroutine svdw_hvp_vjp_eval

   !> Jet-contracted own-atom nuclear Hessian block of one active atom.
   !>
   !> Given the adjoint weights (`w0`, `w1`, `w2`) of one evaluation point,
   !> returns the symmetric 3x3 block
   !>
   !>    vjp_f2_rArA(s, t) = w0 * f2_rArB(s, t) + sum_a w1(a) * f3_r_rArB(a, s, t)
   !>                        + sum_a sum_b w2(a, b) * f4_rr_rArB(a, b, s, t)
   !>
   !> with `B = A`, i.e. the three outputs of `svdw_pair_diag_eval` contracted
   !> against the weights.  `w2` is a *general* 3x3, contracted over all nine
   !> entries.  Two nuclear derivations of the collapsed weighted jet on the same
   !> atom; they commute, so only the six `s <= t` entries are computed.
   !>
   !> Why it exists.  The weighted hvp row of `svdw_hvp_vjp_eval` is linear in
   !> the aggregate `ws` and in the atom's own `aw`.  Moving atom A by `v_A`
   !> changes both: `aw{n}(v_A) = -v_A . at{n+1}` enters `ws` as one summand and
   !> the `aw` inputs directly, and that total is this block times `v_A`.  So
   !> the row of atom A for *any* direction is
   !>
   !>    row_A(v) = hvp_vjp . (ws(v) - aw_A(v)) + vjp_f2_rArA . v_A ,
   !>
   !> the `ws` form applied to the other atoms' share plus this block (feeding
   !> the full `ws(v)` would count atom A twice), and the diagonal block of the
   !> weighted mixed nuclear Hessian is this routine alone, where the
   !> row-per-direction spelling calls `svdw_hvp_eval` three times per atom.
   !>
   !> Reads the power sums to order 2 and the atom's tensors to order 4.
   !>
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps1          Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    ps2          Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all
   !>                            active atoms; last index selects the kind (1 = p1, 2 = p2, 3 = p3,
   !>                            4 = q)
   !> @param[in]    at0          Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at1          Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at2          Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at3          Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    at4          Atom-A kind tensor d^4 u_A^m / dr^4 from svdw_atom_eval; last index
   !>                            selects the kind
   !> @param[in]    w0           Adjoint weight of the level-set value S
   !> @param[in]    w1           Adjoint weights of the spatial gradient dS/dr
   !> @param[in]    w2           Adjoint weights of the spatial Hessian d^2S/dr^2; a general 3x3,
   !>                            contracted over all nine entries
   !> @param[out]   vjp_f2_rArA  w0 f2_rArA + w1 . f3_r_rArA + w2 : f4_rr_rArA of this atom; index
   !>                            order (nuclear s, nuclear t), symmetric
   pure module subroutine svdw_vjp_diag_eval(blend_k, s_1, s_2, s_3, ps0, ps1, ps2, at0, at1, at2, at3, &
                                      at4, w0, w1, w2, vjp_f2_rArA)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Accumulated power-sum tensor d^0 p_kind / dr^0, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(nkind)
      !> Accumulated power-sum tensor d^1 p_kind / dr^1, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(3, nkind)
      !> Accumulated power-sum tensor d^2 p_kind / dr^2, summed over all active atoms; last index
      !> selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(3, 3, nkind)
      !> Atom-A kind tensor d^0 u_A^m / dr^0 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at0(nkind)
      !> Atom-A kind tensor d^1 u_A^m / dr^1 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at1(3, nkind)
      !> Atom-A kind tensor d^2 u_A^m / dr^2 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at2(3, 3, nkind)
      !> Atom-A kind tensor d^3 u_A^m / dr^3 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at3(3, 3, 3, nkind)
      !> Atom-A kind tensor d^4 u_A^m / dr^4 from svdw_atom_eval; last index selects the kind
      real(wp), intent(in) :: at4(3, 3, 3, 3, nkind)
      !> Adjoint weight of the level-set value S
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient dS/dr
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the spatial Hessian d^2S/dr^2; a general 3x3, contracted over all nine
      !> entries
      real(wp), intent(in) :: w2(3, 3)
      !> w0 f2_rArA + w1 . f3_r_rArA + w2 : f4_rr_rArA of this atom; index order (nuclear s, nuclear
      !> t), symmetric
      real(wp), intent(out) :: vjp_f2_rArA(3, 3)
   end subroutine svdw_vjp_diag_eval

   !> Strip power sums through order 0, atom loop outer, point loop inner
   !>
   !> - accumulates `d^n/dr^n sum_A g_kind(r_p - R_A)`, n = 0..0, for the
   !>   `np <= svdw_strip_width(0) = 32` points of one strip over one shared
   !>   atom list, in list order; only orbit representatives (packed slots)
   !> - atom function `g = (w(d_A) u_A)**m`, `d_A = ||r - R_A|| - R_a`,
   !>   `u_A = exp(-(blend_k/3) d_A)` from the inline exp, `w = 1 - S_4(t)`,
   !>   `t = clamp((d_A - fade_start)/fade_width, 0, 1)`
   !> - per atom against the strip box: plateau (`d_max <= fade_start`) runs the
   !>   plain body, beyond (`d_min >= fade_start + fade_width`) is skipped, else
   !>   the faded body; band at infinity: `fade_start = huge(0.0_wp)` with a
   !>   finite `fade_width` (huge + huge would overflow `fade_start + fade_width`)
   !> - hard screening per point: `u = 0` where `||r - R_A||**2 > bound_sq(A)`,
   !>   the gate's own squared compare
   !> - kind-4 (`q`) slots are exact zeros when `include_q` is false
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not
   !>   handled: their sums are not finite (`1/x` is Inf, no FPE trap assumed) and
   !>   the caller re-evaluates the strip on the guarded per-atom path
   !> - requires `blend_k*max(radii)/3 <= 708`
   !> - geometry cache: the gate-masked, unfaded `x`, `1/x`, `u` of every atom
   !>   and point, the same values the sums saw
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    points       Strip points, point-last
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's `cand_screen` row
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    include_q    Accumulate the kind-4 (`q`) power sums; pass `blend_2b /= 0`
   !> @param[in]    fade_start   Fade onset X_on in d_A; huge(0.0_wp) for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, finite and > 0
   !> @param[out]   ps0          Order-0 power sums of the strip, point-first; last index the kind
   !>                            (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   n_zero       Point-atom pairs at zero distance whose geometry was evaluated;
   !>                            nonzero means the sums are not finite
   !> @param[out]   cache_x      Distance ||r_p - R_A|| per point and atom (optional)
   !> @param[out]   cache_inv_x  Reciprocal distance per point and atom (optional)
   !> @param[out]   cache_u      Gate-masked unfaded screening factor u_A per point and atom
   !>                            (optional)
   pure module subroutine svdw_strip_powersums_o0(np, points, nat, centers, radii, bound_sq, blend_k, &
                                           include_q, fade_start, fade_width, ps0, n_zero, &
                                           cache_x, cache_inv_x, cache_u)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Strip points, point-last
      real(wp), intent(in) :: points(3, np)
      !> Atoms in the shared list
      integer, intent(in) :: nat
      !> Atom centers of the shared list
      real(wp), intent(in) :: centers(3, nat)
      !> Atom radii R_a
      real(wp), intent(in) :: radii(nat)
      !> Squared screening reach per atom, the gate's `cand_screen` row
      real(wp), intent(in) :: bound_sq(nat)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Accumulate the kind-4 (`q`) power sums; pass `blend_2b /= 0`
      logical, intent(in) :: include_q
      !> Fade onset X_on in d_A; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, finite and > 0
      real(wp), intent(in) :: fade_width
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 = p2, 3 = p3,
      !> 4 = q)
      real(wp), intent(out) :: ps0(np, nkind)
      !> Point-atom pairs at zero distance whose geometry was evaluated; nonzero means the sums are
      !> not finite
      integer, intent(out) :: n_zero
      !> Distance ||r_p - R_A|| per point and atom (optional)
      real(wp), optional, intent(out) :: cache_x(np, nat)
      !> Reciprocal distance per point and atom (optional)
      real(wp), optional, intent(out) :: cache_inv_x(np, nat)
      !> Gate-masked unfaded screening factor u_A per point and atom (optional)
      real(wp), optional, intent(out) :: cache_u(np, nat)
   end subroutine svdw_strip_powersums_o0

   !> Strip lift through order 0: packed spatial jet of S from strip power sums
   !>
   !> - same algebra as `svdw_spatial_eval` at `max_deriv = 0`; writes `lsf0` point-last, symmetric
   !>   groups packed (sorted tuples, last index slowest)
   !> - `np <= svdw_strip_width(0) = 32`; requires `Z > 0` for every point (escape is the caller's
   !>   test on the power sums)
   !> - weights used as given: no low-atom-count guard
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Order-0 power sums of the strip, point-first; last index the kind
   !>                            (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   lsf0         Level-set value S = -log(Z)/blend_k
   pure module subroutine svdw_strip_lift_o0(np, blend_k, s_1, s_2, s_3, ps0, lsf0)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 = p2, 3 = p3,
      !> 4 = q)
      real(wp), intent(in) :: ps0(np, nkind)
      !> Level-set value S = -log(Z)/blend_k
      real(wp), intent(out) :: lsf0(np)
   end subroutine svdw_strip_lift_o0

   !> Strip power sums through order 1, atom loop outer, point loop inner
   !>
   !> - accumulates `d^n/dr^n sum_A g_kind(r_p - R_A)`, n = 0..1, for the
   !>   `np <= svdw_strip_width(1) = 32` points of one strip over one shared
   !>   atom list, in list order; only orbit representatives (packed slots)
   !> - atom function `g = (w(d_A) u_A)**m`, `d_A = ||r - R_A|| - R_a`,
   !>   `u_A = exp(-(blend_k/3) d_A)` from the inline exp, `w = 1 - S_4(t)`,
   !>   `t = clamp((d_A - fade_start)/fade_width, 0, 1)`
   !> - per atom against the strip box: plateau (`d_max <= fade_start`) runs the
   !>   plain body, beyond (`d_min >= fade_start + fade_width`) is skipped, else
   !>   the faded body; band at infinity: `fade_start = huge(0.0_wp)` with a
   !>   finite `fade_width` (huge + huge would overflow `fade_start + fade_width`)
   !> - hard screening per point: `u = 0` where `||r - R_A||**2 > bound_sq(A)`,
   !>   the gate's own squared compare
   !> - kind-4 (`q`) slots are exact zeros when `include_q` is false
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not
   !>   handled: their sums are not finite (`1/x` is Inf, no FPE trap assumed) and
   !>   the caller re-evaluates the strip on the guarded per-atom path
   !> - requires `blend_k*max(radii)/3 <= 708`
   !> - geometry cache: the gate-masked, unfaded `x`, `1/x`, `u` of every atom
   !>   and point, the same values the sums saw
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    points       Strip points, point-last
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's `cand_screen` row
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    include_q    Accumulate the kind-4 (`q`) power sums; pass `blend_2b /= 0`
   !> @param[in]    fade_start   Fade onset X_on in d_A; huge(0.0_wp) for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, finite and > 0
   !> @param[out]   ps0          Order-0 power sums of the strip, point-first; last index the kind
   !>                            (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   ps1          Order-1 power sums of the strip, point-first; middle index the
   !>                            Cartesian component; last index the kind (1 = p1, 2 = p2, 3 = p3, 4
   !>                            = q)
   !> @param[out]   n_zero       Point-atom pairs at zero distance whose geometry was evaluated;
   !>                            nonzero means the sums are not finite
   !> @param[out]   cache_x      Distance ||r_p - R_A|| per point and atom (optional)
   !> @param[out]   cache_inv_x  Reciprocal distance per point and atom (optional)
   !> @param[out]   cache_u      Gate-masked unfaded screening factor u_A per point and atom
   !>                            (optional)
   pure module subroutine svdw_strip_powersums_o1(np, points, nat, centers, radii, bound_sq, blend_k, &
                                           include_q, fade_start, fade_width, ps0, ps1, n_zero, &
                                           cache_x, cache_inv_x, cache_u)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Strip points, point-last
      real(wp), intent(in) :: points(3, np)
      !> Atoms in the shared list
      integer, intent(in) :: nat
      !> Atom centers of the shared list
      real(wp), intent(in) :: centers(3, nat)
      !> Atom radii R_a
      real(wp), intent(in) :: radii(nat)
      !> Squared screening reach per atom, the gate's `cand_screen` row
      real(wp), intent(in) :: bound_sq(nat)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Accumulate the kind-4 (`q`) power sums; pass `blend_2b /= 0`
      logical, intent(in) :: include_q
      !> Fade onset X_on in d_A; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, finite and > 0
      real(wp), intent(in) :: fade_width
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 = p2, 3 = p3,
      !> 4 = q)
      real(wp), intent(out) :: ps0(np, nkind)
      !> Order-1 power sums of the strip, point-first; middle index the Cartesian component; last
      !> index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(out) :: ps1(np, 3, nkind)
      !> Point-atom pairs at zero distance whose geometry was evaluated; nonzero means the sums are
      !> not finite
      integer, intent(out) :: n_zero
      !> Distance ||r_p - R_A|| per point and atom (optional)
      real(wp), optional, intent(out) :: cache_x(np, nat)
      !> Reciprocal distance per point and atom (optional)
      real(wp), optional, intent(out) :: cache_inv_x(np, nat)
      !> Gate-masked unfaded screening factor u_A per point and atom (optional)
      real(wp), optional, intent(out) :: cache_u(np, nat)
   end subroutine svdw_strip_powersums_o1

   !> Strip lift through order 1: packed spatial jet of S from strip power sums
   !>
   !> - same algebra as `svdw_spatial_eval` at `max_deriv = 1`; writes `lsf0`, `lsf1_r` point-last,
   !>   symmetric groups packed (sorted tuples, last index slowest)
   !> - `np <= svdw_strip_width(1) = 32`; requires `Z > 0` for every point (escape is the caller's
   !>   test on the power sums)
   !> - weights used as given: no low-atom-count guard
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Order-0 power sums of the strip, point-first; last index the kind
   !>                            (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[in]    ps1          Order-1 power sums of the strip, point-first; middle index the
   !>                            Cartesian component; last index the kind (1 = p1, 2 = p2, 3 = p3, 4
   !>                            = q)
   !> @param[out]   lsf0         Level-set value S = -log(Z)/blend_k
   !> @param[out]   lsf1_r       Spatial gradient dS/dr
   pure module subroutine svdw_strip_lift_o1(np, blend_k, s_1, s_2, s_3, ps0, ps1, lsf0, lsf1_r)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 = p2, 3 = p3,
      !> 4 = q)
      real(wp), intent(in) :: ps0(np, nkind)
      !> Order-1 power sums of the strip, point-first; middle index the Cartesian component; last
      !> index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(np, 3, nkind)
      !> Level-set value S = -log(Z)/blend_k
      real(wp), intent(out) :: lsf0(np)
      !> Spatial gradient dS/dr
      real(wp), intent(out) :: lsf1_r(3, np)
   end subroutine svdw_strip_lift_o1

   !> Strip power sums through order 2, atom loop outer, point loop inner
   !>
   !> - accumulates `d^n/dr^n sum_A g_kind(r_p - R_A)`, n = 0..2, for the
   !>   `np <= svdw_strip_width(2) = 32` points of one strip over one shared
   !>   atom list, in list order; only orbit representatives (packed slots)
   !> - atom function `g = (w(d_A) u_A)**m`, `d_A = ||r - R_A|| - R_a`,
   !>   `u_A = exp(-(blend_k/3) d_A)` from the inline exp, `w = 1 - S_4(t)`,
   !>   `t = clamp((d_A - fade_start)/fade_width, 0, 1)`
   !> - per atom against the strip box: plateau (`d_max <= fade_start`) runs the
   !>   plain body, beyond (`d_min >= fade_start + fade_width`) is skipped, else
   !>   the faded body; band at infinity: `fade_start = huge(0.0_wp)` with a
   !>   finite `fade_width` (huge + huge would overflow `fade_start + fade_width`)
   !> - hard screening per point: `u = 0` where `||r - R_A||**2 > bound_sq(A)`,
   !>   the gate's own squared compare
   !> - kind-4 (`q`) slots are exact zeros when `include_q` is false
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not
   !>   handled: their sums are not finite (`1/x` is Inf, no FPE trap assumed) and
   !>   the caller re-evaluates the strip on the guarded per-atom path
   !> - requires `blend_k*max(radii)/3 <= 708`
   !> - geometry cache: the gate-masked, unfaded `x`, `1/x`, `u` of every atom
   !>   and point, the same values the sums saw
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    points       Strip points, point-last
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's `cand_screen` row
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    include_q    Accumulate the kind-4 (`q`) power sums; pass `blend_2b /= 0`
   !> @param[in]    fade_start   Fade onset X_on in d_A; huge(0.0_wp) for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, finite and > 0
   !> @param[out]   ps0          Order-0 power sums of the strip, point-first; last index the kind
   !>                            (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   ps1          Order-1 power sums of the strip, point-first; middle index the
   !>                            Cartesian component; last index the kind (1 = p1, 2 = p2, 3 = p3, 4
   !>                            = q)
   !> @param[out]   ps2          Order-2 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-2 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   n_zero       Point-atom pairs at zero distance whose geometry was evaluated;
   !>                            nonzero means the sums are not finite
   !> @param[out]   cache_x      Distance ||r_p - R_A|| per point and atom (optional)
   !> @param[out]   cache_inv_x  Reciprocal distance per point and atom (optional)
   !> @param[out]   cache_u      Gate-masked unfaded screening factor u_A per point and atom
   !>                            (optional)
   pure module subroutine svdw_strip_powersums_o2(np, points, nat, centers, radii, bound_sq, blend_k, &
                                           include_q, fade_start, fade_width, ps0, ps1, ps2, &
                                           n_zero, cache_x, cache_inv_x, cache_u)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Strip points, point-last
      real(wp), intent(in) :: points(3, np)
      !> Atoms in the shared list
      integer, intent(in) :: nat
      !> Atom centers of the shared list
      real(wp), intent(in) :: centers(3, nat)
      !> Atom radii R_a
      real(wp), intent(in) :: radii(nat)
      !> Squared screening reach per atom, the gate's `cand_screen` row
      real(wp), intent(in) :: bound_sq(nat)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Accumulate the kind-4 (`q`) power sums; pass `blend_2b /= 0`
      logical, intent(in) :: include_q
      !> Fade onset X_on in d_A; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, finite and > 0
      real(wp), intent(in) :: fade_width
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 = p2, 3 = p3,
      !> 4 = q)
      real(wp), intent(out) :: ps0(np, nkind)
      !> Order-1 power sums of the strip, point-first; middle index the Cartesian component; last
      !> index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(out) :: ps1(np, 3, nkind)
      !> Order-2 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-2 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(out) :: ps2(np, 6, nkind)
      !> Point-atom pairs at zero distance whose geometry was evaluated; nonzero means the sums are
      !> not finite
      integer, intent(out) :: n_zero
      !> Distance ||r_p - R_A|| per point and atom (optional)
      real(wp), optional, intent(out) :: cache_x(np, nat)
      !> Reciprocal distance per point and atom (optional)
      real(wp), optional, intent(out) :: cache_inv_x(np, nat)
      !> Gate-masked unfaded screening factor u_A per point and atom (optional)
      real(wp), optional, intent(out) :: cache_u(np, nat)
   end subroutine svdw_strip_powersums_o2

   !> Strip lift through order 2: packed spatial jet of S from strip power sums
   !>
   !> - same algebra as `svdw_spatial_eval` at `max_deriv = 2`; writes `lsf0`, `lsf1_r`, `lsf2_rr`
   !>   point-last, symmetric groups packed (sorted tuples, last index slowest)
   !> - `np <= svdw_strip_width(2) = 32`; requires `Z > 0` for every point (escape is the caller's
   !>   test on the power sums)
   !> - weights used as given: no low-atom-count guard
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Order-0 power sums of the strip, point-first; last index the kind
   !>                            (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[in]    ps1          Order-1 power sums of the strip, point-first; middle index the
   !>                            Cartesian component; last index the kind (1 = p1, 2 = p2, 3 = p3, 4
   !>                            = q)
   !> @param[in]    ps2          Order-2 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-2 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   lsf0         Level-set value S = -log(Z)/blend_k
   !> @param[out]   lsf1_r       Spatial gradient dS/dr
   !> @param[out]   lsf2_rr      Spatial Hessian, packed (6)
   pure module subroutine svdw_strip_lift_o2(np, blend_k, s_1, s_2, s_3, ps0, ps1, ps2, lsf0, lsf1_r, &
                                      lsf2_rr)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 = p2, 3 = p3,
      !> 4 = q)
      real(wp), intent(in) :: ps0(np, nkind)
      !> Order-1 power sums of the strip, point-first; middle index the Cartesian component; last
      !> index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(np, 3, nkind)
      !> Order-2 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-2 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(np, 6, nkind)
      !> Level-set value S = -log(Z)/blend_k
      real(wp), intent(out) :: lsf0(np)
      !> Spatial gradient dS/dr
      real(wp), intent(out) :: lsf1_r(3, np)
      !> Spatial Hessian, packed (6)
      real(wp), intent(out) :: lsf2_rr(6, np)
   end subroutine svdw_strip_lift_o2

   !> Strip power sums through order 3, atom loop outer, point loop inner
   !>
   !> - accumulates `d^n/dr^n sum_A g_kind(r_p - R_A)`, n = 0..3, for the
   !>   `np <= svdw_strip_width(3) = 16` points of one strip over one shared
   !>   atom list, in list order; only orbit representatives (packed slots)
   !> - atom function `g = (w(d_A) u_A)**m`, `d_A = ||r - R_A|| - R_a`,
   !>   `u_A = exp(-(blend_k/3) d_A)` from the inline exp, `w = 1 - S_4(t)`,
   !>   `t = clamp((d_A - fade_start)/fade_width, 0, 1)`
   !> - per atom against the strip box: plateau (`d_max <= fade_start`) runs the
   !>   plain body, beyond (`d_min >= fade_start + fade_width`) is skipped, else
   !>   the faded body; band at infinity: `fade_start = huge(0.0_wp)` with a
   !>   finite `fade_width` (huge + huge would overflow `fade_start + fade_width`)
   !> - hard screening per point: `u = 0` where `||r - R_A||**2 > bound_sq(A)`,
   !>   the gate's own squared compare
   !> - kind-4 (`q`) slots are exact zeros when `include_q` is false
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not
   !>   handled: their sums are not finite (`1/x` is Inf, no FPE trap assumed) and
   !>   the caller re-evaluates the strip on the guarded per-atom path
   !> - requires `blend_k*max(radii)/3 <= 708`
   !> - geometry cache: the gate-masked, unfaded `x`, `1/x`, `u` of every atom
   !>   and point, the same values the sums saw
   !>
   !> @param[in]    np           Live points in the strip, 1..16
   !> @param[in]    points       Strip points, point-last
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's `cand_screen` row
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    include_q    Accumulate the kind-4 (`q`) power sums; pass `blend_2b /= 0`
   !> @param[in]    fade_start   Fade onset X_on in d_A; huge(0.0_wp) for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, finite and > 0
   !> @param[out]   ps0          Order-0 power sums of the strip, point-first; last index the kind
   !>                            (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   ps1          Order-1 power sums of the strip, point-first; middle index the
   !>                            Cartesian component; last index the kind (1 = p1, 2 = p2, 3 = p3, 4
   !>                            = q)
   !> @param[out]   ps2          Order-2 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-2 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   ps3          Order-3 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-3 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   n_zero       Point-atom pairs at zero distance whose geometry was evaluated;
   !>                            nonzero means the sums are not finite
   !> @param[out]   cache_x      Distance ||r_p - R_A|| per point and atom (optional)
   !> @param[out]   cache_inv_x  Reciprocal distance per point and atom (optional)
   !> @param[out]   cache_u      Gate-masked unfaded screening factor u_A per point and atom
   !>                            (optional)
   pure module subroutine svdw_strip_powersums_o3(np, points, nat, centers, radii, bound_sq, blend_k, &
                                           include_q, fade_start, fade_width, ps0, ps1, ps2, ps3, &
                                           n_zero, cache_x, cache_inv_x, cache_u)
      !> Live points in the strip, 1..16
      integer, intent(in) :: np
      !> Strip points, point-last
      real(wp), intent(in) :: points(3, np)
      !> Atoms in the shared list
      integer, intent(in) :: nat
      !> Atom centers of the shared list
      real(wp), intent(in) :: centers(3, nat)
      !> Atom radii R_a
      real(wp), intent(in) :: radii(nat)
      !> Squared screening reach per atom, the gate's `cand_screen` row
      real(wp), intent(in) :: bound_sq(nat)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Accumulate the kind-4 (`q`) power sums; pass `blend_2b /= 0`
      logical, intent(in) :: include_q
      !> Fade onset X_on in d_A; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, finite and > 0
      real(wp), intent(in) :: fade_width
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 = p2, 3 = p3,
      !> 4 = q)
      real(wp), intent(out) :: ps0(np, nkind)
      !> Order-1 power sums of the strip, point-first; middle index the Cartesian component; last
      !> index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(out) :: ps1(np, 3, nkind)
      !> Order-2 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-2 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(out) :: ps2(np, 6, nkind)
      !> Order-3 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-3 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(out) :: ps3(np, 10, nkind)
      !> Point-atom pairs at zero distance whose geometry was evaluated; nonzero means the sums are
      !> not finite
      integer, intent(out) :: n_zero
      !> Distance ||r_p - R_A|| per point and atom (optional)
      real(wp), optional, intent(out) :: cache_x(np, nat)
      !> Reciprocal distance per point and atom (optional)
      real(wp), optional, intent(out) :: cache_inv_x(np, nat)
      !> Gate-masked unfaded screening factor u_A per point and atom (optional)
      real(wp), optional, intent(out) :: cache_u(np, nat)
   end subroutine svdw_strip_powersums_o3

   !> Strip lift through order 3: packed spatial jet of S from strip power sums
   !>
   !> - same algebra as `svdw_spatial_eval` at `max_deriv = 3`; writes `lsf0`, `lsf1_r`, `lsf2_rr`,
   !>   `lsf3_rrr` point-last, symmetric groups packed (sorted tuples, last index slowest)
   !> - `np <= svdw_strip_width(3) = 16`; requires `Z > 0` for every point (escape is the caller's
   !>   test on the power sums)
   !> - weights used as given: no low-atom-count guard
   !>
   !> @param[in]    np           Live points in the strip, 1..16
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Order-0 power sums of the strip, point-first; last index the kind
   !>                            (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[in]    ps1          Order-1 power sums of the strip, point-first; middle index the
   !>                            Cartesian component; last index the kind (1 = p1, 2 = p2, 3 = p3, 4
   !>                            = q)
   !> @param[in]    ps2          Order-2 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-2 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[in]    ps3          Order-3 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-3 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   lsf0         Level-set value S = -log(Z)/blend_k
   !> @param[out]   lsf1_r       Spatial gradient dS/dr
   !> @param[out]   lsf2_rr      Spatial Hessian, packed (6)
   !> @param[out]   lsf3_rrr     Third spatial derivative, packed (10)
   pure module subroutine svdw_strip_lift_o3(np, blend_k, s_1, s_2, s_3, ps0, ps1, ps2, ps3, lsf0, &
                                      lsf1_r, lsf2_rr, lsf3_rrr)
      !> Live points in the strip, 1..16
      integer, intent(in) :: np
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 = p2, 3 = p3,
      !> 4 = q)
      real(wp), intent(in) :: ps0(np, nkind)
      !> Order-1 power sums of the strip, point-first; middle index the Cartesian component; last
      !> index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(np, 3, nkind)
      !> Order-2 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-2 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(np, 6, nkind)
      !> Order-3 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-3 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps3(np, 10, nkind)
      !> Level-set value S = -log(Z)/blend_k
      real(wp), intent(out) :: lsf0(np)
      !> Spatial gradient dS/dr
      real(wp), intent(out) :: lsf1_r(3, np)
      !> Spatial Hessian, packed (6)
      real(wp), intent(out) :: lsf2_rr(6, np)
      !> Third spatial derivative, packed (10)
      real(wp), intent(out) :: lsf3_rrr(10, np)
   end subroutine svdw_strip_lift_o3

   !> Strip power sums through order 4, atom loop outer, point loop inner
   !>
   !> - accumulates `d^n/dr^n sum_A g_kind(r_p - R_A)`, n = 0..4, for the
   !>   `np <= svdw_strip_width(4) = 8` points of one strip over one shared
   !>   atom list, in list order; only orbit representatives (packed slots)
   !> - atom function `g = (w(d_A) u_A)**m`, `d_A = ||r - R_A|| - R_a`,
   !>   `u_A = exp(-(blend_k/3) d_A)` from the inline exp, `w = 1 - S_4(t)`,
   !>   `t = clamp((d_A - fade_start)/fade_width, 0, 1)`
   !> - per atom against the strip box: plateau (`d_max <= fade_start`) runs the
   !>   plain body, beyond (`d_min >= fade_start + fade_width`) is skipped, else
   !>   the faded body; band at infinity: `fade_start = huge(0.0_wp)` with a
   !>   finite `fade_width` (huge + huge would overflow `fade_start + fade_width`)
   !> - hard screening per point: `u = 0` where `||r - R_A||**2 > bound_sq(A)`,
   !>   the gate's own squared compare
   !> - kind-4 (`q`) slots are exact zeros when `include_q` is false
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not
   !>   handled: their sums are not finite (`1/x` is Inf, no FPE trap assumed) and
   !>   the caller re-evaluates the strip on the guarded per-atom path
   !> - requires `blend_k*max(radii)/3 <= 708`
   !> - geometry cache: the gate-masked, unfaded `x`, `1/x`, `u` of every atom
   !>   and point, the same values the sums saw
   !>
   !> @param[in]    np           Live points in the strip, 1..8
   !> @param[in]    points       Strip points, point-last
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's `cand_screen` row
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    include_q    Accumulate the kind-4 (`q`) power sums; pass `blend_2b /= 0`
   !> @param[in]    fade_start   Fade onset X_on in d_A; huge(0.0_wp) for no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, finite and > 0
   !> @param[out]   ps0          Order-0 power sums of the strip, point-first; last index the kind
   !>                            (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   ps1          Order-1 power sums of the strip, point-first; middle index the
   !>                            Cartesian component; last index the kind (1 = p1, 2 = p2, 3 = p3, 4
   !>                            = q)
   !> @param[out]   ps2          Order-2 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-2 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   ps3          Order-3 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-3 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   ps4          Order-4 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-4 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   n_zero       Point-atom pairs at zero distance whose geometry was evaluated;
   !>                            nonzero means the sums are not finite
   !> @param[out]   cache_x      Distance ||r_p - R_A|| per point and atom (optional)
   !> @param[out]   cache_inv_x  Reciprocal distance per point and atom (optional)
   !> @param[out]   cache_u      Gate-masked unfaded screening factor u_A per point and atom
   !>                            (optional)
   pure module subroutine svdw_strip_powersums_o4(np, points, nat, centers, radii, bound_sq, blend_k, &
                                           include_q, fade_start, fade_width, ps0, ps1, ps2, ps3, &
                                           ps4, n_zero, cache_x, cache_inv_x, cache_u)
      !> Live points in the strip, 1..8
      integer, intent(in) :: np
      !> Strip points, point-last
      real(wp), intent(in) :: points(3, np)
      !> Atoms in the shared list
      integer, intent(in) :: nat
      !> Atom centers of the shared list
      real(wp), intent(in) :: centers(3, nat)
      !> Atom radii R_a
      real(wp), intent(in) :: radii(nat)
      !> Squared screening reach per atom, the gate's `cand_screen` row
      real(wp), intent(in) :: bound_sq(nat)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Accumulate the kind-4 (`q`) power sums; pass `blend_2b /= 0`
      logical, intent(in) :: include_q
      !> Fade onset X_on in d_A; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, finite and > 0
      real(wp), intent(in) :: fade_width
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 = p2, 3 = p3,
      !> 4 = q)
      real(wp), intent(out) :: ps0(np, nkind)
      !> Order-1 power sums of the strip, point-first; middle index the Cartesian component; last
      !> index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(out) :: ps1(np, 3, nkind)
      !> Order-2 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-2 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(out) :: ps2(np, 6, nkind)
      !> Order-3 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-3 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(out) :: ps3(np, 10, nkind)
      !> Order-4 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-4 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(out) :: ps4(np, 15, nkind)
      !> Point-atom pairs at zero distance whose geometry was evaluated; nonzero means the sums are
      !> not finite
      integer, intent(out) :: n_zero
      !> Distance ||r_p - R_A|| per point and atom (optional)
      real(wp), optional, intent(out) :: cache_x(np, nat)
      !> Reciprocal distance per point and atom (optional)
      real(wp), optional, intent(out) :: cache_inv_x(np, nat)
      !> Gate-masked unfaded screening factor u_A per point and atom (optional)
      real(wp), optional, intent(out) :: cache_u(np, nat)
   end subroutine svdw_strip_powersums_o4

   !> Strip lift through order 4: packed spatial jet of S from strip power sums
   !>
   !> - same algebra as `svdw_spatial_eval` at `max_deriv = 4`; writes `lsf0`, `lsf1_r`, `lsf2_rr`,
   !>   `lsf3_rrr`, `lsf4_rrrr` point-last, symmetric groups packed (sorted tuples, last index
   !>   slowest)
   !> - `np <= svdw_strip_width(4) = 8`; requires `Z > 0` for every point (escape is the caller's
   !>   test on the power sums)
   !> - weights used as given: no low-atom-count guard
   !>
   !> @param[in]    np           Live points in the strip, 1..8
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Order-0 power sums of the strip, point-first; last index the kind
   !>                            (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[in]    ps1          Order-1 power sums of the strip, point-first; middle index the
   !>                            Cartesian component; last index the kind (1 = p1, 2 = p2, 3 = p3, 4
   !>                            = q)
   !> @param[in]    ps2          Order-2 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-2 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[in]    ps3          Order-3 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-3 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[in]    ps4          Order-4 power sums of the strip, point-first; middle index the
   !>                            packed slot of the sorted rank-4 tuple; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[out]   lsf0         Level-set value S = -log(Z)/blend_k
   !> @param[out]   lsf1_r       Spatial gradient dS/dr
   !> @param[out]   lsf2_rr      Spatial Hessian, packed (6)
   !> @param[out]   lsf3_rrr     Third spatial derivative, packed (10)
   !> @param[out]   lsf4_rrrr    Fourth spatial derivative, packed (15)
   pure module subroutine svdw_strip_lift_o4(np, blend_k, s_1, s_2, s_3, ps0, ps1, ps2, ps3, ps4, lsf0, &
                                      lsf1_r, lsf2_rr, lsf3_rrr, lsf4_rrrr)
      !> Live points in the strip, 1..8
      integer, intent(in) :: np
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 = p2, 3 = p3,
      !> 4 = q)
      real(wp), intent(in) :: ps0(np, nkind)
      !> Order-1 power sums of the strip, point-first; middle index the Cartesian component; last
      !> index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(np, 3, nkind)
      !> Order-2 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-2 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(np, 6, nkind)
      !> Order-3 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-3 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps3(np, 10, nkind)
      !> Order-4 power sums of the strip, point-first; middle index the packed slot of the sorted
      !> rank-4 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps4(np, 15, nkind)
      !> Level-set value S = -log(Z)/blend_k
      real(wp), intent(out) :: lsf0(np)
      !> Spatial gradient dS/dr
      real(wp), intent(out) :: lsf1_r(3, np)
      !> Spatial Hessian, packed (6)
      real(wp), intent(out) :: lsf2_rr(6, np)
      !> Third spatial derivative, packed (10)
      real(wp), intent(out) :: lsf3_rrr(10, np)
      !> Fourth spatial derivative, packed (15)
      real(wp), intent(out) :: lsf4_rrrr(15, np)
   end subroutine svdw_strip_lift_o4

   !> Strip adjoint lift: lambda of the packed weighted jet w.r.t. the order-2 power sums
   !>
   !> - per point, `lam = d/d ps` of `w0 S + w1 . dS/dr + W : d2S/dr2` with `S =
   !>   -log(Z)/blend_k` and `Z` the power-sum polynomial; the same algebra as
   !>   `svdw_strip_lift_o2`, differentiated once in the power sums
   !> - reads `ps0`, `ps1`, `ps2` of `svdw_strip_powersums_o2` only; `lam` has their
   !>   layout, one entry per orbit representative and kind
   !> - seeds hold the symmetric part of the (3, 3) Hessian seed packed as `ps2` (xx, xy,
   !>   yy, xz, yz, zz), entries not doubled; the factor 2 of the off-diagonal pairs is
   !>   applied here
   !> - feeds `svdw_strip_vjp`; a `lam` row of a kind with a zero blending weight (`q` at
   !>   `blend_2b = 0`) is exactly zero
   !> - `np <= svdw_strip_width(2) = 32`; requires `Z > 0` for every point (escape is the
   !>   caller's test on the power sums)
   !> - weights used as given: no low-atom-count guard
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    s_1          One-body blending weight blend_1b
   !> @param[in]    s_2          Two-body blending weight blend_2b
   !> @param[in]    s_3          Three-body blending weight blend_3b
   !> @param[in]    ps0          Order-0 power sums of the strip, point-first; last index
   !>                            the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[in]    ps1          Order-1 power sums of the strip, point-first; middle
   !>                            index the Cartesian component; last index the kind (1 =
   !>                            p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[in]    ps2          Order-2 power sums of the strip, point-first; middle
   !>                            index the packed slot of the sorted rank-2 tuple; last
   !>                            index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
   !> @param[in]    w0           Seed of the level-set value S
   !> @param[in]    w1           Seed of the spatial gradient dS/dr, point-last
   !> @param[in]    w2           Symmetric part of the (3, 3) seed of d2S/dr2, packed (xx,
   !>                            xy, yy, xz, yz, zz), point-last; entries not doubled
   !> @param[out]   lam0         Adjoint of the order-0 power sums, layout of `ps0`
   !> @param[out]   lam1         Adjoint of the order-1 power sums, layout of `ps1`
   !> @param[out]   lam2         Adjoint of the order-2 power sums, layout of `ps2`
   pure module subroutine svdw_strip_lift_vjp(np, blend_k, s_1, s_2, s_3, ps0, ps1, ps2, w0, &
                                       w1, w2, lam0, lam1, lam2)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> One-body blending weight blend_1b
      real(wp), intent(in) :: s_1
      !> Two-body blending weight blend_2b
      real(wp), intent(in) :: s_2
      !> Three-body blending weight blend_3b
      real(wp), intent(in) :: s_3
      !> Order-0 power sums of the strip, point-first; last index the kind (1 = p1, 2 =
      !> p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps0(np, nkind)
      !> Order-1 power sums of the strip, point-first; middle index the Cartesian
      !> component; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps1(np, 3, nkind)
      !> Order-2 power sums of the strip, point-first; middle index the packed slot of the
      !> sorted rank-2 tuple; last index the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp), intent(in) :: ps2(np, 6, nkind)
      !> Seed of the level-set value S
      real(wp), intent(in) :: w0(np)
      !> Seed of the spatial gradient dS/dr, point-last
      real(wp), intent(in) :: w1(3, np)
      !> Symmetric part of the (3, 3) seed of d2S/dr2, packed (xx, xy, yy, xz, yz, zz),
      !> point-last; entries not doubled
      real(wp), intent(in) :: w2(6, np)
      !> Adjoint of the order-0 power sums, layout of `ps0`
      real(wp), intent(out) :: lam0(np, nkind)
      !> Adjoint of the order-1 power sums, layout of `ps1`
      real(wp), intent(out) :: lam1(np, 3, nkind)
      !> Adjoint of the order-2 power sums, layout of `ps2`
      real(wp), intent(out) :: lam2(np, 6, nkind)
   end subroutine svdw_strip_lift_vjp

   !> Strip nuclear vjp: gradient w.r.t. every atom center of the weighted jet
   !>
   !> - `res(:, a, p) = d/dR_a [w0 S + w1 . dS/dr + W : d2S/dr2](r_p)`, which is
   !>   `-sum_kind [lam0 at1 + lam1 : at2 + lam2 : at3]` contracted over the jet indices,
   !>   `lam` from `svdw_strip_lift_vjp`; `at3` is never formed: the radial chain `at(J) =
   !>   sum_j c_J,j(d, 1/x) g^(j)(d_A)` is contracted with `lam` first and the four kinds
   !>   are summed on the radial derivatives `g^(j)`
   !> - reads the geometry cache of `svdw_strip_powersums_o2` for the same strip
   !>   (`cache_x`, `cache_inv_x`, `cache_u`, all `(np, nat)`), so there is no second exp;
   !>   `points`, `centers`, `radii`, `blend_k`, `include_q`, `fade_start`, `fade_width`
   !>   must be the stage-1 arguments
   !> - atom function `g = (w(d_A) u_A)**m` and the classification against the strip box
   !>   as in stage 1: plateau (`d_max <= fade_start`) runs the plain body, beyond the
   !>   band (`d_min >= fade_start + fade_width`) writes exact zeros, else the faded body;
   !>   the `q` kind is skipped when `include_q` is false
   !> - on-nucleus pairs (`cache_x == 0`) are counted in `n_zero`, not handled:
   !>   `cache_inv_x` is Inf there and the `res` lane is not finite; the caller
   !>   re-evaluates such a strip on the guarded per-atom path
   !> - a point whose stage-1 column was unreached (`Z == 0`) or failed has no `lam`; its
   !>   `res` column is the caller's to discard, nothing here tests it
   !> - atoms beyond a point's reach inside a strip that meets the band get exact zeros
   !>   from the fade (`s = 0`) or from the gate-masked `cache_u`
   !> - atom-outer, point-inner; the point loops carry no branch and no call; scratch is
   !>   point-first, `res` is written point-last
   !> - `np <= svdw_strip_width(2) = 32`, the width stage 1 runs at in this pipeline
   !>   (order 2)
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    points       Strip points, point-last; those of stage 1
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    blend_k      Blending sharpness k of u = exp(-(k/3) d)
   !> @param[in]    include_q    Contract the kind-4 (`q`) adjoints; pass `blend_2b /= 0`
   !> @param[in]    fade_start   Fade onset X_on in d_A, as in stage 1; huge(0.0_wp) for
   !>                            no fade
   !> @param[in]    fade_width   Fade band width Delta in d_A, as in stage 1; > 0
   !> @param[in]    lam0         Adjoint of the order-0 power sums from
   !>                            `svdw_strip_lift_vjp`
   !> @param[in]    lam1         Adjoint of the order-1 power sums from
   !>                            `svdw_strip_lift_vjp`
   !> @param[in]    lam2         Adjoint of the order-2 power sums from
   !>                            `svdw_strip_lift_vjp`
   !> @param[in]    cache_x      Distance ||r_p - R_A|| per point and atom, stage-1 cache
   !> @param[in]    cache_inv_x  Reciprocal distance per point and atom, stage-1 cache
   !> @param[in]    cache_u      Gate-masked unfaded screening factor u_A per point and
   !>                            atom, stage-1 cache
   !> @param[out]   res          Nuclear gradient of the weighted jet per atom and point;
   !>                            exact zeros out of reach
   !> @param[out]   n_zero       Point-atom pairs at zero distance; nonzero means `res` is
   !>                            not finite there
   pure module subroutine svdw_strip_vjp(np, points, nat, centers, radii, blend_k, include_q, &
                                  fade_start, fade_width, lam0, lam1, lam2, cache_x, &
                                  cache_inv_x, cache_u, res, n_zero)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Strip points, point-last; those of stage 1
      real(wp), intent(in) :: points(3, np)
      !> Atoms in the shared list
      integer, intent(in) :: nat
      !> Atom centers of the shared list
      real(wp), intent(in) :: centers(3, nat)
      !> Atom radii R_a
      real(wp), intent(in) :: radii(nat)
      !> Blending sharpness k of u = exp(-(k/3) d)
      real(wp), intent(in) :: blend_k
      !> Contract the kind-4 (`q`) adjoints; pass `blend_2b /= 0`
      logical, intent(in) :: include_q
      !> Fade onset X_on in d_A, as in stage 1; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width Delta in d_A, as in stage 1; > 0
      real(wp), intent(in) :: fade_width
      !> Adjoint of the order-0 power sums from `svdw_strip_lift_vjp`
      real(wp), intent(in) :: lam0(np, nkind)
      !> Adjoint of the order-1 power sums from `svdw_strip_lift_vjp`
      real(wp), intent(in) :: lam1(np, 3, nkind)
      !> Adjoint of the order-2 power sums from `svdw_strip_lift_vjp`
      real(wp), intent(in) :: lam2(np, 6, nkind)
      !> Distance ||r_p - R_A|| per point and atom, stage-1 cache
      real(wp), intent(in) :: cache_x(np, nat)
      !> Reciprocal distance per point and atom, stage-1 cache
      real(wp), intent(in) :: cache_inv_x(np, nat)
      !> Gate-masked unfaded screening factor u_A per point and atom, stage-1 cache
      real(wp), intent(in) :: cache_u(np, nat)
      !> Nuclear gradient of the weighted jet per atom and point; exact zeros out of reach
      real(wp), intent(out) :: res(3, nat, np)
      !> Point-atom pairs at zero distance; nonzero means `res` is not finite there
      integer, intent(out) :: n_zero
   end subroutine svdw_strip_vjp

   !> Normalized level set S/||grad S|| and its nuclear gradient, for one atom.
   !>
   !> `norm0 = f0 / ||f1_r||` and
   !>
   !>    norm1_rA(s) = f1_rA(s)/g - f0 * (sum_j f1_r(j) f2_r_rA(j, s)) / g**3 ,
   !>
   !> `g = ||f1_r||`.  All four inputs come from `svdw_spatial_eval` and
   !> `svdw_nuclear_eval`, so this adds no derivative work of its own.  A
   !> vanishing gradient (only reachable at a stationary point of the blend)
   !> returns zeros, matching the hand-written implementation.
   !>
   !> @param[in]  f0        Level-set value S
   !> @param[in]  f1_r      Spatial gradient dS/dr
   !> @param[in]  f1_rA     Nuclear gradient dS/dR_A of this atom
   !> @param[in]  f2_r_rA   Mixed derivative d^2S/(dr dR_A) of this atom
   !> @param[out] norm0     S/||grad S||
   !> @param[out] norm1_rA  d/dR_A of S/||grad S||
   pure module subroutine svdw_normalized_eval(f0, f1_r, f1_rA, f2_r_rA, norm0, norm1_rA)
      !> Level-set value
      real(wp), intent(in) :: f0
      !> Spatial gradient
      real(wp), intent(in) :: f1_r(3)
      !> Nuclear gradient of this atom
      real(wp), intent(in) :: f1_rA(3)
      !> Mixed spatial/nuclear derivative of this atom
      real(wp), intent(in) :: f2_r_rA(3, 3)
      !> Normalized level-set value
      real(wp), intent(out) :: norm0
      !> Nuclear gradient of the normalized level-set value
      real(wp), intent(out) :: norm1_rA(3)
   end subroutine svdw_normalized_eval

!* ================================================================================= *!
!*                            Reference orchestration                                *!
!* ================================================================================= *!
!
! Plain-loop drivers wiring the kernel routines together over an active-atom list.
! They exist so the kernel can be exercised end to end and so a caller has a
! worked example of the intended loop structure; a production caller that already
! owns an SSD/screening system will inline these loops and share one power-sum
! accumulation across several outputs instead of rebuilding it per driver.
!
! Every driver is O(n_active) except `svdw_lsf_pair`, whose *result* is
! O(n_active**2) by construction.
!
! The drivers run the band at infinity with libm `exp`: `svdw_atom_geometry` stands
! in for the strip kernel's geometry cache, so the validation compares the
! generated algebra against the unfaded model independently of the inline exp.

   !> Libm reference of the stage-1 geometry cache of one atom
   !>
   !> @param[in]  d_a      Displacement r - R_A of the atom
   !> @param[in]  R_a      Radius of the atom
   !> @param[in]  blend_k  Blending sharpness k
   !> @param[out] x_a      Distance ||d_a||
   !> @param[out] inv_x_a  1/x_a; huge(1.0_wp) on the nucleus, where the strip cache
   !>                      holds Inf: both undefined there, and only order 0 is read
   !> @param[out] u_a      exp(-(blend_k/3)(x_a - R_a)), libm; not gate-masked, the
   !>                      drivers have no screening gate
   pure module subroutine svdw_atom_geometry(d_a, R_a, blend_k, x_a, inv_x_a, u_a)
      !> Displacement r - R_A of the atom
      real(wp), intent(in) :: d_a(3)
      !> Radius of the atom
      real(wp), intent(in) :: R_a
      !> Blending sharpness k
      real(wp), intent(in) :: blend_k
      !> Distance
      real(wp), intent(out) :: x_a
      !> Reciprocal distance
      real(wp), intent(out) :: inv_x_a
      !> Unfaded screening factor
      real(wp), intent(out) :: u_a
   end subroutine svdw_atom_geometry

   !> Accumulate the direction-contracted power-sum tensors of one point.
   !>
   !> @param[in]  d          Displacements r - R_A, shape (3, n_active)
   !> @param[in]  radii      Radii of the active atoms, shape (n_active)
   !> @param[in]  v          Nuclear displacement directions, shape (3, n_active)
   !> @param[in]  blend_k    Blending sharpness k
   !> @param[in]  max_deriv  Highest spatial-derivative order to accumulate (0..3)
   !> @param[out] ws0        Order-0 contracted tensors
   !> @param[out] ws1        Order-1 contracted tensors
   !> @param[out] ws2        Order-2 contracted tensors
   !> @param[out] ws3        Order-3 contracted tensors
   pure module subroutine svdw_powersums_tangent(d, radii, v, blend_k, max_deriv, &
                                          ws0, ws1, ws2, ws3, vrad)
      !> Displacements r - R_A
      real(wp), intent(in) :: d(:, :)
      !> Radii of the active atoms
      real(wp), intent(in) :: radii(:)
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Blending sharpness k
      real(wp), intent(in) :: blend_k
      !> Highest spatial-derivative order to accumulate
      integer, intent(in) :: max_deriv
      !> Order-0 contracted tensors
      real(wp), intent(out) :: ws0(nkind)
      !> Order-1 contracted tensors
      real(wp), intent(out) :: ws1(3, nkind)
      !> Order-2 contracted tensors
      real(wp), intent(out) :: ws2(3, 3, nkind)
      !> Order-3 contracted tensors
      real(wp), intent(out) :: ws3(3, 3, 3, nkind)
      !> Radius half of the joint direction, vr_A = sum_B (dR_A/dX_B) . v_B.
      !> Present only for a geometry-dependent radius model; when it is absent
      !> not a single floating-point operation below changes, which is what keeps
      !> a static-radius run bit-for-bit on its current results
      real(wp), intent(in), optional :: vrad(:)
   end subroutine svdw_powersums_tangent

   !> Level-set value and pure spatial derivatives over an active-atom list.
   !>
   !> @param[in]  d          Displacements r - R_A, shape (3, n_active)
   !> @param[in]  radii      Radii of the active atoms, shape (n_active)
   !> @param[in]  blend_k    Blending sharpness k
   !> @param[in]  s_1        One-body weight
   !> @param[in]  s_2        Two-body weight
   !> @param[in]  s_3        Three-body weight
   !> @param[in]  max_deriv  Highest spatial-derivative order (0..4)
   !> @param[out] f0         Level-set value
   !> @param[out] f1_r       Spatial gradient
   !> @param[out] f2_rr      Spatial Hessian
   !> @param[out] f3_rrr     Third spatial derivative
   !> @param[out] f4_rrrr    Fourth spatial derivative
   pure module subroutine svdw_lsf_spatial(d, radii, blend_k, s_1, s_2, s_3, max_deriv, &
                                    f0, f1_r, f2_rr, f3_rrr, f4_rrrr)
      !> Displacements r - R_A
      real(wp), intent(in) :: d(:, :)
      !> Radii of the active atoms
      real(wp), intent(in) :: radii(:)
      !> Blending sharpness k
      real(wp), intent(in) :: blend_k
      !> One-body weight
      real(wp), intent(in) :: s_1
      !> Two-body weight
      real(wp), intent(in) :: s_2
      !> Three-body weight
      real(wp), intent(in) :: s_3
      !> Highest spatial-derivative order
      integer, intent(in) :: max_deriv
      !> Level-set value
      real(wp), intent(out) :: f0
      !> Spatial gradient
      real(wp), intent(out) :: f1_r(3)
      !> Spatial Hessian
      real(wp), intent(out) :: f2_rr(3, 3)
      !> Third spatial derivative
      real(wp), intent(out) :: f3_rrr(3, 3, 3)
      !> Fourth spatial derivative
      real(wp), intent(out) :: f4_rrrr(3, 3, 3, 3)
   end subroutine svdw_lsf_spatial

   !> One-nuclear-index derivatives for every active atom.
   !>
   !> @param[in]  d           Displacements r - R_A, shape (3, n_active)
   !> @param[in]  radii       Radii of the active atoms, shape (n_active)
   !> @param[in]  blend_k     Blending sharpness k
   !> @param[in]  s_1         One-body weight
   !> @param[in]  s_2         Two-body weight
   !> @param[in]  s_3         Three-body weight
   !> @param[in]  max_deriv   Highest *spatial* derivative order (0..3)
   !> @param[out] f1_rA       dS/dR_A, shape (3, n_active)
   !> @param[out] f2_r_rA     d^2S/(dr dR_A), shape (3, 3, n_active)
   !> @param[out] f3_rr_rA    d^3S/(dr^2 dR_A), shape (3, 3, 3, n_active)
   !> @param[out] f4_rrr_rA   d^4S/(dr^3 dR_A), shape (3, 3, 3, 3, n_active)
   pure module subroutine svdw_lsf_nuclear(d, radii, blend_k, s_1, s_2, s_3, max_deriv, &
                                    f1_rA, f2_r_rA, f3_rr_rA, f4_rrr_rA)
      !> Displacements r - R_A
      real(wp), intent(in) :: d(:, :)
      !> Radii of the active atoms
      real(wp), intent(in) :: radii(:)
      !> Blending sharpness k
      real(wp), intent(in) :: blend_k
      !> One-body weight
      real(wp), intent(in) :: s_1
      !> Two-body weight
      real(wp), intent(in) :: s_2
      !> Three-body weight
      real(wp), intent(in) :: s_3
      !> Highest spatial-derivative order
      integer, intent(in) :: max_deriv
      !> dS/dR_A
      real(wp), intent(out) :: f1_rA(:, :)
      !> d^2S/(dr dR_A)
      real(wp), intent(out) :: f2_r_rA(:, :, :)
      !> d^3S/(dr^2 dR_A)
      real(wp), intent(out) :: f3_rr_rA(:, :, :, :)
      !> d^4S/(dr^3 dR_A)
      real(wp), intent(out) :: f4_rrr_rA(:, :, :, :, :)
   end subroutine svdw_lsf_nuclear

   !> Radius derivatives over an active-atom list.
   !>
   !> The radius twin of `svdw_lsf_nuclear`, and deliberately identical to it in
   !> shape: `d/dR_a` costs no derivative order, so this driver reads the same
   !> `ps*` with the radius tensors `atd*` in place of `at*`, runs the same single
   !> O(n_active) loop and reaches the same spatial order.  The results carry no
   !> trailing derivative index.
   !>
   !> @param[in]  d          Displacements r - R_A, shape (3, n_active)
   !> @param[in]  radii      Radii of the active atoms
   !> @param[in]  blend_k    Blending sharpness k
   !> @param[in]  s_1        One-body weight
   !> @param[in]  s_2        Two-body weight
   !> @param[in]  s_3        Three-body weight
   !> @param[in]  max_deriv  Highest spatial-derivative order (0..3)
   !> @param[out] f1_rad     dS/dR_a per atom
   !> @param[out] f2_r_rad   d^2S/(dr dR_a) per atom
   !> @param[out] f3_rr_rad  d^3S/(dr^2 dR_a) per atom
   !> @param[out] f4_rrr_rad d^4S/(dr^3 dR_a) per atom
   pure module subroutine svdw_lsf_radius(d, radii, blend_k, s_1, s_2, s_3, max_deriv, &
                                   f1_rad, f2_r_rad, f3_rr_rad, f4_rrr_rad)
      !> Displacements r - R_A
      real(wp), intent(in) :: d(:, :)
      !> Radii of the active atoms
      real(wp), intent(in) :: radii(:)
      !> Blending sharpness k
      real(wp), intent(in) :: blend_k
      !> One-body weight
      real(wp), intent(in) :: s_1
      !> Two-body weight
      real(wp), intent(in) :: s_2
      !> Three-body weight
      real(wp), intent(in) :: s_3
      !> Highest spatial-derivative order
      integer, intent(in) :: max_deriv
      !> dS/dR_a
      real(wp), intent(out) :: f1_rad(:)
      !> d^2S/(dr dR_a)
      real(wp), intent(out) :: f2_r_rad(:, :)
      !> d^3S/(dr^2 dR_a)
      real(wp), intent(out) :: f3_rr_rad(:, :, :)
      !> d^4S/(dr^3 dR_a)
      real(wp), intent(out) :: f4_rrr_rad(:, :, :, :)
   end subroutine svdw_lsf_radius

   !> Two-nuclear-index derivatives for every ordered pair of active atoms.
   !>
   !> The per-atom kind tensors are built once into a local cache, so the pair
   !> loop costs no extra exponentials; the O(n_active**2) here is the result
   !> tensor itself, which is why `svdw_lsf_hvp` exists.
   !>
   !> @param[in]  d             Displacements r - R_A, shape (3, n_active)
   !> @param[in]  radii         Radii of the active atoms, shape (n_active)
   !> @param[in]  blend_k       Blending sharpness k
   !> @param[in]  s_1           One-body weight
   !> @param[in]  s_2           Two-body weight
   !> @param[in]  s_3           Three-body weight
   !> @param[in]  max_deriv     Highest *spatial* derivative order (0..2)
   !> @param[out] f2_rArB       d^2S/(dR_A dR_B), shape (3, n, 3, n)
   !> @param[out] f3_r_rArB     d^3S/(dr dR_A dR_B), shape (3, 3, n, 3, n)
   !> @param[out] f4_rr_rArB    d^4S/(dr^2 dR_A dR_B), shape (3, 3, 3, n, 3, n)
   pure module subroutine svdw_lsf_pair(d, radii, blend_k, s_1, s_2, s_3, max_deriv, &
                                 f2_rArB, f3_r_rArB, f4_rr_rArB)
      !> Displacements r - R_A
      real(wp), intent(in) :: d(:, :)
      !> Radii of the active atoms
      real(wp), intent(in) :: radii(:)
      !> Blending sharpness k
      real(wp), intent(in) :: blend_k
      !> One-body weight
      real(wp), intent(in) :: s_1
      !> Two-body weight
      real(wp), intent(in) :: s_2
      !> Three-body weight
      real(wp), intent(in) :: s_3
      !> Highest spatial-derivative order
      integer, intent(in) :: max_deriv
      !> d^2S/(dR_A dR_B)
      real(wp), intent(out) :: f2_rArB(:, :, :, :)
      !> d^3S/(dr dR_A dR_B)
      real(wp), intent(out) :: f3_r_rArB(:, :, :, :, :)
      !> d^4S/(dr^2 dR_A dR_B)
      real(wp), intent(out) :: f4_rr_rArB(:, :, :, :, :, :)
   end subroutine svdw_lsf_pair

   !> Direction-contracted nuclear derivatives of the spatial tensors.
   !>
   !> @param[in]  d              Displacements r - R_A, shape (3, n_active)
   !> @param[in]  radii          Radii of the active atoms, shape (n_active)
   !> @param[in]  v              Displacement directions, shape (3, n_active)
   !> @param[in]  blend_k        Blending sharpness k
   !> @param[in]  s_1            One-body weight
   !> @param[in]  s_2            Two-body weight
   !> @param[in]  s_3            Three-body weight
   !> @param[in]  max_deriv      Highest spatial-derivative order (0..3)
   !> @param[out] tangent_f0     sum_B v_B . dS/dR_B
   !> @param[out] tangent_f1_r   sum_B v_B . d/dR_B of dS/dr
   !> @param[out] tangent_f2_rr  sum_B v_B . d/dR_B of d^2S/dr^2
   !> @param[out] tangent_f3_rrr sum_B v_B . d/dR_B of d^3S/dr^3
   pure module subroutine svdw_lsf_tangent(d, radii, v, blend_k, s_1, s_2, s_3, max_deriv, &
                                    tangent_f0, tangent_f1_r, tangent_f2_rr, &
                                    tangent_f3_rrr, vrad)
      !> Displacements r - R_A
      real(wp), intent(in) :: d(:, :)
      !> Radii of the active atoms
      real(wp), intent(in) :: radii(:)
      !> Displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Blending sharpness k
      real(wp), intent(in) :: blend_k
      !> One-body weight
      real(wp), intent(in) :: s_1
      !> Two-body weight
      real(wp), intent(in) :: s_2
      !> Three-body weight
      real(wp), intent(in) :: s_3
      !> Highest spatial-derivative order
      integer, intent(in) :: max_deriv
      !> Contracted value derivative
      real(wp), intent(out) :: tangent_f0
      !> Contracted gradient derivative
      real(wp), intent(out) :: tangent_f1_r(3)
      !> Contracted Hessian derivative
      real(wp), intent(out) :: tangent_f2_rr(3, 3)
      !> Contracted third-derivative derivative
      real(wp), intent(out) :: tangent_f3_rrr(3, 3, 3)
      !> Radius half of the joint direction; present for a geometry-dependent
      !> radius model, absent leaves this a pure nuclear tangent
      real(wp), intent(in), optional :: vrad(:)
   end subroutine svdw_lsf_tangent

   !> Nuclear Hessian-vector products for every active atom.
   !>
   !> @param[in]  d             Displacements r - R_A, shape (3, n_active)
   !> @param[in]  radii         Radii of the active atoms, shape (n_active)
   !> @param[in]  v             Displacement directions, shape (3, n_active)
   !> @param[in]  blend_k       Blending sharpness k
   !> @param[in]  s_1           One-body weight
   !> @param[in]  s_2           Two-body weight
   !> @param[in]  s_3           Three-body weight
   !> @param[in]  max_deriv     Highest spatial-derivative order (0..2)
   !> @param[out] hvp_f1_rA     shape (3, n_active)
   !> @param[out] hvp_f2_r_rA   shape (3, 3, n_active)
   !> @param[out] hvp_f3_rr_rA  shape (3, 3, 3, n_active)
   pure module subroutine svdw_lsf_hvp(d, radii, v, blend_k, s_1, s_2, s_3, max_deriv, &
                                hvp_f1_rA, hvp_f2_r_rA, hvp_f3_rr_rA, vrad)
      !> Displacements r - R_A
      real(wp), intent(in) :: d(:, :)
      !> Radii of the active atoms
      real(wp), intent(in) :: radii(:)
      !> Displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Blending sharpness k
      real(wp), intent(in) :: blend_k
      !> One-body weight
      real(wp), intent(in) :: s_1
      !> Two-body weight
      real(wp), intent(in) :: s_2
      !> Three-body weight
      real(wp), intent(in) :: s_3
      !> Highest spatial-derivative order
      integer, intent(in) :: max_deriv
      !> Contracted nuclear Hessian
      real(wp), intent(out) :: hvp_f1_rA(:, :)
      !> Contracted mixed third derivative
      real(wp), intent(out) :: hvp_f2_r_rA(:, :, :)
      !> Contracted mixed fourth derivative
      real(wp), intent(out) :: hvp_f3_rr_rA(:, :, :, :)
      !> Radius half of the joint direction; present for a geometry-dependent
      !> radius model.  Supplying it turns these three results from the pure
      !> nuclear-nuclear Hessian block into the nuclear row of the *joint*
      !> position/radius Hessian; `svdw_lsf_radius_hvp` returns the radius row
      real(wp), intent(in), optional :: vrad(:)
   end subroutine svdw_lsf_hvp

   !> Radius row of the joint position/radius Hessian, over an active-atom list.
   !>
   !> Returns `sum_B (v_B . d/dR_B + vr_B d/dR_b) d/dR_a d^n S / dr^n`, n = 0..2,
   !> one value per active atom A and per spatial component.  Together with
   !> `svdw_lsf_hvp` called with the same `vrad` this is the complete Hessian-
   !> vector product of the level set with respect to `(X, R)`.
   !>
   !> `vrad` is *not* optional here: without a radius direction there is no
   !> radius row to compute, and a caller with a geometry-independent radius
   !> model should not be calling this at all.
   !>
   !> @param[in]  d              Displacements r - R_A, shape (3, n_active)
   !> @param[in]  radii          Radii of the active atoms
   !> @param[in]  v              Nuclear displacement directions
   !> @param[in]  vrad           Radius directions, vr_A = sum_B (dR_A/dX_B) . v_B
   !> @param[in]  blend_k        Blending sharpness k
   !> @param[in]  s_1            One-body weight
   !> @param[in]  s_2            Two-body weight
   !> @param[in]  s_3            Three-body weight
   !> @param[in]  max_deriv      Highest spatial-derivative order (0..2)
   !> @param[out] hvp_f1_rad     Contracted radius Hessian row
   !> @param[out] hvp_f2_r_rad   Contracted mixed third derivative
   !> @param[out] hvp_f3_rr_rad  Contracted mixed fourth derivative
   pure module subroutine svdw_lsf_radius_hvp(d, radii, v, vrad, blend_k, s_1, s_2, s_3, &
                                       max_deriv, hvp_f1_rad, hvp_f2_r_rad, &
                                       hvp_f3_rr_rad)
      !> Displacements r - R_A
      real(wp), intent(in) :: d(:, :)
      !> Radii of the active atoms
      real(wp), intent(in) :: radii(:)
      !> Displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Blending sharpness k
      real(wp), intent(in) :: blend_k
      !> One-body weight
      real(wp), intent(in) :: s_1
      !> Two-body weight
      real(wp), intent(in) :: s_2
      !> Three-body weight
      real(wp), intent(in) :: s_3
      !> Highest spatial-derivative order
      integer, intent(in) :: max_deriv
      !> Contracted radius Hessian row
      real(wp), intent(out) :: hvp_f1_rad(:)
      !> Contracted mixed third derivative
      real(wp), intent(out) :: hvp_f2_r_rad(:, :)
      !> Contracted mixed fourth derivative
      real(wp), intent(out) :: hvp_f3_rr_rad(:, :, :)
   end subroutine svdw_lsf_radius_hvp

   !> Normalized level set and its nuclear gradient over an active-atom list.
   !>
   !> @param[in]  d          Displacements r - R_A, shape (3, n_active)
   !> @param[in]  radii      Radii of the active atoms, shape (n_active)
   !> @param[in]  blend_k    Blending sharpness k
   !> @param[in]  s_1        One-body weight
   !> @param[in]  s_2        Two-body weight
   !> @param[in]  s_3        Three-body weight
   !> @param[out] norm0      S/||grad S||
   !> @param[out] norm1_rA   d/dR_A of S/||grad S||, shape (3, n_active)
   pure module subroutine svdw_lsf_normalized(d, radii, blend_k, s_1, s_2, s_3, &
                                       norm0, norm1_rA)
      !> Displacements r - R_A
      real(wp), intent(in) :: d(:, :)
      !> Radii of the active atoms
      real(wp), intent(in) :: radii(:)
      !> Blending sharpness k
      real(wp), intent(in) :: blend_k
      !> One-body weight
      real(wp), intent(in) :: s_1
      !> Two-body weight
      real(wp), intent(in) :: s_2
      !> Three-body weight
      real(wp), intent(in) :: s_3
      !> Normalized level-set value
      real(wp), intent(out) :: norm0
      !> Nuclear gradient of the normalized level-set value
      real(wp), intent(out) :: norm1_rA(:, :)
   end subroutine svdw_lsf_normalized

   end interface

end module moist_cavity_drop_lsf_svdw_kernel
