!> CFC (COSMO Fine Cavity) pseudo-density kernel
!>
!> Model
!> -----
!> Diedenhofen-Klamt 2018 pseudo-density PD and level set S, interior negative as in SvdW
!>
!>    PD(r) = sum_a exp{a1 (s_a - 1)}
!>          + sum_{a<b} c (1 - vec(s_a).vec(s_b))**4 exp{a2 (s_a + s_b - 2)}
!>    S(r)  = -log PD(r)
!>
!> - `s_a = ||r - r_a||/R_a`, `vec(s_a) = (r - r_a)/R_a`
!> - r enters only through the displacements `d_a = r - r_a`:
!>   `d/dr_i = d/d(d_a)_i + d/d(d_b)_i`, `d/dr_A,i = -d/d(d_A)_i`
!> - pair term does not factorize: O(n**2) work in the pair loop
!> - nuclear index contracted inside the kernel: O(n) storage, the pairwise nuclear
!>   Hessian is never formed
!>
!> Two layers
!> ----------
!> The log needs the complete sum, so a caller runs two layers in this order:
!>
!> 1. Accumulate PD, per atom and per pair
!>    - `*_term_eval`, `cfc_pair_spatial_eval`: spatial derivatives
!>    - `*_nuclear_eval`: one retained nuclear index, per atom
!>    - `*_tangent_eval`: one contracted nuclear index, aggregate
!>    - `*_hvp_eval`: one retained and one contracted nuclear index, per atom
!>    - nuclear signs already applied
!> 2. Lift to S
!>    - `cfc_spatial_eval`, `cfc_tangent_eval`: O(1) per point
!>    - `cfc_hvp_eval`: O(1) per atom
!>    - no chain rule, nuclear sign or split-to-spatial recombination left to the caller
!>
!> Buffers and orders
!> - accumulators `pd*`, `tg*` aggregate; `qn*`, `hvp*` with a trailing atom index
!> - accumulators `intent(inout)`: zeroed by the caller once per point, then only added
!>   to, orders above `max_deriv` included
!> - lift results `intent(out)`: all zeroed before the dispatch, so orders above
!>   `max_deriv` are exact zeros, never undefined
!> - dispatch on `max_deriv`: no high-order temporary for a low-order request; a request
!>   above the top order runs the top branch
!>
!> Fade
!> ----
!> Every layer-1 routine takes the band `fade_start`, `fade_width` (fs, fw below) in
!> `sigma_A = ||r - R_A||/R_a - 1`, as `cfc_strip_powersums_o*`, and accumulates
!>
!>    PD_f = sum_a w_a T_a + sum_{a<b} w_a w_b T_ab
!>    w    = 1 - S_4(t) = S_4(s),   t = clamp((sigma - fade_start)/fade_width, 0, 1),
!>    s    = 1 - t
!>
!> - w evaluated as `s**(n+1) P_n`, `P_n = sum_k C(n+k, k) t**k`, n = 4: positive
!>   coefficients, exact zero beyond the band
!> - band in e-folds of the screening factor `u_A = exp(a2 sigma_A)`:
!>   `fade_start = a/|a2|`, `fade_width = b/|a2|`, the same e-fold count for every radius
!> - `u_A u_B = exp(a2 (s_a + s_b - 2))` is the pair envelope; the atomic term
!>   `exp(a1 sigma_A) = u_A**(a1/a2)` decays faster
!> - sigma recomputed from the displacements in every routine
!> - every center on the plateau (`sigma <= fade_start`): plain body
!> - `fade_start = huge(0.0_wp)`, finite positive `fade_width`: unfaded kernel to the bit
!> - center beyond its band (`sigma >= fade_start + fade_width`): no contribution
!> - lift routines linear in the accumulators, no band
!>
!> Per-point recipe
!> ----------------
!> Zero every accumulator, then:
!>
!>   do A = 1, n_active
!>      call cfc_atomic_term_eval   (d_A, R_A, a1, fs, fw, n,   pd0..pd4_rrrr)
!>      call cfc_atomic_nuclear_eval(d_A, R_A, a1, fs, fw, n-1, qn*(.., A))
!>      call cfc_atomic_hvp_eval    (d_A, R_A, a1, fs, fw, v_A, n-1, tg*, hvp*(.., A))
!>      do B = A+1, n_active
!>         call cfc_pair_spatial_eval(d_A, d_B, .., fs, fw, n,   pd0..pd4_rrrr)
!>         call cfc_pair_nuclear_eval(d_A, d_B, .., fs, fw, n-1, qn*_a(.., A),
!>                                    qn*_b(.., B))
!>         call cfc_pair_hvp_eval    (d_A, d_B, .., fs, fw, v_A, v_B, n-1,
!>                                    tg*, hvp*_a(.., A), hvp*_b(.., B))
!>      end do
!>   end do
!>   call cfc_spatial_eval(pd*, n, f0, f1_r, f2_rr, f3_rrr, f4_rrrr)
!>   call cfc_tangent_eval(pd*, tg*, n-1, tangent_f0 .. tangent_f3_rrr)
!>   do A = 1, n_active
!>      call cfc_nuclear_eval(pd*, qn*(.., A), n-1, f1_rA .. f4_rrr_rA)
!>      call cfc_hvp_eval(pd*, tg*, qn*(.., A), hvp*(.., A), n-2,
!>                        hvp_f1_rA .. hvp_f3_rr_rA)
!>   end do
!>
!> - drop every call the requested order does not need
!> - values and spatial derivatives: the `*_term_eval` and `cfc_pair_spatial_eval` lines
!>   plus `cfc_spatial_eval`
!> - nuclear gradient: add `*_nuclear_eval` and `cfc_nuclear_eval`, no direction vector
!> - directional derivative without the HVP: `*_tangent_eval`; `*_hvp_eval` fills the
!>   same `tg*`, so never call both for one atom or pair
!>
!> Reverse mode, per-point adjoint jet (w0, w1, w2), in place of `cfc_nuclear_eval`:
!>
!>   call cfc_vjp_eval(pd*, qn*(.., A), w0, w1, w2, vjp_f1_rA)
!>
!> - same two families to the same order; 3 numbers per atom instead of 3 + 9 + 27
!> - jet contracted before the common subexpressions are taken, so the rank-3 tensor is
!>   never built
!>
!> Uncontracted nuclear Hessian: zero a per-pair `qq*`, fill it with `*_hessian_eval`
!> (one call for an off-diagonal pair, a sweep over the partners for a diagonal one),
!> then, once per ordered pair (A, B):
!>
!>   call cfc_hessian_eval(pd*, qn*(.., A), qn*(.., B), qq*, n-2,
!>                         f2_rArB, f3_r_rArB, f4_rr_rArB)
!>
!> - O(n**2) calls, the size of the result
!> - Hessian applied to a displacement field: `cfc_hvp_eval` instead, O(n)
!>
!> Cavity radii
!> ------------
!> - radii R_a are level-set parameters; a model that ties them to the geometry
!>   differentiates through them
!> - separate channel: a caller with fixed radii calls none of its routines and the
!>   others are unchanged to the bit
!> - R_a is a scalar parameter, not a coordinate: no sign flip, no trailing Cartesian
!>   index, no extra derivative order consumed
!> - `_rad` outputs have the rank of the spatial tensor they come from; `_rad`, not
!>   `_Ra`, since Fortran names are case-insensitive and `f3_rr_rA` would collide
!>
!>   do A = 1, n_active
!>      call cfc_atomic_radius_eval(d_A, R_A, a1, fs, fw, n-1, rd*(.., A))
!>      do B = A+1, n_active
!>         call cfc_pair_radius_eval(d_A, d_B, .., fs, fw, n-1, rd*_a(.., A),
!>                                   rd*_b(.., B))
!>      end do
!>   end do
!>   do A = 1, n_active
!>      call cfc_radius_eval(pd*, rd*(.., A), n-1, f1_rad .. f4_rrr_rad)
!>   end do
!>
!> Reverse mode, in place of `cfc_radius_eval`:
!>
!>   call cfc_radius_vjp_eval(pd*, rd*(.., A), w0, w1, w2, vjp_f1_rad)
!>
!> - same two families to the same order; 1 number per atom instead of 1 + 3 + 9
!> - jet contracted before the common subexpressions are taken, so the rank-2 tensor is
!>   never built
!>
!> Joint position/radius Hessian: the direction gains radius components vr, and the
!> contracted index of `tg*` and `hvp*` runs over radii too; one more call per atom and
!> per pair:
!>
!>   call cfc_atomic_radius_hvp_eval(d_A, R_A, a1, fs, fw, v_A, vr_A, n-1,
!>                                   tg*, hvp*(.., A), rh*(.., A))
!>
!> - `*_radius_hvp_eval` adds the radius half to `tg*`, `hvp*` and fills `rh*` in the
!>   same pass
!> - `cfc_hvp_eval` then gives the position row, `cfc_radius_hvp_eval` the radius row;
!>   together all four blocks (nuclear-nuclear, nuclear-radius, radius-nuclear,
!>   radius-radius)
!> - `*_radius_tangent_eval`: tangent-only counterpart; never call it and
!>   `*_radius_hvp_eval` for one atom or pair
!>
!> Uncontracted radius blocks, as `cfc_hessian_eval`: zero a per-pair buffer, fill it
!> with `*_radius2_eval` or `*_nucrad_eval` (one call for an off-diagonal pair, a sweep
!> over the partners for a diagonal one), then, once per ordered pair (A, B):
!>
!>   call cfc_radius2_eval(pd*, rd*(.., A), rd*(.., B), rr*, n-2,
!>                         f2_radrad, f3_r_radrad, f4_rr_radrad)
!>   call cfc_nucrad_eval(pd*, qn*(.., A), rd*(.., B), mx*, n-2,
!>                        f2_rA_rad, f3_r_rA_rad, f4_rr_rA_rad)
!>
!> - `cfc_radius2_eval` reads `rd*` in the A and the B slot; `cfc_nucrad_eval` reads
!>   `qn*` for the position of A and `rd*` for the radius of B
!> - radius-radius block symmetric, like `f2_rArB`; the nuclear-radius block is not, so
!>   the full block needs both orderings
!>
!> Strip kernel
!> ------------
!> - `cfc_strip_powersums_o{L}` (stage 1) and `cfc_strip_lift_o{L}` (stage 2),
!>   L = 0..4: spatial jet of S through order L for `np <= cfc_strip_width(L)` points
!>   sharing one atom list
!> - strip widths by order 0..4: 32, 32, 32, 16, 8
!> - stage 1: atom pass (atom outer, point inner: geometry cache, jets, atomic terms),
!>   then pair pass (pair outer, point inner) over the atoms that survive the strip's
!>   classification
!> - no branch inside a point loop: skipped, plateau and faded atoms decided per atom
!>   against the strip box; the gate is a separate select loop
!> - pair body: Leibniz product of the `f(q)` jet with both atoms' cached Cartesian jets
!>   of `g_A = w_A u_A`, one body for plateau and faded pairs
!> - faded atomic term `w_A exp(a1 sigma_A)`
!> - per-point path: a strip of one
!> - on-nucleus pairs counted in `n_zero`, not guarded; the caller falls back to the
!>   per-point path for such a strip
!>
!> Layouts
!> - points `(3, np)`
!> - stage-1 sums point-first: `pd0(np)`, `pd1(np, 3)`, `pd{n}(np, cfc_npack(n))`
!> - lift outputs point-last: `lsf0(np)`, `lsf1_r(3, np)`, `lsf2_rr(6, np)`,
!>   `lsf3_rrr(10, np)`, `lsf4_rrrr(15, np)`
!> - symmetric groups packed as sorted tuples a <= b <= c <= d, last index slowest,
!>   as `moist_math_packed`:
!>
!>    slot = a + b(b-1)/2 + (c-1)c(c+1)/6 + (d-1)d(d+1)(d+2)/24
!>
!> - geometry cache `geo(cfc_strip_width(L), cfc_strip_geo_slots(L), nat)`: slots 1-3
!>   `r - R_A`, 4 `1/||r - R_A||`, then the gate-masked Cartesian jet
!>   `d^J (w_A u_A)/dr^J`, |J| = 0..L, packed, by order
!> - one caller buffer of `cfc_strip_geo_len` doubles per atom serves every order by
!>   sequence association
!>
!> Inline exp
!> - `exp(x) = 2**k (1 + r + r**2 q(r))`, `k = anint(x/ln 2)`, two-constant Cody-Waite
!>   reduction to |r| <= ln2/2, q of degree 10 (total degree 12)
!> - `2**k` from its exponent bits; same constants as the SvdW kernel
!> - argument clamped below at -708; above, `max(-a1, -a2) <= 708`
!>   checked per strip (the argument is `a sigma_A <= -a` for `a <= 0`)
!> - <= 1 ulp, measured by the SvdW generator's `--check` sweep with and without FMA
!>
!> Strip nuclear vjp
!> -----------------
!> - `cfc_strip_lift_vjp` and `cfc_strip_vjp`, after stage 1 at order 2 on the same
!>   strip (`np <= cfc_strip_width(2) = 32`)
!> - batched `res(:, a, p) = d/dR_a [w0 S + w1 . dS/dr + W : d2S/dr2](r_p)` over the
!>   shared atom list
!> - adjoint lift: seeds `w0(np)`, `w1(3, np)`, `w2(6, np)` (symmetric part, packed,
!>   not doubled) and the order-2 sums to `lam[I] = d(w . jet)/d pd[I]`, in the layout
!>   of the sums
!> - contraction: atom outer (atomic term) and pair outer (`A < B`), point inner, over
!>   the stage-1 cache `geo` (read, slots 1-14)
!> - exp and fade recomputed per survivor; the rank-3 Cartesian jet of `g = w u` goes to
!>   the workspace `jet3(nps, 10, nat)`, since `d/dR_A` of the order-2 jet is rank 3
!> - `lam` contracted into the gradients symbolically before emission: six numbers per
!>   (pair, point), for centers A and B; the `qn*` tensors of the per-point path are
!>   never formed
!> - classification as in stage 1 (skip, plateau, band, pair reach test); beyond reach
!>   exact zeros; on-nucleus lanes counted in `n_zero`
!> - unreached columns (`pd0 <= 0` or not finite) have no valid `lam`; the caller
!>   discards their `res` columns
!>
module moist_cavity_drop_lsf_cfc_kernel
   use mctc_env_accuracy, only: wp
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none (type, external)
   private

   !> Layer 1: pseudo-density accumulators (per atom / per pair)
   public :: cfc_atomic_term_eval
   public :: cfc_pair_term_eval
   public :: cfc_pair_spatial_eval
   public :: cfc_atomic_nuclear_eval
   public :: cfc_pair_nuclear_eval
   public :: cfc_atomic_hessian_eval
   public :: cfc_pair_hessian_eval
   public :: cfc_atomic_tangent_eval
   public :: cfc_pair_tangent_eval
   public :: cfc_atomic_hvp_eval
   public :: cfc_pair_hvp_eval

   !> Layer 1: cavity-radius accumulators (per atom / per pair)
   public :: cfc_atomic_radius_eval
   public :: cfc_pair_radius_eval
   public :: cfc_atomic_radius2_eval
   public :: cfc_pair_radius2_eval
   public :: cfc_atomic_nucrad_eval
   public :: cfc_pair_nucrad_eval
   public :: cfc_atomic_radius_tangent_eval
   public :: cfc_pair_radius_tangent_eval
   public :: cfc_atomic_radius_hvp_eval
   public :: cfc_pair_radius_hvp_eval

   !> Layer 2: level-set lift (aggregate / per atom)
   public :: cfc_spatial_eval
   public :: cfc_nuclear_eval
   public :: cfc_hessian_eval
   public :: cfc_tangent_eval
   public :: cfc_hvp_eval
   public :: cfc_vjp_eval
   public :: cfc_radius_eval
   public :: cfc_radius2_eval
   public :: cfc_nucrad_eval
   public :: cfc_radius_hvp_eval
   public :: cfc_radius_vjp_eval

   !> Points per strip of the strip kernels, by derivative order 0..4
   integer, parameter, public :: cfc_strip_width(0:4) = [32, 32, 32, 16, 8]
   !> Packed components of a symmetric rank-n spatial tensor, n = 0..4
   integer, parameter, public :: cfc_npack(0:4) = [1, 3, 6, 10, 15]
   !> Slots per point and atom of the strip geometry cache, by order 0..4
   integer, parameter, public :: cfc_strip_geo_slots(0:4) = [5, 8, 14, 24, 39]
   !> Doubles per atom of a geometry-cache buffer that serves every order
   integer, parameter, public :: cfc_strip_geo_len = 448

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

   !> Strip kernel: pseudo-density jet (stage 1) and lift (stage 2)
   public :: cfc_strip_powersums_o0
   public :: cfc_strip_lift_o0
   public :: cfc_strip_powersums_o1
   public :: cfc_strip_lift_o1
   public :: cfc_strip_powersums_o2
   public :: cfc_strip_lift_o2
   public :: cfc_strip_powersums_o3
   public :: cfc_strip_lift_o3
   public :: cfc_strip_powersums_o4
   public :: cfc_strip_lift_o4

   !> Strip nuclear vjp: adjoint lift and per-atom / per-pair contraction
   public :: cfc_strip_lift_vjp
   public :: cfc_strip_vjp

   !* ================================================================================= *!
   !*                                     Interfaces                                    *!
   !* ================================================================================= *!

   interface

   !> Accumulate one atom's contribution to PD and its spatial derivatives.
   !>
   !> Computes ``exp(a1 (s_a - 1))`` with ``s_a = ||d_a|| / R_a`` and adds the
   !> value plus ``d^k PD / dr^k`` (k = 1..max_deriv) into the accumulators.
   !> The orchestrator zeroes the accumulators once per evaluation point and
   !> then calls this routine for every active atom.
   !>
   !> The expression depends on r only through ``d_a = r - r_a``, so spatial
   !> and displacement derivatives coincide, and the nuclear derivatives are
   !> ``d^k/dr_a^k = (-1)^k d^k/d(d_a)^k``.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Maximum derivative order to accumulate (0..4); a larger value is
   !>                          treated as 4
   !> @param[inout] pd0        Pseudo-density value accumulator
   !> @param[inout] pd1_r      d^1 K / d(d_a)^1 accumulator (fully symmetric)
   !> @param[inout] pd2_rr     d^2 K / d(d_a)^2 accumulator (fully symmetric)
   !> @param[inout] pd3_rrr    d^3 K / d(d_a)^3 accumulator (fully symmetric)
   !> @param[inout] pd4_rrrr   d^4 K / d(d_a)^4 accumulator (fully symmetric)
   pure module subroutine cfc_atomic_term_eval(d_a, R_a, a1, fade_start, fade_width, max_deriv, pd0, &
                                        pd1_r, pd2_rr, pd3_rrr, pd4_rrrr)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Maximum derivative order to accumulate (0..4); a larger value is treated as 4
      integer, intent(in) :: max_deriv
      !> Pseudo-density value accumulator
      real(wp), intent(inout) :: pd0
      !> d^1 K / d(d_a)^1 accumulator (fully symmetric)
      real(wp), intent(inout) :: pd1_r(3)
      !> d^2 K / d(d_a)^2 accumulator (fully symmetric)
      real(wp), intent(inout) :: pd2_rr(3, 3)
      !> d^3 K / d(d_a)^3 accumulator (fully symmetric)
      real(wp), intent(inout) :: pd3_rrr(3, 3, 3)
      !> d^4 K / d(d_a)^4 accumulator (fully symmetric)
      real(wp), intent(inout) :: pd4_rrrr(3, 3, 3, 3)
   end subroutine cfc_atomic_term_eval

   !> Faded term of `cfc_atomic_term_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Maximum derivative order to accumulate (0..4); a larger value is
   !>                          treated as 4
   !> @param[inout] pd0        Pseudo-density value accumulator
   !> @param[inout] pd1_r      d^1 K / d(d_a)^1 accumulator (fully symmetric)
   !> @param[inout] pd2_rr     d^2 K / d(d_a)^2 accumulator (fully symmetric)
   !> @param[inout] pd3_rrr    d^3 K / d(d_a)^3 accumulator (fully symmetric)
   !> @param[inout] pd4_rrrr   d^4 K / d(d_a)^4 accumulator (fully symmetric)
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   pure module subroutine cfc_atomic_term_eval_faded(d_a, R_a, a1, fade_start, fade_width, max_deriv, &
                                              pd0, pd1_r, pd2_rr, pd3_rrr, pd4_rrrr, sg_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Maximum derivative order to accumulate (0..4); a larger value is treated as 4
      integer, intent(in) :: max_deriv
      !> Pseudo-density value accumulator
      real(wp), intent(inout) :: pd0
      !> d^1 K / d(d_a)^1 accumulator (fully symmetric)
      real(wp), intent(inout) :: pd1_r(3)
      !> d^2 K / d(d_a)^2 accumulator (fully symmetric)
      real(wp), intent(inout) :: pd2_rr(3, 3)
      !> d^3 K / d(d_a)^3 accumulator (fully symmetric)
      real(wp), intent(inout) :: pd3_rrr(3, 3, 3)
      !> d^4 K / d(d_a)^4 accumulator (fully symmetric)
      real(wp), intent(inout) :: pd4_rrrr(3, 3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
   end subroutine cfc_atomic_term_eval_faded

   !> Accumulate one pair's contribution to PD, split by displacement vector.
   !>
   !> Computes ``c (1 - vec(s_a).vec(s_b))**m exp(a2 (s_a + s_b - 2))`` with
   !> ``m = 4`` baked into the symbolic differentiation, and accumulates the
   !> value together with every partial derivative taken separately with
   !> respect to ``d_a = r - r_a`` and ``d_b = r - r_b``, up to total order
   !> 4.
   !>
   !> Spatial and nuclear derivatives follow from the splits.  With
   !> ``d/dr = d/d(d_a) + d/d(d_b)`` and ``d/dr_A = -d/d(d_A)``:
   !>
   !>    dK/dr_i             = pd1_a_i + pd1_b_i
   !>    d2K/dr_i dr_j       = pd2_aa_ij + pd2_ab_ij + pd2_ab_ji + pd2_bb_ij
   !>    d3K/dr_i dr_j dr_k  = pd3_aaa + sym(pd3_aab) + sym(pd3_abb) + pd3_bbb
   !>    d4K/dr^4            = pd4_aaaa + sym(pd4_aaab) + sym(pd4_aabb)
   !>                          + sym(pd4_abbb) + pd4_bbbb
   !>
   !> where ``sym`` sums the distinct index assignments of the spatial indices
   !> to the a/b slots (4 for aaab and abbb, 6 for aabb, 3 for aab and abb).
   !> Nuclear derivatives pick up one minus sign per nuclear index.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Maximum derivative order to accumulate (0..4); a larger value is
   !>                          treated as 4
   !> @param[inout] pd0        Pseudo-density value accumulator
   !> @param[inout] pd1_a      d^1 K / d(d_a)^1 d(d_b)^0 accumulator (no index symmetry)
   !> @param[inout] pd1_b      d^1 K / d(d_a)^0 d(d_b)^1 accumulator (no index symmetry)
   !> @param[inout] pd2_aa     d^2 K / d(d_a)^2 d(d_b)^0 accumulator (symmetric in the first 2
   !>                          indices)
   !> @param[inout] pd2_ab     d^2 K / d(d_a)^1 d(d_b)^1 accumulator (no index symmetry)
   !> @param[inout] pd2_bb     d^2 K / d(d_a)^0 d(d_b)^2 accumulator (symmetric in the last 2
   !>                          indices)
   !> @param[inout] pd3_aaa    d^3 K / d(d_a)^3 d(d_b)^0 accumulator (symmetric in the first 3
   !>                          indices)
   !> @param[inout] pd3_aab    d^3 K / d(d_a)^2 d(d_b)^1 accumulator (symmetric in the first 2
   !>                          indices)
   !> @param[inout] pd3_abb    d^3 K / d(d_a)^1 d(d_b)^2 accumulator (symmetric in the last 2
   !>                          indices)
   !> @param[inout] pd3_bbb    d^3 K / d(d_a)^0 d(d_b)^3 accumulator (symmetric in the last 3
   !>                          indices)
   !> @param[inout] pd4_aaaa   d^4 K / d(d_a)^4 d(d_b)^0 accumulator (symmetric in the first 4
   !>                          indices)
   !> @param[inout] pd4_aaab   d^4 K / d(d_a)^3 d(d_b)^1 accumulator (symmetric in the first 3
   !>                          indices)
   !> @param[inout] pd4_aabb   d^4 K / d(d_a)^2 d(d_b)^2 accumulator (symmetric within each index
   !>                          group)
   !> @param[inout] pd4_abbb   d^4 K / d(d_a)^1 d(d_b)^3 accumulator (symmetric in the last 3
   !>                          indices)
   !> @param[inout] pd4_bbbb   d^4 K / d(d_a)^0 d(d_b)^4 accumulator (symmetric in the last 4
   !>                          indices)
   pure module subroutine cfc_pair_term_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, fade_width, &
                                      max_deriv, pd0, pd1_a, pd1_b, pd2_aa, pd2_ab, pd2_bb, &
                                      pd3_aaa, pd3_aab, pd3_abb, pd3_bbb, pd4_aaaa, pd4_aaab, &
                                      pd4_aabb, pd4_abbb, pd4_bbbb)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Maximum derivative order to accumulate (0..4); a larger value is treated as 4
      integer, intent(in) :: max_deriv
      !> Pseudo-density value accumulator
      real(wp), intent(inout) :: pd0
      !> d^1 K / d(d_a)^1 d(d_b)^0 accumulator (no index symmetry)
      real(wp), intent(inout) :: pd1_a(3)
      !> d^1 K / d(d_a)^0 d(d_b)^1 accumulator (no index symmetry)
      real(wp), intent(inout) :: pd1_b(3)
      !> d^2 K / d(d_a)^2 d(d_b)^0 accumulator (symmetric in the first 2 indices)
      real(wp), intent(inout) :: pd2_aa(3, 3)
      !> d^2 K / d(d_a)^1 d(d_b)^1 accumulator (no index symmetry)
      real(wp), intent(inout) :: pd2_ab(3, 3)
      !> d^2 K / d(d_a)^0 d(d_b)^2 accumulator (symmetric in the last 2 indices)
      real(wp), intent(inout) :: pd2_bb(3, 3)
      !> d^3 K / d(d_a)^3 d(d_b)^0 accumulator (symmetric in the first 3 indices)
      real(wp), intent(inout) :: pd3_aaa(3, 3, 3)
      !> d^3 K / d(d_a)^2 d(d_b)^1 accumulator (symmetric in the first 2 indices)
      real(wp), intent(inout) :: pd3_aab(3, 3, 3)
      !> d^3 K / d(d_a)^1 d(d_b)^2 accumulator (symmetric in the last 2 indices)
      real(wp), intent(inout) :: pd3_abb(3, 3, 3)
      !> d^3 K / d(d_a)^0 d(d_b)^3 accumulator (symmetric in the last 3 indices)
      real(wp), intent(inout) :: pd3_bbb(3, 3, 3)
      !> d^4 K / d(d_a)^4 d(d_b)^0 accumulator (symmetric in the first 4 indices)
      real(wp), intent(inout) :: pd4_aaaa(3, 3, 3, 3)
      !> d^4 K / d(d_a)^3 d(d_b)^1 accumulator (symmetric in the first 3 indices)
      real(wp), intent(inout) :: pd4_aaab(3, 3, 3, 3)
      !> d^4 K / d(d_a)^2 d(d_b)^2 accumulator (symmetric within each index group)
      real(wp), intent(inout) :: pd4_aabb(3, 3, 3, 3)
      !> d^4 K / d(d_a)^1 d(d_b)^3 accumulator (symmetric in the last 3 indices)
      real(wp), intent(inout) :: pd4_abbb(3, 3, 3, 3)
      !> d^4 K / d(d_a)^0 d(d_b)^4 accumulator (symmetric in the last 4 indices)
      real(wp), intent(inout) :: pd4_bbbb(3, 3, 3, 3)
   end subroutine cfc_pair_term_eval

   !> Faded term of `cfc_pair_term_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Maximum derivative order to accumulate (0..4); a larger value is
   !>                          treated as 4
   !> @param[inout] pd0        Pseudo-density value accumulator
   !> @param[inout] pd1_a      d^1 K / d(d_a)^1 d(d_b)^0 accumulator (no index symmetry)
   !> @param[inout] pd1_b      d^1 K / d(d_a)^0 d(d_b)^1 accumulator (no index symmetry)
   !> @param[inout] pd2_aa     d^2 K / d(d_a)^2 d(d_b)^0 accumulator (symmetric in the first 2
   !>                          indices)
   !> @param[inout] pd2_ab     d^2 K / d(d_a)^1 d(d_b)^1 accumulator (no index symmetry)
   !> @param[inout] pd2_bb     d^2 K / d(d_a)^0 d(d_b)^2 accumulator (symmetric in the last 2
   !>                          indices)
   !> @param[inout] pd3_aaa    d^3 K / d(d_a)^3 d(d_b)^0 accumulator (symmetric in the first 3
   !>                          indices)
   !> @param[inout] pd3_aab    d^3 K / d(d_a)^2 d(d_b)^1 accumulator (symmetric in the first 2
   !>                          indices)
   !> @param[inout] pd3_abb    d^3 K / d(d_a)^1 d(d_b)^2 accumulator (symmetric in the last 2
   !>                          indices)
   !> @param[inout] pd3_bbb    d^3 K / d(d_a)^0 d(d_b)^3 accumulator (symmetric in the last 3
   !>                          indices)
   !> @param[inout] pd4_aaaa   d^4 K / d(d_a)^4 d(d_b)^0 accumulator (symmetric in the first 4
   !>                          indices)
   !> @param[inout] pd4_aaab   d^4 K / d(d_a)^3 d(d_b)^1 accumulator (symmetric in the first 3
   !>                          indices)
   !> @param[inout] pd4_aabb   d^4 K / d(d_a)^2 d(d_b)^2 accumulator (symmetric within each index
   !>                          group)
   !> @param[inout] pd4_abbb   d^4 K / d(d_a)^1 d(d_b)^3 accumulator (symmetric in the last 3
   !>                          indices)
   !> @param[inout] pd4_bbbb   d^4 K / d(d_a)^0 d(d_b)^4 accumulator (symmetric in the last 4
   !>                          indices)
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_term_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                            fade_width, max_deriv, pd0, pd1_a, pd1_b, pd2_aa, &
                                            pd2_ab, pd2_bb, pd3_aaa, pd3_aab, pd3_abb, pd3_bbb, &
                                            pd4_aaaa, pd4_aaab, pd4_aabb, pd4_abbb, pd4_bbbb, &
                                            sg_a, sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Maximum derivative order to accumulate (0..4); a larger value is treated as 4
      integer, intent(in) :: max_deriv
      !> Pseudo-density value accumulator
      real(wp), intent(inout) :: pd0
      !> d^1 K / d(d_a)^1 d(d_b)^0 accumulator (no index symmetry)
      real(wp), intent(inout) :: pd1_a(3)
      !> d^1 K / d(d_a)^0 d(d_b)^1 accumulator (no index symmetry)
      real(wp), intent(inout) :: pd1_b(3)
      !> d^2 K / d(d_a)^2 d(d_b)^0 accumulator (symmetric in the first 2 indices)
      real(wp), intent(inout) :: pd2_aa(3, 3)
      !> d^2 K / d(d_a)^1 d(d_b)^1 accumulator (no index symmetry)
      real(wp), intent(inout) :: pd2_ab(3, 3)
      !> d^2 K / d(d_a)^0 d(d_b)^2 accumulator (symmetric in the last 2 indices)
      real(wp), intent(inout) :: pd2_bb(3, 3)
      !> d^3 K / d(d_a)^3 d(d_b)^0 accumulator (symmetric in the first 3 indices)
      real(wp), intent(inout) :: pd3_aaa(3, 3, 3)
      !> d^3 K / d(d_a)^2 d(d_b)^1 accumulator (symmetric in the first 2 indices)
      real(wp), intent(inout) :: pd3_aab(3, 3, 3)
      !> d^3 K / d(d_a)^1 d(d_b)^2 accumulator (symmetric in the last 2 indices)
      real(wp), intent(inout) :: pd3_abb(3, 3, 3)
      !> d^3 K / d(d_a)^0 d(d_b)^3 accumulator (symmetric in the last 3 indices)
      real(wp), intent(inout) :: pd3_bbb(3, 3, 3)
      !> d^4 K / d(d_a)^4 d(d_b)^0 accumulator (symmetric in the first 4 indices)
      real(wp), intent(inout) :: pd4_aaaa(3, 3, 3, 3)
      !> d^4 K / d(d_a)^3 d(d_b)^1 accumulator (symmetric in the first 3 indices)
      real(wp), intent(inout) :: pd4_aaab(3, 3, 3, 3)
      !> d^4 K / d(d_a)^2 d(d_b)^2 accumulator (symmetric within each index group)
      real(wp), intent(inout) :: pd4_aabb(3, 3, 3, 3)
      !> d^4 K / d(d_a)^1 d(d_b)^3 accumulator (symmetric in the last 3 indices)
      real(wp), intent(inout) :: pd4_abbb(3, 3, 3, 3)
      !> d^4 K / d(d_a)^0 d(d_b)^4 accumulator (symmetric in the last 4 indices)
      real(wp), intent(inout) :: pd4_bbbb(3, 3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_term_eval_faded

   !> Accumulate the nuclear directional derivative of the spatial derivatives.
   !>
   !> Given one displacement direction per atom, this returns
   !>
   !>    tg{{k}}(i1..ik) = sum_A sum_alpha v_A,alpha
   !>                      d/dr_A,alpha  d^k K / dr_i1..dr_ik
   !>
   !> i.e. the single nuclear index is contracted away and the result is a
   !> purely spatial tensor.  Summed over atoms (atomic term) and pairs (pair
   !> term) this is the directional derivative of the whole pseudo-density
   !> tensor, which is what the caller wants; the full per-atom tensor is never
   !> formed.  Signs of the nuclear derivatives are already applied.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A v_A . d/dr_A of the order-0 spatial derivative
   !> @param[inout] tg1_r      sum_A v_A . d/dr_A of the order-1 spatial derivative
   !> @param[inout] tg2_rr     sum_A v_A . d/dr_A of the order-2 spatial derivative
   !> @param[inout] tg3_rrr    sum_A v_A . d/dr_A of the order-3 spatial derivative
   pure module subroutine cfc_atomic_tangent_eval(d_a, R_a, a1, fade_start, fade_width, v_a, max_deriv, &
                                           tg0, tg1_r, tg2_rr, tg3_rrr)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A v_A . d/dr_A of the order-0 spatial derivative
      real(wp), intent(inout) :: tg0
      !> sum_A v_A . d/dr_A of the order-1 spatial derivative
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A v_A . d/dr_A of the order-2 spatial derivative
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A v_A . d/dr_A of the order-3 spatial derivative
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
   end subroutine cfc_atomic_tangent_eval

   !> Faded term of `cfc_atomic_tangent_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A v_A . d/dr_A of the order-0 spatial derivative
   !> @param[inout] tg1_r      sum_A v_A . d/dr_A of the order-1 spatial derivative
   !> @param[inout] tg2_rr     sum_A v_A . d/dr_A of the order-2 spatial derivative
   !> @param[inout] tg3_rrr    sum_A v_A . d/dr_A of the order-3 spatial derivative
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   pure module subroutine cfc_atomic_tangent_eval_faded(d_a, R_a, a1, fade_start, fade_width, v_a, &
                                                 max_deriv, tg0, tg1_r, tg2_rr, tg3_rrr, sg_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A v_A . d/dr_A of the order-0 spatial derivative
      real(wp), intent(inout) :: tg0
      !> sum_A v_A . d/dr_A of the order-1 spatial derivative
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A v_A . d/dr_A of the order-2 spatial derivative
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A v_A . d/dr_A of the order-3 spatial derivative
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
   end subroutine cfc_atomic_tangent_eval_faded

   !> Accumulate the nuclear directional derivative of the spatial derivatives.
   !>
   !> Given one displacement direction per atom, this returns
   !>
   !>    tg{{k}}(i1..ik) = sum_A sum_alpha v_A,alpha
   !>                      d/dr_A,alpha  d^k K / dr_i1..dr_ik
   !>
   !> i.e. the single nuclear index is contracted away and the result is a
   !> purely spatial tensor.  Summed over atoms (atomic term) and pairs (pair
   !> term) this is the directional derivative of the whole pseudo-density
   !> tensor, which is what the caller wants; the full per-atom tensor is never
   !> formed.  Signs of the nuclear derivatives are already applied.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    v_b        Nuclear displacement direction of atom b
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A v_A . d/dr_A of the order-0 spatial derivative
   !> @param[inout] tg1_r      sum_A v_A . d/dr_A of the order-1 spatial derivative
   !> @param[inout] tg2_rr     sum_A v_A . d/dr_A of the order-2 spatial derivative
   !> @param[inout] tg3_rrr    sum_A v_A . d/dr_A of the order-3 spatial derivative
   pure module subroutine cfc_pair_tangent_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, fade_width, &
                                         v_a, v_b, max_deriv, tg0, tg1_r, tg2_rr, tg3_rrr)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Nuclear displacement direction of atom b
      real(wp), intent(in) :: v_b(3)
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A v_A . d/dr_A of the order-0 spatial derivative
      real(wp), intent(inout) :: tg0
      !> sum_A v_A . d/dr_A of the order-1 spatial derivative
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A v_A . d/dr_A of the order-2 spatial derivative
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A v_A . d/dr_A of the order-3 spatial derivative
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
   end subroutine cfc_pair_tangent_eval

   !> Faded term of `cfc_pair_tangent_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    v_b        Nuclear displacement direction of atom b
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A v_A . d/dr_A of the order-0 spatial derivative
   !> @param[inout] tg1_r      sum_A v_A . d/dr_A of the order-1 spatial derivative
   !> @param[inout] tg2_rr     sum_A v_A . d/dr_A of the order-2 spatial derivative
   !> @param[inout] tg3_rrr    sum_A v_A . d/dr_A of the order-3 spatial derivative
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_tangent_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                               fade_width, v_a, v_b, max_deriv, tg0, tg1_r, &
                                               tg2_rr, tg3_rrr, sg_a, sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Nuclear displacement direction of atom b
      real(wp), intent(in) :: v_b(3)
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A v_A . d/dr_A of the order-0 spatial derivative
      real(wp), intent(inout) :: tg0
      !> sum_A v_A . d/dr_A of the order-1 spatial derivative
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A v_A . d/dr_A of the order-2 spatial derivative
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A v_A . d/dr_A of the order-3 spatial derivative
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_tangent_eval_faded

   !> Accumulate a nuclear Hessian-vector product, plus the tangent ladder.
   !>
   !> Given one displacement direction per atom, this returns, for each center
   !> A of this atom,
   !>
   !>    hvp{k}_A(i1..ik, beta) = sum_B sum_alpha v_B,alpha
   !>                               d^2/(dr_A,beta dr_B,alpha)
   !>                               d^k K / dr_i1..dr_ik
   !>
   !> i.e. exactly one of the two nuclear indices is contracted, the other is
   !> retained as the trailing array index.  The caller scatters the ``_a`` and
   !> ``_b`` blocks into the columns of its per-atom result, turning what would
   !> be an O(n**2) tensor into O(n) storage.  Signs are already applied (the
   !> two nuclear sign flips cancel).
   !>
   !> The `tg*` accumulators of `*_tangent_eval` are filled as well.  They cost
   !> nothing: the two ladders share almost all of their derivative structure,
   !> and CSE'd together the top branch needs 1752 temporaries against 3298 for
   !> the two computed apart.  CALL EITHER THIS ROUTINE OR `*_tangent_eval` FOR
   !> A GIVEN PAIR, NEVER BOTH -- the tangent contribution would be counted
   !> twice.  The HVP ladder stops at spatial order 2 (its total order is one
   !> higher than the tangent's); the tangent ladder runs to 3.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A v_A . d/dr_A of the order-0 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg1_r      sum_A v_A . d/dr_A of the order-1 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg2_rr     sum_A v_A . d/dr_A of the order-2 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg3_rrr    sum_A v_A . d/dr_A of the order-3 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] hvp0_a     sum_B v_B . d2/(dr_a dr_B) of the order-0 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[inout] hvp1_r_a   sum_B v_B . d2/(dr_a dr_B) of the order-1 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[inout] hvp2_rr_a  sum_B v_B . d2/(dr_a dr_B) of the order-2 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   pure module subroutine cfc_atomic_hvp_eval(d_a, R_a, a1, fade_start, fade_width, v_a, max_deriv, tg0, &
                                       tg1_r, tg2_rr, tg3_rrr, hvp0_a, hvp1_r_a, hvp2_rr_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A v_A . d/dr_A of the order-0 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg0
      !> sum_A v_A . d/dr_A of the order-1 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A v_A . d/dr_A of the order-2 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A v_A . d/dr_A of the order-3 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-0 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp0_a(3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-1 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp1_r_a(3, 3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-2 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp2_rr_a(3, 3, 3)
   end subroutine cfc_atomic_hvp_eval

   !> Faded term of `cfc_atomic_hvp_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A v_A . d/dr_A of the order-0 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg1_r      sum_A v_A . d/dr_A of the order-1 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg2_rr     sum_A v_A . d/dr_A of the order-2 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg3_rrr    sum_A v_A . d/dr_A of the order-3 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] hvp0_a     sum_B v_B . d2/(dr_a dr_B) of the order-0 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[inout] hvp1_r_a   sum_B v_B . d2/(dr_a dr_B) of the order-1 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[inout] hvp2_rr_a  sum_B v_B . d2/(dr_a dr_B) of the order-2 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   pure module subroutine cfc_atomic_hvp_eval_faded(d_a, R_a, a1, fade_start, fade_width, v_a, &
                                             max_deriv, tg0, tg1_r, tg2_rr, tg3_rrr, hvp0_a, &
                                             hvp1_r_a, hvp2_rr_a, sg_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A v_A . d/dr_A of the order-0 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg0
      !> sum_A v_A . d/dr_A of the order-1 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A v_A . d/dr_A of the order-2 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A v_A . d/dr_A of the order-3 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-0 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp0_a(3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-1 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp1_r_a(3, 3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-2 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp2_rr_a(3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
   end subroutine cfc_atomic_hvp_eval_faded

   !> Accumulate a nuclear Hessian-vector product, plus the tangent ladder.
   !>
   !> Given one displacement direction per atom, this returns, for each center
   !> A of this pair,
   !>
   !>    hvp{k}_A(i1..ik, beta) = sum_B sum_alpha v_B,alpha
   !>                               d^2/(dr_A,beta dr_B,alpha)
   !>                               d^k K / dr_i1..dr_ik
   !>
   !> i.e. exactly one of the two nuclear indices is contracted, the other is
   !> retained as the trailing array index.  The caller scatters the ``_a`` and
   !> ``_b`` blocks into the columns of its per-atom result, turning what would
   !> be an O(n**2) tensor into O(n) storage.  Signs are already applied (the
   !> two nuclear sign flips cancel).
   !>
   !> The `tg*` accumulators of `*_tangent_eval` are filled as well.  They cost
   !> nothing: the two ladders share almost all of their derivative structure,
   !> and CSE'd together the top branch needs 1752 temporaries against 3298 for
   !> the two computed apart.  CALL EITHER THIS ROUTINE OR `*_tangent_eval` FOR
   !> A GIVEN PAIR, NEVER BOTH -- the tangent contribution would be counted
   !> twice.  The HVP ladder stops at spatial order 2 (its total order is one
   !> higher than the tangent's); the tangent ladder runs to 3.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    v_b        Nuclear displacement direction of atom b
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A v_A . d/dr_A of the order-0 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg1_r      sum_A v_A . d/dr_A of the order-1 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg2_rr     sum_A v_A . d/dr_A of the order-2 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg3_rrr    sum_A v_A . d/dr_A of the order-3 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] hvp0_a     sum_B v_B . d2/(dr_a dr_B) of the order-0 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[inout] hvp0_b     sum_B v_B . d2/(dr_b dr_B) of the order-0 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom b
   !> @param[inout] hvp1_r_a   sum_B v_B . d2/(dr_a dr_B) of the order-1 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[inout] hvp1_r_b   sum_B v_B . d2/(dr_b dr_B) of the order-1 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom b
   !> @param[inout] hvp2_rr_a  sum_B v_B . d2/(dr_a dr_B) of the order-2 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[inout] hvp2_rr_b  sum_B v_B . d2/(dr_b dr_B) of the order-2 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom b
   pure module subroutine cfc_pair_hvp_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, fade_width, v_a, &
                                     v_b, max_deriv, tg0, tg1_r, tg2_rr, tg3_rrr, hvp0_a, hvp0_b, &
                                     hvp1_r_a, hvp1_r_b, hvp2_rr_a, hvp2_rr_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Nuclear displacement direction of atom b
      real(wp), intent(in) :: v_b(3)
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A v_A . d/dr_A of the order-0 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg0
      !> sum_A v_A . d/dr_A of the order-1 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A v_A . d/dr_A of the order-2 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A v_A . d/dr_A of the order-3 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-0 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp0_a(3)
      !> sum_B v_B . d2/(dr_b dr_B) of the order-0 spatial derivative; trailing index is the
      !> retained nuclear component of atom b
      real(wp), intent(inout) :: hvp0_b(3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-1 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp1_r_a(3, 3)
      !> sum_B v_B . d2/(dr_b dr_B) of the order-1 spatial derivative; trailing index is the
      !> retained nuclear component of atom b
      real(wp), intent(inout) :: hvp1_r_b(3, 3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-2 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp2_rr_a(3, 3, 3)
      !> sum_B v_B . d2/(dr_b dr_B) of the order-2 spatial derivative; trailing index is the
      !> retained nuclear component of atom b
      real(wp), intent(inout) :: hvp2_rr_b(3, 3, 3)
   end subroutine cfc_pair_hvp_eval

   !> Faded term of `cfc_pair_hvp_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    v_b        Nuclear displacement direction of atom b
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A v_A . d/dr_A of the order-0 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg1_r      sum_A v_A . d/dr_A of the order-1 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg2_rr     sum_A v_A . d/dr_A of the order-2 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] tg3_rrr    sum_A v_A . d/dr_A of the order-3 spatial derivative (free by-product
   !>                          of the HVP; see the note above)
   !> @param[inout] hvp0_a     sum_B v_B . d2/(dr_a dr_B) of the order-0 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[inout] hvp0_b     sum_B v_B . d2/(dr_b dr_B) of the order-0 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom b
   !> @param[inout] hvp1_r_a   sum_B v_B . d2/(dr_a dr_B) of the order-1 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[inout] hvp1_r_b   sum_B v_B . d2/(dr_b dr_B) of the order-1 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom b
   !> @param[inout] hvp2_rr_a  sum_B v_B . d2/(dr_a dr_B) of the order-2 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom a
   !> @param[inout] hvp2_rr_b  sum_B v_B . d2/(dr_b dr_B) of the order-2 spatial derivative;
   !>                          trailing index is the retained nuclear component of atom b
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_hvp_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, fade_width, &
                                           v_a, v_b, max_deriv, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                           hvp0_a, hvp0_b, hvp1_r_a, hvp1_r_b, hvp2_rr_a, &
                                           hvp2_rr_b, sg_a, sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Nuclear displacement direction of atom b
      real(wp), intent(in) :: v_b(3)
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A v_A . d/dr_A of the order-0 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg0
      !> sum_A v_A . d/dr_A of the order-1 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A v_A . d/dr_A of the order-2 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A v_A . d/dr_A of the order-3 spatial derivative (free by-product of the HVP; see the
      !> note above)
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-0 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp0_a(3)
      !> sum_B v_B . d2/(dr_b dr_B) of the order-0 spatial derivative; trailing index is the
      !> retained nuclear component of atom b
      real(wp), intent(inout) :: hvp0_b(3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-1 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp1_r_a(3, 3)
      !> sum_B v_B . d2/(dr_b dr_B) of the order-1 spatial derivative; trailing index is the
      !> retained nuclear component of atom b
      real(wp), intent(inout) :: hvp1_r_b(3, 3)
      !> sum_B v_B . d2/(dr_a dr_B) of the order-2 spatial derivative; trailing index is the
      !> retained nuclear component of atom a
      real(wp), intent(inout) :: hvp2_rr_a(3, 3, 3)
      !> sum_B v_B . d2/(dr_b dr_B) of the order-2 spatial derivative; trailing index is the
      !> retained nuclear component of atom b
      real(wp), intent(inout) :: hvp2_rr_b(3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_hvp_eval_faded

   !> Accumulate one pair's contribution to PD and its spatial derivatives.
   !>
   !> The assembled form of `cfc_pair_term_eval`: the sum over the 2**k
   !> assignments of the spatial indices to the two centers has already been
   !> carried out symbolically, so the caller adds nothing up.  Prefer this
   !> routine whenever only spatial derivatives are wanted -- it is both simpler
   !> to call and cheaper (1079 temporaries at order 4 against 1440 for the 210
   !> split components).
   !>
   !> `cfc_pair_term_eval` remains the routine to call when the *split* d_a/d_b
   !> derivatives themselves are wanted, i.e. for the uncontracted pairwise
   !> nuclear tensors.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Maximum derivative order to accumulate (0..4); a larger value is
   !>                          treated as 4
   !> @param[inout] pd0        Accumulated order-0 spatial derivative of PD
   !> @param[inout] pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[inout] pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[inout] pd3_rrr    Accumulated order-3 spatial derivative of PD
   !> @param[inout] pd4_rrrr   Accumulated order-4 spatial derivative of PD
   pure module subroutine cfc_pair_spatial_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, fade_width, &
                                         max_deriv, pd0, pd1_r, pd2_rr, pd3_rrr, pd4_rrrr)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Maximum derivative order to accumulate (0..4); a larger value is treated as 4
      integer, intent(in) :: max_deriv
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(inout) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(inout) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(inout) :: pd2_rr(3, 3)
      !> Accumulated order-3 spatial derivative of PD
      real(wp), intent(inout) :: pd3_rrr(3, 3, 3)
      !> Accumulated order-4 spatial derivative of PD
      real(wp), intent(inout) :: pd4_rrrr(3, 3, 3, 3)
   end subroutine cfc_pair_spatial_eval

   !> Faded term of `cfc_pair_spatial_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Maximum derivative order to accumulate (0..4); a larger value is
   !>                          treated as 4
   !> @param[inout] pd0        Accumulated order-0 spatial derivative of PD
   !> @param[inout] pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[inout] pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[inout] pd3_rrr    Accumulated order-3 spatial derivative of PD
   !> @param[inout] pd4_rrrr   Accumulated order-4 spatial derivative of PD
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_spatial_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                               fade_width, max_deriv, pd0, pd1_r, pd2_rr, &
                                               pd3_rrr, pd4_rrrr, sg_a, sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Maximum derivative order to accumulate (0..4); a larger value is treated as 4
      integer, intent(in) :: max_deriv
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(inout) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(inout) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(inout) :: pd2_rr(3, 3)
      !> Accumulated order-3 spatial derivative of PD
      real(wp), intent(inout) :: pd3_rrr(3, 3, 3)
      !> Accumulated order-4 spatial derivative of PD
      real(wp), intent(inout) :: pd4_rrrr(3, 3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_spatial_eval_faded

   !> Accumulate the uncontracted nuclear derivative of the spatial derivatives.
   !>
   !> Given no direction vector, this returns, for each center A of this atom,
   !>
   !>    qn{k}_A(i1..ik, beta) = d/dr_A,beta  d^k K / dr_i1..dr_ik
   !>
   !> i.e. one nuclear index, kept.  Summed over an atom's own term and every
   !> pair it belongs to, this is the per-atom family `cfc_hvp_eval` needs (and
   !> the same object the uncontracted `f3_rr_rA` accessor is built from).
   !> Nuclear signs are already applied.
   !>
   !> A caller wanting nuclear derivatives wants the spatial ones too, and runs
   !> `*_spatial_eval` alongside this: the two are separate routines on purpose,
   !> see the note on the builder.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] qn0_a      d/dr_a of the order-0 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn1_r_a    d/dr_a of the order-1 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn2_rr_a   d/dr_a of the order-2 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn3_rrr_a  d/dr_a of the order-3 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   pure module subroutine cfc_atomic_nuclear_eval(d_a, R_a, a1, fade_start, fade_width, max_deriv, &
                                           qn0_a, qn1_r_a, qn2_rr_a, qn3_rrr_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> d/dr_a of the order-0 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn0_a(3)
      !> d/dr_a of the order-1 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn1_r_a(3, 3)
      !> d/dr_a of the order-2 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn2_rr_a(3, 3, 3)
      !> d/dr_a of the order-3 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn3_rrr_a(3, 3, 3, 3)
   end subroutine cfc_atomic_nuclear_eval

   !> Faded term of `cfc_atomic_nuclear_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] qn0_a      d/dr_a of the order-0 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn1_r_a    d/dr_a of the order-1 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn2_rr_a   d/dr_a of the order-2 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn3_rrr_a  d/dr_a of the order-3 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   pure module subroutine cfc_atomic_nuclear_eval_faded(d_a, R_a, a1, fade_start, fade_width, max_deriv, &
                                                 qn0_a, qn1_r_a, qn2_rr_a, qn3_rrr_a, sg_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> d/dr_a of the order-0 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn0_a(3)
      !> d/dr_a of the order-1 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn1_r_a(3, 3)
      !> d/dr_a of the order-2 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn2_rr_a(3, 3, 3)
      !> d/dr_a of the order-3 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn3_rrr_a(3, 3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
   end subroutine cfc_atomic_nuclear_eval_faded

   !> Accumulate the uncontracted nuclear derivative of the spatial derivatives.
   !>
   !> Given no direction vector, this returns, for each center A of this pair,
   !>
   !>    qn{k}_A(i1..ik, beta) = d/dr_A,beta  d^k K / dr_i1..dr_ik
   !>
   !> i.e. one nuclear index, kept.  Summed over an atom's own term and every
   !> pair it belongs to, this is the per-atom family `cfc_hvp_eval` needs (and
   !> the same object the uncontracted `f3_rr_rA` accessor is built from).
   !> Nuclear signs are already applied.
   !>
   !> A caller wanting nuclear derivatives wants the spatial ones too, and runs
   !> `*_spatial_eval` alongside this: the two are separate routines on purpose,
   !> see the note on the builder.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] qn0_a      d/dr_a of the order-0 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn0_b      d/dr_b of the order-0 spatial derivative; trailing index is the
   !>                          nuclear component of atom b
   !> @param[inout] qn1_r_a    d/dr_a of the order-1 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn1_r_b    d/dr_b of the order-1 spatial derivative; trailing index is the
   !>                          nuclear component of atom b
   !> @param[inout] qn2_rr_a   d/dr_a of the order-2 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn2_rr_b   d/dr_b of the order-2 spatial derivative; trailing index is the
   !>                          nuclear component of atom b
   !> @param[inout] qn3_rrr_a  d/dr_a of the order-3 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn3_rrr_b  d/dr_b of the order-3 spatial derivative; trailing index is the
   !>                          nuclear component of atom b
   pure module subroutine cfc_pair_nuclear_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, fade_width, &
                                         max_deriv, qn0_a, qn0_b, qn1_r_a, qn1_r_b, qn2_rr_a, &
                                         qn2_rr_b, qn3_rrr_a, qn3_rrr_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> d/dr_a of the order-0 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn0_a(3)
      !> d/dr_b of the order-0 spatial derivative; trailing index is the nuclear component of atom b
      real(wp), intent(inout) :: qn0_b(3)
      !> d/dr_a of the order-1 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn1_r_a(3, 3)
      !> d/dr_b of the order-1 spatial derivative; trailing index is the nuclear component of atom b
      real(wp), intent(inout) :: qn1_r_b(3, 3)
      !> d/dr_a of the order-2 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn2_rr_a(3, 3, 3)
      !> d/dr_b of the order-2 spatial derivative; trailing index is the nuclear component of atom b
      real(wp), intent(inout) :: qn2_rr_b(3, 3, 3)
      !> d/dr_a of the order-3 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn3_rrr_a(3, 3, 3, 3)
      !> d/dr_b of the order-3 spatial derivative; trailing index is the nuclear component of atom b
      real(wp), intent(inout) :: qn3_rrr_b(3, 3, 3, 3)
   end subroutine cfc_pair_nuclear_eval

   !> Faded term of `cfc_pair_nuclear_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] qn0_a      d/dr_a of the order-0 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn0_b      d/dr_b of the order-0 spatial derivative; trailing index is the
   !>                          nuclear component of atom b
   !> @param[inout] qn1_r_a    d/dr_a of the order-1 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn1_r_b    d/dr_b of the order-1 spatial derivative; trailing index is the
   !>                          nuclear component of atom b
   !> @param[inout] qn2_rr_a   d/dr_a of the order-2 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn2_rr_b   d/dr_b of the order-2 spatial derivative; trailing index is the
   !>                          nuclear component of atom b
   !> @param[inout] qn3_rrr_a  d/dr_a of the order-3 spatial derivative; trailing index is the
   !>                          nuclear component of atom a
   !> @param[inout] qn3_rrr_b  d/dr_b of the order-3 spatial derivative; trailing index is the
   !>                          nuclear component of atom b
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_nuclear_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                               fade_width, max_deriv, qn0_a, qn0_b, qn1_r_a, &
                                               qn1_r_b, qn2_rr_a, qn2_rr_b, qn3_rrr_a, qn3_rrr_b, &
                                               sg_a, sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> d/dr_a of the order-0 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn0_a(3)
      !> d/dr_b of the order-0 spatial derivative; trailing index is the nuclear component of atom b
      real(wp), intent(inout) :: qn0_b(3)
      !> d/dr_a of the order-1 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn1_r_a(3, 3)
      !> d/dr_b of the order-1 spatial derivative; trailing index is the nuclear component of atom b
      real(wp), intent(inout) :: qn1_r_b(3, 3)
      !> d/dr_a of the order-2 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn2_rr_a(3, 3, 3)
      !> d/dr_b of the order-2 spatial derivative; trailing index is the nuclear component of atom b
      real(wp), intent(inout) :: qn2_rr_b(3, 3, 3)
      !> d/dr_a of the order-3 spatial derivative; trailing index is the nuclear component of atom a
      real(wp), intent(inout) :: qn3_rrr_a(3, 3, 3, 3)
      !> d/dr_b of the order-3 spatial derivative; trailing index is the nuclear component of atom b
      real(wp), intent(inout) :: qn3_rrr_b(3, 3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_nuclear_eval_faded

   !> Accumulate the uncontracted two-nucleus derivative of the spatial ones.
   !>
   !> Given no direction vector, this returns, for each ordered pair of centers
   !> (A, B) this atom can supply,
   !>
   !>    qq{k}_AB(i1..ik, beta, gamma) = d2/(dr_A,beta dr_B,gamma)
   !>                                      d^k K / dr_i1..dr_ik
   !>
   !> i.e. two nuclear indices, both kept.  This is the uncontracted sibling of
   !> `*_hvp_eval`, and the object `cfc_hessian_eval` lifts into `f2_rArB`,
   !> `f3_r_rArB` and `f4_rr_rArB`.  Nuclear signs are already applied (the two
   !> flips cancel).
   !>
   !> Which block a caller wants follows from the atom pair it is filling.  For
   !> an off-diagonal (A, B) only the pair term containing *both* survives two
   !> different nuclear derivatives, so a single call with A in slot a and B in
   !> slot b supplies everything through the `_ab` block.  For a diagonal
   !> (A, A) the contributions are atom A's own term plus the `_aa`/`_bb` block
   !> of every pair containing A, which is one sweep over the partners.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] qq0_aa     d2/(dr_a dr_a) of the order-0 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[inout] qq1_r_aa   d2/(dr_a dr_a) of the order-1 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[inout] qq2_rr_aa  d2/(dr_a dr_a) of the order-2 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   pure module subroutine cfc_atomic_hessian_eval(d_a, R_a, a1, fade_start, fade_width, max_deriv, &
                                           qq0_aa, qq1_r_aa, qq2_rr_aa)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dr_a dr_a) of the order-0 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq0_aa(3, 3)
      !> d2/(dr_a dr_a) of the order-1 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq1_r_aa(3, 3, 3)
      !> d2/(dr_a dr_a) of the order-2 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq2_rr_aa(3, 3, 3, 3)
   end subroutine cfc_atomic_hessian_eval

   !> Faded term of `cfc_atomic_hessian_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] qq0_aa     d2/(dr_a dr_a) of the order-0 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[inout] qq1_r_aa   d2/(dr_a dr_a) of the order-1 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[inout] qq2_rr_aa  d2/(dr_a dr_a) of the order-2 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   pure module subroutine cfc_atomic_hessian_eval_faded(d_a, R_a, a1, fade_start, fade_width, max_deriv, &
                                                 qq0_aa, qq1_r_aa, qq2_rr_aa, sg_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dr_a dr_a) of the order-0 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq0_aa(3, 3)
      !> d2/(dr_a dr_a) of the order-1 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq1_r_aa(3, 3, 3)
      !> d2/(dr_a dr_a) of the order-2 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq2_rr_aa(3, 3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
   end subroutine cfc_atomic_hessian_eval_faded

   !> Accumulate the uncontracted two-nucleus derivative of the spatial ones.
   !>
   !> Given no direction vector, this returns, for each ordered pair of centers
   !> (A, B) this pair can supply,
   !>
   !>    qq{k}_AB(i1..ik, beta, gamma) = d2/(dr_A,beta dr_B,gamma)
   !>                                      d^k K / dr_i1..dr_ik
   !>
   !> i.e. two nuclear indices, both kept.  This is the uncontracted sibling of
   !> `*_hvp_eval`, and the object `cfc_hessian_eval` lifts into `f2_rArB`,
   !> `f3_r_rArB` and `f4_rr_rArB`.  Nuclear signs are already applied (the two
   !> flips cancel).
   !>
   !> Which block a caller wants follows from the atom pair it is filling.  For
   !> an off-diagonal (A, B) only the pair term containing *both* survives two
   !> different nuclear derivatives, so a single call with A in slot a and B in
   !> slot b supplies everything through the `_ab` block.  For a diagonal
   !> (A, A) the contributions are atom A's own term plus the `_aa`/`_bb` block
   !> of every pair containing A, which is one sweep over the partners.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] qq0_aa     d2/(dr_a dr_a) of the order-0 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[inout] qq0_ab     d2/(dr_a dr_b) of the order-0 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and b
   !> @param[inout] qq0_bb     d2/(dr_b dr_b) of the order-0 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms b and b
   !> @param[inout] qq1_r_aa   d2/(dr_a dr_a) of the order-1 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[inout] qq1_r_ab   d2/(dr_a dr_b) of the order-1 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and b
   !> @param[inout] qq1_r_bb   d2/(dr_b dr_b) of the order-1 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms b and b
   !> @param[inout] qq2_rr_aa  d2/(dr_a dr_a) of the order-2 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[inout] qq2_rr_ab  d2/(dr_a dr_b) of the order-2 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and b
   !> @param[inout] qq2_rr_bb  d2/(dr_b dr_b) of the order-2 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms b and b
   pure module subroutine cfc_pair_hessian_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, fade_width, &
                                         max_deriv, qq0_aa, qq0_ab, qq0_bb, qq1_r_aa, qq1_r_ab, &
                                         qq1_r_bb, qq2_rr_aa, qq2_rr_ab, qq2_rr_bb)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dr_a dr_a) of the order-0 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq0_aa(3, 3)
      !> d2/(dr_a dr_b) of the order-0 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and b
      real(wp), intent(inout) :: qq0_ab(3, 3)
      !> d2/(dr_b dr_b) of the order-0 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms b and b
      real(wp), intent(inout) :: qq0_bb(3, 3)
      !> d2/(dr_a dr_a) of the order-1 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq1_r_aa(3, 3, 3)
      !> d2/(dr_a dr_b) of the order-1 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and b
      real(wp), intent(inout) :: qq1_r_ab(3, 3, 3)
      !> d2/(dr_b dr_b) of the order-1 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms b and b
      real(wp), intent(inout) :: qq1_r_bb(3, 3, 3)
      !> d2/(dr_a dr_a) of the order-2 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq2_rr_aa(3, 3, 3, 3)
      !> d2/(dr_a dr_b) of the order-2 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and b
      real(wp), intent(inout) :: qq2_rr_ab(3, 3, 3, 3)
      !> d2/(dr_b dr_b) of the order-2 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms b and b
      real(wp), intent(inout) :: qq2_rr_bb(3, 3, 3, 3)
   end subroutine cfc_pair_hessian_eval

   !> Faded term of `cfc_pair_hessian_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] qq0_aa     d2/(dr_a dr_a) of the order-0 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[inout] qq0_ab     d2/(dr_a dr_b) of the order-0 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and b
   !> @param[inout] qq0_bb     d2/(dr_b dr_b) of the order-0 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms b and b
   !> @param[inout] qq1_r_aa   d2/(dr_a dr_a) of the order-1 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[inout] qq1_r_ab   d2/(dr_a dr_b) of the order-1 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and b
   !> @param[inout] qq1_r_bb   d2/(dr_b dr_b) of the order-1 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms b and b
   !> @param[inout] qq2_rr_aa  d2/(dr_a dr_a) of the order-2 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and a
   !> @param[inout] qq2_rr_ab  d2/(dr_a dr_b) of the order-2 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms a and b
   !> @param[inout] qq2_rr_bb  d2/(dr_b dr_b) of the order-2 spatial derivative; the two trailing
   !>                          indices are the nuclear components of atoms b and b
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_hessian_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                               fade_width, max_deriv, qq0_aa, qq0_ab, qq0_bb, &
                                               qq1_r_aa, qq1_r_ab, qq1_r_bb, qq2_rr_aa, &
                                               qq2_rr_ab, qq2_rr_bb, sg_a, sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dr_a dr_a) of the order-0 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq0_aa(3, 3)
      !> d2/(dr_a dr_b) of the order-0 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and b
      real(wp), intent(inout) :: qq0_ab(3, 3)
      !> d2/(dr_b dr_b) of the order-0 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms b and b
      real(wp), intent(inout) :: qq0_bb(3, 3)
      !> d2/(dr_a dr_a) of the order-1 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq1_r_aa(3, 3, 3)
      !> d2/(dr_a dr_b) of the order-1 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and b
      real(wp), intent(inout) :: qq1_r_ab(3, 3, 3)
      !> d2/(dr_b dr_b) of the order-1 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms b and b
      real(wp), intent(inout) :: qq1_r_bb(3, 3, 3)
      !> d2/(dr_a dr_a) of the order-2 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and a
      real(wp), intent(inout) :: qq2_rr_aa(3, 3, 3, 3)
      !> d2/(dr_a dr_b) of the order-2 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms a and b
      real(wp), intent(inout) :: qq2_rr_ab(3, 3, 3, 3)
      !> d2/(dr_b dr_b) of the order-2 spatial derivative; the two trailing indices are the nuclear
      !> components of atoms b and b
      real(wp), intent(inout) :: qq2_rr_bb(3, 3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_hessian_eval_faded

   !> Accumulate the radius derivative of the spatial derivatives.
   !>
   !> Given no direction vector, this returns, for each center A of this atom,
   !>
   !>    rd{k}_A(i1..ik) = d/dR_A  d^k K / dr_i1..dr_ik
   !>
   !> where ``R_A`` is atom A's *cavity radius*, a parameter of the level set
   !> rather than a coordinate.  Three consequences, all visible in the
   !> signature: no direction vector, no sign flip (the displacement is
   !> ``r - r_a``, so positions flip and radii do not), and no trailing index --
   !> a radius is a scalar, so the result has exactly the rank of the spatial
   !> tensor it came from.  Contrast `*_nuclear_eval`, whose output grows one
   !> Cartesian index.
   !>
   !> Summed over an atom's own term and every pair it belongs to, this is the
   !> per-atom family `cfc_radius_eval` lifts into `f1_rad` .. `f4_rrr_rad`.
   !> A caller wanting radius derivatives wants the spatial ones too, and runs
   !> `*_spatial_eval` alongside this.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] rd0_a      d/dR_a of the order-0 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd1_r_a    d/dR_a of the order-1 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd2_rr_a   d/dR_a of the order-2 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd3_rrr_a  d/dR_a of the order-3 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   pure module subroutine cfc_atomic_radius_eval(d_a, R_a, a1, fade_start, fade_width, max_deriv, rd0_a, &
                                          rd1_r_a, rd2_rr_a, rd3_rrr_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> d/dR_a of the order-0 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd0_a
      !> d/dR_a of the order-1 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd1_r_a(3)
      !> d/dR_a of the order-2 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd2_rr_a(3, 3)
      !> d/dR_a of the order-3 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd3_rrr_a(3, 3, 3)
   end subroutine cfc_atomic_radius_eval

   !> Faded term of `cfc_atomic_radius_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] rd0_a      d/dR_a of the order-0 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd1_r_a    d/dR_a of the order-1 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd2_rr_a   d/dR_a of the order-2 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd3_rrr_a  d/dR_a of the order-3 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   pure module subroutine cfc_atomic_radius_eval_faded(d_a, R_a, a1, fade_start, fade_width, max_deriv, &
                                                rd0_a, rd1_r_a, rd2_rr_a, rd3_rrr_a, sg_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> d/dR_a of the order-0 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd0_a
      !> d/dR_a of the order-1 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd1_r_a(3)
      !> d/dR_a of the order-2 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd2_rr_a(3, 3)
      !> d/dR_a of the order-3 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd3_rrr_a(3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
   end subroutine cfc_atomic_radius_eval_faded

   !> Accumulate the radius derivative of the spatial derivatives.
   !>
   !> Given no direction vector, this returns, for each center A of this pair,
   !>
   !>    rd{k}_A(i1..ik) = d/dR_A  d^k K / dr_i1..dr_ik
   !>
   !> where ``R_A`` is atom A's *cavity radius*, a parameter of the level set
   !> rather than a coordinate.  Three consequences, all visible in the
   !> signature: no direction vector, no sign flip (the displacement is
   !> ``r - r_a``, so positions flip and radii do not), and no trailing index --
   !> a radius is a scalar, so the result has exactly the rank of the spatial
   !> tensor it came from.  Contrast `*_nuclear_eval`, whose output grows one
   !> Cartesian index.
   !>
   !> Summed over an atom's own term and every pair it belongs to, this is the
   !> per-atom family `cfc_radius_eval` lifts into `f1_rad` .. `f4_rrr_rad`.
   !> A caller wanting radius derivatives wants the spatial ones too, and runs
   !> `*_spatial_eval` alongside this.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] rd0_a      d/dR_a of the order-0 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd0_b      d/dR_b of the order-0 spatial derivative; R_b is atom b's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd1_r_a    d/dR_a of the order-1 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd1_r_b    d/dR_b of the order-1 spatial derivative; R_b is atom b's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd2_rr_a   d/dR_a of the order-2 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd2_rr_b   d/dR_b of the order-2 spatial derivative; R_b is atom b's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd3_rrr_a  d/dR_a of the order-3 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd3_rrr_b  d/dR_b of the order-3 spatial derivative; R_b is atom b's radius, so
   !>                          there is no trailing nuclear index
   pure module subroutine cfc_pair_radius_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, fade_width, &
                                        max_deriv, rd0_a, rd0_b, rd1_r_a, rd1_r_b, rd2_rr_a, &
                                        rd2_rr_b, rd3_rrr_a, rd3_rrr_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> d/dR_a of the order-0 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd0_a
      !> d/dR_b of the order-0 spatial derivative; R_b is atom b's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd0_b
      !> d/dR_a of the order-1 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd1_r_a(3)
      !> d/dR_b of the order-1 spatial derivative; R_b is atom b's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd1_r_b(3)
      !> d/dR_a of the order-2 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd2_rr_a(3, 3)
      !> d/dR_b of the order-2 spatial derivative; R_b is atom b's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd2_rr_b(3, 3)
      !> d/dR_a of the order-3 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd3_rrr_a(3, 3, 3)
      !> d/dR_b of the order-3 spatial derivative; R_b is atom b's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd3_rrr_b(3, 3, 3)
   end subroutine cfc_pair_radius_eval

   !> Faded term of `cfc_pair_radius_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] rd0_a      d/dR_a of the order-0 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd0_b      d/dR_b of the order-0 spatial derivative; R_b is atom b's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd1_r_a    d/dR_a of the order-1 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd1_r_b    d/dR_b of the order-1 spatial derivative; R_b is atom b's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd2_rr_a   d/dR_a of the order-2 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd2_rr_b   d/dR_b of the order-2 spatial derivative; R_b is atom b's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd3_rrr_a  d/dR_a of the order-3 spatial derivative; R_a is atom a's radius, so
   !>                          there is no trailing nuclear index
   !> @param[inout] rd3_rrr_b  d/dR_b of the order-3 spatial derivative; R_b is atom b's radius, so
   !>                          there is no trailing nuclear index
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_radius_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                              fade_width, max_deriv, rd0_a, rd0_b, rd1_r_a, &
                                              rd1_r_b, rd2_rr_a, rd2_rr_b, rd3_rrr_a, rd3_rrr_b, &
                                              sg_a, sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> d/dR_a of the order-0 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd0_a
      !> d/dR_b of the order-0 spatial derivative; R_b is atom b's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd0_b
      !> d/dR_a of the order-1 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd1_r_a(3)
      !> d/dR_b of the order-1 spatial derivative; R_b is atom b's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd1_r_b(3)
      !> d/dR_a of the order-2 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd2_rr_a(3, 3)
      !> d/dR_b of the order-2 spatial derivative; R_b is atom b's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd2_rr_b(3, 3)
      !> d/dR_a of the order-3 spatial derivative; R_a is atom a's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd3_rrr_a(3, 3, 3)
      !> d/dR_b of the order-3 spatial derivative; R_b is atom b's radius, so there is no trailing
      !> nuclear index
      real(wp), intent(inout) :: rd3_rrr_b(3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_radius_eval_faded

   !> Accumulate the uncontracted two-*radius* derivative of the spatial ones.
   !>
   !> The radius twin of `*_hessian_eval`.  For each pair of centers (A, B) this
   !> atom can supply,
   !>
   !>    rr{k}_AB(i1..ik) = d2/(dRad_A dRad_B)  d^k K / dr_i1..dr_ik
   !>
   !> Two retained radii, both kept, neither adding an index: the results have the
   !> same rank as the spatial tensors they come from, two below the `qq*` family,
   !> and no nuclear sign is involved because a radius is not a coordinate.
   !>
   !> `cfc_radius2_eval` lifts these into `f2_radrad`, `f3_r_radrad` and
   !> `f4_rr_radrad`.  Which block a caller wants follows from the atom pair it is
   !> filling, exactly as for `*_hessian_eval`: `_ab` for an off-diagonal (A, B),
   !> `_aa`/`_bb` swept over the partners for a diagonal (A, A).
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] rr0_aa     d2/(dR_a dR_a) of the order-0 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr1_r_aa   d2/(dR_a dR_a) of the order-1 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr2_rr_aa  d2/(dR_a dR_a) of the order-2 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   pure module subroutine cfc_atomic_radius2_eval(d_a, R_a, a1, fade_start, fade_width, max_deriv, &
                                           rr0_aa, rr1_r_aa, rr2_rr_aa)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dR_a dR_a) of the order-0 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr0_aa
      !> d2/(dR_a dR_a) of the order-1 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr1_r_aa(3)
      !> d2/(dR_a dR_a) of the order-2 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr2_rr_aa(3, 3)
   end subroutine cfc_atomic_radius2_eval

   !> Faded term of `cfc_atomic_radius2_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] rr0_aa     d2/(dR_a dR_a) of the order-0 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr1_r_aa   d2/(dR_a dR_a) of the order-1 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr2_rr_aa  d2/(dR_a dR_a) of the order-2 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   pure module subroutine cfc_atomic_radius2_eval_faded(d_a, R_a, a1, fade_start, fade_width, max_deriv, &
                                                 rr0_aa, rr1_r_aa, rr2_rr_aa, sg_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dR_a dR_a) of the order-0 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr0_aa
      !> d2/(dR_a dR_a) of the order-1 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr1_r_aa(3)
      !> d2/(dR_a dR_a) of the order-2 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr2_rr_aa(3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
   end subroutine cfc_atomic_radius2_eval_faded

   !> Accumulate the uncontracted two-*radius* derivative of the spatial ones.
   !>
   !> The radius twin of `*_hessian_eval`.  For each pair of centers (A, B) this
   !> pair can supply,
   !>
   !>    rr{k}_AB(i1..ik) = d2/(dRad_A dRad_B)  d^k K / dr_i1..dr_ik
   !>
   !> Two retained radii, both kept, neither adding an index: the results have the
   !> same rank as the spatial tensors they come from, two below the `qq*` family,
   !> and no nuclear sign is involved because a radius is not a coordinate.
   !>
   !> `cfc_radius2_eval` lifts these into `f2_radrad`, `f3_r_radrad` and
   !> `f4_rr_radrad`.  Which block a caller wants follows from the atom pair it is
   !> filling, exactly as for `*_hessian_eval`: `_ab` for an off-diagonal (A, B),
   !> `_aa`/`_bb` swept over the partners for a diagonal (A, A).
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] rr0_aa     d2/(dR_a dR_a) of the order-0 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr0_ab     d2/(dR_a dR_b) of the order-0 spatial derivative, where R_a and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr0_bb     d2/(dR_b dR_b) of the order-0 spatial derivative, where R_b and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr1_r_aa   d2/(dR_a dR_a) of the order-1 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr1_r_ab   d2/(dR_a dR_b) of the order-1 spatial derivative, where R_a and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr1_r_bb   d2/(dR_b dR_b) of the order-1 spatial derivative, where R_b and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr2_rr_aa  d2/(dR_a dR_a) of the order-2 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr2_rr_ab  d2/(dR_a dR_b) of the order-2 spatial derivative, where R_a and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr2_rr_bb  d2/(dR_b dR_b) of the order-2 spatial derivative, where R_b and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   pure module subroutine cfc_pair_radius2_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, fade_width, &
                                         max_deriv, rr0_aa, rr0_ab, rr0_bb, rr1_r_aa, rr1_r_ab, &
                                         rr1_r_bb, rr2_rr_aa, rr2_rr_ab, rr2_rr_bb)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dR_a dR_a) of the order-0 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr0_aa
      !> d2/(dR_a dR_b) of the order-0 spatial derivative, where R_a and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr0_ab
      !> d2/(dR_b dR_b) of the order-0 spatial derivative, where R_b and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr0_bb
      !> d2/(dR_a dR_a) of the order-1 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr1_r_aa(3)
      !> d2/(dR_a dR_b) of the order-1 spatial derivative, where R_a and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr1_r_ab(3)
      !> d2/(dR_b dR_b) of the order-1 spatial derivative, where R_b and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr1_r_bb(3)
      !> d2/(dR_a dR_a) of the order-2 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr2_rr_aa(3, 3)
      !> d2/(dR_a dR_b) of the order-2 spatial derivative, where R_a and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr2_rr_ab(3, 3)
      !> d2/(dR_b dR_b) of the order-2 spatial derivative, where R_b and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr2_rr_bb(3, 3)
   end subroutine cfc_pair_radius2_eval

   !> Faded term of `cfc_pair_radius2_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] rr0_aa     d2/(dR_a dR_a) of the order-0 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr0_ab     d2/(dR_a dR_b) of the order-0 spatial derivative, where R_a and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr0_bb     d2/(dR_b dR_b) of the order-0 spatial derivative, where R_b and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr1_r_aa   d2/(dR_a dR_a) of the order-1 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr1_r_ab   d2/(dR_a dR_b) of the order-1 spatial derivative, where R_a and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr1_r_bb   d2/(dR_b dR_b) of the order-1 spatial derivative, where R_b and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr2_rr_aa  d2/(dR_a dR_a) of the order-2 spatial derivative, where R_a and R_a
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr2_rr_ab  d2/(dR_a dR_b) of the order-2 spatial derivative, where R_a and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[inout] rr2_rr_bb  d2/(dR_b dR_b) of the order-2 spatial derivative, where R_b and R_b
   !>                          are the atoms' radii; no trailing index, radii are scalars
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_radius2_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                               fade_width, max_deriv, rr0_aa, rr0_ab, rr0_bb, &
                                               rr1_r_aa, rr1_r_ab, rr1_r_bb, rr2_rr_aa, &
                                               rr2_rr_ab, rr2_rr_bb, sg_a, sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dR_a dR_a) of the order-0 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr0_aa
      !> d2/(dR_a dR_b) of the order-0 spatial derivative, where R_a and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr0_ab
      !> d2/(dR_b dR_b) of the order-0 spatial derivative, where R_b and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr0_bb
      !> d2/(dR_a dR_a) of the order-1 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr1_r_aa(3)
      !> d2/(dR_a dR_b) of the order-1 spatial derivative, where R_a and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr1_r_ab(3)
      !> d2/(dR_b dR_b) of the order-1 spatial derivative, where R_b and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr1_r_bb(3)
      !> d2/(dR_a dR_a) of the order-2 spatial derivative, where R_a and R_a are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr2_rr_aa(3, 3)
      !> d2/(dR_a dR_b) of the order-2 spatial derivative, where R_a and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr2_rr_ab(3, 3)
      !> d2/(dR_b dR_b) of the order-2 spatial derivative, where R_b and R_b are the atoms' radii;
      !> no trailing index, radii are scalars
      real(wp), intent(inout) :: rr2_rr_bb(3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_radius2_eval_faded

   !> Accumulate the uncontracted mixed nuclear-radius derivative.
   !>
   !> For each ordered pair of centers (A, B) this atom can supply,
   !>
   !>    mx{k}_AB(i1..ik, beta) = d2/(dr_A,beta dRad_B)  d^k K / dr_i1..dr_ik
   !>
   !> The first slot letter names the center of the *position* derivative, the
   !> second the center of the *radius* derivative.  Those are different kinds of
   !> derivative, so the pair term emits all four combinations rather than three:
   !> `_ab` and `_ba` are distinct expressions, not transposes of one another.
   !>
   !> `cfc_nucrad_eval` lifts these into `f2_rA_rad`, `f3_r_rA_rad` and
   !> `f4_rr_rA_rad`.  The nuclear sign is already applied.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] mx0_aa     d2/(dr_a dR_a) of the order-0 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[inout] mx1_r_aa   d2/(dr_a dR_a) of the order-1 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[inout] mx2_rr_aa  d2/(dr_a dR_a) of the order-2 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   pure module subroutine cfc_atomic_nucrad_eval(d_a, R_a, a1, fade_start, fade_width, max_deriv, &
                                          mx0_aa, mx1_r_aa, mx2_rr_aa)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dr_a dR_a) of the order-0 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx0_aa(3)
      !> d2/(dr_a dR_a) of the order-1 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx1_r_aa(3, 3)
      !> d2/(dr_a dR_a) of the order-2 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx2_rr_aa(3, 3, 3)
   end subroutine cfc_atomic_nucrad_eval

   !> Faded term of `cfc_atomic_nucrad_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] mx0_aa     d2/(dr_a dR_a) of the order-0 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[inout] mx1_r_aa   d2/(dr_a dR_a) of the order-1 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[inout] mx2_rr_aa  d2/(dr_a dR_a) of the order-2 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   pure module subroutine cfc_atomic_nucrad_eval_faded(d_a, R_a, a1, fade_start, fade_width, max_deriv, &
                                                mx0_aa, mx1_r_aa, mx2_rr_aa, sg_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dr_a dR_a) of the order-0 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx0_aa(3)
      !> d2/(dr_a dR_a) of the order-1 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx1_r_aa(3, 3)
      !> d2/(dr_a dR_a) of the order-2 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx2_rr_aa(3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
   end subroutine cfc_atomic_nucrad_eval_faded

   !> Accumulate the uncontracted mixed nuclear-radius derivative.
   !>
   !> For each ordered pair of centers (A, B) this pair can supply,
   !>
   !>    mx{k}_AB(i1..ik, beta) = d2/(dr_A,beta dRad_B)  d^k K / dr_i1..dr_ik
   !>
   !> The first slot letter names the center of the *position* derivative, the
   !> second the center of the *radius* derivative.  Those are different kinds of
   !> derivative, so the pair term emits all four combinations rather than three:
   !> `_ab` and `_ba` are distinct expressions, not transposes of one another.
   !>
   !> `cfc_nucrad_eval` lifts these into `f2_rA_rad`, `f3_r_rA_rad` and
   !> `f4_rr_rA_rad`.  The nuclear sign is already applied.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] mx0_aa     d2/(dr_a dR_a) of the order-0 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[inout] mx0_ab     d2/(dr_a dR_b) of the order-0 spatial derivative, position on atom a
   !>                          and radius on atom b; trailing index is a's nuclear component
   !> @param[inout] mx0_ba     d2/(dr_b dR_a) of the order-0 spatial derivative, position on atom b
   !>                          and radius on atom a; trailing index is b's nuclear component
   !> @param[inout] mx0_bb     d2/(dr_b dR_b) of the order-0 spatial derivative, position on atom b
   !>                          and radius on atom b; trailing index is b's nuclear component
   !> @param[inout] mx1_r_aa   d2/(dr_a dR_a) of the order-1 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[inout] mx1_r_ab   d2/(dr_a dR_b) of the order-1 spatial derivative, position on atom a
   !>                          and radius on atom b; trailing index is a's nuclear component
   !> @param[inout] mx1_r_ba   d2/(dr_b dR_a) of the order-1 spatial derivative, position on atom b
   !>                          and radius on atom a; trailing index is b's nuclear component
   !> @param[inout] mx1_r_bb   d2/(dr_b dR_b) of the order-1 spatial derivative, position on atom b
   !>                          and radius on atom b; trailing index is b's nuclear component
   !> @param[inout] mx2_rr_aa  d2/(dr_a dR_a) of the order-2 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[inout] mx2_rr_ab  d2/(dr_a dR_b) of the order-2 spatial derivative, position on atom a
   !>                          and radius on atom b; trailing index is a's nuclear component
   !> @param[inout] mx2_rr_ba  d2/(dr_b dR_a) of the order-2 spatial derivative, position on atom b
   !>                          and radius on atom a; trailing index is b's nuclear component
   !> @param[inout] mx2_rr_bb  d2/(dr_b dR_b) of the order-2 spatial derivative, position on atom b
   !>                          and radius on atom b; trailing index is b's nuclear component
   pure module subroutine cfc_pair_nucrad_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, fade_width, &
                                        max_deriv, mx0_aa, mx0_ab, mx0_ba, mx0_bb, mx1_r_aa, &
                                        mx1_r_ab, mx1_r_ba, mx1_r_bb, mx2_rr_aa, mx2_rr_ab, &
                                        mx2_rr_ba, mx2_rr_bb)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dr_a dR_a) of the order-0 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx0_aa(3)
      !> d2/(dr_a dR_b) of the order-0 spatial derivative, position on atom a and radius on atom b;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx0_ab(3)
      !> d2/(dr_b dR_a) of the order-0 spatial derivative, position on atom b and radius on atom a;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx0_ba(3)
      !> d2/(dr_b dR_b) of the order-0 spatial derivative, position on atom b and radius on atom b;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx0_bb(3)
      !> d2/(dr_a dR_a) of the order-1 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx1_r_aa(3, 3)
      !> d2/(dr_a dR_b) of the order-1 spatial derivative, position on atom a and radius on atom b;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx1_r_ab(3, 3)
      !> d2/(dr_b dR_a) of the order-1 spatial derivative, position on atom b and radius on atom a;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx1_r_ba(3, 3)
      !> d2/(dr_b dR_b) of the order-1 spatial derivative, position on atom b and radius on atom b;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx1_r_bb(3, 3)
      !> d2/(dr_a dR_a) of the order-2 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx2_rr_aa(3, 3, 3)
      !> d2/(dr_a dR_b) of the order-2 spatial derivative, position on atom a and radius on atom b;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx2_rr_ab(3, 3, 3)
      !> d2/(dr_b dR_a) of the order-2 spatial derivative, position on atom b and radius on atom a;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx2_rr_ba(3, 3, 3)
      !> d2/(dr_b dR_b) of the order-2 spatial derivative, position on atom b and radius on atom b;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx2_rr_bb(3, 3, 3)
   end subroutine cfc_pair_nucrad_eval

   !> Faded term of `cfc_pair_nucrad_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..2); a larger value
   !>                          is treated as 2
   !> @param[inout] mx0_aa     d2/(dr_a dR_a) of the order-0 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[inout] mx0_ab     d2/(dr_a dR_b) of the order-0 spatial derivative, position on atom a
   !>                          and radius on atom b; trailing index is a's nuclear component
   !> @param[inout] mx0_ba     d2/(dr_b dR_a) of the order-0 spatial derivative, position on atom b
   !>                          and radius on atom a; trailing index is b's nuclear component
   !> @param[inout] mx0_bb     d2/(dr_b dR_b) of the order-0 spatial derivative, position on atom b
   !>                          and radius on atom b; trailing index is b's nuclear component
   !> @param[inout] mx1_r_aa   d2/(dr_a dR_a) of the order-1 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[inout] mx1_r_ab   d2/(dr_a dR_b) of the order-1 spatial derivative, position on atom a
   !>                          and radius on atom b; trailing index is a's nuclear component
   !> @param[inout] mx1_r_ba   d2/(dr_b dR_a) of the order-1 spatial derivative, position on atom b
   !>                          and radius on atom a; trailing index is b's nuclear component
   !> @param[inout] mx1_r_bb   d2/(dr_b dR_b) of the order-1 spatial derivative, position on atom b
   !>                          and radius on atom b; trailing index is b's nuclear component
   !> @param[inout] mx2_rr_aa  d2/(dr_a dR_a) of the order-2 spatial derivative, position on atom a
   !>                          and radius on atom a; trailing index is a's nuclear component
   !> @param[inout] mx2_rr_ab  d2/(dr_a dR_b) of the order-2 spatial derivative, position on atom a
   !>                          and radius on atom b; trailing index is a's nuclear component
   !> @param[inout] mx2_rr_ba  d2/(dr_b dR_a) of the order-2 spatial derivative, position on atom b
   !>                          and radius on atom a; trailing index is b's nuclear component
   !> @param[inout] mx2_rr_bb  d2/(dr_b dR_b) of the order-2 spatial derivative, position on atom b
   !>                          and radius on atom b; trailing index is b's nuclear component
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_nucrad_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                              fade_width, max_deriv, mx0_aa, mx0_ab, mx0_ba, &
                                              mx0_bb, mx1_r_aa, mx1_r_ab, mx1_r_ba, mx1_r_bb, &
                                              mx2_rr_aa, mx2_rr_ab, mx2_rr_ba, mx2_rr_bb, sg_a, &
                                              sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Highest spatial-derivative order to accumulate (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dr_a dR_a) of the order-0 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx0_aa(3)
      !> d2/(dr_a dR_b) of the order-0 spatial derivative, position on atom a and radius on atom b;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx0_ab(3)
      !> d2/(dr_b dR_a) of the order-0 spatial derivative, position on atom b and radius on atom a;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx0_ba(3)
      !> d2/(dr_b dR_b) of the order-0 spatial derivative, position on atom b and radius on atom b;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx0_bb(3)
      !> d2/(dr_a dR_a) of the order-1 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx1_r_aa(3, 3)
      !> d2/(dr_a dR_b) of the order-1 spatial derivative, position on atom a and radius on atom b;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx1_r_ab(3, 3)
      !> d2/(dr_b dR_a) of the order-1 spatial derivative, position on atom b and radius on atom a;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx1_r_ba(3, 3)
      !> d2/(dr_b dR_b) of the order-1 spatial derivative, position on atom b and radius on atom b;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx1_r_bb(3, 3)
      !> d2/(dr_a dR_a) of the order-2 spatial derivative, position on atom a and radius on atom a;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx2_rr_aa(3, 3, 3)
      !> d2/(dr_a dR_b) of the order-2 spatial derivative, position on atom a and radius on atom b;
      !> trailing index is a's nuclear component
      real(wp), intent(inout) :: mx2_rr_ab(3, 3, 3)
      !> d2/(dr_b dR_a) of the order-2 spatial derivative, position on atom b and radius on atom a;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx2_rr_ba(3, 3, 3)
      !> d2/(dr_b dR_b) of the order-2 spatial derivative, position on atom b and radius on atom b;
      !> trailing index is b's nuclear component
      real(wp), intent(inout) :: mx2_rr_bb(3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_nucrad_eval_faded

   !> Accumulate the *radius half* of the nuclear directional derivative.
   !>
   !> Given one radius direction per atom, this ADDS
   !>
   !>    tg{k}(i1..ik) += sum_A vr_A  d/dR_A  d^k K / dr_i1..dr_ik
   !>
   !> into the same `tg*` accumulators `*_tangent_eval` fills with the position
   !> half.  Together the two make the joint contraction
   !> ``sum_A (v_A . d/dr_A + vr_A d/dR_A)``, which is what `cfc_tangent_eval`
   !> and `cfc_hvp_eval` consume; the lift algebra is identical either way,
   !> because both halves land in the same family.
   !>
   !> Split rather than fused on purpose.  This layer is an O(n**2) pair loop
   !> run at every grid point, and the common caller holds its radii fixed:
   !> fusing would charge that caller for terms it discards on every pair, and
   !> would perturb its results in the last bits by regrouping the CSE block.
   !> A fixed-radius caller therefore simply never calls this routine.
   !>
   !> CALL EITHER THIS ROUTINE OR `*_radius_hvp_eval` FOR A GIVEN ATOM,
   !> NEVER BOTH -- that one fills the same `tg*` accumulators.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    vr_a       Radius direction of atom a
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg1_r      sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg2_rr     sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg3_rrr    sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   pure module subroutine cfc_atomic_radius_tangent_eval(d_a, R_a, a1, fade_start, fade_width, vr_a, &
                                                  max_deriv, tg0, tg1_r, tg2_rr, tg3_rrr)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Radius direction of atom a
      real(wp), intent(in) :: vr_a
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg0
      !> sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
   end subroutine cfc_atomic_radius_tangent_eval

   !> Faded term of `cfc_atomic_radius_tangent_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    vr_a       Radius direction of atom a
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg1_r      sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg2_rr     sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg3_rrr    sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   pure module subroutine cfc_atomic_radius_tangent_eval_faded(d_a, R_a, a1, fade_start, fade_width, &
                                                        vr_a, max_deriv, tg0, tg1_r, tg2_rr, &
                                                        tg3_rrr, sg_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Radius direction of atom a
      real(wp), intent(in) :: vr_a
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg0
      !> sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
   end subroutine cfc_atomic_radius_tangent_eval_faded

   !> Accumulate the *radius half* of the nuclear directional derivative.
   !>
   !> Given one radius direction per atom, this ADDS
   !>
   !>    tg{k}(i1..ik) += sum_A vr_A  d/dR_A  d^k K / dr_i1..dr_ik
   !>
   !> into the same `tg*` accumulators `*_tangent_eval` fills with the position
   !> half.  Together the two make the joint contraction
   !> ``sum_A (v_A . d/dr_A + vr_A d/dR_A)``, which is what `cfc_tangent_eval`
   !> and `cfc_hvp_eval` consume; the lift algebra is identical either way,
   !> because both halves land in the same family.
   !>
   !> Split rather than fused on purpose.  This layer is an O(n**2) pair loop
   !> run at every grid point, and the common caller holds its radii fixed:
   !> fusing would charge that caller for terms it discards on every pair, and
   !> would perturb its results in the last bits by regrouping the CSE block.
   !> A fixed-radius caller therefore simply never calls this routine.
   !>
   !> CALL EITHER THIS ROUTINE OR `*_radius_hvp_eval` FOR A GIVEN PAIR,
   !> NEVER BOTH -- that one fills the same `tg*` accumulators.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    vr_a       Radius direction of atom a
   !> @param[in]    vr_b       Radius direction of atom b
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg1_r      sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg2_rr     sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg3_rrr    sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   pure module subroutine cfc_pair_radius_tangent_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                                fade_width, vr_a, vr_b, max_deriv, tg0, tg1_r, &
                                                tg2_rr, tg3_rrr)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Radius direction of atom a
      real(wp), intent(in) :: vr_a
      !> Radius direction of atom b
      real(wp), intent(in) :: vr_b
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg0
      !> sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
   end subroutine cfc_pair_radius_tangent_eval

   !> Faded term of `cfc_pair_radius_tangent_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    vr_a       Radius direction of atom a
   !> @param[in]    vr_b       Radius direction of atom b
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg1_r      sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg2_rr     sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[inout] tg3_rrr    sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the
   !>                          position half accumulated by the tangent/HVP routines
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_radius_tangent_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                                      fade_width, vr_a, vr_b, max_deriv, tg0, &
                                                      tg1_r, tg2_rr, tg3_rrr, sg_a, sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Radius direction of atom a
      real(wp), intent(in) :: vr_a
      !> Radius direction of atom b
      real(wp), intent(in) :: vr_b
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg0
      !> sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the position half
      !> accumulated by the tangent/HVP routines
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_radius_tangent_eval_faded

   !> Accumulate everything the joint position/radius Hessian adds to `*_hvp_eval`.
   !>
   !> Three ladders, sharing one CSE block because they share their radius
   !> derivatives.  Two are *corrections* -- they add the radius half of a
   !> contraction whose position half `*_hvp_eval` already accumulated:
   !>
   !>    tg{k}(i1..ik)         += sum_B vr_B d/dR_B  d^k K / dr^k
   !>    hvp{k}_A(i1..ik, beta) += sum_B vr_B d^2/(dr_A,beta dR_B)  d^k K / dr^k
   !>
   !> The third is complete in itself, the retained-radius row of the joint
   !> Hessian, and so carries both halves of its own contraction:
   !>
   !>    rh{k}_A(i1..ik) = sum_B (v_B . d/dr_B + vr_B d/dR_B)
   !>                        d/dR_A  d^k K / dr_i1..dr_ik
   !>
   !> A joint-Hessian caller therefore leaves its existing per-ATOM calls
   !> untouched and adds exactly one more.  The `tg*` ladder runs to spatial
   !> order 3, the `hvp*` and `rh*` ladders to 2 (their total order is one
   !> higher), all capped by ``max_deriv``.
   !>
   !> Fusing the three was measured, not assumed: the top pair branch needs 6111
   !> temporaries together against 4427 + 2230 for the corrections and the `rh*`
   !> ladder computed apart.  A modest 8%, but it goes the right way, and every
   !> caller of this routine wants all three.
   !>
   !> CALL EITHER THIS ROUTINE OR `*_radius_tangent_eval` FOR A GIVEN ATOM,
   !> NEVER BOTH -- both fill the `tg*` accumulators.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    vr_a       Radius direction of atom a
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_atomic_hvp_eval
   !> @param[inout] tg1_r      sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_atomic_hvp_eval
   !> @param[inout] tg2_rr     sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_atomic_hvp_eval
   !> @param[inout] tg3_rrr    sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_atomic_hvp_eval
   !> @param[inout] hvp0_a     sum_B vr_B d2/(dr_a dR_B) of the order-0 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] hvp1_r_a   sum_B vr_B d2/(dr_a dR_B) of the order-1 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] hvp2_rr_a  sum_B vr_B d2/(dr_a dR_B) of the order-2 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] rh0_a      sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-0 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh1_r_a    sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-1 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh2_rr_a   sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-2 spatial
   !>                          derivative; complete, not a correction
   pure module subroutine cfc_atomic_radius_hvp_eval(d_a, R_a, a1, fade_start, fade_width, v_a, vr_a, &
                                              max_deriv, tg0, tg1_r, tg2_rr, tg3_rrr, hvp0_a, &
                                              hvp1_r_a, hvp2_rr_a, rh0_a, rh1_r_a, rh2_rr_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Radius direction of atom a
      real(wp), intent(in) :: vr_a
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the position half
      !> accumulated by cfc_atomic_hvp_eval
      real(wp), intent(inout) :: tg0
      !> sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the position half
      !> accumulated by cfc_atomic_hvp_eval
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the position half
      !> accumulated by cfc_atomic_hvp_eval
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the position half
      !> accumulated by cfc_atomic_hvp_eval
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-0 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp0_a(3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-1 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp1_r_a(3, 3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-2 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp2_rr_a(3, 3, 3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-0 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh0_a
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-1 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh1_r_a(3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-2 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh2_rr_a(3, 3)
   end subroutine cfc_atomic_radius_hvp_eval

   !> Faded term of `cfc_atomic_radius_hvp_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    a1         Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    vr_a       Radius direction of atom a
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_atomic_hvp_eval
   !> @param[inout] tg1_r      sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_atomic_hvp_eval
   !> @param[inout] tg2_rr     sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_atomic_hvp_eval
   !> @param[inout] tg3_rrr    sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_atomic_hvp_eval
   !> @param[inout] hvp0_a     sum_B vr_B d2/(dr_a dR_B) of the order-0 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] hvp1_r_a   sum_B vr_B d2/(dr_a dR_B) of the order-1 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] hvp2_rr_a  sum_B vr_B d2/(dr_a dR_B) of the order-2 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] rh0_a      sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-0 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh1_r_a    sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-1 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh2_rr_a   sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-2 spatial
   !>                          derivative; complete, not a correction
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   pure module subroutine cfc_atomic_radius_hvp_eval_faded(d_a, R_a, a1, fade_start, fade_width, v_a, &
                                                    vr_a, max_deriv, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                                    hvp0_a, hvp1_r_a, hvp2_rr_a, rh0_a, rh1_r_a, &
                                                    rh2_rr_a, sg_a)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Atomic-term exponent prefactor a1 (Diedenhofen-Klamt 2018: -15)
      real(wp), intent(in) :: a1
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Radius direction of atom a
      real(wp), intent(in) :: vr_a
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the position half
      !> accumulated by cfc_atomic_hvp_eval
      real(wp), intent(inout) :: tg0
      !> sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the position half
      !> accumulated by cfc_atomic_hvp_eval
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the position half
      !> accumulated by cfc_atomic_hvp_eval
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the position half
      !> accumulated by cfc_atomic_hvp_eval
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-0 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp0_a(3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-1 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp1_r_a(3, 3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-2 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp2_rr_a(3, 3, 3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-0 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh0_a
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-1 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh1_r_a(3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-2 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh2_rr_a(3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
   end subroutine cfc_atomic_radius_hvp_eval_faded

   !> Accumulate everything the joint position/radius Hessian adds to `*_hvp_eval`.
   !>
   !> Three ladders, sharing one CSE block because they share their radius
   !> derivatives.  Two are *corrections* -- they add the radius half of a
   !> contraction whose position half `*_hvp_eval` already accumulated:
   !>
   !>    tg{k}(i1..ik)         += sum_B vr_B d/dR_B  d^k K / dr^k
   !>    hvp{k}_A(i1..ik, beta) += sum_B vr_B d^2/(dr_A,beta dR_B)  d^k K / dr^k
   !>
   !> The third is complete in itself, the retained-radius row of the joint
   !> Hessian, and so carries both halves of its own contraction:
   !>
   !>    rh{k}_A(i1..ik) = sum_B (v_B . d/dr_B + vr_B d/dR_B)
   !>                        d/dR_A  d^k K / dr_i1..dr_ik
   !>
   !> A joint-Hessian caller therefore leaves its existing per-PAIR calls
   !> untouched and adds exactly one more.  The `tg*` ladder runs to spatial
   !> order 3, the `hvp*` and `rh*` ladders to 2 (their total order is one
   !> higher), all capped by ``max_deriv``.
   !>
   !> Fusing the three was measured, not assumed: the top pair branch needs 6111
   !> temporaries together against 4427 + 2230 for the corrections and the `rh*`
   !> ladder computed apart.  A modest 8%, but it goes the right way, and every
   !> caller of this routine wants all three.
   !>
   !> CALL EITHER THIS ROUTINE OR `*_radius_tangent_eval` FOR A GIVEN PAIR,
   !> NEVER BOTH -- both fill the `tg*` accumulators.
   !>
   !> Fade: the term enters as `w_a T_a` (atomic) or `w_a w_b T_ab` (pair),
   !> `w = 1 - S_4(t)`, `t = clamp((sigma - fade_start)/fade_width, 0, 1)`,
   !> `sigma = ||d||/R - 1` per center, the model of `cfc_strip_powersums_o*`
   !>
   !> - plateau (every center at `sigma <= fade_start`): the plain term, the
   !>   unfaded arithmetic; band at infinity: `fade_start = huge(0.0_wp)`
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded term; radius derivatives chain through
   !>   `dsigma/dR = -||d||/R**2`, positions through `dsigma/dr = d/(R ||d||)`
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    v_b        Nuclear displacement direction of atom b
   !> @param[in]    vr_a       Radius direction of atom a
   !> @param[in]    vr_b       Radius direction of atom b
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_pair_hvp_eval
   !> @param[inout] tg1_r      sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_pair_hvp_eval
   !> @param[inout] tg2_rr     sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_pair_hvp_eval
   !> @param[inout] tg3_rrr    sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_pair_hvp_eval
   !> @param[inout] hvp0_a     sum_B vr_B d2/(dr_a dR_B) of the order-0 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] hvp0_b     sum_B vr_B d2/(dr_b dR_B) of the order-0 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom b
   !> @param[inout] hvp1_r_a   sum_B vr_B d2/(dr_a dR_B) of the order-1 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] hvp1_r_b   sum_B vr_B d2/(dr_b dR_B) of the order-1 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom b
   !> @param[inout] hvp2_rr_a  sum_B vr_B d2/(dr_a dR_B) of the order-2 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] hvp2_rr_b  sum_B vr_B d2/(dr_b dR_B) of the order-2 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom b
   !> @param[inout] rh0_a      sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-0 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh0_b      sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-0 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh1_r_a    sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-1 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh1_r_b    sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-1 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh2_rr_a   sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-2 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh2_rr_b   sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-2 spatial
   !>                          derivative; complete, not a correction
   pure module subroutine cfc_pair_radius_hvp_eval(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                            fade_width, v_a, v_b, vr_a, vr_b, max_deriv, tg0, &
                                            tg1_r, tg2_rr, tg3_rrr, hvp0_a, hvp0_b, hvp1_r_a, &
                                            hvp1_r_b, hvp2_rr_a, hvp2_rr_b, rh0_a, rh0_b, &
                                            rh1_r_a, rh1_r_b, rh2_rr_a, rh2_rr_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Nuclear displacement direction of atom b
      real(wp), intent(in) :: v_b(3)
      !> Radius direction of atom a
      real(wp), intent(in) :: vr_a
      !> Radius direction of atom b
      real(wp), intent(in) :: vr_b
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the position half
      !> accumulated by cfc_pair_hvp_eval
      real(wp), intent(inout) :: tg0
      !> sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the position half
      !> accumulated by cfc_pair_hvp_eval
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the position half
      !> accumulated by cfc_pair_hvp_eval
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the position half
      !> accumulated by cfc_pair_hvp_eval
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-0 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp0_a(3)
      !> sum_B vr_B d2/(dr_b dR_B) of the order-0 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom b
      real(wp), intent(inout) :: hvp0_b(3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-1 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp1_r_a(3, 3)
      !> sum_B vr_B d2/(dr_b dR_B) of the order-1 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom b
      real(wp), intent(inout) :: hvp1_r_b(3, 3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-2 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp2_rr_a(3, 3, 3)
      !> sum_B vr_B d2/(dr_b dR_B) of the order-2 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom b
      real(wp), intent(inout) :: hvp2_rr_b(3, 3, 3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-0 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh0_a
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-0 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh0_b
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-1 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh1_r_a(3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-1 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh1_r_b(3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-2 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh2_rr_a(3, 3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-2 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh2_rr_b(3, 3)
   end subroutine cfc_pair_radius_hvp_eval

   !> Faded term of `cfc_pair_radius_hvp_eval`: some center past `fade_start`
   !>
   !> - a center beyond its band (`sigma >= fade_start + fade_width`): adds nothing
   !> - otherwise the faded dispatch; arguments as there, plus sigma of each center
   !>
   !> @param[in]    d_a        Displacement r - r_a
   !> @param[in]    d_b        Displacement r - r_b
   !> @param[in]    R_a        Radius of atom a
   !> @param[in]    R_b        Radius of atom b
   !> @param[in]    a2         Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
   !> @param[in]    c_par      Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
   !> @param[in]    fade_start Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel;
   !>                          huge(0.0_wp) for no fade
   !> @param[in]    fade_width Fade band width in sigma, finite and > 0 (read only when a center is
   !>                          past fade_start)
   !> @param[in]    v_a        Nuclear displacement direction of atom a
   !> @param[in]    v_b        Nuclear displacement direction of atom b
   !> @param[in]    vr_a       Radius direction of atom a
   !> @param[in]    vr_b       Radius direction of atom b
   !> @param[in]    max_deriv  Highest spatial-derivative order to accumulate (0..3); a larger value
   !>                          is treated as 3
   !> @param[inout] tg0        sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_pair_hvp_eval
   !> @param[inout] tg1_r      sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_pair_hvp_eval
   !> @param[inout] tg2_rr     sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_pair_hvp_eval
   !> @param[inout] tg3_rrr    sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the
   !>                          position half accumulated by cfc_pair_hvp_eval
   !> @param[inout] hvp0_a     sum_B vr_B d2/(dr_a dR_B) of the order-0 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] hvp0_b     sum_B vr_B d2/(dr_b dR_B) of the order-0 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom b
   !> @param[inout] hvp1_r_a   sum_B vr_B d2/(dr_a dR_B) of the order-1 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] hvp1_r_b   sum_B vr_B d2/(dr_b dR_B) of the order-1 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom b
   !> @param[inout] hvp2_rr_a  sum_B vr_B d2/(dr_a dR_B) of the order-2 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom a
   !> @param[inout] hvp2_rr_b  sum_B vr_B d2/(dr_b dR_B) of the order-2 spatial derivative, *added*
   !>                          to the position half; trailing index is the retained nuclear
   !>                          component of atom b
   !> @param[inout] rh0_a      sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-0 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh0_b      sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-0 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh1_r_a    sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-1 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh1_r_b    sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-1 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh2_rr_a   sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-2 spatial
   !>                          derivative; complete, not a correction
   !> @param[inout] rh2_rr_b   sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-2 spatial
   !>                          derivative; complete, not a correction
   !> @param[in]    sg_a       Reduced distance sigma = ||d_a||/R_a - 1 of center a
   !> @param[in]    sg_b       Reduced distance sigma = ||d_b||/R_b - 1 of center b
   pure module subroutine cfc_pair_radius_hvp_eval_faded(d_a, d_b, R_a, R_b, a2, c_par, fade_start, &
                                                  fade_width, v_a, v_b, vr_a, vr_b, max_deriv, &
                                                  tg0, tg1_r, tg2_rr, tg3_rrr, hvp0_a, hvp0_b, &
                                                  hvp1_r_a, hvp1_r_b, hvp2_rr_a, hvp2_rr_b, &
                                                  rh0_a, rh0_b, rh1_r_a, rh1_r_b, rh2_rr_a, &
                                                  rh2_rr_b, sg_a, sg_b)
      !> Displacement r - r_a
      real(wp), intent(in) :: d_a(3)
      !> Displacement r - r_b
      real(wp), intent(in) :: d_b(3)
      !> Radius of atom a
      real(wp), intent(in) :: R_a
      !> Radius of atom b
      real(wp), intent(in) :: R_b
      !> Pair-term exponent prefactor a2 (Diedenhofen-Klamt 2018: -9)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c (Diedenhofen-Klamt 2018: 5)
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma = ||d||/R - 1, as passed to the strip kernel; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma, finite and > 0 (read only when a center is past fade_start)
      real(wp), intent(in) :: fade_width
      !> Nuclear displacement direction of atom a
      real(wp), intent(in) :: v_a(3)
      !> Nuclear displacement direction of atom b
      real(wp), intent(in) :: v_b(3)
      !> Radius direction of atom a
      real(wp), intent(in) :: vr_a
      !> Radius direction of atom b
      real(wp), intent(in) :: vr_b
      !> Highest spatial-derivative order to accumulate (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_A vr_A d/dR_A of the order-0 spatial derivative, *added* to the position half
      !> accumulated by cfc_pair_hvp_eval
      real(wp), intent(inout) :: tg0
      !> sum_A vr_A d/dR_A of the order-1 spatial derivative, *added* to the position half
      !> accumulated by cfc_pair_hvp_eval
      real(wp), intent(inout) :: tg1_r(3)
      !> sum_A vr_A d/dR_A of the order-2 spatial derivative, *added* to the position half
      !> accumulated by cfc_pair_hvp_eval
      real(wp), intent(inout) :: tg2_rr(3, 3)
      !> sum_A vr_A d/dR_A of the order-3 spatial derivative, *added* to the position half
      !> accumulated by cfc_pair_hvp_eval
      real(wp), intent(inout) :: tg3_rrr(3, 3, 3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-0 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp0_a(3)
      !> sum_B vr_B d2/(dr_b dR_B) of the order-0 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom b
      real(wp), intent(inout) :: hvp0_b(3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-1 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp1_r_a(3, 3)
      !> sum_B vr_B d2/(dr_b dR_B) of the order-1 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom b
      real(wp), intent(inout) :: hvp1_r_b(3, 3)
      !> sum_B vr_B d2/(dr_a dR_B) of the order-2 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom a
      real(wp), intent(inout) :: hvp2_rr_a(3, 3, 3)
      !> sum_B vr_B d2/(dr_b dR_B) of the order-2 spatial derivative, *added* to the position half;
      !> trailing index is the retained nuclear component of atom b
      real(wp), intent(inout) :: hvp2_rr_b(3, 3, 3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-0 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh0_a
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-0 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh0_b
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-1 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh1_r_a(3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-1 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh1_r_b(3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_a of the order-2 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh2_rr_a(3, 3)
      !> sum_B (v_B . d/dr_B + vr_B d/dR_B) of d/dR_b of the order-2 spatial derivative; complete,
      !> not a correction
      real(wp), intent(inout) :: rh2_rr_b(3, 3)
      !> Reduced distance sigma = ||d_a||/R_a - 1 of center a
      real(wp), intent(in) :: sg_a
      !> Reduced distance sigma = ||d_b||/R_b - 1 of center b
      real(wp), intent(in) :: sg_b
   end subroutine cfc_pair_radius_hvp_eval_faded

   !> Lift the accumulated pseudo-density tensors to the level set.
   !>
   !> The CFC level set is ``S = -log PD``.  PD is a sum over atoms and pairs,
   !> so the log can only be taken once everything has been accumulated; this
   !> routine is that step, and it is O(1) -- call it once per evaluation point
   !> after the `*_term_eval` loops.
   !>
   !> Outputs are derivatives of ``S`` itself, sign included: the caller negates
   !> nothing.  (The previous `cfc_log_lift` returned derivatives of ``+log PD``
   !> and required the caller to negate.  The rename is deliberate, so that a
   !> stale call site fails to compile instead of flipping sign silently.)
   !>
   !> Every result is zeroed on entry, so a fully screened-out point (PD == 0)
   !> and the orders above `max_deriv` both come back as exact zeros; nothing
   !> this routine returns is ever left undefined.
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    pd3_rrr    Accumulated order-3 spatial derivative of PD
   !> @param[in]    pd4_rrrr   Accumulated order-4 spatial derivative of PD
   !> @param[in]    max_deriv  Highest derivative order to lift (0..4); a larger value is treated as
   !>                          4
   !> @param[out]   f0         Order-0 spatial derivative of the level set
   !> @param[out]   f1_r       Order-1 spatial derivative of the level set
   !> @param[out]   f2_rr      Order-2 spatial derivative of the level set
   !> @param[out]   f3_rrr     Order-3 spatial derivative of the level set
   !> @param[out]   f4_rrrr    Order-4 spatial derivative of the level set
   pure module subroutine cfc_spatial_eval(pd0, pd1_r, pd2_rr, pd3_rrr, pd4_rrrr, max_deriv, f0, f1_r, &
                                    f2_rr, f3_rrr, f4_rrrr)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Accumulated order-3 spatial derivative of PD
      real(wp), intent(in) :: pd3_rrr(3, 3, 3)
      !> Accumulated order-4 spatial derivative of PD
      real(wp), intent(in) :: pd4_rrrr(3, 3, 3, 3)
      !> Highest derivative order to lift (0..4); a larger value is treated as 4
      integer, intent(in) :: max_deriv
      !> Order-0 spatial derivative of the level set
      real(wp), intent(out) :: f0
      !> Order-1 spatial derivative of the level set
      real(wp), intent(out) :: f1_r(3)
      !> Order-2 spatial derivative of the level set
      real(wp), intent(out) :: f2_rr(3, 3)
      !> Order-3 spatial derivative of the level set
      real(wp), intent(out) :: f3_rrr(3, 3, 3)
      !> Order-4 spatial derivative of the level set
      real(wp), intent(out) :: f4_rrrr(3, 3, 3, 3)
   end subroutine cfc_spatial_eval

   !> Lift the accumulated tensors to the level-set nuclear gradient.
   !>
   !> With ``S = -log PD`` and one retained nucleus A,
   !>
   !>    f{k+1}_..._rA(i1..ik, beta) = d/dR_A,beta  d^k S / dr_i1..dr_ik
   !>
   !> Two accumulated families feed it: `pd*` (`*_term_eval` /
   !> `cfc_pair_spatial_eval`, an aggregate) and `qn*` (`*_nuclear_eval`, atom
   !> A's own share).  No direction vector is involved, so a caller that wants
   !> only the nuclear gradient never touches the contracted families.
   !>
   !> Call once per atom A per evaluation point, i.e. O(n_active) in total.
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    pd3_rrr    Accumulated order-3 spatial derivative of PD
   !> @param[in]    qn0        Atom-A d/dR_A of the order-0 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qn1_r      Atom-A d/dR_A of the order-1 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qn2_rr     Atom-A d/dR_A of the order-2 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qn3_rrr    Atom-A d/dR_A of the order-3 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    max_deriv  Highest spatial-derivative order to lift (0..3); a larger value is
   !>                          treated as 3
   !> @param[out]   f1_rA      d/dR_A of the order-0 spatial derivative of the level set; trailing
   !>                          index is the nuclear component
   !> @param[out]   f2_r_rA    d/dR_A of the order-1 spatial derivative of the level set; trailing
   !>                          index is the nuclear component
   !> @param[out]   f3_rr_rA   d/dR_A of the order-2 spatial derivative of the level set; trailing
   !>                          index is the nuclear component
   !> @param[out]   f4_rrr_rA  d/dR_A of the order-3 spatial derivative of the level set; trailing
   !>                          index is the nuclear component
   pure module subroutine cfc_nuclear_eval(pd0, pd1_r, pd2_rr, pd3_rrr, qn0, qn1_r, qn2_rr, qn3_rrr, &
                                    max_deriv, f1_rA, f2_r_rA, f3_rr_rA, f4_rrr_rA)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Accumulated order-3 spatial derivative of PD
      real(wp), intent(in) :: pd3_rrr(3, 3, 3)
      !> Atom-A d/dR_A of the order-0 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qn0(3)
      !> Atom-A d/dR_A of the order-1 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qn1_r(3, 3)
      !> Atom-A d/dR_A of the order-2 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qn2_rr(3, 3, 3)
      !> Atom-A d/dR_A of the order-3 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qn3_rrr(3, 3, 3, 3)
      !> Highest spatial-derivative order to lift (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> d/dR_A of the order-0 spatial derivative of the level set; trailing index is the nuclear
      !> component
      real(wp), intent(out) :: f1_rA(3)
      !> d/dR_A of the order-1 spatial derivative of the level set; trailing index is the nuclear
      !> component
      real(wp), intent(out) :: f2_r_rA(3, 3)
      !> d/dR_A of the order-2 spatial derivative of the level set; trailing index is the nuclear
      !> component
      real(wp), intent(out) :: f3_rr_rA(3, 3, 3)
      !> d/dR_A of the order-3 spatial derivative of the level set; trailing index is the nuclear
      !> component
      real(wp), intent(out) :: f4_rrr_rA(3, 3, 3, 3)
   end subroutine cfc_nuclear_eval

   !> Lift the accumulated tensors to the uncontracted level-set nuclear Hessian.
   !>
   !> With ``S = -log PD`` and two retained nuclei A and B,
   !>
   !>    f{k+2}_..._rArB(i1..ik, beta, gamma) = d2/(dR_A,beta dR_B,gamma)
   !>                                           d^k S / dr_i1..dr_ik
   !>
   !> Four accumulated families feed it.  `pd*` is the aggregate; `qa*` and
   !> `qb*` are the per-atom nuclear derivatives of A and of B (both from
   !> `*_nuclear_eval`, just evaluated for different atoms); `qq*` is the
   !> two-nucleus family of `*_hessian_eval` for this particular pair.
   !>
   !> Call once per ordered atom pair per evaluation point, i.e. O(n_active**2)
   !> in total -- which is the size of the result, not an avoidable cost.  A
   !> caller that only needs the Hessian contracted against a displacement
   !> field should use `cfc_hvp_eval` instead and stay O(n_active).
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    qa0        Atom-A d/dR_A of the order-0 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qa1_r      Atom-A d/dR_A of the order-1 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qa2_rr     Atom-A d/dR_A of the order-2 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qb0        Atom-B d/dR_B of the order-0 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qb1_r      Atom-B d/dR_B of the order-1 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qb2_rr     Atom-B d/dR_B of the order-2 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qq0        d2/(dR_A dR_B) of the order-0 spatial derivative of PD; the two
   !>                          trailing indices are the nuclear components of A and B
   !> @param[in]    qq1_r      d2/(dR_A dR_B) of the order-1 spatial derivative of PD; the two
   !>                          trailing indices are the nuclear components of A and B
   !> @param[in]    qq2_rr     d2/(dR_A dR_B) of the order-2 spatial derivative of PD; the two
   !>                          trailing indices are the nuclear components of A and B
   !> @param[in]    max_deriv  Highest spatial-derivative order to lift (0..2); a larger value is
   !>                          treated as 2
   !> @param[out]   f2_rArB    d2/(dR_A dR_B) of the order-0 spatial derivative of the level set;
   !>                          the two trailing indices are the nuclear components of A and B
   !> @param[out]   f3_r_rArB  d2/(dR_A dR_B) of the order-1 spatial derivative of the level set;
   !>                          the two trailing indices are the nuclear components of A and B
   !> @param[out]   f4_rr_rArB d2/(dR_A dR_B) of the order-2 spatial derivative of the level set;
   !>                          the two trailing indices are the nuclear components of A and B
   pure module subroutine cfc_hessian_eval(pd0, pd1_r, pd2_rr, qa0, qa1_r, qa2_rr, qb0, qb1_r, qb2_rr, &
                                    qq0, qq1_r, qq2_rr, max_deriv, f2_rArB, f3_r_rArB, f4_rr_rArB)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Atom-A d/dR_A of the order-0 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qa0(3)
      !> Atom-A d/dR_A of the order-1 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qa1_r(3, 3)
      !> Atom-A d/dR_A of the order-2 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qa2_rr(3, 3, 3)
      !> Atom-B d/dR_B of the order-0 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qb0(3)
      !> Atom-B d/dR_B of the order-1 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qb1_r(3, 3)
      !> Atom-B d/dR_B of the order-2 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qb2_rr(3, 3, 3)
      !> d2/(dR_A dR_B) of the order-0 spatial derivative of PD; the two trailing indices are the
      !> nuclear components of A and B
      real(wp), intent(in) :: qq0(3, 3)
      !> d2/(dR_A dR_B) of the order-1 spatial derivative of PD; the two trailing indices are the
      !> nuclear components of A and B
      real(wp), intent(in) :: qq1_r(3, 3, 3)
      !> d2/(dR_A dR_B) of the order-2 spatial derivative of PD; the two trailing indices are the
      !> nuclear components of A and B
      real(wp), intent(in) :: qq2_rr(3, 3, 3, 3)
      !> Highest spatial-derivative order to lift (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dR_A dR_B) of the order-0 spatial derivative of the level set; the two trailing indices
      !> are the nuclear components of A and B
      real(wp), intent(out) :: f2_rArB(3, 3)
      !> d2/(dR_A dR_B) of the order-1 spatial derivative of the level set; the two trailing indices
      !> are the nuclear components of A and B
      real(wp), intent(out) :: f3_r_rArB(3, 3, 3)
      !> d2/(dR_A dR_B) of the order-2 spatial derivative of the level set; the two trailing indices
      !> are the nuclear components of A and B
      real(wp), intent(out) :: f4_rr_rArB(3, 3, 3, 3)
   end subroutine cfc_hessian_eval

   !> Lift the accumulated tensors to the level-set directional derivative.
   !>
   !> With ``S = -log PD`` and one displacement direction per atom,
   !>
   !>    tangent_f{k}(i1..ik) = sum_B sum_beta v_B,beta
   !>                           d/dR_B,beta  d^k S / dr_i1..dr_ik
   !>
   !> Both inputs are aggregates: `pd*` from the `*_term_eval` loops and `tg*`
   !> from the `*_tangent_eval` loops, each summed over every atom and pair.
   !> The lift itself is O(1) -- one call per evaluation point.
   !>
   !> This is the SvdW `svdw_tangent_eval` contract, term for term: the caller
   !> receives derivatives of the level set and never applies a chain rule.
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    pd3_rrr    Accumulated order-3 spatial derivative of PD
   !> @param[in]    tg0        Accumulated sum_B v_B . d/dR_B of the order-0 derivative of PD
   !> @param[in]    tg1_r      Accumulated sum_B v_B . d/dR_B of the order-1 derivative of PD
   !> @param[in]    tg2_rr     Accumulated sum_B v_B . d/dR_B of the order-2 derivative of PD
   !> @param[in]    tg3_rrr    Accumulated sum_B v_B . d/dR_B of the order-3 derivative of PD
   !> @param[in]    max_deriv  Highest spatial-derivative order to lift (0..3); a larger value is
   !>                          treated as 3
   !> @param[out]   tangent_f0 sum_B v_B . d/dR_B of the order-0 spatial derivative of the level set
   !> @param[out]   tangent_f1_r sum_B v_B . d/dR_B of the order-1 spatial derivative of the level
   !>                            set
   !> @param[out]   tangent_f2_rr sum_B v_B . d/dR_B of the order-2 spatial derivative of the level
   !>                             set
   !> @param[out]   tangent_f3_rrr sum_B v_B . d/dR_B of the order-3 spatial derivative of the level
   !>                              set
   pure module subroutine cfc_tangent_eval(pd0, pd1_r, pd2_rr, pd3_rrr, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                    max_deriv, tangent_f0, tangent_f1_r, tangent_f2_rr, &
                                    tangent_f3_rrr)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Accumulated order-3 spatial derivative of PD
      real(wp), intent(in) :: pd3_rrr(3, 3, 3)
      !> Accumulated sum_B v_B . d/dR_B of the order-0 derivative of PD
      real(wp), intent(in) :: tg0
      !> Accumulated sum_B v_B . d/dR_B of the order-1 derivative of PD
      real(wp), intent(in) :: tg1_r(3)
      !> Accumulated sum_B v_B . d/dR_B of the order-2 derivative of PD
      real(wp), intent(in) :: tg2_rr(3, 3)
      !> Accumulated sum_B v_B . d/dR_B of the order-3 derivative of PD
      real(wp), intent(in) :: tg3_rrr(3, 3, 3)
      !> Highest spatial-derivative order to lift (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> sum_B v_B . d/dR_B of the order-0 spatial derivative of the level set
      real(wp), intent(out) :: tangent_f0
      !> sum_B v_B . d/dR_B of the order-1 spatial derivative of the level set
      real(wp), intent(out) :: tangent_f1_r(3)
      !> sum_B v_B . d/dR_B of the order-2 spatial derivative of the level set
      real(wp), intent(out) :: tangent_f2_rr(3, 3)
      !> sum_B v_B . d/dR_B of the order-3 spatial derivative of the level set
      real(wp), intent(out) :: tangent_f3_rrr(3, 3, 3)
   end subroutine cfc_tangent_eval

   !> Lift the accumulated tensors to the level-set Hessian-vector product.
   !>
   !> With ``S = -log PD``, one displacement direction per atom and one retained
   !> nucleus A,
   !>
   !>    hvp_f{k+1}_..._rA(i1..ik, beta) = sum_B sum_gamma v_B,gamma
   !>                                      d^2/(dR_A,beta dR_B,gamma)
   !>                                      d^k S / dr_i1..dr_ik
   !>
   !> Four accumulated families feed it.  Two are aggregates, summed over all
   !> atoms and pairs: `pd*` (`*_term_eval`) and `tg*` (`*_tangent_eval`).  Two
   !> are per-atom, summed over atom A's own term and every pair containing A:
   !> `qn*` (`*_nuclear_eval`) and `hvp*` (`*_hvp_eval`).  Nothing here is
   !> O(n**2) in storage; the pairwise nuclear Hessian is never formed.
   !>
   !> Call once per atom A per evaluation point, i.e. O(n_active) in total.
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    tg0        Accumulated sum_B v_B . d/dR_B of the order-0 derivative of PD
   !> @param[in]    tg1_r      Accumulated sum_B v_B . d/dR_B of the order-1 derivative of PD
   !> @param[in]    tg2_rr     Accumulated sum_B v_B . d/dR_B of the order-2 derivative of PD
   !> @param[in]    qn0        Atom-A d/dR_A of the order-0 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qn1_r      Atom-A d/dR_A of the order-1 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qn2_rr     Atom-A d/dR_A of the order-2 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    hvp0       Atom-A sum_B v_B . d2/(dR_A dR_B) of the order-0 spatial derivative
   !>                          of PD; trailing index is the retained nuclear component
   !> @param[in]    hvp1_r     Atom-A sum_B v_B . d2/(dR_A dR_B) of the order-1 spatial derivative
   !>                          of PD; trailing index is the retained nuclear component
   !> @param[in]    hvp2_rr    Atom-A sum_B v_B . d2/(dR_A dR_B) of the order-2 spatial derivative
   !>                          of PD; trailing index is the retained nuclear component
   !> @param[in]    max_deriv  Highest spatial-derivative order to lift (0..2); a larger value is
   !>                          treated as 2
   !> @param[out]   hvp_f1_rA  sum_B v_B . d2/(dR_A dR_B) of the order-0 spatial derivative of the
   !>                          level set; trailing index is the retained nuclear component
   !> @param[out]   hvp_f2_r_rA sum_B v_B . d2/(dR_A dR_B) of the order-1 spatial derivative of the
   !>                           level set; trailing index is the retained nuclear component
   !> @param[out]   hvp_f3_rr_rA sum_B v_B . d2/(dR_A dR_B) of the order-2 spatial derivative of the
   !>                            level set; trailing index is the retained nuclear component
   pure module subroutine cfc_hvp_eval(pd0, pd1_r, pd2_rr, tg0, tg1_r, tg2_rr, qn0, qn1_r, qn2_rr, hvp0, &
                                hvp1_r, hvp2_rr, max_deriv, hvp_f1_rA, hvp_f2_r_rA, hvp_f3_rr_rA)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Accumulated sum_B v_B . d/dR_B of the order-0 derivative of PD
      real(wp), intent(in) :: tg0
      !> Accumulated sum_B v_B . d/dR_B of the order-1 derivative of PD
      real(wp), intent(in) :: tg1_r(3)
      !> Accumulated sum_B v_B . d/dR_B of the order-2 derivative of PD
      real(wp), intent(in) :: tg2_rr(3, 3)
      !> Atom-A d/dR_A of the order-0 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qn0(3)
      !> Atom-A d/dR_A of the order-1 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qn1_r(3, 3)
      !> Atom-A d/dR_A of the order-2 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qn2_rr(3, 3, 3)
      !> Atom-A sum_B v_B . d2/(dR_A dR_B) of the order-0 spatial derivative of PD; trailing index
      !> is the retained nuclear component
      real(wp), intent(in) :: hvp0(3)
      !> Atom-A sum_B v_B . d2/(dR_A dR_B) of the order-1 spatial derivative of PD; trailing index
      !> is the retained nuclear component
      real(wp), intent(in) :: hvp1_r(3, 3)
      !> Atom-A sum_B v_B . d2/(dR_A dR_B) of the order-2 spatial derivative of PD; trailing index
      !> is the retained nuclear component
      real(wp), intent(in) :: hvp2_rr(3, 3, 3)
      !> Highest spatial-derivative order to lift (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> sum_B v_B . d2/(dR_A dR_B) of the order-0 spatial derivative of the level set; trailing
      !> index is the retained nuclear component
      real(wp), intent(out) :: hvp_f1_rA(3)
      !> sum_B v_B . d2/(dR_A dR_B) of the order-1 spatial derivative of the level set; trailing
      !> index is the retained nuclear component
      real(wp), intent(out) :: hvp_f2_r_rA(3, 3)
      !> sum_B v_B . d2/(dR_A dR_B) of the order-2 spatial derivative of the level set; trailing
      !> index is the retained nuclear component
      real(wp), intent(out) :: hvp_f3_rr_rA(3, 3, 3)
   end subroutine cfc_hvp_eval

   !> Lift the accumulated tensors to the adjoint-contracted nuclear gradient.
   !>
   !> The reverse-mode mirror of `cfc_hvp_eval`.  With ``S = -log PD``, one
   !> retained nucleus A and a per-point adjoint jet ``(w0, w1, w2)``,
   !>
   !>    vjp_f1_rA(beta) = w0 f1_rA(beta)
   !>                    + sum_a w1(a) f2_r_rA(a, beta)
   !>                    + sum_a sum_b w2(a, b) f3_rr_rA(a, b, beta)
   !>
   !> i.e. the `f1_rA` rung with the jet indices contracted away, exactly as
   !> `hvp_f1_rA` is that same rung contracted with a nuclear direction.
   !>
   !> Same two accumulated families as `cfc_nuclear_eval` and the same orders:
   !> `pd*` (`*_term_eval` / `cfc_pair_spatial_eval`, an aggregate) and `qn*`
   !> (`*_nuclear_eval`, atom A's own share), both up to order 2.  A caller that
   !> has prepared enough for `f3_rr_rA` has prepared enough for this.
   !>
   !> The contraction happens here rather than at the call site, and that is the
   !> point: the weights enter before the common-subexpression pass, so this
   !> routine evaluates the tree that produces three numbers instead of the one
   !> that produces 3 + 9 + 27, and no caller ever forms the rank-3 nuclear
   !> tensor.  `w2` is a general 3x3 -- every one of the nine entries is read,
   !> no symmetry is assumed and no factor of two is folded in.
   !>
   !> There is no `max_deriv` argument: the jet spans a closed set of orders, so
   !> there is no lower-order subset to dispatch on.
   !>
   !> Call once per atom A per evaluation point, i.e. O(n_active) in total.
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    qn0        Atom-A d/dR_A of the order-0 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qn1_r      Atom-A d/dR_A of the order-1 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qn2_rr     Atom-A d/dR_A of the order-2 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    w0         Adjoint weight of the value rung `f1_rA`
   !> @param[in]    w1         Adjoint weights of the gradient rung `f2_r_rA`
   !> @param[in]    w2         Adjoint weights of the Hessian rung `f3_rr_rA`; a general 3x3, read
   !>                          entry by entry
   !> @param[out]   vjp_f1_rA  Adjoint jet contracted onto the nuclear gradient of the level set;
   !>                          the only index is the retained nuclear component
   pure module subroutine cfc_vjp_eval(pd0, pd1_r, pd2_rr, qn0, qn1_r, qn2_rr, w0, w1, w2, vjp_f1_rA)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Atom-A d/dR_A of the order-0 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qn0(3)
      !> Atom-A d/dR_A of the order-1 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qn1_r(3, 3)
      !> Atom-A d/dR_A of the order-2 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qn2_rr(3, 3, 3)
      !> Adjoint weight of the value rung `f1_rA`
      real(wp), intent(in) :: w0
      !> Adjoint weights of the gradient rung `f2_r_rA`
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the Hessian rung `f3_rr_rA`; a general 3x3, read entry by entry
      real(wp), intent(in) :: w2(3, 3)
      !> Adjoint jet contracted onto the nuclear gradient of the level set; the only index is the
      !> retained nuclear component
      real(wp), intent(out) :: vjp_f1_rA(3)
   end subroutine cfc_vjp_eval

   !> Lift the accumulated tensors to the level-set radius gradient.
   !>
   !> With ``S = -log PD`` and one retained atom A,
   !>
   !>    f{k+1}_..._rad(i1..ik) = d/dRad_A  d^k S / dr_i1..dr_ik
   !>
   !> where ``Rad_A`` is atom A's cavity radius.  The radius sibling of
   !> `cfc_nuclear_eval`, and structurally the same lift; the outputs simply do
   !> not grow a nuclear index, because a radius is a scalar.
   !>
   !> Two accumulated families feed it: `pd*` (`*_term_eval` /
   !> `cfc_pair_spatial_eval`, an aggregate) and `rd*` (`*_radius_eval`, atom
   !> A's own share).  No direction vector is involved.
   !>
   !> Spelled `_rad` and not `_rA` because Fortran identifiers are
   !> case-insensitive: `f3_rr_rA` and `f3_rr_Ra` would be one name.
   !>
   !> Call once per atom A per evaluation point, i.e. O(n_active) in total.
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    pd3_rrr    Accumulated order-3 spatial derivative of PD
   !> @param[in]    rd0        Atom-A d/dRad_A of the order-0 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rd1_r      Atom-A d/dRad_A of the order-1 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rd2_rr     Atom-A d/dRad_A of the order-2 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rd3_rrr    Atom-A d/dRad_A of the order-3 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    max_deriv  Highest spatial-derivative order to lift (0..3); a larger value is
   !>                          treated as 3
   !> @param[out]   f1_rad     d/dRad_A of the order-0 spatial derivative of the level set
   !> @param[out]   f2_r_rad   d/dRad_A of the order-1 spatial derivative of the level set
   !> @param[out]   f3_rr_rad  d/dRad_A of the order-2 spatial derivative of the level set
   !> @param[out]   f4_rrr_rad d/dRad_A of the order-3 spatial derivative of the level set
   pure module subroutine cfc_radius_eval(pd0, pd1_r, pd2_rr, pd3_rrr, rd0, rd1_r, rd2_rr, rd3_rrr, &
                                   max_deriv, f1_rad, f2_r_rad, f3_rr_rad, f4_rrr_rad)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Accumulated order-3 spatial derivative of PD
      real(wp), intent(in) :: pd3_rrr(3, 3, 3)
      !> Atom-A d/dRad_A of the order-0 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rd0
      !> Atom-A d/dRad_A of the order-1 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rd1_r(3)
      !> Atom-A d/dRad_A of the order-2 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rd2_rr(3, 3)
      !> Atom-A d/dRad_A of the order-3 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rd3_rrr(3, 3, 3)
      !> Highest spatial-derivative order to lift (0..3); a larger value is treated as 3
      integer, intent(in) :: max_deriv
      !> d/dRad_A of the order-0 spatial derivative of the level set
      real(wp), intent(out) :: f1_rad
      !> d/dRad_A of the order-1 spatial derivative of the level set
      real(wp), intent(out) :: f2_r_rad(3)
      !> d/dRad_A of the order-2 spatial derivative of the level set
      real(wp), intent(out) :: f3_rr_rad(3, 3)
      !> d/dRad_A of the order-3 spatial derivative of the level set
      real(wp), intent(out) :: f4_rrr_rad(3, 3, 3)
   end subroutine cfc_radius_eval

   !> Lift the accumulated tensors to the uncontracted level-set radius Hessian.
   !>
   !> With ``S = -log PD`` and two retained radii,
   !>
   !>    f{k+2}_..._radrad(i1..ik) = d2/(dRad_A dRad_B)  d^k S / dr_i1..dr_ik
   !>
   !> Four accumulated families feed it, mirroring `cfc_hessian_eval` slot for
   !> slot: `pd*` is the aggregate; `ra*` and `rb*` are the per-atom radius
   !> derivatives of A and of B (both from `*_radius_eval`, just evaluated for
   !> different atoms); `rr*` is the two-radius family of `*_radius2_eval` for
   !> this particular pair.
   !>
   !> Call once per atom pair per evaluation point.  A caller that only needs the
   !> radius block contracted against a direction should use `cfc_radius_hvp_eval`
   !> instead and stay O(n_active).
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    ra0        Atom-A d/dRad_A of the order-0 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    ra1_r      Atom-A d/dRad_A of the order-1 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    ra2_rr     Atom-A d/dRad_A of the order-2 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rb0        Atom-B d/dRad_B of the order-0 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rb1_r      Atom-B d/dRad_B of the order-1 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rb2_rr     Atom-B d/dRad_B of the order-2 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rr0        d2/(dRad_A dRad_B) of the order-0 spatial derivative of PD
   !> @param[in]    rr1_r      d2/(dRad_A dRad_B) of the order-1 spatial derivative of PD
   !> @param[in]    rr2_rr     d2/(dRad_A dRad_B) of the order-2 spatial derivative of PD
   !> @param[in]    max_deriv  Highest spatial-derivative order to lift (0..2); a larger value is
   !>                          treated as 2
   !> @param[out]   f2_radrad  d2/(dRad_A dRad_B) of the order-0 spatial derivative of the level set
   !> @param[out]   f3_r_radrad d2/(dRad_A dRad_B) of the order-1 spatial derivative of the level
   !>                           set
   !> @param[out]   f4_rr_radrad d2/(dRad_A dRad_B) of the order-2 spatial derivative of the level
   !>                            set
   pure module subroutine cfc_radius2_eval(pd0, pd1_r, pd2_rr, ra0, ra1_r, ra2_rr, rb0, rb1_r, rb2_rr, &
                                    rr0, rr1_r, rr2_rr, max_deriv, f2_radrad, f3_r_radrad, &
                                    f4_rr_radrad)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Atom-A d/dRad_A of the order-0 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: ra0
      !> Atom-A d/dRad_A of the order-1 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: ra1_r(3)
      !> Atom-A d/dRad_A of the order-2 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: ra2_rr(3, 3)
      !> Atom-B d/dRad_B of the order-0 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rb0
      !> Atom-B d/dRad_B of the order-1 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rb1_r(3)
      !> Atom-B d/dRad_B of the order-2 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rb2_rr(3, 3)
      !> d2/(dRad_A dRad_B) of the order-0 spatial derivative of PD
      real(wp), intent(in) :: rr0
      !> d2/(dRad_A dRad_B) of the order-1 spatial derivative of PD
      real(wp), intent(in) :: rr1_r(3)
      !> d2/(dRad_A dRad_B) of the order-2 spatial derivative of PD
      real(wp), intent(in) :: rr2_rr(3, 3)
      !> Highest spatial-derivative order to lift (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dRad_A dRad_B) of the order-0 spatial derivative of the level set
      real(wp), intent(out) :: f2_radrad
      !> d2/(dRad_A dRad_B) of the order-1 spatial derivative of the level set
      real(wp), intent(out) :: f3_r_radrad(3)
      !> d2/(dRad_A dRad_B) of the order-2 spatial derivative of the level set
      real(wp), intent(out) :: f4_rr_radrad(3, 3)
   end subroutine cfc_radius2_eval

   !> Lift the accumulated tensors to the uncontracted mixed nuclear-radius block.
   !>
   !> With ``S = -log PD``, a retained nucleus A and a retained radius of B,
   !>
   !>    f{k+2}_..._rA_rad(i1..ik, beta) = d2/(dR_A,beta dRad_B)
   !>                                      d^k S / dr_i1..dr_ik
   !>
   !> Four accumulated families feed it: the aggregate `pd*`, atom A's nuclear
   !> derivative `qa*` (from `*_nuclear_eval`), atom B's radius derivative `rb*`
   !> (from `*_radius_eval`), and `mx*` of `*_nucrad_eval` -- the two-letter stem
   !> `mx` names the nuclear x radius family throughout, alongside `pd`/`tg` for
   !> the aggregates, `qn`/`qa`/`qb`/`qq` for retained nuclei and
   !> `rd`/`ra`/`rb`/`rr` for retained radii.
   !>
   !> The result is not symmetric in the two atoms; call the routine again with
   !> the roles of A and B swapped to obtain the other half of the block.
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    qa0        Atom-A d/dR_A of the order-0 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qa1_r      Atom-A d/dR_A of the order-1 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    qa2_rr     Atom-A d/dR_A of the order-2 spatial derivative of PD; trailing index
   !>                          is the nuclear component
   !> @param[in]    rb0        Atom-B d/dRad_B of the order-0 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rb1_r      Atom-B d/dRad_B of the order-1 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rb2_rr     Atom-B d/dRad_B of the order-2 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    mx0        d2/(dR_A dRad_B) of the order-0 spatial derivative of PD; trailing
   !>                          index is atom A's nuclear component
   !> @param[in]    mx1_r      d2/(dR_A dRad_B) of the order-1 spatial derivative of PD; trailing
   !>                          index is atom A's nuclear component
   !> @param[in]    mx2_rr     d2/(dR_A dRad_B) of the order-2 spatial derivative of PD; trailing
   !>                          index is atom A's nuclear component
   !> @param[in]    max_deriv  Highest spatial-derivative order to lift (0..2); a larger value is
   !>                          treated as 2
   !> @param[out]   f2_rA_rad  d2/(dR_A dRad_B) of the order-0 spatial derivative of the level set;
   !>                          trailing index is atom A's nuclear component
   !> @param[out]   f3_r_rA_rad d2/(dR_A dRad_B) of the order-1 spatial derivative of the level set;
   !>                           trailing index is atom A's nuclear component
   !> @param[out]   f4_rr_rA_rad d2/(dR_A dRad_B) of the order-2 spatial derivative of the level
   !>                            set; trailing index is atom A's nuclear component
   pure module subroutine cfc_nucrad_eval(pd0, pd1_r, pd2_rr, qa0, qa1_r, qa2_rr, rb0, rb1_r, rb2_rr, &
                                   mx0, mx1_r, mx2_rr, max_deriv, f2_rA_rad, f3_r_rA_rad, &
                                   f4_rr_rA_rad)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Atom-A d/dR_A of the order-0 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qa0(3)
      !> Atom-A d/dR_A of the order-1 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qa1_r(3, 3)
      !> Atom-A d/dR_A of the order-2 spatial derivative of PD; trailing index is the nuclear
      !> component
      real(wp), intent(in) :: qa2_rr(3, 3, 3)
      !> Atom-B d/dRad_B of the order-0 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rb0
      !> Atom-B d/dRad_B of the order-1 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rb1_r(3)
      !> Atom-B d/dRad_B of the order-2 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rb2_rr(3, 3)
      !> d2/(dR_A dRad_B) of the order-0 spatial derivative of PD; trailing index is atom A's
      !> nuclear component
      real(wp), intent(in) :: mx0(3)
      !> d2/(dR_A dRad_B) of the order-1 spatial derivative of PD; trailing index is atom A's
      !> nuclear component
      real(wp), intent(in) :: mx1_r(3, 3)
      !> d2/(dR_A dRad_B) of the order-2 spatial derivative of PD; trailing index is atom A's
      !> nuclear component
      real(wp), intent(in) :: mx2_rr(3, 3, 3)
      !> Highest spatial-derivative order to lift (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> d2/(dR_A dRad_B) of the order-0 spatial derivative of the level set; trailing index is atom
      !> A's nuclear component
      real(wp), intent(out) :: f2_rA_rad(3)
      !> d2/(dR_A dRad_B) of the order-1 spatial derivative of the level set; trailing index is atom
      !> A's nuclear component
      real(wp), intent(out) :: f3_r_rA_rad(3, 3)
      !> d2/(dR_A dRad_B) of the order-2 spatial derivative of the level set; trailing index is atom
      !> A's nuclear component
      real(wp), intent(out) :: f4_rr_rA_rad(3, 3, 3)
   end subroutine cfc_nucrad_eval

   !> Lift the accumulated tensors to the radius row of the joint Hessian.
   !>
   !> With ``S = -log PD``, one retained atom A and a direction that moves both
   !> the nuclei and the radii,
   !>
   !>    hvp_f{k+1}_..._rad(i1..ik) = sum_B (sum_gamma v_B,gamma d/dR_B,gamma
   !>                                        + vr_B d/dRad_B)
   !>                                 d/dRad_A  d^k S / dr_i1..dr_ik
   !>
   !> `cfc_hvp_eval` returns the position row of the same Hessian; between them
   !> the two cover all four blocks (nuclear-nuclear, nuclear-radius,
   !> radius-nuclear, radius-radius), because each contracted index runs over
   !> positions *and* radii.
   !>
   !> Four accumulated families feed it.  Two are aggregates: `pd*`
   !> (`*_term_eval`) and `tg*` (`*_tangent_eval` plus the radius half from
   !> `*_radius_tangent_eval` / `*_radius_hvp_eval`).  Two are per-atom, summed
   !> over atom A's own term and every pair containing A: `rd*`
   !> (`*_radius_eval`) and `rh*` (`*_radius_hvp_eval`).
   !>
   !> Call once per atom A per evaluation point, i.e. O(n_active) in total.
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    tg0        Accumulated sum_B (v_B . d/dR_B + vr_B d/dRad_B) of the order-0
   !>                          derivative of PD
   !> @param[in]    tg1_r      Accumulated sum_B (v_B . d/dR_B + vr_B d/dRad_B) of the order-1
   !>                          derivative of PD
   !> @param[in]    tg2_rr     Accumulated sum_B (v_B . d/dR_B + vr_B d/dRad_B) of the order-2
   !>                          derivative of PD
   !> @param[in]    rd0        Atom-A d/dRad_A of the order-0 spatial derivative of PD
   !> @param[in]    rd1_r      Atom-A d/dRad_A of the order-1 spatial derivative of PD
   !> @param[in]    rd2_rr     Atom-A d/dRad_A of the order-2 spatial derivative of PD
   !> @param[in]    rh0        Atom-A sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the
   !>                          order-0 spatial derivative of PD
   !> @param[in]    rh1_r      Atom-A sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the
   !>                          order-1 spatial derivative of PD
   !> @param[in]    rh2_rr     Atom-A sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the
   !>                          order-2 spatial derivative of PD
   !> @param[in]    max_deriv  Highest spatial-derivative order to lift (0..2); a larger value is
   !>                          treated as 2
   !> @param[out]   hvp_f1_rad sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the order-0
   !>                          spatial derivative of the level set
   !> @param[out]   hvp_f2_r_rad sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the order-1
   !>                            spatial derivative of the level set
   !> @param[out]   hvp_f3_rr_rad sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the order-2
   !>                             spatial derivative of the level set
   pure module subroutine cfc_radius_hvp_eval(pd0, pd1_r, pd2_rr, tg0, tg1_r, tg2_rr, rd0, rd1_r, &
                                       rd2_rr, rh0, rh1_r, rh2_rr, max_deriv, hvp_f1_rad, &
                                       hvp_f2_r_rad, hvp_f3_rr_rad)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Accumulated sum_B (v_B . d/dR_B + vr_B d/dRad_B) of the order-0 derivative of PD
      real(wp), intent(in) :: tg0
      !> Accumulated sum_B (v_B . d/dR_B + vr_B d/dRad_B) of the order-1 derivative of PD
      real(wp), intent(in) :: tg1_r(3)
      !> Accumulated sum_B (v_B . d/dR_B + vr_B d/dRad_B) of the order-2 derivative of PD
      real(wp), intent(in) :: tg2_rr(3, 3)
      !> Atom-A d/dRad_A of the order-0 spatial derivative of PD
      real(wp), intent(in) :: rd0
      !> Atom-A d/dRad_A of the order-1 spatial derivative of PD
      real(wp), intent(in) :: rd1_r(3)
      !> Atom-A d/dRad_A of the order-2 spatial derivative of PD
      real(wp), intent(in) :: rd2_rr(3, 3)
      !> Atom-A sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the order-0 spatial derivative
      !> of PD
      real(wp), intent(in) :: rh0
      !> Atom-A sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the order-1 spatial derivative
      !> of PD
      real(wp), intent(in) :: rh1_r(3)
      !> Atom-A sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the order-2 spatial derivative
      !> of PD
      real(wp), intent(in) :: rh2_rr(3, 3)
      !> Highest spatial-derivative order to lift (0..2); a larger value is treated as 2
      integer, intent(in) :: max_deriv
      !> sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the order-0 spatial derivative of the
      !> level set
      real(wp), intent(out) :: hvp_f1_rad
      !> sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the order-1 spatial derivative of the
      !> level set
      real(wp), intent(out) :: hvp_f2_r_rad(3)
      !> sum_B (v_B . d/dR_B + vr_B d/dRad_B) of d/dRad_A of the order-2 spatial derivative of the
      !> level set
      real(wp), intent(out) :: hvp_f3_rr_rad(3, 3)
   end subroutine cfc_radius_hvp_eval

   !> Lift the accumulated tensors to the adjoint-contracted radius gradient.
   !>
   !> The radius twin of `cfc_vjp_eval`.  With ``S = -log PD``, one retained
   !> atom A and a per-point adjoint jet ``(w0, w1, w2)``,
   !>
   !>    vjp_f1_rad = w0 f1_rad
   !>               + sum_a w1(a) f2_r_rad(a)
   !>               + sum_a sum_b w2(a, b) f3_rr_rad(a, b)
   !>
   !> i.e. the `f1_rad` rung with the jet indices contracted away.  A radius is
   !> a scalar parameter, so the result is a single number per atom, where the
   !> nuclear twin returns one Cartesian component -- the same difference
   !> `cfc_radius_eval` has from `cfc_nuclear_eval`.
   !>
   !> Same two accumulated families as `cfc_radius_eval` and the same orders:
   !> `pd*` (`*_term_eval` / `cfc_pair_spatial_eval`, an aggregate) and `rd*`
   !> (`*_radius_eval`, atom A's own share), both up to order 2.  A caller that
   !> has prepared enough for `f3_rr_rad` has prepared enough for this; no
   !> direction vector is involved.
   !>
   !> The contraction happens here rather than at the call site, and that is the
   !> point: the weights enter before the common-subexpression pass, so this
   !> routine evaluates the tree that produces one number instead of the one
   !> that produces 1 + 3 + 9, and no caller ever forms the rank-2 radius
   !> tensor.  `w2` is a general 3x3 -- every one of the nine entries is read,
   !> no symmetry is assumed and no factor of two is folded in.
   !>
   !> There is no `max_deriv` argument: the jet spans a closed set of orders, so
   !> there is no lower-order subset to dispatch on.
   !>
   !> Call once per atom A per evaluation point, i.e. O(n_active) in total.
   !>
   !> @param[in]    pd0        Accumulated order-0 spatial derivative of PD
   !> @param[in]    pd1_r      Accumulated order-1 spatial derivative of PD
   !> @param[in]    pd2_rr     Accumulated order-2 spatial derivative of PD
   !> @param[in]    rd0        Atom-A d/dRad_A of the order-0 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rd1_r      Atom-A d/dRad_A of the order-1 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    rd2_rr     Atom-A d/dRad_A of the order-2 spatial derivative of PD (no trailing
   !>                          index: a radius is a scalar parameter)
   !> @param[in]    w0         Adjoint weight of the value rung `f1_rad`
   !> @param[in]    w1         Adjoint weights of the gradient rung `f2_r_rad`
   !> @param[in]    w2         Adjoint weights of the Hessian rung `f3_rr_rad`; a general 3x3, read
   !>                          entry by entry
   !> @param[out]   vjp_f1_rad Adjoint jet contracted onto the radius gradient of the level set; a
   !>                          scalar, because a radius carries no index
   pure module subroutine cfc_radius_vjp_eval(pd0, pd1_r, pd2_rr, rd0, rd1_r, rd2_rr, w0, w1, w2, &
                                       vjp_f1_rad)
      !> Accumulated order-0 spatial derivative of PD
      real(wp), intent(in) :: pd0
      !> Accumulated order-1 spatial derivative of PD
      real(wp), intent(in) :: pd1_r(3)
      !> Accumulated order-2 spatial derivative of PD
      real(wp), intent(in) :: pd2_rr(3, 3)
      !> Atom-A d/dRad_A of the order-0 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rd0
      !> Atom-A d/dRad_A of the order-1 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rd1_r(3)
      !> Atom-A d/dRad_A of the order-2 spatial derivative of PD (no trailing index: a radius is a
      !> scalar parameter)
      real(wp), intent(in) :: rd2_rr(3, 3)
      !> Adjoint weight of the value rung `f1_rad`
      real(wp), intent(in) :: w0
      !> Adjoint weights of the gradient rung `f2_r_rad`
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the Hessian rung `f3_rr_rad`; a general 3x3, read entry by entry
      real(wp), intent(in) :: w2(3, 3)
      !> Adjoint jet contracted onto the radius gradient of the level set; a scalar, because a
      !> radius carries no index
      real(wp), intent(out) :: vjp_f1_rad
   end subroutine cfc_radius_vjp_eval

   !> Strip pseudo-density jet through order 0, point loop innermost
   !>
   !> - accumulates `d^k/dr^k PD(r_p)`, k = 0..0, for the
   !>   `np <= cfc_strip_width(0) = 32` points of one strip over one shared
   !>   atom list: atomic terms in list order, then pairs `A < B` in list order;
   !>   only orbit representatives (packed slots)
   !> - per-atom factors in the reduced distance `sigma_A = ||r - R_A||/R_a - 1`:
   !>   `u_A = exp(a2 sigma_A)` (pair envelope, the screening factor) and
   !>   `e_A = exp(a1 sigma_A)` (atomic term), both from the inline exp
   !> - fade `w_A = 1 - S_4(t)`, `t = clamp((sigma_A - fade_start)/fade_width, 0, 1)`:
   !>   atomic term `w_A e_A`, pair term `f(q_AB) (w_A u_A)(w_B u_B)`,
   !>   `f(q) = c (1 - q)**4`, `q_AB = (r - R_A).(r - R_B)/(R_a R_b)`
   !> - per atom against the strip box: outside the gate sphere or beyond the band
   !>   (`sigma_min >= fade_start + fade_width`) is skipped, plateau
   !>   (`sigma_max <= fade_start`) runs the plain jets, else the faded jets; the
   !>   pair pass sees only survivors and skips pairs whose reach spheres are
   !>   disjoint; band at infinity: `fade_start = huge(0.0_wp)` with a finite
   !>   `fade_width`
   !> - pair body: Leibniz product of the `f(q)` jet (Faa di Bruno, `q` is
   !>   quadratic in r) and both atoms' cached Cartesian jets of `w u`; plateau
   !>   and faded pairs share it
   !> - each CSE temporary emitted just before its first reader
   !> - hard screening per point: `u = e = 0` where `||r - R_A||**2 > bound_sq(A)`,
   !>   the gate's own squared compare
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not
   !>   handled: their sums are not finite (`1/x` is Inf, no FPE trap assumed) and
   !>   the caller re-evaluates the strip on the per-point path
   !> - requires `a1 <= 0`, `a2 <= 0`, `max(-a1, -a2) <= 708`
   !> - geometry cache `geo(:, s, a)`, rows `1:np`, written for every list atom
   !>   not skipped (rows of skipped atoms are not written): s = 1-3 `r - R_A`,
   !>   4 `1/||r - R_A||`, 5..5 the gate-masked Cartesian jet
   !>   `d^J (w_A u_A)/dr^J`, |J| = 0..0, packed, by order; the caller's buffer
   !>   holds `cfc_strip_geo_len` doubles per atom for every order
   !> - survivor and reach workspace `sv(nat)`, `reach(nat)` from the caller, not
   !>   read on entry: the routine owns no `nat`-sized automatic array
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    points       Strip points, point-last
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's `cand_screen` row
   !> @param[in]    a1           Atomic-term exponent a1 (<= 0)
   !> @param[in]    a2           Pair-term exponent a2 (<= 0)
   !> @param[in]    c_par        Pair-term coupling constant c
   !> @param[in]    fade_start   Fade onset in sigma_A = ||r - R_A||/R_a - 1; huge(0.0_wp) for no
   !>                            fade
   !> @param[in]    fade_width   Fade band width in sigma_A, finite and > 0
   !> @param[out]   pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[out]   n_zero       Point-atom pairs at zero distance among the evaluated atoms;
   !>                            nonzero means the sums are not finite
   !> @param[out]   geo          Point-first geometry cache and pair-pass workspace, one block per
   !>                            list atom (slots in the routine notes)
   !> @param[out]   sv           Survivor workspace: list indices of the atoms kept by the strip
   !>                            classification, leading entries written
   !> @param[out]   reach        Reach workspace: pair-test radius of each survivor, leading entries
   !>                            written
   pure module subroutine cfc_strip_powersums_o0(np, points, nat, centers, radii, bound_sq, a1, a2, &
                                          c_par, fade_start, fade_width, pd0, n_zero, geo, sv, &
                                          reach)
      !> Strip width of this order
      integer, parameter :: nps = cfc_strip_width(0)
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
      !> Atomic-term exponent a1 (<= 0)
      real(wp), intent(in) :: a1
      !> Pair-term exponent a2 (<= 0)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma_A = ||r - R_A||/R_a - 1; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma_A, finite and > 0
      real(wp), intent(in) :: fade_width
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(out) :: pd0(np)
      !> Point-atom pairs at zero distance among the evaluated atoms; nonzero means the sums are not
      !> finite
      integer, intent(out) :: n_zero
      !> Point-first geometry cache and pair-pass workspace, one block per list atom (slots in the
      !> routine notes)
      real(wp), intent(out) :: geo(nps, 5, nat)
      !> Survivor workspace: list indices of the atoms kept by the strip classification, leading
      !> entries written
      integer, intent(out) :: sv(nat)
      !> Reach workspace: pair-test radius of each survivor, leading entries written
      real(wp), intent(out) :: reach(nat)
   end subroutine cfc_strip_powersums_o0

   !> Strip lift through order 0: packed spatial jet of S = -log PD
   !>
   !> - same algebra as `cfc_spatial_eval` at `max_deriv = 0`; writes `lsf0` point-last, symmetric
   !>   groups packed (sorted tuples, last index slowest)
   !> - `np <= cfc_strip_width(0) = 32`; requires `pd0 > 0` for every point (escape is the caller's
   !>   test on the stage-1 sums), no zero fill
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[out]   lsf0         Level-set value S = -log(PD)
   pure module subroutine cfc_strip_lift_o0(np, pd0, lsf0)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(in) :: pd0(np)
      !> Level-set value S = -log(PD)
      real(wp), intent(out) :: lsf0(np)
   end subroutine cfc_strip_lift_o0

   !> Strip pseudo-density jet through order 1, point loop innermost
   !>
   !> - accumulates `d^k/dr^k PD(r_p)`, k = 0..1, for the
   !>   `np <= cfc_strip_width(1) = 32` points of one strip over one shared
   !>   atom list: atomic terms in list order, then pairs `A < B` in list order;
   !>   only orbit representatives (packed slots)
   !> - per-atom factors in the reduced distance `sigma_A = ||r - R_A||/R_a - 1`:
   !>   `u_A = exp(a2 sigma_A)` (pair envelope, the screening factor) and
   !>   `e_A = exp(a1 sigma_A)` (atomic term), both from the inline exp
   !> - fade `w_A = 1 - S_4(t)`, `t = clamp((sigma_A - fade_start)/fade_width, 0, 1)`:
   !>   atomic term `w_A e_A`, pair term `f(q_AB) (w_A u_A)(w_B u_B)`,
   !>   `f(q) = c (1 - q)**4`, `q_AB = (r - R_A).(r - R_B)/(R_a R_b)`
   !> - per atom against the strip box: outside the gate sphere or beyond the band
   !>   (`sigma_min >= fade_start + fade_width`) is skipped, plateau
   !>   (`sigma_max <= fade_start`) runs the plain jets, else the faded jets; the
   !>   pair pass sees only survivors and skips pairs whose reach spheres are
   !>   disjoint; band at infinity: `fade_start = huge(0.0_wp)` with a finite
   !>   `fade_width`
   !> - pair body: Leibniz product of the `f(q)` jet (Faa di Bruno, `q` is
   !>   quadratic in r) and both atoms' cached Cartesian jets of `w u`; plateau
   !>   and faded pairs share it
   !> - each CSE temporary emitted just before its first reader
   !> - hard screening per point: `u = e = 0` where `||r - R_A||**2 > bound_sq(A)`,
   !>   the gate's own squared compare
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not
   !>   handled: their sums are not finite (`1/x` is Inf, no FPE trap assumed) and
   !>   the caller re-evaluates the strip on the per-point path
   !> - requires `a1 <= 0`, `a2 <= 0`, `max(-a1, -a2) <= 708`
   !> - geometry cache `geo(:, s, a)`, rows `1:np`, written for every list atom
   !>   not skipped (rows of skipped atoms are not written): s = 1-3 `r - R_A`,
   !>   4 `1/||r - R_A||`, 5..8 the gate-masked Cartesian jet
   !>   `d^J (w_A u_A)/dr^J`, |J| = 0..1, packed, by order; the caller's buffer
   !>   holds `cfc_strip_geo_len` doubles per atom for every order
   !> - survivor and reach workspace `sv(nat)`, `reach(nat)` from the caller, not
   !>   read on entry: the routine owns no `nat`-sized automatic array
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    points       Strip points, point-last
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's `cand_screen` row
   !> @param[in]    a1           Atomic-term exponent a1 (<= 0)
   !> @param[in]    a2           Pair-term exponent a2 (<= 0)
   !> @param[in]    c_par        Pair-term coupling constant c
   !> @param[in]    fade_start   Fade onset in sigma_A = ||r - R_A||/R_a - 1; huge(0.0_wp) for no
   !>                            fade
   !> @param[in]    fade_width   Fade band width in sigma_A, finite and > 0
   !> @param[out]   pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[out]   pd1          Order-1 pseudo-density jet of the strip, point-first; second index
   !>                            the Cartesian component
   !> @param[out]   n_zero       Point-atom pairs at zero distance among the evaluated atoms;
   !>                            nonzero means the sums are not finite
   !> @param[out]   geo          Point-first geometry cache and pair-pass workspace, one block per
   !>                            list atom (slots in the routine notes)
   !> @param[out]   sv           Survivor workspace: list indices of the atoms kept by the strip
   !>                            classification, leading entries written
   !> @param[out]   reach        Reach workspace: pair-test radius of each survivor, leading entries
   !>                            written
   pure module subroutine cfc_strip_powersums_o1(np, points, nat, centers, radii, bound_sq, a1, a2, &
                                          c_par, fade_start, fade_width, pd0, pd1, n_zero, geo, &
                                          sv, reach)
      !> Strip width of this order
      integer, parameter :: nps = cfc_strip_width(1)
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
      !> Atomic-term exponent a1 (<= 0)
      real(wp), intent(in) :: a1
      !> Pair-term exponent a2 (<= 0)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma_A = ||r - R_A||/R_a - 1; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma_A, finite and > 0
      real(wp), intent(in) :: fade_width
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(out) :: pd0(np)
      !> Order-1 pseudo-density jet of the strip, point-first; second index the Cartesian component
      real(wp), intent(out) :: pd1(np, 3)
      !> Point-atom pairs at zero distance among the evaluated atoms; nonzero means the sums are not
      !> finite
      integer, intent(out) :: n_zero
      !> Point-first geometry cache and pair-pass workspace, one block per list atom (slots in the
      !> routine notes)
      real(wp), intent(out) :: geo(nps, 8, nat)
      !> Survivor workspace: list indices of the atoms kept by the strip classification, leading
      !> entries written
      integer, intent(out) :: sv(nat)
      !> Reach workspace: pair-test radius of each survivor, leading entries written
      real(wp), intent(out) :: reach(nat)
   end subroutine cfc_strip_powersums_o1

   !> Strip lift through order 1: packed spatial jet of S = -log PD
   !>
   !> - same algebra as `cfc_spatial_eval` at `max_deriv = 1`; writes `lsf0`, `lsf1_r` point-last,
   !>   symmetric groups packed (sorted tuples, last index slowest)
   !> - `np <= cfc_strip_width(1) = 32`; requires `pd0 > 0` for every point (escape is the caller's
   !>   test on the stage-1 sums), no zero fill
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[in]    pd1          Order-1 pseudo-density jet of the strip, point-first; second index
   !>                            the Cartesian component
   !> @param[out]   lsf0         Level-set value S = -log(PD)
   !> @param[out]   lsf1_r       Spatial gradient dS/dr
   pure module subroutine cfc_strip_lift_o1(np, pd0, pd1, lsf0, lsf1_r)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(in) :: pd0(np)
      !> Order-1 pseudo-density jet of the strip, point-first; second index the Cartesian component
      real(wp), intent(in) :: pd1(np, 3)
      !> Level-set value S = -log(PD)
      real(wp), intent(out) :: lsf0(np)
      !> Spatial gradient dS/dr
      real(wp), intent(out) :: lsf1_r(3, np)
   end subroutine cfc_strip_lift_o1

   !> Strip pseudo-density jet through order 2, point loop innermost
   !>
   !> - accumulates `d^k/dr^k PD(r_p)`, k = 0..2, for the
   !>   `np <= cfc_strip_width(2) = 32` points of one strip over one shared
   !>   atom list: atomic terms in list order, then pairs `A < B` in list order;
   !>   only orbit representatives (packed slots)
   !> - per-atom factors in the reduced distance `sigma_A = ||r - R_A||/R_a - 1`:
   !>   `u_A = exp(a2 sigma_A)` (pair envelope, the screening factor) and
   !>   `e_A = exp(a1 sigma_A)` (atomic term), both from the inline exp
   !> - fade `w_A = 1 - S_4(t)`, `t = clamp((sigma_A - fade_start)/fade_width, 0, 1)`:
   !>   atomic term `w_A e_A`, pair term `f(q_AB) (w_A u_A)(w_B u_B)`,
   !>   `f(q) = c (1 - q)**4`, `q_AB = (r - R_A).(r - R_B)/(R_a R_b)`
   !> - per atom against the strip box: outside the gate sphere or beyond the band
   !>   (`sigma_min >= fade_start + fade_width`) is skipped, plateau
   !>   (`sigma_max <= fade_start`) runs the plain jets, else the faded jets; the
   !>   pair pass sees only survivors and skips pairs whose reach spheres are
   !>   disjoint; band at infinity: `fade_start = huge(0.0_wp)` with a finite
   !>   `fade_width`
   !> - pair body: Leibniz product of the `f(q)` jet (Faa di Bruno, `q` is
   !>   quadratic in r) and both atoms' cached Cartesian jets of `w u`; plateau
   !>   and faded pairs share it
   !> - each CSE temporary emitted just before its first reader
   !> - hard screening per point: `u = e = 0` where `||r - R_A||**2 > bound_sq(A)`,
   !>   the gate's own squared compare
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not
   !>   handled: their sums are not finite (`1/x` is Inf, no FPE trap assumed) and
   !>   the caller re-evaluates the strip on the per-point path
   !> - requires `a1 <= 0`, `a2 <= 0`, `max(-a1, -a2) <= 708`
   !> - geometry cache `geo(:, s, a)`, rows `1:np`, written for every list atom
   !>   not skipped (rows of skipped atoms are not written): s = 1-3 `r - R_A`,
   !>   4 `1/||r - R_A||`, 5..14 the gate-masked Cartesian jet
   !>   `d^J (w_A u_A)/dr^J`, |J| = 0..2, packed, by order; the caller's buffer
   !>   holds `cfc_strip_geo_len` doubles per atom for every order
   !> - survivor and reach workspace `sv(nat)`, `reach(nat)` from the caller, not
   !>   read on entry: the routine owns no `nat`-sized automatic array
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    points       Strip points, point-last
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's `cand_screen` row
   !> @param[in]    a1           Atomic-term exponent a1 (<= 0)
   !> @param[in]    a2           Pair-term exponent a2 (<= 0)
   !> @param[in]    c_par        Pair-term coupling constant c
   !> @param[in]    fade_start   Fade onset in sigma_A = ||r - R_A||/R_a - 1; huge(0.0_wp) for no
   !>                            fade
   !> @param[in]    fade_width   Fade band width in sigma_A, finite and > 0
   !> @param[out]   pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[out]   pd1          Order-1 pseudo-density jet of the strip, point-first; second index
   !>                            the Cartesian component
   !> @param[out]   pd2          Order-2 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-2 tuple
   !> @param[out]   n_zero       Point-atom pairs at zero distance among the evaluated atoms;
   !>                            nonzero means the sums are not finite
   !> @param[out]   geo          Point-first geometry cache and pair-pass workspace, one block per
   !>                            list atom (slots in the routine notes)
   !> @param[out]   sv           Survivor workspace: list indices of the atoms kept by the strip
   !>                            classification, leading entries written
   !> @param[out]   reach        Reach workspace: pair-test radius of each survivor, leading entries
   !>                            written
   pure module subroutine cfc_strip_powersums_o2(np, points, nat, centers, radii, bound_sq, a1, a2, &
                                          c_par, fade_start, fade_width, pd0, pd1, pd2, n_zero, &
                                          geo, sv, reach)
      !> Strip width of this order
      integer, parameter :: nps = cfc_strip_width(2)
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
      !> Atomic-term exponent a1 (<= 0)
      real(wp), intent(in) :: a1
      !> Pair-term exponent a2 (<= 0)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma_A = ||r - R_A||/R_a - 1; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma_A, finite and > 0
      real(wp), intent(in) :: fade_width
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(out) :: pd0(np)
      !> Order-1 pseudo-density jet of the strip, point-first; second index the Cartesian component
      real(wp), intent(out) :: pd1(np, 3)
      !> Order-2 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-2 tuple
      real(wp), intent(out) :: pd2(np, 6)
      !> Point-atom pairs at zero distance among the evaluated atoms; nonzero means the sums are not
      !> finite
      integer, intent(out) :: n_zero
      !> Point-first geometry cache and pair-pass workspace, one block per list atom (slots in the
      !> routine notes)
      real(wp), intent(out) :: geo(nps, 14, nat)
      !> Survivor workspace: list indices of the atoms kept by the strip classification, leading
      !> entries written
      integer, intent(out) :: sv(nat)
      !> Reach workspace: pair-test radius of each survivor, leading entries written
      real(wp), intent(out) :: reach(nat)
   end subroutine cfc_strip_powersums_o2

   !> Strip lift through order 2: packed spatial jet of S = -log PD
   !>
   !> - same algebra as `cfc_spatial_eval` at `max_deriv = 2`; writes `lsf0`, `lsf1_r`, `lsf2_rr`
   !>   point-last, symmetric groups packed (sorted tuples, last index slowest)
   !> - `np <= cfc_strip_width(2) = 32`; requires `pd0 > 0` for every point (escape is the caller's
   !>   test on the stage-1 sums), no zero fill
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[in]    pd1          Order-1 pseudo-density jet of the strip, point-first; second index
   !>                            the Cartesian component
   !> @param[in]    pd2          Order-2 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-2 tuple
   !> @param[out]   lsf0         Level-set value S = -log(PD)
   !> @param[out]   lsf1_r       Spatial gradient dS/dr
   !> @param[out]   lsf2_rr      Spatial Hessian, packed (6)
   pure module subroutine cfc_strip_lift_o2(np, pd0, pd1, pd2, lsf0, lsf1_r, lsf2_rr)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(in) :: pd0(np)
      !> Order-1 pseudo-density jet of the strip, point-first; second index the Cartesian component
      real(wp), intent(in) :: pd1(np, 3)
      !> Order-2 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-2 tuple
      real(wp), intent(in) :: pd2(np, 6)
      !> Level-set value S = -log(PD)
      real(wp), intent(out) :: lsf0(np)
      !> Spatial gradient dS/dr
      real(wp), intent(out) :: lsf1_r(3, np)
      !> Spatial Hessian, packed (6)
      real(wp), intent(out) :: lsf2_rr(6, np)
   end subroutine cfc_strip_lift_o2

   !> Strip pseudo-density jet through order 3, point loop innermost
   !>
   !> - accumulates `d^k/dr^k PD(r_p)`, k = 0..3, for the
   !>   `np <= cfc_strip_width(3) = 16` points of one strip over one shared
   !>   atom list: atomic terms in list order, then pairs `A < B` in list order;
   !>   only orbit representatives (packed slots)
   !> - per-atom factors in the reduced distance `sigma_A = ||r - R_A||/R_a - 1`:
   !>   `u_A = exp(a2 sigma_A)` (pair envelope, the screening factor) and
   !>   `e_A = exp(a1 sigma_A)` (atomic term), both from the inline exp
   !> - fade `w_A = 1 - S_4(t)`, `t = clamp((sigma_A - fade_start)/fade_width, 0, 1)`:
   !>   atomic term `w_A e_A`, pair term `f(q_AB) (w_A u_A)(w_B u_B)`,
   !>   `f(q) = c (1 - q)**4`, `q_AB = (r - R_A).(r - R_B)/(R_a R_b)`
   !> - per atom against the strip box: outside the gate sphere or beyond the band
   !>   (`sigma_min >= fade_start + fade_width`) is skipped, plateau
   !>   (`sigma_max <= fade_start`) runs the plain jets, else the faded jets; the
   !>   pair pass sees only survivors and skips pairs whose reach spheres are
   !>   disjoint; band at infinity: `fade_start = huge(0.0_wp)` with a finite
   !>   `fade_width`
   !> - pair body: Leibniz product of the `f(q)` jet (Faa di Bruno, `q` is
   !>   quadratic in r) and both atoms' cached Cartesian jets of `w u`; plateau
   !>   and faded pairs share it
   !> - each CSE temporary emitted just before its first reader
   !> - hard screening per point: `u = e = 0` where `||r - R_A||**2 > bound_sq(A)`,
   !>   the gate's own squared compare
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not
   !>   handled: their sums are not finite (`1/x` is Inf, no FPE trap assumed) and
   !>   the caller re-evaluates the strip on the per-point path
   !> - requires `a1 <= 0`, `a2 <= 0`, `max(-a1, -a2) <= 708`
   !> - geometry cache `geo(:, s, a)`, rows `1:np`, written for every list atom
   !>   not skipped (rows of skipped atoms are not written): s = 1-3 `r - R_A`,
   !>   4 `1/||r - R_A||`, 5..24 the gate-masked Cartesian jet
   !>   `d^J (w_A u_A)/dr^J`, |J| = 0..3, packed, by order; the caller's buffer
   !>   holds `cfc_strip_geo_len` doubles per atom for every order
   !> - survivor and reach workspace `sv(nat)`, `reach(nat)` from the caller, not
   !>   read on entry: the routine owns no `nat`-sized automatic array
   !>
   !> @param[in]    np           Live points in the strip, 1..16
   !> @param[in]    points       Strip points, point-last
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's `cand_screen` row
   !> @param[in]    a1           Atomic-term exponent a1 (<= 0)
   !> @param[in]    a2           Pair-term exponent a2 (<= 0)
   !> @param[in]    c_par        Pair-term coupling constant c
   !> @param[in]    fade_start   Fade onset in sigma_A = ||r - R_A||/R_a - 1; huge(0.0_wp) for no
   !>                            fade
   !> @param[in]    fade_width   Fade band width in sigma_A, finite and > 0
   !> @param[out]   pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[out]   pd1          Order-1 pseudo-density jet of the strip, point-first; second index
   !>                            the Cartesian component
   !> @param[out]   pd2          Order-2 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-2 tuple
   !> @param[out]   pd3          Order-3 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-3 tuple
   !> @param[out]   n_zero       Point-atom pairs at zero distance among the evaluated atoms;
   !>                            nonzero means the sums are not finite
   !> @param[out]   geo          Point-first geometry cache and pair-pass workspace, one block per
   !>                            list atom (slots in the routine notes)
   !> @param[out]   sv           Survivor workspace: list indices of the atoms kept by the strip
   !>                            classification, leading entries written
   !> @param[out]   reach        Reach workspace: pair-test radius of each survivor, leading entries
   !>                            written
   pure module subroutine cfc_strip_powersums_o3(np, points, nat, centers, radii, bound_sq, a1, a2, &
                                          c_par, fade_start, fade_width, pd0, pd1, pd2, pd3, &
                                          n_zero, geo, sv, reach)
      !> Strip width of this order
      integer, parameter :: nps = cfc_strip_width(3)
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
      !> Atomic-term exponent a1 (<= 0)
      real(wp), intent(in) :: a1
      !> Pair-term exponent a2 (<= 0)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma_A = ||r - R_A||/R_a - 1; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma_A, finite and > 0
      real(wp), intent(in) :: fade_width
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(out) :: pd0(np)
      !> Order-1 pseudo-density jet of the strip, point-first; second index the Cartesian component
      real(wp), intent(out) :: pd1(np, 3)
      !> Order-2 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-2 tuple
      real(wp), intent(out) :: pd2(np, 6)
      !> Order-3 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-3 tuple
      real(wp), intent(out) :: pd3(np, 10)
      !> Point-atom pairs at zero distance among the evaluated atoms; nonzero means the sums are not
      !> finite
      integer, intent(out) :: n_zero
      !> Point-first geometry cache and pair-pass workspace, one block per list atom (slots in the
      !> routine notes)
      real(wp), intent(out) :: geo(nps, 24, nat)
      !> Survivor workspace: list indices of the atoms kept by the strip classification, leading
      !> entries written
      integer, intent(out) :: sv(nat)
      !> Reach workspace: pair-test radius of each survivor, leading entries written
      real(wp), intent(out) :: reach(nat)
   end subroutine cfc_strip_powersums_o3

   !> Strip lift through order 3: packed spatial jet of S = -log PD
   !>
   !> - same algebra as `cfc_spatial_eval` at `max_deriv = 3`; writes `lsf0`, `lsf1_r`, `lsf2_rr`,
   !>   `lsf3_rrr` point-last, symmetric groups packed (sorted tuples, last index slowest)
   !> - `np <= cfc_strip_width(3) = 16`; requires `pd0 > 0` for every point (escape is the caller's
   !>   test on the stage-1 sums), no zero fill
   !>
   !> @param[in]    np           Live points in the strip, 1..16
   !> @param[in]    pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[in]    pd1          Order-1 pseudo-density jet of the strip, point-first; second index
   !>                            the Cartesian component
   !> @param[in]    pd2          Order-2 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-2 tuple
   !> @param[in]    pd3          Order-3 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-3 tuple
   !> @param[out]   lsf0         Level-set value S = -log(PD)
   !> @param[out]   lsf1_r       Spatial gradient dS/dr
   !> @param[out]   lsf2_rr      Spatial Hessian, packed (6)
   !> @param[out]   lsf3_rrr     Third spatial derivative, packed (10)
   pure module subroutine cfc_strip_lift_o3(np, pd0, pd1, pd2, pd3, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      !> Live points in the strip, 1..16
      integer, intent(in) :: np
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(in) :: pd0(np)
      !> Order-1 pseudo-density jet of the strip, point-first; second index the Cartesian component
      real(wp), intent(in) :: pd1(np, 3)
      !> Order-2 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-2 tuple
      real(wp), intent(in) :: pd2(np, 6)
      !> Order-3 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-3 tuple
      real(wp), intent(in) :: pd3(np, 10)
      !> Level-set value S = -log(PD)
      real(wp), intent(out) :: lsf0(np)
      !> Spatial gradient dS/dr
      real(wp), intent(out) :: lsf1_r(3, np)
      !> Spatial Hessian, packed (6)
      real(wp), intent(out) :: lsf2_rr(6, np)
      !> Third spatial derivative, packed (10)
      real(wp), intent(out) :: lsf3_rrr(10, np)
   end subroutine cfc_strip_lift_o3

   !> Strip pseudo-density jet through order 4, point loop innermost
   !>
   !> - accumulates `d^k/dr^k PD(r_p)`, k = 0..4, for the
   !>   `np <= cfc_strip_width(4) = 8` points of one strip over one shared
   !>   atom list: atomic terms in list order, then pairs `A < B` in list order;
   !>   only orbit representatives (packed slots)
   !> - per-atom factors in the reduced distance `sigma_A = ||r - R_A||/R_a - 1`:
   !>   `u_A = exp(a2 sigma_A)` (pair envelope, the screening factor) and
   !>   `e_A = exp(a1 sigma_A)` (atomic term), both from the inline exp
   !> - fade `w_A = 1 - S_4(t)`, `t = clamp((sigma_A - fade_start)/fade_width, 0, 1)`:
   !>   atomic term `w_A e_A`, pair term `f(q_AB) (w_A u_A)(w_B u_B)`,
   !>   `f(q) = c (1 - q)**4`, `q_AB = (r - R_A).(r - R_B)/(R_a R_b)`
   !> - per atom against the strip box: outside the gate sphere or beyond the band
   !>   (`sigma_min >= fade_start + fade_width`) is skipped, plateau
   !>   (`sigma_max <= fade_start`) runs the plain jets, else the faded jets; the
   !>   pair pass sees only survivors and skips pairs whose reach spheres are
   !>   disjoint; band at infinity: `fade_start = huge(0.0_wp)` with a finite
   !>   `fade_width`
   !> - pair body: Leibniz product of the `f(q)` jet (Faa di Bruno, `q` is
   !>   quadratic in r) and both atoms' cached Cartesian jets of `w u`; plateau
   !>   and faded pairs share it
   !> - each CSE temporary emitted just before its first reader
   !> - hard screening per point: `u = e = 0` where `||r - R_A||**2 > bound_sq(A)`,
   !>   the gate's own squared compare
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not
   !>   handled: their sums are not finite (`1/x` is Inf, no FPE trap assumed) and
   !>   the caller re-evaluates the strip on the per-point path
   !> - requires `a1 <= 0`, `a2 <= 0`, `max(-a1, -a2) <= 708`
   !> - geometry cache `geo(:, s, a)`, rows `1:np`, written for every list atom
   !>   not skipped (rows of skipped atoms are not written): s = 1-3 `r - R_A`,
   !>   4 `1/||r - R_A||`, 5..39 the gate-masked Cartesian jet
   !>   `d^J (w_A u_A)/dr^J`, |J| = 0..4, packed, by order; the caller's buffer
   !>   holds `cfc_strip_geo_len` doubles per atom for every order
   !> - survivor and reach workspace `sv(nat)`, `reach(nat)` from the caller, not
   !>   read on entry: the routine owns no `nat`-sized automatic array
   !>
   !> @param[in]    np           Live points in the strip, 1..8
   !> @param[in]    points       Strip points, point-last
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's `cand_screen` row
   !> @param[in]    a1           Atomic-term exponent a1 (<= 0)
   !> @param[in]    a2           Pair-term exponent a2 (<= 0)
   !> @param[in]    c_par        Pair-term coupling constant c
   !> @param[in]    fade_start   Fade onset in sigma_A = ||r - R_A||/R_a - 1; huge(0.0_wp) for no
   !>                            fade
   !> @param[in]    fade_width   Fade band width in sigma_A, finite and > 0
   !> @param[out]   pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[out]   pd1          Order-1 pseudo-density jet of the strip, point-first; second index
   !>                            the Cartesian component
   !> @param[out]   pd2          Order-2 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-2 tuple
   !> @param[out]   pd3          Order-3 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-3 tuple
   !> @param[out]   pd4          Order-4 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-4 tuple
   !> @param[out]   n_zero       Point-atom pairs at zero distance among the evaluated atoms;
   !>                            nonzero means the sums are not finite
   !> @param[out]   geo          Point-first geometry cache and pair-pass workspace, one block per
   !>                            list atom (slots in the routine notes)
   !> @param[out]   sv           Survivor workspace: list indices of the atoms kept by the strip
   !>                            classification, leading entries written
   !> @param[out]   reach        Reach workspace: pair-test radius of each survivor, leading entries
   !>                            written
   pure module subroutine cfc_strip_powersums_o4(np, points, nat, centers, radii, bound_sq, a1, a2, &
                                          c_par, fade_start, fade_width, pd0, pd1, pd2, pd3, pd4, &
                                          n_zero, geo, sv, reach)
      !> Strip width of this order
      integer, parameter :: nps = cfc_strip_width(4)
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
      !> Atomic-term exponent a1 (<= 0)
      real(wp), intent(in) :: a1
      !> Pair-term exponent a2 (<= 0)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma_A = ||r - R_A||/R_a - 1; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma_A, finite and > 0
      real(wp), intent(in) :: fade_width
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(out) :: pd0(np)
      !> Order-1 pseudo-density jet of the strip, point-first; second index the Cartesian component
      real(wp), intent(out) :: pd1(np, 3)
      !> Order-2 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-2 tuple
      real(wp), intent(out) :: pd2(np, 6)
      !> Order-3 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-3 tuple
      real(wp), intent(out) :: pd3(np, 10)
      !> Order-4 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-4 tuple
      real(wp), intent(out) :: pd4(np, 15)
      !> Point-atom pairs at zero distance among the evaluated atoms; nonzero means the sums are not
      !> finite
      integer, intent(out) :: n_zero
      !> Point-first geometry cache and pair-pass workspace, one block per list atom (slots in the
      !> routine notes)
      real(wp), intent(out) :: geo(nps, 39, nat)
      !> Survivor workspace: list indices of the atoms kept by the strip classification, leading
      !> entries written
      integer, intent(out) :: sv(nat)
      !> Reach workspace: pair-test radius of each survivor, leading entries written
      real(wp), intent(out) :: reach(nat)
   end subroutine cfc_strip_powersums_o4

   !> Strip lift through order 4: packed spatial jet of S = -log PD
   !>
   !> - same algebra as `cfc_spatial_eval` at `max_deriv = 4`; writes `lsf0`, `lsf1_r`, `lsf2_rr`,
   !>   `lsf3_rrr`, `lsf4_rrrr` point-last, symmetric groups packed (sorted tuples, last index
   !>   slowest)
   !> - `np <= cfc_strip_width(4) = 8`; requires `pd0 > 0` for every point (escape is the caller's
   !>   test on the stage-1 sums), no zero fill
   !>
   !> @param[in]    np           Live points in the strip, 1..8
   !> @param[in]    pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[in]    pd1          Order-1 pseudo-density jet of the strip, point-first; second index
   !>                            the Cartesian component
   !> @param[in]    pd2          Order-2 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-2 tuple
   !> @param[in]    pd3          Order-3 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-3 tuple
   !> @param[in]    pd4          Order-4 pseudo-density jet of the strip, point-first; second index
   !>                            the packed slot of the sorted rank-4 tuple
   !> @param[out]   lsf0         Level-set value S = -log(PD)
   !> @param[out]   lsf1_r       Spatial gradient dS/dr
   !> @param[out]   lsf2_rr      Spatial Hessian, packed (6)
   !> @param[out]   lsf3_rrr     Third spatial derivative, packed (10)
   !> @param[out]   lsf4_rrrr    Fourth spatial derivative, packed (15)
   pure module subroutine cfc_strip_lift_o4(np, pd0, pd1, pd2, pd3, pd4, lsf0, lsf1_r, lsf2_rr, &
                                     lsf3_rrr, lsf4_rrrr)
      !> Live points in the strip, 1..8
      integer, intent(in) :: np
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(in) :: pd0(np)
      !> Order-1 pseudo-density jet of the strip, point-first; second index the Cartesian component
      real(wp), intent(in) :: pd1(np, 3)
      !> Order-2 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-2 tuple
      real(wp), intent(in) :: pd2(np, 6)
      !> Order-3 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-3 tuple
      real(wp), intent(in) :: pd3(np, 10)
      !> Order-4 pseudo-density jet of the strip, point-first; second index the packed slot of the
      !> sorted rank-4 tuple
      real(wp), intent(in) :: pd4(np, 15)
      !> Level-set value S = -log(PD)
      real(wp), intent(out) :: lsf0(np)
      !> Spatial gradient dS/dr
      real(wp), intent(out) :: lsf1_r(3, np)
      !> Spatial Hessian, packed (6)
      real(wp), intent(out) :: lsf2_rr(6, np)
      !> Third spatial derivative, packed (10)
      real(wp), intent(out) :: lsf3_rrr(10, np)
      !> Fourth spatial derivative, packed (15)
      real(wp), intent(out) :: lsf4_rrrr(15, np)
   end subroutine cfc_strip_lift_o4

   !> Strip adjoint lift: lambda of the weighted jet w.r.t. the order-2 sums
   !>
   !> - per point, `lam[I] = d/d pd[I]` of `w0 S + w1 . dS/dr + W : d2S/dr2` with `S =
   !>   -log PD`: the algebra of `cfc_strip_lift_o2` differentiated once in the sums, `I`
   !>   over the orbit representatives (1 + 3 + 6)
   !> - reads `pd0`, `pd1`, `pd2` of `cfc_strip_powersums_o2` only; `lam` has their
   !>   layout, point-first, packed slots
   !> - seeds hold the symmetric part of the (3, 3) Hessian seed, packed (xx, xy, yy, xz,
   !>   yz, zz) point-last, entries not doubled; the factor 2 of the off-diagonal pairs is
   !>   applied here
   !> - feeds `cfc_strip_vjp`; a non-symmetric seed is the caller's to symmetrize (`W = (W
   !>   + W^T)/2`), the second derivative of `S` is symmetric
   !> - `np <= cfc_strip_width(2) = 32`; requires `pd0 > 0` for every point (escape is the
   !>   caller's test on the stage-1 sums), no zero fill
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    pd0          Order-0 pseudo-density jet of the strip, point-first
   !> @param[in]    pd1          Order-1 pseudo-density jet of the strip, point-first;
   !>                            second index the Cartesian component
   !> @param[in]    pd2          Order-2 pseudo-density jet of the strip, point-first;
   !>                            second index the packed slot of the sorted rank-2 tuple
   !> @param[in]    w0           Seed of the level-set value S
   !> @param[in]    w1           Seed of the spatial gradient dS/dr, point-last
   !> @param[in]    w2           Symmetric part of the (3, 3) seed of d2S/dr2, packed (xx,
   !>                            xy, yy, xz, yz, zz), point-last; entries not doubled
   !> @param[out]   lam0         Adjoint of the order-0 sums, layout of `pd0`
   !> @param[out]   lam1         Adjoint of the order-1 sums, layout of `pd1`
   !> @param[out]   lam2         Adjoint of the order-2 sums, layout of `pd2`
   pure module subroutine cfc_strip_lift_vjp(np, pd0, pd1, pd2, w0, w1, w2, lam0, lam1, lam2)
      !> Live points in the strip, 1..32
      integer, intent(in) :: np
      !> Order-0 pseudo-density jet of the strip, point-first
      real(wp), intent(in) :: pd0(np)
      !> Order-1 pseudo-density jet of the strip, point-first; second index the Cartesian
      !> component
      real(wp), intent(in) :: pd1(np, 3)
      !> Order-2 pseudo-density jet of the strip, point-first; second index the packed
      !> slot of the sorted rank-2 tuple
      real(wp), intent(in) :: pd2(np, 6)
      !> Seed of the level-set value S
      real(wp), intent(in) :: w0(np)
      !> Seed of the spatial gradient dS/dr, point-last
      real(wp), intent(in) :: w1(3, np)
      !> Symmetric part of the (3, 3) seed of d2S/dr2, packed (xx, xy, yy, xz, yz, zz),
      !> point-last; entries not doubled
      real(wp), intent(in) :: w2(6, np)
      !> Adjoint of the order-0 sums, layout of `pd0`
      real(wp), intent(out) :: lam0(np)
      !> Adjoint of the order-1 sums, layout of `pd1`
      real(wp), intent(out) :: lam1(np, 3)
      !> Adjoint of the order-2 sums, layout of `pd2`
      real(wp), intent(out) :: lam2(np, 6)
   end subroutine cfc_strip_lift_vjp

   !> Strip nuclear vjp: gradient w.r.t. every atom center of the weighted jet
   !>
   !> - `res(:, a, p) = d/dR_a [w0 S + w1 . dS/dr + W : d2S/dr2](r_p)` over the shared
   !>   atom list, which is `sum_I lam[I] d pd[I]/dR_a` with `lam` from
   !>   `cfc_strip_lift_vjp`; the contraction with `lam` is done before emission, so the
   !>   nuclear tensors `qn*` of the per-point path are never formed: a (pair, point)
   !>   yields the six numbers of its two centers
   !> - atomic term `w_A e_A`: `-sum_I lam[I] d^(I+s) (w e)(||d_A||)` over the radial
   !>   chain of the SvdW generator, one set of three per atom; pair term `f(q_AB) g_A
   !>   g_B`, `g = w u`: `d/dR_A = -d/dd_A`, `d/dr = d/dd_A + d/dd_B`, contracted over the
   !>   Cartesian jets of `g_A`, `g_B` (rank <= 3) and the `q` derivatives `f_k`, k <= 3
   !> - reads of the stage-1 cache `geo` of `cfc_strip_powersums_o2` for the same strip:
   !>   slots 1-4 (`r - R_A`, `1/||r - R_A||`) and the gate-masked Cartesian jets of `g_A`
   !>   through rank 2 (slots 5-14); recomputed here: the inline exp and the fade of every
   !>   survivor (the cache keeps neither), the rank-3 jet of `g` (into `jet3`, the cache
   !>   stops at rank 2) and the atomic-term jets, which stage 1 never caches
   !> - `points`, `centers`, `radii`, `bound_sq`, `a1`, `a2`, `c_par`, `fade_start`,
   !>   `fade_width` must be the stage-1 arguments of that strip, so the classification,
   !>   the gate mask and the survivor list are the same
   !> - classification as stage 1: per atom against the strip box, outside the gate sphere
   !>   or beyond the band is skipped (`res` exactly 0, no cache row is read), plateau
   !>   (`sigma_max <= fade_start`) runs the plain radial jets, else the faded jets (`w =
   !>   1 - S_4(t)`); pairs of survivors whose reach spheres are disjoint are skipped
   !>   (exact zeros); plateau and faded pairs share one body
   !> - a pair with no point inside both gates is skipped too: every lane of it has one
   !>   factor's jet gate-masked to exact zeros; the gate bits of each survivor (one per
   !>   point) are kept in `res(1, i, 1)` until `res` is written at the end
   !> - on-nucleus pairs (`||r_p - R_A|| == 0`) are counted in `n_zero`, not handled:
   !>   `1/x` is Inf there and the `res` lane is not finite (Inf or NaN, also for the
   !>   other atoms of that point); the caller re-evaluates such a strip on the guarded
   !>   per-point path
   !> - a point whose stage-1 column was unreached (`pd0 <= 0`, not finite) has no valid
   !>   `lam`; its `res` column is the caller's to discard, nothing here tests it
   !> - atom pass (atom outer, point inner) then pair pass (pair outer, point inner); the
   !>   point loops carry no branch and no call; the result is accumulated point-first in
   !>   `gacc` (the atom of the outer pair loop in a local) and transposed to `res` once
   !>   at the end
   !> - each pair-pass accumulation is 3 independent partial sums over round-robin terms
   !>   (`ps<n>_<k>`), added to the accumulator at the end of the point iteration: one
   !>   serial chain per output is latency bound
   !> - each CSE temporary emitted just before its first reader
   !> - `np <= cfc_strip_width(2) = 32`, the width stage 1 runs at in this pipeline (order
   !>   2); no `nat`-sized automatic array, `jet3`, `gacc`, `sv`, `reach` come from the
   !>   caller
   !>
   !> @param[in]    np           Live points in the strip, 1..32
   !> @param[in]    points       Strip points, point-last; those of stage 1
   !> @param[in]    nat          Atoms in the shared list
   !> @param[in]    centers      Atom centers of the shared list
   !> @param[in]    radii        Atom radii R_a
   !> @param[in]    bound_sq     Squared screening reach per atom, the gate's
   !>                            `cand_screen` row
   !> @param[in]    a1           Atomic-term exponent a1 (<= 0)
   !> @param[in]    a2           Pair-term exponent a2 (<= 0)
   !> @param[in]    c_par        Pair-term coupling constant c
   !> @param[in]    fade_start   Fade onset in sigma_A, as in stage 1; huge(0.0_wp) for no
   !>                            fade
   !> @param[in]    fade_width   Fade band width in sigma_A, as in stage 1; > 0
   !> @param[in]    lam0         Adjoint of `pd0` from `cfc_strip_lift_vjp`
   !> @param[in]    lam1         Adjoint of `pd1`, layout of `pd1`
   !> @param[in]    lam2         Adjoint of `pd2`, packed, layout of `pd2`
   !> @param[in]    geo          Stage-1 geometry cache of `cfc_strip_powersums_o2` for
   !>                            this strip, read for the survivors only (slots in that
   !>                            routine's notes)
   !> @param[out]   jet3         Workspace: rank-3 Cartesian jet `d^J g_A/dr^J` of every
   !>                            survivor, packed, point-first; a buffer of
   !>                            `cfc_strip_geo_len` doubles per atom serves by sequence
   !>                            association
   !> @param[out]   gacc         Workspace: accumulated gradient of every survivor,
   !>                            point-first; a buffer of `cfc_strip_geo_len` doubles per
   !>                            atom serves by sequence association
   !> @param[out]   sv           Survivor workspace: list indices of the atoms kept by the
   !>                            strip classification, leading entries written
   !> @param[out]   reach        Reach workspace: pair-test radius of each survivor,
   !>                            leading entries written
   !> @param[out]   res          Nuclear gradient of the weighted jet per atom and point,
   !>                            point-last; exact zeros out of reach
   !> @param[out]   n_zero       Point-atom pairs at zero distance among the survivors;
   !>                            nonzero means `res` is not finite there
   pure module subroutine cfc_strip_vjp(np, points, nat, centers, radii, bound_sq, a1, a2, &
                                 c_par, fade_start, fade_width, lam0, lam1, lam2, geo, &
                                 jet3, gacc, sv, reach, res, n_zero)
      !> Strip width of this pipeline (the order-2 width)
      integer, parameter :: nps = cfc_strip_width(2)
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
      !> Squared screening reach per atom, the gate's `cand_screen` row
      real(wp), intent(in) :: bound_sq(nat)
      !> Atomic-term exponent a1 (<= 0)
      real(wp), intent(in) :: a1
      !> Pair-term exponent a2 (<= 0)
      real(wp), intent(in) :: a2
      !> Pair-term coupling constant c
      real(wp), intent(in) :: c_par
      !> Fade onset in sigma_A, as in stage 1; huge(0.0_wp) for no fade
      real(wp), intent(in) :: fade_start
      !> Fade band width in sigma_A, as in stage 1; > 0
      real(wp), intent(in) :: fade_width
      !> Adjoint of `pd0` from `cfc_strip_lift_vjp`
      real(wp), intent(in) :: lam0(np)
      !> Adjoint of `pd1`, layout of `pd1`
      real(wp), intent(in) :: lam1(np, 3)
      !> Adjoint of `pd2`, packed, layout of `pd2`
      real(wp), intent(in) :: lam2(np, 6)
      !> Stage-1 geometry cache of `cfc_strip_powersums_o2` for this strip, read for the
      !> survivors only (slots in that routine's notes)
      real(wp), intent(in) :: geo(nps, 14, nat)
      !> Workspace: rank-3 Cartesian jet `d^J g_A/dr^J` of every survivor, packed,
      !> point-first; a buffer of `cfc_strip_geo_len` doubles per atom serves by sequence
      !> association
      real(wp), intent(out) :: jet3(nps, 10, nat)
      !> Workspace: accumulated gradient of every survivor, point-first; a buffer of
      !> `cfc_strip_geo_len` doubles per atom serves by sequence association
      real(wp), intent(out) :: gacc(nps, 3, nat)
      !> Survivor workspace: list indices of the atoms kept by the strip classification,
      !> leading entries written
      integer, intent(out) :: sv(nat)
      !> Reach workspace: pair-test radius of each survivor, leading entries written
      real(wp), intent(out) :: reach(nat)
      !> Nuclear gradient of the weighted jet per atom and point, point-last; exact zeros
      !> out of reach
      real(wp), intent(out) :: res(3, nat, np)
      !> Point-atom pairs at zero distance among the survivors; nonzero means `res` is not
      !> finite there
      integer, intent(out) :: n_zero
   end subroutine cfc_strip_vjp

   end interface

end module moist_cavity_drop_lsf_cfc_kernel
