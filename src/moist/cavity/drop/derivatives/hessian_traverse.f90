!> Single grid traversal of the DROP surface Hessian
!>
!> The nuclear gradient of `nuclear.f90` is `J^T omega` with `omega` the
!> accumulated surface adjoint. Its directional derivative splits into
!>
!>     d/dv [ J^T omega ]  =  (dJ^T/dv) omega  +  J^T (d omega/dv)
!>                            ^ fixed channel     ^ response channel
!>
!> Both terms are built here, in **one** walk of the projected grid, because
!> everything before their contractions is the same computation: one
!> [[drop_point_prologue]] per point (level-set `prepare`, jets, seed state,
!> bordered KKT factorization), one application of the **same 16 basis seeds**
!> ([[point_seed_basis]]), one materialisation of the point's jet tensors
!> ([[drop_field_jet_point]]) and one collection of the anchor's iSwiG
!> neighbourhood. A weight set, or a direction, enters only at a contraction.
!>
!> The two public entry points and the multi-branch refusal live in
!> `hessian.f90`; the half-accessors [[get_surface_hessian_fixed_drop]],
!> [[get_surface_hessian_fixed_dirs_drop]] and
!> [[get_surface_hessian_response_drop]] are thin wrappers that enable one
!> channel alone, kept because their unit suites verify the halves in isolation.
!>
!>
!> ## Why sharing the seed applications is legitimate
!>
!> [[apply_seed]] reads the per-point forward state and the point's own
!> bordered solve, and nothing else: the seed directions are constant matrices,
!> the induced point motion comes out of the KKT batch, and no surface adjoint
!> reaches it. So one application of the basis serves any number of weight sets
!> at that point, which is why [[seed_jet_basis]] and [[seed_anchor]] are
!> offered split into an `_apply` and a `_contract` half. The fixed channel
!> needs one thing more, the derived-state tangent `dstate` of every seed, and
!> asks for it through the optional argument the `_apply` halves carry.
!>
!>
!> ## The fixed-adjoint channel, and its two modes
!>
!> The second derivative of the surface map contracted against adjoints held
!> *fixed*. Per grid point, the term along one nuclear direction `v` is one
!> pass of the second-order chain, [[fixed_direction_chain]]: the jet tangent
!> along `v` at the frozen point, the projection riding on it through the same
!> bordered system the seeds use, the tangent of the normal fold, the tangent
!> of the 16 seed right-hand sides, 16 [[apply_seed_tangent]] calls, and the
!> field-row tangent [[drop_field_tangent_dir]]. Its output is one gradient
!> column: a field row over the active atoms, three anchor entries on the
!> owner, and -- outside the chain, because it needs no seed at all -- the
!> switching block contracted with `v`. The chain is linear in `v`, and the
!> channel is offered in two modes that differ only in which `v` they run it
!> for:
!>
!>  * **`drop_fixed_rank4`** (the dense Hessian). The projected point depends
!>    on the nuclei only through the level set's active atoms and through the
!>    owner, so the chain is run for the `3 * |active union owner|` Cartesian
!>    unit directions of that local set and the columns are accumulated into
!>    the direction-free rank-4 block `(3, nsph, 3, nsph)`. One traversal serves
!>    every direction; the price is a per-thread accumulator, which is sparse
!>    for the reason [[drop_hess_sparse_type]] gives. The jet tangent of a unit
!>    direction is a *column* of the point's tensors ([[drop_field_jet_column]]),
!>    and the explicit nuclear motion of every column is one block of the
!>    level set, `vjp_f2_rArB`, formed once per point; the chain makes no LSF
!>    call per direction. For SvdW that block is a rank-52 product of per-atom
!>    quantities, `O(n_active)` kernel work and one `gemm`; a level set without
!>    a factorised form inherits the column-by-column default and pays
!>    `3 n_local` passes of `hvp_jet_rA` per point. What is left per direction
!>    is the second-order chain itself, `O(1)` in the active count.
!>
!>  * **`drop_fixed_per_dir`** (Hessian-vector products). The chain is run once
!>    per *supplied* direction, with the jet tangent formed by contracting the
!>    same tensors with the direction ([[drop_field_jet_tangent]]), and the
!>    column lands straight in the response channel's `(3, nsph, ndir)`
!>    per-thread buffer, with the explicit nuclear motion of each direction from
!>    `hvp_jet_rA`. Cost one `O(n_active)` accessor pass per direction and
!>    point, against the rank-4 mode's per-point block and `3 n_local` chains,
!>    so it is the right mode for a few directions and the wrong one for a
!>    basis; `hessian.f90` makes that choice.
!>
!> Both modes read the same weights of the point, [[drop_point_weights_type]]:
!> the normal fold and the 13 jet-seed contractions of the *fixed* adjoints.
!>
!> The nine single-entry Hessian jet seeds are asymmetric, and
!> [[apply_seed_tangent]] is exact only on a symmetric seed: its per-seed error
!> is antisymmetric in `(i, j)`. The chain therefore never reads one entry of
!> `dw_lsf2` on its own; that tangent reaches the field row only through the
!> contraction against `lsf3_rr_rA`, which is symmetric in its two spatial
!> indices, and an antisymmetric error contracted with a symmetric weight
!> cancels to machine zero.
!>
!>
!> ## The adjoint-response channel
!>
!> Here the **primal map is held fixed and the adjoints move**. The raw host
!> adjoints of a `cavity_surface_adjoint_type` never move; what moves is
!> [[prepare_surface_weights]]'s *folding* of them, `eff(R)`. Writing the
!> gradient as `G(R) = Phi(R) . eff(R)`, `G` is linear in `eff`, so the
!> response term is [[get_surface_gradient_drop]] run with `d(eff)/dv` in place
!> of `eff`, seed for seed. Its passes are:
!>
!>  1. **Forward tangent** -- [[get_surface_tangent_drop]] returns `d_a`,
!>     `d_wleb`, `d_xi0` and `d_wbranch`, one column per direction.
!>  2. **Weight tangent** -- [[prepare_surface_weights_tangent]] turns those
!>     into `d(eff)`, which carries `w_xi`, `w_f` and `branch_phi_adj` and
!>     nothing else: `w_xyz`, `w_n`, `w_k1` and `w_k2` are copies of the raw
!>     adjoints, so their tangents vanish identically. Consequently the normal
!>     channel of the contraction sees the hard zero `w_xyz_zero`, and the
!>     curvature channel is off (`deff%have_wk = .false.`). This pass runs
!>     **serially, one call per direction**, because
!>     [[branch_phi_adj_tangent]] reduces over contiguous anchor groups and a
!>     group split across threads would corrupt that reduction silently.
!>  3. **Contraction** -- the grid loop below: the 13 jet contractions of
!>     `deff(jdir)` give a weight set, [[drop_field_jet_contract]] turns it into
!>     the field row off the point's jet tensors, and the anchor and switching
!>     rows follow the gradient path. No LSF call is made per direction.
!>
!> This channel is first order in the seed chain -- it reads `res` and never
!> `dres` -- so the asymmetric-seed caveat above does not touch it; the seeds
!> still leave only through the symmetric `lsf3_rr_rA` contraction.
!>
!>
!> ## Direction blocks
!>
!> Passes 1 and 2 materialise seven `(ngrid, ndir)` grid arrays and the column
!> accumulator is `(3, nsph, ndir, nthreads)`; nothing above this submodule
!> bounds `ndir`. The traversal is therefore blocked over directions -- see
!> `drop_hvp_chunk_dirs` -- and **the whole grid is re-traversed per block**.
!> Two invariants follow:
!>
!>  * the rank-4 fixed channel is direction free and accumulates **in the first
!>    block only**, `fixed_rank4_here = ilo == 1`; every fixed-only step reads
!>    that flag. Later blocks neither accumulate it twice nor pay for the
!>    fourth-order jet buffer, the state tangents or the switching block;
!>  * the per-direction fixed channel is per direction and runs **in every
!>    block**.
!>
!> A direction's column receives the same additions in the same order whatever
!> the blocking, so a blocked run reproduces an unblocked one to the bit;
!> `hvp_direction_chunking` in the end-to-end suite holds that down. That
!> guarantee is why the prologue *configuration* -- the level set's derivative
!> order and the curvature request -- is the same in every block of a
!> traversal, although a block without a fixed channel reads nothing above
!> order 3 and would be cheaper to prepare at that order. Measured: cloning
!> the later blocks' level sets at order 3 moves their response columns by one
!> ulp against a single-block run, because the generated atom kernels
!> (`svdw_atom_eval` and its CFC counterpart) branch on the requested level and
!> the branches schedule the lower orders differently. The saving is only
!> available once those kernels round their lower orders level-independently;
!> until then reproducibility across blocking wins.
!>
!> Passes 1 and 2 cannot be fused into the contraction walk: pass 2's group
!> reduction needs the whole grid's pass-1 output before any point is
!> contracted, so each block is one pass-1 walk and one contraction walk.
!>
!>
!> ## Configuration and failure contract
!>
!> `eff` is a [[drop_surface_weights_type]] that [[prepare_surface_weights]] has
!> already produced; `hessian.f90` folds once and drives both channels off the
!> same object. The response channel needs the raw `acc` as well, because
!> pass 2 differentiates the fold out of the raw channels and the primal `eff`
!> together. What neither channel offers -- the multi-branch second-order term
!> -- is refused at the public entry points, not here: the response channel is
!> correct on a branched grid on its own, and its suite finite-differences one.
!>
!> Both accumulators are *added* to and left untouched on failure, except for
!> the one case the blocking makes unavoidable: a multi-block run failing in a
!> later block has already reduced the earlier blocks. The column accumulator
!> is restored from a copy taken up front, and the rank-4 block is only ever
!> reached through `hessian.f90`, which stages it in a local buffer.
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_hessian_traverse
!$ use omp_lib, only: omp_get_thread_num
   use, intrinsic :: iso_fortran_env, only: int64
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type, &
      & drop_point_scratch_type
   use moist_cavity_drop_derivatives_kernel, only: &
      & drop_seed_state_type, drop_seed_result_type, drop_seed_state_tangent_type, &
      & drop_seed_input_tangent_type, drop_seed_result_tangent_type, &
      & drop_surface_weights_type, apply_seed, apply_seed_tangent, &
      & seed_weight_tol, seed_contribution
   use moist_cavity_drop_derivatives_seeds, only: drop_n_jet_seeds, &
      & drop_n_point_seeds, seed_normal_channel, fill_seed_basis, scatter_jet_weight, &
      & seed_contribution_tangent, seed_dzero1, seed_dzero2, &
      & seed_jet_basis_apply, seed_jet_basis_contract, &
      & seed_anchor_apply, seed_anchor_contract
   use moist_cavity_drop_derivatives_field_tangent, only: drop_field_tangent_point, &
      & drop_field_jet_point, drop_field_jet_contract, drop_field_jet_tangent, &
      & drop_field_jet_column, drop_field_f4_fold, drop_field_tangent_dir, &
      & drop_field_tangent_work_type
   use moist_cavity_drop_derivatives_weights_tangent, only: prepare_surface_weights_tangent
   use moist_cavity_drop_gaussian_scatter, only: scatter_iswig_block_indexed, &
      & contract_iswig_block
   implicit none(type, external)

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Initial entry capacity and bucket count of a per-thread accumulator
   !>
   !> Neither is a bound: both grow on demand, and both are trimmed at birth to
   !> a molecule that cannot fill them. The bucket count stays a power of two --
   !> which is what makes the mask in [[hess_sparse_entry]] a legal index.
   integer, parameter :: hess_init_ent = 256
   integer, parameter :: hess_init_tab = 1024

   !> Load factor above which the bucket table doubles, as `num/den`
   integer, parameter :: hess_load_num = 7, hess_load_den = 10

   !> FNV-1a (32-bit) constants of [[hess_pair_hash]]; every product is masked
   !> back to 32 bits, so no intermediate leaves the signed 64-bit range
   integer(int64), parameter :: fnv_offset = 2166136261_int64
   integer(int64), parameter :: fnv_prime = 16777619_int64
   integer(int64), parameter :: mask32 = int(z'FFFFFFFF', int64)

   !> Sparse per-thread accumulator of the rank-4 fixed-adjoint Hessian
   !>
   !> One thread's share of `(3, nsph, 3, nsph)`, held as the `(3, 3)` blocks of
   !> the atom pairs it actually reaches. A grid point reaches the level set's
   !> active atoms, its anchor's owner and the switching influence set of that
   !> anchor -- all of them local neighbourhoods -- so the pair count grows like
   !> `nsph * k` and not like `nsph^2`, where a dense per-thread slab would be
   !> `9 nsph^2` doubles: 1.15 GB at a thousand atoms on sixteen threads.
   !>
   !> **This container reproduces a dense one to the bit, and that is a
   !> property of how it is used, not of the arithmetic.** Two orders carry it:
   !>
   !>  1. *Accumulation order.* Every contribution is added in place, in grid
   !>     order, exactly where a dense slab's `+=` would stand. One pair has one
   !>     entry -- the hash guarantees that -- so an element receives its
   !>     contributions in sequence. Growth preserves it: the entry arrays are
   !>     copied by `move_alloc` so indices never move, and the bucket table is
   !>     rebuilt from the entries rather than the reverse.
   !>
   !>  2. *Merge order.* [[hess_sparse_reduce]] is called in the fixed
   !>     `1 .. nthreads` sequence and adds each thread's block to the
   !>     destination in one add per element. A pair no thread touched is
   !>     skipped, where a dense reduction would add an exact zero to it.
   !>
   !> The bucket table is open addressing with linear probing, the idiom
   !> [[timer_type]] already uses in this codebase, with entry indices as values
   !> and `0` for an empty slot.
   type :: drop_hess_sparse_type
      !> Entries in use
      integer :: nent = 0
      !> Pairs that exist at all, `nsph^2` clamped to the integer range
      !>
      !> The entry arrays double, so a container that ends up holding every
      !> pair would otherwise carry up to twice the dense slab. Capping the
      !> capacity here bounds the overshoot at the point where the two meet,
      !> and is exact: a distinct pair beyond this one cannot be asked for
      integer :: maxent = 0
      !> Row and column atom of each entry
      integer, allocatable :: pair_i(:), pair_j(:)
      !> Accumulated `(3, 3)` block of each entry, `(row axis, column axis, entry)`
      real(wp), allocatable :: blocks(:, :, :)
      !> Open-addressing bucket table; values are entry indices, `0` is empty
      integer, allocatable :: htab(:)
   end type drop_hess_sparse_type

   !> The 16 basis seeds of one grid point and their responses
   !>
   !> Filled once per point by [[fill_seed_basis]] and [[point_seed_basis]] and
   !> read by both channels; `dstate` is filled for the fixed channel only. A
   !> per-thread object, so that the parallel region names it once.
   type :: drop_point_seeds_type
      !> Seed perturbations of the level-set gradient and Hessian
      real(wp) :: dlsf1(3, drop_n_point_seeds) = 0.0_wp
      real(wp) :: dlsf2(3, 3, drop_n_point_seeds) = 0.0_wp
      !> Induced point motion and multiplier change of each seed
      real(wp) :: x(4, drop_n_point_seeds) = 0.0_wp
      !> Linear response of each seed
      type(drop_seed_result_type) :: res(drop_n_point_seeds)
      !> Tangent of the derived state along each seed (fixed channel only)
      type(drop_seed_state_tangent_type) :: dstate(drop_n_point_seeds)
   end type drop_point_seeds_type

   !> One weight set contracted at a grid point
   !>
   !> The point-local level-set adjoints built from the 13 jet seeds, the
   !> effective position adjoint every seed sees, and the normal fold it was
   !> built from. The fixed channel fills all of it from `eff`; the response
   !> channel fills the three level-set adjoints from `deff(jdir)` and leaves
   !> the rest at zero, which is exact for it.
   type :: drop_point_weights_type
      !> Point-local level-set adjoint weights
      real(wp) :: w_lsf0 = 0.0_wp, w_lsf1(3) = 0.0_wp, w_lsf2(3, 3) = 0.0_wp
      !> Effective position adjoint and the normal fold it is built from
      real(wp) :: w_xyz(3) = 0.0_wp, normal_grad(3) = 0.0_wp, nwn = 0.0_wp
   end type drop_point_weights_type

   !> Per-direction temporaries of [[fixed_direction_chain]]
   !>
   !> Held in one per-thread object rather than as locals of the chain so that
   !> the derived-type temporaries are not default-initialised on every call.
   type :: drop_fixed_dir_scratch_type
      !> Right-hand side and solution of the directional projection response
      real(wp) :: rhs_v(4, 1) = 0.0_wp, dr_v(3) = 0.0_wp, dl_v = 0.0_wp
      !> Response, derived-state tangent and input tangent of the direction
      type(drop_seed_result_type) :: res_v
      type(drop_seed_state_tangent_type) :: dstate_v
      type(drop_seed_input_tangent_type) :: dinp_v
      !> Second-order response of one seed along the direction
      type(drop_seed_result_tangent_type) :: dres
      !> Tangent of the bordered KKT system and of its 16 right-hand sides
      real(wp) :: dH_lag(3, 3, 1) = 0.0_wp, dg_tot(3, 1) = 0.0_wp
      real(wp) :: dseed_x(4, drop_n_point_seeds) = 0.0_wp
      !> Tangents of the normal fold and of the level-set adjoint weights
      real(wp) :: dnormal_grad(3) = 0.0_wp, dw_xyz(3) = 0.0_wp
      real(wp) :: dw_lsf0 = 0.0_wp, dw_lsf1(3) = 0.0_wp, dw_lsf2(3, 3) = 0.0_wp
   end type drop_fixed_dir_scratch_type

   !> The rank-4 mode's local direction set and pair table of one point
   !>
   !> `pair_ent` and `vdir` are grown from what a point actually needs, never
   !> sized to `nsph^2`: a molecule-sized pair map is the quadratic the sparse
   !> accumulator exists to avoid.
   type :: drop_rank4_scratch_type
      !> Atoms of the local direction set, and the owner's position in it
      integer :: ndir_atom = 0, owner_loc = 0
      integer, allocatable :: dir_atoms(:)
      !> Accumulator entry of every `(row, column)` pair the point reaches
      integer, allocatable :: pair_ent(:, :)
      !> The current Cartesian unit direction, molecule sized and otherwise zero
      real(wp), allocatable :: vdir(:, :)
      !> Weighted mixed nuclear Hessian block of the point over its active
      !> slots, `(3 n_active, 3 n_active)` with the Cartesian component fastest.
      !> Sized to the point, so the level set writes it in place
      real(wp), allocatable :: mblk(:, :)
   end type drop_rank4_scratch_type

   !> The anchor's iSwiG neighbourhood as one point sees it
   type :: drop_swi_scratch_type
      !> Whether the fixed and the response channel want the neighbourhood
      logical :: fixed = .false., resp = .false.
      !> Collected switching factor, and the response channel's owner row
      real(wp) :: f0 = 0.0_wp, owner_row(3) = 0.0_wp, dxi = 0.0_wp
      !> Influence-set size of the second-derivative block
      integer :: n = 0
      !> Influence set, its second-derivative block, and its pair entries
      !>
      !> Grown on demand from `n_nb + 1`: the block is quadratic in the
      !> influence set, and a molecule-sized one would put back per thread the
      !> quadratic the sparse accumulator avoids
      integer, allocatable :: idx(:), ent(:, :)
      real(wp), allocatable :: blk(:, :, :, :)
   end type drop_swi_scratch_type

contains

   !* ================================================================================= *!
   !*                            The half-accessor wrappers                             *!
   !* ================================================================================= *!

   !> The rank-4 fixed-adjoint half of the surface Hessian
   !>
   !> Accumulates `(dJ^T/dv) omega` for the energy whose folded surface adjoints
   !> `eff` are, as the direction-free block, and nothing else. The result is
   !> *added* to `hessian`, and left untouched when anything fails.
   !>
   !> Kept as an entry point because `test_cavity_drop_hessian_fixed`
   !> finite-differences this half on its own; the composition in `hessian.f90`
   !> asks for both channels in one traversal instead. It carries no weight
   !> guard: `eff` arrives folded, and the term a live `w_a` or `w_w` adds is
   !> exactly what the response channel computes.
   !>
   !> @param[in]    self    DROP cavity instance (must hold a projected grid)
   !> @param[in]    eff     Folded surface adjoints, held fixed
   !> @param[inout] hessian Nuclear-Hessian accumulator (3, nsph, 3, nsph)
   !> @param[out]   error   Error object, allocated on failure
   module subroutine get_surface_hessian_fixed_drop(self, eff, hessian, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Effective primitive surface adjoints, as [[prepare_surface_weights]]
      !> returned them
      type(drop_surface_weights_type), intent(in) :: eff
      !> Nuclear-Hessian accumulator
      real(wp), intent(inout) :: hessian(:, :, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (any(shape(hessian) /= [ndim, self%nsph, ndim, self%nsph])) then
         call fatal_error(error, "get_surface_hessian_fixed_drop: hessian shape mismatch")
         return
      end if
      if (self%ngrid <= 0) return

      call drop_hessian_traverse(self, eff, drop_fixed_rank4, .false., &
                                 "get_surface_hessian_fixed_drop", &
                                 hess_fixed=hessian, error=error)

   end subroutine get_surface_hessian_fixed_drop

   !> The per-direction fixed-adjoint half of the surface Hessian
   !>
   !> Accumulates `(dJ^T/dv) omega` along each supplied direction, one gradient
   !> column per direction, and nothing else. The result is *added* to `hvp`
   !> and left untouched when anything fails. Mathematically the rank-4 half
   !> contracted with `dirs`; computationally the second-order chain run once
   !> per direction. The fixed suite cross-checks the two.
   !>
   !> @param[in]    self  DROP cavity instance (must hold a projected grid)
   !> @param[in]    eff   Folded surface adjoints, held fixed
   !> @param[in]    dirs  Nuclear directions `(3, nsph, ndir)`
   !> @param[inout] hvp   Hessian-vector accumulator `(3, nsph, ndir)`
   !> @param[out]   error Error object, allocated on failure
   module subroutine get_surface_hessian_fixed_dirs_drop(self, eff, dirs, hvp, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Folded surface adjoints, as [[prepare_surface_weights]] returned them
      type(drop_surface_weights_type), intent(in) :: eff
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Hessian-vector accumulator
      real(wp), intent(inout) :: hvp(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_direction_set(self, dirs, hvp, "get_surface_hessian_fixed_dirs_drop", error)
      if (allocated(error)) return
      if (self%ngrid <= 0) return

      call drop_hessian_traverse(self, eff, drop_fixed_per_dir, .false., &
                                 "get_surface_hessian_fixed_dirs_drop", &
                                 dirs=dirs, hvp=hvp, error=error)

   end subroutine get_surface_hessian_fixed_dirs_drop

   !> The adjoint-response half of the surface Hessian
   !>
   !> Accumulates `J^T (d omega/dv)` for the energy whose raw surface adjoints
   !> `acc` holds, one gradient column per nuclear direction, and nothing else.
   !> The result is *added* to `hvp` and left untouched when anything fails.
   !>
   !> This half needs **both** forms of the adjoints: `eff` is the primal fold,
   !> `acc` is what that fold was built from, and pass 2 differentiates the one
   !> out of the other.
   !>
   !> No branch guard, deliberately: the public accessors refuse a multi-branch
   !> grid because the *composite* is short a second-order branch term; this
   !> half on its own is not, and its suite finite-differences a branched grid.
   !>
   !> @param[in]    self  DROP cavity instance (must hold a projected grid)
   !> @param[in]    acc   Raw surface-observable adjoints, held fixed
   !> @param[in]    eff   Folded surface adjoints of the base geometry
   !> @param[in]    dirs  Nuclear directions `(3, nsph, ndir)`
   !> @param[inout] hvp   Hessian-vector accumulator `(3, nsph, ndir)`
   !> @param[out]   error Error object, allocated on failure
   module subroutine get_surface_hessian_response_drop(self, acc, eff, dirs, hvp, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Raw surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Folded surface adjoints, as [[prepare_surface_weights]] returned them
      type(drop_surface_weights_type), intent(in) :: eff
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Hessian-vector accumulator
      real(wp), intent(inout) :: hvp(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_surface_adjoint(self, acc, "get_surface_hessian_response_drop", error)
      if (allocated(error)) return
      call check_direction_set(self, dirs, hvp, "get_surface_hessian_response_drop", error)
      if (allocated(error)) return
      if (self%ngrid <= 0) return

      call drop_hessian_traverse(self, eff, drop_fixed_none, .true., &
                                 "get_surface_hessian_response_drop", &
                                 acc=acc, dirs=dirs, hvp=hvp, error=error)

   end subroutine get_surface_hessian_response_drop

   !> Refuse a direction set or a column accumulator of the wrong shape
   !>
   !> Every message is prefixed with the caller's name, as
   !> [[check_surface_adjoint]] does, so a failure names the entry point the
   !> user actually called.
   !>
   !> @param[in]  self    DROP cavity instance
   !> @param[in]  dirs    Nuclear directions, expected `(3, nsph, ndir)`
   !> @param[in]  hvp     Column accumulator, expected `(3, nsph, ndir)`
   !> @param[in]  context Calling routine, used to prefix the diagnostics
   !> @param[out] error   Error object, allocated on a shape mismatch
   module subroutine check_direction_set(self, dirs, hvp, context, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Column accumulator
      real(wp), intent(in) :: hvp(:, :, :)
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Direction count of the request
      integer :: ndir

      if (size(dirs, 1) /= ndim .or. size(dirs, 2) /= self%nsph) then
         call fatal_error(error, context//": dirs must be (3, nsph, ndir)")
         return
      end if
      ndir = size(dirs, 3)
      if (ndir <= 0) then
         call fatal_error(error, context//": no direction supplied")
         return
      end if
      if (size(hvp, 1) /= ndim .or. size(hvp, 2) /= self%nsph .or. size(hvp, 3) /= ndir) then
         call fatal_error(error, context//": hvp must be (3, nsph, ndir)")
         return
      end if
   end subroutine check_direction_set

   !* ================================================================================= *!
   !*                             The merged grid traversal                             *!
   !* ================================================================================= *!

   !> Walk the projected grid once, accumulating into the enabled channels
   !>
   !> Mirrors [[get_surface_gradient_drop]] throughout: same thread setup, same
   !> grid loop, same error latching and the same deterministic reduction. What
   !> differs is the derivative order, the contractions the point feeds, and
   !> the direction blocking.
   !>
   !> `fixed_mode` is one of `drop_fixed_none`, `drop_fixed_rank4` and
   !> `drop_fixed_per_dir` (see the module header). The optional arguments are
   !> the channels' own: `hess_fixed` belongs to the rank-4 fixed channel,
   !> `dirs` and `hvp` to the per-direction fixed channel and to the response
   !> channel, `acc` to the response channel alone. Each must be present when
   !> its channel is enabled, which is checked.
   !>
   !> @param[in]    self           DROP cavity instance (must hold a projected grid)
   !> @param[in]    eff            Folded surface adjoints of the base geometry
   !> @param[in]    fixed_mode     Which form of `(dJ^T/dv) omega` to accumulate, if any
   !> @param[in]    want_response  Accumulate `J^T (d omega/dv)`
   !> @param[in]    context        Calling routine, used to prefix the diagnostics
   !> @param[in]    acc            Raw surface adjoints; required by the response channel
   !> @param[in]    dirs           Nuclear directions; required by every per-direction channel
   !> @param[inout] hess_fixed     Rank-4 fixed-adjoint accumulator `(3, nsph, 3, nsph)`
   !> @param[inout] hvp            Column accumulator `(3, nsph, ndir)`
   !> @param[out]   error          Error object, allocated on failure
   module subroutine drop_hessian_traverse(self, eff, fixed_mode, want_response, context, &
                                           acc, dirs, hess_fixed, hvp, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Folded surface adjoints, as [[prepare_surface_weights]] returned them
      type(drop_surface_weights_type), intent(in) :: eff
      !> Mode of the fixed channel
      integer, intent(in) :: fixed_mode
      !> Whether the response channel is accumulated
      logical, intent(in) :: want_response
      !> Calling routine, so a failure names the entry point the user called
      character(len=*), intent(in) :: context
      !> Raw surface-observable adjoints, held fixed
      type(cavity_surface_adjoint_type), intent(in), optional :: acc
      !> Nuclear directions
      real(wp), intent(in), optional :: dirs(:, :, :)
      !> Rank-4 accumulator of the fixed channel
      real(wp), intent(inout), optional :: hess_fixed(:, :, :, :)
      !> Column accumulator of the per-direction channels
      real(wp), intent(inout), optional :: hvp(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type) :: slots
      !> Per-thread rank-4 buffers, summed deterministically after the region
      type(drop_hess_sparse_type), allocatable :: hess_threads(:)
      !> Per-thread column buffers, summed deterministically after the region
      real(wp), allocatable :: hvp_threads(:, :, :, :)
      !> Column accumulator as it was handed in, kept only while a block can
      !> still fail after an earlier one has already been reduced into it
      real(wp), allocatable :: hvp_entry(:, :, :)
      !> Tangent of the folded surface adjoints along each direction of the block
      type(drop_surface_weights_type), allocatable :: deff(:)
      !> Thread bookkeeping
      integer :: thread_slot, ithread
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort
      !> Per-thread failure on its way to the latch
      type(error_type), allocatable :: worker_error

      !> Whether any channel writes gradient columns
      logical :: want_cols
      !> Prologue configuration: jet order and curvature request
      integer :: max_deriv
      logical :: want_curvature
      !> Direction block: its size, its bounds in `dirs` and its own extent
      integer :: ndir, nchunk, ilo, ihi, nblk
      !> Which fixed channel this block runs, if any; see the module header
      logical :: fixed_rank4_here, fixed_per_dir_here, fixed_here

      !> Per-thread state of the grid point being opened
      type(drop_point_scratch_type) :: pt
      !> Whether the shared prologue cleared the point for this traversal
      logical :: point_ok
      !> The 16 basis seeds of the point and their responses
      type(drop_point_seeds_type) :: seeds
      !> The fixed adjoints contracted at the point, and one response weight set
      type(drop_point_weights_type) :: pw, rw
      !> Per-direction temporaries of the second-order chain
      type(drop_fixed_dir_scratch_type) :: sc
      !> Jet tensors and scratch of the field contraction
      type(drop_field_tangent_work_type) :: ft_work
      !> Rank-4 direction set and pair table
      type(drop_rank4_scratch_type) :: r4
      !> The anchor's iSwiG neighbourhood
      type(drop_swi_scratch_type) :: swi

      !> Jet tangent along the current direction at the frozen point
      real(wp) :: dv0, dv1(3), dv2(3, 3), dv3(3, 3, 3)
      !> Gradient column of the current direction: field row and owner entries
      real(wp), allocatable :: field_row(:, :)
      real(wp) :: anchor_row(3)
      !> A direction restricted to the point's active slots
      real(wp), allocatable :: v_act(:, :)
      !> Effective position adjoint of the response channel; identically zero
      real(wp) :: w_xyz_zero(3)

      !> Grid, atom, axis, seed, slot and direction indices
      integer :: igrid, ient, i, iaxis
      integer :: idir, dir_atom, dir_axis, dir_loc, ia, ib
      integer :: jdir, idir_g, iatom, jj, knb

      !> Timer handle
      integer :: h_hess

      !* --------------------------- Channel preconditions ---------------------------- *!
      select case (fixed_mode)
      case (drop_fixed_none, drop_fixed_rank4, drop_fixed_per_dir)
      case default
         call fatal_error(error, context//": unknown fixed-channel mode")
         return
      end select
      if (fixed_mode == drop_fixed_none .and. .not. want_response) return
      if (fixed_mode == drop_fixed_rank4 .and. .not. present(hess_fixed)) then
         call fatal_error(error, context//": the rank-4 fixed Hessian channel needs an"// &
                          " accumulator")
         return
      end if
      if (fixed_mode == drop_fixed_per_dir .and. .not. (present(dirs) .and. present(hvp))) then
         call fatal_error(error, context//": the per-direction fixed Hessian channel"// &
                          " needs a direction set and an accumulator")
         return
      end if
      if (want_response .and. .not. (present(acc) .and. present(dirs) .and. present(hvp))) then
         call fatal_error(error, context//": the response Hessian channel needs the raw"// &
                          " adjoints, a direction set and an accumulator")
         return
      end if
      if (self%ngrid <= 0) return

      want_cols = want_response .or. fixed_mode == drop_fixed_per_dir

      if (fixed_mode == drop_fixed_per_dir) then
         if (want_response) then
            h_hess = self%ctx%timer%resolve("Surface Hessian (per direction)", &
                                            self%ctx%timer%current(), cat_gradient)
         else
            h_hess = self%ctx%timer%resolve("Surface Hessian (fixed adjoint, per direction)", &
                                            self%ctx%timer%current(), cat_gradient)
         end if
      else if (fixed_mode == drop_fixed_rank4) then
         if (want_response) then
            h_hess = self%ctx%timer%resolve("Surface Hessian (both halves)", &
                                            self%ctx%timer%current(), cat_gradient)
         else
            h_hess = self%ctx%timer%resolve("Surface Hessian (fixed adjoint)", &
                                            self%ctx%timer%current(), cat_gradient)
         end if
      else
         h_hess = self%ctx%timer%resolve("Surface Hessian (adjoint response)", &
                                         self%ctx%timer%current(), cat_gradient)
      end if
      call self%ctx%timer%start(h_hess)

      !* -------------------------------- Thread setup -------------------------------- *!
      ! Order 4 rather than the gradient path's 3 whenever a fixed channel is
      ! on: its field tangent reads `f4_rrrr` and `f4_rrr_rA`, and CFC asks for
      ! the highest *total* order. The response channel reads nothing above the
      ! third derivative, so a response-only traversal keeps order 3 and asks
      ! for no curvature. Both settings hold for every block of the traversal;
      ! see the module header for why a later block is not prepared cheaper.
      if (fixed_mode /= drop_fixed_none) then
         max_deriv = 4
         want_curvature = eff%have_wk
      else
         max_deriv = 3
         want_curvature = .false.
      end if
      call slots%init(self%ctx, self%lsf_model, max_deriv, self%param, self%mol, self%radii)

      if (fixed_mode == drop_fixed_rank4) allocate (hess_threads(slots%nthreads))

      !* ------------------------------ Direction blocks ------------------------------ *!
      ! See `drop_hvp_chunk_dirs` for the memory bound and the recompute factor.
      ! A direction set no larger than one block -- the common case -- takes a
      ! single iteration and is exactly the unblocked traversal. A rank-4-only
      ! traversal is a single block with no directions in it at all.
      !
      ! `hvp` is `intent(inout)` and this routine owes the caller an untouched
      ! accumulator when anything fails, which a multi-block run can no longer
      ! promise by construction: the reduction lands in `hvp` block by block, so
      ! that a direction's column sees the same additions in the same order
      ! whatever the blocking. A copy taken up front is what restores the
      ! promise, and it is taken only when a second block can actually fail.
      if (want_cols) then
         ndir = size(dirs, 3)
         nchunk = min(ndir, drop_hvp_chunk_dirs)
         allocate (hvp_threads(ndim, self%nsph, nchunk, slots%nthreads))
         if (nchunk < ndir) allocate (hvp_entry, source=hvp)
      else
         ndir = 1
         nchunk = 1
      end if

      do ilo = 1, ndir, nchunk
         ihi = min(ilo + nchunk - 1, ndir)
         nblk = ihi - ilo + 1

         ! The rank-4 channel is direction free and belongs to exactly one
         ! block; the per-direction channel belongs to every block. Every
         ! fixed-only step of the grid loop reads these rather than the mode.
         fixed_rank4_here = fixed_mode == drop_fixed_rank4 .and. ilo == 1
         fixed_per_dir_here = fixed_mode == drop_fixed_per_dir
         fixed_here = fixed_rank4_here .or. fixed_per_dir_here

         if (want_response) then
            !* ------------------ Passes 1 and 2: the moving weights ------------------ *!
            ! Serial over this block's directions, for the group-reduction
            ! reason the module header documents; `deff` is indexed `1 .. nblk`.
            call weight_tangents(self, acc, eff, dirs(:, :, ilo:ihi), deff, error)
            if (allocated(error)) then
               if (allocated(hvp_entry)) hvp = hvp_entry
               call self%ctx%timer%stop(h_hess)
               return
            end if
         end if
         if (want_cols) hvp_threads = 0.0_wp

         call abort%reset()

         !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, &
         !$omp& igrid, pt, point_ok, seeds, pw, rw, sc, ft_work, r4, swi, &
         !$omp& dv0, dv1, dv2, dv3, field_row, anchor_row, v_act, w_xyz_zero, &
         !$omp& ient, i, iaxis, idir, dir_atom, dir_axis, dir_loc, ia, ib, &
         !$omp& jdir, idir_g, iatom, jj, knb, worker_error)
         thread_slot = 1
!$       thread_slot = omp_get_thread_num() + 1

         ! `want_lsf4` asks the prologue for the fourth-derivative buffer, read
         ! by the fixed channel alone: a later block that allocated it would
         ! have the prologue evaluate a fourth-order kernel per point for
         ! nothing, and `lsf4_rrrr` is a pure output, so leaving it out changes
         ! no other value. `want_vjp` is the response channel's row buffer.
         call pt%init(self%nsph, self%iswig, want_seed_batch=.true., &
                      want_lsf4=fixed_here, want_vjp=want_response)
         ! The response channel's contractions see a hard zero here; see the
         ! module header for why its normal fold is not merely small
         w_xyz_zero = 0.0_wp

         if (fixed_here) then
            allocate (field_row(3, self%nsph), source=0.0_wp)
            ! Grown on demand from `n_nb + 1` rather than sized to `nsph`; the
            ! block is quadratic in the influence set
            allocate (swi%blk(3, 1, 3, 1), swi%idx(1))
         end if
         if (fixed_rank4_here) then
            call hess_sparse_new(hess_threads(thread_slot), self%nsph)
            allocate (r4%dir_atoms(self%nsph))
            ! Grown from the direction set it indexes; both start at one slot
            ! and are grown by the first point
            allocate (r4%pair_ent(1, 1), swi%ent(1, 1))
            allocate (r4%vdir(3, self%nsph), source=0.0_wp)
         end if
         if (fixed_per_dir_here) then
            allocate (v_act(3, self%nsph), source=0.0_wp)
         end if

         !$omp do schedule(static, 8)
         do igrid = 1, self%ngrid
            !* ------------------------ Shared point prologue ------------------------- *!
            ! Point, jets, seed state and the solved jet and anchor seeds. Every
            ! failure path -- a latch already set, a refusing level set, a
            ! degenerate state, a singular bordered system -- has recorded
            ! itself. `pt%lsf4_rrrr` comes back filled when a fixed channel
            ! asked for it.
            call drop_point_prologue(self, slots, thread_slot, igrid, want_curvature, &
                                     context, abort, pt, point_ok)
            if (.not. point_ok) cycle

            !* ------------------- Basis seeds and their responses -------------------- *!
            ! The 16 seeds of [[seed_jet_basis]] and [[seed_anchor]], applied
            ! once for both channels. `fill_seed_basis` reads the same solved
            ! batch through the same [[seed_rhs_column]] map, so it writes
            ! `seeds%x` with the values the applies below write again; what
            ! only it provides are the seed directions the second-order chain
            ! needs.
            if (fixed_here) then
               call fill_seed_basis(pt%kkt_rhs, seeds%dlsf1, seeds%dlsf2, seeds%x)
            end if
            call point_seed_basis(pt%state, pt%kkt_rhs, fixed_here, seeds%res, seeds%x, &
                                  seeds%dstate)

            ! The active atoms of the level set at this point. Both channels
            ! index rows by them, and every direction is gathered onto them.
            pt%n_active = slots%lsf(thread_slot)%lsf%active_count()
            do i = 1, pt%n_active
               pt%active_idx(i) = slots%lsf(thread_slot)%lsf%active_atom(i)
            end do

            !* ----------------------- Jet tensors of the point ----------------------- *!
            ! Direction free, so formed once per point and read by every
            ! contraction below. The fills are unconditional on purpose: the
            ! buffers outlive the point, and a skipped fill would serve the
            ! previous point's tensors. The fixed channel also needs the mixed
            ! fourth derivative; the response channel does not.
            if (fixed_here) then
               call drop_field_tangent_point(slots%lsf(thread_slot)%lsf, ft_work)
            else
               call drop_field_jet_point(slots%lsf(thread_slot)%lsf, ft_work)
            end if

            !* --------------- Shared iSwiG neighbourhood of the anchor --------------- *!
            ! The switching factor's neighbour set depends on the geometry
            ! alone, so one `swi_collect` serves both channels; only what is
            ! read out of it differs. The fixed channel wants the whole
            ! second-derivative block when its own weight is live, the response
            ! channel the first-derivative rows when any direction of the block
            ! carries a live switching tangent. The width outputs of the block
            ! are not asked for: the switching factor is evaluated at the
            ! *anchor* width, which is nuclear-geometry independent.
            swi%fixed = .false.
            if (fixed_here) swi%fixed = abs(eff%w_f(igrid)) > seed_weight_tol
            swi%resp = .false.
            if (want_response) then
               do jdir = 1, nblk
                  if (abs(deff(jdir)%w_f(igrid)) > seed_weight_tol) then
                     swi%resp = .true.
                     exit
                  end if
               end do
            end if
            if (swi%fixed .or. swi%resp) then
               call self%iswig%swi_collect(pt%anchor, pt%owner_idx, self%anchor_xi0(igrid), &
                                           swi%f0, pt%iswig_work)
            end if
            if (swi%resp) then
               call self%iswig%swi1_rA_sparse(pt%iswig_work, pt%swi_rows, &
                                              swi%owner_row, swi%dxi)
            end if
            if (swi%fixed) then
               if (size(swi%idx) < pt%iswig_work%n_nb + 1) then
                  deallocate (swi%blk, swi%idx)
                  allocate (swi%blk(3, pt%iswig_work%n_nb + 1, 3, pt%iswig_work%n_nb + 1))
                  allocate (swi%idx(pt%iswig_work%n_nb + 1))
                  if (fixed_rank4_here) then
                     deallocate (swi%ent)
                     allocate (swi%ent(pt%iswig_work%n_nb + 1, pt%iswig_work%n_nb + 1))
                  end if
               end if
               call self%iswig%swi2_rArB_block(pt%iswig_work, swi%n, swi%idx, swi%blk)
            end if

            !* ======================== Fixed-adjoint channel ========================= *!
            if (fixed_here) then
               call point_weights(pt, eff, igrid, seeds, pw)
               ! The Hessian weights fold into the mixed fourth derivative once
               ! per point; every direction then reads three numbers per slot
               call drop_field_f4_fold(ft_work, pt%n_active, pw%w_lsf2)
            end if

            if (fixed_rank4_here) then
               !* ------------------------ Local direction set ------------------------ *!
               ! The projected point, the jet and the anchor depend on the
               ! nuclei only through the level set's active atoms and through
               ! the owner, so every other column of this point's block is
               ! exactly zero. The active atoms come first, in slot order, so
               ! that a direction's local index below `n_active` *is* its slot.
               r4%ndir_atom = pt%n_active
               r4%dir_atoms(1:pt%n_active) = pt%active_idx(1:pt%n_active)
               if (.not. any(r4%dir_atoms(1:r4%ndir_atom) == pt%owner_idx)) then
                  r4%ndir_atom = r4%ndir_atom + 1
                  r4%dir_atoms(r4%ndir_atom) = pt%owner_idx
               end if
               r4%owner_loc = r4%ndir_atom
               do i = 1, r4%ndir_atom
                  if (r4%dir_atoms(i) == pt%owner_idx) then
                     r4%owner_loc = i
                     exit
                  end if
               end do

               ! Rows are the active atoms and the owner, columns the same set,
               ! so the whole square is resolved here -- once per point rather
               ! than once per direction. Creating an entry commits nothing: an
               ! untouched one holds an exact zero.
               if (size(r4%pair_ent, 1) < r4%ndir_atom) then
                  deallocate (r4%pair_ent)
                  allocate (r4%pair_ent(r4%ndir_atom, r4%ndir_atom))
               end if
               do ib = 1, r4%ndir_atom
                  do ia = 1, r4%ndir_atom
                     call hess_sparse_entry(hess_threads(thread_slot), r4%dir_atoms(ia), &
                                            r4%dir_atoms(ib), r4%pair_ent(ia, ib))
                  end do
               end do

               !* --------------- Explicit nuclear motion of every column -------------- *!
               ! The weighted mixed nuclear Hessian block of the level set over
               ! the active slots: the explicit half of every column of this
               ! point in one accessor call, where the chain would call
               ! `hvp_jet_rA` once per column. Sized to the point rather than
               ! grown, so the level set's product lands in it without a copy.
               if (pt%n_active > 0) then
                  if (allocated(r4%mblk)) then
                     if (size(r4%mblk, 1) /= 3*pt%n_active) deallocate (r4%mblk)
                  end if
                  if (.not. allocated(r4%mblk)) then
                     allocate (r4%mblk(3*pt%n_active, 3*pt%n_active))
                  end if
                  call slots%lsf(thread_slot)%lsf%vjp_f2_rArB(pw%w_lsf0, pw%w_lsf1, &
                                                                 pw%w_lsf2, r4%mblk)
               end if

               !* ------------------- One chain per basis direction ------------------- *!
               do idir = 1, 3*r4%ndir_atom
                  dir_loc = (idir - 1)/3 + 1
                  dir_atom = r4%dir_atoms(dir_loc)
                  dir_axis = mod(idir - 1, 3) + 1
                  r4%vdir(dir_axis, dir_atom) = 1.0_wp

                  ! The jet tangent along a unit direction is a column of the
                  ! point's tensors when the atom is active. An owner outside
                  ! the active set does not enter the level set, so its column
                  ! is zero and only the anchor moves.
                  if (dir_loc <= pt%n_active) then
                     call drop_field_jet_column(ft_work, pt%n_active, dir_axis, dir_loc, &
                                                dv0, dv1, dv2, dv3)
                  else
                     dv0 = 0.0_wp
                     dv1 = 0.0_wp
                     dv2 = 0.0_wp
                     dv3 = 0.0_wp
                  end if

                  call fixed_direction_chain(self, slots%lsf(thread_slot)%lsf, pt, eff, &
                                             igrid, seeds, pw, r4%vdir, dv0, dv1, dv2, dv3, &
                                             .false., context, sc, ft_work, field_row, &
                                             anchor_row, worker_error)
                  r4%vdir(dir_axis, dir_atom) = 0.0_wp
                  if (allocated(worker_error)) then
                     call abort%latch_error(worker_error, igrid)
                     exit
                  end if

                  ! The explicit nuclear motion of this column, read off the
                  ! point's block; an owner outside the active set has none
                  if (dir_loc <= pt%n_active) then
                     do i = 1, pt%n_active
                        field_row(:, i) = field_row(:, i) &
                           + r4%mblk(3*(i - 1) + 1:3*i, 3*(dir_loc - 1) + dir_axis)
                     end do
                  end if

                  ! The anchor entries belong to the owner's row; the field row
                  ! to the active atoms' rows. Both land in column
                  ! `(dir_axis, dir_atom)` of the block.
                  ient = r4%pair_ent(r4%owner_loc, dir_loc)
                  do iaxis = 1, 3
                     hess_threads(thread_slot)%blocks(iaxis, dir_axis, ient) = &
                        hess_threads(thread_slot)%blocks(iaxis, dir_axis, ient) &
                        + anchor_row(iaxis)
                  end do
                  do i = 1, pt%n_active
                     ient = r4%pair_ent(i, dir_loc)
                     hess_threads(thread_slot)%blocks(:, dir_axis, ient) = &
                        hess_threads(thread_slot)%blocks(:, dir_axis, ient) &
                        + field_row(:, i)
                  end do
               end do
               ! A latched direction abandons the whole point, the response
               ! channel below included: the run is failing and nothing will be
               ! reduced out of it.
               if (abort%requested) cycle

               !* ---------------------- iSwiG switching channel ---------------------- *!
               ! `f_i` depends on the nuclear geometry alone and its adjoint is
               ! fixed, so the whole channel is one block over the influence
               ! set, with no loop over directions at all. The influence set is
               ! not the direction set, so its pairs are resolved on their own;
               ! any pair the two share resolves to the one entry the field
               ! channel already wrote, and lands after it.
               if (swi%fixed) then
                  do ib = 1, swi%n
                     do ia = 1, swi%n
                        call hess_sparse_entry(hess_threads(thread_slot), swi%idx(ia), &
                                               swi%idx(ib), swi%ent(ia, ib))
                     end do
                  end do
                  call scatter_iswig_block_indexed(swi%n, swi%ent, swi%blk, eff%w_f(igrid), &
                                                   hess_threads(thread_slot)%blocks)
               end if
            end if

            if (fixed_per_dir_here) then
               !* ----------------- One chain per supplied direction ------------------ *!
               do jdir = 1, nblk
                  idir_g = ilo + jdir - 1

                  ! The jet tangent along the direction: the point's tensors
                  ! contracted with the direction gathered onto the active
                  ! slots. Components on atoms outside the active set do not
                  ! enter the level set; the owner's enters through the chain.
                  do i = 1, pt%n_active
                     v_act(:, i) = dirs(:, pt%active_idx(i), idir_g)
                  end do
                  call drop_field_jet_tangent(ft_work, pt%n_active, v_act, dv0, dv1, dv2, dv3)

                  call fixed_direction_chain(self, slots%lsf(thread_slot)%lsf, pt, eff, &
                                             igrid, seeds, pw, dirs(:, :, idir_g), &
                                             dv0, dv1, dv2, dv3, .true., context, sc, &
                                             ft_work, field_row, anchor_row, worker_error)
                  if (allocated(worker_error)) then
                     call abort%latch_error(worker_error, igrid)
                     exit
                  end if

                  do i = 1, pt%n_active
                     iatom = pt%active_idx(i)
                     hvp_threads(:, iatom, jdir, thread_slot) = &
                        hvp_threads(:, iatom, jdir, thread_slot) + field_row(:, i)
                  end do
                  hvp_threads(:, pt%owner_idx, jdir, thread_slot) = &
                     hvp_threads(:, pt%owner_idx, jdir, thread_slot) + anchor_row

                  ! The switching block of the point contracted with the
                  ! direction over the influence set
                  if (swi%fixed) then
                     call contract_iswig_block(swi%n, swi%idx, swi%blk, eff%w_f(igrid), &
                                               dirs(:, :, idir_g), &
                                               hvp_threads(:, :, jdir, thread_slot))
                  end if
               end do
               if (abort%requested) cycle
            end if

            !* ======================= Adjoint-response channel ======================= *!
            if (want_response) then
               do jdir = 1, nblk

                  !* ----------- Field seeds -> level-set adjoint tangents ------------ *!
                  rw%w_lsf0 = 0.0_wp
                  rw%w_lsf1 = 0.0_wp
                  rw%w_lsf2 = 0.0_wp
                  call seed_jet_basis_contract(deff(jdir), igrid, pt%phi1_r, w_xyz_zero, &
                                               seeds%res(1:drop_n_jet_seeds), &
                                               seeds%x(:, 1:drop_n_jet_seeds), &
                                               rw%w_lsf0, rw%w_lsf1, rw%w_lsf2)

                  !* ------------- Field channel: the weighted row, read off ---------- *!
                  ! the jet tensors formed once above; no LSF call per direction
                  call drop_field_jet_contract(ft_work, pt%n_active, rw%w_lsf0, rw%w_lsf1, &
                                               rw%w_lsf2, pt%vjp_pt)
                  do i = 1, pt%n_active
                     iatom = pt%active_idx(i)
                     hvp_threads(:, iatom, jdir, thread_slot) = &
                        hvp_threads(:, iatom, jdir, thread_slot) + pt%vjp_pt(:, i)
                  end do

                  !* --------------------- Anchor channel (owner) --------------------- *!
                  call seed_anchor_contract(deff(jdir), igrid, pt%phi1_r, w_xyz_zero, &
                                            seeds%res(drop_n_jet_seeds + 1:), &
                                            seeds%x(:, drop_n_jet_seeds + 1:), &
                                            hvp_threads(:, pt%owner_idx, jdir, thread_slot))

                  !* -------------------- iSwiG switching channel --------------------- *!
                  if (abs(deff(jdir)%w_f(igrid)) > seed_weight_tol) then
                     do jj = 1, pt%iswig_work%n_nb
                        knb = pt%iswig_work%idx(jj)
                        hvp_threads(:, knb, jdir, thread_slot) = &
                           hvp_threads(:, knb, jdir, thread_slot) &
                           + deff(jdir)%w_f(igrid)*pt%swi_rows(:, jj)
                     end do
                     hvp_threads(:, pt%owner_idx, jdir, thread_slot) = &
                        hvp_threads(:, pt%owner_idx, jdir, thread_slot) &
                        + deff(jdir)%w_f(igrid)*swi%owner_row
                  end if

               end do
            end if

         end do
         !$omp end do

         if (fixed_here) deallocate (field_row, swi%blk, swi%idx)
         if (fixed_rank4_here) then
            deallocate (r4%dir_atoms, r4%pair_ent, r4%vdir, swi%ent)
            if (allocated(r4%mblk)) deallocate (r4%mblk)
         end if
         if (fixed_per_dir_here) deallocate (v_act)
         call pt%destroy()
         !$omp end parallel

         if (abort%requested) then
            call abort%raise(context, error)
            if (allocated(hvp_entry)) hvp = hvp_entry
            call self%ctx%timer%stop(h_hess)
            return
         end if

         ! Deterministic reductions: fixed thread order, independent of
         ! scheduling. The column channels take one block's columns at a time,
         ! and blocks are disjoint in the direction index, so a direction sees
         ! exactly the same additions in the same order however the blocking
         ! falls.
         if (fixed_rank4_here) then
            do ithread = 1, slots%nthreads
               call hess_sparse_reduce(hess_threads(ithread), hess_fixed)
               call hess_sparse_destroy(hess_threads(ithread))
            end do
         end if
         if (want_cols) then
            do ithread = 1, slots%nthreads
               hvp(:, :, ilo:ihi) = hvp(:, :, ilo:ihi) + hvp_threads(:, :, 1:nblk, ithread)
            end do
         end if

      end do

      call self%ctx%timer%stop(h_hess)

   end subroutine drop_hessian_traverse

   !* ================================================================================= *!
   !*                     Per-point pieces of the fixed channel                         *!
   !* ================================================================================= *!

   !> Apply the 16 basis seeds of one grid point, with or without their tangents
   !>
   !> The single seed-application site of the traversal. The jet and anchor
   !> halves are the shipped ones, written into one `drop_n_point_seeds` array
   !> in the order [[fill_seed_basis]] lays out, so the second-order chain and
   !> the response contractions index the same object.
   !>
   !> `want_dstate` rather than an optional argument on the caller's side: the
   !> derived-state tangents are the fixed channel's alone, and a
   !> response-only traversal must not pay for them.
   !>
   !> @param[in]    state       Per-grid point forward state
   !> @param[in]    kkt         Solved KKT sensitivities of the point
   !> @param[in]    want_dstate Whether the derived-state tangents are needed
   !> @param[out]   res_seed    Linear response of each seed
   !> @param[out]   seed_x      Induced point motion and multiplier change of each seed
   !> @param[inout] dstate_seed Tangent of the derived state, when asked for
   subroutine point_seed_basis(state, kkt, want_dstate, res_seed, seed_x, dstate_seed)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Solved KKT sensitivities
      real(wp), intent(in) :: kkt(:, :)
      !> Whether the derived-state tangents are needed
      logical, intent(in) :: want_dstate
      !> Linear response of each seed
      type(drop_seed_result_type), intent(out) :: res_seed(drop_n_point_seeds)
      !> Induced point motion and multiplier change of each seed
      real(wp), intent(out) :: seed_x(4, drop_n_point_seeds)
      !> Tangent of the derived seed state along each seed
      type(drop_seed_state_tangent_type), intent(inout) :: dstate_seed(drop_n_point_seeds)

      if (want_dstate) then
         call seed_jet_basis_apply(state, kkt, res_seed(1:drop_n_jet_seeds), &
                                   seed_x(:, 1:drop_n_jet_seeds), &
                                   dstate_seed(1:drop_n_jet_seeds))
         call seed_anchor_apply(state, kkt, res_seed(drop_n_jet_seeds + 1:), &
                                seed_x(:, drop_n_jet_seeds + 1:), &
                                dstate_seed(drop_n_jet_seeds + 1:))
      else
         call seed_jet_basis_apply(state, kkt, res_seed(1:drop_n_jet_seeds), &
                                   seed_x(:, 1:drop_n_jet_seeds))
         call seed_anchor_apply(state, kkt, res_seed(drop_n_jet_seeds + 1:), &
                                seed_x(:, drop_n_jet_seeds + 1:))
      end if

   end subroutine point_seed_basis

   !> Contract the fixed adjoints of one grid point into its weight set
   !>
   !> The outward-normal channel, as on the gradient path, and the 13 jet-seed
   !> contributions scattered into the level-set adjoint weights. The chain
   !> below needs `normal_grad` and `nwn` again to build their tangents, so they
   !> are asked of the fold rather than recomputed: a second copy of one
   !> floating-point chain is free to contract differently, and
   !> [[seed_normal_channel]] owns this one. Both come back zero when the
   !> channel is inactive.
   !>
   !> The bound of the seed loop is a *jet* one: only the 13 jet slots have a
   !> weight to land in. An anchor seed's contraction belongs to the owner's
   !> gradient row instead, which the chain handles -- and is why
   !> [[scatter_jet_weight]] deliberately has no anchor case.
   !>
   !> @param[in]  pt    Point scratch, prologue run
   !> @param[in]  eff   Folded surface adjoints, held fixed
   !> @param[in]  igrid Grid point
   !> @param[in]  seeds Applied basis seeds of the point
   !> @param[out] pw    Weight set of the point
   pure subroutine point_weights(pt, eff, igrid, seeds, pw)
      !> Point scratch
      type(drop_point_scratch_type), intent(in) :: pt
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Applied basis seeds
      type(drop_point_seeds_type), intent(in) :: seeds
      !> Weight set of the point
      type(drop_point_weights_type), intent(out) :: pw

      !> One seed's contribution
      real(wp) :: contribution
      !> Seed index
      integer :: ibasis

      pw%w_lsf0 = 0.0_wp
      pw%w_lsf1 = 0.0_wp
      pw%w_lsf2 = 0.0_wp
      call seed_normal_channel(pt%state, eff, igrid, pt%lsf2_rr, pw%w_lsf1, pw%w_xyz, &
                               normal_grad_pt=pw%normal_grad, nwn_pt=pw%nwn)

      do ibasis = 1, drop_n_jet_seeds
         contribution = seed_contribution(eff, igrid, pw%w_xyz, seeds%x(1:3, ibasis), &
                                          seeds%res(ibasis), pt%phi1_r)
         call scatter_jet_weight(ibasis, contribution, pw%w_lsf0, pw%w_lsf1, pw%w_lsf2)
      end do
   end subroutine point_weights

   !> The second-order chain of the fixed channel along one nuclear direction
   !>
   !> Given the jet tangent `(dv0, dv1, dv2, dv3)` of the level set along `v` at
   !> the *frozen* projected point, returns the gradient column of
   !> `(dJ^T/dv) omega` at this point: the field row over the active slots
   !> (`field_row(:, 1:n_active)`, overwritten) and the three anchor-seed
   !> entries of the owner (`anchor_row`). The switching channel is not part
   !> of the chain, because it needs no seed; the caller contracts it.
   !>
   !> Everything the chain reads comes in through its arguments, and that is a
   !> constraint rather than a style: it is called from inside an OpenMP
   !> region, and a procedure that reached the traversal's variables by host
   !> association would read the shared originals of what the threads hold
   !> privately.
   !>
   !> A nuclear direction is just another seed of the same linear map, so
   !> [[apply_seed]] produces the whole forward tangent of the point; the
   !> input tangents are *total* (`state%lsf1_r` is grad S at the projected
   !> point, whose `v`-tangent carries the point motion), which is what
   !> `res_v%dg` and `res_v%dH` already are, and the third spatial derivative
   !> is folded by hand. The seed right-hand sides move as `K dx = db - dK x`
   !> on the same factors, where only the three gradient seeds carry a `db`.
   !> Every basis seed is a constant matrix, so its own tangent vanishes and
   !> only the induced point motion moves.
   !>
   !> `dw_lsf2` is read only inside [[drop_field_tangent_dir]], against
   !> `lsf3_rr_rA`; see the module header for why that is what makes the nine
   !> asymmetric Hessian seeds legitimate. Never read a single off-diagonal
   !> entry of it on its own.
   !>
   !> @param[in]    self         DROP cavity instance, for the objective parameter
   !> @param[in]    lsf          This thread's level set, prepared at the point
   !> @param[in]    pt           Point scratch, prologue run
   !> @param[in]    eff          Folded surface adjoints, held fixed
   !> @param[in]    igrid        Grid point
   !> @param[in]    seeds        Applied basis seeds of the point, with tangents
   !> @param[in]    pw           Weight set of the point
   !> @param[in]    v            Nuclear direction `(3, nsph)`
   !> @param[in]    dv0          Jet tangent of the value along `v` at the frozen point
   !> @param[in]    dv1          Jet tangent of the gradient
   !> @param[in]    dv2          Jet tangent of the Hessian
   !> @param[in]    dv3          Jet tangent of the third derivative
   !> @param[in]    explicit     Include the explicit nuclear motion of the field row through
   !>                            `hvp_jet_rA`; `.false.` when the caller adds it from a block
   !> @param[in]    context      Calling routine, used to prefix the diagnostics
   !> @param[inout] sc           Per-direction temporaries
   !> @param[inout] ft_work      Jet tensors of the point, fills and fold run
   !> @param[inout] field_row    Field row of the column `(3, >= n_active)`
   !> @param[out]   anchor_row   Owner entries of the column
   !> @param[out]   worker_error Failure of the bordered solves, if any
   subroutine fixed_direction_chain(self, lsf, pt, eff, igrid, seeds, pw, v, &
                                    dv0, dv1, dv2, dv3, explicit, context, sc, ft_work, &
                                    field_row, anchor_row, worker_error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> This thread's level set, prepared at the point
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Point scratch
      type(drop_point_scratch_type), intent(in) :: pt
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Applied basis seeds of the point
      type(drop_point_seeds_type), intent(in) :: seeds
      !> Weight set of the point
      type(drop_point_weights_type), intent(in) :: pw
      !> Nuclear direction
      real(wp), intent(in) :: v(:, :)
      !> Jet tangent along `v` at the frozen point
      real(wp), intent(in) :: dv0, dv1(3), dv2(3, 3), dv3(3, 3, 3)
      !> Whether the field row includes the explicit nuclear motion
      logical, intent(in) :: explicit
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Per-direction temporaries
      type(drop_fixed_dir_scratch_type), intent(inout) :: sc
      !> Jet tensors of the point
      type(drop_field_tangent_work_type), intent(inout) :: ft_work
      !> Field row of the column
      real(wp), intent(inout) :: field_row(:, :)
      !> Owner entries of the column
      real(wp), intent(out) :: anchor_row(3)
      !> Failure of the bordered solves
      type(error_type), allocatable, intent(out) :: worker_error

      !> One seed's contribution
      real(wp) :: contribution
      !> Cartesian and seed indices
      integer :: iaxis, ibasis, k

      anchor_row = 0.0_wp

      !* ------------- Directional response of the projected point -------------- *!
      ! The projection rides on the jet tangent through the same bordered
      ! system the seeds use, with the anchor moving rigidly with its owner:
      ! `d^2 phi/(dr dR_owner) = -alpha I`.
      sc%rhs_v = 0.0_wp
      sc%rhs_v(1:3, 1) = pt%lambda_val*dv1 + self%param%phi_alpha*v(:, pt%owner_idx)
      sc%rhs_v(4, 1) = -dv0
      call pt%kkt_fac%solve(sc%rhs_v, context, worker_error, igrid)
      if (allocated(worker_error)) return
      sc%dr_v = sc%rhs_v(1:3, 1)
      sc%dl_v = sc%rhs_v(4, 1)

      !* ---------------- Directional state and input tangents ------------------ *!
      call apply_seed(pt%state, dv1, dv2, sc%dr_v, sc%dl_v, sc%res_v, sc%dstate_v)

      sc%dinp_v%dlsf1_r = sc%res_v%dg
      sc%dinp_v%dlsf2_rr = sc%res_v%dH
      sc%dinp_v%dlsf3_rrr = dv3
      do k = 1, 3
         sc%dinp_v%dlsf3_rrr(:, :, :) = sc%dinp_v%dlsf3_rrr(:, :, :) &
                                        + pt%lsf4_rrrr(:, :, :, k)*sc%dr_v(k)
      end do
      sc%dinp_v%dlambda_val = sc%dl_v
      sc%dinp_v%dcpjac_scal0 = sc%res_v%dJ
      sc%dinp_v%dw_f0 = sc%res_v%dw_f
      sc%dinp_v%dwleb = sc%res_v%dwleb
      sc%dinp_v%dxi0 = sc%res_v%dxi
      ! The anchor's Lebedev weight is a property of the rigid sphere and the
      ! branch weight is one for every group this channel admits, so both
      ! tangents vanish; see the scope limit in `hessian.f90`.
      sc%dinp_v%danchor_wleb0 = 0.0_wp
      sc%dinp_v%dwbranch = 0.0_wp

      !* --------------------- Tangent of the normal fold ----------------------- *!
      sc%dnormal_grad = 0.0_wp
      sc%dw_xyz = 0.0_wp
      if (eff%have_wn) then
         sc%dnormal_grad = (-sc%res_v%dn_surf*pw%nwn &
                            - pt%state%n_surf*dot_product(sc%res_v%dn_surf, eff%w_n(:, igrid))) &
                           /pt%state%g_norm &
                           - pw%normal_grad*sc%res_v%d_gnorm/pt%state%g_norm
         sc%dw_xyz = matmul(sc%res_v%dH, pw%normal_grad) &
                     + matmul(pt%lsf2_rr, sc%dnormal_grad)
      end if

      !* ---------------- Tangent of the seed right-hand sides ------------------ *!
      sc%dH_lag(:, :, 1) = -sc%dl_v*pt%lsf2_rr - pt%lambda_val*sc%res_v%dH
      sc%dg_tot(:, 1) = sc%res_v%dg
      sc%dseed_x = 0.0_wp
      do iaxis = 1, 3
         sc%dseed_x(iaxis, 1 + iaxis) = sc%dl_v
      end do
      call pt%kkt_fac%solve_tangent(sc%dH_lag, sc%dg_tot, seeds%x, sc%dseed_x, &
                                    context, worker_error, igrid)
      if (allocated(worker_error)) return

      !* -------------------------- Weight tangents ----------------------------- *!
      sc%dw_lsf0 = 0.0_wp
      sc%dw_lsf1 = sc%dnormal_grad
      sc%dw_lsf2 = 0.0_wp
      do ibasis = 1, drop_n_point_seeds
         call apply_seed_tangent(pt%state, sc%dstate_v, sc%dinp_v, sc%res_v, &
                                 seeds%dlsf1(:, ibasis), seeds%dlsf2(:, :, ibasis), &
                                 seeds%x(1:3, ibasis), seeds%x(4, ibasis), &
                                 seed_dzero1, seed_dzero2, &
                                 sc%dseed_x(1:3, ibasis), sc%dseed_x(4, ibasis), &
                                 seeds%res(ibasis), seeds%dstate(ibasis), sc%dres)
         contribution = seed_contribution_tangent(eff, igrid, pw%w_xyz, sc%dw_xyz, &
                                                  seeds%x(1:3, ibasis), &
                                                  sc%dseed_x(1:3, ibasis), sc%dres)
         if (ibasis > drop_n_jet_seeds) then
            ! Anchor seed: the gradient row it feeds belongs to the owner
            anchor_row(ibasis - drop_n_jet_seeds) = contribution
         else
            call scatter_jet_weight(ibasis, contribution, sc%dw_lsf0, sc%dw_lsf1, sc%dw_lsf2)
         end if
      end do

      !* --------------------------- Field channel ------------------------------ *!
      if (pt%n_active > 0) then
         call drop_field_tangent_dir(lsf, pw%w_lsf0, pw%w_lsf1, pw%w_lsf2, &
                                     sc%dw_lsf0, sc%dw_lsf1, sc%dw_lsf2, sc%dr_v, v, &
                                     explicit, ft_work, field_row)
      end if

   end subroutine fixed_direction_chain

   !* ================================================================================= *!
   !*                       Passes 1 and 2: d(eff) per direction                        *!
   !* ================================================================================= *!

   !> Build the directional tangent of the folded surface weights
   !>
   !> Runs the forward tangent once for the whole batch and the weight tangent
   !> once per direction. The second loop is **serial on purpose**: it is the
   !> stage that carries [[branch_phi_adj_tangent]]'s group reduction, the only
   !> cross-point coupling in the scheme, and the primitive is serial over the
   !> whole grid, so every contiguous anchor group is seen whole by exactly one
   !> call. Directions are the admissible parallel axis if this ever needs one;
   !> grid points are not, and never will be.
   !>
   !> `deff(idir)` carries the three channels pass 2 emits and nothing else --
   !> see the module header for why the other five are absent rather than zero.
   !>
   !> `dirs` is one *block* of the caller's direction set and both extents are
   !> taken from it, so `deff` is indexed within the block and the seven
   !> `(ngrid, ndir)` arrays this routine holds are block sized. That is where
   !> the working set of the whole traversal is bounded; see
   !> `drop_hvp_chunk_dirs`.
   !>
   !> @param[in]  self  DROP cavity instance
   !> @param[in]  acc   Raw surface adjoints, held fixed
   !> @param[in]  eff   Folded weights, as [[prepare_surface_weights]] returned them
   !> @param[in]  dirs  Nuclear directions of one block, `(3, nsph, ndir)`
   !> @param[out] deff  Tangent of the folded weights, one element per direction
   !> @param[out] error Error object, allocated on failure
   subroutine weight_tangents(self, acc, eff, dirs, deff, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Raw surface adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Folded weights
      type(drop_surface_weights_type), intent(in) :: eff
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Tangent of the folded weights, one element per direction
      type(drop_surface_weights_type), allocatable, intent(out) :: deff(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Directional tangents of the four grid scalars pass 2 consumes
      real(wp), allocatable :: d_a(:, :), d_wleb(:, :), d_xi0(:, :), d_wbranch(:, :)
      !> Branch bookkeeping, defaulted when the cavity carries none
      integer, allocatable :: branch_count(:), anchor_id(:)
      !> Grid extent and direction index
      integer :: ngrid, ndir, idir

      ngrid = self%ngrid
      ndir = size(dirs, 3)

      !* -------------------------- Pass 1: forward tangent --------------------------- *!
      allocate (d_a(ngrid, ndir), d_wleb(ngrid, ndir), d_xi0(ngrid, ndir), &
                d_wbranch(ngrid, ndir))
      call self%get_surface_tangent(dirs, d_a, d_wleb, d_xi0, d_wbranch, error)
      if (allocated(error)) return

      ! `compute_branch_phi_adj` is skipped by the primal when the cavity holds
      ! no branch bookkeeping; a single-branch stand-in reproduces that early
      ! exit in the tangent instead of duplicating the guard.
      allocate (branch_count(ngrid), source=1)
      allocate (anchor_id(ngrid), source=0)
      if (allocated(self%branch_count)) branch_count = self%branch_count(1:ngrid)
      if (allocated(self%anchor_id)) anchor_id = self%anchor_id(1:ngrid)

      !* --------------------------- Pass 2: weight tangent --------------------------- *!
      allocate (deff(ndir))
      do idir = 1, ndir
         allocate (deff(idir)%w_xi(ngrid), deff(idir)%w_f(ngrid), &
                   deff(idir)%branch_phi_adj(ngrid))
         ! `have_wn` and `have_wk` keep their `.false.` default: the normal and
         ! curvature adjoints are copies of the fixed raw channels, so their
         ! tangent is identically zero.
         call prepare_surface_weights_tangent(acc, eff, .true., &
                                              self%a(1:ngrid), self%wleb(1:ngrid), &
                                              self%xi0(1:ngrid), self%wbranch(1:ngrid), &
                                              self%radii, self%owner(1:ngrid), &
                                              branch_count, anchor_id, &
                                              self%branch_weight%s, &
                                              d_a(:, idir), d_wleb(:, idir), &
                                              d_xi0(:, idir), d_wbranch(:, idir), &
                                              deff(idir)%w_xi, deff(idir)%w_f, &
                                              deff(idir)%branch_phi_adj, error)
         if (allocated(error)) return
      end do

   end subroutine weight_tangents

   !* ================================================================================= *!
   !*                       Sparse per-thread Hessian accumulator                       *!
   !* ================================================================================= *!

   !> Give one thread an empty accumulator
   !>
   !> Called from inside the parallel region on the thread that will own it, so
   !> the storage is first touched by the thread that writes it.
   !>
   !> @param[inout] self Accumulator of one thread
   !> @param[in]    nsph Atoms the pair indices range over
   subroutine hess_sparse_new(self, nsph)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self
      !> Atoms the pair indices range over
      integer, intent(in) :: nsph

      integer :: cap, ntab

      self%nent = 0
      self%maxent = int(min(int(nsph, int64)**2, int(huge(self%maxent), int64)))

      cap = max(1, min(hess_init_ent, self%maxent))
      allocate (self%pair_i(cap), self%pair_j(cap))
      allocate (self%blocks(ndim, ndim, cap))

      ! Halve the default table while the one below it would still hold `cap`
      ! entries under the load factor, so a molecule with a handful of pairs
      ! does not carry a kilobyte of empty buckets per thread
      ntab = hess_init_tab
      do while (ntab > 4 .and. cap*hess_load_den < (ntab/2)*hess_load_num)
         ntab = ntab/2
      end do
      allocate (self%htab(ntab), source=0)

   end subroutine hess_sparse_new

   !> Release one thread's accumulator
   !>
   !> @param[inout] self Accumulator of one thread
   subroutine hess_sparse_destroy(self)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self

      self%nent = 0
      self%maxent = 0
      if (allocated(self%pair_i)) deallocate (self%pair_i)
      if (allocated(self%pair_j)) deallocate (self%pair_j)
      if (allocated(self%blocks)) deallocate (self%blocks)
      if (allocated(self%htab)) deallocate (self%htab)

   end subroutine hess_sparse_destroy

   !> FNV-1a hash of an atom pair, masked to 32 bits
   !>
   !> Fed the two indices as whole words rather than byte by byte. Multiplying
   !> by an odd constant is a bijection modulo any power of two, so the low bits
   !> this hash is indexed by are a permutation of a contiguous run of atom ids
   !> -- which is exactly the run one grid point's neighbourhood hands it.
   !>
   !> @param[in] iatom Row atom
   !> @param[in] jatom Column atom
   !> @return          Hash value in `[0, 2^32)`
   pure function hess_pair_hash(iatom, jatom) result(h)
      !> Row and column atom
      integer, intent(in) :: iatom, jatom
      !> Hash value
      integer(int64) :: h

      h = ieor(fnv_offset, iand(int(iatom, int64), mask32))
      h = iand(h*fnv_prime, mask32)
      h = ieor(h, iand(int(jatom, int64), mask32))
      h = iand(h*fnv_prime, mask32)

   end function hess_pair_hash

   !> Index of the `(3, 3)` block of one atom pair, creating it if it is new
   !>
   !> A new entry is zeroed here and nowhere else, so an entry that is created
   !> and then never written contributes an exact zero, and one that is created
   !> twice is impossible. Both are what makes the merge a single add per
   !> element per thread; see [[drop_hess_sparse_type]].
   !>
   !> @param[inout] self  Accumulator of one thread
   !> @param[in]    iatom Row atom
   !> @param[in]    jatom Column atom
   !> @param[out]   ient  Entry index of the pair
   subroutine hess_sparse_entry(self, iatom, jatom, ient)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self
      !> Row and column atom
      integer, intent(in) :: iatom, jatom
      !> Entry index of the pair
      integer, intent(out) :: ient

      integer :: slot

      slot = int(iand(hess_pair_hash(iatom, jatom), &
                      int(size(self%htab) - 1, int64)), kind(slot)) + 1
      do
         ient = self%htab(slot)
         if (ient == 0) exit
         if (self%pair_i(ient) == iatom .and. self%pair_j(ient) == jatom) return
         slot = slot + 1
         if (slot > size(self%htab)) slot = 1
      end do

      call hess_grow_entries(self)
      self%nent = self%nent + 1
      ient = self%nent
      self%pair_i(ient) = iatom
      self%pair_j(ient) = jatom
      self%blocks(:, :, ient) = 0.0_wp

      ! The probe above found this slot on the *current* table, so it is only
      ! valid while that table stands; a growth reinserts every entry instead.
      if (self%nent*hess_load_den >= size(self%htab)*hess_load_num) then
         call hess_grow_table(self)
      else
         self%htab(slot) = ient
      end if

   end subroutine hess_sparse_entry

   !> Double the entry arrays if the next entry would not fit
   !>
   !> `move_alloc` on a copy rather than a reallocation in place: entry indices
   !> are handed out to the caller and must not move, and the accumulated blocks
   !> are copied verbatim, never re-derived.
   !>
   !> @param[inout] self Accumulator of one thread
   subroutine hess_grow_entries(self)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self

      integer, allocatable :: new_i(:), new_j(:)
      real(wp), allocatable :: new_blocks(:, :, :)
      integer :: cap, new_cap

      cap = size(self%pair_i)
      if (self%nent < cap) return

      new_cap = min(2*cap, self%maxent)
      if (new_cap <= cap) then
         ! Unreachable: `maxent` counts every pair that exists, and each is
         ! entered once, so a full container is never asked for another
         error stop "moist DROP Hessian: the sparse accumulator was asked for more "// &
            "pairs than the molecule has"
      end if

      allocate (new_i(new_cap))
      new_i(1:cap) = self%pair_i(1:cap)
      call move_alloc(new_i, self%pair_i)

      allocate (new_j(new_cap))
      new_j(1:cap) = self%pair_j(1:cap)
      call move_alloc(new_j, self%pair_j)

      allocate (new_blocks(ndim, ndim, new_cap))
      new_blocks(:, :, 1:cap) = self%blocks(:, :, 1:cap)
      call move_alloc(new_blocks, self%blocks)

   end subroutine hess_grow_entries

   !> Double the bucket table and reinsert every entry
   !>
   !> Touches no block and no entry index -- only the pair-to-entry lookup is
   !> rebuilt, so neither the accumulated values nor the order they were
   !> accumulated in can move.
   !>
   !> @param[inout] self Accumulator of one thread
   subroutine hess_grow_table(self)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self

      integer, allocatable :: new_tab(:)
      integer :: ient, slot

      allocate (new_tab(size(self%htab)*2), source=0)
      call move_alloc(new_tab, self%htab)

      do ient = 1, self%nent
         slot = int(iand(hess_pair_hash(self%pair_i(ient), self%pair_j(ient)), &
                         int(size(self%htab) - 1, int64)), kind(slot)) + 1
         do
            if (self%htab(slot) == 0) then
               self%htab(slot) = ient
               exit
            end if
            slot = slot + 1
            if (slot > size(self%htab)) slot = 1
         end do
      end do

   end subroutine hess_grow_table

   !> Add one thread's accumulator to the caller's Hessian
   !>
   !> One add per element, on the elements this thread reached. Called in the
   !> fixed `1 .. nthreads` order, which is the second half of the
   !> bit-reproducibility argument in [[drop_hess_sparse_type]].
   !>
   !> @param[in]    self    Accumulator of one thread
   !> @param[inout] hessian Nuclear-Hessian accumulator (3, nsph, 3, nsph)
   subroutine hess_sparse_reduce(self, hessian)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(in) :: self
      !> Nuclear-Hessian accumulator
      real(wp), intent(inout) :: hessian(:, :, :, :)

      integer :: ient, iatom, jatom

      do ient = 1, self%nent
         iatom = self%pair_i(ient)
         jatom = self%pair_j(ient)
         hessian(:, iatom, :, jatom) = hessian(:, iatom, :, jatom) &
                                       + self%blocks(:, :, ient)
      end do

   end subroutine hess_sparse_reduce

end submodule moist_cavity_drop_derivatives_hessian_traverse
