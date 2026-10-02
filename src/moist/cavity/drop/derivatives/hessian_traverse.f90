!> Single grid traversal of the DROP surface Hessian
!>
!> Directional derivative of the nuclear gradient `J^T omega`
!>
!>     d/dv [ J^T omega ]  =  (dJ^T/dv) omega  +  J^T (d omega/dv)
!>                            ^ fixed channel     ^ response channel
!>
!> ## Shared per point
!> - One [[drop_point_prologue]], one application of the 16 basis seeds
!>   ([[point_seed_basis]]), one jet-tensor fill, one iSwiG collection
!> - Weight sets and directions enter only at the contractions
!> - Entry points in `hessian.f90`; half-accessors here enable one channel alone
!>
!> ## Fixed channel, adjoints held fixed
!> - One [[fixed_direction_chain]] per direction, linear in the direction
!> - Column: field row over the active atoms, three owner anchor entries,
!>   switching block added outside the chain
!> - `drop_fixed_rank4`: direction-free `(3, nsph, 3, nsph)` block; chain on
!>   the 23 elements of [[chain_basis_element]], columns of a point are
!>   `L C` plus the level set's `vjp_f2_rArB` block
!> - `drop_fixed_per_dir`: one chain per supplied direction into the
!>   `(3, nsph, ndir)` columns; explicit nuclear motion from `hvp_jet_rA`,
!>   or `vjp_f2_rArB_apply` from `drop_hvp_apply_min` directions on
!> - Few directions: per-direction mode; a basis: rank-4; chosen in `hessian.f90`
!> - Never read a single entry of `dw_lsf2`: Hessian jet seeds are
!>   asymmetric, exact only contracted against the symmetric `lsf3_rr_rA`
!>
!> ## Branch weight in the fixed channel
!> - `d(wbranch)` couples all points of a multi-branch anchor group, the one
!>   chain input that is not point local
!> - Per direction: from the forward tangent, [[weight_tangents]] or
!>   [[branch_weight_tangents]]
!> - Rank-4: grid loop keeps [[drop_branch_point_type]] elements,
!>   [[branch_weight_block]] closes the term after the thread reduction
!> - Unbranched grid: one scan of `branch_count`, columns unchanged to the bit
!>
!> ## Response channel, primal map held fixed
!> - Gradient contraction with `d(eff)/dv` in place of the fold `eff`
!> - Pass 1: forward tangent, [[get_surface_tangent_drop]]
!> - Pass 2: [[prepare_surface_weights_tangent]], serial, one call per
!>   direction; [[branch_phi_adj_tangent]] must see each anchor group whole
!> - Without a model response `d(eff)` carries `w_xi`, `w_f` and
!>   `branch_phi_adj` only: normal channel sees a hard zero, curvature off
!> - Pass 2b, with `omega_v`: raw adjoint response of the model, folded at the
!>   base geometry and added to `d(eff)`; serial, outside the parallel region
!> - Curvature weights the primal `eff` lacks are refused
!> - Pass 3: contraction in the grid loop, no LSF call per direction
!>
!> ## Direction blocks
!> - Blocked over directions, `drop_hvp_chunk_dirs` and `block_bytes`; the
!>   whole grid is re-traversed per block
!> - One pass-1 walk and one contraction walk per block, not fusable
!> - Rank-4 channel accumulates in the first block only, `ilo == 1`
!> - Per-direction channel runs in every block
!> - Blocked equals unblocked to the bit: same additions in the same order
!>   per column, held by `hvp_direction_chunking`
!> - Prologue configuration identical in every block; an order-3 later
!>   block shifts columns by one ulp
!>
!> ## Failure contract
!> - Accumulators are added to and left untouched on failure
!> - Multi-block run: `hvp` restored from a copy taken up front
!> - Rank-4 block staged by `hessian.f90`, not here
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_hessian_traverse
!$ use omp_lib, only: omp_get_thread_num
   use, intrinsic :: iso_fortran_env, only: int64
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type, &
      & drop_point_scratch_type
   use moist_cavity_drop_derivatives_kernel, only: &
      & drop_seed_state_type, drop_seed_result_type, drop_seed_state_tangent_type, &
      & drop_seed_input_tangent_type, drop_n_host_jet, add_host_jet, &
      & drop_surface_weights_type, apply_seed, apply_seed_tangent, &
      & seed_weight_tol, next_branch_group, max_branch_group_size
   use moist_cavity_drop_derivatives_seeds, only: drop_n_jet_seeds, drop_n_anchor_seeds, &
      & drop_n_point_seeds, seed_normal_channel, fill_seed_basis, scatter_jet_weight, &
      & seed_contribution_tangent, seed_dzero1, seed_dzero2, &
      & seed_jet_basis_apply, seed_jet_basis_contract, &
      & seed_anchor_apply, seed_anchor_contract
   use moist_cavity_drop_derivatives_field_tangent, only: drop_field_tangent_point, &
      & drop_field_jet_point, drop_field_jet_contract, drop_field_jet_tangent, &
      & drop_field_jet_column_packed, drop_field_f4_fold, drop_field_tangent_dir, &
      & drop_field_tangent_work_type, drop_n_sym2, drop_n_jet_coef, &
      & drop_sym2_idx, drop_sym3_idx
   use moist_math_blas, only: gemm
   use moist_cavity_drop_derivatives_weights_tangent, only: prepare_surface_weights_tangent
   use moist_cavity_drop_gaussian_scatter, only: scatter_iswig_block_indexed, &
      & contract_iswig_block
   implicit none(type, external)

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Memory cap of one direction block of the response channel, bytes
   !>
   !> - Only shortens a block on a large grid; `drop_hvp_chunk_dirs` stays the
   !>   upper bound
   !> - Blocked and unblocked runs agree to the bit
   integer(int64), parameter :: block_bytes = 512_int64*1024_int64*1024_int64
   !> Doubles per grid point and direction of a block
   !>
   !> - Surface tangent (12), branch-weight tangent (1), moving fold channels (3)
   integer(int64), parameter :: block_doubles_plain = 16_int64
   !> Same with a model adjoint response
   !>
   !> - Adds the raw adjoints (12) and the copied fold channels (8)
   integer(int64), parameter :: block_doubles_omega = 36_int64
   !> Doubles added by the host's jet tangents
   integer(int64), parameter :: block_doubles_jets = 40_int64
   !> Doubles added by the response tangent
   !>
   !> - Position, surface charge and level-set weight tangents
   integer(int64), parameter :: block_doubles_rt = 16_int64

   !> Dimension of the chain's linear input space
   !>
   !> - Packed jet tangent of a direction plus the owner displacement, see
   !>   [[chain_basis_element]]
   integer, parameter :: drop_n_chain_basis = drop_n_jet_coef + ndim

   !> Initial entry capacity of a per-thread accumulator
   !>
   !> - Not a bound: grows on demand, trimmed at birth for a small molecule
   integer, parameter :: hess_init_ent = 256
   !> Initial bucket count of a per-thread accumulator
   !>
   !> - Not a bound: grows on demand, trimmed at birth for a small molecule
   !> - Must stay a power of two, masked as an index in [[hess_probe]]
   integer, parameter :: hess_init_tab = 1024

   !> Load factor above which the bucket table doubles, as `num/den`
   integer, parameter :: hess_load_num = 7, hess_load_den = 10

   !> FNV-1a 32-bit offset basis of [[hess_pair_hash]]
   integer(int64), parameter :: fnv_offset = 2166136261_int64
   !> FNV-1a 32-bit prime
   integer(int64), parameter :: fnv_prime = 16777619_int64
   !> 32-bit mask
   !>
   !> - Applied after every product, intermediates stay in the signed 64-bit range
   integer(int64), parameter :: mask32 = int(z'FFFFFFFF', int64)

   !> Sparse per-thread accumulator of the rank-4 fixed-adjoint Hessian
   !>
   !> - One thread's share of `(3, nsph, 3, nsph)` as the `(3, 3)` blocks of
   !>   the atom pairs it reaches; pair count grows like `nsph * k`
   !> - Reproduces a dense slab to the bit; accumulation and merge order are
   !>   load bearing
   !> - Accumulation: added in place in grid order, one entry per pair; entry
   !>   indices never move on growth, bucket table rebuilt from the entries
   !> - Merge: [[hess_sparse_reduce]] in fixed `1 .. nthreads` order, one add
   !>   per element, untouched pairs skipped
   !> - Open addressing with linear probing, entry indices as values
   type :: drop_hess_sparse_type
      !> Entries in use
      integer :: nent = 0
      !> Pairs that exist at all, `nsph^2` clamped to the integer range
      !>
      !> - Caps the doubling entry arrays at the dense size
      !> - Exact: no distinct pair beyond it can be requested
      integer :: maxent = 0
      !> Row and column atom of each entry
      integer, allocatable :: pair_i(:), pair_j(:)
      !> Accumulated `(3, 3)` block of each entry, `(row axis, column axis, entry)`
      real(wp), allocatable :: blocks(:, :, :)
      !> Open-addressing bucket table; values are entry indices, `0` is empty
      integer, allocatable :: htab(:)
   end type drop_hess_sparse_type

   !> Basis seeds of one grid point and their responses, 16 per point
   !>
   !> - Filled once per point by [[fill_seed_basis]] and [[point_seed_basis]]
   !> - Read by both channels
   !> - One object per thread
   type :: drop_point_seeds_type
      !> Seed perturbations of the level-set gradient
      real(wp) :: dlsf1(3, drop_n_point_seeds) = 0.0_wp
      !> Seed perturbations of the level-set Hessian
      real(wp) :: dlsf2(3, 3, drop_n_point_seeds) = 0.0_wp
      !> Induced point motion and multiplier change of each seed
      real(wp) :: x(4, drop_n_point_seeds) = 0.0_wp
      !> Linear response of each seed
      type(drop_seed_result_type) :: res(drop_n_point_seeds)
      !> Derived-state tangent along each seed, fixed channel only
      type(drop_seed_state_tangent_type) :: dstate(drop_n_point_seeds)
   end type drop_point_seeds_type

   !> One weight set contracted at a grid point
   !>
   !> - Fixed channel: everything, from `eff`
   !> - Response channel: level-set adjoints from `deff(jdir)`, `w_xyz` only
   !>   with a model response; the rest stays zero, exact for that channel
   type :: drop_point_weights_type
      !> Point-local level-set adjoint weights
      real(wp) :: w_lsf0 = 0.0_wp, w_lsf1(3) = 0.0_wp, w_lsf2(3, 3) = 0.0_wp
      !> Effective position adjoint and the normal fold it is built from
      real(wp) :: w_xyz(3) = 0.0_wp, normal_grad(3) = 0.0_wp, nwn = 0.0_wp
   end type drop_point_weights_type

   !> Per-direction temporaries of [[fixed_direction_chain]]
   !>
   !> - One object per thread, not default-initialised on every call
   type :: drop_fixed_dir_scratch_type
      !> Right-hand side and solution of the directional projection response
      real(wp) :: rhs_v(4, 1) = 0.0_wp, dr_v(3) = 0.0_wp, dl_v = 0.0_wp
      !> Response of the direction
      type(drop_seed_result_type) :: res_v
      !> Derived-state tangent of the direction
      type(drop_seed_state_tangent_type) :: dstate_v
      !> Input tangent of the direction
      type(drop_seed_input_tangent_type) :: dinp_v
      !> Second-order response of one seed along the direction
      type(drop_seed_result_type) :: dres
      !> Tangent of the bordered KKT system
      real(wp) :: dH_lag(3, 3, 1) = 0.0_wp, dg_tot(3, 1) = 0.0_wp
      !> Tangent of the 16 seed right-hand sides
      real(wp) :: dseed_x(4, drop_n_point_seeds) = 0.0_wp
      !> Tangents of the normal fold and of the effective position adjoint
      real(wp) :: dnormal_grad(3) = 0.0_wp, dw_xyz(3) = 0.0_wp
      !> Tangents of the level-set adjoint weights
      real(wp) :: dw_lsf0 = 0.0_wp, dw_lsf1(3) = 0.0_wp, dw_lsf2(3, 3) = 0.0_wp
      !> Tangent of the objective gradient, `alpha (dr_v - v_owner)`
      real(wp) :: dphi1_v(3) = 0.0_wp
   end type drop_fixed_dir_scratch_type

   !> Local direction set and pair table of one point, rank-4 mode
   !>
   !> - `pair_ent` grown from what a point needs, never sized to `nsph^2`
   type :: drop_rank4_scratch_type
      !> Atom count of the local direction set, and the owner's position in it
      integer :: ndir_atom = 0, owner_loc = 0
      !> Atoms of the local direction set, active atoms first in slot order
      integer, allocatable :: dir_atoms(:)
      !> Accumulator entry of every `(row, column)` pair the point reaches
      integer, allocatable :: pair_ent(:, :)
      !> Current Cartesian unit direction, `(3, nsph)`, otherwise zero
      real(wp), allocatable :: vdir(:, :)
      !> Weighted mixed nuclear Hessian block of the point over its active slots
      !>
      !> - `(3 n_active, 3 n_active)`, Cartesian component fastest
      !> - Sized exactly to the point, written in place by the level set
      real(wp), allocatable :: mblk(:, :)
      !> Chain on the symmetrised basis, `(3 n_active + 3, 23)`
      !>
      !> - Field rows of the active slots, Cartesian component fastest, then the
      !>   three anchor entries of the owner
      real(wp), allocatable :: lmat(:, :)
      !> Basis coordinates of every unit direction of the local set, `(23, 3 n_local)`
      real(wp), allocatable :: cmat(:, :)
      !> All columns of the point, `lmat * cmat` plus the explicit block
      !>
      !> - Sized exactly, product lands in place
      real(wp), allocatable :: oblk(:, :)
   end type drop_rank4_scratch_type

   !> Anchor's iSwiG neighbourhood as one point sees it
   type :: drop_swi_scratch_type
      !> Whether the fixed and the response channel want the neighbourhood
      logical :: fixed = .false., resp = .false.
      !> Collected switching factor, and the response channel's owner row
      real(wp) :: f0 = 0.0_wp, owner_row(3) = 0.0_wp, dxi = 0.0_wp
      !> Influence-set size of the second-derivative block
      integer :: n = 0
      !> Influence set and its pair entries
      !>
      !> - Grown on demand from `n_nb + 1`, never molecule sized
      integer, allocatable :: idx(:), ent(:, :)
      !> Second-derivative block of the influence set, grown with `idx`
      real(wp), allocatable :: blk(:, :, :, :)
   end type drop_swi_scratch_type

   !> Branch element of one grid point of a multi-branch anchor group
   !>
   !> - Point-local factors of the branch-weight motion, rank-4 mode
   !> - Combined per group by [[branch_weight_block]] after the grid loop
   !> - One element per branched point, written by the one thread walking it
   type :: drop_branch_point_type
      !> Atoms the rows refer to
      !>
      !> - Active atoms in slot order, then the owner
      !> - Owner listed again when it is active as well
      integer, allocatable :: atoms(:)
      !> Chain on a unit branch-weight tangent, `(3, size(atoms))`
      !>
      !> - Jet and owner at rest
      real(wp), allocatable :: lrow(:, :)
      !> Nuclear gradient of the branch objective `Phi`, `(3, size(atoms))`
      real(wp), allocatable :: prow(:, :)
   end type drop_branch_point_type

contains

   !* ================================================================================= *!
   !*                            The half-accessor wrappers                             *!
   !* ================================================================================= *!

   !> Accumulate the rank-4 fixed-adjoint half of the surface Hessian
   !>
   !> - `(dJ^T/dv) omega` as the direction-free block, added to `hessian`
   !> - `hessian` untouched on failure
   !> - No weight guard: `eff` arrives folded, a live `w_a` or `w_w` belongs to
   !>   the response channel
   !> - Entry point of `test_cavity_drop_hessian_fixed`
   !>
   !> @param[in]     self     DROP cavity instance, must hold a projected grid
   !> @param[in]     eff      folded surface adjoints, held fixed
   !> @param[in,out] hessian  nuclear-Hessian accumulator `(3, nsph, 3, nsph)`
   !> @param[out]    error    allocated on failure
   module subroutine get_surface_hessian_fixed_drop(self, eff, hessian, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Folded surface adjoints from [[prepare_surface_weights]]
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

   !> Accumulate the per-direction fixed-adjoint half of the surface Hessian
   !>
   !> - `(dJ^T/dv) omega`, one gradient column per direction, added to `hvp`
   !> - `hvp` untouched on failure
   !> - Equals the rank-4 half contracted with `dirs`, cross-checked by the
   !>   fixed suite
   !> - Multi-branch grid: every block also runs the forward tangent for
   !>   `d(wbranch)`
   !>
   !> @param[in]     self   DROP cavity instance, must hold a projected grid
   !> @param[in]     eff    folded surface adjoints, held fixed
   !> @param[in]     dirs   nuclear directions `(3, nsph, ndir)`
   !> @param[in,out] hvp    Hessian-vector accumulator `(3, nsph, ndir)`
   !> @param[out]    error  allocated on failure
   module subroutine get_surface_hessian_fixed_dirs_drop(self, eff, dirs, hvp, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Folded surface adjoints from [[prepare_surface_weights]]
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

   !> Accumulate the adjoint-response half of the surface Hessian
   !>
   !> - `J^T (d omega/dv)`, one gradient column per direction, added to `hvp`
   !> - `hvp` untouched on failure
   !> - Needs both forms of the adjoints: raw `acc` and its primal fold `eff`
   !>
   !> @param[in]     self   DROP cavity instance, must hold a projected grid
   !> @param[in]     acc    raw surface-observable adjoints, held fixed
   !> @param[in]     eff    folded surface adjoints of the base geometry
   !> @param[in]     dirs   nuclear directions `(3, nsph, ndir)`
   !> @param[in,out] hvp    Hessian-vector accumulator `(3, nsph, ndir)`
   !> @param[out]    error  allocated on failure
   module subroutine get_surface_hessian_response_drop(self, acc, eff, dirs, hvp, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Raw surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Folded surface adjoints from [[prepare_surface_weights]]
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
   !> - Messages prefixed with `context`, as in [[check_surface_adjoint]]
   !>
   !> @param[in]  self     DROP cavity instance
   !> @param[in]  dirs     nuclear directions, expected `(3, nsph, ndir)`, `ndir > 0`
   !> @param[in]  hvp      column accumulator, expected `(3, nsph, ndir)`
   !> @param[in]  context  calling routine, prefixes the diagnostics
   !> @param[out] error    allocated on a shape mismatch
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

   !> Walk the projected grid once and accumulate into the enabled channels
   !>
   !> - Thread setup, grid loop, error latching and deterministic reduction as
   !>   in [[get_surface_gradient_drop]]
   !> - `host`, and `rt` with level-set weights, need the per-direction fixed mode
   !> - Host-defined level set without `host` is refused
   !> - `hvp` untouched on failure; `hess_fixed` is staged by the caller
   !>
   !> @param[in]     self           DROP cavity instance, must hold a projected grid
   !> @param[in]     eff            folded surface adjoints of the base geometry
   !> @param[in]     fixed_mode     fixed-channel mode, one of the `drop_fixed_*` constants
   !> @param[in]     want_response  accumulate `J^T (d omega/dv)`
   !> @param[in]     context        calling routine, prefixes the diagnostics
   !> @param[in]     acc            raw surface adjoints; required by the response channel
   !> @param[in]     dirs           `(3, nsph, ndir)`; required by every per-direction channel
   !> @param[in,out] hess_fixed     `(3, nsph, 3, nsph)`; required by the rank-4 mode
   !> @param[in,out] hvp            `(3, nsph, ndir)`; required with `dirs`
   !> @param[out]    error          allocated on failure
   !> @param[in,out] omega_v        model adjoint response; needs the response channel
   !> @param[in,out] host           host tangents along `dirs`; needs the response channel
   !> @param[in,out] rt             initialised for `(ngrid, ndir)`; needs the response channel
   module subroutine drop_hessian_traverse(self, eff, fixed_mode, want_response, context, &
                                           acc, dirs, hess_fixed, hvp, error, omega_v, host, rt)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Folded surface adjoints from [[prepare_surface_weights]]
      type(drop_surface_weights_type), intent(in) :: eff
      !> Mode of the fixed channel
      integer, intent(in) :: fixed_mode
      !> Whether the response channel is accumulated
      logical, intent(in) :: want_response
      !> Calling routine, named in the diagnostics
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
      !> Surface-adjoint response of the model
      class(surface_adjoint_response_type), intent(inout), optional :: omega_v
      !> Second-order host exchange, host tangents along the directions
      class(coupling_tangent_type), intent(inout), optional :: host
      !> Response tangent of the whole direction set, filled per block
      type(response_tangent_type), intent(inout), optional :: rt

      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type) :: slots
      !> Whether the raw adjoints move through a model response
      logical :: have_omega
      !> Host tangents supplied, host-defined level set, `rt` filled, with level-set weights
      logical :: have_host, have_jets, want_rt, want_rt_lsf
      !> Host's jet tangents of one block, `(40, ngrid, nblk)`
      real(wp), allocatable :: host_jets(:, :, :)
      !> Non-surface columns of the model response of one block, staged outside `hvp`
      real(wp), allocatable :: hvp_direct(:, :, :)
      !> Per-thread rank-4 buffers, reduced in fixed thread order
      type(drop_hess_sparse_type), allocatable :: hess_threads(:)
      !> Per-thread column buffers, reduced in fixed thread order
      real(wp), allocatable :: hvp_threads(:, :, :, :)
      !> Entry copy of `hvp`, restores it when a later block fails
      real(wp), allocatable :: hvp_entry(:, :, :)
      !> Tangent of the folded surface adjoints per direction of the block
      type(drop_surface_weights_type), allocatable :: deff(:)
      !> Thread bookkeeping
      integer :: thread_slot, ithread
      !> First failure seen in the parallel region
      type(drop_abort_latch_type) :: abort
      !> Per-thread failure on its way to the latch
      type(error_type), allocatable :: worker_error

      !> Whether any channel writes gradient columns
      logical :: want_cols
      !> Prologue configuration: jet order and curvature request
      integer :: max_deriv
      logical :: want_curvature
      !> Direction count, block size, block bounds in `dirs`, block extent
      integer :: ndir, nchunk, ilo, ihi, nblk
      !> Per-direction working set of a block, doubles per grid point
      integer(int64) :: block_doubles
      !> Fixed channel this block runs, if any
      logical :: fixed_rank4_here, fixed_per_dir_here, fixed_here

      !> Per-thread state of the grid point
      type(drop_point_scratch_type) :: pt
      !> Whether the shared prologue cleared the point
      logical :: point_ok
      !> Basis seeds of the point and their responses
      type(drop_point_seeds_type) :: seeds
      !> Fixed adjoints contracted at the point, and one response weight set
      type(drop_point_weights_type) :: pw, rw
      !> Per-direction temporaries of the second-order chain
      type(drop_fixed_dir_scratch_type) :: sc
      !> Jet tensors and scratch of the field contraction
      type(drop_field_tangent_work_type) :: ft_work
      !> Rank-4 direction set and pair table
      type(drop_rank4_scratch_type) :: r4
      !> Anchor's iSwiG neighbourhood
      type(drop_swi_scratch_type) :: swi

      !> Whether the grid carries a multi-branch anchor group
      logical :: branched
      !> Per-direction mode: branch-weight tangent of the block, `(ngrid, nblk)`
      real(wp), allocatable :: dwb_blk(:, :)
      !> Its value at the point being walked
      real(wp) :: dwb_dir
      !> Rank-4 mode: element of each multi-branch point in `br_pts`, zero elsewhere
      integer, allocatable :: br_slot(:)
      !> Branch elements left for [[branch_weight_block]]
      type(drop_branch_point_type), allocatable :: br_pts(:)
      !> Element of the point being walked
      integer :: islot

      !> Jet tangent along the current direction at the frozen point
      real(wp) :: dv0, dv1(3), dv2(3, 3), dv3(3, 3, 3)
      !> Owner displacement of a basis element
      real(wp) :: vown(3)
      !> Row count of the point's column block, basis element
      integer :: nrow, ibasis
      !> Gradient column of the current direction: field row and owner entries
      real(wp), allocatable :: field_row(:, :)
      real(wp) :: anchor_row(3)
      !> Direction restricted to the point's active slots
      real(wp), allocatable :: v_act(:, :)
      !> Explicit nuclear motion of the field row per block direction, applied block
      real(wp), allocatable :: xrow(:, :, :)
      !> Whether the explicit nuclear motion comes from the applied block
      logical :: use_apply
      !> Effective position adjoint of the response channel, identically zero
      real(wp) :: w_xyz_zero(3)

      !> Grid, atom, axis, seed, slot and direction indices
      integer :: igrid, ient, i, iaxis
      integer :: idir, dir_axis, dir_loc, ia, ib
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
      have_omega = present(omega_v)
      if (have_omega .and. .not. want_response) then
         call fatal_error(error, context//": a model adjoint response needs the response"// &
                          " Hessian channel")
         return
      end if
      ! Host jet tangents enter both channels
      ! Level-set weights of `rt` need the per-direction chain's `dw_lsf`
      have_host = present(host)
      have_jets = have_host .and. self%lsf_model%host_defined
      ! Host-defined level set without host exchange: surface-motion term silently zero
      if (self%lsf_model%host_defined .and. .not. have_host) then
         call fatal_error(error, context//": this cavity's level set is the host's, so its"// &
                          " tangents along the directions are host data; the nuclear Hessian"// &
                          " needs the second-order host exchange to supply them")
         return
      end if
      want_rt = present(rt)
      want_rt_lsf = .false.
      if (want_rt) want_rt_lsf = rt%want_lsf()
      if ((have_host .or. want_rt) .and. .not. (want_response .and. present(dirs))) then
         call fatal_error(error, context//": the second-order host exchange needs the"// &
                          " response Hessian channel and a direction set")
         return
      end if
      if ((have_host .or. want_rt_lsf) .and. fixed_mode /= drop_fixed_per_dir) then
         call fatal_error(error, context//": the second-order host exchange needs the"// &
                          " per-direction fixed Hessian channel")
         return
      end if
      if (want_rt) then
         if (.not. rt%is_initialized(self%ngrid, size(dirs, 3))) then
            call fatal_error(error, context//": the response tangent is not initialised"// &
                             " for this grid and direction set")
            return
         end if
      end if
      if (self%ngrid <= 0) return

      want_cols = want_response .or. fixed_mode == drop_fixed_per_dir

      ! Multi-branch anchor group: fixed channel carries a branch-weight motion
      branched = .false.
      if (fixed_mode /= drop_fixed_none .and. allocated(self%branch_count) &
          .and. allocated(self%anchor_id)) then
         branched = any(self%branch_count(1:self%ngrid) > 1)
      end if

      if (fixed_mode == drop_fixed_per_dir) then
         if (want_response) then
            h_hess = self%ctx%timer%resolve("Surface Hessian (per direction)", &
                                            self%ctx%timer%current(), cat_hessian)
         else
            h_hess = self%ctx%timer%resolve("Surface Hessian (fixed adjoint, per direction)", &
                                            self%ctx%timer%current(), cat_hessian)
         end if
      else if (fixed_mode == drop_fixed_rank4) then
         if (want_response) then
            h_hess = self%ctx%timer%resolve("Surface Hessian (both halves)", &
                                            self%ctx%timer%current(), cat_hessian)
         else
            h_hess = self%ctx%timer%resolve("Surface Hessian (fixed adjoint)", &
                                            self%ctx%timer%current(), cat_hessian)
         end if
      else
         h_hess = self%ctx%timer%resolve("Surface Hessian (adjoint response)", &
                                         self%ctx%timer%current(), cat_hessian)
      end if
      call self%ctx%timer%start(h_hess)

      !* -------------------------------- Thread setup -------------------------------- *!
      ! Fixed channel: order 4 for `f4_rrrr` and `f4_rrr_rA`, CFC takes the
      ! highest total order
      ! Response only: order 3, no curvature
      ! Same configuration in every block, see the module header
      if (fixed_mode /= drop_fixed_none) then
         max_deriv = 4
         want_curvature = eff%have_wk
      else
         max_deriv = 3
         want_curvature = .false.
      end if
      call slots%init(self%ctx, self%lsf_model, max_deriv, self%param, self%mol, self%radii)

      if (fixed_mode == drop_fixed_rank4) allocate (hess_threads(slots%nthreads))
      if (fixed_mode == drop_fixed_rank4 .and. branched) then
         call branch_point_slots(self, br_slot, br_pts)
      end if

      !* ------------------------------ Direction blocks ------------------------------ *!
      ! Single block when the direction set fits; rank-4 only: one block, no
      ! directions
      ! Entry copy of `hvp` restores it on failure, taken only for a multi-block run
      if (want_cols) then
         ndir = size(dirs, 3)
         nchunk = min(ndir, drop_hvp_chunk_dirs)
         ! Working set kept under `block_bytes`; memory only, results unchanged
         ! to the bit
         if (want_response) then
            block_doubles = merge(block_doubles_omega, block_doubles_plain, have_omega)
            if (have_jets) block_doubles = block_doubles + block_doubles_jets
            if (want_rt) block_doubles = block_doubles + block_doubles_rt
            nchunk = min(nchunk, max(1, int(block_bytes/ &
                                            (8_int64*block_doubles*int(self%ngrid, int64)))))
         else if (branched) then
            ! Fixed-only branched traversal holds the forward tangent's plain working set
            nchunk = min(nchunk, max(1, int(block_bytes/ &
                                            (8_int64*block_doubles_plain*int(self%ngrid, int64)))))
         end if
         allocate (hvp_threads(ndim, self%nsph, nchunk, slots%nthreads))
         if (nchunk < ndir) allocate (hvp_entry, source=hvp)
         if (have_omega) allocate (hvp_direct(ndim, self%nsph, nchunk))
      else
         ndir = 1
         nchunk = 1
      end if
      ! Applied block from `drop_hvp_apply_min` directions on, else `hvp_jet_rA`
      ! per direction
      ! Keyed on the whole set, never on a block: columns independent of blocking
      use_apply = fixed_mode == drop_fixed_per_dir .and. ndir >= drop_hvp_apply_min
      if (have_jets) allocate (host_jets(drop_n_host_jet, self%ngrid, nchunk))

      do ilo = 1, ndir, nchunk
         ihi = min(ilo + nchunk - 1, ndir)
         nblk = ihi - ilo + 1

         !* -------------------- The host's jet tangents of the block -------------------- *!
         ! Once per block, read by the tangent pass and the contraction walk
         if (have_jets) then
            call host%level_set_tangent(ilo, dirs(:, :, ilo:ihi), self%xyz(:, 1:self%ngrid), &
                                        host_jets(:, :, 1:nblk), error)
            if (allocated(error)) then
               if (allocated(hvp_entry)) hvp = hvp_entry
               call self%ctx%timer%stop(h_hess)
               return
            end if
         end if

         ! Grid loop reads these flags, never the mode
         fixed_rank4_here = fixed_mode == drop_fixed_rank4 .and. ilo == 1
         fixed_per_dir_here = fixed_mode == drop_fixed_per_dir
         fixed_here = fixed_rank4_here .or. fixed_per_dir_here

         if (want_response) then
            !* ------------------ Passes 1 and 2: the moving weights ------------------ *!
            ! Serial over the block's directions; `deff` indexed `1 .. nblk`
            if (have_omega) hvp_direct = 0.0_wp
            if (have_jets) then
               call weight_tangents(self, acc, eff, dirs(:, :, ilo:ihi), &
                                    fixed_mode == drop_fixed_per_dir, want_curvature, context, &
                                    deff, dwb_blk, error, omega_v=omega_v, &
                                    hvp_direct=hvp_direct, rt=rt, first=ilo, &
                                    host_jets=host_jets(:, :, 1:nblk))
            else
               call weight_tangents(self, acc, eff, dirs(:, :, ilo:ihi), &
                                    fixed_mode == drop_fixed_per_dir, want_curvature, context, &
                                    deff, dwb_blk, error, omega_v=omega_v, &
                                    hvp_direct=hvp_direct, rt=rt, first=ilo)
            end if
            if (allocated(error)) then
               if (allocated(hvp_entry)) hvp = hvp_entry
               call self%ctx%timer%stop(h_hess)
               return
            end if
         else if (fixed_per_dir_here .and. branched) then
            !* -------------- Branch-weight tangents of a fixed-only block -------------- *!
            ! No response channel to supply `d(wbranch)`, forward tangent run for it
            ! Contracted, so the fixed half is the same column with or without the
            ! response half
            call branch_weight_tangents(self, dirs(:, :, ilo:ihi), dwb_blk, error)
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
         !$omp& dv0, dv1, dv2, dv3, vown, nrow, ibasis, field_row, anchor_row, v_act, xrow, &
         !$omp& w_xyz_zero, ient, i, iaxis, idir, dir_axis, dir_loc, ia, ib, &
         !$omp& jdir, idir_g, iatom, jj, knb, dwb_dir, islot, worker_error)
         thread_slot = 1
!$       thread_slot = omp_get_thread_num() + 1

         ! `want_lsf4`: fourth-derivative buffer, fixed channel only, pure output
         ! `want_vjp`: row buffer of the response channel
         call pt%init(self%nsph, self%iswig, want_seed_batch=.true., &
                      want_lsf4=fixed_here, want_vjp=want_response)
         ! Hard zero for the response channel's contractions
         w_xyz_zero = 0.0_wp

         if (fixed_here) then
            allocate (field_row(3, self%nsph), source=0.0_wp)
            ! Grown on demand from `n_nb + 1`, quadratic in the influence set
            allocate (swi%blk(3, 1, 3, 1), swi%idx(1))
         end if
         if (fixed_rank4_here) then
            call hess_sparse_new(hess_threads(thread_slot), self%nsph)
            allocate (r4%dir_atoms(self%nsph))
            ! One slot each, grown by the first point
            allocate (r4%pair_ent(1, 1), swi%ent(1, 1))
            allocate (r4%vdir(3, self%nsph), source=0.0_wp)
         end if
         if (fixed_per_dir_here) then
            allocate (v_act(3, self%nsph), source=0.0_wp)
            if (use_apply) allocate (xrow(3, self%nsph, nchunk))
         end if

         !$omp do schedule(static, 8)
         do igrid = 1, self%ngrid
            !* ------------------------ Shared point prologue ------------------------- *!
            ! Every failure path has latched itself
            ! `pt%lsf4_rrrr` filled when a fixed channel asked for it
            call drop_point_prologue(self, slots, thread_slot, igrid, want_curvature, &
                                     context, abort, pt, point_ok)
            if (.not. point_ok) cycle

            !* ------------------- Basis seeds and their responses -------------------- *!
            ! Applied once for both channels
            ! `fill_seed_basis` adds the seed directions of the second-order chain;
            ! `seeds%x` is rewritten below with the same values
            if (fixed_here) then
               call fill_seed_basis(pt%kkt_rhs, seeds%dlsf1, seeds%dlsf2, seeds%x)
            end if
            call point_seed_basis(pt%state, pt%kkt_rhs, fixed_here, seeds%res, seeds%x, &
                                  seeds%dstate)

            !* ----------------------- Jet tensors of the point ----------------------- *!
            ! Direction free, formed once per point
            ! Fills are unconditional: readers abort on a stale slot marker
            if (fixed_here) then
               call drop_field_tangent_point(slots%lsf(thread_slot)%lsf, ft_work)
            else
               call drop_field_jet_point(slots%lsf(thread_slot)%lsf, ft_work)
            end if

            !* --------------- Shared iSwiG neighbourhood of the anchor --------------- *!
            ! One `swi_collect` for both channels
            ! Fixed channel: second-derivative block, only with a live weight
            ! Response channel: first-derivative rows, only with a live tangent
            ! Width outputs unused: switching factor taken at the anchor width
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
               call self%iswig%swi_collect(pt%state%anchor, pt%owner_idx, self%anchor_xi0(igrid), &
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
               ! Hessian weights folded into the mixed fourth derivative once per point
               call drop_field_f4_fold(ft_work, pt%n_active, pw%w_lsf2)
               ! Base level-set weights of `rt`, direction free, first block only
               ! Hessian weight symmetrised: only the transpose-pair sum of the
               ! single-entry seeds is meaningful
               if (want_rt_lsf .and. ilo == 1) then
                  rt%w_value(igrid) = pw%w_lsf0
                  rt%w_gradient(:, igrid) = pw%w_lsf1
                  rt%w_hessian(:, :, igrid) = 0.5_wp*(pw%w_lsf2 + transpose(pw%w_lsf2))
               end if
            end if

            if (fixed_rank4_here) then
               !* ------------------------ Local direction set ------------------------ *!
               ! Active atoms plus owner, every other column of the point is zero
               ! Active atoms first: a local index up to `n_active` is the slot
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

               ! Pair square resolved once per point; an untouched entry holds an exact zero
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
               ! Explicit half of every column in one accessor call
               ! Sized exactly to the point, written in place
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

               !* ------------------ The chain on the symmetrised basis ------------------ *!
               ! Elements `1 .. 20`: packed jet classes, owner at rest
               ! Elements `21 .. 23`: owner alone
               ! Field row and anchor entries of an element form one column of `L`
               nrow = 3*pt%n_active + 3
               if (allocated(r4%lmat)) then
                  if (size(r4%lmat, 1) /= nrow .or. size(r4%oblk, 2) /= 3*r4%ndir_atom) then
                     deallocate (r4%lmat, r4%cmat, r4%oblk)
                  end if
               end if
               if (.not. allocated(r4%lmat)) then
                  allocate (r4%lmat(nrow, drop_n_chain_basis))
                  allocate (r4%cmat(drop_n_chain_basis, 3*r4%ndir_atom))
                  allocate (r4%oblk(nrow, 3*r4%ndir_atom))
               end if

               do ibasis = 1, drop_n_chain_basis
                  call chain_basis_element(ibasis, dv0, dv1, dv2, dv3, vown)
                  r4%vdir(:, pt%owner_idx) = vown
                  call fixed_direction_chain(self, slots%lsf(thread_slot)%lsf, pt, eff, &
                                             igrid, seeds, pw, r4%vdir, dv0, dv1, dv2, dv3, &
                                             0.0_wp, .false., context, sc, ft_work, field_row, &
                                             anchor_row, worker_error)
                  r4%vdir(:, pt%owner_idx) = 0.0_wp
                  if (allocated(worker_error)) then
                     call abort%latch_error(worker_error, igrid)
                     exit
                  end if
                  do i = 1, pt%n_active
                     r4%lmat(3*(i - 1) + 1:3*i, ibasis) = field_row(:, i)
                  end do
                  r4%lmat(nrow - 2:nrow, ibasis) = anchor_row
               end do
               ! Latched element abandons the whole point, response channel included
               if (abort%requested) cycle

               !* --------------- Branch element of a multi-branch point --------------- *!
               ! Unit branch-weight tangent, jet and owner at rest: coefficient of
               ! `d(wbranch)` in every column of the point
               ! Closed per group after the loop by [[branch_weight_block]]
               islot = 0
               if (allocated(br_slot)) islot = br_slot(igrid)
               if (islot > 0) then
                  dv0 = 0.0_wp
                  dv1 = 0.0_wp
                  dv2 = 0.0_wp
                  dv3 = 0.0_wp
                  call fixed_direction_chain(self, slots%lsf(thread_slot)%lsf, pt, eff, &
                                             igrid, seeds, pw, r4%vdir, dv0, dv1, dv2, dv3, &
                                             1.0_wp, .false., context, sc, ft_work, field_row, &
                                             anchor_row, worker_error)
                  if (allocated(worker_error)) then
                     call abort%latch_error(worker_error, igrid)
                     cycle
                  end if
                  call fill_branch_point(pt, seeds, ft_work, field_row, anchor_row, &
                                         br_pts(islot))
               end if

               !* ----------------- Coordinates of every unit direction ------------------ *!
               ! Active atom: packed column of the point's tensors
               ! Owner outside the active set: its displacement alone
               do idir = 1, 3*r4%ndir_atom
                  dir_loc = (idir - 1)/3 + 1
                  dir_axis = mod(idir - 1, 3) + 1
                  if (dir_loc <= pt%n_active) then
                     call drop_field_jet_column_packed(ft_work, pt%n_active, dir_axis, dir_loc, &
                                                       r4%cmat(1:drop_n_jet_coef, idir))
                  else
                     r4%cmat(1:drop_n_jet_coef, idir) = 0.0_wp
                  end if
                  r4%cmat(drop_n_jet_coef + 1:drop_n_chain_basis, idir) = 0.0_wp
                  if (dir_loc == r4%owner_loc) then
                     r4%cmat(drop_n_jet_coef + dir_axis, idir) = 1.0_wp
                  end if
               end do

               !* ----------------------- All columns of the point ----------------------- *!
               ! `L C` plus the explicit block on the active columns
               call gemm(r4%lmat, r4%cmat, r4%oblk)
               if (pt%n_active > 0) then
                  r4%oblk(1:3*pt%n_active, 1:3*pt%n_active) = &
                     r4%oblk(1:3*pt%n_active, 1:3*pt%n_active) + r4%mblk
               end if

               ! Anchor entries to the owner's row, field rows to the active atoms'
               ! rows, both in column `(dir_axis, dir_atom)`
               do idir = 1, 3*r4%ndir_atom
                  dir_loc = (idir - 1)/3 + 1
                  dir_axis = mod(idir - 1, 3) + 1
                  ient = r4%pair_ent(r4%owner_loc, dir_loc)
                  do iaxis = 1, 3
                     hess_threads(thread_slot)%blocks(iaxis, dir_axis, ient) = &
                        hess_threads(thread_slot)%blocks(iaxis, dir_axis, ient) &
                        + r4%oblk(nrow - 3 + iaxis, idir)
                  end do
                  do i = 1, pt%n_active
                     ient = r4%pair_ent(i, dir_loc)
                     hess_threads(thread_slot)%blocks(:, dir_axis, ient) = &
                        hess_threads(thread_slot)%blocks(:, dir_axis, ient) &
                        + r4%oblk(3*(i - 1) + 1:3*i, idir)
                  end do
               end do

               !* ---------------------- iSwiG switching channel ---------------------- *!
               ! One block over the influence set, no loop over directions
               ! Pair shared with the field channel: same entry, lands after it
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
               !* ---------- Explicit nuclear motion of the block's directions ---------- *!
               ! Weighted block applied to every direction of the block in one call
               if (use_apply .and. pt%n_active > 0) then
                  call slots%lsf(thread_slot)%lsf%vjp_f2_rArB_apply(pw%w_lsf0, pw%w_lsf1, &
                                                                       pw%w_lsf2, dirs(:, :, ilo:ihi), &
                                                                       xrow(:, :, 1:nblk))
               end if

               !* ----------------- One chain per supplied direction ------------------ *!
               do jdir = 1, nblk
                  idir_g = ilo + jdir - 1

                  ! Direction gathered onto the active slots; the owner's component
                  ! enters through the chain
                  do i = 1, pt%n_active
                     v_act(:, i) = dirs(:, pt%active_idx(i), idir_g)
                  end do
                  call drop_field_jet_tangent(ft_work, pt%n_active, v_act, dv0, dv1, dv2, dv3)
                  ! Host's partial tangents on top, third order included; the whole
                  ! tangent for a host-defined level set
                  if (have_jets) call add_host_jet(host_jets(:, igrid, jdir), dv0, dv1, dv2, dv3)

                  ! Branch-weight tangent of the direction, zero off a multi-branch group
                  dwb_dir = 0.0_wp
                  if (branched) dwb_dir = dwb_blk(igrid, jdir)

                  call fixed_direction_chain(self, slots%lsf(thread_slot)%lsf, pt, eff, &
                                             igrid, seeds, pw, dirs(:, :, idir_g), &
                                             dv0, dv1, dv2, dv3, dwb_dir, .not. use_apply, &
                                             context, sc, ft_work, field_row, anchor_row, &
                                             worker_error)
                  if (allocated(worker_error)) then
                     call abort%latch_error(worker_error, igrid)
                     exit
                  end if
                  if (use_apply) then
                     do i = 1, pt%n_active
                        field_row(:, i) = field_row(:, i) + xrow(:, i, jdir)
                     end do
                  end if
                  ! Fixed half of the level-set weight tangent, Hessian entry symmetrised
                  if (want_rt_lsf) then
                     rt%dw_value(igrid, idir_g) = rt%dw_value(igrid, idir_g) + sc%dw_lsf0
                     rt%dw_gradient(:, igrid, idir_g) = rt%dw_gradient(:, igrid, idir_g) &
                                                        + sc%dw_lsf1
                     rt%dw_hessian(:, :, igrid, idir_g) = rt%dw_hessian(:, :, igrid, idir_g) &
                                                          + 0.5_wp*(sc%dw_lsf2 + transpose(sc%dw_lsf2))
                  end if

                  do i = 1, pt%n_active
                     iatom = pt%active_idx(i)
                     hvp_threads(:, iatom, jdir, thread_slot) = &
                        hvp_threads(:, iatom, jdir, thread_slot) + field_row(:, i)
                  end do
                  hvp_threads(:, pt%owner_idx, jdir, thread_slot) = &
                     hvp_threads(:, pt%owner_idx, jdir, thread_slot) + anchor_row

                  ! Switching block contracted with the direction over the influence set
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
                  ! Normal fold only with a model response, else hard zeros
                  if (have_omega) then
                     call seed_normal_channel(pt%state, deff(jdir), igrid, pt%state%lsf2_rr, &
                                              rw%w_lsf1, rw%w_xyz)
                  else
                     rw%w_xyz = w_xyz_zero
                  end if
                  call seed_jet_basis_contract(deff(jdir), igrid, pt%phi1_r, rw%w_xyz, &
                                               seeds%res(1:drop_n_jet_seeds), &
                                               seeds%x(:, 1:drop_n_jet_seeds), &
                                               rw%w_lsf0, rw%w_lsf1, rw%w_lsf2)
                  ! Response half of the level-set weight tangent
                  if (want_rt_lsf) then
                     idir_g = ilo + jdir - 1
                     rt%dw_value(igrid, idir_g) = rt%dw_value(igrid, idir_g) + rw%w_lsf0
                     rt%dw_gradient(:, igrid, idir_g) = rt%dw_gradient(:, igrid, idir_g) &
                                                        + rw%w_lsf1
                     rt%dw_hessian(:, :, igrid, idir_g) = rt%dw_hessian(:, :, igrid, idir_g) &
                                                          + 0.5_wp*(rw%w_lsf2 + transpose(rw%w_lsf2))
                  end if

                  !* ------------------------- Field channel -------------------------- *!
                  ! Row off the point's jet tensors, no LSF call per direction
                  call drop_field_jet_contract(ft_work, pt%n_active, rw%w_lsf0, rw%w_lsf1, &
                                               rw%w_lsf2, pt%vjp_pt)
                  do i = 1, pt%n_active
                     iatom = pt%active_idx(i)
                     hvp_threads(:, iatom, jdir, thread_slot) = &
                        hvp_threads(:, iatom, jdir, thread_slot) + pt%vjp_pt(:, i)
                  end do

                  !* --------------------- Anchor channel (owner) --------------------- *!
                  call seed_anchor_contract(deff(jdir), igrid, pt%phi1_r, rw%w_xyz, &
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
            if (allocated(r4%lmat)) deallocate (r4%lmat, r4%cmat, r4%oblk)
         end if
         if (fixed_per_dir_here) then
            deallocate (v_act)
            if (use_apply) deallocate (xrow)
         end if
         call pt%destroy()
         !$omp end parallel

         if (abort%requested) then
            call abort%raise(context, error)
            if (allocated(hvp_entry)) hvp = hvp_entry
            call self%ctx%timer%stop(h_hess)
            return
         end if

         ! Deterministic reductions in fixed thread order
         ! Blocks are disjoint in the direction index: columns independent of blocking
         if (fixed_rank4_here) then
            do ithread = 1, slots%nthreads
               call hess_sparse_reduce(hess_threads(ithread), hess_fixed)
               call hess_sparse_destroy(hess_threads(ithread))
            end do
            ! Branch-weight motion of the multi-branch groups, serial, in grid order
            if (branched) call branch_weight_block(self, br_slot, br_pts, hess_fixed)
         end if
         if (want_cols) then
            do ithread = 1, slots%nthreads
               hvp(:, :, ilo:ihi) = hvp(:, :, ilo:ihi) + hvp_threads(:, :, 1:nblk, ithread)
            end do
         end if
         if (have_omega) hvp(:, :, ilo:ihi) = hvp(:, :, ilo:ihi) + hvp_direct(:, :, 1:nblk)

      end do

      call self%ctx%timer%stop(h_hess)

   end subroutine drop_hessian_traverse

   !* ================================================================================= *!
   !*                     Per-point pieces of the fixed channel                         *!
   !* ================================================================================= *!

   !> Apply the 16 basis seeds of one grid point
   !>
   !> - Single seed-application site of the traversal
   !> - Jet seeds, then anchor seeds, in the order of [[fill_seed_basis]]
   !> - `dstate_seed` written only with `want_dstate`, fixed channel only
   !>
   !> @param[in]     state        forward state of the point
   !> @param[in]     kkt          solved KKT sensitivities of the point
   !> @param[in]     want_dstate  request the derived-state tangents
   !> @param[out]    res_seed     linear response per seed
   !> @param[out]    seed_x       point motion and multiplier change per seed
   !> @param[in,out] dstate_seed  derived-state tangent per seed, untouched unless requested
   subroutine point_seed_basis(state, kkt, want_dstate, res_seed, seed_x, dstate_seed)
      !> Forward state of the grid point
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
   !> - Normal fold as on the gradient path, then the 13 jet-seed contributions
   !> - `normal_grad` and `nwn` come from [[seed_normal_channel]], never
   !>   recomputed here; zero when the channel is inactive
   !> - Jet seeds only: anchor seeds feed the owner's row in the chain,
   !>   [[scatter_jet_weight]] has no anchor case
   !>
   !> @param[in]  pt     point scratch, prologue run
   !> @param[in]  eff    folded surface adjoints, held fixed
   !> @param[in]  igrid  grid point
   !> @param[in]  seeds  applied basis seeds of the point
   !> @param[out] pw     weight set of the point
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

      pw%w_lsf0 = 0.0_wp
      pw%w_lsf1 = 0.0_wp
      pw%w_lsf2 = 0.0_wp
      call seed_normal_channel(pt%state, eff, igrid, pt%state%lsf2_rr, pw%w_lsf1, pw%w_xyz, &
                               normal_grad_pt=pw%normal_grad, nwn_pt=pw%nwn)
      call seed_jet_basis_contract(eff, igrid, pt%phi1_r, pw%w_xyz, &
                                   seeds%res(1:drop_n_jet_seeds), &
                                   seeds%x(:, 1:drop_n_jet_seeds), &
                                   pw%w_lsf0, pw%w_lsf1, pw%w_lsf2)
   end subroutine point_weights

   !> Build one element of the symmetrised basis of the chain's direction inputs
   !>
   !>     1        dv0 = 1
   !>     2 .. 4   dv1 = e_a
   !>     5 .. 10  dv2 = 1 at (a, b) and (b, a)
   !>     11 .. 20 dv3 = 1 at every permutation of (a, b, c)
   !>     21 .. 23 owner displacement e_a
   !>
   !> - Jet tangent `(dv0, dv1, dv2, dv3)` plus owner displacement, dimension
   !>   `1 + 3 + 6 + 10 + 3 = 23`
   !> - Order of [[drop_field_jet_column_packed]], `a <= b` and `a <= b <= c`,
   !>   then the three owner axes
   !> - Symmetric class holds a one at every index of the class, so its packed
   !>   coordinate is the coefficient
   !> - No element for the branch-weight tangent, closed per group instead
   !>
   !> @param[in]  ibasis  basis element, `1 .. drop_n_chain_basis`
   !> @param[out] dv0     jet-value component
   !> @param[out] dv1     jet-gradient component `(3)`
   !> @param[out] dv2     jet-Hessian component `(3, 3)`
   !> @param[out] dv3     jet third-derivative component `(3, 3, 3)`
   !> @param[out] vown    owner displacement `(3)`
   pure subroutine chain_basis_element(ibasis, dv0, dv1, dv2, dv3, vown)
      !> Basis element
      integer, intent(in) :: ibasis
      !> Jet components of the element
      real(wp), intent(out) :: dv0, dv1(3), dv2(3, 3), dv3(3, 3, 3)
      !> Owner displacement of the element
      real(wp), intent(out) :: vown(3)

      !> Symmetry class and its representative indices
      integer :: k, a, b, c

      dv0 = 0.0_wp
      dv1 = 0.0_wp
      dv2 = 0.0_wp
      dv3 = 0.0_wp
      vown = 0.0_wp

      if (ibasis == 1) then
         dv0 = 1.0_wp
      else if (ibasis <= 1 + ndim) then
         dv1(ibasis - 1) = 1.0_wp
      else if (ibasis <= 1 + ndim + drop_n_sym2) then
         k = ibasis - 1 - ndim
         a = drop_sym2_idx(1, k)
         b = drop_sym2_idx(2, k)
         dv2(a, b) = 1.0_wp
         dv2(b, a) = 1.0_wp
      else if (ibasis <= drop_n_jet_coef) then
         k = ibasis - 1 - ndim - drop_n_sym2
         a = drop_sym3_idx(1, k)
         b = drop_sym3_idx(2, k)
         c = drop_sym3_idx(3, k)
         ! Assignments, not increments: a repeated index names one entry twice
         dv3(a, b, c) = 1.0_wp
         dv3(a, c, b) = 1.0_wp
         dv3(b, a, c) = 1.0_wp
         dv3(b, c, a) = 1.0_wp
         dv3(c, a, b) = 1.0_wp
         dv3(c, b, a) = 1.0_wp
      else
         vown(ibasis - drop_n_jet_coef) = 1.0_wp
      end if
   end subroutine chain_basis_element

   !> Run the second-order chain of the fixed channel along one nuclear direction
   !>
   !> - Gradient column of `(dJ^T/dv) omega` at the point: field row over the
   !>   active slots and the three anchor entries of the owner
   !> - Switching channel not included, contracted by the caller
   !> - Called inside an OpenMP region: all state through arguments, no host
   !>   association
   !> - Input tangents are total, point motion included
   !> - `dw_lsf2` only against the symmetric `lsf3_rr_rA`, see the module header
   !> - Linear in `dwb` jointly with the jet tangent and the owner displacement
   !>
   !> @param[in]     self          DROP cavity instance, for the objective parameter
   !> @param[in]     lsf           level set of this thread, prepared at the point
   !> @param[in]     pt            point scratch, prologue run
   !> @param[in]     eff           folded surface adjoints, held fixed
   !> @param[in]     igrid         grid point
   !> @param[in]     seeds         applied basis seeds of the point, with tangents
   !> @param[in]     pw            weight set of the point
   !> @param[in]     v             nuclear direction `(3, nsph)`
   !> @param[in]     dv0           jet tangent of the value along `v` at the frozen point
   !> @param[in]     dv1           jet tangent of the gradient
   !> @param[in]     dv2           jet tangent of the Hessian
   !> @param[in]     dv3           jet tangent of the third derivative
   !> @param[in]     dwb           branch-weight tangent of the point along `v`
   !> @param[in]     explicit      include the `hvp_jet_rA` motion, else the caller adds it
   !> @param[in]     context       calling routine, prefixes the diagnostics
   !> @param[in,out] sc            per-direction temporaries
   !> @param[in,out] ft_work       jet tensors of the point, fills and fold run
   !> @param[in,out] field_row     field row `(3, >= n_active)`, active slots overwritten
   !> @param[out]    anchor_row    owner entries of the column
   !> @param[out]    worker_error  failure of the bordered solves, if any
   subroutine fixed_direction_chain(self, lsf, pt, eff, igrid, seeds, pw, v, &
                                    dv0, dv1, dv2, dv3, dwb, explicit, context, sc, ft_work, &
                                    field_row, anchor_row, worker_error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Level set of this thread, prepared at the point
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
      !> Tangent of the point's branch weight along `v`
      real(wp), intent(in) :: dwb
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
      !> Branch motion of the Lebedev weight
      real(wp) :: dwleb_branch
      !> Cartesian and seed indices
      integer :: iaxis, ibasis, k

      anchor_row = 0.0_wp

      !* ------------- Directional response of the projected point -------------- *!
      ! Same bordered system as the seeds, anchor rigid with its owner:
      ! `d^2 phi/(dr dR_owner) = -alpha I`
      sc%rhs_v = 0.0_wp
      sc%rhs_v(1:3, 1) = pt%state%lambda_val*dv1 + self%param%phi_alpha*v(:, pt%owner_idx)
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
      ! Anchor Lebedev weight belongs to the rigid sphere, zero tangent
      sc%dinp_v%danchor_wleb0 = 0.0_wp
      ! `apply_seed` froze the branch weight: its motion is added back to
      ! `d(wleb)` and `d(xi0)`, guards as in `drop_surface_tangent_core`
      ! Zero off a multi-branch anchor group
      sc%dinp_v%dwbranch = dwb
      if (dwb /= 0.0_wp .and. pt%state%wbranch > tiny(1.0_wp)) then
         dwleb_branch = (pt%state%wleb/pt%state%wbranch)*dwb
         sc%dinp_v%dwleb = sc%dinp_v%dwleb + dwleb_branch
         if (pt%state%wleb > seed_weight_tol) then
            sc%dinp_v%dxi0 = sc%dinp_v%dxi0 &
                             - 0.5_wp*pt%state%xi0*dwleb_branch/pt%state%wleb
         end if
      end if

      ! Tangent of `phi1_r = alpha (r* - anchor)`, anchor rigid with its owner;
      ! read by the branch term alone
      sc%dphi1_v = self%param%phi_alpha*(sc%dr_v - v(:, pt%owner_idx))

      !* --------------------- Tangent of the normal fold ----------------------- *!
      sc%dnormal_grad = 0.0_wp
      sc%dw_xyz = 0.0_wp
      if (eff%have_wn) then
         sc%dnormal_grad = (-sc%res_v%dn_surf*pw%nwn &
                            - pt%state%n_surf*dot_product(sc%res_v%dn_surf, eff%w_n(:, igrid))) &
                           /pt%state%g_norm &
                           - pw%normal_grad*sc%res_v%d_gnorm/pt%state%g_norm
         sc%dw_xyz = matmul(sc%res_v%dH, pw%normal_grad) &
                     + matmul(pt%state%lsf2_rr, sc%dnormal_grad)
      end if

      !* ---------------- Tangent of the seed right-hand sides ------------------ *!
      sc%dH_lag(:, :, 1) = -sc%dl_v*pt%state%lsf2_rr - pt%state%lambda_val*sc%res_v%dH
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
                                 seeds%x(1:3, ibasis), seeds%x(4, ibasis), &
                                 seed_dzero1, seed_dzero2, &
                                 sc%dseed_x(1:3, ibasis), sc%dseed_x(4, ibasis), &
                                 seeds%res(ibasis), seeds%dstate(ibasis), sc%dres)
         if (ibasis > drop_n_jet_seeds) then
            ! Anchor seed: feeds the owner's row, rigid-motion shift
            ! `phi1_r(iaxis)` moves with `phi1_r`
            iaxis = ibasis - drop_n_jet_seeds
            contribution = seed_contribution_tangent(eff, igrid, pw%w_xyz, sc%dw_xyz, &
                                                     seeds%x(1:3, ibasis), &
                                                     sc%dseed_x(1:3, ibasis), sc%dres, &
                                                     pt%phi1_r, sc%dphi1_v, sc%dphi1_v(iaxis))
            anchor_row(iaxis) = contribution
         else
            contribution = seed_contribution_tangent(eff, igrid, pw%w_xyz, sc%dw_xyz, &
                                                     seeds%x(1:3, ibasis), &
                                                     sc%dseed_x(1:3, ibasis), sc%dres, &
                                                     pt%phi1_r, sc%dphi1_v)
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
   !> - Pass 1: forward tangent, once for the whole block
   !> - Pass 2: weight tangent, once per direction, serial
   !> - Grid points are never a parallel axis: [[branch_phi_adj_tangent]]
   !>   reduces over whole anchor groups; directions may be
   !> - Pass 2b, with `omega_v`: model response folded at the base geometry and
   !>   added; curvature weights without `want_curvature` are refused
   !> - Without `omega_v`, `deff(idir)` carries `w_xi`, `w_f` and
   !>   `branch_phi_adj` only
   !> - `dirs` is one block; `deff` and the `(ngrid, ndir)` arrays are block sized
   !>
   !> @param[in]     self            DROP cavity instance
   !> @param[in]     acc             raw surface adjoints, held fixed
   !> @param[in]     eff             folded weights from [[prepare_surface_weights]]
   !> @param[in]     dirs            nuclear directions of one block `(3, nsph, ndir)`
   !> @param[in]     contracted      contracted accessor in pass 1; set by mode, not by block
   !> @param[in]     want_curvature  prologue carries the curvature invariants
   !> @param[in]     context         calling routine, prefixes the diagnostics
   !> @param[out]    deff            tangent of the folded weights per direction
   !> @param[out]    d_wbranch       branch-weight tangent `(ngrid, ndir)`
   !> @param[out]    error           allocated on failure
   !> @param[in,out] omega_v         surface-adjoint response of the model, optional
   !> @param[in,out] hvp_direct      non-surface columns, added to; required with `omega_v`
   !> @param[in,out] rt              response tangent of the whole direction set, optional
   !> @param[in]     first           global index of the first direction; required with `rt`
   !> @param[in]     host_jets       host's jet tangents `(40, ngrid, ndir)`, optional
   subroutine weight_tangents(self, acc, eff, dirs, contracted, want_curvature, context, &
                              deff, d_wbranch, error, omega_v, hvp_direct, rt, first, host_jets)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Raw surface adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Folded weights
      type(drop_surface_weights_type), intent(in) :: eff
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Whether pass 1 uses the contracted accessor per direction
      logical, intent(in) :: contracted
      !> Whether the prologue carries the curvature invariants
      logical, intent(in) :: want_curvature
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Tangent of the folded weights, one element per direction
      type(drop_surface_weights_type), allocatable, intent(out) :: deff(:)
      !> Directional tangent of the branch weight
      real(wp), allocatable, intent(out) :: d_wbranch(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Surface-adjoint response of the model
      class(surface_adjoint_response_type), intent(inout), optional :: omega_v
      !> Non-surface columns of the model response
      real(wp), intent(inout), optional :: hvp_direct(:, :, :)
      !> Response tangent of the whole direction set
      type(response_tangent_type), intent(inout), optional :: rt
      !> Global index of the block's first direction
      integer, intent(in), optional :: first
      !> Host's jet tangents of the block
      real(wp), intent(in), optional :: host_jets(:, :, :)

      !> Surface tangent, every channel
      type(cavity_surface_tangent_type) :: tangent
      !> Raw adjoint response of the model, one accumulator per direction
      type(cavity_surface_adjoint_type), allocatable :: dacc(:)
      !> Fold of one raw adjoint response
      type(drop_surface_weights_type) :: eff_v
      !> Branch bookkeeping, single-branch default when the cavity carries none
      integer, allocatable :: branch_count(:), anchor_id(:)
      !> Grid extent and direction index
      integer :: ngrid, ndir, idir

      ngrid = self%ngrid
      ndir = size(dirs, 3)

      !* -------------------------- Pass 1: forward tangent --------------------------- *!
      ! Curvature channels only for a model response whose primal adjoints
      ! carry curvature weights, the configuration of the contraction walk
      call tangent%init(ngrid, ndir, present(omega_v) .and. want_curvature)
      allocate (d_wbranch(ngrid, ndir))
      call drop_surface_tangent_core(self, dirs, tangent%have_curvature, tangent, error, &
                                     d_wbranch=d_wbranch, contracted=contracted, &
                                     host_jets=host_jets)
      if (allocated(error)) return

      ! Surface point motion for the host, every direction of the block
      if (present(rt)) then
         if (.not. present(first)) then
            call fatal_error(error, context//": the response tangent needs the block offset")
            return
         end if
         rt%xyz(:, :, first:first + ndir - 1) = tangent%d_xyz
      end if

      ! Single-branch stand-in reproduces the primal's early exit when the
      ! cavity holds no branch bookkeeping
      allocate (branch_count(ngrid), source=1)
      allocate (anchor_id(ngrid), source=0)
      if (allocated(self%branch_count)) branch_count = self%branch_count(1:ngrid)
      if (allocated(self%anchor_id)) anchor_id = self%anchor_id(1:ngrid)

      !* --------------------------- Pass 2: weight tangent --------------------------- *!
      allocate (deff(ndir))
      do idir = 1, ndir
         allocate (deff(idir)%w_xi(ngrid), deff(idir)%w_f(ngrid), &
                   deff(idir)%branch_phi_adj(ngrid))
         ! `have_wn` and `have_wk` stay `.false.`: normal and curvature adjoints
         ! are copies of the fixed raw channels, zero tangent
         call prepare_surface_weights_tangent(acc, eff, .true., &
                                              self%a(1:ngrid), self%wleb(1:ngrid), &
                                              self%xi0(1:ngrid), self%wbranch(1:ngrid), &
                                              self%radii, self%owner(1:ngrid), &
                                              branch_count, anchor_id, &
                                              self%branch_weight%s, &
                                              tangent%d_a(:, idir), tangent%d_w(:, idir), &
                                              tangent%d_xi(:, idir), d_wbranch(:, idir), &
                                              deff(idir)%w_xi, deff(idir)%w_f, &
                                              deff(idir)%branch_phi_adj, error)
         if (allocated(error)) return
      end do
      if (.not. present(omega_v)) return

      !* ---------------------- Pass 2b: model adjoint response ---------------------- *!
      if (.not. present(hvp_direct)) then
         call fatal_error(error, context//": a model adjoint response needs a buffer for"// &
                          " its non-surface columns")
         return
      end if
      allocate (dacc(ndir))
      do idir = 1, ndir
         call dacc(idir)%init(ngrid)
      end do
      call omega_v%apply(self, dirs, tangent, dacc, hvp_direct(:, :, 1:ndir), error, &
                         first=first, rt=rt)
      if (allocated(error)) return
      call tangent%destroy()

      ! Fold is linear in the raw adjoints: the response folded at the base
      ! geometry adds to `d(eff)`; copied channels come over whole
      do idir = 1, ndir
         call check_surface_adjoint(self, dacc(idir), context//" (adjoint response)", error)
         if (allocated(error)) return
         call prepare_surface_weights(self, dacc(idir), .true., eff_v)
         if (eff_v%have_wk .and. .not. want_curvature) then
            call fatal_error(error, context//": the model's adjoint response carries"// &
                             " curvature weights but its adjoints do not; the traversal"// &
                             " was configured without curvature")
            return
         end if
         deff(idir)%w_xi = deff(idir)%w_xi + eff_v%w_xi
         deff(idir)%w_f = deff(idir)%w_f + eff_v%w_f
         deff(idir)%branch_phi_adj = deff(idir)%branch_phi_adj + eff_v%branch_phi_adj
         call move_alloc(eff_v%w_xyz, deff(idir)%w_xyz)
         call move_alloc(eff_v%w_n, deff(idir)%w_n)
         call move_alloc(eff_v%w_k1, deff(idir)%w_k1)
         call move_alloc(eff_v%w_k2, deff(idir)%w_k2)
         deff(idir)%have_wn = eff_v%have_wn
         deff(idir)%have_wk = eff_v%have_wk
         ! Raw response released once folded, bounds the block's peak memory
         call dacc(idir)%destroy()
      end do

   end subroutine weight_tangents

   !* ================================================================================= *!
   !*                  Branch-weight motion of the fixed channel                        *!
   !* ================================================================================= *!

   !> Compute the branch-weight tangents of one block, fixed-only traversal
   !>
   !> - Forward tangent run for its `d(wbranch)` channel alone
   !> - Contracted, as in [[weight_tangents]] for the per-direction mode: same
   !>   tangent to the bit from both sources
   !>
   !> @param[in]  self       DROP cavity instance
   !> @param[in]  dirs       nuclear directions of one block `(3, nsph, ndir)`
   !> @param[out] d_wbranch  branch-weight tangent `(ngrid, ndir)`
   !> @param[out] error      allocated on failure
   subroutine branch_weight_tangents(self, dirs, d_wbranch, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Directional tangent of the branch weight
      real(wp), allocatable, intent(out) :: d_wbranch(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface tangent, filled but not read
      type(cavity_surface_tangent_type) :: tangent

      call tangent%init(self%ngrid, size(dirs, 3), .false.)
      allocate (d_wbranch(self%ngrid, size(dirs, 3)))
      call drop_surface_tangent_core(self, dirs, .false., tangent, error, &
                                     d_wbranch=d_wbranch, contracted=.true.)
      call tangent%destroy()
   end subroutine branch_weight_tangents

   !> Number the grid points of the multi-branch anchor groups
   !>
   !> - Membership from [[next_branch_group]], not from `branch_count > 1`:
   !>   only a group's head is gated on the count
   !> - Same group walk as [[branch_weight_block]]
   !>
   !> @param[in]  self     DROP cavity instance, branch bookkeeping allocated
   !> @param[out] br_slot  element of every group point in `br_pts`, zero elsewhere
   !> @param[out] br_pts   branch elements, one per group point, unfilled
   subroutine branch_point_slots(self, br_slot, br_pts)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Element of every group point
      integer, allocatable, intent(out) :: br_slot(:)
      !> Branch elements
      type(drop_branch_point_type), allocatable, intent(out) :: br_pts(:)

      !> Group bookkeeping and the running element count
      integer :: cursor, ifirst, ilast, igrid, nslot
      !> Whether a further anchor group exists
      logical :: found

      allocate (br_slot(self%ngrid), source=0)
      nslot = 0
      cursor = 1
      do
         call next_branch_group(self%branch_count(1:self%ngrid), &
                                self%anchor_id(1:self%ngrid), cursor, ifirst, ilast, found)
         if (.not. found) exit
         do igrid = ifirst, ilast
            nslot = nslot + 1
            br_slot(igrid) = nslot
         end do
      end do
      allocate (br_pts(nslot))
   end subroutine branch_point_slots

   !> Store the two point-local factors of the branch-weight motion
   !>
   !> - `lrow`: chain column on a unit branch-weight tangent
   !> - `prow`: nuclear gradient of `Phi = 0.5 alpha |r* - anchor|^2`
   !> - Jet seeds: `phi1_r . dr`, turned into a field row off the jet tensors
   !> - Anchor seeds: same less the rigid-motion shift on the owner, as in
   !>   [[seed_contribution]]
   !>
   !> @param[in]  pt          point scratch, prologue run
   !> @param[in]  seeds       applied basis seeds of the point
   !> @param[in]  ft_work     jet tensors of the point
   !> @param[in]  field_row   field row of the unit branch-weight chain
   !> @param[in]  anchor_row  owner entries of the unit branch-weight chain
   !> @param[out] bp          branch element of the point
   pure subroutine fill_branch_point(pt, seeds, ft_work, field_row, anchor_row, bp)
      !> Point scratch
      type(drop_point_scratch_type), intent(in) :: pt
      !> Applied basis seeds
      type(drop_point_seeds_type), intent(in) :: seeds
      !> Jet tensors of the point
      type(drop_field_tangent_work_type), intent(in) :: ft_work
      !> Column of the unit branch-weight chain
      real(wp), intent(in) :: field_row(:, :), anchor_row(3)
      !> Branch element of the point
      type(drop_branch_point_type), intent(out) :: bp

      !> Level-set adjoints of a unit branch adjoint
      real(wp) :: w0, w1(3), w2(3, 3)
      !> Row count, seed and Cartesian indices
      integer :: nrow, ibasis, iaxis

      nrow = pt%n_active + 1
      allocate (bp%atoms(nrow), bp%lrow(3, nrow), bp%prow(3, nrow))
      bp%atoms(1:pt%n_active) = pt%active_idx(1:pt%n_active)
      bp%atoms(nrow) = pt%owner_idx

      bp%lrow(:, 1:pt%n_active) = field_row(:, 1:pt%n_active)
      bp%lrow(:, nrow) = anchor_row

      w0 = 0.0_wp
      w1 = 0.0_wp
      w2 = 0.0_wp
      do ibasis = 1, drop_n_jet_seeds
         call scatter_jet_weight(ibasis, dot_product(pt%phi1_r, seeds%x(1:3, ibasis)), &
                                 w0, w1, w2)
      end do
      bp%prow = 0.0_wp
      call drop_field_jet_contract(ft_work, pt%n_active, w0, w1, w2, bp%prow)
      do iaxis = 1, drop_n_anchor_seeds
         bp%prow(iaxis, nrow) = &
            dot_product(pt%phi1_r, seeds%x(1:3, drop_n_jet_seeds + iaxis)) - pt%phi1_r(iaxis)
      end do
   end subroutine fill_branch_point

   !> Close the branch-weight motion of the rank-4 fixed channel over its groups
   !>
   !>     H(:, :, beta, B)  +=  sum_m  L_m  d(p_m)/dR_(beta, B)
   !>
   !> - `L_m`: chain of point `m` on a unit branch-weight tangent
   !> - `p_m`: softmax weight of point `m`, derivative from
   !>   [[branch_weight_type:weights_grad]]
   !> - One parameter per Cartesian component of the atoms the group reaches,
   !>   the union of the members' active atoms and owners
   !> - Columns reach atoms outside a point's level set, hence the dense block
   !> - Serial, in grid order, after the thread reduction: same additions in
   !>   the same order whatever the thread count
   !>
   !> @param[in]     self        DROP cavity instance
   !> @param[in]     br_slot     element of every group point in `br_pts`
   !> @param[in]     br_pts      branch elements the grid loop filled
   !> @param[in,out] hess_fixed  rank-4 fixed-adjoint accumulator `(3, nsph, 3, nsph)`
   subroutine branch_weight_block(self, br_slot, br_pts, hess_fixed)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Element of every group point
      integer, intent(in) :: br_slot(:)
      !> Branch elements
      type(drop_branch_point_type), intent(in) :: br_pts(:)
      !> Rank-4 accumulator of the fixed channel
      real(wp), intent(inout) :: hess_fixed(:, :, :, :)

      !> Atoms a group reaches, and each atom's position among them
      integer, allocatable :: union_atoms(:), union_pos(:)
      !> Softmax scratch of one group, parameters first
      real(wp), allocatable :: phi(:), weights(:), dphi(:, :), dweights(:, :)
      !> Group bookkeeping
      integer :: cursor, ifirst, ilast, group_size, nunion, nbranch_max
      !> Grid point, group member, row, union atom, axis and parameter indices
      integer :: igrid, m_branch, irow, katom, kaxis, kparam
      !> Whether a further anchor group exists
      logical :: found

      nbranch_max = max_branch_group_size(self%branch_count(1:self%ngrid), &
                                          self%anchor_id(1:self%ngrid))
      allocate (phi(nbranch_max), weights(nbranch_max))
      allocate (union_atoms(self%nsph))
      allocate (union_pos(self%nsph), source=0)

      cursor = 1
      do
         call next_branch_group(self%branch_count(1:self%ngrid), &
                                self%anchor_id(1:self%ngrid), cursor, ifirst, ilast, found)
         if (.not. found) exit
         group_size = ilast - ifirst + 1

         ! Atoms the group's rows reach, in order of first appearance
         nunion = 0
         do igrid = ifirst, ilast
            associate (bp => br_pts(br_slot(igrid)))
               do irow = 1, size(bp%atoms)
                  if (union_pos(bp%atoms(irow)) == 0) then
                     nunion = nunion + 1
                     union_pos(bp%atoms(irow)) = nunion
                     union_atoms(nunion) = bp%atoms(irow)
                  end if
               end do
            end associate
         end do

         ! Objective gradients over that set; an owner that is active as well
         ! has two rows on the same atom, and they add
         allocate (dphi(ndim*nunion, group_size), source=0.0_wp)
         allocate (dweights(ndim*nunion, group_size))
         do igrid = ifirst, ilast
            m_branch = igrid - ifirst + 1
            phi(m_branch) = self%phi0(igrid)
            associate (bp => br_pts(br_slot(igrid)))
               do irow = 1, size(bp%atoms)
                  kparam = ndim*(union_pos(bp%atoms(irow)) - 1)
                  dphi(kparam + 1:kparam + ndim, m_branch) = &
                     dphi(kparam + 1:kparam + ndim, m_branch) + bp%prow(:, irow)
               end do
            end associate
         end do

         call self%branch_weight%weights_grad(phi(1:group_size), dphi, &
                                              weights=weights(1:group_size), &
                                              dweights=dweights)

         do igrid = ifirst, ilast
            m_branch = igrid - ifirst + 1
            associate (bp => br_pts(br_slot(igrid)))
               do katom = 1, nunion
                  do kaxis = 1, ndim
                     kparam = ndim*(katom - 1) + kaxis
                     do irow = 1, size(bp%atoms)
                        hess_fixed(:, bp%atoms(irow), kaxis, union_atoms(katom)) = &
                           hess_fixed(:, bp%atoms(irow), kaxis, union_atoms(katom)) &
                           + bp%lrow(:, irow)*dweights(kparam, m_branch)
                     end do
                  end do
               end do
            end associate
         end do

         union_pos(union_atoms(1:nunion)) = 0
         deallocate (dphi, dweights)
      end do
   end subroutine branch_weight_block

   !* ================================================================================= *!
   !*                       Sparse per-thread Hessian accumulator                       *!
   !* ================================================================================= *!

   !> Initialise one thread's empty accumulator
   !>
   !> - Call inside the parallel region on the owning thread, first touch
   !>
   !> @param[in,out] self  accumulator of one thread
   !> @param[in]     nsph  atoms the pair indices range over
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

      ! Halve the default table while the smaller one still holds `cap`
      ! entries under the load factor
      ntab = hess_init_tab
      do while (ntab > 4 .and. cap*hess_load_den < (ntab/2)*hess_load_num)
         ntab = ntab/2
      end do
      allocate (self%htab(ntab), source=0)

   end subroutine hess_sparse_new

   !> Release one thread's accumulator
   !>
   !> @param[in,out] self  accumulator of one thread
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
   !> - Indices fed as whole words, not byte by byte
   !> - Odd multiplier is a bijection modulo a power of two: the low bits
   !>   permute a contiguous run of atom ids
   !>
   !> @param[in] iatom  row atom
   !> @param[in] jatom  column atom
   pure function hess_pair_hash(iatom, jatom) result(h)
      !> Row and column atom
      integer, intent(in) :: iatom, jatom
      !> Hash value in `[0, 2^32)`
      integer(int64) :: h

      h = ieor(fnv_offset, iand(int(iatom, int64), mask32))
      h = iand(h*fnv_prime, mask32)
      h = ieor(h, iand(int(jatom, int64), mask32))
      h = iand(h*fnv_prime, mask32)

   end function hess_pair_hash

   !> Find the `(3, 3)` block entry of one atom pair, creating it if new
   !>
   !> - New entry zeroed here and nowhere else; a never-written entry
   !>   contributes an exact zero
   !> - No pair is created twice
   !> - Both make the merge one add per element, see [[drop_hess_sparse_type]]
   !>
   !> @param[in,out] self   accumulator of one thread
   !> @param[in]     iatom  row atom
   !> @param[in]     jatom  column atom
   !> @param[out]    ient   entry index of the pair
   subroutine hess_sparse_entry(self, iatom, jatom, ient)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self
      !> Row and column atom
      integer, intent(in) :: iatom, jatom
      !> Entry index of the pair
      integer, intent(out) :: ient

      integer :: slot

      slot = hess_probe(self, iatom, jatom)
      ient = self%htab(slot)
      if (ient /= 0) return

      call hess_grow_entries(self)
      self%nent = self%nent + 1
      ient = self%nent
      self%pair_i(ient) = iatom
      self%pair_j(ient) = jatom
      self%blocks(:, :, ient) = 0.0_wp

      ! Probed slot is valid on the current table only; growth reinserts every entry
      if (self%nent*hess_load_den >= size(self%htab)*hess_load_num) then
         call hess_grow_table(self)
      else
         self%htab(slot) = ient
      end if

   end subroutine hess_sparse_entry

   !> Double the entry arrays if the next entry would not fit
   !>
   !> - Copy and `move_alloc`: entry indices must not move, blocks copied verbatim
   !>
   !> @param[in,out] self  accumulator of one thread
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
         ! Unreachable: `maxent` counts every pair, each entered once
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
   !> - Only the pair-to-entry lookup is rebuilt; blocks, entry indices and
   !>   accumulation order untouched
   !>
   !> @param[in,out] self  accumulator of one thread
   subroutine hess_grow_table(self)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self

      integer, allocatable :: new_tab(:)
      integer :: ient, slot

      allocate (new_tab(size(self%htab)*2), source=0)
      call move_alloc(new_tab, self%htab)

      ! Keys are unique, the probe stops at an empty slot
      do ient = 1, self%nent
         slot = hess_probe(self, self%pair_i(ient), self%pair_j(ient))
         self%htab(slot) = ient
      end do

   end subroutine hess_grow_table

   !> Slot of one pair in the open-addressing table, by linear probing
   !>
   !> - Slot of the pair's entry, else the empty slot where it belongs on the
   !>   current table
   !> - Single wraparound rule, see [[drop_hess_sparse_type]]
   !>
   !> @param[in] self   accumulator
   !> @param[in] iatom  row atom
   !> @param[in] jatom  column atom
   pure function hess_probe(self, iatom, jatom) result(slot)
      !> Accumulator
      type(drop_hess_sparse_type), intent(in) :: self
      !> Row and column atom
      integer, intent(in) :: iatom, jatom
      !> Slot of the pair, or the empty slot it would take
      integer :: slot

      integer :: ient

      slot = int(iand(hess_pair_hash(iatom, jatom), &
                      int(size(self%htab) - 1, int64)), kind(slot)) + 1
      do
         ient = self%htab(slot)
         if (ient == 0) return
         if (self%pair_i(ient) == iatom .and. self%pair_j(ient) == jatom) return
         slot = slot + 1
         if (slot > size(self%htab)) slot = 1
      end do
   end function hess_probe

   !> Add one thread's accumulator to the caller's Hessian
   !>
   !> - One add per element this thread reached
   !> - Call in the fixed `1 .. nthreads` order, see [[drop_hess_sparse_type]]
   !>
   !> @param[in]     self     accumulator of one thread
   !> @param[in,out] hessian  nuclear-Hessian accumulator `(3, nsph, 3, nsph)`
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
