!> Provides the seeds for the DROP derivatives
!>
!> * A seed is one unit perturbation direction pushed through the linear
!>   per-point map of [[moist_cavity_drop_derivatives_kernel]]
!> * The kernel yields the surface quantity for a given seed
!> * This module sends needed seeds through the kernel and collects
!>   their responses
!>
!> Two bases are seeded:
!> * the 13 level-set jet directions spanning `(S, grad S, grad^2 S)`
!> * the 3 anchor directions spanning the rigid motion of the owner sphere
!>
!> Because the projected point is an implicit minimizer, a seed of the field
!> also affects the point position; the bordered KKT solve acconts for that
!>
!> The outward-normal adjoint is folded into the position and gradient
!> channels (since all bases consume it)
!>
!> Both bases are offered twice over: as one routine that seeds and contracts
!> in a single sweep, and as an `_apply`/`_contract` pair. The split is not
!> cosmetic -- the expensive half, [[apply_seed]], reads the grid point alone,
!> while the surface adjoints reach only the cheap contraction. A caller that
!> contracts one point against many weight sets (the response half of the
!> Hessian does, once per nuclear direction) therefore pays the eigen
!> decomposition, tangent frame, Jacobian and curvature once instead of once
!> per set. The single-sweep routines are thin compositions of the two halves,
!> so there is exactly one copy of every floating-point chain
!>
!> The seed *layout* lives here too, and only here: [[jet_seed_index]] says
!> which jet direction a slot carries, [[seed_rhs_column]] which column of the
!> solved batch it occupies, and [[seed_standard_rhs]] emits that batch's
!> right-hand sides. The grid driver, the two `_apply` halves and
!> [[fill_seed_basis]] are all consumers, so no producer and no reader of a
!> column can drift away from the others
!>
!> Everything here is `self`-free and takes plain arguments, so it can be
!> called from a submodule of `moist_cavity_drop` and unit-tested without a
!> cavity
module moist_cavity_drop_derivatives_seeds
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use moist_math_lapack_getrf, only: lapack_getrf
   use moist_math_lapack_getrs, only: lapack_getrs
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, drop_seed_result_type, &
      & drop_seed_result_tangent_type, drop_surface_weights_type, apply_seed, &
      & seed_status_message, seed_contribution

   implicit none(type, external)
   private

   public :: drop_kkt_factor_type, seed_normal_channel, seed_jet_basis, seed_anchor
   public :: seed_jet_basis_apply, seed_jet_basis_contract
   public :: seed_anchor_apply, seed_anchor_contract
   public :: jet_seed_index, seed_rhs_column, seed_standard_rhs
   public :: fill_seed_basis, scatter_jet_weight, seed_contribution_tangent

   !> Number of level-set jet directions: one value, three gradient, nine Hessian
   integer, parameter, public :: drop_n_jet_seeds = 13

   !> Number of anchor directions: the rigid motion of the owner sphere
   integer, parameter, public :: drop_n_anchor_seeds = 3

   !> Seeds one grid point pushes through the kernel: the `drop_n_jet_seeds`
   !> level-set jet directions of [[seed_jet_basis]] followed by the
   !> `drop_n_anchor_seeds` anchor directions of [[seed_anchor]], in that order
   !>
   !> [[fill_seed_basis]] lays that order out and the second-order chain of
   !> `hessian_fixed.f90` walks the seeds by this index, so the two agree by
   !> construction rather than by two matching literals
   integer, parameter, public :: drop_n_point_seeds = drop_n_jet_seeds + drop_n_anchor_seeds

   !> Jet directions that move the projected point: the value direction and the
   !> three gradient directions
   !>
   !> The nine Hessian directions leave the stationarity conditions alone --
   !> those see only `S` and `grad S` -- so their right-hand side vanishes
   !> identically and they occupy no column of the seed batch
   integer, parameter, public :: drop_n_moving_jet_seeds = 1 + 3

   !> Columns of the standard seed batch [[seed_standard_rhs]] emits: the moving
   !> jet directions followed by the three anchor directions
   integer, parameter, public :: drop_n_seed_columns = drop_n_moving_jet_seeds &
                                                       + drop_n_anchor_seeds

   !> Vanishing tangent of a basis seed
   !>
   !> Every seed [[fill_seed_basis]] emits is a constant matrix, so its own
   !> derivative along a nuclear direction is zero and only the point motion it
   !> induces survives; these are the arguments that carries into
   !> [[apply_seed_tangent]]
   real(wp), parameter, public :: seed_dzero1(3) = 0.0_wp
   real(wp), parameter, public :: seed_dzero2(3, 3) = 0.0_wp

   !> Slot kinds of the jet-seed layout, as classified by [[jet_seed_index]]
   !>
   !> The layout of the `drop_n_jet_seeds` directions is owned here and nowhere
   !> else: every producer and every consumer of a jet slot goes through
   !> [[jet_seed_index]], so symmetrising the nine Hessian seeds is a change to
   !> one routine rather than to three re-derivations of `(ibasis - 5)/3 + 1`.
   !>
   !> Slot outside `1 .. drop_n_jet_seeds` -- an anchor seed, for instance
   integer, parameter, public :: drop_jet_seed_none = 0
   !> Slot 1: the level-set value direction
   integer, parameter, public :: drop_jet_seed_value = 1
   !> Slots 2-4: the three level-set gradient directions
   integer, parameter, public :: drop_jet_seed_grad = 2
   !> Slots 5-13: the nine level-set Hessian directions, in row-major order
   integer, parameter, public :: drop_jet_seed_hess = 3

   !> LU factorization of the bordered KKT sensitivity matrix
   !>
   !> With the Lagrangian Hessian `H_L = phi_rr - lambda S_rr` and the level-set
   !> gradient `g = S_r`,
   !>
   !> ```
   !>   [ H_L  -g ] [ dr/dp      ]   [ b_1:3 ]
   !>   [ g^T   0 ] [ dlambda/dp ] = [ b_4   ]
   !> ```
   !>
   !> The full bordered system is factored rather than eliminating `dlambda`,
   !> so `H_L` itself need not be invertible.
   !>
   !> Factorization is split from the solve because the Hessian passes need the
   !> factors to across solve calls at one grid point:
   !> the primal right-hand sides are known up front, the per-direction tangent
   !> ones (`K dx = db - dK x`) only once the primal solution exists; `solve`
   !> takes the first, `solve_tangent` the second. Components are
   !> fixed size, so an instance is stack-local inside the OpenMP grid loops and
   !> the path stays allocation-free -- the tangent batch buffer belongs to the
   !> caller for that same reason; only a failure allocates, and that through
   !> the error object all three routines report with.
   !>
   !> All three take a `context` -- the API entry point the user called -- and
   !> an optional `igrid`, and prefix their diagnostic with `context//": "`,
   !> matching [[check_surface_adjoint]]. A singular projection is a property of
   !> the geometry a user has to be able to trace back to a call and a grid
   !> point, and by the time the bordered matrix is assembled neither is
   !> recoverable from the arguments.
   type :: drop_kkt_factor_type
      !> LU factors of the bordered matrix, as returned by `getrf`
      real(wp) :: lu(4, 4)
      !> Pivot indices of the factorization
      integer(lapack_ik) :: ipiv(4)
   contains
      !> Assemble and factor the bordered matrix
      procedure :: factor => drop_kkt_factor
      !> Solve a right-hand side batch with the stored factors
      procedure :: solve => drop_kkt_apply
      !> Solve the tangent of a right-hand side batch with the same factors
      procedure :: solve_tangent => drop_kkt_solve_tangent
   end type drop_kkt_factor_type

contains

   !> Assemble and factor the bordered KKT sensitivity matrix
   !>
   !> A singular bordered matrix means a degenerate projected point, which is a
   !> condition of the geometry rather than a programming error. The LAPACK
   !> status is therefore turned into an error object here, so no caller has to
   !> know that a status code is what LAPACK returns.
   !>
   !> @param[out] self         Factorization
   !> @param[in]  H_lagrangian Lagrangian Hessian at the projected point
   !> @param[in]  lsf1_r       Level-set gradient at the projected point
   !> @param[in]  context      Calling routine, used to prefix the diagnostic
   !> @param[out] error        Error object, allocated when the system is singular
   !> @param[in]  igrid        Grid point being factored, when the caller has one
   subroutine drop_kkt_factor(self, H_lagrangian, lsf1_r, context, error, igrid)
      !> Factorization
      class(drop_kkt_factor_type), intent(out) :: self
      !> Lagrangian Hessian
      real(wp), intent(in) :: H_lagrangian(3, 3)
      !> Level-set gradient
      real(wp), intent(in) :: lsf1_r(3)
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Failing grid point, when the caller has one
      integer, intent(in), optional :: igrid

      !> LAPACK status
      integer(lapack_ik) :: info
      !> Rendered status; fixed length, so the message build is thread safe
      character(len=32) :: status
      !> Rendered grid index and its message suffix; fixed length, for the same reason
      character(len=32) :: idx
      character(len=48) :: at_point

      self%lu = 0.0_wp
      self%lu(1:3, 1:3) = H_lagrangian
      self%lu(1:3, 4) = -lsf1_r
      self%lu(4, 1:3) = lsf1_r

      call lapack_getrf(4_lapack_ik, 4_lapack_ik, self%lu, 4_lapack_ik, self%ipiv, info)
      if (info /= 0_lapack_ik) then
         write (status, "(i0)") info
         at_point = ""
         if (present(igrid)) then
            write (idx, "(i0)") igrid
            at_point = " at grid point "//trim(idx)
         end if
         call fatal_error(error, context//": Bordered KKT sensitivity matrix is singular"// &
                          " (getrf status "//trim(status)//")"//trim(at_point))
      end if
   end subroutine drop_kkt_factor

   !> Solve a right-hand side batch with the stored factors
   !>
   !> Requires a successful [[drop_kkt_factor]]. Once the factors exist `getrs`
   !> has no failure mode of its own: only an illegal argument sets its status,
   !> and that is a programming error rather than a degenerate point. So this
   !> path cannot be reached from a correct caller and cannot be exercised by a
   !> test; it reports for consistency with the factorization, not because a
   !> caller is expected to meet it.
   !>
   !> @param[in]    self    Factorization
   !> @param[inout] rhs     `(4, nrhs)`; right-hand sides in, solutions out
   !> @param[in]    context Calling routine, used to prefix the diagnostic
   !> @param[out]   error   Error object, allocated when LAPACK rejects the call
   !> @param[in]    igrid   Grid point being solved, when the caller has one
   subroutine drop_kkt_apply(self, rhs, context, error, igrid)
      !> Factorization
      class(drop_kkt_factor_type), intent(in) :: self
      !> Right-hand sides in, solutions out
      real(wp), contiguous, intent(inout) :: rhs(:, :)
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Failing grid point, when the caller has one
      integer, intent(in), optional :: igrid

      !> LAPACK status
      integer(lapack_ik) :: info
      !> Rendered status; fixed length, so the message build is thread safe
      character(len=32) :: status
      !> Rendered grid index and its message suffix; fixed length, for the same reason
      character(len=32) :: idx
      character(len=48) :: at_point

      call lapack_getrs("n", 4_lapack_ik, int(size(rhs, 2), lapack_ik), self%lu, &
                        4_lapack_ik, self%ipiv, rhs, 4_lapack_ik, info)
      if (info /= 0_lapack_ik) then
         write (status, "(i0)") info
         at_point = ""
         if (present(igrid)) then
            write (idx, "(i0)") igrid
            at_point = " at grid point "//trim(idx)
         end if
         call fatal_error(error, context//": Bordered KKT sensitivity solve failed"// &
                          " (getrs status "//trim(status)//")"//trim(at_point))
      end if
   end subroutine drop_kkt_apply

   !> Solve the tangent of a right-hand side batch with the same factors
   !>
   !> Differentiating `K x = b` along a perturbation gives `K dx = db - dK x`,
   !> where `dK` carries the same bordered layout [[drop_kkt_factor]] assembles:
   !>
   !> ```
   !>   dK = [ dH_L  -dg ]
   !>        [ dg^T    0 ]
   !> ```
   !>
   !> so the two `dlsf1_r` terms enter with *opposite* signs -- the multiplier
   !> term on rows 1-3 adds, the position term on row 4 subtracts.
   !>
   !> `H_L = phi_rr - lambda S_rr` and the objective Hessian `phi_rr = alpha*I`
   !> is constant in both the point and the anchor, so `d(phi_rr)` vanishes
   !> identically and the caller passes `dH_lagrangian = -dlambda S_rr -
   !> lambda d(S_rr)`, with nothing from the objective.
   !>
   !> Batched over directions deliberately: `dK` varies with the direction while
   !> `x` varies with the seed, so one grid point carries `nseed*ndir` tangent
   !> right-hand sides that all share the single 4x4 factorization. Forming them
   !> all and issuing one `getrs` costs a grid point one LAPACK call instead of
   !> `ndir` of them, on a matrix small enough that per-call overhead would
   !> otherwise dominate the solve. Column `(idir-1)*nseed + iseed`, so the
   !> seeds of one direction are contiguous; the driver depends on that.
   !>
   !> The batch is an argument rather than a local because `(4, nseed*ndir)` is
   !> not a fixed size: allocating it here would allocate once per grid point
   !> and lose the allocation-free path this type advertises. The caller
   !> allocates it once per thread outside the grid loop and reuses it for every
   !> point; this routine only checks that its shape is consistent.
   !>
   !> Requires a successful [[drop_kkt_factor]]. As in [[drop_kkt_apply]], once
   !> the factors exist `getrs` has no failure mode a correct caller can reach,
   !> so that status is reported for consistency rather than because it is
   !> expected. The shape mismatch above it is a programming error in the
   !> caller, and unlike the `getrs` status it is reachable and tested.
   !>
   !> @param[in]    self          Factorization, from a prior `factor` call
   !> @param[in]    dH_lagrangian Tangent of the Lagrangian Hessian, `(3, 3, ndir)`
   !> @param[in]    dlsf1_r       Tangent of the level-set gradient, `(3, ndir)`
   !> @param[in]    x             Primal solutions `(4, nseed)`, as `solve` returns them
   !> @param[inout] rhs           `db` in, `dx` out; `(4, nseed*ndir)`
   !> @param[in]    context       Calling routine, used to prefix the diagnostic
   !> @param[out]   error         Error object, allocated on inconsistent shapes
   !> @param[in]    igrid         Grid point being solved, when the caller has one
   subroutine drop_kkt_solve_tangent(self, dH_lagrangian, dlsf1_r, x, rhs, context, &
                                     error, igrid)
      !> Factorization, from a prior `factor` call
      class(drop_kkt_factor_type), intent(in) :: self
      !> Tangent of the Lagrangian Hessian block, `(3, 3, ndir)`
      real(wp), contiguous, intent(in) :: dH_lagrangian(:, :, :)
      !> Tangent of the level-set gradient, `(3, ndir)`
      real(wp), contiguous, intent(in) :: dlsf1_r(:, :)
      !> Primal solutions `x`, `(4, nseed)`, as returned by `solve`; never modified
      real(wp), contiguous, intent(in) :: x(:, :)
      !> `db` in, `dx` out; `(4, nseed*ndir)`, column `(idir-1)*nseed + iseed`
      real(wp), contiguous, intent(inout) :: rhs(:, :)
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Failing grid point, when the caller has one
      integer, intent(in), optional :: igrid

      !> Batch extents, taken from the tangent and primal arguments
      integer :: ndir, nseed
      !> Direction, seed, column and Cartesian indices
      integer :: idir, iseed, icol, iaxis
      !> LAPACK status
      integer(lapack_ik) :: info
      !> Rendered status; fixed length, so the message build is thread safe
      character(len=32) :: status
      !> Rendered shapes; fixed length, for the same reason
      character(len=64) :: shapes
      !> Rendered grid index and its message suffix; fixed length, for the same reason
      character(len=32) :: idx
      character(len=48) :: at_point

      ndir = size(dH_lagrangian, 3)
      nseed = size(x, 2)

      if (size(dH_lagrangian, 1) /= 3 .or. size(dH_lagrangian, 2) /= 3 .or. &
          size(dlsf1_r, 1) /= 3 .or. size(dlsf1_r, 2) /= ndir .or. &
          size(x, 1) /= 4 .or. size(rhs, 1) /= 4 .or. &
          size(rhs, 2) /= nseed*ndir) then
         write (shapes, "(4(a, i0))") "nseed ", nseed, ", ndir ", ndir, &
            ", rhs ", size(rhs, 1), " x ", size(rhs, 2)
         at_point = ""
         if (present(igrid)) then
            write (idx, "(i0)") igrid
            at_point = " at grid point "//trim(idx)
         end if
         call fatal_error(error, context//": Tangent KKT batch has inconsistent shapes ("// &
                          trim(shapes)//"; expected rhs 4 x nseed*ndir)"//trim(at_point))
         return
      end if

      do idir = 1, ndir
         do iseed = 1, nseed
            icol = (idir - 1)*nseed + iseed
            do iaxis = 1, 3
               ! Rows 1-3: the -dg block sits in column 4, so the multiplier
               ! term comes back into the right-hand side with a plus sign.
               rhs(iaxis, icol) = rhs(iaxis, icol) &
                                  - dH_lagrangian(iaxis, 1, idir)*x(1, iseed) &
                                  - dH_lagrangian(iaxis, 2, idir)*x(2, iseed) &
                                  - dH_lagrangian(iaxis, 3, idir)*x(3, iseed) &
                                  + dlsf1_r(iaxis, idir)*x(4, iseed)
            end do
            ! Row 4 carries +dg^T, so its term subtracts.
            rhs(4, icol) = rhs(4, icol) &
                           - dlsf1_r(1, idir)*x(1, iseed) &
                           - dlsf1_r(2, idir)*x(2, iseed) &
                           - dlsf1_r(3, idir)*x(3, iseed)
         end do
      end do

      call lapack_getrs("n", 4_lapack_ik, int(size(rhs, 2), lapack_ik), self%lu, &
                        4_lapack_ik, self%ipiv, rhs, 4_lapack_ik, info)
      if (info /= 0_lapack_ik) then
         write (status, "(i0)") info
         at_point = ""
         if (present(igrid)) then
            write (idx, "(i0)") igrid
            at_point = " at grid point "//trim(idx)
         end if
         call fatal_error(error, context//": Tangent KKT sensitivity solve failed"// &
                          " (getrs status "//trim(status)//")"//trim(at_point))
      end if
   end subroutine drop_kkt_solve_tangent

   !> Fold an outward-normal adjoint into the level-set gradient and position channels
   !>
   !> The normal `n = grad S / |grad S|` depends on the level set twice over: at
   !> a fixed point through `grad S`, which lands on the gradient channel, and
   !> through the point's own motion, which lands on the position channel as
   !> `H @ normal_grad`.
   !>
   !> The two intermediates are handed back on request, because a second-order
   !> caller needs their own tangents and must not rebuild them: a duplicated
   !> floating-point chain contracts differently under `-ffp-contract=fast`, so
   !> this routine is the single author of both. They are zero whenever the
   !> outward-normal channel is inactive, matching the fold itself.
   !>
   !> @param[in]    state          Per-grid point forward state
   !> @param[in]    eff            Folded surface adjoints
   !> @param[in]    igrid          Grid point
   !> @param[in]    lsf2_rr        Level-set Hessian at the projected point
   !> @param[inout] w_lsf1_pt      Point-local level-set gradient adjoint
   !> @param[out]   w_xyz_pt       Effective position adjoint for every seed
   !> @param[out]   normal_grad_pt Tangential normal adjoint over |grad S|, optional
   !> @param[out]   nwn_pt         Normal component of the normal adjoint, optional
   pure subroutine seed_normal_channel(state, eff, igrid, lsf2_rr, w_lsf1_pt, w_xyz_pt, &
                                       normal_grad_pt, nwn_pt)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Level-set Hessian
      real(wp), intent(in) :: lsf2_rr(3, 3)
      !> Point-local level-set gradient adjoint
      real(wp), intent(inout) :: w_lsf1_pt(3)
      !> Effective position adjoint
      real(wp), intent(out) :: w_xyz_pt(3)
      !> Tangential part of the normal adjoint, divided by |grad S|
      real(wp), intent(out), optional :: normal_grad_pt(3)
      !> Normal component of the normal adjoint
      real(wp), intent(out), optional :: nwn_pt

      !> Tangential part of the normal adjoint, divided by |grad S|
      real(wp) :: normal_grad(3)
      !> Normal component of the normal adjoint
      real(wp) :: nwn

      normal_grad = 0.0_wp
      nwn = 0.0_wp
      if (eff%have_wn) then
         nwn = dot_product(state%n_surf, eff%w_n(:, igrid))
         normal_grad = (eff%w_n(:, igrid) - state%n_surf*nwn)/state%g_norm
         w_lsf1_pt = w_lsf1_pt + normal_grad
         w_xyz_pt = eff%w_xyz(:, igrid) + matmul(lsf2_rr, normal_grad)
      else
         w_xyz_pt = eff%w_xyz(:, igrid)
      end if

      if (present(normal_grad_pt)) normal_grad_pt = normal_grad
      if (present(nwn_pt)) nwn_pt = nwn
   end subroutine seed_normal_channel

   !> Classify one seed slot of the level-set jet layout
   !>
   !> The `drop_n_jet_seeds` directions spanning `(S, grad S, grad^2 S)` are laid
   !> out as slot 1 for the value, slots 2-4 for the gradient components and
   !> slots 5-13 for the Hessian components in row-major order. This routine is
   !> the single author of that layout; [[seed_jet_basis]], [[fill_seed_basis]]
   !> and [[scatter_jet_weight]] are all consumers of it, so the emitting and
   !> the reading side of a slot cannot drift apart.
   !>
   !> `iaxis` and `jaxis` are always defined: the value slot and any slot
   !> outside the layout report zero for both, and a gradient slot reports its
   !> component in `iaxis` and zero in `jaxis`.
   !>
   !> @param[in]  ibasis    Seed slot
   !> @param[out] slot_kind One of the `drop_jet_seed_*` kinds
   !> @param[out] iaxis     First Cartesian index of the slot, or zero
   !> @param[out] jaxis     Second Cartesian index of the slot, or zero
   pure subroutine jet_seed_index(ibasis, slot_kind, iaxis, jaxis)
      !> Seed slot
      integer, intent(in) :: ibasis
      !> Kind of the slot
      integer, intent(out) :: slot_kind
      !> Cartesian indices of the slot
      integer, intent(out) :: iaxis, jaxis

      iaxis = 0
      jaxis = 0
      if (ibasis == 1) then
         slot_kind = drop_jet_seed_value
      else if (ibasis <= 4) then
         slot_kind = drop_jet_seed_grad
         iaxis = ibasis - 1
      else if (ibasis <= drop_n_jet_seeds) then
         slot_kind = drop_jet_seed_hess
         iaxis = (ibasis - 5)/3 + 1
         jaxis = mod(ibasis - 5, 3) + 1
      else
         slot_kind = drop_jet_seed_none
      end if
   end subroutine jet_seed_index

   !> Column of the standard seed batch one point seed's solution occupies
   !>
   !> [[seed_standard_rhs]] emits the batch and every consumer reads it back
   !> through this map, so the emitting and the reading side of a column cannot
   !> drift apart -- the same contract [[jet_seed_index]] enforces on the jet
   !> slots themselves, and what this routine is built on.
   !>
   !> The nine level-set Hessian directions have a vanishing right-hand side and
   !> therefore no column at all: they report zero, which a caller reading a seed
   !> by column has to take as "this seed moves neither the point nor the
   !> multiplier". A slot outside the layout reports zero for the same reason.
   !>
   !> @param[in] ibasis Point-seed slot, `1 .. drop_n_point_seeds`
   !> @return    icol   Batch column of the slot, or zero when it has none
   pure function seed_rhs_column(ibasis) result(icol)
      !> Point-seed slot
      integer, intent(in) :: ibasis
      !> Batch column of the slot
      integer :: icol

      !> Kind of the slot and its Cartesian indices
      integer :: slot_kind, iaxis, jaxis

      if (ibasis > drop_n_point_seeds) then
         icol = 0
         return
      end if
      if (ibasis > drop_n_jet_seeds) then
         icol = drop_n_moving_jet_seeds + (ibasis - drop_n_jet_seeds)
         return
      end if

      call jet_seed_index(ibasis, slot_kind, iaxis, jaxis)
      select case (slot_kind)
      case (drop_jet_seed_value)
         icol = 1
      case (drop_jet_seed_grad)
         icol = 1 + iaxis
      case default
         icol = 0
      end select
   end function seed_rhs_column

   !> Build the right-hand sides of the standard jet and anchor seed batch
   !>
   !> The batch every reverse traversal solves once per grid point, on the
   !> factorization [[drop_kkt_factor]] built there. It is written here rather
   !> than at the grid driver because the columns *are* the seed layout: which
   !> column carries which direction is [[seed_rhs_column]]'s to say, and
   !> [[fill_seed_basis]], [[seed_jet_basis_apply]] and [[seed_anchor_apply]]
   !> read the solutions back through that same map.
   !>
   !> The two live channels are the level-set jet and the anchor:
   !>
   !>   * perturbing `S` by one unit drives the stationarity condition `S = 0`
   !>     and lands on the bordered row alone;
   !>   * perturbing a component of `grad S` drives the multiplier term of the
   !>     Lagrangian and carries `lambda`;
   !>   * moving the owner sphere rigidly leaves the level-set field untouched
   !>     and reaches the system only through the objective's mixed derivative
   !>     `-d^2 phi/(dr dR_owner) = +alpha I`.
   !>
   !> @param[in]  lambda_val Lagrange multiplier of the projection
   !> @param[in]  phi_alpha  Quadratic coefficient of the projection objective
   !> @param[out] kkt_rhs    Right-hand sides, `(4, drop_n_seed_columns)`
   pure subroutine seed_standard_rhs(lambda_val, phi_alpha, kkt_rhs)
      !> Lagrange multiplier of the projection
      real(wp), intent(in) :: lambda_val
      !> Quadratic coefficient of the projection objective
      real(wp), intent(in) :: phi_alpha
      !> Right-hand sides of the batch
      real(wp), intent(out) :: kkt_rhs(4, drop_n_seed_columns)

      !> Seed slot, its column, its kind and its Cartesian indices
      integer :: ibasis, icol, slot_kind, iaxis, jaxis

      kkt_rhs = 0.0_wp

      do ibasis = 1, drop_n_jet_seeds
         icol = seed_rhs_column(ibasis)
         if (icol == 0) cycle
         call jet_seed_index(ibasis, slot_kind, iaxis, jaxis)
         if (slot_kind == drop_jet_seed_value) then
            kkt_rhs(4, icol) = -1.0_wp
         else
            kkt_rhs(iaxis, icol) = lambda_val
         end if
      end do

      do iaxis = 1, drop_n_anchor_seeds
         icol = seed_rhs_column(drop_n_jet_seeds + iaxis)
         kkt_rhs(iaxis, icol) = phi_alpha
      end do
   end subroutine seed_standard_rhs

   !> Lay out the 16 basis seeds in the order the second-order chain reads them
   !>
   !> Slots 1-13 are the level-set jet directions of [[seed_jet_basis]] -- the
   !> value, the three gradient components and the nine Hessian components, the
   !> last of which are single-entry matrices and therefore asymmetric -- and
   !> slots 14-16 the three anchor directions of [[seed_anchor]]. `x` collects
   !> the induced point motion and multiplier change of every seed in the
   !> `(4, nseed)` layout [[drop_kkt_solve_tangent]] expects; the nine Hessian
   !> seeds move neither, because the stationarity conditions see only `S` and
   !> `grad S`.
   !>
   !> Which jet slot carries which direction is not decided here:
   !> [[jet_seed_index]] owns that layout and [[seed_rhs_column]] owns the map
   !> onto the solved batch, so a change to the jet basis -- symmetrising the
   !> Hessian seeds, say -- does not have to be repeated in any consumer.
   !>
   !> @param[in]  kkt      Solved KKT sensitivities, as [[seed_standard_rhs]] laid them out
   !> @param[out] dlsf1_r  Gradient perturbation of each seed
   !> @param[out] dlsf2_rr Hessian perturbation of each seed
   !> @param[out] x        Induced point motion and multiplier change of each seed
   pure subroutine fill_seed_basis(kkt, dlsf1_r, dlsf2_rr, x)
      !> Solved KKT sensitivities
      real(wp), intent(in) :: kkt(:, :)
      !> Seed perturbations of the level-set jet
      real(wp), intent(out) :: dlsf1_r(3, drop_n_point_seeds)
      real(wp), intent(out) :: dlsf2_rr(3, 3, drop_n_point_seeds)
      !> Induced point motion and multiplier change
      real(wp), intent(out) :: x(4, drop_n_point_seeds)

      !> Seed index, its column, its kind and its Cartesian indices
      integer :: ibasis, icol, slot_kind, iaxis, jaxis

      dlsf1_r = 0.0_wp
      dlsf2_rr = 0.0_wp
      x = 0.0_wp

      do ibasis = 1, drop_n_point_seeds
         icol = seed_rhs_column(ibasis)
         if (icol > 0) x(:, ibasis) = kkt(1:4, icol)

         ! An anchor seed perturbs no level-set jet component at all: it reaches
         ! the system through the objective and is already fully described by
         ! the point motion its column carries.
         if (ibasis > drop_n_jet_seeds) cycle

         call jet_seed_index(ibasis, slot_kind, iaxis, jaxis)
         if (slot_kind == drop_jet_seed_hess) then
            dlsf2_rr(iaxis, jaxis, ibasis) = 1.0_wp
         else if (slot_kind == drop_jet_seed_grad) then
            dlsf1_r(iaxis, ibasis) = 1.0_wp
         end if
      end do
   end subroutine fill_seed_basis

   !> Place one jet seed's contribution in the level-set adjoint weights
   !>
   !> The same layout [[fill_seed_basis]] writes, read back through the same
   !> classifier: [[jet_seed_index]] says which of the value, gradient and
   !> Hessian channels a slot belongs to, and this routine only places the
   !> number. Anchor slots are not accepted -- they classify as
   !> `drop_jet_seed_none` and land nowhere -- because their contribution
   !> belongs to the owner's gradient row, not to a weight.
   !>
   !> @param[in]    ibasis Seed slot, `1 .. drop_n_jet_seeds`
   !> @param[in]    contrib Contribution of that seed
   !> @param[inout] w0     Level-set value adjoint
   !> @param[inout] w1     Level-set gradient adjoint
   !> @param[inout] w2     Level-set Hessian adjoint
   pure subroutine scatter_jet_weight(ibasis, contrib, w0, w1, w2)
      !> Seed slot
      integer, intent(in) :: ibasis
      !> Contribution of that seed
      real(wp), intent(in) :: contrib
      !> Level-set adjoints
      real(wp), intent(inout) :: w0, w1(3), w2(3, 3)

      !> Kind of the slot and its Cartesian indices
      integer :: slot_kind, iaxis, jaxis

      call jet_seed_index(ibasis, slot_kind, iaxis, jaxis)
      select case (slot_kind)
      case (drop_jet_seed_value)
         w0 = w0 + contrib
      case (drop_jet_seed_grad)
         w1(iaxis) = w1(iaxis) + contrib
      case (drop_jet_seed_hess)
         w2(iaxis, jaxis) = w2(iaxis, jaxis) + contrib
      end select
   end subroutine scatter_jet_weight

   !> Push the 13 level-set jet directions through the kernel
   !>
   !> The per-point map is linear in its seed, so seeding each basis direction
   !> of the jet `(S, grad S, grad^2 S)` once and collecting the responses gives
   !> the adjoint of the whole jet. Only the value and gradient directions move
   !> the point; the nine Hessian directions have `dr/dp = 0`.
   !>
   !> The results are point-local. `potential.f90` scatters them into
   !> `w_lsf*(..., igrid)`; `nuclear.f90` contracts them with `lsf*_rA`.
   !>
   !> A thin composition of [[seed_jet_basis_apply]] and
   !> [[seed_jet_basis_contract]], for the single-weight callers that have no
   !> reason to hold the responses. A caller that contracts the same point
   !> against several weight sets calls the two halves itself and pays the
   !> expensive one once; see [[get_surface_hessian_response_drop]].
   !>
   !> @param[in]    state     Per-grid point forward state
   !> @param[in]    eff       Folded surface adjoints
   !> @param[in]    igrid     Grid point
   !> @param[in]    phi1_r    Objective gradient at the projected point
   !> @param[in]    kkt       Solved KKT sensitivities, read through [[seed_rhs_column]]
   !> @param[in]    w_xyz_pt  Effective position adjoint from [[seed_normal_channel]]
   !> @param[inout] w_lsf0_pt Point-local level-set value adjoint
   !> @param[inout] w_lsf1_pt Point-local level-set gradient adjoint
   !> @param[inout] w_lsf2_pt Point-local level-set Hessian adjoint
   subroutine seed_jet_basis(state, eff, igrid, phi1_r, kkt, w_xyz_pt, &
                             w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Objective gradient
      real(wp), intent(in) :: phi1_r(3)
      !> Solved KKT sensitivities
      real(wp), intent(in) :: kkt(:, :)
      !> Effective position adjoint
      real(wp), intent(in) :: w_xyz_pt(3)
      !> Point-local level-set adjoints
      real(wp), intent(inout) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)

      !> Linear response of each seed, and its induced point motion
      type(drop_seed_result_type) :: res_seed(drop_n_jet_seeds)
      real(wp) :: seed_x(4, drop_n_jet_seeds)

      call seed_jet_basis_apply(state, kkt, res_seed, seed_x)
      call seed_jet_basis_contract(eff, igrid, phi1_r, w_xyz_pt, res_seed, seed_x, &
                                   w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)
   end subroutine seed_jet_basis

   !> Apply half of [[seed_jet_basis]]: the 13 responses of one grid point
   !>
   !> Everything here is a function of the point alone -- the seed directions
   !> are constant matrices, the induced motion comes from the point's own
   !> bordered solve, and [[apply_seed]] reads only `state`. No surface adjoint
   !> enters, which is exactly what lets a caller with several weight sets at
   !> one point run this once and the contraction many times.
   !>
   !> `seed_x` carries the induced point motion and multiplier change in the
   !> `(4, nseed)` layout the second-order chain uses, so the two halves of the
   !> DROP Hessian hold their seed batches the same way round.
   !>
   !> Which slot carries which direction is [[jet_seed_index]]'s to say, so this
   !> routine and [[seed_jet_basis_contract]] classify the slot the same way and
   !> agree by construction.
   !>
   !> @param[in]  state    Per-grid point forward state
   !> @param[in]  kkt      Solved KKT sensitivities, read through [[seed_rhs_column]]
   !> @param[out] res_seed Linear response of each seed, `(drop_n_jet_seeds)`
   !> @param[out] seed_x   Induced point motion and multiplier change, `(4, nseed)`
   subroutine seed_jet_basis_apply(state, kkt, res_seed, seed_x)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Solved KKT sensitivities
      real(wp), intent(in) :: kkt(:, :)
      !> Linear response of each seed
      type(drop_seed_result_type), intent(out) :: res_seed(drop_n_jet_seeds)
      !> Induced point motion and multiplier change of each seed
      real(wp), intent(out) :: seed_x(4, drop_n_jet_seeds)

      !> Seed perturbation of the level-set jet
      real(wp) :: dlsf1_r(3), dlsf2_rr(3, 3)
      !> Seed index, its kind and its Cartesian indices
      integer :: ibasis, slot_kind, iaxis, jaxis

      do ibasis = 1, drop_n_jet_seeds
         dlsf1_r = 0.0_wp
         dlsf2_rr = 0.0_wp
         call jet_seed_index(ibasis, slot_kind, iaxis, jaxis)
         if (slot_kind == drop_jet_seed_hess) then
            dlsf2_rr(iaxis, jaxis) = 1.0_wp
            seed_x(:, ibasis) = 0.0_wp
         else
            if (slot_kind == drop_jet_seed_grad) dlsf1_r(iaxis) = 1.0_wp
            seed_x(:, ibasis) = kkt(1:4, seed_rhs_column(ibasis))
         end if

         call apply_seed(state, dlsf1_r, dlsf2_rr, seed_x(1:3, ibasis), &
                         seed_x(4, ibasis), res_seed(ibasis))
      end do
   end subroutine seed_jet_basis_apply

   !> Contract half of [[seed_jet_basis]]: the 13 responses against one weight set
   !>
   !> The cheap half, and the only one a surface adjoint reaches. Every seed
   !> lands in its own slot of the level-set adjoints, so the loop carries no
   !> accumulation of its own and may be run once per weight set at a point
   !> whose [[seed_jet_basis_apply]] ran once.
   !>
   !> @param[in]    eff       Folded surface adjoints
   !> @param[in]    igrid     Grid point
   !> @param[in]    phi1_r    Objective gradient at the projected point
   !> @param[in]    w_xyz_pt  Effective position adjoint from [[seed_normal_channel]]
   !> @param[in]    res_seed  Linear response of each seed
   !> @param[in]    seed_x    Induced point motion and multiplier change of each seed
   !> @param[inout] w_lsf0_pt Point-local level-set value adjoint
   !> @param[inout] w_lsf1_pt Point-local level-set gradient adjoint
   !> @param[inout] w_lsf2_pt Point-local level-set Hessian adjoint
   subroutine seed_jet_basis_contract(eff, igrid, phi1_r, w_xyz_pt, res_seed, seed_x, &
                                      w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Objective gradient
      real(wp), intent(in) :: phi1_r(3)
      !> Effective position adjoint
      real(wp), intent(in) :: w_xyz_pt(3)
      !> Linear response of each seed
      type(drop_seed_result_type), intent(in) :: res_seed(drop_n_jet_seeds)
      !> Induced point motion and multiplier change of each seed
      real(wp), intent(in) :: seed_x(4, drop_n_jet_seeds)
      !> Point-local level-set adjoints
      real(wp), intent(inout) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)

      !> Adjoint contribution of one seed
      real(wp) :: contribution
      !> Seed index, its kind and its Cartesian indices
      integer :: ibasis, slot_kind, iaxis, jaxis

      do ibasis = 1, drop_n_jet_seeds
         call jet_seed_index(ibasis, slot_kind, iaxis, jaxis)

         ! A field seed leaves the anchor alone, so no rigid-motion shift; the
         ! switching factor is an anchor-only iSwiG overlap and is absent for
         ! the same reason. See [[seed_contribution]].
         contribution = seed_contribution(eff, igrid, w_xyz_pt, seed_x(1:3, ibasis), &
                                          res_seed(ibasis), phi1_r)

         select case (slot_kind)
         case (drop_jet_seed_value)
            w_lsf0_pt = w_lsf0_pt + contribution
         case (drop_jet_seed_grad)
            w_lsf1_pt(iaxis) = w_lsf1_pt(iaxis) + contribution
         case (drop_jet_seed_hess)
            w_lsf2_pt(iaxis, jaxis) = w_lsf2_pt(iaxis, jaxis) + contribution
         end select
      end do
   end subroutine seed_jet_basis_contract

   !> Push the three anchor directions through the kernel
   !>
   !> The anchor moves rigidly with its owner sphere, so `da_i/dR_I = delta`.
   !> The level-set field is untouched; the whole channel enters through the
   !> objective, whose mixed derivative `-d^2 phi / dr dR_I` is `+alpha * I`.
   !> That makes it three extra right-hand sides on the same factorization,
   !> independent of the number of spheres.
   !>
   !> Split the same way as [[seed_jet_basis]], and for the same reason.
   !>
   !> @param[in]    state      Per-grid point forward state
   !> @param[in]    eff        Folded surface adjoints
   !> @param[in]    igrid      Grid point
   !> @param[in]    phi1_r     Objective gradient at the projected point
   !> @param[in]    kkt        Solved KKT sensitivities, read through [[seed_rhs_column]]
   !> @param[in]    w_xyz_pt   Effective position adjoint from [[seed_normal_channel]]
   !> @param[inout] grad_owner Nuclear-gradient accumulator of the owner sphere
   subroutine seed_anchor(state, eff, igrid, phi1_r, kkt, w_xyz_pt, grad_owner)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Objective gradient
      real(wp), intent(in) :: phi1_r(3)
      !> Solved KKT sensitivities
      real(wp), intent(in) :: kkt(:, :)
      !> Effective position adjoint
      real(wp), intent(in) :: w_xyz_pt(3)
      !> Owner-sphere gradient accumulator
      real(wp), intent(inout) :: grad_owner(3)

      !> Linear response of each seed, and its induced point motion
      type(drop_seed_result_type) :: res_seed(drop_n_anchor_seeds)
      real(wp) :: seed_x(4, drop_n_anchor_seeds)

      call seed_anchor_apply(state, kkt, res_seed, seed_x)
      call seed_anchor_contract(eff, igrid, phi1_r, w_xyz_pt, res_seed, seed_x, &
                                grad_owner)
   end subroutine seed_anchor

   !> Apply half of [[seed_anchor]]: the three responses of one grid point
   !>
   !> @param[in]  state    Per-grid point forward state
   !> @param[in]  kkt      Solved KKT sensitivities, read through [[seed_rhs_column]]
   !> @param[out] res_seed Linear response of each seed, `(drop_n_anchor_seeds)`
   !> @param[out] seed_x   Induced point motion and multiplier change, `(4, nseed)`
   subroutine seed_anchor_apply(state, kkt, res_seed, seed_x)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Solved KKT sensitivities
      real(wp), intent(in) :: kkt(:, :)
      !> Linear response of each seed
      type(drop_seed_result_type), intent(out) :: res_seed(drop_n_anchor_seeds)
      !> Induced point motion and multiplier change of each seed
      real(wp), intent(out) :: seed_x(4, drop_n_anchor_seeds)

      !> Zero field perturbation: rigid anchor motion leaves the level set alone
      real(wp) :: dlsf1_r(3), dlsf2_rr(3, 3)
      !> Cartesian index
      integer :: iaxis

      dlsf1_r = 0.0_wp
      dlsf2_rr = 0.0_wp

      do iaxis = 1, drop_n_anchor_seeds
         seed_x(:, iaxis) = kkt(1:4, seed_rhs_column(drop_n_jet_seeds + iaxis))

         call apply_seed(state, dlsf1_r, dlsf2_rr, seed_x(1:3, iaxis), &
                         seed_x(4, iaxis), res_seed(iaxis))
      end do
   end subroutine seed_anchor_apply

   !> Contract half of [[seed_anchor]]: the three responses against one weight set
   !>
   !> @param[in]    eff        Folded surface adjoints
   !> @param[in]    igrid      Grid point
   !> @param[in]    phi1_r     Objective gradient at the projected point
   !> @param[in]    w_xyz_pt   Effective position adjoint from [[seed_normal_channel]]
   !> @param[in]    res_seed   Linear response of each seed
   !> @param[in]    seed_x     Induced point motion and multiplier change of each seed
   !> @param[inout] grad_owner Nuclear-gradient accumulator of the owner sphere
   subroutine seed_anchor_contract(eff, igrid, phi1_r, w_xyz_pt, res_seed, seed_x, &
                                   grad_owner)
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Objective gradient
      real(wp), intent(in) :: phi1_r(3)
      !> Effective position adjoint
      real(wp), intent(in) :: w_xyz_pt(3)
      !> Linear response of each seed
      type(drop_seed_result_type), intent(in) :: res_seed(drop_n_anchor_seeds)
      !> Induced point motion and multiplier change of each seed
      real(wp), intent(in) :: seed_x(4, drop_n_anchor_seeds)
      !> Owner-sphere gradient accumulator
      real(wp), intent(inout) :: grad_owner(3)

      !> Adjoint contribution of one seed
      real(wp) :: contribution
      !> Cartesian index
      integer :: iaxis

      do iaxis = 1, drop_n_anchor_seeds
         ! phi = 0.5*alpha*|r - anchor|^2, so at fixed r the owner's rigid
         ! motion contributes -phi1_r on top of the point-motion term; that is
         ! the shift [[seed_contribution]] takes.
         contribution = seed_contribution(eff, igrid, w_xyz_pt, seed_x(1:3, iaxis), &
                                          res_seed(iaxis), phi1_r, phi1_r(iaxis))

         grad_owner(iaxis) = grad_owner(iaxis) + contribution
      end do
   end subroutine seed_anchor_contract

   !> Directional derivative of [[seed_contribution]]
   !>
   !> Term by term the product rule applied to that contraction, with the
   !> surface adjoints themselves held fixed -- which is the whole premise of
   !> the fixed-adjoint half of the Hessian, this routine's one caller. The
   !> position adjoint still moves, because the normal fold inside it is built
   !> from the level-set gradient at the projected point.
   !>
   !> It lives beside its primal rather than in that caller so that the seed
   !> layout, its contraction and the contraction's tangent stay in one file.
   !> The branch term is absent because a grid carrying a multi-branch anchor
   !> group is refused at the public entry points in `hessian.f90` -- the only
   !> route to this routine's caller -- and the switching term for the reason
   !> [[seed_contribution]] gives.
   !>
   !> @param[in] eff       Folded surface adjoints
   !> @param[in] igrid     Grid point
   !> @param[in] w_xyz_pt  Effective position adjoint
   !> @param[in] dw_xyz_pt Tangent of the effective position adjoint
   !> @param[in] dr        Induced point motion of the seed
   !> @param[in] ddr       Tangent of that point motion
   !> @param[in] dres      Second-order response of the seed
   !> @return              Tangent of the adjoint contribution
   pure function seed_contribution_tangent(eff, igrid, w_xyz_pt, dw_xyz_pt, dr, ddr, dres) &
      result(contribution)
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Effective position adjoint and its tangent
      real(wp), intent(in) :: w_xyz_pt(3), dw_xyz_pt(3)
      !> Induced point motion and its tangent
      real(wp), intent(in) :: dr(3), ddr(3)
      !> Second-order response
      type(drop_seed_result_tangent_type), intent(in) :: dres
      !> Tangent of the adjoint contribution
      real(wp) :: contribution

      contribution = dot_product(dw_xyz_pt, dr) + dot_product(w_xyz_pt, ddr) &
                     + eff%w_xi(igrid)*dres%dxi
      if (eff%have_wk) then
         contribution = contribution + eff%w_k1(igrid)*dres%dk1 + eff%w_k2(igrid)*dres%dk2
      end if
   end function seed_contribution_tangent

end module moist_cavity_drop_derivatives_seeds
