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
      & drop_surface_weights_type, apply_seed, seed_status_message, seed_contribution

   implicit none(type, external)
   private

   public :: drop_kkt_factor_type, seed_normal_channel, seed_jet_basis, seed_anchor
   public :: seed_jet_basis_apply, seed_jet_basis_contract
   public :: seed_anchor_apply, seed_anchor_contract
   public :: jet_seed_index

   !> Number of level-set jet directions: one value, three gradient, nine Hessian
   integer, parameter, public :: drop_n_jet_seeds = 13

   !> Number of anchor directions: the rigid motion of the owner sphere
   integer, parameter, public :: drop_n_anchor_seeds = 3

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
   !> @param[in]    kkt       Solved KKT sensitivities; columns 1-4 are the jet seeds
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
   !> @param[in]  kkt      Solved KKT sensitivities; columns 1-4 are the jet seeds
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
            seed_x(:, ibasis) = kkt(1:4, ibasis)
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
   !> @param[in]    kkt        Solved KKT sensitivities; columns 5-7 are the anchor seeds
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
   !> @param[in]  kkt      Solved KKT sensitivities; columns 5-7 are the anchor seeds
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
         seed_x(:, iaxis) = kkt(1:4, 4 + iaxis)

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

end module moist_cavity_drop_derivatives_seeds
