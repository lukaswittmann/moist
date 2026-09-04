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
      & drop_surface_weights_type, apply_seed, seed_weight_tol, seed_status_message

   implicit none(type, external)
   private

   public :: drop_kkt_factor_type, seed_normal_channel, seed_jet_basis, seed_anchor
   public :: degenerate_point_error

   !> Number of level-set jet directions: one value, three gradient, nine Hessian
   integer, parameter, public :: drop_n_jet_seeds = 13

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
   !> @param[out] error        Error object, allocated when the system is singular
   subroutine drop_kkt_factor(self, H_lagrangian, lsf1_r, error)
      !> Factorization
      class(drop_kkt_factor_type), intent(out) :: self
      !> Lagrangian Hessian
      real(wp), intent(in) :: H_lagrangian(3, 3)
      !> Level-set gradient
      real(wp), intent(in) :: lsf1_r(3)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> LAPACK status
      integer(lapack_ik) :: info
      !> Rendered status; fixed length, so the message build is thread safe
      character(len=32) :: status

      self%lu = 0.0_wp
      self%lu(1:3, 1:3) = H_lagrangian
      self%lu(1:3, 4) = -lsf1_r
      self%lu(4, 1:3) = lsf1_r

      call lapack_getrf(4_lapack_ik, 4_lapack_ik, self%lu, 4_lapack_ik, self%ipiv, info)
      if (info /= 0_lapack_ik) then
         write (status, "(i0)") info
         call fatal_error(error, "Bordered KKT sensitivity matrix is singular"// &
                          " (getrf status "//trim(status)//")")
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
   !> @param[in]    self  Factorization
   !> @param[inout] rhs   `(4, nrhs)`; right-hand sides in, solutions out
   !> @param[out]   error Error object, allocated when LAPACK rejects the call
   subroutine drop_kkt_apply(self, rhs, error)
      !> Factorization
      class(drop_kkt_factor_type), intent(in) :: self
      !> Right-hand sides in, solutions out
      real(wp), contiguous, intent(inout) :: rhs(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> LAPACK status
      integer(lapack_ik) :: info
      !> Rendered status; fixed length, so the message build is thread safe
      character(len=32) :: status

      call lapack_getrs("n", 4_lapack_ik, int(size(rhs, 2), lapack_ik), self%lu, &
                        4_lapack_ik, self%ipiv, rhs, 4_lapack_ik, info)
      if (info /= 0_lapack_ik) then
         write (status, "(i0)") info
         call fatal_error(error, "Bordered KKT sensitivity solve failed"// &
                          " (getrs status "//trim(status)//")")
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
   !> @param[out]   error         Error object, allocated on inconsistent shapes
   subroutine drop_kkt_solve_tangent(self, dH_lagrangian, dlsf1_r, x, rhs, error)
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
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

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

      ndir = size(dH_lagrangian, 3)
      nseed = size(x, 2)

      if (size(dH_lagrangian, 1) /= 3 .or. size(dH_lagrangian, 2) /= 3 .or. &
          size(dlsf1_r, 1) /= 3 .or. size(dlsf1_r, 2) /= ndir .or. &
          size(x, 1) /= 4 .or. size(rhs, 1) /= 4 .or. &
          size(rhs, 2) /= nseed*ndir) then
         write (shapes, "(4(a, i0))") "nseed ", nseed, ", ndir ", ndir, &
            ", rhs ", size(rhs, 1), " x ", size(rhs, 2)
         call fatal_error(error, "Tangent KKT batch has inconsistent shapes ("// &
                          trim(shapes)//"; expected rhs 4 x nseed*ndir)")
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
         call fatal_error(error, "Tangent KKT sensitivity solve failed"// &
                          " (getrs status "//trim(status)//")")
      end if
   end subroutine drop_kkt_solve_tangent

   !> Fold an outward-normal adjoint into the level-set gradient and position channels
   !>
   !> The normal `n = grad S / |grad S|` depends on the level set twice over: at
   !> a fixed point through `grad S`, which lands on the gradient channel, and
   !> through the point's own motion, which lands on the position channel as
   !> `H @ normal_grad`.
   !>
   !> @param[in]    state     Per-grid point forward state
   !> @param[in]    eff       Folded surface adjoints
   !> @param[in]    igrid     Grid point
   !> @param[in]    lsf2_rr   Level-set Hessian at the projected point
   !> @param[inout] w_lsf1_pt Point-local level-set gradient adjoint
   !> @param[out]   w_xyz_pt  Effective position adjoint for every seed
   pure subroutine seed_normal_channel(state, eff, igrid, lsf2_rr, w_lsf1_pt, w_xyz_pt)
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
      real(wp) :: normal_grad(3)
      !> Normal component of the normal adjoint
      real(wp) :: nwn

      if (eff%have_wn) then
         nwn = dot_product(state%n_surf, eff%w_n(:, igrid))
         normal_grad = (eff%w_n(:, igrid) - state%n_surf*nwn)/state%g_norm
         w_lsf1_pt = w_lsf1_pt + normal_grad
         w_xyz_pt = eff%w_xyz(:, igrid) + matmul(lsf2_rr, normal_grad)
      else
         w_xyz_pt = eff%w_xyz(:, igrid)
      end if
   end subroutine seed_normal_channel

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

      !> Linear response to one seed
      type(drop_seed_result_type) :: res
      !> Seed perturbation of the level-set jet
      real(wp) :: dlsf1_r(3), dlsf2_rr(3, 3)
      !> Induced point motion and multiplier change
      real(wp) :: dr_dp(3), dlambda_dp
      !> Adjoint contribution of one seed
      real(wp) :: contribution
      !> Seed and Cartesian indices
      integer :: ibasis, iaxis, jaxis

      do ibasis = 1, drop_n_jet_seeds
         dlsf1_r = 0.0_wp
         dlsf2_rr = 0.0_wp
         iaxis = 0
         jaxis = 0
         if (ibasis <= 4) then
            if (ibasis > 1) dlsf1_r(ibasis - 1) = 1.0_wp
            dr_dp = kkt(1:3, ibasis)
            dlambda_dp = kkt(4, ibasis)
         else
            iaxis = (ibasis - 5)/3 + 1
            jaxis = mod(ibasis - 5, 3) + 1
            dlsf2_rr(iaxis, jaxis) = 1.0_wp
            dr_dp = 0.0_wp
            dlambda_dp = 0.0_wp
         end if

         call apply_seed(state, dlsf1_r, dlsf2_rr, dr_dp, dlambda_dp, res)

         ! No w_f term: the switching factor is an anchor-only iSwig overlap,
         ! so a level-set perturbation at fixed nuclei leaves it alone.
         contribution = dot_product(w_xyz_pt, dr_dp) + eff%w_xi(igrid)*res%dxi
         if (abs(eff%branch_phi_adj(igrid)) > seed_weight_tol) then
            contribution = contribution + eff%branch_phi_adj(igrid)*dot_product(phi1_r, dr_dp)
         end if
         if (eff%have_wk) then
            contribution = contribution + eff%w_k1(igrid)*res%dk1 + eff%w_k2(igrid)*res%dk2
         end if

         if (ibasis == 1) then
            w_lsf0_pt = w_lsf0_pt + contribution
         else if (ibasis <= 4) then
            w_lsf1_pt(ibasis - 1) = w_lsf1_pt(ibasis - 1) + contribution
         else
            w_lsf2_pt(iaxis, jaxis) = w_lsf2_pt(iaxis, jaxis) + contribution
         end if
      end do
   end subroutine seed_jet_basis

   !> Push the three anchor directions through the kernel
   !>
   !> The anchor moves rigidly with its owner sphere, so `da_i/dR_I = delta`.
   !> The level-set field is untouched; the whole channel enters through the
   !> objective, whose mixed derivative `-d^2 phi / dr dR_I` is `+alpha * I`.
   !> That makes it three extra right-hand sides on the same factorization,
   !> independent of the number of spheres.
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

      !> Linear response to one seed
      type(drop_seed_result_type) :: res
      !> Zero field perturbation: rigid anchor motion leaves the level set alone
      real(wp) :: dlsf1_r(3), dlsf2_rr(3, 3)
      !> Induced point motion and multiplier change
      real(wp) :: dr_dp(3), dlambda_dp
      !> Adjoint contribution of one seed
      real(wp) :: contribution
      !> Cartesian index
      integer :: iaxis

      dlsf1_r = 0.0_wp
      dlsf2_rr = 0.0_wp

      do iaxis = 1, 3
         dr_dp = kkt(1:3, 4 + iaxis)
         dlambda_dp = kkt(4, 4 + iaxis)

         call apply_seed(state, dlsf1_r, dlsf2_rr, dr_dp, dlambda_dp, res)

         contribution = dot_product(w_xyz_pt, dr_dp) + eff%w_xi(igrid)*res%dxi
         if (abs(eff%branch_phi_adj(igrid)) > seed_weight_tol) then
            ! phi = 0.5*alpha*|r - anchor|^2, so at fixed r the owner's rigid
            ! motion contributes -phi1_r on top of the point-motion term
            contribution = contribution + eff%branch_phi_adj(igrid) &
                           *(dot_product(phi1_r, dr_dp) - phi1_r(iaxis))
         end if
         if (eff%have_wk) then
            contribution = contribution + eff%w_k1(igrid)*res%dk1 + eff%w_k2(igrid)*res%dk2
         end if

         grad_owner(iaxis) = grad_owner(iaxis) + contribution
      end do
   end subroutine seed_anchor

   !> Report a degenerate projected point from the sensitivity kernel
   !>
   !> @param[in]  context Calling routine, used to prefix the diagnostic
   !> @param[in]  status  Status code returned by `build_seed_state`
   !> @param[in]  igrid   Grid point that failed
   !> @param[out] error   Error object, always allocated on return
   subroutine degenerate_point_error(context, status, igrid, error)
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Kernel status code
      integer, intent(in) :: status
      !> Failing grid point
      integer, intent(in) :: igrid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Rendered grid index
      character(len=32) :: idx

      write (idx, "(i0)") igrid
      call fatal_error(error, "[Error] "//context//": "//seed_status_message(status)// &
                       " at grid point "//trim(idx))
   end subroutine degenerate_point_error

end module moist_cavity_drop_derivatives_seeds
