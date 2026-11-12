!> Provides the seeds for the DROP derivative paths
!>
!> - A seed is one unit perturbation direction pushed through the linear
!>   per-point map of [[moist_cavity_drop_derivatives_kernel]]
!> - The kernel yields the surface observables for a given seed
!> - This module sends needed seeds through the kernel and collects
!>   their responses
!>
!> Two bases are seeded:
!> - the 10 level-set jet directions spanning `(S, grad S, grad^2 S)`, the
!>   Hessian part in symmetric packed form
!> - the 3 anchor directions spanning the rigid motion of the owner sphere
!>
!> Because the projected point is an implicit minimizer, a seed of the field
!> also affects the point position; the bordered KKT solve accounts for that
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
   use moist_math_lapack_gesv, only: lapack_gesv
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_packed, only: packed_contract2
   use moist_math_smoothstep, only: moist_math_smoothstep_type
   use moist_cavity_drop_phi, only: moist_cavity_drop_phi_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, drop_seed_result_type, &
      & drop_surface_weights_type, apply_seed, build_seed_state, seed_weight_tol, &
      & seed_status_message, seed_state_ok, seed_state_singular_kkt

   implicit none(type, external)
   private

   public :: drop_kkt_solve, seed_normal_channel, seed_jet_basis, seed_anchor
   public :: seed_point_jet
   public :: degenerate_point_error

   !> Number of level-set jet directions: one value, three gradient, six packed Hessian
   integer, parameter, public :: drop_n_jet_seeds = 10

contains

   !> Level-set jet adjoint of one projected surface point
   !>
   !> The per-point kernel of both reverse paths: objective jet, seed state,
   !> outward-normal fold, bordered KKT solve of the four jet columns and the
   !> 10 jet seeds
   !>
   !> - the caller fills the cavity-owned inputs of `state` first
   !>   (`fill_seed_state`); the level-set jet and the multiplier are set here
   !> - `kkt(:, 1:4)` is written here; further columns (the nuclear path's
   !>   anchor seeds) are set by the caller and solved on the same factorization
   !> - a degenerate point or a singular KKT system returns a `seed_state_*`
   !>   code other than `seed_state_ok`; the weights are then unspecified and the
   !>   caller must latch the point, never keep it
   !>
   !> @param[inout] phi            Projection objective of this thread
   !> @param[in]    point          Projected surface point
   !> @param[in]    owner_idx      Owner sphere of the anchor
   !> @param[in]    lsf1_r         Level-set gradient at the point
   !> @param[in]    lsf2_rr        Packed level-set Hessian at the point
   !> @param[in]    lsf3_rrr       Packed level-set third derivative at the point
   !> @param[in]    lambda_val     Lagrange multiplier of the projection
   !> @param[in]    f_crit         Critical-gradient switching function
   !> @param[in]    f_foc          Focusing switching function
   !> @param[in]    f_wleb         Lebedev-weight pruning switching function
   !> @param[in]    use_wleb_prune Whether Lebedev-weight pruning is active
   !> @param[in]    eff            Folded surface adjoints
   !> @param[in]    igrid          Grid point
   !> @param[inout] state          Seed state; cavity-owned inputs in, rest written
   !> @param[inout] kkt            `(4, >= 4)`; columns 5+ in, every column solved out
   !> @param[out]   phi1_r         Objective gradient at the point
   !> @param[out]   w_lsf0_pt      Level-set value adjoint
   !> @param[out]   w_lsf1_pt      Level-set gradient adjoint
   !> @param[out]   w_lsf2_pt      Packed symmetric level-set Hessian adjoint
   !> @param[out]   w_xyz_pt       Effective position adjoint for every seed
   !> @param[out]   status         `seed_state_ok` or the failure code
   subroutine seed_point_jet(phi, point, owner_idx, lsf1_r, lsf2_rr, lsf3_rrr, lambda_val, &
                             f_crit, f_foc, f_wleb, use_wleb_prune, eff, igrid, state, kkt, &
                             phi1_r, w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, w_xyz_pt, status)
      !> Projection objective of this thread
      type(moist_cavity_drop_phi_type), intent(inout) :: phi
      !> Projected surface point
      real(wp), intent(in) :: point(3)
      !> Owner sphere of the anchor
      integer, intent(in) :: owner_idx
      !> Level-set gradient
      real(wp), intent(in) :: lsf1_r(3)
      !> Packed level-set Hessian
      real(wp), intent(in) :: lsf2_rr(:)
      !> Packed level-set third derivative
      real(wp), intent(in) :: lsf3_rrr(:)
      !> Lagrange multiplier
      real(wp), intent(in) :: lambda_val
      !> Switching functions owned by the cavity
      class(moist_math_smoothstep_type), intent(in) :: f_crit, f_foc, f_wleb
      !> Whether Lebedev-weight pruning is active
      logical, intent(in) :: use_wleb_prune
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Seed state
      type(drop_seed_state_type), intent(inout) :: state
      !> Bordered KKT RHSs in, solutions out
      real(wp), intent(inout) :: kkt(:, :)
      !> Objective gradient
      real(wp), intent(out) :: phi1_r(3)
      !> Level-set value adjoint
      real(wp), intent(out) :: w_lsf0_pt
      !> Level-set gradient adjoint
      real(wp), intent(out) :: w_lsf1_pt(3)
      !> Packed level-set Hessian adjoint, off-diagonals not doubled
      real(wp), intent(out) :: w_lsf2_pt(6)
      !> Effective position adjoint
      real(wp), intent(out) :: w_xyz_pt(3)
      !> Kernel status
      integer, intent(out) :: status

      !> Objective value and packed Hessian
      real(wp) :: phi0, phi2_rr(6)
      !> LAPACK status of the KKT solve
      integer(lapack_ik) :: info

      w_lsf0_pt = 0.0_wp
      w_lsf1_pt = 0.0_wp
      w_lsf2_pt = 0.0_wp
      w_xyz_pt = 0.0_wp

      call phi%f012_r(point, state%anchor, owner_idx, phi0, phi1_r, phi2_rr)

      state%lsf1_r = lsf1_r
      state%lsf2_rr = lsf2_rr
      state%lsf3_rrr = lsf3_rrr
      state%lambda_val = lambda_val
      call build_seed_state(state, f_crit, f_foc, f_wleb, use_wleb_prune, status)
      if (status /= seed_state_ok) return

      ! Outward-normal channel:
      ! the direct grad-S term relies on w_lsf1 and is picked up by the field contraction;
      ! its point-motion coupling extends the effective position weight used by every seed
      call seed_normal_channel(state, eff, igrid, lsf2_rr, w_lsf1_pt, w_xyz_pt)

      ! Columns 1-4 are the level-set value and gradient seeds;
      ! the six Hessian seeds have a zero RHS and do not move the point
      kkt(:, 1:4) = 0.0_wp
      kkt(4, 1) = -1.0_wp
      kkt(1, 2) = lambda_val
      kkt(2, 3) = lambda_val
      kkt(3, 4) = lambda_val
      call drop_kkt_solve(phi2_rr - lambda_val*lsf2_rr, lsf1_r, kkt, info)
      if (info /= 0_lapack_ik) then
         status = seed_state_singular_kkt
         return
      end if

      call seed_jet_basis(state, eff, igrid, phi1_r, kkt, w_xyz_pt, w_lsf0_pt, w_lsf1_pt, &
                          w_lsf2_pt)
   end subroutine seed_point_jet

   !> Solve the bordered KKT sensitivity system for a batch of RHSs
   !>
   !> With the Lagrangian Hessian `H_L = phi_rr - lambda S_rr` and the level-set
   !> gradient `g = S_r`,
   !>
   !> ```
   !>   [ H_L  -g ] [ dr/dp      ]   [ b_1:3 ]
   !>   [ g^T   0 ] [ dlambda/dp ] = [ b_4   ]
   !> ```
   !>
   !> The full bordered system is solved rather than eliminating `dlambda`,
   !> so `H_L` itself need not be invertible
   !>
   !> `H_L` is packed, order of [[moist_math_packed]]
   !>
   !> @param[in]    H_lagrangian Packed Lagrangian Hessian at the projected point
   !> @param[in]    lsf1_r       Level-set gradient at the projected point
   !> @param[inout] rhs          `(4, nrhs)`; RHSs in, solutions out
   !> @param[out]   info         LAPACK status; nonzero means singular
   subroutine drop_kkt_solve(H_lagrangian, lsf1_r, rhs, info)
      !> Packed Lagrangian Hessian
      real(wp), intent(in) :: H_lagrangian(:)
      !> Level-set gradient
      real(wp), intent(in) :: lsf1_r(3)
      !> RHSs in, solutions out
      real(wp), intent(inout) :: rhs(:, :)
      !> LAPACK status
      integer(lapack_ik), intent(out) :: info

      !> Bordered matrix, destroyed by the factorization
      real(wp) :: kkt_mat(4, 4)
      !> Pivot indices
      integer(lapack_ik) :: ipiv(4)

      kkt_mat = 0.0_wp
      kkt_mat(1, 1) = H_lagrangian(1)
      kkt_mat(1, 2) = H_lagrangian(2)
      kkt_mat(2, 1) = H_lagrangian(2)
      kkt_mat(2, 2) = H_lagrangian(3)
      kkt_mat(1, 3) = H_lagrangian(4)
      kkt_mat(3, 1) = H_lagrangian(4)
      kkt_mat(2, 3) = H_lagrangian(5)
      kkt_mat(3, 2) = H_lagrangian(5)
      kkt_mat(3, 3) = H_lagrangian(6)
      kkt_mat(1:3, 4) = -lsf1_r
      kkt_mat(4, 1:3) = lsf1_r

      call lapack_gesv(4_lapack_ik, int(size(rhs, 2), lapack_ik), kkt_mat, &
                       4_lapack_ik, ipiv, rhs, 4_lapack_ik, info)
   end subroutine drop_kkt_solve

   !> Fold an outward-normal adjoint into the level-set gradient and position channels
   !>
   !> The normal `n = grad S / |grad S|` depends on the level set twice over:
   !>
   !> - at a fixed point through `grad S`, landing on the gradient channel
   !> - through the point's own motion, landing on the position channel as
   !>   `H @ normal_grad`
   !>
   !> Call this only once the point is known to be usable -- a rejected point
   !> must not leave a half-contracted weight behind
   !>
   !> @param[in]    state     Per-grid point forward state
   !> @param[in]    eff       Folded surface adjoints
   !> @param[in]    igrid     Grid point
   !> @param[in]    lsf2_rr   Packed level-set Hessian at the projected point
   !> @param[inout] w_lsf1_pt Point-local level-set gradient adjoint
   !> @param[out]   w_xyz_pt  Effective position adjoint for every seed
   pure subroutine seed_normal_channel(state, eff, igrid, lsf2_rr, w_lsf1_pt, w_xyz_pt)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Packed level-set Hessian
      real(wp), intent(in) :: lsf2_rr(:)
      !> Point-local level-set gradient adjoint
      real(wp), intent(inout) :: w_lsf1_pt(3)
      !> Effective position adjoint
      real(wp), intent(out) :: w_xyz_pt(3)

      !> Tangential part of the normal adjoint, divided by |grad S|
      real(wp) :: normal_grad(3)
      !> Normal component of the normal adjoint
      real(wp) :: nwn
      !> `H normal_grad`
      real(wp) :: H_normal_grad(3)

      if (eff%have_wn) then
         nwn = dot_product(state%n_surf, eff%w_n(:, igrid))
         normal_grad = (eff%w_n(:, igrid) - state%n_surf*nwn)/state%g_norm
         w_lsf1_pt = w_lsf1_pt + normal_grad
         call packed_contract2(lsf2_rr, normal_grad, H_normal_grad, 3)
         w_xyz_pt = eff%w_xyz(:, igrid) + H_normal_grad
      else
         w_xyz_pt = eff%w_xyz(:, igrid)
      end if
   end subroutine seed_normal_channel

   !> Push the 10 level-set jet directions through the kernel
   !>
   !> The per-point map is linear in its seed, so seeding each basis direction of
   !> the jet `(S, grad S, grad^2 S)` once and collecting the responses gives the
   !> adjoint of the whole jet
   !>
   !> - Hessian directions in packed order of [[moist_math_packed]]: `E_aa`
   !>   on the diagonal, `(E_ab + E_ba)/2` off it
   !> - the response to `(E_ab + E_ba)/2` is `(W_ab + W_ba)/2`, so `w_lsf2_pt`
   !>   is the packed symmetric adjoint, off-diagonals not doubled
   !> - only the value and gradient directions move the point, the six Hessian
   !>   directions having `dr/dp = 0`
   !> - the results are point-local: `potential.f90` scatters them into
   !>   `w_lsf*(..., igrid)`, `nuclear.f90` contracts them with `lsf*_rA`
   !>
   !> @param[in]    state     Per-grid point forward state
   !> @param[in]    eff       Folded surface adjoints
   !> @param[in]    igrid     Grid point
   !> @param[in]    phi1_r    Objective gradient at the projected point
   !> @param[in]    kkt       Solved KKT sensitivities; columns 1-4 are the jet seeds
   !> @param[in]    w_xyz_pt  Effective position adjoint from [[seed_normal_channel]]
   !> @param[inout] w_lsf0_pt Point-local level-set value adjoint
   !> @param[inout] w_lsf1_pt Point-local level-set gradient adjoint
   !> @param[inout] w_lsf2_pt Point-local packed level-set Hessian adjoint
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
      !> Point-local level-set adjoints, Hessian packed
      real(wp), intent(inout) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(6)

      !> Linear response to one seed
      type(drop_seed_result_type) :: res
      !> Seed perturbation of the level-set jet
      real(wp) :: dlsf1_r(3), dlsf2_rr(3, 3)
      !> Induced point motion and multiplier change
      real(wp) :: dr_dp(3), dlambda_dp
      !> Adjoint contribution of one seed
      real(wp) :: contribution
      !> Seed and packed Hessian indices
      integer :: ibasis, ipack
      !> Cartesian indices of the packed slot
      integer :: iaxis, jaxis
      !> Cartesian pair of each packed slot: xx, xy, yy, xz, yz, zz
      integer, parameter :: packed_pair(2, 6) = &
         & reshape([1, 1, 1, 2, 2, 2, 1, 3, 2, 3, 3, 3], [2, 6])

      do ibasis = 1, drop_n_jet_seeds
         dlsf1_r = 0.0_wp
         dlsf2_rr = 0.0_wp
         ipack = 0
         if (ibasis <= 4) then
            if (ibasis > 1) dlsf1_r(ibasis - 1) = 1.0_wp
            dr_dp = kkt(1:3, ibasis)
            dlambda_dp = kkt(4, ibasis)
         else
            ipack = ibasis - 4
            iaxis = packed_pair(1, ipack)
            jaxis = packed_pair(2, ipack)
            if (iaxis == jaxis) then
               dlsf2_rr(iaxis, iaxis) = 1.0_wp
            else
               dlsf2_rr(iaxis, jaxis) = 0.5_wp
               dlsf2_rr(jaxis, iaxis) = 0.5_wp
            end if
            dr_dp = 0.0_wp
            dlambda_dp = 0.0_wp
         end if

         call apply_seed(state, dlsf1_r, dlsf2_rr, dr_dp, dlambda_dp, res)

         ! TODO: This will change with the planned branching change
         ! No w_f term: the switching factor is an anchor-only iSwig overlap,
         ! so a level-set perturbation at fixed nuclei leaves it alone
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
            w_lsf2_pt(ipack) = w_lsf2_pt(ipack) + contribution
         end if
      end do
   end subroutine seed_jet_basis

   !> Push the three anchor directions through the kernel
   !>
   !> The anchor moves rigidly with its owner sphere, so `da_i/dR_I = delta`
   !> The level-set field is untouched; the whole channel enters through the
   !> objective, whose mixed derivative `-d^2 phi / dr dR_I` is `+alpha * I`
   !> That makes it three extra RHSs on the same factorization,
   !> independent of the number of spheres
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
      call fatal_error(error, context//": "//seed_status_message(status)// &
                       " at grid point "//trim(idx))
   end subroutine degenerate_point_error

end module moist_cavity_drop_derivatives_seeds
