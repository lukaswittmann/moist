!> Provides the seeds for the DROP derivative paths
!>
!> * A seed is one unit perturbation direction pushed through the linear
!>   per-point map of [[moist_cavity_drop_derivatives_kernel]]
!> * The kernel yields the surface observables for a given seed
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
   use moist_math_lapack_gesv, only: lapack_gesv
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, drop_seed_result_type, &
      & drop_surface_weights_type, apply_seed, seed_weight_tol, seed_status_message

   implicit none (type, external)
   private

   public :: drop_kkt_solve, seed_normal_channel, seed_jet_basis, seed_anchor
   public :: degenerate_point_error

   !> Number of level-set jet directions: one value, three gradient, nine Hessian
   integer, parameter, public :: drop_n_jet_seeds = 13

contains

   !> Solve the bordered KKT sensitivity system for a batch of right-hand sides
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
   !> so `H_L` itself need not be invertible.
   !>
   !> @param[in]    H_lagrangian Lagrangian Hessian at the projected point
   !> @param[in]    lsf1_r       Level-set gradient at the projected point
   !> @param[inout] rhs          `(4, nrhs)`; right-hand sides in, solutions out
   !> @param[out]   info         LAPACK status; nonzero means singular
   subroutine drop_kkt_solve(H_lagrangian, lsf1_r, rhs, info)
      !> Lagrangian Hessian
      real(wp), intent(in) :: H_lagrangian(3, 3)
      !> Level-set gradient
      real(wp), intent(in) :: lsf1_r(3)
      !> Right-hand sides in, solutions out
      real(wp), intent(inout) :: rhs(:, :)
      !> LAPACK status
      integer(lapack_ik), intent(out) :: info

      !> Bordered matrix, destroyed by the factorization
      real(wp) :: kkt_mat(4, 4)
      !> Pivot indices
      integer(lapack_ik) :: ipiv(4)

      kkt_mat = 0.0_wp
      kkt_mat(1:3, 1:3) = H_lagrangian
      kkt_mat(1:3, 4) = -lsf1_r
      kkt_mat(4, 1:3) = lsf1_r

      call lapack_gesv(4_lapack_ik, int(size(rhs, 2), lapack_ik), kkt_mat, &
                       4_lapack_ik, ipiv, rhs, 4_lapack_ik, info)
   end subroutine drop_kkt_solve

   !> Fold an outward-normal adjoint into the level-set gradient and position channels
   !>
   !> The normal `n = grad S / |grad S|` depends on the level set twice over: at
   !> a fixed point through `grad S`, which lands on the gradient channel, and
   !> through the point's own motion, which lands on the position channel as
   !> `H @ normal_grad`.
   !>
   !> Call this only once the point is known to be usable -- a rejected point
   !> must not leave a half-contracted weight behind.
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
