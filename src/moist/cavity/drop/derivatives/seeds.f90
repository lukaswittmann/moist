!> Seed layout and seed pushes of the DROP derivatives
!>
!> - Seed: unit perturbation direction pushed through the linear per-point
!>   map of [[moist_cavity_drop_derivatives_kernel]]
!> - Two bases: 13 level-set jet directions spanning `(S, grad S, grad^2 S)`,
!>   3 anchor directions for the rigid motion of the owner sphere
!> - Projected point is an implicit minimizer; the bordered KKT solve gives
!>   the point motion of each seed
!> - Outward-normal adjoint folded into the position and gradient channels
!> - Each basis as one sweep and as `_apply`/`_contract` pair: apply reads
!>   the grid point only, contract takes the surface adjoints; the sweep
!>   composes both, one copy of every floating-point chain
!> - Seed layout owned here: [[jet_seed_index]], [[seed_rhs_column]],
!>   [[seed_standard_rhs]]
!> - No `self`, plain arguments: callable from a submodule of
!>   `moist_cavity_drop`, testable without a cavity
module moist_cavity_drop_derivatives_seeds
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use moist_math_lapack_getrf, only: lapack_getrf
   use moist_math_lapack_getrs, only: lapack_getrs
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, drop_seed_result_type, &
      & drop_seed_state_tangent_type, &
      & drop_surface_weights_type, apply_seed, &
      & seed_contribution, seed_weight_tol

   implicit none(type, external)
   private

   public :: drop_kkt_factor_type, seed_normal_channel, seed_jet_basis, seed_anchor
   public :: seed_jet_basis_apply, seed_jet_basis_contract
   public :: seed_anchor_apply, seed_anchor_contract
   public :: jet_seed_index, seed_rhs_column, seed_standard_rhs
   public :: fill_seed_basis, scatter_jet_weight, seed_contribution_tangent

   !> Number of level-set jet directions: 1 value, 3 gradient, 9 Hessian
   integer, parameter, public :: drop_n_jet_seeds = 13

   !> Number of anchor directions, rigid motion of the owner sphere
   integer, parameter, public :: drop_n_anchor_seeds = 3

   !> Seeds per grid point: jet directions first, then anchor directions
   !>
   !> - Order laid out by [[fill_seed_basis]], walked by `hessian_traverse.f90`
   integer, parameter, public :: drop_n_point_seeds = drop_n_jet_seeds + drop_n_anchor_seeds

   !> Jet directions moving the projected point: value and three gradient
   !>
   !> - Hessian directions have a zero right-hand side and no batch column
   integer, parameter, public :: drop_n_moving_jet_seeds = 1 + 3

   !> Columns of the [[seed_standard_rhs]] batch: moving jet, then anchor
   integer, parameter, public :: drop_n_seed_columns = drop_n_moving_jet_seeds &
                                                       + drop_n_anchor_seeds

   !> Zero gradient tangent of a basis seed, for [[apply_seed_tangent]]
   !>
   !> - Basis seeds are constant; only their induced point motion has a tangent
   real(wp), parameter, public :: seed_dzero1(3) = 0.0_wp
   !> Zero Hessian tangent of a basis seed, for [[apply_seed_tangent]]
   real(wp), parameter, public :: seed_dzero2(3, 3) = 0.0_wp

   !> Slot kind beyond `drop_n_jet_seeds`, e.g. an anchor seed
   !>
   !> - Slot kinds assigned by [[jet_seed_index]], sole owner of the jet layout
   integer, parameter, public :: drop_jet_seed_none = 0
   !> Slot 1: level-set value direction
   integer, parameter, public :: drop_jet_seed_value = 1
   !> Slots 2-4: level-set gradient directions
   integer, parameter, public :: drop_jet_seed_grad = 2
   !> Slots 5-13: level-set Hessian directions, row-major
   integer, parameter, public :: drop_jet_seed_hess = 3

   !> LU factorization of the bordered KKT sensitivity matrix
   !>
   !> ```
   !>   [ H_L  -g ] [ dr/dp      ]   [ b_1:3 ]
   !>   [ g^T   0 ] [ dlambda/dp ] = [ b_4   ]
   !> ```
   !>
   !> - `H_L = phi_rr - lambda S_rr`, `g = S_r`
   !> - Full bordered system factored, `H_L` need not be invertible
   !> - Factored once per grid point: `solve` for the primal right-hand sides,
   !>   `solve_tangent` for `K dx = db - dK x` once the primal solution exists
   !> - Fixed size: stack-local in the OpenMP grid loops, allocation only for
   !>   the error object
   !> - Tangent batch buffer owned by the caller
   !> - Diagnostics prefixed with `context//": "` plus optional `igrid`, as in
   !>   [[check_surface_adjoint]]
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
   !> - Singular matrix: degenerate projected point, reported through `error`
   !>
   !> @param[out] self          factorization
   !> @param[in]  H_lagrangian  Lagrangian Hessian at the projected point
   !> @param[in]  lsf1_r        level-set gradient at the projected point
   !> @param[in]  context       API entry point, prefixes the diagnostic
   !> @param[out] error         allocated when the system is singular
   !> @param[in]  igrid         grid point for the diagnostic, optional
   subroutine drop_kkt_factor(self, H_lagrangian, lsf1_r, context, error, igrid)
      !> Factorization
      class(drop_kkt_factor_type), intent(out) :: self
      !> Lagrangian Hessian
      real(wp), intent(in) :: H_lagrangian(3, 3)
      !> Level-set gradient
      real(wp), intent(in) :: lsf1_r(3)
      !> Calling API entry point
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Failing grid point
      integer, intent(in), optional :: igrid

      !> LAPACK status
      integer(lapack_ik) :: info

      self%lu = 0.0_wp
      self%lu(1:3, 1:3) = H_lagrangian
      self%lu(1:3, 4) = -lsf1_r
      self%lu(4, 1:3) = lsf1_r

      call lapack_getrf(4_lapack_ik, 4_lapack_ik, self%lu, 4_lapack_ik, self%ipiv, info)
      if (info /= 0_lapack_ik) then
         call kkt_report(context, "Bordered KKT sensitivity matrix is singular", "getrf", &
                         info, error, igrid)
      end if
   end subroutine drop_kkt_factor

   !> Turn a LAPACK status of the bordered system into an error object
   !>
   !> - Fixed-length renderings, thread safe inside parallel regions
   !>
   !> @param[in]  context  API entry point, prefixes the diagnostic
   !> @param[in]  what     failure description
   !> @param[in]  routine  LAPACK routine reporting the status
   !> @param[in]  info     LAPACK status
   !> @param[out] error    error object, always allocated
   !> @param[in]  igrid    grid point for the diagnostic, optional
   subroutine kkt_report(context, what, routine, info, error, igrid)
      !> Calling API entry point
      character(len=*), intent(in) :: context
      !> Failure description and reporting LAPACK routine
      character(len=*), intent(in) :: what, routine
      !> LAPACK status
      integer(lapack_ik), intent(in) :: info
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Failing grid point
      integer, intent(in), optional :: igrid

      !> Rendered status, fixed length for thread safety
      character(len=32) :: status
      !> Rendered grid index and its message suffix, fixed length
      character(len=32) :: idx
      character(len=48) :: at_point

      write (status, "(i0)") info
      at_point = ""
      if (present(igrid)) then
         write (idx, "(i0)") igrid
         at_point = " at grid point "//trim(idx)
      end if
      call fatal_error(error, context//": "//what//" ("//routine//" status "// &
                       trim(status)//")"//trim(at_point))
   end subroutine kkt_report

   !> Solve a right-hand side batch with the stored factors
   !>
   !> - Requires a successful [[drop_kkt_factor]]
   !> - `getrs` status flags an illegal argument only: unreachable from a
   !>   correct caller, not covered by a test
   !>
   !> @param[in]     self     factorization
   !> @param[in,out] rhs      `(4, nrhs)`; right-hand sides in, solutions out
   !> @param[in]     context  API entry point, prefixes the diagnostic
   !> @param[out]    error    allocated when LAPACK rejects the call
   !> @param[in]     igrid    grid point for the diagnostic, optional
   subroutine drop_kkt_apply(self, rhs, context, error, igrid)
      !> Factorization
      class(drop_kkt_factor_type), intent(in) :: self
      !> Right-hand sides in, solutions out
      real(wp), contiguous, intent(inout) :: rhs(:, :)
      !> Calling API entry point
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Failing grid point
      integer, intent(in), optional :: igrid

      !> LAPACK status
      integer(lapack_ik) :: info

      call lapack_getrs("n", 4_lapack_ik, int(size(rhs, 2), lapack_ik), self%lu, &
                        4_lapack_ik, self%ipiv, rhs, 4_lapack_ik, info)
      if (info /= 0_lapack_ik) then
         call kkt_report(context, "Bordered KKT sensitivity solve failed", "getrs", &
                         info, error, igrid)
      end if
   end subroutine drop_kkt_apply

   !> Solve the tangent of a right-hand side batch with the same factors
   !>
   !> ```
   !>   K dx = db - dK x,   dK = [ dH_L  -dg ]
   !>                            [ dg^T    0 ]
   !> ```
   !>
   !> - `dlsf1_r` terms enter with opposite signs: rows 1-3 add the
   !>   multiplier term, row 4 subtracts the position term
   !> - `phi_rr = alpha*I` is constant: caller passes
   !>   `dH_lagrangian = -dlambda S_rr - lambda d(S_rr)`
   !> - One `getrs` for all `nseed*ndir` right-hand sides of a grid point
   !> - Column `(idir-1)*nseed + iseed`, seeds of one direction contiguous;
   !>   the driver depends on that
   !> - `rhs` allocated by the caller, once per thread outside the grid loop
   !> - Requires a successful [[drop_kkt_factor]]
   !> - Shape mismatch: reported before `rhs` is touched, reachable and tested
   !> - `getrs` status unreachable from a correct caller, see [[drop_kkt_apply]]
   !>
   !> @param[in]     self           factorization from a prior `factor` call
   !> @param[in]     dH_lagrangian  tangent of the Lagrangian Hessian, `(3, 3, ndir)`
   !> @param[in]     dlsf1_r        tangent of the level-set gradient, `(3, ndir)`
   !> @param[in]     x              primal solutions from `solve`, `(4, nseed)`
   !> @param[in,out] rhs            `db` in, `dx` out; `(4, nseed*ndir)`
   !> @param[in]     context        API entry point, prefixes the diagnostic
   !> @param[out]    error          allocated on inconsistent shapes or LAPACK status
   !> @param[in]     igrid          grid point for the diagnostic, optional
   subroutine drop_kkt_solve_tangent(self, dH_lagrangian, dlsf1_r, x, rhs, context, &
                                     error, igrid)
      !> Factorization from a prior `factor` call
      class(drop_kkt_factor_type), intent(in) :: self
      !> Tangent of the Lagrangian Hessian block, `(3, 3, ndir)`
      real(wp), contiguous, intent(in) :: dH_lagrangian(:, :, :)
      !> Tangent of the level-set gradient, `(3, ndir)`
      real(wp), contiguous, intent(in) :: dlsf1_r(:, :)
      !> Primal solutions from `solve`, `(4, nseed)`
      real(wp), contiguous, intent(in) :: x(:, :)
      !> `db` in, `dx` out; `(4, nseed*ndir)`, column `(idir-1)*nseed + iseed`
      real(wp), contiguous, intent(inout) :: rhs(:, :)
      !> Calling API entry point
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Failing grid point
      integer, intent(in), optional :: igrid

      !> Batch extents, from the tangent and primal arguments
      integer :: ndir, nseed
      !> Direction, seed, column and Cartesian indices
      integer :: idir, iseed, icol, iaxis
      !> LAPACK status
      integer(lapack_ik) :: info
      !> Rendered shapes, fixed length for thread safety
      character(len=64) :: shapes
      !> Rendered grid index and its message suffix, fixed length
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
               ! Rows 1-3: -dg block in column 4, multiplier term adds
               rhs(iaxis, icol) = rhs(iaxis, icol) &
                                  - dH_lagrangian(iaxis, 1, idir)*x(1, iseed) &
                                  - dH_lagrangian(iaxis, 2, idir)*x(2, iseed) &
                                  - dH_lagrangian(iaxis, 3, idir)*x(3, iseed) &
                                  + dlsf1_r(iaxis, idir)*x(4, iseed)
            end do
            ! Row 4: +dg^T, position term subtracts
            rhs(4, icol) = rhs(4, icol) &
                           - dlsf1_r(1, idir)*x(1, iseed) &
                           - dlsf1_r(2, idir)*x(2, iseed) &
                           - dlsf1_r(3, idir)*x(3, iseed)
         end do
      end do

      call lapack_getrs("n", 4_lapack_ik, int(size(rhs, 2), lapack_ik), self%lu, &
                        4_lapack_ik, self%ipiv, rhs, 4_lapack_ik, info)
      if (info /= 0_lapack_ik) then
         call kkt_report(context, "Tangent KKT sensitivity solve failed", "getrs", &
                         info, error, igrid)
      end if
   end subroutine drop_kkt_solve_tangent

   !> Fold an outward-normal adjoint into the level-set gradient and position channels
   !>
   !> - `n = grad S / |grad S|`: gradient channel gets `normal_grad`, position
   !>   channel gets `lsf2_rr @ normal_grad` from the point motion
   !> - Sole author of `normal_grad` and `nwn`: second-order callers request
   !>   them here and never rebuild them, one floating-point chain
   !> - Both zero when the outward-normal channel is inactive
   !>
   !> @param[in]     state           per-grid point forward state
   !> @param[in]     eff             folded surface adjoints
   !> @param[in]     igrid           grid point
   !> @param[in]     lsf2_rr         level-set Hessian at the projected point
   !> @param[in,out] w_lsf1_pt       point-local level-set gradient adjoint
   !> @param[out]    w_xyz_pt        effective position adjoint for every seed
   !> @param[out]    normal_grad_pt  tangential normal adjoint over |grad S|, optional
   !> @param[out]    nwn_pt          normal component of the normal adjoint, optional
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
   !> - Slot 1 value, slots 2-4 gradient, slots 5-13 Hessian in row-major order
   !> - Sole owner of the layout; [[seed_jet_basis]], [[fill_seed_basis]] and
   !>   [[scatter_jet_weight]] consume it
   !> - `iaxis`, `jaxis` always defined: both zero for the value slot and slots
   !>   beyond the layout, `jaxis` zero for a gradient slot
   !>
   !> @param[in]  ibasis     seed slot, positive
   !> @param[out] slot_kind  one of the `drop_jet_seed_*` kinds
   !> @param[out] iaxis      first Cartesian index of the slot, or zero
   !> @param[out] jaxis      second Cartesian index of the slot, or zero
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

   !> Column of the standard seed batch holding one point seed's solution
   !>
   !> - Sole map between seed slots and [[seed_standard_rhs]] columns, built
   !>   on [[jet_seed_index]]
   !> - Zero for the nine Hessian slots and slots beyond the layout: no
   !>   column, seed moves neither the point nor the multiplier
   !>
   !> @param[in] ibasis  point-seed slot, `1 .. drop_n_point_seeds`
   pure function seed_rhs_column(ibasis) result(icol)
      !> Point-seed slot
      integer, intent(in) :: ibasis
      !> Batch column of the slot, zero when it has none
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
   !> - Solved once per grid point on the [[drop_kkt_factor]] factors
   !> - Columns follow [[seed_rhs_column]]
   !> - Value seed: unit change of `S` drives `S = 0`, row 4 is -1
   !> - Gradient seed: multiplier term of the Lagrangian, row `iaxis` is `lambda`
   !> - Anchor seed: field untouched, enters through
   !>   `-d^2 phi/(dr dR_owner) = +alpha I`, row `iaxis` is `phi_alpha`
   !>
   !> @param[in]  lambda_val  Lagrange multiplier of the projection
   !> @param[in]  phi_alpha   quadratic coefficient of the projection objective
   !> @param[out] kkt_rhs     right-hand sides, `(4, drop_n_seed_columns)`
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
   !> - Slots 1-13: jet directions of [[seed_jet_basis]]; slots 14-16: anchor
   !>   directions of [[seed_anchor]]
   !> - Hessian seeds are single-entry matrices, hence asymmetric; their `x`
   !>   is zero
   !> - `x` in the `(4, nseed)` layout [[drop_kkt_solve_tangent]] expects
   !> - Slot layout from [[jet_seed_index]], batch columns from [[seed_rhs_column]]
   !>
   !> @param[in]  kkt       solved KKT sensitivities, [[seed_standard_rhs]] layout
   !> @param[out] dlsf1_r   gradient perturbation of each seed
   !> @param[out] dlsf2_rr  Hessian perturbation of each seed
   !> @param[out] x         induced point motion and multiplier change of each seed
   pure subroutine fill_seed_basis(kkt, dlsf1_r, dlsf2_rr, x)
      !> Solved KKT sensitivities
      real(wp), intent(in) :: kkt(:, :)
      !> Gradient perturbation of each seed
      real(wp), intent(out) :: dlsf1_r(3, drop_n_point_seeds)
      !> Hessian perturbation of each seed
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

         ! Anchor seed: no jet perturbation, described by its point motion alone
         if (ibasis > drop_n_jet_seeds) cycle

         call jet_seed_index(ibasis, slot_kind, iaxis, jaxis)
         if (slot_kind == drop_jet_seed_hess) then
            dlsf2_rr(iaxis, jaxis, ibasis) = 1.0_wp
         else if (slot_kind == drop_jet_seed_grad) then
            dlsf1_r(iaxis, ibasis) = 1.0_wp
         end if
      end do
   end subroutine fill_seed_basis

   !> Add one jet seed's contribution to the level-set adjoint weights
   !>
   !> - Slot classified by [[jet_seed_index]], same layout as [[fill_seed_basis]]
   !> - Anchor slots classify as `drop_jet_seed_none` and land nowhere; their
   !>   contribution belongs to the owner's gradient row
   !>
   !> @param[in]     ibasis   seed slot, `1 .. drop_n_jet_seeds`
   !> @param[in]     contrib  contribution of that seed
   !> @param[in,out] w0       level-set value adjoint
   !> @param[in,out] w1       level-set gradient adjoint
   !> @param[in,out] w2       level-set Hessian adjoint
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
   !> - Per-point map is linear in its seed: the responses of all basis
   !>   directions give the adjoint of the whole jet
   !> - Only value and gradient directions move the point; Hessian directions
   !>   have `dr/dp = 0`
   !> - Results are point-local: `potential.f90` scatters them into
   !>   `w_lsf*(..., igrid)`, `nuclear.f90` contracts them with `lsf*_rA`
   !> - Composes [[seed_jet_basis_apply]] and [[seed_jet_basis_contract]]; with
   !>   several weight sets per point call the halves directly, as
   !>   [[drop_hessian_traverse]] does
   !>
   !> @param[in]     state      per-grid point forward state
   !> @param[in]     eff        folded surface adjoints
   !> @param[in]     igrid      grid point
   !> @param[in]     phi1_r     objective gradient at the projected point
   !> @param[in]     kkt        solved KKT sensitivities, read through [[seed_rhs_column]]
   !> @param[in]     w_xyz_pt   effective position adjoint from [[seed_normal_channel]]
   !> @param[in,out] w_lsf0_pt  point-local level-set value adjoint
   !> @param[in,out] w_lsf1_pt  point-local level-set gradient adjoint
   !> @param[in,out] w_lsf2_pt  point-local level-set Hessian adjoint
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

      !> Linear response and induced point motion of each seed
      type(drop_seed_result_type) :: res_seed(drop_n_jet_seeds)
      real(wp) :: seed_x(4, drop_n_jet_seeds)

      call seed_jet_basis_apply(state, kkt, res_seed, seed_x)
      call seed_jet_basis_contract(eff, igrid, phi1_r, w_xyz_pt, res_seed, seed_x, &
                                   w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)
   end subroutine seed_jet_basis

   !> Apply the 13 jet seeds at one grid point, first half of [[seed_jet_basis]]
   !>
   !> - Reads `state` and the bordered solve only, no surface adjoint: run
   !>   once per point, contract once per weight set
   !> - `seed_x` in the `(4, nseed)` layout of the second-order chain, zero
   !>   for Hessian seeds
   !> - Slots classified by [[jet_seed_index]], as in [[seed_jet_basis_contract]]
   !> - `dstate` chain of [[apply_seed]] computed only when requested
   !>
   !> @param[in]  state     per-grid point forward state
   !> @param[in]  kkt       solved KKT sensitivities, read through [[seed_rhs_column]]
   !> @param[out] res_seed  linear response of each seed, `(drop_n_jet_seeds)`
   !> @param[out] seed_x    induced point motion and multiplier change, `(4, nseed)`
   !> @param[out] dstate    tangent of the derived seed state, optional
   subroutine seed_jet_basis_apply(state, kkt, res_seed, seed_x, dstate)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Solved KKT sensitivities
      real(wp), intent(in) :: kkt(:, :)
      !> Linear response of each seed
      type(drop_seed_result_type), intent(out) :: res_seed(drop_n_jet_seeds)
      !> Induced point motion and multiplier change of each seed
      real(wp), intent(out) :: seed_x(4, drop_n_jet_seeds)
      !> Tangent of the derived block of `state` along each seed
      type(drop_seed_state_tangent_type), intent(out), optional :: dstate(drop_n_jet_seeds)

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

         ! Element of an absent optional array cannot be forwarded
         if (present(dstate)) then
            call apply_seed(state, dlsf1_r, dlsf2_rr, seed_x(1:3, ibasis), &
                            seed_x(4, ibasis), res_seed(ibasis), dstate(ibasis))
         else
            call apply_seed(state, dlsf1_r, dlsf2_rr, seed_x(1:3, ibasis), &
                            seed_x(4, ibasis), res_seed(ibasis))
         end if
      end do
   end subroutine seed_jet_basis_apply

   !> Contract the 13 jet responses with one weight set
   !>
   !> - Second half of [[seed_jet_basis]], the only one taking surface adjoints
   !> - Each seed lands in its own adjoint slot, no cross-seed accumulation:
   !>   repeatable per weight set after one [[seed_jet_basis_apply]]
   !>
   !> @param[in]     eff        folded surface adjoints
   !> @param[in]     igrid      grid point
   !> @param[in]     phi1_r     objective gradient at the projected point
   !> @param[in]     w_xyz_pt   effective position adjoint from [[seed_normal_channel]]
   !> @param[in]     res_seed   linear response of each seed
   !> @param[in]     seed_x     induced point motion and multiplier change of each seed
   !> @param[in,out] w_lsf0_pt  point-local level-set value adjoint
   !> @param[in,out] w_lsf1_pt  point-local level-set gradient adjoint
   !> @param[in,out] w_lsf2_pt  point-local level-set Hessian adjoint
   pure subroutine seed_jet_basis_contract(eff, igrid, phi1_r, w_xyz_pt, res_seed, seed_x, &
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
      !> Seed index
      integer :: ibasis

      do ibasis = 1, drop_n_jet_seeds
         ! Field seed leaves the anchor alone: no rigid-motion shift, no
         ! switching term, see [[seed_contribution]]
         contribution = seed_contribution(eff, igrid, w_xyz_pt, seed_x(1:3, ibasis), &
                                          res_seed(ibasis), phi1_r)
         call scatter_jet_weight(ibasis, contribution, w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)
      end do
   end subroutine seed_jet_basis_contract

   !> Push the three anchor directions through the kernel
   !>
   !> - Anchor moves rigidly with its owner sphere, `da_i/dR_I = delta`
   !> - Level-set field untouched; the channel enters through the objective,
   !>   `-d^2 phi / dr dR_I = +alpha * I`
   !> - Three extra right-hand sides on the same factorization, independent
   !>   of the number of spheres
   !> - Composes [[seed_anchor_apply]] and [[seed_anchor_contract]], split as
   !>   in [[seed_jet_basis]]
   !>
   !> @param[in]     state       per-grid point forward state
   !> @param[in]     eff         folded surface adjoints
   !> @param[in]     igrid       grid point
   !> @param[in]     phi1_r      objective gradient at the projected point
   !> @param[in]     kkt         solved KKT sensitivities, read through [[seed_rhs_column]]
   !> @param[in]     w_xyz_pt    effective position adjoint from [[seed_normal_channel]]
   !> @param[in,out] grad_owner  nuclear-gradient accumulator of the owner sphere
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

      !> Linear response and induced point motion of each seed
      type(drop_seed_result_type) :: res_seed(drop_n_anchor_seeds)
      real(wp) :: seed_x(4, drop_n_anchor_seeds)

      call seed_anchor_apply(state, kkt, res_seed, seed_x)
      call seed_anchor_contract(eff, igrid, phi1_r, w_xyz_pt, res_seed, seed_x, &
                                grad_owner)
   end subroutine seed_anchor

   !> Apply the three anchor seeds at one grid point, first half of [[seed_anchor]]
   !>
   !> - `dstate` optional as in [[seed_jet_basis_apply]]: one caller collects
   !>   the whole `drop_n_point_seeds` basis in one pass
   !>
   !> @param[in]  state     per-grid point forward state
   !> @param[in]  kkt       solved KKT sensitivities, read through [[seed_rhs_column]]
   !> @param[out] res_seed  linear response of each seed, `(drop_n_anchor_seeds)`
   !> @param[out] seed_x    induced point motion and multiplier change, `(4, nseed)`
   !> @param[out] dstate    tangent of the derived seed state, optional
   subroutine seed_anchor_apply(state, kkt, res_seed, seed_x, dstate)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Solved KKT sensitivities
      real(wp), intent(in) :: kkt(:, :)
      !> Linear response of each seed
      type(drop_seed_result_type), intent(out) :: res_seed(drop_n_anchor_seeds)
      !> Induced point motion and multiplier change of each seed
      real(wp), intent(out) :: seed_x(4, drop_n_anchor_seeds)
      !> Tangent of the derived block of `state` along each seed
      type(drop_seed_state_tangent_type), intent(out), optional :: dstate(drop_n_anchor_seeds)

      !> Zero field perturbation: rigid anchor motion leaves the level set alone
      real(wp) :: dlsf1_r(3), dlsf2_rr(3, 3)
      !> Cartesian index
      integer :: iaxis

      dlsf1_r = 0.0_wp
      dlsf2_rr = 0.0_wp

      do iaxis = 1, drop_n_anchor_seeds
         seed_x(:, iaxis) = kkt(1:4, seed_rhs_column(drop_n_jet_seeds + iaxis))

         if (present(dstate)) then
            call apply_seed(state, dlsf1_r, dlsf2_rr, seed_x(1:3, iaxis), &
                            seed_x(4, iaxis), res_seed(iaxis), dstate(iaxis))
         else
            call apply_seed(state, dlsf1_r, dlsf2_rr, seed_x(1:3, iaxis), &
                            seed_x(4, iaxis), res_seed(iaxis))
         end if
      end do
   end subroutine seed_anchor_apply

   !> Contract the three anchor responses with one weight set
   !>
   !> - Second half of [[seed_anchor]]
   !>
   !> @param[in]     eff         folded surface adjoints
   !> @param[in]     igrid       grid point
   !> @param[in]     phi1_r      objective gradient at the projected point
   !> @param[in]     w_xyz_pt    effective position adjoint from [[seed_normal_channel]]
   !> @param[in]     res_seed    linear response of each seed
   !> @param[in]     seed_x      induced point motion and multiplier change of each seed
   !> @param[in,out] grad_owner  nuclear-gradient accumulator of the owner sphere
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
         ! phi = 0.5*alpha*|r - anchor|^2: rigid owner motion at fixed r adds
         ! -phi1_r, the shift [[seed_contribution]] takes
         contribution = seed_contribution(eff, igrid, w_xyz_pt, seed_x(1:3, iaxis), &
                                          res_seed(iaxis), phi1_r, phi1_r(iaxis))

         grad_owner(iaxis) = grad_owner(iaxis) + contribution
      end do
   end subroutine seed_anchor_contract

   !> Directional derivative of [[seed_contribution]] at fixed surface adjoints
   !>
   !> - Product rule on the primal contraction; sole caller is the
   !>   fixed-adjoint half of the Hessian
   !> - Position adjoint still moves through its normal fold, hence `dw_xyz_pt`
   !> - No switching term, as in [[seed_contribution]]
   !> - Branch term: tangent of `branch_phi_adj * (phi1_r . dr - branch_shift)`
   !>   with `dphi1_r = alpha (dr_v - v_owner)`
   !> - Anchor seed: `dbranch_shift` is the matching component of `dphi1_r`
   !> - Branch term gated on the primal's condition and threshold
   !> - Accumulation order as in the primal: position, width, branch, curvature
   !>
   !> @param[in] eff            folded surface adjoints
   !> @param[in] igrid          grid point
   !> @param[in] w_xyz_pt       effective position adjoint
   !> @param[in] dw_xyz_pt      tangent of the effective position adjoint
   !> @param[in] dr             induced point motion of the seed
   !> @param[in] ddr            tangent of that point motion
   !> @param[in] dres           second-order response of the seed
   !> @param[in] phi1_r         objective gradient at the projected point
   !> @param[in] dphi1_r        tangent of the objective gradient
   !> @param[in] dbranch_shift  tangent of the rigid-motion shift, anchor seeds only
   pure function seed_contribution_tangent(eff, igrid, w_xyz_pt, dw_xyz_pt, dr, ddr, dres, &
                                           phi1_r, dphi1_r, dbranch_shift) &
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
      type(drop_seed_result_type), intent(in) :: dres
      !> Objective gradient at the projected point and its tangent
      real(wp), intent(in) :: phi1_r(3), dphi1_r(3)
      !> Tangent of the rigid-motion shift, present for an anchor seed only
      real(wp), intent(in), optional :: dbranch_shift
      !> Tangent of the adjoint contribution
      real(wp) :: contribution

      !> Tangent of the objective motion the branch adjoint contracts against
      real(wp) :: branch_ddphi

      contribution = dot_product(dw_xyz_pt, dr) + dot_product(w_xyz_pt, ddr) &
                     + eff%w_xi(igrid)*dres%dxi
      if (abs(eff%branch_phi_adj(igrid)) > seed_weight_tol) then
         branch_ddphi = dot_product(dphi1_r, dr) + dot_product(phi1_r, ddr)
         if (present(dbranch_shift)) branch_ddphi = branch_ddphi - dbranch_shift
         contribution = contribution + eff%branch_phi_adj(igrid)*branch_ddphi
      end if
      if (eff%have_wk) then
         contribution = contribution + eff%w_k1(igrid)*dres%dk1 + eff%w_k2(igrid)*dres%dk2
      end if
   end function seed_contribution_tangent

end module moist_cavity_drop_derivatives_seeds
