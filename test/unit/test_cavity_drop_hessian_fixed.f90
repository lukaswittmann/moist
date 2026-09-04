!> Fixed-adjoint half of the DROP surface Hessian
!>
!> [[get_surface_hessian_fixed_drop]] is the first piece of the Hessian scheme
!> that can be checked against production code end to end. Everything below it
!> is verified by finite-differencing a neighbouring primitive; this suite
!> central-differences the *shipped* [[get_surface_gradient_drop]] along a
!> nuclear direction with the surface adjoints held fixed, and compares that
!> against the rank-4 result contracted with the same direction:
!>
!>     sum_{gamma, B} H(:, A, gamma, B) v(gamma, B)
!>         ==  [ g(R + h v) - g(R - h v) ] / (2 h)
!>
!> **The identity holds only where the effective weights are geometry
!> independent.** [[prepare_surface_weights]] folds `w_a` and `w_w` into `w_xi`
!> and `w_f` through geometry-dependent grid data, and derives `branch_phi_adj`
!> from the branch softmax. This fixture therefore never populates `w_a` or
!> `w_w`, and [[assert_frozen_eff]] pins the primal quantities the remaining
!> fold reads -- `branch_count` and `wbranch` -- at every geometry of every
!> stencil, so the suite fails loudly rather than silently if a future fixture
!> starts branching.
!>
!> The frozen adjoints are a pure function of the persistent point id
!> `cavity%numbering`, not of the grid slot, because the grid is filtered and
!> reordered on every rebuild. [[assert_grid_match]] additionally requires the
!> displaced grids to carry exactly the same points with the same owners: a
!> point appearing or vanishing under displacement would put a step
!> discontinuity into the differenced gradient, so that assertion is also the
!> upper bound on the finite-difference step.
!>
!> ## Where the numerical floor comes from
!>
!> `FD_ABS` and `FD_REL` are set to the `1e-10` the project asked for and four
!> of the cases below do not reach it. The reason was measured, not guessed:
!>
!> * It is **not** the projection tolerance. Rebuilding the whole fixture at
!>   `tolerance` = 1e-10, 1e-12 and 1e-14 reproduces the differenced gradient to
!>   three or four significant figures at every step from `3e-3` to `1e-6`
!>   (e.g. 7.93e-8 / 7.81e-8 / 7.82e-8 at `h = 3e-4`). The projection stops
!>   binding somewhere around `1e-10`; at this fixture's `1e-14` it is four
!>   orders past the point where it matters. Loosening it *does* bite -- at
!>   `tolerance = 1e-8` the floor is `5e-7`, and at `1e-6` the identity breaks
!>   outright at `1.7e-3` -- so `PROJ_TOL` earns its value, it just cannot buy
!>   anything more.
!>
!> * It is the **curvature channel**. Driven one at a time, `w_xi`, `w_f`,
!>   `w_xyz` and `w_n` all bottom out at `5e-11 .. 1e-10` absolute near
!>   `h = 1e-4` -- clean `O(h^4)` truncation, `1e-10` is in reach for them.
!>   `w_k1` and `w_k2` never descend at all: `7.4e-8` at `h = 1e-3` rising to
!>   `3.9e-6` at `h = 1e-5`, which is `1/h`, the signature of a fixed noise
!>   floor. Second-differencing the surface gradient at displacements down to
!>   `1e-12` puts that floor at `2e-15` relative for the four clean channels
!>   (a handful of ulp) and `2.5e-11` relative for `w_k1` -- five orders worse.
!>
!> * The amplifier is the near-umbilic point. `compute_curvature` forms the
!>   principal curvatures stably, `disc = hypot(half_diff, S12)`, but the
!>   derivative path re-derives the same discriminant by cancellation,
!>   `disc_curv = sqrt(KM^2 - D)`, and then divides by it:
!>   `d_disc_c = (KM dT - dD) / (2 disc_curv)`. On this fixture the smallest
!>   `|k1 - k2| / |KM|` is `1.3e-6` (SvdW, 10 of 126 points below `1e-4`) and
!>   `9.4e-11` (CFC, one point below the `curv_disc_guard = 1e-10` cutoff
!>   entirely). `eps KM^2 / (2 disc |k1|)` predicts `4e-11` relative for SvdW
!>   and `5e-10` for CFC against `2.5e-11` and `1.9e-9` measured.
!>
!> Reaching `1e-10` on `w_k1` and `w_k2` is therefore a source change -- the
!> stable discriminant in the derivative path, the one `compute_curvature` and
!> `drop_seed_state`'s `sqrt_disc_B` already use -- not a tolerance the test
!> can pick. No parameter setting reaches it.
module test_cavity_drop_hessian_fixed
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_radii, only: default_cpcm_radii
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: fill_legacy_radii

   implicit none(type, external)
   private

   public :: collect_cavity_drop_hessian_fixed

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Level-set model of the fixture
   integer, parameter :: LSF_SVDW = 1, LSF_CFC = 2

   !> Nuclear direction of the contraction
   integer, parameter :: DIR_SINGLE = 1, DIR_MULTI = 2

   !> Surface-adjoint channels this half can differentiate. `w_a` and `w_w` are
   !> absent by construction: their fold into `w_xi` and `w_f` is the one part
   !> of [[prepare_surface_weights]] that moves with the geometry
   integer, parameter :: CH_XI = 1, CH_F = 2, CH_XYZ = 3, CH_N = 4
   integer, parameter :: CH_K1 = 5, CH_K2 = 6, NCHAN = 6

   !> Shared fixture: SvdW blending, Lebedev order and projection settings
   real(wp), parameter :: BLEND_K = 2.5_wp
   real(wp), parameter :: BLEND_3B = 1.0_wp
   integer, parameter :: NUM_LEB = 50
   real(wp), parameter :: PROJ_TOL = 1.0E-14_wp
   integer, parameter :: PROJ_MAXITER = 1000
   integer, parameter :: PROJ_LEVEL = 2
   integer, parameter :: WLEB_PRUNE = 4

   !> Central-difference steps every assertion runs at
   !>
   !> Measured, not assumed. Swept over `3e-3 .. 1e-5` on all four
   !> level-set/direction combinations: the `O(h^4)` stencil below reaches the
   !> round-off floor already at `1e-3`, sits on it through `3e-4`, and starts
   !> climbing again below `1e-4`. Both steps are required, so a value that
   !> agrees at one step only -- the signature of a step on the round-off wall
   !> or above the plateau -- fails.
   !>
   !> There is also a ceiling, and it is enforced rather than assumed:
   !> [[assert_grid_match]] rejects a step large enough to change the grid, and
   !> the stencil reaches out to `2h`.
   real(wp), parameter :: FD_STEPS(2) = [1.0E-3_wp, 3.0E-4_wp]

   !> Finite-difference agreement bounds
   !>
   !> **These are the target, not the measured floor: four of the eight
   !> finite-difference cases below fail at this bound and are meant to.**
   !>
   !> Measured 2026-09-03 on this fixture, worst deviation over both steps.
   !> The first four rows are what the assertion below now reports, absolute and
   !> relative on the *same* component; the six channel rows come from a sweep
   !> and take each maximum independently:
   !>
   !> | case         | absolute | relative | at              |
   !> |--------------|----------|----------|-----------------|
   !> | svdw/single  | 3.2e-08  | 1.6e-08  | atom 2 axis 1, h = 3e-4 |
   !> | svdw/multi   | 2.1e-07  | 1.7e-08  | atom 1 axis 1, h = 1e-3 |
   !> | cfc/single   | 6.7e-07  | 3.7e-06  | atom 2 axis 3, h = 3e-4 |
   !> | cfc/multi    | 3.0e-07  | 7.9e-09  | atom 1 axis 3, h = 1e-3 |
   !> | w_xi         | 1.9e-08  | 7.8e-08  | h = 1e-3 |
   !> | w_f          | 2.1e-07  | 3.0e-08  | h = 1e-3 |
   !> | w_xyz        | 4.0e-10  | 1.8e-10  | h = 1e-3 / 3e-4 |
   !> | w_n          | 3.4e-10  | 4.8e-10  | h = 1e-3 |
   !> | w_k1         | 1.1e-07  | 6.4e-06  | h = 3e-4 |
   !> | w_k2         | 1.7e-07  | 1.2e-07  | h = 3e-4 |
   !>
   !> Everything except `w_k1` and `w_k2` is still on the `O(h^4)` truncation
   !> wing at these two steps and reaches `5e-11 .. 1e-10` absolute at
   !> `h = 1e-4`; `FD_STEPS` is deliberately left where it is so the numbers
   !> above stay comparable with the previous sweep. `w_k1` and `w_k2` grow as
   !> `1/h` all the way from `1e-3` to `1e-5` (7.4e-8 -> 3.9e-6): they are
   !> noise limited, and no step reaches `1e-10`. See the module header.
   real(wp), parameter :: FD_ABS = 1.0E-10_wp
   real(wp), parameter :: FD_REL = 1.0E-10_wp

   !> Below this the differenced gradient carries no information at all and a
   !> comparison against it would pass for the wrong reason
   real(wp), parameter :: VACUITY_THR = 1.0E-4_wp

   !> Hessian symmetry bound. Both triangles are assembled by different code
   !> paths, so this is round-off of the assembly, not of a difference.
   !>
   !> Not a finite-difference tolerance and therefore left where it is. Measured
   !> worst asymmetry on this fixture: `9.0e-12` absolute (SvdW, `max |H| = 12.0`)
   !> and `3.5e-11` absolute / `2.3e-11` scaled (CFC, `max |H| = 17.4`), so it
   !> already carries two to three orders of headroom and would survive `1e-10`.
   real(wp), parameter :: SYM_TOL = 1.0E-8_wp

contains

   !> Collect the suite
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_cavity_drop_hessian_fixed(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("svdw_hvp_single_atom", test_svdw_single), &
                  new_unittest("svdw_hvp_multi_atom", test_svdw_multi), &
                  new_unittest("cfc_hvp_single_atom", test_cfc_single), &
                  new_unittest("cfc_hvp_multi_atom", test_cfc_multi), &
                  new_unittest("single_channels", test_single_channels), &
                  new_unittest("hessian_symmetry", test_symmetry), &
                  new_unittest("frozen_weight_guards", test_frozen_weight_guards), &
                  new_unittest("shape_guard", test_shape_guard) &
                  ]
   end subroutine collect_cavity_drop_hessian_fixed

   !* ================================================================================= *!
   !*                       End-to-end finite-difference tests                          *!
   !* ================================================================================= *!

   !> SvdW, one atom moved along one axis
   !>
   !> @param[out] error Error handle
   subroutine test_svdw_single(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_hvp_fd(LSF_SVDW, DIR_SINGLE, all_channels(), "svdw/single", error)
   end subroutine test_svdw_single

   !> SvdW, every atom moved along a different direction
   !>
   !> @param[out] error Error handle
   subroutine test_svdw_multi(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_hvp_fd(LSF_SVDW, DIR_MULTI, all_channels(), "svdw/multi", error)
   end subroutine test_svdw_multi

   !> CFC, one atom moved along one axis
   !>
   !> @param[out] error Error handle
   subroutine test_cfc_single(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_hvp_fd(LSF_CFC, DIR_SINGLE, all_channels(), "cfc/single", error)
   end subroutine test_cfc_single

   !> CFC, every atom moved along a different direction
   !>
   !> @param[out] error Error handle
   subroutine test_cfc_multi(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_hvp_fd(LSF_CFC, DIR_MULTI, all_channels(), "cfc/multi", error)
   end subroutine test_cfc_multi

   !> Each adjoint channel driven on its own
   !>
   !> A channel the Hessian dropped entirely would still pass the combined
   !> tests above if a larger channel dominated the sum, and the channels of
   !> this fixture do not have comparable magnitudes. Each one therefore gets
   !> its own case with its own anti-vacuity floor.
   !>
   !> @param[out] error Error handle
   subroutine test_single_channels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Channel index
      integer :: ichannel
      !> Channel label
      character(len=12) :: label(NCHAN)

      label = [character(len=12) :: "w_xi", "w_f", "w_xyz", "w_n", "w_k1", "w_k2"]

      do ichannel = 1, NCHAN
         call run_hvp_fd(LSF_SVDW, DIR_MULTI, [ichannel], trim(label(ichannel)), error)
         if (allocated(error)) return
      end do
   end subroutine test_single_channels

   !> Central-difference the shipped surface gradient against the rank-4 result
   !>
   !> @param[in]  lsf_kind Level-set model of the fixture
   !> @param[in]  dir_kind Nuclear direction to contract with
   !> @param[in]  channels Adjoint channels to populate
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine run_hvp_fd(lsf_kind, dir_kind, channels, label, error)
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Nuclear direction
      integer, intent(in) :: dir_kind
      !> Adjoint channels to populate
      integer, intent(in) :: channels(:)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol

      !> Analytic Hessian, its contraction and the differenced reference
      real(wp), allocatable :: hess(:, :, :, :), hv(:, :), fd(:, :, :)
      !> Nuclear direction
      real(wp), allocatable :: vdir(:, :)
      !> Grid and atom extents, loop indices
      integer :: nsph, iatom, iaxis, istep
      !> Deviation bookkeeping
      real(wp) :: diff, worst, worst_rel, ref, bad_ana, bad_ref
      integer :: bad_step, bad_atom, bad_axis

      call fixture_geometry(mol)
      call build_cavity(cavity, ctx, mol, lsf_kind, error)
      if (allocated(error)) return

      nsph = cavity%nsph
      call assert_frozen_eff(cavity, "base geometry", error)
      if (allocated(error)) return

      call frozen_adjoint(cavity, channels, acc, error)
      if (allocated(error)) return

      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      call cavity%get_surface_hessian_fixed(acc, hess, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "fixed-adjoint Hessian failed ("//label//"): "// &
                          cav_error%message)
         return
      end if

      call build_direction(dir_kind, nsph, vdir)

      allocate (hv(ndim, nsph), source=0.0_wp)
      do iatom = 1, nsph
         do iaxis = 1, ndim
            hv(:, :) = hv(:, :) + hess(:, :, iaxis, iatom)*vdir(iaxis, iatom)
         end do
      end do

      !* ------------------------- Central-difference reference ------------------------ *!
      allocate (fd(ndim, nsph, size(FD_STEPS)), source=0.0_wp)
      do istep = 1, size(FD_STEPS)
         call fd_surface_gradient(mol, cavity, lsf_kind, vdir, channels, FD_STEPS(istep), &
                                  fd(:, :, istep), label, error)
         if (allocated(error)) return
      end do

      ! Anti-vacuity: a reference at machine zero would be matched by anything
      if (maxval(abs(fd(:, :, 1))) <= VACUITY_THR) then
         call test_failed(error, "differenced reference is vacuous for "//label// &
                          " (max |dg| = "//to_string(maxval(abs(fd(:, :, 1))))//")")
         return
      end if

      ! The whole grid is scanned before anything is reported: the deviation this
      ! suite exists to expose is the *worst* one, and failing on the first
      ! component over the bound would name an arbitrary early one instead
      worst = 0.0_wp
      worst_rel = 0.0_wp
      bad_step = 0
      do istep = 1, size(FD_STEPS)
         do iatom = 1, nsph
            do iaxis = 1, ndim
               ref = fd(iaxis, iatom, istep)
               diff = abs(hv(iaxis, iatom) - ref)
               if (diff > FD_ABS .and. diff > FD_REL*abs(ref)) then
                  if (diff > worst) then
                     worst = diff
                     worst_rel = diff/max(abs(ref), tiny(1.0_wp))
                     bad_step = istep
                     bad_atom = iatom
                     bad_axis = iaxis
                     bad_ana = hv(iaxis, iatom)
                     bad_ref = ref
                  end if
               end if
            end do
         end do
      end do

      if (bad_step > 0) then
         call test_failed(error, "fixed-adjoint Hessian mismatch for "//label// &
                          ": worst deviation "//to_string(worst)//" absolute, "// &
                          to_string(worst_rel)//" relative, at atom "//to_string(bad_atom)// &
                          " axis "//to_string(bad_axis)//" (h = "// &
                          to_string(FD_STEPS(bad_step))//"): analytic "//to_string(bad_ana)// &
                          " finite difference "//to_string(bad_ref))
         return
      end if

   end subroutine run_hvp_fd

   !> Central difference of the shipped surface gradient along a nuclear direction
   !>
   !> All four displaced geometries carry the *same* frozen adjoint: the weights are
   !> a pure function of `cavity%numbering`, so no explicit mapping is needed
   !> once [[assert_grid_match]] has established that the point set is the same.
   !>
   !> @param[in]  mol      Base structure
   !> @param[in]  ref_cav  Base cavity, used for the grid and weight comparison
   !> @param[in]  lsf_kind Level-set model of the fixture
   !> @param[in]  vdir     Nuclear direction (3, nsph)
   !> @param[in]  channels Adjoint channels to populate
   !> @param[in]  step     Central-difference step
   !> @param[out] deriv    Differenced gradient (3, nsph)
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine fd_surface_gradient(mol, ref_cav, lsf_kind, vdir, channels, step, deriv, &
                                  label, error)
      !> Base structure
      type(structure_type), intent(in) :: mol
      !> Base cavity
      type(cavity_type_drop), intent(in) :: ref_cav
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Nuclear direction
      real(wp), intent(in) :: vdir(:, :)
      !> Adjoint channels to populate
      integer, intent(in) :: channels(:)
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Differenced gradient
      real(wp), intent(out) :: deriv(:, :)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol_disp

      !> Five-point central stencil of the first derivative, `O(h^4)`
      !>
      !> A three-point stencil leaves an `O(h^2)` truncation error that, at a
      !> step small enough to keep the grid stable, is still four orders above
      !> the round-off floor of the gradient -- large enough to set the
      !> comparison tolerance itself. The extra pair of geometries costs two
      !> cavity builds and buys three orders, which is what lets the assertion
      !> below be tight enough to mean something.
      integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
      real(wp), parameter :: COEFF(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp

      real(wp), allocatable :: grad(:, :)
      integer :: iside
      character(len=12) :: side

      deriv = 0.0_wp

      do iside = 1, size(OFFSET)
         write (side, "(a, i0)") "offset ", OFFSET(iside)

         mol_disp = mol
         mol_disp%xyz = mol%xyz + real(OFFSET(iside), wp)*step*vdir

         call build_cavity(cavity, ctx, mol_disp, lsf_kind, error)
         if (allocated(error)) return

         call assert_frozen_eff(cavity, label//" "//trim(side), error)
         if (allocated(error)) return
         call assert_grid_match(ref_cav, cavity, label//" "//trim(side)//" (h = "// &
                                to_string(step)//")", error)
         if (allocated(error)) return

         call frozen_adjoint(cavity, channels, acc, error)
         if (allocated(error)) return

         if (allocated(grad)) deallocate (grad)
         allocate (grad(size(deriv, 1), size(deriv, 2)), source=0.0_wp)
         call cavity%get_surface_gradient(acc, grad, cav_error)
         if (allocated(cav_error)) then
            call test_failed(error, "surface gradient failed ("//label//"): "// &
                             cav_error%message)
            return
         end if

         deriv = deriv + COEFF(iside)*grad/step

         deallocate (cavity)
      end do

   end subroutine fd_surface_gradient

   !* ================================================================================= *!
   !*                              Structural properties                                *!
   !* ================================================================================= *!

   !> Symmetry of the assembled fixed-adjoint block
   !>
   !> With the effective weights frozen this half is the second derivative of
   !> the scalar `sum_p omega_p . Gamma_p`, so the block must be symmetric under
   !> exchanging `(A, alpha)` with `(B, beta)`. Nothing symmetrises it: the two
   !> triangles are produced by different code -- the field row of one atom
   !> against the direction of another -- so this is an independent check on the
   !> assembly rather than a restatement of it.
   !>
   !> @param[out] error Error handle
   subroutine test_symmetry(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol

      real(wp), allocatable :: hess(:, :, :, :)
      integer :: nsph, iatom, jatom, iaxis, jaxis
      real(wp) :: diff, scale

      call fixture_geometry(mol)
      call build_cavity(cavity, ctx, mol, LSF_SVDW, error)
      if (allocated(error)) return

      nsph = cavity%nsph
      call frozen_adjoint(cavity, all_channels(), acc, error)
      if (allocated(error)) return

      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      call cavity%get_surface_hessian_fixed(acc, hess, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "fixed-adjoint Hessian failed: "//cav_error%message)
         return
      end if

      if (maxval(abs(hess)) <= VACUITY_THR) then
         call test_failed(error, "Hessian is vacuous (max |H| = "// &
                          to_string(maxval(abs(hess)))//")")
         return
      end if

      do jatom = 1, nsph
         do jaxis = 1, ndim
            do iatom = 1, nsph
               do iaxis = 1, ndim
                  diff = abs(hess(iaxis, iatom, jaxis, jatom) &
                             - hess(jaxis, jatom, iaxis, iatom))
                  scale = max(abs(hess(iaxis, iatom, jaxis, jatom)), 1.0_wp)
                  if (diff > SYM_TOL*scale) then
                     call test_failed(error, "fixed-adjoint Hessian is asymmetric at ("// &
                                      to_string(iaxis)//", "//to_string(iatom)//") / ("// &
                                      to_string(jaxis)//", "//to_string(jatom)//"): "// &
                                      to_string(hess(iaxis, iatom, jaxis, jatom))//" vs "// &
                                      to_string(hess(jaxis, jatom, iaxis, iatom)))
                     return
                  end if
               end do
            end do
         end do
      end do

   end subroutine test_symmetry

   !* ================================================================================= *!
   !*                                     Guards                                        *!
   !* ================================================================================= *!

   !> A geometry-dependent effective weight must be refused, not dropped
   !>
   !> The area and integration-weight channels fold into `w_xi` and `w_f`
   !> through grid data that moves with the nuclei. Differentiating that fold is
   !> the weight-tangent pass; until it exists the fixed-adjoint half owes the
   !> caller an error rather than a Hessian silently missing those terms.
   !>
   !> @param[out] error Error handle
   subroutine test_frozen_weight_guards(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol

      real(wp), allocatable :: hess(:, :, :, :), ws(:)
      integer :: ichannel

      call fixture_geometry(mol)
      call build_cavity(cavity, ctx, mol, LSF_SVDW, error)
      if (allocated(error)) return

      allocate (hess(ndim, cavity%nsph, ndim, cavity%nsph))
      allocate (ws(cavity%ngrid), source=0.5_wp)

      do ichannel = 1, 2
         call acc%init(cavity%ngrid)
         if (ichannel == 1) then
            call acc%add_surface_weights(cav_error, w_a=ws)
         else
            call acc%add_surface_weights(cav_error, w_w=ws)
         end if
         if (allocated(cav_error)) then
            call test_failed(error, "failed to seed the guard channel: "//cav_error%message)
            return
         end if

         hess = 0.0_wp
         call cavity%get_surface_hessian_fixed(acc, hess, cav_error)
         if (.not. allocated(cav_error)) then
            call test_failed(error, "a geometry-dependent effective weight was accepted"// &
                             " (channel "//to_string(ichannel)//")")
            return
         end if
         deallocate (cav_error)

         ! The accumulator the caller passed in must come back untouched
         if (maxval(abs(hess)) /= 0.0_wp) then
            call test_failed(error, "the Hessian accumulator was written on a failure")
            return
         end if
      end do

   end subroutine test_frozen_weight_guards

   !> A mis-shaped accumulator must be refused
   !>
   !> @param[out] error Error handle
   subroutine test_shape_guard(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol

      real(wp), allocatable :: hess(:, :, :, :)

      call fixture_geometry(mol)
      call build_cavity(cavity, ctx, mol, LSF_SVDW, error)
      if (allocated(error)) return

      call frozen_adjoint(cavity, all_channels(), acc, error)
      if (allocated(error)) return

      allocate (hess(ndim, cavity%nsph + 1, ndim, cavity%nsph), source=0.0_wp)
      call cavity%get_surface_hessian_fixed(acc, hess, cav_error)
      if (.not. allocated(cav_error)) then
         call test_failed(error, "a mis-shaped Hessian accumulator was accepted")
         return
      end if

   end subroutine test_shape_guard

   !* ================================================================================= *!
   !*                             Preconditions of the identity                         *!
   !* ================================================================================= *!

   !> Assert that this geometry's effective weights are geometry independent
   !>
   !> `w_a` and `w_w` are never populated by this fixture, so the only remaining
   !> geometry-dependent output of [[prepare_surface_weights]] is
   !> `branch_phi_adj`. [[compute_branch_phi_adj]] returns identically zero
   !> unless some anchor group carries more than one branch, and `wbranch` is
   !> exactly one for every point outside such a group. Both are asserted, on
   !> the primal data, because `eff` itself is not reachable from a test.
   !>
   !> @param[in]  cavity Cavity to inspect
   !> @param[in]  label  Human-readable geometry description
   !> @param[out] error  Error handle
   subroutine assert_frozen_eff(cavity, label, error)
      !> Cavity to inspect
      type(cavity_type_drop), intent(in) :: cavity
      !> Geometry description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer :: igrid

      if (allocated(cavity%branch_count)) then
         if (any(cavity%branch_count(1:cavity%ngrid) > 1)) then
            call test_failed(error, "the fixture branched at "//label// &
                             "; the branch adjoint is then geometry dependent and the"// &
                             " frozen-weight identity under test does not hold"// &
                             " (max branch_count "// &
                             to_string(maxval(cavity%branch_count(1:cavity%ngrid)))//")")
            return
         end if
      end if

      do igrid = 1, cavity%ngrid
         if (cavity%wbranch(igrid) /= 1.0_wp) then
            call test_failed(error, "wbranch is not exactly one at "//label// &
                             ", grid point "//to_string(igrid)//" ("// &
                             to_string(cavity%wbranch(igrid))//"); the branch fold would"// &
                             " then move with the geometry")
            return
         end if
      end do

   end subroutine assert_frozen_eff

   !> Assert that two geometries carry the very same grid points
   !>
   !> The frozen adjoint is keyed on `numbering`, so a point that appears or
   !> vanishes between the two displaced geometries would silently add a step
   !> to the differenced gradient. Owners are compared as well, because the
   !> switching channel and the anchor channel both key on them. This assertion
   !> is what bounds the finite-difference step from above.
   !>
   !> @param[in]  ref   Reference cavity
   !> @param[in]  cav   Displaced cavity
   !> @param[in]  label Human-readable geometry description
   !> @param[out] error Error handle
   subroutine assert_grid_match(ref, cav, label, error)
      !> Reference cavity
      type(cavity_type_drop), intent(in) :: ref
      !> Displaced cavity
      type(cavity_type_drop), intent(in) :: cav
      !> Geometry description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer :: igrid

      if (cav%ngrid /= ref%ngrid) then
         call test_failed(error, "the grid changed size at "//label//" ("// &
                          to_string(ref%ngrid)//" -> "//to_string(cav%ngrid)// &
                          "); the step is above the branch-enumeration ceiling")
         return
      end if

      do igrid = 1, ref%ngrid
         if (cav%numbering(igrid) /= ref%numbering(igrid)) then
            call test_failed(error, "the grid was reordered or repopulated at "//label// &
                             ", slot "//to_string(igrid)//" ("//to_string(ref%numbering(igrid))// &
                             " -> "//to_string(cav%numbering(igrid))//")")
            return
         end if
         if (cav%owner(igrid) /= ref%owner(igrid)) then
            call test_failed(error, "the owner sphere changed at "//label// &
                             ", slot "//to_string(igrid)//" ("//to_string(ref%owner(igrid))// &
                             " -> "//to_string(cav%owner(igrid))//")")
            return
         end if
      end do

   end subroutine assert_grid_match

   !* ================================================================================= *!
   !*                                    Fixture                                        *!
   !* ================================================================================= *!

   !> Every channel this half can differentiate
   !>
   !> @return Channel identifiers
   pure function all_channels() result(channels)
      !> Channel identifiers
      integer :: channels(NCHAN)

      channels = [CH_XI, CH_F, CH_XYZ, CH_N, CH_K1, CH_K2]
   end function all_channels

   !> Populate the requested adjoint channels
   !>
   !> `w_a` and `w_w` are never reachable from here: they are the two channels
   !> whose fold is geometry dependent. Each weight is a pure function of the
   !> *persistent* point id, so the same adjoint is reproduced on a displaced
   !> grid without an explicit mapping.
   !>
   !> @param[in]  cavity   Cavity supplying the grid
   !> @param[in]  channels Channel identifiers to populate
   !> @param[out] acc      Surface-adjoint accumulator
   !> @param[out] error    Error handle
   subroutine frozen_adjoint(cavity, channels, acc, error)
      !> Cavity supplying the grid
      type(cavity_type_drop), intent(in) :: cavity
      !> Channel identifiers to populate
      integer, intent(in) :: channels(:)
      !> Surface-adjoint accumulator
      type(cavity_surface_adjoint_type), intent(out) :: acc
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: add_error
      real(wp), allocatable :: ws(:), wv(:, :)
      integer :: ngrid, igrid, iaxis, ichannel

      ngrid = cavity%ngrid
      allocate (ws(ngrid), wv(ndim, ngrid))

      call acc%init(ngrid)
      do ichannel = 1, size(channels)
         do igrid = 1, ngrid
            ws(igrid) = point_weight(cavity%numbering(igrid), channels(ichannel))
            do iaxis = 1, ndim
               wv(iaxis, igrid) = point_weight(cavity%numbering(igrid), &
                                               NCHAN*channels(ichannel) + iaxis)
            end do
         end do

         select case (channels(ichannel))
         case (CH_XI)
            call acc%add_surface_weights(add_error, w_xi=ws)
         case (CH_F)
            call acc%add_surface_weights(add_error, w_f=ws)
         case (CH_XYZ)
            call acc%add_surface_weights(add_error, w_xyz=wv)
         case (CH_N)
            call acc%add_surface_weights(add_error, w_n=wv)
         case (CH_K1)
            call acc%add_surface_weights(add_error, w_k1=ws)
         case (CH_K2)
            call acc%add_surface_weights(add_error, w_k2=ws)
         end select
         if (allocated(add_error)) then
            call test_failed(error, "failed to seed the surface adjoint: "//add_error%message)
            return
         end if
      end do

   end subroutine frozen_adjoint

   !> Reproducible adjoint weight of one persistent point and channel
   !>
   !> Varies across the grid so that a term which happens to cancel for uniform
   !> weights still shows up, and is bounded well away from zero so that no
   !> channel is accidentally switched off.
   !>
   !> @param[in] id      Persistent point id, `cavity%numbering`
   !> @param[in] channel Channel selector
   !> @return            Adjoint weight
   pure function point_weight(id, channel) result(w)
      !> Persistent point id
      integer, intent(in) :: id
      !> Channel selector
      integer, intent(in) :: channel
      !> Adjoint weight
      real(wp) :: w

      w = 0.60_wp + 0.35_wp*sin(0.7_wp*real(id, wp) + 1.3_wp*real(channel, wp)) &
          + 0.11_wp*cos(0.23_wp*real(id, wp)*real(channel + 2, wp))
   end function point_weight

   !> Nuclear direction of the contraction
   !>
   !> @param[in]  dir_kind Direction selector
   !> @param[in]  nsph     Number of spheres
   !> @param[out] vdir     Nuclear direction (3, nsph)
   subroutine build_direction(dir_kind, nsph, vdir)
      !> Direction selector
      integer, intent(in) :: dir_kind
      !> Number of spheres
      integer, intent(in) :: nsph
      !> Nuclear direction
      real(wp), allocatable, intent(out) :: vdir(:, :)

      integer :: iatom, iaxis

      allocate (vdir(ndim, nsph), source=0.0_wp)

      select case (dir_kind)
      case (DIR_SINGLE)
         ! One atom, one axis: the sparsest column of the block, and the one a
         ! wrong influence set would zero out entirely
         vdir(3, 2) = 1.0_wp
      case (DIR_MULTI)
         ! Every atom moved along its own direction, none of them a translation
         do iatom = 1, nsph
            do iaxis = 1, ndim
               vdir(iaxis, iatom) = sin(1.1_wp*real(iatom, wp) + 0.6_wp*real(iaxis, wp))
            end do
         end do
      end select

   end subroutine build_direction

   !> Shared fixture geometry
   !>
   !> Asymmetric on purpose: a symmetric geometry drives the multistart
   !> projection into sibling branches, and a branched grid violates the
   !> frozen-weight precondition of the whole suite.
   !>
   !> @param[out] mol Structure
   subroutine fixture_geometry(mol)
      !> Structure
      type(structure_type), intent(out) :: mol

      call new(mol, [8, 6, 1], reshape([ &
                                       0.00_wp, 0.00_wp, 0.00_wp, &
                                       0.00_wp, 0.00_wp, 4.60_wp, &
                                       2.60_wp, 0.40_wp, -1.10_wp], [3, 3]))
   end subroutine fixture_geometry

   !> Build the DROP cavity for a structure and a level-set model
   !>
   !> @param[out]   cavity   Constructed cavity
   !> @param[inout] ctx      Run context borrowed by the cavity; must outlive it
   !> @param[in]    mol      Structure to build on
   !> @param[in]    lsf_kind Level-set model
   !> @param[out]   error    Error handle
   subroutine build_cavity(cavity, ctx, mol, lsf_kind, error)
      !> Constructed cavity
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(inout) :: ctx
      !> Structure to build on
      type(structure_type), intent(in) :: mol
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: radii(:)
      type(mctc_error), allocatable :: cav_error

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      allocate (cavity)
      call new_context(ctx, verbosity=0)
      select case (lsf_kind)
      case (LSF_SVDW)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=BLEND_K, blend_3b=BLEND_3B)
            call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                                 tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                                 proj_level=PROJ_LEVEL, wleb_prune_level=WLEB_PRUNE, &
                                 radius_model=default_cpcm_radii(), &
                                 lsf_model=svdw_template, error=cav_error)
         end block
      case default
         block
            type(moist_cavity_drop_lsf_cfc_type) :: cfc_template
            call cfc_template%new()
            call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                                 tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                                 proj_level=PROJ_LEVEL, wleb_prune_level=WLEB_PRUNE, &
                                 radius_model=default_cpcm_radii(), &
                                 lsf_model=cfc_template, error=cav_error)
         end block
      end select
      if (allocated(cav_error)) then
         call test_failed(error, "failed to initialize cavity: "//cav_error%message)
         return
      end if

      ! Curvature and normals are surface observables this suite drives, so the
      ! cavity has to be asked for them
      call cavity%properties(do_fine=.true.)

      call cavity%update(mol, error=cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "failed to build cavity: "//cav_error%message)
         return
      end if
      if (cavity%ngrid == 0) then
         call test_failed(error, "empty grid")
         return
      end if

   end subroutine build_cavity

end module test_cavity_drop_hessian_fixed
