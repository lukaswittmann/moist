!> Pass 1 of the DROP Hessian: the forward surface tangent
!>
!> [[get_surface_tangent_drop]] returns the directional derivatives of the four
!> grid scalars pass 2 consumes -- `a`, `wleb`, `xi0` and `wbranch` -- along a
!> batch of nuclear directions. The ground truth here is the **primal** surface
!> itself: displace the nuclei along a direction, rebuild the cavity, and
!> central-difference the four cavity arrays per grid point.
!>
!> ## What the fixtures are for
!>
!> Two geometries, because the interesting term only exists on one of them:
!>
!>   * `FIX_PLAIN`, an asymmetric OCH triple, never branches. `wbranch` is
!>     exactly one at every point and `d_wbranch` must be identically zero.
!>     Both level-set models run on it.
!>   * `FIX_CROSS`, the five-carbon cross of [[get_test_cross]] at
!>     `proj_level = 7` with a softened branch softmax (`s = 0.5`; at the
!>     production 0.05 the prune in `filter.f90` collapses every group back to
!>     one branch), does branch. That makes `d_wbranch` nonzero and, through
!>     `d_wleb = res%dwleb + d_wbranch * wleb / wbranch`, makes the branch term
!>     of the tangent load bearing. SvdW only -- the point of the fixture is the
!>     branching, not the level set.
!>
!> Both directions of the batch are passed in **one** call, so the per-point
!> batched KKT solve is exercised rather than a degenerate `ndir = 1` path.
!>
!> ## The `wbranch` trap, and the mutation that proves it was handled
!>
!> [[apply_seed]] freezes `wbranch` in `res%dwleb` on purpose -- in first-order
!> reverse mode the branch weight's motion travels through `branch_phi_adj`
!> instead -- so a forward tangent has to add `d_wbranch * wleb / wbranch` back
!> by hand. That term was deleted from the implementation and the suite re-run;
!> the signature this suite is shaped to catch is what came back:
!>
!>   * `svdw_cross_branching_fd` fails, first on `d_a`: worst deviation
!>     `2.5e-3` absolute and `2.4e-2` relative at grid point 92 (analytic
!>     `-0.1076765`, finite difference `-0.1051619`). Per channel, the worst
!>     absolute deviations are `2.5e-3` on `d_a`, `1.8e-4` on `d_wleb` and
!>     `8.9e+2` on `d_xi0`, at a worst relative deviation of `2.17` -- a 217%
!>     error, seven to twelve orders above the `1e-11` the same channels reach
!>     with the term in place;
!>   * those numbers are *identical at every step of the sweep*, from `6e-4`
!>     down to `1.5e-4`. A step-independent deviation is a missing term, not a
!>     finite-difference artefact, and that is the cleanest part of the
!>     signature;
!>   * `d_wbranch` still passes, unchanged at `4e-13 .. 1.5e-12`: it is
!>     produced by the branch stage, not by the dropped fold;
!>   * both `FIX_PLAIN` cases still pass, with deviations identical to the
!>     unmutated run in every digit, because `d_wbranch` is identically zero
!>     there. So does `single_branch_wbranch_is_zero`.
!>
!> The two exact identities below also still pass with the term dropped: they
!> tie `d_xi0` and `d_a` to whatever `d_wleb` says, and say nothing about
!> whether `d_wleb` is right. Only the finite difference does, and only on a
!> branched fixture -- which is why one exists.
!>
!> ## Step and tolerance
!>
!> `FD_STEPS` is measured, not assumed; the sweep is in the comment on that
!> parameter. The grid is guarded at every stencil geometry
!> ([[assert_grid_match]]): `numbering`, `owner`, `branch_count` and
!> `anchor_id` must all be slot-wise identical to the base cavity, so a step
!> large enough to re-enumerate branches or to move a point through a filter
!> threshold fails loudly instead of quietly putting a step discontinuity into
!> the reference.
module test_cavity_drop_tangent_forward
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_gaussian, only: iswig_workspace_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_radii, only: default_cpcm_radii
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: fill_legacy_radii, get_test_cross

   implicit none(type, external)
   private

   public :: collect_cavity_drop_tangent_forward

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Level-set model of a fixture
   integer, parameter :: LSF_SVDW = 1, LSF_CFC = 2

   !> Geometry of a fixture
   integer, parameter :: FIX_PLAIN = 1, FIX_CROSS = 2

   !> Directions pushed through in one call: one sparse, one dense
   integer, parameter :: NDIR = 2

   !> Output channels, in the order [[compare_channels]] walks them
   integer, parameter :: CH_A = 1, CH_WLEB = 2, CH_XI = 3, CH_WBRANCH = 4
   integer, parameter :: NCHAN = 4

   !> Shared fixture settings
   real(wp), parameter :: BLEND_K = 2.5_wp
   real(wp), parameter :: BLEND_3B = 1.0_wp
   integer, parameter :: NUM_LEB = 50
   real(wp), parameter :: PROJ_TOL = 1.0E-14_wp
   integer, parameter :: PROJ_MAXITER = 1000
   integer, parameter :: PROJ_LEVEL = 2
   integer, parameter :: WLEB_PRUNE = 4

   !> Branching fixture: multistart level and softmax temperature
   real(wp), parameter :: CROSS_BLEND_K = 1.0_wp
   integer, parameter :: CROSS_PROJ_LEVEL = 7
   real(wp), parameter :: CROSS_BRANCH_S = 0.5_wp

   !> Production softmax temperature, used by the non-branching fixture
   real(wp), parameter :: PLAIN_BRANCH_S = 0.0025_wp

   !> Central-difference steps every finite-difference assertion runs at
   !>
   !> Measured, over all four channels and both directions. The pair reported
   !> per case is `worst absolute / worst relative`, the relative one taken over
   !> exactly the components that already exceed `FD_ABS` -- that is, over the
   !> components that decide the test:
   !>
   !> | h      | svdw/plain        | cfc/plain         | svdw/cross        |
   !> |--------|-------------------|-------------------|-------------------|
   !> | 3e-3   | 1.2e-6  / 1.7e-6  | 2.2e-6  / 1.7e-6  | 2.5e-5  / 2.0e-7  |
   !> | 1e-3   | 1.5e-8  / 6.8e-9  | 2.7e-8  / 1.4e-8  | 3.0e-7  / 2.4e-9  |
   !> | 6e-4   | 1.9e-9  / 8.9e-10 | 3.5e-9  / 1.8e-9  | 4.0e-8  / 3.2e-10 |
   !> | 4e-4   | 3.8e-10 / 5.3e-11 | 6.8e-10 / 2.9e-10 | 9.7e-9  / 7.3e-11 |
   !> | 3e-4   | 1.2e-10 / 1.6e-11 | 2.1e-10 / 4.6e-11 | 6.3e-9  / 5.5e-11 |
   !> | 2e-4   | 3.4e-11 / 0       | 4.2e-11 / 0       | 5.4e-9  / 4.9e-11 |
   !> | 1.5e-4 | 1.7e-11 / 0       | 2.0e-11 / 0       | 1.1e-8  / 4.0e-11 |
   !> | 1e-4   | 3.5e-11 / 0       | 5.1e-11 / 0       | 2.2e-8  / 2.5e-10 |
   !> | 3e-5   | 1.4e-10 / 0       | 1.1e-10 / 6.2e-11 | 3.9e-8  / 1.8e-8  |
   !> | 1e-5   | 2.9e-10 / 3.0e-6  | 4.7e-10 / 2.2e-5  | 1.5e-7  / 2.3e-7  |
   !>
   !> A clean bowl: `h^4` truncation on the way down -- the absolute column
   !> falls by a factor of ~123 from `1e-3` to `3e-4`, which is `(10/3)^4` --
   !> and `eps_primal/h` round-off on the way up. Both bounds are met over
   !> `3e-4 .. 1.5e-4`, and the two steps below are the ends of that window, a
   !> factor of two apart. Both are required, so a value that agrees at one step
   !> only -- the signature of a step sitting on the round-off wall -- fails.
   !>
   !> `svdw/cross` is the case whose absolute column never gets small, and the
   !> reason is `d_xi0`, not the tangent: with `xi0 = swx/(R sqrt(wleb))` and a
   !> pruning floor at `wleb_cut = 1e-6`, the smallest-weight points of the
   !> branched grid carry `xi0` three orders above the typical one, and their
   !> `d_xi0` with them. Their *relative* accuracy is the same `5e-11` as
   !> everywhere else.
   real(wp), parameter :: FD_STEPS(2) = [3.0E-4_wp, 1.5E-4_wp]

   !> Finite-difference agreement bounds
   !>
   !> The project target, `1e-10` absolute and `1e-10` relative, held: a
   !> component fails only when it misses *both*, and over the two steps below
   !> nothing does. The margin is not large -- `4.6e-11` against the relative
   !> bound for `cfc/plain` at `3e-4`, `4.0e-11` for `svdw/cross` at `1.5e-4`,
   !> so between two and three times -- and that is a statement about the
   !> reference, not about the tangent.
   !>
   !> What sets the floor is the round-off of the *primal*, not the stencil.
   !> The projected point is the solution of an iterative minimization stopped
   !> at `PROJ_TOL = 1e-14`, so `a`, `wleb` and `xi0` carry a few times `1e-14`
   !> of irreducible noise; a central difference divides it by `h`, which puts
   !> `~3e-11` under every quotient at `h = 3e-4` and grows it as `1/h` below
   !> that. The `h^4` truncation is already under that floor at the same step,
   !> which is why the bowl is flat rather than pointed. Tightening the
   !> tolerance further would need a tighter projection, not a better stencil.
   real(wp), parameter :: FD_ABS = 1.0E-10_wp
   real(wp), parameter :: FD_REL = 1.0E-10_wp

   !> Exactness bound of the two algebraic identities
   !>
   !> These are not finite differences: they relate three outputs of one call
   !> through closed-form relations the implementation is built on, so anything
   !> above assembly round-off is a real defect.
   real(wp), parameter :: ID_TOL = 1.0E-12_wp

   !> Below this a channel carries no information and its comparison would pass
   !> for the wrong reason
   real(wp), parameter :: VACUITY_THR(NCHAN) = &
                          [1.0E-3_wp, 1.0E-4_wp, 1.0E-3_wp, 1.0E-4_wp]

contains

   !> Collect the suite
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_cavity_drop_tangent_forward(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("svdw_plain_fd", test_svdw_plain), &
                  new_unittest("cfc_plain_fd", test_cfc_plain), &
                  new_unittest("svdw_cross_branching_fd", test_svdw_cross), &
                  new_unittest("width_identity", test_width_identity), &
                  new_unittest("area_identity", test_area_identity), &
                  new_unittest("single_branch_wbranch_is_zero", test_single_branch), &
                  new_unittest("shape_guards", test_shape_guards) &
                  ]
   end subroutine collect_cavity_drop_tangent_forward

   !* ================================================================================= *!
   !*                        End-to-end finite-difference tests                         *!
   !* ================================================================================= *!

   !> SvdW on the non-branching fixture
   !>
   !> @param[out] error Error handle
   subroutine test_svdw_plain(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_tangent_fd(FIX_PLAIN, LSF_SVDW, "svdw/plain", error)
   end subroutine test_svdw_plain

   !> CFC on the non-branching fixture
   !>
   !> @param[out] error Error handle
   subroutine test_cfc_plain(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_tangent_fd(FIX_PLAIN, LSF_CFC, "cfc/plain", error)
   end subroutine test_cfc_plain

   !> SvdW on the branching fixture
   !>
   !> The only case in which `d_wbranch` is nonzero, and therefore the only one
   !> that can see the `wbranch` term of `d_wleb`.
   !>
   !> @param[out] error Error handle
   subroutine test_svdw_cross(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_tangent_fd(FIX_CROSS, LSF_SVDW, "svdw/cross", error)
   end subroutine test_svdw_cross

   !> Central-difference the primal surface against the analytic tangent
   !>
   !> @param[in]  fix_kind Geometry of the fixture
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine run_tangent_fd(fix_kind, lsf_kind, label, error)
      !> Geometry of the fixture
      integer, intent(in) :: fix_kind
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol

      !> Analytic tangent, one column per direction
      real(wp), allocatable :: tan_out(:, :, :)
      !> Differenced reference of one direction, `(ngrid, NCHAN)`
      real(wp), allocatable :: fd(:, :)
      !> Nuclear directions
      real(wp), allocatable :: dirs(:, :, :)
      !> Grid extent, direction and step indices
      integer :: ngrid, idir, istep

      call fixture_geometry(fix_kind, mol)
      call build_cavity(cavity, ctx, mol, fix_kind, lsf_kind, error)
      if (allocated(error)) return

      call assert_branching(cavity, fix_kind, "base geometry", error)
      if (allocated(error)) return

      ngrid = cavity%ngrid
      call build_directions(cavity%nsph, dirs)

      allocate (tan_out(ngrid, NDIR, NCHAN), source=0.0_wp)
      call cavity%get_surface_tangent(dirs, tan_out(:, :, CH_A), tan_out(:, :, CH_WLEB), &
                                      tan_out(:, :, CH_XI), tan_out(:, :, CH_WBRANCH), &
                                      cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "surface tangent failed ("//label//"): "//cav_error%message)
         return
      end if

      call assert_live(tan_out, fix_kind, label, error)
      if (allocated(error)) return

      allocate (fd(ngrid, NCHAN))
      do istep = 1, size(FD_STEPS)
         do idir = 1, NDIR
            call fd_primal(mol, cavity, fix_kind, lsf_kind, dirs(:, :, idir), &
                           FD_STEPS(istep), fd, label, error)
            if (allocated(error)) return
            call compare_channels(tan_out(:, idir, :), fd, FD_STEPS(istep), idir, &
                                  label, error)
            if (allocated(error)) return
         end do
      end do

   end subroutine run_tangent_fd

   !> Five-point central difference of the four primal grid scalars
   !>
   !> `O(h^4)`. A three-point stencil leaves a truncation error that, at a step
   !> small enough to keep the branch enumeration stable, still sits three
   !> orders above the round-off floor of the primal -- large enough to set the
   !> comparison tolerance itself and to put the `1e-10` target out of reach.
   !>
   !> @param[in]  mol      Base structure
   !> @param[in]  ref_cav  Base cavity, for the grid comparison
   !> @param[in]  fix_kind Geometry of the fixture
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  vdir     Nuclear direction `(3, nsph)`
   !> @param[in]  step     Central-difference step
   !> @param[out] fd       Differenced primal `(ngrid, NCHAN)`
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine fd_primal(mol, ref_cav, fix_kind, lsf_kind, vdir, step, fd, label, error)
      !> Base structure
      type(structure_type), intent(in) :: mol
      !> Base cavity
      type(cavity_type_drop), intent(in) :: ref_cav
      !> Geometry of the fixture
      integer, intent(in) :: fix_kind
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Nuclear direction
      real(wp), intent(in) :: vdir(:, :)
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Differenced primal
      real(wp), intent(out) :: fd(:, :)
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Five-point central stencil of the first derivative
      integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
      real(wp), parameter :: COEFF(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol_disp
      integer :: iside, ngrid
      character(len=32) :: side

      fd = 0.0_wp
      ngrid = ref_cav%ngrid

      do iside = 1, size(OFFSET)
         write (side, "(a, i0, a, es9.2, a)") "offset ", OFFSET(iside), " (h = ", step, ")"

         mol_disp = mol
         mol_disp%xyz = mol%xyz + real(OFFSET(iside), wp)*step*vdir

         call build_cavity(cavity, ctx, mol_disp, fix_kind, lsf_kind, error)
         if (allocated(error)) return

         call assert_grid_match(ref_cav, cavity, label//" "//trim(side), error)
         if (allocated(error)) return

         fd(:, CH_A) = fd(:, CH_A) + COEFF(iside)*cavity%a(1:ngrid)/step
         fd(:, CH_WLEB) = fd(:, CH_WLEB) + COEFF(iside)*cavity%wleb(1:ngrid)/step
         fd(:, CH_XI) = fd(:, CH_XI) + COEFF(iside)*cavity%xi0(1:ngrid)/step
         fd(:, CH_WBRANCH) = fd(:, CH_WBRANCH) + COEFF(iside)*cavity%wbranch(1:ngrid)/step

         deallocate (cavity)
      end do

   end subroutine fd_primal

   !> Compare one direction's analytic tangent against the differenced primal
   !>
   !> @param[in]  tangent Analytic tangent of this direction `(ngrid, NCHAN)`
   !> @param[in]  fd      Differenced primal `(ngrid, NCHAN)`
   !> @param[in]  step    Central-difference step behind `fd`
   !> @param[in]  idir    Direction index, for the diagnostic
   !> @param[in]  label   Human-readable case description
   !> @param[out] error   Error handle
   subroutine compare_channels(tangent, fd, step, idir, label, error)
      !> Analytic tangent of this direction
      real(wp), intent(in) :: tangent(:, :)
      !> Differenced primal
      real(wp), intent(in) :: fd(:, :)
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Direction index
      integer, intent(in) :: idir
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Worst offender of the channel and its coordinates
      integer :: ichannel, igrid, iworst
      real(wp) :: diff, ref, rel, worst_abs, worst_rel
      logical :: failed

      do ichannel = 1, NCHAN
         failed = .false.
         worst_abs = 0.0_wp
         worst_rel = 0.0_wp
         iworst = 0
         do igrid = 1, size(fd, 1)
            ref = fd(igrid, ichannel)
            diff = abs(tangent(igrid, ichannel) - ref)
            rel = diff/max(abs(ref), tiny(1.0_wp))
            ! A component has to miss both bounds to be a failure: the absolute
            ! one alone would condemn the large `d_xi0` of a low-weight point,
            ! the relative one alone a component that is numerically zero
            if (diff > FD_ABS .and. diff > FD_REL*abs(ref)) then
               failed = .true.
               if (diff > worst_abs) then
                  worst_abs = diff
                  worst_rel = rel
                  iworst = igrid
               end if
            end if
         end do

         if (failed) then
            call test_failed(error, "surface tangent mismatch for "//label//" "// &
                             trim(channel_name(ichannel))//", direction "// &
                             to_string(idir)//" (h = "//to_string(step)// &
                             "): worst deviation "//to_string(worst_abs)// &
                             " absolute, "//to_string(worst_rel)//" relative, at"// &
                             " grid point "//to_string(iworst)//": analytic "// &
                             to_string(tangent(iworst, ichannel))//" finite difference "// &
                             to_string(fd(iworst, ichannel)))
            return
         end if
      end do

   end subroutine compare_channels


   !* ================================================================================= *!
   !*                                 Exact identities                                  *!
   !* ================================================================================= *!

   !> `xi0 = swx/(R sqrt(wleb))`, hence `d_wleb = -2 wleb d_xi0 / xi0`
   !>
   !> Checked on both fixtures, so the relation is asserted with the branch term
   !> present in `d_wleb` as well as without it.
   !>
   !> @param[out] error Error handle
   subroutine test_width_identity(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      real(wp), allocatable :: dirs(:, :, :), tan_out(:, :, :)
      integer :: ifix, idir, igrid, nlive
      real(wp) :: implied, diff

      do ifix = FIX_PLAIN, FIX_CROSS
         call fixture_geometry(ifix, mol)
         call build_cavity(cavity, ctx, mol, ifix, LSF_SVDW, error)
         if (allocated(error)) return
         call build_directions(cavity%nsph, dirs)
         call tangent_of(cavity, dirs, tan_out, error)
         if (allocated(error)) return

         nlive = 0
         do idir = 1, NDIR
            do igrid = 1, cavity%ngrid
               if (cavity%wleb(igrid) <= 0.0_wp .or. cavity%xi0(igrid) <= 0.0_wp) cycle
               implied = -2.0_wp*cavity%wleb(igrid)*tan_out(igrid, idir, CH_XI) &
                         /cavity%xi0(igrid)
               diff = abs(implied - tan_out(igrid, idir, CH_WLEB))
               if (abs(tan_out(igrid, idir, CH_WLEB)) > VACUITY_THR(CH_WLEB)) nlive = nlive + 1
               if (diff > ID_TOL*max(1.0_wp, abs(tan_out(igrid, idir, CH_WLEB)))) then
                  call test_failed(error, "width identity violated on fixture "// &
                                   to_string(ifix)//" at grid point "//to_string(igrid)// &
                                   ", direction "//to_string(idir)//": d_wleb "// &
                                   to_string(tan_out(igrid, idir, CH_WLEB))// &
                                   " vs -2 wleb d_xi0/xi0 "//to_string(implied))
                  return
               end if
            end do
         end do

         if (nlive == 0) then
            call test_failed(error, "width identity is vacuous on fixture "//to_string(ifix))
            return
         end if

         deallocate (cavity, dirs, tan_out)
      end do

   end subroutine test_width_identity

   !> `a = R^2 f wleb`, hence `d_a = R^2 (wleb d_f + f d_wleb)`
   !>
   !> `d_f` is rebuilt here from the shipped iSwiG gradient rather than taken
   !> from the routine under test, so this is a genuine check that the
   !> switching channel is present, uses the owner radius and multiplies the
   !> *final* `wleb`; only the contraction order is shared with the
   !> implementation, which is why it holds at assembly round-off.
   !>
   !> @param[out] error Error handle
   subroutine test_area_identity(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(iswig_workspace_type) :: work
      real(wp), allocatable :: dirs(:, :, :), tan_out(:, :, :), swi_grad(:, :)
      integer :: ifix, idir, igrid, nlive, iown
      real(wp) :: d_f, r_own, implied, diff

      do ifix = FIX_PLAIN, FIX_CROSS
         call fixture_geometry(ifix, mol)
         call build_cavity(cavity, ctx, mol, ifix, LSF_SVDW, error)
         if (allocated(error)) return
         call build_directions(cavity%nsph, dirs)
         call tangent_of(cavity, dirs, tan_out, error)
         if (allocated(error)) return

         call work%init(cavity%iswig)
         allocate (swi_grad(ndim, cavity%nsph))

         nlive = 0
         do igrid = 1, cavity%ngrid
            iown = cavity%owner(igrid)
            r_own = cavity%radii(iown)
            call cavity%iswig%swi1_rA(cavity%anchorxyz(:, igrid), iown, &
                                      cavity%anchor_xi0(igrid), work, swi_grad)
            do idir = 1, NDIR
               d_f = sum(swi_grad*dirs(:, :, idir))
               implied = r_own*r_own*(cavity%wleb(igrid)*d_f &
                                      + cavity%f(igrid)*tan_out(igrid, idir, CH_WLEB))
               diff = abs(implied - tan_out(igrid, idir, CH_A))
               if (abs(tan_out(igrid, idir, CH_A)) > VACUITY_THR(CH_A)) nlive = nlive + 1
               if (diff > ID_TOL*max(1.0_wp, abs(tan_out(igrid, idir, CH_A)))) then
                  call test_failed(error, "area identity violated on fixture "// &
                                   to_string(ifix)//" at grid point "//to_string(igrid)// &
                                   ", direction "//to_string(idir)//": d_a "// &
                                   to_string(tan_out(igrid, idir, CH_A))// &
                                   " vs R^2 (wleb d_f + f d_wleb) "//to_string(implied))
                  return
               end if
            end do
         end do

         if (nlive == 0) then
            call test_failed(error, "area identity is vacuous on fixture "//to_string(ifix))
            return
         end if

         call work%destroy()
         deallocate (cavity, dirs, tan_out, swi_grad)
      end do

   end subroutine test_area_identity

   !* ================================================================================= *!
   !*                              Structural properties                                *!
   !* ================================================================================= *!

   !> Without a branched anchor, `d_wbranch` must be exactly zero
   !>
   !> `wbranch` is the constant one on every unbranched point, so this is an
   !> identity rather than a tolerance. It is the companion of the branching
   !> case: together they are the signature that separates a correct branch
   !> stage from one that leaks into single-branch points.
   !>
   !> @param[out] error Error handle
   subroutine test_single_branch(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      real(wp), allocatable :: dirs(:, :, :), tan_out(:, :, :)

      call fixture_geometry(FIX_PLAIN, mol)
      call build_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error)
      if (allocated(error)) return
      call assert_branching(cavity, FIX_PLAIN, "single-branch fixture", error)
      if (allocated(error)) return

      call build_directions(cavity%nsph, dirs)
      call tangent_of(cavity, dirs, tan_out, error)
      if (allocated(error)) return

      if (maxval(abs(tan_out(:, :, CH_WBRANCH))) /= 0.0_wp) then
         call test_failed(error, "d_wbranch is nonzero on an unbranched grid (max "// &
                          to_string(maxval(abs(tan_out(:, :, CH_WBRANCH))))//")")
         return
      end if

      ! The other three channels must still be alive, or the assertion above
      ! would be satisfied by a routine that returns nothing at all
      if (maxval(abs(tan_out(:, :, CH_WLEB))) <= VACUITY_THR(CH_WLEB)) then
         call test_failed(error, "the unbranched fixture carries no weight tangent")
         return
      end if

   end subroutine test_single_branch

   !> Mis-shaped arguments must be refused
   !>
   !> @param[out] error Error handle
   subroutine test_shape_guards(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol
      !> Four distinct buffers: the four outputs are `intent(out)` dummies, so
      !> passing one array to several of them would alias
      real(wp), allocatable :: dirs(:, :, :), o1(:, :), o2(:, :), o3(:, :), o4(:, :)
      integer :: ngrid

      call fixture_geometry(FIX_PLAIN, mol)
      call build_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error)
      if (allocated(error)) return

      ngrid = cavity%ngrid
      allocate (o1(ngrid, NDIR), o2(ngrid, NDIR), o3(ngrid, NDIR), o4(ngrid, NDIR))

      ! Wrong number of spheres in `dirs`
      allocate (dirs(ndim, cavity%nsph + 1, NDIR), source=0.1_wp)
      call cavity%get_surface_tangent(dirs, o1, o2, o3, o4, cav_error)
      if (.not. allocated(cav_error)) then
         call test_failed(error, "a mis-shaped direction array was accepted")
         return
      end if
      deallocate (cav_error, dirs)

      ! Output with the wrong direction count
      allocate (dirs(ndim, cavity%nsph, NDIR), source=0.1_wp)
      deallocate (o4)
      allocate (o4(ngrid, NDIR - 1))
      call cavity%get_surface_tangent(dirs, o1, o2, o3, o4, cav_error)
      if (.not. allocated(cav_error)) then
         call test_failed(error, "a mis-shaped output array was accepted")
         return
      end if
      deallocate (cav_error)

   end subroutine test_shape_guards

   !* ================================================================================= *!
   !*                             Preconditions of the tests                            *!
   !* ================================================================================= *!

   !> Assert that the fixture branches exactly as much as it is meant to
   !>
   !> @param[in]  cavity   Cavity to inspect
   !> @param[in]  fix_kind Geometry of the fixture
   !> @param[in]  label    Human-readable geometry description
   !> @param[out] error    Error handle
   subroutine assert_branching(cavity, fix_kind, label, error)
      !> Cavity to inspect
      type(cavity_type_drop), intent(in) :: cavity
      !> Geometry of the fixture
      integer, intent(in) :: fix_kind
      !> Geometry description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      logical :: branched

      if (.not. allocated(cavity%branch_count)) then
         call test_failed(error, "no branch bookkeeping at "//label)
         return
      end if
      branched = any(cavity%branch_count(1:cavity%ngrid) > 1)

      select case (fix_kind)
      case (FIX_PLAIN)
         if (branched) then
            call test_failed(error, "the single-branch fixture branched at "//label// &
                             " (max branch_count "// &
                             to_string(maxval(cavity%branch_count(1:cavity%ngrid)))//")")
            return
         end if
         if (any(cavity%wbranch(1:cavity%ngrid) /= 1.0_wp)) then
            call test_failed(error, "wbranch is not exactly one on the single-branch"// &
                             " fixture at "//label)
            return
         end if
      case default
         if (.not. branched) then
            call test_failed(error, "the branching fixture did not branch at "//label// &
                             "; d_wbranch would then be structurally zero and the"// &
                             " branch term of d_wleb untested")
            return
         end if
      end select

   end subroutine assert_branching

   !> Assert that two geometries carry the very same grid points
   !>
   !> The reference is compared slot by slot, so a point that appears, vanishes
   !> or changes branch identity between two stencil geometries would silently
   !> put a step into the differenced primal. `branch_count` and `anchor_id` are
   !> compared as well as `numbering` and `owner`, because a group that
   !> re-enumerates without changing the point set would corrupt the branch
   !> stage alone. This assertion is what bounds the step from above.
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
                             ", slot "//to_string(igrid)//" ("// &
                             to_string(ref%numbering(igrid))//" -> "// &
                             to_string(cav%numbering(igrid))//")")
            return
         end if
         if (cav%owner(igrid) /= ref%owner(igrid)) then
            call test_failed(error, "the owner sphere changed at "//label//", slot "// &
                             to_string(igrid)//" ("//to_string(ref%owner(igrid))//" -> "// &
                             to_string(cav%owner(igrid))//")")
            return
         end if
         if (cav%branch_count(igrid) /= ref%branch_count(igrid)) then
            call test_failed(error, "the branch count changed at "//label//", slot "// &
                             to_string(igrid)//" ("//to_string(ref%branch_count(igrid))// &
                             " -> "//to_string(cav%branch_count(igrid))//")")
            return
         end if
         if (cav%anchor_id(igrid) /= ref%anchor_id(igrid)) then
            call test_failed(error, "the anchor group changed at "//label//", slot "// &
                             to_string(igrid)//" ("//to_string(ref%anchor_id(igrid))// &
                             " -> "//to_string(cav%anchor_id(igrid))//")")
            return
         end if
      end do

   end subroutine assert_grid_match

   !> Anti-vacuity floor on every channel the fixture is supposed to drive
   !>
   !> Without it a fixture that quietly stopped exercising a term -- a switching
   !> function saturated, a group that stopped branching -- would keep passing
   !> every tolerance in the file.
   !>
   !> @param[in]  tan_out  Analytic tangent `(ngrid, NDIR, NCHAN)`
   !> @param[in]  fix_kind Geometry of the fixture
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine assert_live(tan_out, fix_kind, label, error)
      !> Analytic tangent
      real(wp), intent(in) :: tan_out(:, :, :)
      !> Geometry of the fixture
      integer, intent(in) :: fix_kind
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer :: ichannel, idir
      real(wp) :: peak

      do ichannel = 1, NCHAN
         if (ichannel == CH_WBRANCH .and. fix_kind == FIX_PLAIN) cycle
         do idir = 1, NDIR
            peak = maxval(abs(tan_out(:, idir, ichannel)))
            if (peak <= VACUITY_THR(ichannel)) then
               call test_failed(error, "channel "//trim(channel_name(ichannel))// &
                                " is vacuous for "//label//", direction "// &
                                to_string(idir)//" (max "//to_string(peak)//")")
               return
            end if
         end do
      end do

   end subroutine assert_live

   !* ================================================================================= *!
   !*                                     Fixture                                       *!
   !* ================================================================================= *!

   !> Run the tangent for a built cavity
   !>
   !> @param[in]  cavity  Cavity to differentiate
   !> @param[in]  dirs    Nuclear directions `(3, nsph, NDIR)`
   !> @param[out] tan_out Analytic tangent `(ngrid, NDIR, NCHAN)`
   !> @param[out] error   Error handle
   subroutine tangent_of(cavity, dirs, tan_out, error)
      !> Cavity to differentiate
      type(cavity_type_drop), intent(in) :: cavity
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Analytic tangent
      real(wp), allocatable, intent(out) :: tan_out(:, :, :)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: cav_error

      allocate (tan_out(cavity%ngrid, NDIR, NCHAN), source=0.0_wp)
      call cavity%get_surface_tangent(dirs, tan_out(:, :, CH_A), tan_out(:, :, CH_WLEB), &
                                      tan_out(:, :, CH_XI), tan_out(:, :, CH_WBRANCH), &
                                      cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "surface tangent failed: "//cav_error%message)
         return
      end if

   end subroutine tangent_of

   !> Name of an output channel, for diagnostics
   !>
   !> Fixed length on purpose. A deferred-length allocatable result returned
   !> into test-drive's `!$omp parallel do` over the suite's tests picks up a
   !> thread-shared static length temporary on gfortran and corrupts the
   !> message; the callers `trim` this instead.
   !>
   !> @param[in] ichannel Channel index
   !> @return             Channel name
   pure function channel_name(ichannel) result(name)
      !> Channel index
      integer, intent(in) :: ichannel
      !> Channel name
      character(len=10) :: name

      select case (ichannel)
      case (CH_A)
         name = "d_a"
      case (CH_WLEB)
         name = "d_wleb"
      case (CH_XI)
         name = "d_xi0"
      case default
         name = "d_wbranch"
      end select

   end function channel_name

   !> The two nuclear directions pushed through in one call
   !>
   !> Direction 1 moves a single atom along a single axis -- the sparsest
   !> column, and the one a wrong influence set would zero out. Direction 2
   !> moves every atom along its own vector and is not a rigid translation, so
   !> nothing about it cancels.
   !>
   !> @param[in]  nsph Number of spheres
   !> @param[out] dirs Nuclear directions `(3, nsph, NDIR)`
   subroutine build_directions(nsph, dirs)
      !> Number of spheres
      integer, intent(in) :: nsph
      !> Nuclear directions
      real(wp), allocatable, intent(out) :: dirs(:, :, :)

      integer :: iatom, iaxis

      if (allocated(dirs)) deallocate (dirs)
      allocate (dirs(ndim, nsph, NDIR), source=0.0_wp)

      dirs(3, 2, 1) = 1.0_wp

      do iatom = 1, nsph
         do iaxis = 1, ndim
            dirs(iaxis, iatom, 2) = sin(1.1_wp*real(iatom, wp) + 0.6_wp*real(iaxis, wp))
         end do
      end do

   end subroutine build_directions

   !> Fixture geometry
   !>
   !> `FIX_PLAIN` is asymmetric on purpose: a symmetric geometry drives the
   !> multistart projection into sibling branches, which is the other fixture's
   !> job. `FIX_CROSS` is the shared five-carbon cross, whose concave seams give
   !> several minima per anchor.
   !>
   !> @param[in]  fix_kind Geometry selector
   !> @param[out] mol      Structure
   subroutine fixture_geometry(fix_kind, mol)
      !> Geometry selector
      integer, intent(in) :: fix_kind
      !> Structure
      type(structure_type), intent(out) :: mol

      select case (fix_kind)
      case (FIX_PLAIN)
         call new(mol, [8, 6, 1], reshape([ &
                                          0.00_wp, 0.00_wp, 0.00_wp, &
                                          0.00_wp, 0.00_wp, 4.60_wp, &
                                          2.60_wp, 0.40_wp, -1.10_wp], [3, 3]))
      case default
         call get_test_cross(mol)
      end select

   end subroutine fixture_geometry

   !> Build the DROP cavity for a fixture and a level-set model
   !>
   !> @param[out]   cavity   Constructed cavity
   !> @param[inout] ctx      Run context borrowed by the cavity; must outlive it
   !> @param[in]    mol      Structure to build on
   !> @param[in]    fix_kind Geometry of the fixture
   !> @param[in]    lsf_kind Level-set model
   !> @param[out]   error    Error handle
   subroutine build_cavity(cavity, ctx, mol, fix_kind, lsf_kind, error)
      !> Constructed cavity
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(inout) :: ctx
      !> Structure to build on
      type(structure_type), intent(in) :: mol
      !> Geometry of the fixture
      integer, intent(in) :: fix_kind
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: radii(:)
      type(mctc_error), allocatable :: cav_error
      !> Per-fixture settings. Named apart from the module parameters they are
      !> assigned from: Fortran is case insensitive, so a local `blend_k` would
      !> shadow `BLEND_K` and turn the assignment into a silent self-assignment.
      real(wp) :: blend_k_loc, branch_s_loc
      integer :: proj_level_loc

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      ! The softmax scale is only raised on the branching fixture. It also sets
      ! the admissible branch radius (`branch_dphi_max` in `parameters.f90`),
      ! and at `s = 0.5` with this `wleb_cut` that radius is tens of Bohr --
      ! which is what keeps the cross's far siblings alive, and what empties
      ! the grid outright on the compact one.
      if (fix_kind == FIX_CROSS) then
         blend_k_loc = CROSS_BLEND_K
         proj_level_loc = CROSS_PROJ_LEVEL
         branch_s_loc = CROSS_BRANCH_S
      else
         blend_k_loc = BLEND_K
         proj_level_loc = PROJ_LEVEL
         branch_s_loc = PLAIN_BRANCH_S
      end if

      allocate (cavity)
      call new_context(ctx, verbosity=0)
      select case (lsf_kind)
      case (LSF_SVDW)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=blend_k_loc, blend_3b=BLEND_3B)
            call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                                 tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                                 proj_level=proj_level_loc, wleb_prune_level=WLEB_PRUNE, &
                                 branch_weight_s=branch_s_loc, &
                                 radius_model=default_cpcm_radii(), &
                                 lsf_model=svdw_template, error=cav_error)
         end block
      case default
         block
            type(moist_cavity_drop_lsf_cfc_type) :: cfc_template
            call cfc_template%new()
            call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                                 tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                                 proj_level=proj_level_loc, wleb_prune_level=WLEB_PRUNE, &
                                 branch_weight_s=branch_s_loc, &
                                 radius_model=default_cpcm_radii(), &
                                 lsf_model=cfc_template, error=cav_error)
         end block
      end select
      if (allocated(cav_error)) then
         call test_failed(error, "failed to initialize cavity: "//cav_error%message)
         return
      end if

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

end module test_cavity_drop_tangent_forward
