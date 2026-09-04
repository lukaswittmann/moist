!> Test suite for the surface-weight folding tangent
!>
!> [[prepare_surface_weights_tangent]] differentiates the folding of
!> `derivatives/weights.f90` along one direction with the raw host adjoints
!> held fixed. Finite differences of the folding itself are the ground truth.
!>
!> **How the ground truth is obtained, and why.** `prepare_surface_weights` is
!> a plain module procedure of `moist_cavity_drop` -- not type bound, and not in
!> that module's public list -- so a test cannot call it without building a
!> complete DROP cavity, which is far more setup than a pure function of a
!> handful of per-point scalars deserves. The suite therefore splits the
!> reference in two:
!>
!>   * the four fold statements are transcribed into [[reference_fold]] below,
!>     against `weights.f90:100-118`. Four lines, and a transcription slip
!>     shows up immediately because the tangent under test was derived from the
!>     shipped source rather than from the transcription.
!>   * the branch reverse pass -- the only cross-point coupling in the scheme,
!>     and the part where a slip would be hard to see -- is the **shipped**
!>     [[compute_branch_phi_adj]], which is public. So the subtle half is
!>     finite-differenced against production code.
!>
!> The fixture displaces `a`, `wleb`, `xi0` and `wbranch` along an arbitrary
!> linear path, deliberately ignoring the identities that tie them together for
!> a physical direction (`dwleb = -2 wleb dxi0/xi0`, `da = R^2(wleb df + f
!> dwleb)`). The primitive does not impose those identities -- for the reason
!> [[drop_seed_input_tangent_type]] gives -- so an arbitrary path is a stronger
!> test of it than a physical one, which would leave two of the four input
!> channels linearly dependent on the third.
!>
!> Coverage the fixture is built for, point by point:
!>
!> | igrid | group        | `w_a` | `w_w` | what it pins                       |
!> |-------|--------------|-------|-------|------------------------------------|
!> | 1     | A, 3 branches| live  | live  | both folds on                      |
!> | 2     | A            | zero  | live  | area fold gated off                |
!> | 3     | A            | zero  | zero  | `w_xi == 0`, branch gate shut      |
!> | 4     | B, 2 branches| live  | zero  | weight fold gated off              |
!> | 5     | B            | live  | live  | both folds on inside a 2-group     |
!> | 6     | single       | live  | live  | folds on, no branch coupling       |
!> | 7     | single       | live  | live  | group walk resumes after singles   |
!> | 8, 9  | C, 2 branches| live  | live  | a group after the singles          |
!>
!> Point 3 is the load-bearing one: its own branch adjoint is gated off in both
!> the primal and the tangent, so everything its `branch_phi_adj` tangent
!> carries comes through the group reduction. `branch_coupling_is_group_local`
!> uses it to show the coupling is live and does not cross group boundaries.
module test_cavity_drop_weights_tangent
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, test_failed
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_drop_derivatives_kernel, only: drop_surface_weights_type, &
      & compute_branch_phi_adj, seed_weight_tol
   use moist_cavity_drop_derivatives_weights_tangent, only: prepare_surface_weights_tangent, &
      & branch_phi_adj_tangent

   implicit none(type, external)
   private

   public :: collect_cavity_drop_weights_tangent

   !> Fixture extent; nine points over three spheres
   integer, parameter :: ngrid = 9, nsph = 3

   !> Softmax temperature of the fixture, comfortably above `seed_weight_tol`
   real(wp), parameter :: sigma_phi = 0.37_wp

   !> Anchor groups: one triple, one pair, two singles, one pair. The two
   !> singles sit between groups so the walk has to resume correctly.
   integer, parameter :: anchor_id(ngrid) = [1, 1, 1, 2, 2, 3, 4, 5, 5]
   integer, parameter :: branch_count(ngrid) = [3, 3, 3, 2, 2, 1, 1, 2, 2]
   integer, parameter :: owner(ngrid) = [1, 1, 1, 2, 2, 3, 1, 2, 3]

   !> Sphere radii; distinct, so a wrong `owner` lookup cannot pass
   real(wp), parameter :: radii(nsph) = [1.9_wp, 2.1_wp, 1.7_wp]

   !> Which components must be nonzero, and which must be structurally zero
   !>
   !> These double as the anti-vacuity floor and as the guard-branch assertion:
   !> a `.false.` entry is checked to be exactly zero, a `.true.` one to clear
   !> `live_floor` before its value is compared with the difference quotient.
   !> Without them a fixture that silently stopped exercising a term would keep
   !> passing every tolerance in the file.
   logical, parameter :: live_w_xi(ngrid) = &
                         [.true., .true., .false., .true., .true., .true., .true., .true., .true.]
   logical, parameter :: live_w_f(ngrid) = &
                         [.true., .false., .false., .true., .true., .true., .true., .true., .true.]
   logical, parameter :: live_bpa(ngrid) = &
                         [.true., .true., .true., .true., .true., .false., .false., .true., .true.]

   !> Nothing is live; used for the channel that is a verbatim copy
   logical, parameter :: none_live(ngrid) = .false.

   !> Magnitude every live component must clear
   real(wp), parameter :: live_floor = 1.0e-3_wp

   !> Tolerance of the finite-difference comparisons, scaled by the component
   !> magnitude. Observed worst deviation over all three channels is 5.8e-11 at
   !> `h = 1e-5` and 1.2e-10 at `h = 1e-6`, against components of size 0.28 to
   !> 2.2. No monotonicity in `h` is claimed or asserted: the smaller step is
   !> already roundoff dominated.
   !>
   !> Tightened to the project-wide `1e-10` target on 2026-09-03. This sits
   !> right on the measured floor. Confirmed: both `h = 1e-5` cases pass, and
   !> both `h = 1e-6` cases fail marginally -- `dw_f` at point 7 by `1.21e-10`
   !> and `dbranch_phi_adj` at point 1 by `1.16e-10`. That the deviation *rises*
   !> as the step shrinks is the signature of roundoff, so `1e-10` is within a
   !> small factor of the best this comparison can do; `1.5e-10` would hold.
   real(wp), parameter :: fd_tol = 1.0e-10_wp

contains

   !> Register the surface-weight tangent tests
   subroutine collect_cavity_drop_weights_tangent(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fold_switching_fd_h1em5", test_switching_h1em5), &
                  new_unittest("fold_switching_fd_h1em6", test_switching_h1em6), &
                  new_unittest("no_fold_switching_fd_h1em5", test_no_switching_h1em5), &
                  new_unittest("no_fold_switching_fd_h1em6", test_no_switching_h1em6), &
                  new_unittest("branch_coupling_is_group_local", test_branch_coupling), &
                  new_unittest("branch_pass_early_exits", test_branch_early_exits), &
                  new_unittest("shape_mismatch_reported", test_shape_mismatch), &
                  new_unittest("uninitialized_accumulator_reported", test_uninitialized) &
                  ]
   end subroutine collect_cavity_drop_weights_tangent

   !* ================================================================================= *!
   !*                                      Fixture                                      *!
   !* ================================================================================= *!

   !> Build the accumulator, the cavity scalars and the direction tangents
   !>
   !> Every value is O(1) and no two are equal, so a transposed or off-by-one
   !> index shows up as a gross error rather than as a rounding difference.
   !> `w_a` is exactly zero at points 2 and 3 and `w_w` exactly zero at points 3
   !> and 4, which is what shuts the folds; point 3 additionally has
   !> `w_xi == 0`, so its folded width channel is identically zero along the
   !> whole path and the branch gate stays shut at every displacement.
   !>
   !> @param[out] acc      Raw surface adjoints
   !> @param[out] a        Surface areas
   !> @param[out] wleb     Lebedev weights
   !> @param[out] xi0      Gaussian widths
   !> @param[out] wbranch  Branch weights
   !> @param[out] da       Direction tangent of the areas
   !> @param[out] dwleb    Direction tangent of the Lebedev weights
   !> @param[out] dxi0     Direction tangent of the widths
   !> @param[out] dwbranch Direction tangent of the branch weights
   subroutine weights_fixture(acc, a, wleb, xi0, wbranch, da, dwleb, dxi0, dwbranch)
      !> Raw surface adjoints
      type(cavity_surface_adjoint_type), intent(out) :: acc
      !> Cavity grid scalars
      real(wp), intent(out) :: a(ngrid), wleb(ngrid), xi0(ngrid), wbranch(ngrid)
      !> Direction tangents of the same
      real(wp), intent(out) :: da(ngrid), dwleb(ngrid), dxi0(ngrid), dwbranch(ngrid)

      call acc%init(ngrid)

      acc%w_xi = [0.71_wp, -0.43_wp, 0.00_wp, 0.58_wp, -0.92_wp, &
                  0.34_wp, 0.66_wp, -0.27_wp, 0.81_wp]
      acc%w_f = [0.19_wp, 0.28_wp, -0.37_wp, 0.45_wp, 0.13_wp, &
                 -0.22_wp, 0.31_wp, 0.52_wp, -0.16_wp]
      acc%w_a = [0.63_wp, 0.00_wp, 0.00_wp, -0.74_wp, 0.48_wp, &
                 0.29_wp, 0.55_wp, -0.36_wp, 0.42_wp]
      acc%w_w = [-0.51_wp, 0.67_wp, 0.00_wp, 0.00_wp, 0.39_wp, &
                 -0.44_wp, 0.23_wp, 0.58_wp, -0.31_wp]

      a = [0.83_wp, 1.14_wp, 0.71_wp, 1.32_wp, 0.96_wp, &
           1.05_wp, 0.77_wp, 1.21_wp, 0.89_wp]
      wleb = [0.52_wp, 0.61_wp, 0.44_wp, 0.73_wp, 0.58_wp, &
              0.66_wp, 0.49_wp, 0.71_wp, 0.55_wp]
      xi0 = [1.31_wp, 0.94_wp, 1.62_wp, 1.08_wp, 1.47_wp, &
             0.87_wp, 1.23_wp, 1.05_wp, 1.38_wp]
      wbranch = [0.41_wp, 0.35_wp, 0.24_wp, 0.62_wp, 0.38_wp, &
                 1.00_wp, 1.00_wp, 0.57_wp, 0.43_wp]

      da = [0.37_wp, -0.62_wp, 0.45_wp, 0.28_wp, -0.51_wp, &
            0.74_wp, -0.33_wp, 0.59_wp, -0.48_wp]
      dwleb = [-0.24_wp, 0.53_wp, -0.31_wp, 0.66_wp, 0.42_wp, &
               -0.57_wp, 0.39_wp, -0.68_wp, 0.26_wp]
      dxi0 = [0.46_wp, -0.35_wp, 0.58_wp, -0.29_wp, 0.63_wp, &
              0.41_wp, -0.52_wp, 0.34_wp, -0.61_wp]
      ! Nonzero at the two single-branch points on purpose: their
      ! `branch_phi_adj` tangent must stay exactly zero regardless.
      dwbranch = [0.32_wp, -0.27_wp, 0.45_wp, -0.38_wp, 0.21_wp, &
                  0.15_wp, -0.22_wp, 0.49_wp, -0.44_wp]

   end subroutine weights_fixture

   !> Transcription of the folding of [[prepare_surface_weights]]
   !>
   !> `weights.f90:100-118`, statement for statement. The one deviation is
   !> `r_own*r_own` for the primal's `radii(owner(igrid))**2`, which the project
   !> forbids for small integer powers; the two are the same number.
   !>
   !> @param[in]  acc            Raw surface adjoints
   !> @param[in]  fold_switching Whether the area channel folds into `w_f`
   !> @param[in]  a              Surface areas
   !> @param[in]  wleb           Lebedev weights
   !> @param[in]  xi0            Gaussian widths
   !> @param[out] w_xi           Folded width channel
   !> @param[out] w_f            Folded switching channel
   subroutine reference_fold(acc, fold_switching, a, wleb, xi0, w_xi, w_f)
      !> Raw surface adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Whether the area channel folds into the switching channel
      logical, intent(in) :: fold_switching
      !> Cavity grid scalars
      real(wp), intent(in) :: a(ngrid), wleb(ngrid), xi0(ngrid)
      !> Folded channels
      real(wp), intent(out) :: w_xi(ngrid), w_f(ngrid)

      !> Grid index
      integer :: igrid
      !> Owner radius
      real(wp) :: r_own

      w_xi = acc%w_xi
      w_f = acc%w_f

      do igrid = 1, ngrid
         if (abs(acc%w_a(igrid)) > seed_weight_tol) then
            w_xi(igrid) = w_xi(igrid) - 2.0_wp*a(igrid)*acc%w_a(igrid)/xi0(igrid)
         end if
         if (abs(acc%w_w(igrid)) > seed_weight_tol) then
            w_xi(igrid) = w_xi(igrid) - 2.0_wp*wleb(igrid)*acc%w_w(igrid)/xi0(igrid)
         end if
      end do

      if (fold_switching) then
         do igrid = 1, ngrid
            if (abs(acc%w_a(igrid)) > seed_weight_tol) then
               r_own = radii(owner(igrid))
               w_f(igrid) = w_f(igrid) + acc%w_a(igrid)*r_own*r_own*wleb(igrid)
            end if
         end do
      end if

   end subroutine reference_fold

   !> Full primal: the transcribed fold plus the shipped branch reverse pass
   !>
   !> @param[in]  acc            Raw surface adjoints
   !> @param[in]  fold_switching Whether the area channel folds into `w_f`
   !> @param[in]  a              Surface areas
   !> @param[in]  wleb           Lebedev weights
   !> @param[in]  xi0            Gaussian widths
   !> @param[in]  wbranch        Branch weights
   !> @param[out] w_xi           Folded width channel
   !> @param[out] w_f            Folded switching channel
   !> @param[out] bpa            Branch-objective adjoint
   subroutine reference_prepare(acc, fold_switching, a, wleb, xi0, wbranch, w_xi, w_f, bpa)
      !> Raw surface adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Whether the area channel folds into the switching channel
      logical, intent(in) :: fold_switching
      !> Cavity grid scalars
      real(wp), intent(in) :: a(ngrid), wleb(ngrid), xi0(ngrid), wbranch(ngrid)
      !> Folded channels and the branch adjoint
      real(wp), intent(out) :: w_xi(ngrid), w_f(ngrid), bpa(ngrid)

      call reference_fold(acc, fold_switching, a, wleb, xi0, w_xi, w_f)
      call compute_branch_phi_adj(branch_count, anchor_id, wbranch, wleb, xi0, &
                                  sigma_phi, w_xi, bpa)

   end subroutine reference_prepare

   !* ================================================================================= *!
   !*                             Finite-difference driver                              *!
   !* ================================================================================= *!

   !> Compare the tangent with a central difference of the folding
   !>
   !> Displaces all four grid scalars along the fixture direction, which is not
   !> a physical one; see the module header for why that is deliberate.
   !>
   !> @param[out] error          Test error
   !> @param[in]  fold_switching Whether the area channel folds into `w_f`
   !> @param[in]  step           Finite-difference step
   subroutine run_fd_case(error, fold_switching, step)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Whether the area channel folds into the switching channel
      logical, intent(in) :: fold_switching
      !> Finite-difference step
      real(wp), intent(in) :: step

      type(cavity_surface_adjoint_type) :: acc
      type(drop_surface_weights_type) :: eff
      type(mctc_error), allocatable :: merr
      real(wp) :: a(ngrid), wleb(ngrid), xi0(ngrid), wbranch(ngrid)
      real(wp) :: da(ngrid), dwleb(ngrid), dxi0(ngrid), dwbranch(ngrid)
      real(wp) :: dw_xi(ngrid), dw_f(ngrid), dbpa(ngrid)
      real(wp) :: fd_w_xi(ngrid), fd_w_f(ngrid), fd_bpa(ngrid)
      real(wp) :: wp_xi(ngrid), wp_f(ngrid), bp(ngrid)
      real(wp) :: wm_xi(ngrid), wm_f(ngrid), bm(ngrid)
      character(len=8) :: tag

      call weights_fixture(acc, a, wleb, xi0, wbranch, da, dwleb, dxi0, dwbranch)

      write (tag, "(es8.1)") step

      allocate (eff%w_xi(ngrid), eff%w_f(ngrid))
      call reference_fold(acc, fold_switching, a, wleb, xi0, eff%w_xi, eff%w_f)

      call prepare_surface_weights_tangent(acc, eff, fold_switching, &
                                           a, wleb, xi0, wbranch, radii, owner, &
                                           branch_count, anchor_id, sigma_phi, &
                                           da, dwleb, dxi0, dwbranch, &
                                           dw_xi, dw_f, dbpa, merr)
      if (allocated(merr)) then
         call test_failed(error, "tangent rejected the fixture: "//merr%message)
         return
      end if

      call reference_prepare(acc, fold_switching, a + step*da, wleb + step*dwleb, &
                             xi0 + step*dxi0, wbranch + step*dwbranch, wp_xi, wp_f, bp)
      call reference_prepare(acc, fold_switching, a - step*da, wleb - step*dwleb, &
                             xi0 - step*dxi0, wbranch - step*dwbranch, wm_xi, wm_f, bm)

      fd_w_xi = (wp_xi - wm_xi)/(2.0_wp*step)
      fd_w_f = (wp_f - wm_f)/(2.0_wp*step)
      fd_bpa = (bp - bm)/(2.0_wp*step)

      call check_channel(error, dw_xi, fd_w_xi, live_w_xi, "dw_xi at h="//trim(tag))
      if (allocated(error)) return

      ! With the switching fold off, `w_f` is a verbatim copy of the fixed raw
      ! adjoint, so every point is a structural zero rather than a small number.
      if (fold_switching) then
         call check_channel(error, dw_f, fd_w_f, live_w_f, "dw_f at h="//trim(tag))
      else
         call check_channel(error, dw_f, fd_w_f, none_live, "dw_f at h="//trim(tag))
      end if
      if (allocated(error)) return

      call check_channel(error, dbpa, fd_bpa, live_bpa, "dbranch_phi_adj at h="//trim(tag))

   end subroutine run_fd_case

   !> Assert one channel against its difference quotient, with the vacuity floor
   !>
   !> @param[inout] error   Test error
   !> @param[in]    tangent Analytic tangent
   !> @param[in]    fd      Central difference of the same quantity
   !> @param[in]    live    Which components must carry a nonzero value
   !> @param[in]    label   Channel name for the diagnostics
   subroutine check_channel(error, tangent, fd, live, label)
      !> Test error
      type(error_type), allocatable, intent(inout) :: error
      !> Analytic tangent and its difference quotient
      real(wp), intent(in) :: tangent(ngrid), fd(ngrid)
      !> Components that must be nonzero
      logical, intent(in) :: live(ngrid)
      !> Channel name
      character(len=*), intent(in) :: label

      !> Grid index
      integer :: igrid
      !> Allowed deviation at this component
      real(wp) :: thr

      do igrid = 1, ngrid
         if (live(igrid)) then
            call check(error, abs(tangent(igrid)) > live_floor, &
                       label//": point "//to_string(igrid)//" is "// &
                       to_string(tangent(igrid))//", too small to test; the fixture"// &
                       " no longer exercises this term")
            if (allocated(error)) return
         else
            call check(error, tangent(igrid) == 0.0_wp, &
                       label//": point "//to_string(igrid)//" must be a structural"// &
                       " zero but is "//to_string(tangent(igrid)))
            if (allocated(error)) return
            call check(error, fd(igrid) == 0.0_wp, &
                       label//": point "//to_string(igrid)//" is a structural zero"// &
                       " in the tangent but the primal moved by "//to_string(fd(igrid)))
            if (allocated(error)) return
         end if

         thr = fd_tol*max(1.0_wp, abs(tangent(igrid)))
         call check(error, tangent(igrid), fd(igrid), thr=thr, &
                    message=label//": point "//to_string(igrid)//" is "// &
                    to_string(tangent(igrid))//", finite difference "//to_string(fd(igrid)))
         if (allocated(error)) return
      end do

   end subroutine check_channel

   !* ================================================================================= *!
   !*                                       Tests                                       *!
   !* ================================================================================= *!

   !> Nuclear path: the area channel also folds into `w_f`
   !>
   !> @param[out] error Test error
   subroutine test_switching_h1em5(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      call run_fd_case(error, .true., 1.0e-5_wp)
   end subroutine test_switching_h1em5

   !> Nuclear path at the smaller step
   !>
   !> @param[out] error Test error
   subroutine test_switching_h1em6(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      call run_fd_case(error, .true., 1.0e-6_wp)
   end subroutine test_switching_h1em6

   !> Electronic path: `f` is fixed, so the area channel skips `w_f`
   !>
   !> @param[out] error Test error
   subroutine test_no_switching_h1em5(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      call run_fd_case(error, .false., 1.0e-5_wp)
   end subroutine test_no_switching_h1em5

   !> Electronic path at the smaller step
   !>
   !> @param[out] error Test error
   subroutine test_no_switching_h1em6(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      call run_fd_case(error, .false., 1.0e-6_wp)
   end subroutine test_no_switching_h1em6

   !> The branch tangent couples within a group and only within a group
   !>
   !> Point 3's own branch adjoint is gated off, so its tangent is pure group
   !> reduction. Zeroing the tangents of points 1 and 2 -- its two partners --
   !> must therefore move it, and must leave the unrelated group at 8, 9 and the
   !> single-branch points untouched. A finite-difference test cannot separate
   !> these two failures: a driver that reduced globally instead of per group
   !> would still match its own reference.
   !>
   !> @param[out] error Test error
   subroutine test_branch_coupling(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_surface_adjoint_type) :: acc
      type(drop_surface_weights_type) :: eff
      type(mctc_error), allocatable :: merr
      real(wp) :: a(ngrid), wleb(ngrid), xi0(ngrid), wbranch(ngrid)
      real(wp) :: da(ngrid), dwleb(ngrid), dxi0(ngrid), dwbranch(ngrid)
      real(wp) :: dw_xi(ngrid), dw_f(ngrid), dbpa(ngrid), dbpa_cut(ngrid)
      integer :: igrid

      call weights_fixture(acc, a, wleb, xi0, wbranch, da, dwleb, dxi0, dwbranch)
      allocate (eff%w_xi(ngrid), eff%w_f(ngrid))
      call reference_fold(acc, .true., a, wleb, xi0, eff%w_xi, eff%w_f)

      call prepare_surface_weights_tangent(acc, eff, .true., &
                                           a, wleb, xi0, wbranch, radii, owner, &
                                           branch_count, anchor_id, sigma_phi, &
                                           da, dwleb, dxi0, dwbranch, &
                                           dw_xi, dw_f, dbpa, merr)
      if (allocated(merr)) then
         call test_failed(error, "tangent rejected the fixture: "//merr%message)
         return
      end if

      ! Silence the two partners of point 3, nothing else
      da(1:2) = 0.0_wp
      dwleb(1:2) = 0.0_wp
      dxi0(1:2) = 0.0_wp
      dwbranch(1:2) = 0.0_wp

      call prepare_surface_weights_tangent(acc, eff, .true., &
                                           a, wleb, xi0, wbranch, radii, owner, &
                                           branch_count, anchor_id, sigma_phi, &
                                           da, dwleb, dxi0, dwbranch, &
                                           dw_xi, dw_f, dbpa_cut, merr)
      if (allocated(merr)) then
         call test_failed(error, "tangent rejected the cut fixture: "//merr%message)
         return
      end if

      call check(error, abs(dbpa(3) - dbpa_cut(3)) > live_floor, &
                 "Point 3 has no branch adjoint of its own, so its tangent "// &
                 to_string(dbpa(3))//" must come from the group reduction, but "// &
                 "silencing points 1 and 2 left it at "//to_string(dbpa_cut(3)))
      if (allocated(error)) return

      do igrid = 6, ngrid
         call check(error, dbpa(igrid) == dbpa_cut(igrid), &
                    "Point "//to_string(igrid)//" is outside the perturbed anchor "// &
                    "group but moved from "//to_string(dbpa(igrid))//" to "// &
                    to_string(dbpa_cut(igrid))//"; the reduction is not group local")
         if (allocated(error)) return
      end do

   end subroutine test_branch_coupling

   !> The branch tangent takes the same early exits as the primal
   !>
   !> No group wider than one branch, and a vanishing softmax temperature, both
   !> leave the whole channel at zero. The fixture tangents are left live, so a
   !> missing early exit produces nonzero output rather than a silent pass.
   !>
   !> @param[out] error Test error
   subroutine test_branch_early_exits(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_surface_adjoint_type) :: acc
      real(wp) :: a(ngrid), wleb(ngrid), xi0(ngrid), wbranch(ngrid)
      real(wp) :: da(ngrid), dwleb(ngrid), dxi0(ngrid), dwbranch(ngrid)
      real(wp) :: w_xi(ngrid), w_f(ngrid), dw_xi(ngrid), dbpa(ngrid)
      integer :: single(ngrid)

      call weights_fixture(acc, a, wleb, xi0, wbranch, da, dwleb, dxi0, dwbranch)
      call reference_fold(acc, .true., a, wleb, xi0, w_xi, w_f)
      dw_xi = da + dwleb + dxi0

      single = 1
      call branch_phi_adj_tangent(single, anchor_id, wbranch, wleb, xi0, sigma_phi, &
                                  w_xi, dwbranch, dwleb, dxi0, dw_xi, dbpa)
      call check(error, all(dbpa == 0.0_wp), &
                 "Every point is single branched, but the tangent reached "// &
                 to_string(maxval(abs(dbpa))))
      if (allocated(error)) return

      call branch_phi_adj_tangent(branch_count, anchor_id, wbranch, wleb, xi0, 0.0_wp, &
                                  w_xi, dwbranch, dwleb, dxi0, dw_xi, dbpa)
      call check(error, all(dbpa == 0.0_wp), &
                 "A vanishing softmax temperature must skip the branch tangent, "// &
                 "but it reached "//to_string(maxval(abs(dbpa))))

   end subroutine test_branch_early_exits

   !> A grid argument of the wrong length is reported, not folded
   !>
   !> @param[out] error Test error
   subroutine test_shape_mismatch(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_surface_adjoint_type) :: acc
      type(drop_surface_weights_type) :: eff
      type(mctc_error), allocatable :: merr
      real(wp) :: a(ngrid), wleb(ngrid), xi0(ngrid), wbranch(ngrid)
      real(wp) :: da(ngrid), dwleb(ngrid), dxi0(ngrid), dwbranch(ngrid)
      real(wp) :: dw_xi(ngrid), dw_f(ngrid), dbpa(ngrid)

      call weights_fixture(acc, a, wleb, xi0, wbranch, da, dwleb, dxi0, dwbranch)
      allocate (eff%w_xi(ngrid), eff%w_f(ngrid))
      call reference_fold(acc, .true., a, wleb, xi0, eff%w_xi, eff%w_f)

      call prepare_surface_weights_tangent(acc, eff, .true., &
                                           a, wleb, xi0, wbranch, radii, owner, &
                                           branch_count, anchor_id, sigma_phi, &
                                           da, dwleb, dxi0(1:ngrid - 1), dwbranch, &
                                           dw_xi, dw_f, dbpa, merr)
      call check(error, allocated(merr), &
                 "A width tangent of "//to_string(ngrid - 1)//" points was accepted "// &
                 "for a grid of "//to_string(ngrid))

   end subroutine test_shape_mismatch

   !> An accumulator that was never initialized is reported
   !>
   !> @param[out] error Test error
   subroutine test_uninitialized(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      type(cavity_surface_adjoint_type) :: acc, bare
      type(drop_surface_weights_type) :: eff
      type(mctc_error), allocatable :: merr
      real(wp) :: a(ngrid), wleb(ngrid), xi0(ngrid), wbranch(ngrid)
      real(wp) :: da(ngrid), dwleb(ngrid), dxi0(ngrid), dwbranch(ngrid)
      real(wp) :: dw_xi(ngrid), dw_f(ngrid), dbpa(ngrid)

      call weights_fixture(acc, a, wleb, xi0, wbranch, da, dwleb, dxi0, dwbranch)
      allocate (eff%w_xi(ngrid), eff%w_f(ngrid))
      call reference_fold(acc, .true., a, wleb, xi0, eff%w_xi, eff%w_f)

      call prepare_surface_weights_tangent(bare, eff, .true., &
                                           a, wleb, xi0, wbranch, radii, owner, &
                                           branch_count, anchor_id, sigma_phi, &
                                           da, dwleb, dxi0, dwbranch, &
                                           dw_xi, dw_f, dbpa, merr)
      call check(error, allocated(merr), &
                 "An unallocated accumulator was accepted")

   end subroutine test_uninitialized

end module test_cavity_drop_weights_tangent
