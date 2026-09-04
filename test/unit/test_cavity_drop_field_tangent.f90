!> Tests for the DROP field-contraction tangent
!>
!> [[drop_field_tangent]] differentiates the jet-contracted nuclear-gradient
!> row `vjp_f1_rA(w0, w1, w2)` along a nuclear direction `v`, with the adjoint
!> weights carrying their own tangents `(dw0, dw1, dw2)` and the evaluation
!> point riding along on `dr`. That splits into two structurally different
!> halves, and this suite pins each of them by a different mechanism:
!>
!>   * **Term (A)**, the weight tangents. The row is *linear* in the weights,
!>     so this half is exactly `vjp_f1_rA(dw0, dw1, dw2)` -- an algebraic
!>     identity, checked at roundoff rather than at finite-difference accuracy.
!>   * **Term (B)**, the motion of the jet. `vjp_f1_rA` at *fixed* weights is a
!>     function of the nuclear positions and of the evaluation point, so this
!>     half is the ordinary derivative of that function along `(v, dr)` and is
!>     checked against a finite difference of it.
!>
!> Both are FD'd along the *combined* direction rather than atom by atom, so a
!> reference costs four LSF evaluations, not `3 * nat` of them.
!>
!> Five channel combinations run at every sampling point, so that a failure
!> names the term that broke rather than just "the tangent is wrong":
!>
!>   1. `weights` -- `w = 0`, everything else on. Term (B) vanishes identically
!>      and the answer must be exactly `vjp_f1_rA(dw0, dw1, dw2)`. Also proves
!>      that a nonzero `dr`/`v` cannot leak in through a vanishing weight.
!>   2. `nuclear` -- `dr = 0`, `dw = 0`. The three `hvp_*` channels alone.
!>   3. `point`   -- `v = 0`, `dw = 0`. The three point-motion channels alone,
!>      two of which are folded into the `vjp_f1_rA` call and one of which is
!>      not.
!>   4. `f4`      -- `v = 0`, `dw = 0`, `w0 = w1 = 0`. Everything except the
!>      unfoldable `w2_ab T3_abk dr_k` term is zero by construction, so this
!>      case *is* the `f4_rrr_rA` channel. It is the one term with no
!>      contracted accessor behind it, hence the one most easily left out.
!>   5. `full`    -- all channels at once, to catch a term that is right in
!>      isolation and mis-combined in the sum.
!>
!> Every case carries an anti-vacuity floor on the half it is meant to
!> exercise: a reference at machine zero would let a wrong implementation pass,
!> so it fails the test instead.
!>
!> Fixtures, radii, sampling points and the `STEP_SIZE = 1e-3` 4-point stencil
!> are the ones `test_cavity_drop_lsf.f90` uses. The stencil matters: the FD
!> here differentiates a row that already carries three derivatives, so its
!> reference reaches total order four -- the same regime as the `f4_rrr_rA`
!> checks there. The thresholds below are measured at that step rather than
!> inherited; see their comments.
module test_cavity_drop_field_tangent
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use test_helpers, only: get_test_structures, get_test_radii, get_test_points, fd4_scalar
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_derivatives_field_tangent, only: drop_field_tangent, &
                                                         drop_field_tangent_work_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

   implicit none(type, external)
   private

   public :: collect_cavity_drop_field_tangent

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Nuclear/point displacement of the finite-difference stencil. The step and
   !> the 4-point `fd4_scalar` stencil are [[test_cavity_drop_lsf]]'s, and the
   !> choice is not free: the reference here is a fourth total derivative, so
   !> only an `O(h^4)` stencil reaches it at all.
   !>
   !> `1e-3` is the measured optimum rather than an inherited constant. Halving
   !> the step from `2e-3` divides the worst deviation by 16.1 (SvdW) and 16.1
   !> (CFC) -- clean truncation -- while halving it again only buys a further
   !> factor 7.6 and 3.8, cancellation having taken over. See the threshold
   !> comments below for the values at this step.
   real(wp), parameter :: STEP_SIZE = 1.0e-3_wp

   !> Finite-difference thresholds of the SvdW dispatch.
   !>
   !> Tightened to the project-wide `1e-10 / 1e-10` target on 2026-09-03.
   !> The previously measured worst deviation over the whole fixture sweep is
   !> `1.1e-11`, an order below the new bound, so this pair is expected to hold;
   !> the old `2e-10 / 1e-9` was headroom rather than floor. It is also the pair
   !> `svdw_f4_rrr_ra` runs at, which is the expected coincidence -- both
   !> references are a fourth total derivative of the same level set at the same
   !> step.
   real(wp), parameter :: ABS_THR = 1.0e-10_wp
   real(wp), parameter :: REL_THR = 1.0e-10_wp

   !> Finite-difference thresholds of the CFC dispatch.
   !>
   !> CFC's steep exponents (a1 = -15, a2 = -9) make it vary on a length scale
   !> of R/15, so the same stencil carries more truncation error: the worst
   !> deviation here is `5.9e-11`, five times SvdW's. One decade of headroom
   !> above that is `1e-9`.
   !>
   !> That is a decade *tighter* than the `CFC_ABS_THR` the LSF suite uses for
   !> its own order-4 checks, and deliberately so -- this reference is a single
   !> directional difference of an already-contracted row, not a difference of
   !> the raw `f3_rrr` ladder, and it is measurably more accurate. Loosening to
   !> the LSF suite's pair would leave two decades of slack.
   !>
   !> Tightened to `1e-10 / 1e-10` on 2026-09-03 anyway, against the measured
   !> `5.9e-11`, and it holds: both CFC cases pass at the new bound with about a
   !> factor two of margin. The old `1e-9 / 1e-8` was headroom, not floor.
   real(wp), parameter :: CFC_ABS_THR = 1.0e-10_wp
   real(wp), parameter :: CFC_REL_THR = 1.0e-10_wp

   !> Threshold of the algebraic half. Nothing numerical stands between
   !> [[drop_field_tangent]] with zero weights and `vjp_f1_rA(dw0, dw1, dw2)`:
   !> the extra terms are multiplied by an exact zero, so the two agree to the
   !> last bit and this only guards against a compiler reassociating them.
   real(wp), parameter :: EXACT_THR = 1.0e-14_wp

   !> Anti-vacuity floor. Every case must move the half of the tangent it
   !> exists to test by at least this much, or it is testing nothing.
   real(wp), parameter :: VACUITY_FLOOR = 1.0e-6_wp

   !> Sampling points per structure. CFC gets fewer purely because one CFC
   !> evaluation runs an `O(n_active^2)` pair kernel; see the same split in
   !> [[test_cavity_drop_lsf]].
   integer, parameter :: n_points_svdw = 3
   integer, parameter :: n_points_cfc = 1

   !> Channel combinations run at every sampling point
   integer, parameter :: n_subcases = 5

   !> Legacy SvdW smoothing, matching the LSF suite's default dispatch
   real(wp), parameter :: svdw_legacy_blend_k = 3.0_wp
   real(wp), parameter :: svdw_legacy_blend_2b = 1.0_wp
   real(wp), parameter :: svdw_legacy_blend_3b = 1.0_wp

   !> Adjoint weights and their tangents. Deterministic, non-symmetric and
   !> unrelated to one another, so no accidental cancellation hides a term.
   real(wp), parameter :: W0_FULL = 0.37_wp
   real(wp), parameter :: W1_FULL(ndim) = [0.19_wp, -0.53_wp, 0.71_wp]
   real(wp), parameter :: DW0_FULL = -0.61_wp
   real(wp), parameter :: DW1_FULL(ndim) = [0.43_wp, 0.13_wp, -0.29_wp]
   !> Induced motion of the projected point, at unit max-norm for the same
   !> reason the nuclear direction field is normalized: it is a finite-
   !> difference direction, so `STEP_SIZE` must remain the step.
   real(wp), parameter :: DR_FULL(ndim) = [1.0_wp, -0.606557377049180_wp, &
                                           0.377049180327869_wp]

contains

   !> Collect the suite
   subroutine collect_cavity_drop_field_tangent(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("svdw_field_tangent_fd", test_svdw_field_tangent), &
                  new_unittest("cfc_field_tangent_fd", test_cfc_field_tangent) &
                  ]
   end subroutine collect_cavity_drop_field_tangent

   !* ================================================================================= *!
   !*                                    Dispatches                                     *!
   !* ================================================================================= *!

   !> SvdW dispatch of the field-contraction tangent
   subroutine test_svdw_field_tangent(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Concrete level set under test
      type(moist_cavity_drop_lsf_svdw_type) :: prim

      prim%screening_threshold = 0.0_wp
      call prim%new(blend_k=svdw_legacy_blend_k, blend_2b=svdw_legacy_blend_2b, &
                    blend_3b=svdw_legacy_blend_3b)
      call run_field_tangent(prim, n_points_svdw, ABS_THR, REL_THR, "svdw", error)
   end subroutine test_svdw_field_tangent

   !> CFC dispatch of the field-contraction tangent
   subroutine test_cfc_field_tangent(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Concrete level set under test
      type(moist_cavity_drop_lsf_cfc_type) :: prim

      prim%screening_threshold = 0.0_wp
      call prim%new()
      call run_field_tangent(prim, n_points_cfc, CFC_ABS_THR, CFC_REL_THR, "cfc", error)
   end subroutine test_cfc_field_tangent

   !* ================================================================================= *!
   !*                            Shared finite-difference driver                        *!
   !* ================================================================================= *!

   !> Compare [[drop_field_tangent]] against `vjp_f1_rA` and a difference of it
   !>
   !> The reference of every case is `vjp_f1_rA(dw0, dw1, dw2)` -- exact -- plus
   !> a 4-point central difference of `vjp_f1_rA(w0, w1, w2)` taken along the
   !> combined displacement `(t*v, t*dr)` of the nuclei and the evaluation
   !> point. The two halves are added, never fitted: the case table switches
   !> channels off by zeroing their inputs, so each case compares a genuinely
   !> narrower quantity rather than a rescaled one.
   !>
   !> @param[inout] lsf      Level set under test, already constructed
   !> @param[in]    npts     Sampling points per structure
   !> @param[in]    thr_abs  Absolute comparison threshold
   !> @param[in]    thr_rel  Relative comparison threshold
   !> @param[in]    label    Dispatch name, for failure messages
   !> @param[out]   error    Error handle
   subroutine run_field_tangent(lsf, npts, thr_abs, thr_rel, label, error)
      !> Level set under test
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Sampling points per structure
      integer, intent(in) :: npts
      !> Comparison thresholds
      real(wp), intent(in) :: thr_abs, thr_rel
      !> Dispatch name
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Test fixtures
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer, allocatable :: atomic_numbers(:)

      !> Base and displaced nuclear coordinates
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      !> Nuclear direction field and its per-case mask
      real(wp), allocatable :: v_full(:, :), vdir(:, :)
      !> Analytic tangent, its reference, and the two halves of that reference
      real(wp), allocatable :: res(:, :), ref(:, :), term_a(:, :), fd_ref(:, :)
      !> The four stencil rows of the finite difference
      real(wp), allocatable :: stencil(:, :, :)
      !> Active-list identity of the base state
      integer, allocatable :: act_ref(:)

      !> Scratch of the routine under test
      type(drop_field_tangent_work_type) :: work

      !> Weights of the current case
      real(wp) :: w0, w1(ndim), w2(ndim, ndim)
      real(wp) :: dw0, dw1(ndim), dw2(ndim, ndim)
      real(wp) :: dr(ndim)
      !> Full (unmasked) Hessian weights
      real(wp) :: w2_full(ndim, ndim), dw2_full(ndim, ndim)
      !> Evaluation point and its displaced copy
      real(wp) :: point0(ndim), point_local(ndim)
      !> Stencil offsets, in units of `STEP_SIZE`
      real(wp), parameter :: stencil_t(4) = [2.0_wp, 1.0_wp, -1.0_wp, -2.0_wp]

      !> Loop indices and sizes
      integer :: icase, ipt, isub, istep, iat, i, s, j, k, nat, n_active
      !> Case name and failure message
      character(len=16) :: sub_name
      character(len=192) :: msg
      !> Magnitude of the half a case exists to exercise
      real(wp) :: probe

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         nat = mol%nat
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, npts)

         if (allocated(atomic_numbers)) deallocate (atomic_numbers)
         allocate (atomic_numbers(nat))
         do iat = 1, nat
            atomic_numbers(iat) = mol%num(mol%id(iat))
         end do

         if (allocated(centers_base)) deallocate (centers_base, centers_local, v_full, &
                                                  vdir, res, ref, term_a, fd_ref, &
                                                  stencil, act_ref)
         allocate (centers_base(ndim, nat), centers_local(ndim, nat))
         allocate (v_full(ndim, nat), vdir(ndim, nat))
         allocate (res(ndim, nat), ref(ndim, nat), term_a(ndim, nat), fd_ref(ndim, nat))
         allocate (stencil(ndim, nat, size(stencil_t)))
         allocate (act_ref(nat))
         centers_base = mol%xyz

         !* A deterministic, non-symmetric direction field, normalized to unit
         !* max-norm. The normalization is what makes the finite difference
         !* well conditioned: `STEP_SIZE` is the *step*, so with `max|v| = 1`
         !* no atom moves further than `2 * STEP_SIZE` at the outer stencil
         !* point -- exactly the displacement the per-atom nuclear FDs in
         !* [[test_cavity_drop_lsf]] take, and the scale their thresholds were
         !* measured at. An unnormalized field whose components grow with the
         !* atom index would make the effective step grow with the fixture and
         !* drown the check in truncation error on the largest structure.
         do iat = 1, nat
            v_full(1, iat) = 0.31_wp*real(mod(iat, 7), wp) - 0.7_wp
            v_full(2, iat) = -0.17_wp*real(mod(iat, 5), wp) + 0.4_wp
            v_full(3, iat) = 0.23_wp*real(mod(iat, 3) + 1, wp)
         end do
         v_full = v_full/maxval(abs(v_full))

         !* Asymmetric Hessian weights: `vjp_f1_rA` contracts all nine entries
         !* with no symmetry assumption, and the folded `outer(w1, dr)` the
         !* routine builds is asymmetric too, so a symmetric probe would not see
         !* a term that silently symmetrises them
         do k = 1, ndim
            do j = 1, ndim
               w2_full(j, k) = 0.11_wp*real(j, wp) - 0.29_wp*real(k, wp) &
                               + 0.07_wp*real(j*k, wp)
               dw2_full(j, k) = -0.23_wp*real(j, wp) + 0.17_wp*real(k, wp) &
                                + 0.13_wp*real(j*j*k, wp)
            end do
         end do

         call lsf%update(mol, radii)

         do ipt = 1, size(points, 2)
            point0 = points(:, ipt)

            do isub = 1, n_subcases
               w0 = 0.0_wp
               w1 = 0.0_wp
               w2 = 0.0_wp
               dw0 = 0.0_wp
               dw1 = 0.0_wp
               dw2 = 0.0_wp
               dr = 0.0_wp
               vdir = 0.0_wp
               select case (isub)
               case (1)
                  !* Weight tangents alone; `dr`/`v` are on but must not leak
                  dw0 = DW0_FULL
                  dw1 = DW1_FULL
                  dw2 = dw2_full
                  dr = DR_FULL
                  vdir = v_full
                  sub_name = "weights"
               case (2)
                  !* Explicit nuclear motion: the three `hvp_*` channels
                  w0 = W0_FULL
                  w1 = W1_FULL
                  w2 = w2_full
                  vdir = v_full
                  sub_name = "nuclear"
               case (3)
                  !* Point motion: two folded channels plus the `f4` one
                  w0 = W0_FULL
                  w1 = W1_FULL
                  w2 = w2_full
                  dr = DR_FULL
                  sub_name = "point"
               case (4)
                  !* Point motion through `w2` alone: the `f4_rrr_rA` channel
                  w2 = w2_full
                  dr = DR_FULL
                  sub_name = "f4"
               case default
                  !* Everything at once
                  w0 = W0_FULL
                  w1 = W1_FULL
                  w2 = w2_full
                  dw0 = DW0_FULL
                  dw1 = DW1_FULL
                  dw2 = dw2_full
                  dr = DR_FULL
                  vdir = v_full
                  sub_name = "full"
               end select

               !* ------------------------- Analytic tangent ------------------------- *!
               call prepare_at(lsf, atomic_numbers, radii, centers_base, point0, error)
               if (allocated(error)) return
               n_active = lsf%active_count()
               if (n_active == 0) cycle
               do i = 1, n_active
                  act_ref(i) = lsf%active_atom(i)
               end do

               res = 0.0_wp
               call drop_field_tangent(lsf, w0, w1, w2, dw0, dw1, dw2, dr, vdir, work, res)

               !* Term (A): exact, on the very same prepared state
               term_a = 0.0_wp
               call lsf%vjp_f1_rA(dw0, dw1, dw2, term_a)

               !* -------------------- Term (B): finite difference ------------------- *!
               fd_ref = 0.0_wp
               if (isub /= 1) then
                  do istep = 1, size(stencil_t)
                     centers_local = centers_base + (stencil_t(istep)*STEP_SIZE)*vdir
                     point_local = point0 + (stencil_t(istep)*STEP_SIZE)*dr
                     call prepare_at(lsf, atomic_numbers, radii, centers_local, &
                                     point_local, error)
                     if (allocated(error)) return
                     call check_active_list(lsf, n_active, act_ref, label, sub_name, error)
                     if (allocated(error)) return
                     stencil(:, :, istep) = 0.0_wp
                     call lsf%vjp_f1_rA(w0, w1, w2, stencil(:, :, istep))
                  end do
                  do i = 1, n_active
                     do s = 1, ndim
                        fd_ref(s, i) = fd4_scalar(stencil(s, i, 1), stencil(s, i, 2), &
                                                  stencil(s, i, 3), stencil(s, i, 4), &
                                                  STEP_SIZE)
                     end do
                  end do
               end if

               ref = term_a + fd_ref

               !* --------------------------- Anti-vacuity --------------------------- *!
               ! Each case must actually move the half it exists to exercise. A
               ! reference at machine zero would pass against an implementation
               ! that dropped the very term the case isolates.
               if (isub == 1) then
                  probe = maxval(abs(term_a(:, 1:n_active)))
               else
                  probe = maxval(abs(fd_ref(:, 1:n_active)))
               end if
               if (probe < VACUITY_FLOOR) then
                  write (msg, '(a,a,a,a,a,i0,a,i0,a,es12.5)') &
                     "vacuous field-tangent case: ", label, "/", trim(sub_name), &
                     " struct ", icase, " point ", ipt, " reference max ", probe
                  call test_failed(error, trim(msg))
                  return
               end if

               !* --------------------------- Comparison ----------------------------- *!
               do i = 1, n_active
                  do s = 1, ndim
                     write (msg, '(a,a,a,a,a,i0,a,i0,a,i0,a,i0)') &
                        "field tangent ", label, "/", trim(sub_name), &
                        ": struct ", icase, " point ", ipt, &
                        " atom ", act_ref(i), " axis ", s
                     if (isub == 1) then
                        call check(error, res(s, i), ref(s, i), thr_abs=EXACT_THR, &
                                   thr_rel=EXACT_THR, message=trim(msg))
                     else
                        call check(error, res(s, i), ref(s, i), thr_abs=thr_abs, &
                                   thr_rel=thr_rel, message=trim(msg))
                     end if
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine run_field_tangent

   !* ================================================================================= *!
   !*                                     Helpers                                       *!
   !* ================================================================================= *!

   !> Bind a displaced geometry and prepare the level set at `point`
   !>
   !> `update` refreshes the concrete's per-atom caches (and is what a nuclear
   !> finite difference needs); `set_centers` then rebuilds the spatial sort and
   !> screening bounds for exactly those coordinates. `set_max_deriv` is
   !> repeated after every `update` rather than once up front, so no accessor
   !> can silently run at an order a rebuild reset.
   !>
   !> @param[inout] lsf             Level set to prepare
   !> @param[in]    atomic_numbers  Atomic numbers of the fixture [nat]
   !> @param[in]    radii           Per-atom radii [nat]
   !> @param[in]    centers         Nuclear coordinates to bind [3, nat]
   !> @param[in]    point           Evaluation point [3]
   !> @param[out]   error           Error handle
   subroutine prepare_at(lsf, atomic_numbers, radii, centers, point, error)
      !> Level set to prepare
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Atomic numbers of the fixture
      integer, intent(in) :: atomic_numbers(:)
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Nuclear coordinates to bind
      real(wp), intent(in) :: centers(:, :)
      !> Evaluation point
      real(wp), intent(in) :: point(:)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Displaced fixture
      type(structure_type) :: mol_shift
      !> LSF-side failure
      type(mctc_error), allocatable :: lsf_err

      call new(mol_shift, atomic_numbers, centers)
      call lsf%update(mol_shift, radii)
      call lsf%set_centers(centers)
      !* `f4_rrr_rA` is the highest rung reached here: order 3 for SvdW, 4 for
      !* CFC. Ask for 4 unconditionally and let each concrete clamp it.
      call lsf%set_max_deriv(4)
      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "LSF prepare failed: "//lsf_err%message)
         return
      end if
   end subroutine prepare_at

   !> Fail unless the active list is the one the base state had
   !>
   !> Every comparison here is slot by slot, which is only meaningful while the
   !> displaced states keep the base state's active list. With
   !> `screening_threshold = 0` they do, but a silent reordering would turn a
   !> correct implementation into a confusing finite-difference failure, so it
   !> is caught as itself instead.
   !>
   !> @param[in]  lsf      Level set in its displaced state
   !> @param[in]  n_ref    Active count of the base state
   !> @param[in]  act_ref  Active list of the base state [>= n_ref]
   !> @param[in]  label    Dispatch name, for the failure message
   !> @param[in]  sub_name Case name, for the failure message
   !> @param[out] error    Error handle
   subroutine check_active_list(lsf, n_ref, act_ref, label, sub_name, error)
      !> Level set in its displaced state
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Active count of the base state
      integer, intent(in) :: n_ref
      !> Active list of the base state
      integer, intent(in) :: act_ref(:)
      !> Dispatch and case names
      character(len=*), intent(in) :: label, sub_name
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Failure message
      character(len=192) :: msg
      !> Active-slot index
      integer :: i

      if (lsf%active_count() /= n_ref) then
         write (msg, '(a,a,a,a,a,i0,a,i0)') &
            "field tangent ", label, "/", trim(sub_name), &
            ": active count changed under the FD displacement, ", n_ref, " -> ", &
            lsf%active_count()
         call test_failed(error, trim(msg))
         return
      end if
      do i = 1, n_ref
         if (lsf%active_atom(i) /= act_ref(i)) then
            write (msg, '(a,a,a,a,a,i0)') &
               "field tangent ", label, "/", trim(sub_name), &
               ": active list reordered under the FD displacement at slot ", i
            call test_failed(error, trim(msg))
            return
         end if
      end do
   end subroutine check_active_list

end module test_cavity_drop_field_tangent
