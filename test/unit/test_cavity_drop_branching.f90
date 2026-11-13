!> DROP branching tests; contains tests for the branch weights individually,
!> and full DROP nuclear gradients
!>
!> - branch weights are the pure `branch_weights_with_cut` of the projection
!>   finalize, one anchor group at a time, no cavity
!> - [[test_cross_branching]] runs the per-point derivative check of the
!>   gradient suite on the five-carbon cross, whose anchors branch, so the
!>   branch post-pass is exercised
!> - [[test_branching_xyz_totals]] checks the total area and volume gradients
!>   of the same cross against finite differences
module test_cavity_drop_branching
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use test_helpers, only: get_test_cross, fill_legacy_radii, build_numbering_map
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_branching, only: branch_weight_type
   use moist_cavity_drop_projection_finalize, only: branch_weights_with_cut
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_radii, only: default_cpcm_radii
   use moist_context, only: moist_context_type, new_context
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none(type, external)
   private

   public :: collect_cavity_drop_branching

   integer, parameter :: ndim = 3

   real(wp), parameter :: k = 2.5_wp
   real(wp), parameter :: blend_3b = 1.0_wp
   integer, parameter :: NUM_LEB = 50

   real(wp), parameter :: STEP_SIZE = 2.5E-4_wp
   real(wp), parameter :: ABS_THR = 5.0E-9_wp
   real(wp), parameter :: REL_THR = 5.0E-8_wp
   real(wp), parameter :: ABS_THRPRINT = 5.0E-9_wp

   real(wp), parameter :: PROJ_TOL = 1E-14_wp
   integer, parameter :: PROJ_MAXITER = 1000
   integer, parameter :: PROJ_LEVEL = 1

   !> Branch softmax scale for [[test_branching_xyz_totals]] so we have a lot of branching
   real(wp), parameter :: TOTALS_BRANCH_S = 0.5_wp

   !> FD step used by [[test_branching_xyz_totals]] only
   real(wp), parameter :: TOTALS_STEP_SIZE = 1.0E-3_wp

   !> Softmax weights theshold
   real(wp), parameter :: SOFTMAX_THR = 1.0e-14_wp

contains

   !> Register the branch weight and branching gradient tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_cavity_drop_branching(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("branch_partition_of_unity", test_branch_partition_of_unity), &
                  new_unittest("branch_singletons_untouched", test_branch_singletons_untouched), &
                  new_unittest("branch_uniform_phi", test_branch_uniform_phi), &
                  new_unittest("branch_dominant_sibling", test_branch_dominant_sibling), &
                  new_unittest("branch_all_below_cut_fallback", test_branch_fallback), &
                  new_unittest("branch_softmax_width_limits", test_branch_width_limits), &
                  new_unittest("branch_pruning_contract_edges", test_branch_contract_edges), &
                  new_unittest("cross_branching", test_cross_branching), &
                  new_unittest("branching_xyz_totals", test_branching_xyz_totals) &
                  ]

   end subroutine collect_cavity_drop_branching

   !* ================================================================================= *!
   !*                              branch_weights_with_cut                              *!
   !* ================================================================================= *!

   !> Branch weights of the kept siblings sum to one, both for a group that
   !> survives intact and for one that loses a sibling and is re-normalized
   subroutine test_branch_partition_of_unity(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: sigma = 1.0_wp
      real(wp) :: phi3(3), w3(3), wref(3), phi2(2), w2(2)
      integer :: n_kept

      call bw%init(sigma)

      !> Intact group: weights match the plain softmax and sum to one
      phi3 = [0.0_wp, 0.2_wp, 0.5_wp]
      call weigh(bw, phi3, 1.0e-3_wp, w3, n_kept)
      wref = ref_softmax(phi3, sigma)
      call check(error, maxval(abs(w3 - wref)) <= SOFTMAX_THR, "intact group is the plain softmax")
      if (allocated(error)) return
      call check(error, sum(w3), 1.0_wp, thr_abs=SOFTMAX_THR, thr_rel=SOFTMAX_THR)
      if (allocated(error)) return
      call check(error, n_kept, 3, "intact group keeps all siblings")
      if (allocated(error)) return

      !> Pruned group: the survivor carries the whole anchor weight
      phi2 = [0.0_wp, 12.0_wp]
      call weigh(bw, phi2, 1.0e-3_wp, w2, n_kept)
      call check(error, w2(1), 1.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, w2(2), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, n_kept, 1, "collapsed group keeps one sibling")

   end subroutine test_branch_partition_of_unity

   !> A single sibling is weight one whatever its objective
   subroutine test_branch_singletons_untouched(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp) :: phi(4), w(1)
      integer :: i, n_kept

      call bw%init(1.0_wp)
      phi = [0.0_wp, 5.0_wp, -2.0_wp, 100.0_wp]
      do i = 1, size(phi)
         call weigh(bw, phi(i:i), 1.0e-3_wp, w, n_kept)
         call check(error, w(1), 1.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
         if (allocated(error)) return
         call check(error, n_kept, 1, "singleton count")
         if (allocated(error)) return
      end do

   end subroutine test_branch_singletons_untouched

   !> Degenerate siblings share the anchor weight equally
   subroutine test_branch_uniform_phi(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp) :: phi(4), w(4)
      integer :: n_kept

      call bw%init(1.0_wp)
      phi = 0.75_wp
      call weigh(bw, phi, 1.0e-3_wp, w, n_kept)
      call check(error, all(w == 0.25_wp), "equal weights")
      if (allocated(error)) return
      call check(error, n_kept, 4, "all siblings retained")

   end subroutine test_branch_uniform_phi

   !> One clearly better branch takes the whole anchor weight
   subroutine test_branch_dominant_sibling(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp) :: w(3)
      integer :: n_kept

      call bw%init(1.0_wp)
      call weigh(bw, [0.0_wp, 30.0_wp, 30.0_wp], 1.0e-3_wp, w, n_kept)
      call check(error, w(1), 1.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, all(w(2:3) == 0.0_wp), "weak siblings cut")
      if (allocated(error)) return
      call check(error, n_kept, 1, "group collapsed to one branch")

   end subroutine test_branch_dominant_sibling

   !> When the cutoff would throw away every sibling, the strongest one is kept
   !> instead of discarding the anchor entirely
   subroutine test_branch_fallback(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp) :: w(4)
      integer :: n_kept

      call bw%init(1.0_wp)
      !> Softmax weights are 0.087/0.644/0.237/0.032, so a cutoff of 0.7 rejects
      !> all four and forces the safety fallback
      call weigh(bw, [3.0_wp, 1.0_wp, 2.0_wp, 4.0_wp], 0.7_wp, w, n_kept)
      call check(error, w(2), 1.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, count(w == 0.0_wp), 3, "the others are cut")
      if (allocated(error)) return
      call check(error, n_kept, 1, "fallback keeps a single branch")

   end subroutine test_branch_fallback

   !> The softmax width interpolates between argmax and uniform weighting
   subroutine test_branch_width_limits(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp) :: w(2)
      integer :: n_kept

      !> A zero cutoff keeps both siblings alive so the raw weights are visible
      call bw%init(0.01_wp)
      call weigh(bw, [0.0_wp, 1.0_wp], 0.0_wp, w, n_kept)
      call check(error, w(1), 1.0_wp, thr_abs=1.0e-30_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, w(2) < 1.0e-40_wp, "narrow softmax suppresses the loser")
      if (allocated(error)) return
      call check(error, w(2) > 0.0_wp, "the loser keeps a positive weight")
      if (allocated(error)) return

      !> A wide softmax washes the objective difference out
      call bw%init(1.0e6_wp)
      call weigh(bw, [0.0_wp, 1.0_wp], 0.0_wp, w, n_kept)
      call check(error, w(1), 0.5_wp, thr_abs=1.0e-6_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, w(2), 0.5_wp, thr_abs=1.0e-6_wp, thr_rel=0.0_wp)

   end subroutine test_branch_width_limits

   !> Strict cutoff with a leading discarded sibling, and exact ties at the cutoff
   subroutine test_branch_contract_edges(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp) :: w3(3), w4(4), wref(2)
      integer :: n_kept

      call bw%init(1.0_wp)
      call weigh(bw, [20.0_wp, 0.0_wp, 0.3_wp], 0.1_wp, w3, n_kept)
      wref = ref_softmax([0.0_wp, 0.3_wp], 1.0_wp)
      call check(error, w3(1), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(w3(2:3) - wref)) <= SOFTMAX_THR, "kept siblings re-softmaxed")
      if (allocated(error)) return
      call check(error, n_kept, 2, "two retained siblings")
      if (allocated(error)) return

      !> Equal weights exactly at the cutoff are all rejected, so the fallback keeps one
      !> sibling; which of the tied siblings survives is not part of the contract
      call weigh(bw, [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], 0.25_wp, w4, n_kept)
      call check(error, count(w4 == 1.0_wp), 1, "exactly one tied sibling survives")
      if (allocated(error)) return
      call check(error, count(w4 == 0.0_wp), 3, "the other tied siblings are discarded")
      if (allocated(error)) return
      call check(error, n_kept, 1, "tie fallback keeps a single branch")

   end subroutine test_branch_contract_edges

   !* ================================================================================= *!
   !*                               Branch weight helpers                               *!
   !* ================================================================================= *!

   !> Run the sibling-cut softmax with its own workspace
   !>
   !> @param[in]  bw       Softmax model
   !> @param[in]  phi      Objective per sibling
   !> @param[in]  wleb_cut Weight cutoff
   !> @param[out] w        Branch weights
   !> @param[out] n_kept   Kept siblings
   subroutine weigh(bw, phi, wleb_cut, w, n_kept)
      !> Softmax model
      type(branch_weight_type), intent(in) :: bw
      !> Objective per sibling
      real(wp), intent(in) :: phi(:)
      !> Weight cutoff
      real(wp), intent(in) :: wleb_cut
      !> Branch weights
      real(wp), intent(out) :: w(:)
      !> Kept siblings
      integer, intent(out) :: n_kept
      !> Workspace
      real(wp) :: work(2*size(phi))
      logical :: keep(size(phi))

      call branch_weights_with_cut(bw, phi, wleb_cut, w, n_kept, work, keep)

   end subroutine weigh

   !> Softmax written straight from the definition, without the shift by the
   !> minimum that the implementation uses
   pure function ref_softmax(phi, sigma) result(weights)
      !> Objective values of one branch group
      real(wp), intent(in) :: phi(:)
      !> Softmax width
      real(wp), intent(in) :: sigma
      !> Reference weights
      real(wp) :: weights(size(phi))

      weights = exp(-phi/sigma)
      weights = weights/sum(weights)

   end function ref_softmax

   !* ================================================================================= *!
   !*                                Branching gradients                                *!
   !* ================================================================================= *!

   !> Full derivative suite on a geometry whose anchors actually branch
   !>
   !> The five-carbon cross has concave seams where the projection has several
   !> minima per anchor, so the multistart projector returns siblings. Keeping
   !> them alive additionally needs a softer branch softmax: at the production
   !> scale (0.05, set by `new_cavity_drop`) the prune in filter.f90 discards
   !> every sibling but the strongest and `branch_count` collapses back to 1
   !> Measured onset is s ~ 0.2; 0.5 leaves margin so the stencil geometries
   !> branch too
   !>
   !> This is the only fixture that exercises the branch post-pass, which
   !> corrects `wleb1_rA`, `xi1_rA`, `a_i1_rA` and `v1_rA` after the main loop
   !> `xi1_rA` was missing from that list and nothing caught it
   subroutine test_cross_branching(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_test_cross(mol)

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii, proj_level_override=3, nleb_override=50, &
                   blend_k_override=1.0_wp, blend_3b_override=1.0_wp, &
                   branch_weight_s_override=0.5_wp, require_branching=.true.)
   end subroutine test_cross_branching

   !* ================================================================================= *!
   !*                              Finite-difference driver                             *!
   !* ================================================================================= *!

   !> Test gradient of gridpoints w.r.t. atomic positions
   subroutine do_test(error, mol, radii, blend_k_override, blend_3b_override, nleb_override, proj_level_override, &
                      branch_weight_s_override, require_branching)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Atomic radii
      real(wp), intent(in) :: radii(:)
      !> Optional override for blending steepness parameter
      real(wp), intent(in), optional :: blend_k_override
      !> Optional override for blending blend_3b parameter
      real(wp), intent(in), optional :: blend_3b_override
      !> Optional override for Lebedev grid size
      integer, intent(in), optional :: nleb_override
      !> Optional override for projection level
      integer, intent(in), optional :: proj_level_override
      !> Optional override for branch-weight softmax scale (larger = softer)
      real(wp), intent(in), optional :: branch_weight_s_override
      !> Require the reference build to actually produce branched anchors
      logical, intent(in), optional :: require_branching

      type(structure_type) :: mol_fd
      type(cavity_type_drop), allocatable :: cavity
      real(wp), allocatable :: en_xyz1_rA(:, :, :, :)
      real(wp), allocatable :: num_xyz1_rA(:, :, :, :)
      real(wp), allocatable :: nn_xyz(:, :), n_xyz(:, :)
      real(wp), allocatable :: p_xyz(:, :), pp_xyz(:, :)
      real(wp), allocatable :: nnn_xyz(:, :), ppp_xyz(:, :)
      integer, allocatable :: ref_numbering(:), numbering_to_idx(:)
      logical, allocatable :: valid_gridpoint_ref(:)
      logical, allocatable :: valid_gridpoint(:, :, :)  ! Validity per (idir, iat, igrid)
      !> Branch count of each reference point's anchor (ngrid_set)
      integer, allocatable :: bc_ref(:)
      !> Anchor branch count unchanged at every stencil geometry (ndim, nat, ngrid_set)
      logical, allocatable :: bc_kept(:, :, :)
      integer :: nbranched, nflipped
      integer :: iat, idir, jdir, igrid, ngrid_set, jgrid, num_idn, idx_map
      real(wp) :: diff
      real(wp) :: blend_k_local
      real(wp) :: blend_3b_local
      integer :: nleb_local, proj_level_local

      logical, allocatable :: nn_converged(:, :, :), n_converged(:, :, :), &
                              p_converged(:, :, :), pp_converged(:, :, :), &
                              nnn_converged(:, :, :), ppp_converged(:, :, :)
      integer, allocatable :: num_nn(:, :, :), num_n(:, :, :), num_p(:, :, :), num_pp(:, :, :)
      integer, allocatable :: num_nnn(:, :, :), num_ppp(:, :, :)
      type(mctc_error), allocatable :: cavity_error

      ! r_iI gradient arrays (r_iI from gridpoint to its owning atom)
      real(wp), allocatable :: en_r_iI1_rA(:, :, :)
      real(wp), allocatable :: num_r_iI1_rA(:, :, :)
      real(wp), allocatable :: nn_r_iI(:), n_r_iI(:)
      real(wp), allocatable :: p_r_iI(:), pp_r_iI(:)
      real(wp), allocatable :: nnn_r_iI(:), ppp_r_iI(:)

      ! Surface area gradient arrays
      real(wp), allocatable :: en_area1_rA(:, :, :)
      real(wp), allocatable :: num_area1_rA(:, :, :)
      real(wp), allocatable :: nn_area(:), n_area(:)
      real(wp), allocatable :: p_area(:), pp_area(:)
      real(wp), allocatable :: nnn_area(:), ppp_area(:)

      ! Switching function gradient arrays
      real(wp), allocatable :: en_iswig1_rA(:, :, :)
      real(wp), allocatable :: num_iswig1_rA(:, :, :)
      real(wp), allocatable :: nn_iswig(:), n_iswig(:)
      real(wp), allocatable :: p_iswig(:), pp_iswig(:)
      real(wp), allocatable :: nnn_iswig(:), ppp_iswig(:)

      real(wp), allocatable :: num_pou_f1_rA(:, :, :)
      real(wp), allocatable :: nn_pou_f(:), n_pou_f(:)
      real(wp), allocatable :: p_pou_f(:), pp_pou_f(:)
      real(wp), allocatable :: nnn_pou_f(:), ppp_pou_f(:)

      ! Volume gradient arrays
      real(wp), allocatable :: en_volume1_rA(:, :, :)
      real(wp), allocatable :: num_volume1_rA(:, :, :)
      real(wp), allocatable :: nn_volume(:), n_volume(:)
      real(wp), allocatable :: p_volume(:), pp_volume(:)
      real(wp), allocatable :: nnn_volume(:), ppp_volume(:)

      ! Jacobian scaling gradient arrays
      real(wp), allocatable :: en_cpjac1_rA(:, :, :)
      real(wp), allocatable :: num_cpjac1_rA(:, :, :)
      real(wp), allocatable :: nn_cpjac(:), n_cpjac(:)
      real(wp), allocatable :: p_cpjac(:), pp_cpjac(:)
      real(wp), allocatable :: nnn_cpjac(:), ppp_cpjac(:)

      ! Lebedev weight gradient arrays
      real(wp), allocatable :: en_wleb1_rA(:, :, :)
      real(wp), allocatable :: num_wleb1_rA(:, :, :)
      real(wp), allocatable :: nn_wleb(:), n_wleb(:)
      real(wp), allocatable :: p_wleb(:), pp_wleb(:)
      real(wp), allocatable :: nnn_wleb(:), ppp_wleb(:)

      ! Gaussian width (xi) gradient arrays
      real(wp), allocatable :: en_xi1_rA(:, :, :)
      real(wp), allocatable :: num_xi1_rA(:, :, :)
      real(wp), allocatable :: nn_xi(:), n_xi(:)
      real(wp), allocatable :: p_xi(:), pp_xi(:)
      real(wp), allocatable :: nnn_xi(:), ppp_xi(:)

      ! w_f switching function gradient arrays
      real(wp), allocatable :: en_w_f1_rA(:, :, :)
      real(wp), allocatable :: num_w_f1_rA(:, :, :)
      real(wp), allocatable :: nn_w_f0(:), n_w_f0(:)
      real(wp), allocatable :: p_w_f0(:), pp_w_f0(:)
      real(wp), allocatable :: nnn_w_f0(:), ppp_w_f0(:)

      ! Surface normal gradient arrays (normal has 3 components)
      real(wp), allocatable :: en_normal1_rA(:, :, :, :)
      real(wp), allocatable :: num_normal1_rA(:, :, :, :)
      real(wp), allocatable :: nn_normal(:, :), n_normal(:, :)
      real(wp), allocatable :: p_normal(:, :), pp_normal(:, :)
      real(wp), allocatable :: nnn_normal(:, :), ppp_normal(:, :)
      integer :: kdir
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      !> Initialize cavity with configurable blending and Lebedev grid
      call new_context(ctx, nthreads=0, verbosity=0)

      blend_k_local = k
      if (present(blend_k_override)) blend_k_local = blend_k_override
      blend_3b_local = blend_3b
      if (present(blend_3b_override)) blend_3b_local = blend_3b_override
      nleb_local = NUM_LEB
      if (present(nleb_override)) nleb_local = nleb_override
      proj_level_local = PROJ_LEVEL
      if (present(proj_level_override)) proj_level_local = proj_level_override
      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=blend_k_local, &
            blend_3b=blend_3b_local))
         call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
            error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=nleb_local, do_fine=.true., &
            tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=proj_level_local, wleb_prune_level=4), ctx=ctx)
      end block
      if (allocated(cavity_error)) call test_failed(error, "Failed to initialize cavity: "//cavity_error%message)
      if (allocated(error)) return
      ! Raise wleb_cut; with xi~1/sqrt(wleb) and small wleb value and gradient is increased
      ! to magnitudes where FD noise dominate
      if (present(branch_weight_s_override)) then
         cavity%param%branch_weight_s = branch_weight_s_override
         ! The admissible branch set is derived from the softmax scale, so it
         ! has to be recomputed alongside it
         call cavity%param%compute_derived(cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, "Failed to recompute derived parameters: "//cavity_error%message)
            return
         end if
         call cavity%branch_weight%init(branch_weight_s_override)
      end if
      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) call test_failed(error, "Failed to build cavity: "//cavity_error%message)
      if (allocated(error)) return
      ngrid_set = cavity%ngrid

      if (present(require_branching)) then
         if (require_branching) then
            if (.not. allocated(cavity%offset)) then
               call test_failed(error, "Fixture requested branching but has no row layout")
               return
            end if
            if (.not. any(cavity%offset(2:cavity%anchor%ngrid + 1) - cavity%offset(1:cavity%anchor%ngrid) > 1)) then
               call test_failed(error, "Fixture requested branching but no anchor branched; "// &
                                "the branch post-pass would go untested")
               return
            end if
         end if
      end if

      allocate (ref_numbering(ngrid_set))
      if (ngrid_set > 0) ref_numbering = cavity%numbering(1:ngrid_set)
      call build_numbering_map(ref_numbering, numbering_to_idx)

      ! Store reference convergence status before finite differences
      allocate (valid_gridpoint_ref(ngrid_set), source=.false.)
      allocate (bc_ref(ngrid_set), source=0)
      do jgrid = 1, cavity%ngrid
         num_idn = cavity%numbering(jgrid)
         if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
         idx_map = numbering_to_idx(num_idn)
         if (idx_map <= 0) cycle
         valid_gridpoint_ref(idx_map) = .true.
         bc_ref(idx_map) = cavity%offset(cavity%anchor_of_row(jgrid) + 1) &
                           - cavity%offset(cavity%anchor_of_row(jgrid))
      end do

      !> Get analytic gradients
      call cavity%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      if (.not. allocated(cavity%xyz1_rA)) then
         call test_failed(error, "xyz1_rA not allocated after get_gradient")
         return
      end if

      allocate (en_xyz1_rA(ndim, ndim, mol%nat, ngrid_set), source=0.0_wp)
      allocate (en_r_iI1_rA(ndim, mol%nat, ngrid_set), source=0.0_wp)
      allocate (en_iswig1_rA(ndim, mol%nat, ngrid_set), source=0.0_wp)
      allocate (en_area1_rA(ndim, mol%nat, ngrid_set), source=0.0_wp)
      allocate (en_volume1_rA(ndim, mol%nat, ngrid_set), source=0.0_wp)
      allocate (en_cpjac1_rA(ndim, mol%nat, ngrid_set), source=0.0_wp)
      allocate (en_wleb1_rA(ndim, mol%nat, ngrid_set), source=0.0_wp)
      allocate (en_w_f1_rA(ndim, mol%nat, ngrid_set), source=0.0_wp)
      allocate (en_xi1_rA(ndim, mol%nat, ngrid_set), source=0.0_wp)
      allocate (en_normal1_rA(ndim, mol%nat, ndim, ngrid_set), source=0.0_wp)

      ! Fill analytical derivatives mapped onto reference numbering
      do jgrid = 1, cavity%ngrid
         num_idn = cavity%numbering(jgrid)
         if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
         idx_map = numbering_to_idx(num_idn)
         if (idx_map <= 0) cycle
         en_xyz1_rA(:, :, :, idx_map) = cavity%xyz1_rA(:, :, :, jgrid)
         en_r_iI1_rA(:, :, idx_map) = cavity%r_iI1_rA(:, :, jgrid)
         en_iswig1_rA(:, :, idx_map) = cavity%f1_rA(:, :, jgrid)
         en_area1_rA(:, :, idx_map) = cavity%a_i1_rA(:, :, jgrid)
         en_volume1_rA(:, :, idx_map) = cavity%v1_rA(:, :, jgrid)
         en_cpjac1_rA(:, :, idx_map) = cavity%cpjac_scal1_rA(:, :, jgrid)
         en_wleb1_rA(:, :, idx_map) = cavity%wleb1_rA(:, :, jgrid)
         en_w_f1_rA(:, :, idx_map) = cavity%w_f1_rA(:, :, jgrid)
         en_xi1_rA(:, :, idx_map) = cavity%xi1_rA(:, :, jgrid)
         en_normal1_rA(:, :, :, idx_map) = cavity%normal1_rA(:, :, :, jgrid)
      end do

      allocate (num_xyz1_rA(ndim, ndim, mol%nat, ngrid_set))
      allocate (nn_xyz(ndim, ngrid_set))
      allocate (n_xyz(ndim, ngrid_set))
      allocate (p_xyz(ndim, ngrid_set))
      allocate (pp_xyz(ndim, ngrid_set))
      allocate (nnn_xyz(ndim, ngrid_set))
      allocate (ppp_xyz(ndim, ngrid_set))

      allocate (num_r_iI1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_r_iI(ngrid_set))
      allocate (n_r_iI(ngrid_set))
      allocate (p_r_iI(ngrid_set))
      allocate (pp_r_iI(ngrid_set))
      allocate (nnn_r_iI(ngrid_set))
      allocate (ppp_r_iI(ngrid_set))

      allocate (num_iswig1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_iswig(ngrid_set))
      allocate (n_iswig(ngrid_set))
      allocate (p_iswig(ngrid_set))
      allocate (pp_iswig(ngrid_set))
      allocate (nnn_iswig(ngrid_set))
      allocate (ppp_iswig(ngrid_set))

      allocate (num_pou_f1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_pou_f(ngrid_set))
      allocate (n_pou_f(ngrid_set))
      allocate (p_pou_f(ngrid_set))
      allocate (pp_pou_f(ngrid_set))
      allocate (nnn_pou_f(ngrid_set))
      allocate (ppp_pou_f(ngrid_set))

      allocate (num_area1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_area(ngrid_set))
      allocate (n_area(ngrid_set))
      allocate (p_area(ngrid_set))
      allocate (pp_area(ngrid_set))
      allocate (nnn_area(ngrid_set))
      allocate (ppp_area(ngrid_set))

      allocate (num_volume1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_volume(ngrid_set))
      allocate (n_volume(ngrid_set))
      allocate (p_volume(ngrid_set))
      allocate (pp_volume(ngrid_set))
      allocate (nnn_volume(ngrid_set))
      allocate (ppp_volume(ngrid_set))

      allocate (num_cpjac1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_cpjac(ngrid_set))
      allocate (n_cpjac(ngrid_set))
      allocate (p_cpjac(ngrid_set))
      allocate (pp_cpjac(ngrid_set))
      allocate (nnn_cpjac(ngrid_set))
      allocate (ppp_cpjac(ngrid_set))

      allocate (num_wleb1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_wleb(ngrid_set))
      allocate (n_wleb(ngrid_set))
      allocate (p_wleb(ngrid_set))
      allocate (pp_wleb(ngrid_set))
      allocate (nnn_wleb(ngrid_set))
      allocate (ppp_wleb(ngrid_set))

      allocate (num_w_f1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_w_f0(ngrid_set))
      allocate (n_w_f0(ngrid_set))
      allocate (p_w_f0(ngrid_set))
      allocate (pp_w_f0(ngrid_set))
      allocate (nnn_w_f0(ngrid_set))
      allocate (ppp_w_f0(ngrid_set))

      allocate (num_xi1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_xi(ngrid_set))
      allocate (n_xi(ngrid_set))
      allocate (p_xi(ngrid_set))
      allocate (pp_xi(ngrid_set))
      allocate (nnn_xi(ngrid_set))
      allocate (ppp_xi(ngrid_set))

      allocate (num_normal1_rA(ndim, mol%nat, ndim, ngrid_set))
      allocate (nn_normal(ndim, ngrid_set))
      allocate (n_normal(ndim, ngrid_set))
      allocate (p_normal(ndim, ngrid_set))
      allocate (pp_normal(ndim, ngrid_set))
      allocate (nnn_normal(ndim, ngrid_set))
      allocate (ppp_normal(ndim, ngrid_set))

      ! Huge tensors that contain converged and numbering info for all gridpoints
      allocate (nn_converged(ndim, mol%nat, ngrid_set), source=.false.)
      allocate (n_converged(ndim, mol%nat, ngrid_set), source=.false.)
      allocate (p_converged(ndim, mol%nat, ngrid_set), source=.false.)
      allocate (pp_converged(ndim, mol%nat, ngrid_set), source=.false.)
      allocate (nnn_converged(ndim, mol%nat, ngrid_set), source=.false.)
      allocate (ppp_converged(ndim, mol%nat, ngrid_set), source=.false.)
      allocate (bc_kept(ndim, mol%nat, ngrid_set), source=.true.)
      allocate (num_nn(ndim, mol%nat, ngrid_set), source=0)
      allocate (num_n(ndim, mol%nat, ngrid_set), source=0)
      allocate (num_p(ndim, mol%nat, ngrid_set), source=0)
      allocate (num_pp(ndim, mol%nat, ngrid_set), source=0)
      allocate (num_nnn(ndim, mol%nat, ngrid_set), source=0)
      allocate (num_ppp(ndim, mol%nat, ngrid_set), source=0)

      do iat = 1, mol%nat
         do idir = 1, ndim
            ! Initialize arrays to zero
            nn_xyz = 0.0_wp
            n_xyz = 0.0_wp
            p_xyz = 0.0_wp
            pp_xyz = 0.0_wp
            nnn_xyz = 0.0_wp
            ppp_xyz = 0.0_wp
            nn_r_iI = 0.0_wp
            n_r_iI = 0.0_wp
            p_r_iI = 0.0_wp
            pp_r_iI = 0.0_wp
            nnn_r_iI = 0.0_wp
            ppp_r_iI = 0.0_wp
            nn_iswig = 0.0_wp
            n_iswig = 0.0_wp
            p_iswig = 0.0_wp
            pp_iswig = 0.0_wp
            nnn_iswig = 0.0_wp
            ppp_iswig = 0.0_wp
            nn_pou_f = 0.0_wp
            n_pou_f = 0.0_wp
            p_pou_f = 0.0_wp
            pp_pou_f = 0.0_wp
            nnn_pou_f = 0.0_wp
            ppp_pou_f = 0.0_wp
            nn_area = 0.0_wp
            n_area = 0.0_wp
            p_area = 0.0_wp
            pp_area = 0.0_wp
            nnn_area = 0.0_wp
            ppp_area = 0.0_wp
            nn_volume = 0.0_wp
            n_volume = 0.0_wp
            p_volume = 0.0_wp
            pp_volume = 0.0_wp
            nnn_volume = 0.0_wp
            ppp_volume = 0.0_wp
            nn_cpjac = 0.0_wp
            n_cpjac = 0.0_wp
            p_cpjac = 0.0_wp
            pp_cpjac = 0.0_wp
            nnn_cpjac = 0.0_wp
            ppp_cpjac = 0.0_wp
            nn_wleb = 0.0_wp
            n_wleb = 0.0_wp
            p_wleb = 0.0_wp
            pp_wleb = 0.0_wp
            nnn_wleb = 0.0_wp
            ppp_wleb = 0.0_wp
            nn_w_f0 = 0.0_wp
            n_w_f0 = 0.0_wp
            p_w_f0 = 0.0_wp
            pp_w_f0 = 0.0_wp
            nnn_w_f0 = 0.0_wp
            ppp_w_f0 = 0.0_wp
            nn_xi = 0.0_wp
            n_xi = 0.0_wp
            p_xi = 0.0_wp
            pp_xi = 0.0_wp
            nnn_xi = 0.0_wp
            ppp_xi = 0.0_wp
            nn_normal = 0.0_wp
            n_normal = 0.0_wp
            p_normal = 0.0_wp
            pp_normal = 0.0_wp
            nnn_normal = 0.0_wp
            ppp_normal = 0.0_wp

            ! -3h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) - 3.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
            if (allocated(error)) return
            do jgrid = 1, cavity%ngrid
               num_idn = cavity%numbering(jgrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_map = numbering_to_idx(num_idn)
               if (idx_map <= 0) cycle
               nnn_xyz(:, idx_map) = cavity%xyz(:, jgrid)
               nnn_converged(idir, iat, idx_map) = .true.
               bc_kept(idir, iat, idx_map) = bc_kept(idir, iat, idx_map) .and. &
                  cavity%offset(cavity%anchor_of_row(jgrid) + 1) &
                  - cavity%offset(cavity%anchor_of_row(jgrid)) == bc_ref(idx_map)
               num_nnn(idir, iat, idx_map) = num_idn
               nnn_r_iI(idx_map) = cavity%r_iI0(jgrid)
               nnn_iswig(idx_map) = cavity%f(jgrid)
               nnn_pou_f(idx_map) = cavity%anchor%f0(cavity%anchor_of_row(jgrid))
               nnn_area(idx_map) = cavity%a(jgrid)
               nnn_volume(idx_map) = cavity%v(jgrid)
               nnn_cpjac(idx_map) = cavity%cpjac_scal0(jgrid)
               nnn_wleb(idx_map) = cavity%wleb(jgrid)
               nnn_w_f0(idx_map) = cavity%w_f0(jgrid)
               nnn_xi(idx_map) = cavity%xi0(jgrid)
               nnn_normal(:, idx_map) = cavity%normal0(:, jgrid)
            end do

            ! -2h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) - 2.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
            if (allocated(error)) return
            do jgrid = 1, cavity%ngrid
               num_idn = cavity%numbering(jgrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_map = numbering_to_idx(num_idn)
               if (idx_map <= 0) cycle
               nn_xyz(:, idx_map) = cavity%xyz(:, jgrid)
               nn_converged(idir, iat, idx_map) = .true.
               bc_kept(idir, iat, idx_map) = bc_kept(idir, iat, idx_map) .and. &
                  cavity%offset(cavity%anchor_of_row(jgrid) + 1) &
                  - cavity%offset(cavity%anchor_of_row(jgrid)) == bc_ref(idx_map)
               num_nn(idir, iat, idx_map) = num_idn
               nn_r_iI(idx_map) = cavity%r_iI0(jgrid)
               nn_iswig(idx_map) = cavity%f(jgrid)
               nn_pou_f(idx_map) = cavity%anchor%f0(cavity%anchor_of_row(jgrid))
               nn_area(idx_map) = cavity%a(jgrid)
               nn_volume(idx_map) = cavity%v(jgrid)
               nn_cpjac(idx_map) = cavity%cpjac_scal0(jgrid)
               nn_wleb(idx_map) = cavity%wleb(jgrid)
               nn_w_f0(idx_map) = cavity%w_f0(jgrid)
               nn_xi(idx_map) = cavity%xi0(jgrid)
               nn_normal(:, idx_map) = cavity%normal0(:, jgrid)
            end do

            ! -1h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) - 1.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
            if (allocated(error)) return
            do jgrid = 1, cavity%ngrid
               num_idn = cavity%numbering(jgrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_map = numbering_to_idx(num_idn)
               if (idx_map <= 0) cycle
               n_xyz(:, idx_map) = cavity%xyz(:, jgrid)
               n_converged(idir, iat, idx_map) = .true.
               bc_kept(idir, iat, idx_map) = bc_kept(idir, iat, idx_map) .and. &
                  cavity%offset(cavity%anchor_of_row(jgrid) + 1) &
                  - cavity%offset(cavity%anchor_of_row(jgrid)) == bc_ref(idx_map)
               num_n(idir, iat, idx_map) = num_idn
               n_r_iI(idx_map) = cavity%r_iI0(jgrid)
               n_iswig(idx_map) = cavity%f(jgrid)
               n_pou_f(idx_map) = cavity%anchor%f0(cavity%anchor_of_row(jgrid))
               n_area(idx_map) = cavity%a(jgrid)
               n_volume(idx_map) = cavity%v(jgrid)
               n_cpjac(idx_map) = cavity%cpjac_scal0(jgrid)
               n_wleb(idx_map) = cavity%wleb(jgrid)
               n_w_f0(idx_map) = cavity%w_f0(jgrid)
               n_xi(idx_map) = cavity%xi0(jgrid)
               n_normal(:, idx_map) = cavity%normal0(:, jgrid)
            end do

            ! +1h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) + 1.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
            if (allocated(error)) return
            do jgrid = 1, cavity%ngrid
               num_idn = cavity%numbering(jgrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_map = numbering_to_idx(num_idn)
               if (idx_map <= 0) cycle
               p_xyz(:, idx_map) = cavity%xyz(:, jgrid)
               p_converged(idir, iat, idx_map) = .true.
               bc_kept(idir, iat, idx_map) = bc_kept(idir, iat, idx_map) .and. &
                  cavity%offset(cavity%anchor_of_row(jgrid) + 1) &
                  - cavity%offset(cavity%anchor_of_row(jgrid)) == bc_ref(idx_map)
               num_p(idir, iat, idx_map) = num_idn
               p_r_iI(idx_map) = cavity%r_iI0(jgrid)
               p_iswig(idx_map) = cavity%f(jgrid)
               p_pou_f(idx_map) = cavity%anchor%f0(cavity%anchor_of_row(jgrid))
               p_area(idx_map) = cavity%a(jgrid)
               p_volume(idx_map) = cavity%v(jgrid)
               p_cpjac(idx_map) = cavity%cpjac_scal0(jgrid)
               p_wleb(idx_map) = cavity%wleb(jgrid)
               p_w_f0(idx_map) = cavity%w_f0(jgrid)
               p_xi(idx_map) = cavity%xi0(jgrid)
               p_normal(:, idx_map) = cavity%normal0(:, jgrid)
            end do

            ! +2h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) + 2.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
            if (allocated(error)) return
            do jgrid = 1, cavity%ngrid
               num_idn = cavity%numbering(jgrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_map = numbering_to_idx(num_idn)
               if (idx_map <= 0) cycle
               pp_xyz(:, idx_map) = cavity%xyz(:, jgrid)
               pp_converged(idir, iat, idx_map) = .true.
               bc_kept(idir, iat, idx_map) = bc_kept(idir, iat, idx_map) .and. &
                  cavity%offset(cavity%anchor_of_row(jgrid) + 1) &
                  - cavity%offset(cavity%anchor_of_row(jgrid)) == bc_ref(idx_map)
               num_pp(idir, iat, idx_map) = num_idn
               pp_r_iI(idx_map) = cavity%r_iI0(jgrid)
               pp_iswig(idx_map) = cavity%f(jgrid)
               pp_pou_f(idx_map) = cavity%anchor%f0(cavity%anchor_of_row(jgrid))
               pp_area(idx_map) = cavity%a(jgrid)
               pp_volume(idx_map) = cavity%v(jgrid)
               pp_cpjac(idx_map) = cavity%cpjac_scal0(jgrid)
               pp_wleb(idx_map) = cavity%wleb(jgrid)
               pp_w_f0(idx_map) = cavity%w_f0(jgrid)
               pp_xi(idx_map) = cavity%xi0(jgrid)
               pp_normal(:, idx_map) = cavity%normal0(:, jgrid)
            end do

            ! +3h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) + 3.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
            if (allocated(error)) return
            do jgrid = 1, cavity%ngrid
               num_idn = cavity%numbering(jgrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_map = numbering_to_idx(num_idn)
               if (idx_map <= 0) cycle
               ppp_xyz(:, idx_map) = cavity%xyz(:, jgrid)
               ppp_converged(idir, iat, idx_map) = .true.
               bc_kept(idir, iat, idx_map) = bc_kept(idir, iat, idx_map) .and. &
                  cavity%offset(cavity%anchor_of_row(jgrid) + 1) &
                  - cavity%offset(cavity%anchor_of_row(jgrid)) == bc_ref(idx_map)
               num_ppp(idir, iat, idx_map) = num_idn
               ppp_r_iI(idx_map) = cavity%r_iI0(jgrid)
               ppp_iswig(idx_map) = cavity%f(jgrid)
               ppp_pou_f(idx_map) = cavity%anchor%f0(cavity%anchor_of_row(jgrid))
               ppp_area(idx_map) = cavity%a(jgrid)
               ppp_volume(idx_map) = cavity%v(jgrid)
               ppp_cpjac(idx_map) = cavity%cpjac_scal0(jgrid)
               ppp_wleb(idx_map) = cavity%wleb(jgrid)
               ppp_w_f0(idx_map) = cavity%w_f0(jgrid)
               ppp_xi(idx_map) = cavity%xi0(jgrid)
               ppp_normal(:, idx_map) = cavity%normal0(:, jgrid)
            end do

            ! O(h^6) central difference, as hessian_fixed: with the 79 branched anchors of the
            ! level-3 cross the O(h^4) truncation error exceeds REL_THR (cpjac 8.6e-8 relative)
            num_xyz1_rA(:, idir, iat, :) = &
               (45.0_wp*(p_xyz(:, :) - n_xyz(:, :)) - 9.0_wp*(pp_xyz(:, :) - nn_xyz(:, :)) &
                + (ppp_xyz(:, :) - nnn_xyz(:, :))) &
               /(60.0_wp*STEP_SIZE)

            num_r_iI1_rA(idir, iat, :) = &
               (45.0_wp*(p_r_iI(:) - n_r_iI(:)) - 9.0_wp*(pp_r_iI(:) - nn_r_iI(:)) &
                + (ppp_r_iI(:) - nnn_r_iI(:))) &
               /(60.0_wp*STEP_SIZE)

            num_iswig1_rA(idir, iat, :) = &
               (45.0_wp*(p_iswig(:) - n_iswig(:)) - 9.0_wp*(pp_iswig(:) - nn_iswig(:)) &
                + (ppp_iswig(:) - nnn_iswig(:))) &
               /(60.0_wp*STEP_SIZE)

            num_pou_f1_rA(idir, iat, :) = &
               (45.0_wp*(p_pou_f(:) - n_pou_f(:)) - 9.0_wp*(pp_pou_f(:) - nn_pou_f(:)) &
                + (ppp_pou_f(:) - nnn_pou_f(:))) &
               /(60.0_wp*STEP_SIZE)

            num_area1_rA(idir, iat, :) = &
               (45.0_wp*(p_area(:) - n_area(:)) - 9.0_wp*(pp_area(:) - nn_area(:)) &
                + (ppp_area(:) - nnn_area(:))) &
               /(60.0_wp*STEP_SIZE)

            num_volume1_rA(idir, iat, :) = &
               (45.0_wp*(p_volume(:) - n_volume(:)) - 9.0_wp*(pp_volume(:) - nn_volume(:)) &
                + (ppp_volume(:) - nnn_volume(:))) &
               /(60.0_wp*STEP_SIZE)
            num_cpjac1_rA(idir, iat, :) = &
               (45.0_wp*(p_cpjac(:) - n_cpjac(:)) - 9.0_wp*(pp_cpjac(:) - nn_cpjac(:)) &
                + (ppp_cpjac(:) - nnn_cpjac(:))) &
               /(60.0_wp*STEP_SIZE)
            num_wleb1_rA(idir, iat, :) = &
               (45.0_wp*(p_wleb(:) - n_wleb(:)) - 9.0_wp*(pp_wleb(:) - nn_wleb(:)) &
                + (ppp_wleb(:) - nnn_wleb(:))) &
               /(60.0_wp*STEP_SIZE)
            num_w_f1_rA(idir, iat, :) = &
               (45.0_wp*(p_w_f0(:) - n_w_f0(:)) - 9.0_wp*(pp_w_f0(:) - nn_w_f0(:)) &
                + (ppp_w_f0(:) - nnn_w_f0(:))) &
               /(60.0_wp*STEP_SIZE)
            num_xi1_rA(idir, iat, :) = &
               (45.0_wp*(p_xi(:) - n_xi(:)) - 9.0_wp*(pp_xi(:) - nn_xi(:)) &
                + (ppp_xi(:) - nnn_xi(:))) &
               /(60.0_wp*STEP_SIZE)
            num_normal1_rA(idir, iat, :, :) = &
               (45.0_wp*(p_normal(:, :) - n_normal(:, :)) - 9.0_wp*(pp_normal(:, :) - nn_normal(:, :)) &
                + (ppp_normal(:, :) - nnn_normal(:, :))) &
               /(60.0_wp*STEP_SIZE)
         end do
      end do

      ! Build comprehensive validity tracking array
      ! A gridpoint is valid only if:
      ! 1. It converged in the reference configuration
      ! 2. It exists (numbering > 0) in all 6 FD steps
      ! 3. It converged in all 6 FD steps
      allocate (valid_gridpoint(ndim, mol%nat, ngrid_set), source=.false.)

      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Start with reference convergence
               valid_gridpoint(idir, iat, igrid) = valid_gridpoint_ref(igrid)

               ! Check all 6 FD steps: must exist (numbering > 0) AND converged
               if (valid_gridpoint(idir, iat, igrid)) then
                  valid_gridpoint(idir, iat, igrid) = &
                     (num_nn(idir, iat, igrid) > 0 .and. nn_converged(idir, iat, igrid)) .and. &
                     (num_n(idir, iat, igrid) > 0 .and. n_converged(idir, iat, igrid)) .and. &
                     (num_p(idir, iat, igrid) > 0 .and. p_converged(idir, iat, igrid)) .and. &
                     (num_pp(idir, iat, igrid) > 0 .and. pp_converged(idir, iat, igrid)) .and. &
                     (num_nnn(idir, iat, igrid) > 0 .and. nnn_converged(idir, iat, igrid)) .and. &
                     (num_ppp(idir, iat, igrid) > 0 .and. ppp_converged(idir, iat, igrid))
               end if
            end do
         end do
      end do

      !> Compare analytic vs numeric for valid gridpoints only; every nuclear direction must
      !  retain at least a single point each FD step
      !  Without this, no points left would indicate a pass!
      do iat = 1, mol%nat
         do idir = 1, ndim
            call check(error, count(valid_gridpoint(idir, iat, :)) > 0, &
                       "No valid FD comparisons for atom "//to_string(iat)// &
                       " direction "//to_string(idir))
            if (allocated(error)) return
            nflipped = count(valid_gridpoint(idir, iat, :) .and. .not. bc_kept(idir, iat, :))
            if (nflipped > 0) then
               call test_failed(error, to_string(nflipped)//" compared points change their branch "// &
                                "count across the stencil for atom "//to_string(iat)// &
                                " direction "//to_string(idir))
               return
            end if
            if (present(require_branching)) then
               if (require_branching) then
                  nbranched = count(valid_gridpoint(idir, iat, :) .and. bc_ref > 1)
                  if (nbranched < 1) then
                     call test_failed(error, "No branched point survives the stencil for atom "// &
                                      to_string(iat)//" direction "//to_string(idir))
                     return
                  end if
               end if
            end if
            do igrid = 1, ngrid_set
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle
               call check(error, &
                          all(ieee_is_finite(num_xyz1_rA(:, idir, iat, igrid))) .and. &
                          ieee_is_finite(num_r_iI1_rA(idir, iat, igrid)) .and. &
                          all(ieee_is_finite(num_normal1_rA(idir, iat, :, igrid))) .and. &
                          ieee_is_finite(num_cpjac1_rA(idir, iat, igrid)) .and. &
                          ieee_is_finite(num_wleb1_rA(idir, iat, igrid)) .and. &
                          ieee_is_finite(num_w_f1_rA(idir, iat, igrid)) .and. &
                          ieee_is_finite(num_xi1_rA(idir, iat, igrid)) .and. &
                          ieee_is_finite(num_iswig1_rA(idir, iat, igrid)) .and. &
                          ieee_is_finite(num_area1_rA(idir, iat, igrid)) .and. &
                          ieee_is_finite(num_volume1_rA(idir, iat, igrid)), &
                          "Non-finite FD reference for atom "//to_string(iat)// &
                          " direction "//to_string(idir)//" point "//to_string(igrid))
               if (allocated(error)) return
            end do
         end do
      end do

      ! Test 1: Gridpoint positions
      do iat = 1, mol%nat
         do idir = 1, ndim
            do jdir = 1, ndim
               do igrid = 1, ngrid_set
                  ! Skip invalid gridpoints (missing or failed convergence)
                  if (.not. valid_gridpoint(idir, iat, igrid)) cycle

                  diff = abs(en_xyz1_rA(jdir, idir, iat, igrid) &
                             - num_xyz1_rA(jdir, idir, iat, igrid))

                  if (diff > ABS_THRPRINT) then
                     write (*, "(A,I6,I6,I6,I6,4ES15.5)") "gridpoint: ", &
                        iat, idir, jdir, igrid, &
                        en_xyz1_rA(jdir, idir, iat, igrid), &
                        num_xyz1_rA(jdir, idir, iat, igrid), &
                        en_xyz1_rA(jdir, idir, iat, igrid) &
                        - num_xyz1_rA(jdir, idir, iat, igrid), &
                        en_xyz1_rA(jdir, idir, iat, igrid) &
                        /(en_xyz1_rA(jdir, idir, iat, igrid) &
                          - num_xyz1_rA(jdir, idir, iat, igrid))
                  end if

                  call check(error, en_xyz1_rA(jdir, idir, iat, igrid), &
                             num_xyz1_rA(jdir, idir, iat, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                             more="Grid point position gradient mismatch")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do

      ! Test 2: r_iIs
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_r_iI1_rA(idir, iat, igrid) &
                          - num_r_iI1_rA(idir, iat, igrid))

               if (diff > ABS_THRPRINT) then
                  write (*, "(A,I6,I6,I6,4ES15.5)") "r_iI: ", &
                     iat, idir, igrid, &
                     en_r_iI1_rA(idir, iat, igrid), &
                     num_r_iI1_rA(idir, iat, igrid), &
                     en_r_iI1_rA(idir, iat, igrid) &
                     - num_r_iI1_rA(idir, iat, igrid), &
                     en_r_iI1_rA(idir, iat, igrid) &
                     /(en_r_iI1_rA(idir, iat, igrid) &
                       - num_r_iI1_rA(idir, iat, igrid))
               end if

               call check(error, en_r_iI1_rA(idir, iat, igrid), &
                          num_r_iI1_rA(idir, iat, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                          more="r_iI gradient mismatch")
               if (allocated(error)) return
            end do
         end do
      end do

      ! Test 3: Surface normals
      do iat = 1, mol%nat
         do idir = 1, ndim
            do kdir = 1, ndim
               do igrid = 1, ngrid_set
                  ! Skip invalid gridpoints
                  if (.not. valid_gridpoint(idir, iat, igrid)) cycle

                  diff = abs(en_normal1_rA(kdir, iat, idir, igrid) &
                             - num_normal1_rA(idir, iat, kdir, igrid))

                  if (diff > ABS_THRPRINT) then
                     write (*, "(A,I6,I6,I6,I6,3ES15.5)") "normal: ", &
                        iat, idir, kdir, igrid, &
                        en_normal1_rA(kdir, iat, idir, igrid), &
                        num_normal1_rA(idir, iat, kdir, igrid), &
                        en_normal1_rA(kdir, iat, idir, igrid) &
                        - num_normal1_rA(idir, iat, kdir, igrid)
                  end if

                  call check(error, en_normal1_rA(kdir, iat, idir, igrid), &
                             num_normal1_rA(idir, iat, kdir, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                             more="Surface normal gradient mismatch")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do

      ! Test 7: Jacobian scaling
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_cpjac1_rA(idir, iat, igrid) &
                          - num_cpjac1_rA(idir, iat, igrid))

               if (diff > ABS_THRPRINT) then
                  write (*, "(A,I6,I6,I6,4ES15.5)") "cpjac: ", &
                     iat, idir, igrid, &
                     en_cpjac1_rA(idir, iat, igrid), &
                     num_cpjac1_rA(idir, iat, igrid), &
                     en_cpjac1_rA(idir, iat, igrid) &
                     - num_cpjac1_rA(idir, iat, igrid), &
                     en_cpjac1_rA(idir, iat, igrid) &
                     /(en_cpjac1_rA(idir, iat, igrid) &
                       - num_cpjac1_rA(idir, iat, igrid))
               end if

               call check(error, en_cpjac1_rA(idir, iat, igrid), &
                          num_cpjac1_rA(idir, iat, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                          more="Jacobian scaling gradient mismatch")
               if (allocated(error)) return
            end do
         end do
      end do

      ! Test 8: Lebedev weights
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_wleb1_rA(idir, iat, igrid) &
                          - num_wleb1_rA(idir, iat, igrid))

               if (diff > ABS_THR) then
                  write (*, "(A,I6,I6,I6,4ES15.5)") "wleb: ", &
                     iat, idir, igrid, &
                     en_wleb1_rA(idir, iat, igrid), &
                     num_wleb1_rA(idir, iat, igrid), &
                     en_wleb1_rA(idir, iat, igrid) &
                     - num_wleb1_rA(idir, iat, igrid), &
                     en_wleb1_rA(idir, iat, igrid) &
                     /(en_wleb1_rA(idir, iat, igrid) &
                       - num_wleb1_rA(idir, iat, igrid))
               end if

               ! cross_branching fails here by decision since 2026-11-03 (79-branch topology):
               ! finalize cuts rows after the softmax, so the stored branch weights of a
               ! group no longer sum to 1; do not re-pin
               call check(error, en_wleb1_rA(idir, iat, igrid), &
                          num_wleb1_rA(idir, iat, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                          more="Lebedev weight gradient mismatch")
               if (allocated(error)) return
            end do
         end do
      end do

      ! w_f switching function gradient
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_w_f1_rA(idir, iat, igrid) &
                          - num_w_f1_rA(idir, iat, igrid))

               if (diff > ABS_THRPRINT) then
                  write (*, "(A,I6,I6,I6,4ES15.5)") "w_f: ", &
                     iat, idir, igrid, &
                     en_w_f1_rA(idir, iat, igrid), &
                     num_w_f1_rA(idir, iat, igrid), &
                     en_w_f1_rA(idir, iat, igrid) &
                     - num_w_f1_rA(idir, iat, igrid), &
                     en_w_f1_rA(idir, iat, igrid) &
                     /(en_w_f1_rA(idir, iat, igrid) &
                       - num_w_f1_rA(idir, iat, igrid))
               end if

               call check(error, en_w_f1_rA(idir, iat, igrid), &
                          num_w_f1_rA(idir, iat, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                          more="w_f switching function gradient mismatch")
               if (allocated(error)) return
            end do
         end do
      end do

      ! Test 9: Gaussian width (xi) derivatives
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_xi1_rA(idir, iat, igrid) &
                          - num_xi1_rA(idir, iat, igrid))

               if (diff > ABS_THRPRINT) then
                  write (*, "(A,I6,I6,I6,4ES15.5)") "xi: ", &
                     iat, idir, igrid, &
                     en_xi1_rA(idir, iat, igrid), &
                     num_xi1_rA(idir, iat, igrid), &
                     en_xi1_rA(idir, iat, igrid) &
                     - num_xi1_rA(idir, iat, igrid), &
                     en_xi1_rA(idir, iat, igrid) &
                     /(en_xi1_rA(idir, iat, igrid) &
                       - num_xi1_rA(idir, iat, igrid))
               end if

               call check(error, en_xi1_rA(idir, iat, igrid), &
                          num_xi1_rA(idir, iat, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                          more="Gaussian width (xi) gradient mismatch")
               if (allocated(error)) return
            end do
         end do
      end do

      ! Test 3d: iSwiG switching function
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_iswig1_rA(idir, iat, igrid) &
                          - num_iswig1_rA(idir, iat, igrid))

               if (diff > ABS_THRPRINT) then
                  write (*, "(A,I6,I6,I6,4ES15.5)") "iswig: ", &
                     iat, idir, igrid, &
                     en_iswig1_rA(idir, iat, igrid), &
                     num_iswig1_rA(idir, iat, igrid), &
                     en_iswig1_rA(idir, iat, igrid) &
                     - num_iswig1_rA(idir, iat, igrid), &
                     en_iswig1_rA(idir, iat, igrid) &
                     /(en_iswig1_rA(idir, iat, igrid) &
                       - num_iswig1_rA(idir, iat, igrid))
               end if

               call check(error, en_iswig1_rA(idir, iat, igrid), &
                          num_iswig1_rA(idir, iat, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                          more="Switching function gradient mismatch")
               if (allocated(error)) return
            end do
         end do
      end do

      ! Test 4: Surface areas
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_area1_rA(idir, iat, igrid) &
                          - num_area1_rA(idir, iat, igrid))

               if (diff > ABS_THRPRINT) then
                  write (*, "(A,I6,I6,I6,4ES15.5)") "area: ", &
                     iat, idir, igrid, &
                     en_area1_rA(idir, iat, igrid), &
                     num_area1_rA(idir, iat, igrid), &
                     en_area1_rA(idir, iat, igrid) &
                     - num_area1_rA(idir, iat, igrid), &
                     en_area1_rA(idir, iat, igrid) &
                     /(en_area1_rA(idir, iat, igrid) &
                       - num_area1_rA(idir, iat, igrid))
               end if

               call check(error, en_area1_rA(idir, iat, igrid), &
                          num_area1_rA(idir, iat, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                          more="Surface area gradient mismatch")
               if (allocated(error)) return
            end do
         end do
      end do

      ! Test 6: Volumes
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_volume1_rA(idir, iat, igrid) &
                          - num_volume1_rA(idir, iat, igrid))

               if (diff > ABS_THR) then
                  write (*, "(A,I6,I6,I6,4ES15.5)") "volume: ", &
                     iat, idir, igrid, &
                     en_volume1_rA(idir, iat, igrid), &
                     num_volume1_rA(idir, iat, igrid), &
                     en_volume1_rA(idir, iat, igrid) &
                     - num_volume1_rA(idir, iat, igrid), &
                     en_volume1_rA(idir, iat, igrid) &
                     /(en_volume1_rA(idir, iat, igrid) &
                       - num_volume1_rA(idir, iat, igrid))
               end if

               call check(error, en_volume1_rA(idir, iat, igrid), &
                          num_volume1_rA(idir, iat, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                          more="Volume gradient mismatch")
               if (allocated(error)) return
            end do
         end do
      end do

   end subroutine do_test

   !* ================================================================================= *!
   !*                                  Branching totals                                 *!
   !* ================================================================================= *!

   !> Total-area and total-volume gradient test for four carbons in a planar cross plus a central carbon
   !>
   !> Softer branch softmax, so seam anchors keep their siblings; fails if no anchor branches
   subroutine test_branching_xyz_totals(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error
      real(wp), allocatable :: radii(:)
      real(wp) :: en_dA_drA(ndim, 5), en_dV_drA(ndim, 5)
      real(wp) :: num_dA_drA(ndim, 5), num_dV_drA(ndim, 5)
      real(wp) :: A_nn, A_n, A_p, A_pp
      real(wp) :: V_nn, V_n, V_p, V_pp
      integer :: iat, idir
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)

      call get_test_cross(mol)

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=1.0_wp, &
            blend_3b=1.0_wp))
         call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
            error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=110, do_fine=.true., &
            tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=3), ctx=ctx)
      end block
      if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
      if (allocated(error)) return

      !> The cavity must *borrow* the caller-owned run context, not copy it
      call check(error, associated(cavity%ctx, ctx), &
                 "DROP cavity does not borrow the caller-owned run context")
      if (allocated(error)) return
      call check(error, cavity%ctx%verbosity() == 0, &
                 "DROP cavity does not observe the context verbosity")
      if (allocated(error)) return

      cavity%param%branch_weight_s = TOTALS_BRANCH_S
      call cavity%param%compute_derived(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "Failed to recompute derived parameters: "//cavity_error%message)
         return
      end if
      call cavity%branch_weight%init(TOTALS_BRANCH_S)

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
      if (allocated(error)) return

      if (.not. any(cavity%offset(2:cavity%anchor%ngrid + 1) - cavity%offset(1:cavity%anchor%ngrid) > 1)) then
         call test_failed(error, "Fixture requested branching but no anchor branched")
         return
      end if

      call cavity%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      do iat = 1, mol%nat
         do idir = 1, ndim
            en_dA_drA(idir, iat) = sum(cavity%asph1_rA(idir, :, iat))
            en_dV_drA(idir, iat) = sum(cavity%vsph1_rA(idir, :, iat))
         end do
      end do

      do iat = 1, mol%nat
         do idir = 1, ndim
            call build_and_totalize(cavity, mol, iat, idir, -2.0_wp*TOTALS_STEP_SIZE, A_nn, V_nn, error)
            if (allocated(error)) return
            call build_and_totalize(cavity, mol, iat, idir, -1.0_wp*TOTALS_STEP_SIZE, A_n, V_n, error)
            if (allocated(error)) return
            call build_and_totalize(cavity, mol, iat, idir, 1.0_wp*TOTALS_STEP_SIZE, A_p, V_p, error)
            if (allocated(error)) return
            call build_and_totalize(cavity, mol, iat, idir, 2.0_wp*TOTALS_STEP_SIZE, A_pp, V_pp, error)
            if (allocated(error)) return

            num_dA_drA(idir, iat) = (-A_pp + 8.0_wp*A_p - 8.0_wp*A_n + A_nn)/(12.0_wp*TOTALS_STEP_SIZE)
            num_dV_drA(idir, iat) = (-V_pp + 8.0_wp*V_p - 8.0_wp*V_n + V_nn)/(12.0_wp*TOTALS_STEP_SIZE)
         end do
      end do

      do iat = 1, mol%nat
         do idir = 1, ndim
            call check(error, ieee_is_finite(num_dA_drA(idir, iat)) .and. &
                       ieee_is_finite(num_dV_drA(idir, iat)), "Non-finite total area or volume FD reference")
            if (allocated(error)) return
            call check(error, en_dA_drA(idir, iat), num_dA_drA(idir, iat), &
                       thr_abs=ABS_THR, thr_rel=REL_THR, &
                       more="Total area gradient mismatch")
            if (allocated(error)) return
            call check(error, en_dV_drA(idir, iat), num_dV_drA(idir, iat), &
                       thr_abs=ABS_THR, thr_rel=REL_THR, &
                       more="Total volume gradient mismatch")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_branching_xyz_totals

   !> Helper: rebuild cavity with atom (iat, idir) shifted by delta and
   !> return sum(a) and sum(v) over the resulting grid points
   subroutine build_and_totalize(cavity, mol, iat, idir, delta, A_tot, V_tot, error)
      type(cavity_type_drop), intent(inout) :: cavity
      type(structure_type), intent(in) :: mol
      integer, intent(in) :: iat, idir
      real(wp), intent(in) :: delta
      real(wp), intent(out) :: A_tot, V_tot
      !> Error handle: set (and caller should return) if the cavity rebuild fails
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol_fd
      type(mctc_error), allocatable :: cavity_error

      mol_fd = mol
      mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) + delta
      call cavity%update(mol_fd, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      A_tot = sum(cavity%a(1:cavity%ngrid))
      V_tot = sum(cavity%v(1:cavity%ngrid))
   end subroutine build_and_totalize

end module test_cavity_drop_branching
