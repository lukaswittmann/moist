module test_cavity_drop_gradient
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use test_helpers, only: get_test_cross, fill_legacy_radii, build_numbering_map
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_branching, only: branch_weight_type, softmax_weights
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_projector, only: drop_projector_type
   use moist_cavity_drop_types, only: projection_workspace_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_radii, only: default_cpcm_radii
   use moist_data_radii_legacy, only: get_radius_func
   use mstore, only: get_structure
   use moist_context, only: moist_context_type, new_context
   use, intrinsic :: iso_fortran_env, only: error_unit
   implicit none(type, external)
   private

   public :: collect_cavity_drop_gradient

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
   integer, parameter :: PROJ_LEVEL = 2

   !> FD step used by [[test_branching_xyz_totals]] only
   real(wp), parameter :: TOTALS_STEP_SIZE = 1.0E-3_wp

contains

   subroutine collect_cavity_drop_gradient(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  ! Cavity derivative tests
                  new_unittest("single_atom", test_single_atom), &
                  new_unittest("dimer", test_dimer), &
                  new_unittest("amino20x4_gly_xab", test_amino20x4_gly_xab), &
                  new_unittest("mb16_43_01", test_mb16_43_01), &
                  new_unittest("but14diol_1", test_but14diol_1), &
                  new_unittest("il16_008", test_il16_008), &
                  ! Adjoint tests
                  new_unittest("adjoint_channels_fd", test_adjoint_channels_fd), &
                  new_unittest("adjoint_area_channels", test_adjoint_area_channels), &
                  ! Branching tests (more expensive but kept for now)
                  new_unittest("cross_branching", test_cross_branching), &
                  new_unittest("branching_xyz_totals", test_branching_xyz_totals) &
                  ]
   end subroutine collect_cavity_drop_gradient

   !> Test gradient for a single atom
   subroutine test_single_atom(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      ! Create single oxygen atom
      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_single_atom

   !> Test gradient for a dimer
   subroutine test_dimer(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      ! Create dimer (two oxygen atoms)
      call new(mol, [6, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                     6.0_wp, 1.1_wp, 0.0_wp], [3, 2]))

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii, blend_k_override=2.0_wp, nleb_override=194, proj_level_override=2)
   end subroutine test_dimer

   !> Test gradient for Amino20x4 GLY_xab
   subroutine test_amino20x4_gly_xab(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "Amino20x4", "GLY_xab")

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_amino20x4_gly_xab

   !> Test gradient for MB16-43 01
   subroutine test_mb16_43_01(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "MB16-43", "01")

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_mb16_43_01

   !> Test gradient for But14diol 1
   subroutine test_but14diol_1(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "But14diol", "1")

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_but14diol_1

   !> Test gradient for IL16 008
   subroutine test_il16_008(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "IL16", "008")

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_il16_008

   !> Full derivative suite on a geometry whose anchors actually branch
   !>
   !> The five-carbon cross has concave seams where the projection has several
   !> minima per anchor, so the multistart projector returns siblings. Keeping
   !> them alive additionally needs a softer branch softmax: at the production
   !> scale (0.05, set by `new_cavity_drop`) the prune in filter.f90 discards
   !> every sibling but the strongest and `branch_count` collapses back to 1.
   !> Measured onset is s ~ 0.2; 0.5 leaves margin so the stencil geometries
   !> branch too.
   !>
   !> This is the only fixture that exercises the branch post-pass, which
   !> corrects `wleb1_rA`, `xi1_rA`, `a_i1_rA` and `v1_rA` after the main loop.
   !> `xi1_rA` was missing from that list and nothing caught it.
   subroutine test_cross_branching(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_test_cross(mol)

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii, proj_level_override=7, nleb_override=50, &
                   blend_k_override=1.0_wp, blend_3b_override=1.0_wp, &
                   branch_weight_s_override=0.5_wp, require_branching=.true.)
   end subroutine test_cross_branching

   !> Total-area and total-volume gradient test for four carbons in a planar cross plus a central carbon
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

      call new_context(ctx, verbosity=0)

      call get_test_cross(mol)

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=1.0_wp, blend_3b=1.0_wp)
         call new_cavity_drop(cavity, ctx, nleb=110, &
                              do_fine=.true., &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=7, &
                              radius_model=default_cpcm_radii(), &
                              lsf_model=svdw_template, error=cavity_error)
      end block
      if (allocated(cavity_error)) call test_failed(error, cavity_error%message)

      !> The cavity must *borrow* the caller-owned run context, not copy it
      call check(error, associated(cavity%ctx, ctx), &
                 "DROP cavity does not borrow the caller-owned run context")
      if (allocated(error)) return
      call check(error, cavity%ctx%verbosity == 0, &
                 "DROP cavity does not observe the context verbosity")
      if (allocated(error)) return

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) call test_failed(error, cavity_error%message)

      if (.not. allocated(cavity%branch_count)) then
         call test_failed(error, "Fixture requested branching but has no branch_count")
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
   !> return sum(a) and sum(v) over the resulting grid points.
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

   !> Sort the reference branches (and mirror the sort on phi_ref) so
   !> subsequent perturbed projections can be matched branch-by-branch
   !> via the same "nearest original atom" labelling.
   subroutine sort_by_projected_owner(work, n, mol, phi)
      type(projection_workspace_type), intent(inout) :: work
      integer, intent(in) :: n
      type(structure_type), intent(in) :: mol
      real(wp), intent(inout) :: phi(:)

      integer :: i, j, key_owner(n)
      integer :: tmp_i
      real(wp) :: tmp_r
      real(wp) :: tmp_v(3)
      logical :: tmp_l

      do i = 1, n
         key_owner(i) = nearest_atom(work%points(:, i), mol)
      end do

      ! Simple insertion sort by key_owner ascending. n is small (typically 4).
      do i = 2, n
         do j = i, 2, -1
            if (key_owner(j) < key_owner(j - 1)) then
               tmp_i = key_owner(j); key_owner(j) = key_owner(j - 1); key_owner(j - 1) = tmp_i
               tmp_r = phi(j); phi(j) = phi(j - 1); phi(j - 1) = tmp_r
               tmp_v = work%points(:, j); work%points(:, j) = work%points(:, j - 1); work%points(:, j - 1) = tmp_v
               tmp_r = work%rho(j); work%rho(j) = work%rho(j - 1); work%rho(j - 1) = tmp_r
               tmp_r = work%lambda(j); work%lambda(j) = work%lambda(j - 1); work%lambda(j - 1) = tmp_r
               tmp_v = work%normals(:, j); work%normals(:, j) = work%normals(:, j - 1); work%normals(:, j - 1) = tmp_v
               tmp_l = work%converged(j); work%converged(j) = work%converged(j - 1); work%converged(j - 1) = tmp_l
            else
               exit
            end if
         end do
      end do
   end subroutine sort_by_projected_owner

   !> Return the atom whose centers is closest to the given point.
   integer function nearest_atom(point, mol) result(idx)
      real(wp), intent(in) :: point(3)
      type(structure_type), intent(in) :: mol
      integer :: i
      real(wp) :: d, d_best

      idx = 1
      d_best = huge(1.0_wp)
      do i = 1, mol%nat
         d = sum((point - mol%xyz(:, i))**2)
         if (d < d_best) then
            d_best = d
            idx = i
         end if
      end do
   end function nearest_atom

   !> Translate an mctc_error message into the testdrive error_type so the
   !> check assertion framework can pick it up.
   subroutine fatal_error_from_mctc(error, msg)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: msg
      allocate (error)
      error%message = msg
      error%stat = 1
   end subroutine fatal_error_from_mctc

   !> Test gradient of gridpoints w.r.t. atomic positions
   subroutine do_test(error, mol, radii, blend_k_override, blend_3b_override, nleb_override, proj_level_override, &
                      branch_rho_cut_override, branch_weight_s_override, require_branching)
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
      !> Optional override for branch rho cutoff
      real(wp), intent(in), optional :: branch_rho_cut_override
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
      integer, allocatable :: ref_numbering(:), numbering_to_idx(:)
      logical, allocatable :: valid_gridpoint_ref(:)
      logical, allocatable :: valid_gridpoint(:, :, :)  ! Validity per (idir, iat, igrid)
      integer :: iat, idir, jdir, igrid, ngrid_set, jgrid, num_idn, idx_map
      real(wp) :: diff, max_diff
      real(wp) :: blend_k_local
      real(wp) :: blend_3b_local
      integer :: nleb_local, proj_level_local

      logical, allocatable :: nn_converged(:, :, :), n_converged(:, :, :), &
                              p_converged(:, :, :), pp_converged(:, :, :)
      integer, allocatable :: num_nn(:, :, :), num_n(:, :, :), num_p(:, :, :), num_pp(:, :, :)
      type(mctc_error), allocatable :: cavity_error

      ! r_iI gradient arrays (r_iI from gridpoint to its owning atom)
      real(wp), allocatable :: en_r_iI1_rA(:, :, :)
      real(wp), allocatable :: num_r_iI1_rA(:, :, :)
      real(wp), allocatable :: nn_r_iI(:), n_r_iI(:)
      real(wp), allocatable :: p_r_iI(:), pp_r_iI(:)

      ! Surface area gradient arrays
      real(wp), allocatable :: en_area1_rA(:, :, :)
      real(wp), allocatable :: num_area1_rA(:, :, :)
      real(wp), allocatable :: nn_area(:), n_area(:)
      real(wp), allocatable :: p_area(:), pp_area(:)

      ! Switching function gradient arrays
      real(wp), allocatable :: en_iswig1_rA(:, :, :)
      real(wp), allocatable :: num_iswig1_rA(:, :, :)
      real(wp), allocatable :: nn_iswig(:), n_iswig(:)
      real(wp), allocatable :: p_iswig(:), pp_iswig(:)

      real(wp), allocatable :: num_pou_f1_rA(:, :, :)
      real(wp), allocatable :: nn_pou_f(:), n_pou_f(:)
      real(wp), allocatable :: p_pou_f(:), pp_pou_f(:)

      ! Volume gradient arrays
      real(wp), allocatable :: en_volume1_rA(:, :, :)
      real(wp), allocatable :: num_volume1_rA(:, :, :)
      real(wp), allocatable :: nn_volume(:), n_volume(:)
      real(wp), allocatable :: p_volume(:), pp_volume(:)

      ! Jacobian scaling gradient arrays
      real(wp), allocatable :: en_cpjac1_rA(:, :, :)
      real(wp), allocatable :: num_cpjac1_rA(:, :, :)
      real(wp), allocatable :: nn_cpjac(:), n_cpjac(:)
      real(wp), allocatable :: p_cpjac(:), pp_cpjac(:)

      ! Lebedev weight gradient arrays
      real(wp), allocatable :: en_wleb1_rA(:, :, :)
      real(wp), allocatable :: num_wleb1_rA(:, :, :)
      real(wp), allocatable :: nn_wleb(:), n_wleb(:)
      real(wp), allocatable :: p_wleb(:), pp_wleb(:)

      ! Gaussian width (xi) gradient arrays
      real(wp), allocatable :: en_xi1_rA(:, :, :)
      real(wp), allocatable :: num_xi1_rA(:, :, :)
      real(wp), allocatable :: nn_xi(:), n_xi(:)
      real(wp), allocatable :: p_xi(:), pp_xi(:)

      ! w_f switching function gradient arrays
      real(wp), allocatable :: en_w_f1_rA(:, :, :)
      real(wp), allocatable :: num_w_f1_rA(:, :, :)
      real(wp), allocatable :: nn_w_f0(:), n_w_f0(:)
      real(wp), allocatable :: p_w_f0(:), pp_w_f0(:)

      ! Surface normal gradient arrays (normal has 3 components)
      real(wp), allocatable :: en_normal1_rA(:, :, :, :)
      real(wp), allocatable :: num_normal1_rA(:, :, :, :)
      real(wp), allocatable :: nn_normal(:, :), n_normal(:, :)
      real(wp), allocatable :: p_normal(:, :), pp_normal(:, :)
      integer :: kdir
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      !> Initialize cavity with configurable blending and Lebedev grid
      call new_context(ctx, verbosity=0)

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
         call svdw_template%new(blend_k=blend_k_local, blend_3b=blend_3b_local)
         call new_cavity_drop(cavity, ctx, nleb=nleb_local, &
                              do_fine=.true., &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=proj_level_local, &
                              wleb_prune_level=4, &
                              radius_model=default_cpcm_radii(), &
                              lsf_model=svdw_template, error=cavity_error)
      end block
      if (allocated(cavity_error)) call test_failed(error, "Failed to initialize cavity: "//cavity_error%message)
      ! Raise wleb_cut; with xi~1/sqrt(wleb) and small wleb value and gradient is increased
      ! to magnitudes where FD noise dominate
      if (present(branch_rho_cut_override)) then
         cavity%param%branch_rho_cut = branch_rho_cut_override
      end if
      if (present(branch_weight_s_override)) then
         cavity%param%branch_weight_s = branch_weight_s_override
         call cavity%branch_weight%init(branch_weight_s_override)
      end if
      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) call test_failed(error, "Failed to build cavity: "//cavity_error%message)
      ngrid_set = cavity%ngrid

      if (present(require_branching)) then
         if (require_branching) then
            if (.not. allocated(cavity%branch_count)) then
               call test_failed(error, "Fixture requested branching but has no branch_count")
               return
            end if
            if (.not. any(cavity%branch_count(1:ngrid_set) > 1)) then
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
      do jgrid = 1, cavity%ngrid
         num_idn = cavity%numbering(jgrid)
         if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
         idx_map = numbering_to_idx(num_idn)
         if (idx_map <= 0) cycle
         valid_gridpoint_ref(idx_map) = cavity%converged(jgrid)
      end do

      !> Get analytic gradients
      call cavity%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      if (.not. allocated(cavity%xyz1_rA)) then
         call test_failed(error, "xyz1_rA not allocated after get_gradient")
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

      ! Fill analytical derivatives mapped onto reference numbering.
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

      allocate (num_r_iI1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_r_iI(ngrid_set))
      allocate (n_r_iI(ngrid_set))
      allocate (p_r_iI(ngrid_set))
      allocate (pp_r_iI(ngrid_set))

      allocate (num_iswig1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_iswig(ngrid_set))
      allocate (n_iswig(ngrid_set))
      allocate (p_iswig(ngrid_set))
      allocate (pp_iswig(ngrid_set))

      allocate (num_pou_f1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_pou_f(ngrid_set))
      allocate (n_pou_f(ngrid_set))
      allocate (p_pou_f(ngrid_set))
      allocate (pp_pou_f(ngrid_set))

      allocate (num_area1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_area(ngrid_set))
      allocate (n_area(ngrid_set))
      allocate (p_area(ngrid_set))
      allocate (pp_area(ngrid_set))

      allocate (num_volume1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_volume(ngrid_set))
      allocate (n_volume(ngrid_set))
      allocate (p_volume(ngrid_set))
      allocate (pp_volume(ngrid_set))

      allocate (num_cpjac1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_cpjac(ngrid_set))
      allocate (n_cpjac(ngrid_set))
      allocate (p_cpjac(ngrid_set))
      allocate (pp_cpjac(ngrid_set))

      allocate (num_wleb1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_wleb(ngrid_set))
      allocate (n_wleb(ngrid_set))
      allocate (p_wleb(ngrid_set))
      allocate (pp_wleb(ngrid_set))

      allocate (num_w_f1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_w_f0(ngrid_set))
      allocate (n_w_f0(ngrid_set))
      allocate (p_w_f0(ngrid_set))
      allocate (pp_w_f0(ngrid_set))

      allocate (num_xi1_rA(ndim, mol%nat, ngrid_set))
      allocate (nn_xi(ngrid_set))
      allocate (n_xi(ngrid_set))
      allocate (p_xi(ngrid_set))
      allocate (pp_xi(ngrid_set))

      allocate (num_normal1_rA(ndim, mol%nat, ndim, ngrid_set))
      allocate (nn_normal(ndim, ngrid_set))
      allocate (n_normal(ndim, ngrid_set))
      allocate (p_normal(ndim, ngrid_set))
      allocate (pp_normal(ndim, ngrid_set))

      ! Huge tensors that contain converged and numbering info for all gridpoints
      allocate (nn_converged(ndim, mol%nat, ngrid_set), source=.false.)
      allocate (n_converged(ndim, mol%nat, ngrid_set), source=.false.)
      allocate (p_converged(ndim, mol%nat, ngrid_set), source=.false.)
      allocate (pp_converged(ndim, mol%nat, ngrid_set), source=.false.)
      allocate (num_nn(ndim, mol%nat, ngrid_set), source=0)
      allocate (num_n(ndim, mol%nat, ngrid_set), source=0)
      allocate (num_p(ndim, mol%nat, ngrid_set), source=0)
      allocate (num_pp(ndim, mol%nat, ngrid_set), source=0)

      do iat = 1, mol%nat
         do idir = 1, ndim
            ! Initialize arrays to zero
            nn_xyz = 0.0_wp
            n_xyz = 0.0_wp
            p_xyz = 0.0_wp
            pp_xyz = 0.0_wp
            nn_r_iI = 0.0_wp
            n_r_iI = 0.0_wp
            p_r_iI = 0.0_wp
            pp_r_iI = 0.0_wp
            nn_iswig = 0.0_wp
            n_iswig = 0.0_wp
            p_iswig = 0.0_wp
            pp_iswig = 0.0_wp
            nn_pou_f = 0.0_wp
            n_pou_f = 0.0_wp
            p_pou_f = 0.0_wp
            pp_pou_f = 0.0_wp
            nn_area = 0.0_wp
            n_area = 0.0_wp
            p_area = 0.0_wp
            pp_area = 0.0_wp
            nn_volume = 0.0_wp
            n_volume = 0.0_wp
            p_volume = 0.0_wp
            pp_volume = 0.0_wp
            nn_cpjac = 0.0_wp
            n_cpjac = 0.0_wp
            p_cpjac = 0.0_wp
            pp_cpjac = 0.0_wp
            nn_wleb = 0.0_wp
            n_wleb = 0.0_wp
            p_wleb = 0.0_wp
            pp_wleb = 0.0_wp
            nn_w_f0 = 0.0_wp
            n_w_f0 = 0.0_wp
            p_w_f0 = 0.0_wp
            pp_w_f0 = 0.0_wp
            nn_xi = 0.0_wp
            n_xi = 0.0_wp
            p_xi = 0.0_wp
            pp_xi = 0.0_wp
            nn_normal = 0.0_wp
            n_normal = 0.0_wp
            p_normal = 0.0_wp
            pp_normal = 0.0_wp

            ! -2h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) - 2.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
            do jgrid = 1, cavity%ngrid
               num_idn = cavity%numbering(jgrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_map = numbering_to_idx(num_idn)
               if (idx_map <= 0) cycle
               nn_xyz(:, idx_map) = cavity%xyz(:, jgrid)
               nn_converged(idir, iat, idx_map) = cavity%converged(jgrid)
               num_nn(idir, iat, idx_map) = num_idn
               nn_r_iI(idx_map) = cavity%r_iI0(jgrid)
               nn_iswig(idx_map) = cavity%f(jgrid)
               nn_pou_f(idx_map) = cavity%iswig_f0(jgrid)
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
            do jgrid = 1, cavity%ngrid
               num_idn = cavity%numbering(jgrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_map = numbering_to_idx(num_idn)
               if (idx_map <= 0) cycle
               n_xyz(:, idx_map) = cavity%xyz(:, jgrid)
               n_converged(idir, iat, idx_map) = cavity%converged(jgrid)
               num_n(idir, iat, idx_map) = num_idn
               n_r_iI(idx_map) = cavity%r_iI0(jgrid)
               n_iswig(idx_map) = cavity%f(jgrid)
               n_pou_f(idx_map) = cavity%iswig_f0(jgrid)
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
            do jgrid = 1, cavity%ngrid
               num_idn = cavity%numbering(jgrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_map = numbering_to_idx(num_idn)
               if (idx_map <= 0) cycle
               p_xyz(:, idx_map) = cavity%xyz(:, jgrid)
               p_converged(idir, iat, idx_map) = cavity%converged(jgrid)
               num_p(idir, iat, idx_map) = num_idn
               p_r_iI(idx_map) = cavity%r_iI0(jgrid)
               p_iswig(idx_map) = cavity%f(jgrid)
               p_pou_f(idx_map) = cavity%iswig_f0(jgrid)
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
            do jgrid = 1, cavity%ngrid
               num_idn = cavity%numbering(jgrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_map = numbering_to_idx(num_idn)
               if (idx_map <= 0) cycle
               pp_xyz(:, idx_map) = cavity%xyz(:, jgrid)
               pp_converged(idir, iat, idx_map) = cavity%converged(jgrid)
               num_pp(idir, iat, idx_map) = num_idn
               pp_r_iI(idx_map) = cavity%r_iI0(jgrid)
               pp_iswig(idx_map) = cavity%f(jgrid)
               pp_pou_f(idx_map) = cavity%iswig_f0(jgrid)
               pp_area(idx_map) = cavity%a(jgrid)
               pp_volume(idx_map) = cavity%v(jgrid)
               pp_cpjac(idx_map) = cavity%cpjac_scal0(jgrid)
               pp_wleb(idx_map) = cavity%wleb(jgrid)
               pp_w_f0(idx_map) = cavity%w_f0(jgrid)
               pp_xi(idx_map) = cavity%xi0(jgrid)
               pp_normal(:, idx_map) = cavity%normal0(:, jgrid)
            end do

            ! central diff formula: f'(x) ~= [-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)] / (12h)
            num_xyz1_rA(:, idir, iat, :) = &
               (-pp_xyz(:, :) + 8.0_wp*p_xyz(:, :) &
                - 8.0_wp*n_xyz(:, :) + nn_xyz(:, :)) &
               /(12.0_wp*STEP_SIZE)

            num_r_iI1_rA(idir, iat, :) = &
               (-pp_r_iI(:) + 8.0_wp*p_r_iI(:) &
                - 8.0_wp*n_r_iI(:) + nn_r_iI(:)) &
               /(12.0_wp*STEP_SIZE)

            num_iswig1_rA(idir, iat, :) = &
               (-pp_iswig(:) + 8.0_wp*p_iswig(:) &
                - 8.0_wp*n_iswig(:) + nn_iswig(:)) &
               /(12.0_wp*STEP_SIZE)

            num_pou_f1_rA(idir, iat, :) = &
               (-pp_pou_f(:) + 8.0_wp*p_pou_f(:) &
                - 8.0_wp*n_pou_f(:) + nn_pou_f(:)) &
               /(12.0_wp*STEP_SIZE)

            num_area1_rA(idir, iat, :) = &
               (-pp_area(:) + 8.0_wp*p_area(:) &
                - 8.0_wp*n_area(:) + nn_area(:)) &
               /(12.0_wp*STEP_SIZE)

            num_volume1_rA(idir, iat, :) = &
               (-pp_volume(:) + 8.0_wp*p_volume(:) &
                - 8.0_wp*n_volume(:) + nn_volume(:)) &
               /(12.0_wp*STEP_SIZE)
            num_cpjac1_rA(idir, iat, :) = &
               (-pp_cpjac(:) + 8.0_wp*p_cpjac(:) &
                - 8.0_wp*n_cpjac(:) + nn_cpjac(:)) &
               /(12.0_wp*STEP_SIZE)
            num_wleb1_rA(idir, iat, :) = &
               (-pp_wleb(:) + 8.0_wp*p_wleb(:) &
                - 8.0_wp*n_wleb(:) + nn_wleb(:)) &
               /(12.0_wp*STEP_SIZE)
            num_w_f1_rA(idir, iat, :) = &
               (-pp_w_f0(:) + 8.0_wp*p_w_f0(:) &
                - 8.0_wp*n_w_f0(:) + nn_w_f0(:)) &
               /(12.0_wp*STEP_SIZE)
            num_xi1_rA(idir, iat, :) = &
               (-pp_xi(:) + 8.0_wp*p_xi(:) &
                - 8.0_wp*n_xi(:) + nn_xi(:)) &
               /(12.0_wp*STEP_SIZE)
            num_normal1_rA(idir, iat, :, :) = &
               (-pp_normal(:, :) + 8.0_wp*p_normal(:, :) &
                - 8.0_wp*n_normal(:, :) + nn_normal(:, :)) &
               /(12.0_wp*STEP_SIZE)
         end do
      end do

      ! Build comprehensive validity tracking array
      ! A gridpoint is valid only if:
      ! 1. It converged in the reference configuration
      ! 2. It exists (numbering > 0) in all 4 FD steps
      ! 3. It converged in all 4 FD steps
      allocate (valid_gridpoint(ndim, mol%nat, ngrid_set), source=.false.)

      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Start with reference convergence
               valid_gridpoint(idir, iat, igrid) = valid_gridpoint_ref(igrid)

               ! Check all 4 FD steps: must exist (numbering > 0) AND converged
               if (valid_gridpoint(idir, iat, igrid)) then
                  valid_gridpoint(idir, iat, igrid) = &
                     (num_nn(idir, iat, igrid) > 0 .and. nn_converged(idir, iat, igrid)) .and. &
                     (num_n(idir, iat, igrid) > 0 .and. n_converged(idir, iat, igrid)) .and. &
                     (num_p(idir, iat, igrid) > 0 .and. p_converged(idir, iat, igrid)) .and. &
                     (num_pp(idir, iat, igrid) > 0 .and. pp_converged(idir, iat, igrid))
               end if
            end do
         end do
      end do

      !> Compare analytic vs numeric for valid gridpoints only

      ! Test 1: Gridpoint positions
      max_diff = 0.0_wp
      do iat = 1, mol%nat
         do idir = 1, ndim
            do jdir = 1, ndim
               do igrid = 1, ngrid_set
                  ! Skip invalid gridpoints (missing or failed convergence)
                  if (.not. valid_gridpoint(idir, iat, igrid)) cycle

                  diff = abs(en_xyz1_rA(jdir, idir, iat, igrid) &
                             - num_xyz1_rA(jdir, idir, iat, igrid))
                  max_diff = max(max_diff, diff)

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
      max_diff = 0.0_wp
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_r_iI1_rA(idir, iat, igrid) &
                          - num_r_iI1_rA(idir, iat, igrid))
               max_diff = max(max_diff, diff)

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
      max_diff = 0.0_wp
      do iat = 1, mol%nat
         do idir = 1, ndim
            do kdir = 1, ndim
               do igrid = 1, ngrid_set
                  ! Skip invalid gridpoints
                  if (.not. valid_gridpoint(idir, iat, igrid)) cycle

                  diff = abs(en_normal1_rA(kdir, iat, idir, igrid) &
                             - num_normal1_rA(idir, iat, kdir, igrid))
                  max_diff = max(max_diff, diff)

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
      max_diff = 0.0_wp
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_cpjac1_rA(idir, iat, igrid) &
                          - num_cpjac1_rA(idir, iat, igrid))
               max_diff = max(max_diff, diff)

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
      max_diff = 0.0_wp
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_wleb1_rA(idir, iat, igrid) &
                          - num_wleb1_rA(idir, iat, igrid))
               max_diff = max(max_diff, diff)

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

               call check(error, en_wleb1_rA(idir, iat, igrid), &
                          num_wleb1_rA(idir, iat, igrid), thr_abs=ABS_THR, thr_rel=REL_THR, &
                          more="Lebedev weight gradient mismatch")
               if (allocated(error)) return
            end do
         end do
      end do

      ! w_f switching function gradient
      max_diff = 0.0_wp
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_w_f1_rA(idir, iat, igrid) &
                          - num_w_f1_rA(idir, iat, igrid))
               max_diff = max(max_diff, diff)

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
      max_diff = 0.0_wp
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_xi1_rA(idir, iat, igrid) &
                          - num_xi1_rA(idir, iat, igrid))
               max_diff = max(max_diff, diff)

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
      max_diff = 0.0_wp
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_iswig1_rA(idir, iat, igrid) &
                          - num_iswig1_rA(idir, iat, igrid))
               max_diff = max(max_diff, diff)

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
      max_diff = 0.0_wp
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_area1_rA(idir, iat, igrid) &
                          - num_area1_rA(idir, iat, igrid))
               max_diff = max(max_diff, diff)

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
      max_diff = 0.0_wp
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid_set
               ! Skip invalid gridpoints
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               diff = abs(en_volume1_rA(idir, iat, igrid) &
                          - num_volume1_rA(idir, iat, igrid))
               max_diff = max(max_diff, diff)

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

   !> Finite-difference validation of the optional surface-adjoint channels
   !> w_n, w_k1, and w_k2 of contract_surface_lsf_weights.
   !>
   !> For a *non-owner* atom A the projection objective (the anchor term phi)
   !> carries no dependence on R_A, so displacing A perturbs the surface only
   !> through the level set field. The adjoint output (w_lsf0, w_lsf1, w_lsf2)
   !> contracted with the LSF's own nuclear-derivative jet at the projected point
   !> is then the *complete* derivative with respect to R_A of the functional the
   !> adjoint was asked about,
   !>
   !>    L_i = w_n_i . n_i + w_k1_i k1_i + w_k2_i k2_i,
   !>
   !> which is finite-differenced here by rebuilding the cavity at +-h and +-2h
   !> and reading the forward normals and principal curvatures back per grid point
   !> (matched through the persistent grid numbering, as in do_test).
   !>
   !> Every channel is contracted on its own *and* all three together: the normal
   !> channel folds a Hessian coupling (H @ dS-gradient weight) into the effective
   !> position weight that exists only when w_n is present, so a combined-only
   !> check could not separate it from the curvature contributions.
   subroutine test_adjoint_channels_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Adjoint channel combinations probed: w_n, w_k1, w_k2, and all three
      integer, parameter :: NCHAN = 4
      !> Near-umbilic guard. At k1 == k2 the individual principal curvatures are
      !> not differentiable, so such grid points are skipped for the k channels.
      real(wp), parameter :: CURV_GAP_THR = 5.0e-2_wp
      !> FD-vs-analytic thresholds for this test. The measured worst-case
      !> deviations are 6.4e-12 (w_n), 2.1e-12 (w_k1), 2.1e-12 (w_k2), and
      !> 6.3e-12 (all three) against adjoint magnitudes up to ~0.5, i.e. the
      !> agreement is limited by the 4-point FD roundoff floor at STEP_SIZE.
      !> These bounds keep roughly one order of magnitude of headroom, far below
      !> the suite-wide ABS_THR/REL_THR used for the cheaper forward gradients.
      real(wp), parameter :: ADJ_ABS = 5.0e-11_wp
      real(wp), parameter :: ADJ_REL = 5.0e-10_wp

      type(structure_type) :: mol, mol_fd
      real(wp), allocatable :: radii(:)
      type(cavity_type_drop), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error
      type(mctc_error), allocatable :: lsf_err
      type(moist_context_type), target :: ctx
      class(moist_cavity_drop_lsf_type), allocatable :: lsf

      real(wp), allocatable :: w_n(:, :), w_k1(:), w_k2(:)
      !> Per-channel surface-adjoint accumulators (all other channels zero)
      type(cavity_surface_adjoint_type) :: acc
      real(wp), allocatable :: w_lsf0(:, :), w_lsf1(:, :, :), w_lsf2(:, :, :, :)
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :), lsf3_rr_rA(:, :, :, :)
      real(wp), allocatable :: en_adj(:, :, :, :), num_adj(:, :, :, :)
      real(wp), allocatable :: ref_k1(:), ref_k2(:)
      real(wp), allocatable :: st_normal(:, :, :), st_k1(:, :), st_k2(:, :)
      integer, allocatable :: ref_owner(:), numbering_to_idx(:)
      logical, allocatable :: ref_conv(:), st_ok(:, :), valid(:, :, :)

      real(wp) :: fd_coeff(4), fd_delta(4)
      real(wp) :: d_normal(ndim), d_k1, d_k2, rhs
      integer :: iat, idir, igrid, jgrid, ich, istep, ax, a, b
      integer :: ngrid_set, nsph, num_idn, idx_map
      integer :: ncompared(NCHAN)

      fd_coeff = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/(12.0_wp*STEP_SIZE)
      fd_delta = [-2.0_wp, -1.0_wp, 1.0_wp, 2.0_wp]*STEP_SIZE

      !> Small asymmetric cluster: the two heavy centers make an elongated
      !> cavity with well-separated principal curvatures over most of the grid,
      !> and the off-axis hydrogen removes the residual rotational symmetry.
      call new(mol, [8, 6, 1], reshape([ &
                                       0.00_wp, 0.00_wp, 0.00_wp, &
                                       0.00_wp, 0.00_wp, 4.60_wp, &
                                       2.60_wp, 0.40_wp, -1.10_wp], [3, 3]))

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=k, blend_3b=blend_3b)
         call new_context(ctx, verbosity=0)
         call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                              proj_level=PROJ_LEVEL, wleb_prune_level=4, &
                              radius_model=default_cpcm_radii(), &
                              lsf_model=svdw_template, error=cavity_error)
      end block
      if (allocated(cavity_error)) then
         call test_failed(error, "Failed to initialize cavity: "//cavity_error%message)
         return
      end if
      !> The curvature channels need the forward principal curvatures; the
      !> normals are stored by the projection regardless.
      call cavity%properties(do_curvature=.true., do_normal=.true.)

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "Failed to build cavity: "//cavity_error%message)
         return
      end if

      ngrid_set = cavity%ngrid
      nsph = cavity%nsph
      if (ngrid_set == 0) then
         call test_failed(error, "empty grid")
         return
      end if

      !> Reference bookkeeping: persistent numbering -> reference grid index.
      call build_numbering_map(cavity%numbering(1:ngrid_set), numbering_to_idx)
      allocate (ref_conv(ngrid_set), source=.false.)
      allocate (ref_owner(ngrid_set), source=0)
      allocate (ref_k1(ngrid_set), source=0.0_wp)
      allocate (ref_k2(ngrid_set), source=0.0_wp)
      do igrid = 1, ngrid_set
         ref_conv(igrid) = cavity%converged(igrid)
         ref_owner(igrid) = cavity%owner(igrid)
         ref_k1(igrid) = cavity%k1(igrid)
         ref_k2(igrid) = cavity%k2(igrid)
      end do

      !> Arbitrary, distinct, non-proportional surface adjoint weights.
      allocate (w_n(ndim, ngrid_set), w_k1(ngrid_set), w_k2(ngrid_set))
      do igrid = 1, ngrid_set
         do ax = 1, ndim
            w_n(ax, igrid) = 0.4_wp + 0.1_wp*real(ax, wp) &
                             + 0.3_wp*sin(0.7_wp*real(igrid, wp))
         end do
         w_k1(igrid) = 1.0_wp + 0.5_wp*real(igrid, wp)/real(ngrid_set, wp)
         w_k2(igrid) = -0.7_wp + 0.3_wp*real(igrid, wp)/real(ngrid_set, wp)
      end do

      !> All non-probed surface channels are switched off, so w_lsf is the pure
      !> adjoint of the channel combination under test.
      allocate (w_lsf0(ngrid_set, NCHAN), source=0.0_wp)
      allocate (w_lsf1(ndim, ngrid_set, NCHAN), source=0.0_wp)
      allocate (w_lsf2(ndim, ndim, ngrid_set, NCHAN), source=0.0_wp)

      call acc%init(ngrid_set)
      acc%w_n = w_n
      call cavity%contract_surface_lsf_weights(acc, &
                                               w_lsf0(:, 1), w_lsf1(:, :, 1), w_lsf2(:, :, :, 1), &
                                               cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "w_n adjoint failed: "//cavity_error%message)
         return
      end if
      call acc%init(ngrid_set)
      acc%w_k1 = w_k1
      call cavity%contract_surface_lsf_weights(acc, &
                                               w_lsf0(:, 2), w_lsf1(:, :, 2), w_lsf2(:, :, :, 2), &
                                               cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "w_k1 adjoint failed: "//cavity_error%message)
         return
      end if
      call acc%init(ngrid_set)
      acc%w_k2 = w_k2
      call cavity%contract_surface_lsf_weights(acc, &
                                               w_lsf0(:, 3), w_lsf1(:, :, 3), w_lsf2(:, :, :, 3), &
                                               cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "w_k2 adjoint failed: "//cavity_error%message)
         return
      end if
      call acc%init(ngrid_set)
      acc%w_n = w_n
      acc%w_k1 = w_k1
      acc%w_k2 = w_k2
      call cavity%contract_surface_lsf_weights(acc, &
                                               w_lsf0(:, 4), w_lsf1(:, :, 4), w_lsf2(:, :, :, 4), &
                                               cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "combined adjoint failed: "//cavity_error%message)
         return
      end if

      !> Analytic dL/dR_A: the adjoint weights contracted with the LSF nuclear
      !> jet at the reference projected points (before any FD rebuild moves them).
      allocate (en_adj(ndim, nsph, ngrid_set, NCHAN), source=0.0_wp)
      allocate (lsf, source=cavity%lsf_model)
      call lsf%set_max_deriv(3)
      allocate (lsf1_rA(ndim, nsph), lsf2_r_rA(ndim, ndim, nsph))
      do igrid = 1, ngrid_set
         if (.not. ref_conv(igrid)) cycle
         call lsf%prepare(cavity%xyz(:, igrid), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "LSF prepare failed: "//lsf_err%message)
            return
         end if
         call lsf%f3_rr_rA_screened(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
         do ich = 1, NCHAN
            do iat = 1, nsph
               do ax = 1, ndim
                  rhs = w_lsf0(igrid, ich)*lsf1_rA(ax, iat)
                  do a = 1, ndim
                     rhs = rhs + w_lsf1(a, igrid, ich)*lsf2_r_rA(a, ax, iat)
                     do b = 1, ndim
                        rhs = rhs + w_lsf2(a, b, igrid, ich)*lsf3_rr_rA(a, b, ax, iat)
                     end do
                  end do
                  en_adj(ax, iat, igrid, ich) = rhs
               end do
            end do
         end do
      end do

      !> Numeric dL/dR_A from 4-point central differences of the forward surface.
      allocate (num_adj(ndim, nsph, ngrid_set, NCHAN), source=0.0_wp)
      allocate (valid(ndim, nsph, ngrid_set), source=.false.)
      allocate (st_normal(ndim, ngrid_set, 4), st_k1(ngrid_set, 4), st_k2(ngrid_set, 4))
      allocate (st_ok(ngrid_set, 4))

      do iat = 1, nsph
         do idir = 1, ndim
            st_normal = 0.0_wp
            st_k1 = 0.0_wp
            st_k2 = 0.0_wp
            st_ok = .false.

            do istep = 1, 4
               mol_fd = mol
               mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) + fd_delta(istep)
               call cavity%update(mol_fd, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, "FD cavity rebuild failed: "//cavity_error%message)
                  return
               end if
               do jgrid = 1, cavity%ngrid
                  num_idn = cavity%numbering(jgrid)
                  if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
                  idx_map = numbering_to_idx(num_idn)
                  if (idx_map <= 0) cycle
                  st_normal(:, idx_map, istep) = cavity%normal0(:, jgrid)
                  st_k1(idx_map, istep) = cavity%k1(jgrid)
                  st_k2(idx_map, istep) = cavity%k2(jgrid)
                  st_ok(idx_map, istep) = cavity%converged(jgrid)
               end do
            end do

            do igrid = 1, ngrid_set
               valid(idir, iat, igrid) = ref_conv(igrid) .and. all(st_ok(igrid, :))
               if (.not. valid(idir, iat, igrid)) cycle

               d_normal = matmul(st_normal(:, igrid, :), fd_coeff)
               d_k1 = dot_product(st_k1(igrid, :), fd_coeff)
               d_k2 = dot_product(st_k2(igrid, :), fd_coeff)

               num_adj(idir, iat, igrid, 1) = dot_product(w_n(:, igrid), d_normal)
               num_adj(idir, iat, igrid, 2) = w_k1(igrid)*d_k1
               num_adj(idir, iat, igrid, 3) = w_k2(igrid)*d_k2
               num_adj(idir, iat, igrid, 4) = num_adj(idir, iat, igrid, 1) &
                                              + num_adj(idir, iat, igrid, 2) &
                                              + num_adj(idir, iat, igrid, 3)
            end do
         end do
      end do

      !> Compare, skipping the owner atom (whose anchor term is not part of the
      !> level set adjoint) and near-umbilic points for the curvature channels.
      ncompared = 0
      do ich = 1, NCHAN
         do igrid = 1, ngrid_set
            if (ich > 1 .and. abs(ref_k1(igrid) - ref_k2(igrid)) < CURV_GAP_THR) cycle
            do iat = 1, nsph
               if (iat == ref_owner(igrid)) cycle
               do idir = 1, ndim
                  if (.not. valid(idir, iat, igrid)) cycle
                  ncompared(ich) = ncompared(ich) + 1
                  call check(error, en_adj(idir, iat, igrid, ich), &
                             num_adj(idir, iat, igrid, ich), &
                             thr_abs=ADJ_ABS, thr_rel=ADJ_REL, &
                             more="Surface adjoint channel "//to_string(ich)// &
                             " does not match its finite difference")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do

      !> Guard against a vacuously green run: every channel must have been
      !> exercised on a healthy number of grid points.
      do ich = 1, NCHAN
         call check(error, ncompared(ich) > 100, &
                    "Adjoint channel "//to_string(ich)//" was never compared")
         if (allocated(error)) return
      end do
   end subroutine test_adjoint_channels_fd

   !> The area and integration-weight adjoint channels must be equivalent to the
   !> `w_xi` weights they fold into.
   !>
   !> A DROP point area is `a = wleb * f / xi**2`, so a weight on `a` and a
   !> weight on `wleb` both reach the level set only through
   !> `w_xi <- -2*(a*w_a + wleb*w_w)/xi0`. Feeding `w_a`/`w_w` to the accumulator
   !> and feeding the pre-folded `w_xi` directly must therefore produce the same
   !> contraction, at all three LSF orders. Not a finite-difference test: the two
   !> paths run the same kernel on algebraically identical input, so they agree
   !> to roundoff.
   subroutine test_adjoint_area_channels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(cavity_type_drop), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error
      type(moist_context_type), target :: ctx

      !> Accumulator carrying w_a/w_w, and the pre-folded reference
      type(cavity_surface_adjoint_type) :: acc, acc_ref
      real(wp), allocatable :: w_a(:), w_w(:)
      real(wp), allocatable :: w_lsf0(:), w_lsf1(:, :), w_lsf2(:, :, :)
      real(wp), allocatable :: w_lsf0_ref(:), w_lsf1_ref(:, :), w_lsf2_ref(:, :, :)
      integer :: igrid, ngrid_set

      !> The two paths differ only by the order of a handful of flops.
      real(wp), parameter :: ADJ_EQ_THR = 1.0e-13_wp

      !> Same asymmetric cluster as `test_adjoint_channels_fd`.
      call new(mol, [8, 6, 1], reshape([ &
                                       0.00_wp, 0.00_wp, 0.00_wp, &
                                       0.00_wp, 0.00_wp, 4.60_wp, &
                                       2.60_wp, 0.40_wp, -1.10_wp], [3, 3]))

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=k, blend_3b=blend_3b)
         call new_context(ctx, verbosity=0)
         call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                              proj_level=PROJ_LEVEL, wleb_prune_level=4, &
                              radius_model=default_cpcm_radii(), &
                              lsf_model=svdw_template, error=cavity_error)
      end block
      if (allocated(cavity_error)) then
         call test_failed(error, "Failed to initialize cavity: "//cavity_error%message)
         return
      end if
      call cavity%properties(do_curvature=.true., do_normal=.true.)

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "Failed to build cavity: "//cavity_error%message)
         return
      end if

      ngrid_set = cavity%ngrid
      if (ngrid_set == 0) then
         call test_failed(error, "empty grid")
         return
      end if

      call acc%init(ngrid_set)
      call acc_ref%init(ngrid_set)
      allocate (w_a(ngrid_set), w_w(ngrid_set))
      do igrid = 1, ngrid_set
         w_a(igrid) = 0.2_wp + 0.01_wp*real(igrid, wp)
         w_w(igrid) = -0.1_wp + 0.005_wp*real(igrid, wp)
      end do
      call acc%add_surface_weights(cavity_error, w_a=w_a, w_w=w_w)
      if (allocated(cavity_error)) then
         call test_failed(error, "surface-weight accumulation failed: "//cavity_error%message)
         return
      end if
      acc_ref%w_xi = -2.0_wp*(cavity%a*acc%w_a + cavity%wleb*acc%w_w)/cavity%xi0

      allocate (w_lsf0(ngrid_set), w_lsf1(ndim, ngrid_set), &
                w_lsf2(ndim, ndim, ngrid_set))
      allocate (w_lsf0_ref(ngrid_set), w_lsf1_ref(ndim, ngrid_set), &
                w_lsf2_ref(ndim, ndim, ngrid_set))

      call cavity%contract_surface_lsf_weights(acc, w_lsf0, w_lsf1, w_lsf2, &
                                               cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "derived surface-weight contraction failed: " &
                          //cavity_error%message)
         return
      end if
      call cavity%contract_surface_lsf_weights(acc_ref, w_lsf0_ref, w_lsf1_ref, &
                                               w_lsf2_ref, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "canonical surface-weight contraction failed: " &
                          //cavity_error%message)
         return
      end if

      ! Guard against a vacuous comparison: both paths returning zero would pass.
      call check(error, maxval(abs(w_lsf0_ref)) > 1.0e-8_wp, &
                 "derived area-channel contraction is identically zero")
      if (allocated(error)) return

      call check(error, maxval(abs(w_lsf0 - w_lsf0_ref)), 0.0_wp, &
                 thr=ADJ_EQ_THR, more="derived scalar surface-weight channel")
      if (allocated(error)) return
      call check(error, maxval(abs(w_lsf1 - w_lsf1_ref)), 0.0_wp, &
                 thr=ADJ_EQ_THR, more="derived gradient surface-weight channel")
      if (allocated(error)) return
      call check(error, maxval(abs(w_lsf2 - w_lsf2_ref)), 0.0_wp, &
                 thr=ADJ_EQ_THR, more="derived Hessian surface-weight channel")
   end subroutine test_adjoint_area_channels

   !> Fill per-atom CPCM radii, turning a failed lookup into a test failure.

end module test_cavity_drop_gradient
