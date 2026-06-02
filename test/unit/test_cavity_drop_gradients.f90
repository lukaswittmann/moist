module test_cavity_drop_gradients
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_branching, only: branch_weight_type, softmax_weights
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_projector, only: drop_projector_type
   use moist_cavity_drop_types, only: projection_workspace_type
   use moist_radii, only : default_cpcm_radii
   use moist_data_radii_legacy, only: get_radius_func
   use mctc_io_convert, only: aatoau
   use mstore, only: get_structure
   use, intrinsic :: iso_fortran_env, only: error_unit
   implicit none
   private

   public :: collect_cavity_drop_gradients

   integer, parameter :: ndim = 3

   real(wp), parameter :: k = 2.5_wp
   real(wp), parameter :: gamma = 1.0_wp
   integer, parameter :: NUM_LEB = 50

   real(wp), parameter :: STEP_SIZE = 2.5E-4_wp
   real(wp), parameter :: ABS_THR = 5.0E-9_wp
   real(wp), parameter :: REL_THR = 5.0E-8_wp
   real(wp), parameter :: ABS_THRPRINT = 5.0E-9_wp

   real(wp), parameter :: PROJ_TOL = 1E-14_wp
   integer, parameter :: PROJ_MAXITER = 1000
   integer, parameter :: PROJ_LEVEL = 2

contains

   subroutine collect_cavity_drop_gradients(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  ! new_unittest("single_atom", test_single_atom), &
                  ! new_unittest("dimer", test_dimer), &
                  ! new_unittest("ar5_blendk_09", test_ar5_blendk_09), &
                  ! new_unittest("bih3_h2o", test_bih3_h2o), &
                  ! new_unittest("upu_0a", test_upu_0a), &
                  ! new_unittest("upu_1g", test_upu_1g), &
                  ! new_unittest("heavy28_pbh4", test_heavy28_pbh4), &
                  ! new_unittest("amino20x4_GLY_xab", test_amino20x4_gly_xab), &
                  ! new_unittest("amino20x4_TRP_xac", test_amino20x4_trp_xac), &
                  ! new_unittest("mb16-43_01", test_mb16_43_01), &
                  ! new_unittest("mb16-43_H2", test_mb16_43_h2), &
                  ! new_unittest("but14diol_32", test_but14diol_32), &
                  ! new_unittest("il16_231B", test_il16_231b) &
                  new_unittest("but14diol_1", test_but14diol_1), &
                  ! new_unittest("il16_008", test_il16_008), &
                  new_unittest("dimer_branching", test_dimer_branching), &
                  ! This test is expensive but we keep it for now
                  new_unittest("branching_xyz_totals", test_branching_xyz_totals) & 
                  ]
   end subroutine collect_cavity_drop_gradients

   !> Test gradient for a single atom
   subroutine test_single_atom(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      ! Create single oxygen atom
      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_single_atom

   !> Test gradient for a dimer
   subroutine test_dimer(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      ! Create dimer (two oxygen atoms)
      call new(mol, [6, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                     6.0_wp, 1.1_wp, 0.0_wp], [3, 2]))

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii, blend_k_override=2.0_wp, nleb_override=194, proj_level_override=2)
   end subroutine test_dimer

   !> Test gradient for 5-argon geometry with custom blend-k
   subroutine test_ar5_blendk_09(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call new(mol, [18, 18, 18, 18, 18], reshape([ &
                                                  0.2_wp, 0.0_wp, 5.1_wp, &
                                                  -2.2_wp, -2.2_wp, 0.0_wp, &
                                                  2.2_wp, -2.2_wp, 0.0_wp, &
                                                  -2.2_wp, 2.2_wp, 0.0_wp, &
                                                  2.2_wp, 2.2_wp, 0.0_wp], [3, 5]))

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii, blend_k_override=0.9_wp)
   end subroutine test_ar5_blendk_09

   !> Test gradient for bih3_h2o system
   subroutine test_bih3_h2o(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "Heavy28", "bih3_h2o")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_bih3_h2o

   !> Test gradient for UPU23 0a
   subroutine test_upu_0a(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "UPU23", "0a")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_upu_0a

   !> Test gradient for UPU23 1g
   subroutine test_upu_1g(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "UPU23", "1g")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_upu_1g

   !> Test gradient for Heavy28 h2o
   subroutine test_heavy28_h2o(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "Heavy28", "h2o")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_heavy28_h2o

   !> Test gradient for Heavy28 pbh4
   subroutine test_heavy28_pbh4(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "Heavy28", "pbh4")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_heavy28_pbh4

   !> Test gradient for Amino20x4 GLY_xab
   subroutine test_amino20x4_gly_xab(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "Amino20x4", "GLY_xab")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_amino20x4_gly_xab

   !> Test gradient for Amino20x4 TRP_xac
   subroutine test_amino20x4_trp_xac(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "Amino20x4", "TRP_xac")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_amino20x4_trp_xac

   !> Test gradient for MB16-43 01
   subroutine test_mb16_43_01(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "MB16-43", "01")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_mb16_43_01

   !> Test gradient for MB16-43 H2
   subroutine test_mb16_43_h2(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "MB16-43", "H2")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_mb16_43_h2

   !> Test gradient for But14diol 1
   subroutine test_but14diol_1(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "But14diol", "1")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_but14diol_1

   !> Test gradient for But14diol 32
   subroutine test_but14diol_32(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "But14diol", "32")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_but14diol_32

   !> Test gradient for IL16 008
   subroutine test_il16_008(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "IL16", "008")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_il16_008

   !> Branching test: carbon dimer nominally aligned with z-axis, placed
   !> near the dissociation limit (~5.75 bohr between carbons whose default
   !> radii are ~3 bohr) so that the cavity has strong concave pinch near
   !> the midpoint. A 0.01 bohr perturbation in x and y breaks perfect
   !> axial symmetry. proj_level=7 (full Lebedev onion solver) is
   !> the canonical enumerator: it exercises every LSF branch needed by
   !> the branch-weight softmax derivative path without leaning on any
   !> deflation heuristics.
   subroutine test_dimer_branching(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call new(mol, [6, 6], reshape([0.0000_wp, 0.000_wp, 0.00_wp, &
                                     0.0001_wp, 0.003_wp, 6.1_wp], [3, 2]))

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii, proj_level_override=7, nleb_override=50, blend_k_override=1.0_wp)
   end subroutine test_dimer_branching

   !> Total-area and total-volume gradient test for the branching.xyz
   !> geometry (four carbons in a planar cross plus a central carbon,
   !> coordinates in the file are Angstrom). Uses blend_k=1 and
   !> proj_level=7 (full Lebedev onion solver), which guarantees
   !> every concave-middle branch is enumerated.
   !>
   !> Unlike the per-grid-point harness in do_test, this test simply
   !> integrates cavity%a (area) and cavity%v (volume) over all grid
   !> points to get totals, and their analytic gradients by summing
   !> cavity%a_i1_rA and cavity%v_i1_rA over the grid axis. Those sums
   !> are invariant under any reordering of branch siblings or any
   !> relabelling of grid points, which is exactly the failure mode the
   !> per-gridpoint comparison would hit for a branching geometry.
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

      call new(mol, [6, 6, 6, 6, 6], reshape([ &
                                              0.00_wp,  4.21_wp,  0.00_wp, &
                                              0.00_wp,  0.00_wp,  4.22_wp, &
                                              0.00_wp, -4.18_wp,  0.00_wp, &
                                              0.00_wp,  0.00_wp, -4.15_wp, &
                                              0.02_wp,  0.10_wp, -0.20_wp], &
                                              [3, 5]) * aatoau)

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=1.0_wp, blend_3b=1.0_wp)
         call new_cavity_drop(cavity, nleb=110, &
                             do_fine=.true., &
                             tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=7, &
                             debug=.false., verbose=0, &
                             radius_model=default_cpcm_radii(), &
                             lsf_model=svdw_template, error=cavity_error)
      end block
      if (allocated(cavity_error)) call test_failed(error, cavity_error%message)

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
      call cavity%get_gradient()

      do iat = 1, mol%nat
         do idir = 1, ndim
            en_dA_drA(idir, iat) = sum(cavity%asph1_rA(idir, :, iat))
            en_dV_drA(idir, iat) = sum(cavity%vsph1_rA(idir, :, iat))
         end do
      end do

      do iat = 1, mol%nat
         do idir = 1, ndim
            call build_and_totalize(cavity, mol, iat, idir, -2.0_wp*STEP_SIZE, A_nn, V_nn, error)
            if (allocated(error)) return
            call build_and_totalize(cavity, mol, iat, idir, -1.0_wp*STEP_SIZE, A_n,  V_n, error)
            if (allocated(error)) return
            call build_and_totalize(cavity, mol, iat, idir,  1.0_wp*STEP_SIZE, A_p,  V_p, error)
            if (allocated(error)) return
            call build_and_totalize(cavity, mol, iat, idir,  2.0_wp*STEP_SIZE, A_pp, V_pp, error)
            if (allocated(error)) return

            num_dA_drA(idir, iat) = (-A_pp + 8.0_wp*A_p - 8.0_wp*A_n + A_nn) / (12.0_wp*STEP_SIZE)
            num_dV_drA(idir, iat) = (-V_pp + 8.0_wp*V_p - 8.0_wp*V_n + V_nn) / (12.0_wp*STEP_SIZE)
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
               tmp_r = phi(j);        phi(j)       = phi(j - 1);        phi(j - 1)       = tmp_r
               tmp_v = work%points(:, j); work%points(:, j) = work%points(:, j - 1); work%points(:, j - 1) = tmp_v
               tmp_r = work%rho(j);    work%rho(j)    = work%rho(j - 1);    work%rho(j - 1)    = tmp_r
               tmp_r = work%lambda(j); work%lambda(j) = work%lambda(j - 1); work%lambda(j - 1) = tmp_r
               tmp_v = work%normals(:, j); work%normals(:, j) = work%normals(:, j - 1); work%normals(:, j - 1) = tmp_v
               tmp_l = work%converged(j); work%converged(j) = work%converged(j - 1); work%converged(j - 1) = tmp_l
            else
               exit
            end if
         end do
      end do
   end subroutine sort_by_projected_owner

   !> Return the atom whose centre is closest to the given point.
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
      allocate(error)
      error%message = msg
      error%stat = 1
   end subroutine fatal_error_from_mctc

   !> Test gradient for IL16 231B
   subroutine test_il16_231b(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "IL16", "231B")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_il16_231b

   !> Test gradient of gridpoints w.r.t. atomic positions
   subroutine do_test(error, mol, radii, blend_k_override, gamma_override, nleb_override, proj_level_override, &
         branch_rho_cut_override, branch_weight_s_override)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Atomic radii
      real(wp), intent(in) :: radii(:)
      !> Optional override for blending steepness parameter
      real(wp), intent(in), optional :: blend_k_override
      !> Optional override for blending gamma parameter
      real(wp), intent(in), optional :: gamma_override
      !> Optional override for Lebedev grid size
      integer, intent(in), optional :: nleb_override
      !> Optional override for projection level
      integer, intent(in), optional :: proj_level_override
      !> Optional override for branch rho cutoff
      real(wp), intent(in), optional :: branch_rho_cut_override
      !> Optional override for branch-weight softmax scale (larger = softer)
      real(wp), intent(in), optional :: branch_weight_s_override


      type(structure_type) :: mol_fd
      type(cavity_type_drop), allocatable :: cavity
      real(wp), allocatable :: en_xyz1_rA(:, :, :, :)
      real(wp), allocatable :: num_xyz1_rA(:, :, :, :)
      real(wp), allocatable :: nn_xyz(:, :), n_xyz(:, :)
      real(wp), allocatable :: p_xyz(:, :), pp_xyz(:, :)
      integer, allocatable :: ref_numbering(:), numbering_to_idx(:)
      logical, allocatable :: valid_gridpoint_ref(:)
      logical, allocatable :: valid_gridpoint(:, :, :)  ! Validity per (idir, iat, igrid)
      integer :: iat, idir, jdir, igrid, ngrid_set, jgrid, num_idn, idx_map, max_numbering
      real(wp) :: diff, max_diff
      real(wp) :: blend_k_local
      real(wp) :: gamma_local
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

      !> Initialize cavity with configurable blending and Lebedev grid
      blend_k_local = k
      if (present(blend_k_override)) blend_k_local = blend_k_override
      gamma_local = gamma
      if (present(gamma_override)) gamma_local = gamma_override
      nleb_local = NUM_LEB
      if (present(nleb_override)) nleb_local = nleb_override
      proj_level_local = PROJ_LEVEL
      if (present(proj_level_override)) proj_level_local = proj_level_override
      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=blend_k_local, blend_3b=gamma_local)
         call new_cavity_drop(cavity, nleb=nleb_local, &
                             do_fine=.true., &
                             tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=proj_level_local, &
                             debug=.false., verbose=0, &
                             wleb_prune_level=4, &
                             radius_model=default_cpcm_radii(), &
                             lsf_model=svdw_template, error=cavity_error)
      end block
      if (allocated(cavity_error)) call test_failed(error, "Failed to initialize cavity: "//cavity_error%message)
      !> Raise wleb_cut above the tolerance-derived default (5e-16). With xi
      !> ~ 1/sqrt(wleb), points near the cutoff would otherwise inflate xi and
      !> its gradient to magnitudes where FD noise dominates REL_THR.
      if (present(branch_rho_cut_override)) &
         cavity%param%branch_rho_cut = branch_rho_cut_override
      if (present(branch_weight_s_override)) then
         cavity%param%branch_weight_s = branch_weight_s_override
         call cavity%branch_weight%init(branch_weight_s_override)
      end if
      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) call test_failed(error, "Failed to build cavity: "//cavity_error%message)
      ngrid_set = cavity%ngrid

      allocate (ref_numbering(ngrid_set))
      if (ngrid_set > 0) then
         ref_numbering = cavity%numbering(1:ngrid_set)
         max_numbering = maxval(ref_numbering)
      else
         max_numbering = 0
      end if
      allocate (numbering_to_idx(max(1, max_numbering)), source=0)
      do igrid = 1, ngrid_set
         num_idn = ref_numbering(igrid)
         if (num_idn > 0 .and. num_idn <= size(numbering_to_idx)) then
            numbering_to_idx(num_idn) = igrid
         end if
      end do

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
      call cavity%get_gradient()

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
         en_volume1_rA(:, :, idx_map) = cavity%v_i1_rA(:, :, jgrid)
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
                     write (*, '(A,I6,I6,I6,I6,4ES15.5)') 'gridpoint: ', &
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
                  write (*, '(A,I6,I6,I6,4ES15.5)') 'r_iI: ', &
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
                     write (*, '(A,I6,I6,I6,I6,3ES15.5)') 'normal: ', &
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
                  write (*, '(A,I6,I6,I6,4ES15.5)') 'cpjac: ', &
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
                  write (*, '(A,I6,I6,I6,4ES15.5)') 'wleb: ', &
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
                  write (*, '(A,I6,I6,I6,4ES15.5)') 'w_f: ', &
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
                  write (*, '(A,I6,I6,I6,4ES15.5)') 'xi: ', &
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
                  write (*, '(A,I6,I6,I6,4ES15.5)') 'iswig: ', &
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
                  write (*, '(A,I6,I6,I6,4ES15.5)') 'area: ', &
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
                  write (*, '(A,I6,I6,I6,4ES15.5)') 'volume: ', &
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

end module test_cavity_drop_gradients
