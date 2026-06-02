module test_cavity_drop_robustness
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use mctc_io_convert, only: aatoau
   use mctc_io_write, only: write_structure
   use mctc_io_filetype, only: filetype
   use mctc_io_symbols, only: to_symbol
   use mstore, only: get_structure
   use moist_data_radii_legacy, only: get_radius_func
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only : default_cpcm_radii
   use moist_utils_env, only: resolve_dir, ensure_dir
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private

   public :: collect_cavity_drop_robustness

   real(wp), parameter :: k = 2.5_wp
   real(wp), parameter :: gamma = 1.0_wp
   integer, parameter :: NUM_LEB = 110

contains

   subroutine collect_cavity_drop_robustness(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("upu", test_robustness_upu), &
                  new_unittest("heavy28", test_robustness_heavy28), &
                  new_unittest("amino20x4", test_robustness_amino20x4), &
                  new_unittest("mb16-43", test_robustness_mb16_43), &
                  new_unittest("but14diol", test_robustness_but14diol), &
                  new_unittest("il16", test_robustness_il16), &
                  new_unittest("fuzz", test_robustness_fuzz) &
                  ]
   end subroutine collect_cavity_drop_robustness

   !> Timing benchmark for ssd012_r vs separate calls
   subroutine test_robustness_upu(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_drop), allocatable :: cavity
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: i, iid
      type(mctc_error), allocatable :: cavity_error

      ! Array of all UPU23 IDs to iterate over
      character(len=2), parameter :: upu23_ids(24) = [ &
                                     '0a', '0b', '1a', '1b', '1c', '1e', '1f', '1g', &
                                     '1m', '1p', '2a', '2h', '2p', '3a', '3b', '3d', &
                                     '4b', '5z', '6p', '7a', '7p', '8d', '9a', 'aa']

      ! Iterate over all IDs in UPU23
      do iid = 1, size(upu23_ids)

         write (*, '(2x, A, A)') "Testing UPU23 ID: ", upu23_ids(iid)

         call get_structure(mol, "UPU23", upu23_ids(iid))

         ! Get radii
         allocate (radii(mol%nat))
         do i = 1, mol%nat
            radii(i) = get_radius_func(mol%num(mol%id(i)))
         end do

         allocate (cavity)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=k, blend_3b=gamma)
            call new_cavity_drop(cavity, nleb=NUM_LEB, &
                                debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
                                lsf_model=svdw_template, error=cavity_error)
         end block
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         if (allocated(cavity)) deallocate (cavity)

         deallocate (radii)
      end do
   end subroutine test_robustness_upu

   !> Robustness test for heavy28 set
   subroutine test_robustness_heavy28(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_drop), allocatable :: cavity
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: i, iid
      type(mctc_error), allocatable :: cavity_error

      ! Array of all heavy28 IDs to iterate over
      character(len=11), parameter :: heavy28_ids(38) = [ &
                                      'bih3       ', 'bih3_2     ', 'bih3_h2o   ', 'bih3_h2s   ', &
                                      'bih3_hbr   ', 'bih3_hcl   ', 'bih3_hi    ', 'bih3_nh3   ', &
                                      'h2o        ', 'h2s        ', 'hbr        ', 'hcl        ', &
                                      'hi         ', 'nh3        ', 'pbh4       ', 'pbh4_2     ', &
                                      'pbh4_bih3  ', 'pbh4_h2o   ', 'pbh4_hbr   ', 'pbh4_hcl   ', &
                                      'pbh4_hi    ', 'pbh4_teh2  ', 'sbh3       ', 'sbh3_2     ', &
                                      'sbh3_h2o   ', 'sbh3_h2s   ', 'sbh3_hbr   ', 'sbh3_hcl   ', &
                                      'sbh3_hi    ', 'sbh3_nh3   ', 'teh2       ', 'teh2_2     ', &
                                      'teh2_h2o   ', 'teh2_h2s   ', 'teh2_hbr   ', 'teh2_hcl   ', &
                                      'teh2_hi    ', 'teh2_nh3   ']

      ! Iterate over all IDs in heavy28
      do iid = 1, size(heavy28_ids)

         write (*, '(2x, A, A)') "Testing heavy28 ID: ", trim(heavy28_ids(iid))

         call get_structure(mol, "Heavy28", trim(heavy28_ids(iid)))

         ! Get radii
         allocate (radii(mol%nat))
         do i = 1, mol%nat
            radii(i) = get_radius_func(mol%num(mol%id(i)))
         end do

         allocate (cavity)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=k, blend_3b=gamma)
            call new_cavity_drop(cavity, nleb=NUM_LEB, &
                                debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
                                lsf_model=svdw_template, error=cavity_error)
         end block
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         if (allocated(cavity)) deallocate (cavity)

         deallocate (radii)
      end do
   end subroutine test_robustness_heavy28

   !> Robustness test for amino20x4 set
   subroutine test_robustness_amino20x4(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_drop), allocatable :: cavity
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: i, iid
      type(mctc_error), allocatable :: cavity_error

      ! Array of all amino20x4 IDs to iterate over
      character(len=7), parameter :: amino20x4_ids(100) = [ &
                                     'ALA_xab', 'ALA_xac', 'ALA_xag', 'ALA_xai', 'ALA_xak', &
                                     'ARG_xak', 'ARG_xbv', 'ARG_xbx', 'ARG_xby', 'ARG_xci', &
                                     'ASN_xab', 'ASN_xae', 'ASN_xaf', 'ASN_xah', 'ASN_xaj', &
                                     'ASP_xad', 'ASP_xau', 'ASP_xay', 'ASP_xaz', 'ASP_xbc', &
                                     'CYS_xag', 'CYS_xah', 'CYS_xai', 'CYS_xal', 'CYS_xao', &
                                     'GLN_xai', 'GLN_xal', 'GLN_xan', 'GLN_xap', 'GLN_xat', &
                                     'GLU_xad', 'GLU_xal', 'GLU_xar', 'GLU_xav', 'GLU_xbi', &
                                     'GLY_xab', 'GLY_xac', 'GLY_xad', 'GLY_xae', 'GLY_xag', &
                                     'HIS_xah', 'HIS_xam', 'HIS_xaq', 'HIS_xau', 'HIS_xav', &
                                     'ILE_xae', 'ILE_xag', 'ILE_xaj', 'ILE_xak', 'ILE_xaq', &
                                     'LEU_xad', 'LEU_xae', 'LEU_xap', 'LEU_xaq', 'LEU_xbb', &
                                     'LYS_xan', 'LYS_xao', 'LYS_xap', 'LYS_xas', 'LYS_xat', &
                                     'MET_xag', 'MET_xav', 'MET_xbf', 'MET_xbm', 'MET_xbo', &
                                     'PHE_xab', 'PHE_xal', 'PHE_xan', 'PHE_xar', 'PHE_xaw', &
                                     'PRO_xab', 'PRO_xac', 'PRO_xad', 'PRO_xae', 'PRO_xaf', &
                                     'SER_xad', 'SER_xaf', 'SER_xah', 'SER_xak', 'SER_xar', &
                                     'THR_xab', 'THR_xag', 'THR_xah', 'THR_xal', 'THR_xaq', &
                                     'TRP_xac', 'TRP_xaf', 'TRP_xag', 'TRP_xah', 'TRP_xao', &
                                     'TYR_xab', 'TYR_xag', 'TYR_xah', 'TYR_xan', 'TYR_xar', &
                                     'VAL_xad', 'VAL_xaf', 'VAL_xah', 'VAL_xaj', 'VAL_xak']

      ! Iterate over all IDs in amino20x4
      do iid = 1, size(amino20x4_ids)

         write (*, '(2x, A, A)') "Testing amino20x4 ID: ", amino20x4_ids(iid)

         call get_structure(mol, "Amino20x4", amino20x4_ids(iid))

         ! Get radii
         allocate (radii(mol%nat))
         do i = 1, mol%nat
            radii(i) = get_radius_func(mol%num(mol%id(i)))
         end do

         allocate (cavity)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=k, blend_3b=gamma)
            call new_cavity_drop(cavity, nleb=NUM_LEB, &
                                debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
                                lsf_model=svdw_template, error=cavity_error)
         end block
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         if (allocated(cavity)) deallocate (cavity)

         deallocate (radii)
      end do
   end subroutine test_robustness_amino20x4

   !> Robustness test for mb16-43 set
   subroutine test_robustness_mb16_43(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_drop), allocatable :: cavity
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: i, iid
      type(mctc_error), allocatable :: cavity_error

      ! Array of all mb16-43 IDs to iterate over
      character(len=4), parameter :: mb16_43_ids(59) = [ &
                                     '01  ', '02  ', '03  ', '04  ', '05  ', '06  ', '07  ', '08  ', &
                                     '09  ', '10  ', '11  ', '12  ', '13  ', '14  ', '15  ', '16  ', &
                                     '17  ', '18  ', '19  ', '20  ', '21  ', '22  ', '23  ', '24  ', &
                                     '25  ', '26  ', '27  ', '28  ', '29  ', '30  ', '31  ', '32  ', &
                                     '33  ', '34  ', '35  ', '36  ', '37  ', '38  ', '39  ', '40  ', &
                                     '41  ', '42  ', '43  ', 'AlH3', 'BH3 ', 'BeH2', 'CH4 ', 'Cl2 ', &
                                     'F2  ', 'H2  ', 'LiH ', 'MgH2', 'N2  ', 'NaH ', 'O2  ', 'P2  ', &
                                     'S2  ', 'PCl ', 'SiH4']

      ! Iterate over all IDs in mb16-43
      do iid = 1, size(mb16_43_ids)

         write (*, '(2x, A, A)') "Testing MB16-43 ID: ", trim(mb16_43_ids(iid))

         call get_structure(mol, "MB16-43", trim(mb16_43_ids(iid)))

         ! Get radii
         allocate (radii(mol%nat))
         do i = 1, mol%nat
            radii(i) = get_radius_func(mol%num(mol%id(i)))
         end do

         allocate (cavity)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=k, blend_3b=gamma)
            call new_cavity_drop(cavity, nleb=NUM_LEB, &
                                debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
                                lsf_model=svdw_template, error=cavity_error)
         end block
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         if (allocated(cavity)) deallocate (cavity)

         deallocate (radii)
      end do
   end subroutine test_robustness_mb16_43

   !> Robustness test for But14diol set
   subroutine test_robustness_but14diol(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_drop), allocatable :: cavity
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: i, iid
      type(mctc_error), allocatable :: cavity_error

      ! Array of all But14diol IDs to iterate over
      character(len=2), parameter :: but14diol_ids(65) = [ &
                                     '1 ', '2 ', '3 ', '4 ', '5 ', '6 ', '7 ', '8 ', '9 ', '10', &
                                     '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', &
                                     '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', &
                                     '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', &
                                     '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', &
                                     '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', &
                                     '61', '62', '63', '64', '65']

      ! Iterate over all IDs in But14diol
      do iid = 1, size(but14diol_ids)

         write (*, '(2x, A, A)') "Testing But14diol ID: ", trim(but14diol_ids(iid))

         call get_structure(mol, "But14diol", trim(but14diol_ids(iid)))

         ! Get radii
         allocate (radii(mol%nat))
         do i = 1, mol%nat
            radii(i) = get_radius_func(mol%num(mol%id(i)))
         end do

         allocate (cavity)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=k, blend_3b=gamma)
            call new_cavity_drop(cavity, nleb=NUM_LEB, &
                                debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
                                lsf_model=svdw_template, error=cavity_error)
         end block
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         if (allocated(cavity)) deallocate (cavity)

         deallocate (radii)
      end do
   end subroutine test_robustness_but14diol

   !> Robustness test for IL16 set
   subroutine test_robustness_il16(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_drop), allocatable :: cavity
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: i, iid
      type(mctc_error), allocatable :: cavity_error

      ! Array of all IL16 IDs to iterate over
      character(len=4), parameter :: il16_ids(48) = [ &
                                     '008 ', '008A', '008B', '144 ', '144A', '144B', '147 ', '147A', &
                                     '147B', '148 ', '148A', '148B', '150 ', '150A', '150B', '152 ', &
                                     '152A', '152B', '187 ', '187A', '187B', '202 ', '202A', '202B', &
                                     '212 ', '212A', '212B', '213 ', '213A', '213B', '214 ', '214A', &
                                     '214B', '227 ', '227A', '227B', '228 ', '228A', '228B', '229 ', &
                                     '229A', '229B', '230 ', '230A', '230B', '231 ', '231A', '231B']

      ! Iterate over all IDs in IL16
      do iid = 1, size(il16_ids)

         write (*, '(2x, A, A)') "Testing IL16 ID: ", trim(il16_ids(iid))

         call get_structure(mol, "IL16", trim(il16_ids(iid)))

         ! Get radii
         allocate (radii(mol%nat))
         do i = 1, mol%nat
            radii(i) = get_radius_func(mol%num(mol%id(i)))
         end do

         allocate (cavity)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=k, blend_3b=gamma)
            call new_cavity_drop(cavity, nleb=NUM_LEB, &
                                debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
                                lsf_model=svdw_template, error=cavity_error)
         end block
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         if (allocated(cavity)) deallocate (cavity)

         deallocate (radii)
      end do
   end subroutine test_robustness_il16

   !> Fuzz testing with random molecular structures
   subroutine test_robustness_fuzz(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_drop), allocatable :: cavity
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), xyz(:, :)
      integer, allocatable :: num(:)
      integer :: i, itest, nat, iat, iat_ref
      real(wp) :: rand_val, distance, min_dist, vdw_radius
      logical :: too_close
      integer :: seed_array(8)
      character(len=300) :: test_xyz_path
      character(len=:), allocatable :: issues_dir
      type(mctc_error), allocatable :: cavity_error

      ! Fuzz test settings

      !> Number of tests
      integer, parameter :: NUM_TESTS = 100

      !> Min and max number of atoms
      integer, parameter :: MIN_ATOMS = 5
      integer, parameter :: MAX_ATOMS = 175
      !> Min and max atomic numbers
      integer, parameter :: MIN_ATNUM = 1
      integer, parameter :: MAX_ATNUM = 86

      !> Probability to place new atom near existing ones
      real(wp), parameter :: VICINITY_PROB = 0.97_wp
      !> Minimum distance factor to avoid extreme overlaps
      real(wp), parameter :: MIN_DIST_FACTOR = 0.5_wp
      !> Vicinity radius for placing atoms near existing ones
      real(wp), parameter :: VICINITY_RADIUS = 3.0_wp*aatoau
      !> Size of the far away box
      real(wp), parameter :: FAR_BOX_SIZE = 20.0_wp*aatoau

      !> Seed value for rng
      integer, parameter :: SEED_VALUE = 42

      ! Initialize random seed for reproducibility
      seed_array = SEED_VALUE
      call random_seed(put=seed_array)

      write (*, '(2x, A)') "Starting fuzz testing with random structures..."

      ! Resolve directory for problematic-structure dumps
      ! (MOIST_DROP_ISSUES_DIR override, repo-relative default)
      issues_dir = resolve_dir("MOIST_DROP_ISSUES_DIR", "test/outputs/drop_issues")
      call ensure_dir(issues_dir)

      ! Main fuzz test loop
      do itest = 1, NUM_TESTS

         ! Generate random number of atoms
         call random_number(rand_val)
         nat = MIN_ATOMS + int(rand_val*(MAX_ATOMS - MIN_ATOMS + 1))

         ! Allocate arrays
         allocate (xyz(3, nat), num(nat), radii(nat))

         ! Generate random structure
         do iat = 1, nat

            ! Random atomic number
            call random_number(rand_val)
            num(iat) = MIN_ATNUM + int(rand_val*(MAX_ATNUM - MIN_ATNUM + 1))

            ! Get VdW radius for this atom
            vdw_radius = get_radius_func(num(iat))
            radii(iat) = vdw_radius

            ! Generate position: 90% near existing atoms, 10% far away
            if (iat == 1) then
               ! First atom at origin with small random offset
               call random_number(xyz(:, iat))
               xyz(:, iat) = (xyz(:, iat) - 0.5_wp)*2.0_wp*aatoau
            else
               call random_number(rand_val)

               if (rand_val < VICINITY_PROB) then
                  ! Place near an existing atom
                  call random_number(rand_val)
                  iat_ref = 1 + int(rand_val*(iat - 1))

                  ! Random direction
                  call random_number(xyz(:, iat))
                  xyz(:, iat) = xyz(:, iat) - 0.5_wp
                  distance = sqrt(sum(xyz(:, iat)**2))
                  if (distance > 1.0e-10_wp) then
                     xyz(:, iat) = xyz(:, iat)/distance
                  else
                     xyz(:, iat) = [1.0_wp, 0.0_wp, 0.0_wp]
                  end if

                  ! Random distance within vicinity
                  call random_number(rand_val)
                  distance = MIN_DIST_FACTOR*(vdw_radius + radii(iat_ref)) + &
                             rand_val*VICINITY_RADIUS

                  xyz(:, iat) = xyz(:, iat_ref) + xyz(:, iat)*distance
               else
                  ! Place far away in box
                  call random_number(xyz(:, iat))
                  xyz(:, iat) = (xyz(:, iat) - 0.5_wp)*FAR_BOX_SIZE
               end if

               ! Check minimum distance constraint (avoid extreme overlaps)
               too_close = .false.
               do i = 1, iat - 1
                  distance = sqrt(sum((xyz(:, iat) - xyz(:, i))**2))
                  min_dist = MIN_DIST_FACTOR*(radii(iat) + radii(i))
                  if (distance < min_dist) then
                     too_close = .true.
                     exit
                  end if
               end do

               ! If too close, try to place it elsewhere (simple retry)
               if (too_close) then
                  call random_number(xyz(:, iat))
                  xyz(:, iat) = (xyz(:, iat) - 0.5_wp)*FAR_BOX_SIZE
               end if
            end if
         end do

         ! Create structure_type
         mol%nat = nat
         mol%nid = 0
         mol%nbd = 0
         mol%uhf = 0
         mol%charge = 0.0_wp

         ! Allocate and set atomic data
         if (allocated(mol%num)) deallocate (mol%num)
         if (allocated(mol%xyz)) deallocate (mol%xyz)
         if (allocated(mol%id)) deallocate (mol%id)
         if (allocated(mol%sym)) deallocate (mol%sym)
         allocate (mol%num(nat), mol%xyz(3, nat), mol%id(nat), mol%sym(nat))
         mol%num = num
         mol%xyz = xyz
         mol%id = [(i, i=1, nat)]

         ! Set element symbols from atomic numbers
         do i = 1, nat
            mol%sym(i) = to_symbol(num(i))
         end do

         ! Initialize lattice and periodic (non-periodic system)
         if (allocated(mol%lattice)) deallocate (mol%lattice)
         if (allocated(mol%periodic)) deallocate (mol%periodic)
         allocate (mol%lattice(3, 3), mol%periodic(3))
         mol%lattice = 0.0_wp
         mol%periodic = .false.

         ! Write structure before attempting cavity construction
         ! (if it crashes, we have the problematic structure)
         call write_test_xyz(mol, itest, test_xyz_path)

         ! Printout - extract just the filename from the full path
         i = index(test_xyz_path, '/', back=.true.)
         write (*, '(2x, A, I0, A, I0, A, A, A)') "Testing fuzz structure #", itest, " with ", nat, &
            " atoms (", trim(test_xyz_path(i + 1:)), ")"

         ! Try to construct cavity
         allocate (cavity)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=k, blend_3b=gamma)
            call new_cavity_drop(cavity, nleb=NUM_LEB, &
                                debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
                                lsf_model=svdw_template, error=cavity_error)
         end block
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) call test_failed(error, cavity_error%message)

         ! Check if any grid points failed to converge
         if (allocated(cavity%converged)) then
            if (any(.not. cavity%converged)) then
               ! At least one projection failed - write convergence failure file
               call write_convergence_failure_xyz(mol, itest)
               write (*, '(4x, A, I0, A, I0, A)') "⚠ Warning: ", count(.not. cavity%converged), &
                  " of ", size(cavity%converged), " projections failed to converge"
            end if
         end if

         ! If we get here without crashing, deallocate and remove the xyz file
         if (allocated(cavity)) deallocate (cavity)

         ! Clean up for next iteration
         deallocate (xyz, num, radii)
         if (allocated(mol%num)) deallocate (mol%num)
         if (allocated(mol%xyz)) deallocate (mol%xyz)
         if (allocated(mol%id)) deallocate (mol%id)

      end do

      write (*, '(2x, A, I0, A)') "Fuzz testing completed: ", NUM_TESTS, " structures tested successfully"

   contains

      !> Write XYZ file before attempting cavity construction
      subroutine write_test_xyz(mol, test_num, filepath_out)
         type(structure_type), intent(in) :: mol
         integer, intent(in) :: test_num
         character(len=*), intent(out) :: filepath_out
         character(len=20) :: timestamp
         integer :: dt(8)
         type(mctc_error), allocatable :: write_error

         ! Get timestamp
         call date_and_time(values=dt)
         write (timestamp, '(I4.4,I2.2,I2.2,A,I2.2,I2.2,I2.2)') &
            dt(1), dt(2), dt(3), '_', dt(5), dt(6), dt(7)

         write (filepath_out, '(A,A,I4.4,A,A,A)') &
            trim(issues_dir)//'/', &
            'failed_', test_num, '_', trim(timestamp), '.xyz'

         call write_structure(mol, trim(filepath_out), write_error, filetype%xyz)

      end subroutine write_test_xyz

      !> Write XYZ file for structures with convergence failures
      subroutine write_convergence_failure_xyz(mol, test_num)
         type(structure_type), intent(in) :: mol
         integer, intent(in) :: test_num
         character(len=300) :: filepath
         character(len=20) :: timestamp
         integer :: dt(8)
         type(mctc_error), allocatable :: write_error

         ! Get timestamp
         call date_and_time(values=dt)
         write (timestamp, '(I4.4,I2.2,I2.2,A,I2.2,I2.2,I2.2)') &
            dt(1), dt(2), dt(3), '_', dt(5), dt(6), dt(7)

         write (filepath, '(A,A,I4.4,A,A,A)') &
            trim(issues_dir)//'/', &
            'convergence_', test_num, '_', trim(timestamp), '.xyz'

         call write_structure(mol, trim(filepath), write_error, filetype%xyz)

      end subroutine write_convergence_failure_xyz

   end subroutine test_robustness_fuzz

end module test_cavity_drop_robustness
