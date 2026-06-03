!> Regression tests for DROP cavity area/volume against fixed MC references.
module test_cavity_drop_integration_ref
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: default_cpcm_radii

   implicit none
   private

   public :: collect_cavity_drop_integration_ref

   !> One independent integration reference case
   type :: integration_case_type
      !> Blending parameter k.
      real(wp) :: blend_k
      !> Blending parameter beta.
      real(wp) :: blend_2b
      !> Blending parameter gamma.
      real(wp) :: blend_3b
      !> Dataset name used with mstore.
      character(len=12) :: dataset
      !> Structure identifier inside dataset.
      character(len=7) :: structure
      !> Reference marching-cubes area.
      real(wp) :: mc_area
      !> Reference marching-cubes volume.
      real(wp) :: mc_volume
   end type integration_case_type

   integer, parameter :: NUM_LEB = 194
   real(wp), parameter :: PROJ_TOL = 1.0e-12_wp
   integer, parameter :: PROJ_MAXITER = 1000
   integer, parameter :: PROJ_LEVEL = 1

   real(wp), parameter :: AREA_REL_THR = 1.0e-2_wp
   real(wp), parameter :: VOLUME_REL_THR = 1.0e-2_wp

   ! Reference are generated from MC using test/dev/test_cavity_drop_deflation_comparison.f90
   type(integration_case_type), parameter :: cases(*) = [ &
                           integration_case_type(1.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "Ar     ", 221.742887_wp, 310.437355_wp), &
                           integration_case_type(2.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "Ar     ", 221.742887_wp, 310.437355_wp), &
                           integration_case_type(2.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "Ar     ", 221.742887_wp, 310.437355_wp), &
                          integration_case_type(10.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "Ar     ", 221.742887_wp, 310.437355_wp), &
                           integration_case_type(1.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "O2     ", 242.589143_wp, 353.723353_wp), &
                           integration_case_type(2.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "O2     ", 196.161273_wp, 255.138423_wp), &
                           integration_case_type(2.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "O2     ", 196.161273_wp, 255.138423_wp), &
                           integration_case_type(3.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "O2     ", 185.924268_wp, 233.928477_wp), &
                           integration_case_type(3.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "O2     ", 185.924268_wp, 233.928477_wp), &
                          integration_case_type(10.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "O2     ", 179.419066_wp, 218.353773_wp), &
                           integration_case_type(1.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "CH4    ", 396.929665_wp, 742.412003_wp), &
                           integration_case_type(2.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "CH4    ", 255.697218_wp, 380.888334_wp), &
                           integration_case_type(2.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "CH4    ", 270.055161_wp, 414.747388_wp), &
                           integration_case_type(3.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "CH4    ", 229.987373_wp, 323.622720_wp), &
                          integration_case_type(10.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "CH4    ", 201.379772_wp, 261.188771_wp), &
                          integration_case_type(2.0_wp, 1.0_wp, 0.0_wp, "Amino20x4   ", "THR_xab", 871.709782_wp, 2046.069766_wp), &
                          integration_case_type(3.0_wp, 1.0_wp, 0.0_wp, "Amino20x4   ", "THR_xab", 802.552998_wp, 1673.692695_wp), &
                          integration_case_type(3.0_wp, 1.0_wp, 1.0_wp, "Amino20x4   ", "THR_xab", 816.233266_wp, 1801.648358_wp), &
                          integration_case_type(5.0_wp, 1.0_wp, 1.0_wp, "Amino20x4   ", "THR_xab", 778.058686_wp, 1516.976987_wp), &
                         integration_case_type(10.0_wp, 1.0_wp, 0.0_wp, "Amino20x4   ", "THR_xab", 782.802167_wp, 1410.364379_wp), &
                         integration_case_type(1.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "16     ", 1033.798426_wp, 3088.935546_wp), &
                          integration_case_type(2.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "16     ", 739.232911_wp, 1734.806473_wp), &
                          integration_case_type(2.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "16     ", 775.484041_wp, 1921.288922_wp), &
                          integration_case_type(3.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "16     ", 703.771421_wp, 1567.182120_wp), &
                          integration_case_type(5.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "16     ", 690.487755_wp, 1425.662212_wp), &
                          integration_case_type(2.0_wp, 1.0_wp, 0.0_wp, "But14diol   ", "30     ", 577.749734_wp, 1207.837574_wp), &
                          integration_case_type(2.0_wp, 1.0_wp, 1.0_wp, "But14diol   ", "30     ", 634.177253_wp, 1431.003011_wp), &
                           integration_case_type(3.0_wp, 1.0_wp, 0.0_wp, "But14diol   ", "30     ", 519.505143_wp, 980.963503_wp), &
                          integration_case_type(3.0_wp, 1.0_wp, 1.0_wp, "But14diol   ", "30     ", 533.169514_wp, 1051.688680_wp), &
                           integration_case_type(5.0_wp, 1.0_wp, 1.0_wp, "But14diol   ", "30     ", 496.399532_wp, 882.073907_wp), &
                          integration_case_type(10.0_wp, 1.0_wp, 0.0_wp, "But14diol   ", "30     ", 493.931152_wp, 816.379242_wp), &
                           integration_case_type(2.0_wp, 0.0_wp, 1.0_wp, "MB16-43     ", "CH4    ", 243.767818_wp, 354.652628_wp), &
                          integration_case_type(3.0_wp, 0.0_wp, 1.0_wp, "Amino20x4   ", "THR_xab", 789.581576_wp, 1694.721566_wp), &
                           integration_case_type(3.0_wp, 0.0_wp, 1.0_wp, "But14diol   ", "30     ", 512.762446_wp, 981.787682_wp), &
                          integration_case_type(3.0_wp, 0.0_wp, 1.0_wp, "MB16-43     ", "16     ", 690.329407_wp, 1491.585502_wp), &
                            integration_case_type(2.0_wp, 0.0_wp, 1.0_wp, "MB16-43     ", "O2     ", 184.603535_wp, 231.274026_wp) &
                                             ]

contains

   !> Collect all regular integration reference tests.
   subroutine collect_cavity_drop_integration_ref(testsuite)
      !> Collection of unit tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      allocate (testsuite(size(cases)))
      testsuite(1) = new_unittest(case_to_string(cases(1)), test_case_001)
      testsuite(2) = new_unittest(case_to_string(cases(2)), test_case_002)
      testsuite(3) = new_unittest(case_to_string(cases(3)), test_case_003)
      testsuite(4) = new_unittest(case_to_string(cases(4)), test_case_004)
      testsuite(5) = new_unittest(case_to_string(cases(5)), test_case_005)
      testsuite(6) = new_unittest(case_to_string(cases(6)), test_case_006)
      testsuite(7) = new_unittest(case_to_string(cases(7)), test_case_007)
      testsuite(8) = new_unittest(case_to_string(cases(8)), test_case_008)
      testsuite(9) = new_unittest(case_to_string(cases(9)), test_case_009)
      testsuite(10) = new_unittest(case_to_string(cases(10)), test_case_010)
      testsuite(11) = new_unittest(case_to_string(cases(11)), test_case_011)
      testsuite(12) = new_unittest(case_to_string(cases(12)), test_case_012)
      testsuite(13) = new_unittest(case_to_string(cases(13)), test_case_013)
      testsuite(14) = new_unittest(case_to_string(cases(14)), test_case_014)
      testsuite(15) = new_unittest(case_to_string(cases(15)), test_case_015)
      testsuite(16) = new_unittest(case_to_string(cases(16)), test_case_016)
      testsuite(17) = new_unittest(case_to_string(cases(17)), test_case_017)
      testsuite(18) = new_unittest(case_to_string(cases(18)), test_case_018)
      testsuite(19) = new_unittest(case_to_string(cases(19)), test_case_019)
      testsuite(20) = new_unittest(case_to_string(cases(20)), test_case_020)
      testsuite(21) = new_unittest(case_to_string(cases(21)), test_case_021)
      testsuite(22) = new_unittest(case_to_string(cases(22)), test_case_022)
      testsuite(23) = new_unittest(case_to_string(cases(23)), test_case_023)
      testsuite(24) = new_unittest(case_to_string(cases(24)), test_case_024)
      testsuite(25) = new_unittest(case_to_string(cases(25)), test_case_025)
      testsuite(26) = new_unittest(case_to_string(cases(26)), test_case_026)
      testsuite(27) = new_unittest(case_to_string(cases(27)), test_case_027)
      testsuite(28) = new_unittest(case_to_string(cases(28)), test_case_028)
      testsuite(29) = new_unittest(case_to_string(cases(29)), test_case_029)
      testsuite(30) = new_unittest(case_to_string(cases(30)), test_case_030)
      testsuite(31) = new_unittest(case_to_string(cases(31)), test_case_031)
      testsuite(32) = new_unittest(case_to_string(cases(32)), test_case_032)
      testsuite(33) = new_unittest(case_to_string(cases(33)), test_case_033)
      testsuite(34) = new_unittest(case_to_string(cases(34)), test_case_034)
      testsuite(35) = new_unittest(case_to_string(cases(35)), test_case_035)
      testsuite(36) = new_unittest(case_to_string(cases(36)), test_case_036)
   end subroutine collect_cavity_drop_integration_ref

   !> Validate one DROP cavity area/volume case against fixed MC references.
   !> @param[out] error     Test failure state.
   !> @param[in]  case_idx  Index into the `cases` table.
   subroutine run_single_case(error, case_idx)
      !> Test failure state.
      type(error_type), allocatable, intent(out) :: error
      !> Index into the `cases` table.
      integer, intent(in) :: case_idx

      !> Molecular structure for current case.
      type(structure_type) :: mol
      !> DROP cavity instance.
      type(cavity_type_drop), allocatable :: cavity
      !> Error from cavity routines.
      type(mctc_error), allocatable :: cavity_error
      !> Computed total cavity area.
      real(wp) :: cavity_area
      !> Computed total cavity volume.
      real(wp) :: cavity_volume
      !> Computed area ratio (cavity/reference).
      real(wp) :: area_ratio
      !> Computed volume ratio (cavity/reference).
      real(wp) :: volume_ratio

      call load_structure(cases(case_idx)%dataset, cases(case_idx)%structure, mol)

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=cases(case_idx)%blend_k, &
                                blend_2b=cases(case_idx)%blend_2b, &
                                blend_3b=cases(case_idx)%blend_3b)
         call new_cavity_drop(cavity, nleb=NUM_LEB, tolerance=PROJ_TOL, &
                              proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL, &
                              debug=.false., verbose=0, &
                              radius_model=default_cpcm_radii(), &
                              lsf_model=svdw_template, error=cavity_error)
      end block
      if (allocated(cavity_error)) then
         call test_failed(error, "new_cavity_drop failed for "//trim(case_to_string(cases(case_idx)))// &
                          ": "//trim(cavity_error%message))
         return
      end if

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "cavity%update failed for "//trim(case_to_string(cases(case_idx)))// &
                          ": "//trim(cavity_error%message))
         return
      end if

      cavity_area = cavity%total_area
      cavity_volume = cavity%total_volume

      area_ratio = cavity_area/cases(case_idx)%mc_area
      volume_ratio = cavity_volume/cases(case_idx)%mc_volume

      ! ! print out all areas and volumes, differences and ratios
      ! write(*, "(A10, A8, F12.6, A8, F12.6, A8, F12.6, A8, F12.6)") "Area", "DROP:", cavity_area, &
      !       " MC:", cases(case_idx)%mc_area, " diff:", cavity_area - cases(case_idx)%mc_area, &
      !       " ratio:", area_ratio
      ! write(*, "(A10, A8, F12.6, A8, F12.6, A8, F12.6, A8, F12.6)") "Volume", "DROP:", cavity_volume, &
      !       " MC:", cases(case_idx)%mc_volume, " diff:", cavity_volume - cases(case_idx)%mc_volume, &
      !       " ratio:", volume_ratio

      call check(error, area_ratio, 1.0_wp, thr=AREA_REL_THR, &
                 message="Area ratio mismatch for "//trim(case_to_string(cases(case_idx))))
      if (allocated(error)) return

      call check(error, volume_ratio, 1.0_wp, thr=VOLUME_REL_THR, &
                 message="Volume ratio mismatch for "//trim(case_to_string(cases(case_idx))))
      if (allocated(error)) return

   end subroutine run_single_case

   subroutine test_case_001(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 1)
   end subroutine test_case_001

   subroutine test_case_002(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 2)
   end subroutine test_case_002

   subroutine test_case_003(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 3)
   end subroutine test_case_003

   subroutine test_case_004(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 4)
   end subroutine test_case_004

   subroutine test_case_005(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 5)
   end subroutine test_case_005

   subroutine test_case_006(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 6)
   end subroutine test_case_006

   subroutine test_case_007(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 7)
   end subroutine test_case_007

   subroutine test_case_008(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 8)
   end subroutine test_case_008

   subroutine test_case_009(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 9)
   end subroutine test_case_009

   subroutine test_case_010(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 10)
   end subroutine test_case_010

   subroutine test_case_011(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 11)
   end subroutine test_case_011

   subroutine test_case_012(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 12)
   end subroutine test_case_012

   subroutine test_case_013(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 13)
   end subroutine test_case_013

   subroutine test_case_014(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 14)
   end subroutine test_case_014

   subroutine test_case_015(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 15)
   end subroutine test_case_015

   subroutine test_case_016(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 16)
   end subroutine test_case_016

   subroutine test_case_017(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 17)
   end subroutine test_case_017

   subroutine test_case_018(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 18)
   end subroutine test_case_018

   subroutine test_case_019(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 19)
   end subroutine test_case_019

   subroutine test_case_020(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 20)
   end subroutine test_case_020

   subroutine test_case_021(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 21)
   end subroutine test_case_021

   subroutine test_case_022(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 22)
   end subroutine test_case_022

   subroutine test_case_023(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 23)
   end subroutine test_case_023

   subroutine test_case_024(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 24)
   end subroutine test_case_024

   subroutine test_case_025(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 25)
   end subroutine test_case_025

   subroutine test_case_026(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 26)
   end subroutine test_case_026

   subroutine test_case_027(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 27)
   end subroutine test_case_027

   subroutine test_case_028(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 28)
   end subroutine test_case_028

   subroutine test_case_029(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 29)
   end subroutine test_case_029

   subroutine test_case_030(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 30)
   end subroutine test_case_030

   subroutine test_case_031(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 31)
   end subroutine test_case_031

   subroutine test_case_032(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 32)
   end subroutine test_case_032

   subroutine test_case_033(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 33)
   end subroutine test_case_033

   subroutine test_case_034(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 34)
   end subroutine test_case_034

   subroutine test_case_035(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 35)
   end subroutine test_case_035

   subroutine test_case_036(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 36)
   end subroutine test_case_036

   subroutine test_case_037(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 37)
   end subroutine test_case_037

   subroutine test_case_038(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 38)
   end subroutine test_case_038

   !> Load one structure for an integration reference case.
   !> @param[in]  dataset   Dataset name in mstore.
   !> @param[in]  structure Structure ID in mstore.
   !> @param[out] mol       Loaded molecular structure.
   subroutine load_structure(dataset, structure, mol)
      !> Dataset name in mstore.
      character(len=*), intent(in) :: dataset
      !> Structure ID in mstore.
      character(len=*), intent(in) :: structure
      !> Loaded molecular structure.
      type(structure_type), intent(out) :: mol

      if (trim(structure) == "Ar") then
         call new(mol, [18], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      else
         call get_structure(mol, trim(dataset), trim(structure))
      end if
   end subroutine load_structure

   !> Convert one integration case to a compact label.
   !> @param[in] c    Reference case entry.
   !> @return    str  Printable case label without trailing blanks.
   pure function case_to_string(c) result(str)
      !> Reference case entry.
      type(integration_case_type), intent(in) :: c
      !> Printable case label without trailing blanks.
      character(len=:), allocatable :: str
      !> Formatted k value.
      character(len=8) :: k_str
      !> Formatted beta value.
      character(len=8) :: b_str
      !> Formatted gamma value.
      character(len=8) :: g_str

      write (k_str, '(F4.1)') c%blend_k
      write (b_str, '(F4.1)') c%blend_2b
      write (g_str, '(F4.1)') c%blend_3b
      str = trim(c%dataset)//" "//trim(c%structure)//" k="// &
            k_str(1:4)//" b="//b_str(1:4)//" g="//g_str(1:4)
   end function case_to_string

end module test_cavity_drop_integration_ref
