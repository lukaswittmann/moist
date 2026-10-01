!> DROP cavity area and volume regression against marching-cubes references
!>
!> Shared structures and DROP settings for both level-set kinds:
!> - SvdW: `blend_k`, `blend_2b`, `blend_3b`; marching-cubes area and volume
!> - CFC: `a1`, `a2`, `c`; kernel fixes `m = 4`
!>   Marching-cubes cavity via `new_cavity_marchingcubes`, default spacing,
!>   unscreened CFC level set, and `default_cpcm_radii`; matching shape parameters
!>
!> One test per (kind, structure) pair
module test_cavity_drop_integration
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_cfc_param, only: moist_cavity_drop_lsf_cfc_param_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: default_cpcm_radii
   use moist_context, only: moist_context_type, new_context
   implicit none(type, external)
   private

   public :: collect_cavity_drop_integration

   !> Independent integration reference case
   !>
   !> Unused shape parameters default to 0
   type :: integration_case_type
      !> Level-set kind, `"svdw"` or `"cfc"`
      character(len=4) :: lsf_kind
      !> mstore dataset
      character(len=12) :: dataset
      !> Structure ID in dataset
      character(len=7) :: structure
      !> SvdW blending parameter k
      real(wp) :: blend_k = 0.0_wp
      !> SvdW blending parameter beta
      real(wp) :: blend_2b = 0.0_wp
      !> SvdW blending parameter gamma
      real(wp) :: blend_3b = 0.0_wp
      !> CFC atomic-term exponent a1
      real(wp) :: a1 = 0.0_wp
      !> CFC pair-term exponent a2
      real(wp) :: a2 = 0.0_wp
      !> CFC pair-term coupling c
      real(wp) :: c = 0.0_wp
      !> Reference marching-cubes area
      real(wp) :: mc_area
      !> Reference marching-cubes volume
      real(wp) :: mc_volume
   end type integration_case_type

   !> Lebedev points per atomic sphere
   integer, parameter :: NUM_LEB = 194
   !> DROP convergence tolerance
   real(wp), parameter :: PROJ_TOL = 1.0e-12_wp
   !> Projection iteration cap
   integer, parameter :: PROJ_MAXITER = 1000
   !> Projection level, 1 = SLSQP only
   integer, parameter :: PROJ_LEVEL = 2

   !> DROP/reference area-ratio tolerance
   real(wp), parameter :: AREA_REL_THR = 1.0e-2_wp
   !> DROP/reference volume-ratio tolerance
   real(wp), parameter :: VOLUME_REL_THR = 1.0e-2_wp

   !> Reference cases
   type(integration_case_type), parameter :: cases(61) = [ &
      integration_case_type("svdw", "MB16-43", "Ar", blend_k=1.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=221.742887_wp, mc_volume=310.437355_wp), &
      integration_case_type("svdw", "MB16-43", "Ar", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=221.742887_wp, mc_volume=310.437355_wp), &
      integration_case_type("svdw", "MB16-43", "Ar", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=221.742887_wp, mc_volume=310.437355_wp), &
      integration_case_type("svdw", "MB16-43", "Ar", blend_k=10.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=221.742887_wp, mc_volume=310.437355_wp), &
      integration_case_type("svdw", "MB16-43", "O2", blend_k=1.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=242.589143_wp, mc_volume=353.723353_wp), &
      integration_case_type("svdw", "MB16-43", "O2", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=196.161273_wp, mc_volume=255.138423_wp), &
      integration_case_type("svdw", "MB16-43", "O2", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=196.161273_wp, mc_volume=255.138423_wp), &
      integration_case_type("svdw", "MB16-43", "O2", blend_k=3.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=185.924268_wp, mc_volume=233.928477_wp), &
      integration_case_type("svdw", "MB16-43", "O2", blend_k=3.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=185.924268_wp, mc_volume=233.928477_wp), &
      integration_case_type("svdw", "MB16-43", "O2", blend_k=10.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=179.419066_wp, mc_volume=218.353773_wp), &
      integration_case_type("svdw", "MB16-43", "CH4", blend_k=1.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=396.929665_wp, mc_volume=742.412003_wp), &
      integration_case_type("svdw", "MB16-43", "CH4", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=255.697218_wp, mc_volume=380.888334_wp), &
      integration_case_type("svdw", "MB16-43", "CH4", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=270.055161_wp, mc_volume=414.747388_wp), &
      integration_case_type("svdw", "MB16-43", "CH4", blend_k=3.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=229.987373_wp, mc_volume=323.622720_wp), &
      integration_case_type("svdw", "MB16-43", "CH4", blend_k=10.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=201.379772_wp, mc_volume=261.188771_wp), &
      integration_case_type("svdw", "Amino20x4", "THR_xab", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=871.709782_wp, mc_volume=2046.069766_wp), &
      integration_case_type("svdw", "Amino20x4", "THR_xab", blend_k=3.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=802.552998_wp, mc_volume=1673.692695_wp), &
      integration_case_type("svdw", "Amino20x4", "THR_xab", blend_k=3.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=816.233266_wp, mc_volume=1801.648358_wp), &
      integration_case_type("svdw", "Amino20x4", "THR_xab", blend_k=5.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=778.058686_wp, mc_volume=1516.976987_wp), &
      integration_case_type("svdw", "Amino20x4", "THR_xab", blend_k=10.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=782.802167_wp, mc_volume=1410.364379_wp), &
      integration_case_type("svdw", "MB16-43", "16", blend_k=1.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=1033.798426_wp, mc_volume=3088.935546_wp), &
      integration_case_type("svdw", "MB16-43", "16", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=739.232911_wp, mc_volume=1734.806473_wp), &
      integration_case_type("svdw", "MB16-43", "16", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=775.484041_wp, mc_volume=1921.288922_wp), &
      integration_case_type("svdw", "MB16-43", "16", blend_k=3.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=703.771421_wp, mc_volume=1567.182120_wp), &
      integration_case_type("svdw", "MB16-43", "16", blend_k=5.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=690.487755_wp, mc_volume=1425.662212_wp), &
      integration_case_type("svdw", "But14diol", "30", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=577.749734_wp, mc_volume=1207.837574_wp), &
      integration_case_type("svdw", "But14diol", "30", blend_k=2.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=634.177253_wp, mc_volume=1431.003011_wp), &
      integration_case_type("svdw", "But14diol", "30", blend_k=3.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=519.505143_wp, mc_volume=980.963503_wp), &
      integration_case_type("svdw", "But14diol", "30", blend_k=3.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=533.169514_wp, mc_volume=1051.688680_wp), &
      integration_case_type("svdw", "But14diol", "30", blend_k=5.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, &
                            mc_area=496.399532_wp, mc_volume=882.073907_wp), &
      integration_case_type("svdw", "But14diol", "30", blend_k=10.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp, &
                            mc_area=493.931152_wp, mc_volume=816.379242_wp), &
      integration_case_type("svdw", "MB16-43", "CH4", blend_k=2.0_wp, blend_2b=0.0_wp, blend_3b=1.0_wp, &
                            mc_area=243.767818_wp, mc_volume=354.652628_wp), &
      integration_case_type("svdw", "Amino20x4", "THR_xab", blend_k=3.0_wp, blend_2b=0.0_wp, blend_3b=1.0_wp, &
                            mc_area=789.581576_wp, mc_volume=1694.721566_wp), &
      integration_case_type("svdw", "But14diol", "30", blend_k=3.0_wp, blend_2b=0.0_wp, blend_3b=1.0_wp, &
                            mc_area=512.762446_wp, mc_volume=981.787682_wp), &
      integration_case_type("svdw", "MB16-43", "16", blend_k=3.0_wp, blend_2b=0.0_wp, blend_3b=1.0_wp, &
                            mc_area=690.329407_wp, mc_volume=1491.585502_wp), &
      integration_case_type("svdw", "MB16-43", "O2", blend_k=2.0_wp, blend_2b=0.0_wp, blend_3b=1.0_wp, &
                            mc_area=184.603535_wp, mc_volume=231.274026_wp), &
      integration_case_type("cfc", "MB16-43", "Ar", a1=-15.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=221.644539_wp, mc_volume=310.219352_wp), &
      integration_case_type("cfc", "MB16-43", "O2", a1=-15.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=179.585114_wp, mc_volume=219.289158_wp), &
      integration_case_type("cfc", "MB16-43", "CH4", a1=-15.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=202.648865_wp, mc_volume=263.879456_wp), &
      integration_case_type("cfc", "Amino20x4", "THR_xab", a1=-15.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=761.936610_wp, mc_volume=1468.405797_wp), &
      integration_case_type("cfc", "MB16-43", "16", a1=-15.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=683.018739_wp, mc_volume=1408.227683_wp), &
      integration_case_type("cfc", "But14diol", "30", a1=-15.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=486.990697_wp, mc_volume=840.481784_wp), &
      integration_case_type("cfc", "MB16-43", "O2", a1=-10.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=180.166379_wp, mc_volume=221.925529_wp), &
      integration_case_type("cfc", "MB16-43", "CH4", a1=-10.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=207.288656_wp, mc_volume=273.573900_wp), &
      integration_case_type("cfc", "Amino20x4", "THR_xab", a1=-10.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=768.428592_wp, mc_volume=1512.722523_wp), &
      integration_case_type("cfc", "MB16-43", "16", a1=-10.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=684.257010_wp, mc_volume=1435.514879_wp), &
      integration_case_type("cfc", "But14diol", "30", a1=-10.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=492.761391_wp, mc_volume=869.601869_wp), &
      integration_case_type("cfc", "MB16-43", "O2", a1=-20.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=179.307289_wp, mc_volume=218.178104_wp), &
      integration_case_type("cfc", "MB16-43", "CH4", a1=-20.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=200.956384_wp, mc_volume=260.175413_wp), &
      integration_case_type("cfc", "But14diol", "30", a1=-20.0_wp, a2=-9.0_wp, c=5.0_wp, &
                            mc_area=485.343691_wp, mc_volume=829.473708_wp), &
      integration_case_type("cfc", "MB16-43", "O2", a1=-15.0_wp, a2=-6.0_wp, c=5.0_wp, &
                            mc_area=179.668448_wp, mc_volume=219.430137_wp), &
      integration_case_type("cfc", "MB16-43", "CH4", a1=-15.0_wp, a2=-6.0_wp, c=5.0_wp, &
                            mc_area=203.102076_wp, mc_volume=264.778704_wp), &
      integration_case_type("cfc", "MB16-43", "16", a1=-15.0_wp, a2=-6.0_wp, c=5.0_wp, &
                            mc_area=681.158101_wp, mc_volume=1438.670151_wp), &
      integration_case_type("cfc", "MB16-43", "CH4", a1=-15.0_wp, a2=-12.0_wp, c=5.0_wp, &
                            mc_area=202.636397_wp, mc_volume=263.806572_wp), &
      integration_case_type("cfc", "Amino20x4", "THR_xab", a1=-15.0_wp, a2=-12.0_wp, c=5.0_wp, &
                            mc_area=767.613503_wp, mc_volume=1449.274241_wp), &
      integration_case_type("cfc", "MB16-43", "O2", a1=-15.0_wp, a2=-9.0_wp, c=2.5_wp, &
                            mc_area=179.575342_wp, mc_volume=219.270717_wp), &
      integration_case_type("cfc", "MB16-43", "CH4", a1=-15.0_wp, a2=-9.0_wp, c=2.5_wp, &
                            mc_area=202.637462_wp, mc_volume=263.855981_wp), &
      integration_case_type("cfc", "But14diol", "30", a1=-15.0_wp, a2=-9.0_wp, c=2.5_wp, &
                            mc_area=488.094604_wp, mc_volume=836.122085_wp), &
      integration_case_type("cfc", "MB16-43", "CH4", a1=-15.0_wp, a2=-9.0_wp, c=10.0_wp, &
                            mc_area=202.671697_wp, mc_volume=263.926430_wp), &
      integration_case_type("cfc", "Amino20x4", "THR_xab", a1=-15.0_wp, a2=-9.0_wp, c=10.0_wp, &
                            mc_area=759.486532_wp, mc_volume=1483.876765_wp), &
      integration_case_type("cfc", "MB16-43", "16", a1=-15.0_wp, a2=-9.0_wp, c=10.0_wp, &
                            mc_area=680.931308_wp, mc_volume=1418.494143_wp) &
      ]

contains

   !> Collect one test per (level-set kind, structure)
   !>
   !> @param[out] testsuite  Collection of unit tests
   subroutine collect_cavity_drop_integration(testsuite)
      !> Collection of unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("svdw_Ar", test_svdw_ar), &
                  new_unittest("svdw_O2", test_svdw_o2), &
                  new_unittest("svdw_CH4", test_svdw_ch4), &
                  new_unittest("svdw_THR_xab", test_svdw_thr_xab), &
                  new_unittest("svdw_16", test_svdw_16), &
                  new_unittest("svdw_30", test_svdw_30), &
                  new_unittest("cfc_Ar", test_cfc_ar), &
                  new_unittest("cfc_O2", test_cfc_o2), &
                  new_unittest("cfc_CH4", test_cfc_ch4), &
                  new_unittest("cfc_THR_xab", test_cfc_thr_xab), &
                  new_unittest("cfc_16", test_cfc_16), &
                  new_unittest("cfc_30", test_cfc_30) &
                  ]
   end subroutine collect_cavity_drop_integration

   !> SvdW-DROP cases of argon
   subroutine test_svdw_ar(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "svdw", "Ar")
   end subroutine test_svdw_ar

   !> SvdW-DROP cases of O2
   subroutine test_svdw_o2(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "svdw", "O2")
   end subroutine test_svdw_o2

   !> SvdW-DROP cases of CH4
   subroutine test_svdw_ch4(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "svdw", "CH4")
   end subroutine test_svdw_ch4

   !> SvdW-DROP cases of Amino20x4 THR_xab
   subroutine test_svdw_thr_xab(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "svdw", "THR_xab")
   end subroutine test_svdw_thr_xab

   !> SvdW-DROP cases of MB16-43 16
   subroutine test_svdw_16(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "svdw", "16")
   end subroutine test_svdw_16

   !> SvdW-DROP cases of But14diol 30
   subroutine test_svdw_30(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "svdw", "30")
   end subroutine test_svdw_30

   !> CFC-DROP cases of argon
   subroutine test_cfc_ar(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "cfc", "Ar")
   end subroutine test_cfc_ar

   !> CFC-DROP cases of O2
   subroutine test_cfc_o2(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "cfc", "O2")
   end subroutine test_cfc_o2

   !> CFC-DROP cases of CH4
   subroutine test_cfc_ch4(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "cfc", "CH4")
   end subroutine test_cfc_ch4

   !> CFC-DROP cases of Amino20x4 THR_xab
   subroutine test_cfc_thr_xab(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "cfc", "THR_xab")
   end subroutine test_cfc_thr_xab

   !> CFC-DROP cases of MB16-43 16
   subroutine test_cfc_16(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "cfc", "16")
   end subroutine test_cfc_16

   !> CFC-DROP cases of But14diol 30
   subroutine test_cfc_30(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      call run_group(error, "cfc", "30")
   end subroutine test_cfc_30

   !> Run one group until its first failure
   !>
   !> @param[out] error      Test failure state; names the failing case
   !> @param[in]  lsf_kind   `"svdw"` or `"cfc"`
   !> @param[in]  structure  Structure identifier as in the `cases` table
   subroutine run_group(error, lsf_kind, structure)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error
      !> Level-set kind of the group
      character(len=*), intent(in) :: lsf_kind
      !> Structure identifier of the group
      character(len=*), intent(in) :: structure

      !> Case index
      integer :: ic
      !> Number of cases run
      integer :: nrun

      nrun = 0
      do ic = 1, size(cases)
         if (cases(ic)%lsf_kind /= lsf_kind .or. cases(ic)%structure /= structure) cycle
         call run_single_case(error, cases(ic))
         if (allocated(error)) return
         nrun = nrun + 1
      end do

      ! An empty group would pass silently
      if (nrun == 0) then
         call test_failed(error, "No "//lsf_kind//" case for structure "//structure)
      end if
   end subroutine run_group

   !> Check one DROP area and volume case against MC references
   !>
   !> @param[in]  c      Reference case entry
   subroutine run_single_case(error, c)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error
      !> Reference case entry
      type(integration_case_type), intent(in) :: c

      !> Molecular structure for current case
      type(structure_type) :: mol
      !> Level set of the case's kind
      class(moist_cavity_drop_lsf_type), allocatable :: lsf_template
      !> DROP cavity instance
      type(cavity_type_drop), allocatable :: cavity
      !> Error from cavity routines
      type(mctc_error), allocatable :: cavity_error
      !> Computed area ratio (cavity/reference)
      real(wp) :: area_ratio
      !> Computed volume ratio (cavity/reference)
      real(wp) :: volume_ratio
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, verbosity=0)

      call load_structure(c%dataset, c%structure, mol)

      select case (c%lsf_kind)
      case ("svdw")
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=c%blend_k, &
               blend_2b=c%blend_2b, blend_3b=c%blend_3b))
            allocate (lsf_template, source=svdw_template)
         end block
      case ("cfc")
         block
            type(moist_cavity_drop_lsf_cfc_type) :: cfc_template
            ! `m` keeps its default
            call cfc_template%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=c%a1, a2=c%a2, c=c%c))
            allocate (lsf_template, source=cfc_template)
         end block
      case default
         call test_failed(error, "Unknown level-set kind for "//case_to_string(c))
         return
      end select

      allocate (cavity)
      call new_cavity_drop(cavity, ctx, radius_model=default_cpcm_radii(), lsf_model=lsf_template, &
         error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=NUM_LEB, tolerance=PROJ_TOL, &
         proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL))
      if (allocated(cavity_error)) then
         call test_failed(error, "new_cavity_drop failed for "//case_to_string(c)// &
                          ": "//trim(cavity_error%message))
         return
      end if

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "cavity%update failed for "//case_to_string(c)// &
                          ": "//trim(cavity_error%message))
         return
      end if

      area_ratio = cavity%total_area/c%mc_area
      volume_ratio = cavity%total_volume/c%mc_volume

      call check(error, area_ratio, 1.0_wp, thr=AREA_REL_THR, &
                 message="Area ratio mismatch for "//case_to_string(c))
      if (allocated(error)) return

      call check(error, volume_ratio, 1.0_wp, thr=VOLUME_REL_THR, &
                 message="Volume ratio mismatch for "//case_to_string(c))
      if (allocated(error)) return

   end subroutine run_single_case

   !> Load an integration case structure
   !>
   !> @param[in]  dataset   Dataset name in mstore
   !> @param[in]  structure Structure ID in mstore
   !> @param[out] mol       Loaded molecular structure
   subroutine load_structure(dataset, structure, mol)
      !> Dataset name in mstore
      character(len=*), intent(in) :: dataset
      !> Structure ID in mstore
      character(len=*), intent(in) :: structure
      !> Loaded molecular structure
      type(structure_type), intent(out) :: mol

      if (trim(structure) == "Ar") then
         call new(mol, [18], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      else
         call get_structure(mol, trim(dataset), trim(structure))
      end if
   end subroutine load_structure

   !> Unique compact label for an integration case
   !>
   !> @param[in] c  Reference case entry
   pure function case_to_string(c) result(str)
      !> Reference case entry
      type(integration_case_type), intent(in) :: c
      !> Printable case label: kind, dataset, structure, shape parameters
      character(len=:), allocatable :: str
      !> Formatted first shape parameter (k or a1)
      character(len=8) :: p1_str
      !> Formatted second shape parameter (beta or a2)
      character(len=8) :: p2_str
      !> Formatted third shape parameter (gamma or c)
      character(len=8) :: p3_str

      select case (c%lsf_kind)
      case ("svdw")
         write (p1_str, "(F4.1)") c%blend_k
         write (p2_str, "(F4.1)") c%blend_2b
         write (p3_str, "(F4.1)") c%blend_3b
         str = trim(c%lsf_kind)//" "//trim(c%dataset)//" "//trim(c%structure)//" k="// &
               p1_str(1:4)//" b="//p2_str(1:4)//" g="//p3_str(1:4)
      case ("cfc")
         write (p1_str, "(F5.1)") c%a1
         write (p2_str, "(F5.1)") c%a2
         write (p3_str, "(F4.1)") c%c
         str = trim(c%lsf_kind)//" "//trim(c%dataset)//" "//trim(c%structure)//" a1="//trim(adjustl(p1_str))// &
               " a2="//trim(adjustl(p2_str))//" c="//trim(adjustl(p3_str))
      case default
         str = trim(c%lsf_kind)//" "//trim(c%dataset)//" "//trim(c%structure)
      end select
   end function case_to_string

end module test_cavity_drop_integration
