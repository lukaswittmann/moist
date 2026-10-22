!> Energies of the PCM components against analytic results and other codes
!>
!> - `born`: ion of either sign at the center of a spherical cavity, Gaussian and
!>   point surface, against the analytic energy of every model
!> - `reference`: point charges in spheres against exact series, other codes and
!>   the literature; rows without a reference value skip and print the moist
!>   energy (see `reference_table`)
!> - `switched_off_point`: a point switched off on the surface drops out of
!>   the energy of every model, so a cavity may drop it without a jump
!> - `conductor_limit`: eps = +Infinity on Gaussian and point spheres against the
!>   analytic f = 1 limit and a very large finite eps; IEF-PCM and COSMO equal
!>   CPCM there on a molecular surface
!> - CPCM identities exact at any discretization: rigid translation, linearity
!>   in f(eps)
!> - Generated Gaussian IEF-PCM and SS(V)PE kernels against PySCF on an
!>   identical surface
module test_model_component_pcm_references
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_is_finite
   use mctc_env, only: wp
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use mstore, only: get_structure
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed, skip_test
   use test_helpers, only: check_moist_error, component_view, get_test_structures, &
      & center_at_origin, get_test_cavity_iswig, build_test_cavity, fill_point_charge_potential, &
      & stage_point_charge_energy, submit
   use test_model_component_pcm_fixtures, only: pcm_cpcm, pcm_cosmo, pcm_iefpcm, pcm_ssvpe, &
      & pcm_dpcm, pcm_dpcm_original, pcm_lpb, pcm_slic, pcm_nonlocal, pcm_physics, new_pcm, &
      & seed_drop, first_seed, host_energy, stage_energy, fd_cavity, get_seed, nmol, mol_label, &
      & get_solute, host_charges, thr_vacuous
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_type, only: cavity_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, &
      & moist_cavity_iswig_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: radius_type, radius_type_static, new_radii_custom_atoms
   use moist_channels_coupling, only: coupling_type, gaussian_potential_request_type, &
      & point_potential_request_type, nuclear_charge_request_type
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, &
      & moist_pcm_parameters_type, solver_type
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, &
      & new_component_cpcm
   use moist_model_continuum_component_pcm_solvers, only: pcm_dense_solver
   use pcmgen_gauss_ssvpe, only: ssvpe_state, ssvpe_assemble, ssvpe_medium, ssvpe_solve
   use pcmgen_gauss_iefpcm, only: iefpcm_state, iefpcm_assemble, iefpcm_medium, iefpcm_solve
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_references

   !> Tolerance of identities that hold to roundoff
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)

   !> Born sphere radius, bohr
   real(wp), parameter :: born_radius = 4.0_wp
   !> Lebedev order of the Born sphere
   integer, parameter :: born_nleb = 302
   !> Nonlocal correlation length of the Born check, bohr
   !>
   !> nu = sqrt(eps/eps_inf) R/lambda ~ 1, where the nonlocal energy lies far
   !> from both Born limits and Lambda = 3.8 bohr is resolved by the grid
   real(wp), parameter :: born_lambda = 25.0_wp

   !> Very large finite epsilon compared with the conductor limit
   real(wp), parameter :: eps_large = 1.0e12_wp
   !> Relative agreement of eps = +Infinity with `eps_large`
   real(wp), parameter :: thr_large = 1.0e-10_wp
   !> Relative agreement of IEF-PCM and COSMO with CPCM at eps = +Infinity, measured 3e-16
   real(wp), parameter :: thr_conductor_identity = 1.0e-13_wp

   !> Factors that switch a surface point off in `check_switched_off_point`
   !>
   !> The decay is checked on the last two, where the leading order dominates: a
   !> point surface has terms in s^(1/2) and s that can cancel near s = 1e-4
   real(wp), parameter :: switch_off(3) = [1.0e-4_wp, 1.0e-8_wp, 1.0e-12_wp]
   !> Relative tolerances of the strong-screening checks, twice the measured deviations
   !>
   !> LPB at kappa = 40 and 100: 4.7e-5 and 2.6e-5 (Debye-Hueckel); nonlocal at
   !> lambda = 1e-2 and 1e-4 bohr: 2.0e-4 and 2.3e-6 (dense and iterative alike)
   real(wp), parameter :: lpb_strong_tolerance = 1.0e-4_wp, nonlocal_strong_tolerance = 5.0e-4_wp
   !> Roundoff floor of the energy gap of a switched-off point, Eh
   real(wp), parameter :: thr_gap = 1.0e-12_wp

   !> Reference systems of the smoke test: point charges in spheres of explicit radii
   integer, parameter :: ref_sphere = 1, ref_water = 2, ref_ion_pair = 3, ref_ammonium = 4, &
      & ref_guanidinium = 5, ref_acetate = 6, ref_sulfate = 7, ref_glycine = 8, ref_benzene = 9
   !> Number of reference systems
   integer, parameter :: nref_system = 9
   !> Reference system names for messages
   character(len=*), parameter :: ref_system_label(nref_system) = [character(len=24) :: &
      & "off-center pair, sphere", "TIP3P water", "Na+ Cl- contact pair", "ammonium (+1)", &
      & "guanidinium (+1)", "acetate (-1)", "sulfate (-2)", "glycine zwitterion", "benzene"]
   !> Lebedev order of the moist side of the smoke test
   integer, parameter :: ref_nleb = 302
   !> Sentinel of a reference energy not filled in yet
   real(wp), parameter :: ref_pending = huge(1.0_wp)

   !> Row of the reference table: one model on one reference system
   type :: reference_row
      !> Model identifier
      integer :: kind
      !> Reference system
      integer :: system
      !> Reference polarization energy, hartree; `ref_pending` until filled in
      real(wp) :: energy = ref_pending
      !> The reference is the exact series of `sphere_reference` instead of `energy`
      logical :: exact = .false.
      !> Relative tolerance, the expected discretization difference of the two codes
      real(wp) :: rtol = 1.0e-2_wp
      !> Origin of the reference: code and version, or citation
      character(len=160) :: source = ""
   end type reference_row

   !> PySCF reference of the generated layer kernels
   !>
   !> PySCF `solvent/pcm.py` (SWIG surface) driven with a point-charge potential,
   !> both family members on one surface: plain names SS(V)PE, `_ief_` names
   !> IEF-PCM. S, D, A, R and the potential are shared, only K, q, q_sym and the
   !> energy differ. Generated by `moist_dev/tools/pcm_layer_ref/gen_golden.py`,
   !> case hf (H 0 0 0; F 0 0 1.75)
   real(wp), parameter :: golden_eps = 78.3553_wp
   integer, parameter :: hf_n = 27
   real(wp), parameter :: hf_xyz(3, 27) = reshape([ &
      & 2.494438484425882_wp, 0.0_wp, 0.0_wp, &
      & -2.494438484425882_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 2.494438484425882_wp, 0.0_wp, &
      & 0.0_wp, -2.494438484425882_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, -2.494438484425882_wp, &
      & 1.440164730460245_wp, 1.440164730460245_wp, 1.440164730460245_wp, &
      & -1.440164730460245_wp, 1.440164730460245_wp, 1.440164730460245_wp, &
      & 1.440164730460245_wp, -1.440164730460245_wp, 1.440164730460245_wp, &
      & -1.440164730460245_wp, -1.440164730460245_wp, 1.440164730460245_wp, &
      & 1.440164730460245_wp, 1.440164730460245_wp, -1.440164730460245_wp, &
      & -1.440164730460245_wp, 1.440164730460245_wp, -1.440164730460245_wp, &
      & 1.440164730460245_wp, -1.440164730460245_wp, -1.440164730460245_wp, &
      & -1.440164730460245_wp, -1.440164730460245_wp, -1.440164730460245_wp, &
      & 3.333476883732769_wp, 0.0_wp, 1.75_wp, &
      & -3.333476883732769_wp, 0.0_wp, 1.75_wp, &
      & 0.0_wp, 3.333476883732769_wp, 1.75_wp, &
      & 0.0_wp, -3.333476883732769_wp, 1.75_wp, &
      & 0.0_wp, 0.0_wp, 5.083476883732769_wp, &
      & 0.0_wp, 0.0_wp, -1.5834768837327688_wp, &
      & 1.924583776160509_wp, 1.924583776160509_wp, 3.674583776160509_wp, &
      & -1.924583776160509_wp, 1.924583776160509_wp, 3.674583776160509_wp, &
      & 1.924583776160509_wp, -1.924583776160509_wp, 3.674583776160509_wp, &
      & -1.924583776160509_wp, -1.924583776160509_wp, 3.674583776160509_wp, &
      & 1.924583776160509_wp, 1.924583776160509_wp, -0.1745837761605089_wp, &
      & -1.924583776160509_wp, 1.924583776160509_wp, -0.1745837761605089_wp, &
      & 1.924583776160509_wp, -1.924583776160509_wp, -0.1745837761605089_wp, &
      & -1.924583776160509_wp, -1.924583776160509_wp, -0.1745837761605089_wp &
      & ], [3, 27])
   real(wp), parameter :: hf_normal(3, 27) = reshape([ &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & -1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, &
      & 0.0_wp, -1.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, -1.0_wp, &
      & 0.5773502691896257_wp, 0.5773502691896257_wp, 0.5773502691896257_wp, &
      & -0.5773502691896257_wp, 0.5773502691896257_wp, 0.5773502691896257_wp, &
      & 0.5773502691896257_wp, -0.5773502691896257_wp, 0.5773502691896257_wp, &
      & -0.5773502691896257_wp, -0.5773502691896257_wp, 0.5773502691896257_wp, &
      & 0.5773502691896257_wp, 0.5773502691896257_wp, -0.5773502691896257_wp, &
      & -0.5773502691896257_wp, 0.5773502691896257_wp, -0.5773502691896257_wp, &
      & 0.5773502691896257_wp, -0.5773502691896257_wp, -0.5773502691896257_wp, &
      & -0.5773502691896257_wp, -0.5773502691896257_wp, -0.5773502691896257_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & -1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, &
      & 0.0_wp, -1.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.0_wp, &
      & 0.0_wp, 0.0_wp, -1.0_wp, &
      & 0.5773502691896257_wp, 0.5773502691896257_wp, 0.5773502691896257_wp, &
      & -0.5773502691896257_wp, 0.5773502691896257_wp, 0.5773502691896257_wp, &
      & 0.5773502691896257_wp, -0.5773502691896257_wp, 0.5773502691896257_wp, &
      & -0.5773502691896257_wp, -0.5773502691896257_wp, 0.5773502691896257_wp, &
      & 0.5773502691896257_wp, 0.5773502691896257_wp, -0.5773502691896257_wp, &
      & -0.5773502691896257_wp, 0.5773502691896257_wp, -0.5773502691896257_wp, &
      & 0.5773502691896257_wp, -0.5773502691896257_wp, -0.5773502691896257_wp, &
      & -0.5773502691896257_wp, -0.5773502691896257_wp, -0.5773502691896257_wp &
      & ], [3, 27])
   real(wp), parameter :: hf_xi(*) = [ &
      & 2.130658696182589_wp, 2.130658696182589_wp, 2.130658696182589_wp, &
      & 2.130658696182589_wp, 2.130658696182589_wp, 2.0088042832863953_wp, &
      & 2.0088042832863953_wp, 2.0088042832863953_wp, 2.0088042832863953_wp, &
      & 2.0088042832863953_wp, 2.0088042832863953_wp, 2.0088042832863953_wp, &
      & 2.0088042832863953_wp, 1.594370452925747_wp, 1.594370452925747_wp, &
      & 1.594370452925747_wp, 1.594370452925747_wp, 1.594370452925747_wp, &
      & 1.594370452925747_wp, 1.5031868786496838_wp, 1.5031868786496838_wp, &
      & 1.5031868786496838_wp, 1.5031868786496838_wp, 1.5031868786496838_wp, &
      & 1.5031868786496838_wp, 1.5031868786496838_wp, 1.5031868786496838_wp &
      & ]
   real(wp), parameter :: hf_f(*) = [ &
      & 0.3742539186401483_wp, 0.3742539186401483_wp, 0.3742539186401483_wp, &
      & 0.3742539186401483_wp, 0.9351770959980755_wp, 0.02031572973675191_wp, &
      & 0.02031572973675191_wp, 0.02031572973675191_wp, 0.02031572973675191_wp, &
      & 0.7701381620437818_wp, 0.7701381620437818_wp, 0.7701381620437818_wp, &
      & 0.7701381620437818_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp, &
      & 0.02800972227932615_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 0.702062547712545_wp, &
      & 0.702062547712545_wp, 0.702062547712545_wp, 0.702062547712545_wp &
      & ]
   real(wp), parameter :: hf_area(*) = [ &
      & 1.9508800058775093_wp, 1.9508800058775093_wp, 1.9508800058775093_wp, &
      & 1.9508800058775093_wp, 4.874814150687485_wp, 0.11913768325734597_wp, &
      & 0.11913768325734597_wp, 0.11913768325734597_wp, 0.11913768325734597_wp, &
      & 4.516326885761973_wp, 4.516326885761973_wp, 4.516326885761973_wp, &
      & 4.516326885761973_wp, 9.309224431242598_wp, 9.309224431242598_wp, &
      & 9.309224431242598_wp, 9.309224431242598_wp, 9.309224431242598_wp, &
      & 0.26074879095502307_wp, 10.472877485147922_wp, 10.472877485147922_wp, &
      & 10.472877485147922_wp, 10.472877485147922_wp, 7.352615049104301_wp, &
      & 7.352615049104301_wp, 7.352615049104301_wp, 7.352615049104301_wp &
      & ]
   real(wp), parameter :: hf_curv(*) = [ &
      & 0.40089182645454546_wp, 0.40089182645454546_wp, 0.40089182645454546_wp, &
      & 0.40089182645454546_wp, 0.40089182645454546_wp, 0.40089182645454546_wp, &
      & 0.40089182645454546_wp, 0.40089182645454546_wp, 0.40089182645454546_wp, &
      & 0.40089182645454546_wp, 0.40089182645454546_wp, 0.40089182645454546_wp, &
      & 0.40089182645454546_wp, 0.2999870810204082_wp, 0.2999870810204082_wp, &
      & 0.2999870810204082_wp, 0.2999870810204082_wp, 0.2999870810204082_wp, &
      & 0.2999870810204082_wp, 0.2999870810204082_wp, 0.2999870810204082_wp, &
      & 0.2999870810204082_wp, 0.2999870810204082_wp, 0.2999870810204082_wp, &
      & 0.2999870810204082_wp, 0.2999870810204082_wp, 0.2999870810204082_wp &
      & ]
   real(wp), parameter :: hf_v(*) = [ &
      & 0.025448238614779267_wp, 0.025448238614779267_wp, 0.025448238614779267_wp, &
      & 0.025448238614779267_wp, 0.05785129048386302_wp, -0.02957982826240152_wp, &
      & -0.02957982826240152_wp, -0.02957982826240152_wp, -0.02957982826240152_wp, &
      & 0.04783887615082673_wp, 0.04783887615082673_wp, 0.04783887615082673_wp, &
      & 0.04783887615082673_wp, -0.012031837027736608_wp, -0.012031837027736608_wp, &
      & -0.012031837027736608_wp, -0.012031837027736608_wp, -0.036144963639547265_wp, &
      & 0.11603711365325416_wp, -0.028456144301747574_wp, -0.028456144301747574_wp, &
      & -0.028456144301747574_wp, -0.028456144301747574_wp, 0.023333469405338828_wp, &
      & 0.023333469405338828_wp, 0.023333469405338828_wp, 0.023333469405338828_wp &
      & ]
   real(wp), parameter :: hf_q(*) = [ &
      & -0.0030222591944487794_wp, -0.0030222591944487824_wp, -0.0030222591944487807_wp, &
      & -0.0030222591944487794_wp, -0.020014917889904178_wp, 0.0004153630575754697_wp, &
      & 0.0004153630575754696_wp, 0.00041536305757546965_wp, 0.00041536305757546965_wp, &
      & -0.014399179288649411_wp, -0.01439917928864941_wp, -0.01439917928864941_wp, &
      & -0.014399179288649418_wp, 0.009836184967798187_wp, 0.009836184967798187_wp, &
      & 0.009836184967798188_wp, 0.009836184967798185_wp, 0.017904752026660992_wp, &
      & -0.0023264576889079932_wp, 0.017596419368429813_wp, 0.017596419368429813_wp, &
      & 0.017596419368429816_wp, 0.01759641936842981_wp, -0.008255641652360476_wp, &
      & -0.008255641652360476_wp, -0.008255641652360484_wp, -0.008255641652360474_wp &
      & ]
   real(wp), parameter :: hf_qsym(*) = [ &
      & -0.0032327372952371605_wp, -0.0032327372952371623_wp, -0.003232737295237161_wp, &
      & -0.0032327372952371605_wp, -0.020821410844149832_wp, 0.00035344644856540584_wp, &
      & 0.0003534464485654059_wp, 0.00035344644856540584_wp, 0.00035344644856540584_wp, &
      & -0.014990579633768306_wp, -0.014990579633768304_wp, -0.0149905796337683_wp, &
      & -0.01499057963376831_wp, 0.009553971001864694_wp, 0.009553971001864696_wp, &
      & 0.009553971001864696_wp, 0.009553971001864697_wp, 0.017958207139327394_wp, &
      & -0.0017051226551002593_wp, 0.017756326628328063_wp, 0.01775632662832806_wp, &
      & 0.017756326628328066_wp, 0.01775632662832806_wp, -0.007094466255850573_wp, &
      & -0.007094466255850574_wp, -0.007094466255850578_wp, -0.007094466255850572_wp &
      & ]
   real(wp), parameter :: hf_energy = -0.0042169943992978725_wp
   real(wp), parameter :: hf_ief_q(*) = [ &
      & -0.0032738612301247533_wp, -0.0032738612301247594_wp, -0.003273861230124754_wp, &
      & -0.003273861230124756_wp, -0.020699921252084356_wp, 0.0003427729916246641_wp, &
      & 0.00034277299162466387_wp, 0.00034277299162466403_wp, 0.00034277299162466414_wp, &
      & -0.015001116086919793_wp, -0.015001116086919791_wp, -0.015001116086919794_wp, &
      & -0.015001116086919793_wp, 0.009462005663253314_wp, 0.009462005663253325_wp, &
      & 0.00946200566325332_wp, 0.009462005663253318_wp, 0.01794860966088233_wp, &
      & -0.0016234647836306064_wp, 0.017735249811301163_wp, 0.017735249811301163_wp, &
      & 0.017735249811301167_wp, 0.017735249811301163_wp, -0.007057594108837933_wp, &
      & -0.007057594108837931_wp, -0.007057594108837933_wp, -0.007057594108837932_wp &
      & ]
   real(wp), parameter :: hf_ief_qsym(*) = [ &
      & -0.0032782894188899004_wp, -0.003278289418889903_wp, -0.003278289418889901_wp, &
      & -0.0032782894188899026_wp, -0.020716605599694737_wp, 0.00034146899742661026_wp, &
      & 0.0003414689974266102_wp, 0.0003414689974266103_wp, 0.0003414689974266103_wp, &
      & -0.01501294277751896_wp, -0.015012942777518962_wp, -0.015012942777518964_wp, &
      & -0.015012942777518964_wp, 0.009456721728217085_wp, 0.009456721728217085_wp, &
      & 0.009456721728217087_wp, 0.009456721728217087_wp, 0.01794921911411844_wp, &
      & -0.0016106734867473806_wp, 0.017738162363294348_wp, 0.017738162363294348_wp, &
      & 0.017738162363294348_wp, 0.017738162363294348_wp, -0.007034668180173885_wp, &
      & -0.007034668180173882_wp, -0.00703466818017388_wp, -0.007034668180173883_wp &
      & ]
   real(wp), parameter :: hf_ief_energy = -0.004205905577815418_wp

contains

   !> Collect the reference tests, three per model
   !>
   !> @param[out] testsuite Collection of tests
   subroutine collect_model_component_pcm_references(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("cpcm_born", test_cpcm_born), &
         & new_unittest("cpcm_reference", test_cpcm_reference), &
         & new_unittest("cpcm_switched_off_point", test_cpcm_switched_off_point), &
         & new_unittest("cosmo_born", test_cosmo_born), &
         & new_unittest("cosmo_reference", test_cosmo_reference), &
         & new_unittest("cosmo_switched_off_point", test_cosmo_switched_off_point), &
         & new_unittest("iefpcm_born", test_iefpcm_born), &
         & new_unittest("iefpcm_reference", test_iefpcm_reference), &
         & new_unittest("iefpcm_switched_off_point", test_iefpcm_switched_off_point), &
         & new_unittest("ssvpe_born", test_ssvpe_born), &
         & new_unittest("ssvpe_reference", test_ssvpe_reference), &
         & new_unittest("ssvpe_switched_off_point", test_ssvpe_switched_off_point), &
         & new_unittest("dpcm_born", test_dpcm_born), &
         & new_unittest("dpcm_reference", test_dpcm_reference), &
         & new_unittest("dpcm_switched_off_point", test_dpcm_switched_off_point), &
         & new_unittest("dpcm_original_born", test_dpcm_original_born), &
         & new_unittest("dpcm_original_reference", test_dpcm_original_reference), &
         & new_unittest("dpcm_original_switched_off_point", test_dpcm_original_switched_off_point), &
         & new_unittest("lpb_born", test_lpb_born), &
         & new_unittest("lpb_reference", test_lpb_reference), &
         & new_unittest("lpb_switched_off_point", test_lpb_switched_off_point), &
         & new_unittest("slic_born", test_slic_born), &
         & new_unittest("slic_reference", test_slic_reference), &
         & new_unittest("slic_switched_off_point", test_slic_switched_off_point), &
         & new_unittest("slic_static_potential", test_slic_static_potential), &
         & new_unittest("nonlocal_born", test_nonlocal_born), &
         & new_unittest("nonlocal_reference", test_nonlocal_reference), &
         & new_unittest("nonlocal_switched_off_point", test_nonlocal_switched_off_point), &
         & new_unittest("cpcm_conductor_limit", test_cpcm_conductor_limit), &
         & new_unittest("cosmo_conductor_limit", test_cosmo_conductor_limit), &
         & new_unittest("iefpcm_conductor_limit", test_iefpcm_conductor_limit), &
         & new_unittest("ssvpe_conductor_limit", test_ssvpe_conductor_limit), &
         & new_unittest("dpcm_conductor_limit", test_dpcm_conductor_limit), &
         & new_unittest("dpcm_original_conductor_limit", test_dpcm_original_conductor_limit), &
         & new_unittest("lpb_conductor_limit", test_lpb_conductor_limit), &
         & new_unittest("slic_conductor_limit", test_slic_conductor_limit), &
         & new_unittest("iefpcm_equals_cpcm_conductor", test_iefpcm_equals_cpcm_conductor), &
         & new_unittest("cosmo_equals_cpcm_conductor", test_cosmo_equals_cpcm_conductor), &
         & new_unittest("nonlocal_interior_permittivity", test_nonlocal_interior_permittivity), &
         & new_unittest("nonlocal_local_limit_lambda", test_nonlocal_local_limit_lambda), &
         & new_unittest("nonlocal_escaped_charge", test_nonlocal_escaped_charge), &
         & new_unittest("lpb_strong_screening", test_lpb_strong_screening), &
         & new_unittest("nonlocal_strong_screening", test_nonlocal_strong_screening), &
         & new_unittest("cpcm_translation_invariance", test_cpcm_translation_invariance), &
         & new_unittest("cpcm_dielectric_scaling", test_cpcm_dielectric_scaling), &
         & new_unittest("layer_kernels_vs_pyscf", test_layer_goldens) &
         & ]
   end subroutine collect_model_component_pcm_references

   !* ================================================================================= *!
   !*                                 Per-model entries                                 *!
   !* ================================================================================= *!

   !> CPCM: Born ion
   subroutine test_cpcm_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_cpcm)
   end subroutine test_cpcm_born

   !> CPCM: reference energies of other codes and the literature
   subroutine test_cpcm_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_cpcm)
   end subroutine test_cpcm_reference

   !> CPCM: a switched-off point drops out of the energy
   subroutine test_cpcm_switched_off_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_switched_off_point(error, pcm_cpcm)
   end subroutine test_cpcm_switched_off_point

   !> COSMO: Born ion
   subroutine test_cosmo_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_cosmo)
   end subroutine test_cosmo_born

   !> COSMO: reference energies of other codes and the literature
   subroutine test_cosmo_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_cosmo)
   end subroutine test_cosmo_reference

   !> COSMO: a switched-off point drops out of the energy
   subroutine test_cosmo_switched_off_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_switched_off_point(error, pcm_cosmo)
   end subroutine test_cosmo_switched_off_point

   !> IEF-PCM: Born ion
   subroutine test_iefpcm_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_iefpcm)
   end subroutine test_iefpcm_born

   !> IEF-PCM: reference energies of other codes and the literature
   subroutine test_iefpcm_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_iefpcm)
   end subroutine test_iefpcm_reference

   !> IEF-PCM: a switched-off point drops out of the energy
   subroutine test_iefpcm_switched_off_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_switched_off_point(error, pcm_iefpcm)
   end subroutine test_iefpcm_switched_off_point

   !> SS(V)PE: Born ion
   subroutine test_ssvpe_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_ssvpe)
   end subroutine test_ssvpe_born

   !> SS(V)PE: reference energies of other codes and the literature
   subroutine test_ssvpe_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_ssvpe)
   end subroutine test_ssvpe_reference

   !> SS(V)PE: a switched-off point drops out of the energy
   subroutine test_ssvpe_switched_off_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_switched_off_point(error, pcm_ssvpe)
   end subroutine test_ssvpe_switched_off_point

   !> DPCM: Born ion
   subroutine test_dpcm_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_dpcm)
   end subroutine test_dpcm_born

   !> DPCM: reference energies of other codes and the literature
   subroutine test_dpcm_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_dpcm)
   end subroutine test_dpcm_reference

   !> DPCM: a switched-off point drops out of the energy
   subroutine test_dpcm_switched_off_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_switched_off_point(error, pcm_dpcm)
   end subroutine test_dpcm_switched_off_point

   !> Original DPCM: Born ion
   subroutine test_dpcm_original_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_born

   !> Original DPCM: reference energies of other codes and the literature
   subroutine test_dpcm_original_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_reference

   !> DPCM original: a switched-off point drops out of the energy
   subroutine test_dpcm_original_switched_off_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_switched_off_point(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_switched_off_point

   !> LPB: Born ion
   subroutine test_lpb_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_lpb)
   end subroutine test_lpb_born

   !> LPB: reference energies of other codes and the literature
   subroutine test_lpb_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_lpb)
   end subroutine test_lpb_reference

   !> LPB: a switched-off point drops out of the energy
   subroutine test_lpb_switched_off_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_switched_off_point(error, pcm_lpb)
   end subroutine test_lpb_switched_off_point

   !> SLIC: Born ion
   subroutine test_slic_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_slic)
   end subroutine test_slic_born

   !> SLIC: reference energies of other codes and the literature
   subroutine test_slic_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_slic)
   end subroutine test_slic_reference

   !> SLIC: the static term adds phi_static times the net solute charge
   !>
   !> Born ion of either sign; the static term is a constant, so raising
   !> phi_static by 0.1 Eh/e shifts the energy by exactly 0.1*Q
   subroutine test_slic_static_potential(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the cavity and the components
      type(moist_context_type), target :: ctx
      !> Components without and with the static potential
      class(model_continuum_component_pcm), allocatable :: plain, static
      !> Single-sphere solute
      type(structure_type) :: mol
      !> Spherical cavity
      class(cavity_type), allocatable :: cavity
      !> Physical parameters of the shifted component
      type(pcm_physics) :: physics
      !> Static potential step, hartree per elementary charge
      real(wp), parameter :: dphi_static = 0.1_wp
      !> Charge and energies without and with the static potential
      real(wp) :: charge, e0, e1
      !> Sign index
      integer :: isign

      call new_context(ctx, verbosity=0)
      call new_pcm(pcm_slic, ctx, plain, error)
      if (allocated(error)) return
      physics%slic_phi_static = dphi_static
      call new_pcm(pcm_slic, ctx, static, error, physics)
      if (allocated(error)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_sphere_cavity(pcm_slic, [born_radius], born_nleb, ctx, cavity, error)
      if (allocated(error)) return
      call static%configure_cavity(cavity, err)
      call check_moist_error(error, err, "cavity configuration")
      if (allocated(error)) return
      call cavity%update(mol, error=err)
      call check_moist_error(error, err, "cavity update")
      if (allocated(error)) return

      do isign = 1, 2
         charge = merge(1.0_wp, -1.0_wp, isign == 1)
         mol%charge = charge
         call host_energy(plain, mol, cavity, mol, [charge], e0, error)
         if (allocated(error)) return
         call host_energy(static, mol, cavity, mol, [charge], e1, error)
         if (allocated(error)) return
         call check(error, e1 - e0, dphi_static*charge, thr=1.0e-12_wp, &
            & message="SLIC static term is not phi_static times the net charge")
         if (allocated(error)) return
      end do
   end subroutine test_slic_static_potential

   !> SLIC: a switched-off point drops out of the energy
   subroutine test_slic_switched_off_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_switched_off_point(error, pcm_slic)
   end subroutine test_slic_switched_off_point

   !> Nonlocal: Born ion
   subroutine test_nonlocal_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_nonlocal)
   end subroutine test_nonlocal_born

   !> Nonlocal: reference energies of other codes and the literature
   subroutine test_nonlocal_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_nonlocal)
   end subroutine test_nonlocal_reference

   !> Nonlocal: a switched-off point drops out of the energy
   subroutine test_nonlocal_switched_off_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_switched_off_point(error, pcm_nonlocal)
   end subroutine test_nonlocal_switched_off_point

   !> CPCM: conductor limit
   subroutine test_cpcm_conductor_limit(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_conductor_limit(error, pcm_cpcm)
   end subroutine test_cpcm_conductor_limit

   !> COSMO: conductor limit
   subroutine test_cosmo_conductor_limit(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_conductor_limit(error, pcm_cosmo)
   end subroutine test_cosmo_conductor_limit

   !> IEFPCM: conductor limit
   subroutine test_iefpcm_conductor_limit(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_conductor_limit(error, pcm_iefpcm)
   end subroutine test_iefpcm_conductor_limit

   !> SSVPE: conductor limit
   subroutine test_ssvpe_conductor_limit(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_conductor_limit(error, pcm_ssvpe)
   end subroutine test_ssvpe_conductor_limit

   !> DPCM: conductor limit
   subroutine test_dpcm_conductor_limit(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_conductor_limit(error, pcm_dpcm)
   end subroutine test_dpcm_conductor_limit

   !> DPCM_ORIGINAL: conductor limit
   subroutine test_dpcm_original_conductor_limit(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_conductor_limit(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_conductor_limit

   !> LPB: conductor limit
   subroutine test_lpb_conductor_limit(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_conductor_limit(error, pcm_lpb)
   end subroutine test_lpb_conductor_limit

   !> SLIC: conductor limit
   subroutine test_slic_conductor_limit(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_conductor_limit(error, pcm_slic)
   end subroutine test_slic_conductor_limit

   !> IEF-PCM equals CPCM at eps = +Infinity
   subroutine test_iefpcm_equals_cpcm_conductor(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_conductor_identity(error, pcm_iefpcm)
   end subroutine test_iefpcm_equals_cpcm_conductor

   !> COSMO equals CPCM at eps = +Infinity
   subroutine test_cosmo_equals_cpcm_conductor(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_conductor_identity(error, pcm_cosmo)
   end subroutine test_cosmo_equals_cpcm_conductor

   !* ================================================================================= *!
   !*                                   CPCM identities                                 *!
   !* ================================================================================= *!

!> The CPCM energy is invariant under rigid translation of the solute
!>
!> Every iSwiG grid point is placed relative to an atomic center, so shifting the
!> whole structure shifts the surface with it and leaves every distance in both
!> the matrix and the potential trace unchanged. Unlike rotation -- which
!> reorients the molecule against a lab-fixed Lebedev grid and therefore only
!> holds up to discretization error -- this is an exact identity, and it fails
!> the moment any part of the assembly picks up an absolute coordinate
   subroutine test_cpcm_translation_invariance(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: shifted
      type(model_continuum_component_cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      real(wp), allocatable :: qat(:)

      !> Dielectric constant
      real(wp), parameter :: epsilon = 78.4_wp
      !> Lebedev grid size
      integer, parameter :: nleb = 26
      !> Rigid shift applied to the whole structure, bohr
      real(wp), parameter :: shift(3) = [1.7_wp, -2.3_wp, 0.9_wp]

      integer :: imol, iat
      real(wp) :: energy, energy_shifted
      character(len=64) :: context

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      call get_test_structures(mols, 5)

      do imol = 1, size(mols)
         call center_at_origin(mols(imol))
         call make_neutral_charges(mols(imol)%nat, qat)

         call get_test_cavity_iswig(mols(imol), cavity, err, nleb=nleb)
         if (allocated(err)) then
            call test_failed(error, "Cavity setup failed: "//err%message)
            return
         end if
         call new_component_cpcm(pcm_model, epsilon=epsilon, error=err, &
                                 param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
         if (allocated(err)) then
            call test_failed(error, "CPCM initialization failed: "//err%message)
            return
         end if
         call pcm_model%update(mols(imol), cavity, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM update failed: "//err%message)
            return
         end if
         call stage_point_charge_energy(error, pcm_model, cavity, qat, mols(imol), coupling)
         if (allocated(error)) return
         energy = 0.0_wp
         call pcm_model%get_energy(component_view(coupling), cavity, energy, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed: "//err%message)
            return
         end if

         shifted = mols(imol)
         do iat = 1, shifted%nat
            shifted%xyz(:, iat) = shifted%xyz(:, iat) + shift
         end do

         call get_test_cavity_iswig(shifted, cavity, err, nleb=nleb)
         if (allocated(err)) then
            call test_failed(error, "Shifted cavity setup failed: "//err%message)
            return
         end if
         call pcm_model%update(shifted, cavity, err)
         if (allocated(err)) then
            call test_failed(error, "Shifted CPCM update failed: "//err%message)
            return
         end if
         ! The potential trace moves with the surface: re-stage and refill it
         ! on the shifted cavity
         call pcm_model%prepare_energy(cavity, coupling, err)
         if (allocated(err)) then
            call test_failed(error, "Shifted CPCM staging failed: "//err%message)
            return
         end if
         call fill_point_charge_potential(cavity, coupling, qat, shifted)
         energy_shifted = 0.0_wp
         call pcm_model%get_energy(component_view(coupling), cavity, energy_shifted, err)
         if (allocated(err)) then
            call test_failed(error, "Shifted CPCM energy failed: "//err%message)
            return
         end if

         ! Guard against a vacuous comparison: a structure whose polarization
         ! energy happened to vanish would pass any invariance test
         write (context, "(a,i0)") "structure ", imol
         call check(error, abs(energy) > 1.0e-6_wp, &
            & more=trim(context)//": CPCM energy is zero, invariance is vacuous")
         if (allocated(error)) return
         call check(error, energy_shifted, energy, thr=1.0e-13_wp*abs(energy), &
            & more=trim(context)//": CPCM energy is not translation invariant")
         if (allocated(error)) return

         deallocate (qat)
      end do

   end subroutine test_cpcm_translation_invariance

!> The CPCM energy is exactly linear in the dielectric factor f(eps)
!>
!> CPCM solves A q = -f(eps) phi and reports E = 1/2 q.phi, so at fixed geometry
!> and fixed potential the whole dielectric dependence is the single scalar
!> f(eps) = (eps-1)/eps. Energies computed at different dielectrics must
!> therefore rescale exactly, with no discretization error of their own: the
!> matrix, its factorization and the trace are identical in all of them. This
!> pins down the one place the dielectric enters and would catch, for instance,
!> f(eps) being applied twice or to the wrong side of the solve
   subroutine test_cpcm_dielectric_scaling(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type), allocatable :: mols(:)
      type(model_continuum_component_cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      real(wp), allocatable :: qat(:)

      !> Dielectrics spanning the weakly polar to the conductor limit
      real(wp), parameter :: epsilons(4) = [2.0_wp, 8.93_wp, 78.4_wp, 1.0e6_wp]
      !> Lebedev grid size
      integer, parameter :: nleb = 26

      integer :: imol, ieps
      real(wp) :: energy, e_ref, feps, feps_ref, predicted
      character(len=96) :: context

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      call get_test_structures(mols, 5)

      do imol = 1, size(mols)
         call center_at_origin(mols(imol))
         call make_neutral_charges(mols(imol)%nat, qat)

         call get_test_cavity_iswig(mols(imol), cavity, err, nleb=nleb)
         if (allocated(err)) then
            call test_failed(error, "Cavity setup failed: "//err%message)
            return
         end if

         e_ref = 0.0_wp
         feps_ref = 0.0_wp
         do ieps = 1, size(epsilons)
            feps = (epsilons(ieps) - 1.0_wp)/epsilons(ieps)
            call new_component_cpcm(pcm_model, epsilon=epsilons(ieps), error=err, &
                                    param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
            if (allocated(err)) then
               call test_failed(error, "CPCM initialization failed: "//err%message)
               return
            end if
            call pcm_model%update(mols(imol), cavity, err)
            if (allocated(err)) then
               call test_failed(error, "CPCM update failed: "//err%message)
               return
            end if
            call stage_point_charge_energy(error, pcm_model, cavity, qat, mols(imol), coupling)
            if (allocated(error)) return
            energy = 0.0_wp
            call pcm_model%get_energy(component_view(coupling), cavity, energy, err)
            if (allocated(err)) then
               call test_failed(error, "CPCM energy failed: "//err%message)
               return
            end if

            if (ieps == 1) then
               e_ref = energy
               feps_ref = feps
               write (context, "(a,i0)") "structure ", imol
               call check(error, abs(e_ref) > 1.0e-6_wp, &
                  & more=trim(context)//": CPCM energy is zero, scaling is vacuous")
               if (allocated(error)) return
               cycle
            end if

            predicted = e_ref*feps/feps_ref
            write (context, "(a,i0,a,f0.2)") "structure ", imol, ", eps = ", epsilons(ieps)
            call check(error, energy, predicted, thr=1.0e-14_wp*abs(predicted), &
               & more=trim(context)//": CPCM energy does not scale with f(eps)")
            if (allocated(error)) return
         end do

         deallocate (qat)
      end do

   end subroutine test_cpcm_dielectric_scaling

!> Deterministic, charge-neutral atomic charges for a structure
!>
!> @param[in]  nat  Number of atoms
!> @param[out] qat  Atomic charges summing to zero
   subroutine make_neutral_charges(nat, qat)

      !> Number of atoms
      integer, intent(in) :: nat

      !> Atomic charges summing to zero
      real(wp), allocatable, intent(out) :: qat(:)

      integer :: iat

      allocate (qat(nat))
      do iat = 1, nat
         qat(iat) = 0.2_wp*sin(1.3_wp*real(iat, wp))
      end do
      qat = qat - sum(qat)/real(nat, wp)

   end subroutine make_neutral_charges

   !* ================================================================================= *!
   !*                                  PySCF layer kernels                              *!
   !* ================================================================================= *!

   !> Generated Gaussian IEF-PCM and SS(V)PE kernels reproduce PySCF charges
   subroutine test_layer_goldens(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err
      type(ssvpe_state) :: ssv
      type(iefpcm_state) :: ief
      type(pcm_dense_solver) :: solver
      real(wp) :: q(hf_n), qsym(hf_n), energy

      call ssvpe_assemble(ssv, hf_xyz, hf_normal, hf_curv, hf_curv, hf_xi, hf_f, &
         & hf_area, err, solver)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call ssvpe_medium(ssv, golden_eps, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call ssvpe_solve(ssv, hf_v, energy, q, qsym, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, maxval(abs(q - hf_q)), 0.0_wp, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(qsym - hf_qsym)), 0.0_wp, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, energy, hf_energy, thr=1.0e-10_wp)
      if (allocated(error)) return

      call iefpcm_assemble(ief, hf_xyz, hf_normal, hf_curv, hf_curv, hf_xi, hf_f, &
         & hf_area, err, solver)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call iefpcm_medium(ief, golden_eps, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call iefpcm_solve(ief, hf_v, energy, q, qsym, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, maxval(abs(q - hf_ief_q)), 0.0_wp, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(qsym - hf_ief_qsym)), 0.0_wp, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, energy, hf_ief_energy, thr=1.0e-10_wp)
   end subroutine test_layer_goldens

   !* ================================================================================= *!
   !*                                       Born ion                                    *!
   !* ================================================================================= *!

   !> Ion of either sign at the center of a spherical cavity against the analytic energy
   !>
   !> The host potential is that of the central charge, Gaussian-probed on the
   !> iSwiG surface (SLIC: SvdW-DROP, see `first_seed`); the even models give
   !> the same energy for both signs, SLIC solvates the anion more strongly for
   !> beta < 0. The Gaussian surface runs the default parameters, the same
   !> sphere stripped of its Gaussian fields (point surface, point branch of
   !> component and host) every case of `born_cases`: dielectric values,
   !> screening and correlation lengths. Finite differences cannot see a
   !> point branch with the wrong dielectric, only these references do
   !> @param[in]  kind  Model identifier
   subroutine check_born(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the cavity and the component
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Single-sphere solute
      type(structure_type) :: mol
      !> Spherical cavity
      class(cavity_type), allocatable :: cavity
      !> Physical parameters of every case
      type(pcm_physics), allocatable :: cases(:)
      !> Charge, energy, reference and relative deviation per sign
      real(wp) :: charge, energy(2), reference, deviation
      !> Surface, case and sign indices
      integer :: isurf, icase, ncase, isign
      !> Gaussian surface (else point surface)
      logical :: smeared
      !> Failure context
      character(len=200) :: context

      call new_context(ctx, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_sphere_cavity(kind, [born_radius], born_nleb, ctx, cavity, error)
      if (allocated(error)) return
      cases = born_cases(kind)

      do isurf = 1, 2
         smeared = isurf == 1
         call new_pcm(kind, ctx, pcm, error, cases(1))
         if (allocated(error)) return
         call pcm%configure_cavity(cavity, err)
         call check_moist_error(error, err, "cavity configuration")
         if (allocated(error)) return
         call cavity%update(mol, error=err)
         call check_moist_error(error, err, "cavity update")
         if (allocated(error)) return
         if (.not. smeared) then
            if (allocated(cavity%xi0)) deallocate (cavity%xi0)
            if (allocated(cavity%f)) deallocate (cavity%f)
         end if
         ncase = merge(1, size(cases), smeared)
         do icase = 1, ncase
            call new_pcm(kind, ctx, pcm, error, cases(icase))
            if (allocated(error)) return
            do isign = 1, 2
               charge = merge(1.0_wp, -1.0_wp, isign == 1)
               mol%charge = charge
               call host_energy(pcm, mol, cavity, mol, [charge], energy(isign), error)
               if (allocated(error)) return
               reference = born_reference(kind, cases(icase), charge, born_radius)
               deviation = abs(energy(isign) - reference)/abs(reference)
               write (context, '(a,", ",a," surface, eps ",f6.1,", kappa ",f4.2,", lambda ",f5.1, &
                  & ", charge ",f4.1,": energy ",es16.8," reference ",es16.8,", rel. dev ",es9.2)') &
                  & trim(pcm%name), trim(merge("Gaussian", "point   ", smeared)), cases(icase)%eps, &
                  & cases(icase)%kappa, cases(icase)%lambda, charge, energy(isign), reference, deviation
               call check(error, deviation <= born_tolerance(kind, smeared), &
                  & message="Born energy deviates from the analytic result", more=trim(context))
               if (allocated(error)) return
            end do
            if (kind == pcm_slic) then
               call check(error, energy(2) < energy(1), &
                  & message="SLIC must solvate the anion more strongly for beta < 0", &
                  & more=trim(context))
            else
               call check(error, abs(energy(2) - energy(1)) <= 1.0e-12_wp*abs(energy(1)), &
                  & message=trim(pcm%name)//": Born energy is not even in the charge", &
                  & more=trim(context))
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine check_born

   !> Parameter sets of the Born check, the first one the Gaussian default
   !>
   !> - All models: eps = 78.4 and 2
   !> - LPB: also kappa = 0 (pure dielectric limit) at both eps, and eps = 1
   !>   with kappa > 0 (screening alone)
   !> - Nonlocal: also a second correlation length at eps = 78.4
   !> @param[in] kind Model identifier
   function born_cases(kind) result(cases)
      !> Model identifier
      integer, intent(in) :: kind
      !> Parameter sets
      type(pcm_physics), allocatable :: cases(:)

      type(pcm_physics) :: base

      base%lambda = born_lambda
      cases = [base, base]
      cases(2)%eps = 2.0_wp
      select case (kind)
      case (pcm_lpb)
         cases = [cases, base, base, base]
         cases(3)%kappa = 0.0_wp
         cases(4)%kappa = 0.0_wp
         cases(4)%eps = 2.0_wp
         cases(5)%eps = 1.0_wp
      case (pcm_nonlocal)
         cases = [cases, base]
         cases(3)%lambda = 5.0_wp
      end select
   end function born_cases

   !> Energy of a central charge in the Born sphere, Gaussian or point surface
   !>
   !> The point surface is the Gaussian one without its widths: the same 302
   !> Lebedev points and areas, collocated with point charges
   !>
   !> @param[in]  kind          Model identifier
   !> @param[in]  physics       Physical parameters
   !> @param[in]  point_surface Point instead of Gaussian surface charges
   !> @param[in]  solver        Linear solver of the component
   !> @param[in]  charge        Central charge, e
   !> @param[out] energy        Polarization energy, hartree
   !> @param[out] error         Test failure
   !> @param[in]  offset        Charge position relative to the center, bohr (default 0)
   !> @param[in]  alpha         Exponent of a Gaussian density at the center instead of the
   !>                           point charge, 1/bohr^2
   subroutine central_charge_energy(kind, physics, point_surface, solver, charge, energy, error, &
      & offset, alpha)
      !> Model identifier
      integer, intent(in) :: kind
      !> Physical parameters
      type(pcm_physics), intent(in) :: physics
      !> Point instead of Gaussian surface charges
      logical, intent(in) :: point_surface
      !> Linear solver of the component
      integer, intent(in) :: solver
      !> Central charge
      real(wp), intent(in) :: charge
      !> Polarization energy
      real(wp), intent(out) :: energy
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Charge position relative to the center
      real(wp), intent(in), optional :: offset(3)
      !> Exponent of a Gaussian density instead of the point charge
      real(wp), intent(in), optional :: alpha

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the cavity and the component
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Single-sphere solute and the position of its charge
      type(structure_type) :: mol, src
      !> Spherical cavity
      class(cavity_type), allocatable :: cavity
      !> Host coupling of the Gaussian density
      type(coupling_type), target :: coupling

      energy = 0.0_wp
      call new_context(ctx, verbosity=0)
      call new_pcm(kind, ctx, pcm, error, physics)
      if (allocated(error)) return
      pcm%param%solver = solver
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      mol%charge = charge
      call new_sphere_cavity(kind, [born_radius], born_nleb, ctx, cavity, error)
      if (allocated(error)) return
      call pcm%configure_cavity(cavity, err)
      call check_moist_error(error, err, "cavity configuration")
      if (allocated(error)) return
      call cavity%update(mol, error=err)
      call check_moist_error(error, err, "cavity update")
      if (allocated(error)) return
      if (point_surface) deallocate (cavity%xi0)
      if (present(alpha)) then
         call stage_energy(pcm, mol, cavity, coupling, error)
         if (allocated(error)) return
         call answer_gaussian_density(cavity, coupling, charge, alpha)
         call pcm%get_energy(component_view(coupling), cavity, energy, err)
         call check_moist_error(error, err, "energy")
         return
      end if
      src = mol
      if (present(offset)) src%xyz(:, 1) = offset
      call host_energy(pcm, mol, cavity, src, [charge], energy, error)
   end subroutine central_charge_energy

   !> Answer the potential requests of a Gaussian density at the origin
   !>
   !> charge (alpha/pi)^(3/2) exp(-alpha r^2) has the potential charge erf(beta r)/r with
   !> beta^2 = alpha, probed by a Gaussian of width xi0 beta^2 = alpha xi0^2/(alpha + xi0^2);
   !> no nuclei
   !>
   !> @param[in]     cavity   Cavity the coupling was staged on
   !> @param[in,out] coupling Staged coupling
   !> @param[in]     charge   Charge of the density, e
   !> @param[in]     alpha    Exponent of the density, 1/bohr^2
   subroutine answer_gaussian_density(cavity, coupling, charge, alpha)
      !> Cavity the coupling was staged on
      class(cavity_type), intent(in) :: cavity
      !> Staged coupling
      type(coupling_type), intent(inout), target :: coupling
      !> Charge of the density
      real(wp), intent(in) :: charge
      !> Exponent of the density
      real(wp), intent(in) :: alpha

      !> Exponents of the probed density, potential and its gradient on the grid
      real(wp), allocatable :: beta(:), phi(:), dphi(:, :)
      !> Distance of a grid point from the center
      real(wp) :: r
      !> Grid point index
      integer :: i

      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            beta = sqrt(alpha*cavity%xi0**2/(alpha + cavity%xi0**2))
         type is (point_potential_request_type)
            allocate (beta(cavity%ngrid), source=sqrt(alpha))
         type is (nuclear_charge_request_type)
            call submit(coupling, "z", [0.0_wp])
            cycle
         class default
            cycle
         end select
         allocate (phi(cavity%ngrid), dphi(3, cavity%ngrid))
         do i = 1, cavity%ngrid
            r = norm2(cavity%xyz(:, i))
            phi(i) = charge*erf(beta(i)*r)/r
            dphi(:, i) = (2.0_wp*charge*beta(i)/sqrt(acos(-1.0_wp))*exp(-(beta(i)*r)**2) &
               & - phi(i))*cavity%xyz(:, i)/r**2
         end do
         call submit(coupling, "phi", phi)
         call submit(coupling, "dphi_dr", dphi)
         deallocate (beta, phi, dphi)
      end do
   end subroutine answer_gaussian_density

   !> Nonlocal: interior permittivity, local limit and full nonlocal energy
   !>
   !> The host supplies the potentials of the charge in vacuum; the energy is the
   !> physical reaction energy for eps_solute = eta, on Gaussian and point
   !> surfaces: the local limit eps = eps_inf = 80 and the full nonlocal energy at
   !> lambda = 25 bohr, deviations in units of Q**2/(2 R) (tolerances at `interior_tolerance`)
   subroutine test_nonlocal_interior_permittivity(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(pcm_physics) :: physics
      real(wp) :: energy, reference, deviation
      real(wp), parameter :: eta(2) = [2.0_wp, 4.0_wp]
      integer :: ieta, isurf, icase
      character(len=160) :: context

      do ieta = 1, size(eta)
         do isurf = 1, 2
            do icase = 1, 2
               physics = pcm_physics(eps_solute=eta(ieta), lambda=born_lambda)
               ! local limit eps = eps_inf: Q**2/(2 R) (1/eps - 1/eta) whatever lambda
               if (icase == 1) physics = pcm_physics(eps=80.0_wp, eps_inf=80.0_wp, lambda=3.0_wp, &
                  & eps_solute=eta(ieta))
               call central_charge_energy(pcm_nonlocal, physics, isurf == 2, solver_type%dense, &
                  & 1.0_wp, energy, error)
               if (allocated(error)) return
               reference = born_reference(pcm_nonlocal, physics, 1.0_wp, born_radius)
               ! relative to the Born scale Q**2/(2 R): at eta = 4 the full energy is
               ! 3e-3 Eh after a cancellation of 0.25 against 0.26 (in units of Q**2/R)
               deviation = abs(energy - reference)*2.0_wp*born_radius
               write (context, '("eta ",f4.1,1x,a,1x,a,", energy ",es16.8,", reference ",es16.8, &
                  & ", dev/(Q^2/2R) ",es9.2)') eta(ieta), merge("point   ", "Gaussian", isurf == 2), &
                  & merge("local limit  ", "full nonlocal", icase == 1), energy, reference, deviation
               call check(error, deviation <= interior_tolerance(isurf == 2, icase == 1), &
                  & message="Nonlocal energy deviates from the interior-permittivity result", &
                  & more=trim(context))
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_nonlocal_interior_permittivity

   !> Tolerance of the interior-permittivity check in units of Q**2/(2 R)
   !>
   !> Twice the measured deviations of the 302-point sphere at eta = 2 and 4:
   !> - Gaussian: 1.0e-8 in the local limit, 1.3e-4 full nonlocal (the
   !>   discretization error of the eta = 1 check, the same at every eta)
   !> - point: 2.9e-3 and 2.1e-3
   !>
   !> @param[in] point_surface Point instead of Gaussian surface charges
   !> @param[in] local_limit   Local limit eps = eps_inf
   pure function interior_tolerance(point_surface, local_limit) result(tol)
      !> Point instead of Gaussian surface charges
      logical, intent(in) :: point_surface
      !> Local limit eps = eps_inf
      logical, intent(in) :: local_limit
      !> Tolerance
      real(wp) :: tol

      if (point_surface) then
         tol = merge(6.0e-3_wp, 5.0e-3_wp, local_limit)
      else
         tol = merge(3.0e-8_wp, 3.0e-4_wp, local_limit)
      end if
   end function interior_tolerance

   !> Nonlocal: the local limit eps = eps_inf does not depend on lambda
   !>
   !> lambda enters only through the Yukawa double layer, which cancels exactly on
   !> the discrete surface; an off-center charge, as the central one hides part of
   !> a leak by symmetry. Gaussian and point surfaces, interior permittivity 1 and 2
   subroutine test_nonlocal_local_limit_lambda(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(pcm_physics) :: physics
      real(wp), parameter :: lambda(3) = [0.5_wp, 3.0_wp, 25.0_wp], eta(2) = [1.0_wp, 2.0_wp]
      real(wp) :: energy(size(lambda)), spread
      integer :: ieta, isurf, ilambda
      character(len=160) :: context

      do ieta = 1, size(eta)
         do isurf = 1, 2
            do ilambda = 1, size(lambda)
               physics = pcm_physics(eps=80.0_wp, eps_inf=80.0_wp, lambda=lambda(ilambda), &
                  & eps_solute=eta(ieta))
               call central_charge_energy(pcm_nonlocal, physics, isurf == 2, solver_type%dense, &
                  & 1.0_wp, energy(ilambda), error, offset=[0.0_wp, 0.0_wp, 1.5_wp])
               if (allocated(error)) return
            end do
            spread = (maxval(energy) - minval(energy))/abs(energy(1))
            write (context, '("eta ",f4.1,1x,a,", energies ",3es20.12,", rel. spread ",es9.2)') &
               & eta(ieta), merge("point   ", "Gaussian", isurf == 2), energy, spread
            call check(error, spread <= 1.0e-12_wp, &
               & message="Nonlocal local limit depends on lambda", more=trim(context))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_nonlocal_local_limit_lambda

   !> Nonlocal: solute charge outside the cavity is polarized as in IEF-PCM
   !>
   !> A unit Gaussian density exp(-0.2 r^2) at the center of the sphere, 9.4 % of it
   !> outside, in the local limit eps = eps_inf = 80. The potential right-hand side and
   !> the surface-charge energy recover its volume polarization, 1/2 int rho dphi with
   !> dphi(r) = (1/eps - 1) int Q(s)/s^2 ds from max(r, R) on (Q the enclosed charge,
   !> quadrature to 1e-13). Written with the normal field the model would see only the
   !> enclosed charge, (1/eps - 1) Q(R)^2/(2 R) = -0.1014 Eh, 16 % off
   subroutine test_nonlocal_escaped_charge(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(pcm_physics) :: physics
      real(wp), parameter :: alpha = 0.2_wp, exact = -1.206972975269118e-1_wp
      real(wp) :: energy, deviation
      integer :: isurf
      character(len=160) :: context

      physics = pcm_physics(eps=80.0_wp, eps_inf=80.0_wp, lambda=3.0_wp)
      do isurf = 1, 2
         call central_charge_energy(pcm_nonlocal, physics, isurf == 2, solver_type%dense, &
            & 1.0_wp, energy, error, alpha=alpha)
         if (allocated(error)) return
         deviation = abs(energy - exact)/abs(exact)
         write (context, '(a,", energy ",es20.12,", exact ",es20.12,", rel. dev ",es9.2)') &
            & merge("point   ", "Gaussian", isurf == 2), energy, exact, deviation
         call check(error, deviation <= escaped_tolerance(isurf == 2), &
            & message="Nonlocal energy of a density outside the cavity deviates", &
            & more=trim(context))
         if (allocated(error)) return
      end do
   end subroutine test_nonlocal_escaped_charge

   !> Tolerance of the escaped-charge check, twice the measured deviation
   !>
   !> 9.8e-4 Gaussian (with the Gaussian-probed potential), 5.5e-3 point
   !>
   !> @param[in] point_surface Point instead of Gaussian surface charges
   pure function escaped_tolerance(point_surface) result(tol)
      !> Point instead of Gaussian surface charges
      logical, intent(in) :: point_surface
      !> Relative tolerance
      real(wp) :: tol

      tol = merge(1.2e-2_wp, 2.0e-3_wp, point_surface)
   end function escaped_tolerance

   !> LPB at inverse Debye lengths of 40 and 100 per bohr stays finite and Debye-Hueckel
   !>
   !> The exponentials of the smeared screened Coulomb kernel overflow at
   !> beta**2 = (kappa/(2 z))**2 > 700 unless evaluated as exp(-x**2) erfcx(beta - x)
   subroutine test_lpb_strong_screening(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(pcm_physics) :: physics
      real(wp) :: energy, reference, deviation
      real(wp), parameter :: kappa(2) = [40.0_wp, 100.0_wp]
      integer :: ikappa, isolver
      character(len=160) :: context

      do ikappa = 1, size(kappa)
         do isolver = solver_type%dense, solver_type%iterative
            physics%kappa = kappa(ikappa)
            call central_charge_energy(pcm_lpb, physics, .false., isolver, 1.0_wp, energy, error)
            if (allocated(error)) return
            reference = born_reference(pcm_lpb, physics, 1.0_wp, born_radius)
            deviation = abs(energy - reference)/abs(reference)
            write (context, '("kappa ",f6.1," solver ",i1," energy ",es16.8," reference ",es16.8, &
               & ", rel. dev ",es9.2)') kappa(ikappa), isolver, energy, reference, deviation
            call check(error, ieee_is_finite(energy), message="LPB energy is not finite", &
               & more=trim(context))
            if (allocated(error)) return
            call check(error, deviation <= lpb_strong_tolerance, &
               & message="LPB energy deviates from Debye-Hueckel", more=trim(context))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_lpb_strong_screening

   !> Nonlocal at correlation lengths of 0.01 and 1e-4 bohr stays finite
   !>
   !> The reference is the exact monopole energy of `born_reference`; its values
   !> for R = 4 bohr, eps = 78.4, eps_inf = 1.8 are pinned to the audit values
   subroutine test_nonlocal_strong_screening(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(pcm_physics) :: physics
      real(wp) :: energy, reference, deviation
      real(wp), parameter :: lambda(2) = [1.0e-2_wp, 1.0e-4_wp]
      real(wp), parameter :: audit(2) = [-0.12337991988706842_wp, -0.12340535522496877_wp]
      integer :: ilambda, isolver
      character(len=160) :: context

      do ilambda = 1, size(lambda)
         physics%lambda = lambda(ilambda)
         reference = born_reference(pcm_nonlocal, physics, 1.0_wp, born_radius)
         call check(error, reference, audit(ilambda), thr=1.0e-13_wp)
         if (allocated(error)) return
         do isolver = solver_type%dense, solver_type%iterative
            call central_charge_energy(pcm_nonlocal, physics, .false., isolver, 1.0_wp, energy, error)
            if (allocated(error)) return
            deviation = abs(energy - reference)/abs(reference)
            write (context, '("lambda ",es8.1," solver ",i1," energy ",es16.8," reference ",es16.8, &
               & ", rel. dev ",es9.2)') lambda(ilambda), isolver, energy, reference, deviation
            call check(error, ieee_is_finite(energy), message="Nonlocal energy is not finite", &
               & more=trim(context))
            if (allocated(error)) return
            call check(error, deviation <= nonlocal_strong_tolerance, &
               & message="Nonlocal energy deviates from the monopole result", more=trim(context))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_nonlocal_strong_screening

   !> Analytic energy of a central charge in a sphere
   !>
   !> - Conductor: -feps q**2/(2 R) with the model's dielectric factor
   !> - IEF-PCM, SS(V)PE and both DPCM forms: Born, feps = (eps - 1)/eps
   !> - LPB: Debye-Hueckel, feps = 1 - 1/(eps (1 + kappa R))
   !> - SLIC: feps = (eps - 1)/(eps - (eps - 1) h(q/R**2))
   !> - Nonlocal: exact monopole solution of the local reformulation with
   !>   Hildebrandt's approximate interface conditions (Bardhan, Knepley,
   !>   Brune, arXiv:1208.3866, eqs. 21-26) for interior permittivity eta,
   !>   feps = 1/eta - 1/eps - (1/eps_inf - 1/eps)/(1 + nu), nu = sqrt(eps/eps_inf) R/lambda;
   !>   the exterior sees only the displacement flux of the charge, so eta enters
   !>   through the interior potential 1/(eta r) alone (checked numerically
   !>   against the generated system at eta = 1, 2, 4 on 110 to 600 points)
   !> - The factor sinh(nu) exp(-nu)/nu in place of 1/(1 + nu) is not the
   !>   implemented model: it shares both lambda limits but differs by 5 % at nu ~ 1
   !>
   !> @param[in] kind    Model identifier
   !> @param[in] physics Physical parameters
   !> @param[in] charge  Central charge, e
   !> @param[in] radius  Sphere radius, bohr
   pure function born_reference(kind, physics, charge, radius) result(energy)
      !> Model identifier
      integer, intent(in) :: kind
      !> Physical parameters
      type(pcm_physics), intent(in) :: physics
      !> Central charge
      real(wp), intent(in) :: charge
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Polarization energy, hartree
      real(wp) :: energy
      !> Dielectric factor, SLIC perturbation and nonlocal ratio R/Lambda
      real(wp) :: feps, h, nu

      associate (eps => physics%eps)
         select case (kind)
         case (pcm_cosmo)
            feps = merge((eps - 1.0_wp)/(eps + 0.5_wp), 1.0_wp, ieee_is_finite(eps))
         case (pcm_lpb)
            feps = 1.0_wp - 1.0_wp/(eps*(1.0_wp + physics%kappa*radius))
         case (pcm_slic)
            h = physics%slic_alpha*tanh(physics%slic_beta*charge/radius**2 - physics%slic_gamma) &
               & + physics%slic_mu
            if (ieee_is_finite(eps)) then
               feps = (eps - 1.0_wp)/(eps - (eps - 1.0_wp)*h)
            else
               feps = 1.0_wp/(1.0_wp - h)
            end if
         case (pcm_nonlocal)
            nu = sqrt(eps/physics%eps_inf)*radius/physics%lambda
            feps = 1.0_wp/physics%eps_solute - 1.0_wp/eps &
               & - (1.0_wp/physics%eps_inf - 1.0_wp/eps)/(1.0_wp + nu)
         case default
            feps = merge((eps - 1.0_wp)/eps, 1.0_wp, ieee_is_finite(eps))
         end select
      end associate
      energy = -0.5_wp*feps*charge**2/radius
   end function born_reference

   !> Relative tolerance of the Born check on the 302-point sphere
   !>
   !> Point surface, all models and parameter sets of `born_cases`: 1 %, the
   !> discretization error of the point quadrature (the self terms) measured
   !> as 0.2 % for C-PCM, COSMO, IEF-PCM, SS(V)PE, 0.1 % DPCM, 0.13 % SLIC,
   !> 0.2-0.45 % LPB (largest at eps = 1, kappa = 0.3); nonlocal 2.5 %, measured
   !> 0.45-1.1 % (largest at lambda = 5 bohr). A point branch with the wrong
   !> dielectric is off by tens of percent at eps = 2
   !>
   !> Gaussian surface:
   !> Measured deviations on 110, 302 and 1202 points:
   !> - C-PCM and COSMO: 2e-13, 1e-13 and 3e-14; the Born-fitted zeta
   !>   (bornzeta.f90) makes the conductor sphere exact, so this guards the
   !>   table
   !> - IEF-PCM and SS(V)PE: 3e-8, 1e-8 and 3e-9; DPCM: 5e-7, 5e-8 and 1.4e-7
   !> - LPB: 9e-6, 4e-6 and 1e-6, the screened kernels
   !> - Nonlocal (no curvature terms): O(h) convergence, 0.7 % at 302 points
   !>   and lambda = 25 bohr, Richardson extrapolation to the reference within
   !>   0.1 %; the sinh form of the reference is 5 % away and fails
   !> - SLIC on the SvdW-DROP sphere: 5.6e-7, 5.7e-8 and 1.8e-7
   !>
   !> @param[in] kind    Model identifier
   !> @param[in] smeared Gaussian surface (else point surface)
   pure function born_tolerance(kind, smeared) result(tol)
      !> Model identifier
      integer, intent(in) :: kind
      !> Gaussian surface (else point surface)
      logical, intent(in) :: smeared
      !> Relative tolerance
      real(wp) :: tol

      if (.not. smeared) then
         tol = merge(2.5e-2_wp, 1.0e-2_wp, kind == pcm_nonlocal)
         return
      end if
      select case (kind)
      case (pcm_cpcm, pcm_cosmo)
         tol = 1.0e-11_wp
      case (pcm_lpb)
         tol = 1.0e-5_wp
      case (pcm_nonlocal)
         tol = 1.0e-2_wp
      case default
         tol = 1.0e-6_wp
      end select
   end function born_tolerance

   !* ================================================================================= *!
   !*                                  Conductor limit                                  *!
   !* ================================================================================= *!

   !> Energy of a unit charge at the center of the Born sphere
   !>
   !> @param[in]  kind    Model identifier
   !> @param[in]  physics Physical parameters
   !> @param[in]  smeared Gaussian surface (else point surface: widths and switching factors dropped)
   !> @param[out] energy  Polarization energy, hartree
   !> @param[out] error   Test failure
   subroutine sphere_energy(kind, physics, smeared, energy, error)
      !> Model identifier
      integer, intent(in) :: kind
      !> Physical parameters
      type(pcm_physics), intent(in) :: physics
      !> Gaussian surface
      logical, intent(in) :: smeared
      !> Polarization energy
      real(wp), intent(out) :: energy
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the cavity and the component
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Single-sphere solute
      type(structure_type) :: mol
      !> Spherical cavity
      class(cavity_type), allocatable :: cavity

      energy = 0.0_wp
      call new_context(ctx, verbosity=0)
      call new_pcm(kind, ctx, pcm, error, physics)
      if (allocated(error)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), charge=1.0_wp)
      call new_sphere_cavity(kind, [born_radius], born_nleb, ctx, cavity, error)
      if (allocated(error)) return
      call pcm%configure_cavity(cavity, err)
      call check_moist_error(error, err, "cavity configuration")
      if (allocated(error)) return
      call cavity%update(mol, error=err)
      call check_moist_error(error, err, "cavity update")
      if (allocated(error)) return
      if (.not. smeared) then
         if (allocated(cavity%xi0)) deallocate (cavity%xi0)
         if (allocated(cavity%f)) deallocate (cavity%f)
      end if
      call host_energy(pcm, mol, cavity, mol, [1.0_wp], energy, error)
   end subroutine sphere_energy

   !> Relative tolerance of the conductor-limit Born check
   !>
   !> Measured deviations on the 302-point sphere at eps = +Infinity:
   !> - Gaussian: C-PCM, COSMO, IEF-PCM and LPB 1e-13, SS(V)PE 3e-10, DPCM
   !>   (both forms) 5e-8, SLIC 6e-8 against its own limit 1/(1 - h)
   !> - Point: C-PCM, COSMO, IEF-PCM, SS(V)PE and LPB 2.0e-3, DPCM 1.0e-3,
   !>   SLIC 1.1e-3, the point quadrature error of the sphere
   !>
   !> @param[in] kind    Model identifier
   !> @param[in] smeared Gaussian surface
   pure function conductor_tolerance(kind, smeared) result(tol)
      !> Model identifier
      integer, intent(in) :: kind
      !> Gaussian surface
      logical, intent(in) :: smeared
      !> Relative tolerance
      real(wp) :: tol

      select case (kind)
      case (pcm_cpcm, pcm_cosmo, pcm_iefpcm, pcm_lpb)
         tol = merge(1.0e-11_wp, 3.0e-3_wp, smeared)
      case (pcm_ssvpe)
         tol = merge(1.0e-8_wp, 3.0e-3_wp, smeared)
      case default
         ! DPCM, both forms, and SLIC
         tol = merge(1.0e-6_wp, 1.5e-3_wp, smeared)
      end select
   end function conductor_tolerance

   !> Conductor limit eps = +Infinity of a central charge in a sphere
   !>
   !> - The energy is finite and equals the analytic f = 1 limit -q**2/(2R),
   !>   never an inf/inf evaluation (SLIC: the limit of its own factor)
   !> - It agrees with a very large finite eps, which catches a special case
   !>   that differs from the continuous limit (f = 0.99, an eps fallback, a
   !>   NaN from (eps - 1)/(eps + 1))
   !>
   !> @param[out] error Test failure
   !> @param[in]  kind  Model identifier
   subroutine check_conductor_limit(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind

      !> Physical parameters at infinite and at very large epsilon
      type(pcm_physics) :: inf_physics, big_physics
      !> Energy at infinite epsilon, at very large epsilon, analytic limit
      real(wp) :: energy, energy_big, reference, deviation, gap
      !> Surface index
      integer :: isurf
      !> Gaussian surface (else point surface)
      logical :: smeared
      !> Failure context
      character(len=200) :: context

      inf_physics%eps = ieee_value(0.0_wp, ieee_positive_inf)
      big_physics%eps = eps_large
      reference = born_reference(kind, inf_physics, 1.0_wp, born_radius)
      do isurf = 1, 2
         smeared = isurf == 1
         call sphere_energy(kind, inf_physics, smeared, energy, error)
         if (allocated(error)) return
         call sphere_energy(kind, big_physics, smeared, energy_big, error)
         if (allocated(error)) return
         deviation = abs(energy - reference)/abs(reference)
         gap = abs(energy - energy_big)/abs(energy_big)
         write (context, '(a," surface: energy ",es16.8," reference ",es16.8, &
            & ", rel. dev ",es9.2,", vs eps = 1e12 ",es9.2)') &
            & trim(merge("Gaussian", "point   ", smeared)), energy, reference, deviation, gap
         call check(error, ieee_is_finite(energy), &
            & message="energy at infinite epsilon is not finite", more=trim(context))
         if (allocated(error)) return
         call check(error, deviation <= conductor_tolerance(kind, smeared), &
            & message="energy at infinite epsilon deviates from the conductor limit", &
            & more=trim(context))
         if (allocated(error)) return
         call check(error, gap <= thr_large, &
            & message="energy at infinite epsilon differs from that at very large epsilon", &
            & more=trim(context))
         if (allocated(error)) return
      end do
   end subroutine check_conductor_limit

   !> A conductor-limit model equals CPCM at eps = +Infinity on the same molecular surface
   !>
   !> - Both solve the conductor boundary condition S q = -phi: IEF-PCM through
   !>   (2 pi - D A) (S q + phi) = 0, COSMO with f = 1 at its conductor limit
   !> - The identity holds on Gaussian and point surfaces to roundoff, on
   !>   asymmetric solutes without symmetry
   !>
   !> @param[out] error Test failure
   !> @param[in]  kind  Model identifier, IEF-PCM or COSMO
   subroutine check_conductor_identity(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind

      !> Run context borrowed by the components
      type(moist_context_type), target :: ctx
      !> Components under test
      class(model_continuum_component_pcm), allocatable :: pcm, ref
      !> Solute
      type(structure_type) :: mol
      !> Synthetic cavity
      type(fd_cavity) :: cavity
      !> Host charges
      real(wp), allocatable :: q(:)
      !> Physical parameters
      type(pcm_physics) :: physics
      !> Energies of the model and of CPCM
      real(wp) :: energy, energy_cpcm
      !> Solute and surface indices
      integer :: imol, isurf
      !> Gaussian surface (else point surface)
      logical :: smeared
      !> Failure context
      character(len=200) :: context

      physics%eps = ieee_value(0.0_wp, ieee_positive_inf)
      call new_context(ctx, verbosity=0)
      do imol = 1, nmol
         call get_solute(imol, mol)
         q = host_charges(mol)
         do isurf = 1, 2
            smeared = isurf == 1
            call get_seed(first_seed(kind), smeared, mol, cavity, error)
            if (allocated(error)) return
            call new_pcm(kind, ctx, pcm, error, physics)
            if (allocated(error)) return
            call host_energy(pcm, mol, cavity, mol, q, energy, error)
            if (allocated(error)) return
            call new_pcm(pcm_cpcm, ctx, ref, error, physics)
            if (allocated(error)) return
            call host_energy(ref, mol, cavity, mol, q, energy_cpcm, error)
            if (allocated(error)) return
            write (context, '(a,", ",a," surface: ",a," ",es22.14," CPCM ",es22.14," rel. dev ",es9.2)') &
               & trim(mol_label(imol)), trim(merge("Gaussian", "point   ", smeared)), &
               & trim(pcm%name), energy, energy_cpcm, abs(energy - energy_cpcm)/abs(energy_cpcm)
            call check(error, abs(energy_cpcm) > thr_vacuous, &
               & message="vacuous conductor energy", more=trim(context))
            if (allocated(error)) return
            call check(error, abs(energy - energy_cpcm) <= thr_conductor_identity*abs(energy_cpcm), &
               & message="energy at infinite epsilon differs from CPCM", more=trim(context))
            if (allocated(error)) return
         end do
      end do
   end subroutine check_conductor_identity

   !* ================================================================================= *!
   !*                                 Reference energies                                *!
   !* ================================================================================= *!

   !> Energies of the reference systems against the rows of the reference table
   !>
   !> Every row of the model is evaluated. A filled row must agree within its
   !> tolerance; pending rows skip the test and list the moist energies
   !> @param[in]  kind  Model identifier
   subroutine check_reference(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind

      !> Reference table
      type(reference_row), allocatable :: rows(:)
      !> Failure of one row's evaluation, collected and then cleared
      type(error_type), allocatable :: row_error
      !> Sphere centers and point-charge positions of an exact row
      type(structure_type) :: mol, src
      !> Point charges and sphere radii of an exact row
      real(wp), allocatable :: q(:), radii(:)
      !> Moist energy, reference energy and relative deviation
      real(wp) :: energy, reference, deviation
      !> Row index
      integer :: irow
      !> Context of one row, the failed rows and the pending rows
      character(len=320) :: context
      character(len=:), allocatable :: failed, pending

      rows = reference_table()
      failed = ""
      pending = ""
      do irow = 1, size(rows)
         if (rows(irow)%kind /= kind) cycle
         call reference_energy(kind, rows(irow)%system, energy, row_error)
         if (allocated(row_error)) then
            call add_line(failed, trim(ref_system_label(rows(irow)%system))//": "// &
               & row_error%message)
            ! handled: reported with the other rows below
            row_error%stat = 0
            deallocate (row_error)
            cycle
         end if
         reference = rows(irow)%energy
         if (rows(irow)%exact) then
            call reference_system(rows(irow)%system, mol, src, q, radii)
            reference = sphere_reference(kind, reference_physics(), radii(1), src%xyz, q)
         end if
         if (reference == ref_pending) then
            write (context, '(a,": moist ",es20.12," Eh, reference from ",a)') &
               & trim(ref_system_label(rows(irow)%system)), energy, trim(rows(irow)%source)
            call add_line(pending, trim(context))
            cycle
         end if
         deviation = abs(energy - reference)/abs(reference)
         write (context, '(a,": moist ",es20.12," vs ",es20.12," (",a,"), rel. dev ",es9.2)') &
            & trim(ref_system_label(rows(irow)%system)), energy, reference, &
            & trim(rows(irow)%source), deviation
         if (.not. deviation <= rows(irow)%rtol) call add_line(failed, trim(context))
      end do
      if (len(failed) > 0) then
         if (len(pending) > 0) call add_line(failed, "pending: "//pending)
         call test_failed(error, "energies deviate from the reference or failed", more=failed)
      else if (len(pending) > 0) then
         call skip_test(error, "reference energies pending", more=pending)
      end if
   end subroutine check_reference

   !> Append a line to a message, indented under the testdrive message
   !>
   !> @param[in,out] text Message
   !> @param[in]     line Line to append
   pure subroutine add_line(text, line)
      !> Message
      character(len=:), allocatable, intent(inout) :: text
      !> Line to append
      character(len=*), intent(in) :: line

      if (len(text) > 0) text = text//new_line("a")//"           "
      text = text//line
   end subroutine add_line

   !> Moist energy of a reference system on a Gaussian iSwiG surface
   !>
   !> @param[in]  kind   Model identifier
   !> @param[in]  system Reference system
   !> @param[out] energy Polarization energy, hartree
   !> @param[out] error  Test failure
   subroutine reference_energy(kind, system, energy, error)
      !> Model identifier
      integer, intent(in) :: kind
      !> Reference system
      integer, intent(in) :: system
      !> Polarization energy
      real(wp), intent(out) :: energy
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the cavity and the component
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Sphere centers and point-charge positions
      type(structure_type) :: mol, src
      !> Point charges and sphere radii
      real(wp), allocatable :: q(:), radii(:)
      !> Cavity
      class(cavity_type), allocatable :: cavity

      energy = 0.0_wp
      call reference_system(system, mol, src, q, radii)
      call new_context(ctx, verbosity=0)
      call new_pcm(kind, ctx, pcm, error, reference_physics())
      if (allocated(error)) return
      call new_sphere_cavity(kind, radii, ref_nleb, ctx, cavity, error)
      if (allocated(error)) return
      call pcm%configure_cavity(cavity, err)
      call check_moist_error(error, err, "cavity configuration")
      if (allocated(error)) return
      call cavity%update(mol, error=err)
      call check_moist_error(error, err, "cavity update")
      if (allocated(error)) return
      call host_energy(pcm, mol, cavity, src, q, energy, error)
   end subroutine reference_energy

   !> Rows of the smoke test: energies of other codes and exact series
   !>
   !> Every row reproduces outside moist from `reference_system` and
   !> `reference_physics`: point charges (no nuclei, no electrons) in the
   !> union of spheres of the listed radii, no radius scaling, energy
   !> E = 1/2 sum q_i phi_reaction(r_i) in hartree. A `ref_pending` row lists
   !> the moist energy; `rtol` is the expected difference of the two
   !> discretizations, `source` names the code and version or the series
   !>
   !> - Single sphere: exact series (`sphere_reference`) for every model but
   !>   SLIC. Measured deviation on 302 points: 1.8e-4 (CPCM, COSMO, IEF-PCM,
   !>   SS(V)PE, LPB) and 4e-5 (DPCM), falling to 2.4e-5 and 5e-6 on 1202
   !> - Nonlocal sphere: 2.5e-3 on 302 points and not monotone in the grid
   !>   (7.6e-4 on 110, 2.6e-3 on 434, 2.1e-3 on 1202). Without curvature
   !>   terms every mode converges as O(h) with its own sign; the monopole
   !>   alone (central ion, lambda = 10 bohr) is 1.6 % off on 302 points and
   !>   0.85 % on 1202. A gross-error guard only
   !> - CPCM, COSMO, IEF-PCM, SS(V)PE on molecules: PySCF 2.13.1 on the same
   !>   iSwiG surface (302 points, ISWIG requested explicitly), driven with
   !>   point charges through its operator intermediates. Both use the same
   !>   XI (bornzeta.f90); the difference of at most 5e-10 is the pruning
   !>   threshold, and with it matched the energies agree to 3e-13
   !> - DPCM (both forms) on molecules: the dielectric continuum limit, shared
   !>   with IEF-PCM and SS(V)PE (PySCF IEF-PCM on up to 5810 points,
   !>   Lebedev-extrapolated; APBS 3.4.1 at kappa = 0 agrees to 3.4e-4)
   !> - LPB on molecules: APBS 3.4.1, sharp sphere union (`srfm mol`,
   !>   `srad 0`), ion radius 0, Debye length 7.85129 A; extrapolated from grid
   !>   spacings 0.025 to 0.15 A over 3 grid shifts (3.7e-4 from the exact
   !>   sphere series)
   !> - DPCM and LPB tolerances: 1.5 times the measured 302-point
   !>   discretization error of moist plus twice the reference uncertainty,
   !>   0.15 % to 3.6 %; benzene is the loosest, every model there is 2 % from
   !>   the continuum on 302 points
   !> - Not used: domain-decomposition codes (ddX), whose discretization of the
   !>   boundary problem differs. PCMSolver 1.2.3 (GePol) confirms the
   !>   continuum values where its S matrix stays positive definite
   !> - SLIC: evaluated on the SvdW-DROP cavity (default level set), where
   !>   iSwiG seams fold its interface condition. Sphere: Bardhan's
   !>   testasymmetry/pointbem converged to 1e-10 and extrapolated in the point
   !>   density, with the normal field just inside (his code takes the principal
   !>   value; one jump line enabled, plus his K' sigma weight fixed for unequal
   !>   weights); uncertainty 2.9e-5, 3.3e-5 from moist on 302 points.
   !>   Molecules: the same code on moist's own DROP points (positions,
   !>   normals, areas; the analytic self-patch also takes moist's curvatures)
   !>   at Lebedev 302 to 5810, so it checks the electrostatics and the
   !>   interface condition but shares the surface quadrature. Reference: the
   !>   finest order where his Newton converges (3074 to 5810; on anions it
   !>   stalls at some orders, at 1e-2 residual). There it is 2.3e-4 to 1.1e-3
   !>   from moist on the same points, as in the linear model (7e-4 to 3.8e-3),
   !>   and the SLIC part (SLIC minus IEF-PCM) agrees to 7e-4 to 1.7e-2.
   !>   Tolerance: 1.5 times moist's 302-point error on DROP (its own finest
   !>   order) plus twice that gap; benzene, 2 % off on 302 points, is a
   !>   gross check
   !> - Nonlocal molecules: no reference, rows commented out until one exists.
   function reference_table() result(rows)
      !> Rows of every model
      type(reference_row), allocatable :: rows(:)

      !> Tolerance of the PySCF rows, 3 times the largest difference (PySCF
      !> zeroes switching factors below 1e-8, moist cuts points at 1e-10)
      real(wp), parameter :: pyscf_rtol(ref_water:nref_system) = 1.5e-9_wp
      !> Dielectric continuum energies, the DPCM references
      real(wp), parameter :: dielectric_energy(ref_water:nref_system) = [ &
         & -9.247244539606e-3_wp, -5.986466050790e-2_wp, -1.247026222991e-1_wp, &
         & -9.797658606323e-2_wp, -1.158654839645e-1_wp, -3.812554143673e-1_wp, &
         & -7.832255674277e-2_wp, -3.115387503725e-3_wp]
      !> Tolerance of the DPCM rows
      real(wp), parameter :: dpcm_rtol(ref_water:nref_system) = [ &
         & 1.2e-2_wp, 4.3e-3_wp, 1.9e-2_wp, 2.9e-2_wp, 1.5e-2_wp, 1.8e-3_wp, 2.2e-2_wp, 3.6e-2_wp]
      !> Source of the DPCM rows
      character(len=*), parameter :: dielectric_source = &
         & "dielectric continuum: PySCF 2.13.1 IEF-PCM, Lebedev-extrapolated; APBS 3.4.1 kappa = 0"

      allocate (rows(0))
      rows = [rows, reference_row(pcm_cpcm, ref_sphere, exact=.true., rtol=5.0e-4_wp, &
         & source="exact: Kelvin images times (eps - 1)/eps")]
      call add_reference_rows(rows, pcm_cpcm, &
         & "PySCF 2.13.1 C-PCM, iSwiG 302 (point-charge driver, same XI)", [ &
         & -9.282588447688e-3_wp, -6.002336394632e-2_wp, -1.244313700754e-1_wp, &
         & -9.787375677457e-2_wp, -1.157020380190e-1_wp, -3.806497039596e-1_wp, &
         & -7.824217756579e-2_wp, -3.088155427039e-3_wp], pyscf_rtol)

      rows = [rows, reference_row(pcm_cosmo, ref_sphere, exact=.true., rtol=5.0e-4_wp, &
         & source="exact: Kelvin images times (eps - 1)/(eps + 1/2)")]
      call add_reference_rows(rows, pcm_cosmo, &
         & "PySCF 2.13.1 COSMO, iSwiG 302 (point-charge driver, same XI)", [ &
         & -9.223755969251e-3_wp, -5.964293953292e-2_wp, -1.236427316544e-1_wp, &
         & -9.725343888399e-2_wp, -1.149687255712e-1_wp, -3.782371693928e-1_wp, &
         & -7.774628342480e-2_wp, -3.068582886622e-3_wp], pyscf_rtol)

      rows = [rows, reference_row(pcm_iefpcm, ref_sphere, exact=.true., rtol=5.0e-4_wp, &
         & source="exact: Kirkwood series")]
      call add_reference_rows(rows, pcm_iefpcm, &
         & "PySCF 2.13.1 IEF-PCM, iSwiG 302 (point-charge driver, same XI)", [ &
         & -9.230680472207e-3_wp, -5.983386654374e-2_wp, -1.244240617773e-1_wp, &
         & -9.783600680854e-2_wp, -1.156218073482e-1_wp, -3.805956547272e-1_wp, &
         & -7.795047884686e-2_wp, -3.052371993565e-3_wp], pyscf_rtol)

      rows = [rows, reference_row(pcm_ssvpe, ref_sphere, exact=.true., rtol=5.0e-4_wp, &
         & source="exact: Kirkwood series")]
      call add_reference_rows(rows, pcm_ssvpe, &
         & "PySCF 2.13.1 SS(V)PE, iSwiG 302 (point-charge driver, same XI)", [ &
         & -9.213745328787e-3_wp, -5.983850556443e-2_wp, -1.215136755362e-1_wp, &
         & -9.463881031330e-2_wp, -1.133424390666e-1_wp, -3.704190466665e-1_wp, &
         & -7.728286938272e-2_wp, -2.989230591327e-3_wp], pyscf_rtol)

      rows = [rows, reference_row(pcm_dpcm, ref_sphere, exact=.true., rtol=5.0e-4_wp, &
         & source="exact: Kirkwood series")]
      call add_reference_rows(rows, pcm_dpcm, dielectric_source, dielectric_energy, dpcm_rtol)

      rows = [rows, reference_row(pcm_dpcm_original, ref_sphere, exact=.true., rtol=5.0e-4_wp, &
         & source="exact: Kirkwood series")]
      call add_reference_rows(rows, pcm_dpcm_original, dielectric_source, dielectric_energy, &
         & dpcm_rtol)

      rows = [rows, reference_row(pcm_lpb, ref_sphere, exact=.true., rtol=5.0e-4_wp, &
         & source="exact: Kirkwood series with salt")]
      call add_reference_rows(rows, pcm_lpb, &
         & "APBS 3.4.1, grid-extrapolated (srfm mol, srad 0, ion radius 0)", [ &
         & -9.249964431513e-3_wp, -5.989640506049e-2_wp, -1.250557810616e-1_wp, &
         & -9.827887845846e-2_wp, -1.161912979478e-1_wp, -3.825071303454e-1_wp, &
         & -7.840535566853e-2_wp, -3.118313818774e-3_wp], [ &
         & 3.6e-3_wp, 1.5e-3_wp, 4.8e-3_wp, 3.4e-3_wp, 4.3e-3_wp, 3.7e-3_wp, 8.8e-3_wp, 3.6e-2_wp])

      rows = [rows, reference_row(pcm_slic, ref_sphere, energy=-3.997005781069e-2_wp, rtol=1.1e-4_wp, &
         & source="Bardhan testasymmetry/pointbem, E_n inside, converged and extrapolated")]
      call add_reference_rows(rows, pcm_slic, &
         & "Bardhan testasymmetry/pointbem on moist's SvdW-DROP points, E_n inside", [ &
         & -9.543944718520e-3_wp, -6.744635159891e-2_wp, -8.305701409185e-2_wp, &
         & -7.074330147437e-2_wp, -1.505648524855e-1_wp, -4.945038921462e-1_wp, &
         & -7.632376814487e-2_wp, -3.061163654352e-3_wp], [ &
         & 4.7e-3_wp, 3.0e-3_wp, 6.0e-3_wp, 3.9e-3_wp, 3.3e-3_wp, 1.2e-2_wp, 2.9e-3_wp, 3.3e-2_wp])

      rows = [rows, reference_row(pcm_nonlocal, ref_sphere, exact=.true., rtol=5.0e-3_wp, &
         & source="exact: per-mode series, Hildebrandt interface conditions")]

      ! No reference for nonlocal yet
      ! call add_reference_rows(rows, pcm_nonlocal, "none (NESSie.jl: check its Born limit)")

   end function reference_table

   !> Append the rows of one model on the molecular systems, water to benzene
   !>
   !> @param[in,out] rows     Reference table
   !> @param[in]     kind     Model identifier
   !> @param[in]     source   Origin of the references
   !> @param[in]     energies Reference energies, hartree, in system order; absent: pending rows
   !> @param[in]     rtol     Relative tolerance of each row, present with `energies`
   subroutine add_reference_rows(rows, kind, source, energies, rtol)
      !> Reference table
      type(reference_row), allocatable, intent(inout) :: rows(:)
      !> Model identifier
      integer, intent(in) :: kind
      !> Origin of the references
      character(len=*), intent(in) :: source
      !> Reference energies
      real(wp), intent(in), optional :: energies(ref_water:)
      !> Relative tolerances
      real(wp), intent(in), optional :: rtol(ref_water:)
      !> Row of one system
      type(reference_row) :: row
      !> System index
      integer :: system

      do system = ref_water, nref_system
         row = reference_row(kind, system, source=source)
         if (present(energies)) then
            row%energy = energies(system)
            row%rtol = rtol(system)
         end if
         rows = [rows, row]
      end do
   end subroutine add_reference_rows

   !> Exact energy of point charges inside a sphere centered at the origin
   !>
   !> Multipole series, solute permittivity 1:
   !> E = 1/2 sum_ij q_i q_j sum_n c_n (r_i r_j)**n P_n(cos g_ij)/a**(2n + 1)
   !>
   !> - Conductor (Kelvin images): c_n = -1, times the dielectric factor of
   !>   CPCM or COSMO
   !> - Dielectric (Kirkwood): c_n = (n + 1)(1 - eps)/((n + 1) eps + n)
   !> - LPB with the ions excluded from the sphere (Kirkwood): c_n = (g_n + n + 1)/(n - g_n),
   !>   g_n = eps x k_n'(x)/k_n(x), x = kappa a, k_n the modified spherical
   !>   Bessel function of the second kind
   !> - Nonlocal (local reformulation, Hildebrandt's interface conditions, as
   !>   in `born_reference`): c_n = (n + 1)(F_n - 1)/(n + 1 + n F_n),
   !>   F_n = 1/eps - (n + 1)(1/eps_inf - 1/eps)/g_n, g_n = x k_n'(x)/k_n(x),
   !>   x = a sqrt(eps/eps_inf)/lambda; n = 0 is the Born form
   !> - SLIC: `ref_pending`
   !>
   !> @param[in] kind    Model identifier
   !> @param[in] physics Physical parameters
   !> @param[in] radius  Sphere radius, bohr
   !> @param[in] xyz     Charge positions inside the sphere (3, n), bohr
   !> @param[in] q       Point charges (n)
   function sphere_reference(kind, physics, radius, xyz, q) result(energy)
      !> Model identifier
      integer, intent(in) :: kind
      !> Physical parameters
      type(pcm_physics), intent(in) :: physics
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Charge positions
      real(wp), intent(in) :: xyz(:, :)
      !> Point charges
      real(wp), intent(in) :: q(:)
      !> Polarization energy, hartree
      real(wp) :: energy

      !> Highest multipole order, far past convergence for charges 1.5 bohr inside 4 bohr
      integer, parameter :: nmax = 60
      !> Reaction coefficients
      real(wp) :: c(0:nmax)
      !> Legendre polynomials of the pair angle
      real(wp) :: p(0:nmax)
      !> Dielectric factor, x = kappa a, ratio k_(n-1)/k_n and x k_n'/k_n
      real(wp) :: feps, x, ratio, dlog
      !> Nonlocal: 1/eps_inf - 1/eps and the mode factor F_n
      real(wp) :: delta, f
      !> Pair cosine and radial factor
      real(wp) :: cosg, t
      !> Charge and order indices
      integer :: i, j, n

      associate (eps => physics%eps)
         feps = 1.0_wp
         select case (kind)
         case (pcm_cpcm, pcm_cosmo)
            c = -1.0_wp
            feps = merge((eps - 1.0_wp)/eps, (eps - 1.0_wp)/(eps + 0.5_wp), kind == pcm_cpcm)
         case (pcm_iefpcm, pcm_ssvpe, pcm_dpcm, pcm_dpcm_original)
            do n = 0, nmax
               c(n) = real(n + 1, wp)*(1.0_wp - eps)/(real(n + 1, wp)*eps + real(n, wp))
            end do
         case (pcm_lpb)
            ! k_0 = exp(-x)/x, k_1 = exp(-x)(1 + x)/x**2, k_(n+1) = k_(n-1) + (2n + 1)/x k_n,
            ! k_n' = -k_(n-1) - (n + 1)/x k_n, carried as the ratio k_(n-1)/k_n
            x = physics%kappa*radius
            dlog = -(1.0_wp + x)
            c(0) = (eps*dlog + 1.0_wp)/(-eps*dlog)
            ratio = x/(1.0_wp + x)
            do n = 1, nmax
               dlog = -x*ratio - real(n + 1, wp)
               c(n) = (eps*dlog + real(n + 1, wp))/(real(n, wp) - eps*dlog)
               ratio = 1.0_wp/(ratio + real(2*n + 1, wp)/x)
            end do
         case (pcm_nonlocal)
            ! Same recursion as LPB with x = a/Lambda; dlog and ratio advance to n + 1
            x = radius*sqrt(eps/physics%eps_inf)/physics%lambda
            delta = 1.0_wp/physics%eps_inf - 1.0_wp/eps
            dlog = -(1.0_wp + x)
            ratio = x/(1.0_wp + x)
            do n = 0, nmax
               f = 1.0_wp/eps - real(n + 1, wp)*delta/dlog
               c(n) = real(n + 1, wp)*(f - 1.0_wp)/(real(n + 1, wp) + real(n, wp)*f)
               dlog = -x*ratio - real(n + 2, wp)
               ratio = 1.0_wp/(ratio + real(2*n + 3, wp)/x)
            end do
         case default
            energy = ref_pending
            return
         end select
      end associate

      energy = 0.0_wp
      do i = 1, size(q)
         do j = 1, size(q)
            t = norm2(xyz(:, i))*norm2(xyz(:, j))/radius**2
            cosg = 1.0_wp
            if (t > 0.0_wp) cosg = dot_product(xyz(:, i), xyz(:, j))/(t*radius**2)
            p(0) = 1.0_wp
            p(1) = cosg
            do n = 1, nmax - 1
               p(n + 1) = (real(2*n + 1, wp)*cosg*p(n) - real(n, wp)*p(n - 1))/real(n + 1, wp)
            end do
            energy = energy + 0.5_wp*q(i)*q(j)*sum(c*t**[(n, n=0, nmax)]*p)/radius
         end do
      end do
      energy = feps*energy
   end function sphere_reference

   !> Physical parameters of the smoke test
   !>
   !> - Water at 298.15 K: eps = 78.39
   !> - LPB: 0.15 M 1:1 salt, Debye length 7.855 A, kappa = 0.0674 / bohr
   !> - Nonlocal: eps_inf = 1.8, lambda = 10 bohr
   !> - SLIC: Bardhan and Knepley, J. Chem. Phys. 141, 131103 (2014), water:
   !>   alpha 0.5, gamma -0.5, mu = -alpha tanh(-gamma) (h(0) = 0) and
   !>   beta -60 per e/(4 pi A**2), his Green's function being 1/(4 pi r),
   !>   = -60/(4 pi) (A/bohr)**2 per e/bohr**2
   pure function reference_physics() result(physics)
      !> Physical parameters
      type(pcm_physics) :: physics

      physics%eps = 78.39_wp
      physics%kappa = 0.0674_wp
      physics%eps_inf = 1.8_wp
      physics%lambda = 10.0_wp
      physics%slic_alpha = 0.5_wp
      physics%slic_beta = -1.705057857513347e1_wp
      physics%slic_gamma = -0.5_wp
      physics%slic_mu = -2.310585786300049e-1_wp
   end function reference_physics

   !> Point charges and spheres of a reference system
   !>
   !> Charges on the sphere centers except for the single sphere; radii Bondi
   !> times 1.2 (`reference_radius`); partial charges are fixed numbers summing
   !> to the net charge, not fitted to any electron density
   !>
   !> - Off-center pair: one sphere of radius 4 bohr at the origin, +1 at
   !>   (0, 0, 1.5) bohr and -0.5 at (1, 0.5, -1) bohr
   !> - TIP3P water: r(OH) = 0.9572 A, HOH = 104.52 deg, charges -0.834 and +0.417
   !> - Na+ Cl- contact pair: 2.8 A apart, charges +1 and -1
   !> - Ammonium: tetrahedral, N-H 1.03 A
   !> - Guanidinium: planar D3h, C-N 1.33 A, N-H 1.01 A
   !> - Acetate: C-O 1.26 A, OCO 125 deg, C-C 1.52 A, staggered methyl
   !> - Sulfate: tetrahedral, S-O 1.49 A
   !> - Glycine zwitterion: heavy atoms planar, N-CA 1.49, CA-C 1.53, C-O 1.25 A
   !> - Benzene: C-C 1.39 A, C-H 1.09 A
   !>
   !> @param[in]  system Reference system
   !> @param[out] mol    Sphere centers
   !> @param[out] src    Point-charge positions
   !> @param[out] q      Point charges
   !> @param[out] radii  Sphere radii, bohr
   subroutine reference_system(system, mol, src, q, radii)
      !> Reference system
      integer, intent(in) :: system
      !> Sphere centers
      type(structure_type), intent(out) :: mol
      !> Point-charge positions
      type(structure_type), intent(out) :: src
      !> Point charges
      real(wp), allocatable, intent(out) :: q(:)
      !> Sphere radii
      real(wp), allocatable, intent(out) :: radii(:)

      !> TIP3P bond length (bohr) and half angle (rad)
      real(wp), parameter :: r_oh = 0.9572_wp*aatoau, half_hoh = 0.5_wp*104.52_wp*acos(-1.0_wp)/180.0_wp
      !> Atomic numbers of the sphere centers
      integer, allocatable :: num(:)
      !> Sphere centers, bohr
      real(wp), allocatable :: xyz(:, :)
      !> Atom index
      integer :: i

      select case (system)
      case (ref_sphere)
         call new(mol, [11], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
         call new(src, [1, 1], reshape([0.0_wp, 0.0_wp, 1.5_wp, 1.0_wp, 0.5_wp, -1.0_wp], [3, 2]))
         q = [1.0_wp, -0.5_wp]
         radii = [4.0_wp]
         mol%charge = sum(q)
         return
      case (ref_water)
         num = [8, 1, 1]
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, &
            & r_oh*sin(half_hoh), 0.0_wp, r_oh*cos(half_hoh), &
            & -r_oh*sin(half_hoh), 0.0_wp, r_oh*cos(half_hoh)], [3, 3])
         q = [-0.834_wp, 0.417_wp, 0.417_wp]
      case (ref_ion_pair)
         num = [11, 17]
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.8_wp*aatoau], [3, 2])
         q = [1.0_wp, -1.0_wp]
      case (ref_ammonium)
         num = [7, 1, 1, 1, 1]
         xyz = reshape([0.000000_wp, 0.000000_wp, 0.000000_wp, &
            & 0.594671_wp, 0.594671_wp, 0.594671_wp, &
            & 0.594671_wp, -0.594671_wp, -0.594671_wp, &
            & -0.594671_wp, 0.594671_wp, -0.594671_wp, &
            & -0.594671_wp, -0.594671_wp, 0.594671_wp], [3, 5])*aatoau
         q = [-0.4_wp, 0.35_wp, 0.35_wp, 0.35_wp, 0.35_wp]
      case (ref_guanidinium)
         num = [6, 7, 7, 7, 1, 1, 1, 1, 1, 1]
         xyz = reshape([0.000000_wp, 0.000000_wp, 0.000000_wp, &
            & 0.000000_wp, 1.330000_wp, 0.000000_wp, &
            & -1.151814_wp, -0.665000_wp, 0.000000_wp, &
            & 1.151814_wp, -0.665000_wp, 0.000000_wp, &
            & -0.874686_wp, 1.835000_wp, 0.000000_wp, &
            & 0.874686_wp, 1.835000_wp, 0.000000_wp, &
            & -1.151814_wp, -1.675000_wp, 0.000000_wp, &
            & -2.026499_wp, -0.160000_wp, 0.000000_wp, &
            & 2.026499_wp, -0.160000_wp, 0.000000_wp, &
            & 1.151814_wp, -1.675000_wp, 0.000000_wp], [3, 10])*aatoau
         q = [0.7_wp, -0.8_wp, -0.8_wp, -0.8_wp, 0.45_wp, &
            & 0.45_wp, 0.45_wp, 0.45_wp, 0.45_wp, 0.45_wp]
      case (ref_acetate)
         num = [6, 8, 8, 6, 1, 1, 1]
         xyz = reshape([0.000000_wp, 0.000000_wp, 0.000000_wp, &
            & 1.117634_wp, 0.000000_wp, 0.581803_wp, &
            & -1.117634_wp, 0.000000_wp, 0.581803_wp, &
            & 0.000000_wp, 0.000000_wp, -1.520000_wp, &
            & 0.513740_wp, -0.889823_wp, -1.883849_wp, &
            & -1.027479_wp, 0.000000_wp, -1.883849_wp, &
            & 0.513740_wp, 0.889823_wp, -1.883849_wp], [3, 7])*aatoau
         q = [0.7_wp, -0.8_wp, -0.8_wp, -0.28_wp, 0.06_wp, 0.06_wp, 0.06_wp]
      case (ref_sulfate)
         num = [16, 8, 8, 8, 8]
         xyz = reshape([0.000000_wp, 0.000000_wp, 0.000000_wp, &
            & 0.860252_wp, 0.860252_wp, 0.860252_wp, &
            & 0.860252_wp, -0.860252_wp, -0.860252_wp, &
            & -0.860252_wp, 0.860252_wp, -0.860252_wp, &
            & -0.860252_wp, -0.860252_wp, 0.860252_wp], [3, 5])*aatoau
         q = [1.2_wp, -0.8_wp, -0.8_wp, -0.8_wp, -0.8_wp]
      case (ref_glycine)
         num = [7, 6, 6, 8, 8, 1, 1, 1, 1, 1]
         xyz = reshape([1.490000_wp, 0.000000_wp, 0.000000_wp, &
            & 0.000000_wp, 0.000000_wp, 0.000000_wp, &
            & -0.548303_wp, -1.428378_wp, 0.000000_wp, &
            & 0.288110_wp, -2.357309_wp, 0.000000_wp, &
            & -1.791455_wp, -1.559039_wp, 0.000000_wp, &
            & 1.842281_wp, 0.967883_wp, 0.000000_wp, &
            & 1.842281_wp, -0.483942_wp, 0.838212_wp, &
            & 1.842281_wp, -0.483942_wp, -0.838212_wp, &
            & -0.353907_wp, 0.515969_wp, -0.892539_wp, &
            & -0.353907_wp, 0.515969_wp, 0.892539_wp], [3, 10])*aatoau
         q = [-0.3_wp, 0.06_wp, 0.72_wp, -0.795_wp, -0.795_wp, &
            & 0.33_wp, 0.33_wp, 0.33_wp, 0.06_wp, 0.06_wp]
      case (ref_benzene)
         num = [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1]
         xyz = reshape([1.390000_wp, 0.000000_wp, 0.000000_wp, &
            & 0.695000_wp, 1.203775_wp, 0.000000_wp, &
            & -0.695000_wp, 1.203775_wp, 0.000000_wp, &
            & -1.390000_wp, 0.000000_wp, 0.000000_wp, &
            & -0.695000_wp, -1.203775_wp, 0.000000_wp, &
            & 0.695000_wp, -1.203775_wp, 0.000000_wp, &
            & 2.480000_wp, 0.000000_wp, 0.000000_wp, &
            & 1.240000_wp, 2.147743_wp, 0.000000_wp, &
            & -1.240000_wp, 2.147743_wp, 0.000000_wp, &
            & -2.480000_wp, 0.000000_wp, 0.000000_wp, &
            & -1.240000_wp, -2.147743_wp, 0.000000_wp, &
            & 1.240000_wp, -2.147743_wp, 0.000000_wp], [3, 12])*aatoau
         q = [-0.115_wp, -0.115_wp, -0.115_wp, -0.115_wp, -0.115_wp, -0.115_wp, &
            & 0.115_wp, 0.115_wp, 0.115_wp, 0.115_wp, 0.115_wp, 0.115_wp]
      case default
         error stop "reference_system: unknown system"
      end select
      call new(mol, num, xyz)
      src = mol
      radii = [(reference_radius(num(i)), i=1, size(num))]*aatoau
      mol%charge = sum(q)
   end subroutine reference_system

   !> Sphere radius of an element in the reference systems: Bondi times 1.2
   !>
   !> @param[in] num Atomic number
   pure function reference_radius(num) result(radius)
      !> Atomic number
      integer, intent(in) :: num
      !> Radius, angstrom
      real(wp) :: radius

      select case (num)
      case (1)
         radius = 1.44_wp
      case (6)
         radius = 2.04_wp
      case (7)
         radius = 1.86_wp
      case (8)
         radius = 1.824_wp
      case (11)
         radius = 2.724_wp
      case (16)
         radius = 2.16_wp
      case (17)
         radius = 2.10_wp
      case default
         error stop "reference_radius: element without a reference radius"
      end select
   end function reference_radius

   !> Cavity of spheres with explicit radii: SvdW-DROP (default level set) for SLIC, else iSwiG
   !>
   !> On a single sphere both are the sphere itself
   !>
   !> @param[in]  kind    Model identifier
   !> @param[in]  radii   Sphere radii, bohr
   !> @param[in]  num_leb Lebedev order per sphere
   !> @param[in]  ctx     Run context borrowed by the cavity
   !> @param[out] cavity  Cavity, not yet configured
   !> @param[out] error   Test failure
   subroutine new_sphere_cavity(kind, radii, num_leb, ctx, cavity, error)
      !> Model identifier
      integer, intent(in) :: kind
      !> Sphere radii
      real(wp), intent(in) :: radii(:)
      !> Lebedev order
      integer, intent(in) :: num_leb
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(in) :: ctx
      !> Cavity
      class(cavity_type), allocatable, intent(out) :: cavity
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Radius model of the explicit radii
      class(radius_type), allocatable :: radius_model
      !> Concrete cavities
      type(cavity_type_iswig), allocatable :: iswig
      type(cavity_type_drop), allocatable :: drop
      !> Level set of the DROP cavity
      type(moist_cavity_drop_lsf_svdw_type) :: svdw

      call new_radii_custom_atoms(radii, radius_model, err)
      call check_moist_error(error, err, "radius model")
      if (allocated(error)) return
      if (first_seed(kind) == seed_drop) then
         allocate (drop)
         call svdw%new()
         call new_cavity_drop(drop, radius_model=radius_model, lsf_model=svdw, error=err, &
            & param=moist_cavity_drop_parameters_type(num_leb=num_leb), ctx=ctx)
         call check_moist_error(error, err, "DROP setup")
         if (allocated(error)) return
         call move_alloc(drop, cavity)
      else
         allocate (iswig)
         call new_cavity_iswig(iswig, radius_model=radius_model, error=err, &
            & param=moist_cavity_iswig_parameters_type(num_leb=num_leb), ctx=ctx)
         call check_moist_error(error, err, "iSwiG setup")
         if (allocated(error)) return
         call move_alloc(iswig, cavity)
      end if
   end subroutine new_sphere_cavity

   !* ================================================================================= *!
   !*                                 Switched-off point                                *!
   !* ================================================================================= *!

   !> A point switched off on the surface drops out of the energy
   !>
   !> One point of each solute's surface (SLIC: SvdW-DROP, see `first_seed`) is
   !> switched off by a factor s: on a Gaussian surface its switching factor
   !> and its area (a = F R^2 w as on iSwiG, the width kept), on a point surface
   !> its area. The energy must approach that of the surface without the point,
   !> the gap shrinking at least like s^(1/2) between the two smallest factors, the
   !> order of the point-surface self terms 1.0694 sqrt(4 pi/a). A gap that stays
   !> finite is a jump in the energy wherever a cavity drops a point that is nearly
   !> switched off
   !> @param[in]  kind  Model identifier
   subroutine check_switched_off_point(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind

      !> Run context borrowed by the components
      type(moist_context_type), target :: ctx
      !> Component under test, fresh for every surface
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Solute
      type(structure_type) :: mol
      !> Synthetic cavity and the surface without or with the switched-off point
      type(fd_cavity) :: cavity, trial
      !> Host charges on the solute atoms
      real(wp), allocatable :: q(:)
      !> Energy without the point, energies and gaps with the point switched off
      real(wp) :: e_removed, e_off(3), gap(3)
      !> Solute, surface kind, point and factor indices
      integer :: imol, isurf, ipt, k
      !> Gaussian surface (else point surface)
      logical :: smeared
      !> Failure context
      character(len=200) :: context

      call new_context(ctx, verbosity=0)
      do imol = 1, nmol
         call get_solute(imol, mol)
         q = host_charges(mol)
         do isurf = 1, 2
            smeared = isurf == 1
            call get_seed(first_seed(kind), smeared, mol, cavity, error)
            if (allocated(error)) return
            ipt = 1 + cavity%ngrid/3
            call remove_point(cavity, ipt, trial)
            call new_pcm(kind, ctx, pcm, error)
            if (allocated(error)) return
            call host_energy(pcm, mol, trial, mol, q, e_removed, error)
            if (allocated(error)) return
            do k = 1, size(switch_off)
               trial = cavity
               trial%a(ipt) = switch_off(k)*cavity%a(ipt)
               if (smeared) trial%f(ipt) = switch_off(k)*cavity%f(ipt)
               call new_pcm(kind, ctx, pcm, error)
               if (allocated(error)) return
               call host_energy(pcm, mol, trial, mol, q, e_off(k), error)
               if (allocated(error)) return
            end do
            gap = e_off - e_removed
            write (context, '(a,", ",a,", ",a," surface, point ",i0,": gap ",3(es9.2, &
               & " at s = ",es7.1,:,", "))') trim(pcm%name), trim(mol_label(imol)), &
               & trim(merge("Gaussian", "point   ", smeared)), ipt, &
               & (gap(k), switch_off(k), k=1, size(switch_off))
            call check(error, abs(gap(1)) > thr_gap, message="vacuous switched-off point", &
               & more=trim(context))
            if (allocated(error)) return
            call check(error, abs(gap(3)) <= 2*sqrt(switch_off(3)/switch_off(2))*abs(gap(2)) &
               & + thr_gap, message="switched-off point does not drop out of the energy", &
               & more=trim(context))
            if (allocated(error)) return
         end do
      end do
   end subroutine check_switched_off_point

   !> Copy of a synthetic cavity without one of its points
   !>
   !> @param[in]  cavity  Synthetic cavity
   !> @param[in]  ipt     Point to remove
   !> @param[out] removed Copy without the point
   subroutine remove_point(cavity, ipt, removed)
      !> Synthetic cavity
      type(fd_cavity), intent(in) :: cavity
      !> Point to remove
      integer, intent(in) :: ipt
      !> Copy without the point
      type(fd_cavity), intent(out) :: removed

      !> Points kept
      integer, allocatable :: keep(:)
      !> Grid-point index
      integer :: i

      keep = pack([(i, i=1, cavity%ngrid)], [(i, i=1, cavity%ngrid)] /= ipt)
      removed = cavity
      removed%ngrid = size(keep)
      removed%xyz = cavity%xyz(:, keep)
      removed%normal0 = cavity%normal0(:, keep)
      removed%a = cavity%a(keep)
      removed%k1 = cavity%k1(keep)
      removed%k2 = cavity%k2(keep)
      if (allocated(cavity%xi0)) removed%xi0 = cavity%xi0(keep)
      if (allocated(cavity%f)) removed%f = cavity%f(keep)
   end subroutine remove_point

end module test_model_component_pcm_references
