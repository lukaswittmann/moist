!> Generic checks run on every concrete PCM component
module test_model_component_pcm_general
   use mctc_env, only: wp
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use mstore, only: get_structure
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed, skip_test
   use test_helpers, only: fd4_scalar, fd4_offsets, fd6_scalar, fd6_offsets, check_moist_error, &
      & component_view, copy_potential_adjoint, point_charge_potential, point_charge_field, submit, &
      & get_test_cross
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, &
      & moist_cavity_iswig_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_radii, only: default_cpcm_radii, radius_type, new_radii_custom_atoms
   use moist_channels_coupling, only: coupling_type, gaussian_potential_request_type, &
      & point_potential_request_type, nuclear_charge_request_type
   use moist_channels_response, only: response_type, potential_adjoint_response_type
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, &
      & moist_pcm_parameters_type, solver_type
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, &
      & new_component_cpcm, model_continuum_component_cosmo, new_component_cosmo
   use moist_model_continuum_component_pcm_iefpcm, only: model_continuum_component_iefpcm, &
      & new_component_iefpcm
   use moist_model_continuum_component_pcm_ssvpe, only: model_continuum_component_ssvpe, &
      & new_component_ssvpe
   use moist_model_continuum_component_pcm_dpcm, only: model_continuum_component_dpcm, &
      & new_component_dpcm
   use moist_model_continuum_component_pcm_dpcm_original, only: &
      & model_continuum_component_dpcm_original, new_component_dpcm_original
   use moist_model_continuum_component_pcm_lpb, only: model_continuum_component_lpb, &
      & new_component_lpb
   use moist_model_continuum_component_pcm_slic, only: model_continuum_component_slic, &
      & new_component_slic
   use moist_model_continuum_component_pcm_nonlocal, only: model_continuum_component_nonlocal, &
      & new_component_nonlocal
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_general

   !> Model identifiers of the component factory
   integer, parameter :: pcm_cpcm = 1, pcm_cosmo = 2, pcm_iefpcm = 3, pcm_ssvpe = 4, &
      & pcm_dpcm = 5, pcm_dpcm_original = 6, pcm_lpb = 7, pcm_slic = 8, pcm_nonlocal = 9
   !> Number of models
   integer, parameter :: npcm = 9

   !> Seeds of the synthetic cavity and cavities of the nuclear check
   integer, parameter :: seed_iswig = 1, seed_drop = 2
   !> Seed names for messages
   character(len=5), parameter :: seed_name(2) = ["iSwiG", "DROP "]

   !> Surface channels: positions, normals, curvatures, widths, switching, areas
   integer, parameter :: ch_r = 1, ch_n = 2, ch_k = 3, ch_xi = 4, ch_f = 5, ch_a = 6
   !> Number of surface channels
   integer, parameter :: nchannel = 6
   !> Channel names, as in the expected-channel lists
   character(len=2), parameter :: channel_name(nchannel) = ["r ", "n ", "k ", "xi", "f ", "a "]
   !> Weight channels of a surface accumulator
   character(len=5), parameter :: weight_name(7) = ["w_xyz", "w_n  ", "w_k1 ", "w_k2 ", &
      & "w_xi ", "w_f  ", "w_a  "]
   !> Surface channel of each weight channel
   integer, parameter :: weight_channel(7) = [ch_r, ch_n, ch_k, ch_k, ch_xi, ch_f, ch_a]

   !> Physical parameters of a component, defaults those of the difference checks
   type :: pcm_physics
      !> Static permittivity of every model, bulk permittivity of the nonlocal one
      real(wp) :: eps = 78.4_wp
      !> High-frequency permittivity of the nonlocal model
      real(wp) :: eps_inf = 1.8_wp
      !> Correlation length of the nonlocal model, bohr
      real(wp) :: lambda = 3.0_wp
      !> Inverse Debye length of the LPB model, 1/bohr
      real(wp) :: kappa = 0.3_wp
      !> SLIC h(E_n) = alpha tanh(beta E_n - gamma) + mu, unsaturated at the test fields
      real(wp) :: slic_alpha = 0.3_wp
      !> SLIC beta, 1/(atomic field unit)
      real(wp) :: slic_beta = -2.0_wp
      !> SLIC gamma
      real(wp) :: slic_gamma = -0.4_wp
      !> SLIC mu
      real(wp) :: slic_mu = 0.05_wp
   end type pcm_physics
   !> Newton tolerance of SLIC, far below the finite-difference tolerances
   real(wp), parameter :: newton_tol = 1.0e-13_wp

   !> Lebedev order of the seeding and nuclear-check cavities
   integer, parameter :: nleb = 26
   !> Switching cutoff of the seeds, keeping f well away from zero
   real(wp), parameter :: cut_f = 1.0e-2_wp
   !> Minimum point separation of the point seeds, bohr
   !>
   !> A point surface keeps seam points whatever their f; collocation of
   !> full-area points closer than their own extent is nearly singular, which
   !> spoils the differences (truncation or roundoff) on larger solutes
   real(wp), parameter :: min_sep = 1.0_wp
   !> Step along the surface and trace directions (relative to the field scale)
   !>
   !> With the six-point stencil the differences are roundoff-limited; the
   !> short nonlocal length makes a larger step truncation-limited
   real(wp), parameter :: step = 5.0e-4_wp
   !> Tangent step on the DROP seed, whose points lie closer
   real(wp), parameter :: step_drop = 2.5e-4_wp
   !> Surface-gradient step on the DROP seed (SLIC only)
   !>
   !> Roundoff-limited at about 1e-12 on every step tried; at 1e-3 the worst
   !> case (C5 cross) is 4 (point) to 6 (Gaussian) times below the tolerance
   real(wp), parameter :: step_drop_surface = 1.0e-3_wp
   !> Nuclear-displacement step, bohr
   real(wp), parameter :: step_nuc = 1.0e-4_wp
   !> Energy differences: |an - fd| <= atol + rtol |fd|
   !>
   !> All tolerances are about ten times the largest deviation over every
   !> model and solute (six-point stencil, roundoff-limited)
   real(wp), parameter :: grad_atol = 3.0e-12_wp, grad_rtol = 1.0e-9_wp
   !> Nuclear gradient (four-point stencil), with cavity rebuild roundoff
   !>
   !> The SLIC Newton tolerance sets the largest deviation, 9e-10 relative
   real(wp), parameter :: nuc_atol = 1.0e-12_wp, nuc_rtol = 1.0e-8_wp
   !> Weight differences, per channel relative to its largest entry
   real(wp), parameter :: hess_atol = 3.0e-12_wp, hess_rtol = 1.0e-9_wp
   !> Magnitude below which an expected-active quantity counts as vacuous
   real(wp), parameter :: thr_vacuous = 1.0e-7_wp
   !> Energy derivative relative to |E| below which an expected-active direction is vacuous
   !>
   !> The difference noise is about 1e-12 |E|
   real(wp), parameter :: thr_active = 1.0e-8_wp

   !> Solutes of the difference checks
   integer, parameter :: nmol = 9
   !> Hand-placed bent triatomic without symmetry, the solute of the lifecycle check
   integer, parameter :: mol_triatomic = 1
   !> Five-carbon cross of `get_test_cross`
   integer, parameter :: mol_cross = 9
   !> Solute names for messages
   !>
   !> - triatomic: bent O-H-H, three spheres and no symmetry
   !> - IL16 008B: chloride, a single charged sphere (k1 = k2 on DROP too)
   !> - MB16-43 LiH: two spheres of very different radii
   !> - Amino20x4 GLY_xab: glycine, polar organic
   !> - But14diol 1: intramolecular hydrogen bond, concave crevices
   !> - MB16-43 01: mindless molecule with unusual elements and bonding
   !> - IL16 008A: dimethylimidazolium, net charge +1
   !> - Heavy28 bih3_h2o: bismuth (large radius) and a weakly bound water
   !> - C5 cross: the five-carbon cross of the DROP gradient tests, concave seams
   character(len=*), parameter :: mol_label(nmol) = [character(len=17) :: "triatomic", &
      & "IL16 008B", "MB16-43 LiH", "Amino20x4 GLY_xab", "But14diol 1", "MB16-43 01", &
      & "IL16 008A", "Heavy28 bih3_h2o", "C5 cross"]
   !> mstore collection of each solute (blank for the hand-placed ones)
   character(len=*), parameter :: mol_collection(nmol) = [character(len=9) :: "", "IL16", &
      & "MB16-43", "Amino20x4", "But14diol", "MB16-43", "IL16", "Heavy28", ""]
   !> mstore record of each solute (blank for the hand-placed ones)
   character(len=*), parameter :: mol_record(nmol) = [character(len=8) :: "", "008B", "LiH", &
      & "GLY_xab", "1", "01", "008A", "bih3_h2o", ""]
   !> Atomic numbers of the triatomic
   integer, parameter :: solute_num(3) = [8, 1, 1]
   !> Coordinates of the triatomic, bohr
   real(wp), parameter :: solute_xyz(3, 3) = reshape([0.0_wp, 0.0_wp, 0.1_wp, &
      & 1.45_wp, 1.05_wp, 0.0_wp, -1.40_wp, 1.15_wp, -0.2_wp], [3, 3])
   !> Largest solute whose nuclear coordinates and tangent directions are swept one by one
   !>
   !> Larger solutes get `ndir_nuc` directions over all nuclei and the joint
   !> tangent direction only
   integer, parameter :: max_nat_sweep = 3
   !> Directions of the nuclear check on the larger solutes
   integer, parameter :: ndir_nuc = 3

   !> Born sphere radius, bohr
   real(wp), parameter :: born_radius = 4.0_wp
   !> Lebedev order of the Born sphere
   integer, parameter :: born_nleb = 302
   !> Nonlocal correlation length of the Born check, bohr
   !>
   !> nu = sqrt(eps/eps_inf) R/lambda ~ 1, where the nonlocal energy lies far
   !> from both Born limits and Lambda = 3.8 bohr is resolved by the grid
   real(wp), parameter :: born_lambda = 25.0_wp

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

   !> Synthetic cavity whose surface fields are set directly by the test
   type, extends(cavity_type) :: fd_cavity
   contains
      !> No-op: the fields are set by the test
      procedure :: update => fd_cavity_update
      !> No-op: the synthetic surface has no nuclear derivatives
      procedure :: get_gradient => fd_cavity_get_gradient
   end type fd_cavity

contains

   !> Collect the generic PCM tests, five per model
   !>
   !> @param[out] testsuite Collection of tests
   subroutine collect_model_component_pcm_general(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("cpcm_born", test_cpcm_born), &
         & new_unittest("cpcm_surface_gradient", test_cpcm_surface_gradient), &
         & new_unittest("cpcm_tangent", test_cpcm_tangent), &
         & new_unittest("cpcm_nuclear_gradient", test_cpcm_nuclear_gradient), &
         & new_unittest("cpcm_reference", test_cpcm_reference), &
         & new_unittest("cosmo_born", test_cosmo_born), &
         & new_unittest("cosmo_surface_gradient", test_cosmo_surface_gradient), &
         & new_unittest("cosmo_tangent", test_cosmo_tangent), &
         & new_unittest("cosmo_nuclear_gradient", test_cosmo_nuclear_gradient), &
         & new_unittest("cosmo_reference", test_cosmo_reference), &
         & new_unittest("iefpcm_born", test_iefpcm_born), &
         & new_unittest("iefpcm_surface_gradient", test_iefpcm_surface_gradient), &
         & new_unittest("iefpcm_tangent", test_iefpcm_tangent), &
         & new_unittest("iefpcm_nuclear_gradient", test_iefpcm_nuclear_gradient), &
         & new_unittest("iefpcm_reference", test_iefpcm_reference), &
         & new_unittest("ssvpe_born", test_ssvpe_born), &
         & new_unittest("ssvpe_surface_gradient", test_ssvpe_surface_gradient), &
         & new_unittest("ssvpe_tangent", test_ssvpe_tangent), &
         & new_unittest("ssvpe_nuclear_gradient", test_ssvpe_nuclear_gradient), &
         & new_unittest("ssvpe_reference", test_ssvpe_reference), &
         & new_unittest("dpcm_born", test_dpcm_born), &
         & new_unittest("dpcm_surface_gradient", test_dpcm_surface_gradient), &
         & new_unittest("dpcm_tangent", test_dpcm_tangent), &
         & new_unittest("dpcm_nuclear_gradient", test_dpcm_nuclear_gradient), &
         & new_unittest("dpcm_reference", test_dpcm_reference), &
         & new_unittest("dpcm_original_born", test_dpcm_original_born), &
         & new_unittest("dpcm_original_surface_gradient", test_dpcm_original_surface_gradient), &
         & new_unittest("dpcm_original_tangent", test_dpcm_original_tangent), &
         & new_unittest("dpcm_original_nuclear_gradient", test_dpcm_original_nuclear_gradient), &
         & new_unittest("dpcm_original_reference", test_dpcm_original_reference), &
         & new_unittest("lpb_born", test_lpb_born), &
         & new_unittest("lpb_surface_gradient", test_lpb_surface_gradient), &
         & new_unittest("lpb_tangent", test_lpb_tangent), &
         & new_unittest("lpb_nuclear_gradient", test_lpb_nuclear_gradient), &
         & new_unittest("lpb_reference", test_lpb_reference), &
         & new_unittest("slic_born", test_slic_born), &
         & new_unittest("slic_surface_gradient", test_slic_surface_gradient), &
         & new_unittest("slic_tangent", test_slic_tangent), &
         & new_unittest("slic_nuclear_gradient", test_slic_nuclear_gradient), &
         & new_unittest("slic_reference", test_slic_reference), &
         & new_unittest("nonlocal_born", test_nonlocal_born), &
         & new_unittest("nonlocal_surface_gradient", test_nonlocal_surface_gradient), &
         & new_unittest("nonlocal_tangent", test_nonlocal_tangent), &
         & new_unittest("nonlocal_nuclear_gradient", test_nonlocal_nuclear_gradient), &
         & new_unittest("nonlocal_reference", test_nonlocal_reference), &
         & new_unittest("tangent_requires_solution", test_tangent_requires_solution) &
         & ]
   end subroutine collect_model_component_pcm_general

   !* ================================================================================= *!
   !*                                 Per-model entries                                 *!
   !* ================================================================================= *!

   !> CPCM: Born ion
   !>
   !> @param[out] error Test failure
   subroutine test_cpcm_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_cpcm)
   end subroutine test_cpcm_born

   !> CPCM: surface and trace gradients
   !>
   !> @param[out] error Test failure
   subroutine test_cpcm_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_cpcm)
   end subroutine test_cpcm_surface_gradient

   !> CPCM: surface weight tangent
   !>
   !> @param[out] error Test failure
   subroutine test_cpcm_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_cpcm)
   end subroutine test_cpcm_tangent

   !> CPCM: nuclear gradient
   !>
   !> @param[out] error Test failure
   subroutine test_cpcm_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_cpcm)
   end subroutine test_cpcm_nuclear_gradient

   !> CPCM: reference energies of other codes and the literature
   !>
   !> @param[out] error Test failure
   subroutine test_cpcm_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_cpcm)
   end subroutine test_cpcm_reference

   !> COSMO: Born ion
   !>
   !> @param[out] error Test failure
   subroutine test_cosmo_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_cosmo)
   end subroutine test_cosmo_born

   !> COSMO: surface and trace gradients
   !>
   !> @param[out] error Test failure
   subroutine test_cosmo_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_cosmo)
   end subroutine test_cosmo_surface_gradient

   !> COSMO: surface weight tangent
   !>
   !> @param[out] error Test failure
   subroutine test_cosmo_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_cosmo)
   end subroutine test_cosmo_tangent

   !> COSMO: nuclear gradient
   !>
   !> @param[out] error Test failure
   subroutine test_cosmo_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_cosmo)
   end subroutine test_cosmo_nuclear_gradient

   !> COSMO: reference energies of other codes and the literature
   !>
   !> @param[out] error Test failure
   subroutine test_cosmo_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_cosmo)
   end subroutine test_cosmo_reference

   !> IEF-PCM: Born ion
   !>
   !> @param[out] error Test failure
   subroutine test_iefpcm_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_iefpcm)
   end subroutine test_iefpcm_born

   !> IEF-PCM: surface and trace gradients
   !>
   !> @param[out] error Test failure
   subroutine test_iefpcm_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_iefpcm)
   end subroutine test_iefpcm_surface_gradient

   !> IEF-PCM: surface weight tangent
   !>
   !> @param[out] error Test failure
   subroutine test_iefpcm_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_iefpcm)
   end subroutine test_iefpcm_tangent

   !> IEF-PCM: nuclear gradient
   !>
   !> @param[out] error Test failure
   subroutine test_iefpcm_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_iefpcm)
   end subroutine test_iefpcm_nuclear_gradient

   !> IEF-PCM: reference energies of other codes and the literature
   !>
   !> @param[out] error Test failure
   subroutine test_iefpcm_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_iefpcm)
   end subroutine test_iefpcm_reference

   !> SS(V)PE: Born ion
   !>
   !> @param[out] error Test failure
   subroutine test_ssvpe_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_ssvpe)
   end subroutine test_ssvpe_born

   !> SS(V)PE: surface and trace gradients
   !>
   !> @param[out] error Test failure
   subroutine test_ssvpe_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_ssvpe)
   end subroutine test_ssvpe_surface_gradient

   !> SS(V)PE: surface weight tangent
   !>
   !> @param[out] error Test failure
   subroutine test_ssvpe_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_ssvpe)
   end subroutine test_ssvpe_tangent

   !> SS(V)PE: nuclear gradient
   !>
   !> @param[out] error Test failure
   subroutine test_ssvpe_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_ssvpe)
   end subroutine test_ssvpe_nuclear_gradient

   !> SS(V)PE: reference energies of other codes and the literature
   !>
   !> @param[out] error Test failure
   subroutine test_ssvpe_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_ssvpe)
   end subroutine test_ssvpe_reference

   !> DPCM: Born ion
   !>
   !> @param[out] error Test failure
   subroutine test_dpcm_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_dpcm)
   end subroutine test_dpcm_born

   !> DPCM: surface and trace gradients
   !>
   !> @param[out] error Test failure
   subroutine test_dpcm_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_dpcm)
   end subroutine test_dpcm_surface_gradient

   !> DPCM: surface weight tangent
   !>
   !> @param[out] error Test failure
   subroutine test_dpcm_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_dpcm)
   end subroutine test_dpcm_tangent

   !> DPCM: nuclear gradient
   !>
   !> @param[out] error Test failure
   subroutine test_dpcm_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_dpcm)
   end subroutine test_dpcm_nuclear_gradient

   !> DPCM: reference energies of other codes and the literature
   !>
   !> @param[out] error Test failure
   subroutine test_dpcm_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_dpcm)
   end subroutine test_dpcm_reference

   !> Original DPCM: Born ion
   !>
   !> @param[out] error Test failure
   subroutine test_dpcm_original_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_born

   !> Original DPCM: surface and trace gradients
   !>
   !> @param[out] error Test failure
   subroutine test_dpcm_original_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_surface_gradient

   !> Original DPCM: surface weight tangent
   !>
   !> @param[out] error Test failure
   subroutine test_dpcm_original_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_tangent

   !> Original DPCM: nuclear gradient
   !>
   !> @param[out] error Test failure
   subroutine test_dpcm_original_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_nuclear_gradient

   !> Original DPCM: reference energies of other codes and the literature
   !>
   !> @param[out] error Test failure
   subroutine test_dpcm_original_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_reference

   !> LPB: Born ion
   !>
   !> @param[out] error Test failure
   subroutine test_lpb_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_lpb)
   end subroutine test_lpb_born

   !> LPB: surface and trace gradients
   !>
   !> @param[out] error Test failure
   subroutine test_lpb_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_lpb)
   end subroutine test_lpb_surface_gradient

   !> LPB: surface weight tangent
   !>
   !> @param[out] error Test failure
   subroutine test_lpb_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_lpb)
   end subroutine test_lpb_tangent

   !> LPB: nuclear gradient
   !>
   !> @param[out] error Test failure
   subroutine test_lpb_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_lpb)
   end subroutine test_lpb_nuclear_gradient

   !> LPB: reference energies of other codes and the literature
   !>
   !> @param[out] error Test failure
   subroutine test_lpb_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_lpb)
   end subroutine test_lpb_reference

   !> SLIC: Born ion
   !>
   !> @param[out] error Test failure
   subroutine test_slic_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_slic)
   end subroutine test_slic_born

   !> SLIC: surface and trace gradients
   !>
   !> @param[out] error Test failure
   subroutine test_slic_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_slic)
   end subroutine test_slic_surface_gradient

   !> SLIC: surface weight tangent
   !>
   !> @param[out] error Test failure
   subroutine test_slic_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_slic)
   end subroutine test_slic_tangent

   !> SLIC: nuclear gradient
   !>
   !> @param[out] error Test failure
   subroutine test_slic_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_slic)
   end subroutine test_slic_nuclear_gradient

   !> SLIC: reference energies of other codes and the literature
   !>
   !> @param[out] error Test failure
   subroutine test_slic_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_slic)
   end subroutine test_slic_reference

   !> Nonlocal: Born ion
   !>
   !> @param[out] error Test failure
   subroutine test_nonlocal_born(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_born(error, pcm_nonlocal)
   end subroutine test_nonlocal_born

   !> Nonlocal: surface and trace gradients
   !>
   !> @param[out] error Test failure
   subroutine test_nonlocal_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_nonlocal)
   end subroutine test_nonlocal_surface_gradient

   !> Nonlocal: surface weight tangent
   !>
   !> @param[out] error Test failure
   subroutine test_nonlocal_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_nonlocal)
   end subroutine test_nonlocal_tangent

   !> Nonlocal: nuclear gradient
   !>
   !> @param[out] error Test failure
   subroutine test_nonlocal_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_nonlocal)
   end subroutine test_nonlocal_nuclear_gradient

   !> Nonlocal: reference energies of other codes and the literature
   !>
   !> @param[out] error Test failure
   subroutine test_nonlocal_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_reference(error, pcm_nonlocal)
   end subroutine test_nonlocal_reference

   !> The tangent needs current charges and nothing else, in every model
   !>
   !> - After `update` the operators are new and no charges exist yet: refused
   !> - Directly after the energy, with no surface weights requested: accepted
   !> - After a second `update` the charges belong to the previous operators: refused
   !>
   !> @param[out] error Test failure
   subroutine test_tangent_requires_solution(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the components
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Solute and synthetic cavity
      type(structure_type) :: mol
      type(fd_cavity) :: cavity
      !> Zero direction and tangent accumulator
      type(cavity_surface_adjoint_type) :: direction, tangent
      !> Potential direction, adjoint tangent and energy
      real(wp), allocatable :: dphi(:), dw_phi(:)
      real(wp) :: energy
      !> Model index
      integer :: kind

      call new_context(ctx, verbosity=0)
      call get_solute(mol_triatomic, mol)
      call get_seed(seed_iswig, .true., mol, cavity, error)
      if (allocated(error)) return
      allocate (dphi(cavity%ngrid), source=0.0_wp)
      do kind = 1, npcm
         call new_pcm(kind, ctx, pcm, error)
         if (allocated(error)) return
         call pcm%update(mol, cavity, err)
         call check_moist_error(error, err, "update")
         if (allocated(error)) return
         call tangent%init(cavity%ngrid)
         call pcm%surface_weight_tangent(cavity, direction, dphi, tangent, dw_phi, err)
         call check(error, allocated(err), &
            & message=trim(pcm%name)//": tangent accepted operators without charges")
         if (allocated(error)) return
         call host_energy(pcm, mol, cavity, mol, host_charges(mol), energy, error)
         if (allocated(error)) return
         call pcm%surface_weight_tangent(cavity, direction, dphi, tangent, dw_phi, err)
         call check_moist_error(error, err, trim(pcm%name)//": tangent directly after the energy")
         if (allocated(error)) return
         call pcm%update(mol, cavity, err)
         call check_moist_error(error, err, "update")
         if (allocated(error)) return
         call pcm%surface_weight_tangent(cavity, direction, dphi, tangent, dw_phi, err)
         call check(error, allocated(err), &
            & message=trim(pcm%name)//": tangent accepted charges of stale operators")
         if (allocated(error)) return
         deallocate (err)
      end do
   end subroutine test_tangent_requires_solution

   !* ================================================================================= *!
   !*                                      Drivers                                      *!
   !* ================================================================================= *!

   !> Surface and trace gradients on Gaussian and point surfaces of every solute
   !>
   !> @param[out] error Test failure
   !> @param[in]  kind  Model identifier
   subroutine run_surface_gradient(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind
      !> Solute index
      integer :: imol

      do imol = 1, nmol
         call check_surface_gradient(error, kind, imol, .true.)
         if (allocated(error)) return
         call check_surface_gradient(error, kind, imol, .false.)
         if (allocated(error)) return
      end do
   end subroutine run_surface_gradient

   !> Surface weight tangent on both seeds (SLIC: DROP only), Gaussian and point surfaces of every solute
   !>
   !> @param[out] error Test failure
   !> @param[in]  kind  Model identifier
   subroutine run_tangent(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind
      !> Solute and seed indices
      integer :: imol, seed

      do imol = 1, nmol
         do seed = first_seed(kind), seed_drop
            call check_tangent(error, kind, imol, seed, .true.)
            if (allocated(error)) return
            call check_tangent(error, kind, imol, seed, .false.)
            if (allocated(error)) return
         end do
      end do
   end subroutine run_tangent

   !> Nuclear gradient on rebuilt iSwiG and SvdW-DROP cavities of every solute (SLIC: DROP only)
   !>
   !> @param[out] error Test failure
   !> @param[in]  kind  Model identifier
   subroutine run_nuclear_gradient(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind
      !> Solute and seed indices
      integer :: imol, seed

      do imol = 1, nmol
         do seed = first_seed(kind), seed_drop
            call check_nuclear_gradient(error, kind, imol, seed)
            if (allocated(error)) return
         end do
      end do
   end subroutine run_nuclear_gradient

   !> Ion of either sign at the center of a spherical cavity against the analytic energy
   !>
   !> The host potential is that of the central charge, Gaussian-probed on the
   !> iSwiG surface (SLIC: SvdW-DROP, see `first_seed`); the even models give
   !> the same energy for both signs, SLIC solvates the anion more strongly for
   !> beta < 0
   !>
   !> @param[out] error Test failure
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
      !> Physical parameters, nonlocal correlation length of the Born check
      type(pcm_physics) :: physics
      !> Charge, energy, reference and relative deviation per sign
      real(wp) :: charge, energy(2), reference, deviation
      !> Sign index
      integer :: isign
      !> Failure context
      character(len=160) :: context

      physics%lambda = born_lambda
      call new_context(ctx, verbosity=0)
      call new_pcm(kind, ctx, pcm, error, physics)
      if (allocated(error)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_sphere_cavity(kind, [born_radius], born_nleb, ctx, cavity, error)
      if (allocated(error)) return
      call pcm%configure_cavity(cavity, err)
      call check_moist_error(error, err, "cavity configuration")
      if (allocated(error)) return
      call cavity%update(mol, error=err)
      call check_moist_error(error, err, "cavity update")
      if (allocated(error)) return

      do isign = 1, 2
         charge = merge(1.0_wp, -1.0_wp, isign == 1)
         mol%charge = charge
         call host_energy(pcm, mol, cavity, mol, [charge], energy(isign), error)
         if (allocated(error)) return
         reference = born_reference(kind, physics, charge, born_radius)
         deviation = abs(energy(isign) - reference)/abs(reference)
         write (context, '(a,": charge ",f4.1," energy ",es16.8," reference ",es16.8, &
            & ", rel. dev ",es9.2)') trim(pcm%name), charge, energy(isign), reference, deviation
         call check(error, deviation <= born_tolerance(kind), &
            & message="Born energy deviates from the analytic result", more=trim(context))
         if (allocated(error)) return
      end do
      if (kind == pcm_slic) then
         call check(error, energy(2) < energy(1), &
            & message="SLIC must solvate the anion more strongly for beta < 0")
      else
         call check(error, abs(energy(2) - energy(1)) <= 1.0e-12_wp*abs(energy(1)), &
            & message=trim(pcm%name)//": Born energy is not even in the charge")
      end if
   end subroutine check_born

   !> Analytic energy of a central charge in a sphere
   !>
   !> - Conductor: -feps q**2/(2 R) with the model's dielectric factor
   !> - IEF-PCM, SS(V)PE and both DPCM forms: Born, feps = (eps - 1)/eps
   !> - LPB: Debye-Hueckel, feps = 1 - 1/(eps (1 + kappa R))
   !> - SLIC: feps = (eps - 1)/(eps - (eps - 1) h(q/R**2))
   !> - Nonlocal: exact monopole solution of the local reformulation with
   !>   Hildebrandt's approximate interface conditions (Bardhan, Knepley,
   !>   Brune, arXiv:1208.3866, eqs. 21-26), solute permittivity 1,
   !>   feps = 1 - 1/eps - (1/eps_inf - 1/eps)/(1 + nu), nu = sqrt(eps/eps_inf) R/lambda
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
            feps = (eps - 1.0_wp)/(eps + 0.5_wp)
         case (pcm_lpb)
            feps = 1.0_wp - 1.0_wp/(eps*(1.0_wp + physics%kappa*radius))
         case (pcm_slic)
            h = physics%slic_alpha*tanh(physics%slic_beta*charge/radius**2 - physics%slic_gamma) &
               & + physics%slic_mu
            feps = (eps - 1.0_wp)/(eps - (eps - 1.0_wp)*h)
         case (pcm_nonlocal)
            nu = sqrt(eps/physics%eps_inf)*radius/physics%lambda
            feps = 1.0_wp - 1.0_wp/eps - (1.0_wp/physics%eps_inf - 1.0_wp/eps)/(1.0_wp + nu)
         case default
            feps = (eps - 1.0_wp)/eps
         end select
      end associate
      energy = -0.5_wp*feps*charge**2/radius
   end function born_reference

   !> Relative tolerance of the Born check on the 302-point Gaussian sphere
   !>
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
   !> @param[in] kind Model identifier
   pure function born_tolerance(kind) result(tol)
      !> Model identifier
      integer, intent(in) :: kind
      !> Relative tolerance
      real(wp) :: tol

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

   !> Energies of the reference systems against the rows of the reference table
   !>
   !> Every row of the model is evaluated. A filled row must agree within its
   !> tolerance; pending rows skip the test and list the moist energies
   !>
   !> @param[out] error Test failure
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

   !> Surface weights and trace adjoints against energy differences
   !>
   !> The host is fixed point charges on the solute atoms, re-evaluated on the
   !> displaced surface, so a surface step moves the traces as well and the
   !> difference is the total derivative `get_gradient_surface_weights` returns
   !>
   !> @param[out] error   Test failure
   !> @param[in]  kind    Model identifier
   !> @param[in]  imol    Solute index
   !> @param[in]  smeared Gaussian surface (else point surface)
   subroutine check_surface_gradient(error, kind, imol, smeared)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind
      !> Solute index
      integer, intent(in) :: imol
      !> Gaussian surface
      logical, intent(in) :: smeared

      !> Run context borrowed by the component
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Solute
      type(structure_type) :: mol
      !> Synthetic cavity and its displacement
      type(fd_cavity) :: cavity, trial
      !> Total surface weights and a one-channel direction
      type(cavity_surface_adjoint_type) :: acc, direction
      !> Host charges on the solute atoms
      real(wp), allocatable :: q(:)
      !> Host traces, their directions and adjoints
      real(wp), allocatable :: phi(:), dphi_dr(:, :), du(:), dg(:, :), w_phi(:), w_dphi_dr(:, :)
      !> Reference energy, stencil energies, analytic and FD directional derivatives
      real(wp) :: energy, values(6), analytic, fd
      !> Difference step of the seed
      real(wp) :: h
      !> Channel and stencil indices, cavity seed
      integer :: ich, k, seed
      !> Failure context
      character(len=200) :: context, prefix

      call new_context(ctx, verbosity=0)
      call new_pcm(kind, ctx, pcm, error)
      if (allocated(error)) return
      call get_solute(imol, mol)
      q = host_charges(mol)
      seed = first_seed(kind)
      call get_seed(seed, smeared, mol, cavity, error)
      if (allocated(error)) return
      h = merge(step_drop_surface, step, seed == seed_drop)
      prefix = trim(pcm%name)//", "//trim(mol_label(imol))//", "//trim(seed_name(seed))// &
         & merge(", Gaussian", ", point   ", smeared)
      call host_surface_weights(pcm, mol, cavity, q, acc, energy, error)
      if (allocated(error)) return

      do ich = 1, nchannel
         if (.not. provides(cavity, ich)) cycle
         call channel_direction(cavity, ich, direction)
         analytic = contract(acc, direction)
         do k = 1, 6
            call displace(cavity, direction, fd6_offsets(k)*h, trial)
            call host_energy(pcm, mol, trial, mol, q, values(k), error)
            if (allocated(error)) return
         end do
         call fd6(values, h, fd, error)
         if (allocated(error)) return
         write (context, '(a,": channel ",a," dE/dt ",es12.4," vs FD ",es12.4,", dev ",es9.2)') &
            & trim(prefix), trim(channel_name(ich)), analytic, fd, abs(analytic - fd)
         call check(error, abs(analytic - fd) <= grad_atol + grad_rtol*abs(fd), &
            & message="surface weights disagree with the energy difference", more=trim(context))
         if (allocated(error)) return
         if (expects(active_channels(kind, smeared), ich)) then
            call check(error, abs(fd) > thr_active*abs(energy), message="vacuous channel", &
               & more=trim(context))
            if (allocated(error)) return
         end if
      end do

      call host_traces(mol, q, cavity, smeared, phi, dphi_dr)
      call trace_adjoints(pcm, mol, cavity, phi, dphi_dr, w_phi, w_dphi_dr, error)
      if (allocated(error)) return
      call trace_direction(cavity%ngrid, 1, du)
      call field_direction(cavity%ngrid, dg)

      ! potential trace
      analytic = dot_product(w_phi, du)
      do k = 1, 6
         call trace_energy(pcm, mol, cavity, phi + fd6_offsets(k)*h*du, dphi_dr, values(k), error)
         if (allocated(error)) return
      end do
      call fd6(values, h, fd, error)
      if (allocated(error)) return
      write (context, '(a,": phi dE/dt ",es12.4," vs FD ",es12.4,", dev ",es9.2)') &
         & trim(prefix), analytic, fd, abs(analytic - fd)
      call check(error, abs(analytic - fd) <= grad_atol + grad_rtol*abs(fd) .and. &
         & abs(fd) > thr_active*abs(energy), message="dE/dphi disagrees with the energy difference", &
         & more=trim(context))
      if (allocated(error)) return

      ! spatial derivative of the potential, read through n . dphi_dr
      analytic = 0.0_wp
      if (allocated(w_dphi_dr)) analytic = sum(w_dphi_dr*dg)
      do k = 1, 6
         call trace_energy(pcm, mol, cavity, phi, dphi_dr + fd6_offsets(k)*h*dg, values(k), error)
         if (allocated(error)) return
      end do
      call fd6(values, h, fd, error)
      if (allocated(error)) return
      write (context, '(a,": dphi_dr dE/dt ",es12.4," vs FD ",es12.4,", dev ",es9.2)') &
         & trim(prefix), analytic, fd, abs(analytic - fd)
      call check(error, abs(analytic - fd) <= grad_atol + grad_rtol*abs(fd), &
         & message="dE/d(dphi_dr) disagrees with the energy difference", more=trim(context))
      if (allocated(error)) return
      call check(error, allocated(w_dphi_dr) .eqv. pcm%uses_normal_trace(), &
         & message="dE/d(dphi_dr) is returned exactly for normal-trace models", more=trim(context))
      if (allocated(error)) return
      if (pcm%uses_normal_trace()) then
         call check(error, abs(fd) > thr_active*abs(energy), message="vacuous normal trace", &
            & more=trim(context))
      end if
   end subroutine check_surface_gradient

   !> Surface weight tangent against differences of the weights at fixed traces
   !>
   !> Directions: all channels and traces jointly; on solutes up to
   !> `max_nat_sweep` atoms also each surface channel the cavity provides alone
   !> and each trace alone. The weights are those of `get_surface_weights` on
   !> a field-independent cavity (operator part only), the adjoints those of
   !> `potential_adjoint` and `normal_trace_adjoint`
   !>
   !> @param[out] error   Test failure
   !> @param[in]  kind    Model identifier
   !> @param[in]  imol    Solute index
   !> @param[in]  seed    Seed of the synthetic cavity
   !> @param[in]  smeared Gaussian surface (else point surface)
   subroutine check_tangent(error, kind, imol, seed, smeared)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind
      !> Solute index
      integer, intent(in) :: imol
      !> Seed of the synthetic cavity
      integer, intent(in) :: seed
      !> Gaussian surface
      logical, intent(in) :: smeared

      !> Run context borrowed by the component
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Solute
      type(structure_type) :: mol
      !> Synthetic cavity
      type(fd_cavity) :: cavity
      !> Direction of the surface fields
      type(cavity_surface_adjoint_type) :: direction
      !> Fixed traces, their directions and a zero direction
      real(wp), allocatable :: phi(:), dphi_dr(:, :), vn(:), du(:), dq(:), zero(:)
      !> Step of the differences
      real(wp) :: h
      !> Channel and grid-point indices
      integer :: ich, i
      !> Prefix of the failure context
      character(len=128) :: prefix

      call new_context(ctx, verbosity=0)
      call new_pcm(kind, ctx, pcm, error)
      if (allocated(error)) return
      call get_solute(imol, mol)
      call get_seed(seed, smeared, mol, cavity, error)
      if (allocated(error)) return
      h = merge(step_drop, step, seed == seed_drop)
      prefix = trim(pcm%name)//", "//trim(mol_label(imol))//", "//trim(seed_name(seed))// &
         & merge(" Gaussian", " point   ", smeared)
      call host_traces(mol, host_charges(mol), cavity, smeared, phi, dphi_dr)
      allocate (vn(cavity%ngrid))
      do i = 1, cavity%ngrid
         vn(i) = dot_product(cavity%normal0(:, i), dphi_dr(:, i))
      end do
      call trace_direction(cavity%ngrid, 1, du)
      call trace_direction(cavity%ngrid, 2, dq)
      allocate (zero(cavity%ngrid), source=0.0_wp)

      if (mol%nat <= max_nat_sweep) then
         do ich = 1, nchannel
            if (.not. provides(cavity, ich)) cycle
            call channel_direction(cavity, ich, direction)
            call tangent_case(error, pcm, mol, cavity, phi, vn, direction, zero, zero, h, &
               & trim(prefix)//", channel "//trim(channel_name(ich)), "")
            if (allocated(error)) return
         end do

         call clear_direction(direction)
         call tangent_case(error, pcm, mol, cavity, phi, vn, direction, du, zero, h, &
            & trim(prefix)//", phi", "")
         if (allocated(error)) return
         if (pcm%uses_normal_trace()) then
            call tangent_case(error, pcm, mol, cavity, phi, vn, direction, zero, dq, h, &
               & trim(prefix)//", vn", "")
            if (allocated(error)) return
         end if
      end if

      call clear_direction(direction)
      do ich = 1, nchannel
         if (provides(cavity, ich)) call add_channel_direction(cavity, ich, direction)
      end do
      call tangent_case(error, pcm, mol, cavity, phi, vn, direction, du, dq, h, &
         & trim(prefix)//", joint", active_channels(kind, smeared))
   end subroutine check_tangent

   !> One direction of the tangent check
   !>
   !> @param[out]   error     Test failure
   !> @param[in,out] pcm      Component under test
   !> @param[in]    mol       Solute
   !> @param[in,out] cavity   Reference cavity
   !> @param[in]    phi       Fixed potential trace (ngrid)
   !> @param[in]    vn        Fixed normal trace (ngrid)
   !> @param[in]    direction Direction of the surface fields
   !> @param[in]    du        Direction of the potential trace (ngrid)
   !> @param[in]    dq        Direction of the normal trace (ngrid)
   !> @param[in]    h         Step of the differences
   !> @param[in]    prefix    Failure context prefix
   !> @param[in]    channels  Channels whose tangent must not vanish (blank for none)
   subroutine tangent_case(error, pcm, mol, cavity, phi, vn, direction, du, dq, h, prefix, channels)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Reference cavity
      type(fd_cavity), intent(inout) :: cavity
      !> Fixed potential trace
      real(wp), intent(in) :: phi(:)
      !> Fixed normal trace
      real(wp), intent(in) :: vn(:)
      !> Direction of the surface fields
      type(cavity_surface_adjoint_type), intent(in) :: direction
      !> Direction of the potential trace
      real(wp), intent(in) :: du(:)
      !> Direction of the normal trace
      real(wp), intent(in) :: dq(:)
      !> Step of the differences
      real(wp), intent(in) :: h
      !> Failure context prefix
      character(len=*), intent(in) :: prefix
      !> Channels whose tangent must not vanish
      character(len=*), intent(in) :: channels

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Displaced cavity
      type(fd_cavity) :: trial
      !> Weights and adjoints at the reference point, adjoint tangents
      real(wp), allocatable :: w_phi(:), w_vn(:), dw_phi(:), dw_vn(:)
      !> Adjoints on the stencil
      real(wp), allocatable :: s_phi(:, :), s_vn(:, :)
      !> Tangent accumulator, reference and stencil weights
      type(cavity_surface_adjoint_type) :: tangent, acc0, acc(6)
      !> Stencil index, channel index and number of grid points
      integer :: k, ich, ngrid
      !> Whether the expected channels must be active
      logical :: active

      ngrid = cavity%ngrid
      active = len_trim(channels) > 0
      call fixed_trace_state(pcm, mol, cavity, phi, vn, acc0, w_phi, w_vn, error)
      if (allocated(error)) return
      call tangent%init(ngrid)
      call pcm%surface_weight_tangent(cavity, direction, du, tangent, dw_phi, err, &
         & dvn=dq, dw_vn=dw_vn)
      call check_moist_error(error, err, prefix//": tangent")
      if (allocated(error)) return
      call check(error, allocated(dw_phi), message=prefix//": dw_phi is missing")
      if (allocated(error)) return
      call check(error, allocated(dw_vn) .eqv. pcm%uses_normal_trace(), &
         & message=prefix//": dw_vn is returned exactly for normal-trace models")
      if (allocated(error)) return

      allocate (s_phi(ngrid, 6), s_vn(ngrid, 6))
      do k = 1, 6
         call displace(cavity, direction, fd6_offsets(k)*h, trial)
         call fixed_trace_state(pcm, mol, trial, phi + fd6_offsets(k)*h*du, &
            & vn + fd6_offsets(k)*h*dq, acc(k), w_phi, w_vn, error)
         if (allocated(error)) return
         s_phi(:, k) = w_phi
         s_vn(:, k) = w_vn
      end do

      do ich = 1, size(weight_name)
         call compare_channel(error, prefix, trim(weight_name(ich)), &
            & weights_of(tangent, weight_name(ich)), stencil_weights(acc, weight_name(ich)), h, &
            & active .and. expects(channels, weight_channel(ich)))
         if (allocated(error)) return
      end do
      call compare_channel(error, prefix, "dE/dphi", dw_phi, s_phi, h, active)
      if (allocated(error)) return
      if (pcm%uses_normal_trace()) then
         call compare_channel(error, prefix, "dE/dvn", dw_vn, s_vn, h, active)
      end if
   end subroutine tangent_case

   !> Compare an analytic tangent entrywise with the six-point difference
   !>
   !> Passes when max |an - fd| <= atol + rtol max |fd| over the channel
   !>
   !> @param[out] error    Test failure
   !> @param[in]  prefix   Failure context prefix
   !> @param[in]  name     Channel name
   !> @param[in]  analytic Analytic tangent
   !> @param[in]  values   Values on the stencil (n, 6), in `fd6_offsets` order
   !> @param[in]  h        Step of the differences
   !> @param[in]  active   The tangent must not vanish
   subroutine compare_channel(error, prefix, name, analytic, values, h, active)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Failure context prefix
      character(len=*), intent(in) :: prefix
      !> Channel name
      character(len=*), intent(in) :: name
      !> Analytic tangent
      real(wp), intent(in) :: analytic(:)
      !> Values on the stencil
      real(wp), intent(in) :: values(:, :)
      !> Step of the differences
      real(wp), intent(in) :: h
      !> The tangent must not vanish
      logical, intent(in) :: active

      !> Entrywise six-point differences
      real(wp), allocatable :: fd(:)
      !> Largest deviation and largest difference
      real(wp) :: dev, scale
      !> Entry index
      integer :: i
      !> Failure context
      character(len=200) :: context

      allocate (fd(size(analytic)))
      do i = 1, size(analytic)
         call fd6(values(i, :), h, fd(i), error)
         if (allocated(error)) return
      end do
      dev = maxval(abs(analytic - fd))
      scale = maxval(abs(fd))
      write (context, '(a,", ",a,": max dev ",es10.2," scale ",es10.2)') prefix, name, dev, scale
      call check(error, dev <= hess_atol + hess_rtol*scale, &
         & message="tangent disagrees with the weight difference", more=trim(context))
      if (allocated(error)) return
      if (active) then
         call check(error, maxval(abs(analytic)) > thr_vacuous, message="vacuous tangent", &
            & more=trim(context))
      end if
   end subroutine compare_channel

   !> Six-point central difference of values in `fd6_offsets` order
   !>
   !> @param[in]  values Values on the stencil (6)
   !> @param[in]  h      Step
   !> @param[out] df     Derivative
   !> @param[out] error  Test failure on a nonfinite derivative
   subroutine fd6(values, h, df, error)
      !> Values on the stencil
      real(wp), intent(in) :: values(:)
      !> Step
      real(wp), intent(in) :: h
      !> Derivative
      real(wp), intent(out) :: df
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call fd6_scalar(values(1), values(2), values(3), values(4), values(5), values(6), h, df, error)
   end subroutine fd6

   !> One weight channel of an accumulator, flattened
   !>
   !> @param[in] acc  Surface weights
   !> @param[in] name Weight channel, one of `weight_name`
   function weights_of(acc, name) result(w)
      !> Surface weights
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Weight channel
      character(len=*), intent(in) :: name
      !> Flattened weights
      real(wp), allocatable :: w(:)

      select case (trim(name))
      case ("w_xyz")
         w = reshape(acc%w_xyz, [size(acc%w_xyz)])
      case ("w_n")
         w = reshape(acc%w_n, [size(acc%w_n)])
      case ("w_k1")
         w = acc%w_k1
      case ("w_k2")
         w = acc%w_k2
      case ("w_xi")
         w = acc%w_xi
      case ("w_f")
         w = acc%w_f
      case ("w_a")
         w = acc%w_a
      case default
         error stop "weights_of: unknown weight channel"
      end select
   end function weights_of

   !> One weight channel of the stencil accumulators, (n, 6)
   !>
   !> @param[in] acc  Surface weights on the stencil (6)
   !> @param[in] name Weight channel, one of `weight_name`
   function stencil_weights(acc, name) result(values)
      !> Surface weights on the stencil
      type(cavity_surface_adjoint_type), intent(in) :: acc(:)
      !> Weight channel
      character(len=*), intent(in) :: name
      !> Flattened weights per stencil point
      real(wp), allocatable :: values(:, :)
      !> Stencil index
      integer :: k

      allocate (values(size(weights_of(acc(1), name)), size(acc)))
      do k = 1, size(acc)
         values(:, k) = weights_of(acc(k), name)
      end do
   end function stencil_weights

   !> Nuclear gradient against energy differences with the cavity rebuilt
   !>
   !> The host is the solute nuclei (charges Z, answered as the nuclear charges)
   !> plus electron charges fixed in space. The electrons do not move with the
   !> nuclei, so the energy difference over a nuclear displacement is the
   !> complete gradient. Small solutes difference every coordinate, larger ones
   !> move all nuclei at once along a few fixed directions
   !>
   !> @param[out] error Test failure
   !> @param[in]  kind  Model identifier
   !> @param[in]  imol  Solute index
   !> @param[in]  seed  Cavity kind rebuilt at every geometry
   subroutine check_nuclear_gradient(error, kind, imol, seed)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind
      !> Solute index
      integer, intent(in) :: imol
      !> Cavity kind
      integer, intent(in) :: seed

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the cavity and the component
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Solute and displaced solute
      type(structure_type) :: mol, trial
      !> Molecular cavity, rebuilt at every geometry
      class(cavity_type), allocatable :: cavity
      !> Host coupling of the reference geometry
      type(coupling_type), target :: coupling
      !> Host part of the gradient phase, unused
      type(response_type) :: response
      !> Fixed electron positions (3, nat) and charges (nat)
      real(wp), allocatable :: el_xyz(:, :), el_q(:)
      !> Analytic nuclear gradient (3, nat) and displacement directions (3, nat, ndir)
      real(wp), allocatable :: gradient(:, :), dirs(:, :, :)
      !> Energy, stencil energies, analytic and FD derivatives, scale of the terms
      real(wp) :: energy, values(4), analytic, fd, scale
      !> Direction and stencil indices
      integer :: idir, k
      !> Failure context
      character(len=240) :: context

      call new_context(ctx, verbosity=0)
      call new_pcm(kind, ctx, pcm, error)
      if (allocated(error)) return
      call get_solute(imol, mol)
      call electron_host(mol, el_xyz, el_q)
      call new_host_cavity(seed, pcm, ctx, cavity, error)
      if (allocated(error)) return

      call nuclear_state(pcm, mol, el_xyz, el_q, cavity, coupling, energy, error)
      if (allocated(error)) return
      call pcm%prepare_gradient(cavity, coupling, err)
      call check_moist_error(error, err, "gradient staging")
      if (allocated(error)) return
      call answer_nuclear_host(cavity, coupling, mol, el_xyz, el_q)
      allocate (gradient(3, mol%nat), source=0.0_wp)
      call pcm%get_gradient(component_view(coupling), cavity, response, gradient, err)
      call check_moist_error(error, err, trim(pcm%name)//", "//trim(mol_label(imol))//" on "// &
         & trim(seed_name(seed))//": nuclear gradient")
      if (allocated(error)) return
      call check(error, maxval(abs(gradient)) > thr_vacuous, message=trim(pcm%name)//", "// &
         & trim(mol_label(imol))//": nuclear gradient is vacuously zero")
      if (allocated(error)) return

      call nuclear_directions(mol%nat, dirs)
      do idir = 1, size(dirs, 3)
         do k = 1, 4
            trial = mol
            trial%xyz = mol%xyz + fd4_offsets(k)*step_nuc*dirs(:, :, idir)
            call nuclear_state(pcm, trial, el_xyz, el_q, cavity, coupling, values(k), error)
            if (allocated(error)) return
         end do
         call fd4_scalar(values(1), values(2), values(3), values(4), step_nuc, fd, error)
         if (allocated(error)) return
         analytic = sum(gradient*dirs(:, :, idir))
         scale = sum(abs(gradient*dirs(:, :, idir)))
         write (context, '(a,", ",a,", ",a,": direction ",i0," dE/dt ",es12.4," vs FD ",es12.4, &
            & ", dev ",es9.2)') trim(pcm%name), trim(mol_label(imol)), trim(seed_name(seed)), &
            & idir, analytic, fd, abs(analytic - fd)
         call check(error, abs(analytic - fd) <= nuc_atol + nuc_rtol*scale, &
            & message="nuclear gradient disagrees with the energy difference", &
            & more=trim(context))
         if (allocated(error)) return
      end do
   end subroutine check_nuclear_gradient

   !> Displacement directions of the nuclear check
   !>
   !> Every Cartesian coordinate up to `max_nat_sweep` atoms, else `ndir_nuc`
   !> fixed directions moving all nuclei by up to one unit per coordinate
   !>
   !> @param[in]  nat  Number of atoms
   !> @param[out] dirs Directions (3, nat, ndir)
   subroutine nuclear_directions(nat, dirs)
      !> Number of atoms
      integer, intent(in) :: nat
      !> Directions
      real(wp), allocatable, intent(out) :: dirs(:, :, :)
      !> Atom, axis and direction indices
      integer :: iatom, iaxis, idir

      if (nat <= max_nat_sweep) then
         allocate (dirs(3, nat, 3*nat), source=0.0_wp)
         do iatom = 1, nat
            do iaxis = 1, 3
               dirs(iaxis, iatom, 3*(iatom - 1) + iaxis) = 1.0_wp
            end do
         end do
      else
         allocate (dirs(3, nat, ndir_nuc))
         do idir = 1, ndir_nuc
            do iatom = 1, nat
               do iaxis = 1, 3
                  dirs(iaxis, iatom, idir) = sin(0.83_wp*real(3*(iatom - 1) + iaxis, wp) &
                     & + 1.37_wp*real(idir, wp))
               end do
            end do
         end do
      end if
   end subroutine nuclear_directions

   !* ================================================================================= *!
   !*                                 Reference energies                                *!
   !* ================================================================================= *!

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
   !> - Nonlocal molecules: no reference. NESSie.jl is tested only against the
   !>   sinh-form Born energy; check its Born limit first
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
      call add_reference_rows(rows, pcm_nonlocal, "none (NESSie.jl: check its Born limit)")
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

   !* ================================================================================= *!
   !*                               Component evaluation                                *!
   !* ================================================================================= *!

   !> Component factory
   !>
   !> Dense solves; SLIC with a Newton tolerance far below the differences
   !>
   !> @param[in]  kind    Model identifier
   !> @param[in]  ctx     Run context borrowed by the component
   !> @param[out] pcm     Constructed component
   !> @param[out] error   Test failure
   !> @param[in]  physics Physical parameters, those of the difference checks by default
   subroutine new_pcm(kind, ctx, pcm, error, physics)
      !> Model identifier
      integer, intent(in) :: kind
      !> Run context borrowed by the component
      type(moist_context_type), target, intent(in) :: ctx
      !> Constructed component
      class(model_continuum_component_pcm), allocatable, intent(out) :: pcm
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Physical parameters
      type(pcm_physics), intent(in), optional :: physics

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Solver controls
      type(moist_pcm_parameters_type) :: param
      !> Physical parameters in effect
      type(pcm_physics) :: phys

      param%solver = solver_type%dense
      param%newton_tol = newton_tol
      if (present(physics)) phys = physics
      select case (kind)
      case (pcm_cpcm)
         allocate (model_continuum_component_cpcm :: pcm)
      case (pcm_cosmo)
         allocate (model_continuum_component_cosmo :: pcm)
      case (pcm_iefpcm)
         allocate (model_continuum_component_iefpcm :: pcm)
      case (pcm_ssvpe)
         allocate (model_continuum_component_ssvpe :: pcm)
      case (pcm_dpcm)
         allocate (model_continuum_component_dpcm :: pcm)
      case (pcm_dpcm_original)
         allocate (model_continuum_component_dpcm_original :: pcm)
      case (pcm_lpb)
         allocate (model_continuum_component_lpb :: pcm)
      case (pcm_slic)
         allocate (model_continuum_component_slic :: pcm)
      case (pcm_nonlocal)
         allocate (model_continuum_component_nonlocal :: pcm)
      case default
         call test_failed(error, "unknown PCM model")
         return
      end select
      associate (eps => phys%eps)
         select type (pcm)
         type is (model_continuum_component_cpcm)
            call new_component_cpcm(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_cosmo)
            call new_component_cosmo(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_iefpcm)
            call new_component_iefpcm(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_ssvpe)
            call new_component_ssvpe(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_dpcm)
            call new_component_dpcm(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_dpcm_original)
            call new_component_dpcm_original(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_lpb)
            call new_component_lpb(pcm, eps, phys%kappa, err, param, ctx=ctx)
         type is (model_continuum_component_slic)
            call new_component_slic(pcm, eps, phys%slic_alpha, phys%slic_beta, phys%slic_gamma, &
               & err, param, mu=phys%slic_mu, ctx=ctx)
         type is (model_continuum_component_nonlocal)
            call new_component_nonlocal(pcm, eps, phys%eps_inf, phys%lambda, err, param, ctx=ctx)
         end select
      end associate
      call check_moist_error(error, err, "component setup")
   end subroutine new_pcm

   !> Surface channels whose weights and tangents must not vanish
   !>
   !> @param[in] kind    Model identifier
   !> @param[in] smeared Gaussian surface
   function active_channels(kind, smeared) result(channels)
      !> Model identifier
      integer, intent(in) :: kind
      !> Gaussian surface
      logical, intent(in) :: smeared
      !> Blank-separated channel names
      character(len=:), allocatable :: channels

      select case (kind)
      case (pcm_cpcm, pcm_cosmo)
         channels = merge("r xi f", "r a   ", smeared)
      case (pcm_iefpcm, pcm_ssvpe, pcm_lpb)
         channels = merge("r n k xi f a", "r n k a     ", smeared)
      case (pcm_dpcm, pcm_dpcm_original, pcm_slic)
         channels = merge("r n k xi a", "r n k a   ", smeared)
      case default
         channels = merge("r n xi a", "r n a   ", smeared)
      end select
      channels = trim(channels)
   end function active_channels

   !> Energy with the traces of host point charges
   !>
   !> @param[in,out] pcm    Component under test
   !> @param[in]     mol    Solute
   !> @param[in,out] cavity Cavity, current for `mol`
   !> @param[in]     src    Structure carrying the host charge positions
   !> @param[in]     qsrc   Host charges
   !> @param[out]    energy Polarization energy
   !> @param[out]    error  Test failure
   subroutine host_energy(pcm, mol, cavity, src, qsrc, energy, error)
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Cavity
      class(cavity_type), intent(inout) :: cavity
      !> Host charge positions
      type(structure_type), intent(in) :: src
      !> Host charges
      real(wp), intent(in) :: qsrc(:)
      !> Polarization energy
      real(wp), intent(out) :: energy
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Host coupling
      type(coupling_type), target :: coupling

      energy = 0.0_wp
      call stage_energy(pcm, mol, cavity, coupling, error)
      if (allocated(error)) return
      call answer_sources(cavity, coupling, src, qsrc, real(mol%num(mol%id), wp))
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err, "energy")
   end subroutine host_energy

   !> Total surface weights of the gradient phase with host point charges
   !>
   !> @param[in,out] pcm    Component under test
   !> @param[in]     mol    Solute
   !> @param[in,out] cavity Cavity
   !> @param[in]     q      Host charges on the solute atoms
   !> @param[out]    acc    Total surface weights
   !> @param[out]    energy Polarization energy
   !> @param[out]    error  Test failure
   subroutine host_surface_weights(pcm, mol, cavity, q, acc, energy, error)
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Cavity
      class(cavity_type), intent(inout) :: cavity
      !> Host charges on the solute atoms
      real(wp), intent(in) :: q(:)
      !> Total surface weights
      type(cavity_surface_adjoint_type), intent(out) :: acc
      !> Polarization energy
      real(wp), intent(out) :: energy
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Host coupling
      type(coupling_type), target :: coupling
      !> Nuclear charges, unused by the surface weights
      real(wp), allocatable :: z(:)

      energy = 0.0_wp
      z = real(mol%num(mol%id), wp)
      call stage_energy(pcm, mol, cavity, coupling, error)
      if (allocated(error)) return
      call answer_sources(cavity, coupling, mol, q, z)
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err, "energy")
      if (allocated(error)) return
      call pcm%prepare_gradient(cavity, coupling, err)
      call check_moist_error(error, err, "gradient staging")
      if (allocated(error)) return
      call answer_sources(cavity, coupling, mol, q, z)
      call acc%init(cavity%ngrid)
      call pcm%get_gradient_surface_weights(component_view(coupling), cavity, acc, err)
      call check_moist_error(error, err, "surface weights")
   end subroutine host_surface_weights

   !> Energy with the potential and its spatial derivative given directly
   !>
   !> @param[in,out] pcm     Component under test
   !> @param[in]     mol     Solute
   !> @param[in,out] cavity  Cavity
   !> @param[in]     phi     Potential trace (ngrid)
   !> @param[in]     dphi_dr Spatial derivative of the potential (3, ngrid)
   !> @param[out]    energy  Polarization energy
   !> @param[out]    error   Test failure
   subroutine trace_energy(pcm, mol, cavity, phi, dphi_dr, energy, error)
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Cavity
      class(cavity_type), intent(inout) :: cavity
      !> Potential trace
      real(wp), intent(in) :: phi(:)
      !> Spatial derivative of the potential
      real(wp), intent(in) :: dphi_dr(:, :)
      !> Polarization energy
      real(wp), intent(out) :: energy
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Host coupling
      type(coupling_type), target :: coupling

      energy = 0.0_wp
      call stage_energy(pcm, mol, cavity, coupling, error)
      if (allocated(error)) return
      call answer_traces(coupling, phi, dphi_dr)
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err, "energy")
   end subroutine trace_energy

   !> Trace adjoints of `get_trace_response` with the traces given directly
   !>
   !> @param[in,out] pcm       Component under test
   !> @param[in]     mol       Solute
   !> @param[in,out] cavity    Cavity
   !> @param[in]     phi       Potential trace (ngrid)
   !> @param[in]     dphi_dr   Spatial derivative of the potential (3, ngrid)
   !> @param[out]    w_phi     dE/dphi (ngrid)
   !> @param[out]    w_dphi_dr dE/d(dphi_dr) (3, ngrid), normal-trace models only
   !> @param[out]    error     Test failure
   subroutine trace_adjoints(pcm, mol, cavity, phi, dphi_dr, w_phi, w_dphi_dr, error)
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Cavity
      class(cavity_type), intent(inout) :: cavity
      !> Potential trace
      real(wp), intent(in) :: phi(:)
      !> Spatial derivative of the potential
      real(wp), intent(in) :: dphi_dr(:, :)
      !> Potential adjoint
      real(wp), allocatable, intent(out) :: w_phi(:)
      !> Adjoint of the spatial derivative
      real(wp), allocatable, intent(out) :: w_dphi_dr(:, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Host coupling
      type(coupling_type), target :: coupling
      !> Response accumulator and its potential adjoint item
      type(response_type) :: response
      type(potential_adjoint_response_type), allocatable :: item

      call stage_energy(pcm, mol, cavity, coupling, error)
      if (allocated(error)) return
      call answer_traces(coupling, phi, dphi_dr)
      call pcm%get_trace_response(component_view(coupling), cavity, response, err)
      call check_moist_error(error, err, "trace response")
      if (allocated(error)) return
      call copy_potential_adjoint(response, item)
      call check(error, allocated(item), message="trace response has no potential adjoint")
      if (allocated(error)) return
      call check(error, allocated(item%w_phi), message="potential adjoint has no w_phi")
      if (allocated(error)) return
      call move_alloc(item%w_phi, w_phi)
      if (allocated(item%w_dphi_dr)) call move_alloc(item%w_dphi_dr, w_dphi_dr)
   end subroutine trace_adjoints

   !> Operator weights and trace adjoints with the potential and normal trace fixed
   !>
   !> The spatial derivative is answered as vn n/|n|**2, so n . dphi_dr = vn
   !> holds on displaced (non-unit) normals too
   !>
   !> @param[in,out] pcm    Component under test
   !> @param[in]     mol    Solute
   !> @param[in,out] cavity Cavity
   !> @param[in]     phi    Potential trace (ngrid)
   !> @param[in]     vn     Normal trace (ngrid)
   !> @param[out]    acc    Operator surface weights
   !> @param[out]    w_phi  Potential adjoint (ngrid)
   !> @param[out]    w_vn   Normal-trace adjoint (ngrid), zero for potential-only models
   !> @param[out]    error  Test failure
   subroutine fixed_trace_state(pcm, mol, cavity, phi, vn, acc, w_phi, w_vn, error)
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Cavity
      class(cavity_type), intent(inout) :: cavity
      !> Potential trace
      real(wp), intent(in) :: phi(:)
      !> Normal trace
      real(wp), intent(in) :: vn(:)
      !> Operator surface weights
      type(cavity_surface_adjoint_type), intent(out) :: acc
      !> Potential adjoint
      real(wp), allocatable, intent(out) :: w_phi(:)
      !> Normal-trace adjoint
      real(wp), allocatable, intent(out) :: w_vn(:)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Host coupling
      type(coupling_type), target :: coupling
      !> Spatial derivative reproducing the normal trace
      real(wp), allocatable :: dphi_dr(:, :)
      !> Energy of the solve
      real(wp) :: energy
      !> Grid-point index
      integer :: i

      allocate (dphi_dr(3, cavity%ngrid))
      do i = 1, cavity%ngrid
         dphi_dr(:, i) = vn(i)*cavity%normal0(:, i)/sum(cavity%normal0(:, i)**2)
      end do
      call stage_energy(pcm, mol, cavity, coupling, error)
      if (allocated(error)) return
      call answer_traces(coupling, phi, dphi_dr)
      energy = 0.0_wp
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err, "energy")
      if (allocated(error)) return
      call acc%init(cavity%ngrid)
      call pcm%get_surface_weights(component_view(coupling), cavity, acc, err)
      call check_moist_error(error, err, "surface weights")
      if (allocated(error)) return
      call pcm%potential_adjoint(w_phi, err)
      call check_moist_error(error, err, "potential adjoint")
      if (allocated(error)) return
      call pcm%normal_trace_adjoint(w_vn, err)
      call check_moist_error(error, err, "normal-trace adjoint")
   end subroutine fixed_trace_state

   !> Rebuild the cavity and the component at a geometry and solve with the nuclear host
   !>
   !> @param[in,out] pcm      Component under test
   !> @param[in]     mol      Solute
   !> @param[in]     el_xyz   Fixed electron positions (3, nat)
   !> @param[in]     el_q     Fixed electron charges (nat)
   !> @param[in,out] cavity   Cavity to rebuild
   !> @param[out]    coupling Coupling staged for the energy phase and answered
   !> @param[out]    energy   Polarization energy
   !> @param[out]    error    Test failure
   subroutine nuclear_state(pcm, mol, el_xyz, el_q, cavity, coupling, energy, error)
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Fixed electron positions
      real(wp), intent(in) :: el_xyz(:, :)
      !> Fixed electron charges
      real(wp), intent(in) :: el_q(:)
      !> Cavity to rebuild
      class(cavity_type), intent(inout) :: cavity
      !> Host coupling
      type(coupling_type), intent(out), target :: coupling
      !> Polarization energy
      real(wp), intent(out) :: energy
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err

      energy = 0.0_wp
      call cavity%update(mol, err)
      call check_moist_error(error, err, "cavity rebuild")
      if (allocated(error)) return
      call stage_energy(pcm, mol, cavity, coupling, error)
      if (allocated(error)) return
      call answer_nuclear_host(cavity, coupling, mol, el_xyz, el_q)
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err, "energy")
   end subroutine nuclear_state

   !> Update the component and stage a new coupling for the energy phase
   !>
   !> @param[in,out] pcm      Component under test
   !> @param[in]     mol      Solute
   !> @param[in,out] cavity   Cavity, current for `mol`
   !> @param[out]    coupling Coupling staged for the energy phase
   !> @param[out]    error    Test failure
   subroutine stage_energy(pcm, mol, cavity, coupling, error)
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Cavity
      class(cavity_type), intent(inout) :: cavity
      !> Coupling staged for the energy phase
      type(coupling_type), intent(out), target :: coupling
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err

      call pcm%update(mol, cavity, err)
      call check_moist_error(error, err, "update")
      if (allocated(error)) return
      call pcm%new_coupling(cavity, coupling, err)
      call check_moist_error(error, err, "coupling setup")
      if (allocated(error)) return
      call pcm%prepare_energy(cavity, coupling, err)
      call check_moist_error(error, err, "energy staging")
   end subroutine stage_energy

   !* ================================================================================= *!
   !*                                   Host answers                                    *!
   !* ================================================================================= *!

   !> Answer every request of the staged phase from host point charges
   !>
   !> Potential and all its derivatives, Gaussian-probed on a Gaussian request;
   !> `z` as the nuclear charges
   !>
   !> @param[in]     cavity   Cavity the coupling was staged on
   !> @param[in,out] coupling Staged coupling
   !> @param[in]     src      Structure carrying the host charge positions
   !> @param[in]     qsrc     Host charges
   !> @param[in]     z        Nuclear charges of the solute
   subroutine answer_sources(cavity, coupling, src, qsrc, z)
      !> Cavity the coupling was staged on
      class(cavity_type), intent(in) :: cavity
      !> Staged coupling
      type(coupling_type), intent(inout), target :: coupling
      !> Host charge positions
      type(structure_type), intent(in) :: src
      !> Host charges
      real(wp), intent(in) :: qsrc(:)
      !> Nuclear charges of the solute
      real(wp), intent(in) :: z(:)

      !> Raw potential derivatives
      real(wp), allocatable :: g(:, :), gxi(:), hrr(:, :, :), hrxi(:, :)

      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            call point_charge_field(cavity, qsrc, src, .true., g, gxi, hrr, hrxi)
            call submit(coupling, "phi", point_charge_potential(cavity, qsrc, src, .true.))
            call submit(coupling, "dphi_dr", g)
            call submit(coupling, "dphi_dxi", gxi)
            call submit(coupling, "d2phi_dr2", hrr)
            call submit(coupling, "d2phi_dr_dxi", hrxi)
         type is (point_potential_request_type)
            call point_charge_field(cavity, qsrc, src, .false., g, gxi, hrr, hrxi)
            call submit(coupling, "phi", point_charge_potential(cavity, qsrc, src, .false.))
            call submit(coupling, "dphi_dr", g)
            call submit(coupling, "d2phi_dr2", hrr)
         type is (nuclear_charge_request_type)
            call submit(coupling, "z", z)
         end select
      end do
   end subroutine answer_sources

   !> Answer the potential and its spatial derivative of the energy phase directly
   !>
   !> @param[in,out] coupling Coupling staged for the energy phase
   !> @param[in]     phi      Potential trace (ngrid)
   !> @param[in]     dphi_dr  Spatial derivative of the potential (3, ngrid)
   subroutine answer_traces(coupling, phi, dphi_dr)
      !> Staged coupling
      type(coupling_type), intent(inout), target :: coupling
      !> Potential trace
      real(wp), intent(in) :: phi(:)
      !> Spatial derivative of the potential
      real(wp), intent(in) :: dphi_dr(:, :)

      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            call submit(coupling, "phi", phi)
            call submit(coupling, "dphi_dr", dphi_dr)
         type is (point_potential_request_type)
            call submit(coupling, "phi", phi)
            call submit(coupling, "dphi_dr", dphi_dr)
         end select
      end do
   end subroutine answer_traces

   !> Answer the host of moving nuclei and fixed electrons
   !>
   !> @param[in]     cavity   Cavity the coupling was staged on
   !> @param[in,out] coupling Staged coupling
   !> @param[in]     mol      Solute, its nuclei carry the charges Z
   !> @param[in]     el_xyz   Fixed electron positions (3, nat)
   !> @param[in]     el_q     Fixed electron charges (nat)
   subroutine answer_nuclear_host(cavity, coupling, mol, el_xyz, el_q)
      !> Cavity the coupling was staged on
      class(cavity_type), intent(in) :: cavity
      !> Staged coupling
      type(coupling_type), intent(inout), target :: coupling
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Fixed electron positions
      real(wp), intent(in) :: el_xyz(:, :)
      !> Fixed electron charges
      real(wp), intent(in) :: el_q(:)

      !> Host charge positions: nuclei, then electrons
      type(structure_type) :: src
      !> Nuclear charges
      real(wp), allocatable :: z(:)

      z = real(mol%num(mol%id), wp)
      call new(src, [mol%num(mol%id), spread(1, 1, size(el_q))], &
         & reshape([mol%xyz, el_xyz], [3, mol%nat + size(el_q)]))
      call answer_sources(cavity, coupling, src, [z, el_q], z)
   end subroutine answer_nuclear_host

   !> Host point charges on the solute atoms, summing to the solute charge
   !>
   !> @param[in] mol Solute
   pure function host_charges(mol) result(q)
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Host charges (nat)
      real(wp), allocatable :: q(:)
      !> Atom index
      integer :: iatom

      allocate (q(mol%nat))
      do iatom = 1, mol%nat
         q(iatom) = 0.6_wp*sin(1.7_wp*real(iatom, wp) + 0.3_wp)
      end do
      q = q + (mol%charge - sum(q))/real(mol%nat, wp)
   end function host_charges

   !> Electron charges of the nuclear check, fixed near the reference atoms
   !>
   !> Each atom keeps its host charge: the electron charge is -(Z - q), placed
   !> up to 0.25 bohr off the nucleus so that the host carries dipoles
   !>
   !> @param[in]  mol    Solute at the reference geometry
   !> @param[out] el_xyz Electron positions (3, nat)
   !> @param[out] el_q   Electron charges (nat)
   subroutine electron_host(mol, el_xyz, el_q)
      !> Solute at the reference geometry
      type(structure_type), intent(in) :: mol
      !> Electron positions
      real(wp), allocatable, intent(out) :: el_xyz(:, :)
      !> Electron charges
      real(wp), allocatable, intent(out) :: el_q(:)
      !> Atom index
      integer :: iatom

      el_q = host_charges(mol) - real(mol%num(mol%id), wp)
      allocate (el_xyz(3, mol%nat))
      do iatom = 1, mol%nat
         el_xyz(:, iatom) = mol%xyz(:, iatom) + 0.25_wp*sin(1.9_wp*real(iatom, wp) &
            & + [0.1_wp, 2.2_wp, 4.3_wp])
      end do
   end subroutine electron_host

   !> Potential and spatial derivative of host point charges on a cavity
   !>
   !> @param[in]  src     Structure carrying the host charge positions
   !> @param[in]  qsrc    Host charges
   !> @param[in]  cavity  Cavity
   !> @param[in]  smeared Gaussian-probed traces
   !> @param[out] phi     Potential trace (ngrid)
   !> @param[out] dphi_dr Spatial derivative of the potential (3, ngrid)
   subroutine host_traces(src, qsrc, cavity, smeared, phi, dphi_dr)
      !> Host charge positions
      type(structure_type), intent(in) :: src
      !> Host charges
      real(wp), intent(in) :: qsrc(:)
      !> Cavity
      class(cavity_type), intent(in) :: cavity
      !> Gaussian-probed traces
      logical, intent(in) :: smeared
      !> Potential trace
      real(wp), allocatable, intent(out) :: phi(:)
      !> Spatial derivative of the potential
      real(wp), allocatable, intent(out) :: dphi_dr(:, :)

      !> Unused higher derivatives
      real(wp), allocatable :: gxi(:), hrr(:, :, :), hrxi(:, :)

      phi = point_charge_potential(cavity, qsrc, src, smeared)
      call point_charge_field(cavity, qsrc, src, smeared, dphi_dr, gxi, hrr, hrxi)
   end subroutine host_traces

   !* ================================================================================= *!
   !*                                     Cavities                                      *!
   !* ================================================================================= *!

   !> Solute of the difference checks
   !>
   !> @param[in]  imol Solute index
   !> @param[out] mol  Solute
   subroutine get_solute(imol, mol)
      !> Solute index
      integer, intent(in) :: imol
      !> Solute
      type(structure_type), intent(out) :: mol

      select case (imol)
      case (mol_triatomic)
         call new(mol, solute_num, solute_xyz)
      case (mol_cross)
         call get_test_cross(mol)
      case default
         call get_structure(mol, trim(mol_collection(imol)), trim(mol_record(imol)))
      end select
   end subroutine get_solute

   !> Empty molecular cavity of the nuclear check, configured for the component
   !>
   !> @param[in]  seed   Cavity kind
   !> @param[in,out] pcm Component whose geometry needs are requested
   !> @param[in]  ctx    Run context borrowed by the cavity
   !> @param[out] cavity Cavity, built by `nuclear_state`
   !> @param[out] error  Test failure
   subroutine new_host_cavity(seed, pcm, ctx, cavity, error)
      !> Cavity kind
      integer, intent(in) :: seed
      !> Component whose geometry needs are requested
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(in) :: ctx
      !> Cavity
      class(cavity_type), allocatable, intent(out) :: cavity
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err

      call new_seed_cavity(seed, ctx, cavity, error)
      if (allocated(error)) return
      call pcm%configure_cavity(cavity, err)
      call check_moist_error(error, err, "cavity configuration")
   end subroutine new_host_cavity

   !> First cavity seed of a model's checks: SvdW-DROP for SLIC, else iSwiG
   !>
   !> Concave iSwiG seams can fold the SLIC interface condition (the component
   !> warns there), so SLIC is checked on the smooth cavity only
   !>
   !> @param[in] kind Model identifier
   pure function first_seed(kind) result(seed)
      !> Model identifier
      integer, intent(in) :: kind
      !> Cavity seed
      integer :: seed

      seed = merge(seed_drop, seed_iswig, kind == pcm_slic)
   end function first_seed

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

   !> Empty iSwiG or SvdW-DROP cavity with CPCM radii
   !>
   !> DROP projects to 1e-14 so that its surface is smooth at the difference steps
   !>
   !> @param[in]  seed   Cavity kind
   !> @param[in]  ctx    Run context borrowed by the cavity
   !> @param[out] cavity Cavity
   !> @param[out] error  Test failure
   subroutine new_seed_cavity(seed, ctx, cavity, error)
      !> Cavity kind
      integer, intent(in) :: seed
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(in) :: ctx
      !> Cavity
      class(cavity_type), allocatable, intent(out) :: cavity
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Concrete cavities
      type(cavity_type_iswig), allocatable :: iswig
      type(cavity_type_drop), allocatable :: drop
      !> Level set of the DROP cavity
      type(moist_cavity_drop_lsf_svdw_type) :: svdw

      select case (seed)
      case (seed_iswig)
         allocate (iswig)
         call new_cavity_iswig(iswig, radius_model=default_cpcm_radii(), error=err, &
            & param=moist_cavity_iswig_parameters_type(num_leb=nleb), ctx=ctx)
         call check_moist_error(error, err, "iSwiG setup")
         if (allocated(error)) return
         call move_alloc(iswig, cavity)
      case (seed_drop)
         allocate (drop)
         call svdw%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=1.5_wp, blend_3b=1.0_wp))
         call new_cavity_drop(drop, radius_model=default_cpcm_radii(), lsf_model=svdw, error=err, &
            & param=moist_cavity_drop_parameters_type(num_leb=nleb, tolerance=1.0e-14_wp, &
            & proj_maxiter=150, proj_level=2), ctx=ctx)
         call check_moist_error(error, err, "DROP setup")
         if (allocated(error)) return
         call move_alloc(drop, cavity)
      case default
         call test_failed(error, "unknown cavity seed")
      end select
   end subroutine new_seed_cavity

   !> Synthetic cavity seeded from the surface of the solute
   !>
   !> - iSwiG: k1 = k2, modulated around 1/R so that the curvature is not
   !>   piecewise constant
   !> - DROP: its principal curvatures, k1 /= k2
   !> - Only switched-on points (f > cut_f); widths and switching factors only
   !>   for a Gaussian surface
   !> - Point surface: no two points closer than `min_sep`, the most
   !>   switched-on point of a close group kept
   !>
   !> @param[in]  seed    Seed cavity
   !> @param[in]  smeared Keep the widths and switching factors
   !> @param[in]  mol     Solute
   !> @param[out] cavity  Synthetic cavity
   !> @param[out] error   Test failure
   subroutine get_seed(seed, smeared, mol, cavity, error)
      !> Seed cavity
      integer, intent(in) :: seed
      !> Keep the widths and switching factors
      logical, intent(in) :: smeared
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Synthetic cavity
      type(fd_cavity), intent(out) :: cavity
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the seeding cavity
      type(moist_context_type), target :: ctx
      !> Seeding cavity
      class(cavity_type), allocatable :: source
      !> Seed points kept
      integer, allocatable :: keep(:)
      !> Grid-point index
      integer :: i

      call new_context(ctx, verbosity=0)
      call new_seed_cavity(seed, ctx, source, error)
      if (allocated(error)) return
      call source%configure(normal=.true., curvature=.true., error=err)
      call check_moist_error(error, err, "seed configuration")
      if (allocated(error)) return
      call source%update(mol, err)
      call check_moist_error(error, err, "seed update")
      if (allocated(error)) return

      keep = pack([(i, i=1, source%ngrid)], source%f(1:source%ngrid) > cut_f)
      if (.not. smeared) call separate_points(source%xyz, source%f, keep)
      cavity%ngrid = size(keep)
      cavity%nsph = source%nsph
      cavity%xyz = source%xyz(:, keep)
      cavity%normal0 = source%normal0(:, keep)
      cavity%a = source%a(keep)
      cavity%k1 = source%k1(keep)
      cavity%k2 = source%k2(keep)
      if (seed == seed_iswig) then
         do i = 1, cavity%ngrid
            cavity%k1(i) = cavity%k1(i)*(1.0_wp + 0.2_wp*sin(0.37_wp*real(i, wp)))
         end do
         cavity%k2 = cavity%k1
      end if
      if (smeared) then
         cavity%xi0 = source%xi0(keep)
         cavity%f = source%f(keep)
      end if
      cavity%updated = .true.
      cavity%built = .true.
      call check(error, cavity%ngrid > 10, message="unexpected seed size")
      if (allocated(error)) return
      if (seed == seed_drop .and. mol%nat > 1) then
         call check(error, maxval(abs(cavity%k1 - cavity%k2)) > thr_vacuous, &
            & message="DROP seed: principal curvatures coincide")
      end if
   end subroutine get_seed

   !> Thin a point selection so that no two points lie closer than `min_sep`
   !>
   !> Greedy in order of decreasing switching factor; the kept indices stay in
   !> ascending order
   !>
   !> @param[in]     xyz  Point positions (3, n)
   !> @param[in]     f    Switching factors (n)
   !> @param[in,out] keep Selected point indices
   subroutine separate_points(xyz, f, keep)
      !> Point positions
      real(wp), intent(in) :: xyz(:, :)
      !> Switching factors
      real(wp), intent(in) :: f(:)
      !> Selected point indices
      integer, allocatable, intent(inout) :: keep(:)

      !> Candidates by decreasing f and the accepted mask
      integer, allocatable :: order(:)
      logical, allocatable :: accepted(:)
      !> Candidate indices
      integer :: i, j, ibest

      order = keep
      do i = 1, size(order) - 1
         ibest = i - 1 + maxloc(f(order(i:)), dim=1)
         j = order(i)
         order(i) = order(ibest)
         order(ibest) = j
      end do
      allocate (accepted(size(xyz, 2)), source=.false.)
      do i = 1, size(order)
         accepted(order(i)) = all(norm2(xyz(:, order(:i - 1)) &
            & - spread(xyz(:, order(i)), 2, i - 1), dim=1) >= min_sep .or. &
            & .not. accepted(order(:i - 1)))
      end do
      keep = pack(keep, accepted(keep))
   end subroutine separate_points

   !> No-op update: the fields are set by the test
   !>
   !> @param[in,out] self  Synthetic cavity
   !> @param[in]     mol   Solute, unused
   !> @param[out]    error Never allocated
   subroutine fd_cavity_update(self, mol, error)
      !> Synthetic cavity
      class(fd_cavity), intent(inout) :: self
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      associate (unused_self => self, unused_mol => mol)
      end associate
   end subroutine fd_cavity_update

   !> No-op: the synthetic surface has no nuclear derivatives
   !>
   !> @param[in,out] self  Synthetic cavity
   !> @param[out]    error Never allocated
   subroutine fd_cavity_get_gradient(self, error)
      !> Synthetic cavity
      class(fd_cavity), intent(inout) :: self
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      associate (unused_self => self)
      end associate
   end subroutine fd_cavity_get_gradient

   !* ================================================================================= *!
   !*                                    Directions                                     *!
   !* ================================================================================= *!

   !> Whether the cavity provides a surface channel
   !>
   !> @param[in] cavity Cavity
   !> @param[in] ich    Channel index
   pure function provides(cavity, ich) result(has)
      !> Cavity
      type(fd_cavity), intent(in) :: cavity
      !> Channel index
      integer, intent(in) :: ich
      !> The channel can be perturbed
      logical :: has

      select case (ich)
      case (ch_xi)
         has = allocated(cavity%xi0)
      case (ch_f)
         has = allocated(cavity%f)
      case default
         has = .true.
      end select
   end function provides

   !> Whether a channel is in a blank-separated channel list
   !>
   !> @param[in] channels Channel list
   !> @param[in] ich      Channel index
   pure function expects(channels, ich) result(listed)
      !> Channel list
      character(len=*), intent(in) :: channels
      !> Channel index
      integer, intent(in) :: ich
      !> The channel is listed
      logical :: listed

      listed = index(" "//channels//" ", " "//trim(channel_name(ich))//" ") > 0
   end function expects

   !> Deterministic direction of a trace
   !>
   !> @param[in]  n    Number of grid points
   !> @param[in]  seed Variant of the direction
   !> @param[out] d    Direction (n)
   subroutine trace_direction(n, seed, d)
      !> Number of grid points
      integer, intent(in) :: n
      !> Variant of the direction
      integer, intent(in) :: seed
      !> Direction
      real(wp), allocatable, intent(out) :: d(:)
      !> Grid-point index
      integer :: i

      allocate (d(n))
      do i = 1, n
         d(i) = 0.1_wp*sin(1.31_wp*real(i, wp) + 0.77_wp*real(seed, wp))
      end do
   end subroutine trace_direction

   !> Deterministic direction of the spatial derivative of the potential
   !>
   !> @param[in]  n Number of grid points
   !> @param[out] d Direction (3, n)
   subroutine field_direction(n, d)
      !> Number of grid points
      integer, intent(in) :: n
      !> Direction
      real(wp), allocatable, intent(out) :: d(:, :)
      !> Entry index
      integer :: i

      allocate (d(3, n))
      d = reshape([(0.01_wp*cos(0.53_wp*real(i, wp) + 0.2_wp), i=1, 3*n)], [3, n])
   end subroutine field_direction

   !> Deallocate every channel of a direction (a zero direction)
   !>
   !> @param[in,out] direction Direction
   subroutine clear_direction(direction)
      !> Direction
      type(cavity_surface_adjoint_type), intent(inout) :: direction

      if (allocated(direction%w_xyz)) deallocate (direction%w_xyz)
      if (allocated(direction%w_n)) deallocate (direction%w_n)
      if (allocated(direction%w_k1)) deallocate (direction%w_k1)
      if (allocated(direction%w_k2)) deallocate (direction%w_k2)
      if (allocated(direction%w_xi)) deallocate (direction%w_xi)
      if (allocated(direction%w_f)) deallocate (direction%w_f)
      if (allocated(direction%w_a)) deallocate (direction%w_a)
      if (allocated(direction%w_w)) deallocate (direction%w_w)
      if (allocated(direction%w_c)) deallocate (direction%w_c)
   end subroutine clear_direction

   !> Direction along a single surface channel, every other channel unallocated
   !>
   !> @param[in]  cavity    Cavity setting the scale of the fields
   !> @param[in]  ich       Channel index
   !> @param[out] direction Direction
   subroutine channel_direction(cavity, ich, direction)
      !> Cavity setting the scale of the fields
      type(fd_cavity), intent(in) :: cavity
      !> Channel index
      integer, intent(in) :: ich
      !> Direction
      type(cavity_surface_adjoint_type), intent(out) :: direction

      call add_channel_direction(cavity, ich, direction)
   end subroutine channel_direction

   !> Set one channel of a direction, deterministic and scaled to the field
   !>
   !> Positions and normals move by O(1), the scalar fields relative to their
   !> own size; the principal curvatures move along independent patterns where
   !> they differ (DROP seed), else along the same one
   !>
   !> @param[in]     cavity    Cavity setting the scale of the fields
   !> @param[in]     ich       Channel index
   !> @param[in,out] direction Direction
   subroutine add_channel_direction(cavity, ich, direction)
      !> Cavity setting the scale of the fields
      type(fd_cavity), intent(in) :: cavity
      !> Channel index
      integer, intent(in) :: ich
      !> Direction
      type(cavity_surface_adjoint_type), intent(inout) :: direction

      !> Deterministic unit-scale pattern
      real(wp), allocatable :: s(:)
      !> Grid size and entry index
      integer :: ngrid, i

      ngrid = cavity%ngrid
      allocate (s(3*ngrid))
      do i = 1, 3*ngrid
         s(i) = sin(0.91_wp*real(i, wp) + 1.7_wp*real(ich, wp))
      end do
      select case (ich)
      case (ch_r)
         direction%w_xyz = reshape(s, [3, ngrid])
      case (ch_n)
         direction%w_n = reshape(s, [3, ngrid])
      case (ch_k)
         direction%w_k1 = 0.5_wp*cavity%k1*s(1:ngrid)
         if (any(cavity%k1 /= cavity%k2)) then
            direction%w_k2 = 0.5_wp*cavity%k2*s(ngrid + 1:2*ngrid)
         else
            direction%w_k2 = direction%w_k1
         end if
      case (ch_xi)
         direction%w_xi = 0.5_wp*cavity%xi0*s(1:ngrid)
      case (ch_f)
         direction%w_f = 0.5_wp*cavity%f*s(1:ngrid)
      case (ch_a)
         direction%w_a = 0.5_wp*cavity%a*s(1:ngrid)
      case default
         error stop "add_channel_direction: unknown channel"
      end select
   end subroutine add_channel_direction

   !> Directional derivative of the energy from surface weights
   !>
   !> @param[in] acc       Surface weights
   !> @param[in] direction Direction, unallocated channels are zero
   pure function contract(acc, direction) result(dedt)
      !> Surface weights
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Direction
      type(cavity_surface_adjoint_type), intent(in) :: direction
      !> Directional derivative
      real(wp) :: dedt

      dedt = 0.0_wp
      if (allocated(direction%w_xyz)) dedt = dedt + sum(acc%w_xyz*direction%w_xyz)
      if (allocated(direction%w_n)) dedt = dedt + sum(acc%w_n*direction%w_n)
      if (allocated(direction%w_k1)) dedt = dedt + sum(acc%w_k1*direction%w_k1)
      if (allocated(direction%w_k2)) dedt = dedt + sum(acc%w_k2*direction%w_k2)
      if (allocated(direction%w_xi)) dedt = dedt + sum(acc%w_xi*direction%w_xi)
      if (allocated(direction%w_f)) dedt = dedt + sum(acc%w_f*direction%w_f)
      if (allocated(direction%w_a)) dedt = dedt + sum(acc%w_a*direction%w_a)
   end function contract

   !> Cavity moved by t along a direction
   !>
   !> @param[in]  cavity    Reference cavity
   !> @param[in]  direction Direction
   !> @param[in]  t         Step
   !> @param[out] trial     Displaced cavity
   subroutine displace(cavity, direction, t, trial)
      !> Reference cavity
      type(fd_cavity), intent(in) :: cavity
      !> Direction
      type(cavity_surface_adjoint_type), intent(in) :: direction
      !> Step
      real(wp), intent(in) :: t
      !> Displaced cavity
      type(fd_cavity), intent(out) :: trial

      trial = cavity
      if (allocated(direction%w_xyz)) trial%xyz = trial%xyz + t*direction%w_xyz
      if (allocated(direction%w_n)) trial%normal0 = trial%normal0 + t*direction%w_n
      if (allocated(direction%w_k1)) trial%k1 = trial%k1 + t*direction%w_k1
      if (allocated(direction%w_k2)) trial%k2 = trial%k2 + t*direction%w_k2
      if (allocated(direction%w_xi)) trial%xi0 = trial%xi0 + t*direction%w_xi
      if (allocated(direction%w_f)) trial%f = trial%f + t*direction%w_f
      if (allocated(direction%w_a)) trial%a = trial%a + t*direction%w_a
   end subroutine displace

end module test_model_component_pcm_general
