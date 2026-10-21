!> Derivatives of the PCM components against finite differences
!>
!> Generic checks through `model_continuum_component_pcm`, every model on every
!> solute of `mol_label`:
!>
!> - `surface_gradient`: `get_gradient_surface_weights` with a point-charge host
!>   against six-point differences of the energy along every surface channel,
!>   and the trace adjoints of `get_trace_response` against differences over
!>   the potential and its spatial derivative; Gaussian and point surfaces
!> - `tangent`: `surface_weight_tangent` against six-point differences of the
!>   operator weights and trace adjoints at fixed traces; iSwiG seed (k1 = k2)
!>   and SvdW-DROP seed (k1 /= k2), Gaussian and point surfaces
!> - `mf_tangent`: the matrix-free GMRES tangent against the dense analytic
!>   tangent of the same state, the dense one being the one `tangent` validates
!> - `nuclear_gradient`: `get_gradient` with moving nuclei and fixed electron
!>   charges against energy differences, the cavity rebuilt at every geometry
!>
!> The surface checks set the fields of a synthetic cavity seeded from an iSwiG
!> or SvdW-DROP surface; SLIC runs on SvdW-DROP only (its interface condition
!> can fold on iSwiG seams). Beyond the generic checks:
!>
!> - `component_paths`: every model on the shared seven-point fixture (k1 /= k2)
!>   with a quadratic host, `get_response` adjoints and each geometry channel
!>   separately differenced, dense against iterative solves
!> - CPCM host features: an independently coded energy on the fixture, an ECP
!>   host with forward/reverse gradient parity, a supplied external potential
!> - `cpcm_iterative_surface_gradient`: `surface_gradient` with GMRES at the
!>   default tolerance, whose residual-corrected energy the weights differentiate
!> - `point_potential_*`: `surface_gradient`, `nuclear_gradient` and
!>   `component_paths` of every model on Gaussian surfaces coupled through
!>   point potential integrals; every check also confirms the staged request
module test_model_component_pcm_derivatives
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: fd4_scalar, fd4_offsets, fd6_scalar, fd6_offsets, check_moist_error, &
      & component_view, copy_potential_adjoint, point_charge_potential, point_charge_field, submit, &
      & get_test_structures, center_at_origin, get_test_cavity_iswig, &
      & fill_point_charge_potential, fill_point_charge_field, set_host_potential, &
      & fill_missing_with_zeros, stage_point_charge_energy
   use test_model_component_helper, only: surface_fixture, new_surface_fixture, &
      & check_surface_weights, nfixture => fixture_ngrid_param, fixture_xyz => fixture_xyz_param, &
      & fixture_normals => fixture_normals_param, fixture_areas => fixture_areas_param, &
      & fixture_xis => fixture_xis_param, fixture_fs => fixture_fs_param
   use test_pcm_conductor_probe, only: probe_conductor_bilinear_weights, &
      & probe_conductor_nuclear_gradient
   use test_model_component_pcm_fixtures, only: pcm_cpcm, pcm_cosmo, pcm_iefpcm, pcm_ssvpe, &
      & pcm_dpcm, pcm_dpcm_original, pcm_lpb, pcm_slic, pcm_nonlocal, npcm, new_pcm, make_component, &
      & seed_drop, seed_name, first_seed, new_seed_cavity, fd_cavity, get_seed, nmol, &
      & mol_label, get_solute, thr_vacuous, pcm_physics, stage_energy, answer_sources, host_energy, host_charges
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_iswig, only: cavity_type_iswig
   use moist_cavity_drop, only: cavity_type_drop
   use moist_cavity_drop_lsf_isodensity_internal, only: moist_cavity_drop_lsf_isodensity_internal_type
   use moist_channels_coupling, only: coupling_type, gaussian_potential_request_type, &
      & point_potential_request_type, nuclear_charge_request_type
   use moist_channels_response, only: response_type, potential_adjoint_response_type, &
      & parameter_adjoint_response_type
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, &
      & moist_pcm_parameters_type, solver_type
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, &
      & new_component_cpcm
   use moist_model_continuum_component_pcm_solvers, only: solve_pcm_lu, pcm_storage_matrix_free
   use moist_model_continuum_component_pcm_slic, only: model_continuum_component_slic, &
      & new_component_slic
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_derivatives

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
   !> Energy differences with GMRES at the default tolerance, absolute part
   !>
   !> The iterate's error enters to first order in the residual: measured 1.4e-10
   !> with the residual-corrected energy, up to 1.8e-7 with 0.5*dot(q, phi)
   real(wp), parameter :: gmres_atol = 1.5e-9_wp
   !> Nuclear gradient (four-point stencil), with cavity rebuild roundoff
   !>
   !> The SLIC Newton tolerance sets the largest deviation, 9e-10 relative
   real(wp), parameter :: nuc_atol = 1.0e-12_wp, nuc_rtol = 1.0e-8_wp
   !> Parameter derivatives against six-point energy differences
   !>
   !> Absolute part relative to max(1, |E|), relative part of the difference
   !>
   !> Measured over every model and parameter, Gaussian and point: at most
   !> 3.5e-14 absolute and 1.4e-11 relative (SS(V)PE eps), the six-point
   !> stencil being roundoff-limited; tolerances are about ten times that
   real(wp), parameter :: param_atol = 1.0e-15_wp, param_rtol = 1.5e-10_wp
   !> Weight differences, per channel relative to its largest entry
   real(wp), parameter :: hess_atol = 3.0e-12_wp, hess_rtol = 1.0e-9_wp
   !> Matrix-free tangent against the dense analytic tangent, per channel
   !>
   !> Absolute part, and relative part of the largest dense entry; the dense
   !> solves are direct and the GMRES tolerance is 1e-13
   !>
   !> Measured over every model, solute, seed, surface and direction: at most
   !> 4.8e-14 absolute (CPCM dE/dphi) and 1.5e-12 relative to the channel scale
   !> (entries above 1e-3); tolerances are about ten times that
   real(wp), parameter :: mf_atol = 5.0e-13_wp, mf_rtol = 1.5e-11_wp
   !> Energy derivative relative to |E| below which an expected-active direction is vacuous
   !>
   !> The difference noise is about 1e-12 |E|
   real(wp), parameter :: thr_active = 1.0e-8_wp

   !> Largest solute whose nuclear coordinates and tangent directions are swept one by one
   !>
   !> Larger solutes get `ndir_nuc` directions over all nuclei and the joint
   !> tangent direction only
   integer, parameter :: max_nat_sweep = 3
   !> Directions of the nuclear check on the larger solutes
   integer, parameter :: ndir_nuc = 3

   !> Tolerance of identities that hold to roundoff
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   !> Absolute tolerance for charge and energy identities
   real(wp), parameter :: thr2 = 1.0e-12_wp

   !> Host surface potential driving the surface-weight test
   real(wp), parameter :: sw_phi(nfixture) = [-0.31_wp, 0.22_wp, -0.17_wp, 0.41_wp, &
                                              -0.28_wp, 0.13_wp, -0.36_wp]

   !> iSwiG fixture with manufactured nuclear width derivatives
   !>
   !> Ordinary iSwiG widths are fixed under nuclear motion. The width-channel
   !> regression supplies its own xi1_rA, contracted through the cavity callback.
   type, extends(cavity_type_iswig) :: width_gradient_cavity
   contains
      procedure :: contract_surface_gradient => fixture_surface_gradient
   end type width_gradient_cavity

contains

   !> Collect the derivative tests: per model surface gradient, tangents, nuclear gradient and
   !> parameter derivatives; then the fixture and CPCM checks
   !>
   !> @param[out] testsuite Collection of tests
   subroutine collect_model_component_pcm_derivatives(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("cpcm_surface_gradient", test_cpcm_surface_gradient), &
         & new_unittest("cpcm_iterative_surface_gradient", test_cpcm_iterative_surface_gradient), &
         & new_unittest("cpcm_tangent", test_cpcm_tangent), &
         & new_unittest("cpcm_mf_tangent", test_cpcm_mf_tangent), &
         & new_unittest("cpcm_nuclear_gradient", test_cpcm_nuclear_gradient), &
         & new_unittest("cpcm_parameter_derivatives_gaussian", test_cpcm_parameter_derivatives_gaussian), &
         & new_unittest("cpcm_parameter_derivatives_point", test_cpcm_parameter_derivatives_point), &
         & new_unittest("cpcm_mf_parameter_derivatives_gaussian", test_cpcm_mf_parameter_derivatives_gaussian), &
         & new_unittest("cpcm_mf_parameter_derivatives_point", test_cpcm_mf_parameter_derivatives_point), &
         & new_unittest("cosmo_surface_gradient", test_cosmo_surface_gradient), &
         & new_unittest("cosmo_tangent", test_cosmo_tangent), &
         & new_unittest("cosmo_mf_tangent", test_cosmo_mf_tangent), &
         & new_unittest("cosmo_nuclear_gradient", test_cosmo_nuclear_gradient), &
         & new_unittest("cosmo_parameter_derivatives_gaussian", test_cosmo_parameter_derivatives_gaussian), &
         & new_unittest("cosmo_parameter_derivatives_point", test_cosmo_parameter_derivatives_point), &
         & new_unittest("cosmo_mf_parameter_derivatives_gaussian", test_cosmo_mf_parameter_derivatives_gaussian), &
         & new_unittest("cosmo_mf_parameter_derivatives_point", test_cosmo_mf_parameter_derivatives_point), &
         & new_unittest("iefpcm_surface_gradient", test_iefpcm_surface_gradient), &
         & new_unittest("iefpcm_tangent", test_iefpcm_tangent), &
         & new_unittest("iefpcm_mf_tangent", test_iefpcm_mf_tangent), &
         & new_unittest("iefpcm_nuclear_gradient", test_iefpcm_nuclear_gradient), &
         & new_unittest("iefpcm_parameter_derivatives_gaussian", test_iefpcm_parameter_derivatives_gaussian), &
         & new_unittest("iefpcm_parameter_derivatives_point", test_iefpcm_parameter_derivatives_point), &
         & new_unittest("iefpcm_mf_parameter_derivatives_gaussian", test_iefpcm_mf_parameter_derivatives_gaussian), &
         & new_unittest("iefpcm_mf_parameter_derivatives_point", test_iefpcm_mf_parameter_derivatives_point), &
         & new_unittest("ssvpe_surface_gradient", test_ssvpe_surface_gradient), &
         & new_unittest("ssvpe_tangent", test_ssvpe_tangent), &
         & new_unittest("ssvpe_mf_tangent", test_ssvpe_mf_tangent), &
         & new_unittest("ssvpe_nuclear_gradient", test_ssvpe_nuclear_gradient), &
         & new_unittest("ssvpe_parameter_derivatives_gaussian", test_ssvpe_parameter_derivatives_gaussian), &
         & new_unittest("ssvpe_parameter_derivatives_point", test_ssvpe_parameter_derivatives_point), &
         & new_unittest("ssvpe_mf_parameter_derivatives_gaussian", test_ssvpe_mf_parameter_derivatives_gaussian), &
         & new_unittest("ssvpe_mf_parameter_derivatives_point", test_ssvpe_mf_parameter_derivatives_point), &
         & new_unittest("dpcm_surface_gradient", test_dpcm_surface_gradient), &
         & new_unittest("dpcm_tangent", test_dpcm_tangent), &
         & new_unittest("dpcm_mf_tangent", test_dpcm_mf_tangent), &
         & new_unittest("dpcm_nuclear_gradient", test_dpcm_nuclear_gradient), &
         & new_unittest("dpcm_parameter_derivatives_gaussian", test_dpcm_parameter_derivatives_gaussian), &
         & new_unittest("dpcm_parameter_derivatives_point", test_dpcm_parameter_derivatives_point), &
         & new_unittest("dpcm_mf_parameter_derivatives_gaussian", test_dpcm_mf_parameter_derivatives_gaussian), &
         & new_unittest("dpcm_mf_parameter_derivatives_point", test_dpcm_mf_parameter_derivatives_point), &
         & new_unittest("dpcm_original_surface_gradient", test_dpcm_original_surface_gradient), &
         & new_unittest("dpcm_original_tangent", test_dpcm_original_tangent), &
         & new_unittest("dpcm_original_mf_tangent", test_dpcm_original_mf_tangent), &
         & new_unittest("dpcm_original_nuclear_gradient", test_dpcm_original_nuclear_gradient), &
         & new_unittest("dpcm_original_parameter_derivatives_gaussian", test_dpcm_original_parameter_derivatives_gaussian), &
         & new_unittest("dpcm_original_parameter_derivatives_point", test_dpcm_original_parameter_derivatives_point), &
         & new_unittest("dpcm_original_mf_parameter_derivatives_gaussian", test_dpcm_original_mf_parameter_derivatives_gaussian), &
         & new_unittest("dpcm_original_mf_parameter_derivatives_point", test_dpcm_original_mf_parameter_derivatives_point), &
         & new_unittest("lpb_surface_gradient", test_lpb_surface_gradient), &
         & new_unittest("lpb_tangent", test_lpb_tangent), &
         & new_unittest("lpb_mf_tangent", test_lpb_mf_tangent), &
         & new_unittest("lpb_nuclear_gradient", test_lpb_nuclear_gradient), &
         & new_unittest("lpb_parameter_derivatives_gaussian", test_lpb_parameter_derivatives_gaussian), &
         & new_unittest("lpb_parameter_derivatives_point", test_lpb_parameter_derivatives_point), &
         & new_unittest("lpb_mf_parameter_derivatives_gaussian", test_lpb_mf_parameter_derivatives_gaussian), &
         & new_unittest("lpb_mf_parameter_derivatives_point", test_lpb_mf_parameter_derivatives_point), &
         & new_unittest("slic_surface_gradient", test_slic_surface_gradient), &
         & new_unittest("slic_tangent", test_slic_tangent), &
         & new_unittest("slic_mf_tangent", test_slic_mf_tangent), &
         & new_unittest("slic_nuclear_gradient", test_slic_nuclear_gradient), &
         & new_unittest("slic_parameter_derivatives_gaussian", test_slic_parameter_derivatives_gaussian), &
         & new_unittest("slic_parameter_derivatives_point", test_slic_parameter_derivatives_point), &
         & new_unittest("slic_mf_parameter_derivatives_gaussian", test_slic_mf_parameter_derivatives_gaussian), &
         & new_unittest("slic_mf_parameter_derivatives_point", test_slic_mf_parameter_derivatives_point), &
         & new_unittest("slic_default_mu_parameter_derivatives_gaussian", &
         & test_slic_default_mu_parameter_derivatives_gaussian), &
         & new_unittest("slic_default_mu_parameter_derivatives_point", &
         & test_slic_default_mu_parameter_derivatives_point), &
         & new_unittest("slic_default_mu_mf_parameter_derivatives_gaussian", &
         & test_slic_default_mu_mf_parameter_derivatives_gaussian), &
         & new_unittest("slic_default_mu_mf_parameter_derivatives_point", &
         & test_slic_default_mu_mf_parameter_derivatives_point), &
         & new_unittest("nonlocal_surface_gradient", test_nonlocal_surface_gradient), &
         & new_unittest("nonlocal_tangent", test_nonlocal_tangent), &
         & new_unittest("nonlocal_mf_tangent", test_nonlocal_mf_tangent), &
         & new_unittest("nonlocal_nuclear_gradient", test_nonlocal_nuclear_gradient), &
         & new_unittest("nonlocal_parameter_derivatives_gaussian", test_nonlocal_parameter_derivatives_gaussian), &
         & new_unittest("nonlocal_parameter_derivatives_point", test_nonlocal_parameter_derivatives_point), &
         & new_unittest("nonlocal_mf_parameter_derivatives_gaussian", test_nonlocal_mf_parameter_derivatives_gaussian), &
         & new_unittest("nonlocal_mf_parameter_derivatives_point", test_nonlocal_mf_parameter_derivatives_point), &
         & new_unittest("component_paths_gaussian", test_gaussian_component_paths), &
         & new_unittest("component_paths_point", test_point_component_paths), &
         & new_unittest("point_potential_surface_gradient", test_point_potential_surface_gradient), &
         & new_unittest("point_potential_nuclear_gradient", test_point_potential_nuclear_gradient), &
         & new_unittest("point_potential_component_paths", test_point_potential_component_paths), &
         & new_unittest("cpcm_surface_weights", test_cpcm_surface_weights), &
         & new_unittest("cpcm_nuclear_gradient_ecp", test_cpcm_nuclear_gradient_ecp), &
         & new_unittest("cpcm_external_nuclear_gradient", test_cpcm_external_nuclear_gradient) &
         & ]
   end subroutine collect_model_component_pcm_derivatives

   !* ================================================================================= *!
   !*                                 Per-model entries                                 *!
   !* ================================================================================= *!

   !> CPCM: surface and trace gradients
   subroutine test_cpcm_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_cpcm)
   end subroutine test_cpcm_surface_gradient

   !> CPCM: surface and trace gradients with inexact GMRES solves
   !>
   !> Default solver tolerance; the energy carries the residual correction of
   !> the generated solve, so the analytic weights stay its exact derivative
   subroutine test_cpcm_iterative_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_cpcm, iterative=.true.)
   end subroutine test_cpcm_iterative_surface_gradient

   !> CPCM: surface weight tangent
   subroutine test_cpcm_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_cpcm)
   end subroutine test_cpcm_tangent

   !> CPCM: surface weight tangent from matrix-free GMRES, against dense
   subroutine test_cpcm_mf_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_cpcm, mf=.true.)
   end subroutine test_cpcm_mf_tangent

   !> CPCM: nuclear gradient
   subroutine test_cpcm_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_cpcm)
   end subroutine test_cpcm_nuclear_gradient

   !> COSMO: surface and trace gradients
   subroutine test_cosmo_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_cosmo)
   end subroutine test_cosmo_surface_gradient

   !> COSMO: surface weight tangent
   subroutine test_cosmo_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_cosmo)
   end subroutine test_cosmo_tangent

   !> COSMO: surface weight tangent from matrix-free GMRES, against dense
   subroutine test_cosmo_mf_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_cosmo, mf=.true.)
   end subroutine test_cosmo_mf_tangent

   !> COSMO: nuclear gradient
   subroutine test_cosmo_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_cosmo)
   end subroutine test_cosmo_nuclear_gradient

   !> IEF-PCM: surface and trace gradients
   subroutine test_iefpcm_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_iefpcm)
   end subroutine test_iefpcm_surface_gradient

   !> IEF-PCM: surface weight tangent
   subroutine test_iefpcm_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_iefpcm)
   end subroutine test_iefpcm_tangent

   !> IEFPCM: surface weight tangent from matrix-free GMRES, against dense
   subroutine test_iefpcm_mf_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_iefpcm, mf=.true.)
   end subroutine test_iefpcm_mf_tangent

   !> IEF-PCM: nuclear gradient
   subroutine test_iefpcm_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_iefpcm)
   end subroutine test_iefpcm_nuclear_gradient

   !> SS(V)PE: surface and trace gradients
   subroutine test_ssvpe_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_ssvpe)
   end subroutine test_ssvpe_surface_gradient

   !> SS(V)PE: surface weight tangent
   subroutine test_ssvpe_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_ssvpe)
   end subroutine test_ssvpe_tangent

   !> SSVPE: surface weight tangent from matrix-free GMRES, against dense
   subroutine test_ssvpe_mf_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_ssvpe, mf=.true.)
   end subroutine test_ssvpe_mf_tangent

   !> SS(V)PE: nuclear gradient
   subroutine test_ssvpe_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_ssvpe)
   end subroutine test_ssvpe_nuclear_gradient

   !> DPCM: surface and trace gradients
   subroutine test_dpcm_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_dpcm)
   end subroutine test_dpcm_surface_gradient

   !> DPCM: surface weight tangent
   subroutine test_dpcm_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_dpcm)
   end subroutine test_dpcm_tangent

   !> DPCM: surface weight tangent from matrix-free GMRES, against dense
   subroutine test_dpcm_mf_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_dpcm, mf=.true.)
   end subroutine test_dpcm_mf_tangent

   !> DPCM: nuclear gradient
   subroutine test_dpcm_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_dpcm)
   end subroutine test_dpcm_nuclear_gradient

   !> Original DPCM: surface and trace gradients
   subroutine test_dpcm_original_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_surface_gradient

   !> Original DPCM: surface weight tangent
   subroutine test_dpcm_original_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_tangent

   !> DPCM ORIGINAL: surface weight tangent from matrix-free GMRES, against dense
   subroutine test_dpcm_original_mf_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_dpcm_original, mf=.true.)
   end subroutine test_dpcm_original_mf_tangent

   !> Original DPCM: nuclear gradient
   subroutine test_dpcm_original_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_dpcm_original)
   end subroutine test_dpcm_original_nuclear_gradient

   !> LPB: surface and trace gradients
   subroutine test_lpb_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_lpb)
   end subroutine test_lpb_surface_gradient

   !> LPB: surface weight tangent
   subroutine test_lpb_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_lpb)
   end subroutine test_lpb_tangent

   !> LPB: surface weight tangent from matrix-free GMRES, against dense
   subroutine test_lpb_mf_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_lpb, mf=.true.)
   end subroutine test_lpb_mf_tangent

   !> LPB: nuclear gradient
   subroutine test_lpb_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_lpb)
   end subroutine test_lpb_nuclear_gradient

   !> SLIC: surface and trace gradients
   subroutine test_slic_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_slic)
   end subroutine test_slic_surface_gradient

   !> SLIC: surface weight tangent
   subroutine test_slic_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_slic)
   end subroutine test_slic_tangent

   !> SLIC: surface weight tangent from matrix-free GMRES, against dense
   subroutine test_slic_mf_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_slic, mf=.true.)
   end subroutine test_slic_mf_tangent

   !> SLIC: nuclear gradient
   subroutine test_slic_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_slic)
   end subroutine test_slic_nuclear_gradient

   !> Nonlocal: surface and trace gradients
   subroutine test_nonlocal_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_surface_gradient(error, pcm_nonlocal)
   end subroutine test_nonlocal_surface_gradient

   !> Nonlocal: surface weight tangent
   subroutine test_nonlocal_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_nonlocal)
   end subroutine test_nonlocal_tangent

   !> Nonlocal: surface weight tangent from matrix-free GMRES, against dense
   subroutine test_nonlocal_mf_tangent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_tangent(error, pcm_nonlocal, mf=.true.)
   end subroutine test_nonlocal_mf_tangent

   !> Nonlocal: nuclear gradient
   subroutine test_nonlocal_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_nuclear_gradient(error, pcm_nonlocal)
   end subroutine test_nonlocal_nuclear_gradient

   !* ================================================================================= *!
   !*                                      Drivers                                      *!
   !* ================================================================================= *!

   !> Surface and trace gradients on Gaussian and point surfaces of every solute
   !>
   !> @param[out] error     Test failure
   !> @param[in]  kind      Model identifier
   !> @param[in]  iterative Solve with GMRES at the default tolerance (default dense)
   subroutine run_surface_gradient(error, kind, iterative)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind
      !> Solve with GMRES at the default tolerance
      logical, intent(in), optional :: iterative
      !> Solute index
      integer :: imol

      do imol = 1, nmol
         call check_surface_gradient(error, kind, imol, .true., iterative)
         if (allocated(error)) return
         call check_surface_gradient(error, kind, imol, .false., iterative)
         if (allocated(error)) return
      end do
   end subroutine run_surface_gradient

   !> Matrix-free storage with tight GMRES
   !>
   !> Tolerance 1e-13, restart 100; a surface within one 256-point preconditioner
   !> block converges in few iterations, larger ones need the iteration budget
   !>
   !> @param[in,out] pcm Component under test
   subroutine set_matrix_free(pcm)
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm

      pcm%param%solver = solver_type%iterative
      pcm%param%storage = pcm_storage_matrix_free
      pcm%param%solver_tol = 1.0e-13_wp
      pcm%param%solver_maxiter = 2000
      pcm%param%solver_restart = 100
   end subroutine set_matrix_free

   !> Surface weight tangent on both seeds (SLIC: DROP only), Gaussian and point surfaces of every solute
   !> @param[in]  kind  Model identifier
   !> @param[in]  mf    Tangent from the matrix-free GMRES path (default dense)
   subroutine run_tangent(error, kind, mf)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind
      !> Matrix-free GMRES tangent
      logical, intent(in), optional :: mf
      !> Solute and seed indices
      integer :: imol, seed

      do imol = 1, nmol
         do seed = first_seed(kind), seed_drop
            call check_tangent(error, kind, imol, seed, .true., mf)
            if (allocated(error)) return
            call check_tangent(error, kind, imol, seed, .false., mf)
            if (allocated(error)) return
         end do
      end do
   end subroutine run_tangent

   !> Nuclear gradient on rebuilt iSwiG and SvdW-DROP cavities of every solute (SLIC: DROP only)
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

   !> Surface weights and trace adjoints against energy differences
   !>
   !> The host is fixed point charges on the solute atoms, re-evaluated on the
   !> displaced surface, so a surface step moves the traces as well and the
   !> difference is the total derivative `get_gradient_surface_weights` returns
   !>
   !> @param[out] error     Test failure
   !> @param[in]  kind      Model identifier
   !> @param[in]  imol      Solute index
   !> @param[in]  smeared   Gaussian surface (else point surface)
   !> @param[in]  iterative Solve with GMRES at the default tolerance (default dense)
   !> @param[in]  point_potential Couple through point potential integrals (default off)
   subroutine check_surface_gradient(error, kind, imol, smeared, iterative, point_potential)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind
      !> Solute index
      integer, intent(in) :: imol
      !> Gaussian surface
      logical, intent(in) :: smeared
      !> Solve with GMRES at the default tolerance
      logical, intent(in), optional :: iterative
      !> Couple through point potential integrals
      logical, intent(in), optional :: point_potential

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
      !> Absolute part of the difference tolerance
      real(wp) :: atol
      !> Channel and stencil indices, cavity seed
      integer :: ich, k, seed
      !> Failure context
      character(len=200) :: context, prefix
      !> Library error handling
      type(moist_error_type), allocatable :: err

      call new_context(ctx, verbosity=0)
      call new_pcm(kind, ctx, pcm, error)
      if (allocated(error)) return
      atol = grad_atol
      if (present(iterative)) then
         if (iterative) then
            pcm%param%solver = solver_type%iterative
            atol = gmres_atol
         end if
      end if
      if (present(point_potential)) pcm%param%point_potential = point_potential
      call get_solute(imol, mol)
      q = host_charges(mol)
      seed = first_seed(kind)
      call get_seed(seed, smeared, mol, cavity, error)
      if (allocated(error)) return
      h = merge(step_drop_surface, step, seed == seed_drop)
      prefix = trim(pcm%name)//", "//trim(mol_label(imol))//", "//trim(seed_name(seed))// &
         & merge(", Gaussian", ", point   ", smeared)
      if (pcm%param%point_potential) prefix = trim(prefix)//", point potential"
      if (pcm%param%solver == solver_type%iterative) prefix = trim(prefix)//", GMRES"
      call check_potential_request(error, pcm, mol, cavity, smeared .and. .not. pcm%param%point_potential, &
         & trim(prefix))
      if (allocated(error)) return
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
         call check(error, abs(analytic - fd) <= atol + grad_rtol*abs(fd), &
            & message="surface weights disagree with the energy difference", more=trim(context))
         if (allocated(error)) return
         if (expects(active_channels(kind, smeared), ich)) then
            call check(error, abs(fd) > thr_active*abs(energy), message="vacuous channel", &
               & more=trim(context))
            if (allocated(error)) return
         end if
      end do

      ! operators of the unchanged cavity, assembled once for every trace
      call pcm%update(mol, cavity, err)
      call check_moist_error(error, err, "update")
      if (allocated(error)) return
      call host_traces(mol, q, cavity, smeared .and. .not. pcm%param%point_potential, phi, dphi_dr)
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
      call check(error, abs(analytic - fd) <= atol + grad_rtol*abs(fd) .and. &
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
      call check(error, abs(analytic - fd) <= atol + grad_rtol*abs(fd), &
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
   !> @param[in]  mf      Matrix-free tangent against the dense analytic one (default dense against differences)
   subroutine check_tangent(error, kind, imol, seed, smeared, mf)
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
      !> Matrix-free GMRES tangent
      logical, intent(in), optional :: mf

      !> Run context borrowed by the components
      type(moist_context_type), target :: ctx
      !> Component under test and the dense component of the stencil
      class(model_continuum_component_pcm), allocatable :: pcm, ref
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
      !> Matrix-free GMRES tangent
      logical :: matrix_free

      call new_context(ctx, verbosity=0)
      call new_pcm(kind, ctx, pcm, error)
      if (allocated(error)) return
      call new_pcm(kind, ctx, ref, error)
      if (allocated(error)) return
      matrix_free = .false.
      if (present(mf)) matrix_free = mf
      if (matrix_free) call set_matrix_free(pcm)
      call get_solute(imol, mol)
      call get_seed(seed, smeared, mol, cavity, error)
      if (allocated(error)) return
      h = merge(step_drop, step, seed == seed_drop)
      prefix = trim(pcm%name)//", "//trim(mol_label(imol))//", "//trim(seed_name(seed))// &
         & merge(" Gaussian", " point   ", smeared)//merge(" MF", "   ", matrix_free)
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
            call tangent_case(error, pcm, ref, mol, cavity, phi, vn, direction, zero, zero, h, &
               & trim(prefix)//", channel "//trim(channel_name(ich)), "")
            if (allocated(error)) return
         end do

         call clear_direction(direction)
         call tangent_case(error, pcm, ref, mol, cavity, phi, vn, direction, du, zero, h, &
            & trim(prefix)//", phi", "")
         if (allocated(error)) return
         if (pcm%uses_normal_trace()) then
            call tangent_case(error, pcm, ref, mol, cavity, phi, vn, direction, zero, dq, h, &
               & trim(prefix)//", vn", "")
            if (allocated(error)) return
         end if
      end if

      call clear_direction(direction)
      do ich = 1, nchannel
         if (provides(cavity, ich)) call add_channel_direction(cavity, ich, direction)
      end do
      call tangent_case(error, pcm, ref, mol, cavity, phi, vn, direction, du, dq, h, &
         & trim(prefix)//", joint", active_channels(kind, smeared))
   end subroutine check_tangent

   !> One direction of the tangent check
   !>
   !> Dense components are checked against the six-point difference of the
   !> weights at fixed traces. A matrix-free component is checked against the
   !> dense analytic tangent of `ref` at the same state, every weight channel
   !> and adjoint tangent entrywise. A direction without surface motion keeps
   !> the operators of the reference cavity, the stencil only changes the traces
   !>
   !> @param[out]   error     Test failure
   !> @param[in,out] pcm      Component under test
   !> @param[in,out] ref      Dense component, the oracle (stencil or analytic tangent)
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
   subroutine tangent_case(error, pcm, ref, mol, cavity, phi, vn, direction, du, dq, h, prefix, channels)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Dense component, the oracle
      class(model_continuum_component_pcm), intent(inout) :: ref
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
      !> Adjoints on the stencil, adjoint tangents of the dense oracle
      real(wp), allocatable :: s_phi(:, :), s_vn(:, :), dense_phi(:), dense_vn(:)
      !> Tangent accumulator, dense tangent, reference and stencil weights
      type(cavity_surface_adjoint_type) :: tangent, dense, acc0, acc(6)
      !> Stencil index, channel index and number of grid points
      integer :: k, ich, ngrid
      !> Whether the expected channels must be active
      logical :: active
      !> The surface moves along the direction
      logical :: moves

      ngrid = cavity%ngrid
      active = len_trim(channels) > 0
      moves = allocated(direction%w_xyz) .or. allocated(direction%w_n) .or. allocated(direction%w_k1) &
         & .or. allocated(direction%w_k2) .or. allocated(direction%w_xi) .or. allocated(direction%w_f) &
         & .or. allocated(direction%w_a)
      call pcm%update(mol, cavity, err)
      call check_moist_error(error, err, prefix//": update")
      if (allocated(error)) return
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

      if (pcm%param%storage == pcm_storage_matrix_free .or. .not. moves) then
         call ref%update(mol, cavity, err)
         call check_moist_error(error, err, prefix//": update")
         if (allocated(error)) return
      end if
      if (pcm%param%storage == pcm_storage_matrix_free) then
         call fixed_trace_state(ref, mol, cavity, phi, vn, acc0, w_phi, w_vn, error)
         if (allocated(error)) return
         call dense%init(ngrid)
         call ref%surface_weight_tangent(cavity, direction, du, dense, dense_phi, err, &
            & dvn=dq, dw_vn=dense_vn)
         call check_moist_error(error, err, prefix//": dense tangent")
         if (allocated(error)) return
         do ich = 1, size(weight_name)
            call compare_dense(error, prefix, trim(weight_name(ich)), weights_of(tangent, weight_name(ich)), &
               & weights_of(dense, weight_name(ich)), active .and. expects(channels, weight_channel(ich)))
            if (allocated(error)) return
         end do
         call compare_dense(error, prefix, "dE/dphi", dw_phi, dense_phi, active)
         if (allocated(error)) return
         if (pcm%uses_normal_trace()) call compare_dense(error, prefix, "dE/dvn", dw_vn, dense_vn, active)
         return
      end if

      allocate (s_phi(ngrid, 6), s_vn(ngrid, 6))
      trial = cavity
      do k = 1, 6
         if (moves) then
            call displace(cavity, direction, fd6_offsets(k)*h, trial)
            call ref%update(mol, trial, err)
            call check_moist_error(error, err, prefix//": update")
            if (allocated(error)) return
         end if
         call fixed_trace_state(ref, mol, trial, phi + fd6_offsets(k)*h*du, &
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

   !> Compare a matrix-free tangent entrywise with the dense analytic tangent
   !>
   !> Passes when max |mf - dense| <= mf_atol + mf_rtol max |dense| over the
   !> channel; an active channel must not vanish in either
   !>
   !> @param[out] error    Test failure
   !> @param[in]  prefix   Failure context prefix
   !> @param[in]  name     Channel name
   !> @param[in]  analytic Matrix-free tangent
   !> @param[in]  dense    Dense analytic tangent
   !> @param[in]  active   The tangent must not vanish
   subroutine compare_dense(error, prefix, name, analytic, dense, active)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Failure context prefix
      character(len=*), intent(in) :: prefix
      !> Channel name
      character(len=*), intent(in) :: name
      !> Matrix-free tangent
      real(wp), intent(in) :: analytic(:)
      !> Dense analytic tangent
      real(wp), intent(in) :: dense(:)
      !> The tangent must not vanish
      logical, intent(in) :: active

      !> Largest deviation and largest dense entry
      real(wp) :: dev, scale
      !> Failure context
      character(len=200) :: context

      dev = maxval(abs(analytic - dense))
      scale = maxval(abs(dense))
      write (context, '(a,", ",a,": max dev ",es10.2," scale ",es10.2)') prefix, name, dev, scale
      call check(error, dev <= mf_atol + mf_rtol*scale, &
         & message="matrix-free tangent disagrees with the dense tangent", more=trim(context))
      if (allocated(error)) return
      if (active) then
         call check(error, scale > thr_vacuous .and. maxval(abs(analytic)) > thr_vacuous, &
            & message="vacuous tangent", more=trim(context))
      end if
   end subroutine compare_dense

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
   !> @param[in]  kind  Model identifier
   !> @param[in]  imol  Solute index
   !> @param[in]  seed  Cavity kind rebuilt at every geometry
   !> @param[in]  point_potential Couple through point potential integrals (default off)
   subroutine check_nuclear_gradient(error, kind, imol, seed, point_potential)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier
      integer, intent(in) :: kind
      !> Solute index
      integer, intent(in) :: imol
      !> Cavity kind
      integer, intent(in) :: seed
      !> Couple through point potential integrals
      logical, intent(in), optional :: point_potential

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
      !> Failure context prefix
      character(len=:), allocatable :: prefix

      call new_context(ctx, verbosity=0)
      call new_pcm(kind, ctx, pcm, error)
      if (allocated(error)) return
      if (present(point_potential)) pcm%param%point_potential = point_potential
      call get_solute(imol, mol)
      call electron_host(mol, el_xyz, el_q)
      call new_host_cavity(seed, pcm, ctx, cavity, error)
      if (allocated(error)) return
      prefix = trim(pcm%name)//", "//trim(mol_label(imol))//", "//trim(seed_name(seed))
      if (pcm%param%point_potential) prefix = prefix//", point potential"

      call nuclear_state(pcm, mol, el_xyz, el_q, cavity, coupling, energy, error)
      if (allocated(error)) return
      call check_potential_request(error, pcm, mol, cavity, .not. pcm%param%point_potential, prefix)
      if (allocated(error)) return
      call pcm%prepare_gradient(cavity, coupling, err)
      call check_moist_error(error, err, "gradient staging")
      if (allocated(error)) return
      call answer_nuclear_host(cavity, coupling, mol, el_xyz, el_q)
      allocate (gradient(3, mol%nat), source=0.0_wp)
      call pcm%get_gradient(component_view(coupling), cavity, response, gradient, err)
      call check_moist_error(error, err, prefix//": nuclear gradient")
      if (allocated(error)) return
      call check(error, maxval(abs(gradient)) > thr_vacuous, message=prefix// &
         & ": nuclear gradient is vacuously zero")
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
         write (context, '(a,": direction ",i0," dE/dt ",es12.4," vs FD ",es12.4, &
            & ", dev ",es9.2)') prefix, idir, analytic, fd, abs(analytic - fd)
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
   !*                               Component evaluation                                *!
   !* ================================================================================= *!

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
      case (pcm_dpcm, pcm_dpcm_original, pcm_slic)
         channels = merge("r n k xi a", "r n k a   ", smeared)
      case default
         channels = merge("r n k xi f a", "r n k a     ", smeared)
      end select
      channels = trim(channels)
   end function active_channels

   !> Potential request the component stages for the energy phase
   !>
   !> The fixture hosts answer whichever request is staged, so the energy
   !> differences alone cannot tell whether `point_potential` took effect
   !>
   !> @param[out]    error    Test failure
   !> @param[in,out] pcm      Component under test, updated on `cavity`
   !> @param[in]     mol      Solute
   !> @param[in,out] cavity   Cavity, current for `mol`
   !> @param[in]     gaussian Expect `gaussian_potential` (else `point_potential`)
   !> @param[in]     context  Failure context
   subroutine check_potential_request(error, pcm, mol, cavity, gaussian, context)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Cavity
      class(cavity_type), intent(inout) :: cavity
      !> Expect Gaussian potential integrals
      logical, intent(in) :: gaussian
      !> Failure context
      character(len=*), intent(in) :: context

      !> Unanswered energy-phase coupling
      type(coupling_type), target :: coupling
      !> Staged potential requests
      logical :: found_gaussian, found_point

      call stage_energy(pcm, mol, cavity, coupling, error)
      if (allocated(error)) return
      found_gaussian = .false.
      found_point = .false.
      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            found_gaussian = .true.
         type is (point_potential_request_type)
            found_point = .true.
         end select
      end do
      call check(error, (found_gaussian .eqv. gaussian) .and. (found_point .neqv. gaussian), &
         & message="unexpected potential request", more=context// &
         & merge(": expected gaussian_potential", ": expected point_potential   ", gaussian))
   end subroutine check_potential_request

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
   !> The component stays updated on `cavity`: only the traces change
   !>
   !> @param[in,out] pcm     Component under test, updated on `cavity`
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
      call stage_energy(pcm, mol, cavity, coupling, error, keep=.true.)
      if (allocated(error)) return
      call answer_traces(coupling, phi, dphi_dr)
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err, "energy")
   end subroutine trace_energy

   !> Trace adjoints of `get_trace_response` with the traces given directly
   !>
   !> The component stays updated on `cavity`: only the traces change
   !>
   !> @param[in,out] pcm       Component under test, updated on `cavity`
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

      call stage_energy(pcm, mol, cavity, coupling, error, keep=.true.)
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
   !> The component stays updated on `cavity`: only the traces change
   !>
   !> @param[in,out] pcm    Component under test, updated on `cavity`
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
      call stage_energy(pcm, mol, cavity, coupling, error, keep=.true.)
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

   !* ================================================================================= *!
   !*                              Seven-point fixture paths                            *!
   !* ================================================================================= *!

   !> Gaussian component solves, host responses and geometry adjoints
   subroutine test_gaussian_component_paths(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: kind

      do kind = 1, 9
         call check_component_path(kind, .true., error)
         if (allocated(error)) return
      end do
   end subroutine test_gaussian_component_paths

   !> Point component solves, host responses and geometry adjoints
   subroutine test_point_component_paths(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: kind

      do kind = 1, 9
         call check_component_path(kind, .false., error)
         if (allocated(error)) return
      end do
   end subroutine test_point_component_paths

   !* ================================================================================= *!
   !*                        Point potential on Gaussian surfaces                       *!
   !* ================================================================================= *!

   !> Every model: `surface_gradient` on Gaussian surfaces with a point-potential host
   subroutine test_point_potential_surface_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model and solute indices
      integer :: kind, imol

      do kind = 1, npcm
         do imol = 1, nmol
            call check_surface_gradient(error, kind, imol, .true., point_potential=.true.)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_point_potential_surface_gradient

   !> Every model: `nuclear_gradient` on rebuilt Gaussian cavities with a point-potential host
   subroutine test_point_potential_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model, solute and seed indices
      integer :: kind, imol, seed

      do kind = 1, npcm
         do imol = 1, nmol
            do seed = first_seed(kind), seed_drop
               call check_nuclear_gradient(error, kind, imol, seed, point_potential=.true.)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_point_potential_nuclear_gradient

   !> Every model: `component_paths` on the Gaussian fixture with a point-potential host,
   !> including the density-following surface of the response phase
   subroutine test_point_potential_component_paths(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model index
      integer :: kind

      do kind = 1, npcm
         call check_component_path(kind, .true., error, point_potential=.true.)
         if (allocated(error)) return
      end do
   end subroutine test_point_potential_component_paths

   !> Check a component through its public coupling and surface interfaces
   !>
   !> @param[in] kind Component constructor index
   !> @param[in] gaussian Whether the surface has Gaussian widths
   !> @param[out] error Test failure
   !> @param[in] point_potential Couple through point potential integrals (default off)
   subroutine check_component_path(kind, gaussian, error, point_potential)
      !> Component constructor index
      integer, intent(in) :: kind
      !> Whether the surface has Gaussian widths
      logical, intent(in) :: gaussian
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Couple through point potential integrals
      logical, intent(in), optional :: point_potential
      !> The host potential depends on the widths
      logical :: host_gaussian
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(model_continuum_component_pcm), allocatable :: component
      type(cavity_type_drop) :: cavity, displaced
      type(cavity_surface_adjoint_type) :: weights, response_weights
      type(coupling_type), target :: coupling
      type(response_type) :: response
      type(potential_adjoint_response_type) :: trace
      real(wp), parameter :: h = 1.0e-4_wp
      real(wp), parameter :: direction(nfixture) = &
         & [0.02_wp, -0.03_wp, 0.04_wp, -0.01_wp, 0.05_wp, -0.02_wp, 0.03_wp]
      real(wp) :: energy, cached, iterative, samples(4), expected, numerical, shift
      real(wp) :: dxyz(3, nfixture), dn(3, nfixture), dg(3, nfixture)
      real(wp), allocatable :: dense_phi(:), dense_normal(:), iterative_normal(:)
      integer :: channel, offset, i
      character(len=:), allocatable :: label

      call new_context(ctx, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call make_component(kind, component, err, ctx)
      call require_success("constructor")
      if (allocated(error)) return
      label = trim(component%name)//merge(" Gaussian", " point   ", gaussian)
      component%param%solver_tol = 1.0e-12_wp
      component%param%newton_tol = 1.0e-13_wp
      if (present(point_potential)) component%param%point_potential = point_potential
      if (component%param%point_potential) label = label//", point potential"
      host_gaussian = gaussian .and. .not. component%param%point_potential
      cavity%ngrid = nfixture
      cavity%nsph = 1
      cavity%updated = .true.
      cavity%built = .true.
      cavity%xyz = fixture_xyz
      cavity%normal0 = fixture_normals
      cavity%a = fixture_areas
      allocate (cavity%k1(nfixture), source=0.4_wp)
      allocate (cavity%k2(nfixture), source=0.5_wp)
      if (gaussian) then
         cavity%xi0 = fixture_xis
         cavity%f = fixture_fs
      end if
      do i = 1, nfixture
         dxyz(:, i) = direction(i)*[1.0_wp, -0.7_wp, 0.3_wp]
         dn(:, i) = direction(i)*[-0.4_wp, 0.2_wp, 0.8_wp]
         dg(:, i) = direction(i)*[0.3_wp, 0.6_wp, -0.2_wp]
      end do

      call check_potential_request(error, component, mol, cavity, host_gaussian, label)
      if (allocated(error)) return
      call component%update(mol, cavity, err)
      call require_success("update")
      if (allocated(error)) return
      call sample_energy(cavity, 0.0_wp, energy)
      if (allocated(error)) return
      call check(error, ieee_is_finite(energy) .and. abs(energy) > 1.0e-10_wp, &
         & message=label//" nontrivial finite energy")
      if (allocated(error)) return
      call sample_energy(cavity, 0.0_wp, cached)
      if (allocated(error)) return
      call check(error, cached, energy, thr=1.0e-14_wp, message=label//" cached energy")
      if (allocated(error)) return

      call component%prepare_response(cavity, coupling, err)
      call require_success("prepare response")
      if (allocated(error)) return
      call answer_host(cavity, 0.0_wp)
      if (allocated(error)) return
      call component%get_response(component_view(coupling), cavity, response, err)
      call require_success("response")
      if (allocated(error)) return
      call check(error, response%next(), message=label//" potential response exists")
      if (allocated(error)) return
      select type (item => response%item())
      type is (potential_adjoint_response_type)
         trace = item
      class default
         call test_failed(error, label//" unexpected response item")
         return
      end select
      call check(error, allocated(trace%w_phi), message=label//" potential adjoint exists")
      if (allocated(error)) return
      call check(error, allocated(trace%w_dphi_dr) .eqv. component%uses_normal_trace(), &
         & message=label//" normal-trace response requirement")
      if (allocated(error)) return
      dense_phi = trace%w_phi
      call component%normal_trace_adjoint(dense_normal, err)
      call require_success("normal adjoint")
      if (allocated(error)) return
      expected = dot_product(trace%w_phi, direction)
      if (allocated(trace%w_dphi_dr)) expected = expected + sum(trace%w_dphi_dr*dg)
      do offset = 1, 4
         call sample_energy(cavity, h*fd4_offsets(offset), samples(offset))
         if (allocated(error)) return
      end do
      call fd4_scalar(samples(1), samples(2), samples(3), samples(4), h, numerical, error)
      if (allocated(error)) return
      call check(error, numerical, expected, thr=2.0e-8_wp, message=label//" host-trace adjoint")
      if (allocated(error)) return
      expected = 0.0_wp
      if (allocated(trace%w_dphi_dr)) expected = sum(trace%w_dphi_dr*dg)
      do offset = 1, 4
         call sample_energy(cavity, h*fd4_offsets(offset), samples(offset), normal_only=.true.)
         if (allocated(error)) return
      end do
      call fd4_scalar(samples(1), samples(2), samples(3), samples(4), h, numerical, error)
      if (allocated(error)) return
      call check(error, numerical, expected, thr=2.0e-8_wp, message=label//" changed normal trace")
      if (allocated(error)) return

      call sample_energy(cavity, 0.0_wp, cached)
      if (allocated(error)) return
      call component%prepare_gradient(cavity, coupling, err)
      call require_success("prepare gradient")
      if (allocated(error)) return
      call answer_host(cavity, 0.0_wp)
      if (allocated(error)) return
      call weights%init(nfixture)
      call component%get_gradient_surface_weights(component_view(coupling), cavity, weights, err)
      call require_success("surface adjoints")
      if (allocated(error)) return
      if (gaussian) then
         allocate (moist_cavity_drop_lsf_isodensity_internal_type :: cavity%lsf_model)
         call component%prepare_response(cavity, coupling, err)
         call require_success("prepare field-dependent response")
         if (allocated(error)) return
         call answer_host(cavity, 0.0_wp)
         if (allocated(error)) return
         call response_weights%init(nfixture)
         call component%get_surface_weights(component_view(coupling), cavity, response_weights, err)
         call require_success("field-dependent surface adjoints")
         if (allocated(error)) return
         call check(error, maxval(abs(response_weights%w_xyz - weights%w_xyz)), &
            & 0.0_wp, thr=1.0e-12_wp, message=label//" field-dependent position adjoints")
         if (allocated(error)) return
         call check(error, maxval(abs(response_weights%w_xi - weights%w_xi)), &
            & 0.0_wp, thr=1.0e-12_wp, message=label//" field-dependent width adjoints")
         if (allocated(error)) return
         call check(error, maxval(abs(response_weights%w_n - weights%w_n)), &
            & 0.0_wp, thr=1.0e-12_wp, message=label//" field-dependent normal adjoints")
         if (allocated(error)) return
      end if
      do channel = 1, merge(7, 5, gaussian)
         select case (channel)
         case (1)
            expected = sum(weights%w_xyz*dxyz)
         case (2)
            expected = sum(weights%w_n*dn)
         case (3)
            expected = dot_product(weights%w_a, direction)
         case (4)
            expected = dot_product(weights%w_k1, direction)
         case (5)
            expected = dot_product(weights%w_k2, direction)
         case (6)
            expected = dot_product(weights%w_xi, direction)
         case (7)
            expected = dot_product(weights%w_f, direction)
         case default
            call test_failed(error, label//" unknown geometry channel")
            return
         end select
         do offset = 1, 4
            displaced = cavity
            shift = h*fd4_offsets(offset)
            select case (channel)
            case (1)
               displaced%xyz = cavity%xyz + shift*dxyz
            case (2)
               displaced%normal0 = cavity%normal0 + shift*dn
            case (3)
               displaced%a = cavity%a + shift*direction
            case (4)
               displaced%k1 = cavity%k1 + shift*direction
            case (5)
               displaced%k2 = cavity%k2 + shift*direction
            case (6)
               displaced%xi0 = cavity%xi0 + shift*direction
            case (7)
               displaced%f = cavity%f + shift*direction
            case default
               call test_failed(error, label//" unknown geometry channel")
               return
            end select
            call component%update(mol, displaced, err)
            call require_success("displaced update")
            if (allocated(error)) return
            call sample_energy(displaced, 0.0_wp, samples(offset))
            if (allocated(error)) return
         end do
         call fd4_scalar(samples(1), samples(2), samples(3), samples(4), h, numerical, error)
         if (allocated(error)) return
         call check(error, numerical, expected, thr=2.0e-8_wp, &
            & message=label//" geometry channel "//achar(iachar("0") + channel))
         if (allocated(error)) return
      end do

      component%param%solver = solver_type%iterative
      component%param%storage = pcm_storage_matrix_free
      call component%update(mol, cavity, err)
      call require_success("iterative update")
      if (allocated(error)) return
      call sample_energy(cavity, 0.0_wp, iterative)
      if (allocated(error)) return
      call check(error, iterative, energy, thr=2.0e-9_wp, message=label//" dense/iterative energy")
      if (allocated(error)) return
      call component%potential_adjoint(trace%w_phi, err)
      call require_success("iterative adjoint")
      if (allocated(error)) return
      call check(error, maxval(abs(trace%w_phi - dense_phi)), 0.0_wp, thr=2.0e-9_wp, &
         & message=label//" dense/iterative potential adjoint")
      if (allocated(error)) return
      call component%normal_trace_adjoint(iterative_normal, err)
      call require_success("iterative normal adjoint")
      if (allocated(error)) return
      call check(error, maxval(abs(iterative_normal - dense_normal)), 0.0_wp, thr=2.0e-9_wp, &
         & message=label//" dense/iterative normal adjoint")
      if (allocated(error)) return
      call component%prepare_gradient(cavity, coupling, err)
      call require_success("prepare iterative gradient")
      if (allocated(error)) return
      call answer_host(cavity, 0.0_wp)
      if (allocated(error)) return
      call response_weights%init(nfixture)
      call component%get_gradient_surface_weights(component_view(coupling), cavity, response_weights, err)
      call require_success("iterative surface adjoints")
      if (allocated(error)) return
      call check(error, max(maxval(abs(response_weights%w_xyz - weights%w_xyz)), &
         & maxval(abs(response_weights%w_n - weights%w_n)), &
         & maxval(abs(response_weights%w_a - weights%w_a)), &
         & maxval(abs(response_weights%w_k1 - weights%w_k1)), &
         & maxval(abs(response_weights%w_k2 - weights%w_k2)), &
         & maxval(abs(response_weights%w_xi - weights%w_xi)), &
         & maxval(abs(response_weights%w_f - weights%w_f))), &
         & 0.0_wp, thr=2.0e-9_wp, message=label//" dense/iterative surface adjoints")

   contains

      !> Turn a library diagnostic into a model-labelled test failure
      !>
      !> @param[in] phase Operation that produced the diagnostic
      subroutine require_success(phase)
         !> Operation that produced the diagnostic
         character(len=*), intent(in) :: phase
         if (.not. allocated(err)) return
         if (allocated(label)) then
            call test_failed(error, label//" "//phase//": "//err%message)
         else
            call test_failed(error, "PCM component "//phase//": "//err%message)
         end if
      end subroutine require_success

      !> Solve a staged potential, preserving the current assembled operators
      !>
      !> @param[inout] surface Synthetic cavity surface
      !> @param[in] trace_shift Independent host-potential perturbation
      !> @param[out] value Component energy
      !> @param[in] normal_only Change the normal trace while keeping phi fixed
      subroutine sample_energy(surface, trace_shift, value, normal_only)
         !> Synthetic cavity surface
         type(cavity_type_drop), intent(inout) :: surface
         !> Independent host-potential perturbation
         real(wp), intent(in) :: trace_shift
         !> Component energy
         real(wp), intent(out) :: value
         !> Change the normal trace while keeping phi fixed
         logical, intent(in), optional :: normal_only

         call component%new_coupling(surface, coupling, err)
         call require_success("new coupling")
         if (allocated(error)) return
         call component%prepare_energy(surface, coupling, err)
         call require_success("prepare energy")
         if (allocated(error)) return
         call answer_host(surface, trace_shift, normal_only)
         if (allocated(error)) return
         value = 0.0_wp
         call component%get_energy(component_view(coupling), surface, value, err)
         call require_success("energy")
      end subroutine sample_energy

      !> Answer a quadratic host potential and its analytic raw derivatives
      !>
      !> @param[in] surface Synthetic cavity surface
      !> @param[in] trace_shift Independent host-potential perturbation
      !> @param[in] normal_only Change the normal trace while keeping phi fixed
      subroutine answer_host(surface, trace_shift, normal_only)
         !> Synthetic cavity surface
         type(cavity_type_drop), intent(in) :: surface
         !> Independent host-potential perturbation
         real(wp), intent(in) :: trace_shift
         !> Change the normal trace while keeping phi fixed
         logical, intent(in), optional :: normal_only
         real(wp) :: phi(nfixture), grad(3, nfixture), width(nfixture)
         real(wp) :: hessian(3, 3, nfixture), mixed(3, nfixture), dx(3), dxi, base_grad(3), phi_shift
         integer :: point, axis

         hessian = 0.0_wp
         mixed = 0.003_wp
         phi_shift = trace_shift
         if (present(normal_only)) then
            if (normal_only) phi_shift = 0.0_wp
         end if
         do point = 1, nfixture
            dx = surface%xyz(:, point) - fixture_xyz(:, point)
            dxi = 0.0_wp
            if (host_gaussian) dxi = surface%xi0(point) - fixture_xis(point)
            base_grad = 0.01_wp*real(point, wp)*[1.0_wp, -0.5_wp, 0.7_wp]
            do axis = 1, 3
               hessian(axis, axis, point) = 0.002_wp*real(axis + point, wp)
            end do
            width(point) = 0.013_wp*real(point, wp) + dot_product(mixed(:, point), dx)
            phi(point) = (-1.0_wp)**point*0.05_wp*real(point, wp) + &
               & dot_product(base_grad, dx) + 0.5_wp*dot_product(dx, matmul(hessian(:, :, point), dx)) + &
               & width(point)*dxi + phi_shift*direction(point)
            grad(:, point) = base_grad + matmul(hessian(:, :, point), dx) + &
               & mixed(:, point)*dxi + trace_shift*dg(:, point)
         end do
         do while (coupling%next())
            select type (request => coupling%request())
            type is (gaussian_potential_request_type)
               call coupling%answer("phi", phi, err)
               call require_success("submit phi")
               if (allocated(error)) return
               call coupling%answer("dphi_dr", grad, err)
               call require_success("submit dphi_dr")
               if (allocated(error)) return
               call coupling%answer("dphi_dxi", width, err)
               call require_success("submit dphi_dxi")
               if (allocated(error)) return
               call coupling%answer("d2phi_dr2", hessian, err)
               call require_success("submit d2phi_dr2")
               if (allocated(error)) return
               call coupling%answer("d2phi_dr_dxi", mixed, err)
               call require_success("submit d2phi_dr_dxi")
               if (allocated(error)) return
            type is (point_potential_request_type)
               call coupling%answer("phi", phi, err)
               call require_success("submit phi")
               if (allocated(error)) return
               call coupling%answer("dphi_dr", grad, err)
               call require_success("submit dphi_dr")
               if (allocated(error)) return
               call coupling%answer("d2phi_dr2", hessian, err)
               call require_success("submit d2phi_dr2")
               if (allocated(error)) return
            type is (nuclear_charge_request_type)
               call coupling%answer("z", [1.0_wp], err)
               call require_success("submit nuclear charge")
               if (allocated(error)) return
            end select
         end do
      end subroutine answer_host

   end subroutine check_component_path

   !* ================================================================================= *!
   !*                              Parameter derivatives                                *!
   !* ================================================================================= *!

   !> CPCM: parameter derivatives on the Gaussian fixture
   subroutine test_cpcm_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_cpcm, .true., .false., error)
   end subroutine test_cpcm_parameter_derivatives_gaussian

   !> CPCM: parameter derivatives on the point fixture
   subroutine test_cpcm_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_cpcm, .false., .false., error)
   end subroutine test_cpcm_parameter_derivatives_point

   !> COSMO: parameter derivatives on the Gaussian fixture
   subroutine test_cosmo_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_cosmo, .true., .false., error)
   end subroutine test_cosmo_parameter_derivatives_gaussian

   !> COSMO: parameter derivatives on the point fixture
   subroutine test_cosmo_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_cosmo, .false., .false., error)
   end subroutine test_cosmo_parameter_derivatives_point

   !> IEF-PCM: parameter derivatives on the Gaussian fixture
   subroutine test_iefpcm_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_iefpcm, .true., .false., error)
   end subroutine test_iefpcm_parameter_derivatives_gaussian

   !> IEF-PCM: parameter derivatives on the point fixture
   subroutine test_iefpcm_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_iefpcm, .false., .false., error)
   end subroutine test_iefpcm_parameter_derivatives_point

   !> SS(V)PE: parameter derivatives on the Gaussian fixture
   subroutine test_ssvpe_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_ssvpe, .true., .false., error)
   end subroutine test_ssvpe_parameter_derivatives_gaussian

   !> SS(V)PE: parameter derivatives on the point fixture
   subroutine test_ssvpe_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_ssvpe, .false., .false., error)
   end subroutine test_ssvpe_parameter_derivatives_point

   !> DPCM: parameter derivatives on the Gaussian fixture
   subroutine test_dpcm_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_dpcm, .true., .false., error)
   end subroutine test_dpcm_parameter_derivatives_gaussian

   !> DPCM: parameter derivatives on the point fixture
   subroutine test_dpcm_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_dpcm, .false., .false., error)
   end subroutine test_dpcm_parameter_derivatives_point

   !> DPCM original: parameter derivatives on the Gaussian fixture
   subroutine test_dpcm_original_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_dpcm_original, .true., .false., error)
   end subroutine test_dpcm_original_parameter_derivatives_gaussian

   !> DPCM original: parameter derivatives on the point fixture
   subroutine test_dpcm_original_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_dpcm_original, .false., .false., error)
   end subroutine test_dpcm_original_parameter_derivatives_point

   !> LPB: parameter derivatives on the Gaussian fixture
   subroutine test_lpb_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_lpb, .true., .false., error)
   end subroutine test_lpb_parameter_derivatives_gaussian

   !> LPB: parameter derivatives on the point fixture
   subroutine test_lpb_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_lpb, .false., .false., error)
   end subroutine test_lpb_parameter_derivatives_point

   !> SLIC: parameter derivatives on the Gaussian fixture
   subroutine test_slic_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_slic, .true., .false., error)
   end subroutine test_slic_parameter_derivatives_gaussian

   !> SLIC: parameter derivatives on the point fixture
   subroutine test_slic_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_slic, .false., .false., error)
   end subroutine test_slic_parameter_derivatives_point

   !> Nonlocal: parameter derivatives on the Gaussian fixture
   subroutine test_nonlocal_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_nonlocal, .true., .false., error)
   end subroutine test_nonlocal_parameter_derivatives_gaussian

   !> Nonlocal: parameter derivatives on the point fixture
   subroutine test_nonlocal_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_nonlocal, .false., .false., error)
   end subroutine test_nonlocal_parameter_derivatives_point

   !> SLIC with `mu` left at its default, which follows alpha and gamma, Gaussian fixture
   subroutine test_slic_default_mu_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_slic, .true., .true., error)
   end subroutine test_slic_default_mu_parameter_derivatives_gaussian

   !> SLIC with `mu` left at its default, which follows alpha and gamma, point fixture
   subroutine test_slic_default_mu_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_slic, .false., .true., error)
   end subroutine test_slic_default_mu_parameter_derivatives_point

   !> CPCM matrix-free parameter derivatives, Gaussian fixture
   subroutine test_cpcm_mf_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_cpcm, .true., .false., error, mf=.true.)
   end subroutine test_cpcm_mf_parameter_derivatives_gaussian

   !> CPCM matrix-free parameter derivatives, point fixture
   subroutine test_cpcm_mf_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_cpcm, .false., .false., error, mf=.true.)
   end subroutine test_cpcm_mf_parameter_derivatives_point

   !> COSMO matrix-free parameter derivatives, Gaussian fixture
   subroutine test_cosmo_mf_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_cosmo, .true., .false., error, mf=.true.)
   end subroutine test_cosmo_mf_parameter_derivatives_gaussian

   !> COSMO matrix-free parameter derivatives, point fixture
   subroutine test_cosmo_mf_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_cosmo, .false., .false., error, mf=.true.)
   end subroutine test_cosmo_mf_parameter_derivatives_point

   !> IEFPCM matrix-free parameter derivatives, Gaussian fixture
   subroutine test_iefpcm_mf_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_iefpcm, .true., .false., error, mf=.true.)
   end subroutine test_iefpcm_mf_parameter_derivatives_gaussian

   !> IEFPCM matrix-free parameter derivatives, point fixture
   subroutine test_iefpcm_mf_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_iefpcm, .false., .false., error, mf=.true.)
   end subroutine test_iefpcm_mf_parameter_derivatives_point

   !> SSVPE matrix-free parameter derivatives, Gaussian fixture
   subroutine test_ssvpe_mf_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_ssvpe, .true., .false., error, mf=.true.)
   end subroutine test_ssvpe_mf_parameter_derivatives_gaussian

   !> SSVPE matrix-free parameter derivatives, point fixture
   subroutine test_ssvpe_mf_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_ssvpe, .false., .false., error, mf=.true.)
   end subroutine test_ssvpe_mf_parameter_derivatives_point

   !> DPCM matrix-free parameter derivatives, Gaussian fixture
   subroutine test_dpcm_mf_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_dpcm, .true., .false., error, mf=.true.)
   end subroutine test_dpcm_mf_parameter_derivatives_gaussian

   !> DPCM matrix-free parameter derivatives, point fixture
   subroutine test_dpcm_mf_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_dpcm, .false., .false., error, mf=.true.)
   end subroutine test_dpcm_mf_parameter_derivatives_point

   !> DPCM_ORIGINAL matrix-free parameter derivatives, Gaussian fixture
   subroutine test_dpcm_original_mf_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_dpcm_original, .true., .false., error, mf=.true.)
   end subroutine test_dpcm_original_mf_parameter_derivatives_gaussian

   !> DPCM_ORIGINAL matrix-free parameter derivatives, point fixture
   subroutine test_dpcm_original_mf_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_dpcm_original, .false., .false., error, mf=.true.)
   end subroutine test_dpcm_original_mf_parameter_derivatives_point

   !> LPB matrix-free parameter derivatives, Gaussian fixture
   subroutine test_lpb_mf_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_lpb, .true., .false., error, mf=.true.)
   end subroutine test_lpb_mf_parameter_derivatives_gaussian

   !> LPB matrix-free parameter derivatives, point fixture
   subroutine test_lpb_mf_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_lpb, .false., .false., error, mf=.true.)
   end subroutine test_lpb_mf_parameter_derivatives_point

   !> SLIC matrix-free parameter derivatives, Gaussian fixture
   subroutine test_slic_mf_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_slic, .true., .false., error, mf=.true.)
   end subroutine test_slic_mf_parameter_derivatives_gaussian

   !> SLIC matrix-free parameter derivatives, point fixture
   subroutine test_slic_mf_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_slic, .false., .false., error, mf=.true.)
   end subroutine test_slic_mf_parameter_derivatives_point

   !> NONLOCAL matrix-free parameter derivatives, Gaussian fixture
   subroutine test_nonlocal_mf_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_nonlocal, .true., .false., error, mf=.true.)
   end subroutine test_nonlocal_mf_parameter_derivatives_gaussian

   !> NONLOCAL matrix-free parameter derivatives, point fixture
   subroutine test_nonlocal_mf_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_nonlocal, .false., .false., error, mf=.true.)
   end subroutine test_nonlocal_mf_parameter_derivatives_point

   !> SLIC default `mu` matrix-free parameter derivatives, Gaussian fixture
   subroutine test_slic_default_mu_mf_parameter_derivatives_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_slic, .true., .true., error, mf=.true.)
   end subroutine test_slic_default_mu_mf_parameter_derivatives_gaussian

   !> SLIC default `mu` matrix-free parameter derivatives, point fixture
   subroutine test_slic_default_mu_mf_parameter_derivatives_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_parameter_derivatives(pcm_slic, .false., .true., error, mf=.true.)
   end subroutine test_slic_default_mu_mf_parameter_derivatives_point

   !> Parameter derivatives of the response against energy differences
   !>
   !> Six-point central differences of the energy at fixed cavity and host
   !> potential in every physical parameter the model has; the response item
   !> holds exactly those parameters and no other
   !>
   !> @param[in] kind Component constructor index
   !> @param[in] gaussian Whether the surface has Gaussian widths
   !> @param[in] mu_default SLIC only: construct without `mu`
   !> @param[out] error Test failure
   !> @param[in] mf Analytic side on matrix-free GMRES storage, energy differences stay dense
   subroutine check_parameter_derivatives(kind, gaussian, mu_default, error, mf)
      !> Component constructor index
      integer, intent(in) :: kind
      !> Whether the surface has Gaussian widths
      logical, intent(in) :: gaussian
      !> SLIC without explicit mu
      logical, intent(in) :: mu_default
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Matrix-free analytic side
      logical, intent(in), optional :: mf
      type(moist_error_type), allocatable :: err
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(model_continuum_component_pcm), allocatable :: component
      type(cavity_type_drop) :: cavity
      type(coupling_type), target :: coupling
      type(response_type) :: response
      type(parameter_adjoint_response_type) :: item
      type(pcm_physics) :: base, phys
      real(wp) :: samples(6), analytic, numerical, h, energy0
      logical :: present_ok, found
      integer :: ip, offset
      character(len=:), allocatable :: label

      label = "PCM parameter derivatives"
      call new_context(ctx, verbosity=0)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      mol%charge = 1.0_wp
      cavity%ngrid = nfixture
      cavity%nsph = 1
      cavity%updated = .true.
      cavity%built = .true.
      cavity%xyz = fixture_xyz
      cavity%normal0 = fixture_normals
      cavity%a = fixture_areas
      allocate (cavity%k1(nfixture), source=0.4_wp)
      allocate (cavity%k2(nfixture), source=0.5_wp)
      if (gaussian) then
         cavity%xi0 = fixture_xis
         cavity%f = fixture_fs
      end if
      base%eps_solute = 1.5_wp

      ! Analytic derivatives, requested through the parameter
      call build(base, .true., component)
      if (allocated(error)) return
      if (present(mf)) then
         if (mf) call set_matrix_free(component)
      end if
      label = trim(component%name)//merge(" Gaussian", " point   ", gaussian)// &
         & merge(" default mu", "           ", mu_default)
      if (present(mf)) then
         if (mf) label = label//" matrix-free"
      end if
      call component%update(mol, cavity, err)
      call require_success("update")
      if (allocated(error)) return
      call component%prepare_response(cavity, coupling, err)
      call require_success("prepare response")
      if (allocated(error)) return
      call answer_host()
      if (allocated(error)) return
      call component%get_response(component_view(coupling), cavity, response, err)
      call require_success("response")
      if (allocated(error)) return
      found = .false.
      do while (response%next())
         select type (it => response%item())
         type is (parameter_adjoint_response_type)
            item = it
            found = .true.
         end select
      end do
      call check(error, found, message=label//" parameter item present")
      if (allocated(error)) return
      call energy_at(base, energy0)
      if (allocated(error)) return

      do ip = 1, 10
         present_ok = applicable(ip)
         call check(error, item_has(ip), present_ok, message=label//" parameter set: "//pname(ip))
         if (allocated(error)) return
         if (.not. present_ok) cycle
         h = 1.0e-3_wp*max(1.0_wp, abs(value_of(base, ip)))
         do offset = 1, 6
            phys = base
            call shift(phys, ip, h*fd6_offsets(offset))
            call energy_at(phys, samples(offset))
            if (allocated(error)) return
         end do
         call fd6_scalar(samples(1), samples(2), samples(3), samples(4), samples(5), samples(6), &
            & h, numerical, error)
         if (allocated(error)) return
         analytic = item_value(ip)
         call check(error, analytic, numerical, thr=param_atol*max(1.0_wp, abs(energy0)) + &
            & param_rtol*abs(numerical), message=label//" d/d"//pname(ip))
         if (allocated(error)) return
      end do

   contains

      !> Parameter names, in the order of the response item
      function pname(i_) result(name)
         !> Parameter number
         integer, intent(in) :: i_
         character(len=10) :: name
         character(len=10), parameter :: names(10) = [character(len=10) :: "eps", "eps_inf", &
            & "lambda", "eps_solute", "kappa", "alpha", "beta", "gamma", "mu", "phi_static"]

         name = names(i_)
      end function pname

      !> Whether the model has the parameter as an independent constructor argument
      function applicable(i_) result(yes)
         !> Parameter number
         integer, intent(in) :: i_
         logical :: yes

         select case (kind)
         case (pcm_lpb)
            yes = i_ == 1 .or. i_ == 5
         case (pcm_slic)
            yes = i_ == 1 .or. (i_ >= 6 .and. i_ <= 8) .or. i_ == 10 .or. (i_ == 9 .and. .not. mu_default)
         case (pcm_nonlocal)
            yes = i_ <= 4
         case default
            yes = i_ == 1
         end select
      end function applicable

      !> Whether the response item carries the parameter
      function item_has(i_) result(yes)
         !> Parameter number
         integer, intent(in) :: i_
         logical :: yes

         select case (i_)
         case (1); yes = allocated(item%dg_deps)
         case (2); yes = allocated(item%dg_deps_inf)
         case (3); yes = allocated(item%dg_dlambda)
         case (4); yes = allocated(item%dg_deps_solute)
         case (5); yes = allocated(item%dg_dkappa)
         case (6); yes = allocated(item%dg_dalpha)
         case (7); yes = allocated(item%dg_dbeta)
         case (8); yes = allocated(item%dg_dgamma)
         case (9); yes = allocated(item%dg_dmu)
         case default; yes = allocated(item%dg_dphi_static)
         end select
      end function item_has

      !> Value of a parameter in the response item
      function item_value(i_) result(v)
         !> Parameter number
         integer, intent(in) :: i_
         real(wp) :: v

         select case (i_)
         case (1); v = item%dg_deps
         case (2); v = item%dg_deps_inf
         case (3); v = item%dg_dlambda
         case (4); v = item%dg_deps_solute
         case (5); v = item%dg_dkappa
         case (6); v = item%dg_dalpha
         case (7); v = item%dg_dbeta
         case (8); v = item%dg_dgamma
         case (9); v = item%dg_dmu
         case default; v = item%dg_dphi_static
         end select
      end function item_value

      !> Value of a parameter of the physics record
      function value_of(ph, i_) result(v)
         !> Physics
         type(pcm_physics), intent(in) :: ph
         !> Parameter number
         integer, intent(in) :: i_
         real(wp) :: v

         select case (i_)
         case (1); v = ph%eps
         case (2); v = ph%eps_inf
         case (3); v = ph%lambda
         case (4); v = ph%eps_solute
         case (5); v = ph%kappa
         case (6); v = ph%slic_alpha
         case (7); v = ph%slic_beta
         case (8); v = ph%slic_gamma
         case (9); v = ph%slic_mu
         case default; v = ph%slic_phi_static
         end select
      end function value_of

      !> Shift one parameter of the physics record
      subroutine shift(ph, i_, d)
         !> Physics
         type(pcm_physics), intent(inout) :: ph
         !> Parameter number
         integer, intent(in) :: i_
         !> Shift
         real(wp), intent(in) :: d

         select case (i_)
         case (1); ph%eps = ph%eps + d
         case (2); ph%eps_inf = ph%eps_inf + d
         case (3); ph%lambda = ph%lambda + d
         case (4); ph%eps_solute = ph%eps_solute + d
         case (5); ph%kappa = ph%kappa + d
         case (6); ph%slic_alpha = ph%slic_alpha + d
         case (7); ph%slic_beta = ph%slic_beta + d
         case (8); ph%slic_gamma = ph%slic_gamma + d
         case (9); ph%slic_mu = ph%slic_mu + d
         case default; ph%slic_phi_static = ph%slic_phi_static + d
         end select
      end subroutine shift

      !> Component for a physics record
      subroutine build(ph, derivatives, comp)
         !> Physics
         type(pcm_physics), intent(in) :: ph
         !> Request the parameter derivatives
         logical, intent(in) :: derivatives
         !> Component
         class(model_continuum_component_pcm), allocatable, intent(out) :: comp
         type(moist_pcm_parameters_type) :: param

         if (kind == pcm_slic .and. mu_default) then
            param%solver = solver_type%dense
            param%newton_tol = 1.0e-13_wp
            allocate (model_continuum_component_slic :: comp)
            select type (comp)
            type is (model_continuum_component_slic)
               call new_component_slic(comp, ph%eps, ph%slic_alpha, ph%slic_beta, ph%slic_gamma, &
                  & err, param, phi_static=ph%slic_phi_static, ctx=ctx)
            end select
            call require_success("constructor")
         else
            call new_pcm(kind, ctx, comp, error, ph)
         end if
         if (allocated(error)) return
         comp%param%parameter_derivatives = derivatives
      end subroutine build

      !> Energy of a physics record at the fixed cavity and host potential
      subroutine energy_at(ph, value)
         !> Physics
         type(pcm_physics), intent(in) :: ph
         !> Energy
         real(wp), intent(out) :: value
         class(model_continuum_component_pcm), allocatable :: comp

         value = 0.0_wp
         call build(ph, .false., comp)
         if (allocated(error)) return
         call comp%update(mol, cavity, err)
         call require_success("update")
         if (allocated(error)) return
         call comp%new_coupling(cavity, coupling, err)
         call require_success("new coupling")
         if (allocated(error)) return
         call comp%prepare_energy(cavity, coupling, err)
         call require_success("prepare energy")
         if (allocated(error)) return
         call answer_host()
         if (allocated(error)) return
         call comp%get_energy(component_view(coupling), cavity, value, err)
         call require_success("energy")
      end subroutine energy_at

      !> Answer a smooth host potential and its first derivatives
      subroutine answer_host()
         real(wp) :: phi(nfixture), grad(3, nfixture), width(nfixture)
         integer :: point

         do point = 1, nfixture
            phi(point) = (-1.0_wp)**point*0.05_wp*real(point, wp)
            grad(:, point) = 0.01_wp*real(point, wp)*[1.0_wp, -0.5_wp, 0.7_wp]
            width(point) = 0.013_wp*real(point, wp)
         end do
         do while (coupling%next())
            select type (request => coupling%request())
            type is (gaussian_potential_request_type)
               call coupling%answer("phi", phi, err)
               call require_success("submit phi")
               if (allocated(error)) return
               call coupling%answer("dphi_dr", grad, err)
               call require_success("submit dphi_dr")
               if (allocated(error)) return
               call coupling%answer("dphi_dxi", width, err)
               call require_success("submit dphi_dxi")
               if (allocated(error)) return
            type is (point_potential_request_type)
               call coupling%answer("phi", phi, err)
               call require_success("submit phi")
               if (allocated(error)) return
               call coupling%answer("dphi_dr", grad, err)
               call require_success("submit dphi_dr")
               if (allocated(error)) return
            type is (nuclear_charge_request_type)
               call coupling%answer("z", [1.0_wp], err)
               call require_success("submit nuclear charge")
               if (allocated(error)) return
            end select
         end do
      end subroutine answer_host

      !> Turn a library error into a test failure
      subroutine require_success(stage)
         !> Stage name
         character(len=*), intent(in) :: stage

         if (allocated(err)) then
            call test_failed(error, label//" "//stage//": "//err%message)
            deallocate (err)
         end if
      end subroutine require_success

   end subroutine check_parameter_derivatives

   !* ================================================================================= *!
   !*                                 CPCM host features                                *!
   !* ================================================================================= *!

   !> Contract the usual iSwiG adjoint and the manufactured width channel
   !>
   !> @param[in] self Test cavity, optionally carrying xi1_rA
   !> @param[in] acc Surface weights (ngrid)
   !> @param[in,out] gradient Nuclear-gradient accumulator (3, nat)
   !> @param[out] error Allocated on failure
   subroutine fixture_surface_gradient(self, acc, gradient, error)
      !> Test cavity
      class(width_gradient_cavity), intent(in) :: self
      !> Accumulated surface adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error
      integer :: igrid

      call self%cavity_type_iswig%contract_surface_gradient(acc, gradient, error)
      if (allocated(error)) return
      if (.not. allocated(self%xi1_rA)) return
      do igrid = 1, self%ngrid
         gradient = gradient + acc%w_xi(igrid)*self%xi1_rA(:, :, igrid)
      end do
   end subroutine fixture_surface_gradient

!> Component surface adjoints on the shared synthetic seven-point surface
!>
!> Drives `get_surface_weights` on a hand-built DROP surface and differences it
!> with the shared 4-point harness against an independently coded CPCM energy
   subroutine test_cpcm_surface_weights(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      type(structure_type) :: mol
      type(model_continuum_component_cpcm) :: pcm_model
      type(cavity_type_drop) :: cavity
      type(cavity_surface_adjoint_type) :: weights
      type(surface_fixture) :: surface
      type(coupling_type) :: coupling
      !> Energy of the fixture, solved to stage the response phase
      real(wp) :: energy
      !> Nonconstant host potential and its derivatives
      real(wp) :: dphi_dr(3, nfixture), dphi_dxi(nfixture)
      logical :: moving_potential

      !> Dielectric constant
      real(wp), parameter :: epsilon = 78.4_wp
      !> Dielectric scaling shared with the reference energy
      real(wp), parameter :: feps_ref = (epsilon - 1.0_wp)/epsilon

      !> Dummy molecular geometry
      real(wp) :: xyz_mol(3, 1)
      integer :: igrid, iaxis, icoord
      real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :)
      real(wp) :: grad_rA(3, 1), grad_ref(3, 1), gradient_deviation

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      moving_potential = .false.
      call new_context(ctx, nthreads=0)

      ! A single dummy center: with a host-supplied potential the molecular
      ! geometry never enters the energy, but update() stores it
      xyz_mol(:, 1) = 0.0_wp
      call new(mol, [1], xyz_mol)

      ! Synthetic surface carrying only the fields the CPCM matrix reads
      cavity%nsph = mol%nat
      cavity%ngrid = nfixture
      allocate (cavity%a, source=fixture_areas)
      ! The hand-built surface counts as updated and built
      cavity%updated = .true.
      cavity%built = .true.
      allocate (cavity%xi0, source=fixture_xis)
      allocate (cavity%f, source=fixture_fs)
      allocate (cavity%xyz, source=fixture_xyz)
      allocate (cavity%normal0, source=fixture_normals)

      call new_component_cpcm(pcm_model, epsilon=epsilon, error=err, &
                              param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      if (allocated(err)) then
         call test_failed(error, "CPCM initialization failed: "//err%message)
         return
      end if

      call pcm_model%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if

      ! The hand-built DROP surface responds to the field, so the response
      ! phase asks for the host surface weights; a fixed host potential has
      ! none, and the zeros are answered explicitly
      call pcm_model%new_coupling(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Coupling setup failed: "//err%message)
         return
      end if
      call pcm_model%prepare_energy(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Energy staging failed: "//err%message)
         return
      end if
      call set_host_potential(coupling, sw_phi)
      energy = 0.0_wp
      call pcm_model%get_energy(component_view(coupling), cavity, energy, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM energy failed: "//err%message)
         return
      end if
      call pcm_model%prepare_response(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Response staging failed: "//err%message)
         return
      end if
      call fill_missing_with_zeros(cavity, coupling)

      call weights%init(nfixture)
      call pcm_model%get_surface_weights(component_view(coupling), cavity, weights, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM surface-weight assembly failed: "//err%message)
         return
      end if

      call new_surface_fixture(surface, fixture_areas, fixture_xis, fixture_fs, fixture_xyz, fixture_normals)
      call check_surface_weights(error, surface, weights, cpcm_surface_energy, &
         & "cpcm", step=1.0e-4_wp, thr_abs=2.0e-11_wp, thr_rel=2.0e-9_wp, &
         & check_normal=.false.)
      if (allocated(error)) return

      allocate (w_xi(nfixture), w_f(nfixture), w_xyz(3, nfixture))
      call probe_conductor_bilinear_weights(cavity%xi0, cavity%f, cavity%xyz, &
         & pcm_model%q, pcm_model%q, w_xi, w_f, w_xyz, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM A-matrix surface weights failed: "//err%message)
         return
      end if
      call check(error, maxval(abs(weights%w_xi - 0.5_wp*w_xi/feps_ref)), &
         & 0.0_wp, thr=thr, message="type-bound xi weights disagree with PCM weights")
      if (allocated(error)) return
      call check(error, maxval(abs(weights%w_f - 0.5_wp*w_f/feps_ref)), &
         & 0.0_wp, thr=thr, message="type-bound f weights disagree with PCM weights")
      if (allocated(error)) return
      call check(error, maxval(abs(weights%w_xyz - 0.5_wp*w_xyz/feps_ref)), &
         & 0.0_wp, thr=thr, message="type-bound xyz weights disagree with PCM weights")
      if (allocated(error)) return

      allocate (cavity%xi1_rA(3, 1, nfixture))
      allocate (cavity%f1_rA(3, 1, nfixture))
      allocate (cavity%xyz1_rA(3, 3, 1, nfixture))
      do igrid = 1, nfixture
         do iaxis = 1, 3
            cavity%xi1_rA(iaxis, 1, igrid) = 0.01_wp*real(iaxis + igrid, wp)
            cavity%f1_rA(iaxis, 1, igrid) = -0.02_wp*real(2*iaxis + igrid, wp)
            do icoord = 1, 3
               cavity%xyz1_rA(icoord, iaxis, 1, igrid) = &
                  & 0.005_wp*real(icoord - iaxis + igrid, wp)
            end do
         end do
      end do

      call probe_conductor_nuclear_gradient(cavity%xi1_rA, cavity%f1_rA, cavity%xyz1_rA, &
         & w_xi, w_f, w_xyz, grad_rA, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM A-matrix nuclear gradient failed: "//err%message)
         return
      end if

      grad_ref = 0.0_wp
      do igrid = 1, nfixture
         do iaxis = 1, 3
            grad_ref(iaxis, 1) = grad_ref(iaxis, 1) &
               & + w_xi(igrid)*cavity%xi1_rA(iaxis, 1, igrid) &
               & + w_f(igrid)*cavity%f1_rA(iaxis, 1, igrid)
            do icoord = 1, 3
               grad_ref(iaxis, 1) = grad_ref(iaxis, 1) &
                  & + w_xyz(icoord, igrid)*cavity%xyz1_rA(icoord, iaxis, 1, igrid)
            end do
         end do
      end do
      gradient_deviation = maxval(abs(grad_rA - grad_ref))
      call check(error, gradient_deviation, 0.0_wp, thr=thr, &
         & message="type-bound A-matrix nuclear gradient disagrees with surface contraction")

      if (allocated(error)) return

      ! Test raw host derivatives in the response phase on a density-dependent surface
      allocate (moist_cavity_drop_lsf_isodensity_internal_type :: cavity%lsf_model)
      call pcm_model%prepare_response(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      do igrid = 1, nfixture
         dphi_dr(:, igrid) = [0.03_wp, -0.02_wp, 0.04_wp]*real(igrid, wp)
         dphi_dxi(igrid) = -0.01_wp*real(igrid, wp)
      end do
      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            call coupling%answer("dphi_dr", dphi_dr, err)
            if (.not. allocated(err)) call coupling%answer("dphi_dxi", dphi_dxi, err)
         end select
         if (allocated(err)) then
            call test_failed(error, "Answering the raw derivatives failed: "//err%message)
            return
         end if
      end do
      call weights%init(nfixture)
      call pcm_model%get_surface_weights(component_view(coupling), cavity, weights, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      moving_potential = .true.
      call check_surface_weights(error, surface, weights, cpcm_surface_energy, &
         & "cpcm raw potential response", step=1.0e-4_wp, thr_abs=2.0e-11_wp, thr_rel=2.0e-9_wp, &
         & check_normal=.false.)

   contains

      !> Rebuild the CPCM energy independently on a perturbed surface
      function cpcm_surface_energy(trial) result(energy)
         !> Perturbed surface fixture
         type(surface_fixture), intent(in) :: trial
         !> CPCM polarization energy
         real(wp) :: energy

         real(wp) :: amat_trial(nfixture, nfixture)
         real(wp) :: q_trial(nfixture), rhs(nfixture), phi_trial(nfixture)
         real(wp) :: xi_i, xi_j, xi_ij, r_ij
         integer :: i, j
         type(moist_error_type), allocatable :: solve_error

         amat_trial = 0.0_wp
         do i = 1, nfixture
            amat_trial(i, i) = trial%xi(i)*sqrt(2.0_wp/pi)/trial%f(i)
            do j = 1, i - 1
               xi_i = trial%xi(i)
               xi_j = trial%xi(j)
               xi_ij = xi_i*xi_j/sqrt(xi_i*xi_i + xi_j*xi_j)
               r_ij = norm2(trial%xyz(:, i) - trial%xyz(:, j))
               amat_trial(i, j) = erf(xi_ij*r_ij)/r_ij
               amat_trial(j, i) = amat_trial(i, j)
            end do
         end do

         phi_trial = sw_phi
         if (moving_potential) then
            do i = 1, nfixture
               phi_trial(i) = phi_trial(i) &
                  & + dot_product(dphi_dr(:, i), trial%xyz(:, i) - fixture_xyz(:, i)) &
                  & + dphi_dxi(i)*(trial%xi(i) - fixture_xis(i))
            end do
         end if
         rhs(:) = -feps_ref*phi_trial
         call solve_pcm_lu(amat_trial, rhs, q_trial, solve_error)
         if (allocated(solve_error)) error stop "cpcm_surface_energy: CPCM solve failed"

         energy = 0.5_wp*dot_product(q_trial, phi_trial)
      end function cpcm_surface_energy

   end subroutine test_cpcm_surface_weights

!> Check the CPCM nuclear gradient of an ECP host against finite differences
!>
!> The host potential is that of the reduced charges `z`, so the direct term
!> must use the `z` the host answers in `nuclear_charges`. The driver also pins
!> forward against reverse gradients and the raw host width contraction on a
!> manufactured xi1_rA; the all-electron differences are in `cpcm_nuclear_gradient`
   subroutine test_cpcm_nuclear_gradient_ecp(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_cpcm_nuclear_gradient(error, core=2)
   end subroutine test_cpcm_nuclear_gradient_ecp

!> Check the complete CPCM nuclear gradient against finite differences
!>
!> This exercises the component entry point end to end with the point-charge
!> potential of the nuclear charges supplied through the coupling trace and
!> re-evaluated on every displaced cavity. Both terms in
!>
!>    dE/dR_A = q^T dphi/dR_A
!>              + 1/(2 f(eps)) q^T (dA/dR_A) q
!>
!> are active because rebuilding the iSwiG cavity moves its  grid points and changes
!> its switching factors while the source charges move with the nuclei
!>
!> @param[out] error Test failure
!> @param[in]  core  ECP core electrons on one displaced atom (0: all-electron)
   subroutine check_cpcm_nuclear_gradient(error, core)

      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> ECP core electrons on one displaced atom
      integer, intent(in) :: core

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Reference and displaced molecular structures
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: trial
      !> PCM component and molecular cavity
      type(model_continuum_component_cpcm) :: pcm_model
      type(width_gradient_cavity) :: cavity
      !> Point-charge potential trace, refilled on every displaced cavity
      type(coupling_type) :: coupling
      !> Point charges generating the potential: the host's nuclear charges, so
      !> the analytic gradient's moving sources coincide with the potential's
      real(wp), allocatable :: qat(:)
      !> Nuclear charges of the ECP variant
      integer, allocatable :: zat(:)
      !> Atom index used to fill the charges
      integer :: iat
      !> Analytic PCM nuclear gradient
      real(wp), allocatable :: gradient(:, :), forward_gradient(:, :)
      !> Manufactured width-gradient probes and raw host width derivative
      real(wp), allocatable :: width_gradient(:, :, :), width_derivative(:), width_forward(:, :)
      !> Surface point and probe index for the width contraction
      integer :: width_point, iprobe
      !> Saved cavity width derivative at the probe point and expected probe difference
      real(wp), allocatable :: saved_xi1(:, :), width_expected(:, :)
      !> Probe-difference residual and its tolerance
      real(wp) :: width_residual, width_thr
      !> Energy of the reference geometry, solved to stage the gradient phase
      real(wp) :: energy
      !> Host part of the gradient phase, unused here
      type(response_type) :: response
      !> Four-point finite-difference energies
      real(wp) :: values(4)
      !> Finite-difference derivative and saved coordinate
      real(wp) :: fd, saved
      !> Atom, Cartesian-axis, and stencil indices
      integer :: iatom, iaxis, k
      !> Failure context
      character(len=96) :: context

      !> Dielectric constant
      real(wp), parameter :: epsilon = 78.4_wp
      !> Lebedev grid size
      integer, parameter :: nleb = 26
      !> Nuclear-displacement step in bohr
      real(wp), parameter :: step = 1.0e-4_wp
      !> Tight tolerances accounting for cavity and linear-solve roundoff
      real(wp), parameter :: fd_atol = 5.0e-9_wp
      real(wp), parameter :: fd_rtol = 5.0e-8_wp

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      call get_test_structures(mols, 5)
      call center_at_origin(mols(1))
      zat = mols(1)%num(mols(1)%id)
      if (core > 0) then
         ! The core goes on an atom the finite differences below displace: only
         ! its own direct term carries its charge
         iat = findloc(zat(:min(2, mols(1)%nat)) > core, .true., dim=1)
         call check(error, iat > 0, more="ECP variant needs a displaced atom heavier than its core")
         if (allocated(error)) return
         zat(iat) = zat(iat) - core
      end if
      qat = real(zat, wp)

      call get_test_cavity_iswig(mols(1), cavity%cavity_type_iswig, err, nleb=nleb)
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
      call pcm_model%update(mols(1), cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if
      call stage_point_charge_energy(error, pcm_model, cavity, qat, mols(1), coupling)
      if (allocated(error)) return
      energy = 0.0_wp
      call pcm_model%get_energy(component_view(coupling), cavity, energy, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM energy failed: "//err%message)
         return
      end if
      ! The gradient phase asks for the host total position weight
      ! `q_i grad phi(r_i)` of the point-charge potential (the width weight,
      ! zero for a width-independent potential, is answered with zeros)
      call pcm_model%prepare_gradient(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Gradient staging failed: "//err%message)
         return
      end if
      call fill_missing_with_zeros(cavity, coupling)
      call fill_point_charge_field(cavity, coupling, qat, mols(1))

      allocate (gradient(3, mols(1)%nat), source=0.0_wp)
      call pcm_model%get_gradient(component_view(coupling), cavity, response, gradient, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM nuclear gradient failed: "//err%message)
         return
      end if
      call check(error, maxval(abs(gradient)) > 1.0e-8_wp, &
         & more="CPCM nuclear gradient is zero, the test is vacuous")
      if (allocated(error)) return

      call check(error,.not. allocated(cavity%xyz1_rA), &
                 more="reverse CPCM gradient built a surface nuclear Jacobian")
      if (allocated(error)) return
      allocate (forward_gradient(3, mols(1)%nat), source=0.0_wp)
      call pcm_model%get_forward_gradient(component_view(coupling), cavity, response, forward_gradient, err)
      call check(error,.not. allocated(err), more="forward CPCM nuclear gradient failed")
      if (allocated(error)) return
      call check(error, allocated(cavity%xyz1_rA), more="forward CPCM gradient did not build its Jacobian")
      if (allocated(error)) return
      call check(error, maxval(abs(forward_gradient - gradient)), 0.0_wp, thr=thr2, &
                 more="forward and reverse CPCM nuclear gradients differ")
      if (allocated(error)) return

      ! Isolate the raw host width contraction using a manufactured cavity
      ! derivative; the identical matrix and direct terms cancel between probes,
      ! leaving q_i dphi/dxi_i xi1_rA(:, :, i) at the probe point
      call cavity%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, "Width-gradient fixture setup failed: "//err%message)
         return
      end if
      width_point = maxloc(abs(pcm_model%q), dim=1)
      allocate (width_gradient(3, mols(1)%nat, 2), source=0.0_wp)
      allocate (width_forward(3, mols(1)%nat), source=0.0_wp)
      allocate (width_derivative(cavity%ngrid), source=0.0_wp)
      saved_xi1 = cavity%xi1_rA(:, :, width_point)
      cavity%xi1_rA(2, 1, width_point) = 0.37_wp
      width_expected = pcm_model%q(width_point)*0.23_wp*cavity%xi1_rA(:, :, width_point)
      call check(error, abs(width_expected(2, 1)) > 1.0e-8_wp, &
         & more="manufactured CPCM width gradient is zero")
      if (allocated(error)) return
      do iprobe = 1, 2
         call stage_point_charge_energy(error, pcm_model, cavity, qat, mols(1), coupling)
         if (allocated(error)) return
         call pcm_model%prepare_gradient(cavity, coupling, err)
         if (allocated(err)) then
            call test_failed(error, "Width-gradient staging failed: "//err%message)
            return
         end if
         width_derivative(width_point) = real(iprobe - 1, wp)*0.23_wp
         do while (coupling%next())
            select type (request => coupling%request())
            type is (gaussian_potential_request_type)
               call coupling%answer("dphi_dr", spread(spread(0.0_wp, 1, 3), 2, cavity%ngrid), err)
               if (.not. allocated(err)) call coupling%answer("dphi_dxi", width_derivative, err)
            type is (nuclear_charge_request_type)
               call coupling%answer("z", qat, err)
            end select
            if (allocated(err)) then
               call test_failed(error, "Width-gradient input failed: "//err%message)
               return
            end if
         end do
         call pcm_model%get_gradient(component_view(coupling), cavity, response, &
            & width_gradient(:, :, iprobe), err)
         if (allocated(err)) then
            call test_failed(error, "Width-gradient probe failed: "//err%message)
            return
         end if
         width_forward = 0.0_wp
         call pcm_model%get_forward_gradient(component_view(coupling), cavity, response, width_forward, err)
         call check(error,.not. allocated(err), more="forward manufactured width gradient failed")
         if (allocated(error)) return
         call check(error, maxval(abs(width_forward - width_gradient(:, :, iprobe))), 0.0_wp, thr=thr2, &
                    more="forward and reverse manufactured width gradients differ")
         if (allocated(error)) return
      end do
      ! Absolute 5e-12 with a 10x looser fallback relative to max|gradient|: the
      ! charges are bitwise identical between probes, but the adjoint contraction
      ! reduces over OpenMP threads in no fixed order, so a few ulp of |gradient|
      ! may differ; measured 5.6e-17 at max|gradient| = 10 (macOS arm64 gfortran
      ! 14.3, two and four threads)
      width_residual = maxval(abs(width_gradient(:, :, 2) - width_gradient(:, :, 1) &
         & - width_expected))
      width_thr = max(5.0e-12_wp, 5.0e-11_wp*maxval(abs(width_gradient)))
      call check(error, width_residual, 0.0_wp, thr=width_thr, &
         & more="CPCM nuclear raw width contraction")
      if (allocated(error)) return
      cavity%xi1_rA(:, :, width_point) = saved_xi1

      do iatom = 1, min(2, mols(1)%nat)
         do iaxis = 1, 3
            saved = mols(1)%xyz(iaxis, iatom)
            do k = 1, 4
               trial = mols(1)
               trial%xyz(iaxis, iatom) = saved + fd4_offsets(k)*step
               call displaced_energy(trial, values(k))
               if (allocated(error)) return
            end do
            call fd4_scalar(values(1), values(2), values(3), values(4), step, fd, error)
            if (allocated(error)) return
            write (context, "(a,i0,a,i0)") "CPCM gradient atom ", iatom, &
               & ", axis ", iaxis
            call check(error, gradient(iaxis, iatom), fd, &
               & thr_abs=fd_atol, thr_rel=fd_rtol, more=trim(context))
            if (allocated(error)) return
            call check(error, forward_gradient(iaxis, iatom), fd, &
               & thr_abs=fd_atol, thr_rel=fd_rtol, more="forward "//trim(context))
            if (allocated(error)) return
         end do
      end do

   contains

      !> Evaluate the CPCM energy for one displaced molecular structure
      !> The point-charge potential is re-evaluated on the displaced cavity
      !>
      !> @param[in]  displaced_mol Displaced molecular structure
      !> @param[out] energy        CPCM polarization energy
      subroutine displaced_energy(displaced_mol, energy)
         !> Displaced molecular structure
         type(structure_type), intent(in) :: displaced_mol
         !> CPCM polarization energy
         real(wp), intent(out) :: energy

         !> Local library error handling
         type(moist_error_type), allocatable :: local_err

         call get_test_cavity_iswig(displaced_mol, cavity%cavity_type_iswig, local_err, nleb=nleb)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced cavity setup failed: "//local_err%message)
            return
         end if
         call pcm_model%update(displaced_mol, cavity, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced CPCM update failed: "//local_err%message)
            return
         end if
         call pcm_model%prepare_energy(cavity, coupling, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced CPCM staging failed: "//local_err%message)
            return
         end if
         call fill_point_charge_potential(cavity, coupling, qat, displaced_mol)
         energy = 0.0_wp
         call pcm_model%get_energy(component_view(coupling), cavity, energy, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced CPCM energy failed: "//local_err%message)
         end if
      end subroutine displaced_energy

   end subroutine check_cpcm_nuclear_gradient

!> Check the external-potential CPCM nuclear gradient against finite differences
!>
!> The supplied potential is the sum of the solute nuclear potential and a
!> linear electronic potential. Both gradients are known in closed form, so
!> the host total position weight `w_xyz(:,i) = q_i*dphi/dr_i` (nuclear plus
!> electronic) is exact. This exercises the component's nuclear-charge
!> construction (the direct term) and the host position weight through the
!> complete component entry point, and pins that a mis-shaped weight answer
!> is refused before the accumulator is touched
   subroutine test_cpcm_external_nuclear_gradient(error)

      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Reference and displaced molecular structures
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: trial
      !> PCM component and molecular cavity
      type(model_continuum_component_cpcm) :: pcm_model
      type(cavity_type_iswig) :: cavity
      !> External coupling data
      type(coupling_type), target :: coupling
      !> External potential on the current cavity
      real(wp), allocatable :: phi(:)
      !> Host total position weight at the grid points
      real(wp), allocatable :: w_xyz(:, :)
      !> Analytic gradient and rejected-call sentinel
      real(wp), allocatable :: gradient(:, :), bad_gradient(:, :)
      !> Host part of the gradient phase, unused here
      type(response_type) :: response
      !> Four-point finite-difference energies
      real(wp) :: values(4)
      !> Finite-difference derivative and saved coordinate
      real(wp) :: fd, saved
      !> Surface, Cartesian-axis, and stencil indices
      integer :: igrid, iaxis, k
      !> Failure context
      character(len=96) :: context

      !> Dielectric constant
      real(wp), parameter :: epsilon = 78.4_wp
      !> Lebedev grid size
      integer, parameter :: nleb = 26
      !> Constant gradient of the synthetic electronic potential
      real(wp), parameter :: electronic_field(3) = [0.013_wp, -0.017_wp, 0.009_wp]
      !> Nuclear-displacement step in bohr
      real(wp), parameter :: step = 1.0e-4_wp
      !> Tight tolerances accounting for cavity and linear-solve roundoff
      real(wp), parameter :: fd_atol = 5.0e-9_wp
      real(wp), parameter :: fd_rtol = 5.0e-8_wp

      !> Run context owned here and borrowed by the component
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0)

      call get_test_structures(mols, 5)
      call center_at_origin(mols(1))
      call get_test_cavity_iswig(mols(1), cavity, err, nleb=nleb)
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
      call pcm_model%update(mols(1), cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if

      call pcm_model%new_coupling(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Coupling setup failed: "//err%message)
         return
      end if
      call pcm_model%prepare_energy(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Energy staging failed: "//err%message)
         return
      end if
      call build_external_potential(mols(1), cavity, phi)
      call set_host_potential(coupling, phi)
      call pcm_model%ensure_charges(component_view(coupling), cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM charge solution failed: "//err%message)
         return
      end if

      call pcm_model%prepare_gradient(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, "Gradient staging failed: "//err%message)
         return
      end if
      ! The nuclear charges of the external potential in a pass of their own
      do while (coupling%next())
         select type (request => coupling%request())
         type is (nuclear_charge_request_type)
            call coupling%answer("z", real(mols(1)%num(mols(1)%id), wp), err)
            if (allocated(err)) then
               call test_failed(error, "Answering the nuclear charges failed: "//err%message)
               return
            end if
         end select
      end do
      ! One step of the walk: the potential request stays current below, so
      ! the mis-shaped answer further down replaces this one
      if (.not. coupling%next()) then
         call test_failed(error, "CPCM asked for no potential derivatives in the gradient phase")
         return
      end if
      allocate (w_xyz(3, cavity%ngrid))
      do igrid = 1, cavity%ngrid
         w_xyz(:, igrid) = external_potential_gradient(mols(1), &
            & cavity%xyz(:, igrid)) + electronic_field
      end do
      select type (position => coupling%request())
      type is (gaussian_potential_request_type)
         call coupling%answer("dphi_dr", w_xyz, err)
         if (.not. allocated(err)) then
            call coupling%answer("dphi_dxi", spread(0.0_wp, 1, size(w_xyz, 2)), err)
         end if
      class default
         call test_failed(error, "CPCM asked for "//trim(position%name())//", not a Gaussian potential")
         return
      end select
      if (allocated(err)) then
         call test_failed(error, "Answering the position weights failed: "//err%message)
         return
      end if
      allocate (gradient(3, mols(1)%nat), source=0.0_wp)
      call pcm_model%get_gradient(component_view(coupling), cavity, response, gradient, err)
      if (allocated(err)) then
         call test_failed(error, "External CPCM nuclear gradient failed: "//err%message)
         return
      end if
      call check(error, maxval(abs(w_xyz)) > 1.0e-8_wp, &
         & more="external position-weight contribution is zero, the test is vacuous")
      if (allocated(error)) return

      ! A mis-shaped answer is rejected on submission and leaves the output
      ! unanswered, so the next accessor fails before it touches the
      ! accumulator; the next `prepare_energy` below stages afresh
      deallocate (w_xyz)
      allocate (w_xyz(2, cavity%ngrid), source=0.0_wp)
      call coupling%answer("dphi_dr", w_xyz, err)
      call check(error, allocated(err), more="mis-shaped dphi_dr was accepted")
      if (allocated(error)) return
      call coupling%answer("dphi_dxi", spread(0.0_wp, 1, size(w_xyz, 2)), err)
      allocate (bad_gradient(3, mols(1)%nat), source=1.0_wp)
      call pcm_model%get_gradient(component_view(coupling), cavity, response, bad_gradient, err)
      call check(error, allocated(err), &
         & more="external CPCM gradient accepted a malformed position weight")
      if (allocated(error)) return
      call check(error, maxval(abs(bad_gradient - 1.0_wp)), 0.0_wp, thr=0.0_wp, &
         & more="rejected external CPCM gradient modified the accumulator")
      if (allocated(error)) return
      deallocate (err)
      call pcm_model%get_forward_gradient(component_view(coupling), cavity, response, bad_gradient, err)
      call check(error, allocated(err), more="forward CPCM accepted a malformed position weight")
      if (allocated(error)) return
      call check(error, maxval(abs(bad_gradient - 1.0_wp)), 0.0_wp, thr=0.0_wp, &
                 more="rejected forward CPCM gradient modified the accumulator")
      if (allocated(error)) return
      deallocate (err)

      do iaxis = 1, 3
         saved = mols(1)%xyz(iaxis, 1)
         do k = 1, 4
            trial = mols(1)
            trial%xyz(iaxis, 1) = saved + fd4_offsets(k)*step
            call displaced_external_energy(trial, values(k))
            if (allocated(error)) return
         end do
         call fd4_scalar(values(1), values(2), values(3), values(4), step, fd, error)
         if (allocated(error)) return
         write (context, "(a,i0)") "external CPCM gradient atom 1, axis ", iaxis
         call check(error, gradient(iaxis, 1), fd, &
            & thr_abs=fd_atol, thr_rel=fd_rtol, more=trim(context))
         if (allocated(error)) return
      end do

   contains

      !> Build the synthetic nuclear-plus-electronic external potential
      !>
      !> @param[in]  potential_mol    Molecular structure providing nuclear sources
      !> @param[in]  potential_cavity Cavity carrying evaluation points
      !> @param[out] potential        External molecular potential
      subroutine build_external_potential(potential_mol, potential_cavity, potential)
         !> Molecular structure providing nuclear sources
         type(structure_type), intent(in) :: potential_mol
         !> Cavity carrying evaluation points
         type(cavity_type_iswig), intent(in) :: potential_cavity
         !> External molecular potential
         real(wp), allocatable, intent(out) :: potential(:)

         !> Surface and atom indices
         integer :: i, iatom
         !> Nuclear charge
         real(wp) :: za

         allocate (potential(potential_cavity%ngrid), source=0.0_wp)
         do i = 1, potential_cavity%ngrid
            do iatom = 1, potential_mol%nat
               za = real(potential_mol%num(potential_mol%id(iatom)), wp)
               potential(i) = potential(i) + za/norm2( &
                  & potential_cavity%xyz(:, i) - potential_mol%xyz(:, iatom))
            end do
            potential(i) = potential(i) + &
               & dot_product(electronic_field, potential_cavity%xyz(:, i))
         end do
      end subroutine build_external_potential

      !> Gradient of the nuclear part of the external potential at one point
      !>
      !> @param[in] potential_mol Molecular structure providing nuclear sources
      !> @param[in] point         Evaluation point
      !> @return grad_phi Gradient of the nuclear point-charge potential
      function external_potential_gradient(potential_mol, point) result(grad_phi)
         !> Molecular structure providing nuclear sources
         type(structure_type), intent(in) :: potential_mol
         !> Evaluation point
         real(wp), intent(in) :: point(3)
         !> Gradient of the nuclear point-charge potential
         real(wp) :: grad_phi(3)

         !> Atom index
         integer :: iatom
         !> Nuclear charge, separation vector and distance
         real(wp) :: za, r_vec(3), r_dist

         grad_phi = 0.0_wp
         do iatom = 1, potential_mol%nat
            za = real(potential_mol%num(potential_mol%id(iatom)), wp)
            r_vec = point - potential_mol%xyz(:, iatom)
            r_dist = norm2(r_vec)
            grad_phi = grad_phi - za*r_vec/(r_dist*r_dist*r_dist)
         end do
      end function external_potential_gradient

      !> Evaluate the external-potential CPCM energy for one displaced structure
      !>
      !> @param[in]  displaced_mol Displaced molecular structure
      !> @param[out] energy        CPCM polarization energy
      subroutine displaced_external_energy(displaced_mol, energy)
         !> Displaced molecular structure
         type(structure_type), intent(in) :: displaced_mol
         !> CPCM polarization energy
         real(wp), intent(out) :: energy

         !> Local library error handling
         type(moist_error_type), allocatable :: local_err
         !> External potential on the displaced cavity
         real(wp), allocatable :: displaced_phi(:)

         call get_test_cavity_iswig(displaced_mol, cavity, local_err, nleb=nleb)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced cavity setup failed: "//local_err%message)
            return
         end if
         call pcm_model%update(displaced_mol, cavity, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced CPCM update failed: "//local_err%message)
            return
         end if
         call pcm_model%prepare_energy(cavity, coupling, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced CPCM staging failed: "//local_err%message)
            return
         end if
         call build_external_potential(displaced_mol, cavity, displaced_phi)
         call set_host_potential(coupling, displaced_phi)
         energy = 0.0_wp
         call pcm_model%get_energy(component_view(coupling), cavity, energy, local_err)
         if (allocated(local_err)) then
            call test_failed(error, "Displaced CPCM energy failed: "//local_err%message)
         end if
      end subroutine displaced_external_energy

   end subroutine test_cpcm_external_nuclear_gradient

end module test_model_component_pcm_derivatives
