!> Shared fixtures of the PCM component tests, NOT a test suite
!>
!> - Component factories of all nine models: `new_pcm` (difference-check physics,
!>   dense solves) and `make_component` (constructor defaults, optional context)
!> - Solutes of the difference checks and the iSwiG/SvdW-DROP seed cavities
!> - `fd_cavity`, a synthetic cavity whose surface fields the test sets directly
!> - Host of point charges answering every request of a staged phase
module test_model_component_pcm_fixtures
   use mctc_env, only: wp, fatal_error
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: error_type, check, test_failed
   use test_helpers, only: check_moist_error, component_view, point_charge_potential, &
      & point_charge_field, submit, get_test_cross
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_type, only: cavity_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, &
      & moist_cavity_iswig_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_radii, only: default_cpcm_radii
   use moist_channels_coupling, only: coupling_type, gaussian_potential_request_type, &
      & point_potential_request_type, nuclear_charge_request_type
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

   public :: pcm_cpcm, pcm_cosmo, pcm_iefpcm, pcm_ssvpe, pcm_dpcm, pcm_dpcm_original, pcm_lpb, &
      & pcm_slic, pcm_nonlocal, npcm
   public :: pcm_physics, new_pcm, make_component
   public :: seed_iswig, seed_drop, seed_name, first_seed, new_seed_cavity
   public :: fd_cavity, get_seed
   public :: nmol, mol_triatomic, mol_label, get_solute
   public :: thr_vacuous
   public :: host_energy, stage_energy, answer_sources, host_charges

   !> Model identifiers of the component factory
   integer, parameter :: pcm_cpcm = 1, pcm_cosmo = 2, pcm_iefpcm = 3, pcm_ssvpe = 4, &
      & pcm_dpcm = 5, pcm_dpcm_original = 6, pcm_lpb = 7, pcm_slic = 8, pcm_nonlocal = 9
   !> Number of models
   integer, parameter :: npcm = 9

   !> Seeds of the synthetic cavity and cavities of the nuclear check
   integer, parameter :: seed_iswig = 1, seed_drop = 2
   !> Seed names for messages
   character(len=5), parameter :: seed_name(2) = ["iSwiG", "DROP "]

   !> Physical parameters of a component, defaults those of the difference checks
   type :: pcm_physics
      !> Static permittivity of every model, bulk permittivity of the nonlocal one
      real(wp) :: eps = 78.4_wp
      !> High-frequency permittivity of the nonlocal model
      real(wp) :: eps_inf = 1.8_wp
      !> Correlation length of the nonlocal model, bohr
      real(wp) :: lambda = 3.0_wp
      !> Relative permittivity inside the cavity of the nonlocal model
      real(wp) :: eps_solute = 1.0_wp
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
      !> SLIC static potential, hartree per elementary charge
      real(wp) :: slic_phi_static = 0.0_wp
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
   !> Magnitude below which an expected-active quantity counts as vacuous
   real(wp), parameter :: thr_vacuous = 1.0e-7_wp

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

   !> Synthetic cavity whose surface fields are set directly by the test
   type, extends(cavity_type) :: fd_cavity
   contains
      !> No-op: the fields are set by the test
      procedure :: update => fd_cavity_update
      !> No-op: the synthetic surface has no nuclear derivatives
      procedure :: get_gradient => fd_cavity_get_gradient
   end type fd_cavity

contains

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
               & err, param, mu=phys%slic_mu, phi_static=phys%slic_phi_static, ctx=ctx)
         type is (model_continuum_component_nonlocal)
            call new_component_nonlocal(pcm, eps, phys%eps_inf, phys%lambda, err, param, &
               & eps_solute=phys%eps_solute, ctx=ctx)
         end select
      end associate
      call check_moist_error(error, err, "component setup")
   end subroutine new_pcm

   !> Construct each PCM family with optional context forwarding
   !>
   !> @param[in] kind Constructor index
   !> @param[out] component Constructed component
   !> @param[out] err Library diagnostic
   !> @param[in] ctx Optional borrowed context
   subroutine make_component(kind, component, err, ctx)
      !> Constructor index
      integer, intent(in) :: kind
      !> Constructed component
      class(model_continuum_component_pcm), allocatable, intent(out) :: component
      !> Library diagnostic
      type(moist_error_type), allocatable, intent(out) :: err
      !> Optional borrowed context
      type(moist_context_type), target, optional, intent(in) :: ctx

      select case (kind)
      case (1)
         allocate (model_continuum_component_cpcm :: component)
         select type (component)
         type is (model_continuum_component_cpcm)
            call new_component_cpcm(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (2)
         allocate (model_continuum_component_cosmo :: component)
         select type (component)
         type is (model_continuum_component_cosmo)
            call new_component_cosmo(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (3)
         allocate (model_continuum_component_iefpcm :: component)
         select type (component)
         type is (model_continuum_component_iefpcm)
            call new_component_iefpcm(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (4)
         allocate (model_continuum_component_ssvpe :: component)
         select type (component)
         type is (model_continuum_component_ssvpe)
            call new_component_ssvpe(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (5)
         allocate (model_continuum_component_dpcm :: component)
         select type (component)
         type is (model_continuum_component_dpcm)
            call new_component_dpcm(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (6)
         allocate (model_continuum_component_dpcm_original :: component)
         select type (component)
         type is (model_continuum_component_dpcm_original)
            call new_component_dpcm_original(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (7)
         allocate (model_continuum_component_lpb :: component)
         select type (component)
         type is (model_continuum_component_lpb)
            call new_component_lpb(component, 80.0_wp, 0.1_wp, error=err, ctx=ctx)
         end select
      case (8)
         allocate (model_continuum_component_slic :: component)
         select type (component)
         type is (model_continuum_component_slic)
            call new_component_slic(component, 80.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, error=err, ctx=ctx)
         end select
      case (9)
         allocate (model_continuum_component_nonlocal :: component)
         select type (component)
         type is (model_continuum_component_nonlocal)
            call new_component_nonlocal(component, 80.0_wp, 2.0_wp, 1.0_wp, error=err, ctx=ctx)
         end select
      case default
         call fatal_error(err, "Unknown PCM component fixture")
      end select
   end subroutine make_component

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

   !> Update the component and stage a new coupling for the energy phase
   !>
   !> With `keep` the operators of the last update stay, for new host traces on
   !> an unchanged cavity; only the charge cache sees the change
   !>
   !> @param[in,out] pcm      Component under test
   !> @param[in]     mol      Solute
   !> @param[in,out] cavity   Cavity, current for `mol`
   !> @param[out]    coupling Coupling staged for the energy phase
   !> @param[out]    error    Test failure
   !> @param[in]     keep     Skip the update (default update)
   subroutine stage_energy(pcm, mol, cavity, coupling, error, keep)
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
      !> Skip the update
      logical, intent(in), optional :: keep

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Update skipped
      logical :: skip

      skip = .false.
      if (present(keep)) skip = keep
      if (.not. skip) then
         call pcm%update(mol, cavity, err)
         call check_moist_error(error, err, "update")
         if (allocated(error)) return
      end if
      call pcm%new_coupling(cavity, coupling, err)
      call check_moist_error(error, err, "coupling setup")
      if (allocated(error)) return
      call pcm%prepare_energy(cavity, coupling, err)
      call check_moist_error(error, err, "energy staging")
   end subroutine stage_energy

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

end module test_model_component_pcm_fixtures
