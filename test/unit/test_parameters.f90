!> Native parameter construction, JSON/TOML input/output, and copy ownership
module test_parameters
   use, intrinsic :: iso_c_binding, only: c_null_funptr, c_null_ptr
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_model_parameters, only: moist_model_parameters_type
   use moist, only: moist_cavity_drop_parameters_type, moist_cavity_iswig_parameters_type, &
                    moist_cavity_numsa_parameters_type, moist_cavity_mc_parameters_type, &
                    moist_cavity_drop_lsf_svdw_param_type, moist_cavity_drop_lsf_cfc_param_type, &
                    moist_cavity_drop_lsf_isodensity_param_type, moist_pcm_parameters_type, moist_gostshyp_parameters_type, &
                    moist_cavity_mt_parameters_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_branching, only: branch_weight_type
   use moist_cavity_drop_parameters, only: drop_seed_parameters_type, drop_branch_parameters_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_cavity_numsa, only: cavity_type_numsa, new_cavity_numsa
   use moist_cavity_mt, only: cavity_type_mt, new_cavity_mt
   use moist_cavity_mc, only: cavity_type_mc, new_cavity_mc
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_isodensity_callback, only: moist_cavity_drop_lsf_isodensity_callback_type
   use moist_model_continuum_component_type, only: model_continuum_component_type
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, new_component_cpcm
   use moist_model_continuum_component_gostshyp, only: model_continuum_component_gostshyp, new_component_gostshyp
   use moist_model_continuum_component_pcm_type, only: solver_type
   use moist_model_continuum_component_pcm_solvers, only: pcm_storage_dense
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cosmo, new_component_cosmo
   use moist_model_continuum_component_pcm_iefpcm, only: model_continuum_component_iefpcm, new_component_iefpcm
   use moist_model_continuum_component_pcm_lpb, only: model_continuum_component_lpb, new_component_lpb
   use moist_model_continuum_component_pcm_slic, only: model_continuum_component_slic, new_component_slic
   use moist_model_continuum_component_pcm_nonlocal, only: model_continuum_component_nonlocal, new_component_nonlocal
   use moist_context, only: moist_context_type, new_context
   use moist_radii, only: default_cpcm_radii
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use test_helpers, only: read_printout, printed_entry
   implicit none(type, external)
   private
   public :: collect_parameters

   !> Temporary path of the copy test
   character(len=*), parameter :: parameter_file = "moist-test-parameters-copies.json"

   !> Exercise less common registration helpers through the abstract interface
   type, extends(moist_model_parameters_type) :: vector_parameters_type
      real(wp) :: weights(2) = [1.0_wp, 2.0_wp]
      character(len=16) :: label = "default"
   contains
      procedure :: init_defaults => init_vector_defaults
      procedure :: register_entries => register_vector_entries
   end type vector_parameters_type

contains

   !> Check model-specific physical inputs before the first update
   subroutine test_pcm_physical_inputs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Construction error
      type(moist_error), allocatable :: err
      !> Components under test
      type(model_continuum_component_lpb) :: lpb
      type(model_continuum_component_slic) :: slic
      type(model_continuum_component_nonlocal) :: nonlocal
      type(model_continuum_component_cpcm) :: cpcm
      type(model_continuum_component_cosmo) :: cosmo
      type(model_continuum_component_iefpcm) :: iefpcm
      !> Captured parameter output
      character(len=160) :: lines(80)
      !> Scratch unit and line count
      integer :: unit, nline

      call new_component_lpb(lpb, 80.0_wp, 0.125_wp, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      open (newunit=unit, status="scratch", action="readwrite")
      call lpb%print_parameters(unit, level=2)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, printed_entry(lines, "Epsilon", "80.000000") .and. &
         & printed_entry(lines, "Kappa", "0.125000 bohr^-1"), "LPB screening input missing")
      if (allocated(error)) return
      call check(error,.not. printed_entry(lines, "f(eps)", ""), "LPB prints an internal factor")
      if (allocated(error)) return

      call new_component_slic(slic, 80.0_wp, 0.125_wp, 0.25_wp, 0.5_wp, err, &
         & mu=0.0625_wp, phi_static=1.25_wp)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      open (newunit=unit, status="scratch", action="readwrite")
      call slic%print_parameters(unit, level=2)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, printed_entry(lines, "Epsilon", "80.000000") .and. &
         & printed_entry(lines, "Alpha", "0.125000") .and. &
         & printed_entry(lines, "Beta", "0.250000 a.u.^-1") .and. &
         & printed_entry(lines, "Gamma", "0.500000") .and. &
         & printed_entry(lines, "Mu", "0.062500") .and. &
         & printed_entry(lines, "Static potential", "1.250000 Eh/e"), "SLIC physical inputs missing")
      if (allocated(error)) return
      call check(error,.not. printed_entry(lines, "f(eps)", ""), "SLIC prints an internal factor")
      if (allocated(error)) return

      call new_component_nonlocal(nonlocal, 80.0_wp, 2.0_wp, 10.0_wp, err, epsilon_solute=1.5_wp)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      open (newunit=unit, status="scratch", action="readwrite")
      call nonlocal%print_parameters(unit, level=2)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, printed_entry(lines, "Epsilon", "80.000000") .and. &
         & printed_entry(lines, "Epsilon (high-frequency)", "2.000000") .and. &
         & printed_entry(lines, "Epsilon (solute)", "1.500000") .and. &
         & printed_entry(lines, "Correlation length", "10.000000 bohr"), "Nonlocal physical inputs missing")
      if (allocated(error)) return
      call check(error,.not. printed_entry(lines, "f(eps)", ""), "Nonlocal prints an internal factor")
      if (allocated(error)) return

      call new_component_iefpcm(iefpcm, 80.0_wp, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      open (newunit=unit, status="scratch", action="readwrite")
      call iefpcm%print_parameters(unit, level=2)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, printed_entry(lines, "Epsilon", "80.000000") .and. &
         & .not. printed_entry(lines, "f(eps)", ""), "Layer PCM prints an internal factor")
      if (allocated(error)) return

      call new_component_cpcm(cpcm, 2.0_wp, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      open (newunit=unit, status="scratch", action="readwrite")
      call cpcm%print_parameters(unit, level=2)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, printed_entry(lines, "f(eps)", "0.500000"), "CPCM scaling factor missing")
      if (allocated(error)) return

      call new_component_cosmo(cosmo, 2.0_wp, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      open (newunit=unit, status="scratch", action="readwrite")
      call cosmo%print_parameters(unit, level=2)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, printed_entry(lines, "f(eps)", "0.400000"), "COSMO scaling factor missing")
   end subroutine test_pcm_physical_inputs

   !> Reset the helper parameter values
   !>
   !> @param[inout] self Test parameters
   subroutine init_vector_defaults(self)
      class(vector_parameters_type), intent(inout) :: self
      self%weights = [1.0_wp, 2.0_wp]
      self%label = "default"
   end subroutine init_vector_defaults

   !> Register vector and fixed-length string fields
   !>
   !> @param[inout] self Test parameters
   subroutine register_vector_entries(self)
      class(vector_parameters_type), intent(inout), target :: self
      call self%register_real_vector("weights", self%weights)
      call self%register_string("label", self%label, level=2)
   end subroutine register_vector_entries

   !> Print a parameter table into a scratch unit and read its lines back
   !>
   !> @param[in] param Parameter values
   !> @param[out] lines Printed lines
   !> @param[out] nline Number of printed lines
   !> @param[in] level Verbosity level of the printout
   subroutine capture_table(param, lines, nline, level)
      class(moist_model_parameters_type), intent(in) :: param
      character(len=*), intent(out) :: lines(:)
      integer, intent(out) :: nline
      integer, intent(in), optional :: level
      type(prettyprinter) :: pp
      integer :: unit

      open (newunit=unit, status="scratch", action="readwrite")
      pp = new_prettyprinter(unit=unit)
      call param%print_table(pp, level)
      call read_printout(unit, lines, nline)
      close (unit)
   end subroutine capture_table

   !> A table printed for a verbosity keeps the fields at or below it
   subroutine test_print_levels(error)
      type(error_type), allocatable, intent(out) :: error
      type(vector_parameters_type) :: vector
      type(moist_gostshyp_parameters_type) :: gostshyp
      type(moist_pcm_parameters_type) :: pcm
      character(len=128) :: lines(64)
      integer :: nline, unit
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(cavity_type_iswig) :: iswig
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(cavity_type_drop) :: drop

      call capture_table(vector, lines, nline)
      call check(error, printed_entry(lines, "label", "default") .and. printed_entry(lines, "weights(1)", "1."), &
         & "a table without a level lists every field")
      if (allocated(error)) return

      call capture_table(vector, lines, nline, level=2)
      call check(error, printed_entry(lines, "label", "default"), "level 2 keeps the normal field")
      if (allocated(error)) return
      call check(error,.not. printed_entry(lines, "weights(1)", "1."), "level 2 drops the detailed field")
      if (allocated(error)) return

      ! Every GOSTSHYP setting is detailed: the normal table is empty
      call capture_table(gostshyp, lines, nline, level=2)
      call check(error, nline, 0, "a table without normal fields prints nothing")
      if (allocated(error)) return
      call capture_table(gostshyp, lines, nline, level=3)
      call check(error, printed_entry(lines, "regularization_start", ""), "level 3 lists the detailed fields")
      if (allocated(error)) return

      call check(error, pcm%solver, solver_type%iterative, "default iterative solver")
      if (allocated(error)) return
      call check(error, pcm%storage, pcm_storage_dense, "default dense storage")
      if (allocated(error)) return
      call capture_table(pcm, lines, nline, level=2)
      call check(error, nline, 2, "normal PCM settings contain only the named solver and storage")
      if (allocated(error)) return
      call check(error, printed_entry(lines, "Solver", "Iterative") .and. printed_entry(lines, "Storage", "Dense"), &
         & "level 2 names the solver and storage")
      if (allocated(error)) return
      call check(error, .not. (printed_entry(lines, "solver", "") .or. printed_entry(lines, "storage", "")), &
         & "level 2 drops the raw solver and storage")
      if (allocated(error)) return
      call capture_table(pcm, lines, nline, level=3)
      call check(error, nline, 11, "detailed PCM settings include the names and all solver controls")
      if (allocated(error)) return
      call check(error, printed_entry(lines, "solver", "") .and. printed_entry(lines, "storage", "") .and. &
         & printed_entry(lines, "solver_tol", "") .and. &
         & printed_entry(lines, "solver_maxiter", "") .and. printed_entry(lines, "solver_restart", "") .and. &
         & printed_entry(lines, "newton_tol", "") .and. printed_entry(lines, "newton_maxiter", "") .and. &
         & printed_entry(lines, "point_potential", ""))
      if (allocated(error)) return

      ! Storage only applies to the iterative solver
      pcm%solver = solver_type%dense
      call capture_table(pcm, lines, nline, level=2)
      call check(error, nline, 1, "the dense solver omits the storage")
      if (allocated(error)) return
      call check(error, printed_entry(lines, "Solver", "Dense"), "level 2 names the dense solver")
      if (allocated(error)) return

      ! Nothing prints at verbosity 0
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_cavity_iswig(iswig, default_cpcm_radii(), err, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      open (newunit=unit, status="scratch", action="readwrite")
      call iswig%print_parameters(unit)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, nline, 0, "verbosity 0 prints settings")
      if (allocated(error)) return
      call ctx%delete()

      ! From verbosity 1 an explicit request prints the normal settings
      call new_context(ctx, nthreads=0, verbosity=1)
      call new_cavity_iswig(iswig, default_cpcm_radii(), err, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      open (newunit=unit, status="scratch", action="readwrite")
      call iswig%print_parameters(unit)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, printed_entry(lines, "num_leb", "110"), "verbosity 1 omits the normal settings")
      if (allocated(error)) return
      call check(error, .not. printed_entry(lines, "cut_a", ""), "verbosity 1 lists the detailed settings")
      if (allocated(error)) return

      call lsf%new()
      call new_cavity_drop(drop, default_cpcm_radii(), lsf, err, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      open (newunit=unit, status="scratch", action="readwrite")
      call drop%print_parameters(unit)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, printed_entry(lines, "Number of Leb. points", "194"), &
         & "verbosity 1 omits the normal DROP settings")
      if (allocated(error)) return
      call check(error, .not. printed_entry(lines, "Projection tol.", ""), "verbosity 1 lists the detailed DROP settings")
      call ctx%delete()
   end subroutine test_print_levels

   !> Level 3 adds the detailed settings of every cavity, level set and component
   !>
   !> Each detailed entry is also required to be absent at level 2, so the
   !> check fails if the level gate stops working in either direction
   subroutine test_print_detailed(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Run context borrowed by every part
      type(moist_context_type), target :: ctx
      !> DROP settings reaching the octree and switching rows
      type(moist_cavity_drop_parameters_type) :: param
      !> Level sets
      type(moist_cavity_drop_lsf_svdw_type) :: svdw
      type(moist_cavity_drop_lsf_cfc_type) :: cfc
      type(moist_cavity_drop_lsf_isodensity_callback_type) :: rho
      !> Cavities
      type(cavity_type_drop) :: drop
      type(cavity_type_iswig) :: iswig
      type(cavity_type_numsa) :: numsa
      type(cavity_type_mt) :: mt
      type(cavity_type_mc) :: mc
      !> Components
      type(model_continuum_component_cpcm) :: pcm
      type(model_continuum_component_gostshyp) :: gostshyp

      call new_context(ctx, nthreads=0, verbosity=0)
      call param%new(proj_level=4, wleb_prune_level=2, error=err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return

      call svdw%new()
      call new_cavity_drop(drop, default_cpcm_radii(), svdw, err, param, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_cavity(drop, [character(len=32) :: "Projection tol.", "Wleb switch off", "xi_born", &
         "Octree seed size", "Max. distance excess", "Cell grid full-scan below", "Surface normals"], error)
      if (allocated(error)) return

      call cfc%new()
      call new_cavity_drop(drop, default_cpcm_radii(), cfc, err, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_cavity(drop, [character(len=32) :: "Projection tol."], error, &
         normal=[character(len=32) :: "Atomic exponent (a1)"])
      if (allocated(error)) return

      call rho%new(c_null_funptr, c_null_ptr, moist_cavity_drop_lsf_isodensity_param_type(log_grad_out=6.0_wp))
      call new_cavity_drop(drop, default_cpcm_radii(), rho, err, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_cavity(drop, [character(len=32) :: "ECP core fill radius (of R_A)", &
         "Exclusion |grad ln rho| (out)", "Exclusion radius cap"], error, &
         normal=[character(len=32) :: "Isovalue (rho_iso)", "ECP core fill"])
      if (allocated(error)) return

      call new_cavity_iswig(iswig, default_cpcm_radii(), err, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_cavity(iswig, [character(len=32) :: "cut_a", "cut_f"], error)
      if (allocated(error)) return

      call new_cavity_numsa(numsa, default_cpcm_radii(), err, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_cavity(numsa, [character(len=32) :: "offset", "smoothing", "tolsesp"], error, &
         normal=[character(len=32) :: "probe"])
      if (allocated(error)) return

      call new_cavity_mt(mt, default_cpcm_radii(), svdw, err, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_cavity(mt, [character(len=32) :: "area_tolerance", "w_0ls_from", "w_0ls_to", "w_0ls_p", &
         "w_0ls_a", "weight_cutoff"], error, normal=[character(len=32) :: "spacing", "target_area", "gaussian_scale"])
      if (allocated(error)) return

      ! Unset export paths print nothing
      call new_cavity_mc(mc, default_cpcm_radii(), svdw, err, &
         moist_cavity_mc_parameters_type(obj_file="mesh.obj", pqr_file="mesh.pqr"), ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_cavity(mc, [character(len=32) :: "obj_file", "pqr_file"], error, &
         normal=[character(len=32) :: "spacing"])
      if (allocated(error)) return

      call new_component_cpcm(pcm, 80.0_wp, error=err, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_component(pcm, [character(len=32) :: "solver", "storage", "solver_tol", "solver_maxiter"], error, &
         normal=[character(len=32) :: "Solver", "Storage"])
      if (allocated(error)) return

      call new_component_gostshyp(gostshyp, 1.0e-6_wp, ctx=ctx, error=err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_component(gostshyp, [character(len=32) :: "regularization_start", "regularization_end", &
         "suppress_negative_amplitudes"], error, normal=[character(len=32) :: "Pressure"])
      call ctx%delete()

   contains

      !> Print a cavity at levels 2 and 3 and compare the detailed entries
      !>
      !> @param[in] cavity Constructed cavity
      !> @param[in] detailed Entries printed only from level 3
      !> @param[out] error Test failure
      !> @param[in] normal Entries printed at both levels
      subroutine require_cavity(cavity, detailed, error, normal)
         !> Constructed cavity
         class(cavity_type), intent(in) :: cavity
         !> Entries printed only from level 3
         character(len=*), intent(in) :: detailed(:)
         !> Test failure
         type(error_type), allocatable, intent(out) :: error
         !> Entries printed at both levels
         character(len=*), intent(in), optional :: normal(:)
         !> Printed lines per level
         character(len=128) :: lines2(128), lines3(128)
         !> Line count and scratch unit
         integer :: nline, unit

         open (newunit=unit, status="scratch", action="readwrite")
         call cavity%print_parameters(unit, level=2)
         call read_printout(unit, lines2, nline)
         close (unit)
         open (newunit=unit, status="scratch", action="readwrite")
         call cavity%print_parameters(unit, level=3)
         call read_printout(unit, lines3, nline)
         close (unit)
         call compare_levels(lines2, lines3, detailed, error, normal)
      end subroutine require_cavity

      !> Print a component at levels 2 and 3 and compare the detailed entries
      !>
      !> @param[in] component Constructed component
      !> @param[in] detailed Entries printed only from level 3
      !> @param[out] error Test failure
      !> @param[in] normal Entries printed at both levels
      subroutine require_component(component, detailed, error, normal)
         !> Constructed component
         class(model_continuum_component_type), intent(in) :: component
         !> Entries printed only from level 3
         character(len=*), intent(in) :: detailed(:)
         !> Test failure
         type(error_type), allocatable, intent(out) :: error
         !> Entries printed at both levels
         character(len=*), intent(in), optional :: normal(:)
         !> Printed lines per level
         character(len=128) :: lines2(64), lines3(64)
         !> Line count and scratch unit
         integer :: nline, unit

         open (newunit=unit, status="scratch", action="readwrite")
         call component%print_parameters(unit, level=2)
         call read_printout(unit, lines2, nline)
         close (unit)
         open (newunit=unit, status="scratch", action="readwrite")
         call component%print_parameters(unit, level=3)
         call read_printout(unit, lines3, nline)
         close (unit)
         call compare_levels(lines2, lines3, detailed, error, normal)
      end subroutine require_component

      !> Detailed entries appear only at level 3, normal entries at both
      !>
      !> @param[in] lines2 Level-2 printout
      !> @param[in] lines3 Level-3 printout
      !> @param[in] detailed Entries printed only from level 3
      !> @param[out] error Test failure
      !> @param[in] normal Entries printed at both levels
      subroutine compare_levels(lines2, lines3, detailed, error, normal)
         !> Level-2 printout
         character(len=*), intent(in) :: lines2(:)
         !> Level-3 printout
         character(len=*), intent(in) :: lines3(:)
         !> Entries printed only from level 3
         character(len=*), intent(in) :: detailed(:)
         !> Test failure
         type(error_type), allocatable, intent(out) :: error
         !> Entries printed at both levels
         character(len=*), intent(in), optional :: normal(:)
         !> Entry index
         integer :: i

         do i = 1, size(detailed)
            call check(error, printed_entry(lines3, trim(detailed(i)), ""), &
               more="level 3 omits "//trim(detailed(i)))
            if (allocated(error)) return
            call check(error, .not. printed_entry(lines2, trim(detailed(i)), ""), &
               more="level 2 prints "//trim(detailed(i)))
            if (allocated(error)) return
         end do
         if (.not. present(normal)) return
         do i = 1, size(normal)
            call check(error, printed_entry(lines2, trim(normal(i)), "") .and. &
               printed_entry(lines3, trim(normal(i)), ""), more="a level omits "//trim(normal(i)))
            if (allocated(error)) return
         end do
      end subroutine compare_levels

   end subroutine test_print_detailed

   !> DROP settings load from TOML and JSON with their derived values recomputed
   !>
   !> References for the derived values come from the in-memory constructor,
   !> not from restated formulas
   subroutine test_drop_load_file(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Loaded settings and the constructor reference
      type(moist_cavity_drop_parameters_type) :: loaded, reference
      !> Owned I/O state
      integer :: unit, stat, close_stat
      !> Owned files
      character(len=*), parameter :: toml_path = "moist-test-parameters-drop-load.toml"
      character(len=*), parameter :: json_path = "moist-test-parameters-drop-load.json"

      call reference%new(nleb=110, tolerance=3.0e-8_wp, proj_level=4, wleb_prune_level=3, error=err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return

      open (newunit=unit, file=toml_path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)", iostat=stat) "tolerance = 3e-8", "[grid]", "num_leb = 110", &
         "[projection]", "level = 4", "[switching]", "wleb_prune_level = 3"
      close (unit, iostat=close_stat)
      call check(error, stat == 0 .and. close_stat == 0)
      if (allocated(error)) return
      call loaded%load_file(toml_path, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call same_settings(loaded, reference, error)
      if (allocated(error)) return

      ! A written file loads back to the same settings
      call loaded%write_file(json_path, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call loaded%init_defaults()
      call loaded%load_file(json_path, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call same_settings(loaded, reference, error)
      if (allocated(error)) return

      open (newunit=unit, file=toml_path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      open (newunit=unit, file=json_path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete", iostat=stat)
      call check(error, stat, 0)

   contains

      !> Compare the read and the derived settings with the reference
      !>
      !> @param[in] loaded Settings read from a file
      !> @param[in] reference Settings from the constructor
      !> @param[out] error Test failure
      subroutine same_settings(loaded, reference, error)
         !> Settings read from a file
         type(moist_cavity_drop_parameters_type), intent(in) :: loaded
         !> Settings from the constructor
         type(moist_cavity_drop_parameters_type), intent(in) :: reference
         !> Test failure
         type(error_type), allocatable, intent(out) :: error

         call check(error, loaded%num_leb == reference%num_leb .and. loaded%proj_level == reference%proj_level &
            .and. loaded%wleb_prune_level == reference%wleb_prune_level)
         if (allocated(error)) return
         call check(error, maxval(abs([loaded%tolerance, loaded%solve%tol, loaded%iswig_xi_born, &
            loaded%finalize%wleb_prune_off, loaded%finalize%wleb_prune_on, loaded%finalize%wleb_cut, loaded%branch%dpdist_max] - &
            [reference%tolerance, reference%solve%tol, reference%iswig_xi_born, reference%finalize%wleb_prune_off, &
            reference%finalize%wleb_prune_on, reference%finalize%wleb_cut, reference%branch%dpdist_max])), 0.0_wp, thr=1.0e-14_wp)
      end subroutine same_settings

   end subroutine test_drop_load_file

   !> Check vector values and escaped fixed-length strings survive serialization
   subroutine test_vector_and_string(error)
      type(error_type), allocatable, intent(out) :: error
      type(vector_parameters_type) :: param
      !> Library error and owned fixture I/O
      type(moist_error), allocatable :: err
      !> File unit and I/O statuses
      integer :: unit, stat, close_stat, format_index
      !> Independent vector input
      character(len=*), parameter :: path = "moist-test-parameters-vector-input.json"

      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)", iostat=stat) '{"weights":[0.3,-5.0],"label":"input"}'
      close (unit, iostat=close_stat)
      call check(error, stat == 0 .and. close_stat == 0)
      if (allocated(error)) return
      call param%read_file(path, err)
      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, maxval(abs(param%weights - [0.3_wp, -5.0_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, trim(param%label) == "input")
      if (allocated(error)) return
      param%weights = [0.3_wp, -5.0_wp]
      param%label = 'a "quoted" name'
      do format_index = 1, 2
         call roundtrip(param, "moist-test-parameters-vector", error, format_index)
         if (allocated(error)) return
         call check(error, maxval(abs(param%weights - [0.3_wp, -5.0_wp])), 0.0_wp, thr=1.0e-14_wp)
         if (allocated(error)) return
         call check(error, trim(param%label) == 'a "quoted" name')
         if (allocated(error)) return
      end do
   end subroutine test_vector_and_string

   !> Structurally wrong fields are refused by name
   !>
   !> A group key holding a scalar, a vector of the wrong length and an
   !> oversized fixed-length string each fail in their own registration helper
   subroutine test_invalid_fields(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Parameter set with dotted group keys
      type(moist_cavity_drop_parameters_type) :: drop
      !> Fresh vector/string sets, one per refused field
      type(vector_parameters_type) :: vector, label
      !> Temporary path owned by this test
      character(len=*), parameter :: path = "moist-test-parameters-invalid.json"

      call read_invalid(drop, '{"grid": 110}', "Invalid parameter group: grid.", error)
      if (allocated(error)) return
      call read_invalid(vector, '{"weights": [1.0, 2.0, 3.0]}', &
                        "Invalid parameter vector: weights", error)
      if (allocated(error)) return
      call read_invalid(label, '{"label": "seventeen letters!"}', &
                        "Parameter string is too long: label", error)

   contains

      !> Write one document, read it back and require the named refusal
      !>
      !> @param[inout] param Parameter set reading the document
      !> @param[in] text JSON document
      !> @param[in] expected Distinctive substring of the expected message
      !> @param[out] error Test failure
      subroutine read_invalid(param, text, expected, error)
         !> Parameter set reading the document
         class(moist_model_parameters_type), intent(inout) :: param
         !> JSON document
         character(len=*), intent(in) :: text
         !> Distinctive substring of the expected message
         character(len=*), intent(in) :: expected
         !> Test failure
         type(error_type), allocatable, intent(out) :: error
         !> Library error handling
         type(moist_error), allocatable :: err
         !> File unit and I/O status
         integer :: unit, stat

         open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write (unit, "(a)", iostat=stat) text
         close (unit)
         call check(error, stat, 0)
         if (allocated(error)) return
         call param%read_file(path, err)
         open (newunit=unit, file=path, status="old", iostat=stat)
         if (stat == 0) close (unit, status="delete", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call check(error, allocated(err), more="accepted "//text)
         if (allocated(error)) return
         call check(error, index(err%message, expected) > 0, more=err%message)
      end subroutine read_invalid

   end subroutine test_invalid_fields

   !> Fail an expected-failure test only on the targeted library error
   !>
   !> - Used by `should_fail=.true.` tests, where test-drive inverts the
   !>   verdict: a raised test failure passes, a clean return fails
   !> - Fails the test only when `err` names `expected`; no error or a
   !>   different one returns cleanly, so the test is reported as failed
   !>
   !> @param[out] error     test failure, set only on the expected error
   !> @param[in]  err       library error, possibly unallocated
   !> @param[in]  expected  substring the expected error message must contain
   subroutine expect_error(error, err, expected)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error, possibly unallocated
      type(moist_error), allocatable, intent(in) :: err
      !> Substring the expected error message must contain
      character(len=*), intent(in) :: expected

      if (.not. allocated(err)) return
      if (index(err%message, expected) > 0) call test_failed(error, err%message)
   end subroutine expect_error

   !> A grid without a fitted Born zeta is refused at load
   subroutine test_bad_drop_load_unfitted_grid(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      type(moist_cavity_drop_parameters_type) :: loaded
      integer :: unit, stat, close_stat
      character(len=*), parameter :: toml_path = "moist-test-parameters-drop-unfitted.toml"

      open (newunit=unit, file=toml_path, status="replace", action="write", iostat=stat)
      ! Setup failures return cleanly, so the test is reported as failed
      if (stat /= 0) return
      write (unit, "(a)", iostat=stat) "[grid]", "num_leb = 1"
      close (unit, iostat=close_stat)
      if (stat /= 0 .or. close_stat /= 0) return
      call loaded%load_file(toml_path, err)
      open (newunit=unit, file=toml_path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete", iostat=stat)
      call expect_error(error, err, "No fitted DROP Born zeta")
   end subroutine test_bad_drop_load_unfitted_grid

   !> A missing parameter file is an error, not a silent default
   subroutine test_bad_drop_load_missing_file(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      type(moist_cavity_drop_parameters_type) :: loaded
      character(len=*), parameter :: toml_path = "moist-test-parameters-drop-missing.toml"

      call loaded%load_file(toml_path, err)
      call expect_error(error, err, "Cannot read parameter file")
   end subroutine test_bad_drop_load_missing_file

   !> A CFC pair power other than 4 is refused when read from JSON
   subroutine test_bad_cfc_pair_power(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      !> Fixture file unit and I/O status
      integer :: unit, stat
      !> Fixture path
      character(len=*), parameter :: path = "moist-test-parameters-cfc-pair-power.json"

      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      if (stat /= 0) return
      write (unit, "(a)", iostat=stat) '{"m":6}'
      close (unit, iostat=stat)
      call param%read_file(path, err)
      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete", iostat=stat)
      call expect_error(error, err, "CFC pair power m must be 4")
   end subroutine test_bad_cfc_pair_power

   !> A non-positive projection batch size is refused
   subroutine test_bad_drop_batch_size(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%batch_size = 0
      call param%compute_derived(err)
      call expect_error(error, err, "Projection batch size must be positive")
   end subroutine test_bad_drop_batch_size

   !> A negative Newton stage-switch residual is refused
   subroutine test_bad_drop_newton_switch_fnorm(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%solve%switch_fnorm = -1.0e-3_wp
      call param%compute_derived(err)
      call expect_error(error, err, "Newton stage-switch residual must be nonnegative")
   end subroutine test_bad_drop_newton_switch_fnorm


   !> A horizon band whose on-edge lies below its off-edge is refused
   subroutine test_bad_drop_horizon_bounds(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%w_0hor_off = 0.2_wp
      param%w_0hor_on = 0.1_wp
      call param%compute_derived(err)
      call expect_error(error, err, "Horizon switch bounds must satisfy 0 <= off < on <= 1")
   end subroutine test_bad_drop_horizon_bounds

   !> A zero on-edge of the horizon band is accepted and leaves no cosine inside
   !>
   !> - The default band is handed on unchanged
   subroutine test_drop_horizon_off(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      call param%compute_derived(err)
      call check(error,.not. allocated(err), "default horizon band is refused")
      if (allocated(error)) return
      call check(error, param%hor_off == param%w_0hor_off .and. param%hor_on == param%w_0hor_on &
                 .and. param%hor_on > 0.0_wp, "default horizon band is not handed on")
      if (allocated(error)) return

      param%w_0hor_on = 0.0_wp
      call param%compute_derived(err)
      call check(error,.not. allocated(err), "zero horizon on-edge is refused")
      if (allocated(error)) return
      call check(error, param%hor_off < param%hor_on .and. param%hor_on <= -1.0_wp, &
                 "disabled horizon band still holds cosines")
   end subroutine test_drop_horizon_off

   !> Reversed critical-gradient edges are refused: they would open the gate at
   !> critical points and close it on the regular surface
   subroutine test_bad_drop_critical_bounds(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%w_0ls_off = 0.6_wp
      param%w_0ls_on = 0.25_wp
      call param%compute_derived(err)
      call expect_error(error, err, "Critical-gradient switch bounds must satisfy 0 <= off < on")
   end subroutine test_bad_drop_critical_bounds

   !> A negative focal off-edge is refused: it would open the gate past the focal point
   subroutine test_bad_drop_focal_bounds(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%w_0tra_off = -0.1_wp
      call param%compute_derived(err)
      call expect_error(error, err, "Focal switch bounds must satisfy 0 <= off < on")
   end subroutine test_bad_drop_focal_bounds

   !> A zero sigmoid exponent is refused: the switch would lose its flat ends
   subroutine test_bad_drop_switch_shape(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%w_0ls_p = 0.0_wp
      call param%compute_derived(err)
      call expect_error(error, err, "Switch shape parameters w_0ls_p and w_0ls_a must be positive")
   end subroutine test_bad_drop_switch_shape

   !> A key no field reads is refused, here the name before the pdist rename
   subroutine test_bad_unknown_key(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param
      !> Whether the fixture was written and removed
      logical :: ok
      !> Fixture path
      character(len=*), parameter :: path = "moist-test-parameters-unknown-key.json"

      call read_fixture(param, path, [character(len=64) :: &
                        '{"grid":{"num_leb":110},"branching":{"drho_max":0.5}}'], err, ok)
      if (.not. ok) return
      call expect_error(error, err, "Unknown parameter: branching.drho_max")
   end subroutine test_bad_unknown_key

   !> Every unknown key is named, nested ones with their group; known keys
   !> next to them are no complaint
   subroutine test_unknown_keys_named(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param
      !> Whether the fixture was written and removed
      logical :: ok
      !> Fixture path
      character(len=*), parameter :: path = "moist-test-parameters-unknown-keys.toml"

      call read_fixture(param, path, [character(len=32) :: "do_rho = true", "[grid]", "num_leb = 110", &
                        "[projection.newton]", "max_cheap_steps = 3", "cheap_steps = 3"], err, ok)
      call check(error, ok, "fixture file")
      if (allocated(error)) return
      call check(error, allocated(err), "unknown keys are accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Unknown parameters: do_rho, projection.newton.cheap_steps") > 0, &
                 more=err%message)
      if (allocated(error)) return
      ! Control: the same file without the unknown keys is read
      call read_fixture(param, path, [character(len=32) :: "[grid]", "num_leb = 110", "[projection.newton]", &
                        "max_cheap_steps = 3"], err, ok)
      call check(error, ok, "fixture file")
      if (allocated(error)) return
      call check(error, .not. allocated(err), "known keys are refused")
      if (allocated(error)) return
      call check(error, param%num_leb == 110 .and. param%solve%max_cheap_steps == 3, "known keys are not read")
   end subroutine test_unknown_keys_named

   !> Write a parameter file, read it and remove it again
   !>
   !> @param[inout] param Settings that read the file
   !> @param[in]    path  Fixture path, owned by the calling test
   !> @param[in]    lines File lines, trailing blanks dropped
   !> @param[out]   err   Library error of the read
   !> @param[out]   ok    Whether the file was written and removed
   subroutine read_fixture(param, path, lines, err, ok)
      class(moist_model_parameters_type), intent(inout) :: param
      character(len=*), intent(in) :: path
      character(len=*), intent(in) :: lines(:)
      type(moist_error), allocatable, intent(out) :: err
      logical, intent(out) :: ok
      integer :: unit, stat, close_stat, i

      ok = .false.
      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      if (stat /= 0) return
      do i = 1, size(lines)
         write (unit, "(a)", iostat=stat) trim(lines(i))
         if (stat /= 0) exit
      end do
      close (unit, iostat=close_stat)
      if (stat /= 0 .or. close_stat /= 0) return
      call param%read_file(path, err)
      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat /= 0) return
      close (unit, status="delete", iostat=stat)
      ok = stat == 0
   end subroutine read_fixture

   !> A reach band whose open edge lies beyond its closed edge is refused
   subroutine test_bad_drop_reach_bounds(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%branch%reach_on = 3.0_wp
      param%branch%reach_off = 2.5_wp
      call param%compute_derived(err)
      call expect_error(error, err, "Reach bounds must satisfy 0 < on < off")
   end subroutine test_bad_drop_reach_bounds

   !> A reach band of zero width is refused: equal edges would close the gate at the anchor
   subroutine test_bad_drop_reach_zero_width(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%branch%reach_on = 2.5_wp
      param%branch%reach_off = 2.5_wp
      call param%compute_derived(err)
      call expect_error(error, err, "Reach bounds must satisfy 0 < on < off")
   end subroutine test_bad_drop_reach_zero_width

   !> The default reach is an accepted band, an explicit one is kept as it is,
   !> and the defaults restore the band and switch the full-reach search off
   subroutine test_drop_reach_defaults(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param
      !> Default band
      real(wp) :: reach_on, reach_off

      call param%init_defaults()
      call param%compute_derived(err)
      call check(error, .not. allocated(err) .and. param%branch%reach_on > 0.0_wp &
                 .and. param%branch%reach_on < param%branch%reach_off, "the default reach is no band")
      if (allocated(error)) return
      reach_on = param%branch%reach_on
      reach_off = param%branch%reach_off
      param%branch%reach_on = 1.5_wp
      param%branch%reach_off = 2.5_wp
      call param%compute_derived(err)
      call check(error, .not. allocated(err) .and. param%branch%reach_on == 1.5_wp &
                 .and. param%branch%reach_off == 2.5_wp, "an explicit reach is not kept")
      if (allocated(error)) return
      call check(error, .not. param%branch%full_reach, "full-reach search is on by default")
      if (allocated(error)) return
      param%branch%full_reach = .true.
      call param%init_defaults()
      call check(error, .not. param%branch%full_reach, "defaults leave the full-reach search on")
      if (allocated(error)) return
      call check(error, param%branch%reach_on == reach_on .and. param%branch%reach_off == reach_off, &
                 "defaults do not restore the reach band")
   end subroutine test_drop_reach_defaults

   !> The octree preset follows the projection level: the coarse octree at
   !> levels 1 and 3 drops a leaf without sign change, the fine one of level 4
   !> has smaller leaves and refines them
   subroutine test_drop_seed_follows_level(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param
      !> Preset of each level
      type(drop_seed_parameters_type) :: seed1, seed3, seed4

      call param%init_defaults()
      param%proj_level = 1
      call param%compute_derived(err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      seed1 = param%seed
      param%proj_level = 4
      call param%compute_derived(err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      seed4 = param%seed
      param%proj_level = 3
      call param%compute_derived(err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      seed3 = param%seed
      call check(error, seed1%seed_size == seed3%seed_size .and. seed1%refine_levels == seed3%refine_levels, &
                 "the level-1 fallback is not the coarse octree")
      if (allocated(error)) return
      call check(error, seed3%refine_levels == 0 .and. seed4%refine_levels > 0, &
                 "only the fine octree refines a leaf without sign change")
      if (allocated(error)) return
      call check(error, seed4%seed_size < seed3%seed_size, "the fine octree has no smaller leaves")
   end subroutine test_drop_seed_follows_level

   !> A nonpositive octree box budget is refused
   subroutine test_bad_drop_octree_max_boxes(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%seed%max_boxes = 0
      call param%compute_derived(err)
      call expect_error(error, err, "Octree box budget must be positive")
   end subroutine test_bad_drop_octree_max_boxes

   !> A nonpositive octree survivor budget is refused
   subroutine test_bad_drop_octree_max_survivors(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%seed%max_survivors = -1
      call param%compute_derived(err)
      call expect_error(error, err, "Octree survivor budget must be positive")
   end subroutine test_bad_drop_octree_max_survivors

   !> An octree depth limit past the addressable 20 levels is refused
   subroutine test_bad_drop_octree_max_depth(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%seed%max_depth = 21
      call param%compute_derived(err)
      call expect_error(error, err, "Octree depth limit must lie in 1..20")
   end subroutine test_bad_drop_octree_max_depth

   !> A non-finite SvdW fade band is refused while the fade is on
   subroutine test_bad_svdw_fade_nonfinite(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_lsf_svdw_param_type) :: param

      param%fade_enabled = .true.
      param%fade_width = ieee_value(1.0_wp, ieee_positive_inf)
      call param%validate(err)
      call expect_error(error, err, "SvdW fade_start and fade_width must be finite")
   end subroutine test_bad_svdw_fade_nonfinite

   !> A zero SvdW fade onset is refused while the fade is on
   subroutine test_bad_svdw_fade_nonpositive(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_lsf_svdw_param_type) :: param

      param%fade_enabled = .true.
      param%fade_start = 0.0_wp
      call param%validate(err)
      call expect_error(error, err, "SvdW fade_start and fade_width must be positive")
   end subroutine test_bad_svdw_fade_nonpositive

   !> A non-finite CFC fade band is refused while the fade is on
   subroutine test_bad_cfc_fade_nonfinite(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_lsf_cfc_param_type) :: param

      param%fade_enabled = .true.
      param%fade_width = ieee_value(1.0_wp, ieee_positive_inf)
      call param%validate(err)
      call expect_error(error, err, "CFC fade_start and fade_width must be finite")
   end subroutine test_bad_cfc_fade_nonfinite

   !> A zero CFC fade onset is refused while the fade is on
   subroutine test_bad_cfc_fade_nonpositive(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_lsf_cfc_param_type) :: param

      param%fade_enabled = .true.
      param%fade_start = 0.0_wp
      call param%validate(err)
      call expect_error(error, err, "CFC fade_start and fade_width must be positive")
   end subroutine test_bad_cfc_fade_nonpositive

   !> A CFC fade with `a2 = 0` is refused: the band counts e-folds of `exp(a2 sigma)`
   subroutine test_bad_cfc_fade_a2(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_lsf_cfc_param_type) :: param

      param%fade_enabled = .true.
      param%a2 = 0.0_wp
      call param%validate(err)
      call expect_error(error, err, "CFC fade needs a2 < 0")
   end subroutine test_bad_cfc_fade_a2

   !> A negative Newton cheap-step budget is refused
   subroutine test_bad_drop_newton_max_cheap_steps(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%solve%max_cheap_steps = -1
      call param%compute_derived(err)
      call expect_error(error, err, "Newton cheap-step budget must not be negative")
   end subroutine test_bad_drop_newton_max_cheap_steps

   !> A nonpositive Newton step cap is refused
   subroutine test_bad_drop_newton_max_step(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%solve%max_step = 0.0_wp
      call param%compute_derived(err)
      call expect_error(error, err, "Newton step cap must be positive")
   end subroutine test_bad_drop_newton_max_step

   !> A nonpositive Newton merit margin is refused
   subroutine test_bad_drop_newton_merit_margin(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      param%solve%merit_margin = -1.0_wp
      call param%compute_derived(err)
      call expect_error(error, err, "Newton merit margin must be positive")
   end subroutine test_bad_drop_newton_merit_margin

   !> A NaN tolerance is refused; a `<= 0` test let it through
   subroutine test_bad_drop_tolerance_nan(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      ! Control: a finite positive value is accepted; a refusal returns cleanly,
      ! so the test is reported as failed
      param%tolerance = 1.0e-8_wp
      call param%compute_derived(err)
      if (allocated(err)) return
      param%tolerance = ieee_value(1.0_wp, ieee_quiet_nan)
      call param%compute_derived(err)
      call expect_error(error, err, "tolerance must be finite and positive")
   end subroutine test_bad_drop_tolerance_nan

   !> An infinite tolerance is refused
   subroutine test_bad_drop_tolerance_inf(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      ! Control, as in the NaN case
      param%tolerance = 1.0e-8_wp
      call param%compute_derived(err)
      if (allocated(err)) return
      param%tolerance = ieee_value(1.0_wp, ieee_positive_inf)
      call param%compute_derived(err)
      call expect_error(error, err, "tolerance must be finite and positive")
   end subroutine test_bad_drop_tolerance_inf

   !> A NaN density kernel length is refused; a `<= 0` test let it through
   subroutine test_bad_drop_rho_grid_h_nan(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      ! Control: a finite positive value is accepted; a refusal returns cleanly,
      ! so the test is reported as failed
      param%rho_grid_h = 0.5_wp
      call param%compute_derived(err)
      if (allocated(err)) return
      param%rho_grid_h = ieee_value(1.0_wp, ieee_quiet_nan)
      call param%compute_derived(err)
      call expect_error(error, err, "rho_grid_h must be finite and positive")
   end subroutine test_bad_drop_rho_grid_h_nan

   !> An infinite density kernel length is refused
   subroutine test_bad_drop_rho_grid_h_inf(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      ! Control, as in the NaN case
      param%rho_grid_h = 0.5_wp
      call param%compute_derived(err)
      if (allocated(err)) return
      param%rho_grid_h = ieee_value(1.0_wp, ieee_positive_inf)
      call param%compute_derived(err)
      call expect_error(error, err, "rho_grid_h must be finite and positive")
   end subroutine test_bad_drop_rho_grid_h_inf

   !> A NaN objective alpha is refused; a `<= 0` test let it through
   subroutine test_bad_drop_phi_alpha_nan(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      ! Control: a finite positive value is accepted; a refusal returns cleanly,
      ! so the test is reported as failed
      param%phi_alpha = 0.25_wp
      call param%compute_derived(err)
      if (allocated(err)) return
      param%phi_alpha = ieee_value(1.0_wp, ieee_quiet_nan)
      call param%compute_derived(err)
      call expect_error(error, err, "Objective alpha must be finite and positive")
   end subroutine test_bad_drop_phi_alpha_nan

   !> An infinite objective alpha is refused
   subroutine test_bad_drop_phi_alpha_inf(error)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param

      call param%init_defaults()
      ! Control, as in the NaN case
      param%phi_alpha = 0.25_wp
      call param%compute_derived(err)
      if (allocated(err)) return
      param%phi_alpha = ieee_value(1.0_wp, ieee_positive_inf)
      call param%compute_derived(err)
      call expect_error(error, err, "Objective alpha must be finite and positive")
   end subroutine test_bad_drop_phi_alpha_inf

   !> Collect parameter API regressions
   subroutine collect_parameters(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("constructors", test_constructors), &
         new_unittest("file_roundtrip", test_roundtrip), &
         new_unittest("copies_and_errors", test_copies_and_errors), &
         new_unittest("vector_and_string", test_vector_and_string), &
         new_unittest("print_levels", test_print_levels), &
         new_unittest("print_detailed", test_print_detailed), &
         new_unittest("drop_load_file", test_drop_load_file), &
         new_unittest("bad_drop_load_unfitted_grid", test_bad_drop_load_unfitted_grid, should_fail=.true.), &
         new_unittest("bad_drop_load_missing_file", test_bad_drop_load_missing_file, should_fail=.true.), &
         new_unittest("bad_cfc_pair_power", test_bad_cfc_pair_power, should_fail=.true.), &
         new_unittest("bad_drop_batch_size", test_bad_drop_batch_size, should_fail=.true.), &
         new_unittest("bad_drop_newton_switch_fnorm", test_bad_drop_newton_switch_fnorm, should_fail=.true.), &
         new_unittest("bad_drop_newton_max_cheap_steps", test_bad_drop_newton_max_cheap_steps, should_fail=.true.), &
         new_unittest("bad_drop_horizon_bounds", test_bad_drop_horizon_bounds, should_fail=.true.), &
         new_unittest("drop_horizon_off", test_drop_horizon_off), &
         new_unittest("bad_drop_critical_bounds", test_bad_drop_critical_bounds, should_fail=.true.), &
         new_unittest("bad_drop_focal_bounds", test_bad_drop_focal_bounds, should_fail=.true.), &
         new_unittest("bad_drop_switch_shape", test_bad_drop_switch_shape, should_fail=.true.), &
         new_unittest("bad_unknown_key", test_bad_unknown_key, should_fail=.true.), &
         new_unittest("unknown_keys_named", test_unknown_keys_named), &
         new_unittest("bad_drop_reach_bounds", test_bad_drop_reach_bounds, should_fail=.true.), &
         new_unittest("bad_drop_reach_zero_width", test_bad_drop_reach_zero_width, should_fail=.true.), &
         new_unittest("drop_reach_defaults", test_drop_reach_defaults), &
         new_unittest("drop_seed_follows_level", test_drop_seed_follows_level), &
         new_unittest("bad_drop_octree_max_boxes", test_bad_drop_octree_max_boxes, should_fail=.true.), &
         new_unittest("bad_drop_octree_max_survivors", test_bad_drop_octree_max_survivors, should_fail=.true.), &
         new_unittest("bad_drop_octree_max_depth", test_bad_drop_octree_max_depth, should_fail=.true.), &
         new_unittest("bad_svdw_fade_nonfinite", test_bad_svdw_fade_nonfinite, should_fail=.true.), &
         new_unittest("bad_svdw_fade_nonpositive", test_bad_svdw_fade_nonpositive, should_fail=.true.), &
         new_unittest("bad_cfc_fade_nonfinite", test_bad_cfc_fade_nonfinite, should_fail=.true.), &
         new_unittest("bad_cfc_fade_nonpositive", test_bad_cfc_fade_nonpositive, should_fail=.true.), &
         new_unittest("bad_cfc_fade_a2", test_bad_cfc_fade_a2, should_fail=.true.), &
         new_unittest("bad_drop_newton_max_step", test_bad_drop_newton_max_step, should_fail=.true.), &
         new_unittest("bad_drop_newton_merit_margin", test_bad_drop_newton_merit_margin, should_fail=.true.), &
         new_unittest("bad_drop_tolerance_nan", test_bad_drop_tolerance_nan, should_fail=.true.), &
         new_unittest("bad_drop_tolerance_inf", test_bad_drop_tolerance_inf, should_fail=.true.), &
         new_unittest("bad_drop_rho_grid_h_nan", test_bad_drop_rho_grid_h_nan, should_fail=.true.), &
         new_unittest("bad_drop_rho_grid_h_inf", test_bad_drop_rho_grid_h_inf, should_fail=.true.), &
         new_unittest("bad_drop_phi_alpha_nan", test_bad_drop_phi_alpha_nan, should_fail=.true.), &
         new_unittest("bad_drop_phi_alpha_inf", test_bad_drop_phi_alpha_inf, should_fail=.true.), &
         new_unittest("pcm_physical_inputs", test_pcm_physical_inputs), &
         new_unittest("invalid_fields", test_invalid_fields), &
         new_unittest("file_formats", test_file_formats), &
         new_unittest("drop_contracts", test_drop_contracts), &
         new_unittest("configuration_fields", test_configuration_fields), &
         new_unittest("parameter_validation", test_parameter_validation), &
         new_unittest("optional_overrides", test_optional_overrides)]
   end subroutine collect_parameters

   !> Check copied settings, derived values, and resetting to defaults
   subroutine test_constructors(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moist_cavity_drop_parameters_type) :: param
      type(moist_cavity_drop_lsf_svdw_param_type) :: shape
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_isodensity_callback_type) :: rho
      type(cavity_type_drop) :: drop
      type(cavity_type_iswig) :: iswig
      type(model_continuum_component_cpcm) :: pcm

      call new_context(ctx, nthreads=0, verbosity=0)
      shape%blend_k = 7.0_wp
      call lsf%new(shape)
      shape%blend_k = 9.0_wp
      call check(error, lsf%param%blend_k, 7.0_wp)
      if (allocated(error)) return
      param%num_leb = 110
      param%tolerance = 1.0e-8_wp
      param%solve%tol = -1.0_wp
      param%do_fine = .true.
      call new_cavity_drop(drop, default_cpcm_radii(), lsf, err, param, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, drop%param%solve%tol, param%tolerance)
      if (allocated(error)) return
      call check(error, param%solve%tol, -1.0_wp)
      if (allocated(error)) return
      call check(error, drop%param%do_grid_density .and. drop%param%do_curvature .and. &
                 drop%param%do_normal .and. drop%param%do_atomarea, &
                 more="do_fine did not switch every optional property on")
      if (allocated(error)) return
      call check(error,.not. param%do_curvature, more="constructor changed the caller's settings")
      if (allocated(error)) return
      call check(error, drop%lsf_model%screening_threshold, 1.0e-9_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call new_cavity_drop(drop, default_cpcm_radii(), lsf, err, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, drop%param%num_leb, 194)
      if (allocated(error)) return
      call check(error,.not. (drop%param%do_grid_density .or. drop%param%do_curvature .or. &
                              drop%param%do_normal .or. drop%param%do_atomarea), &
                 more="default DROP settings switch an optional property on")
      if (allocated(error)) return
      call new_cavity_iswig(iswig, default_cpcm_radii(), err, &
                            moist_cavity_iswig_parameters_type(num_leb=194, cut_a=0.3_wp, cut_f=0.2_wp), ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, iswig%param%num_leb, 194)
      if (allocated(error)) return
      call check(error, iswig%param%cut_a, 0.3_wp)
      if (allocated(error)) return
      call check(error, iswig%param%cut_f, 0.2_wp)
      if (allocated(error)) return
      call new_cavity_iswig(iswig, default_cpcm_radii(), err, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, iswig%param%num_leb, 110)
      if (allocated(error)) return
      call check(error, iswig%param%cut_f, 1.0e-10_wp)
      if (allocated(error)) return
      call lsf%new()
      call check(error, lsf%param%blend_k, 5.5_wp)
      if (allocated(error)) return
      call rho%new(c_null_funptr, c_null_ptr, &
                   moist_cavity_drop_lsf_isodensity_param_type(rho_iso=0.004_wp, exclusion_cap=4.0_wp))
      call check(error, rho%param%exclusion_cap, 4.0_wp)
      if (allocated(error)) return
      call new_component_cpcm(pcm, 80.0_wp, error=err, &
         & param=moist_pcm_parameters_type(solver=solver_type%iterative, &
         & solver_tol=1.0e-8_wp, solver_maxiter=71), ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, pcm%param%solver, solver_type%iterative)
      if (allocated(error)) return
      call check(error, pcm%param%solver_tol, 1.0e-8_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call check(error, pcm%param%solver_maxiter, 71)
      if (allocated(error)) return
      call new_component_cpcm(pcm, 80.0_wp, error=err, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, pcm%param%solver, solver_type%iterative)
      if (allocated(error)) return
      param%num_leb = 1
      call new_cavity_drop(drop, default_cpcm_radii(), lsf, err, param, ctx=ctx)
      call check(error, allocated(err))
      call ctx%delete()
   end subroutine test_constructors

   !> All concrete parameter sets support the same file and print operations
   subroutine test_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_parameters_type) :: drop
      type(moist_cavity_iswig_parameters_type) :: iswig
      type(moist_cavity_numsa_parameters_type) :: numsa
      type(moist_cavity_mc_parameters_type) :: mc
      type(moist_cavity_drop_lsf_svdw_param_type) :: svdw
      type(moist_cavity_drop_lsf_cfc_param_type) :: cfc
      type(moist_cavity_drop_lsf_isodensity_param_type) :: rho
      type(moist_pcm_parameters_type) :: pcm
      !> Registered hydrostatic-pressure settings
      type(moist_gostshyp_parameters_type) :: gostshyp
      character(len=*), parameter :: stem = "moist-test-parameters-roundtrip"
      !> Each input format must preserve logical state after a single pass
      integer :: format_index

      drop%tolerance = 1.0e-8_wp
      drop%do_fine = .true.
      drop%branch%dpdist_max = 0.3_wp
      do format_index = 1, 2
         call roundtrip(drop, stem, error, format_index)
         if (allocated(error)) return
         call check(error, drop%solve%tol, 1.0e-8_wp)
         if (allocated(error)) return
         call check(error, drop%do_fine)
         if (allocated(error)) return
         call check(error, drop%branch%dpdist_max, 0.3_wp)
         if (allocated(error)) return
      end do
      iswig%cut_f = 0.02_wp
      call roundtrip(iswig, stem, error)
      if (allocated(error)) return
      call check(error, iswig%cut_f, 0.02_wp)
      if (allocated(error)) return
      numsa%smoothing = 0.4_wp
      call roundtrip(numsa, stem, error)
      if (allocated(error)) return
      call check(error, numsa%smoothing, 0.4_wp)
      if (allocated(error)) return
      mc%obj_file = 'path with "quotes" and \slash/mesh.obj'
      mc%spacing = 0.4_wp
      call roundtrip(mc, stem, error)
      if (allocated(error)) return
      call check(error, allocated(mc%obj_file))
      if (allocated(error)) return
      call check(error, mc%obj_file == 'path with "quotes" and \slash/mesh.obj')
      if (allocated(error)) return
      call check(error,.not. allocated(mc%pqr_file))
      if (allocated(error)) return
      svdw%blend_k = 6.0_wp
      call roundtrip(svdw, stem, error)
      if (allocated(error)) return
      call check(error, svdw%blend_k, 6.0_wp)
      if (allocated(error)) return
      cfc%a1 = -12.0_wp
      call roundtrip(cfc, stem, error)
      if (allocated(error)) return
      call check(error, cfc%a1, -12.0_wp)
      if (allocated(error)) return
      rho%rho_iso = 0.005_wp
      rho%exclusion_cap = 3.0_wp
      call roundtrip(rho, stem, error)
      if (allocated(error)) return
      call check(error, rho%exclusion_cap, 3.0_wp)
      if (allocated(error)) return
      pcm%solver = solver_type%dense
      pcm%solver_maxiter = 42
      call roundtrip(pcm, stem, error)
      if (allocated(error)) return
      call check(error, pcm%solver_maxiter, 42)
      if (allocated(error)) return
      gostshyp%regularization_start = 2.0e-9_wp
      gostshyp%regularization_end = 2.0e-8_wp
      gostshyp%suppress_negative_amplitudes = .true.
      do format_index = 1, 2
         call roundtrip(gostshyp, stem, error, format_index)
         if (allocated(error)) return
         call check(error, gostshyp%regularization_start, 2.0e-9_wp, thr=1.0e-21_wp)
         if (allocated(error)) return
         call check(error, gostshyp%regularization_end, 2.0e-8_wp, thr=1.0e-20_wp)
         if (allocated(error)) return
         call check(error, gostshyp%suppress_negative_amplitudes)
         if (allocated(error)) return
      end do
      call gostshyp%init_defaults()
      call check(error, gostshyp%regularization_start, 1.0e-12_wp, thr=1.0e-24_wp)
      if (allocated(error)) return
      call check(error, gostshyp%regularization_end, 1.0e-10_wp, thr=1.0e-22_wp)
      if (allocated(error)) return
      call check(error,.not. gostshyp%suppress_negative_amplitudes)
   end subroutine test_roundtrip

   !> Exercise the abstract file/printing interface without knowing the concrete type
   !>
   !> @param[inout] param Parameter values
   !> @param[in] stem File name without extension, owned by the calling test
   !> @param[out] error Test failure
   !> @param[in] format_index Select JSON or TOML; omission exercises both
   subroutine roundtrip(param, stem, error, format_index)
      class(moist_model_parameters_type), intent(inout) :: param
      character(len=*), intent(in) :: stem
      type(error_type), allocatable, intent(out) :: error
      !> Select one format for assertions between passes
      integer, intent(in), optional :: format_index
      type(moist_error), allocatable :: err
      integer :: unit, stat, i, first, last
      character(len=5), parameter :: extensions(2) = [".json", ".toml"]

      first = 1
      last = size(extensions)
      if (present(format_index)) then
         first = format_index
         last = format_index
      end if
      do i = first, last
         call param%write_file(stem//extensions(i), err)
         call check(error,.not. allocated(err))
         if (allocated(error)) return
         call param%init_defaults()
         call param%read_file(stem//extensions(i), err)
         call check(error,.not. allocated(err))
         if (allocated(error)) return
         open (newunit=unit, status="scratch", action="readwrite", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call param%print_parameters(err, unit)
         close (unit, iostat=stat)
         call check(error,.not. allocated(err))
         if (allocated(error)) return
         call check(error, stat, 0)
         if (allocated(error)) return
         open (newunit=unit, file=stem//extensions(i), status="old", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         close (unit, status="delete", iostat=stat)
         call check(error, stat, 0)
      end do
   end subroutine roundtrip

   !> Read handwritten TOML, reject malformed input and unsupported extensions
   subroutine test_file_formats(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_cavity_drop_parameters_type) :: param
      type(moist_cavity_iswig_parameters_type) :: iswig
      integer :: unit, stat, close_stat, i
      character(len=40), parameter :: bad_paths(3) = [character(len=40) :: &
         & "moist-test-parameters-formats.txt", "moist-test-parameters-formats", &
         & "directory.json/no-extension"]
      logical :: exists
      !> Uppercase extension; must not match another test's file on case-insensitive systems
      character(len=*), parameter :: toml_file = "moist-test-parameters-formats.TOML"

      open (newunit=unit, file=toml_file, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)", iostat=stat) "# Independently authored TOML", &
         "tolerance = 2e-8", "[grid]", "num_leb = 110"
      close (unit, iostat=close_stat)
      call check(error, stat == 0 .and. close_stat == 0)
      if (allocated(error)) return
      call param%read_file(toml_file, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, param%num_leb, 110)
      if (allocated(error)) return
      call check(error, param%solve%tol, 2.0e-8_wp)
      if (allocated(error)) return
      call check(error, param%solve%max_iter, 50)
      if (allocated(error)) return
      ! Uppercase output extensions select TOML too
      call param%write_file(toml_file, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call param%read_file(toml_file, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      open (newunit=unit, file=toml_file, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)", iostat=stat) "num_leb = ["
      close (unit, iostat=close_stat)
      call check(error, stat == 0 .and. close_stat == 0)
      if (allocated(error)) return
      call iswig%read_file(toml_file, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      open (newunit=unit, file=toml_file, status="old", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      close (unit, status="delete", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      call iswig%read_file(toml_file, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      do i = 1, size(bad_paths)
         call iswig%write_file(trim(bad_paths(i)), err)
         call check(error, allocated(err))
         if (allocated(error)) return
         inquire (file=trim(bad_paths(i)), exist=exists)
         call check(error,.not. exists)
         if (allocated(error)) return
         call iswig%read_file(trim(bad_paths(i)), err)
         call check(error, allocated(err))
         if (allocated(error)) return
      end do
   end subroutine test_file_formats

   !> Copies never share field bindings; malformed input returns an error
   subroutine test_copies_and_errors(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_cavity_drop_parameters_type) :: source, copied
      type(moist_cavity_iswig_parameters_type) :: iswig
      integer :: unit, stat
      integer :: i
      character(len=32), parameter :: invalid(3) = [character(len=32) :: &
                                                    "[]", "{broken", '{"num_leb": "bad"}']

      source%tolerance = 1.0e-8_wp
      call source%write_file(parameter_file, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      copied = source
      copied%tolerance = 1.0e-6_wp
      call copied%write_file(parameter_file, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call copied%read_file(parameter_file, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, copied%tolerance, 1.0e-6_wp)
      if (allocated(error)) return
      call check(error, source%tolerance, 1.0e-8_wp)
      if (allocated(error)) return
      do i = 1, size(invalid)
         open (newunit=unit, file=parameter_file, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write (unit, "(a)", iostat=stat) trim(invalid(i))
         close (unit)
         call check(error, stat, 0)
         if (allocated(error)) return
         call iswig%read_file(parameter_file, err)
         call check(error, allocated(err))
         if (allocated(error)) return
      end do
      open (newunit=unit, file=parameter_file, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)", iostat=stat) '{"num_leb": 194}'
      close (unit)
      call check(error, stat, 0)
      if (allocated(error)) return
      iswig%cut_f = 0.5_wp
      call iswig%read_file(parameter_file, err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, iswig%num_leb, 194)
      if (allocated(error)) return
      call check(error, iswig%cut_f, 1.0e-10_wp)
      if (allocated(error)) return
      call iswig%write_file("missing-moist-parameter-directory/config.json", err)
      call check(error, allocated(err))
      if (allocated(error)) return
      open (newunit=unit, file=parameter_file, status="old", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      close (unit, status="delete", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      call iswig%read_file(parameter_file, err)
      call check(error, allocated(err))
   end subroutine test_copies_and_errors

   !> Exercise DROP derivation from physical bounds and invalid boundary inputs
   subroutine test_drop_contracts(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Settings under test
      type(moist_cavity_drop_parameters_type) :: param
      !> Pruning and invalid-input cases
      integer :: level, invalid
      !> Independent pruning decade
      real(wp) :: lower

      call param%new(nleb=110, tolerance=2.0e-7_wp, proj_maxiter=81, proj_level=4, &
                     branch=drop_branch_parameters_type(dpdist_max=0.03_wp), rho_grid_h=0.75_wp, wleb_prune_level=0, error=err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, param%num_leb == 110 .and. param%solve%max_iter == 81 .and. param%proj_level == 4)
      if (allocated(error)) return
      call check(error, param%finalize%wleb_cut, 1.0e-8_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call check(error, param%screening_threshold, 2.0e-8_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call check(error, param%solve%tol, 2.0e-7_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call check(error, param%branch%sep_cut, 2.0e-6_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      ! Adjacency cutoff is four kernel lengths; the soft density screens with
      ! cutoff*rho_grid_h as 4 h^2
      call check(error, param%adj_list_grid_cutoff, 4.0_wp*0.75_wp)
      if (allocated(error)) return
      call check(error, param%iswig_xi_born > 4.8_wp .and. param%iswig_xi_born < 5.0_wp)
      if (allocated(error)) return
      ! The branch window is an input and survives the derivation untouched
      call check(error, param%branch%dpdist_max, 0.03_wp, thr=0.0_wp)
      if (allocated(error)) return
      param%phi_alpha = 3.0_wp
      call param%compute_derived(err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      ! The admissible radius is where the branch weight vanishes, whatever alpha
      block
         type(branch_weight_type) :: bw
         real(wp) :: pdist_max, w(2)

         call bw%init(param%branch%dpdist_max, param%phi_alpha)
         pdist_max = bw%pdist_max(2.0_wp)
         call bw%weights(0.5_wp*param%phi_alpha*[2.0_wp, pdist_max*(1.0_wp + 1.0e-12_wp)]**2, &
                         [1.0_wp, 1.0_wp], w)
         call check(error, w(2) == 0.0_wp, "a branch at the admissible radius has no weight")
         if (allocated(error)) return
         call bw%weights(0.5_wp*param%phi_alpha*[2.0_wp, 0.5_wp*(2.0_wp + pdist_max)]**2, &
                         [1.0_wp, 1.0_wp], w)
         call check(error, w(2) > 0.0_wp, "a branch inside the admissible radius has weight")
         if (allocated(error)) return
         call check(error, pdist_max - 2.0_wp, param%branch%dpdist_max, thr=1.0e-2_wp*param%branch%dpdist_max)
         if (allocated(error)) return
      end block
      do level = 0, 6
         call param%new(wleb_prune_level=level, error=err)
         call check(error,.not. allocated(err))
         if (allocated(error)) return
         if (level == 0) then
            call check(error, param%finalize%wleb_prune_off == 0.0_wp .and. param%finalize%wleb_prune_on == 0.0_wp)
         else
            lower = 10.0_wp**(-14 + 2*level)
            call check(error, param%finalize%wleb_prune_off/lower, 1.0_wp, thr=1.0e-14_wp)
            if (allocated(error)) return
            call check(error, param%finalize%wleb_prune_on/lower, 100.0_wp, thr=1.0e-12_wp)
            if (allocated(error)) return
            if (allocated(error)) return
            ! The cut does not follow the prune onset
            call check(error, param%finalize%wleb_cut, 0.05_wp*param%tolerance, thr=1.0e-20_wp)
         end if
         if (allocated(error)) return
      end do
      do invalid = 1, 11
         call param%init_defaults()
         select case (invalid)
         case (1)
            param%tolerance = 0.0_wp
         case (2)
            param%rho_grid_h = 0.0_wp
         case (3)
            param%solve%max_iter = 0
         case (4)
            param%proj_level = 0
         case (5)
            param%proj_level = 5
         case (6)
            param%proj_level = 2
         case (7)
            param%branch%dpdist_max = ieee_value(1.0_wp, ieee_quiet_nan)
         case (8)
            param%phi_alpha = 0.0_wp
         case (9)
            param%branch%dpdist_max = -0.1_wp
         case (10)
            param%wleb_prune_level = 7
         case (11)
            param%num_leb = 1
         case default
            error stop "invalid DROP test case"
         end select
         call param%compute_derived(err)
         call check(error, allocated(err), more="accepted invalid DROP setting")
         if (allocated(error)) return
      end do
   end subroutine test_drop_contracts

   !> Independently authored configuration inputs exercise public field names
   subroutine test_configuration_fields(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Concrete parameter sets
      type(moist_cavity_drop_parameters_type) :: drop
      type(moist_cavity_iswig_parameters_type) :: iswig
      type(moist_cavity_numsa_parameters_type) :: numsa
      type(moist_cavity_mc_parameters_type) :: mc
      type(moist_cavity_mt_parameters_type) :: mt
      type(moist_cavity_drop_lsf_svdw_param_type) :: svdw
      type(moist_cavity_drop_lsf_cfc_param_type) :: cfc
      type(moist_cavity_drop_lsf_isodensity_param_type) :: rho
      type(moist_pcm_parameters_type) :: pcm

      call configured(drop, '{"tolerance":2e-7,"do_fine":true,"grid":{"num_leb":110},"objective":{"alpha":2},'// &
                      '"projection":{"maxiter":72,"level":4,"batch_size":64,'// &
                      '"newton":{"switch_fnorm":1e-2,"max_cheap_steps":2,'// &
                      '"max_step":0.6,"merit_margin":2.5},'// &
                      '"octree":{"max_boxes":123,"max_survivors":43,'// &
                      '"max_depth":8}},'// &
                      '"screening":{"cell_grid_full_scan_below":19},'// &
                      '"switching":{"w_0ls_off":0.2,"w_0ls_on":0.7,"w_0ls_p":0.9,"w_0ls_a":1.7,'// &
                      '"w_0tra_off":0.15,"w_0tra_on":0.4,"w_0hor_off":0.02,"w_0hor_on":0.3,'// &
                      '"wleb_prune_level":2},'// &
                      '"density":{"rho_grid_h":0.8},"branching":{"dpdist_max":0.04,"reach_on":1.5,'// &
                      '"reach_off":2.25}}', error)
      if (allocated(error)) return
      call check(error, drop%num_leb == 110 .and. drop%do_fine .and. &
                 drop%solve%max_iter == 72 .and. drop%proj_level == 4 .and. &
                 drop%seed%max_boxes == 123 .and. drop%seed%max_survivors == 43 .and. &
                 drop%seed%max_depth == 8 .and. &
                 drop%cell_grid_full_scan_below == 19 .and. drop%wleb_prune_level == 2)
      if (allocated(error)) return
      call check(error, drop%batch_size == 64 .and. drop%solve%max_cheap_steps == 2, &
                 more="batched projection integer keys")
      if (allocated(error)) return
      call check(error, maxval(abs([drop%solve%switch_fnorm, drop%solve%max_step, drop%solve%merit_margin] - &
                                   [1.0e-2_wp, 0.6_wp, 2.5_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, maxval(abs([drop%tolerance, drop%phi_alpha, drop%w_0ls_off, &
                                    drop%w_0ls_on, drop%w_0ls_p, drop%w_0ls_a, &
                                    drop%w_0tra_off, drop%w_0tra_on, drop%w_0hor_off, drop%w_0hor_on, &
                                    drop%rho_grid_h, drop%branch%dpdist_max, &
                                    drop%branch%reach_on, drop%branch%reach_off] - &
                                   [2.0e-7_wp, 2.0_wp, 0.2_wp, 0.7_wp, 0.9_wp, 1.7_wp, &
                                    0.15_wp, 0.4_wp, 0.02_wp, 0.3_wp, &
                                    0.8_wp, 0.04_wp, 1.5_wp, 2.25_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, drop%do_grid_density .and. drop%do_curvature .and. drop%do_normal .and. &
                 drop%do_atomarea, &
                 more="do_fine from a file did not switch every property on")
      if (allocated(error)) return
      ! A single property flag stays single
      call configured(drop, '{"do_curvature":true}', error)
      if (allocated(error)) return
      call check(error, drop%do_curvature .and. .not. (drop%do_fine .or. drop%do_grid_density .or. &
                                                       drop%do_normal .or. drop%do_atomarea), &
                 more="do_curvature from a file switched other properties")
      if (allocated(error)) return
      call configured(iswig, '{"num_leb":194,"cut_a":0.3,"cut_f":0.4}', error)
      if (allocated(error)) return
      call check(error, iswig%num_leb, 194)
      if (allocated(error)) return
      call check(error, maxval(abs([iswig%cut_a, iswig%cut_f] - [0.3_wp, 0.4_wp])), &
                 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call configured(numsa, '{"num_leb":194,"probe":0.8,"offset":1.1,"smoothing":0.4,"tolsesp":0.02}', error)
      if (allocated(error)) return
      call check(error, numsa%num_leb, 194)
      if (allocated(error)) return
      call check(error, maxval(abs([numsa%probe, numsa%offset, numsa%smoothing, numsa%tolsesp] - &
                                   [0.8_wp, 1.1_wp, 0.4_wp, 0.02_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call configured(mc, '{"spacing":0.3,"obj_file":"one.obj","pqr_file":"two.pqr"}', error)
      if (allocated(error)) return
      call check(error, mc%spacing, 0.3_wp)
      if (allocated(error)) return
      call check(error, allocated(mc%obj_file) .and. allocated(mc%pqr_file))
      if (allocated(error)) return
      call check(error, mc%obj_file == "one.obj" .and. mc%pqr_file == "two.pqr")
      if (allocated(error)) return
      call configured(mt, '{"spacing":0.3,"target_area":0.5,"area_tolerance":0.2,"gaussian_scale":4.5,'// &
                      '"w_0ls_from":0.2,"w_0ls_to":0.7,"w_0ls_p":0.9,"w_0ls_a":1.7,"weight_cutoff":1e-6}', error)
      if (allocated(error)) return
      call check(error, maxval(abs([mt%spacing, mt%target_area, mt%area_tolerance, mt%gaussian_scale, &
                                    mt%w_0ls_from, mt%w_0ls_to, mt%w_0ls_p, mt%w_0ls_a, mt%weight_cutoff] - &
                                   [0.3_wp, 0.5_wp, 0.2_wp, 4.5_wp, 0.2_wp, 0.7_wp, 0.9_wp, 1.7_wp, 1.0e-6_wp])), &
                 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call configured(svdw, '{"blend_k":7,"blend_1b":0.6,"blend_2b":0.4,"blend_3b":2,'// &
                      '"fade_enabled":true,"fade_start":12,"fade_width":5}', error)
      if (allocated(error)) return
      call check(error, maxval(abs([svdw%blend_k, svdw%blend_1b, svdw%blend_2b, svdw%blend_3b, &
                                    svdw%fade_start, svdw%fade_width] - &
                                   [7.0_wp, 0.6_wp, 0.4_wp, 2.0_wp, 12.0_wp, 5.0_wp])), 0.0_wp, &
                 thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, svdw%fade_enabled, more="fade_enabled read from JSON")
      if (allocated(error)) return
      call configured(cfc, '{"a1":-12,"a2":-7,"c":3,"m":4,'// &
                      '"fade_enabled":true,"fade_start":20,"fade_width":4}', error)
      if (allocated(error)) return
      call check(error, maxval(abs([cfc%a1, cfc%a2, cfc%c, cfc%fade_start, cfc%fade_width] - &
                                   [-12.0_wp, -7.0_wp, 3.0_wp, 20.0_wp, 4.0_wp])), &
                 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, cfc%m, 4)
      if (allocated(error)) return
      call check(error, cfc%fade_enabled, more="CFC fade_enabled read from JSON")
      if (allocated(error)) return
      call configured(rho, '{"rho_iso":0.004,"scale":120,"log_grad_out":6,"log_grad_cusp":3,"exclusion_cap":4,' &
         //'"core_fill":false,"core_radius_fraction":0.3}', error)
      if (allocated(error)) return
      call check(error, maxval(abs([rho%rho_iso, rho%scale, rho%log_grad_out, rho%log_grad_cusp, rho%exclusion_cap, &
         rho%core_radius_fraction] - [0.004_wp, 120.0_wp, 6.0_wp, 3.0_wp, 4.0_wp, 0.3_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, .not. rho%core_fill, more="core_fill read from JSON")
      if (allocated(error)) return
      call configured(pcm, '{"solver":2,"solver_tol":0.00002,"solver_maxiter":71}', error)
      if (allocated(error)) return
      call check(error, pcm%solver == solver_type%iterative .and. pcm%solver_maxiter == 71)
      if (allocated(error)) return
      call check(error, pcm%solver_tol, 2.0e-5_wp, thr=1.0e-14_wp)

   contains

      !> Read a fixture, then exercise serialization with the supplied values
      !>
      !> @param[inout] param Settings being loaded
      !> @param[in] document Independent input fixture
      !> @param[out] error Test failure
      subroutine configured(param, document, error)
         !> Settings being loaded
         class(moist_model_parameters_type), intent(inout) :: param
         !> Independent input fixture
         character(len=*), intent(in) :: document
         !> Test failure
         type(error_type), allocatable, intent(out) :: error
         !> Library error
         type(moist_error), allocatable :: err
         !> Owned file unit and I/O status
         integer :: unit, stat, close_stat
         !> Fixture path
         character(len=*), parameter :: path = "moist-test-parameters-field-fixture.json"

         open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write (unit, "(a)", iostat=stat) document
         close (unit, iostat=close_stat)
         call check(error, stat == 0 .and. close_stat == 0)
         if (allocated(error)) return
         call param%read_file(path, err)
         open (newunit=unit, file=path, status="old", iostat=stat)
         if (stat == 0) close (unit, status="delete", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call check(error,.not. allocated(err))
         if (allocated(error)) return
         call roundtrip(param, "moist-test-parameters-all-fields", error)
         if (allocated(error)) return
         call check_reset(param, error)
      end subroutine configured

      !> Reset a populated copy and compare it with a freshly allocated instance
      !>
      !> @param[in] param Populated concrete settings
      !> @param[out] error Test failure
      subroutine check_reset(param, error)
         !> Populated concrete settings
         class(moist_model_parameters_type), intent(in) :: param
         !> Test failure
         type(error_type), allocatable, intent(out) :: error
         !> Independent objects with the same dynamic type
         class(moist_model_parameters_type), allocatable :: reset, fresh
         !> Library error
         type(moist_error), allocatable :: err
         !> Scratch units and read statuses
         integer :: reset_unit, fresh_unit, reset_stat, fresh_stat
         !> Serialized lines; no copied field/default tables
         character(len=4096) :: reset_line, fresh_line

         allocate (reset, source=param)
         allocate (fresh, mold=param)
         call reset%init_defaults()
         call reset%validate(err)
         call check(error,.not. allocated(err))
         if (allocated(error)) return
         call fresh%validate(err)
         call check(error,.not. allocated(err))
         if (allocated(error)) return
         open (newunit=reset_unit, status="scratch", action="readwrite", iostat=reset_stat)
         call check(error, reset_stat, 0)
         if (allocated(error)) return
         open (newunit=fresh_unit, status="scratch", action="readwrite", iostat=fresh_stat)
         if (fresh_stat /= 0) close (reset_unit)
         call check(error, fresh_stat, 0)
         if (allocated(error)) return
         ! Leave through the block so both scratch units are closed on failure
         compare: block
            call reset%print_parameters(err, reset_unit)
            call check(error,.not. allocated(err))
            if (allocated(error)) exit compare
            call fresh%print_parameters(err, fresh_unit)
            call check(error,.not. allocated(err))
            if (allocated(error)) exit compare
            rewind (reset_unit)
            rewind (fresh_unit)
            do
               read (reset_unit, "(a)", iostat=reset_stat) reset_line
               read (fresh_unit, "(a)", iostat=fresh_stat) fresh_line
               if (reset_stat < 0 .and. fresh_stat < 0) exit
               call check(error, reset_stat == 0 .and. fresh_stat == 0)
               if (allocated(error)) exit compare
               call check(error, trim(reset_line) == trim(fresh_line), &
                          more="default reset retained a configured field")
               if (allocated(error)) exit compare
            end do
         end block compare
         close (reset_unit, iostat=reset_stat)
         close (fresh_unit, iostat=fresh_stat)
         if (allocated(error)) return
         call check(error, reset_stat == 0 .and. fresh_stat == 0)
      end subroutine check_reset

   end subroutine test_configuration_fields

   !> Reject structurally valid files whose numerical settings are invalid
   subroutine test_parameter_validation(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(moist_error), allocatable :: err
      !> Parameter classes with numerical validation
      type(moist_cavity_drop_parameters_type) :: drop
      type(moist_cavity_mc_parameters_type) :: mc
      type(moist_cavity_mt_parameters_type) :: mt
      type(moist_pcm_parameters_type) :: pcm
      type(moist_cavity_drop_lsf_svdw_param_type) :: svdw
      type(moist_cavity_drop_lsf_cfc_param_type) :: cfc
      !> Invalid solver cases
      integer :: invalid

      ! SvdW fade band: read only while the fade is on (default on); the refusals
      ! are the bad_svdw_fade_* tests
      call check(error, svdw%fade_enabled, more="fade on by default")
      if (allocated(error)) return
      svdw%fade_enabled = .false.
      svdw%fade_start = 0.0_wp
      call svdw%validate(err)
      call check(error, .not. allocated(err), more="fade off must not validate the band")
      if (allocated(error)) return
      svdw%fade_enabled = .true.
      svdw%fade_start = 18.0_wp
      svdw%fade_width = 7.0_wp
      call svdw%validate(err)
      call check(error, .not. allocated(err), more="valid fade band refused")
      if (allocated(error)) return
      ! CFC fade band: on by default, read only while on; the refusals are the
      ! bad_cfc_fade_* tests
      call check(error, cfc%fade_enabled, more="CFC fade on by default")
      if (allocated(error)) return
      cfc%fade_enabled = .false.
      cfc%fade_start = 0.0_wp
      call cfc%validate(err)
      call check(error, .not. allocated(err), more="CFC fade off must not validate the band")
      if (allocated(error)) return
      cfc%fade_enabled = .true.
      cfc%fade_start = 27.0_wp
      call cfc%validate(err)
      call check(error, .not. allocated(err), more="CFC default fade band refused")
      if (allocated(error)) return
      mc%spacing = 0.0_wp
      call mc%validate(err)
      call check(error, allocated(err))
      if (allocated(error)) return
      do invalid = 1, 4
         call pcm%init_defaults()
         select case (invalid)
         case (1)
            pcm%solver = solver_type%dense - 1
         case (2)
            pcm%solver = solver_type%iterative + 1
         case (3)
            pcm%solver_tol = 0.0_wp
         case (4)
            pcm%solver_maxiter = 0
         case default
            error stop "invalid PCM test case"
         end select
         call pcm%validate(err)
         call check(error, allocated(err))
         if (allocated(error)) return
      end do
      ! File loading must dispatch the concrete validation hook
      call read_bad(drop, '{"tolerance":0.0}', error)
      if (allocated(error)) return
      call read_bad(mc, '{"spacing":0.0}', error)
      if (allocated(error)) return
      call read_bad(mt, '{"area_tolerance":1.5}', error)
      if (allocated(error)) return
      call read_bad(pcm, '{"solver":0}', error)
      if (allocated(error)) return
      call read_bad(svdw, '{"fade_enabled":true,"fade_width":-1}', error)
      if (allocated(error)) return
      call read_bad(cfc, '{"fade_enabled":true,"fade_start":-1}', error)

   contains

      !> Require an error when loading a numerically invalid document
      !>
      !> @param[inout] param Parameter instance
      !> @param[in] document Invalid input
      !> @param[out] error Test failure
      subroutine read_bad(param, document, error)
         !> Parameter instance
         class(moist_model_parameters_type), intent(inout) :: param
         !> Invalid input
         character(len=*), intent(in) :: document
         !> Test failure
         type(error_type), allocatable, intent(out) :: error
         !> Owned I/O state
         integer :: unit, stat, close_stat
         !> Owned file
         character(len=*), parameter :: path = "moist-test-parameters-validation.json"

         open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write (unit, "(a)", iostat=stat) document
         close (unit, iostat=close_stat)
         call check(error, stat == 0 .and. close_stat == 0)
         if (allocated(error)) return
         call param%read_file(path, err)
         open (newunit=unit, file=path, status="old", iostat=stat)
         if (stat == 0) close (unit, status="delete", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call check(error, allocated(err))
      end subroutine read_bad

   end subroutine test_parameter_validation

   !> Partial override constructors preserve fields omitted by their callers
   subroutine test_optional_overrides(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Shape parameter sets
      type(moist_cavity_drop_lsf_svdw_param_type) :: svdw
      type(moist_cavity_drop_lsf_cfc_param_type) :: cfc
      type(moist_cavity_drop_lsf_isodensity_param_type) :: rho

      call svdw%new(blend_k=7.0_wp, blend_1b=0.6_wp, blend_2b=0.4_wp, blend_3b=2.0_wp, &
                    fade_enabled=.true., fade_start=12.0_wp, fade_width=5.0_wp)
      call svdw%new(blend_k=8.0_wp)
      call check(error, maxval(abs([svdw%blend_k, svdw%blend_1b, svdw%blend_2b, svdw%blend_3b, &
                                    svdw%fade_start, svdw%fade_width] - &
                                   [8.0_wp, 0.6_wp, 0.4_wp, 2.0_wp, 12.0_wp, 5.0_wp])), 0.0_wp, &
                 thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, svdw%fade_enabled, more="partial override dropped fade_enabled")
      if (allocated(error)) return
      call cfc%new(a1=-12.0_wp, a2=-7.0_wp, c=3.0_wp, m=4, fade_enabled=.true., &
                   fade_start=20.0_wp, fade_width=4.0_wp)
      call cfc%new(a1=-11.0_wp)
      call check(error, maxval(abs([cfc%a1, cfc%a2, cfc%c, cfc%fade_start, cfc%fade_width] - &
                                   [-11.0_wp, -7.0_wp, 3.0_wp, 20.0_wp, 4.0_wp])), &
                 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, cfc%m, 4)
      if (allocated(error)) return
      call check(error, cfc%fade_enabled, more="CFC partial override dropped fade_enabled")
      if (allocated(error)) return
      call rho%new(rho_iso=0.004_wp, scale=120.0_wp, log_grad_out=6.0_wp, log_grad_cusp=3.0_wp, exclusion_cap=4.0_wp)
      call rho%new(scale=150.0_wp, core_fill=.false., core_radius_fraction=0.3_wp)
      call check(error, maxval(abs([rho%rho_iso, rho%scale, rho%log_grad_out, rho%log_grad_cusp, rho%exclusion_cap, &
         rho%core_radius_fraction] - [0.004_wp, 150.0_wp, 6.0_wp, 3.0_wp, 4.0_wp, 0.3_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, .not. rho%core_fill, more="core_fill override")
   end subroutine test_optional_overrides

end module test_parameters
