!> Native parameter construction, JSON/TOML input/output, and copy ownership
module test_parameters
   use, intrinsic :: iso_c_binding, only: c_null_funptr, c_null_ptr
   use mctc_env, only: wp, moist_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_model_parameters, only: moist_model_parameters_type
   use moist, only: moist_cavity_drop_parameters_type, moist_cavity_iswig_parameters_type, &
                    moist_cavity_numsa_parameters_type, moist_cavity_marchingcubes_parameters_type, &
                    moist_cavity_drop_lsf_svdw_param_type, moist_cavity_drop_lsf_cfc_param_type, &
                    moist_cavity_drop_lsf_isodensity_param_type, moist_pcm_parameters_type, moist_gostshyp_parameters_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_cavity_numsa, only: cavity_type_numsa, new_cavity_numsa
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_isodensity_callback, only: moist_cavity_drop_lsf_isodensity_callback_type
   use moist_model_continuum_component_type, only: model_continuum_component_type
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, new_component_cpcm
   use moist_model_continuum_component_gostshyp, only: model_continuum_component_gostshyp, new_component_gostshyp
   use moist_model_continuum_component_pcm_type, only: solver_type
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
      !> DROP settings reaching the conditional octree and switching rows
      type(moist_cavity_drop_parameters_type) :: param
      !> Level sets
      type(moist_cavity_drop_lsf_svdw_type) :: svdw
      type(moist_cavity_drop_lsf_cfc_type) :: cfc
      type(moist_cavity_drop_lsf_isodensity_callback_type) :: rho
      !> Cavities
      type(cavity_type_drop) :: drop
      type(cavity_type_iswig) :: iswig
      type(cavity_type_numsa) :: numsa
      !> Components
      type(model_continuum_component_cpcm) :: pcm
      type(model_continuum_component_gostshyp) :: gostshyp

      call new_context(ctx, nthreads=0, verbosity=0)
      call param%new(proj_level=9, wleb_prune_level=2, error=err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return

      call svdw%new()
      call new_cavity_drop(drop, default_cpcm_radii(), svdw, err, param, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_cavity(drop, [character(len=32) :: "Projection tol.", "Wleb switch from", "xi_born", &
         "Octree seed size", "Softmax scale", "Cell grid fraction", "Surface normals"], error)
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

      call new_component_cpcm(pcm, 80.0_wp, error=err, ctx=ctx)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call require_component(pcm, [character(len=32) :: "solver_maxiter"], error, &
         normal=[character(len=32) :: "solver_tol"])
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

      call reference%new(nleb=110, tolerance=3.0e-8_wp, proj_level=9, wleb_prune_level=3, error=err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return

      open (newunit=unit, file=toml_path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)", iostat=stat) "tolerance = 3e-8", "[grid]", "num_leb = 110", &
         "[projection]", "level = 9", "[switching]", "wleb_prune_level = 3"
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
         call check(error, maxval(abs([loaded%tolerance, loaded%proj_tol, loaded%iswig_xi_born, &
            loaded%wleb_prune_from, loaded%wleb_prune_to, loaded%wleb_cut, loaded%branch_dphi_max] - &
            [reference%tolerance, reference%proj_tol, reference%iswig_xi_born, reference%wleb_prune_from, &
            reference%wleb_prune_to, reference%wleb_cut, reference%branch_dphi_max])), 0.0_wp, thr=1.0e-14_wp)
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

   !> Collect parameter API regressions
   !>
   !> @param[out] testsuite Tests to run
   subroutine collect_parameters(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("constructors", test_constructors), &
         new_unittest("file-roundtrip", test_roundtrip), &
         new_unittest("memory-budget-gib", test_memory_budget_gib), &
         new_unittest("copies-and-errors", test_copies_and_errors), &
         new_unittest("vector-and-string", test_vector_and_string), &
         new_unittest("print-levels", test_print_levels), &
         new_unittest("print-detailed", test_print_detailed), &
         new_unittest("drop-load-file", test_drop_load_file), &
         new_unittest("bad_drop_load_unfitted_grid", test_bad_drop_load_unfitted_grid, should_fail=.true.), &
         new_unittest("bad_drop_load_missing_file", test_bad_drop_load_missing_file, should_fail=.true.), &
         new_unittest("invalid-fields", test_invalid_fields), &
         new_unittest("file-formats", test_file_formats), &
         new_unittest("drop-contracts", test_drop_contracts), &
         new_unittest("configuration-fields", test_configuration_fields), &
         new_unittest("parameter-validation", test_parameter_validation), &
         new_unittest("optional-overrides", test_optional_overrides)]
   end subroutine collect_parameters

   !> Integer GiB budgets survive table output and JSON/TOML round trips
   !>
   !> @param[out] error Test failure
   subroutine test_memory_budget_gib(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> PCM settings exercised through the shared parameter interface
      type(moist_pcm_parameters_type) :: pcm
      !> Library validation error
      type(moist_error), allocatable :: err
      !> Table printer
      type(prettyprinter) :: pp
      !> Practical GiB budgets from the minimum through one PiB
      integer, parameter :: budgets(*) = [1, 4, 5, 64, 1024, 1048576]
      !> Captured table lines and the current line without indentation
      character(len=256) :: lines(32), line
      !> Scratch I/O, indices, parsed value, and matching field count
      integer :: unit, stat, nline, icase, iline, pos, shown, matches, format_index

      call check(error, pcm%memory_budget_gib, 4)
      if (allocated(error)) return
      do icase = 1, size(budgets)
         pcm%memory_budget_gib = budgets(icase)
         open (newunit=unit, status="scratch", action="readwrite", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         pp = new_prettyprinter(unit=unit)
         call pcm%print_table(pp)
         call read_printout(unit, lines, nline)
         close (unit, iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call check(error, pcm%memory_budget_gib, budgets(icase))
         if (allocated(error)) return
         matches = 0
         do iline = 1, nline
            line = adjustl(lines(iline))
            if (index(line, "memory_budget_gib ") /= 1) cycle
            matches = matches + 1
            pos = index(line, "...", back=.true.)
            call check(error, pos > 0)
            if (allocated(error)) return
            read (line(pos + 3:), *, iostat=stat) shown
            call check(error, stat, 0)
            if (allocated(error)) return
            call check(error, shown, budgets(icase))
            if (allocated(error)) return
         end do
         call check(error, matches, 1)
         if (allocated(error)) return
         do format_index = 1, 2
            call roundtrip(pcm, "moist-test-memory-budget", error, format_index)
            if (allocated(error)) return
            call check(error, pcm%memory_budget_gib, budgets(icase))
            if (allocated(error)) return
         end do
      end do
      ! Exercise integer limits in serialization without relying on printer width
      pcm%memory_budget_gib = huge(0)
      do format_index = 1, 2
         call roundtrip(pcm, "moist-test-memory-budget", error, format_index)
         if (allocated(error)) return
         call check(error, pcm%memory_budget_gib, huge(0))
         if (allocated(error)) return
      end do
      do icase = -1, 0
         pcm%memory_budget_gib = icase
         call pcm%validate(err)
         call check(error, allocated(err))
         if (allocated(error)) return
      end do
      call pcm%init_defaults()
      call check(error, pcm%memory_budget_gib, 4)
   end subroutine test_memory_budget_gib

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
      param%proj_tol = -1.0_wp
      param%do_fine = .true.
      call new_cavity_drop(drop, default_cpcm_radii(), lsf, err, param, ctx=ctx)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, drop%param%proj_tol, param%tolerance)
      if (allocated(error)) return
      call check(error, param%proj_tol, -1.0_wp)
      if (allocated(error)) return
      call check(error, drop%param%do_grid_density .and. drop%param%do_curvature .and. &
                 drop%param%do_normal .and. drop%param%do_r_iI .and. drop%param%do_rho, &
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
                              drop%param%do_normal .or. drop%param%do_r_iI .or. drop%param%do_rho), &
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
      call check(error, pcm%param%solver, solver_type%dense)
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
      type(moist_cavity_marchingcubes_parameters_type) :: mc
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
      do format_index = 1, 2
         call roundtrip(drop, stem, error, format_index)
         if (allocated(error)) return
         call check(error, drop%proj_tol, 1.0e-8_wp)
         if (allocated(error)) return
         call check(error, drop%do_fine)
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
      pcm%memory_budget_gib = 5
      call roundtrip(pcm, stem, error)
      if (allocated(error)) return
      call check(error, pcm%solver_maxiter, 42)
      if (allocated(error)) return
      call check(error, pcm%memory_budget_gib, 5)
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
      call check(error, param%proj_tol, 2.0e-8_wp)
      if (allocated(error)) return
      call check(error, param%proj_maxiter, 150)
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

      call param%new(nleb=110, tolerance=2.0e-7_wp, proj_maxiter=81, proj_level=9, &
                     branch_weight_s=0.03_wp, rho_grid_h=0.75_wp, wleb_prune_level=0, error=err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, param%num_leb == 110 .and. param%proj_maxiter == 81 .and. param%proj_level == 9)
      if (allocated(error)) return
      call check(error, param%wleb_cut, 1.0e-8_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call check(error, param%screening_threshold, 2.0e-8_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call check(error, param%proj_tol, 2.0e-7_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call check(error, param%branch_sep_cut, 2.0e-6_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      ! Adjacency cutoff is four kernel lengths; the soft density screens with
      ! cutoff*rho_grid_h as 4 h^2
      call check(error, param%adj_list_grid_cutoff, 4.0_wp*0.75_wp)
      if (allocated(error)) return
      call check(error, param%iswig_xi_born > 4.8_wp .and. param%iswig_xi_born < 5.0_wp)
      if (allocated(error)) return
      ! The unnormalized softmax reaches the fallback floor at the bound
      call check(error, exp(-param%branch_dphi_max/0.03_wp), 1.0e-8_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call check(error, param%branch_weight_floor, 0.0_wp)
      if (allocated(error)) return
      param%branch_weight_floor = 0.25_wp
      param%phi_alpha = 3.0_wp
      call param%compute_derived(err)
      call check(error,.not. allocated(err))
      if (allocated(error)) return
      call check(error, exp(-param%branch_dphi_max/0.03_wp), 0.25_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      ! Quadratic objective excess equals the admissible objective bound
      call check(error, 1.5_wp*(param%rho_max_from(2.0_wp)**2 - 4.0_wp), &
                 param%branch_dphi_max, thr=1.0e-14_wp)
      if (allocated(error)) return
      do level = 0, 6
         call param%new(wleb_prune_level=level, error=err)
         call check(error,.not. allocated(err))
         if (allocated(error)) return
         if (level == 0) then
            call check(error, param%wleb_prune_from == 0.0_wp .and. param%wleb_prune_to == 0.0_wp)
         else
            lower = 10.0_wp**(-14 + 2*level)
            call check(error, param%wleb_prune_from/lower, 1.0_wp, thr=1.0e-14_wp)
            if (allocated(error)) return
            call check(error, param%wleb_prune_to/lower, 100.0_wp, thr=1.0e-12_wp)
            if (allocated(error)) return
            call check(error, param%wleb_cut/lower, 1.0_wp, thr=1.0e-14_wp)
         end if
         if (allocated(error)) return
      end do
      do invalid = 1, 10
         call param%init_defaults()
         select case (invalid)
         case (1)
            param%tolerance = 0.0_wp
         case (2)
            param%rho_grid_h = 0.0_wp
         case (3)
            param%proj_maxiter = 0
         case (4)
            param%proj_level = 0
         case (5)
            param%proj_level = 10
         case (6)
            param%branch_weight_floor = 1.0_wp
         case (7)
            param%phi_alpha = 0.0_wp
         case (8)
            param%branch_weight_s = -0.1_wp
         case (9)
            param%wleb_prune_level = 7
         case (10)
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
      type(moist_cavity_marchingcubes_parameters_type) :: mc
      type(moist_cavity_drop_lsf_svdw_param_type) :: svdw
      type(moist_cavity_drop_lsf_cfc_param_type) :: cfc
      type(moist_cavity_drop_lsf_isodensity_param_type) :: rho
      type(moist_pcm_parameters_type) :: pcm

      call configured(drop, '{"tolerance":2e-7,"do_fine":true,"grid":{"num_leb":110},"objective":{"alpha":2},'// &
                      '"projection":{"maxiter":72,"level":9,'// &
                      '"octree":{"seed_size":0.7,"max_boxes":123,"max_survivors":43,"max_depth":8,"seed_mode":2}},'// &
                      '"screening":{"cell_grid_full_scan_below":19,"cell_grid_fraction":0.4},'// &
                      '"switching":{"w_0ls_from":0.2,"w_0ls_to":0.7,"w_0ls_p":0.9,"w_0ls_a":1.7,'// &
                      '"w_0tra_from":0.15,"w_0tra_to":0.4,"wleb_prune_level":2},'// &
                      '"density":{"rho_grid_h":0.8},"branching":{"softmax_scale":0.04,"weight_floor":0.2}}', error)
      if (allocated(error)) return
      call check(error, drop%num_leb == 110 .and. drop%do_fine .and. &
                 drop%proj_maxiter == 72 .and. drop%proj_level == 9 .and. &
                 drop%octree_max_boxes == 123 .and. drop%octree_max_survivors == 43 .and. &
                 drop%octree_max_depth == 8 .and. drop%octree_seed_mode == 2 .and. &
                 drop%cell_grid_full_scan_below == 19 .and. drop%wleb_prune_level == 2)
      if (allocated(error)) return
      call check(error, maxval(abs([drop%tolerance, drop%phi_alpha, drop%octree_seed_size, &
                                    drop%cell_grid_fraction, drop%w_0ls_from, drop%w_0ls_to, drop%w_0ls_p, drop%w_0ls_a, &
                                    drop%w_0tra_from, drop%w_0tra_to, drop%rho_grid_h, drop%branch_weight_s, &
                                    drop%branch_weight_floor] - &
                                   [2.0e-7_wp, 2.0_wp, 0.7_wp, 0.4_wp, 0.2_wp, 0.7_wp, 0.9_wp, 1.7_wp, 0.15_wp, 0.4_wp, &
                                    0.8_wp, 0.04_wp, 0.2_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, drop%do_grid_density .and. drop%do_curvature .and. drop%do_normal .and. &
                 drop%do_r_iI .and. drop%do_rho, more="do_fine from a file did not switch every property on")
      if (allocated(error)) return
      ! A single property flag stays single
      call configured(drop, '{"do_curvature":true}', error)
      if (allocated(error)) return
      call check(error, drop%do_curvature .and. .not. (drop%do_fine .or. drop%do_grid_density .or. &
                                                       drop%do_normal .or. drop%do_r_iI .or. drop%do_rho), &
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
      call configured(svdw, '{"blend_k":7,"blend_1b":0.6,"blend_2b":0.4,"blend_3b":2}', error)
      if (allocated(error)) return
      call check(error, maxval(abs([svdw%blend_k, svdw%blend_1b, svdw%blend_2b, svdw%blend_3b] - &
                                   [7.0_wp, 0.6_wp, 0.4_wp, 2.0_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call configured(cfc, '{"a1":-12,"a2":-7,"c":3,"m":6}', error)
      if (allocated(error)) return
      call check(error, maxval(abs([cfc%a1, cfc%a2, cfc%c] - [-12.0_wp, -7.0_wp, 3.0_wp])), &
                 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, cfc%m, 6)
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
      type(moist_cavity_marchingcubes_parameters_type) :: mc
      type(moist_pcm_parameters_type) :: pcm
      !> Invalid solver cases
      integer :: invalid

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
      call read_bad(pcm, '{"solver":0}', error)

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

      call svdw%new(blend_k=7.0_wp, blend_1b=0.6_wp, blend_2b=0.4_wp, blend_3b=2.0_wp)
      call svdw%new(blend_k=8.0_wp)
      call check(error, maxval(abs([svdw%blend_k, svdw%blend_1b, svdw%blend_2b, svdw%blend_3b] - &
                                   [8.0_wp, 0.6_wp, 0.4_wp, 2.0_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call cfc%new(a1=-12.0_wp, a2=-7.0_wp, c=3.0_wp, m=6)
      call cfc%new(a1=-11.0_wp)
      call check(error, maxval(abs([cfc%a1, cfc%a2, cfc%c] - [-11.0_wp, -7.0_wp, 3.0_wp])), &
                 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, cfc%m, 6)
      if (allocated(error)) return
      call rho%new(rho_iso=0.004_wp, scale=120.0_wp, log_grad_out=6.0_wp, log_grad_cusp=3.0_wp, exclusion_cap=4.0_wp)
      call rho%new(scale=150.0_wp, core_fill=.false., core_radius_fraction=0.3_wp)
      call check(error, maxval(abs([rho%rho_iso, rho%scale, rho%log_grad_out, rho%log_grad_cusp, rho%exclusion_cap, &
         rho%core_radius_fraction] - [0.004_wp, 150.0_wp, 6.0_wp, 3.0_wp, 4.0_wp, 0.3_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, .not. rho%core_fill, more="core_fill override")
   end subroutine test_optional_overrides

end module test_parameters
