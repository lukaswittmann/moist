!> Native parameter construction, JSON/TOML input/output, and copy ownership
module test_parameters
   use, intrinsic :: iso_c_binding, only: c_null_funptr, c_null_ptr
   use mctc_env, only: wp, moist_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_model_parameters, only: moist_model_parameters_type
   use moist, only: moist_cavity_drop_parameters_type, moist_cavity_iswig_parameters_type, &
      moist_cavity_numsa_parameters_type, moist_cavity_marchingcubes_parameters_type, &
      moist_cavity_drop_lsf_svdw_param_type, moist_cavity_drop_lsf_cfc_param_type, &
      moist_cavity_drop_lsf_isodensity_param_type, moist_pcm_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_isodensity_callback, only: moist_cavity_drop_lsf_isodensity_callback_type
   use moist_model_continuum_component_pcm_cpcm, only: model_continuum_component_cpcm, new_component_cpcm
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
      call self%register_string("label", self%label)
   end subroutine register_vector_entries

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

      open(newunit=unit, file=path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write(unit, "(a)", iostat=stat) '{"weights":[0.3,-5.0],"label":"input"}'
      close(unit, iostat=close_stat)
      call check(error, stat == 0 .and. close_stat == 0)
      if (allocated(error)) return
      call param%read_file(path, err)
      open(newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close(unit, status="delete", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      call check(error, .not. allocated(err))
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

         open(newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write(unit, "(a)", iostat=stat) text
         close(unit)
         call check(error, stat, 0)
         if (allocated(error)) return
         call param%read_file(path, err)
         open(newunit=unit, file=path, status="old", iostat=stat)
         if (stat == 0) close(unit, status="delete", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call check(error, allocated(err), more="accepted "//text)
         if (allocated(error)) return
         call check(error, index(err%message, expected) > 0, more=err%message)
      end subroutine read_invalid

   end subroutine test_invalid_fields

   !> Collect parameter API regressions
   !>
   !> @param[out] testsuite Tests to run
   subroutine collect_parameters(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("constructors", test_constructors), &
         new_unittest("file-roundtrip", test_roundtrip), &
         new_unittest("copies-and-errors", test_copies_and_errors), &
         new_unittest("vector-and-string", test_vector_and_string), &
         new_unittest("invalid-fields", test_invalid_fields), &
         new_unittest("file-formats", test_file_formats), &
         new_unittest("drop-contracts", test_drop_contracts), &
         new_unittest("configuration-fields", test_configuration_fields), &
         new_unittest("parameter-validation", test_parameter_validation), &
         new_unittest("optional-overrides", test_optional_overrides)]
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

      call new_context(ctx, verbosity=0)
      shape%blend_k = 7.0_wp
      call lsf%new(shape)
      shape%blend_k = 9.0_wp
      call check(error, lsf%param%blend_k, 7.0_wp)
      if (allocated(error)) return
      param%num_leb = 110
      param%tolerance = 1.0e-8_wp
      param%proj_tol = -1.0_wp
      param%do_fine = .true.
      call new_cavity_drop(drop, ctx, default_cpcm_radii(), lsf, err, param)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, drop%param%proj_tol, param%tolerance)
      if (allocated(error)) return
      call check(error, param%proj_tol, -1.0_wp)
      if (allocated(error)) return
      call check(error, drop%param%do_grid_density .and. drop%param%do_curvature .and. &
         drop%param%do_normal .and. drop%param%do_r_iI .and. drop%param%do_rho, &
         more="do_fine did not switch every optional property on")
      if (allocated(error)) return
      call check(error, .not. param%do_curvature, more="constructor changed the caller's settings")
      if (allocated(error)) return
      call check(error, drop%lsf_model%screening_threshold, 1.0e-9_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call new_cavity_drop(drop, ctx, default_cpcm_radii(), lsf, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, drop%param%num_leb, 194)
      if (allocated(error)) return
      call check(error, .not. (drop%param%do_grid_density .or. drop%param%do_curvature .or. &
         drop%param%do_normal .or. drop%param%do_r_iI .or. drop%param%do_rho), &
         more="default DROP settings switch an optional property on")
      if (allocated(error)) return
      call new_cavity_iswig(iswig, ctx, default_cpcm_radii(), err, &
         moist_cavity_iswig_parameters_type(num_leb=194, cut_a=0.3_wp, cut_f=0.2_wp))
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, iswig%param%num_leb, 194)
      if (allocated(error)) return
      call check(error, iswig%param%cut_a, 0.3_wp)
      if (allocated(error)) return
      call check(error, iswig%param%cut_f, 0.2_wp)
      if (allocated(error)) return
      call new_cavity_iswig(iswig, ctx, default_cpcm_radii(), err)
      call check(error, .not. allocated(err))
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
      call new_component_cpcm(pcm, ctx, 80.0_wp, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%iterative, solver_tol=1.0e-8_wp, solver_maxiter=71))
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, pcm%param%solver, solver_type%iterative)
      if (allocated(error)) return
      call check(error, pcm%param%solver_tol, 1.0e-8_wp, thr=1.0e-15_wp)
      if (allocated(error)) return
      call check(error, pcm%param%solver_maxiter, 71)
      if (allocated(error)) return
      call new_component_cpcm(pcm, ctx, 80.0_wp, error=err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, pcm%param%solver, solver_type%cholesky)
      if (allocated(error)) return
      param%num_leb = 1
      call new_cavity_drop(drop, ctx, default_cpcm_radii(), lsf, err, param)
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
      call check(error, .not. allocated(mc%pqr_file))
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
      pcm%solver = solver_type%lu
      pcm%solver_maxiter = 42
      call roundtrip(pcm, stem, error)
      if (allocated(error)) return
      call check(error, pcm%solver_maxiter, 42)
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
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call param%init_defaults()
         call param%read_file(stem//extensions(i), err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         open(newunit=unit, status="scratch", action="readwrite", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call param%print_parameters(err, unit)
         close(unit, iostat=stat)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call check(error, stat, 0)
         if (allocated(error)) return
         open(newunit=unit, file=stem//extensions(i), status="old", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         close(unit, status="delete", iostat=stat)
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
         "moist-test-parameters-formats.txt", "moist-test-parameters-formats", "directory.json/no-extension"]
      logical :: exists
      !> Uppercase extension; must not match another test's file on case-insensitive systems
      character(len=*), parameter :: toml_file = "moist-test-parameters-formats.TOML"

      open(newunit=unit, file=toml_file, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write(unit, "(a)", iostat=stat) "# Independently authored TOML", &
         "tolerance = 2e-8", "[grid]", "num_leb = 110"
      close(unit, iostat=close_stat)
      call check(error, stat == 0 .and. close_stat == 0)
      if (allocated(error)) return
      call param%read_file(toml_file, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, param%num_leb, 110)
      if (allocated(error)) return
      call check(error, param%proj_tol, 2.0e-8_wp)
      if (allocated(error)) return
      call check(error, param%proj_maxiter, 150)
      if (allocated(error)) return
      ! Uppercase output extensions select TOML too
      call param%write_file(toml_file, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call param%read_file(toml_file, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      open(newunit=unit, file=toml_file, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write(unit, "(a)", iostat=stat) "num_leb = ["
      close(unit, iostat=close_stat)
      call check(error, stat == 0 .and. close_stat == 0)
      if (allocated(error)) return
      call iswig%read_file(toml_file, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      open(newunit=unit, file=toml_file, status="old", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      close(unit, status="delete", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      call iswig%read_file(toml_file, err)
      call check(error, allocated(err))
      if (allocated(error)) return
      do i = 1, size(bad_paths)
         call iswig%write_file(trim(bad_paths(i)), err)
         call check(error, allocated(err))
         if (allocated(error)) return
         inquire(file=trim(bad_paths(i)), exist=exists)
         call check(error, .not. exists)
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
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      copied = source
      copied%tolerance = 1.0e-6_wp
      call copied%write_file(parameter_file, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call copied%read_file(parameter_file, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, copied%tolerance, 1.0e-6_wp)
      if (allocated(error)) return
      call check(error, source%tolerance, 1.0e-8_wp)
      if (allocated(error)) return
      do i = 1, size(invalid)
         open(newunit=unit, file=parameter_file, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write(unit, "(a)", iostat=stat) trim(invalid(i))
         close(unit)
         call check(error, stat, 0)
         if (allocated(error)) return
         call iswig%read_file(parameter_file, err)
         call check(error, allocated(err))
         if (allocated(error)) return
      end do
      open(newunit=unit, file=parameter_file, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write(unit, "(a)", iostat=stat) '{"num_leb": 194}'
      close(unit)
      call check(error, stat, 0)
      if (allocated(error)) return
      iswig%cut_f = 0.5_wp
      call iswig%read_file(parameter_file, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, iswig%num_leb, 194)
      if (allocated(error)) return
      call check(error, iswig%cut_f, 1.0e-10_wp)
      if (allocated(error)) return
      call iswig%write_file("missing-moist-parameter-directory/config.json", err)
      call check(error, allocated(err))
      if (allocated(error)) return
      open(newunit=unit, file=parameter_file, status="old", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      close(unit, status="delete", iostat=stat)
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
      call check(error, .not. allocated(err))
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
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, exp(-param%branch_dphi_max/0.03_wp), 0.25_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      ! Quadratic objective excess equals the admissible objective bound
      call check(error, 1.5_wp*(param%rho_max_from(2.0_wp)**2 - 4.0_wp), &
         param%branch_dphi_max, thr=1.0e-14_wp)
      if (allocated(error)) return
      do level = 0, 6
         call param%new(wleb_prune_level=level, error=err)
         call check(error, .not. allocated(err))
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
      call configured(rho, '{"rho_iso":0.004,"scale":120,"log_grad_out":6,"log_grad_cusp":3,"exclusion_cap":4}', error)
      if (allocated(error)) return
      call check(error, maxval(abs([rho%rho_iso, rho%scale, rho%log_grad_out, rho%log_grad_cusp, rho%exclusion_cap] - &
         [0.004_wp, 120.0_wp, 6.0_wp, 3.0_wp, 4.0_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call configured(pcm, '{"solver":2,"solver_tol":0.00002,"solver_maxiter":71}', error)
      if (allocated(error)) return
      call check(error, pcm%solver == solver_type%lu .and. pcm%solver_maxiter == 71)
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

         open(newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write(unit, "(a)", iostat=stat) document
         close(unit, iostat=close_stat)
         call check(error, stat == 0 .and. close_stat == 0)
         if (allocated(error)) return
         call param%read_file(path, err)
         open(newunit=unit, file=path, status="old", iostat=stat)
         if (stat == 0) close(unit, status="delete", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call check(error, .not. allocated(err))
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

         allocate(reset, source=param)
         allocate(fresh, mold=param)
         call reset%init_defaults()
         call reset%validate(err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         call fresh%validate(err)
         call check(error, .not. allocated(err))
         if (allocated(error)) return
         open(newunit=reset_unit, status="scratch", action="readwrite", iostat=reset_stat)
         call check(error, reset_stat, 0)
         if (allocated(error)) return
         open(newunit=fresh_unit, status="scratch", action="readwrite", iostat=fresh_stat)
         if (fresh_stat /= 0) close(reset_unit)
         call check(error, fresh_stat, 0)
         if (allocated(error)) return
         ! Leave through the block so both scratch units are closed on failure
         compare: block
            call reset%print_parameters(err, reset_unit)
            call check(error, .not. allocated(err))
            if (allocated(error)) exit compare
            call fresh%print_parameters(err, fresh_unit)
            call check(error, .not. allocated(err))
            if (allocated(error)) exit compare
            rewind(reset_unit)
            rewind(fresh_unit)
            do
               read(reset_unit, "(a)", iostat=reset_stat) reset_line
               read(fresh_unit, "(a)", iostat=fresh_stat) fresh_line
               if (reset_stat < 0 .and. fresh_stat < 0) exit
               call check(error, reset_stat == 0 .and. fresh_stat == 0)
               if (allocated(error)) exit compare
               call check(error, trim(reset_line) == trim(fresh_line), &
                  more="default reset retained a configured field")
               if (allocated(error)) exit compare
            end do
         end block compare
         close(reset_unit, iostat=reset_stat)
         close(fresh_unit, iostat=fresh_stat)
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
            pcm%solver = solver_type%inversion - 1
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

         open(newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write(unit, "(a)", iostat=stat) document
         close(unit, iostat=close_stat)
         call check(error, stat == 0 .and. close_stat == 0)
         if (allocated(error)) return
         call param%read_file(path, err)
         open(newunit=unit, file=path, status="old", iostat=stat)
         if (stat == 0) close(unit, status="delete", iostat=stat)
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
      call rho%new(scale=150.0_wp)
      call check(error, maxval(abs([rho%rho_iso, rho%scale, rho%log_grad_out, rho%log_grad_cusp, rho%exclusion_cap] - &
         [0.004_wp, 150.0_wp, 6.0_wp, 3.0_wp, 4.0_wp])), 0.0_wp, thr=1.0e-14_wp)
   end subroutine test_optional_overrides

end module test_parameters
