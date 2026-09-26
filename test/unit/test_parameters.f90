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
   use moist_model_component_pcm_cpcm, only: solvation_model_component_cpcm, new_component_cpcm
   use moist_model_component_pcm_type, only: solver_type
   use moist_context, only: moist_context_type, new_context
   use moist_radii, only: default_cpcm_radii
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
   !>
   !> @param[out] error Test failure
   subroutine test_vector_and_string(error)
      type(error_type), allocatable, intent(out) :: error
      type(vector_parameters_type) :: param

      param%weights = [0.3_wp, -5.0_wp]
      param%label = 'a "quoted" name'
      call roundtrip(param, "moist-test-parameters-vector", error)
      if (allocated(error)) return
      call check(error, maxval(abs(param%weights - [0.3_wp, -5.0_wp])), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, trim(param%label) == 'a "quoted" name')
   end subroutine test_vector_and_string

   !> Structurally wrong fields are refused by name
   !>
   !> A group key holding a scalar, a vector of the wrong length and an
   !> oversized fixed-length string each fail in their own registration helper
   !>
   !> @param[out] error Test failure
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
         new_unittest("file-formats", test_file_formats)]
   end subroutine collect_parameters

   !> Check copied settings, derived values, and resetting to defaults
   !>
   !> @param[out] error Test failure
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
      type(solvation_model_component_cpcm) :: pcm

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
      call check(error, drop%request%curvature)
      if (allocated(error)) return
      call new_cavity_drop(drop, ctx, default_cpcm_radii(), lsf, err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, drop%param%num_leb, 194)
      if (allocated(error)) return
      call check(error, .not. drop%request%curvature)
      if (allocated(error)) return
      call new_cavity_iswig(iswig, ctx, default_cpcm_radii(), err, &
         moist_cavity_iswig_parameters_type(num_leb=194, cut_f=0.2_wp))
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call new_cavity_iswig(iswig, ctx, default_cpcm_radii(), err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, iswig%num_leb, 110)
      if (allocated(error)) return
      call check(error, iswig%cut_f, 1.0e-10_wp)
      if (allocated(error)) return
      call lsf%new()
      call check(error, lsf%param%blend_k, 5.5_wp)
      if (allocated(error)) return
      call rho%new(c_null_funptr, c_null_ptr, &
         moist_cavity_drop_lsf_isodensity_param_type(rho_iso=0.004_wp, exclusion_cap=4.0_wp))
      call check(error, rho%param%exclusion_cap, 4.0_wp)
      if (allocated(error)) return
      call new_component_cpcm(pcm, ctx, 80.0_wp, error=err, &
         param=moist_pcm_parameters_type(solver=solver_type%iterative, solver_tol=1.0e-8_wp, solver_maxiter=50))
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, pcm%solver_maxiter, 50)
      if (allocated(error)) return
      call new_component_cpcm(pcm, ctx, 80.0_wp, error=err)
      call check(error, .not. allocated(err))
      if (allocated(error)) return
      call check(error, pcm%solver, solver_type%cholesky)
      if (allocated(error)) return
      param%num_leb = 1
      call new_cavity_drop(drop, ctx, default_cpcm_radii(), lsf, err, param)
      call check(error, allocated(err))
      call ctx%delete()
   end subroutine test_constructors

   !> All concrete parameter sets support the same file and print operations
   !>
   !> @param[out] error Test failure
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

      drop%tolerance = 1.0e-8_wp
      drop%do_fine = .true.
      call roundtrip(drop, stem, error)
      if (allocated(error)) return
      call check(error, drop%proj_tol, 1.0e-8_wp)
      if (allocated(error)) return
      call check(error, drop%do_fine)
      if (allocated(error)) return
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
   subroutine roundtrip(param, stem, error)
      class(moist_model_parameters_type), intent(inout) :: param
      character(len=*), intent(in) :: stem
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      integer :: unit, stat, i
      character(len=5), parameter :: extensions(2) = [".json", ".toml"]

      do i = 1, size(extensions)
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
   !>
   !> @param[out] error Test failure
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
   !>
   !> @param[out] error Test failure
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

end module test_parameters
