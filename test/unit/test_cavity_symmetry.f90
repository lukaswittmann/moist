!> Tests for the symmetry model and the cavities oriented by it (iSwiG, NumSA, DROP)
!>
!> Covers
!>
!> - settings and model: the top-level `symmetry` table (TOML, defaults,
!>   validation, thresholds), species from radii, tolerance and errors, reuse
!>   of a host's update, timers, `symmetry_*` fields, report and warnings
!> - orientation: a cavity without a model or in the lab orientation is
!>   bitwise c + r u; molecular points are c + r R u, R from the math layer
!> - exact symmetry: with a group the points map onto themselves
!> - frame gradient: analytic gradients against central differences at the
!>   frame levels eigen, plane and sphere, and with a group
!> - rigid motion: gradients along finite rigid paths, and cavities of
!>   rotated, translated and renumbered copies compared end to end
!> - special cases: linear molecules, a single atom, tied frames, frames
!>   near a threshold, further point groups and subgroups
!>
!> Conventions
!>
!> - cavity kinds `kind_iswig`, `kind_numsa`, `kind_drop`; tests named
!>   `drop_*` and `numsa_*` cover that kind, the others iSwiG or several
!> - shared checks take the kind or a fixture: `check_frame_gradient`
!>   (`frame_fixture`), `check_rigid_paths`, `check_e2e_rotation`
!>   (`e2e_state`), `check_exact_symmetry`
!> - cavities are built like a host builds them: `host_model`, then the
!>   public constructor (`updated_cavity`)
!> - controls that must fail: the lab orientation end to end
!>   (`*_lab_sensitive`), a gradient without its frame term (`frame_off`)
!> - tolerances: `TOL` (symmetry tolerance of the settings, bohr), `FD_*`
!>   (differences), `E2E_TOL` and `E2E_TOL_DROP` (end to end, relative)
!> - structures and rigid motions come from `test_cavity_symmetry_fixtures`;
!>   math-level tests (detection, embedding, frame) are in `test_math_symmetry`
module test_cavity_symmetry
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use mctc_io_constants, only: pi
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_grid_frame, only: grid_frame_type
   use moist_cavity_adjoint, only: cavity_adjoint_type
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm
   use moist_math_symmetry_frame, only: frame_level_none, frame_level_eigen, frame_level_plane, frame_level_sphere
   use moist_cavity_grid_model_symmetry, only: grid_model_symmetry_type, new_grid_model_symmetry, &
      & moist_cavity_grid_model_symmetry_parameters_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, moist_cavity_iswig_parameters_type
   use moist_cavity_numsa, only: cavity_type_numsa, new_cavity_numsa, moist_cavity_numsa_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_channels_fields, only: field_query_type
   use moist_context, only: moist_context_type, new_context
   use moist_radii, only: default_cpcm_radii, new_radii_custom_atoms, radius_type
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_type, new_lebedev_grid
   use moist_math_symmetry_orient, only: sym_orient_type, new_sym_orient, sym_orient_vjp, &
      & orient_lab, orient_molecular
   use test_helpers, only: read_printout, printed_entry, component_view, check_moist_error, &
      & fd4_scalar, fd4_offsets
   use test_model_component_pcm_fixtures, only: pcm_cpcm, new_pcm, stage_energy, answer_sources
   use test_cavity_drop_fd, only: check_drop_row_gradient, drop_fd_param, drop_fd_svdw
   use test_cavity_symmetry_fixtures, only: model_nrot, model_rot_index, model_rotation, model_rotate_grid, &
      & model_add_torque, model_frame_jacobian, &
      & symmetry_of, structure_of, test_rotation, test_permutation, rigid_motion, &
      & rotation_matrix, new_tilted, &
      & surface_points_type, surface_points_of, match_surface_points, molecule_symmetry, same_bits, &
      & gradient_covariance_error, sphere_gradient_covariance_error
   use moist_math_linalg_geometry, only: identity3, cross_product
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: collect_cavity_symmetry

   !> Exactness tolerance of the symmetric tests, bohr
   real(wp), parameter :: TOL = 1.0e-8_wp
   !> Length of a captured printout line
   integer, parameter :: line_len = 200

   !> Cavities with symmetry handling
   integer, parameter :: kind_iswig = 1, kind_numsa = 2, kind_drop = 3
   !> Their names in messages
   character(len=*), parameter :: kind_name(3) = [character(len=5) :: "iSwiG", "NumSA", "DROP"]

   !> Quantities of the frame-gradient checks
   !>
   !> Area, volume, generic surface adjoint energy, CPCM energy
   integer, parameter :: q_area = 1, q_volume = 2, q_adjoint = 3, q_cpcm = 4, nquant = 4
   !> Analytic gradients of the frame-gradient checks
   integer, parameter :: ngrad = 8
   !> Quantity each analytic gradient differentiates
   integer, parameter :: grad_quantity(ngrad) = [q_area, q_area, q_volume, q_volume, q_volume, &
      & q_adjoint, q_cpcm, q_cpcm]
   !> Their names in messages
   character(len=*), parameter :: grad_name(ngrad) = [character(len=26) :: "area_grad", &
      & "area reverse", "volume_grad", "sum of v1_rA", "volume reverse (PV)", &
      & "surface adjoint reverse", "CPCM reverse", "CPCM forward (xyz1/f1_rA)"]
   !> Raw adjoint channels: xi, f, a, w, k1, k2, then xyz and n (three each)
   integer, parameter :: nchannel = 12
   !> Central-difference step of the frame-gradient checks, bohr
   real(wp), parameter :: FD_STEP = 1.0e-4_wp
   !> Tolerances of the frame-gradient checks, as in the iSwiG and DROP suites
   real(wp), parameter :: FD_ABS = 5.0e-9_wp, FD_REL = 5.0e-8_wp
   !> Central-difference step of the DROP frame-gradient checks, bohr
   real(wp), parameter :: FD_STEP_DROP = 1.0e-4_wp
   !> Absolute tolerance of the NumSA frame-gradient checks, that of its suite
   real(wp), parameter :: FD_ABS_NUMSA = sqrt(epsilon(1.0_wp))
   !> Rounding of a differenced value, in units of epsilon times its size
   !>
   !> - a sum over the points of the surface, another set and order at every
   !>   stencil point
   !> - differences of the rigid paths miss by up to 20 of these units at
   !>   the step used and by up to 44 at half and twice that step (GCC,
   !>   macOS arm64): a factor of 13 and of 6 below this
   real(wp), parameter :: FD_ROUND = 256.0_wp
   !> Sum of the absolute weights of the four-point stencil, (1 + 8 + 8 + 1)/12
   !>
   !> - rounding `d` of every value costs the difference at most this times
   !>   `d` over the step
   real(wp), parameter :: FD4_WEIGHT_SUM = 1.5_wp
   !> Step of the rigid paths of NumSA where spheres turn against each other, rad
   !>
   !> - the 13-ring with C2v: the nearest point that enters or leaves a
   !>   smoothing band lies between 7e-5 and 8e-5 along the first path, the
   !>   stencil reaches 4e-5
   real(wp), parameter :: FD_STEP_NUMSA_BAND = 2.0e-5_wp
   !> Smallest relative area of a DROP point with surface adjoint weights
   !>
   !> - relative to the largest area; a point near the weight cut may leave the
   !>   grid within the stencil
   real(wp), parameter :: drop_pin_floor = 1.0e-6_wp
   !> Relative tolerance of the end-to-end rotation checks (iSwiG)
   real(wp), parameter :: E2E_TOL = 1.0e-12_wp
   !> Relative tolerance of the end-to-end rotation checks of SvdW-DROP
   !>
   !> - its points are converged projections
   real(wp), parameter :: E2E_TOL_DROP = 1.0e-10_wp
   !> Projection level of the DROP end-to-end checks: local Newton
   integer, parameter :: e2e_drop_level = 1
   !> Rigid motions per molecule in the end-to-end rotation checks
   integer, parameter :: e2e_nmotion = 6

   !> Molecule, settings and host of a frame-gradient check
   type :: frame_fixture
      !> Cavity kind, iSwiG, DROP or NumSA
      integer :: kind = kind_iswig
      !> Name in messages
      character(len=:), allocatable :: name
      !> Reference structure
      type(structure_type) :: mol
      !> Symmetry settings
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      !> Fixed host charges beside the nuclei: positions (3, nel), charges (nel)
      real(wp), allocatable :: el_xyz(:, :), el_q(:)
      !> Whether the surface adjoint energy pins its weights to raw points
      !>
      !> Off with a group, whose ties may permute equivalent points
      logical :: pinned = .true.
      !> Whether both principal curvatures of a point share one weight
      !>
      !> - a node on a (near) threefold axis is an umbilic point: k1 and k2
      !>   cross there and neither has a derivative, their sum does
      logical :: shared_curvature = .false.
      !> Lebedev points per sphere
      integer :: num_leb = 110
      !> Switching cutoff, away from every switching value of the stencil (iSwiG)
      real(wp) :: cut_f = 1.0e-7_wp
      !> Central-difference step, bohr
      real(wp) :: step = FD_STEP
      !> Surviving raw points of the reference
      integer, allocatable :: numbering(:)
      !> Raw points that carry the surface adjoint energy's weights
      !>
      !> - iSwiG: every surviving point, the set held fixed by the cutoff
      !> - DROP: those whose area is above `drop_pin_floor` of the largest; a
      !>   point near the weight cut may leave the grid within the stencil
      integer, allocatable :: pinned_ids(:)
      !> Frame level of the reference
      integer :: level = 0
      !> Group used at the reference
      character(len=:), allocatable :: group_used
      !> Number of distinct rotations at the reference
      integer :: nrot = 0
      !> Rotation of every sphere at the reference (3, 3, nat)
      real(wp), allocatable :: rot(:, :, :)
   end type frame_fixture

   !> End-to-end state of a cavity built like a host builds it
   type :: e2e_state
      !> Number of points
      integer :: ngrid = 0
      !> Total area, bohr**2, and volume, bohr**3
      real(wp) :: area = 0.0_wp, volume = 0.0_wp
      !> Total-area gradient of the cavity (3, nat)
      real(wp), allocatable :: grad_area(:, :)
      !> Energy of the cavity
      !>
      !> - CPCM energy of the nuclear point charges
      !> - NumSA, which has no surface to polarize: an element-weighted sum of
      !>   the sphere areas
      real(wp) :: energy = 0.0_wp
      !> Its nuclear gradient, reverse and forward mode (3, nat)
      !>
      !> - NumSA has no reverse mode; its forward one runs through `dsdr`
      real(wp), allocatable :: grad_reverse(:, :), grad_forward(:, :)
      !> NumSA: area per sphere (nat), bohr**2
      real(wp), allocatable :: asph(:)
      !> NumSA: its gradients (3, nat, nat), `dsdr`
      real(wp), allocatable :: sphere_grad(:, :, :)
      !> Surface points (NumSA: quadrature nodes)
      type(surface_points_type) :: pts
      !> Group used
      character(len=:), allocatable :: group_used
      !> Its operations (3, 3, order); unallocated without a group
      real(wp), allocatable :: ops(:, :, :)
   end type e2e_state

contains

   !> Collect all cavity symmetry tests
   !>
   !> @param[out] testsuite Collection of tests
   subroutine collect_cavity_symmetry(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("settings_toml", test_settings_toml), &
                  new_unittest("settings_invalid", test_settings_invalid), &
                  new_unittest("settings_not_in_cavities", test_settings_not_in_cavities), &
                  new_unittest("defaults_lab_identity", test_defaults_lab), &
                  new_unittest("no_model_bitwise", test_no_model_bitwise), &
                  new_unittest("molecular_matches_orient", test_molecular_matches_orient), &
                  new_unittest("species_from_radii", test_species_from_radii), &
                  new_unittest("explicit_inexact_fails", test_explicit_inexact), &
                  new_unittest("lab_with_group_fails", test_lab_with_group), &
                  new_unittest("group_implies_molecular", test_group_implies_molecular), &
                  new_unittest("fingerprint_reuse", test_fingerprint_reuse), &
                  new_unittest("timer_sections", test_timer_sections), &
                  new_unittest("tolerance_fallback", test_tolerance_fallback), &
                  new_unittest("fields_numeric", test_fields), &
                  new_unittest("report", test_report), &
                  new_unittest("iswig_lab_bitwise", test_iswig_lab_bitwise), &
                  new_unittest("iswig_molecular_points", test_iswig_molecular_points), &
                  new_unittest("iswig_invariance", test_iswig_invariance), &
                  new_unittest("iswig_exact_symmetry", test_iswig_exact_symmetry), &
                  new_unittest("iswig_linear_anywhere", test_iswig_linear), &
                  new_unittest("iswig_frame_gradient_eigen", test_iswig_gradient_eigen), &
                  new_unittest("iswig_frame_gradient_plane", test_iswig_gradient_plane), &
                  new_unittest("iswig_frame_gradient_symmetric", test_iswig_gradient_symmetric), &
                  new_unittest("frame_helpers", test_frame_helpers), &
                  new_unittest("e2e_rotation_frame", test_e2e_frame), &
                  new_unittest("e2e_rotation_symmetric", test_e2e_symmetric), &
                  new_unittest("e2e_rotation_lab_sensitive", test_e2e_lab_sensitive), &
                  new_unittest("drop_anchors_rotated", test_drop_anchors), &
                  new_unittest("drop_frame_gradient_eigen", test_drop_gradient_eigen), &
                  new_unittest("drop_frame_gradient_plane", test_drop_gradient_plane), &
                  new_unittest("drop_frame_gradient_symmetric", test_drop_gradient_symmetric), &
                  new_unittest("drop_setup_timer", test_drop_setup_timer), &
                  new_unittest("drop_e2e_rotation_frame", test_drop_e2e_frame), &
                  new_unittest("drop_e2e_rotation_symmetric", test_drop_e2e_symmetric), &
                  new_unittest("drop_e2e_rotation_lab_sensitive", test_drop_e2e_lab_sensitive), &
                  new_unittest("numsa_rotated_nodes", test_numsa_nodes), &
                  new_unittest("numsa_frame_gradient_eigen", test_numsa_gradient_eigen), &
                  new_unittest("numsa_frame_gradient_plane", test_numsa_gradient_plane), &
                  new_unittest("numsa_frame_gradient_symmetric", test_numsa_gradient_symmetric), &
                  new_unittest("numsa_e2e_rotation_frame", test_numsa_e2e_frame), &
                  new_unittest("numsa_e2e_rotation_symmetric", test_numsa_e2e_symmetric), &
                  new_unittest("numsa_e2e_rotation_lab_sensitive", test_numsa_e2e_lab_sensitive), &
                  new_unittest("rigid_rotation_frame", test_rigid_rotation_frame), &
                  new_unittest("rigid_rotation_group", test_rigid_rotation_group), &
                  new_unittest("drop_rigid_rotation_group", test_drop_rigid_rotation_group), &
                  new_unittest("sphere_frame_and_atom", test_sphere_frame_and_atom), &
                  new_unittest("single_frame_tops", test_single_frame_tops), &
                  new_unittest("near_switch_warning", test_near_switch_warning), &
                  new_unittest("tied_frames_warning", test_tied_frames_warning), &
                  new_unittest("threshold_settings", test_threshold_settings), &
                  new_unittest("harmonic_and_linear_thresholds", test_harmonic_linear_thresholds), &
                  new_unittest("numsa_exact_symmetry", test_numsa_exact_symmetry), &
                  new_unittest("drop_exact_symmetry", test_drop_exact_symmetry), &
                  new_unittest("e2e_rotation_subgroups", test_e2e_subgroups), &
                  new_unittest("drop_e2e_rotation_subgroups", test_drop_e2e_subgroups), &
                  new_unittest("numsa_e2e_rotation_subgroups", test_numsa_e2e_subgroups), &
                  new_unittest("drop_rigid_rotation_more", test_drop_rigid_rotation_more), &
                  new_unittest("linear_cavities", test_linear_cavities), &
                  new_unittest("lab_angle_rigid_rotation", test_lab_angle_rigid_rotation), &
                  new_unittest("point_groups_iswig", test_point_groups_iswig), &
                  new_unittest("point_groups_numsa_drop", test_point_groups_numsa_drop), &
                  new_unittest("iswig_frame_gradient_sphere", test_iswig_gradient_sphere), &
                  new_unittest("numsa_frame_gradient_sphere", test_numsa_gradient_sphere), &
                  new_unittest("drop_frame_gradient_sphere", test_drop_gradient_sphere), &
                  new_unittest("single_atom_group", test_single_atom_group) &
                  ]

   end subroutine collect_cavity_symmetry

   !* ================================================================================= *!
   !*                                     Settings                                      *!
   !* ================================================================================= *!

   !> The top-level `symmetry` table reads, defaults and round-trips
   !>
   !> @param[out] error Error handle
   subroutine test_settings_toml(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_grid_model_symmetry_parameters_type) :: param, copy, defaults
      type(mctc_error), allocatable :: err
      character(len=*), parameter :: path = "moist-test-cavity-symmetry-settings.toml"
      character(len=*), parameter :: copy_path = "moist-test-cavity-symmetry-copy.json"
      character(len=line_len) :: lines(16)
      integer :: nline, unit
      logical :: ok

      ! Absent table: compiled defaults
      call read_fixture(param, path, [character(len=1) :: ""], err, ok)
      call check(error, ok .and. .not. allocated(err), "empty file not read")
      if (allocated(error)) return
      call check(error, trim(param%orientation) == trim(defaults%orientation) .and. &
         & trim(param%group) == trim(defaults%group) .and. &
         & abs(param%tolerance - defaults%tolerance) < 1.0e-20_wp, "defaults changed")
      if (allocated(error)) return
      call check(error, trim(param%orientation) == "auto" .and. trim(param%group) == "off" .and. &
         & param%effective_orientation() == "lab" .and. abs(param%tolerance - 1.0e-8_wp) < 1.0e-20_wp, &
         & "default is not auto, off, lab, 1e-8")
      if (allocated(error)) return

      ! A full table
      call read_fixture(param, path, [character(len=32) :: "[symmetry]", 'orientation = "molecular"', &
         & 'group = "D3h"', "tolerance = 1e-6"], err, ok)
      call check(error, ok, "fixture file")
      if (allocated(error)) return
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, trim(param%orientation) == "molecular" .and. trim(param%group) == "D3h" &
         & .and. abs(param%tolerance - 1.0e-6_wp) < 1.0e-20_wp, "table not read")
      if (allocated(error)) return

      ! Written and read back
      call param%write_file(copy_path, err)
      call check(error,.not. allocated(err), "write failed")
      if (allocated(error)) return
      call copy%read_file(copy_path, err)
      call check(error, delete_file(copy_path), "written file not removed")
      if (allocated(error)) return
      call check(error,.not. allocated(err), "written file not read")
      if (allocated(error)) return
      call check(error, trim(copy%orientation) == "molecular" .and. trim(copy%group) == "D3h" &
         & .and. abs(copy%tolerance - 1.0e-6_wp) < 1.0e-20_wp, "round trip")
      if (allocated(error)) return

      ! The table printout: one symmetry section with the normal keys
      open (newunit=unit, status="scratch", action="readwrite")
      call copy%print_parameters(err, unit)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error,.not. allocated(err) .and. any(index(lines, '"symmetry"') > 0) .and. &
         & any(index(lines, '"group"') > 0), "JSON printout lacks the symmetry table")
      if (allocated(error)) return

      ! The model takes the settings
      block
         type(grid_model_symmetry_type) :: model

         call new_grid_model_symmetry(model, err, param=copy)
         call check(error,.not. allocated(err) .and. trim(model%param%group) == "D3h" .and. &
            & .not. model%computed, "model does not take the settings")
      end block
   end subroutine test_settings_toml

   !> Unknown keys, an unknown orientation and a bad tolerance are refused
   !>
   !> @param[out] error Error handle
   subroutine test_settings_invalid(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_grid_model_symmetry_parameters_type) :: param
      type(mctc_error), allocatable :: err
      character(len=*), parameter :: path = "moist-test-cavity-symmetry-invalid.toml"
      logical :: ok

      call read_fixture(param, path, [character(len=32) :: "[symmetry]", 'grp = "auto"'], err, ok)
      call check(error, ok .and. allocated(err), "unknown key accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Unknown parameter: symmetry.grp") > 0, err%message)
      if (allocated(error)) return

      ! The keys live in the symmetry table, not at the top level
      call read_fixture(param, path, [character(len=32) :: 'group = "auto"'], err, ok)
      call check(error, ok .and. allocated(err), "top-level key accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Unknown parameter: group") > 0, err%message)
      if (allocated(error)) return

      call read_fixture(param, path, [character(len=32) :: "[symmetry]", 'orientation = "sideways"'], &
         & err, ok)
      call check(error, ok .and. allocated(err), "unknown orientation accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "orientation must be") > 0, err%message)
      if (allocated(error)) return

      call read_fixture(param, path, [character(len=32) :: "[symmetry]", 'orientation = ""'], err, ok)
      call check(error, ok .and. allocated(err), "blank orientation accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "orientation must be") > 0, err%message)
      if (allocated(error)) return

      call read_fixture(param, path, [character(len=32) :: "[symmetry]", "tolerance = -1.0"], err, ok)
      call check(error, ok .and. allocated(err), "negative tolerance accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "tolerance must be finite and positive") > 0, err%message)
      if (allocated(error)) return

      ! The constructor checks the same
      block
         type(grid_model_symmetry_type) :: model

         param%orientation = "auto"
         param%group = "off"
         param%tolerance = 0.0_wp
         call new_grid_model_symmetry(model, err, param=param)
         call check(error, allocated(err), "zero tolerance constructed")
         if (allocated(error)) return
         call check(error, index(err%message, "tolerance must be finite and positive") > 0, err%message)
      end block
   end subroutine test_settings_invalid

   !> Cavity parameter files have no `symmetry` table anymore
   !>
   !> @param[out] error Error handle
   subroutine test_settings_not_in_cavities(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      class(moist_model_parameters_type), allocatable :: param
      type(mctc_error), allocatable :: err
      character(len=*), parameter :: path = "moist-test-cavity-symmetry-cavity.toml"
      character(len=*), parameter :: dump = "moist-test-cavity-symmetry-cavity-dump.toml"
      character(len=line_len) :: lines(256)
      integer :: kind, nline, unit
      logical :: ok

      do kind = 1, size(kind_name)
         call new_param(kind, param)
         call read_fixture(param, path, [character(len=32) :: "[symmetry]", 'group = "auto"'], err, ok)
         call check(error, ok .and. allocated(err), trim(kind_name(kind))//": symmetry table accepted")
         if (allocated(error)) return
         call check(error, index(err%message, "Unknown parameter: symmetry.group") > 0, err%message)
         if (allocated(error)) return

         call new_param(kind, param)
         call param%init_defaults()
         call param%write_file(dump, err)
         call check(error,.not. allocated(err), trim(kind_name(kind))//": dump failed")
         if (allocated(error)) return
         open (newunit=unit, file=dump, status="old", action="read")
         call read_printout(unit, lines, nline)
         close (unit)
         call check(error, delete_file(dump), "dump not removed")
         if (allocated(error)) return
         call check(error,.not. any(index(lines(1:nline), "symmetry") > 0), &
            & trim(kind_name(kind))//": parameters still write a symmetry table")
         if (allocated(error)) return
      end do
   end subroutine test_settings_not_in_cavities

   !* ================================================================================= *!
   !*                                     Defaults                                      *!
   !* ================================================================================= *!

   !> Lab and off: exact identity, bitwise unit-grid copies, no frame gradient
   !>
   !> @param[out] error Error handle
   subroutine test_defaults_lab(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      type(moist_math_grid_angular_type) :: leb
      real(wp), allocatable :: rotated(:, :, :), gradient(:, :), torque(:, :)
      real(wp) :: rot(3, 3)
      integer :: kind, isph, ipt, ileb, iat

      call new_context(ctx, nthreads=0, verbosity=0)
      call structure_of("water", mol)
      call new_lebedev_grid(leb, err, npts=110)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if

      do kind = 1, size(kind_name)
         call updated_cavity(kind, mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, trim(kind_name(kind))//": "//err%message)
            return
         end if
         call check(error, allocated(cav%grid_model) .and. cav%plain_grid(), &
            & trim(kind_name(kind))//": default model missing or not in the lab frame")
         if (allocated(error)) return
         associate (sym => symmetry_of(cav))
            call check(error, sym%computed .and. sym%is_lab() .and. model_nrot(sym) == 1, &
               & trim(kind_name(kind))//": defaults are not one lab rotation")
            if (allocated(error)) return
            call check(error, sym%orient%group_used == "C1" .and. len(sym%orient%group_detected) == 0, &
               & trim(kind_name(kind))//": group off ran a detection")
            if (allocated(error)) return
            do isph = 1, mol%nat
               rot = model_rotation(sym, isph)
               call check(error, same_bits(rot, identity3) .and. model_rot_index(sym, isph) == 1, &
                  & trim(kind_name(kind))//": a sphere rotation is not the exact identity")
               if (allocated(error)) return
            end do
            call model_rotate_grid(sym, leb%points, rotated)
            call check(error, all(shape(rotated) == [3, 110, 1]), trim(kind_name(kind))//": rotated grid shape")
            if (allocated(error)) return
            call check(error, same_bits(rotated(:, :, 1), leb%points), &
               & trim(kind_name(kind))//": lab grid is not a bitwise copy")
            if (allocated(error)) return
            allocate (gradient(3, mol%nat), source=0.5_wp)
            allocate (torque(3, mol%nat), source=1.0_wp)
            call model_add_torque(sym, torque, gradient)
            call check(error, same_bits(gradient, spread(spread(0.5_wp, 1, 3), 2, mol%nat)), &
               & trim(kind_name(kind))//": lab mode has a frame gradient")
            if (allocated(error)) return
            deallocate (gradient, torque)
         end associate

         ! iSwiG points are c + r R u with the cached unit grid, bitwise
         select type (cav)
         type is (cavity_type_iswig)
            associate (gm_ => symmetry_of(cav))
               call model_rotate_grid(gm_, cav%ang_grid, rotated)
            end associate
            do ipt = 1, cav%ngrid
               iat = cav%owner(ipt)
               ileb = cav%numbering(ipt) - (iat - 1)*cav%param%num_leb
               associate (gm_ => symmetry_of(cav))
                  call check(error, same_bits(cav%xyz(:, ipt:ipt), spread(cav%sphxyz(:, iat) + &
                     & cav%radii(iat)*rotated(:, ileb, model_rot_index(gm_, iat)), 2, 1)), &
                     & "iSwiG: a grid point is not c + r R u bitwise")
               end associate
               if (allocated(error)) return
            end do
         end select
      end do
   end subroutine test_defaults_lab

   !> No model passed: no symmetry state, results bitwise those of a lab model
   !>
   !> - iSwiG points, areas, switching, normals, totals and gradients; NumSA
   !>   and DROP totals and points; NumSA sphere areas and every area
   !>   gradient; DROP anchors, every forward array and a reverse contraction
   !>   of every adjoint channel
   !> - no `symmetry_*` field and no report without a model
   !>
   !> @param[out] error Error handle
   subroutine test_no_model_bitwise(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: plain, lab
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      type(field_query_type) :: query
      character(len=line_len) :: lines(64)
      integer :: kind, unit, nline

      call new_context(ctx, nthreads=0, verbosity=0)
      call get_structure(mol, "MB16-43", "01")
      settings%orientation = "lab"
      do kind = 1, size(kind_name)
         call updated_cavity(kind, mol, settings, ctx, plain, err, with_model=.false.)
         if (.not. allocated(err)) call updated_cavity(kind, mol, settings, ctx, lab, err)
         if (allocated(err)) then
            call test_failed(error, trim(kind_name(kind))//": "//err%message)
            return
         end if
         call check(error,.not. allocated(plain%grid_model) .and. plain%plain_grid(), &
            & trim(kind_name(kind))//": a cavity without a model has symmetry state")
         if (allocated(error)) return
         call check(error, plain%ngrid == lab%ngrid, trim(kind_name(kind))//": point counts differ")
         if (allocated(error)) return
         call check(error, same_bits(reshape([plain%total_area, plain%total_volume], [2, 1]), &
            & reshape([lab%total_area, lab%total_volume], [2, 1])), &
            & trim(kind_name(kind))//": area or volume not bitwise")
         if (allocated(error)) return
         if (plain%ngrid > 0 .and. allocated(plain%xyz)) then
            call check(error, same_bits(plain%xyz(:, 1:plain%ngrid), lab%xyz(:, 1:lab%ngrid)), &
               & trim(kind_name(kind))//": points not bitwise")
            if (allocated(error)) return
         end if

         select type (plain)
         type is (cavity_type_iswig)
            select type (lab)
            type is (cavity_type_iswig)
               call check(error, same_bits(reshape(plain%a, [1, plain%ngrid]), reshape(lab%a, [1, lab%ngrid])) &
                  & .and. same_bits(reshape(plain%f, [1, plain%ngrid]), reshape(lab%f, [1, lab%ngrid])) &
                  & .and. same_bits(plain%normal0, lab%normal0), "iSwiG: areas, switching or normals not bitwise")
               if (allocated(error)) return
               call plain%get_gradient(err)
               if (.not. allocated(err)) call lab%get_gradient(err)
               if (allocated(err)) then
                  call test_failed(error, err%message)
                  return
               end if
               call check(error, same_bits(plain%area_grad, lab%area_grad) .and. &
                  & same_bits(plain%volume_grad, lab%volume_grad), "iSwiG: gradients not bitwise")
               if (allocated(error)) return
            end select
         type is (cavity_type_drop)
            select type (lab)
            type is (cavity_type_drop)
               call check_drop_bitwise(plain, lab, error)
               if (allocated(error)) return
            end select
         type is (cavity_type_numsa)
            select type (lab)
            type is (cavity_type_numsa)
               call plain%get_gradient(err)
               if (.not. allocated(err)) call lab%get_gradient(err)
               if (allocated(err)) then
                  call test_failed(error, err%message)
                  return
               end if
               call check(error, same_bits(reshape(plain%asph, [1, plain%nsph]), reshape(lab%asph, [1, lab%nsph])) &
                  & .and. same_bits(plain%area_grad, lab%area_grad) .and. same_bits(reshape(plain%dsdr, &
                  & [3, plain%nsph**2]), reshape(lab%dsdr, [3, lab%nsph**2])), &
                  & "NumSA: sphere areas or area gradients not bitwise")
               if (allocated(error)) return
               call check(error,.not. ((.not. lab%grid%plain) .or. allocated(lab%frame_torque)), &
                  & "NumSA: the lab orientation built rotated nodes or a torque")
               if (allocated(error)) return
            end select
         end select

         call query%enumerate()
         call plain%list_fields(query)
         call check(error, query%index_of("symmetry_orientation") == 0, &
            & trim(kind_name(kind))//": symmetry fields without a model")
         if (allocated(error)) return
         open (newunit=unit, status="scratch", action="readwrite")
         call plain%print_radii(unit, 3)
         call read_printout(unit, lines, nline)
         close (unit)
         call check(error, any(index(lines(1:nline), "Used radii") > 0) .and. &
            & .not. any(index(lines(1:nline), "Symmetry:") > 0), &
            & trim(kind_name(kind))//": radii printout without a model has a symmetry report")
         if (allocated(error)) return
      end do
   end subroutine test_no_model_bitwise

   !> Two DROP cavities of one structure agree bit for bit
   !>
   !> Anchors, the forward Jacobians and totals, and the reverse contraction of
   !> a surface adjoint on every channel
   !>
   !> @param[in,out] plain Cavity without a symmetry model
   !> @param[in,out] lab   Cavity with a lab-frame model
   !> @param[out]    error Test failure
   subroutine check_drop_bitwise(plain, lab, error)
      !> Cavity without a symmetry model
      type(cavity_type_drop), intent(inout) :: plain
      !> Cavity with a lab-frame model
      type(cavity_type_drop), intent(inout) :: lab
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: grad_plain(:, :), grad_lab(:, :)
      integer :: ip, n

      call check(error, same_bits(plain%anchor%xyz, lab%anchor%xyz) .and. .not. (.not. lab%grid%plain), &
         & "DROP: anchors not bitwise")
      if (allocated(error)) return
      call plain%get_gradient(err)
      if (.not. allocated(err)) call lab%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      n = plain%ngrid
      call check(error, same_bits(reshape(plain%xyz1_rA, [9, plain%nsph*n]), reshape(lab%xyz1_rA, [9, lab%nsph*n])) &
         & .and. same_bits(reshape(plain%f1_rA, [3, plain%nsph*n]), reshape(lab%f1_rA, [3, lab%nsph*n])) &
         & .and. same_bits(reshape(plain%xi1_rA, [3, plain%nsph*n]), reshape(lab%xi1_rA, [3, lab%nsph*n])) &
         & .and. same_bits(reshape(plain%v1_rA, [3, plain%nsph*n]), reshape(lab%v1_rA, [3, lab%nsph*n])) &
         & .and. same_bits(plain%A_tot1_rA, lab%A_tot1_rA) .and. same_bits(plain%V_tot1_rA, lab%V_tot1_rA), &
         & "DROP: forward gradients not bitwise")
      if (allocated(error)) return

      call acc%init(n)
      do ip = 1, n
         acc%w_xi(ip) = sin(0.3_wp*real(ip, wp))
         acc%w_f(ip) = cos(0.7_wp*real(ip, wp))
         acc%w_a(ip) = sin(1.1_wp*real(ip, wp) + 0.2_wp)
         acc%w_w(ip) = cos(0.5_wp*real(ip, wp) + 0.4_wp)
         acc%w_xyz(:, ip) = sin(0.9_wp*real(ip, wp) + [0.1_wp, 1.2_wp, 2.3_wp])
         acc%w_n(:, ip) = cos(0.4_wp*real(ip, wp) + [0.5_wp, 1.6_wp, 2.7_wp])
      end do
      allocate (grad_plain(3, plain%nsph), grad_lab(3, lab%nsph), source=0.0_wp)
      call plain%get_surface_gradient(acc, grad_plain, err)
      if (.not. allocated(err)) call lab%get_surface_gradient(acc, grad_lab, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, same_bits(grad_plain, grad_lab), "DROP: reverse gradient not bitwise")
   end subroutine check_drop_bitwise

   !* ================================================================================= *!
   !*                                   Orientations                                    *!
   !* ================================================================================= *!

   !> Molecular orientation and auto group on a cavity: the math layer's result
   !>
   !> @param[out] error Error handle
   subroutine test_molecular_matches_orient(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(sym_orient_type) :: ref
      type(mctc_error), allocatable :: err
      type(moist_math_grid_angular_type) :: leb
      real(wp), allocatable :: rotated(:, :, :), torque(:, :), gradient(:, :), expected(:, :)
      integer :: kind, isph, irot, imode
      character(len=*), parameter :: groups(2) = [character(len=4) :: "auto", "off"]

      call new_context(ctx, nthreads=0, verbosity=0)
      call structure_of("ammonia", mol)
      call new_lebedev_grid(leb, err, npts=26)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      allocate (torque(3, mol%nat))
      do isph = 1, mol%nat
         torque(:, isph) = [0.3_wp*real(isph, wp), -0.2_wp, 0.1_wp + 0.05_wp*real(isph, wp)]
      end do

      do imode = 1, size(groups)
         settings%orientation = "molecular"
         settings%group = groups(imode)
         call new_sym_orient(ref, mol%xyz, mol%num(mol%id), orient_molecular, trim(groups(imode)), TOL, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         do kind = 1, size(kind_name)
            call updated_cavity(kind, mol, settings, ctx, cav, err)
            if (allocated(err)) then
               call test_failed(error, trim(kind_name(kind))//": "//err%message)
               return
            end if
            associate (sym => symmetry_of(cav))
               call check(error,.not. sym%is_lab() .and. sym%orient%orientation == orient_molecular, &
                  & trim(kind_name(kind))//": molecular orientation not used")
               if (allocated(error)) return
               call check(error, sym%orient%group_used == ref%group_used .and. sym%orient%nrot == ref%nrot .and. &
                  & all(sym%orient%rot_of_atom == ref%rot_of_atom), &
                  & trim(kind_name(kind))//" ("//trim(groups(imode))//"): group or rotation map differs")
               if (allocated(error)) return
               call check(error, maxval(abs(sym%orient%rotation - ref%rotation)) < 1.0e-14_wp, &
                  & trim(kind_name(kind))//": rotations differ from the math layer")
               if (allocated(error)) return
               if (imode == 1) then
                  call check(error, sym%orient%group_used == "C3v", trim(kind_name(kind))//": ammonia is not C3v, got "// &
                     & sym%orient%group_used)
               else
                  ! Ammonia's three harmonic candidates give one grid, its
                  ! threefold axis on the grid diagonal
                  call check(error, sym%orient%nrot == 1 .and. &
                     & maxval(abs(sym%orient%rotation(:, :, 1) - sym%orient%frame%axes)) < 1.0e-14_wp, &
                     & trim(kind_name(kind))//": molecular off is not the frame")
               end if
               if (allocated(error)) return

               ! Accessors agree with the stored data
               call model_rotate_grid(sym, leb%points, rotated)
               call check(error, size(rotated, 3), sym%orient%nrot, trim(kind_name(kind))//": one copy per rotation")
               if (allocated(error)) return
               do irot = 1, sym%orient%nrot
                  call check(error, maxval(abs(rotated(:, :, irot) - &
                     & matmul(sym%orient%rotation(:, :, irot), leb%points))) < 1.0e-15_wp, &
                     & trim(kind_name(kind))//": rotated copy is not R u")
                  if (allocated(error)) return
               end do
               do isph = 1, mol%nat
                  call check(error, same_bits(model_rotation(sym, isph), &
                     & sym%orient%rotation(:, :, sym%orient%rot_of_atom(isph))) .and. &
                     & model_rot_index(sym, isph) == sym%orient%rot_of_atom(isph), trim(kind_name(kind))//": sphere accessors")
                  if (allocated(error)) return
               end do

               ! Frame gradient: the math layer's VJP, added
               allocate (gradient(3, mol%nat), source=1.0_wp)
               allocate (expected(3, mol%nat))
               call sym_orient_vjp(ref, mol%xyz, torque, expected)
               call model_add_torque(sym, torque, gradient)
               call check(error, maxval(abs(gradient - 1.0_wp - expected)) < 1.0e-12_wp, &
                  & trim(kind_name(kind))//": frame gradient differs from the math layer")
               if (allocated(error)) return
               call check(error, maxval(abs(expected)) > 0.0_wp, trim(kind_name(kind))//": vacuous frame gradient")
               if (allocated(error)) return
               deallocate (gradient, expected)
            end associate
         end do
      end do
   end subroutine test_molecular_matches_orient

   !> Same element, different radius: not equivalent; roundoff does not split
   !>
   !> @param[out] error Error handle
   subroutine test_species_from_radii(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      class(radius_type), allocatable :: radius_model
      type(mctc_error), allocatable :: err
      real(wp), parameter :: r_o = 3.2_wp, r_h = 2.3_wp
      real(wp) :: radii(3, 3)
      character(len=3) :: want(3)
      integer :: kind, icase

      call new_context(ctx, nthreads=0, verbosity=0)
      call structure_of("water", mol)
      radii(:, 1) = [r_o, r_h, r_h]
      want(1) = "C2v"
      radii(:, 2) = [r_o, r_h, r_h*(1.0_wp + 1.0e-13_wp)]
      want(2) = "C2v"
      radii(:, 3) = [r_o, r_h, 1.1_wp*r_h]
      want(3) = "Cs"
      settings%orientation = "molecular"
      settings%group = "auto"

      do icase = 1, size(want)
         call new_radii_custom_atoms(radii(:, icase), radius_model, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         do kind = 1, size(kind_name)
            call updated_cavity(kind, mol, settings, ctx, cav, err, radius_model)
            if (allocated(err)) then
               call test_failed(error, trim(kind_name(kind))//": "//err%message)
               return
            end if
            associate (gm_ => symmetry_of(cav))
               call check(error, gm_%orient%group_used == trim(want(icase)), &
                  & trim(kind_name(kind))//": expected "//trim(want(icase))//", got "//gm_%orient%group_used)
            end associate
            if (allocated(error)) return
         end do
      end do
   end subroutine test_species_from_radii

   !> An explicit group the geometry misses fails the update
   !>
   !> - auto lowers the group instead; a bad orientation fails the constructor
   !>
   !> @param[out] error Error handle
   subroutine test_explicit_inexact(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      integer :: kind, unit

      ! Auto on a distortion warns; keep it off the test output
      open (newunit=unit, status="scratch", action="readwrite")
      call new_context(ctx, nthreads=0, verbosity=0, unit=unit)
      call structure_of("water_distorted", mol)
      do kind = 1, size(kind_name)
         settings%orientation = "molecular"
         settings%group = "C2v"
         call updated_cavity(kind, mol, settings, ctx, cav, err)
         call check(error, allocated(err), trim(kind_name(kind))//": inexact C2v accepted")
         if (allocated(error)) exit
         call check(error, index(err%message, "deviates from C2v") > 0, err%message)
         if (allocated(error)) exit
         associate (gm_ => symmetry_of(cav))
            call check(error,.not. cav%is_updated() .and. .not. gm_%computed, &
               & trim(kind_name(kind))//": failed update left the cavity updated")
         end associate
         if (allocated(error)) exit

         ! Auto lowers the group instead
         settings%group = "auto"
         call updated_cavity(kind, mol, settings, ctx, cav, err)
         call check(error,.not. allocated(err), trim(kind_name(kind))//": auto failed on a distortion")
         if (allocated(error)) exit
         associate (gm_ => symmetry_of(cav))
            call check(error, gm_%orient%group_used == "Cs", trim(kind_name(kind))//": distorted water not Cs")
         end associate
         if (allocated(error)) exit

         settings%group = "off"
         settings%orientation = "sideways"
         call updated_cavity(kind, mol, settings, ctx, cav, err)
         call check(error, allocated(err), trim(kind_name(kind))//": unknown orientation accepted")
         if (allocated(error)) exit
         call check(error, index(err%message, "orientation must be") > 0, err%message)
         if (allocated(error)) exit
      end do
      close (unit)
   end subroutine test_explicit_inexact

   !> A group other than off with an explicit lab orientation is refused
   !>
   !> - by the constructor, by an update of a model changed afterwards, and in
   !>   a settings file
   !>
   !> @param[out] error Error handle
   subroutine test_lab_with_group(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings, param
      type(grid_model_symmetry_type) :: model
      type(mctc_error), allocatable :: err
      character(len=*), parameter :: path = "moist-test-cavity-symmetry-lab-group.toml"
      character(len=*), parameter :: groups(3) = [character(len=4) :: "auto", "C2v", "C1"]
      integer :: kind, igroup
      logical :: ok

      call new_context(ctx, nthreads=0, verbosity=0)
      call structure_of("water", mol)
      do igroup = 1, size(groups)
         settings%orientation = "lab"
         settings%group = groups(igroup)
         do kind = 1, size(kind_name)
            call updated_cavity(kind, mol, settings, ctx, cav, err)
            call check(error, allocated(err), trim(kind_name(kind))//": lab with "//trim(groups(igroup))// &
               & " accepted")
            if (allocated(error)) return
            call check(error, index(err%message, 'symmetry.group = "'//trim(groups(igroup))// &
               & '" requires symmetry.orientation = "molecular"') > 0, err%message)
            if (allocated(error)) return
         end do

         ! A model changed after construction is checked at update
         call new_grid_model_symmetry(model, err)
         model%param = settings
         call model%update(mol, [3.0_wp, 2.0_wp, 2.0_wp], err)
         call check(error, allocated(err) .and. .not. model%computed, "update of lab with a group accepted")
         if (allocated(error)) return
      end do

      ! The same in a file
      call read_fixture(param, path, [character(len=32) :: "[symmetry]", 'orientation = "lab"', &
         & 'group = "auto"'], err, ok)
      call check(error, ok .and. allocated(err), "lab with auto read")
      if (allocated(error)) return
      call check(error, index(err%message, 'symmetry.group = "auto" requires symmetry.orientation = '// &
         & '"molecular"') > 0, err%message)
   end subroutine test_lab_with_group

   !> Auto orientation is molecular with a group and lab without
   !>
   !> - files preserve "auto" across a round trip; the report marks it as implied
   !>
   !> @param[out] error Error handle
   subroutine test_group_implies_molecular(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings, param, read_back
      type(mctc_error), allocatable :: err
      character(len=line_len) :: lines(64)
      character(len=*), parameter :: path = "moist-test-cavity-symmetry-implied.toml"
      character(len=*), parameter :: dump = "moist-test-cavity-symmetry-implied-dump.toml"
      integer :: kind, nline, unit
      logical :: ok

      call new_context(ctx, nthreads=0, verbosity=0)
      call structure_of("water", mol)
      do kind = 1, size(kind_name)
         ! Auto alone: molecular, marked as implied
         settings%orientation = "auto"
         settings%group = "auto"
         call updated_cavity(kind, mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, trim(kind_name(kind))//": "//err%message)
            return
         end if
         associate (gm_ => symmetry_of(cav))
            call check(error, gm_%orient%orientation == orient_molecular .and. &
               & gm_%orient%group_used == "C2v", trim(kind_name(kind))//": auto alone is not molecular C2v")
         end associate
         if (allocated(error)) return
         open (newunit=unit, status="scratch", action="readwrite")
         call cav%print_grid_model(unit, 1)
         call read_printout(unit, lines, nline)
         close (unit)
         call check(error, printed_entry(lines, "Orientation", "molecular (implied by group)"), &
            & trim(kind_name(kind))//": implied orientation not marked")
         if (allocated(error)) return

         ! Off with the default: lab, not marked
         settings%group = "off"
         call updated_cavity(kind, mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, trim(kind_name(kind))//": "//err%message)
            return
         end if
         associate (gm_ => symmetry_of(cav))
            call check(error, gm_%is_lab() .and. gm_%orient%orientation == orient_lab, &
               & trim(kind_name(kind))//": off with the default is not lab")
         end associate
         if (allocated(error)) return
         open (newunit=unit, status="scratch", action="readwrite")
         call cav%print_grid_model(unit, 1)
         call read_printout(unit, lines, nline)
         close (unit)
         call check(error, printed_entry(lines, "Orientation", "lab") .and. &
            & .not. any(index(lines, "implied") > 0), trim(kind_name(kind))//": off printout")
         if (allocated(error)) return
      end do

      ! A file with a group only: orientation defaults to auto and round-trips as auto
      call read_fixture(param, path, [character(len=32) :: "[symmetry]", 'group = "auto"'], err, ok)
      call check(error, ok .and. .not. allocated(err), "group alone refused")
      if (allocated(error)) return
      call check(error, trim(param%orientation) == "auto" .and. &
         & param%effective_orientation() == "molecular" .and. param%orientation_implied(), &
         & "group alone not implied")
      if (allocated(error)) return
      call param%write_file(dump, err)
      call check(error,.not. allocated(err), "dump failed")
      if (allocated(error)) return
      open (newunit=unit, file=dump, status="old", action="read")
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, delete_file(dump), "dump not removed")
      if (allocated(error)) return
      call check(error, any(index(lines, "[symmetry]") > 0) .and. &
         & any(index(lines, 'orientation = "auto"') > 0), &
         & "dump is no symmetry table with auto orientation")
      if (allocated(error)) return
      call read_fixture(read_back, path, lines(1:nline), err, ok)
      call check(error, ok .and. .not. allocated(err), "dump not read back")
      if (allocated(error)) return
      call check(error, trim(read_back%orientation) == "auto" .and. &
         & read_back%effective_orientation() == "molecular" .and. read_back%orientation_implied() .and. &
         & trim(read_back%group) == "auto", "dump round trip")
      if (allocated(error)) return

      ! Explicit auto without a group resolves to lab
      call read_fixture(param, path, [character(len=32) :: "[symmetry]", 'orientation = "auto"', &
         & 'group = "off"'], err, ok)
      call check(error, ok .and. .not. allocated(err), "explicit auto with off refused")
      if (allocated(error)) return
      call check(error, trim(param%orientation) == "auto" .and. &
         & param%effective_orientation() == "lab" .and. .not. param%orientation_implied(), &
         & "explicit auto with off is not lab")
   end subroutine test_group_implies_molecular

   !* ================================================================================= *!
   !*                              Ownership and lifetime                               *!
   !* ================================================================================= *!

   !> One detection per geometry: the host's update is reused by the cavity
   !>
   !> - host update, then a cavity constructed with the model: the cavity's
   !>   update finds it current, the timer shows one detection
   !> - the same structure again: still one; a moved structure: two, the
   !>   second under the cavity's "Setup"; reused again after that
   !> - a changed radius or setting is no longer current
   !>
   !> @param[out] error Error handle
   subroutine test_fingerprint_reuse(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, moved
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(grid_model_symmetry_type) :: model
      class(radius_type), allocatable :: radius_model
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: radii(:)
      integer :: kind, istep
      !> Structure of every update and the detections expected after it
      integer, parameter :: nstep = 4
      integer, parameter :: which(nstep) = [1, 1, 2, 2], expect(nstep) = [1, 1, 2, 2]

      call structure_of("ammonia", mol)
      call rigid_motion(mol, test_rotation(1), [0.4_wp, -1.1_wp, 0.7_wp], test_permutation(mol%nat, 1), moved)
      settings%group = "auto"
      do kind = 1, size(kind_name)
         call new_context(ctx, nthreads=0, verbosity=0)
         call host_model(mol, settings, ctx, model, radius_model, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         call check(error, model%computed .and. model%is_current(mol, radius_model%f0), &
            & "host model not current after its update")
         if (allocated(error)) return
         call check(error, detections(ctx), 1, trim(kind_name(kind))//": host update not one detection")
         if (allocated(error)) return

         call new_cavity_of(kind, radius_model, ctx, cav, err, model)
         do istep = 1, nstep
            if (.not. allocated(err)) then
               if (which(istep) == 1) then
                  call cav%update(mol, err)
               else
                  call cav%update(moved, err)
               end if
            end if
            if (allocated(err)) then
               call test_failed(error, trim(kind_name(kind))//": "//err%message)
               return
            end if
            call check(error, detections(ctx), expect(istep), trim(kind_name(kind))// &
               & ": unexpected number of detections after an update")
            if (allocated(error)) return
         end do
         associate (gm_ => symmetry_of(cav))
            call check(error, gm_%orient%group_used == "C3v", trim(kind_name(kind))//": moved not C3v")
         end associate
         if (allocated(error)) return
         call check(error, timer_path(ctx, "Setup/Symmetry/Detection"), &
            & trim(kind_name(kind))//": cavity's own computation not under Setup")
         if (allocated(error)) return
         call ctx%delete()
      end do

      ! Radii and settings are part of the fingerprint
      radii = radius_model%f0
      call check(error, model%is_current(mol, radii), "model not current")
      if (allocated(error)) return
      radii(2) = radii(2)*(1.0_wp + epsilon(1.0_wp))
      call check(error,.not. model%is_current(mol, radii), "changed radius still current")
      if (allocated(error)) return
      call check(error,.not. model%is_current(moved, radius_model%f0), "moved structure still current")
      if (allocated(error)) return
      model%param%tolerance = 2.0e-8_wp
      call check(error,.not. model%is_current(mol, radius_model%f0), "changed tolerance still current")
      if (allocated(error)) return
      model%param%tolerance = settings%tolerance
      model%param%group = "C3v"
      call check(error,.not. model%is_current(mol, radius_model%f0), "changed group still current")
   end subroutine test_fingerprint_reuse

   !> Timer layout of the symmetry computation
   !>
   !> - host update with a group: "Symmetry" with "Molecular frame",
   !>   "Detection", "Atom matching" and "Orientation" under the open node
   !> - molecular orientation without a group: the frame only
   !> - lab orientation without a group: no section, on the host or the cavity
   !>
   !> @param[out] error Error handle
   subroutine test_timer_sections(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(grid_model_symmetry_type) :: model
      class(radius_type), allocatable :: radius_model
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      character(len=*), parameter :: phases(4) = [character(len=15) :: "Molecular frame", "Detection", &
         & "Atom matching", "Orientation"]
      integer :: k

      call structure_of("benzene", mol)

      call new_context(ctx, nthreads=0, verbosity=0)
      call ctx%timer%start("host")
      settings%group = "auto"
      call host_model(mol, settings, ctx, model, radius_model, err)
      call ctx%timer%stop("host")
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, timer_path(ctx, "host/Symmetry"), "no Symmetry section under the open node")
      if (allocated(error)) return
      do k = 1, size(phases)
         call check(error, timer_path(ctx, "host/Symmetry/"//trim(phases(k))), "Symmetry lacks "//trim(phases(k)))
         if (allocated(error)) return
      end do
      call check(error, ctx%timer%current_depth(), 0, "symmetry update left a timer open")
      if (allocated(error)) return
      call ctx%delete()

      ! Molecular orientation, no group: the frame only
      call new_context(ctx, nthreads=0, verbosity=0)
      settings%group = "off"
      settings%orientation = "molecular"
      call host_model(mol, settings, ctx, model, radius_model, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, timer_path(ctx, "Symmetry/Molecular frame") .and. &
         & .not. timer_has(ctx, "Detection"), "molecular off is not the frame alone")
      if (allocated(error)) return
      call ctx%delete()

      ! Lab orientation, no group: nothing on the host, nor on a cavity that
      ! computes the model itself
      call new_context(ctx, nthreads=0, verbosity=0)
      settings%orientation = "lab"
      call host_model(mol, settings, ctx, model, radius_model, err)
      if (.not. allocated(err)) then
         call model%clear()
         call new_cavity_of(kind_iswig, radius_model, ctx, cav, err, model)
      end if
      if (.not. allocated(err)) call cav%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      associate (gm_ => symmetry_of(cav))
         call check(error, gm_%computed .and. .not. timer_has(ctx, "Symmetry") .and. &
            & .not. timer_has(ctx, "Setup"), "lab orientation without a group has a timer section")
      end associate
      if (allocated(error)) return
      call ctx%delete()
   end subroutine test_timer_sections

   !> A detected group the geometry misses only by the tolerance
   !>
   !> - noisy benzene: every operation but the identity fails, C1 is used;
   !>   the report shows the detected group with its deviation and the reason
   !>   with both numbers; the context warns
   !> - a larger tolerance uses the group, no warning
   !> - in-plane distorted water: Cs remains, the reason and warning name C2v
   !>
   !> @param[out] error Error handle
   subroutine test_tolerance_fallback(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(grid_model_symmetry_type) :: model
      class(radius_type), allocatable :: radius_model
      type(mctc_error), allocatable :: err
      character(len=line_len) :: lines(64), out(64)
      character(len=:), allocatable :: dev
      integer :: unit, nline, nout, ctx_unit

      call structure_of("benzene_noisy", mol)
      settings%group = "auto"
      open (newunit=ctx_unit, status="scratch", action="readwrite")
      call new_context(ctx, nthreads=0, verbosity=0, unit=ctx_unit)
      call host_model(mol, settings, ctx, model, radius_model, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         close (ctx_unit)
         return
      end if
      call read_printout(ctx_unit, out, nout)
      close (ctx_unit)
      call ctx%delete()

      associate (orient => model%orient)
         call check(error, orient%group_detected == "D6h" .and. orient%group_used == "C1" .and. &
            & orient%nrot == 1 .and. size(orient%ops, 3) == 1 .and. orient%norbit == mol%nat, &
            & "noisy benzene: expected D6h detected, C1 used, got "//orient%group_detected//", "// &
            & orient%group_used)
         if (allocated(error)) return
         call check(error, orient%deviation_detected > TOL .and. orient%deviation_detected < 1.0e-4_wp, &
            & "noisy benzene: detected deviation out of range")
         if (allocated(error)) return
         call check(error, model%missed_tolerance(), "noisy benzene: tolerance miss not flagged")
         if (allocated(error)) return
         dev = to_string(orient%deviation_detected, "(es9.2)")
      end associate

      open (newunit=unit, status="scratch", action="readwrite")
      call model%print(unit, 1)
      call read_printout(unit, lines, nline)
      close (unit)
      call check(error, printed_entry(lines, "Group detected", "D6h (deviation "//dev//" bohr)") .and. &
         & printed_entry(lines, "Group used", "C1") .and. .not. printed_entry(lines, "Group used", "order") &
         & .and. printed_entry(lines, "Reason", "geometry matches D6h only to "//dev// &
         & " bohr, above the tolerance 1.00E-08 bohr; no symmetry used") .and. &
         & printed_entry(lines, "Tolerance", "1.00E-08 bohr"), "noisy benzene report lacks a key line")
      if (allocated(error)) return
      call check(error, nout >= 1 .and. any(index(out(1:nout), "[Warning] Symmetry: the geometry matches D6h only to "// &
         & dev//" bohr") > 0 .and. index(out(1:nout), "symmetrize the input geometry or raise symmetry.tolerance") > 0), &
         & "noisy benzene: no warning")
      if (allocated(error)) return

      ! A tolerance above the noise: D6h itself, no warning
      settings%tolerance = 1.0e-3_wp
      open (newunit=ctx_unit, status="scratch", action="readwrite")
      call new_context(ctx, nthreads=0, verbosity=0, unit=ctx_unit)
      call host_model(mol, settings, ctx, model, radius_model, err)
      call read_printout(ctx_unit, out, nout)
      close (ctx_unit)
      call ctx%delete()
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, model%orient%group_used == "D6h" .and. len(model%orient%reason) == 0 .and. &
         & .not. model%missed_tolerance() .and. .not. any(index(out(1:nout), "only to") > 0), &
         & "raised tolerance: expected D6h without a tolerance warning")
      if (allocated(error)) return

      ! A subgroup survives the tolerance: reason and warning name the miss
      call structure_of("water_distorted", mol)
      settings%tolerance = TOL
      open (newunit=ctx_unit, status="scratch", action="readwrite")
      call new_context(ctx, nthreads=0, verbosity=0, unit=ctx_unit)
      call host_model(mol, settings, ctx, model, radius_model, err)
      call read_printout(ctx_unit, out, nout)
      close (ctx_unit)
      call ctx%delete()
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      dev = to_string(model%orient%deviation_detected, "(es9.2)")
      call check(error, model%orient%group_used == "Cs" .and. index(model%orient%reason, &
         & "detected C2v, using Cs (C2v matches only to "//dev//" bohr, above the tolerance 1.00E-08 bohr)") > 0, &
         & "distorted water: reason '"//model%orient%reason//"'")
      if (allocated(error)) return
      call check(error, nout >= 1 .and. any(index(out(1:nout), "matches C2v only to "//dev//" bohr") > 0 .and. &
         & index(out(1:nout), "so Cs is used") > 0), "distorted water: no warning")
   end subroutine test_tolerance_fallback

   !* ================================================================================= *!
   !*                                      Outputs                                      *!
   !* ================================================================================= *!

   !> The `symmetry_*` fields: numeric only, from the update on
   !>
   !> @param[out] error Error handle
   subroutine test_fields(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      type(field_query_type) :: query
      integer :: kind, pos, iname
      character(len=*), parameter :: names(13) = [character(len=24) :: "symmetry_orientation", &
         & "symmetry_order", "symmetry_deviation", "symmetry_frame_level", "symmetry_frame_gap", &
         & "symmetry_nrot", "symmetry_near_gap", "symmetry_near_linear", "symmetry_near_harmonic", &
         & "symmetry_tied_frames", "symmetry_rot_index", "symmetry_orbit", "symmetry_rotation"]
      character(len=*), parameter :: removed(4) = [character(len=24) :: "symmetry_group_requested", &
         & "symmetry_group_detected", "symmetry_group_used", "symmetry_reason"]

      call new_context(ctx, nthreads=0, verbosity=0)
      call structure_of("water", mol)
      do kind = 1, size(kind_name)
         settings%orientation = "molecular"
         settings%group = "auto"
         call updated_cavity(kind, mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, trim(kind_name(kind))//": "//err%message)
            return
         end if
         call query%enumerate()
         call cav%list_fields(query)
         do iname = 1, size(names)
            call check(error, query%index_of(trim(names(iname))) > 0, &
               & trim(kind_name(kind))//": "//trim(names(iname))//" not listed")
            if (allocated(error)) return
         end do
         do iname = 1, size(removed)
            call check(error, query%index_of(trim(removed(iname))) == 0, &
               & trim(kind_name(kind))//": text field "//trim(removed(iname))//" still listed")
            if (allocated(error)) return
         end do
         pos = query%index_of("symmetry_rotation")
         call check(error, query%info(pos)%rank == 3 .and. all(query%info(pos)%dims == [3, 3, mol%nat]), &
            & trim(kind_name(kind))//": rotation shape")
         if (allocated(error)) return

         call query%fetch("symmetry_order")
         call cav%list_fields(query)
         call check(error, query%found .and. all(query%ivals == [4]), trim(kind_name(kind))//": order field")
         if (allocated(error)) return
         call query%fetch("symmetry_orbit")
         call cav%list_fields(query)
         call check(error, query%found .and. all(query%ivals == [0, 1, 1]), trim(kind_name(kind))//": orbit field")
         if (allocated(error)) return
         call query%fetch("symmetry_orientation")
         call cav%list_fields(query)
         call check(error, query%found .and. all(query%ivals == [1]), trim(kind_name(kind))//": orientation field")
         if (allocated(error)) return
         call query%fetch("symmetry_deviation")
         call cav%list_fields(query)
         call check(error, query%found .and. query%rvals(1) <= TOL, trim(kind_name(kind))//": deviation field")
         if (allocated(error)) return
         call query%fetch("symmetry_rotation")
         call cav%list_fields(query)
         associate (gm_ => symmetry_of(cav))
            call check(error, query%found .and. &
               & maxval(abs(query%rvals(1:9) - reshape(model_rotation(gm_, 1), [9]))) < 1.0e-15_wp, &
               & trim(kind_name(kind))//": rotation field")
         end associate
         if (allocated(error)) return

         ! Off: lab frame, order 1
         settings%orientation = "lab"
         settings%group = "off"
         call updated_cavity(kind, mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, trim(kind_name(kind))//": "//err%message)
            return
         end if
         call query%fetch("symmetry_orientation")
         call cav%list_fields(query)
         call check(error, query%found .and. all(query%ivals == [0]), trim(kind_name(kind))//": lab orientation field")
         if (allocated(error)) return
         call query%fetch("symmetry_order")
         call cav%list_fields(query)
         call check(error, query%found .and. all(query%ivals == [1]), trim(kind_name(kind))//": off order field")
         if (allocated(error)) return
      end do
   end subroutine test_fields

   !> The symmetry report for off, auto and explicit groups, and its place
   !>
   !> - the model prints it; a cavity prints it right after its radii, never
   !>   in its results
   !>
   !> @param[out] error Error handle
   subroutine test_report(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      character(len=line_len) :: lines(64)
      integer :: nline, iradii, isym

      call new_context(ctx, nthreads=0, verbosity=0)

      ! Off: orientation and group only
      call structure_of("water", mol)
      call updated_cavity(kind_iswig, mol, settings, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call capture(cav, 3, lines, nline)
      call check(error, printed_entry(lines, "Orientation", "lab") .and. printed_entry(lines, "Group", "off"), &
         & "off report lacks orientation or group")
      if (allocated(error)) return
      call check(error,.not. printed_entry(lines, "Group used", "") .and. &
         & .not. any(index(lines, "Symmetry-equivalent atoms") > 0), "off report prints a group or table")
      if (allocated(error)) return

      ! Auto: groups, frame and the orbit table from level 2
      settings%orientation = "molecular"
      settings%group = "auto"
      call updated_cavity(kind_iswig, mol, settings, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call capture(cav, 2, lines, nline)
      call check(error, printed_entry(lines, "Orientation", "molecular") .and. &
         & printed_entry(lines, "Group requested", "auto") .and. &
         & printed_entry(lines, "Group detected", "C2v (deviation ") .and. &
         & printed_entry(lines, "Group used", "C2v (order 4)") .and. printed_entry(lines, "Orbits", "2") .and. &
         & printed_entry(lines, "Molecular frame", "eigen") .and. &
         & printed_entry(lines, "Max. deviation", "bohr (tolerance 1.00E-08 bohr)"), &
         & "auto report lacks a key line")
      if (allocated(error)) return
      call check(error,.not. printed_entry(lines, "Reason", ""), "auto water reports a reason")
      if (allocated(error)) return
      call check(error, any(index(lines, "Symmetry-equivalent atoms") > 0) .and. &
         & count(index(lines, " O ") > 0) == 1 .and. count(index(lines, "2,3") > 0) == 1, &
         & "auto report lacks the orbit table")
      if (allocated(error)) return
      call capture(cav, 1, lines, nline)
      call check(error, printed_entry(lines, "Group used", "C2v") .and. &
         & .not. any(index(lines, "Symmetry-equivalent atoms") > 0), "level 1 prints the orbit table")
      if (allocated(error)) return

      ! The radii printout carries the report right after the radii table
      open (newunit=isym, status="scratch", action="readwrite")
      call cav%print_radii(isym, 2)
      call read_printout(isym, lines, nline)
      close (isym)
      iradii = findloc(index(lines(1:nline), "Used radii") > 0, .true., dim=1)
      isym = findloc(index(lines(1:nline), "Symmetry:") > 0, .true., dim=1)
      call check(error, iradii > 0 .and. isym > iradii .and. printed_entry(lines, "Group used", "C2v"), &
         & "radii printout lacks the report after the radii")
      if (allocated(error)) return

      ! Explicit D3h on BF3: used as C3v, with the reason
      call structure_of("bf3", mol)
      settings%group = "D3h"
      call updated_cavity(kind_drop, mol, settings, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call capture(cav, 3, lines, nline)
      call check(error, printed_entry(lines, "Group requested", "D3h") .and. &
         & printed_entry(lines, "Group detected", "D3h (deviation") .and. &
         & printed_entry(lines, "Group used", "C3v") .and. printed_entry(lines, "Reason", "requested D3h"), &
         & "explicit BF3 report lacks a key line")
      if (allocated(error)) return
      call check(error, any(index(lines, "Symmetry-equivalent atoms") > 0) .and. &
         & count(index(lines, "2-4") > 0) == 1, "explicit BF3 report lacks the fluorine orbit")
      if (allocated(error)) return

      ! Explicit C1: no detection, no table
      settings%group = "C1"
      call updated_cavity(kind_numsa, mol, settings, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call capture(cav, 3, lines, nline)
      call check(error, printed_entry(lines, "Group used", "C1") .and. printed_entry(lines, "Group detected", "-") &
         & .and. printed_entry(lines, "Tolerance", "1.00E-08 bohr") .and. &
         & .not. any(index(lines, "Symmetry-equivalent atoms") > 0), "explicit C1 report")
      if (allocated(error)) return

      ! Nothing at level 0; the results printout carries no report anymore
      call capture(cav, 0, lines, nline)
      call check(error, nline, 0, "level 0 prints the report")
      if (allocated(error)) return
      call capture_info(cav, lines, nline)
      call check(error, any(index(lines, "Results:") > 0) .and. .not. any(index(lines, "Symmetry:") > 0), &
         & "cavity results still carry the symmetry report")

   contains

      !> Print the report at a level into a scratch unit
      !>
      !> @param[in]  cav   Updated cavity
      !> @param[in]  level Verbosity level
      !> @param[out] lines Printed lines
      !> @param[out] nline Number of lines
      subroutine capture(cav, level, lines, nline)
         !> Updated cavity
         class(cavity_type), intent(in) :: cav
         !> Verbosity level
         integer, intent(in) :: level
         !> Printed lines
         character(len=line_len), intent(out) :: lines(:)
         !> Number of lines
         integer, intent(out) :: nline

         integer :: unit

         open (newunit=unit, status="scratch", action="readwrite")
         call cav%print_grid_model(unit, level)
         call read_printout(unit, lines, nline)
         close (unit)
      end subroutine capture

      !> Print the cavity results into a scratch unit
      !>
      !> @param[in]  cav   Updated cavity
      !> @param[out] lines Printed lines
      !> @param[out] nline Number of lines
      subroutine capture_info(cav, lines, nline)
         !> Updated cavity
         class(cavity_type), intent(in) :: cav
         !> Printed lines
         character(len=line_len), intent(out) :: lines(:)
         !> Number of lines
         integer, intent(out) :: nline

         integer :: unit

         open (newunit=unit, status="scratch", action="readwrite")
         call cav%print(unit)
         call read_printout(unit, lines, nline)
         close (unit)
      end subroutine capture_info

   end subroutine test_report

   !* ================================================================================= *!
   !*                                 iSwiG orientation                                 *!
   !* ================================================================================= *!

   !> Lab orientation: points and normals bitwise those of c + r u
   !>
   !> The construction before sphere orientations, component by component
   !> from the cached nodes, for an explicit and the default lab orientation;
   !> no rotated copy of the nodes is made
   !>
   !> @param[out] error Error handle
   subroutine test_iswig_lab_bitwise(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      real(wp) :: p(3, 1), n(3, 1)
      integer :: ipt, iat, ileb, imode

      call new_context(ctx, nthreads=0, verbosity=0)
      call get_structure(mol, "MB16-43", "01")
      do imode = 1, 2
         settings%orientation = merge("lab ", "auto", imode == 1)
         call new_iswig(mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         associate (gm_ => symmetry_of(cav))
            call check(error, gm_%is_lab() .and. .not. (.not. cav%grid%plain), &
               & "lab orientation built rotated nodes")
         end associate
         if (allocated(error)) return
         do ipt = 1, cav%ngrid
            iat = cav%owner(ipt)
            ileb = cav%numbering(ipt) - (iat - 1)*cav%param%num_leb
            p(1, 1) = cav%sphxyz(1, iat) + cav%radii(iat)*cav%ang_grid(1, ileb)
            p(2, 1) = cav%sphxyz(2, iat) + cav%radii(iat)*cav%ang_grid(2, ileb)
            p(3, 1) = cav%sphxyz(3, iat) + cav%radii(iat)*cav%ang_grid(3, ileb)
            n(:, 1) = (p(:, 1) - cav%sphxyz(:, iat))/cav%radii(iat)
            call check(error, same_bits(cav%xyz(:, ipt:ipt), p) .and. same_bits(cav%normal0(:, ipt:ipt), n), &
               & "lab point or normal is not the construction from c + r u")
            if (allocated(error)) return
            call check(error, cav%wleb(ipt) == cav%ang_weight(ileb), "lab weight changed")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_iswig_lab_bitwise

   !> Base frame helpers against the dense contraction of the rotation Jacobians
   !>
   !> Benzene with its group (several distinct rotations): every add procedure
   !> of [[grid_frame_type]] matches sum_k jac(k, i, j, rot(a)) t(k), and the
   !> base node cache matches the model's rotation indices
   !>
   !> @param[out] error Error handle
   subroutine test_frame_helpers(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: tol = 1.0e-13_wp
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      type(grid_frame_type) :: frame
      real(wp), allocatable :: jac(:, :, :, :), t(:, :), tv(:, :, :), s_sph(:, :)
      real(wp), allocatable :: q(:, :), q_ref(:, :), qv(:, :, :), qv_ref(:, :, :)
      real(wp), allocatable :: qs(:, :, :), qs_ref(:, :, :)
      integer :: nat, a, i, j, k, c, irot

      call new_context(ctx, nthreads=0, verbosity=0)
      call structure_of("benzene", mol)
      settings%group = "auto"
      call new_iswig(mol, settings, ctx, cav, err, num_leb=26)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      nat = mol%nat
      associate (gm_ => symmetry_of(cav))
         call check(error, model_nrot(gm_) > 1, "benzene should need several rotations")
      end associate
      if (allocated(error)) return
      call check(error, allocated(cav%grid%comp_rot) .and. (.not. cav%grid%plain), "node cache missing")
      if (allocated(error)) return
      associate (gm_ => symmetry_of(cav))
         call check(error, all([(cav%grid%comp_rot(a) == model_rot_index(gm_, a), a=1, nat)]), &
            & "node cache rotation indices")
      end associate
      if (allocated(error)) return

      call cav%get_frame(frame)
      associate (gm_ => symmetry_of(cav))
         call model_frame_jacobian(gm_, jac)
      end associate
      allocate (t(3, nat), tv(3, 3, nat), s_sph(3, nat))
      do a = 1, nat
         t(:, a) = sin(0.7_wp*real(a, wp) + [0.1_wp, 1.3_wp, 2.6_wp])
         s_sph(:, a) = cos(0.4_wp*real(a, wp) + [0.5_wp, 1.1_wp, 2.2_wp])
         do c = 1, 3
            tv(:, c, a) = sin(0.3_wp*real(a + 7*c, wp) + [0.2_wp, 0.9_wp, 1.7_wp])
         end do
      end do

      ! One point per owner sphere, scalar and vector
      allocate (q(3, nat), q_ref(3, nat), qv(3, 3, nat), qv_ref(3, 3, nat))
      do a = 1, nat
         associate (gm_ => symmetry_of(cav))
            irot = model_rot_index(gm_, a)
         end associate
         q = 1.0_wp
         qv = 1.0_wp
         q_ref = 1.0_wp
         qv_ref = 1.0_wp
         do j = 1, nat
            do i = 1, 3
               do k = 1, 3
                  q_ref(i, j) = q_ref(i, j) + jac(k, i, j, irot)*t(k, a)
                  do c = 1, 3
                     qv_ref(c, i, j) = qv_ref(c, i, j) + jac(k, i, j, irot)*tv(k, c, a)
                  end do
               end do
            end do
         end do
         call frame%add_scalar(a, t(:, a), q)
         call frame%add_vector(a, tv(:, :, a), qv)
         call check(error, maxval(abs(q - q_ref)) < tol, "add_scalar")
         if (allocated(error)) return
         call check(error, maxval(abs(qv - qv_ref)) < tol, "add_vector")
         if (allocated(error)) return
      end do

      ! Per-sphere sums and the total
      allocate (qs(3, nat, nat), qs_ref(3, nat, nat))
      qs = 0.0_wp
      qs_ref = 0.0_wp
      q = 0.0_wp
      q_ref = 0.0_wp
      do a = 1, nat
         associate (gm_ => symmetry_of(cav))
            irot = model_rot_index(gm_, a)
         end associate
         do j = 1, nat
            do i = 1, 3
               do k = 1, 3
                  qs_ref(i, a, j) = qs_ref(i, a, j) + jac(k, i, j, irot)*s_sph(k, a)
                  q_ref(i, j) = q_ref(i, j) + jac(k, i, j, irot)*s_sph(k, a)
               end do
            end do
         end do
      end do
      call frame%add_spheres(s_sph, qs)
      call frame%add_total(s_sph, q)
      call check(error, maxval(abs(qs - qs_ref)) < tol, "add_spheres")
      if (allocated(error)) return
      call check(error, maxval(abs(q - q_ref)) < tol, "add_total")
   end subroutine test_frame_helpers

   !> Molecular orientation: sphere a carries the nodes R_a u, weights unchanged
   !>
   !> One rotated copy of the nodes per distinct rotation; frame only (C1) and
   !> with a group (ammonia, C3v)
   !>
   !> @param[out] error Error handle
   subroutine test_iswig_molecular_points(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      character(len=*), parameter :: names(2) = [character(len=7) :: "c1", "ammonia"]
      character(len=*), parameter :: groups(2) = [character(len=4) :: "off", "auto"]
      real(wp) :: rot(3, 3), u(3), zeta
      integer :: im, ipt, iat, ileb

      call new_context(ctx, nthreads=0, verbosity=0)
      do im = 1, size(names)
         call structure_of(trim(names(im)), mol)
         settings%orientation = "molecular"
         settings%group = groups(im)
         call new_iswig(mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         associate (gm_ => symmetry_of(cav))
            call check(error,.not. gm_%is_lab() .and. (.not. cav%grid%plain), &
               & trim(names(im))//": no rotated nodes")
         end associate
         if (allocated(error)) return
         associate (gm_ => symmetry_of(cav))
            call check(error, all(shape(cav%grid%node) == [3, cav%param%num_leb*model_nrot(gm_)]), &
               & trim(names(im))//": one rotated copy per distinct rotation")
         end associate
         if (allocated(error)) return
         zeta = cav%cached_swx
         do ipt = 1, cav%ngrid
            iat = cav%owner(ipt)
            ileb = cav%numbering(ipt) - (iat - 1)*cav%param%num_leb
            associate (gm_ => symmetry_of(cav))
               rot = model_rotation(gm_, iat)
            end associate
            u = matmul(rot, cav%ang_grid(:, ileb))
            call check(error, maxval(abs(cav%xyz(:, ipt) - cav%sphxyz(:, iat) - cav%radii(iat)*u)) < 1.0e-13_wp &
               & .and. maxval(abs(cav%normal0(:, ipt) - u)) < 1.0e-13_wp, &
               & trim(names(im))//": point is not c + r R u")
            if (allocated(error)) return
            call check(error, cav%wleb(ipt) == cav%ang_weight(ileb) .and. &
               & cav%xi0(ipt) == zeta/(cav%radii(iat)*sqrt(cav%ang_weight(ileb))), &
               & trim(names(im))//": weight or width changed by the rotation")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_iswig_molecular_points

   !> Molecular orientation: area, volume and a CPCM energy follow the molecule
   !>
   !> - rotated, translated and listed in another order, to 1e-12 relative
   !> - frame only on asymmetric and exactly symmetric molecules: a tie of the
   !>   frame's discrete choices gives the image of the cavity under a
   !>   symmetry of the atoms, so the host charges depend on the element only
   !> - with a group, including sites whose admissible orientations tie in
   !>   their closeness to the frame (the C3v sites of ammonia, CH3F and
   !>   ethane): the neighbour form breaks the tie
   !> - the lab orientation does change, so the check has teeth
   !>
   !> @param[out] error Error handle
   subroutine test_iswig_invariance(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(14) = [character(len=7) :: "c1", "mb16", "water", "ammonia", &
         & "bf3", "benzene", "water", "ammonia", "methane", "benzene", "sf6", "hocl", "ch3f", "ethane"]
      character(len=*), parameter :: groups(14) = [character(len=4) :: "off", "off", "off", "off", &
         & "off", "off", "auto", "auto", "auto", "auto", "auto", "auto", "auto", "auto"]
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, moved
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      real(wp) :: ref(3), val(3), q(3, 3)
      integer :: im

      call new_context(ctx, nthreads=0, verbosity=0)
      q = test_rotation(1)
      do im = 1, size(names)
         call structure_of(trim(names(im)), mol)
         call moved_structure(mol, q, [1.7_wp, -0.6_wp, 2.9_wp], moved)
         settings%orientation = "molecular"
         settings%group = groups(im)
         call invariance_values(mol, settings, ctx, ref, error)
         if (allocated(error)) return
         call invariance_values(moved, settings, ctx, val, error)
         if (allocated(error)) return
         call check(error, all(abs(val - ref) <= 1.0e-12_wp*abs(ref)), trim(names(im))//" ("// &
            & trim(groups(im))//"): area, volume or CPCM energy not invariant")
         if (allocated(error)) return
      end do

      ! The lab orientation is not invariant
      call structure_of("c1", mol)
      call moved_structure(mol, q, [1.7_wp, -0.6_wp, 2.9_wp], moved)
      settings%orientation = "lab"
      settings%group = "off"
      call invariance_values(mol, settings, ctx, ref, error)
      if (allocated(error)) return
      call invariance_values(moved, settings, ctx, val, error)
      if (allocated(error)) return
      call check(error, abs(val(1) - ref(1)) > 1.0e-9_wp*abs(ref(1)), "lab orientation looks invariant")
   end subroutine test_iswig_invariance

   !> With a group, the point set maps onto itself under every operation used
   !>
   !> The image of every point lands on a point of the image sphere with the
   !> same weight, area, switching and width and the image normal
   !>
   !> @param[out] error Error handle
   subroutine test_iswig_exact_symmetry(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(7) = [character(len=7) :: "water", "ammonia", "methane", &
         & "bf3", "benzene", "sf6", "hocl"]
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      real(wp) :: g(3, 3), image(3), dist, best
      integer :: im, ig, ip, jp, jbest, ib
      integer, allocatable :: count_a(:)

      call new_context(ctx, nthreads=0, verbosity=0)
      settings%group = "auto"
      do im = 1, size(names)
         call structure_of(trim(names(im)), mol)
         call new_iswig(mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, trim(names(im))//": "//err%message)
            return
         end if
         associate (sym_ => symmetry_of(cav))
            call check(error, size(sym_%orient%ops, 3) > 1, trim(names(im))//": no symmetry used")
            if (allocated(error)) return
            count_a = [(count(cav%owner == ib), ib=1, mol%nat)]
            do ig = 1, size(sym_%orient%ops, 3)
               g = sym_%orient%ops(:, :, ig)
               do ip = 1, cav%ngrid
                  ib = sym_%orient%perm(cav%owner(ip), ig)
                  call check(error, count_a(ib) == count_a(cav%owner(ip)), &
                     & trim(names(im))//": equivalent spheres keep different numbers of points")
                  if (allocated(error)) return
                  image = sym_%orient%center + matmul(g, cav%xyz(:, ip) - sym_%orient%center)
                  best = huge(1.0_wp)
                  jbest = 0
                  do jp = 1, cav%ngrid
                     if (cav%owner(jp) /= ib) cycle
                     dist = norm2(cav%xyz(:, jp) - image)
                     if (dist < best) then
                        best = dist
                        jbest = jp
                     end if
                  end do
                  call check(error, best < 1.0e-9_wp, trim(names(im))//": image of a point is no point")
                  if (allocated(error)) return
                  call check(error, maxval(abs(cav%normal0(:, jbest) - matmul(g, cav%normal0(:, ip)))) &
                     & < 1.0e-9_wp, trim(names(im))//": image normal differs")
                  if (allocated(error)) return
                  call check(error, cav%wleb(jbest) == cav%wleb(ip) .and. &
                     & abs(cav%a(jbest) - cav%a(ip)) <= 1.0e-10_wp*cav%a(ip) .and. &
                     & abs(cav%f(jbest) - cav%f(ip)) <= 1.0e-10_wp*cav%f(ip) .and. &
                     & cav%xi0(jbest) == cav%xi0(ip), &
                     & trim(names(im))//": image point has another weight, area, switching or width")
                  if (allocated(error)) return
               end do
            end do
         end associate
      end do
   end subroutine test_iswig_exact_symmetry

   !> Frame gradient at the eigen level: an asymmetric molecule, frame only
   !>
   !> Every forward array and the reverse contraction against central
   !> differences along each Cartesian coordinate; area and volume carry no
   !> net force or torque, which without the frame term they would
   !>
   !> @param[out] error Error handle
   subroutine test_iswig_gradient_eigen(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: grads(:, :, :)
      real(wp) :: torque(3)
      integer :: iat

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_frame_fixture("c1", "off", fix, ctx, error)
      if (allocated(error)) return
      call check(error, fix%level, frame_level_eigen, "c1: frame not at the eigen level")
      if (allocated(error)) return
      call check_frame_gradient(fix, ctx, cartesian_directions(fix%mol%nat), FD_STEP, FD_ABS, FD_REL, &
         & error, grads)
      if (allocated(error)) return
      call check_rigid(fix, grads, error)
      if (allocated(error)) return

      ! Lab orientation: the area does feel a torque
      fix%settings%orientation = "lab"
      call new_iswig(fix%mol, fix%settings, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call cav%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      torque = 0.0_wp
      do iat = 1, fix%mol%nat
         torque = torque + cross_product(fix%mol%xyz(:, iat), cav%area_grad(:, iat))
      end do
      call check(error, norm2(torque) > 1.0e-6_wp, "lab area has no torque, the rigid check is vacuous")
   end subroutine test_iswig_gradient_eigen

   !> Frame gradient at the plane level
   !>
   !> An asymmetric molecule whose second moment has two equal eigenvalues:
   !> axis from the separated eigenvalue, in-plane angle from a harmonic;
   !> every stencil point stays at that level and branch (no symmetry, so no
   !> tie between the in-plane candidates)
   !>
   !> @param[out] error Error handle
   subroutine test_iswig_gradient_plane(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      real(wp), allocatable :: grads(:, :, :)

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_frame_fixture("c1_plane", "off", fix, ctx, error)
      if (allocated(error)) return
      call check(error, fix%level, frame_level_plane, "c1_plane: frame not at the plane level")
      if (allocated(error)) return
      call check_frame_gradient(fix, ctx, cartesian_directions(fix%mol%nat), 2.0_wp*FD_STEP, FD_ABS, FD_REL, &
         & error, grads)
      if (allocated(error)) return
      call check_rigid(fix, grads, error)
   end subroutine test_iswig_gradient_plane

   !> Frame gradient with a group: directions that keep it
   !>
   !> - ammonia (C3v): fluorine-free family sites, the frame does not turn
   !> - HOCl (Cs): every sphere turns about the plane normal with the frame
   !> - benzene (D6h): entries shared by spheres modulo Oh
   !>
   !> @param[out] error Error handle
   subroutine test_iswig_gradient_symmetric(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(3) = [character(len=7) :: "ammonia", "hocl", "benzene"]
      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: grads(:, :, :), dirs(:, :, :)
      integer :: im

      call new_context(ctx, nthreads=0, verbosity=0)
      do im = 1, size(names)
         call new_frame_fixture(trim(names(im)), "auto", fix, ctx, error)
         if (allocated(error)) return
         call new_iswig(fix%mol, fix%settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         associate (gm_ => symmetry_of(cav))
            dirs = group_directions(gm_%orient%ops, gm_%orient%perm, 3)
         end associate
         call check_frame_gradient(fix, ctx, dirs, FD_STEP, FD_ABS, FD_REL, error, grads)
         if (allocated(error)) return
      end do
   end subroutine test_iswig_gradient_symmetric

   !> Pseudo-random displacement directions averaged over a group
   !>
   !> Atom a moves along the average of g^T v(perm(a, g)) over the operations g,
   !> so every direction keeps the group; scaled to a largest entry of one
   !>
   !> @param[in] ops  Operations of the group (3, 3, order)
   !> @param[in] perm Image atom of every atom under every operation (nat, order)
   !> @param[in] ndir Number of directions
   pure function group_directions(ops, perm, ndir) result(dirs)
      !> Operations of the group
      real(wp), intent(in) :: ops(:, :, :)
      !> Image atoms
      integer, intent(in) :: perm(:, :)
      !> Number of directions
      integer, intent(in) :: ndir
      !> Directions (3, nat, ndir)
      real(wp), allocatable :: dirs(:, :, :)

      integer :: idir, iat, ig

      allocate (dirs(3, size(perm, 1), ndir), source=0.0_wp)
      do idir = 1, ndir
         do iat = 1, size(perm, 1)
            do ig = 1, size(ops, 3)
               dirs(:, iat, idir) = dirs(:, iat, idir) + matmul(transpose(ops(:, :, ig)), &
                  & sin(0.9_wp*real(3*perm(iat, ig) + idir, wp) + [0.3_wp, 1.4_wp, 2.6_wp]))
            end do
         end do
         dirs(:, :, idir) = dirs(:, :, idir)/maxval(abs(dirs(:, :, idir)))
      end do
   end function group_directions

   !* ================================================================================= *!
   !*                                 DROP orientation                                  *!
   !* ================================================================================= *!

   !> DROP anchors: sphere a carries the Lebedev nodes R_a u
   !>
   !> - molecular orientation (frame only, and C3v ammonia): one rotated copy
   !>   of the nodes per distinct rotation, R_k u to round-off; every anchor is
   !>   exactly c_a + r_a times its sphere's copy, weight and width unchanged
   !> - lab orientation: no copy, every anchor bitwise c + r u
   !>
   !> @param[out] error Error handle
   subroutine test_drop_anchors(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      character(len=*), parameter :: names(3) = [character(len=7) :: "c1", "ammonia", "c1"]
      character(len=*), parameter :: groups(3) = [character(len=4) :: "off", "auto", "off"]
      character(len=*), parameter :: orientations(3) = [character(len=9) :: "molecular", "molecular", "lab"]
      real(wp) :: rot(3, 3), p(3, 1), u(3)
      integer :: im, ia, iat, ileb, irot, k

      call new_context(ctx, nthreads=0, verbosity=0)
      do im = 1, size(names)
         call structure_of(trim(names(im)), mol)
         settings%orientation = trim(orientations(im))
         settings%group = trim(groups(im))
         call updated_cavity(kind_drop, mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         select type (cav)
         type is (cavity_type_drop)
            if (settings%orientation == "lab") then
               associate (gm_ => symmetry_of(cav))
                  call check(error, gm_%is_lab() .and. .not. (.not. cav%grid%plain), &
                     & "DROP lab orientation built rotated nodes")
               end associate
               if (allocated(error)) return
            else
               associate (gm_ => symmetry_of(cav))
                  call check(error,.not. gm_%is_lab() .and. (.not. cav%grid%plain), &
                     & trim(names(im))//": no rotated DROP nodes")
               end associate
               if (allocated(error)) return
               associate (gm_ => symmetry_of(cav))
                  call check(error, all(shape(cav%grid%node) == [3, cav%param%num_leb*model_nrot(gm_)]), &
                     & trim(names(im))//": one rotated copy of the nodes per distinct rotation")
               end associate
               if (allocated(error)) return
               do irot = 1, cav%grid%nset
                  associate (gm_ => symmetry_of(cav))
                     rot = gm_%orient%rotation(:, :, irot)
                  end associate
                  do ileb = 1, cav%param%num_leb
                     u = matmul(rot, cav%ang_grid(:, ileb))
                     call check(error, maxval(abs(cav%grid%node(:, (irot - 1)*size(cav%ang_grid, 2) + ileb) - u)) &
                        & <= 4.0_wp*epsilon(1.0_wp), trim(names(im))//": rotated node is not R u")
                     if (allocated(error)) return
                  end do
               end do
            end if
            call check(error, cav%anchor%ngrid > 0, trim(names(im))//": no DROP anchors")
            if (allocated(error)) return
            do ia = 1, cav%anchor%ngrid
               iat = cav%anchor%owner(ia)
               ileb = cav%anchor%id(ia) - (iat - 1)*cav%param%num_leb
               if (settings%orientation == "lab") then
                  do k = 1, 3
                     p(k, 1) = cav%sphxyz(k, iat) + cav%radii(iat)*cav%ang_grid(k, ileb)
                  end do
               else
                  associate (gm_ => symmetry_of(cav))
                     irot = model_rot_index(gm_, iat)
                  end associate
                  p(:, 1) = cav%sphxyz(:, iat) + cav%radii(iat)*cav%grid%node(:, (irot - 1)*size(cav%ang_grid, 2) + ileb)
               end if
               call check(error, same_bits(cav%anchor%xyz(:, ia:ia), p), &
                  & trim(names(im))//" ("//trim(orientations(im))//"): anchor is not c + r R u")
               if (allocated(error)) return
               call check(error, cav%anchor%wleb0(ia) == cav%ang_weight(ileb), &
                  & trim(names(im))//": anchor weight changed by the rotation")
               if (allocated(error)) return
            end do
         class default
            call test_failed(error, "not a DROP cavity")
            return
         end select
      end do
   end subroutine test_drop_anchors

   !> DROP frame gradient at the eigen level: an asymmetric molecule, frame only
   !>
   !> - area, volume, a generic surface adjoint energy and a CPCM energy, every
   !>   forward total and reverse contraction, against central differences
   !>   along every Cartesian coordinate; area and volume carry no net force
   !>   or torque
   !> - every per-point forward array against the 6-point differences of the
   !>   DROP per-point check, with its settings and tolerances
   !> - without the frame term every analytic gradient misses the differences
   !>
   !> @param[out] error Error handle
   subroutine test_drop_gradient_eigen(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_drop_frame_gradient("c1", frame_level_eigen, error)
   end subroutine test_drop_gradient_eigen

   !> DROP frame gradient at the plane level
   !>
   !> As the eigen level, for the near-symmetric propeller; the per-point
   !> stencil (6 points, step 2.5e-5 bohr) stays at the plane level
   !>
   !> - the surface adjoint energy weights both principal curvatures of a
   !>   point alike: the axis of the propeller passes through a node (a lab
   !>   axis or the body diagonal of the grid), an umbilic point of the
   !>   surface up to the distortion, where k1 and k2 cross within the
   !>   stencil; in the lab orientation they miss the differences as well
   !>
   !> @param[out] error Error handle
   subroutine test_drop_gradient_plane(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_drop_frame_gradient("c1_plane", frame_level_plane, error, row_step=2.5e-5_wp, &
         & shared_curvature=.true.)
   end subroutine test_drop_gradient_plane

   !> DROP frame gradient of a molecule without a group, totals and per point
   !>
   !> @param[in]  name     Molecule
   !> @param[in]  level    Frame level the reference must have
   !> @param[out] error    Test failure
   !> @param[in]  row_step Step of the per-point check, bohr (default its own)
   !> @param[in]  shared_curvature Whether the surface adjoint energy weights
   !>                      both principal curvatures of a point alike (default no)
   subroutine check_drop_frame_gradient(name, level, error, row_step, shared_curvature)
      !> Molecule
      character(len=*), intent(in) :: name
      !> Frame level
      integer, intent(in) :: level
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Step of the per-point check
      real(wp), intent(in), optional :: row_step
      !> Whether both principal curvatures share one weight
      logical, intent(in), optional :: shared_curvature

      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      type(grid_model_symmetry_type) :: model
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: grads(:, :, :)

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_frame_fixture(name, "off", fix, ctx, error, kind=kind_drop)
      if (allocated(error)) return
      if (present(shared_curvature)) fix%shared_curvature = shared_curvature
      call check(error, fix%level, level, fix%name//": frame not at the expected level")
      if (allocated(error)) return
      call check_frame_gradient(fix, ctx, cartesian_directions(fix%mol%nat), fix%step, FD_ABS, FD_REL, &
         & error, grads)
      if (allocated(error)) return
      call check_rigid(fix, grads, error)
      if (allocated(error)) return

      ! The frame term is needed
      call check_frame_gradient(fix, ctx, cartesian_directions(fix%mol%nat), fix%step, FD_ABS, FD_REL, &
         & error, grads, frame_off=.true.)
      if (allocated(error)) return

      ! Every per-point forward array
      call new_grid_model_symmetry(model, err, param=fix%settings)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check_drop_row_gradient(error, fix%mol, step=row_step, symmetry=model)
   end subroutine check_drop_frame_gradient

   !> DROP frame gradient with a group: directions that keep it
   !>
   !> - ammonia (C3v): the frame does not turn along such directions
   !> - HOCl (Cs): every sphere turns about the plane normal with the frame,
   !>   so without the frame term every analytic gradient misses
   !>
   !> @param[out] error Error handle
   subroutine test_drop_gradient_symmetric(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(2) = [character(len=7) :: "ammonia", "hocl"]
      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: grads(:, :, :), dirs(:, :, :)
      integer :: im

      call new_context(ctx, nthreads=0, verbosity=0)
      do im = 1, size(names)
         call new_frame_fixture(trim(names(im)), "auto", fix, ctx, error, kind=kind_drop)
         if (allocated(error)) return
         call new_frame_cavity(fix, fix%mol, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         associate (gm_ => symmetry_of(cav))
            dirs = group_directions(gm_%orient%ops, gm_%orient%perm, 3)
         end associate
         call check_frame_gradient(fix, ctx, dirs, fix%step, FD_ABS, FD_REL, error, grads)
         if (allocated(error)) return
         if (trim(names(im)) == "hocl") then
            call check_frame_gradient(fix, ctx, dirs, fix%step, FD_ABS, FD_REL, error, grads, frame_off=.true.)
            if (allocated(error)) return
         end if
      end do
   end subroutine test_drop_gradient_symmetric

   !> DROP computing its own symmetry model times it once per update
   !>
   !> - "Setup" counts one call per update and holds "Symmetry"; a model the
   !>   host already updated opens no "Symmetry" under the cavity
   !>
   !> @param[out] error Error handle
   subroutine test_drop_setup_timer(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, moved
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      character(len=80) :: label
      integer :: istep

      call structure_of("ammonia", mol)
      call rigid_motion(mol, test_rotation(1), [0.4_wp, -1.1_wp, 0.7_wp], test_permutation(mol%nat, 1), moved)
      settings%group = "auto"
      call new_context(ctx, nthreads=0, verbosity=0)
      call updated_cavity(kind_drop, mol, settings, ctx, cav, err)
      do istep = 1, 2
         if (allocated(err)) exit
         call cav%update(moved, err)
         if (allocated(err)) exit
         call cav%update(mol, err)
      end do
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, timer_path(ctx, "Setup/Symmetry/Detection"), "DROP: symmetry not under Setup")
      if (allocated(error)) return
      write (label, '("DROP: Setup ",i0,"x, Symmetry ",i0,"x in 5 updates")') timer_calls(ctx, "Setup"), &
         & timer_calls(ctx, "Symmetry")
      call check(error, timer_calls(ctx, "Setup") == 5 .and. timer_calls(ctx, "Symmetry") == 5, &
         & trim(label)//", not one call per update")
   end subroutine test_drop_setup_timer

   !* ================================================================================= *!
   !*                                 NumSA orientation                                 *!
   !* ================================================================================= *!

   !> NumSA nodes: sphere a integrates over R_a u
   !>
   !> - lab orientation (explicit and default): no rotated copy, no torque
   !> - molecular orientation (frame only, and C3v ammonia): one rotated copy
   !>   of the nodes per distinct rotation, R_k u to round-off, and a torque
   !>   per sphere; the sphere areas add up to the total
   !>
   !> @param[out] error Error handle
   subroutine test_numsa_nodes(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(mctc_error), allocatable :: err
      character(len=*), parameter :: names(4) = [character(len=7) :: "c1", "c1", "c1", "ammonia"]
      character(len=*), parameter :: orientations(4) = [character(len=9) :: "lab", "auto", "molecular", &
         & "molecular"]
      character(len=*), parameter :: groups(4) = [character(len=4) :: "off", "off", "off", "auto"]
      character(len=40) :: label
      integer :: im, iat, irot

      call new_context(ctx, nthreads=0, verbosity=0)
      do im = 1, size(names)
         call structure_of(trim(names(im)), mol)
         settings%orientation = orientations(im)
         settings%group = groups(im)
         write (label, '("NumSA ",a," (",a,", ",a,")")') trim(names(im)), trim(orientations(im)), &
            & trim(groups(im))
         call updated_cavity(kind_numsa, mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, trim(label)//": "//err%message)
            return
         end if
         select type (cav)
         type is (cavity_type_numsa)
            call check(error, abs(sum(cav%asph) - cav%total_area) <= 1.0e-12_wp*cav%total_area, &
               & trim(label)//": sphere areas do not add up")
            if (allocated(error)) return
            if (cav%plain_grid()) then
               call check(error,.not. ((.not. cav%grid%plain) .or. allocated(cav%frame_torque)), &
                  & trim(label)//": the lab orientation built rotated nodes or a torque")
               if (allocated(error)) return
               cycle
            end if
            call check(error, (.not. cav%grid%plain) .and. allocated(cav%frame_torque), &
               & trim(label)//": no rotated nodes or torque")
            if (allocated(error)) return
            associate (gm_ => symmetry_of(cav))
               call check(error, all(shape(cav%grid%node) == [3, cav%param%num_leb*model_nrot(gm_)]) .and. &
                  & all(shape(cav%frame_torque) == [3, mol%nat]), &
                  & trim(label)//": not one rotated copy per distinct rotation, or no torque per sphere")
            end associate
            if (allocated(error)) return
            do iat = 1, mol%nat
               associate (gm_ => symmetry_of(cav))
                  irot = model_rot_index(gm_, iat)
               end associate
               associate (gm_ => symmetry_of(cav))
                  call check(error, maxval(abs(cav%grid%node(:, (irot - 1)*size(cav%ang_grid, 2) + 1: &
                     & irot*size(cav%ang_grid, 2)) - matmul(model_rotation(gm_, iat), cav%ang_grid))) < 1.0e-14_wp, &
                     & trim(label)//": nodes of a sphere are not R u")
               end associate
               if (allocated(error)) return
            end do
         class default
            call test_failed(error, trim(label)//": not a NumSA cavity")
            return
         end select
      end do
   end subroutine test_numsa_nodes

   !> NumSA frame gradient at the eigen level: an asymmetric molecule
   !>
   !> - total area (`area_grad`) and a sphere-weighted area (`dsdr`) against
   !>   central differences along every Cartesian coordinate, at the accuracy
   !>   of the NumSA suite; neither carries a net force or torque
   !> - without the frame term both miss the differences
   !>
   !> @param[out] error Error handle
   subroutine test_numsa_gradient_eigen(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_numsa_frame_cartesian("c1", frame_level_eigen, error)
   end subroutine test_numsa_gradient_eigen

   !> NumSA frame gradient at the plane level: two equal second moments
   !>
   !> @param[out] error Error handle
   subroutine test_numsa_gradient_plane(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_numsa_frame_cartesian("c1_plane", frame_level_plane, error)
   end subroutine test_numsa_gradient_plane

   !> NumSA frame gradient of a molecule without a group, every coordinate
   !>
   !> @param[in]  name  Molecule
   !> @param[in]  level Frame level the reference must have
   !> @param[out] error Test failure
   subroutine check_numsa_frame_cartesian(name, level, error)
      !> Molecule
      character(len=*), intent(in) :: name
      !> Frame level
      integer, intent(in) :: level
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      real(wp), allocatable :: grads(:, :, :)

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_frame_fixture(name, "off", fix, ctx, error, kind=kind_numsa)
      if (allocated(error)) return
      call check(error, fix%level, level, fix%name//": frame not at the expected level")
      if (allocated(error)) return
      call check_numsa_frame_gradient(fix, ctx, cartesian_directions(fix%mol%nat), error, grads)
      if (allocated(error)) return
      call check_no_net_force(fix, grads(:, :, 1), "area_grad", error)
      if (allocated(error)) return
      call check_no_net_force(fix, grads(:, :, 2), "weighted dsdr", error)
      if (allocated(error)) return

      ! The frame term is needed
      call check_numsa_frame_gradient(fix, ctx, cartesian_directions(fix%mol%nat), error, grads, &
         & frame_off=.true.)
   end subroutine check_numsa_frame_cartesian

   !> NumSA frame gradient with a group: directions that keep it
   !>
   !> - ammonia (C3v): the frame does not turn along such directions
   !> - HOCl (Cs): every sphere turns about the plane normal with the frame,
   !>   so without the frame term both gradients miss
   !>
   !> @param[out] error Error handle
   subroutine test_numsa_gradient_symmetric(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(2) = [character(len=7) :: "ammonia", "hocl"]
      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: grads(:, :, :), dirs(:, :, :)
      integer :: im

      call new_context(ctx, nthreads=0, verbosity=0)
      do im = 1, size(names)
         call new_frame_fixture(trim(names(im)), "auto", fix, ctx, error, kind=kind_numsa)
         if (allocated(error)) return
         call new_frame_cavity(fix, fix%mol, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         associate (gm_ => symmetry_of(cav))
            dirs = group_directions(gm_%orient%ops, gm_%orient%perm, 3)
         end associate
         call check_numsa_frame_gradient(fix, ctx, dirs, error, grads)
         if (allocated(error)) return
         if (trim(names(im)) == "hocl") then
            call check_numsa_frame_gradient(fix, ctx, dirs, error, grads, frame_off=.true.)
            if (allocated(error)) return
         end if
      end do
   end subroutine test_numsa_gradient_symmetric

   !> Gradients along rigid motions, no group: the frame turns with the molecule
   !>
   !> @param[out] error Error handle
   subroutine test_rigid_rotation_frame(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(3) = [character(len=5) :: "c1", "hocl", "water"]
      integer :: im

      do im = 1, size(names)
         call check_rigid_paths(kind_iswig, trim(names(im)), "off", error)
         if (allocated(error)) return
         call check_rigid_paths(kind_numsa, trim(names(im)), "off", error)
         if (allocated(error)) return
      end do
      call check_rigid_paths(kind_iswig, "co2", "off", error)
   end subroutine test_rigid_rotation_frame

   !> Gradients along rigid motions with a group: the group turns with the
   !> molecule, so every sphere turns rigidly
   !>
   !> - quantities that are not totally symmetric: the CPCM energy of fixed
   !>   off-nucleus charges (reverse and forward), sphere-weighted NumSA areas
   !>   through `dsdr`
   !> - water (C2v and C2), ammonia (C3v), HOCl (Cs), benzene (D6h and its
   !>   subgroup C6h, whose rotation about the axis keeps the group), methane
   !>   (Td and C3v), CO2 (linear)
   !> - no subgroup with conjugate copies in the full group that the fixed
   !>   charges tell apart (D2h in D6h): which copy is used is a tie the
   !>   molecule cannot break, and it may change along a rotation
   !>
   !> @param[out] error Error handle
   subroutine test_rigid_rotation_group(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(6) = [character(len=7) :: "water", "ammonia", "hocl", "benzene", &
         & "methane", "co2"]
      integer :: im

      do im = 1, size(names)
         call check_rigid_paths(kind_iswig, trim(names(im)), "auto", error)
         if (allocated(error)) return
         call check_rigid_paths(kind_numsa, trim(names(im)), "auto", error)
         if (allocated(error)) return
      end do
      ! Subgroups, rotations only: a displacement that keeps the subgroup
      ! alone lowers the symmetry of the geometry and so leaves the branch
      call check_rigid_paths(kind_iswig, "benzene", "C6h", error, npath=4)
      if (allocated(error)) return
      call check_rigid_paths(kind_numsa, "benzene", "C6h", error, npath=4)
      if (allocated(error)) return
      call check_rigid_paths(kind_iswig, "benzene", "C2h", error, tilted=.true.)
      if (allocated(error)) return
      call check_rigid_paths(kind_numsa, "benzene", "C2h", error, tilted=.true.)
      if (allocated(error)) return
      call check_rigid_paths(kind_iswig, "water", "C2", error, npath=4)
      if (allocated(error)) return
      call check_rigid_paths(kind_numsa, "methane", "C3v", error, npath=4)
   end subroutine test_rigid_rotation_group

   !> DROP gradients along rigid motions with a group
   !>
   !> @param[out] error Error handle
   subroutine test_drop_rigid_rotation_group(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_rigid_paths(kind_drop, "water", "auto", error)
      if (allocated(error)) return
      call check_rigid_paths(kind_drop, "ammonia", "auto", error)
   end subroutine test_drop_rigid_rotation_group

   !> Analytic gradients against differences along finite rigid motions
   !>
   !> Paths x(t) = c + t s + Q(t u) (x + t a d - c) with exact rotations Q
   !> about the centroid c, which keep the group at every stencil point:
   !> rotations about x, y, z and a generic axis, and a rotation composed
   !> with a translation s and a displacement d that keeps the group
   !>
   !> - iSwiG and DROP: every gradient of [[frame_gradients]]
   !> - NumSA: total area and the sphere-weighted area through `dsdr`
   !> - step `fix%step` in radians, tolerances of the frame-gradient checks
   !> - on top of the absolute tolerance, the rounding of the differences:
   !>   `FD4_WEIGHT_SUM*FD_ROUND*epsilon*max|value|/step`. Area and volume of
   !>   a large cavity do not change along a rotation, and their differences
   !>   are rounding alone
   !> - iSwiG and DROP: fails when the number of points differs between the
   !>   reference and a stencil point, where the differences would cross a
   !>   jump of the value and miss a correct gradient
   !> - `keep_all` (iSwiG): no point of weight is dropped, `cut_f` at 1e-20.
   !>   Where spheres turn against each other, points sweep the switching
   !>   region and cross a larger cutoff within the stencil: a jump of the
   !>   differenced value, not of the gradient
   !>
   !> @param[in]  kind     Cavity kind
   !> @param[in]  name     Molecule ("co2" or one of [[structure_of]])
   !> @param[in]  group    Symmetry group setting
   !> @param[out] error    Test failure
   !> @param[in]  npath    Number of paths (default all five)
   !> @param[in]  tilted   Start from a generic pose
   !> @param[in]  keep_all Whether no point of weight is dropped (default off)
   !> @param[in]  step     Step, rad (default that of the fixture)
   subroutine check_rigid_paths(kind, name, group, error, npath, tilted, keep_all, step)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Molecule
      character(len=*), intent(in) :: name
      !> Symmetry group setting
      character(len=*), intent(in) :: group
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Number of paths
      integer, intent(in), optional :: npath
      !> Start from a generic pose
      logical, intent(in), optional :: tilted
      !> Whether no point of weight is dropped
      logical, intent(in), optional :: keep_all
      !> Step
      real(wp), intent(in), optional :: step

      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      type(structure_type) :: mol, moved
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: grads(:, :, :), keep(:, :, :), tangent(:, :), xyz(:, :), vals(:, :), numeric(:)
      real(wp), allocatable :: w(:)
      real(wp) :: axis(3), shift(3), center(3), q(3, 3), theta, amp, analytic, thr_abs, thr_rel
      character(len=240) :: context
      integer :: nat, iat, ipath, k, igrad, iq, nq, np, ngrid(4)

      np = 5
      if (present(npath)) np = npath
      call new_context(ctx, nthreads=0, verbosity=0)
      if (name == "co2") then
         call new(mol, [6, 8, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.16_wp, &
            & 0.0_wp, 0.0_wp, -1.16_wp], [3, 3])*aatoau)
         ! Off the lab axes, where the lab reference of the line is generic
         call rigid_motion(mol, test_rotation(1), [0.0_wp, 0.0_wp, 0.0_wp], [1, 2, 3], moved)
         call new_frame_fixture(name, group, fix, ctx, error, kind=kind, mol=moved)
      else if (present(tilted)) then
         call structure_of(name, mol)
         call rigid_motion(mol, test_rotation(1), [0.0_wp, 0.0_wp, 0.0_wp], [(iat, iat=1, mol%nat)], moved)
         call new_frame_fixture(name, group, fix, ctx, error, kind=kind, mol=moved)
      else
         call new_frame_fixture(name, group, fix, ctx, error, kind=kind)
      end if
      if (allocated(error)) return
      fix%name = fix%name//" ("//group//")"
      ! Numbering of the nodes aside: among the frames that give one sphere
      ! the numbering follows a lab direction, which a rotation may cross
      fix%pinned = .false.
      if (present(keep_all)) then
         if (keep_all) then
            ! A point that crosses it changes nothing to rounding; at zero
            ! points with f = 0 stay, where the CPCM gradient is not finite.
            ! Unpinned: the surviving points of the fixture are not used
            fix%cut_f = 1.0e-20_wp
         end if
      end if
      if (present(step)) fix%step = step
      nat = fix%mol%nat
      center = sum(fix%mol%xyz, dim=2)/real(nat, wp)
      call new_frame_cavity(fix, fix%mol, ctx, cav, err)
      if (.not. allocated(err)) call cav%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, fix%name//": "//err%message)
         return
      end if
      associate (gm_ => symmetry_of(cav))
         keep = group_directions(gm_%orient%ops, gm_%orient%perm, 1)
      end associate

      if (kind == kind_numsa) then
         nq = 2
         thr_abs = FD_ABS_NUMSA
         thr_rel = 0.0_wp
         w = numsa_sphere_weights(nat)
         allocate (grads(3, nat, 2), source=0.0_wp)
         select type (cav)
         type is (cavity_type_numsa)
            grads(:, :, 1) = cav%area_grad
            do iat = 1, nat
               grads(:, :, 2) = grads(:, :, 2) + w(iat)*cav%dsdr(:, :, iat)
            end do
         end select
      else
         nq = nquant
         ! Area and volume do not change: their differences are rounding
         thr_abs = 4.0_wp*FD_ABS
         thr_rel = FD_REL
         call frame_gradients(fix, ctx, grads, error)
         if (allocated(error)) return
      end if
      allocate (tangent(3, nat), xyz(3, nat), vals(nq, 4), numeric(nq))

      do ipath = 1, np
         shift = 0.0_wp
         amp = 0.0_wp
         select case (ipath)
         case (1, 2, 3)
            axis = 0.0_wp
            axis(ipath) = 1.0_wp
         case (4)
            axis = [0.48_wp, -0.6_wp, 0.64_wp]
         case default
            axis = [-0.36_wp, 0.8_wp, 0.48_wp]
            shift = [0.3_wp, -0.2_wp, 0.5_wp]
            amp = 0.4_wp
         end select
         do iat = 1, nat
            tangent(:, iat) = shift + cross_product(axis, fix%mol%xyz(:, iat) - center) + amp*keep(:, iat, 1)
         end do
         do k = 1, 4
            theta = fd4_offsets(k)*fix%step
            q = axis_rotation(axis, theta)
            do iat = 1, nat
               xyz(:, iat) = center + theta*shift &
                  & + matmul(q, fix%mol%xyz(:, iat) + theta*amp*keep(:, iat, 1) - center)
            end do
            if (kind == kind_numsa) then
               call numsa_frame_values(fix, xyz, ctx, w, vals(:, k), error)
            else
               call frame_values(fix, xyz, ctx, vals(:, k), error, ngrid=ngrid(k))
            end if
            if (allocated(error)) return
         end do
         if (kind /= kind_numsa) then
            ! Unpinned, the surviving points are not compared one by one
            write (context, '(a,", path ",i0,": ",i0," points, ",4(1x,i0)," at the stencil points")') &
               & fix%name, ipath, cav%ngrid, ngrid
            call check(error, all(ngrid == cav%ngrid), &
               & "a stencil point of a rigid motion changes the number of points: differences across a jump", &
               & more=trim(context))
            if (allocated(error)) return
         end if
         do iq = 1, nq
            call fd4_scalar(vals(iq, 1), vals(iq, 2), vals(iq, 3), vals(iq, 4), fix%step, numeric(iq), error)
            if (allocated(error)) return
         end do
         do igrad = 1, size(grads, 3)
            if (kind == kind_numsa) then
               iq = igrad
               write (context, '(a,", gradient ",i0,", path ",i0)') fix%name, igrad, ipath
            else
               iq = grad_quantity(igrad)
               if (iq == q_adjoint .and. .not. fix%pinned) cycle
               write (context, '(a,", ",a,", path ",i0)') fix%name, trim(grad_name(igrad)), ipath
            end if
            analytic = sum(grads(:, :, igrad)*tangent)
            write (context(len_trim(context) + 1:), '(": ",es14.6," vs FD ",es14.6)') analytic, numeric(iq)
            ! Rounding of the differenced values on top of the tolerance
            call check(error, analytic, numeric(iq), thr_abs=thr_abs &
               & + FD4_WEIGHT_SUM*FD_ROUND*epsilon(1.0_wp)*maxval(abs(vals(iq, :)))/fix%step, thr_rel=thr_rel, &
               & message="gradient along a rigid motion disagrees with central differences", &
               & more=trim(context))
            if (allocated(error)) return
         end do
      end do
   end subroutine check_rigid_paths

   !> Rotation by an angle about a unit axis (Rodrigues)
   !>
   !> @param[in] axis  Unit axis
   !> @param[in] theta Angle, rad
   pure function axis_rotation(axis, theta) result(q)
      !> Unit axis
      real(wp), intent(in) :: axis(3)
      !> Angle
      real(wp), intent(in) :: theta
      !> Rotation
      real(wp) :: q(3, 3)

      real(wp) :: k(3, 3)

      k = reshape([0.0_wp, axis(3), -axis(2), -axis(3), 0.0_wp, axis(1), axis(2), -axis(1), 0.0_wp], [3, 3])
      q = identity3 + sin(theta)*k + (1.0_wp - cos(theta))*matmul(k, k)
   end function axis_rotation

   !> NumSA area gradients of a fixture against central differences
   !>
   !> - total area through `area_grad`, and sum_a w_a A_a through `dsdr` with
   !>   pseudo-random sphere weights w_a ([[numsa_sphere_weights]])
   !> - step `fix%step`, tolerance `FD_ABS_NUMSA` (absolute), as the NumSA suite
   !> - with `frame_off`, the frame term is taken back out of both, mapped from
   !>   `frame_torque` by `add_gradient` (the reverse map, independent
   !>   of the forward `jacobian` the cavity uses): both must then miss
   !>   in some direction
   !>
   !> @param[in]  fix       Fixture of kind NumSA
   !> @param[in]  ctx       Run context
   !> @param[in]  dirs      Displacement directions (3, nat, ndir)
   !> @param[out] error     Test failure
   !> @param[out] grads     Analytic gradients (3, nat, 2): area, weighted area
   !> @param[in]  frame_off Leave the frame term out and expect failures
   !> @param[in]  thr       Absolute tolerance (default that of the NumSA suite)
   subroutine check_numsa_frame_gradient(fix, ctx, dirs, error, grads, frame_off, thr)
      !> Fixture
      type(frame_fixture), intent(in) :: fix
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Displacement directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Analytic gradients
      real(wp), allocatable, intent(out) :: grads(:, :, :)
      !> Leave the frame term out
      logical, intent(in), optional :: frame_off
      !> Absolute tolerance
      real(wp), intent(in), optional :: thr

      character(len=*), parameter :: names(2) = [character(len=13) :: "area_grad", "weighted dsdr"]
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: frame(:, :)
      !> Sphere weights of the weighted area
      real(wp) :: w(fix%mol%nat)
      real(wp) :: vals(2, 4), numeric(2), analytic, thr_abs
      character(len=200) :: context
      logical :: off, missed(2)
      integer :: nat, iat, idir, k, iq

      off = .false.
      if (present(frame_off)) off = frame_off
      thr_abs = FD_ABS_NUMSA
      if (present(thr)) thr_abs = thr
      nat = fix%mol%nat
      w = numsa_sphere_weights(nat)
      allocate (grads(3, nat, 2), source=0.0_wp)
      call new_frame_cavity(fix, fix%mol, ctx, cav, err)
      if (.not. allocated(err)) call cav%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, fix%name//": "//err%message)
         return
      end if
      select type (cav)
      type is (cavity_type_numsa)
         grads(:, :, 1) = cav%area_grad
         do iat = 1, nat
            grads(:, :, 2) = grads(:, :, 2) + w(iat)*cav%dsdr(:, :, iat)
         end do
         if (off) then
            allocate (frame(3, nat), source=0.0_wp)
            associate (gm_ => symmetry_of(cav))
               call model_add_torque(gm_, cav%frame_torque, frame)
            end associate
            grads(:, :, 1) = grads(:, :, 1) - frame
            frame = 0.0_wp
            associate (gm_ => symmetry_of(cav))
               call model_add_torque(gm_, cav%frame_torque*spread(w, 1, 3), frame)
            end associate
            grads(:, :, 2) = grads(:, :, 2) - frame
         end if
      class default
         call test_failed(error, fix%name//": not a NumSA cavity")
         return
      end select
      do iq = 1, 2
         call check(error, maxval(abs(grads(:, :, iq))) > 1.0e-8_wp, fix%name//": "//trim(names(iq))//" vanishes")
         if (allocated(error)) return
      end do

      missed = .false.
      do idir = 1, size(dirs, 3)
         do k = 1, 4
            call numsa_frame_values(fix, fix%mol%xyz + fd4_offsets(k)*fix%step*dirs(:, :, idir), ctx, w, &
               & vals(:, k), error)
            if (allocated(error)) return
         end do
         do iq = 1, 2
            call fd4_scalar(vals(iq, 1), vals(iq, 2), vals(iq, 3), vals(iq, 4), fix%step, numeric(iq), error)
            if (allocated(error)) return
            analytic = sum(grads(:, :, iq)*dirs(:, :, idir))
            if (off) then
               missed(iq) = missed(iq) .or. abs(analytic - numeric(iq)) > FD_ABS_NUMSA
               cycle
            end if
            write (context, '(a,", ",a,", direction ",i0,": ",es14.6," vs FD ",es14.6)') fix%name, &
               & trim(names(iq)), idir, analytic, numeric(iq)
            call check(error, analytic, numeric(iq), thr=thr_abs, &
               & message="frame gradient disagrees with central differences", more=trim(context))
            if (allocated(error)) return
         end do
      end do
      if (.not. off) return
      do iq = 1, 2
         call check(error, missed(iq), fix%name//": "//trim(names(iq))// &
            & " agrees with central differences without the frame term")
         if (allocated(error)) return
      end do
   end subroutine check_numsa_frame_gradient

   !> NumSA areas of a fixture at a displaced geometry
   !>
   !> Fails when a stencil point leaves the reference branch
   !> ([[check_stencil_branch]])
   !>
   !> @param[in]  fix   Fixture of kind NumSA
   !> @param[in]  xyz   Displaced coordinates (3, nat)
   !> @param[in]  ctx   Run context
   !> @param[in]  w     Sphere weights (nat)
   !> @param[out] vals  Total area, sum_a w_a A_a
   !> @param[out] error Test failure
   subroutine numsa_frame_values(fix, xyz, ctx, w, vals, error)
      !> Fixture
      type(frame_fixture), intent(in) :: fix
      !> Displaced coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Sphere weights
      real(wp), intent(in) :: w(:)
      !> Total area, weighted area
      real(wp), intent(out) :: vals(2)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err

      vals = 0.0_wp
      mol = fix%mol
      mol%xyz = xyz
      call new_frame_cavity(fix, mol, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, fix%name//": "//err%message)
         return
      end if
      call check_stencil_branch(fix, cav, error)
      if (allocated(error)) return
      vals(1) = cav%total_area
      vals(2) = sum(w*cav%asph)
   end subroutine numsa_frame_values

   !> Pseudo-random sphere weights of the NumSA weighted area
   !>
   !> @param[in] nat Number of spheres
   pure function numsa_sphere_weights(nat) result(w)
      !> Number of spheres
      integer, intent(in) :: nat
      !> Weights (nat), between 0.5 and 1.5
      real(wp) :: w(nat)

      integer :: iat

      do iat = 1, nat
         w(iat) = 1.0_wp + 0.5_wp*sin(1.7_wp*real(iat, wp) + 0.3_wp)
      end do
   end function numsa_sphere_weights

   !* ================================================================================= *!
   !*                              iSwiG orientation helpers                            *!
   !* ================================================================================= *!

   !> Construct and update an iSwiG cavity with a symmetry model
   !>
   !> The model is built from the settings and computed by the cavity
   !>
   !> @param[in]  mol      Structure
   !> @param[in]  settings Symmetry settings
   !> @param[in]  ctx      Run context
   !> @param[out] cav      Updated cavity
   !> @param[out] err      Construction or update error
   !> @param[in]  num_leb  Lebedev points per sphere (default 110)
   !> @param[in]  cut_f    Switching cutoff (default the compiled one)
   subroutine new_iswig(mol, settings, ctx, cav, err, num_leb, cut_f)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Symmetry settings
      type(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: settings
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Cavity
      type(cavity_type_iswig), allocatable, intent(out) :: cav
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: err
      !> Lebedev points per sphere
      integer, intent(in), optional :: num_leb
      !> Switching cutoff
      real(wp), intent(in), optional :: cut_f

      type(moist_cavity_iswig_parameters_type) :: param
      type(grid_model_symmetry_type) :: model

      if (present(num_leb)) param%num_leb = num_leb
      if (present(cut_f)) param%cut_f = cut_f
      call new_grid_model_symmetry(model, err, param=settings)
      if (allocated(err)) return
      allocate (cav)
      call new_cavity_iswig(cav, default_cpcm_radii(), err, param=param, ctx=ctx, grid_model=model)
      if (allocated(err)) return
      call cav%update(mol, err)
   end subroutine new_iswig

   !* ================================================================================= *!
   !*                                 Frame thresholds                                  *!
   !* ================================================================================= *!

   !> Structures for the single-frame checks, in a generic position
   !>
   !> - `ammonia`, `benzene`, `methane`: exact tops
   !> - `rod`: thick helix 36 bohr long, no line: one eigenvector frame
   !> - `ammonia_near`: slightly distorted, eigen-gap near its threshold
   !> - `ring_weak`: planar rings with one weak in-plane harmonic
   !> - `co2_near`: bent by 0.01 rad, near the line threshold
   !> - `co2_bent`: bent by 0.3 rad
   !>
   !> @param[in]  name Case
   !> @param[out] mol  Structure
   subroutine frame_structure(name, mol)
      !> Case
      character(len=*), intent(in) :: name
      !> Structure
      type(structure_type), intent(out) :: mol

      type(structure_type) :: base
      real(wp) :: xyz(3, 3), rod(3, 24), t
      integer :: k

      select case (name)
      case ("ammonia", "benzene", "methane", "water")
         call structure_of(name, base)
      case ("rod")
         do k = 1, 24
            t = 1.7_wp*real(k, wp)
            rod(:, k) = [0.63_wp*cos(t), 0.45_wp*sin(t) + 0.02_wp*sin(0.3_wp*real(k*k, wp)), 1.5_wp*real(k, wp)]
         end do
         call new(base, [(6, k=1, 24)], rod)
      case ("ring_weak")
         ! Planar 13-ring (no harmonic up to order 12) around a small
         ! triangle: a lone third harmonic of relative size 5e-4
         block
            real(wp) :: ring(3, 16), rho, phi

            rho = 2.0_wp*(5.0e-4_wp*13.0_wp/(3.0_wp*(1.0_wp - 5.0e-4_wp)))**(1.0_wp/3.0_wp)
            do k = 1, 13
               phi = 2.0_wp*pi*real(k - 1, wp)/13.0_wp
               ring(:, k) = [2.0_wp*cos(phi), 2.0_wp*sin(phi), 0.0_wp]
            end do
            do k = 1, 3
               phi = 0.3_wp + 2.0_wp*pi*real(k - 1, wp)/3.0_wp
               ring(:, 13 + k) = [rho*cos(phi), rho*sin(phi), 0.0_wp]
            end do
            call new(base, [(6, k=1, 16)], ring)
         end block
      case ("ammonia_near")
         call structure_of("ammonia", base)
         base%xyz(:, 2) = base%xyz(:, 2) + [0.0_wp, 2.0e-3_wp, 1.0e-3_wp]
      case default
         t = merge(0.01_wp, 0.3_wp, name == "co2_near")
         xyz(:, 1) = 0.0_wp
         xyz(:, 2) = 2.2_wp*[sin(t), 0.0_wp, cos(t)]
         xyz(:, 3) = 2.2_wp*[sin(t), 0.0_wp, -cos(t)]
         call new(base, [6, 8, 8], xyz)
      end select
      call rigid_motion(base, rotation_matrix(0.5_wp, 1.2_wp, -0.4_wp), [0.4_wp, -1.1_wp, 2.2_wp], &
         & [(k, k=1, base%nat)], mol)
   end subroutine frame_structure

   !> Total area of an iSwiG or NumSA cavity with its nuclear gradient
   !>
   !> @param[in]  kind     Cavity kind, iSwiG or NumSA
   !> @param[in]  mol      Structure
   !> @param[in]  settings Symmetry settings
   !> @param[in]  ctx      Run context
   !> @param[out] area     Total area
   !> @param[out] error    Test failure
   !> @param[out] grad     Forward gradient of the area (3, nat)
   subroutine frame_area(kind, mol, settings, ctx, area, error, grad)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Symmetry settings
      type(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: settings
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Total area
      real(wp), intent(out) :: area
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Forward gradient
      real(wp), intent(out), optional :: grad(:, :)

      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err

      area = 0.0_wp
      call updated_cavity(kind, mol, settings, ctx, cav, err)
      if (.not. allocated(err) .and. present(grad)) call cav%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, trim(kind_name(kind))//": "//err%message)
         return
      end if
      area = cav%total_area
      if (present(grad)) then
         select type (cav)
         type is (cavity_type_iswig)
            grad = cav%area_grad
         type is (cavity_type_numsa)
            grad = cav%area_grad
         class default
            call test_failed(error, "frame_area: unsupported cavity kind")
            return
         end select
      end if
   end subroutine frame_area

   !> iSwiG and NumSA on tops and a rod without a group: one grid per sphere
   !>
   !> - one component per sphere of weight one
   !> - the area follows a rotation, shift and renumbering of the molecule,
   !>   also where the candidates of a harmonic tie by symmetry
   !>
   !> @param[out] error Error handle
   subroutine test_single_frame_tops(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(4) = [character(len=8) :: "ammonia", "benzene", "methane", "rod"]
      type(moist_context_type), target :: ctx
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(structure_type) :: mol, moved
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp) :: area, area_moved
      integer :: im, kind, k

      call new_context(ctx, nthreads=0, verbosity=0)
      settings%orientation = "molecular"
      do im = 1, size(names)
         call frame_structure(trim(names(im)), mol)
         call rigid_motion(mol, rotation_matrix(1.1_wp, -0.7_wp, 2.3_wp), [3.0_wp, -1.0_wp, 0.5_wp], &
            & [(k, k=mol%nat, 1, -1)], moved)
         do kind = kind_iswig, kind_numsa
            call updated_cavity(kind, mol, settings, ctx, cav, err)
            if (allocated(err)) then
               call test_failed(error, trim(names(im))//": "//err%message)
               return
            end if
            call check(error, cav%grid%ncomp == mol%nat .and. all(cav%grid%comp_weight == 1.0_wp) .and. &
               & .not. symmetry_near(cav), trim(names(im))//": not one plain grid per sphere")
            if (allocated(error)) return
            deallocate (cav)
            call frame_area(kind, mol, settings, ctx, area, error)
            if (allocated(error)) return
            call frame_area(kind, moved, settings, ctx, area_moved, error)
            if (allocated(error)) return
            call check(error, abs(area_moved - area) < 1.0e-10_wp*area, &
               & trim(names(im))//" ("//trim(kind_name(kind))//"): area does not follow the molecule")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_single_frame_tops

   !> Whether the frame of a cavity's symmetry model is near a switch or tied
   !>
   !> @param[in] cav Updated cavity
   function symmetry_near(cav) result(near)
      !> Cavity
      class(cavity_type), intent(in) :: cav
      !> Whether any flag of the frame is set
      logical :: near

      associate (sym => symmetry_of(cav))
         near = sym%near_switch() .or. sym%tied_frames()
      end associate
   end function symmetry_near

   !> A frame near a switch is flagged, printed and warned about
   !>
   !> - slightly distorted ammonia (eigen-gap) and slightly bent CO2 (line)
   !>   warn; water, exact ammonia and bent CO2 do not
   !> - with the group of exact ammonia in use nothing is flagged
   !> - the lab orientation never warns
   !>
   !> @param[out] error Error handle
   subroutine test_near_switch_warning(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(6) = [character(len=12) :: "ammonia_near", "co2_near", &
         & "water", "ammonia", "co2_bent", "ammonia"]
      !> Expected flags: nearly degenerate eigenvalues, nearly linear
      logical, parameter :: gap(6) = [.true., .false., .false., .false., .false., .false.]
      logical, parameter :: line(6) = [.false., .true., .false., .false., .false., .false.]
      character(len=*), parameter :: fields(4) = [character(len=24) :: "symmetry_near_gap", &
         & "symmetry_near_linear", "symmetry_near_harmonic", "symmetry_tied_frames"]
      logical :: flags(4)
      integer :: k
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(grid_model_symmetry_type) :: model
      class(radius_type), allocatable :: radius_model
      type(mctc_error), allocatable :: err
      type(field_query_type) :: query
      character(len=line_len) :: lines(64), out(64)
      integer :: im, unit, nline, nout, ctx_unit
      logical :: warned

      do im = 1, size(names) + 1
         settings = moist_cavity_grid_model_symmetry_parameters_type()
         settings%orientation = "molecular"
         if (im == 6) settings%group = "auto"
         if (im == 7) settings%orientation = "lab"
         call frame_structure(trim(names(min(im, 6))), mol)
         if (im == 7) call frame_structure("ammonia_near", mol)
         open (newunit=ctx_unit, status="scratch", action="readwrite")
         call new_context(ctx, nthreads=0, verbosity=0, unit=ctx_unit)
         call host_model(mol, settings, ctx, model, radius_model, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            close (ctx_unit)
            return
         end if
         call read_printout(ctx_unit, out, nout)
         close (ctx_unit)
         call ctx%delete()

         flags = [model%near_gap(), model%near_linear(), model%near_harmonic(), model%tied_frames()]
         call check(error, all(flags .eqv. [gap(min(im, 6)) .and. im /= 7, line(min(im, 6)) .and. im /= 7, &
            & .false., .false.]), trim(names(min(im, 6)))//": unexpected near-switch flags")
         if (allocated(error)) return
         call check(error, model%near_switch() .eqv. any(flags(1:3)), &
            & trim(names(min(im, 6)))//": near_switch differs from its three flags")
         if (allocated(error)) return
         warned = .false.
         if (nout >= 1) warned = any(index(out(1:nout), "[Warning] Symmetry:") > 0 .and. &
            & index(out(1:nout), "close to a switch") > 0)
         call check(error, warned .eqv. model%near_switch(), &
            & trim(names(min(im, 6)))//": warning does not follow the flags")
         if (allocated(error)) return
         if (im == 1) then
            call check(error, any(index(out(1:nout), "nearly degenerate second-moment eigenvalues") > 0), &
               & "distorted ammonia: warning does not name the eigenvalues")
            if (allocated(error)) return
         else if (im == 2) then
            call check(error, any(index(out(1:nout), "nearly linear molecule") > 0), &
               & "slightly bent CO2: warning does not name the line")
            if (allocated(error)) return
         end if

         open (newunit=unit, status="scratch", action="readwrite")
         call model%print(unit, 1)
         call read_printout(unit, lines, nline)
         close (unit)
         call check(error, any(index(lines(1:nline), "Near a switch") > 0) .eqv. model%near_switch(), &
            & trim(names(min(im, 6)))//": report does not follow the flags")
         if (allocated(error)) return

         do k = 1, size(fields)
            call query%fetch(trim(fields(k)))
            call model%list_fields(query)
            call check(error, query%found .and. allocated(query%lvals), trim(fields(k))//": field missing")
            if (allocated(error)) return
            call check(error, size(query%lvals) == 1 .and. all(query%lvals .eqv. flags(k)), &
               & trim(fields(k))//": field differs from the flag")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_near_switch_warning

   !> Frames that tie exactly and the species tell apart are reported
   !>
   !> The frame follows the atom positions alone. Where these have a
   !> symmetry the cube lacks and the species do not share it, the frames
   !> that tie give different cavities, and the pose decides
   !>
   !> - a regular pentagon and a regular heptagon with one N, an icosahedron
   !>   with one N: `tied_frames`, a warning, a line in the
   !>   report; also with the detected group in use
   !> - benzene, a ring of alternating B and N, ammonia, methane, the plain
   !>   icosahedron and the ci structure: tied frames, if any, show the same
   !>   molecule; no flag, and the iSwiG area follows a rotation, shift and
   !>   renumbering
   !> - no flag in the lab orientation
   !>
   !> @param[out] error Error handle
   subroutine test_tied_frames_warning(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ncase = 9
      character(len=*), parameter :: names(ncase) = [character(len=11) :: "pentagon_n", "heptagon_n", &
         & "icosa_n", "benzene", "ring_bn", "ammonia", "methane", "icosahedron", "ci"]
      logical, parameter :: tied(ncase) = [.true., .true., .true., .false., .false., .false., .false., &
         & .false., .false.]
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, moved
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(grid_model_symmetry_type) :: model
      class(radius_type), allocatable :: radius_model
      type(mctc_error), allocatable :: err
      character(len=line_len) :: lines(64), out(64)
      real(wp) :: area0, area
      integer :: im, ig, it, unit, nline, nout, ctx_unit
      logical :: warned

      do im = 1, ncase
         call tied_structure(trim(names(im)), mol)
         do ig = 1, 3
            settings = moist_cavity_grid_model_symmetry_parameters_type()
            settings%orientation = "molecular"
            if (ig == 2) settings%group = "auto"
            if (ig == 3) settings%orientation = "lab"
            open (newunit=ctx_unit, status="scratch", action="readwrite")
            call new_context(ctx, nthreads=0, verbosity=0, unit=ctx_unit)
            call host_model(mol, settings, ctx, model, radius_model, err)
            if (allocated(err)) then
               call test_failed(error, trim(names(im))//": "//err%message)
               close (ctx_unit)
               return
            end if
            call read_printout(ctx_unit, out, nout)
            close (ctx_unit)
            call ctx%delete()

            call check(error, model%tied_frames() .eqv. (tied(im) .and. ig < 3), &
               & trim(names(im))//": tied frames flagged against the species")
            if (allocated(error)) return
            warned = .false.
            if (nout >= 1) warned = any(index(out(1:nout), "[Warning] Symmetry:") > 0 .and. &
               & index(out(1:nout), "do not determine the grid orientation uniquely") > 0)
            call check(error, warned .eqv. (tied(im) .and. ig < 3), &
               & trim(names(im))//": warning does not follow the flag")
            if (allocated(error)) return
            open (newunit=unit, status="scratch", action="readwrite")
            call model%print(unit, 1)
            call read_printout(unit, lines, nline)
            close (unit)
            call check(error, any(index(lines(1:nline), "Tied frames") > 0) .eqv. (tied(im) .and. ig < 3), &
               & trim(names(im))//": report does not follow the flag")
            if (allocated(error)) return
         end do

         ! No flag: the pose does not matter
         if (tied(im)) cycle
         settings = moist_cavity_grid_model_symmetry_parameters_type()
         settings%orientation = "molecular"
         call new_context(ctx, nthreads=0, verbosity=0)
         do it = 0, 4
            if (it == 0) then
               moved = mol
            else
               call rigid_motion(mol, test_rotation(it), 5.0_wp*sin(real(it, wp) + [0.5_wp, 2.1_wp, 3.7_wp]), &
                  & test_permutation(mol%nat, it), moved)
            end if
            call frame_area(kind_iswig, moved, settings, ctx, area, error)
            if (allocated(error)) return
            if (it == 0) area0 = area
            call check(error, abs(area - area0) <= 1.0e-10_wp*area0, &
               & trim(names(im))//": area depends on the pose without a flag")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_tied_frames_warning

   !> Structures of the tied-frame test
   !>
   !> pentagon_n, heptagon_n: regular rings of C with one N; ring_bn: a
   !> regular hexagon of alternating B and N; icosa_n: the icosahedron with
   !> one N; any other name from [[structure_of]]
   !>
   !> @param[in]  name Structure
   !> @param[out] mol  Structure in bohr
   subroutine tied_structure(name, mol)
      !> Structure
      character(len=*), intent(in) :: name
      !> Structure
      type(structure_type), intent(out) :: mol

      type(structure_type) :: base
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: num(:)
      integer :: n, k

      select case (name)
      case ("pentagon_n", "heptagon_n", "ring_bn")
         n = 5
         if (name == "heptagon_n") n = 7
         if (name == "ring_bn") n = 6
         allocate (xyz(3, n), num(n))
         do k = 1, n
            xyz(:, k) = 1.3_wp*aatoau*[cos(2.0_wp*pi*real(k - 1, wp)/real(n, wp)), &
               & sin(2.0_wp*pi*real(k - 1, wp)/real(n, wp)), 0.0_wp]
         end do
         num = 6
         num(1) = 7
         if (name == "ring_bn") num = [(merge(5, 7, mod(k, 2) == 1), k=1, n)]
         call new_tilted(mol, num, xyz)
      case ("icosa_n")
         call structure_of("icosahedron", base)
         num = base%num(base%id)
         num(1) = 7
         call new(mol, num, base%xyz)
      case default
         call structure_of(name, mol)
      end select
   end subroutine tied_structure

   !> Thresholds of the frame choices are settings
   !>
   !> - a gap threshold far below the eigen-gap of distorted ammonia gives the
   !>   eigenvector frame without a flag; the default gives the harmonic frame
   !> - a threshold that is not positive is refused
   !>
   !> @param[out] error Error handle
   subroutine test_threshold_settings(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(grid_model_symmetry_type) :: model
      class(radius_type), allocatable :: radius_model
      type(mctc_error), allocatable :: err

      call new_context(ctx, nthreads=0, verbosity=0)
      settings%orientation = "molecular"
      call frame_structure("ammonia_near", mol)
      settings%gap_threshold = 1.0e-2_wp
      call host_model(mol, settings, ctx, model, radius_model, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, model%orient%frame%level == frame_level_plane, &
         & "distorted ammonia: a wide gap threshold does not give the harmonic frame")
      if (allocated(error)) return

      settings%gap_threshold = 1.0e-9_wp
      call host_model(mol, settings, ctx, model, radius_model, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, model%orient%frame%level == frame_level_eigen .and. .not. model%near_switch(), &
         & "distorted ammonia: a tiny gap threshold still switches or flags")
      if (allocated(error)) return

      settings%gap_threshold = -1.0_wp
      call new_grid_model_symmetry(model, err, param=settings)
      call check(error, allocated(err), "negative threshold accepted")
   end subroutine test_threshold_settings

   !> Frame of an isotropic second moment and of a single atom
   !>
   !> - methane: one `sphere` frame, its threefold axis on the body diagonal
   !>   of the grid, so the hydrogens lie on body diagonals
   !> - distorted methane, group off: iSwiG and NumSA area gradients match
   !>   central differences through the `sphere` frame
   !> - one atom: the identity, and a full sphere
   !>
   !> @param[out] error Error handle
   subroutine test_sphere_frame_and_atom(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: step = 1.0e-6_wp
      type(moist_context_type), target :: ctx
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(structure_type) :: mol, shifted, atom
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: grad(:, :), fd(:, :)
      real(wp) :: area, ap, am, u(3), rot(3, 3)
      integer :: kind, iat, ic

      call new_context(ctx, nthreads=0, verbosity=0)
      settings%orientation = "molecular"

      call frame_structure("methane", mol)
      call updated_cavity(kind_iswig, mol, settings, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, "methane: "//err%message)
         return
      end if
      associate (sym => symmetry_of(cav))
         call check(error, sym%orient%frame%level == frame_level_sphere, "methane: not a sphere-level frame")
         if (allocated(error)) return
         rot = model_rotation(sym, 1)
      end associate
      do iat = 2, mol%nat
         u = matmul(transpose(rot), mol%xyz(:, iat) - mol%xyz(:, 1))
         u = abs(u)/norm2(u)
         call check(error, maxval(abs(u - 1.0_wp/sqrt(3.0_wp))) < 1.0e-10_wp, &
            & "methane: a hydrogen is off the body diagonals of the grid")
         if (allocated(error)) return
      end do
      deallocate (cav)

      ! A distortion well inside the sphere level, away from its ties
      mol%xyz(:, 2) = mol%xyz(:, 2) + [3.0e-5_wp, -1.0e-5_wp, 2.0e-5_wp]
      mol%xyz(:, 4) = mol%xyz(:, 4) + [-1.0e-5_wp, 2.0e-5_wp, 1.0e-5_wp]
      allocate (grad(3, mol%nat), fd(3, mol%nat))
      do kind = kind_iswig, kind_numsa
         call frame_area(kind, mol, settings, ctx, area, error, grad=grad)
         if (allocated(error)) return
         call updated_cavity(kind, mol, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, "distorted methane: "//err%message)
            return
         end if
         associate (sym => symmetry_of(cav))
            call check(error, sym%orient%frame%level == frame_level_sphere, &
               & "distorted methane: frame is not at the sphere level")
            if (allocated(error)) return
         end associate
         deallocate (cav)
         do iat = 1, mol%nat
            do ic = 1, 3
               shifted = mol
               shifted%xyz(ic, iat) = mol%xyz(ic, iat) + step
               call frame_area(kind, shifted, settings, ctx, ap, error)
               if (allocated(error)) return
               shifted%xyz(ic, iat) = mol%xyz(ic, iat) - step
               call frame_area(kind, shifted, settings, ctx, am, error)
               if (allocated(error)) return
               fd(ic, iat) = (ap - am)/(2.0_wp*step)
            end do
         end do
         call check(error, maxval(abs(grad - fd)) < 5.0e-6_wp*max(1.0_wp, maxval(abs(fd))), &
            & "distorted methane ("//trim(kind_name(kind))//"): area gradient differs from finite differences")
         if (allocated(error)) return
      end do

      call new(atom, [8], reshape([0.3_wp, -0.2_wp, 0.5_wp], [3, 1]))
      do kind = kind_iswig, kind_numsa
         call frame_area(kind, atom, settings, ctx, area, error, grad=grad(:, 1:1))
         if (allocated(error)) return
         call updated_cavity(kind, atom, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, "atom: "//err%message)
            return
         end if
         associate (sym => symmetry_of(cav))
            call check(error, sym%orient%frame%level == frame_level_none, "atom: not the identity frame")
            if (allocated(error)) return
         end associate
         ! NumSA adds its probe to the radius
         call check(error, (kind == kind_numsa .or. abs(area - 4.0_wp*pi*cav%radii(1)**2) < 1.0e-10_wp*area) &
            & .and. area > 0.0_wp .and. maxval(abs(grad(:, 1))) < 1.0e-12_wp, &
            & "atom ("//trim(kind_name(kind))//"): not a full sphere at rest")
         if (allocated(error)) return
         deallocate (cav)
      end do
   end subroutine test_sphere_frame_and_atom

   !> The harmonic and line thresholds move their switches and flags
   !>
   !> - a lone weak harmonic: flagged, warned with its own text, and neither
   !>   once the threshold is far below or far above it
   !> - slightly bent CO2: flagged at the default line threshold, not with a
   !>   threshold far below or far above its distance from the line
   !>
   !> @param[out] error Error handle
   subroutine test_harmonic_linear_thresholds(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(grid_model_symmetry_type) :: model
      class(radius_type), allocatable :: radius_model
      type(mctc_error), allocatable :: err
      character(len=line_len) :: out(64)
      real(wp), parameter :: harmonic(3) = [1.0e-4_wp, 1.0e-7_wp, 2.0e-2_wp]
      real(wp), parameter :: linear(3) = [1.0e-4_wp, 1.0e-9_wp, 1.0e-1_wp]
      integer :: k, nout, ctx_unit

      call frame_structure("ring_weak", mol)
      do k = 1, 3
         settings = moist_cavity_grid_model_symmetry_parameters_type()
         settings%orientation = "molecular"
         settings%harmonic_threshold = harmonic(k)
         open (newunit=ctx_unit, status="scratch", action="readwrite")
         call new_context(ctx, nthreads=0, verbosity=0, unit=ctx_unit)
         call host_model(mol, settings, ctx, model, radius_model, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            close (ctx_unit)
            return
         end if
         call read_printout(ctx_unit, out, nout)
         close (ctx_unit)
         call ctx%delete()
         call check(error, (model%near_harmonic() .eqv. k == 1) .and. .not. (model%near_gap() .or. &
            & model%near_linear() .or. model%tied_frames()), "weak harmonic: unexpected flags")
         if (allocated(error)) return
         call check(error, model%orient%frame%harmonic, merge(0, 3, k == 3), &
            & "weak harmonic: used or skipped against its threshold")
         if (allocated(error)) return
         if (k == 1) then
            call check(error, nout >= 1, "weak harmonic: no warning")
            if (allocated(error)) return
            call check(error, any(index(out(1:nout), "weak or nearly tied in-plane harmonic") > 0), &
               & "weak harmonic: warning does not name the harmonic")
            if (allocated(error)) return
         end if
      end do

      call frame_structure("co2_near", mol)
      do k = 1, 3
         settings = moist_cavity_grid_model_symmetry_parameters_type()
         settings%orientation = "molecular"
         settings%linear_threshold = linear(k)
         call new_context(ctx, nthreads=0, verbosity=0)
         call host_model(mol, settings, ctx, model, radius_model, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         call check(error, model%near_linear() .eqv. k == 1, &
            & "slightly bent CO2: unexpected line flag")
         if (allocated(error)) return
         call ctx%delete()
      end do

      settings%harmonic_threshold = -1.0_wp
      call new_grid_model_symmetry(model, err, param=settings)
      call check(error, allocated(err), "negative harmonic threshold accepted")
      if (allocated(error)) return
      deallocate (err)
      settings%harmonic_threshold = 1.0e-4_wp
      settings%linear_threshold = -1.0_wp
      call new_grid_model_symmetry(model, err, param=settings)
      call check(error, allocated(err), "negative line threshold accepted")
   end subroutine test_harmonic_linear_thresholds

   !> The surface of a cavity maps onto itself under every operation used
   !>
   !> Every point (NumSA: quadrature node) has an image point on the image
   !> sphere with the image normal and the same area, switching, width and
   !> volume element; NumSA: equivalent spheres have equal areas
   !>
   !> - iSwiG and NumSA to `E2E_TOL`, DROP to `E2E_TOL_DROP`
   !>
   !> @param[in]  kind  Cavity kind
   !> @param[in]  name  Molecule
   !> @param[in]  group Symmetry group setting
   !> @param[out] error Test failure
   !> @param[in]  used  Group the cavity must use
   subroutine check_exact_symmetry(kind, name, group, error, used)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Molecule
      character(len=*), intent(in) :: name
      !> Symmetry group setting
      character(len=*), intent(in) :: group
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Group the cavity must use
      character(len=*), intent(in), optional :: used

      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      type(surface_points_type) :: pts
      real(wp) :: g(3, 3), shift(3), tol
      integer, allocatable :: back(:)
      character(len=:), allocatable :: message
      character(len=8) :: label
      logical :: ok
      integer :: ig, ia

      tol = E2E_TOL
      if (kind == kind_drop) tol = E2E_TOL_DROP
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_frame_fixture(name, group, fix, ctx, error, kind=kind)
      if (allocated(error)) return
      call new_frame_cavity(fix, fix%mol, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, fix%name//": "//err%message)
         return
      end if
      call surface_points_of(cav, pts)
      allocate (back(fix%mol%nat))
      associate (sym_ => symmetry_of(cav))
         call check(error, size(sym_%orient%ops, 3) > 1, fix%name//" ("//group//"): no symmetry used")
         if (allocated(error)) return
         if (present(used)) then
            call check(error, sym_%orient%group_used == used, fix%name//" ("//group//"): expected "// &
               & used//", got "//sym_%orient%group_used)
            if (allocated(error)) return
         end if
         do ig = 1, size(sym_%orient%ops, 3)
            g = sym_%orient%ops(:, :, ig)
            shift = sym_%orient%center - matmul(g, sym_%orient%center)
            ! Atom i is the image of atom back(i)
            do ia = 1, fix%mol%nat
               back(sym_%orient%perm(ia, ig)) = ia
            end do
            call match_surface_points(pts, pts, g, shift, back, tol, ok, message)
            write (label, "(i0)") ig
            call check(error, ok, fix%name//" ("//group//"): surface is not invariant under operation "// &
               & trim(label)//": "//message)
            if (allocated(error)) return
            select type (cav)
            type is (cavity_type_numsa)
               call check(error, maxval(abs(cav%asph(sym_%orient%perm(:, ig)) - cav%asph)) <= &
                  & E2E_TOL*maxval(cav%asph), fix%name//" ("//group//"): equivalent spheres differ in area")
               if (allocated(error)) return
            end select
         end do
      end associate
   end subroutine check_exact_symmetry

   !> NumSA nodes and sphere areas are exactly symmetric with a group
   !>
   !> @param[out] error Error handle
   subroutine test_numsa_exact_symmetry(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(7) = [character(len=7) :: "water", "ammonia", "methane", &
         & "bf3", "benzene", "sf6", "hocl"]
      integer :: im

      do im = 1, size(names)
         call check_exact_symmetry(kind_numsa, trim(names(im)), "auto", error)
         if (allocated(error)) return
      end do
   end subroutine test_numsa_exact_symmetry

   !> The projected DROP surface is symmetric with a group, to its convergence
   !>
   !> @param[out] error Error handle
   subroutine test_drop_exact_symmetry(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(7) = [character(len=7) :: "water", "ammonia", "methane", &
         & "bf3", "benzene", "sf6", "hocl"]
      integer :: im

      do im = 1, size(names)
         call check_exact_symmetry(kind_drop, trim(names(im)), "auto", error)
         if (allocated(error)) return
      end do
   end subroutine test_drop_exact_symmetry

   !> Rotation, shift and renumbering with explicit groups and subgroups
   !>
   !> - the full group by name, and subgroups; where the molecule has a
   !>   symmetry outside the subgroup (conjugate copies: D2h in D6h, D2d in
   !>   Td, Cs in C3v) the copy may take the image of the cavity under it,
   !>   as one whole
   !>
   !> @param[in]  kind  Cavity kind
   !> @param[out] error Test failure
   subroutine check_e2e_subgroups(kind, error)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ncase = 12
      character(len=7), parameter :: names(ncase) = [character(len=7) :: "water", "water", "water", &
         & "ammonia", "ammonia", "benzene", "benzene", "benzene", "methane", "methane", "benzene", "benzene"]
      character(len=3), parameter :: groups(ncase) = [character(len=3) :: "C2v", "C2", "Cs", "C3v", "Cs", &
         & "D6h", "D2h", "D3d", "Td", "D2d", "C2h", "C6h"]
      logical, parameter :: image(ncase) = [.false., .false., .true., .false., .true., .false., .true., &
         & .true., .false., .true., .true., .true.]
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      integer :: ic

      do ic = 1, ncase
         settings%group = trim(groups(ic))
         call check_e2e_rotation(kind, trim(names(ic)), settings, error, group=trim(groups(ic)), &
            & allow_image=image(ic))
         if (allocated(error)) return
      end do
   end subroutine check_e2e_subgroups

   !> iSwiG follows rigid motions with explicit groups and subgroups
   !>
   !> @param[out] error Error handle
   subroutine test_e2e_subgroups(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_e2e_subgroups(kind_iswig, error)
   end subroutine test_e2e_subgroups

   !> DROP follows rigid motions with explicit groups and subgroups
   !>
   !> @param[out] error Error handle
   subroutine test_drop_e2e_subgroups(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_e2e_subgroups(kind_drop, error)
   end subroutine test_drop_e2e_subgroups

   !> NumSA follows rigid motions with explicit groups and subgroups
   !>
   !> @param[out] error Error handle
   subroutine test_numsa_e2e_subgroups(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_e2e_subgroups(kind_numsa, error)
   end subroutine test_numsa_e2e_subgroups

   !> DROP gradients along rigid motions: further groups, and no group
   !>
   !> @param[out] error Error handle
   subroutine test_drop_rigid_rotation_more(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(4) = [character(len=7) :: "hocl", "benzene", "methane", "co2"]
      integer :: im

      do im = 1, size(names)
         call check_rigid_paths(kind_drop, trim(names(im)), "auto", error)
         if (allocated(error)) return
      end do
      call check_rigid_paths(kind_drop, "c1", "off", error)
      if (allocated(error)) return
      call check_rigid_paths(kind_drop, "hocl", "off", error)
      if (allocated(error)) return
      call check_rigid_paths(kind_drop, "water", "C2", error, npath=4)
   end subroutine test_drop_rigid_rotation_more

   !> Linear molecules in every cavity: poses, rigid motions and stretches
   !>
   !> Tilted CO2 and HCN without a group, with the detected group and with
   !> finite subgroups
   !>
   !> - the area and the volume do not depend on the pose or the group
   !> - gradients along rigid rotations, and along a rotation with a shift
   !>   and a displacement that keeps the group (a stretch along the line),
   !>   match central differences
   !>
   !> @param[out] error Error handle
   subroutine test_linear_cavities(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=3), parameter :: names(2) = ["co2", "hcn"]
      character(len=4), parameter :: groups(5) = [character(len=4) :: "off", "auto", "C2v", "C4", "D2h"]
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, moved
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp) :: area0, volume0
      integer :: im, kind, ig, it

      call new_context(ctx, nthreads=0, verbosity=0)
      do im = 1, 2
         call structure_of(names(im), mol)
         do kind = kind_iswig, kind_drop
            do ig = 1, size(groups)
               ! HCN has no operation that turns the line over
               if (im == 2 .and. ig == 5) cycle
               settings = moist_cavity_grid_model_symmetry_parameters_type()
               settings%orientation = "molecular"
               settings%group = trim(groups(ig))
               do it = 0, 3
                  if (it == 0) then
                     moved = mol
                  else
                     call rigid_motion(mol, test_rotation(it), 3.0_wp*sin(real(it, wp) + [0.5_wp, 2.1_wp, 3.7_wp]), &
                        & test_permutation(3, it), moved)
                  end if
                  call updated_cavity(kind, moved, settings, ctx, cav, err)
                  if (allocated(err)) then
                     call test_failed(error, names(im)//" ("//trim(groups(ig))//"): "//err%message)
                     return
                  end if
                  if (it == 0 .and. ig == 1) then
                     area0 = cav%total_area
                     volume0 = cav%total_volume
                  end if
                  call check(error, abs(cav%total_area - area0) <= 1.0e-10_wp*area0 .and. &
                     & abs(cav%total_volume - volume0) <= 1.0e-10_wp*max(volume0, 1.0_wp), &
                     & trim(kind_name(kind))//" "//names(im)//" ("//trim(groups(ig))// &
                     & "): area or volume depends on the pose or the group")
                  if (allocated(error)) return
               end do
               if (ig <= 2) then
                  call check_rigid_paths(kind, names(im), trim(groups(ig)), error)
               else
                  call check_rigid_paths(kind, names(im), trim(groups(ig)), error, npath=4)
               end if
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_linear_cavities

   !> Gradients along rigid rotations where the frame angle follows the lab
   !>
   !> The frame does not turn about its axis with the molecule: lines with a
   !> subgroup that keeps no rotation about the line (HCN Cs, CO2 C2 and
   !> C2h), CO2 bent by 0.005 bohr (C2v) and a regular 13-ring (D13h)
   !>
   !> - iSwiG through the reverse map, NumSA through the Jacobian
   !> - the 13-ring with its subgroup C2v, iSwiG: spheres off the sites
   !>   follow the frame, the others the group, and the cavity changes along
   !>   the rotation; every point kept, see [[check_rigid_paths]]
   !> - NumSA there at a step of `FD_STEP_NUMSA_BAND`: its switching
   !>   polynomial is C1, and a point that enters or leaves a smoothing band
   !>   within the stencil (seen at steps from 4e-5 rad, not up to 3.5e-5)
   !>   costs the differences an error of order h. Rounding allows 4e-6 at
   !>   that step, on gradients of 0.2 to 1.7 whose frame term is of order 1
   !>
   !> @param[out] error Error handle
   subroutine test_lab_angle_rigid_rotation(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=8), parameter :: names(6) = [character(len=8) :: "hcn", "co2", "co2_bent", "co2", "ring13", &
         & "ring13"]
      character(len=4), parameter :: groups(6) = [character(len=4) :: "Cs", "C2", "auto", "C2h", "auto", "C2v"]
      integer :: im, kind, last

      do im = 1, size(names)
         ! DROP for the triatomics that showed the error
         last = merge(kind_drop, kind_numsa, im <= 3)
         do kind = kind_iswig, last
            if (im == 6 .and. kind == kind_numsa) then
               call check_rigid_paths(kind, trim(names(im)), trim(groups(im)), error, npath=4, &
                  & step=FD_STEP_NUMSA_BAND)
            else
               call check_rigid_paths(kind, trim(names(im)), trim(groups(im)), error, npath=4, keep_all=im == 6)
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_lab_angle_rigid_rotation

   !> iSwiG with further point groups, detected: exact, covariant, differentiable
   !>
   !> C2, C2h, Ci, D2, S4, D2d, C4v, C6v, Th, Ih, T, O and a centred cube (Oh)
   !>
   !> - the surface maps onto itself under every operation used
   !> - it follows a rotation, shift and renumbering of the molecule; where a
   !>   subgroup is used (C6v, Ih) up to a symmetry of the molecule
   !> - gradients along rigid rotations (also about the axis of C2, C2h, S4
   !>   and through Ci, where the rotation keeps the group) and along a
   !>   displacement that keeps the group match central differences; not
   !>   where a subgroup is used, whose conjugate copies the molecule cannot
   !>   tell apart
   !>
   !> @param[out] error Error handle
   subroutine test_point_groups_iswig(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ncase = 13
      character(len=14), parameter :: names(ncase) = [character(len=14) :: "h2o2", "h2o2_trans", "ci", &
         & "ethene_twisted", "s4", "allene", "pyramid", "c6v", "pyrite", "icosahedron", "t12", "o24", "cube"]
      character(len=3), parameter :: used(ncase) = [character(len=3) :: "C2", "C2h", "Ci", "D2", "S4", &
         & "D2d", "C4v", "C3v", "Th", "Th", "T", "O", "Oh"]
      logical, parameter :: reduced(ncase) = [.false., .false., .false., .false., .false., .false., .false., &
         & .true., .false., .true., .false., .false., .false.]
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      integer :: im

      settings%group = "auto"
      do im = 1, ncase
         call check_exact_symmetry(kind_iswig, trim(names(im)), "auto", error, used=trim(used(im)))
         if (allocated(error)) return
         ! Up to 25 atoms: rounding of the gradients reaches 1e-12 relative
         call check_e2e_rotation(kind_iswig, trim(names(im)), settings, error, group=trim(used(im)), &
            & allow_image=reduced(im), tol_factor=10.0_wp)
         if (allocated(error)) return
         if (reduced(im)) cycle
         call check_rigid_paths(kind_iswig, trim(names(im)), "auto", error)
         if (allocated(error)) return
      end do
   end subroutine test_point_groups_iswig

   !> NumSA and DROP with further point groups: exact surfaces, rigid motions
   !>
   !> @param[out] error Error handle
   subroutine test_point_groups_numsa_drop(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ncase = 8
      character(len=14), parameter :: names(ncase) = [character(len=14) :: "h2o2", "h2o2_trans", "ci", &
         & "ethene_twisted", "s4", "allene", "pyramid", "pyrite"]
      integer :: im

      do im = 1, ncase
         call check_exact_symmetry(kind_numsa, trim(names(im)), "auto", error)
         if (allocated(error)) return
         call check_rigid_paths(kind_numsa, trim(names(im)), "auto", error)
         if (allocated(error)) return
         call check_exact_symmetry(kind_drop, trim(names(im)), "auto", error)
         if (allocated(error)) return
      end do
      do im = 1, 6
         call check_rigid_paths(kind_drop, trim(names(im)), "auto", error)
         if (allocated(error)) return
      end do
   end subroutine test_point_groups_numsa_drop

   !> iSwiG frame gradient of a molecule without a group, every coordinate
   !>
   !> Every forward array and reverse contraction of [[frame_gradients]]
   !>
   !> @param[in]  name  Molecule
   !> @param[in]  level Frame level the reference must have
   !> @param[out] error Test failure
   subroutine check_iswig_frame_cartesian(name, level, error)
      !> Molecule
      character(len=*), intent(in) :: name
      !> Frame level
      integer, intent(in) :: level
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      real(wp), allocatable :: grads(:, :, :)

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_frame_fixture(name, "off", fix, ctx, error)
      if (allocated(error)) return
      call check(error, fix%level, level, fix%name//": frame not at the expected level")
      if (allocated(error)) return
      call check_frame_gradient(fix, ctx, cartesian_directions(fix%mol%nat), FD_STEP, FD_ABS, FD_REL, &
         & error, grads)
      if (allocated(error)) return
      call check_rigid(fix, grads, error)
      if (allocated(error)) return
      ! The frame term is needed
      call check_frame_gradient(fix, ctx, cartesian_directions(fix%mol%nat), FD_STEP, FD_ABS, FD_REL, &
         & error, grads, frame_off=.true.)
   end subroutine check_iswig_frame_cartesian

   !> iSwiG frame gradient at the sphere level
   !>
   !> A molecule without symmetry at a gap threshold no pair of eigenvalues
   !> passes: the axis from the cubic moment, every coordinate and rigid
   !> motions
   !>
   !> - atoms in pairs through a centre (the `ci` structure): the axis from
   !>   the quartic moment, along displacements that keep the centre (any
   !>   other brings the cubic moment back) and rigid motions. The in-plane
   !>   harmonics tie two frames there; the harmonic odd in the height
   !>   tells them apart, so the frame does not jump between stencil points
   !>
   !> @param[out] error Error handle
   subroutine test_iswig_gradient_sphere(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      real(wp), allocatable :: dirs(:, :, :), grads(:, :, :)
      integer :: idir, ia

      call check_sphere_axis("c1_sphere", 3, error)
      if (allocated(error)) return
      call check_iswig_frame_cartesian("c1_sphere", frame_level_sphere, error)
      if (allocated(error)) return
      call check_rigid_paths(kind_iswig, "c1_sphere", "off", error)
      if (allocated(error)) return

      call check_sphere_axis("ci_sphere", 4, error)
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_frame_fixture("ci_sphere", "off", fix, ctx, error)
      if (allocated(error)) return
      allocate (dirs(3, fix%mol%nat, 6))
      do idir = 1, 6
         do ia = 1, fix%mol%nat, 2
            dirs(:, ia, idir) = sin(1.1_wp*real(ia + 7*idir, wp) + [0.4_wp, 1.6_wp, 3.1_wp])
            dirs(:, ia + 1, idir) = -dirs(:, ia, idir)
         end do
         dirs(:, :, idir) = dirs(:, :, idir)/norm2(dirs(:, :, idir))
      end do
      call check_frame_gradient(fix, ctx, dirs, FD_STEP, FD_ABS, FD_REL, error, grads)
      if (allocated(error)) return
      call check_frame_gradient(fix, ctx, dirs, FD_STEP, FD_ABS, FD_REL, error, grads, frame_off=.true.)
      if (allocated(error)) return
      call check_rigid_paths(kind_iswig, "ci_sphere", "off", error, npath=4)
   end subroutine test_iswig_gradient_sphere

   !> Fail unless a sphere-level fixture takes its axis from the given moment
   !>
   !> @param[in]  name  Molecule
   !> @param[in]  order Order of the moment, 3 or 4
   !> @param[out] error Test failure
   subroutine check_sphere_axis(name, order, error)
      !> Molecule
      character(len=*), intent(in) :: name
      !> Order of the moment
      integer, intent(in) :: order
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(frame_fixture) :: fix
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_frame_fixture(name, "off", fix, ctx, error)
      if (allocated(error)) return
      call new_frame_cavity(fix, fix%mol, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, fix%name//": "//err%message)
         return
      end if
      associate (sym_ => symmetry_of(cav))
         call check(error, sym_%orient%frame%level == frame_level_sphere .and. &
            & sym_%orient%frame%axis_order == order, fix%name//": axis not from the expected moment")
      end associate
   end subroutine check_sphere_axis

   !> NumSA frame gradient at the sphere level, area and `dsdr`
   !>
   !> @param[out] error Error handle
   subroutine test_numsa_gradient_sphere(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_numsa_frame_cartesian("c1_sphere", frame_level_sphere, error)
   end subroutine test_numsa_gradient_sphere

   !> DROP frame gradient at the sphere level, totals and per point
   !>
   !> @param[out] error Error handle
   subroutine test_drop_gradient_sphere(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_drop_frame_gradient("c1_sphere", frame_level_sphere, error)
   end subroutine test_drop_gradient_sphere

   !> One atom with a group setting: accepted without symmetry, or refused
   !>
   !> @param[out] error Error handle
   subroutine test_single_atom_group(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: atom
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      integer :: kind

      call new_context(ctx, nthreads=0, verbosity=0)
      call new(atom, [8], reshape([0.3_wp, -0.2_wp, 0.5_wp], [3, 1]))
      do kind = kind_iswig, kind_drop
         settings%group = "D6h"
         call updated_cavity(kind, atom, settings, ctx, cav, err)
         if (allocated(err)) then
            call test_failed(error, trim(kind_name(kind))//": one atom with D6h: "//err%message)
            return
         end if
         associate (sym => symmetry_of(cav))
            call check(error, sym%orient%group_used == "C1" .and. index(sym%orient%reason, "single atom") > 0, &
               & trim(kind_name(kind))//": one atom with D6h: no reason given")
            if (allocated(error)) return
         end associate
         deallocate (cav)
         settings%group = "Bogus"
         call updated_cavity(kind, atom, settings, ctx, cav, err)
         call check(error, allocated(err), trim(kind_name(kind))//": one atom with an unknown group accepted")
         if (allocated(error)) return
         call check(error, index(err%message, "unknown point group Bogus") > 0, &
            & trim(kind_name(kind))//": wrong error: "//err%message)
         if (allocated(error)) return
         deallocate (err)
         if (allocated(cav)) deallocate (cav)
      end do
   end subroutine test_single_atom_group

   !> Rotated, translated copy of a structure with its atoms in another order
   !>
   !> @param[in]  mol   Structure
   !> @param[in]  q     Rotation
   !> @param[in]  shift Translation, bohr
   !> @param[out] moved Moved copy, atoms reversed and cycled by one
   subroutine moved_structure(mol, q, shift, moved)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Rotation
      real(wp), intent(in) :: q(3, 3)
      !> Translation
      real(wp), intent(in) :: shift(3)
      !> Moved copy
      type(structure_type), intent(out) :: moved

      call rigid_motion(mol, q, shift, test_permutation(mol%nat, 1), moved)
   end subroutine moved_structure

   !> Area, volume and CPCM energy of host charges that depend on the element only
   !>
   !> @param[in]  mol      Structure
   !> @param[in]  settings Symmetry settings
   !> @param[in]  ctx      Run context
   !> @param[out] vals     Area, volume, energy
   !> @param[out] error    Test failure
   subroutine invariance_values(mol, settings, ctx, vals, error)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Symmetry settings
      type(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: settings
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Area, volume, energy
      real(wp), intent(out) :: vals(3)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: q(:)

      vals = 0.0_wp
      call new_iswig(mol, settings, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      vals(1) = cav%total_area
      vals(2) = cav%total_volume
      ! Net charge q(Z) on every nucleus: electrons of charge q - Z on the nuclei
      q = 0.4_wp*sin(1.3_wp*real(mol%num(mol%id), wp))
      call cpcm_host(cav, mol, mol%xyz, q - real(mol%num(mol%id), wp), ctx, vals(3), error)
   end subroutine invariance_values

   !> CPCM energy of nuclei moving with the atoms and fixed point charges
   !>
   !> Optionally its nuclear gradient in reverse mode (cavity contraction)
   !> and forward mode (`xyz1_rA`, `f1_rA`, `xi1_rA`)
   !>
   !> @param[in,out] cav          Updated cavity
   !> @param[in]     mol          Structure of the update; nuclei carry Z
   !> @param[in]     el_xyz       Fixed charge positions (3, nel)
   !> @param[in]     el_q         Fixed charges (nel)
   !> @param[in]     ctx          Run context
   !> @param[out]    energy       Polarization energy
   !> @param[out]    error        Test failure
   !> @param[out]    grad_reverse Reverse-mode nuclear gradient (3, nat)
   !> @param[out]    grad_forward Forward-mode nuclear gradient (3, nat)
   subroutine cpcm_host(cav, mol, el_xyz, el_q, ctx, energy, error, grad_reverse, grad_forward)
      !> Updated cavity
      class(cavity_type), intent(inout) :: cav
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Fixed charge positions
      real(wp), intent(in) :: el_xyz(:, :)
      !> Fixed charges
      real(wp), intent(in) :: el_q(:)
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Polarization energy
      real(wp), intent(out) :: energy
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Reverse-mode gradient
      real(wp), allocatable, intent(out), optional :: grad_reverse(:, :)
      !> Forward-mode gradient
      real(wp), allocatable, intent(out), optional :: grad_forward(:, :)

      class(model_continuum_component_pcm), allocatable :: pcm
      type(coupling_type), target :: coupling
      type(response_type) :: response
      type(structure_type) :: src
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: z(:)
      integer :: ipass

      energy = 0.0_wp
      call new_pcm(pcm_cpcm, ctx, pcm, error)
      if (allocated(error)) return
      call stage_energy(pcm, mol, cav, coupling, error)
      if (allocated(error)) return
      z = real(mol%num(mol%id), wp)
      call new(src, [mol%num(mol%id), spread(1, 1, size(el_q))], &
         & reshape([mol%xyz, el_xyz], [3, mol%nat + size(el_q)]))
      call answer_sources(cav, coupling, src, [z, el_q], z)
      call pcm%get_energy(component_view(coupling), cav, energy, err)
      call check_moist_error(error, err, "CPCM energy")
      if (allocated(error)) return

      do ipass = 1, 2
         if (ipass == 1 .and. .not. present(grad_reverse)) cycle
         if (ipass == 2 .and. .not. present(grad_forward)) cycle
         call pcm%prepare_gradient(cav, coupling, err)
         call check_moist_error(error, err, "CPCM gradient staging")
         if (allocated(error)) return
         call answer_sources(cav, coupling, src, [z, el_q], z)
         if (ipass == 1) then
            allocate (grad_reverse(3, mol%nat), source=0.0_wp)
            call pcm%get_gradient(component_view(coupling), cav, response, grad_reverse, err)
         else
            allocate (grad_forward(3, mol%nat), source=0.0_wp)
            call pcm%get_forward_gradient(component_view(coupling), cav, response, grad_forward, err)
         end if
         call check_moist_error(error, err, "CPCM gradient")
         if (allocated(error)) return
      end do
   end subroutine cpcm_host

   !> Reference state of a frame-gradient check
   !>
   !> Fixed host electrons up to 0.25 bohr off the reference nuclei
   !>
   !> - DROP: the settings of the per-point DROP checks ([[drop_fd_param]]:
   !>   50 Lebedev points, projection level 1, prune level 4, tolerance
   !>   1e-14) with normals and curvatures, and the step of the totals check
   !>
   !> @param[in]  name  Molecule
   !> @param[in]  group Symmetry group setting; the orientation is molecular
   !> @param[out] fix   Fixture
   !> @param[in]  ctx   Run context
   !> @param[out] error Test failure
   !> @param[in]  cut_f Switching cutoff (default 1e-7, iSwiG)
   !> @param[in]  kind  Cavity kind (default iSwiG)
   !> @param[in]  mol   Structure (default the molecule of that name)
   subroutine new_frame_fixture(name, group, fix, ctx, error, cut_f, kind, mol)
      !> Molecule
      character(len=*), intent(in) :: name
      !> Symmetry group setting
      character(len=*), intent(in) :: group
      !> Fixture
      type(frame_fixture), intent(out) :: fix
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Switching cutoff
      real(wp), intent(in), optional :: cut_f
      !> Cavity kind
      integer, intent(in), optional :: kind
      !> Structure
      type(structure_type), intent(in), optional :: mol

      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      integer :: iat

      if (present(kind)) fix%kind = kind
      fix%name = trim(kind_name(fix%kind))//" "//name
      if (present(cut_f)) fix%cut_f = cut_f
      if (fix%kind == kind_drop) then
         fix%num_leb = 50
         fix%step = FD_STEP_DROP
      end if
      if (present(mol)) then
         fix%mol = mol
      else
         call structure_of(name, fix%mol)
      end if
      fix%settings%orientation = "molecular"
      fix%settings%group = group
      ! The stencil moves the eigen-gap of the near-symmetric propeller by
      ! about 1e-4: a threshold well above it keeps every point at one level
      if (name == "c1_plane") fix%settings%gap_threshold = 1.0e-2_wp
      ! No pair of eigenvalues counts as separated: the sphere level, on a
      ! molecule whose moments are far from any tie
      if (index(name, "_sphere") > 0) fix%settings%gap_threshold = 0.99_wp
      fix%pinned = group == "off"
      allocate (fix%el_xyz(3, fix%mol%nat), fix%el_q(fix%mol%nat))
      do iat = 1, fix%mol%nat
         fix%el_xyz(:, iat) = fix%mol%xyz(:, iat) + 0.25_wp*sin(1.9_wp*real(iat, wp) + [0.1_wp, 2.2_wp, 4.3_wp])
         fix%el_q(iat) = 0.6_wp*sin(1.7_wp*real(iat, wp) + 0.3_wp) - real(fix%mol%num(fix%mol%id(iat)), wp)
      end do
      call new_frame_cavity(fix, fix%mol, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, fix%name//": "//err%message)
         return
      end if
      fix%numbering = cavity_numbering(cav)
      if (fix%kind == kind_drop) then
         fix%pinned_ids = pack(fix%numbering, cav%a(1:cav%ngrid) > drop_pin_floor*maxval(cav%a(1:cav%ngrid)))
      else
         fix%pinned_ids = fix%numbering
      end if
      associate (gm_ => symmetry_of(cav))
         fix%level = gm_%orient%frame%level
      end associate
      associate (gm_ => symmetry_of(cav))
         fix%group_used = gm_%orient%group_used
      end associate
      associate (gm_ => symmetry_of(cav))
         fix%nrot = model_nrot(gm_)
      end associate
      allocate (fix%rot(3, 3, fix%mol%nat))
      do iat = 1, fix%mol%nat
         associate (gm_ => symmetry_of(cav))
            fix%rot(:, :, iat) = model_rotation(gm_, iat)
         end associate
      end do
   end subroutine new_frame_fixture

   !> Updated cavity of a frame-gradient fixture at a structure
   !>
   !> The symmetry model is built from the fixture's settings and computed by
   !> the cavity
   !>
   !> @param[in]  fix Fixture
   !> @param[in]  mol Structure
   !> @param[in]  ctx Run context
   !> @param[out] cav Updated cavity
   !> @param[out] err Construction or update error
   subroutine new_frame_cavity(fix, mol, ctx, cav, err)
      !> Fixture
      type(frame_fixture), intent(in) :: fix
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Updated cavity
      class(cavity_type), allocatable, intent(out) :: cav
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: err

      if (fix%kind == kind_drop) then
         block
            type(cavity_type_drop), allocatable :: drop
            type(moist_cavity_drop_parameters_type) :: param
            type(grid_model_symmetry_type) :: model

            param = drop_fd_param()
            param%num_leb = fix%num_leb
            param%do_normal = .true.
            param%do_curvature = .true.
            call new_grid_model_symmetry(model, err, param=fix%settings)
            if (allocated(err)) return
            allocate (drop)
            call new_cavity_drop(drop, default_cpcm_radii(), drop_fd_svdw(), err, param=param, ctx=ctx, &
               & grid_model=model)
            if (allocated(err)) return
            call drop%update(mol, err)
            call move_alloc(drop, cav)
         end block
      else if (fix%kind == kind_numsa) then
         block
            type(cavity_type_numsa), allocatable :: numsa
            type(grid_model_symmetry_type) :: model

            call new_grid_model_symmetry(model, err, param=fix%settings)
            if (allocated(err)) return
            allocate (numsa)
            call new_cavity_numsa(numsa, default_cpcm_radii(), err, &
               & param=moist_cavity_numsa_parameters_type(num_leb=fix%num_leb), ctx=ctx, grid_model=model)
            if (allocated(err)) return
            call numsa%update(mol, err)
            call move_alloc(numsa, cav)
         end block
      else
         block
            type(cavity_type_iswig), allocatable :: iswig

            call new_iswig(mol, fix%settings, ctx, iswig, err, num_leb=fix%num_leb, cut_f=fix%cut_f)
            call move_alloc(iswig, cav)
         end block
      end if
   end subroutine new_frame_cavity

   !> Every analytic gradient of a fixture against central differences
   !>
   !> With `frame_off`, the gradients leave the frame term out (the model
   !> switched to the lab orientation after the update, the anchors kept):
   !> every analytic gradient must then miss the differences in some direction
   !>
   !> @param[in]  fix       Fixture
   !> @param[in]  ctx       Run context
   !> @param[in]  dirs      Displacement directions (3, nat, ndir)
   !> @param[in]  step      Central-difference step, bohr
   !> @param[in]  thr_abs   Absolute tolerance
   !> @param[in]  thr_rel   Relative tolerance
   !> @param[out] error     Test failure
   !> @param[out] grads     Analytic gradients (3, nat, ngrad)
   !> @param[in]  frame_off Leave the frame term out and expect failures
   subroutine check_frame_gradient(fix, ctx, dirs, step, thr_abs, thr_rel, error, grads, frame_off)
      !> Fixture
      type(frame_fixture), intent(in) :: fix
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Displacement directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Absolute tolerance
      real(wp), intent(in) :: thr_abs
      !> Relative tolerance
      real(wp), intent(in) :: thr_rel
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Analytic gradients
      real(wp), allocatable, intent(out) :: grads(:, :, :)
      !> Leave the frame term out
      logical, intent(in), optional :: frame_off

      real(wp) :: vals(nquant, 4), numeric(nquant), analytic
      character(len=200) :: context
      logical :: off
      !> Per gradient, whether some direction misses the differences
      logical :: missed(ngrad)
      integer :: idir, k, igrad

      off = .false.
      if (present(frame_off)) off = frame_off
      call frame_gradients(fix, ctx, grads, error, frame_off=off)
      if (allocated(error)) return
      do igrad = 1, ngrad
         if (grad_quantity(igrad) == q_adjoint .and. .not. fix%pinned) cycle
         call check(error, maxval(abs(grads(:, :, igrad))) > 1.0e-8_wp, &
            & fix%name//": "//trim(grad_name(igrad))//" vanishes")
         if (allocated(error)) return
      end do

      missed = .false.
      do idir = 1, size(dirs, 3)
         do k = 1, 4
            call frame_values(fix, fix%mol%xyz + fd4_offsets(k)*step*dirs(:, :, idir), ctx, vals(:, k), error)
            if (allocated(error)) return
         end do
         do k = 1, nquant
            call fd4_scalar(vals(k, 1), vals(k, 2), vals(k, 3), vals(k, 4), step, numeric(k), error)
            if (allocated(error)) return
         end do
         do igrad = 1, ngrad
            if (grad_quantity(igrad) == q_adjoint .and. .not. fix%pinned) cycle
            analytic = sum(grads(:, :, igrad)*dirs(:, :, idir))
            if (off) then
               ! The criterion of the check below
               missed(igrad) = missed(igrad) .or. abs(analytic - numeric(grad_quantity(igrad))) &
                  & > max(thr_abs, abs(thr_rel*numeric(grad_quantity(igrad))))
               cycle
            end if
            write (context, '(a,", ",a,", direction ",i0,": ",es14.6," vs FD ",es14.6)') fix%name, &
               & trim(grad_name(igrad)), idir, analytic, numeric(grad_quantity(igrad))
            call check(error, analytic, numeric(grad_quantity(igrad)), thr_abs=thr_abs, thr_rel=thr_rel, &
               & message="frame gradient disagrees with central differences", more=trim(context))
            if (allocated(error)) return
         end do
      end do
      if (.not. off) return
      do igrad = 1, ngrad
         if (grad_quantity(igrad) == q_adjoint .and. .not. fix%pinned) cycle
         call check(error, missed(igrad), fix%name//": "//trim(grad_name(igrad))// &
            & " agrees with central differences without the frame term")
         if (allocated(error)) return
      end do
   end subroutine check_frame_gradient

   !> Analytic gradients of the fixture quantities at the reference
   !>
   !> @param[in]  fix       Fixture
   !> @param[in]  ctx       Run context
   !> @param[out] grads     Gradients (3, nat, ngrad)
   !> @param[out] error     Test failure
   !> @param[in]  frame_off Switch the model to the lab orientation after the
   !>                       update: rotated spheres, no frame term
   subroutine frame_gradients(fix, ctx, grads, error, frame_off)
      !> Fixture
      type(frame_fixture), intent(in) :: fix
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Gradients
      real(wp), allocatable, intent(out) :: grads(:, :, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Leave the frame term out
      logical, intent(in), optional :: frame_off

      class(cavity_type), allocatable :: cav
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: grad_reverse(:, :), grad_forward(:, :), w(:, :)
      integer, allocatable :: numbering(:)
      real(wp) :: origin(3), w_c(3), energy
      integer :: nat, ip

      nat = fix%mol%nat
      allocate (grads(3, nat, ngrad), source=0.0_wp)
      call new_frame_cavity(fix, fix%mol, ctx, cav, err)
      if (.not. allocated(err) .and. present(frame_off)) then
         if (frame_off) then
            block
               type(grid_model_symmetry_type), pointer :: sym_ptr
               sym_ptr => symmetry_of(cav)
               sym_ptr%orient%orientation = orient_lab
            end block
         end if
      end if
      if (.not. allocated(err)) call cav%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, fix%name//": "//err%message)
         return
      end if
      select type (cav)
      type is (cavity_type_iswig)
         grads(:, :, 1) = cav%area_grad
         grads(:, :, 3) = cav%volume_grad
      type is (cavity_type_drop)
         grads(:, :, 1) = cav%A_tot1_rA
         grads(:, :, 3) = cav%V_tot1_rA
      end select
      grads(:, :, 4) = sum(cav%v1_rA, dim=3)

      ! Area: a unit area adjoint
      call acc%init(cav%ngrid)
      acc%w_a = 1.0_wp
      call cav%get_surface_gradient(acc, grads(:, :, 2), err)
      if (allocated(err)) then
         call test_failed(error, fix%name//": "//err%message)
         return
      end if

      ! Volume: the adjoints of the pressure-volume term at unit pressure
      call acc%init(cav%ngrid)
      origin = cav%get_volume_origin()
      w_c = 0.0_wp
      do ip = 1, cav%ngrid
         acc%w_a(ip) = dot_product(cav%xyz(:, ip) - origin, cav%normal0(:, ip))/3.0_wp
         acc%w_xyz(:, ip) = cav%a(ip)/3.0_wp*cav%normal0(:, ip)
         acc%w_n(:, ip) = cav%a(ip)/3.0_wp*(cav%xyz(:, ip) - origin)
         w_c = w_c - acc%w_xyz(:, ip)
      end do
      acc%w_c = w_c
      call cav%get_surface_gradient(acc, grads(:, :, 5), err)
      if (allocated(err)) then
         call test_failed(error, fix%name//": "//err%message)
         return
      end if

      ! Generic surface adjoint energy, every channel
      if (fix%pinned) then
         call acc%init(cav%ngrid)
         numbering = cavity_numbering(cav)
         allocate (w(nchannel, cav%ngrid))
         do ip = 1, cav%ngrid
            w(:, ip) = pinned_weights(fix, numbering(ip))
         end do
         acc%w_xi = w(1, :)
         acc%w_f = w(2, :)
         acc%w_a = w(3, :)
         acc%w_w = w(4, :)
         acc%w_k1 = w(5, :)
         acc%w_k2 = w(6, :)
         acc%w_xyz = w(7:9, :)
         acc%w_n = w(10:12, :)
         call cav%get_surface_gradient(acc, grads(:, :, 6), err)
         if (allocated(err)) then
            call test_failed(error, fix%name//": "//err%message)
            return
         end if
      end if

      call cpcm_host(cav, fix%mol, fix%el_xyz, fix%el_q, ctx, energy, error, grad_reverse, grad_forward)
      if (allocated(error)) return
      grads(:, :, 7) = grad_reverse
      grads(:, :, 8) = grad_forward
   end subroutine frame_gradients

   !> Fixture quantities at a displaced geometry
   !>
   !> Fails when a stencil point leaves the reference branch: another frame
   !> level, group or number of rotations, a jump of a sphere rotation, or
   !> (pinned weights) another surviving point set
   !>
   !> - curvatures: the cavity's own (DROP), 1/r on a sphere (iSwiG)
   !>
   !> @param[in]  fix   Fixture
   !> @param[in]  xyz   Displaced coordinates (3, nat)
   !> @param[in]  ctx   Run context
   !> @param[out] vals  Area, volume, surface adjoint energy, CPCM energy
   !> @param[out] error Test failure
   !> @param[out] ngrid Number of points
   subroutine frame_values(fix, xyz, ctx, vals, error, ngrid)
      !> Fixture
      type(frame_fixture), intent(in) :: fix
      !> Displaced coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Quantities
      real(wp), intent(out) :: vals(nquant)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Number of points
      integer, intent(out), optional :: ngrid

      type(structure_type) :: mol
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      integer, allocatable :: numbering(:)
      real(wp), allocatable :: wleb(:)
      real(wp) :: w(nchannel), k1, k2
      integer :: ip

      vals = 0.0_wp
      if (present(ngrid)) ngrid = 0
      mol = fix%mol
      mol%xyz = xyz
      call new_frame_cavity(fix, mol, ctx, cav, err)
      if (allocated(err)) then
         call test_failed(error, fix%name//": "//err%message)
         return
      end if
      if (present(ngrid)) ngrid = cav%ngrid
      call check_stencil_branch(fix, cav, error)
      if (allocated(error)) return
      numbering = cavity_numbering(cav)
      if (fix%pinned) then
         if (fix%kind == kind_iswig) then
            call check(error, size(numbering) == size(fix%numbering), &
               & fix%name//": a stencil point changes the surviving points")
            if (allocated(error)) return
            call check(error, all(numbering == fix%numbering), &
               & fix%name//": a stencil point changes the surviving points")
            if (allocated(error)) return
         end if
         do ip = 1, size(fix%pinned_ids)
            call check(error, any(numbering == fix%pinned_ids(ip)), &
               & fix%name//": a stencil point loses a point of the surface adjoint energy")
            if (allocated(error)) return
         end do
      end if

      vals(q_area) = cav%total_area
      vals(q_volume) = cav%total_volume
      if (fix%pinned) then
         select type (cav)
         type is (cavity_type_iswig)
            wleb = cav%wleb
         type is (cavity_type_drop)
            wleb = cav%wleb
         class default
            allocate (wleb(cav%ngrid), source=0.0_wp)
         end select
         do ip = 1, cav%ngrid
            w = pinned_weights(fix, numbering(ip))
            if (allocated(cav%k1) .and. allocated(cav%k2)) then
               k1 = cav%k1(ip)
               k2 = cav%k2(ip)
            else
               k1 = 1.0_wp/cav%radii(cav%owner(ip))
               k2 = k1
            end if
            vals(q_adjoint) = vals(q_adjoint) + w(1)*cav%xi0(ip) + w(2)*cav%f(ip) + w(3)*cav%a(ip) &
               & + w(4)*wleb(ip) + w(5)*k1 + w(6)*k2 &
               & + dot_product(w(7:9), cav%xyz(:, ip)) + dot_product(w(10:12), cav%normal0(:, ip))
         end do
      end if
      call cpcm_host(cav, mol, fix%el_xyz, fix%el_q, ctx, vals(q_cpcm), error)
   end subroutine frame_values

   !> Fail when a stencil point leaves the reference branch
   !>
   !> Another frame level, group or number of rotations; without a group
   !> (pinned), a jump of a sphere rotation
   !>
   !> @param[in]  fix   Fixture
   !> @param[in]  cav   Cavity at the stencil point
   !> @param[out] error Test failure
   subroutine check_stencil_branch(fix, cav, error)
      !> Fixture
      type(frame_fixture), intent(in) :: fix
      !> Cavity at the stencil point
      class(cavity_type), intent(in) :: cav
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer :: iat

      associate (sym_ => symmetry_of(cav))
         call check(error, sym_%orient%frame%level == fix%level .and. sym_%orient%group_used == fix%group_used &
            & .and. sym_%orient%nrot == fix%nrot, fix%name//": a stencil point changes the frame level or group")
         if (allocated(error)) return
      end associate
      if (.not. fix%pinned) return
      do iat = 1, fix%mol%nat
         associate (gm_ => symmetry_of(cav))
            call check(error, maxval(abs(model_rotation(gm_, iat) - fix%rot(:, :, iat))) < 0.05_wp, &
               & fix%name//": a stencil point jumps to another sphere rotation")
         end associate
         if (allocated(error)) return
      end do
   end subroutine check_stencil_branch

   !> Surface adjoint weights of one raw point, zero unless pinned
   !>
   !> @param[in] fix Fixture
   !> @param[in] id  Raw point index
   pure function pinned_weights(fix, id) result(w)
      !> Fixture
      type(frame_fixture), intent(in) :: fix
      !> Raw point index
      integer, intent(in) :: id
      !> Weights: xi, f, a, w, k1, k2, xyz (3), n (3)
      real(wp) :: w(nchannel)

      w = 0.0_wp
      if (any(fix%pinned_ids == id)) w = raw_weights(id, fix%mol%nat*fix%num_leb)
      if (fix%shared_curvature) w(6) = w(5)
   end function pinned_weights

   !> Raw point numbering of an iSwiG or DROP cavity
   !>
   !> @param[in] cav Updated cavity
   function cavity_numbering(cav) result(numbering)
      !> Updated cavity
      class(cavity_type), intent(in) :: cav
      !> Raw index of every point (ngrid)
      integer, allocatable :: numbering(:)

      select type (cav)
      type is (cavity_type_iswig)
         numbering = cav%numbering(1:cav%ngrid)
      type is (cavity_type_drop)
         numbering = cav%numbering(1:cav%ngrid)
      class default
         allocate (numbering(0))
      end select
   end function cavity_numbering

   !> No net force or torque on the rotation-invariant quantities
   !>
   !> Area and volume of a cavity that follows the molecule: every gradient
   !> of them sums to zero and has no moment about the origin
   !>
   !> @param[in]  fix   Fixture
   !> @param[in]  grads Analytic gradients (3, nat, ngrad)
   !> @param[out] error Test failure
   subroutine check_rigid(fix, grads, error)
      !> Fixture
      type(frame_fixture), intent(in) :: fix
      !> Analytic gradients
      real(wp), intent(in) :: grads(:, :, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer :: igrad

      do igrad = 1, ngrad
         if (grad_quantity(igrad) /= q_area .and. grad_quantity(igrad) /= q_volume) cycle
         call check_no_net_force(fix, grads(:, :, igrad), trim(grad_name(igrad)), error)
         if (allocated(error)) return
      end do
   end subroutine check_rigid

   !> No net force or torque on the gradient of a rotation-invariant quantity
   !>
   !> @param[in]  fix   Fixture
   !> @param[in]  grad  Analytic gradient (3, nat)
   !> @param[in]  name  Gradient in messages
   !> @param[out] error Test failure
   subroutine check_no_net_force(fix, grad, name, error)
      !> Fixture
      type(frame_fixture), intent(in) :: fix
      !> Analytic gradient
      real(wp), intent(in) :: grad(:, :)
      !> Gradient in messages
      character(len=*), intent(in) :: name
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: force(3), torque(3), scale
      integer :: iat

      force = sum(grad, dim=2)
      torque = 0.0_wp
      do iat = 1, fix%mol%nat
         torque = torque + cross_product(fix%mol%xyz(:, iat), grad(:, iat))
      end do
      scale = maxval(abs(grad))*max(1.0_wp, maxval(abs(fix%mol%xyz)))
      call check(error, norm2(force) <= 1.0e-11_wp*scale .and. norm2(torque) <= 1.0e-11_wp*scale, &
         & fix%name//": "//name//" has a net force or torque")
   end subroutine check_no_net_force

   !> Unit displacement of every Cartesian coordinate
   !>
   !> @param[in] nat Number of atoms
   pure function cartesian_directions(nat) result(dirs)
      !> Number of atoms
      integer, intent(in) :: nat
      !> Directions (3, nat, 3 nat)
      real(wp), allocatable :: dirs(:, :, :)

      integer :: iat, k

      allocate (dirs(3, nat, 3*nat), source=0.0_wp)
      do iat = 1, nat
         do k = 1, 3
            dirs(k, iat, 3*(iat - 1) + k) = 1.0_wp
         end do
      end do
   end function cartesian_directions

   !> Surface adjoint weights of one raw point in every channel
   !>
   !> Pinned to the raw index, so a rebuilt grid weights the same point alike
   !>
   !> @param[in] iraw Raw point index
   !> @param[in] nraw Number of raw points
   pure function raw_weights(iraw, nraw) result(w)
      !> Raw point index
      integer, intent(in) :: iraw
      !> Number of raw points
      integer, intent(in) :: nraw
      !> Weights: xi, f, a, w, k1, k2, xyz (3), n (3)
      real(wp) :: w(nchannel)

      integer :: ich

      do ich = 1, nchannel
         w(ich) = sin(0.37_wp*real(iraw, wp) + 1.7_wp*real(ich, wp))*10.0_wp/real(nraw, wp)
      end do
   end function raw_weights
   !* ================================================================================= *!
   !*                           End-to-end rotation (iSwiG)                             *!
   !* ================================================================================= *!
   !
   ! Cavities built exactly like a host builds them: a radius model and a
   ! symmetry model updated for the structure, passed into the public
   ! constructor, then update, cavity gradient and a CPCM energy with its
   ! nuclear gradient. The same for rigidly rotated, translated and permuted
   ! copies; everything is compared after mapping the copy back (fixtures in
   ! `test_cavity_symmetry_fixtures`, meant for SvdW-DROP as well)

   !> Molecular frame without a group: covariant end to end
   !>
   !> - C1 molecule: the frame has no ties
   !> - ammonia and benzene with the group off: a tie of the frame's discrete
   !>   choices may give the image of the cavity under a symmetry of the
   !>   molecule, which only permutes equivalent grids
   !>
   !> @param[out] error Error handle
   subroutine test_e2e_frame(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(3) = [character(len=7) :: "c1", "ammonia", "benzene"]
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      integer :: im

      settings%orientation = "molecular"
      settings%group = "off"
      do im = 1, size(names)
         call check_e2e_rotation(kind_iswig, trim(names(im)), settings, error)
         if (allocated(error)) return
      end do
   end subroutine test_e2e_frame

   !> Symmetric frame (group auto): covariant end to end
   !>
   !> @param[out] error Error handle
   subroutine test_e2e_symmetric(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(4) = [character(len=7) :: "water", "ammonia", "methane", &
         & "benzene"]
      character(len=*), parameter :: groups(4) = [character(len=3) :: "C2v", "C3v", "Td", "D6h"]
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      integer :: im

      settings%group = "auto"
      do im = 1, size(names)
         call check_e2e_rotation(kind_iswig, trim(names(im)), settings, error, trim(groups(im)))
         if (allocated(error)) return
      end do
   end subroutine test_e2e_symmetric

   !> Lab orientation: the same rigid motion changes the cavity measurably
   !>
   !> Proves the end-to-end checks can fail: the point sets do not match,
   !> modulo any symmetry, and area and energy move beyond 1e-9 relative
   !>
   !> @param[out] error Error handle
   subroutine test_e2e_lab_sensitive(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_e2e_lab_sensitive(kind_iswig, 1.0e-9_wp, error)
   end subroutine test_e2e_lab_sensitive

   !> The lab orientation fails the end-to-end comparison of one cavity kind
   !>
   !> - C1 molecule and benzene, one rigid motion: the point sets do not match
   !>   to the kind's tolerance, modulo any symmetry of the molecule, and area
   !>   and energy (CPCM; NumSA: weighted sphere areas) move beyond `thr`
   !>   relative
   !>
   !> @param[in]  kind  Cavity kind
   !> @param[in]  thr   Smallest relative change of area and energy
   !> @param[out] error Test failure
   subroutine check_e2e_lab_sensitive(kind, thr, error)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Smallest relative change
      real(wp), intent(in) :: thr
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, moved
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      type(e2e_state) :: ref, val
      real(wp), allocatable :: ops(:, :, :)
      integer, allocatable :: atom_perm(:, :), perm(:)
      real(wp) :: q(3, 3), center(3), tol
      character(len=:), allocatable :: message
      character(len=200) :: label
      logical :: ok
      integer :: im
      character(len=*), parameter :: names(2) = [character(len=7) :: "c1", "benzene"]

      tol = E2E_TOL
      if (kind == kind_drop) tol = E2E_TOL_DROP
      call new_context(ctx, nthreads=0, verbosity=0)
      settings%orientation = "lab"
      do im = 1, size(names)
         call structure_of(trim(names(im)), mol)
         call molecule_symmetry(mol, ops, atom_perm, center)
         call e2e_build(kind, mol, settings, ctx, ref, error)
         if (allocated(error)) return
         q = test_rotation(1)
         perm = test_permutation(mol%nat, 1)
         call rigid_motion(mol, q, [1.7_wp, -0.6_wp, 2.9_wp], perm, moved)
         call e2e_build(kind, moved, settings, ctx, val, error)
         if (allocated(error)) return
         write (label, '(a,": ",a,", area ",es10.3,", energy ",es10.3," relative")') trim(kind_name(kind)), &
            & trim(names(im)), abs(val%area - ref%area)/ref%area, abs(val%energy - ref%energy)/abs(ref%energy)
         call match_surface_points(ref%pts, val%pts, q, [1.7_wp, -0.6_wp, 2.9_wp], perm, tol, ok, message, &
            & ops, atom_perm, center)
         call check(error,.not. ok, trim(label)//": lab point sets match after a rotation")
         if (allocated(error)) return
         call check(error, abs(val%area - ref%area) > thr*ref%area .and. &
            & abs(val%energy - ref%energy) > thr*abs(ref%energy), &
            & trim(label)//": lab area or energy unchanged by a rotation")
         if (allocated(error)) return
      end do
   end subroutine check_e2e_lab_sensitive

   !> SvdW-DROP, molecular frame without a group: covariant end to end
   !>
   !> Projection level 1; C1 molecule, ammonia and benzene (frame ties give a
   !> symmetry image, as for iSwiG)
   !>
   !> @param[out] error Error handle
   subroutine test_drop_e2e_frame(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(3) = [character(len=7) :: "c1", "ammonia", "benzene"]
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      integer :: im

      settings%orientation = "molecular"
      settings%group = "off"
      do im = 1, size(names)
         call check_e2e_rotation(kind_drop, trim(names(im)), settings, error)
         if (allocated(error)) return
      end do
   end subroutine test_drop_e2e_frame

   !> SvdW-DROP, symmetric frame (group auto): covariant end to end
   !>
   !> @param[out] error Error handle
   subroutine test_drop_e2e_symmetric(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(4) = [character(len=7) :: "water", "ammonia", "methane", &
         & "benzene"]
      character(len=*), parameter :: groups(4) = [character(len=3) :: "C2v", "C3v", "Td", "D6h"]
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      integer :: im

      settings%group = "auto"
      do im = 1, size(names)
         call check_e2e_rotation(kind_drop, trim(names(im)), settings, error, trim(groups(im)))
         if (allocated(error)) return
      end do
   end subroutine test_drop_e2e_symmetric

   !> SvdW-DROP, lab orientation: the end-to-end comparison fails
   !>
   !> @param[out] error Error handle
   subroutine test_drop_e2e_lab_sensitive(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_e2e_lab_sensitive(kind_drop, 1.0e-9_wp, error)
   end subroutine test_drop_e2e_lab_sensitive

   !> NumSA, molecular frame without a group: covariant end to end
   !>
   !> C1 molecule, ammonia and benzene (frame ties give a symmetry image, as
   !> for iSwiG); quadrature nodes, total and sphere areas and their gradients
   !>
   !> @param[out] error Error handle
   subroutine test_numsa_e2e_frame(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(3) = [character(len=7) :: "c1", "ammonia", "benzene"]
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      integer :: im

      settings%orientation = "molecular"
      settings%group = "off"
      do im = 1, size(names)
         call check_e2e_rotation(kind_numsa, trim(names(im)), settings, error)
         if (allocated(error)) return
      end do
   end subroutine test_numsa_e2e_frame

   !> NumSA, symmetric frame (group auto): covariant end to end
   !>
   !> @param[out] error Error handle
   subroutine test_numsa_e2e_symmetric(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: names(4) = [character(len=7) :: "water", "ammonia", "methane", &
         & "benzene"]
      character(len=*), parameter :: groups(4) = [character(len=3) :: "C2v", "C3v", "Td", "D6h"]
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      integer :: im

      settings%group = "auto"
      do im = 1, size(names)
         call check_e2e_rotation(kind_numsa, trim(names(im)), settings, error, trim(groups(im)))
         if (allocated(error)) return
      end do
   end subroutine test_numsa_e2e_symmetric

   !> NumSA, lab orientation: the end-to-end comparison fails
   !>
   !> @param[out] error Error handle
   subroutine test_numsa_e2e_lab_sensitive(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_e2e_lab_sensitive(kind_numsa, 1.0e-9_wp, error)
   end subroutine test_numsa_e2e_lab_sensitive

   !> Compare a molecule's cavity with those of rigidly moved, permuted copies
   !>
   !> Rotations `test_rotation(1..e2e_nmotion)` (half turns included),
   !> translations up to 15 bohr and permutations `test_permutation`; per copy
   !>
   !> - number of points, total area and volume, energy (relative)
   !> - point sets mapped back, modulo a symmetry of the molecule: positions,
   !>   normals, areas, switching, widths and volume elements (NumSA: its
   !>   quadrature nodes)
   !> - total-area gradient, reverse and forward energy gradients:
   !>   grad'(:, i) = Q grad(:, perm(i))
   !> - NumSA: sphere areas and their gradients, both atom indices permuted
   !> - tolerance `E2E_TOL` for iSwiG and NumSA, `E2E_TOL_DROP` for SvdW-DROP
   !>
   !> @param[in]  kind     Cavity kind
   !> @param[in]  name     Molecule
   !> @param[in]  settings Symmetry settings
   !> @param[out] error    Test failure
   !> @param[in]  group    Group the reference must use
   !> @param[in]  allow_image Whether the copy may be a symmetry image
   !> @param[in]  tol_factor  Factor on the tolerance (default one)
   subroutine check_e2e_rotation(kind, name, settings, error, group, allow_image, tol_factor)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Molecule
      character(len=*), intent(in) :: name
      !> Symmetry settings
      type(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: settings
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Group the reference must use
      character(len=*), intent(in), optional :: group
      !> Whether the copy may take the image of the cavity under a symmetry
      !> of the molecule (an explicit subgroup with conjugate copies)
      logical, intent(in), optional :: allow_image
      !> Factor on the tolerance
      real(wp), intent(in), optional :: tol_factor

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, moved
      type(e2e_state) :: ref, val
      logical :: image_ok
      real(wp), allocatable :: ops(:, :, :)
      integer, allocatable :: atom_perm(:, :), perm(:), gperm(:)
      real(wp) :: q(3, 3), shift(3), center(3), dev, g(3, 3), tol
      character(len=:), allocatable :: message
      character(len=200) :: label
      logical :: ok
      integer :: imotion, iop, ia
      !> Atom of the original matched by every atom of the copy
      integer, allocatable :: back(:)
      !> Operations of the molecule that carry the group used onto the copy's
      logical, allocatable :: allowed(:)

      tol = E2E_TOL
      if (kind == kind_drop) tol = E2E_TOL_DROP
      if (present(tol_factor)) tol = tol_factor*tol
      call new_context(ctx, nthreads=0, verbosity=0)
      call structure_of(name, mol)
      call molecule_symmetry(mol, ops, atom_perm, center)
      call e2e_build(kind, mol, settings, ctx, ref, error)
      if (allocated(error)) return
      if (present(group)) then
         call check(error, ref%group_used == group, name//": expected "//group//", got "//ref%group_used)
         if (allocated(error)) return
      end if

      do imotion = 1, e2e_nmotion
         q = test_rotation(imotion)
         perm = test_permutation(mol%nat, imotion)
         shift = 15.0_wp*sin(1.3_wp*real(imotion, wp) + [0.5_wp, 2.1_wp, 3.7_wp])
         call rigid_motion(mol, q, shift, perm, moved)
         call e2e_build(kind, moved, settings, ctx, val, error)
         if (allocated(error)) return
         write (label, '(a,": ",a," (",a,", motion ",i0,")")') trim(kind_name(kind)), name, &
            & trim(settings%group), imotion

         call check(error, val%ngrid, ref%ngrid, trim(label)//": number of points")
         if (allocated(error)) return
         call check(error, val%group_used == ref%group_used, trim(label)//": group used changed")
         if (allocated(error)) return
         call check(error, val%area, ref%area, thr=tol, rel=.true., message=trim(label)//": total area")
         if (allocated(error)) return
         ! Scaled threshold, not rel: the volume of NumSA is exactly zero
         call check(error, val%volume, ref%volume, thr=tol*abs(ref%volume), message=trim(label)//": total volume")
         if (allocated(error)) return
         call check(error, val%energy, ref%energy, thr=tol, rel=.true., message=trim(label)//": energy")
         if (allocated(error)) return

         ! Only an operation that carries the group used by the original
         ! onto the one used by the copy: a point set of higher symmetry
         ! than the group matches under others too, whose gradients differ
         allowed = [(carries_group(iop), iop=0, size(ops, 3))]
         call match_surface_points(ref%pts, val%pts, q, shift, perm, tol, ok, message, ops, atom_perm, &
            & center, iop, allowed)
         call check(error, ok, trim(label)//": point sets differ: "//message)
         if (allocated(error)) return
         ! A group used that is the full group never needs another operation
         image_ok = settings%group == "off"
         if (present(allow_image)) image_ok = image_ok .or. allow_image
         call check(error, iop == 0 .or. image_ok, trim(label)//": cavity is a symmetry image")
         if (allocated(error)) return

         ! The gradient of the branch the copy took: that of g(original)
         if (iop == 0) then
            g = identity3
            gperm = [(ia, ia=1, mol%nat)]
         else
            g = ops(:, :, iop)
            gperm = atom_perm(:, iop)
         end if
         dev = gradient_covariance_error(ref%grad_area, val%grad_area, q, perm, g, gperm)
         call check_deviation("area gradient")
         if (allocated(error)) return
         if (allocated(ref%grad_reverse)) then
            dev = gradient_covariance_error(ref%grad_reverse, val%grad_reverse, q, perm, g, gperm)
            call check_deviation("reverse energy gradient")
            if (allocated(error)) return
         end if
         dev = gradient_covariance_error(ref%grad_forward, val%grad_forward, q, perm, g, gperm)
         call check_deviation("forward energy gradient")
         if (allocated(error)) return

         ! NumSA: per sphere
         if (allocated(ref%asph)) then
            back = [(findloc(gperm, perm(ia), dim=1), ia=1, mol%nat)]
            dev = maxval(abs(val%asph - ref%asph(back)))/maxval(abs(ref%asph))
            call check_deviation("sphere areas")
            if (allocated(error)) return
            dev = sphere_gradient_covariance_error(ref%sphere_grad, val%sphere_grad, q, perm, g, gperm)
            call check_deviation("sphere area gradients")
            if (allocated(error)) return
         end if
      end do

   contains

      !> Whether g G g^T is the group of the copy turned back, q^T G' q
      !>
      !> @param[in] jop Operation of the molecule, 0 the identity
      function carries_group(jop) result(carries)
         !> Operation of the molecule
         integer, intent(in) :: jop
         !> Whether the groups agree
         logical :: carries

         real(wp) :: gg(3, 3), a(3, 3), b(3, 3)
         integer :: i, j
         logical :: found

         carries = .true.
         if (.not. (allocated(ref%ops) .and. allocated(val%ops))) return
         gg = identity3
         if (jop > 0) gg = ops(:, :, jop)
         carries = size(ref%ops, 3) == size(val%ops, 3)
         if (.not. carries) return
         do i = 1, size(ref%ops, 3)
            a = matmul(gg, matmul(ref%ops(:, :, i), transpose(gg)))
            found = .false.
            do j = 1, size(val%ops, 3)
               b = matmul(transpose(q), matmul(val%ops(:, :, j), q))
               found = maxval(abs(a - b)) < 1.0e-6_wp
               if (found) exit
            end do
            carries = found
            if (.not. carries) return
         end do
      end function carries_group

      !> Fail unless `dev` is within the tolerance
      !>
      !> @param[in] name Compared quantity
      subroutine check_deviation(name)
         !> Compared quantity
         character(len=*), intent(in) :: name

         character(len=240) :: text

         write (text, '(a,": ",a," not covariant, relative deviation ",es10.3)') trim(label), name, dev
         call check(error, dev <= tol, trim(text))
      end subroutine check_deviation

   end subroutine check_e2e_rotation

   !> iSwiG of a linear molecule in the molecular orientation, wherever it lies
   !>
   !> - CO2 along z, shifted, and tilted and shifted: same point count and
   !>   area, area gradient equal to central differences in every coordinate;
   !>   the in-plane angle about the line comes from the lab frame, which
   !>   the axially symmetric spheres do not see
   !> - group off and auto
   !>
   !> @param[out] error Error handle
   subroutine test_iswig_linear(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, moved
      type(moist_cavity_grid_model_symmetry_parameters_type) :: settings
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      type(structure_type) :: shifted
      real(wp) :: xyz(3, 3), area0, ap, am
      integer :: ngrid0, it, ig, iat, ic

      xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.16_wp, 0.0_wp, 0.0_wp, -1.16_wp], &
         & [3, 3])*aatoau
      call new(mol, [6, 8, 8], xyz)
      call new_context(ctx, nthreads=0, verbosity=0)
      settings%orientation = "molecular"
      do ig = 1, 2
         if (ig == 2) settings%group = "auto"
         do it = 1, 3
            select case (it)
            case (1)
               moved = mol
            case (2)
               call rigid_motion(mol, identity3, [0.1_wp, 0.2_wp, -0.3_wp], [1, 2, 3], moved)
            case default
               call rigid_motion(mol, test_rotation(1), [0.4_wp, -1.1_wp, 0.7_wp], [1, 2, 3], moved)
            end select
            call updated_cavity(kind_iswig, moved, settings, ctx, cav, err)
            if (allocated(err)) then
               call test_failed(error, err%message)
               return
            end if
            if (it == 1) then
               ngrid0 = cav%ngrid
               area0 = cav%total_area
            end if
            call check(error, cav%ngrid, ngrid0, "linear molecule: point count depends on the pose")
            if (allocated(error)) return
            call check(error, abs(cav%total_area - area0) <= 1.0e-11_wp*area0, &
               & "linear molecule: area depends on the pose")
            if (allocated(error)) return
            ! A harmonic of roundoff would scale the frame term by 1/roundoff
            call cav%get_gradient(err)
            if (allocated(err)) then
               call test_failed(error, err%message)
               return
            end if
            select type (cav)
            type is (cavity_type_iswig)
               call check(error, all(abs(sum(cav%area_grad, dim=2)) <= 1.0e-9_wp), &
                  & "linear molecule: area gradient not translation invariant")
               if (allocated(error)) return
               ! Every coordinate against central differences: a harmonic of
               ! roundoff would scale the frame term by 1/roundoff
               do iat = 1, 3
                  do ic = 1, 3
                     shifted = moved
                     shifted%xyz(ic, iat) = moved%xyz(ic, iat) + 1.0e-5_wp
                     call frame_area(kind_iswig, shifted, settings, ctx, ap, error)
                     if (allocated(error)) return
                     shifted%xyz(ic, iat) = moved%xyz(ic, iat) - 1.0e-5_wp
                     call frame_area(kind_iswig, shifted, settings, ctx, am, error)
                     if (allocated(error)) return
                     call check(error, abs(cav%area_grad(ic, iat) - (ap - am)/2.0e-5_wp) <= 5.0e-6_wp* &
                        & max(1.0_wp, maxval(abs(cav%area_grad))), &
                        & "linear molecule: area gradient differs from central differences")
                     if (allocated(error)) return
                  end do
               end do
            end select
         end do
      end do
   end subroutine test_iswig_linear

   !> Build a cavity the way a host does and evaluate it
   !>
   !> Radius model and symmetry model updated for the structure (the symmetry
   !> report's state), both passed into the public constructor, update, cavity
   !> gradient, then the CPCM energy of a net charge q(Z) on every nucleus with
   !> its reverse and forward nuclear gradient
   !>
   !> - DROP: SvdW level set at projection level `e2e_drop_level`
   !> - NumSA: the energy is sum_i s(Z_i) A_i over the sphere areas, its
   !>   gradient through `dsdr`
   !>
   !> @param[in]     kind     Cavity kind
   !> @param[in]     mol      Structure
   !> @param[in]     settings Symmetry settings
   !> @param[in,out] ctx      Run context
   !> @param[out]    state    Results
   !> @param[out]    error    Test failure
   subroutine e2e_build(kind, mol, settings, ctx, state, error)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Symmetry settings
      type(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: settings
      !> Run context
      type(moist_context_type), intent(inout), target :: ctx
      !> Results
      type(e2e_state), intent(out) :: state
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(grid_model_symmetry_type) :: model
      class(radius_type), allocatable :: radius_model
      class(cavity_type), allocatable :: cav
      type(mctc_error), allocatable :: err
      type(moist_cavity_drop_parameters_type) :: drop_param
      real(wp), allocatable :: q(:)
      integer :: iat

      drop_param%proj_level = e2e_drop_level
      call host_model(mol, settings, ctx, model, radius_model, err)
      if (.not. allocated(err)) call new_cavity_of(kind, radius_model, ctx, cav, err, model, drop_param)
      if (.not. allocated(err)) call cav%update(mol, err)
      if (.not. allocated(err)) call cav%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      state%ngrid = cav%ngrid
      state%area = cav%total_area
      state%volume = cav%total_volume
      associate (gm_ => symmetry_of(cav))
         state%group_used = gm_%orient%group_used
         if (allocated(gm_%orient%ops)) state%ops = gm_%orient%ops
      end associate
      call surface_points_of(cav, state%pts)
      q = 0.4_wp*sin(1.3_wp*real(mol%num(mol%id), wp))
      select type (cav)
      type is (cavity_type_iswig)
         state%grad_area = cav%area_grad
      type is (cavity_type_drop)
         state%grad_area = cav%A_tot1_rA
      type is (cavity_type_numsa)
         state%grad_area = cav%area_grad
         state%asph = cav%asph
         state%sphere_grad = cav%dsdr
         ! Element weights s(Z) = 1 + q(Z) keep the energy permutation invariant
         state%energy = sum((1.0_wp + q)*cav%asph)
         allocate (state%grad_forward(3, mol%nat), source=0.0_wp)
         do iat = 1, mol%nat
            state%grad_forward = state%grad_forward + (1.0_wp + q(iat))*cav%dsdr(:, :, iat)
         end do
         return
      end select
      call cpcm_host(cav, mol, mol%xyz, q - real(mol%num(mol%id), wp), ctx, state%energy, error, &
         & state%grad_reverse, state%grad_forward)
   end subroutine e2e_build

   !* ================================================================================= *!
   !*                                      Helpers                                      *!
   !* ================================================================================= *!

   !> Construct and update a cavity of one kind with a symmetry model
   !>
   !> The model is built from the settings and computed by the cavity itself
   !>
   !> @param[in]  kind         Cavity kind
   !> @param[in]  mol          Structure
   !> @param[in]  settings     Symmetry settings
   !> @param[in]  ctx          Run context
   !> @param[out] cav          Cavity; allocated unless the model is refused
   !> @param[out] err          Model, construction or update error
   !> @param[in]  radius_model Radii; CPCM by default
   !> @param[in]  with_model   Pass the model (default); false passes none
   subroutine updated_cavity(kind, mol, settings, ctx, cav, err, radius_model, with_model)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Symmetry settings
      type(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: settings
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Cavity
      class(cavity_type), allocatable, intent(out) :: cav
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: err
      !> Radii
      class(radius_type), intent(in), optional :: radius_model
      !> Whether the model is passed
      logical, intent(in), optional :: with_model

      class(radius_type), allocatable :: radii
      type(grid_model_symmetry_type) :: model
      logical :: pass

      if (present(radius_model)) then
         allocate (radii, source=radius_model)
      else
         allocate (radii, source=default_cpcm_radii())
      end if
      pass = .true.
      if (present(with_model)) pass = with_model

      if (pass) then
         call new_grid_model_symmetry(model, err, param=settings)
         if (allocated(err)) return
         call new_cavity_of(kind, radii, ctx, cav, err, model)
      else
         call new_cavity_of(kind, radii, ctx, cav, err)
      end if
      if (allocated(err)) return
      call cav%update(mol, err)
   end subroutine updated_cavity

   !> Construct a cavity of one kind through its public constructor
   !>
   !> @param[in]  kind         Cavity kind
   !> @param[in]  radius_model Radii
   !> @param[in]  ctx          Run context
   !> @param[out] cav          Constructed cavity
   !> @param[out] err          Construction error
   !> @param[in]  model        Symmetry model; absent, none
   !> @param[in]  drop_param   DROP settings; absent, the compiled defaults
   subroutine new_cavity_of(kind, radius_model, ctx, cav, err, model, drop_param)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Radii
      class(radius_type), intent(in) :: radius_model
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Cavity
      class(cavity_type), allocatable, intent(out) :: cav
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: err
      !> Symmetry model
      type(grid_model_symmetry_type), intent(in), optional :: model
      !> DROP settings
      type(moist_cavity_drop_parameters_type), intent(in), optional :: drop_param

      select case (kind)
      case (kind_iswig)
         block
            type(cavity_type_iswig), allocatable :: iswig

            allocate (iswig)
            call new_cavity_iswig(iswig, radius_model, err, ctx=ctx, grid_model=model)
            call move_alloc(iswig, cav)
         end block
      case (kind_numsa)
         block
            type(cavity_type_numsa), allocatable :: numsa

            allocate (numsa)
            call new_cavity_numsa(numsa, radius_model, err, ctx=ctx, grid_model=model)
            call move_alloc(numsa, cav)
         end block
      case default
         block
            type(cavity_type_drop), allocatable :: drop
            type(moist_cavity_drop_lsf_svdw_type) :: lsf

            call lsf%new()
            allocate (drop)
            call new_cavity_drop(drop, radius_model, lsf, err, param=drop_param, ctx=ctx, grid_model=model)
            call move_alloc(drop, cav)
         end block
      end select
   end subroutine new_cavity_of

   !> The host's part: radii and the symmetry model for a structure
   !>
   !> CPCM radii updated for `mol`, then a model built from the settings and
   !> updated with those radii on the context, as the CLI driver does
   !>
   !> @param[in]     mol          Structure
   !> @param[in]     settings     Symmetry settings
   !> @param[in,out] ctx          Run context: timer and warnings
   !> @param[out]    model        Updated symmetry model
   !> @param[out]    radius_model Updated radii
   !> @param[out]    err          Radii, model or update error
   subroutine host_model(mol, settings, ctx, model, radius_model, err)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Symmetry settings
      type(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: settings
      !> Run context
      type(moist_context_type), intent(inout), target :: ctx
      !> Symmetry model
      type(grid_model_symmetry_type), intent(out) :: model
      !> Radii
      class(radius_type), allocatable, intent(out) :: radius_model
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: err

      allocate (radius_model, source=default_cpcm_radii())
      call radius_model%update(mol, err)
      if (allocated(err)) return
      call new_grid_model_symmetry(model, err, param=settings)
      if (allocated(err)) return
      call model%update(mol, radius_model%f0, err, ctx)
   end subroutine host_model

   !> Fresh parameter set of one cavity kind
   !>
   !> @param[in]  kind  Cavity kind
   !> @param[out] param Parameters with compiled defaults
   subroutine new_param(kind, param)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Parameters
      class(moist_model_parameters_type), allocatable, intent(out) :: param

      select case (kind)
      case (kind_iswig)
         allocate (moist_cavity_iswig_parameters_type :: param)
      case (kind_numsa)
         allocate (moist_cavity_numsa_parameters_type :: param)
      case default
         allocate (moist_cavity_drop_parameters_type :: param)
      end select
   end subroutine new_param

   !> Total calls of every "Detection" timer of a context
   !>
   !> @param[in,out] ctx Run context
   function detections(ctx) result(ncall)
      !> Run context
      type(moist_context_type), intent(inout) :: ctx
      !> Number of detections
      integer :: ncall

      ncall = timer_calls(ctx, "Detection")
   end function detections

   !> Total calls of every timer of a name in a context
   !>
   !> Read from the timer report, whose rows end in the call count ("12x")
   !>
   !> @param[in,out] ctx  Run context
   !> @param[in]     name Node name
   function timer_calls(ctx, name) result(ncall)
      !> Run context
      type(moist_context_type), intent(inout) :: ctx
      !> Node name
      character(len=*), intent(in) :: name
      !> Number of calls
      integer :: ncall

      character(len=line_len) :: lines(200)
      integer :: unit, nline, i, ix, n, stat

      ncall = 0
      open (newunit=unit, status="scratch", action="readwrite")
      call ctx%timer%write(unit)
      call read_printout(unit, lines, nline)
      close (unit)
      do i = 1, nline
         if (index(lines(i), " "//name//" ") == 0) cycle
         ix = index(lines(i), "x", back=.true.)
         if (ix < 10) cycle
         read (lines(i) (ix - 9:ix - 1), *, iostat=stat) n
         if (stat == 0) ncall = ncall + n
      end do
   end function timer_calls

   !> Whether a timer node exists at a path from the root
   !>
   !> - every segment must name a node at its own depth; the measured time is
   !>   not used, a fast run may round it to zero
   !>
   !> @param[in] ctx  Run context
   !> @param[in] path Slash-separated node names from the top level
   function timer_path(ctx, path) result(found)
      !> Run context
      type(moist_context_type), intent(in) :: ctx
      !> Path
      character(len=*), intent(in) :: path
      !> Whether every segment exists at its depth
      logical :: found

      integer :: id, first, last, depth, top

      top = huge(1)
      do id = 1, ctx%timer%num_nodes()
         top = min(top, ctx%timer%node_depth(id))
      end do
      first = 1
      depth = top
      found = len(path) > 0
      do while (found .and. first <= len(path))
         last = index(path(first:), "/")
         last = merge(len(path), first + last - 2, last == 0)
         found = .false.
         do id = 1, ctx%timer%num_nodes()
            if (ctx%timer%node_depth(id) /= depth) cycle
            if (ctx%timer%node_name(id) == path(first:last)) found = .true.
         end do
         first = last + 2
         depth = depth + 1
      end do
   end function timer_path

   !> Whether a context's timer has a node of a name
   !>
   !> @param[in] ctx  Run context
   !> @param[in] name Node name
   function timer_has(ctx, name) result(found)
      !> Run context
      type(moist_context_type), intent(in) :: ctx
      !> Node name
      character(len=*), intent(in) :: name
      !> Whether the node exists
      logical :: found

      integer :: id

      found = .false.
      do id = 1, ctx%timer%num_nodes()
         if (ctx%timer%node_name(id) == name) found = .true.
      end do
   end function timer_has

   !> Write a parameter file, read it and remove it again
   !>
   !> @param[in,out] param Settings that read the file
   !> @param[in]     path  Fixture path, owned by the calling test
   !> @param[in]     lines File lines, trailing blanks dropped
   !> @param[out]    err   Library error of the read
   !> @param[out]    ok    Whether the file was written and removed
   subroutine read_fixture(param, path, lines, err, ok)
      !> Settings that read the file
      class(moist_model_parameters_type), intent(inout) :: param
      !> Fixture path
      character(len=*), intent(in) :: path
      !> File lines
      character(len=*), intent(in) :: lines(:)
      !> Library error of the read
      type(mctc_error), allocatable, intent(out) :: err
      !> Whether the file was written and removed
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
      ok = delete_file(path)
   end subroutine read_fixture

   !> Remove an existing file
   !>
   !> @param[in] path File path
   function delete_file(path) result(deleted)
      !> File path
      character(len=*), intent(in) :: path
      !> Whether the file existed and was removed
      logical :: deleted

      integer :: unit, stat

      deleted = .false.
      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat /= 0) return
      close (unit, status="delete", iostat=stat)
      deleted = stat == 0
   end function delete_file

end module test_cavity_symmetry
