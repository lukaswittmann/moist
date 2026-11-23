!> Tests for cavity symmetry handling, shared by every Lebedev-based cavity
!>
!> - detection: libmsym groups of reference molecules, linear groups
!> - group used: largest subgroup whose site groups admit a Lebedev-invariant
!>   sphere (BF3 D3h -> C3v, icosahedron Ih -> Th, CO2 D_inf_h -> D4h)
!> - exactness: every site sphere invariant under its site group, every
!>   transported sphere the image of its representative, modulo Oh
!> - storage: orientations equal modulo Oh with equal derivative maps share
!>   one entry, minimal and covariant (SF6 and a 48-atom Oh cluster)
!> - covariance: rotating, translating and permuting the atoms rotates the
!>   spheres with them, up to one symmetry of the molecule, also where the
!>   discrete choices tie (neighbour forms break those ties)
!> - errors: absent or inexact requested groups
!> - derivatives: reverse-mode frame and orientation torques against central
!>   differences of an Oh-invariant function of the sphere orientations; the
!>   forward-mode orientation Jacobian against the reverse mode
module test_math_symmetry
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, i8
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io_constants, only: pi
   use mctc_io_convert, only: aatoau
   use moist_math_symmetry_msym, only: msym_group_type, find_point_group
   use moist_math_symmetry_frame, only: sym_frame_type, sym_frame_vjp, &
      & frame_level_none, frame_level_eigen, frame_level_plane, frame_level_sphere, frame_level_axial, &
      & frame_level_name, frame_level_codes, sym_threshold_type, new_sym_frame
   use moist_math_symmetry_embed, only: is_oh_element, site_candidates, neighbour_form_type, &
      & new_neighbour_form, compare_forms
   use moist_math_symmetry_orient, only: sym_orient_type, new_sym_orient, sym_orient_vjp, &
      & sym_orient_jacobian, orient_lab, orient_molecular, orient_codes
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_io, only: structure_type
   use test_cavity_symmetry_fixtures, only: structure_of, molecule_symmetry
   use moist_math_linalg_geometry, only: identity3, cross_product
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: collect_math_symmetry

   !> Exactness tolerance of the symmetric tests, bohr
   real(wp), parameter :: TOL = 1.0e-8_wp
   !> Central-difference step, bohr
   real(wp), parameter :: STEP = 1.0e-5_wp

contains

   !> Collect all cavity symmetry tests
   subroutine collect_math_symmetry(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("detect_reference_groups", test_detect_reference_groups), &
                  new_unittest("detect_linear_groups", test_detect_linear), &
                  new_unittest("detect_linear_candidates", test_detect_linear_candidates), &
                  new_unittest("detect_nearly_linear", test_detect_nearly_linear), &
                  new_unittest("group_used_auto", test_group_used_auto), &
                  new_unittest("group_used_requested", test_group_used_requested), &
                  new_unittest("site_spheres_invariant", test_site_spheres_invariant), &
                  new_unittest("transported_spheres_consistent", test_transport_consistent), &
                  new_unittest("rotations_deduplicated", test_rotations_deduplicated), &
                  new_unittest("rotations_modulo_oh", test_rotations_modulo_oh), &
                  new_unittest("rotations_modulo_oh_covariant", test_rotations_modulo_oh_covariant), &
                  new_unittest("lab_and_molecular_modes", test_lab_and_molecular), &
                  new_unittest("explicit_c1_like_off", test_explicit_c1), &
                  new_unittest("orbits", test_orbits), &
                  new_unittest("frame_levels", test_frame_levels), &
                  new_unittest("frame_covariant", test_frame_covariant), &
                  new_unittest("frame_linear_anywhere", test_frame_linear), &
                  new_unittest("frame_vjp_axial_stretch", test_frame_vjp_axial_stretch), &
                  new_unittest("frame_vjp_axial_rotation", test_frame_vjp_axial_rotation), &
                  new_unittest("symmetric_covariant", test_symmetric_covariant), &
                  new_unittest("neighbour_forms", test_neighbour_forms), &
                  new_unittest("requested_group_absent", test_requested_absent), &
                  new_unittest("requested_group_inexact", test_requested_inexact), &
                  new_unittest("requested_subgroup_messages", test_requested_messages), &
                  new_unittest("single_atom_groups", test_single_atom_groups), &
                  new_unittest("linear_subgroups", test_linear_subgroups), &
                  new_unittest("point_group_coverage", test_point_groups), &
                  new_unittest("frame_harmonic_flags", test_frame_harmonic_flags), &
                  new_unittest("auto_inexact_falls_back", test_auto_inexact), &
                  new_unittest("auto_loose_tolerance", test_auto_loose), &
                  new_unittest("refine_large_clusters", test_refine_large), &
                  new_unittest("detect_large_clusters", test_detect_large), &
                  new_unittest("frame_vjp_eigen", test_frame_vjp_eigen), &
                  new_unittest("frame_vjp_plane", test_frame_vjp_plane), &
                  new_unittest("orient_vjp_symmetric", test_orient_vjp_symmetric), &
                  new_unittest("orient_jacobian_matches_vjp", test_orient_jacobian), &
                  new_unittest("orient_jacobian_lab_angle", test_orient_lab_angle), &
                  new_unittest("frame_gap_line", test_frame_gap_line), &
                  new_unittest("frame_single_tops", test_frame_single_tops), &
                  new_unittest("frame_near_switch", test_frame_near_switch), &
                  new_unittest("frame_tied_stable", test_frame_tied_stable), &
                  new_unittest("orbit_ties_follow", test_orbit_ties_follow) &
                  ]

   end subroutine collect_math_symmetry

   !* ================================================================================= *!
   !*                                     Detection                                     *!
   !* ================================================================================= *!

   !> libmsym finds the textbook groups of reference molecules
   subroutine test_detect_reference_groups(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call expect_detected(error, "water", "C2v")
      if (allocated(error)) return
      call expect_detected(error, "ammonia", "C3v")
      if (allocated(error)) return
      call expect_detected(error, "methane", "Td")
      if (allocated(error)) return
      call expect_detected(error, "benzene", "D6h")
      if (allocated(error)) return
      call expect_detected(error, "bf3", "D3h")
      if (allocated(error)) return
      call expect_detected(error, "sf6", "Oh")
      if (allocated(error)) return
      call expect_detected(error, "icosahedron", "Ih")
      if (allocated(error)) return
      call expect_detected(error, "c1", "C1")
   end subroutine test_detect_reference_groups

   !> Linear molecules are flagged linear and carry no subgroup list
   subroutine test_detect_linear(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(msym_group_type) :: group
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)

      call molecule("co2", xyz, z)
      call find_point_group(xyz, z, group, merr)
      call check(error, .not. allocated(merr), "CO2 detection failed")
      if (allocated(error)) return
      call check(error, group%linear, "CO2 not flagged linear, got "//group%name)
      if (allocated(error)) return
      call check(error, size(group%subgroups), 0)
      if (allocated(error)) return

      call molecule("hcn", xyz, z)
      call find_point_group(xyz, z, group, merr)
      call check(error, .not. allocated(merr), "HCN detection failed")
      if (allocated(error)) return
      call check(error, group%linear, "HCN not flagged linear, got "//group%name)
   end subroutine test_detect_linear

   !> Detection at a tolerance ranks linear groups by their operation count
   !>
   !> - an exact linear molecule keeps its group at the tighter threshold
   !> - one CO bond stretched beyond the tolerance: C_inf_v at the tighter
   !>   threshold, D_inf_h at the default one, which is the larger candidate
   subroutine test_detect_linear_candidates(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(msym_group_type) :: group
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)

      call molecule("co2", xyz, z)
      call find_point_group(xyz, z, group, merr, TOL)
      call check(error, .not. allocated(merr), "CO2 detection at a tolerance failed")
      if (allocated(error)) return
      call check(error, group%linear .and. group%name == "D0h" .and. size(group%ops, 3) == 4, &
         & "CO2 at a tolerance: detected "//group%name)
      if (allocated(error)) return

      call molecule("hcn", xyz, z)
      call find_point_group(xyz, z, group, merr, TOL)
      call check(error, .not. allocated(merr), "HCN detection at a tolerance failed")
      if (allocated(error)) return
      call check(error, group%linear .and. group%name == "C0v" .and. size(group%ops, 3) == 2, &
         & "HCN at a tolerance: detected "//group%name)
      if (allocated(error)) return

      call molecule("co2", xyz, z)
      xyz(3, 2) = xyz(3, 2) + 1.0e-5_wp*aatoau
      call find_point_group(xyz, z, group, merr, TOL)
      call check(error, .not. allocated(merr), "stretched CO2: detection failed")
      if (allocated(error)) return
      call check(error, group%linear .and. group%name == "D0h" .and. size(group%ops, 3) == 4, &
         & "stretched CO2: detected "//group%name)
   end subroutine test_detect_linear_candidates

   !> A bent chain libmsym reads as linear and then rejects is detected as bent
   !>
   !> - libmsym merges the axes through the atoms into one C_inf axis up to
   !>   0.045 rad, then fails its permutation test: "Invalid equivalence set"
   !> - an unsymmetric triatomic bent by 0.005 bohr and a chain of five atoms
   !>   bent by 0.002 bohr; each fails at the default angle threshold
   !> - both are exactly Cs: the group is found, requested or not
   subroutine test_detect_nearly_linear(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(msym_group_type) :: group
      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp) :: xyz(3, 3), chain(3, 5)
      character(len=4), parameter :: requests(2) = [character(len=4) :: "Cs", "auto"]
      integer :: k

      xyz = reshape([-2.2_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.005_wp, 0.0_wp, 2.5_wp, 0.0_wp, 0.0_wp], [3, 3])
      call find_point_group(xyz, [8, 6, 8], group, merr)
      call check(error, .not. allocated(merr), "bent triatomic: detection failed")
      if (allocated(error)) return
      call check(error, group%name == "Cs", "bent triatomic: detected "//group%name)
      if (allocated(error)) return
      do k = 1, size(requests)
         call new_sym_orient(orient, xyz, [8, 6, 8], orient_molecular, trim(requests(k)), TOL, merr)
         call check(error, .not. allocated(merr), "bent triatomic with "//trim(requests(k))//" refused")
         if (allocated(error)) return
         call check(error, orient%group_used == "Cs" .and. orient%deviation <= TOL .and. &
            & len(orient%reason) == 0, "bent triatomic with "//trim(requests(k))//": used "// &
            & orient%group_used//", reason '"//orient%reason//"'")
         if (allocated(error)) return
      end do

      chain = reshape([-4.4_wp, 0.0_wp, 0.0_wp, -2.2_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.002_wp, 0.0_wp, &
         & 2.3_wp, 0.0_wp, 0.0_wp, 4.6_wp, 0.0_wp, 0.0_wp], [3, 5])
      call new_sym_orient(orient, chain, [8, 6, 6, 6, 8], orient_molecular, "auto", TOL, merr)
      call check(error, .not. allocated(merr), "bent chain refused")
      if (allocated(error)) return
      call check(error, orient%group_detected == "Cs" .and. orient%group_used == "Cs", &
         & "bent chain: detected "//orient%group_detected//", reason '"//orient%reason//"'")
   end subroutine test_detect_nearly_linear

   !* ================================================================================= *!
   !*                                    Group used                                     *!
   !* ================================================================================= *!

   !> Auto mode keeps every group whose sites embed in Oh and reduces the rest
   subroutine test_group_used_auto(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call expect_used(error, "water", "auto", "C2v", .false.)
      if (allocated(error)) return
      call expect_used(error, "ammonia", "auto", "C3v", .false.)
      if (allocated(error)) return
      call expect_used(error, "methane", "auto", "Td", .false.)
      if (allocated(error)) return
      ! Every benzene atom has site group C2v, so D6h itself is admissible
      call expect_used(error, "benzene", "auto", "D6h", .false.)
      if (allocated(error)) return
      call expect_used(error, "sf6", "auto", "Oh", .false.)
      if (allocated(error)) return
      ! Boron's D3h site has an S3: C3v keeps the planes through the B-F bonds
      call expect_used(error, "bf3", "auto", "C3v", .true.)
      if (allocated(error)) return
      ! Vertex sites C5v: Th is the largest subgroup of Ih inside Oh
      call expect_used(error, "icosahedron", "auto", "Th", .true.)
      if (allocated(error)) return
      call expect_used(error, "co2", "auto", "D4h", .true.)
      if (allocated(error)) return
      call expect_used(error, "hcn", "auto", "C4v", .true.)
      if (allocated(error)) return
      call expect_used(error, "c1", "auto", "C1", .false.)
   end subroutine test_group_used_auto

   !> A requested symbol selects that group or the best subgroup of it
   subroutine test_group_used_requested(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call expect_used(error, "water", "C2v", "C2v", .false.)
      if (allocated(error)) return
      call expect_used(error, "water", "cs", "Cs", .false.)
      if (allocated(error)) return
      call expect_used(error, "benzene", "D2h", "D2h", .false.)
      if (allocated(error)) return
      call expect_used(error, "bf3", "D3h", "C3v", .true.)
      if (allocated(error)) return
      call expect_used(error, "co2", "Dinfh", "D4h", .true.)
   end subroutine test_group_used_requested

   !* ================================================================================= *!
   !*                                     Exactness                                     *!
   !* ================================================================================= *!

   !> R_a^T s R_a is in Oh for every operation s fixing atom a
   subroutine test_site_spheres_invariant(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=12), parameter :: names(10) = [character(len=12) :: "water", "ammonia", &
         & "methane", "benzene", "bf3", "sf6", "icosahedron", "co2", "hcn", "oh48"]
      type(sym_orient_type) :: orient
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      real(wp) :: r(3, 3)
      integer :: im, ia, ig

      do im = 1, size(names)
         call orient_of(error, names(im), "auto", xyz, z, orient)
         if (allocated(error)) return
         do ia = 1, size(xyz, 2)
            r = orient%rotation(:, :, orient%rot_of_atom(ia))
            do ig = 1, size(orient%ops, 3)
               if (orient%perm(ia, ig) /= ia) cycle
               call check(error, is_oh_element(matmul(transpose(r), matmul(orient%ops(:, :, ig), r))), &
                  & trim(names(im))//": sphere not invariant under its site group")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_site_spheres_invariant

   !> The sphere of g(a) is det(g) g times the sphere of a, modulo Oh
   subroutine test_transport_consistent(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=12), parameter :: names(9) = [character(len=12) :: "water", "ammonia", &
         & "methane", "benzene", "bf3", "sf6", "icosahedron", "co2", "oh48"]
      type(sym_orient_type) :: orient
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      real(wp) :: ra(3, 3), rb(3, 3), g(3, 3)
      integer :: im, ia, ig

      do im = 1, size(names)
         call orient_of(error, names(im), "auto", xyz, z, orient)
         if (allocated(error)) return
         do ia = 1, size(xyz, 2)
            ra = orient%rotation(:, :, orient%rot_of_atom(ia))
            do ig = 1, size(orient%ops, 3)
               g = orient%ops(:, :, ig)
               rb = orient%rotation(:, :, orient%rot_of_atom(orient%perm(ia, ig)))
               call check(error, is_oh_element(matmul(transpose(rb), matmul(g, ra))), &
                  & trim(names(im))//": transported sphere is not the image of its source")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_transport_consistent

   !> Equal orientations are stored once: benzene's 12 atoms share 6
   subroutine test_rotations_deduplicated(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: orient
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)

      call orient_of(error, "benzene", "auto", xyz, z, orient)
      if (allocated(error)) return
      call check(error, orient%nrot <= 6, "benzene needs at most six orientations")
      if (allocated(error)) return

      call orient_of(error, "c1", "off", xyz, z, orient)
      if (allocated(error)) return
      call check(error, orient%nrot, 1)
   end subroutine test_rotations_deduplicated

   !> Orientations equal modulo Oh with equal derivative maps share one entry
   !>
   !> - a Lebedev sphere only sees R modulo Oh, and R_b = R_a o with constant o
   !>   has the left angular velocity of R_a: same points, moving alike
   !> - SF6: seven cube-aligned spheres, one entry
   !> - oh48: 29 exactly distinct orientations, at most 21 modulo Oh; the
   !>   24 atoms in the diagonal mirrors turn about their mirror normals and
   !>   stay apart by their derivative maps
   !> - minimal: no two entries agree modulo Oh with equal maps
   !>
   !> @param[out] error Error handle
   subroutine test_rotations_modulo_oh(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=12), parameter :: names(6) = [character(len=12) :: "water", "benzene", &
         & "sf6", "icosahedron", "co2", "oh48"]
      type(sym_orient_type) :: orient
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      integer :: im, ir, jr

      call orient_of(error, "sf6", "auto", xyz, z, orient)
      if (allocated(error)) return
      call check(error, orient%group_used == "Oh" .and. orient%nrot == 1, &
         & "sf6: cube-aligned spheres are not one entry")
      if (allocated(error)) return

      call orient_of(error, "oh48", "auto", xyz, z, orient)
      if (allocated(error)) return
      call check(error, orient%group_used == "Oh" .and. orient%nrot <= 21, &
         & "oh48: more than 21 entries modulo Oh")
      if (allocated(error)) return
      call check(error, any(orient%dw /= 0.0_wp), "oh48: no entry carries a derivative map")
      if (allocated(error)) return

      do im = 1, size(names)
         call orient_of(error, names(im), "auto", xyz, z, orient)
         if (allocated(error)) return
         call check(error, all(orient%rot_of_atom >= 1 .and. orient%rot_of_atom <= orient%nrot) .and. &
            & size(orient%rotation, 3) == orient%nrot .and. size(orient%dw, 3) == orient%nrot, &
            & trim(names(im))//": entry indices out of range")
         if (allocated(error)) return
         do ir = 1, orient%nrot
            call check(error, any(orient%rot_of_atom == ir), trim(names(im))//": unused entry")
            if (allocated(error)) return
            do jr = ir + 1, orient%nrot
               call check(error, .not. (is_oh_element(matmul(transpose(orient%rotation(:, :, ir)), &
                  & orient%rotation(:, :, jr))) .and. &
                  & all(abs(orient%dw(:, :, ir) - orient%dw(:, :, jr)) < 1.0e-6_wp)), &
                  & trim(names(im))//": two entries agree modulo Oh")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_rotations_modulo_oh

   !> Which atoms share an entry follows the atoms, up to an operation
   !>
   !> oh48 rotated, translated and listed in reverse: atoms a and b share an
   !> entry exactly when their images do, for one operation g of the group
   !> (a tie of the discrete choices gives the image under g)
   !>
   !> - the sharing only; [[test_symmetric_covariant]] checks the orientations
   !> - a failure lists the entries and orbits of both copies, the second
   !>   through the operation with the fewest differing pairs
   !>
   !> @param[out] error Error handle
   subroutine test_rotations_modulo_oh_covariant(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: o1, o2
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :), moved(:, :), ops(:, :, :)
      integer, allocatable :: z(:), perm(:, :)
      real(wp) :: q(3, 3)
      integer :: nat, ig, ia, ib, ga, gb, nbad, best, best_bad
      !> Diagnostic on failure
      character(len=:), allocatable :: more

      q = rotation_matrix(0.7_wp, -1.3_wp, 0.25_wp)
      call orient_of(error, "oh48", "auto", xyz, z, o1)
      if (allocated(error)) return
      call symmetry_of(xyz, z, ops, perm)
      nat = size(xyz, 2)
      moved = transformed(xyz, q, [0.9_wp, -1.4_wp, 0.3_wp])
      call new_sym_orient(o2, moved(:, nat:1:-1), z(nat:1:-1), orient_molecular, "auto", TOL, merr)
      call check(error, .not. allocated(merr), "oh48: rotated copy failed")
      if (allocated(error)) return
      call check(error, o2%group_used == o1%group_used .and. o2%nrot == o1%nrot, &
         & "oh48: group or number of entries changed")
      if (allocated(error)) return

      ! Operation with the fewest atom pairs whose sharing differs
      best = 0
      best_bad = huge(best_bad)
      do ig = 1, size(ops, 3)
         nbad = 0
         do ia = 1, nat
            ga = nat + 1 - perm(ia, ig)
            do ib = ia + 1, nat
               gb = nat + 1 - perm(ib, ig)
               if ((o1%rot_of_atom(ia) == o1%rot_of_atom(ib)) .neqv. &
                  & (o2%rot_of_atom(ga) == o2%rot_of_atom(gb))) nbad = nbad + 1
            end do
         end do
         if (nbad < best_bad) then
            best = ig
            best_bad = nbad
         end if
      end do
      if (best_bad == 0) return

      ! Entries of both copies, the second through the closest operation
      more = "closest operation "//to_string(best)//" of "//to_string(size(ops, 3))//": " &
         & //to_string(best_bad)//" atom pairs differ"//new_line("a")//"entry (o1|o2 g):"
      do ia = 1, nat
         more = more//" "//to_string(o1%rot_of_atom(ia))//"|"//to_string(o2%rot_of_atom(nat + 1 - perm(ia, best)))
      end do
      more = more//new_line("a")//"orbit (o1|o2 g):"
      do ia = 1, nat
         more = more//" "//to_string(o1%orbit_of_atom(ia))//"|"//to_string(o2%orbit_of_atom(nat + 1 - perm(ia, best)))
      end do
      more = more//new_line("a")//"frame level "//to_string(o1%frame%level)//"|"//to_string(o2%frame%level) &
         & //", gap "//to_string(o1%frame%gap, "(es10.3)")//"|"//to_string(o2%frame%gap, "(es10.3)") &
         & //", ambiguous "//merge("T", "F", o1%ambiguous)//"|"//merge("T", "F", o2%ambiguous) &
         & //", deviation "//to_string(o1%deviation, "(es10.3)")//"|"//to_string(o2%deviation, "(es10.3)")
      call test_failed(error, "oh48: shared entries do not follow the atoms", more)
   end subroutine test_rotations_modulo_oh_covariant

   !> Lab mode keeps the identity, molecular mode the frame, both one rotation
   subroutine test_lab_and_molecular(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)

      call molecule("c1", xyz, z)
      call new_sym_orient(orient, xyz, z, orient_lab, "off", TOL, merr)
      call check(error, .not. allocated(merr), "lab mode failed")
      if (allocated(error)) return
      call check(error, all(abs(orient%rotation(:, :, 1) - identity3) < 1.0e-14_wp), &
         & "lab mode is not the identity")
      if (allocated(error)) return
      call check(error, all(orient%dw == 0.0_wp), "lab mode has a frame derivative")
      if (allocated(error)) return

      call new_sym_orient(orient, xyz, z, orient_molecular, "off", TOL, merr)
      call check(error, .not. allocated(merr), "molecular mode failed")
      if (allocated(error)) return
      call check(error, all(abs(orient%rotation(:, :, 1) - orient%frame%axes) < 1.0e-14_wp), &
         & "molecular mode is not the frame")
      if (allocated(error)) return
      call check(error, abs(det(orient%rotation(:, :, 1)) - 1.0_wp) < 1.0e-12_wp, &
         & "frame is not a proper rotation")
   end subroutine test_lab_and_molecular

   !* ================================================================================= *!
   !*                                 Molecular frame                                   *!
   !* ================================================================================= *!

   !> The frame resolves at the level its second moment allows
   subroutine test_frame_levels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_frame_type) :: frame
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      !> Listing rebuilt from the codes
      character(len=:), allocatable :: expected
      integer :: level

      call molecule("c1", xyz, z)
      call default_frame(frame, xyz)
      call check(error, frame%level, frame_level_eigen)
      if (allocated(error)) return

      call molecule("ammonia", xyz, z)
      call default_frame(frame, xyz)
      call check(error, frame%level, frame_level_plane)
      if (allocated(error)) return
      call check(error, frame%harmonic, 3)
      if (allocated(error)) return

      call molecule("co2", xyz, z)
      call default_frame(frame, xyz)
      call check(error, frame%level, frame_level_axial)
      if (allocated(error)) return

      call check(error, frame_level_name(frame_level_none) == "none" .and. &
         & frame_level_name(frame_level_eigen) == "eigen" .and. &
         & frame_level_name(frame_level_plane) == "plane" .and. &
         & frame_level_name(frame_level_sphere) == "sphere" .and. &
         & frame_level_name(frame_level_axial) == "axial" .and. &
         & frame_level_name(-1) == "unknown", "unexpected frame_level_name")
      if (allocated(error)) return

      ! The listings in the field descriptions follow the codes and names
      expected = ""
      do level = frame_level_none, frame_level_axial
         if (level > frame_level_none) expected = expected//", "
         expected = expected//to_string(level)//" "//trim(frame_level_name(level))
      end do
      call check(error, frame_level_codes, expected, "frame_level_codes out of step with the levels")
      if (allocated(error)) return
      call check(error, orient_codes, to_string(orient_lab)//" lab frame, "// &
         & to_string(orient_molecular)//" molecular frame", "orient_codes out of step with the orientations")
   end subroutine test_frame_levels

   !> Linear molecules reach `axial` wherever they lie
   !>
   !> - shifted and tilted copies: the transverse coordinates are roundoff,
   !>   not structure, so no harmonic may fix the in-plane angle
   !> - the third axis is the line; the frame derivative stays of normal size
   !> - a bend of 1e-6 rad stays below the absolute line threshold: `axial`; a
   !>   bend of 0.3 rad is structure: `eigen`
   !>
   !> @param[out] error Error handle
   subroutine test_frame_linear(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=9), parameter :: names(3) = [character(len=9) :: "co2", "hcn", "acetylene"]
      type(sym_frame_type) :: frame
      real(wp), allocatable :: xyz(:, :), moved(:, :), grad(:, :)
      integer, allocatable :: z(:)
      real(wp) :: q(3, 3), shift(3), line(3), w(3, 3)
      integer :: im, it

      w = weight_matrix()
      do im = 1, size(names)
         call molecule(trim(names(im)), xyz, z)
         do it = 1, 4
            select case (it)
            case (1)
               q = identity3
               shift = 0.0_wp
            case (2)
               q = identity3
               shift = [0.1_wp, 0.2_wp, -0.3_wp]
            case (3)
               q = rotation_matrix(0.3_wp, 1.1_wp, -0.7_wp)
               shift = [0.1_wp, 0.2_wp, -0.3_wp]
            case default
               q = rotation_matrix(-2.1_wp, 0.4_wp, 1.9_wp)
               shift = [1.0e4_wp, -3.0e3_wp, 7.0e3_wp]
            end select
            moved = transformed(xyz, q, shift)
            call default_frame(frame, moved)
            call check(error, frame%level, frame_level_axial, trim(names(im))//": not axial")
            if (allocated(error)) return
            call check(error, frame%harmonic, 0, trim(names(im))//": a harmonic of roundoff")
            if (allocated(error)) return
            line = moved(:, 2) - moved(:, 1)
            line = line/norm2(line)
            call check(error, abs(abs(dot_product(frame%axes(:, 3), line)) - 1.0_wp) <= 1.0e-10_wp, &
               & trim(names(im))//": third axis is not the line")
            if (allocated(error)) return
            allocate (grad(3, size(moved, 2)))
            call sym_frame_vjp(frame, moved, oh_torque(frame%axes, w), grad)
            call check(error, all(abs(grad) <= 1.0e2_wp), trim(names(im))//": frame derivative blows up")
            if (allocated(error)) return
            deallocate (grad)
         end do
      end do

      ! Bent by 1e-6 rad: below the line threshold, which is absolute
      call molecule("co2", xyz, z)
      xyz(1, 2) = xyz(1, 2) + 1.0e-6_wp*xyz(3, 2)
      call default_frame(frame, xyz)
      call check(error, frame%level, frame_level_axial, "CO2 bent by 1e-6 rad is no line")
      if (allocated(error)) return
      ! Bent by 0.3 rad: structure, three separated eigenvalues
      xyz(1, 2) = xyz(1, 2) + 0.3_wp*xyz(3, 2)
      call default_frame(frame, xyz)
      call check(error, frame%level, frame_level_eigen, "bent CO2 is still a line")
   end subroutine test_frame_linear

   !> Axial frame derivative under stretches, which keep the molecule linear
   !>
   !> The frame does not move, so both sides vanish
   !>
   !> @param[out] error Error handle
   subroutine test_frame_vjp_axial_stretch(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_axial_vjp(error, .false.)
   end subroutine test_frame_vjp_axial_stretch

   !> Axial frame derivative under rigid rotations
   !>
   !> Covers the twist of the lab reference direction of the in-plane angle
   !>
   !> @param[out] error Error handle
   subroutine test_frame_vjp_axial_rotation(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_axial_vjp(error, .true.)
   end subroutine test_frame_vjp_axial_rotation

   !> Axial frame VJP against central differences on tilted HCN
   !>
   !> Only displacements that keep the molecule linear: any bend crosses the
   !> axial/plane switch
   !>
   !> @param[out] error    Error handle
   !> @param[in]  rotation Rigid rotations if true, stretches along the line else
   subroutine check_axial_vjp(error, rotation)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Rigid rotations or stretches
      logical, intent(in) :: rotation

      type(sym_frame_type) :: frame
      real(wp), allocatable :: xyz(:, :), grad(:, :), mode(:, :), d(:, :)
      integer, allocatable :: z(:)
      real(wp) :: w(3, 3), line(3), omega(3), ep, em, analytic, numeric
      real(wp), parameter :: coef(3, 3) = reshape([1.0_wp, 0.0_wp, 0.0_wp, 0.3_wp, -0.8_wp, 0.5_wp, &
         & -0.2_wp, 0.4_wp, 0.9_wp], [3, 3])
      integer :: imode, ia

      w = weight_matrix()
      call molecule("hcn", xyz, z)
      xyz = transformed(xyz, rotation_matrix(0.3_wp, 1.1_wp, -0.7_wp), [0.1_wp, 0.2_wp, -0.3_wp])
      call default_frame(frame, xyz)
      call check(error, frame%level, frame_level_axial)
      if (allocated(error)) return
      allocate (grad(3, 3), mode(3, 3))
      call sym_frame_vjp(frame, xyz, oh_torque(frame%axes, w), grad)
      line = (xyz(:, 3) - xyz(:, 1))/norm2(xyz(:, 3) - xyz(:, 1))
      d = xyz - spread(sum(xyz, dim=2)/3.0_wp, 2, 3)

      do imode = 1, 3
         do ia = 1, 3
            if (rotation) then
               omega = coef(:, imode)
               mode(:, ia) = cross_product(omega, d(:, ia))
            else
               mode(:, ia) = coef(ia, imode)*line
            end if
         end do
         analytic = sum(grad*mode)
         call default_frame(frame, xyz + STEP*mode)
         call check(error, frame%level, frame_level_axial)
         if (allocated(error)) return
         ep = oh_value(frame%axes, w)
         call default_frame(frame, xyz - STEP*mode)
         em = oh_value(frame%axes, w)
         numeric = (ep - em)/(2.0_wp*STEP)
         call check(error, analytic, numeric, thr=1.0e-7_wp)
         if (allocated(error)) return
      end do
   end subroutine check_axial_vjp

   !> Rotating, translating and permuting the atoms rotates the frame with them
   !>
   !> Up to a symmetry operation of the atoms: at exact symmetry, the frame's
   !> discrete choices tie between symmetry-related candidates
   subroutine test_frame_covariant(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=12), parameter :: names(5) = [character(len=12) :: "c1", "water", "ammonia", &
         & "benzene", "methane"]
      type(sym_frame_type) :: f1, f2
      real(wp), allocatable :: xyz(:, :), moved(:, :), ops(:, :, :)
      integer, allocatable :: z(:), perm(:, :)
      real(wp) :: q(3, 3)
      integer :: im, ig
      logical :: follows

      q = rotation_matrix(0.3_wp, -1.1_wp, 2.4_wp)
      do im = 1, size(names)
         call molecule(names(im), xyz, z)
         call symmetry_of(xyz, z, ops, perm)
         call default_frame(f1, xyz)
         moved = transformed(xyz, q, [1.5_wp, -0.7_wp, 2.2_wp])
         call default_frame(f2, moved(:, size(xyz, 2):1:-1))
         follows = .false.
         do ig = 1, size(ops, 3)
            follows = follows .or. is_oh_element(matmul(transpose(f2%axes), &
               & matmul(q, det(ops(:, :, ig))*matmul(ops(:, :, ig), f1%axes))))
         end do
         call check(error, follows, trim(names(im))//": frame does not follow the molecule")
         if (allocated(error)) return
      end do
   end subroutine test_frame_covariant

   !> Symmetric orientations follow the molecule, atom by atom
   !>
   !> Every molecule rotated, translated and listed in another order, in
   !> several trials: one symmetry M of the atoms, the same for every atom,
   !> maps the original spheres onto the new ones modulo Oh,
   !> R'_k = Q M R_a o_a for the atom a listed at k
   !>
   !> - M any symmetry of the full molecule, not only of the group used: a
   !>   different choice between images under M gives a congruent cavity
   !> - linear molecules: M includes the rotations about the axis (their
   !>   perpendicular axes come from the lab frame)
   !> - the discrete choices tie at exact symmetry: C3v sites (ammonia,
   !>   CH3F, ethane, PF5), the C2v sites of oh48, the Cs sites of the
   !>   icosahedron (Th used), the C3v subgroups of C6v and the D3d subgroups
   !>   of D6h (a centred sandwich)
   !> - symmetry-related ties left on the main axis, several of them with
   !>   different site groups (a sandwich with axial ligands): one symmetry
   !>   must relate all of them
   !>
   !> @param[out] error Error handle
   subroutine test_symmetric_covariant(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=12), parameter :: names(21) = [character(len=12) :: "water", "ammonia", &
         & "methane", "benzene", "bf3", "sf6", "icosahedron", "hocl", "co2", "hcn", "oh48", "ch3f", &
         & "ethane", "allene", "cyclopropane", "pf5", "xef4", "benzene_li", "benzene_lif", "sandwich", &
         & "sandwich_x2"]
      !> Trials: rotation, translation and atom order
      integer, parameter :: ntrial = 12
      type(sym_orient_type) :: o1, o2
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :), moved(:, :)
      integer, allocatable :: z(:), order(:)
      character(len=:), allocatable :: failed
      character(len=8) :: label
      real(wp) :: q(3, 3), shift(3)
      integer :: im, it, nat

      failed = ""
      do im = 1, size(names)
         call orient_of(error, names(im), "auto", xyz, z, o1)
         if (allocated(error)) return
         nat = size(xyz, 2)
         do it = 1, ntrial
            q = rotation_matrix(1.7_wp*real(it, wp) - 0.89_wp, 0.45_wp - 1.72_wp*real(it, wp), &
               & 2.33_wp - 0.83_wp*real(it**2, wp))
            shift = [-2.0_wp, 0.4_wp, 1.0_wp]*real(it, wp)
            order = trial_order(nat, it)
            moved = transformed(xyz(:, order), q, shift)
            call new_sym_orient(o2, moved, z(order), orient_molecular, "auto", TOL, merr)
            call check(error, .not. allocated(merr), trim(names(im))//": moved copy failed")
            if (allocated(error)) return
            call check(error, o2%group_used == o1%group_used, trim(names(im))//": group changed")
            if (allocated(error)) return
            if (.not. spheres_follow(xyz, z, o1, o2, q, order)) then
               write (label, "(i0)") it
               failed = failed//" "//trim(names(im))//"#"//trim(label)
            end if
         end do
      end do
      call check(error, len(failed) == 0, "spheres do not follow their atoms:"//failed)
   end subroutine test_symmetric_covariant

   !> Tied site orientations have equal neighbour forms when symmetry-related
   !>
   !> Equal forms exactly when a symmetry of the molecule fixing the site
   !> relates the two orientations
   !>
   !> - C3v sites with an exact tie in closeness to the frame, R and C2(u) R
   !>   about the axis u: no symmetry relates them on ammonia's N and CH3F's
   !>   C and F, sigma_h does on the centers of BF3 and PF5 (C3v used)
   !> - a form follows rigid motion and atom order
   !>
   !> @param[out] error Error handle
   subroutine test_neighbour_forms(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=12), parameter :: names(5) = [character(len=12) :: "ammonia", "ch3f", "ch3f", &
         & "bf3", "pf5"]
      !> Site atom of every case
      integer, parameter :: sites(5) = [1, 1, 2, 1, 1]
      !> Whether a symmetry of the molecule relates the tied pair
      logical, parameter :: related(5) = [.false., .false., .false., .true., .true.]
      type(sym_orient_type) :: orient
      type(neighbour_form_type) :: f1, f2
      real(wp), allocatable :: xyz(:, :), site(:, :, :), rot(:, :, :), dw(:, :, :), moved(:, :)
      integer, allocatable :: z(:), order(:)
      real(wp) :: value, q(3, 3), oh(3, 3, 48), m(3, 3)
      integer, allocatable :: image(:)
      integer :: im, ia, ig, io, isg, k
      logical :: found, by_symmetry

      call oh_group(oh)
      do im = 1, size(names)
         call orient_of(error, names(im), "auto", xyz, z, orient)
         if (allocated(error)) return
         ia = sites(im)
         site = orient%ops(:, :, pack([(ig, ig=1, size(orient%ops, 3))], orient%perm(ia, :) == ia))
         call site_candidates(site, orient%frame%axes, found, value, rot, dw)
         call check(error, found .and. size(rot, 3) == 2, trim(names(im))//": expected two tied orientations")
         if (allocated(error)) return

         ! Some orthogonal map R2 o R1^T (o in Oh, either sign) a symmetry fixing the site
         by_symmetry = .false.
         allocate (image(size(z)))
         do io = 1, 48
            do isg = -1, 1, 2
               m = real(isg, wp)*matmul(rot(:, :, 2), matmul(oh(:, :, io), transpose(rot(:, :, 1))))
               if (is_molecular_symmetry(xyz, z, m, image)) by_symmetry = by_symmetry .or. image(ia) == ia
            end do
         end do
         deallocate (image)
         call check(error, by_symmetry .eqv. related(im), trim(names(im))//": unexpected symmetry relation")
         if (allocated(error)) return

         call new_neighbour_form(f1, xyz, z, ia, rot(:, :, 1))
         call new_neighbour_form(f2, xyz, z, ia, rot(:, :, 2))
         call check(error, (compare_forms(f1, f2) == 0) .eqv. related(im), &
            & trim(names(im))//": forms equal without a symmetry, or differ with one")
         if (allocated(error)) return
         call check(error, compare_forms(f1, f2) == -compare_forms(f2, f1), &
            & trim(names(im))//": form order not antisymmetric")
         if (allocated(error)) return

         ! The same sphere on a moved, reordered copy
         q = rotation_matrix(0.4_wp, -2.1_wp, 1.3_wp)
         order = trial_order(size(z), 4)
         moved = transformed(xyz(:, order), q, [0.3_wp, 1.7_wp, -0.9_wp])
         k = findloc(order, ia, dim=1)
         call new_neighbour_form(f2, moved, z(order), k, matmul(q, rot(:, :, 1)))
         call check(error, compare_forms(f1, f2) == 0, trim(names(im))//": form does not follow the molecule")
         if (allocated(error)) return
      end do
   end subroutine test_neighbour_forms

   !> Whether the spheres of a moved copy are the images of the original ones
   !>
   !> Some orthogonal M maps the atoms onto atoms of their species (about the
   !> centroid), a onto M(a), and every sphere with them: R'_k = Q M R_a
   !> modulo Oh for the copy k of M(a); M is tried from every Oh image of the
   !> first atom's sphere onto every sphere of its species
   !>
   !> @param[in] xyz   Original coordinates (3, nat)
   !> @param[in] z     Species (nat)
   !> @param[in] o1    Orientations of the original
   !> @param[in] o2    Orientations of the moved copy
   !> @param[in] q     Rotation of the copy
   !> @param[in] order Original atom listed at each position of the copy (nat)
   function spheres_follow(xyz, z, o1, o2, q, order) result(follows)
      !> Original coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species
      integer, intent(in) :: z(:)
      !> Orientations of the original
      type(sym_orient_type), intent(in) :: o1
      !> Orientations of the moved copy
      type(sym_orient_type), intent(in) :: o2
      !> Rotation of the copy
      real(wp), intent(in) :: q(3, 3)
      !> Atom order of the copy
      integer, intent(in) :: order(:)
      !> Whether one symmetry maps every sphere
      logical :: follows

      real(wp) :: oh(3, 3, 48), m(3, 3), ra(3, 3), rb(3, 3)
      integer, allocatable :: at(:), image(:)
      integer :: nat, io, isg, ia, ib, k

      nat = size(xyz, 2)
      allocate (at(nat), image(nat))
      do k = 1, nat
         at(order(k)) = k
      end do
      call oh_group(oh)
      follows = .false.
      do ib = 1, nat
         if (z(ib) /= z(1)) cycle
         do io = 1, 48
            do isg = -1, 1, 2
               ra = o1%rotation(:, :, o1%rot_of_atom(1))
               rb = o2%rotation(:, :, o2%rot_of_atom(at(ib)))
               m = real(isg, wp)*matmul(transpose(q), matmul(rb, matmul(oh(:, :, io), transpose(ra))))
               if (.not. is_molecular_symmetry(xyz, z, m, image)) cycle
               if (image(1) /= ib) cycle
               follows = .true.
               do ia = 1, nat
                  ra = o1%rotation(:, :, o1%rot_of_atom(ia))
                  rb = o2%rotation(:, :, o2%rot_of_atom(at(image(ia))))
                  follows = follows .and. is_oh_element(matmul(transpose(rb), matmul(q, matmul(m, ra))))
               end do
               if (follows) return
            end do
         end do
      end do
   end function spheres_follow

   !> Whether an orthogonal map about the centroid is a symmetry of the atoms
   !>
   !> Every atom onto an atom of its species, to 1e-6 bohr
   !>
   !> @param[in]  xyz   Coordinates (3, nat)
   !> @param[in]  z     Species (nat)
   !> @param[in]  m     Orthogonal map
   !> @param[out] image Image atom of every atom (nat), when a symmetry
   function is_molecular_symmetry(xyz, z, m, image) result(yes)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species
      integer, intent(in) :: z(:)
      !> Orthogonal map
      real(wp), intent(in) :: m(3, 3)
      !> Image atoms
      integer, intent(out) :: image(:)
      !> Whether m is a symmetry
      logical :: yes

      real(wp) :: center(3), target(3)
      integer :: ia, ib

      center = sum(xyz, dim=2)/real(size(xyz, 2), wp)
      yes = .true.
      do ia = 1, size(xyz, 2)
         target = center + matmul(m, xyz(:, ia) - center)
         image(ia) = 0
         do ib = 1, size(xyz, 2)
            if (z(ib) == z(ia) .and. norm2(xyz(:, ib) - target) < 1.0e-6_wp) image(ia) = ib
         end do
         if (image(ia) == 0) then
            yes = .false.
            return
         end if
      end do
   end function is_molecular_symmetry

   !> Atom order of a covariance trial
   !>
   !> Reversed, cycled by one, then pseudo-random shuffles
   !>
   !> @param[in] nat   Number of atoms
   !> @param[in] trial Trial number
   pure function trial_order(nat, trial) result(order)
      !> Number of atoms
      integer, intent(in) :: nat
      !> Trial number
      integer, intent(in) :: trial
      !> Original atom at every position
      integer, allocatable :: order(:)

      integer(int64) :: state
      integer :: i, j, tmp

      order = [(i, i=1, nat)]
      select case (trial)
      case (1)
         order = order(nat:1:-1)
      case (2)
         order = cshift(order, 1)
      case default
         ! Fisher-Yates with a linear congruential generator
         state = int(trial, int64)*7919_int64
         do i = nat, 2, -1
            state = modulo(1103515245_int64*state + 12345_int64, 2147483648_int64)
            j = int(modulo(state, int(i, int64))) + 1
            tmp = order(i)
            order(i) = order(j)
            order(j) = tmp
         end do
      end select
   end function trial_order

   !> The 48 signed permutation matrices
   !>
   !> @param[out] oh Elements (3, 3, 48)
   pure subroutine oh_group(oh)
      !> Elements
      real(wp), intent(out) :: oh(3, 3, 48)

      integer, parameter :: perms(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, &
         & 3, 2, 1], [3, 6])
      integer :: ip, is, i, n

      n = 0
      oh = 0.0_wp
      do ip = 1, 6
         do is = 0, 7
            n = n + 1
            do i = 1, 3
               oh(perms(i, ip), i, n) = merge(-1.0_wp, 1.0_wp, btest(is, i - 1))
            end do
         end do
      end do
   end subroutine oh_group

   !> An explicit C1 is off: trivial orientation, no detection, never an error
   subroutine test_explicit_c1(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      character(len=*), parameter :: names(3) = [character(len=8) :: "water", "bf3", "c1"]
      integer :: k

      do k = 1, size(names)
         call molecule(trim(names(k)), xyz, z)
         call new_sym_orient(orient, xyz, z, orient_lab, "C1", TOL, merr)
         call check(error, .not. allocated(merr), trim(names(k))//": explicit C1 failed")
         if (allocated(error)) return
         call check(error, orient%group_used == "C1" .and. orient%group_requested == "C1" .and. &
            & len(orient%group_detected) == 0 .and. len(orient%reason) == 0, &
            & trim(names(k))//": explicit C1 is not reported as C1 without detection")
         if (allocated(error)) return
         call check(error, orient%orientation == orient_lab .and. orient%nrot == 1, &
            & trim(names(k))//": explicit C1 changed the orientation")
         if (allocated(error)) return
         call check(error, all(orient%rotation(:, :, 1) == identity3), &
            & trim(names(k))//": explicit C1 in lab mode is not the exact identity")
         if (allocated(error)) return
      end do

      ! Lower case, molecular orientation: the frame, as with off
      call molecule("water", xyz, z)
      call new_sym_orient(orient, xyz, z, orient_molecular, "c1", TOL, merr)
      call check(error, .not. allocated(merr), "explicit c1 in molecular mode failed")
      if (allocated(error)) return
      call check(error, orient%group_used == "C1" .and. orient%nrot == 1, "explicit c1 is not trivial")
      if (allocated(error)) return
      call check(error, all(abs(orient%rotation(:, :, 1) - orient%frame%axes) < 1.0e-14_wp), &
         & "explicit c1 in molecular mode is not the frame")
   end subroutine test_explicit_c1

   !> Orbits: the atoms the group used maps onto each other
   subroutine test_orbits(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      integer :: k

      call orient_of(error, "water", "auto", xyz, z, orient)
      if (allocated(error)) return
      call check(error, orient%norbit, 2, "water: two orbits")
      if (allocated(error)) return
      call check(error, all(orient%orbit_of_atom == [1, 2, 2]), "water: orbits {O}, {H1, H2}")
      if (allocated(error)) return
      call check_orbits(error, "water", orient)
      if (allocated(error)) return

      call orient_of(error, "benzene", "auto", xyz, z, orient)
      if (allocated(error)) return
      call check(error, orient%norbit, 2, "benzene: two orbits")
      if (allocated(error)) return
      call check(error, all(orient%orbit_of_atom == [(1, k=1, 6), (2, k=1, 6)]), &
         & "benzene: orbits {6 C}, {6 H}")
      if (allocated(error)) return
      call check_orbits(error, "benzene", orient)
      if (allocated(error)) return

      ! BF3 used as C3v: boron alone, the three fluorines together
      call orient_of(error, "bf3", "auto", xyz, z, orient)
      if (allocated(error)) return
      call check(error, all(orient%orbit_of_atom == [1, 2, 2, 2]), "bf3: orbits {B}, {3 F}")
      if (allocated(error)) return
      call check_orbits(error, "bf3", orient)
      if (allocated(error)) return

      ! Off: every atom its own orbit
      call molecule("benzene", xyz, z)
      call new_sym_orient(orient, xyz, z, orient_lab, "off", TOL, merr)
      call check(error, .not. allocated(merr), "benzene: off failed")
      if (allocated(error)) return
      call check(error, orient%norbit, size(z), "off: one orbit per atom")
      if (allocated(error)) return
      call check(error, all(orient%orbit_of_atom == [(k, k=1, size(z))]) .and. &
         & all(orient%orbit_rep == [(k, k=1, size(z))]), "off: atom k is orbit k")

   contains

      !> Every orbit is closed and its representative reaches every member
      !>
      !> @param[out] error  Error handle
      !> @param[in]  name   Molecule
      !> @param[in]  orient Orientations
      subroutine check_orbits(error, name, orient)
         !> Error handle
         type(error_type), allocatable, intent(out) :: error
         !> Molecule
         character(len=*), intent(in) :: name
         !> Orientations
         type(sym_orient_type), intent(in) :: orient

         integer :: ia, ig, rep

         call check(error, size(orient%orbit_rep), orient%norbit, name//": one representative per orbit")
         if (allocated(error)) return
         do ia = 1, size(orient%orbit_of_atom)
            rep = orient%orbit_rep(orient%orbit_of_atom(ia))
            call check(error, orient%orbit_of_atom(rep) == orient%orbit_of_atom(ia) .and. &
               & any(orient%perm(rep, :) == ia), name//": representative does not reach its orbit")
            if (allocated(error)) return
            do ig = 1, size(orient%perm, 2)
               call check(error, orient%orbit_of_atom(orient%perm(ia, ig)) == orient%orbit_of_atom(ia), &
                  & name//": an operation leaves the orbit")
               if (allocated(error)) return
            end do
         end do
      end subroutine check_orbits

   end subroutine test_orbits

   !* ================================================================================= *!
   !*                                      Errors                                       *!
   !* ================================================================================= *!

   !> Water has no C3v subgroup: the error names both groups
   subroutine test_requested_absent(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)

      call molecule("water", xyz, z)
      call new_sym_orient(orient, xyz, z, orient_molecular, "C3v", TOL, merr)
      call check(error, allocated(merr), "C3v accepted for water")
      if (allocated(error)) return
      call check(error, index(merr%message, "requested group C3v is not a subgroup of the detected group C2v") > 0, &
         & "wrong error for an absent group: "//merr%message)
   end subroutine test_requested_absent

   !> A distortion above the tolerance is an error for a requested group
   !>
   !> - the message carries the deviation of the requested group and the tolerance
   subroutine test_requested_inexact(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      real(wp) :: dev

      ! An in-plane stretch of one O-H bond; any triatomic stays exactly Cs
      call molecule("water", xyz, z)
      xyz(2, 2) = xyz(2, 2) + 1.0e-5_wp
      call new_sym_orient(orient, xyz, z, orient_molecular, "C2v", TOL, merr)
      call check(error, allocated(merr), "inexact C2v accepted for water")
      if (allocated(error)) return
      call check(error, index(merr%message, "geometry deviates from C2v by") > 0 .and. &
         & index(merr%message, "above the tolerance") > 0, "wrong error for an inexact group: "//merr%message)
      if (allocated(error)) return
      dev = message_deviation(merr%message)
      call check(error, dev > 1.0e-6_wp .and. dev < 2.0e-5_wp, "deviation of the stretch not reported: "// &
         & merr%message)
   end subroutine test_requested_inexact

   !> Deviation in bohr an error message or reason reports after " by " or " to "
   !>
   !> @param[in] message Message
   function message_deviation(message) result(dev)
      !> Message
      character(len=*), intent(in) :: message
      !> Deviation, negative when none is found
      real(wp) :: dev

      integer :: i, stat

      dev = -1.0_wp
      i = index(message, " by ")
      if (i == 0) i = index(message, " to ")
      if (i == 0) return
      read (message(i + 4:), *, iostat=stat) dev
      if (stat /= 0) dev = -1.0_wp
   end function message_deviation

   !> Explicit subgroups report their own deviation and the true reason
   !>
   !> - benzene stretched along one C2 axis by 1e-4 bohr (D2h kept) with one
   !>   atom moved by 1e-6 bohr: requesting D2h fails with a deviation of the
   !>   order of 1e-6, not that of D6h
   !> - benzene with an axial atom and the hydrogen ring turned by 1e-5 rad
   !>   (C6 kept, the mirror planes lost): requesting C6 succeeds with a
   !>   subgroup, and the reason is the site of the axial atom, not a deviation
   !>
   !> @param[out] error Error handle
   subroutine test_requested_messages(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      real(wp) :: dev, c, s
      integer :: k

      call molecule("benzene", xyz, z)
      xyz(1, [1, 7]) = xyz(1, [1, 7]) + 1.0e-4_wp
      xyz(1, [4, 10]) = xyz(1, [4, 10]) - 1.0e-4_wp
      xyz(:, 2) = xyz(:, 2) + [0.0_wp, 0.0_wp, 1.0e-6_wp]
      call new_sym_orient(orient, xyz, z, orient_molecular, "D2h", TOL, merr)
      call check(error, allocated(merr), "inexact D2h accepted for stretched benzene")
      if (allocated(error)) return
      dev = message_deviation(merr%message)
      call check(error, index(merr%message, "geometry deviates from D2h by") > 0 .and. &
         & dev > 1.0e-7_wp .and. dev < 1.0e-5_wp, "D2h reported with another group's deviation: "// &
         & merr%message)
      if (allocated(error)) return
      deallocate (merr)

      call molecule("benzene_li", xyz, z)
      c = cos(1.0e-5_wp)
      s = sin(1.0e-5_wp)
      do k = 7, 12
         xyz(1:2, k) = [c*xyz(1, k) - s*xyz(2, k), s*xyz(1, k) + c*xyz(2, k)]
      end do
      call new_sym_orient(orient, xyz, z, orient_molecular, "C6", TOL, merr)
      if (allocated(merr)) then
         call test_failed(error, "C6 of a twisted C6v molecule: "//merr%message)
         return
      end if
      call check(error, orient%group_used /= "C6" .and. orient%group_used /= "C1", &
         & "twisted C6v: expected a subgroup of C6, got "//orient%group_used)
      if (allocated(error)) return
      call check(error, index(orient%reason, "requested C6, using "//orient%group_used) > 0 .and. &
         & index(orient%reason, "site group") > 0 .and. index(orient%reason, "matches only") == 0, &
         & "twisted C6v: wrong reason: "//orient%reason)
   end subroutine test_requested_messages

   !> One atom: any point group is accepted and needs nothing, a symbol that
   !> names none is an error
   !>
   !> @param[out] error Error handle
   subroutine test_single_atom_groups(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=5), parameter :: good(13) = [character(len=5) :: "auto", "off", "C1", "D6h", "c2v", &
         & "Td", "S4", "Dinfh", "C0v", "d0h", "S10", "C3h", "D2"]
      ! No order, a zero or leading-zero order, odd or small Sn, C1v, D1
      character(len=5), parameter :: bad(15) = [character(len=5) :: "Bogus", "C", "D2x", "Sv", "3", &
         & "C0", "S3", "S2", "S0", "D0", "C0h", "D0d", "C02v", "C1v", "D1h"]
      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp) :: xyz(3, 1)
      integer :: k

      xyz(:, 1) = [0.3_wp, -0.2_wp, 0.5_wp]
      do k = 1, size(good)
         call new_sym_orient(orient, xyz, [8], orient_molecular, trim(good(k)), TOL, merr)
         call check(error, .not. allocated(merr), "one atom: "//trim(good(k))//" refused")
         if (allocated(error)) return
         call check(error, orient%group_used == "C1" .and. orient%nrot == 1, "one atom: not the trivial result")
         if (allocated(error)) return
         call check(error, (len(orient%reason) > 0) .eqv. (k > 3), &
            & "one atom: "//trim(good(k))//" has the wrong reason: "//orient%reason)
         if (allocated(error)) return
      end do
      do k = 1, size(bad)
         call new_sym_orient(orient, xyz, [8], orient_molecular, trim(bad(k)), TOL, merr)
         call check(error, allocated(merr), "one atom: "//trim(bad(k))//" accepted")
         if (allocated(error)) return
         call check(error, index(merr%message, "unknown point group") > 0, &
            & "one atom: wrong error: "//merr%message)
         if (allocated(error)) return
         deallocate (merr)
      end do
   end subroutine test_single_atom_groups

   !> Linear molecules take the finite subgroups of their stand-in group
   !>
   !> - CO2 (D4h): every subgroup named; HCN (C4v): those without an
   !>   operation that turns the axis over
   !> - the operations used match the name in order, map the atoms exactly
   !>   and leave every site sphere invariant
   !> - a finite group used as requested reports no reason; the linear
   !>   symbol reports its stand-in
   !>
   !> @param[out] error Error handle
   subroutine test_linear_subgroups(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=3), parameter :: groups(15) = [character(len=3) :: "D4h", "D2h", "D2d", "C4v", "C4h", &
         & "D4", "C2v", "C2h", "D2", "C4", "S4", "Cs", "Ci", "C2", "C3v"]
      integer, parameter :: order(15) = [16, 8, 8, 8, 8, 8, 4, 4, 4, 4, 4, 2, 2, 2, 0]
      !> Whether HCN has the group
      logical, parameter :: polar(15) = [.false., .false., .false., .true., .false., .false., .true., &
         & .false., .false., .true., .false., .true., .false., .true., .false.]
      character(len=3), parameter :: names(2) = ["co2", "hcn"]
      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      type(structure_type) :: mol
      integer, allocatable :: z(:)
      real(wp) :: image(3), center(3)
      integer :: im, ig, iop, ia
      logical :: expected

      do im = 1, 2
         call structure_of(names(im), mol)
         z = mol%num(mol%id)
         do ig = 1, size(groups)
            expected = order(ig) > 0 .and. (im == 1 .or. polar(ig))
            call new_sym_orient(orient, mol%xyz, z, orient_molecular, groups(ig), TOL, merr)
            call check(error, allocated(merr) .neqv. expected, names(im)//" with "//trim(groups(ig))// &
               & merge(": refused ", ": accepted", expected))
            if (allocated(error)) return
            if (.not. expected) then
               call check(error, index(merr%message, "is not a subgroup of the detected group") > 0, &
                  & names(im)//" with "//trim(groups(ig))//": wrong error: "//merr%message)
               if (allocated(error)) return
               deallocate (merr)
               cycle
            end if
            call check(error, orient%group_used == trim(groups(ig)) .and. size(orient%ops, 3) == order(ig), &
               & names(im)//" with "//trim(groups(ig))//": used "//orient%group_used)
            if (allocated(error)) return
            ! The request is honoured as it stands: nothing to report
            call check(error, len(orient%reason) == 0, names(im)//" with "//trim(groups(ig))// &
               & ": unexpected reason '"//orient%reason//"'")
            if (allocated(error)) return
            center = orient%center
            do iop = 1, size(orient%ops, 3)
               do ia = 1, mol%nat
                  image = center + matmul(orient%ops(:, :, iop), mol%xyz(:, ia) - center)
                  call check(error, norm2(image - mol%xyz(:, orient%perm(ia, iop))) < 1.0e-10_wp .and. &
                     & z(orient%perm(ia, iop)) == z(ia), names(im)//" with "//trim(groups(ig))// &
                     & ": an operation does not map the atoms")
                  if (allocated(error)) return
                  ! The sphere of the image atom is the image sphere, modulo Oh
                  call check(error, is_oh_element(matmul(transpose( &
                     & orient%rotation(:, :, orient%rot_of_atom(orient%perm(ia, iop)))), &
                     & matmul(orient%ops(:, :, iop), orient%rotation(:, :, orient%rot_of_atom(ia))))), &
                     & names(im)//" with "//trim(groups(ig))//": a sphere is not carried to its image")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do

      ! The linear symbol itself names no finite group: the stand-in is reported
      call structure_of("co2", mol)
      call new_sym_orient(orient, mol%xyz, mol%num(mol%id), orient_molecular, "D0h", TOL, merr)
      call check(error, .not. allocated(merr), "co2 with D0h refused")
      if (allocated(error)) return
      call check(error, index(orient%reason, "requested D0h, using "//orient%group_used// &
         & " (finite subgroup of a linear group)") == 1, "co2 with D0h: reason '"//orient%reason//"'")
   end subroutine test_linear_subgroups

   !> Further point groups: covariance, exact spheres and derivatives
   !>
   !> C2, C2h, Ci, D2, S4, D2d, C4v, C6v, Th, Ih, T, O, Oh and the linear
   !> groups, detected with "auto"
   !>
   !> - the group used is the expected one (the largest whose site groups
   !>   fit a cube)
   !> - rotating, shifting and renumbering the atoms carries the spheres along
   !> - every operation used carries the sphere of an atom to that of its image
   !> - the orientation Jacobian matches the VJP; a rigid rotation turns every
   !>   sphere with the molecule, also in the groups that keep a rotation
   !>   (C2, C2h, S4, Ci about their axis)
   !> - along a displacement that keeps the group, the VJP matches central
   !>   differences of an Oh-invariant function of the spheres
   !>
   !> @param[out] error Error handle
   subroutine test_point_groups(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ncase = 15
      character(len=14), parameter :: names(ncase) = [character(len=14) :: "h2o2", "h2o2_trans", "ci", &
         & "ethene_twisted", "s4", "allene", "pyramid", "c6v", "pyrite", "icosahedron", "t12", "o24", &
         & "cube", "co2", "hcn"]
      character(len=3), parameter :: detected(ncase) = [character(len=3) :: "C2", "C2h", "Ci", "D2", "S4", &
         & "D2d", "C4v", "C6v", "Th", "Ih", "T", "O", "Oh", "D0h", "C0v"]
      character(len=3), parameter :: used(ncase) = [character(len=3) :: "C2", "C2h", "Ci", "D2", "S4", &
         & "D2d", "C4v", "C3v", "Th", "Th", "T", "O", "Oh", "D4h", "C4v"]
      integer, parameter :: ntrial = 4
      type(sym_orient_type) :: o1, o2
      type(moist_error_type), allocatable :: merr
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :), moved(:, :), torque(:, :), grad(:, :), contracted(:, :), &
         & jac(:, :, :, :), mode(:, :), sym_mode(:, :)
      integer, allocatable :: z(:), order(:)
      !> Operations of the detected group and their atom images
      real(wp), allocatable :: full_ops(:, :, :)
      integer, allocatable :: full_perm(:, :)
      real(wp) :: q(3, 3), shift(3), w(3, 3), center(3), omega(3), rigid, ep, em
      integer :: im, it, nat, ia, ib, k, ig
      character(len=40) :: text

      w = weight_matrix()
      do im = 1, ncase
         call structure_of(trim(names(im)), mol)
         xyz = mol%xyz
         z = mol%num(mol%id)
         nat = mol%nat
         call new_sym_orient(o1, xyz, z, orient_molecular, "auto", TOL, merr)
         if (allocated(merr)) then
            call test_failed(error, trim(names(im))//": "//merr%message)
            return
         end if
         call check(error, o1%group_detected == trim(detected(im)) .and. o1%group_used == trim(used(im)), &
            & trim(names(im))//": detected "//o1%group_detected//", used "//o1%group_used)
         if (allocated(error)) return

         ! Exact spheres under every operation used
         do ig = 1, size(o1%ops, 3)
            do ia = 1, nat
               call check(error, norm2(o1%center + matmul(o1%ops(:, :, ig), xyz(:, ia) - o1%center) &
                  & - xyz(:, o1%perm(ia, ig))) < 1.0e-9_wp .and. &
                  & is_oh_element(matmul(transpose(o1%rotation(:, :, o1%rot_of_atom(o1%perm(ia, ig)))), &
                  & matmul(o1%ops(:, :, ig), o1%rotation(:, :, o1%rot_of_atom(ia))))), &
                  & trim(names(im))//": a sphere is not carried to its image")
               if (allocated(error)) return
            end do
         end do

         ! Covariance; the lab-derived angle of a line is no property of it
         if (o1%group_detected /= "D0h" .and. o1%group_detected /= "C0v") then
            do it = 1, ntrial
               q = rotation_matrix(1.7_wp*real(it, wp) - 0.89_wp, 0.45_wp - 1.72_wp*real(it, wp), &
                  & 2.33_wp - 0.83_wp*real(it**2, wp))
               shift = [-2.0_wp, 0.4_wp, 1.0_wp]*real(it, wp)
               order = trial_order(nat, it)
               moved = transformed(xyz(:, order), q, shift)
               call new_sym_orient(o2, moved, z(order), orient_molecular, "auto", TOL, merr)
               call check(error, .not. allocated(merr), trim(names(im))//": moved copy failed")
               if (allocated(error)) return
               call check(error, o2%group_used == o1%group_used .and. spheres_follow(xyz, z, o1, o2, q, order), &
                  & trim(names(im))//": spheres do not follow their atoms")
               if (allocated(error)) return
            end do
         end if

         ! Jacobian against the VJP, rigid rotations and translations
         allocate (torque(3, nat), grad(3, nat), contracted(3, nat), jac(3, 3, nat, o1%nrot), &
            & mode(3, nat), sym_mode(3, nat))
         do ia = 1, nat
            torque(:, ia) = sin(1.3_wp*real(ia, wp) + [0.2_wp, 1.9_wp, 3.7_wp])
            mode(:, ia) = sin(0.7_wp*real(ia, wp) + [1.1_wp, 2.9_wp, 0.3_wp])
         end do
         ! A line: the angle about it follows the lab frame, so only the
         ! torque across the line obeys the rule of rigid rotations
         if (o1%group_detected == "D0h" .or. o1%group_detected == "C0v") then
            do ia = 1, nat
               torque(:, ia) = torque(:, ia) - dot_product(torque(:, ia), line_axis(xyz))*line_axis(xyz)
            end do
         end if
         call sym_orient_vjp(o1, xyz, torque, grad)
         call sym_orient_jacobian(o1, xyz, jac)
         contracted = 0.0_wp
         do ia = 1, nat
            do ib = 1, nat
               do k = 1, 3
                  contracted(k, ib) = contracted(k, ib) + dot_product(torque(:, ia), jac(:, k, ib, o1%rot_of_atom(ia)))
               end do
            end do
         end do
         call check(error, maxval(abs(contracted - grad)) <= 1.0e-11_wp*max(1.0_wp, maxval(abs(grad))), &
            & trim(names(im))//": Jacobian disagrees with the VJP")
         if (allocated(error)) return
         call check(error, maxval(abs(sum(grad, dim=2))) <= 1.0e-11_wp, trim(names(im))//": net force")
         if (allocated(error)) return
         center = sum(xyz, dim=2)/real(nat, wp)
         do k = 1, 3
            omega = 0.0_wp
            omega(k) = 1.0_wp
            rigid = 0.0_wp
            do ia = 1, nat
               rigid = rigid + dot_product(grad(:, ia), cross_product(omega, xyz(:, ia) - center))
            end do
            ! About the axis of a line no sphere turns: the molecule stays put
            if (o1%group_detected == "D0h" .or. o1%group_detected == "C0v") then
               omega = omega - dot_product(omega, line_axis(xyz))*line_axis(xyz)
               call check(error, abs(rigid - dot_product(omega, sum(torque, dim=2))) <= 1.0e-9_wp, &
                  & trim(names(im))//": a rigid rotation of the line does not turn its spheres")
            else
               call check(error, abs(rigid - sum(torque(k, :))) <= 1.0e-9_wp, &
                  & trim(names(im))//": a rigid rotation does not turn every sphere with the molecule")
            end if
            if (allocated(error)) return
         end do

         ! A displacement that keeps the whole detected group (one that
         ! keeps a subgroup used alone would lower the symmetry): the
         ! average of its images; a line has the operations used.
         ! Not where a subgroup is used: the molecule cannot tell the
         ! conjugate choices apart (the two angle classes of the sixfold
         ! harmonic of C6v), so rounding picks one image of the spheres,
         ! and the function summed here is no invariant of the group
         if (.not. (used(im) == detected(im) .or. detected(im) == "D0h" .or. detected(im) == "C0v")) then
            deallocate (torque, grad, contracted, jac, mode, sym_mode)
            cycle
         end if
         if (o1%group_detected == "D0h" .or. o1%group_detected == "C0v") then
            full_ops = o1%ops
            full_perm = o1%perm
         else
            call molecule_symmetry(mol, full_ops, full_perm, center)
         end if
         sym_mode = 0.0_wp
         do ig = 1, size(full_ops, 3)
            do ia = 1, nat
               sym_mode(:, full_perm(ia, ig)) = sym_mode(:, full_perm(ia, ig)) &
                  & + matmul(full_ops(:, :, ig), mode(:, ia))
            end do
         end do
         sym_mode = sym_mode/real(size(full_ops, 3), wp)
         do ia = 1, nat
            torque(:, ia) = oh_torque(o1%rotation(:, :, o1%rot_of_atom(ia)), w)
         end do
         call sym_orient_vjp(o1, xyz, torque, grad)
         call orient_value(xyz + STEP*sym_mode, z, "auto", w, ep, error)
         if (allocated(error)) return
         call orient_value(xyz - STEP*sym_mode, z, "auto", w, em, error)
         if (allocated(error)) return
         write (text, '(es14.6," vs FD ",es14.6)') sum(grad*sym_mode), (ep - em)/(2.0_wp*STEP)
         call check(error, sum(grad*sym_mode), (ep - em)/(2.0_wp*STEP), thr=2.0e-7_wp, &
            & message=trim(names(im))//": VJP along a displacement that keeps the group: "//text)
         if (allocated(error)) return
         deallocate (torque, grad, contracted, jac, mode, sym_mode)
      end do
   end subroutine test_point_groups

   !> Unit vector along a linear molecule
   !>
   !> @param[in] xyz Coordinates (3, nat)
   pure function line_axis(xyz) result(axis)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Axis
      real(wp) :: axis(3)

      axis = xyz(:, size(xyz, 2)) - xyz(:, 1)
      axis = axis/norm2(axis)
   end function line_axis

   !> Flags of the in-plane harmonics
   !>
   !> Planar rings about a centre; a regular 13-ring has no harmonic up to
   !> the largest order used
   !>
   !> - a lone harmonic well above its threshold, of a high order: no flag
   !>   (two 10-rings of slightly different size, turned by 18 degrees)
   !> - a lone weak harmonic (13-ring and a small triangle): flagged, and no
   !>   longer once the threshold moves away in either direction
   !> - two harmonics of nearly equal preference (triangle and square): flagged
   !>
   !> @param[out] error Error handle
   subroutine test_frame_harmonic_flags(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_frame_type) :: frame
      type(sym_threshold_type) :: thr
      real(wp), allocatable :: xyz(:, :)
      real(wp) :: rho

      ! rel(10) = (r1^10 - r2^10)/(r1^10 + r2^10) = 0.02
      xyz = rings([10, 10], [2.0_wp*(1.02_wp/0.98_wp)**0.1_wp, 2.0_wp], [0.0_wp, pi/10.0_wp])
      call new_sym_frame(frame, xyz, thr)
      call check(error, frame%level == frame_level_plane .and. frame%harmonic == 10, &
         & "ten-fold rings: tenth harmonic not used")
      if (allocated(error)) return
      call check(error, .not. frame%near_harmonic, "a lone harmonic far above its threshold is flagged")
      if (allocated(error)) return

      ! rel(3) = 3 rho^3/(13 R^3 + 3 rho^3) = 5e-4
      rho = 2.0_wp*(5.0e-4_wp*13.0_wp/(3.0_wp*(1.0_wp - 5.0e-4_wp)))**(1.0_wp/3.0_wp)
      xyz = rings([13, 3], [2.0_wp, rho], [0.0_wp, 0.3_wp])
      call new_sym_frame(frame, xyz, thr)
      call check(error, frame%harmonic == 3 .and. frame%near_harmonic .and. &
         & abs(frame%near_harmonic_value - 5.0e-4_wp) < 1.0e-6_wp, "a weak harmonic is not flagged")
      if (allocated(error)) return
      thr%harmonic = 1.0e-7_wp
      call new_sym_frame(frame, xyz, thr)
      call check(error, frame%harmonic == 3 .and. .not. frame%near_harmonic, &
         & "a harmonic far above a lowered threshold is flagged")
      if (allocated(error)) return
      thr%harmonic = 2.0e-2_wp
      call new_sym_frame(frame, xyz, thr)
      call check(error, frame%harmonic == 0 .and. .not. frame%near_harmonic, &
         & "a harmonic far below a raised threshold is used or flagged")
      if (allocated(error)) return

      ! Preferences rel(3) and 0.9 rel(4) within 1e-3 of each other
      thr = sym_threshold_type()
      xyz = rings([13, 3, 4], [2.0_wp, 0.9_wp, 1.0285_wp], [0.0_wp, 0.3_wp, 0.5_wp])
      call new_sym_frame(frame, xyz, thr)
      call check(error, frame%near_harmonic .and. frame%near_harmonic_value < 5.0e-3_wp, &
         & "two harmonics of nearly equal preference are not flagged")
   end subroutine test_frame_harmonic_flags

   !> Concentric planar regular rings in the plane z = 0
   !>
   !> @param[in] n      Atoms per ring
   !> @param[in] radius Ring radii, bohr
   !> @param[in] phase  Angle of the first atom of every ring, rad
   pure function rings(n, radius, phase) result(xyz)
      !> Atoms per ring
      integer, intent(in) :: n(:)
      !> Radii
      real(wp), intent(in) :: radius(:)
      !> Phases
      real(wp), intent(in) :: phase(:)
      !> Coordinates (3, sum(n))
      real(wp), allocatable :: xyz(:, :)

      real(wp) :: phi
      integer :: ir, k, m

      allocate (xyz(3, sum(n)))
      m = 0
      do ir = 1, size(n)
         do k = 1, n(ir)
            m = m + 1
            phi = phase(ir) + 2.0_wp*pi*real(k - 1, wp)/real(n(ir), wp)
            xyz(:, m) = [radius(ir)*cos(phi), radius(ir)*sin(phi), 0.0_wp]
         end do
      end do
   end function rings

   !> In auto mode the same distortion only lowers the group
   !>
   !> - the detected group keeps its deviation, the reason names it with the
   !>   tolerance
   !> - every atom off: no operation but the identity passes, C1 with the
   !>   tolerance reason; exact geometries detect a deviation below it
   !>
   !> @param[out] error Error handle
   subroutine test_auto_inexact(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      integer :: ia

      call molecule("water", xyz, z)
      xyz(2, 2) = xyz(2, 2) + 1.0e-5_wp
      call new_sym_orient(orient, xyz, z, orient_molecular, "auto", TOL, merr)
      call check(error, .not. allocated(merr), "auto mode must not fail on a distortion")
      if (allocated(error)) return
      call check(error, orient%group_used == "Cs", "distorted water: expected Cs, got "// &
         & orient%group_used)
      if (allocated(error)) return
      call check(error, orient%deviation <= TOL, "group used exceeds the tolerance")
      if (allocated(error)) return
      call check(error, orient%group_detected == "C2v" .and. orient%deviation_detected > TOL .and. &
         & orient%deviation_detected < 1.0e-4_wp, "distorted water: detected deviation")
      if (allocated(error)) return
      call check(error, index(orient%reason, "detected C2v, using Cs (C2v matches only to ") == 1 .and. &
         & index(orient%reason, "bohr, above the tolerance 1.00E-08 bohr)") > 0, &
         & "distorted water: reason '"//orient%reason//"'")
      if (allocated(error)) return

      ! Every atom off: C1, the tolerance named in the reason
      call molecule("benzene", xyz, z)
      do ia = 1, size(xyz, 2)
         xyz(:, ia) = xyz(:, ia) + 2.0e-5_wp*sin(1.7_wp*real(ia, wp) + [0.3_wp, 1.9_wp, 4.1_wp])
      end do
      call new_sym_orient(orient, xyz, z, orient_molecular, "auto", TOL, merr)
      call check(error, .not. allocated(merr), "auto mode must not fail on noise")
      if (allocated(error)) return
      call check(error, orient%group_detected == "D6h" .and. orient%group_used == "C1" .and. &
         & orient%nrot == 1 .and. orient%deviation == 0.0_wp .and. orient%deviation_detected > TOL, &
         & "noisy benzene: expected D6h detected, C1 used, got "//orient%group_detected//", "// &
         & orient%group_used)
      if (allocated(error)) return
      call check(error, index(orient%reason, "geometry matches D6h only to ") == 1 .and. &
         & index(orient%reason, " bohr, above the tolerance 1.00E-08 bohr; no symmetry used") > 0, &
         & "noisy benzene: reason '"//orient%reason//"'")
      if (allocated(error)) return

      ! Exact: the detected deviation is roundoff
      call molecule("benzene", xyz, z)
      call new_sym_orient(orient, xyz, z, orient_molecular, "auto", TOL, merr)
      call check(error, .not. allocated(merr) .and. orient%deviation_detected <= TOL .and. &
         & abs(orient%deviation_detected - orient%deviation) <= TOL, "exact benzene: detected deviation")
   end subroutine test_auto_inexact

   !> A tolerance wider than an atom pair: the reason names the true cause
   !>
   !> - two hydrogens 2.1e-4 bohr apart, unevenly about the mirror plane the
   !>   other atoms fix: within 1e-3 bohr both reflect onto the same atom, so
   !>   the plane is no permutation although no atom is off by the tolerance
   !>
   !> @param[out] error Error handle
   subroutine test_auto_loose(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: loose = 1.0e-3_wp
      integer, parameter :: z(7) = [8, 6, 1, 1, 7, 9, 9]
      real(wp), parameter :: xyz(3, 7) = reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 1.5_wp, 2.0e-4_wp, 1.0_wp, &
         & 1.5_wp, -1.0e-5_wp, 1.0_wp, 2.5_wp, 0.0_wp, 0.3_wp, 2.0_wp, 2.0_wp, -0.5_wp, &
         & 2.0_wp, -2.0_wp, -0.5_wp], [3, 7])
      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr

      call new_sym_orient(orient, xyz, z, orient_molecular, "auto", loose, merr)
      call check(error, .not. allocated(merr), "auto mode must not fail on a loose tolerance")
      if (allocated(error)) return
      call check(error, orient%group_detected == "Cs" .and. orient%group_used == "C1" .and. &
         & orient%deviation_detected <= loose, "close pair: detected "//orient%group_detected// &
         & ", used "//orient%group_used)
      if (allocated(error)) return
      call check(error, index(orient%reason, "the tolerance 1.00E-03 bohr is too loose for Cs") == 1 .and. &
         & index(orient%reason, "do not map the atoms one to one; no symmetry used") > 0 .and. &
         & index(orient%reason, "above the tolerance") == 0, "close pair: reason '"//orient%reason//"'")
   end subroutine test_auto_loose

   !> Exactly symmetric clusters of about a thousand atoms keep their group
   !>
   !> - libmsym places the operations of such clusters only to about 1e-4
   !>   bohr (seen for these seeds); the rigid fit brings them to roundoff
   !> - an atom moved by 1e-6 bohr still loses the full group at the default
   !>   tolerance and keeps it at a loosened one
   subroutine test_refine_large(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=3), parameter :: groups(2) = ["C2v", "D2h"]
      integer, parameter :: nasym(2) = [300, 120]
      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      integer :: ig

      do ig = 1, size(groups)
         call symmetric_cluster(groups(ig), nasym(ig), 1, xyz, z)
         call new_sym_orient(orient, xyz, z, orient_molecular, "auto", TOL, merr)
         call check(error, .not. allocated(merr), groups(ig)//" cluster: detection failed")
         if (allocated(error)) return
         call check(error, orient%group_detected == groups(ig) .and. orient%group_used == groups(ig), &
            & groups(ig)//" cluster: detected "//orient%group_detected//", used "//orient%group_used)
         if (allocated(error)) return
         call check(error, orient%deviation_detected <= 1.0e-11_wp .and. orient%deviation <= 1.0e-11_wp, &
            & groups(ig)//" cluster: deviation above roundoff")
         if (allocated(error)) return
         call check(error, len(orient%reason) == 0, groups(ig)//" cluster: reason '"//orient%reason//"'")
         if (allocated(error)) return

         ! One atom off by 1e-6 bohr
         xyz(1, 1) = xyz(1, 1) + 1.0e-6_wp
         call new_sym_orient(orient, xyz, z, orient_molecular, "auto", TOL, merr)
         call check(error, .not. allocated(merr), groups(ig)//" cluster, distorted: detection failed")
         if (allocated(error)) return
         call check(error, orient%group_detected == groups(ig) .and. orient%group_used /= groups(ig) .and. &
            & orient%deviation_detected > TOL .and. orient%deviation_detected < 2.0e-6_wp .and. &
            & orient%deviation <= TOL, groups(ig)//" cluster, distorted: used "//orient%group_used)
         if (allocated(error)) return
         call new_sym_orient(orient, xyz, z, orient_molecular, "auto", 1.0e-5_wp, merr)
         call check(error, .not. allocated(merr), groups(ig)//" cluster, loose: detection failed")
         if (allocated(error)) return
         call check(error, orient%group_used == groups(ig) .and. orient%deviation < 2.0e-6_wp, &
            & groups(ig)//" cluster, loose: used "//orient%group_used)
         if (allocated(error)) return
      end do
   end subroutine test_refine_large

   !> Large exactly symmetric clusters are detected in full
   !>
   !> - at libmsym's default equivalence threshold unrelated orbits of these
   !>   clusters merge: detection fails for the random ones and reads the
   !>   lattice cut as Th; the threshold tied to the tolerance finds the group
   !> - the cluster of 1000 atoms also runs the second detection at the
   !>   default threshold, whose failure must not replace the group found
   !> - a loosened tolerance returns to the default threshold; the failure is
   !>   retried at tighter ones
   subroutine test_detect_large(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=3), parameter :: groups(3) = ["C2v", "D2h", "D2h"]
      integer, parameter :: nasym(3) = [250, 125, 300], seeds(3) = [2, 15, 2]
      type(sym_orient_type) :: orient
      type(msym_group_type) :: group
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)
      integer :: ic

      do ic = 1, size(groups)
         call symmetric_cluster(groups(ic), nasym(ic), seeds(ic), xyz, z)
         ! The premise: one detection at the default threshold fails
         call find_point_group(xyz, z, group, merr)
         call check(error, allocated(merr), groups(ic)//" cluster: the default threshold no longer fails")
         if (allocated(error)) return
         deallocate (merr)
         call new_sym_orient(orient, xyz, z, orient_molecular, "auto", TOL, merr)
         call check(error, .not. allocated(merr), groups(ic)//" cluster: detection failed")
         if (allocated(error)) return
         call check(error, orient%group_detected == groups(ic) .and. orient%group_used == groups(ic) .and. &
            & orient%deviation <= 1.0e-11_wp .and. len(orient%reason) == 0, groups(ic)// &
            & " cluster: detected "//orient%group_detected//", used "//orient%group_used//" "//orient%reason)
         if (allocated(error)) return
      end do

      ! Loosened tolerance on the last cluster: the retry at tighter thresholds
      call new_sym_orient(orient, xyz, z, orient_molecular, "D2h", 1.0e-4_wp, merr)
      call check(error, .not. allocated(merr), "D2h cluster, loose tolerance: detection failed")
      if (allocated(error)) return
      call check(error, orient%group_used == "D2h" .and. orient%deviation <= 1.0e-4_wp, &
         & "D2h cluster, loose tolerance: used "//orient%group_used)
      if (allocated(error)) return

      ! Sphere cut from an fcc lattice: Oh, read as Th at the default threshold
      call lattice_cluster(9.5_wp, xyz, z)
      call find_point_group(xyz, z, group, merr)
      call check(error, .not. allocated(merr), "lattice cluster: default detection failed")
      if (allocated(error)) return
      call check(error, group%name /= "Oh", "lattice cluster: the default threshold no longer misses Oh")
      if (allocated(error)) return
      call find_point_group(xyz, z, group, merr, TOL)
      call check(error, .not. allocated(merr), "lattice cluster: detection failed")
      if (allocated(error)) return
      call check(error, group%name == "Oh" .and. size(group%ops, 3) == 48, &
         & "lattice cluster: detected "//group%name)
   end subroutine test_detect_large

   !> Atoms of an fcc lattice within a sphere about one of them
   !>
   !> @param[in]  radius Sphere radius in half lattice constants
   !> @param[out] xyz    Coordinates (3, nat), bohr
   !> @param[out] z      Species ids (nat), all one
   subroutine lattice_cluster(radius, xyz, z)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Coordinates
      real(wp), allocatable, intent(out) :: xyz(:, :)
      !> Species ids
      integer, allocatable, intent(out) :: z(:)

      !> Half the lattice constant, bohr
      real(wp), parameter :: half = 1.9_wp
      integer :: i, j, k, m, n, pass

      m = int(radius) + 1
      do pass = 1, 2
         n = 0
         do i = -m, m
            do j = -m, m
               do k = -m, m
                  if (mod(abs(i + j + k), 2) /= 0) cycle
                  if (real(i*i + j*j + k*k, wp) > radius**2) cycle
                  n = n + 1
                  if (pass == 2) xyz(:, n) = half*real([i, j, k], wp)
               end do
            end do
         end do
         if (pass == 1) then
            allocate (xyz(3, n))
            allocate (z(n), source=1)
         end if
      end do
   end subroutine lattice_cluster

   !> Cluster of random points and their images under C2v or D2h
   !>
   !> Points in general position from a linear congruential generator; the
   !> orbits alternate between two species
   !>
   !> @param[in]  group Schoenflies symbol, "C2v" or "D2h"
   !> @param[in]  nasym Number of orbits
   !> @param[in]  seed  Generator seed
   !> @param[out] xyz   Coordinates (3, nat), bohr
   !> @param[out] z     Species ids (nat)
   subroutine symmetric_cluster(group, nasym, seed, xyz, z)
      !> Schoenflies symbol
      character(len=*), intent(in) :: group
      !> Number of orbits
      integer, intent(in) :: nasym
      !> Generator seed
      integer, intent(in) :: seed
      !> Coordinates
      real(wp), allocatable, intent(out) :: xyz(:, :)
      !> Species ids
      integer, allocatable, intent(out) :: z(:)

      real(wp) :: sign(3, 8), p(3), scale
      integer(i8) :: state
      integer :: nop, i, k

      sign(:, 1) = [1.0_wp, 1.0_wp, 1.0_wp]
      sign(:, 2) = [-1.0_wp, -1.0_wp, 1.0_wp]
      sign(:, 3) = [-1.0_wp, 1.0_wp, 1.0_wp]
      sign(:, 4) = [1.0_wp, -1.0_wp, 1.0_wp]
      sign(:, 5:8) = -sign(:, 1:4)
      nop = merge(8, 4, group == "D2h")
      allocate (xyz(3, nop*nasym), z(nop*nasym))
      scale = 2.2_wp*real(nop*nasym, wp)**(1.0_wp/3.0_wp)
      state = int(seed, i8)*7919_i8 + 12345_i8
      do i = 1, nasym
         do k = 1, 3
            state = mod(state*1103515245_i8 + 12345_i8, 2147483648_i8)
            p(k) = scale*(0.1_wp + 0.9_wp*real(state, wp)/2147483648.0_wp)
         end do
         do k = 1, nop
            xyz(:, (i - 1)*nop + k) = sign(:, k)*p
            z((i - 1)*nop + k) = 1 + mod(i, 2)
         end do
      end do
   end subroutine symmetric_cluster

   !* ================================================================================= *!
   !*                                    Derivatives                                    *!
   !* ================================================================================= *!

   !> Eigen-level frame torque against central differences
   subroutine test_frame_vjp_eigen(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call check_frame_vjp(error, "c1")
   end subroutine test_frame_vjp_eigen

   !> Plane-level frame torque against central differences
   !>
   !> Ammonia distorted by 1e-4 bohr: the in-plane eigenvalues stay within
   !> the degeneracy gap, the third harmonic fixes the angle, and the
   !> distortion separates its three candidates (at exact C3v they tie, and
   !> differences across the tie are meaningless)
   subroutine test_frame_vjp_plane(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_frame_type) :: frame
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)

      call molecule("ammonia_distorted", xyz, z)
      call default_frame(frame, xyz)
      call check(error, frame%level, frame_level_plane)
      if (allocated(error)) return
      call check(error, frame%harmonic, 3)
      if (allocated(error)) return
      call check_frame_vjp(error, "ammonia_distorted", fd_step=1.0e-7_wp)
   end subroutine test_frame_vjp_plane

   !> Symmetric orientations: torque through the family angles of Cs sites
   !>
   !> Planar HOCl in Cs: every site sphere may turn about the plane normal,
   !> and in-plane displacements (which keep Cs) turn the frame within the
   !> plane, so the derivative is not zero by symmetry
   subroutine test_orient_vjp_symmetric(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :), torque(:, :), grad(:, :), mode(:, :)
      integer, allocatable :: z(:)
      real(wp) :: w(3, 3), ep, em, analytic, numeric
      integer :: imode, ia

      w = weight_matrix()
      call molecule("hocl", xyz, z)
      allocate (torque(3, 3), grad(3, 3))
      do imode = 1, 2
         if (imode == 1) then
            mode = reshape([0.1_wp, 0.2_wp, 0.0_wp, 0.3_wp, -0.7_wp, 0.0_wp, -0.4_wp, 0.5_wp, 0.0_wp], [3, 3])
         else
            mode = reshape([-0.6_wp, 0.1_wp, 0.0_wp, 0.2_wp, 0.4_wp, 0.0_wp, 0.3_wp, -0.2_wp, 0.0_wp], [3, 3])
         end if

         call new_sym_orient(orient, xyz, z, orient_molecular, "Cs", TOL, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, any(orient%dw /= 0.0_wp), "Cs sites carry no frame derivative")
         if (allocated(error)) return
         do ia = 1, 3
            torque(:, ia) = oh_torque(orient%rotation(:, :, orient%rot_of_atom(ia)), w)
         end do
         call sym_orient_vjp(orient, xyz, torque, grad)
         analytic = sum(grad*mode)

         call orient_value(xyz + STEP*mode, z, "Cs", w, ep, error)
         if (allocated(error)) return
         call orient_value(xyz - STEP*mode, z, "Cs", w, em, error)
         if (allocated(error)) return
         numeric = (ep - em)/(2.0_wp*STEP)
         call check(error, analytic, numeric, thr=1.0e-7_wp)
         if (allocated(error)) return
      end do
   end subroutine test_orient_vjp_symmetric

   !> Forward-mode orientation Jacobian against the reverse mode
   !>
   !> Contracting jac(:, :, :, rot_of_atom(a)) with the torque on every atom a
   !> must reproduce [[sym_orient_vjp]]: frame levels eigen and plane, family
   !> site groups (Cs, C2), a cluster whose group (Oh) keeps no rotation of the
   !> frame; lab orientation all zero
   !>
   !> - rotational sum rule: a rigid rotation turns every sphere with the
   !>   molecule, whatever the group, so the gradient along Omega x d is
   !>   Omega . (total torque), also for torques the group does not keep
   !>
   !> @param[out] error Error handle
   subroutine test_orient_jacobian(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=17), parameter :: names(5) = [character(len=17) :: "c1", "ammonia_distorted", &
         & "hocl", "water", "oh48"]
      character(len=4), parameter :: groups(5) = [character(len=4) :: "off", "off", "Cs", "C2", "auto"]
      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :), torque(:, :), grad(:, :), contracted(:, :), jac(:, :, :, :)
      integer, allocatable :: z(:)
      real(wp) :: center(3), omega(3), rigid
      integer :: im, ia, ib, k

      do im = 1, size(names)
         call molecule(trim(names(im)), xyz, z)
         call new_sym_orient(orient, xyz, z, orient_molecular, trim(groups(im)), TOL, merr)
         if (allocated(merr)) then
            call test_failed(error, trim(names(im))//": "//merr%message)
            return
         end if
         allocate (torque(3, size(z)), grad(3, size(z)), contracted(3, size(z)), &
            & jac(3, 3, size(z), orient%nrot))
         do ia = 1, size(z)
            torque(:, ia) = sin(1.3_wp*real(ia, wp) + [0.2_wp, 1.9_wp, 3.7_wp])
         end do
         call sym_orient_vjp(orient, xyz, torque, grad)
         call sym_orient_jacobian(orient, xyz, jac)
         contracted = 0.0_wp
         do ia = 1, size(z)
            do ib = 1, size(z)
               do k = 1, 3
                  contracted(k, ib) = contracted(k, ib) &
                     & + dot_product(torque(:, ia), jac(:, k, ib, orient%rot_of_atom(ia)))
               end do
            end do
         end do
         center = sum(xyz, dim=2)/real(size(z), wp)
         do k = 1, 3
            omega = 0.0_wp
            omega(k) = 1.0_wp
            rigid = 0.0_wp
            do ia = 1, size(z)
               rigid = rigid + dot_product(grad(:, ia), cross_product(omega, xyz(:, ia) - center))
            end do
            call check(error, abs(rigid - sum(torque(k, :))) <= 1.0e-10_wp, &
               & trim(names(im))//": a rigid rotation does not turn every sphere with the molecule")
            if (allocated(error)) return
         end do
         call check(error, maxval(abs(sum(grad, dim=2))) <= 1.0e-12_wp, &
            & trim(names(im))//": torque gradient has a net force")
         if (allocated(error)) return
         call check(error, maxval(abs(contracted - grad)) <= 1.0e-12_wp*max(1.0_wp, maxval(abs(grad))), &
            & trim(names(im))//": Jacobian disagrees with the VJP")
         if (allocated(error)) return
         deallocate (torque, grad, contracted, jac)
      end do

      ! Lab orientation: no derivative
      call molecule("c1", xyz, z)
      call new_sym_orient(orient, xyz, z, orient_lab, "off", TOL, merr)
      allocate (jac(3, 3, size(z), orient%nrot), source=1.0_wp)
      call sym_orient_jacobian(orient, xyz, jac)
      call check(error, all(jac == 0.0_wp), "lab orientation has a Jacobian")
   end subroutine test_orient_jacobian

   !> Orientation derivative where the in-plane angle of the frame follows the lab
   !>
   !> The `axial` frame does not turn about its axis with the molecule: lines
   !> with a subgroup (HCN Cs, CO2 C2 and C2h), whose operations the frame
   !> places, and a triatomic bent by 0.005 bohr (C2v) and a regular 13-ring
   !> (D13h, and its subgroup C2v with spheres off the sites), whose
   !> operations turn with the atoms
   !>
   !> - angular velocity of every sphere along rigid rotations: the Jacobian
   !>   against central differences
   !> - the VJP is the contracted Jacobian
   !>
   !> @param[out] error Error handle
   subroutine test_orient_lab_angle(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=6), parameter :: names(6) = [character(len=6) :: "hcn", "co2", "co2", "bent", "ring13", "ring13"]
      character(len=4), parameter :: groups(6) = [character(len=4) :: "Cs", "C2", "C2h", "auto", "auto", "C2v"]
      character(len=4), parameter :: used(6) = [character(len=4) :: "Cs", "C2", "C2h", "C2v", "D13h", "C2v"]
      type(sym_orient_type) :: o0, shifted(2)
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :), d(:, :), moved(:, :), jac(:, :, :, :), torque(:, :), grad(:, :)
      integer, allocatable :: z(:)
      real(wp) :: center(3), axis(3), numeric(3), analytic(3), theta, contracted
      character(len=80) :: text
      integer :: im, iaxis, ia, ib, k, nat, is

      do im = 1, size(names)
         select case (trim(names(im)))
         case ("bent")
            z = [8, 6, 8]
            xyz = reshape([-2.2_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.005_wp, 0.0_wp, 2.2_wp, 0.0_wp, 0.0_wp], [3, 3])
         case ("ring13")
            z = [(6, k=1, 13)]
            xyz = reshape([(4.0_wp*cos(2.0_wp*pi*real(k, wp)/13.0_wp), 4.0_wp*sin(2.0_wp*pi*real(k, wp)/13.0_wp), &
               & 0.0_wp, k=1, 13)], [3, 13])
         case default
            call molecule(trim(names(im)), xyz, z)
         end select
         nat = size(z)
         xyz = transformed(xyz, rotation_matrix(0.5_wp, 1.2_wp, -0.4_wp), [0.4_wp, -1.1_wp, 2.2_wp])
         center = sum(xyz, dim=2)/real(nat, wp)
         d = xyz - spread(center, 2, nat)
         text = trim(names(im))//" ("//trim(groups(im))//")"
         call new_sym_orient(o0, xyz, z, orient_molecular, trim(groups(im)), TOL, merr)
         if (allocated(merr)) then
            call test_failed(error, trim(text)//": "//merr%message)
            return
         end if
         call check(error, o0%group_used == trim(used(im)) .and. o0%frame%level == frame_level_axial, &
            & trim(text)//": group "//o0%group_used//" or frame level not as expected")
         if (allocated(error)) return
         allocate (jac(3, 3, nat, o0%nrot), torque(3, nat), grad(3, nat), moved(3, nat))
         call sym_orient_jacobian(o0, xyz, jac)

         do iaxis = 1, 3
            axis = sin(1.7_wp*real(iaxis, wp) + [0.3_wp, 1.5_wp, 2.9_wp])
            axis = axis/norm2(axis)
            do is = 1, 2
               theta = merge(STEP, -STEP, is == 1)
               do ia = 1, nat
                  moved(:, ia) = center + cos(theta)*d(:, ia) + sin(theta)*cross_product(axis, d(:, ia)) &
                     & + (1.0_wp - cos(theta))*dot_product(axis, d(:, ia))*axis
               end do
               call new_sym_orient(shifted(is), moved, z, orient_molecular, trim(groups(im)), TOL, merr)
               if (allocated(merr)) then
                  call test_failed(error, trim(text)//": "//merr%message)
                  return
               end if
            end do
            do ia = 1, nat
               numeric = turn_rate(o0%rotation(:, :, o0%rot_of_atom(ia)), &
                  & shifted(1)%rotation(:, :, shifted(1)%rot_of_atom(ia)), &
                  & shifted(2)%rotation(:, :, shifted(2)%rot_of_atom(ia)))
               analytic = 0.0_wp
               do ib = 1, nat
                  analytic = analytic + matmul(jac(:, :, ib, o0%rot_of_atom(ia)), cross_product(axis, d(:, ib)))
               end do
               call check(error, maxval(abs(analytic - numeric)) <= 1.0e-7_wp, &
                  & trim(text)//": a sphere does not turn as its Jacobian along a rigid rotation")
               if (allocated(error)) return
            end do
         end do

         do ia = 1, nat
            torque(:, ia) = sin(1.3_wp*real(ia, wp) + [0.2_wp, 1.9_wp, 3.7_wp])
         end do
         call sym_orient_vjp(o0, xyz, torque, grad)
         do ib = 1, nat
            do k = 1, 3
               contracted = 0.0_wp
               do ia = 1, nat
                  contracted = contracted + dot_product(torque(:, ia), jac(:, k, ib, o0%rot_of_atom(ia)))
               end do
               call check(error, abs(contracted - grad(k, ib)) <= 1.0e-9_wp*max(1.0_wp, maxval(abs(grad))), &
                  & trim(text)//": Jacobian disagrees with the VJP")
               if (allocated(error)) return
            end do
         end do
         deallocate (jac, torque, grad, moved)
      end do

   contains

      !> Left angular velocity of an orientation by central differences
      !>
      !> The shifted orientations are compared modulo Oh
      !>
      !> @param[in] r0 Orientation
      !> @param[in] rp Orientation a step ahead
      !> @param[in] rm Orientation a step behind
      pure function turn_rate(r0, rp, rm) result(omega)
         !> Orientation
         real(wp), intent(in) :: r0(3, 3)
         !> Orientation a step ahead
         real(wp), intent(in) :: rp(3, 3)
         !> Orientation a step behind
         real(wp), intent(in) :: rm(3, 3)
         !> Angular velocity
         real(wp) :: omega(3)

         real(wp) :: a(3, 3)

         a = matmul(oh_nearest(r0, rp) - oh_nearest(r0, rm), transpose(r0))/(2.0_wp*STEP)
         omega = 0.5_wp*[a(3, 2) - a(2, 3), a(1, 3) - a(3, 1), a(2, 1) - a(1, 2)]
      end function turn_rate

      !> Image r o, o in Oh, of an orientation nearest to a reference
      !>
      !> @param[in] r0 Reference
      !> @param[in] r  Orientation
      pure function oh_nearest(r0, r) result(image)
         !> Reference
         real(wp), intent(in) :: r0(3, 3)
         !> Orientation
         real(wp), intent(in) :: r(3, 3)
         !> Image
         real(wp) :: image(3, 3)

         real(wp) :: m(3, 3), o(3, 3)
         integer :: i, j

         m = matmul(transpose(r0), r)
         o = 0.0_wp
         do j = 1, 3
            i = maxloc(abs(m(:, j)), dim=1)
            o(i, j) = sign(1.0_wp, m(i, j))
         end do
         image = matmul(r, transpose(o))
      end function oh_nearest
   end subroutine test_orient_lab_angle

   !> Frame torque of an Oh-invariant function against central differences
   !>
   !> @param[out] error   Error handle
   !> @param[in]  name    Molecule
   !> @param[in]  fd_step Central-difference step, default `STEP`
   subroutine check_frame_vjp(error, name, fd_step)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Molecule
      character(len=*), intent(in) :: name
      !> Central-difference step
      real(wp), intent(in), optional :: fd_step

      type(sym_frame_type) :: frame
      real(wp), allocatable :: xyz(:, :), grad(:, :), disp(:, :)
      integer, allocatable :: z(:)
      real(wp) :: w(3, 3), ep, em, h, numeric
      integer :: ia, k

      h = STEP
      if (present(fd_step)) h = fd_step

      w = weight_matrix()
      call molecule(name, xyz, z)
      call default_frame(frame, xyz)
      allocate (grad(3, size(xyz, 2)))
      call sym_frame_vjp(frame, xyz, oh_torque(frame%axes, w), grad)

      do ia = 1, size(xyz, 2)
         do k = 1, 3
            disp = xyz
            disp(k, ia) = disp(k, ia) + h
            call default_frame(frame, disp)
            ep = oh_value(frame%axes, w)
            disp(k, ia) = disp(k, ia) - 2.0_wp*h
            call default_frame(frame, disp)
            em = oh_value(frame%axes, w)
            numeric = (ep - em)/(2.0_wp*h)
            call check(error, abs(grad(k, ia) - numeric) <= 1.0e-6_wp*max(1.0_wp, abs(numeric)), &
               & name//": frame torque disagrees with central differences")
            if (allocated(error)) return
         end do
      end do
   end subroutine check_frame_vjp

   !> Sum over atoms of the Oh-invariant function of their sphere orientation
   !>
   !> @param[in]  xyz   Coordinates
   !> @param[in]  z     Species
   !> @param[in]  group Requested group
   !> @param[in]  w     Weight matrix
   !> @param[out] value Function value
   !> @param[out] error Error handle
   subroutine orient_value(xyz, z, group, w, value, error)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species
      integer, intent(in) :: z(:)
      !> Requested group
      character(len=*), intent(in) :: group
      !> Weight matrix
      real(wp), intent(in) :: w(3, 3)
      !> Function value
      real(wp), intent(out) :: value
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_orient_type) :: orient
      type(moist_error_type), allocatable :: merr
      integer :: ia

      value = 0.0_wp
      call new_sym_orient(orient, xyz, z, orient_molecular, group, TOL, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      do ia = 1, size(xyz, 2)
         value = value + oh_value(orient%rotation(:, :, orient%rot_of_atom(ia)), w)
      end do
   end subroutine orient_value

   !> f(R) = sum_k (r_k^T W r_k)^2 over the columns r_k: invariant under R -> R o, o in Oh
   !>
   !> @param[in] r Orientation
   !> @param[in] w Symmetric weight matrix
   pure function oh_value(r, w) result(f)
      !> Orientation
      real(wp), intent(in) :: r(3, 3)
      !> Weight matrix
      real(wp), intent(in) :: w(3, 3)
      !> Value
      real(wp) :: f

      integer :: k

      f = 0.0_wp
      do k = 1, 3
         f = f + dot_product(r(:, k), matmul(w, r(:, k)))**2
      end do
   end function oh_value

   !> Torque of `oh_value`: delta f = tau . delta omega for delta R = [delta omega]_x R
   !>
   !> @param[in] r Orientation
   !> @param[in] w Symmetric weight matrix
   pure function oh_torque(r, w) result(tau)
      !> Orientation
      real(wp), intent(in) :: r(3, 3)
      !> Weight matrix
      real(wp), intent(in) :: w(3, 3)
      !> Torque
      real(wp) :: tau(3)

      real(wp) :: s, wr(3)
      integer :: k

      tau = 0.0_wp
      do k = 1, 3
         wr = matmul(w, r(:, k))
         s = dot_product(r(:, k), wr)
         tau = tau + 4.0_wp*s*cross_product(r(:, k), wr)
      end do
   end function oh_torque

   !> A fixed symmetric weight matrix without special axes
   pure function weight_matrix() result(w)
      !> Weight matrix
      real(wp) :: w(3, 3)

      w = reshape([1.3_wp, 0.2_wp, -0.4_wp, 0.2_wp, 0.7_wp, 0.5_wp, -0.4_wp, 0.5_wp, -0.9_wp], [3, 3])
   end function weight_matrix

   !* ================================================================================= *!
   !*                                      Helpers                                      *!
   !* ================================================================================= *!

   !> Detected group of a reference molecule
   !>
   !> @param[out] error Error handle
   !> @param[in]  name  Molecule
   !> @param[in]  want  Expected Schoenflies symbol
   subroutine expect_detected(error, name, want)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Molecule
      character(len=*), intent(in) :: name
      !> Expected symbol
      character(len=*), intent(in) :: want

      type(msym_group_type) :: group
      type(moist_error_type), allocatable :: merr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)

      call molecule(name, xyz, z)
      call find_point_group(xyz, z, group, merr)
      if (allocated(merr)) then
         call test_failed(error, name//": "//merr%message)
         return
      end if
      call check(error, group%name == want, name//": expected "//want//", got "//group%name)
   end subroutine expect_detected

   !> Group used for a reference molecule
   !>
   !> @param[out] error   Error handle
   !> @param[in]  name    Molecule
   !> @param[in]  request Requested group
   !> @param[in]  want    Expected group used
   !> @param[in]  reduced Whether a reason must be reported
   subroutine expect_used(error, name, request, want, reduced)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Molecule
      character(len=*), intent(in) :: name
      !> Requested group
      character(len=*), intent(in) :: request
      !> Expected group used
      character(len=*), intent(in) :: want
      !> Whether a reason must be reported
      logical, intent(in) :: reduced

      type(sym_orient_type) :: orient
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)

      call orient_of(error, name, request, xyz, z, orient)
      if (allocated(error)) return
      call check(error, orient%group_used == want, name//" ("//request//"): expected "//want// &
         & ", got "//orient%group_used)
      if (allocated(error)) return
      call check(error, (len(orient%reason) > 0) .eqv. reduced, name//": unexpected reason '"// &
         & orient%reason//"'")
   end subroutine expect_used

   !> Orientations of a reference molecule
   !>
   !> @param[out] error   Error handle
   !> @param[in]  name    Molecule
   !> @param[in]  request Requested group
   !> @param[out] xyz     Coordinates
   !> @param[out] z       Species
   !> @param[out] orient  Orientations
   subroutine orient_of(error, name, request, xyz, z, orient)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Molecule
      character(len=*), intent(in) :: name
      !> Requested group
      character(len=*), intent(in) :: request
      !> Coordinates
      real(wp), allocatable, intent(out) :: xyz(:, :)
      !> Species
      integer, allocatable, intent(out) :: z(:)
      !> Orientations
      type(sym_orient_type), intent(out) :: orient

      type(moist_error_type), allocatable :: merr

      call molecule(trim(name), xyz, z)
      call new_sym_orient(orient, xyz, z, orient_molecular, request, TOL, merr)
      if (allocated(merr)) call test_failed(error, trim(name)//": "//merr%message)
   end subroutine orient_of

   !> The reported eigen-gap of a line is zero wherever the line points
   !>
   !> Off the lab axes the two small eigenvalues are rounding of the
   !> eigensolver, not zero
   !>
   !> @param[out] error Error handle
   subroutine test_frame_gap_line(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_frame_type) :: frame
      real(wp) :: xyz(3, 3), dir(3)
      integer :: k

      do k = 1, 3
         select case (k)
         case (1)
            dir = [0.0_wp, 0.0_wp, 1.0_wp]
         case (2)
            dir = [1.0_wp, 1.0_wp, 1.0_wp]/sqrt(3.0_wp)
         case default
            dir = [0.3_wp, -0.8_wp, 0.52_wp]/norm2([0.3_wp, -0.8_wp, 0.52_wp])
         end select
         xyz(:, 1) = [0.4_wp, -1.1_wp, 2.2_wp]
         xyz(:, 2) = xyz(:, 1) + 2.2_wp*dir
         xyz(:, 3) = xyz(:, 1) - 2.2_wp*dir
         call default_frame(frame, xyz)
         call check(error, frame%level == frame_level_axial .and. frame%gap == 0.0_wp .and. &
            & .not. frame%near_linear, "frame of a line reports a gap or a near switch")
         if (allocated(error)) return
         ! Two atoms
         call default_frame(frame, xyz(:, 1:2))
         call check(error, frame%level == frame_level_axial .and. frame%gap == 0.0_wp, &
            & "frame of two atoms reports a gap")
         if (allocated(error)) return
      end do
   end subroutine test_frame_gap_line

   !> One frame for tops, rods and near-isotropic molecules
   !>
   !> - level, harmonic and no flag for exact tops and a thick rod
   !> - the frame follows a rotation, shift and renumbering modulo Oh, also
   !>   where the candidates of a harmonic tie by symmetry
   !> - the frame torque against central differences for the cases without
   !>   an exact tie
   !>
   !> @param[out] error Error handle
   subroutine test_frame_single_tops(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=16), parameter :: names(9) = [character(len=16) :: "ammonia", "benzene", "methane", &
         & "rod", "bipyramid", "bipyramid_off", "methane_off", "ammonia_off", "c1_round"]
      integer, parameter :: levels(9) = [frame_level_plane, frame_level_plane, frame_level_sphere, &
         & frame_level_eigen, frame_level_plane, frame_level_eigen, frame_level_eigen, frame_level_eigen, &
         & frame_level_eigen]
      !> Whether the geometry has no exact tie, so differences are meaningful
      logical, parameter :: generic(9) = [.false., .false., .false., .true., .false., .true., .true., &
         & .true., .true.]
      type(sym_frame_type) :: f1, f2
      real(wp), allocatable :: xyz(:, :), moved(:, :), grad(:, :), x(:, :)
      real(wp) :: q(3, 3), w(3, 3), ep, em, numeric
      integer :: im, nat, ia, k

      w = weight_matrix()
      q = rotation_matrix(1.1_wp, -0.7_wp, 2.3_wp)
      do im = 1, size(names)
         call frame_case(names(im), xyz)
         nat = size(xyz, 2)
         call default_frame(f1, xyz)
         call check(error, f1%level, levels(im), trim(names(im))//": unexpected frame level")
         if (allocated(error)) return
         if (im <= 5) then
            call check(error, .not. (f1%near_gap .or. f1%near_linear .or. f1%near_harmonic), &
               & trim(names(im))//": flagged near a switch")
            if (allocated(error)) return
         end if

         moved = transformed(xyz, q, [3.0_wp, -1.0_wp, 0.5_wp])
         call default_frame(f2, moved(:, nat:1:-1))
         if (generic(im)) then
            call check(error, is_oh_element(matmul(transpose(matmul(q, f1%axes)), f2%axes)), &
               & trim(names(im))//": frame does not follow the molecule")
         else
            ! A tie by symmetry: rounding picks one of the equivalent frames
            call check(error, f2%level == f1%level .and. f2%harmonic == f1%harmonic, &
               & trim(names(im))//": frame level changed")
         end if
         if (allocated(error)) return
         if (.not. generic(im)) cycle

         allocate (grad(3, nat))
         call sym_frame_vjp(f1, xyz, oh_torque(f1%axes, w), grad)
         do ia = 1, nat
            do k = 1, 3
               x = xyz
               x(k, ia) = x(k, ia) + 1.0e-6_wp
               call default_frame(f2, x)
               ep = oh_value(f2%axes, w)
               x(k, ia) = x(k, ia) - 2.0e-6_wp
               call default_frame(f2, x)
               em = oh_value(f2%axes, w)
               numeric = (ep - em)/2.0e-6_wp
               call check(error, abs(grad(k, ia) - numeric) <= 2.0e-5_wp*max(1.0_wp, abs(numeric)), &
                  & trim(names(im))//": frame torque disagrees with central differences")
               if (allocated(error)) return
            end do
         end do
         deallocate (grad)
      end do
   end subroutine test_frame_single_tops

   !> Exactly tied frames are listed, and rounding does not pick among them
   !>
   !> - benzene: the two classes of the sixfold harmonic tie, one tied
   !>   frame; a regular pentagon: five classes, four tied frames; water: none
   !> - benzene at a generic pose, turned by multiples of 1e-4 rad about the
   !>   lab axes: the frame turns along, modulo Oh (a choice by rounding
   !>   would switch the class between neighbouring poses)
   !> - atoms in pairs through a centre at the sphere level: the in-plane
   !>   harmonics tie two frames a sixth of a turn apart, the harmonic odd in
   !>   the height does not; the frame follows a rotation, modulo Oh
   !>
   !> @param[out] error Error handle
   subroutine test_frame_tied_stable(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_frame_type) :: frame, frame0
      type(sym_threshold_type) :: thr
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :)
      real(wp) :: ring(3, 5), q(3, 3), q0(3, 3), t
      integer :: k, it, axis

      call structure_of("benzene", mol)
      call new_sym_frame(frame, mol%xyz, thr)
      call check(error, allocated(frame%tied), "benzene: no tied frame listed")
      if (allocated(error)) return
      call check(error, size(frame%tied, 3), 1, "benzene: number of tied frames")
      if (allocated(error)) return
      do k = 1, 5
         ring(:, k) = 2.4_wp*[cos(2.0_wp*pi*real(k - 1, wp)/5.0_wp), sin(2.0_wp*pi*real(k - 1, wp)/5.0_wp), 0.0_wp]
      end do
      call new_sym_frame(frame, ring, thr)
      call check(error, allocated(frame%tied), "pentagon: no tied frame listed")
      if (allocated(error)) return
      call check(error, size(frame%tied, 3), 4, "pentagon: number of tied frames")
      if (allocated(error)) return
      call structure_of("water", mol)
      call new_sym_frame(frame, mol%xyz, thr)
      call check(error, .not. allocated(frame%tied), "water: tied frames listed")
      if (allocated(error)) return

      call structure_of("benzene", mol)
      q0 = axis_turn(3, 0.7_wp)
      q0 = matmul(q0, matmul(axis_turn(1, 1.3_wp), axis_turn(3, 0.4_wp)))
      call new_sym_frame(frame0, matmul(q0, mol%xyz), thr)
      do axis = 1, 3
         do it = -3, 3
            t = 1.0e-4_wp*real(it, wp)
            q = axis_turn(axis, t)
            call new_sym_frame(frame, matmul(q, matmul(q0, mol%xyz)), thr)
            call check(error, is_oh_element(matmul(transpose(matmul(q, frame0%axes)), frame%axes)), &
               & "benzene: the tied frame switches under a small rotation")
            if (allocated(error)) return
         end do
      end do

      call structure_of("ci", mol)
      thr%gap = 0.99_wp
      xyz = mol%xyz
      call new_sym_frame(frame0, xyz, thr)
      call check(error, frame0%level == frame_level_sphere .and. frame0%axis_order == 4, &
         & "ci: not a sphere-level frame from the quartic moment")
      if (allocated(error)) return
      call check(error, .not. allocated(frame0%tied), "ci: tied frames listed")
      if (allocated(error)) return
      do it = 1, 6
         q = matmul(axis_turn(3, 0.9_wp*real(it, wp)), matmul(axis_turn(1, 1.7_wp*real(it, wp)), &
            & axis_turn(2, 0.3_wp + real(it, wp))))
         call new_sym_frame(frame, matmul(q, xyz), thr)
         call check(error, is_oh_element(matmul(transpose(matmul(q, frame0%axes)), frame%axes)), &
            & "ci: the frame does not follow a rotation")
         if (allocated(error)) return
      end do

   contains

      !> Rotation about a lab axis
      !>
      !> @param[in] k     Lab axis
      !> @param[in] angle Angle, rad
      pure function axis_turn(k, angle) result(r)
         !> Lab axis
         integer, intent(in) :: k
         !> Angle
         real(wp), intent(in) :: angle
         !> Rotation
         real(wp) :: r(3, 3)

         integer :: i, j

         i = mod(k, 3) + 1
         j = mod(k + 1, 3) + 1
         r = 0.0_wp
         r(k, k) = 1.0_wp
         r(i, i) = cos(angle)
         r(j, j) = cos(angle)
         r(j, i) = sin(angle)
         r(i, j) = -sin(angle)
      end function axis_turn

   end subroutine test_frame_tied_stable

   !> Orbits that tie under a symmetry outside the group follow one choice
   !>
   !> Benzene in C6h, D6, D3h, C3h, C6, C2v and C2h: the spheres of the
   !> carbon and of the hydrogen ring each tie between orientations that a
   !> symmetry of the molecule outside the group relates
   !>
   !> - a rotated and renumbered copy gets the spheres of the reference
   !>   under one operation of D6h, all twelve of them
   !>
   !> @param[out] error Error handle
   subroutine test_orbit_ties_follow(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=3), parameter :: groups(7) = [character(len=3) :: "C6h", "D6", "D3h", "C3h", "C6", "C2v", "C2h"]
      integer, parameter :: n = 12
      type(sym_orient_type) :: full, ref, copy
      type(structure_type) :: mol
      type(moist_error_type), allocatable :: err
      real(wp) :: q(3, 3), x(3, n), a, b, c
      integer :: z(n), p(n), ig, it, iop, ia, ja, nok, best, k
      character(len=80) :: label

      call structure_of("benzene", mol)
      z = mol%num(mol%id)
      call new_sym_orient(full, mol%xyz, z, orient_molecular, "D6h", 1.0e-8_wp, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      do ig = 1, size(groups)
         call new_sym_orient(ref, mol%xyz, z, orient_molecular, trim(groups(ig)), 1.0e-8_wp, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         do it = 1, 8
            a = 1.7_wp*real(it, wp)
            b = 0.9_wp*real(it, wp) + 0.3_wp
            c = 2.3_wp*real(it, wp)
            q = matmul(reshape([cos(a), sin(a), 0.0_wp, -sin(a), cos(a), 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3]), &
               & matmul(reshape([cos(b), 0.0_wp, -sin(b), 0.0_wp, 1.0_wp, 0.0_wp, sin(b), 0.0_wp, cos(b)], [3, 3]), &
               & reshape([cos(c), sin(c), 0.0_wp, -sin(c), cos(c), 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])))
            p = cshift([(k, k=1, n)], 5*it)
            if (mod(it, 2) == 0) p = p(n:1:-1)
            if (mod(it, 3) == 0) p([1, 7]) = p([7, 1])
            x = matmul(q, mol%xyz(:, p))
            call new_sym_orient(copy, x, z(p), orient_molecular, trim(groups(ig)), 1.0e-8_wp, err)
            if (allocated(err)) then
               call test_failed(error, err%message)
               return
            end if
            best = 0
            do iop = 1, size(full%ops, 3)
               nok = 0
               do ia = 1, n
                  ! Moved atom ia is original p(ia), the image of atom ja
                  ja = findloc(full%perm(:, iop), p(ia), dim=1)
                  if (is_oh_element(matmul(transpose(copy%rotation(:, :, copy%rot_of_atom(ia))), &
                     & matmul(q, matmul(full%ops(:, :, iop), ref%rotation(:, :, ref%rot_of_atom(ja))))))) nok = nok + 1
               end do
               best = max(best, nok)
            end do
            write (label, '("benzene in ",a,", motion ",i0,": ",i0," of 12 spheres follow")') &
               & trim(groups(ig)), it, best
            call check(error, best, n, trim(label))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_orbit_ties_follow

   !> Flags of a frame near a switch of its choices
   !>
   !> - eigen-gap: two second moments 1e-3 of their sum apart
   !> - line: CO2 bent by 0.01 rad, mean square distance about 1e-4 bohr^2
   !> - not for water and CO2 bent by 0.3 rad
   !> - moving a threshold away clears the flag
   !>
   !> @param[out] error Error handle
   subroutine test_frame_near_switch(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(sym_frame_type) :: frame
      type(sym_threshold_type) :: thr
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: z(:)

      call frame_case("near_top", xyz)
      call default_frame(frame, xyz)
      call check(error, frame%near_gap .and. abs(frame%near_gap_value - 1.0e-3_wp) < 1.0e-4_wp, &
         & "near-symmetric top: eigen-gap not flagged")
      if (allocated(error)) return
      thr%gap = 1.0e-6_wp
      call new_sym_frame(frame, xyz, thr)
      call check(error, .not. frame%near_gap .and. frame%level == frame_level_eigen, &
         & "near-symmetric top: flagged far above the threshold")
      if (allocated(error)) return

      call frame_case("co2_bent_small", xyz)
      call default_frame(frame, xyz)
      call check(error, frame%near_linear .and. frame%tau > 1.0e-5_wp .and. frame%tau < 1.0e-2_wp, &
         & "slightly bent CO2: line not flagged")
      if (allocated(error)) return
      call frame_case("co2_bent", xyz)
      call default_frame(frame, xyz)
      call check(error, .not. frame%near_linear, "bent CO2: line flagged")
      if (allocated(error)) return

      call frame_case("methane_off", xyz)
      call default_frame(frame, xyz)
      call check(error, frame%near_gap, "distorted tetrahedron: eigen-gap not flagged")
      if (allocated(error)) return

      call molecule("water", xyz, z)
      call default_frame(frame, xyz)
      call check(error, .not. (frame%near_gap .or. frame%near_linear .or. frame%near_harmonic), &
         & "water: flagged near a switch")
   end subroutine test_frame_near_switch

   !> Geometries for the single-frame tests, in a generic position
   !>
   !> @param[in]  name Case
   !> @param[out] xyz  Coordinates (3, nat)
   subroutine frame_case(name, xyz)
      !> Case
      character(len=*), intent(in) :: name
      !> Coordinates
      real(wp), allocatable, intent(out) :: xyz(:, :)

      integer, allocatable :: z(:)
      type(sym_frame_type) :: f
      real(wp), allocatable :: d(:, :), loc(:, :)
      real(wp) :: t
      integer :: nat, k

      select case (trim(name))
      case ("ammonia", "benzene", "methane")
         call molecule(trim(name), xyz, z)
      case ("ammonia_off")
         call molecule("ammonia", xyz, z)
         xyz(:, 2) = xyz(:, 2) + [0.0_wp, 0.2_wp, 0.1_wp]
      case ("rod")
         ! Thick helix 60 bohr long with an elliptic cross-section: about
         ! 1e-3 of its second moment is transverse, yet it is no line
         nat = 40
         allocate (xyz(3, nat))
         do k = 1, nat
            t = 1.7_wp*real(k, wp)
            xyz(:, k) = [0.45_wp*cos(t), 0.45_wp*sin(t), 1.5_wp*real(k, wp)]
            xyz(1, k) = 1.4_wp*xyz(1, k)
            xyz(2, k) = xyz(2, k) + 0.02_wp*sin(0.3_wp*real(k*k, wp))
         end do
      case ("co2_bent_small", "co2_bent")
         t = merge(0.01_wp, 0.3_wp, trim(name) == "co2_bent_small")
         allocate (xyz(3, 3))
         xyz(:, 1) = 0.0_wp
         xyz(:, 2) = 2.2_wp*[sin(t), 0.0_wp, cos(t)]
         xyz(:, 3) = 2.2_wp*[sin(t), 0.0_wp, -cos(t)]
      case ("bipyramid", "bipyramid_off")
         ! Trigonal bipyramid with an isotropic second moment: the fourth
         ! moment gives the axis; one atom moved for `bipyramid_off`
         allocate (xyz(3, 6), source=0.0_wp)
         do k = 1, 3
            t = 2.0_wp*pi*real(k, wp)/3.0_wp
            xyz(:, k + 1) = 2.0_wp*[cos(t), sin(t), 0.0_wp]
         end do
         xyz(3, 5) = 2.0_wp*sqrt(0.75_wp)
         xyz(3, 6) = -xyz(3, 5)
         if (trim(name) == "bipyramid_off") xyz(:, 2) = xyz(:, 2) + [0.1_wp, 0.2_wp, -0.15_wp]
      case ("methane_off")
         ! Distorted tetrahedron
         call molecule("methane", xyz, z)
         xyz(:, 2) = xyz(:, 2) + [0.04_wp, -0.02_wp, 0.03_wp]
         xyz(:, 4) = xyz(:, 4) + [0.01_wp, 0.03_wp, -0.02_wp]
      case ("near_top", "c1_round")
         ! No symmetry, the two smallest second moments 1e-3 (`near_top`)
         ! or 1e-6 of their sum apart
         t = merge(1.0e-3_wp, 1.0e-6_wp, trim(name) == "near_top")
         call molecule("c1", xyz, z)
         nat = size(xyz, 2)
         call default_frame(f, xyz)
         d = xyz - spread(sum(xyz, dim=2)/real(nat, wp), 2, nat)
         loc = matmul(transpose(f%evec), d)
         loc(2, :) = loc(2, :)*sqrt(f%lambda(1)*(1.0_wp + t)/((1.0_wp - t)*f%lambda(2)))
         xyz = matmul(f%evec, loc)
      case default
         error stop "frame_case: unknown case"
      end select
      xyz = transformed(xyz, rotation_matrix(0.5_wp, 1.2_wp, -0.4_wp), [0.4_wp, -1.1_wp, 2.2_wp])
   end subroutine frame_case

   !> Molecular frame under the default thresholds
   !>
   !> @param[out] frame Frame
   !> @param[in]  xyz   Coordinates (3, nat)
   subroutine default_frame(frame, xyz)
      !> Frame
      type(sym_frame_type), intent(out) :: frame
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)

      call new_sym_frame(frame, xyz, sym_threshold_type())
   end subroutine default_frame

   !> @param[out] xyz  Coordinates (3, nat)
   !> @param[out] z    Atomic numbers
   recursive subroutine molecule(name, xyz, z)
      !> Molecule
      character(len=*), intent(in) :: name
      !> Coordinates
      real(wp), allocatable, intent(out) :: xyz(:, :)
      !> Atomic numbers
      integer, allocatable, intent(out) :: z(:)

      real(wp) :: phi, c, s, a, gr
      integer :: k

      select case (trim(name))
      case ("water")
         z = [8, 1, 1]
         xyz = reshape([0.0_wp, 0.0_wp, 0.1173_wp, 0.0_wp, 0.7572_wp, -0.4692_wp, &
            & 0.0_wp, -0.7572_wp, -0.4692_wp], [3, 3])*aatoau
      case ("ammonia")
         z = [7, 1, 1, 1]
         allocate (xyz(3, 4))
         xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k + 1) = [0.9377_wp*cos(phi), 0.9377_wp*sin(phi), -0.3816_wp]
         end do
         xyz = xyz*aatoau
      case ("ammonia_distorted")
         call molecule("ammonia", xyz, z)
         xyz(:, 2) = xyz(:, 2) + [1.0e-4_wp, 0.4e-4_wp, 0.0_wp]
      case ("methane")
         z = [6, 1, 1, 1, 1]
         a = 0.629_wp*aatoau
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, a, a, a, a, -a, -a, -a, a, -a, -a, -a, a], [3, 5])
      case ("benzene")
         allocate (xyz(3, 12), z(12))
         do k = 1, 6
            phi = 2.0_wp*pi*real(k - 1, wp)/6.0_wp
            c = cos(phi)
            s = sin(phi)
            xyz(:, k) = [1.397_wp*c, 1.397_wp*s, 0.0_wp]*aatoau
            xyz(:, k + 6) = [2.481_wp*c, 2.481_wp*s, 0.0_wp]*aatoau
            z(k) = 6
            z(k + 6) = 1
         end do
      case ("bf3")
         z = [5, 9, 9, 9]
         allocate (xyz(3, 4))
         xyz(:, 1) = 0.0_wp
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k + 1) = [1.30_wp*cos(phi), 1.30_wp*sin(phi), 0.0_wp]*aatoau
         end do
      case ("sf6")
         z = [16, 9, 9, 9, 9, 9, 9]
         a = 1.56_wp*aatoau
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, a, 0.0_wp, 0.0_wp, -a, 0.0_wp, 0.0_wp, &
            & 0.0_wp, a, 0.0_wp, 0.0_wp, -a, 0.0_wp, 0.0_wp, 0.0_wp, a, 0.0_wp, 0.0_wp, -a], [3, 7])
      case ("icosahedron")
         gr = (1.0_wp + sqrt(5.0_wp))/2.0_wp
         z = [(6, k=1, 12)]
         xyz = reshape([0.0_wp, 1.0_wp, gr, 0.0_wp, 1.0_wp, -gr, 0.0_wp, -1.0_wp, gr, &
            & 0.0_wp, -1.0_wp, -gr, 1.0_wp, gr, 0.0_wp, 1.0_wp, -gr, 0.0_wp, -1.0_wp, gr, 0.0_wp, &
            & -1.0_wp, -gr, 0.0_wp, gr, 0.0_wp, 1.0_wp, gr, 0.0_wp, -1.0_wp, -gr, 0.0_wp, 1.0_wp, &
            & -gr, 0.0_wp, -1.0_wp], [3, 12])*1.4_wp*aatoau
      case ("co2")
         z = [6, 8, 8]
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.16_wp, 0.0_wp, 0.0_wp, -1.16_wp], &
            & [3, 3])*aatoau
      case ("hocl")
         z = [8, 1, 17]
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.966_wp, 0.0_wp, 0.0_wp, -0.383_wp, 1.643_wp, 0.0_wp], &
            & [3, 3])*aatoau
      case ("oh48")
         ! Oh cluster: six on the C4 axes twice, twelve on the C2 axes and 24
         ! in the diagonal mirror planes
         xyz = reshape([oh_orbit([2.0_wp, 0.0_wp, 0.0_wp]), oh_orbit([3.5_wp, 0.0_wp, 0.0_wp]), &
            & oh_orbit([2.2_wp, 2.2_wp, 0.0_wp]), oh_orbit([1.5_wp, 1.5_wp, 3.1_wp])], [3, 48])
         z = [(16, k=1, 6), (9, k=1, 6), (1, k=1, 12), (6, k=1, 24)]
      case ("acetylene")
         z = [1, 6, 6, 1]
         xyz = reshape([0.0_wp, 0.0_wp, -1.664_wp, 0.0_wp, 0.0_wp, -0.601_wp, 0.0_wp, 0.0_wp, 0.601_wp, &
            & 0.0_wp, 0.0_wp, 1.664_wp], [3, 4])*aatoau
      case ("hcn")
         z = [1, 6, 7]
         xyz = reshape([0.0_wp, 0.0_wp, -1.066_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.156_wp], &
            & [3, 3])*aatoau
      case ("ch3f")
         ! C3v: carbon and fluorine on the axis
         z = [6, 9, 1, 1, 1]
         allocate (xyz(3, 5), source=0.0_wp)
         xyz(:, 2) = [0.0_wp, 0.0_wp, 1.383_wp]
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k + 2) = [1.029_wp*cos(phi), 1.029_wp*sin(phi), -0.350_wp]
         end do
         xyz = xyz*aatoau
      case ("ethane")
         ! Staggered, D3d
         z = [6, 6, 1, 1, 1, 1, 1, 1]
         allocate (xyz(3, 8), source=0.0_wp)
         xyz(:, 1) = [0.0_wp, 0.0_wp, 0.7655_wp]
         xyz(:, 2) = [0.0_wp, 0.0_wp, -0.7655_wp]
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k + 2) = [1.020_wp*cos(phi), 1.020_wp*sin(phi), 1.1606_wp]
            phi = phi + pi/3.0_wp
            xyz(:, k + 5) = [1.020_wp*cos(phi), 1.020_wp*sin(phi), -1.1606_wp]
         end do
         xyz = xyz*aatoau
      case ("allene")
         ! D2d: the two CH2 planes perpendicular
         z = [6, 6, 6, 1, 1, 1, 1]
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.308_wp, 0.0_wp, 0.0_wp, -1.308_wp, &
            & 0.9327_wp, 0.0_wp, 1.866_wp, -0.9327_wp, 0.0_wp, 1.866_wp, &
            & 0.0_wp, 0.9327_wp, -1.866_wp, 0.0_wp, -0.9327_wp, -1.866_wp], [3, 7])*aatoau
      case ("cyclopropane")
         ! D3h without an atom on the axis
         z = [6, 6, 6, 1, 1, 1, 1, 1, 1]
         allocate (xyz(3, 9))
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k) = [0.8718_wp*cos(phi), 0.8718_wp*sin(phi), 0.0_wp]
            xyz(:, 2*k + 2) = [1.4553_wp*cos(phi), 1.4553_wp*sin(phi), 0.9124_wp]
            xyz(:, 2*k + 3) = [1.4553_wp*cos(phi), 1.4553_wp*sin(phi), -0.9124_wp]
         end do
         xyz = xyz*aatoau
      case ("pf5")
         ! D3h trigonal bipyramid
         z = [15, 9, 9, 9, 9, 9]
         allocate (xyz(3, 6), source=0.0_wp)
         xyz(:, 2) = [0.0_wp, 0.0_wp, 1.577_wp]
         xyz(:, 3) = [0.0_wp, 0.0_wp, -1.577_wp]
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k + 3) = [1.534_wp*cos(phi), 1.534_wp*sin(phi), 0.0_wp]
         end do
         xyz = xyz*aatoau
      case ("xef4")
         ! D4h square plane
         z = [54, 9, 9, 9, 9]
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.95_wp, 0.0_wp, 0.0_wp, -1.95_wp, 0.0_wp, 0.0_wp, &
            & 0.0_wp, 1.95_wp, 0.0_wp, 0.0_wp, -1.95_wp, 0.0_wp], [3, 5])*aatoau
      case ("benzene_li")
         ! C6v: benzene with an atom on the axis; its site C6v leaves two C3v
         ! subgroups that the molecular frame cannot tell apart
         call molecule("benzene", xyz, z)
         xyz = reshape([xyz, [0.0_wp, 0.0_wp, 1.9_wp*aatoau]], [3, 13])
         z = [z, 3]
      case ("benzene_lif")
         ! C6v with two different atoms on the axis
         call molecule("benzene", xyz, z)
         xyz = reshape([xyz, [0.0_wp, 0.0_wp, 1.9_wp*aatoau], [0.0_wp, 0.0_wp, -2.6_wp*aatoau]], [3, 14])
         z = [z, 3, 9]
      case ("sandwich_x2")
         ! D6h with axial ligands: the metal (site D3d of the D3d used) and the
         ! ligands (site C3v) both tie under C2 about the axis
         call molecule("sandwich", xyz, z)
         xyz = reshape([xyz, [0.0_wp, 0.0_wp, 3.9_wp*aatoau], [0.0_wp, 0.0_wp, -3.9_wp*aatoau]], [3, 15])
         z = [z, 17, 17]
      case ("sandwich")
         ! D6h: an atom between two eclipsed six-rings; its site D6h leaves two
         ! D3d subgroups that the molecular frame cannot tell apart
         z = [24, (6, k=1, 12)]
         allocate (xyz(3, 13), source=0.0_wp)
         do k = 1, 6
            phi = 2.0_wp*pi*real(k - 1, wp)/6.0_wp
            xyz(:, k + 1) = [1.42_wp*cos(phi), 1.42_wp*sin(phi), 1.61_wp]
            xyz(:, k + 7) = [1.42_wp*cos(phi), 1.42_wp*sin(phi), -1.61_wp]
         end do
         xyz = xyz*aatoau
      case default
         ! Bromochlorofluoromethane-like, no symmetry
         z = [6, 1, 9, 17, 35]
         xyz = reshape([0.02_wp, -0.01_wp, 0.03_wp, 0.65_wp, 0.71_wp, 0.48_wp, -0.81_wp, 0.52_wp, &
            & -0.69_wp, 0.93_wp, -1.24_wp, -0.61_wp, -1.12_wp, -0.88_wp, 1.37_wp], [3, 5])*aatoau
      end select
   end subroutine molecule

   !> Distinct signed permutations of a point: its orbit under Oh
   !>
   !> @param[in] p Point
   pure function oh_orbit(p) result(pts)
      !> Point
      real(wp), intent(in) :: p(3)
      !> Orbit (3, n), in a fixed order
      real(wp), allocatable :: pts(:, :)

      integer, parameter :: perms(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, &
         & 3, 2, 1], [3, 6])
      real(wp) :: all_pts(3, 48), cand(3)
      integer :: ip, is, n, k
      logical :: seen

      n = 0
      do ip = 1, 6
         do is = 0, 7
            cand = p(perms(:, ip))*[merge(-1.0_wp, 1.0_wp, btest(is, 0)), &
               & merge(-1.0_wp, 1.0_wp, btest(is, 1)), merge(-1.0_wp, 1.0_wp, btest(is, 2))]
            seen = .false.
            do k = 1, n
               seen = seen .or. all(abs(all_pts(:, k) - cand) < 1.0e-12_wp)
            end do
            if (seen) cycle
            n = n + 1
            all_pts(:, n) = cand
         end do
      end do
      pts = all_pts(:, 1:n)
   end function oh_orbit

   !> Operations of the detected group and the atom permutations they induce
   !>
   !> The identity alone when libmsym finds nothing
   !>
   !> @param[in]  xyz  Coordinates
   !> @param[in]  z    Species
   !> @param[out] ops  Operations (3, 3, order), about the centroid
   !> @param[out] perm Image atoms (nat, order)
   subroutine symmetry_of(xyz, z, ops, perm)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species
      integer, intent(in) :: z(:)
      !> Operations
      real(wp), allocatable, intent(out) :: ops(:, :, :)
      !> Image atoms
      integer, allocatable, intent(out) :: perm(:, :)

      type(msym_group_type) :: group
      type(moist_error_type), allocatable :: merr
      real(wp) :: center(3), target(3)
      integer :: ig, ia, nat

      nat = size(xyz, 2)
      call find_point_group(xyz, z, group, merr)
      if (allocated(merr) .or. group%linear) then
         ops = reshape(identity3, [3, 3, 1])
      else
         ops = group%ops
      end if
      center = sum(xyz, dim=2)/real(nat, wp)
      allocate (perm(nat, size(ops, 3)))
      do ig = 1, size(ops, 3)
         do ia = 1, nat
            target = center + matmul(ops(:, :, ig), xyz(:, ia) - center)
            perm(ia, ig) = minloc(norm2(xyz - spread(target, 2, nat), dim=1), dim=1)
         end do
      end do
   end subroutine symmetry_of

   !> Rotated and translated coordinates
   !>
   !> @param[in] xyz   Coordinates
   !> @param[in] q     Rotation
   !> @param[in] shift Translation
   pure function transformed(xyz, q, shift) result(moved)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Rotation
      real(wp), intent(in) :: q(3, 3)
      !> Translation
      real(wp), intent(in) :: shift(3)
      !> Moved coordinates
      real(wp), allocatable :: moved(:, :)

      moved = matmul(q, xyz) + spread(shift, 2, size(xyz, 2))
   end function transformed

   !> Rotation from z-y-z Euler angles
   !>
   !> @param[in] a First angle
   !> @param[in] b Second angle
   !> @param[in] c Third angle
   pure function rotation_matrix(a, b, c) result(q)
      !> First angle
      real(wp), intent(in) :: a
      !> Second angle
      real(wp), intent(in) :: b
      !> Third angle
      real(wp), intent(in) :: c
      !> Rotation
      real(wp) :: q(3, 3)

      real(wp) :: ra(3, 3), rb(3, 3), rc(3, 3)

      ra = reshape([cos(a), sin(a), 0.0_wp, -sin(a), cos(a), 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      rb = reshape([cos(b), 0.0_wp, -sin(b), 0.0_wp, 1.0_wp, 0.0_wp, sin(b), 0.0_wp, cos(b)], [3, 3])
      rc = reshape([cos(c), sin(c), 0.0_wp, -sin(c), cos(c), 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      q = matmul(ra, matmul(rb, rc))
   end function rotation_matrix

   !> Determinant of a 3x3 matrix
   !>
   !> @param[in] m Matrix
   pure function det(m) result(d)
      !> Matrix
      real(wp), intent(in) :: m(3, 3)
      !> Determinant
      real(wp) :: d

      d = dot_product(m(:, 1), cross_product(m(:, 2), m(:, 3)))
   end function det

end module test_math_symmetry
