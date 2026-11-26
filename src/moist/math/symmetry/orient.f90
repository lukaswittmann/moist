!> Sphere orientation of every atom: lab, molecular frame or point-group symmetric
!>
!> Provides
!>
!> - [[sym_orient_type]]: orientations, orbits, the group used and why
!> - [[new_sym_orient]]: orientations of a set of atoms
!> - [[sym_orient_vjp]], [[sym_orient_jacobian]]: their derivative in reverse
!>   and forward mode
!>
!> | orientation | group        | sphere of atom a                                       |
!> |-------------|--------------|--------------------------------------------------------|
!> | lab         | off          | identity                                               |
!> | molecular   | off          | molecular frame F ([[moist_math_symmetry_frame]])      |
!> | (implied)   | auto, symbol | invariant under the site group of a, closest to F      |
!>
!> - grid convention: the Lebedev node u of atom a sits at c_a + r_a R u,
!>   R = `rotation(:, :, rot_of_atom(a))`
!>
!> With a group (`symmetric_orient`):
!>
!>  1. libmsym detects a candidate group G, its equivalence threshold
!>     following `tolerance` ([[moist_math_symmetry_msym:find_point_group]])
!>  2. every operation is checked against the atoms at `tolerance`, after a
!>     rigid fit of the operations where one misses (`refine_operations`)
!>  3. the group used, H, is the best of G and its subgroups that passes the
!>     check and admits an orientation at every site
!>     ([[moist_math_symmetry_embed]]); largest order first (`candidate_key`)
!>  4. a linear G is replaced by D4h (D_inf_h) or C4v (C_inf_v) about the
!>     molecular axis (`linear_group`)
!>  5. per orbit, one representative takes the admissible orientation
!>     closest to F; atom b = g(rep) gets det(g) g R_rep, an improper g
!>     acting through -g since the grid has inversion symmetry
!>     (`orient_orbits`)
!>
!> Exact ties of these choices are broken by neighbour forms
!> ([[moist_math_symmetry_embed:neighbour_form_type]]), functions of the
!> molecule alone, never by the lab frame or the atom order
!>
!> | choice               | tie (example)                          | broken by            | in               |
!> |----------------------|----------------------------------------|----------------------|------------------|
!> | site orientation     | ammonia N: R and C2(u) R               | form at the site     | `orient_orbits`  |
!> | orbit representative | members whose best orientations tie    | form at one member   | `orient_orbits`  |
!> | group H              | the two C3v of C6v, the two D3d of D6h | sorted orbit forms   | `config_order`   |
!>
!> - `species` ids must not depend on the atom order (atomic numbers, say):
!>   forms sort atoms by them
module moist_math_symmetry_orient
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matdet_3x3
   use mctc_io_utils, only: to_lower
   use moist_utils_timer, only: timer_type
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_lapack_syev, only: lapack_syev
   use moist_math_linalg_geometry, only: identity3, cross_product
   use moist_math_symmetry_msym, only: msym_group_type, find_point_group, msym_op_proper, &
      & msym_type_dnh
   use moist_math_symmetry_frame, only: sym_frame_type, sym_frame_vjp, sym_threshold_type, &
      & new_sym_frame
   use moist_math_symmetry_embed, only: site_candidates, smallest_forms, oh_closeness, is_oh_element, sort_keys, &
      & neighbour_form_type, new_neighbour_form, compare_forms, tie_tol
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: sym_orient_type, new_sym_orient, sym_orient_vjp, sym_orient_jacobian
   public :: orient_lab, orient_molecular, orient_codes

   !> Spheres in the lab frame
   integer, parameter :: orient_lab = 0
   !> Spheres in the molecular frame
   integer, parameter :: orient_molecular = 1
   !> Codes and names of the orientations, for descriptions of the orientation
   character(len=*), parameter :: orient_codes = "0 lab frame, 1 molecular frame"

   !> Tolerance for matching orientations (modulo Oh) and their derivative maps
   real(wp), parameter :: dedupe_tol = 1.0e-12_wp
   !> Length of a `subgroup_name` result
   integer, parameter :: subgroup_name_len = 3
   !> Tolerance of the candidate group keys (integers and a sum of closenesses)
   real(wp), parameter :: key_tol = 1.0e-8_wp
   !> Smallest inertia eigenvalue with a rotation, relative to the largest
   real(wp), parameter :: inertia_tol = 1.0e-12_wp

   !> Per-atom sphere orientations and the data of their derivative
   type :: sym_orient_type
      !> `orient_lab` or `orient_molecular`
      integer :: orientation = orient_lab
      !> Requested group: "off", "auto" or a Schoenflies symbol
      character(len=:), allocatable :: group_requested
      !> Group libmsym detected, empty when not asked or failed
      character(len=:), allocatable :: group_detected
      !> Group whose operations the orientations respect, "C1" without symmetry
      character(len=:), allocatable :: group_used
      !> Why the group used differs from the detected or requested one, or empty
      character(len=:), allocatable :: reason
      !> Largest atom deviation under the operations used
      real(wp) :: deviation = 0.0_wp
      !> Largest atom deviation under the operations of the detected group
      !>
      !> - 0 without a detection
      real(wp) :: deviation_detected = 0.0_wp
      !> Molecular frame
      type(sym_frame_type) :: frame
      !> Whether tied frames differ in species
      !>
      !> - the positions leave the frame open; orientations then depend on the pose
      logical :: ambiguous = .false.
      !> Number of distinct orientations
      integer :: nrot = 0
      !> Distinct orientations modulo Oh (3, 3, nrot)
      !>
      !> - atoms share an entry when R_b = R_a o, o in Oh, and their `dw`
      !>   agree: the Lebedev sphere is Oh-invariant, so both give the same
      !>   points in another order, and a constant o leaves the angular
      !>   velocity unchanged, so the points also move alike
      !> - the entry keeps the orientation of the first atom that needs it
      !> - the exact identity in the lab orientation
      real(wp), allocatable :: rotation(:, :, :)
      !> Derivative maps (3, 3, nrot): delta omega_R = dw delta omega_F
      !>
      !> - delta omega_F the rotation of the frame, delta R = [delta omega_R]_x R
      real(wp), allocatable :: dw(:, :, :)
      !> Orientation index of every atom (nat)
      integer, allocatable :: rot_of_atom(:)
      !> Number of orbits
      !>
      !> - orbit: atoms the operations used map onto each other; without
      !>   symmetry every atom is its own orbit
      integer :: norbit = 0
      !> Orbit index of every atom (nat), orbits numbered by their lowest atom
      integer, allocatable :: orbit_of_atom(:)
      !> Representative of every orbit (norbit), its orientation transported to the rest
      integer, allocatable :: orbit_rep(:)
      !> Operations of the group used (3, 3, order), about `center`
      real(wp), allocatable :: ops(:, :, :)
      !> Image atom of every atom under every operation (nat, order)
      integer, allocatable :: perm(:, :)
      !> Fixed point of the operations
      real(wp) :: center(3) = 0.0_wp
      !> Whether the operations are placed by the frame (a line)
      !>
      !> - they then turn with the frame, not with the atoms
      logical :: ops_on_frame = .false.
      !> Axes of a cubic lattice about the atom centroid, as columns (3, 3)
      real(wp) :: lattice(3, 3) = identity3
   end type sym_orient_type

   !> Candidate group: indices into an operation list
   type :: candidate_type
      !> Schoenflies symbol
      character(len=:), allocatable :: name
      !> Operation indices
      integer, allocatable :: ops(:)
   end type candidate_type

   !> Orientations an orbit member offers: its admissible ones tied closest to F
   type :: offer_type
      !> Closeness of the best one to the frame
      real(wp) :: value = 0.0_wp
      !> Tied orientations (3, 3, n)
      real(wp), allocatable :: rot(:, :, :)
      !> Their derivative maps (3, 3, n)
      real(wp), allocatable :: dw(:, :, :)
   end type offer_type

contains

   !> Sphere orientations of a set of atoms
   !>
   !> - a group other than "off" and "C1" implies the molecular orientation
   !> - "C1" is handled like "off", without libmsym
   !> - a single atom needs no symmetry: any valid symbol is accepted
   !> - timer phases under the open node: "Molecular frame", "Detection"
   !>   (libmsym), "Atom matching" (image atoms and the exactness check),
   !>   "Orientation" (group choice and orbit orientations)
   !>
   !> @param[out] self        Orientations
   !> @param[in]  xyz         Coordinates (3, nat), bohr
   !> @param[in]  species     Species id per atom (nat), positive; atoms of one
   !>                         id must be interchangeable for the cavity, and
   !>                         ids must not depend on the atom order
   !> @param[in]  orientation `orient_lab` or `orient_molecular`
   !> @param[in]  group       "off", "auto" or a Schoenflies symbol
   !> @param[in]  tolerance   Largest atom deviation under an operation, bohr
   !> @param[out] error       Explicit group absent or not exact, bad input
   !> @param[in,out] timer    Records the phases under its open node
   !> @param[in]  threshold   Thresholds of the frame choices; defaults
   subroutine new_sym_orient(self, xyz, species, orientation, group, tolerance, error, timer, threshold)
      !> Orientations
      type(sym_orient_type), intent(out) :: self
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Orientation mode
      integer, intent(in) :: orientation
      !> Requested group
      character(len=*), intent(in) :: group
      !> Exactness tolerance
      real(wp), intent(in) :: tolerance
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Phase timer
      type(timer_type), intent(inout), optional :: timer
      !> Thresholds of the frame choices
      type(sym_threshold_type), intent(in), optional :: threshold

      type(sym_threshold_type) :: thr
      integer :: nat
      logical :: off, auto
      !> Requested group, lower case
      character(len=:), allocatable :: req

      nat = size(xyz, 2)
      if (size(species) /= nat) then
         call fatal_error(error, "Symmetry: species must have one entry per atom")
         return
      end if
      ! C1 has no operation to respect; libmsym may not list it as a subgroup
      self%group_requested = trim(adjustl(group))
      req = to_lower(self%group_requested)
      off = req == "off" .or. req == "c1"
      auto = req == "auto"
      self%group_detected = ""
      self%group_used = "C1"
      self%reason = ""
      self%orientation = orientation
      if (.not. off) self%orientation = orient_molecular

      if (self%orientation == orient_molecular) then
         call time_phase(timer, "Molecular frame", .true.)
         if (present(threshold)) thr = threshold
         call new_sym_frame(self%frame, xyz, thr)
         self%ambiguous = frame_ambiguous(self%frame, xyz, species)
         call time_phase(timer, "Molecular frame", .false.)
      end if
      call set_trivial(self, nat)
      call set_lattice(self, xyz, species)
      if (off) return
      if (nat < 2) then
         ! One atom has every point group and no operation to respect
         if (.not. auto) then
            if (.not. valid_symbol(req)) then
               call fatal_error(error, "Symmetry: unknown point group "//trim(group))
               return
            end if
            self%reason = "requested "//trim(group)//", a single atom needs no symmetry"
         end if
         return
      end if

      call symmetric_orient(self, xyz, species, auto, tolerance, error, timer)
      if (allocated(error)) return
      call set_lattice(self, xyz, species)
   end subroutine new_sym_orient

   !> Axes of a cubic lattice about the atom centroid
   !>
   !> - lab orientation: the identity; otherwise the molecular frame F
   !> - with operations: the lattice is a sphere whose site group is the whole
   !>   group, at the centroid (their fixed point); the admissible orientation
   !>   closest to F ([[site_candidates]]), ties broken by the neighbour form of
   !>   the centroid, so the operations map the lattice onto itself
   !> - F when the group fits no cube (C5, C6, S3 or higher axes)
   !>
   !> @param[in,out] self    Orientations, frame and operations set
   !> @param[in]     xyz     Coordinates (3, nat), bohr
   !> @param[in]     species Species ids (nat)
   subroutine set_lattice(self, xyz, species)
      !> Orientations
      type(sym_orient_type), intent(inout) :: self
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)

      !> Admissible orientations tied closest to the frame, and their maps
      real(wp), allocatable :: rot(:, :, :), dw(:, :, :)
      !> Atoms and their centroid, the origin of the forms
      real(wp), allocatable :: points(:, :)
      !> Whether each tied orientation has the smallest form
      logical, allocatable :: smallest(:)
      real(wp) :: value
      integer :: nat
      logical :: found

      self%lattice = self%frame%axes
      if (self%orientation /= orient_molecular) return
      if (.not. allocated(self%ops)) return
      if (size(self%ops, 3) <= 1) return
      call site_candidates(self%ops, self%frame%axes, found, value, rot, dw)
      if (.not. found .or. size(rot, 3) == 0) return
      if (size(rot, 3) == 1) then
         self%lattice = rot(:, :, 1)
         return
      end if
      nat = size(xyz, 2)
      allocate (points(3, nat + 1))
      points(:, 1:nat) = xyz
      points(:, nat + 1) = sum(xyz, dim=2)/real(nat, wp)
      smallest = smallest_forms(points, [species, 0], nat + 1, rot)
      self%lattice = rot(:, :, findloc(smallest, .true., dim=1))
   end subroutine set_lattice

   !> Whether a frame ties exactly with one the species tell apart from it
   !>
   !> - the frame follows the atom positions alone
   !> - frames that tie exactly are related by a symmetry of the positions,
   !>   or by one the frame rules do not see
   !> - the molecule in two such frames looks the same only when the species
   !>   agree as well: its form about the centroid, modulo Oh
   !>
   !> @param[in] frame   Molecular frame with its tied frames
   !> @param[in] xyz     Coordinates (3, nat)
   !> @param[in] species Species ids (nat)
   pure function frame_ambiguous(frame, xyz, species) result(ambiguous)
      !> Molecular frame
      type(sym_frame_type), intent(in) :: frame
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Whether a tied frame shows another molecule
      logical :: ambiguous

      type(neighbour_form_type) :: chosen, other
      real(wp), allocatable :: points(:, :)
      integer :: nat, k

      ambiguous = .false.
      if (.not. allocated(frame%tied)) return
      nat = size(xyz, 2)
      ! The centroid as one more point, the origin of the forms
      allocate (points(3, nat + 1))
      points(:, 1:nat) = xyz
      points(:, nat + 1) = sum(xyz, dim=2)/real(nat, wp)
      call new_neighbour_form(chosen, points, [species, 0], nat + 1, frame%axes)
      do k = 1, size(frame%tied, 3)
         call new_neighbour_form(other, points, [species, 0], nat + 1, frame%tied(:, :, k))
         ambiguous = compare_forms(chosen, other) /= 0
         if (ambiguous) return
      end do
   end function frame_ambiguous

   !> Start or stop a phase timer, nothing without a timer
   !>
   !> @param[in,out] timer Phase timer
   !> @param[in]     name  Phase name
   !> @param[in]     start Whether the phase starts, else stops
   subroutine time_phase(timer, name, start)
      !> Phase timer
      type(timer_type), intent(inout), optional :: timer
      !> Phase name
      character(len=*), intent(in) :: name
      !> Whether the phase starts
      logical, intent(in) :: start

      if (.not. present(timer)) return
      if (start) then
         call timer%start(name)
      else
         call timer%stop(name)
      end if
   end subroutine time_phase

   !> No symmetry: every atom shares the frame (or the lab frame)
   !>
   !> @param[in,out] self Orientations
   !> @param[in]     nat  Number of atoms
   subroutine set_trivial(self, nat)
      !> Orientations
      type(sym_orient_type), intent(inout) :: self
      !> Number of atoms
      integer, intent(in) :: nat

      integer :: ia

      self%nrot = 1
      if (allocated(self%rotation)) deallocate (self%rotation)
      if (allocated(self%dw)) deallocate (self%dw)
      allocate (self%rotation(3, 3, 1), self%dw(3, 3, 1))
      self%rotation(:, :, 1) = self%frame%axes
      self%dw(:, :, 1) = 0.0_wp
      if (self%orientation == orient_molecular) self%dw(:, :, 1) = identity3
      self%rot_of_atom = [(1, ia=1, nat)]
      self%norbit = nat
      self%orbit_of_atom = [(ia, ia=1, nat)]
      self%orbit_rep = [(ia, ia=1, nat)]
      self%ops = reshape(identity3, [3, 3, 1])
      self%perm = reshape([(ia, ia=1, nat)], [nat, 1])
      self%deviation = 0.0_wp
   end subroutine set_trivial

   !> Symmetric mode: detect, check, choose H and orient every orbit
   !>
   !> - check: every atom lands within `tolerance` of an atom of its species
   !>   under every operation, one to one; the geometry is never symmetrized
   !> - candidates: the detected group and its subgroups, restricted to the
   !>   requested group and its subgroups when a symbol is given
   !> - H: the best candidate (`candidate_key`) that passes the check and is
   !>   embeddable (`orient_orbits`)
   !> - candidates with equal keys (the two C3v of C6v, the two D3d of D6h,
   !>   related by an operation that keeps F modulo Oh): `config_order`
   !> - no candidate left: no symmetry, molecular frame, `reason` says why;
   !>   an error only with a requested group (detection failed, group
   !>   absent or not exact)
   !>
   !> @param[in,out] self      Orientations, molecular frame already set
   !> @param[in]     xyz       Coordinates (3, nat)
   !> @param[in]     species   Species ids (nat)
   !> @param[in]     auto      Whether the group is detected, not requested
   !> @param[in]     tolerance Exactness tolerance
   !> @param[out]    error     Requested group absent or not exact
   !> @param[in,out] timer     Phase timer; a failed call leaves its phase open
   subroutine symmetric_orient(self, xyz, species, auto, tolerance, error, timer)
      !> Orientations
      type(sym_orient_type), intent(inout) :: self
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Detected, not requested
      logical, intent(in) :: auto
      !> Exactness tolerance
      real(wp), intent(in) :: tolerance
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Phase timer
      type(timer_type), intent(inout), optional :: timer

      type(msym_group_type) :: detected
      type(error_type), allocatable :: detect_error
      !> Operation list the candidates index
      real(wp), allocatable :: gops(:, :, :)
      type(candidate_type), allocatable :: cands(:)
      !> Per operation: image atoms, deviation, pass
      integer, allocatable :: perm_op(:, :)
      real(wp), allocatable :: dev_op(:)
      logical, allocatable :: pass_op(:), in_pool(:), tried(:), tie(:)
      real(wp), allocatable :: key(:, :)
      character(len=:), allocatable :: req
      !> Orientations of a tied candidate and of the best one so far
      type(sym_orient_type) :: trial, chosen
      integer :: nat, iop, ic, jc, best
      logical :: any_req, embeddable
      !> Smallest deviation among the candidates of the requested name
      real(wp) :: dev_req
      !> Distance of every atom from the fixed point, the atoms in its order
      !> and the place of every atom in it
      real(wp), allocatable :: radius(:)
      integer, allocatable :: by_radius(:), rank(:)
      integer, parameter :: no_key(1, 1) = 0

      nat = size(xyz, 2)
      req = to_lower(self%group_requested)

      call time_phase(timer, "Detection", .true.)
      call find_point_group(xyz, species, detected, detect_error, tolerance)
      call time_phase(timer, "Detection", .false.)
      if (allocated(detect_error)) then
         if (auto) then
            self%reason = "detection failed ("//detect_error%message//"); no symmetry used"
            call set_trivial(self, nat)
            return
         end if
         call fatal_error(error, "Symmetry: "//detect_error%message)
         return
      end if
      self%group_detected = detected%name
      self%center = detected%center

      call time_phase(timer, "Atom matching", .true.)
      if (detected%linear) then
         call linear_group(self, detected, gops, cands)
      else
         gops = detected%ops
         allocate (cands(size(detected%subgroups) + 1))
         cands(1)%name = detected%name
         cands(1)%ops = [(iop, iop=1, size(gops, 3))]
         do ic = 1, size(detected%subgroups)
            cands(ic + 1)%name = detected%subgroups(ic)%name
            cands(ic + 1)%ops = detected%subgroups(ic)%ops
         end do
      end if

      ! Exactness of every operation at the caller's tolerance; an image keeps
      ! its distance from the fixed point, so the search runs over the atoms
      ! ordered by that distance
      allocate (radius(nat), by_radius(nat), rank(nat))
      do ic = 1, nat
         radius(ic) = norm2(xyz(:, ic) - self%center)
      end do
      call sort_keys(no_key, by_radius, radius)
      rank(by_radius) = [(ic, ic=1, nat)]
      allocate (perm_op(nat, size(gops, 3)), dev_op(size(gops, 3)))
      do iop = 1, size(gops, 3)
         call atom_images(xyz, species, self%center, gops(:, :, iop), radius, by_radius, rank, &
            & perm_op(:, iop), dev_op(iop))
      end do
      ! A loose tolerance can match two atoms to the same image; such an
      ! operation is no permutation and does not count as a symmetry
      allocate (pass_op(size(gops, 3)))
      do iop = 1, size(gops, 3)
         pass_op(iop) = dev_op(iop) <= tolerance .and. is_permutation(perm_op(:, iop))
      end do
      ! libmsym places the operations at its own thresholds; when one misses
      ! the tolerance the whole set is fitted to the atoms as one rigid body
      if (.not. all(pass_op)) then
         call refine_operations(xyz, species, self%center, radius, by_radius, rank, tolerance, &
            & gops, perm_op, dev_op, pass_op)
      end if
      self%deviation_detected = maxval(dev_op)
      call time_phase(timer, "Atom matching", .false.)

      call time_phase(timer, "Orientation", .true.)
      ! Pool of admissible candidates
      allocate (in_pool(size(cands)))
      if (auto) then
         do ic = 1, size(cands)
            in_pool(ic) = all(pass_op(cands(ic)%ops))
         end do
      else
         in_pool = .false.
         any_req = .false.
         dev_req = huge(1.0_wp)
         do ic = 1, size(cands)
            if (.not. same_group_name(cands(ic)%name, req, detected%linear)) cycle
            any_req = .true.
            dev_req = min(dev_req, maxval(dev_op(cands(ic)%ops)))
            if (.not. all(pass_op(cands(ic)%ops))) cycle
            ! The requested group and every subgroup of it
            do jc = 1, size(cands)
               if (all(pass_op(cands(jc)%ops)) .and. is_subset(cands(jc)%ops, cands(ic)%ops)) then
                  in_pool(jc) = .true.
               end if
            end do
         end do
         if (.not. any_req) then
            call fatal_error(error, "Symmetry: requested group "//self%group_requested// &
               & " is not a subgroup of the detected group "//detected%name)
            return
         end if
         if (.not. any(in_pool)) then
            if (dev_req > tolerance) then
               call fatal_error(error, "Symmetry: geometry deviates from "//self%group_requested// &
                  & " by "//to_string(dev_req, "(es9.2)")// &
                  & " bohr, above the tolerance "//to_string(tolerance, "(es9.2)")//" bohr")
            else
               ! Within the tolerance, but two atoms share an image
               call fatal_error(error, "Symmetry: the tolerance "// &
                  & to_string(tolerance, "(es9.2)")// &
                  & " bohr is too loose for "//self%group_requested// &
                  & ", its operations do not map the atoms one to one")
            end if
            return
         end if
      end if

      ! Best embeddable candidate; candidates whose keys tie (the two C3v
      ! subgroups of C6v, related by an operation that keeps F but is no
      ! symmetry) go to the smaller orbit forms, not to libmsym's order
      allocate (key(4, size(cands)), tried(size(cands)), tie(size(cands)))
      do ic = 1, size(cands)
         key(:, ic) = candidate_key(gops(:, :, cands(ic)%ops), self%frame%axes)
      end do
      tried = .not. in_pool
      best = 0
      do while (.not. all(tried) .and. best == 0)
         ic = best_candidate(key, tried)
         tie = .not. tried
         do jc = 1, size(cands)
            if (tie(jc)) tie(jc) = all(abs(key(:, jc) - key(:, ic)) <= key_tol)
         end do
         tried = tried .or. tie
         do jc = 1, size(cands)
            if (.not. tie(jc)) cycle
            trial = self
            call orient_orbits(trial, xyz, species, gops(:, :, cands(jc)%ops), perm_op(:, cands(jc)%ops), &
               & embeddable)
            if (.not. embeddable) cycle
            if (best /= 0) then
               if (config_order(trial, chosen, xyz, species) >= 0) cycle
            end if
            best = jc
            chosen = trial
         end do
      end do
      if (best /= 0) self = chosen

      call time_phase(timer, "Orientation", .false.)
      if (best == 0) then
         call set_trivial(self, nat)
         if (.not. any(in_pool) .and. self%deviation_detected > tolerance) then
            self%reason = "geometry matches "//detected%name//" only to "// &
               & to_string(self%deviation_detected, "(es9.2)")// &
               & " bohr, above the tolerance "//to_string(tolerance, "(es9.2)")// &
               & " bohr; no symmetry used"
         else if (.not. any(in_pool)) then
            ! Within the tolerance, but two atoms share an image
            self%reason = "the tolerance "//to_string(tolerance, "(es9.2)")// &
               & " bohr is too loose for "//detected%name// &
               & ", its operations do not map the atoms one to one; no symmetry used"
         else
            self%reason = "no subgroup of "//detected%name//" admits Lebedev-invariant site "// &
               & "orientations; no symmetry used"
         end if
         return
      end if

      self%group_used = cands(best)%name
      self%deviation = maxval(dev_op(cands(best)%ops))
      if (auto .and. (best /= 1 .or. detected%linear)) then
         self%reason = "detected "//detected%name//", using "//cands(best)%name
      else if (.not. auto .and. to_lower(cands(best)%name) /= req) then
         self%reason = "requested "//self%group_requested//", using "//cands(best)%name
      end if
      if (len(self%reason) > 0) then
         if (detected%linear) then
            self%reason = self%reason//" (finite subgroup of a linear group)"
         else if (auto .and. any(dev_op(cands(1)%ops) > tolerance)) then
            self%reason = self%reason//" ("//detected%name//" matches only to "// &
               & to_string(self%deviation_detected, "(es9.2)")// &
               & " bohr, above the tolerance "//to_string(tolerance, "(es9.2)")//" bohr)"
         else if (auto .and. any(.not. pass_op(cands(1)%ops))) then
            self%reason = self%reason//" (the tolerance "//to_string(tolerance, "(es9.2)")// &
               & " bohr is too loose for "//detected%name// &
               & ", its operations do not map the atoms one to one)"
         else
            self%reason = self%reason//" (a site group has no Lebedev-invariant orientation)"
         end if
      end if
   end subroutine symmetric_orient

   !> Orient every orbit for a group, or report that a site is not embeddable
   !>
   !> - every member offers its admissible orientations tied closest to F;
   !>   the offers within `tie_tol` of the orbit's best, carried to the first
   !>   such member (the representative), are the orbit's candidates (the
   !>   lowest index alone would make the result depend on the atom order,
   !>   since F itself is not invariant under the group)
   !> - several candidates: the smallest neighbour form at the representative
   !> - equal smallest forms (a symmetry of the molecule outside the group
   !>   used relates the candidates and fixes the site: BF3's B, the metal of
   !>   a sandwich and its axial ligands): the orbits are settled in an order
   !>   that does not depend
   !>   on the atoms' numbering (species, site-group order, distance from the
   !>   centroid), and the tied candidates are compared again by the form
   !>   that also holds the spheres already settled. The first such orbit
   !>   has nothing to tell its candidates apart and takes any; every later
   !>   one follows it, so the cavity is that of any other numbering or pose
   !>   up to one symmetry of the molecule
   !> - orbits that agree in species, site-group order and distance are
   !>   settled in atom order
   !>
   !> @param[in,out] self       Orientations, receives rotations on success
   !> @param[in]     xyz        Coordinates (3, nat)
   !> @param[in]     species    Species ids (nat)
   !> @param[in]     ops        Operations (3, 3, order)
   !> @param[in]     perm       Image atoms (nat, order)
   !> @param[out]    embeddable Whether every site group admits an orientation
   subroutine orient_orbits(self, xyz, species, ops, perm, embeddable)
      !> Orientations
      type(sym_orient_type), intent(inout) :: self
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Operations
      real(wp), intent(in) :: ops(:, :, :)
      !> Image atoms
      integer, intent(in) :: perm(:, :)
      !> Whether every site is embeddable
      logical, intent(out) :: embeddable

      real(wp), allocatable :: rots(:, :, :), dws(:, :, :)
      integer, allocatable :: rot_of_atom(:), site(:), orbit_of_atom(:), orbit_rep(:), member(:)
      !> Species and site-group order of every atom in one number
      integer, allocatable :: marked(:)
      logical, allocatable :: done(:)
      !> Offers of the members of one orbit
      type(offer_type), allocatable :: offer(:)
      !> Candidates of every orbit at its representative
      type(offer_type), allocatable :: cand(:)
      !> First atom of every orbit, and the order the orbits are settled in
      integer, allocatable :: first(:), order(:)
      !> Orientation of every orbit at its representative
      real(wp), allocatable :: rep_rot(:, :, :), rep_dw(:, :, :)
      !> Atoms and the axis ends of the spheres settled so far, with marks
      real(wp), allocatable :: points(:, :)
      integer, allocatable :: point_mark(:)
      !> Candidates of smallest form, and the tied ones among them
      logical, allocatable :: smallest(:), finer(:)
      integer, allocatable :: tied(:)
      real(wp) :: r_rep(3, 3), dw_rep(3, 3), r_b(3, 3), dw_b(3, 3), sgn, best_value
      integer :: nat, norder, ia, ib, ig, jg, ir, nrot, nsite, irep, norbit, nmember, k, ic, jc, ncand
      integer :: io, npoint, i
      logical :: found, known

      nat = size(perm, 1)
      norder = size(ops, 3)
      allocate (rots(3, 3, nat), dws(3, 3, nat), rot_of_atom(nat), done(nat), site(norder), &
         & orbit_of_atom(nat), orbit_rep(nat), member(norder), offer(norder), cand(nat), first(nat))
      done = .false.
      nrot = 0
      norbit = 0
      embeddable = .false.
      ! Species marked by the order of the site group under the operations
      ! used: a symmetry of the molecule outside the group (the other mirror
      ! planes of ammonia in Cs) then no longer makes two candidates look
      ! alike that differ relative to the group
      allocate (marked(nat))
      do ia = 1, nat
         marked(ia) = species(ia)*(norder + 1) + count(perm(ia, :) == ia)
      end do

      ! Candidates of every orbit
      do ia = 1, nat
         if (done(ia)) cycle
         norbit = norbit + 1
         first(norbit) = ia
         nmember = 0
         do ig = 1, norder
            ib = perm(ia, ig)
            if (any(perm(ia, 1:ig - 1) == ib)) cycle
            done(ib) = .true.
            orbit_of_atom(ib) = norbit
            nsite = 0
            do jg = 1, norder
               if (perm(ib, jg) == ib) then
                  nsite = nsite + 1
                  site(nsite) = jg
               end if
            end do
            nmember = nmember + 1
            member(nmember) = ib
            call site_candidates(ops(:, :, site(1:nsite)), self%frame%axes, found, offer(nmember)%value, &
               & offer(nmember)%rot, offer(nmember)%dw)
            if (.not. found) return
         end do

         ! Tied offers, carried to the representative by det(g) g, one per
         ! class modulo Oh; the representative's own come first
         best_value = maxval(offer(1:nmember)%value)
         ncand = 0
         do k = 1, nmember
            if (offer(k)%value >= best_value - tie_tol) ncand = ncand + size(offer(k)%rot, 3)
         end do
         allocate (cand(norbit)%rot(3, 3, ncand), cand(norbit)%dw(3, 3, ncand))
         irep = 0
         ncand = 0
         do k = 1, nmember
            if (offer(k)%value < best_value - tie_tol) cycle
            ib = member(k)
            if (irep == 0) irep = ib
            do ig = 1, norder
               if (perm(ib, ig) == irep) exit
            end do
            ! No operation back to the representative: the images do not
            ! form a group action (tolerance too loose)
            if (ig > norder) return
            sgn = matdet_3x3(ops(:, :, ig))
            do ic = 1, size(offer(k)%rot, 3)
               r_b = sgn*matmul(ops(:, :, ig), offer(k)%rot(:, :, ic))
               known = .false.
               do jc = 1, ncand
                  known = is_oh_element(matmul(transpose(cand(norbit)%rot(:, :, jc)), r_b))
                  if (known) exit
               end do
               if (known) cycle
               ncand = ncand + 1
               cand(norbit)%rot(:, :, ncand) = r_b
               cand(norbit)%dw(:, :, ncand) = sgn*matmul(ops(:, :, ig), offer(k)%dw(:, :, ic))
            end do
         end do
         cand(norbit)%rot = cand(norbit)%rot(:, :, 1:ncand)
         cand(norbit)%dw = cand(norbit)%dw(:, :, 1:ncand)
         orbit_rep(norbit) = irep
         ! No operation from the representative to the first atom: no group action
         if (.not. any(perm(irep, :) == ia)) return
      end do

      ! One candidate per orbit, in an order free of the atoms' numbering
      allocate (rep_rot(3, 3, norbit), rep_dw(3, 3, norbit), order(norbit))
      allocate (points(3, 7*nat), point_mark(7*nat))
      points(:, 1:nat) = xyz
      point_mark(1:nat) = marked
      point_mark(nat + 1:) = maxval(marked) + 1
      npoint = nat
      call settle_order(xyz, marked, first(1:norbit), order)
      do i = 1, norbit
         io = order(i)
         irep = orbit_rep(io)
         ncand = size(cand(io)%rot, 3)
         ic = 1
         if (ncand > 1) then
            smallest = smallest_forms(xyz, marked, irep, cand(io)%rot)
            ic = findloc(smallest, .true., dim=1)
            ! Equal forms: a symmetry of the molecule relates the candidates
            ! (BF3's B: R and sigma_h R; the rings of benzene in C6h). The
            ! spheres settled so far tell them apart, unless that symmetry
            ! keeps those too: then either candidate gives the same cavity
            if (count(smallest) > 1 .and. npoint > nat) then
               tied = pack([(jc, jc=1, ncand)], smallest)
               finer = smallest_forms(points(:, 1:npoint), point_mark(1:npoint), irep, cand(io)%rot(:, :, tied))
               ic = tied(findloc(finer, .true., dim=1))
            end if
         end if
         rep_rot(:, :, io) = cand(io)%rot(:, :, ic)
         rep_dw(:, :, io) = cand(io)%dw(:, :, ic)
         ! The axis ends of the orbit's spheres join the points: as a set they
         ! hold a sphere's orientation modulo Oh
         do ig = 1, norder
            ib = perm(irep, ig)
            if (any(perm(irep, 1:ig - 1) == ib)) cycle
            ! Orbits that overlap hold more than nat members: no group action
            if (npoint + 6 > size(points, 2)) return
            r_b = matdet_3x3(ops(:, :, ig))*matmul(ops(:, :, ig), rep_rot(:, :, io))
            do k = 1, 3
               points(:, npoint + 1) = xyz(:, ib) + r_b(:, k)
               points(:, npoint + 2) = xyz(:, ib) - r_b(:, k)
               npoint = npoint + 2
            end do
         end do
      end do

      ! Rotations of all members, numbered in atom order
      done = .false.
      do io = 1, norbit
         ia = first(io)
         irep = orbit_rep(io)
         ig = findloc(perm(irep, :), ia, dim=1)
         ! g_ia maps the representative onto ia; every member b = g(ia) is
         ! reached from the representative by g g_ia
         r_rep = matdet_3x3(ops(:, :, ig))*matmul(ops(:, :, ig), rep_rot(:, :, io))
         dw_rep = matdet_3x3(ops(:, :, ig))*matmul(ops(:, :, ig), rep_dw(:, :, io))
         do ig = 1, norder
            ib = perm(ia, ig)
            if (done(ib)) cycle
            done(ib) = .true.
            sgn = matdet_3x3(ops(:, :, ig))
            r_b = sgn*matmul(ops(:, :, ig), r_rep)
            dw_b = sgn*matmul(ops(:, :, ig), dw_rep)
            rot_of_atom(ib) = 0
            do ir = 1, nrot
               if (same_modulo_oh(rots(:, :, ir), r_b) .and. &
                  & all(abs(dws(:, :, ir) - dw_b) < dedupe_tol)) then
                  rot_of_atom(ib) = ir
                  exit
               end if
            end do
            if (rot_of_atom(ib) == 0) then
               nrot = nrot + 1
               rots(:, :, nrot) = r_b
               dws(:, :, nrot) = dw_b
               rot_of_atom(ib) = nrot
            end if
         end do
      end do

      embeddable = .true.
      self%nrot = nrot
      self%rotation = rots(:, :, 1:nrot)
      self%dw = dws(:, :, 1:nrot)
      self%rot_of_atom = rot_of_atom
      self%norbit = norbit
      self%orbit_of_atom = orbit_of_atom
      self%orbit_rep = orbit_rep(1:norbit)
      self%ops = ops
      self%perm = perm
   end subroutine orient_orbits

   !> Order of the orbits that does not depend on the atoms' numbering
   !>
   !> By mark, then by the distance from the centroid in clusters of
   !> `radial_tol` relative width; orbits equal in both keep their atom order
   !>
   !> @param[in]  xyz    Coordinates (3, nat)
   !> @param[in]  marked Species and site-group order of every atom (nat)
   !> @param[in]  first  First atom of every orbit (norbit)
   !> @param[out] order  Orbits in the order they are settled (norbit)
   pure subroutine settle_order(xyz, marked, first, order)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Marks
      integer, intent(in) :: marked(:)
      !> First atoms
      integer, intent(in) :: first(:)
      !> Order
      integer, intent(out) :: order(:)

      real(wp), parameter :: radial_tol = 1.0e-7_wp
      real(wp), allocatable :: radius(:)
      integer, allocatable :: key(:, :), by_radius(:)
      real(wp) :: center(3), eps
      integer :: n, i

      n = size(first)
      center = sum(xyz, dim=2)/real(size(xyz, 2), wp)
      allocate (radius(n), key(2, n), by_radius(n))
      do i = 1, n
         radius(i) = norm2(xyz(:, first(i)) - center)
         key(1, i) = marked(first(i))
      end do
      eps = radial_tol*max(maxval(radius), 1.0_wp)
      call sort_keys(key(1:1, :), by_radius, radius)
      do i = 1, n
         if (i == 1) then
            key(2, by_radius(i)) = 1
         else if (radius(by_radius(i)) - radius(by_radius(i - 1)) > eps) then
            key(2, by_radius(i)) = key(2, by_radius(i - 1)) + 1
         else
            key(2, by_radius(i)) = key(2, by_radius(i - 1))
         end if
      end do
      call sort_keys(key, order)
   end subroutine settle_order

   !> D4h or C4v standing in for a linear group
   !>
   !> Axis from libmsym's infinite rotation; the perpendicular axes from the
   !> molecular frame (any choice is a rigid rotation about the molecular axis)
   !>
   !> @param[in,out] self     Orientations
   !> @param[in]     detected Linear group
   !> @param[out]    gops     Operations of the finite group (3, 3, 16 or 8)
   !> @param[out]    cands    The finite group and its subgroups
   subroutine linear_group(self, detected, gops, cands)
      !> Orientations
      type(sym_orient_type), intent(inout) :: self
      !> Linear group
      type(msym_group_type), intent(in) :: detected
      !> Operations
      real(wp), allocatable, intent(out) :: gops(:, :, :)
      !> Candidates
      type(candidate_type), allocatable, intent(out) :: cands(:)

      real(wp) :: axis(3), x(3), y(3), basis(3, 3), local(3, 3, 16), c4(3, 3), c2x(3, 3), sx(3, 3)
      integer :: iop, n, k

      self%ops_on_frame = .true.
      axis = self%frame%axes(:, 3)
      do iop = 1, size(detected%info, 2)
         if (detected%info(1, iop) == msym_op_proper .and. detected%info(2, iop) == 0) then
            axis = detected%axis(:, iop)
            exit
         end if
      end do
      x = self%frame%axes(:, 1) - axis*dot_product(axis, self%frame%axes(:, 1))
      if (norm2(x) < 1.0e-8_wp) x = self%frame%axes(:, 2) - axis*dot_product(axis, self%frame%axes(:, 2))
      x = x/norm2(x)
      y = cross_product(axis, x)
      basis(:, 1) = x
      basis(:, 2) = y
      basis(:, 3) = axis

      ! Generators in the local frame (z the molecular axis)
      c4 = reshape([0.0_wp, 1.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      c2x = reshape([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp], [3, 3])
      sx = reshape([-1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      n = 0
      do k = 0, 3
         n = n + 1
         local(:, :, n) = matpow(c4, k)
      end do
      if (detected%kind == msym_type_dnh) then
         ! D4h = C4 x {E, C2'} x {E, i}
         do k = 1, 4
            local(:, :, 4 + k) = matmul(local(:, :, k), c2x)
         end do
         do k = 1, 8
            local(:, :, 8 + k) = -local(:, :, k)
         end do
         n = 16
      else
         ! C4v = C4 x {E, sigma_v}
         do k = 1, 4
            local(:, :, 4 + k) = matmul(local(:, :, k), sx)
         end do
         n = 8
      end if
      allocate (gops(3, 3, n))
      do k = 1, n
         gops(:, :, k) = matmul(basis, matmul(local(:, :, k), transpose(basis)))
      end do
      call finite_subgroups(local(:, :, 1:n), cands)
   end subroutine linear_group

   !> A subgroup of D4h and all its subgroups but C1 as candidates
   !>
   !> - closures of every one, two and three elements, on the product table
   !> - the group itself comes first, the rest by falling order
   !>
   !> @param[in]  ops   Operations (3, 3, order), z the fourfold axis
   !> @param[out] cands The group and its subgroups
   subroutine finite_subgroups(ops, cands)
      !> Operations
      real(wp), intent(in) :: ops(:, :, :)
      !> Candidates
      type(candidate_type), allocatable, intent(out) :: cands(:)

      integer :: n, i, j, k, l, ng, ig, count_new
      integer, allocatable :: table(:, :)
      logical, allocatable :: member(:), groups(:, :)
      logical :: grew
      real(wp) :: prod(3, 3)

      n = size(ops, 3)
      allocate (table(n, n), member(n), groups(n, n*n*n + 1))
      do i = 1, n
         do j = 1, n
            prod = matmul(ops(:, :, i), ops(:, :, j))
            table(i, j) = 0
            do k = 1, n
               if (all(abs(prod - ops(:, :, k)) < 1.0e-8_wp)) table(i, j) = k
            end do
         end do
      end do

      ng = 1
      groups(:, 1) = .true.
      do i = 1, n
         do j = i, n
            do k = j, n
               member = .false.
               member([i, j, k]) = .true.
               grew = .true.
               do while (grew)
                  grew = .false.
                  do l = 1, n
                     if (.not. member(l)) cycle
                     do ig = 1, n
                        if (.not. member(ig)) cycle
                        if (.not. member(table(l, ig))) then
                           member(table(l, ig)) = .true.
                           grew = .true.
                        end if
                     end do
                  end do
               end do
               if (count(member) < 2) cycle
               count_new = 0
               do ig = 1, ng
                  if (all(groups(:, ig) .eqv. member)) count_new = 1
               end do
               if (count_new == 1) cycle
               ng = ng + 1
               groups(:, ng) = member
            end do
         end do
      end do

      allocate (cands(ng))
      l = 0
      do k = n, 2, -1
         do ig = 1, ng
            if (count(groups(:, ig)) /= k) cycle
            l = l + 1
            cands(l)%ops = pack([(i, i=1, n)], groups(:, ig))
            cands(l)%name = trim(subgroup_name(ops(:, :, cands(l)%ops)))
         end do
      end do
   end subroutine finite_subgroups

   !> Schoenflies symbol of a subgroup of D4h
   !>
   !> @param[in] ops Operations (3, 3, order)
   function subgroup_name(ops) result(name)
      !> Operations
      real(wp), intent(in) :: ops(:, :, :)
      !> Symbol, blank-padded
      character(len=subgroup_name_len) :: name

      integer :: k, nproper
      logical :: inversion, c4, s4, square(size(ops, 3))

      nproper = 0
      inversion = .false.
      do k = 1, size(ops, 3)
         if (matdet_3x3(ops(:, :, k)) > 0.0_wp) nproper = nproper + 1
         if (all(abs(ops(:, :, k) + identity3) < 1.0e-8_wp)) inversion = .true.
         square(k) = all(abs(matmul(ops(:, :, k), ops(:, :, k)) - identity3) < 1.0e-8_wp)
      end do
      c4 = .false.
      s4 = .false.
      do k = 1, size(ops, 3)
         if (square(k)) cycle
         if (matdet_3x3(ops(:, :, k)) > 0.0_wp) then
            c4 = .true.
         else
            s4 = .true.
         end if
      end do

      select case (size(ops, 3))
      case (2)
         if (inversion) then
            name = "Ci"
         else if (nproper == 2) then
            name = "C2"
         else
            name = "Cs"
         end if
      case (4)
         if (c4) then
            name = "C4"
         else if (s4) then
            name = "S4"
         else if (nproper == 4) then
            name = "D2"
         else if (inversion) then
            name = "C2h"
         else
            name = "C2v"
         end if
      case (8)
         if (nproper == 8) then
            name = "D4"
         else if (inversion .and. c4) then
            name = "C4h"
         else if (inversion) then
            name = "D2h"
         else if (c4) then
            name = "C4v"
         else
            name = "D2d"
         end if
      case default
         name = "D4h"
      end select
   end function subgroup_name

   !> Image atoms of an operation and the largest deviation
   !>
   !> - the image of atom a is the nearest atom of its species to
   !>   center + op (x_a - center), the lowest index among equally near ones
   !> - an image lies as far from the center as a, and |x_b - target| >=
   !>   |r_b - r_a|, so the atoms are visited by increasing |r_b - r_a| and
   !>   the search stops once that bound exceeds the nearest distance found
   !>   (a rounding margin keeps it exact): about the atoms of a's shell
   !>   instead of all atoms
   !>
   !> @param[in]  xyz       Coordinates (3, nat)
   !> @param[in]  species   Species ids (nat)
   !> @param[in]  center    Fixed point
   !> @param[in]  op        Operation
   !> @param[in]  radius    Distance of every atom from the center (nat)
   !> @param[in]  by_radius Atoms by increasing radius (nat)
   !> @param[in]  rank      Place of every atom in `by_radius` (nat)
   !> @param[out] image     Nearest atom of the same species to each image (nat)
   !> @param[out] deviation Largest distance of an image to its atom
   pure subroutine atom_images(xyz, species, center, op, radius, by_radius, rank, image, deviation)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Fixed point
      real(wp), intent(in) :: center(3)
      !> Operation
      real(wp), intent(in) :: op(3, 3)
      !> Distance from the center
      real(wp), intent(in) :: radius(:)
      !> Atoms by radius
      integer, intent(in) :: by_radius(:)
      !> Place in `by_radius`
      integer, intent(in) :: rank(:)
      !> Image atoms
      integer, intent(out) :: image(:)
      !> Largest deviation
      real(wp), intent(out) :: deviation

      !> Rounding margin of the radial bound, bohr
      real(wp), parameter :: margin = 1.0e-9_wp
      real(wp) :: target(3), dist, best, r_a, gap_lo, gap_hi
      integer :: ia, ib, nat, lo, hi

      nat = size(xyz, 2)
      deviation = 0.0_wp
      do ia = 1, nat
         target = center + matmul(op, xyz(:, ia) - center)
         r_a = radius(ia)
         best = huge(1.0_wp)
         image(ia) = ia
         lo = rank(ia)
         hi = rank(ia) + 1
         do
            ! Next atom by radial gap, from below or above
            gap_lo = huge(1.0_wp)
            gap_hi = huge(1.0_wp)
            if (lo >= 1) gap_lo = r_a - radius(by_radius(lo))
            if (hi <= nat) gap_hi = radius(by_radius(hi)) - r_a
            if (min(gap_lo, gap_hi) > best + margin) exit
            if (gap_lo <= gap_hi) then
               ib = by_radius(lo)
               lo = lo - 1
            else
               ib = by_radius(hi)
               hi = hi + 1
            end if
            if (species(ib) /= species(ia)) cycle
            dist = norm2(xyz(:, ib) - target)
            if (dist < best .or. (dist == best .and. ib < image(ia))) then
               best = dist
               image(ia) = ib
            end if
         end do
         deviation = max(deviation, best)
      end do
   end subroutine atom_images

   !> Fit the operations to the atoms as one rigid body
   !>
   !> - libmsym builds an exact group in its own frame and places it by axes
   !>   found at its thresholds, so the set as a whole can be turned off the
   !>   atoms (1e-7 rad seen for thousands of atoms) although the geometry is
   !>   symmetric; one rotation Q of the set, g -> Q g Q^T, keeps the group
   !>   closed
   !> - Q minimizes the squared distances of all images to their atoms over
   !>   the fitted operations (Gauss-Newton on the rotation vector; directions
   !>   the operations do not fix, such as the turn about the axis of Cn, stay
   !>   untouched)
   !> - operations that still miss the tolerance leave the fit, the worst
   !>   ones first, so an exact subgroup is not traded against operations the
   !>   geometry does not have
   !> - the fitted set is taken when no operation that passed is lost and
   !>   either more pass or the largest deviation did not grow; the caller's
   !>   set is kept otherwise
   !>
   !> @param[in]     xyz       Coordinates (3, nat)
   !> @param[in]     species   Species ids (nat)
   !> @param[in]     center    Fixed point
   !> @param[in]     radius    Distance of every atom from the center (nat)
   !> @param[in]     by_radius Atoms by increasing radius (nat)
   !> @param[in]     rank      Place of every atom in `by_radius` (nat)
   !> @param[in]     tolerance Exactness tolerance, bohr
   !> @param[in,out] ops       Operations (3, 3, nop)
   !> @param[in,out] perm      Image atoms per operation (nat, nop)
   !> @param[in,out] dev       Largest deviation per operation (nop)
   !> @param[in,out] pass      Whether an operation is exact (nop)
   subroutine refine_operations(xyz, species, center, radius, by_radius, rank, tolerance, ops, perm, dev, pass)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Fixed point
      real(wp), intent(in) :: center(3)
      !> Distance from the center
      real(wp), intent(in) :: radius(:)
      !> Atoms by radius
      integer, intent(in) :: by_radius(:)
      !> Place in `by_radius`
      integer, intent(in) :: rank(:)
      !> Exactness tolerance
      real(wp), intent(in) :: tolerance
      !> Operations
      real(wp), intent(inout) :: ops(:, :, :)
      !> Image atoms
      integer, intent(inout) :: perm(:, :)
      !> Largest deviations
      real(wp), intent(inout) :: dev(:)
      !> Exact operations
      logical, intent(inout) :: pass(:)

      !> Rounds of dropping operations and Gauss-Newton steps per round
      integer, parameter :: max_round = 8, max_step = 4
      !> Rotation below which the next step is rounding, rad
      real(wp), parameter :: step_tol = 1.0e-8_wp
      real(wp), allocatable :: trial(:, :, :), fit_dev(:), new_dev(:)
      integer, allocatable :: new_perm(:, :)
      logical, allocatable :: active(:), new_pass(:)
      real(wp) :: w(3), q(3, 3), worst
      integer :: nop, iop, round, step

      nop = size(ops, 3)
      allocate (active(nop), fit_dev(nop))
      do iop = 1, nop
         active(iop) = is_permutation(perm(:, iop))
      end do
      trial = ops
      do round = 1, max_round
         if (.not. any(active)) exit
         do step = 1, max_step
            call rigid_step(xyz, center, trial, perm, active, w)
            q = rotation_vector(w)
            do iop = 1, nop
               trial(:, :, iop) = matmul(q, matmul(trial(:, :, iop), transpose(q)))
            end do
            if (norm2(w) <= step_tol) exit
         end do
         fit_dev = 0.0_wp
         do iop = 1, nop
            if (active(iop)) fit_dev(iop) = image_deviation(xyz, center, trial(:, :, iop), perm(:, iop))
         end do
         worst = maxval(fit_dev)
         if (worst <= tolerance) exit
         active = active .and. fit_dev <= max(tolerance, 0.5_wp*worst)
      end do

      allocate (new_perm(size(perm, 1), nop), new_dev(nop), new_pass(nop))
      do iop = 1, nop
         call atom_images(xyz, species, center, trial(:, :, iop), radius, by_radius, rank, &
            & new_perm(:, iop), new_dev(iop))
         new_pass(iop) = new_dev(iop) <= tolerance .and. is_permutation(new_perm(:, iop))
      end do
      if (any(pass .and. .not. new_pass)) return
      if (count(new_pass) == count(pass) .and. maxval(new_dev) > maxval(dev)) return
      ops = trial
      perm = new_perm
      dev = new_dev
      pass = new_pass
   end subroutine refine_operations

   !> Gauss-Newton step of the rigid fit of an operation set
   !>
   !> With Q = exp([w]x), the image of a under Q g Q^T moves by
   !> w x (g a) - g (w x a) to first order; w solves the normal equations of
   !> all pairs, by the pseudo-inverse so that an undetermined direction
   !> stays zero
   !>
   !> @param[in]  xyz    Coordinates (3, nat)
   !> @param[in]  center Fixed point
   !> @param[in]  ops    Operations (3, 3, nop)
   !> @param[in]  perm   Image atoms per operation (nat, nop)
   !> @param[in]  active Operations in the fit (nop)
   !> @param[out] w      Rotation vector, rad
   subroutine rigid_step(xyz, center, ops, perm, active, w)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Fixed point
      real(wp), intent(in) :: center(3)
      !> Operations
      real(wp), intent(in) :: ops(:, :, :)
      !> Image atoms
      integer, intent(in) :: perm(:, :)
      !> Operations in the fit
      logical, intent(in) :: active(:)
      !> Rotation vector
      real(wp), intent(out) :: w(3)

      !> Eigenvalues below this fraction of the summed squared radii are
      !> undetermined (the identity and the inversion fix no direction)
      real(wp), parameter :: rank_tol = 1.0e-10_wp
      real(wp) :: hess(3, 3), rhs(3), a(3), ga(3), jac(3, 3), eval(3), work(102), scale
      integer(lapack_ik) :: info
      integer :: iop, ia, k

      hess = 0.0_wp
      rhs = 0.0_wp
      scale = 0.0_wp
      do iop = 1, size(ops, 3)
         if (.not. active(iop)) cycle
         do ia = 1, size(xyz, 2)
            a = xyz(:, ia) - center
            scale = scale + dot_product(a, a)
            ga = matmul(ops(:, :, iop), a)
            jac = matmul(ops(:, :, iop), skew(a)) - skew(ga)
            hess = hess + matmul(transpose(jac), jac)
            rhs = rhs - matmul(transpose(jac), ga - (xyz(:, perm(ia, iop)) - center))
         end do
      end do
      w = 0.0_wp
      call lapack_syev("V", "U", 3_lapack_ik, hess, 3_lapack_ik, eval, work, &
         & int(size(work), lapack_ik), info)
      if (info /= 0) return
      do k = 1, 3
         if (eval(k) > rank_tol*scale) w = w + hess(:, k)*dot_product(hess(:, k), rhs)/eval(k)
      end do
   end subroutine rigid_step

   !> Largest distance of the images under an operation to given atoms
   !>
   !> @param[in] xyz    Coordinates (3, nat)
   !> @param[in] center Fixed point
   !> @param[in] op     Operation
   !> @param[in] image  Image atoms (nat)
   pure function image_deviation(xyz, center, op, image) result(deviation)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Fixed point
      real(wp), intent(in) :: center(3)
      !> Operation
      real(wp), intent(in) :: op(3, 3)
      !> Image atoms
      integer, intent(in) :: image(:)
      !> Largest distance
      real(wp) :: deviation

      integer :: ia

      deviation = 0.0_wp
      do ia = 1, size(xyz, 2)
         deviation = max(deviation, norm2(center + matmul(op, xyz(:, ia) - center) - xyz(:, image(ia))))
      end do
   end function image_deviation

   !> Matrix of the cross product with a vector
   !>
   !> @param[in] v Vector
   pure function skew(v) result(s)
      !> Vector
      real(wp), intent(in) :: v(3)
      !> Matrix with s x = v cross x
      real(wp) :: s(3, 3)

      s = reshape([0.0_wp, v(3), -v(2), -v(3), 0.0_wp, v(1), v(2), -v(1), 0.0_wp], [3, 3])
   end function skew

   !> Rotation by a rotation vector (Rodrigues)
   !>
   !> @param[in] w Rotation vector, rad
   pure function rotation_vector(w) result(q)
      !> Rotation vector
      real(wp), intent(in) :: w(3)
      !> Rotation matrix
      real(wp) :: q(3, 3)

      real(wp) :: angle, k(3, 3)

      q = identity3
      angle = norm2(w)
      if (angle <= tiny(1.0_wp)) return
      k = skew(w/angle)
      q = q + sin(angle)*k + (1.0_wp - cos(angle))*matmul(k, k)
   end function rotation_vector

   !> Whether image atoms are a permutation of 1, ..., nat
   !>
   !> @param[in] image Image atoms (nat)
   pure function is_permutation(image) result(valid)
      !> Image atoms
      integer, intent(in) :: image(:)
      !> Whether every atom is hit exactly once
      logical :: valid

      logical, allocatable :: hit(:)
      integer :: ia

      allocate (hit(size(image)), source=.false.)
      valid = .false.
      do ia = 1, size(image)
         if (image(ia) < 1 .or. image(ia) > size(image)) return
         if (hit(image(ia))) return
         hit(image(ia)) = .true.
      end do
      valid = .true.
   end function is_permutation

   !> Ranking key of a candidate: order, principal order, mirrors, F alignment
   !>
   !> - larger is better, compared in this order: number of operations,
   !>   highest rotation order, number of mirror planes, summed closeness of
   !>   g F to F modulo Oh
   !>
   !> @param[in] ops  Operations (3, 3, order)
   !> @param[in] frame Molecular frame
   function candidate_key(ops, frame) result(key)
      !> Operations
      real(wp), intent(in) :: ops(:, :, :)
      !> Molecular frame
      real(wp), intent(in) :: frame(3, 3)
      !> Key, compared lexicographically
      real(wp) :: key(4)

      real(wp) :: c, theta
      integer :: ig

      key = 0.0_wp
      key(1) = real(size(ops, 3), wp)
      do ig = 1, size(ops, 3)
         if (matdet_3x3(ops(:, :, ig)) > 0.0_wp) then
            c = min(1.0_wp, max(-1.0_wp, (trace3(ops(:, :, ig)) - 1.0_wp)/2.0_wp))
            theta = acos(c)
            if (theta > 1.0e-6_wp) key(2) = max(key(2), real(nint(2.0_wp*pi/theta), wp))
         else if (abs(trace3(ops(:, :, ig)) - 1.0_wp) < 1.0e-6_wp) then
            key(3) = key(3) + 1.0_wp
         end if
         key(4) = key(4) + oh_closeness(frame, matmul(ops(:, :, ig), frame))
      end do
   end function candidate_key

   !> Untried candidate with the largest key
   !>
   !> @param[in] key   Keys (4, ncand)
   !> @param[in] tried Candidates already tried or excluded
   function best_candidate(key, tried) result(best)
      !> Keys
      real(wp), intent(in) :: key(:, :)
      !> Excluded candidates
      logical, intent(in) :: tried(:)
      !> Best candidate
      integer :: best

      integer :: ic, k
      logical :: better

      best = 0
      do ic = 1, size(tried)
         if (tried(ic)) cycle
         if (best == 0) then
            best = ic
            cycle
         end if
         better = .false.
         do k = 1, 4
            if (key(k, ic) > key(k, best) + key_tol) then
               better = .true.
               exit
            else if (key(k, ic) < key(k, best) - key_tol) then
               exit
            end if
         end do
         if (better) best = ic
      end do
   end function best_candidate

   !> Order of two orientations of the same atoms by their orbits
   !>
   !> - per orbit: species of the representative, size, neighbour form of the
   !>   representative's sphere; orbits sorted by these
   !> - fewer orbits first, then the sorted orbit lists lexicographically
   !> - a function of the molecule and the spheres: equal for congruent
   !>   orientations
   !>
   !> @param[in] a       First orientations
   !> @param[in] b       Second orientations
   !> @param[in] xyz     Coordinates (3, nat)
   !> @param[in] species Species ids (nat)
   function config_order(a, b, xyz, species) result(order)
      !> First orientations
      type(sym_orient_type), intent(in) :: a
      !> Second orientations
      type(sym_orient_type), intent(in) :: b
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> -1 when a comes first, 1 when b does, 0 when equal
      integer :: order

      type(neighbour_form_type), allocatable :: fa(:), fb(:)
      integer, allocatable :: sa(:), sb(:)
      integer :: k

      order = 0
      if (a%norbit /= b%norbit) then
         order = merge(-1, 1, a%norbit < b%norbit)
         return
      end if
      call sorted_orbits(a, xyz, species, fa, sa)
      call sorted_orbits(b, xyz, species, fb, sb)
      do k = 1, a%norbit
         order = orbit_order(a, fa, sa(k), b, fb, sb(k), species)
         if (order /= 0) return
      end do
   end function config_order

   !> Neighbour forms of the orbit representatives, and the orbits sorted
   !>
   !> @param[in]  self    Orientations
   !> @param[in]  xyz     Coordinates (3, nat)
   !> @param[in]  species Species ids (nat)
   !> @param[out] forms   Form per orbit (norbit)
   !> @param[out] sorted  Orbits in ascending order (norbit)
   subroutine sorted_orbits(self, xyz, species, forms, sorted)
      !> Orientations
      type(sym_orient_type), intent(in) :: self
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Form per orbit
      type(neighbour_form_type), allocatable, intent(out) :: forms(:)
      !> Orbits in ascending order
      integer, allocatable, intent(out) :: sorted(:)

      integer :: k, j, rep, tmp

      allocate (forms(self%norbit))
      sorted = [(k, k=1, self%norbit)]
      do k = 1, self%norbit
         rep = self%orbit_rep(k)
         call new_neighbour_form(forms(k), xyz, species, rep, self%rotation(:, :, self%rot_of_atom(rep)))
      end do
      ! Insertion sort, stable
      do k = 2, self%norbit
         j = k
         do while (j > 1)
            if (orbit_order(self, forms, sorted(j), self, forms, sorted(j - 1), species) >= 0) exit
            tmp = sorted(j)
            sorted(j) = sorted(j - 1)
            sorted(j - 1) = tmp
            j = j - 1
         end do
      end do
   end subroutine sorted_orbits

   !> Order of two orbits: species of the representative, size, form
   !>
   !> @param[in] a       Orientations of the first orbit
   !> @param[in] fa      Their orbit forms
   !> @param[in] ka      First orbit
   !> @param[in] b       Orientations of the second orbit
   !> @param[in] fb      Their orbit forms
   !> @param[in] kb      Second orbit
   !> @param[in] species Species ids (nat)
   function orbit_order(a, fa, ka, b, fb, kb, species) result(order)
      !> Orientations of the first orbit
      type(sym_orient_type), intent(in) :: a
      !> Their orbit forms
      type(neighbour_form_type), intent(in) :: fa(:)
      !> First orbit
      integer, intent(in) :: ka
      !> Orientations of the second orbit
      type(sym_orient_type), intent(in) :: b
      !> Their orbit forms
      type(neighbour_form_type), intent(in) :: fb(:)
      !> Second orbit
      integer, intent(in) :: kb
      !> Species ids
      integer, intent(in) :: species(:)
      !> -1 when the first comes first, 1 when the second does, 0 when equal
      integer :: order

      integer :: na, nb, sa, sb

      sa = species(a%orbit_rep(ka))
      sb = species(b%orbit_rep(kb))
      na = count(a%orbit_of_atom == ka)
      nb = count(b%orbit_of_atom == kb)
      if (sa /= sb) then
         order = merge(-1, 1, sa < sb)
      else if (na /= nb) then
         order = merge(-1, 1, na < nb)
      else
         order = compare_forms(fa(ka), fb(kb))
      end if
   end function orbit_order

   !> Reverse-mode derivative of the orientations
   !>
   !> Given the torque on every sphere, tau_a = sum_p (p - c_a) x dE/dp over
   !> its points, returns dE/dxyz through the orientations
   !>
   !> With a group a sphere follows the frame within its site (dw) and the
   !> operations for the rest. A displacement dx splits into P dx, its
   !> average over the group, and the rest, of which only the rigid rotation
   !> Omega' = I^+ C (1 - P) dx counts:
   !>
   !> delta omega_r = dw_r J_F (P dx + C^T Omega') + (1 - dw_r) omega_G
   !>
   !> - J_F: frame Jacobian ([[sym_frame_vjp]]); C dx = sum_a d_a x dx_a about
   !>   the centroid; I = C C^T, the inertia tensor of unit masses
   !> - omega_G, the rotation of the operations: Omega', they turn with the
   !>   atoms (P dx keeps them); placed by the frame (`ops_on_frame`, a
   !>   line) they turn with it instead, J_F C^T I^+ C dx
   !> - a frame that turns with the molecule, J_F C^T Omega = Omega, leaves
   !>   J_F P dx within the sites and Omega' for every sphere; a frame whose
   !>   in-plane angle follows the lab (no harmonic) does not
   !> - exact along rigid motions and displacements that keep the group, for
   !>   torques of any symmetry; other displacements break the group and
   !>   contribute their rigid rotation only
   !> - without a group P = 1, dw = 1 and J_F dx remains
   !>
   !> @param[in]  self   Orientations from [[new_sym_orient]] at xyz
   !> @param[in]  xyz    Coordinates (3, nat)
   !> @param[in]  torque Torque per sphere (3, nat)
   !> @param[out] grad   Contribution to dE/dxyz (3, nat)
   subroutine sym_orient_vjp(self, xyz, torque, grad)
      !> Orientations
      type(sym_orient_type), intent(in) :: self
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Torque per sphere
      real(wp), intent(in) :: torque(:, :)
      !> Gradient contribution
      real(wp), intent(out) :: grad(:, :)

      real(wp), allocatable :: tau(:, :)
      !> Coordinates about the centroid, rigid term and its part the group keeps
      real(wp), allocatable :: d(:, :), rigid(:, :), kept(:, :)
      !> Torque on the frame and on the operations, rigid rotation
      real(wp) :: tau_f(3), tau_g(3), w(3), inv(3, 3)
      integer :: ia, ir

      grad = 0.0_wp
      if (self%orientation /= orient_molecular) return

      allocate (tau(3, self%nrot), source=0.0_wp)
      do ia = 1, size(torque, 2)
         ir = self%rot_of_atom(ia)
         tau(:, ir) = tau(:, ir) + torque(:, ia)
      end do
      tau_f = 0.0_wp
      do ir = 1, self%nrot
         tau_f = tau_f + matmul(transpose(self%dw(:, :, ir)), tau(:, ir))
      end do
      call sym_frame_vjp(self%frame, xyz, tau_f, grad)
      if (size(self%ops, 3) <= 1) return

      ! Rigid part of the displacements the group does not keep: the frame
      ! turns by J_F C^T Omega', the operations by Omega'
      allocate (d(3, size(xyz, 2)), rigid(3, size(xyz, 2)), kept(3, size(xyz, 2)))
      d = xyz - spread(sum(xyz, dim=2)/real(size(xyz, 2), wp), 2, size(xyz, 2))
      inv = inertia_pinv(d)
      tau_g = sum(tau, dim=2) - tau_f
      w = moment(d, grad)
      if (.not. self%ops_on_frame) w = w + tau_g
      w = matmul(inv, w)
      call symmetrize_gradient(self, grad)
      do ia = 1, size(xyz, 2)
         rigid(:, ia) = cross_product(w, d(:, ia))
      end do
      kept = rigid
      call symmetrize_gradient(self, kept)
      grad = grad + rigid - kept
      if (.not. self%ops_on_frame) return

      ! Operations placed by the frame turn with it along every rigid rotation
      call sym_frame_vjp(self%frame, xyz, tau_g, rigid)
      w = matmul(inv, moment(d, rigid))
      do ia = 1, size(xyz, 2)
         grad(:, ia) = grad(:, ia) + cross_product(w, d(:, ia))
      end do
   end subroutine sym_orient_vjp

   !> Forward-mode derivative of the orientations, one Jacobian per entry
   !>
   !> jac(k, :, :, r) = d omega_r(k)/d xyz, the left angular velocity of
   !> `rotation(:, :, r)` per atom displacement, split with a group like
   !> [[sym_orient_vjp]]: contracting it with the torques summed per entry
   !> reproduces that VJP
   !>
   !> - zero in the lab orientation
   !>
   !> @param[in]  self Orientations from [[new_sym_orient]] at xyz
   !> @param[in]  xyz  Coordinates (3, nat)
   !> @param[out] jac  Jacobians (3, 3, nat, nrot)
   subroutine sym_orient_jacobian(self, xyz, jac)
      !> Orientations
      type(sym_orient_type), intent(in) :: self
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Jacobians
      real(wp), intent(out) :: jac(:, :, :, :)

      !> Frame Jacobian, d omega_F(m)/d xyz in (:, :, m); with a group less
      !> that of the operations
      real(wp), allocatable :: frame_jac(:, :, :)
      !> Coordinates about the centroid, Jacobian of the operations per
      !> component and a rigid term
      real(wp), allocatable :: d(:, :), ops_jac(:, :, :), rigid(:, :)
      !> Frame rotation per rigid rotation, times the inverse inertia tensor
      real(wp) :: turn(3, 3), inv(3, 3), w(3)
      real(wp) :: unit_torque(3)
      logical :: group
      integer :: ia, ir, k, m

      jac = 0.0_wp
      if (self%orientation /= orient_molecular) return

      group = size(self%ops, 3) > 1
      allocate (frame_jac(3, size(xyz, 2), 3))
      if (group) then
         allocate (d(3, size(xyz, 2)), ops_jac(3, size(xyz, 2), 3), rigid(3, size(xyz, 2)))
         d = xyz - spread(sum(xyz, dim=2)/real(size(xyz, 2), wp), 2, size(xyz, 2))
      end if
      do m = 1, 3
         unit_torque = 0.0_wp
         unit_torque(m) = 1.0_wp
         call sym_frame_vjp(self%frame, xyz, unit_torque, frame_jac(:, :, m))
         if (group) turn(m, :) = moment(d, frame_jac(:, :, m))
      end do

      if (group) then
         ! omega_r = omega_G + dw_r (omega_F - omega_G), with the frame along
         ! P dx + C^T Omega', Omega' = I^+ C (1 - P) dx, and the operations
         ! along Omega' or, placed by the frame, with it along C^T I^+ C dx
         inv = inertia_pinv(d)
         turn = matmul(turn, inv)
         do m = 1, 3
            w = turn(m, :)
            if (.not. self%ops_on_frame) w = w - inv(m, :)
            do ia = 1, size(xyz, 2)
               rigid(:, ia) = cross_product(w, d(:, ia))
            end do
            frame_jac(:, :, m) = frame_jac(:, :, m) - rigid
            call symmetrize_gradient(self, frame_jac(:, :, m))
            if (self%ops_on_frame) then
               ops_jac(:, :, m) = rigid
               cycle
            end if
            frame_jac(:, :, m) = frame_jac(:, :, m) + rigid
            do ia = 1, size(xyz, 2)
               rigid(:, ia) = cross_product(inv(m, :), d(:, ia))
            end do
            ops_jac(:, :, m) = rigid
            call symmetrize_gradient(self, rigid)
            ops_jac(:, :, m) = ops_jac(:, :, m) - rigid
         end do
      end if
      do ir = 1, self%nrot
         do k = 1, 3
            if (group) jac(k, :, :, ir) = ops_jac(:, :, k)
            do m = 1, 3
               if (self%dw(k, m, ir) == 0.0_wp) cycle
               jac(k, :, :, ir) = jac(k, :, :, ir) + self%dw(k, m, ir)*frame_jac(:, :, m)
            end do
         end do
      end do
   end subroutine sym_orient_jacobian

   !> Moment of a set of vectors about the centroid: sum_a d_a x v_a
   !>
   !> @param[in] d Coordinates about the centroid (3, nat)
   !> @param[in] v Vector per atom (3, nat)
   pure function moment(d, v) result(m)
      !> Coordinates about the centroid
      real(wp), intent(in) :: d(:, :)
      !> Vector per atom
      real(wp), intent(in) :: v(:, :)
      !> Moment
      real(wp) :: m(3)

      integer :: ia

      m = 0.0_wp
      do ia = 1, size(d, 2)
         m = m + cross_product(d(:, ia), v(:, ia))
      end do
   end function moment

   !> Pseudo-inverse of the uniform-weight inertia tensor
   !>
   !> I = sum_a (|d_a|^2 1 - d_a d_a^T); Omega = I^+ sum_a d_a x dx_a is the
   !> rigid rotation of a displacement (Eckart, uniform weights like the frame)
   !>
   !> - a line has no rotation about itself and a point none at all: those
   !>   directions (eigenvalues below `inertia_tol` of the largest) are left out
   !>
   !> @param[in] d Coordinates about the centroid (3, nat)
   pure function inertia_pinv(d) result(inv)
      !> Coordinates about the centroid
      real(wp), intent(in) :: d(:, :)
      !> Pseudo-inverse
      real(wp) :: inv(3, 3)

      real(wp) :: evec(3, 3), lambda(3), work(102)
      integer(lapack_ik) :: info
      integer :: ia, i, j

      evec = 0.0_wp
      do ia = 1, size(d, 2)
         do i = 1, 3
            evec(i, i) = evec(i, i) + sum(d(:, ia)**2)
            evec(:, i) = evec(:, i) - d(i, ia)*d(:, ia)
         end do
      end do
      inv = 0.0_wp
      call lapack_syev("V", "U", 3_lapack_ik, evec, 3_lapack_ik, lambda, work, int(size(work), lapack_ik), info)
      if (info /= 0) return
      do j = 1, 3
         if (lambda(j) <= inertia_tol*lambda(3)) cycle
         do i = 1, 3
            inv(:, i) = inv(:, i) + evec(:, j)*evec(i, j)/lambda(j)
         end do
      end do
   end function inertia_pinv

   !> Average a gradient over the group used: (1/|G|) sum_g g^T grad(g(a))
   !>
   !> Projects onto the displacements that keep the group; nothing for C1
   !>
   !> @param[in]     self Orientations
   !> @param[in,out] grad Gradient (3, nat)
   subroutine symmetrize_gradient(self, grad)
      !> Orientations
      type(sym_orient_type), intent(in) :: self
      !> Gradient
      real(wp), intent(inout) :: grad(:, :)

      real(wp), allocatable :: sym(:, :)
      integer :: ia, ig, norder

      norder = size(self%ops, 3)
      if (norder <= 1) return
      allocate (sym(3, size(grad, 2)), source=0.0_wp)
      do ia = 1, size(grad, 2)
         do ig = 1, norder
            sym(:, ia) = sym(:, ia) + matmul(transpose(self%ops(:, :, ig)), grad(:, self%perm(ia, ig)))
         end do
      end do
      grad = sym/real(norder, wp)
   end subroutine symmetrize_gradient

   !> Whether a candidate's name matches the requested symbol
   !>
   !> @param[in] name   Candidate symbol
   !> @param[in] req    Requested symbol, lower case
   !> @param[in] linear Whether the detected group is linear
   pure function same_group_name(name, req, linear) result(same)
      !> Candidate symbol
      character(len=*), intent(in) :: name
      !> Requested symbol
      character(len=*), intent(in) :: req
      !> Linear detected group
      logical, intent(in) :: linear
      !> Whether they match
      logical :: same

      same = to_lower(name) == req
      if (linear .and. .not. same) then
         select case (req)
         case ("dinfh", "dooh", "d*h", "d0h")
            same = name == "D4h"
         case ("cinfv", "coov", "c*v", "c0v")
            same = name == "C4v"
         case default
            same = .false.
         end select
      end if
   end function same_group_name

   !> Whether a lower-case string is a Schoenflies symbol
   !>
   !> - cs, ci, cn, cnv, cnh, dn, dnh, dnd, sn, t, td, th, o, oh, i, ih and
   !>   the linear symbols (c0v, d0h and their spellings)
   !> - n >= 1 for cn, n >= 2 for cnv, cnh and the dihedral groups, n even
   !>   and >= 4 for sn; no leading zero
   !>
   !> @param[in] req Symbol, lower case
   pure function valid_symbol(req) result(valid)
      !> Symbol
      character(len=*), intent(in) :: req
      !> Whether it names a point group
      logical :: valid

      integer :: n, last, order, stat
      character(len=1) :: tail

      valid = .false.
      n = len_trim(req)
      if (n < 1) return
      select case (req(:n))
      case ("cs", "ci", "t", "td", "th", "o", "oh", "i", "ih", "dinfh", "dooh", "d*h", "d0h", &
         & "cinfv", "coov", "c*v", "c0v")
         valid = .true.
         return
      case default
         valid = .false.
      end select
      if (n < 2) return
      if (verify(req(1:1), "cds") /= 0) return
      ! Digits without a leading zero, then at most one letter
      last = verify(req(2:n), "0123456789")
      if (last == 0) last = n
      if (last == 1 .or. last < n - 1 .or. req(2:2) == "0") return
      read (req(2:last), *, iostat=stat) order
      if (stat /= 0) return
      tail = req(last + 1:n)
      select case (req(1:1))
      case ("c")
         ! C1 alone, Cnv and Cnh from n = 2 (C1v and C1h are Cs)
         valid = tail == " " .or. (order >= 2 .and. (tail == "v" .or. tail == "h"))
      case ("d")
         valid = order >= 2 .and. (tail == " " .or. tail == "h" .or. tail == "d")
      case default
         ! Sn: n even from 4 (S2 is Ci), as libmsym
         valid = tail == " " .and. order >= 4 .and. mod(order, 2) == 0
      end select
   end function valid_symbol

   !> Whether two orientations give the same Lebedev sphere: r2 = r1 o, o in Oh
   !>
   !> r1^T r2 a signed permutation matrix to `dedupe_tol`
   !>
   !> @param[in] r1 First orientation
   !> @param[in] r2 Second orientation
   pure function same_modulo_oh(r1, r2) result(same)
      !> First orientation
      real(wp), intent(in) :: r1(3, 3)
      !> Second orientation
      real(wp), intent(in) :: r2(3, 3)
      !> Whether r1^T r2 is in Oh
      logical :: same

      real(wp) :: m(3, 3)

      m = matmul(transpose(r1), r2)
      same = all(abs(m) < dedupe_tol .or. abs(abs(m) - 1.0_wp) < dedupe_tol)
      if (same) same = all(count(abs(m) > 0.5_wp, dim=1) == 1) .and. &
         & all(count(abs(m) > 0.5_wp, dim=2) == 1)
   end function same_modulo_oh

   !> Whether every element of a is in b
   !>
   !> @param[in] a Indices
   !> @param[in] b Indices
   pure function is_subset(a, b) result(yes)
      !> Indices
      integer, intent(in) :: a(:)
      !> Indices
      integer, intent(in) :: b(:)
      !> Whether a is a subset of b
      logical :: yes

      integer :: i

      yes = .true.
      do i = 1, size(a)
         if (.not. any(b == a(i))) then
            yes = .false.
            return
         end if
      end do
   end function is_subset

   !> Integer power of a 3x3 matrix, identity for k = 0
   !>
   !> @param[in] m Matrix
   !> @param[in] k Power, non-negative
   pure function matpow(m, k) result(p)
      !> Matrix
      real(wp), intent(in) :: m(3, 3)
      !> Power
      integer, intent(in) :: k
      !> m^k
      real(wp) :: p(3, 3)

      integer :: i

      p = identity3
      do i = 1, k
         p = matmul(p, m)
      end do
   end function matpow

   !> Trace of a 3x3 matrix
   !>
   !> @param[in] m Matrix
   pure function trace3(m) result(t)
      !> Matrix
      real(wp), intent(in) :: m(3, 3)
      !> Trace
      real(wp) :: t

      t = m(1, 1) + m(2, 2) + m(3, 3)
   end function trace3

end module moist_math_symmetry_orient
