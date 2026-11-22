!> Point-group symmetry and sphere orientations of a structure
!>
!> - host: builds the model from its settings
!>   ([[moist_cavity_grid_model_symmetry_parameters]]) like a radius model,
!>   updates it with the structure and the sphere radii, prints its report
!>   after the radii
!> - cavity: a Lebedev-based one (iSwiG, NumSA, DROP) takes the model
!>   optionally, keeps a copy and updates it per geometry; MT and MC take none
!> - every in-scope cavity uses Lebedev grids, so the model chooses the group
!>   used and the per-atom sphere orientations once the radii are known
!>   ([[moist_math_symmetry_orient]]) and its report is complete
!> - species: atoms are equivalent only with the same atomic number and the
!>   same sphere radius (relative `radius_tol`); ids follow (atomic number,
!>   radius), not the atom order
!> - fingerprint: an update with bitwise the coordinates, atomic numbers,
!>   radii and settings of the latest computation reuses it, so a cavity
!>   updating its copy right after the host never detects twice
!> - lab orientation with group off computes nothing but the identity: no
!>   frame, no libmsym, no timer section
!> - grid convention: the unit point u (a lab-frame Lebedev node) of sphere a
!>   sits at p = c_a + r_a R u with R the rotation of its component; lab
!>   mode stores R = I exactly, so lab points are bitwise unchanged
!> - frame torque: tau = sum_p (p - c_a) x dE/dp over the points of a
!>   component; `add_gradient` maps it to dE/dc through the orientations
!>   (reverse mode), `jacobian` gives d omega/dc per distinct rotation for
!>   forward-mode Jacobians, delta p = delta omega x (p - c_a)
!> - timer: a "Symmetry" section (category setup) under the open node, with
!>   "Molecular frame", "Detection", "Atom matching" and "Orientation"
module moist_cavity_grid_model_symmetry_type
   use, intrinsic :: iso_fortran_env, only: output_unit, int64
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_symbols, only: to_symbol
   use mctc_io_utils, only: to_lower, to_string
   use moist_context, only: moist_context_type
   use moist_utils_timer, only: cat_setup
   use moist_cavity_grid_model_symmetry_parameters, only: moist_cavity_grid_model_symmetry_parameters_type, &
      & check_symmetry_parameters
   use moist_math_symmetry_orient, only: sym_orient_type, new_sym_orient, sym_orient_vjp, &
      & sym_orient_jacobian, orient_lab, orient_molecular
   use moist_math_symmetry_frame, only: sym_threshold_type, frame_level_name, frame_level_list
   use moist_channels_fields, only: field_query_type
   use moist_cavity_grid_model_type, only: grid_model_type, grid_components_type, grid_size_default
   use moist_utils_format, only: format_string
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   implicit none(type, external)
   private

   public :: grid_model_symmetry_type, new_grid_model_symmetry

   !> Relative radius difference below which two spheres count as equal
   real(wp), parameter :: radius_tol = 1.0e-10_wp
   !> Longest atom or rotation list printed in one table cell
   integer, parameter :: max_list_len = 40
   !> What an ambiguous frame means for the result
   character(len=*), parameter :: ambiguous_text = "the atom positions alone do not determine the"// &
      & " grid orientation uniquely (frames tie that the species tell apart), so the result can"// &
      & " depend on the orientation of the molecule"

   !> Point-group symmetry and sphere orientations of the latest structure
   !>
   !> - as a grid model: one component of weight one per sphere, of the
   !>   cavity's Lebedev size, turned by the sphere's rotation
   type, extends(grid_model_type) :: grid_model_symmetry_type
      !> Settings
      type(moist_cavity_grid_model_symmetry_parameters_type) :: param
      !> Whether the state below belongs to the fingerprint
      logical :: computed = .false.
      !> Fingerprint: coordinates of the computation (3, nat), bohr
      real(wp), allocatable :: xyz(:, :)
      !> Fingerprint: atomic number per atom (nat)
      integer, allocatable :: num(:)
      !> Fingerprint: sphere radius per atom (nat), bohr
      real(wp), allocatable :: radii(:)
      !> Fingerprint: settings of the computation
      type(moist_cavity_grid_model_symmetry_parameters_type) :: used
      !> Orientations, groups and orbits
      type(sym_orient_type) :: orient
   contains
      !> Compute the orientations of a structure unless already current
      procedure :: update => update_symmetry
      !> Whether a structure and radii match the latest computation
      procedure :: is_current => symmetry_is_current
      !> Whether the settings ask for nothing: lab orientation, group off
      procedure :: is_trivial => symmetry_is_trivial
      !> Whether the detected group failed the tolerance in auto mode
      procedure :: missed_tolerance => symmetry_missed_tolerance
      !> Whether the molecular frame lies near any switch of its choices
      procedure :: near_switch => symmetry_near_switch
      !> Whether two eigenvalues of the frame tensor are nearly degenerate
      procedure :: near_gap => symmetry_near_gap
      !> Whether the molecule is nearly linear
      procedure :: near_linear => symmetry_near_linear
      !> Whether the in-plane harmonic is weak or nearly tied
      procedure :: near_harmonic => symmetry_near_harmonic
      !> Whether frames tie that the species tell apart
      procedure :: tied_frames => symmetry_tied_frames
      !> Forget the computation: every sphere back in the lab frame
      procedure :: clear => clear_symmetry
      !> Whether every sphere uses the exact identity
      procedure :: is_lab => symmetry_is_lab
      !> Grid model: every sphere keeps its plain grid (lab orientation)
      procedure :: is_plain => symmetry_is_lab
      !> Grid model: the rotated components of every sphere
      procedure :: components => symmetry_components
      !> Grid model: map the component torques and weight adjoints to the atoms
      procedure :: add_gradient => symmetry_add_gradient
      !> Grid model: rotation and weight Jacobians
      procedure :: jacobian => symmetry_jacobian
      !> Declare the readable numeric symmetry results
      procedure :: list_fields => list_symmetry_fields
      !> Print the symmetry report
      procedure :: print => print_symmetry_report
   end type grid_model_symmetry_type

contains

   !* ================================================================================= *!
   !*                                   Construction                                    *!
   !* ================================================================================= *!

   !> Construct from settings; omitted, lab orientation and group off
   !>
   !> @param[out] self  Symmetry model, nothing computed
   !> @param[out] error Invalid setting
   !> @param[in]  param Settings copied by value
   subroutine new_grid_model_symmetry(self, error, param)
      !> Symmetry model
      type(grid_model_symmetry_type), intent(out) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Settings; omitted means compiled defaults
      type(moist_cavity_grid_model_symmetry_parameters_type), intent(in), optional :: param

      if (present(param)) self%param = param
      call check_symmetry_parameters(self%param, error)
   end subroutine new_grid_model_symmetry

   !* ================================================================================= *!
   !*                                    Computation                                    *!
   !* ================================================================================= *!

   !> Orientations of a structure unless the latest computation is current
   !>
   !> - atoms are equivalent only with the same atomic number and radius
   !> - timed as "Symmetry" under the context's open node unless the settings
   !>   are trivial (lab orientation, group off)
   !> - warns through the context when auto mode detects a group the geometry
   !>   misses only by the tolerance, and when the molecular frame lies near a
   !>   switch of its choices
   !> - a failed update leaves nothing computed
   !>
   !> @param[in,out] self  Symmetry model
   !> @param[in]     mol   Structure
   !> @param[in]     radii Sphere radius per atom (nat), bohr
   !> @param[out]    error Invalid setting, explicit group absent or inexact
   !> @param[in,out] ctx   Run context: timer and warnings
   subroutine update_symmetry(self, mol, radii, error, ctx)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(inout) :: self
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Sphere radii
      real(wp), intent(in) :: radii(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Run context
      type(moist_context_type), intent(inout), optional :: ctx

      !> Timer stack depth at entry; error paths unwind back to it
      integer :: d0

      if (self%is_current(mol, radii)) return
      call self%clear()
      call check_symmetry_parameters(self%param, error)
      if (allocated(error)) return
      if (size(radii) /= mol%nat) then
         call fatal_error(error, "Symmetry: radii must have one entry per atom")
         return
      end if

      if (self%is_trivial() .or. .not. present(ctx)) then
         call compute_symmetry(self, mol, radii, error)
         if (allocated(error)) return
      else
         d0 = ctx%timer%current_depth()
         call ctx%timer%start("Symmetry", category=cat_setup)
         call compute_symmetry(self, mol, radii, error, ctx)
         if (allocated(error)) then
            call ctx%timer%unwind(d0)
            return
         end if
         call ctx%timer%stop("Symmetry")
      end if

      if (present(ctx) .and. self%missed_tolerance()) then
         call ctx%warning("Symmetry: the geometry matches "//self%orient%group_detected// &
            & " only to "//format_string(self%orient%deviation_detected, "(es9.2)")// &
            & " bohr, above the tolerance "//format_string(self%used%tolerance, "(es9.2)")// &
            & " bohr, so "//self%orient%group_used// &
            & " is used; symmetrize the input geometry or raise symmetry.tolerance")
      end if
      if (present(ctx) .and. self%tied_frames()) then
         call ctx%warning("Symmetry: "//ambiguous_text//", as in the lab frame; a point group"// &
            & " (symmetry.group) removes the choice only where the grid can carry it")
      end if
      if (present(ctx) .and. self%near_switch()) then
         call ctx%warning("Symmetry: "//near_text(self%orient)//"; the molecular frame is close"// &
            & " to a switch and the grid orientation can change abruptly under a small"// &
            & " displacement, which can cause problems in geometry optimizations and other"// &
            & " calculations; symmetrize the geometry and specify its point group"// &
            & " (symmetry.group), or use the lab frame")
      end if
   end subroutine update_symmetry

   !> Run the computation and record its fingerprint
   !>
   !> @param[in,out] self  Symmetry model, cleared
   !> @param[in]     mol   Structure
   !> @param[in]     radii Sphere radii (nat)
   !> @param[out]    error Explicit group absent or inexact
   !> @param[in,out] ctx   Run context whose timer records the phases
   subroutine compute_symmetry(self, mol, radii, error, ctx)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(inout) :: self
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Sphere radii
      real(wp), intent(in) :: radii(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Run context
      type(moist_context_type), intent(inout), optional :: ctx

      !> Atomic number and species id per atom
      integer, allocatable :: num(:), species(:)
      !> Requested group, lower case
      character(len=:), allocatable :: group
      integer :: orientation
      type(sym_threshold_type) :: threshold

      threshold%gap = self%param%gap_threshold
      threshold%harmonic = self%param%harmonic_threshold
      threshold%linear = self%param%linear_threshold
      orientation = orient_lab
      if (self%param%effective_orientation() == "molecular") orientation = orient_molecular
      num = mol%num(mol%id)
      group = to_lower(trim(adjustl(self%param%group)))
      ! Without a group only the molecular frame asks, for its tied choices
      if ((group == "off" .or. group == "c1") .and. orientation /= orient_molecular) then
         allocate (species(mol%nat), source=1)
      else
         call sphere_species(num, radii, species)
      end if

      if (present(ctx)) then
         call new_sym_orient(self%orient, mol%xyz, species, orientation, trim(adjustl(self%param%group)), &
            & self%param%tolerance, error, ctx%timer, threshold=threshold)
      else
         call new_sym_orient(self%orient, mol%xyz, species, orientation, trim(adjustl(self%param%group)), &
            & self%param%tolerance, error, threshold=threshold)
      end if
      if (allocated(error)) return
      self%xyz = mol%xyz
      self%num = num
      self%radii = radii
      self%used = self%param
      self%computed = .true.
   end subroutine compute_symmetry

   !> Species ids 1..n of the distinct (atomic number, radius) pairs
   !>
   !> - radii closer than `radius_tol` relative count as equal, so roundoff in
   !>   geometry-dependent radii of equivalent atoms does not split them
   !> - numbered by atomic number, then radius, never by the atom order: the
   !>   orientation tie-breaker sorts atoms by species
   !>
   !> @param[in]  num     Atomic numbers (nat)
   !> @param[in]  radii   Sphere radii (nat)
   !> @param[out] species Species id per atom (nat)
   pure subroutine sphere_species(num, radii, species)
      !> Atomic numbers
      integer, intent(in) :: num(:)
      !> Sphere radii
      real(wp), intent(in) :: radii(:)
      !> Species ids
      integer, allocatable, intent(out) :: species(:)

      !> Atoms by atomic number, then radius
      integer, allocatable :: order(:), work(:)
      integer :: n, width, lo, mid, hi, i, j, k, prev, cur

      n = size(num)
      allocate (species(n), order(n), work(n))
      order = [(i, i=1, n)]
      width = 1
      do while (width < n)
         lo = 1
         do while (lo <= n)
            mid = min(lo + width - 1, n)
            hi = min(lo + 2*width - 1, n)
            i = lo
            j = mid + 1
            k = lo
            do while (i <= mid .and. j <= hi)
               if (num(order(j)) < num(order(i)) .or. (num(order(j)) == num(order(i)) &
                  & .and. radii(order(j)) < radii(order(i)))) then
                  work(k) = order(j)
                  j = j + 1
               else
                  work(k) = order(i)
                  i = i + 1
               end if
               k = k + 1
            end do
            work(k:k + mid - i) = order(i:mid)
            k = k + mid - i + 1
            work(k:k + hi - j) = order(j:hi)
            lo = lo + 2*width
         end do
         order = work
         width = 2*width
      end do

      ! A new species where the atomic number changes or the radius steps
      ! beyond roundoff from its neighbour (single linkage, so the grouping
      ! does not depend on the atom order)
      do k = 1, n
         cur = order(k)
         if (k == 1) then
            species(cur) = 1
         else if (num(cur) /= num(prev) .or. abs(radii(cur) - radii(prev)) > &
            & radius_tol*max(abs(radii(cur)), abs(radii(prev)))) then
            species(cur) = species(prev) + 1
         else
            species(cur) = species(prev)
         end if
         prev = cur
      end do
   end subroutine sphere_species

   !> Whether a structure and radii match the latest computation bit for bit
   !>
   !> Coordinates, atomic numbers, radii and the settings
   !>
   !> @param[in] self  Symmetry model
   !> @param[in] mol   Structure
   !> @param[in] radii Sphere radii (nat)
   function symmetry_is_current(self, mol, radii) result(current)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Sphere radii
      real(wp), intent(in) :: radii(:)
      !> Whether the latest computation applies
      logical :: current

      current = .false.
      if (.not. self%computed) return
      if (.not. (allocated(self%xyz) .and. allocated(self%num) .and. allocated(self%radii))) return
      if (size(self%num) /= mol%nat .or. size(radii) /= mol%nat) return
      if (any(shape(mol%xyz) /= shape(self%xyz))) return
      if (any(self%num /= mol%num(mol%id))) return
      if (.not. same_bits(reshape(self%xyz, [size(self%xyz)]), reshape(mol%xyz, [size(mol%xyz)]))) return
      if (.not. same_bits(self%radii, radii)) return
      if (self%used%orientation /= self%param%orientation .or. self%used%group /= self%param%group) return
      current = same_bits([self%used%tolerance, self%used%gap_threshold, self%used%harmonic_threshold, &
         & self%used%linear_threshold], [self%param%tolerance, self%param%gap_threshold, &
         & self%param%harmonic_threshold, self%param%linear_threshold])
   end function symmetry_is_current

   !> Whether the settings ask for nothing: lab orientation and group off
   !>
   !> @param[in] self Symmetry model
   pure function symmetry_is_trivial(self) result(trivial)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Whether only the identity results
      logical :: trivial

      trivial = self%param%effective_orientation() == "lab" .and. .not. self%param%has_group()
   end function symmetry_is_trivial

   !> Whether auto mode detected a group the geometry misses by the tolerance
   !>
   !> The group used is then a proper subgroup or C1
   !>
   !> @param[in] self Symmetry model
   pure function symmetry_missed_tolerance(self) result(missed)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Whether the detected group exceeds the tolerance
      logical :: missed

      missed = .false.
      if (.not. self%computed) return
      if (.not. (allocated(self%orient%group_requested) .and. allocated(self%orient%group_detected))) return
      if (to_lower(self%orient%group_requested) /= "auto" .or. len(self%orient%group_detected) == 0) return
      missed = self%orient%deviation_detected > self%used%tolerance
   end function symmetry_missed_tolerance

   !> Whether the molecular frame lies near any switch of its choices
   !>
   !> - nearly degenerate eigenvalues of the frame tensor, a nearly linear
   !>   molecule, or a weak or nearly tied in-plane harmonic
   !> - an exact degeneracy of a symmetric molecule lies far from the switch
   !>   and is not counted; tied frames are reported by `tied_frames`
   !> - false in the lab orientation and before a computation
   !>
   !> @param[in] self Symmetry model
   pure function symmetry_near_switch(self) result(near)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Whether any choice lies near its switch
      logical :: near

      near = self%near_gap() .or. self%near_linear() .or. self%near_harmonic()
   end function symmetry_near_switch

   !> Whether the molecular frame of the latest computation is in use
   !>
   !> @param[in] self Symmetry model
   pure function has_frame(self) result(has)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Whether the frame flags apply
      logical :: has

      has = self%computed
      if (has) has = self%orient%orientation == orient_molecular
   end function has_frame

   !> Whether two eigenvalues of the frame tensor are nearly degenerate
   !>
   !> @param[in] self Symmetry model
   pure function symmetry_near_gap(self) result(near)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Whether the eigen-gap lies near its threshold
      logical :: near

      near = has_frame(self)
      if (near) near = self%orient%frame%near_gap
   end function symmetry_near_gap

   !> Whether the molecule is nearly linear
   !>
   !> @param[in] self Symmetry model
   pure function symmetry_near_linear(self) result(near)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Whether the distance from the axis lies near its threshold
      logical :: near

      near = has_frame(self)
      if (near) near = self%orient%frame%near_linear
   end function symmetry_near_linear

   !> Whether the in-plane harmonic is weak or nearly tied
   !>
   !> @param[in] self Symmetry model
   pure function symmetry_near_harmonic(self) result(near)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Whether the harmonic lies near its threshold or a second one
      logical :: near

      near = has_frame(self)
      if (near) near = self%orient%frame%near_harmonic
   end function symmetry_near_harmonic

   !> Whether frames tie exactly that the species tell apart
   !>
   !> The atom positions alone then do not fix the grid orientation
   !>
   !> @param[in] self Symmetry model
   pure function symmetry_tied_frames(self) result(tied)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Whether the orientation is ambiguous
      logical :: tied

      tied = has_frame(self)
      if (tied) tied = self%orient%ambiguous
   end function symmetry_tied_frames

   !> Text naming the choices of the frame near their switch
   !>
   !> @param[in] orient Orientations with a flagged frame
   function near_text(orient) result(text)
      !> Orientations
      type(sym_orient_type), intent(in) :: orient
      !> Description
      character(len=:), allocatable :: text

      text = ""
      if (orient%frame%near_gap) text = text//", nearly degenerate second-moment eigenvalues"// &
         & " (relative gap "//format_string(orient%frame%near_gap_value, "(es9.2)")//")"
      if (orient%frame%near_linear) text = text//", nearly linear molecule (mean square distance"// &
         & " from the axis "//format_string(orient%frame%tau, "(es9.2)")//" bohr^2)"
      if (orient%frame%near_harmonic) text = text//", weak or nearly tied in-plane harmonic ("// &
         & format_string(orient%frame%near_harmonic_value, "(es9.2)")//")"
      text = text(3:)
   end function near_text

   !> Forget the computation: every sphere back in the lab frame
   !>
   !> @param[in,out] self Symmetry model
   subroutine clear_symmetry(self)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(inout) :: self

      !> Orientations as constructed
      type(sym_orient_type) :: empty

      self%computed = .false.
      self%orient = empty
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%num)) deallocate (self%num)
      if (allocated(self%radii)) deallocate (self%radii)
   end subroutine clear_symmetry

   !* ================================================================================= *!
   !*                                     Accessors                                     *!
   !* ================================================================================= *!

   !> Whether every sphere uses the exact identity
   !>
   !> - lab orientation, or nothing computed
   !>
   !> @param[in] self Symmetry model
   pure function symmetry_is_lab(self) result(lab)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Whether every rotation is the identity
      logical :: lab

      lab = .true.
      if (self%computed) lab = self%orient%orientation == orient_lab
   end function symmetry_is_lab

   !> Grid components of the cavity's Lebedev size: one per sphere, weight one
   !>
   !> @param[in]  self  Updated symmetry model
   !> @param[in]  nsph  Number of spheres
   !> @param[out] comps Components of every sphere and the distinct rotations
   subroutine symmetry_components(self, nsph, comps)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Number of spheres
      integer, intent(in) :: nsph
      !> Components
      type(grid_components_type), intent(out) :: comps

      integer :: isph

      ! TODO: choose the grid family per sphere by its site symmetry;
      ! - Lebedev grids carry only the octahedral group, so sites with a fivefold or
      !   sixfold axis or a mirror the cube lacks are reduced (BF3 D3h -> C3v)
      ! - a product grid (Gauss-Legendre in the polar angle times n equidistant
      !   azimuths, any Cnv/Dnh) or an
      ! - icosahedral grid on those spheres would carry the site group; needs a
      !   per-size and per-family node provider in place of the single Lebedev set
      !   of build_grid
      allocate (comps%sph_first(nsph + 1), comps%rot(nsph))
      allocate (comps%weight(nsph), source=1.0_wp)
      allocate (comps%size(nsph), source=grid_size_default)
      comps%sph_first = [(isph, isph=1, nsph + 1)]
      if (self%computed) then
         comps%rot = self%orient%rot_of_atom(1:nsph)
         comps%rotation = self%orient%rotation(:, :, 1:self%orient%nrot)
      else
         comps%rot = 1
         allocate (comps%rotation(3, 3, 1))
         comps%rotation(:, :, 1) = reshape([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
            & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      end if
   end subroutine symmetry_components

   !> Grid model gradient: the frame torques
   !>
   !> Nothing in lab mode; otherwise dE/dc through the molecular frame and
   !> the orientations; with a group, along displacements that keep it and
   !> along rigid rotations
   !>
   !> @param[in]     self     Symmetry model
   !> @param[in]     torque   Torque per component, one per sphere (3, nat)
   !> @param[in]     wbar     Weight adjoint per component (nat); unused, the
   !>                         weights are constant
   !> @param[in,out] gradient Nuclear gradient (3, nat), accumulated
   subroutine symmetry_add_gradient(self, torque, wbar, gradient)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Torque per component
      real(wp), intent(in) :: torque(:, :)
      !> Weight adjoint per component
      real(wp), intent(in) :: wbar(:)
      !> Nuclear gradient
      real(wp), intent(inout) :: gradient(:, :)

      real(wp), allocatable :: grad(:, :)

      if (self%is_lab()) return
      allocate (grad(3, size(gradient, 2)))
      call sym_orient_vjp(self%orient, self%xyz, torque, grad)
      gradient = gradient + grad
   end subroutine symmetry_add_gradient

   !> Grid model Jacobians: angular velocities of the rotations
   !>
   !> - `jac_rot(k, :, :, r)` = d omega_r(k)/dc for the distinct rotation R_r,
   !>   delta R_r = [delta omega_r]_x R_r; a point of a component with that
   !>   rotation moves by delta omega_r x (p - c_a)
   !> - the same split with a group as `add_gradient`, which it
   !>   reproduces when contracted with the torques; zero in lab mode
   !>
   !> @param[in]  self         Symmetry model
   !> @param[out] jac_rot      d omega_r(k)/dc_j(i) as (k, i, j, r) (3, 3, nat, nrot)
   !> @param[out] jac_weight   Left unallocated: the weights are constant
   !> @param[out] weight_class Left unallocated
   subroutine symmetry_jacobian(self, jac_rot, jac_weight, weight_class)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Rotation Jacobians
      real(wp), allocatable, intent(out) :: jac_rot(:, :, :, :)
      !> Weight Jacobians
      real(wp), allocatable, intent(out) :: jac_weight(:, :, :)
      !> Weight class of every component
      integer, allocatable, intent(out) :: weight_class(:)

      integer :: nat

      if (.not. self%computed) then
         allocate (jac_rot(3, 3, 0, 1))
         return
      end if
      nat = size(self%xyz, 2)
      allocate (jac_rot(3, 3, nat, self%orient%nrot), source=0.0_wp)
      if (self%is_lab()) return
      call sym_orient_jacobian(self%orient, self%xyz, jac_rot)
   end subroutine symmetry_jacobian

   !* ================================================================================= *!
   !*                                      Output                                       *!
   !* ================================================================================= *!

   !> Declare the readable numeric symmetry results; nothing before a computation
   !>
   !> Groups and the reason are text, which the field query has no type for;
   !> the report prints them
   !>
   !> @param[in]     self  Symmetry model
   !> @param[in,out] query Walker collecting or fetching the declarations
   subroutine list_symmetry_fields(self, query)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Walker collecting or fetching the declarations
      type(field_query_type), intent(inout) :: query

      !> Per-sphere rotations and indices
      real(wp), allocatable :: rotation(:, :, :)
      integer, allocatable :: ids(:)
      !> Scalars in the allocatable form the query takes
      real(wp), allocatable :: scalar
      !> One flag in the array form the query takes
      logical, allocatable :: flag(:)
      integer :: isph

      if (.not. self%computed) return
      associate (orient => self%orient)
         call query%add_int_value("symmetry_orientation", &
            & "Sphere orientation used: "//to_string(orient_lab)//" lab frame, "// &
            & to_string(orient_molecular)//" molecular frame", orient%orientation)
         call query%add_int_value("symmetry_order", "Order of the group used", size(orient%ops, 3))
         scalar = orient%deviation
         call query%add_real_scalar("symmetry_deviation", &
            & "Largest atom deviation under the operations used, bohr", scalar)
         call query%add_int_value("symmetry_frame_level", &
            & "Molecular frame level: "//frame_level_list(), orient%frame%level)
         scalar = orient%frame%gap
         call query%add_real_scalar("symmetry_frame_gap", &
            & "Smallest gap of two neighbouring eigenvalues of the frame tensor relative to their sum", scalar)
         call query%add_int_value("symmetry_nrot", "Number of distinct sphere rotations", orient%nrot)
         flag = [self%near_gap()]
         call query%add_bool("symmetry_near_gap", &
            & "Whether two eigenvalues of the frame tensor are nearly degenerate (1)", flag)
         flag = [self%near_linear()]
         call query%add_bool("symmetry_near_linear", "Whether the molecule is nearly linear (1)", flag)
         flag = [self%near_harmonic()]
         call query%add_bool("symmetry_near_harmonic", &
            & "Whether the in-plane harmonic of the frame is weak or nearly tied (1)", flag)
         flag = [self%tied_frames()]
         call query%add_bool("symmetry_tied_frames", &
            & "Whether frames tie that the species tell apart (1)", flag)
         ids = orient%rot_of_atom
         call query%add_int("symmetry_rot_index", "Distinct rotation of each sphere, 0-based (nsph)", &
            & ids, zero_based=.true.)
         ids = orient%orbit_of_atom
         call query%add_int("symmetry_orbit", "Orbit of each sphere under the group used, 0-based (nsph)", &
            & ids, zero_based=.true.)
         allocate (rotation(3, 3, size(orient%rot_of_atom)))
         do isph = 1, size(orient%rot_of_atom)
            rotation(:, :, isph) = orient%rotation(:, :, orient%rot_of_atom(isph))
         end do
         call query%add_real3("symmetry_rotation", &
            & "Sphere rotation R, unit point u at center + radius R u (3, 3, nsph)", rotation)
      end associate
   end subroutine list_symmetry_fields

   !> Print the symmetry report of the latest computation
   !>
   !> - group off: the orientation and "off"
   !> - otherwise requested, detected (with its deviation) and used group, the
   !>   reason when they differ, deviation or tolerance, molecular frame,
   !>   orbits and rotations
   !> - from `level` 2, a table of the orbits when the group used has any
   !>   operation besides the identity
   !> - nothing before a computation or at `level` 0
   !>
   !> @param[in] self  Symmetry model
   !> @param[in] unit  Output unit (default standard output)
   !> @param[in] level Verbosity level (default everything)
   subroutine print_symmetry_report(self, unit, level)
      !> Symmetry model
      class(grid_model_symmetry_type), intent(in) :: self
      !> Output unit
      integer, intent(in), optional :: unit
      !> Verbosity level
      integer, intent(in), optional :: level

      type(prettyprinter) :: pp
      character(len=:), allocatable :: orientation, used
      integer :: iu, lvl, order

      if (.not. self%computed) return
      iu = output_unit
      if (present(unit)) iu = unit
      lvl = huge(1)
      if (present(level)) lvl = level
      if (lvl < 1) return

      associate (orient => self%orient)
         orientation = "lab"
         if (orient%orientation == orient_molecular) orientation = "molecular"
         if (self%used%orientation_implied()) orientation = orientation//" (implied by group)"

         pp = new_prettyprinter(unit=iu, fmt_len=20)
         call pp%push("Symmetry:")
         call pp%kv("Orientation", orientation)
         if (to_lower(orient%group_requested) == "off") then
            call pp%kv("Group", "off")
            if (orient%orientation == orient_molecular) call print_frames(pp, orient)
            call pp%pop()
            call pp%blank()
            return
         end if

         order = size(orient%ops, 3)
         call pp%kv("Group requested", orient%group_requested)
         if (len(orient%group_detected) > 0) then
            call pp%kv("Group detected", orient%group_detected//" (deviation "// &
               & format_string(orient%deviation_detected, "(es9.2)")//" bohr)")
         else
            call pp%kv("Group detected", "-")
         end if
         used = orient%group_used
         if (order > 1) used = used//" (order "//to_string(order)//")"
         call pp%kv("Group used", used)
         if (len(orient%reason) > 0) call pp%kv("Reason", orient%reason)
         if (order > 1) then
            call pp%kv("Max. deviation", format_string(orient%deviation, "(es9.2)")// &
               & " bohr (tolerance "//format_string(self%used%tolerance, "(es9.2)")//" bohr)")
         else
            call pp%kv("Tolerance", format_string(self%used%tolerance, "(es9.2)")//" bohr")
         end if
         if (orient%orientation == orient_molecular) call print_frames(pp, orient)
         call pp%kv("Orbits", to_string(orient%norbit))
         call pp%kv("Grid rotations", to_string(orient%nrot))
         call pp%pop()
         call pp%blank()
         if (lvl >= 2 .and. order > 1) call print_orbit_table(orient, iu, self%num, self%radii)
      end associate
   end subroutine print_symmetry_report

   !> Report lines of the molecular frame
   !>
   !> - the level and eigen-gap of the frame
   !> - the choices near their switch, when any
   !>
   !> @param[in,out] pp     Printer inside the symmetry section
   !> @param[in]     orient Orientations
   subroutine print_frames(pp, orient)
      !> Printer
      type(prettyprinter), intent(inout) :: pp
      !> Orientations
      type(sym_orient_type), intent(in) :: orient

      call pp%kv("Molecular frame", frame_level_name(orient%frame%level)//" (eigen-gap "// &
         & format_string(orient%frame%gap, "(es9.2)")//")")
      if (orient%frame%near_gap .or. orient%frame%near_linear .or. orient%frame%near_harmonic) then
         call pp%kv("Near a switch", near_text(orient))
      end if
      if (orient%ambiguous) call pp%kv("Tied frames", "the species tell them apart; orientation dependent")
   end subroutine print_frames

   !> Table of the orbits: one row per set of symmetry-equivalent atoms
   !>
   !> @param[in] orient Orientations
   !> @param[in] unit   Output unit
   !> @param[in] num    Atomic numbers (nat)
   !> @param[in] radii  Sphere radii (nat), bohr
   subroutine print_orbit_table(orient, unit, num, radii)
      !> Orientations
      type(sym_orient_type), intent(in) :: orient
      !> Output unit
      integer, intent(in) :: unit
      !> Atomic numbers
      integer, intent(in) :: num(:)
      !> Sphere radii
      real(wp), intent(in) :: radii(:)

      type(prettylistprinter) :: plp
      !> Atoms sorted by orbit, ascending within one, and the orbit offsets
      integer, allocatable :: members(:), start(:)
      !> Distinct rotations of one orbit
      integer, allocatable :: rots(:)
      integer :: nat, iorb, ia, irep, nmem, nrots, k

      nat = size(orient%orbit_of_atom)
      allocate (start(orient%norbit + 1), source=0)
      do ia = 1, nat
         start(orient%orbit_of_atom(ia) + 1) = start(orient%orbit_of_atom(ia) + 1) + 1
      end do
      start(1) = 1
      do iorb = 1, orient%norbit
         start(iorb + 1) = start(iorb + 1) + start(iorb)
      end do
      allocate (members(nat), rots(nat))
      block
         integer, allocatable :: fill(:)
         fill = start(1:orient%norbit)
         do ia = 1, nat
            members(fill(orient%orbit_of_atom(ia))) = ia
            fill(orient%orbit_of_atom(ia)) = fill(orient%orbit_of_atom(ia)) + 1
         end do
      end block

      plp = new_prettylistprinter([6, 4, 10, 5, 5, 10, 20], &
         & [character(len=9) :: "Orbit", "Sym", "R (bohr)", "Mult", "Site", "Rotations", "Atoms"], &
         & unit=unit)
      write (plp%unit, "(a)") "Symmetry-equivalent atoms:"
      call plp%separator()
      call plp%print_header()
      call plp%separator()
      do iorb = 1, orient%norbit
         irep = orient%orbit_rep(iorb)
         nmem = start(iorb + 1) - start(iorb)
         nrots = 0
         do k = start(iorb), start(iorb + 1) - 1
            call insert_sorted(rots, nrots, orient%rot_of_atom(members(k)))
         end do
         call plp%begin_row()
         call plp%add(iorb)
         call plp%add(trim(to_symbol(num(irep))))
         call plp%add(radii(irep), fmt="f10.4")
         call plp%add(nmem)
         call plp%add(count(orient%perm(irep, :) == irep))
         call plp%add(compact_ranges(rots(1:nrots), max_list_len))
         call plp%add(compact_ranges(members(start(iorb):start(iorb + 1) - 1), max_list_len))
         call plp%end_row()
      end do
      call plp%blank()
   end subroutine print_orbit_table

   !> Insert a value into a sorted list unless already present
   !>
   !> @param[in,out] list Sorted values, room for one more
   !> @param[in,out] n    Number of values in `list`
   !> @param[in]     val  Value to insert
   pure subroutine insert_sorted(list, n, val)
      !> Sorted values
      integer, intent(inout) :: list(:)
      !> Number of values
      integer, intent(inout) :: n
      !> Value to insert
      integer, intent(in) :: val

      integer :: pos

      pos = n + 1
      do while (pos > 1)
         if (list(pos - 1) < val) exit
         if (list(pos - 1) == val) return
         pos = pos - 1
      end do
      list(pos + 1:n + 1) = list(pos:n)
      list(pos) = val
      n = n + 1
   end subroutine insert_sorted

   !> Ascending distinct integers as ranges, "1-6", "1,3,5", "2,7-9"
   !>
   !> Cut with ",..." once longer than `max_len`
   !>
   !> @param[in] ids     Ascending distinct integers
   !> @param[in] max_len Longest rendering before the cut
   pure function compact_ranges(ids, max_len) result(text)
      !> Ascending distinct integers
      integer, intent(in) :: ids(:)
      !> Longest rendering before the cut
      integer, intent(in) :: max_len
      !> Rendering
      character(len=:), allocatable :: text

      character(len=:), allocatable :: item
      integer :: i, j

      text = ""
      i = 1
      do while (i <= size(ids))
         j = i
         do while (j < size(ids))
            if (ids(j + 1) /= ids(j) + 1) exit
            j = j + 1
         end do
         ! Runs of three or more as a range, shorter ones listed
         if (j >= i + 2) then
            item = to_string(ids(i))//"-"//to_string(ids(j))
         else
            j = i
            item = to_string(ids(i))
         end if
         if (len(text) > 0) then
            if (len(text) + 1 + len(item) > max_len) then
               text = text//",..."
               return
            end if
            text = text//","
         end if
         text = text//item
         i = j + 1
      end do
   end function compact_ranges

   !> Whether two real arrays agree bit for bit
   !>
   !> @param[in] a First array
   !> @param[in] b Second array
   pure function same_bits(a, b) result(same)
      !> First array
      real(wp), intent(in) :: a(:)
      !> Second array
      real(wp), intent(in) :: b(:)
      !> Whether sizes and bit patterns agree
      logical :: same

      same = size(a) == size(b)
      if (same) same = all(transfer(a, [0_int64]) == transfer(b, [0_int64]))
   end function same_bits

end module moist_cavity_grid_model_symmetry_type
