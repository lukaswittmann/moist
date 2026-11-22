!> Orientations of a Lebedev sphere that keep it invariant under a site group
!>
!> Provides
!>
!> - [[site_candidates]]: admissible orientations closest to the molecular frame
!> - [[neighbour_form_type]], [[smallest_forms]]: tie-break among them
!> - [[operation_axis]], [[oh_closeness]], [[is_oh_element]]: building blocks
!>
!> Concepts
!>
!> - Oh: symmetry group of the cube and of a Lebedev grid, the 48 signed
!>   permutation matrices
!> - site group: operations of the molecule that fix one atom, as matrices s
!> - admissible R: the sphere rotated by R is invariant under the site
!>   group, R^T s R in Oh for every s
!> - one exists only when every proper part det(s) s is the identity or a
!>   rotation by 90, 120 or 180 degrees (no C5, C6, S3 or higher): the site
!>   group is then a subgroup of Oh up to conjugation
!> - closeness of R to the molecular frame F: max over o in Oh of
!>   tr(F^T R o); the grid sees R modulo Oh only
!> - tie: admissible orientations within `tie_tol` of the best closeness
!> - neighbour form: the other atoms in the grid frame of R, reduced modulo
!>   Oh; the smallest form wins a tie, equal forms give congruent cavities
!>
!> Cases by the rotation axes of the site group, and why ties arise:
!> [[site_candidates]]
module moist_math_symmetry_embed
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matdet_3x3
   use moist_math_linalg_geometry, only: cross_product
   implicit none(type, external)
   private

   public :: site_candidates, operation_axis, oh_closeness, is_oh_element
   public :: neighbour_form_type, new_neighbour_form, compare_forms, smallest_forms, sort_keys
   public :: axis_class_none, axis_class_invalid, tie_tol

   !> No rotation axis: identity, inversion
   integer, parameter :: axis_class_none = 1
   !> Proper part not a rotation by 90, 120 or 180 degrees
   integer, parameter :: axis_class_invalid = 0

   !> Tolerance of matrix entries and angles
   real(wp), parameter :: tol = 1.0e-6_wp
   !> Closeness difference within which two orientations tie (closeness is at most 3)
   real(wp), parameter :: tie_tol = 1.0e-9_wp
   !> Coordinate tolerance of neighbour forms, relative to the farthest atom
   real(wp), parameter :: form_tol = 1.0e-7_wp

   !> Canonical form of the other atoms seen from a sphere
   !>
   !> - the vectors d_b = x_b - x_a to every other atom b in the grid frame,
   !>   o R^T d_b, sorted by coordinates (then species), for the o in Oh whose
   !>   sorted list comes first; compared lexicographically by
   !>   [[compare_forms]], coordinates to `tol`
   !> - a function of the molecule and of R modulo Oh: unchanged by rigid
   !>   motion and atom order, given species ids that do not depend on the
   !>   order (atomic numbers, say)
   !> - equal for R and R' exactly when R' = g R o, o in Oh, for an
   !>   orthogonal g mapping every d_b onto a d_b' of its species: a symmetry
   !>   of the molecule fixing atom a, so both spheres give congruent cavities
   !> - sees the odd part of the neighbour set too; a sum of Oh-invariant
   !>   harmonics is even and misses e.g. a planar ring with one axial atom
   type :: neighbour_form_type
      !> Sorted vectors in the grid frame (3, nat - 1)
      real(wp), allocatable :: xyz(:, :)
      !> Their species (nat - 1)
      integer, allocatable :: species(:)
      !> Comparison tolerance of coordinates, `form_tol` times the largest |d_b|
      real(wp) :: tol = 0.0_wp
   end type neighbour_form_type

contains

   !> Admissible sphere orientations of a site group that tie closest to the frame
   !>
   !> Every admissible R whose closeness to F is within `tie_tol` of the
   !> best, one per class modulo Oh, in enumeration order (which depends on
   !> the lab frame: break ties with [[smallest_forms]])
   !>
   !> | site operations span | admissible R                      | candidates               | dw              |
   !> |----------------------|-----------------------------------|--------------------------|-----------------|
   !> | no axis (E, i)       | any                               | R = F                    | identity        |
   !> | one axis line u      | one family, free rotation about u | best angle per target, o | u (dtheta/dw_F) |
   !> | two or more lines    | finitely many                     | all, by enumeration      | zero            |
   !>
   !> Why orientations tie exactly: F is only defined modulo Oh, so every
   !> k in F Oh F^T keeps every closeness, and admissible R and k R tie
   !>
   !> - k (or -k) a symmetry of the molecule fixing the site: congruent
   !>   cavities, either will do
   !> - otherwise the cavities differ and the neighbour form must decide
   !>
   !> | site                         | tied pair                         | k a symmetry |
   !> |------------------------------|-----------------------------------|--------------|
   !> | ammonia N, CH3F C and F, C3v | R, C2(u) R about the C3 axis u    | no           |
   !> | BF3 B, C3v used of D3h       | R, C2(u) R = -sigma_h R           | yes          |
   !> | C3 family, u an F axis       | angles 90 degrees apart about u   | no           |
   !> | Cs family, u normal to F z   | mirror images under C2 about F z  | no           |
   !>
   !> @param[in]  site  Site operations (3, 3, ns), linear parts about the atom
   !> @param[in]  frame Molecular frame F, axes as columns
   !> @param[out] found Whether the site group admits an orientation
   !> @param[out] value Closeness of the best orientation to F
   !> @param[out] rot   Tied orientations (3, 3, ntie), ntie = 0 when not found
   !> @param[out] dw    Their maps delta omega_R = dw delta omega_F (3, 3, ntie)
   subroutine site_candidates(site, frame, found, value, rot, dw)
      !> Site operations
      real(wp), intent(in) :: site(:, :, :)
      !> Molecular frame
      real(wp), intent(in) :: frame(3, 3)
      !> Whether an orientation exists
      logical, intent(out) :: found
      !> Closeness of the best orientation
      real(wp), intent(out) :: value
      !> Tied orientations
      real(wp), allocatable, intent(out) :: rot(:, :, :)
      !> Their derivative maps
      real(wp), allocatable, intent(out) :: dw(:, :, :)

      !> Distinct axis lines and their highest class (4, 3 or 2)
      real(wp) :: lines(3, size(site, 3))
      integer :: classes(size(site, 3))
      real(wp) :: axis(3)
      integer :: nline, is, il, cls
      logical :: known

      found = .false.
      value = -huge(1.0_wp)
      allocate (rot(3, 3, 0), dw(3, 3, 0))

      nline = 0
      do is = 1, size(site, 3)
         call operation_axis(site(:, :, is), axis, cls)
         if (cls == axis_class_invalid) return
         if (cls == axis_class_none) cycle
         known = .false.
         do il = 1, nline
            if (abs(dot_product(lines(:, il), axis)) > 1.0_wp - tol) then
               classes(il) = max(classes(il), cls)
               known = .true.
               exit
            end if
         end do
         if (.not. known) then
            nline = nline + 1
            lines(:, nline) = axis
            classes(nline) = cls
         end if
      end do

      select case (nline)
      case (0)
         found = .true.
         value = 3.0_wp
         rot = reshape(frame, [3, 3, 1])
         dw = reshape([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3, 1])
      case (1)
         call family_candidates(site, frame, lines(:, 1), classes(1), found, value, rot, dw)
      case default
         ! The highest class first: it has the fewest targets
         il = maxloc(classes(1:nline), dim=1)
         is = merge(2, 1, il == 1)
         call finite_candidates(site, frame, lines(:, il), classes(il), lines(:, is), classes(is), &
            & found, value, rot)
         deallocate (dw)
         allocate (dw(3, 3, size(rot, 3)), source=0.0_wp)
      end select
   end subroutine site_candidates

   !> One-parameter family: R = R0 Rot(U, theta), R0 U = u
   !>
   !> - per target U and o in Oh, the best angle in closed form; tr(o F^T R)
   !>   = U^T a U + A cos(theta) + B sin(theta), a = o F^T R0
   !> - A = B = 0 makes the angle free but is never the best: its value is
   !>   at most 1, below the closeness of any orientation
   !> - pairs within `tie_tol` of the best give the candidates; each keeps
   !>   the angle derivative of its own pair
   !>
   !> @param[in]  site  Site operations (3, 3, ns)
   !> @param[in]  frame Molecular frame
   !> @param[in]  u     The axis line
   !> @param[in]  cls   Its class
   !> @param[out] found Whether an orientation exists
   !> @param[out] best  Closeness of the best orientation
   !> @param[out] rot   Tied orientations (3, 3, ntie)
   !> @param[out] dw    Their derivative maps (3, 3, ntie)
   subroutine family_candidates(site, frame, u, cls, found, best, rot, dw)
      !> Site operations
      real(wp), intent(in) :: site(:, :, :)
      !> Molecular frame
      real(wp), intent(in) :: frame(3, 3)
      !> Axis line
      real(wp), intent(in) :: u(3)
      !> Axis class
      integer, intent(in) :: cls
      !> Whether an orientation exists
      logical, intent(out) :: found
      !> Closeness of the best orientation
      real(wp), intent(out) :: best
      !> Tied orientations
      real(wp), allocatable, intent(out) :: rot(:, :, :)
      !> Their derivative maps
      real(wp), allocatable, intent(out) :: dw(:, :, :)

      real(wp) :: targets(3, 2), r0s(3, 3, 2), a(3, 3), da(3, 3), k(3, 3), ek(3, 3), oh(3, 3, 48)
      real(wp) :: big_u(3), ac, b, theta, dac, db, dtheta(3), cand(3, 3)
      !> Value of every (target, o) pair, -huge for invalid targets
      real(wp) :: values(48, 2)
      real(wp) :: rots(3, 3, 96), dws(3, 3, 96)
      integer :: ntarget, it, io, ik, is, ntie, jt
      logical :: valid, known

      found = .false.
      best = -huge(1.0_wp)
      call class_targets(cls, targets, ntarget)
      call oh_elements(oh)
      values = -huge(1.0_wp)

      do it = 1, ntarget
         big_u = targets(:, it)
         r0s(:, :, it) = map_line(big_u, u)
         ! Every member of the family is admissible once one is
         cand = r0s(:, :, it)
         valid = .true.
         do is = 1, size(site, 3)
            if (.not. is_oh_element(matmul(transpose(cand), matmul(site(:, :, is), cand)))) then
               valid = .false.
               exit
            end if
         end do
         if (.not. valid) cycle
         k = skew(big_u)
         do io = 1, 48
            a = matmul(oh(:, :, io), matmul(transpose(frame), r0s(:, :, it)))
            ac = trace3(a) - dot_product(big_u, matmul(a, big_u))
            b = trace3(matmul(a, k))
            values(io, it) = dot_product(big_u, matmul(a, big_u)) + sqrt(ac**2 + b**2)
         end do
         found = .true.
      end do
      if (.not. found) then
         allocate (rot(3, 3, 0), dw(3, 3, 0))
         return
      end if
      best = maxval(values)

      ntie = 0
      do it = 1, ntarget
         big_u = targets(:, it)
         k = skew(big_u)
         do io = 1, 48
            if (values(io, it) < best - tie_tol) cycle
            a = matmul(oh(:, :, io), matmul(transpose(frame), r0s(:, :, it)))
            ac = trace3(a) - dot_product(big_u, matmul(a, big_u))
            b = trace3(matmul(a, k))
            theta = atan2(b, ac)
            cand = matmul(r0s(:, :, it), rotation(big_u, theta))
            known = .false.
            do jt = 1, ntie
               known = is_oh_element(matmul(transpose(rots(:, :, jt)), cand))
               if (known) exit
            end do
            if (known) cycle
            ntie = ntie + 1
            rots(:, :, ntie) = cand

            ! a = o F^T R0, so delta a = -o F^T [delta omega]_x R0 = -a R0^T [dw]_x R0
            dws(:, :, ntie) = 0.0_wp
            if (ac**2 + b**2 > tiny(1.0_wp)) then
               do ik = 1, 3
                  ek = skew(unit3(ik))
                  da = -matmul(a, matmul(transpose(r0s(:, :, it)), matmul(ek, r0s(:, :, it))))
                  dac = trace3(da) - dot_product(big_u, matmul(da, big_u))
                  db = trace3(matmul(da, k))
                  dtheta(ik) = (ac*db - b*dac)/(ac**2 + b**2)
               end do
               do ik = 1, 3
                  dws(:, ik, ntie) = u*dtheta(ik)
               end do
            end if
         end do
      end do
      rot = rots(:, :, 1:ntie)
      dw = dws(:, :, 1:ntie)
   end subroutine family_candidates

   !> Finitely many orientations: map (u1, u2) onto pairs of Oh lines
   !>
   !> @param[in]  site  Site operations (3, 3, ns)
   !> @param[in]  frame Molecular frame
   !> @param[in]  u1    First axis line
   !> @param[in]  c1    Its class
   !> @param[in]  u2    Second axis line, not parallel to u1
   !> @param[in]  c2    Its class
   !> @param[out] found Whether an orientation exists
   !> @param[out] best  Closeness of the best orientation
   !> @param[out] rot   Tied orientations (3, 3, ntie), distinct modulo Oh
   subroutine finite_candidates(site, frame, u1, c1, u2, c2, found, best, rot)
      !> Site operations
      real(wp), intent(in) :: site(:, :, :)
      !> Molecular frame
      real(wp), intent(in) :: frame(3, 3)
      !> First axis line
      real(wp), intent(in) :: u1(3)
      !> Its class
      integer, intent(in) :: c1
      !> Second axis line
      real(wp), intent(in) :: u2(3)
      !> Its class
      integer, intent(in) :: c2
      !> Whether an orientation exists
      logical, intent(out) :: found
      !> Closeness of the best orientation
      real(wp), intent(out) :: best
      !> Tied orientations
      real(wp), allocatable, intent(out) :: rot(:, :, :)

      !> Every admissible orientation (at most 2 x 9 x 4) and its closeness
      real(wp) :: cands(3, 3, 72), values(72)
      real(wp) :: t1(3, 2), all2(3, 13), cand(3, 3), cosa, v2(3), s1
      integer :: n1, n2, i1, i2, is, isign, jsign, ncand, ic, jc, ntie
      logical :: valid, known
      logical :: keep(72)

      found = .false.
      best = -huge(1.0_wp)
      ncand = 0
      call class_targets(c1, t1, n1)
      call class_lines(c2, all2, n2)
      cosa = dot_product(u1, u2)

      ! Lines carry no sign: R takes U1 onto +-u1 and V onto u2 at the same angle
      do i1 = 1, n1
         do i2 = 1, n2
            if (abs(abs(dot_product(t1(:, i1), all2(:, i2))) - 1.0_wp) < tol) cycle
            do isign = 0, 1
               s1 = merge(-1.0_wp, 1.0_wp, isign == 1)
               do jsign = 0, 1
                  v2 = all2(:, i2)*merge(-1.0_wp, 1.0_wp, jsign == 1)
                  if (abs(dot_product(t1(:, i1), v2) - s1*cosa) > tol) cycle
                  cand = map_pair(t1(:, i1), v2, s1*u1, u2)
                  valid = .true.
                  do is = 1, size(site, 3)
                     if (.not. is_oh_element(matmul(transpose(cand), matmul(site(:, :, is), cand)))) then
                        valid = .false.
                        exit
                     end if
                  end do
                  if (.not. valid) cycle
                  ncand = ncand + 1
                  cands(:, :, ncand) = cand
                  values(ncand) = oh_closeness(frame, cand)
               end do
            end do
         end do
      end do
      found = ncand > 0
      if (found) best = maxval(values(1:ncand))

      ! Tied orientations, one per class modulo Oh
      ntie = 0
      keep = .false.
      do ic = 1, ncand
         if (values(ic) < best - tie_tol) cycle
         known = .false.
         do jc = 1, ic - 1
            known = keep(jc) .and. is_oh_element(matmul(transpose(cands(:, :, jc)), cands(:, :, ic)))
            if (known) exit
         end do
         keep(ic) = .not. known
         if (keep(ic)) ntie = ntie + 1
      end do
      allocate (rot(3, 3, ntie))
      ntie = 0
      do ic = 1, ncand
         if (.not. keep(ic)) cycle
         ntie = ntie + 1
         rot(:, :, ntie) = cands(:, :, ic)
      end do
   end subroutine finite_candidates

   !> The orientations with the smallest neighbour form at a site
   !>
   !> Equal forms give congruent cavities (a symmetry of the molecule fixing
   !> the site relates them)
   !>
   !> @param[in] xyz     Coordinates (3, nat)
   !> @param[in] species Species ids (nat), labels independent of the atom order
   !> @param[in] iat     Site atom
   !> @param[in] rots    Orientations (3, 3, n), n >= 1
   function smallest_forms(xyz, species, iat, rots) result(smallest)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Site atom
      integer, intent(in) :: iat
      !> Orientations
      real(wp), intent(in) :: rots(:, :, :)
      !> Whether each orientation has the smallest form (n)
      logical, allocatable :: smallest(:)

      type(neighbour_form_type), allocatable :: forms(:)
      integer :: ir, ibest

      allocate (forms(size(rots, 3)), smallest(size(rots, 3)))
      ibest = 1
      do ir = 1, size(rots, 3)
         call new_neighbour_form(forms(ir), xyz, species, iat, rots(:, :, ir))
         if (compare_forms(forms(ir), forms(ibest)) < 0) ibest = ir
      end do
      do ir = 1, size(rots, 3)
         smallest(ir) = compare_forms(forms(ir), forms(ibest)) == 0
      end do
   end function smallest_forms

   !> Neighbour form of a sphere: the other atoms in its grid frame, modulo Oh
   !>
   !> @param[out] self    Form
   !> @param[in]  xyz     Coordinates (3, nat)
   !> @param[in]  species Species ids (nat), labels independent of the atom order
   !> @param[in]  iat     Atom of the sphere
   !> @param[in]  rot     Sphere orientation R
   pure subroutine new_neighbour_form(self, xyz, species, iat, rot)
      !> Form
      type(neighbour_form_type), intent(out) :: self
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Atom of the sphere
      integer, intent(in) :: iat
      !> Sphere orientation
      real(wp), intent(in) :: rot(3, 3)

      !> Vectors to the other atoms in the grid frame, and their species
      real(wp), allocatable :: d(:, :)
      integer, allocatable :: sp(:), idx(:)
      type(neighbour_form_type) :: trial
      real(wp) :: oh(3, 3, 48)
      integer :: nat, ib, n, io

      nat = size(xyz, 2)
      allocate (d(3, max(nat - 1, 0)), sp(max(nat - 1, 0)), idx(max(nat - 1, 0)))
      n = 0
      do ib = 1, nat
         if (ib == iat) cycle
         n = n + 1
         d(:, n) = matmul(transpose(rot), xyz(:, ib) - xyz(:, iat))
         sp(n) = species(ib)
      end do
      self%tol = 0.0_wp
      if (n > 0) self%tol = form_tol*maxval(norm2(d, dim=1))
      trial%tol = self%tol

      call oh_elements(oh)
      do io = 1, 48
         trial%xyz = matmul(oh(:, :, io), d)
         call sort_columns(trial%xyz, sp, self%tol, idx)
         trial%xyz = trial%xyz(:, idx)
         trial%species = sp(idx)
         if (io == 1) then
            self%xyz = trial%xyz
            self%species = trial%species
         else if (compare_forms(trial, self) < 0) then
            self%xyz = trial%xyz
            self%species = trial%species
         end if
      end do
   end subroutine new_neighbour_form

   !> Lexicographic order of two neighbour forms
   !>
   !> Coordinates column by column to the larger tolerance, then species
   !>
   !> @param[in] a First form
   !> @param[in] b Second form
   pure function compare_forms(a, b) result(order)
      !> First form
      type(neighbour_form_type), intent(in) :: a
      !> Second form
      type(neighbour_form_type), intent(in) :: b
      !> -1 when a comes first, 1 when b does, 0 when equal
      integer :: order

      real(wp) :: eps
      integer :: i, j, n

      order = 0
      eps = max(a%tol, b%tol)
      n = min(size(a%species), size(b%species))
      do j = 1, n
         do i = 1, 3
            if (a%xyz(i, j) < b%xyz(i, j) - eps) then
               order = -1
               return
            else if (a%xyz(i, j) > b%xyz(i, j) + eps) then
               order = 1
               return
            end if
         end do
      end do
      do j = 1, n
         if (a%species(j) /= b%species(j)) then
            order = merge(-1, 1, a%species(j) < b%species(j))
            return
         end if
      end do
      if (size(a%species) /= size(b%species)) order = merge(-1, 1, size(a%species) < size(b%species))
   end function compare_forms

   !> Stable merge sort of columns by coordinates (to a tolerance), then species
   !>
   !> - every coordinate is first reduced to the index of its cluster of
   !>   values closer than the tolerance (single linkage), so the order is
   !>   transitive and does not depend on the order of the columns
   !>
   !> @param[in]  x   Columns (3, n)
   !> @param[in]  sp  Species (n)
   !> @param[in]  eps Coordinate tolerance
   !> @param[out] idx Sorting permutation (n)
   pure subroutine sort_columns(x, sp, eps, idx)
      !> Columns
      real(wp), intent(in) :: x(:, :)
      !> Species
      integer, intent(in) :: sp(:)
      !> Coordinate tolerance
      real(wp), intent(in) :: eps
      !> Sorting permutation
      integer, intent(out) :: idx(:)

      integer, allocatable :: key(:, :), order(:)
      integer :: n, i, k

      n = size(sp)
      allocate (key(4, n), order(n))
      do k = 1, 3
         ! Exact order of one coordinate, then its clusters
         key(1, :) = 0
         call sort_keys(key(1:1, :), order, x(k, :))
         do i = 1, n
            if (i == 1) then
               key(k + 1, order(i)) = 1
            else if (x(k, order(i)) - x(k, order(i - 1)) > eps) then
               key(k + 1, order(i)) = key(k + 1, order(i - 1)) + 1
            else
               key(k + 1, order(i)) = key(k + 1, order(i - 1))
            end if
         end do
      end do
      key(1:3, :) = key(2:4, :)
      key(4, :) = sp
      call sort_keys(key, idx)
   end subroutine sort_columns

   !> Stable merge sort by integer keys, or by a real value
   !>
   !> @param[in]  key   Keys (nkey, n), compared in turn
   !> @param[out] idx   Sorting permutation (n)
   !> @param[in]  value Sort by this value (n) in place of the keys
   pure subroutine sort_keys(key, idx, value)
      !> Keys
      integer, intent(in) :: key(:, :)
      !> Sorting permutation
      integer, intent(out) :: idx(:)
      !> Value to sort by
      real(wp), intent(in), optional :: value(:)

      integer, allocatable :: work(:)
      integer :: n, width, lo, mid, hi, i, j, k
      logical :: before

      n = size(idx)
      idx = [(i, i=1, n)]
      allocate (work(n))
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
               if (present(value)) then
                  before = value(idx(j)) < value(idx(i))
               else
                  before = key_order(key(:, idx(j)), key(:, idx(i))) < 0
               end if
               if (before) then
                  work(k) = idx(j)
                  j = j + 1
               else
                  work(k) = idx(i)
                  i = i + 1
               end if
               k = k + 1
            end do
            do while (i <= mid)
               work(k) = idx(i)
               i = i + 1
               k = k + 1
            end do
            do while (j <= hi)
               work(k) = idx(j)
               j = j + 1
               k = k + 1
            end do
            lo = lo + 2*width
         end do
         idx = work
         width = 2*width
      end do
   end subroutine sort_keys

   !> Lexicographic order of two integer keys
   !>
   !> @param[in] a First key
   !> @param[in] b Second key
   pure function key_order(a, b) result(order)
      !> First key
      integer, intent(in) :: a(:)
      !> Second key
      integer, intent(in) :: b(:)
      !> -1 when a comes first, 1 when b does, 0 when equal
      integer :: order

      integer :: k

      order = 0
      do k = 1, size(a)
         if (a(k) /= b(k)) then
            order = merge(-1, 1, a(k) < b(k))
            return
         end if
      end do
   end function key_order

   !> Axis line and class of an operation
   !>
   !> Class of the proper part det(s) s: `axis_class_none` for the identity,
   !> 4, 3 or 2 for a rotation by 90, 120 or 180 degrees, `axis_class_invalid`
   !> otherwise
   !>
   !> @param[in]  s    Orthogonal matrix
   !> @param[out] axis Rotation axis of the proper part, unit (sign arbitrary)
   !> @param[out] cls  Class
   pure subroutine operation_axis(s, axis, cls)
      !> Orthogonal matrix
      real(wp), intent(in) :: s(3, 3)
      !> Rotation axis
      real(wp), intent(out) :: axis(3)
      !> Class
      integer, intent(out) :: cls

      real(wp) :: sp(3, 3), c, theta
      integer :: k

      sp = s*matdet_3x3(s)
      c = min(1.0_wp, max(-1.0_wp, (trace3(sp) - 1.0_wp)/2.0_wp))
      axis = [0.0_wp, 0.0_wp, 1.0_wp]
      if (c > 1.0_wp - tol) then
         cls = axis_class_none
         return
      end if
      theta = acos(c)
      if (abs(theta - pi/2.0_wp) < tol) then
         cls = 4
      else if (abs(theta - 2.0_wp*pi/3.0_wp) < tol) then
         cls = 3
      else if (abs(theta - pi) < tol) then
         cls = 2
      else
         cls = axis_class_invalid
      end if

      if (theta < pi - tol) then
         axis = [sp(3, 2) - sp(2, 3), sp(1, 3) - sp(3, 1), sp(2, 1) - sp(1, 2)]
      else
         ! (sp + I)/2 = u u^T for a half turn
         k = maxloc([sp(1, 1), sp(2, 2), sp(3, 3)], dim=1)
         axis = (sp(:, k) + unit3(k))/2.0_wp
      end if
      axis = axis/norm2(axis)
   end subroutine operation_axis

   !> Closeness of R to F modulo Oh: max over o in Oh of tr(F^T R o)
   !>
   !> Equals the largest sum of |(F^T R)_{i, pi(i)}| over permutations pi
   !>
   !> @param[in] frame Molecular frame F
   !> @param[in] rot   Orientation R
   pure function oh_closeness(frame, rot) result(value)
      !> Molecular frame
      real(wp), intent(in) :: frame(3, 3)
      !> Orientation
      real(wp), intent(in) :: rot(3, 3)
      !> Closeness, at most 3
      real(wp) :: value

      integer, parameter :: perms(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, &
         & 3, 2, 1], [3, 6])
      real(wp) :: b(3, 3)
      integer :: ip

      b = abs(matmul(transpose(frame), rot))
      value = -huge(1.0_wp)
      do ip = 1, 6
         value = max(value, b(1, perms(1, ip)) + b(2, perms(2, ip)) + b(3, perms(3, ip)))
      end do
   end function oh_closeness

   !> Whether a matrix is a signed permutation matrix
   !>
   !> @param[in] m Matrix
   pure function is_oh_element(m) result(yes)
      !> Matrix
      real(wp), intent(in) :: m(3, 3)
      !> Whether m is in Oh
      logical :: yes

      yes = all(abs(m) < tol .or. abs(abs(m) - 1.0_wp) < tol) .and. &
         & all(count(abs(m) > 0.5_wp, dim=1) == 1) .and. all(count(abs(m) > 0.5_wp, dim=2) == 1)
   end function is_oh_element

   !> The 48 elements of Oh, signed permutation matrices
   !>
   !> @param[out] oh Elements (3, 3, 48)
   pure subroutine oh_elements(oh)
      !> Elements
      real(wp), intent(out) :: oh(3, 3, 48)

      integer, parameter :: perms(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, &
         & 3, 2, 1], [3, 6])
      integer :: ip, is, i, n
      real(wp) :: sgn(3)

      n = 0
      oh = 0.0_wp
      do ip = 1, 6
         do is = 0, 7
            sgn = [merge(-1.0_wp, 1.0_wp, btest(is, 0)), merge(-1.0_wp, 1.0_wp, btest(is, 1)), &
               & merge(-1.0_wp, 1.0_wp, btest(is, 2))]
            n = n + 1
            do i = 1, 3
               oh(perms(i, ip), i, n) = sgn(i)
            end do
         end do
      end do
   end subroutine oh_elements

   !> One representative Oh line per Oh orbit admitting a class
   !>
   !> @param[in]  cls     Class 4, 3 or 2
   !> @param[out] targets Representative lines (3, 2)
   !> @param[out] n       Number of representatives
   pure subroutine class_targets(cls, targets, n)
      !> Class
      integer, intent(in) :: cls
      !> Representative lines
      real(wp), intent(out) :: targets(3, 2)
      !> Number of representatives
      integer, intent(out) :: n

      targets = 0.0_wp
      select case (cls)
      case (4)
         n = 1
         targets(:, 1) = [0.0_wp, 0.0_wp, 1.0_wp]
      case (3)
         n = 1
         targets(:, 1) = [1.0_wp, 1.0_wp, 1.0_wp]/sqrt(3.0_wp)
      case default
         n = 2
         targets(:, 1) = [0.0_wp, 0.0_wp, 1.0_wp]
         targets(:, 2) = [1.0_wp, 1.0_wp, 0.0_wp]/sqrt(2.0_wp)
      end select
   end subroutine class_targets

   !> Every Oh line admitting a class
   !>
   !> @param[in]  cls   Class 4, 3 or 2
   !> @param[out] lines Lines (3, 13)
   !> @param[out] n     Number of lines
   pure subroutine class_lines(cls, lines, n)
      !> Class
      integer, intent(in) :: cls
      !> Lines
      real(wp), intent(out) :: lines(3, 13)
      !> Number of lines
      integer, intent(out) :: n

      real(wp) :: r2, r3

      r2 = 1.0_wp/sqrt(2.0_wp)
      r3 = 1.0_wp/sqrt(3.0_wp)
      lines = 0.0_wp
      n = 0
      if (cls == 4 .or. cls == 2) then
         lines(:, 1) = [1.0_wp, 0.0_wp, 0.0_wp]
         lines(:, 2) = [0.0_wp, 1.0_wp, 0.0_wp]
         lines(:, 3) = [0.0_wp, 0.0_wp, 1.0_wp]
         n = 3
      end if
      if (cls == 2) then
         lines(:, 4) = [r2, r2, 0.0_wp]
         lines(:, 5) = [r2, -r2, 0.0_wp]
         lines(:, 6) = [r2, 0.0_wp, r2]
         lines(:, 7) = [r2, 0.0_wp, -r2]
         lines(:, 8) = [0.0_wp, r2, r2]
         lines(:, 9) = [0.0_wp, r2, -r2]
         n = 9
      end if
      if (cls == 3) then
         lines(:, 1) = [r3, r3, r3]
         lines(:, 2) = [r3, r3, -r3]
         lines(:, 3) = [r3, -r3, r3]
         lines(:, 4) = [-r3, r3, r3]
         n = 4
      end if
   end subroutine class_lines

   !> A rotation taking the unit vector a onto the unit vector b
   !>
   !> @param[in] a Source
   !> @param[in] b Target
   pure function map_line(a, b) result(r)
      !> Source
      real(wp), intent(in) :: a(3)
      !> Target
      real(wp), intent(in) :: b(3)
      !> Rotation with r a = b
      real(wp) :: r(3, 3)

      r = map_pair(a, perpendicular(a), b, perpendicular(b))
   end function map_line

   !> The rotation taking the orthonormalized pair (a1, a2) onto (b1, b2)
   !>
   !> @param[in] a1 First source direction, unit
   !> @param[in] a2 Second source direction, not parallel to a1
   !> @param[in] b1 First target direction, unit
   !> @param[in] b2 Second target direction, at the angle of (a1, a2) to b1
   pure function map_pair(a1, a2, b1, b2) result(r)
      !> First source direction
      real(wp), intent(in) :: a1(3)
      !> Second source direction
      real(wp), intent(in) :: a2(3)
      !> First target direction
      real(wp), intent(in) :: b1(3)
      !> Second target direction
      real(wp), intent(in) :: b2(3)
      !> Rotation
      real(wp) :: r(3, 3)

      real(wp) :: fa(3, 3), fb(3, 3)

      fa(:, 1) = a1
      fa(:, 2) = a2 - a1*dot_product(a1, a2)
      fa(:, 2) = fa(:, 2)/norm2(fa(:, 2))
      fa(:, 3) = cross_product(fa(:, 1), fa(:, 2))
      fb(:, 1) = b1
      fb(:, 2) = b2 - b1*dot_product(b1, b2)
      fb(:, 2) = fb(:, 2)/norm2(fb(:, 2))
      fb(:, 3) = cross_product(fb(:, 1), fb(:, 2))
      r = matmul(fb, transpose(fa))
   end function map_pair

   !> A unit vector normal to the unit vector a
   !>
   !> @param[in] a Unit vector
   pure function perpendicular(a) result(p)
      !> Unit vector
      real(wp), intent(in) :: a(3)
      !> Normal unit vector
      real(wp) :: p(3)

      integer :: k

      k = minloc(abs(a), dim=1)
      p = 0.0_wp
      p(k) = 1.0_wp
      p = p - a*dot_product(a, p)
      p = p/norm2(p)
   end function perpendicular

   !> Rotation by theta about the unit vector u
   !>
   !> @param[in] u     Axis, unit
   !> @param[in] theta Angle
   pure function rotation(u, theta) result(r)
      !> Axis
      real(wp), intent(in) :: u(3)
      !> Angle
      real(wp), intent(in) :: theta
      !> Rotation matrix
      real(wp) :: r(3, 3)

      integer :: i

      r = sin(theta)*skew(u)
      do i = 1, 3
         r(:, i) = r(:, i) + (1.0_wp - cos(theta))*u*u(i)
         r(i, i) = r(i, i) + cos(theta)
      end do
   end function rotation

   !> Cross-product matrix, skew(u) v = u x v
   !>
   !> @param[in] u Vector
   pure function skew(u) result(k)
      !> Vector
      real(wp), intent(in) :: u(3)
      !> Skew-symmetric matrix
      real(wp) :: k(3, 3)

      k = 0.0_wp
      k(2, 1) = u(3)
      k(1, 2) = -u(3)
      k(3, 1) = -u(2)
      k(1, 3) = u(2)
      k(3, 2) = u(1)
      k(2, 3) = -u(1)
   end function skew

   !> Unit vector along a Cartesian axis
   !>
   !> @param[in] k Axis index
   pure function unit3(k) result(e)
      !> Axis index
      integer, intent(in) :: k
      !> Unit vector
      real(wp) :: e(3)

      e = 0.0_wp
      e(k) = 1.0_wp
   end function unit3

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

end module moist_math_symmetry_embed
