!> Point-group detection through libmsym
!>
!> Provides
!>
!> - [[find_point_group]]: candidate group of a set of atoms, with thresholds
!>   and retries chosen for the caller's tolerance
!> - [[msym_group_type]], [[msym_subgroup_type]]: operations and subgroups
!> - `msym_type_*`, `msym_op_*`: libmsym codes the caller needs
!>
!> Conventions
!>
!> - libmsym is vendored under `src/moist/math/symmetry/libmsym/` (MIT; see
!>   its `PROVENANCE.md`) and reached through the C shim `msym_shim.c`
!> - atoms are told apart by the caller's species id alone
!> - operations act on coordinates relative to `center`, the plain centroid
!> - linear groups (n = 0) carry the molecular axis and no subgroup list;
!>   their infinite rotation is not an operation matrix the caller can use
!> - libmsym keeps its error details in a static buffer: the calls run in
!>   an OpenMP critical section
module moist_math_symmetry_msym
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_ptr, &
      & c_null_char, c_associated
   use mctc_env, only: error_type, fatal_error, wp
   implicit none(type, external)
   private

   public :: msym_group_type, msym_subgroup_type, find_point_group
   public :: msym_type_cnv, msym_type_dnh, msym_op_identity, msym_op_proper, &
      & msym_op_improper, msym_op_reflection, msym_op_inversion

   !> libmsym point group type of C_nv (and C_inf_v for n = 0)
   integer, parameter :: msym_type_cnv = 6
   !> libmsym point group type of D_nh (and D_inf_h for n = 0)
   integer, parameter :: msym_type_dnh = 8

   !> libmsym operation type of the identity
   integer, parameter :: msym_op_identity = 0
   !> libmsym operation type of a proper rotation
   integer, parameter :: msym_op_proper = 1
   !> libmsym operation type of an improper rotation
   integer, parameter :: msym_op_improper = 2
   !> libmsym operation type of a reflection
   integer, parameter :: msym_op_reflection = 3
   !> libmsym operation type of the inversion
   integer, parameter :: msym_op_inversion = 4

   !> Length of the error message buffer
   integer, parameter :: msg_len = 256

   !> libmsym's default equivalence threshold, relative
   !>
   !> - merges unrelated orbits of a large molecule into one set, which
   !>   libmsym then rejects or reads as a lower group: the threshold used
   !>   follows the caller's tolerance instead
   real(wp), parameter :: equivalence_default = 5.0e-4_wp
   !> Smallest equivalence threshold used; below it rounding splits orbits
   real(wp), parameter :: equivalence_floor = 1.0e-10_wp
   !> Equivalence threshold per tolerance in bohr
   real(wp), parameter :: equivalence_per_tolerance = 10.0_wp
   !> Tighter equivalence thresholds tried in turn when libmsym fails
   real(wp), parameter :: equivalence_ladder(3) = [5.0e-5_wp, 1.0e-7_wp, 1.0e-9_wp]
   !> Angle threshold of libmsym's defaults, as the shim reads "not positive"
   real(wp), parameter :: angle_default = 0.0_wp
   !> Angle threshold of the last retry, on 1 - |cos|
   !>
   !> - the default 1e-3 calls axes up to 0.045 rad apart parallel, wider than
   !>   the permutation test (5e-3) accepts for the group it then builds: the
   !>   axes through the atoms of a nearly linear molecule merge into one
   !>   C_inf axis that is then rejected
   !> - 1e-6 keeps parallel axes within 1.4e-3 rad
   real(wp), parameter :: angle_tight = 1.0e-6_wp
   !> Most atoms for the second detection at the default threshold
   integer, parameter :: hint_max_atoms = 1000
   !> Offset ranking a linear group above every finite one
   integer, parameter :: linear_rank = huge(1)/2

   !> Subgroup of a detected point group
   type :: msym_subgroup_type
      !> Schoenflies symbol, e.g. "C2v"
      character(len=:), allocatable :: name
      !> Indices of its operations into the parent's operation list (order)
      integer, allocatable :: ops(:)
   end type msym_subgroup_type

   !> Detected point group with its operations about `center`
   type :: msym_group_type
      !> Schoenflies symbol as libmsym spells it, e.g. "D6h", "D0h" for D_inf_h
      character(len=:), allocatable :: name
      !> libmsym point group type
      integer :: kind = 0
      !> Principal order, 0 for linear groups
      integer :: n = 0
      !> Whether the group is C_inf_v or D_inf_h
      logical :: linear = .false.
      !> Operation matrices (3, 3, order), acting on coordinates relative to `center`
      real(wp), allocatable :: ops(:, :, :)
      !> Operation type, order, power and orientation (4, order), libmsym codes
      integer, allocatable :: info(:, :)
      !> Axis or plane normal of each operation (3, order), zero for E and i
      real(wp), allocatable :: axis(:, :)
      !> Fixed point of every operation, the centroid
      real(wp) :: center(3) = 0.0_wp
      !> Subgroups, empty for linear groups
      type(msym_subgroup_type), allocatable :: subgroups(:)
   end type msym_group_type

   interface

      !> Create a libmsym context
      !>
      !> @param[out] handle Opaque context, released with `moist_msym_release`
      function moist_msym_create(handle) result(status) bind(C, name="moist_msym_create")
         import :: c_int, c_ptr
         implicit none(type, external)
         !> Opaque context
         type(c_ptr), intent(out) :: handle
         !> Zero on success
         integer(c_int) :: status
      end function moist_msym_create

      !> Release a libmsym context
      !>
      !> @param[in] handle Opaque context
      subroutine moist_msym_release(handle) bind(C, name="moist_msym_release")
         import :: c_ptr
         implicit none(type, external)
         !> Opaque context
         type(c_ptr), value, intent(in) :: handle
      end subroutine moist_msym_release

      !> Set the atoms and detect their point group
      !>
      !> @param[in] handle  Opaque context
      !> @param[in] nat     Number of atoms
      !> @param[in] xyz     Coordinates (3, nat)
      !> @param[in] species Species id per atom, positive
      !> @param[in] equivalence Equivalence threshold; not positive keeps libmsym's
      !> @param[in] angle   Angle threshold; not positive keeps libmsym's
      function moist_msym_find(handle, nat, xyz, species, equivalence, angle) result(status) &
         & bind(C, name="moist_msym_find")
         import :: c_int, c_ptr, c_double
         implicit none(type, external)
         !> Opaque context
         type(c_ptr), value, intent(in) :: handle
         !> Number of atoms
         integer(c_int), value, intent(in) :: nat
         !> Coordinates
         real(c_double), intent(in) :: xyz(*)
         !> Species ids
         integer(c_int), intent(in) :: species(*)
         !> Equivalence threshold
         real(c_double), value, intent(in) :: equivalence
         !> Angle threshold
         real(c_double), value, intent(in) :: angle
         !> Zero on success
         integer(c_int) :: status
      end function moist_msym_find

      !> Name, type, principal order and order of the detected group
      !>
      !> @param[in]  handle Opaque context
      !> @param[out] name   Schoenflies symbol, NUL padded
      !> @param[out] kind   libmsym point group type
      !> @param[out] n      Principal order, 0 for linear groups
      !> @param[out] order  Number of operations
      function moist_msym_group(handle, name, kind, n, order) result(status) &
         & bind(C, name="moist_msym_group")
         import :: c_int, c_ptr, c_char
         implicit none(type, external)
         !> Opaque context
         type(c_ptr), value, intent(in) :: handle
         !> Schoenflies symbol
         character(kind=c_char), intent(out) :: name(8)
         !> Point group type
         integer(c_int), intent(out) :: kind
         !> Principal order
         integer(c_int), intent(out) :: n
         !> Number of operations
         integer(c_int), intent(out) :: order
         !> Zero on success
         integer(c_int) :: status
      end function moist_msym_group

      !> Operation matrices, codes and the center they act about
      !>
      !> @param[in]  handle Opaque context
      !> @param[in]  order  Number of operations
      !> @param[out] mats   Matrices (3, 3, order), column-major
      !> @param[out] info   Type, order, power, orientation (4, order)
      !> @param[out] center Centroid
      function moist_msym_operations(handle, order, mats, info, center) result(status) &
         & bind(C, name="moist_msym_operations")
         import :: c_int, c_ptr, c_double
         implicit none(type, external)
         !> Opaque context
         type(c_ptr), value, intent(in) :: handle
         !> Number of operations
         integer(c_int), value, intent(in) :: order
         !> Operation matrices
         real(c_double), intent(out) :: mats(*)
         !> Operation codes
         integer(c_int), intent(out) :: info(*)
         !> Centroid
         real(c_double), intent(out) :: center(3)
         !> Zero on success
         integer(c_int) :: status
      end function moist_msym_operations

      !> Axis or plane normal of one operation
      !>
      !> @param[in]  handle Opaque context
      !> @param[in]  iop    Operation index, 0-based
      !> @param[out] axis   Axis or normal
      function moist_msym_axis(handle, iop, axis) result(status) bind(C, name="moist_msym_axis")
         import :: c_int, c_ptr, c_double
         implicit none(type, external)
         !> Opaque context
         type(c_ptr), value, intent(in) :: handle
         !> Operation index
         integer(c_int), value, intent(in) :: iop
         !> Axis or normal
         real(c_double), intent(out) :: axis(3)
         !> Zero on success
         integer(c_int) :: status
      end function moist_msym_axis

      !> Number of subgroups of the detected group
      !>
      !> @param[in]  handle Opaque context
      !> @param[out] nsub   Number of subgroups
      function moist_msym_subgroup_count(handle, nsub) result(status) &
         & bind(C, name="moist_msym_subgroup_count")
         import :: c_int, c_ptr
         implicit none(type, external)
         !> Opaque context
         type(c_ptr), value, intent(in) :: handle
         !> Number of subgroups
         integer(c_int), intent(out) :: nsub
         !> Zero on success
         integer(c_int) :: status
      end function moist_msym_subgroup_count

      !> One subgroup as indices into the group's operation list
      !>
      !> @param[in]  handle   Opaque context
      !> @param[in]  isub     Subgroup index, 0-based
      !> @param[out] name     Schoenflies symbol, NUL padded
      !> @param[out] order    Number of operations
      !> @param[in]  maxorder Capacity of `ops`
      !> @param[out] ops      Operation indices, 0-based
      function moist_msym_subgroup(handle, isub, name, order, maxorder, ops) result(status) &
         & bind(C, name="moist_msym_subgroup")
         import :: c_int, c_ptr, c_char
         implicit none(type, external)
         !> Opaque context
         type(c_ptr), value, intent(in) :: handle
         !> Subgroup index
         integer(c_int), value, intent(in) :: isub
         !> Schoenflies symbol
         character(kind=c_char), intent(out) :: name(8)
         !> Number of operations
         integer(c_int), intent(out) :: order
         !> Capacity of ops
         integer(c_int), value, intent(in) :: maxorder
         !> Operation indices
         integer(c_int), intent(out) :: ops(*)
         !> Zero on success
         integer(c_int) :: status
      end function moist_msym_subgroup

      !> Render a status and libmsym's last error details
      !>
      !> @param[in]  status Status of a failed call
      !> @param[out] buf    NUL-terminated message
      !> @param[in]  len    Capacity of buf
      subroutine moist_msym_error_message(status, buf, len) &
         & bind(C, name="moist_msym_error_message")
         import :: c_int, c_char
         implicit none(type, external)
         !> Status of a failed call
         integer(c_int), value, intent(in) :: status
         !> Message
         character(kind=c_char), intent(out) :: buf(*)
         !> Capacity of buf
         integer(c_int), value, intent(in) :: len
      end subroutine moist_msym_error_message

   end interface

contains

   !> Detect the point group of a set of atoms
   !>
   !> - `species` decides which atoms are interchangeable; it must be positive
   !> - a single atom or an empty set is an error: the caller has no frame to
   !>   take from a point
   !> - the group is a candidate: the caller checks it at its own tolerance
   !> - without `tolerance`: libmsym's default thresholds; on failure once
   !>   more at `angle_tight`
   !>
   !> With `tolerance`, in this order:
   !>
   !>  1. equivalence threshold min(default, max(floor, 10 tolerance)); a
   !>     failure of libmsym is retried at the tighter thresholds of
   !>     `equivalence_ladder`
   !>  2. up to `hint_max_atoms` atoms, a second detection at the default
   !>     equivalence threshold unless step 1 used it; its group is taken
   !>     when step 1 failed or found a smaller one (`group_size`), so a
   !>     geometry symmetric only beyond the tolerance keeps its group as a
   !>     candidate, which the caller reports with its deviation
   !>  3. when all of that failed, step 1 once more at `angle_tight`: a nearly
   !>     linear molecule that libmsym reads as linear and then rejects
   !>     ("Invalid equivalence set") gets the group of its bent shape
   !>
   !> @param[in]  xyz       Coordinates (3, nat), any length unit
   !> @param[in]  species   Species id per atom (nat), positive
   !> @param[out] group     Detected group, operations and subgroups
   !> @param[out] error     Set when libmsym fails
   !> @param[in]  tolerance Exactness the caller checks the group at, bohr
   subroutine find_point_group(xyz, species, group, error, tolerance)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Detected group
      type(msym_group_type), intent(out) :: group
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Exactness tolerance of the caller
      real(wp), intent(in), optional :: tolerance

      !> Detection at the default threshold
      type(msym_group_type) :: loose
      type(error_type), allocatable :: loose_error
      !> Equivalence threshold following the tolerance
      real(wp) :: equivalence
      !> Whether the default threshold was tried already
      logical :: default_tried

      if (.not. present(tolerance)) then
         call detect_point_group(xyz, species, equivalence_default, angle_default, group, error)
         if (allocated(error)) then
            call detect_point_group(xyz, species, equivalence_default, angle_tight, loose, loose_error)
            if (.not. allocated(loose_error)) then
               deallocate (error)
               group = loose
            end if
         end if
         return
      end if

      equivalence = min(equivalence_default, max(equivalence_floor, equivalence_per_tolerance*tolerance))
      default_tried = equivalence >= equivalence_default
      call detect_ladder(angle_default, group, error)

      if (size(xyz, 2) <= hint_max_atoms .and. .not. default_tried) then
         call detect_point_group(xyz, species, equivalence_default, angle_default, loose, loose_error)
         if (.not. allocated(loose_error)) then
            if (allocated(error)) then
               deallocate (error)
               group = loose
            else if (group_size(loose) > group_size(group)) then
               group = loose
            end if
         end if
      end if
      if (.not. allocated(error)) return

      ! Every threshold failed; that error stays when this fails as well
      call detect_ladder(angle_tight, loose, loose_error)
      if (allocated(loose_error)) return
      deallocate (error)
      group = loose

   contains

      !> Detection at the caller's equivalence threshold, then at the tighter ones
      !>
      !> @param[in]  angle libmsym angle threshold; not positive keeps libmsym's
      !> @param[out] found Detected group
      !> @param[out] err   Set when every threshold fails
      subroutine detect_ladder(angle, found, err)
         !> Angle threshold
         real(wp), intent(in) :: angle
         !> Detected group
         type(msym_group_type), intent(out) :: found
         !> Error handling
         type(error_type), allocatable, intent(out) :: err

         !> Equivalence threshold of the current detection
         real(wp) :: eq
         integer :: k

         eq = equivalence
         call detect_point_group(xyz, species, eq, angle, found, err)
         do k = 1, size(equivalence_ladder)
            if (.not. allocated(err)) exit
            if (equivalence_ladder(k) >= eq) cycle
            eq = equivalence_ladder(k)
            deallocate (err)
            call detect_point_group(xyz, species, eq, angle, found, err)
         end do
      end subroutine detect_ladder

   end subroutine find_point_group

   !> Size of a group for comparing two detections
   !>
   !> - a finite group counts its operations
   !> - a linear group is larger than any finite one; libmsym lists only a
   !>   few of its operations (2 for C_inf_v, 4 for D_inf_h), which still
   !>   rank D_inf_h above C_inf_v
   !>
   !> @param[in] group Detected group
   pure function group_size(group) result(rank)
      !> Detected group
      type(msym_group_type), intent(in) :: group
      !> Rank, larger for the larger group
      integer :: rank

      rank = size(group%ops, 3)
      if (group%linear) rank = rank + linear_rank
   end function group_size

   !> Detect the point group at one equivalence threshold
   !>
   !> @param[in]  xyz         Coordinates (3, nat)
   !> @param[in]  species     Species id per atom (nat), positive
   !> @param[in]  equivalence libmsym equivalence threshold, relative
   !> @param[in]  angle       libmsym angle threshold; not positive keeps libmsym's
   !> @param[out] group       Detected group, operations and subgroups
   !> @param[out] error       Set when libmsym fails
   subroutine detect_point_group(xyz, species, equivalence, angle, group, error)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Species ids
      integer, intent(in) :: species(:)
      !> Equivalence threshold
      real(wp), intent(in) :: equivalence
      !> Angle threshold
      real(wp), intent(in) :: angle
      !> Detected group
      type(msym_group_type), intent(out) :: group
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> libmsym context
      type(c_ptr) :: handle
      !> Shim status
      integer(c_int) :: status
      !> Group name, NUL padded
      character(kind=c_char) :: cname(8)
      !> Group type, principal order and order
      integer(c_int) :: kind, n, order
      !> Number of subgroups and order of one
      integer(c_int) :: nsub, sub_order
      !> Operation indices of a subgroup, 0-based
      integer(c_int), allocatable :: sub_ops(:)
      !> Flat operation matrices and codes
      real(c_double), allocatable :: mats(:)
      integer(c_int), allocatable :: info(:)
      !> Centroid
      real(c_double) :: center(3)
      !> Axis of one operation
      real(c_double) :: axis(3)
      integer :: nat, iop, isub

      nat = size(xyz, 2)
      if (size(xyz, 1) /= 3 .or. size(species) /= nat) then
         call fatal_error(error, "find_point_group: xyz must be (3, nat) and species (nat)")
         return
      end if
      if (nat < 2) then
         call fatal_error(error, "find_point_group: at least two atoms are needed")
         return
      end if
      if (any(species < 1)) then
         call fatal_error(error, "find_point_group: species ids must be positive")
         return
      end if

      ! libmsym keeps its error details in a global buffer
      !$omp critical (moist_msym)
      call query()
      !$omp end critical (moist_msym)

   contains

      !> All libmsym work: create, detect, read back, release
      subroutine query()
         status = moist_msym_create(handle)
         if (status /= 0 .or. .not. c_associated(handle)) then
            call fail(status, "creating the libmsym context")
            return
         end if

         status = moist_msym_find(handle, int(nat, c_int), real(reshape(xyz, [3*nat]), c_double), &
            & int(species, c_int), real(equivalence, c_double), real(angle, c_double))
         if (status /= 0) then
            call fail(status, "detecting the point group")
            call moist_msym_release(handle)
            return
         end if

         status = moist_msym_group(handle, cname, kind, n, order)
         if (status /= 0) then
            call fail(status, "reading the point group")
            call moist_msym_release(handle)
            return
         end if
         call from_c(cname, group%name)
         group%kind = kind
         group%n = n
         group%linear = n == 0 .and. (kind == msym_type_cnv .or. kind == msym_type_dnh)

         allocate (mats(9*order), info(4*order))
         status = moist_msym_operations(handle, order, mats, info, center)
         if (status /= 0) then
            call fail(status, "reading the symmetry operations")
            call moist_msym_release(handle)
            return
         end if
         group%ops = reshape(real(mats, wp), [3, 3, int(order)])
         group%info = reshape(int(info), [4, int(order)])
         group%center = real(center, wp)

         allocate (group%axis(3, order), source=0.0_wp)
         do iop = 1, order
            if (group%info(1, iop) == msym_op_identity .or. group%info(1, iop) == msym_op_inversion) cycle
            status = moist_msym_axis(handle, int(iop - 1, c_int), axis)
            if (status /= 0) then
               call fail(status, "reading an operation axis")
               call moist_msym_release(handle)
               return
            end if
            group%axis(:, iop) = real(axis, wp)
            if (norm2(group%axis(:, iop)) > 0.0_wp) then
               group%axis(:, iop) = group%axis(:, iop)/norm2(group%axis(:, iop))
            end if
         end do

         if (group%linear) then
            allocate (group%subgroups(0))
            call moist_msym_release(handle)
            return
         end if

         status = moist_msym_subgroup_count(handle, nsub)
         if (status /= 0) then
            ! C1 and other groups without proper subgroups report none
            nsub = 0
         end if
         allocate (group%subgroups(nsub), sub_ops(order))
         do isub = 1, nsub
            status = moist_msym_subgroup(handle, int(isub - 1, c_int), cname, sub_order, order, sub_ops)
            if (status /= 0) then
               call fail(status, "reading a subgroup")
               call moist_msym_release(handle)
               return
            end if
            call from_c(cname, group%subgroups(isub)%name)
            group%subgroups(isub)%ops = int(sub_ops(1:sub_order)) + 1
         end do

         call moist_msym_release(handle)
      end subroutine query

      !> Turn a shim status into an error
      !>
      !> @param[in] code Shim status
      !> @param[in] what Step that failed
      subroutine fail(code, what)
         !> Shim status
         integer(c_int), intent(in) :: code
         !> Step that failed
         character(len=*), intent(in) :: what

         !> Rendered message
         character(kind=c_char) :: buf(msg_len)
         !> Fortran copy of the message
         character(len=:), allocatable :: msg

         call moist_msym_error_message(code, buf, int(msg_len, c_int))
         call from_c(buf, msg)
         call fatal_error(error, "libmsym failed "//what//": "//msg)
      end subroutine fail

   end subroutine detect_point_group

   !> Fortran string from a NUL-terminated C character array
   !>
   !> @param[in]  buf C characters
   !> @param[out] str Trimmed string up to the first NUL
   pure subroutine from_c(buf, str)
      !> C characters
      character(kind=c_char), intent(in) :: buf(:)
      !> Trimmed string up to the first NUL
      character(len=:), allocatable, intent(out) :: str

      integer :: i, n

      n = size(buf)
      do i = 1, size(buf)
         if (buf(i) == c_null_char) then
            n = i - 1
            exit
         end if
      end do
      allocate (character(len=n) :: str)
      do i = 1, n
         str(i:i) = buf(i)
      end do
      str = trim(str)
   end subroutine from_c

end module moist_math_symmetry_msym
