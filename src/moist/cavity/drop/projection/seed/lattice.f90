!> Lattice helpers of the lattice seeder ([[projection_seeder_type]])
!>
!> - [[key_table_type]], [[lattice_key]]: open-addressing table from a 64-bit
!>   cell key to an index, and the key of three lattice coordinates
!> - [[evaluate_slices]]: level set jets of a build, in fixed slices over
!>   per-thread clones ([[drop_worker_slots_type]])
!> - [[segment_crossing]]: crossing of a sign-change segment
!> - [[lattice_frame_type]]: origin and axes of the lattice, maps between the
!>   lab and the lattice frame
module moist_cavity_drop_projection_seed_lattice
   !$ use omp_lib, only: omp_get_thread_num
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_math_linalg_geometry, only: identity3
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: key_bits, key_bias, key_table_type, table_init, table_find, table_insert, lattice_key
   public :: evaluate_slices, segment_crossing, lattice_frame_type

   !> Bits of each of the three fields of a cell key
   integer, parameter :: key_bits = 21
   !> Offset making a lattice coordinate non-negative
   integer, parameter :: key_bias = 2**(key_bits - 1)
   !> Spatial hash primes of the key tables
   !>
   !> - mix all three 21-bit coordinates into the low slot bits without int64 overflow
   integer(int64), parameter :: hash_mult(3) = [73856093_int64, 19349663_int64, 83492791_int64]

   !> Open-addressing table from a cell key to an index
   !>
   !> - a key holds three non-negative fields of `key_bits` bits
   type :: key_table_type
      !> Stored keys
      integer(int64), allocatable :: key(:)
      !> Index of each slot, 0 when empty
      integer, allocatable :: idx(:)
      !> Entries
      integer :: count = 0
   end type key_table_type

   !> Lattice frame: lattice-frame positions `x_l = o + A^T (x - o)`
   !>
   !> - lab axes (`A = I`) skip the maps, positions are returned bitwise
   type :: lattice_frame_type
      !> Lattice origin `o` (Bohr)
      real(wp) :: origin(3) = 0.0_wp
      !> Lattice axes `A` as columns
      real(wp) :: axes(3, 3) = identity3
      !> Whether the axes differ from the lab axes
      logical :: rotated = .false.
   contains
      !> Set origin and axes
      procedure :: init => frame_init
      !> Lab position of a lattice-frame position
      procedure :: to_lab => frame_to_lab
      !> Lattice-frame position of a lab position
      procedure :: to_lattice => frame_to_lattice
   end type lattice_frame_type

contains

   !* ================================================================================= *!
   !*                                    Key tables                                     *!
   !* ================================================================================= *!

   !> Empty a table for about `n` entries; unallocated when `n` is not positive
   !>
   !> @param[inout] table Table
   !> @param[in]    n     Expected entries
   pure subroutine table_init(table, n)
      !> Table
      type(key_table_type), intent(inout) :: table
      !> Expected entries
      integer, intent(in) :: n

      integer :: cap

      cap = 16
      do while (cap < 2*n)
         cap = 2*cap
      end do
      if (allocated(table%key)) deallocate (table%key, table%idx)
      table%count = 0
      if (n <= 0) return
      allocate (table%key(cap), table%idx(cap))
      table%idx = 0
   end subroutine table_init

   !> First slot of a key
   !>
   !> @param[in] table Table
   !> @param[in] key   Key
   pure function table_slot(table, key) result(slot)
      !> Table
      type(key_table_type), intent(in) :: table
      !> Key
      integer(int64), intent(in) :: key
      !> Slot, 1..size(table%idx)
      integer :: slot

      integer(int64), parameter :: mask = 2_int64**key_bits - 1_int64
      integer(int64) :: h

      h = ieor(ieor(iand(key, mask)*hash_mult(1), iand(shiftr(key, key_bits), mask)*hash_mult(2)), &
               iand(shiftr(key, 2*key_bits), mask)*hash_mult(3))
      slot = int(iand(h, int(size(table%idx) - 1, int64))) + 1
   end function table_slot

   !> Index stored under a key, 0 when absent
   !>
   !> @param[in] table Table
   !> @param[in] key   Key
   pure function table_find(table, key) result(idx)
      !> Table
      type(key_table_type), intent(in) :: table
      !> Key
      integer(int64), intent(in) :: key
      !> Stored index, 0 when absent
      integer :: idx

      integer :: slot

      idx = 0
      if (.not. allocated(table%idx)) return
      slot = table_slot(table, key)
      do while (table%idx(slot) /= 0)
         if (table%key(slot) == key) then
            idx = table%idx(slot)
            return
         end if
         slot = slot + 1
         if (slot > size(table%idx)) slot = 1
      end do
   end function table_find

   !> Key of three lattice coordinates, each field `c(i) + key_bias` in `0 .. 2^key_bits - 1`
   !>
   !> @param[in] c Coordinates, the first in the lowest field
   pure function lattice_key(c) result(key)
      !> Coordinates, the first in the lowest field
      integer, intent(in) :: c(3)
      !> Key
      integer(int64) :: key

      key = ior(ior(int(c(1) + key_bias, int64), shiftl(int(c(2) + key_bias, int64), key_bits)), &
                shiftl(int(c(3) + key_bias, int64), 2*key_bits))
   end function lattice_key

   !> Store an absent key, growing the table at half load
   !>
   !> - growth rehashes the stored keys in slot order into the doubled table
   !>
   !> @param[inout] table Table
   !> @param[in]    key   Key
   !> @param[in]    idx   Index, positive
   pure subroutine table_insert(table, key, idx)
      !> Table
      type(key_table_type), intent(inout) :: table
      !> Key
      integer(int64), intent(in) :: key
      !> Index, positive
      integer, intent(in) :: idx

      integer(int64), allocatable :: old_key(:)
      integer, allocatable :: old_idx(:)
      integer :: k

      if (.not. allocated(table%idx)) call table_init(table, 16)
      if (2*(table%count + 1) > size(table%idx)) then
         call move_alloc(table%key, old_key)
         call move_alloc(table%idx, old_idx)
         allocate (table%key(2*size(old_idx)), table%idx(2*size(old_idx)))
         table%idx = 0
         table%count = 0
         do k = 1, size(old_idx)
            if (old_idx(k) /= 0) call table_place(table, old_key(k), old_idx(k))
         end do
      end if
      call table_place(table, key, idx)
   end subroutine table_insert

   !> Store an absent key in the first free slot from its own, no growth
   !>
   !> @param[inout] table Table with a free slot
   !> @param[in]    key   Key
   !> @param[in]    idx   Index, positive
   pure subroutine table_place(table, key, idx)
      !> Table with a free slot
      type(key_table_type), intent(inout) :: table
      !> Key
      integer(int64), intent(in) :: key
      !> Index, positive
      integer, intent(in) :: idx

      integer :: slot

      slot = table_slot(table, key)
      do while (table%idx(slot) /= 0)
         slot = slot + 1
         if (slot > size(table%idx)) slot = 1
      end do
      table%key(slot) = key
      table%idx(slot) = idx
      table%count = table%count + 1
   end subroutine table_place

   !* ================================================================================= *!
   !*                                      Probes                                       *!
   !* ================================================================================= *!

   !> Level set jets at a batch of points, in fixed slices over the thread clones
   !>
   !> - slice `i` holds points `(i - 1) chunk + 1 .. i chunk`, whatever the thread
   !>   count; the threads take the slices in any order, each writes only its own
   !> - a thread clones `lsf` into its slot on first use ([[worker_slots_fill]]);
   !>   the slot keeps the clone for later calls until the caller releases it
   !> - every slice starts on a clone without evaluation history
   !>   (`reset_history`), so no value depends on the slices a thread did before
   !> - `reached` false where the level set does not reach the point
   !> - on failure no slice behind the failing one evaluates (a host callback is
   !>   not called again after it failed); the error of the first failing slice
   !>
   !> @param[inout] slots   Clones of the build, sized; filled on first use
   !> @param[in]    lsf     Level set, the source of the clones; may be borrowed
   !>                       from until the slots are released
   !> @param[in]    x       Points (3, n)
   !> @param[in]    order   Derivative order, 0..2
   !> @param[in]    chunk   Points per slice
   !> @param[out]   value   Values (n)
   !> @param[out]   reached Whether the level set reaches the point (n)
   !> @param[out]   error   LSF failure, order out of range or its output missing
   !> @param[out]   grad    Gradients (3, n), required for order >= 1
   !> @param[out]   hess    Packed Hessians (6, n), required for order 2
   subroutine evaluate_slices(slots, lsf, x, order, chunk, value, reached, error, grad, hess)
      !> Clones of the build, sized; filled on first use
      type(drop_worker_slots_type), intent(inout) :: slots
      !> Level set, the source of the clones; may be borrowed
      !> from until the slots are released
      class(moist_cavity_drop_lsf_type), intent(in), target :: lsf
      !> Points (3, n)
      real(wp), intent(in), contiguous :: x(:, :)
      !> Derivative order, 0..2
      integer, intent(in) :: order
      !> Points per slice
      integer, intent(in) :: chunk
      !> Values (n)
      real(wp), intent(out), contiguous :: value(:)
      !> Whether the level set reaches the point (n)
      logical, intent(out), contiguous :: reached(:)
      !> LSF failure, order out of range or its output missing
      type(error_type), allocatable, intent(out) :: error
      !> Gradients (3, n), required for order >= 1
      real(wp), intent(out), contiguous, optional :: grad(:, :)
      !> Packed Hessians (6, n), required for order 2
      real(wp), intent(out), contiguous, optional :: hess(:, :)

      type(drop_abort_latch_type) :: latch
      type(error_type), allocatable :: failed
      integer :: n, nslice, islice, i0, i1, k, tid

      if (order < 0 .or. order > 2) then
         call fatal_error(error, "Lattice seeder: level set order "//to_string(order)//" outside 0..2")
         return
      end if
      if (order >= 1 .and. .not. present(grad)) then
         call fatal_error(error, "Lattice seeder: level set order "//to_string(order)//" without gradients")
         return
      end if
      if (order == 2 .and. .not. present(hess)) then
         call fatal_error(error, "Lattice seeder: level set order 2 without Hessians")
         return
      end if
      n = size(x, 2)
      if (n == 0) return
      nslice = (n + chunk - 1)/chunk
      call latch%reset()
      !$omp parallel do num_threads(slots%nthreads) schedule(dynamic) default(shared) &
      !$omp& private(islice, i0, i1, k, tid, failed)
      do islice = 1, nslice
         ! Lowest failing slice wins, none behind it calls the level set
         if (latch%skips(islice)) cycle
         tid = 1
         !$ tid = omp_get_thread_num() + 1
         if (.not. allocated(slots%lsf(tid)%lsf)) call slots%fill(tid, lsf, lsf%max_deriv)
         call slots%lsf(tid)%lsf%reset_history()
         i0 = (islice - 1)*chunk + 1
         i1 = min(n, islice*chunk)
         select case (order)
         case (0)
            call slots%lsf(tid)%lsf%evaluate_batch(x(:, i0:i1), 0, value(i0:i1), error=failed)
         case (1)
            call slots%lsf(tid)%lsf%evaluate_batch(x(:, i0:i1), 1, value(i0:i1), lsf1_r=grad(:, i0:i1), &
                                                   error=failed)
         case default
            call slots%lsf(tid)%lsf%evaluate_batch(x(:, i0:i1), 2, value(i0:i1), lsf1_r=grad(:, i0:i1), &
                                                   lsf2_rr_packed=hess(:, i0:i1), error=failed)
         end select
         if (allocated(failed)) then
            call latch%latch_error(failed, islice)
            cycle
         end if
         do k = i0, i1
            reached(k) = .not. slots%lsf(tid)%lsf%unreached(k - i0 + 1)
         end do
      end do
      !$omp end parallel do
      if (latch%is_set()) call move_alloc(latch%error, error)
   end subroutine evaluate_slices

   !> Crossing of a sign-change segment by linear interpolation of its end values
   !>
   !> @param[in] a  First end point
   !> @param[in] b  Second end point
   !> @param[in] va Value at `a`
   !> @param[in] vb Value at `b`, on the other side of the surface
   pure function segment_crossing(a, b, va, vb) result(x)
      !> First end point
      real(wp), intent(in) :: a(3)
      !> Second end point
      real(wp), intent(in) :: b(3)
      !> Value at `a`
      real(wp), intent(in) :: va
      !> Value at `b`, on the other side of the surface
      real(wp), intent(in) :: vb
      !> Crossing
      real(wp) :: x(3)

      x = a + (va/(va - vb))*(b - a)
   end function segment_crossing

   !* ================================================================================= *!
   !*                                   Lattice frame                                   *!
   !* ================================================================================= *!

   !> Set origin and axes of the frame
   !>
   !> @param[inout] self   Frame
   !> @param[in]    origin Lattice origin `o`
   !> @param[in]    axes   Lattice axes `A` as columns, a proper rotation
   pure subroutine frame_init(self, origin, axes)
      !> Frame
      class(lattice_frame_type), intent(inout) :: self
      !> Lattice origin `o`
      real(wp), intent(in) :: origin(3)
      !> Lattice axes `A` as columns, a proper rotation
      real(wp), intent(in) :: axes(3, 3)

      self%origin = origin
      self%axes = axes
      self%rotated = any(axes /= identity3)
   end subroutine frame_init

   !> Lab position of a lattice-frame position: `o + A (x_l - o)`
   !>
   !> @param[in] self Frame
   !> @param[in] xl   Lattice-frame position
   pure function frame_to_lab(self, xl) result(x)
      !> Frame
      class(lattice_frame_type), intent(in) :: self
      !> Lattice-frame position
      real(wp), intent(in) :: xl(3)
      !> Lab position
      real(wp) :: x(3)

      x = xl
      if (self%rotated) x = self%origin + matmul(self%axes, xl - self%origin)
   end function frame_to_lab

   !> Lattice-frame position of a lab position: `o + A^T (x - o)`
   !>
   !> @param[in] self Frame
   !> @param[in] x    Lab position
   pure function frame_to_lattice(self, x) result(xl)
      !> Frame
      class(lattice_frame_type), intent(in) :: self
      !> Lab position
      real(wp), intent(in) :: x(3)
      !> Lattice-frame position
      real(wp) :: xl(3)

      xl = x
      if (self%rotated) xl = self%origin + matmul(x - self%origin, self%axes)
   end function frame_to_lattice

end module moist_cavity_drop_projection_seed_lattice
