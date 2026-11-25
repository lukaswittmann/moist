!> Lattice helpers of the surface mesh seeder ([[projection_seeder_type]])
!>
!> - [[key_table_type]], [[lattice_key]]: open-addressing table from a 64-bit
!>   cell key to an index, and the key of three lattice coordinates
!> - [[build_probe_type]], [[evaluate_slices]]: level set jets of a build, in
!>   fixed slices over per-thread clones
!> - [[segment_crossing]]: crossing of a sign-change segment
module moist_cavity_drop_projection_seed_lattice
   !$ use omp_lib, only: omp_get_thread_num
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   implicit none(type, external)
   private

   public :: key_bits, key_table_type, table_init, table_find, table_insert, lattice_key
   public :: build_probe_type, evaluate_slices, segment_crossing

   !> Bits of each of the three fields of a cell key
   integer, parameter :: key_bits = 21
   !> Spatial hash primes of the key tables
   !>
   !> - mix all three 21-bit coordinates into the low slot bits without int64 overflow
   integer(int64), parameter :: hash_mult(3) = [73856093_int64, 19349663_int64, 83492791_int64]

   !> Failure of one slice of a build's level set call
   type :: slice_error_type
      !> Failure of the slice, unallocated on success
      type(error_type), allocatable :: error
   end type slice_error_type

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

   !> Level set clones of a build, one per thread
   type :: build_probe_type
      !> Clones, allocated by the thread that uses them first
      type(lsf_thread_slot), allocatable :: slot(:)
      !> Threads of the probe batches
      integer :: nthreads = 1
   end type build_probe_type

contains

   !* ================================================================================= *!
   !*                                    Key tables                                     *!
   !* ================================================================================= *!

   !> Empty a table for about `n` entries; unallocated when `n` is not positive
   !>
   !> @param[inout] table Table
   !> @param[in]    n     Expected entries
   pure subroutine table_init(table, n)
      type(key_table_type), intent(inout) :: table
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
      type(key_table_type), intent(in) :: table
      integer(int64), intent(in) :: key
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
      type(key_table_type), intent(in) :: table
      integer(int64), intent(in) :: key
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

   !> Key of three lattice coordinates, each field `c(i) + bias` in `0 .. 2^key_bits - 1`
   !>
   !> @param[in] c    Coordinates, the first in the lowest field
   !> @param[in] bias Offset making every field non-negative
   pure function lattice_key(c, bias) result(key)
      integer, intent(in) :: c(3)
      integer, intent(in) :: bias
      integer(int64) :: key

      key = ior(ior(int(c(1) + bias, int64), shiftl(int(c(2) + bias, int64), key_bits)), &
                shiftl(int(c(3) + bias, int64), 2*key_bits))
   end function lattice_key

   !> Store an absent key, growing the table at half load
   !>
   !> @param[inout] table Table
   !> @param[in]    key   Key
   !> @param[in]    idx   Index, positive
   pure recursive subroutine table_insert(table, key, idx)
      type(key_table_type), intent(inout) :: table
      integer(int64), intent(in) :: key
      integer, intent(in) :: idx

      integer(int64), allocatable :: old_key(:)
      integer, allocatable :: old_idx(:)
      integer :: slot, k

      if (.not. allocated(table%idx)) call table_init(table, 16)
      if (2*(table%count + 1) > size(table%idx)) then
         call move_alloc(table%key, old_key)
         call move_alloc(table%idx, old_idx)
         allocate (table%key(2*size(old_idx)), table%idx(2*size(old_idx)))
         table%idx = 0
         table%count = 0
         do k = 1, size(old_idx)
            if (old_idx(k) /= 0) call table_insert(table, old_key(k), old_idx(k))
         end do
      end if
      slot = table_slot(table, key)
      do while (table%idx(slot) /= 0)
         slot = slot + 1
         if (slot > size(table%idx)) slot = 1
      end do
      table%key(slot) = key
      table%idx(slot) = idx
      table%count = table%count + 1
   end subroutine table_insert

   !* ================================================================================= *!
   !*                                      Probes                                       *!
   !* ================================================================================= *!

   !> Level set jets at a batch of points, in fixed slices over the thread clones
   !>
   !> - slice `i` holds points `(i - 1) chunk + 1 .. i chunk`, whatever the thread
   !>   count; the threads take the slices in any order, each writes only its own
   !> - every slice starts on a clone without evaluation history
   !>   (`reset_history`), so no value depends on the slices a thread did before
   !> - `reached` false where the level set does not reach the point
   !> - on failure the error of the first failing slice
   !>
   !> @param[inout] probe   Build probe (clones, allocated by the thread that uses them first)
   !> @param[in]    lsf     Level set, the source of the clones
   !> @param[in]    x       Points (3, n)
   !> @param[in]    order   Derivative order, 0..2
   !> @param[in]    chunk   Points per slice
   !> @param[out]   value   Values (n)
   !> @param[out]   reached Whether the level set reaches the point (n)
   !> @param[out]   error   LSF failure
   !> @param[out]   grad    Gradients (3, n), order >= 1
   !> @param[out]   hess    Packed Hessians (6, n), order 2
   subroutine evaluate_slices(probe, lsf, x, order, chunk, value, reached, error, grad, hess)
      type(build_probe_type), intent(inout) :: probe
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      real(wp), intent(in), contiguous :: x(:, :)
      integer, intent(in) :: order
      integer, intent(in) :: chunk
      real(wp), intent(out), contiguous :: value(:)
      logical, intent(out), contiguous :: reached(:)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(out), contiguous, optional :: grad(:, :)
      real(wp), intent(out), contiguous, optional :: hess(:, :)

      type(slice_error_type), allocatable :: failed(:)
      integer :: n, nslice, islice, i0, i1, k, tid

      n = size(x, 2)
      if (n == 0) return
      nslice = (n + chunk - 1)/chunk
      allocate (failed(nslice))
      !$omp parallel do num_threads(probe%nthreads) schedule(dynamic) default(shared) &
      !$omp& private(islice, i0, i1, k, tid)
      do islice = 1, nslice
         tid = 1
         !$ tid = omp_get_thread_num() + 1
         if (.not. allocated(probe%slot(tid)%lsf)) allocate (probe%slot(tid)%lsf, source=lsf)
         call probe%slot(tid)%lsf%reset_history()
         i0 = (islice - 1)*chunk + 1
         i1 = min(n, islice*chunk)
         select case (order)
         case (0)
            call probe%slot(tid)%lsf%evaluate_batch(x(:, i0:i1), 0, value(i0:i1), error=failed(islice)%error)
         case (1)
            call probe%slot(tid)%lsf%evaluate_batch(x(:, i0:i1), 1, value(i0:i1), lsf1_r=grad(:, i0:i1), &
                                                    error=failed(islice)%error)
         case default
            call probe%slot(tid)%lsf%evaluate_batch(x(:, i0:i1), 2, value(i0:i1), lsf1_r=grad(:, i0:i1), &
                                                    lsf2_rr_packed=hess(:, i0:i1), error=failed(islice)%error)
         end select
         if (.not. allocated(failed(islice)%error)) then
            do k = i0, i1
               reached(k) = .not. probe%slot(tid)%lsf%unreached(k - i0 + 1)
            end do
         end if
      end do
      !$omp end parallel do
      do islice = 1, nslice
         if (allocated(failed(islice)%error)) then
            call move_alloc(failed(islice)%error, error)
            return
         end if
      end do
   end subroutine evaluate_slices

   !> Crossing of a sign-change segment by linear interpolation of its end values
   !>
   !> @param[in] a  First end point
   !> @param[in] b  Second end point
   !> @param[in] va Value at `a`
   !> @param[in] vb Value at `b`, on the other side of the surface
   pure function segment_crossing(a, b, va, vb) result(x)
      real(wp), intent(in) :: a(3)
      real(wp), intent(in) :: b(3)
      real(wp), intent(in) :: va
      real(wp), intent(in) :: vb
      real(wp) :: x(3)

      x = a + (va/(va - vb))*(b - a)
   end function segment_crossing

end module moist_cavity_drop_projection_seed_lattice
