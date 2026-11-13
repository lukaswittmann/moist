!> Tiled level-set evaluation at projected surface points
!>
!> Tile: a batch of surface_tile_width surface points
!>
!> Shared by the surface properties, the reverse nuclear gradient and the
!> response contraction: one shared candidate list per tile, the order-1 to
!> order-3 jet through `evaluate_batch`, the seed-contracted nuclear partials
!> through `vjp_f1_rA_batch`
!>
module moist_cavity_drop_tile
   use mctc_env, only: error_type
   use mctc_env_accuracy, only: wp
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_threads, only: drop_abort_latch_type
   implicit none(type, external)
   private

   public :: surface_tile_width, surface_tile_candidates, surface_jet_tile, surface_vjp_tile

   !> Points per surface tile of the derivative loops, scheduled `dynamic, 1`
   integer, parameter :: surface_tile_width = 64

contains

   !> Shared candidate list of a tile: every atom whose screening reach meets
   !> the tile's bounding box
   !>
   !> Superset of the per-point gate for every point of the tile, ascending in
   !> the LSF's `candidate_space`
   !>
   !> @param[in]  lsf    LSF clone, after `update`
   !> @param[in]  points Tile points (3, np), np >= 1
   !> @param[out] cand   Candidate ids in `cand(1:nc)`, size >= ncenters
   !> @param[out] nc     Number of candidates
   subroutine surface_tile_candidates(lsf, points, cand, nc)
      !> LSF clone
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Tile points
      real(wp), intent(in) :: points(:, :)
      !> Candidate ids
      integer, intent(out) :: cand(:)
      !> Number of candidates
      integer, intent(out) :: nc

      !> Bounding box of the tile
      real(wp) :: lo(3), hi(3)
      integer :: p

      lo = points(:, 1)
      hi = points(:, 1)
      do p = 2, size(points, 2)
         lo = min(lo, points(:, p))
         hi = max(hi, points(:, p))
      end do
      call lsf%candidates_in_box(lo, hi, cand, nc)
   end subroutine surface_tile_candidates

   !> Level-set jet of one surface tile, order 1 to 3
   !>
   !> - packed columns as `evaluate_batch`: Hessian (6, np), third (10, np)
   !> - `lsf2_rr` required and written for order >= 2, `lsf3_rrr` for order 3
   !> - an `error`, or an unreached column (`lsf%unreached(column)`), is latched
   !>   with its grid point and `ok` is false; the tile must then be dropped whole
   !>
   !> @param[inout] lsf      LSF clone of this thread
   !> @param[in]    points   Tile points (3, np)
   !> @param[in]    cand     Shared candidate list of the tile, ascending
   !> @param[in]    igrid0   Grid index of the first tile column
   !> @param[in]    who      Calling routine, prefixes the diagnostic
   !> @param[in]    order    Highest spatial derivative, 1 to 3
   !> @param[out]   lsf0     Values (>= np)
   !> @param[out]   lsf1_r   Gradients (3, >= np)
   !> @param[out]   lsf2_rr  Packed Hessians (6, >= np), order >= 2
   !> @param[out]   lsf3_rrr Packed third derivatives (10, >= np), order 3
   !> @param[inout] abort    Shared abort latch
   !> @returns      ok       Every column evaluated
   function surface_jet_tile(lsf, points, cand, igrid0, who, order, lsf0, lsf1_r, lsf2_rr, &
                             lsf3_rrr, abort) result(ok)
      !> LSF clone of this thread
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Tile points
      real(wp), intent(in), contiguous :: points(:, :)
      !> Shared candidate list
      integer, intent(in) :: cand(:)
      !> Grid index of the first tile column
      integer, intent(in) :: igrid0
      !> Calling routine
      character(len=*), intent(in) :: who
      !> Highest spatial derivative
      integer, intent(in) :: order
      !> Values
      real(wp), intent(out), contiguous :: lsf0(:)
      !> Gradients
      real(wp), intent(out), contiguous :: lsf1_r(:, :)
      !> Packed Hessians
      real(wp), intent(out), contiguous, optional :: lsf2_rr(:, :)
      !> Packed third derivatives
      real(wp), intent(out), contiguous, optional :: lsf3_rrr(:, :)
      !> Shared abort latch
      type(drop_abort_latch_type), intent(inout) :: abort
      !> Every column evaluated
      logical :: ok

      type(error_type), allocatable :: error

      select case (order)
      case (1)
         call lsf%evaluate_batch(points, 1, lsf0, lsf1_r=lsf1_r, candidate_indices=cand, &
                                 error=error)
      case (2)
         call lsf%evaluate_batch(points, 2, lsf0, lsf1_r=lsf1_r, lsf2_rr_packed=lsf2_rr, &
                                 candidate_indices=cand, error=error)
      case (3)
         call lsf%evaluate_batch(points, 3, lsf0, lsf1_r=lsf1_r, lsf2_rr_packed=lsf2_rr, &
                                 lsf3_rrr=lsf3_rrr, candidate_indices=cand, error=error)
      case default
         error stop who//": surface_jet_tile order must be 1 to 3"
      end select
      ok = tile_ok(lsf, error, size(points, 2), igrid0, who, "level-set evaluation", abort)
   end function surface_jet_tile

   !> Seed-contracted nuclear partials of the level-set jet of one surface tile
   !>
   !> - `w2` the packed symmetric part of the Hessian seed, not doubled
   !> - `res(3, nc, np)`: column k is candidate `cand(k)`, user atom
   !>   `cand_to_user(cand(k))`
   !> - failures latched as [[surface_jet_tile]], including a non-finite
   !>   contraction of a finite jet
   !>
   !> @param[inout] lsf    LSF clone of this thread
   !> @param[in]    points Tile points (3, np)
   !> @param[in]    cand   Shared candidate list of the tile, ascending
   !> @param[in]    w0     Value seeds (>= np)
   !> @param[in]    w1     Gradient seeds (3, >= np)
   !> @param[in]    w2     Packed symmetric Hessian seeds (6, >= np)
   !> @param[in]    igrid0 Grid index of the first tile column
   !> @param[in]    who    Calling routine, prefixes the diagnostic
   !> @param[out]   res    Contracted nuclear gradient (3, nc, np)
   !> @param[inout] abort  Shared abort latch
   !> @returns      ok     Every column contracted
   function surface_vjp_tile(lsf, points, cand, w0, w1, w2, igrid0, who, res, abort) &
      result(ok)
      !> LSF clone of this thread
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Tile points
      real(wp), intent(in), contiguous :: points(:, :)
      !> Shared candidate list
      integer, intent(in) :: cand(:)
      !> Value seeds
      real(wp), intent(in), contiguous :: w0(:)
      !> Gradient seeds
      real(wp), intent(in), contiguous :: w1(:, :)
      !> Packed symmetric Hessian seeds
      real(wp), intent(in), contiguous :: w2(:, :)
      !> Grid index of the first tile column
      integer, intent(in) :: igrid0
      !> Calling routine
      character(len=*), intent(in) :: who
      !> Contracted nuclear gradient
      real(wp), intent(out), contiguous :: res(:, :, :)
      !> Shared abort latch
      type(drop_abort_latch_type), intent(inout) :: abort
      !> Every column contracted
      logical :: ok

      type(error_type), allocatable :: error

      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=error, candidate_indices=cand)
      ok = tile_ok(lsf, error, size(points, 2), igrid0, who, "nuclear level-set contraction", &
                   abort)
   end function surface_vjp_tile

   !> Latch the first failure of a batched tile call
   !>
   !> - an `error` first: its message names the tile column
   !> - else the first unreached column (`lsf%unreached(column)`), by grid point
   !>
   !> @param[in]    lsf    LSF clone right after the batch call
   !> @param[inout] error  Batch error, consumed
   !> @param[in]    np     Tile size
   !> @param[in]    igrid0 Grid index of the first tile column
   !> @param[in]    who    Calling routine
   !> @param[in]    what   Failed stage
   !> @param[inout] abort  Shared abort latch
   !> @returns      ok     No failure
   function tile_ok(lsf, error, np, igrid0, who, what, abort) result(ok)
      !> LSF clone right after the batch call
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Batch error
      type(error_type), allocatable, intent(inout) :: error
      !> Tile size
      integer, intent(in) :: np
      !> Grid index of the first tile column
      integer, intent(in) :: igrid0
      !> Calling routine
      character(len=*), intent(in) :: who
      !> Failed stage
      character(len=*), intent(in) :: what
      !> Shared abort latch
      type(drop_abort_latch_type), intent(inout) :: abort
      !> No failure
      logical :: ok

      character(len=16) :: idx
      integer :: p

      ok = .false.
      if (allocated(error)) then
         write (idx, "(i0)") igrid0
         error%message = who//": "//what//" failed on the surface tile from grid point "// &
                         trim(idx)//": "//error%message
         call abort%latch_error(error, igrid0)
         return
      end if
      do p = 1, np
         if (.not. lsf%unreached(p)) cycle
         write (idx, "(i0)") igrid0 + p - 1
         call abort%latch_message(who//": "//what//" at surface grid point "//trim(idx)// &
                                  " is unreached or not finite", igrid0 + p - 1)
         return
      end do
      ok = .true.
   end function tile_ok

end module moist_cavity_drop_tile
