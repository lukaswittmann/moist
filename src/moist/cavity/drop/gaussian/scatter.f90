!> iSwiG influence-set block scatter
!>
!> [[iswig_swi_f2_rArB_block]] returns one point's second-derivative block over
!> its influence set - a square `(3, n, 3, n)` with `n = work%n_nb + 1` and
!> `idx` **owner first** - while the buffers a derivative pass reduces into are
!> indexed by the whole molecule or by a resolved pair table. This module is
!> that bridge, and it lives beside the producer because what it encodes is the
!> producer's slot order.
!>
!> It is a module of its own so the zero-cost contract stays auditable: nothing
!> allocates, no local array is declared, and every routine is `pure`, so it may
!> be called on a thread-local slice inside an `!$omp` region.
!>
!> Both scatters accumulate, and both address each element under its own scalar
!> index: `acc(:, idx(1:n))` would be a many-one section as soon as an influence
!> set repeats an id, and the repeated writes would be lost.
module moist_cavity_drop_gaussian_scatter
   use mctc_env_accuracy, only: wp

   implicit none(type, external)
   private

   public :: scatter_iswig_block, scatter_iswig_block_indexed

   !> Cartesian dimension
   integer, parameter :: ndim = 3

contains

   !* ================================================================================= *!
   !*                            Influence-set block scatter                            *!
   !* ================================================================================= *!

   !> Scatter one point's iSwiG second-derivative block into a dense Hessian
   !>
   !> `blk(i, ia, j, ib)` is the second derivative w.r.t. `idx(ia)` and
   !> `idx(ib)`, `idx` owner first; only `(:, 1:n, :, 1:n)` is read. Accumulated
   !> under scalar ids, so an influence set naming one atom twice sums all four
   !> of its `(ia, ib)` contributions instead of overwriting three.
   !>
   !> [[scatter_iswig_block_indexed]] mirrors this loop nest and the two are
   !> required to agree bit for bit; keep them in step.
   !>
   !> @param[in]    n      Influence-set size, `work%n_nb + 1`
   !> @param[in]    idx    Atom ids of the influence set, owner first; only `1:n` read
   !> @param[in]    blk    Local second-derivative block; only `(:, 1:n, :, 1:n)` read
   !> @param[in]    weight Scalar the whole point is weighted by
   !> @param[inout] acc    Hessian accumulator (3, nsph, 3, nsph)
   pure subroutine scatter_iswig_block(n, idx, blk, weight, acc)
      !> Influence-set size
      integer, intent(in) :: n
      !> Atom ids of the influence set, owner first; only `1:n` is read
      integer, intent(in) :: idx(:)
      !> Local second-derivative block (3, >= n, 3, >= n)
      real(wp), intent(in) :: blk(:, :, :, :)
      !> Scalar the whole point is weighted by
      real(wp), intent(in) :: weight
      !> Hessian accumulator (3, nsph, 3, nsph)
      real(wp), intent(inout) :: acc(:, :, :, :)

      integer :: ia, ib, katom, latom, jaxis

      do ib = 1, n
         latom = idx(ib)
         do jaxis = 1, ndim
            do ia = 1, n
               katom = idx(ia)
               acc(:, katom, jaxis, latom) = acc(:, katom, jaxis, latom) &
                                             + weight*blk(:, ia, jaxis, ib)
            end do
         end do
      end do

   end subroutine scatter_iswig_block

   !> Scatter the local iSwiG block into a pair-indexed `(3, 3)`-block accumulator
   !>
   !> For a caller whose Hessian is not a molecule-sized rank-4 slab: instead of
   !> atom ids it hands in `ent(ia, ib)`, the index of the `(3, 3)` block that
   !> `(idx(ia), idx(ib))` accumulates into. Resolving the pairs is the caller's
   !> business because this module must not know what the container is.
   !>
   !> The loop nest is identical to [[scatter_iswig_block]]'s so the two agree to
   !> the bit on the same input.
   !>
   !> @param[in]    n      Influence-set size, `work%n_nb + 1`
   !> @param[in]    ent    Entry index of each pair; only `(1:n, 1:n)` read
   !> @param[in]    blk    Local second-derivative block; only `(:, 1:n, :, 1:n)` read
   !> @param[in]    weight Scalar the whole point is weighted by
   !> @param[inout] blocks Pair-indexed accumulator (3, 3, nentries)
   pure subroutine scatter_iswig_block_indexed(n, ent, blk, weight, blocks)
      !> Influence-set size
      integer, intent(in) :: n
      !> Entry index of each influence-set pair; only (1:n, 1:n) is read
      integer, intent(in) :: ent(:, :)
      !> Local second-derivative block (3, >= n, 3, >= n)
      real(wp), intent(in) :: blk(:, :, :, :)
      !> Scalar the whole point is weighted by
      real(wp), intent(in) :: weight
      !> Pair-indexed accumulator (3, 3, nentries)
      real(wp), intent(inout) :: blocks(:, :, :)

      integer :: ia, ib, ient, jaxis

      do ib = 1, n
         do jaxis = 1, ndim
            do ia = 1, n
               ient = ent(ia, ib)
               blocks(:, jaxis, ient) = blocks(:, jaxis, ient) &
                                        + weight*blk(:, ia, jaxis, ib)
            end do
         end do
      end do

   end subroutine scatter_iswig_block_indexed

end module moist_cavity_drop_gaussian_scatter
