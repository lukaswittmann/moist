!> Sparse iSwiG row and block scatter for the DROP Hessian
!>
!> The iSwiG switching factor of one surface point depends on the owner atom
!> and on the neighbours [[iswig_swi_collect]] cached for it, and on nothing
!> else. Its derivative rows are therefore returned in a compact index space of
!> size `work%n_nb + 1`, while the buffers a derivative pass reduces into are
!> indexed by the whole molecule. This module is the bridge, and the reason it
!> is a module of its own is that the bridge must not cost anything:
!>
!> * nothing here allocates, and no local array is declared at all, so no
!>   temporary can scale with `nsph`;
!> * the compact rows are read where they are and written straight into the
!>   caller's accumulator, never through a dense `(3, nsph)` intermediate;
!> * every routine is `pure` and holds no module state, so it may be called
!>   from inside an `!$omp` region on a thread-local slice.
!>
!> The two producers in [[moist_cavity_drop_gaussian]] use *different* index
!> conventions, and each scatter here matches one of them:
!>
!> * [[iswig_swi_f1_rA_sparse]] and [[iswig_swi_f2_rArB_sparse]] write
!>   `rows(:, jj)` for the neighbour `work%idx(jj)`, `jj = 1 .. work%n_nb`, and
!>   deliver the owner separately in `owner_row`. The owner is *not* a member of
!>   `work%idx`. [[scatter_iswig_rows]] consumes that pair.
!> * [[iswig_swi_f2_rArB_block]] writes a square `(3, n, 3, n)` block over the
!>   influence set with `n = work%n_nb + 1` and `idx` **owner first**, so the
!>   owner is slot 1 of the same array rather than a separate argument.
!>   [[scatter_iswig_block]] consumes that.
!>
!> Both scatters accumulate. A pass reduces many grid points into one buffer, so
!> assignment would drop every contribution but the last, and an atom that
!> appears twice in one influence set must sum rather than overwrite. That is
!> also why neither routine ever writes through a vector subscript: `acc(:,
!> idx(1:n))` is a many-one section as soon as `idx` repeats an id, and the
!> repeated writes would be lost.
module moist_cavity_drop_derivatives_iswig_scatter
   use mctc_env_accuracy, only: wp

   implicit none(type, external)
   private

   public :: scatter_iswig_rows, scatter_iswig_block

   !> Cartesian dimension of a scatter row
   integer, parameter :: ndim = 3

contains

   !* ================================================================================= *!
   !*                          Compact row scatter                                      *!
   !* ================================================================================= *!

   !> Scatter the compact iSwiG derivative rows of one point into a gradient buffer
   !>
   !> Consumes the `(rows, owner_row)` pair of either sparse producer -- the two
   !> share a layout, so this serves the first-order rows of
   !> [[iswig_swi_f1_rA_sparse]] weighted by an adjoint and the directional
   !> second-order rows of [[iswig_swi_f2_rArB_sparse]] weighted by the same
   !> adjoint alike. Only `1:n_nb` of `rows` and `nb_idx` is read, because the
   !> producers size those from the workspace capacity and write no further.
   !>
   !> `acc` is accumulated into, never assigned, and each row is written under
   !> its own scalar index, so a repeated id in `nb_idx` sums correctly.
   !>
   !> @param[in]    n_nb      Number of cached neighbours, `work%n_nb`
   !> @param[in]    nb_idx    Neighbour atom ids, `work%idx`; only `1:n_nb` read
   !> @param[in]    rows      Derivative rows of the neighbours; only `1:n_nb` read
   !> @param[in]    owner     Owner atom id, `work%owner`
   !> @param[in]    owner_row Derivative row of the owner atom
   !> @param[in]    weight    Scalar the whole point is weighted by
   !> @param[inout] acc       Gradient accumulator (3, nsph)
   pure subroutine scatter_iswig_rows(n_nb, nb_idx, rows, owner, owner_row, weight, acc)
      !> Number of cached neighbours
      integer, intent(in) :: n_nb
      !> Neighbour atom ids; only `1:n_nb` is read
      integer, intent(in) :: nb_idx(:)
      !> Derivative rows of the neighbours (3, >= n_nb)
      real(wp), intent(in) :: rows(:, :)
      !> Owner atom id
      integer, intent(in) :: owner
      !> Derivative row of the owner atom
      real(wp), intent(in) :: owner_row(ndim)
      !> Scalar the whole point is weighted by
      real(wp), intent(in) :: weight
      !> Gradient accumulator (3, nsph)
      real(wp), intent(inout) :: acc(:, :)

      integer :: jj, katom

      do jj = 1, n_nb
         katom = nb_idx(jj)
         acc(:, katom) = acc(:, katom) + weight*rows(:, jj)
      end do

      acc(:, owner) = acc(:, owner) + weight*owner_row(:)

   end subroutine scatter_iswig_rows

   !* ================================================================================= *!
   !*                          Influence-set block scatter                               *!
   !* ================================================================================= *!

   !> Scatter the local iSwiG second-derivative block of one point into a Hessian
   !>
   !> Consumes the `(n, idx, blk)` triple of [[iswig_swi_f2_rArB_block]]: `idx`
   !> is owner first, so the owner needs no separate path here, and `blk(i, ia,
   !> j, ib)` is the second derivative with respect to `idx(ia)` and `idx(ib)`.
   !> Only `(:, 1:n, :, 1:n)` is read, matching what the producer writes.
   !>
   !> This is the direction-free half of the Hessian: the block is assembled
   !> once per point and lands in the assembled matrix in `9 n^2` scattered
   !> writes, with no loop over nuclear directions and no dense intermediate.
   !> The loop order runs the influence-set column outermost and the length-3
   !> Cartesian row innermost, so both arrays are walked with unit stride.
   !>
   !> `acc` is accumulated into and every element is addressed under its own
   !> scalar pair of ids, so an influence set that names the same atom twice
   !> sums all four of its `(ia, ib)` contributions into the one diagonal block
   !> instead of overwriting three of them.
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

end module moist_cavity_drop_derivatives_iswig_scatter
