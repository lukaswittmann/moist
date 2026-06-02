!> Outer-product constructors, from rank-2 dyadics through rank-4 tensors.
!>
!> Direct small-size implementations that avoid BLAS/temporary overhead.
module moist_math_linalg_outer
   use mctc_env_accuracy, only: wp
   implicit none
   private

   public :: outer_matrix
   public :: outer3
   public :: outer3_linear
   public :: outer4

contains

   !> Compute outer product of two vectors
   !> Computes the dyadic (outer) product:
   !> $$
   !> M_{ij} = \ell_i \, r_j
   !> $$
   !> Direct implementation is used for 3x3 matrices as it is faster than
   !> BLAS dger due to lower function call overhead at this scale.
   !> @param[in] left   Left vector [3]
   !> @param[in] right  Right vector [3]
   !> @returns   mat    Outer product matrix [3, 3]
   pure function outer_matrix(left, right) result(mat)
      !> Left vector
      real(wp), intent(in) :: left(3)
      !> Right vector
      real(wp), intent(in) :: right(3)
      !> Outer product matrix
      real(wp) :: mat(3, 3)

      mat(:, 1) = left*right(1)
      mat(:, 2) = left*right(2)
      mat(:, 3) = left*right(3)
   end function outer_matrix

   !> Add a rank-1 outer product to a 3x3 matrix.
   !>
   !> Computes the in-place update
   !> \f[
   !>   A \leftarrow A + l\,r^{T},
   !> \f]
   !> where \c l = \c left and \c r = \c right are 3-vectors and \c A is a 3x3 matrix.
   !> This is equivalent to (but avoids forming the temporary)
   !> \code
   !>   A = A + outer_matrix(left, right)
   !> \endcode
   !>
   !> Notes:
   !> - Implemented via three column updates \c A(:,j) += left * right(j).
   !> - Avoids allocating/storing a temporary 3x3 matrix; typically faster in tight loops.
   !>
   !> \param[in,out] A     3x3 matrix updated in-place.
   !> \param[in]     left  Left vector (length 3).
   !> \param[in]     right Right vector (length 3).
   pure subroutine outer_add(A, left, right)
      real(wp), intent(inout) :: A(3, 3)
      real(wp), intent(in)    :: left(3), right(3)

      A(:, 1) = A(:, 1) + left*right(1)
      A(:, 2) = A(:, 2) + left*right(2)
      A(:, 3) = A(:, 3) + left*right(3)
   end subroutine outer_add

   !> Compute triple outer product of a vector with itself
   !> Computes the rank-3 tensor:
   !> $$
   !> T_{ijk} = v_i \, v_j \, v_k
   !> $$
   !> This creates a rank-3 tensor representing the third-order self outer product.
   !> Used in fourth derivatives and tensor contractions.
   !> @param[in] vec    Input vector [3]
   !> @returns   tensor Rank-3 tensor [3, 3, 3]
   pure function outer3(vec) result(tensor)
      !> Input vector
      real(wp), intent(in) :: vec(3)
      !> Triple outer product tensor
      real(wp) :: tensor(3, 3, 3)
      integer :: i, j, k

      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               tensor(i, j, k) = vec(i)*vec(j)*vec(k)
            end do
         end do
      end do

   end function outer3

   !> Compute derivative of triple outer product
   !> Computes the linearization of the triple outer product:
   !> $$
   !> \frac{d}{dt}\left[\mathbf{v}(t) \otimes \mathbf{v}(t) \otimes \mathbf{v}(t)\right]
   !> = \delta\mathbf{v} \otimes \mathbf{v} \otimes \mathbf{v}
   !> + \mathbf{v} \otimes \delta\mathbf{v} \otimes \mathbf{v}
   !> + \mathbf{v} \otimes \mathbf{v} \otimes \delta\mathbf{v}
   !> $$
   !> In index notation:
   !> $$
   !> T_{ijk} = \delta v_i \, v_j \, v_k + v_i \, \delta v_j \, v_k + v_i \, v_j \, \delta v_k
   !> $$
   !> This appears in sensitivities of fourth derivatives.
   !> @param[in] vec    Base vector [3]
   !> @param[in] dvec   Derivative/perturbation of vector [3]
   !> @returns   tensor Linearized rank-3 tensor [3, 3, 3]
   pure function outer3_linear(vec, dvec) result(tensor)
      !> Base vector
      real(wp), intent(in) :: vec(3)
      !> Derivative vector
      real(wp), intent(in) :: dvec(3)
      !> Linearized tensor
      real(wp) :: tensor(3, 3, 3)
      integer :: i, j, k

      tensor = 0.0_wp
      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               tensor(i, j, k) = dvec(i)*vec(j)*vec(k) &
                                 + vec(i)*dvec(j)*vec(k) &
                                 + vec(i)*vec(j)*dvec(k)
            end do
         end do
      end do

   end function outer3_linear

   !> Self outer product of rank 4: T_{ijkl} = v_i v_j v_k v_l
   !> @param[in] v   Input vector [3]
   !> @returns   t   Rank-4 tensor [3, 3, 3, 3]
   pure function outer4(v) result(t)
      !> Input vector
      real(wp), intent(in) :: v(3)
      !> Rank-4 self outer product
      real(wp) :: t(3, 3, 3, 3)
      integer :: i, j, k, l
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               do i = 1, 3
                  t(i, j, k, l) = v(i)*v(j)*v(k)*v(l)
               end do
            end do
         end do
      end do
   end function outer4

end module moist_math_linalg_outer
