!> Symmetric tensor combinations used in third- and fourth-derivative assembly.
!>
!> Symmetrized contractions of Hessians, gradients and rank-3 tensors into
!> rank-3 and rank-4 tensors.
module moist_math_linalg_symmetrize
   use mctc_env_accuracy, only: wp
   implicit none
   private

   public :: sym3_21
   public :: sym4_31
   public :: sym4_22
   public :: sym4_211

contains

   !> Symmetrically combine Hessian and gradient to form a rank-3 tensor
   !> Computes the symmetric tensor contraction:
   !> $$
   !> T_{ijk} = H_{ij} g_k + H_{ik} g_j + H_{jk} g_i
   !> $$
   !> This operation appears in third derivatives where the Hessian is contracted
   !> with a gradient vector in a symmetric way.
   !> @param[in] hess   Hessian matrix [3, 3]
   !> @param[in] grad   Gradient vector [3]
   !> @returns   tensor Rank-3 tensor [3, 3, 3]
   pure function sym3_21(hess, grad) result(tensor)
      !> Hessian matrix
      real(wp), intent(in) :: hess(3, 3)
      !> Gradient vector
      real(wp), intent(in) :: grad(3)
      !> Resulting rank-3 tensor
      real(wp) :: tensor(3, 3, 3)
      !> Loop indices
      integer :: i, j, k

      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               tensor(i, j, k) = hess(i, j)*grad(k) + hess(i, k)*grad(j) + hess(j, k)*grad(i)
            end do
         end do
      end do

   end function sym3_21

   !> Symmetrized rank-4 tensor: 4-term sym {grad x third}
   !> T_{ijkl} = g_i T_{jkl} + g_j T_{ikl} + g_k T_{ijl} + g_l T_{ijk}
   !> @param[in] g   Gradient vector [3]
   !> @param[in] h3  Rank-3 tensor [3, 3, 3]
   !> @returns   t   Rank-4 tensor [3, 3, 3, 3]
   pure function sym4_31(g, h3) result(t)
      !> Gradient vector
      real(wp), intent(in) :: g(3)
      !> Rank-3 tensor
      real(wp), intent(in) :: h3(3, 3, 3)
      !> Symmetrized rank-4 tensor
      real(wp) :: t(3, 3, 3, 3)
      integer :: i, j, k, l
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               do i = 1, 3
                  t(i, j, k, l) = g(i)*h3(j, k, l) + g(j)*h3(i, k, l) &
                                  + g(k)*h3(i, j, l) + g(l)*h3(i, j, k)
               end do
            end do
         end do
      end do
   end function sym4_31

   !> Symmetrized rank-4 tensor: 3 unordered pair-partitions of 4 indices
   !> T_{ijkl} = H_{ij} K_{kl} + H_{ik} K_{jl} + H_{il} K_{jk}
   !>
   !> NOTE: This produces only 3 of the 6 ordered 2+2 pair-partition terms
   !> (those with index i on the left factor). For the full 6-term symmetric
   !> version, call twice with swapped arguments: sym4_22(A,B) + sym4_22(B,A).
   !> @param[in] h2a  First Hessian [3, 3]
   !> @param[in] h2b  Second Hessian [3, 3]
   !> @returns   t    Rank-4 tensor [3, 3, 3, 3]
   pure function sym4_22(h2a, h2b) result(t)
      !> First Hessian matrix
      real(wp), intent(in) :: h2a(3, 3)
      !> Second Hessian matrix
      real(wp), intent(in) :: h2b(3, 3)
      !> Symmetrized rank-4 tensor
      real(wp) :: t(3, 3, 3, 3)
      integer :: i, j, k, l
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               do i = 1, 3
                  t(i, j, k, l) = h2a(i, j)*h2b(k, l) &
                                  + h2a(i, k)*h2b(j, l) &
                                  + h2a(i, l)*h2b(j, k)
               end do
            end do
         end do
      end do
   end function sym4_22

   !> Symmetrized rank-4 tensor: 6-term sym {(g g) x hess}
   !> Sum over all 6 unordered pairs {a,b} chosen from {i,j,k,l}:
   !>   g_a g_b * H_{cd}, where {c,d} = complement
   !> @param[in] g   Gradient vector [3]
   !> @param[in] h2  Hessian matrix [3, 3]
   !> @returns   t   Rank-4 tensor [3, 3, 3, 3]
   pure function sym4_211(g, h2) result(t)
      !> Gradient vector
      real(wp), intent(in) :: g(3)
      !> Hessian matrix
      real(wp), intent(in) :: h2(3, 3)
      !> Symmetrized rank-4 tensor
      real(wp) :: t(3, 3, 3, 3)
      integer :: i, j, k, l
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               do i = 1, 3
                  t(i, j, k, l) = g(i)*g(j)*h2(k, l) &
                                  + g(i)*g(k)*h2(j, l) &
                                  + g(i)*g(l)*h2(j, k) &
                                  + g(j)*g(k)*h2(i, l) &
                                  + g(j)*g(l)*h2(i, k) &
                                  + g(k)*g(l)*h2(i, j)
               end do
            end do
         end do
      end do
   end function sym4_211

end module moist_math_linalg_symmetrize
