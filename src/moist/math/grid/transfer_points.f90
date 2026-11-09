!> Interpolation from scattered source points by moving least squares
!>
!> - A target fits an affine function of [1, dx, dy, dz] to its k nearest source points,
!>   weighted by the Wendland C2 function of radius `radius_factor` times the distance of
!>   the k-th point; the fit, solved by SVD of the weighted design matrix, gives a
!>   stencil: the value is sum_j a_j f_j over the neighbours
!> - k starts at `k_first` and doubles up to `k_last` while a stencil is not accepted
!> - Accepted: the nearest source point lies within its reach (`reach_factor` times its
!>   distance to its nearest distinct source point), the singular values satisfy
!>   s_min >= `min_conditioning` s_max, and the amplification sum_j |a_j| is at most
!>   `max_amplification`; this rejects targets beyond the reach and coplanar, collinear or
!>   coincident neighbourhoods
!> - Extrapolation is bounded, not excluded: a target within the reach but outside its
!>   source points is extrapolated linearly while the amplification stays within its bound
!> - Affine fields are reproduced at every accepted target up to roundoff; quadrature
!>   weights play no part
!> - Neighbours are ordered by distance, then index, so stencils do not depend on threads
module moist_math_grid_transfer_points
   use mctc_env, only: wp
   use moist_math_kdtree, only: moist_kdtree_type
   use moist_math_lapack_gesvd, only: lapack_gesvd
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none(type, external)
   private

   public :: point_stencil, point_reach
   public :: k_first, k_last
   public :: stencil_fitted, stencil_far, stencil_conditioning, stencil_amplification

   !> Neighbours of the first fit and the most neighbours of a fit
   integer, parameter :: k_first = 16, k_last = 64
   !> Weight radius in units of the distance of the k-th neighbour
   real(wp), parameter :: radius_factor = 1.1_wp
   !> Smallest ratio of the smallest to the largest singular value
   real(wp), parameter :: min_conditioning = 1.0e-3_wp
   !> Largest sum of the absolute stencil weights
   real(wp), parameter :: max_amplification = 4.0_wp
   !> Reach of a source point in units of the distance to its nearest source point
   real(wp), parameter :: reach_factor = 4.0_wp

   !> Outcome of a stencil
   integer, parameter :: stencil_fitted = 0, stencil_far = 1, stencil_conditioning = 2, stencil_amplification = 3

   !> Work space of the SVD
   integer, parameter :: lwork = 512

contains

   !> Reach of each source point: `reach_factor` times the distance to its nearest distinct
   !> source point, farther than `floor`; at least `floor`
   !>
   !> - Coinciding source points do not shrink the reach of their neighbourhood; a point
   !>   whose `k_last` nearest points all coincide with it keeps the reach `floor`
   !>
   !> @param[in]  tree      k-d tree of the source points
   !> @param[in]  floor     smallest reach and largest distance of coinciding points, bohr
   !> @param[out] reach     reach of each source point, bohr (n)
   !> @param[in]  nthreads  threads
   subroutine point_reach(tree, floor, reach, nthreads)
      !> Tree of the source points
      type(moist_kdtree_type), intent(in) :: tree
      !> Smallest reach
      real(wp), intent(in) :: floor
      !> Reach of each source point
      real(wp), allocatable, intent(out) :: reach(:)
      !> Threads
      integer, intent(in) :: nthreads
      integer :: nearest(k_last), j, k, m, n
      real(wp) :: dist2(k_last)

      n = size(tree%x, 2)
      allocate (reach(n), source=floor)
      if (n < 2) return
      !$omp parallel do schedule(dynamic, 256) num_threads(nthreads) private(nearest, dist2, k, m)
      do j = 1, n
         k = 2
         search: do
            call tree%knn(tree%x(:, j), k, nearest, dist2)
            do m = 2, k
               if (dist2(m) > floor**2) then
                  reach(j) = max(floor, reach_factor*sqrt(dist2(m)))
                  exit search
               end if
            end do
            if (k >= min(k_last, n)) exit search
            k = min(2*k, k_last, n)
         end do search
      end do
   end subroutine point_reach

   !> Stencil of one target from its nearest source points
   !>
   !> @param[in]  tree     k-d tree of the source points
   !> @param[in]  reach    reach of each source point, bohr
   !> @param[in]  y        target in the source frame, bohr
   !> @param[out] nearest  source points of the stencil (at least `k_last`)
   !> @param[out] weight   stencil weights (at least `k_last`)
   !> @param[out] k        points of the stencil; 0 without one
   !> @param[out] status   `stencil_fitted`, or why there is no stencil
   pure subroutine point_stencil(tree, reach, y, nearest, weight, k, status)
      !> Tree of the source points
      type(moist_kdtree_type), intent(in) :: tree
      !> Reach of each source point
      real(wp), intent(in) :: reach(:)
      !> Target
      real(wp), intent(in) :: y(3)
      !> Source points of the stencil
      integer, intent(out) :: nearest(:)
      !> Stencil weights
      real(wp), intent(out) :: weight(:)
      !> Points of the stencil
      integer, intent(out) :: k
      !> Outcome
      integer, intent(out) :: status
      real(wp) :: dist2(k_last), b(k_last, 4), u(k_last, 4), s(4), vt(4, 4), work(lwork), sw(k_last), d(3)
      real(wp) :: radius, r, amplification
      integer(lapack_ik) :: info
      integer :: kk, navail, j, m

      k = 0
      weight = 0.0_wp
      status = stencil_far
      navail = size(tree%x, 2)
      kk = min(k_first, navail)
      do while (kk >= 1)
         call tree%knn(y, kk, nearest, dist2)
         if (sqrt(dist2(1)) > reach(nearest(1))) then
            status = stencil_far
            return
         end if
         status = stencil_conditioning
         if (kk < 4) return
         if (.not. (dist2(kk) > 0.0_wp)) then
            ! Every neighbour coincides with the target: no spread to fit; more neighbours
            if (kk >= min(k_last, navail)) exit
            kk = min(2*kk, k_last, navail)
            cycle
         end if
         radius = radius_factor*sqrt(dist2(kk))
         do j = 1, kk
            r = sqrt(dist2(j))/radius
            sw(j) = sqrt((1.0_wp - r)**4*(4.0_wp*r + 1.0_wp))
            d = (tree%x(:, nearest(j)) - y)/radius
            b(j, :) = sw(j)*[1.0_wp, d(1), d(2), d(3)]
         end do
         call lapack_gesvd("S", "A", int(kk, lapack_ik), 4_lapack_ik, b, int(k_last, lapack_ik), s, u, &
            & int(k_last, lapack_ik), vt, 4_lapack_ik, work, int(lwork, lapack_ik), info)
         if (info == 0 .and. s(4) >= min_conditioning*s(1) .and. s(1) > 0.0_wp) then
            ! Weights of the constant term of the fit: a = W^(1/2) U S^-1 V^T e_1
            do j = 1, kk
               weight(j) = 0.0_wp
               do m = 1, 4
                  weight(j) = weight(j) + vt(m, 1)*u(j, m)/s(m)
               end do
               weight(j) = sw(j)*weight(j)
            end do
            amplification = sum(abs(weight(1:kk)))
            status = stencil_amplification
            if (amplification <= max_amplification) then
               status = stencil_fitted
               k = kk
               return
            end if
         end if
         ! A larger neighbourhood, while there are more points
         if (kk >= min(k_last, navail)) exit
         kk = min(2*kk, k_last, navail)
      end do
      weight = 0.0_wp
   end subroutine point_stencil

end module moist_math_grid_transfer_points
