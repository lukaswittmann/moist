!> Provides the Signed Sphere Distance (SSD) function and its derivatives
!>
!> Screening System:
!>   Also provides screening types for optimizing triplet summations by pre-computing derivatives only
!>   for atoms that contribute significantly (exponential k3f0 > threshold).
!>
!> Layout:
!>   Uses Structure-of-Arrays (SoA) for active atoms. All per-atom derivative data is stored in
!>   contiguous arrays (k3f0_arr, f1_r_arr, f2_rr_arr, etc.)
!>
module moist_cavity_drop_lsf_svdw_ssd
   use mctc_env_accuracy, only: wp
   use moist_math_sorter_counting_sort, only: counting_argsort
   implicit none
   private

   integer, parameter :: ndim = 3

   !> System manager for collection of active atoms (SoA layout)
   ! TODO: Feed in neighbour list into the build/update functions to avoid full rebuilds
   type :: moist_cavity_drop_lsf_svdw_ssd_type
      !> Number of active atoms (valid entries in SoA arrays)
      integer :: n_active = 0
      !> Sharpness parameter k
      real(wp) :: k = 0.0_wp
      !> Screening threshold for exp(-k/3 * d_I)
      real(wp) :: threshold = 0.0_wp
      !> Precomputed distance cutoff: (x_i - radius) > screening_cutoff
      !> implies exp(-k/3 * (x_i - radius)) < threshold, so the atom fails
      !> screening and can be skipped without evaluating exp().
      real(wp) :: screening_cutoff = huge(0.0_wp)
      !> Maximum derivative order computed (0-4)
      integer :: max_deriv = 0
      !> Stored atom centers in user-space (caller) order [ndim, n_atoms].
      !> Kept in caller ordering so that external consumers (e.g.
      !> projector.f90) can index with their atom ids directly.
      real(wp), allocatable :: centers(:, :)
      !> Stored atom radii in user-space (caller) order [n_atoms].
      real(wp), allocatable :: radii(:)
      !> Spatially sorted copy of centers [ndim, n_atoms]. sorted_centers(:, j)
      !> is the center of the atom whose user-space index is sorted_to_orig(j).
      !> Used exclusively inside the compute hot loop so that gathers driven
      !> by cell-local candidate index lists become near-sequential in memory.
      real(wp), allocatable :: sorted_centers(:, :)
      !> Spatially sorted copy of radii [n_atoms]. sorted_radii(j) is the
      !> radius of the atom whose user-space index is sorted_to_orig(j).
      real(wp), allocatable :: sorted_radii(:)
      !> Permutation mapping a user-space atom index i to its internal
      !> spatially-sorted position. Populated in %update; size = n_atoms.
      integer, allocatable :: orig_to_sorted(:)
      !> Inverse of orig_to_sorted. sorted_to_orig(j) gives the user-space
      !> atom index of the atom stored at internal sorted position j.
      integer, allocatable :: sorted_to_orig(:)

      !> --- SoA active-atom derivative storage ---
      !> All arrays pre-allocated to capacity = size(radii) in %update;
      !> only the first n_active entries are valid after each compute call.
      !> Original atom indices in user-space ordering [n_alloc]
      integer, allocatable :: atom_indices(:)
      !> Cached screening weight exp(-(k/3) * f0) [n_alloc]
      real(wp), allocatable :: k3f0_arr(:)
      !> Cached sqrt(k3f0), used by LSF for e2 exponential kind [n_alloc]
      real(wp), allocatable :: sqrt_k3f0_arr(:)
      !> Gradient of f0 w.r.t. field point r [ndim, n_alloc]
      real(wp), allocatable :: f1_r_arr(:, :)
      !> Reciprocal distance 1/x per active atom [n_alloc] (if max_deriv >= 2)
      !> The Hessian f2_rr = (I - n n^T)/x is reconstruted on the fly from f1_r (= n) and inv_x in z012/pou f012
      real(wp), allocatable :: inv_x_arr(:)
      !> Hessian of f0 w.r.t. r [ndim, ndim, n_alloc]; only saved if max_deriv >= 3
      !> (where the higher-order / nuclear routines require it; see inv_x_arr)
      real(wp), allocatable :: f2_rr_arr(:, :, :)
      !> Third derivative w.r.t. r [ndim, ndim, ndim, n_alloc] (max_deriv >= 3)
      real(wp), allocatable :: f3_rrr_arr(:, :, :, :)
      !> Fourth derivative w.r.t. r [ndim, ndim, ndim, ndim, n_alloc] (max_deriv >= 4)
      real(wp), allocatable :: f4_rrrr_arr(:, :, :, :, :)

   contains
      !> Initialize system with parameters
      procedure :: new => new_ssd_system
      !> Update geometry (centers and radii)
      procedure :: update => ssd_system_update
      !> Compute derivatives at point for a caller-provided candidate subset
      procedure :: compute => ssd_system_compute
      !> Get number of active nodes
      procedure :: get_n_active => ssd_system_get_n_active
   end type moist_cavity_drop_lsf_svdw_ssd_type

   ! SSD derivative functions
   public :: ssd0
   public :: ssd1_r, ssd2_rr, ssd3_rrr, ssd4_rrrr
   public :: ssd2_r_rA
   public :: ssd1_rA, ssd2_rArB
   public :: ssd012_r

   ! Screening types
   public :: moist_cavity_drop_lsf_svdw_ssd_type

contains

   !* ================================================================================= *!
   !*                  SSD functions and derivatives with 3d positions                  *!
   !* ================================================================================= *!

   !> Combined computation of value, gradient, and hessian for SSD.
   !> Wrapper routine for backward compatibility that calls individual functions.
   !> @param[in]  point     Evaluation coordinates [ndim]
   !> @param[in]  center    Sphere center coordinates [ndim]
   !> @param[in]  radius    Sphere radius
   !> @param[out] value     SSD value (signed distance)
   !> @param[out] gradient  First derivative w.r.t. point (optional) [ndim]
   !> @param[out] hessian   Second derivative w.r.t. point (optional) [ndim, ndim]
   pure subroutine ssd012_r(point, center, radius, value, gradient, hessian)
      !> Cartesian coordinates
      real(wp), intent(in) :: point(ndim)
      !> Sphere center
      real(wp), intent(in) :: center(ndim)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> SSD value
      real(wp), intent(out) :: value
      !> First derivative
      real(wp), intent(out), optional :: gradient(ndim)
      !> Second derivative
      real(wp), intent(out), optional :: hessian(ndim, ndim)

      value = ssd0(point, center, radius)
      if (present(gradient)) then
         gradient = ssd1_r(point, center, radius)
         if (present(hessian)) then
            hessian = ssd2_rr(point, center, radius)
         end if
      end if

   end subroutine ssd012_r

   !> Signed distance from a point to a sphere surface
   !>
   !> @param[in] point    Evaluation coordinates [ndim]
   !> @param[in] center   Sphere center coordinates [ndim]
   !> @param[in] radius   Sphere radius
   !> @returns   distance Signed distance
   pure function ssd0(point, center, radius) result(distance)
      !> Cartesian coordinates
      real(wp), intent(in) :: point(ndim)
      !> Sphere center
      real(wp), intent(in) :: center(ndim)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Signed distance
      real(wp) :: distance

      distance = norm2(point - center) - radius
   end function ssd0
   !> First derivative of SSD with respect to point coordinates
   !>
   !> @param[in] point    Evaluation coordinates [ndim]
   !> @param[in] center   Sphere center coordinates [ndim]
   !> @param[in] radius   Sphere radius
   !> @returns   gradient First derivative vector [ndim]
   pure function ssd1_r(point, center, radius) result(gradient)
      !> Cartesian coordinates
      real(wp), intent(in) :: point(ndim)
      !> Sphere center
      real(wp), intent(in) :: center(ndim)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Gradient vector
      real(wp) :: gradient(ndim)
      !> Displacement vector and its norm
      real(wp) :: diff(ndim), x

      diff = point - center
      x = norm2(diff)
      if (x > 0.0_wp) then
         gradient = diff/x
      else
         gradient = 0.0_wp
      end if
   end function ssd1_r

   !> Second derivative of SSD with respect to point coordinates
   !>
   !> @param[in] point   Evaluation coordinates [ndim]
   !> @param[in] center  Sphere center coordinates [ndim]
   !> @param[in] radius  Sphere radius
   !> @returns   hessian Second derivative matrix [ndim, ndim]
   pure function ssd2_rr(point, center, radius) result(hessian)
      !> Cartesian coordinates
      real(wp), intent(in) :: point(ndim)
      !> Sphere center
      real(wp), intent(in) :: center(ndim)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Hessian matrix
      real(wp) :: hessian(ndim, ndim)
      !> Displacement, distance, unit vector
      real(wp) :: diff(ndim), x, n(ndim)
      integer :: i, j

      diff = point - center
      x = norm2(diff)
      if (x > 0.0_wp) then
         n = diff/x
         do j = 1, ndim
            do i = 1, ndim
               hessian(i, j) = -n(i)*n(j)/x
            end do
            hessian(j, j) = hessian(j, j) + 1.0_wp/x
         end do
      else
         hessian = 0.0_wp
      end if
   end function ssd2_rr

   !> Third derivative of SSD with respect to point coordinates
   !>
   !> @param[in] point  Evaluation coordinates [ndim]
   !> @param[in] center Sphere center coordinates [ndim]
   !> @param[in] radius Sphere radius
   !> @returns   third  Third derivative tensor [ndim, ndim, ndim]
   pure function ssd3_rrr(point, center, radius) result(third)
      !> Cartesian coordinates
      real(wp), intent(in) :: point(ndim)
      !> Sphere center
      real(wp), intent(in) :: center(ndim)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Third derivative tensor
      real(wp) :: third(ndim, ndim, ndim)
      !> Displacement, distance squared, unit vector
      real(wp) :: diff(ndim), x2, n(ndim)
      real(wp) :: dij, dik, djk
      integer :: i, j, k

      diff = point - center
      x2 = dot_product(diff, diff)
      if (x2 > 0.0_wp) then
         n = diff/sqrt(x2)
         do k = 1, ndim
            do j = 1, ndim
               djk = merge(1.0_wp, 0.0_wp, j == k)
               do i = 1, ndim
                  dij = merge(1.0_wp, 0.0_wp, i == j)
                  dik = merge(1.0_wp, 0.0_wp, i == k)
                  third(i, j, k) = -(dij*n(k) + dik*n(j) + djk*n(i) &
                                     - 3.0_wp*n(i)*n(j)*n(k))/x2
               end do
            end do
         end do
      else
         third = 0.0_wp
      end if
   end function ssd3_rrr

   !> Fourth derivative of SSD with respect to point coordinates.
   !>
   !> @param[in] point  Evaluation coordinates [ndim]
   !> @param[in] center Sphere center coordinates [ndim]
   !> @param[in] radius Sphere radius
   !> @returns   fourth Fourth derivative tensor [ndim, ndim, ndim, ndim]
   pure function ssd4_rrrr(point, center, radius) result(fourth)
      !> Cartesian coordinates
      real(wp), intent(in) :: point(ndim)
      !> Sphere center
      real(wp), intent(in) :: center(ndim)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Fourth derivative tensor
      real(wp) :: fourth(ndim, ndim, ndim, ndim)
      !> Displacement, x^3, unit vector
      real(wp) :: diff(ndim), x, x3, n(ndim)
      real(wp) :: dij, dik, dil, djk, djl, dkl, nn
      integer :: i, j, k, l

      diff = point - center
      x = norm2(diff)
      if (x > 0.0_wp) then
         x3 = x*x*x
         n = diff/x
         do l = 1, ndim
            do k = 1, ndim
               dkl = merge(1.0_wp, 0.0_wp, k == l)
               do j = 1, ndim
                  djk = merge(1.0_wp, 0.0_wp, j == k)
                  djl = merge(1.0_wp, 0.0_wp, j == l)
                  do i = 1, ndim
                     dij = merge(1.0_wp, 0.0_wp, i == j)
                     dik = merge(1.0_wp, 0.0_wp, i == k)
                     dil = merge(1.0_wp, 0.0_wp, i == l)
                     nn = n(i)*n(j)*n(k)*n(l)
                     fourth(i, j, k, l) = ( &
                                          -(dij*dkl + dik*djl + dil*djk) &
                                          + 3.0_wp*(dij*n(k)*n(l) + dik*n(j)*n(l) + dil*n(j)*n(k) &
                                                    + djk*n(i)*n(l) + djl*n(i)*n(k) + dkl*n(i)*n(j)) &
                                          - 15.0_wp*nn)/x3
                  end do
               end do
            end do
         end do
      else
         fourth = 0.0_wp
      end if
   end function ssd4_rrrr

   !> Mixed second derivative of SSD with respect to point and center
   !>
   !> @param[in] point         Evaluation coordinates [ndim]
   !> @param[in] center        Sphere center coordinates [ndim]
   !> @param[in] radius        Sphere radius
   !> @returns   mixed_second  Mixed second derivative [ndim, ndim]
   pure function ssd2_r_rA(point, center, radius) result(mixed_second)
      !> Cartesian coordinates
      real(wp), intent(in) :: point(ndim)
      !> Sphere center
      real(wp), intent(in) :: center(ndim)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Mixed derivative matrix
      real(wp) :: mixed_second(ndim, ndim)

      mixed_second = -ssd2_rr(point, center, radius)
   end function ssd2_r_rA

   !> First derivative of SSD with respect to center (nuclear) coordinates
   !>
   !> @param[in] point            Evaluation coordinates [ndim]
   !> @param[in] center           Sphere center coordinates [ndim]
   !> @param[in] radius           Sphere radius
   !> @returns   nuclear_gradient Nuclear gradient vector [ndim]
   pure function ssd1_rA(point, center, radius) result(nuclear_gradient)
      !> Cartesian coordinates
      real(wp), intent(in) :: point(ndim)
      !> Sphere center
      real(wp), intent(in) :: center(ndim)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Nuclear gradient
      real(wp) :: nuclear_gradient(ndim)

      nuclear_gradient = -ssd1_r(point, center, radius)
   end function ssd1_rA

   !> Second derivative of SSD with respect to center coordinates
   !>
   !> @param[in] point           Evaluation coordinates [ndim]
   !> @param[in] center          Sphere center coordinates [ndim]
   !> @param[in] radius          Sphere radius
   !> @returns   nuclear_hessian Nuclear Hessian matrix [ndim, ndim]
   pure function ssd2_rArB(point, center, radius) result(nuclear_hessian)
      !> Cartesian coordinates
      real(wp), intent(in) :: point(ndim)
      !> Sphere center
      real(wp), intent(in) :: center(ndim)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Nuclear Hessian
      real(wp) :: nuclear_hessian(ndim, ndim)

      nuclear_hessian = ssd2_rr(point, center, radius)
   end function ssd2_rArB

   !* ================================================================================= *!
   !*            SoA population for a single active atom slot                          *!
   !* ================================================================================= *!

   !> Populate SoA slot idx with derivative data for one active atom
   !>
   !> @param[inout] self      System whose SoA arrays to populate
   !> @param[in]    idx       Slot index (1..n_active)
   !> @param[in]    n         Unit direction vector (r-A)/||r-A|| [ndim]
   !> @param[in]    x         Distance ||r - A||
   !> @param[in]    radius    Sphere radius
   !> @param[in]    atom_idx  Original atom index in user-space ordering
   !> @param[in]    ef0       Pre-computed exp(-(k/3) * f0) from screening pass
   pure subroutine ssd_soa_populate(self, idx, n, x, radius, atom_idx, ef0)
      class(moist_cavity_drop_lsf_svdw_ssd_type), intent(inout) :: self
      !> Slot index in SoA arrays
      integer, intent(in) :: idx
      !> Unit direction vector from atom to field point [ndim]
      real(wp), intent(in) :: n(ndim)
      !> Distance ||r - A||
      real(wp), intent(in) :: x
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Original atom index
      integer, intent(in) :: atom_idx
      !> Pre-computed exp(-(k/3) * f0) from the screening pass
      real(wp), intent(in) :: ef0
      !> Reciprocal distance and its powers
      real(wp) :: inv_x, inv_x2, inv_x3
      !> Kronecker delta temporaries
      real(wp) :: dij, dik, dil, djk, djl, dkl
      !> Loop indices (k_ avoids shadowing the parameter k)
      integer :: i, j, k_, l

      self%atom_indices(idx) = atom_idx
      self%k3f0_arr(idx) = ef0
      self%sqrt_k3f0_arr(idx) = sqrt(ef0)

      if (self%max_deriv >= 1) then
         self%f1_r_arr(:, idx) = n
      else
         self%f1_r_arr(:, idx) = 0.0_wp
      end if

      if (self%max_deriv >= 2 .and. x > 0.0_wp) then
         inv_x = 1.0_wp/x
         self%inv_x_arr(idx) = inv_x

         if (self%max_deriv >= 3) then
            ! Materialize the full Hessian
            do j = 1, ndim
               do i = 1, ndim
                  self%f2_rr_arr(i, j, idx) = -n(i)*n(j)*inv_x
               end do
               self%f2_rr_arr(j, j, idx) = self%f2_rr_arr(j, j, idx) + inv_x
            end do

            inv_x2 = inv_x*inv_x
            do k_ = 1, ndim
               do j = 1, ndim
                  djk = merge(1.0_wp, 0.0_wp, j == k_)
                  do i = 1, ndim
                     dij = merge(1.0_wp, 0.0_wp, i == j)
                     dik = merge(1.0_wp, 0.0_wp, i == k_)
                     self%f3_rrr_arr(i, j, k_, idx) = -(dij*n(k_) + dik*n(j) + djk*n(i) &
                                                        - 3.0_wp*n(i)*n(j)*n(k_))*inv_x2
                  end do
               end do
            end do

            if (self%max_deriv >= 4) then
               inv_x3 = inv_x2*inv_x
               do l = 1, ndim
                  do k_ = 1, ndim
                     dkl = merge(1.0_wp, 0.0_wp, k_ == l)
                     do j = 1, ndim
                        djk = merge(1.0_wp, 0.0_wp, j == k_)
                        djl = merge(1.0_wp, 0.0_wp, j == l)
                        do i = 1, ndim
                           dij = merge(1.0_wp, 0.0_wp, i == j)
                           dik = merge(1.0_wp, 0.0_wp, i == k_)
                           dil = merge(1.0_wp, 0.0_wp, i == l)
                           self%f4_rrrr_arr(i, j, k_, l, idx) = ( &
                                                                -(dij*dkl + dik*djl + dil*djk) &
                                                                + 3.0_wp*(dij*n(k_)*n(l) + dik*n(j)*n(l) + dil*n(j)*n(k_) &
                                                                          + djk*n(i)*n(l) + djl*n(i)*n(k_) + dkl*n(i)*n(j)) &
                                                                - 15.0_wp*n(i)*n(j)*n(k_)*n(l))*inv_x3
                        end do
                     end do
                  end do
               end do
            end if
         end if
      else if (self%max_deriv >= 2) then
         ! x == 0: zero second-order field
         self%inv_x_arr(idx) = 0.0_wp
         if (self%max_deriv >= 3) self%f2_rr_arr(:, :, idx) = 0.0_wp
      end if
   end subroutine ssd_soa_populate

   !* ================================================================================= *!
   !*                 Initializers for the whole ssd system (molecule)                 *!
   !* ================================================================================= *!

   !> Initialize system with parameters
   !>
   !> @param[inout] self      System to initialize
   !> @param[in]    k         Sharpness parameter
   !> @param[in]    threshold Screening cutoff for exponential k3f0 (optional, default=tiny)
   !> @param[in]    max_deriv Maximum derivative order (0-4, default=2)
   subroutine new_ssd_system(self, k, threshold, max_deriv)
      class(moist_cavity_drop_lsf_svdw_ssd_type), intent(inout) :: self
      real(wp), intent(in) :: k
      real(wp), intent(in), optional :: threshold
      integer, intent(in), optional :: max_deriv
      real(wp) :: thresh
      integer :: deriv_order

      ! Clean up any existing data
      if (allocated(self%centers)) deallocate (self%centers)
      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%sorted_centers)) deallocate (self%sorted_centers)
      if (allocated(self%sorted_radii)) deallocate (self%sorted_radii)
      if (allocated(self%orig_to_sorted)) deallocate (self%orig_to_sorted)
      if (allocated(self%sorted_to_orig)) deallocate (self%sorted_to_orig)
      if (allocated(self%atom_indices)) deallocate (self%atom_indices)
      if (allocated(self%k3f0_arr)) deallocate (self%k3f0_arr)
      if (allocated(self%sqrt_k3f0_arr)) deallocate (self%sqrt_k3f0_arr)
      if (allocated(self%f1_r_arr)) deallocate (self%f1_r_arr)
      if (allocated(self%inv_x_arr)) deallocate (self%inv_x_arr)
      if (allocated(self%f2_rr_arr)) deallocate (self%f2_rr_arr)
      if (allocated(self%f3_rrr_arr)) deallocate (self%f3_rrr_arr)
      if (allocated(self%f4_rrrr_arr)) deallocate (self%f4_rrrr_arr)
      self%n_active = 0

      ! Store parameters
      self%k = k
      thresh = tiny(0.0_wp)
      if (present(threshold)) thresh = threshold
      self%threshold = thresh
      deriv_order = 2
      if (present(max_deriv)) deriv_order = max_deriv
      self%max_deriv = deriv_order

      ! Precompute screening distance cutoff.
      ! From exp(-(k/3) * (x - R)) >= threshold we have
      ! (x - R) <= -3*log(threshold)/k.
      ! If threshold <= 0 or k <= 0 the check would be ill-defined, so disable it.
      if (self%threshold > 0.0_wp .and. self%k > 0.0_wp) then
         self%screening_cutoff = -3.0_wp*log(self%threshold)/self%k
      else
         self%screening_cutoff = huge(0.0_wp)
      end if
   end subroutine new_ssd_system

   !> Update geometry (centers and radii)
   !>
   !> @param[inout] self    System to update
   !> @param[in]    centers Atom center coordinates [ndim, n_atoms]
   !> @param[in]    radii   Atom radii [n_atoms]
   pure subroutine ssd_system_update(self, centers, radii)
      class(moist_cavity_drop_lsf_svdw_ssd_type), intent(inout) :: self
      real(wp), intent(in) :: centers(:, :)
      real(wp), intent(in) :: radii(:)
      !> Loop index / atom count
      integer :: i, j, n
      !> Per-axis minima / extents used to discretize coordinates for bucket assignment
      real(wp) :: xmin(ndim), xmax(ndim), span(ndim)
      !> Discretized per-axis bucket index
      integer :: bx, by, bz
      !> Buckets per axis for spatial grouping (6^3 = 216 buckets, fits in 0..255)
      integer, parameter :: nbx = 6, nby = 6, nbz = 6
      !> Total number of spatial buckets
      integer, parameter :: n_buckets = nbx*nby*nbz
      !> Per-atom spatial bucket (0..n_buckets-1)
      integer, allocatable :: buckets(:)

      ! Store geometry in user-space order (caller ordering) so that external
      ! consumers that index %centers with user-space atom ids keep working.
      if (allocated(self%centers)) deallocate (self%centers)
      if (allocated(self%radii)) deallocate (self%radii)
      allocate (self%centers, source=centers)
      allocate (self%radii, source=radii)

      n = size(self%radii)

      ! (Re)allocate permutation maps.
      if (allocated(self%orig_to_sorted)) deallocate (self%orig_to_sorted)
      if (allocated(self%sorted_to_orig)) deallocate (self%sorted_to_orig)
      allocate (self%orig_to_sorted(n))
      allocate (self%sorted_to_orig(n))

      ! Spatial bucket sort. Assign each atom to a coarse 3D bucket
      ! (nbx * nby * nbz = 216 buckets) and use a single-pass counting
      ! sort to build the sorted_to_orig permutation. O(N) time, inspired
      ! by stdlib's int8 radix_sort (which degenerates to counting sort).
      if (n <= 0) then
         ! Nothing to sort; skip permutation population.
      else if (n == 1) then
         self%sorted_to_orig(1) = 1
         self%orig_to_sorted(1) = 1
      else
         allocate (buckets(n))
         xmin = minval(centers, dim=2)
         xmax = maxval(centers, dim=2)
         span = xmax - xmin

         ! Assign each atom to a spatial bucket
         do i = 1, n
            if (span(1) > 0.0_wp) then
               bx = min(nbx - 1, int((centers(1, i) - xmin(1))/span(1)*nbx))
            else
               bx = 0
            end if
            if (span(2) > 0.0_wp) then
               by = min(nby - 1, int((centers(2, i) - xmin(2))/span(2)*nby))
            else
               by = 0
            end if
            if (span(3) > 0.0_wp) then
               bz = min(nbz - 1, int((centers(3, i) - xmin(3))/span(3)*nbz))
            else
               bz = 0
            end if
            buckets(i) = bx + by*nbx + bz*nbx*nby
         end do

         ! Single-pass counting-sort argsort (O(N), stable)
         call counting_argsort(buckets, n_buckets - 1, self%sorted_to_orig)
         do j = 1, n
            self%orig_to_sorted(self%sorted_to_orig(j)) = j
         end do
         deallocate (buckets)
      end if

      ! Populate the spatially-sorted mirror of centers/radii used inside the
      ! compute hot loop.
      if (allocated(self%sorted_centers)) deallocate (self%sorted_centers)
      if (allocated(self%sorted_radii)) deallocate (self%sorted_radii)
      allocate (self%sorted_centers(ndim, n))
      allocate (self%sorted_radii(n))
      do j = 1, n
         self%sorted_centers(:, j) = centers(:, self%sorted_to_orig(j))
         self%sorted_radii(j) = radii(self%sorted_to_orig(j))
      end do

      ! (Re)allocate SoA arrays sized to the full atom count
      ! Only the first self%n_active entries are valid after each
      ! compute call
      ! sequential access / BLAS for higher-order derivatives.
      if (allocated(self%atom_indices)) deallocate (self%atom_indices)
      if (allocated(self%k3f0_arr)) deallocate (self%k3f0_arr)
      if (allocated(self%sqrt_k3f0_arr)) deallocate (self%sqrt_k3f0_arr)
      if (allocated(self%f1_r_arr)) deallocate (self%f1_r_arr)
      if (allocated(self%inv_x_arr)) deallocate (self%inv_x_arr)
      if (allocated(self%f2_rr_arr)) deallocate (self%f2_rr_arr)
      if (allocated(self%f3_rrr_arr)) deallocate (self%f3_rrr_arr)
      if (allocated(self%f4_rrrr_arr)) deallocate (self%f4_rrrr_arr)
      allocate (self%atom_indices(n))
      allocate (self%k3f0_arr(n))
      allocate (self%sqrt_k3f0_arr(n))
      ! f1_r and inv_x are always allocated; z012_rr / pou f012 reconstruct the Hessian (I - n n^T)/x on the fly
      allocate (self%f1_r_arr(ndim, n))
      allocate (self%inv_x_arr(n))
      if (self%max_deriv >= 3) allocate (self%f2_rr_arr(ndim, ndim, n))
      if (self%max_deriv >= 3) allocate (self%f3_rrr_arr(ndim, ndim, ndim, n))
      if (self%max_deriv >= 4) allocate (self%f4_rrrr_arr(ndim, ndim, ndim, ndim, n))
   end subroutine ssd_system_update

   !> Compute derivatives at point
   !>
   !> Derivative orders:
   !>  - 0: Only f0, k3f0 (for value-only computations)
   !>  - 1: + f1_r (for gradients)
   !>  - 2: + inv_x (for Hessians; the full f2_rr is reconstructed on the fly
   !>       in the consumer from f1_r and inv_x) [DEFAULT]
   !>  - 3: + f2_rr, f3_rrr (for 3rd derivatives, nuclear Hessians)
   !>  - 4: + f4_rrrr (for 4th derivatives, nuclear 3rd derivatives)
   !>
   !> @param[inout] self              System to compute
   !> @param[in]    point             Evaluation point coordinates [ndim]
   !> @param[in]    candidate_indices Sorted-index atom ids to screen
   subroutine ssd_system_compute(self, point, candidate_indices)
      class(moist_cavity_drop_lsf_svdw_ssd_type), intent(inout) :: self
      real(wp), intent(in) :: point(ndim)
      integer, intent(in) :: candidate_indices(:)
      integer :: n_total, n_active_local

      self%n_active = 0
      if (.not. allocated(self%sorted_centers)) return

      n_total = size(candidate_indices)
      if (n_total == 0) return

      call ssd_screen_subset(self, point, candidate_indices, n_total, n_active_local)
      self%n_active = n_active_local
   end subroutine ssd_system_compute

   !> Screen a caller-provided subset of atoms and populate SoA arrays
   !>
   !> @param[inout] self              System
   !> @param[in]    point             Evaluation point [ndim]
   !> @param[in]    candidate_indices Spatially-sorted atom ids to screen
   !>                                 (already in `sorted_*` index space; the
   !>                                 caller remaps via `remap_candidate_grid`)
   !> @param[in]    n_total           Size of candidate_indices
   !> @param[out]   n_active_local    Number of atoms that passed screening
   pure subroutine ssd_screen_subset(self, point, candidate_indices, n_total, n_active_local)
      class(moist_cavity_drop_lsf_svdw_ssd_type), intent(inout) :: self
      real(wp), intent(in) :: point(ndim)
      integer, intent(in) :: candidate_indices(:)
      integer, intent(in) :: n_total
      integer, intent(out) :: n_active_local
      integer :: i, sorted_idx, user_idx
      real(wp) :: x_i, w_i, diff(ndim), n_vec(ndim), bound

      n_active_local = 0
      do i = 1, n_total
         ! Candidates already carry sorted-index ids
         sorted_idx = candidate_indices(i)

         diff = point - self%sorted_centers(:, sorted_idx)
         ! Conservative squared-distance pre-filter
         bound = self%sorted_radii(sorted_idx) + self%screening_cutoff
         if (bound > 0.0_wp) then
            if (sum(diff*diff) > bound*bound*(1.0_wp + 1.0e-12_wp)) cycle
         end if
         x_i = norm2(diff)
         if (x_i - self%sorted_radii(sorted_idx) > self%screening_cutoff) cycle
         w_i = exp(-(self%k/3.0_wp)*(x_i - self%sorted_radii(sorted_idx)))

         if (w_i >= self%threshold) then
            n_active_local = n_active_local + 1
            if (x_i > 0.0_wp) then
               n_vec = diff/x_i
            else
               n_vec = 0.0_wp
            end if
            ! Translate back to user space only for the (far fewer) kept atoms
            user_idx = self%sorted_to_orig(sorted_idx)
            call ssd_soa_populate(self, n_active_local, &
                                  n_vec, x_i, self%radii(user_idx), user_idx, w_i)
         end if
      end do
   end subroutine ssd_screen_subset

   !> Get number of active nodes
   !>
   !> @param[in] self      System instance
   !> @returns   n_active  Number of active atoms
   pure function ssd_system_get_n_active(self) result(n_active)
      class(moist_cavity_drop_lsf_svdw_ssd_type), intent(in) :: self
      integer :: n_active

      n_active = self%n_active
   end function ssd_system_get_n_active

end module moist_cavity_drop_lsf_svdw_ssd
