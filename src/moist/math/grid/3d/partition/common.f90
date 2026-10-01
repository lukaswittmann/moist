!> Batched normalized pair-product partitions
module moist_math_grid_3d_partition_common
   use mctc_env, only: wp
   use moist_data_atomicrad, only: covalent_rad
!$ use omp_lib, only: omp_get_max_threads, omp_in_parallel
   implicit none(type, external)
   private

   public :: pair_partition_weights

   !> Callback contract for pair-cell switches
   abstract interface
      !> Pair-cell weight, decreasing from one to zero
      !>
      !> @param[in] x Pair coordinate
      pure function partition_cell(x) result(s)
         import :: wp
         implicit none(type, external)
         !> Pair coordinate
         real(wp), intent(in) :: x
         !> Weight in [0, 1]
         real(wp) :: s
      end function partition_cell
   end interface

contains

   !> Owner weights from pair cells, normalized in log space to avoid product underflow
   !>
   !> - Unchecked shapes: points (3,npts), xyz (3,nat), numbers (nat), w (npts)
   !> - Valid owner, atomic numbers and positive optional widths required
   !> - Default: size-adjusted elliptic coordinates, clipped to [-1,1]
   !> - Power width selects (q_i-q_j)/width, q_i=|r-R_i|**2-radius_i**2
   !> - Power coordinates have no distance norms or geometry-dependent minima
   !> - Point threading suppressed inside an enclosing OpenMP team
   !>
   !> @param[in]  owner       Owner index in [1,nat]
   !> @param[in]  points      Sample points in bohr
   !> @param[in]  xyz         Atomic positions in bohr
   !> @param[in]  numbers     Atomic numbers
   !> @param[out] w           Normalized owner weights
   !> @param[in]  cell        Pure pair-cell function
   !> @param[in]  power_width Optional power-coordinate half-width in bohr**2
   !> @param[in]  radii       Optional power radii in bohr, shape (nat); default covalent radii
   subroutine pair_partition_weights(owner, points, xyz, numbers, w, cell, power_width, radii)
      !> Owner index
      integer, intent(in) :: owner
      !> Sample points
      real(wp), intent(in) :: points(:, :)
      !> Atomic positions
      real(wp), intent(in) :: xyz(:, :)
      !> Atomic numbers
      integer, intent(in) :: numbers(:)
      !> Owner weights
      real(wp), intent(out) :: w(:)
      !> Pair-cell evaluator
      procedure(partition_cell) :: cell
      !> Power-coordinate half-width
      real(wp), intent(in), optional :: power_width
      !> Power radii
      real(wp), intent(in), optional :: radii(:)

      real(wp), allocatable :: pair_r(:, :), pair_a(:, :), radius(:), dist(:), logs(:)
      real(wp) :: chi, u, width, peak
      integer :: nat, npts, i, j, ip, nthreads
      logical :: power

      nat = size(xyz, 2)
      npts = size(points, 2)
      if (npts == 0) return
      if (nat == 1) then
         w = 1.0_wp
         return
      end if
      power = present(power_width)
      width = 1.0_wp
      if (power) width = power_width
      radius = covalent_rad(numbers)
      if (present(radii)) radius = radii
      allocate (pair_r(nat, nat), pair_a(nat, nat))
      pair_r = 0.0_wp
      pair_a = 0.0_wp
      do i = 1, nat
         do j = 1, i - 1
            if (power) then
               pair_a(j, i) = (radius(j) - radius(i))*(radius(j) + radius(i))
            else
               pair_r(j, i) = norm2(xyz(:, i) - xyz(:, j))
               chi = radius(i)/radius(j)
               u = (chi - 1.0_wp)/(chi + 1.0_wp)
               pair_a(j, i) = max(-0.5_wp, min(0.5_wp, u/(u*u - 1.0_wp)))
            end if
            pair_r(i, j) = pair_r(j, i)
            pair_a(i, j) = -pair_a(j, i)
         end do
      end do

      nthreads = 1
!$    if (.not. omp_in_parallel()) nthreads = min(omp_get_max_threads(), npts)
      !$omp parallel num_threads(nthreads) default(none) &
      !$omp shared(owner, points, xyz, pair_r, pair_a, width, power, w, nat, npts) &
      !$omp private(ip, i, dist, logs, peak)
      allocate (dist(nat), logs(nat))
      !$omp do schedule(dynamic, 64)
      do ip = 1, npts
         dist = 0.0_wp
         if (.not. power) then
            do i = 1, nat
               dist(i) = norm2(points(:, ip) - xyz(:, i))
            end do
         end if
         logs(owner) = cell_log(owner, points(:, ip), xyz, dist, pair_r, pair_a, cell, width, power)
         if (logs(owner) == -huge(1.0_wp)) then
            w(ip) = 0.0_wp
            cycle
         end if
         do i = 1, nat
            if (i /= owner) then
               logs(i) = cell_log(i, points(:, ip), xyz, dist, pair_r, pair_a, cell, width, power)
            end if
         end do
         peak = maxval(logs)
         w(ip) = exp(logs(owner) - peak)/sum(exp(logs - peak))
      end do
      !$omp end do
      deallocate (dist, logs)
      !$omp end parallel

   end subroutine pair_partition_weights

   !> Logarithm of one cell product
   !>
   !> @param[in] ia      Atom index
   !> @param[in] point   Sample point in bohr, shape (3)
   !> @param[in] xyz     Atomic positions in bohr, shape (3,nat)
   !> @param[in] dist    Point-to-atom distances in bohr, shape (nat)
   !> @param[in] pair_r  Pair distances in bohr, shape (nat,nat)
   !> @param[in] pair_a  Pair size adjustments or power offsets, shape (nat,nat)
   !> @param[in] cell    Pair-cell evaluator
   !> @param[in] width   Power-coordinate half-width in bohr**2
   !> @param[in] power   Use power coordinates
   pure function cell_log(ia, point, xyz, dist, pair_r, pair_a, cell, width, power) result(value)
      !> Atom index
      integer, intent(in) :: ia
      !> Sample point
      real(wp), intent(in) :: point(:)
      !> Atomic positions
      real(wp), intent(in) :: xyz(:, :)
      !> Point-to-atom distances
      real(wp), intent(in) :: dist(:)
      !> Pair distances
      real(wp), intent(in) :: pair_r(:, :)
      !> Pair size adjustments or power offsets
      real(wp), intent(in) :: pair_a(:, :)
      !> Pair-cell evaluator
      procedure(partition_cell) :: cell
      !> Power-coordinate half-width
      real(wp), intent(in) :: width
      !> Power-coordinate selection
      logical, intent(in) :: power
      !> Log cell product or -huge for zero
      real(wp) :: value

      integer :: ja
      real(wp) :: mu, nu, s

      value = 0.0_wp
      do ja = 1, size(dist)
         if (ja == ia) cycle
         if (power) then
            ! Difference of powers without cancellation between squared distances
            nu = (dot_product((point - xyz(:, ia)) + (point - xyz(:, ja)), &
               & xyz(:, ja) - xyz(:, ia)) + pair_a(ja, ia))/width
         else
            if (pair_r(ja, ia) <= 0.0_wp) cycle
            mu = max(-1.0_wp, min(1.0_wp, (dist(ia) - dist(ja))/pair_r(ja, ia)))
            nu = mu + pair_a(ja, ia)*(1.0_wp - mu*mu)
         end if
         s = cell(nu)
         if (s <= 0.0_wp) then
            value = -huge(1.0_wp)
            return
         end if
         value = value + log(s)
      end do
   end function cell_log

end module moist_math_grid_3d_partition_common
