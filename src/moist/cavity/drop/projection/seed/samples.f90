!> Lattice seeder build, samples: retraction onto the surface, bins
submodule(moist_cavity_drop_projection_seed) moist_cavity_drop_projection_seed_samples
   implicit none(type, external)

contains

   !* ================================================================================= *!
   !*                                      Samples                                      *!
   !* ================================================================================= *!

   !> Retract the crossings onto `S = 0` and store normal and second fundamental form
   !>
   !> - `retract_steps` normal Newton steps of at most one edge length each, then
   !>   the order-2 jet and one more step without evaluation when it is shorter
   !>   than the edge
   !> - a sample whose level set is unreached or flat keeps its position and gets
   !>   no normal
   !> - the covering radius is a cell diagonal plus the longest retraction
   !>
   !> @param[inout] self  Seeder, samples written
   !> @param[in]    lsf   Level set, the source of the clones
   !> @param[inout] work  Build scratch
   !> @param[out]   error Level set failure
   module subroutine retract_samples(self, lsf, work, error)
      type(projection_seeder_type), intent(inout) :: self
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      type(build_type), intent(inout) :: work
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: y(:, :), val(:), grad(:, :), hess(:, :)
      logical, allocatable :: reached(:), ok(:)
      integer :: n, k, it
      real(wp) :: gn, step(3), slen, q1(3), q2(3), h3(3, 3), dmax

      n = work%ncross
      self%nsample = n
      allocate (self%y(3, n), self%normal(3, n), self%form(3, n), self%kmax(n), self%spacing(n))
      self%edge_mid = work%cross_mid(:, 1:n)
      self%edge_axis = work%cross_axis(1:n)
      self%spacing_max = 0.0_wp
      if (n > 0) self%spacing_max = maxval(work%cross_h(1:n))
      self%cover = 0.0_wp
      if (n == 0) return

      allocate (y(3, n), val(n), grad(3, n), hess(6, n), reached(n), ok(n))
      y = work%cross(:, 1:n)
      ok = .true.
      do it = 1, retract_steps
         call evaluate_slices(work%probe, lsf, y, 1, build_chunk, val, reached, error, grad=grad)
         if (allocated(error)) return
         self%n_order1 = self%n_order1 + n
         do k = 1, n
            if (.not. ok(k)) cycle
            gn = norm2(grad(:, k))
            if (.not. reached(k) .or. .not. (gn > 0.0_wp)) then
               ok(k) = .false.
               cycle
            end if
            step = -(val(k)/(gn*gn))*grad(:, k)
            slen = norm2(step)
            if (slen > work%cross_h(k)) step = step*(work%cross_h(k)/slen)
            y(:, k) = y(:, k) + step
         end do
      end do
      call evaluate_slices(work%probe, lsf, y, 2, build_chunk, val, reached, error, grad=grad, hess=hess)
      if (allocated(error)) return
      self%n_order2 = self%n_order2 + n

      dmax = 0.0_wp
      do k = 1, n
         self%normal(:, k) = 0.0_wp
         self%form(:, k) = 0.0_wp
         self%kmax(k) = 0.0_wp
         self%spacing(k) = work%cross_h(k)
         gn = norm2(grad(:, k))
         if (ok(k) .and. reached(k) .and. gn > 0.0_wp) then
            self%normal(:, k) = grad(:, k)/gn
            call setup_tangent_frame(self%normal(:, k), q1, q2)
            call unpack_sym_matrix(hess(:, k), h3, 3)
            self%form(1, k) = dot_product(q1, matmul(h3, q1))/gn
            self%form(2, k) = dot_product(q1, matmul(h3, q2))/gn
            self%form(3, k) = dot_product(q2, matmul(h3, q2))/gn
            self%kmax(k) = 0.5_wp*abs(self%form(1, k) + self%form(3, k)) &
                           + hypot(0.5_wp*(self%form(1, k) - self%form(3, k)), self%form(2, k))
            step = -(val(k)/(gn*gn))*grad(:, k)
            if (norm2(step) <= work%cross_h(k)) y(:, k) = y(:, k) + step
         end if
         self%y(:, k) = y(:, k)
         dmax = max(dmax, norm2(y(:, k) - work%cross(:, k)))
      end do
      self%cover = sqrt(3.0_wp)*self%spacing_max + dmax
   end subroutine retract_samples

   !> Sort the samples into cubic bins of edge `edge` by their edge midpoint
   !>
   !> @param[inout] self Seeder
   !> @param[in]    edge Bin edge (Bohr)
   module subroutine build_bins(self, edge)
      type(projection_seeder_type), intent(inout) :: self
      real(wp), intent(in) :: edge

      integer, allocatable :: bin_of(:), order(:), count_of(:)
      integer(int64) :: key
      integer :: k, nbin, b, pos

      self%bin_size = edge
      self%bin_origin = 0.0_wp
      if (self%nsample > 0) self%bin_origin = minval(self%edge_mid, dim=2)
      call table_init(self%bins, max(16, self%nsample/4))
      allocate (bin_of(self%nsample))
      nbin = 0
      do k = 1, self%nsample
         key = lattice_key(bin_coords(self, self%edge_mid(:, k)), key_bias)
         b = table_find(self%bins, key)
         if (b == 0) then
            nbin = nbin + 1
            b = nbin
            call table_insert(self%bins, key, b)
         end if
         bin_of(k) = b
      end do

      ! Counting sort by bin, stable in sample order
      allocate (count_of(nbin + 1), source=0)
      do k = 1, self%nsample
         count_of(bin_of(k) + 1) = count_of(bin_of(k) + 1) + 1
      end do
      allocate (self%bin_first(nbin + 1))
      self%bin_first(1) = 1
      do b = 1, nbin
         self%bin_first(b + 1) = self%bin_first(b) + count_of(b + 1)
      end do
      allocate (order(self%nsample))
      count_of(1:nbin) = self%bin_first(1:nbin)
      do k = 1, self%nsample
         pos = count_of(bin_of(k))
         order(pos) = k
         count_of(bin_of(k)) = pos + 1
      end do
      self%y = self%y(:, order)
      self%normal = self%normal(:, order)
      self%form = self%form(:, order)
      self%kmax = self%kmax(order)
      self%spacing = self%spacing(order)
      self%edge_mid = self%edge_mid(:, order)
      self%edge_axis = self%edge_axis(order)
   end subroutine build_bins

   !> Bin coordinates of a point
   !>
   !> @param[in] self Seeder
   !> @param[in] x    Point
   pure module function bin_coords(self, x) result(j)
      type(projection_seeder_type), intent(in) :: self
      real(wp), intent(in) :: x(3)
      integer :: j(3)

      j = floor((x - self%bin_origin)/self%bin_size)
   end function bin_coords

end submodule moist_cavity_drop_projection_seed_samples
