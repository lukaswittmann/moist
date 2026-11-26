!> Lattice seeder, seeds of one anchor: ball query, filters, pre-step
submodule(moist_cavity_drop_projection_seed) moist_cavity_drop_projection_seed_query
   implicit none(type, external)

contains

   !* ================================================================================= *!
   !*                                      Seeding                                      *!
   !* ================================================================================= *!

   !> Seeds of one anchor from the built lattice, no level set evaluation
   !>
   !> - the ball is the request's `pdist_start`; when `pdist_inner >= 0`, only the
   !>   outer shell `(pdist_inner, pdist_start]` is emitted and
   !>   `scratch%shell_filtered` is set
   !> - the filters switched on drop samples ([[step_candidate]]), the pre-step
   !>   moves the kept ones; the back-facing filter needs the anchor's sphere
   !>   and is skipped for a request without an owner (`owner == 0`)
   !> - `reach`: the covering radius plus the pre-step's longest move
   !> - `nearest_first` orders the seeds by distance to the anchor, ties in
   !>   sample order
   !>
   !> @param[in]    self    Seeder (built)
   !> @param[inout] scratch Calling thread's scratch
   !> @param[in]    lsf     Level set; its stationarity certificate only
   !> @param[in]    request Anchor request, radius positive and finite
   !> @param[out]   error   Invalid radius, seeder used before its build
   module subroutine seeder_seed(self, scratch, lsf, request, error)
      !> Seeder (built)
      class(projection_seeder_type), intent(in) :: self
      !> Calling thread's scratch
      type(seeder_scratch_type), intent(inout) :: scratch
      !> Level set; its stationarity certificate only
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Anchor request, radius positive and finite
      type(seed_request_type), intent(in) :: request
      !> Invalid radius, seeder used before its build
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: dx(3), x0(3), n_sph(3), nsph_len
      integer :: i, n, nkeep
      logical :: keep, do_hor

      ! The lattice answers without evaluation; `lsf` only certifies bins
      scratch%n_points = 0
      scratch%reach = 0.0_wp
      scratch%shell_filtered = .false.
      ! `points` is allocated even when no seed follows
      if (.not. allocated(scratch%points)) allocate (scratch%points(3, 32))
      if (.not. (request%pdist_start > 0.0_wp)) then
         call fatal_error(error, "Lattice seeder: request radius must be positive")
         return
      end if
      if (.not. ieee_is_finite(request%pdist_start)) then
         call fatal_error(error, "Lattice seeder: request radius must be finite")
         return
      end if
      if (.not. self%built) then
         call fatal_error(error, "Lattice seeder: seed before build")
         return
      end if
      scratch%shell_filtered = request%pdist_inner >= 0.0_wp
      if (self%nsample == 0) return
      scratch%reach = self%cover + self%spacing_max

      ! Samples of the cells meeting the ball (or outer shell)
      call ball_query(self, scratch, request%anchor, request%pdist_start, request%pdist_inner, lsf)
      n = scratch%ncand
      if (n == 0) return
      n_sph = request%anchor - request%owner_center
      nsph_len = norm2(n_sph)
      ! Without an owner `owner_center` is not the anchor's sphere center
      do_hor = self%filter_back_facing .and. request%has_accepted .and. request%owner > 0 &
               .and. (nsph_len > 0.0_wp)
      if (do_hor) n_sph = n_sph/nsph_len
      nkeep = 0
      do i = 1, n
         call step_candidate(self, scratch%cand(i), request%anchor, request%has_accepted, do_hor, n_sph, keep, x0)
         if (.not. keep) cycle
         nkeep = nkeep + 1
         scratch%cand_x0(:, nkeep) = x0
         if (request%nearest_first) then
            dx = x0 - request%anchor
            scratch%cand_dist(nkeep) = dot_product(dx, dx)
         end if
      end do
      if (nkeep == 0) return
      call reserve_points(scratch, nkeep)
      if (request%nearest_first) then
         call sort_index(scratch%cand_dist(1:nkeep), scratch%cand_perm(1:nkeep))
         do i = 1, nkeep
            scratch%points(:, i) = scratch%cand_x0(:, scratch%cand_perm(i))
         end do
      else
         scratch%points(:, 1:nkeep) = scratch%cand_x0(:, 1:nkeep)
      end if
      scratch%n_points = nkeep
   end subroutine seeder_seed

   !> Samples whose 4-cell box meets the ball `radius` (and lies strictly outside
   !> `inner` when non-negative), appended to the candidates
   !>
   !> - samples are binned by `edge_mid`, so every point of a bin's 4-cell boxes
   !>   lies within `0.5 bin_size + spacing_max` of the bin center (independent of
   !>   the retraction distance `dmax`)
   !> - bins and boxes in the lattice frame; the stationarity certificate sees
   !>   lab boxes, the axis-aligned hulls of the turned ones (larger, so it
   !>   excludes less, never more)
   !> - a bin is skipped when its inflated box misses the outer ball, lies
   !>   strictly inside `inner`, or is certified stationary-free by `lsf`
   !>
   !> @param[in]    self    Seeder
   !> @param[inout] scratch Scratch; kept candidates written
   !> @param[in]    x_lab   Point (lab)
   !> @param[in]    radius  Outer ball radius (Bohr)
   !> @param[in]    inner   Inner excluded ball radius (Bohr), negative for full ball
   !> @param[in]    lsf     Level set whose stationarity certificate prunes bins
   subroutine ball_query(self, scratch, x_lab, radius, inner, lsf)
      !> Seeder
      type(projection_seeder_type), intent(in) :: self
      !> Scratch; kept candidates written
      type(seeder_scratch_type), intent(inout) :: scratch
      !> Point (lab)
      real(wp), intent(in) :: x_lab(3)
      !> Outer ball radius (Bohr)
      real(wp), intent(in) :: radius
      !> Inner excluded ball radius (Bohr), negative for full ball
      real(wp), intent(in) :: inner
      !> Level set whose stationarity certificate prunes bins
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf

      integer :: lo(3), hi(3), ix, iy, iz, b, k
      real(wp) :: limit, half_box, half_box_pad, half_bin_pad, r2, r2_hi, r_inner2, r_inner2_lo
      real(wp) :: center(3), reg_lo(3), reg_hi(3), atom_lo(3), atom_hi(3)
      !> Point in the lattice frame, lab hull widths of a unit box, a box's lab center and half width
      real(wp) :: x(3), stretch(3), center_lab(3), half_lab
      real(wp) :: adz, ady, adx, dz_min, dy_min, dx_min, dz_min2, dy_min2, d_min2
      real(wp) :: dz_max, dy_max, dx_max, dz_max2, dy_max2, d_max2, d2
      logical :: can_exclude, use_inner, cull_inner

      scratch%ncand = 0
      if (self%nsample == 0) return
      x = self%frame%to_lattice(x_lab)
      limit = nearest(radius + self%spacing_max, 1.0_wp)
      lo = max(bin_coords(self, x - limit), -key_bias + 1)
      hi = min(bin_coords(self, x + limit), key_bias - 1)
      half_box = nearest(0.5_wp*self%bin_size + self%spacing_max, 1.0_wp)
      half_box_pad = half_box*(1.0_wp + 1.0e-12_wp)
      half_bin_pad = 0.5_wp*self%bin_size*(1.0_wp + 1.0e-12_wp)
      r2 = radius*radius
      r2_hi = r2*(1.0_wp + 1.0e-12_wp) + 1.0e-24_wp
      use_inner = inner >= 0.0_wp
      cull_inner = inner > 0.0_wp
      r_inner2 = 0.0_wp
      r_inner2_lo = 0.0_wp
      if (use_inner) then
         r_inner2 = inner*inner
         r_inner2_lo = r_inner2*(1.0_wp - 1.0e-12_wp)
      end if
      reg_lo = nearest(nearest(self%bin_origin + self%bin_size*real(lo, wp), -1.0_wp) - self%spacing_max, -1.0_wp)
      reg_hi = nearest(nearest(self%bin_origin + self%bin_size*real(hi + 1, wp), 1.0_wp) + self%spacing_max, 1.0_wp)
      half_lab = half_box
      if (self%frame%rotated) then
         ! Lab hull of the region, widened by that of a tested box beyond its own
         stretch = sum(abs(self%frame%axes), dim=2)
         half_lab = maxval(stretch)*half_box
         center_lab = self%frame%to_lab(0.5_wp*(reg_lo + reg_hi))
         reg_hi = matmul(abs(self%frame%axes), 0.5_wp*(reg_hi - reg_lo)) + half_lab
         reg_lo = nearest(center_lab - reg_hi, -1.0_wp)
         reg_hi = nearest(center_lab + reg_hi, 1.0_wp)
      end if
      call lsf%stationary_region_bounds(x_lab, reg_lo, reg_hi, atom_lo, atom_hi, can_exclude)
      do iz = lo(3), hi(3)
         center(3) = self%bin_origin(3) + self%bin_size*(real(iz, wp) + 0.5_wp)
         adz = abs(x(3) - center(3))
         dz_min = max(adz - half_box_pad, 0.0_wp)
         dz_min2 = dz_min*dz_min
         if (dz_min2 > r2_hi) cycle
         dz_max = adz + half_bin_pad
         dz_max2 = dz_max*dz_max
         do iy = lo(2), hi(2)
            center(2) = self%bin_origin(2) + self%bin_size*(real(iy, wp) + 0.5_wp)
            ady = abs(x(2) - center(2))
            dy_min = max(ady - half_box_pad, 0.0_wp)
            dy_min2 = dz_min2 + dy_min*dy_min
            if (dy_min2 > r2_hi) cycle
            dy_max = ady + half_bin_pad
            dy_max2 = dz_max2 + dy_max*dy_max
            do ix = lo(1), hi(1)
               center(1) = self%bin_origin(1) + self%bin_size*(real(ix, wp) + 0.5_wp)
               adx = abs(x(1) - center(1))
               dx_min = max(adx - half_box_pad, 0.0_wp)
               d_min2 = dy_min2 + dx_min*dx_min
               if (d_min2 > r2_hi) cycle
               if (cull_inner) then
                  dx_max = adx + half_bin_pad
                  d_max2 = dy_max2 + dx_max*dx_max
                  if (d_max2 < r_inner2_lo) cycle
               end if
               b = table_find(self%bins, lattice_key([ix, iy, iz]))
               if (b == 0) cycle
               if (can_exclude) then
                  if (lsf%stationary_box_excluded(x_lab, self%frame%to_lab(center), half_lab, atom_lower=atom_lo, &
                                                  atom_upper=atom_hi)) cycle
               end if
               do k = self%bin_first(b), self%bin_first(b + 1) - 1
                  d2 = edge_box_distance_sq(self, k, x)
                  if (d2 > r2) cycle
                  if (use_inner .and. d2 <= r_inner2) cycle
                  call push_candidate(scratch, k)
               end do
            end do
         end do
      end do
   end subroutine ball_query

   !> Append one candidate sample
   !>
   !> @param[inout] self Scratch
   !> @param[in]    k    Sample
   subroutine push_candidate(self, k)
      !> Scratch
      type(seeder_scratch_type), intent(inout) :: self
      !> Sample
      integer, intent(in) :: k

      integer :: cap

      if (.not. allocated(self%cand)) then
         allocate (self%cand(64), self%cand_perm(64), self%cand_x0(3, 64), self%cand_dist(64))
      end if
      if (self%ncand + 1 > size(self%cand)) then
         ! Only the candidates matter; their starts are written per request
         cap = 2*size(self%cand)
         call grow_array(self%cand, cap)
         call grow_array(self%cand_perm, cap)
         call grow_array(self%cand_x0, 3, cap)
         call grow_array(self%cand_dist, cap)
      end if
      self%ncand = self%ncand + 1
      self%cand(self%ncand) = k
   end subroutine push_candidate

   !> Ensure `scratch%points` holds at least `n` columns; contents are not kept
   !>
   !> @param[inout] scratch Scratch, `points` allocated
   !> @param[in]    n       Required columns
   subroutine reserve_points(scratch, n)
      !> Scratch, `points` allocated
      type(seeder_scratch_type), intent(inout) :: scratch
      !> Required columns
      integer, intent(in) :: n

      if (size(scratch%points, 2) < n) call grow_array(scratch%points, 3, max(n, 2*size(scratch%points, 2)))
   end subroutine reserve_points

   !> Squared distance from a point to the four cells around a sample's edge
   !>
   !> @param[in] self Seeder
   !> @param[in] v    Sample
   !> @param[in] x    Point
   pure function edge_box_distance_sq(self, v, x) result(dist2)
      !> Seeder
      type(projection_seeder_type), intent(in) :: self
      !> Sample
      integer, intent(in) :: v
      !> Point
      real(wp), intent(in) :: x(3)
      !> Squared distance (Bohr^2)
      real(wp) :: dist2

      real(wp) :: half(3), d(3)

      half = self%spacing(v)
      half(self%edge_axis(v)) = 0.5_wp*self%spacing(v)
      d = max(abs(x - self%edge_mid(:, v)) - half, 0.0_wp)
      dist2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
   end function edge_box_distance_sq

   !> Filter one sample for an anchor and move it by the osculating pre-step
   !>
   !> - a sample without a normal is kept and not moved
   !> - each filter drops a sample on its own, when switched on:
   !>   - `filter_back_facing` (`do_hor`, an anchor with an accepted local point):
   !>     `n . n_sph <= back_facing_cos`; the horizon switch closes a point there
   !>   - `filter_grazing`: `r > h` and `|d . n| / r < 0.20`, a patch seen edge-on
   !>   - `filter_cone` (an anchor with an accepted local point):
   !>     `||g|| > cone_scale h (1 + r k_max)`, `g` the tangential part of `d`;
   !>     at a point `g = 0`, and it grows by at most `1 + r k_max` per distance
   !>     along the surface, so a sample outside the cone is no point's nearest one
   !>   - `filter_focal`: `I - sigma II` is negative definite, the sample lies
   !>     beyond both focal centers; a point there is a maximum of the distance,
   !>     which the focal switch closes. Beyond one center only the sample is
   !>     kept: Newton leaves it along the crease for an open point
   !> - a kept sample with `I - sigma II` positive definite steps towards the point
   !>   on its osculating quadric, the tangent step `t` solving `(I - sigma II) t = -g`,
   !>   capped at `h`, and the normal offset `z = -t^T II t / 2`; `t` and `z` are then
   !>   scaled together (`z` with the square) so the move `|t|^2 + z^2` stays within
   !>   `h^2`, the bound the seeder's `reach` counts
   !>
   !> @param[in]  self         Seeder
   !> @param[in]  v            Sample index
   !> @param[in]  anchor       Anchor position (3)
   !> @param[in]  has_accepted Whether the anchor already has an accepted local point
   !> @param[in]  do_hor       Whether to drop back-facing samples
   !> @param[in]  n_sph        Unit outward anchor-sphere direction (`(anchor - owner_center)/r`)
   !> @param[out] keep         Whether the sample survives
   !> @param[out] x0           Starting point (valid when `keep` is `.true.`)
   pure subroutine step_candidate(self, v, anchor, has_accepted, do_hor, n_sph, keep, x0)
      !> Seeder
      type(projection_seeder_type), intent(in) :: self
      !> Sample index
      integer, intent(in) :: v
      !> Anchor position (3)
      real(wp), intent(in) :: anchor(3)
      !> Whether the anchor already has an accepted local point
      logical, intent(in) :: has_accepted
      !> Whether to drop back-facing samples
      logical, intent(in) :: do_hor
      !> Unit outward anchor-sphere direction (`(anchor - owner_center)/r`)
      real(wp), intent(in) :: n_sph(3)
      !> Whether the sample survives
      logical, intent(out) :: keep
      !> Starting point (valid when `keep` is `.true.`)
      real(wp), intent(out) :: x0(3)

      real(wp) :: q1(3), q2(3), d(3), nrm(3), sigma, g1, g2, m11, m12, m22, det, t(2), tn(2), z
      real(wp) :: tlen2, h, r2, g2_tan, thr, scale2

      keep = .true.
      x0 = self%y(:, v)
      nrm = self%normal(:, v)
      if (.not. (dot_product(nrm, nrm) > 0.0_wp)) return
      if (do_hor .and. dot_product(nrm, n_sph) <= back_facing_cos) then
         keep = .false.
         return
      end if
      d = x0 - anchor
      sigma = dot_product(d, nrm)
      h = self%spacing(v)

      r2 = dot_product(d, d)
      ! Grazing: r > h and |d . n| / r < 0.20
      if (self%filter_grazing .and. r2 > h*h .and. sigma*sigma < grazing_cos2*r2) then
         keep = .false.
         return
      end if
      ! Outside the tangent cone, for an anchor with an accepted local point
      if (self%filter_cone .and. has_accepted) then
         g2_tan = max(0.0_wp, r2 - sigma*sigma)
         thr = self%cone_scale*h*(1.0_wp + sqrt(r2)*self%kmax(v))
         if (g2_tan > thr*thr) then
            keep = .false.
            return
         end if
      end if

      m11 = 1.0_wp - sigma*self%form(1, v)
      m12 = -sigma*self%form(2, v)
      m22 = 1.0_wp - sigma*self%form(3, v)
      det = m11*m22 - m12*m12

      ! Beyond both focal centers: I - sigma II negative definite <=> det > 0 and m11 < 0
      if (self%filter_focal .and. det > 0.0_wp .and. m11 < 0.0_wp) then
         keep = .false.
         return
      end if

      if (.not. (m11 > 0.0_wp .and. det > 0.0_wp)) return
      call setup_tangent_frame(nrm, q1, q2)
      g1 = dot_product(d, q1)
      g2 = dot_product(d, q2)
      t(1) = -(m22*g1 - m12*g2)/det
      t(2) = -(-m12*g1 + m11*g2)/det
      tlen2 = t(1)*t(1) + t(2)*t(2)
      if (.not. ieee_is_finite(tlen2)) then
         ! Overflow at a tiny `det`: the direction from the numerators, length `h`
         tn = [m22*g1 - m12*g2, -m12*g1 + m11*g2]
         t = -tn*(h/norm2(tn))
         tlen2 = t(1)*t(1) + t(2)*t(2)
      end if
      ! The tangent step alone first, so `z` stays finite near a focal point
      if (tlen2 > h*h) then
         t = t*(h/sqrt(tlen2))
         tlen2 = t(1)*t(1) + t(2)*t(2)
      end if
      z = -0.5_wp*(self%form(1, v)*t(1)*t(1) + 2.0_wp*self%form(2, v)*t(1)*t(2) + self%form(3, v)*t(2)*t(2))
      if (tlen2 + z*z > h*h) then
         ! Scale s: s^2 |t|^2 + s^4 z^2 = h^2, the root in s^2 written without cancellation
         scale2 = 2.0_wp*h*h/(tlen2 + sqrt(tlen2*tlen2 + 4.0_wp*z*z*h*h))
         t = t*sqrt(scale2)
         z = z*scale2
      end if
      x0 = x0 + t(1)*q1 + t(2)*q2 + z*nrm
   end subroutine step_candidate

end submodule moist_cavity_drop_projection_seed_query
