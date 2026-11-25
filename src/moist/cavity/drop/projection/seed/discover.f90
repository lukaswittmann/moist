!> Surface mesh build, discovery: coarse nodes of the search region
submodule(moist_cavity_drop_projection_seed) moist_cavity_drop_projection_seed_discover
   implicit none(type, external)

contains

   !* ================================================================================= *!
   !*                                     Discovery                                     *!
   !* ================================================================================= *!

   !> Probe every coarse node within a cell diagonal of the searched anchor balls
   !>
   !> - every coarse cell meeting any anchor's search ball `B(a_k, r_k)` has all
   !>   eight corners within `r_k + sqrt(3) h_0` of `a_k`
   !> - stamps the union of balls into a compact 3D byte mask when the bounding
   !>   box fits [[max_mask_cells]] and the lattice key range, emitting nodes once
   !>   in deterministic `(z, y, x)` order; falls back to direct insertion otherwise
   !>
   !> @param[inout] self     Seeder (tallies)
   !> @param[in]    lsf      Level set (certificates)
   !> @param[in]    geometry Build input
   !> @param[inout] work     Build scratch
   !> @param[out]   error    Level set failure or lattice out of range
   module subroutine discover_volume(self, lsf, geometry, work, error)
      type(projection_seeder_type), intent(inout) :: self
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      type(seed_build_type), intent(in) :: geometry
      type(build_type), intent(inout) :: work
      type(error_type), allocatable, intent(out) :: error

      integer(int8), allocatable :: mask(:, :, :)
      integer(int64) :: glo(3), ghi(3), lo64(3), hi64(3), nx, ny, nz, step64, key_lim
      integer(int64) :: ix64, iy64, iz64
      integer :: k, nanc, step, idx, ix, iy, iz, ix_lo, ix_hi, lo(3), hi(3), ox, oy, oz, nmark, cap, ijk(3)
      real(wp) :: pad, r_k, rho, rho2, c(3), dz, dz2, dy, dy2, rem, span_x, inv_h0, x(3)
      logical :: has_rmax, use_mask

      nanc = size(geometry%anchors, 2)
      if (nanc == 0) return
      has_rmax = allocated(geometry%radius_max)
      step = level_step(work, 0)
      step64 = int(step, int64)
      key_lim = int(key_bias/2, int64)
      pad = sqrt(3.0_wp)*work%h0*(1.0_wp + 1.0e-12_wp)
      inv_h0 = 1.0_wp/work%h0

      glo = huge(1_int64)
      ghi = -huge(1_int64)
      do k = 1, nanc
         if (has_rmax) then
            r_k = geometry%radius_max(k)
         else
            r_k = geometry%reach
         end if
         rho = r_k + pad
         lo64 = floor((geometry%anchors(:, k) - rho - work%origin)*inv_h0, kind=int64)
         hi64 = ceiling((geometry%anchors(:, k) + rho - work%origin)*inv_h0, kind=int64)
         glo = min(glo, lo64)
         ghi = max(ghi, hi64)
      end do

      nx = ghi(1) - glo(1) + 1_int64
      ny = ghi(2) - glo(2) + 1_int64
      nz = ghi(3) - glo(3) + 1_int64
      use_mask = nx > 0_int64 .and. ny > 0_int64 .and. nz > 0_int64 &
                 .and. nx <= 4096_int64 .and. ny <= 4096_int64 .and. nz <= 4096_int64 &
                 .and. all(abs(step64*glo) < key_lim) .and. all(abs(step64*ghi) < key_lim)
      if (use_mask) use_mask = nx*ny*nz <= max_mask_cells

      if (use_mask) then
         allocate (mask(int(nx), int(ny), int(nz)), source=0_int8)
         ox = int(glo(1)) - 1
         oy = int(glo(2)) - 1
         oz = int(glo(3)) - 1
         do k = 1, nanc
            if (has_rmax) then
               r_k = geometry%radius_max(k)
            else
               r_k = geometry%reach
            end if
            rho = r_k + pad
            rho2 = rho*rho
            c = geometry%anchors(:, k)
            lo = floor((c - rho - work%origin)*inv_h0)
            hi = ceiling((c + rho - work%origin)*inv_h0)
            do iz = lo(3), hi(3)
               dz = work%origin(3) + work%h0*real(iz, wp) - c(3)
               dz2 = dz*dz
               if (dz2 > rho2) cycle
               do iy = lo(2), hi(2)
                  dy = work%origin(2) + work%h0*real(iy, wp) - c(2)
                  dy2 = dy*dy
                  rem = rho2 - dz2 - dy2
                  if (rem < 0.0_wp) cycle
                  span_x = sqrt(rem)
                  ix_lo = max(lo(1), floor((c(1) - span_x - work%origin(1))*inv_h0))
                  ix_hi = min(hi(1), ceiling((c(1) + span_x - work%origin(1))*inv_h0))
                  if (ix_lo <= ix_hi) mask(ix_lo - ox:ix_hi - ox, iy - oy, iz - oz) = 1_int8
               end do
            end do
         end do

         nmark = 0
         do iz = 1, int(nz)
            do iy = 1, int(ny)
               do ix = 1, int(nx)
                  if (mask(ix, iy, iz) /= 0_int8) nmark = nmark + 1
               end do
            end do
         end do
         cap = max(1024, nmark)
         call table_init(work%nodes%table, max(4096, nmark))
         allocate (work%nodes%ijk(3, cap), work%nodes%val(cap), work%nodes%free(cap), work%nodes%reached(cap))
         do iz = 1, int(nz)
            do iy = 1, int(ny)
               do ix = 1, int(nx)
                  if (mask(ix, iy, iz) == 0_int8) cycle
                  work%nodes%n = work%nodes%n + 1
                  idx = work%nodes%n
                  ijk = step*[ix + ox, iy + oy, iz + oz]
                  work%nodes%ijk(:, idx) = ijk
                  call table_insert(work%nodes%table, lattice_key(ijk, key_bias), idx)
               end do
            end do
         end do
      else
         call table_init(work%nodes%table, 4096)
         do k = 1, nanc
            if (has_rmax) then
               r_k = geometry%radius_max(k)
            else
               r_k = geometry%reach
            end if
            rho = r_k + pad
            rho2 = rho*rho
            c = geometry%anchors(:, k)
            lo64 = floor((c - rho - work%origin)*inv_h0, kind=int64)
            hi64 = ceiling((c + rho - work%origin)*inv_h0, kind=int64)
            do iz64 = lo64(3), hi64(3)
               x(3) = work%origin(3) + work%h0*real(iz64, wp)
               dz2 = (x(3) - c(3))**2
               if (dz2 > rho2) cycle
               do iy64 = lo64(2), hi64(2)
                  x(2) = work%origin(2) + work%h0*real(iy64, wp)
                  dy2 = (x(2) - c(2))**2
                  if (dz2 + dy2 > rho2) cycle
                  do ix64 = lo64(1), hi64(1)
                     x(1) = work%origin(1) + work%h0*real(ix64, wp)
                     if (dz2 + dy2 + (x(1) - c(1))**2 > rho2) cycle
                     ! Every later node lies between coarse ones, so this keeps all of them in range
                     if (any(abs(step64*[ix64, iy64, iz64]) >= key_lim)) then
                        call range_error(error)
                        return
                     end if
                     call add_node(work%nodes, step*[int(ix64), int(iy64), int(iz64)])
                  end do
               end do
            end do
         end do
      end if
      call evaluate_nodes(self, lsf, work, error)
   end subroutine discover_volume

end submodule moist_cavity_drop_projection_seed_discover
