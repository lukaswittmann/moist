module moist_cavity_iswig
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use mctc_io_structure, only: structure_type
   use mctc_io, only: new
   use mctc_env, only: error_type, fatal_error, wp
   use iso_fortran_env, only: error_unit, output_unit

   use moist_math_grid_lebedev, only: get_angular_grid, grid_size, lebedev_order_from_num
   use moist_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_context, only: moist_context_type
   use moist_radius_type, only: radius_type

   implicit none
   private
   public :: cavity_type_iswig, new_cavity_iswig

   ! iSwiG implementation of cavity
   type, extends(cavity_type) :: cavity_type_iswig

      !> Number of Lebedev points per sphere
      integer :: num_leb = 110
      !> Default area cutoff
      real(wp) :: cut_a = 0.0_wp
      !> Default iSwiG value cutoff
      real(wp) :: cut_f = 1.0E-10_wp

      !> Raw Lebedev weights (ngrid)
      real(wp), allocatable :: wleb(:)
      !> Accumulated gradient of total area w.r.t. atom positions (3, nsph)
      real(wp), allocatable :: area_grad(:, :)
      !> Accumulated gradient of total volume w.r.t. atom positions (3, nsph)
      real(wp), allocatable :: volume_grad(:, :)
      !> Grid point number (not changed after removing points)
      integer, allocatable :: numbering(:)

      ! Cached Lebedev data (reused across updates)
      integer :: cached_num_leb = 0
      integer :: cached_oleb = 0
      real(wp) :: cached_swx = 0.0_wp
      real(wp), allocatable :: ang_grid(:, :) ! (3, num_leb)
      real(wp), allocatable :: ang_weight(:)  ! (num_leb)

   contains
      procedure :: update => update_cavity_iswig
      procedure :: get_gradient => compute_gradient_iswig
      !> Contract surface-observable adjoints into the nuclear gradient
      procedure :: get_surface_gradient => get_surface_gradient_iswig
      procedure :: write_csv_debug => write_cavity_csv_debug
   end type cavity_type_iswig

contains

   !> Constructor for iSwiG cavity
   !> Initialize an already-declared object; no allocation of the object itself.
   subroutine new_cavity_iswig(self, ctx, nleb, cut_a, cut_f, radius_model, error)
      !> Cavity type instance to initialize
      type(cavity_type_iswig), intent(inout) :: self
      !> Shared run context (verbosity/debug/timer); borrowed, must outlive self
      type(moist_context_type), intent(in), target :: ctx
      !> Number of lebedev grid points per unit sphere
      integer, intent(in), optional :: nleb
      !> Settings for iSwiG cavity
      real(wp), intent(in), optional :: cut_a, cut_f
      !> Enable simplified mode
      !> Optional radii model
      class(radius_type), intent(in) :: radius_model
      !> Constructor error
      type(error_type), allocatable, intent(out) :: error

      !> Borrow the shared run context (owns verbosity/debug/timer)
      self%ctx => ctx

      !> Set configuration values (leave previously allocated buffers untouched)
      if (present(nleb)) self%num_leb = nleb
      if (present(cut_a)) self%cut_a = cut_a
      if (present(cut_f)) self%cut_f = cut_f
      if (allocated(self%radius_model)) deallocate (self%radius_model)
      allocate (self%radius_model, source=radius_model)

   end subroutine new_cavity_iswig

   !> Write grid to CSV, including numbering, Lebedev weight, and switching value
   subroutine write_cavity_csv_debug(self, filename, error)
      class(cavity_type_iswig), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, ipt

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open CSV file for writing: "//trim(filename))
         return
      end if

      write (unit, "(a)") "ngrid,numbering,x,y,z,owner,radius,area,w_leb,f"

      do ipt = 1, self%ngrid
         write (unit, '(i0,10('','',g0))') ipt, self%numbering(ipt), &
            self%xyz(1, ipt), self%xyz(2, ipt), self%xyz(3, ipt), &
            self%owner(ipt), self%radii(self%owner(ipt)), &
            self%a(ipt), self%wleb(ipt), self%f(ipt)
      end do
      close (unit)

      write (output_unit, "(a,1x,a)") "[Info] Wrote cavity grid to", trim(filename)

   end subroutine write_cavity_csv_debug

   !> Update surface with the current new geometry (iSwiG implementation)
   subroutine update_cavity_iswig(self, mol, error)
      class(cavity_type_iswig), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error

      !> Set number of spheres
      self%nsph = mol%nat

      ! Nuclear derivative arrays belong to the previous geometry until the
      ! caller explicitly requests a fresh gradient build.
      if (allocated(self%area_grad)) deallocate (self%area_grad)
      if (allocated(self%volume_grad)) deallocate (self%volume_grad)
      if (allocated(self%xi1_rA)) deallocate (self%xi1_rA)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      if (allocated(self%xyz1_rA)) deallocate (self%xyz1_rA)
      if (allocated(self%v1_rA)) deallocate (self%v1_rA)

      call self%radius_model%update(mol, error)
      if (allocated(error)) return
      if (allocated(self%radii)) deallocate (self%radii)
      self%radii = self%radius_model%f0

      !> Set centers of spheres
      if (allocated(self%sphxyz)) deallocate (self%sphxyz)
      allocate (self%sphxyz(3, self%nsph))
      self%sphxyz = mol%xyz

      !> Allocation
      if (allocated(self%asph)) deallocate (self%asph)
      allocate (self%asph(self%nsph))
      if (allocated(self%total_area)) deallocate (self%total_area)
      if (allocated(self%total_volume)) deallocate (self%total_volume)
      allocate (self%total_area)
      allocate (self%total_volume)

      !> Ensure Lebedev cache for current num_leb
      call ensure_lebedev_cache(self, error)
      if (allocated(error)) return

      !> Construct full cavity surface
      call setup_iswig_surface( &
         nsph=self%nsph, &
         centers=self%sphxyz, &
         radii=self%radii, &
         cut_a=self%cut_a, &
         cut_f=self%cut_f, &
         oleb=self%cached_oleb, &
         zeta_born=self%cached_swx, &
         ang_grid=self%ang_grid, &
         ang_weight=self%ang_weight, &
         ngrid=self%ngrid, &
         owner=self%owner, &
         grid_xyz=self%xyz, &
         xi=self%xi0, &
         f=self%f, &
         wleb=self%wleb, &
         a=self%a, &
         normal0=self%normal0, &
         v=self%v, &
         numbering=self%numbering, &
         asph=self%asph, &
         total_area=self%total_area, &
         total_volume=self%total_volume, &
         error=error &
         )
      if (allocated(error)) return

   end subroutine update_cavity_iswig

   !> Unified gradient computation for the iSwiG cavity.
   !> Populates self%area_grad(3, nsph) and self%volume_grad(3, nsph)
   !> in a single pass over grid points and switching function derivatives.
   !>
   !> Both gradients share the same expensive inner loop over pairs
   !> (grid point ip, atom jat) for the switching function derivative df/ds.
   !> The area and volume gradients differ only in the weight applied:
   !>   area:   dA/ds  = sum_p  R^2 * w * (df/ds)
   !>   volume: dV/ds  = sum_p  R * w * r_dot_p/3 * (df/ds)  + geom. term
   subroutine compute_gradient_iswig(self, error)
      class(cavity_type_iswig), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: nsph, ip, iat, jat, iaxis
      real(wp) :: weight, zeta, r_own
      real(wp) :: px, py, pz, rx, ry, rz, r_dot_p
      real(wp) :: dx, dy, dz, dvec(3)
      real(wp) :: dfdR, dswitch
      real(wp) :: area_weight, vol_weight

      nsph = size(self%radii)

      if (allocated(self%area_grad)) deallocate (self%area_grad)
      if (allocated(self%volume_grad)) deallocate (self%volume_grad)
      if (allocated(self%xi1_rA)) deallocate (self%xi1_rA)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      if (allocated(self%xyz1_rA)) deallocate (self%xyz1_rA)
      if (allocated(self%v1_rA)) deallocate (self%v1_rA)
      allocate (self%area_grad(3, nsph), source=0.0_wp)
      allocate (self%volume_grad(3, nsph), source=0.0_wp)
      allocate (self%xi1_rA(3, nsph, self%ngrid), source=0.0_wp)
      allocate (self%f1_rA(3, nsph, self%ngrid), source=0.0_wp)
      allocate (self%xyz1_rA(3, 3, nsph, self%ngrid), source=0.0_wp)
      allocate (self%v1_rA(3, nsph, self%ngrid), source=0.0_wp)

      if (self%ngrid <= 0) return

      do ip = 1, self%ngrid
         iat = self%owner(ip)
         r_own = self%radii(iat)
         do iaxis = 1, 3
            self%xyz1_rA(iaxis, iaxis, iat, ip) = 1.0_wp
         end do

         ! Per-point quantities
         weight = self%wleb(ip)
         zeta = self%xi0(ip)
         px = self%xyz(1, ip)
         py = self%xyz(2, ip)
         pz = self%xyz(3, ip)

         ! Outward normal (point relative to owner center) and volume integrand
         rx = px - self%sphxyz(1, iat)
         ry = py - self%sphxyz(2, iat)
         rz = pz - self%sphxyz(3, iat)
         r_dot_p = rx*px + ry*py + rz*pz

         ! Per-grid-point weights for each gradient type:
         !   area:   R^2 w df/dR
         !   volume: R w (n dot c)/3 df/dR
         area_weight = 1.0_wp
         vol_weight = r_dot_p/(3.0_wp*r_own)

         ! Switching function derivative contributions
         do jat = 1, nsph
            if (jat == iat .or. self%radii(jat) == 0.0_wp) cycle

            call swi_pair_dfdr([px, py, pz], self%sphxyz(:, jat), zeta, self%radii(jat), &
                               self%f(ip), dvec, dfdR)
            dx = dvec(1)
            dy = dvec(2)
            dz = dvec(3)

            self%f1_rA(:, iat, ip) = self%f1_rA(:, iat, ip) + dfdR*dvec
            self%f1_rA(:, jat, ip) = self%f1_rA(:, jat, ip) - dfdR*dvec
            dswitch = r_own*r_own*weight*dfdR

            ! Area gradient
            self%area_grad(1, iat) = self%area_grad(1, iat) + dswitch*area_weight*dx
            self%area_grad(2, iat) = self%area_grad(2, iat) + dswitch*area_weight*dy
            self%area_grad(3, iat) = self%area_grad(3, iat) + dswitch*area_weight*dz
            self%area_grad(1, jat) = self%area_grad(1, jat) - dswitch*area_weight*dx
            self%area_grad(2, jat) = self%area_grad(2, jat) - dswitch*area_weight*dy
            self%area_grad(3, jat) = self%area_grad(3, jat) - dswitch*area_weight*dz

            ! Volume gradient (switching function part), accumulated both per
            ! grid point and into the total; contracting v1_rA over the grid
            ! reproduces volume_grad up to summation order.
            self%volume_grad(1, iat) = self%volume_grad(1, iat) + dswitch*vol_weight*dx
            self%volume_grad(2, iat) = self%volume_grad(2, iat) + dswitch*vol_weight*dy
            self%volume_grad(3, iat) = self%volume_grad(3, iat) + dswitch*vol_weight*dz
            self%volume_grad(1, jat) = self%volume_grad(1, jat) - dswitch*vol_weight*dx
            self%volume_grad(2, jat) = self%volume_grad(2, jat) - dswitch*vol_weight*dy
            self%volume_grad(3, jat) = self%volume_grad(3, jat) - dswitch*vol_weight*dz

            self%v1_rA(:, iat, ip) = self%v1_rA(:, iat, ip) &
               & + dswitch*vol_weight*dvec
            self%v1_rA(:, jat, ip) = self%v1_rA(:, jat, ip) &
               & - dswitch*vol_weight*dvec
         end do

         ! Volume geometric term (owner atom only)
         self%volume_grad(1, iat) = self%volume_grad(1, iat) &
            & + r_own*weight*self%f(ip)/3.0_wp*rx
         self%volume_grad(2, iat) = self%volume_grad(2, iat) &
            & + r_own*weight*self%f(ip)/3.0_wp*ry
         self%volume_grad(3, iat) = self%volume_grad(3, iat) &
            & + r_own*weight*self%f(ip)/3.0_wp*rz

         self%v1_rA(:, iat, ip) = self%v1_rA(:, iat, ip) &
            & + r_own*weight*self%f(ip)/3.0_wp*[rx, ry, rz]
      end do

   end subroutine compute_gradient_iswig

   !> Contract a surface adjoint into the nuclear gradient (reverse mode)
   !>
   !> Accumulates `dE/dR_A` for the energy whose surface adjoints `acc` holds,
   !> without ever forming the `(3, nsph, ngrid)` forward Jacobians
   !>
   !> The result is *added* to `gradient`, so several passes can share one 
   !> accumulator
   !>
   !> @param[in]    self     iSwiG cavity instance (must hold an updated grid)
   !> @param[in]    acc      Accumulated surface-observable adjoints
   !> @param[inout] gradient Nuclear-gradient accumulator (3, nsph)
   !> @param[out]   error    Error object, allocated on failure
   subroutine get_surface_gradient_iswig(self, acc, gradient, error)
      !> iSwiG cavity instance
      class(cavity_type_iswig), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Grid, owner and neighbour sphere indices
      integer :: igrid, iat, jat
      !> Effective switching adjoint, area channel folded in
      real(wp) :: w_f_eff
      !> Point-to-center separation and the radial switching derivative
      real(wp) :: dvec(3), dfdr
      !> Gradient contribution of one (point, sphere) pair
      real(wp) :: contribution(3)

      if (.not. acc%is_initialized()) then
         call fatal_error(error, "get_surface_gradient_iswig: accumulator is not initialized")
         return
      end if
      if (size(acc%w_xi) /= self%ngrid) then
         call fatal_error(error, "get_surface_gradient_iswig: accumulator grid size mismatch")
         return
      end if
      if (any(shape(gradient) /= [3, self%nsph])) then
         call fatal_error(error, "get_surface_gradient_iswig: gradient shape mismatch")
         return
      end if
      if (.not. allocated(self%xyz) .or. .not. allocated(self%owner) .or. &
          .not. allocated(self%sphxyz) .or. .not. allocated(self%radii) .or. &
          .not. allocated(self%wleb) .or. .not. allocated(self%xi0) .or. &
          .not. allocated(self%f)) then
         call fatal_error(error, "get_surface_gradient_iswig: cavity surface data are incomplete")
         return
      end if

      ! Geometry-dependent radii would move the width, weight, normal and
      ! curvature channels this routine drops, and the forward path in
      ! [[compute_gradient_iswig]] ignores them just as completely.
      if (allocated(self%radius_model)) then
         if (allocated(self%radius_model%f1_rA)) then
            if (any(self%radius_model%f1_rA /= 0.0_wp)) then
               call fatal_error(error, "get_surface_gradient_iswig: radii models with a "// &
                                "nuclear dependence are not supported")
               return
            end if
         end if
      end if

      if (self%ngrid <= 0) return

      do igrid = 1, self%ngrid
         iat = self%owner(igrid)

         !* --------------------------- Position channel --------------------------- *!
         gradient(:, iat) = gradient(:, iat) + acc%w_xyz(:, igrid)

         !* -------------------------- Switching channel -------------------------- *!
         ! a_i = R_I^2 wleb_i f_i, and neither radius nor Lebedev weight moves,
         ! so the area adjoint enters purely through df_i/dR_A.
         w_f_eff = acc%w_f(igrid) &
                   + acc%w_a(igrid)*self%radii(iat)**2*self%wleb(igrid)
         if (w_f_eff == 0.0_wp) cycle

         do jat = 1, self%nsph
            if (jat == iat .or. self%radii(jat) == 0.0_wp) cycle

            call swi_pair_dfdr(self%xyz(:, igrid), self%sphxyz(:, jat), self%xi0(igrid), &
                               self%radii(jat), self%f(igrid), dvec, dfdr)

            contribution = w_f_eff*dfdr*dvec
            gradient(:, iat) = gradient(:, iat) + contribution
            gradient(:, jat) = gradient(:, jat) - contribution
         end do
      end do

   end subroutine get_surface_gradient_iswig

   !> Ensure Lebedev grid cache is initialized and matches the requested size
   subroutine ensure_lebedev_cache(self, error)
      class(cavity_type_iswig), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      integer :: isize, oleb

      !> iSwiG-supported Lebedev orders (indexing into swig_xi_tab)
      integer :: iswig_order

      ! Precompute constant swig_xi value for this Lebedev order
      real(wp), parameter :: swig_xi_tab(11) = [ &
                             4.865_wp, 4.855_wp, 4.893_wp, 4.901_wp, 4.903_wp, &
                             4.905_wp, 4.906_wp, 4.905_wp, 4.899_wp, 4.907_wp, 4.907_wp]
      integer, parameter :: iswig_grid_sizes(11) = [ &
                            14, 26, 50, 110, 194, 302, 434, 590, 770, 974, 1202]

      ! Map requested num_leb to Lebedev order index
      call lebedev_order_from_num(self%num_leb, oleb, error)
      if (allocated(error)) return

      !> Check if the self%num_leb is available for iswig (xi)
      iswig_order = -1
      do isize = 1, size(iswig_grid_sizes)
         if (self%num_leb == iswig_grid_sizes(isize)) iswig_order = isize
      end do
      if (iswig_order < 0) then
         write (error_unit, "(a,i0)") "[ERROR] Unsupported Lebedev size in iSwiG: ", self%num_leb
         write (error_unit, "(a)") "Supported sizes:"
         write (error_unit, "(8i10)") iswig_grid_sizes(1:8)
         write (error_unit, "(8i10)") iswig_grid_sizes(9:)
         call fatal_error(error, "Unsupported Lebedev size in iSwiG")
         return
      end if

      if (self%cached_num_leb == self%num_leb .and. allocated(self%ang_grid) .and. allocated(self%ang_weight)) then
         return
      end if

      self%cached_num_leb = self%num_leb
      self%cached_oleb = oleb
      self%cached_swx = swig_xi_tab(iswig_order)

      if (allocated(self%ang_grid)) deallocate (self%ang_grid)
      if (allocated(self%ang_weight)) deallocate (self%ang_weight)

      allocate (self%ang_grid(3, self%num_leb))
      allocate (self%ang_weight(self%num_leb))
      call get_angular_grid(self%cached_oleb, self%ang_grid, self%ang_weight, error)
      if (allocated(error)) return

   end subroutine ensure_lebedev_cache

   !> Compute the iswig surface using cached Lebedev data and precomputed factors
   subroutine setup_iswig_surface( &
      nsph, centers, radii, &
      cut_a, cut_f, &
      oleb, zeta_born, ang_grid, ang_weight, &
      ngrid, owner, grid_xyz, xi, f, wleb, a, normal0, v, numbering, asph, &
      total_area, total_volume, error)

      integer, intent(in) :: nsph
      real(wp), intent(in) :: centers(3, nsph)
      real(wp), intent(in) :: radii(nsph)
      real(wp), intent(in) :: cut_a, cut_f
      integer, intent(in) :: oleb
      real(wp), intent(in) :: zeta_born
      real(wp), intent(in) :: ang_grid(:, :)
      real(wp), intent(in) :: ang_weight(:)

      integer, intent(out) :: ngrid
      integer, allocatable, intent(out) :: owner(:)
      real(wp), allocatable, intent(out) :: grid_xyz(:, :)
      real(wp), allocatable, intent(out) :: a(:)
      real(wp), allocatable, intent(out) :: xi(:)
      real(wp), allocatable, intent(out) :: f(:)
      real(wp), allocatable, intent(out) :: wleb(:)
      !> Outward unit normals (3, ngrid)
      real(wp), allocatable, intent(out) :: normal0(:, :)
      !> Point volume elements (ngrid)
      real(wp), allocatable, intent(out) :: v(:)
      integer, allocatable, intent(out) :: numbering(:)
      real(wp), intent(out) :: asph(:)
      real(wp), intent(out) :: total_area
      real(wp), intent(out) :: total_volume
      type(error_type), allocatable, intent(out) :: error

      integer :: iraw, ipt, num_leb, nraw
      real(wp), allocatable :: xyz_raw(:, :), area_raw(:)
      integer, allocatable  :: owner_raw(:)
      real(wp), allocatable :: zeta_raw(:), weight_raw(:), switch_raw(:)
      real(wp) :: rx, ry, rz

      ! Allocate raw (pre-filter) arrays of total size
      num_leb = grid_size(oleb)
      nraw = nsph*num_leb
      allocate (xyz_raw(3, nraw), source=0.0_wp)
      allocate (area_raw(nraw), source=0.0_wp)
      allocate (owner_raw(nraw), source=-1)
      allocate (zeta_raw(nraw), source=0.0_wp)
      allocate (weight_raw(nraw), source=0.0_wp)
      allocate (switch_raw(nraw), source=0.0_wp)

      ! Fill raw arrays
      call fill_intermediate_arrays(nsph, centers, radii, num_leb, ang_grid, &
                                    ang_weight, zeta_born, nraw, xyz_raw, area_raw, owner_raw, &
                                    zeta_raw, weight_raw, switch_raw)

      ! Compute switch_raw(iraw) = product_{j /= owner} [1 - 0.5*(erf(arg_plus)+erf(arg_minus))]
      call compute_switching_function(nraw, nsph, owner_raw, xyz_raw, centers, &
                                      zeta_raw, radii, switch_raw)

      ! Filter out points according to cut_a or cut_f
      ngrid = 0
      do iraw = 1, nraw
         if (cut_a > 0.0_wp) then
            if (switch_raw(iraw)*area_raw(iraw) > cut_a) ngrid = ngrid + 1
         else
            if (switch_raw(iraw) > cut_f) ngrid = ngrid + 1
         end if
      end do

      allocate (grid_xyz(3, ngrid), source=0.0_wp)
      allocate (a(ngrid), source=0.0_wp)
      allocate (owner(ngrid), source=-1)
      allocate (xi(ngrid), source=0.0_wp)
      allocate (f(ngrid), source=0.0_wp)
      allocate (wleb(ngrid), source=0.0_wp)
      allocate (normal0(3, ngrid), source=0.0_wp)
      allocate (v(ngrid), source=0.0_wp)
      allocate (numbering(ngrid), source=-1)

      if (ngrid == 0) then
         call fatal_error(error, "iSwiG: no points left after filtering.")
         return
      end if

      total_area = 0.0_wp
      total_volume = 0.0_wp
      asph = 0.0_wp

      ipt = 0
      do iraw = 1, nraw
         if (cut_a > 0.0_wp) then
            if ((switch_raw(iraw)*area_raw(iraw)) <= cut_a) cycle
         else
            if (switch_raw(iraw) <= cut_f) cycle
         end if

         ipt = ipt + 1
         numbering(ipt) = iraw

         grid_xyz(1, ipt) = xyz_raw(1, iraw)
         grid_xyz(2, ipt) = xyz_raw(2, iraw)
         grid_xyz(3, ipt) = xyz_raw(3, iraw)

         a(ipt) = area_raw(iraw)*switch_raw(iraw)
         total_area = total_area + a(ipt)

         owner(ipt) = owner_raw(iraw)
         xi(ipt) = zeta_raw(iraw)
         f(ipt) = switch_raw(iraw)
         wleb(ipt) = weight_raw(iraw)

         asph(owner(ipt)) = asph(owner(ipt)) + a(ipt)

         rx = grid_xyz(1, ipt) - centers(1, owner(ipt))
         ry = grid_xyz(2, ipt) - centers(2, owner(ipt))
         rz = grid_xyz(3, ipt) - centers(3, owner(ipt))

         ! Points sit on their owner sphere, so the outward unit normal is the
         ! radial direction and the volume element is the divergence-theorem
         ! contribution a_i (r_i . n_i)/3 that the total below accumulates.
         normal0(1, ipt) = rx/radii(owner(ipt))
         normal0(2, ipt) = ry/radii(owner(ipt))
         normal0(3, ipt) = rz/radii(owner(ipt))

         v(ipt) = a(ipt)* &
                  (rx*grid_xyz(1, ipt) + ry*grid_xyz(2, ipt) + rz*grid_xyz(3, ipt))/ &
                  (3.0_wp*radii(owner(ipt)))
         total_volume = total_volume + v(ipt)
      end do

   end subroutine setup_iswig_surface

   !> Fill raw (pre-filter) arrays with surface points and initial values
   subroutine fill_intermediate_arrays( &
      nsph, centers, radii, num_leb, ang_grid, ang_weight, zeta_born, &
      nraw, xyz_raw, area_raw, owner_raw, zeta_raw, weight_raw, switch_raw)
      implicit none

      !> Number of spheres
      integer, intent(in) :: nsph
      !> Coordinates of sphere centers in bohr
      real(wp), intent(in) :: centers(3, nsph)
      !> Radii of spheres in bohr
      real(wp), intent(in) :: radii(nsph)
      !> Number of Lebedev points per sphere
      integer, intent(in) :: num_leb
      !> Unit vectors for Lebedev grid (3, num_leb)
      real(wp), intent(in) :: ang_grid(3, num_leb)
      !> Angular weights (unitless, sum = 4*pi) (num_leb)
      real(wp), intent(in) :: ang_weight(num_leb)
      !> Gaussian-width scale factor for this Lebedev order
      real(wp), intent(in) :: zeta_born
      !> Total raw points = nsph * num_leb
      integer, intent(in) :: nraw

      !> Raw coords (3, nraw)
      real(wp), intent(out) :: xyz_raw(3, nraw)
      !> Raw area before switching (nraw)
      real(wp), intent(out) :: area_raw(nraw)
      !> Raw owner indices (nraw)
      integer, intent(out) :: owner_raw(nraw)
      !> Raw Gaussian widths (nraw)
      real(wp), intent(out) :: zeta_raw(nraw)
      !> Raw Lebedev weights (nraw)
      real(wp), intent(out) :: weight_raw(nraw)
      !> Raw switching function (nraw)
      real(wp), intent(out) :: switch_raw(nraw)

      !> Loop variables
      integer :: iat, ileb, iraw

      iraw = 0
      do iat = 1, nsph
         do ileb = 1, num_leb
            iraw = iraw + 1

            ! Construct raw Lebedev weight from ang_weight(ileb)
            weight_raw(iraw) = ang_weight(ileb)*(4.0_wp*pi)

            ! Cartesian location of point on sphere iat:
            xyz_raw(1, iraw) = centers(1, iat) + radii(iat)*ang_grid(1, ileb)
            xyz_raw(2, iraw) = centers(2, iat) + radii(iat)*ang_grid(2, ileb)
            xyz_raw(3, iraw) = centers(3, iat) + radii(iat)*ang_grid(3, ileb)

            ! "area before switching" = r**2 * raw_weight
            area_raw(iraw) = radii(iat)**2*weight_raw(iraw)

            ! owner-atom index:
            owner_raw(iraw) = iat

            ! zeta = zeta_born / (radii(iat) * sqrt(raw_weight))
            zeta_raw(iraw) = zeta_born/(radii(iat)*sqrt(weight_raw(iraw)))

            ! initialize switching-value to 1.0:
            switch_raw(iraw) = 1.0_wp
         end do
      end do

   end subroutine fill_intermediate_arrays

   !> Compute switching function values for all surface points
   subroutine compute_switching_function( &
      nraw, nsph, owner_raw, xyz_raw, centers, zeta_raw, radii, switch_raw)
      implicit none

      !> Total number of raw points
      integer, intent(in) :: nraw
      !> Number of spheres
      integer, intent(in) :: nsph
      !> Owner indices for each point (nraw)
      integer, intent(in) :: owner_raw(nraw)
      !> Cartesian coordinates of points (3, nraw)
      real(wp), intent(in) :: xyz_raw(3, nraw)
      !> Coordinates of sphere centers (3, nsph)
      real(wp), intent(in) :: centers(3, nsph)
      !> Gaussian widths at each point (nraw)
      real(wp), intent(in) :: zeta_raw(nraw)
      !> Radii of spheres (nsph)
      real(wp), intent(in) :: radii(nsph)
      !> Switching function values (nraw)
      real(wp), intent(inout) :: switch_raw(nraw)

      !> Loop variables
      integer :: iraw, iat
      !> Coordinate differences
      real(wp) :: dx, dy, dz
      !> Distance between point and sphere center
      real(wp) :: dist
      !> Error function arguments
      real(wp) :: arg_plus, arg_minus
      !> Pairwise switching value
      real(wp) :: switch_pair

      do iraw = 1, nraw
         do iat = 1, nsph
            if (iat == owner_raw(iraw)) cycle ! skip self

            dx = xyz_raw(1, iraw) - centers(1, iat)
            dy = xyz_raw(2, iraw) - centers(2, iat)
            dz = xyz_raw(3, iraw) - centers(3, iat)
            dist = sqrt(dx*dx + dy*dy + dz*dz)

            arg_plus = zeta_raw(iraw)*(radii(iat) + dist)
            arg_minus = zeta_raw(iraw)*(radii(iat) - dist)
            switch_pair = 1.0_wp - 0.5_wp*(erf(arg_plus) + erf(arg_minus))

            switch_raw(iraw) = switch_raw(iraw)*switch_pair
         end do
      end do

   end subroutine compute_switching_function

   !> Calculate prefactors needed for switching function derivatives
   !> (for the iSwiG elementary switching function)
   pure subroutine factors_swi_derivs(point, center, zeta, radius, &
                                      dx, dy, dz, dist, arg_plus_sq, arg_minus_sq)

      !> Charge (grid point) position
      real(wp), intent(in) :: point(3)
      !> Sphere center position
      real(wp), intent(in) :: center(3)
      !> Gaussian charge exponent
      real(wp), intent(in) :: zeta
      !> Sphere radius
      real(wp), intent(in) :: radius

      !> Distance components (can be negative)
      real(wp), intent(out) :: dx, dy, dz
      !> Distance between charge and sphere center
      real(wp), intent(out) :: dist
      !> Squared error-function arguments
      real(wp), intent(out) :: arg_plus_sq, arg_minus_sq

      !> Intermediate values
      real(wp) :: arg_plus, arg_minus

      dx = point(1) - center(1)
      dy = point(2) - center(2)
      dz = point(3) - center(3)

      dist = sqrt(dx*dx + dy*dy + dz*dz)

      arg_plus = zeta*(radius + dist)
      arg_minus = zeta*(radius - dist)

      arg_plus_sq = arg_plus*arg_plus
      arg_minus_sq = arg_minus*arg_minus

   end subroutine factors_swi_derivs

   !> Derivative of a switching factor with respect to one sphere center
   !>
   !> `f_i` is the product of the elementary switching factors over all spheres
   !> other than the owner, so its derivative with respect to the separation
   !> from sphere `j` factorizes into `f_i` times the logarithmic derivative of
   !> that one factor. The result is returned as the radial coefficient `dfdr`
   !> and the separation vector `dvec`, with
   !>
   !>   d f_i / d r  =  dfdr * dvec  =  -d f_i / d c_j
   !>
   !> Shared by the forward Jacobian in [[compute_gradient_iswig]] and the
   !> reverse contraction in [[get_surface_gradient_iswig]] so that both paths
   !> evaluate the identical expression.
   !>
   !> @param[in]  point  Grid point position
   !> @param[in]  center Sphere center position
   !> @param[in]  zeta   Gaussian width of the grid point
   !> @param[in]  radius Radius of the sphere
   !> @param[in]  f_val  Switching factor of the grid point
   !> @param[out] dvec   Separation of point and sphere center
   !> @param[out] dfdr   Radial derivative coefficient
   pure subroutine swi_pair_dfdr(point, center, zeta, radius, f_val, dvec, dfdr)
      !> Grid point position
      real(wp), intent(in) :: point(3)
      !> Sphere center position
      real(wp), intent(in) :: center(3)
      !> Gaussian width of the grid point
      real(wp), intent(in) :: zeta
      !> Radius of the sphere
      real(wp), intent(in) :: radius
      !> Switching factor of the grid point
      real(wp), intent(in) :: f_val

      !> Separation of point and sphere center
      real(wp), intent(out) :: dvec(3)
      !> Radial derivative coefficient
      real(wp), intent(out) :: dfdr

      !> Separation components and distance
      real(wp) :: dx, dy, dz, dist
      !> Error-function arguments and the elementary switching factor
      real(wp) :: arg_plus, arg_minus, arg_plus_sq, arg_minus_sq, switch_pair

      call factors_swi_derivs(point, center, zeta, radius, &
                              dx, dy, dz, dist, arg_plus_sq, arg_minus_sq)

      arg_plus = zeta*(radius + dist)
      arg_minus = zeta*(radius - dist)
      switch_pair = 1.0_wp - 0.5_wp*(erf(arg_plus) + erf(arg_minus))

      dfdr = -f_val*zeta/(sqrt(pi)*switch_pair*dist) &
            & *(exp(-arg_plus_sq) - exp(-arg_minus_sq))
      dvec = [dx, dy, dz]

   end subroutine swi_pair_dfdr

end module moist_cavity_iswig
