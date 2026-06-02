module moist_cavity_iswig
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use mctc_io_structure, only: structure_type
   use mctc_io, only: new
   use mctc_env, only: error_type, fatal_error, wp
   use iso_fortran_env, only: error_unit, output_unit

   use moist_math_grid_lebedev, only: get_angular_grid, grid_size, lebedev_order_from_num
   use moist_type, only: cavity_type
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

      ! Spheres of system
      !> Number of spheres (atoms) in the cavity
      integer :: nsph
      !> Coordinates of sphere centers (3, nsph), in Bohr
      real(wp), allocatable :: sphxyz(:, :)

      ! Surface properties
      !> Gaussian-widths at each point (ngrid)
      real(wp), allocatable :: xi(:)
      !> Switching function (ngrid)
      real(wp), allocatable :: f(:)
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
      procedure :: contract_amat1_q1q2_rA => contract_amat1_q1q2_rA_iswig
      procedure :: write_csv_debug => write_cavity_csv_debug
      !> ISWIG-specific matrix assembly (symmetric, diagonally dominant)
      procedure :: get_amat => get_amat_iswig
   end type cavity_type_iswig

contains

   !> Constructor for iSwiG cavity
   !> Initialize an already-declared object; no allocation of the object itself.
   subroutine new_cavity_iswig(self, nleb, cut_a, cut_f, radius_model, error)
      !> Cavity type instance to initialize
      type(cavity_type_iswig), intent(inout) :: self
      !> Number of lebedev grid points per unit sphere
      integer, intent(in), optional :: nleb
      !> Settings for iSwiG cavity
      real(wp), intent(in), optional :: cut_a, cut_f
      !> Enable simplified mode
      !> Optional radii model
      class(radius_type), intent(in) :: radius_model
      !> Constructor error
      type(error_type), allocatable, intent(out) :: error

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

      open (file=filename, newunit=unit, status='replace', action='write', iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, 'Could not open CSV file for writing: '//trim(filename))
         return
      end if

      write (unit, '(a)') 'ngrid,numbering,x,y,z,owner,radius,area,w_leb,f'

      do ipt = 1, self%ngrid
         write (unit, '(i0,10('','',g0))') ipt, self%numbering(ipt), &
            self%xyz(1, ipt), self%xyz(2, ipt), self%xyz(3, ipt), &
            self%owner(ipt), self%radii(self%owner(ipt)), &
            self%a(ipt), self%wleb(ipt), self%f(ipt)
      end do
      close (unit)

      write (output_unit, '(a,1x,a)') '[Info] Wrote cavity grid to', trim(filename)

   end subroutine write_cavity_csv_debug

   !> Update surface with the current new geometry (iSwiG implementation)
   subroutine update_cavity_iswig(self, mol, error)
      class(cavity_type_iswig), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error

      !> Set number of spheres
      self%nsph = mol%nat

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
         xi=self%xi, &
         f=self%f, &
         wleb=self%wleb, &
         a=self%a, &
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
   subroutine compute_gradient_iswig(self)
      class(cavity_type_iswig), intent(inout) :: self

      integer :: nsph, ip, iat, jat
      real(wp) :: weight, zeta, r_own
      real(wp) :: px, py, pz, rx, ry, rz, r_dot_p
      real(wp) :: dx, dy, dz, dist, arg_plus, arg_minus, arg_plus_sq, arg_minus_sq
      real(wp) :: switch_pair
      real(wp) :: pref, pref_zeta, pref_pair, dswitch
      real(wp) :: area_weight, vol_weight

      nsph = size(self%radii)

      if (.not. allocated(self%area_grad)) allocate (self%area_grad(3, nsph))
      if (.not. allocated(self%volume_grad)) allocate (self%volume_grad(3, nsph))
      self%area_grad = 0.0_wp
      self%volume_grad = 0.0_wp

      if (self%ngrid <= 0) return

      do ip = 1, self%ngrid
         iat = self%owner(ip)
         r_own = self%radii(iat)

         ! Per-point quantities
         weight = self%wleb(ip)
         zeta = self%xi(ip)
         px = self%xyz(1, ip)
         py = self%xyz(2, ip)
         pz = self%xyz(3, ip)

         ! Outward normal (point relative to owner center) and volume integrand
         rx = px - self%sphxyz(1, iat)
         ry = py - self%sphxyz(2, iat)
         rz = pz - self%sphxyz(3, iat)
         r_dot_p = rx*px + ry*py + rz*pz

         ! Common base factor for switching function derivative
         pref = -r_own*r_own*weight*self%f(ip)/sqrt(pi)
         pref_zeta = pref*zeta

         ! Per-grid-point weights for each gradient type
         area_weight = 1.0_wp
         vol_weight = r_dot_p/(3.0_wp*r_own)

         ! Switching function derivative contributions
         do jat = 1, nsph
            if (jat == iat .or. self%radii(jat) == 0.0_wp) cycle

            call factors_swi_derivs([px, py, pz], self%sphxyz(:, jat), zeta, self%radii(jat), &
                                    dx, dy, dz, dist, arg_plus_sq, arg_minus_sq)

            arg_plus = zeta*(self%radii(jat) + dist)
            arg_minus = zeta*(self%radii(jat) - dist)
            switch_pair = 1.0_wp - 0.5_wp*(erf(arg_plus) + erf(arg_minus))

            pref_pair = pref_zeta/(switch_pair*dist)
            dswitch = pref_pair*(exp(-arg_plus_sq) - exp(-arg_minus_sq))

            ! Area gradient
            self%area_grad(1, iat) = self%area_grad(1, iat) + dswitch*area_weight*dx
            self%area_grad(2, iat) = self%area_grad(2, iat) + dswitch*area_weight*dy
            self%area_grad(3, iat) = self%area_grad(3, iat) + dswitch*area_weight*dz
            self%area_grad(1, jat) = self%area_grad(1, jat) - dswitch*area_weight*dx
            self%area_grad(2, jat) = self%area_grad(2, jat) - dswitch*area_weight*dy
            self%area_grad(3, jat) = self%area_grad(3, jat) - dswitch*area_weight*dz

            ! Volume gradient (switching function part)
            self%volume_grad(1, iat) = self%volume_grad(1, iat) + dswitch*vol_weight*dx
            self%volume_grad(2, iat) = self%volume_grad(2, iat) + dswitch*vol_weight*dy
            self%volume_grad(3, iat) = self%volume_grad(3, iat) + dswitch*vol_weight*dz
            self%volume_grad(1, jat) = self%volume_grad(1, jat) - dswitch*vol_weight*dx
            self%volume_grad(2, jat) = self%volume_grad(2, jat) - dswitch*vol_weight*dy
            self%volume_grad(3, jat) = self%volume_grad(3, jat) - dswitch*vol_weight*dz
         end do

         ! Volume geometric term (owner atom only)
         self%volume_grad(1, iat) = self%volume_grad(1, iat) &
            & + r_own*weight*self%f(ip)/3.0_wp*rx
         self%volume_grad(2, iat) = self%volume_grad(2, iat) &
            & + r_own*weight*self%f(ip)/3.0_wp*ry
         self%volume_grad(3, iat) = self%volume_grad(3, iat) &
            & + r_own*weight*self%f(ip)/3.0_wp*rz
      end do

   end subroutine compute_gradient_iswig

   !> Compute the contracted A-matrix derivative: grad = q1^T (dA/dR) q2.
   !> Returns grad(3, nsph) without forming the full derivative tensor.
   !>
   !> The A-matrix depends on positions through:
   !>  (a) area elements a_p = R^2 * w * f_p (switching function)
   !>  (b) grid point positions c_p (move with owner atom)
   !>
   !> @param[in]  q1    first charge vector (ngrid)
   !> @param[in]  q2    second charge vector (ngrid)
   !> @param[out] grad  contracted gradient (3, nsph)
   !> @param[out] error error handling
   subroutine contract_amat1_q1q2_rA_iswig(self, q1, q2, grad, error)
      class(cavity_type_iswig), intent(in) :: self
      !> First charge vector (ngrid)
      real(wp), intent(in) :: q1(:)
      !> Second charge vector (ngrid)
      real(wp), intent(in) :: q2(:)
      !> Contracted gradient (3, nsph)
      real(wp), intent(out) :: grad(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: nsph, ngrid, ip, jp, iat, jat
      real(wp) :: px, py, pz, zeta, weight, r_own, switch
      real(wp) :: dx, dy, dz, dist
      real(wp) :: arg_plus, arg_minus, arg_plus_sq, arg_minus_sq, switch_pair
      real(wp) :: pref, pref_zeta, pref_pair, dswitch
      real(wp) :: r_vec(3), r_dist, r_dist3
      real(wp) :: sqrt_ap, sqrt_aq, dA_dc, pair_weight
      real(wp), parameter :: fourpi = 4.0_wp*pi
      real(wp), parameter :: iswig_factor = 5.0_wp
      real(wp), allocatable :: dform_da(:)
      real(wp) :: dA_da

      nsph = size(self%radii)
      ngrid = self%ngrid

      grad = 0.0_wp
      if (ngrid <= 0) return

      ! Part 1: Precompute, for each point, d(q1^T A q2) / d(a_p)
      allocate (dform_da(ngrid))
      dform_da = 0.0_wp

      do ip = 1, ngrid
         sqrt_ap = sqrt(self%a(ip))
         ! Off-diagonal contributions to dform_da(ip)
         do jp = 1, ngrid
            if (jp == ip) cycle
            sqrt_aq = sqrt(self%a(jp))
            r_vec = self%xyz(:, ip) - self%xyz(:, jp)
            r_dist = sqrt(sum(r_vec**2))
            ! A(p,q) = -sqrt_ap * sqrt_aq / (fourpi * r_dist)
            ! dA(p,q)/da_p = A(p,q) / (2*a_p)
            !              = -sqrt_aq / (2 * sqrt_ap * fourpi * r_dist)
            dA_da = -sqrt_aq/(2.0_wp*sqrt_ap*fourpi*r_dist)
            dform_da(ip) = dform_da(ip) + q1(ip)*dA_da*q2(jp) + q1(jp)*dA_da*q2(ip)
         end do
         ! Diagonal contribution to dform_da(ip)
         do jp = 1, ngrid
            if (jp == ip) cycle
            sqrt_aq = sqrt(self%a(jp))
            r_vec = self%xyz(:, ip) - self%xyz(:, jp)
            r_dist = sqrt(sum(r_vec**2))
            dA_da = -sqrt_aq/(2.0_wp*sqrt_ap*fourpi*r_dist)
            ! d(iswig_factor * |A(p,q)|)/da_p contributes to A(p,p) diagonal
            dform_da(ip) = dform_da(ip) + q1(ip)*(-iswig_factor*dA_da)*q2(ip)
            ! A(q,q) also depends on a_p through iswig row-sum: |A(q,p)|
            dform_da(ip) = dform_da(ip) + q1(jp)*(-iswig_factor*dA_da)*q2(jp)
         end do
         ! d(2*pi/a_p)/da_p = -2*pi/a_p^2
         dform_da(ip) = dform_da(ip) &
            & + q1(ip)*(-2.0_wp*pi/self%a(ip)**2)*q2(ip)
      end do

      do ip = 1, ngrid
         iat = self%owner(ip)
         r_own = self%radii(iat)
         weight = self%wleb(ip)
         zeta = self%xi(ip)
         switch = self%f(ip)
         px = self%xyz(1, ip)
         py = self%xyz(2, ip)
         pz = self%xyz(3, ip)

         if (abs(dform_da(ip)) < epsilon(1.0_wp)) cycle

         ! pref = -R^2 * w * f / sqrt(pi), pref_zeta = pref * zeta
         pref = -r_own*r_own*weight*switch/sqrt(pi)
         pref_zeta = pref*zeta

         do jat = 1, nsph
            if (jat == iat .or. self%radii(jat) == 0.0_wp) cycle

            call factors_swi_derivs([px, py, pz], self%sphxyz(:, jat), zeta, self%radii(jat), &
                                    dx, dy, dz, dist, arg_plus_sq, arg_minus_sq)

            arg_plus = zeta*(self%radii(jat) + dist)
            arg_minus = zeta*(self%radii(jat) - dist)
            switch_pair = 1.0_wp - 0.5_wp*(erf(arg_plus) + erf(arg_minus))

            pref_pair = pref_zeta/(switch_pair*dist)
            dswitch = pref_pair*(exp(-arg_plus_sq) - exp(-arg_minus_sq))

            ! dswitch * direction = R^2 * w * (partial df_p / partial s)
            ! Weight by dform_da(ip) to get contribution to grad
            grad(1, iat) = grad(1, iat) + dswitch*dform_da(ip)*dx
            grad(2, iat) = grad(2, iat) + dswitch*dform_da(ip)*dy
            grad(3, iat) = grad(3, iat) + dswitch*dform_da(ip)*dz

            grad(1, jat) = grad(1, jat) - dswitch*dform_da(ip)*dx
            grad(2, jat) = grad(2, jat) - dswitch*dform_da(ip)*dy
            grad(3, jat) = grad(3, jat) - dswitch*dform_da(ip)*dz
         end do
      end do

      deallocate (dform_da)

      ! Part 2: Derivative through grid point positions c_p
      do ip = 1, ngrid
         iat = self%owner(ip)
         sqrt_ap = sqrt(self%a(ip))
         do jp = ip + 1, ngrid
            sqrt_aq = sqrt(self%a(jp))

            r_vec = self%xyz(:, ip) - self%xyz(:, jp)
            r_dist = sqrt(sum(r_vec**2))
            r_dist3 = r_dist*r_dist*r_dist

            dA_dc = sqrt_ap*sqrt_aq/(fourpi*r_dist3)

            pair_weight = (q1(ip)*q2(jp) + q1(jp)*q2(ip)) &
               & - iswig_factor*(q1(ip)*q2(ip) + q1(jp)*q2(jp))

            grad(:, iat) = grad(:, iat) + pair_weight*dA_dc*r_vec
            grad(:, self%owner(jp)) = grad(:, self%owner(jp)) - pair_weight*dA_dc*r_vec
         end do
      end do

   end subroutine contract_amat1_q1q2_rA_iswig

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
         write (error_unit, '(a,i0)') "[ERROR] Unsupported Lebedev size in iSwiG: ", self%num_leb
         write (error_unit, '(a)') "Supported sizes:"
         write (error_unit, '(8i10)') iswig_grid_sizes(1:8)
         write (error_unit, '(8i10)') iswig_grid_sizes(9:)
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
      ngrid, owner, grid_xyz, xi, f, wleb, a, numbering, asph, total_area, total_volume, error)

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
         total_volume = total_volume + a(ipt)* &
                        (rx*grid_xyz(1, ipt) + ry*grid_xyz(2, ipt) + rz*grid_xyz(3, ipt))/ &
                        (3.0_wp*radii(owner(ipt)))
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
   !> (for the ISWIG elementary switching function)
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

   !> Assemble PCM interaction matrix using ISWIG-style assembly.
   !> Builds a symmetric, diagonally dominant matrix suitable for iterative solvers.
   !> Off-diagonal: symmetric area weighting with geometric mean sqrt(a_i * a_j).
   !> Diagonal: self-potential + row-sum enhancement for positive definiteness.
   !> @param[out] amat  Interaction matrix (ngrid, ngrid)
   !> @param[out] error Error handling
   subroutine get_amat_iswig(self, amat, error)
      !> iSwiG cavity instance
      class(cavity_type_iswig), intent(in) :: self
      !> Output: assembled matrix (ngrid, ngrid)
      real(wp), intent(out) :: amat(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: ip, jp, ngrid
      real(wp) :: r_vec(3), r_dist, row_sum
      real(wp), parameter :: fourpi = 4.0_wp*pi
      real(wp), parameter :: iswig_factor = 5.0_wp

      ngrid = self%ngrid

      ! Check dimensions
      if (size(amat, 1) /= ngrid .or. size(amat, 2) /= ngrid) then
         call fatal_error(error, &
            & "[get_amat_iswig] Matrix dimension mismatch")
         return
      end if

      ! Build ISWIG-style matrix (symmetric, diagonally dominant)
      do ip = 1, ngrid
         row_sum = 0.0_wp
         do jp = 1, ngrid
            if (ip /= jp) then
               r_vec(:) = self%xyz(:, ip) - self%xyz(:, jp)
               r_dist = sqrt(sum(r_vec**2))
               amat(ip, jp) = -sqrt(self%a(ip)*self%a(jp))/(fourpi*r_dist)
               row_sum = row_sum + abs(amat(ip, jp))
            end if
         end do
         amat(ip, ip) = (2.0_wp*pi)/self%a(ip) + iswig_factor*row_sum
      end do

   end subroutine get_amat_iswig

end module moist_cavity_iswig
