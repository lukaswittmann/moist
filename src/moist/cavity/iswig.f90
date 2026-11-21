module moist_cavity_iswig
   use moist_output_citations, only: citation
   use moist_model_parameters, only: moist_model_parameters_type
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use mctc_io_structure, only: structure_type
   use mctc_io, only: new
   use mctc_env, only: error_type, fatal_error, wp
   use, intrinsic :: iso_fortran_env, only: output_unit

   use moist_math_grid_angular_grid, only: moist_math_grid_angular_type, new_lebedev_grid
   use moist_math_linalg_geometry, only: cross_product
   use moist_cavity_type, only: cavity_type, list_cavity_fields_base, &
      & add_volume_origin_gradient, update_cavity_base
   use moist_cavity_adjoint, only: cavity_adjoint_type
   use moist_cavity_bornzeta, only: born_zeta_nleb, get_born_zeta
   use moist_cavity_grid_model_type, only: grid_model_type
   use moist_cavity_grid_type, only: cavity_grid_type
   use moist_cavity_grid_frame, only: grid_frame_type
   use moist_context, only: moist_context_type
   use moist_radius_type, only: radius_type
   use moist_channels_fields, only: field_query_type

   implicit none(type, external)
   private
   public :: cavity_type_iswig, new_cavity_iswig

   public :: moist_cavity_iswig_parameters_type

   !> iSwiG construction parameters
   type, extends(moist_model_parameters_type) :: moist_cavity_iswig_parameters_type
      !> Lebedev points per sphere
      integer :: num_leb = 110
      !> Area cutoff
      real(wp) :: cut_a = 0.0_wp
      !> Switching-factor cutoff
      real(wp) :: cut_f = 1.0e-10_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Declare fields for JSON input, output, and printing
      procedure :: register_entries => register_parameter_entries
   end type moist_cavity_iswig_parameters_type

   !> iSwiG cavity state
   type, extends(cavity_type) :: cavity_type_iswig

      !> Construction settings: Lebedev points per sphere, area and switching cutoffs
      type(moist_cavity_iswig_parameters_type) :: param

      !> Raw Lebedev weights (ngrid)
      real(wp), allocatable :: wleb(:)
      !> Accumulated gradient of total area w.r.t. atom positions (3, nsph)
      real(wp), allocatable :: area_grad(:, :)
      !> Accumulated gradient of total volume w.r.t. atom positions (3, nsph)
      real(wp), allocatable :: volume_grad(:, :)
      !> Grid point number (not changed after removing points)
      integer, allocatable :: numbering(:)

      ! Cached Lebedev data (reused across updates)
      !> Lebedev size of the cache, 0 before the first update
      integer :: cached_num_leb = 0
      !> Born-fitted Gaussian width scale of the cached size
      real(wp) :: cached_swx = 0.0_wp
      !> Lebedev nodes in the lab frame (3, num_leb)
      real(wp), allocatable :: ang_grid(:, :)
      !> Solid-angle weights summing to 4*pi (num_leb)
      real(wp), allocatable :: ang_weight(:)
      !> Whether the next update should expose principal curvatures
      logical :: need_curvature = .false.

   contains
      !> Request surface geometry for the next update
      procedure :: configure => configure_cavity_iswig
      procedure :: update => update_cavity_iswig
      procedure :: get_gradient => compute_gradient_iswig
      !> Contract surface-observable adjoints into the nuclear gradient
      procedure :: contract_surface_gradient => get_surface_gradient_iswig
      procedure :: write_csv_debug => write_cavity_csv_debug
      !> Declare the readable results an iSwiG cavity holds
      procedure :: list_fields => list_cavity_fields_iswig
      !> Construction settings
      procedure :: parameters => iswig_parameters
   end type cavity_type_iswig

contains

   !> Enable optional surface geometry for the next update
   !>
   !> @param[in,out] self iSwiG cavity
   !> @param[in] normal Request normals, which iSwiG already stores
   !> @param[in] curvature Request principal curvatures
   !> @param[out] error Configuration error
   subroutine configure_cavity_iswig(self, normal, curvature, error)
      !> iSwiG cavity
      class(cavity_type_iswig), intent(inout) :: self
      !> Request normals, which iSwiG already stores
      logical, intent(in), optional :: normal
      !> Request principal curvatures
      logical, intent(in), optional :: curvature
      !> Configuration error
      type(error_type), allocatable, intent(out) :: error

      if (present(curvature)) self%need_curvature = self%need_curvature .or. curvature
   end subroutine configure_cavity_iswig

   !> Restore compiled parameter defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_parameter_defaults(self)
      class(moist_cavity_iswig_parameters_type), intent(inout) :: self
      type(moist_cavity_iswig_parameters_type) :: defaults

      self%num_leb = defaults%num_leb
      self%cut_a = defaults%cut_a
      self%cut_f = defaults%cut_f
   end subroutine init_parameter_defaults

   !> Declare parameter fields for JSON input, output, and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_parameter_entries(self)
      class(moist_cavity_iswig_parameters_type), intent(inout), target :: self

      call self%register_int_scalar("num_leb", self%num_leb, level=2)
      call self%register_real_scalar("cut_a", self%cut_a)
      call self%register_real_scalar("cut_f", self%cut_f)
   end subroutine register_parameter_entries

   !> Declare the results an iSwiG cavity holds, on top of the generic ones
   !>
   !> @param[in]    self   iSwiG cavity instance
   !> @param[inout] query  Walker collecting or fetching the declarations
   subroutine list_cavity_fields_iswig(self, query)
      !> iSwiG cavity instance
      class(cavity_type_iswig), intent(in) :: self
      !> Walker collecting or fetching the declarations
      type(field_query_type), intent(inout) :: query

      call list_cavity_fields_base(self, query)

      ! Everything below belongs to the surface of the latest build
      if (.not. self%built) return
      call query%add_int_value("num_leb", "Lebedev points per sphere the grid was built with", &
         & self%param%num_leb)
      call query%add_int("numbering", "Stable point id, kept when points are removed (ngrid)", &
         & self%numbering)
      call query%add_real("wleb", "Lebedev quadrature weight per point (ngrid)", self%wleb)

   end subroutine list_cavity_fields_iswig

   !> Construct from parameter values; omission uses compiled defaults
   !>
   !> @param[in,out] self Object to initialize
   !> @param[in] radius_model Atomic radius model to copy
   !> @param[out] error Construction error
   !> @param[in] param Configuration copied by value
   !> @param[in] ctx Borrowed run context; omitted, a model supplies its own
   !> @param[in] grid_model Grid model to copy, e.g. a symmetry model; omitted,
   !>                       every sphere keeps its plain grid
   subroutine new_cavity_iswig(self, radius_model, error, param, ctx, grid_model)
      !> Cavity to initialize
      type(cavity_type_iswig), intent(inout) :: self
      !> Radius model to copy
      class(radius_type), intent(in) :: radius_model
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Configuration; omitted means compiled defaults
      type(moist_cavity_iswig_parameters_type), intent(in), optional :: param
      !> Borrowed run context; omitted, a model supplies its own
      type(moist_context_type), intent(in), target, optional :: ctx
      !> Grid model to copy; omitted means the plain grid
      class(grid_model_type), intent(in), optional :: grid_model
      !> Resolved configuration
      type(moist_cavity_iswig_parameters_type) :: settings

      if (present(param)) settings = param
      nullify (self%ctx)
      if (present(ctx)) self%ctx => ctx
      self%label = "vdW iSwiG"
      self%param = settings
      if (allocated(self%radius_model)) deallocate (self%radius_model)
      allocate (self%radius_model, source=radius_model)
      call self%set_grid_model(grid_model)
      ! A constructed cavity has no structure, no surface and no coupling yet
      if (allocated(self%mol)) deallocate (self%mol)
      self%updated = .false.
      call self%clear_surface()
      call self%clear_couplings()

   end subroutine new_cavity_iswig

   !> Construction settings of the cavity
   !>
   !> @param[in] self iSwiG cavity instance
   function iswig_parameters(self) result(param)
      !> iSwiG cavity instance
      class(cavity_type_iswig), intent(in), target :: self
      !> Settings, valid while the cavity is
      class(moist_model_parameters_type), pointer :: param

      param => self%param

   end function iswig_parameters

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
         write (unit, "(i0,10(',',g0))") ipt, self%numbering(ipt), &
            self%xyz(1, ipt), self%xyz(2, ipt), self%xyz(3, ipt), &
            self%owner(ipt), self%radii(self%owner(ipt)), &
            self%a(ipt), self%wleb(ipt), self%f(ipt)
      end do
      close (unit)

      write (output_unit, "(a,1x,a)") "[Info] Wrote cavity grid to", trim(filename)

   end subroutine write_cavity_csv_debug

   !> Update surface with the current new geometry (iSwiG implementation)
   !>
   !> - every sphere carries the nodes of its grid components (`grid`); the
   !>   plain grid builds bitwise the points c + r u
   !>
   !> @param[in,out] self  iSwiG cavity
   !> @param[in]     mol   Molecular structure
   !> @param[out]    error Error handling
   subroutine update_cavity_iswig(self, mol, error)
      !> iSwiG cavity
      class(cavity_type_iswig), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call update_cavity_base(self, mol, error)
      if (allocated(error)) return
      call self%require_context(error)
      if (allocated(error)) return
      call self%ctx%cite(citation%iswig)
      call self%ctx%cite(citation%iswig_adaptive)

      !> Set number of spheres
      self%nsph = mol%nat

      ! Nuclear derivative arrays belong to the previous geometry until the
      ! caller explicitly requests a fresh gradient build
      if (allocated(self%area_grad)) deallocate (self%area_grad)
      if (allocated(self%volume_grad)) deallocate (self%volume_grad)
      if (allocated(self%xi1_rA)) deallocate (self%xi1_rA)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      if (allocated(self%xyz1_rA)) deallocate (self%xyz1_rA)
      if (allocated(self%v1_rA)) deallocate (self%v1_rA)
      if (allocated(self%k1)) deallocate (self%k1)
      if (allocated(self%k2)) deallocate (self%k2)

      call self%radius_model%update(mol, error)
      if (allocated(error)) return
      if (allocated(self%radii)) deallocate (self%radii)
      self%radii = self%radius_model%f0
      call self%update_grid_model(error)
      if (allocated(error)) return

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

      ! Grid components of every sphere: the plain nodes, or those of the
      ! grid model (rotated once per distinct rotation)
      call self%build_grid(self%ang_grid, self%ang_weight, self%cached_swx, error)
      if (allocated(error)) return

      !> Construct full cavity surface
      call setup_iswig_surface( &
         nsph=self%nsph, &
         centers=self%sphxyz, &
         origin=self%get_volume_origin(), &
         radii=self%radii, &
         cut_a=self%param%cut_a, &
         cut_f=self%param%cut_f, &
         grid=self%grid, &
         ngrid=self%ngrid, &
         owner=self%owner, &
         comp=self%comp, &
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
      self%updated = .true.
      self%built = .true.

      if (self%need_curvature) then
         self%k1 = 1.0_wp/self%radii(self%owner)
         self%k2 = self%k1
      end if

   end subroutine update_cavity_iswig

   !> Unified gradient computation for the iSwiG cavity
   !>
   !> - populates self%area_grad(3, nsph) and self%volume_grad(3, nsph) in a
   !>   single pass over grid points and switching function derivatives
   !> - both gradients share the same expensive inner loop over pairs
   !>   (grid point ip, atom jat) for the switching function derivative df/ds
   !> - they differ only in the weight applied:
   !>   area:   dA/ds  = sum_p  R^2 * w * (df/ds)
   !>   volume: dV/ds  = sum_p  R * w * r_dot_p/3 * (df/ds)  + geom. term
   !> - r_dot_p is measured from the volume origin, whose own motion adds
   !>   dv_p/dc = -a_p n_p/3 through `add_volume_origin_gradient`
   !> - molecular orientation: the sphere rotations follow every atom, so
   !>   every forward array and total gains the frame term
   !>   ([[grid_frame_type]]); a point's arrays take it at the end of its
   !>   own iteration, while they are in cache, the totals once from the
   !>   per-sphere sums. Per point, the derivatives with respect to the
   !>   rotation omega of its sphere (d = p - c_a = r_a n):
   !>   - position: delta p = omega x d
   !>   - switching: delta f = omega . (d x g), g = df/dp, the owner column
   !>     of `f1_rA` before the frame term
   !>   - area: a = r^2 w f, so r^2 w times the switching one
   !>   - volume: v = a n.(p - c)/3, so n.(p - c)/3 times the area one plus
   !>     omega . (a/3) n x (p - c) (n . delta p = 0 since d is parallel to n)
   !>   - widths do not move: xi depends on the radius and weight only
   !>
   !> @param[in,out] self  iSwiG cavity with a built surface
   !> @param[out]    error Error handling
   subroutine compute_gradient_iswig(self, error)
      !> iSwiG cavity
      class(cavity_type_iswig), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: nsph, ip, iat, jat, iaxis, icomp
      real(wp) :: weight, zeta, r_own
      real(wp) :: px, py, pz, rx, ry, rz, r_dot_p
      real(wp) :: dx, dy, dz, dvec(3)
      real(wp) :: dfdR, dswitch
      real(wp) :: area_weight, vol_weight
      !> Volume origin and the per-point volume derivative wrt. it
      real(wp) :: origin(3), dv_dc(3)
      !> Whether the spheres follow a molecular frame, and its Jacobian
      logical :: frame
      type(grid_frame_type) :: frame_jac
      !> Point offset and normal, offset from the volume origin
      real(wp) :: d(3), n(3), po(3)
      !> Changes per rotation of the owner's grid: position (axis, component),
      !> switching, area, volume
      real(wp) :: t_xyz(3, 3), t_f(3), t_a(3), t_v(3)
      !> Area and volume changes per rotation summed per sphere (3, nsph)
      real(wp), allocatable :: s_a(:, :), s_v(:, :)
      !> Area and volume per unit component weight, summed per component
      real(wp), allocatable :: s_wa(:), s_wv(:)
      !> Inverse component weight and the volume element per unit weight
      real(wp) :: w_inv, v_w

      call self%require_built(error)
      if (allocated(error)) return
      nsph = size(self%radii)
      origin = self%get_volume_origin()

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

      frame = .not. self%plain_grid()
      if (frame) then
         call self%get_frame(frame_jac)
         allocate (s_a(3, self%grid%ncomp), s_v(3, self%grid%ncomp), source=0.0_wp)
         allocate (s_wa(self%grid%ncomp), s_wv(self%grid%ncomp), source=0.0_wp)
      end if

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
         r_dot_p = rx*(px - origin(1)) + ry*(py - origin(2)) + rz*(pz - origin(3))

         ! Per-grid point weights for each gradient type:
         !   area:   R^2 w df/dR
         !   volume: R w (n . (p - c))/3 df/dR, c the volume origin
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
            ! reproduces volume_grad up to summation order
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

         ! Volume-origin term, reaching every atom
         dv_dc = -self%a(ip)*self%normal0(:, ip)/3.0_wp
         call add_volume_origin_gradient(dv_dc, self%volume_grad)
         call add_volume_origin_gradient(dv_dc, self%v1_rA(:, :, ip))

         ! Frame term of the point, dense over the atoms
         if (frame) then
            d = [rx, ry, rz]
            n = self%normal0(:, ip)
            po = self%xyz(:, ip) - origin
            ! delta p_c = omega . (d x e_c)
            t_xyz(:, 1) = [0.0_wp, d(3), -d(2)]
            t_xyz(:, 2) = [-d(3), 0.0_wp, d(1)]
            t_xyz(:, 3) = [d(2), -d(1), 0.0_wp]
            t_f = cross_product(d, self%f1_rA(:, iat, ip))
            t_a = r_own*r_own*weight*t_f
            t_v = dot_product(n, po)/3.0_wp*t_a + self%a(ip)/3.0_wp*cross_product(n, po)
            icomp = self%comp(ip)
            call frame_jac%add_vector(icomp, t_xyz, self%xyz1_rA(:, :, :, ip))
            call frame_jac%add_scalar(icomp, t_f, self%f1_rA(:, :, ip))
            call frame_jac%add_scalar(icomp, t_v, self%v1_rA(:, :, ip))
            s_a(:, icomp) = s_a(:, icomp) + t_a
            s_v(:, icomp) = s_v(:, icomp) + t_v
            ! Weight term: f, a and v carry the component weight as a factor
            ! (a point of a weightless component carries nothing to divide)
            if (frame_jac%has_weight .and. self%grid%comp_weight(icomp) > 0.0_wp) then
               w_inv = 1.0_wp/self%grid%comp_weight(icomp)
               v_w = self%a(ip)*dot_product(n, po)/3.0_wp*w_inv
               call frame_jac%add_weight(icomp, self%f(ip)*w_inv, self%f1_rA(:, :, ip))
               call frame_jac%add_weight(icomp, v_w, self%v1_rA(:, :, ip))
               s_wa(icomp) = s_wa(icomp) + self%a(ip)*w_inv
               s_wv(icomp) = s_wv(icomp) + v_w
            end if
         end if
      end do

      if (frame) then
         call frame_jac%add_total(s_a, self%area_grad)
         call frame_jac%add_total(s_v, self%volume_grad)
         if (frame_jac%has_weight) then
            do icomp = 1, self%grid%ncomp
               call frame_jac%add_weight(icomp, s_wa(icomp), self%area_grad)
               call frame_jac%add_weight(icomp, s_wv(icomp), self%volume_grad)
            end do
         end if
      end if

   end subroutine compute_gradient_iswig

   !> Contract a surface adjoint into the nuclear gradient (reverse mode)
   !>
   !> Accumulates `dE/dR_A` for the energy whose surface adjoints `acc` holds,
   !> without ever forming the `(3, nsph, ngrid)` forward Jacobians
   !>
   !> The result is *added* to `gradient`, so several passes can share one
   !> accumulator
   !>
   !> Molecular orientation: the sphere rotations follow the atoms; the same
   !> sweep sums the torque on every sphere, tau_a = sum_p d_p x pbar_p +
   !> n_p x nbar_p (d_p = p - c_a, pbar the total position adjoint including
   !> the switching channels, nbar = `w_n`) per grid component, together
   !> with the weight adjoint wbar = sum_p fbar_p f_p/w (the switching value
   !> is linear in the component weight w); `add_grid_gradient` maps both to
   !> the atoms once; nothing else of a point depends on its component
   !>
   !> @param[in]    self     iSwiG cavity instance (must hold an updated grid)
   !> @param[in]    acc      Accumulated surface-observable adjoints
   !> @param[in,out] gradient Nuclear-gradient accumulator (3, nsph)
   !> @param[out]   error    Error object, allocated on failure
   subroutine get_surface_gradient_iswig(self, acc, gradient, error)
      !> iSwiG cavity instance
      class(cavity_type_iswig), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Grid, owner and neighbour sphere indices, grid component
      integer :: igrid, iat, jat, icomp
      !> Effective switching adjoint, area channel folded in
      real(wp) :: w_f_eff
      !> Point-to-center separation and the radial switching derivative
      real(wp) :: dvec(3), dfdr
      !> Gradient contribution of one (point, sphere) pair
      real(wp) :: contribution(3)
      !> Whether the spheres follow a molecular frame
      logical :: frame
      !> Total position adjoint of a point
      real(wp) :: pbar(3)
      !> Torque per grid component (3, ncomp), not for the plain grid
      real(wp), allocatable :: torque(:, :)
      !> Weight adjoint per grid component (ncomp), not for the plain grid
      real(wp), allocatable :: wbar(:)

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
      ! [[compute_gradient_iswig]] ignores them just as completely
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

      frame = .not. self%plain_grid()
      if (frame) then
         if (.not. allocated(self%normal0)) then
            call fatal_error(error, "get_surface_gradient_iswig: cavity surface data are incomplete")
            return
         end if
         allocate (torque(3, self%grid%ncomp), wbar(self%grid%ncomp), source=0.0_wp)
      end if

      do igrid = 1, self%ngrid
         iat = self%owner(igrid)

         !* --------------------------- Position channel --------------------------- *!
         gradient(:, iat) = gradient(:, iat) + acc%w_xyz(:, igrid)
         if (frame) pbar = acc%w_xyz(:, igrid)

         !* -------------------------- Switching channel -------------------------- *!
         ! a_i = R_I^2 wleb_i f_i, and neither radius nor Lebedev weight moves,
         ! so the area adjoint enters purely through df_i/dR_A
         w_f_eff = acc%w_f(igrid) &
                   + acc%w_a(igrid)*self%radii(iat)**2*self%wleb(igrid)
         if (w_f_eff /= 0.0_wp) then
            do jat = 1, self%nsph
               if (jat == iat .or. self%radii(jat) == 0.0_wp) cycle

               call swi_pair_dfdr(self%xyz(:, igrid), self%sphxyz(:, jat), self%xi0(igrid), &
                                  self%radii(jat), self%f(igrid), dvec, dfdr)

               contribution = w_f_eff*dfdr*dvec
               gradient(:, iat) = gradient(:, iat) + contribution
               gradient(:, jat) = gradient(:, jat) - contribution
               ! The point itself moves the same way as its owner here
               if (frame) pbar = pbar + contribution
            end do
         end if

         !* ---------------------------- Frame torque ----------------------------- *!
         if (frame) then
            icomp = self%comp(igrid)
            torque(:, icomp) = torque(:, icomp) &
               & + cross_product(self%xyz(:, igrid) - self%sphxyz(:, iat), pbar) &
               & + cross_product(self%normal0(:, igrid), acc%w_n(:, igrid))
            if (self%grid%comp_weight(icomp) > 0.0_wp) then
               wbar(icomp) = wbar(icomp) + w_f_eff*self%f(igrid)/self%grid%comp_weight(icomp)
            end if
         end if
      end do

      if (frame) call self%add_grid_gradient(torque, wbar, gradient)

   end subroutine get_surface_gradient_iswig

   !> Ensure Lebedev grid cache is initialized and matches the requested size
   !>
   !> - Sizes without a fitted zeta (bornzeta.f90) are errors here, checked first
   !> - None of the fitted sizes has negative weights; the generator would
   !>   refuse such a rule
   !>
   !> @param[in,out] self   Cavity whose angular cache is refreshed
   !> @param[out]    error  Set for an unsupported or rejected Lebedev size
   subroutine ensure_lebedev_cache(self, error)
      !> Cavity instance
      class(cavity_type_iswig), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Born-fitted Gaussian width scale of this Lebedev size
      real(wp) :: zeta

      !> Whether this Lebedev size has a fitted zeta
      logical :: found

      !> Angular quadrature filling the cache, solid-angle weights summing to 4*pi
      type(moist_math_grid_angular_type) :: leb

      !> Error message buffer
      character(len=512) :: msg

      if (self%cached_num_leb == self%param%num_leb .and. allocated(self%ang_grid) .and. allocated(self%ang_weight)) then
         return
      end if

      ! Gaussian width scale shared with DROP
      call get_born_zeta(self%param%num_leb, zeta, found)
      if (.not. found) then
         write (msg, "(a,i0,a,*(i0,:,', '))") "Unsupported Lebedev size in iSwiG ", self%param%num_leb, &
            & "; supported sizes: ", born_zeta_nleb
         call fatal_error(error, trim(msg))
         return
      end if

      ! Negative weights would give negative areas and imaginary widths
      call new_lebedev_grid(leb, error, npts=self%param%num_leb, positive_weights_only=.true.)
      if (allocated(error)) return

      self%cached_num_leb = self%param%num_leb
      self%cached_swx = zeta

      ! Cache nodes and weights in the plain arrays the integrators consume
      call move_alloc(leb%points, self%ang_grid)
      call move_alloc(leb%weights, self%ang_weight)

   end subroutine ensure_lebedev_cache

   !> Compute the iswig surface using cached Lebedev data and precomputed factors
   !>
   !> - every sphere carries the nodes of its grid components
   !>
   !> @param[in]  nsph         Number of spheres
   !> @param[in]  centers      Sphere centers (3, nsph), bohr
   !> @param[in]  origin       Volume origin, bohr
   !> @param[in]  radii        Sphere radii (nsph), bohr
   !> @param[in]  cut_a        Area cutoff, used when positive
   !> @param[in]  cut_f        Switching cutoff, used otherwise
   !> @param[in]  grid         Grid components of every sphere
   !> @param[out] ngrid        Number of surviving points
   !> @param[out] owner        Sphere of every point (ngrid)
   !> @param[out] comp         Grid component of every point (ngrid)
   !> @param[out] grid_xyz     Points (3, ngrid), bohr
   !> @param[out] xi           Gaussian widths (ngrid)
   !> @param[out] f            Switching factors (ngrid)
   !> @param[out] wleb         Lebedev weights (ngrid)
   !> @param[out] a            Areas (ngrid), bohr**2
   !> @param[out] normal0      Outward unit normals (3, ngrid)
   !> @param[out] v            Volume elements (ngrid), bohr**3
   !> @param[out] numbering    Raw index of every point (ngrid)
   !> @param[out] asph         Area per sphere (nsph), bohr**2
   !> @param[out] total_area   Total area, bohr**2
   !> @param[out] total_volume Total volume, bohr**3
   !> @param[out] error        No point survives the cutoff
   subroutine setup_iswig_surface( &
      nsph, centers, origin, radii, &
      cut_a, cut_f, &
      grid, &
      ngrid, owner, comp, grid_xyz, xi, f, wleb, a, normal0, v, numbering, asph, &
      total_area, total_volume, error)

      !> Number of spheres
      integer, intent(in) :: nsph
      !> Sphere centers, bohr
      real(wp), intent(in) :: centers(3, nsph)
      !> Volume origin, bohr
      real(wp), intent(in) :: origin(3)
      !> Sphere radii, bohr
      real(wp), intent(in) :: radii(nsph)
      !> Area cutoff
      real(wp), intent(in) :: cut_a
      !> Switching cutoff
      real(wp), intent(in) :: cut_f
      !> Grid components of every sphere
      type(cavity_grid_type), intent(in) :: grid

      !> Number of surviving points
      integer, intent(out) :: ngrid
      !> Sphere of every point (ngrid)
      integer, allocatable, intent(out) :: owner(:)
      !> Grid component of every point (ngrid)
      integer, allocatable, intent(out) :: comp(:)
      !> Points (3, ngrid), bohr
      real(wp), allocatable, intent(out) :: grid_xyz(:, :)
      !> Areas (ngrid), bohr**2
      real(wp), allocatable, intent(out) :: a(:)
      !> Gaussian widths (ngrid)
      real(wp), allocatable, intent(out) :: xi(:)
      !> Switching factors (ngrid)
      real(wp), allocatable, intent(out) :: f(:)
      !> Lebedev weights (ngrid)
      real(wp), allocatable, intent(out) :: wleb(:)
      !> Outward unit normals (3, ngrid)
      real(wp), allocatable, intent(out) :: normal0(:, :)
      !> Point volume elements (ngrid)
      real(wp), allocatable, intent(out) :: v(:)
      !> Raw index of every point (ngrid)
      integer, allocatable, intent(out) :: numbering(:)
      !> Area per sphere (nsph), bohr**2
      real(wp), intent(out) :: asph(:)
      !> Total area, bohr**2
      real(wp), intent(out) :: total_area
      !> Total volume, bohr**3
      real(wp), intent(out) :: total_volume
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: iraw, ipt, nraw
      real(wp), allocatable :: xyz_raw(:, :), area_raw(:)
      integer, allocatable  :: owner_raw(:), comp_raw(:)
      real(wp), allocatable :: zeta_raw(:), weight_raw(:), switch_raw(:)
      real(wp) :: rx, ry, rz

      ! Allocate raw (pre-filter) arrays of total size
      nraw = grid%nnode
      allocate (xyz_raw(3, nraw), source=0.0_wp)
      allocate (area_raw(nraw), source=0.0_wp)
      allocate (owner_raw(nraw), source=-1)
      allocate (comp_raw(nraw), source=-1)
      allocate (zeta_raw(nraw), source=0.0_wp)
      allocate (weight_raw(nraw), source=0.0_wp)
      allocate (switch_raw(nraw), source=0.0_wp)

      ! Fill raw arrays
      call fill_intermediate_arrays(nsph, centers, radii, grid, nraw, xyz_raw, area_raw, &
                                    owner_raw, comp_raw, zeta_raw, weight_raw, switch_raw)

      ! Compute switch_raw(iraw) = w_comp product_{j /= owner} 0.5*(erfc(arg_plus)+erfc(arg_minus))
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
      allocate (comp(ngrid), source=-1)
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
         comp(ipt) = comp_raw(iraw)
         xi(ipt) = zeta_raw(iraw)
         f(ipt) = switch_raw(iraw)
         wleb(ipt) = weight_raw(iraw)

         asph(owner(ipt)) = asph(owner(ipt)) + a(ipt)

         rx = grid_xyz(1, ipt) - centers(1, owner(ipt))
         ry = grid_xyz(2, ipt) - centers(2, owner(ipt))
         rz = grid_xyz(3, ipt) - centers(3, owner(ipt))

         ! Points sit on their owner sphere, so the outward unit normal is the
         ! radial direction and the volume element is the divergence-theorem
         ! contribution a_i ((r_i - c) . n_i)/3 that the total below accumulates
         normal0(1, ipt) = rx/radii(owner(ipt))
         normal0(2, ipt) = ry/radii(owner(ipt))
         normal0(3, ipt) = rz/radii(owner(ipt))

         v(ipt) = a(ipt)* &
                  (rx*(grid_xyz(1, ipt) - origin(1)) + ry*(grid_xyz(2, ipt) - origin(2)) &
                   + rz*(grid_xyz(3, ipt) - origin(3)))/ &
                  (3.0_wp*radii(owner(ipt)))
         total_volume = total_volume + v(ipt)
      end do

   end subroutine setup_iswig_surface

   !> Fill raw (pre-filter) arrays with surface points and initial values
   !>
   !> - node j of component c of sphere a sits at c_a + r_a u_j, u_j the
   !>   node of the component's node set
   !> - the switching value starts at the component weight, so the weight
   !>   scales the point like one more switching factor
   !>
   !> @param[in]  nsph       Number of spheres
   !> @param[in]  centers    Sphere centers (3, nsph), bohr
   !> @param[in]  radii      Sphere radii (nsph), bohr
   !> @param[in]  grid       Grid components of every sphere
   !> @param[in]  nraw       Number of nodes over all components
   !> @param[out] xyz_raw    Points (3, nraw), bohr
   !> @param[out] area_raw   Areas before switching (nraw), bohr**2
   !> @param[out] owner_raw  Sphere of every point (nraw)
   !> @param[out] comp_raw   Grid component of every point (nraw)
   !> @param[out] zeta_raw   Gaussian widths (nraw)
   !> @param[out] weight_raw Lebedev weights (nraw)
   !> @param[out] switch_raw Switching factors, the component weights (nraw)
   subroutine fill_intermediate_arrays( &
      nsph, centers, radii, grid, nraw, xyz_raw, area_raw, owner_raw, comp_raw, &
      zeta_raw, weight_raw, switch_raw)
      implicit none(type, external)

      !> Number of spheres
      integer, intent(in) :: nsph
      !> Coordinates of sphere centers in bohr
      real(wp), intent(in) :: centers(3, nsph)
      !> Radii of spheres in bohr
      real(wp), intent(in) :: radii(nsph)
      !> Grid components of every sphere
      type(cavity_grid_type), intent(in) :: grid
      !> Total raw points
      integer, intent(in) :: nraw

      !> Raw coords (3, nraw)
      real(wp), intent(out) :: xyz_raw(3, nraw)
      !> Raw area before switching (nraw)
      real(wp), intent(out) :: area_raw(nraw)
      !> Raw owner indices (nraw)
      integer, intent(out) :: owner_raw(nraw)
      !> Raw component indices (nraw)
      integer, intent(out) :: comp_raw(nraw)
      !> Raw Gaussian widths (nraw)
      real(wp), intent(out) :: zeta_raw(nraw)
      !> Raw Lebedev weights (nraw)
      real(wp), intent(out) :: weight_raw(nraw)
      !> Raw switching function (nraw)
      real(wp), intent(out) :: switch_raw(nraw)

      !> Loop variables
      integer :: iat, icomp, inode, iraw, n0

      do iat = 1, nsph
         do icomp = grid%sph_first(iat), grid%sph_first(iat + 1) - 1
            n0 = grid%set_first(grid%comp_set(icomp)) - grid%comp_first(icomp)
            do iraw = grid%comp_first(icomp), grid%comp_first(icomp + 1) - 1
               inode = n0 + iraw

               ! Raw Lebedev weight, solid angle of the node, rotated or not
               weight_raw(iraw) = grid%wang(inode)

               ! Cartesian location of point on sphere iat:
               xyz_raw(1, iraw) = centers(1, iat) + radii(iat)*grid%node(1, inode)
               xyz_raw(2, iraw) = centers(2, iat) + radii(iat)*grid%node(2, inode)
               xyz_raw(3, iraw) = centers(3, iat) + radii(iat)*grid%node(3, inode)

               ! "area before switching" = r**2 * raw_weight
               area_raw(iraw) = radii(iat)**2*weight_raw(iraw)

               ! owner-atom and component index:
               owner_raw(iraw) = iat
               comp_raw(iraw) = icomp

               ! zeta = zeta_comp / (radii(iat) * sqrt(raw_weight))
               zeta_raw(iraw) = grid%comp_zeta(icomp)/(radii(iat)*sqrt(weight_raw(iraw)))

               ! initialize the switching value to the component weight (1.0
               ! for a single component):
               switch_raw(iraw) = grid%comp_weight(icomp)
            end do
         end do
      end do

   end subroutine fill_intermediate_arrays

   !> Compute switching function values for all surface points
   subroutine compute_switching_function( &
      nraw, nsph, owner_raw, xyz_raw, centers, zeta_raw, radii, switch_raw)
      implicit none(type, external)

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
            switch_pair = 0.5_wp*(erfc(arg_plus) + erfc(arg_minus))

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
   !> that one factor, returned as the radial coefficient `dfdr` and the
   !> separation vector `dvec`, with
   !>
   !>   d f_i / d r  =  dfdr * dvec  =  -d f_i / d c_j
   !>
   !> Shared by the forward Jacobian in [[compute_gradient_iswig]] and the
   !> reverse contraction in [[get_surface_gradient_iswig]] so that both paths
   !> evaluate the identical expression
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
      switch_pair = 0.5_wp*(erfc(arg_plus) + erfc(arg_minus))

      dfdr = -f_val*zeta/(sqrt(pi)*switch_pair*dist) &
            & *(exp(-arg_plus_sq) - exp(-arg_minus_sq))
      dvec = [dx, dy, dz]

   end subroutine swi_pair_dfdr

end module moist_cavity_iswig
