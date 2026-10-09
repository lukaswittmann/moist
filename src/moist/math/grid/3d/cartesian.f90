!> Cartesian 3D grid and 3D Fourier transforms for 3D-RISM
module moist_math_grid_3d_cartesian
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type, moist_math_grid_3d_trafo_type, &
                                      integrand_3d, check_trafo_blocks
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_fft, only: moist_fft_r2c_3d, moist_fft_c2r_3d, &
      & moist_fft_r2c_3d_batch, moist_fft_c2r_3d_batch
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
   implicit none(type, external)
   private

   public :: moist_math_grid_3d_cartesian_type
   public :: new_cartesian_point_grid, new_cartesian_gaussian_grid
   public :: moist_math_grid_3d_cartesian_trafo_type
   public :: nk_x_from_nx

   !> Default real-space points per axis
   integer, parameter :: default_box_points = 128
   !> Default real-space spacing, bohr
   real(wp), parameter :: default_box_spacing = 0.6_wp
   !> Default minimum distance from every nucleus to the box faces, bohr
   real(wp), parameter :: default_box_margin = 10.0_wp

   !> Cartesian volume domain with a centroid-following origin
   !>
   !> No transform handles; destroy borrowed engines before geometry mutation
   type, extends(moist_math_grid_3d_type) :: moist_math_grid_3d_cartesian_type
      !> Publish Gaussian widths
      logical, private :: gaussian = .false.
      !> Molecular centroid of the last successful update, bohr
      real(wp), private :: center(3) = 0.0_wp
      !> Real-space size along x
      integer :: nx = default_box_points
      !> Real-space size along y
      integer :: ny = default_box_points
      !> Real-space size along z
      integer :: nz = default_box_points
      !> First reciprocal-space dimension (R2C halved axis)
      integer, private :: nkx = 0
      !> Real-space spacing (bohr)
      real(wp) :: dr = default_box_spacing
      !> Minimum distance from every nucleus to the box faces (bohr); keeps the
      !> solute away from its periodic images
      real(wp) :: margin = default_box_margin
      !> Gaussian exponent scale, xi0 = xi0_factor/dr
      real(wp) :: xi0_factor = 1.0_wp
      !> Real-space cell volume (dr**3)
      real(wp), private :: dv = 0.0_wp
      !> Total box volume
      real(wp), private :: vbox = 0.0_wp
      !> Reciprocal-space spacing along x, inverse bohr
      real(wp), private :: dkx = 0.0_wp
      !> Reciprocal-space spacing along y, inverse bohr
      real(wp), private :: dky = 0.0_wp
      !> Reciprocal-space spacing along z, inverse bohr
      real(wp), private :: dkz = 0.0_wp
      !> Real-space coordinate of grid point (1, 1, 1)
      real(wp), private :: origin(3) = 0.0_wp
      !> Reciprocal frequencies along x, inverse bohr, length nkx
      real(wp), allocatable, private :: kx(:)
      !> Reciprocal frequencies along y, inverse bohr, length ny, FFT-ordered
      real(wp), allocatable, private :: ky(:)
      !> Reciprocal frequencies along z, inverse bohr, length nz, FFT-ordered
      real(wp), allocatable, private :: kz(:)
   contains
      procedure :: validate => validate_cartesian_grid
      procedure :: update => update_cartesian
      procedure :: get_volume_gradient => get_volume_gradient_cartesian
      procedure :: get_volume_hessian_vector => get_volume_hessian_vector_cartesian
      procedure :: get_volume_hessian => get_volume_hessian_cartesian
      procedure :: kind_name => cartesian_kind_name
      procedure :: kpoint => cartesian_grid_3d_kpoint
      procedure :: destroy => cartesian_grid_3d_dealloc
      procedure :: integrate_field => cartesian_grid_3d_integrate_field
      procedure :: new_trafo => cartesian_grid_3d_new_trafo
      !> Real-space coordinate of grid point (1, 1, 1), bohr
      procedure :: get_origin => cartesian_get_origin
      !> Real-space cell volume `dr**3`, bohr**3
      procedure :: get_cell_volume => cartesian_get_cell_volume
      !> Periodic box volume, bohr**3
      procedure :: get_box_volume => cartesian_get_box_volume
   end type moist_math_grid_3d_cartesian_type

   !> 3D Fourier transform engine for a Cartesian grid
   !>
   !> Stateless apart from the grid pointer; the backend transforms the
   !> caller's buffers directly, which is thread-safe as long as the buffers
   !> are distinct; a single instance may therefore be shared across
   !> threads, but one per thread is equally fine
   type, extends(moist_math_grid_3d_trafo_type) :: moist_math_grid_3d_cartesian_trafo_type
      !> Grid this trafo transforms on (not owned; must outlive the trafo)
      class(moist_math_grid_3d_cartesian_type), pointer, private :: grid => null()
   contains
      procedure :: fft_r2k => cartesian_trafo_fft_r2k
      procedure :: fft_k2r => cartesian_trafo_fft_k2r
      procedure :: destroy => cartesian_trafo_destroy
   end type moist_math_grid_3d_cartesian_trafo_type

contains

   !> Configure a Cartesian point grid without realizing geometry
   !>
   !> @param[in,out] self Grid configuration
   !> @param[out] error Invalid settings
   !> @param[in] nx Optional x point count, default 128
   !> @param[in] ny Optional y point count, default 128
   !> @param[in] nz Optional z point count, default 128
   !> @param[in] dr Optional spacing, default 0.6 bohr
   !> @param[in] margin Optional nucleus-to-face distance, default 10 bohr
   subroutine new_cartesian_point_grid(self, error, nx, ny, nz, dr, margin)
      !> Grid configuration
      type(moist_math_grid_3d_cartesian_type), intent(inout) :: self
      !> Configuration error
      type(error_type), allocatable, intent(out) :: error
      !> X point count
      integer, intent(in), optional :: nx
      !> Y point count
      integer, intent(in), optional :: ny
      !> Z point count
      integer, intent(in), optional :: nz
      !> Spacing
      real(wp), intent(in), optional :: dr
      !> Nucleus-to-face distance
      real(wp), intent(in), optional :: margin

      call configure_cartesian_grid(self, error, .false., nx, ny, nz, dr, margin=margin)
   end subroutine new_cartesian_point_grid

   !> Configure a Cartesian gaussian grid without realizing geometry
   !>
   !> @param[in,out] self Grid configuration
   !> @param[out] error Invalid settings
   !> @param[in] nx Optional x point count, default 128
   !> @param[in] ny Optional y point count, default 128
   !> @param[in] nz Optional z point count, default 128
   !> @param[in] dr Optional spacing, default 0.6 bohr
   !> @param[in] xi0_factor Optional width scale, default 1
   !> @param[in] margin Optional nucleus-to-face distance, default 10 bohr
   subroutine new_cartesian_gaussian_grid(self, error, nx, ny, nz, dr, xi0_factor, margin)
      !> Grid configuration
      type(moist_math_grid_3d_cartesian_type), intent(inout) :: self
      !> Configuration error
      type(error_type), allocatable, intent(out) :: error
      !> X point count
      integer, intent(in), optional :: nx
      !> Y point count
      integer, intent(in), optional :: ny
      !> Z point count
      integer, intent(in), optional :: nz
      !> Spacing
      real(wp), intent(in), optional :: dr
      !> Gaussian width scale
      real(wp), intent(in), optional :: xi0_factor
      !> Nucleus-to-face distance
      real(wp), intent(in), optional :: margin

      call configure_cartesian_grid(self, error, .true., nx, ny, nz, dr, xi0_factor, margin)
   end subroutine new_cartesian_gaussian_grid

   !> Shared configuration reset for Cartesian variants
   !>
   !> @param[in,out] self Grid configuration
   !> @param[out] error Invalid settings
   !> @param[in] gaussian Publish Gaussian widths
   !> @param[in] nx Optional x point count
   !> @param[in] ny Optional y point count
   !> @param[in] nz Optional z point count
   !> @param[in] dr Optional spacing, bohr
   !> @param[in] xi0_factor Optional Gaussian width scale
   !> @param[in] margin Optional nucleus-to-face distance, bohr
   subroutine configure_cartesian_grid(self, error, gaussian, nx, ny, nz, dr, xi0_factor, margin)
      !> Grid configuration
      type(moist_math_grid_3d_cartesian_type), intent(inout) :: self
      !> Configuration error
      type(error_type), allocatable, intent(out) :: error
      !> Gaussian variant
      logical, intent(in) :: gaussian
      !> X point count
      integer, intent(in), optional :: nx
      !> Y point count
      integer, intent(in), optional :: ny
      !> Z point count
      integer, intent(in), optional :: nz
      !> Spacing
      real(wp), intent(in), optional :: dr
      !> Gaussian width scale
      real(wp), intent(in), optional :: xi0_factor
      !> Nucleus-to-face distance
      real(wp), intent(in), optional :: margin

      call self%destroy()
      self%nx = default_box_points
      self%ny = default_box_points
      self%nz = default_box_points
      self%dr = default_box_spacing
      self%margin = default_box_margin
      self%xi0_factor = 1.0_wp
      self%gaussian = gaussian
      self%center = 0.0_wp
      if (present(nx)) self%nx = nx
      if (present(ny)) self%ny = ny
      if (present(nz)) self%nz = nz
      if (present(dr)) self%dr = dr
      if (present(xi0_factor)) self%xi0_factor = xi0_factor
      if (present(margin)) self%margin = margin
      call check_cartesian_settings(self, error)
   end subroutine configure_cartesian_grid

   !> Contract centroid-following point adjoints into the nuclear gradient
   !>
   !> Every point moves with the centroid and weights and widths are constant,
   !> so each atom receives `sum(w_xyz, dim=2)/natom`
   !>
   !> @param[in] self Successfully updated grid
   !> @param[in] acc Volume-observable adjoints
   !> @param[in,out] gradient Nuclear-gradient accumulator (3, natom)
   !> @param[out] error Invalid geometry or adjoints
   subroutine get_volume_gradient_cartesian(self, acc, gradient, error)
      !> Updated grid
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Volume adjoints
      type(volume_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Invalid input
      type(error_type), allocatable, intent(out) :: error

      !> Contribution of every atom
      real(wp) :: contribution(3)
      !> Atom index
      integer :: iat

      call self%check_volume_adjoint(acc, error, vector=gradient)
      if (allocated(error)) return
      contribution = sum(acc%w_xyz, dim=2)/real(self%natom, wp)
      if (.not. all(ieee_is_finite(contribution))) then
         call fatal_error(error, "cartesian domain: non-finite gradient contraction")
         return
      end if
      do iat = 1, self%natom
         gradient(:, iat) = gradient(:, iat) + contribution
      end do
   end subroutine get_volume_gradient_cartesian

   !> Zero curvature of linear centroid motion and constant quadrature weights
   !>
   !> @param[in] self Successfully updated grid
   !> @param[in] acc Fixed volume adjoints
   !> @param[in] direction Nuclear displacement direction (3, natom)
   !> @param[in,out] hessian_vector Nuclear Hessian-vector accumulator (3, natom)
   !> @param[out] error Invalid geometry, direction or adjoints
   subroutine get_volume_hessian_vector_cartesian(self, acc, direction, hessian_vector, error)
      !> Updated grid
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Fixed volume adjoints
      type(volume_adjoint_type), intent(in) :: acc
      !> Nuclear displacement direction
      real(wp), intent(in) :: direction(:, :)
      !> Nuclear Hessian-vector accumulator
      real(wp), intent(inout) :: hessian_vector(:, :)
      !> Invalid input
      type(error_type), allocatable, intent(out) :: error

      ! Exact zero contribution: the accumulator stays as it is
      call self%check_volume_adjoint(acc, error, vector=hessian_vector, direction=direction)
   end subroutine get_volume_hessian_vector_cartesian

   !> Zero curvature of linear centroid motion and constant quadrature weights
   !>
   !> @param[in] self Successfully updated grid
   !> @param[in] acc Fixed volume adjoints
   !> @param[in,out] hessian Nuclear Hessian accumulator (3, natom, 3, natom)
   !> @param[out] error Invalid geometry, direction or adjoints
   subroutine get_volume_hessian_cartesian(self, acc, hessian, error)
      !> Updated grid
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Fixed volume adjoints
      type(volume_adjoint_type), intent(in) :: acc
      !> Nuclear Hessian accumulator
      real(wp), intent(inout) :: hessian(:, :, :, :)
      !> Invalid input
      type(error_type), allocatable, intent(out) :: error

      ! Exact zero contribution: the accumulator stays as it is
      call self%check_volume_adjoint(acc, error, hessian=hessian)
   end subroutine get_volume_hessian_cartesian

   !> Validate a configured Cartesian grid
   !>
   !> @param[in] self Grid configuration
   !> @param[out] error Invalid settings
   subroutine validate_cartesian_grid(self, error)
      !> Grid configuration
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Configuration error
      type(error_type), allocatable, intent(out) :: error

      call check_cartesian_settings(self, error)
   end subroutine validate_cartesian_grid

   !> Reject invalid sizes and non-finite scales before grid allocation
   !>
   !> @param[in] self Grid construction settings
   !> @param[out] error Invalid setting
   subroutine check_cartesian_settings(self, error)
      !> Grid construction settings
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Invalid setting
      type(error_type), allocatable, intent(out) :: error
      !> Overflow-safe running point count
      integer(int64) :: count

      if (min(self%nx, self%ny, self%nz) < 1) then
         call fatal_error(error, "cartesian domain: nx, ny and nz must be positive")
         return
      end if
      count = int(self%nx, int64)*int(self%ny, int64)
      if (count > int(huge(0), int64)/int(self%nz, int64)) then
         call fatal_error(error, "cartesian domain: point count exceeds the integer range")
         return
      end if
      if (.not. all(ieee_is_finite([self%dr, self%xi0_factor]))) then
         call fatal_error(error, "cartesian domain: dr and xi0_factor must be finite")
         return
      end if
      if (self%dr <= 0.0_wp) then
         call fatal_error(error, "cartesian domain: dr must be finite and positive")
      else if (self%xi0_factor <= 0.0_wp) then
         call fatal_error(error, "cartesian domain: xi0_factor must be finite and positive")
      else if (.not. ieee_is_finite(self%margin) .or. self%margin < 0.0_wp) then
         call fatal_error(error, "cartesian domain: margin must be finite and nonnegative")
      end if
   end subroutine check_cartesian_settings

   !> Require every nucleus to stay `margin` away from the box faces
   !>
   !> - The box spans `center +- n*dr/2` per axis and repeats periodically
   !> - `reach(d)` is the largest nuclear distance from the centroid along axis d
   !>
   !> @param[in]  self   Grid settings
   !> @param[in]  reach  Largest nuclear distance from the centroid per axis, bohr
   !> @param[out] error  Names the first axis without room
   subroutine check_cartesian_fit(self, reach, error)
      !> Grid settings
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Largest nuclear distance from the centroid per axis
      real(wp), intent(in) :: reach(3)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=1), parameter :: axis(3) = ["x", "y", "z"]
      character(len=24) :: needed
      real(wp) :: half(3)
      integer :: d

      half = 0.5_wp*self%dr*real([self%nx, self%ny, self%nz], wp)
      do d = 1, 3
         if (reach(d) + self%margin > half(d)) then
            write (needed, "(f0.3)") 2.0_wp*(reach(d) + self%margin)
            call fatal_error(error, "cartesian domain: solute and margin do not fit the box along "// &
               & axis(d)//"; n"//axis(d)//"*dr must be at least "//trim(needed)//" bohr")
            return
         end if
      end do
   end subroutine check_cartesian_fit

   !> Commit a successfully constructed Cartesian discretization
   !>
   !> @param[in,out] source Candidate grid
   !> @param[in,out] dest Destination retaining settings
   subroutine move_cartesian_geometry(source, dest)
      !> Candidate grid
      type(moist_math_grid_3d_cartesian_type), intent(inout) :: source
      !> Destination grid
      class(moist_math_grid_3d_cartesian_type), intent(inout) :: dest
      dest%ngrid = source%ngrid
      dest%npts_k = source%npts_k
      dest%nx = source%nx
      dest%ny = source%ny
      dest%nz = source%nz
      dest%nkx = source%nkx
      dest%dr = source%dr
      dest%dv = source%dv
      dest%vbox = source%vbox
      dest%dkx = source%dkx
      dest%dky = source%dky
      dest%dkz = source%dkz
      dest%origin = source%origin
      call move_alloc(source%xyz, dest%xyz)
      call move_alloc(source%w, dest%w)
      call move_alloc(source%owner, dest%owner)
      call move_alloc(source%xi0, dest%xi0)
      call move_alloc(source%kx, dest%kx)
      call move_alloc(source%ky, dest%ky)
      call move_alloc(source%kz, dest%kz)
   end subroutine move_cartesian_geometry

   !> Center the fixed box on the current solute centroid
   !>
   !> Every nucleus must stay `margin` away from the box faces
   !>
   !> Failure keeps the committed box and center but zeroes `natom`; the grid
   !> stays unusable (contractions refuse) until an update succeeds
   !>
   !> @param[in,out] self Domain instance
   !> @param[in] mol Solute structure
   !> @param[out] error Invalid geometry or allocation failure
   subroutine update_cartesian(self, mol, error)
      !> Domain instance
      class(moist_math_grid_3d_cartesian_type), intent(inout) :: self
      !> Solute structure
      type(structure_type), intent(in) :: mol
      !> Invalid geometry or allocation failure
      type(error_type), allocatable, intent(out) :: error
      !> Requested box center
      real(wp) :: center(3)
      !> Largest nuclear distance from the requested center per axis
      real(wp) :: reach(3)
      !> Atom index
      integer :: iat

      call self%check_update_input(mol, error)
      if (allocated(error)) return
      center = sum(mol%xyz, dim=2)/real(mol%nat, wp)
      reach = 0.0_wp
      do iat = 1, mol%nat
         reach = max(reach, abs(mol%xyz(:, iat) - center))
      end do
      call check_cartesian_fit(self, reach, error)
      if (allocated(error)) return
      call realize_cartesian(self, center, error)
      if (allocated(error)) return
      self%center = center
      self%natom = mol%nat
   end subroutine update_cartesian

   !> Construct and commit a Cartesian grid from its validated settings
   !>
   !> @param[in,out] self Grid retaining atom count
   !> @param[in] center Box center, bohr
   !> @param[out] error Invalid settings or allocation failure
   subroutine realize_cartesian(self, center, error)
      !> Grid retaining atom count
      class(moist_math_grid_3d_cartesian_type), intent(inout) :: self
      !> Box center, bohr
      real(wp), intent(in) :: center(3)
      !> Invalid settings or allocation failure
      type(error_type), allocatable, intent(out) :: error
      !> Candidate grid
      type(moist_math_grid_3d_cartesian_type), allocatable :: grid
      !> Box origin
      real(wp) :: origin(3)
      !> Allocation status
      integer :: stat

      allocate (grid, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "cartesian domain: cannot allocate grid")
         return
      end if
      origin = center - 0.5_wp*self%dr*real([self%nx, self%ny, self%nz], wp)
      grid%gaussian = self%gaussian
      grid%xi0_factor = self%xi0_factor
      call build_cartesian_grid(grid, self%nx, self%ny, self%nz, self%dr, origin, error)
      if (allocated(error)) return
      call move_cartesian_geometry(grid, self)
   end subroutine realize_cartesian

   !> Domain name for diagnostics
   !>
   !> @param[in] self Domain instance
   function cartesian_kind_name(self) result(name)
      !> Domain instance
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Diagnostic name
      character(len=:), allocatable :: name
      name = "cartesian"
   end function cartesian_kind_name

   !> R2C reciprocal-space extent along x from the real-space extent
   !>
   !> @param[in] nx  real-space size along x
   pure function nk_x_from_nx(nx) result(nkx)
      !> Real-space size along x
      integer, intent(in) :: nx
      !> R2C-halved reciprocal-space size along x
      integer :: nkx
      nkx = nx/2 + 1
   end function nk_x_from_nx

   !> Initialise Cartesian grid storage from concrete arguments
   !>
   !> @param[in,out] grid Configured candidate receiving geometry
   !> @param[in]  nx     Real-space size along x
   !> @param[in]  ny     Real-space size along y
   !> @param[in]  nz     Real-space size along z
   !> @param[in]  dr     Real-space spacing, bohr
   !> @param[in]  origin Box origin, bohr
   !> @param[out] error  Invalid settings or allocation failure
   subroutine build_cartesian_grid(grid, nx, ny, nz, dr, origin, error)
      !> Grid to construct
      type(moist_math_grid_3d_cartesian_type), intent(inout) :: grid
      !> Real-space x size
      integer, intent(in) :: nx
      !> Real-space y size
      integer, intent(in) :: ny
      !> Real-space z size
      integer, intent(in) :: nz
      !> Real-space spacing, bohr
      real(wp), intent(in) :: dr
      !> Box origin, bohr
      real(wp), intent(in) :: origin(3)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i, ix, iy, iz, stat
      real(wp), parameter :: two_pi = 8.0_wp*atan(1.0_wp)

      grid%nx = nx
      grid%ny = ny
      grid%nz = nz
      grid%dr = dr
      call grid%validate(error)
      if (allocated(error)) return
      if (.not. all(ieee_is_finite(origin))) then
         call fatal_error(error, "cartesian domain: origin must be finite")
         return
      end if
      grid%ngrid = nx*ny*nz
      grid%nkx = nk_x_from_nx(nx)
      grid%npts_k = grid%nkx*ny*nz
      grid%dr = dr
      grid%dv = dr*dr*dr
      grid%vbox = real(grid%ngrid, wp)*grid%dv

      grid%dkx = two_pi/(real(nx, wp)*dr)
      grid%dky = two_pi/(real(ny, wp)*dr)
      grid%dkz = two_pi/(real(nz, wp)*dr)

      grid%origin = origin

      allocate (grid%xyz(3, grid%ngrid), stat=stat)
      if (stat == 0) allocate (grid%w(grid%ngrid), stat=stat)
      if (stat == 0) allocate (grid%owner(grid%ngrid), stat=stat)
      if (stat == 0 .and. grid%gaussian) allocate (grid%xi0(grid%ngrid), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "cartesian domain: cannot allocate geometry")
         return
      end if
      i = 0
      do iz = 1, nz
         do iy = 1, ny
            do ix = 1, nx
               i = i + 1
               grid%xyz(:, i) = origin + (real([ix, iy, iz], wp) - 0.5_wp)*dr
            end do
         end do
      end do
      grid%w = grid%dv
      grid%owner = 0
      if (grid%gaussian) grid%xi0 = grid%xi0_factor/dr

      allocate (grid%kx(grid%nkx), stat=stat)
      if (stat == 0) allocate (grid%ky(ny), stat=stat)
      if (stat == 0) allocate (grid%kz(nz), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate Cartesian reciprocal axes")
         return
      end if
      do i = 1, grid%nkx
         grid%kx(i) = real(i - 1, wp)*grid%dkx
      end do

      do i = 1, ny
         if (i - 1 <= ny/2) then
            grid%ky(i) = real(i - 1, wp)*grid%dky
         else
            grid%ky(i) = -real(ny - (i - 1), wp)*grid%dky
         end if
      end do
      do i = 1, nz
         if (i - 1 <= nz/2) then
            grid%kz(i) = real(i - 1, wp)*grid%dkz
         else
            grid%kz(i) = -real(nz - (i - 1), wp)*grid%dkz
         end if
      end do

   end subroutine build_cartesian_grid

   !> Reciprocal-space coordinate of flat k-point j (R2C layout, kx fastest)
   !>
   !> j-1 = (ikx-1) + nkx*((iky-1) + ny*(ikz-1))
   !>
   !> @param[in]  self  Grid instance
   !> @param[in]  j     Flat reciprocal-space index (1..npts_k)
   pure function cartesian_grid_3d_kpoint(self, j) result(k)
      !> Grid instance
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Flat reciprocal-space index
      integer, intent(in) :: j
      !> Coordinate of k-point j
      real(wp) :: k(3)

      integer :: j0, ikx, iky, ikz

      ! Out-of-range or absent reciprocal grid: no silent wraparound or divide by zero
      if (j < 1 .or. j > self%npts_k) error stop "cartesian domain: kpoint index outside the reciprocal grid"
      j0 = j - 1
      ikx = mod(j0, self%nkx) + 1
      iky = mod(j0/self%nkx, self%ny) + 1
      ikz = j0/(self%nkx*self%ny) + 1
      k = [self%kx(ikx), self%ky(iky), self%kz(ikz)]
   end function cartesian_grid_3d_kpoint

   !> Allocate a `moist_math_grid_3d_cartesian_trafo_type` bound to `self`
   !>
   !> @param[in]  self   Grid instance (must have the target attribute)
   !> @param[out] trafo  Allocated, bound transform engine
   !> @param[out] error  Set on allocation failure
   subroutine cartesian_grid_3d_new_trafo(self, trafo, error)
      !> Grid instance
      class(moist_math_grid_3d_cartesian_type), intent(in), target :: self
      !> Allocated transform engine
      class(moist_math_grid_3d_trafo_type), allocatable, intent(out) :: trafo
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_cartesian_trafo_type), allocatable :: t
      integer :: stat

      if (self%ngrid < 1) then
         call fatal_error(error, "cartesian domain: initialize or update before new_trafo")
         return
      end if
      allocate (t, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate Cartesian 3D trafo")
         return
      end if
      t%grid => self
      call move_alloc(t, trafo)
   end subroutine cartesian_grid_3d_new_trafo

   !> Release Cartesian-only storage, then the common geometry
   !>
   !> @param[in,out] self  Grid instance
   subroutine cartesian_grid_3d_dealloc(self)
      !> Grid instance
      class(moist_math_grid_3d_cartesian_type), intent(inout) :: self

      if (allocated(self%kx)) deallocate (self%kx)
      if (allocated(self%ky)) deallocate (self%ky)
      if (allocated(self%kz)) deallocate (self%kz)
      call self%clear_geometry()
      self%nkx = 0
      self%dv = 0.0_wp
      self%vbox = 0.0_wp
      self%dkx = 0.0_wp
      self%dky = 0.0_wp
      self%dkz = 0.0_wp
      self%origin = 0.0_wp
   end subroutine cartesian_grid_3d_dealloc

   !> Volume integral of a field already tabulated on the grid points
   !>
   !> Every Cartesian cell carries the same volume element `dV = dr**3`, so
   !> the weighted sum collapses to `dV * sum(f)`; compensated summation
   !> preserves small field contributions on fine grids
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Per-point field values (length ngrid); any other length is an error
   !> @param[out] result  Quadrature result
   !> @param[out] error   Error handling
   subroutine cartesian_grid_3d_integrate_field(self, f, result, error)
      !> Grid instance
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Per-point field values (length ngrid)
      real(wp), intent(in) :: f(:)
      !> Quadrature result
      real(wp), intent(out) :: result
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: correction, term, subtotal
      integer :: i

      result = 0.0_wp
      if (size(f) /= self%ngrid) then
         call fatal_error(error, "cartesian domain: integrate_field needs one value per grid point")
         return
      end if
      correction = 0.0_wp
      do i = 1, self%ngrid
         term = f(i) - correction
         subtotal = result + term
         correction = (subtotal - result) - term
         result = subtotal
      end do
      result = self%dv*result
   end subroutine cartesian_grid_3d_integrate_field

   !> Real-space coordinate of grid point (1, 1, 1)
   !>
   !> @param[in] self Grid instance
   pure function cartesian_get_origin(self) result(origin)
      !> Grid instance
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Box origin (bohr)
      real(wp) :: origin(3)
      origin = self%origin
   end function cartesian_get_origin

   !> Real-space cell volume, zero before a successful update
   !>
   !> @param[in] self Grid instance
   pure function cartesian_get_cell_volume(self) result(dv)
      !> Grid instance
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> `dr**3` (bohr**3)
      real(wp) :: dv
      dv = self%dv
   end function cartesian_get_cell_volume

   !> Periodic box volume, zero before a successful update
   !>
   !> @param[in] self Grid instance
   pure function cartesian_get_box_volume(self) result(vbox)
      !> Grid instance
      class(moist_math_grid_3d_cartesian_type), intent(in) :: self
      !> Box volume (bohr**3)
      real(wp) :: vbox
      vbox = self%vbox
   end function cartesian_get_box_volume

   !> Detach the trafo from its grid, idempotent
   !>
   !> @param[in,out] self  Trafo instance
   subroutine cartesian_trafo_destroy(self)
      !> Trafo instance
      class(moist_math_grid_3d_cartesian_trafo_type), intent(inout) :: self

      self%grid => null()
   end subroutine cartesian_trafo_destroy

   !> Forward 3D FFT of all `nv` sites, real-space to reciprocal-space
   !>
   !> - Real-space block `f_r(ngrid, nv)` -> reciprocal-space block
   !>   `f_k(npts_k, nv)`
   !> - Batched transform distributes spatial lines and sites across one
   !>   backend thread pool, using the OpenMP thread budget
   !> - Input preserved despite its `intent(inout)`; callers may treat it as
   !>   read-only on return
   !> - Each result multiplied by `dV = dr**3` to match the continuous
   !>   Fourier convention used elsewhere in moist
   !> - Nonzero backend status reported through `error`
   !>
   !> @param[in,out] self  Trafo instance
   !> @param[in,out] f_r   Real-space field block, shape (ngrid, nv); preserved on return
   !> @param[out]    f_k   Reciprocal-space field block, shape (npts_k, nv)
   !> @param[out]    error Error handling
   subroutine cartesian_trafo_fft_r2k(self, f_r, f_k, error)
      !> Trafo instance
      class(moist_math_grid_3d_cartesian_trafo_type), intent(inout) :: self
      !> Real-space field block, shape (ngrid, nv); preserved on return
      real(wp), intent(inout), contiguous, target :: f_r(:, :)
      !> Reciprocal-space field block, shape (npts_k, nv)
      complex(wp), intent(out), contiguous, target :: f_k(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Use one backend team across all FFT lines, including small site counts
      integer :: nthreads, iv

      call cartesian_trafo_check_blocks(self, shape(f_r), shape(f_k), error)
      if (allocated(error)) return
      nthreads = self%grid%team_size()
      ! Rank-three calls avoid batched-view overhead when only one worker is available
      if (nthreads == 1) then
         do iv = 1, size(f_r, 2)
            call cartesian_fft_r2k_one(self%grid, f_r(:, iv), f_k(:, iv), error)
            if (allocated(error)) return
         end do
         return
      end if
      if (moist_fft_r2c_3d_batch(int(self%grid%nz, c_int), int(self%grid%ny, c_int), &
         & int(self%grid%nx, c_int), int(size(f_r, 2), c_int), f_r, f_k, &
         & self%grid%dv, int(nthreads, c_int)) /= 0) then
         call fatal_error(error, "cartesian trafo: forward FFT backend failed")
      end if
   end subroutine cartesian_trafo_fft_r2k

   !> Inverse 3D FFT of all `nv` sites, reciprocal-space to real-space
   !>
   !> - Reciprocal-space block `f_k(npts_k, nv)` -> real-space block
   !>   `f_r(ngrid, nv)`
   !> - Backend uses `f_k` as its transform scratch, so `f_k` is destroyed
   !> - Spatial lines and sites share one backend thread pool
   !> - Each output normalised by `1/Vbox` to match the continuous Fourier
   !>   convention
   !> - Nonzero backend status reported through `error`
   !>
   !> @param[in,out] self  Trafo instance
   !> @param[in,out] f_k   Reciprocal-space field block, shape (npts_k, nv); destroyed on return
   !> @param[out]    f_r   Real-space field block, shape (ngrid, nv)
   !> @param[out]    error Error handling
   subroutine cartesian_trafo_fft_k2r(self, f_k, f_r, error)
      !> Trafo instance
      class(moist_math_grid_3d_cartesian_trafo_type), intent(inout) :: self
      !> Reciprocal-space field block, shape (npts_k, nv); destroyed on return
      complex(wp), intent(inout), contiguous, target :: f_k(:, :)
      !> Real-space field block, shape (ngrid, nv)
      real(wp), intent(out), contiguous, target :: f_r(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Use one backend team across all FFT lines, including small site counts
      integer :: nthreads, iv

      call cartesian_trafo_check_blocks(self, shape(f_r), shape(f_k), error)
      if (allocated(error)) return
      nthreads = self%grid%team_size()
      ! Rank-three calls avoid batched-view overhead when only one worker is available
      if (nthreads == 1) then
         do iv = 1, size(f_k, 2)
            call cartesian_fft_k2r_one(self%grid, f_k(:, iv), f_r(:, iv), error)
            if (allocated(error)) return
         end do
         return
      end if
      if (moist_fft_c2r_3d_batch(int(self%grid%nz, c_int), int(self%grid%ny, c_int), &
         & int(self%grid%nx, c_int), int(size(f_k, 2), c_int), f_k, f_r, &
         & 1.0_wp/self%grid%vbox, int(nthreads, c_int)) /= 0) then
         call fatal_error(error, "cartesian trafo: backward FFT backend failed")
      end if
   end subroutine cartesian_trafo_fft_k2r

   !> Guard - trafo bound to a grid and block shapes match it
   !>
   !> @param[in]  self     Trafo instance
   !> @param[in]  shape_r  Shape of the caller's real-space block
   !> @param[in]  shape_k  Shape of the caller's reciprocal-space block
   !> @param[out] error    Set when the trafo cannot transform these blocks
   subroutine cartesian_trafo_check_blocks(self, shape_r, shape_k, error)
      !> Trafo instance
      class(moist_math_grid_3d_cartesian_trafo_type), intent(in) :: self
      !> Shape of the real-space block
      integer, intent(in) :: shape_r(2)
      !> Shape of the reciprocal-space block
      integer, intent(in) :: shape_k(2)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. associated(self%grid)) then
         call fatal_error(error, "cartesian trafo: not bound to a grid")
         return
      end if
      call check_trafo_blocks(self%grid%ngrid, self%grid%npts_k, shape_r, shape_k, error)
   end subroutine cartesian_trafo_check_blocks

   !> Single-column forward R2C transform
   !>
   !> - Backend takes logical extents in the order (n0, n1, n2) with n0 the
   !>   slowest dimension; flat memory is Fortran column-major with x
   !>   fastest, so extents are passed as (nz, ny, nx)
   !> - `dV = dr**3` factor matching moist's continuous Fourier convention
   !>   folded into the transform rather than applied in a second pass
   !>
   !> @param[in]     grid  Grid supplying the transform geometry
   !> @param[in,out] f_r   One real-space field (ngrid); preserved on return
   !> @param[out]    f_k   One reciprocal-space field (npts_k)
   !> @param[out]    error A nonzero backend status
   subroutine cartesian_fft_r2k_one(grid, f_r, f_k, error)
      !> Grid supplying the transform geometry
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
      !> Real-space field column
      real(wp), intent(inout), contiguous, target :: f_r(:)
      !> Reciprocal-space field column
      complex(wp), intent(out), contiguous, target :: f_k(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (moist_fft_r2c_3d(int(grid%nz, c_int), int(grid%ny, c_int), &
         & int(grid%nx, c_int), f_r, f_k, grid%dv) /= 0) then
         call fatal_error(error, "cartesian trafo: forward FFT backend failed")
      end if
   end subroutine cartesian_fft_r2k_one

   !> Single-column inverse C2R transform, normalised by `1/Vbox`
   !>
   !> The backend transforms the slow axes in place in `f_k`, so `f_k` is
   !> destroyed (no per-call scratch allocation)
   !>
   !> @param[in]     grid  Grid supplying the transform geometry
   !> @param[in,out] f_k   One reciprocal-space field (npts_k); destroyed
   !> @param[out]    f_r   One real-space field (ngrid)
   !> @param[out]    error A nonzero backend status
   subroutine cartesian_fft_k2r_one(grid, f_k, f_r, error)
      !> Grid supplying the transform geometry
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
      !> Reciprocal-space field column
      complex(wp), intent(inout), contiguous, target :: f_k(:)
      !> Real-space field column
      real(wp), intent(out), contiguous, target :: f_r(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (moist_fft_c2r_3d(int(grid%nz, c_int), int(grid%ny, c_int), &
         & int(grid%nx, c_int), f_k, f_r, 1.0_wp/grid%vbox) /= 0) then
         call fatal_error(error, "cartesian trafo: backward FFT backend failed")
      end if
   end subroutine cartesian_fft_k2r_one

end module moist_math_grid_3d_cartesian
