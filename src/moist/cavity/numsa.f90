!> Numerical solvent accessible surface area (NUMSA) integrator
!>
!> Implementation based on:
!> - Original code: [github.com/grimme-lab/numsa](https://github.com/grimme-lab/numsa)
!> - Theory: Im, W., Lee, M. S., Brooks, C. L. (2003)
!>   "Generalized Born Model with a Simple Smoothing Function."
!>   *J. Comput. Chem.*, **24**(14), 1691-1702
!>   [DOI: 10.1002/jcc.10321](https://doi.org/10.1002/jcc.10321)
!>
!> Grid components (optional grid model, e.g. symmetry): node u of a
!> component of sphere a sits at c_a + r_a R u, R the component's rotation,
!> and counts with the component's weight; the rotations follow the atoms,
!> so every area derivative (`dsdr`, `area_grad`) carries the frame term
!> through the per-component torque (`frame_torque`) and the grid Jacobians
!> (`get_frame`)
module moist_cavity_numsa
   use moist_output_citations, only: citation
   use moist_model_parameters, only: moist_model_parameters_type
   use mctc_env, only: wp
   use mctc_env, only: error_type, fatal_error, get_argument, wp
   use mctc_io_convert, only: aatoau
   use moist_math_linalg_geometry, only: cross_product
   use moist_cavity_type, only: cavity_type, update_cavity_base
   use moist_cavity_grid_model_type, only: grid_model_type
   use moist_cavity_grid_frame, only: grid_frame_type
   use moist_context, only: moist_context_type
   use moist_radius_type, only: radius_type
   use mctc_io, only: structure_type
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_type, new_lebedev_grid

   implicit none(type, external)
   private

   public :: cavity_type_numsa
   public :: new_cavity_numsa

   public :: moist_cavity_numsa_parameters_type

   !> NUMSA construction parameters
   type, extends(moist_model_parameters_type) :: moist_cavity_numsa_parameters_type
      !> Lebedev points per sphere
      integer :: num_leb = 110
      !> Probe radius in bohr
      real(wp) :: probe = 0.0_wp
      !> Neighbor cutoff offset in bohr
      real(wp) :: offset = 2.0_wp*aatoau
      !> Smoothing width in bohr
      real(wp) :: smoothing = 0.3_wp*aatoau
      !> Surface exclusion tolerance
      real(wp) :: tolsesp = 1.0e-6_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Declare fields for JSON input, output, and printing
      procedure :: register_entries => register_parameter_entries
   end type moist_cavity_numsa_parameters_type

   !> NUMSA cavity state
   type, extends(cavity_type) :: cavity_type_numsa

      !> Construction settings: Lebedev points, probe radius, neighbor-list
      !> offset, smoothing width and surface exclusion tolerance
      type(moist_cavity_numsa_parameters_type) :: param

      !> Atomic numbers (nat)
      integer, allocatable :: at(:)
      !> Number of unique atom pairs
      integer :: ntpair
      !> Pair indices for neighbor list construction (2, ntpair)
      integer, allocatable :: ppind(:, :)
      !> Lebedev angular grid points (3, num_leb)
      real(wp), allocatable :: ang_grid(:, :)
      !> Lebedev solid-angle weights (num_leb), summing to 4*pi
      real(wp), allocatable :: ang_weight(:)
      !> Torque of a sphere's area on the rotation of each of its grid
      !> components (3, ncomp)
      !>
      !> - tau_c = sum_p (p - c_i) x dA_i/dp over the nodes p of component c
      !>   of sphere i, at fixed sphere centers; the frame term of `dsdr` is
      !>   its contraction with the grid Jacobians ([[grid_frame_type]])
      !> - not for the plain grid
      real(wp), allocatable :: frame_torque(:, :)
      !> Area of every grid component per unit component weight (ncomp), with it
      real(wp), allocatable :: comp_area(:)
      !> Neighbor-list cutoff radius (bohr)
      real(wp) :: srcut
      !> Number of neighbors for each atom (nat)
      integer, allocatable :: nnsas(:)
      !> Neighbor lists (nat, nat)
      integer, allocatable :: nnlists(:, :)
      !> SASA sphere radii: $r_{\mathrm{vdW}} + r_{\mathrm{probe}}$ (nat)
      real(wp), allocatable :: vdwsa(:)
      !> Precomputed radial integral weights (nat)
      real(wp), allocatable :: wrp(:)
      !> Squared smoothing boundaries: $(R_i \pm w)^2$ (2, nat)
      real(wp), allocatable :: trj2(:, :)
      !> Switching function coefficients: $a_0, a_1, a_3$
      real(wp) :: ah0, ah1, ah3

      !> Total area gradient w.r.t. atom positions (3, nat)
      real(wp), allocatable :: area_grad(:, :)
      !> Atomic area gradients (unused, kept for compatibility)
      real(wp), allocatable :: atomic_area_grad(:, :, :)
      !> Atomic surface gradient w.r.t. atom positions (3, nat, nat)
      !>
      !> - dsdr(:, j, i) = dA_i/dR_j, with the frame term in the molecular
      !>   orientation
      real(wp), allocatable :: dsdr(:, :, :)

   contains
      procedure :: update => update_cavity_numsa
      procedure :: get_gradient => compute_area_gradient_numsa
      !> Construction settings
      procedure :: parameters => numsa_parameters
   end type cavity_type_numsa

contains

   !> Restore compiled parameter defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_parameter_defaults(self)
      class(moist_cavity_numsa_parameters_type), intent(inout) :: self
      type(moist_cavity_numsa_parameters_type) :: defaults

      self%num_leb = defaults%num_leb
      self%probe = defaults%probe
      self%offset = defaults%offset
      self%smoothing = defaults%smoothing
      self%tolsesp = defaults%tolsesp
   end subroutine init_parameter_defaults

   !> Declare parameter fields for JSON input, output, and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_parameter_entries(self)
      class(moist_cavity_numsa_parameters_type), intent(inout), target :: self

      call self%register_int_scalar("num_leb", self%num_leb, level=2)
      call self%register_real_scalar("probe", self%probe, level=2)
      call self%register_real_scalar("offset", self%offset)
      call self%register_real_scalar("smoothing", self%smoothing)
      call self%register_real_scalar("tolsesp", self%tolsesp)
   end subroutine register_parameter_entries

   !> Construct from parameter values; omission uses compiled defaults
   !>
   !> @param[in,out] self Object to initialize
   !> @param[in] radii Atomic radius model to copy
   !> @param[out] error Construction error
   !> @param[in] param Configuration copied by value
   !> @param[in] ctx Borrowed run context; omitted, a model supplies its own
   !> @param[in] grid_model Grid model to copy, e.g. a symmetry model; omitted,
   !>                       every sphere keeps its plain grid
   subroutine new_cavity_numsa(self, radii, error, param, ctx, grid_model)
      !> Cavity to initialize
      type(cavity_type_numsa), intent(inout) :: self
      !> Radius model to copy
      class(radius_type), intent(in) :: radii
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Configuration; omitted means compiled defaults
      type(moist_cavity_numsa_parameters_type), intent(in), optional :: param
      !> Borrowed run context; omitted, a model supplies its own
      type(moist_context_type), intent(in), target, optional :: ctx
      !> Grid model to copy; omitted means the plain grid
      class(grid_model_type), intent(in), optional :: grid_model
      !> Resolved configuration
      type(moist_cavity_numsa_parameters_type) :: settings

      if (present(param)) settings = param
      nullify (self%ctx)
      if (present(ctx)) self%ctx => ctx
      self%label = "NUMSA"
      self%param = settings
      if (allocated(self%radius_model)) deallocate(self%radius_model)
      allocate(self%radius_model, source=radii)
      call self%set_grid_model(grid_model)
      ! A constructed cavity has no structure, no surface and no coupling yet
      if (allocated(self%mol)) deallocate (self%mol)
      self%updated = .false.
      call self%clear_surface()
      call self%clear_couplings()

   end subroutine new_cavity_numsa

   !> Construction settings of the cavity
   !>
   !> @param[in] self NUMSA cavity instance
   function numsa_parameters(self) result(param)
      !> NUMSA cavity instance
      class(cavity_type_numsa), intent(in), target :: self
      !> Settings, valid while the cavity is
      class(moist_model_parameters_type), pointer :: param

      param => self%param

   end function numsa_parameters

   !> Update cavity surface and gradients for current molecular geometry
   !>
   !> - sphere a integrates over the nodes of its grid components (`grid`);
   !>   the plain grid integrates bitwise over c + r u
   !> - molecular orientation: `dsdr` and `area_grad` carry the frame term
   !>   ([[add_frame_jacobians]])
   !>
   !> @param[in,out] self  The cavity object to update
   !> @param[in]     mol   Molecular structure with coordinates
   !> @param[out]    error Error handling
   subroutine update_cavity_numsa(self, mol, error)
      !> NumSA cavity
      class(cavity_type_numsa), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: nat
      real(wp), allocatable :: surface(:)
      real(wp), allocatable :: derivs(:, :, :)

      call update_cavity_base(self, mol, error)
      if (allocated(error)) return
      call self%require_context(error)
      if (allocated(error)) return
      call self%ctx%cite(citation%numsa)

      nat = mol%nat
      self%nsph = nat
      if (allocated(self%sphxyz)) deallocate (self%sphxyz)
      allocate (self%sphxyz(3, self%nsph), source=mol%xyz)

      call self%radius_model%update(mol, error)
      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(error)) return
      allocate (self%radii(size(self%radius_model%f0)))
      self%radii = self%radius_model%f0
      call self%update_grid_model(error)
      if (allocated(error)) return

      if (.not. allocated(self%asph)) allocate (self%asph(nat))

      ! initialize the internal numsa state and neighbour list
      call init_numsa(self, mol%num(mol%id), self%radii, self%param%probe, self%param%num_leb, &
                      self%param%offset, self%param%smoothing, error)
      if (allocated(error)) return
      call update_nnlist(self, mol%xyz)

      ! Grid components of every sphere: the plain nodes and no torque, or
      ! those of the grid model; NumSA has no Gaussian widths
      call self%build_grid(self%ang_grid, self%ang_weight, 0.0_wp, error)
      if (allocated(error)) return
      if (self%plain_grid()) then
         if (allocated(self%frame_torque)) deallocate (self%frame_torque)
         if (allocated(self%comp_area)) deallocate (self%comp_area)
      end if

      allocate (surface(nat))
      allocate (derivs(3, nat, nat))
      if (allocated(self%area_grad)) deallocate (self%area_grad)
      allocate (self%area_grad(3, nat))

      ! Surface, its derivatives with the frame term, and the total-area
      ! gradient (the sum of dsdr over the spheres)
      call compute_numsa(self, mol%xyz, surface, derivs, self%area_grad)

      call move_alloc(derivs, self%dsdr)
      self%asph = surface

      if (.not. allocated(self%total_area)) allocate (self%total_area)
      if (.not. allocated(self%total_volume)) allocate (self%total_volume)
      self%total_area = sum(self%asph)
      ! numsa integrator does not provide cavity volume; set to zero for now
      self%total_volume = 0.0_wp

      deallocate (surface)
      self%updated = .true.
      self%built = .true.

   end subroutine update_cavity_numsa

!> Retrieve gradient of total surface area
!>
!> Ensures the area gradient was computed during the last
!> [[update_cavity_numsa]] call; the gradient is stored in
!> `self%area_grad(3, nat)` and represents:
!>
!> $$
!> \frac{\partial A_{\mathrm{tot}}}{\partial {\bf R}_j} =
!> \sum_i \frac{\partial A_i}{\partial {\bf R}_j}
!> $$
!>
!> where $A_{\mathrm{tot}} = \sum_i A_i$ is the total surface area; in the
!> molecular orientation it includes the turn of the sphere grids with the
!> atoms (frame term)
!>
!> @param[in,out] self  The cavity object with cached gradients
!> @param[out]    error Error handling
   subroutine compute_area_gradient_numsa(self, error)
      !> NumSA cavity
      class(cavity_type_numsa), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%require_built(error)
      if (allocated(error)) return
      ! area_grad already prepared during update; ensure allocated
      if (.not. allocated(self%area_grad)) then
         call fatal_error(error, "Area gradient requested before cavity update in numsa cavity.")
      end if
   end subroutine compute_area_gradient_numsa

   !> Initialize internal NUMSA integrator state
   !>
   !> Prepares all geometry-independent quantities needed for surface integration:
   !>
   !> - **Switching function coefficients**: Computes $a_0, a_1, a_3$ for the cubic
   !>   polynomial $s(u) = a_0 + (a_1 + a_3 u^2) u$ that smoothly transitions from
   !>   0 to 1 over the interval $u \in [-w, w]$
   !>
   !> - **Radial integral weights**: Precomputes the analytical integral
   !>   $$
   !>   w_r = \int_{R-w}^{R+w} \left(\frac{1}{4w} + 3a_3 f(r,R,w)\right) r^3 \, dr
   !>   $$
   !>   to eliminate the need for radial quadrature
   !>
   !> - **Angular grid**: Sets up Lebedev quadrature points and weights on the unit sphere
   !>
   !> - **Neighbor pair indices**: Prepares all unique (i,j) pairs for neighbor list construction
   !>
   !> @param[inout] self      The cavity object to initialize
   !> @param[in]    num       Atomic numbers (nat)
   !> @param[in]    rad       Atomic radii in bohr (nat)
   !> @param[in]    probe     Probe radius in bohr
   !> @param[in]    nang      Number of Lebedev angular grid points
   !> @param[in]    offset    Optional cutoff radius offset (bohr)
   !> @param[in]    smoothing Optional smoothing width $w$ (bohr)
   !> @param[out]   error     Error handling
   subroutine init_numsa(self, num, rad, probe, nang, offset, smoothing, error)

      !> The cavity object to initialize
      class(cavity_type_numsa), intent(inout) :: self

      !> Atomic numbers of length nat
      integer, intent(in) :: num(:)

      !> Radii of length nat
      real(wp), intent(in) :: rad(:)

      !> Probe radius
      real(wp), intent(in) :: probe

      !> Number of Lebedev points to use
      integer, intent(in) :: nang

      !> Offset to add to the cutoff radius
      real(wp), intent(in), optional :: offset

      !> Smoothing width for soft cavity
      real(wp), intent(in), optional :: smoothing

      !> Error handling
      type(error_type), intent(out), allocatable :: error

      integer :: iat, jat, ij, izp
      !> Angular quadrature filling the cached grid, solid-angle weights summing to 4*pi
      type(moist_math_grid_angular_type) :: leb
      real(wp) :: ws, rr

      ! Set number of atoms
      self%nsph = size(num)
      if (allocated(self%at)) deallocate (self%at)
      allocate (self%at(self%nsph))
      self%at = num

      ! Allocate pair indices for all unique (i,j) combinations
      self%ntpair = self%nsph*(self%nsph - 1)/2
      if (allocated(self%ppind)) deallocate (self%ppind)
      allocate (self%ppind(2, self%ntpair))
      if (allocated(self%nnsas)) deallocate (self%nnsas)
      allocate (self%nnsas(self%nsph))
      if (allocated(self%nnlists)) deallocate (self%nnlists)
      allocate (self%nnlists(self%nsph, self%nsph))

      ! Build list of unique atom pairs for neighbor detection
      ij = 0
      do iat = 1, self%nsph
         do jat = 1, iat - 1
            ij = ij + 1
            self%ppind(1, ij) = iat
            self%ppind(2, ij) = jat
         end do
      end do

      if (allocated(self%vdwsa)) deallocate (self%vdwsa)
      if (allocated(self%trj2)) deallocate (self%trj2)
      if (allocated(self%wrp)) deallocate (self%wrp)
      allocate (self%vdwsa(self%nsph))
      allocate (self%trj2(2, self%nsph))
      allocate (self%wrp(self%nsph))

      ! Set smoothing width parameter
      if (present(smoothing)) then
         ws = smoothing
      else
         ws = 0.3_wp*aatoau
      end if

      ! Switching function polynomial coefficients
      !   s(u) = a_0 + (a_1 + a_3*u^2)*u
      !   cubic, rising smoothly from 0 to 1 over u in [-w, w]
      !   following Im, Lee, Brooks (2003) eq. (11)
      self%ah0 = 0.5_wp
      self%ah1 = 3._wp/(4.0_wp*ws)
      self%ah3 = -1._wp/(4.0_wp*(ws*(ws*ws)))

      ! Compute atom-specific quantities
      do iat = 1, self%nsph
         izp = num(iat)
         ! SASA sphere radius: vdW + probe
         self%vdwsa(iat) = rad(iat) + probe
         ! Squared boundaries of smoothing region: (R-w)^2 and (R+w)^2
         self%trj2(1, iat) = (self%vdwsa(iat) - ws)**2
         self%trj2(2, iat) = (self%vdwsa(iat) + ws)**2
         ! Precomputed radial integral weight (analytical primitive), so no
         ! radial quadrature is needed
         rr = self%vdwsa(iat) + ws
         self%wrp(iat) = (0.25_wp/ws + &
            & 3.0_wp*self%ah3*(0.2_wp*rr*rr - 0.5_wp*rr*self%vdwsa(iat) + &
            &  self%vdwsa(iat)*self%vdwsa(iat)/3.0_wp))*rr*rr*rr
         rr = self%vdwsa(iat) - ws
         self%wrp(iat) = self%wrp(iat) - (0.25_wp/ws + &
            & 3.0_wp*self%ah3*(0.2_wp*rr*rr - 0.5_wp*rr*self%vdwsa(iat) + &
            &  self%vdwsa(iat)*self%vdwsa(iat)/3.0_wp))*rr*rr*rr
      end do

      ! Neighbor-list cutoff: large enough to capture all overlapping smoothing regions
      self%srcut = 2*(ws + maxval(self%vdwsa))
      if (present(offset)) then
         self%srcut = self%srcut + offset
      else
         self%srcut = self%srcut + 2.0_wp*aatoau
      end if

      ! Set up the angular quadrature grid through the S2 grid type and cache
      ! its nodes and weights in the plain arrays the integrators consume
      call new_lebedev_grid(leb, error, npts=nang)
      if (allocated(error)) return

      call move_alloc(leb%points, self%ang_grid)
      call move_alloc(leb%weights, self%ang_weight)

   end subroutine init_numsa

   !> Update neighbor list for current molecular geometry
   !>
   !> Builds an efficient neighbor list containing only atom pairs within
   !> the cutoff distance `self%srcut`. This avoids checking all atom pairs
   !> during surface integration
   !>
   !> For each atom $i$, stores indices of neighbors $j$ such that
   !> $|{\bf R}_i - {\bf R}_j| < r_{\mathrm{cut}}$
   !>
   !> @param[inout] self The cavity object with neighbor list storage
   !> @param[in]    xyz  Atomic coordinates (3, nat) in bohr
   subroutine update_nnlist(self, xyz)
      class(cavity_type_numsa), intent(inout) :: self
      real(wp), intent(in) :: xyz(:, :)

      !> Loop indices and atom pair indices
      integer :: kk, i1, i2
      !> Squared cutoff radius
      real(wp) :: srcut2
      !> Cartesian components of interatomic vector
      real(wp) :: x, y, z
      !> Squared interatomic distance
      real(wp) :: dr2
      !> Temporary neighbor count for atom i
      integer :: nntmp_i
      !> Temporary neighbor counts (nat)
      integer, allocatable :: nntmp(:)
      !> Temporary neighbor lists (nat, nat)
      integer, allocatable :: nnls(:, :)

      srcut2 = self%srcut*self%srcut
      allocate (nnls(self%nsph, self%nsph))
      allocate (nntmp(self%nsph))
      nntmp = 0
      nnls = 0
      self%nnsas = 0
      self%nnlists = 0

      ! Check all unique pairs and add to neighbor lists if within cutoff
      do kk = 1, self%ntpair
         i1 = self%ppind(1, kk)
         i2 = self%ppind(2, kk)
         x = xyz(1, i1) - xyz(1, i2)
         y = xyz(2, i1) - xyz(2, i2)
         z = xyz(3, i1) - xyz(3, i2)
         dr2 = x*x + y*y + z*z
         if (dr2 < srcut2) then
            nntmp(i1) = nntmp(i1) + 1
            nntmp(i2) = nntmp(i2) + 1
            nnls(nntmp(i1), i1) = i2
            nnls(nntmp(i2), i2) = i1
         end if
      end do

      ! Copy temporary neighbor lists to persistent storage
      do i1 = 1, self%nsph
         nntmp_i = nntmp(i1)
         if (nntmp_i > 0) then
            do i2 = 1, nntmp_i
               self%nnlists(self%nnsas(i1) + i2, i1) = nnls(i2, i1)
            end do
            self%nnsas(i1) = self%nnsas(i1) + nntmp_i
         end if
      end do

      deallocate (nnls, nntmp)

   end subroutine update_nnlist

   !> Compute NUMSA surface areas and gradients via Lebedev quadrature
   !>
   !> Core integration routine; for each atom $i$:
   !>
   !> 1. Places Lebedev grid points on a sphere of radius
   !>    $R_i = r_{\mathrm{vdW},i} + r_{\mathrm{probe}}$
   !> 2. At each point ${\bf x}_p$, computes the accessibility weight
   !>    $s_p = \prod_j H_j({\bf x}_p)$ where $H_j$ is the switching function
   !>    describing exclusion by neighbor $j$
   !> 3. Integrates: $A_i = \sum_p w_p w_r s_p$ where $w_p$ is the Lebedev weight
   !>    and $w_r$ is the precomputed radial weight
   !> 4. Accumulates gradients $\partial A_i / \partial {\bf R}_j$ using the chain rule
   !>
   !> The switching function product ensures that buried points (inside neighbors)
   !> contribute zero, while exposed points contribute their full weight
   !>
   !> Sphere i sums over the nodes of its grid components (`grid`), each
   !> node counting with its component's weight; unless the grid is plain,
   !> `frame_torque(:, c)` takes the sum over the nodes of component c of
   !> (x_p - c_i) x dA_i/dx_p (c_i the center); `dsdrt` gains the frame term
   !> of sphere i ([[add_frame_jacobians]]) while its row is in cache
   !>
   !> @param[in,out] self      The cavity object with grid and parameters
   !> @param[in]     xyz       Atomic coordinates (3, nat) in bohr
   !> @param[out]    surface   Atomic surface areas (nat) in bohr^2
   !> @param[out]    dsdrt     Surface gradients (3, nat, nat) in bohr^2/bohr
   !> @param[out]    area_grad Total-area gradient (3, nat), the sum of dsdrt over the spheres
   subroutine compute_numsa(self, xyz, surface, dsdrt, area_grad)
      !> NumSA cavity
      class(cavity_type_numsa), intent(inout) :: self
      !> Atomic coordinates; contiguous, passed on per node without a copy
      real(wp), intent(in), contiguous :: xyz(:, :)
      !> Atomic surface areas
      real(wp), intent(out) :: surface(:)
      !> Surface gradients: dsdrt(:, j, i) = dA_i/dR_j
      real(wp), intent(out), contiguous :: dsdrt(:, :, :)
      !> Total-area gradient
      real(wp), intent(out), contiguous :: area_grad(:, :)

      !> Current atom index
      integer :: iat
      !> Lebedev grid point index
      integer :: ip
      !> Neighbor loop index and neighbor atom index
      integer :: jj, nnj
      !> Number of neighbors affecting current point
      integer :: nni
      !> Total number of neighbors for current atom
      integer :: nno
      !> SASA sphere radius for current atom
      real(wp) :: rsas
      !> Accumulated surface area for current atom
      real(wp) :: sasai
      !> Precomputed radial weight for current atom
      real(wp) :: wr
      !> Combined weight (angular * radial * accessibility)
      real(wp) :: wsa
      !> Accessibility weight at current grid point
      real(wp) :: sasap
      !> Coordinates of current atom center
      real(wp) :: xyza(3)
      !> Coordinates of current surface point
      real(wp) :: xyzp(3)
      !> Gradient contribution from neighbor jj
      real(wp) :: drjj(3)
      !> Point gradients w.r.t. neighbors (3, maxneigh)
      real(wp), allocatable :: grds(:, :)
      !> Atom gradients w.r.t. all atoms (3, nat)
      real(wp), allocatable :: grads(:, :)
      !> Indices of affecting neighbors (maxneigh)
      integer, allocatable :: grdi(:)
      !> Whether the spheres follow a molecular frame
      logical :: frame
      !> Grid component, the offset of its nodes and its weighted radial weight
      integer :: icomp, n0, nnode
      real(wp) :: wrc
      !> Node offset from the atom center, molecular orientation only
      real(wp) :: dnode(3)
      !> Derivative of a node's contribution in its position
      real(wp) :: dpoint(3)
      !> Grid Jacobians, molecular orientation only
      type(grid_frame_type) :: jacobians

      ! Every sphere assigns its whole row of dsdrt below
      surface(:) = 0.0_wp
      area_grad(:, :) = 0.0_wp

      frame = .not. self%grid%plain
      if (frame) then
         call self%get_frame(jacobians)
         if (allocated(self%frame_torque)) deallocate (self%frame_torque)
         allocate (self%frame_torque(3, self%grid%ncomp), source=0.0_wp)
         if (allocated(self%comp_area)) deallocate (self%comp_area)
         allocate (self%comp_area(self%grid%ncomp), source=0.0_wp)
      end if

      allocate (grads(3, self%nsph))
      grads = 0.0_wp
      allocate (grds(3, maxval(self%nnsas)))
      allocate (grdi(maxval(self%nnsas)))

      ! Loop over all atoms
      do iat = 1, self%nsph
         rsas = self%vdwsa(iat)
         nno = self%nnsas(iat)
         grads = 0.0_wp
         sasai = 0.0_wp
         xyza(:) = xyz(:, iat)
         wr = self%wrp(iat)

         do icomp = self%grid%sph_first(iat), self%grid%sph_first(iat + 1) - 1
            n0 = self%grid%set_first(self%grid%comp_set(icomp)) - 1
            ! The component weight scales every node of the component (1.0
            ! for a single component)
            wrc = wr*self%grid%comp_weight(icomp)

            nnode = self%grid%comp_first(icomp + 1) - self%grid%comp_first(icomp)

            ! Lebedev quadrature over the angular grid of the component; the
            ! associate names keep the node loop off the container's index
            ! arithmetic (measurably faster)
            associate (node => self%grid%node(:, n0 + 1:n0 + nnode), wang => self%grid%wang(n0 + 1:n0 + nnode))
            do ip = 1, nnode
               ! Place grid point on SASA sphere, in the component's orientation
               if (frame) then
                  dnode(:) = rsas*node(:, ip)
                  xyzp(:) = xyza(:) + dnode(:)
               else
                  xyzp(:) = xyza(:) + rsas*node(:, ip)
               end if
               ! Compute accessibility weight and gradients at this point
               call compute_w_sp(self, self%nsph, self%nnlists(:nno, iat), self%trj2, &
                                 self%vdwsa, xyz, nno, xyzp, self%ah0, self%ah1, self%ah3, sasap, &
                                 grds, nni, grdi)

               ! Accumulate surface contribution if point is accessible
               if (sasap > self%param%tolsesp) then
                  wsa = wang(ip)*wrc*sasap
                  sasai = sasai + wsa
                  ! Accumulate gradient contributions
                  ! grds contains d(sasap)/d(R_j) for each affecting neighbor
                  do jj = 1, nni
                     nnj = grdi(jj)
                     drjj(:) = wsa*grds(:, jj)
                     grads(:, iat) = grads(:, iat) + drjj(:)
                     grads(:, nnj) = grads(:, nnj) - drjj(:)
                  end do
                  ! The node also turns with its component: delta x_p = omega x dnode
                  if (frame) then
                     dpoint(:) = 0.0_wp
                     do jj = 1, nni
                        dpoint(:) = dpoint(:) + grds(:, jj)
                     end do
                     dpoint(:) = wsa*dpoint(:)
                     self%frame_torque(:, icomp) = self%frame_torque(:, icomp) + cross_product(dnode, dpoint)
                     self%comp_area(icomp) = self%comp_area(icomp) + wang(ip)*wr*sasap
                  end if
               end if
            end do
            end associate
         end do

         surface(iat) = sasai
         if (frame) call add_frame_jacobians(self, jacobians, iat, grads)
         dsdrt(:, :, iat) = grads
         area_grad(:, :) = area_grad(:, :) + grads(:, :)
      end do

      deallocate (grads, grds, grdi)

   end subroutine compute_numsa

   !> Frame term of the area derivatives: component rotations that follow the atoms
   !>
   !> The nodes of component c of sphere i turn with the angular velocity
   !> omega of its rotation, which reaches every atom j ([[grid_frame_type]]),
   !> so dA_i/dR_j gains tau_c . d omega_c/dR_j with tau_c =
   !> `frame_torque(:, c)`; only sphere i's own nodes enter A_i
   !>
   !> - called per sphere from [[compute_numsa]], on the row it just summed
   !>
   !> @param[in]     self  Cavity with the torques and areas of sphere i's components
   !> @param[in]     frame Grid Jacobians
   !> @param[in]     isph  Sphere i
   !> @param[in,out] dsdr  Surface gradients of sphere i (3, nat), dsdr(:, j) = dA_i/dR_j
   subroutine add_frame_jacobians(self, frame, isph, dsdr)
      !> NumSA cavity
      class(cavity_type_numsa), intent(in) :: self
      !> Grid Jacobians
      type(grid_frame_type), intent(in) :: frame
      !> Sphere
      integer, intent(in) :: isph
      !> Surface gradients
      real(wp), intent(inout), contiguous :: dsdr(:, :)

      integer :: icomp

      do icomp = self%grid%sph_first(isph), self%grid%sph_first(isph + 1) - 1
         call frame%add_scalar(icomp, self%frame_torque(:, icomp), dsdr)
         ! A_i is the weighted sum of its components' areas
         call frame%add_weight(icomp, self%comp_area(icomp), dsdr)
      end do

   end subroutine add_frame_jacobians

   !> Compute switching weight at a single surface point
   !>
   !> For a given point ${\bf x}_p$ on atom $i$'s SASA sphere, computes the
   !> accessibility weight as a product of switching functions from all neighbors:
   !>
   !> $$
   !> s_p = \prod_{j \in \text{neighbors}} H_j({\bf x}_p)
   !> $$
   !>
   !> where $H_j$ is the smooth switching polynomial that depends on the distance
   !> from ${\bf x}_p$ to neighbor $j$'s surface:
   !>
   !> $$
   !> u_j = |{\bf x}_p - {\bf R}_j| - R_j, \quad
   !> H_j(u_j) = \begin{cases}
   !> 0 & u_j < -w
   !> \\ a_0 + (a_1 + a_3 u_j^2) u_j & -w \le u_j \le w
   !> \\ 1 & u_j > w
   !> \end{cases}
   !> $$
   !>
   !> Simultaneously computes gradients using the chain rule:
   !>
   !> $$
   !> \frac{\partial s_p}{\partial {\bf R}_j} = s_p \frac{H_j'(u_j)}{H_j(u_j)}
   !> \frac{{\bf x}_p - {\bf R}_j}{|{\bf x}_p - {\bf R}_j|}
   !> $$
   !>
   !> Product rule applied to the logarithmic derivative
   !>
   !> @param[in]  self     The cavity object with parameters
   !> @param[in]  nat      Total number of atoms
   !> @param[in]  nnlists  Neighbor indices for this atom (nno)
   !> @param[in]  trj2     Squared smoothing boundaries (2, nat)
   !> @param[in]  vdwsa    SASA sphere radii (nat)
   !> @param[in]  xyza     All atomic coordinates (3, nat)
   !> @param[in]  nno      Number of neighbors to check
   !> @param[in]  xyzp     Surface point coordinates (3)
   !> @param[in]  ah0      Switching polynomial coefficient $a_0$
   !> @param[in]  ah1      Switching polynomial coefficient $a_1$
   !> @param[in]  ah3      Switching polynomial coefficient $a_3$
   !> @param[out] sasap    Accessibility weight $s_p$
   !> @param[out] grds     Gradient contributions (3, nno)
   !> @param[out] nni      Number of neighbors affecting this point
   !> @param[out] grdi     Indices of affecting neighbors (nno)
   pure subroutine compute_w_sp(self, nat, nnlists, trj2, vdwsa, xyza, &
                                nno, xyzp, ah0, ah1, ah3, sasap, grds, nni, grdi)
      class(cavity_type_numsa), intent(in) :: self
      integer, intent(in) :: nat
      integer, intent(in) :: nno
      integer, intent(in) :: nnlists(nno)
      integer, intent(out) :: nni
      real(wp), intent(in) :: xyza(3, nat)
      real(wp), intent(in) :: xyzp(3)
      real(wp), intent(in) :: ah0, ah1, ah3
      real(wp), intent(out) :: sasap
      real(wp), intent(out) :: grds(3, nno)
      integer, intent(out) :: grdi(nno)
      real(wp), intent(in) :: trj2(2, nat)
      real(wp), intent(in) :: vdwsa(nat)

      !> Neighbor loop index
      integer :: i
      !> Current neighbor atom index
      integer :: ia
      !> Vector from surface point to neighbor center
      real(wp) :: tj(3)
      !> Squared distance from surface point to neighbor center
      real(wp) :: tj2
      !> Distance from surface point to neighbor center
      real(wp) :: sqtj
      !> Distance from surface point to neighbor's SASA sphere
      real(wp) :: uj
      !> Intermediate value: $a_3 u^2$
      real(wp) :: ah3uj2
      !> Derivative of switching function: $s'(u)$
      real(wp) :: dsasaij
      !> Switching function value: $s(u)$
      real(wp) :: sasaij

      nni = 0
      sasap = 1.0_wp  ! Start with fully accessible

      ! Loop over neighbors and accumulate switching weights
      do i = 1, nno
         ia = nnlists(i)
         tj(:) = xyzp(:) - xyza(:, ia)
         tj2 = dot_product(tj, tj)
         ! Check if point is within outer smoothing boundary
         if (tj2 < trj2(2, ia)) then
            ! Check if point is fully buried (inside inner boundary)
            if (tj2 <= trj2(1, ia)) then
               sasap = 0.0_wp
               return  ! Completely inaccessible, skip this point
            else
               ! Point is in smoothing region: apply switching function
               sqtj = sqrt(tj2)
               uj = sqtj - vdwsa(ia)  ! Distance from neighbor's SASA sphere
               ah3uj2 = ah3*uj*uj
               dsasaij = ah1 + 3.0_wp*ah3uj2  ! Derivative s'(u)
               sasaij = ah0 + (ah1 + ah3uj2)*uj  ! Switching value s(u)
               sasap = sasap*sasaij  ! Accumulate product
               ! Compute gradient: chain rule for log derivative
               ! d(ln s)/dR = s'(u)/s(u) * du/dR where du/dR = (x_p - R_j)/|x_p - R_j|
               dsasaij = dsasaij/(sasaij*sqtj)
               nni = nni + 1
               grdi(nni) = ia
               grds(:, nni) = dsasaij*tj(:)
            end if
         end if
      end do

   end subroutine compute_w_sp

end module moist_cavity_numsa
