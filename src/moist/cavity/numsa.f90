!> Numerical solvent accessible surface area (NUMSA) integrator
!>
!> This module implements the NUMSA method for computing solvent accessible
!> surface area (SASA) and its derivatives using Lebedev angular quadrature
!> with smooth switching functions for neighbor exclusion.
!>
!> ## References
!>
!> Implementation based on:
!> - Original code: [github.com/grimme-lab/numsa](https://github.com/grimme-lab/numsa)
!> - Theory: Im, W., Lee, M. S., Brooks, C. L. (2003).
!>   "Generalized Born Model with a Simple Smoothing Function."
!>   *J. Comput. Chem.*, **24**(14), 1691-1702.
!>   [DOI: 10.1002/jcc.10321](https://doi.org/10.1002/jcc.10321)
!>
!> ## Mathematical Background
!>
!> The method computes SASA via angular integration on Lebedev spheres with
!> smooth exclusion weights describing overlaps with neighboring atoms.
!>
!> For each atom $i$, surface points are placed on a sphere of radius
!> $R_i = r_{\mathrm{vdW},i} + r_{\mathrm{probe}}$ and weighted by a product
!> of switching functions from all neighbors:
!>
!> $$
!> A_i = \int_{\Omega} w_r \prod_{j \in \text{neighbors}} H_j({\bf x}) \, d\Omega
!> $$
!>
!> where $w_r$ is a precomputed radial integral and $H_j$ is a smooth switching
!> polynomial that transitions from 0 (fully buried) to 1 (fully exposed).
!>
module moist_cavity_numsa
   use mctc_env, only: wp
   use mctc_env, only: error_type, fatal_error, get_argument, wp
   use mctc_io_convert, only: aatoau
   use moist_type, only: cavity_type
   use moist_radius_type, only: radius_type
   use mctc_io, only: structure_type
   use moist_math_grid_lebedev, only: get_angular_grid, lebedev_order_from_num
   use mctc_io_constants, only: pi

   implicit none
   private

   public :: cavity_type_numsa
   public :: new_cavity_numsa

   !> NUMSA cavity integrator type
   !>
   !> Extends the base cavity_type with NUMSA-specific state for computing
   !> solvent accessible surface area and its gradients.
   type, extends(cavity_type) :: cavity_type_numsa

      !* Configuration parameters

      !> Number of Lebedev angular grid points
      integer :: num_leb = 110
      !> Probe radius for solvent sphere (bohr)
      real(wp) :: probe = 0.0_wp*aatoau
      !> Offset added to neighbor-list cutoff radius (bohr)
      real(wp) :: offset = 2.0_wp*aatoau
      !> Smoothing width $w$ for switching function (bohr)
      real(wp) :: smoothing = 0.3_wp*aatoau
      !> Tolerance for surface point exclusion
      real(wp) :: tolsesp = 1.e-6_wp

      !* Internal integrator state

      !> Number of atoms
      integer :: nat
      !> Atomic numbers (nat)
      integer, allocatable :: at(:)
      !> Number of unique atom pairs
      integer :: ntpair
      !> Pair indices for neighbor list construction (2, ntpair)
      integer, allocatable :: ppind(:, :)
      !> Lebedev angular grid points (3, num_leb)
      real(wp), allocatable :: ang_grid(:, :)
      !> Lebedev quadrature weights (num_leb)
      real(wp), allocatable :: ang_weight(:)
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

      !* Cached gradients

      !> Total area gradient w.r.t. atom positions (3, nat)
      real(wp), allocatable :: area_grad(:, :)
      !> Atomic area gradients (unused, kept for compatibility)
      real(wp), allocatable :: atomic_area_grad(:, :, :)
      !> Atomic surface gradient w.r.t. atom positions (3, nat, nat)
      real(wp), allocatable :: dsdr(:, :, :)

   contains
      procedure :: update => update_cavity_numsa
      procedure :: get_gradient => compute_area_gradient_numsa
   end type cavity_type_numsa

contains

   ! !> Legacy numsa wrapper
   ! subroutine calc_numsa(mol, total_area, total_volume, asph, error)

   !    type(structure_type), intent(in) :: mol
   !    real(wp), intent(out) :: total_area, total_volume
   !    real(wp), allocatable, intent(out) :: asph(:)
   !    type(error_type), allocatable, intent(out) :: error
   !    integer :: ngrid, num_leb, nsph
   !    real(wp) :: cut_a, cut_f
   !    real(wp), allocatable :: c(:, :), a(:), xi(:), f(:), wleb(:)
   !    integer, allocatable :: owner(:)
   !    type(cavity_type_numsa), allocatable :: cavity
   !    class(radius_type), allocatable :: radius_model

   !    nsph = mol%nat
   !    num_leb = 110
   !    cut_a = 0.0_wp
   !    cut_f = 1.0e-7_wp

   !    allocate(asph(nsph))

   !    allocate(cavity)
   !    call new_radii("cpcm", radius_model, error)
   !    if (allocated(error)) return

   !    call new_cavity_iswig(cavity, num_leb, cut_a, cut_f, radius_model=radius_model, error=error)
   !    if (allocated(error)) return

   !    call cavity%update(mol, error=error)
   !    if (allocated(error)) return

   !    total_area = cavity%total_area
   !    total_volume = cavity%total_volume
   !    asph = cavity%asph

   ! end subroutine calc_numsa

!> Constructor for NUMSA cavity integrator
!>
!> Initializes a NUMSA cavity object with optional configuration parameters.
!> The actual surface computation happens later during [[update_cavity_numsa]].
!>
!> @param[inout] self       The cavity object to initialize
!> @param[in]    nleb       Number of Lebedev angular grid points (optional)
!> @param[in]    probe_r    Probe radius in bohr (optional)
!> @param[in]    offset_r   Cutoff offset in bohr (optional)
!> @param[in]    smoothing_r Smoothing width $w$ in bohr (optional)
   subroutine new_cavity_numsa(self, nleb, probe_r, offset_r, smoothing_r, &
                               radii, error)
      type(cavity_type_numsa), intent(inout) :: self
      integer, intent(in), optional :: nleb
      real(wp), intent(in), optional :: probe_r, offset_r, smoothing_r
      class(radius_type), intent(in) :: radii
      type(error_type), allocatable, intent(out) :: error

      if (present(nleb)) self%num_leb = nleb
      if (present(probe_r)) self%probe = probe_r
      if (present(offset_r)) self%offset = offset_r
      if (present(smoothing_r)) self%smoothing = smoothing_r
      if (allocated(self%radius_model)) deallocate (self%radius_model)
      allocate (self%radius_model, source=radii)

   end subroutine new_cavity_numsa

!> Update cavity surface and gradients for current molecular geometry
!>
!> This is the main entry point for NUMSA surface computation. It:
!>
!> 1. Initializes internal state (angular grid, switching parameters)
!> 2. Builds neighbor lists for efficient overlap detection
!> 3. Computes atomic surface areas via Lebedev quadrature
!> 4. Computes surface gradients $\partial A_i / \partial {\bf R}_j$
!> 5. Caches total area gradient for later retrieval
!>
!> The surface area for each atom is stored in `self%asph`, and the
!> full gradient tensor $\partial A_i / \partial {\bf R}_j$ is stored
!> in `self%dsdr(3, j, i)`.
!>
!> @param[inout] self   The cavity object to update
!> @param[in]    mol    Molecular structure with coordinates
!> @param[in]    radii  Atomic radii (typically vdW radii) in bohr (nat)
   subroutine update_cavity_numsa(self, mol, error)
      class(cavity_type_numsa), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error

      integer :: nat
      real(wp), allocatable :: surface(:)
      real(wp), allocatable :: dsdr(:, :, :)
      integer :: iat, jatom, i

      nat = mol%nat

      call self%radius_model%update(mol, error)
      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(error)) return
      allocate (self%radii(size(self%radius_model%f0)))
      self%radii = self%radius_model%f0

      if (.not. allocated(self%asph)) allocate (self%asph(nat))

      ! initialize the internal numsa state and neighbour list
      call init_numsa(self, mol%num(mol%id), self%radii, self%probe, self%num_leb, self%offset, self%smoothing, error)
      call update_nnlist(self, mol%xyz)

      allocate (surface(nat))
      allocate (dsdr(3, nat, nat))

      ! compute surface and gradients via internal implementation
      call compute_numsa(self, mol%xyz, surface, dsdr)

      self%dsdr = dsdr
      self%asph = surface

      if (.not. allocated(self%total_area)) allocate (self%total_area)
      if (.not. allocated(self%total_volume)) allocate (self%total_volume)
      self%total_area = sum(self%asph)
      ! numsa integrator does not provide cavity volume; set to zero for now
      self%total_volume = 0.0_wp

      ! cache dsdr inside area_grad storage for later gradient computation
      if (allocated(self%area_grad)) deallocate (self%area_grad)
      allocate (self%area_grad(3, nat))
      self%area_grad = 0.0_wp
      ! accumulate total-area gradient by summing contributions of each sphere
      do jatom = 1, nat
         do iat = 1, nat
            self%area_grad(:, jatom) = self%area_grad(:, jatom) + dsdr(:, jatom, iat)
         end do
      end do

      deallocate (surface, dsdr)

   end subroutine update_cavity_numsa

!> Retrieve gradient of total surface area
!>
!> This routine ensures that the area gradient has been computed during
!> the last [[update_cavity_numsa]] call. The gradient is already stored
!> in `self%area_grad(3, nat)` and represents:
!>
!> $$
!> \frac{\partial A_{\mathrm{tot}}}{\partial {\bf R}_j} =
!> \sum_i \frac{\partial A_i}{\partial {\bf R}_j}
!> $$
!>
!> where $A_{\mathrm{tot}} = \sum_i A_i$ is the total surface area.
!>
!> @param[inout] self The cavity object with cached gradients
   subroutine compute_area_gradient_numsa(self)
      class(cavity_type_numsa), intent(inout) :: self
      type(error_type), allocatable :: error
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
!>   0 to 1 over the interval $u \in [-w, w]$.
!>
!> - **Radial integral weights**: Precomputes the analytical integral
!>   $$
!>   w_r = \int_{R-w}^{R+w} \left(\frac{1}{4w} + 3a_3 f(r,R,w)\right) r^3 \, dr
!>   $$
!>   to eliminate the need for radial quadrature.
!>
!> - **Angular grid**: Sets up Lebedev quadrature points and weights on the unit sphere.
!>
!> - **Neighbor pair indices**: Prepares all unique (i,j) pairs for neighbor list construction.
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

      integer :: iat, jat, ij, oleb, izp
      real(wp) :: ws, rr

      ! Set number of atoms
      self%nat = size(num)
      if (allocated(self%at)) deallocate (self%at)
      allocate (self%at(self%nat))
      self%at = num

      ! Set number of Lebedev points
      self%num_leb = nang

      ! Allocate pair indices for all unique (i,j) combinations
      self%ntpair = self%nat*(self%nat - 1)/2
      if (allocated(self%ppind)) deallocate (self%ppind)
      allocate (self%ppind(2, self%ntpair))
      if (allocated(self%nnsas)) deallocate (self%nnsas)
      allocate (self%nnsas(self%nat))
      if (allocated(self%nnlists)) deallocate (self%nnlists)
      allocate (self%nnlists(self%nat, self%nat))

      ! Build list of unique atom pairs for neighbor detection
      ij = 0
      do iat = 1, self%nat
         do jat = 1, iat - 1
            ij = ij + 1
            self%ppind(1, ij) = iat
            self%ppind(2, ij) = jat
         end do
      end do

      if (allocated(self%vdwsa)) deallocate (self%vdwsa)
      if (allocated(self%trj2)) deallocate (self%trj2)
      if (allocated(self%wrp)) deallocate (self%wrp)
      allocate (self%vdwsa(self%nat))
      allocate (self%trj2(2, self%nat))
      allocate (self%wrp(self%nat))

      ! Set smoothing width parameter
      if (present(smoothing)) then
         ws = smoothing
      else
         ws = 0.3_wp*aatoau
      end if

      ! Switching function polynomial coefficients for s(u) = a_0 + (a_1 + a_3*u^2)*u
      ! This cubic polynomial smoothly transitions from 0 to 1 over u in [-w, w]
      ! following Im, Lee, Brooks (2003) eq. (11)
      self%ah0 = 0.5_wp
      self%ah1 = 3._wp/(4.0_wp*ws)
      self%ah3 = -1._wp/(4.0_wp*(ws*(ws*ws)))

      ! Compute atom-specific quantities
      do iat = 1, self%nat
         izp = num(iat)
         ! SASA sphere radius: vdW + probe
         self%vdwsa(iat) = rad(iat) + probe
         ! Squared boundaries of smoothing region: (R-w)^2 and (R+w)^2
         self%trj2(1, iat) = (self%vdwsa(iat) - ws)**2
         self%trj2(2, iat) = (self%vdwsa(iat) + ws)**2
         ! Precomputed radial integral weight (analytical primitive)
         ! This eliminates the need for radial quadrature
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

      ! Set up Lebedev angular quadrature grid
      ! Map requested num_leb to Lebedev order index
      call lebedev_order_from_num(nang, oleb, error)
      if (allocated(error)) return

      if (allocated(self%ang_grid)) deallocate (self%ang_grid)
      if (allocated(self%ang_weight)) deallocate (self%ang_weight)
      allocate (self%ang_grid(3, nang))
      allocate (self%ang_weight(nang))
      call get_angular_grid(oleb, self%ang_grid, self%ang_weight, error)
      if (allocated(error)) return

      ! Scale weights for full sphere (Lebedev weights integrate to 1)
      self%ang_weight(:) = self%ang_weight*4.0_wp*pi

   end subroutine init_numsa

!> Update neighbor list for current molecular geometry
!>
!> Builds an efficient neighbor list containing only atom pairs within
!> the cutoff distance `self%srcut`. This avoids checking all atom pairs
!> during surface integration.
!>
!> For each atom $i$, stores indices of neighbors $j$ such that
!> $|{\bf R}_i - {\bf R}_j| < r_{\mathrm{cut}}$.
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
      allocate (nnls(self%nat, self%nat))
      allocate (nntmp(self%nat))
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
      do i1 = 1, self%nat
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
!> This is the core integration routine. For each atom $i$:
!>
!> 1. Places Lebedev grid points on a sphere of radius $R_i = r_{\mathrm{vdW},i} + r_{\mathrm{probe}}$
!> 2. At each point ${\bf x}_p$, computes the accessibility weight
!>    $s_p = \prod_j H_j({\bf x}_p)$ where $H_j$ is the switching function
!>    describing exclusion by neighbor $j$
!> 3. Integrates: $A_i = \sum_p w_p w_r s_p$ where $w_p$ is the Lebedev weight
!>    and $w_r$ is the precomputed radial weight
!> 4. Accumulates gradients $\partial A_i / \partial {\bf R}_j$ using the chain rule
!>
!> The switching function product ensures that buried points (inside neighbors)
!> contribute zero, while exposed points contribute their full weight.
!>
!> @param[inout] self    The cavity object with grid and parameters
!> @param[in]    xyz     Atomic coordinates (3, nat) in bohr
!> @param[out]   surface Atomic surface areas (nat) in bohr^2
!> @param[out]   dsdrt   Surface gradients (3, nat, nat) in bohr^2/bohr
   subroutine compute_numsa(self, xyz, surface, dsdrt)
      class(cavity_type_numsa), intent(inout) :: self
      real(wp), intent(in) :: xyz(:, :)
      real(wp), intent(out) :: surface(:)
      real(wp), intent(out) :: dsdrt(:, :, :)

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

      surface(:) = 0.0_wp
      dsdrt(:, :, :) = 0.0_wp

      allocate (grads(3, self%nat))
      grads = 0.0_wp
      allocate (grds(3, maxval(self%nnsas)))
      allocate (grdi(maxval(self%nnsas)))

      ! Loop over all atoms
      do iat = 1, self%nat
         rsas = self%vdwsa(iat)
         nno = self%nnsas(iat)
         grads = 0.0_wp
         sasai = 0.0_wp
         xyza(:) = xyz(:, iat)
         wr = self%wrp(iat)

         ! Lebedev quadrature over angular grid
         do ip = 1, size(self%ang_grid, 2)
            ! Place grid point on SASA sphere
            xyzp(:) = xyza(:) + rsas*self%ang_grid(:, ip)
            ! Compute accessibility weight and gradients at this point
            call compute_w_sp(self, self%nat, self%nnlists(:nno, iat), self%trj2, self%vdwsa, xyz, nno, xyzp, &
               & self%ah0, self%ah1, self%ah3, sasap, grds, nni, grdi)

            ! Accumulate surface contribution if point is accessible
            if (sasap > self%tolsesp) then
               wsa = self%ang_weight(ip)*wr*sasap
               sasai = sasai + wsa
               ! Accumulate gradient contributions
               ! grds contains d(sasap)/d(R_j) for each affecting neighbor
               do jj = 1, nni
                  nnj = grdi(jj)
                  drjj(:) = wsa*grds(:, jj)
                  grads(:, iat) = grads(:, iat) + drjj(:)
                  grads(:, nnj) = grads(:, nnj) - drjj(:)
               end do
            end if
         end do

         surface(iat) = sasai
         dsdrt(:, :, iat) = grads
      end do

      deallocate (grads, grds, grdi)

   end subroutine compute_numsa

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
!> 0 & u_j < -w \\
!> a_0 + (a_1 + a_3 u_j^2) u_j & -w \le u_j \le w \\
!> 1 & u_j > w
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
!> This is the product rule applied to the logarithmic derivative.
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
      integer, intent(in) :: nnlists(nno)
      integer, intent(in) :: nno
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
