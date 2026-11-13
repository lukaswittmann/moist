!> Gaussian-based switching functions for DROP cavities
!>
!> For efficient evaluation, a sorted per-atom neighbor list is built during
!> `set_input` using [[adjacency_list_type]], an O(N) cell-grid build
!>
!> - enables O(n_neighbors) switching function evaluation with early exit,
!>   instead of O(nsph) per point
module moist_cavity_drop_gaussian
   use mctc_env_accuracy, only: wp
   use mctc_io_constants, only: pi
   use mctc_io_structure, only: structure_type
   use moist_math_adjacency, only: adjacency_list_type

   implicit none(type, external)
   private

   public :: moist_cavity_drop_iswig, new_iswig, iswig_workspace_type

   !> Erf argument beyond which erf(x) = 1 in double precision (erf(6) = 1 - 2.2e-17)
   real(wp), parameter :: erf_cutoff = 6.0_wp

   !> Per-thread scratch of the switching gradient: one row per neighbour
   type :: iswig_workspace_type
      !> Neighbours collected by the last call
      integer :: n = 0
      !> Neighbour atoms
      integer, allocatable :: k(:)
      !> Point minus neighbour centre (3, n)
      real(wp), allocatable :: d(:, :)
      !> Distance, scaled erfc arguments and pair factor
      real(wp), allocatable :: rij(:), rplus(:), rminus(:), fij(:)
   end type iswig_workspace_type

   !> iSwig switching function type
   type :: moist_cavity_drop_iswig
      !> Gaussian width parameter (swx)
      real(wp) :: swx = 0.0_wp
      !> Number of atomic spheres
      integer :: nsph = 0
      !> Atomic positions (3, nsph)
      real(wp), allocatable :: xyz(:, :)
      !> Atomic radii (nsph)
      real(wp), allocatable :: radii(:)
      !> Maximum atomic radius (for conservative break distance)
      real(wp) :: R_max = 0.0_wp

      !> Atom-atom adjacency list (CSR, sorted by distance via cell grid in O(N))
      type(adjacency_list_type) :: adj_list
   contains
      procedure :: update => iswig_set_input
      procedure :: xi0 => iswig_xi0
      procedure :: xi1_rA => iswig_xi1_rA
      procedure :: swi0 => iswig_swi_f0
      procedure :: new_workspace => iswig_new_workspace
      procedure :: add_f1_rA => iswig_add_f1_rA
   end type moist_cavity_drop_iswig

contains

   !> Constructor for iSwig switching function
   subroutine new_iswig(self, swx)
      type(moist_cavity_drop_iswig), intent(out) :: self
      !> Gaussian width parameter
      real(wp), intent(in) :: swx

      self%swx = swx

   end subroutine new_iswig

   !> Set molecular geometry and radii for iSwig switching function
   !>
   !> @param[in] mol      Molecular structure
   !> @param[in] radii    Atomic radii (bohr)
   !> @param[in] wleb_max Maximum Lebedev quadrature weight (optional)
   !> @param[in] nthreads OpenMP team size of the neighbour-list build (optional)
   subroutine iswig_set_input(self, mol, radii, wleb_max, nthreads)
      class(moist_cavity_drop_iswig), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Atomic radii (bohr)
      real(wp), intent(in) :: radii(:)
      !> Maximum Lebedev weight (enables sorted neighbor list construction)
      real(wp), intent(in), optional :: wleb_max
      !> OpenMP team size of the neighbour-list build
      integer, intent(in), optional :: nthreads

      self%nsph = mol%nat

      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%radii)) deallocate (self%radii)

      allocate (self%xyz(3, self%nsph))
      allocate (self%radii(self%nsph))

      self%xyz = mol%xyz
      self%radii = radii

      ! Build sorted neighbor list when Lebedev weight info is available
      if (present(wleb_max)) then
         call iswig_build_neighbors(self, wleb_max, nthreads)
      else
         call self%adj_list%destroy()
      end if

   end subroutine iswig_set_input

   !> Build per-atom neighbor list for iSwiG screening
   !>
   !> Uses [[adjacency_list_type]] for O(N) construction via its internal
   !> cell grid; the adjacency list sorts each atom's neighbors by
   !> center-center distance, enabling early exit in the switching function
   !>
   !> The global cutoff is the maximum per-atom break distance of
   !>   cutoff = R_max * (2 + erf_cutoff * sqrt(wleb_max) / swx)
   !>
   !> @param[inout] self     iSwig instance with xyz and radii set
   !> @param[in]    wleb_max Maximum Lebedev quadrature weight
   !> @param[in]    nthreads OpenMP team size (optional)
   subroutine iswig_build_neighbors(self, wleb_max, nthreads)
      class(moist_cavity_drop_iswig), intent(inout) :: self
      !> Maximum Lebedev weight (for computing minimum xi per sphere)
      real(wp), intent(in) :: wleb_max
      !> OpenMP team size
      integer, intent(in), optional :: nthreads

      real(wp) :: cutoff_global

      self%R_max = maxval(self%radii)

      ! Global cutoff: conservative bound that includes all relevant pairs
      ! Per-atom break distance is R_i + R_max + erf_cutoff / xi_min(i)
      ! where xi_min(i) = swx / (R_i * sqrt(wleb_max))
      ! Maximum over all atoms occurs at R_i = R_max
      cutoff_global = self%R_max*(2.0_wp + erf_cutoff*sqrt(wleb_max)/self%swx)

      ! Build adjacency list (distances and sorting handled internally)
      call self%adj_list%init(cutoff=cutoff_global, sorted=.true.)
      call self%adj_list%update(self%xyz, nthreads)

   end subroutine iswig_build_neighbors

   !> Compute xi (Gaussian width) for a single grid point
   !>
   !> Computes xi = swx / (R_owner * sqrt(wleb))
   !>
   !> @param[in] owner Owner atom index
   !> @param[in] wleb  Lebedev weight (raw weight, not normalized)
   !> @return    xi    Gaussian width parameter
   pure function iswig_xi0(self, owner, wleb) result(xi)
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Owner atom index
      integer, intent(in) :: owner
      !> Lebedev weight (raw weight, not normalized)
      real(wp), intent(in) :: wleb
      !> Gaussian width
      real(wp) :: xi

      if (wleb == 0.0_wp) then
         xi = 0.0_wp
      else
         xi = self%swx/(self%radii(owner)*sqrt(wleb))
      end if

   end function iswig_xi0

   !> Compute derivative of xi w.r.t. atomic positions
   !>
   !> Computes dxi / dR_A where xi = swx / (R_owner * sqrt(wleb))
   !>
   !> When an optional active index list is provided, only the derivative
   !> components for those atoms are computed (all others remain zero)
   !>
   !> @param[in] owner     Owner atom index
   !> @param[in] wleb      Lebedev weight
   !> @param[in] wleb1_rA  Derivative of wleb w.r.t. atomic positions (3, nsph)
   !> @param[in] active    Optional list of active atom indices for screening
   !> @return    xi1_rA    Derivative of xi w.r.t. atomic positions (3, nsph)
   pure function iswig_xi1_rA(self, owner, wleb, wleb1_rA, active) result(xi1_rA)
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Owner atom index
      integer, intent(in) :: owner
      !> Lebedev weight
      real(wp), intent(in) :: wleb
      !> Derivative of wleb w.r.t. atomic positions (3, nsph)
      real(wp), intent(in) :: wleb1_rA(3, self%nsph)
      !> Optional list of active atom indices for screening
      integer, intent(in), optional :: active(:)
      !> Derivative of xi w.r.t. atomic positions (3, nsph)
      real(wp) :: xi1_rA(3, self%nsph)

      real(wp) :: factor
      integer :: ii, iatom

      factor = -self%swx/(2.0_wp*self%radii(owner)*sqrt(wleb)*wleb)

      if (present(active)) then
         xi1_rA = 0.0_wp
         do ii = 1, size(active)
            iatom = active(ii)
            xi1_rA(:, iatom) = factor*wleb1_rA(:, iatom)
         end do
      else
         xi1_rA(:, :) = factor*wleb1_rA(:, :)
      end if

   end function iswig_xi1_rA

   !> Compute iSwig switching function value for a single surface point
   !>
   !> The switching function is computed as:
   !>   f = prod 0.5 * (erfc(xi*(R_j+r_ij)) + erfc(xi*(R_j-r_ij)))
   !>
   !> When a sorted per-atom neighbor list is available (built in set_input),
   !> the loop iterates over sorted neighbors and exits early once the
   !> center-center distance exceeds the break threshold
   !>
   !> @param[in] pos   Position of surface point (3, bohr)
   !> @param[in] owner Owner atom index
   !> @param[in] xi    Gaussian width parameter (precomputed)
   function iswig_swi_f0(self, pos, owner, xi) result(f)
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Position of surface point (3, bohr)
      real(wp), intent(in) :: pos(3)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Gaussian width parameter (precomputed)
      real(wp), intent(in) :: xi
      !> Switching function value
      real(wp) :: f

      integer :: i, ii, start, count
      real(wp) :: rij, rplus, rminus, f_tmp, break_thresh

      f = 1.0_wp

      if (self%adj_list%n > 0) then
         ! Sorted per-atom neighbor list: iterate and exit early
         start = self%adj_list%inl(owner)
         count = self%adj_list%nnl(owner)
         break_thresh = self%radii(owner) + self%R_max + erf_cutoff/xi

         do ii = 1, count
            ! Sorted by d_ij ascending: all remaining neighbors are farther
            if (self%adj_list%dist(start + ii) > break_thresh) exit

            i = self%adj_list%nlat(start + ii)
            rij = norm2(pos(:) - self%xyz(:, i))

            ! Per-atom skip: avoid erfc for atoms beyond individual cutoff
            if (xi*(rij - self%radii(i)) > erf_cutoff) cycle

            rplus = xi*(self%radii(i) + rij)
            rminus = xi*(self%radii(i) - rij)
            f_tmp = 0.5_wp*(erfc(rplus) + erfc(rminus))
            f = f*f_tmp

            if (f < 1.0e-14_wp) then
               f = 0.0_wp
               return
            end if
         end do
      else
         ! Fallback: loop over all atoms
         do i = 1, self%nsph
            if (i == owner) cycle

            rij = norm2(pos(:) - self%xyz(:, i))

            rplus = xi*(self%radii(i) + rij)
            rminus = xi*(self%radii(i) - rij)
            f_tmp = 0.5_wp*(erfc(rplus) + erfc(rminus))
            f = f*f_tmp

            if (f < 1.0e-14_wp) then
               f = 0.0_wp
               return
            end if
         end do
      end if

   end function iswig_swi_f0

   !> Workspace of [[iswig_add_f1_rA]] for the largest neighbour list
   !>
   !> Allocate once per thread, outside the point loop
   !>
   !> @param[out] work Workspace sized for every owner
   pure subroutine iswig_new_workspace(self, work)
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Workspace sized for every owner
      type(iswig_workspace_type), intent(out) :: work

      integer :: nmax

      if (self%adj_list%n > 0) then
         nmax = maxval(self%adj_list%nnl)
      else
         nmax = self%nsph
      end if
      nmax = max(1, nmax)
      allocate (work%k(nmax), work%d(3, nmax), work%rij(nmax), work%rplus(nmax), &
                work%rminus(nmax), work%fij(nmax))
      work%n = 0
   end subroutine iswig_new_workspace

   !> Add the weighted nuclear gradient of the switching function at fixed width
   !>
   !> Adds `weight * df/dR_k`, with
   !>    df/dR_k = f * sum_j (1/f_j) * df_j/dR_k,
   !> to the owner and to every neighbour inside the erf cutoff
   !>
   !> - `xi` is held fixed: DROP evaluates `f` at the anchor, whose width has no
   !>   nuclear dependence
   !> - O(n_nb) and allocation-free; the atoms touched besides the owner are
   !>   left in `work%k(1:work%n)`
   !> - the sorted per-atom neighbour list, when built, is walked with early exit
   !>
   !> @param[in]    pos    Position of surface point (3, bohr)
   !> @param[in]    owner  Owner atom index
   !> @param[in]    xi     Gaussian width parameter (precomputed)
   !> @param[in]    weight Adjoint weight of `f`
   !> @param[inout] work   Workspace from [[iswig_new_workspace]]
   !> @param[inout] grad   Gradient accumulator (3, nsph)
   pure subroutine iswig_add_f1_rA(self, pos, owner, xi, weight, work, grad)
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Position of surface point (3, bohr)
      real(wp), intent(in) :: pos(3)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Gaussian width parameter (precomputed)
      real(wp), intent(in) :: xi
      !> Adjoint weight of `f`
      real(wp), intent(in) :: weight
      !> Workspace
      type(iswig_workspace_type), intent(inout) :: work
      !> Gradient accumulator (3, nsph)
      real(wp), intent(inout) :: grad(:, :)

      integer :: k, ii, start, count, jj
      real(wp) :: xdif, ydif, zdif, rij, rplus, rminus, fij, f_val, prefac, coeff
      real(wp) :: break_thresh, exp_rp, exp_rm
      !> Walk the sorted neighbour list
      logical :: use_list
      !> Owner row before the weight
      real(wp) :: g_owner(3)

      ! Pass 1: f_val product, neighbour intermediates cached in `work`; the
      ! sorted neighbour list exits early and skips atoms beyond their erf
      ! cutoff, the fallback walks every atom without that cutoff
      use_list = self%adj_list%n > 0
      if (use_list) then
         start = self%adj_list%inl(owner)
         count = self%adj_list%nnl(owner)
         break_thresh = self%radii(owner) + self%R_max + erf_cutoff/xi
      else
         start = 0
         count = self%nsph
         break_thresh = huge(1.0_wp)
      end if
      f_val = 1.0_wp
      work%n = 0
      do ii = 1, count
         if (use_list) then
            if (self%adj_list%dist(start + ii) > break_thresh) exit
            k = self%adj_list%nlat(start + ii)
         else
            k = ii
            if (k == owner) cycle
         end if
         xdif = pos(1) - self%xyz(1, k)
         ydif = pos(2) - self%xyz(2, k)
         zdif = pos(3) - self%xyz(3, k)
         rij = sqrt(xdif*xdif + ydif*ydif + zdif*zdif)
         if (use_list) then
            if (xi*(rij - self%radii(k)) > erf_cutoff) cycle
         end if

         rplus = xi*(self%radii(k) + rij)
         rminus = xi*(self%radii(k) - rij)
         fij = 0.5_wp*(erfc(rplus) + erfc(rminus))
         f_val = f_val*fij

         work%n = work%n + 1
         work%k(work%n) = k
         work%d(:, work%n) = [xdif, ydif, zdif]
         work%rij(work%n) = rij
         work%rplus(work%n) = rplus
         work%rminus(work%n) = rminus
         work%fij(work%n) = fij
      end do

      ! Pass 2: spatial gradient; the owner row collects the opposite of every
      ! neighbour row
      prefac = -f_val*xi/sqrt(pi)
      g_owner = 0.0_wp
      do jj = 1, work%n
         k = work%k(jj)
         rij = work%rij(jj)
         fij = work%fij(jj)
         if (rij < 1.0e-30_wp) cycle
         if (abs(fij) < 1.0e-30_wp) cycle

         exp_rp = exp(-work%rplus(jj)*work%rplus(jj))
         exp_rm = exp(-work%rminus(jj)*work%rminus(jj))
         coeff = prefac/(fij*rij)*(exp_rp - exp_rm)

         grad(:, k) = grad(:, k) + weight*(-coeff*work%d(:, jj))
         g_owner = g_owner + coeff*work%d(:, jj)
      end do
      grad(:, owner) = grad(:, owner) + weight*g_owner

   end subroutine iswig_add_f1_rA

end module moist_cavity_drop_gaussian
