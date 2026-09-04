!> Gaussian-based switching functions for DROP cavities
!>
!> For efficient evaluation, a sorted per-atom neighbor list is built
!> during set_input using [[adjacency_list_type]] (O(N) cell-grid build).
!> This enables O(n_neighbors) switching function evaluation with early
!> exit instead of O(nsph) per point
module moist_cavity_drop_gaussian
   use mctc_env_accuracy, only: wp
   use mctc_io_constants, only: pi
   use mctc_io_structure, only: structure_type
   use moist_math_adjacency_list, only: adjacency_list_type
   use moist_cavity_drop_gaussian_kernel, only: iswig_pair_f0, iswig_pair_coeffs

   implicit none
   private

   public :: moist_cavity_drop_iswig, new_iswig, iswig_workspace_type

   !> Erf argument threshold beyond which erf(x) = 1 within double precision (erf(6) = 1 - 2.2e-17)
   real(wp), parameter :: erf_cutoff = 6.0_wp

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
      procedure :: swi_collect => iswig_swi_collect
      procedure :: swi1_rA_sparse => iswig_swi_f1_rA_sparse
      procedure :: swi2_rArB_sparse => iswig_swi_f2_rArB_sparse
      procedure :: swi2_rArB_block => iswig_swi_f2_rArB_block
      procedure :: swi1_rA => iswig_swi_f1_rA
   end type moist_cavity_drop_iswig

   !> Per-thread cache of the contributing neighbours at one surface point
   !>
   !> Filled once per point by [[iswig_swi_collect]]; every derivative row of
   !> that point is then a transcendental-free O(n_nb) contraction over the
   !> cache, so a single fill serves all 3*nsph directions of a Hessian pass.
   !>
   !> Allocate one per thread with `init` and reuse it: the whole point of the
   !> type is that the neighbour arrays are never allocated inside the grid
   !> loop. `init` sizes it from the adjacency list, so `reserve` only ever
   !> fires if the geometry changed underneath it.
   type :: iswig_workspace_type
      !> Neighbour slots the arrays are sized for
      integer :: capacity = 0
      !> Owner atom the cache was filled for
      integer :: owner = 0
      !> Gaussian width the cache was filled for
      real(wp) :: xi = 0.0_wp
      !> Switching function value at the cached point
      real(wp) :: f_val = 1.0_wp
      !> Contributing neighbours cached
      integer :: n_nb = 0
      !> Neighbour atom indices (capacity)
      integer, allocatable :: idx(:)
      !> Unit vectors (pos - xyz(:, k)) / r_k (3, capacity)
      real(wp), allocatable :: nhat(:, :)
      !> Reciprocal neighbour distances 1 / r_k (capacity)
      real(wp), allocatable :: rinv(:)
      !> Radial log-derivatives d ln f_k / dr and d2 ln f_k / dr2 (capacity)
      real(wp), allocatable :: c(:), e(:)
      !> Width log-derivatives d ln f_k / dxi, d2 / dxi2 and d2 / (dr dxi) (capacity)
      real(wp), allocatable :: a(:), b(:), dc(:)
   contains
      procedure :: init => iswig_workspace_init
      procedure :: reserve => iswig_workspace_reserve
      procedure :: destroy => iswig_workspace_destroy
   end type iswig_workspace_type

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
   subroutine iswig_set_input(self, mol, radii, wleb_max)
      class(moist_cavity_drop_iswig), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Atomic radii (bohr)
      real(wp), intent(in) :: radii(:)
      !> Maximum Lebedev weight (enables sorted neighbor list construction)
      real(wp), intent(in), optional :: wleb_max

      self%nsph = mol%nat

      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%radii)) deallocate (self%radii)

      allocate (self%xyz(3, self%nsph))
      allocate (self%radii(self%nsph))

      self%xyz = mol%xyz
      self%radii = radii

      ! Build sorted neighbor list when Lebedev weight info is available
      if (present(wleb_max)) then
         call iswig_build_neighbors(self, wleb_max)
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
   subroutine iswig_build_neighbors(self, wleb_max)
      class(moist_cavity_drop_iswig), intent(inout) :: self
      !> Maximum Lebedev weight (for computing minimum xi per sphere)
      real(wp), intent(in) :: wleb_max

      real(wp) :: cutoff_global

      self%R_max = maxval(self%radii)

      ! Global cutoff: conservative bound that includes all relevant pairs
      cutoff_global = self%R_max*(2.0_wp + erf_cutoff*sqrt(wleb_max)/self%swx)

      ! Build adjacency list (distances and sorting handled internally)
      call self%adj_list%init(cutoff=cutoff_global, sorted=.true.)
      call self%adj_list%update(self%xyz)

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

   !* ================================================================================= *!
   !*                          Neighbour workspace                                      *!
   !* ================================================================================= *!

   !> Size a workspace for the switching function it will be used with
   !>
   !> The adjacency list bounds the neighbour count by `maxval(nnl)`, which is
   !> O(1) in the system size; without one the fallback traversal can reach
   !> every atom. Call this once per thread, after [[iswig_set_input]]
   !>
   !> @param[out] self  Workspace to size
   !> @param[in]  iswig Switching function whose geometry it will cache
   subroutine iswig_workspace_init(self, iswig)
      !> Workspace to size
      class(iswig_workspace_type), intent(out) :: self
      !> Switching function whose geometry it will cache
      class(moist_cavity_drop_iswig), intent(in) :: iswig

      integer :: capacity

      if (iswig%adj_list%n > 0) then
         capacity = maxval(iswig%adj_list%nnl)
      else
         capacity = iswig%nsph
      end if

      call self%reserve(capacity)

   end subroutine iswig_workspace_init

   !> Ensure the workspace holds at least `n` neighbour slots
   !>
   !> @param[inout] self Workspace to grow
   !> @param[in]    n    Required neighbour slots
   pure subroutine iswig_workspace_reserve(self, n)
      !> Workspace to grow
      class(iswig_workspace_type), intent(inout) :: self
      !> Required neighbour slots
      integer, intent(in) :: n

      integer :: want

      want = max(n, 1)
      if (allocated(self%idx) .and. self%capacity >= want) return

      call self%destroy()
      self%capacity = want
      allocate (self%idx(want))
      allocate (self%nhat(3, want))
      allocate (self%rinv(want), self%c(want), self%e(want))
      allocate (self%a(want), self%b(want), self%dc(want))

   end subroutine iswig_workspace_reserve

   !> Release the neighbour arrays
   !>
   !> @param[inout] self Workspace to empty
   pure subroutine iswig_workspace_destroy(self)
      !> Workspace to empty
      class(iswig_workspace_type), intent(inout) :: self

      if (allocated(self%idx)) deallocate (self%idx)
      if (allocated(self%nhat)) deallocate (self%nhat)
      if (allocated(self%rinv)) deallocate (self%rinv)
      if (allocated(self%c)) deallocate (self%c)
      if (allocated(self%e)) deallocate (self%e)
      if (allocated(self%a)) deallocate (self%a)
      if (allocated(self%b)) deallocate (self%b)
      if (allocated(self%dc)) deallocate (self%dc)
      self%capacity = 0
      self%n_nb = 0

   end subroutine iswig_workspace_destroy

   !* ================================================================================= *!
   !*                          Neighbour collection                                     *!
   !* ================================================================================= *!

   !> Collect the contributing neighbours of one surface point
   !>
   !> The switching function is the product
   !>   f = prod [1 - 0.5 * (erf(xi*(R_j+r_ij)) + erf(xi*(R_j-r_ij)))]
   !> and this is the single place its neighbour enumeration is written down
   !>
   !> Both traversals are implemented: the sorted adjacency list with early exit, and
   !> the `adj_list%n == 0` fallback over all atoms
   !>
   !> @param[in]    self  iSwig instance
   !> @param[in]    pos   Position of surface point (3, bohr)
   !> @param[in]    owner Owner atom index
   !> @param[in]    xi    Gaussian width parameter (precomputed)
   !> @param[out]   f_val Switching function value
   !> @param[inout] work  Neighbour cache; absent for a value-only evaluation
   pure subroutine iswig_swi_collect(self, pos, owner, xi, f_val, work)
      !> iSwig instance
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Position of surface point (3, bohr)
      real(wp), intent(in) :: pos(3)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Gaussian width parameter (precomputed)
      real(wp), intent(in) :: xi
      !> Switching function value
      real(wp), intent(out) :: f_val
      !> Neighbour cache; absent for a value-only evaluation
      type(iswig_workspace_type), intent(inout), optional :: work

      logical :: have_adj, want_rows
      integer :: k, ii, start, count, n_nb
      real(wp) :: dif(3), rij, fij, break_thresh, rad_own
      real(wp) :: cc, ee, aa, bb, dcc

      f_val = 1.0_wp
      n_nb = 0
      want_rows = present(work)
      have_adj = self%adj_list%n > 0
      rad_own = self%radii(owner)

      if (have_adj) then
         start = self%adj_list%inl(owner)
         count = self%adj_list%nnl(owner)
         break_thresh = rad_own + self%R_max + erf_cutoff/xi
      else
         start = 0
         count = self%nsph
      end if

      if (want_rows) then
         call work%reserve(count)
         work%owner = owner
         work%xi = xi
      end if

      do ii = 1, count
         if (have_adj) then
            ! Sorted by d_ij ascending: all remaining neighbors are farther
            if (self%adj_list%dist(start + ii) > break_thresh) exit

            k = self%adj_list%nlat(start + ii)

            ! Centre-distance pre-screen
            if (xi*(self%adj_list%dist(start + ii) - rad_own - self%radii(k)) &
                > erf_cutoff) cycle
         else
            k = ii
            if (k == owner) cycle
         end if

         dif = pos(:) - self%xyz(:, k)
         rij = sqrt(dif(1)*dif(1) + dif(2)*dif(2) + dif(3)*dif(3))

         ! Per-atom skip: avoid erf for atoms beyond individual cutoff
         if (have_adj) then
            if (xi*(rij - self%radii(k)) > erf_cutoff) cycle
         end if

         call iswig_pair_f0(xi, self%radii(k), rij, fij)
         f_val = f_val*fij

         if (.not. want_rows) then
            if (f_val < 1.0e-14_wp) then
               f_val = 0.0_wp
               return
            end if
            cycle
         end if

         n_nb = n_nb + 1
         work%idx(n_nb) = k

         if (rij < 1.0e-30_wp .or. abs(fij) < 1.0e-30_wp) then
            ! Degenerate neighbour: it still multiplies into the product, but no
            ! derivative row can be formed for it, because every coefficient is a
            ! log-derivative dividing by `fij` and the row assembly divides by
            ! `rij`. Zero the whole slot so the contractions see exact zeros
            ! rather than an Inf that would poison them through 0 * Inf.
            work%nhat(:, n_nb) = 0.0_wp
            work%rinv(n_nb) = 0.0_wp
            work%c(n_nb) = 0.0_wp
            work%e(n_nb) = 0.0_wp
            work%a(n_nb) = 0.0_wp
            work%b(n_nb) = 0.0_wp
            work%dc(n_nb) = 0.0_wp
         else
            work%rinv(n_nb) = 1.0_wp/rij
            work%nhat(:, n_nb) = dif*work%rinv(n_nb)
            call iswig_pair_coeffs(xi, self%radii(k), rij, fij, cc, ee, aa, bb, dcc)
            work%c(n_nb) = cc
            work%e(n_nb) = ee
            work%a(n_nb) = aa
            work%b(n_nb) = bb
            work%dc(n_nb) = dcc
         end if
      end do

      if (want_rows) then
         work%n_nb = n_nb
         work%f_val = f_val
      end if

   end subroutine iswig_swi_collect

   !* ================================================================================= *!
   !*                          Switching function rows                                  *!
   !* ================================================================================= *!

   !> Compute iSwig switching function value for a single surface point.
   !>
   !> @param[in] self  iSwig instance
   !> @param[in] pos   Position of surface point (3, bohr)
   !> @param[in] owner Owner atom index
   !> @param[in] xi    Gaussian width parameter (precomputed)
   !> @return    f     Switching function value
   pure function iswig_swi_f0(self, pos, owner, xi) result(f)
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Position of surface point (3, bohr)
      real(wp), intent(in) :: pos(3)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Gaussian width parameter (precomputed)
      real(wp), intent(in) :: xi
      !> Switching function value
      real(wp) :: f

      call iswig_swi_collect(self, pos, owner, xi, f)

   end function iswig_swi_f0

   !> Sparse gradient of the iSwig switching function w.r.t. atomic positions
   !>
   !> Only the owner atom and its cached neighbours are nonzero, so the row is
   !> returned in the cache's compact index space: `rows(:, jj)` is the gradient
   !> w.r.t. atom `work%idx(jj)`, and `owner_row` the one w.r.t. the owner
   !>
   !> @param[in]  self      iSwig instance
   !> @param[in]  work      Neighbour cache filled by [[iswig_swi_collect]]
   !> @param[out] rows      d f / d r_k for the cached neighbours (3, >= n_nb)
   !> @param[out] owner_row d f / d r_owner (3)
   !> @param[out] dxi       d f / d xi at fixed geometry
   pure subroutine iswig_swi_f1_rA_sparse(self, work, rows, owner_row, dxi)
      !> iSwig instance; every input is read from `work`
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Neighbour cache filled by [[iswig_swi_collect]]
      type(iswig_workspace_type), intent(in) :: work
      !> Gradient w.r.t. the cached neighbours; only 1:work%n_nb is written
      real(wp), intent(out) :: rows(:, :)
      !> Gradient w.r.t. the owner atom
      real(wp), intent(out) :: owner_row(3)
      !> Derivative w.r.t. the Gaussian width
      real(wp), intent(out) :: dxi

      integer :: jj
      real(wp) :: coeff

      owner_row = 0.0_wp
      dxi = 0.0_wp

      do jj = 1, work%n_nb
         coeff = work%f_val*work%c(jj)
         rows(:, jj) = -coeff*work%nhat(:, jj)
         owner_row = owner_row - rows(:, jj)
         dxi = dxi + work%a(jj)
      end do

      dxi = work%f_val*dxi

   end subroutine iswig_swi_f1_rA_sparse

   !> Directional second derivative of the iSwig switching function
   !>
   !> Returns `sum_B v_B . d2 f / (d r_A d r_B)` in the same sparse layout as
   !> [[iswig_swi_f1_rA_sparse]]; supplying `vxi` changes the contraction from
   !> the position-position block to the joint `(v, v_xi)` direction, exactly as
   !> `vrad` does for the level-set `hvp_*` accessors
   !>
   !> @param[in]  self       iSwig instance
   !> @param[in]  work       Neighbour cache filled by [[iswig_swi_collect]]
   !> @param[in]  v          Nuclear displacement directions (3, nsph)
   !> @param[out] rows2      Contracted Hessian rows of the neighbours (3, >= n_nb)
   !> @param[out] owner_row2 Contracted Hessian row of the owner (3)
   !> @param[out] dxi2       Directional derivative of `dxi`
   !> @param[in]  vxi        d xi / dt along `v`; taken as zero when absent
   pure subroutine iswig_swi_f2_rArB_sparse(self, work, v, rows2, owner_row2, dxi2, vxi)
      !> iSwig instance; every geometric input is read from `work`
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Neighbour cache filled by [[iswig_swi_collect]]
      type(iswig_workspace_type), intent(in) :: work
      !> Nuclear displacement directions (3, nsph)
      real(wp), intent(in) :: v(:, :)
      !> Contracted Hessian rows of the neighbours; only 1:work%n_nb is written
      real(wp), intent(out) :: rows2(:, :)
      !> Contracted Hessian row of the owner
      real(wp), intent(out) :: owner_row2(3)
      !> Directional derivative of the width derivative
      real(wp), intent(out) :: dxi2
      !> Directional derivative of the Gaussian width; zero when absent
      real(wp), intent(in), optional :: vxi

      integer :: jj, k, owner
      real(wp) :: dr(work%n_nb)
      real(wp) :: wvec(3), mw(3), crinv, vxi_val
      real(wp) :: s_sum, a_sum, t_sum, b_sum, dl

      owner = work%owner
      vxi_val = 0.0_wp
      if (present(vxi)) vxi_val = vxi

      ! Pass 1: the four scalar reductions the rows share
      s_sum = 0.0_wp
      a_sum = 0.0_wp
      t_sum = 0.0_wp
      b_sum = 0.0_wp
      do jj = 1, work%n_nb
         k = work%idx(jj)
         wvec = v(:, owner) - v(:, k)
         dr(jj) = dot_product(work%nhat(:, jj), wvec)
         s_sum = s_sum + work%c(jj)*dr(jj)
         a_sum = a_sum + work%a(jj)
         t_sum = t_sum + work%dc(jj)*dr(jj)
         b_sum = b_sum + work%b(jj)
      end do
      dl = s_sum + a_sum*vxi_val

      ! Pass 2: the rows themselves. The owner row is the running negation, so
      ! a uniform translation - which makes every `w_k`, and hence `dl`, vanish -
      ! returns exact zeros.
      owner_row2 = 0.0_wp
      do jj = 1, work%n_nb
         k = work%idx(jj)
         wvec = v(:, owner) - v(:, k)
         crinv = work%c(jj)*work%rinv(jj)
         mw = (work%e(jj) - crinv)*dr(jj)*work%nhat(:, jj) + crinv*wvec
         rows2(:, jj) = work%f_val*(-(work%c(jj)*dl + work%dc(jj)*vxi_val) &
                                    *work%nhat(:, jj) - mw)
         owner_row2 = owner_row2 - rows2(:, jj)
      end do

      dxi2 = work%f_val*(a_sum*dl + t_sum + b_sum*vxi_val)

   end subroutine iswig_swi_f2_rArB_sparse

   !> Local second-derivative block of the iSwig switching function
   !>
   !> Returns the whole `(3, n, 3, n)` position-position block of one surface
   !> point over its influence set - the owner atom and its cached neighbours -
   !> together with the mixed width rows and the pure width curvature
   !>
   !> @param[in]  self  iSwig instance; every input is read from `work`
   !> @param[in]  work  Neighbour cache filled by [[iswig_swi_collect]]
   !> @param[out] n     Influence-set size, `work%n_nb + 1`
   !> @param[out] idx   Atom ids of the influence set, owner first
   !> @param[out] blk   Second derivative w.r.t. the influence-set positions
   !> @param[out] mix   Mixed position-width second derivative
   !> @param[out] d2xi  Second derivative w.r.t. the Gaussian width
   pure subroutine iswig_swi_f2_rArB_block(self, work, n, idx, blk, mix, d2xi)
      !> iSwig instance; every input is read from `work`
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Neighbour cache filled by [[iswig_swi_collect]]
      type(iswig_workspace_type), intent(in) :: work
      !> Influence-set size, n_nb + 1
      integer, intent(out) :: n
      !> Atom ids of the influence set; idx(1) is the owner. Only 1:n written
      integer, intent(out) :: idx(:)
      !> d2f/dR_A dR_B over the influence set. Only (:, 1:n, :, 1:n) written
      real(wp), intent(out) :: blk(:, :, :, :)
      !> d2f/dR_A dxi. Only (:, 1:n) written
      real(wp), intent(out) :: mix(:, :)
      !> d2f/dxi2
      real(wp), intent(out) :: d2xi

      integer :: ii, jj, i, j, n_nb
      real(wp) :: cn(3, work%n_nb), fcn(3, work%n_nb)
      real(wp) :: gl1(3), fgl1(3), msum(3, 3), mk(3, 3), acc(3, 3)
      real(wp) :: fval, crinv, dxi_l, b_sum, coeff

      n_nb = work%n_nb
      n = n_nb + 1
      fval = work%f_val
      idx(1) = work%owner

      ! Pass 1: the per-neighbour log-gradient vectors and the two width sums.
      ! `cn` is `c_k n_k`, the neighbour's log-gradient up to its sign, and
      ! `fcn` its `f`-scaled twin, so every rank-one entry below is a single
      ! `fcn * cn` product: `f` is applied once, and never divided out.
      gl1(:) = 0.0_wp
      dxi_l = 0.0_wp
      b_sum = 0.0_wp
      do jj = 1, n_nb
         idx(1 + jj) = work%idx(jj)
         cn(:, jj) = work%c(jj)*work%nhat(:, jj)
         fcn(:, jj) = fval*cn(:, jj)
         gl1(:) = gl1(:) + cn(:, jj)
         dxi_l = dxi_l + work%a(jj)
         b_sum = b_sum + work%b(jj)
      end do
      fgl1(:) = fval*gl1(:)

      ! Pass 2: the neighbour-neighbour blocks, upper triangle only. The
      ! diagonal carries that neighbour's own curvature `M_k`; every
      ! off-diagonal pair is pure rank one, which is why a pair of neighbours
      ! never needs a coefficient of its own.
      msum(:, :) = 0.0_wp
      do jj = 1, n_nb
         crinv = work%c(jj)*work%rinv(jj)
         do j = 1, 3
            do i = 1, j
               mk(i, j) = (work%e(jj) - crinv)*work%nhat(i, jj)*work%nhat(j, jj)
               mk(j, i) = mk(i, j)
            end do
            mk(j, j) = mk(j, j) + crinv
         end do
         msum(:, :) = msum(:, :) + mk(:, :)

         do j = 1, 3
            do i = 1, j
               blk(i, 1 + jj, j, 1 + jj) = fcn(i, jj)*cn(j, jj) + fval*mk(i, j)
               blk(j, 1 + jj, i, 1 + jj) = blk(i, 1 + jj, j, 1 + jj)
            end do
         end do

         do ii = 1, jj - 1
            do j = 1, 3
               do i = 1, 3
                  blk(i, 1 + ii, j, 1 + jj) = fcn(i, ii)*cn(j, jj)
                  blk(j, 1 + jj, i, 1 + ii) = blk(i, 1 + ii, j, 1 + jj)
               end do
            end do
         end do
      end do

      ! Pass 3: the owner strip, as the running negation of the neighbour
      ! blocks. `sum_A d2f/(dR_A dR_B) = 0` then holds to the bit for every
      ! neighbour column, and by the mirror for every neighbour row, so a
      ! uniform translation summed neighbours first contracts to exact zeros.
      do jj = 1, n_nb
         acc(:, :) = 0.0_wp
         do ii = 1, n_nb
            acc(:, :) = acc(:, :) + blk(:, 1 + ii, :, 1 + jj)
         end do
         do j = 1, 3
            do i = 1, 3
               blk(i, 1, j, 1 + jj) = -acc(i, j)
               blk(j, 1 + jj, i, 1) = -acc(i, j)
            end do
         end do
      end do

      ! Pass 4: the owner's own block, mirrored rather than negated - the
      ! header says why those two cannot both be exact in this one place.
      do j = 1, 3
         do i = 1, j
            blk(i, 1, j, 1) = fgl1(i)*gl1(j) + fval*msum(i, j)
            blk(j, 1, i, 1) = blk(i, 1, j, 1)
         end do
      end do

      ! Pass 5: the width rows, the same running negation. `dxi_L` multiplies
      ! the log-gradient, so a saturated neighbour - `a = dc = 0` bitwise -
      ! contributes nothing here either.
      mix(:, 1) = 0.0_wp
      do jj = 1, n_nb
         coeff = fval*(work%c(jj)*dxi_l + work%dc(jj))
         mix(:, 1 + jj) = -coeff*work%nhat(:, jj)
         mix(:, 1) = mix(:, 1) - mix(:, 1 + jj)
      end do

      d2xi = fval*(dxi_l*dxi_l + b_sum)

   end subroutine iswig_swi_f2_rArB_block

   !> Dense gradient of the iSwig switching function w.r.t. atomic positions
   !>
   !> @param[in]    self  iSwig instance
   !> @param[in]    pos   Position of surface point (3, bohr)
   !> @param[in]    owner Owner atom index
   !> @param[in]    xi    Gaussian width parameter (precomputed)
   !> @param[inout] work  Caller-owned neighbour cache, sized by `init`
   !> @param[out]   grad  Gradient w.r.t. atomic positions (3, nsph)
   !> @param[out]   dxi   Derivative w.r.t. the Gaussian width (optional)
   pure subroutine iswig_swi_f1_rA(self, pos, owner, xi, work, grad, dxi)
      !> iSwig instance
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Position of surface point (3, bohr)
      real(wp), intent(in) :: pos(3)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Gaussian width parameter (precomputed)
      real(wp), intent(in) :: xi
      !> Caller-owned neighbour cache
      type(iswig_workspace_type), intent(inout) :: work
      !> Gradient w.r.t. atomic positions (3, nsph)
      real(wp), intent(out) :: grad(3, self%nsph)
      !> Derivative w.r.t. the Gaussian width
      real(wp), intent(out), optional :: dxi

      ! `n_nb` never exceeds `nsph`, on either traversal, so this bound holds
      ! whatever `work%reserve` decides.
      real(wp) :: rows(3, self%nsph), owner_row(3), f_val, dxi_local
      integer :: jj

      call iswig_swi_collect(self, pos, owner, xi, f_val, work)
      call iswig_swi_f1_rA_sparse(self, work, rows, owner_row, dxi_local)

      grad = 0.0_wp
      do jj = 1, work%n_nb
         grad(:, work%idx(jj)) = rows(:, jj)
      end do
      grad(:, owner) = grad(:, owner) + owner_row

      if (present(dxi)) dxi = dxi_local

   end subroutine iswig_swi_f1_rA

end module moist_cavity_drop_gaussian
