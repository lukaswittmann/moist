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

   implicit none
   private

   public :: moist_cavity_drop_iswig, new_iswig

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
      procedure :: swi1_rA => iswig_swi_f1_rA
   end type moist_cavity_drop_iswig

contains

   !> Constructor for iSwig switching function
   subroutine new_iswig(self, swx)
      type(moist_cavity_drop_iswig), intent(out) :: self
      !> Gaussian width parameter
      real(wp), intent(in) :: swx

      self%swx = swx

   end subroutine new_iswig

   !> Set molecular geometry and radii for iSwig switching function.
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

      ! Global cutoff: conservative bound that includes all relevant pairs.
      ! Per-atom break distance is R_i + R_max + erf_cutoff / xi_min(i)
      ! where xi_min(i) = swx / (R_i * sqrt(wleb_max)).
      ! Maximum over all atoms occurs at R_i = R_max.
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

   !> Compute iSwig switching function value for a single surface point.
   !>
   !> The switching function is computed as:
   !>   f = prod [1 - 0.5 * (erf(xi*(R_j+r_ij)) + erf(xi*(R_j-r_ij)))]
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

            ! Per-atom skip: avoid erf for atoms beyond individual cutoff
            if (xi*(rij - self%radii(i)) > erf_cutoff) cycle

            rplus = xi*(self%radii(i) + rij)
            rminus = xi*(self%radii(i) - rij)
            f_tmp = 1.0_wp - 0.5_wp*(erf(rplus) + erf(rminus))
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
            f_tmp = 1.0_wp - 0.5_wp*(erf(rplus) + erf(rminus))
            f = f*f_tmp

            if (f < 1.0e-14_wp) then
               f = 0.0_wp
               return
            end if
         end do
      end if

   end function iswig_swi_f0

   !> Compute gradient of iSwig switching function w.r.t. atomic positions
   !>
   !> The gradient is computed as:
   !>    df/dr_k = f * sum_j (1/f_j) * df_j/dr_k
   !>
   !> When a sorted per-atom neighbor list is available, all three inner
   !> loops use it with early exit
   !>
   !> @param[in]  pos    Position of surface point (3, bohr)
   !> @param[in]  owner  Owner atom index
   !> @param[in]  xi     Gaussian width parameter (precomputed)
   !> @param[in]  xi1_rA Derivative of xi w.r.t. atomic positions (3, nsph)
   !> @param[in]  active Optional list of active atom indices for output screening
   pure function iswig_swi_f1_rA(self, pos, owner, xi, xi1_rA, active) result(grad)
      class(moist_cavity_drop_iswig), intent(in) :: self
      !> Position of surface point (3, bohr)
      real(wp), intent(in) :: pos(3)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Gaussian width parameter (precomputed)
      real(wp), intent(in) :: xi
      !> Derivative of xi w.r.t. atomic positions (3, nsph)
      real(wp), intent(in) :: xi1_rA(3, self%nsph)
      !> Optional list of active atom indices for output screening
      integer, intent(in), optional :: active(:)
      !> Gradient w.r.t. atomic positions (3, nsph)
      real(wp) :: grad(3, self%nsph)

      integer :: k, ii, iatom, start, count, n_nb, jj
      real(wp) :: xdif, ydif, zdif, rij, rplus, rminus, rplus2, rminus2
      real(wp) :: fij, f_val, prefac, coeff, total_coeff, break_thresh
      real(wp) :: exp_rp, exp_rm
      integer, allocatable :: nb_k(:)
      real(wp), allocatable :: nb_xdif(:), nb_ydif(:), nb_zdif(:)
      real(wp), allocatable :: nb_rij(:), nb_rplus(:), nb_rminus(:), nb_fij(:)

      grad = 0.0_wp

      if (self%adj_list%n > 0) then
         ! Sorted per-atom neighbor list: iterate with early exit
         start = self%adj_list%inl(owner)
         count = self%adj_list%nnl(owner)
         break_thresh = self%radii(owner) + self%R_max + erf_cutoff/xi

         allocate (nb_k(count))
         allocate (nb_xdif(count), nb_ydif(count), nb_zdif(count), &
                   nb_rij(count), nb_rplus(count), nb_rminus(count), nb_fij(count))

         ! Pass 1: f_val product and cache intermediates
         f_val = 1.0_wp
         n_nb = 0
         do ii = 1, count
            if (self%adj_list%dist(start + ii) > break_thresh) exit
            k = self%adj_list%nlat(start + ii)

            xdif = pos(1) - self%xyz(1, k)
            ydif = pos(2) - self%xyz(2, k)
            zdif = pos(3) - self%xyz(3, k)
            rij = sqrt(xdif*xdif + ydif*ydif + zdif*zdif)

            if (xi*(rij - self%radii(k)) > erf_cutoff) cycle

            rplus = xi*(self%radii(k) + rij)
            rminus = xi*(self%radii(k) - rij)
            fij = 1.0_wp - 0.5_wp*(erf(rplus) + erf(rminus))
            f_val = f_val*fij

            n_nb = n_nb + 1
            nb_k(n_nb) = k
            nb_xdif(n_nb) = xdif
            nb_ydif(n_nb) = ydif
            nb_zdif(n_nb) = zdif
            nb_rij(n_nb) = rij
            nb_rplus(n_nb) = rplus
            nb_rminus(n_nb) = rminus
            nb_fij(n_nb) = fij
         end do

         ! Pass 2: fused spatial gradient + xi chain-rule coefficient
         prefac = -f_val*xi/sqrt(pi)
         total_coeff = 0.0_wp

         do jj = 1, n_nb
            k = nb_k(jj)
            rij = nb_rij(jj)
            fij = nb_fij(jj)

            if (rij < 1.0e-30_wp) cycle
            if (abs(fij) < 1.0e-30_wp) cycle

            rplus2 = nb_rplus(jj)*nb_rplus(jj)
            rminus2 = nb_rminus(jj)*nb_rminus(jj)
            exp_rp = exp(-rplus2)
            exp_rm = exp(-rminus2)

            coeff = prefac/(fij*rij)*(exp_rp - exp_rm)

            grad(1, k) = grad(1, k) - coeff*nb_xdif(jj)
            grad(2, k) = grad(2, k) - coeff*nb_ydif(jj)
            grad(3, k) = grad(3, k) - coeff*nb_zdif(jj)

            grad(1, owner) = grad(1, owner) + coeff*nb_xdif(jj)
            grad(2, owner) = grad(2, owner) + coeff*nb_ydif(jj)
            grad(3, owner) = grad(3, owner) + coeff*nb_zdif(jj)

            if (self%radii(k) /= 0.0_wp) then
               total_coeff = total_coeff - f_val/(fij*sqrt(pi))* &
                             ((self%radii(k) + rij)*exp_rp + (self%radii(k) - rij)*exp_rm)
            end if
         end do
      else
         ! Fallback: full loop over all atoms

         allocate (nb_k(self%nsph))
         allocate (nb_xdif(self%nsph), nb_ydif(self%nsph), nb_zdif(self%nsph), &
                   nb_rij(self%nsph), nb_rplus(self%nsph), nb_rminus(self%nsph), &
                   nb_fij(self%nsph))

         ! Pass 1: f_val product and cache intermediates
         f_val = 1.0_wp
         n_nb = 0
         do k = 1, self%nsph
            if (k == owner) cycle
            xdif = pos(1) - self%xyz(1, k)
            ydif = pos(2) - self%xyz(2, k)
            zdif = pos(3) - self%xyz(3, k)
            rij = sqrt(xdif*xdif + ydif*ydif + zdif*zdif)
            rplus = xi*(self%radii(k) + rij)
            rminus = xi*(self%radii(k) - rij)
            fij = 1.0_wp - 0.5_wp*(erf(rplus) + erf(rminus))
            f_val = f_val*fij

            n_nb = n_nb + 1
            nb_k(n_nb) = k
            nb_xdif(n_nb) = xdif
            nb_ydif(n_nb) = ydif
            nb_zdif(n_nb) = zdif
            nb_rij(n_nb) = rij
            nb_rplus(n_nb) = rplus
            nb_rminus(n_nb) = rminus
            nb_fij(n_nb) = fij
         end do

         ! Pass 2: fused spatial gradient + xi chain-rule coefficient
         prefac = -f_val*xi/sqrt(pi)
         total_coeff = 0.0_wp

         do jj = 1, n_nb
            k = nb_k(jj)
            rij = nb_rij(jj)
            fij = nb_fij(jj)
            if (rij < 1.0e-30_wp) cycle
            if (abs(fij) < 1.0e-30_wp) cycle
            rplus2 = nb_rplus(jj)*nb_rplus(jj)
            rminus2 = nb_rminus(jj)*nb_rminus(jj)
            exp_rp = exp(-rplus2)
            exp_rm = exp(-rminus2)
            coeff = prefac/(fij*rij)*(exp_rp - exp_rm)
            grad(1, k) = grad(1, k) - coeff*nb_xdif(jj)
            grad(2, k) = grad(2, k) - coeff*nb_ydif(jj)
            grad(3, k) = grad(3, k) - coeff*nb_zdif(jj)
            grad(1, owner) = grad(1, owner) + coeff*nb_xdif(jj)
            grad(2, owner) = grad(2, owner) + coeff*nb_ydif(jj)
            grad(3, owner) = grad(3, owner) + coeff*nb_zdif(jj)
            if (self%radii(k) /= 0.0_wp) then
               total_coeff = total_coeff - f_val/(fij*sqrt(pi))* &
                             ((self%radii(k) + rij)*exp_rp + (self%radii(k) - rij)*exp_rm)
            end if
         end do
      end if

      ! Apply accumulated xi chain-rule coefficient (screened by active list)
      if (present(active)) then
         do ii = 1, size(active)
            iatom = active(ii)
            grad(:, iatom) = grad(:, iatom) + total_coeff*xi1_rA(:, iatom)
         end do
      else
         do iatom = 1, self%nsph
            grad(:, iatom) = grad(:, iatom) + total_coeff*xi1_rA(:, iatom)
         end do
      end if

   end function iswig_swi_f1_rA

end module moist_cavity_drop_gaussian
