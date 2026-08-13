!> Cartesian-monomial Gaussian basis + density evaluator for the internal
!> isodensity DROP level set function.
!>
!> The internal isodensity LSF needs the electron density rho(r) and its spatial
!> derivatives (up to third order) at arbitrary points.  Rather than reproduce a
!> host code's spherical-harmonic / normalization conventions, moist evaluates
!> only bare cartesian monomial Gaussians
!>
!>    g_c(r) = (x-Rx)^lx (y-Ry)^ly (z-Rz)^lz * sum_p coeff_p exp(-a_p |r-R|^2),
!>
!> and consumes a density matrix that has already been transformed into this
!> cartesian-monomial basis (D_cart).  The physical density
!>
!>    rho(r) = sum_{c,c'} D_cart(c,c') g_c(r) g_c'(r)
!>
!> is representation-invariant, so the host builds the (fixed) transform once
!> from its own basis metadata and passes D_cart each SCF step.  The cartesian
!> component ordering used here is reported through the API so the host can match
!> it exactly.
!>
!> A single cartesian primitive factorizes as g = X(dx) Y(dy) Z(dz), so every
!> spatial derivative is a product of the exp-stripped 1D factor derivatives
!> tabulated by [[moist_iso_gto_poly1d]] (sympy-generated).
module moist_cavity_drop_lsf_isodensity_gto
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_math_blas, only: symv, gemm
   implicit none
   private

   public :: moist_iso_gto_type
   public :: moist_iso_gto_ncart

   !> Highest supported angular momentum (s=0 .. l=8).  The parameter and the
   !> [[moist_iso_gto_poly1d]] routine below are emitted by
   !> config/gen_isodensity_gto.py into the marked regions -- DO NOT EDIT THEM BY
   !> HAND; regenerate with `python3 config/gen_isodensity_gto.py`.
   ! >>> GENERATED lmax >>>
   integer, parameter :: moist_iso_gto_lmax = 8
   ! <<< GENERATED lmax <<<

   !> Cartesian-component count above which the dense density contraction routes
   !> through BLAS dsymv rather than the matmul intrinsic.  For small solutes the
   !> per-call BLAS overhead outweighs the tuned kernel; for large ones (the
   !> compact-solute case where the dense contraction dominates) dsymv wins ~3x.
   integer, parameter :: iso_blas_ncart_min = 96

   !> Symmetric Hessian slot index (4..9) for spatial axis pair (i, j).
   integer, parameter :: hess_slot(3, 3) = reshape([ &
                                                   4, 5, 6, &
                                                   5, 7, 8, &
                                                   6, 8, 9], [3, 3])

   !> Symmetric third-derivative slot index (10..19) for axis triple (i, j, k).
   !> Built column-major as third_slot(k, j, i); symmetric under any permutation.
   integer, parameter :: third_slot(3, 3, 3) = reshape([ &
                                                       10, 11, 12, 11, 13, 14, 12, 14, 15, &
                                                       11, 13, 14, 13, 16, 17, 14, 17, 18, &
                                                       12, 14, 15, 14, 17, 18, 15, 18, 19], [3, 3, 3])

   !> Cumulative number of packed derivative slots through each order, 0..4.
   integer, parameter :: nslot_by_order(0:4) = [1, 4, 10, 20, 35]

   !> Number of x derivatives represented by each packed slot, 0..34.
   integer, parameter :: deriv_x(0:34) = [ &
                         0, 1, 0, 0, 2, 1, 1, 0, 0, 0, &
                         3, 2, 2, 1, 1, 1, 0, 0, 0, 0, &
                         4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0]

   !> Number of y derivatives represented by each packed slot, 0..34.
   integer, parameter :: deriv_y(0:34) = [ &
                         0, 0, 1, 0, 0, 1, 0, 2, 1, 0, &
                         0, 1, 0, 2, 1, 0, 3, 2, 1, 0, &
                         0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0]

   !> Number of z derivatives represented by each packed slot, 0..34.
   integer, parameter :: deriv_z(0:34) = [ &
                         0, 0, 0, 1, 0, 0, 1, 0, 1, 2, &
                         0, 0, 1, 0, 1, 2, 0, 1, 2, 3, &
                         0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4]

   !> Cartesian-monomial Gaussian basis together with a density matrix expressed
   !> in that basis.  Read-only after setup / set_density, so it is safe to share
   !> across the cavity's per-thread LSF clones.
   type :: moist_iso_gto_type
      !> Number of shells
      integer :: nshell = 0
      !> Total number of cartesian components (sum of ncart(l) over shells)
      integer :: ncart = 0
      !> Highest angular momentum present
      integer :: lmax = 0
      !> Per-shell owner atom (1-based), size nshell
      integer, allocatable :: sh_atom(:)
      !> Per-shell angular momentum, size nshell
      integer, allocatable :: sh_l(:)
      !> Per-shell primitive CSR offsets into exps/coeffs, size nshell+1
      integer, allocatable :: sh_poff(:)
      !> Per-shell cartesian-component CSR offsets, size nshell+1; the global
      !> component indices of shell s are sh_coff(s)+1 .. sh_coff(s+1)
      integer, allocatable :: sh_coff(:)
      !> Primitive Gaussian exponents, size nprim_tot
      real(wp), allocatable :: exps(:)
      !> Primitive contraction coefficients (host-normalized), size nprim_tot
      real(wp), allocatable :: coeffs(:)
      !> Current shell centers in Bohr (refreshed from the molecule), (3, nshell)
      real(wp), allocatable :: center(:, :)
      !> Cartesian monomial powers per global component, (3, ncart); this defines
      !> the ordering the host must match when building the density transform
      integer, allocatable :: comp_l(:, :)
      !> Density matrix in the cartesian-monomial basis, symmetric (ncart, ncart)
      real(wp), allocatable :: dcart(:, :)
      !> Highest owner-atom index present in the basis (= size of the CSR below)
      integer :: natom_grid = 0
      !> Atom -> shell CSR offsets, size natom_grid+1; atom A owns shells
      !> atom_shell(atom_soff(A)+1 : atom_soff(A+1)).  Lets a candidate-atom list
      !> (from the cavity cell grid) be turned into a shell list without a scan.
      integer, allocatable :: atom_soff(:)
      !> Flat-packed shell indices grouped by owner atom, size nshell
      integer, allocatable :: atom_shell(:)
      !> Per-shell squared radial cutoff: the shell contributes nothing (to value
      !> or derivatives through the order requested from [[gto_build_screening]])
      !> beyond |r - center|^2 > sh_rcut2. Set by
      !> [[gto_build_screening]]; huge (no screening) until then / when disabled.
      real(wp), allocatable :: sh_rcut2(:)
      !> Largest shell reach sqrt(max(sh_rcut2)); the radius the cavity cell grid
      !> must span so no contributing atom is missed as a candidate.
      real(wp) :: max_rcut = 0.0_wp
   contains
      !> Configure the basis from per-shell primitive data
      procedure :: init => gto_init
      !> Refresh shell centers from a molecular structure
      procedure :: refresh_centers => gto_refresh_centers
      !> Install the cartesian-monomial density matrix
      procedure :: set_density => gto_set_density
      !> Whether a density matrix has been installed
      procedure :: has_density => gto_has_density
      !> (Re)compute the per-shell radial screening cutoffs for a threshold
      procedure :: build_screening => gto_build_screening
      !> Conservative global shell reach (Bohr) for a screening threshold
      procedure :: reach => gto_reach
      !> Evaluate the density and the requested spatial derivatives at one point
      procedure :: eval => gto_eval
   end type moist_iso_gto_type

contains

   !> Number of cartesian components for angular momentum l.
   !>
   !> @param[in] l  Angular momentum
   !> @returns      (l+1)(l+2)/2
   elemental function moist_iso_gto_ncart(l) result(n)
      !> Angular momentum
      integer, intent(in) :: l
      integer :: n

      n = (l + 1)*(l + 2)/2
   end function moist_iso_gto_ncart

   !> Symmetric fourth-derivative slot index (20..34) for a spatial axis 4-tuple.
   !>
   !> The fourth-derivative table has the 15 unique cartesian components with
   !> kx+ky+kz = 4, laid out (like the lower orders) in descending-kx-then-ky
   !> order at phi columns 20..34.  This maps any axis quadruple (i, j, k, l),
   !> each in 1..3, to its column by counting the per-axis derivative orders --
   !> cheaper and less error-prone than a 3x3x3x3 lookup table, and only used on
   !> the (non-hot) fourth-derivative path.
   !>
   !> @param[in] i First spatial axis (1..3)
   !> @param[in] j Second spatial axis (1..3)
   !> @param[in] k Third spatial axis (1..3)
   !> @param[in] l Fourth spatial axis (1..3)
   !> @returns     phi column index in 20..34
   pure function fourth_slot(i, j, k, l) result(slot)
      integer, intent(in) :: i, j, k, l
      integer :: slot

      integer :: ax(4), nx, ny, kx

      ax = [i, j, k, l]
      nx = count(ax == 1)
      ny = count(ax == 2)
      !> 0-based position of (kx=nx, ky=ny) in the descending-kx-then-ky order.
      slot = 4 - nx - ny
      do kx = nx + 1, 4
         slot = slot + (5 - kx)
      end do
      slot = 20 + slot
   end function fourth_slot

   !> Configure the basis and derive the cartesian-component layout.
   !>
   !> @param[inout] self      Basis instance
   !> @param[in]    sh_atom   Per-shell owner atom index (1-based), size nshell
   !> @param[in]    sh_l      Per-shell angular momentum, size nshell
   !> @param[in]    sh_nprim  Per-shell primitive count, size nshell
   !> @param[in]    exps      Primitive exponents (concatenated per shell)
   !> @param[in]    coeffs    Primitive contraction coefficients (host-normalized)
   !> @param[out]   error     Set on invalid input
   subroutine gto_init(self, sh_atom, sh_l, sh_nprim, exps, coeffs, error)
      !> Basis instance
      class(moist_iso_gto_type), intent(inout) :: self
      !> Per-shell owner atom index (1-based)
      integer, intent(in) :: sh_atom(:)
      !> Per-shell angular momentum
      integer, intent(in) :: sh_l(:)
      !> Per-shell primitive count
      integer, intent(in) :: sh_nprim(:)
      !> Primitive exponents
      real(wp), intent(in) :: exps(:)
      !> Primitive contraction coefficients
      real(wp), intent(in) :: coeffs(:)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer :: s, l, c, lx, ly, lz, nprim_tot

      self%nshell = size(sh_l)
      if (size(sh_atom) /= self%nshell .or. size(sh_nprim) /= self%nshell) then
         call fatal_error(error, "isodensity basis shell arrays have inconsistent length")
         return
      end if
      if (self%nshell < 1) then
         call fatal_error(error, "isodensity basis has no shells")
         return
      end if
      if (any(sh_l < 0) .or. any(sh_l > moist_iso_gto_lmax)) then
         call fatal_error(error, "isodensity basis angular momentum out of supported range")
         return
      end if
      if (any(sh_nprim < 1)) then
         call fatal_error(error, "isodensity basis shell has no primitives")
         return
      end if

      nprim_tot = sum(sh_nprim)
      if (size(exps) /= nprim_tot .or. size(coeffs) /= nprim_tot) then
         call fatal_error(error, "isodensity basis primitive arrays have inconsistent length")
         return
      end if

      ! Reinitialization is supported: discard all storage derived from the
      ! previous basis only after the new input has passed validation.
      if (allocated(self%sh_poff)) deallocate (self%sh_poff)
      if (allocated(self%sh_coff)) deallocate (self%sh_coff)
      if (allocated(self%center)) deallocate (self%center)
      if (allocated(self%comp_l)) deallocate (self%comp_l)
      if (allocated(self%dcart)) deallocate (self%dcart)
      if (allocated(self%atom_soff)) deallocate (self%atom_soff)
      if (allocated(self%atom_shell)) deallocate (self%atom_shell)
      if (allocated(self%sh_rcut2)) deallocate (self%sh_rcut2)

      self%lmax = maxval(sh_l)
      self%sh_atom = sh_atom
      self%sh_l = sh_l
      self%exps = exps
      self%coeffs = coeffs

      !> Primitive CSR offsets (1-based, exclusive upper bound at s+1)
      allocate (self%sh_poff(self%nshell + 1))
      self%sh_poff(1) = 1
      do s = 1, self%nshell
         self%sh_poff(s + 1) = self%sh_poff(s) + sh_nprim(s)
      end do

      !> Cartesian-component CSR offsets
      allocate (self%sh_coff(self%nshell + 1))
      self%sh_coff(1) = 0
      do s = 1, self%nshell
         self%sh_coff(s + 1) = self%sh_coff(s) + moist_iso_gto_ncart(sh_l(s))
      end do
      self%ncart = self%sh_coff(self%nshell + 1)

      !> Canonical cartesian monomial ordering per shell: lx from l down to 0,
      !> ly from (l-lx) down to 0, lz = l-lx-ly.  Matches the common CINT order.
      allocate (self%comp_l(3, self%ncart))
      do s = 1, self%nshell
         l = sh_l(s)
         c = self%sh_coff(s)
         do lx = l, 0, -1
            do ly = l - lx, 0, -1
               lz = l - lx - ly
               c = c + 1
               self%comp_l(:, c) = [lx, ly, lz]
            end do
         end do
      end do

      allocate (self%center(3, self%nshell), source=0.0_wp)

      !> Atom -> shell CSR (counting sort on the owner atom).  Candidate atom
      !> lists from the cavity cell grid index straight into this, so the hot
      !> projection path visits only the shells on nearby atoms.
      self%natom_grid = maxval(sh_atom)
      allocate (self%atom_soff(self%natom_grid + 1), source=0)
      do s = 1, self%nshell
         self%atom_soff(sh_atom(s) + 1) = self%atom_soff(sh_atom(s) + 1) + 1
      end do
      do s = 1, self%natom_grid
         self%atom_soff(s + 1) = self%atom_soff(s + 1) + self%atom_soff(s)
      end do
      allocate (self%atom_shell(self%nshell))
      block
         integer :: fill(self%natom_grid)
         fill = 0
         do s = 1, self%nshell
            fill(sh_atom(s)) = fill(sh_atom(s)) + 1
            self%atom_shell(self%atom_soff(sh_atom(s)) + fill(sh_atom(s))) = s
         end do
      end block

      !> Default to "no screening" (huge cutoff) until build_screening runs, so a
      !> direct user (no cavity / zero threshold) evaluates every shell exactly.
      allocate (self%sh_rcut2(self%nshell), source=huge(1.0_wp))
      self%max_rcut = 0.0_wp
   end subroutine gto_init

   !> Refresh the shell centers from a molecular structure.
   !>
   !> @param[inout] self  Basis instance
   !> @param[in]    mol   Molecular structure (positions in Bohr)
   subroutine gto_refresh_centers(self, mol)
      !> Basis instance
      class(moist_iso_gto_type), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      integer :: s

      do s = 1, self%nshell
         self%center(:, s) = mol%xyz(:, self%sh_atom(s))
      end do
   end subroutine gto_refresh_centers

   !> Install the cartesian-monomial density matrix.
   !>
   !> @param[inout] self   Basis instance
   !> @param[in]    dcart  Density matrix in the cartesian-monomial basis
   !> @param[out]   error  Set on a size mismatch
   subroutine gto_set_density(self, dcart, error)
      !> Basis instance
      class(moist_iso_gto_type), intent(inout) :: self
      !> Cartesian-monomial density matrix (ncart, ncart)
      real(wp), intent(in) :: dcart(:, :)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      if (size(dcart, 1) /= self%ncart .or. size(dcart, 2) /= self%ncart) then
         call fatal_error(error, "isodensity density matrix has the wrong shape")
         return
      end if
      if (allocated(self%dcart)) deallocate (self%dcart)
      self%dcart = dcart
   end subroutine gto_set_density

   !> Whether a density matrix has been installed.
   !>
   !> @param[in] self  Basis instance
   !> @returns         .true. once set_density has been called
   pure function gto_has_density(self) result(ok)
      !> Basis instance
      class(moist_iso_gto_type), intent(in) :: self
      logical :: ok

      ok = allocated(self%dcart)
   end function gto_has_density

   !> Conservative squared amplitude bound for one primitive at radius ``r``.
   !>
   !> Bounds the magnitude of any cartesian component of ``x^lx y^ly z^lz
   !> exp(-a r^2)`` and of its spatial derivatives up to order ``k`` by
   !> ``|c| (r+1)^l (2 a r + l + 1)^k exp(-a r^2)`` -- an over-estimate (the
   !> monomial is <= r^l <= (r+1)^l and each differentiation multiplies the
   !> bound by at most 2 a r + l), so screening on it never drops a shell that
   !> still contributes above the threshold.
   !>
   !> @param[in] c  Absolute contraction coefficient
   !> @param[in] a  Primitive exponent
   !> @param[in] l  Shell angular momentum
   !> @param[in] k  Highest derivative order to cover
   !> @param[in] r  Radius (Bohr)
   !> @returns      Amplitude upper bound
   pure function prim_amp_bound(c, a, l, k, r) result(b)
      real(wp), intent(in) :: c, a, r
      integer, intent(in) :: l, k
      real(wp) :: b

      b = c*(r + 1.0_wp)**l*(2.0_wp*a*r + real(l + 1, wp))**k*exp(-a*r*r)
   end function prim_amp_bound

   !> Radius beyond which a whole shell's amplitude bound falls below ``tol``.
   !>
   !> Sums [[prim_amp_bound]] over the shell's primitives at ``max_deriv`` and
   !> bisects for the single outer crossing of
   !> ``tol``; the bound is large at r=0 and decays monotonically past its peak,
   !> so the crossing is unique.  ``tol <= 0`` disables screening (returns a huge
   !> reach).
   !>
   !> @param[in] l      Shell angular momentum
   !> @param[in] exps   Primitive exponents
   !> @param[in] coeffs Primitive contraction coefficients
   !> @param[in] tol    Amplitude screening threshold
   !> @param[in] max_deriv Highest spatial derivative order covered (0..4)
   !> @returns          Radial cutoff (Bohr)
   pure function iso_shell_rcut(l, exps, coeffs, tol, max_deriv) result(rcut)
      integer, intent(in) :: l
      real(wp), intent(in) :: exps(:), coeffs(:)
      real(wp), intent(in) :: tol
      integer, intent(in) :: max_deriv
      real(wp) :: rcut

      integer, parameter :: nbisect = 60
      real(wp), parameter :: r_cap = 1.0e3_wp
      real(wp) :: lo, hi, mid
      integer :: it

      if (max_deriv < 0 .or. max_deriv > 4) then
         error stop "Gaussian screening derivative order must be 0..4"
      end if

      if (tol <= 0.0_wp) then
         rcut = r_cap
         return
      end if

      !> Expand hi until the bound is below tol (or the cap is hit).
      hi = 1.0_wp
      do
         if (shell_amp(hi) <= tol .or. hi >= r_cap) exit
         hi = 2.0_wp*hi
      end do
      if (hi >= r_cap) then
         rcut = r_cap
         return
      end if

      !> Bisect [lo, hi] with amp(lo) > tol >= amp(hi) for the outer crossing.
      lo = 0.0_wp
      do it = 1, nbisect
         mid = 0.5_wp*(lo + hi)
         if (shell_amp(mid) > tol) then
            lo = mid
         else
            hi = mid
         end if
      end do
      rcut = hi
   contains
      pure function shell_amp(r) result(amp)
         real(wp), intent(in) :: r
         real(wp) :: amp
         integer :: p

         amp = 0.0_wp
         do p = 1, size(exps)
            amp = amp + prim_amp_bound(abs(coeffs(p)), exps(p), l, max_deriv, r)
         end do
      end function shell_amp
   end function iso_shell_rcut

   !> Recompute the per-shell squared radial cutoffs for a screening threshold
   !>
   !> ``threshold`` is an absolute amplitude cutoff on the (bounded) AO product contribution;
   !> ``threshold <= 0`` disables screening
   !>
   !> @param[inout] self      Basis instance
   !> @param[in]    threshold Amplitude screening threshold
   !> @param[in]    max_deriv Optional highest derivative order covered; default 3
   subroutine gto_build_screening(self, threshold, max_deriv)
      class(moist_iso_gto_type), intent(inout) :: self
      real(wp), intent(in) :: threshold
      integer, intent(in), optional :: max_deriv

      integer :: s, p0, p1, derivative_order
      real(wp) :: rcut

      derivative_order = 3
      if (present(max_deriv)) derivative_order = max_deriv
      self%max_rcut = 0.0_wp
      do s = 1, self%nshell
         p0 = self%sh_poff(s)
         p1 = self%sh_poff(s + 1) - 1
         rcut = iso_shell_rcut(self%sh_l(s), self%exps(p0:p1), self%coeffs(p0:p1), &
                               threshold, derivative_order)
         self%sh_rcut2(s) = rcut*rcut
         self%max_rcut = max(self%max_rcut, rcut)
      end do
   end subroutine gto_build_screening

   !> Conservative global shell reach (Bohr) for a screening threshold
   !>
   !> Independent of prior [[gto_build_screening]] calls
   !>
   !> @param[in] self      Basis instance
   !> @param[in] threshold Amplitude screening threshold
   !> @param[in] max_deriv Optional highest derivative order covered; default 3
   !> @returns             Largest shell cutoff radius over the whole basis
   pure function gto_reach(self, threshold, max_deriv) result(r)
      class(moist_iso_gto_type), intent(in) :: self
      real(wp), intent(in) :: threshold
      integer, intent(in), optional :: max_deriv
      real(wp) :: r

      integer :: s, p0, p1, derivative_order

      derivative_order = 3
      if (present(max_deriv)) derivative_order = max_deriv
      r = 0.0_wp
      do s = 1, self%nshell
         p0 = self%sh_poff(s)
         p1 = self%sh_poff(s + 1) - 1
         r = max(r, iso_shell_rcut(self%sh_l(s), self%exps(p0:p1), &
                                   self%coeffs(p0:p1), threshold, derivative_order))
      end do
   end function gto_reach

   !> Evaluate the electron density and its spatial derivatives at one point.
   !>
   !> The derivative order is set by *which* outputs the caller passes: only the
   !> orders whose output argument is present are computed, so the hot
   !> value+gradient path never pays for the Hessian or the third derivative.
   !> The caller supplies persistent per-thread scratch so that path performs no
   !> heap allocation either.
   !>
   !> Two levels of screening keep the cost proportional to the *nearby* basis
   !> rather than the whole molecule: when ``cand_atoms`` is present only the
   !> shells owned by those atoms (the cavity cell grid's candidate list) are
   !> visited, and every shell is further skipped when the point lies beyond its
   !> radial cutoff [[gto_type:sh_rcut2]].  The density contraction then runs over
   !> only the surviving ("active") cartesian components -- ``nact`` of them --
   !> so it is O(nact^2) instead of dense O(ncart^2).  With screening disabled
   !> (huge cutoffs) every component is active and the result is the exact dense
   !> contraction.
   !>
   !> @param[in]    self       Basis instance (read-only)
   !> @param[in]    point      Evaluation point in Bohr
   !> @param[inout] phi        Scratch AO-derivative table, shape (ncart, 0:34);
   !>                          only the slots of the requested orders are touched
   !> @param[inout] t0         Scratch, shape (ncart): dcart . phi(:,value)
   !> @param[inout] tm         Scratch, shape (ncart, 3): dcart . phi(:,grad)
   !> @param[inout] act        Scratch active-component index list, shape (ncart)
   !> @param[out]   rho        Electron density
   !> @param[out]   drho       Density gradient (3)
   !> @param[out]   d2rho      Optional density Hessian (3, 3)
   !> @param[out]   d3rho      Optional density third derivative (3, 3, 3)
   !> @param[in]    cand_atoms Optional candidate owner-atom list; when present,
   !>                          only shells on these atoms are considered
   !> @param[out]   d4rho      Optional density fourth derivative (3, 3, 3, 3);
   !>                          requires the ``tmm`` scratch
   !> @param[inout] tmm        Optional scratch, shape (ncart, 6): dcart . phi(:,Hess);
   !>                          required when ``d4rho`` is present
   subroutine gto_eval(self, point, phi, t0, tm, act, &
                       rho, drho, d2rho, d3rho, &
                       cand_atoms, d4rho, tmm)
      !> Basis instance
      class(moist_iso_gto_type), intent(in) :: self
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Scratch AO-derivative table (ncart, 0:34)
      real(wp), intent(inout) :: phi(:, 0:)
      !> Scratch density-weighted value vector (ncart)
      real(wp), intent(inout) :: t0(:)
      !> Scratch density-weighted gradient vectors (ncart, 3)
      real(wp), intent(inout) :: tm(:, :)
      !> Scratch active-component index list (ncart)
      integer, intent(inout) :: act(:)
      !> Electron density
      real(wp), intent(out) :: rho
      !> Density gradient
      real(wp), intent(out) :: drho(3)
      !> Optional density Hessian
      real(wp), intent(out), optional :: d2rho(3, 3)
      !> Optional density third derivative
      real(wp), intent(out), optional :: d3rho(3, 3, 3)
      !> Optional candidate owner-atom list (cavity cell-grid screening)
      integer, intent(in), optional :: cand_atoms(:)
      !> Optional density fourth derivative
      real(wp), intent(out), optional :: d4rho(3, 3, 3, 3)
      !> Optional density-weighted Hessian vectors (ncart, 6)
      real(wp), intent(inout), optional :: tmm(:, :)

      integer :: s, ka, aidx, idx, nslot, nact, ndloc
      integer :: i, j, k, l, m, ii, jj, ci
      real(wp) :: acc, s1, s2
      logical :: use_blas, want_d2, want_d3, want_d4

      want_d2 = present(d2rho)
      want_d3 = present(d3rho)
      want_d4 = present(d4rho)

      ! The requested order is exactly the highest order the caller asked back.
      ! The fourth derivative additionally needs the extra density-weighted
      ! Hessian scratch.
      ndloc = 1
      if (want_d2) ndloc = 2
      if (want_d3) ndloc = 3
      if (want_d4) then
         ndloc = 4
         if (.not. present(tmm)) error stop "gto_eval: d4rho requires the tmm scratch"
      end if
      nslot = nslot_by_order(ndloc)
      if (size(phi, 2) < nslot) error stop "gto_eval: phi scratch has too few derivative slots"
      rho = 0.0_wp
      drho = 0.0_wp
      if (want_d2) d2rho = 0.0_wp
      if (want_d3) d3rho = 0.0_wp
      if (want_d4) d4rho = 0.0_wp
      nact = 0

      !> Zero the whole used table once
      phi(:, 0:nslot - 1) = 0.0_wp

      ! Assemble the AO-derivative table for the surviving shells only
      !+ record their cartesian components in the active list
      if (present(cand_atoms)) then
         do ka = 1, size(cand_atoms)
            aidx = cand_atoms(ka)
            if (aidx < 1 .or. aidx > self%natom_grid) cycle
            do idx = self%atom_soff(aidx) + 1, self%atom_soff(aidx + 1)
               call assemble_shell(self%atom_shell(idx))
            end do
         end do
      else
         do s = 1, self%nshell
            call assemble_shell(s)
         end do
      end if

      if (nact == 0) return

      ! Contract with the cartesian-monomial density matrix.  rho = a^T D a,
      ! d_i rho = 2 (D a) . a^i, and the higher orders follow by the product
      ! rule (see module header), tm(:,j) = D a^j.
      !
      ! Two regimes: when screening left most components active a dense matmul
      ! over the (zeroed) full table is more cache-friendly than gathering the
      ! active block; once screening removes the bulk (nact well below ncart)
      ! the O(nact^2) active-only contraction wins.  Both give the same result.
      if (3*nact > self%ncart) then
         ! Dense contraction (little screened).  dcart is symmetric and the
         ! inactive components are zero, so the full matvec ``D a`` is exact;
         ! [[dcart_matvec]] routes it through BLAS dsymv for large solutes (where
         ! this path dominates) and the matmul intrinsic for small ones.
         use_blas = self%ncart >= iso_blas_ncart_min
         call dcart_matvec(phi(:, 0), t0)
         rho = dot_product(phi(:, 0), t0)
         drho(1) = 2.0_wp*dot_product(phi(:, 1), t0)
         drho(2) = 2.0_wp*dot_product(phi(:, 2), t0)
         drho(3) = 2.0_wp*dot_product(phi(:, 3), t0)

         ! The gradient weights tm are needed by every order above the first,
         ! so they follow ndloc rather than the presence of d2rho itself.
         if (ndloc >= 2) then
            call dcart_matvec(phi(:, 1), tm(:, 1))
            call dcart_matvec(phi(:, 2), tm(:, 2))
            call dcart_matvec(phi(:, 3), tm(:, 3))
            if (want_d2) then
               do j = 1, 3
                  do i = 1, 3
                     d2rho(i, j) = 2.0_wp*(dot_product(phi(:, hess_slot(i, j)), t0) &
                                           + dot_product(phi(:, i), tm(:, j)))
                  end do
               end do
            end if
         end if

         if (want_d3) then
            do k = 1, 3
               do j = 1, 3
                  do i = 1, 3
                     d3rho(i, j, k) = 2.0_wp*( &
                                      dot_product(phi(:, third_slot(i, j, k)), t0) &
                                      + dot_product(phi(:, hess_slot(j, k)), tm(:, i)) &
                                      + dot_product(phi(:, hess_slot(i, k)), tm(:, j)) &
                                      + dot_product(phi(:, hess_slot(i, j)), tm(:, k)))
                  end do
               end do
            end do
         end if
      else
         do ii = 1, nact
            ci = act(ii)
            acc = 0.0_wp
            do jj = 1, nact
               acc = acc + self%dcart(ci, act(jj))*phi(act(jj), 0)
            end do
            t0(ci) = acc
         end do
         do ii = 1, nact
            ci = act(ii)
            rho = rho + phi(ci, 0)*t0(ci)
            drho(1) = drho(1) + phi(ci, 1)*t0(ci)
            drho(2) = drho(2) + phi(ci, 2)*t0(ci)
            drho(3) = drho(3) + phi(ci, 3)*t0(ci)
         end do
         drho = 2.0_wp*drho

         if (ndloc >= 2) then
            do j = 1, 3
               do ii = 1, nact
                  ci = act(ii)
                  acc = 0.0_wp
                  do jj = 1, nact
                     acc = acc + self%dcart(ci, act(jj))*phi(act(jj), j)
                  end do
                  tm(ci, j) = acc
               end do
            end do
            if (want_d2) then
               do j = 1, 3
                  do i = 1, 3
                     s1 = 0.0_wp
                     s2 = 0.0_wp
                     do ii = 1, nact
                        ci = act(ii)
                        s1 = s1 + phi(ci, hess_slot(i, j))*t0(ci)
                        s2 = s2 + phi(ci, i)*tm(ci, j)
                     end do
                     d2rho(i, j) = 2.0_wp*(s1 + s2)
                  end do
               end do
            end if
         end if

         if (want_d3) then
            do k = 1, 3
               do j = 1, 3
                  do i = 1, 3
                     s1 = 0.0_wp
                     do ii = 1, nact
                        ci = act(ii)
                        s1 = s1 + phi(ci, third_slot(i, j, k))*t0(ci) &
                             + phi(ci, hess_slot(j, k))*tm(ci, i) &
                             + phi(ci, hess_slot(i, k))*tm(ci, j) &
                             + phi(ci, hess_slot(i, j))*tm(ci, k)
                     end do
                     d3rho(i, j, k) = 2.0_wp*s1
                  end do
               end do
            end do
         end if
      end if

      ! Fourth derivative (optional).  Reuses the active list, so it is correct
      ! after either contraction branch (inactive phi are zero).  First form the
      ! density-weighted Hessian vectors tmm(:,m) = D phi(:,Hess m), then the
      ! Leibniz expansion of the fourth derivative of rho = phi^T D phi: one
      ! 4th-order term, four (3rd x 1st) terms, and three (Hess x Hess) terms.
      if (want_d4) then
         do m = 1, 6
            do ii = 1, nact
               ci = act(ii)
               acc = 0.0_wp
               do jj = 1, nact
                  acc = acc + self%dcart(ci, act(jj))*phi(act(jj), 3 + m)
               end do
               tmm(ci, m) = acc
            end do
         end do
         do l = 1, 3
            do k = 1, 3
               do j = 1, 3
                  do i = 1, 3
                     s1 = 0.0_wp
                     do ii = 1, nact
                        ci = act(ii)
                        s1 = s1 + phi(ci, fourth_slot(i, j, k, l))*t0(ci) &
                             + phi(ci, third_slot(j, k, l))*tm(ci, i) &
                             + phi(ci, third_slot(i, k, l))*tm(ci, j) &
                             + phi(ci, third_slot(i, j, l))*tm(ci, k) &
                             + phi(ci, third_slot(i, j, k))*tm(ci, l) &
                             + phi(ci, hess_slot(i, j))*tmm(ci, hess_slot(k, l) - 3) &
                             + phi(ci, hess_slot(i, k))*tmm(ci, hess_slot(j, l) - 3) &
                             + phi(ci, hess_slot(i, l))*tmm(ci, hess_slot(j, k) - 3)
                     end do
                     d4rho(i, j, k, l) = 2.0_wp*s1
                  end do
               end do
            end do
         end do
      end if

   contains

      !> Accumulate one shell's cartesian AO-derivative contributions into phi
      !> (skipping shells beyond their radial cutoff) and append its components
      !> to the active list.  Host-associated: reads point/ndloc/self, writes
      !> phi/act/nact.  Delegates the actual slot fills to the shared module
      !> routine [[gto_assemble_shell]] so the batched evaluator stays in sync.
      subroutine assemble_shell(sh)
         integer, intent(in) :: sh
         integer :: c, cbase, ncc
         logical :: contributed

         call gto_assemble_shell(self, sh, point, ndloc, phi, contributed)
         if (.not. contributed) return

         cbase = self%sh_coff(sh)
         ncc = moist_iso_gto_ncart(self%sh_l(sh))
         do c = 1, ncc
            nact = nact + 1
            act(nact) = cbase + c
         end do
      end subroutine assemble_shell

      !> Dense density matvec ``y = dcart . x`` over all ncart components, via
      !> BLAS dsymv (large solutes) or the matmul intrinsic (small).  Host
      !> associated: reads self%dcart and the host-set ``use_blas`` flag.
      subroutine dcart_matvec(x, y)
         real(wp), intent(in) :: x(:)
         real(wp), intent(inout) :: y(:)

         if (use_blas) then
            call symv(self%dcart, x, y)
         else
            y(1:self%ncart) = matmul(self%dcart, x)
         end if
      end subroutine dcart_matvec

   end subroutine gto_eval

   !> Accumulate one shell's cartesian AO-derivative contributions for a single
   !> point into the (single-point) derivative table ``phi(ncart, 0:nslot-1)``
   !>
   !> Shared by the single-point [[gto_eval]] and the batched [[gto_eval_batch]]
   !> so the (screening + primitive loop + per-slot factor products) live in one
   !> place.  Skips the shell when the point is beyond its squared radial cutoff,
   !> reporting that through ``contributed`` so the caller can decide whether to
   !> record the shell's components in an active list.
   !>
   !> @param[in]    self        Basis instance
   !> @param[in]    sh          Shell index
   !> @param[in]    point       Evaluation point in Bohr
   !> @param[in]    ndloc       Highest spatial derivative order to fill (1..4)
   !> @param[inout] phi         Derivative table for this point (accumulated into)
   !> @param[out]   contributed .true. if the shell was inside its radial cutoff
   subroutine gto_assemble_shell(self, sh, point, ndloc, phi, contributed)
      class(moist_iso_gto_type), intent(in) :: self
      integer, intent(in) :: sh
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: ndloc
      real(wp), intent(inout) :: phi(:, 0:)
      logical, intent(out) :: contributed

      integer :: p, c, cbase, l, ncc, lx, ly, lz
      real(wp) :: dx, dy, dz, u, a, ce
      real(wp) :: px(0:moist_iso_gto_lmax, 0:4)
      real(wp) :: py(0:moist_iso_gto_lmax, 0:4)
      real(wp) :: pz(0:moist_iso_gto_lmax, 0:4)

      contributed = .false.

      dx = point(1) - self%center(1, sh)
      dy = point(2) - self%center(2, sh)
      dz = point(3) - self%center(3, sh)
      u = dx*dx + dy*dy + dz*dz
      if (u > self%sh_rcut2(sh)) return

      contributed = .true.
      l = self%sh_l(sh)
      ncc = moist_iso_gto_ncart(l)
      cbase = self%sh_coff(sh)
      do p = self%sh_poff(sh), self%sh_poff(sh + 1) - 1
         a = self%exps(p)
         ce = self%coeffs(p)*exp(-a*u)
         if (ce == 0.0_wp) cycle
         call moist_iso_gto_poly1d(dx, a, ndloc, l, px)
         call moist_iso_gto_poly1d(dy, a, ndloc, l, py)
         call moist_iso_gto_poly1d(dz, a, ndloc, l, pz)
         do c = 1, ncc
            lx = self%comp_l(1, cbase + c)
            ly = self%comp_l(2, cbase + c)
            lz = self%comp_l(3, cbase + c)
            phi(cbase + c, 0) = phi(cbase + c, 0) + ce*px(lx, 0)*py(ly, 0)*pz(lz, 0)
            phi(cbase + c, 1) = phi(cbase + c, 1) + ce*px(lx, 1)*py(ly, 0)*pz(lz, 0)
            phi(cbase + c, 2) = phi(cbase + c, 2) + ce*px(lx, 0)*py(ly, 1)*pz(lz, 0)
            phi(cbase + c, 3) = phi(cbase + c, 3) + ce*px(lx, 0)*py(ly, 0)*pz(lz, 1)
            if (ndloc >= 2) then
               phi(cbase + c, 4) = phi(cbase + c, 4) + ce*px(lx, 2)*py(ly, 0)*pz(lz, 0)
               phi(cbase + c, 5) = phi(cbase + c, 5) + ce*px(lx, 1)*py(ly, 1)*pz(lz, 0)
               phi(cbase + c, 6) = phi(cbase + c, 6) + ce*px(lx, 1)*py(ly, 0)*pz(lz, 1)
               phi(cbase + c, 7) = phi(cbase + c, 7) + ce*px(lx, 0)*py(ly, 2)*pz(lz, 0)
               phi(cbase + c, 8) = phi(cbase + c, 8) + ce*px(lx, 0)*py(ly, 1)*pz(lz, 1)
               phi(cbase + c, 9) = phi(cbase + c, 9) + ce*px(lx, 0)*py(ly, 0)*pz(lz, 2)
            end if
            if (ndloc >= 3) then
               phi(cbase + c, 10) = phi(cbase + c, 10) + ce*px(lx, 3)*py(ly, 0)*pz(lz, 0)
               phi(cbase + c, 11) = phi(cbase + c, 11) + ce*px(lx, 2)*py(ly, 1)*pz(lz, 0)
               phi(cbase + c, 12) = phi(cbase + c, 12) + ce*px(lx, 2)*py(ly, 0)*pz(lz, 1)
               phi(cbase + c, 13) = phi(cbase + c, 13) + ce*px(lx, 1)*py(ly, 2)*pz(lz, 0)
               phi(cbase + c, 14) = phi(cbase + c, 14) + ce*px(lx, 1)*py(ly, 1)*pz(lz, 1)
               phi(cbase + c, 15) = phi(cbase + c, 15) + ce*px(lx, 1)*py(ly, 0)*pz(lz, 2)
               phi(cbase + c, 16) = phi(cbase + c, 16) + ce*px(lx, 0)*py(ly, 3)*pz(lz, 0)
               phi(cbase + c, 17) = phi(cbase + c, 17) + ce*px(lx, 0)*py(ly, 2)*pz(lz, 1)
               phi(cbase + c, 18) = phi(cbase + c, 18) + ce*px(lx, 0)*py(ly, 1)*pz(lz, 2)
               phi(cbase + c, 19) = phi(cbase + c, 19) + ce*px(lx, 0)*py(ly, 0)*pz(lz, 3)
            end if
            if (ndloc >= 4) then
               phi(cbase + c, 20) = phi(cbase + c, 20) + ce*px(lx, 4)*py(ly, 0)*pz(lz, 0)
               phi(cbase + c, 21) = phi(cbase + c, 21) + ce*px(lx, 3)*py(ly, 1)*pz(lz, 0)
               phi(cbase + c, 22) = phi(cbase + c, 22) + ce*px(lx, 3)*py(ly, 0)*pz(lz, 1)
               phi(cbase + c, 23) = phi(cbase + c, 23) + ce*px(lx, 2)*py(ly, 2)*pz(lz, 0)
               phi(cbase + c, 24) = phi(cbase + c, 24) + ce*px(lx, 2)*py(ly, 1)*pz(lz, 1)
               phi(cbase + c, 25) = phi(cbase + c, 25) + ce*px(lx, 2)*py(ly, 0)*pz(lz, 2)
               phi(cbase + c, 26) = phi(cbase + c, 26) + ce*px(lx, 1)*py(ly, 3)*pz(lz, 0)
               phi(cbase + c, 27) = phi(cbase + c, 27) + ce*px(lx, 1)*py(ly, 2)*pz(lz, 1)
               phi(cbase + c, 28) = phi(cbase + c, 28) + ce*px(lx, 1)*py(ly, 1)*pz(lz, 2)
               phi(cbase + c, 29) = phi(cbase + c, 29) + ce*px(lx, 1)*py(ly, 0)*pz(lz, 3)
               phi(cbase + c, 30) = phi(cbase + c, 30) + ce*px(lx, 0)*py(ly, 4)*pz(lz, 0)
               phi(cbase + c, 31) = phi(cbase + c, 31) + ce*px(lx, 0)*py(ly, 3)*pz(lz, 1)
               phi(cbase + c, 32) = phi(cbase + c, 32) + ce*px(lx, 0)*py(ly, 2)*pz(lz, 2)
               phi(cbase + c, 33) = phi(cbase + c, 33) + ce*px(lx, 0)*py(ly, 1)*pz(lz, 3)
               phi(cbase + c, 34) = phi(cbase + c, 34) + ce*px(lx, 0)*py(ly, 0)*pz(lz, 4)
            end if
         end do
      end do
   end subroutine gto_assemble_shell

   ! >>> GENERATED poly1d >>>
   !> Fill the exp-stripped 1D Gaussian-factor derivative table
   !>
   !> GENERATED by gen_isodensity_gto.py
   !>
   !> ``p(n, k) = [ d^k/dt^k ( t**n exp(-a t**2) ) ] / exp(-a t**2)`` for factor
   !> power ``n`` (0..nmax) and derivative order ``k`` (0..nderiv)
   !> Rows above ``nmax`` and columns above ``nderiv`` are untouched
   !>
   !> @param[in]  t       Displacement component (Bohr)
   !> @param[in]  a       Primitive Gaussian exponent
   !> @param[in]  nderiv  Highest derivative order to fill (1..4)
   !> @param[in]  nmax    Highest factor power (row) to fill = the shell's l
   !> @param[out] p       Factor-derivative table, shape (0:lmax, 0:4)
   pure subroutine moist_iso_gto_poly1d(t, a, nderiv, nmax, p)
      !> Displacement component
      real(wp), intent(in) :: t
      !> Primitive Gaussian exponent
      real(wp), intent(in) :: a
      !> Highest derivative order to fill
      integer, intent(in) :: nderiv
      !> Highest factor power to fill (the shell's angular momentum)
      integer, intent(in) :: nmax
      !> Exp-stripped factor-derivative table
      real(wp), intent(out) :: p(0:8, 0:4)

      select case (nmax)
      case (0)
         block
            ! value (k=0)
            p(0, 0) = 1.0_wp
            ! first derivative (k=1)
            p(0, 1) = -2*a*t
            if (nderiv >= 2) then
               p(0, 2) = 4*a**2*t**2 - 2*a
            end if
            if (nderiv >= 3) then
               p(0, 3) = -8*a**3*t**3 + 12*a**2*t
            end if
            if (nderiv >= 4) then
               p(0, 4) = 16*a**4*t**4 - 48*a**3*t**2 + 12*a**2
            end if
         end block
      case (1)
         block
            real(wp) :: c1_1_0, c1_2_0, c1_3_0, c1_3_1, c1_4_0, c1_4_1, c1_4_2
            ! value (k=0)
            p(0, 0) = 1.0_wp
            p(1, 0) = t
            ! first derivative (k=1)
            c1_1_0 = 2*a
            p(0, 1) = -c1_1_0*t
            p(1, 1) = -c1_1_0*t**2 + 1
            if (nderiv >= 2) then
               c1_2_0 = 4*a**2
               p(0, 2) = -2*a + c1_2_0*t**2
               p(1, 2) = -6*a*t + c1_2_0*t**3
            end if
            if (nderiv >= 3) then
               c1_3_0 = a**2
               c1_3_1 = 8*a**3
               p(0, 3) = 12*c1_3_0*t - c1_3_1*t**3
               p(1, 3) = -6*a + 24*c1_3_0*t**2 - c1_3_1*t**4
            end if
            if (nderiv >= 4) then
               c1_4_0 = a**2
               c1_4_1 = a**3
               c1_4_2 = 16*a**4
               p(0, 4) = 12*c1_4_0 - 48*c1_4_1*t**2 + c1_4_2*t**4
               p(1, 4) = 60*c1_4_0*t - 80*c1_4_1*t**3 + c1_4_2*t**5
            end if
         end block
      case (2)
         block
            real(wp) :: c2_1_0, c2_2_0, c2_2_1, c2_3_0, c2_3_1, c2_3_2, c2_4_0, c2_4_1
            real(wp) :: c2_4_2, c2_4_3, c2_4_4
            ! value (k=0)
            p(0, 0) = 1.0_wp
            p(1, 0) = t
            p(2, 0) = t**2
            ! first derivative (k=1)
            c2_1_0 = 2*a
            p(0, 1) = -2*a*t
            p(1, 1) = -c2_1_0*t**2 + 1
            p(2, 1) = -c2_1_0*t**3 + 2*t
            if (nderiv >= 2) then
               c2_2_0 = t**2
               c2_2_1 = 4*a**2
               p(0, 2) = -2*a + c2_2_0*c2_2_1
               p(1, 2) = -6*a*t + c2_2_1*t**3
               p(2, 2) = -10*a*c2_2_0 + c2_2_1*t**4 + 2
            end if
            if (nderiv >= 3) then
               c2_3_0 = a**2
               c2_3_1 = t**3
               c2_3_2 = 8*a**3
               p(0, 3) = 12*c2_3_0*t - c2_3_1*c2_3_2
               p(1, 3) = -6*a + 24*c2_3_0*t**2 - c2_3_2*t**4
               p(2, 3) = -24*a*t + 36*c2_3_0*c2_3_1 - c2_3_2*t**5
            end if
            if (nderiv >= 4) then
               c2_4_0 = a**2
               c2_4_1 = a**3
               c2_4_2 = t**2
               c2_4_3 = t**4
               c2_4_4 = 16*a**4
               p(0, 4) = 12*c2_4_0 - 48*c2_4_1*c2_4_2 + c2_4_3*c2_4_4
               p(1, 4) = 60*c2_4_0*t - 80*c2_4_1*t**3 + c2_4_4*t**5
               p(2, 4) = -24*a + 156*c2_4_0*c2_4_2 - 112*c2_4_1*c2_4_3 + c2_4_4*t**6
            end if
         end block
      case (3)
         block
            real(wp) :: c3_1_0, c3_1_1, c3_2_0, c3_2_1, c3_2_2, c3_2_3, c3_3_0, c3_3_1
            real(wp) :: c3_3_2, c3_3_3, c3_3_4, c3_4_0, c3_4_1, c3_4_2, c3_4_3, c3_4_4
            real(wp) :: c3_4_5, c3_4_6
            ! value (k=0)
            p(0, 0) = 1.0_wp
            p(1, 0) = t
            p(2, 0) = t**2
            p(3, 0) = t**3
            ! first derivative (k=1)
            c3_1_0 = t**2
            c3_1_1 = 2*a
            p(0, 1) = -2*a*t
            p(1, 1) = -c3_1_0*c3_1_1 + 1
            p(2, 1) = -c3_1_1*t**3 + 2*t
            p(3, 1) = 3*c3_1_0 - c3_1_1*t**4
            if (nderiv >= 2) then
               c3_2_0 = t**2
               c3_2_1 = 4*a**2
               c3_2_2 = 6*t
               c3_2_3 = t**3
               p(0, 2) = -2*a + c3_2_0*c3_2_1
               p(1, 2) = -a*c3_2_2 + c3_2_1*c3_2_3
               p(2, 2) = -10*a*c3_2_0 + c3_2_1*t**4 + 2
               p(3, 2) = -14*a*c3_2_3 + c3_2_1*t**5 + c3_2_2
            end if
            if (nderiv >= 3) then
               c3_3_0 = a**2
               c3_3_1 = t**3
               c3_3_2 = 8*a**3
               c3_3_3 = t**2
               c3_3_4 = t**4
               p(0, 3) = 12*c3_3_0*t - c3_3_1*c3_3_2
               p(1, 3) = -6*a + 24*c3_3_0*c3_3_3 - c3_3_2*c3_3_4
               p(2, 3) = -24*a*t + 36*c3_3_0*c3_3_1 - c3_3_2*t**5
               p(3, 3) = -54*a*c3_3_3 + 48*c3_3_0*c3_3_4 - c3_3_2*t**6 + 6
            end if
            if (nderiv >= 4) then
               c3_4_0 = a**2
               c3_4_1 = a**3
               c3_4_2 = t**2
               c3_4_3 = t**4
               c3_4_4 = 16*a**4
               c3_4_5 = t**3
               c3_4_6 = t**5
               p(0, 4) = 12*c3_4_0 - 48*c3_4_1*c3_4_2 + c3_4_3*c3_4_4
               p(1, 4) = 60*c3_4_0*t - 80*c3_4_1*c3_4_5 + c3_4_4*c3_4_6
               p(2, 4) = -24*a + 156*c3_4_0*c3_4_2 - 112*c3_4_1*c3_4_3 + c3_4_4*t**6
               p(3, 4) = -120*a*t + 300*c3_4_0*c3_4_5 - 144*c3_4_1*c3_4_6 + c3_4_4*t**7
            end if
         end block
      case (4)
         block
            real(wp) :: c4_1_0, c4_1_1, c4_1_2, c4_2_0, c4_2_1, c4_2_2, c4_2_3, c4_2_4
            real(wp) :: c4_3_0, c4_3_1, c4_3_2, c4_3_3, c4_3_4, c4_3_5, c4_4_0, c4_4_1
            real(wp) :: c4_4_2, c4_4_3, c4_4_4, c4_4_5, c4_4_6, c4_4_7
            ! value (k=0)
            p(0, 0) = 1.0_wp
            p(1, 0) = t
            p(2, 0) = t**2
            p(3, 0) = t**3
            p(4, 0) = t**4
            ! first derivative (k=1)
            c4_1_0 = t**2
            c4_1_1 = 2*a
            c4_1_2 = t**3
            p(0, 1) = -2*a*t
            p(1, 1) = -c4_1_0*c4_1_1 + 1
            p(2, 1) = -c4_1_1*c4_1_2 + 2*t
            p(3, 1) = 3*c4_1_0 - c4_1_1*t**4
            p(4, 1) = -c4_1_1*t**5 + 4*c4_1_2
            if (nderiv >= 2) then
               c4_2_0 = t**2
               c4_2_1 = 4*a**2
               c4_2_2 = 6*t
               c4_2_3 = t**3
               c4_2_4 = t**4
               p(0, 2) = -2*a + c4_2_0*c4_2_1
               p(1, 2) = -a*c4_2_2 + c4_2_1*c4_2_3
               p(2, 2) = -10*a*c4_2_0 + c4_2_1*c4_2_4 + 2
               p(3, 2) = -14*a*c4_2_3 + c4_2_1*t**5 + c4_2_2
               p(4, 2) = -18*a*c4_2_4 + 12*c4_2_0 + c4_2_1*t**6
            end if
            if (nderiv >= 3) then
               c4_3_0 = a**2
               c4_3_1 = t**3
               c4_3_2 = 8*a**3
               c4_3_3 = t**2
               c4_3_4 = t**4
               c4_3_5 = t**5
               p(0, 3) = 12*c4_3_0*t - c4_3_1*c4_3_2
               p(1, 3) = -6*a + 24*c4_3_0*c4_3_3 - c4_3_2*c4_3_4
               p(2, 3) = -24*a*t + 36*c4_3_0*c4_3_1 - c4_3_2*c4_3_5
               p(3, 3) = -54*a*c4_3_3 + 48*c4_3_0*c4_3_4 - c4_3_2*t**6 + 6
               p(4, 3) = -96*a*c4_3_1 + 60*c4_3_0*c4_3_5 - c4_3_2*t**7 + 24*t
            end if
            if (nderiv >= 4) then
               c4_4_0 = a**2
               c4_4_1 = a**3
               c4_4_2 = t**2
               c4_4_3 = t**4
               c4_4_4 = 16*a**4
               c4_4_5 = t**3
               c4_4_6 = t**5
               c4_4_7 = t**6
               p(0, 4) = 12*c4_4_0 - 48*c4_4_1*c4_4_2 + c4_4_3*c4_4_4
               p(1, 4) = 60*c4_4_0*t - 80*c4_4_1*c4_4_5 + c4_4_4*c4_4_6
               p(2, 4) = -24*a + 156*c4_4_0*c4_4_2 - 112*c4_4_1*c4_4_3 + c4_4_4*c4_4_7
               p(3, 4) = -120*a*t + 300*c4_4_0*c4_4_5 - 144*c4_4_1*c4_4_6 + c4_4_4*t**7
               p(4, 4) = -336*a*c4_4_2 + 492*c4_4_0*c4_4_3 - 176*c4_4_1*c4_4_7 + c4_4_4*t**8 + 24
            end if
         end block
      case (5)
         block
            real(wp) :: c5_1_0, c5_1_1, c5_1_2, c5_1_3, c5_2_0, c5_2_1, c5_2_2, c5_2_3
            real(wp) :: c5_2_4, c5_2_5, c5_3_0, c5_3_1, c5_3_2, c5_3_3, c5_3_4, c5_3_5
            real(wp) :: c5_3_6, c5_4_0, c5_4_1, c5_4_2, c5_4_3, c5_4_4, c5_4_5, c5_4_6
            real(wp) :: c5_4_7, c5_4_8, c5_4_9
            ! value (k=0)
            p(0, 0) = 1.0_wp
            p(1, 0) = t
            p(2, 0) = t**2
            p(3, 0) = t**3
            p(4, 0) = t**4
            p(5, 0) = t**5
            ! first derivative (k=1)
            c5_1_0 = t**2
            c5_1_1 = 2*a
            c5_1_2 = t**3
            c5_1_3 = t**4
            p(0, 1) = -2*a*t
            p(1, 1) = -c5_1_0*c5_1_1 + 1
            p(2, 1) = -c5_1_1*c5_1_2 + 2*t
            p(3, 1) = 3*c5_1_0 - c5_1_1*c5_1_3
            p(4, 1) = -c5_1_1*t**5 + 4*c5_1_2
            p(5, 1) = -c5_1_1*t**6 + 5*c5_1_3
            if (nderiv >= 2) then
               c5_2_0 = t**2
               c5_2_1 = 4*a**2
               c5_2_2 = 6*t
               c5_2_3 = t**3
               c5_2_4 = t**4
               c5_2_5 = t**5
               p(0, 2) = -2*a + c5_2_0*c5_2_1
               p(1, 2) = -a*c5_2_2 + c5_2_1*c5_2_3
               p(2, 2) = -10*a*c5_2_0 + c5_2_1*c5_2_4 + 2
               p(3, 2) = -14*a*c5_2_3 + c5_2_1*c5_2_5 + c5_2_2
               p(4, 2) = -18*a*c5_2_4 + 12*c5_2_0 + c5_2_1*t**6
               p(5, 2) = -22*a*c5_2_5 + c5_2_1*t**7 + 20*c5_2_3
            end if
            if (nderiv >= 3) then
               c5_3_0 = a**2
               c5_3_1 = t**3
               c5_3_2 = 8*a**3
               c5_3_3 = t**2
               c5_3_4 = t**4
               c5_3_5 = t**5
               c5_3_6 = t**6
               p(0, 3) = 12*c5_3_0*t - c5_3_1*c5_3_2
               p(1, 3) = -6*a + 24*c5_3_0*c5_3_3 - c5_3_2*c5_3_4
               p(2, 3) = -24*a*t + 36*c5_3_0*c5_3_1 - c5_3_2*c5_3_5
               p(3, 3) = -54*a*c5_3_3 + 48*c5_3_0*c5_3_4 - c5_3_2*c5_3_6 + 6
               p(4, 3) = -96*a*c5_3_1 + 60*c5_3_0*c5_3_5 - c5_3_2*t**7 + 24*t
               p(5, 3) = -150*a*c5_3_4 + 72*c5_3_0*c5_3_6 - c5_3_2*t**8 + 60*c5_3_3
            end if
            if (nderiv >= 4) then
               c5_4_0 = a**2
               c5_4_1 = a**3
               c5_4_2 = t**2
               c5_4_3 = t**4
               c5_4_4 = 16*a**4
               c5_4_5 = t**3
               c5_4_6 = t**5
               c5_4_7 = t**6
               c5_4_8 = 120*t
               c5_4_9 = t**7
               p(0, 4) = 12*c5_4_0 - 48*c5_4_1*c5_4_2 + c5_4_3*c5_4_4
               p(1, 4) = 60*c5_4_0*t - 80*c5_4_1*c5_4_5 + c5_4_4*c5_4_6
               p(2, 4) = -24*a + 156*c5_4_0*c5_4_2 - 112*c5_4_1*c5_4_3 + c5_4_4*c5_4_7
               p(3, 4) = -a*c5_4_8 + 300*c5_4_0*c5_4_5 - 144*c5_4_1*c5_4_6 + c5_4_4*c5_4_9
               p(4, 4) = -336*a*c5_4_2 + 492*c5_4_0*c5_4_3 - 176*c5_4_1*c5_4_7 + c5_4_4*t**8 + 24
               p(5, 4) = -720*a*c5_4_5 + 732*c5_4_0*c5_4_6 - 208*c5_4_1*c5_4_9 + c5_4_4*t**9 + c5_4_8
            end if
         end block
      case (6)
         block
            real(wp) :: c6_1_0, c6_1_1, c6_1_2, c6_1_3, c6_1_4, c6_2_0, c6_2_1, c6_2_2
            real(wp) :: c6_2_3, c6_2_4, c6_2_5, c6_2_6, c6_3_0, c6_3_1, c6_3_2, c6_3_3
            real(wp) :: c6_3_4, c6_3_5, c6_3_6, c6_3_7, c6_4_0, c6_4_1, c6_4_2, c6_4_3
            real(wp) :: c6_4_4, c6_4_5, c6_4_6, c6_4_7, c6_4_8, c6_4_9, c6_4_10
            ! value (k=0)
            p(0, 0) = 1.0_wp
            p(1, 0) = t
            p(2, 0) = t**2
            p(3, 0) = t**3
            p(4, 0) = t**4
            p(5, 0) = t**5
            p(6, 0) = t**6
            ! first derivative (k=1)
            c6_1_0 = t**2
            c6_1_1 = 2*a
            c6_1_2 = t**3
            c6_1_3 = t**4
            c6_1_4 = t**5
            p(0, 1) = -2*a*t
            p(1, 1) = -c6_1_0*c6_1_1 + 1
            p(2, 1) = -c6_1_1*c6_1_2 + 2*t
            p(3, 1) = 3*c6_1_0 - c6_1_1*c6_1_3
            p(4, 1) = -c6_1_1*c6_1_4 + 4*c6_1_2
            p(5, 1) = -c6_1_1*t**6 + 5*c6_1_3
            p(6, 1) = -c6_1_1*t**7 + 6*c6_1_4
            if (nderiv >= 2) then
               c6_2_0 = t**2
               c6_2_1 = 4*a**2
               c6_2_2 = 6*t
               c6_2_3 = t**3
               c6_2_4 = t**4
               c6_2_5 = t**5
               c6_2_6 = t**6
               p(0, 2) = -2*a + c6_2_0*c6_2_1
               p(1, 2) = -a*c6_2_2 + c6_2_1*c6_2_3
               p(2, 2) = -10*a*c6_2_0 + c6_2_1*c6_2_4 + 2
               p(3, 2) = -14*a*c6_2_3 + c6_2_1*c6_2_5 + c6_2_2
               p(4, 2) = -18*a*c6_2_4 + 12*c6_2_0 + c6_2_1*c6_2_6
               p(5, 2) = -22*a*c6_2_5 + c6_2_1*t**7 + 20*c6_2_3
               p(6, 2) = -26*a*c6_2_6 + c6_2_1*t**8 + 30*c6_2_4
            end if
            if (nderiv >= 3) then
               c6_3_0 = a**2
               c6_3_1 = t**3
               c6_3_2 = 8*a**3
               c6_3_3 = t**2
               c6_3_4 = t**4
               c6_3_5 = t**5
               c6_3_6 = t**6
               c6_3_7 = t**7
               p(0, 3) = 12*c6_3_0*t - c6_3_1*c6_3_2
               p(1, 3) = -6*a + 24*c6_3_0*c6_3_3 - c6_3_2*c6_3_4
               p(2, 3) = -24*a*t + 36*c6_3_0*c6_3_1 - c6_3_2*c6_3_5
               p(3, 3) = -54*a*c6_3_3 + 48*c6_3_0*c6_3_4 - c6_3_2*c6_3_6 + 6
               p(4, 3) = -96*a*c6_3_1 + 60*c6_3_0*c6_3_5 - c6_3_2*c6_3_7 + 24*t
               p(5, 3) = -150*a*c6_3_4 + 72*c6_3_0*c6_3_6 - c6_3_2*t**8 + 60*c6_3_3
               p(6, 3) = -216*a*c6_3_5 + 84*c6_3_0*c6_3_7 + 120*c6_3_1 - c6_3_2*t**9
            end if
            if (nderiv >= 4) then
               c6_4_0 = a**2
               c6_4_1 = a**3
               c6_4_2 = t**2
               c6_4_3 = t**4
               c6_4_4 = 16*a**4
               c6_4_5 = t**3
               c6_4_6 = t**5
               c6_4_7 = t**6
               c6_4_8 = 120*t
               c6_4_9 = t**7
               c6_4_10 = t**8
               p(0, 4) = 12*c6_4_0 - 48*c6_4_1*c6_4_2 + c6_4_3*c6_4_4
               p(1, 4) = 60*c6_4_0*t - 80*c6_4_1*c6_4_5 + c6_4_4*c6_4_6
               p(2, 4) = -24*a + 156*c6_4_0*c6_4_2 - 112*c6_4_1*c6_4_3 + c6_4_4*c6_4_7
               p(3, 4) = -a*c6_4_8 + 300*c6_4_0*c6_4_5 - 144*c6_4_1*c6_4_6 + c6_4_4*c6_4_9
               p(4, 4) = -336*a*c6_4_2 + 492*c6_4_0*c6_4_3 - 176*c6_4_1*c6_4_7 + c6_4_10*c6_4_4 + 24
               p(5, 4) = -720*a*c6_4_5 + 732*c6_4_0*c6_4_6 - 208*c6_4_1*c6_4_9 + c6_4_4*t**9 + c6_4_8
               p(6, 4) = -1320*a*c6_4_3 + 1020*c6_4_0*c6_4_7 - 240*c6_4_1*c6_4_10 + 360*c6_4_2 + c6_4_4*t**10
            end if
         end block
      case (7)
         block
            real(wp) :: c7_1_0, c7_1_1, c7_1_2, c7_1_3, c7_1_4, c7_1_5, c7_2_0, c7_2_1
            real(wp) :: c7_2_2, c7_2_3, c7_2_4, c7_2_5, c7_2_6, c7_2_7, c7_3_0, c7_3_1
            real(wp) :: c7_3_2, c7_3_3, c7_3_4, c7_3_5, c7_3_6, c7_3_7, c7_3_8, c7_4_0
            real(wp) :: c7_4_1, c7_4_2, c7_4_3, c7_4_4, c7_4_5, c7_4_6, c7_4_7, c7_4_8
            real(wp) :: c7_4_9, c7_4_10, c7_4_11
            ! value (k=0)
            p(0, 0) = 1.0_wp
            p(1, 0) = t
            p(2, 0) = t**2
            p(3, 0) = t**3
            p(4, 0) = t**4
            p(5, 0) = t**5
            p(6, 0) = t**6
            p(7, 0) = t**7
            ! first derivative (k=1)
            c7_1_0 = t**2
            c7_1_1 = 2*a
            c7_1_2 = t**3
            c7_1_3 = t**4
            c7_1_4 = t**5
            c7_1_5 = t**6
            p(0, 1) = -2*a*t
            p(1, 1) = -c7_1_0*c7_1_1 + 1
            p(2, 1) = -c7_1_1*c7_1_2 + 2*t
            p(3, 1) = 3*c7_1_0 - c7_1_1*c7_1_3
            p(4, 1) = -c7_1_1*c7_1_4 + 4*c7_1_2
            p(5, 1) = -c7_1_1*c7_1_5 + 5*c7_1_3
            p(6, 1) = -c7_1_1*t**7 + 6*c7_1_4
            p(7, 1) = -c7_1_1*t**8 + 7*c7_1_5
            if (nderiv >= 2) then
               c7_2_0 = t**2
               c7_2_1 = 4*a**2
               c7_2_2 = 6*t
               c7_2_3 = t**3
               c7_2_4 = t**4
               c7_2_5 = t**5
               c7_2_6 = t**6
               c7_2_7 = t**7
               p(0, 2) = -2*a + c7_2_0*c7_2_1
               p(1, 2) = -a*c7_2_2 + c7_2_1*c7_2_3
               p(2, 2) = -10*a*c7_2_0 + c7_2_1*c7_2_4 + 2
               p(3, 2) = -14*a*c7_2_3 + c7_2_1*c7_2_5 + c7_2_2
               p(4, 2) = -18*a*c7_2_4 + 12*c7_2_0 + c7_2_1*c7_2_6
               p(5, 2) = -22*a*c7_2_5 + c7_2_1*c7_2_7 + 20*c7_2_3
               p(6, 2) = -26*a*c7_2_6 + c7_2_1*t**8 + 30*c7_2_4
               p(7, 2) = -30*a*c7_2_7 + c7_2_1*t**9 + 42*c7_2_5
            end if
            if (nderiv >= 3) then
               c7_3_0 = a**2
               c7_3_1 = t**3
               c7_3_2 = 8*a**3
               c7_3_3 = t**2
               c7_3_4 = t**4
               c7_3_5 = t**5
               c7_3_6 = t**6
               c7_3_7 = t**7
               c7_3_8 = t**8
               p(0, 3) = 12*c7_3_0*t - c7_3_1*c7_3_2
               p(1, 3) = -6*a + 24*c7_3_0*c7_3_3 - c7_3_2*c7_3_4
               p(2, 3) = -24*a*t + 36*c7_3_0*c7_3_1 - c7_3_2*c7_3_5
               p(3, 3) = -54*a*c7_3_3 + 48*c7_3_0*c7_3_4 - c7_3_2*c7_3_6 + 6
               p(4, 3) = -96*a*c7_3_1 + 60*c7_3_0*c7_3_5 - c7_3_2*c7_3_7 + 24*t
               p(5, 3) = -150*a*c7_3_4 + 72*c7_3_0*c7_3_6 - c7_3_2*c7_3_8 + 60*c7_3_3
               p(6, 3) = -216*a*c7_3_5 + 84*c7_3_0*c7_3_7 + 120*c7_3_1 - c7_3_2*t**9
               p(7, 3) = -294*a*c7_3_6 + 96*c7_3_0*c7_3_8 - c7_3_2*t**10 + 210*c7_3_4
            end if
            if (nderiv >= 4) then
               c7_4_0 = a**2
               c7_4_1 = a**3
               c7_4_2 = t**2
               c7_4_3 = t**4
               c7_4_4 = 16*a**4
               c7_4_5 = t**3
               c7_4_6 = t**5
               c7_4_7 = t**6
               c7_4_8 = 120*t
               c7_4_9 = t**7
               c7_4_10 = t**8
               c7_4_11 = t**9
               p(0, 4) = 12*c7_4_0 - 48*c7_4_1*c7_4_2 + c7_4_3*c7_4_4
               p(1, 4) = 60*c7_4_0*t - 80*c7_4_1*c7_4_5 + c7_4_4*c7_4_6
               p(2, 4) = -24*a + 156*c7_4_0*c7_4_2 - 112*c7_4_1*c7_4_3 + c7_4_4*c7_4_7
               p(3, 4) = -a*c7_4_8 + 300*c7_4_0*c7_4_5 - 144*c7_4_1*c7_4_6 + c7_4_4*c7_4_9
               p(4, 4) = -336*a*c7_4_2 + 492*c7_4_0*c7_4_3 - 176*c7_4_1*c7_4_7 + c7_4_10*c7_4_4 + 24
               p(5, 4) = -720*a*c7_4_5 + 732*c7_4_0*c7_4_6 - 208*c7_4_1*c7_4_9 + c7_4_11*c7_4_4 + c7_4_8
               p(6, 4) = -1320*a*c7_4_3 + 1020*c7_4_0*c7_4_7 - 240*c7_4_1*c7_4_10 + 360*c7_4_2 + c7_4_4*t**10
               p(7, 4) = -2184*a*c7_4_6 + 1356*c7_4_0*c7_4_9 - 272*c7_4_1*c7_4_11 + c7_4_4*t**11 + 840*c7_4_5
            end if
         end block
      case default
         block
            real(wp) :: c8_1_0, c8_1_1, c8_1_2, c8_1_3, c8_1_4, c8_1_5, c8_1_6, c8_2_0
            real(wp) :: c8_2_1, c8_2_2, c8_2_3, c8_2_4, c8_2_5, c8_2_6, c8_2_7, c8_2_8
            real(wp) :: c8_3_0, c8_3_1, c8_3_2, c8_3_3, c8_3_4, c8_3_5, c8_3_6, c8_3_7
            real(wp) :: c8_3_8, c8_3_9, c8_4_0, c8_4_1, c8_4_2, c8_4_3, c8_4_4, c8_4_5
            real(wp) :: c8_4_6, c8_4_7, c8_4_8, c8_4_9, c8_4_10, c8_4_11, c8_4_12
            ! value (k=0)
            p(0, 0) = 1.0_wp
            p(1, 0) = t
            p(2, 0) = t**2
            p(3, 0) = t**3
            p(4, 0) = t**4
            p(5, 0) = t**5
            p(6, 0) = t**6
            p(7, 0) = t**7
            p(8, 0) = t**8
            ! first derivative (k=1)
            c8_1_0 = t**2
            c8_1_1 = 2*a
            c8_1_2 = t**3
            c8_1_3 = t**4
            c8_1_4 = t**5
            c8_1_5 = t**6
            c8_1_6 = t**7
            p(0, 1) = -2*a*t
            p(1, 1) = -c8_1_0*c8_1_1 + 1
            p(2, 1) = -c8_1_1*c8_1_2 + 2*t
            p(3, 1) = 3*c8_1_0 - c8_1_1*c8_1_3
            p(4, 1) = -c8_1_1*c8_1_4 + 4*c8_1_2
            p(5, 1) = -c8_1_1*c8_1_5 + 5*c8_1_3
            p(6, 1) = -c8_1_1*c8_1_6 + 6*c8_1_4
            p(7, 1) = -c8_1_1*t**8 + 7*c8_1_5
            p(8, 1) = -c8_1_1*t**9 + 8*c8_1_6
            if (nderiv >= 2) then
               c8_2_0 = t**2
               c8_2_1 = 4*a**2
               c8_2_2 = 6*t
               c8_2_3 = t**3
               c8_2_4 = t**4
               c8_2_5 = t**5
               c8_2_6 = t**6
               c8_2_7 = t**7
               c8_2_8 = t**8
               p(0, 2) = -2*a + c8_2_0*c8_2_1
               p(1, 2) = -a*c8_2_2 + c8_2_1*c8_2_3
               p(2, 2) = -10*a*c8_2_0 + c8_2_1*c8_2_4 + 2
               p(3, 2) = -14*a*c8_2_3 + c8_2_1*c8_2_5 + c8_2_2
               p(4, 2) = -18*a*c8_2_4 + 12*c8_2_0 + c8_2_1*c8_2_6
               p(5, 2) = -22*a*c8_2_5 + c8_2_1*c8_2_7 + 20*c8_2_3
               p(6, 2) = -26*a*c8_2_6 + c8_2_1*c8_2_8 + 30*c8_2_4
               p(7, 2) = -30*a*c8_2_7 + c8_2_1*t**9 + 42*c8_2_5
               p(8, 2) = -34*a*c8_2_8 + c8_2_1*t**10 + 56*c8_2_6
            end if
            if (nderiv >= 3) then
               c8_3_0 = a**2
               c8_3_1 = t**3
               c8_3_2 = 8*a**3
               c8_3_3 = t**2
               c8_3_4 = t**4
               c8_3_5 = t**5
               c8_3_6 = t**6
               c8_3_7 = t**7
               c8_3_8 = t**8
               c8_3_9 = t**9
               p(0, 3) = 12*c8_3_0*t - c8_3_1*c8_3_2
               p(1, 3) = -6*a + 24*c8_3_0*c8_3_3 - c8_3_2*c8_3_4
               p(2, 3) = -24*a*t + 36*c8_3_0*c8_3_1 - c8_3_2*c8_3_5
               p(3, 3) = -54*a*c8_3_3 + 48*c8_3_0*c8_3_4 - c8_3_2*c8_3_6 + 6
               p(4, 3) = -96*a*c8_3_1 + 60*c8_3_0*c8_3_5 - c8_3_2*c8_3_7 + 24*t
               p(5, 3) = -150*a*c8_3_4 + 72*c8_3_0*c8_3_6 - c8_3_2*c8_3_8 + 60*c8_3_3
               p(6, 3) = -216*a*c8_3_5 + 84*c8_3_0*c8_3_7 + 120*c8_3_1 - c8_3_2*c8_3_9
               p(7, 3) = -294*a*c8_3_6 + 96*c8_3_0*c8_3_8 - c8_3_2*t**10 + 210*c8_3_4
               p(8, 3) = -384*a*c8_3_7 + 108*c8_3_0*c8_3_9 - c8_3_2*t**11 + 336*c8_3_5
            end if
            if (nderiv >= 4) then
               c8_4_0 = a**2
               c8_4_1 = a**3
               c8_4_2 = t**2
               c8_4_3 = t**4
               c8_4_4 = 16*a**4
               c8_4_5 = t**3
               c8_4_6 = t**5
               c8_4_7 = t**6
               c8_4_8 = 120*t
               c8_4_9 = t**7
               c8_4_10 = t**8
               c8_4_11 = t**9
               c8_4_12 = t**10
               p(0, 4) = 12*c8_4_0 - 48*c8_4_1*c8_4_2 + c8_4_3*c8_4_4
               p(1, 4) = 60*c8_4_0*t - 80*c8_4_1*c8_4_5 + c8_4_4*c8_4_6
               p(2, 4) = -24*a + 156*c8_4_0*c8_4_2 - 112*c8_4_1*c8_4_3 + c8_4_4*c8_4_7
               p(3, 4) = -a*c8_4_8 + 300*c8_4_0*c8_4_5 - 144*c8_4_1*c8_4_6 + c8_4_4*c8_4_9
               p(4, 4) = -336*a*c8_4_2 + 492*c8_4_0*c8_4_3 - 176*c8_4_1*c8_4_7 + c8_4_10*c8_4_4 + 24
               p(5, 4) = -720*a*c8_4_5 + 732*c8_4_0*c8_4_6 - 208*c8_4_1*c8_4_9 + c8_4_11*c8_4_4 + c8_4_8
               p(6, 4) = -1320*a*c8_4_3 + 1020*c8_4_0*c8_4_7 - 240*c8_4_1*c8_4_10 + c8_4_12*c8_4_4 + 360*c8_4_2
               p(7, 4) = -2184*a*c8_4_6 + 1356*c8_4_0*c8_4_9 - 272*c8_4_1*c8_4_11 + c8_4_4*t**11 + 840*c8_4_5
               p(8, 4) = -3360*a*c8_4_7 + 1740*c8_4_0*c8_4_10 - 304*c8_4_1*c8_4_12 + 1680*c8_4_3 + c8_4_4*t**12
            end if
         end block
      end select

   end subroutine moist_iso_gto_poly1d

end module moist_cavity_drop_lsf_isodensity_gto
