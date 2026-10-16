!> Cartesian-monomial Gaussian basis + density evaluator for the internal
!> isodensity DROP level set function
!>
!> The internal isodensity LSF needs the electron density rho(r) and its spatial
!> derivatives, up to fourth order, at arbitrary points; rather than reproduce a
!> host code's spherical-harmonic / normalization conventions, moist evaluates
!> only bare cartesian monomial Gaussians
!>
!>    g_c(r) = (x-Rx)^lx (y-Ry)^ly (z-Rz)^lz * sum_p coeff_p exp(-a_p |r-R|^2),
!>
!> and consumes a density matrix already transformed into this
!> cartesian-monomial basis (D_cart), the physical density
!>
!>    rho(r) = sum_{c,c'} D_cart(c,c') g_c(r) g_c'(r)
!>
!> being representation-invariant
!>
!> - the host builds the fixed transform once from its own basis metadata and
!>   passes D_cart each SCF step
!> - the cartesian component ordering used here is reported through the API, so
!>   the host can match it exactly
!>
!> A single cartesian primitive factorizes as g = X(dx) Y(dy) Z(dz), so every
!> spatial derivative is a product of the exp-stripped 1D factor derivatives
!> tabulated by [[moist_iso_gto_poly1d]] (sympy-generated)
!>
!> The density is contracted over internal basis functions, chosen at
!> [[gto_init]] by [[moist_iso_basis]]:
!>
!> - cartesian: the monomials themselves
!> - spherical: per shell with l >= 2, the 2l+1 harmonic polynomials h = U^T g,
!>   each the harmonic extension of a leading monomial (lx <= 1)
!> - spherical is exact whenever D_cart = U D_sph U^T; D_sph is then the leading
!>   sub-block of D_cart, rho = h^T D_sph h, and [[gto_set_density]] rejects any
!>   other density
!> - functions are ordered atom by atom, most diffuse shell first, so the shells
!>   of an atom surviving radial screening form one contiguous run
module moist_cavity_drop_lsf_isodensity_gto
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_math_blas, only: gemm
   implicit none(type, external)
   private

   public :: moist_iso_gto_type
   public :: moist_iso_gto_work_type
   public :: moist_iso_gto_ncart
   public :: moist_iso_gto_nslot
   public :: moist_iso_gto_nrhs
   public :: moist_iso_basis_type
   public :: moist_iso_basis

   !> Highest supported angular momentum (s=0 .. l=8)
   !>
   !> - emitted, with the [[moist_iso_gto_poly1d]] and [[moist_iso_gto_fold]]
   !>   routines and the harmonic tables below, by
   !>   `moist_dev/tools/gen_isodensity_gto.py` into the marked regions
   !>
   integer, parameter :: moist_iso_gto_lmax = 8

   !> Harmonic-extension coefficients of moist's spherical basis, l = 2..lmax
   !>
   !> - entries harm_off(l)+1 .. harm_off(l+1) hold the non-leading terms of the
   !>   2l+1 functions of angular momentum l
   !> - entry e adds harm_coef(e) times cartesian component harm_cart(e) to
   !>   spherical function harm_sph(e), both 1-based within the shell
   integer, parameter :: harm_off(2:9) = [0, 2, 8, 21, 44, 81, 136, 214]
   !> Non-leading cartesian component of each entry
   integer, parameter :: harm_cart(214) = [ &
                         1, 1, 1, 1, 2, 3, 2, 3, 2, 3, 2, 3, 1, 4, 5, 1, 4, 6, 5, 1, &
                         6, 1, 4, 5, 1, 4, 6, 5, 1, 6, 2, 7, 3, 8, 2, 7, 9, 3, 8, 10, &
                         2, 9, 3, 10, 2, 7, 3, 8, 2, 7, 9, 3, 8, 10, 2, 9, 3, 10, 1, 4, &
                         11, 5, 12, 1, 4, 6, 11, 13, 5, 12, 14, 1, 4, 6, 13, 15, 5, 14, 1, 6, &
                         15, 1, 4, 11, 5, 12, 1, 4, 6, 11, 13, 5, 12, 14, 1, 4, 6, 13, 15, 5, &
                         14, 1, 6, 15, 2, 7, 16, 3, 8, 17, 2, 7, 9, 16, 18, 3, 8, 10, 17, 19, &
                         2, 7, 9, 18, 20, 3, 8, 10, 19, 21, 2, 9, 20, 3, 10, 21, 2, 7, 16, 3, &
                         8, 17, 2, 7, 9, 16, 18, 3, 8, 10, 17, 19, 2, 7, 9, 18, 20, 3, 8, 10, &
                         19, 21, 2, 9, 20, 3, 10, 21, 1, 4, 11, 22, 5, 12, 23, 1, 4, 6, 11, 13, &
                         22, 24, 5, 12, 14, 23, 25, 1, 4, 6, 11, 13, 15, 24, 26, 5, 12, 14, 25, 27, &
                         1, 4, 6, 13, 15, 26, 28, 5, 14, 27, 1, 6, 15, 28]
   !> Spherical function of each entry
   integer, parameter :: harm_sph(214) = [ &
                         3, 5, 1, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 5, 6, 7, 7, 7, 8, 9, &
                         9, 1, 1, 2, 3, 3, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9, 9, &
                         10, 10, 11, 11, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 7, &
                         7, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 13, 13, &
                         13, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, &
                         6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, &
                         12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 15, 15, 15, 1, 1, 1, 2, &
                         2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, &
                         6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 11, &
                         11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, &
                         15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 17, 17, 17, 17]
   !> Coefficient U(c, m) of each entry
   real(wp), parameter :: harm_coef(214) = [ &
                          -1.0_wp, -1.0_wp, -0.3333333333333333_wp, -0.3333333333333333_wp, &
                          -3.0_wp, -1.0_wp, -1.0_wp, -3.0_wp, &
                          -1.0_wp, -0.3333333333333333_wp, -0.3333333333333333_wp, -1.0_wp, &
                          1.0_wp, -6.0_wp, -3.0_wp, 0.3333333333333333_wp, &
                          -1.0_wp, -1.0_wp, -3.0_wp, 1.0_wp, &
                          -6.0_wp, 0.2_wp, -2.0_wp, -1.0_wp, &
                          0.06666666666666667_wp, -0.3333333333333333_wp, -0.3333333333333333_wp, -1.0_wp, &
                          0.2_wp, -2.0_wp, 5.0_wp, -10.0_wp, &
                          1.0_wp, -6.0_wp, 1.0_wp, -1.0_wp, &
                          -3.0_wp, 1.0_wp, -3.0_wp, -1.0_wp, &
                          1.0_wp, -6.0_wp, 5.0_wp, -10.0_wp, &
                          1.0_wp, -3.3333333333333335_wp, 0.2_wp, -2.0_wp, &
                          0.2_wp, -0.3333333333333333_wp, -1.0_wp, 0.2_wp, &
                          -1.0_wp, -0.3333333333333333_wp, 0.2_wp, -2.0_wp, &
                          1.0_wp, -3.3333333333333335_wp, -1.0_wp, 15.0_wp, &
                          -15.0_wp, 5.0_wp, -10.0_wp, -0.2_wp, &
                          2.0_wp, 1.0_wp, -1.0_wp, -6.0_wp, &
                          3.0_wp, -3.0_wp, -3.0_wp, -0.2_wp, &
                          1.0_wp, 2.0_wp, -6.0_wp, -1.0_wp, &
                          5.0_wp, -10.0_wp, -1.0_wp, 15.0_wp, &
                          -15.0_wp, -0.14285714285714285_wp, 3.0_wp, -5.0_wp, &
                          1.0_wp, -3.3333333333333335_wp, -0.02857142857142857_wp, 0.4_wp, &
                          0.2_wp, -0.3333333333333333_wp, -2.0_wp, 0.6_wp, &
                          -1.0_wp, -1.0_wp, -0.02857142857142857_wp, 0.2_wp, &
                          0.4_wp, -2.0_wp, -0.3333333333333333_wp, 1.0_wp, &
                          -3.3333333333333335_wp, -0.14285714285714285_wp, 3.0_wp, -5.0_wp, &
                          -7.0_wp, 35.0_wp, -21.0_wp, -1.0_wp, &
                          15.0_wp, -15.0_wp, -1.0_wp, 3.3333333333333335_wp, &
                          5.0_wp, -1.0_wp, -10.0_wp, -0.6_wp, &
                          6.0_wp, 1.0_wp, -3.0_wp, -6.0_wp, &
                          -0.6_wp, 1.0_wp, 6.0_wp, -6.0_wp, &
                          -3.0_wp, -1.0_wp, 5.0_wp, 3.3333333333333335_wp, &
                          -10.0_wp, -1.0_wp, -1.0_wp, 15.0_wp, &
                          -15.0_wp, -7.0_wp, 35.0_wp, -21.0_wp, &
                          -1.0_wp, 7.0_wp, -7.0_wp, -0.14285714285714285_wp, &
                          3.0_wp, -5.0_wp, -0.14285714285714285_wp, 0.6666666666666666_wp, &
                          1.0_wp, -0.3333333333333333_wp, -3.3333333333333335_wp, -0.08571428571428572_wp, &
                          1.2_wp, 0.2_wp, -1.0_wp, -2.0_wp, &
                          -0.08571428571428572_wp, 0.2_wp, 1.2_wp, -2.0_wp, &
                          -1.0_wp, -0.14285714285714285_wp, 1.0_wp, 0.6666666666666666_wp, &
                          -3.3333333333333335_wp, -0.3333333333333333_wp, -0.14285714285714285_wp, 3.0_wp, &
                          -5.0_wp, -1.0_wp, 7.0_wp, -7.0_wp, &
                          1.0_wp, -28.0_wp, 70.0_wp, -28.0_wp, &
                          -7.0_wp, 35.0_wp, -21.0_wp, 0.14285714285714285_wp, &
                          -3.0_wp, -1.0_wp, 5.0_wp, 15.0_wp, &
                          -1.0_wp, -15.0_wp, -3.0_wp, 10.0_wp, &
                          5.0_wp, -3.0_wp, -10.0_wp, 0.08571428571428572_wp, &
                          -1.2_wp, -1.2_wp, 1.0_wp, 12.0_wp, &
                          1.0_wp, -6.0_wp, -6.0_wp, -3.0_wp, &
                          5.0_wp, 10.0_wp, -10.0_wp, -3.0_wp, &
                          0.14285714285714285_wp, -1.0_wp, -3.0_wp, 15.0_wp, &
                          5.0_wp, -15.0_wp, -1.0_wp, -7.0_wp, &
                          35.0_wp, -21.0_wp, 1.0_wp, -28.0_wp, &
                          70.0_wp, -28.0_wp]

   !> Enumerator for the density representation of [[moist_iso_gto_type]]
   type :: moist_iso_basis_type
      !> Cartesian monomials, valid for any host basis
      integer :: cartesian = 1
      !> Harmonic polynomials, valid for a host with spherical basis functions
      integer :: spherical = 2
   end type moist_iso_basis_type

   !> Global instance of the density-representation enum
   type(moist_iso_basis_type), parameter :: moist_iso_basis = moist_iso_basis_type()

   !> Width of the packed density panel of the screened multi-column contraction
   integer, parameter :: iso_panel_width = 256

   !> Largest out-of-subspace residual of a spherical density, relative to max |D|
   real(wp), parameter :: iso_harmonic_tol = 1.0e-10_wp

   !> Symmetric Hessian slot index (4..9) for spatial axis pair (i, j)
   integer, parameter :: hess_slot(3, 3) = reshape([ &
                                                   4, 5, 6, &
                                                   5, 7, 8, &
                                                   6, 8, 9], [3, 3])

   !> Symmetric third-derivative slot index (10..19) for axis triple (i, j, k)
   !> Built column-major as third_slot(k, j, i); symmetric under any permutation
   integer, parameter :: third_slot(3, 3, 3) = reshape([ &
                                                       10, 11, 12, 11, 13, 14, 12, 14, 15, &
                                                       11, 13, 14, 13, 16, 17, 14, 17, 18, &
                                                       12, 14, 15, 14, 17, 18, 15, 18, 19], [3, 3, 3])

   !> Cumulative number of packed derivative slots through each order, 0..4
   integer, parameter :: moist_iso_gto_nslot(0:4) = [1, 4, 10, 20, 35]

   !> Number of density-weighted AO-derivative columns per derivative order, 0..4
   !>
   !> - value and gradient need D phi; Hessian and third order add D phi_i;
   !>   fourth order adds D phi_ij
   integer, parameter :: moist_iso_gto_nrhs(0:4) = [1, 1, 4, 4, 10]

   !> Cartesian-monomial Gaussian basis together with a density matrix expressed
   !> over the internal basis functions
   !>
   !> - read-only after `init` / `set_density`, hence safe to share across the
   !>   cavity's per-thread LSF clones
   type :: moist_iso_gto_type
      !> Number of shells
      integer :: nshell = 0
      !> Total number of cartesian components (sum of ncart(l) over shells)
      integer :: ncart = 0
      !> Highest angular momentum present
      integer :: lmax = 0
      !> Density representation, a [[moist_iso_basis]] value
      integer :: basis = moist_iso_basis%cartesian
      !> Number of internal basis functions the density is contracted over
      integer :: nfun = 0
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
      !> Per-shell offset into the internal function order, size nshell; the
      !> functions of shell s are sh_foff(s)+1 .. sh_foff(s)+nfun(l)
      integer, allocatable :: sh_foff(:)
      !> Host cartesian component behind each internal function, size nfun
      !>
      !> - the monomial itself, or the leading monomial of a harmonic polynomial
      integer, allocatable :: fun_src(:)
      !> Density matrix over the internal functions, symmetric (nfun, nfun)
      real(wp), allocatable :: dens(:, :)
      !> Highest owner-atom index present in the basis (= size of the CSR below)
      integer :: natom_grid = 0
      !> Atom -> shell CSR offsets, size natom_grid+1; atom A owns shells
      !> `atom_shell(atom_soff(A)+1 : atom_soff(A+1))`
      !>
      !> - turns a candidate-atom list from the cavity cell grid into a shell list
      !>   without a scan
      integer, allocatable :: atom_soff(:)
      !> Flat-packed shell indices grouped by owner atom, most diffuse first
      !> within an atom; this is the internal function order, size nshell
      integer, allocatable :: atom_shell(:)
      !> Per-shell squared radial cutoff, beyond which the shell contributes
      !> nothing to the value or to derivatives through the order requested from
      !> [[gto_build_screening]], i.e. for `|r - center|^2 > sh_rcut2`
      !>
      !> - set by [[gto_build_screening]]
      !> - huge, meaning no screening, until then and when disabled
      real(wp), allocatable :: sh_rcut2(:)
      !> Largest shell reach sqrt(max(sh_rcut2)); the radius the cavity cell grid
      !> must span so no contributing atom is missed as a candidate
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

   !> Scratch of [[gto_eval]], one per evaluating thread
   !>
   !> - sized on first use for the basis and derivative order, then reused, so
   !>   evaluation performs no heap allocation
   !> - rows of `phi` are compact: row q holds the q-th active function in visiting
   !>   order, the active runs listing which internal functions those are
   !> - `dphi`, `panel` and `rows` are storage only, viewed per point as contiguous
   !>   blocks of leading dimension nact (or the panel width), so every gemm
   !>   operand is contiguous and no array temporary is ever created
   type :: moist_iso_gto_work_type
      !> AO-derivative table of the active functions, (nfun, 0:nslot-1)
      real(wp), allocatable :: phi(:, :)
      !> Storage of the density-weighted AO derivatives D phi(:, 0:nrhs-1),
      !> (nfun, 0:nrhs-1), viewed as (nact, 0:nrhs-1)
      real(wp), allocatable :: dphi(:, :)
      !> Storage of the packed density panel of the screened multi-column
      !> contraction, flat (nfun*iso_panel_width)
      real(wp), allocatable :: panel(:)
      !> Storage of the phi rows matching one panel, flat (iso_panel_width*10)
      real(wp), allocatable :: rows(:)
      !> Internal offset of each active run, size nshell
      integer, allocatable :: run_off(:)
      !> Function count of each active run, size nshell
      integer, allocatable :: run_len(:)
   end type moist_iso_gto_work_type

contains

   !> Number of cartesian components for angular momentum l
   !>
   !> @param[in] l  angular momentum
   !> @returns      (l+1)(l+2)/2
   elemental function moist_iso_gto_ncart(l) result(n)
      !> Angular momentum
      integer, intent(in) :: l
      integer :: n

      n = (l + 1)*(l + 2)/2
   end function moist_iso_gto_ncart

   !> Number of internal functions of a shell
   !>
   !> @param[in] basis  density representation, a [[moist_iso_basis]] value
   !> @param[in] l      angular momentum
   !> @returns          2l+1 for spherical, (l+1)(l+2)/2 for cartesian
   elemental function shell_nfun(basis, l) result(n)
      !> Density representation
      integer, intent(in) :: basis
      !> Angular momentum
      integer, intent(in) :: l
      !> Function count
      integer :: n

      if (basis == moist_iso_basis%spherical) then
         n = 2*l + 1
      else
         n = moist_iso_gto_ncart(l)
      end if
   end function shell_nfun

   !> Symmetric fourth-derivative slot index (20..34) for a spatial axis 4-tuple
   !>
   !> The fourth-derivative table has the 15 unique cartesian components with
   !> kx+ky+kz = 4, laid out like the lower orders in descending-kx-then-ky order
   !> at phi columns 20..34
   !>
   !> - maps any axis quadruple (i, j, k, l), each in 1..3, to its column by
   !>   counting the per-axis derivative orders
   !> - cheaper and less error-prone than a 3x3x3x3 lookup table, and only used on
   !>   the non-hot fourth-derivative path
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
      !> 0-based position of (kx=nx, ky=ny) in the descending-kx-then-ky order
      slot = 4 - nx - ny
      do kx = nx + 1, 4
         slot = slot + (5 - kx)
      end do
      slot = 20 + slot
   end function fourth_slot

   !> Configure the basis, derive the cartesian-component layout and the internal
   !> function order
   !>
   !> @param[inout] self      Basis instance
   !> @param[in]    sh_atom   Per-shell owner atom index (1-based), size nshell
   !> @param[in]    sh_l      Per-shell angular momentum, size nshell
   !> @param[in]    sh_nprim  Per-shell primitive count, size nshell
   !> @param[in]    exps      Primitive exponents (concatenated per shell)
   !> @param[in]    coeffs    Primitive contraction coefficients (host-normalized)
   !> @param[out]   error     Set on invalid input
   !> @param[in]    basis     Density representation, a [[moist_iso_basis]] value;
   !>                         cartesian when absent
   subroutine gto_init(self, sh_atom, sh_l, sh_nprim, exps, coeffs, error, basis)
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
      !> Density representation
      integer, intent(in), optional :: basis

      integer :: s, l, c, lx, ly, lz, nprim_tot, idx, nf, first, selected

      selected = moist_iso_basis%cartesian
      if (present(basis)) selected = basis
      if (selected /= moist_iso_basis%cartesian .and. selected /= moist_iso_basis%spherical) then
         call fatal_error(error, "isodensity basis representation must be cartesian or spherical")
         return
      end if

      self%nshell = size(sh_l)
      if (size(sh_atom) /= self%nshell .or. size(sh_nprim) /= self%nshell) then
         call fatal_error(error, "isodensity basis shell arrays have inconsistent length")
         return
      end if
      if (self%nshell < 1) then
         call fatal_error(error, "isodensity basis has no shells")
         return
      end if
      if (any(sh_atom < 1)) then
         call fatal_error(error, "isodensity basis shell atom index must be positive")
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
      ! previous basis only after the new input has passed validation
      if (allocated(self%sh_poff)) deallocate (self%sh_poff)
      if (allocated(self%sh_coff)) deallocate (self%sh_coff)
      if (allocated(self%center)) deallocate (self%center)
      if (allocated(self%comp_l)) deallocate (self%comp_l)
      if (allocated(self%sh_foff)) deallocate (self%sh_foff)
      if (allocated(self%fun_src)) deallocate (self%fun_src)
      if (allocated(self%dens)) deallocate (self%dens)
      if (allocated(self%atom_soff)) deallocate (self%atom_soff)
      if (allocated(self%atom_shell)) deallocate (self%atom_shell)
      if (allocated(self%sh_rcut2)) deallocate (self%sh_rcut2)

      self%basis = selected
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
      !> ly from (l-lx) down to 0, lz = l-lx-ly, matching the common CINT order
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

      !> Atom -> shell CSR, a counting sort on the owner atom
      !>
      !> - candidate atom lists from the cavity cell grid index straight into it,
      !>   so the hot projection path visits only the shells on nearby atoms
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
      do s = 1, self%natom_grid
         call sort_diffuse_first(self%atom_shell(self%atom_soff(s) + 1:self%atom_soff(s + 1)))
      end do

      !> Internal function order: atom by atom along atom_shell
      allocate (self%sh_foff(self%nshell))
      nf = 0
      do idx = 1, self%nshell
         s = self%atom_shell(idx)
         self%sh_foff(s) = nf
         nf = nf + shell_nfun(self%basis, self%sh_l(s))
      end do
      self%nfun = nf

      !> Spherical functions take the leading monomials, the last 2l+1 of a shell
      allocate (self%fun_src(self%nfun))
      do s = 1, self%nshell
         l = self%sh_l(s)
         first = self%sh_coff(s) + moist_iso_gto_ncart(l) - shell_nfun(self%basis, l)
         do c = 1, shell_nfun(self%basis, l)
            self%fun_src(self%sh_foff(s) + c) = first + c
         end do
      end do

      !> Default to "no screening" (huge cutoff) until build_screening runs, so a
      !> direct user (no cavity / zero threshold) evaluates every shell exactly
      allocate (self%sh_rcut2(self%nshell), source=huge(1.0_wp))
      self%max_rcut = 0.0_wp

   contains

      !> Stable insertion sort of one atom's shells by smallest exponent
      !>
      !> - the most diffuse shell reaches furthest, so the shells inside their
      !>   radial cutoff at any point are, up to coefficient effects, a prefix
      !>
      !> @param[inout] shells  shell indices of one atom
      subroutine sort_diffuse_first(shells)
         !> Shell indices of one atom
         integer, intent(inout) :: shells(:)

         integer :: i, j, moving
         real(wp) :: key

         do i = 2, size(shells)
            moving = shells(i)
            key = smallest_exponent(moving)
            j = i - 1
            do while (j >= 1)
               if (smallest_exponent(shells(j)) <= key) exit
               shells(j + 1) = shells(j)
               j = j - 1
            end do
            shells(j + 1) = moving
         end do
      end subroutine sort_diffuse_first

      !> Smallest primitive exponent of a shell
      !>
      !> @param[in] sh  shell index
      !> @returns       smallest exponent
      pure function smallest_exponent(sh) result(amin)
         !> Shell index
         integer, intent(in) :: sh
         !> Smallest exponent
         real(wp) :: amin

         amin = minval(self%exps(self%sh_poff(sh):self%sh_poff(sh + 1) - 1))
      end function smallest_exponent

   end subroutine gto_init

   !> Refresh the shell centers from a molecular structure
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

   !> Install the cartesian-monomial density matrix
   !>
   !> Gathers the density over the internal functions; in the spherical
   !> representation first checks that it lies in the span of the spherical
   !> functions, see [[check_spherical_density]]
   !>
   !> - unchanged on error
   !>
   !> @param[inout] self   Basis instance
   !> @param[in]    dcart  Density matrix in the cartesian-monomial basis
   !> @param[out]   error  Set on a size mismatch, a non-spherical density in the
   !>                      spherical representation, or a failed allocation
   subroutine gto_set_density(self, dcart, error)
      !> Basis instance
      class(moist_iso_gto_type), intent(inout) :: self
      !> Cartesian-monomial density matrix (ncart, ncart)
      real(wp), intent(in) :: dcart(:, :)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer :: i, j, src, stat

      if (size(dcart, 1) /= self%ncart .or. size(dcart, 2) /= self%ncart) then
         call fatal_error(error, "isodensity density matrix has the wrong shape")
         return
      end if
      if (self%basis == moist_iso_basis%spherical) then
         call check_spherical_density(self, dcart, error)
         if (allocated(error)) return
      end if

      if (allocated(self%dens)) then
         if (size(self%dens, 1) /= self%nfun) deallocate (self%dens)
      end if
      if (.not. allocated(self%dens)) then
         allocate (self%dens(self%nfun, self%nfun), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "cannot allocate the isodensity density matrix")
            return
         end if
      end if
      do j = 1, self%nfun
         src = self%fun_src(j)
         do i = 1, self%nfun
            self%dens(i, j) = dcart(self%fun_src(i), src)
         end do
      end do
   end subroutine gto_set_density

   !> Check that a cartesian density lies in the span of the spherical functions
   !>
   !> Per shell with l >= 2, each column of D_cart must be U times its leading
   !> entries, hence D_cart = U D_sph U^T with D_sph the leading block
   !>
   !> - a host with cartesian basis functions carries weight on the r^2 l-2
   !>   contaminants there and fails
   !> - the largest deviation is measured relative to max |D_cart| and must stay
   !>   below [[iso_harmonic_tol]], far above host round-off
   !>
   !> @param[in]  self   Basis instance
   !> @param[in]  dcart  Density matrix in the cartesian-monomial basis
   !> @param[out] error  Set when the density leaves the spherical subspace
   subroutine check_spherical_density(self, dcart, error)
      !> Basis instance
      class(moist_iso_gto_type), intent(in) :: self
      !> Cartesian-monomial density matrix (ncart, ncart)
      real(wp), intent(in) :: dcart(:, :)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: pred((moist_iso_gto_lmax + 1)*(moist_iso_gto_lmax + 2)/2)
      real(wp) :: dmax, resid
      integer :: j, s, l, cb, nlead0, e, c
      character(len=16) :: buffer

      dmax = maxval(abs(dcart))
      resid = 0.0_wp
      do j = 1, self%ncart
         do s = 1, self%nshell
            l = self%sh_l(s)
            if (l < 2) cycle
            cb = self%sh_coff(s)
            ! Non-leading components come first in a shell, the 2l+1 leading last
            nlead0 = moist_iso_gto_ncart(l) - (2*l + 1)
            pred(1:nlead0) = 0.0_wp
            do e = harm_off(l) + 1, harm_off(l + 1)
               pred(harm_cart(e)) = pred(harm_cart(e)) &
                                    + harm_coef(e)*dcart(cb + nlead0 + harm_sph(e), j)
            end do
            do c = 1, nlead0
               resid = max(resid, abs(dcart(cb + c, j) - pred(c)))
            end do
         end do
      end do

      if (resid > iso_harmonic_tol*dmax) then
         write (buffer, "(es9.2)") resid/dmax
         call fatal_error(error, "isodensity density matrix has weight outside the spherical "// &
                          "subspace (relative residual "//trim(adjustl(buffer))//"); a host with "// &
                          "cartesian basis functions needs the cartesian representation")
      end if
   end subroutine check_spherical_density

   !> Whether a density matrix has been installed
   !>
   !> @param[in] self  Basis instance
   !> @returns         .true. once set_density has been called
   pure function gto_has_density(self) result(ok)
      !> Basis instance
      class(moist_iso_gto_type), intent(in) :: self
      logical :: ok

      ok = allocated(self%dens)
   end function gto_has_density

   !> Conservative squared amplitude bound for one primitive at radius `r`
   !>
   !> Bounds the magnitude of any cartesian component of `x^lx y^ly z^lz
   !> exp(-a r^2)`, and of its spatial derivatives up to order `k`, by
   !> `|c| (r+1)^l (2 a r + l + 1)^k exp(-a r^2)`
   !>
   !> - an over-estimate: the monomial is <= r^l <= (r+1)^l, and each
   !>   differentiation multiplies the bound by at most 2 a r + l
   !> - screening on it therefore never drops a shell still contributing above
   !>   the threshold
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

   !> Radius beyond which a whole shell's amplitude bound falls below `tol`
   !>
   !> Sums [[prim_amp_bound]] over the shell's primitives at `max_deriv` and
   !> bisects for the single outer crossing of
   !> `tol`; the bound is large at r=0 and decays monotonically past its peak,
   !> so the crossing is unique.  `tol <= 0` disables screening (returns a huge
   !> reach)
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

      !> Expand hi until the bound is below tol (or the cap is hit)
      hi = 1.0_wp
      do
         if (shell_amp(hi) <= tol .or. hi >= r_cap) exit
         hi = 2.0_wp*hi
      end do
      if (hi >= r_cap) then
         rcut = r_cap
         return
      end if

      !> Bisect [lo, hi] with amp(lo) > tol >= amp(hi) for the outer crossing
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
   !> `threshold` is an absolute amplitude cutoff on the (bounded) AO product contribution;
   !> `threshold <= 0` disables screening
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

   !> Evaluate the electron density and its spatial derivatives at one point
   !>
   !> The derivative order is set by *which* outputs the caller passes: only the
   !> orders whose output argument is present are computed, so the hot
   !> value+gradient path never pays for the higher orders
   !>
   !> Three stages, each proportional to the *nearby* basis rather than the whole
   !> molecule:
   !>
   !> - assembly: with `cand_atoms` present, only the shells owned by those atoms
   !>   (the cavity cell grid's candidate list) are visited, and every shell beyond
   !>   its radial cutoff [[gto_type:sh_rcut2]] is skipped; the survivors fill
   !>   compact rows of `work%phi`, folded to the internal functions
   !> - contraction: `work%dphi = D phi` for every column the requested order needs,
   !>   in one pass over the active block of the density, see [[gto_contract]]
   !> - derivatives: all from the small gram matrix phi^T D phi, the unique
   !>   components only, then symmetrized
   !> - with screening disabled (huge cutoffs) every function is active and the
   !>   result is the exact dense contraction
   !>
   !> @param[in]    self       Basis instance (read-only)
   !> @param[in]    point      Evaluation point in Bohr
   !> @param[inout] work       Per-thread scratch, sized on first use
   !> @param[out]   rho        Electron density
   !> @param[out]   drho       Density gradient (3)
   !> @param[out]   d2rho      Optional density Hessian (3, 3)
   !> @param[out]   d3rho      Optional density third derivative (3, 3, 3)
   !> @param[in]    cand_atoms Optional candidate owner-atom list; when present,
   !>                          only shells on these atoms are considered
   !> @param[out]   d4rho      Optional density fourth derivative (3, 3, 3, 3)
   subroutine gto_eval(self, point, work, rho, drho, d2rho, d3rho, cand_atoms, d4rho)
      !> Basis instance
      class(moist_iso_gto_type), intent(in) :: self
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Per-thread scratch
      type(moist_iso_gto_work_type), intent(inout), target :: work
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

      !> Density-weighted AO derivatives of the active functions, contiguous
      real(wp), pointer, contiguous :: dphi(:, :)
      real(wp) :: gram(0:34, 0:9)
      integer :: ka, aidx, idx, nslot, nrhs, nact, nrun, ndloc

      ! The requested order is exactly the highest order the caller asked back
      ndloc = 1
      if (present(d2rho)) ndloc = 2
      if (present(d3rho)) ndloc = 3
      if (present(d4rho)) ndloc = 4
      nslot = moist_iso_gto_nslot(ndloc)
      nrhs = moist_iso_gto_nrhs(ndloc)
      call size_work(self, work, nslot, nrhs)

      rho = 0.0_wp
      drho = 0.0_wp
      if (present(d2rho)) d2rho = 0.0_wp
      if (present(d3rho)) d3rho = 0.0_wp
      if (present(d4rho)) d4rho = 0.0_wp

      ! Assemble the surviving shells along the internal function order
      nact = 0
      nrun = 0
      if (present(cand_atoms)) then
         do ka = 1, size(cand_atoms)
            aidx = cand_atoms(ka)
            if (aidx < 1 .or. aidx > self%natom_grid) cycle
            do idx = self%atom_soff(aidx) + 1, self%atom_soff(aidx + 1)
               call add_shell(self%atom_shell(idx))
            end do
         end do
      else
         do idx = 1, self%nshell
            call add_shell(self%atom_shell(idx))
         end do
      end if
      if (nact == 0) return

      dphi(1:nact, 0:nrhs - 1) => work%dphi
      call gto_contract(self, work%phi, dphi, work%panel, work%rows, &
                        work%run_off(1:nrun), work%run_len(1:nrun))
      call gto_gram(work%phi, dphi, nslot, gram)

      rho = gram(0, 0)
      drho = 2.0_wp*gram(1:3, 0)
      if (present(d2rho)) call density_hessian(gram, d2rho)
      if (present(d3rho)) call density_third(gram, d3rho)
      if (present(d4rho)) call density_fourth(gram, d4rho)

   contains

      !> Assemble one shell into the next compact rows of `work%phi` and extend
      !> the active runs, skipping it beyond its radial cutoff
      !>
      !> - host-associated: reads point/ndloc/nslot/self, writes work/nact/nrun
      !> - in the spherical representation, shells with l >= 2 are assembled in a
      !>   local cartesian block and folded by [[moist_iso_gto_fold]]
      !>
      !> @param[in] sh  shell index
      subroutine add_shell(sh)
         !> Shell index
         integer, intent(in) :: sh

         real(wp) :: cart((moist_iso_gto_lmax + 1)*(moist_iso_gto_lmax + 2)/2, 0:34)
         real(wp) :: d(3), u
         integer :: l, nc, nf, off

         d = point - self%center(:, sh)
         u = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
         if (u > self%sh_rcut2(sh)) return

         l = self%sh_l(sh)
         nc = moist_iso_gto_ncart(l)
         nf = shell_nfun(self%basis, l)
         if (nf < nc) then
            cart(1:nc, 0:nslot - 1) = 0.0_wp
            call gto_assemble_shell(self, sh, d, u, ndloc, cart)
            call moist_iso_gto_fold(l, nslot, cart, work%phi(nact + 1:nact + nf, :))
         else
            work%phi(nact + 1:nact + nf, 0:nslot - 1) = 0.0_wp
            call gto_assemble_shell(self, sh, d, u, ndloc, work%phi(nact + 1:nact + nf, :))
         end if

         off = self%sh_foff(sh)
         nact = nact + nf
         if (nrun > 0) then
            if (work%run_off(nrun) + work%run_len(nrun) == off) then
               work%run_len(nrun) = work%run_len(nrun) + nf
               return
            end if
         end if
         nrun = nrun + 1
         work%run_off(nrun) = off
         work%run_len(nrun) = nf
      end subroutine add_shell

   end subroutine gto_eval

   !> Size the evaluation scratch for this basis and derivative order
   !>
   !> - reallocates only when the basis changed or a higher order is requested
   !>
   !> @param[in]    self   Basis instance
   !> @param[inout] work   Evaluation scratch
   !> @param[in]    nslot  number of AO-derivative slots
   !> @param[in]    nrhs   number of density-weighted columns
   subroutine size_work(self, work, nslot, nrhs)
      !> Basis instance
      class(moist_iso_gto_type), intent(in) :: self
      !> Evaluation scratch
      type(moist_iso_gto_work_type), intent(inout) :: work
      !> Number of AO-derivative slots
      integer, intent(in) :: nslot
      !> Number of density-weighted columns
      integer, intent(in) :: nrhs

      if (allocated(work%phi)) then
         if (size(work%phi, 1) /= self%nfun .or. size(work%phi, 2) < nslot) deallocate (work%phi)
      end if
      if (.not. allocated(work%phi)) allocate (work%phi(self%nfun, 0:nslot - 1))

      if (allocated(work%dphi)) then
         if (size(work%dphi, 1) /= self%nfun .or. size(work%dphi, 2) < nrhs) deallocate (work%dphi)
      end if
      if (.not. allocated(work%dphi)) allocate (work%dphi(self%nfun, 0:nrhs - 1))

      if (nrhs > 1) then
         if (allocated(work%panel)) then
            if (size(work%panel) /= self%nfun*iso_panel_width) deallocate (work%panel)
         end if
         if (.not. allocated(work%panel)) allocate (work%panel(self%nfun*iso_panel_width))
         if (.not. allocated(work%rows)) allocate (work%rows(iso_panel_width*moist_iso_gto_nrhs(4)))
      end if

      if (allocated(work%run_off)) then
         if (size(work%run_off) /= self%nshell) deallocate (work%run_off, work%run_len)
      end if
      if (.not. allocated(work%run_off)) allocate (work%run_off(self%nshell), work%run_len(self%nshell))
   end subroutine size_work

   !> Accumulate one shell's cartesian AO-derivative contributions at a point
   !>
   !> Fills the shell's components as rows 1..ncart(l) of `buf`, so screening, the
   !> primitive loop and the per-slot factor products live in one place
   !>
   !> @param[in]    self   Basis instance
   !> @param[in]    sh     Shell index
   !> @param[in]    d      Displacement from the shell center in Bohr
   !> @param[in]    u      Squared distance from the shell center
   !> @param[in]    ndloc  Highest spatial derivative order to fill (1..4)
   !> @param[inout] buf    Derivative block of the shell (accumulated into)
   pure subroutine gto_assemble_shell(self, sh, d, u, ndloc, buf)
      !> Basis instance
      class(moist_iso_gto_type), intent(in) :: self
      !> Shell index
      integer, intent(in) :: sh
      !> Displacement from the shell center
      real(wp), intent(in) :: d(3)
      !> Squared distance from the shell center
      real(wp), intent(in) :: u
      !> Highest spatial derivative order to fill
      integer, intent(in) :: ndloc
      !> Derivative block of the shell
      real(wp), intent(inout) :: buf(:, 0:)

      integer :: p, c, cbase, l, ncc, lx, ly, lz
      real(wp) :: a, ce
      real(wp) :: px(0:moist_iso_gto_lmax, 0:4)
      real(wp) :: py(0:moist_iso_gto_lmax, 0:4)
      real(wp) :: pz(0:moist_iso_gto_lmax, 0:4)

      l = self%sh_l(sh)
      ncc = moist_iso_gto_ncart(l)
      cbase = self%sh_coff(sh)
      do p = self%sh_poff(sh), self%sh_poff(sh + 1) - 1
         a = self%exps(p)
         ce = self%coeffs(p)*exp(-a*u)
         if (ce == 0.0_wp) cycle
         call moist_iso_gto_poly1d(d(1), a, ndloc, l, px)
         call moist_iso_gto_poly1d(d(2), a, ndloc, l, py)
         call moist_iso_gto_poly1d(d(3), a, ndloc, l, pz)
         do c = 1, ncc
            lx = self%comp_l(1, cbase + c)
            ly = self%comp_l(2, cbase + c)
            lz = self%comp_l(3, cbase + c)
            buf(c, 0) = buf(c, 0) + ce*px(lx, 0)*py(ly, 0)*pz(lz, 0)
            buf(c, 1) = buf(c, 1) + ce*px(lx, 1)*py(ly, 0)*pz(lz, 0)
            buf(c, 2) = buf(c, 2) + ce*px(lx, 0)*py(ly, 1)*pz(lz, 0)
            buf(c, 3) = buf(c, 3) + ce*px(lx, 0)*py(ly, 0)*pz(lz, 1)
            if (ndloc >= 2) then
               buf(c, 4) = buf(c, 4) + ce*px(lx, 2)*py(ly, 0)*pz(lz, 0)
               buf(c, 5) = buf(c, 5) + ce*px(lx, 1)*py(ly, 1)*pz(lz, 0)
               buf(c, 6) = buf(c, 6) + ce*px(lx, 1)*py(ly, 0)*pz(lz, 1)
               buf(c, 7) = buf(c, 7) + ce*px(lx, 0)*py(ly, 2)*pz(lz, 0)
               buf(c, 8) = buf(c, 8) + ce*px(lx, 0)*py(ly, 1)*pz(lz, 1)
               buf(c, 9) = buf(c, 9) + ce*px(lx, 0)*py(ly, 0)*pz(lz, 2)
            end if
            if (ndloc >= 3) then
               buf(c, 10) = buf(c, 10) + ce*px(lx, 3)*py(ly, 0)*pz(lz, 0)
               buf(c, 11) = buf(c, 11) + ce*px(lx, 2)*py(ly, 1)*pz(lz, 0)
               buf(c, 12) = buf(c, 12) + ce*px(lx, 2)*py(ly, 0)*pz(lz, 1)
               buf(c, 13) = buf(c, 13) + ce*px(lx, 1)*py(ly, 2)*pz(lz, 0)
               buf(c, 14) = buf(c, 14) + ce*px(lx, 1)*py(ly, 1)*pz(lz, 1)
               buf(c, 15) = buf(c, 15) + ce*px(lx, 1)*py(ly, 0)*pz(lz, 2)
               buf(c, 16) = buf(c, 16) + ce*px(lx, 0)*py(ly, 3)*pz(lz, 0)
               buf(c, 17) = buf(c, 17) + ce*px(lx, 0)*py(ly, 2)*pz(lz, 1)
               buf(c, 18) = buf(c, 18) + ce*px(lx, 0)*py(ly, 1)*pz(lz, 2)
               buf(c, 19) = buf(c, 19) + ce*px(lx, 0)*py(ly, 0)*pz(lz, 3)
            end if
            if (ndloc >= 4) then
               buf(c, 20) = buf(c, 20) + ce*px(lx, 4)*py(ly, 0)*pz(lz, 0)
               buf(c, 21) = buf(c, 21) + ce*px(lx, 3)*py(ly, 1)*pz(lz, 0)
               buf(c, 22) = buf(c, 22) + ce*px(lx, 3)*py(ly, 0)*pz(lz, 1)
               buf(c, 23) = buf(c, 23) + ce*px(lx, 2)*py(ly, 2)*pz(lz, 0)
               buf(c, 24) = buf(c, 24) + ce*px(lx, 2)*py(ly, 1)*pz(lz, 1)
               buf(c, 25) = buf(c, 25) + ce*px(lx, 2)*py(ly, 0)*pz(lz, 2)
               buf(c, 26) = buf(c, 26) + ce*px(lx, 1)*py(ly, 3)*pz(lz, 0)
               buf(c, 27) = buf(c, 27) + ce*px(lx, 1)*py(ly, 2)*pz(lz, 1)
               buf(c, 28) = buf(c, 28) + ce*px(lx, 1)*py(ly, 1)*pz(lz, 2)
               buf(c, 29) = buf(c, 29) + ce*px(lx, 1)*py(ly, 0)*pz(lz, 3)
               buf(c, 30) = buf(c, 30) + ce*px(lx, 0)*py(ly, 4)*pz(lz, 0)
               buf(c, 31) = buf(c, 31) + ce*px(lx, 0)*py(ly, 3)*pz(lz, 1)
               buf(c, 32) = buf(c, 32) + ce*px(lx, 0)*py(ly, 2)*pz(lz, 2)
               buf(c, 33) = buf(c, 33) + ce*px(lx, 0)*py(ly, 1)*pz(lz, 3)
               buf(c, 34) = buf(c, 34) + ce*px(lx, 0)*py(ly, 0)*pz(lz, 4)
            end if
         end do
      end do
   end subroutine gto_assemble_shell

   !> Density-weighted AO derivatives of the active functions, dphi = D phi
   !>
   !> Reads the active block of the density once for all columns of `dphi`:
   !>
   !> - every function active, in internal order: one gemm on the stored density
   !> - value and gradient (one column): run-pair loop over the upper block
   !>   triangle, each density block feeding both halves of the product
   !> - several columns: the active block packed panel by panel, one gemm each
   !>
   !> @param[in]    self     Basis instance
   !> @param[in]    phi      compact AO-derivative table, (nfun, 0:nslot-1)
   !> @param[out]   dphi     D phi of the first columns, (nact, 0:nrhs-1)
   !> @param[inout] panel    panel storage, (nfun*iso_panel_width)
   !> @param[inout] rows     panel-row storage of phi, (iso_panel_width*nrhs)
   !> @param[in]    run_off  internal offset of each active run
   !> @param[in]    run_len  function count of each active run
   subroutine gto_contract(self, phi, dphi, panel, rows, run_off, run_len)
      !> Basis instance
      class(moist_iso_gto_type), intent(in) :: self
      !> Compact AO-derivative table
      real(wp), intent(in), contiguous :: phi(:, 0:)
      !> Density-weighted AO derivatives
      real(wp), intent(out), contiguous :: dphi(:, 0:)
      !> Panel storage
      real(wp), intent(inout), contiguous, optional :: panel(:)
      !> Panel-row storage
      real(wp), intent(inout), contiguous, optional :: rows(:)
      !> Internal offset of each active run
      integer, intent(in) :: run_off(:)
      !> Function count of each active run
      integer, intent(in) :: run_len(:)

      integer :: nact, nrhs

      nact = size(dphi, 1)
      nrhs = size(dphi, 2)
      if (size(run_off) == 1 .and. run_off(1) == 0 .and. nact == self%nfun) then
         ! TODO: this gemm (and the one in contract_panels) runs per point inside
         ! the OpenMP projection; pthreaded OpenBLAS can fight with that
         ! (~2x slower for ~600 functions)
         call gemm(self%dens, phi(:, 0:nrhs - 1), dphi)
      else if (nrhs == 1) then
         call contract_runs(self, phi, dphi, run_off, run_len)
      else
         call contract_panels(self, phi, dphi, panel, rows, run_off, run_len)
      end if
   end subroutine gto_contract

   !> One-column contraction over the active runs, reading each off-diagonal run
   !> block of the density once for both of its symmetric halves
   !>
   !> @param[in]  self     Basis instance
   !> @param[in]  phi      compact AO-derivative table, (nfun, 0:nslot-1)
   !> @param[out] dphi     D phi, (nact, 0:0)
   !> @param[in]  run_off  internal offset of each active run
   !> @param[in]  run_len  function count of each active run
   subroutine contract_runs(self, phi, dphi, run_off, run_len)
      !> Basis instance
      class(moist_iso_gto_type), intent(in) :: self
      !> Compact AO-derivative table
      real(wp), intent(in) :: phi(:, 0:)
      !> Density-weighted AO derivatives
      real(wp), intent(out) :: dphi(:, 0:)
      !> Internal offset of each active run
      integer, intent(in) :: run_off(:)
      !> Function count of each active run
      integer, intent(in) :: run_len(:)

      integer :: ir, jr, oi, oj, ni, nj, pi, pj, c, i
      real(wp) :: p, acc

      dphi(:, 0) = 0.0_wp
      pj = 0
      do jr = 1, size(run_off)
         oj = run_off(jr)
         nj = run_len(jr)
         pi = 0
         do ir = 1, jr - 1
            oi = run_off(ir)
            ni = run_len(ir)
            do c = 1, nj
               p = phi(pj + c, 0)
               acc = 0.0_wp
               do i = 1, ni
                  dphi(pi + i, 0) = dphi(pi + i, 0) + self%dens(oi + i, oj + c)*p
                  acc = acc + self%dens(oi + i, oj + c)*phi(pi + i, 0)
               end do
               dphi(pj + c, 0) = dphi(pj + c, 0) + acc
            end do
            pi = pi + ni
         end do
         do c = 1, nj
            p = phi(pj + c, 0)
            do i = 1, nj
               dphi(pj + i, 0) = dphi(pj + i, 0) + self%dens(oj + i, oj + c)*p
            end do
         end do
         pj = pj + nj
      end do
   end subroutine contract_runs

   !> Several-column contraction: the active block of the density packed panel by
   !> panel of columns, each panel multiplied by its rows of phi in one gemm
   !>
   !> @param[in]    self     Basis instance
   !> @param[in]    phi      compact AO-derivative table, (nfun, 0:nslot-1)
   !> @param[out]   dphi     D phi of the first columns, (nact, 0:nrhs-1)
   !> @param[inout] panel    panel storage, (nfun*iso_panel_width)
   !> @param[inout] rows     panel-row storage of phi, (iso_panel_width*nrhs)
   !> @param[in]    run_off  internal offset of each active run
   !> @param[in]    run_len  function count of each active run
   subroutine contract_panels(self, phi, dphi, panel, rows, run_off, run_len)
      !> Basis instance
      class(moist_iso_gto_type), intent(in) :: self
      !> Compact AO-derivative table
      real(wp), intent(in) :: phi(:, 0:)
      !> Density-weighted AO derivatives
      real(wp), intent(out), contiguous :: dphi(:, 0:)
      !> Panel storage
      real(wp), intent(inout), contiguous, target :: panel(:)
      !> Panel-row storage
      real(wp), intent(inout), contiguous, target :: rows(:)
      !> Internal offset of each active run
      integer, intent(in) :: run_off(:)
      !> Function count of each active run
      integer, intent(in) :: run_len(:)

      !> Packed density columns of one panel, (nact, w)
      real(wp), pointer, contiguous :: dpan(:, :)
      !> Rows of phi matching the panel columns, (w, 0:nrhs-1)
      real(wp), pointer, contiguous :: ppan(:, :)
      integer :: nact, nrhs, q0, jr, joff, w, take, c, cj, pos, ir, oi, ni
      real(wp) :: beta

      nact = size(dphi, 1)
      nrhs = size(dphi, 2)
      q0 = 0
      jr = 1
      joff = 0
      beta = 0.0_wp
      do while (q0 < nact)
         ! Width of the next panel: up to iso_panel_width columns, splitting runs
         w = min(iso_panel_width, nact - q0)
         dpan(1:nact, 1:w) => panel
         ppan(1:w, 0:nrhs - 1) => rows
         ppan(:, :) = phi(q0 + 1:q0 + w, 0:nrhs - 1)
         c = 0
         do while (c < w)
            take = min(w - c, run_len(jr) - joff)
            do cj = run_off(jr) + joff + 1, run_off(jr) + joff + take
               c = c + 1
               pos = 0
               do ir = 1, size(run_off)
                  oi = run_off(ir)
                  ni = run_len(ir)
                  dpan(pos + 1:pos + ni, c) = self%dens(oi + 1:oi + ni, cj)
                  pos = pos + ni
               end do
            end do
            joff = joff + take
            if (joff == run_len(jr)) then
               jr = jr + 1
               joff = 0
            end if
         end do
         call gemm(dpan, ppan, dphi, beta=beta)
         beta = 1.0_wp
         q0 = q0 + w
      end do
   end subroutine contract_panels

   !> Gram entries phi(:, a) . dphi(:, b) over the active rows that the
   !> derivative assembly reads
   !>
   !> - column 0 against every slot; columns 1..3 against the gradient through
   !>   third-order slots; columns 4..9 against the Hessian slots
   !>
   !> @param[in]  phi    compact AO-derivative table, (nfun, 0:nslot-1)
   !> @param[in]  dphi   density-weighted AO derivatives, (nact, 0:nrhs-1)
   !> @param[in]  nslot  number of AO-derivative slots
   !> @param[out] gram   gram entries, (0:34, 0:9); unread entries undefined
   pure subroutine gto_gram(phi, dphi, nslot, gram)
      !> Compact AO-derivative table
      real(wp), intent(in) :: phi(:, 0:)
      !> Density-weighted AO derivatives
      real(wp), intent(in) :: dphi(:, 0:)
      !> Number of AO-derivative slots
      integer, intent(in) :: nslot
      !> Gram entries
      real(wp), intent(out) :: gram(0:34, 0:9)

      integer :: a, b, alo, ahi, nact

      nact = size(dphi, 1)
      do b = 0, size(dphi, 2) - 1
         select case (b)
         case (0)
            alo = 0
            ahi = nslot - 1
         case (1:3)
            alo = 1
            ahi = min(nslot - 1, 19)
         case default
            alo = 4
            ahi = 9
         end select
         do a = alo, ahi
            gram(a, b) = dot_product(phi(1:nact, a), dphi(:, b))
         end do
      end do
   end subroutine gto_gram

   !> Density Hessian from the gram entries
   !>
   !> d_ij rho = 2 (phi_ij . D phi + phi_i . D phi_j)
   !>
   !> @param[in]  gram   gram entries
   !> @param[out] d2rho  density Hessian, exactly symmetric
   pure subroutine density_hessian(gram, d2rho)
      !> Gram entries
      real(wp), intent(in) :: gram(0:34, 0:9)
      !> Density Hessian
      real(wp), intent(out) :: d2rho(3, 3)

      integer :: i, j

      do j = 1, 3
         do i = 1, j
            d2rho(i, j) = 2.0_wp*(gram(hess_slot(i, j), 0) + gram(i, j))
            d2rho(j, i) = d2rho(i, j)
         end do
      end do
   end subroutine density_hessian

   !> Density third derivative from the gram entries
   !>
   !> d_ijk rho = 2 (phi_ijk . D phi + phi_jk . D phi_i + phi_ik . D phi_j
   !> + phi_ij . D phi_k), the 10 unique components spread to all 27
   !>
   !> @param[in]  gram   gram entries
   !> @param[out] d3rho  density third derivative, exactly symmetric
   pure subroutine density_third(gram, d3rho)
      !> Gram entries
      real(wp), intent(in) :: gram(0:34, 0:9)
      !> Density third derivative
      real(wp), intent(out) :: d3rho(3, 3, 3)

      real(wp) :: unique(10:19)
      integer :: i, j, k

      do k = 1, 3
         do j = 1, k
            do i = 1, j
               unique(third_slot(i, j, k)) = 2.0_wp*(gram(third_slot(i, j, k), 0) &
                                                     + gram(hess_slot(j, k), i) &
                                                     + gram(hess_slot(i, k), j) &
                                                     + gram(hess_slot(i, j), k))
            end do
         end do
      end do
      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               d3rho(i, j, k) = unique(third_slot(i, j, k))
            end do
         end do
      end do
   end subroutine density_third

   !> Density fourth derivative from the gram entries
   !>
   !> Leibniz expansion of rho = phi^T D phi: one 4th-order term, four
   !> (3rd x 1st) terms and three (Hess x Hess) terms; the 15 unique components
   !> spread to all 81
   !>
   !> @param[in]  gram   gram entries
   !> @param[out] d4rho  density fourth derivative, exactly symmetric
   pure subroutine density_fourth(gram, d4rho)
      !> Gram entries
      real(wp), intent(in) :: gram(0:34, 0:9)
      !> Density fourth derivative
      real(wp), intent(out) :: d4rho(3, 3, 3, 3)

      real(wp) :: unique(20:34)
      integer :: i, j, k, l

      do l = 1, 3
         do k = 1, l
            do j = 1, k
               do i = 1, j
                  unique(fourth_slot(i, j, k, l)) = 2.0_wp*(gram(fourth_slot(i, j, k, l), 0) &
                                                            + gram(third_slot(j, k, l), i) &
                                                            + gram(third_slot(i, k, l), j) &
                                                            + gram(third_slot(i, j, l), k) &
                                                            + gram(third_slot(i, j, k), l) &
                                                            + gram(hess_slot(i, j), hess_slot(k, l)) &
                                                            + gram(hess_slot(i, k), hess_slot(j, l)) &
                                                            + gram(hess_slot(i, l), hess_slot(j, k)))
               end do
            end do
         end do
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               do i = 1, 3
                  d4rho(i, j, k, l) = unique(fourth_slot(i, j, k, l))
               end do
            end do
         end do
      end do
   end subroutine density_fourth

   ! >>> GENERATED poly1d >>>
   !> Fill the exp-stripped 1D Gaussian-factor derivative table
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

   ! >>> GENERATED fold >>>
   !> Fold one shell's cartesian AO-derivative block into moist's spherical basis
   !>
   !> ``sph(m, k) = sum_c U(c, m) cart(c, k)``, where spherical function m is the
   !> harmonic extension of the m-th leading monomial (``lx <= 1``): that monomial
   !> plus terms with ``lx >= 2``
   !>
   !> @param[in]    l      shell angular momentum, 2..lmax
   !> @param[in]    nslot  number of derivative slots to fold
   !> @param[in]    cart   cartesian block, shape (ncart(l), 0:nslot-1)
   !> @param[inout] sph    spherical block, shape (2l+1, 0:nslot-1)
   pure subroutine moist_iso_gto_fold(l, nslot, cart, sph)
      !> Shell angular momentum
      integer, intent(in) :: l
      !> Number of derivative slots
      integer, intent(in) :: nslot
      !> Cartesian AO-derivative block
      real(wp), intent(in) :: cart(:, 0:)
      !> Spherical AO-derivative block
      real(wp), intent(inout) :: sph(:, 0:)

      integer :: k

      select case (l)
      case (2)
         do k = 0, nslot - 1
            sph(1, k) = cart(2, k)
            sph(2, k) = cart(3, k)
            sph(3, k) = cart(4, k) - cart(1, k)
            sph(4, k) = cart(5, k)
            sph(5, k) = cart(6, k) - cart(1, k)
         end do
      case (3)
         do k = 0, nslot - 1
            sph(1, k) = cart(4, k) - 0.3333333333333333_wp*cart(1, k)
            sph(2, k) = cart(5, k)
            sph(3, k) = cart(6, k) - 0.3333333333333333_wp*cart(1, k)
            sph(4, k) = cart(7, k) - 3.0_wp*cart(2, k)
            sph(5, k) = cart(8, k) - cart(3, k)
            sph(6, k) = cart(9, k) - cart(2, k)
            sph(7, k) = cart(10, k) - 3.0_wp*cart(3, k)
         end do
      case (4)
         do k = 0, nslot - 1
            sph(1, k) = cart(7, k) - cart(2, k)
            sph(2, k) = cart(8, k) - 0.3333333333333333_wp*cart(3, k)
            sph(3, k) = cart(9, k) - 0.3333333333333333_wp*cart(2, k)
            sph(4, k) = cart(10, k) - cart(3, k)
            sph(5, k) = cart(11, k) + cart(1, k) - 6.0_wp*cart(4, k)
            sph(6, k) = cart(12, k) - 3.0_wp*cart(5, k)
            sph(7, k) = cart(13, k) + 0.3333333333333333_wp*cart(1, k) - cart(4, k) - cart(6, k)
            sph(8, k) = cart(14, k) - 3.0_wp*cart(5, k)
            sph(9, k) = cart(15, k) + cart(1, k) - 6.0_wp*cart(6, k)
         end do
      case (5)
         do k = 0, nslot - 1
            sph(1, k) = cart(11, k) + 0.2_wp*cart(1, k) - 2.0_wp*cart(4, k)
            sph(2, k) = cart(12, k) - cart(5, k)
            sph(3, k) = &
               cart(13, k) + 0.06666666666666667_wp*cart(1, k) - 0.3333333333333333_wp*cart(4, k) &
               - 0.3333333333333333_wp*cart(6, k)
            sph(4, k) = cart(14, k) - cart(5, k)
            sph(5, k) = cart(15, k) + 0.2_wp*cart(1, k) - 2.0_wp*cart(6, k)
            sph(6, k) = cart(16, k) + 5.0_wp*cart(2, k) - 10.0_wp*cart(7, k)
            sph(7, k) = cart(17, k) + cart(3, k) - 6.0_wp*cart(8, k)
            sph(8, k) = cart(18, k) + cart(2, k) - cart(7, k) - 3.0_wp*cart(9, k)
            sph(9, k) = cart(19, k) + cart(3, k) - 3.0_wp*cart(8, k) - cart(10, k)
            sph(10, k) = cart(20, k) + cart(2, k) - 6.0_wp*cart(9, k)
            sph(11, k) = cart(21, k) + 5.0_wp*cart(3, k) - 10.0_wp*cart(10, k)
         end do
      case (6)
         do k = 0, nslot - 1
            sph(1, k) = cart(16, k) + cart(2, k) - 3.3333333333333335_wp*cart(7, k)
            sph(2, k) = cart(17, k) + 0.2_wp*cart(3, k) - 2.0_wp*cart(8, k)
            sph(3, k) = cart(18, k) + 0.2_wp*cart(2, k) - 0.3333333333333333_wp*cart(7, k) - cart(9, k)
            sph(4, k) = cart(19, k) + 0.2_wp*cart(3, k) - cart(8, k) - 0.3333333333333333_wp*cart(10, k)
            sph(5, k) = cart(20, k) + 0.2_wp*cart(2, k) - 2.0_wp*cart(9, k)
            sph(6, k) = cart(21, k) + cart(3, k) - 3.3333333333333335_wp*cart(10, k)
            sph(7, k) = cart(22, k) - cart(1, k) + 15.0_wp*cart(4, k) - 15.0_wp*cart(11, k)
            sph(8, k) = cart(23, k) + 5.0_wp*cart(5, k) - 10.0_wp*cart(12, k)
            sph(9, k) = &
               cart(24, k) - 0.2_wp*cart(1, k) + 2.0_wp*cart(4, k) + cart(6, k) - cart(11, k) - 6.0_wp*cart(13, k)
            sph(10, k) = cart(25, k) + 3.0_wp*cart(5, k) - 3.0_wp*cart(12, k) - 3.0_wp*cart(14, k)
            sph(11, k) = &
               cart(26, k) - 0.2_wp*cart(1, k) + cart(4, k) + 2.0_wp*cart(6, k) - 6.0_wp*cart(13, k) - cart(15, k)
            sph(12, k) = cart(27, k) + 5.0_wp*cart(5, k) - 10.0_wp*cart(14, k)
            sph(13, k) = cart(28, k) - cart(1, k) + 15.0_wp*cart(6, k) - 15.0_wp*cart(15, k)
         end do
      case (7)
         do k = 0, nslot - 1
            sph(1, k) = cart(22, k) - 0.14285714285714285_wp*cart(1, k) + 3.0_wp*cart(4, k) - 5.0_wp*cart(11, k)
            sph(2, k) = cart(23, k) + cart(5, k) - 3.3333333333333335_wp*cart(12, k)
            sph(3, k) = &
               cart(24, k) - 0.02857142857142857_wp*cart(1, k) + 0.4_wp*cart(4, k) + 0.2_wp*cart(6, k) &
               - 0.3333333333333333_wp*cart(11, k) - 2.0_wp*cart(13, k)
            sph(4, k) = cart(25, k) + 0.6_wp*cart(5, k) - cart(12, k) - cart(14, k)
            sph(5, k) = &
               cart(26, k) - 0.02857142857142857_wp*cart(1, k) + 0.2_wp*cart(4, k) + 0.4_wp*cart(6, k) &
               - 2.0_wp*cart(13, k) - 0.3333333333333333_wp*cart(15, k)
            sph(6, k) = cart(27, k) + cart(5, k) - 3.3333333333333335_wp*cart(14, k)
            sph(7, k) = cart(28, k) - 0.14285714285714285_wp*cart(1, k) + 3.0_wp*cart(6, k) - 5.0_wp*cart(15, k)
            sph(8, k) = cart(29, k) - 7.0_wp*cart(2, k) + 35.0_wp*cart(7, k) - 21.0_wp*cart(16, k)
            sph(9, k) = cart(30, k) - cart(3, k) + 15.0_wp*cart(8, k) - 15.0_wp*cart(17, k)
            sph(10, k) = &
               cart(31, k) - cart(2, k) + 3.3333333333333335_wp*cart(7, k) + 5.0_wp*cart(9, k) &
               - cart(16, k) - 10.0_wp*cart(18, k)
            sph(11, k) = &
               cart(32, k) - 0.6_wp*cart(3, k) + 6.0_wp*cart(8, k) + cart(10, k) - 3.0_wp*cart(17, k) &
               - 6.0_wp*cart(19, k)
            sph(12, k) = &
               cart(33, k) - 0.6_wp*cart(2, k) + cart(7, k) + 6.0_wp*cart(9, k) - 6.0_wp*cart(18, k) &
               - 3.0_wp*cart(20, k)
            sph(13, k) = &
               cart(34, k) - cart(3, k) + 5.0_wp*cart(8, k) + 3.3333333333333335_wp*cart(10, k) &
               - 10.0_wp*cart(19, k) - cart(21, k)
            sph(14, k) = cart(35, k) - cart(2, k) + 15.0_wp*cart(9, k) - 15.0_wp*cart(20, k)
            sph(15, k) = cart(36, k) - 7.0_wp*cart(3, k) + 35.0_wp*cart(10, k) - 21.0_wp*cart(21, k)
         end do
      case (8)
         do k = 0, nslot - 1
            sph(1, k) = cart(29, k) - cart(2, k) + 7.0_wp*cart(7, k) - 7.0_wp*cart(16, k)
            sph(2, k) = cart(30, k) - 0.14285714285714285_wp*cart(3, k) + 3.0_wp*cart(8, k) - 5.0_wp*cart(17, k)
            sph(3, k) = &
               cart(31, k) - 0.14285714285714285_wp*cart(2, k) + 0.6666666666666666_wp*cart(7, k) + cart(9, k) - &
               0.3333333333333333_wp*cart(16, k) - 3.3333333333333335_wp*cart(18, k)
            sph(4, k) = &
               cart(32, k) - 0.08571428571428572_wp*cart(3, k) + 1.2_wp*cart(8, k) &
               + 0.2_wp*cart(10, k) - cart(17, k) - 2.0_wp*cart(19, k)
            sph(5, k) = &
               cart(33, k) - 0.08571428571428572_wp*cart(2, k) + 0.2_wp*cart(7, k) + 1.2_wp*cart(9, k) &
               - 2.0_wp*cart(18, k) - &
               cart(20, k)
            sph(6, k) = &
               cart(34, k) - 0.14285714285714285_wp*cart(3, k) + cart(8, k) + 0.6666666666666666_wp*cart(10, k) - &
               3.3333333333333335_wp*cart(19, k) - 0.3333333333333333_wp*cart(21, k)
            sph(7, k) = cart(35, k) - 0.14285714285714285_wp*cart(2, k) + 3.0_wp*cart(9, k) - 5.0_wp*cart(20, k)
            sph(8, k) = cart(36, k) - cart(3, k) + 7.0_wp*cart(10, k) - 7.0_wp*cart(21, k)
            sph(9, k) = cart(37, k) + cart(1, k) - 28.0_wp*cart(4, k) + 70.0_wp*cart(11, k) - 28.0_wp*cart(22, k)
            sph(10, k) = cart(38, k) - 7.0_wp*cart(5, k) + 35.0_wp*cart(12, k) - 21.0_wp*cart(23, k)
            sph(11, k) = &
               cart(39, k) + 0.14285714285714285_wp*cart(1, k) - 3.0_wp*cart(4, k) - cart(6, k) + 5.0_wp*cart(11, k) + &
               15.0_wp*cart(13, k) - cart(22, k) - 15.0_wp*cart(24, k)
            sph(12, k) = &
               cart(40, k) - 3.0_wp*cart(5, k) + 10.0_wp*cart(12, k) + 5.0_wp*cart(14, k) &
               - 3.0_wp*cart(23, k) - 10.0_wp*cart(25, k)
            sph(13, k) = &
               cart(41, k) + 0.08571428571428572_wp*cart(1, k) - 1.2_wp*cart(4, k) - 1.2_wp*cart(6, k) &
               + cart(11, k) + 12.0_wp*cart(13, k) + cart(15, k) - 6.0_wp*cart(24, k) - 6.0_wp*cart(26, k)
            sph(14, k) = &
               cart(42, k) - 3.0_wp*cart(5, k) + 5.0_wp*cart(12, k) + 10.0_wp*cart(14, k) &
               - 10.0_wp*cart(25, k) - 3.0_wp*cart(27, k)
            sph(15, k) = &
               cart(43, k) + 0.14285714285714285_wp*cart(1, k) - cart(4, k) - 3.0_wp*cart(6, k) &
                + 15.0_wp*cart(13, k) + 5.0_wp*cart(15, k) - 15.0_wp*cart(26, k) - cart(28, k)
            sph(16, k) = cart(44, k) - 7.0_wp*cart(5, k) + 35.0_wp*cart(14, k) - 21.0_wp*cart(27, k)
            sph(17, k) = cart(45, k) + cart(1, k) - 28.0_wp*cart(6, k) + 70.0_wp*cart(15, k) - 28.0_wp*cart(28, k)
         end do
      case default
         error stop "moist_iso_gto_fold: angular momentum out of range"
      end select
   end subroutine moist_iso_gto_fold

end module moist_cavity_drop_lsf_isodensity_gto
