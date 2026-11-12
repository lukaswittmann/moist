!> Internal isodensity level set function for DROP
!>
!> Internal variant of [[moist_cavity_drop_lsf_isodensity_callback_type]]
!>
!> - owns the cartesian-monomial Gaussian basis and the density matrix, and
!>   evaluates the level set internally
!>
!> Relies on the GTO code in [[moist_cavity_drop_lsf_isodensity_gto]]
!>
!> Level set function
!>
!>    S(r) = scale * (rho_iso - rho(r))
!>
!> The host code has to provide
!> - Basis set (once, at construction)
!> - Density matrix (once per SCF step), as the `density_matrix` answer of the
!>   build phase, read by `read_build` before the projection
!> - Nuclear charges of the core fill, as the `nuclear_charges` answer, while
!>   the fill is on
!>
!> The isodensity surface is a level set of the density itself, so the level set
!> function carries no *explicit* nuclear-position dependence for the cavity
!> chain rule, and the mixed spatial/nuclear derivatives vanish here
module moist_cavity_drop_lsf_isodensity_internal
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, &
                                         lsf_base_update, lsf_candidate_space_user, &
                                         lsf_batch_require_outputs, lsf_sorted_candidates, &
                                         lsf_vjp_batch_shape_message
   use moist_cavity_drop_lsf_isodensity_gto, only: moist_iso_gto_type, moist_iso_gto_work_type, &
                                                   moist_iso_gto_tile_type, moist_iso_gto_tile_width
   use moist_math_packed_sym, only: pack_sym_matrix, pack_sym_tensor3, pack_sym_tensor4
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use moist_cavity_drop_lsf_isodensity_param, only: &
      moist_cavity_drop_lsf_isodensity_param_type, isodensity_exclusion_radius, &
      isodensity_max_nuclear_charge
   use moist_cavity_drop_lsf_isodensity_core, only: moist_iso_core_fill_type, iso_check_nuclear_charges
   use moist_channels_coupling, only: coupling_type, coupling_view_type, coupling_register, &
      & request_require, moist_phase_build, nuclear_charge_request_type, density_matrix_request_type
   implicit none(type, external)
   private

   integer, parameter :: ndim = 3

   !> Widest tile of [[lsf_evaluate_batch]]
   integer, parameter :: max_tile = maxval(moist_iso_gto_tile_width)
   !> Bits per axis of the Morton cells ordering the batch columns into tiles
   integer, parameter :: order_bits = 4

   public :: moist_cavity_drop_lsf_isodensity_internal_type

   !> Isodensity LSF backed by an internal cartesian-monomial GTO evaluator
   type, extends(moist_cavity_drop_lsf_type) :: moist_cavity_drop_lsf_isodensity_internal_type
      !> Cartesian-monomial Gaussian basis and density matrix
      type(moist_iso_gto_type) :: gto
      !> Level set parameters (isovalue and constant multiplier)
      type(moist_cavity_drop_lsf_isodensity_param_type) :: param
      !> Cached evaluation point in Bohr
      real(wp) :: point(ndim) = 0.0_wp
      !> Cached LSF value
      real(wp) :: value = 0.0_wp
      !> Cached spatial gradient
      real(wp) :: grad(ndim) = 0.0_wp
      !> Cached spatial Hessian
      real(wp) :: hess(ndim, ndim) = 0.0_wp
      !> Cached third spatial derivative
      real(wp) :: third(ndim, ndim, ndim) = 0.0_wp
      !> Cached fourth spatial derivative
      real(wp) :: fourth(ndim, ndim, ndim, ndim) = 0.0_wp
      !> Highest requested derivative order
      integer :: max_deriv = 0
      !> Largest nuclear charge in the bound structure
      real(wp) :: zmax = 0.0_wp
      !> Core fill of the ECP nuclei in the bound structure
      type(moist_iso_core_fill_type) :: core
      !> Nuclear charges of the host Hamiltonian (nat); unallocated means no ECP nucleus
      real(wp), allocatable :: z(:)
      !> Per-instance evaluation scratch, sized on first use
      type(moist_iso_gto_work_type) :: work
      !> Per-instance tile scratch of `evaluate_batch`, grown on demand by `eval_tile`
      !>
      !> - empty until the first batch, so a clone of a fresh model copies nothing
      type(moist_iso_gto_tile_type) :: tile
      !> Ascending copy of the shared candidate list of `evaluate_batch`
      integer, allocatable :: bat_cand(:)
      !> Morton cell of each batch column, grow-only
      integer, allocatable :: bat_key(:)
      !> Batch columns in tile order, grow-only
      integer, allocatable :: bat_perm(:)
   contains
      procedure, public :: new => lsf_new
      procedure, public :: set_density => lsf_set_density
      procedure, public :: set_nuclear_charges => lsf_set_nuclear_charges
      procedure, public :: declare_coupling => lsf_declare_coupling
      procedure, public :: read_build => lsf_read_build
      procedure, public :: requires_build_inputs => lsf_requires_build_inputs
      procedure, public :: update => lsf_update
      procedure, public :: prepare => lsf_prepare
      procedure, public :: prepare_subset => lsf_prepare_subset
      procedure, public :: set_max_deriv => lsf_set_max_deriv
      procedure, public :: active_count => lsf_active_count
      procedure, public :: active_atom => lsf_active_atom
      procedure, public :: f0 => lsf_f0
      procedure, public :: f012_r => lsf_f012_r
      procedure, public :: f3_rrr => lsf_f3_rrr
      procedure, public :: f4_rrrr => lsf_f4_rrrr
      procedure, public :: evaluate_batch => lsf_evaluate_batch
      procedure, public :: f3_rr_rA => lsf_f3_rr_rA
      procedure, public :: vjp_f1_rA => lsf_vjp_f1_rA
      procedure, public :: vjp_f1_rA_batch => lsf_vjp_f1_rA_batch
      procedure, public :: screening_offset => lsf_screening_offset
      procedure, public :: exclusion_radius => lsf_exclusion_radius
      procedure, public :: density_adjoint_factor => lsf_density_adjoint_factor
      procedure, public :: check_surface => lsf_check_surface
      procedure, public :: has_atom_partition => lsf_has_atom_partition
      procedure, public :: atom_partition => lsf_atom_partition
   end type moist_cavity_drop_lsf_isodensity_internal_type

contains

   !> Construct from parameter values; omission uses compiled defaults
   !>
   !> @param[inout] self Object to initialize
   !> @param[in]    sh_atom Shell atom indices
   !> @param[in]    sh_l Shell angular momenta
   !> @param[in]    sh_nprim Primitive count per shell
   !> @param[in]    exps Gaussian exponents
   !> @param[in]    coeffs Gaussian coefficients
   !> @param[out]   error Construction error
   !> @param[in]    param Configuration copied by value
   subroutine lsf_new(self, sh_atom, sh_l, sh_nprim, exps, coeffs, error, param)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      integer, intent(in) :: sh_atom(:)
      integer, intent(in) :: sh_l(:)
      integer, intent(in) :: sh_nprim(:)
      real(wp), intent(in) :: exps(:)
      real(wp), intent(in) :: coeffs(:)
      type(error_type), allocatable, intent(out) :: error

      !> Configuration; omitted means compiled defaults
      type(moist_cavity_drop_lsf_isodensity_param_type), intent(in), optional :: param

      !> Candidate ids index this LSF's per-atom GTO shells
      self%candidate_space = lsf_candidate_space_user
      self%name = "Isodensity"

      self%radius_dependent = .false.

      call self%param%init_defaults()
      if (present(param)) self%param = param
      call self%gto%init(sh_atom, sh_l, sh_nprim, exps, coeffs, error, basis=self%param%basis)
      if (allocated(error)) return
      !> Scratch sized for the previous basis; resized on first use
      self%work = moist_iso_gto_work_type()
      self%tile = moist_iso_gto_tile_type()
   end subroutine lsf_new

   !> Install the cartesian-monomial density matrix for the current SCF step
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    dcart Density matrix in the cartesian-monomial basis
   !> @param[out]   error Set on a size mismatch or, in the spherical
   !>                     representation, a density outside the spherical subspace
   subroutine lsf_set_density(self, dcart, error)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: dcart(:, :)
      type(error_type), allocatable, intent(out) :: error

      call self%gto%set_density(dcart, error)
   end subroutine lsf_set_density

   !> Declare the build-phase requests at the cavity scope
   !>
   !> - `density_matrix`, required: the density this level set is built from
   !> - `nuclear_charges`, required while the core fill is on
   !>
   !> @param[in]    self     LSF instance
   !> @param[inout] coupling Coupling being declared
   !> @param[out]   error    Error handling
   subroutine lsf_declare_coupling(self, coupling, error)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      type(coupling_type), intent(inout) :: coupling
      type(error_type), allocatable, intent(out) :: error

      type(density_matrix_request_type) :: density
      type(nuclear_charge_request_type) :: charges

      ! The density matrix runs over the cartesian-monomial basis
      density%nbasis = self%gto%ncart
      call request_require(density, moist_phase_build, "dcart", error)
      if (allocated(error)) return
      call coupling_register(coupling, "density_matrix", density, error)
      if (allocated(error)) return
      call request_require(charges, moist_phase_build, "z", self%param%fills_core(), error)
      if (allocated(error)) return
      call coupling_register(coupling, "nuclear_charges", charges, error)
   end subroutine lsf_declare_coupling

   !> Read the density matrix and, with the core fill on, the nuclear charges
   !>
   !> @param[inout] self     LSF instance
   !> @param[in]    coupling Cavity-scoped view of the build phase
   !> @param[in]    mol      Structure of the cavity's latest update
   !> @param[out]   error    Missing answer, a density the basis rejects or a
   !>                        charge outside (0, Z]
   subroutine lsf_read_build(self, coupling, mol, error)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: dcart(:, :), z(:)

      call coupling%read("density_matrix", "dcart", dcart, error)
      if (allocated(error)) return
      call self%set_density(dcart, error)
      if (allocated(error)) return
      if (.not. self%param%fills_core()) return
      call coupling%read("nuclear_charges", "z", z, error)
      if (allocated(error)) return
      call iso_check_nuclear_charges(mol, z, error)
      if (allocated(error)) return
      call self%set_nuclear_charges(z)
   end subroutine lsf_read_build

   !> Whether the projection waits for build-phase answers
   !>
   !> The density matrix, always
   !>
   !> @param[in] self LSF instance
   pure function lsf_requires_build_inputs(self) result(requires)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      logical :: requires

      requires = .true.
   end function lsf_requires_build_inputs

   !> Store the nuclear charges of the host Hamiltonian for the core fill
   !>
   !> Read by the next `update`; atoms with `z < Z` are ECP nuclei
   !>
   !> @param[inout] self LSF instance
   !> @param[in]    z    Nuclear charges (nat)
   subroutine lsf_set_nuclear_charges(self, z)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: z(:)

      self%z = z
   end subroutine lsf_set_nuclear_charges

   !> Bind molecular geometry and refresh the shell centers
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    mol   Molecular structure
   !> @param[in]    radii Per-atom radii
   subroutine lsf_update(self, mol, radii)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)

      call lsf_base_update(self, mol, radii)
      call self%gto%refresh_centers(mol)
      !> Size the per-shell radial screening cutoffs to the cavity's screening threshold (inherited from the base)
      call self%gto%build_screening(self%screening_threshold)
      self%zmax = isodensity_max_nuclear_charge(mol)
      call self%core%update(mol, radii, self%param, self%z)
      if (allocated(self%bat_cand)) then
         if (size(self%bat_cand) < mol%nat) deallocate (self%bat_cand)
      end if
      if (.not. allocated(self%bat_cand)) allocate (self%bat_cand(mol%nat))
   end subroutine lsf_update

   !> Evaluate and cache the level set at one point via the internal GTO evaluator
   !>
   !> Only the derivative orders demanded by `max_deriv` are computed
   !>
   !> - never unreached: `status = 0`, the density with no shell in reach being
   !>   a defined zero
   !>
   !> @param[inout] self   LSF instance
   !> @param[in]    point  Evaluation point in Bohr
   !> @param[out]   error  Set when the point lies in an ECP core fill outside the cavity
   !> @param[out]   status Always 0 (optional)
   subroutine lsf_prepare(self, point, error, status)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(out), optional :: status

      if (present(status)) status = 0
      call lsf_prepare_impl(self, point, error)
   end subroutine lsf_prepare

   !> Shared prepare body: evaluate the level set at `point` and cache it
   !>
   !> When `cand_atoms` is present only the shells owned by those atoms are computed
   !>
   !> @param[inout] self       LSF instance
   !> @param[in]    point      Evaluation point in Bohr
   !> @param[out]   error      Set when the point lies in an ECP core fill outside the cavity
   !> @param[in]    cand_atoms Optional candidate owner-atom list
   subroutine lsf_prepare_impl(self, point, error, cand_atoms)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: cand_atoms(:)

      real(wp) :: rho, drho(ndim), d2rho(ndim, ndim), d3rho(ndim, ndim, ndim)
      real(wp) :: d4rho(ndim, ndim, ndim, ndim)
      integer :: nderiv, inside

      nderiv = min(max(self%max_deriv, 1), 4)

      select case (nderiv)
      case (:1)
         call self%gto%eval(point, self%work, rho, drho, cand_atoms=cand_atoms)
      case (2)
         call self%gto%eval(point, self%work, rho, drho, d2rho=d2rho, cand_atoms=cand_atoms)
      case (3)
         call self%gto%eval(point, self%work, rho, drho, d2rho=d2rho, d3rho=d3rho, &
                            cand_atoms=cand_atoms)
      case default
         call self%gto%eval(point, self%work, rho, drho, d2rho=d2rho, d3rho=d3rho, &
                            cand_atoms=cand_atoms, d4rho=d4rho)
      end select
      ! Every candidate list holds the shells of a point's nearby nuclei, so the
      ! fill is added whatever the screening
      call self%core%add(point, nderiv, rho, drho, d2rho, d3rho, d4rho, inside)
      call self%core%check(inside, rho, error)

      self%point = point

      self%prepared_deriv = nderiv
      self%value = self%param%scale*(self%param%rho_iso - rho)
      self%grad = -self%param%scale*drho
      if (nderiv >= 2) then
         self%hess = -self%param%scale*d2rho
      else
         self%hess = 0.0_wp
      end if
      if (nderiv >= 3) then
         self%third = -self%param%scale*d3rho
      else
         self%third = 0.0_wp
      end if
      if (nderiv >= 4) self%fourth = -self%param%scale*d4rho
   end subroutine lsf_prepare_impl

   !> Evaluate and cache the level set over the shells of a candidate atom list
   !>
   !> - the list is honoured: only the shells owned by the listed atoms are
   !>   evaluated (user-space ids), plus the core fill whatever the list
   !> - never unreached, as [[lsf_prepare]]
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point in Bohr
   !> @param[in]    candidate_indices Candidate owner atoms, user space
   !> @param[out]   error             Set when the point lies in an ECP core fill
   !>                                 outside the cavity
   !> @param[out]   status            Always 0 (optional)
   subroutine lsf_prepare_subset(self, point, candidate_indices, error, status)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: candidate_indices(:)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(out), optional :: status

      if (present(status)) status = 0
      call lsf_prepare_impl(self, point, error, cand_atoms=candidate_indices)
   end subroutine lsf_prepare_subset

   !> Record requested derivative order
   !>
   !> @param[inout] self LSF instance
   !> @param[in]    n    Requested max derivative order
   subroutine lsf_set_max_deriv(self, n)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      integer, intent(in) :: n

      self%max_deriv = max(0, n)
   end subroutine lsf_set_max_deriv

   !> Number of active atoms, true-density LSFs not being atom screened
   !>
   !> @param[in] self LSF instance
   pure function lsf_active_count(self) result(n)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      integer :: n

      if (self%ncenters < 0) then
         n = -1
         return
      end if
      n = 0
   end function lsf_active_count

   !> Active atom lookup, undefined for zero active atoms, returning a zero sentinel
   !>
   !> @param[in] self LSF instance
   !> @param[in] i    Active-list index
   pure function lsf_active_atom(self, i) result(idx)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      integer, intent(in) :: i
      integer :: idx

      if (self%ncenters < 0 .or. i < 0) then
         idx = -1
         return
      end if
      idx = 0
   end function lsf_active_atom

   !> Return cached LSF value
   !>
   !> @param[in]  self LSF instance
   !> @param[out] val  LSF value
   subroutine lsf_f0(self, val)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out) :: val

      val = self%value
   end subroutine lsf_f0

   !> Return cached LSF value, gradient, and Hessian
   !>
   !> @param[in]  self    LSF instance
   !> @param[out] lsf0    Optional LSF value
   !> @param[out] lsf1_r  Optional spatial gradient
   !> @param[out] lsf2_rr Optional packed spatial Hessian [6]
   subroutine lsf_f012_r(self, lsf0, lsf1_r, lsf2_rr)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:)

      if (present(lsf0)) lsf0 = self%value
      if (present(lsf1_r)) lsf1_r(:) = self%grad(:)
      if (present(lsf2_rr)) then
         call self%require_deriv(2, "f012_r(lsf2_rr)")
         call pack_sym_matrix(self%hess, lsf2_rr, ndim)
      end if
   end subroutine lsf_f012_r

   !> Return cached lower derivatives and the third spatial derivative
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] lsf0     Optional LSF value
   !> @param[out] lsf1_r   Optional spatial gradient
   !> @param[out] lsf2_rr  Optional packed spatial Hessian [6]
   !> @param[out] lsf3_rrr Packed spatial third derivative [10]
   subroutine lsf_f3_rrr(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:)
      real(wp), intent(out) :: lsf3_rrr(:)

      call self%require_deriv(3, "f3_rrr")
      call self%f012_r(lsf0, lsf1_r, lsf2_rr)
      call pack_sym_tensor3(self%third, lsf3_rrr, ndim)
   end subroutine lsf_f3_rrr

   !> Return the cached fourth spatial derivative
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf4_rrrr Packed fourth spatial derivative [15]
   subroutine lsf_f4_rrrr(self, lsf4_rrrr)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rrrr(:)

      call self%require_deriv(4, "f4_rrrr")
      call pack_sym_tensor4(self%fourth, lsf4_rrrr, ndim)
   end subroutine lsf_f4_rrrr

   !* ================================================================================= *!
   !*                                Batched evaluation                                 *!
   !* ================================================================================= *!

   !> Level-set jet at a batch of points through the tile contraction
   !>
   !> Override of [[lsf_base_evaluate_batch]], same contract:
   !>
   !> - columns in the caller's order, or by Morton cell over the batch's
   !>   bounding box when that makes the tiles more compact ([[tile_order]]),
   !>   then cut into tiles of `moist_iso_gto_tile_width(order)`, so a scattered
   !>   batch still shares shells within a tile; per tile the AO
   !>   derivatives of every point, one level-3 contraction against the
   !>   density, then per point the core fill and the lift
   !>   ([[moist_iso_gto_type:eval_tile]]); `order` alone governs (no
   !>   `require_deriv`, `max_deriv` untouched); orders 0 and 1 both evaluate
   !>   the gradient, as `prepare` does
   !> - screening is that of `prepare` / `prepare_subset` per point: the
   !>   per-(shell, point) radial cutoff, over the shells of the listed atoms
   !>   when `candidate_indices` is present (unique ids in 1..ncenters, sorted
   !>   ascending before use; atoms without shells add nothing), of every atom
   !>   when absent; every
   !>   column equals the accessors after `prepare` up to the summation order
   !>   of the contraction
   !> - `status` never 1: an isodensity point is never unreached, the density
   !>   with no shell in reach being a defined zero
   !> - `error`: a point inside an ECP core fill that is not inside the cavity,
   !>   or no density installed; names the first failing column, every failing
   !>   column's outputs unspecified and its `status` 2
   !> - the single-point cache is invalidated (`prepared_deriv = -1`): accessors
   !>   abort until the next `prepare`
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    points            Evaluation points (3, n)
   !> @param[in]    order             Highest spatial derivative to write (0..4)
   !> @param[out]   lsf0              Values (n)
   !> @param[out]   lsf1_r            Gradients (3, n), written for order >= 1
   !> @param[out]   lsf2_rr_packed    Packed Hessians (6, n), written for order >= 2
   !> @param[out]   lsf3_rrr          Packed third derivatives (10, n), written for order >= 3
   !> @param[out]   lsf4_rrrr         Packed fourth derivatives (15, n), written for order 4
   !> @param[out]   status            Per column 0 ok, 2 refused by the core fill (optional)
   !> @param[in]    candidate_indices Shared atom list in `candidate_space` (optional)
   !> @param[out]   error             First failing column
   subroutine lsf_evaluate_batch(self, points, order, lsf0, lsf1_r, lsf2_rr_packed, lsf3_rrr, &
                                 lsf4_rrrr, status, candidate_indices, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      !> Evaluation points
      real(wp), intent(in), contiguous :: points(:, :)
      !> Highest spatial derivative to write
      integer, intent(in) :: order
      !> Values
      real(wp), intent(out), contiguous :: lsf0(:)
      !> Gradients
      real(wp), intent(out), optional, contiguous :: lsf1_r(:, :)
      !> Packed Hessians
      real(wp), intent(out), optional, contiguous :: lsf2_rr_packed(:, :)
      !> Packed third derivatives
      real(wp), intent(out), optional, contiguous :: lsf3_rrr(:, :)
      !> Packed fourth derivatives
      real(wp), intent(out), optional, contiguous :: lsf4_rrrr(:, :)
      !> Per-column status
      integer, intent(out), optional, contiguous :: status(:)
      !> Shared atom list
      integer, intent(in), optional :: candidate_indices(:)
      !> First failing column
      type(error_type), allocatable, intent(out) :: error

      !> Number of points and list entries
      integer :: n, nc

      n = size(points, 2)
      call lsf_batch_require_outputs(self%name, points, order, lsf0, lsf1_r, lsf2_rr_packed, &
                                     lsf3_rrr, lsf4_rrrr, status)
      if (n == 0) return
      self%prepared_deriv = -1
      if (.not. self%gto%has_density()) then
         call fatal_error(error, "evaluate_batch failed at point 1: no density matrix installed ("// &
                          self%name//")")
         return
      end if

      if (allocated(self%bat_perm)) then
         if (size(self%bat_perm) < n) deallocate (self%bat_perm, self%bat_key)
      end if
      if (.not. allocated(self%bat_perm)) allocate (self%bat_perm(n), self%bat_key(n))
      call tile_order(points, moist_iso_gto_tile_width(order), self%bat_key(1:n), self%bat_perm(1:n))

      if (present(candidate_indices)) then
         if (allocated(self%bat_cand)) then
            if (size(self%bat_cand) < size(candidate_indices)) deallocate (self%bat_cand)
         end if
         if (.not. allocated(self%bat_cand)) allocate (self%bat_cand(size(candidate_indices)))
         call lsf_sorted_candidates(candidate_indices, self%bat_cand, nc, self%ncenters)
         call batch_tiles(self, points, order, lsf0, lsf1_r, lsf2_rr_packed, lsf3_rrr, &
                          lsf4_rrrr, status, error, self%bat_cand(1:nc))
      else
         call batch_tiles(self, points, order, lsf0, lsf1_r, lsf2_rr_packed, lsf3_rrr, &
                          lsf4_rrrr, status, error)
      end if
   end subroutine lsf_evaluate_batch

   !> Tile loop of [[lsf_evaluate_batch]] over one optional shared list
   !>
   !> - lift expressions those of `prepare`, packing that of the looping default
   !>
   !> @param[inout] self      LSF instance
   !> @param[in]    points    evaluation points (3, n)
   !> @param[in]    order     highest spatial derivative to write (0..4)
   !> @param[out]   lsf0      values (n)
   !> @param[out]   lsf1_r    gradients (3, n)
   !> @param[out]   lsf2_rr_packed packed Hessians (6, n)
   !> @param[out]   lsf3_rrr  packed third derivatives (10, n)
   !> @param[out]   lsf4_rrrr packed fourth derivatives (15, n)
   !> @param[out]   status    per-column status (n), 2 on a core-fill refusal
   !> @param[out]   error     first failing column
   !> @param[in]    cand      ascending shared atom list (optional)
   subroutine batch_tiles(self, points, order, lsf0, lsf1_r, lsf2_rr_packed, lsf3_rrr, &
                          lsf4_rrrr, status, error, cand)
      !> LSF instance
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      !> Evaluation points
      real(wp), intent(in), contiguous :: points(:, :)
      !> Highest spatial derivative to write
      integer, intent(in) :: order
      !> Values
      real(wp), intent(inout), contiguous :: lsf0(:)
      !> Gradients
      real(wp), intent(inout), optional, contiguous :: lsf1_r(:, :)
      !> Packed Hessians
      real(wp), intent(inout), optional, contiguous :: lsf2_rr_packed(:, :)
      !> Packed third derivatives
      real(wp), intent(inout), optional, contiguous :: lsf3_rrr(:, :)
      !> Packed fourth derivatives
      real(wp), intent(inout), optional, contiguous :: lsf4_rrrr(:, :)
      !> Per-column status
      integer, intent(inout), optional, contiguous :: status(:)
      !> First failing column
      type(error_type), allocatable, intent(inout) :: error
      !> Ascending shared atom list
      integer, intent(in), optional :: cand(:)

      !> Density and its derivatives over one tile
      real(wp) :: rho(max_tile), drho(ndim, max_tile), d2rho(ndim, ndim, max_tile)
      real(wp) :: d3rho(ndim, ndim, ndim, max_tile), d4rho(ndim, ndim, ndim, ndim, max_tile)
      !> Points of one tile
      real(wp) :: tp(ndim, max_tile)
      !> Second, third and fourth level-set derivative of one point before packing
      real(wp) :: t2(ndim, ndim), t3(ndim, ndim, ndim), t4(ndim, ndim, ndim, ndim)
      !> Point count, derivative order evaluated, tile width, bounds and size
      integer :: n, nderiv, width, s0, s1, np
      !> Point and column index, filled atom holding the point, first failing column
      integer :: p, i, inside, fail
      !> Refusal of one column
      type(error_type), allocatable :: refusal
      !> Column label
      character(len=16) :: col

      n = size(points, 2)
      nderiv = max(order, 1)
      width = moist_iso_gto_tile_width(order)
      fail = n + 1
      do s0 = 1, n, width
         s1 = min(n, s0 + width - 1)
         np = s1 - s0 + 1
         do p = 1, np
            tp(:, p) = points(:, self%bat_perm(s0 + p - 1))
         end do
         select case (nderiv)
         case (1)
            call self%gto%eval_tile(tp(:, 1:np), self%tile, rho(1:np), drho(:, 1:np), &
                                    cand_atoms=cand)
         case (2)
            call self%gto%eval_tile(tp(:, 1:np), self%tile, rho(1:np), drho(:, 1:np), &
                                    d2rho=d2rho(:, :, 1:np), cand_atoms=cand)
         case (3)
            call self%gto%eval_tile(tp(:, 1:np), self%tile, rho(1:np), drho(:, 1:np), &
                                    d2rho=d2rho(:, :, 1:np), d3rho=d3rho(:, :, :, 1:np), &
                                    cand_atoms=cand)
         case default
            call self%gto%eval_tile(tp(:, 1:np), self%tile, rho(1:np), drho(:, 1:np), &
                                    d2rho=d2rho(:, :, 1:np), d3rho=d3rho(:, :, :, 1:np), &
                                    d4rho=d4rho(:, :, :, :, 1:np), cand_atoms=cand)
         end select

         do p = 1, np
            i = self%bat_perm(s0 + p - 1)
            ! The fill is added whatever the screening, as in prepare
            call self%core%add(tp(:, p), nderiv, rho(p), drho(:, p), d2rho(:, :, p), &
                               d3rho(:, :, :, p), d4rho(:, :, :, :, p), inside)
            call self%core%check(inside, rho(p), refusal)
            if (allocated(refusal)) then
               if (present(status)) status(i) = 2
               ! Tiles run in Morton order: keep the refusal of the lowest column
               if (i < fail) then
                  fail = i
                  call move_alloc(refusal, error)
               else
                  deallocate (refusal)
               end if
               cycle
            end if

            lsf0(i) = self%param%scale*(self%param%rho_iso - rho(p))
            if (order >= 1) lsf1_r(:, i) = -self%param%scale*drho(:, p)
            if (order >= 2) then
               t2 = -self%param%scale*d2rho(:, :, p)
               call pack_sym_matrix(t2, lsf2_rr_packed(:, i), ndim)
            end if
            if (order >= 3) then
               t3 = -self%param%scale*d3rho(:, :, :, p)
               call pack_sym_tensor3(t3, lsf3_rrr(:, i), ndim)
            end if
            if (order >= 4) then
               t4 = -self%param%scale*d4rho(:, :, :, :, p)
               call pack_sym_tensor4(t4, lsf4_rrrr(:, i), ndim)
            end if
            if (present(status)) status(i) = 0
         end do
      end do
      if (fail <= n) then
         write (col, "(i0)") fail
         error%message = "evaluate_batch failed at point "//trim(col)//": "//error%message
      end if
   end subroutine batch_tiles

   !> Batch columns in tile order: the caller's, or the Morton order of the
   !> cells of a 2^order_bits grid over the bounding box of the finite points
   !>
   !> - Morton order only when it lowers the summed bounding-box diagonal of the
   !>   tiles, so a batch the caller already grouped keeps its tiles
   !> - stable, so a batch inside one cell keeps the caller's order
   !> - a point with a non-finite coordinate goes to cell 0
   !>
   !> @param[in]  points  evaluation points (3, n)
   !> @param[in]  width   points per tile
   !> @param[out] key     Morton cell of each column (n)
   !> @param[out] perm    columns in tile order (n)
   subroutine tile_order(points, width, key, perm)
      !> Evaluation points
      real(wp), intent(in) :: points(:, :)
      !> Points per tile
      integer, intent(in) :: width
      !> Morton cell of each column
      integer, intent(out) :: key(:)
      !> Columns in tile order
      integer, intent(out) :: perm(:)

      !> Cells per axis
      integer, parameter :: ncell = 2**order_bits
      !> Bounding box of the finite points
      real(wp) :: lo(ndim), hi(ndim)
      !> Cell of one point per axis
      integer :: cell(ndim)
      integer :: i, ax, b, n

      n = size(points, 2)
      do i = 1, n
         perm(i) = i
      end do
      if (n <= width) return

      lo = huge(1.0_wp)
      hi = -huge(1.0_wp)
      do i = 1, n
         if (.not. all(ieee_is_finite(points(:, i)))) cycle
         lo = min(lo, points(:, i))
         hi = max(hi, points(:, i))
      end do
      do i = 1, n
         key(i) = 0
         if (.not. all(ieee_is_finite(points(:, i)))) cycle
         do ax = 1, ndim
            cell(ax) = 0
            if (hi(ax) > lo(ax)) then
               cell(ax) = min(ncell - 1, int(real(ncell, wp)*((points(ax, i) - lo(ax))/(hi(ax) - lo(ax)))))
            end if
         end do
         do b = 0, order_bits - 1
            do ax = 1, ndim
               if (btest(cell(ax), b)) key(i) = ibset(key(i), ndim*b + ax - 1)
            end do
         end do
      end do
      call cell_argsort(key, perm)
      if (.not. (tile_spread(points, width, perm) < tile_spread(points, width))) then
         do i = 1, n
            perm(i) = i
         end do
      end if
   end subroutine tile_order

   !> Stable counting argsort of the Morton cells, fixed-size counters
   !>
   !> - `perm(j)` is the column of the j-th point in cell order
   !> - one counter per cell of the `2^order_bits` grid on the stack, nothing
   !>   allocated per batch
   !>
   !> @param[in]  key   Morton cell of each column, 0..ncell**3 - 1 (n)
   !> @param[out] perm  Columns in cell order (n)
   pure subroutine cell_argsort(key, perm)
      !> Morton cell of each column
      integer, intent(in) :: key(:)
      !> Columns in cell order
      integer, intent(out) :: perm(:)

      !> Cells per axis
      integer, parameter :: ncell = 2**order_bits
      !> Scatter offset of each cell
      integer :: offset(0:ncell**ndim - 1)
      integer :: i, j, c, nxt

      offset = 0
      do i = 1, size(key)
         offset(key(i)) = offset(key(i)) + 1
      end do
      c = 0
      do j = 0, ncell**ndim - 1
         nxt = c + offset(j)
         offset(j) = c
         c = nxt
      end do
      do i = 1, size(key)
         offset(key(i)) = offset(key(i)) + 1
         perm(offset(key(i))) = i
      end do
   end subroutine cell_argsort

   !> Summed bounding-box diagonal of the tiles of a column order
   !>
   !> @param[in] points  evaluation points (3, n)
   !> @param[in] width   points per tile
   !> @param[in] perm    columns in tile order (n); the caller's order when absent
   !> @returns   spread  summed tile diagonal in Bohr
   function tile_spread(points, width, perm) result(spread)
      !> Evaluation points
      real(wp), intent(in) :: points(:, :)
      !> Points per tile
      integer, intent(in) :: width
      !> Columns in tile order
      integer, intent(in), optional :: perm(:)
      !> Summed tile diagonal
      real(wp) :: spread

      real(wp) :: lo(ndim), hi(ndim), x(ndim)
      integer :: s0, i, j

      spread = 0.0_wp
      do s0 = 1, size(points, 2), width
         lo = huge(1.0_wp)
         hi = -huge(1.0_wp)
         do i = s0, min(size(points, 2), s0 + width - 1)
            j = i
            if (present(perm)) j = perm(i)
            x = points(:, j)
            if (.not. all(ieee_is_finite(x))) cycle
            lo = min(lo, x)
            hi = max(hi, x)
         end do
         if (all(hi >= lo)) spread = spread + norm2(hi - lo)
      end do
   end function tile_spread

   !> Return zero nuclear-derivative placeholders
   !>
   !> The isodensity surface tracks the density itself, so the level set field's
   !> explicit nuclear derivatives vanish (they are carried by the density)
   !>
   !> @param[in]  self       LSF instance
   !> @param[out] lsf1_rA    Optional nuclear gradient placeholder
   !> @param[out] lsf2_r_rA  Optional mixed second derivative placeholder
   !> @param[out] lsf3_rr_rA Packed mixed third derivative placeholder [6, 3, n_active]
   subroutine lsf_f3_rr_rA(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      real(wp), intent(out) :: lsf3_rr_rA(:, :, :)

      if (present(lsf1_rA)) lsf1_rA(:, :) = 0.0_wp
      if (present(lsf2_r_rA)) lsf2_r_rA(:, :, :) = 0.0_wp
      lsf3_rr_rA(:, :, :) = 0.0_wp
   end subroutine lsf_f3_rr_rA

   !> Return a zero jet-contracted nuclear vector-Jacobian product
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of the value
   !> @param[in]  w1   Adjoint weights of the spatial gradient [3]
   !> @param[in]  w2   Packed adjoint weights of the spatial Hessian [6]
   !> @param[out] res  Contracted nuclear gradient placeholder [3, >= n_active]
   subroutine lsf_vjp_f1_rA(self, w0, w1, w2, res)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(in) :: w0
      real(wp), intent(in) :: w1(3)
      real(wp), intent(in) :: w2(:)
      real(wp), intent(out) :: res(:, :)

      res(:, 1:self%active_count()) = 0.0_wp
   end subroutine lsf_vjp_f1_rA

   !> Batched nuclear vector-Jacobian product, identically zero
   !>
   !> Override of [[lsf_base_vjp_f1_rA_batch]], same contract:
   !>
   !> - zero rows for every column without any evaluation: no GTO evaluator call, no
   !>   `prepare`, the per-point cache (`prepared_deriv`) is left alone
   !> - the shape and `candidate_indices` checks of the default, so the same bad
   !>   arguments `error stop` ([[lsf_vjp_batch_shape_message]],
   !>   [[lsf_sorted_candidates]]); the list only sizes `res`
   !> - never unreached: `status = 0` for every column
   !> - the one per-point failure of this LSF (an ECP core-fill refusal) is not
   !>   raised, as no density is evaluated; a consumer sees it from the
   !>   `evaluate_batch` / `prepare` it runs on the same points
   !> - no density installed: `error` "vjp_f1_rA_batch failed at point 1: ...", no
   !>   column written; `np = 0` returns after the argument checks
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    points            Evaluation points (3, np)
   !> @param[in]    w0                Adjoint weights of the value (np)
   !> @param[in]    w1                Adjoint weights of the gradient (3, np)
   !> @param[in]    w2                Packed symmetric part of the Hessian seed (6, np)
   !> @param[out]   res               Contracted nuclear gradient (3, ncand, np), zero
   !> @param[out]   status            Per-column status (np), always 0 (optional)
   !> @param[out]   error             No density installed
   !> @param[in]    candidate_indices Shared atom list in `candidate_space` (optional)
   subroutine lsf_vjp_f1_rA_batch(self, points, w0, w1, w2, res, status, error, &
                                  candidate_indices)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in), contiguous :: points(:, :)
      real(wp), intent(in), contiguous :: w0(:)
      real(wp), intent(in), contiguous :: w1(:, :)
      real(wp), intent(in), contiguous :: w2(:, :)
      real(wp), intent(out), contiguous :: res(:, :, :)
      integer, intent(out), optional, contiguous :: status(:)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: candidate_indices(:)

      !> Ascending copy of `candidate_indices`, validated only
      integer, allocatable :: sorted(:)
      character(len=:), allocatable :: msg
      !> Point count, list size, validated list length
      integer :: np, ncand, nc

      np = size(points, 2)
      ncand = self%ncenters
      if (present(candidate_indices)) ncand = size(candidate_indices)
      if (present(status)) then
         msg = lsf_vjp_batch_shape_message(self%name, points, w0, w1, w2, res, ncand, &
                                           status)
      else
         msg = lsf_vjp_batch_shape_message(self%name, points, w0, w1, w2, res, ncand)
      end if
      if (len(msg) > 0) error stop msg
      if (present(candidate_indices)) then
         allocate (sorted(ncand))
         call lsf_sorted_candidates(candidate_indices, sorted, nc, self%ncenters)
      end if
      if (np == 0) return
      if (.not. self%gto%has_density()) then
         call fatal_error(error, "vjp_f1_rA_batch failed at point 1: "// &
                          "no density matrix installed ("//self%name//")")
         return
      end if

      res(:, :, 1:np) = 0.0_wp
      if (present(status)) status(1:np) = 0
   end subroutine lsf_vjp_f1_rA_batch

   !> Radial offset (from the atom surface) beyond which no shell of this atom
   !> still contributes at the current screening threshold
   !>
   !> Reports the global shell reach for the current screening threshold, minus
   !> the atom radius, so `radius + screening_offset = max(radius, reach)`
   !> With screening disabled (threshold <= 0) the reach is huge and the candidate
   !> cell list degrades to a full scan -- every atom is a candidate, i.e. exact
   !> evaluation
   !>
   !> @param[in] self   LSF instance
   !> @param[in] radius Atom radius (Bohr)
   !> @returns          Radial offset from the atom surface (Bohr)
   pure function lsf_screening_offset(self, radius) result(d)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      !> Atom radius (Bohr)
      real(wp), intent(in) :: radius
      !> Radial offset from the atom surface (Bohr)
      real(wp) :: d

      d = max(0.0_wp, self%gto%reach(self%screening_threshold) - radius)
   end function lsf_screening_offset

   !> Radius of a ball around the evaluation point free of surface
   !>
   !> Delegates to [[isodensity_exclusion_radius]], which carries the derivation
   !> and the caveats; the two isodensity variants share it because they share
   !> the level set, differing only in where `rho` comes from
   !>
   !> @param[in] self  LSF instance
   !> @param[in] lsf0  LSF value at the evaluation point
   !> @returns   r     Surface-free radius (zero when uncertified)
   pure function lsf_exclusion_radius(self, lsf0) result(r)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      !> LSF value at the evaluation point
      real(wp), intent(in) :: lsf0
      !> Surface-free radius
      real(wp) :: r

      r = isodensity_exclusion_radius(self%param, self%zmax, lsf0, self%core%log_grad_bound())
   end function lsf_exclusion_radius

   !> Chain-rule factor between the level set and the density it is built from
   !>
   !> `S = scale * (rho_iso - rho)`, so `dS/drho = -scale` for the value and,
   !> `scale` and `rho_iso` being constants, for every spatial derivative of it
   !> as well; one factor therefore converts the whole level-set adjoint jet
   !> into a density adjoint jet
   !>
   !> @param[in]  self    LSF instance
   !> @param[out] factor  Chain-rule factor `dS/drho`
   !> @returns            Always `.true.`: this level set is a density
   function lsf_density_adjoint_factor(self, factor) result(available)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      !> Chain-rule factor `dS/drho`
      real(wp), intent(out) :: factor
      !> Whether this LSF is backed by a density
      logical :: available

      factor = -self%param%scale
      available = .true.
   end function lsf_density_adjoint_factor

   !> The density splits into atomic shares
   !>
   !> @param[in] self LSF instance
   pure function lsf_has_atom_partition(self) result(available)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      !> Always available
      logical :: available

      available = .true.
   end function lsf_has_atom_partition

   !> Add the Mulliken shares of the density at the prepared point
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    scale Weight of the point
   !> @param[inout] acc   Per-atom accumulator, user space [ncenters]
   subroutine lsf_atom_partition(self, scale, acc)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      !> Weight of the point
      real(wp), intent(in) :: scale
      !> Per-atom accumulator
      real(wp), intent(inout) :: acc(:)

      !> Density, its unused gradient and its atomic shares
      real(wp) :: rho, drho(ndim), rho_atom(self%ncenters)

      rho_atom = 0.0_wp
      call self%gto%eval(self%point, self%work, rho, drho, rho_atom=rho_atom)
      rho = sum(rho_atom)
      if (rho <= 0.0_wp) return
      acc = acc + rho_atom*(scale/rho)
   end subroutine lsf_atom_partition

   !> Reject projected surface points inside an ECP core fill
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  xyz       Cavity points in Bohr (3, n)
   !> @param[in]  projected Whether each point is a converged surface projection (n)
   !> @param[out] error     Set at the first projected point inside a fill ball
   subroutine lsf_check_surface(self, xyz, projected, error)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(in) :: xyz(:, :)
      logical, intent(in) :: projected(:)
      type(error_type), allocatable, intent(out) :: error

      call self%core%check_surface(xyz, projected, error)
   end subroutine lsf_check_surface

end module moist_cavity_drop_lsf_isodensity_internal
