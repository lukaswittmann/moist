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
   use mctc_env, only: error_type
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, &
                                         lsf_base_update, lsf_candidate_space_user
   use moist_cavity_drop_lsf_isodensity_gto, only: moist_iso_gto_type, moist_iso_gto_work_type
   use moist_cavity_drop_lsf_isodensity_param, only: &
      moist_cavity_drop_lsf_isodensity_param_type, isodensity_exclusion_radius, &
      isodensity_max_nuclear_charge
   use moist_cavity_drop_lsf_isodensity_core, only: moist_iso_core_fill_type, iso_check_nuclear_charges
   use moist_channels_coupling, only: coupling_type, coupling_view_type, coupling_register, &
      & request_require, moist_phase_build, nuclear_charge_request_type, density_matrix_request_type
   implicit none(type, external)
   private

   integer, parameter :: ndim = 3

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
      procedure, public :: f3_rr_rA => lsf_f3_rr_rA
      procedure, public :: vjp_f1_rA => lsf_vjp_f1_rA
      procedure, public :: screening_offset => lsf_screening_offset
      procedure, public :: exclusion_radius => lsf_exclusion_radius
      procedure, public :: density_adjoint_factor => lsf_density_adjoint_factor
      procedure, public :: check_surface => lsf_check_surface
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

      self%param = moist_cavity_drop_lsf_isodensity_param_type()
      if (present(param)) self%param = param
      call self%gto%init(sh_atom, sh_l, sh_nprim, exps, coeffs, error, basis=self%param%basis)
      if (allocated(error)) return
      !> Scratch sized for the previous basis; resized on first use
      self%work = moist_iso_gto_work_type()
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
   end subroutine lsf_update

   !> Evaluate and cache the level set at one point via the internal GTO evaluator
   !>
   !> Only the derivative orders demanded by `max_deriv` are computed
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    point Evaluation point in Bohr
   subroutine lsf_prepare(self, point, error)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      type(error_type), allocatable, intent(out) :: error
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

   !> Evaluate cached data; candidate lists are ignored for true density LSFs
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point in Bohr
   !> @param[in]    candidate_indices Ignored atom candidates
   subroutine lsf_prepare_subset(self, point, candidate_indices, error)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: candidate_indices(:)
      type(error_type), allocatable, intent(out) :: error

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
   !> @param[out] lsf2_rr Optional spatial Hessian
   subroutine lsf_f012_r(self, lsf0, lsf1_r, lsf2_rr)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)

      if (present(lsf0)) lsf0 = self%value
      if (present(lsf1_r)) lsf1_r(:) = self%grad(:)
      if (present(lsf2_rr)) then
         call self%require_deriv(2, "f012_r(lsf2_rr)")
         lsf2_rr(:, :) = self%hess(:, :)
      end if
   end subroutine lsf_f012_r

   !> Return cached lower derivatives and the third spatial derivative
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] lsf0     Optional LSF value
   !> @param[out] lsf1_r   Optional spatial gradient
   !> @param[out] lsf2_rr  Optional spatial Hessian
   !> @param[out] lsf3_rrr Spatial third derivative
   subroutine lsf_f3_rrr(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)
      real(wp), intent(out) :: lsf3_rrr(:, :, :)

      call self%require_deriv(3, "f3_rrr")
      call self%f012_r(lsf0, lsf1_r, lsf2_rr)
      lsf3_rrr(:, :, :) = self%third
   end subroutine lsf_f3_rrr

   !> Return the cached fourth spatial derivative
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf4_rrrr Fourth spatial derivative [3, 3, 3, 3]
   subroutine lsf_f4_rrrr(self, lsf4_rrrr)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rrrr(:, :, :, :)

      call self%require_deriv(4, "f4_rrrr")
      lsf4_rrrr(:, :, :, :) = self%fourth
   end subroutine lsf_f4_rrrr

   !> Return zero nuclear-derivative placeholders
   !>
   !> The isodensity surface tracks the density itself, so the level set field's
   !> explicit nuclear derivatives vanish (they are carried by the density)
   !>
   !> @param[in]  self       LSF instance
   !> @param[out] lsf1_rA    Optional nuclear gradient placeholder
   !> @param[out] lsf2_r_rA  Optional mixed second derivative placeholder
   !> @param[out] lsf3_rr_rA Mixed third derivative placeholder
   subroutine lsf_f3_rr_rA(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      real(wp), intent(out) :: lsf3_rr_rA(:, :, :, :)

      if (present(lsf1_rA)) lsf1_rA(:, :) = 0.0_wp
      if (present(lsf2_r_rA)) lsf2_r_rA(:, :, :) = 0.0_wp
      lsf3_rr_rA(:, :, :, :) = 0.0_wp
   end subroutine lsf_f3_rr_rA

   !> Return a zero jet-contracted nuclear vector-Jacobian product
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of the value
   !> @param[in]  w1   Adjoint weights of the spatial gradient [3]
   !> @param[in]  w2   Adjoint weights of the spatial Hessian [3, 3]
   !> @param[out] res  Contracted nuclear gradient placeholder [3, >= n_active]
   subroutine lsf_vjp_f1_rA(self, w0, w1, w2, res)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(in) :: w0
      real(wp), intent(in) :: w1(3)
      real(wp), intent(in) :: w2(3, 3)
      real(wp), intent(out) :: res(:, :)

      res(:, 1:self%active_count()) = 0.0_wp
   end subroutine lsf_vjp_f1_rA

   !> Radial offset (from the atom surface) beyond which no shell of this atom
   !> still contributes at the current screening threshold
   !>
   !> Reports the global shell reach for the current screening threshold, minus
   !> the atom radius, so `radius + screening_offset = max(radius, reach)`
   !> With screening disabled (threshold <= 0) the reach is huge and the cavity
   !> cell grid degrades to a full scan -- every atom is a candidate, i.e. exact
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
