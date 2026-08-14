!> Internal isodensity level set function for DROP
!>
!> Internal variant of [[moist_cavity_drop_lsf_isodensity_callback_type]]
!> This LSF owns the cartesian-monomial Gaussian basis and the density matrix
!> and evaluates the level set internally
!>
!> Relies on the GTO code in [[moist_cavity_drop_lsf_isodensity_gto]]
!>
!> Level set function
!>
!>    S(r) = scale * (rho_iso - rho(r))
!>
!> The host code has to provide
!> - Basis set (once)
!> - Density matrix (once per SCF step)
!>
!> The isodensity surface is a level set of the density itself -> level set function carries no
!> *explicit* nuclear-position dependence for the cavity chain rule
!> -> The mixed spatial/nuclear derivatives therefore vanish here
module moist_cavity_drop_lsf_isodensity_internal
   use mctc_env, only: error_type
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, &
                                         lsf_base_update, lsf_candidate_space_user
   use moist_cavity_drop_lsf_isodensity_gto, only: moist_iso_gto_type, moist_iso_gto_nslot
   implicit none (type, external)
   private

   integer, parameter :: ndim = 3

   public :: moist_cavity_drop_lsf_isodensity_internal_type

   !> Isodensity LSF backed by an internal cartesian-monomial GTO evaluator.
   type, extends(moist_cavity_drop_lsf_type) :: moist_cavity_drop_lsf_isodensity_internal_type
      !> Cartesian-monomial Gaussian basis and density matrix
      type(moist_iso_gto_type) :: gto
      !> Constant multiplier applied to the level set value and derivatives
      real(wp) :: scale = 1.0_wp
      !> Density isovalue defining the surface
      real(wp) :: rho_iso = 0.0_wp
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
      !> Per-instance scratch AO-derivative table (ncart, 0:nslot-1), sized to the
      !> highest order this instance has been asked for so far
      real(wp), allocatable :: phi(:, :)
      !> Per-instance scratch density-weighted value vector (ncart)
      real(wp), allocatable :: t0(:)
      !> Per-instance scratch density-weighted gradient vectors (ncart, 3)
      real(wp), allocatable :: tm(:, :)
      !> Per-instance scratch density-weighted Hessian vectors (ncart, 6)
      real(wp), allocatable :: tmm(:, :)
      !> Per-instance scratch active-component index list (ncart)
      integer, allocatable :: act(:)
   contains
      procedure, public :: new => lsf_new
      procedure, public :: set_density => lsf_set_density
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
      procedure, public :: screening_offset => lsf_screening_offset
   end type moist_cavity_drop_lsf_isodensity_internal_type

contains

   !> Configure the basis and level set parameters
   !>
   !> @param[inout] self     LSF instance
   !> @param[in]    sh_atom  Per-shell owner atom index (1-based)
   !> @param[in]    sh_l     Per-shell angular momentum
   !> @param[in]    sh_nprim Per-shell primitive count
   !> @param[in]    exps     Primitive exponents
   !> @param[in]    coeffs   Primitive contraction coefficients (host-normalized)
   !> @param[in]    rho_iso  Density isovalue defining the surface
   !> @param[in]    scale    Constant level set multiplier
   !> @param[out]   error    Set on invalid basis input
   subroutine lsf_new(self, sh_atom, sh_l, sh_nprim, exps, coeffs, rho_iso, scale, error)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      integer, intent(in) :: sh_atom(:)
      integer, intent(in) :: sh_l(:)
      integer, intent(in) :: sh_nprim(:)
      real(wp), intent(in) :: exps(:)
      real(wp), intent(in) :: coeffs(:)
      real(wp), intent(in) :: rho_iso
      real(wp), intent(in), optional :: scale
      type(error_type), allocatable, intent(out) :: error

      !> Candidate ids index this LSF's per-atom GTO shells
      self%candidate_space = lsf_candidate_space_user

      call self%gto%init(sh_atom, sh_l, sh_nprim, exps, coeffs, error)
      if (allocated(error)) return
      !> The per-instance scratch is sized from ``gto%ncart``, so a re-configured
      !> basis invalidates it. Dropping it here makes ``lsf_prepare_impl`` size it
      !> again; keeping it would let a larger basis write past its end.
      if (allocated(self%phi)) deallocate (self%phi)
      if (allocated(self%t0)) deallocate (self%t0)
      if (allocated(self%tm)) deallocate (self%tm)
      if (allocated(self%tmm)) deallocate (self%tmm)
      if (allocated(self%act)) deallocate (self%act)
      self%rho_iso = rho_iso
      if (present(scale)) self%scale = scale
   end subroutine lsf_new

   !> Install the cartesian-monomial density matrix for the current SCF step
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    dcart Density matrix in the cartesian-monomial basis
   !> @param[out]   error Set on a size mismatch
   subroutine lsf_set_density(self, dcart, error)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: dcart(:, :)
      type(error_type), allocatable, intent(out) :: error

      call self%gto%set_density(dcart, error)
   end subroutine lsf_set_density

   !> Bind molecular geometry and refresh the shell centers.
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
   end subroutine lsf_update

   !> Evaluate and cache the level set at one point via the internal GTO evaluator
   !>
   !> Only the derivative orders demanded by ``max_deriv`` are computed
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    point Evaluation point in Bohr
   subroutine lsf_prepare(self, point, error)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      type(error_type), allocatable, intent(out) :: error
      call lsf_prepare_impl(self, point)
   end subroutine lsf_prepare

   !> Shared prepare body: evaluate the level set at ``point`` and cache it
   !>
   !> When ``cand_atoms`` is present only the shells owned by those atoms are computed
   !>
   !> @param[inout] self       LSF instance
   !> @param[in]    point      Evaluation point in Bohr
   !> @param[in]    cand_atoms Optional candidate owner-atom list
   subroutine lsf_prepare_impl(self, point, cand_atoms)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      integer, intent(in), optional :: cand_atoms(:)

      real(wp) :: rho, drho(ndim), d2rho(ndim, ndim), d3rho(ndim, ndim, ndim)
      real(wp) :: d4rho(ndim, ndim, ndim, ndim)
      integer :: nderiv, nslot

      nderiv = min(max(self%max_deriv, 1), 4)
      nslot = moist_iso_gto_nslot(nderiv)

      if (.not. allocated(self%phi)) then
         allocate (self%phi(self%gto%ncart, 0:nslot - 1))
         allocate (self%t0(self%gto%ncart))
         allocate (self%tm(self%gto%ncart, 3))
         allocate (self%act(self%gto%ncart))
      else if (size(self%phi, 2) < nslot) then
         deallocate (self%phi)
         allocate (self%phi(self%gto%ncart, 0:nslot - 1))
      end if
      if (nderiv >= 4 .and. .not. allocated(self%tmm)) allocate (self%tmm(self%gto%ncart, 6))

      ! The evaluator derives the derivative order from the outputs it is handed,
      ! so the value+gradient phase passes neither the Hessian nor the third
      ! derivative and never pays for them.  An absent ``cand_atoms`` stays absent
      ! when forwarded, so both screening modes share one call per order
      select case (nderiv)
      case (:1)
         call self%gto%eval(point, self%phi, self%t0, self%tm, self%act, &
                            rho, drho, cand_atoms=cand_atoms)
      case (2)
         call self%gto%eval(point, self%phi, self%t0, self%tm, self%act, &
                            rho, drho, d2rho=d2rho, cand_atoms=cand_atoms)
      case (3)
         call self%gto%eval(point, self%phi, self%t0, self%tm, self%act, &
                            rho, drho, d2rho=d2rho, d3rho=d3rho, cand_atoms=cand_atoms)
      case default
         call self%gto%eval(point, self%phi, self%t0, self%tm, self%act, &
                            rho, drho, d2rho=d2rho, d3rho=d3rho, cand_atoms=cand_atoms, &
                            d4rho=d4rho, tmm=self%tmm)
      end select

      self%point = point
      ! Record what this point's cache actually holds, so an accessor asked for
      ! a higher order aborts instead of returning the zeros below
      self%prepared_deriv = nderiv
      self%value = self%scale*(self%rho_iso - rho)
      self%grad = -self%scale*drho
      if (nderiv >= 2) then
         self%hess = -self%scale*d2rho
      else
         self%hess = 0.0_wp
      end if
      if (nderiv >= 3) then
         self%third = -self%scale*d3rho
      else
         self%third = 0.0_wp
      end if
      if (nderiv >= 4) self%fourth = -self%scale*d4rho
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
      !> The candidate atoms come from the cavity's molecular cell grid, whose
      !> per-atom reach is sized by lsf_neighbor_cutoff to the shell reach, so no
      !> contributing shell is missed.  Forward them straight to the evaluator.
      call lsf_prepare_impl(self, point, cand_atoms=candidate_indices)
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

   !> Number of active atoms. True-density LSFs are not atom screened
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

   !> Active atom lookup. Undefined for zero active atoms, returns zero sentinel
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

   !> Return cached LSF value.
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
   !> explicit nuclear derivatives vanish (they are carried by the density).
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

   !> Radial offset (from the atom surface) beyond which no shell of this atom
   !> still contributes at the current screening threshold
   !>
   !> Reports the global shell reach for the current screening threshold, minus
   !> the atom radius, so ``radius + screening_offset = max(radius, reach)``.
   !> With screening disabled (threshold <= 0) the reach is huge and the cavity
   !> cell grid degrades to a full scan -- every atom is a candidate, i.e. exact
   !> evaluation.
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

end module moist_cavity_drop_lsf_isodensity_internal
