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
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_isodensity_gto, only: moist_iso_gto_type
   implicit none
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
      !> Highest requested derivative order
      integer :: max_deriv = 0
      !> Per-instance scratch AO-derivative table (ncart, 0:19)
      real(wp), allocatable :: phi(:, :)
      !> Per-instance scratch density-weighted value vector (ncart)
      real(wp), allocatable :: t0(:)
      !> Per-instance scratch density-weighted gradient vectors (ncart, 3)
      real(wp), allocatable :: tm(:, :)
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
      procedure, public :: f0_screened => lsf_f0_screened
      procedure, public :: f012_r_screened => lsf_f012_r_screened
      procedure, public :: f3_rrr_screened => lsf_f3_rrr_screened
      procedure, public :: f3_rr_rA_screened => lsf_f3_rr_rA_screened
      procedure, public :: neighbor_cutoff => lsf_neighbor_cutoff
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

      call self%gto%init(sh_atom, sh_l, sh_nprim, exps, coeffs, error)
      if (allocated(error)) return
      !> The per-instance scratch is sized from ``gto%ncart``, so a re-configured
      !> basis invalidates it. Dropping it here makes ``lsf_prepare_impl`` size it
      !> again; keeping it would let a larger basis write past its end.
      if (allocated(self%phi)) deallocate (self%phi)
      if (allocated(self%t0)) deallocate (self%t0)
      if (allocated(self%tm)) deallocate (self%tm)
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

      self%mol = mol
      self%radii = radii
      self%ncenters = mol%nat
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
   subroutine lsf_prepare(self, point)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)

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
      integer :: nderiv

      nderiv = min(max(self%max_deriv, 1), 3)

      if (.not. allocated(self%phi)) then
         allocate (self%phi(self%gto%ncart, 0:19))
         allocate (self%t0(self%gto%ncart))
         allocate (self%tm(self%gto%ncart, 3))
         allocate (self%act(self%gto%ncart))
      end if

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
      case default
         call self%gto%eval(point, self%phi, self%t0, self%tm, self%act, &
                            rho, drho, d2rho=d2rho, d3rho=d3rho, cand_atoms=cand_atoms)
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
   end subroutine lsf_prepare_impl

   !> Evaluate cached data; candidate lists are ignored for true density LSFs
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point in Bohr
   !> @param[in]    candidate_indices Ignored atom candidates
   subroutine lsf_prepare_subset(self, point, candidate_indices)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: candidate_indices(:)

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
   subroutine lsf_f0_screened(self, val)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out) :: val

      val = self%value
   end subroutine lsf_f0_screened

   !> Return cached LSF value, gradient, and Hessian
   !>
   !> @param[in]  self    LSF instance
   !> @param[out] lsf0    Optional LSF value
   !> @param[out] lsf1_r  Optional spatial gradient
   !> @param[out] lsf2_rr Optional spatial Hessian
   subroutine lsf_f012_r_screened(self, lsf0, lsf1_r, lsf2_rr)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)

      if (present(lsf0)) lsf0 = self%value
      if (present(lsf1_r)) lsf1_r(:) = self%grad(:)
      if (present(lsf2_rr)) then
         call self%require_deriv(2, "f012_r_screened(lsf2_rr)")
         lsf2_rr(:, :) = self%hess(:, :)
      end if
   end subroutine lsf_f012_r_screened

   !> Return cached lower derivatives and the third spatial derivative
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] lsf0     Optional LSF value
   !> @param[out] lsf1_r   Optional spatial gradient
   !> @param[out] lsf2_rr  Optional spatial Hessian
   !> @param[out] lsf3_rrr Spatial third derivative
   subroutine lsf_f3_rrr_screened(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)
      real(wp), allocatable, intent(out) :: lsf3_rrr(:, :, :)

      call self%require_deriv(3, "f3_rrr_screened")
      call self%f012_r_screened(lsf0, lsf1_r, lsf2_rr)
      allocate (lsf3_rrr(ndim, ndim, ndim), source=self%third)
   end subroutine lsf_f3_rrr_screened

   !> Return zero nuclear-derivative placeholders
   !>
   !> The isodensity surface tracks the density itself, so the level set field's
   !> explicit nuclear derivatives vanish (they are carried by the density).
   !>
   !> @param[in]  self       LSF instance
   !> @param[out] lsf1_rA    Optional nuclear gradient placeholder
   !> @param[out] lsf2_r_rA  Optional mixed second derivative placeholder
   !> @param[out] lsf3_rr_rA Mixed third derivative placeholder
   subroutine lsf_f3_rr_rA_screened(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      real(wp), allocatable, intent(out) :: lsf3_rr_rA(:, :, :, :)

      if (present(lsf1_rA)) lsf1_rA(:, :) = 0.0_wp
      if (present(lsf2_r_rA)) lsf2_r_rA(:, :, :) = 0.0_wp
      allocate (lsf3_rr_rA(ndim, ndim, ndim, self%ncenters), source=0.0_wp)
   end subroutine lsf_f3_rr_rA_screened

   !> Radial offset (from the atom surface) the cavity cell grid must span so
   !> that every atom whose shells still contribute at a point is returned as a
   !> candidate for that point
   !>
   !> Reports the global shell reach for the current screening threshold, minus
   !> the atom radius (the cell-grid reach is ``radius + neighbor_cutoff``), so
   !> ``radius + neighbor_cutoff = max(radius, reach)``.  With screening disabled
   !> (threshold <= 0) the reach is huge and the grid degrades to a full scan --
   !> every atom is a candidate, i.e. exact evaluation.
   !>
   !> @param[in] self   LSF instance
   !> @param[in] radius Atom radius (Bohr)
   !> @returns          Radial offset from the atom surface (Bohr)
   pure function lsf_neighbor_cutoff(self, radius) result(d)
      class(moist_cavity_drop_lsf_isodensity_internal_type), intent(in) :: self
      real(wp), intent(in) :: radius

      real(wp) :: d

      d = max(0.0_wp, self%gto%reach(self%screening_threshold) - radius)
   end function lsf_neighbor_cutoff

end module moist_cavity_drop_lsf_isodensity_internal
