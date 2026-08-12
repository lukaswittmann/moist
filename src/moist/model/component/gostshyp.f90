!> GOSTSHYP (Gaussians On Surface Tesserae Simulating HYdrostatic Pressure)
!>
!> An unnormalized Gaussian is placed on every cavity grid point and its
!> amplitude is fixed by the constraint that the force the Gaussian exerts
!> on the electron density matches the applied pressure times the area,
!>
!>    w_i = pi ln2 / a_i
!>    G_i = exp(-w_i |r - r_i|^2)
!>    p_i = p_inp a_i / ftilde_i
!>    E   = sum_i p_i gtilde_i
!>
!> where `gtilde_i = <G_i>` and `ftilde_i = <n_i . grad_r G_i>` are traces
!> against the solute density
!>
!>
!> Host supplies
!> * the AO three-center integrals
!> * the two higher Gaussian moments, in `coupling%gauss_gt/pt/mt/rt`
!>
!> Only the relative s/p/d/f angular normalization matters for the moments;
!> the per-point Gaussian normalization cancels between energy and amplitudes
module moist_model_component_gostshyp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_type, only: solvation_model_component, cavity_type, coupling_type, potential_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type

   implicit none (type, external)
   private

   public :: solvation_model_component_gostshyp, new_component_gostshyp

   !> Relative floor on `|ftilde_i|` below which a grid point is inactive
   !>
   !> Every GOSTSHYP quantity is the ratio of two exponentially small numbers,
   !> so grid points that have "left the density" do not contribute
   !>
   !> FIXME: This is a pragmatic solution for now
   real(wp), parameter :: overlap_floor = 1.0e-9_wp

   !> Pi times ln 2, the numerator of the Gaussian width
   real(wp), parameter :: pi_ln2 = 3.14159265358979323846_wp*0.69314718055994530942_wp

   !> GOSTSHYP hydrostatic pressure contribution
   type, extends(solvation_model_component) :: solvation_model_component_gostshyp
      !> Applied hydrostatic pressure in atomic units, Hartree/bohr**3
      real(wp) :: pressure = 0.0_wp
   contains
      procedure :: update => gostshyp_update
      procedure :: get_energy => gostshyp_get_energy
      procedure :: get_potential => gostshyp_get_potential
      procedure :: get_gradient => gostshyp_get_gradient
      procedure :: get_surface_weights => gostshyp_get_surface_weights
   end type solvation_model_component_gostshyp

contains

   !* ================================================================================= *!
   !*                                    Constructor                                    *!
   !* ================================================================================= *!

   !> Construct a GOSTSHYP pressure component
   !>
   !> @param[out] self     Component instance
   !> @param[in]  pressure Applied pressure in Hartree/bohr**3
   subroutine new_component_gostshyp(self, pressure)
      !> Component instance
      type(solvation_model_component_gostshyp), intent(out) :: self
      !> Applied pressure
      real(wp), intent(in) :: pressure

      self%name = "GOSTSHYP"
      self%pressure = pressure

   end subroutine new_component_gostshyp

   !> Bind the current molecular structure
   !>
   !> @param[inout] self   Component instance
   !> @param[in]    mol    Molecular structure
   !> @param[inout] cavity Live model cavity
   !> @param[out]   error  Error handling
   subroutine gostshyp_update(self, mol, cavity, error)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      self%mol_solu = mol
      if (.not. allocated(cavity%a) .or. .not. allocated(cavity%xyz) .or. &
          .not. allocated(cavity%normal0)) then
         call fatal_error(error, "GOSTSHYP component requires an updated cavity")
      end if

   end subroutine gostshyp_update

   !* ================================================================================= *!
   !*                              Shared state evaluation                              *!
   !* ================================================================================= *!

   !> Evaluate the per-grid-point amplitudes shared by energy, potential and weights
   !>
   !> Returns the Gaussian widths, the normal-projected gradient trace, and the
   !> two amplitudes, all already zeroed on the inactive grid points and already
   !> scaled by the component's linear `scale` factor.
   !>
   !> @param[in]  self     Component instance
   !> @param[in]  coupling Host coupling data carrying the Gaussian moments
   !> @param[in]  cavity   Live model cavity
   !> @param[out] omega    Gaussian widths, bohr**-2 (ngrid)
   !> @param[out] ftilde   Normal-projected Gaussian gradient trace (ngrid)
   !> @param[out] alpha    Amplitude conjugate to `g_uv,i` (ngrid)
   !> @param[out] beta      Amplitude conjugate to `-f_uv,i` (ngrid)
   !> @param[out] error     Error handling
   !> @param[out] ninactive Grid points switched off by the floor (optional)
   subroutine gostshyp_amplitudes(self, coupling, cavity, omega, ftilde, alpha, beta, error, &
                                  ninactive)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(in) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(in) :: cavity
      !> Gaussian widths and the traces built from the supplied moments
      real(wp), allocatable, intent(out) :: omega(:), ftilde(:), alpha(:), beta(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Grid points switched off by the activity floor
      integer, intent(out), optional :: ninactive

      !> Grid-point index
      integer :: igrid
      !> Activity threshold on the normal-projected trace
      real(wp) :: floor
      !> Candidate amplitudes, kept only once they are known to be representable
      real(wp) :: alpha_i, beta_i
      !> Whether this grid point survives both the floor and the divisions
      logical :: active

      call require_moments(coupling, cavity, error)
      if (allocated(error)) return

      allocate (omega(cavity%ngrid), ftilde(cavity%ngrid))
      allocate (alpha(cavity%ngrid), source=0.0_wp)
      allocate (beta(cavity%ngrid), source=0.0_wp)

      do igrid = 1, cavity%ngrid
         !> A degenerate zero-area grid point is inert rather than infinite.
         if (cavity%a(igrid) > 0.0_wp) then
            omega(igrid) = pi_ln2/cavity%a(igrid)
         else
            omega(igrid) = 0.0_wp
         end if
         !> ftilde is derived rather than supplied: it must share one normal
         !> convention with the derivatives built from the same moments.
         ftilde(igrid) = -2.0_wp*omega(igrid) &
            & *dot_product(cavity%normal0(:, igrid), coupling%gauss_pt(:, igrid))
      end do

      floor = overlap_floor*maxval(abs(ftilde))
      if (present(ninactive)) ninactive = 0
      do igrid = 1, cavity%ngrid
         active = abs(ftilde(igrid)) > floor
         if (active) then
            alpha_i = self%scale*self%pressure*cavity%a(igrid)/ftilde(igrid)
            beta_i = coupling%gauss_gt(igrid)*alpha_i/ftilde(igrid)
            !> The floor is *relative*, so it says nothing about the absolute
            !> size of `ftilde`: a grid uniformly down at the denormals clears it
            !> intact and the two divisions above then leave the reals. No SCF
            !> density gets within 300 decades of that, but the moments arrive
            !> over the C API from an arbitrary host, and an infinity here would
            !> spread silently through the energy, the host's Fock matrix and
            !> every surface weight. A point that cannot be divided is treated
            !> like a point that has left the density.
            active = ieee_is_finite(alpha_i) .and. ieee_is_finite(beta_i)
         end if

         if (active) then
            alpha(igrid) = alpha_i
            beta(igrid) = beta_i
         else if (present(ninactive)) then
            ninactive = ninactive + 1
         end if
      end do

   end subroutine gostshyp_amplitudes

   !> Verify that the host supplied every Gaussian moment at the cavity size
   !>
   !> @param[in]  coupling Host coupling data
   !> @param[in]  cavity   Live model cavity
   !> @param[out] error    Error handling
   subroutine require_moments(coupling, cavity, error)
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(in) :: cavity
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(coupling%gauss_gt) .or. .not. allocated(coupling%gauss_pt) .or. &
          .not. allocated(coupling%gauss_mt) .or. .not. allocated(coupling%gauss_rt)) then
         call fatal_error(error, "GOSTSHYP requires the host to supply Gaussian density "// &
            & "moments for the current cavity")
         return
      end if

      if (size(coupling%gauss_gt) /= cavity%ngrid .or. &
          any(shape(coupling%gauss_pt) /= [3, cavity%ngrid]) .or. &
          any(shape(coupling%gauss_mt) /= [3, 3, cavity%ngrid]) .or. &
          any(shape(coupling%gauss_rt) /= [3, cavity%ngrid])) then
         call fatal_error(error, "GOSTSHYP Gaussian moments do not match the cavity grid; "// &
            & "they must be rebuilt after every cavity update")
      end if

      !> Only the shapes are checked. The Gaussians sit *on* the grid points, so
      !> every cavity update invalidates every moment -- but a geometry step
      !> usually keeps the point count, so moments from the previous surface fit
      !> perfectly and cannot be told apart here. Supplying current moments after
      !> every update is the host's responsibility.

   end subroutine require_moments

   !* ================================================================================= *!
   !*                             Energy, potential, weights                            *!
   !* ================================================================================= *!

   !> Add the GOSTSHYP pressure energy
   !>
   !> @param[inout] self     Component instance
   !> @param[in]    coupling Host coupling data
   !> @param[inout] cavity   Live model cavity
   !> @param[inout] energy   Energy accumulator
   !> @param[out]   error    Error handling
   subroutine gostshyp_get_energy(self, coupling, cavity, energy, error)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Energy accumulator
      real(wp), intent(inout) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Gaussian widths, traces and amplitudes
      real(wp), allocatable :: omega(:), ftilde(:), alpha(:), beta(:)
      !> Grid points switched off by the activity floor
      integer :: ninactive
      !> Diagnostic line
      character(len=80) :: report

      if (self%pressure == 0.0_wp .or. self%scale == 0.0_wp) return
      call gostshyp_amplitudes(self, coupling, cavity, omega, ftilde, alpha, beta, error, &
         & ninactive=ninactive)
      if (allocated(error)) return

      energy = energy + dot_product(alpha, coupling%gauss_gt)

      !> Reported unconditionally rather than above a threshold: a sizeable
      !> inactive fraction is *normal* (15% for fluoroacetate/STO-3G at 50 GPa),
      !> so any threshold loose enough to stay quiet would also stay quiet for
      !> the failure worth catching -- a systematically wrong `gauss_pt`, which
      !> shrinks every `ftilde` and pushes points under the floor. The number is
      !> the diagnostic; what counts as too many is the reader's call.
      !> Reported here rather than in `gostshyp_amplitudes` so one energy
      !> evaluation produces one line, not three.
      if (associated(self%ctx) .and. ninactive > 0) then
         write (report, '(a,i0,a,i0,a)') &
            & "GOSTSHYP: ", ninactive, " of ", cavity%ngrid, &
            & " grid points below the density-overlap floor"
         call self%ctx%message(trim(report), level=2)
      end if

   end subroutine gostshyp_get_energy

   !> Hand the host the amplitudes conjugate to its Gaussian integral blocks
   !>
   !> The host rebuilds its Fock contribution as
   !> `F_uv += sum_i [w_gauss_g(i) g_uv,i + w_gauss_f(i) f_uv,i]`; both signs
   !> are folded in here.
   !>
   !> @param[inout] self      Component instance
   !> @param[in]    coupling  Host coupling data
   !> @param[inout] cavity    Live model cavity
   !> @param[inout] potential Potential accumulator
   !> @param[out]   error     Error handling
   subroutine gostshyp_get_potential(self, coupling, cavity, potential, error)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Potential accumulator
      type(potential_type), intent(inout) :: potential
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Gaussian widths, traces and amplitudes
      real(wp), allocatable :: omega(:), ftilde(:), alpha(:), beta(:)

      if (self%pressure == 0.0_wp .or. self%scale == 0.0_wp) return
      call gostshyp_amplitudes(self, coupling, cavity, omega, ftilde, alpha, beta, error)
      if (allocated(error)) return

      if (.not. allocated(potential%w_gauss_g)) then
         allocate (potential%w_gauss_g(cavity%ngrid), source=0.0_wp)
      end if
      if (.not. allocated(potential%w_gauss_f)) then
         allocate (potential%w_gauss_f(cavity%ngrid), source=0.0_wp)
      end if

      potential%w_gauss_g = potential%w_gauss_g + alpha
      potential%w_gauss_f = potential%w_gauss_f - beta

   end subroutine gostshyp_get_potential

   !> Add the GOSTSHYP surface adjoints
   !>
   !> With `G = exp(-w |r - C|^2)` and the outward normal held fixed, the
   !> supplied moments give every parameter derivative of the two traces,
   !>
   !>    dgtilde/dC_a = 2 w Pt_a          dftilde/dC_b = 2 w n_b gt - 4 w^2 (n.Mt)_b
   !>    dgtilde/dw   = -tr(Mt)           dftilde/dw   = -2 (n.Pt) + 2 w (n.Rt)
   !>    dftilde/dn   = -2 w Pt
   !>
   !> The area enters twice: explicitly through the amplitude `p_i`, and through
   !> the Gaussian width `w_i = pi ln2 / a_i`, whence the `-w_i/a_i` chain
   !> factor on the width route. The switching factor carries no dependence at
   !> all -- `w_f` is exactly zero -- because the Gaussian width is the only
   !> route by which the area reaches the level set.
   !>
   !> These same weights serve the nuclear gradient: the base class points
   !> `get_gradient_surface_weights` here, and the cavity contracts them once in
   !> reverse mode.
   !>
   !> @param[inout] self     Component instance
   !> @param[in]    coupling Host coupling data
   !> @param[in]    cavity   Live model cavity
   !> @param[inout] acc      Surface-adjoint accumulator
   !> @param[out]   error    Error handling
   subroutine gostshyp_get_surface_weights(self, coupling, cavity, acc, error)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(in) :: cavity
      !> Surface accumulator
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Gaussian widths, traces and amplitudes
      real(wp), allocatable :: omega(:), ftilde(:), alpha(:), beta(:)
      !> Adjoints of the grid-point areas, positions and normals
      real(wp), allocatable :: w_a(:), w_xyz(:, :), w_n(:, :)
      !> Parameter derivatives of the two traces at one grid point
      real(wp) :: dgdr(3), dfdr(3), dgdw, dfdw
      !> Normal projections of the supplied moments at one grid point
      real(wp) :: n_pt, n_mt(3), n_rt
      !> Grid-point index
      integer :: igrid

      if (self%pressure == 0.0_wp .or. self%scale == 0.0_wp) return
      call gostshyp_amplitudes(self, coupling, cavity, omega, ftilde, alpha, beta, error)
      if (allocated(error)) return

      allocate (w_a(cavity%ngrid), source=0.0_wp)
      allocate (w_xyz(3, cavity%ngrid), source=0.0_wp)
      allocate (w_n(3, cavity%ngrid), source=0.0_wp)

      do igrid = 1, cavity%ngrid
         !> Inactive grid points were zeroed in both amplitudes, and the width
         !> route below would divide by a trace that carries only round-off.
         if (alpha(igrid) == 0.0_wp) cycle

         n_pt = dot_product(cavity%normal0(:, igrid), coupling%gauss_pt(:, igrid))
         n_mt = matmul(cavity%normal0(:, igrid), coupling%gauss_mt(:, :, igrid))
         n_rt = dot_product(cavity%normal0(:, igrid), coupling%gauss_rt(:, igrid))

         dgdr = 2.0_wp*omega(igrid)*coupling%gauss_pt(:, igrid)
         dfdr = 2.0_wp*omega(igrid)*cavity%normal0(:, igrid)*coupling%gauss_gt(igrid) &
            & - 4.0_wp*omega(igrid)**2*n_mt
         dgdw = -(coupling%gauss_mt(1, 1, igrid) + coupling%gauss_mt(2, 2, igrid) &
            &     + coupling%gauss_mt(3, 3, igrid))
         dfdw = -2.0_wp*n_pt + 2.0_wp*omega(igrid)*n_rt

         w_xyz(:, igrid) = alpha(igrid)*dgdr - beta(igrid)*dfdr
         !> Only ftilde depends on the normal, through ftilde = n . (-2 w Pt).
         w_n(:, igrid) = 2.0_wp*omega(igrid)*beta(igrid)*coupling%gauss_pt(:, igrid)
         w_a(igrid) = self%scale*self%pressure*coupling%gauss_gt(igrid)/ftilde(igrid) &
            & + (alpha(igrid)*dgdw - beta(igrid)*dfdw)*(-omega(igrid)/cavity%a(igrid))
      end do

      call acc%add_surface_weights(error, w_a=w_a, w_xyz=w_xyz, w_n=w_n)

   end subroutine gostshyp_get_surface_weights

   !* ================================================================================= *!
   !*                                  Nuclear gradient                                 *!
   !* ================================================================================= *!

   !> Forward-mode nuclear gradient, which GOSTSHYP does not provide
   !>
   !> The energy reaches the nuclei only through cavity-surface quantities, all
   !> of which `get_surface_weights` already states. The reverse-mode path
   !> contracts those once against the cavity's own nuclear derivatives; a
   !> forward-mode implementation would be a second, independently maintained
   !> derivation of the same numbers, so it is refused rather than approximated.
   !>
   !> @param[inout] self     Component instance
   !> @param[in]    coupling Host coupling data, unused
   !> @param[inout] cavity   Live model cavity, unused
   !> @param[inout] gradient Nuclear-gradient accumulator, unchanged
   !> @param[out]   error    Error handling
   subroutine gostshyp_get_gradient(self, coupling, cavity, gradient, error)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (self%pressure == 0.0_wp .or. self%scale == 0.0_wp) return
      call fatal_error(error, "GOSTSHYP has no forward-mode nuclear gradient; use the "// &
         & "reverse-mode surface path (force_forward_gradient must stay disabled)")

   end subroutine gostshyp_get_gradient

end module moist_model_component_gostshyp
