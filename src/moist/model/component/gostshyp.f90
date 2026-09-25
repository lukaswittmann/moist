!> GOSTSHYP (Gaussians On Surface Tesserae Simulating HYdrostatic Pressure)
!>
!> An unnormalized Gaussian is placed on every cavity grid point and its
!> amplitude is fixed by the constraint that the force the Gaussian exerts
!> on the electron density matches the applied pressure times the area,
!>    w_i = pi ln2 / a_i
!>
!>    G_i = exp(-w_i |r - r_i|^2)
!>    p_i = p_inp a_i / ftilde_i
!>
!>    E   = sum_i p_i gtilde_i
!> where `gtilde_i = <G_i>` and `ftilde_i = <n_i . grad_r G_i>` are traces
!> against the solute density
!>
!>
!> Host supplies
!>
!> * the AO three-center integrals
!> * the Gaussian moments `gt/pt/mt/rt`, answered on the coupling's
!>
!>   `gaussian_moment_request_type`
!>
!> Only the relative s/p/d/f angular normalization matters for the moments;
!> the per-point Gaussian normalization cancels between energy and amplitudes
module moist_model_component_gostshyp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_cavity_type, only: cavity_type
   use moist_model_type, only: solvation_model_component_type
   use moist_channels_response, only: response_type, gostshyp_amplitude_response_type, &
      & response_accumulate
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & gaussian_moment_request_type, moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, coupling_register, request_require
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type

   implicit none(type, external)
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
   real(wp), parameter :: pi_ln2 = acos(-1.0_wp)*log(2.0_wp)

   !> GOSTSHYP hydrostatic pressure contribution
   type, extends(solvation_model_component_type) :: solvation_model_component_gostshyp
      !> Applied hydrostatic pressure in atomic units, Hartree/bohr**3
      real(wp) :: pressure = 0.0_wp
   contains
      procedure :: update => gostshyp_update
      procedure :: get_energy => gostshyp_get_energy
      procedure :: get_response => gostshyp_get_response
      procedure :: get_gradient => gostshyp_get_gradient
      procedure :: get_surface_weights => gostshyp_get_surface_weights
      !> Declare the Gaussian-moment request
      procedure :: declare_coupling => gostshyp_declare_coupling
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
   !> @param[in,out] self   Component instance
   !> @param[in]    mol    Molecular structure
   !> @param[in,out] cavity Live model cavity
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
   !> scaled by the component's linear `scale` factor
   !>
   !> @param[in]  self     Component instance
   !> @param[in]  cavity   Live model cavity
   !> @param[in]  gt       Host-supplied `<G_i>` (ngrid)
   !> @param[in]  pt       Host-supplied `<(r - r_i) G_i>` (3, ngrid)
   !> @param[out] omega    Gaussian widths, bohr**-2 (ngrid)
   !> @param[out] ftilde   Normal-projected Gaussian gradient trace (ngrid)
   !> @param[out] alpha    Amplitude conjugate to `g_uv,i` (ngrid)
   !> @param[out] beta      Amplitude conjugate to `-f_uv,i` (ngrid)
   !> @param[out] ninactive Grid points switched off by the floor (optional)
   subroutine gostshyp_amplitudes(self, cavity, gt, pt, omega, ftilde, alpha, beta, ninactive)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(in) :: self
      !> Live model cavity
      class(cavity_type), intent(in) :: cavity
      !> Host-supplied Gaussian moments
      real(wp), intent(in) :: gt(:)
      !> Pt
      real(wp), intent(in) :: pt(:, :)
      !> Gaussian widths and the traces built from the supplied moments
      real(wp), allocatable, intent(out) :: omega(:)
      !> Ftilde
      real(wp), allocatable, intent(out) :: ftilde(:)
      !> Alpha
      real(wp), allocatable, intent(out) :: alpha(:)
      !> Beta
      real(wp), allocatable, intent(out) :: beta(:)
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

      allocate (omega(cavity%ngrid), ftilde(cavity%ngrid))
      allocate (alpha(cavity%ngrid), source=0.0_wp)
      allocate (beta(cavity%ngrid), source=0.0_wp)

      do igrid = 1, cavity%ngrid
         !> A degenerate zero-area grid point is inert rather than infinite
         if (cavity%a(igrid) > 0.0_wp) then
            omega(igrid) = pi_ln2/cavity%a(igrid)
         else
            omega(igrid) = 0.0_wp
         end if
         !> ftilde is derived rather than supplied: it must share one normal
         !>
         !> convention with the derivatives built from the same moments
         ftilde(igrid) = -2.0_wp*omega(igrid) &
                        & *dot_product(cavity%normal0(:, igrid), pt(:, igrid))
      end do

      floor = overlap_floor*maxval(abs(ftilde))
      if (present(ninactive)) ninactive = 0
      do igrid = 1, cavity%ngrid
         active = abs(ftilde(igrid)) > floor
         if (active) then
            alpha_i = self%scale*self%pressure*cavity%a(igrid)/ftilde(igrid)
            beta_i = gt(igrid)*alpha_i/ftilde(igrid)
            !> The floor is *relative*, so it says nothing about the absolute
            !>
            !> size of `ftilde`: a grid uniformly down at the denormals clears it
            !> intact and the two divisions above then leave the reals; no SCF
            !> density gets within 300 decades of that, but the moments arrive
            !> over the C API from an arbitrary host, and an infinity here would
            !> spread silently through the energy, the host's Fock matrix and
            !> every surface weight, and a point that cannot be divided is treated
            !> like a point that has left the density
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

   !> Read only the Gaussian moments requested by the caller
   !>
   !> Model updates invalidate all previously supplied moments
   !>
   !> @param[in]  coupling Host coupling data
   !> @param[in]  cavity   Live model cavity
   !> @param[out] gt       `<G_i>` (ngrid)
   !> @param[out] pt       `<(r - r_i) G_i>` (3, ngrid)
   !> @param[out] mt       `<(r - r_i)(r - r_i) G_i>` (3, 3, ngrid)
   !> @param[out] rt       `<(r - r_i) |r - r_i|^2 G_i>` (3, ngrid)
   !> @param[out] error    Error handling
   subroutine read_moments(coupling, cavity, gt, pt, mt, rt, error)
      !> Scoped host answers
      class(coupling_view_type), intent(in) :: coupling
      !> Live cavity
      class(cavity_type), intent(in) :: cavity
      !> Scalar moments
      real(wp), allocatable, intent(out) :: gt(:)
      !> First moments
      real(wp), allocatable, intent(out) :: pt(:, :)
      !> Second moments, only read when needed
      real(wp), allocatable, optional, intent(out) :: mt(:, :, :)
      !> Third moments, only read when needed
      real(wp), allocatable, optional, intent(out) :: rt(:, :)
      !> Read error
      type(error_type), allocatable, intent(out) :: error
      call coupling%read("moments", "gt", gt, error)
      if (allocated(error)) return
      call coupling%read("moments", "pt", pt, error)
      if (allocated(error)) return
      if (present(mt)) then
         call coupling%read("moments", "mt", mt, error)
         if (allocated(error)) return
      end if
      if (present(rt)) then
         call coupling%read("moments", "rt", rt, error)
         if (allocated(error)) return
      end if
      if (size(gt) /= cavity%ngrid .or. any(shape(pt) /= [3, cavity%ngrid])) then
         call fatal_error(error, "GOSTSHYP: moment extents do not match the live cavity")
         return
      end if
      if (present(mt)) then
         if (any(shape(mt) /= [3, 3, cavity%ngrid])) then
            call fatal_error(error, "GOSTSHYP: second-moment extents do not match the live cavity")
            return
         end if
      end if
      if (present(rt)) then
         if (any(shape(rt) /= [3, cavity%ngrid])) then
            call fatal_error(error, "GOSTSHYP: third-moment extents do not match the live cavity")
            return
         end if
      end if
   end subroutine read_moments

   !* ================================================================================= *!
   !*                            Coupling declaration and staging                       *!
   !* ================================================================================= *!

   !> Declare the Gaussian-moment request
   !>
   !> Energy and amplitudes consume gt and pt; geometry derivatives additionally
   !> need mt and rt; disabled pressure or scale declares no host work
   !>
   !> @param[in]    self     Component instance
   !> @param[in]    cavity   Cavity the model is built on, unused
   !> @param[in,out] coupling Coupling being declared
   !> @param[out]   error    Error handling
   subroutine gostshyp_declare_coupling(self, cavity, coupling, error)
      !> Pressure component
      class(solvation_model_component_gostshyp), intent(in) :: self
      !> Live cavity supplying moment centers and areas
      class(cavity_type), intent(in) :: cavity
      !> Component registration context
      type(coupling_type), intent(inout) :: coupling
      !> Registration error
      type(error_type), allocatable, intent(out) :: error
      type(gaussian_moment_request_type) :: moments
      integer :: i
      if (self%pressure == 0.0_wp .or. self%scale == 0.0_wp) return
      if (.not. allocated(cavity%a)) then
         call fatal_error(error, "GOSTSHYP requires cavity areas")
         return
      end if
      allocate (moments%width(cavity%ngrid), source=0.0_wp)
      do i = 1, cavity%ngrid
         if (cavity%a(i) > 0.0_wp) moments%width(i) = pi_ln2/cavity%a(i)
      end do
      call request_require(moments, moist_phase_energy, "gt", error)
      if (allocated(error)) return
      call request_require(moments, moist_phase_energy, "pt", error)
      if (allocated(error)) return
      call request_require(moments, moist_phase_response, "gt", error)
      if (allocated(error)) return
      call request_require(moments, moist_phase_response, "pt", error)
      if (allocated(error)) return
      call request_require(moments, moist_phase_response, "mt", cavity%has_field_dependent_geometry(), error)
      if (allocated(error)) return
      call request_require(moments, moist_phase_response, "rt", cavity%has_field_dependent_geometry(), error)
      if (allocated(error)) return
      call request_require(moments, moist_phase_gradient, "gt", error)
      if (allocated(error)) return
      call request_require(moments, moist_phase_gradient, "pt", error)
      if (allocated(error)) return
      call request_require(moments, moist_phase_gradient, "mt", error)
      if (allocated(error)) return
      call request_require(moments, moist_phase_gradient, "rt", error)
      if (allocated(error)) return
      call coupling_register(coupling, "moments", moments, error)
   end subroutine gostshyp_declare_coupling

   !* ================================================================================= *!
   !*                             Energy, potential, weights                            *!
   !* ================================================================================= *!

   !> Add the GOSTSHYP pressure energy
   !>
   !> @param[in,out] self     Component instance
   !> @param[in]    coupling Host coupling data
   !> @param[in,out] cavity   Live model cavity
   !> @param[in,out] energy   Energy accumulator
   !> @param[out]   error    Error handling
   subroutine gostshyp_get_energy(self, coupling, cavity, energy, error)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(inout) :: self
      !> Host coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Energy accumulator
      real(wp), intent(inout) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Host-supplied Gaussian moments
      real(wp), allocatable :: gt(:), pt(:, :)
      !> Gaussian widths, traces and amplitudes
      real(wp), allocatable :: omega(:), ftilde(:), alpha(:), beta(:)
      !> Grid points switched off by the activity floor
      integer :: ninactive
      !> Diagnostic line
      character(len=80) :: report

      ! A switched-off component asks nothing of the host, so the mandatory
      ! check comes after the short circuit
      if (self%pressure == 0.0_wp .or. self%scale == 0.0_wp) return
      call coupling%check_mandatory(coupling%phase, error)
      if (allocated(error)) return
      call read_moments(coupling, cavity, gt, pt, error=error)
      if (allocated(error)) return
      call gostshyp_amplitudes(self, cavity, gt, pt, omega, ftilde, alpha, beta, &
         & ninactive=ninactive)

      energy = energy + dot_product(alpha, gt)

      !> Reported unconditionally rather than above a threshold: a sizeable
      !>
      !> inactive fraction is *normal* (15% for fluoroacetate/STO-3G at 50 GPa),
      !> so any threshold loose enough to stay quiet would also stay quiet for
      !> the failure worth catching -- a systematically wrong `gostshyp%pt`, which
      !> shrinks every `ftilde` and pushes points under the floor; the number
      !> is the diagnostic, and what counts as too many is the reader's call
      !>
      !> Reported here rather than in `gostshyp_amplitudes`, so one energy
      !> evaluation produces one line, not three
      if (associated(self%ctx) .and. ninactive > 0) then
         write (report, "(a,i0,a,i0,a)") &
            & "GOSTSHYP: ", ninactive, " of ", cavity%ngrid, &
            & " grid points below the density-overlap floor"
         call self%ctx%message(trim(report), level=2)
      end if

   end subroutine gostshyp_get_energy

   !> Hand the host the amplitudes conjugate to its Gaussian integral blocks
   !>
   !> The host rebuilds its Fock contribution as
   !>
   !> `F_uv += sum_i [w_overlap(i) g_uv,i + w_normal_deriv(i) f_uv,i]`; both signs
   !> are folded in here
   !>
   !> @param[in,out] self      Component instance
   !> @param[in]    coupling  Host coupling data
   !> @param[in,out] cavity    Live model cavity
   !> @param[in,out] response  Response accumulator
   !> @param[out]   error     Error handling
   subroutine gostshyp_get_response(self, coupling, cavity, response, error)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(inout) :: self
      !> Host coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Potential accumulator
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Host-supplied Gaussian moments
      real(wp), allocatable :: gt(:), pt(:, :)
      !> Gaussian widths, traces and amplitudes
      real(wp), allocatable :: omega(:), ftilde(:), alpha(:), beta(:)
      !> Amplitude item of this component
      type(gostshyp_amplitude_response_type) :: item

      ! A component switched off by zero pressure or zero scale still publishes
      ! its item, filled with zeros: present and contributing nothing, which is
      ! a different statement from having no GOSTSHYP component at all, and
      ! only the latter may leave the item absent
      allocate (item%w_overlap(cavity%ngrid), source=0.0_wp)
      allocate (item%w_normal_deriv(cavity%ngrid), source=0.0_wp)

      if (self%pressure /= 0.0_wp .and. self%scale /= 0.0_wp) then
         call coupling%check_mandatory(coupling%phase, error)
         if (allocated(error)) return
         call read_moments(coupling, cavity, gt, pt, error=error)
         if (allocated(error)) return
         call gostshyp_amplitudes(self, cavity, gt, pt, omega, ftilde, alpha, beta)

         item%w_overlap = item%w_overlap + alpha
         item%w_normal_deriv = item%w_normal_deriv - beta
      end if

      call response_accumulate(response, item, error)

   end subroutine gostshyp_get_response

   !> Add the GOSTSHYP surface adjoints
   !>
   !> With `G = exp(-w |r - C|^2)` and the outward normal held fixed, the
   !> supplied moments give every parameter derivative of the two traces,
   !>    dgtilde/dC_a = 2 w Pt_a          dftilde/dC_b = 2 w n_b gt - 4 w^2 (n.Mt)_b
   !>    dgtilde/dw   = -tr(Mt)           dftilde/dw   = -2 (n.Pt) + 2 w (n.Rt)
   !>    dftilde/dn   = -2 w Pt
   !>
   !> The area enters twice: explicitly through the amplitude `p_i`, and through
   !> the Gaussian width `w_i = pi ln2 / a_i`, whence the `-w_i/a_i` chain
   !> factor on the width route; the switching factor carries no dependence at
   !> all -- `w_f` is exactly zero -- because the Gaussian width is the only
   !> route by which the area reaches the level set
   !>
   !> These same weights serve the nuclear gradient: the base class points
   !>
   !> `get_gradient_surface_weights` here, and the cavity contracts them once in
   !> reverse mode
   !>
   !> @param[in,out] self     Component instance
   !> @param[in]    coupling Host coupling data
   !> @param[in]    cavity   Live model cavity
   !> @param[in,out] acc      Surface-adjoint accumulator
   !> @param[out]   error    Error handling
   subroutine gostshyp_get_surface_weights(self, coupling, cavity, acc, error)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(inout) :: self
      !> Host coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(in) :: cavity
      !> Surface accumulator
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Host-supplied Gaussian moments
      real(wp), allocatable :: gt(:), pt(:, :), mt(:, :, :), rt(:, :)
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
      if (coupling%phase == moist_phase_response) then
         if (.not. cavity%has_field_dependent_geometry()) return
      end if
      call read_moments(coupling, cavity, gt, pt, mt, rt, error)
      if (allocated(error)) return
      call gostshyp_amplitudes(self, cavity, gt, pt, omega, ftilde, alpha, beta)

      allocate (w_a(cavity%ngrid), source=0.0_wp)
      allocate (w_xyz(3, cavity%ngrid), source=0.0_wp)
      allocate (w_n(3, cavity%ngrid), source=0.0_wp)

      do igrid = 1, cavity%ngrid
         !> Inactive grid points were zeroed in both amplitudes, and the width
         !>
         !> route below would divide by a trace that carries only round-off
         if (alpha(igrid) == 0.0_wp) cycle

         n_pt = dot_product(cavity%normal0(:, igrid), pt(:, igrid))
         n_mt = matmul(cavity%normal0(:, igrid), mt(:, :, igrid))
         n_rt = dot_product(cavity%normal0(:, igrid), rt(:, igrid))

         dgdr = 2.0_wp*omega(igrid)*pt(:, igrid)
         dfdr = 2.0_wp*omega(igrid)*cavity%normal0(:, igrid)*gt(igrid) &
            & - 4.0_wp*omega(igrid)**2*n_mt
         dgdw = -(mt(1, 1, igrid) + mt(2, 2, igrid) + mt(3, 3, igrid))
         dfdw = -2.0_wp*n_pt + 2.0_wp*omega(igrid)*n_rt

         w_xyz(:, igrid) = alpha(igrid)*dgdr - beta(igrid)*dfdr
         !> Only ftilde depends on the normal, through ftilde = n . (-2 w Pt)
         w_n(:, igrid) = 2.0_wp*omega(igrid)*beta(igrid)*pt(:, igrid)
         w_a(igrid) = self%scale*self%pressure*gt(igrid)/ftilde(igrid) &
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
   !> of which `get_surface_weights` already states; the reverse-mode path
   !> contracts those once against the cavity's own nuclear derivatives; a
   !> forward-mode implementation would be a second, independently maintained
   !> derivation of the same numbers, so it is refused rather than approximated
   !>
   !> A component switched off by zero pressure or zero scale still publishes
   !> its zero amplitudes, as in the response phase
   !>
   !> @param[in,out] self     Component instance
   !> @param[in]    coupling Host coupling data, unused
   !> @param[in,out] cavity   Live model cavity, unused
   !> @param[in,out] response Host part of the gradient phase (amplitudes)
   !> @param[in,out] gradient Nuclear-gradient accumulator, unchanged
   !> @param[out]   error    Error handling
   subroutine gostshyp_get_gradient(self, coupling, cavity, response, gradient, error)
      !> Component instance
      class(solvation_model_component_gostshyp), intent(inout) :: self
      !> Host coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Host part of the gradient phase
      type(response_type), intent(inout) :: response
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (self%pressure == 0.0_wp .or. self%scale == 0.0_wp) then
         call self%get_response(coupling, cavity, response, error)
         return
      end if
      call coupling%check_mandatory(coupling%phase, error)
      if (allocated(error)) return
      call fatal_error(error, "GOSTSHYP has no forward-mode nuclear gradient; use the "// &
         & "reverse-mode surface path (force_forward_gradient must stay disabled)")

   end subroutine gostshyp_get_gradient

end module moist_model_component_gostshyp
