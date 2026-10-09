!> GOSTSHYP hydrostatic pressure with unit-integral surface Gaussians
!>
!> G_i = (omega_i/pi)**1.5 exp(-omega_i |r-C_i|**2), omega_i = pi ln2/a_i
!> E_i = scale pressure a_i gt_i R(ftilde_i), ftilde_i = -2 omega_i n_i.pt_i
!> R(f) = S(|f|)/f, S a C-infinity switch from 0 at regularization_start to 1 at regularization_end
!> Narrow points are switched off smoothly by area before host integrals lose precision
!> Host moments and their AO integral blocks must use the same normalization
module moist_model_continuum_component_gostshyp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_context, only: moist_context_type
   use moist_cavity_type, only: cavity_type
   use moist_model_continuum_component_type, only: model_continuum_component_type, autogpa
   use moist_channels_response, only: response_type, gaussian_amplitude_response_type, &
      & response_accumulate
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & gaussian_moment_request_type, moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, coupling_register, request_require
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_drop_switching, only: moist_cavity_drop_swif_sigmoid_bump_type, new_swif_sigmoid_bump
   use moist_utils_prettyprint, only: prettyprinter

   implicit none(type, external)
   private

   public :: model_continuum_component_gostshyp, new_component_gostshyp, moist_gostshyp_parameters_type

   !> Fixed numerical regularization of normalized Gaussian gradient traces
   type, extends(moist_model_parameters_type) :: moist_gostshyp_parameters_type
      !> Trace magnitude at and below which a point is off, bohr**-4
      real(wp) :: regularization_start = 1.0e-12_wp
      !> Trace magnitude at and above which the reciprocal is exact, bohr**-4
      real(wp) :: regularization_end = 1.0e-10_wp
      !> Remove negative pressure amplitudes; the switch makes the cut smooth at zero
      logical :: suppress_negative_amplitudes = .false.
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Register fields for JSON, TOML and printing
      procedure :: register_entries => register_parameter_entries
      !> Validate the fixed transition window
      procedure :: validate => validate_gostshyp_parameters
   end type moist_gostshyp_parameters_type

   !> Pi times ln 2, the numerator of the Gaussian width
   real(wp), parameter :: pi_ln2 = acos(-1.0_wp)*log(2.0_wp)

   !> Gaussian widths, bohr**-2, bracketing the C2 switch-off of narrow points
   !>
   !> Host integral codes lose the higher moments of narrow Gaussians; libcint's
   !> rt error grows about omega**2, 2e-9 relative at 1e3 and 3e-7 at 1e4
   real(wp), parameter :: switch_width_on = 1.0e3_wp, switch_width_off = 1.0e4_wp

   !> GOSTSHYP hydrostatic pressure contribution
   type, extends(model_continuum_component_type) :: model_continuum_component_gostshyp
      !> Applied hydrostatic pressure in atomic units, Hartree/bohr**3
      real(wp) :: pressure = 0.0_wp
      !> Numerical settings, constant throughout the calculation
      type(moist_gostshyp_parameters_type) :: param
   contains
      procedure :: update => gostshyp_update
      procedure :: get_energy => gostshyp_get_energy
      procedure :: get_response => gostshyp_get_response
      procedure :: get_gradient => gostshyp_get_gradient
      procedure :: get_surface_weights => gostshyp_get_surface_weights
      !> Declare the Gaussian-moment request
      procedure :: declare_coupling => gostshyp_declare_coupling
      !> Print the applied pressure
      procedure :: print_inputs => gostshyp_print_inputs
      !> Numerical settings of the component
      procedure :: parameters => gostshyp_component_parameters
   end type model_continuum_component_gostshyp

contains


   !* ================================================================================= *!
   !*                                    Constructor                                    *!
   !* ================================================================================= *!

   !> Construct a GOSTSHYP pressure component
   !>
   !> @param[out] self     Component instance
   !> @param[in]  pressure Applied pressure in Hartree/bohr**3
   !> @param[in]  ctx      Borrowed run context; omitted, a model supplies its own
   !> @param[in]  param    Fixed regularization settings
   !> @param[out] error    Optional construction diagnostic
   subroutine new_component_gostshyp(self, pressure, ctx, param, error)
      !> Component instance
      type(model_continuum_component_gostshyp), intent(out) :: self
      !> Applied pressure
      real(wp), intent(in) :: pressure
      !> Borrowed run context; omitted, a model supplies its own
      type(moist_context_type), intent(in), target, optional :: ctx
      !> Numerical settings; omitted selects compiled defaults
      type(moist_gostshyp_parameters_type), intent(in), optional :: param
      !> Optional construction diagnostic; evaluation also validates settings
      type(error_type), allocatable, intent(out), optional :: error

      self%name = "GOSTSHYP"
      self%description = "Gaussians On Surface Tesserae To Simulate HYdrostatic Pressure"
      self%pressure = pressure
      if (present(param)) self%param = param
      if (present(ctx)) self%ctx => ctx
      if (present(error)) call validate_settings(self, error)

   end subroutine new_component_gostshyp

   !> Print the applied pressure in Eh/bohr^3 and GPa
   !>
   !> @param[in]    self Component instance
   !> @param[in,out] pp   Pretty printer inside the component section
   subroutine gostshyp_print_inputs(self, pp)
      !> Component instance
      class(model_continuum_component_gostshyp), intent(in) :: self
      !> Pretty printer inside the component section
      type(prettyprinter), intent(inout) :: pp

      call pp%kv2("Pressure", self%pressure, "Eh/bohr^3", self%pressure*autogpa, "GPa", use_exp1=.true.)

   end subroutine gostshyp_print_inputs

   !> Bind the current molecular structure
   !>
   !> @param[in,out] self   Component instance
   !> @param[in]    mol    Molecular structure
   !> @param[in,out] cavity Live model cavity
   !> @param[out]   error  Error handling
   subroutine gostshyp_update(self, mol, cavity, error)
      !> Component instance
      class(model_continuum_component_gostshyp), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call validate_settings(self, error)
      if (allocated(error)) return
      call self%require_context(error)
      if (allocated(error)) return

      self%mol_solu = mol
      if (.not. allocated(cavity%a) .or. .not. allocated(cavity%xyz) .or. &
          .not. allocated(cavity%normal0)) then
         call fatal_error(error, "GOSTSHYP component requires an updated cavity")
      end if

   end subroutine gostshyp_update

   !* ================================================================================= *!
   !*                              Shared state evaluation                              *!
   !* ================================================================================= *!

   !> Restore compiled GOSTSHYP parameter defaults
   !>
   !> @param[in,out] self Numerical settings
   subroutine init_parameter_defaults(self)
      !> Numerical settings
      class(moist_gostshyp_parameters_type), intent(inout) :: self
      !> Compiled defaults
      type(moist_gostshyp_parameters_type) :: defaults
      self%regularization_start = defaults%regularization_start
      self%regularization_end = defaults%regularization_end
      self%suppress_negative_amplitudes = defaults%suppress_negative_amplitudes
   end subroutine init_parameter_defaults

   !> Declare GOSTSHYP fields for JSON, TOML and printing
   !>
   !> @param[in,out] self Numerical settings
   subroutine register_parameter_entries(self)
      !> Numerical settings
      class(moist_gostshyp_parameters_type), intent(inout), target :: self
      call self%register_real_scalar("regularization_start", self%regularization_start)
      call self%register_real_scalar("regularization_end", self%regularization_end)
      call self%register_logical("suppress_negative_amplitudes", self%suppress_negative_amplitudes)
   end subroutine register_parameter_entries

   !> Validate the fixed reciprocal transition window
   !>
   !> The switch degrades to a step for windows narrower than epsilon, so those are refused
   !>
   !> @param[in,out] self Numerical settings
   !> @param[out] error Invalid window diagnostic
   subroutine validate_gostshyp_parameters(self, error)
      !> Numerical settings
      class(moist_gostshyp_parameters_type), intent(inout) :: self
      !> Invalid settings diagnostic
      type(error_type), allocatable, intent(out) :: error
      if (.not. ieee_is_finite(self%regularization_start) .or. .not. ieee_is_finite(self%regularization_end)) then
         call fatal_error(error, "GOSTSHYP: regularization start and end must be finite")
      else if (self%regularization_start < 0.0_wp) then
         call fatal_error(error, "GOSTSHYP: regularization start must be nonnegative")
      else if (self%regularization_end - self%regularization_start <= epsilon(1.0_wp)) then
         call fatal_error(error, "GOSTSHYP: regularization end must exceed start by more than epsilon")
      end if
   end subroutine validate_gostshyp_parameters

   !> Registered numerical settings of the component
   !>
   !> @param[in] self Pressure component
   function gostshyp_component_parameters(self) result(param)
      !> Pressure component
      class(model_continuum_component_gostshyp), intent(in), target :: self
      !> Numerical settings, valid while the component is
      class(moist_model_parameters_type), pointer :: param
      param => self%param
   end function gostshyp_component_parameters

   !> Validate settings even when the pressure component is disabled
   !>
   !> @param[in] self Component settings
   !> @param[out] error Invalid settings diagnostic
   subroutine validate_settings(self, error)
      !> Component settings
      class(model_continuum_component_gostshyp), intent(in) :: self
      !> Invalid settings diagnostic
      type(error_type), allocatable, intent(out) :: error

      !> Copy permits the parameter validation interface on an immutable component
      type(moist_gostshyp_parameters_type) :: settings

      settings = self%param
      call settings%validate(error)
      if (allocated(error)) return
      if (.not. ieee_is_finite(self%pressure) .or. .not. ieee_is_finite(self%scale)) then
         call fatal_error(error, "GOSTSHYP: pressure and scale must be finite")
         return
      end if
      !> Divide rather than multiply, so the check itself cannot overflow
      if (abs(self%scale) > 1.0_wp) then
         if (abs(self%pressure) > huge(1.0_wp)/abs(self%scale)) then
            call fatal_error(error, "GOSTSHYP: pressure times scale is not representable")
         end if
      end if
   end subroutine validate_settings

   !> Build finite Gaussian exponents; points below the switch-off area are inert
   !>
   !> @param[in] cavity Live cavity
   !> @param[out] omega Finite Gaussian exponents, bohr**-2
   !> @param[out] error Geometry or arithmetic diagnostic
   subroutine gostshyp_widths(cavity, omega, error)
      !> Live cavity areas
      class(cavity_type), intent(in) :: cavity
      !> Gaussian exponents
      real(wp), allocatable, intent(out) :: omega(:)
      !> Geometry or arithmetic diagnostic
      type(error_type), allocatable, intent(out) :: error
      !> Point index
      integer :: i

      if (.not. allocated(cavity%a)) then
         call fatal_error(error, "GOSTSHYP requires cavity areas")
         return
      end if
      if (size(cavity%a) /= cavity%ngrid .or. .not. all(ieee_is_finite(cavity%a))) then
         call fatal_error(error, "GOSTSHYP: invalid cavity areas")
         return
      end if
      if (any(cavity%a < 0.0_wp)) then
         call fatal_error(error, "GOSTSHYP: cavity areas must be nonnegative")
         return
      end if
      allocate (omega(cavity%ngrid), source=0.0_wp)
      do i = 1, cavity%ngrid
         if (cavity%a(i) > pi_ln2/switch_width_off) omega(i) = pi_ln2/cavity%a(i)
      end do
      if (.not. all(ieee_is_finite(omega))) then
         call fatal_error(error, "GOSTSHYP: Gaussian width is not representable")
      end if
   end subroutine gostshyp_widths

   !> C2 switch of a point by its area, between the switch-off and switch-on widths
   !>
   !> @param[in] area Point area, bohr**2
   !> @param[out] switch One above pi ln2/switch_width_on, zero below pi ln2/switch_width_off
   !> @param[out] dswitch Area derivative of the switch
   pure subroutine area_switch(area, switch, dswitch)
      !> Point area
      real(wp), intent(in) :: area
      !> Switch value and its area derivative
      real(wp), intent(out) :: switch, dswitch
      !> Areas bracketing the transition
      real(wp), parameter :: area_on = pi_ln2/switch_width_on, area_off = pi_ln2/switch_width_off
      !> Reduced area inside the transition
      real(wp) :: t

      if (area >= area_on) then
         switch = 1.0_wp
         dswitch = 0.0_wp
      else if (area <= area_off) then
         switch = 0.0_wp
         dswitch = 0.0_wp
      else
         t = (area - area_off)/(area_on - area_off)
         switch = t*t*t*(10.0_wp + t*(-15.0_wp + 6.0_wp*t))
         dswitch = 30.0_wp*t*t*(1.0_wp - t)*(1.0_wp - t)/(area_on - area_off)
      end if
   end subroutine area_switch

   !> Evaluate the two derivatives of E = s(a) H a gt R(ftilde)
   !>
   !> alpha = s H a R, beta = -s H a gt R', R(f) = S(|f|)/f, R'(f) = S'(|f|)/|f| - S(|f|)/f**2
   !> S vanishes with all derivatives at the start, so R is smooth through zero
   !>
   !> @param[in] self Pressure component
   !> @param[in] cavity Live surface
   !> @param[in] gt Normalized scalar moments (ngrid)
   !> @param[in] pt Normalized first moments (3, ngrid)
   !> @param[out] omega Gaussian exponents, bohr**-2
   !> @param[out] ftilde Normal gradient traces, bohr**-4
   !> @param[out] alpha Overlap conjugates
   !> @param[out] beta Negative normal-gradient conjugates
   !> @param[out] error Invalid or unrepresentable result
   !> @param[out] counts Optional regularized, negative, suppressed and switched-off counts
   !> @param[out] w_switch Optional area adjoint of the switch, s'(a) H a gt R
   subroutine gostshyp_amplitudes(self, cavity, gt, pt, omega, ftilde, alpha, beta, error, counts, w_switch)
      !> Pressure component
      class(model_continuum_component_gostshyp), intent(in) :: self
      !> Live surface
      class(cavity_type), intent(in) :: cavity
      !> Normalized zeroth and first moments
      real(wp), intent(in) :: gt(:)
      !> Normalized first moments
      real(wp), intent(in) :: pt(:, :)
      !> Exponents and normal gradient traces
      real(wp), allocatable, intent(out) :: omega(:)
      !> Normal gradient traces
      real(wp), allocatable, intent(out) :: ftilde(:)
      !> Conjugates of overlap and negative normal derivative blocks
      real(wp), allocatable, intent(out) :: alpha(:)
      !> Negative normal-gradient conjugates
      real(wp), allocatable, intent(out) :: beta(:)
      !> Invalid or unrepresentable inputs/results
      type(error_type), allocatable, intent(out) :: error
      !> Regularized, negative candidate, suppressed and switched-off point counts
      integer, intent(out), optional :: counts(4)
      !> Area adjoint of the narrow-point switch
      real(wp), allocatable, intent(out), optional :: w_switch(:)
      !> Trace window, trace magnitude and scaled pressure
      real(wp) :: band_start, band_end, magnitude, pressure
      !> Trace switch S(|f|), its derivative, R(f) and R'(f) at one point
      real(wp) :: trace_switch, dtrace_switch, reciprocal, dreciprocal
      !> C-infinity trace switch rising from band_start to band_end
      type(moist_cavity_drop_swif_sigmoid_bump_type) :: bump
      !> Area switch value and area derivative at one point
      real(wp) :: switch, dswitch
      !> Point index and evaluation diagnostics
      integer :: i, stats(4)
      !> Sign before optional suppression, independent of arithmetic underflow
      logical :: negative

      call validate_settings(self, error)
      if (allocated(error)) return
      call gostshyp_widths(cavity, omega, error)
      if (allocated(error)) return
      if (.not. allocated(cavity%normal0)) then
         call fatal_error(error, "GOSTSHYP requires cavity normals")
         return
      end if
      if (any(shape(cavity%normal0) /= [3, cavity%ngrid]) &
         & .or. .not. all(ieee_is_finite(cavity%normal0))) then
         call fatal_error(error, "GOSTSHYP: invalid cavity normals")
         return
      end if
      allocate (ftilde(cavity%ngrid), alpha(cavity%ngrid), beta(cavity%ngrid))
      ftilde = 0.0_wp
      alpha = 0.0_wp
      beta = 0.0_wp
      if (present(w_switch)) allocate (w_switch(cavity%ngrid), source=0.0_wp)
      stats = 0
      band_start = self%param%regularization_start
      band_end = self%param%regularization_end
      call new_swif_sigmoid_bump(bump, band_start, band_end)
      pressure = self%scale*self%pressure
      do i = 1, cavity%ngrid
         if (cavity%a(i) > 0.0_wp .and. omega(i) == 0.0_wp) stats(4) = stats(4) + 1
         if (omega(i) == 0.0_wp) cycle
         ftilde(i) = -2.0_wp*omega(i)*dot_product(cavity%normal0(:, i), pt(:, i))
         if (.not. ieee_is_finite(ftilde(i))) cycle
         magnitude = abs(ftilde(i))
         if (magnitude < band_end) stats(1) = stats(1) + 1
         !> Off below the start, which also covers a vanishing trace
         if (magnitude <= band_start) cycle
         negative = ((pressure < 0.0_wp .and. ftilde(i) > 0.0_wp) .or. &
            & (pressure > 0.0_wp .and. ftilde(i) < 0.0_wp))
         if (negative) stats(2) = stats(2) + 1
         if (negative .and. self%param%suppress_negative_amplitudes) then
            stats(3) = stats(3) + 1
            cycle
         end if
         if (magnitude < band_end) then
            call bump%eval(magnitude, trace_switch, dtrace_switch)
            reciprocal = trace_switch/ftilde(i)
            dreciprocal = dtrace_switch/magnitude - reciprocal/ftilde(i)
            alpha(i) = pressure*cavity%a(i)*reciprocal
            beta(i) = -pressure*cavity%a(i)*gt(i)*dreciprocal
         else
            alpha(i) = pressure*cavity%a(i)/ftilde(i)
            beta(i) = (gt(i)/ftilde(i))*alpha(i)
         end if
         call area_switch(cavity%a(i), switch, dswitch)
         if (present(w_switch)) w_switch(i) = dswitch*alpha(i)*gt(i)
         alpha(i) = switch*alpha(i)
         beta(i) = switch*beta(i)
      end do
      if (.not. all(ieee_is_finite(ftilde)) &
         & .or. .not. all(ieee_is_finite(alpha)) .or. .not. all(ieee_is_finite(beta))) then
         call fatal_error(error, "GOSTSHYP: gradient trace or amplitude is not representable")
         return
      end if
      if (present(counts)) counts = stats
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
      if (.not. all(ieee_is_finite(gt)) .or. .not. all(ieee_is_finite(pt))) then
         call fatal_error(error, "GOSTSHYP: nonfinite host gt or pt moments")
         return
      end if
      if (present(mt)) then
         if (.not. all(ieee_is_finite(mt))) then
            call fatal_error(error, "GOSTSHYP: nonfinite host mt moments")
            return
         end if
      end if
      if (present(rt)) then
         if (.not. all(ieee_is_finite(rt))) call fatal_error(error, "GOSTSHYP: nonfinite host rt moments")
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
      class(model_continuum_component_gostshyp), intent(in) :: self
      !> Live cavity supplying moment centers and areas
      class(cavity_type), intent(in) :: cavity
      !> Component registration context
      type(coupling_type), intent(inout) :: coupling
      !> Registration error
      type(error_type), allocatable, intent(out) :: error
      type(gaussian_moment_request_type) :: moments
      call validate_settings(self, error)
      if (allocated(error)) return
      if (self%pressure == 0.0_wp .or. self%scale == 0.0_wp) return
      call gostshyp_widths(cavity, moments%width, error)
      if (allocated(error)) return
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
      class(model_continuum_component_gostshyp), intent(inout) :: self
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
      !> Regularized, negative candidate, suppressed and switched-off point counts
      integer :: counts(4)
      !> Diagnostic line
      character(len=160) :: report

      call validate_settings(self, error)
      if (allocated(error)) return
      if (self%pressure == 0.0_wp .or. self%scale == 0.0_wp) return
      call read_moments(coupling, cavity, gt, pt, error=error)
      if (allocated(error)) return
      call gostshyp_amplitudes(self, cavity, gt, pt, omega, ftilde, alpha, beta, error, counts)
      if (allocated(error)) return
      energy = energy + dot_product(alpha, gt)

      if (associated(self%ctx)) then
         if (counts(2) > 0) then
            write (report, "(a,i0,a,i0,a)") "GOSTSHYP warning: ", counts(2), &
               & " negative pressure amplitudes; ", counts(3), " suppressed"
            call self%ctx%message_flush(trim(report), level=1)
         end if
         write (report, "(a,i0,a,i0,a,i0,a,i0,a,i0)") "GOSTSHYP: ", counts(1), &
            & " regularized of ", cavity%ngrid, "; negative candidates: ", counts(2), &
            & "; suppressed: ", counts(3), "; switched off: ", counts(4)
         call self%ctx%message_flush(trim(report), level=2)
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
      class(model_continuum_component_gostshyp), intent(inout) :: self
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
      type(gaussian_amplitude_response_type) :: item

      call validate_settings(self, error)
      if (allocated(error)) return

      ! A component switched off by zero pressure or zero scale still publishes
      ! its item, filled with zeros: present and contributing nothing, which is
      ! a different statement from having no GOSTSHYP component at all, and
      ! only the latter may leave the item absent
      allocate (item%w_overlap(cavity%ngrid), source=0.0_wp)
      allocate (item%w_normal_deriv(cavity%ngrid), source=0.0_wp)

      if (self%pressure /= 0.0_wp .and. self%scale /= 0.0_wp) then
         call read_moments(coupling, cavity, gt, pt, error=error)
         if (allocated(error)) return
         call gostshyp_amplitudes(self, cavity, gt, pt, omega, ftilde, alpha, beta, error)
         if (allocated(error)) return

         item%w_overlap = item%w_overlap + alpha
         item%w_normal_deriv = item%w_normal_deriv - beta
      end if

      call response_accumulate(response, item, error)

   end subroutine gostshyp_get_response

   !> Add the GOSTSHYP surface adjoints
   !>
   !> With `G = (w/pi)**1.5 exp(-w |r - C|^2)` and the outward normal held fixed, the
   !> supplied moments give every parameter derivative of the two traces,
   !>    dgtilde/dC_a = 2 w Pt_a          dftilde/dC_b = 2 w n_b gt - 4 w^2 (n.Mt)_b
   !>    dgtilde/dw   = -tr(Mt) + 3 gt/(2w)
   !>    dftilde/dw   = -2 (n.Pt) + 2 w (n.Rt) + 3 ftilde/(2w)
   !>    dftilde/dn   = -2 w Pt
   !>
   !> The area enters three times: explicitly through the amplitude `p_i`, through
   !> the narrow-point switch `s(a_i)`, and through the Gaussian width and its
   !> normalization, whence the `-w_i/a_i` chain factor on the width route; the
   !> cavity switching factor carries no dependence at
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
      class(model_continuum_component_gostshyp), intent(inout) :: self
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
      !> Area adjoint of the narrow-point switch
      real(wp), allocatable :: w_switch(:)
      !> Adjoints of the grid-point areas, positions and normals
      real(wp), allocatable :: w_a(:), w_xyz(:, :), w_n(:, :)
      !> Parameter derivatives of the two traces at one grid point
      real(wp) :: dgdr(3), dfdr(3), dgdw, dfdw
      !> Normal projections of the supplied moments at one grid point
      real(wp) :: n_pt, n_mt(3), n_rt
      !> Grid-point index
      integer :: igrid

      call validate_settings(self, error)
      if (allocated(error)) return
      if (self%pressure == 0.0_wp .or. self%scale == 0.0_wp) return
      if (coupling%phase == moist_phase_response) then
         if (.not. cavity%has_field_dependent_geometry()) return
      end if
      call read_moments(coupling, cavity, gt, pt, mt, rt, error)
      if (allocated(error)) return
      call gostshyp_amplitudes(self, cavity, gt, pt, omega, ftilde, alpha, beta, error, w_switch=w_switch)
      if (allocated(error)) return

      allocate (w_a(cavity%ngrid), source=0.0_wp)
      allocate (w_xyz(3, cavity%ngrid), source=0.0_wp)
      allocate (w_n(3, cavity%ngrid), source=0.0_wp)

      do igrid = 1, cavity%ngrid
         if (alpha(igrid) == 0.0_wp .and. beta(igrid) == 0.0_wp) cycle

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
         w_a(igrid) = alpha(igrid)*gt(igrid)/cavity%a(igrid) &
            & + (alpha(igrid)*dgdw - beta(igrid)*dfdw)*(-omega(igrid)/cavity%a(igrid))
         !> dN/da contribution; the reciprocal branch cancels it exactly
         if (abs(ftilde(igrid)) < self%param%regularization_end) then
            w_a(igrid) = w_a(igrid) - 1.5_wp &
               & *(alpha(igrid)*gt(igrid) - beta(igrid)*ftilde(igrid))/cavity%a(igrid)
         end if
         w_a(igrid) = w_a(igrid) + w_switch(igrid)
      end do

      if (.not. all(ieee_is_finite(w_a)) .or. .not. all(ieee_is_finite(w_xyz)) &
         & .or. .not. all(ieee_is_finite(w_n))) then
         call fatal_error(error, "GOSTSHYP: surface weights are not representable")
         return
      end if
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
      class(model_continuum_component_gostshyp), intent(inout) :: self
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
      call fatal_error(error, "GOSTSHYP has no forward-mode nuclear gradient; use the "// &
         & "reverse-mode surface path (force_forward_gradient must stay disabled)")

   end subroutine gostshyp_get_gradient

end module moist_model_continuum_component_gostshyp
