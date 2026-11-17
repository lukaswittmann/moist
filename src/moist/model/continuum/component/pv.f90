!> Pressure-volume energy contribution
module moist_model_continuum_component_pv
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_output_citations, only: citation
   use moist_context, only: moist_context_type
   use moist_cavity_type, only: cavity_type
   use moist_model_continuum_component_type, only: model_continuum_component_type, autogpa
   use moist_channels_coupling, only: coupling_type, coupling_view_type
   use moist_channels_response, only: response_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_surface_tangent, only: cavity_surface_tangent_type
   use moist_utils_prettyprint, only: prettyprinter

   implicit none(type, external)
   private

   public :: model_continuum_component_pv, new_component_pv

   !> Pressure-volume energy contribution `pressure * cavity_volume`
   type, extends(model_continuum_component_type) :: model_continuum_component_pv
      !> Pressure multiplying the cavity volume in atomic units
      real(wp) :: pressure = 0.0_wp
   contains
      procedure :: update => pv_update
      procedure :: get_energy => pv_get_energy
      procedure :: get_response => pv_get_response
      procedure :: get_gradient => pv_get_gradient
      procedure :: get_surface_weights => pv_get_surface_weights
      procedure :: supports_hessian => pv_supports_hessian
      procedure :: get_hessian_surface_weights => pv_get_hessian_surface_weights
      !> Print the pressure
      procedure :: print_inputs => pv_print_inputs
   end type model_continuum_component_pv

contains


   !> Construct a pressure-volume energy component
   !>
   !> @param[out] self     Component instance
   !> @param[in]  pressure Pressure multiplying the cavity volume
   !> @param[in]  ctx      Borrowed run context; omitted, a model supplies its own
   subroutine new_component_pv(self, pressure, ctx)
      !> Component instance
      type(model_continuum_component_pv), intent(out) :: self
      !> Pressure multiplying the cavity volume
      real(wp), intent(in) :: pressure
      !> Borrowed run context; omitted, a model supplies its own
      type(moist_context_type), intent(in), target, optional :: ctx

      self%name = "PV"
      self%description = "Pressure-volume work, pressure times cavity volume"
      self%pressure = pressure
      if (present(ctx)) self%ctx => ctx

   end subroutine new_component_pv

   !> Print the pressure in Eh/bohr^3 and GPa
   !>
   !> @param[in]    self Component instance
   !> @param[inout] pp   Pretty printer inside the component section
   subroutine pv_print_inputs(self, pp)
      !> Component instance
      class(model_continuum_component_pv), intent(in) :: self
      !> Pretty printer inside the component section
      type(prettyprinter), intent(inout) :: pp

      call pp%kv2("Pressure", self%pressure, "Eh/bohr^3", self%pressure*autogpa, "GPa", use_exp1=.true.)

   end subroutine pv_print_inputs

   !> Bind the current molecular structure
   !>
   !> @param[inout] self   Component instance
   !> @param[in]    mol    Molecular structure
   !> @param[inout] cavity Live model cavity
   !> @param[out]   error  Error handling
   subroutine pv_update(self, mol, cavity, error)
      !> Component instance
      class(model_continuum_component_pv), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%require_context(error)
      if (allocated(error)) return

      call self%ctx%cite(citation%pv)

      self%mol_solu = mol
      if (.not. allocated(cavity%total_volume)) then
         call fatal_error(error, "PV component requires an updated cavity")
      end if

   end subroutine pv_update

   !> Add the linear cavity-volume energy
   !>
   !> @param[inout] self     Component instance
   !> @param[in]    coupling Host coupling data, unused
   !> @param[inout] cavity   Live model cavity
   !> @param[inout] energy   Energy accumulator
   !> @param[out]   error    Error handling
   subroutine pv_get_energy(self, coupling, cavity, energy, error)
      !> Component instance
      class(model_continuum_component_pv), intent(inout) :: self
      !> Host coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Energy accumulator
      real(wp), intent(inout) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(cavity%total_volume)) then
         call fatal_error(error, "Cavity volume is unavailable")
         return
      end if
      energy = energy + self%pressure*cavity%total_volume

   end subroutine pv_get_energy

   !> No direct host-trace response is produced by a volume contribution
   !>
   !> @param[inout] self      Component instance
   !> @param[in]    coupling  Host coupling data, unused
   !> @param[inout] cavity    Live model cavity, unused
   !> @param[inout] response  Response accumulator, unchanged
   !> @param[out]   error     Error handling
   subroutine pv_get_response(self, coupling, cavity, response, error)
      !> Component instance
      class(model_continuum_component_pv), intent(inout) :: self
      !> Host coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Potential accumulator
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine pv_get_response

   !> Add pressure-scaled total-volume surface adjoints
   !>
   !> The enclosed volume is the divergence-theorem sum over the grid points,
   !> `V = sum_i a_i ((r_i - c) . n_i)/3` with `c` the cavity's volume origin, so
   !> its surface adjoints follow directly from the per-point geometry; the area
   !> channel carries `dV/da_i = ((r_i - c) . n_i)/3`, the cavity folds it into
   !> whichever primitive channels it is built from, and the origin channel
   !> carries `dV/dc = -sum_i a_i n_i/3`
   !>
   !> @param[inout] self     Component instance
   !> @param[in]    coupling Host coupling data, unused
   !> @param[in]    cavity   Live model cavity
   !> @param[inout] acc      Surface-adjoint accumulator
   !> @param[out]   error    Error handling
   subroutine pv_get_surface_weights(self, coupling, cavity, acc, error)
      !> Component instance
      class(model_continuum_component_pv), intent(inout) :: self
      !> Host coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(in) :: cavity
      !> Surface accumulator
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Pressure-scaled adjoints of the grid point areas, positions, and normals
      real(wp), allocatable :: w_a(:), w_xyz(:, :), w_n(:, :)
      !> Volume origin and its pressure-scaled adjoint
      real(wp) :: origin(3), w_c(3)
      !> Grid point index
      integer :: igrid

      if (self%pressure == 0.0_wp) return
      if (.not. allocated(cavity%a) .or. .not. allocated(cavity%xyz) .or. &
          .not. allocated(cavity%normal0) .or. .not. allocated(cavity%sphxyz)) then
         call fatal_error(error, "Cavity volume surface data are incomplete")
         return
      end if

      origin = cavity%get_volume_origin()
      allocate (w_a(cavity%ngrid))
      allocate (w_xyz(3, cavity%ngrid))
      allocate (w_n(3, cavity%ngrid))
      w_c = 0.0_wp
      do igrid = 1, cavity%ngrid
         w_a(igrid) = self%pressure/3.0_wp &
            & *dot_product(cavity%xyz(:, igrid) - origin, cavity%normal0(:, igrid))
         w_xyz(:, igrid) = (self%pressure*cavity%a(igrid)/3.0_wp)*cavity%normal0(:, igrid)
         w_n(:, igrid) = (self%pressure*cavity%a(igrid)/3.0_wp)*(cavity%xyz(:, igrid) - origin)
         ! The volume reads positions only as r_i - c
         w_c = w_c - w_xyz(:, igrid)
      end do
      call acc%add_surface_weights(error, w_a=w_a, w_xyz=w_xyz, w_n=w_n, w_c=w_c)

   end subroutine pv_get_surface_weights

   !> Report second-order support for the pressure-volume component
   !>
   !> @param[in] self  component instance
   pure function pv_supports_hessian(self) result(supported)
      !> Component instance
      class(model_continuum_component_pv), intent(in) :: self
      !> Second-order support flag
      logical :: supported

      supported = .true.

   end function pv_supports_hessian

   !> Add directional tangents of the pressure-scaled total-volume surface adjoints
   !>
   !> @param[in,out] self      component instance
   !> @param[in]     coupling  scoped host data, unused
   !> @param[in]     cavity    live model cavity
   !> @param[in]     dirs      nuclear directions `(3, nsph, nblk)`
   !> @param[in]     tangent   surface tangent of the block
   !> @param[in,out] dacc      surface-adjoint tangent per direction `(nblk)`
   !> @param[out]    error     error handling
   !> @param[in]     first     global index of the block's first direction, optional
   subroutine pv_get_hessian_surface_weights(self, coupling, cavity, dirs, tangent, dacc, &
                                             error, first)
      !> Component instance
      class(model_continuum_component_pv), intent(inout) :: self
      !> Scoped host data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(in) :: cavity
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Surface tangent of the block
      type(cavity_surface_tangent_type), intent(in) :: tangent
      !> Surface-adjoint tangent per direction
      type(cavity_surface_adjoint_type), intent(inout) :: dacc(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Global index of the block's first direction
      integer, intent(in), optional :: first

      !> Tangents of grid point area, position, normal, and volume-origin adjoints
      real(wp), allocatable :: dw_a(:), dw_xyz(:, :), dw_n(:, :)
      real(wp) :: origin(3), dc(3), dr(3), r_rel(3), dw_c(3), p3
      integer :: ngrid, ndir, igrid, idir

      if (self%pressure == 0.0_wp) return
      if (.not. allocated(cavity%a) .or. .not. allocated(cavity%xyz) .or. &
          .not. allocated(cavity%normal0) .or. .not. allocated(cavity%sphxyz)) then
         call fatal_error(error, "Cavity volume surface data are incomplete")
         return
      end if
      ngrid = cavity%ngrid
      ndir = size(dirs, 3)
      if (size(dacc) /= ndir .or. .not. tangent%is_initialized() .or. tangent%ndir() /= ndir) then
         call fatal_error(error, "PV Hessian surface-tangent shape mismatch")
         return
      end if
      if (size(dirs, 1) /= 3 .or. size(dirs, 2) /= cavity%nsph) then
         call fatal_error(error, "PV Hessian direction shape mismatch")
         return
      end if

      origin = cavity%get_volume_origin()
      p3 = self%pressure/3.0_wp
      allocate (dw_a(ngrid), dw_xyz(3, ngrid), dw_n(3, ngrid))
      do idir = 1, ndir
         dc = sum(dirs(:, :cavity%nsph, idir), dim=2)/real(max(1, cavity%nsph), wp)
         dw_c = 0.0_wp
         do igrid = 1, ngrid
            r_rel = cavity%xyz(:, igrid) - origin
            dr = tangent%d_xyz(:, igrid, idir) - dc
            dw_a(igrid) = p3*(dot_product(dr, cavity%normal0(:, igrid)) + &
                              dot_product(r_rel, tangent%d_n(:, igrid, idir)))
            dw_xyz(:, igrid) = p3*(tangent%d_a(igrid, idir)*cavity%normal0(:, igrid) + &
                                   cavity%a(igrid)*tangent%d_n(:, igrid, idir))
            dw_n(:, igrid) = p3*(tangent%d_a(igrid, idir)*r_rel + cavity%a(igrid)*dr)
            dw_c = dw_c - dw_xyz(:, igrid)
         end do
         call dacc(idir)%add_surface_weights(error, w_a=dw_a, w_xyz=dw_xyz, w_n=dw_n, w_c=dw_c)
         if (allocated(error)) return
      end do

   end subroutine pv_get_hessian_surface_weights

   !> Add the pressure-scaled cavity-volume nuclear gradient
   !>
   !> The total volume is the grid sum of the per-grid point volume elements, so
   !> its nuclear derivative is the contraction of `v1_rA` over the grid
   !>
   !> @param[inout] self     Component instance
   !> @param[in]    coupling Host coupling data, unused
   !> @param[inout] cavity   Live model cavity
   !> @param[inout] response Host part of the gradient phase, unchanged
   !> @param[inout] gradient Nuclear-gradient accumulator
   !> @param[out]   error    Error handling
   subroutine pv_get_gradient(self, coupling, cavity, response, gradient, error)
      !> Component instance
      class(model_continuum_component_pv), intent(inout) :: self
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

      if (self%pressure == 0.0_wp) return
      if (any(shape(gradient) /= [3, cavity%nsph])) then
         call fatal_error(error, "Cavity-volume gradient shape mismatch")
         return
      end if

      !> The per-grid point volume derivatives belong to the current geometry only
      !> once a gradient run has produced them
      if (.not. allocated(cavity%v1_rA)) then
         call cavity%get_gradient(error)
         if (allocated(error)) return
      end if
      if (.not. allocated(cavity%v1_rA)) then
         call fatal_error(error, "Cavity does not provide volume nuclear derivatives")
         return
      end if
      if (any(shape(cavity%v1_rA) /= [3, cavity%nsph, cavity%ngrid])) then
         call fatal_error(error, "Cavity volume-derivative shape mismatch")
         return
      end if

      gradient = gradient + self%pressure*sum(cavity%v1_rA, dim=3)

   end subroutine pv_get_gradient

end module moist_model_continuum_component_pv
