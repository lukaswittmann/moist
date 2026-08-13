!> Pressure-volume energy contribution
module moist_model_component_pv
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_type, only: solvation_model_component_type, cavity_type
   use moist_channels, only: coupling_type, response_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type

   implicit none (type, external)
   private

   public :: solvation_model_component_pv, new_component_pv

   !> Pressure-volume energy contribution `pressure * cavity_volume`
   type, extends(solvation_model_component_type) :: solvation_model_component_pv
      !> Pressure multiplying the cavity volume in atomic units
      real(wp) :: pressure = 0.0_wp
   contains
      procedure :: update => pv_update
      procedure :: get_energy => pv_get_energy
      procedure :: get_response => pv_get_response
      procedure :: get_gradient => pv_get_gradient
      procedure :: get_surface_weights => pv_get_surface_weights
   end type solvation_model_component_pv

contains

   !> Construct a pressure-volume energy component
   !>
   !> @param[out] self     Component instance
   !> @param[in]  pressure Pressure multiplying the cavity volume
   subroutine new_component_pv(self, pressure)
      !> Component instance
      type(solvation_model_component_pv), intent(out) :: self
      !> Pressure multiplying the cavity volume
      real(wp), intent(in) :: pressure

      self%name = "PV"
      self%pressure = pressure

   end subroutine new_component_pv

   !> Bind the current molecular structure
   !>
   !> @param[inout] self   Component instance
   !> @param[in]    mol    Molecular structure
   !> @param[inout] cavity Live model cavity
   !> @param[out]   error  Error handling
   subroutine pv_update(self, mol, cavity, error)
      !> Component instance
      class(solvation_model_component_pv), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

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
      class(solvation_model_component_pv), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
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
      class(solvation_model_component_pv), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
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
   !> `V = sum_i a_i (r_i . n_i)/3`, so its surface adjoints follow directly from
   !> the per-point geometry. The area channel carries `dV/da_i = (r_i . n_i)/3`;
   !> the cavity folds it into whichever primitive channels it is built from.
   !>
   !> @param[inout] self     Component instance
   !> @param[in]    coupling Host coupling data, unused
   !> @param[in]    cavity   Live model cavity
   !> @param[inout] acc      Surface-adjoint accumulator
   !> @param[out]   error    Error handling
   subroutine pv_get_surface_weights(self, coupling, cavity, acc, error)
      !> Component instance
      class(solvation_model_component_pv), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(in) :: cavity
      !> Surface accumulator
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Pressure-scaled adjoints of the grid point areas, positions, and normals
      real(wp), allocatable :: w_a(:), w_xyz(:, :), w_n(:, :)
      !> Grid point index
      integer :: igrid

      if (self%pressure == 0.0_wp) return
      if (.not. allocated(cavity%a) .or. .not. allocated(cavity%xyz) .or. &
          .not. allocated(cavity%normal0)) then
         call fatal_error(error, "Cavity volume surface data are incomplete")
         return
      end if

      allocate (w_a(cavity%ngrid))
      allocate (w_xyz(3, cavity%ngrid))
      allocate (w_n(3, cavity%ngrid))
      do igrid = 1, cavity%ngrid
         w_a(igrid) = self%pressure/3.0_wp &
            & *dot_product(cavity%xyz(:, igrid), cavity%normal0(:, igrid))
         w_xyz(:, igrid) = (self%pressure*cavity%a(igrid)/3.0_wp)*cavity%normal0(:, igrid)
         w_n(:, igrid) = (self%pressure*cavity%a(igrid)/3.0_wp)*cavity%xyz(:, igrid)
      end do
      call acc%add_surface_weights(error, w_a=w_a, w_xyz=w_xyz, w_n=w_n)

   end subroutine pv_get_surface_weights

   !> Add the pressure-scaled cavity-volume nuclear gradient
   !>
   !> The total volume is the grid sum of the per-grid point volume elements, so
   !> its nuclear derivative is the contraction of `v1_rA` over the grid.
   !>
   !> @param[inout] self     Component instance
   !> @param[in]    coupling Host coupling data, unused
   !> @param[inout] cavity   Live model cavity
   !> @param[inout] gradient Nuclear-gradient accumulator
   !> @param[out]   error    Error handling
   subroutine pv_get_gradient(self, coupling, cavity, gradient, error)
      !> Component instance
      class(solvation_model_component_pv), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
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
      !> once a gradient run has produced them.
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

end module moist_model_component_pv
