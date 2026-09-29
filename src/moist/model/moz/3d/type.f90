!> Three-dimensional MOZ model state
module moist_model_moz_3d
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_model_type, only: solvation_model_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type
   implicit none(type, external)
   private
   public :: model_moz_3d_type, new_moz_3d_model

   !> Direct 3D MOZ model; correlation theory is pending
   type, extends(solvation_model_type) :: model_moz_3d_type
      private
      !> Owned spatial grid
      class(moist_math_grid_3d_type), allocatable :: grid
      !> Whether the grid reflects the latest successful update
      logical :: updated = .false.
   contains
      procedure :: update => moz_3d_update
      procedure :: get_energy => moz_3d_energy
      procedure :: get_response => moz_3d_response
      procedure :: get_gradient => moz_3d_gradient
      procedure :: get_grid => moz_3d_get_grid
      procedure :: atom_count => moz_3d_atom_count
      procedure :: is_updated => moz_3d_is_updated
   end type model_moz_3d_type
contains
   !> Copy the configured spatial grid
   !>
   !> @param[out] self Model
   !> @param[in] grid Spatial grid template
   !> @param[in] ctx Run context
   !> @param[out] error Allocation failure
   subroutine new_moz_3d_model(self, grid, ctx, error)
      !> Self
      type(model_moz_3d_type), intent(out) :: self
      !> Grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Ctx
      type(moist_context_type), intent(in), target :: ctx
      !> Error
      type(error_type), allocatable, intent(out) :: error
      !> Stat
      integer :: stat
      allocate (self%grid, source=grid, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to copy 3D MOZ grid")
         return
      end if
      self%ctx => ctx
      self%grid%ctx => ctx
   end subroutine new_moz_3d_model

   !> Update owned grid geometry
   !>
   !> @param[in,out] self Model
   !> @param[in] mol Solute structure
   !> @param[out] error Grid error
   subroutine moz_3d_update(self, mol, error)
      !> Self
      class(model_moz_3d_type), intent(inout) :: self
      !> Mol
      class(structure_type), intent(in) :: mol
      !> Error
      type(error_type), allocatable, intent(out) :: error
      self%updated = .false.
      if (.not. allocated(self%grid)) then
         call fatal_error(error, "Construct the 3D MOZ model before update")
         return
      end if
      call self%grid%update(mol, error)
      if (allocated(error)) return
      self%updated = .true.
   end subroutine moz_3d_update

   !> Borrow the owned grid; model must outlive the pointer
   !>
   !> @param[in] self Self
   function moz_3d_get_grid(self) result(grid)
      !> Self
      class(model_moz_3d_type), intent(in), target :: self
      !> Grid
      class(moist_math_grid_3d_type), pointer :: grid
      nullify (grid)
      if (allocated(self%grid)) grid => self%grid
   end function moz_3d_get_grid

   !> Number of solute atoms in the grid
   !>
   !> @param[in] self Self
   function moz_3d_atom_count(self) result(natom)
      !> Self
      class(model_moz_3d_type), intent(in) :: self
      !> Natom
      integer :: natom
      natom = 0
      if (allocated(self%grid)) natom = self%grid%natom
   end function moz_3d_atom_count

   !> Whether the grid update succeeded
   !>
   !> @param[in] self Self
   function moz_3d_is_updated(self) result(updated)
      !> Self
      class(model_moz_3d_type), intent(in) :: self
      !> Updated
      logical :: updated
      updated = self%updated
   end function moz_3d_is_updated

   !> Report pending 3D MOZ energy theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] energy Energy
   !> @param[out] error Error
   subroutine moz_3d_energy(self, coupling, energy, error)
      !> Self
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Energy
      real(wp), intent(inout) :: energy
      !> Error
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "3D MOZ energy is not implemented")
   end subroutine moz_3d_energy

   !> Report pending 3D MOZ response theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] response Response
   !> @param[out] error Error
   subroutine moz_3d_response(self, coupling, response, error)
      !> Self
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Response
      type(response_type), intent(inout) :: response
      !> Error
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "3D MOZ response is not implemented")
   end subroutine moz_3d_response

   !> Report pending 3D MOZ gradient theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] response Response
   !> @param[in,out] gradient Gradient
   !> @param[out] error Error
   subroutine moz_3d_gradient(self, coupling, response, gradient, error)
      !> Self
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Response
      type(response_type), intent(inout) :: response
      !> Gradient
      real(wp), intent(inout) :: gradient(:, :)
      !> Error
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "3D MOZ gradient is not implemented")
   end subroutine moz_3d_gradient
end module moist_model_moz_3d
