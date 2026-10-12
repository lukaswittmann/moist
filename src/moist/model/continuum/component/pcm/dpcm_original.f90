!> DPCM original component backed by generated Gaussian and point kernels
module moist_model_continuum_component_pcm_dpcm_original
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_continuum_component_pcm_type, only: moist_pcm_parameters_type
   use moist_model_continuum_component_pcm_layer, only: pcm_layer_base, pcm_layer_feps
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver
   use moist_model_continuum_component_pcm_kernel_gauss_dpcm_original, only: dpcm_original_complex_state, &
      & dpcm_original_assemble, dpcm_original_medium, dpcm_original_solve, dpcm_original_weights
   use moist_model_continuum_component_pcm_kernel_point_dpcm_original, only: dpcm_original_point_complex_state, &
      & dpcm_original_point_assemble, dpcm_original_point_medium, dpcm_original_point_solve, &
      & dpcm_original_point_weights
   implicit none(type, external)
   private

   public :: model_continuum_component_dpcm_original, new_component_dpcm_original

   !> Concrete DPCM original component
   type, extends(pcm_layer_base) :: model_continuum_component_dpcm_original
      type(dpcm_original_complex_state) :: gauss
      type(dpcm_original_point_complex_state) :: point
   contains
      procedure :: assemble_layer => dpcm_original_assemble_layer
      procedure :: solve_layer => dpcm_original_solve_layer
      procedure :: layer_weights => dpcm_original_layer_weights
   end type model_continuum_component_dpcm_original

contains

   !> Construct DPCM original from physical and numerical parameters
   subroutine new_component_dpcm_original(self, epsilon, error, param, ctx)
      type(model_continuum_component_dpcm_original), intent(out) :: self
      !> Borrowed run context; omitted, a model supplies its own
      type(moist_context_type), intent(in), target, optional :: ctx
      real(wp), intent(in) :: epsilon
      type(error_type), allocatable, intent(out) :: error
      type(moist_pcm_parameters_type), intent(in), optional :: param
      type(moist_pcm_parameters_type) :: settings

      if (epsilon < 1.0_wp .or. epsilon /= epsilon) then
         call fatal_error(error, "DPCM original epsilon must be >= 1")
         return
      end if
      if (present(param)) settings = param
      call settings%validate(error)
      if (allocated(error)) return
      if (present(ctx)) self%ctx => ctx
      self%epsilon = epsilon
      self%feps = pcm_layer_feps(epsilon)
      self%param = settings
      self%name = "DPCM original"
      self%needs_switching = .false.
      self%needs_normal = .true.
   end subroutine new_component_dpcm_original

   !> Assemble the DPCM original operators
   subroutine dpcm_original_assemble_layer(self, cavity, solver, error)
      class(model_continuum_component_dpcm_original), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(pcm_solver), intent(in) :: solver
      type(error_type), allocatable, intent(out) :: error
      integer :: ngrid

      ngrid = cavity%ngrid
      associate (r => cavity%xyz(:, 1:ngrid), n => cavity%normal0(:, 1:ngrid), &
         & k1 => cavity%k1(1:ngrid), k2 => cavity%k2(1:ngrid), &
         & a => cavity%a(1:ngrid))
         if (self%smeared) then
            associate (xi => cavity%xi0(1:ngrid))
               call dpcm_original_assemble(self%gauss%dpcm_original_state, r, n, k1, k2, xi, a, error, solver)
               if (allocated(error)) return
               call dpcm_original_medium(self%gauss%dpcm_original_state, self%epsilon, error)
            end associate
         else
            call dpcm_original_point_assemble(self%point%dpcm_original_point_state, r, n, k1, k2, a, error, solver)
            if (allocated(error)) return
            call dpcm_original_point_medium(self%point%dpcm_original_point_state, self%epsilon, error)
         end if
      end associate
   end subroutine dpcm_original_assemble_layer

   !> Solve the current DPCM original host traces
   subroutine dpcm_original_solve_layer(self, error)
      class(model_continuum_component_dpcm_original), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      if (self%smeared) then
         call dpcm_original_solve(self%gauss%dpcm_original_state, self%phi, self%vn, &
            & self%kernel_energy, self%q, self%qsym, self%qn, error)
      else
         call dpcm_original_point_solve(self%point%dpcm_original_point_state, self%phi, self%vn, &
            & self%kernel_energy, self%q, self%qsym, self%qn, error)
      end if
   end subroutine dpcm_original_solve_layer

   !> Accumulate DPCM original operator surface weights
   subroutine dpcm_original_layer_weights(self, cavity, acc, error)
      class(model_continuum_component_dpcm_original), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_xyz(:, :), w_n(:, :), w_k1(:), w_k2(:)
      real(wp), allocatable :: w_a(:), w_xi(:), w_f(:)
      real(wp) :: w_eps
      integer :: ngrid

      ngrid = cavity%ngrid
      allocate (w_xyz(3, ngrid), w_n(3, ngrid), w_k1(ngrid), w_k2(ngrid), w_a(ngrid))
      if (self%smeared) then
         allocate (w_xi(ngrid), source=0.0_wp)
         call dpcm_original_weights(self%gauss%dpcm_original_state, w_xyz, w_n, w_k1, w_k2, w_xi, w_a, w_eps, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, &
            & w_k2=w_k2, w_a=w_a, w_xi=w_xi)
      else
         call dpcm_original_point_weights(self%point%dpcm_original_point_state, w_xyz, w_n, w_k1, w_k2, w_a, w_eps, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, &
            & w_k2=w_k2, w_a=w_a)
      end if
   end subroutine dpcm_original_layer_weights

end module moist_model_continuum_component_pcm_dpcm_original
