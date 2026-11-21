!> Grid model: decides the angular grid components of every atomic sphere
!>
!> A Lebedev-based cavity (iSwiG, NumSA, DROP) builds the angular grid of a
!> sphere as a weighted sum of components ([[moist_cavity_grid_type]]); a grid model
!> is the provider of those components, passed optionally into the cavity
!> constructor like a radius model
!>
!> - component: a Lebedev size, a rotation R and a smooth weight w; its unit
!>   node u sits at c_a + r_a R u, and the switching value of its points is
!>   scaled by w
!> - the model is updated once per geometry (structure and sphere radii) and
!>   then hands over the components of every sphere ([[grid_components_type]]);
!>   no model procedure is called inside a point loop
!> - gradients: a cavity returns per component the torque
!>   tau = sum_p (p - c_a) x dE/dp and the weight adjoint wbar = dE/dw; the
!>   model maps both to the atoms (`add_gradient`, reverse mode) or gives
!>   d omega/dc and dw/dc (`jacobian`, forward mode)
!> - a plain model (`is_plain`) leaves every sphere with its single unrotated
!>   default component of weight one: the cavity then uses its own nodes,
!>   bitwise as without a model
!> - concrete models: [[grid_model_symmetry_type]] (orientations by molecular frame and
!>   point group); further ones get a folder below `src/moist/cavity/grid/model`
module moist_cavity_grid_model_type
   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_channels_fields, only: field_query_type
   implicit none(type, external)
   private

   public :: grid_model_type, grid_components_type, grid_size_default

   !> Component size meaning the Lebedev size of the cavity itself
   integer, parameter :: grid_size_default = 0

   !> Components of every sphere, as a grid model hands them over
   type :: grid_components_type
      !> First component of each sphere (nsph + 1)
      !>
      !> - sphere a owns `sph_first(a):sph_first(a + 1) - 1`
      integer, allocatable :: sph_first(:)
      !> Lebedev size of each component, `grid_size_default` for the cavity's (ncomp)
      integer, allocatable :: size(:)
      !> Distinct rotation of each component, into `rotation` (ncomp)
      integer, allocatable :: rot(:)
      !> Smooth weight of each component; sums to one per sphere (ncomp)
      real(wp), allocatable :: weight(:)
      !> Distinct rotations R, unit node u at c + r R u (3, 3, nrot)
      real(wp), allocatable :: rotation(:, :, :)
   end type grid_components_type

   !> Provider of the angular grid components of every sphere
   type, abstract :: grid_model_type
   contains
      !> Decide the components for a structure unless already current
      procedure(update_iface), deferred :: update
      !> Whether a structure and radii match the latest update
      procedure(current_iface), deferred :: is_current
      !> Whether the settings ask for nothing, so no work is timed
      procedure(query_iface), deferred :: is_trivial
      !> Whether every sphere keeps its single plain default component
      procedure(query_iface), deferred :: is_plain
      !> Components of every sphere after an update
      procedure(components_iface), deferred :: components
      !> Map component torques and weight adjoints to the atoms (reverse mode)
      procedure(gradient_iface), deferred :: add_gradient
      !> Rotation and weight Jacobians per atom displacement (forward mode)
      procedure(jacobian_iface), deferred :: jacobian
      !> Declare the readable results of the model
      procedure(fields_iface), deferred :: list_fields
      !> Print the report of the latest update
      procedure(print_iface), deferred :: print
   end type grid_model_type

   !> Deferred procedures of a grid model
   abstract interface

      !> Decide the components of every sphere unless already current
      !>
      !> @param[in,out] self  Grid model
      !> @param[in]     mol   Structure
      !> @param[in]     radii Sphere radius per atom (nat), bohr
      !> @param[out]    error Invalid setting or failed update
      !> @param[in,out] ctx   Run context: timer and warnings
      subroutine update_iface(self, mol, radii, error, ctx)
         import :: grid_model_type, structure_type, wp, error_type, moist_context_type
         implicit none(type, external)
         !> Grid model
         class(grid_model_type), intent(inout) :: self
         !> Structure
         type(structure_type), intent(in) :: mol
         !> Sphere radii
         real(wp), intent(in) :: radii(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Run context
         type(moist_context_type), intent(inout), optional :: ctx
      end subroutine update_iface

      !> Whether a structure and radii match the latest update
      !>
      !> @param[in] self  Grid model
      !> @param[in] mol   Structure
      !> @param[in] radii Sphere radius per atom (nat), bohr
      function current_iface(self, mol, radii) result(current)
         import :: grid_model_type, structure_type, wp
         implicit none(type, external)
         !> Grid model
         class(grid_model_type), intent(in) :: self
         !> Structure
         type(structure_type), intent(in) :: mol
         !> Sphere radii
         real(wp), intent(in) :: radii(:)
         !> Whether the latest update applies
         logical :: current
      end function current_iface

      !> Yes/no property of a grid model
      !>
      !> @param[in] self Grid model
      pure function query_iface(self) result(yes)
         import :: grid_model_type
         implicit none(type, external)
         !> Grid model
         class(grid_model_type), intent(in) :: self
         !> Answer
         logical :: yes
      end function query_iface

      !> Components of every sphere after an update
      !>
      !> @param[in]  self  Updated grid model
      !> @param[in]  nsph  Number of spheres
      !> @param[out] comps Components per sphere and the distinct rotations
      subroutine components_iface(self, nsph, comps)
         import :: grid_model_type, grid_components_type
         implicit none(type, external)
         !> Grid model
         class(grid_model_type), intent(in) :: self
         !> Number of spheres
         integer, intent(in) :: nsph
         !> Components
         type(grid_components_type), intent(out) :: comps
      end subroutine components_iface

      !> Add the component adjoints to a nuclear gradient (reverse mode)
      !>
      !> @param[in]     self     Updated grid model
      !> @param[in]     torque   Torque per component (3, ncomp): sum over its
      !>                         points of (p - c_a) x dE/dp
      !> @param[in]     wbar     Weight adjoint dE/dw per component (ncomp)
      !> @param[in,out] gradient Nuclear gradient (3, nat), accumulated
      subroutine gradient_iface(self, torque, wbar, gradient)
         import :: grid_model_type, wp
         implicit none(type, external)
         !> Grid model
         class(grid_model_type), intent(in) :: self
         !> Torque per component
         real(wp), intent(in) :: torque(:, :)
         !> Weight adjoint per component
         real(wp), intent(in) :: wbar(:)
         !> Nuclear gradient
         real(wp), intent(inout) :: gradient(:, :)
      end subroutine gradient_iface

      !> Rotation and weight Jacobians per atom displacement (forward mode)
      !>
      !> @param[in]  self       Updated grid model
      !> @param[out] jac_rot    d omega_r(k)/dc_j(i) as (k, i, j, r), the angular
      !>                        velocity of every distinct rotation (3, 3, nat, nrot)
      !> @param[out] jac_weight   dw/dc_j(i) of every distinct weight as (i, j, class)
      !>                          (3, nat, nclass); left unallocated when no
      !>                          weight depends on the geometry
      !> @param[out] weight_class Distinct weight of every component, into
      !>                          `jac_weight` (ncomp); allocated with it
      subroutine jacobian_iface(self, jac_rot, jac_weight, weight_class)
         import :: grid_model_type, wp
         implicit none(type, external)
         !> Grid model
         class(grid_model_type), intent(in) :: self
         !> Rotation Jacobians
         real(wp), allocatable, intent(out) :: jac_rot(:, :, :, :)
         !> Weight Jacobians
         real(wp), allocatable, intent(out) :: jac_weight(:, :, :)
         !> Distinct weight of every component
         integer, allocatable, intent(out) :: weight_class(:)
      end subroutine jacobian_iface

      !> Declare the readable results of the model
      !>
      !> @param[in]     self  Grid model
      !> @param[in,out] query Walker collecting or fetching the declarations
      subroutine fields_iface(self, query)
         import :: grid_model_type, field_query_type
         implicit none(type, external)
         !> Grid model
         class(grid_model_type), intent(in) :: self
         !> Walker collecting or fetching the declarations
         type(field_query_type), intent(inout) :: query
      end subroutine fields_iface

      !> Print the report of the latest update
      !>
      !> @param[in] self  Grid model
      !> @param[in] unit  Output unit (default standard output)
      !> @param[in] level Verbosity level (default everything)
      subroutine print_iface(self, unit, level)
         import :: grid_model_type
         implicit none(type, external)
         !> Grid model
         class(grid_model_type), intent(in) :: self
         !> Output unit
         integer, intent(in), optional :: unit
         !> Verbosity level
         integer, intent(in), optional :: level
      end subroutine print_iface

   end interface

end module moist_cavity_grid_model_type
