!> Parameters of the DROP projection stages
!>
!> - pure parameters with their compiled defaults, held by
!>   [[moist_cavity_drop_parameters_type]] as `seed`, `solve`, `polish`,
!>   `finalize` and `branch`
!>   and read by the stage kernels as they are
!> - named by pipeline step (seed, solve, polish, finalize), not by algorithm;
!>   `branch` holds the branch model shared by the search and the weights
!> - a user setting translated into several values (`tolerance`,
!>   `wleb_prune_level`, `proj_level`) stays on the central type; its derived
!>   values land here
!> - the objective weight `phi_alpha`, the run context and values derived at
!>   worker init are passed to the kernels separately
module moist_cavity_drop_projection_parameters
   use mctc_env_accuracy, only: wp
   implicit none(type, external)
   private

   public :: drop_seed_parameters_type, drop_branch_parameters_type
   public :: drop_solve_parameters_type, drop_polish_parameters_type, drop_finalize_parameters_type

   !> Seed step: the octree surface coverage (levels 3, 4; fallback of 1)
   type :: drop_seed_parameters_type
      !> Leaf edge length of the octree (Bohr), derived from `proj_level`
      real(wp) :: seed_size = 1.0_wp
      !> Refinement levels of a surviving leaf without sign change, derived from
      !> `proj_level`; 0 drops such a leaf unrefined
      integer :: refine_levels = 0
      !> Tree-descent box budget per anchor; the level-synchronous search counts every
      !> box probed at every depth, so it trips earlier than the old lazy best-first
      !> search did for the same value; face-neighbour and refinement probes not counted
      integer :: max_boxes = 500000
      !> Surviving-leaf budget per anchor
      integer :: max_survivors = 200000
      !> Octree depth limit
      integer :: max_depth = 12
      !> Largest accepted `max_depth`, keeping the seeder's integer lattice
      !> coordinates representable (`max_addressable_depth` of the octree seeder)
      !>
      !> - not a configuration key
      integer :: max_addressable_depth = 20
   contains
      !> Coarse octree: levels 1 (fallback) and 3
      procedure :: init_coarse => seed_init_coarse
      !> Fine octree: level 4
      procedure :: init_fine => seed_init_fine
   end type drop_seed_parameters_type

   !> Branch model: the distance window of the branch weights, which also bounds
   !> the branch search behind a fully open point, and the reach
   type :: drop_branch_parameters_type
      !> Projection-distance excess over a fully open branch at which a
      !> branch loses all weight; width of the branch partition of unity (Bohr)
      !>
      !> - softened below `dpdist_max/2` from the surface, where the window
      !>   widens to `sqrt(2)*dpdist_max`
      real(wp) :: dpdist_max = 0.2_wp
      !> Projection distance up to which the reach gate of a point is open
      !> (Bohr; about 2 Angstrom, every level set)
      real(wp) :: reach_on = 3.78_wp
      !> Projection distance from which the reach gate is closed (Bohr; about
      !> 3 Angstrom)
      !>
      !> - no point beyond it carries weight, so it bounds the search of every anchor
      real(wp) :: reach_off = 5.67_wp
      !> Search every anchor to the full reach instead of the adaptive radius
      !>
      !> - debug switch, not a configuration key; same rows as the adaptive search
      logical :: full_reach = .false.
      !> Separation below which two points are one branch (Bohr), derived from `tolerance`
      real(wp) :: sep_cut = 1.0e-8_wp
      !> Lower bound of the derived `sep_cut` (Bohr)
      !>
      !> - not a configuration key
      real(wp) :: sep_floor = 1.0e-6_wp
   end type drop_branch_parameters_type

   !> Solve step: the batched KKT Newton
   type :: drop_solve_parameters_type
      !> Residual tolerance, derived from `tolerance`; `0.1 tol` is the step tolerance
      real(wp) :: tol = 1.0e-10_wp
      !> Step budget per point
      integer :: max_iter = 50
      !> `||F||` below which a point switches to full Newton
      real(wp) :: switch_fnorm = 1.0e-3_wp
      !> Cheap steps after which a point switches to full Newton
      integer :: max_cheap_steps = 2
      !> Smallest line-search step, accepted unconditionally
      real(wp) :: ls_floor = 0.1_wp
      !> Longest spatial Newton step before the line search (Bohr)
      real(wp) :: max_step = 1.0_wp
      !> Margin of the merit penalty over the multiplier, `mu >= |lambda+| + margin`
      real(wp) :: merit_margin = 1.0_wp
      !> Stop a column whose prediction lands on an earlier accepted column of its run
      logical :: capture_duplicates = .false.
   end type drop_solve_parameters_type

   !> Optimality polish settings
   type :: drop_polish_parameters_type
      !> Iterations of the Riemannian saddle escape
      integer :: escape_max_iter = 150
      !> Line-search trials per escape iteration
      integer :: escape_max_ls_iter = 20
      !> Initial line-search step of a non-saddle escape update; 1 is the full Newton step
      real(wp) :: escape_alpha_init = 1.0_wp
      !> Reduction of a rejected line-search step
      real(wp) :: escape_alpha_reduce = 0.75_wp
      !> Armijo coefficient of the escape line search
      real(wp) :: escape_armijo_c = 1.0e-4_wp
      !> Fixed step along negative curvature when leaving a saddle
      real(wp) :: escape_step = 0.5_wp
      !> Reduced-Hessian eigenvalue, relative to `||H_L||_F`, below whose negative a
      !> point is a saddle and above which a minimum; one threshold for check and escape
      real(wp) :: saddle_eig_tol = 1.0e-8_wp
      !> Retraction iterations
      integer :: retract_max_iter = 50
      !> Report a saddle as an error when the caller does not ask for `is_saddle`
      logical :: strict_minimum_required = .true.
   end type drop_polish_parameters_type

   !> Per-anchor finalize settings
   type :: drop_finalize_parameters_type
      !> Quadrature weight cutoff, derived from `tolerance`
      real(wp) :: wleb_cut = 5.0e-12_wp
      !> Prune switch on the geometric weight, derived from `wleb_prune_level`
      logical :: wleb_prune = .false.
      !> Off-edge of the prune switch (below: fully off), derived from `wleb_prune_level`
      real(wp) :: wleb_prune_off = 0.0_wp
      !> On-edge of the prune switch (above: fully on), derived from `wleb_prune_level`
      real(wp) :: wleb_prune_on = 0.0_wp
   end type drop_finalize_parameters_type

contains

   !> Coarse octree preset: 1 Bohr leaves, a leaf without sign change is dropped
   !>
   !> - sets the leaf size and the refinement only; the budgets stay
   !>
   !> @param[inout] self Seed-step settings
   pure subroutine seed_init_coarse(self)
      class(drop_seed_parameters_type), intent(inout) :: self

      self%seed_size = 1.0_wp
      self%refine_levels = 0
   end subroutine seed_init_coarse

   !> Fine octree preset: 0.2 Bohr leaves, a leaf without sign change has its
   !> children probed once
   !>
   !> - sets the leaf size and the refinement only; the budgets stay
   !>
   !> @param[inout] self Seed-step settings
   pure subroutine seed_init_fine(self)
      class(drop_seed_parameters_type), intent(inout) :: self

      self%seed_size = 0.2_wp
      self%refine_levels = 1
   end subroutine seed_init_fine

end module moist_cavity_drop_projection_parameters
