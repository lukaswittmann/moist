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

   !> Seed step: the surface mesh seeder's settings, derived from `proj_level`
   !> ([[seed_init_level]]; levels 2, 3, 4; failed local solves of 1)
   type :: drop_seed_parameters_type
      !> Sample spacing of the mesh (Bohr)
      real(wp) :: seed_size = 1.0_wp
      !> Lattice levels below the sample spacing for a cell without sign change;
      !> 0 leaves such a cell unrefined
      integer :: refine_levels = 0
      !> Seed filters: a seed is dropped
      !> outside the tangent cone of an anchor with an accepted local point
      logical :: filter_cone = .false.
      !> on a patch the anchor sees edge-on
      logical :: filter_grazing = .false.
      !> with a normal turned away from the anchor's sphere direction, for an
      !> anchor with an accepted local point
      logical :: filter_back_facing = .false.
      !> beyond both focal centers of its patch, where a point is a maximum of
      !> the distance
      logical :: filter_focal = .false.
      !> Tangent-cone scale in units of `h (1 + r k_max)`, the same at every level
      real(wp) :: cone_scale = 1.5_wp
   contains
      !> Spacing, refinement and seed filters of a projection level
      procedure :: init_level => seed_init_level
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
      !> Longest spatial Newton step before the line search (Bohr), also of the
      !> polish's Riemannian Newton steps
      real(wp) :: max_step = 1.0_wp
      !> Margin of the merit penalty over the multiplier, `mu >= |lambda+| + margin`
      real(wp) :: merit_margin = 1.0_wp
      !> Stop duplicates of a point the run already holds
      !>
      !> - Newton: a column whose prediction (full step, or flat step at
      !>   `lambda = 0`, unscaled) lands on a captor of its run, an accepted
      !>   column or a strict minimum converged in the same solve
      !> - polish: the escape walk stops on an accepted point of its run, and a
      !>   repeated visit of a saddle is skipped once the walks of its directions
      !>   are confirmed
      logical :: capture_duplicates = .false.
   end type drop_solve_parameters_type

   !> Optimality polish settings
   type :: drop_polish_parameters_type
      !> Iterations of the Riemannian saddle escape
      integer :: escape_max_iter = 150
      !> Initial line-search step of a non-saddle escape update; 1 is the full Newton step
      real(wp) :: escape_alpha_init = 1.0_wp
      !> Feasibility `|S|` for a KKT point
      real(wp) :: kkt_tol_feas = 1.0e-6_wp
      !> Stationarity `||w(x - a) - lambda grad S||` for a KKT point
      real(wp) :: kkt_tol_stat = 1.0e-6_wp
      !> Gradient norm floor for a well-defined tangent plane
      real(wp) :: kkt_tol_grad = 1.0e-10_wp
      !> Reduced-Hessian eigenvalue, relative to `||H_L||_F`, below whose negative a
      !> point is a saddle and above which a minimum; one threshold for check and
      !> escape, and, tenfold, the Newton captor test
      real(wp) :: saddle_eig_tol = 1.0e-8_wp
      !> Retraction iterations
      integer :: retract_max_iter = 50
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

   !> Seed-step preset of a projection level
   !>
   !> | level | spacing (Bohr) | refinement | cone | grazing | back-facing | focal |
   !> |-------|----------------|-----------|------|---------|-------------|-------|
   !> | 1     | 1.0            | 0         | no   | no      | no          | no    |
   !> | 2     | 0.5            | 0         | yes  | yes     | yes         | yes   |
   !> | 3     | 0.3            | 0         | yes  | no      | no          | no    |
   !> | 4     | 0.2            | 1         | no   | no      | no          | no    |
   !>
   !> - level 1 seeds only its failed local solves
   !> - level 2 seeds failed local solves, partial/closed gates, and fully open
   !>   local points beyond `0.50` Bohr
   !> - a back-facing seed leads to a point the horizon switch closes: with the
   !>   switch off, level 2 leaves the filter out
   !>
   !> @param[inout] self                Seed-step settings
   !> @param[in]    level               Projection level
   !> @param[in]    horizon_closes_back Whether the horizon switch closes back-facing points
   pure subroutine seed_init_level(self, level, horizon_closes_back)
      class(drop_seed_parameters_type), intent(inout) :: self
      integer, intent(in) :: level
      logical, intent(in) :: horizon_closes_back

      select case (level)
      case (4)
         self%seed_size = 0.2_wp
         self%refine_levels = 1
      case (3)
         self%seed_size = 0.3_wp
         self%refine_levels = 0
         self%filter_cone = .true.
      case (2)
         self%seed_size = 0.5_wp
         self%refine_levels = 0
         self%filter_cone = .true.
         self%filter_grazing = .true.
         self%filter_back_facing = .true. .and. horizon_closes_back
         self%filter_focal = .true.
      case default
         self%seed_size = 1.0_wp
         self%refine_levels = 0
      end select
      self%cone_scale = 1.5_wp
   end subroutine seed_init_level

end module moist_cavity_drop_projection_parameters
