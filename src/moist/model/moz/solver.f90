!> MOZ solver controls, construction and iteration interfaces
!>
!> JSON/TOML keys: name, picard.*, diis.*, hybrid.*, jfnk.*
module moist_model_moz_solver
   use mctc_env, only: error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_context, only: moist_context_type
   use moist_model_parameters, only: moist_model_parameters_type, print_hidden
   use moist_model_moz_solver_core_base, only: &
      & moist_model_moz_solver_type, moist_model_moz_solver_statistics_type
   use moist_model_moz_solver_core_evaluator, only: moist_model_moz_solver_evaluator_type
   use moist_model_moz_solver_core_iteration, only: moist_model_moz_solver_iteration_type, run_iterations
   use moist_model_moz_solver_core_printer, only: moist_model_moz_solver_printer_type
   use moist_model_moz_solver_core_preconditioner, only: &
      & moist_model_moz_solver_preconditioner_type, check_preconditioner_name
   use moist_model_moz_solver_jfnk_deflation, only: moist_model_moz_solver_deflation_type
   use moist_model_moz_solver_picard, only: &
      & moist_model_moz_solver_picard_parameters_type, new_picard_solver
   use moist_model_moz_solver_diis, only: moist_model_moz_solver_diis_type, &
      & moist_model_moz_solver_diis_parameters_type, new_diis_solver
   use moist_model_moz_solver_jfnk, only: moist_model_moz_solver_jfnk_type, &
      & moist_model_moz_solver_jfnk_parameters_type, moist_model_moz_solver_jfnk_adjoint_parameters_type, &
      & new_jfnk_solver
   use moist_model_moz_solver_hybrid, only: &
      & moist_model_moz_solver_hybrid_parameters_type, new_hybrid_solver
   use moist_model_moz_solver_jfnk_adjoint, only: solve_adjoint
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_type, moist_model_moz_solver_statistics_type
   public :: moist_model_moz_solver_evaluator_type, moist_model_moz_solver_iteration_type
   public :: moist_model_moz_solver_printer_type, moist_model_moz_solver_preconditioner_type
   public :: moist_model_moz_solver_deflation_type, moist_model_moz_solver_parameters_type
   public :: moist_model_moz_solver_picard_parameters_type, moist_model_moz_solver_diis_parameters_type
   public :: moist_model_moz_solver_hybrid_parameters_type, moist_model_moz_solver_jfnk_parameters_type
   public :: moist_model_moz_solver_jfnk_adjoint_parameters_type
   public :: new_solver, run_iterations, solve_adjoint, check_preconditioner_name

   !> Selected solver and the controls of every MOZ solver, one nested group per solver
   !>
   !> Files keep every group; printed tables only the groups the selected solver runs
   type, extends(moist_model_parameters_type) :: moist_model_moz_solver_parameters_type
      !> Solver: PICARD, DIIS, JFNK or HYBRID
      character(len=8) :: name = "JFNK"
      !> Picard damping, also used while DIIS fills its history
      type(moist_model_moz_solver_picard_parameters_type) :: picard
      !> DIIS history and mixing
      type(moist_model_moz_solver_diis_parameters_type) :: diis
      !> JFNK inner solve, forcing and line search
      type(moist_model_moz_solver_jfnk_parameters_type) :: jfnk
      !> Hybrid DIIS-to-JFNK switch
      type(moist_model_moz_solver_hybrid_parameters_type) :: hybrid
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Declare fields for JSON/TOML input, output and printing
      procedure :: register_entries => register_parameter_entries
      !> Check the solver name and the ranges of all controls
      procedure :: validate => validate_solver_parameters
      !> Whether the selected solver runs the controls of one group
      procedure :: runs => solver_runs
      !> Print level of settings the selected solver may not run
      procedure :: group_level => solver_group_level
   end type moist_model_moz_solver_parameters_type

contains

   !> Restore compiled defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_parameter_defaults(self)
      !> Parameter values
      class(moist_model_moz_solver_parameters_type), intent(inout) :: self
      type(moist_model_moz_solver_parameters_type) :: defaults

      self%name = defaults%name
      call self%picard%init_defaults()
      call self%diis%init_defaults()
      call self%jfnk%init_defaults()
      call self%hybrid%init_defaults()
   end subroutine init_parameter_defaults

   !> Declare parameter fields for JSON/TOML input, output and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_parameter_entries(self)
      !> Parameter values
      class(moist_model_moz_solver_parameters_type), intent(inout), target :: self

      call self%register_string("name", self%name, level=2)
      call self%register_group("picard", self%picard, level=self%group_level("picard", 2))
      call self%register_group("diis", self%diis, level=self%group_level("diis", 2))
      call self%register_group("jfnk", self%jfnk, level=self%group_level("jfnk", 2))
      call self%register_group("hybrid", self%hybrid, level=self%group_level("hybrid", 2))
   end subroutine register_parameter_entries

   !> Check the solver name and the ranges of all controls
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error Unknown solver or first invalid control
   subroutine validate_solver_parameters(self, error)
      !> Parameter values
      class(moist_model_moz_solver_parameters_type), intent(inout) :: self
      !> Unknown solver or first invalid control
      type(error_type), allocatable, intent(out) :: error

      if (len(solver_kind(self%name)) == 0) then
         call fatal_error(error, "Unknown MOZ solver: "//trim(self%name)//" (choose PICARD, DIIS, JFNK or HYBRID)")
         return
      end if
      call self%picard%validate(error)
      if (allocated(error)) return
      call self%diis%validate(error)
      if (allocated(error)) return
      call self%hybrid%validate(error)
      if (allocated(error)) return
      call self%jfnk%validate(error)
   end subroutine validate_solver_parameters

   !> Whether the selected solver runs the controls of one group
   !>
   !> - DIIS also runs the Picard damping while it fills its history
   !> - HYBRID runs every group
   !>
   !> @param[in] self Parameter values
   !> @param[in] group Group name: picard, diis, jfnk or hybrid
   pure logical function solver_runs(self, group) result(runs)
      !> Parameter values
      class(moist_model_moz_solver_parameters_type), intent(in) :: self
      !> Group name: picard, diis, jfnk or hybrid
      character(len=*), intent(in) :: group

      select case (solver_kind(self%name))
      case ("picard")
         runs = group == "picard"
      case ("diis")
         runs = group == "picard" .or. group == "diis"
      case ("jfnk")
         runs = group == "jfnk"
      case default
         ! HYBRID; unknown names are rejected by validate
         runs = .true.
      end select
   end function solver_runs

   !> Canonical solver kind of a name: picard, diis, jfnk or hybrid
   !>
   !> - Case-insensitive; MDIIS and Anderson are DIIS, Newton and Newton-Krylov are JFNK
   !> - Empty for an unknown name
   !>
   !> @param[in] name Solver name
   pure function solver_kind(name) result(kind)
      !> Solver name
      character(len=*), intent(in) :: name
      !> Canonical kind; empty when unknown
      character(len=:), allocatable :: kind

      kind = to_lower(trim(adjustl(name)))
      select case (kind)
      case ("picard", "diis", "jfnk", "hybrid")
         continue
      case ("mdiis", "anderson")
         kind = "diis"
      case ("newton", "newton-krylov")
         kind = "jfnk"
      case default
         kind = ""
      end select
   end function solver_kind

   !> Print level of settings that act only when the selected solver runs a group
   !>
   !> @param[in] self Parameter values
   !> @param[in] group Group name: picard, diis, jfnk or hybrid
   !> @param[in] level Print level while the group runs; `print_hidden` otherwise
   pure integer function solver_group_level(self, group, level) result(group_level)
      !> Parameter values
      class(moist_model_moz_solver_parameters_type), intent(in) :: self
      !> Group name: picard, diis, jfnk or hybrid
      character(len=*), intent(in) :: group
      !> Print level while the group runs
      integer, intent(in) :: level

      group_level = merge(level, print_hidden, self%runs(group))
   end function solver_group_level

   !> Allocate the solver named by its controls
   !>
   !> - Verbosity and OpenMP team size come from the borrowed run context
   !>
   !> @param[out] solver Allocated solver; unallocated on failure
   !> @param[in] ctx Run context, borrowed for the lifetime of the solver
   !> @param[out] error Unknown solver or invalid controls
   !> @param[in] parameters Solver name and numerical controls; JFNK defaults when absent
   subroutine new_solver(solver, ctx, error, parameters)
      !> Allocated solver; unallocated on failure
      class(moist_model_moz_solver_type), allocatable, intent(out) :: solver
      !> Run context, borrowed for the lifetime of the solver
      type(moist_context_type), intent(in), target :: ctx
      !> Unknown solver or invalid controls
      type(error_type), allocatable, intent(out) :: error
      !> Optional solver name and numerical controls
      type(moist_model_moz_solver_parameters_type), optional, intent(in) :: parameters
      !> Effective controls
      type(moist_model_moz_solver_parameters_type) :: p
      !> Hybrid components
      type(moist_model_moz_solver_diis_type) :: mdiis
      type(moist_model_moz_solver_jfnk_type) :: jfnk

      if (present(parameters)) p = parameters
      call p%validate(error)
      if (allocated(error)) return

      select case (solver_kind(p%name))
      case ("picard")
         allocate (solver, source=new_picard_solver(ctx, p%picard))
      case ("diis")
         allocate (solver, source=new_diis_solver(ctx, p%diis, picard=p%picard))
      case ("jfnk")
         allocate (solver, source=new_jfnk_solver(ctx, p%jfnk))
      case default
         ! HYBRID; validate rejected unknown names
         mdiis = new_diis_solver(ctx, p%diis, picard=p%picard)
         jfnk = new_jfnk_solver(ctx, p%jfnk)
         allocate (solver, source=new_hybrid_solver(ctx, mdiis, jfnk, p%hybrid))
      end select
   end subroutine new_solver

end module moist_model_moz_solver
