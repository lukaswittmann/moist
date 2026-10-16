!> Construction of Picard, DIIS, JFNK and hybrid solvers
module moist_model_moz_solver_factory
   use mctc_env, only: error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_moz_solver_base, only: moist_model_moz_solver_type
   use moist_model_moz_solver_parameters, only: moist_model_moz_solver_parameters_type
   use moist_model_moz_solver_picard, only: new_picard_solver
   use moist_model_moz_solver_diis, only: moist_model_moz_solver_diis_type, new_diis_solver
   use moist_model_moz_solver_jfnk, only: moist_model_moz_solver_jfnk_type, new_jfnk_solver
   use moist_model_moz_solver_hybrid, only: new_hybrid_solver
   implicit none(type, external)
   private

   public :: new_solver

contains

   !> Allocate a solver from its name and optional numerical controls
   !>
   !> @param[out] solver Allocated solver; unallocated on failure
   !> @param[in] name Solver identifier: PICARD, DIIS, JFNK or HYBRID
   !> @param[out] error Unknown identifier or invalid controls
   !> @param[in] parameters Numerical controls; defaults when absent
   recursive subroutine new_solver(solver, name, error, parameters)
      !> Allocated solver; unallocated on failure
      class(moist_model_moz_solver_type), allocatable, intent(out) :: solver
      !> Solver identifier
      character(len=*), intent(in) :: name
      !> Unknown identifier or invalid controls
      type(error_type), allocatable, intent(out) :: error
      !> Optional numerical controls
      type(moist_model_moz_solver_parameters_type), optional, intent(in) :: parameters
      !> Effective controls
      type(moist_model_moz_solver_parameters_type) :: p
      !> Configure hybrid components through the same factory paths
      class(moist_model_moz_solver_type), allocatable :: mdiis, jfnk

      if (present(parameters)) p = parameters
      call p%validate(error)
      if (allocated(error)) return

      select case (to_lower(trim(adjustl(name))))
      case ("picard")
         allocate (solver, source=new_picard_solver(damping=p%picard_damping, verbosity=p%verbosity))
      case ("diis", "mdiis")
         allocate (solver, source=new_diis_solver(depth=p%diis_depth, damping=p%diis_damping, &
            & picard_damping=p%picard_damping, backtracking=p%diis_backtracking, &
            & backtrack_tau=p%diis_backtrack_tau, min_alpha=p%diis_min_alpha, &
            & max_step_ratio=p%diis_max_step_ratio, verbosity=p%verbosity))
      case ("jfnk")
         allocate (solver, source=new_jfnk_solver(restart=p%jfnk_restart, max_krylov=p%jfnk_max_krylov, &
            & eta0=p%jfnk_eta0, eta_min=p%jfnk_eta_min, eta_max=p%jfnk_eta_max, fd_eps=p%jfnk_fd_eps, &
            & verbosity=p%verbosity, restart_stagnation_tol=p%jfnk_stagnation_tol, cgs2=p%jfnk_cgs2, &
            & exact_jvp=p%jfnk_exact_jvp, recycle=p%jfnk_recycle, deflation=p%jfnk_deflation, &
            & deflation_rmse=p%jfnk_deflation_rmse, alpha_max=p%jfnk_alpha_max))
      case ("hybrid")
         call new_solver(mdiis, "DIIS", error, p)
         if (allocated(error)) return
         call new_solver(jfnk, "JFNK", error, p)
         if (allocated(error)) return
         select type (mdiis)
         type is (moist_model_moz_solver_diis_type)
            select type (jfnk)
            type is (moist_model_moz_solver_jfnk_type)
               allocate (solver, source=new_hybrid_solver(mdiis, jfnk, p%hybrid_switch_rmse))
            end select
         end select
      case default
         call fatal_error(error, "Unknown MOZ solver: "//trim(name)//" (choose PICARD, DIIS, JFNK or HYBRID)")
      end select
      ! Hybrid components already carry the count from their own factory calls
      if (allocated(solver)) solver%nthreads = p%nthreads
   end subroutine new_solver

end module moist_model_moz_solver_factory
