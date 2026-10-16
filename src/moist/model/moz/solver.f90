!> MOZ solver construction and iteration interfaces
module moist_model_moz_solver
   use moist_model_moz_solver_base, only: &
      & moist_model_moz_solver_type, moist_model_moz_solver_statistics_type
   use moist_model_moz_solver_evaluator, only: moist_model_moz_solver_evaluator_type
   use moist_model_moz_solver_iteration, only: moist_model_moz_solver_iteration_type, run_iterations
   use moist_model_moz_solver_printer, only: moist_model_moz_solver_printer_type
   use moist_model_moz_solver_preconditioner, only: moist_model_moz_solver_preconditioner_type
   use moist_model_moz_solver_deflation, only: moist_model_moz_solver_deflation_type
   use moist_model_moz_solver_parameters, only: moist_model_moz_solver_parameters_type
   use moist_model_moz_solver_factory, only: new_solver
   use moist_model_moz_solver_adjoint, only: solve_adjoint
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_type, moist_model_moz_solver_statistics_type
   public :: moist_model_moz_solver_evaluator_type, moist_model_moz_solver_iteration_type
   public :: moist_model_moz_solver_printer_type, moist_model_moz_solver_preconditioner_type
   public :: moist_model_moz_solver_deflation_type, moist_model_moz_solver_parameters_type
   public :: new_solver, run_iterations, solve_adjoint

end module moist_model_moz_solver
