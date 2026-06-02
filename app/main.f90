

!> Entry point for the command line interface of moist
program driver
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : error_type
   use moist_cli, only : run_config, get_arguments
   use moist_driver, only : main
   implicit none (type, external)
   !> Configuration data deteriming the driver behaviour
   type(run_config) :: config
   !> Error handling
   type(error_type), allocatable :: error

   call get_arguments(config, error)
   if (allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call main(config, error)
   if (allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

end program driver
