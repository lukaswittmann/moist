!> Environment-variable lookup and output-directory resolution helpers
module moist_utils_env

   implicit none
   private

   public :: get_env
   public :: resolve_dir
   public :: ensure_dir

contains

   !> Read an environment variable into an allocatable string
   !>
   !> @param[in]  name     Environment variable name to query
   !> @param[in]  default  Optional: value used when the variable is unset or empty
   !> @return              Trimmed value, `default`, or a zero-length string
   function get_env(name, default) result(val)
      !> Environment variable name to query
      character(len=*), intent(in) :: name
      !> Value substituted when the variable is unset or empty
      character(len=*), intent(in), optional :: default
      !> Resolved value; always allocated
      character(len=:), allocatable :: val
      !> Length of the variable value as reported by the runtime
      integer :: length
      !> Query status (0 if the variable exists)
      integer :: status

      call get_environment_variable(name, length=length, status=status)
      if (status == 0 .and. length > 0) then
         allocate (character(len=length) :: val)
         call get_environment_variable(name, value=val)
         val = trim(val)
         if (len(val) > 0) return
      end if

      val = ""
      if (present(default)) val = default
   end function get_env

   !> Resolve a directory from an environment variable, falling back to a default.
   !> Does not touch the filesystem; pair with `ensure_dir` to create it.
   !> @param[in]  env_name  Environment variable consulted first
   !> @param[in]  fallback  Default used when the variable is unset or empty
   !> @return               Resolved directory path (no trailing slash)
   function resolve_dir(env_name, fallback) result(dir)
      !> Environment variable consulted first
      character(len=*), intent(in) :: env_name
      !> Default used when the variable is unset or empty
      character(len=*), intent(in) :: fallback
      !> Resolved directory path
      character(len=:), allocatable :: dir

      dir = get_env(env_name, default=fallback)
   end function resolve_dir

   !> Create a directory (and any missing parents) if it does not already exist.
   !> @param[in]  path  Directory to create
   !> @param[out] stat  Optional status; 0 on success, non-zero on failure
   subroutine ensure_dir(path, stat)
      !> Directory to create
      character(len=*), intent(in) :: path
      !> Optional status (0 on success)
      integer, intent(out), optional :: stat
      !> Status of launching the shell command
      integer :: cmdstat
      !> Exit status of the mkdir command
      integer :: exitstat

      call execute_command_line("mkdir -p '"//trim(path)//"'", &
                                cmdstat=cmdstat, exitstat=exitstat)
      if (present(stat)) then
         if (cmdstat /= 0) then
            stat = cmdstat
         else
            stat = exitstat
         end if
      end if
   end subroutine ensure_dir

end module moist_utils_env
