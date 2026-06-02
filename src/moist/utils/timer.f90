!> Adapted from the timer implementation in xtb
!> Copyright (C) 2017-2020 Stefan Grimme.

module moist_utils_timer
   use mctc_env, only: wp, int64 => i8
   implicit none

   public :: timer_type
   private

   type :: timer_type

      !> number of timers
      integer, private :: n = 0

      !> verbosity
      logical, private :: verbose = .false.

      real(wp), private :: totwall = 0.0_wp
      real(wp), private :: totcpu = 0.0_wp
      logical, private, allocatable :: running(:)
      real(wp), private, allocatable :: twall(:)
      real(wp), private, allocatable :: tcpu(:)
      character(len=40), private, allocatable :: tag(:)
      integer, private, allocatable :: parent(:)

   contains

      procedure :: new => allocate_timer
      procedure :: delete => deallocate_timer
      procedure :: register => register_timer
      procedure :: register_parent => register_parent_timer
      procedure :: measure => timer
      procedure :: write_timing
      procedure :: write => write_all_timings
      procedure :: write_small => write_all_timings_small
      procedure :: get => get_timer
      procedure, private :: start_timing
      procedure, private :: stop_timing

   end type timer_type

contains

!> To initialize timer
   subroutine allocate_timer(self, n, verbose)

      implicit none

      !> instance of timer
      class(timer_type), intent(inout) :: self

      !> number of timers
      integer, intent(in)           :: n

      !> if verbose
      logical, intent(in), optional :: verbose

      real(wp) :: time_cpu
      real(wp) :: time_wall

      call self%delete

      ! capture negative values !
      if (n < 1) return

      self%n = n
      if (present(verbose)) self%verbose = verbose
      allocate (self%twall(0:n), source=0.0_wp)
      allocate (self%tcpu(0:n), source=0.0_wp)
      allocate (self%running(n), source=.false.)
      allocate (self%tag(n)); self%tag = ' '
      allocate (self%parent(n), source=0)

      ! launch timer !
      call self%start_timing(0)

   end subroutine allocate_timer

!> To deallocate memory
   subroutine deallocate_timer(self)

      implicit none

      !> instance of timer
      class(timer_type), intent(inout) :: self

      self%n = 0
      self%totwall = 0
      self%totcpu = 0
      self%verbose = .false.
      if (allocated(self%twall)) deallocate (self%twall)
      if (allocated(self%tcpu)) deallocate (self%tcpu)
      if (allocated(self%running)) deallocate (self%running)
      if (allocated(self%parent)) deallocate (self%parent)

   end subroutine deallocate_timer

!> To register a timer with a name and optional parent
!> If no parent is specified, timer is registered as a top-level (parent) timer
   subroutine register_timer(self, i, name, parent)

      implicit none

      !> instance of timer
      class(timer_type), intent(inout) :: self

      !> timer index
      integer, intent(in) :: i

      !> timer name
      character(len=*), intent(in) :: name

      !> parent timer index (0 = top-level, default if not provided)
      integer, intent(in), optional :: parent

      ! check bounds
      if (i < 1 .or. i > self%n) return

      ! set name
      self%tag(i) = trim(name)

      ! set parent: defaults to 0 (top-level) if not provided
      if (present(parent)) then
         if (parent >= 0 .and. parent < i) then
            self%parent(i) = parent
         end if
      else
         self%parent(i) = 0
      end if

   end subroutine register_timer

!> To register a parent (top-level) timer with a name
   subroutine register_parent_timer(self, i, name)

      implicit none

      !> instance of timer
      class(timer_type), intent(inout) :: self

      !> timer index
      integer, intent(in) :: i

      !> timer name
      character(len=*), intent(in) :: name

      ! Register with parent=0 (top-level)
      call self%register(i, name, parent=0)

   end subroutine register_parent_timer

!> To obtain current elapsed time
   function get_timer(self, i) result(time)

      !> instance of timer
      class(timer_type), intent(inout) :: self

      !> if specific timer
      integer, intent(in), optional :: i

      integer  :: it
      real(wp) :: tcpu, twall
      real(wp) :: time
      logical  :: running

      ! if i is not given, calculate overall elapsed time !
      if (present(i)) then
         it = i
      else
         it = 0
      end if

      if (it > 0) then
         running = self%running(it)
      else
         running = .true.
      end if

      if (running) then
         call timing(tcpu, twall)
         time = self%twall(it) + twall
      else
         time = self%twall(it)
      end if

   end function get_timer

!> To write timing for specific timer
   subroutine write_timing(self, iunit, i, inmsg, verbose)

      implicit none

      !> instance of timer
      class(timer_type), intent(inout) :: self

      !> I/O unit
      integer, intent(in) :: iunit

      !> index
      integer, intent(in) :: i

      !> raw message text
      character(len=*), intent(in), optional :: inmsg

      !> if verbose
      logical, intent(in), optional :: verbose

      character(len=26) :: msg
      real(wp) :: cputime, walltime
      integer(int64) ::  cpudays, cpuhours, cpumins
      integer(int64) :: walldays, wallhours, wallmins
      logical :: lverbose

!  '(1x,a,1x,"time:",1x,a)'
      ! check if tag should be added !
      if (present(inmsg)) then
         msg = inmsg
      else
         msg = self%tag(i)
      end if

      ! verbosity settings !
      if (present(verbose)) then
         lverbose = verbose
      else
         lverbose = self%verbose
      end if

      !           DAYS   HOURS   MINUTES   SECONDS
      ! DAYS        1     1/24    1/1440   1/86400
      ! HOURS      24      1       1/60     1/3600
      ! MINUTES   1440    60        1        1/60
      ! SECONDS  86400   3600      60         1

      ! convert elapsed CPU time into days, hours, minutes !
      cputime = self%tcpu(i)
      cpudays = int(cputime/86400._wp)
      cputime = cputime - cpudays*86400._wp
      cpuhours = int(cputime/3600._wp)
      cputime = cputime - cpuhours*3600._wp
      cpumins = int(cputime/60._wp)
      cputime = cputime - cpumins*60._wp

      ! convert elapsed wall time into days, hours, minutes !
      walltime = self%twall(i)
      walldays = int(walltime/86400._wp)
      walltime = walltime - walldays*86400._wp
      wallhours = int(walltime/3600._wp)
      walltime = walltime - wallhours*3600._wp
      wallmins = int(walltime/60._wp)
      walltime = walltime - wallmins*60._wp

      !----------!
      ! printout !
      !----------!

      if (lverbose) then
         write (iunit, '(1x,a)') msg
         write (iunit, '(" * wall-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
            walldays, wallhours, wallmins, walltime
         write (iunit, '(" *  cpu-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
            cpudays, cpuhours, cpumins, cputime
         write (iunit, '(1x,"*",1x,"ratio c/w:",1x,f9.3,1x,"speedup")') self%tcpu(i)/self%twall(i)
      else
         write (iunit, '(1x,a30,1x,"...",i9," min, ",f6.3," sec")') &
            msg, wallmins, walltime
      end if

   end subroutine write_timing

!> To write timing for all timers
   subroutine write_all_timings(self, iunit, inmsg)

      implicit none

      !> instance of timer
      class(timer_type), intent(inout) :: self

      !> I/O unit
      integer, intent(in) :: iunit

      !> raw message
      character(len=*), intent(in), optional :: inmsg

      character(len=26) :: msg
      real(wp) :: cputime, walltime
      integer  :: i
      integer(int64) ::  cpudays, cpuhours, cpumins
      integer(int64) :: walldays, wallhours, wallmins

      call self%stop_timing(0)

!  '(1x,a,1x,"time:",1x,a)'
      ! check if an external message should be added !
      if (present(inmsg)) then
         msg = inmsg//" timings"
      else
         msg = "total time"
      end if

      !           DAYS   HOURS   MINUTES   SECONDS
      ! DAYS        1     1/24    1/1440   1/86400
      ! HOURS      24      1       1/60     1/3600
      ! MINUTES   1440    60        1        1/60
      ! SECONDS  86400   3600      60         1

      ! convert overall elapsed CPU time into days, hours, minutes !
      cputime = self%tcpu(0)
      cpudays = int(cputime/86400._wp)
      cputime = cputime - cpudays*86400._wp
      cpuhours = int(cputime/3600._wp)
      cputime = cputime - cpuhours*3600._wp
      cpumins = int(cputime/60._wp)
      cputime = cputime - cpumins*60._wp

      ! convert overall elapsed wall time into days, hours, minutes !
      walltime = self%twall(0)
      walldays = int(walltime/86400._wp)
      walltime = walltime - walldays*86400._wp
      wallhours = int(walltime/3600._wp)
      walltime = walltime - wallhours*3600._wp
      wallmins = int(walltime/60._wp)
      walltime = walltime - wallmins*60._wp

      !----------!
      ! printout !
      !----------!
      write (iunit, '(a)')
      write (iunit, '(x,a)') "======================== T I M I N G S ======================="
      ! write(iunit,'(2x,a)') msg
      write (iunit, '(a)')

      ! printout every timer hierarchically
      do i = 1, self%n
         if (self%parent(i) == 0) then
            call print_timer_tree(self, iunit, i, 0)
         end if
      end do

      write (iunit, '(x,a)') "=============================================================="

      if (self%verbose) then
         write (iunit, '(2x, "Wall-time", 4x,i8," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
            walldays, wallhours, wallmins, walltime
         write (iunit, '(2x, "CPU-time", 5x,i8," d, ",i2," h, ",i2," min, ",f6.3," sec", 2x, "(",f0.1,"x)")') &
            cpudays, cpuhours, cpumins, cputime, self%tcpu(0)/self%twall(0)
         ! write(iunit,'(2x,"*",1x,"ratio c/w:",1x,f9.3,1x,"speedup")') self%tcpu (0)/self%twall(0)
      else
         write (iunit, '(2x,a26,i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
            msg, walldays, wallhours, wallmins, walltime
      end if

      write (iunit, '(x,a)') "=============================================================="

      write (iunit, '(a)')

   end subroutine write_all_timings

!> Helper to print timer tree recursively
   recursive subroutine print_timer_tree(self, iunit, i, depth)

      implicit none

      !> instance of timer
      type(timer_type), intent(in) :: self

      !> I/O unit
      integer, intent(in) :: iunit

      !> timer index
      integer, intent(in) :: i

      !> indentation depth
      integer, intent(in) :: depth

      real(wp) :: walltime, parent_time, percent
      integer(int64) :: wallmins
      integer :: j
      character(len=40) :: indent
      logical :: has_children

      ! Skip printing if timer was not activated (has no accumulated time)
      if (self%twall(i) < 1.0e-10_wp) return

      ! calculate indentation
      indent = repeat('  ', depth)

      ! calculate wall time
      walltime = self%twall(i)
      wallmins = int(walltime/60._wp)
      walltime = walltime - wallmins*60._wp

      ! check if this timer has children
      has_children = .false.
      do j = 1, self%n
         if (self%parent(j) == i) then
            has_children = .true.
            exit
         end if
      end do

      ! calculate percentage
      if (self%parent(i) == 0) then
         ! top-level timer: percentage of total
         parent_time = self%twall(0)
         if (parent_time > 1.0e-10_wp) then
            percent = 100*self%twall(i)/parent_time
         else
            percent = 0.0_wp
         end if
         write (iunit, '(x,a)') "=============================================================="
         write (iunit, '(2x,a,a20,4x,i9," min, ",f6.3," sec ",f7.3,"%")') &
            trim(indent), self%tag(i), wallmins, walltime, percent
         if (has_children) then
            write (iunit, '(x,a)') "--------------------------------------------------------------"
         end if
      else
         ! child timer: percentage of parent with parent name
         parent_time = self%twall(self%parent(i))
         if (parent_time > 1.0e-10_wp) then
            percent = 100*self%twall(i)/parent_time
         else
            percent = 0.0_wp
         end if
         write (iunit, '(3x,"-",x,a,a17,15x,f10.3," sec ",f7.3,"%")') &
            trim(indent), self%tag(i), wallmins*60.0_wp + walltime, percent!, trim(self%tag(self%parent(i)))
      end if

      ! print children
      do j = 1, self%n
         if (self%parent(j) == i) then
            call print_timer_tree(self, iunit, j, depth + 1)
         end if
      end do

   end subroutine print_timer_tree

!> To write timing for parent timers only (no hierarchy)
   subroutine write_all_timings_small(self, iunit, inmsg)

      implicit none

      !> instance of timer
      class(timer_type), intent(inout) :: self

      !> I/O unit
      integer, intent(in) :: iunit

      !> raw message
      character(len=*), intent(in), optional :: inmsg

      character(len=26) :: msg
      real(wp) :: cputime, walltime
      integer  :: i
      integer(int64) ::  cpudays, cpuhours, cpumins
      integer(int64) :: walldays, wallhours, wallmins
      real(wp) :: percent

      call self%stop_timing(0)

      ! check if an external message should be added !
      if (present(inmsg)) then
         msg = inmsg//" (total)"
      else
         msg = "total time"
      end if

      ! convert overall elapsed wall time into days, hours, minutes !
      walltime = self%twall(0)
      walldays = int(walltime/86400._wp)
      walltime = walltime - walldays*86400._wp
      wallhours = int(walltime/3600._wp)
      walltime = walltime - wallhours*3600._wp
      wallmins = int(walltime/60._wp)
      walltime = walltime - wallmins*60._wp

      !----------!
      ! printout !
      !----------!
      write (iunit, '(a)')
      write (iunit, '(1x,a26,i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
         msg, walldays, wallhours, wallmins, walltime

      ! printout parent timers only
      do i = 1, self%n
         if (self%parent(i) == 0) then
            walltime = self%twall(i)
            wallmins = int(walltime/60._wp)
            walltime = walltime - wallmins*60._wp
            percent = 100*self%twall(i)/self%twall(0)
            write (iunit, '(1x,a30,4x,i9," min, ",f6.3," sec (",f7.3,"%)")') &
               self%tag(i), wallmins, walltime, percent
         end if
      end do
      write (iunit, '(a)')

   end subroutine write_all_timings_small

!> start/stop button
   subroutine timer(self, i, inmsg)

      implicit none

      !> instance of timer
      class(timer_type), intent(inout) :: self

      !> index
      integer, intent(in) :: i

      !> raw message text
      character(len=*), intent(in), optional :: inmsg

      ! check if appropriate index is given !
      if (i > self%n .or. i < 1) return

      ! switcher between start/stop status !
      if (self%running(i)) then
         call self%stop_timing(i)
      else
         call self%start_timing(i)
      end if

      ! update status !
      self%running(i) = .not. self%running(i)

      ! assign tag to specific timer !
      if (present(inmsg)) self%tag(i) = trim(inmsg)

   end subroutine timer

!> To start counting
   subroutine start_timing(self, i)

      implicit none

      !> instance of timer
      class(timer_type), intent(inout) :: self

      !> index
      integer, intent(in) :: i

      real(wp) :: time_cpu
      real(wp) :: time_wall

      call timing(time_cpu, time_wall)
      self%tcpu(i) = self%tcpu(i) - time_cpu
      self%twall(i) = self%twall(i) - time_wall

   end subroutine start_timing

!> To stop counting
   subroutine stop_timing(self, i)

      implicit none

      !> instance of timer
      class(timer_type), intent(inout) :: self

      !> index
      integer, intent(in) :: i

      real(wp) :: time_cpu
      real(wp) :: time_wall

      call timing(time_cpu, time_wall)
      self%tcpu(i) = self%tcpu(i) + time_cpu
      self%twall(i) = self%twall(i) + time_wall

   end subroutine stop_timing

!> To retrieve the current CPU and wall time
   subroutine timing(time_cpu, time_wall)

      implicit none

      real(wp), intent(out) :: time_cpu
      real(wp), intent(out) :: time_wall

      !> current value of system clock (time passed from arbitary point)
      integer(int64) :: time_count

      !> number of clock ticks per second (conversion factor b/n ticks and seconds)
      integer(int64) :: time_rate
      integer(int64) :: time_max

      call system_clock(time_count, time_rate, time_max)
      call cpu_time(time_cpu)

      ! elapsed time in seconds !
      time_wall = real(time_count, wp)/real(time_rate, wp)

   end subroutine timing

end module moist_utils_timer
