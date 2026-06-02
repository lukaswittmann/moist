!> General-purpose CLI histogram for diagnostic visualization.
!>
!> Provides a derived type that bins real-valued data and renders
!> ASCII histograms to the terminal. Supports optional log-scaled
!> Y-axis, configurable dimensions, and summary statistics.
module moist_utils_histogram
   use mctc_env, only: wp, error_type
   use moist_math_sorter_quicksort, only: qsort
   use iso_fortran_env, only: output_unit
   implicit none
   private

   public :: histogram_type

   type :: histogram_type
      !> Number of bins
      integer :: nbins = 0
      !> Lower bound of the histogram range
      real(wp) :: xmin = 0.0_wp
      !> Upper bound of the histogram range
      real(wp) :: xmax = 0.0_wp
      !> Width of each bin
      real(wp) :: bin_width = 0.0_wp
      !> Bin counts (nbins)
      integer, allocatable :: counts(:)
      !> Bin edges (nbins + 1)
      real(wp), allocatable :: bin_edges(:)
      !> Total number of data points inserted
      integer :: n_total = 0
      !> Number of points below xmin
      integer :: n_underflow = 0
      !> Number of points above xmax
      integer :: n_overflow = 0
      !> Minimum value in data
      real(wp) :: data_min = 0.0_wp
      !> Maximum value in data
      real(wp) :: data_max = 0.0_wp
      !> Mean value of data
      real(wp) :: data_mean = 0.0_wp
      !> Median value of data
      real(wp) :: data_median = 0.0_wp
      !> Optional label / title for the histogram
      character(:), allocatable :: label
      !> Whether data has been inserted and binned
      logical :: filled = .false.
   contains
      !> Insert data, compute bins and statistics
      procedure :: insert => histogram_insert
      !> Render ASCII histogram to a Fortran unit
      procedure :: print => histogram_print
      !> Reset the histogram to its initial state
      procedure :: clear => histogram_clear
   end type histogram_type

contains

   !> Insert data into the histogram, computing bins and statistics.
   !>
   !> Clears any previous state, computes summary statistics (min, max,
   !> mean, median), builds uniform bin edges over [xmin, xmax], and
   !> bins each data value.
   !>
   !> @param[inout] self  Histogram instance
   !> @param[in]    data  Array of values to bin
   !> @param[in]    nbins Number of bins (default: 20)
   !> @param[in]    xmin  Lower range bound (default: data minimum)
   !> @param[in]    xmax  Upper range bound (default: data maximum)
   !> @param[in]    label Title string for display
   subroutine histogram_insert(self, data, nbins, xmin, xmax, label)
      !> Histogram instance
      class(histogram_type), intent(inout) :: self
      !> Array of values to bin
      real(wp), intent(in) :: data(:)
      !> Number of bins (default: 20)
      integer, intent(in), optional :: nbins
      !> Lower range bound (default: data minimum)
      real(wp), intent(in), optional :: xmin
      !> Upper range bound (default: data maximum)
      real(wp), intent(in), optional :: xmax
      !> Title string for display
      character(*), intent(in), optional :: label

      integer :: n, nb, i, ibin
      real(wp) :: x, range_min, range_max, bw
      real(wp), allocatable :: sorted(:)
      type(error_type), allocatable :: error

      call self%clear()

      n = size(data)
      if (n == 0) return

      ! Resolve number of bins
      nb = 20
      if (present(nbins)) nb = max(1, nbins)

      ! Compute data statistics
      self%data_min = minval(data)
      self%data_max = maxval(data)
      self%data_mean = sum(data)/real(n, wp)

      ! Compute median via sorted copy
      allocate (sorted(n))
      sorted = data
      call qsort(sorted, error=error)
      if (mod(n, 2) == 1) then
         self%data_median = sorted((n + 1)/2)
      else
         self%data_median = 0.5_wp*(sorted(n/2) + sorted(n/2 + 1))
      end if
      deallocate (sorted)

      ! Resolve range
      range_min = self%data_min
      if (present(xmin)) range_min = xmin
      range_max = self%data_max
      if (present(xmax)) range_max = xmax

      ! Guard against degenerate range (all values identical)
      if (range_max <= range_min) then
         if (range_min == 0.0_wp) then
            range_min = -0.5_wp
            range_max = 0.5_wp
         else
            range_min = range_min - 0.5_wp*abs(range_min)
            range_max = range_max + 0.5_wp*abs(range_max)
         end if
      end if

      ! Store parameters
      self%nbins = nb
      self%xmin = range_min
      self%xmax = range_max
      self%n_total = n
      self%bin_width = (range_max - range_min)/real(nb, wp)
      if (present(label)) self%label = trim(label)

      ! Allocate bin arrays
      allocate (self%bin_edges(nb + 1))
      allocate (self%counts(nb), source=0)

      ! Compute bin edges
      bw = self%bin_width
      do i = 1, nb + 1
         self%bin_edges(i) = range_min + real(i - 1, wp)*bw
      end do

      ! Bin the data: [edge_i, edge_{i+1}) except last bin which is closed
      do i = 1, n
         x = data(i)
         if (x < range_min) then
            self%n_underflow = self%n_underflow + 1
         else if (x > range_max) then
            self%n_overflow = self%n_overflow + 1
         else if (x == range_max) then
            ! Place exactly-at-max in last bin
            self%counts(nb) = self%counts(nb) + 1
         else
            ibin = int(floor((x - range_min)/bw)) + 1
            ibin = max(1, min(nb, ibin))
            self%counts(ibin) = self%counts(ibin) + 1
         end if
      end do

      self%filled = .true.

   end subroutine histogram_insert

   !> Render the histogram as ASCII art to a Fortran output unit.
   !>
   !> Draws a vertical bar chart with configurable width, height,
   !> optional log-scaled Y-axis, and summary statistics.
   !>
   !> @param[in] self       Histogram instance (after insert)
   !> @param[in] unit       Fortran unit number (default: output_unit)
   !> @param[in] width      Total character width (default: 60)
   !> @param[in] height     Number of bar rows (default: 20)
   !> @param[in] log_y      Use log10 scale on Y-axis (default: .false.)
   !> @param[in] show_stats Print summary statistics (default: .true.)
   subroutine histogram_print(self, unit, width, height, log_y, show_stats)
      !> Histogram instance
      class(histogram_type), intent(in) :: self
      !> Fortran unit number (default: output_unit)
      integer, intent(in), optional :: unit
      !> Total character width (default: 60)
      integer, intent(in), optional :: width
      !> Number of bar rows (default: 20)
      integer, intent(in), optional :: height
      !> Use log10 scale on Y-axis (default: .false.)
      logical, intent(in), optional :: log_y
      !> Print summary statistics below (default: .true.)
      logical, intent(in), optional :: show_stats

      integer :: iu, plot_width, plot_height
      logical :: use_log, print_stats
      integer :: max_count, y_label_width, bar_area_width
      integer :: irow, ibin, chars_per_bin
      real(wp) :: y_max, threshold
      real(wp), allocatable :: display_vals(:)
      character(32) :: buf
      character(256) :: line
      integer :: pos

      ! X-axis tick layout
      integer :: n_ticks, itick, tick_pos, label_len, label_start
      real(wp) :: tick_val, mag
      character(12) :: tick_fmt
      character(16) :: tick_label
      integer, allocatable :: tick_positions(:)

      ! Resolve optionals
      iu = output_unit
      if (present(unit)) iu = unit
      plot_width = 60
      if (present(width)) plot_width = max(20, width)
      plot_height = 20
      if (present(height)) plot_height = max(5, height)
      use_log = .false.
      if (present(log_y)) use_log = log_y
      print_stats = .true.
      if (present(show_stats)) print_stats = show_stats

      ! Guard: no data
      if (.not. self%filled) then
         write (iu, '(2x,a)') '[Histogram] No data inserted.'
         return
      end if

      ! Compute display values (raw or log-scaled)
      max_count = maxval(self%counts)
      allocate (display_vals(self%nbins))

      if (use_log) then
         do ibin = 1, self%nbins
            if (self%counts(ibin) > 0) then
               display_vals(ibin) = log10(real(self%counts(ibin), wp))
            else
               display_vals(ibin) = 0.0_wp
            end if
         end do
         if (max_count > 0) then
            y_max = log10(real(max_count, wp))
         else
            y_max = 0.0_wp
         end if
      else
         do ibin = 1, self%nbins
            display_vals(ibin) = real(self%counts(ibin), wp)
         end do
         y_max = real(max_count, wp)
      end if

      if (y_max <= 0.0_wp) then
         write (iu, '(2x,a)') '[Histogram] All bins empty.'
         return
      end if

      ! Compute layout dimensions
      write (buf, '(i0)') max_count
      y_label_width = len_trim(buf) + 1

      bar_area_width = plot_width - y_label_width - 2
      bar_area_width = max(self%nbins, bar_area_width)
      chars_per_bin = max(1, bar_area_width/self%nbins)
      bar_area_width = chars_per_bin*self%nbins

      ! === Title ===
      if (allocated(self%label)) then
         write (iu, '(2x,a)') self%label
      end if
      if (use_log) then
         write (iu, '(2x,a)') 'Count (log10)'
      else
         write (iu, '(2x,a)') 'Count'
      end if

      ! === Histogram rows (top to bottom) ===
      do irow = plot_height, 1, -1
         threshold = y_max*real(irow, wp)/real(plot_height, wp)

         ! Y-axis label on selected rows
         line = ''
         if (irow == plot_height .or. irow == 1 .or. &
             irow == (plot_height + 1)/2) then
            if (use_log) then
               write (buf, '(f6.1)') threshold
               buf = adjustr(buf(1:y_label_width))
            else
               write (buf, '(i0)') nint(threshold)
               buf = adjustr(buf(1:y_label_width))
            end if
            line(1:y_label_width) = buf(1:y_label_width)
         else
            line(1:y_label_width) = repeat(' ', y_label_width)
         end if

         ! Separator
         pos = y_label_width + 1
         line(pos:pos + 1) = '| '
         pos = pos + 2

         ! Bar characters
         do ibin = 1, self%nbins
            if (display_vals(ibin) >= threshold) then
               line(pos:pos + chars_per_bin - 1) = repeat('#', chars_per_bin)
            else
               line(pos:pos + chars_per_bin - 1) = repeat(' ', chars_per_bin)
            end if
            pos = pos + chars_per_bin
         end do

         write (iu, '(a)') trim(line)
      end do

      ! === X-axis with tick marks ===

      ! Choose a consistent format based on the magnitude of the range
      mag = max(abs(self%xmin), abs(self%xmax))
      if (mag >= 0.01_wp .and. mag < 1000.0_wp) then
         tick_fmt = '(f8.2)'
         label_len = 8
      else if (mag >= 1000.0_wp .and. mag < 100000.0_wp) then
         tick_fmt = '(f8.0)'
         label_len = 8
      else
         tick_fmt = '(es9.1)'
         label_len = 9
      end if

      ! Determine how many ticks fit (need at least label_len+1 chars between)
      n_ticks = max(2, bar_area_width/(label_len + 1))
      ! Cap at a reasonable number
      n_ticks = min(n_ticks, 8)

      ! Compute tick character positions within the bar area (0-based from bar start)
      allocate (tick_positions(n_ticks))
      do itick = 1, n_ticks
         tick_positions(itick) = nint(real((itick - 1)*bar_area_width, wp) &
                                      /real(n_ticks - 1, wp))
      end do
      ! Ensure last tick is exactly at the end
      tick_positions(n_ticks) = bar_area_width

      ! Print the axis line: '-' with '+' at tick positions
      line = ''
      line(1:y_label_width) = repeat(' ', y_label_width)
      pos = y_label_width + 1
      line(pos:pos) = '+'
      pos = pos + 1
      do ibin = 0, bar_area_width
         if (any(tick_positions == ibin)) then
            line(pos + ibin:pos + ibin) = '+'
         else
            line(pos + ibin:pos + ibin) = '-'
         end if
      end do
      write (iu, '(a)') trim(line)

      ! Print tick labels below, centered on each tick position
      line = ''
      do itick = 1, n_ticks
         tick_val = self%xmin + real(tick_positions(itick), wp) &
                    /real(bar_area_width, wp)*(self%xmax - self%xmin)
         write (tick_label, tick_fmt) tick_val
         tick_label = adjustl(tick_label)

         ! Center the label on the tick position
         tick_pos = y_label_width + 2 + tick_positions(itick)
         label_start = tick_pos - len_trim(tick_label)/2
         label_start = max(1, label_start)

         ! Write label if it fits without overwriting previous labels
         if (label_start + len_trim(tick_label) - 1 <= 256 .and. &
             line(label_start:label_start) == ' ') then
            line(label_start:label_start + len_trim(tick_label) - 1) = &
               trim(tick_label)
         end if
      end do
      write (iu, '(a)') trim(line)

      deallocate (tick_positions)

      ! === Statistics ===
      if (print_stats) then
         write (iu, '(2x,a,i0,2x,a,es10.3,2x,a,es10.3)') &
            'N=', self%n_total, 'min=', self%data_min, 'max=', self%data_max
         write (iu, '(2x,a,es10.3,2x,a,es10.3)') &
            'mean=', self%data_mean, 'median=', self%data_median
      end if

      ! === Overflow / underflow ===
      if (self%n_underflow > 0 .or. self%n_overflow > 0) then
         write (iu, '(2x,a,i0,2x,a,i0)') &
            'underflow=', self%n_underflow, 'overflow=', self%n_overflow
      end if

      write (iu, '(a)') ''
      deallocate (display_vals)

   end subroutine histogram_print

   !> Reset the histogram to its initial empty state.
   !>
   !> Deallocates all arrays and resets counters.
   !>
   !> @param[inout] self Histogram instance
   subroutine histogram_clear(self)
      !> Histogram instance
      class(histogram_type), intent(inout) :: self

      self%nbins = 0
      self%xmin = 0.0_wp
      self%xmax = 0.0_wp
      self%bin_width = 0.0_wp
      self%n_total = 0
      self%n_underflow = 0
      self%n_overflow = 0
      self%data_min = 0.0_wp
      self%data_max = 0.0_wp
      self%data_mean = 0.0_wp
      self%data_median = 0.0_wp
      self%filled = .false.
      if (allocated(self%counts)) deallocate (self%counts)
      if (allocated(self%bin_edges)) deallocate (self%bin_edges)
      if (allocated(self%label)) deallocate (self%label)

   end subroutine histogram_clear

end module moist_utils_histogram
