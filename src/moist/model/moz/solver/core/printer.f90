!> Shared iteration table for MOZ solvers
!>
!> - Drivers supply outer diagnostics
!> - JFNK writes the inner GMRES rows
!> - Each solver selects its own columns through `set_columns` when attached
!> - The header row prints before the first row of a stage and after a column change
module moist_model_moz_solver_core_printer
   use mctc_env, only: wp
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_printer_type

   !> Columns of every solver: Iter, t (s), lambda, F, dF, RMSR, max|dg|, max(h), alpha
   !>
   !> - lambda only for staged solves; F and dF only with energy columns
   integer, parameter :: NCOLS_COMMON = 9
   !> Position of the lambda column
   integer, parameter :: COL_LAMBDA = 3
   !> Position of the energy column; the energy-change column follows
   integer, parameter :: COL_ENERGY = 4
   !> Largest energy magnitude in fixed-point format; larger values switch to ES
   real(wp), parameter :: max_fixed_energy = 1.0e7_wp
   !> DIIS column: ndiis
   integer, parameter :: NCOLS_DIIS = 1
   !> GMRES block columns: gmres, krylov, res, tol
   !>
   !> - krylov, res and tol only while inner rows are printed
   integer, parameter :: NCOLS_GMRES = 4
   !> Total table columns
   integer, parameter :: NCOLS_TOTAL = NCOLS_COMMON + NCOLS_DIIS + NCOLS_GMRES
   !> Widths of every column, in table order
   integer, parameter :: widths(NCOLS_TOTAL) = [6, 7, 9, 12, 10, 11, 11, 8, 7, 6, 6, 7, 11, 11]

   !> Iteration table writer shared by the outer driver and the inner GMRES
   !>
   !> - Picard shows the common columns, DIIS adds ndiis, JFNK the GMRES block
   type :: moist_model_moz_solver_printer_type
      !> Underlying table printer, set up by `setup` and rebuilt by `set_columns`
      type(prettylistprinter) :: plp
      !> Headers of every column, with the energy labels of `setup`
      character(8) :: headers(NCOLS_TOTAL) = ""
      !> Whether the DIIS column and the GMRES block are shown
      logical :: diis_column = .true., gmres_block = .true.
      !> Whether inner GMRES rows are printed; shows the krylov, res and tol columns
      logical :: inner_rows = .false.
      !> Whether the lambda column of staged solves is shown
      logical :: lambda_column = .false.
      !> Whether the energy and energy-change columns are shown
      logical :: energy_columns = .true.
      !> Whether the header row prints before the next row
      logical :: header_due = .true.
      !> Energy of the previous row of the stage; huge when there is none
      real(wp) :: energy_prev = huge(1.0_wp)
   contains
      procedure :: setup => printer_setup
      procedure :: set_columns => printer_set_columns
      procedure :: begin_stage => printer_begin_stage
      procedure :: header => printer_header
      procedure :: note => printer_note
      procedure :: sep_note => printer_sep_note
      procedure :: print_header => printer_print_header
      procedure :: separator => printer_separator
      procedure :: blank => printer_blank
      procedure :: stat_line => printer_stat_line
      procedure :: print_outer_row => printer_print_outer_row
      procedure :: print_inner_gmres_row => printer_print_inner_gmres_row
   end type moist_model_moz_solver_printer_type

contains

   !> Initialise the printer with every column, required before any print method
   !>
   !> - Optional energy_label/denergy_label override the default "F"/"dF" headers
   !>
   !> @param[out] self           printer to set up
   !> @param[in]  unit           output unit; absent takes the table printer default
   !> @param[in]  energy_label   header of the energy column; default "F"
   !> @param[in]  denergy_label  header of the energy-change column; default "dF"
   !> @param[in]  inner_rows     print inner GMRES rows and their columns; default off
   !> @param[in]  lambda_column  show the lambda column of staged solves; default off
   !> @param[in]  energy_columns show the energy and energy-change columns; default on
   subroutine printer_setup(self, unit, energy_label, denergy_label, inner_rows, lambda_column, energy_columns)
      !> Printer to set up
      class(moist_model_moz_solver_printer_type), intent(out) :: self
      !> Output unit; absent takes the table printer default
      integer, intent(in), optional :: unit
      !> Headers of the energy and energy-change columns; defaults "F" and "dF"
      character(*), intent(in), optional :: energy_label, denergy_label
      !> Print inner GMRES rows and their columns
      logical, intent(in), optional :: inner_rows
      !> Show the lambda column of staged solves
      logical, intent(in), optional :: lambda_column
      !> Show the energy and energy-change columns
      logical, intent(in), optional :: energy_columns

      self%headers = [character(8) :: "Iter", "t (s)", "lambda", "F", "dF", "RMSR", &
         & "max|dg|", "max(h)", "alpha", "ndiis", "gmres", "krylov", "res", "tol"]
      if (present(energy_label)) self%headers(4) = energy_label
      if (present(denergy_label)) self%headers(5) = denergy_label
      if (present(unit)) self%plp%unit = unit
      if (present(inner_rows)) self%inner_rows = inner_rows
      if (present(lambda_column)) self%lambda_column = lambda_column
      if (present(energy_columns)) self%energy_columns = energy_columns
      call build_table(self)
   end subroutine printer_setup

   !> Show the DIIS column and the GMRES block of the solver in use
   !>
   !> - Solvers call it when a printer is attached and when their columns change
   !> - A change rebuilds the table on the same unit; the header row prints again
   !>
   !> @param[inout] self   printer set up before
   !> @param[in]    diis   show the ndiis column
   !> @param[in]    gmres  show the GMRES block
   subroutine printer_set_columns(self, diis, gmres)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      !> Show the ndiis column and the GMRES block
      logical, intent(in) :: diis, gmres

      if ((diis .eqv. self%diis_column) .and. (gmres .eqv. self%gmres_block)) return
      self%diis_column = diis
      self%gmres_block = gmres
      call build_table(self)
      self%header_due = .true.
   end subroutine printer_set_columns

   !> Start a stage: the header row prints before the next row, dF restarts
   !>
   !> @param[inout] self  printer set up before
   subroutine printer_begin_stage(self)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self

      self%header_due = .true.
      self%energy_prev = huge(1.0_wp)
   end subroutine printer_begin_stage

   !> Print the header row and its separator when due
   !>
   !> @param[inout] self  printer set up before
   subroutine print_due_header(self)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self

      if (.not. self%header_due) return
      call self%plp%print_header()
      call self%plp%separator()
      self%header_due = .false.
   end subroutine print_due_header

   !> Build the table printer from the shown columns, on the current unit
   !>
   !> @param[inout] self  printer with its headers set
   subroutine build_table(self)
      !> Printer with its headers set
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      !> Whether each column is shown
      logical :: shown(NCOLS_TOTAL)

      shown = .true.
      shown(COL_LAMBDA) = self%lambda_column
      shown(COL_ENERGY:COL_ENERGY + 1) = self%energy_columns
      shown(NCOLS_COMMON + 1:NCOLS_COMMON + NCOLS_DIIS) = self%diis_column
      shown(NCOLS_COMMON + NCOLS_DIIS + 1) = self%gmres_block
      shown(NCOLS_COMMON + NCOLS_DIIS + 2:) = self%gmres_block .and. self%inner_rows
      self%plp = new_prettylistprinter(pack(widths, shown), pack(self%headers, shown), unit=self%plp%unit)
   end subroutine build_table

   !> Print the table title
   !>
   !> @param[inout] self   printer set up before
   !> @param[in]    title  title text
   subroutine printer_header(self, title)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      !> Title text
      character(*), intent(in) :: title
      call self%plp%header(title)
   end subroutine printer_header

   !> Print an indented setup detail or solver trace message
   !>
   !> - Empty text prints nothing
   !>
   !> @param[inout] self  printer set up before
   !> @param[in]    text  note text
   subroutine printer_note(self, text)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      !> Note text
      character(*), intent(in) :: text
      if (len_trim(text) == 0) return
      write (self%plp%unit, "(2x,a)") trim(text)
   end subroutine printer_note

   !> Print a separator line carrying a note, e.g. the solver in use
   !>
   !> @param[inout] self  printer set up before
   !> @param[in]    text  note text
   subroutine printer_sep_note(self, text)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      !> Note text
      character(*), intent(in) :: text
      call self%plp%separator()
      call self%note(text)
   end subroutine printer_sep_note

   !> Print the column header row
   !>
   !> @param[inout] self  printer set up before
   subroutine printer_print_header(self)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      call self%plp%print_header()
   end subroutine printer_print_header

   !> Print a table separator line
   !>
   !> @param[inout] self  printer set up before
   subroutine printer_separator(self)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      call self%plp%separator()
   end subroutine printer_separator

   !> Print an empty table line
   !>
   !> @param[inout] self  printer set up before
   subroutine printer_blank(self)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      call self%plp%blank()
   end subroutine printer_blank

   !> Print one labelled count below the table
   !>
   !> @param[inout] self   printer set up before
   !> @param[in]    label  count label
   !> @param[in]    n      count
   subroutine printer_stat_line(self, label, n)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      !> Count label
      character(*), intent(in) :: label
      !> Count
      integer, intent(in) :: n
      !> Label padded to the column of the dots
      character(len=26) :: padded

      padded = label
      write (self%plp%unit, "(2x,a,'...',i10)") padded, n
   end subroutine printer_stat_line

   !> Print one outer iteration row
   !>
   !> - Model diagnostics are optional; absent values leave their columns blank
   !> - Inner residual and tolerance columns are blank for outer rows
   !> - ndiis and gmres_total print only when their solver's columns are shown
   !> - lambda prints only when its column is shown
   !> - The energy is fixed-point below 1e7 in magnitude, else exponential
   !> - dF is the change from the previous row of the stage; a row without an
   !>   energy leaves both blank and restarts dF
   !>
   !> @param[inout] self         printer set up before
   !> @param[in]    iter         iteration number
   !> @param[in]    time         wall time column, seconds
   !> @param[in]    rmse         RMS residual (RMSR column)
   !> @param[in]    maxdx        maximum absolute change (max|dg| column)
   !> @param[in]    hmax         max(h) column; absent leaves it blank
   !> @param[in]    lambda       lambda column; absent leaves it blank
   !> @param[in]    energy       energy column; absent leaves it and dF blank
   !> @param[in]    alpha        step fraction; absent leaves it blank
   !> @param[in]    ndiis        stored DIIS vectors; absent or nonpositive leaves it blank
   !> @param[in]    gmres_total  inner-solver iteration count; absent or nonpositive leaves it blank
   subroutine printer_print_outer_row(self, iter, time, rmse, maxdx, hmax, &
                                  lambda, energy, alpha, ndiis, gmres_total)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      !> Iteration number
      integer, intent(in) :: iter
      !> Wall time in seconds, RMS residual and maximum absolute change
      real(wp), intent(in) :: time, rmse, maxdx
      !> Optional max(h), lambda, energy and step fraction columns
      real(wp), intent(in), optional :: hmax, lambda, energy, alpha
      !> Optional stored DIIS vector count and inner-solver iteration count
      integer, intent(in), optional :: ndiis, gmres_total

      call print_due_header(self)
      call self%plp%begin_row()
      call self%plp%add(iter)
      call self%plp%add(time, fmt="F7.2")
      if (self%lambda_column) then
         if (present(lambda)) then
            call self%plp%add(lambda, fmt="F9.3")
         else
            call self%plp%skip()
         end if
      end if
      if (self%energy_columns) then
         if (present(energy)) then
            if (abs(energy) < max_fixed_energy) then
               call self%plp%add(energy, fmt="F12.4")
            else
               call self%plp%add(energy, fmt="ES12.4")
            end if
            if (self%energy_prev < huge(1.0_wp)) then
               call self%plp%add(energy - self%energy_prev, fmt="ES10.2")
            else
               call self%plp%skip()
            end if
            self%energy_prev = energy
         else
            call self%plp%skip()
            call self%plp%skip()
            self%energy_prev = huge(1.0_wp)
         end if
      end if
      call self%plp%add(rmse, fmt="ES11.2")
      call self%plp%add(maxdx, fmt="ES11.2")
      if (present(hmax)) then
         call self%plp%add(hmax, fmt="F8.2")
      else
         call self%plp%skip()
      end if
      if (present(alpha)) then
         call self%plp%add(alpha, fmt="F7.4")
      else
         call self%plp%skip()
      end if
      if (self%diis_column) then
         if (present(ndiis)) then
            if (ndiis > 0) then
               call self%plp%add(ndiis)
            else
               call self%plp%skip()
            end if
         else
            call self%plp%skip()
         end if
      end if
      if (self%gmres_block) then
         ! Inner-solver iteration count
         if (present(gmres_total)) then
            if (gmres_total > 0) then
               call self%plp%add(gmres_total)
            else
               call self%plp%skip()
            end if
         else
            call self%plp%skip()
         end if
         if (self%inner_rows) then
            call self%plp%skip() ! krylov
            call self%plp%skip() ! res
            call self%plp%skip() ! tol
         end if
      end if
      call self%plp%end_row()
   end subroutine printer_print_outer_row

   !> Print one inner Krylov iteration row
   !>
   !> - Outer iteration columns are blank; inner counts, residual and tolerance are filled
   !> - Nothing prints while the GMRES block is hidden or inner rows are off
   !>
   !> @param[inout] self        printer set up before
   !> @param[in]    gmres_iter  GMRES iteration count
   !> @param[in]    krylov      Krylov dimension count
   !> @param[in]    rho         residual column value
   !> @param[in]    tol_val     tolerance column value
   subroutine printer_print_inner_gmres_row(self, gmres_iter, krylov, rho, tol_val)
      !> Printer set up before
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      !> GMRES and Krylov counts
      integer, intent(in) :: gmres_iter, krylov
      !> Residual and tolerance column values
      real(wp), intent(in) :: rho, tol_val

      !> Blank column counter
      integer :: i

      if (.not. (self%gmres_block .and. self%inner_rows)) return
      call print_due_header(self)
      call self%plp%begin_row()
      do i = 1, NCOLS_COMMON - merge(0, 1, self%lambda_column) - merge(0, 2, self%energy_columns) &
         & + merge(NCOLS_DIIS, 0, self%diis_column)
         call self%plp%skip()
      end do
      call self%plp%add(gmres_iter)
      call self%plp%add(krylov)
      call self%plp%add(rho, fmt="ES11.2")
      call self%plp%add(tol_val, fmt="ES11.2")
      call self%plp%end_row()
   end subroutine printer_print_inner_gmres_row

end module moist_model_moz_solver_core_printer
