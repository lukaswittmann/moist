!> Shared iteration table for MOZ solvers
!> Drivers supply outer diagnostics; JFNK writes the inner GMRES rows
module moist_model_moz_solver_printer
   use mctc_env, only: wp
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_printer_type

   integer, parameter :: NCOLS_OUTER = 10  !< columns before GMRES block
   integer, parameter :: NCOLS_GMRES = 4   !< gmres, krylov, res, tol
   integer, parameter :: NCOLS_TOTAL = NCOLS_OUTER + NCOLS_GMRES

   type :: moist_model_moz_solver_printer_type
      type(prettylistprinter) :: plp
   contains
      procedure :: setup => printer_setup
      procedure :: header => printer_header
      procedure :: note => printer_note
      procedure :: sep_note => printer_sep_note
      procedure :: print_header => printer_print_header
      procedure :: separator => printer_separator
      procedure :: blank => printer_blank
      procedure :: print_outer_row => printer_print_outer_row
      procedure :: print_inner_gmres_row => printer_print_inner_gmres_row
   end type moist_model_moz_solver_printer_type

contains

   !> Initialise the printer, required before any print method
   !> Optional energy_label/denergy_label override the default "F"/"dF" headers
   subroutine printer_setup(self, unit, energy_label, denergy_label)
      class(moist_model_moz_solver_printer_type), intent(out) :: self
      integer, intent(in), optional :: unit
      character(*), intent(in), optional :: energy_label, denergy_label

      integer, parameter :: widths(NCOLS_TOTAL) = &
                            [6, 7, 9, 12, 10, 11, 11, 8, 7, 6, 6, 7, 11, 11]
      character(8) :: hdrs(NCOLS_TOTAL)

      hdrs = [character(8) :: "Iter", "t (s)", "lambda", "F", "dF", "RMSR", &
              "max|dg|", "max(h)", "alpha", "ndiis", "gmres", "krylov", "res", "tol"]
      if (present(energy_label)) hdrs(4) = energy_label
      if (present(denergy_label)) hdrs(5) = denergy_label

      self%plp = new_prettylistprinter(widths, hdrs, unit=unit)
   end subroutine printer_setup

   subroutine printer_header(self, title)
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      character(*), intent(in) :: title
      call self%plp%header(title)
   end subroutine printer_header

   !> Print an indented setup detail or solver trace message
   !>
   !> @param[inout] self Printer set up before
   !> @param[in]    text note text
   subroutine printer_note(self, text)
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      character(*), intent(in) :: text
      write (self%plp%unit, "(2x,a)") trim(text)
   end subroutine printer_note

   !> Print a separator line carrying a note, e.g. the solver in use
   !>
   !> @param[inout] self Printer set up before
   !> @param[in]    text note text
   subroutine printer_sep_note(self, text)
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      character(*), intent(in) :: text
      call self%plp%separator()
      call self%note(text)
   end subroutine printer_sep_note

   subroutine printer_print_header(self)
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      call self%plp%print_header()
   end subroutine printer_print_header

   subroutine printer_separator(self)
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      call self%plp%separator()
   end subroutine printer_separator

   subroutine printer_blank(self)
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      call self%plp%blank()
   end subroutine printer_blank

   !> Print one outer iteration row
   !> Model diagnostics are optional; absent values leave their columns blank
   !> Inner residual and tolerance columns are blank for outer rows
   subroutine printer_print_outer_row(self, iter, time, rmse, maxdx, hmax, &
                                  lambda, free_energy, delta_f, alpha, ndiis, gmres_total)
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      integer, intent(in) :: iter
      real(wp), intent(in) :: time, rmse, maxdx
      real(wp), intent(in), optional :: hmax, lambda, free_energy, delta_f, alpha
      integer, intent(in), optional :: ndiis, gmres_total

      call self%plp%begin_row()
      call self%plp%add(iter)
      call self%plp%add(time, fmt="F7.2")
      if (present(lambda)) then
         call self%plp%add(lambda, fmt="F9.3")
      else
         call self%plp%skip()
      end if
      if (present(free_energy)) then
         call self%plp%add(free_energy, fmt="ES12.4")
      else
         call self%plp%skip()
      end if
      if (present(delta_f)) then
         call self%plp%add(delta_f, fmt="ES10.2")
      else
         call self%plp%skip()
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
      if (present(ndiis)) then
         if (ndiis > 0) then
            call self%plp%add(ndiis)
         else
            call self%plp%skip()
         end if
      else
         call self%plp%skip()
      end if
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
      call self%plp%skip() ! krylov
      call self%plp%skip() ! res
      call self%plp%skip() ! tol
      call self%plp%end_row()
   end subroutine printer_print_outer_row

   !> Print one inner Krylov iteration row
   !> Outer iteration columns are blank; inner counts, residual and tolerance are filled
   subroutine printer_print_inner_gmres_row(self, gmres_iter, krylov, rho, tol_val)
      class(moist_model_moz_solver_printer_type), intent(inout) :: self
      integer, intent(in) :: gmres_iter, krylov
      real(wp), intent(in) :: rho, tol_val

      integer :: i

      call self%plp%begin_row()
      do i = 1, NCOLS_OUTER
         call self%plp%skip()
      end do
      call self%plp%add(gmres_iter)
      call self%plp%add(krylov)
      call self%plp%add(rho, fmt="ES11.2")
      call self%plp%add(tol_val, fmt="ES11.2")
      call self%plp%end_row()
   end subroutine printer_print_inner_gmres_row

end module moist_model_moz_solver_printer
