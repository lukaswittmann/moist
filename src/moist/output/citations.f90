
!> Citation registry for moist and its sub-models.
!> Each citation is stored as a typed entry with category, label,
!> authors, title, journal line, and DOI. The registry can be
!> printed in full or filtered by category/label.
module moist_output_citations
   use moist_output_format, only: print_wrapped
   implicit none
   private

   public :: citation_entry, moist_citations, num_citations
   public :: print_citations, print_citations_by_category, print_citations_by_label

   !> A single literature reference.
   type :: citation_entry
      !> Broad grouping: "General", "Cavities", "Models", or "Solvers"
      character(len=:), allocatable :: category
      !> Short identifier, e.g. "iSwiG", "GEMS", "ALPB"
      character(len=:), allocatable :: label
      !> Author list (single line)
      character(len=:), allocatable :: authors
      !> Full title (may be long)
      character(len=:), allocatable :: title
      !> Journal, year, volume (formatted for display)
      character(len=:), allocatable :: journal
      !> DOI URL
      character(len=:), allocatable :: doi
   end type citation_entry

   !> Number of entries in the registry
   integer, parameter :: num_citations = 17

   !> The global citation registry, populated in init_citations().
   type(citation_entry), target :: moist_citations(num_citations)

   !> Guard against repeated initialisation
   logical :: initialised = .false.

   !> Maximum line width for word-wrapped output (excluding indent).
   integer, parameter :: wrap_width = 57

contains

   !> Populate the registry (called lazily on first access).
   subroutine init_citations()
      if (initialised) return

      ! No dedicated toolkit paper yet; cite the DROP cavity work for now.
      moist_citations(1) = citation_entry( &
                           category="General", &
                           label="moist", &
                           authors="Wittmann, L., Pausch, A.", &
                           title="A Smooth and Fully Differentiable Molecular Cavity Based on "// &
                           " Discretization via Reference-Onto-Surface Projection", &
                           journal="ChemRxiv 2026", &
                           doi="https://doi.org/10.26434/chemrxiv.15003893/v2")

      moist_citations(2) = citation_entry( &
                           category="Cavities", &
                           label="DROP/SvdW", &
                           authors="Wittmann, L., Pausch, A.", &
                           title="A Smooth and Fully Differentiable Molecular Cavity Based on "// &
                           " Discretization via Reference-Onto-Surface Projection", &
                           journal="ChemRxiv 2026", &
                           doi="https://doi.org/10.26434/chemrxiv.15003893/v2")

      moist_citations(3) = citation_entry( &
                           category="Cavities", &
                           label="Improved Switching Gaussian Approach (iSwiG)", &
                           authors="Lange, A.W., Herbert, J.M.", &
                           title="A smooth, nonsingular, and faithful discretization scheme "// &
                           "for PCM: the switching/Gaussian approach.", &
                           journal="J. Chem. Phys. 2010, 133", &
                           doi="https://doi.org/10.1063/1.3511297")

      moist_citations(4) = citation_entry( &
                           category="Cavities", &
                           label="Improved Switching Gaussian Approach (iSwiG) with adaptive radii", &
                           authors="Wittmann, L., Garcia-Rates, M., Riplinger, C.", &
                           title="Analytical first derivatives of the SCF energy for the "// &
                           "conductor-like PCM with non-static radii.", &
                           journal="J. Comput. Chem. 2025, 46", &
                           doi="https://doi.org/10.1002/jcc.70099")

      moist_citations(5) = citation_entry( &
                           category="Cavities", &
                           label="Numerical Surface Area (numSA)", &
                           authors="Im, W., Lee, M.S., Brooks, C.L.", &
                           title="Generalized born model with a simple smoothing function.", &
                           journal="J. Comput. Chem. 2003, 24", &
                           doi="https://doi.org/10.1002/jcc.10321")

      moist_citations(6) = citation_entry( &
                           category="Cavities", &
                           label="COSMO Fine Cavity (CFC)", &
                           authors="Klamt, A., Diedenhofen, M.", &
                           title="A refined cavity construction algorithm for the conductor-like screening model.", &
                           journal="J. Comput. Chem. 2018, 39, 1648-1655", &
                           doi="https://doi.org/10.1002/jcc.25342")

      moist_citations(7) = citation_entry( &
                           category="Models", &
                           label="ALPB", &
                           authors="Ehlert, S., Stahn, M., Spicher, S., Grimme, S.", &
                           title="Robust and Efficient Implicit Solvation Model for "// &
                           "Fast Semiempirical Methods.", &
                           journal="J. Chem. Theory Comput. 2021, 17", &
                           doi="https://doi.org/10.1021/acs.jctc.1c00471")

      moist_citations(8) = citation_entry( &
                           category="Models", &
                           label="SMD", &
                           authors="Marenich, A.V., Cramer, C.J., Truhlar, D.G.", &
                           title="Universal Solvation Model Based on Solute Electron Density "// &
                           "and on a Continuum Model of the Solvent Defined by the "// &
                           "Bulk Dielectric Constant and Atomic Surface Tensions.", &
                           journal="J. Phys. Chem. B 2009, 113, 18", &
                           doi="https://doi.org/10.1021/jp810292n")

      moist_citations(9) = citation_entry( &
                           category="Solvers", &
                           label="SLSQP", &
                           authors="Kraft, D.", &
                           title="A software package for sequential quadratic programming.", &
                           journal="Tech. Rep. DFVLR-FB 88-28, DLR German Aerospace Center, 1988", &
                           doi="")

      moist_citations(10) = citation_entry( &
                            category="Solvers", &
                            label="SLSQP", &
                            authors="Kraft, D.", &
                            title="Algorithm 733: TOMP - Fortran modules for "// &
                            "optimal control calculations.", &
                            journal="ACM Trans. Math. Softw. 1994, 20, 262-281", &
                            doi="https://doi.org/10.1145/192115.192124")

      moist_citations(11) = citation_entry( &
                            category="Solvers", &
                            label="L-BFGS-B", &
                            authors="Byrd, R.H., Lu, P., Nocedal, J., Zhu, C.", &
                            title="A limited memory algorithm for bound constrained optimization.", &
                            journal="SIAM J. Sci. Comput. 1995, 16, 1190-1208", &
                            doi="https://doi.org/10.1137/0916069")

      moist_citations(12) = citation_entry( &
                            category="Solvers", &
                            label="L-BFGS-B", &
                            authors="Zhu, C., Byrd, R.H., Lu, P., Nocedal, J.", &
                            title="Algorithm 778: L-BFGS-B: Fortran subroutines for "// &
                            "large-scale bound-constrained optimization.", &
                            journal="ACM Trans. Math. Softw. 1997, 23, 550-560", &
                            doi="https://doi.org/10.1145/279232.279236")

      moist_citations(13) = citation_entry( &
                            category="Solvers", &
                            label="L-BFGS-B", &
                            authors="Morales, J.L., Nocedal, J.", &
                            title="Remark on Algorithm 778: L-BFGS-B: Fortran subroutines "// &
                            "for large-scale bound constrained optimization.", &
                            journal="ACM Trans. Math. Softw. 2011, 38, 7", &
                            doi="https://doi.org/10.1145/2049662.2049669")

      moist_citations(14) = citation_entry( &
                            category="Solvers", &
                            label="fmin", &
                            authors="Brent, R.P.", &
                            title="Algorithms for Minimization Without Derivatives.", &
                            journal="Prentice-Hall, Englewood Cliffs, NJ, 1973", &
                            doi="https://maths-people.anu.edu.au/~brent/pub/pub011.html")

      moist_citations(15) = citation_entry( &
                            category="Solvers", &
                            label="LSQR", &
                            authors="Paige, C.C., Saunders, M.A.", &
                            title="LSQR: An algorithm for sparse linear equations and "// &
                            "sparse least squares.", &
                            journal="ACM Trans. Math. Softw. 1982, 8, 43-71", &
                            doi="https://doi.org/10.1145/355984.355989")

      moist_citations(16) = citation_entry( &
                            category="Solvers", &
                            label="LSMR", &
                            authors="Fong, D.C.-L., Saunders, M.A.", &
                            title="LSMR: An iterative algorithm for sparse "// &
                            "least-squares problems.", &
                            journal="SIAM J. Sci. Comput. 2011, 33, 2950-2971", &
                            doi="https://doi.org/10.1137/10079687X")

      moist_citations(17) = citation_entry( &
                            category="Solvers", &
                            label="LUSOL", &
                            authors="Gill, P.E., Murray, W., Saunders, M.A., Wright, M.H.", &
                            title="Maintaining LU factors of a general sparse matrix.", &
                            journal="Linear Algebra Appl. 1987, 88-89, 239-270", &
                            doi="https://doi.org/10.1016/0024-3795(87)90112-1")

      initialised = .true.
   end subroutine init_citations

   !> Print category header
   subroutine print_category_header(unit, category)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Category name to print
      character(len=*), intent(in) :: category

      write (unit, '(a,a,":")') "", category
      write (unit, '(a)') ""

   end subroutine print_category_header

   !> Print a single citation entry with consistent formatting.
   !> Long lines are word-wrapped at wrap_width characters.
   subroutine print_entry(unit, entry)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Citation to print
      type(citation_entry), intent(in) :: entry

      call print_wrapped(unit, entry%authors, "    ", wrap_width)
      call print_wrapped(unit, entry%title, "    ", wrap_width)
      if (len_trim(entry%doi) > 0) then
         call print_wrapped(unit, entry%journal//", "//entry%doi, "    ", wrap_width)
      else
         call print_wrapped(unit, entry%journal, "    ", wrap_width)
      end if
      write (unit, '(a)') ""
   end subroutine print_entry

   !> Print all citations, grouped by category.
   subroutine print_citations(unit)
      !> Fortran I/O unit
      integer, intent(in) :: unit

      call init_citations()

      write (unit, '(a)') "Please include the appropriate citations when using our work:"
      write (unit, '(a)') ""

      call print_citations_by_category(unit, "General")
      call print_citations_by_category(unit, "Cavities")
      call print_citations_by_category(unit, "Models")
      call print_citations_by_category(unit, "Solvers")

   end subroutine print_citations

   !> Print all citations matching a given category.
   !> Entries are grouped under their label as sub-header.
   subroutine print_citations_by_category(unit, category)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Category to filter by (e.g. "General", "Cavities", "Models")
      character(len=*), intent(in) :: category

      integer :: i
      character(len=:), allocatable :: last_label

      call init_citations()

      call print_category_header(unit, category)

      last_label = ""
      do i = 1, num_citations
         if (moist_citations(i)%category /= category) cycle
         ! Print label sub-header when it changes
         if (moist_citations(i)%label /= last_label) then
            write (unit, '(2x,a,a,":")') "", moist_citations(i)%label
            write (unit, '(a)') ""
            last_label = moist_citations(i)%label
         end if
         call print_entry(unit, moist_citations(i))
      end do

   end subroutine print_citations_by_category

   !> Print all citations matching a given label (e.g. "iSwiG", "GEMS").
   subroutine print_citations_by_label(unit, label)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Label to filter by
      character(len=*), intent(in) :: label

      integer :: i
      logical :: found

      call init_citations()

      found = .false.
      do i = 1, num_citations
         if (moist_citations(i)%label /= label) cycle
         if (.not. found) then
            write (unit, '(a,a,":")') "", label
            write (unit, '(a)') ""
            found = .true.
         end if
         call print_entry(unit, moist_citations(i))
      end do

   end subroutine print_citations_by_label

end module moist_output_citations
