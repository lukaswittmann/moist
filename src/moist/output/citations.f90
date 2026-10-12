
!> Citation registry for moist and its sub-models
!>
!> Each citation is stored as a typed entry with category, label,
!> authors, title, journal line, and DOI. The registry can be
!> printed in full or filtered by category/label
module moist_output_citations
   use moist_output_format, only: print_wrapped
   implicit none(type, external)
   private

   public :: citation_entry_type, moist_citations, num_citations, citation_label_len
   public :: print_citations, print_citations_by_category, print_citations_by_label
   public :: citation_moist, citation_drop, citation_isodensity, citation_iswig
   public :: citation_iswig_adaptive, citation_numsa, citation_cfc, citation_smd
   public :: citation_slsqp, citation_lbfgsb, citation_fmin, citation_lsqr
   public :: citation_lsmr, citation_lusol, citation_pv, citation_gostshyp
   public :: citation_cpcm, citation_cosmo, citation_iefpcm, citation_ssvpe
   public :: citation_dpcm, citation_dpcm_original, citation_lpb, citation_slic, citation_nonlocal

   !> Registry labels; a part records the ones it used through `ctx%cite`
   character(len=*), parameter :: citation_moist = "moist"
   character(len=*), parameter :: citation_drop = "DROP/SvdW"
   character(len=*), parameter :: citation_isodensity = "DROP/Isodensity"
   character(len=*), parameter :: citation_iswig = "Improved Switching Gaussian Approach (iSwiG)"
   character(len=*), parameter :: citation_iswig_adaptive = "Improved Switching Gaussian Approach (iSwiG) " &
      & // "with dynamic radii"
   character(len=*), parameter :: citation_numsa = "Numerical Surface Area (numSA)"
   character(len=*), parameter :: citation_cfc = "COSMO Fine Cavity (CFC)"
   character(len=*), parameter :: citation_smd = "SMD"
   character(len=*), parameter :: citation_slsqp = "SLSQP"
   character(len=*), parameter :: citation_lbfgsb = "L-BFGS-B"
   character(len=*), parameter :: citation_fmin = "fmin"
   character(len=*), parameter :: citation_lsqr = "LSQR"
   character(len=*), parameter :: citation_lsmr = "LSMR"
   character(len=*), parameter :: citation_lusol = "LUSOL"
   character(len=*), parameter :: citation_pv = "Pressure-volume term (PV)"
   character(len=*), parameter :: citation_gostshyp = "GOSTSHYP"
   character(len=*), parameter :: citation_cpcm = "CPCM"
   character(len=*), parameter :: citation_cosmo = "COSMO"
   character(len=*), parameter :: citation_iefpcm = "IEF-PCM"
   character(len=*), parameter :: citation_ssvpe = "SS(V)PE"
   character(len=*), parameter :: citation_dpcm = "DPCM"
   character(len=*), parameter :: citation_dpcm_original = "DPCM original"
   character(len=*), parameter :: citation_lpb = "LPB"
   character(len=*), parameter :: citation_slic = "SLIC"
   character(len=*), parameter :: citation_nonlocal = "Nonlocal electrostatics"

   !> Length that holds every label
   integer, parameter :: citation_label_len = 80

   !> A single literature reference
   type :: citation_entry_type
      !> Broad grouping: "General", "Cavities", "Models", or "Solvers"
      character(len=:), allocatable :: category
      !> Short identifier, e.g. "iSwiG", "GEMS", "SMD"
      character(len=:), allocatable :: label
      !> Author list (single line)
      character(len=:), allocatable :: authors
      !> Full title (may be long)
      character(len=:), allocatable :: title
      !> Journal, year, volume (formatted for display)
      character(len=:), allocatable :: journal
      !> DOI URL
      character(len=:), allocatable :: doi
   end type citation_entry_type

   !> Number of entries in the registry
   integer, parameter :: num_citations = 31

   !> The global citation registry, populated in init_citations()
   type(citation_entry_type), target :: moist_citations(num_citations)

   !> Guard against repeated initialisation
   logical :: initialised = .false.

   !> Maximum line width for word-wrapped output (excluding indent)
   integer, parameter :: wrap_width = 57

contains

   !> Populate the registry once (called lazily on first access)
   !>
   !> - contexts on different threads may print concurrently, so the first
   !>   fill is serialized; the critical construct also publishes the entries
   subroutine init_citations()
      !$omp critical(moist_citations_init)
      if (.not. initialised) call fill_citations()
      !$omp end critical(moist_citations_init)
   end subroutine init_citations

   !> Fill the registry entries
   subroutine fill_citations()

      ! Concept DOI: resolves to the latest release
      moist_citations(1) = citation_entry_type( &
                           category="General", &
                           label=citation_moist, &
                           authors="Wittmann, L.", &
                           title="MOIST: The Modular and Open-source Implicit Solvation Toolkit.", &
                           journal="Computer software, Zenodo 2026", &
                           doi="https://doi.org/10.5281/zenodo.23473416")

      moist_citations(2) = citation_entry_type( &
                           category="Cavities", &
                           label=citation_drop, &
                           authors="Wittmann, L., Pausch, A.", &
                           title="A Smooth and Fully Differentiable Molecular Cavity Based on "// &
                           "Discretization via Reference-Onto-Surface Projection", &
                           journal="ChemRxiv 2026", &
                           doi="https://doi.org/10.26434/chemrxiv.15003893/v2")

      moist_citations(3) = citation_entry_type( &
                           category="Cavities", &
                           label=citation_isodensity, &
                           authors="Wittmann, L.", &
                           title="From Density to Boundary and Back: A Fully Differentiable, "// &
                           "Self-Consistent Isodensity Cavity.", &
                           journal="ChemRxiv, 06 August 2026", &
                           doi="https://doi.org/10.26434/chemrxiv.15007095/v1")

      moist_citations(4) = citation_entry_type( &
                           category="Cavities", &
                           label=citation_iswig, &
                           authors="Lange, A.W., Herbert, J.M.", &
                           title="A smooth, nonsingular, and faithful discretization scheme "// &
                           "for PCM: the switching/Gaussian approach.", &
                           journal="J. Chem. Phys. 2010, 133", &
                           doi="https://doi.org/10.1063/1.3511297")

      moist_citations(5) = citation_entry_type( &
                           category="Cavities", &
                           label=citation_iswig_adaptive, &
                           authors="Wittmann, L., Garcia-Rates, M., Riplinger, C.", &
                           title="Analytical first derivatives of the SCF energy for the "// &
                           "conductor-like PCM with non-static radii.", &
                           journal="J. Comput. Chem. 2025, 46", &
                           doi="https://doi.org/10.1002/jcc.70099")

      moist_citations(6) = citation_entry_type( &
                           category="Cavities", &
                           label=citation_numsa, &
                           authors="Im, W., Lee, M.S., Brooks, C.L.", &
                           title="Generalized born model with a simple smoothing function.", &
                           journal="J. Comput. Chem. 2003, 24", &
                           doi="https://doi.org/10.1002/jcc.10321")

      moist_citations(7) = citation_entry_type( &
                           category="Cavities", &
                           label=citation_cfc, &
                           authors="Klamt, A., Diedenhofen, M.", &
                           title="A refined cavity construction algorithm for the conductor-like screening model.", &
                           journal="J. Comput. Chem. 2018, 39, 1648-1655", &
                           doi="https://doi.org/10.1002/jcc.25342")

      moist_citations(8) = citation_entry_type( &
                           category="Models", &
                           label=citation_smd, &
                           authors="Marenich, A.V., Cramer, C.J., Truhlar, D.G.", &
                           title="Universal Solvation Model Based on Solute Electron Density "// &
                           "and on a Continuum Model of the Solvent Defined by the "// &
                           "Bulk Dielectric Constant and Atomic Surface Tensions.", &
                           journal="J. Phys. Chem. B 2009, 113, 18", &
                           doi="https://doi.org/10.1021/jp810292n")

      moist_citations(9) = citation_entry_type( &
                           category="Solvers", &
                           label=citation_slsqp, &
                           authors="Kraft, D.", &
                           title="A software package for sequential quadratic programming.", &
                           journal="Tech. Rep. DFVLR-FB 88-28, DLR German Aerospace Center, 1988", &
                           doi="")

      moist_citations(10) = citation_entry_type( &
                            category="Solvers", &
                            label=citation_slsqp, &
                            authors="Kraft, D.", &
                            title="Algorithm 733: TOMP - Fortran modules for "// &
                            "optimal control calculations.", &
                            journal="ACM Trans. Math. Softw. 1994, 20, 262-281", &
                            doi="https://doi.org/10.1145/192115.192124")

      moist_citations(11) = citation_entry_type( &
                            category="Solvers", &
                            label=citation_lbfgsb, &
                            authors="Byrd, R.H., Lu, P., Nocedal, J., Zhu, C.", &
                            title="A limited memory algorithm for bound constrained optimization.", &
                            journal="SIAM J. Sci. Comput. 1995, 16, 1190-1208", &
                            doi="https://doi.org/10.1137/0916069")

      moist_citations(12) = citation_entry_type( &
                            category="Solvers", &
                            label=citation_lbfgsb, &
                            authors="Zhu, C., Byrd, R.H., Lu, P., Nocedal, J.", &
                            title="Algorithm 778: L-BFGS-B: Fortran subroutines for "// &
                            "large-scale bound-constrained optimization.", &
                            journal="ACM Trans. Math. Softw. 1997, 23, 550-560", &
                            doi="https://doi.org/10.1145/279232.279236")

      moist_citations(13) = citation_entry_type( &
                            category="Solvers", &
                            label=citation_lbfgsb, &
                            authors="Morales, J.L., Nocedal, J.", &
                            title="Remark on Algorithm 778: L-BFGS-B: Fortran subroutines "// &
                            "for large-scale bound constrained optimization.", &
                            journal="ACM Trans. Math. Softw. 2011, 38, 7", &
                            doi="https://doi.org/10.1145/2049662.2049669")

      moist_citations(14) = citation_entry_type( &
                            category="Solvers", &
                            label=citation_fmin, &
                            authors="Brent, R.P.", &
                            title="Algorithms for Minimization Without Derivatives.", &
                            journal="Prentice-Hall, Englewood Cliffs, NJ, 1973", &
                            doi="https://maths-people.anu.edu.au/~brent/pub/pub011.html")

      moist_citations(15) = citation_entry_type( &
                            category="Solvers", &
                            label=citation_lsqr, &
                            authors="Paige, C.C., Saunders, M.A.", &
                            title="LSQR: An algorithm for sparse linear equations and "// &
                            "sparse least squares.", &
                            journal="ACM Trans. Math. Softw. 1982, 8, 43-71", &
                            doi="https://doi.org/10.1145/355984.355989")

      moist_citations(16) = citation_entry_type( &
                            category="Solvers", &
                            label=citation_lsmr, &
                            authors="Fong, D.C.-L., Saunders, M.A.", &
                            title="LSMR: An iterative algorithm for sparse "// &
                            "least-squares problems.", &
                            journal="SIAM J. Sci. Comput. 2011, 33, 2950-2971", &
                            doi="https://doi.org/10.1137/10079687X")

      moist_citations(17) = citation_entry_type( &
                            category="Solvers", &
                            label=citation_lusol, &
                            authors="Gill, P.E., Murray, W., Saunders, M.A., Wright, M.H.", &
                            title="Maintaining LU factors of a general sparse matrix.", &
                            journal="Linear Algebra Appl. 1987, 88-89, 239-270", &
                            doi="https://doi.org/10.1016/0024-3795(87)90112-1")

      moist_citations(18) = citation_entry_type( &
                            category="Models", &
                            label=citation_pv, &
                            authors="Spooner, J., Yanciw, B., Wiebe, B., Weinberg, N.", &
                            title="Reaction Profiles and Energy Surfaces of Compressed Species.", &
                            journal="J. Phys. Chem. A 2014, 118, 765-777", &
                            doi="https://doi.org/10.1021/jp410496t")

      moist_citations(19) = citation_entry_type( &
                            category="Models", &
                            label=citation_pv, &
                            authors="Zeller, F., Pracht, P., Neudecker, T.", &
                            title="Using Conformational Sampling to Model Spectral and Structural "// &
                            "Changes of Molecules at Elevated Pressures.", &
                            journal="J. Phys. Chem. A 2025, 129, 2108-2116", &
                            doi="https://doi.org/10.1021/acs.jpca.4c08065")

      moist_citations(20) = citation_entry_type( &
                            category="Models", &
                            label=citation_gostshyp, &
                            authors="Scheurer, M., Dreuw, A., Epifanovsky, E., Head-Gordon, M., Stauch, T.", &
                            title="Modeling Molecules under Pressure with Gaussian Potentials.", &
                            journal="J. Chem. Theory Comput. 2021, 17, 583-597", &
                            doi="https://doi.org/10.1021/acs.jctc.0c01212")

      moist_citations(21) = citation_entry_type( &
                            category="Models", &
                            label=citation_gostshyp, &
                            authors="Pausch, A., Zeller, F., Neudecker, T.", &
                            title="An Integral-Direct GOSTSHYP Algorithm for the Computation of "// &
                            "High Pressure Effects on Molecular and Electronic Structure.", &
                            journal="J. Chem. Theory Comput. 2025, 21, 747-761", &
                            doi="https://doi.org/10.1021/acs.jctc.4c01502")

      moist_citations(22) = citation_entry_type( &
                            category="Models", &
                            label=citation_cpcm, &
                            authors="Barone, V., Cossi, M.", &
                            title="Quantum Calculation of Molecular Energies and Energy Gradients "// &
                            "in Solution by a Conductor Solvent Model.", &
                            journal="J. Phys. Chem. A 1998, 102, 1995-2001", &
                            doi="https://doi.org/10.1021/jp9716997")

      moist_citations(23) = citation_entry_type( &
                            category="Models", &
                            label=citation_cpcm, &
                            authors="Cossi, M., Rega, N., Scalmani, G., Barone, V.", &
                            title="Energies, structures, and electronic properties of molecules "// &
                            "in solution with the C-PCM solvation model.", &
                            journal="J. Comput. Chem. 2003, 24, 669-681", &
                            doi="https://doi.org/10.1002/jcc.10189")

      moist_citations(24) = citation_entry_type( &
                            category="Models", &
                            label=citation_cosmo, &
                            authors="Klamt, A., Schuurmann, G.", &
                            title="COSMO: a new approach to dielectric screening in solvents "// &
                            "with explicit expressions for the screening energy and its gradient.", &
                            journal="J. Chem. Soc., Perkin Trans. 2 1993, 799-805", &
                            doi="https://doi.org/10.1039/P29930000799")

      moist_citations(25) = citation_entry_type( &
                            category="Models", &
                            label=citation_iefpcm, &
                            authors="Cances, E., Mennucci, B., Tomasi, J.", &
                            title="A new integral equation formalism for the polarizable continuum model: "// &
                            "Theoretical background and applications to isotropic and anisotropic dielectrics.", &
                            journal="J. Chem. Phys. 1997, 107, 3032-3041", &
                            doi="https://doi.org/10.1063/1.474659")

      moist_citations(26) = citation_entry_type( &
                            category="Models", &
                            label=citation_ssvpe, &
                            authors="Chipman, D.M.", &
                            title="Reaction field treatment of charge penetration.", &
                            journal="J. Chem. Phys. 2000, 112, 5558-5565", &
                            doi="https://doi.org/10.1063/1.481133")

      moist_citations(27) = citation_entry_type( &
                            category="Models", &
                            label=citation_dpcm, &
                            authors="Miertus, S., Scrocco, E., Tomasi, J.", &
                            title="Electrostatic interaction of a solute with a continuum. "// &
                            "A direct utilization of ab initio molecular potentials for the prevision of solvent effects.", &
                            journal="Chem. Phys. 1981, 55, 117-129", &
                            doi="https://doi.org/10.1016/0301-0104(81)85090-2")

      moist_citations(28) = moist_citations(27)
      moist_citations(28)%label = citation_dpcm_original

      moist_citations(29) = citation_entry_type( &
                            category="Models", &
                            label=citation_lpb, &
                            authors="Juffer, A.H., Botta, E.F.F., van Keulen, B.A.M., "// &
                            "van der Ploeg, A., Berendsen, H.J.C.", &
                            title="The electric potential of a macromolecule in a solvent: "// &
                            "A fundamental approach.", &
                            journal="J. Comput. Phys. 1991, 97, 144-171", &
                            doi="https://doi.org/10.1016/0021-9991(91)90043-K")

      moist_citations(30) = citation_entry_type( &
                            category="Models", &
                            label=citation_slic, &
                            authors="Molavi Tabrizi, A., Goossens, S., Mehdizadeh Rahimi, A., "// &
                            "Knepley, M., Bardhan, J.P.", &
                            title="Predicting solvation free energies and thermodynamics in polar "// &
                            "solvents and mixtures using a solvation-layer interface condition.", &
                            journal="J. Chem. Phys. 2017, 146, 094103", &
                            doi="https://doi.org/10.1063/1.4977037")

      moist_citations(31) = citation_entry_type( &
                            category="Models", &
                            label=citation_nonlocal, &
                            authors="Hildebrandt, A., Blossey, R., Rjasanow, S., "// &
                            "Kohlbacher, O., Lenhof, H.-P.", &
                            title="Novel Formulation of Nonlocal Electrostatics.", &
                            journal="Phys. Rev. Lett. 2004, 93, 108104", &
                            doi="https://doi.org/10.1103/PhysRevLett.93.108104")

      initialised = .true.
   end subroutine fill_citations

   !> Print category header
   subroutine print_category_header(unit, category)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Category name to print
      character(len=*), intent(in) :: category

      write (unit, '(a,a,":")') "", category
      write (unit, "(a)") ""

   end subroutine print_category_header

   !> Print a single citation entry with consistent formatting
   !>
   !> - long lines are word-wrapped at wrap_width characters
   subroutine print_entry(unit, entry, indent)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Citation to print
      type(citation_entry_type), intent(in) :: entry
      !> Prefix of every line
      character(len=*), intent(in) :: indent

      call print_wrapped(unit, entry%authors, indent, wrap_width)
      call print_wrapped(unit, entry%title, indent, wrap_width)
      if (len_trim(entry%doi) > 0) then
         call print_wrapped(unit, entry%journal//", "//entry%doi, indent, wrap_width)
      else
         call print_wrapped(unit, entry%journal, indent, wrap_width)
      end if
      write (unit, "(a)") ""
   end subroutine print_entry

   !> Print the General entries (the moist package), which always print, then
   !> the other citations grouped by category: all of them, or with `labels`
   !> only those labels
   subroutine print_citations(unit, labels)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Labels to print besides the General entries; absent prints all
      character(len=*), intent(in), optional :: labels(:)

      integer :: i

      call init_citations()

      write (unit, "(a)") "Please include the appropriate citations when using our work:"
      write (unit, "(a)") ""

      ! The package itself, outside any category
      do i = 1, num_citations
         if (moist_citations(i)%category == "General") call print_entry(unit, moist_citations(i), "  ")
      end do

      call print_citations_by_category(unit, "Cavities", labels)
      call print_citations_by_category(unit, "Models", labels)
      call print_citations_by_category(unit, "Solvers", labels)

   end subroutine print_citations

   !> Print all citations matching a given category
   !>
   !> - absent `labels`, entries are grouped under their label as sub-header
   !> - with `labels`, only entries carrying one of them, without sub-headers;
   !>   the category header is skipped when none match
   subroutine print_citations_by_category(unit, category, labels)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Category to filter by (e.g. "General", "Cavities", "Models")
      character(len=*), intent(in) :: category
      !> Labels to keep; absent keeps all
      character(len=*), intent(in), optional :: labels(:)

      integer :: i
      logical :: header_done
      character(len=:), allocatable :: last_label

      call init_citations()

      header_done = .false.
      last_label = ""
      do i = 1, num_citations
         if (moist_citations(i)%category /= category) cycle
         if (present(labels)) then
            if (.not. any(labels == moist_citations(i)%label)) cycle
         end if
         if (.not. header_done) then
            call print_category_header(unit, category)
            header_done = .true.
         end if
         if (present(labels)) then
            call print_entry(unit, moist_citations(i), "  ")
            cycle
         end if
         ! Print label sub-header when it changes
         if (moist_citations(i)%label /= last_label) then
            write (unit, '(2x,a,a,":")') "", moist_citations(i)%label
            write (unit, "(a)") ""
            last_label = moist_citations(i)%label
         end if
         call print_entry(unit, moist_citations(i), "    ")
      end do

   end subroutine print_citations_by_category

   !> Print all citations matching a given label (e.g. "iSwiG", "GEMS")
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
            write (unit, "(a)") ""
            found = .true.
         end if
         call print_entry(unit, moist_citations(i), "  ")
      end do

   end subroutine print_citations_by_label

end module moist_output_citations
