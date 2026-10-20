
!> Citation registry for moist and its sub-models
!>
!> Each citation is stored as a typed entry with category, the labels of
!> the parts citing it, authors, title, journal line, and DOI. The
!> registry can be printed in full or filtered by label
module moist_output_citations
   use moist_output_format, only: print_wrapped
   implicit none(type, external)
   private

   public :: citation_entry_type, citation_registry_type, moist_citations
   public :: print_citations, print_citations_by_category
   public :: citation

   !> Registry labels; a part records the ones it used through `ctx%cite`
   type :: citation_label_type
      integer :: moist = 1
      integer :: drop = 2
      integer :: isodensity = 3
      integer :: iswig = 4
      integer :: iswig_adaptive = 5
      integer :: numsa = 6
      integer :: cfc = 7
      integer :: smd = 8
      integer :: slsqp = 9
      integer :: lbfgsb = 10
      integer :: fmin = 11
      integer :: lsqr = 12
      integer :: lsmr = 13
      integer :: lusol = 14
      integer :: gmres = 15
      integer :: pv = 16
      integer :: gostshyp = 17
      integer :: cpcm = 18
      integer :: cosmo = 19
      integer :: iefpcm = 20
      integer :: ssvpe = 21
      integer :: dpcm = 22
      integer :: dpcm_original = 23
      integer :: lpb = 24
      integer :: slic = 25
      integer :: nonlocal = 26
   end type citation_label_type

   !> Global instance of the label enum
   type(citation_label_type), parameter :: citation = citation_label_type()

   !> A single literature reference
   type :: citation_entry_type
      !> Broad grouping: "General", "Cavities & Discretization", "Models", or "Solvers"
      character(len=:), allocatable :: category
      !> Labels of the parts citing it
      integer, allocatable :: labels(:)
      !> Author list (single line)
      character(len=:), allocatable :: authors
      !> Full title (may be long)
      character(len=:), allocatable :: title
      !> Journal, year, volume (formatted for display)
      character(len=:), allocatable :: journal
      !> DOI URL
      character(len=:), allocatable :: doi
   end type citation_entry_type

   !> Growable list of citations, printed in insertion order
   type :: citation_registry_type
      !> Entries; only the first `n` are set
      type(citation_entry_type), allocatable :: entries(:)
      !> Number of entries
      integer :: n = 0
   contains
      !> Append one entry
      procedure :: add => registry_add
      !> Add a label to an already added entry, found by its DOI
      procedure :: alias => registry_alias
   end type citation_registry_type

   !> The global citation registry, populated in init_citations()
   type(citation_registry_type), target :: moist_citations

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
      call moist_citations%add( &
         category="General", &
         label=citation%moist, &
         authors="Wittmann, L.", &
         title="MOIST: The Modular and Open-source Implicit Solvation Toolkit.", &
         journal="Computer software, Zenodo 2026", &
         doi="https://doi.org/10.5281/zenodo.23473416")

      call moist_citations%add( &
         category="Cavities & Discretization", &
         label=citation%drop, &
         authors="Wittmann, L., Pausch, A.", &
         title="A Smooth and Fully Differentiable Molecular Cavity Based on "// &
         "Discretization via Reference-Onto-Surface Projection", &
         journal="ChemRxiv 2026", &
         doi="https://doi.org/10.26434/chemrxiv.15003893/v2")

      call moist_citations%add( &
         category="Cavities & Discretization", &
         label=citation%isodensity, &
         authors="Wittmann, L.", &
         title="From Density to Boundary and Back: A Fully Differentiable, "// &
         "Self-Consistent Isodensity Cavity.", &
         journal="ChemRxiv, 06 August 2026", &
         doi="https://doi.org/10.26434/chemrxiv.15007095/v1")

      call moist_citations%add( &
         category="Cavities & Discretization", &
         label=citation%iswig, &
         authors="Lange, A.W., Herbert, J.M.", &
         title="A smooth, nonsingular, and faithful discretization scheme "// &
         "for polarizable continuum models: The switching/Gaussian approach.", &
         journal="J. Chem. Phys. 2010, 133, 244111", &
         doi="https://doi.org/10.1063/1.3511297")

      call moist_citations%add( &
         category="Cavities & Discretization", &
         label=citation%iswig, &
         authors="York, D.M., Karplus, M.", &
         title="A Smooth Solvation Potential Based on the Conductor-Like Screening Model.", &
         journal="J. Phys. Chem. A 1999, 103, 11060-11079", &
         doi="https://doi.org/10.1021/jp992097l")

      call moist_citations%add( &
         category="Cavities & Discretization", &
         label=citation%iswig_adaptive, &
         authors="Wittmann, L., Garcia-Rates, M., Riplinger, C.", &
         title="Analytical First Derivatives of the SCF Energy for the "// &
         "Conductor-Like Polarizable Continuum Model With Non-Static Radii.", &
         journal="J. Comput. Chem. 2025, 46, e70099", &
         doi="https://doi.org/10.1002/jcc.70099")

      call moist_citations%add( &
         category="Cavities & Discretization", &
         label=citation%numsa, &
         authors="Im, W., Lee, M.S., Brooks, C.L.", &
         title="Generalized born model with a simple smoothing function.", &
         journal="J. Comput. Chem. 2003, 24, 1691-1702", &
         doi="https://doi.org/10.1002/jcc.10321")

      call moist_citations%add( &
         category="Cavities & Discretization", &
         label=citation%cfc, &
         authors="Klamt, A., Diedenhofen, M.", &
         title="A refined cavity construction algorithm for the conductor-like"// &
         " screening model.", &
         journal="J. Comput. Chem. 2018, 39, 1648-1655", &
         doi="https://doi.org/10.1002/jcc.25342")

      call moist_citations%add( &
         category="Models", &
         label=citation%smd, &
         authors="Marenich, A.V., Cramer, C.J., Truhlar, D.G.", &
         title="Universal Solvation Model Based on Solute Electron Density "// &
         "and on a Continuum Model of the Solvent Defined by the "// &
         "Bulk Dielectric Constant and Atomic Surface Tensions.", &
         journal="J. Phys. Chem. B 2009, 113, 6378-6396", &
         doi="https://doi.org/10.1021/jp810292n")

      call moist_citations%add( &
         category="Solvers", &
         label=citation%slsqp, &
         authors="Kraft, D.", &
         title="A software package for sequential quadratic programming.", &
         journal="Tech. Rep. DFVLR-FB 88-28, DLR German Aerospace Center, 1988", &
         doi="")

      call moist_citations%add( &
         category="Solvers", &
         label=citation%slsqp, &
         authors="Kraft, D.", &
         title="Algorithm 733: TOMP - Fortran modules for "// &
         "optimal control calculations.", &
         journal="ACM Trans. Math. Softw. 1994, 20, 262-281", &
         doi="https://doi.org/10.1145/192115.192124")

      call moist_citations%add( &
         category="Solvers", &
         label=citation%lbfgsb, &
         authors="Byrd, R.H., Lu, P., Nocedal, J., Zhu, C.", &
         title="A limited memory algorithm for bound constrained optimization.", &
         journal="SIAM J. Sci. Comput. 1995, 16, 1190-1208", &
         doi="https://doi.org/10.1137/0916069")

      call moist_citations%add( &
         category="Solvers", &
         label=citation%lbfgsb, &
         authors="Zhu, C., Byrd, R.H., Lu, P., Nocedal, J.", &
         title="Algorithm 778: L-BFGS-B: Fortran subroutines for "// &
         "large-scale bound-constrained optimization.", &
         journal="ACM Trans. Math. Softw. 1997, 23, 550-560", &
         doi="https://doi.org/10.1145/279232.279236")

      call moist_citations%add( &
         category="Solvers", &
         label=citation%lbfgsb, &
         authors="Morales, J.L., Nocedal, J.", &
         title="Remark on Algorithm 778: L-BFGS-B: Fortran subroutines "// &
         "for large-scale bound constrained optimization.", &
         journal="ACM Trans. Math. Softw. 2011, 38, Article 7", &
         doi="https://doi.org/10.1145/2049662.2049669")

      call moist_citations%add( &
         category="Solvers", &
         label=citation%fmin, &
         authors="Brent, R.P.", &
         title="Algorithms for Minimization Without Derivatives.", &
         journal="Prentice-Hall, Englewood Cliffs, NJ, 1973", &
         doi="https://maths-people.anu.edu.au/~brent/pub/pub011.html")

      call moist_citations%add( &
         category="Solvers", &
         label=citation%lsqr, &
         authors="Paige, C.C., Saunders, M.A.", &
         title="LSQR: An algorithm for sparse linear equations and "// &
         "sparse least squares.", &
         journal="ACM Trans. Math. Softw. 1982, 8, 43-71", &
         doi="https://doi.org/10.1145/355984.355989")

      call moist_citations%add( &
         category="Solvers", &
         label=citation%lsmr, &
         authors="Fong, D.C.-L., Saunders, M.A.", &
         title="LSMR: An iterative algorithm for sparse "// &
         "least-squares problems.", &
         journal="SIAM J. Sci. Comput. 2011, 33, 2950-2971", &
         doi="https://doi.org/10.1137/10079687X")

      call moist_citations%add( &
         category="Solvers", &
         label=citation%lusol, &
         authors="Gill, P.E., Murray, W., Saunders, M.A., Wright, M.H.", &
         title="Maintaining LU factors of a general sparse matrix.", &
         journal="Linear Algebra Appl. 1987, 88-89, 239-270", &
         doi="https://doi.org/10.1016/0024-3795(87)90112-1")

      call moist_citations%add( &
         category="Solvers", &
         label=citation%gmres, &
         authors="Saad, Y., Schultz, M.H.", &
         title="GMRES: A Generalized Minimal Residual Algorithm for Solving "// &
         "Nonsymmetric Linear Systems.", &
         journal="SIAM J. Sci. Stat. Comput. 1986, 7, 856-869", &
         doi="https://doi.org/10.1137/0907058")

      call moist_citations%add( &
         category="Models", &
         label=citation%pv, &
         authors="Spooner, J., Yanciw, B., Wiebe, B., Weinberg, N.", &
         title="Reaction Profiles and Energy Surfaces of Compressed Species.", &
         journal="J. Phys. Chem. A 2014, 118, 765-777", &
         doi="https://doi.org/10.1021/jp410496t")

      call moist_citations%add( &
         category="Models", &
         label=citation%pv, &
         authors="Zeller, F., Pracht, P., Neudecker, T.", &
         title="Using Conformational Sampling to Model Spectral and Structural "// &
         "Changes of Molecules at Elevated Pressures.", &
         journal="J. Phys. Chem. A 2025, 129, 2108-2116", &
         doi="https://doi.org/10.1021/acs.jpca.4c08065")

      call moist_citations%add( &
         category="Models", &
         label=citation%gostshyp, &
         authors="Scheurer, M., Dreuw, A., Epifanovsky, E., Head-Gordon, M., Stauch, T.", &
         title="Modeling Molecules under Pressure with Gaussian Potentials.", &
         journal="J. Chem. Theory Comput. 2021, 17, 583-597", &
         doi="https://doi.org/10.1021/acs.jctc.0c01212")

      call moist_citations%add( &
         category="Models", &
         label=citation%gostshyp, &
         authors="Pausch, A., Zeller, F., Neudecker, T.", &
         title="An Integral-Direct GOSTSHYP Algorithm for the Computation of "// &
         "High Pressure Effects on Molecular and Electronic Structure.", &
         journal="J. Chem. Theory Comput. 2025, 21, 747-761", &
         doi="https://doi.org/10.1021/acs.jctc.4c01502")

      call moist_citations%add( &
         category="Models", &
         label=citation%cpcm, &
         authors="Barone, V., Cossi, M.", &
         title="Quantum Calculation of Molecular Energies and Energy Gradients "// &
         "in Solution by a Conductor Solvent Model.", &
         journal="J. Phys. Chem. A 1998, 102, 1995-2001", &
         doi="https://doi.org/10.1021/jp9716997")

      call moist_citations%add( &
         category="Models", &
         label=citation%cpcm, &
         authors="Cossi, M., Rega, N., Scalmani, G., Barone, V.", &
         title="Energies, structures, and electronic properties of molecules "// &
         "in solution with the C-PCM solvation model.", &
         journal="J. Comput. Chem. 2003, 24, 669-681", &
         doi="https://doi.org/10.1002/jcc.10189")

      call moist_citations%add( &
         category="Models", &
         label=citation%cosmo, &
         authors="Klamt, A., Schuurmann, G.", &
         title="COSMO: a new approach to dielectric screening in solvents "// &
         "with explicit expressions for the screening energy and its gradient.", &
         journal="J. Chem. Soc., Perkin Trans. 2 1993, 799-805", &
         doi="https://doi.org/10.1039/P29930000799")

      call moist_citations%add( &
         category="Models", &
         label=citation%iefpcm, &
         authors="Cances, E., Mennucci, B., Tomasi, J.", &
         title="A new integral equation formalism for the polarizable continuum model: "// &
         "Theoretical background and applications to isotropic and anisotropic dielectrics.", &
         journal="J. Chem. Phys. 1997, 107, 3032-3041", &
         doi="https://doi.org/10.1063/1.474659")

      call moist_citations%add( &
         category="Models", &
         label=citation%iefpcm, &
         authors="Mennucci, B., Cances, E., Tomasi, J.", &
         title="Evaluation of Solvent Effects in Isotropic and Anisotropic Dielectrics "// &
         "and in Ionic Solutions with a Unified Integral Equation Method: "// &
         "Theoretical Bases, Computational Implementation, and Numerical Applications.", &
         journal="J. Phys. Chem. B 1997, 101, 10506-10517", &
         doi="https://doi.org/10.1021/jp971959k")

      call moist_citations%add( &
         category="Models", &
         label=citation%ssvpe, &
         authors="Chipman, D.M.", &
         title="Reaction field treatment of charge penetration.", &
         journal="J. Chem. Phys. 2000, 112, 5558-5565", &
         doi="https://doi.org/10.1063/1.481133")

      call moist_citations%add( &
         category="Models", &
         label=citation%ssvpe, &
         authors="Lange, A.W., Herbert, J.M.", &
         title="Symmetric versus asymmetric discretization of the integral equations "// &
         "in polarizable continuum solvation models.", &
         journal="Chem. Phys. Lett. 2011, 509, 77-87", &
         doi="https://doi.org/10.1016/j.cplett.2011.04.092")

      call moist_citations%add( &
         category="Models", &
         label=citation%dpcm, &
         authors="Miertus, S., Scrocco, E., Tomasi, J.", &
         title="Electrostatic interaction of a solute with a continuum. "// &
         "A direct utilization of ab initio molecular potentials for the "// &
         "prevision of solvent effects.", &
         journal="Chem. Phys. 1981, 55, 117-129", &
         doi="https://doi.org/10.1016/0301-0104(81)85090-2")

      call moist_citations%alias(citation%dpcm_original, "https://doi.org/10.1016/0301-0104(81)85090-2")

      ! LPB is the IEF form for ionic solutions
      call moist_citations%alias(citation%lpb, "https://doi.org/10.1021/jp971959k")

      call moist_citations%add( &
         category="Models", &
         label=citation%lpb, &
         authors="Cances, E., Mennucci, B.", &
         title="New applications of integral equations methods for solvation "// &
         "continuum models: ionic solutions and liquid crystals.", &
         journal="J. Math. Chem. 1998, 23, 309-326", &
         doi="https://doi.org/10.1023/A:1019133611148")

      call moist_citations%add( &
         category="Models", &
         label=citation%lpb, &
         authors="Juffer, A.H., Botta, E.F.F., van Keulen, B.A.M., "// &
         "van der Ploeg, A., Berendsen, H.J.C.", &
         title="The electric potential of a macromolecule in a solvent: "// &
         "A fundamental approach.", &
         journal="J. Comput. Phys. 1991, 97, 144-171", &
         doi="https://doi.org/10.1016/0021-9991(91)90043-K")

      call moist_citations%add( &
         category="Models", &
         label=citation%slic, &
         authors="Bardhan, J.P., Knepley, M.G.", &
         title="Communication: Modeling charge-sign asymmetric solvation free "// &
         "energies with nonlinear boundary conditions.", &
         journal="J. Chem. Phys. 2014, 141, 131103", &
         doi="https://doi.org/10.1063/1.4897324")

      call moist_citations%add( &
         category="Models", &
         label=citation%slic, &
         authors="Molavi Tabrizi, A., Goossens, S., Mehdizadeh Rahimi, A., "// &
         "Knepley, M., Bardhan, J.P.", &
         title="Predicting solvation free energies and thermodynamics in polar "// &
         "solvents and mixtures using a solvation-layer interface condition.", &
         journal="J. Chem. Phys. 2017, 146, 094103", &
         doi="https://doi.org/10.1063/1.4977037")

      call moist_citations%add( &
         category="Models", &
         label=citation%nonlocal, &
         authors="Hildebrandt, A., Blossey, R., Rjasanow, S., "// &
         "Kohlbacher, O., Lenhof, H.-P.", &
         title="Novel Formulation of Nonlocal Electrostatics.", &
         journal="Phys. Rev. Lett. 2004, 93, 108104", &
         doi="https://doi.org/10.1103/PhysRevLett.93.108104")

      call moist_citations%add( &
         category="Models", &
         label=citation%nonlocal, &
         authors="Hildebrandt, A., Blossey, R., Rjasanow, S., "// &
         "Kohlbacher, O., Lenhof, H.-P.", &
         title="Electrostatic potentials of proteins in water: "// &
         "a structured continuum approach.", &
         journal="Bioinformatics 2007, 23, e99-e103", &
         doi="https://doi.org/10.1093/bioinformatics/btl312")

      call moist_citations%add( &
         category="Models", &
         label=citation%nonlocal, &
         authors="Kemmer, T., Rjasanow, S., Hildebrandt, A.", &
         title="NESSie.jl - Efficient and intuitive finite element and boundary "// &
         "element methods for nonlocal protein electrostatics in the Julia language.", &
         journal="J. Comput. Sci. 2018, 28, 193-203", &
         doi="https://doi.org/10.1016/j.jocs.2018.08.008")

      initialised = .true.
   end subroutine fill_citations

   !> Append one entry, growing the storage by doubling
   subroutine registry_add(self, category, label, authors, title, journal, doi)
      !> Registry instance
      class(citation_registry_type), intent(inout) :: self
      !> Broad grouping: "General", "Cavities & Discretization", "Models", or "Solvers"
      character(len=*), intent(in) :: category
      !> Registry label of the part citing it
      integer, intent(in) :: label
      !> Author list (single line)
      character(len=*), intent(in) :: authors
      !> Full title
      character(len=*), intent(in) :: title
      !> Journal, year, volume (formatted for display)
      character(len=*), intent(in) :: journal
      !> DOI URL, empty when there is none
      character(len=*), intent(in) :: doi

      type(citation_entry_type), allocatable :: entries(:)
      integer :: i

      if (.not. allocated(self%entries)) allocate (self%entries(32))
      if (self%n == size(self%entries)) then
         allocate (entries(2*self%n))
         do i = 1, self%n
            entries(i) = self%entries(i)
         end do
         call move_alloc(entries, self%entries)
      end if
      self%n = self%n + 1
      self%entries(self%n) = citation_entry_type(category=category, labels=[label], authors=authors, &
                                                 title=title, journal=journal, doi=doi)
   end subroutine registry_add

   !> Add a label to an already added entry, so one paper serves two parts
   !> and still prints once
   !>
   !> - the entry is found by its DOI; a DOI that was never added is a
   !>   programming error in fill_citations
   subroutine registry_alias(self, label, doi)
      !> Registry instance
      class(citation_registry_type), intent(inout) :: self
      !> Label to add
      integer, intent(in) :: label
      !> DOI URL of the entry
      character(len=*), intent(in) :: doi

      integer :: i

      do i = 1, self%n
         if (self%entries(i)%doi /= doi) cycle
         self%entries(i)%labels = [self%entries(i)%labels, label]
         return
      end do
      error stop "moist citations: alias of an unknown DOI "//doi
   end subroutine registry_alias

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
   !> only those cited by one of these labels
   subroutine print_citations(unit, labels)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Labels to print besides the General entries; absent prints all
      integer, intent(in), optional :: labels(:)

      integer :: i

      call init_citations()

      write (unit, "(a)") "Please include the appropriate citations when using our work:"
      write (unit, "(a)") ""

      ! The package itself, outside any category
      do i = 1, moist_citations%n
         if (moist_citations%entries(i)%category == "General") call print_entry(unit, moist_citations%entries(i), "  ")
      end do

      call print_citations_by_category(unit, "Cavities & Discretization", labels)
      call print_citations_by_category(unit, "Models", labels)
      call print_citations_by_category(unit, "Solvers", labels)

   end subroutine print_citations

   !> Print all citations matching a given category, or with `labels` only
   !> those cited by one of these labels
   !>
   !> - the category header is skipped when no entry matches
   subroutine print_citations_by_category(unit, category, labels)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Category to filter by (e.g. "General", "Cavities & Discretization", "Models")
      character(len=*), intent(in) :: category
      !> Labels to keep; absent keeps all
      integer, intent(in), optional :: labels(:)

      integer :: i
      logical :: header_done

      call init_citations()

      header_done = .false.
      do i = 1, moist_citations%n
         if (moist_citations%entries(i)%category /= category) cycle
         if (present(labels)) then
            if (.not. cited_by(moist_citations%entries(i), labels)) cycle
         end if
         if (.not. header_done) then
            call print_category_header(unit, category)
            header_done = .true.
         end if
         call print_entry(unit, moist_citations%entries(i), "  ")
      end do

   end subroutine print_citations_by_category

   !> Whether one of the labels citing the entry is among `labels`
   pure logical function cited_by(entry, labels)
      !> Citation to test
      type(citation_entry_type), intent(in) :: entry
      !> Labels to look for
      integer, intent(in) :: labels(:)

      integer :: j

      cited_by = .false.
      do j = 1, size(entry%labels)
         if (any(labels == entry%labels(j))) cited_by = .true.
      end do
   end function cited_by

end module moist_output_citations
