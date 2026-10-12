!> Unit tests for the isodensity DROP level set functions and their internal
!> cartesian-monomial GTO density evaluator
!>
!> Three parts are covered:
!>
!>   * the bare evaluator [[moist_iso_gto_type]]: the assembled density against an
!>     independent direct monomial evaluation, and each analytic derivative order
!>     against a 4-point central FD of the analytic previous order (Richardson
!>     extrapolated for the third order); plus the radial screening, which must not
!>     perturb any evaluated quantity beyond its own threshold
!>
!>   * the two isodensity LSF backends: the internal
!>     [[moist_cavity_drop_lsf_isodensity_internal_type]] (moist evaluates the
!>     density itself) and the callback
!>     [[moist_cavity_drop_lsf_isodensity_callback_type]] (a host callback does)
!>
!>   * the max_deriv zeroing contract, which both backends must honour: with
!>     set_max_deriv(1) the reported Hessian is *exactly* zero (the expensive
!>     density Hessian is never computed), with set_max_deriv(2) it is the real
!>     one
!>
!> Every layer runs on real molecules: the mstore records But14diol/1,
!> MB16-43/01, and Heavy28/h2o in STO-3G, def2-SVP, and def2-TZVP
!> Geoms from mstore, the bases and converged density matrices from the reference
!> files in `test/unit/data` (generated using PySCF)
!>
!> - Cartesian-basis SCF densities run in the cartesian representation
!> - spherical-basis SCF densities (But14diol/1 def2-SVP, Heavy28/h2o def2-TZVP)
!>   run in the spherical representation, so every per-test check covers both
module test_cavity_drop_isodensity
   use moist, only: moist_iso_basis
   use moist_cavity_drop_lsf_isodensity_param, only: moist_cavity_drop_lsf_isodensity_param_type
   use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_null_ptr, c_funloc, &
                                          c_associated, c_f_pointer, c_null_funptr
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use test_helpers, only: fd4_scalar, fd4_offsets, get_test_points, center_at_origin, rel_deviation
   use moist_utils_env, only: get_env
   use moist_model_continuum_gems_utils, only: BuildSuperStructure
   use moist_cavity_drop_lsf_isodensity_gto, only: moist_iso_gto_type, moist_iso_gto_ncart, &
                                                   moist_iso_gto_work_type, &
                                                   moist_iso_gto_nslot, moist_iso_gto_nrhs
   use moist_cavity_drop_lsf_isodensity_internal, only: &
      moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_drop_lsf_isodensity_callback, only: &
      moist_cavity_drop_lsf_isodensity_callback_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_radii, only: default_cpcm_radii, new_radii_custom_atoms, radius_type
   use moist_context, only: moist_context_type, new_context
   use testdrive, only: new_unittest, unittest_type, check, test_failed, error_type
   implicit none(type, external)
   private

   public :: collect_cavity_drop_isodensity

   integer, parameter :: ndim = 3

   !> Finite-difference step
   real(wp), parameter :: STEP_SIZE = 1.0e-3_wp
   !> Absolute finite-difference comparison threshold
   real(wp), parameter :: ABS_THR = 1.0e-10_wp
   !> Relative finite-difference comparison threshold
   real(wp), parameter :: REL_THR = 1.0e-9_wp
   !> Relative agreement of one quantity evaluated at two derivative orders; the
   !> orders contract the density along different paths, so they agree to round-off
   !> - the value is measured against [[value_magnitude]], not against itself
   real(wp), parameter :: ORDER_THR = 1.0e-13_wp

   !> Every test pairs one mstore record with one basis set
   integer, parameter :: nmolecules = 3
   integer, parameter :: nbases = 3
   integer, parameter :: ntests = 10

   !> mstore collection of each test molecule
   character(len=*), parameter :: mol_collection(nmolecules) = &
                                  [character(len=9) :: "But14diol", "MB16-43", "Heavy28"]
   !> mstore record of each test molecule
   character(len=*), parameter :: mol_record(nmolecules) = &
                                  [character(len=3) :: "1", "01", "h2o"]
   !> Reference file name stem of each test molecule
   character(len=*), parameter :: mol_tag(nmolecules) = &
                                  [character(len=11) :: "but14diol_1", "mb16_43_01", "heavy28_h2o"]
   !> Reference file name stem of each basis set
   character(len=*), parameter :: basis_tag(nbases) = &
                                  [character(len=9) :: "sto_3g", "def2_svp", "def2_tzvp"]
   !> PySCF basis-set name of each basis set
   character(len=*), parameter :: basis_name(nbases) = &
                                  [character(len=9) :: "sto-3g", "def2-svp", "def2-tzvp"]

   !> Named positions in the molecule and basis tables
   integer, parameter :: imol_but14diol = 1, imol_mb16_43 = 2, imol_heavy28 = 3
   integer, parameter :: ibasis_sto_3g = 1, ibasis_def2_svp = 2, ibasis_def2_tzvp = 3

   !> Molecule of each test
   integer, parameter :: test_mol(ntests) = &
                         [imol_but14diol, imol_but14diol, imol_but14diol, &
                          imol_mb16_43, imol_mb16_43, &
                          imol_heavy28, imol_heavy28, imol_heavy28, &
                          imol_but14diol, imol_heavy28]
   !> Basis set of each test
   integer, parameter :: test_basis(ntests) = &
                         [ibasis_sto_3g, ibasis_def2_svp, ibasis_def2_tzvp, &
                          ibasis_sto_3g, ibasis_def2_svp, &
                          ibasis_sto_3g, ibasis_def2_svp, ibasis_def2_tzvp, &
                          ibasis_def2_svp, ibasis_def2_tzvp]
   !> Whether the SCF of a test used spherical basis functions; such a density
   !> runs in the spherical representation
   logical, parameter :: test_spherical(ntests) = &
                         [.false., .false., .false., .false., .false., .false., .false., .false., &
                          .true., .true.]

   !> test used wherever one reference density is enough. Heavy28/h2o in
   !> def2-TZVP is the cheapest test that still carries f functions
   integer, parameter :: test_reference = 8
   !> Spherical-basis counterpart of [[test_reference]]
   integer, parameter :: test_reference_sph = 10

   !> Shell and Cartesian-function counts of the checked-in reference files
   !> Asserted so a regenerated file cannot silently change basis set
   integer, parameter :: expected_nshell(ntests) = &
                         [28, 66, 106, 42, 83, 5, 12, 19, 66, 19]
   integer, parameter :: expected_ncart(ntests) = &
                         [40, 140, 276, 68, 189, 7, 25, 48, 140, 48]
   !> Converged RHF energies of the checked-in reference files, in Hartree
   real(wp), parameter :: expected_energy(ntests) = [ &
                          -3.03120634898054163e+02_wp, -3.06792977091366652e+02_wp, &
                          -3.07129116858719954e+02_wp, -1.25920490619214843e+03_wp, &
                          -1.27392053922109562e+03_wp, -7.49631306328200253e+01_wp, &
                          -7.59622034508552559e+01_wp, -7.60593754022876709e+01_wp, &
                          -3.06776230758649149e+02_wp, -7.60589697028347302e+01_wp]

   !> Density isovalue defining the reference surface (Bohr^-3)
   real(wp), parameter :: rho_iso_ref = 4.0e-4_wp

   !> Level set multiplier used by both backends in the agreement tests
   real(wp), parameter :: lsf_scale = 1.0_wp/rho_iso_ref

   !> Reference density shared with the C-interoperable callback
   type(moist_iso_gto_type), save :: cb_gto
   !> Host cartesian density matrix behind [[cb_gto]]
   real(wp), allocatable, save :: cb_dcart(:, :)

   !> Packed spatial derivative multi-indices through third order
   integer, parameter :: ref_deriv_x(0:19) = [ &
                         0, 1, 0, 0, 2, 1, 1, 0, 0, 0, &
                         3, 2, 2, 1, 1, 1, 0, 0, 0, 0]
   integer, parameter :: ref_deriv_y(0:19) = [ &
                         0, 0, 1, 0, 0, 1, 0, 2, 1, 0, &
                         0, 1, 0, 2, 1, 0, 3, 2, 1, 0]
   integer, parameter :: ref_deriv_z(0:19) = [ &
                         0, 0, 0, 1, 0, 0, 1, 0, 1, 2, &
                         0, 0, 1, 0, 1, 2, 0, 1, 2, 3]

contains

   !> Register the isodensity tests
   !>
   !> @param[out] testsuite Collected unit tests
   subroutine collect_cavity_drop_isodensity(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("gto_value_reference", test_gto_value_reference), &
                  new_unittest("gto_grad_fd", test_gto_grad_fd), &
                  new_unittest("gto_hess_fd", test_gto_hess_fd), &
                  new_unittest("gto_third_fd", test_gto_third_fd), &
                  new_unittest("gto_fourth_fd", test_gto_fourth_fd), &
                  new_unittest("gto_screening_separation", test_gto_screening_separation), &
                  new_unittest("gto_screened_paths", test_gto_screened_paths), &
                  new_unittest("gto_spherical_rejects_cartesian", test_gto_spherical_rejects_cartesian), &
                  new_unittest("gto_spherical_residual_tolerance", test_gto_spherical_tolerance), &
                  new_unittest("internal_vs_callback", test_internal_vs_callback), &
                  new_unittest("max_deriv_gating_internal", test_max_deriv_internal), &
                  new_unittest("max_deriv_gating_callback", test_max_deriv_callback), &
                  new_unittest("internal_fourth_fd", test_internal_fourth_fd), &
                  new_unittest("internal_fourth_gating", test_internal_fourth_gating), &
                  new_unittest("internal_screening_equivalence", test_internal_screening_equivalence), &
                  new_unittest("exclusion_radius_never_overclaims", test_exclusion_radius), &
                  new_unittest("exclusion_radius_gated_off", test_exclusion_radius_gate), &
                  new_unittest("exclusion_radius_ignores_zeff", test_exclusion_radius_zeff), &
                  new_unittest("core_fill_removes_dip", test_core_fill_dip), &
                  new_unittest("core_fill_fd", test_core_fill_fd), &
                  new_unittest("core_fill_callback_matches_internal", test_core_fill_callback), &
                  new_unittest("core_fill_guard", test_core_fill_guard), &
                  new_unittest("core_fill_certificate", test_core_fill_certificate), &
                  new_unittest("core_fill_cavity_unchanged", test_core_fill_cavity), &
                  new_unittest("core_fill_surface_check", test_core_fill_surface), &
                  new_unittest("core_fill_surface_check_callback", test_core_fill_surface_callback), &
                  new_unittest("core_fill_skipped_anchor", test_core_fill_skipped_anchor), &
                  new_unittest("core_fill_certificate_bound", test_core_fill_certificate_bound) &
                  ]
   end subroutine collect_cavity_drop_isodensity

   !> Build one realistic molecular basis test and bind its geometry
   !>
   !> @param[out] gto       Configured basis + density
   !> @param[in]  test      test identifier
   !> @param[out] error     Set if the evaluator rejects the test
   !> @param[out] mol_out   Optional copy of the test molecule
   !> @param[out] dcart_out Optional copy of the host cartesian density
   subroutine build_test(gto, test, error, mol_out, dcart_out)
      !> Configured basis + density
      type(moist_iso_gto_type), intent(out) :: gto
      !> test identifier
      integer, intent(in) :: test
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Optional copy of the test molecule
      type(structure_type), optional, intent(out) :: mol_out
      !> Optional copy of the host cartesian density
      real(wp), allocatable, optional, intent(out) :: dcart_out(:, :)

      type(mctc_error), allocatable :: merr
      type(structure_type) :: mol
      integer, allocatable :: sh_atom(:), sh_l(:), sh_nprim(:)
      real(wp), allocatable :: exps(:), coeffs(:)
      real(wp), allocatable :: dcart(:, :)

      call test_molecule(test, mol)
      call load_test_reference(test, mol, sh_atom, sh_l, sh_nprim, exps, coeffs, &
                               dcart, error)
      if (allocated(error)) return

      call gto%init(sh_atom, sh_l, sh_nprim, exps, coeffs, merr, basis=test_representation(test))
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call gto%refresh_centers(mol)
      call check_cart_layout(gto, error)
      if (allocated(error)) return

      call check(error, gto%nshell, expected_nshell(test), &
                 more=test_label(test)//" shell count")
      if (allocated(error)) return
      call check(error, gto%ncart, expected_ncart(test), &
                 more=test_label(test)//" cartesian count")
      if (allocated(error)) return
      call check(error, gto%ncart, size(dcart, 1), &
                 more=test_label(test)//" density dimension")
      if (allocated(error)) return

      call gto%set_density(dcart, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
      if (present(mol_out)) mol_out = mol
      if (present(dcart_out)) call move_alloc(dcart, dcart_out)
   end subroutine build_test

   !> Density representation a test runs in
   pure integer function test_representation(test) result(basis)
      !> test identifier
      integer, intent(in) :: test

      basis = merge(moist_iso_basis%spherical, moist_iso_basis%cartesian, test_spherical(test))
   end function test_representation

   !> Index of the molecule a test belongs to
   pure integer function test_molecule_index(test) result(imol)
      !> test identifier
      integer, intent(in) :: test

      imol = test_mol(test)
   end function test_molecule_index

   !> Index of the basis set a test belongs to
   pure integer function test_basis_index(test) result(ibasis)
      !> test identifier
      integer, intent(in) :: test

      ibasis = test_basis(test)
   end function test_basis_index

   !> Human-readable test name used in assertion messages
   pure function test_label(test) result(label)
      !> test identifier
      integer, intent(in) :: test
      character(len=:), allocatable :: label

      label = trim(mol_collection(test_molecule_index(test)))//"/" &
              //trim(mol_record(test_molecule_index(test)))//"/" &
              //trim(basis_name(test_basis_index(test)))
   end function test_label

   !> Fetch the mstore record a test is built on
   !>
   !> @param[in]  test test identifier
   !> @param[out] mol     test geometry
   subroutine test_molecule(test, mol)
      !> test identifier
      integer, intent(in) :: test
      !> test geometry
      type(structure_type), intent(out) :: mol

      integer :: imol

      imol = test_molecule_index(test)
      call get_structure(mol, trim(mol_collection(imol)), trim(mol_record(imol)))
   end subroutine test_molecule

   !> Check the Cartesian component order assumed by the PySCF density files
   subroutine check_cart_layout(gto, error)
      type(moist_iso_gto_type), intent(in) :: gto
      type(error_type), allocatable, intent(out) :: error

      integer :: shell, component, lx, ly, lz

      do shell = 1, gto%nshell
         component = gto%sh_coff(shell)
         do lx = gto%sh_l(shell), 0, -1
            do ly = gto%sh_l(shell) - lx, 0, -1
               lz = gto%sh_l(shell) - lx - ly
               component = component + 1
               call check(error, gto%comp_l(1, component), lx, more="PySCF Cartesian lx order")
               if (allocated(error)) return
               call check(error, gto%comp_l(2, component), ly, more="PySCF Cartesian ly order")
               if (allocated(error)) return
               call check(error, gto%comp_l(3, component), lz, more="PySCF Cartesian lz order")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine check_cart_layout

   !> Read the checked-in PySCF reference file belonging to a test
   !>
   !> @param[in]  test  test identifier
   !> @param[in]  mol      test geometry from mstore
   !> @param[out] sh_atom  Per-shell atom index
   !> @param[out] sh_l     Per-shell angular momentum
   !> @param[out] sh_nprim Per-shell primitive count
   !> @param[out] exps     Primitive exponents
   !> @param[out] coeffs   Primitive contraction coefficients
   !> @param[out] dcart    Cartesian density matrix
   !> @param[out] error    Set if the file is missing, truncated, or mismatched
   subroutine load_test_reference(test, mol, sh_atom, sh_l, sh_nprim, exps, &
                                  coeffs, dcart, error)
      !> test identifier
      integer, intent(in) :: test
      !> test geometry from mstore
      type(structure_type), intent(in) :: mol
      !> Per-shell atom index
      integer, allocatable, intent(out) :: sh_atom(:)
      !> Per-shell angular momentum
      integer, allocatable, intent(out) :: sh_l(:)
      !> Per-shell primitive count
      integer, allocatable, intent(out) :: sh_nprim(:)
      !> Primitive exponents
      real(wp), allocatable, intent(out) :: exps(:)
      !> Primitive contraction coefficients
      real(wp), allocatable, intent(out) :: coeffs(:)
      !> Cartesian density matrix
      real(wp), allocatable, intent(out) :: dcart(:, :)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Tolerance on the reference geometry, in Bohr
      real(wp), parameter :: THR_GEOMETRY = 1.0e-12_wp

      character(len=:), allocatable :: filename, label
      character(len=128) :: header
      character(len=32) :: collection, record, basis
      integer, allocatable :: number(:)
      real(wp), allocatable :: xyz(:, :), triangle(:)
      real(wp) :: energy
      integer :: unit, stat, nat, ncart, nelectron, nshell, nprim
      integer :: i, j, iat, shell, offset, ientry

      filename = test_filename(test)
      label = test_label(test)
      ! Claim and release the unit under a critical section
      !$omp critical(moist_test_scratch_unit)
      open (newunit=unit, file=filename, status="old", action="read", iostat=stat)
      !$omp end critical(moist_test_scratch_unit)
      if (stat /= 0) then
         call test_failed(error, "cannot open PySCF density test: "//filename)
         return
      end if

      read (unit, "(a)", iostat=stat) header
      if (stat == 0) read (unit, *, iostat=stat) collection, record, basis
      if (stat == 0) read (unit, *, iostat=stat) nat, ncart, nelectron, energy, nshell, nprim
      if (stat == 0 .and. (nat < 1 .or. ncart < 1 .or. nshell < 1 .or. nprim < nshell)) stat = 1
      if (stat /= 0) then
         !$omp critical(moist_test_scratch_unit)
         close (unit)
         !$omp end critical(moist_test_scratch_unit)
         call test_failed(error, "cannot read PySCF density test header: "//filename)
         return
      end if

      allocate (number(nat), xyz(3, nat))
      allocate (sh_atom(nshell), sh_l(nshell), sh_nprim(nshell))
      allocate (exps(nprim), coeffs(nprim))
      allocate (triangle(ncart*(ncart + 1)/2))

      do iat = 1, nat
         read (unit, *, iostat=stat) number(iat), xyz(:, iat)
         if (stat /= 0) exit
      end do
      offset = 0
      if (stat == 0) then
         do shell = 1, nshell
            read (unit, *, iostat=stat) sh_atom(shell), sh_l(shell), sh_nprim(shell)
            if (stat == 0 .and. offset + sh_nprim(shell) > nprim) stat = 1
            if (stat /= 0) exit
            read (unit, *, iostat=stat) &
               (exps(offset + i), coeffs(offset + i), i=1, sh_nprim(shell))
            if (stat /= 0) exit
            offset = offset + sh_nprim(shell)
         end do
      end if
      if (stat == 0) read (unit, *, iostat=stat) triangle
      !$omp critical(moist_test_scratch_unit)
      close (unit)
      !$omp end critical(moist_test_scratch_unit)
      if (stat /= 0) then
         call test_failed(error, "truncated PySCF density test: "//filename)
         return
      end if

      call check(error, trim(header), "PySCF 2.14.0 RHF gas-phase " &
                 //trim(merge("spherical-basis", "Cartesian      ", test_spherical(test)))//" density")
      if (allocated(error)) return
      call check(error, trim(collection)//"/"//trim(record)//"/"//trim(basis), label, &
                 more="PySCF test provenance")
      if (allocated(error)) return
      call check(error, nat, mol%nat, more=label//" atom count")
      if (allocated(error)) return
      call check(error, ncart, expected_ncart(test), more=label//" density dimension")
      if (allocated(error)) return
      call check(error, nshell, expected_nshell(test), more=label//" shell count")
      if (allocated(error)) return
      call check(error, offset, nprim, more=label//" primitive count")
      if (allocated(error)) return
      call check(error, nelectron, sum(mol%num(mol%id)), more=label//" electron count")
      if (allocated(error)) return
      call check(error, energy, expected_energy(test), thr=1.0e-10_wp, &
                 more=label//" RHF energy")
      if (allocated(error)) return

      !> The geometry the density was converged for must be the mstore record
      !> the tests evaluate on, otherwise every reference below is meaningless
      do iat = 1, nat
         call check(error, number(iat), mol%num(mol%id(iat)), &
                    more=label//" atomic number")
         if (allocated(error)) return
      end do
      call check(error, maxval(abs(xyz - mol%xyz)), 0.0_wp, thr=THR_GEOMETRY, &
                 more=label//" geometry deviates from its mstore record")
      if (allocated(error)) return

      allocate (dcart(ncart, ncart))
      ientry = 0
      do j = 1, ncart
         do i = 1, j
            ientry = ientry + 1
            dcart(i, j) = triangle(ientry)
            dcart(j, i) = triangle(ientry)
         end do
      end do
   end subroutine load_test_reference

   !> Path of the reference file belonging to a test
   !>
   !> @param[in] test test identifier
   !> @returns           Path to the checked-in PySCF reference file
   function test_filename(test) result(filename)
      !> test identifier
      integer, intent(in) :: test
      character(len=:), allocatable :: filename

      character(len=:), allocatable :: source_root

      ! get_env always returns an allocated string, so an `allocated` guard would
      ! never fire; the default is what covers an unset variable. meson exports
      ! the source root, a tester started from the project root gets the default
      source_root = get_env("MOIST_SOURCE_ROOT", default=".")
      filename = source_root//"/test/unit/data/" &
                 //trim(mol_tag(test_molecule_index(test)))//"_" &
                 //trim(basis_tag(test_basis_index(test))) &
                 //trim(merge("_sph", "    ", test_spherical(test)))//".density"
   end function test_filename

   !> Evaluate the density and its spatial derivatives at a point
   !>
   !> @param[in]  gto    Basis + density
   !> @param[in]  point  Evaluation point in Bohr
   !> @param[in]  nderiv Highest derivative order to assemble
   !> @param[out] rho    Density
   !> @param[out] drho   Density gradient
   !> @param[out] d2rho  Density Hessian; zero when nderiv < 2
   !> @param[out] d3rho  Third density derivative; zero when nderiv < 3
   !> @param[out] d4rho  Fourth density derivative; zero when nderiv < 4
   subroutine eval_at(gto, point, nderiv, rho, drho, d2rho, d3rho, d4rho)
      !> Basis + density
      type(moist_iso_gto_type), intent(in) :: gto
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Highest derivative order to assemble
      integer, intent(in) :: nderiv
      !> Density
      real(wp), intent(out) :: rho
      !> Density gradient
      real(wp), intent(out) :: drho(3)
      !> Density Hessian
      real(wp), intent(out) :: d2rho(3, 3)
      !> Third density derivative
      real(wp), intent(out) :: d3rho(3, 3, 3)
      !> Fourth density derivative
      real(wp), intent(out), optional :: d4rho(3, 3, 3, 3)

      type(moist_iso_gto_work_type) :: work

      ! The evaluator computes exactly the orders it is asked to return, so the
      ! orders this helper does not request are zeroed here for its own callers
      d2rho = 0.0_wp
      d3rho = 0.0_wp
      if (present(d4rho)) d4rho = 0.0_wp
      select case (nderiv)
      case (:1)
         call gto%eval(point, work, rho, drho)
      case (2)
         call gto%eval(point, work, rho, drho, d2rho=d2rho)
      case (3)
         call gto%eval(point, work, rho, drho, d2rho=d2rho, d3rho=d3rho)
      case default
         if (.not. present(d4rho)) error stop "eval_at: nderiv 4 requires d4rho"
         call gto%eval(point, work, rho, drho, d2rho=d2rho, d3rho=d3rho, d4rho=d4rho)
      end select
   end subroutine eval_at

   !> Independent direct density: rho = sum_cc' D_cc' g_c g_c'
   !>
   !> @param[in] gto   Basis
   !> @param[in] dcart Host cartesian density matrix
   !> @param[in] point Evaluation point in Bohr
   !> @returns         Density at the point
   function rho_reference(gto, dcart, point) result(rho)
      !> Basis
      type(moist_iso_gto_type), intent(in) :: gto
      !> Host cartesian density matrix
      real(wp), intent(in) :: dcart(:, :)
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      real(wp) :: rho

      real(wp), allocatable :: g(:)
      integer :: s, p, c, cbase, l, ncc, lx, ly, lz, i, j
      real(wp) :: d(3), u, radial

      allocate (g(gto%ncart), source=0.0_wp)
      do s = 1, gto%nshell
         l = gto%sh_l(s)
         ncc = moist_iso_gto_ncart(l)
         cbase = gto%sh_coff(s)
         d = point - gto%center(:, s)
         u = dot_product(d, d)
         radial = 0.0_wp
         do p = gto%sh_poff(s), gto%sh_poff(s + 1) - 1
            radial = radial + gto%coeffs(p)*exp(-gto%exps(p)*u)
         end do
         do c = 1, ncc
            lx = gto%comp_l(1, cbase + c)
            ly = gto%comp_l(2, cbase + c)
            lz = gto%comp_l(3, cbase + c)
            g(cbase + c) = d(1)**lx*d(2)**ly*d(3)**lz*radial
         end do
      end do

      rho = 0.0_wp
      do j = 1, gto%ncart
         do i = 1, gto%ncart
            rho = rho + dcart(i, j)*g(i)*g(j)
         end do
      end do
   end function rho_reference

   !> Analytic density value matches the independent direct evaluation
   subroutine test_gto_value_reference(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :), dcart(:, :)
      real(wp) :: rho, drho(3), d2rho(3, 3), d3rho(3, 3, 3)
      integer :: test, ip

      do test = 1, ntests
         call build_test(gto, test, error, mol, dcart)
         if (allocated(error)) return
         call get_test_points(mol, pts, 12)
         do ip = 1, size(pts, 2)
            call eval_at(gto, pts(:, ip), 1, rho, drho, d2rho, d3rho)
            call check(error, rho, rho_reference(gto, dcart, pts(:, ip)), thr=1.0e-12_wp)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_gto_value_reference

   !> Analytic gradient matches a 4-point central FD of the density
   subroutine test_gto_grad_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rho, drho(3), d2(3, 3), d3(3, 3, 3)
      real(wp) :: dg(3), dd2(3, 3), dd3(3, 3, 3)
      real(wp) :: pp(3), rpp, rp, rm, rmm, fd
      integer :: test, ip, ax

      do test = 1, ntests
         call build_test(gto, test, error, mol)
         if (allocated(error)) return
         call get_test_points(mol, pts, 12)
         do ip = 1, size(pts, 2)
            call eval_at(gto, pts(:, ip), 1, rho, drho, d2, d3)
            do ax = 1, ndim
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + 2*STEP_SIZE
               call eval_at(gto, pp, 1, rpp, dg, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + STEP_SIZE
               call eval_at(gto, pp, 1, rp, dg, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - STEP_SIZE
               call eval_at(gto, pp, 1, rm, dg, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - 2*STEP_SIZE
               call eval_at(gto, pp, 1, rmm, dg, dd2, dd3)
               call fd4_scalar(rpp, rp, rm, rmm, STEP_SIZE, fd, error)
               if (allocated(error)) return
               call check(error, drho(ax), fd, thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_gto_grad_fd

   !> Analytic Hessian matches a 4-point central FD of the gradient
   subroutine test_gto_hess_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rho, drho(3), d2(3, 3), d3(3, 3, 3)
      real(wp) :: drr, dd2(3, 3), dd3(3, 3, 3)
      real(wp) :: pp(3), gpp(3), gp(3), gm(3), gmm(3), fd
      integer :: test, ip, ax, jx

      do test = 1, ntests
         call build_test(gto, test, error, mol)
         if (allocated(error)) return
         call get_test_points(mol, pts, 12)
         do ip = 1, size(pts, 2)
            call eval_at(gto, pts(:, ip), 2, rho, drho, d2, d3)
            do ax = 1, ndim
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + 2*STEP_SIZE
               call eval_at(gto, pp, 1, drr, gpp, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + STEP_SIZE
               call eval_at(gto, pp, 1, drr, gp, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - STEP_SIZE
               call eval_at(gto, pp, 1, drr, gm, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - 2*STEP_SIZE
               call eval_at(gto, pp, 1, drr, gmm, dd2, dd3)
               do jx = 1, ndim
                  call fd4_scalar(gpp(jx), gp(jx), gm(jx), gmm(jx), STEP_SIZE, fd, error)
                  if (allocated(error)) return
                  call check(error, d2(ax, jx), fd, thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_gto_hess_fd

   !> Analytic third derivative matches a Richardson-extrapolated 4-point central
   !> FD of the Hessian (steps h and h/2), which cancels the h^4 truncation term
   subroutine test_gto_third_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rho, drho(3), d2(3, 3), d3(3, 3, 3)
      real(wp) :: drr, dg(3), dd3(3, 3, 3)
      real(wp) :: pp(3), hpp(3, 3), hp(3, 3), hm(3, 3), hmm(3, 3)
      real(wp) :: hph(3, 3), hmh(3, 3), fd, fd_full, fd_half
      integer :: test, ip, ax, jx, kx

      do test = 1, ntests
         call build_test(gto, test, error, mol)
         if (allocated(error)) return
         call get_test_points(mol, pts, 12)
         do ip = 1, size(pts, 2)
            call eval_at(gto, pts(:, ip), 3, rho, drho, d2, d3)
            do ax = 1, ndim
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + 2*STEP_SIZE
               call eval_at(gto, pp, 2, drr, dg, hpp, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + STEP_SIZE
               call eval_at(gto, pp, 2, drr, dg, hp, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - STEP_SIZE
               call eval_at(gto, pp, 2, drr, dg, hm, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - 2*STEP_SIZE
               call eval_at(gto, pp, 2, drr, dg, hmm, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + 0.5_wp*STEP_SIZE
               call eval_at(gto, pp, 2, drr, dg, hph, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - 0.5_wp*STEP_SIZE
               call eval_at(gto, pp, 2, drr, dg, hmh, dd3)
               do jx = 1, ndim
                  do kx = 1, ndim
                     ! Richardson extrapolation cancels the leading fourth-order stencil error
                     call fd4_scalar(hpp(jx, kx), hp(jx, kx), hm(jx, kx), hmm(jx, kx), STEP_SIZE, &
                                     fd_full, error)
                     if (allocated(error)) return
                     call fd4_scalar(hp(jx, kx), hph(jx, kx), hmh(jx, kx), hm(jx, kx), &
                                     0.5_wp*STEP_SIZE, fd_half, error)
                     if (allocated(error)) return
                     fd = (16.0_wp*fd_half - fd_full)/15.0_wp
                     call check(error, d3(ax, jx, kx), fd, thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_gto_third_fd

   !> Analytic fourth derivative matches a 4-point central FD of the third, in
   !> both density representations
   subroutine test_gto_fourth_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rho, drho(3), d2(3, 3), d3(3, 3, 3), d4(3, 3, 3, 3)
      real(wp) :: drr, dg(3), dh(3, 3)
      real(wp) :: pp(3), dpp(3, 3, 3), dp(3, 3, 3), dm(3, 3, 3), dmm(3, 3, 3), fd
      real(wp) :: dev
      integer :: ip, ax, ix, jx, kx, itest
      integer, parameter :: tests(2) = [test_reference, test_reference_sph]

      do itest = 1, size(tests)
         call build_test(gto, tests(itest), error, mol)
         if (allocated(error)) return
         call get_test_points(mol, pts, 8)
         do ip = 1, size(pts, 2)
            call eval_at(gto, pts(:, ip), 4, rho, drho, d2, d3, d4)

            dev = 0.0_wp
            do ax = 1, ndim
               do ix = 1, ndim
                  do jx = 1, ndim
                     do kx = 1, ndim
                        dev = max(dev, rel_deviation(d4(ax, ix, jx, kx), d4(ix, ax, jx, kx)))
                        dev = max(dev, rel_deviation(d4(ax, ix, jx, kx), d4(ax, jx, ix, kx)))
                        dev = max(dev, rel_deviation(d4(ax, ix, jx, kx), d4(ax, ix, kx, jx)))
                     end do
                  end do
               end do
            end do
            call check(error, dev, 0.0_wp, thr=1.0e-14_wp, more="d4 permutation symmetry")
            if (allocated(error)) return

            do ax = 1, ndim
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + 2*STEP_SIZE
               call eval_at(gto, pp, 3, drr, dg, dh, dpp)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + STEP_SIZE
               call eval_at(gto, pp, 3, drr, dg, dh, dp)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - STEP_SIZE
               call eval_at(gto, pp, 3, drr, dg, dh, dm)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - 2*STEP_SIZE
               call eval_at(gto, pp, 3, drr, dg, dh, dmm)
               do ix = 1, ndim
                  do jx = 1, ndim
                     do kx = 1, ndim
                        call fd4_scalar(dpp(ix, jx, kx), dp(ix, jx, kx), dm(ix, jx, kx), &
                                        dmm(ix, jx, kx), STEP_SIZE, fd, error)
                        if (allocated(error)) return
                        call check(error, d4(ax, ix, jx, kx), fd, thr_abs=ABS_THR, thr_rel=REL_THR)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_gto_fourth_fd

   !> Exercise screening while a translated molecular copy crosses the cutoff
   !>
   !> Points remain near a stationary molecule while an identical, block-diagonal
   !> density is translated away. The moving copy must contribute when adjacent,
   !> become fully inactive at large separation, and stay within the screening
   !> error bound throughout the transition
   subroutine test_gto_screening_separation(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: base_gto, pair_gto, exact_gto, screened_gto
      type(structure_type) :: center_mol, pair_mol
      real(wp), allocatable :: points(:, :), base_dcart(:, :)
      real(wp) :: rho_e, grad_e(3), hess_e(3, 3), third_e(3, 3, 3)
      real(wp) :: rho_s, grad_s(3), hess_s(3, 3), third_s(3, 3, 3)
      real(wp) :: rho_center, grad_center(3), hess_center(3, 3), third_center(3, 3, 3)
      real(wp) :: threshold, gap, gap_step, bound
      logical :: adjacent_active, adjacent_effect, far_inactive
      integer :: test, ithreshold, istep, ip

      real(wp), parameter :: thresholds(6) = [ &
                             1.0e-8_wp, 1.0e-9_wp, 1.0e-10_wp, &
                             1.0e-11_wp, 1.0e-12_wp, 1.0e-13_wp]
      integer, parameter :: nsteps = 140
      !> Clearance the widest shell must still gain at the end of the scan, in
      !> Bohr; the test points themselves already sit up to `pad` outside the box
      real(wp), parameter :: gap_margin = 10.0_wp

      do test = 1, ntests
         call build_test(base_gto, test, error, center_mol, base_dcart)
         if (allocated(error)) return
         call center_at_origin(center_mol)
         call base_gto%refresh_centers(center_mol)
         call get_test_points(center_mol, points, 12)

         !> The scan has to end beyond the widest shell reach of this test,
         !> which spans a factor of four between STO-3G carbon and the diffuse
         !> def2-TZVP sodium shells of MB16-43/01
         gap_step = (base_gto%reach(minval(thresholds)) + gap_margin)/real(nsteps, wp)

         do ithreshold = 1, size(thresholds)
            threshold = thresholds(ithreshold)
            adjacent_active = .false.
            adjacent_effect = .false.
            far_inactive = .true.

            do istep = 0, nsteps
               gap = real(istep, wp)*gap_step
               call build_density_pair(base_gto, base_dcart, center_mol, gap, pair_gto, pair_mol, error)
               if (allocated(error)) return
               exact_gto = pair_gto
               screened_gto = pair_gto
               call screened_gto%build_screening(threshold)
               bound = real(pair_gto%ncart, wp)**2*threshold

               do ip = 1, size(points, 2)
                  call eval_at(exact_gto, points(:, ip), 3, rho_e, grad_e, hess_e, third_e)
                  call eval_at(screened_gto, points(:, ip), 3, rho_s, grad_s, hess_s, third_s)
                  call check(error, rho_s, rho_e, thr=bound, &
                             more="screened density diverged during molecular separation")
                  if (allocated(error)) return
                  call check(error, maxval(abs(grad_s - grad_e)), 0.0_wp, thr=bound, &
                             more="screened density gradient diverged during separation")
                  if (allocated(error)) return
                  call check(error, maxval(abs(hess_s - hess_e)), 0.0_wp, thr=bound, &
                             more="screened density Hessian diverged during separation")
                  if (allocated(error)) return
                  call check(error, maxval(abs(third_s - third_e)), 0.0_wp, thr=bound, &
                             more="screened density third derivative diverged during separation")
                  if (allocated(error)) return

                  if (istep == 0) then
                     adjacent_active = adjacent_active .or. &
                                       any_moving_shell_active(screened_gto, points(:, ip), center_mol%nat)
                     call eval_at(base_gto, points(:, ip), 1, rho_center, grad_center, &
                                  hess_center, third_center)
                     adjacent_effect = adjacent_effect .or. abs(rho_e - rho_center) > 100.0_wp*threshold
                  else if (istep == nsteps) then
                     far_inactive = far_inactive .and. &
                                    .not. any_moving_shell_active(screened_gto, points(:, ip), center_mol%nat)
                  end if
               end do
            end do

            call check(error, adjacent_active, &
                       "moving density has no active shells at adjacent separation")
            if (allocated(error)) return
            call check(error, adjacent_effect, &
                       "adjacent moving density does not affect the reference density")
            if (allocated(error)) return
            call check(error, far_inactive, &
                       "moving density still has active shells at maximum separation")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_gto_screening_separation

   !> Screened contraction paths reproduce the exact dense contraction
   !>
   !> - order 1 runs the one-column run loop, order 4 the ten-column panels
   !> - candidate atoms arrive in reverse, so the active runs leave the internal
   !>   function order
   !> - one cartesian and one spherical test; the exact side runs unscreened, the
   !>   whole density in one gemm
   subroutine test_gto_screened_paths(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: exact, screened
      type(moist_iso_gto_work_type) :: work_e, work_s
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rho_e, grad_e(3), hess_e(3, 3), third_e(3, 3, 3), fourth_e(3, 3, 3, 3)
      real(wp) :: rho_s, grad_s(3), hess_s(3, 3), third_s(3, 3, 3), fourth_s(3, 3, 3, 3)
      real(wp) :: bound
      integer, allocatable :: reversed(:)
      integer :: itest, ip, iat, nscreened
      !> Screening threshold, small enough that the bound below stays tight
      real(wp), parameter :: threshold = 1.0e-14_wp
      !> But14diol/1 in def2-TZVP (cartesian) and def2-SVP (spherical)
      integer, parameter :: tests(2) = [3, 9]

      do itest = 1, size(tests)
         call build_test(exact, tests(itest), error, mol)
         if (allocated(error)) return
         screened = exact
         call screened%build_screening(threshold, max_deriv=4)
         reversed = [(iat, iat=mol%nat, 1, -1)]
         call get_test_points(mol, pts, 12)

         nscreened = 0
         do ip = 1, size(pts, 2)
            if (any(sum((spread(pts(:, ip), 2, screened%nshell) - screened%center)**2, dim=1) &
                    > screened%sh_rcut2)) nscreened = nscreened + 1

            call exact%eval(pts(:, ip), work_e, rho_e, grad_e)
            call screened%eval(pts(:, ip), work_s, rho_s, grad_s, cand_atoms=reversed)
            bound = real(exact%nfun, wp)**2*threshold
            call check(error, rho_s, rho_e, thr=bound + ORDER_THR*abs(rho_e), &
                       more="screened run-loop density")
            if (allocated(error)) return
            call check(error, maxval(abs(grad_s - grad_e)), 0.0_wp, &
                       thr=bound + ORDER_THR*maxval(abs(grad_e)), more="screened run-loop gradient")
            if (allocated(error)) return

            call exact%eval(pts(:, ip), work_e, rho_e, grad_e, d2rho=hess_e, d3rho=third_e, &
                            d4rho=fourth_e)
            call screened%eval(pts(:, ip), work_s, rho_s, grad_s, d2rho=hess_s, d3rho=third_s, &
                               cand_atoms=reversed, d4rho=fourth_s)
            call check(error, maxval(abs(hess_s - hess_e)), 0.0_wp, &
                       thr=bound + ORDER_THR*maxval(abs(hess_e)), more="screened panel Hessian")
            if (allocated(error)) return
            call check(error, maxval(abs(third_s - third_e)), 0.0_wp, &
                       thr=bound + ORDER_THR*maxval(abs(third_e)), more="screened panel third derivative")
            if (allocated(error)) return
            call check(error, maxval(abs(fourth_s - fourth_e)), 0.0_wp, &
                       thr=bound + ORDER_THR*maxval(abs(fourth_e)), more="screened panel fourth derivative")
            if (allocated(error)) return
         end do
         ! Positive control: the screened side really dropped shells
         call check(error, nscreened > 0, test_label(tests(itest))//" screening dropped no shell")
         if (allocated(error)) return
      end do
   end subroutine test_gto_screened_paths

   !> The spherical representation rejects a density from cartesian basis
   !> functions and keeps the density it had
   !>
   !> - a Cartesian-basis SCF density carries weight on the r^2 contaminants
   !> - after the rejection the installed spherical-basis density still
   !>   evaluates bitwise as before
   subroutine test_gto_spherical_rejects_cartesian(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: sph, cart
      type(moist_iso_gto_work_type) :: work
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: dcart(:, :), pts(:, :)
      real(wp) :: rho_before, rho_after, grad(3)

      call build_test(sph, test_reference_sph, error, mol)
      if (allocated(error)) return
      ! Water in def2-TZVP: 48 monomials, 43 spherical functions
      call check(error, sph%nfun, 43, more="spherical function count")
      if (allocated(error)) return
      call get_test_points(mol, pts, 4)
      call sph%eval(pts(:, 1), work, rho_before, grad)

      ! Same basis, density of the Cartesian-basis SCF
      call build_test(cart, test_reference, error, dcart_out=dcart)
      if (allocated(error)) return
      call sph%set_density(dcart, merr)
      call check(error, allocated(merr), "spherical representation accepted a cartesian-basis density")
      if (allocated(error)) return
      call check(error, index(merr%message, "outside the spherical subspace") > 0, &
                 "unexpected diagnostic: "//merr%message)
      if (allocated(error)) return

      call check(error, sph%has_density(), "a rejected density must keep the installed one")
      if (allocated(error)) return
      call sph%eval(pts(:, 1), work, rho_after, grad)
      call check(error, rho_after, rho_before, thr=0.0_wp, &
                 more="a rejected density must leave the installed one unchanged")
   end subroutine test_gto_spherical_rejects_cartesian

   !> A spherical-basis density with an r^2 contaminant added on one d shell is
   !> rejected above the residual tolerance and accepted below it
   subroutine test_gto_spherical_tolerance(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: dcart(:, :), contaminated(:, :), r2(:)
      integer :: s, c, ishell
      !> Contaminant weights relative to max |D|: far above and far below 1e-10
      real(wp), parameter :: weights(2) = [1.0e-7_wp, 1.0e-13_wp]
      logical, parameter :: rejected(2) = [.true., .false.]
      integer :: iw

      call build_test(gto, test_reference_sph, error, dcart_out=dcart)
      if (allocated(error)) return

      ! r^2 = xx + yy + zz on the first d shell, in the host monomial layout
      ishell = 0
      do s = 1, gto%nshell
         if (gto%sh_l(s) == 2) then
            ishell = s
            exit
         end if
      end do
      call check(error, ishell > 0, "reference basis has no d shell")
      if (allocated(error)) return
      allocate (r2(gto%ncart), source=0.0_wp)
      do c = gto%sh_coff(ishell) + 1, gto%sh_coff(ishell + 1)
         if (maxval(gto%comp_l(:, c)) == 2) r2(c) = 1.0_wp
      end do

      do iw = 1, size(weights)
         contaminated = dcart + weights(iw)*maxval(abs(dcart)) &
                        *spread(r2, 2, gto%ncart)*spread(r2, 1, gto%ncart)
         call gto%set_density(contaminated, merr)
         call check(error, allocated(merr) .eqv. rejected(iw), &
                    "r^2 contaminant handled against the residual tolerance")
         if (allocated(error)) return
      end do
   end subroutine test_gto_spherical_tolerance

   !> Build two translated molecular copies with a block-diagonal density
   subroutine build_density_pair(base_gto, base_dcart, center_mol, gap, pair_gto, pair_mol, error)
      type(moist_iso_gto_type), intent(in) :: base_gto
      real(wp), intent(in) :: base_dcart(:, :)
      type(structure_type), intent(in) :: center_mol
      real(wp), intent(in) :: gap
      type(moist_iso_gto_type), intent(out) :: pair_gto
      type(structure_type), intent(out) :: pair_mol
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: moving_mol
      type(mctc_error), allocatable :: merr
      integer, allocatable :: sh_atom(:), sh_l(:), sh_nprim(:)
      real(wp), allocatable :: exps(:), coeffs(:), density(:, :)
      real(wp) :: shift
      integer :: n

      moving_mol = center_mol
      shift = maxval(center_mol%xyz(1, :)) - minval(center_mol%xyz(1, :)) + gap
      moving_mol%xyz(1, :) = moving_mol%xyz(1, :) + shift
      call BuildSuperStructure(center_mol, moving_mol, pair_mol, no_displacement=.true.)

      sh_nprim = base_gto%sh_poff(2:) - base_gto%sh_poff(:base_gto%nshell)
      sh_atom = [base_gto%sh_atom, base_gto%sh_atom + center_mol%nat]
      sh_l = [base_gto%sh_l, base_gto%sh_l]
      sh_nprim = [sh_nprim, sh_nprim]
      exps = [base_gto%exps, base_gto%exps]
      coeffs = [base_gto%coeffs, base_gto%coeffs]
      call pair_gto%init(sh_atom, sh_l, sh_nprim, exps, coeffs, merr, basis=base_gto%basis)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call pair_gto%refresh_centers(pair_mol)

      n = base_gto%ncart
      allocate (density(2*n, 2*n), source=0.0_wp)
      density(1:n, 1:n) = base_dcart
      density(n + 1:, n + 1:) = base_dcart
      call pair_gto%set_density(density, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine build_density_pair

   !> Whether any shell on the translated molecular copy survives screening
   logical function any_moving_shell_active(gto, point, center_natom) result(active)
      type(moist_iso_gto_type), intent(in) :: gto
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: center_natom

      integer :: shell

      active = .false.
      do shell = 1, gto%nshell
         if (gto%sh_atom(shell) <= center_natom) cycle
         if (sum((point - gto%center(:, shell))**2) <= gto%sh_rcut2(shell)) then
            active = .true.
            return
         end if
      end do
   end function any_moving_shell_active

   !> Install the module-level reference density used by the callback backend
   subroutine set_reference_density(mol, error)
      type(structure_type), intent(out) :: mol
      type(error_type), allocatable, intent(out) :: error

      call build_test(cb_gto, test_reference, error, mol, cb_dcart)
   end subroutine set_reference_density

   !> Independent rho and spatial derivatives from direct monomial products
   !>
   !> Basis-function derivatives are assembled from the Leibniz rule and the
   !> closed-form derivatives of x^l exp(-a*x^2)
   !> Density derivatives then use a second, explicit multidimensional Leibniz
   !> expansion. This covers contracted s, p, d, and f functions without calling
   !> the production GTO evaluator
   !>
   !> @param[in]  point Evaluation point in Bohr
   !> @param[out] rho   Density
   !> @param[out] drho  Density gradient
   !> @param[out] d2rho Density Hessian
   !> @param[out] d3rho Third density derivative
   subroutine rho_product_rule(point, rho, drho, d2rho, d3rho)
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Density
      real(wp), intent(out) :: rho
      !> Density gradient
      real(wp), intent(out) :: drho(3)
      !> Density Hessian
      real(wp), intent(out) :: d2rho(3, 3)
      !> Third density derivative
      real(wp), intent(out) :: d3rho(3, 3, 3)

      real(wp), allocatable :: phi(:, :)
      integer :: a, b, c

      call build_reference_phi(point, phi)
      rho = density_derivative(phi, 0, 0, 0)
      do a = 1, ndim
         drho(a) = density_derivative(phi, merge(1, 0, a == 1), &
                                      merge(1, 0, a == 2), merge(1, 0, a == 3))
         do b = 1, ndim
            d2rho(a, b) = density_derivative(phi, count([a, b] == 1), &
                                             count([a, b] == 2), count([a, b] == 3))
            do c = 1, ndim
               d3rho(a, b, c) = density_derivative(phi, count([a, b, c] == 1), &
                                                   count([a, b, c] == 2), &
                                                   count([a, b, c] == 3))
            end do
         end do
      end do
   end subroutine rho_product_rule

   !> Assemble all Cartesian basis-function derivatives through third order
   subroutine build_reference_phi(point, phi)
      real(wp), intent(in) :: point(3)
      real(wp), allocatable, intent(out) :: phi(:, :)

      real(wp) :: d(3), primitive
      integer :: s, p, c, slot, cbase, lx, ly, lz

      allocate (phi(cb_gto%ncart, 0:19), source=0.0_wp)
      do s = 1, cb_gto%nshell
         d = point - cb_gto%center(:, s)
         cbase = cb_gto%sh_coff(s)
         do c = 1, moist_iso_gto_ncart(cb_gto%sh_l(s))
            lx = cb_gto%comp_l(1, cbase + c)
            ly = cb_gto%comp_l(2, cbase + c)
            lz = cb_gto%comp_l(3, cbase + c)
            do slot = 0, 19
               do p = cb_gto%sh_poff(s), cb_gto%sh_poff(s + 1) - 1
                  primitive = monomial_gaussian_deriv(d(1), lx, cb_gto%exps(p), ref_deriv_x(slot)) &
                              *monomial_gaussian_deriv(d(2), ly, cb_gto%exps(p), ref_deriv_y(slot)) &
                              *monomial_gaussian_deriv(d(3), lz, cb_gto%exps(p), ref_deriv_z(slot))
                  phi(cbase + c, slot) = phi(cbase + c, slot) + cb_gto%coeffs(p)*primitive
               end do
            end do
         end do
      end do
   end subroutine build_reference_phi

   !> Derivative of x^l exp(-a*x^2), orders zero through three
   pure real(wp) function monomial_gaussian_deriv(x, l, exponent, order) result(value)
      real(wp), intent(in) :: x, exponent
      integer, intent(in) :: l, order

      real(wp) :: gaussian(0:3), falling
      integer :: k, m

      gaussian(0) = exp(-exponent*x*x)
      gaussian(1) = -2.0_wp*exponent*x*gaussian(0)
      gaussian(2) = (4.0_wp*exponent**2*x*x - 2.0_wp*exponent)*gaussian(0)
      gaussian(3) = (-8.0_wp*exponent**3*x**3 + 12.0_wp*exponent**2*x)*gaussian(0)
      value = 0.0_wp
      do k = 0, min(order, l)
         falling = 1.0_wp
         do m = 0, k - 1
            falling = falling*real(l - m, wp)
         end do
         value = value + real(binomial_small(order, k), wp)*falling*x**(l - k)*gaussian(order - k)
      end do
   end function monomial_gaussian_deriv

   !> Direct multidimensional Leibniz expansion of one density derivative
   function density_derivative(phi, nx, ny, nz) result(value)
      real(wp), intent(in) :: phi(:, 0:)
      integer, intent(in) :: nx, ny, nz
      real(wp) :: value

      integer :: i, j, kx, ky, kz, slot_i, slot_j, weight

      value = 0.0_wp
      do j = 1, cb_gto%ncart
         do i = 1, cb_gto%ncart
            do kz = 0, nz
               do ky = 0, ny
                  do kx = 0, nx
                     slot_i = reference_slot(kx, ky, kz)
                     slot_j = reference_slot(nx - kx, ny - ky, nz - kz)
                     weight = binomial_small(nx, kx)*binomial_small(ny, ky) &
                              *binomial_small(nz, kz)
                     value = value + real(weight, wp)*cb_dcart(i, j) &
                             *phi(i, slot_i)*phi(j, slot_j)
                  end do
               end do
            end do
         end do
      end do
   end function density_derivative

   pure integer function reference_slot(nx, ny, nz) result(slot)
      integer, intent(in) :: nx, ny, nz
      integer :: candidate

      slot = -1
      do candidate = 0, 19
         if (ref_deriv_x(candidate) == nx .and. ref_deriv_y(candidate) == ny &
             .and. ref_deriv_z(candidate) == nz) then
            slot = candidate
            return
         end if
      end do
      error stop "reference_slot: derivative order is unavailable"
   end function reference_slot

   pure integer function binomial_small(n, k) result(value)
      integer, intent(in) :: n, k

      if (k < 0 .or. k > n) then
         value = 0
      else if (k == 0 .or. k == n) then
         value = 1
      else if (n == 2) then
         value = 2
      else if (n == 3 .and. (k == 1 .or. k == 2)) then
         value = 3
      else
         value = 1
      end if
   end function binomial_small

   !> C-interoperable isodensity density callback
   !>
   !> Returns the bare density rho and its spatial derivatives; the LSF applies
   !> the isovalue, the sign convention and its multiplier
   !>
   !> d2rho/d3rho are NULL when the cavity did not request that order
   !>
   !> @param[in]  context Unused callback context
   !> @param[in]  point   Evaluation point in Bohr
   !> @param[out] rho     Electron density
   !> @param[out] drho    Density spatial gradient
   !> @param[out] d2rho   Density spatial Hessian (Fortran (3,3), or NULL)
   !> @param[out] d3rho   Density third spatial derivative (Fortran (3,3,3), or NULL)
   !> @returns            0 (this reference implementation never fails)
   function iso_reference_callback(context, point, rho, drho, d2rho, d3rho) result(status) bind(C)
      !> Unused callback context (the test keeps its data in module state)
      type(c_ptr), value :: context
      !> Evaluation point in Bohr
      real(c_double), intent(in) :: point(3)
      !> Electron density
      real(c_double), intent(out) :: rho
      !> Density spatial gradient
      real(c_double), intent(out) :: drho(3)
      !> Density spatial Hessian (Fortran (3,3), or NULL)
      type(c_ptr), value :: d2rho
      !> Density third spatial derivative (Fortran (3,3,3), or NULL)
      type(c_ptr), value :: d3rho
      !> Zero on success
      integer(c_int) :: status

      real(c_double), pointer :: hptr(:, :), tptr(:, :, :)
      real(wp) :: rho_w, drho_w(3), d2rho_w(3, 3), d3rho_w(3, 3, 3)
      logical :: want_hess, want_third

      status = 0_c_int
      if (c_associated(context)) return
      want_hess = c_associated(d2rho)
      want_third = c_associated(d3rho)

      call rho_product_rule(real(point, wp), rho_w, drho_w, d2rho_w, d3rho_w)

      rho = real(rho_w, c_double)
      drho = real(drho_w, c_double)
      if (want_hess) then
         call c_f_pointer(d2rho, hptr, [3, 3])
         hptr = real(d2rho_w, c_double)
      end if
      if (want_third) then
         call c_f_pointer(d3rho, tptr, [3, 3, 3])
         tptr = real(d3rho_w, c_double)
      end if
   end function iso_reference_callback

   !> The log-ratio exclusion radius must never over-claim
   subroutine test_exclusion_radius(error)
      type(error_type), allocatable, intent(out) :: error

      !> Fraction of the claimed radius actually probed
      real(wp), parameter :: probe_frac = 0.999_wp

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: points(:, :)
      type(mctc_error), allocatable :: lsf_err
      real(wp) :: center(ndim), probe(ndim), lsf_center, lsf_probe, r
      integer :: icase, ipoint, idir, ithr, n_certified

      ! Unscreened, and at the cavity's production screening threshold. Screening
      ! only ever drops density, so the exterior branch -- the one that divides by
      ! a *smaller* rho -- is the one that could be pushed into over-claiming
      real(wp), parameter :: thresholds(2) = [0.0_wp, 1.0e-11_wp]

      ! A deterministic spread of probe directions; no RNG, so a failure here
      ! reproduces exactly
      integer, parameter :: ndir = 14
      real(wp), parameter :: dirs(ndim, ndir) = reshape([ &
                             1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, &
                             0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
                             0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, &
                             1.0_wp, 1.0_wp, 1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, &
                             1.0_wp, -1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp, -1.0_wp, &
                             1.0_wp, 1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, 1.0_wp, &
                             1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, 1.0_wp, 1.0_wp], &
                             [ndim, ndir])

      n_certified = 0
      do icase = 1, ntests
      do ithr = 1, size(thresholds)
         call build_molecular_internal_lsf(lsf, icase, thresholds(ithr), mol, error)
         if (allocated(error)) return
         call get_test_points(mol, points)

         do ipoint = 1, size(points, 2)
            center = points(:, ipoint)

            call lsf%prepare(center, lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "LSF prepare failed: "//lsf_err%message)
               return
            end if
            call lsf%f0(lsf_center)

            r = lsf%exclusion_radius(lsf_center)
            if (r <= 0.0_wp) cycle
            n_certified = n_certified + 1

            do idir = 1, ndir
               probe = center + probe_frac*r*dirs(:, idir)/norm2(dirs(:, idir))

               call lsf%prepare(probe, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               call lsf%f0(lsf_probe)

               if (lsf_probe*lsf_center <= 0.0_wp) then
                  call test_failed(error, "exclusion radius over-claims on "// &
                                   test_label(icase)//": S changes sign inside "// &
                                   "the ball it certifies as surface-free")
                  return
               end if
            end do
         end do
      end do
      end do

      if (n_certified == 0) then
         call test_failed(error, "the isodensity LSF certified no ball at all: "// &
                          "the probe loop never ran")
      end if
   end subroutine test_exclusion_radius

   !> `log_grad_out = 0` is the off switch for a certificate that is a model
   subroutine test_exclusion_radius_gate(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf_cb
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      real(wp) :: r_out, r_in, unit_s

      call build_molecular_internal_lsf(lsf, test_reference, 0.0_wp, mol, error)
      if (allocated(error)) return

      ! Exterior (S > 0) and interior (S < 0) values of equal density contrast:
      ! rho = rho_iso/e outside, rho = e*rho_iso inside. Built from the LSF's own
      ! scale, since `S = scale (rho_iso - rho)` and the certificate divides that
      ! scale back out
      allocate (radii(mol%nat), source=2.0_wp)
      unit_s = lsf%param%scale*lsf%param%rho_iso
      r_out = lsf%exclusion_radius(unit_s*(1.0_wp - exp(-1.0_wp)))
      r_in = lsf%exclusion_radius(unit_s*(1.0_wp - exp(1.0_wp)))

      if (.not. (r_out > 0.0_wp .and. r_in > 0.0_wp)) then
         call test_failed(error, "the default parameters must certify something")
         return
      end if
      if (r_in >= r_out) then
         call test_failed(error, "the interior branch must be the more "// &
                          "conservative one: its ball can reach a nuclear cusp")
         return
      end if

      ! Same contrast, so the same |ln(rho/rho_iso)| = 1: the radii are the
      ! reciprocals of the two bounds in force
      call check(error, r_out, 1.0_wp/lsf%param%log_grad_out, thr=1.0e-12_wp, &
                 more="exterior radius is 1/log_grad_out at unit log contrast")
      if (allocated(error)) return

      ! The callback variant shares the certificate, so on the same geometry and
      ! the same parameters it must return the same radius -- it differs only in
      ! where `rho` comes from, and the certificate never asks for `rho` itself
      ! No `prepare` here, so the shared reference density is not touched and
      ! this stays independent of the `cb_gto` lock
      call lsf_cb%new(callback_ptr=c_null_funptr, context=c_null_ptr, &
         param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=lsf%param%rho_iso, &
         scale=lsf%param%scale))
      call lsf_cb%update(mol, radii)
      call check(error, lsf_cb%exclusion_radius(unit_s*(1.0_wp - exp(-1.0_wp))), &
                 r_out, thr=0.0_wp, &
                 more="both isodensity variants must certify the same radius")
      if (allocated(error)) return
      call check(error, lsf_cb%exclusion_radius(unit_s*(1.0_wp - exp(1.0_wp))), &
                 r_in, thr=0.0_wp, &
                 more="both isodensity variants must certify the same radius")
      if (allocated(error)) return

      call lsf%param%new(log_grad_out=0.0_wp)
      call check(error, lsf%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 more="log_grad_out = 0 must certify nothing")
   end subroutine test_exclusion_radius_gate

   !> The interior certificate reads the atomic numbers, never `z_eff`
   !>
   !> ECP pseudo-densities exceed the cusp bound of `z_eff` but stay below that
   !> of `Z`, so a structure with ECP charges certifies what one without does;
   !> the core fill, which does depend on `z_eff`, is off
   subroutine test_exclusion_radius_zeff(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf_cb
      type(structure_type) :: mol, mol_ecp
      integer, allocatable :: numbers(:), z_eff(:)
      real(wp), allocatable :: radii(:)
      real(wp) :: s_in, r_bare, r_ecp, r_ecp_cb

      call build_molecular_internal_lsf(lsf, test_reference, 0.0_wp, mol, error)
      if (allocated(error)) return
      allocate (radii(mol%nat), source=2.0_wp)
      lsf%param%core_fill = .false.

      ! Interior value at unit log contrast, rho = e*rho_iso
      s_in = lsf%param%scale*lsf%param%rho_iso*(1.0_wp - exp(1.0_wp))
      r_bare = lsf%exclusion_radius(s_in)

      ! The heaviest element behind a two-electron core
      numbers = mol%num(mol%id)
      z_eff = numbers
      where (numbers == maxval(numbers)) z_eff = numbers - 2
      call new(mol_ecp, numbers, mol%xyz, z_eff=z_eff)

      call lsf%update(mol_ecp, radii)
      r_ecp = lsf%exclusion_radius(s_in)

      ! No `prepare`, so the null callback is never called
      call lsf_cb%new(callback_ptr=c_null_funptr, context=c_null_ptr, param=lsf%param)
      call lsf_cb%update(mol_ecp, radii)
      r_ecp_cb = lsf_cb%exclusion_radius(s_in)

      call check(error, r_ecp, r_bare, thr=0.0_wp, &
                 more="z_eff must not enter the internal cusp bound")
      if (allocated(error)) return
      call check(error, r_ecp_cb, r_bare, thr=0.0_wp, &
                 more="z_eff must not enter the callback cusp bound")
   end subroutine test_exclusion_radius_zeff

   !* =================================================================================== *!
   !*                                  ECP core fill                                      *!
   !* =================================================================================== *!

   !> Pseudo-atom LSF whose density vanishes at its nucleus
   !>
   !> One d shell of exponent 1/2 with every cartesian component weighted once,
   !> so `rho ~ r^4` dips below `rho_iso` at the nucleus as an ECP pseudo-density
   !> does; at the default isovalue the outer surface lies near 3.4 Bohr
   !>
   !> @param[inout] lsf       Internal isodensity LSF
   !> @param[in]    number    Atomic number
   !> @param[in]    z_eff     Effective nuclear charge; below `number` marks an ECP
   !> @param[in]    core_fill Whether the ECP core fill is on
   !> @param[out]   mol       Single atom at the origin
   !> @param[out]   error     Set when the basis or density is rejected
   !> @param[in]    rho_iso   Isovalue (default 1e-3)
   !> @param[in]    core_radius_fraction Fill radius as a fraction of the atom's radius
   subroutine build_pseudo_atom_lsf(lsf, number, z_eff, core_fill, mol, error, rho_iso, core_radius_fraction)
      type(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: lsf
      integer, intent(in) :: number, z_eff
      logical, intent(in) :: core_fill
      type(structure_type), intent(out) :: mol
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: rho_iso
      real(wp), intent(in), optional :: core_radius_fraction

      type(mctc_error), allocatable :: merr
      type(moist_cavity_drop_lsf_isodensity_param_type) :: param
      real(wp) :: dcart(6, 6)
      integer :: i

      call new(mol, [number], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), z_eff=[z_eff])
      param = moist_cavity_drop_lsf_isodensity_param_type(core_fill=core_fill)
      if (present(rho_iso)) then
         param%rho_iso = rho_iso
         param%scale = 1.0_wp/rho_iso
      end if
      if (present(core_radius_fraction)) param%core_radius_fraction = core_radius_fraction
      call lsf%new(sh_atom=[1], sh_l=[2], sh_nprim=[1], exps=[0.5_wp], coeffs=[1.0_wp], &
                   error=merr, param=param)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      dcart = 0.0_wp
      do i = 1, 6
         dcart(i, i) = 1.0_wp
      end do
      call lsf%set_density(dcart, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call lsf%update(mol, [4.0_wp])
   end subroutine build_pseudo_atom_lsf

   !> Unit probe directions, axes and diagonals
   pure function core_fill_dirs() result(dirs)
      real(wp) :: dirs(ndim, 7)
      integer :: i

      dirs = reshape([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
                      1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp, 0.3_wp, &
                      0.2_wp, -0.7_wp, -1.0_wp], [ndim, 7])
      do i = 1, 7
         dirs(:, i) = dirs(:, i)/norm2(dirs(:, i))
      end do
   end function core_fill_dirs

   !> The fill closes the dip of an ECP nucleus and leaves everything else alone
   !>
   !> - without the fill the nucleus lies outside the cavity: the dip itself
   !> - with it, the whole fill ball lies inside
   !> - beyond the ball, value and derivatives are bitwise those without the fill
   !> - an all-electron nucleus (`z_eff = Z`) is never filled
   subroutine test_core_fill_dip(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf_off, lsf_on, lsf_ae
      type(structure_type) :: mol
      real(wp) :: dirs(ndim, 7), point(ndim), r, s_off, s_on
      integer :: idir, istep

      call build_pseudo_atom_lsf(lsf_off, 54, 8, .false., mol, error)
      if (allocated(error)) return
      call build_pseudo_atom_lsf(lsf_on, 54, 8, .true., mol, error)
      if (allocated(error)) return
      call build_pseudo_atom_lsf(lsf_ae, 54, 54, .true., mol, error)
      if (allocated(error)) return
      call lsf_off%set_max_deriv(4)
      call lsf_on%set_max_deriv(4)
      call lsf_ae%set_max_deriv(4)

      call prepare_internal(lsf_off, [0.0_wp, 0.0_wp, 0.0_wp], error)
      if (allocated(error)) return
      call lsf_off%f0(s_off)
      call check(error, s_off > 0.0_wp, more="the pseudo-atom must dip below rho_iso at its nucleus")
      if (allocated(error)) return

      dirs = core_fill_dirs()
      do idir = 1, size(dirs, 2)
         do istep = 0, 100
            r = 0.05_wp*real(istep, wp)
            point = r*dirs(:, idir)
            call prepare_internal(lsf_on, point, error)
            if (allocated(error)) return
            call lsf_on%f0(s_on)
            ! r = 2 itself rounds to either side of the cutoff
            if (r < 1.99_wp) then
               if (.not. (s_on < 0.0_wp)) then
                  call test_failed(error, "the core fill ball must lie inside the cavity")
                  return
               end if
            else if (r > 2.01_wp) then
               call check_lsf_identical(error, lsf_on, lsf_off, point, "beyond the fill cutoff")
               if (allocated(error)) return
            end if
            call check_lsf_identical(error, lsf_ae, lsf_off, point, "all-electron nucleus")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_core_fill_dip

   !> Value and every spatial derivative of two LSFs agree bitwise at a point
   subroutine check_lsf_identical(error, lsf_a, lsf_b, point, what)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: lsf_a, lsf_b
      real(wp), intent(in) :: point(ndim)
      character(len=*), intent(in) :: what

      real(wp) :: v(2), g(ndim, 2), h(ndim, ndim, 2), t(ndim, ndim, ndim, 2)
      real(wp) :: f4(ndim, ndim, ndim, ndim, 2)

      call prepare_internal(lsf_a, point, error)
      if (allocated(error)) return
      call lsf_a%f3_rrr(v(1), g(:, 1), h(:, :, 1), t(:, :, :, 1))
      call lsf_a%f4_rrrr(f4(:, :, :, :, 1))
      call prepare_internal(lsf_b, point, error)
      if (allocated(error)) return
      call lsf_b%f3_rrr(v(2), g(:, 2), h(:, :, 2), t(:, :, :, 2))
      call lsf_b%f4_rrrr(f4(:, :, :, :, 2))
      if (v(1) /= v(2) .or. any(g(:, 1) /= g(:, 2)) .or. any(h(:, :, 1) /= h(:, :, 2)) &
          .or. any(t(:, :, :, 1) /= t(:, :, :, 2)) .or. any(f4(:, :, :, :, 1) /= f4(:, :, :, :, 2))) then
         call test_failed(error, "core fill changed the level set: "//what)
      end if
   end subroutine check_lsf_identical

   !> Every derivative of the filled level set matches a 4-point central FD of
   !> the order below, where the fill dominates and where it fades out
   subroutine test_core_fill_fd(error)
      type(error_type), allocatable, intent(out) :: error

      !> FD step, small against the 2 Bohr scale of the fill
      real(wp), parameter :: h = 1.0e-4_wp
      !> Relative tolerance against the largest entry of each derivative
      real(wp), parameter :: rtol = 1.0e-8_wp
      real(wp), parameter :: radii(4) = [0.25_wp, 0.6_wp, 0.95_wp, 1.4_wp]

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf, lsf_off
      type(structure_type) :: mol
      real(wp) :: dirs(ndim, 7), x0(ndim), x(ndim), fd, ref, v_off
      real(wp) :: v0, g0(ndim), h0(ndim, ndim), t0(ndim, ndim, ndim), q0(ndim, ndim, ndim, ndim)
      real(wp) :: v(4), g(ndim, 4), hh(ndim, ndim, 4), t(ndim, ndim, ndim, 4)
      integer :: ir, ax, k, i, j, l

      call build_pseudo_atom_lsf(lsf, 54, 8, .true., mol, error)
      if (allocated(error)) return
      call build_pseudo_atom_lsf(lsf_off, 54, 8, .false., mol, error)
      if (allocated(error)) return
      call lsf%set_max_deriv(4)
      dirs = core_fill_dirs()

      do ir = 1, size(radii)
         x0 = radii(ir)*dirs(:, 7)
         call prepare_internal(lsf_off, x0, error)
         if (allocated(error)) return
         call lsf_off%f0(v_off)
         call prepare_internal(lsf, x0, error)
         if (allocated(error)) return
         call lsf%f3_rrr(v0, g0, h0, t0)
         call lsf%f4_rrrr(q0)
         call check(error, v0 /= v_off, more="no fill at the probe: the derivatives checked are the density's alone")
         if (allocated(error)) return
         do ax = 1, ndim
            do k = 1, 4
               x = x0
               x(ax) = x0(ax) + fd4_offsets(k)*h
               call prepare_internal(lsf, x, error)
               if (allocated(error)) return
               call lsf%f3_rrr(v(k), g(:, k), hh(:, :, k), t(:, :, :, k))
            end do
            ref = maxval(abs(g0))
            call fd4_scalar(v(1), v(2), v(3), v(4), h, fd, error)
            if (allocated(error)) return
            call check(error, g0(ax), fd, thr=rtol*ref, more="core fill gradient")
            if (allocated(error)) return
            ref = maxval(abs(h0))
            do i = 1, ndim
               call fd4_scalar(g(i, 1), g(i, 2), g(i, 3), g(i, 4), h, fd, error)
               if (allocated(error)) return
               call check(error, h0(i, ax), fd, thr=rtol*ref, more="core fill Hessian")
               if (allocated(error)) return
            end do
            ref = maxval(abs(t0))
            do j = 1, ndim
               do i = 1, ndim
                  call fd4_scalar(hh(i, j, 1), hh(i, j, 2), hh(i, j, 3), hh(i, j, 4), h, fd, error)
                  if (allocated(error)) return
                  call check(error, t0(i, j, ax), fd, thr=rtol*ref, more="core fill third derivative")
                  if (allocated(error)) return
               end do
            end do
            ref = maxval(abs(q0))
            do l = 1, ndim
               do j = 1, ndim
                  do i = 1, ndim
                     call fd4_scalar(t(i, j, l, 1), t(i, j, l, 2), t(i, j, l, 3), t(i, j, l, 4), &
                                     h, fd, error)
                     if (allocated(error)) return
                     call check(error, q0(i, j, l, ax), fd, thr=rtol*ref, more="core fill fourth derivative")
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_core_fill_fd

   !> Constant density with vanishing derivatives, so the callback LSF shows the fill alone
   function constant_density_callback(context, point, rho, drho, d2rho, d3rho) result(status) bind(C)
      type(c_ptr), value :: context
      real(c_double), intent(in) :: point(3)
      real(c_double), intent(out) :: rho
      real(c_double), intent(out) :: drho(3)
      type(c_ptr), value :: d2rho
      type(c_ptr), value :: d3rho
      integer(c_int) :: status

      real(c_double), pointer :: hptr(:, :), tptr(:, :, :)

      status = 0_c_int
      if (c_associated(context) .or. any(point /= point)) status = 1_c_int
      rho = 1.0_c_double
      drho = 0.0_c_double
      if (c_associated(d2rho)) then
         call c_f_pointer(d2rho, hptr, [3, 3])
         hptr = 0.0_c_double
      end if
      if (c_associated(d3rho)) then
         call c_f_pointer(d3rho, tptr, [3, 3, 3])
         tptr = 0.0_c_double
      end if
   end function constant_density_callback

   !> The callback LSF adds the same fill as the internal one
   !>
   !> Each variant's fill is its filled minus its unfilled level set; compared
   !> where the fill dominates, so the internal difference does not cancel
   subroutine test_core_fill_callback(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: int_on, int_off
      type(moist_cavity_drop_lsf_isodensity_callback_type) :: cb_on, cb_off
      type(structure_type) :: mol
      type(mctc_error), allocatable :: lsf_err
      real(wp) :: dirs(ndim, 7), point(ndim), radii(4)
      real(wp) :: vi(2), gi(ndim, 2), hi(ndim, ndim, 2), ti(ndim, ndim, ndim, 2)
      real(wp) :: vc(2), gc(ndim, 2), hc(ndim, ndim, 2), tc(ndim, ndim, ndim, 2)
      integer :: ir, idir

      call build_pseudo_atom_lsf(int_on, 54, 8, .true., mol, error)
      if (allocated(error)) return
      call build_pseudo_atom_lsf(int_off, 54, 8, .false., mol, error)
      if (allocated(error)) return
      call cb_on%new(callback_ptr=c_funloc(constant_density_callback), context=c_null_ptr, &
                     param=moist_cavity_drop_lsf_isodensity_param_type(core_fill=.true.))
      call cb_off%new(callback_ptr=c_funloc(constant_density_callback), context=c_null_ptr, &
                      param=moist_cavity_drop_lsf_isodensity_param_type(core_fill=.false.))
      call cb_on%update(mol, [4.0_wp])
      call cb_off%update(mol, [4.0_wp])
      call int_on%set_max_deriv(3)
      call int_off%set_max_deriv(3)
      call cb_on%set_max_deriv(3)
      call cb_off%set_max_deriv(3)

      radii = [0.15_wp, 0.4_wp, 0.65_wp, 0.85_wp]
      dirs = core_fill_dirs()
      do idir = 1, size(dirs, 2)
         do ir = 1, size(radii)
            point = radii(ir)*dirs(:, idir)
            call prepare_internal(int_on, point, error)
            if (allocated(error)) return
            call int_on%f3_rrr(vi(1), gi(:, 1), hi(:, :, 1), ti(:, :, :, 1))
            call prepare_internal(int_off, point, error)
            if (allocated(error)) return
            call int_off%f3_rrr(vi(2), gi(:, 2), hi(:, :, 2), ti(:, :, :, 2))
            call cb_on%prepare(point, lsf_err)
            if (.not. allocated(lsf_err)) call cb_on%f3_rrr(vc(1), gc(:, 1), hc(:, :, 1), tc(:, :, :, 1))
            if (.not. allocated(lsf_err)) call cb_off%prepare(point, lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "callback LSF prepare failed: "//lsf_err%message)
               return
            end if
            call cb_off%f3_rrr(vc(2), gc(:, 2), hc(:, :, 2), tc(:, :, :, 2))

            call check(error, vi(1) /= vi(2), more="no fill at the probe: both variants would agree on zero")
            if (allocated(error)) return
            call check(error, vc(1) - vc(2), vi(1) - vi(2), thr=1.0e-10_wp*abs(vi(1) - vi(2)), &
                       more="callback and internal core fill values")
            if (allocated(error)) return
            call check(error, maxval(abs((gc(:, 1) - gc(:, 2)) - (gi(:, 1) - gi(:, 2)))) <= &
                       1.0e-10_wp*maxval(abs(gi(:, 1) - gi(:, 2))), more="core fill gradients")
            if (allocated(error)) return
            call check(error, maxval(abs((hc(:, :, 1) - hc(:, :, 2)) - (hi(:, :, 1) - hi(:, :, 2)))) <= &
                       1.0e-10_wp*maxval(abs(hi(:, :, 1) - hi(:, :, 2))), more="core fill Hessians")
            if (allocated(error)) return
            call check(error, maxval(abs((tc(:, :, :, 1) - tc(:, :, :, 2)) - (ti(:, :, :, 1) - ti(:, :, :, 2)))) <= &
                       1.0e-10_wp*maxval(abs(ti(:, :, :, 1) - ti(:, :, :, 2))), more="core fill third derivatives")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_core_fill_callback

   !> A surface entering the fill ball aborts the evaluation
   !>
   !> At a large isovalue the pseudo-atom's outer surface lies inside the fill
   !> cutoff, where the fill would move it
   subroutine test_core_fill_guard(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(structure_type) :: mol
      type(mctc_error), allocatable :: lsf_err

      call build_pseudo_atom_lsf(lsf, 54, 8, .true., mol, error, rho_iso=0.4_wp)
      if (allocated(error)) return
      call lsf%prepare([1.95_wp, 0.0_wp, 0.0_wp], lsf_err)
      call check(error, allocated(lsf_err), more="a surface inside the fill ball must be refused")
      if (allocated(error)) return
      call check(error, index(lsf_err%message, "ECP core fill of atom 1") > 0, &
                 more="the refusal names the filled atom: "//lsf_err%message)
      if (allocated(error)) return

      call build_pseudo_atom_lsf(lsf, 54, 8, .false., mol, error, rho_iso=0.4_wp)
      if (allocated(error)) return
      call lsf%prepare([1.95_wp, 0.0_wp, 0.0_wp], lsf_err)
      call check(error, .not. allocated(lsf_err), more="without the fill the point is ordinary")
   end subroutine test_core_fill_guard

   !> The exclusion certificate covers the steepness of the fill
   !>
   !> A light pseudo-atom's cusp bound is weaker than the fill's log-gradient, so
   !> the fill must shrink the interior ball; every certified ball stays free of
   !> surface
   subroutine test_core_fill_certificate(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: probe_frac = 0.999_wp

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf_on, lsf_off
      type(structure_type) :: mol
      real(wp) :: dirs(ndim, 7), center(ndim), probe(ndim), s_in, s_c, s_p, r
      integer :: idir, jdir, istep, n_certified

      call build_pseudo_atom_lsf(lsf_on, 10, 2, .true., mol, error)
      if (allocated(error)) return
      call build_pseudo_atom_lsf(lsf_off, 10, 2, .false., mol, error)
      if (allocated(error)) return

      s_in = lsf_on%param%scale*lsf_on%param%rho_iso*(1.0_wp - exp(1.0_wp))
      if (.not. (lsf_on%exclusion_radius(s_in) < lsf_off%exclusion_radius(s_in))) then
         call test_failed(error, "the fill's log-gradient must bound the interior ball")
         return
      end if

      dirs = core_fill_dirs()
      n_certified = 0
      do idir = 1, size(dirs, 2)
         do istep = 0, 45
            center = 0.1_wp*real(istep, wp)*dirs(:, idir)
            call prepare_internal(lsf_on, center, error)
            if (allocated(error)) return
            call lsf_on%f0(s_c)
            r = lsf_on%exclusion_radius(s_c)
            if (r <= 0.0_wp) cycle
            n_certified = n_certified + 1
            do jdir = 1, size(dirs, 2)
               probe = center + probe_frac*r*dirs(:, jdir)
               call prepare_internal(lsf_on, probe, error)
               if (allocated(error)) return
               call lsf_on%f0(s_p)
               if (s_p*s_c <= 0.0_wp) then
                  call test_failed(error, "exclusion radius over-claims around a filled nucleus")
                  return
               end if
            end do
         end do
      end do
      call check(error, n_certified > 0, more="no ball was certified: the probe loop never ran")
   end subroutine test_core_fill_certificate

   !> A DROP cavity is bitwise unchanged when a nucleus is marked as ECP
   !>
   !> Water's surface stays beyond the fill cutoff around its oxygen, so neither
   !> the surface nor any projection iterate may see the fill
   subroutine test_core_fill_cavity(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(structure_type) :: mol, mol_ecp
      type(moist_context_type), target :: ctx
      type(cavity_type_drop), allocatable :: cavity(:)
      type(mctc_error), allocatable :: cavity_error
      integer, allocatable :: numbers(:), z_eff(:)
      integer :: ivar
      !> Level set at the oxygen nucleus, unmarked and marked
      real(wp) :: s_nucleus(2)

      call build_molecular_internal_lsf(lsf, test_reference, 0.0_wp, mol, error)
      if (allocated(error)) return
      ! The production scale: DROP's tolerances are absolute in the level set
      lsf%param%scale = lsf_scale
      numbers = mol%num(mol%id)
      z_eff = numbers
      where (numbers == maxval(numbers)) z_eff = numbers - 2
      call new(mol_ecp, numbers, mol%xyz, z_eff=z_eff)

      call new_context(ctx, nthreads=0, verbosity=0)
      allocate (cavity(2))
      do ivar = 1, 2
         call new_cavity_drop(cavity(ivar), radius_model=default_cpcm_radii(), lsf_model=lsf, &
                              error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=110), ctx=ctx)
         if (.not. allocated(cavity_error)) then
            if (ivar == 1) call cavity(ivar)%update(mol, cavity_error)
            if (ivar == 2) call cavity(ivar)%update(mol_ecp, cavity_error)
         end if
         if (allocated(cavity_error)) then
            call test_failed(error, "DROP isodensity cavity failed: "//cavity_error%message)
            return
         end if
      end do

      call check(error, cavity(1)%ngrid > 0 .and. cavity(1)%ngrid == cavity(2)%ngrid, &
                 more="the ECP mark changed the number of surface points")
      if (allocated(error)) return
      if (any(cavity(1)%xyz /= cavity(2)%xyz) .or. any(cavity(1)%a /= cavity(2)%a) &
          .or. any(cavity(1)%f /= cavity(2)%f) .or. any(cavity(1)%normal0 /= cavity(2)%normal0)) then
         call test_failed(error, "the ECP core fill changed the outer surface")
         return
      end if

      ! Positive control: the marked oxygen is filled inside the cavities' own level sets
      do ivar = 1, 2
         call cavity(ivar)%lsf_model%prepare(mol%xyz(:, maxloc(numbers, 1)), cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, "prepare at the oxygen failed: "//cavity_error%message)
            return
         end if
         call cavity(ivar)%lsf_model%f0(s_nucleus(ivar))
      end do
      call check(error, s_nucleus(2) < s_nucleus(1), more="the ECP mark never reached the fill")
   end subroutine test_core_fill_cavity

   !> A converged surface point inside a fill ball aborts the cavity
   !>
   !> At `rho_iso = 0.4` the pseudo-atom's outer surface lies at 1.83 Bohr, inside
   !> the 2 Bohr fill; the filled density matches `rho_iso` there only to the
   !> projection tolerance, so the per-point guard cannot see it. A smaller
   !> `core_radius_fraction` puts the surface back outside the ball, and then the
   !> fill may not touch it
   subroutine test_core_fill_surface(error)
      type(error_type), allocatable, intent(out) :: error

      !> Atom radius, and a fill radius below the surface distance
      real(wp), parameter :: r_atom = 4.0_wp, x_small = 1.0_wp

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(structure_type) :: mol
      type(moist_context_type), target :: ctx
      type(cavity_type_drop) :: cavity(3)
      class(radius_type), allocatable :: radii
      type(mctc_error), allocatable :: cavity_error
      integer :: ivar

      call new_context(ctx, nthreads=0, verbosity=0)
      call new_radii_custom_atoms([r_atom], radii, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ! 1: default fill, 2: no fill, 3: fill shrunk below the surface
      do ivar = 1, 3
         if (ivar == 3) then
            call build_pseudo_atom_lsf(lsf, 54, 8, .true., mol, error, rho_iso=0.4_wp, &
                                       core_radius_fraction=x_small/r_atom)
         else
            call build_pseudo_atom_lsf(lsf, 54, 8, ivar == 1, mol, error, rho_iso=0.4_wp)
         end if
         if (allocated(error)) return
         call new_cavity_drop(cavity(ivar), radius_model=radii, lsf_model=lsf, error=cavity_error, &
                              param=moist_cavity_drop_parameters_type(num_leb=6, proj_level=1), ctx=ctx)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if
         call cavity(ivar)%update(mol, cavity_error)
         if (ivar == 1) then
            call check(error, allocated(cavity_error), more="a surface inside the fill ball must be refused")
            if (allocated(error)) return
            ! The per-point guard may see an iterate inside the ball first
            call check(error, index(cavity_error%message, "ECP core fill") > 0 .and. &
                       index(cavity_error%message, "atom 1") > 0, &
                       more="the refusal names the filled atom: "//cavity_error%message)
            if (allocated(error)) return
            deallocate (cavity_error)
         else if (allocated(cavity_error)) then
            call test_failed(error, "DROP isodensity cavity failed: "//cavity_error%message)
            return
         end if
      end do

      call check(error, cavity(2)%ngrid > 0 .and. cavity(2)%ngrid == cavity(3)%ngrid, &
                 more="the shrunk fill changed the number of surface points")
      if (allocated(error)) return
      call check(error, minval(norm2(cavity(3)%xyz(:, 1:cavity(3)%ngrid), dim=1)) > x_small, &
                 more="a surface point lies inside the shrunk fill")
      if (allocated(error)) return
      if (any(cavity(2)%xyz /= cavity(3)%xyz) .or. any(cavity(2)%a /= cavity(3)%a)) then
         call test_failed(error, "a fill clear of the surface changed it")
      end if
   end subroutine test_core_fill_surface

   !> The callback LSF rejects converged surface points inside its fill too
   !>
   !> Only converged points count, and `core_radius_fraction` sizes the ball; a
   !> fraction <= 0 leaves no fill
   subroutine test_core_fill_surface_callback(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf
      type(structure_type) :: mol
      type(mctc_error), allocatable :: lsf_err
      real(wp) :: xyz(ndim, 2)

      call new(mol, [54], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), z_eff=[8])
      ! One point inside the default 2 Bohr fill, one beyond it
      xyz = reshape([1.5_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.5_wp, 0.0_wp], [ndim, 2])

      call lsf%new(callback_ptr=c_funloc(constant_density_callback), context=c_null_ptr)
      call lsf%update(mol, [4.0_wp])
      call lsf%check_surface(xyz, [.false., .true.], lsf_err)
      call check(error, .not. allocated(lsf_err), more="an unconverged point is not on the surface")
      if (allocated(error)) return
      call lsf%check_surface(xyz, [.true., .true.], lsf_err)
      call check(error, allocated(lsf_err), more="a converged point inside the fill must be refused")
      if (allocated(error)) return
      call check(error, index(lsf_err%message, "surface point 1 lies 1.500 Bohr from the nucleus of atom 1") > 0, &
                 more="the refusal names the point and the atom: "//lsf_err%message)
      if (allocated(error)) return

      call lsf%new(callback_ptr=c_funloc(constant_density_callback), context=c_null_ptr, &
                   param=moist_cavity_drop_lsf_isodensity_param_type(core_radius_fraction=0.25_wp))
      call lsf%update(mol, [4.0_wp])
      call lsf%check_surface(xyz, [.true., .true.], lsf_err)
      call check(error, .not. allocated(lsf_err), more="core_radius_fraction must shrink the fill ball")
      if (allocated(error)) return

      call lsf%new(callback_ptr=c_funloc(constant_density_callback), context=c_null_ptr, &
                   param=moist_cavity_drop_lsf_isodensity_param_type(core_radius_fraction=0.0_wp))
      call lsf%update(mol, [4.0_wp])
      xyz(:, 1) = 0.0_wp
      call lsf%check_surface(xyz, [.true., .true.], lsf_err)
      call check(error, .not. allocated(lsf_err), more="a fraction of zero leaves no fill")
   end subroutine test_core_fill_surface_callback

   !> An anchor skipped for a small switching weight is no surface point
   !>
   !> Such an anchor stays at its Lebedev position until the postfilter drops it;
   !> next to a hydrogen it lies 0.64 Bohr from the marked xenon, deep inside the
   !> fill, and must neither abort the build nor change the cavity
   subroutine test_core_fill_skipped_anchor(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(structure_type) :: mol
      type(moist_context_type), target :: ctx
      type(cavity_type_drop) :: cavity(2)
      type(mctc_error), allocatable :: cavity_error
      real(wp) :: dcart(7, 7)
      !> Level set at the xenon nucleus, without and with the fill
      real(wp) :: s_nucleus(2)
      integer :: ivar, i

      call new(mol, [54, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.10_wp], [3, 2]), z_eff=[8, 1])
      dcart = 0.0_wp
      do i = 1, 7
         dcart(i, i) = 1.0_wp
      end do
      call new_context(ctx, nthreads=0, verbosity=0)

      ! 1: no fill, 2: fill
      do ivar = 1, 2
         call lsf%new(sh_atom=[1, 2], sh_l=[2, 0], sh_nprim=[1, 1], exps=[0.5_wp, 0.4_wp], &
                      coeffs=[1.0_wp, 1.0_wp], error=cavity_error, &
                      param=moist_cavity_drop_lsf_isodensity_param_type(core_fill=ivar == 2))
         if (.not. allocated(cavity_error)) call lsf%set_density(dcart, cavity_error)
         if (.not. allocated(cavity_error)) then
            call new_cavity_drop(cavity(ivar), radius_model=default_cpcm_radii(), lsf_model=lsf, &
                                 error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=6), ctx=ctx)
         end if
         if (.not. allocated(cavity_error)) call cavity(ivar)%update(mol, cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, "DROP isodensity cavity failed: "//cavity_error%message)
            return
         end if
      end do

      call check(error, cavity(1)%ngrid > 0 .and. cavity(1)%ngrid == cavity(2)%ngrid, &
                 more="the fill changed the number of surface points")
      if (allocated(error)) return
      if (any(cavity(1)%xyz /= cavity(2)%xyz) .or. any(cavity(1)%a /= cavity(2)%a)) then
         call test_failed(error, "a fill clear of the surface changed it")
         return
      end if

      ! Positive control: the second cavity's level set fills the xenon nucleus
      do ivar = 1, 2
         call cavity(ivar)%lsf_model%prepare(mol%xyz(:, 1), cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, "prepare at the xenon nucleus failed: "//cavity_error%message)
            return
         end if
         call cavity(ivar)%lsf_model%f0(s_nucleus(ivar))
      end do
      call check(error, s_nucleus(2) < s_nucleus(1), more="the fill never reached the xenon nucleus")
   end subroutine test_core_fill_skipped_anchor

   !> The fill's term in the exclusion certificate bounds `|grad g| / max(g, rho_iso)`
   !>
   !> `g` and its gradient are read off the callback level set over a constant
   !> density, fill on minus fill off; the sampled supremum must stay below the
   !> bound and come close to it for one ball. Two balls 1 Bohr apart exceed the
   !> one-ball bound, which therefore has to count the overlap
   subroutine test_core_fill_certificate_bound(error)
      type(error_type), allocatable, intent(out) :: error

      !> Samples per line
      integer, parameter :: nsample = 4000

      type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf_on, lsf_off
      type(structure_type) :: mol
      real(wp) :: sup(2), bound(2), dirs(ndim, 7), point(ndim), s_in
      integer :: icase, idir, isample

      dirs = core_fill_dirs()
      do icase = 1, 2
         if (icase == 1) then
            call new(mol, [54], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), z_eff=[8])
         else
            call new(mol, [54, 54], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]), &
                     z_eff=[8, 8])
         end if
         call lsf_on%new(callback_ptr=c_funloc(constant_density_callback), context=c_null_ptr)
         call lsf_off%new(callback_ptr=c_funloc(constant_density_callback), context=c_null_ptr, &
                          param=moist_cavity_drop_lsf_isodensity_param_type(core_fill=.false.))
         call lsf_on%update(mol, spread(4.0_wp, 1, mol%nat))
         call lsf_off%update(mol, spread(4.0_wp, 1, mol%nat))
         call lsf_on%set_max_deriv(1)
         call lsf_off%set_max_deriv(1)
         bound(icase) = lsf_on%core%log_grad_bound()

         ! Rays from the first nucleus through both balls (X = 2 Bohr)
         sup(icase) = 0.0_wp
         do idir = 1, size(dirs, 2)
            do isample = 1, nsample
               point = (3.5_wp*real(isample, wp)/real(nsample, wp))*dirs(:, idir)
               sup(icase) = max(sup(icase), fill_log_grad(lsf_on, lsf_off, point, error))
               if (allocated(error)) return
            end do
         end do
         if (icase == 1) then
            ! The callback carries the fill into its certificate
            s_in = lsf_on%param%scale*lsf_on%param%rho_iso*(1.0_wp - exp(1.0_wp))
            call check(error, lsf_on%exclusion_radius(s_in) < lsf_off%exclusion_radius(s_in), &
                       more="the callback certificate ignores the fill")
            if (allocated(error)) return
         end if
      end do

      call check(error, sup(1) <= bound(1) .and. sup(1) >= 0.95_wp*bound(1), &
                 more="one ball: sampled sup outside [0.95, 1] of the bound")
      if (allocated(error)) return
      call check(error, sup(2) <= bound(2), more="overlapping balls exceed the bound")
      if (allocated(error)) return
      call check(error, sup(2) > bound(1), more="overlapping balls never exceeded one ball's bound")
   end subroutine test_core_fill_certificate_bound

   !> `|grad g| / max(g, rho_iso)` at one point, from filled minus unfilled level sets
   function fill_log_grad(lsf_on, lsf_off, point, error) result(ratio)
      type(moist_cavity_drop_lsf_isodensity_callback_type), intent(inout) :: lsf_on, lsf_off
      real(wp), intent(in) :: point(ndim)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: ratio

      type(mctc_error), allocatable :: lsf_err
      real(wp) :: s_on, s_off, g_on(ndim), g_off(ndim), g, grad_g(ndim)

      ratio = 0.0_wp
      call lsf_on%prepare(point, lsf_err)
      if (.not. allocated(lsf_err)) call lsf_off%prepare(point, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "callback prepare failed: "//lsf_err%message)
         return
      end if
      call lsf_on%f012_r(lsf0=s_on, lsf1_r=g_on)
      call lsf_off%f012_r(lsf0=s_off, lsf1_r=g_off)
      g = -(s_on - s_off)/lsf_on%param%scale
      grad_g = -(g_on - g_off)/lsf_on%param%scale
      ratio = norm2(grad_g)/max(g, lsf_on%param%rho_iso)
   end function fill_log_grad

   !> Build the internal isodensity LSF over the module reference density
   !>
   !> @param[inout] lsf       Internal isodensity LSF
   !> @param[in]    mol       Molecular structure
   !> @param[in]    threshold Screening threshold pushed into the LSF
   !> @param[out]   error     Set if the basis or density is rejected
   subroutine build_internal_lsf(lsf, mol, threshold, error)
      !> Internal isodensity LSF
      type(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: lsf
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Screening threshold pushed into the LSF
      real(wp), intent(in) :: threshold
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: merr
      integer, allocatable :: sh_nprim(:)
      real(wp), allocatable :: radii(:)

      sh_nprim = cb_gto%sh_poff(2:) - cb_gto%sh_poff(:cb_gto%nshell)
      allocate (radii(mol%nat), source=2.0_wp)

      call lsf%new(sh_atom=cb_gto%sh_atom, sh_l=cb_gto%sh_l, sh_nprim=sh_nprim, exps=cb_gto%exps, &
         coeffs=cb_gto%coeffs, error=merr, &
         param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=rho_iso_ref, scale=lsf_scale))
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call lsf%set_density(cb_dcart, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      lsf%screening_threshold = threshold
      call lsf%update(mol, radii)
   end subroutine build_internal_lsf

   !> Build an internal LSF from one of the realistic molecular basis tests
   subroutine build_molecular_internal_lsf(lsf, test, threshold, mol, error)
      type(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: lsf
      integer, intent(in) :: test
      real(wp), intent(in) :: threshold
      type(structure_type), intent(out) :: mol
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(mctc_error), allocatable :: merr
      integer, allocatable :: sh_nprim(:)
      real(wp), allocatable :: radii(:), dcart(:, :)

      call build_test(gto, test, error, mol, dcart)
      if (allocated(error)) return
      sh_nprim = gto%sh_poff(2:) - gto%sh_poff(:gto%nshell)
      allocate (radii(mol%nat), source=2.0_wp)

      call lsf%new(sh_atom=gto%sh_atom, sh_l=gto%sh_l, sh_nprim=sh_nprim, exps=gto%exps, coeffs=gto%coeffs, &
         error=merr, param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=rho_iso_ref, scale=1.0_wp, &
         basis=test_representation(test)))
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call lsf%set_density(dcart, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      lsf%screening_threshold = threshold
      call lsf%update(mol, radii)
   end subroutine build_molecular_internal_lsf

   !> Build the callback isodensity LSF over the module reference density
   !>
   !> @param[inout] lsf LSF instance
   !> @param[in]    mol Molecular structure
   subroutine build_callback_lsf(lsf, mol)
      !> LSF instance
      type(moist_cavity_drop_lsf_isodensity_callback_type), intent(inout) :: lsf
      !> Molecular structure
      type(structure_type), intent(in) :: mol

      real(wp), allocatable :: radii(:)

      allocate (radii(mol%nat), source=2.0_wp)
      call lsf%new(callback_ptr=c_funloc(iso_reference_callback), context=c_null_ptr, &
         param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=rho_iso_ref, scale=lsf_scale))
      call lsf%update(mol, radii)
   end subroutine build_callback_lsf

   !> Both isodensity backends must describe the same level set
   subroutine test_internal_vs_callback(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !$omp critical(iso_reference_density)
      call run_internal_vs_callback(error)
      !$omp end critical(iso_reference_density)
   end subroutine test_internal_vs_callback

   !> Body of [[test_internal_vs_callback]], run under the `cb_gto` lock
   subroutine run_internal_vs_callback(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf_int
      type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf_cb
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: v0_i, v0_c, v01_i, v01_c, g01_i(3), g01_c(3)
      real(wp) :: v012_i, v012_c, g012_i(3), g012_c(3), h012_i(3, 3), h012_c(3, 3)
      real(wp) :: v_i, g_i(3), h_i(3, 3)
      real(wp) :: v_c, g_c(3), h_c(3, 3)
      real(wp) :: t_i(3, 3, 3), t_c(3, 3, 3)
      real(wp) :: dev_val, dev_grad, dev_hess, dev_third, scale_ref
      type(mctc_error), allocatable :: lsf_err
      integer :: ip

      !> Roundoff-limited thresholds include headroom for the different summation
      !> orders of the dense production contraction and direct Leibniz reference
      real(wp), parameter :: THR_VAL = 1.0e-12_wp
      real(wp), parameter :: THR_GRAD = 1.0e-12_wp
      real(wp), parameter :: THR_HESS = 1.0e-11_wp
      real(wp), parameter :: THR_THIRD = 1.0e-10_wp

      call set_reference_density(mol, error)
      if (allocated(error)) return
      call build_internal_lsf(lsf_int, mol, 0.0_wp, error)
      if (allocated(error)) return
      call build_callback_lsf(lsf_cb, mol)

      call lsf_int%set_max_deriv(3)
      call lsf_cb%set_max_deriv(3)

      call get_test_points(mol, pts, 12)
      dev_val = 0.0_wp
      dev_grad = 0.0_wp
      dev_hess = 0.0_wp
      dev_third = 0.0_wp
      scale_ref = 0.0_wp
      do ip = 1, size(pts, 2)
         call lsf_int%prepare(pts(:, ip), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "internal LSF prepare failed: "//lsf_err%message)
            return
         end if
         call lsf_int%f0(v0_i)
         call lsf_int%f012_r(v01_i, g01_i)
         call lsf_int%f012_r(v012_i, g012_i, h012_i)
         call lsf_int%f3_rrr(v_i, g_i, h_i, t_i)
         call lsf_cb%prepare(pts(:, ip), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "callback LSF prepare failed: "//lsf_err%message)
            return
         end if
         call lsf_cb%f0(v0_c)
         call lsf_cb%f012_r(v01_c, g01_c)
         call lsf_cb%f012_r(v012_c, g012_c, h012_c)
         call lsf_cb%f3_rrr(v_c, g_c, h_c, t_c)

         call check(error, v0_i, v0_c, thr=THR_VAL, &
                    more="internal vs callback f0 value")
         if (allocated(error)) return
         call check(error, v01_i, v01_c, thr=THR_VAL, &
                    more="internal vs callback f012_r value without Hessian")
         if (allocated(error)) return
         call check(error, maxval(abs(g01_i - g01_c)), 0.0_wp, thr=THR_GRAD, &
                    more="internal vs callback f012_r gradient without Hessian")
         if (allocated(error)) return
         call check(error, v012_i, v012_c, thr=THR_VAL, &
                    more="internal vs callback f012_r value")
         if (allocated(error)) return
         call check(error, maxval(abs(g012_i - g012_c)), 0.0_wp, thr=THR_GRAD, &
                    more="internal vs callback f012_r gradient")
         if (allocated(error)) return
         call check(error, maxval(abs(h012_i - h012_c)), 0.0_wp, thr=THR_HESS, &
                    more="internal vs callback f012_r Hessian")
         if (allocated(error)) return

         dev_val = max(dev_val, abs(v_i - v_c))
         dev_grad = max(dev_grad, maxval(abs(g_i - g_c)))
         dev_hess = max(dev_hess, maxval(abs(h_i - h_c)))
         dev_third = max(dev_third, maxval(abs(t_i - t_c)))
         scale_ref = max(scale_ref, maxval(abs(t_i)))
      end do

      !> The test must actually produce a non-trivial level set; otherwise the
      !> agreement above would be vacuous
      call check(error, scale_ref > 1.0e-3_wp, &
                 "isodensity test produced a vanishing third derivative")
      if (allocated(error)) return

      call check(error, dev_val, 0.0_wp, thr=THR_VAL, &
                 more="internal vs callback level set value")
      if (allocated(error)) return
      call check(error, dev_grad, 0.0_wp, thr=THR_GRAD, &
                 more="internal vs callback level set gradient")
      if (allocated(error)) return
      call check(error, dev_hess, 0.0_wp, thr=THR_HESS, &
                 more="internal vs callback level set Hessian")
      if (allocated(error)) return
      call check(error, dev_third, 0.0_wp, thr=THR_THIRD, &
                 more="internal vs callback level set third derivative")
   end subroutine run_internal_vs_callback

   !> max_deriv must gate what the internal backend caches, without disturbing
   !> the orders that are still requested
   subroutine test_max_deriv_internal(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: v1, g1(3)
      real(wp) :: v2, g2(3), h2(3, 3)
      type(mctc_error), allocatable :: lsf_err
      integer :: test, ip
      logical :: any_nonzero

      do test = 1, ntests
         call build_molecular_internal_lsf(lsf, test, 0.0_wp, mol, error)
         if (allocated(error)) return
         call get_test_points(mol, pts, 12)
         call append_surface_points(lsf, lsf%param, mol, pts, error)
         if (allocated(error)) return

         any_nonzero = .false.
         do ip = 1, size(pts, 2)
            !> Value+gradient only: the Hessian is deliberately not requested
            call lsf%set_max_deriv(1)
            call lsf%prepare(pts(:, ip), lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "internal LSF prepare failed: "//lsf_err%message)
               return
            end if
            call check(error, lsf%prepared_deriv, 1, &
                       more="internal molecular test cached the wrong order at max_deriv=1")
            if (allocated(error)) return
            call lsf%f012_r(v1, g1)

            call lsf%set_max_deriv(2)
            call lsf%prepare(pts(:, ip), lsf_err)
            call check(error, lsf%prepared_deriv, 2, &
                       more="internal molecular test cached the wrong order at max_deriv=2")
            if (allocated(error)) return
            call lsf%f012_r(v2, g2, h2)

            call check_max_deriv_pair(error, lsf%param, v1, g1, v2, g2, h2, &
                                      "internal molecular test", any_nonzero)
            if (allocated(error)) return
         end do

         call check(error, any_nonzero, &
                    "internal molecular test returned a zero Hessian at max_deriv=2")
         if (allocated(error)) return
      end do
   end subroutine test_max_deriv_internal

   !> max_deriv must gate what the callback backend requests through the ABI
   subroutine test_max_deriv_callback(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !$omp critical(iso_reference_density)
      call run_max_deriv_callback(error)
      !$omp end critical(iso_reference_density)
   end subroutine test_max_deriv_callback

   !> Body of [[test_max_deriv_callback]], run under the `cb_gto` lock
   subroutine run_max_deriv_callback(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: v1, g1(3)
      real(wp) :: v2, g2(3), h2(3, 3)
      type(mctc_error), allocatable :: lsf_err
      integer :: ip
      logical :: any_nonzero

      call set_reference_density(mol, error)
      if (allocated(error)) return
      call build_callback_lsf(lsf, mol)
      call get_test_points(mol, pts, 12)
      call append_surface_points(lsf, lsf%param, mol, pts, error)
      if (allocated(error)) return

      any_nonzero = .false.
      do ip = 1, size(pts, 2)
         !> Value+gradient only: the Hessian is deliberately not requested
         call lsf%set_max_deriv(1)
         call lsf%prepare(pts(:, ip), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "callback LSF prepare failed: "//lsf_err%message)
            return
         end if
         call check(error, lsf%prepared_deriv, 1, &
                    more="callback isodensity LSF cached the wrong order at max_deriv=1")
         if (allocated(error)) return
         call lsf%f012_r(v1, g1)

         call lsf%set_max_deriv(2)
         call lsf%prepare(pts(:, ip), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "callback LSF prepare failed: "//lsf_err%message)
            return
         end if
         call check(error, lsf%prepared_deriv, 2, &
                    more="callback isodensity LSF cached the wrong order at max_deriv=2")
         if (allocated(error)) return
         call lsf%f012_r(v2, g2, h2)

         call check_max_deriv_pair(error, lsf%param, v1, g1, v2, g2, h2, "callback", any_nonzero)
         if (allocated(error)) return
      end do

      call check(error, any_nonzero, &
                 "callback isodensity LSF returned a zero Hessian at max_deriv=2")
   end subroutine run_max_deriv_callback

   !> Shared assertions for one max_deriv=1 / max_deriv=2 evaluation pair
   !>
   !> @param[out]   error       Set on contract violation
   !> @param[in]    param       Isodensity settings of the evaluated LSF
   !> @param[in]    v1          Value at max_deriv = 1
   !> @param[in]    g1          Gradient at max_deriv = 1
   !> @param[in]    v2          Value at max_deriv = 2
   !> @param[in]    g2          Gradient at max_deriv = 2
   !> @param[in]    h2          Hessian at max_deriv = 2
   !> @param[in]    label       Backend name used in failure messages
   !> @param[inout] any_nonzero Set once a non-zero max_deriv=2 Hessian is seen
   subroutine check_max_deriv_pair(error, param, v1, g1, v2, g2, h2, label, any_nonzero)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Isodensity settings of the evaluated LSF
      type(moist_cavity_drop_lsf_isodensity_param_type), intent(in) :: param
      !> Value at max_deriv = 1
      real(wp), intent(in) :: v1
      !> Gradient at max_deriv = 1
      real(wp), intent(in) :: g1(3)
      !> Value at max_deriv = 2
      real(wp), intent(in) :: v2
      !> Gradient at max_deriv = 2
      real(wp), intent(in) :: g2(3)
      !> Hessian at max_deriv = 2
      real(wp), intent(in) :: h2(3, 3)
      !> Backend name used in failure messages
      character(len=*), intent(in) :: label
      !> Set once a non-zero max_deriv=2 Hessian is seen
      logical, intent(inout) :: any_nonzero

      !> Lowering max_deriv must not disturb the orders that are still requested,
      !> beyond the round-off of a different contraction path
      call check(error, v1, v2, thr=ORDER_THR*value_magnitude(param, v2), &
                 more=label//" isodensity LSF value depends on max_deriv")
      if (allocated(error)) return
      call check(error, maxval(abs(g1 - g2)), 0.0_wp, thr=ORDER_THR*maxval(abs(g2)), &
                 more=label//" isodensity LSF gradient depends on max_deriv")
      if (allocated(error)) return

      if (maxval(abs(h2)) > 0.0_wp) any_nonzero = .true.
   end subroutine check_max_deriv_pair

   !> Round-off magnitude of a level-set value, scale*max(rho, rho_iso)
   !>
   !> - S = scale*(rho_iso - rho) vanishes on the surface; contraction round-off
   !>   scales with rho itself
   !>
   !> @param[in] param  isodensity settings of the evaluated LSF
   !> @param[in] value  level-set value S
   pure function value_magnitude(param, value) result(magnitude)
      !> Isodensity settings of the evaluated LSF
      type(moist_cavity_drop_lsf_isodensity_param_type), intent(in) :: param
      !> Level-set value S
      real(wp), intent(in) :: value
      !> scale*max(rho, rho_iso), with rho recovered from S
      real(wp) :: magnitude

      magnitude = param%scale*max(param%rho_iso - value/param%scale, param%rho_iso)
   end function value_magnitude

   !> Append isosurface points, one per atom ray that crosses S = 0
   !>
   !> - Rays leave the first atoms away from the centroid; the LSF's own sign is bisected
   !> - Random box points miss the surface, where S vanishes and only rho sets round-off
   !> - Fails when no ray crosses or a point misses S = 0 beyond round-off
   !>
   !> @param[inout] lsf    isodensity LSF, left at max_deriv = 1
   !> @param[in]    param  isodensity settings of `lsf`
   !> @param[in]    mol    molecule whose atoms seed the rays
   !> @param[inout] pts    points in Bohr, (3, n), grown by the surface points
   !> @param[out]   error  set when no surface point is found
   subroutine append_surface_points(lsf, param, mol, pts, error)
      !> Isodensity LSF, left at max_deriv = 1
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Isodensity settings of `lsf`
      type(moist_cavity_drop_lsf_isodensity_param_type), intent(in) :: param
      !> Molecule whose atoms seed the rays
      type(structure_type), intent(in) :: mol
      !> Points in Bohr, grown by the surface points
      real(wp), allocatable, intent(inout) :: pts(:, :)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Rays per molecule, at most one per atom
      integer, parameter :: nray = 4
      !> Bisection steps; 60 halvings of the bracket reach position round-off
      integer, parameter :: nbisect = 60
      !> Ray bracket in Bohr from the atom: inside and outside the isosurface
      real(wp), parameter :: r_in = 0.3_wp, r_out = 10.0_wp
      real(wp) :: centroid(3), dir(3), lo, hi, mid, s_lo, s_hi, s_mid
      real(wp) :: found(3, nray)
      real(wp), allocatable :: grown(:, :)
      integer :: iat, it, nfound, npts

      centroid = sum(mol%xyz, dim=2)/real(mol%nat, wp)
      call lsf%set_max_deriv(1)
      nfound = 0
      do iat = 1, min(nray, mol%nat)
         dir = mol%xyz(:, iat) - centroid
         if (norm2(dir) < 1.0e-8_wp) dir = [0.0_wp, 0.0_wp, 1.0_wp]
         dir = dir/norm2(dir)
         call lsf_value_at(lsf, mol%xyz(:, iat) + r_in*dir, s_lo, error)
         if (allocated(error)) return
         call lsf_value_at(lsf, mol%xyz(:, iat) + r_out*dir, s_hi, error)
         if (allocated(error)) return
         if (s_lo >= 0.0_wp .or. s_hi <= 0.0_wp) cycle
         lo = r_in
         hi = r_out
         do it = 1, nbisect
            mid = 0.5_wp*(lo + hi)
            call lsf_value_at(lsf, mol%xyz(:, iat) + mid*dir, s_mid, error)
            if (allocated(error)) return
            if (s_mid < 0.0_wp) then
               lo = mid
            else
               hi = mid
            end if
         end do
         call lsf_value_at(lsf, mol%xyz(:, iat) + lo*dir, s_mid, error)
         if (allocated(error)) return
         call check(error, abs(s_mid), 0.0_wp, thr=1.0e-10_wp*param%scale*param%rho_iso, &
                    more="bisection did not reach the isosurface")
         if (allocated(error)) return
         nfound = nfound + 1
         found(:, nfound) = mol%xyz(:, iat) + lo*dir
      end do
      call check(error, nfound > 0, "no atom ray crosses the isosurface")
      if (allocated(error)) return

      npts = size(pts, 2)
      allocate (grown(3, npts + nfound))
      grown(:, :npts) = pts
      grown(:, npts + 1:) = found(:, :nfound)
      call move_alloc(grown, pts)
   end subroutine append_surface_points

   !> Level-set value of an LSF at one point
   !>
   !> @param[inout] lsf    level set, prepared at `point`
   !> @param[in]    point  evaluation point in Bohr
   !> @param[out]   value  level-set value S
   !> @param[out]   error  set when `prepare` reports a failure
   subroutine lsf_value_at(lsf, point, value, error)
      !> Level set, prepared at `point`
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Level-set value S
      real(wp), intent(out) :: value
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: lsf_err
      real(wp) :: grad(3)

      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "LSF prepare failed: "//lsf_err%message)
         return
      end if
      call lsf%f012_r(value, grad)
   end subroutine lsf_value_at

   !> A non-positive screening threshold disables screening exactly
   !>
   !> With screening_threshold <= 0 the per-shell reach saturates at the
   !> evaluator's cap, so every shell contributes at every point (exact
   !> evaluation) and neighbor_cutoff degrades the cavity cell grid to a full
   !> scan. A tight positive threshold must reproduce that exact result at both
   !> near and far points, otherwise the screening bound is not conservative
   subroutine test_internal_screening_equivalence(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf_exact, lsf_screened
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: v_e, g_e(3), h_e(3, 3)
      real(wp) :: v_s, g_s(3), h_s(3, 3)
      real(wp) :: t_e(3, 3, 3), t_s(3, 3, 3)
      real(wp) :: far_pts(3, 4), centroid(3), extent
      type(mctc_error), allocatable :: lsf_err
      real(wp), parameter :: THR_SCREEN = 1.0e-12_wp
      real(wp), parameter :: SCREEN_TOL = 1.0e-14_wp
      !> Signed clearance of each far point beyond the molecular extent, in Bohr
      real(wp), parameter :: far_clearance(4) = &
                             [-20.0_wp, -5.0_wp, 5.0_wp, 15.0_wp]
      integer :: test, ip

      do test = 1, ntests
         call build_molecular_internal_lsf(lsf_exact, test, 0.0_wp, mol, error)
         if (allocated(error)) return
         call build_molecular_internal_lsf(lsf_screened, test, SCREEN_TOL, mol, error)
         if (allocated(error)) return

         !> Screening off -> reach saturates at the evaluator cap, so the cell-grid
         !> reach swallows the whole molecule and no shell is ever dropped
         call check(error, lsf_exact%gto%reach(0.0_wp) > 500.0_wp, &
                    "screening_threshold <= 0 must give an unbounded shell reach")
         if (allocated(error)) return
         call check(error, lsf_exact%neighbor_cutoff(2.0_wp) > 500.0_wp, &
                    "screening_threshold <= 0 must disable the cell-grid cutoff")
         if (allocated(error)) return

         !> Screening on -> a finite reach that actually drops shells
         call check(error, lsf_screened%gto%reach(SCREEN_TOL) < 50.0_wp, &
                    "positive screening threshold must give a finite shell reach")
         if (allocated(error)) return

         call lsf_exact%set_max_deriv(3)
         call lsf_screened%set_max_deriv(3)

         call get_test_points(mol, pts, 12)
         do ip = 1, size(pts, 2)
            call lsf_exact%prepare(pts(:, ip), lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "internal LSF prepare failed: "//lsf_err%message)
               return
            end if
            call lsf_exact%f3_rrr(v_e, g_e, h_e, t_e)
            call lsf_screened%prepare(pts(:, ip), lsf_err)
            call lsf_screened%f3_rrr(v_s, g_s, h_s, t_s)

            call check(error, abs(v_e - v_s), 0.0_wp, thr=THR_SCREEN, &
                       more="screened vs exact molecular isodensity value")
            if (allocated(error)) return
            call check(error, maxval(abs(g_e - g_s)), 0.0_wp, thr=THR_SCREEN, &
                       more="screened vs exact molecular isodensity gradient")
            if (allocated(error)) return
            call check(error, maxval(abs(h_e - h_s)), 0.0_wp, thr=THR_SCREEN, &
                       more="screened vs exact molecular isodensity Hessian")
            if (allocated(error)) return
            call check(error, maxval(abs(t_e - t_s)), 0.0_wp, thr=THR_SCREEN, &
                       more="screened vs exact molecular isodensity third derivative")
            if (allocated(error)) return
         end do

         ! Far points exercise outright shell drops beyond the diffuse reach
         ! They are placed relative to the molecule so every test, from a
         ! water molecule to a 16-atom structure, gets the same clearances
         centroid = sum(mol%xyz, dim=2)/real(mol%nat, wp)
         extent = maxval(norm2(mol%xyz - spread(centroid, 2, mol%nat), dim=1))
         do ip = 1, size(far_pts, 2)
            far_pts(:, ip) = centroid
            far_pts(3, ip) = centroid(3) &
                             + sign(extent + abs(far_clearance(ip)), far_clearance(ip))
         end do
         do ip = 1, size(far_pts, 2)
            call lsf_exact%prepare(far_pts(:, ip), lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "internal LSF prepare failed: "//lsf_err%message)
               return
            end if
            call lsf_exact%f3_rrr(v_e, g_e, h_e, t_e)
            call lsf_screened%prepare(far_pts(:, ip), lsf_err)
            call lsf_screened%f3_rrr(v_s, g_s, h_s, t_s)

            call check(error, abs(v_e - v_s), 0.0_wp, thr=THR_SCREEN, &
                       more="screened vs exact molecular isodensity value (far point)")
            if (allocated(error)) return
            call check(error, maxval(abs(g_e - g_s)), 0.0_wp, thr=THR_SCREEN, &
                       more="screened vs exact molecular isodensity gradient (far point)")
            if (allocated(error)) return
            call check(error, maxval(abs(h_e - h_s)), 0.0_wp, thr=THR_SCREEN, &
                       more="screened vs exact molecular isodensity Hessian (far point)")
            if (allocated(error)) return
            call check(error, maxval(abs(t_e - t_s)), 0.0_wp, thr=THR_SCREEN, &
                       more="screened vs exact molecular isodensity third derivative (far point)")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_internal_screening_equivalence

   !* ------------ *                  Fourth spatial derivative (internal) ------------ *!

   !> Prepare the internal LSF and translate an evaluator failure into a test one
   !>
   !> @param[inout] lsf   Internal isodensity LSF
   !> @param[in]    point Evaluation point in Bohr
   !> @param[out]   error Set when `prepare` reports a failure
   subroutine prepare_internal(lsf, point, error)
      !> Internal isodensity LSF
      type(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: lsf
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: lsf_err

      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) call test_failed(error, "internal LSF prepare failed: " &
                                               //lsf_err%message)
   end subroutine prepare_internal

   !> Third spatial derivative of the internal LSF at one point
   !>
   !> The caller must have set `max_deriv` to at least 3 already
   !>
   !> @param[inout] lsf   Internal isodensity LSF
   !> @param[in]    point Evaluation point in Bohr
   !> @param[out]   t3    Third spatial derivative
   !> @param[out]   error Set when `prepare` reports a failure
   subroutine third_at(lsf, point, t3, error)
      !> Internal isodensity LSF
      type(moist_cavity_drop_lsf_isodensity_internal_type), intent(inout) :: lsf
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Third spatial derivative
      real(wp), intent(out) :: t3(ndim, ndim, ndim)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call prepare_internal(lsf, point, error)
      if (allocated(error)) return
      call lsf%f3_rrr(lsf3_rrr=t3)
   end subroutine third_at

   !> The internal backend's analytic `f4_rrrr` matches a 4-point central FD of
   !> its own `f3_rrr`, and is fully index symmetric
   !>
   !> This is the LSF-level counterpart of [[test_gto_fourth_fd]]: it covers the
   !> sign/scale lift, the `(3,3,3,3)` cache and the per-instance scratch that the
   !> bare evaluator test does not touch
   subroutine test_internal_fourth_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: f4(ndim, ndim, ndim, ndim)
      real(wp) :: tpp(ndim, ndim, ndim), tp(ndim, ndim, ndim)
      real(wp) :: tm(ndim, ndim, ndim), tmm(ndim, ndim, ndim)
      real(wp) :: pp(ndim), fd, dev, scale_ref
      integer :: ip, ax, ix, jx, kx

      call build_molecular_internal_lsf(lsf, test_reference, 0.0_wp, mol, error)
      if (allocated(error)) return
      ! Non-unit scale, so an f4_rrrr that drops the scale factor cannot match the FD of
      ! f3_rrr; scale is only read at prepare, so setting it after update is safe
      lsf%param%scale = 2.3_wp
      call get_test_points(mol, pts, 8)

      scale_ref = 0.0_wp
      do ip = 1, size(pts, 2)
         call lsf%set_max_deriv(4)
         call prepare_internal(lsf, pts(:, ip), error)
         if (allocated(error)) return
         call check(error, lsf%prepared_deriv, 4, &
                    more="internal LSF cached the wrong order at max_deriv=4")
         if (allocated(error)) return
         call lsf%f4_rrrr(f4)
         scale_ref = max(scale_ref, maxval(abs(f4)))

         ! A fourth derivative of a smooth field is symmetric under every index
         ! permutation; the packed slot table must reproduce that exactly
         dev = 0.0_wp
         do ax = 1, ndim
            do ix = 1, ndim
               do jx = 1, ndim
                  do kx = 1, ndim
                     dev = max(dev, rel_deviation(f4(ax, ix, jx, kx), f4(ix, ax, jx, kx)))
                     dev = max(dev, rel_deviation(f4(ax, ix, jx, kx), f4(ax, jx, ix, kx)))
                     dev = max(dev, rel_deviation(f4(ax, ix, jx, kx), f4(ax, ix, kx, jx)))
                  end do
               end do
            end do
         end do
         call check(error, dev, 0.0_wp, thr=1.0e-14_wp, &
                    more="internal f4_rrrr permutation symmetry")
         if (allocated(error)) return

         ! Dropping back to order 3 for the FD stencil also exercises the cache
         ! being refilled at a lower order between two order-4 evaluations
         call lsf%set_max_deriv(3)
         do ax = 1, ndim
            pp = pts(:, ip); pp(ax) = pts(ax, ip) + 2*STEP_SIZE
            call third_at(lsf, pp, tpp, error)
            if (allocated(error)) return
            pp = pts(:, ip); pp(ax) = pts(ax, ip) + STEP_SIZE
            call third_at(lsf, pp, tp, error)
            if (allocated(error)) return
            pp = pts(:, ip); pp(ax) = pts(ax, ip) - STEP_SIZE
            call third_at(lsf, pp, tm, error)
            if (allocated(error)) return
            pp = pts(:, ip); pp(ax) = pts(ax, ip) - 2*STEP_SIZE
            call third_at(lsf, pp, tmm, error)
            if (allocated(error)) return

            do ix = 1, ndim
               do jx = 1, ndim
                  do kx = 1, ndim
                     call fd4_scalar(tpp(ix, jx, kx), tp(ix, jx, kx), tm(ix, jx, kx), &
                                     tmm(ix, jx, kx), STEP_SIZE, fd, error)
                     if (allocated(error)) return
                     call check(error, f4(ax, ix, jx, kx), fd, thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do

      !> Otherwise the agreement above would be a comparison of zeros
      call check(error, scale_ref > 1.0e-6_wp, &
                 "internal f4_rrrr vanished at every test point")
   end subroutine test_internal_fourth_fd

   !> Order 4 is reached only when asked for, and asking late is safe
   !>
   !> Three separate contracts:
   !>
   !>   * the per-instance scratch is sized for the order actually requested, and
   !>     raising `max_deriv` *after* a low-order `prepare` already sized it must
   !>     grow it rather than let `gto_eval` write past its end;
   !>   * the multi-column scratch (density panel, D phi_ij columns) stays
   !>     unallocated or narrow until a higher-order `prepare` needs it;
   !>   * raising the order must not disturb the lower orders' values beyond
   !>     round-off
   !>
   !> The callback backend is checked in the same test for the opposite property:
   !> its published C ABI stops at the third derivative, so `set_max_deriv(4)`
   !> must still leave `prepared_deriv` at 3 rather than claim an order it cannot
   !> deliver
   subroutine test_internal_fourth_gating(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: v3, g3(ndim), h3(ndim, ndim), t3(ndim, ndim, ndim)
      real(wp) :: v4, g4(ndim), h4(ndim, ndim), t4(ndim, ndim, ndim)
      real(wp) :: f4(ndim, ndim, ndim, ndim)
      integer :: ip
      logical :: any_nonzero

      call build_molecular_internal_lsf(lsf, test_reference, 0.0_wp, mol, error)
      if (allocated(error)) return
      call get_test_points(mol, pts, 6)
      call append_surface_points(lsf, lsf%param, mol, pts, error)
      if (allocated(error)) return

      ! The value+gradient projection phase sizes the scratch first ...
      call lsf%set_max_deriv(1)
      call prepare_internal(lsf, pts(:, 1), error)
      if (allocated(error)) return
      call check(error, lsf%prepared_deriv, 1, &
                 more="internal LSF cached the wrong order at max_deriv=1")
      if (allocated(error)) return
      call check(error, size(lsf%work%phi, 2), moist_iso_gto_nslot(1), &
                 more="max_deriv=1 must not size the AO table for a higher order")
      if (allocated(error)) return
      call check(error, size(lsf%work%dphi, 2), moist_iso_gto_nrhs(1), &
                 more="max_deriv=1 must not size the density-weighted columns for a higher order")
      if (allocated(error)) return
      call check(error, .not. allocated(lsf%work%panel), &
                 "the multi-column density panel must not be allocated at max_deriv=1")
      if (allocated(error)) return

      any_nonzero = .false.
      do ip = 1, size(pts, 2)
         call lsf%set_max_deriv(3)
         call prepare_internal(lsf, pts(:, ip), error)
         if (allocated(error)) return
         call lsf%f3_rrr(v3, g3, h3, t3)

         ! ... and only then is order 4 asked for, on an already-sized instance
         call lsf%set_max_deriv(4)
         call prepare_internal(lsf, pts(:, ip), error)
         if (allocated(error)) return
         call check(error, lsf%prepared_deriv, 4, &
                    more="internal LSF cached the wrong order at max_deriv=4")
         if (allocated(error)) return
         call check(error, size(lsf%work%phi, 2), moist_iso_gto_nslot(4), &
                    more="max_deriv=4 must grow the AO table to the fourth-order slots")
         if (allocated(error)) return
         call check(error, size(lsf%work%dphi, 2), moist_iso_gto_nrhs(4), &
                    more="max_deriv=4 must grow the density-weighted columns to D phi_ij")
         if (allocated(error)) return
         call lsf%f3_rrr(v4, g4, h4, t4)
         call lsf%f4_rrrr(f4)
         if (maxval(abs(f4)) > 0.0_wp) any_nonzero = .true.

         call check(error, v3, v4, thr=ORDER_THR*value_magnitude(lsf%param, v4), &
                    more="internal isodensity LSF value depends on max_deriv")
         if (allocated(error)) return
         call check(error, maxval(abs(g3 - g4)), 0.0_wp, thr=ORDER_THR*maxval(abs(g4)), &
                    more="internal isodensity LSF gradient depends on max_deriv")
         if (allocated(error)) return
         call check(error, maxval(abs(h3 - h4)), 0.0_wp, thr=ORDER_THR*maxval(abs(h4)), &
                    more="internal isodensity LSF Hessian depends on max_deriv")
         if (allocated(error)) return
         call check(error, maxval(abs(t3 - t4)), 0.0_wp, thr=ORDER_THR*maxval(abs(t4)), &
                    more="internal isodensity LSF third derivative depends on max_deriv")
         if (allocated(error)) return
      end do

      call check(error, any_nonzero, &
                 "internal isodensity LSF returned a zero f4_rrrr at max_deriv=4")
      if (allocated(error)) return

      ! The callback ABI stops at the third derivative and must say so
      !$omp critical(iso_reference_density)
      call run_callback_fourth_cap(error)
      !$omp end critical(iso_reference_density)
   end subroutine test_internal_fourth_gating

   !> Body of the callback half of [[test_internal_fourth_gating]], run under the
   !> `cb_gto` lock
   subroutine run_callback_fourth_cap(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      type(mctc_error), allocatable :: lsf_err
      integer :: ip

      call set_reference_density(mol, error)
      if (allocated(error)) return
      call build_callback_lsf(lsf, mol)
      call get_test_points(mol, pts, 4)

      call lsf%set_max_deriv(4)
      do ip = 1, size(pts, 2)
         call lsf%prepare(pts(:, ip), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "callback LSF prepare failed: "//lsf_err%message)
            return
         end if
         !> Not 4: the C ABI carries no fourth-derivative buffer, so claiming
         !> order 4 here would make `require_deriv` wave through an `f4_rrrr`
         !> read that the base type is supposed to abort
         call check(error, lsf%prepared_deriv, 3, &
                    more="callback isodensity LSF claimed an order its ABI cannot deliver")
         if (allocated(error)) return
      end do
   end subroutine run_callback_fourth_cap

end module test_cavity_drop_isodensity
