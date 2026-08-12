!> Unit tests for the isodensity DROP level set functions and their internal
!> cartesian-monomial GTO density evaluator
!>
!> Three parts are covered:
!>
!>   * the bare evaluator [[moist_iso_gto_type]]: the assembled density against an
!>     independent direct monomial evaluation, and each analytic derivative order
!>     against a 4-point central FD of the analytic previous order; plus the
!>     radial screening, which must not perturb any evaluated quantity beyond its
!>     own threshold.
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
!> MB16-43/01, and Heavy28/h2o in STO-3G, def2-SVP, and def2-TZVP.
!> Geoms from mstore, the bases and converged density matrices from the reference
!> files in `test/unit/data` (generated using PySCF)
module test_cavity_drop_isodensity
   use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_null_ptr, c_funloc, &
                                          c_associated, c_f_pointer
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use test_helpers, only: fd4_scalar, get_test_points, center_at_origin, rel_deviation
   use moist_utils_env, only: get_env
   use moist_model_gems_utils, only: BuildSuperStructure
   use moist_cavity_drop_lsf_isodensity_gto, only: moist_iso_gto_type, moist_iso_gto_ncart
   use moist_cavity_drop_lsf_isodensity_internal, only: &
      moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_drop_lsf_isodensity_callback, only: &
      moist_cavity_drop_lsf_isodensity_callback_type
   use testdrive, only: new_unittest, unittest_type, check, test_failed, error_type
   implicit none(type, external)
   private

   public :: collect_cavity_drop_isodensity

   integer, parameter :: ndim = 3

   !> Every test pairs one mstore record with one basis set
   integer, parameter :: nmolecules = 3
   integer, parameter :: nbases = 3
   integer, parameter :: ntests = 8

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
                          imol_heavy28, imol_heavy28, imol_heavy28]
   !> Basis set of each test
   integer, parameter :: test_basis(ntests) = &
                         [ibasis_sto_3g, ibasis_def2_svp, ibasis_def2_tzvp, &
                          ibasis_sto_3g, ibasis_def2_svp, &
                          ibasis_sto_3g, ibasis_def2_svp, ibasis_def2_tzvp]

   !> test used wherever one reference density is enough. Heavy28/h2o in
   !> def2-TZVP is the cheapest test that still carries f functions.
   integer, parameter :: test_reference = 8

   !> Shell and Cartesian-function counts of the checked-in reference files.
   !> Asserted so a regenerated file cannot silently change basis set.
   integer, parameter :: expected_nshell(ntests) = &
                         [28, 66, 106, 42, 83, 5, 12, 19]
   integer, parameter :: expected_ncart(ntests) = &
                         [40, 140, 276, 68, 189, 7, 25, 48]
   !> Converged RHF energies of the checked-in reference files, in Hartree
   real(wp), parameter :: expected_energy(ntests) = [ &
                          -3.03120634898054163e+02_wp, -3.06792977091366652e+02_wp, &
                          -3.07129116858719954e+02_wp, -1.25920490619214843e+03_wp, &
                          -1.27392053922109562e+03_wp, -7.49631306328200253e+01_wp, &
                          -7.59622034508552559e+01_wp, -7.60593754022876709e+01_wp]

   !> Density isovalue defining the reference surface (Bohr^-3)
   real(wp), parameter :: rho_iso_ref = 4.0e-4_wp

   !> Level set multiplier used by both backends in the agreement tests
   real(wp), parameter :: lsf_scale = 1.0_wp/rho_iso_ref

   !> Reference density shared with the C-interoperable callback
   type(moist_iso_gto_type), save :: cb_gto

   !> Packed spatial derivative multi-indices through third order.
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
                  new_unittest("internal_vs_callback", test_internal_vs_callback), &
                  new_unittest("max_deriv_gating_internal", test_max_deriv_internal), &
                  new_unittest("max_deriv_gating_callback", test_max_deriv_callback), &
                  new_unittest("internal_screening_equivalence", test_internal_screening_equivalence) &
                  ]
   end subroutine collect_cavity_drop_isodensity

   !> Build one realistic molecular basis test and bind its geometry
   !>
   !> @param[out] gto     Configured basis + density
   !> @param[in]  test test identifier
   !> @param[out] error   Set if the evaluator rejects the test
   !> @param[out] mol_out Optional copy of the test molecule
   subroutine build_test(gto, test, error, mol_out)
      !> Configured basis + density
      type(moist_iso_gto_type), intent(out) :: gto
      !> test identifier
      integer, intent(in) :: test
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Optional copy of the test molecule
      type(structure_type), optional, intent(out) :: mol_out

      type(mctc_error), allocatable :: merr
      type(structure_type) :: mol
      integer, allocatable :: sh_atom(:), sh_l(:), sh_nprim(:)
      real(wp), allocatable :: exps(:), coeffs(:)
      real(wp), allocatable :: dcart(:, :)

      call test_molecule(test, mol)
      call load_test_reference(test, mol, sh_atom, sh_l, sh_nprim, exps, coeffs, &
                               dcart, error)
      if (allocated(error)) return

      call gto%init(sh_atom, sh_l, sh_nprim, exps, coeffs, merr)
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
   end subroutine build_test

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

   !> Fetch the mstore record a test is built on.
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

   !> Check the Cartesian component order assumed by the PySCF density files.
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

      call check(error, trim(header), "PySCF 2.14.0 RHF gas-phase Cartesian density")
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
      !> the tests evaluate on, otherwise every reference below is meaningless.
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
      ! the source root, fpm runs the tester from the project root.
      source_root = get_env("MOIST_SOURCE_ROOT", default=".")
      filename = source_root//"/test/unit/data/" &
                 //trim(mol_tag(test_molecule_index(test)))//"_" &
                 //trim(basis_tag(test_basis_index(test)))//".density"
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

      real(wp), allocatable :: phi(:, :), t0(:), tm(:, :), tmm(:, :)
      integer, allocatable :: act(:)

      ! The fourth order packs 35 derivative slots and needs the extra
      ! density-weighted Hessian scratch; the lower orders fit in 20
      if (nderiv >= 4) then
         allocate (phi(gto%ncart, 0:34), tmm(gto%ncart, 6))
      else
         allocate (phi(gto%ncart, 0:19))
      end if
      allocate (t0(gto%ncart), tm(gto%ncart, 3), act(gto%ncart))

      ! The evaluator computes exactly the orders it is asked to return, so the
      ! orders this helper does not request are zeroed here for its own callers
      d2rho = 0.0_wp
      d3rho = 0.0_wp
      if (present(d4rho)) d4rho = 0.0_wp
      select case (nderiv)
      case (:1)
         call gto%eval(point, phi, t0, tm, act, rho, drho)
      case (2)
         call gto%eval(point, phi, t0, tm, act, rho, drho, d2rho=d2rho)
      case (3)
         call gto%eval(point, phi, t0, tm, act, rho, drho, d2rho=d2rho, d3rho=d3rho)
      case default
         if (.not. present(d4rho)) error stop "eval_at: nderiv 4 requires d4rho"
         call gto%eval(point, phi, t0, tm, act, rho, drho, d2rho=d2rho, d3rho=d3rho, &
                       d4rho=d4rho, tmm=tmm)
      end select
   end subroutine eval_at

   !> Independent direct density: rho = sum_cc' D_cc' g_c g_c'
   !>
   !> @param[in] gto   Basis + density
   !> @param[in] point Evaluation point in Bohr
   !> @returns         Density at the point
   function rho_reference(gto, point) result(rho)
      !> Basis + density
      type(moist_iso_gto_type), intent(in) :: gto
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
            rho = rho + gto%dcart(i, j)*g(i)*g(j)
         end do
      end do
   end function rho_reference

   !> Analytic density value matches the independent direct evaluation
   !>
   !> @param[out] error Set on mismatch
   subroutine test_gto_value_reference(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rho, drho(3), d2rho(3, 3), d3rho(3, 3, 3)
      integer :: test, ip

      do test = 1, ntests
         call build_test(gto, test, error, mol)
         if (allocated(error)) return
         call get_test_points(mol, pts, 12)
         do ip = 1, size(pts, 2)
            call eval_at(gto, pts(:, ip), 1, rho, drho, d2rho, d3rho)
            call check(error, rho, rho_reference(gto, pts(:, ip)), thr=1.0e-12_wp)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_gto_value_reference

   !> Analytic gradient matches a 4-point central FD of the density
   !>
   !> @param[out] error Set on mismatch
   subroutine test_gto_grad_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rho, drho(3), d2(3, 3), d3(3, 3, 3)
      real(wp) :: dg(3), dd2(3, 3), dd3(3, 3, 3)
      real(wp) :: pp(3), rpp, rp, rm, rmm, fd
      real(wp), parameter :: h = 1.0e-3_wp
      integer :: test, ip, ax

      do test = 1, ntests
         call build_test(gto, test, error, mol)
         if (allocated(error)) return
         call get_test_points(mol, pts, 12)
         do ip = 1, size(pts, 2)
            call eval_at(gto, pts(:, ip), 1, rho, drho, d2, d3)
            do ax = 1, ndim
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + 2*h
               call eval_at(gto, pp, 1, rpp, dg, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + h
               call eval_at(gto, pp, 1, rp, dg, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - h
               call eval_at(gto, pp, 1, rm, dg, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - 2*h
               call eval_at(gto, pp, 1, rmm, dg, dd2, dd3)
               fd = fd4_scalar(rpp, rp, rm, rmm, h)
               call check(error, drho(ax), fd, thr=1.0e-8_wp)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_gto_grad_fd

   !> Analytic Hessian matches a 4-point central FD of the gradient
   !>
   !> @param[out] error Set on mismatch
   subroutine test_gto_hess_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rho, drho(3), d2(3, 3), d3(3, 3, 3)
      real(wp) :: drr, dd2(3, 3), dd3(3, 3, 3)
      real(wp) :: pp(3), gpp(3), gp(3), gm(3), gmm(3), fd
      real(wp), parameter :: h = 1.0e-3_wp
      integer :: test, ip, ax, jx

      do test = 1, ntests
         call build_test(gto, test, error, mol)
         if (allocated(error)) return
         call get_test_points(mol, pts, 12)
         do ip = 1, size(pts, 2)
            call eval_at(gto, pts(:, ip), 2, rho, drho, d2, d3)
            do ax = 1, ndim
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + 2*h
               call eval_at(gto, pp, 1, drr, gpp, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + h
               call eval_at(gto, pp, 1, drr, gp, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - h
               call eval_at(gto, pp, 1, drr, gm, dd2, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - 2*h
               call eval_at(gto, pp, 1, drr, gmm, dd2, dd3)
               do jx = 1, ndim
                  fd = fd4_scalar(gpp(jx), gp(jx), gm(jx), gmm(jx), h)
                  call check(error, d2(ax, jx), fd, thr=1.0e-7_wp)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_gto_hess_fd

   !> Analytic third derivative matches a 4-point central FD of the Hessian
   !>
   !> @param[out] error Set on mismatch
   subroutine test_gto_third_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: gto
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rho, drho(3), d2(3, 3), d3(3, 3, 3)
      real(wp) :: drr, dg(3), dd3(3, 3, 3)
      real(wp) :: pp(3), hpp(3, 3), hp(3, 3), hm(3, 3), hmm(3, 3), fd
      real(wp), parameter :: h = 2.0e-3_wp
      integer :: test, ip, ax, jx, kx

      do test = 1, ntests
         call build_test(gto, test, error, mol)
         if (allocated(error)) return
         call get_test_points(mol, pts, 12)
         do ip = 1, size(pts, 2)
            call eval_at(gto, pts(:, ip), 3, rho, drho, d2, d3)
            do ax = 1, ndim
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + 2*h
               call eval_at(gto, pp, 2, drr, dg, hpp, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) + h
               call eval_at(gto, pp, 2, drr, dg, hp, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - h
               call eval_at(gto, pp, 2, drr, dg, hm, dd3)
               pp = pts(:, ip); pp(ax) = pts(ax, ip) - 2*h
               call eval_at(gto, pp, 2, drr, dg, hmm, dd3)
               do jx = 1, ndim
                  do kx = 1, ndim
                     fd = fd4_scalar(hpp(jx, kx), hp(jx, kx), hm(jx, kx), hmm(jx, kx), h)
                     call check(error, d3(ax, jx, kx), fd, thr=1.0e-6_wp)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_gto_third_fd

   !> Analytic fourth derivative matches a 4-point central FD of the third
   !>
   !> @param[out] error Set on mismatch
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
      real(wp), parameter :: h = 2.0e-3_wp
      integer :: ip, ax, ix, jx, kx

      call build_test(gto, test_reference, error, mol)
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
            pp = pts(:, ip); pp(ax) = pts(ax, ip) + 2*h
            call eval_at(gto, pp, 3, drr, dg, dh, dpp)
            pp = pts(:, ip); pp(ax) = pts(ax, ip) + h
            call eval_at(gto, pp, 3, drr, dg, dh, dp)
            pp = pts(:, ip); pp(ax) = pts(ax, ip) - h
            call eval_at(gto, pp, 3, drr, dg, dh, dm)
            pp = pts(:, ip); pp(ax) = pts(ax, ip) - 2*h
            call eval_at(gto, pp, 3, drr, dg, dh, dmm)
            do ix = 1, ndim
               do jx = 1, ndim
                  do kx = 1, ndim
                     fd = fd4_scalar(dpp(ix, jx, kx), dp(ix, jx, kx), &
                                     dm(ix, jx, kx), dmm(ix, jx, kx), h)
                     call check(error, d4(ax, ix, jx, kx), fd, thr=1.0e-5_wp)
                     if (allocated(error)) return
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
   !> error bound throughout the transition.
   subroutine test_gto_screening_separation(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_iso_gto_type) :: base_gto, pair_gto, exact_gto, screened_gto
      type(structure_type) :: center_mol, pair_mol
      real(wp), allocatable :: points(:, :)
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
      !> Bohr; the test points themselves already sit up to `pad` outside the box.
      real(wp), parameter :: gap_margin = 10.0_wp

      do test = 1, ntests
         call build_test(base_gto, test, error, center_mol)
         if (allocated(error)) return
         call center_at_origin(center_mol)
         call base_gto%refresh_centers(center_mol)
         call get_test_points(center_mol, points, 12)

         !> The scan has to end beyond the widest shell reach of this test,
         !> which spans a factor of four between STO-3G carbon and the diffuse
         !> def2-TZVP sodium shells of MB16-43/01.
         gap_step = (base_gto%reach(minval(thresholds)) + gap_margin)/real(nsteps, wp)

         do ithreshold = 1, size(thresholds)
            threshold = thresholds(ithreshold)
            adjacent_active = .false.
            adjacent_effect = .false.
            far_inactive = .true.

            do istep = 0, nsteps
               gap = real(istep, wp)*gap_step
               call build_density_pair(base_gto, center_mol, gap, pair_gto, pair_mol, error)
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

   !> Build two translated molecular copies with a block-diagonal density
   subroutine build_density_pair(base_gto, center_mol, gap, pair_gto, pair_mol, error)
      type(moist_iso_gto_type), intent(in) :: base_gto
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
      call pair_gto%init(sh_atom, sh_l, sh_nprim, exps, coeffs, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call pair_gto%refresh_centers(pair_mol)

      n = base_gto%ncart
      allocate (density(2*n, 2*n), source=0.0_wp)
      density(1:n, 1:n) = base_gto%dcart
      density(n + 1:, n + 1:) = base_gto%dcart
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

      call build_test(cb_gto, test_reference, error, mol)
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
                     value = value + real(weight, wp)*cb_gto%dcart(i, j) &
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

   !> C-interoperable isodensity level set callback
   !>
   !> Returns the *unscaled* level set S = rho_iso - rho and its spatial
   !> derivatives
   !>
   !> The LSF applies its own multiplier
   !>
   !> hess/third are NULL when the cavity did not request that order
   !>
   !> @param[in]  context Unused callback context
   !> @param[in]  point   Evaluation point in Bohr
   !> @param[out] value   Level set value
   !> @param[out] grad    Level set spatial gradient
   !> @param[out] hess    Level set spatial Hessian (Fortran (3,3), or NULL)
   !> @param[out] third   Level set third spatial derivative (Fortran (3,3,3), or NULL)
   !> @returns            0 (this reference implementation never fails)
   function iso_reference_callback(context, point, value, grad, hess, third) result(status) bind(C)
      !> Unused callback context (the test keeps its data in module state)
      type(c_ptr), value :: context
      !> Evaluation point in Bohr
      real(c_double), intent(in) :: point(3)
      !> Level set value
      real(c_double), intent(out) :: value
      !> Level set spatial gradient
      real(c_double), intent(out) :: grad(3)
      !> Level set spatial Hessian (Fortran (3,3), or NULL)
      type(c_ptr), value :: hess
      !> Level set third spatial derivative (Fortran (3,3,3), or NULL)
      type(c_ptr), value :: third
      !> Zero on success
      integer(c_int) :: status

      real(c_double), pointer :: hptr(:, :), tptr(:, :, :)
      real(wp) :: rho, drho(3), d2rho(3, 3), d3rho(3, 3, 3)
      logical :: want_hess, want_third

      status = 0_c_int
      if (c_associated(context)) return
      want_hess = c_associated(hess)
      want_third = c_associated(third)

      call rho_product_rule(real(point, wp), rho, drho, d2rho, d3rho)

      value = real(rho_iso_ref - rho, c_double)
      grad = real(-drho, c_double)
      if (want_hess) then
         call c_f_pointer(hess, hptr, [3, 3])
         hptr = real(-d2rho, c_double)
      end if
      if (want_third) then
         call c_f_pointer(third, tptr, [3, 3, 3])
         tptr = real(-d3rho, c_double)
      end if
   end function iso_reference_callback

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

      call lsf%new(cb_gto%sh_atom, cb_gto%sh_l, sh_nprim, cb_gto%exps, cb_gto%coeffs, rho_iso_ref, &
                   lsf_scale, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call lsf%set_density(cb_gto%dcart, merr)
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
      real(wp), allocatable :: radii(:)

      call build_test(gto, test, error, mol)
      if (allocated(error)) return
      sh_nprim = gto%sh_poff(2:) - gto%sh_poff(:gto%nshell)
      allocate (radii(mol%nat), source=2.0_wp)

      call lsf%new(gto%sh_atom, gto%sh_l, sh_nprim, gto%exps, gto%coeffs, &
                   rho_iso_ref, 1.0_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call lsf%set_density(gto%dcart, merr)
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
      call lsf%new(c_funloc(iso_reference_callback), c_null_ptr, lsf_scale)
      call lsf%update(mol, radii)
   end subroutine build_callback_lsf

   !> Both isodensity backends must describe the same level set
   !>
   !> @param[out] error Set on mismatch
   subroutine test_internal_vs_callback(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !$omp critical(iso_reference_density)
      call run_internal_vs_callback(error)
      !$omp end critical(iso_reference_density)
   end subroutine test_internal_vs_callback

   !> Body of [[test_internal_vs_callback]], run under the `cb_gto` lock
   !>
   !> @param[out] error Set on mismatch
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
      real(wp), allocatable :: t_i(:, :, :), t_c(:, :, :)
      real(wp) :: dev_val, dev_grad, dev_hess, dev_third, scale_ref
      type(mctc_error), allocatable :: lsf_err
      integer :: ip

      !> Roundoff-limited thresholds include headroom for the different summation
      !> orders of the dense production contraction and direct Leibniz reference.
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
         call lsf_int%f0_screened(v0_i)
         call lsf_int%f012_r_screened(v01_i, g01_i)
         call lsf_int%f012_r_screened(v012_i, g012_i, h012_i)
         call lsf_int%f3_rrr_screened(v_i, g_i, h_i, t_i)
         call lsf_cb%prepare(pts(:, ip), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "callback LSF prepare failed: "//lsf_err%message)
            return
         end if
         call lsf_cb%f0_screened(v0_c)
         call lsf_cb%f012_r_screened(v01_c, g01_c)
         call lsf_cb%f012_r_screened(v012_c, g012_c, h012_c)
         call lsf_cb%f3_rrr_screened(v_c, g_c, h_c, t_c)

         call check(error, v0_i, v0_c, thr=THR_VAL, &
                    more="internal vs callback f0_screened value")
         if (allocated(error)) return
         call check(error, v01_i, v01_c, thr=THR_VAL, &
                    more="internal vs callback f012_r_screened value without Hessian")
         if (allocated(error)) return
         call check(error, maxval(abs(g01_i - g01_c)), 0.0_wp, thr=THR_GRAD, &
                    more="internal vs callback f012_r_screened gradient without Hessian")
         if (allocated(error)) return
         call check(error, v012_i, v012_c, thr=THR_VAL, &
                    more="internal vs callback f012_r_screened value")
         if (allocated(error)) return
         call check(error, maxval(abs(g012_i - g012_c)), 0.0_wp, thr=THR_GRAD, &
                    more="internal vs callback f012_r_screened gradient")
         if (allocated(error)) return
         call check(error, maxval(abs(h012_i - h012_c)), 0.0_wp, thr=THR_HESS, &
                    more="internal vs callback f012_r_screened Hessian")
         if (allocated(error)) return

         dev_val = max(dev_val, abs(v_i - v_c))
         dev_grad = max(dev_grad, maxval(abs(g_i - g_c)))
         dev_hess = max(dev_hess, maxval(abs(h_i - h_c)))
         dev_third = max(dev_third, maxval(abs(t_i - t_c)))
         scale_ref = max(scale_ref, maxval(abs(t_i)))
      end do

      !> The test must actually produce a non-trivial level set; otherwise the
      !> agreement above would be vacuous.
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
   !>
   !> @param[out] error Set on contract violation
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

         any_nonzero = .false.
         do ip = 1, size(pts, 2)
            !> Value+gradient only: the Hessian is deliberately not requested.
            call lsf%set_max_deriv(1)
            call lsf%prepare(pts(:, ip), lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "internal LSF prepare failed: "//lsf_err%message)
               return
            end if
            call check(error, lsf%prepared_deriv, 1, &
                       more="internal molecular test cached the wrong order at max_deriv=1")
            if (allocated(error)) return
            call lsf%f012_r_screened(v1, g1)

            call lsf%set_max_deriv(2)
            call lsf%prepare(pts(:, ip), lsf_err)
            call check(error, lsf%prepared_deriv, 2, &
                       more="internal molecular test cached the wrong order at max_deriv=2")
            if (allocated(error)) return
            call lsf%f012_r_screened(v2, g2, h2)

            call check_max_deriv_pair(error, v1, g1, v2, g2, h2, &
                                      "internal molecular test", any_nonzero)
            if (allocated(error)) return
         end do

         call check(error, any_nonzero, &
                    "internal molecular test returned a zero Hessian at max_deriv=2")
         if (allocated(error)) return
      end do
   end subroutine test_max_deriv_internal

   !> max_deriv must gate what the callback backend requests through the ABI
   !>
   !> @param[out] error Set on contract violation
   subroutine test_max_deriv_callback(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !$omp critical(iso_reference_density)
      call run_max_deriv_callback(error)
      !$omp end critical(iso_reference_density)
   end subroutine test_max_deriv_callback

   !> Body of [[test_max_deriv_callback]], run under the `cb_gto` lock
   !>
   !> @param[out] error Set on contract violation
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

      any_nonzero = .false.
      do ip = 1, size(pts, 2)
         !> Value+gradient only: the Hessian is deliberately not requested.
         call lsf%set_max_deriv(1)
         call lsf%prepare(pts(:, ip), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "callback LSF prepare failed: "//lsf_err%message)
            return
         end if
         call check(error, lsf%prepared_deriv, 1, &
                    more="callback isodensity LSF cached the wrong order at max_deriv=1")
         if (allocated(error)) return
         call lsf%f012_r_screened(v1, g1)

         call lsf%set_max_deriv(2)
         call lsf%prepare(pts(:, ip), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "callback LSF prepare failed: "//lsf_err%message)
            return
         end if
         call check(error, lsf%prepared_deriv, 2, &
                    more="callback isodensity LSF cached the wrong order at max_deriv=2")
         if (allocated(error)) return
         call lsf%f012_r_screened(v2, g2, h2)

         call check_max_deriv_pair(error, v1, g1, v2, g2, h2, "callback", any_nonzero)
         if (allocated(error)) return
      end do

      call check(error, any_nonzero, &
                 "callback isodensity LSF returned a zero Hessian at max_deriv=2")
   end subroutine run_max_deriv_callback

   !> Shared assertions for one max_deriv=1 / max_deriv=2 evaluation pair
   !>
   !> @param[out]   error       Set on contract violation
   !> @param[in]    v1          Value at max_deriv = 1
   !> @param[in]    g1          Gradient at max_deriv = 1
   !> @param[in]    v2          Value at max_deriv = 2
   !> @param[in]    g2          Gradient at max_deriv = 2
   !> @param[in]    h2          Hessian at max_deriv = 2
   !> @param[in]    label       Backend name used in failure messages
   !> @param[inout] any_nonzero Set once a non-zero max_deriv=2 Hessian is seen
   subroutine check_max_deriv_pair(error, v1, g1, v2, g2, h2, label, any_nonzero)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
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

      !> Lowering max_deriv must not disturb the orders that are still requested.
      call check(error, v1, v2, thr=0.0_wp, &
                 more=label//" isodensity LSF value depends on max_deriv")
      if (allocated(error)) return
      call check(error, maxval(abs(g1 - g2)), 0.0_wp, thr=0.0_wp, &
                 more=label//" isodensity LSF gradient depends on max_deriv")
      if (allocated(error)) return

      if (maxval(abs(h2)) > 0.0_wp) any_nonzero = .true.
   end subroutine check_max_deriv_pair

   !> A non-positive screening threshold disables screening exactly
   !>
   !> With screening_threshold <= 0 the per-shell reach saturates at the
   !> evaluator's cap, so every shell contributes at every point (exact
   !> evaluation) and neighbor_cutoff degrades the cavity cell grid to a full
   !> scan. A tight positive threshold must reproduce that exact result at both
   !> near and far points, otherwise the screening bound is not conservative.
   !>
   !> @param[out] error Set on mismatch
   subroutine test_internal_screening_equivalence(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf_exact, lsf_screened
      type(structure_type) :: mol
      real(wp), allocatable :: pts(:, :)
      real(wp) :: v_e, g_e(3), h_e(3, 3)
      real(wp) :: v_s, g_s(3), h_s(3, 3)
      real(wp), allocatable :: t_e(:, :, :), t_s(:, :, :)
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
         !> reach swallows the whole molecule and no shell is ever dropped.
         call check(error, lsf_exact%gto%reach(0.0_wp) > 500.0_wp, &
                    "screening_threshold <= 0 must give an unbounded shell reach")
         if (allocated(error)) return
         call check(error, lsf_exact%neighbor_cutoff(2.0_wp) > 500.0_wp, &
                    "screening_threshold <= 0 must disable the cell-grid cutoff")
         if (allocated(error)) return

         !> Screening on -> a finite reach that actually drops shells.
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
            call lsf_exact%f3_rrr_screened(v_e, g_e, h_e, t_e)
            call lsf_screened%prepare(pts(:, ip), lsf_err)
            call lsf_screened%f3_rrr_screened(v_s, g_s, h_s, t_s)

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

         ! Far points exercise outright shell drops beyond the diffuse reach.
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
            call lsf_exact%f3_rrr_screened(v_e, g_e, h_e, t_e)
            call lsf_screened%prepare(far_pts(:, ip), lsf_err)
            call lsf_screened%f3_rrr_screened(v_s, g_s, h_s, t_s)

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

end module test_cavity_drop_isodensity
