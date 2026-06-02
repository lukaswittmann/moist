!> Dev test for fitting DROP Gaussian widths to analytical Born energies.
!>
!> The geometry used here is a *vestigial structural carrier*: the
!> physically meaningful object is a single sphere at the origin with
!> a swept radius. To honour the project convention that test molecules
!> come from mstore, we load MB16-43/H2 and park its second atom 100 bohr
!> away with a 1e-3 bohr radius. That atom's surface area (~1e-5 bohr^2)
!> and Born contribution are O(1e-7) - well below the test's
!> max_rel_energy_error = 1e-4 tolerance.
module test_cavity_drop_born_fit
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env_accuracy, only : wp
   use mctc_env_error, only : mctc_error => error_type
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use moist_cavity_drop, only : cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only : moist_cavity_drop_lsf_svdw_type
   use moist_math_grid_lebedev, only : get_angular_grid, grid_size, &
      lebedev_order_from_num
   use moist_model_component_pcm_solvers, only : solve_pcm_cholesky
   use moist_radii, only : new_radii_custom_atoms, radius_type
   use moist_utils_prettylistprint, only : prettylistprinter, new_prettylistprinter
   use testdrive, only : new_unittest, unittest_type, error_type, test_failed
   implicit none
   private

   public :: collect_cavity_drop_born_fit

   !> Dielectric constant used for the CPCM/Born comparison.
   real(wp), parameter :: epsilon = 100.0_wp
   !> CPCM dielectric prefactor.
   real(wp), parameter :: feps = (epsilon - 1.0_wp) / epsilon
   !> Central point charge used for the spherical Born reference.
   real(wp), parameter :: source_charge = 1.0_wp
   !> Radii used to verify the fitted zeta reproduces 1/R Born scaling.
   real(wp), parameter :: sphere_radii(*) = [1.0_wp, 2.0_wp, 3.0_wp, 5.0_wp, 7.0_wp]
   !> Fitting lower bound for zeta.
   real(wp), parameter :: zeta_lower = 1.0_wp
   !> Fitting upper bound for zeta.
   real(wp), parameter :: zeta_upper = 10.0_wp
   !> Maximum tolerated relative energy error after fitting.
   real(wp), parameter :: max_rel_energy_error = 1.0e-4_wp
   !> Maximum tolerated difference to the constructor-selected zeta.
   real(wp), parameter :: max_zeta_param_error = 2.0e-4_wp
   !> Stop threshold for the golden-section fit interval.
   real(wp), parameter :: fit_interval_tol = 1.0e-4_wp
   !> Maximum golden-section iterations.
   integer, parameter :: fit_maxiter = 40

contains

   !> Collect all DROP Born fit tests.
   subroutine collect_cavity_drop_born_fit(testsuite)
      !> Collection of tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("all_lebedev_grids", test_all_lebedev_grids) &
         ]
   end subroutine collect_cavity_drop_born_fit

   !> Fit zeta_Born for all supported Lebedev grids.
   subroutine test_all_lebedev_grids(error)
      !> Test-drive error object.
      type(error_type), allocatable, intent(out) :: error

      integer :: igrid
      integer :: ntested
      type(prettylistprinter) :: plp

      ntested = 0
      plp = new_prettylistprinter([8, 12, 14, 14, 42], &
         [character(len=14) :: "Lebedev", "Status", "zeta", "max rel err", "Note"], &
         unit=error_unit, offset=4, column_gap=2)

      call plp%blank()
      call plp%header("DROP Born zeta fit")
      call plp%blank()
      call plp%print_header()
      call plp%separator()

      do igrid = 1, size(grid_size)
         if (.not. drop_supports_lebedev_grid(grid_size(igrid))) then
            call plp%begin_row()
            call plp%add(grid_size(igrid))
            call plp%add("skipped")
            call plp%add("-")
            call plp%add("-")
            call plp%add("negative weights unsupported by DROP")
            call plp%end_row()
            cycle
         end if

         call check_lebedev_grid(grid_size(igrid), plp, error)
         if (allocated(error)) return
         ntested = ntested + 1
      end do

      call plp%separator()
      call plp%blank()

      if (ntested == 0) then
         call test_failed(error, "No positive-weight Lebedev grids were tested")
      end if
   end subroutine test_all_lebedev_grids

   !> Fit and validate one Lebedev grid size.
   subroutine check_lebedev_grid(nleb, plp, error)
      !> Number of Lebedev points per sphere.
      integer, intent(in) :: nleb
      !> Results table printer.
      type(prettylistprinter), intent(inout) :: plp
      !> Test-drive error object.
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop) :: cavities(size(sphere_radii))
      real(wp) :: zeta_fit
      real(wp) :: max_error
      character(len=:), allocatable :: message

      call build_spherical_cavities(nleb, cavities, error)
      if (allocated(error)) return

      call fit_zeta(cavities, zeta_fit, message)
      if (allocated(message)) then
         call test_failed(error, message)
         return
      end if

      if (.not. ieee_is_finite(zeta_fit) .or. zeta_fit < zeta_lower &
         .or. zeta_fit > zeta_upper) then
         call test_failed(error, "Invalid fitted zeta for nleb=" // int_string(nleb) &
            // ": " // real_string(zeta_fit))
         return
      end if

      call validate_fit(cavities, zeta_fit, max_error, message)
      if (allocated(message)) then
         call test_failed(error, message)
         return
      end if

      if (max_error > max_rel_energy_error) then
         call test_failed(error, "DROP Born fit too inaccurate for nleb=" // int_string(nleb) &
            // ", zeta=" // real_string(zeta_fit) &
            // ", max relative error=" // real_string(max_error))
         return
      end if

      if (abs(zeta_fit - cavities(1)%param%iswig_xi_born) > max_zeta_param_error) then
         call test_failed(error, "Constructor Born zeta mismatch for nleb=" // int_string(nleb) &
            // ", fitted=" // real_string(zeta_fit) &
            // ", constructor=" // real_string(cavities(1)%param%iswig_xi_born))
         return
      end if

      call plp%begin_row()
      call plp%add(nleb)
      call plp%add("fitted")
      call plp%add(zeta_fit, fmt="es14.6")
      call plp%add(max_error, fmt="es14.6")
      call plp%add("validated against Born radii set")
      call plp%end_row()
   end subroutine check_lebedev_grid

   !> Build single-atom DROP cavities for all reference radii.
   subroutine build_spherical_cavities(nleb, cavities, error)
      !> Number of Lebedev points per sphere.
      integer, intent(in) :: nleb
      !> Output DROP cavities.
      type(cavity_type_drop), intent(out) :: cavities(:)
      !> Test-drive error object.
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp) :: radii(2)
      integer :: ir

      !* Load MB16-43/H2 and reshape it into a sphere-at-origin carrier:
      !* atom 1 sits at the origin (centre of the test sphere); atom 2 is
      !* pushed to a far corner of space with a sub-millibohr radius so it
      !* contributes negligibly to either the surface area or the CPCM
      !* solve. See the module header for the tolerance budget.
      call get_structure(mol, "MB16-43", "H2")
      mol%xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      mol%xyz(:, 2) = [100.0_wp, 100.0_wp, 100.0_wp]

      do ir = 1, size(cavities)
         radii(1) = sphere_radii(ir)
         radii(2) = 1.0e-3_wp
         call new_radii_custom_atoms(radii, radius_model, cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, "new_radii_custom_atoms failed for nleb=" &
               // int_string(nleb) // ": " // trim(cavity_error%message))
            return
         end if

         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=3.0_wp, blend_1b=1.0_wp, &
               blend_2b=1.0_wp, blend_3b=1.0_wp)
            call new_cavity_drop(cavities(ir), nleb=nleb, &
               tolerance=1.0e-10_wp, proj_maxiter=150, proj_level=2, &
               radius_model=radius_model, verbose=0, debug=.false., &
               lsf_model=svdw_template, error=cavity_error)
         end block
         if (allocated(cavity_error)) then
            call test_failed(error, "new_cavity_drop failed for nleb=" &
               // int_string(nleb) // ": " // trim(cavity_error%message))
            return
         end if

         call cavities(ir)%update(mol, error=cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, "DROP cavity update failed for nleb=" &
               // int_string(nleb) // ", radius=" // real_string(radii(1)) &
               // ": " // trim(cavity_error%message))
            return
         end if

         if (allocated(radius_model)) deallocate(radius_model)
      end do
   end subroutine build_spherical_cavities

   !> Check whether DROP can use a Lebedev grid without negative weights.
   function drop_supports_lebedev_grid(nleb) result(supported)
      !> Number of Lebedev points.
      integer, intent(in) :: nleb
      !> Support flag.
      logical :: supported

      type(mctc_error), allocatable :: grid_error
      real(wp), allocatable :: xyz(:, :)
      real(wp), allocatable :: weights(:)
      integer :: order

      supported = .false.

      call lebedev_order_from_num(nleb, order, grid_error)
      if (allocated(grid_error)) return

      allocate(xyz(3, nleb), weights(nleb))
      call get_angular_grid(order, xyz, weights, grid_error)
      if (allocated(grid_error)) return

      supported = .not. any(weights < 0.0_wp)
   end function drop_supports_lebedev_grid

   !> Fit zeta by minimizing summed squared relative Born-energy errors.
   subroutine fit_zeta(cavities, zeta_fit, message)
      !> Spherical DROP cavities.
      type(cavity_type_drop), intent(inout) :: cavities(:)
      !> Fitted zeta value.
      real(wp), intent(out) :: zeta_fit
      !> Error message, allocated on failure.
      character(len=:), allocatable, intent(out) :: message

      real(wp) :: a, b, c, d, fc, fd, gr
      integer :: iter

      gr = (sqrt(5.0_wp) - 1.0_wp) / 2.0_wp

      a = zeta_lower
      b = zeta_upper
      c = b - gr * (b - a)
      d = a + gr * (b - a)

      call born_fit_objective(cavities, c, fc, message)
      if (allocated(message)) return
      call born_fit_objective(cavities, d, fd, message)
      if (allocated(message)) return

      do iter = 1, fit_maxiter
         if (abs(b - a) < fit_interval_tol) exit

         if (fc < fd) then
            b = d
            d = c
            fd = fc
            c = b - gr * (b - a)
            call born_fit_objective(cavities, c, fc, message)
            if (allocated(message)) return
         else
            a = c
            c = d
            fc = fd
            d = a + gr * (b - a)
            call born_fit_objective(cavities, d, fd, message)
            if (allocated(message)) return
         end if
      end do

      zeta_fit = 0.5_wp * (a + b)
   end subroutine fit_zeta

   !> Compute summed squared relative Born-energy errors for one zeta.
   subroutine born_fit_objective(cavities, zeta, objective, message)
      !> Spherical DROP cavities.
      type(cavity_type_drop), intent(inout) :: cavities(:)
      !> Candidate zeta value.
      real(wp), intent(in) :: zeta
      !> Summed squared relative error.
      real(wp), intent(out) :: objective
      !> Error message, allocated on failure.
      character(len=:), allocatable, intent(out) :: message

      real(wp) :: energy, reference, rel_error

      objective = 0.0_wp
      ! For a single centered charge in a one-atom spherical cavity, the CPCM
      ! energy scales exactly as 1/R. One radius is sufficient for fitting;
      ! all radii are checked after the fitted zeta is found.
      call compute_drop_born_energy(cavities(1), zeta, energy, message)
      if (allocated(message)) return

      reference = analytical_born_energy(sphere_radii(1))
      rel_error = (energy - reference) / reference
      objective = rel_error**2
   end subroutine born_fit_objective

   !> Validate the final fit and return the maximum relative error.
   subroutine validate_fit(cavities, zeta, max_error, message)
      !> Spherical DROP cavities.
      type(cavity_type_drop), intent(inout) :: cavities(:)
      !> Fitted zeta value.
      real(wp), intent(in) :: zeta
      !> Maximum relative error over radii.
      real(wp), intent(out) :: max_error
      !> Error message, allocated on failure.
      character(len=:), allocatable, intent(out) :: message

      real(wp) :: energy, reference, rel_error
      integer :: ir

      max_error = 0.0_wp
      do ir = 1, size(cavities)
         call compute_drop_born_energy(cavities(ir), zeta, energy, message)
         if (allocated(message)) return

         reference = analytical_born_energy(sphere_radii(ir))
         rel_error = abs((energy - reference) / reference)
         max_error = max(max_error, rel_error)
      end do
   end subroutine validate_fit

   !> Compute DROP CPCM solvation energy for a central point charge.
   subroutine compute_drop_born_energy(cavity, zeta, energy, message)
      !> Spherical DROP cavity.
      type(cavity_type_drop), intent(inout) :: cavity
      !> Candidate Gaussian width scale.
      real(wp), intent(in) :: zeta
      !> CPCM solvation energy.
      real(wp), intent(out) :: energy
      !> Error message, allocated on failure.
      character(len=:), allocatable, intent(out) :: message

      type(mctc_error), allocatable :: cavity_error
      real(wp), allocatable :: amat(:, :)
      real(wp), allocatable :: phi(:)
      real(wp), allocatable :: rhs(:)
      real(wp), allocatable :: sigma(:)
      integer :: igrid

      cavity%iswig%swx = zeta
      call cavity%compute_gaussians(cavity_error)
      if (allocated(cavity_error)) then
         message = "compute_gaussians failed: " // trim(cavity_error%message)
         return
      end if

      call cavity%Amat012_rA(amat, error=cavity_error)
      if (allocated(cavity_error)) then
         message = "DROP A-matrix assembly failed: " // trim(cavity_error%message)
         return
      end if

      allocate(phi(cavity%ngrid), rhs(cavity%ngrid), sigma(cavity%ngrid))
      do igrid = 1, cavity%ngrid
         phi(igrid) = source_charge / norm2(cavity%xyz(:, igrid))
      end do

      rhs = -feps * phi
      call solve_pcm_cholesky(amat, rhs, sigma, error=cavity_error)
      if (allocated(cavity_error)) then
         message = "CPCM Cholesky solve failed: " // trim(cavity_error%message)
         return
      end if

      energy = 0.5_wp * dot_product(sigma, phi)
   end subroutine compute_drop_born_energy

   !> Analytical Born solvation energy for a central unit charge in a sphere.
   pure function analytical_born_energy(radius) result(energy)
      !> Sphere radius in bohr.
      real(wp), intent(in) :: radius
      !> Analytical Born energy in hartree.
      real(wp) :: energy

      energy = -0.5_wp * feps * source_charge**2 / radius
   end function analytical_born_energy

   !> Convert integer to allocatable string.
   function int_string(value) result(string)
      !> Integer value.
      integer, intent(in) :: value
      !> Formatted string.
      character(len=:), allocatable :: string
      character(len=32) :: buffer

      write(buffer, '(i0)') value
      string = trim(buffer)
   end function int_string

   !> Convert real to allocatable string.
   function real_string(value) result(string)
      !> Real value.
      real(wp), intent(in) :: value
      !> Formatted string.
      character(len=:), allocatable :: string
      character(len=48) :: buffer

      write(buffer, '(es22.14)') value
      string = trim(adjustl(buffer))
   end function real_string

end module test_cavity_drop_born_fit
