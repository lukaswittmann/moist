!> DROP optional surface properties against independent geometric references
!>
!> - grid point density of isolated spheres against `num_leb/(4 pi R^2)`
!> - `rho_scal0` of isolated and overlapping atoms
!> - curvatures, normals, point-to-center distances, area and volume of isolated spheres
!> - `analyze_cavity` report against the cavity's own areas and curvatures
module test_cavity_drop_properties
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed, to_string
   use moist_cavity_drop, only: cavity_type_drop
   use test_helpers, only: build_test_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: new_radii_custom_atoms, default_cpcm_radii, radius_type
   use moist_context, only: moist_context_type, new_context
   implicit none(type, external)
   private

   public :: collect_cavity_drop_properties

   !> Projection tolerance
   real(wp), parameter :: PROJ_TOL = 1.0e-12_wp
   !> Projection iteration cap
   integer, parameter :: PROJ_MAXITER = 1000
   !> Projection level, 1 = SLSQP only
   integer, parameter :: PROJ_LEVEL = 2
   !> Sphere radii of the isolated-atom fixtures, bohr
   real(wp), parameter :: SPHERE_RADII(3) = [2.0_wp, 3.5_wp, 5.0_wp]
   !> Density kernel lengths paired with `SPHERE_RADII`, bohr
   real(wp), parameter :: SPHERE_KERNEL_H(3) = [0.5_wp, 0.8_wp, 1.0_wp]
   !> Lebedev points of the density fixtures
   integer, parameter :: DENSITY_NLEB = 974
   !> Relative per-point tolerance of the isolated-sphere density
   !>
   !> Observed worst case 3.1e-3 (R = 5, h = 1, fewest points under the kernel)
   real(wp), parameter :: DENSITY_REL_THR = 5.0e-3_wp
   !> Absolute tolerance of exact sphere identities
   real(wp), parameter :: SPHERE_ABS_THR = 1.0e-9_wp
   !> Absolute tolerance of the isolated-sphere curvature identities
   !>
   !> Observed worst case 6e-16
   real(wp), parameter :: CURVATURE_ABS_THR = 1.0e-12_wp
   !> Tolerance of the point-wise `rho_scal0 = 1` and `cpjac_scal0 = 1` identities for isolated spheres
   !>
   !> Observed worst case 1.5e-12 for two atoms 60 bohr apart, where the level-set tails remain
   real(wp), parameter :: SCAL_ISOLATED_THR = 1.0e-9_wp
   !> Bound of the mean `|rho_scal0 - 1|` of a soft-blend molecule, observed 0.22
   real(wp), parameter :: SCAL_MEAN_DEV_THR = 0.3_wp
   !> Bound of the area-weighted mean `rho_scal0 - 1` of a soft-blend molecule, observed 0.047
   real(wp), parameter :: SCAL_AREA_MEAN_THR = 0.1_wp
   !> Tolerance of numbers read back from the six-decimal report
   real(wp), parameter :: REPORT_ABS_THR = 2.0e-6_wp
   !> Lebedev sizes of the integral-identity fixtures
   integer, parameter :: IDENT_NLEB(2) = [194, 590]
   !> Soft-blend parameters of the integral-identity fixtures
   real(wp), parameter :: IDENT_BLEND_K(2) = [1.0_wp, 2.0_wp]
   !> Number of integral-identity fixtures
   integer, parameter :: IDENT_NFIX = 4
   !> Relative tolerance of the closure and the divergence theorem about another point
   !>
   !> Observed worst case 3.1e-3 (closure, MB16-43 `16`, blend_k 1, 194 points), no systematic gain with
   !> the Lebedev size because the switched and projected points limit the quadrature
   real(wp), parameter :: IDENT_VECTOR_THR = 5.0e-3_wp
   !> Relative tolerance of the Gauss-Bonnet identity per `IDENT_NLEB`
   !>
   !> Observed worst case 6.8e-3 with 194 points and 1.5e-3 with 590 points (MB16-43 `16`, blend_k 2)
   real(wp), parameter :: IDENT_GAUSS_THR(2) = [1.0e-2_wp, 3.0e-3_wp]
   !> Relative tolerance of the Minkowski identity per `IDENT_NLEB`
   !>
   !> Observed worst case 3.3e-3 with 194 points (MB16-43 `16`, blend_k 2) and 9.0e-4 with 590 points
   !> (MB16-43 `16`, blend_k 2); a uniform relative error in H shows up at its own size
   real(wp), parameter :: IDENT_MINK_THR(2) = [5.0e-3_wp, 2.0e-3_wp]
   !> Tolerance of all identities on the separated spheres, which have no switched points
   !>
   !> Observed worst case 2.5e-14 (Minkowski)
   real(wp), parameter :: IDENT_EXACT_THR = 1.0e-12_wp
   !> Tolerance of the divergence theorem about the volume origin, which is the definition of the volume
   !>
   !> Observed worst case 4.2e-15
   real(wp), parameter :: IDENT_VOLUME_THR = 1.0e-10_wp

contains

   !> Collect all exported unit tests
   !>
   !> @param[out] testsuite  Collection of unit tests
   subroutine collect_cavity_drop_properties(testsuite)
      !> Collection of unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("density_single_sphere", test_density_single_sphere), &
                  new_unittest("rho_scal0_molecule", test_rho_scal0_molecule), &
                  new_unittest("curvature_single_sphere", test_curvature_single_sphere), &
                  new_unittest("sphere_geometry", test_sphere_geometry), &
                  new_unittest("analyze_cavity_output", test_analyze_cavity_output), &
                  new_unittest("integral_identities", test_integral_identities) &
                  ]
   end subroutine collect_cavity_drop_properties

   !> SvdW level set of the fixtures
   subroutine new_svdw_lsf(lsf, blend_k)
      type(moist_cavity_drop_lsf_svdw_type), intent(out) :: lsf
      real(wp), intent(in) :: blend_k

      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=blend_k, blend_2b=1.0_wp, blend_3b=0.0_wp))
   end subroutine new_svdw_lsf

   !> DROP settings of the fixtures
   function fixture_param(num_leb, h, do_density, do_curv) result(param)
      integer, intent(in) :: num_leb
      real(wp), intent(in) :: h
      logical, intent(in) :: do_density, do_curv
      type(moist_cavity_drop_parameters_type) :: param

      param = moist_cavity_drop_parameters_type(num_leb=num_leb, tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
         proj_level=PROJ_LEVEL, rho_grid_h=h, do_grid_density=do_density, do_curvature=do_curv)
   end function fixture_param

   !> Turn a cavity build error into a test failure
   subroutine fail_on_build(error, build_error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable, intent(in) :: build_error

      if (allocated(build_error)) call test_failed(error, "DROP cavity build failed: "//trim(build_error%message))
   end subroutine fail_on_build

   !> Isolated-sphere DROP cavity with a custom radius
   !>
   !> @param[out]   error   Test failure state
   !> @param[out]   cavity  Updated cavity
   !> @param[inout] ctx     Context borrowed by the cavity
   !> @param[in]    radius  Sphere radius in bohr
   !> @param[in]    center  Atom position
   !> @param[in]    num_leb Lebedev points
   !> @param[in]    h       Density kernel length
   !> @param[in]    do_density Switch on the grid point density
   !> @param[in]    do_curv    Switch on the curvatures
   subroutine build_sphere(error, cavity, ctx, radius, center, num_leb, h, do_density, do_curv)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      type(moist_context_type), target, intent(inout) :: ctx
      real(wp), intent(in) :: radius
      real(wp), intent(in) :: center(3)
      integer, intent(in) :: num_leb
      real(wp), intent(in) :: h
      logical, intent(in) :: do_density, do_curv

      type(structure_type) :: mol
      class(radius_type), allocatable :: radius_model
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(mctc_error), allocatable :: radius_error, cavity_error

      call new(mol, [18], reshape(center, [3, 1]))
      call new_radii_custom_atoms([radius], radius_model, radius_error)
      if (allocated(radius_error)) then
         call test_failed(error, "new_radii_custom_atoms failed: "//trim(radius_error%message))
         return
      end if
      call new_svdw_lsf(lsf, 1.0_wp)
      allocate (cavity)
      call build_test_cavity_drop(mol, num_leb, ctx, radius_model, lsf, cavity, cavity_error, &
         param=fixture_param(num_leb, h, do_density, do_curv))
      call fail_on_build(error, cavity_error)
   end subroutine build_sphere

   !> Per-point density of an isolated sphere against the uniform-grid value
   !>
   !> - the density is a 2D kernel sum over surface weights, so an isolated sphere
   !>   of radius R carries `num_leb/(4 pi R^2)` points per area everywhere
   !> - the kernel acts on chord distances, for which the sphere integral of a radial
   !>   kernel equals the planar one, so the continuum value is exact for `2h <= 2R`
   !>   and only the Lebedev quadrature of the kernel remains
   !> - an isolated sphere has no projection and no switching, so the soft and
   !>   hard-sphere densities coincide and `rho_scal0 = 1`
   subroutine test_density_single_sphere(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      real(wp) :: expect
      integer :: ir

      call new_context(ctx, nthreads=0, verbosity=0)
      do ir = 1, size(SPHERE_RADII)
         call build_sphere(error, cavity, ctx, SPHERE_RADII(ir), [0.0_wp, 0.0_wp, 0.0_wp], DENSITY_NLEB, &
            SPHERE_KERNEL_H(ir), .true., .false.)
         if (allocated(error)) return
         call check(error, cavity%ngrid, DENSITY_NLEB, message="Isolated sphere lost grid points for R = "// &
                    to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return

         expect = real(DENSITY_NLEB, wp)/(4.0_wp*pi*SPHERE_RADII(ir)**2)
         call check(error, maxval(abs(cavity%rho_grid/expect - 1.0_wp)) < DENSITY_REL_THR, &
            "Isolated sphere density deviates from num_leb/(4 pi R^2) for R = "//to_string(SPHERE_RADII(ir))// &
            ": expected "//to_string(expect)//", range "//to_string(minval(cavity%rho_grid))//" to "// &
            to_string(maxval(cavity%rho_grid)))
         if (allocated(error)) return

         call check(error, maxval(abs(cavity%rho_scal0 - 1.0_wp)) < SCAL_ISOLATED_THR, &
            "rho_scal0 of an isolated sphere is not 1 for R = "//to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return
         deallocate (cavity)
      end do
   end subroutine test_density_single_sphere

   !> `rho_scal0` of non-overlapping and strongly overlapping atoms
   !>
   !> - `rho_scal0` is the hard-sphere over the soft density; the soft density
   !>   weights its neighbours with the Jacobian-scaled weights, so it measures
   !>   surface area and not point count, and `rho_scal0` stays near 1 instead of
   !>   tracking `1/cpjac_scal0`
   !> - two atoms 60 bohr apart keep every point on its sphere, so `rho_scal0 = 1`
   !>   point-wise
   !> - MB16-43 `16` at `blend_k = 1` has `cpjac_scal0` from 0.009 to 5.6; the bounds
   !>   on the mean and the area-weighted mean come from the observed spread
   subroutine test_rho_scal0_molecule(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      class(radius_type), allocatable :: radius_model
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(mctc_error), allocatable :: cavity_error
      real(wp) :: xyz(3, 2), mean_dev, area_mean

      radius_model = default_cpcm_radii()
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_svdw_lsf(lsf, 1.0_wp)

      xyz = 0.0_wp
      xyz(3, 2) = 60.0_wp
      call new(mol, [8, 8], xyz)
      allocate (cavity)
      call build_test_cavity_drop(mol, 194, ctx, radius_model, lsf, cavity, cavity_error, &
         param=fixture_param(194, 1.0_wp, .true., .false.))
      call fail_on_build(error, cavity_error)
      if (allocated(error)) return
      call check(error, maxval(abs(cavity%rho_scal0 - 1.0_wp)) < SCAL_ISOLATED_THR, &
         "rho_scal0 of two distant atoms is not 1")
      if (allocated(error)) return
      call check(error, maxval(abs(cavity%cpjac_scal0 - 1.0_wp)) < SCAL_ISOLATED_THR, &
         "cpjac_scal0 of two distant atoms is not 1, the fixture is not isolated")
      if (allocated(error)) return
      deallocate (cavity)

      call get_structure(mol, "MB16-43", "16")
      allocate (cavity)
      call build_test_cavity_drop(mol, 194, ctx, radius_model, lsf, cavity, cavity_error, &
         param=fixture_param(194, 1.0_wp, .true., .false.))
      call fail_on_build(error, cavity_error)
      if (allocated(error)) return
      call check(error, all(cavity%rho_scal0 > 0.0_wp), "rho_scal0 is not positive everywhere")
      if (allocated(error)) return
      call check(error, maxval(cavity%cpjac_scal0) > 2.0_wp .and. minval(cavity%cpjac_scal0) < 0.5_wp, &
         "MB16-43 16 at blend_k = 1 no longer spans a wide Jacobian range")
      if (allocated(error)) return

      mean_dev = sum(abs(cavity%rho_scal0 - 1.0_wp))/real(cavity%ngrid, wp)
      call check(error, mean_dev < SCAL_MEAN_DEV_THR, "Mean |rho_scal0 - 1| is "//to_string(mean_dev))
      if (allocated(error)) return
      area_mean = sum(cavity%a*cavity%rho_scal0)/sum(cavity%a)
      call check(error, abs(area_mean - 1.0_wp) < SCAL_AREA_MEAN_THR, &
         "Area-weighted mean rho_scal0 is "//to_string(area_mean))
      if (allocated(error)) return
   end subroutine test_rho_scal0_molecule

   !> Principal, mean and Gaussian curvature of isolated spheres
   !>
   !> - a sphere of radius R has `k1 = k2 = 1/R`, `K_M = 1/R`, `K_G = 1/R^2`
   !> - the convention is positive curvature for a convex surface, `k1 >= k2`
   !> - the atom sits off the origin
   subroutine test_curvature_single_sphere(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      real(wp) :: kappa
      integer :: ir

      call new_context(ctx, nthreads=0, verbosity=0)
      do ir = 1, size(SPHERE_RADII)
         call build_sphere(error, cavity, ctx, SPHERE_RADII(ir), [0.3_wp, -0.7_wp, 1.1_wp], &
          194, 1.0_wp, .false., .true.)
         if (allocated(error)) return
         kappa = 1.0_wp/SPHERE_RADII(ir)
         call check(error, maxval(abs(cavity%k1 - kappa)) < CURVATURE_ABS_THR, &
            "k1 differs from 1/R for R = "//to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return
         call check(error, maxval(abs(cavity%k2 - kappa)) < CURVATURE_ABS_THR, &
            "k2 differs from 1/R for R = "//to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return
         call check(error, maxval(abs(cavity%KM - kappa)) < CURVATURE_ABS_THR, &
            "Mean curvature differs from 1/R for R = "//to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return
         call check(error, maxval(abs(cavity%KG - kappa**2)) < CURVATURE_ABS_THR, &
            "Gaussian curvature differs from 1/R^2 for R = "//to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return
         call check(error, all(cavity%k1 >= cavity%k2), "k1 < k2 at some point")
         if (allocated(error)) return
         deallocate (cavity)
      end do
   end subroutine test_curvature_single_sphere

   !> Normals, center distances, displacements, area and volume of isolated spheres
   !>
   !> - the outward normal is the unit radial vector, `r_iI0 = R`, and the
   !>   anchor-to-surface displacement `rho` vanishes
   !> - area `4 pi R^2` and volume `4 pi R^3/3` follow from the weights alone, and
   !>   the volume does not depend on the atom position
   subroutine test_sphere_geometry(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      real(wp), parameter :: CENTER(3) = [0.3_wp, -0.7_wp, 1.1_wp]
      real(wp) :: radial(3), dev
      integer :: ir, i

      call new_context(ctx, nthreads=0, verbosity=0)
      do ir = 1, size(SPHERE_RADII)
         call build_sphere(error, cavity, ctx, SPHERE_RADII(ir), CENTER, 194, 1.0_wp, .false., .false.)
         if (allocated(error)) return

         dev = 0.0_wp
         do i = 1, cavity%ngrid
            radial = (cavity%xyz(:, i) - CENTER)/SPHERE_RADII(ir)
            dev = max(dev, norm2(cavity%normal0(:, i) - radial))
         end do
         call check(error, dev < SPHERE_ABS_THR, "Normal is not the outward radial unit vector for R = "// &
                    to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return
         call check(error, maxval(abs(cavity%r_iI0 - SPHERE_RADII(ir))) < SPHERE_ABS_THR, &
            "r_iI0 differs from R for R = "//to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return
         call check(error, maxval(abs(cavity%rho)) < SPHERE_ABS_THR, &
            "Anchor displacement rho of an isolated sphere is not zero for R = "//to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return
         call check(error, cavity%total_area, 4.0_wp*pi*SPHERE_RADII(ir)**2, thr=SPHERE_ABS_THR, &
                    message="Sphere area differs from 4 pi R^2 for R = "//to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return
         call check(error, cavity%total_volume, 4.0_wp*pi*SPHERE_RADII(ir)**3/3.0_wp, thr=SPHERE_ABS_THR, &
                    message="Sphere volume differs from 4 pi R^3/3 for R = "//to_string(SPHERE_RADII(ir)))
         if (allocated(error)) return
         deallocate (cavity)
      end do
   end subroutine test_sphere_geometry

   !> Read the numbers behind the `...` filler of a report line
   !>
   !> @param[in]  line  Report line
   !> @param[out] vals  Leading numeric tokens of the value field, units skipped
   !> @param[out] found Line carries the filler and `n` numbers
   subroutine report_numbers(line, n, vals, found)
      character(len=*), intent(in) :: line
      integer, intent(in) :: n
      real(wp), intent(out) :: vals(n)
      logical, intent(out) :: found

      integer :: ipos, ival, istart, iend, stat
      character(len=len(line)) :: rest

      vals = 0.0_wp
      found = .false.
      ipos = index(line, "...")
      if (ipos == 0) return
      rest = adjustl(line(ipos + 3:))
      ival = 0
      istart = 1
      do while (ival < n .and. istart <= len_trim(rest))
         iend = index(rest(istart:), " ")
         if (iend == 0) iend = len_trim(rest) - istart + 2
         iend = istart + iend - 2
         if (iend >= istart) then
            read (rest(istart:iend), *, iostat=stat) vals(ival + 1)
            if (stat == 0) ival = ival + 1
         end if
         istart = iend + 2
         do while (istart <= len_trim(rest))
            if (rest(istart:istart) /= " ") exit
            istart = istart + 1
         end do
      end do
      found = ival == n
   end subroutine report_numbers

   !> Report of `analyze_cavity`, run through verbosity 3 of the update
   !>
   !> - convex, saddle and pit areas partition the cavity, so they sum to
   !>   `total_area`, and their percentages to 100
   !> - the printed k1 and k2 ranges equal those of the cavity arrays
   !> - an isolated sphere is fully convex with area `4 pi R^2` in the convex class
   !>   and `k1 = k2 = 1/R`
   !> - the watershed patches cover part of the surface
   subroutine test_analyze_cavity_output(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      class(radius_type), allocatable :: radius_model
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(mctc_error), allocatable :: cavity_error
      real(wp) :: v(3), area_sum, pct_sum, patch_area, patch_pct
      logical :: found, got_convex, got_saddle, got_pit, got_k1, got_k2, got_patch, got_status, got_probe
      integer :: iu, stat
      character(len=256) :: line

      radius_model = default_cpcm_radii()
      call get_structure(mol, "MB16-43", "16")
      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat, 0, more="could not open the context scratch unit")
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, verbosity=3, unit=iu)
      call new_svdw_lsf(lsf, 2.0_wp)
      allocate (cavity)
      call build_test_cavity_drop(mol, 194, ctx, radius_model, lsf, cavity, cavity_error, &
         param=fixture_param(194, 1.0_wp, .false., .true.))
      call fail_on_build(error, cavity_error)
      if (allocated(error)) then
         close (iu)
         return
      end if

      area_sum = 0.0_wp
      pct_sum = 0.0_wp
      patch_area = -1.0_wp
      got_convex = .false.; got_saddle = .false.; got_pit = .false.
      got_k1 = .false.; got_k2 = .false.; got_patch = .false.
      rewind (iu)
      do
         read (iu, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (index(line, "Convex  (k1>0, k2>0)") > 0) then
            call report_numbers(line, 2, v, found)
            got_convex = found
            area_sum = area_sum + v(1)
            pct_sum = pct_sum + v(2)
         else if (index(line, "Saddle  (k1>0, k2<0)") > 0) then
            call report_numbers(line, 2, v, found)
            got_saddle = found
            area_sum = area_sum + v(1)
            pct_sum = pct_sum + v(2)
         else if (index(line, "Pit     (k1<0, k2<0)") > 0) then
            call report_numbers(line, 2, v, found)
            got_pit = found
            area_sum = area_sum + v(1)
            pct_sum = pct_sum + v(2)
         else if (index(line, "k1 [min/mean/max]") > 0) then
            call report_numbers(line, 3, v, found)
            got_k1 = found
            call check(error, abs(v(1) - minval(cavity%k1)) < REPORT_ABS_THR .and. &
               abs(v(3) - maxval(cavity%k1)) < REPORT_ABS_THR .and. &
               abs(v(2) - sum(cavity%k1)/real(cavity%ngrid, wp)) < REPORT_ABS_THR, &
               "Reported k1 min/mean/max differ from the cavity: "//trim(line))
         else if (index(line, "k2 [min/mean/max]") > 0) then
            call report_numbers(line, 3, v, found)
            got_k2 = found
            call check(error, abs(v(1) - minval(cavity%k2)) < REPORT_ABS_THR .and. &
               abs(v(3) - maxval(cavity%k2)) < REPORT_ABS_THR .and. &
               abs(v(2) - sum(cavity%k2)/real(cavity%ngrid, wp)) < REPORT_ABS_THR, &
               "Reported k2 min/mean/max differ from the cavity: "//trim(line))
         else if (index(line, "Patch area total") > 0) then
            call report_numbers(line, 2, v, found)
            got_patch = found
            patch_area = v(1)
            patch_pct = v(2)
         end if
         if (allocated(error)) exit
      end do
      close (iu)
      if (allocated(error)) return

      call check(error, got_convex .and. got_saddle .and. got_pit .and. got_k1 .and. got_k2 .and. got_patch, &
         "analyze_cavity report lacks a regime, curvature or patch line")
      if (allocated(error)) return
      call check(error, abs(area_sum - cavity%total_area) < 3.0_wp*REPORT_ABS_THR, &
         "Curvature regime areas sum to "//to_string(area_sum)//" instead of "//to_string(cavity%total_area))
      if (allocated(error)) return
      call check(error, abs(pct_sum - 100.0_wp) < 3.0_wp*REPORT_ABS_THR, &
         "Curvature regime percentages sum to "//to_string(pct_sum))
      if (allocated(error)) return
      call check(error, patch_area > 0.0_wp .and. patch_area < cavity%total_area, &
         "Watershed patch area "//to_string(patch_area)//" is not within the cavity area")
      if (allocated(error)) return
      call check(error, abs(patch_pct - 100.0_wp*patch_area/cavity%total_area) < 1.0e-4_wp, &
         "Watershed patch percentage is not its area over the total")
      if (allocated(error)) return
      deallocate (cavity)

      ! Isolated sphere
      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat, 0, more="could not open the context scratch unit")
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, verbosity=3, unit=iu)
      call build_sphere(error, cavity, ctx, SPHERE_RADII(2), [0.0_wp, 0.0_wp, 0.0_wp], 194, 1.0_wp, .false., .true.)
      if (allocated(error)) then
         close (iu)
         return
      end if
      got_convex = .false.; got_status = .false.; got_probe = .false.
      v = 0.0_wp
      rewind (iu)
      do
         read (iu, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (index(line, "Convex  (k1>0, k2>0)") > 0) then
            call report_numbers(line, 2, v(1:2), found)
            got_convex = found
         end if
         if (index(line, "fully convex") > 0) got_status = .true.
         if (index(line, "Probe-region   ") > 0) got_probe = .true.
      end do
      close (iu)
      call check(error, got_convex .and. got_status, "Sphere report lacks the convex line or the fully convex status")
      if (allocated(error)) return
      call check(error, .not. got_probe, "Sphere report claims a probe region")
      if (allocated(error)) return
      call check(error, abs(v(1) - 4.0_wp*pi*SPHERE_RADII(2)**2) < REPORT_ABS_THR &
         .and. abs(v(2) - 100.0_wp) < REPORT_ABS_THR, &
         "Sphere convex area is "//to_string(v(1))//" ("//to_string(v(2))//" %)")
      if (allocated(error)) return
   end subroutine test_analyze_cavity_output

   !> Fixture of the integral identities
   !>
   !> 1 two overlapping spheres, 2 three overlapping spheres, 3 MB16-43 `16`, 4 two separated spheres
   subroutine build_identity_fixture(error, cavity, ctx, ifix, blend_k, num_leb)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      type(moist_context_type), target, intent(inout) :: ctx
      integer, intent(in) :: ifix
      real(wp), intent(in) :: blend_k
      integer, intent(in) :: num_leb

      type(structure_type) :: mol
      class(radius_type), allocatable :: radius_model
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(mctc_error), allocatable :: radius_error, cavity_error

      select case (ifix)
      case (1)
         call new(mol, [18, 18], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp, 0.3_wp, 4.0_wp], [3, 2]))
         call new_radii_custom_atoms([3.5_wp, 3.0_wp], radius_model, radius_error)
      case (2)
         call new(mol, [18, 18, 18], &
            reshape([0.0_wp, 0.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, 0.5_wp, 2.0_wp, 3.5_wp, -0.4_wp], [3, 3]))
         call new_radii_custom_atoms([3.0_wp, 3.2_wp, 2.8_wp], radius_model, radius_error)
      case (3)
         call get_structure(mol, "MB16-43", "16")
         radius_model = default_cpcm_radii()
      case default
         call new(mol, [18, 18], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 11.0_wp], [3, 2]))
         call new_radii_custom_atoms([2.0_wp, 2.5_wp], radius_model, radius_error)
      end select
      if (allocated(radius_error)) then
         call test_failed(error, "new_radii_custom_atoms failed: "//trim(radius_error%message))
         return
      end if
      call new_svdw_lsf(lsf, blend_k)
      allocate (cavity)
      call build_test_cavity_drop(mol, num_leb, ctx, radius_model, lsf, cavity, cavity_error, &
         param=fixture_param(num_leb, 1.0_wp, .false., .true.))
      call fail_on_build(error, cavity_error)
   end subroutine build_identity_fixture

   !> Surface integrals of the closed-surface identities
   !>
   !> @param[in]  cavity   Updated cavity with curvatures
   !> @param[in]  c        Reference point of the position vector
   !> @param[out] closure  `sum w n`
   !> @param[out] div      `sum w (r - c).n`
   !> @param[out] gauss    `sum w k1 k2`
   !> @param[out] mink     `sum w (1 - H (r - c).n)`
   !> @param[out] imc      `sum w H`
   subroutine surface_integrals(cavity, c, closure, div, gauss, mink, imc)
      type(cavity_type_drop), intent(in) :: cavity
      real(wp), intent(in) :: c(3)
      real(wp), intent(out) :: closure(3), div, gauss, mink, imc

      real(wp) :: rn
      integer :: i

      closure = 0.0_wp
      div = 0.0_wp
      gauss = 0.0_wp
      mink = 0.0_wp
      imc = 0.0_wp
      do i = 1, cavity%ngrid
         rn = dot_product(cavity%xyz(:, i) - c, cavity%normal0(:, i))
         closure = closure + cavity%a(i)*cavity%normal0(:, i)
         div = div + cavity%a(i)*rn
         gauss = gauss + cavity%a(i)*cavity%k1(i)*cavity%k2(i)
         mink = mink + cavity%a(i)*(1.0_wp - cavity%KM(i)*rn)
         imc = imc + cavity%a(i)*cavity%KM(i)
      end do
   end subroutine surface_integrals

   !> Integral identities of smooth closed surfaces on overlapping-sphere and molecular cavities
   !>
   !> - closure `sum w n = 0`, divergence theorem `sum w (r - c).n = 3 V` for the volume origin and for
   !>   another point, Gauss-Bonnet `sum w k1 k2 = 2 pi chi` and Minkowski `sum w (1 - H (r - c).n) = 0`
   !> - with outward normals and positive curvature of a convex surface, a sphere has `H = 1/R` and
   !>   `(r - c).n = R`, so the Minkowski integrand vanishes there point by point
   !> - `chi = 2` for the connected fixtures and `chi = 4` for two separated spheres
   !> - the Gauss-Bonnet and Minkowski errors fall with the Lebedev degree, so the larger
   !>   size is held to tighter bounds; the closure is limited by the switched points instead
   !> - the separated spheres have no switched points and hold every identity to round-off
   subroutine test_integral_identities(error)
      !> Test failure state
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      real(wp), parameter :: OTHER(3) = [2.0_wp, -1.5_wp, 3.0_wp]
      real(wp) :: closure(3), div, gauss, mink, imc, chi, v3, area, err(5)
      real(wp) :: thr_vector, thr_gauss, thr_mink
      integer :: ifix, ik, il

      call new_context(ctx, nthreads=0, verbosity=0)
      do ifix = 1, IDENT_NFIX
         chi = merge(4.0_wp, 2.0_wp, ifix == 4)
         do ik = 1, size(IDENT_BLEND_K)
            do il = 1, size(IDENT_NLEB)
               if (ifix == 4) then
                  thr_vector = IDENT_EXACT_THR
                  thr_gauss = IDENT_EXACT_THR
                  thr_mink = IDENT_EXACT_THR
               else
                  thr_vector = IDENT_VECTOR_THR
                  thr_gauss = IDENT_GAUSS_THR(il)
                  thr_mink = IDENT_MINK_THR(il)
               end if
               call build_identity_fixture(error, cavity, ctx, ifix, IDENT_BLEND_K(ik), IDENT_NLEB(il))
               if (allocated(error)) return
               area = sum(cavity%a)
               v3 = 3.0_wp*cavity%total_volume
               call surface_integrals(cavity, cavity%get_volume_origin(), closure, div, gauss, mink, imc)
               err(1) = norm2(closure)/area
               err(2) = abs(div - v3)/v3
               err(3) = abs(gauss - 2.0_wp*pi*chi)/(2.0_wp*pi*chi)
               err(4) = abs(mink)/area
               call surface_integrals(cavity, OTHER, closure, div, gauss, mink, imc)
               err(5) = abs(div - v3)/v3
               call check(error, err(1) < thr_vector, "Closure sum w n is not zero for fixture " &
                  //to_string(ifix)// &
                  ", blend_k "//to_string(IDENT_BLEND_K(ik))//", num_leb "//to_string(IDENT_NLEB(il))//": " &
                  //to_string(err(1)))
               if (allocated(error)) return
               call check(error, err(2) < IDENT_VOLUME_THR, "Divergence theorem fails at the volume origin for fixture "// &
                  to_string(ifix)//", blend_k "//to_string(IDENT_BLEND_K(ik))//", num_leb " &
                  //to_string(IDENT_NLEB(il))//": "// &
                  to_string(err(2)))
               if (allocated(error)) return
               call check(error, err(5) < thr_vector, "Divergence theorem fails at another point for fixture "// &
                  to_string(ifix)//", blend_k "//to_string(IDENT_BLEND_K(ik))//", num_leb " &
                  //to_string(IDENT_NLEB(il))//": "// &
                  to_string(err(5)))
               if (allocated(error)) return
               call check(error, err(3) < thr_gauss, "Gauss-Bonnet fails for fixture "//to_string(ifix) &
                  //", blend_k "// &
                  to_string(IDENT_BLEND_K(ik))//", num_leb "//to_string(IDENT_NLEB(il))//": integral " &
                  //to_string(gauss)// &
                  " against "//to_string(2.0_wp*pi*chi))
               if (allocated(error)) return
               call check(error, err(4) < thr_mink, "Minkowski formula fails for fixture "//to_string(ifix)// &
                  ", blend_k "//to_string(IDENT_BLEND_K(ik))//", num_leb "//to_string(IDENT_NLEB(il)) &
                  //": "//to_string(err(4)))
               if (allocated(error)) return
               deallocate (cavity)
            end do
         end do
      end do
   end subroutine test_integral_identities

end module test_cavity_drop_properties
