!> Tests for the marching-cubes cavity
module test_cavity_marchingcubes
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use moist_cavity, only: cavity_type_marchingcubes, new_cavity_marchingcubes
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: new_radii_custom_atoms, radius_type, default_cpcm_radii
   use moist_context, only: moist_context_type, new_context

   implicit none
   private

   public :: collect_cavity_marchingcubes

   !> Finest marching-cubes spacing used throughout; the adaptive refinement
   !> reaches roughly 0.2 % on area and 0.4 % on volume for a sphere at this
   !> setting, so the checks below allow 1 %.
   real(wp), parameter :: MC_SPACING = 0.2_wp
   !> Relative tolerance on the integrated area and volume
   real(wp), parameter :: REL_THR = 1.0E-2_wp
   !> Sphere radius of the single-atom fixture, bohr
   real(wp), parameter :: SPHERE_RADIUS = 3.0_wp

   !> Radii of the analytic multi-sphere fixtures, bohr
   real(wp), parameter :: MULTI_RADII(3) = [2.0_wp, 2.5_wp, 3.0_wp]
   !> Radius of the enclosing sphere of the nested fixture, bohr
   real(wp), parameter :: ENCLOSING_RADIUS = 6.0_wp

   !> One molecular reference case, integrated at the reference spacing of
   !> 0.01 bohr by test/dev/test_cavity_drop_integration.f90
   type :: integration_case_type
      !> Blending sharpness k
      real(wp) :: blend_k
      !> Two-body blending weight
      real(wp) :: blend_2b
      !> Three-body blending weight
      real(wp) :: blend_3b
      !> Dataset name used with mstore
      character(len=12) :: dataset
      !> Structure identifier inside the dataset
      character(len=7) :: structure
      !> Reference marching-cubes area, bohr^2
      real(wp) :: mc_area
      !> Reference marching-cubes volume, bohr^3
      real(wp) :: mc_volume
   end type integration_case_type

   !> Relative tolerance of the molecular cases
   real(wp), parameter :: MOL_REL_THR = 2.0E-3_wp

   !> References taken from test/unit/test_cavity/drop/integration.f90, where they
   !> serve as the target the DROP cavity is validated against. Here they check the
   !> integrator that produced them, so a coarser grid must still reproduce them.
   type(integration_case_type), parameter :: cases(*) = [ &
      integration_case_type(1.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "16     ", 1033.798426_wp, 3088.935546_wp), &
      integration_case_type(2.0_wp, 1.0_wp, 0.0_wp, "MB16-43     ", "16     ", 739.232911_wp, 1734.806473_wp), &
      integration_case_type(2.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "16     ", 775.484041_wp, 1921.288922_wp), &
      integration_case_type(3.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "16     ", 703.771421_wp, 1567.182120_wp), &
      integration_case_type(5.0_wp, 1.0_wp, 1.0_wp, "MB16-43     ", "16     ", 690.487755_wp, 1425.662212_wp), &
      integration_case_type(2.0_wp, 1.0_wp, 0.0_wp, "But14diol   ", "30     ", 577.749734_wp, 1207.837574_wp), &
      integration_case_type(2.0_wp, 1.0_wp, 1.0_wp, "But14diol   ", "30     ", 634.177253_wp, 1431.003011_wp), &
      integration_case_type(3.0_wp, 1.0_wp, 0.0_wp, "But14diol   ", "30     ", 519.505143_wp, 980.963503_wp), &
      integration_case_type(3.0_wp, 1.0_wp, 1.0_wp, "But14diol   ", "30     ", 533.169514_wp, 1051.688680_wp), &
      integration_case_type(5.0_wp, 1.0_wp, 1.0_wp, "But14diol   ", "30     ", 496.399532_wp, 882.073907_wp), &
      integration_case_type(10.0_wp, 1.0_wp, 0.0_wp, "But14diol   ", "30     ", 493.931152_wp, 816.379242_wp), &
      integration_case_type(2.0_wp, 0.0_wp, 1.0_wp, "MB16-43     ", "CH4    ", 243.767818_wp, 354.652628_wp), &
      integration_case_type(3.0_wp, 0.0_wp, 1.0_wp, "Amino20x4   ", "THR_xab", 789.581576_wp, 1694.721566_wp), &
      integration_case_type(3.0_wp, 0.0_wp, 1.0_wp, "But14diol   ", "30     ", 512.762446_wp, 981.787682_wp), &
      integration_case_type(3.0_wp, 0.0_wp, 1.0_wp, "MB16-43     ", "16     ", 690.329407_wp, 1491.585502_wp), &
      integration_case_type(2.0_wp, 0.0_wp, 1.0_wp, "MB16-43     ", "O2     ", 184.603535_wp, 231.274026_wp) &
      ]

contains

   !> Collect all exported unit tests
   subroutine collect_cavity_marchingcubes(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("spherical_isosurface", test_spherical_isosurface), &
         & new_unittest("sphere_radii", test_sphere_radii), &
         & new_unittest("disjoint_spheres", test_disjoint_spheres), &
         & new_unittest("nested_spheres", test_nested_spheres), &
         & new_unittest("no_grid_points", test_no_grid_points), &
         & new_unittest("no_gradient", test_gradient_unavailable), &
         & new_unittest("rejects_bad_spacing", test_rejects_bad_spacing), &
         & new_unittest(case_to_string(cases(1)), test_case_001), &
         & new_unittest(case_to_string(cases(2)), test_case_002), &
         & new_unittest(case_to_string(cases(3)), test_case_003), &
         & new_unittest(case_to_string(cases(4)), test_case_004), &
         & new_unittest(case_to_string(cases(5)), test_case_005), &
         & new_unittest(case_to_string(cases(6)), test_case_006), &
         & new_unittest(case_to_string(cases(7)), test_case_007), &
         & new_unittest(case_to_string(cases(8)), test_case_008), &
         & new_unittest(case_to_string(cases(9)), test_case_009), &
         & new_unittest(case_to_string(cases(10)), test_case_010), &
         & new_unittest(case_to_string(cases(11)), test_case_011), &
         & new_unittest(case_to_string(cases(12)), test_case_012), &
         & new_unittest(case_to_string(cases(13)), test_case_013), &
         & new_unittest(case_to_string(cases(14)), test_case_014), &
         & new_unittest(case_to_string(cases(15)), test_case_015), &
         & new_unittest(case_to_string(cases(16)), test_case_016) &
         & ]

   end subroutine collect_cavity_marchingcubes

   !> Build a marching-cubes cavity over hydrogen centers with custom radii.
   !>
   !> @param[out]   error   Test failure state
   !> @param[inout] ctx     Run context borrowed by the cavity (caller-owned)
   !> @param[out]   mol     Structure carrying one center per entry of `radii`
   !> @param[out]   cav     Constructed cavity
   !> @param[in]    radii   Per-center radii, bohr
   !> @param[in]    xyz     Center coordinates (3, size(radii)), bohr
   !> @param[in]    spacing Finest grid spacing override (optional)
   subroutine build_custom_cavity(error, ctx, mol, cav, radii, xyz, spacing)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), intent(inout), target :: ctx
      type(structure_type), intent(out) :: mol
      type(cavity_type_marchingcubes), allocatable, intent(out) :: cav
      real(wp), intent(in) :: radii(:)
      real(wp), intent(in) :: xyz(:, :)
      real(wp), intent(in), optional :: spacing

      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
      real(wp) :: grid_spacing
      integer :: iat

      grid_spacing = MC_SPACING
      if (present(spacing)) grid_spacing = spacing

      call new_context(ctx, verbosity=0)

      call new(mol, [(1, iat=1, size(radii))], xyz)

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call svdw_template%new()

      allocate (cav)
      call new_cavity_marchingcubes(cav, ctx, radius_model=radius_model, &
         & lsf_model=svdw_template, spacing=grid_spacing, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

   end subroutine build_custom_cavity

   !> Build a single-atom marching-cubes cavity of radius `SPHERE_RADIUS`.
   !>
   !> The SvdW level set of an isolated atom reduces to its one-body term, so its
   !> zero isosurface is exactly the sphere of that radius -- the one geometry with
   !> a closed-form area and volume to check the integrator against.
   !>
   !> @param[out]   error   Test failure state
   !> @param[inout] ctx     Run context borrowed by the cavity (caller-owned)
   !> @param[out]   mol     Single-atom structure
   !> @param[out]   cav     Constructed cavity
   !> @param[in]    spacing Finest grid spacing override (optional)
   subroutine build_sphere_cavity(error, ctx, mol, cav, spacing)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), intent(inout), target :: ctx
      type(structure_type), intent(out) :: mol
      type(cavity_type_marchingcubes), allocatable, intent(out) :: cav
      real(wp), intent(in), optional :: spacing

      real(wp) :: xyz(3, 1)

      xyz(:, 1) = 0.0_wp
      call build_custom_cavity(error, ctx, mol, cav, [SPHERE_RADIUS], xyz, spacing)

   end subroutine build_sphere_cavity

   !> A single atom must integrate to the analytic sphere area and volume
   subroutine test_spherical_isosurface(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_marchingcubes), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavity built here
      type(moist_context_type), target :: ctx
      real(wp) :: area_ref, volume_ref

      call build_sphere_cavity(error, ctx, mol, cav)
      if (allocated(error)) return

      !> The cavity must *borrow* the caller-owned run context, not copy it
      call check(error, associated(cav%ctx, ctx), &
         & more="Marching-cubes cavity does not borrow the caller-owned run context")
      if (allocated(error)) return

      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      area_ref = 4.0_wp*pi*SPHERE_RADIUS**2
      call check(error, cav%total_area/area_ref, 1.0_wp, thr=REL_THR, &
         & more="Single-atom isosurface area does not match the analytic sphere")
      if (allocated(error)) return

      volume_ref = 4.0_wp/3.0_wp*pi*SPHERE_RADIUS**3
      call check(error, cav%total_volume/volume_ref, 1.0_wp, thr=REL_THR, &
         & more="Single-atom isosurface volume does not match the analytic sphere")

   end subroutine test_spherical_isosurface

   !> The integrated totals must follow the analytic r^2 / r^3 scaling
   !>
   !> The adaptive refinement stops on an effectively absolute floor -- its change
   !> test divides by `max(1, |value|)` -- so the relative error grows as the spheres
   !> shrink: at `MC_SPACING` the volume is low by 0.08 % at r = 5 and by 0.8 % at
   !> r = 1.5, crossing the 1 % checked here only near r = 1.4 (2.1 % at r = 1).
   !> The scan therefore starts at 1.5 bohr, below every atomic radius in use.
   subroutine test_sphere_radii(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Radii scanned, bohr
      real(wp), parameter :: radii(*) = [1.5_wp, 2.0_wp, 3.0_wp, 4.5_wp, 6.0_wp]

      type(structure_type) :: mol
      type(cavity_type_marchingcubes), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavity built here
      type(moist_context_type), target :: ctx
      real(wp) :: xyz(3, 1)
      character(len=16) :: label
      integer :: ir

      xyz(:, 1) = 0.0_wp

      do ir = 1, size(radii)
         write (label, '(f6.2)') radii(ir)

         call build_custom_cavity(error, ctx, mol, cav, radii(ir:ir), xyz)
         if (allocated(error)) return

         call cav%update(mol, error=cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if

         call check(error, cav%total_area/(4.0_wp*pi*radii(ir)**2), 1.0_wp, thr=REL_THR, &
            & more="Sphere area does not match the analytic value at r ="//trim(label))
         if (allocated(error)) return

         call check(error, cav%total_volume/(4.0_wp/3.0_wp*pi*radii(ir)**3), 1.0_wp, thr=REL_THR, &
            & more="Sphere volume does not match the analytic value at r ="//trim(label))
         if (allocated(error)) return

         deallocate (cav)
      end do

   end subroutine test_sphere_radii

   !> Well-separated spheres must integrate to the sum of their analytic totals
   subroutine test_disjoint_spheres(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_marchingcubes), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavity built here
      type(moist_context_type), target :: ctx
      real(wp) :: xyz(3, 3)
      real(wp) :: area_ref, volume_ref

      !> Center separations of 13 bohr leave surface-to-surface gaps of at least
      !> 7.5 bohr, far outside the reach of the SvdW blending: what is left is the
      !> pure discretization error of three independent spheres
      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [13.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [0.0_wp, 13.0_wp, 0.0_wp]

      call build_custom_cavity(error, ctx, mol, cav, MULTI_RADII, xyz)
      if (allocated(error)) return

      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call check(error, cav%nsph, size(MULTI_RADII), &
         & more="Disjoint fixture must report one sphere per center")
      if (allocated(error)) return

      area_ref = sum(4.0_wp*pi*MULTI_RADII**2)
      call check(error, cav%total_area/area_ref, 1.0_wp, thr=REL_THR, &
         & more="Disjoint spheres do not integrate to the summed analytic area")
      if (allocated(error)) return

      volume_ref = sum(4.0_wp/3.0_wp*pi*MULTI_RADII**3)
      call check(error, cav%total_volume/volume_ref, 1.0_wp, thr=REL_THR, &
         & more="Disjoint spheres do not integrate to the summed analytic volume")

   end subroutine test_disjoint_spheres

   !> Spheres buried inside a larger one must not contribute anything
   subroutine test_nested_spheres(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_marchingcubes), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavity built here
      type(moist_context_type), target :: ctx
      real(wp) :: radii(4)
      real(wp) :: xyz(3, 4)
      real(wp) :: area_ref, volume_ref

      !> One enclosing sphere plus three small ones fully inside it: the union is
      !> the enclosing sphere alone
      radii = [ENCLOSING_RADIUS, 1.0_wp, 1.5_wp, 2.0_wp]
      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [2.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [0.0_wp, -1.5_wp, 1.0_wp]
      xyz(:, 4) = [-1.0_wp, 1.0_wp, -2.0_wp]

      call build_custom_cavity(error, ctx, mol, cav, radii, xyz)
      if (allocated(error)) return

      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      area_ref = 4.0_wp*pi*ENCLOSING_RADIUS**2
      call check(error, cav%total_area/area_ref, 1.0_wp, thr=REL_THR, &
         & more="Buried spheres must leave the enclosing sphere area unchanged")
      if (allocated(error)) return

      volume_ref = 4.0_wp/3.0_wp*pi*ENCLOSING_RADIUS**3
      call check(error, cav%total_volume/volume_ref, 1.0_wp, thr=REL_THR, &
         & more="Buried spheres must leave the enclosing sphere volume unchanged")

   end subroutine test_nested_spheres

   !> Marching cubes produces no surface discretization, only the two totals
   subroutine test_no_grid_points(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_marchingcubes), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavity built here
      type(moist_context_type), target :: ctx

      call build_sphere_cavity(error, ctx, mol, cav)
      if (allocated(error)) return

      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call check(error, cav%ngrid, 0, &
         & more="Marching cubes must not report grid points")
      if (allocated(error)) return

      call check(error, .not. allocated(cav%xyz), &
         & more="Marching cubes must not allocate grid coordinates")
      if (allocated(error)) return

      call check(error, cav%nsph, 1, &
         & more="Marching cubes must still report the sphere count")

   end subroutine test_no_grid_points

   !> Requesting nuclear derivatives must fail loudly, not return zeros
   subroutine test_gradient_unavailable(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_marchingcubes), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavity built here
      type(moist_context_type), target :: ctx

      call build_sphere_cavity(error, ctx, mol, cav)
      if (allocated(error)) return

      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call cav%get_gradient()

      call check(error, allocated(cav%error), &
         & more="Marching cubes must report that it has no analytic gradients")
      if (allocated(error)) return

      call check(error, .not. allocated(cav%xyz1_rA), &
         & more="Marching cubes must not fabricate nuclear derivatives")

   end subroutine test_gradient_unavailable

   !> A non-positive grid spacing must be rejected at construction
   subroutine test_rejects_bad_spacing(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_marchingcubes), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
      !> Local run context borrowed by the cavity built here
      type(moist_context_type), target :: ctx
      real(wp) :: radii(1)
      real(wp) :: xyz(3, 1)

      call new_context(ctx, verbosity=0)

      xyz(:, 1) = 0.0_wp
      call new(mol, [1], xyz)

      radii = SPHERE_RADIUS
      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call svdw_template%new()

      allocate (cav)
      call new_cavity_marchingcubes(cav, ctx, radius_model=radius_model, &
         & lsf_model=svdw_template, spacing=0.0_wp, error=cavity_error)

      call check(error, allocated(cavity_error), &
         & more="A zero marching-cubes spacing must be rejected")

   end subroutine test_rejects_bad_spacing

   !> Integrate one molecular case and compare against its fine-grid reference
   !>
   !> @param[out] error     Test failure state
   !> @param[in]  case_idx  Index into the `cases` table
   subroutine run_single_case(error, case_idx)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Index into the `cases` table
      integer, intent(in) :: case_idx

      type(structure_type) :: mol
      type(cavity_type_marchingcubes), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
      !> Local run context borrowed by the cavity built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, verbosity=0)

      call get_structure(mol, trim(cases(case_idx)%dataset), trim(cases(case_idx)%structure))

      call svdw_template%new(blend_k=cases(case_idx)%blend_k, &
         & blend_2b=cases(case_idx)%blend_2b, &
         & blend_3b=cases(case_idx)%blend_3b)

      allocate (cav)
      call new_cavity_marchingcubes(cav, ctx, radius_model=default_cpcm_radii(), &
         & lsf_model=svdw_template, spacing=MC_SPACING, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "new_cavity_marchingcubes failed for "// &
            & case_to_string(cases(case_idx))//": "//cavity_error%message)
         return
      end if

      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "cav%update failed for "// &
            & case_to_string(cases(case_idx))//": "//cavity_error%message)
         return
      end if

      call check(error, cav%total_area/cases(case_idx)%mc_area, 1.0_wp, thr=MOL_REL_THR, &
         & more="Area mismatch for "//case_to_string(cases(case_idx)))
      if (allocated(error)) return

      call check(error, cav%total_volume/cases(case_idx)%mc_volume, 1.0_wp, thr=MOL_REL_THR, &
         & more="Volume mismatch for "//case_to_string(cases(case_idx)))

   end subroutine run_single_case

   !> Compact label of one molecular reference case
   !>
   !> @param[in] c    Reference case entry
   !> @return    str  Printable case label without trailing blanks
   pure function case_to_string(c) result(str)

      !> Reference case entry
      type(integration_case_type), intent(in) :: c

      !> Printable case label without trailing blanks
      character(len=:), allocatable :: str

      character(len=8) :: k_str, b_str, g_str

      write (k_str, '(f4.1)') c%blend_k
      write (b_str, '(f4.1)') c%blend_2b
      write (g_str, '(f4.1)') c%blend_3b
      str = trim(c%dataset)//" "//trim(c%structure)//" k="// &
         & k_str(1:4)//" b="//b_str(1:4)//" g="//g_str(1:4)

   end function case_to_string

   subroutine test_case_001(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 1)
   end subroutine test_case_001

   subroutine test_case_002(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 2)
   end subroutine test_case_002

   subroutine test_case_003(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 3)
   end subroutine test_case_003

   subroutine test_case_004(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 4)
   end subroutine test_case_004

   subroutine test_case_005(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 5)
   end subroutine test_case_005

   subroutine test_case_006(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 6)
   end subroutine test_case_006

   subroutine test_case_007(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 7)
   end subroutine test_case_007

   subroutine test_case_008(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 8)
   end subroutine test_case_008

   subroutine test_case_009(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 9)
   end subroutine test_case_009

   subroutine test_case_010(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 10)
   end subroutine test_case_010

   subroutine test_case_011(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 11)
   end subroutine test_case_011

   subroutine test_case_012(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 12)
   end subroutine test_case_012

   subroutine test_case_013(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 13)
   end subroutine test_case_013

   subroutine test_case_014(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 14)
   end subroutine test_case_014

   subroutine test_case_015(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 15)
   end subroutine test_case_015

   subroutine test_case_016(error)
      type(error_type), allocatable, intent(out) :: error
      call run_single_case(error, 16)
   end subroutine test_case_016

end module test_cavity_marchingcubes
