!> Analytic geometry, topology and failure checks for marching tetrahedra
module test_cavity_mt_mesh
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, mctc_error => error_type, fatal_error
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_mt_mesh, only: mt_mesh_type, build_mt_mesh
   implicit none(type, external)
   private
   public :: collect_cavity_mt_mesh

   !> Analytic union or spherical shell with the production LSF clone contract
   type, extends(moist_cavity_drop_lsf_svdw_type) :: analytic_lsf_type
      !> Inner void radius, Bohr; zero selects a union of spheres
      real(wp) :: inner = 0.0_wp
      !> Cached analytic value
      real(wp) :: value = 0.0_wp
      !> Radius multiplier to exercise exterior-bound expansion
      real(wp) :: scale = 1.0_wp
      !> Inject a recoverable preparation failure
      logical :: fail = .false.
   contains
      procedure :: prepare => analytic_prepare
      procedure :: f0 => analytic_f0
   end type analytic_lsf_type

contains

   !> Collect geometric extraction tests
   subroutine collect_cavity_mt_mesh(testsuite)
      !> Exported tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("sphere-convergence", test_sphere_convergence), &
                  new_unittest("disconnected-components", test_disconnected), &
                  new_unittest("exact-grid-zero", test_exact_zero), &
                  new_unittest("enclosed-void", test_void), &
                  new_unittest("thread-determinism", test_threads), &
                  new_unittest("exterior-bound-expansion", test_bounds), &
                  new_unittest("disconnected-cfc-underflow", test_cfc_underflow), &
                  new_unittest("invalid-input-and-evaluation", test_errors)]
   end subroutine collect_cavity_mt_mesh

   !> Verify decreasing planar area and volume errors under lattice refinement
   subroutine test_sphere_convergence(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(analytic_lsf_type) :: lsf
      type(mt_mesh_type) :: coarse, fine
      type(mctc_error), allocatable :: mesh_error
      real(wp) :: xyz(3, 1), area, volume
      xyz = 0.0_wp
      call init_analytic(lsf, xyz, [1.0_wp])
      call build_mt_mesh(lsf, xyz, 0.4_wp, coarse, mesh_error, nthreads=1)
      call require_mesh(mesh_error, error)
      if (allocated(error)) return
      call build_mt_mesh(lsf, xyz, 0.2_wp, fine, mesh_error, nthreads=1)
      call require_mesh(mesh_error, error)
      if (allocated(error)) return
      area = 4.0_wp*pi
      volume = area/3.0_wp
      call check(error, abs(fine%total_area - area) < abs(coarse%total_area - area), "area must converge")
      if (allocated(error)) return
      call check(error, abs(fine%total_volume - volume) < abs(coarse%total_volume - volume), "volume must converge")
      if (allocated(error)) return
      call check(error, abs(fine%total_area/area - 1.0_wp) < 0.02_wp, "sphere area within two percent")
      if (allocated(error)) return
      call check(error, abs(fine%total_volume/volume - 1.0_wp) < 0.03_wp, "sphere volume within three percent")
      if (allocated(error)) return
      call check_closed(fine, 1, error)
   end subroutine test_sphere_convergence

   !> Retain both separated components with correct normals and total volume
   subroutine test_disconnected(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(analytic_lsf_type) :: lsf
      type(mt_mesh_type) :: mesh
      type(mctc_error), allocatable :: mesh_error
      real(wp) :: xyz(3, 2), center(3)
      integer :: it, atom
      xyz = 0.0_wp
      xyz(1, :) = [-2.5_wp, 2.5_wp]
      call init_analytic(lsf, xyz, [1.0_wp, 1.0_wp])
      call build_mt_mesh(lsf, xyz, 0.25_wp, mesh, mesh_error, nthreads=1)
      call require_mesh(mesh_error, error)
      if (allocated(error)) return
      call check_closed(mesh, 2, error)
      if (allocated(error)) return
      call check(error, abs(mesh%total_volume/(8.0_wp*pi/3.0_wp) - 1.0_wp) < 0.05_wp, &
                 "volume must include both disconnected spheres")
      if (allocated(error)) return
      do it = 1, size(mesh%area)
         atom = 1
         if (mesh%centroid(1, it) > 0.0_wp) atom = 2
         center = mesh%centroid(:, it) - xyz(:, atom)
         call check(error, dot_product(center, mesh%normal(:, it)) > 0.0_wp, "disconnected face points inward")
         if (allocated(error)) return
      end do
   end subroutine test_disconnected

   !> Binary spacing places sphere poles exactly on shared lattice vertices
   subroutine test_exact_zero(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(analytic_lsf_type) :: lsf
      type(mt_mesh_type) :: mesh
      type(mctc_error), allocatable :: mesh_error
      real(wp) :: xyz(3, 1)
      integer :: it, poles
      xyz = 0.0_wp
      call init_analytic(lsf, xyz, [1.0_wp])
      call build_mt_mesh(lsf, xyz, 0.25_wp, mesh, mesh_error, nthreads=1)
      call require_mesh(mesh_error, error)
      if (allocated(error)) return
      call check_closed(mesh, 1, error)
      if (allocated(error)) return
      poles = count(sum(abs(mesh%vertices - spread([1.0_wp, 0.0_wp, 0.0_wp], 2, &
                                                   size(mesh%vertices, 2))), dim=1) < epsilon(1.0_wp))
      call check(error, poles, 1, "exact pole must have one shared index")
      if (allocated(error)) return
      do it = 1, size(mesh%area)
         call check(error, dot_product(mesh%centroid(:, it), mesh%normal(:, it)) > 0.0_wp, &
                    "sphere normal points inward")
         if (allocated(error)) return
      end do
   end subroutine test_exact_zero

   !> Preserve the inner boundary and subtract its enclosed volume
   subroutine test_void(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(analytic_lsf_type) :: lsf
      type(mt_mesh_type) :: mesh
      type(mctc_error), allocatable :: mesh_error
      real(wp) :: xyz(3, 1), direction
      integer :: it
      xyz = 0.0_wp
      call init_analytic(lsf, xyz, [2.0_wp])
      lsf%inner = 1.0_wp
      call build_mt_mesh(lsf, xyz, 0.25_wp, mesh, mesh_error, nthreads=1)
      call require_mesh(mesh_error, error)
      if (allocated(error)) return
      call check_closed(mesh, 2, error)
      if (allocated(error)) return
      call check(error, abs(mesh%total_volume/(28.0_wp*pi/3.0_wp) - 1.0_wp) < 0.03_wp, &
                 "inner void volume must be subtracted")
      if (allocated(error)) return
      do it = 1, size(mesh%area)
         direction = dot_product(mesh%centroid(:, it), mesh%normal(:, it))
         if (norm2(mesh%centroid(:, it)) < 1.5_wp) direction = -direction
         call check(error, direction > 0.0_wp, "shell boundary normal has wrong sign")
         if (allocated(error)) return
      end do
   end subroutine test_void

   !> Parallel sampling leaves vertex order, connectivity and measures unchanged
   subroutine test_threads(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(analytic_lsf_type) :: lsf
      type(mt_mesh_type) :: serial, threaded
      type(mctc_error), allocatable :: mesh_error
      real(wp) :: xyz(3, 1)
      xyz(:, 1) = [4.0_wp, -2.0_wp, 7.0_wp]
      call init_analytic(lsf, xyz, [1.0_wp])
      call build_mt_mesh(lsf, xyz, 0.25_wp, serial, mesh_error, nthreads=1)
      call require_mesh(mesh_error, error)
      if (allocated(error)) return
      call build_mt_mesh(lsf, xyz, 0.25_wp, threaded, mesh_error, nthreads=2)
      call require_mesh(mesh_error, error)
      if (allocated(error)) return
      call check(error, size(serial%area), size(threaded%area))
      if (allocated(error)) return
      call check(error, size(serial%vertices, 2), size(threaded%vertices, 2))
      if (allocated(error)) return
      call check(error, all(serial%triangles == threaded%triangles) .and. &
                 all(serial%neighbors == threaded%neighbors), "thread count changes topology")
      if (allocated(error)) return
      call check(error, maxval(abs(serial%vertices - threaded%vertices)), 0.0_wp, thr=0.0_wp)
      if (allocated(error)) return
      call check(error, serial%total_volume, threaded%total_volume, thr=0.0_wp)
      if (allocated(error)) return
      call serial%clear()
      call check(error,.not. allocated(serial%vertices) .and. .not. allocated(serial%triangles), "clear retains storage")
   end subroutine test_threads

   !> Expand a radii-based box until its lattice boundary lies outside
   subroutine test_bounds(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(analytic_lsf_type) :: lsf
      type(mt_mesh_type) :: mesh
      type(mctc_error), allocatable :: mesh_error
      real(wp) :: xyz(3, 1)
      xyz = 0.0_wp
      call init_analytic(lsf, xyz, [1.0_wp])
      lsf%scale = 5.0_wp
      call build_mt_mesh(lsf, xyz, 0.5_wp, mesh, mesh_error, nthreads=1)
      call require_mesh(mesh_error, error)
      if (allocated(error)) return
      call check(error, abs(mesh%total_volume/(500.0_wp*pi/3.0_wp) - 1.0_wp) < 0.02_wp, &
                 "expanded lattice must enclose the full surface")
      if (allocated(error)) return
      call check_closed(mesh, 1, error)
   end subroutine test_bounds

   !> Retain separated CFC spheres despite underflowed exterior pseudo-density
   subroutine test_cfc_underflow(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(structure_type) :: mol
      type(mt_mesh_type) :: mesh
      type(mctc_error), allocatable :: mesh_error
      real(wp) :: xyz(3, 2)
      xyz = 0.0_wp
      xyz(1, :) = [-100.0_wp, 100.0_wp]
      call new(mol, [1, 1], xyz)
      call lsf%new()
      call lsf%update(mol, [1.0_wp, 1.0_wp])
      call build_mt_mesh(lsf, xyz, 0.25_wp, mesh, mesh_error, nthreads=2)
      call require_mesh(mesh_error, error)
      if (allocated(error)) return
      call check_closed(mesh, 2, error)
      if (allocated(error)) return
      call check(error, abs(mesh%total_volume/(8.0_wp*pi/3.0_wp) - 1.0_wp) < 0.1_wp, &
                 "underflowed exterior must not remove disconnected CFC spheres")
   end subroutine test_cfc_underflow

   !> Reject invalid parameters and propagate recoverable LSF errors
   subroutine test_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(analytic_lsf_type) :: lsf
      type(mt_mesh_type) :: mesh
      type(mctc_error), allocatable :: mesh_error
      real(wp) :: xyz(3, 1), spacing(5)
      integer :: i
      xyz = 0.0_wp
      call init_analytic(lsf, xyz, [1.0_wp])
      spacing = [0.0_wp, -1.0_wp, tiny(1.0_wp), huge(1.0_wp), ieee_value(0.0_wp, ieee_quiet_nan)]
      do i = 1, size(spacing)
         call build_mt_mesh(lsf, xyz, spacing(i), mesh, mesh_error, nthreads=1)
         call check(error, allocated(mesh_error), "invalid spacing must return an error")
         if (allocated(error)) return
         call check(error,.not. allocated(mesh%vertices), "failed mesh retains geometry")
         if (allocated(error)) return
      end do
      call build_mt_mesh(lsf, xyz, 0.25_wp, mesh, mesh_error, nthreads=0)
      call check(error, allocated(mesh_error), "invalid thread count must return an error")
      if (allocated(error)) return
      lsf%fail = .true.
      call build_mt_mesh(lsf, xyz, 0.25_wp, mesh, mesh_error, nthreads=2)
      call check(error, allocated(mesh_error), "LSF failure must propagate")
   end subroutine test_errors

   !> Check reciprocal adjacency, positive areas and connected components
   !>
   !> @param[in] mesh Surface geometry
   !> @param[in] expected Expected boundary components
   !> @param[out] error Test failure
   subroutine check_closed(mesh, expected, error)
      !> Surface geometry
      type(mt_mesh_type), intent(in) :: mesh
      !> Expected boundary components
      integer, intent(in) :: expected
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      logical, allocatable :: visited(:)
      integer, allocatable :: queue(:)
      integer :: it, k, neighbor, components, head, tail, jt, side, a, b, c, d
      call check(error, all(mesh%area > 0.0_wp), "zero or negative triangle area")
      if (allocated(error)) return
      call check(error, all(mesh%neighbors > 0) .and. all(mesh%neighbors <= size(mesh%area)), "open or invalid adjacency")
      if (allocated(error)) return
      do it = 1, size(mesh%area)
         do k = 1, 3
            neighbor = mesh%neighbors(k, it)
            side = findloc(mesh%neighbors(:, neighbor), it, dim=1)
            call check(error, side > 0, "nonreciprocal adjacency")
            if (allocated(error)) return
            a = mesh%triangles(mod(k, 3) + 1, it)
            b = mesh%triangles(mod(k + 1, 3) + 1, it)
            c = mesh%triangles(mod(side, 3) + 1, neighbor)
            d = mesh%triangles(mod(side + 1, 3) + 1, neighbor)
            call check(error, a == d .and. b == c, "adjacent faces must share oppositely directed edges")
            if (allocated(error)) return
         end do
      end do
      allocate (visited(size(mesh%area)), queue(size(mesh%area)))
      visited = .false.
      components = 0
      do it = 1, size(mesh%area)
         if (visited(it)) cycle
         components = components + 1
         head = 1
         tail = 1
         queue(1) = it
         visited(it) = .true.
         do while (head <= tail)
            jt = queue(head)
            head = head + 1
            do k = 1, 3
               neighbor = mesh%neighbors(k, jt)
               if (visited(neighbor)) cycle
               visited(neighbor) = .true.
               tail = tail + 1
               queue(tail) = neighbor
            end do
         end do
      end do
      call check(error, components, expected, "wrong boundary component count")
   end subroutine check_closed

   !> Bind an analytic fixture through the production LSF update API
   !>
   !> @param[out] lsf Analytic LSF fixture
   !> @param[in] xyz Atomic centers (3, nat), Bohr
   !> @param[in] radii Sphere radii (nat), Bohr
   subroutine init_analytic(lsf, xyz, radii)
      !> Analytic fixture
      type(analytic_lsf_type), intent(out) :: lsf
      !> Atomic centers
      real(wp), intent(in) :: xyz(:, :)
      !> Sphere radii
      real(wp), intent(in) :: radii(:)
      type(structure_type) :: mol
      integer :: numbers(size(radii))
      numbers = 1
      call new(mol, numbers, xyz)
      call lsf%new()
      call lsf%update(mol, radii)
   end subroutine init_analytic

   !> Cache an exact union distance or spherical-shell level set
   !>
   !> @param[in,out] self Analytic fixture
   !> @param[in] point Evaluation point (3), Bohr
   !> @param[out] error Injected evaluation failure
   subroutine analytic_prepare(self, point, error)
      !> Analytic fixture
      class(analytic_lsf_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Evaluation failure
      type(mctc_error), allocatable, intent(out) :: error
      integer :: i
      if (self%fail) then
         call fatal_error(error, "analytic fixture preparation failure")
         return
      end if
      self%value = huge(1.0_wp)
      do i = 1, self%ncenters
         self%value = min(self%value, norm2(point - self%mol%xyz(:, i)) - self%scale*self%radii(i))
      end do
      if (self%inner > 0.0_wp) self%value = max(self%value, self%inner - norm2(point - self%mol%xyz(:, 1)))
   end subroutine analytic_prepare

   !> Return the cached analytic value
   !>
   !> @param[in] self Analytic fixture
   !> @param[out] val Level-set value
   subroutine analytic_f0(self, val)
      !> Analytic fixture
      class(analytic_lsf_type), intent(in) :: self
      !> Level-set value
      real(wp), intent(out) :: val
      val = self%value
   end subroutine analytic_f0

   !> Convert a geometry failure into a test failure
   !>
   !> @param[in] mesh_error Optional geometry failure
   !> @param[out] error Test failure
   subroutine require_mesh(mesh_error, error)
      !> Optional geometry failure
      type(mctc_error), allocatable, intent(in) :: mesh_error
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      if (allocated(mesh_error)) call test_failed(error, mesh_error%message)
   end subroutine require_mesh

end module test_cavity_mt_mesh
