!> Area conservation and intrinsic separation of planar integration patches
module test_cavity_mt_patches
   use mctc_env, only: wp, mctc_error => error_type
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_cavity_mt_mesh, only: mt_mesh_type
   use moist_cavity_mt_patches, only: mt_patch_type, build_mt_patches
   implicit none(type, external)
   private
   public :: collect_cavity_mt_patches

contains

   !> Collect patch tests
   subroutine collect_cavity_mt_patches(testsuite)
      !> Unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("planar_equal_area", test_planar_equal_area), &
                  new_unittest("close_disconnected_sheets", test_close_sheets), &
                  new_unittest("small_component_and_tiny_face", test_small_component), &
                  new_unittest("component_rounding", test_rounding), &
                  new_unittest("invalid_controls", test_invalid_controls), &
                  new_unittest("centroid_projection", test_centroid_projection), &
                  new_unittest("narrow_neck_unequal_components", test_narrow_neck), &
                  new_unittest("folded_connected_sheet", test_folded_sheet), &
                  new_unittest("unbalanceable_fallback", test_unbalanceable), &
                  new_unittest("lloyd_regular_spacing", test_regular_spacing) &
      ]
   end subroutine collect_cavity_mt_patches

   !> Coarse unequal triangles retain their full area in nearly uniform patches
   subroutine test_planar_equal_area(error)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(mt_mesh_type) :: mesh
      type(mt_patch_type) :: patches
      type(mctc_error), allocatable :: patch_error
      real(wp) :: separation
      integer :: p, q

      call make_sheets(mesh, 1, 4.0_wp)
      ! Move a shared corner to make the adjacent subdivision counts unequal.
      mesh%vertices(:, 3) = [1.0_wp, 4.0_wp, 0.0_wp]
      call face_geometry(mesh)
      call build_mt_patches(mesh, 1.0_wp, 0.1_wp, patches, patch_error)
      if (allocated(patch_error)) then
         call test_failed(error, patch_error%message)
         return
      end if
      call check(error, abs(sum(patches%area) - sum(mesh%area)) < 1.e-12_wp)
      if (allocated(error)) return
      call check(error, all(abs(patches%area - 1.0_wp) <= 0.1_wp))
      if (allocated(error)) return
      call check(error,.not. any(patches%area_exception))
      if (allocated(error)) return
      call verify_on_mesh(error, mesh, patches)
      if (allocated(error)) return
      separation = huge(1.0_wp)
      do p = 1, size(patches%area)
         do q = 1, p - 1
            separation = min(separation, norm2(patches%xyz(:, p) - patches%xyz(:, q)))
         end do
      end do
      call check(error, separation > 0.3_wp, more="Patch representatives clustered on a planar surface")
   end subroutine test_planar_equal_area

   !> Nearby parallel sheets retain separate intrinsic components
   subroutine test_close_sheets(error)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(mt_mesh_type) :: mesh
      type(mt_patch_type) :: patches
      type(mctc_error), allocatable :: patch_error
      integer :: c

      call make_sheets(mesh, 2, 2.0_wp)
      call build_mt_patches(mesh, 0.25_wp, 0.1_wp, patches, patch_error)
      if (allocated(patch_error)) then
         call test_failed(error, patch_error%message)
         return
      end if
      call check(error, size(patches%area) == 32)
      if (allocated(error)) return
      call check(error, all(abs(patches%area/0.25_wp - 1.0_wp) <= 0.1_wp))
      if (allocated(error)) return
      do c = 1, 2
         call check(error, count(patches%component == c) == 16)
         if (allocated(error)) return
         call check(error, abs(sum(patches%area, mask=patches%component == c) - 4.0_wp) < 1.e-12_wp)
         if (allocated(error)) return
         call check(error, all(abs(pack(patches%xyz(3, :), patches%component == c) &
                                   - real(c - 1, wp)*1.e-6_wp) < 1.e-14_wp))
         if (allocated(error)) return
      end do
      call verify_on_mesh(error, mesh, patches)
   end subroutine test_close_sheets

   !> A tiny positive-area component remains present and explicitly exceptional
   subroutine test_small_component(error)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(mt_mesh_type) :: mesh
      type(mt_patch_type) :: patches
      type(mctc_error), allocatable :: patch_error
      integer :: c

      call make_sheets(mesh, 2, 2.0_wp)
      mesh%vertices(1:2, 5:8) = mesh%vertices(1:2, 5:8)*1.e-10_wp
      call face_geometry(mesh)
      call build_mt_patches(mesh, 1.0_wp, 0.1_wp, patches, patch_error)
      if (allocated(patch_error)) then
         call test_failed(error, patch_error%message)
         return
      end if
      call check(error, size(patches%area) == 5)
      if (allocated(error)) return
      call check(error, count(patches%area_exception) == 1)
      if (allocated(error)) return
      c = findloc(patches%component, 2, dim=1)
      call check(error, c > 0)
      if (allocated(error)) return
      call check(error, abs(patches%area(c)/4.e-20_wp - 1.0_wp) < 1.e-12_wp)
      if (allocated(error)) return
      call check(error, patches%area_exception(c))
      if (allocated(error)) return
      call check(error, abs(patches%component_target(c)/4.e-20_wp - 1.0_wp) < 1.e-12_wp)
      if (allocated(error)) return
      call verify_on_mesh(error, mesh, patches)
      if (allocated(error)) return
      call patches%clear()
      call check(error,.not. allocated(patches%xyz) .and. .not. allocated(patches%area_exception))
   end subroutine test_small_component

   !> Integer rounding exceptions are distinguished from attainable point counts
   subroutine test_rounding(error)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(mt_mesh_type) :: mesh
      type(mt_patch_type) :: patches
      type(mctc_error), allocatable :: patch_error

      call make_sheets(mesh, 1, sqrt(1.45_wp))
      call build_mt_patches(mesh, 1.0_wp, 0.1_wp, patches, patch_error)
      if (allocated(patch_error)) then
         call test_failed(error, patch_error%message)
         return
      end if
      call check(error, size(patches%area) == 1 .and. all(patches%area_exception))
      if (allocated(error)) return
      call check(error, abs(sum(patches%area) - 1.45_wp) < 1.e-12_wp)
      if (allocated(error)) return
      call build_mt_patches(mesh, 1.0_wp, 0.3_wp, patches, patch_error)
      if (allocated(patch_error)) then
         call test_failed(error, patch_error%message)
         return
      end if
      call check(error, size(patches%area) == 2 .and. .not. any(patches%area_exception))
      if (allocated(error)) return
      call check(error, all(abs(patches%area - 1.0_wp) <= 0.3_wp))
      if (allocated(error)) return
      call verify_on_mesh(error, mesh, patches)
   end subroutine test_rounding

   !> Invalid target areas and tolerances produce recoverable errors
   subroutine test_invalid_controls(error)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(mt_mesh_type) :: mesh
      type(mt_patch_type) :: patches
      type(mctc_error), allocatable :: patch_error
      call make_sheets(mesh, 1, 1.0_wp)
      call build_mt_patches(mesh, 0.0_wp, 0.1_wp, patches, patch_error)
      call check(error, allocated(patch_error))
      if (allocated(error)) return
      call build_mt_patches(mesh, 1.0_wp, 0.0_wp, patches, patch_error)
      call check(error, allocated(patch_error))
      if (allocated(error)) return
      call build_mt_patches(mesh, 1.0_wp, 1.0_wp, patches, patch_error)
      call check(error, allocated(patch_error))
   end subroutine test_invalid_controls

   !> A planar patch uses its exact area centroid across subdivision boundaries
   subroutine test_centroid_projection(error)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(mt_mesh_type) :: mesh
      type(mt_patch_type) :: patches
      type(mctc_error), allocatable :: patch_error
      call make_sheets(mesh, 1, 2.0_wp)
      call build_mt_patches(mesh, 10.0_wp, 0.1_wp, patches, patch_error)
      if (allocated(patch_error)) then
         call test_failed(error, patch_error%message)
         return
      end if
      call check(error, size(patches%area) == 1)
      if (allocated(error)) return
      call check(error, norm2(patches%xyz(:, 1) - [1.0_wp, 1.0_wp, 0.0_wp]) < 1.e-12_wp, &
                 more="Patch representative did not use the planar area centroid")
      if (allocated(error)) return
      call verify_on_mesh(error, mesh, patches)
   end subroutine test_centroid_projection

   !> A thin bridge and unequal disconnected components preserve capacities
   subroutine test_narrow_neck(error)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(mt_mesh_type) :: mesh
      type(mt_patch_type) :: patches
      type(mctc_error), allocatable :: patch_error
      real(wp) :: separation
      integer :: p, q
      call make_narrow_neck(mesh)
      call build_mt_patches(mesh, 1.0_wp, 0.1_wp, patches, patch_error)
      if (allocated(patch_error)) then
         call test_failed(error, patch_error%message)
         return
      end if
      call check(error, all(abs(patches%area - 1.0_wp) <= 0.1_wp))
      if (allocated(error)) return
      call check(error,.not. any(patches%area_exception))
      if (allocated(error)) return
      call check(error, abs(sum(patches%area, mask=patches%component == 1) - 42.5_wp) < 1.e-11_wp)
      if (allocated(error)) return
      call check(error, abs(sum(patches%area, mask=patches%component == 2) - 9.0_wp) < 1.e-11_wp)
      if (allocated(error)) return
      separation = huge(1.0_wp)
      do p = 1, size(patches%area)
         do q = 1, p - 1
            if (patches%component(p) /= patches%component(q)) cycle
            separation = min(separation, norm2(patches%xyz(:, p) - patches%xyz(:, q)))
         end do
      end do
      call check(error, separation > 0.15_wp, more="Narrow-neck representatives nearly coincide")
      if (allocated(error)) return
      call verify_on_mesh(error, mesh, patches)
   end subroutine test_narrow_neck

   !> A nearly closed fold retains separation after intrinsic unfolding
   subroutine test_folded_sheet(error)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(mt_mesh_type) :: mesh
      type(mt_patch_type) :: patches
      type(mctc_error), allocatable :: patch_error
      real(wp), allocatable :: unfolded(:, :)
      real(wp) :: separation
      integer :: p, q
      call make_sheets(mesh, 1, 4.0_wp)
      mesh%vertices(:, 4) = [4.0_wp, 0.0_wp, 1.e-6_wp]
      call face_geometry(mesh)
      call build_mt_patches(mesh, 1.0_wp, 0.1_wp, patches, patch_error)
      if (allocated(patch_error)) then
         call test_failed(error, patch_error%message)
         return
      end if
      call check(error, all(abs(patches%area - 1.0_wp) <= 0.1_wp))
      if (allocated(error)) return
      call check(error, all(patches%component == 1))
      if (allocated(error)) return
      call check(error, abs(sum(patches%area) - sum(mesh%area)) < 1.e-11_wp)
      if (allocated(error)) return
      call verify_on_mesh(error, mesh, patches)
      if (allocated(error)) return
      unfolded = patches%xyz(1:2, :)
      do p = 1, size(patches%area)
         if (patches%face(p) == 2) unfolded(:, p) = patches%xyz([2, 1], p)
      end do
      separation = huge(1.0_wp)
      do p = 1, size(patches%area)
         do q = 1, p - 1
            separation = min(separation, norm2(unfolded(:, p) - unfolded(:, q)))
         end do
      end do
      call check(error, separation > 0.25_wp, more="Connected folded sheet used a through-space shortcut")
   end subroutine test_folded_sheet

   !> Three lobes joined only through a tiny hub keep their best connected cells
   !>
   !> Two connected cells cannot split the lobes evenly: the cell without the
   !> hub lies inside one lobe, so balancing must fall back instead of failing.
   subroutine test_unbalanceable(error)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(mt_mesh_type) :: mesh
      type(mt_patch_type) :: patches
      type(mctc_error), allocatable :: patch_error

      call make_star(mesh)
      call build_mt_patches(mesh, 0.9_wp, 0.1_wp, patches, patch_error)
      if (allocated(patch_error)) then
         call test_failed(error, patch_error%message)
         return
      end if
      call check(error, size(patches%area) == 2 .and. .not. any(patches%area_exception))
      if (allocated(error)) return
      call check(error,.not. patches%balanced, more="Impossible bound reported as balanced")
      if (allocated(error)) return
      call check(error, patches%max_deviation > 0.1_wp, more="Fallback deviation within tolerance")
      if (allocated(error)) return
      call check(error, maxval(abs(patches%area/0.9_wp - 1.0_wp)), patches%max_deviation, thr=1.e-12_wp)
      if (allocated(error)) return
      call check(error, abs(sum(patches%area) - sum(mesh%area)) < 1.e-12_wp)
      if (allocated(error)) return
      call verify_on_mesh(error, mesh, patches)
   end subroutine test_unbalanceable

   !> Centroidal relaxation spaces planar representatives nearly evenly
   subroutine test_regular_spacing(error)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(mt_mesh_type) :: mesh
      type(mt_patch_type) :: patches
      type(mctc_error), allocatable :: patch_error
      real(wp), allocatable :: nearest(:)
      real(wp) :: mean
      integer :: p, q, np

      call make_sheets(mesh, 1, 6.0_wp)
      call build_mt_patches(mesh, 1.0_wp, 0.1_wp, patches, patch_error)
      if (allocated(patch_error)) then
         call test_failed(error, patch_error%message)
         return
      end if
      np = size(patches%area)
      call check(error, np, 36)
      if (allocated(error)) return
      allocate (nearest(np))
      nearest = huge(1.0_wp)
      do p = 1, np
         do q = 1, np
            if (q == p) cycle
            nearest(p) = min(nearest(p), norm2(patches%xyz(:, p) - patches%xyz(:, q)))
         end do
      end do
      ! Unit-area cells: 1 on a square lattice; without relaxation 0.84 and 0.63
      mean = sum(nearest)/real(np, wp)
      call check(error, mean > 0.875_wp, more="Representatives crowd on a planar sheet")
      if (allocated(error)) return
      call check(error, minval(nearest) > 0.69_wp, more="Two representatives nearly coincide")
      if (allocated(error)) return
      call verify_on_mesh(error, mesh, patches)
   end subroutine test_regular_spacing

   !> Three long planar lobes of 0.6 bohr^2 on the edges of a tiny hub triangle
   !>
   !> @param[out] mesh Retained planar geometry
   subroutine make_star(mesh)
      !> Retained geometry
      type(mt_mesh_type), intent(out) :: mesh
      real(wp), parameter :: side = 0.1_wp, reach = 12.0_wp
      real(wp) :: center(3), middle(3), outward(3)
      integer :: k, a, b
      allocate (mesh%vertices(3, 6), mesh%triangles(3, 4), mesh%neighbors(3, 4), &
                mesh%area(4), mesh%centroid(3, 4), mesh%normal(3, 4))
      mesh%vertices(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      mesh%vertices(:, 2) = [side, 0.0_wp, 0.0_wp]
      mesh%vertices(:, 3) = [0.5_wp*side, 0.5_wp*sqrt(3.0_wp)*side, 0.0_wp]
      center = sum(mesh%vertices(:, 1:3), dim=2)/3.0_wp
      mesh%triangles(:, 1) = [1, 2, 3]
      mesh%neighbors = 0
      do k = 1, 3
         ! The hub edge opposite local vertex k carries lobe k + 1
         a = mesh%triangles(mod(k, 3) + 1, 1)
         b = mesh%triangles(mod(k + 1, 3) + 1, 1)
         middle = 0.5_wp*(mesh%vertices(:, a) + mesh%vertices(:, b))
         outward = (middle - center)/norm2(middle - center)
         mesh%vertices(:, 3 + k) = middle + reach*outward
         mesh%triangles(:, k + 1) = [b, a, 3 + k]
         mesh%neighbors(k, 1) = k + 1
         mesh%neighbors(3, k + 1) = 1
      end do
      call face_geometry(mesh)
   end subroutine make_star

   !> Two broad planar lobes joined by a one-cell bridge, plus a smaller square
   !>
   !> @param[out] mesh Retained planar geometry
   subroutine make_narrow_neck(mesh)
      !> Retained geometry
      type(mt_mesh_type), intent(out) :: mesh
      integer, parameter :: nx = 26, ny = 10, nv = (nx + 1)*(ny + 1)
      integer :: triangles(3, 2*nx*ny + 2), ix, iy, v, f, g, k, l, nt, edge(2), other(2)
      allocate (mesh%vertices(3, nv + 4))
      mesh%vertices = 0.0_wp
      do iy = 0, ny
         do ix = 0, nx
            v = 1 + ix + (nx + 1)*iy
            mesh%vertices(:, v) = [0.5_wp*real(ix, wp), 0.5_wp*real(iy, wp), 0.0_wp]
         end do
      end do
      nt = 0
      do iy = 0, ny - 1
         do ix = 0, nx - 1
            if (ix >= 8 .and. ix < 18 .and. iy /= 4) cycle
            v = 1 + ix + (nx + 1)*iy
            triangles(:, nt + 1) = [v, v + 1, v + nx + 2]
            triangles(:, nt + 2) = [v, v + nx + 2, v + nx + 1]
            nt = nt + 2
         end do
      end do
      mesh%vertices(:, nv + 1) = [20.0_wp, 0.0_wp, 0.0_wp]
      mesh%vertices(:, nv + 2) = [23.0_wp, 0.0_wp, 0.0_wp]
      mesh%vertices(:, nv + 3) = [23.0_wp, 3.0_wp, 0.0_wp]
      mesh%vertices(:, nv + 4) = [20.0_wp, 3.0_wp, 0.0_wp]
      triangles(:, nt + 1) = [nv + 1, nv + 2, nv + 3]
      triangles(:, nt + 2) = [nv + 1, nv + 3, nv + 4]
      nt = nt + 2
      mesh%triangles = triangles(:, :nt)
      allocate (mesh%neighbors(3, nt), mesh%area(nt), mesh%centroid(3, nt), mesh%normal(3, nt))
      mesh%neighbors = 0
      do f = 1, nt
         do k = 1, 3
            if (mesh%neighbors(k, f) /= 0) cycle
            edge = mesh%triangles([mod(k, 3) + 1, mod(k + 1, 3) + 1], f)
            do g = f + 1, nt
               do l = 1, 3
                  other = mesh%triangles([mod(l, 3) + 1, mod(l + 1, 3) + 1], g)
                  if (minval(edge) /= minval(other) .or. maxval(edge) /= maxval(other)) cycle
                  mesh%neighbors(k, f) = g
                  mesh%neighbors(l, g) = f
               end do
            end do
         end do
      end do
      call face_geometry(mesh)
   end subroutine make_narrow_neck

   !> Build separate square sheets from two triangles each
   !>
   !> @param[out] mesh Retained geometry
   !> @param[in] count Number of sheets
   !> @param[in] side Square side length, bohr
   subroutine make_sheets(mesh, count, side)
      !> Retained geometry
      type(mt_mesh_type), intent(out) :: mesh
      !> Number of sheets
      integer, intent(in) :: count
      !> Side length, bohr
      real(wp), intent(in) :: side
      integer :: c, v, f
      allocate (mesh%vertices(3, 4*count), mesh%triangles(3, 2*count), mesh%neighbors(3, 2*count), &
                mesh%area(2*count), mesh%centroid(3, 2*count), mesh%normal(3, 2*count))
      mesh%neighbors = 0
      do c = 1, count
         v = 4*(c - 1)
         f = 2*(c - 1)
         mesh%vertices(:, v + 1) = [0.0_wp, 0.0_wp, real(c - 1, wp)*1.e-6_wp]
         mesh%vertices(:, v + 2) = [side, 0.0_wp, real(c - 1, wp)*1.e-6_wp]
         mesh%vertices(:, v + 3) = [side, side, real(c - 1, wp)*1.e-6_wp]
         mesh%vertices(:, v + 4) = [0.0_wp, side, real(c - 1, wp)*1.e-6_wp]
         mesh%triangles(:, f + 1) = [v + 1, v + 2, v + 3]
         mesh%triangles(:, f + 2) = [v + 1, v + 3, v + 4]
         mesh%neighbors(2, f + 1) = f + 2
         mesh%neighbors(3, f + 2) = f + 1
      end do
      call face_geometry(mesh)
   end subroutine make_sheets

   !> Update geometric face quantities after moving fixture vertices
   !>
   !> @param[in,out] mesh Retained geometry
   subroutine face_geometry(mesh)
      !> Retained geometry
      type(mt_mesh_type), intent(inout) :: mesh
      real(wp) :: a(3), b(3), normal(3)
      integer :: f, ids(3)
      do f = 1, size(mesh%area)
         ids = mesh%triangles(:, f)
         a = mesh%vertices(:, ids(2)) - mesh%vertices(:, ids(1))
         b = mesh%vertices(:, ids(3)) - mesh%vertices(:, ids(1))
         normal = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
         mesh%area(f) = 0.5_wp*norm2(normal)
         mesh%normal(:, f) = normal/norm2(normal)
         mesh%centroid(:, f) = sum(mesh%vertices(:, ids), dim=2)/3.0_wp
      end do
      mesh%total_area = sum(mesh%area)
   end subroutine face_geometry

   !> Check that every representative lies inside its recorded planar triangle
   !> @param[in] mesh Retained geometry
   !> @param[in] patches Integration patches
   subroutine verify_on_mesh(error, mesh, patches)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      !> Retained geometry
      type(mt_mesh_type), intent(in) :: mesh
      !> Integration patches
      type(mt_patch_type), intent(in) :: patches
      real(wp) :: a(3), b(3), r(3), aa, ab, bb, ar, br, denominator, x, y
      integer :: p, ids(3), f
      do p = 1, size(patches%area)
         f = patches%face(p)
         ids = mesh%triangles(:, f)
         a = mesh%vertices(:, ids(2)) - mesh%vertices(:, ids(1))
         b = mesh%vertices(:, ids(3)) - mesh%vertices(:, ids(1))
         r = patches%xyz(:, p) - mesh%vertices(:, ids(1))
         aa = dot_product(a, a)
         ab = dot_product(a, b)
         bb = dot_product(b, b)
         ar = dot_product(a, r)
         br = dot_product(b, r)
         denominator = aa*bb - ab*ab
         x = (ar*bb - br*ab)/denominator
         y = (br*aa - ar*ab)/denominator
         call check(error, x >= -1.e-12_wp .and. y >= -1.e-12_wp .and. x + y <= 1.0_wp + 1.e-12_wp)
         if (allocated(error)) return
         call check(error, norm2(r - x*a - y*b) < 1.e-12_wp*max(norm2(a), norm2(b)))
         if (allocated(error)) return
      end do
      call verify_members(error, mesh, patches)
   end subroutine verify_on_mesh

   !> Check that member faces partition the mesh area into face-connected patches
   !>
   !> Connected pieces imply connected member faces, so a patch whose faces
   !> split into separate face-adjacency groups is disconnected.
   !> @param[in] mesh Retained geometry
   !> @param[in] patches Integration patches
   subroutine verify_members(error, mesh, patches)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      !> Retained geometry
      type(mt_mesh_type), intent(in) :: mesh
      !> Integration patches
      type(mt_patch_type), intent(in) :: patches
      real(wp), allocatable :: face_area(:)
      integer, allocatable :: mark(:), queue(:)
      integer :: np, nf, p, m, f, g, k, head, tail, faces
      np = size(patches%area)
      nf = size(mesh%area)
      call check(error, size(patches%member_offset) == np + 1 .and. patches%member_offset(1) == 1 &
                 .and. all(patches%member_offset(2:) > patches%member_offset(:np)), &
                 more="Member offsets are not increasing from one")
      if (allocated(error)) return
      call check(error, patches%member_offset(np + 1) - 1 == size(patches%member_face) &
                 .and. size(patches%member_face) == size(patches%member_area), &
                 more="Member arrays do not match their offsets")
      if (allocated(error)) return
      call check(error, all(patches%member_face >= 1 .and. patches%member_face <= nf) &
                 .and. all(patches%member_area > 0.0_wp), more="Member face or area out of range")
      if (allocated(error)) return
      allocate (face_area(nf), mark(nf), queue(nf))
      face_area = 0.0_wp
      mark = 0
      do p = 1, np
         faces = 0
         do m = patches%member_offset(p), patches%member_offset(p + 1) - 1
            f = patches%member_face(m)
            face_area(f) = face_area(f) + patches%member_area(m)
            if (mark(f) == p) cycle
            mark(f) = p
            faces = faces + 1
         end do
         call check(error, abs(sum(patches%member_area(patches%member_offset(p):patches%member_offset(p + 1) - 1)) &
                               - patches%area(p)) <= 1.e-12_wp*patches%area(p), &
                    more="Member areas do not sum to the patch area")
         if (allocated(error)) return
         call check(error, mark(patches%face(p)) == p, more="Representative face is not a patch member")
         if (allocated(error)) return
         ! Breadth-first search over member faces; visited faces are marked -p
         head = 1
         tail = 1
         queue(1) = patches%face(p)
         mark(queue(1)) = -p
         do while (head <= tail)
            f = queue(head)
            head = head + 1
            do k = 1, 3
               g = mesh%neighbors(k, f)
               if (g <= 0) cycle
               if (mark(g) /= p) cycle
               mark(g) = -p
               tail = tail + 1
               queue(tail) = g
            end do
         end do
         call check(error, tail == faces, more="Patch member faces are disconnected")
         if (allocated(error)) return
      end do
      call check(error, all(abs(face_area - mesh%area) <= 1.e-12_wp*mesh%area), &
                 more="Member areas do not partition each face")
   end subroutine verify_members

end module test_cavity_mt_patches
