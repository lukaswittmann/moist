!> Gaussian width factors that make the MT quadrature exact for a uniform layer
!>
!> On the unswitched retained surface, with geometric patch areas A, each point
!> i solves one equation for its width xi_i:
!>
!>   sqrt(2/pi) xi_i A_i + sum_(j /= i) erf(xi_ij d_ij)/d_ij A_j = integral 1/|x_i - y| dA
!>
!> The right-hand side is the potential of the unit layer on the retained
!> planar mesh. Pair terms hardly depend on the widths, so one diagonal Newton
!> step from xi = scale/sqrt(A) solves the system. Both sides are restricted to
!> the same smooth window around x_i, since quadrature and integral agree in
!> the far field.
!>
!> The resulting factors correct the patch arrangement only. Applied to the
!> switched widths, they keep the larger self terms by which the critical
!> switch suppresses charge at low-gradient points; a layer consistent with the
!> switched areas would cancel that suppression.
module moist_cavity_mt_widths
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_cavity_mt_mesh, only: mt_mesh_type
   use moist_cavity_mt_patches, only: mt_patch_type
   use moist_math_adjacency_list, only: adjacency_list_type
   implicit none(type, external)
   private
   public :: width_factors

   !> sqrt(pi/2), converting a self potential to a Gaussian width
   real(wp), parameter :: width_per_self = 1.2533141373155003_wp
   !> Largest relative change of one width
   real(wp), parameter :: max_ratio = 1.25_wp
   !> Triangles farther than this many diameters use the centroid rule
   real(wp), parameter :: far_factor = 8.0_wp
   !> Subdivision depth limit for close off-plane triangles
   integer, parameter :: max_depth = 9

contains

   !> One Newton step on the unswitched widths against the uniform-layer potential
   !>
   !> Factors are limited to `max_ratio` in either direction; limited points are
   !> counted.
   !>
   !> @param[in] mesh Retained planar mesh
   !> @param[in] patches Integration patches with face members
   !> @param[in] patch_index Patch of each grid point
   !> @param[in] xyz Grid points (3, n), bohr
   !> @param[in] area Geometric patch areas of the points (n), bohr^2
   !> @param[in] scale Initial width coefficient, xi = scale/sqrt(area)
   !> @param[in] radius Window radius; the window vanishes beyond twice this, bohr
   !> @param[in] nthreads Thread count
   !> @param[out] factor Width correction of each point (n)
   !> @param[out] limited Number of limited corrections
   !> @param[out] error Allocation or evaluation failure
   subroutine width_factors(mesh, patches, patch_index, xyz, area, scale, radius, nthreads, factor, limited, error)
      !> Retained planar mesh
      type(mt_mesh_type), intent(in) :: mesh
      !> Integration patches
      type(mt_patch_type), intent(in) :: patches
      !> Patch of each point
      integer, intent(in) :: patch_index(:)
      !> Grid points
      real(wp), intent(in) :: xyz(:, :)
      !> Geometric areas
      real(wp), intent(in) :: area(:)
      !> Initial width coefficient
      real(wp), intent(in) :: scale
      !> Window radius
      real(wp), intent(in) :: radius
      !> Thread count
      integer, intent(in) :: nthreads
      !> Width corrections
      real(wp), intent(out) :: factor(:)
      !> Number of limited corrections
      integer, intent(out) :: limited
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: neighbours
      real(wp), allocatable :: xi(:)
      real(wp) :: reach, extent, layer, pairs, d, xij
      integer :: n, i, j, k, p, m, stat

      n = size(xyz, 2)
      limited = 0
      factor = 1.0_wp
      if (n == 0) return
      allocate (xi(n), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate marching tetrahedra width correction")
         return
      end if
      ! Patches reach their window faces from at most their own extent away
      extent = 0.0_wp
      do i = 1, n
         p = patch_index(i)
         do m = patches%member_offset(p), patches%member_offset(p + 1) - 1
            extent = max(extent, norm2(mesh%centroid(:, patches%member_face(m)) - xyz(:, i)))
         end do
      end do
      xi = scale/sqrt(area)
      reach = 2.0_wp*radius + extent
      call neighbours%init(reach)
      call neighbours%update(xyz)
!$omp parallel do default(none) schedule(dynamic, 16) num_threads(nthreads) &
!$omp shared(mesh, patches, patch_index, xyz, area, xi, radius, neighbours, factor, n) &
!$omp private(i, j, k, layer, pairs, d, xij) reduction(+:limited)
      do i = 1, n
         layer = patch_layer(mesh, patches, patch_index(i), xyz(:, i), radius)
         pairs = 0.0_wp
         do k = neighbours%inl(i) + 1, neighbours%inl(i) + neighbours%nnl(i)
            j = neighbours%nlat(k)
            layer = layer + patch_layer(mesh, patches, patch_index(j), xyz(:, i), radius)
            d = neighbours%dist(k)
            xij = xi(i)*xi(j)/sqrt(xi(i)**2 + xi(j)**2)
            pairs = pairs + window(d, radius)*erf(xij*d)/d*area(j)
         end do
         ! Self potential the point needs, as a fraction of its initial one
         factor(i) = width_per_self*(layer - pairs)/area(i)/xi(i)
         ! Nonfinite factors are evaluation failures, reported below
         if (.not. ieee_is_finite(factor(i))) cycle
         if (factor(i) < 1.0_wp/max_ratio .or. factor(i) > max_ratio) then
            limited = limited + 1
            factor(i) = merge(max_ratio, 1.0_wp/max_ratio, factor(i) > max_ratio)
         end if
      end do
!$omp end parallel do
      if (.not. all(ieee_is_finite(factor))) then
         call fatal_error(error, "Nonfinite marching tetrahedra width correction")
      end if
   end subroutine width_factors

   !> Windowed potential at x of the unit layer on one patch
   !>
   !> Faces split between patches contribute in proportion to their member areas.
   !>
   !> @param[in] mesh Retained planar mesh
   !> @param[in] patches Integration patches
   !> @param[in] p Patch
   !> @param[in] x Observation point, bohr
   !> @param[in] radius Window radius, bohr
   function patch_layer(mesh, patches, p, x, radius) result(layer)
      !> Retained planar mesh
      type(mt_mesh_type), intent(in) :: mesh
      !> Integration patches
      type(mt_patch_type), intent(in) :: patches
      !> Patch
      integer, intent(in) :: p
      !> Observation point
      real(wp), intent(in) :: x(3)
      !> Window radius
      real(wp), intent(in) :: radius
      !> Windowed layer potential, bohr
      real(wp) :: layer
      real(wp) :: w
      integer :: m, f
      layer = 0.0_wp
      do m = patches%member_offset(p), patches%member_offset(p + 1) - 1
         f = patches%member_face(m)
         w = window(norm2(mesh%centroid(:, f) - x), radius)
         if (w <= 0.0_wp) cycle
         layer = layer + w*patches%member_area(m)/mesh%area(f)*triangle_potential( &
            & mesh%vertices(:, mesh%triangles(1, f)), mesh%vertices(:, mesh%triangles(2, f)), &
            & mesh%vertices(:, mesh%triangles(3, f)), x, 0)
      end do
   end function patch_layer

   !> Smooth cutoff, one inside the radius and zero beyond twice the radius
   !>
   !> @param[in] r Distance, bohr
   !> @param[in] radius Window radius, bohr
   pure function window(r, radius) result(w)
      !> Distance
      real(wp), intent(in) :: r
      !> Window radius
      real(wp), intent(in) :: radius
      !> Window value
      real(wp) :: w
      real(wp) :: s
      s = (r - radius)/radius
      w = 1.0_wp
      if (s <= 0.0_wp) return
      w = 0.0_wp
      if (s >= 1.0_wp) return
      w = 1.0_wp - s*s*(3.0_wp - 2.0_wp*s)
   end function window

   !> Integral of 1/|x - y| over a planar triangle
   !>
   !> Far triangles use the centroid rule. Close triangles containing the
   !> projection plane of x use the analytic in-plane formula; close off-plane
   !> triangles are subdivided until either applies.
   !>
   !> @param[in] a First vertex, bohr
   !> @param[in] b Second vertex, bohr
   !> @param[in] c Third vertex, bohr
   !> @param[in] x Observation point, bohr
   !> @param[in] depth Subdivision level
   pure recursive function triangle_potential(a, b, c, x, depth) result(potential)
      !> First vertex
      real(wp), intent(in) :: a(3)
      !> Second vertex
      real(wp), intent(in) :: b(3)
      !> Third vertex
      real(wp), intent(in) :: c(3)
      !> Observation point
      real(wp), intent(in) :: x(3)
      !> Subdivision level
      integer, intent(in) :: depth
      !> Integral, bohr
      real(wp) :: potential
      real(wp) :: normal(3), centroid(3), ab(3), bc(3), ca(3), area, diameter, height
      normal = cross(b - a, c - a)
      area = 0.5_wp*norm2(normal)
      centroid = (a + b + c)/3.0_wp
      diameter = max(norm2(b - a), norm2(c - b), norm2(a - c))
      ! Slivers from surfaces through lattice nodes have no reliable plane
      potential = 0.0_wp
      if (area <= sqrt(epsilon(1.0_wp))*diameter**2) return
      if (norm2(x - centroid) > far_factor*diameter) then
         potential = area/norm2(x - centroid)
         return
      end if
      normal = normal/(2.0_wp*area)
      height = dot_product(x - a, normal)
      if (abs(height) <= 1.0e-3_wp*diameter .or. depth >= max_depth) then
         ! Neglecting a height h changes the integral by O((h/diameter)**2)
         potential = abs(edge_fan(a, b, x - height*normal, normal) + edge_fan(b, c, x - height*normal, normal) &
                         + edge_fan(c, a, x - height*normal, normal))
         return
      end if
      ab = 0.5_wp*(a + b)
      bc = 0.5_wp*(b + c)
      ca = 0.5_wp*(c + a)
      potential = triangle_potential(a, ab, ca, x, depth + 1) + triangle_potential(ab, b, bc, x, depth + 1) &
                  + triangle_potential(ca, bc, c, x, depth + 1) + triangle_potential(ab, bc, ca, x, depth + 1)
   end function triangle_potential

   !> Signed integral of 1/|x - y| over the triangle (x, a, b) in its plane
   !>
   !> In polar coordinates around x the integrand reduces to the distance of
   !> the edge, d/cos(phi), whose antiderivative is d asinh(tan(phi)).
   !>
   !> @param[in] a Edge start, bohr
   !> @param[in] b Edge end, bohr
   !> @param[in] x In-plane observation point, bohr
   !> @param[in] normal Unit plane normal fixing the orientation
   pure function edge_fan(a, b, x, normal) result(value)
      !> Edge start
      real(wp), intent(in) :: a(3)
      !> Edge end
      real(wp), intent(in) :: b(3)
      !> Observation point
      real(wp), intent(in) :: x(3)
      !> Plane normal
      real(wp), intent(in) :: normal(3)
      !> Signed integral, bohr
      real(wp) :: value
      real(wp) :: tangent(3), length, ta, d
      value = 0.0_wp
      length = norm2(b - a)
      if (length <= 0.0_wp) return
      tangent = (b - a)/length
      ta = dot_product(a - x, tangent)
      d = norm2((a - x) - ta*tangent)
      if (d <= 1.0e-14_wp*length) return
      value = sign(d*(asinh((ta + length)/d) - asinh(ta/d)), dot_product(cross(a - x, b - x), normal))
   end function edge_fan

   !> Cross product
   !>
   !> @param[in] u First vector
   !> @param[in] v Second vector
   pure function cross(u, v) result(w)
      !> First vector
      real(wp), intent(in) :: u(3)
      !> Second vector
      real(wp), intent(in) :: v(3)
      !> Cross product
      real(wp) :: w(3)
      w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]
   end function cross

end module moist_cavity_mt_widths
