!> Cavity contracts and electrostatic accuracy of the MT surface quadrature
module test_cavity_mt
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   use moist_radii, only: radius_type, new_radii_custom_atoms
   use moist_cavity_mt, only: cavity_type_mt, &
      & new_cavity_mt, moist_cavity_mt_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_switching, only: moist_cavity_drop_swif_sigmoid_bump_type, new_swif_sigmoid_bump
   use moist_channels_fields, only: field_query_type
   use pcmgen_gauss_conductor, only: conductor_state, conductor_assemble, conductor_medium, conductor_solve
   implicit none(type, external)
   private
   public :: collect_cavity_mt

contains

   !> Collect cavity and CPCM checks
   !>
   !> @param[out] testsuite Test collection
   subroutine collect_cavity_mt(testsuite)
      !> Test collection
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("sphere_grid", test_sphere_grid), &
                  new_unittest("seam_switch_integration", test_seam_switch), &
                  new_unittest("fully_suppressed", test_fully_suppressed), &
                  new_unittest("frozen_geometry", test_frozen_geometry), &
                  new_unittest("invalid_parameters", test_invalid_parameters), &
                  new_unittest("sphere_cpcm", test_sphere_cpcm), &
                  new_unittest("sliver_widths", test_sliver_widths)]
   end subroutine collect_cavity_mt

   !> Build a custom-radius SvdW fixture with an explicit thread count
   !> @param[in,out] ctx Borrowed run context
   !> @param[out] cav Constructed cavity
   !> @param[in] xyz Atomic coordinates, bohr
   !> @param[in] radii Atomic radii, bohr
   !> @param[in] param Optional settings
   subroutine fixture(error, ctx, cav, xyz, radii, param)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Borrowed run context
      type(moist_context_type), intent(inout), target :: ctx
      !> Constructed cavity
      type(cavity_type_mt), intent(out) :: cav
      !> Atomic coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Atomic radii
      real(wp), intent(in) :: radii(:)
      !> Optional settings
      type(moist_cavity_mt_parameters_type), intent(in), optional :: param
      type(structure_type) :: mol
      class(radius_type), allocatable :: radius_model
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(mctc_error), allocatable :: local_error
      integer :: i

      call new_context(ctx, verbosity=0, nthreads=1)
      call new(mol, [(1, i=1, size(radii))], xyz)
      call new_radii_custom_atoms(radii, radius_model, local_error)
      if (allocated(local_error)) then
         call test_failed(error, local_error%message)
         return
      end if
      call lsf%new()
      call new_cavity_mt(cav, radius_model, lsf, local_error, param, ctx)
      if (.not. allocated(local_error)) call cav%update(mol, local_error)
      if (allocated(local_error)) call test_failed(error, local_error%message)
   end subroutine fixture

   !> Verify the complete common surface contract on an analytic sphere
   subroutine test_sphere_grid(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(field_query_type) :: fields
      real(wp) :: grad(3), xyz(3, 1), radial, d
      real(wp), allocatable :: layer(:)
      type(mctc_error), allocatable :: local_error
      integer :: i, j

      xyz = 0.0_wp
      call fixture(error, ctx, cav, xyz, [3.0_wp])
      if (allocated(error)) return
      call check(error, cav%built .and. cav%updated .and. cav%ngrid > 50)
      if (allocated(error)) return
      call check(error, abs(cav%mesh_area/(36.0_wp*pi) - 1.0_wp) < 0.005_wp)
      if (allocated(error)) return
      call check(error, abs(cav%mesh_volume/(36.0_wp*pi) - 1.0_wp) < 0.008_wp)
      if (allocated(error)) return
      call check(error, cav%total_area, cav%mesh_area, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(cav%patch_area/cav%param%target_area - 1.0_wp)) <= cav%param%area_tolerance)
      if (allocated(error)) return
      ! The width correction stays small on a well-spaced unswitched surface
      call check(error, all(abs(cav%xi0*sqrt(cav%a)/cav%param%gaussian_scale - 1.0_wp) < 0.15_wp), &
                 more="Width correction is unexpectedly large")
      if (allocated(error)) return
      ! A uniform layer has the potential 4 pi R everywhere on a sphere
      allocate (layer(cav%ngrid))
      do i = 1, cav%ngrid
         layer(i) = sqrt(2.0_wp/pi)*cav%xi0(i)*cav%a(i)
         do j = 1, cav%ngrid
            if (j == i) cycle
            d = norm2(cav%xyz(:, i) - cav%xyz(:, j))
            layer(i) = layer(i) + erf(cav%xi0(i)*cav%xi0(j)/sqrt(cav%xi0(i)**2 + cav%xi0(j)**2)*d)/d*cav%a(j)
         end do
      end do
      ! Corrected 0.08%, uncorrected 0.49%
      call check(error, maxval(abs(layer/(12.0_wp*pi) - 1.0_wp)) < 2.0e-3_wp, &
                 more="Uniform layer potential varies over the sphere")
      if (allocated(error)) return
      call check(error, all(cav%owner == 1) .and. all(abs(cav%f - 1.0_wp) < epsilon(1.0_wp)))
      if (allocated(error)) return
      call check(error, cav%asph(1), cav%total_area, thr=1.0e-11_wp)
      if (allocated(error)) return
      call cav%lsf_model%set_max_deriv(1)
      do i = 1, cav%ngrid
         radial = norm2(cav%xyz(:, i))
         call check(error, radial <= 3.0_wp + 1.0e-12_wp .and. radial > 2.98_wp)
         if (allocated(error)) return
         call cav%lsf_model%prepare(cav%xyz(:, i), local_error)
         if (allocated(local_error)) then
            call test_failed(error, local_error%message)
            return
         end if
         call cav%lsf_model%f012_r(lsf1_r=grad)
         call check(error, maxval(abs(cav%normal0(:, i) - grad/norm2(grad))) < 1.0e-13_wp)
         if (allocated(error)) return
      end do
      call cav%list_fields(fields)
      call check(error, fields%index_of("mesh_area") > 0 .and. fields%index_of("critical_weight") > 0)
   end subroutine test_sphere_grid

   !> Seam switching is integrated per face, independent of the patch partition
   subroutine test_seam_switch(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: fine, coarse
      type(moist_cavity_mt_parameters_type) :: param
      type(moist_cavity_drop_swif_sigmoid_bump_type) :: bump
      real(wp) :: xyz(3, 2), area
      integer :: i

      xyz = 0.0_wp
      xyz(1, :) = [-2.05_wp, 2.05_wp]
      param%spacing = 0.25_wp
      param%target_area = 0.5_wp
      call fixture(error, ctx, fine, xyz, [2.0_wp, 2.0_wp], param)
      if (allocated(error)) return
      param%target_area = 2.0_wp
      call fixture(error, ctx, coarse, xyz, [2.0_wp, 2.0_wp], param)
      if (allocated(error)) return
      call check(error, count(fine%critical_weight < 0.99_wp) > 0, more="Fixture has no switched seam")
      if (allocated(error)) return
      call check(error, fine%total_area < fine%mesh_area, more="Seam switch removed no area")
      if (allocated(error)) return
      call check(error, abs(fine%total_area/coarse%total_area - 1.0_wp) < 1.0e-10_wp, &
         & more="Switched area depends on the patch partition")
      if (allocated(error)) return
      call new_swif_sigmoid_bump(bump, param%w_0ls_from, param%w_0ls_to, &
         & p_hi=param%w_0ls_p, p_lo=param%w_0ls_p, a_hi=param%w_0ls_a, a_lo=param%w_0ls_a)
      do i = 1, coarse%ngrid
         call switched_area(error, coarse, bump, coarse%patch_index(i), area)
         if (allocated(error)) return
         call check(error, coarse%a(i), area, thr=1.0e-12_wp)
         if (allocated(error)) return
      end do
   end subroutine test_seam_switch

   !> Reference switched area of one patch from its face members
   !> @param[in,out] cav Built cavity
   !> @param[in] bump Critical switch
   !> @param[in] patch Index in the complete patch collection
   !> @param[out] area Sum of face areas times the switch at face centroids
   subroutine switched_area(error, cav, bump, patch, area)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Built cavity
      type(cavity_type_mt), intent(inout) :: cav
      !> Critical switch
      type(moist_cavity_drop_swif_sigmoid_bump_type), intent(in) :: bump
      !> Patch index
      integer, intent(in) :: patch
      !> Switched area
      real(wp), intent(out) :: area
      type(mctc_error), allocatable :: local_error
      real(wp) :: grad(3)
      integer :: j, f

      area = 0.0_wp
      call cav%lsf_model%set_max_deriv(1)
      do j = cav%patches%member_offset(patch), cav%patches%member_offset(patch + 1) - 1
         f = cav%patches%member_face(j)
         call cav%lsf_model%prepare(cav%mesh%centroid(:, f), local_error)
         if (allocated(local_error)) then
            call test_failed(error, local_error%message)
            return
         end if
         call cav%lsf_model%f012_r(lsf1_r=grad)
         area = area + bump%f0(norm2(grad))*cav%patches%member_area(j)
      end do
   end subroutine switched_area

   !> Suppression removes quadrature points while preserving the full mesh
   subroutine test_fully_suppressed(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(moist_cavity_mt_parameters_type) :: param

      param%spacing = 0.35_wp
      param%target_area = 2.0_wp
      param%w_0ls_from = 10.0_wp
      param%w_0ls_to = 20.0_wp
      call fixture(error, ctx, cav, reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [3.0_wp], param)
      if (allocated(error)) return
      call check(error, cav%ngrid == 0 .and. cav%built .and. size(cav%mesh%area) > 0)
      if (allocated(error)) return
      call check(error, abs(cav%total_area) < tiny(1.0_wp) .and. abs(sum(cav%asph)) < tiny(1.0_wp))
      if (allocated(error)) return
      call check(error, sum(cav%patches%area), cav%mesh_area, thr=1.0e-10_wp)
   end subroutine test_fully_suppressed

   !> Frozen geometry accepts normals and rejects unsupported derivative requests
   subroutine test_frozen_geometry(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_mt) :: cav
      type(mctc_error), allocatable :: local_error

      call cav%configure(normal=.true., error=local_error)
      call check(error,.not. allocated(local_error))
      if (allocated(error)) return
      call check(error,.not. cav%has_field_dependent_geometry())
      if (allocated(error)) return
      call cav%get_gradient(local_error)
      call check(error, allocated(local_error))
      if (allocated(error)) return
      call cav%configure(curvature=.true., error=local_error)
      call check(error, allocated(local_error))
   end subroutine test_frozen_geometry

   !> Reject nonfinite and invalid native settings
   subroutine test_invalid_parameters(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_mt_parameters_type) :: param
      type(mctc_error), allocatable :: local_error

      param%spacing = ieee_value(0.0_wp, ieee_quiet_nan)
      call param%validate(local_error)
      call check(error, allocated(local_error))
      if (allocated(error)) return
      call param%init_defaults()
      param%w_0ls_to = param%w_0ls_from
      call param%validate(local_error)
      call check(error, allocated(local_error))
      if (allocated(error)) return
      call param%init_defaults()
      param%target_area = 0.0_wp
      call param%validate(local_error)
      call check(error, allocated(local_error))
   end subroutine test_invalid_parameters

   !> Check centered and off-center conductor energies and local charge density
   subroutine test_sphere_cpcm(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(moist_cavity_mt_parameters_type) :: param
      type(conductor_state) :: state
      type(mctc_error), allocatable :: local_error
      real(wp), allocatable :: v(:), q(:), gv(:), sigma(:)
      real(wp) :: shift(3), r, energy, reference, density_error
      integer :: i, k

      param%spacing = 0.15_wp
      param%target_area = 0.6_wp
      call fixture(error, ctx, cav, reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [3.0_wp], param)
      if (allocated(error)) return
      allocate (v(cav%ngrid), q(cav%ngrid), gv(cav%ngrid), sigma(cav%ngrid))
      call conductor_assemble(state, cav%xyz, cav%xi0, cav%f, local_error)
      if (.not. allocated(local_error)) call conductor_medium(state, 1.0_wp, local_error)
      if (allocated(local_error)) then
         call test_failed(error, local_error%message)
         return
      end if
      do k = 0, 2
         shift = real(k, wp)*[0.3_wp, 0.4_wp, 0.2_wp]
         reference = -3.0_wp/(2.0_wp*(9.0_wp - sum(shift**2)))
         do i = 1, cav%ngrid
            r = norm2(cav%xyz(:, i) - shift)
            v(i) = erf(cav%xi0(i)*r)/r
            sigma(i) = -(9.0_wp - sum(shift**2))/(12.0_wp*pi*r**3)
         end do
         call conductor_solve(state, v, energy, q, gv, local_error)
         if (allocated(local_error)) then
            call test_failed(error, local_error%message)
            return
         end if
         ! Errors with gfortran: 0.05% energy, 0.015% charge, 0.38% density;
         ! without the width correction 0.07%, 0.08% and 2.7%
         call check(error, abs(energy/reference - 1.0_wp) < 1.5e-3_wp, more="MT sphere conductor energy")
         if (allocated(error)) return
         call check(error, abs(sum(q) + 1.0_wp) < 5.e-4_wp, more="MT induced total charge")
         if (allocated(error)) return
         density_error = sqrt(sum(cav%a*(q/cav%a - sigma)**2)/sum(cav%a*sigma**2))
         call check(error, density_error < 0.01_wp, more="MT local charge density")
         if (allocated(error)) return
      end do
   end subroutine test_sphere_cpcm

   !> Near-zero-area slivers from a surface through lattice nodes leave the widths finite
   subroutine test_sliver_widths(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(moist_cavity_mt_parameters_type) :: param
      real(wp) :: diameter
      integer :: f, slivers, ids(3)

      param%target_area = 0.4_wp
      call fixture(error, ctx, cav, reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [3.0_wp], param)
      if (allocated(error)) return
      slivers = 0
      do f = 1, size(cav%mesh%area)
         ids = cav%mesh%triangles(:, f)
         diameter = max(norm2(cav%mesh%vertices(:, ids(2)) - cav%mesh%vertices(:, ids(1))), &
                        norm2(cav%mesh%vertices(:, ids(3)) - cav%mesh%vertices(:, ids(2))), &
                        norm2(cav%mesh%vertices(:, ids(1)) - cav%mesh%vertices(:, ids(3))))
         if (cav%mesh%area(f) <= sqrt(epsilon(1.0_wp))*diameter**2) slivers = slivers + 1
      end do
      call check(error, slivers > 0, more="Fixture no longer contains a sliver face")
      if (allocated(error)) return
      call check(error, all(abs(cav%xi0*sqrt(cav%a)/cav%param%gaussian_scale - 1.0_wp) < 0.15_wp), &
                 more="Sliver faces distort the width correction")
   end subroutine test_sliver_widths

end module test_cavity_mt
