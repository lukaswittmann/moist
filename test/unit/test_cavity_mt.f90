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
   use moist_cavity_mc, only: cavity_type_mc, new_cavity_mc
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_isodensity_internal, only: moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_drop_lsf_isodensity_param, only: moist_cavity_drop_lsf_isodensity_param_type
   use moist_cavity_drop_switching, only: moist_cavity_drop_swif_sigmoid_bump_type, new_swif_sigmoid_bump
   use moist_cavity_type, only: cavity_type
   use moist_channels_coupling, only: coupling_type, density_matrix_request_type, nuclear_charge_request_type
   use moist_channels_fields, only: field_query_type
   use pcmgen_gauss_conductor, only: conductor_state, conductor_assemble, conductor_medium, conductor_solve
   use test_helpers, only: check_moist_error, read_printout
   implicit none(type, external)
   private
   public :: collect_cavity_mt

   !> SvdW level set with a constant spatial gradient
   type, extends(moist_cavity_drop_lsf_svdw_type) :: constant_gradient_lsf_type
      !> Value of every gradient component
      real(wp) :: gradient = 0.0_wp
   contains
      procedure :: f012_r => constant_gradient_f012_r
   end type constant_gradient_lsf_type

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
                  new_unittest("sliver_widths", test_sliver_widths), &
                  new_unittest("limited_widths_warning", test_limited_widths), &
                  new_unittest("thread_determinism", test_threads), &
                  new_unittest("deferred_build", test_deferred_build), &
                  new_unittest("build_keeps_surface", test_build_keeps_surface), &
                  new_unittest("bad_update_without_context", test_bad_update_without_context, should_fail=.true.), &
                  new_unittest("bad_update_without_lsf", test_bad_update_without_lsf, should_fail=.true.), &
                  new_unittest("bad_update_no_atoms", test_bad_update_no_atoms, should_fail=.true.), &
                  new_unittest("bad_update_area_tolerance", test_bad_update_area_tolerance, should_fail=.true.), &
                  new_unittest("bad_gaussian_scale", test_bad_gaussian_scale, should_fail=.true.), &
                  new_unittest("bad_spacing", test_bad_spacing, should_fail=.true.), &
                  new_unittest("bad_zero_area_tolerance", test_bad_zero_area_tolerance, should_fail=.true.), &
                  new_unittest("bad_negative_switch_start", test_bad_switch_start, should_fail=.true.), &
                  new_unittest("bad_critical_prefactor", test_bad_critical_prefactor, should_fail=.true.), &
                  new_unittest("bad_negative_weight_cutoff", test_bad_negative_weight_cutoff, should_fail=.true.), &
                  new_unittest("bad_critical_bump", test_bad_critical_bump, should_fail=.true.), &
                  new_unittest("bad_weight_cutoff", test_bad_weight_cutoff, should_fail=.true.), &
                  new_unittest("bad_nan_gradient", test_bad_nan_gradient, should_fail=.true.), &
                  new_unittest("bad_overflow_gradient", test_bad_overflow_gradient, should_fail=.true.)]
   end subroutine collect_cavity_mt

   !> Build a custom-radius SvdW fixture with an explicit thread count
   !>
   !> @param[out] error Test failure
   !> @param[in,out] ctx Borrowed run context
   !> @param[out] cav Constructed cavity
   !> @param[in] xyz Atomic coordinates, bohr
   !> @param[in] radii Atomic radii, bohr
   !> @param[in] param Optional settings
   !> @param[in] unit Optional unit receiving copies of warnings
   !> @param[in] nthreads Optional thread count, default one
   subroutine fixture(error, ctx, cav, xyz, radii, param, unit, nthreads)
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
      !> Optional unit receiving copies of warnings
      integer, intent(in), optional :: unit
      !> Optional thread count
      integer, intent(in), optional :: nthreads
      type(structure_type) :: mol
      class(radius_type), allocatable :: radius_model
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(mctc_error), allocatable :: local_error
      integer :: i, nt

      nt = 1
      if (present(nthreads)) nt = nthreads
      call new_context(ctx, verbosity=0, nthreads=nt, unit=unit)
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
      call cav%configure(curvature=.false., error=local_error)
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

   !> Overly narrow initial widths are clamped at every point and reported
   subroutine test_limited_widths(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(moist_cavity_mt_parameters_type) :: param
      character(len=256) :: lines(8)
      character(len=16) :: total
      integer :: unit, nline

      param%spacing = 0.3_wp
      param%target_area = 0.6_wp
      ! The uniform layer needs about 4.9/2 times the initial width
      param%gaussian_scale = 2.0_wp
      open (newunit=unit, status="scratch", action="readwrite")
      call fixture(error, ctx, cav, reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [3.0_wp], param, unit=unit)
      if (.not. allocated(error)) call read_printout(unit, lines, nline)
      close (unit)
      if (allocated(error)) return
      call check(error, all(abs(cav%xi0*sqrt(cav%a)/param%gaussian_scale - 1.25_wp) < 1.0e-12_wp), &
         & more="Width factors not clamped at the upper limit")
      if (allocated(error)) return
      write (total, "(i0)") cav%ngrid
      call check(error, nline, 1, more="Expected exactly one warning line")
      if (allocated(error)) return
      call check(error, trim(lines(1)) == "[Warning] Marching tetrahedra width correction reached its limit at " &
         & //trim(total)//" of "//trim(total)//" points; crowded points remain less accurate", more=trim(lines(1)))
   end subroutine test_limited_widths

   !> Two threads reproduce the serial quadrature exactly
   subroutine test_threads(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: serial_ctx, threaded_ctx
      type(cavity_type_mt) :: serial, threaded
      type(moist_cavity_mt_parameters_type) :: param
      real(wp) :: xyz(3, 2)

      xyz = 0.0_wp
      xyz(1, :) = [-1.6_wp, 1.6_wp]
      param%spacing = 0.3_wp
      param%target_area = 0.6_wp
      call fixture(error, serial_ctx, serial, xyz, [2.0_wp, 2.5_wp], param)
      if (allocated(error)) return
      call fixture(error, threaded_ctx, threaded, xyz, [2.0_wp, 2.5_wp], param, nthreads=2)
      if (allocated(error)) return
      call check(error, serial%ngrid, threaded%ngrid)
      if (allocated(error)) return
      call check(error, all(serial%patch_index == threaded%patch_index), more="Thread count changes the points")
      if (allocated(error)) return
      call check(error, maxval(abs(serial%xyz - threaded%xyz)), 0.0_wp, thr=0.0_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(serial%a - threaded%a)), 0.0_wp, thr=0.0_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(serial%xi0 - threaded%xi0)), 0.0_wp, thr=0.0_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(serial%normal0 - threaded%normal0)), 0.0_wp, thr=0.0_wp)
   end subroutine test_threads

   !> MT and MC defer an isodensity sphere to the build phase and integrate it there
   subroutine test_deferred_build(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      class(radius_type), allocatable :: radius_model
      type(cavity_type_mt), target :: mt
      type(cavity_type_mc), target :: mc
      type(moist_cavity_mt_parameters_type) :: param
      type(field_query_type) :: fields
      real(wp) :: radius
      type(mctc_error), allocatable :: err

      call new_context(ctx, verbosity=0, nthreads=1)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      ! One s Gaussian, exp(-r**2/2): spherical density with unit LSF gradient near 1.7 bohr
      call lsf%new(sh_atom=[1], sh_l=[0], sh_nprim=[1], exps=[0.5_wp], coeffs=[1.0_wp], error=err, &
                   param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=0.01_wp, scale=30.0_wp))
      if (.not. allocated(err)) call new_radii_custom_atoms([4.0_wp], radius_model, err)
      param%spacing = 0.15_wp
      param%target_area = 0.5_wp
      if (.not. allocated(err)) call new_cavity_mt(mt, radius_model, lsf, err, param, ctx)
      if (.not. allocated(err)) call new_cavity_mc(mc, radius_model, lsf, err, ctx=ctx)
      if (.not. allocated(err)) call mt%update(mol, err)
      if (.not. allocated(err)) call mc%update(mol, err)
      call check_moist_error(error, err, "deferred update")
      if (allocated(error)) return
      call mt%list_fields(fields)
      call check(error, fields%index_of("mesh_area") == 0, more="Unbuilt cavity lists mesh fields")
      if (allocated(error)) return
      call check(error, mt%is_updated() .and. .not. mt%is_built() .and. mt%ngrid == 0, &
         & more="MT built a density surface before the build phase")
      if (allocated(error)) return
      call check(error, mc%is_updated() .and. .not. mc%is_built() .and. .not. allocated(mc%total_area), &
         & more="MC integrated a density surface before the build phase")
      if (allocated(error)) return
      call build_with_density(mt, err)
      if (.not. allocated(err)) call build_with_density(mc, err)
      call check_moist_error(error, err, "build")
      if (allocated(error)) return
      call check(error, mt%is_built() .and. mt%ngrid > 0 .and. mc%is_built(), more="Build left no surface")
      if (allocated(error)) return
      call zero_radius(error, mt%lsf_model, radius)
      if (allocated(error)) return
      ! The critical switch is one on this unit-gradient surface
      call check(error, mt%total_area, mt%mesh_area, rel=.true., thr=1.0e-12_wp)
      if (allocated(error)) return
      ! Linear edge interpolation of the concave density LSF places vertices outside:
      ! MT 1.5 % area and 0.8 % volume high, MC 0.4 % and 0.4 %
      call check(error, mt%mesh_area/(4.0_wp*pi*radius**2), 1.0_wp, thr=2.0e-2_wp, more="MT density sphere area")
      if (allocated(error)) return
      call check(error, mt%mesh_volume/(4.0_wp*pi*radius**3/3.0_wp), 1.0_wp, thr=1.0e-2_wp, &
         & more="MT density sphere volume")
      if (allocated(error)) return
      call check(error, mc%total_area/(4.0_wp*pi*radius**2), 1.0_wp, thr=1.0e-2_wp, more="MC density sphere area")
      if (allocated(error)) return
      call check(error, mc%total_volume/(4.0_wp*pi*radius**3/3.0_wp), 1.0_wp, thr=1.0e-2_wp, &
         & more="MC density sphere volume")
   end subroutine test_deferred_build

   !> Bisect the zero of a centered spherical level set along the x axis
   !>
   !> @param[out] error Test failure
   !> @param[in,out] lsf Level set, negative at the origin
   !> @param[out] radius Zero-set radius, bohr
   subroutine zero_radius(error, lsf, radius)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Level set
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Zero-set radius
      real(wp), intent(out) :: radius
      type(mctc_error), allocatable :: err
      real(wp) :: lower, upper, value
      integer :: iter

      lower = 0.0_wp
      upper = 4.0_wp
      call lsf%set_max_deriv(0)
      do iter = 1, 60
         radius = 0.5_wp*(lower + upper)
         call lsf%prepare([radius, 0.0_wp, 0.0_wp], err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         call lsf%f0(value)
         if (value < 0.0_wp) then
            lower = radius
         else
            upper = radius
         end if
      end do
   end subroutine zero_radius

   !> Run the build phase with a unit s density and Z = 1
   !>
   !> @param[in,out] cav Updated deferring cavity
   !> @param[out] err Staging, answer or build failure
   subroutine build_with_density(cav, err)
      !> Updated deferring cavity
      class(cavity_type), intent(inout), target :: cav
      !> Staging, answer or build failure
      type(mctc_error), allocatable, intent(out) :: err
      type(coupling_type), pointer :: coupling
      real(wp) :: dcart(1, 1)

      dcart = 1.0_wp
      call cav%new_coupling(coupling, err)
      if (allocated(err)) return
      call cav%prepare_build(coupling, err)
      if (.not. allocated(err)) then
         do while (coupling%next())
            select type (request => coupling%request())
            type is (density_matrix_request_type)
               call coupling%answer("dcart", dcart, err)
            type is (nuclear_charge_request_type)
               call coupling%answer("z", [1.0_wp], err)
            end select
            if (allocated(err)) exit
         end do
      end if
      if (.not. allocated(err)) call cav%build(coupling, err)
      call cav%release_coupling(coupling)
   end subroutine build_with_density

   !> The build phase of a geometric level set keeps the surface of the update
   subroutine test_build_keeps_surface(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt), target :: cav
      type(coupling_type), pointer :: coupling
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: area(:)

      call fixture(error, ctx, cav, reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [3.0_wp], coarse_parameters())
      if (allocated(error)) return
      area = cav%a
      call cav%new_coupling(coupling, err)
      if (.not. allocated(err)) then
         call cav%prepare_build(coupling, err)
         if (.not. allocated(err)) call cav%build(coupling, err)
         call cav%release_coupling(coupling)
      end if
      call check_moist_error(error, err, "build")
      if (allocated(error)) return
      call check(error, cav%is_built() .and. cav%ngrid == size(area), more="Build changed the point count")
      if (allocated(error)) return
      call check(error, maxval(abs(cav%a - area)), 0.0_wp, thr=0.0_wp)
   end subroutine test_build_keeps_surface

   !> Reject an update of a cavity constructed without a run context
   subroutine test_bad_update_without_context(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      class(radius_type), allocatable :: radius_model
      type(mctc_error), allocatable :: err

      call lsf%new()
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_radii_custom_atoms([3.0_wp], radius_model, err)
      if (.not. allocated(err)) call new_cavity_mt(cav, radius_model, lsf, err)
      if (.not. allocated(err)) call cav%update(mol, err)
      call check_moist_error(error, err)
   end subroutine test_bad_update_without_context

   !> Reject an update after the level set was removed from a constructed cavity
   subroutine test_bad_update_without_lsf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(mctc_error), allocatable :: err

      call lsf%new()
      call new_sphere_cavity(ctx, cav, mol, lsf, err)
      if (.not. allocated(err)) then
         deallocate (cav%lsf_model)
         call cav%update(mol, err)
      end if
      call check_moist_error(error, err)
   end subroutine test_bad_update_without_lsf

   !> Reject a structure without atoms
   subroutine test_bad_update_no_atoms(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol, empty
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(mctc_error), allocatable :: err
      integer :: none(0)
      real(wp) :: nowhere(3, 0)

      call lsf%new()
      call new_sphere_cavity(ctx, cav, mol, lsf, err)
      if (.not. allocated(err)) then
         call new(empty, none, nowhere)
         call cav%update(empty, err)
      end if
      call check_moist_error(error, err)
   end subroutine test_bad_update_no_atoms

   !> Settings changed after construction are validated again on update
   subroutine test_bad_update_area_tolerance(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(mctc_error), allocatable :: err

      call lsf%new()
      call new_sphere_cavity(ctx, cav, mol, lsf, err)
      if (.not. allocated(err)) then
         cav%param%area_tolerance = 1.0_wp
         call cav%update(mol, err)
      end if
      call check_moist_error(error, err)
   end subroutine test_bad_update_area_tolerance

   !> Reject a nonpositive Gaussian width coefficient
   subroutine test_bad_gaussian_scale(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_mt_parameters_type) :: param
      type(mctc_error), allocatable :: err

      call lsf%new()
      param%gaussian_scale = 0.0_wp
      call new_sphere_cavity(ctx, cav, mol, lsf, err, param)
      call check_moist_error(error, err)
   end subroutine test_bad_gaussian_scale

   !> Reject a nonpositive lattice spacing
   subroutine test_bad_spacing(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_mt_parameters_type) :: param
      type(mctc_error), allocatable :: err

      call lsf%new()
      param%spacing = 0.0_wp
      call new_sphere_cavity(ctx, cav, mol, lsf, err, param)
      call check_moist_error(error, err)
   end subroutine test_bad_spacing

   !> Reject a zero area tolerance
   subroutine test_bad_zero_area_tolerance(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_mt_parameters_type) :: param
      type(mctc_error), allocatable :: err

      call lsf%new()
      param%area_tolerance = 0.0_wp
      call new_sphere_cavity(ctx, cav, mol, lsf, err, param)
      call check_moist_error(error, err)
   end subroutine test_bad_zero_area_tolerance

   !> Reject a negative critical-switch start
   subroutine test_bad_switch_start(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_mt_parameters_type) :: param
      type(mctc_error), allocatable :: err

      call lsf%new()
      param%w_0ls_from = -0.1_wp
      call new_sphere_cavity(ctx, cav, mol, lsf, err, param)
      call check_moist_error(error, err)
   end subroutine test_bad_switch_start

   !> Reject a nonpositive critical-bump prefactor
   subroutine test_bad_critical_prefactor(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_mt_parameters_type) :: param
      type(mctc_error), allocatable :: err

      call lsf%new()
      param%w_0ls_a = 0.0_wp
      call new_sphere_cavity(ctx, cav, mol, lsf, err, param)
      call check_moist_error(error, err)
   end subroutine test_bad_critical_prefactor

   !> Reject a negative weight cutoff
   subroutine test_bad_negative_weight_cutoff(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_mt_parameters_type) :: param
      type(mctc_error), allocatable :: err

      call lsf%new()
      param%weight_cutoff = -1.0e-3_wp
      call new_sphere_cavity(ctx, cav, mol, lsf, err, param)
      call check_moist_error(error, err)
   end subroutine test_bad_negative_weight_cutoff

   !> Reject a nonpositive critical-bump exponent
   subroutine test_bad_critical_bump(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_mt_parameters_type) :: param
      type(mctc_error), allocatable :: err

      call lsf%new()
      param%w_0ls_p = 0.0_wp
      call new_sphere_cavity(ctx, cav, mol, lsf, err, param)
      call check_moist_error(error, err)
   end subroutine test_bad_critical_bump

   !> Reject a weight cutoff that suppresses every point
   subroutine test_bad_weight_cutoff(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_mt_parameters_type) :: param
      type(mctc_error), allocatable :: err

      call lsf%new()
      param%weight_cutoff = 1.0_wp
      call new_sphere_cavity(ctx, cav, mol, lsf, err, param)
      call check_moist_error(error, err)
   end subroutine test_bad_weight_cutoff

   !> Reject a level set with a NaN spatial gradient
   subroutine test_bad_nan_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(constant_gradient_lsf_type) :: lsf
      type(mctc_error), allocatable :: err

      call lsf%new()
      lsf%gradient = ieee_value(1.0_wp, ieee_quiet_nan)
      call new_sphere_cavity(ctx, cav, mol, lsf, err, coarse_parameters())
      if (.not. allocated(err)) call cav%update(mol, err)
      call check_moist_error(error, err)
   end subroutine test_bad_nan_gradient

   !> Reject a finite level-set gradient whose norm overflows
   subroutine test_bad_overflow_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(cavity_type_mt) :: cav
      type(structure_type) :: mol
      type(constant_gradient_lsf_type) :: lsf
      type(mctc_error), allocatable :: err

      call lsf%new()
      lsf%gradient = huge(1.0_wp)
      call new_sphere_cavity(ctx, cav, mol, lsf, err, coarse_parameters())
      if (.not. allocated(err)) call cav%update(mol, err)
      call check_moist_error(error, err)
   end subroutine test_bad_overflow_gradient

   !> Coarse settings for cavities that fail after meshing
   function coarse_parameters() result(param)
      !> Coarse mesh and patch settings
      type(moist_cavity_mt_parameters_type) :: param

      param%spacing = 0.5_wp
      param%target_area = 2.0_wp
   end function coarse_parameters

   !> Construct a single-sphere cavity of radius 3 bohr without updating it
   !>
   !> @param[in,out] ctx Borrowed run context
   !> @param[out] cav Constructed cavity
   !> @param[out] mol Single atom at the origin
   !> @param[in] lsf Level set copied into the cavity
   !> @param[out] err Construction failure
   !> @param[in] param Optional settings
   subroutine new_sphere_cavity(ctx, cav, mol, lsf, err, param)
      !> Borrowed run context
      type(moist_context_type), intent(inout), target :: ctx
      !> Constructed cavity
      type(cavity_type_mt), intent(out) :: cav
      !> Single atom at the origin
      type(structure_type), intent(out) :: mol
      !> Level set copied into the cavity
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Construction failure
      type(mctc_error), allocatable, intent(out) :: err
      !> Optional settings
      type(moist_cavity_mt_parameters_type), intent(in), optional :: param
      class(radius_type), allocatable :: radius_model

      call new_context(ctx, verbosity=0, nthreads=1)
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_radii_custom_atoms([3.0_wp], radius_model, err)
      if (allocated(err)) return
      call new_cavity_mt(cav, radius_model, lsf, err, param, ctx)
   end subroutine new_sphere_cavity

   !> Parent value and Hessian with a constant spatial gradient
   !>
   !> @param[in] self Level set
   !> @param[out] lsf0 Level-set value
   !> @param[out] lsf1_r Spatial gradient (3)
   !> @param[out] lsf2_rr Spatial Hessian (3, 3)
   subroutine constant_gradient_f012_r(self, lsf0, lsf1_r, lsf2_rr)
      !> Level set
      class(constant_gradient_lsf_type), intent(in) :: self
      !> Level-set value
      real(wp), intent(out), optional :: lsf0
      !> Spatial gradient
      real(wp), intent(out), optional :: lsf1_r(:)
      !> Spatial Hessian
      real(wp), intent(out), optional :: lsf2_rr(:, :)

      call self%moist_cavity_drop_lsf_svdw_type%f012_r(lsf0, lsf1_r, lsf2_rr)
      if (present(lsf1_r)) lsf1_r = self%gradient
   end subroutine constant_gradient_f012_r

end module test_cavity_mt
