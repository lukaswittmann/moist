!> Test suite for the constant level-set multiplier `scale` of the DROP LSFs
!>
!> - every value and derivative an accessor returns is `scale` times the one at `scale = 1`
!> - scale-invariant: `normalized_f01_rA` and the surface-free radius of the exclusion
!>   certificate, which are read at the scaled value
module test_cavity_drop_lsf_scale
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: expect_error, get_test_radii, get_test_points
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_cfc_param, only: moist_cavity_drop_lsf_cfc_param_type
   use moist_cavity_drop_lsf_sses, only: moist_cavity_drop_lsf_sses_type
   use moist_cavity_drop_lsf_sses_param, only: moist_cavity_drop_lsf_sses_param_type
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: collect_cavity_drop_lsf_scale

   !> Spatial dimension
   integer, parameter :: ndim = 3
   !> Non-trivial multiplier of the scaled instance
   real(wp), parameter :: alpha = 0.37_wp
   !> Relative tolerance of a scaled output against the array's largest entry
   real(wp), parameter :: rel_tol = 1.0e-14_wp
   !> Evaluation points per structure
   integer, parameter :: npoint = 4

contains

   !> Collect the scale tests
   !>
   !> @param[out] testsuite Unit tests of this suite
   subroutine collect_cavity_drop_lsf_scale(testsuite)
      !> Unit tests of this suite
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("svdw_scale_linear", test_svdw_scale_linear), &
                  new_unittest("cfc_scale_linear", test_cfc_scale_linear), &
                  new_unittest("sses_scale_linear", test_sses_scale_linear), &
                  new_unittest("family_scale_defaults", test_family_scale_defaults), &
                  new_unittest("bad_svdw_scale", test_bad_svdw_scale, should_fail=.true.), &
                  new_unittest("bad_cfc_scale", test_bad_cfc_scale, should_fail=.true.)]
   end subroutine collect_cavity_drop_lsf_scale

   !> SvdW accessors at `scale = alpha` are `alpha` times those at one
   !>
   !> @param[out] error Test failure
   subroutine test_svdw_scale_linear(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_svdw_type) :: one, scaled

      call one%new(moist_cavity_drop_lsf_svdw_param_type(scale=1.0_wp))
      call scaled%new(moist_cavity_drop_lsf_svdw_param_type(scale=alpha))
      call check_family(error, one, scaled)
   end subroutine test_svdw_scale_linear

   !> CFC accessors at `scale = alpha` are `alpha` times those at one
   !>
   !> @param[out] error Test failure
   subroutine test_cfc_scale_linear(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_cfc_type) :: one, scaled

      call one%new(moist_cavity_drop_lsf_cfc_param_type(scale=1.0_wp))
      call scaled%new(moist_cavity_drop_lsf_cfc_param_type(scale=alpha))
      call check_family(error, one, scaled)
   end subroutine test_cfc_scale_linear

   !> SSES accessors at `scale = alpha` are `alpha` times those at one
   !>
   !> @param[out] error Test failure
   subroutine test_sses_scale_linear(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: one, scaled

      call one%new(moist_cavity_drop_lsf_sses_param_type(scale=1.0_wp))
      call scaled%new(moist_cavity_drop_lsf_sses_param_type(scale=alpha))
      call check_family(error, one, scaled)
   end subroutine test_sses_scale_linear

   !> The constructors copy each family's default scale into the base
   !>
   !> @param[out] error Test failure
   subroutine test_family_scale_defaults(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_svdw_type) :: svdw
      type(moist_cavity_drop_lsf_cfc_type) :: cfc
      type(moist_cavity_drop_lsf_sses_type) :: sses

      call svdw%new()
      call cfc%new()
      call sses%new()
      call check(error, svdw%scale, 1.0_wp, "SvdW default scale")
      if (allocated(error)) return
      call check(error, cfc%scale, 0.25_wp, "CFC default scale")
      if (allocated(error)) return
      call check(error, sses%scale, 0.088_wp, "SSES default scale")
   end subroutine test_family_scale_defaults

   !> SvdW refuses a non-positive scale at evaluation
   !>
   !> @param[out] error Test failure, set only on the expected refusal
   subroutine test_bad_svdw_scale(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(mctc_error), allocatable :: err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call lsf%new(moist_cavity_drop_lsf_svdw_param_type(scale=0.0_wp))
      call lsf%update(mol, radii)
      call lsf%prepare(mol%xyz(:, 1) + 1.0_wp, err)
      call expect_error(error, err, "SvdW scale must be finite and positive")
   end subroutine test_bad_svdw_scale

   !> CFC refuses a non-positive scale at evaluation
   !>
   !> @param[out] error Test failure, set only on the expected refusal
   subroutine test_bad_cfc_scale(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(mctc_error), allocatable :: err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call lsf%new(moist_cavity_drop_lsf_cfc_param_type(scale=-1.0_wp))
      call lsf%update(mol, radii)
      call lsf%prepare(mol%xyz(:, 1) + 1.0_wp, err)
      call expect_error(error, err, "CFC scale must be finite and positive")
   end subroutine test_bad_cfc_scale

   !> Compare every accessor family of two instances that differ only in `scale`
   !>
   !> @param[out]   error  Test failure
   !> @param[inout] one    Instance at `scale = 1`
   !> @param[inout] scaled Instance at `scale = alpha`
   subroutine check_family(error, one, scaled)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Instance at `scale = 1`
      class(moist_cavity_drop_lsf_type), intent(inout) :: one
      !> Instance at `scale = alpha`
      class(moist_cavity_drop_lsf_type), intent(inout) :: scaled

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      character(len=*), parameter :: records(2) = ["01", "02"]
      integer :: irec, ip

      do irec = 1, size(records)
         call get_structure(mol, "MB16-43", records(irec))
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, npoint)
         call one%update(mol, radii)
         call scaled%update(mol, radii)
         call check(error, scaled%scale, alpha, "base scale copied from the parameter")
         if (allocated(error)) return
         call one%set_max_deriv(4)
         call scaled%set_max_deriv(4)
         do ip = 1, npoint
            call check_point(error, one, scaled, mol%nat, points(:, ip))
            if (allocated(error)) then
               error%message = "MB16-43/"//records(irec)//" point "//to_string(ip)//": "// &
                               error%message
               return
            end if
         end do
         call check_batch(error, one, scaled, mol%nat, points)
         if (allocated(error)) then
            error%message = "MB16-43/"//records(irec)//": "//error%message
            return
         end if
      end do
   end subroutine check_family

   !> Compare the per-point accessors at one point
   !>
   !> @param[out]   error  Test failure
   !> @param[inout] one    Instance at `scale = 1`
   !> @param[inout] scaled Instance at `scale = alpha`
   !> @param[in]    nat    Atom count
   !> @param[in]    point  Evaluation point
   subroutine check_point(error, one, scaled, nat, point)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Instance at `scale = 1`
      class(moist_cavity_drop_lsf_type), intent(inout) :: one
      !> Instance at `scale = alpha`
      class(moist_cavity_drop_lsf_type), intent(inout) :: scaled
      !> Atom count
      integer, intent(in) :: nat
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)

      type(mctc_error), allocatable :: err
      real(wp) :: v(ndim, nat), vrad(nat), w1(ndim), w2(6), s0(2), r0(2)
      real(wp) :: g1(ndim, 2), g2(6, 2), g3(10, 2), g4(15, 2)
      real(wp) :: t1(ndim, nat, 2), t2(ndim, ndim, nat, 2), t3(6, ndim, nat, 2), &
                  t4(10, ndim, nat, 2)
      real(wp) :: q1(nat, 2), q2(ndim, nat, 2), q3(6, nat, 2)
      real(wp) :: p2(ndim, nat, ndim, nat, 2), p3(ndim, ndim, nat, ndim, nat, 2), &
                  p4(6, ndim, nat, ndim, nat, 2)
      real(wp) :: rr2(nat, nat, 2), rr3(ndim, nat, nat, 2), rr4(6, nat, nat, 2)
      real(wp) :: ar2(ndim, nat, nat, 2), ar3(ndim, ndim, nat, nat, 2), ar4(6, ndim, nat, nat, 2)
      real(wp) :: tan0(2), tan1(ndim, 2), tan2(6, 2), tan3(10, 2)
      real(wp) :: hvp1(ndim, nat, 2), hvp2(ndim, ndim, nat, 2), hvp3(6, ndim, nat, 2)
      real(wp) :: hr1(nat, 2), hr2(ndim, nat, 2), hr3(6, nat, 2)
      real(wp) :: vj(ndim, nat, 2), vjr(nat, 2), nval(2), nder(ndim, nat, 2)
      integer :: k, a

      do a = 1, nat
         do k = 1, ndim
            v(k, a) = sin(real(3*a + k, wp))
         end do
         vrad(a) = cos(real(a, wp))
      end do
      w1 = [0.3_wp, -0.7_wp, 0.2_wp]
      w2 = [0.5_wp, -0.1_wp, 0.4_wp, 0.2_wp, -0.3_wp, 0.6_wp]

      g1 = 0.0_wp; g2 = 0.0_wp; g3 = 0.0_wp; g4 = 0.0_wp
      t1 = 0.0_wp; t2 = 0.0_wp; t3 = 0.0_wp; t4 = 0.0_wp
      q1 = 0.0_wp; q2 = 0.0_wp; q3 = 0.0_wp
      p2 = 0.0_wp; p3 = 0.0_wp; p4 = 0.0_wp
      rr2 = 0.0_wp; rr3 = 0.0_wp; rr4 = 0.0_wp
      ar2 = 0.0_wp; ar3 = 0.0_wp; ar4 = 0.0_wp
      hvp1 = 0.0_wp; hvp2 = 0.0_wp; hvp3 = 0.0_wp
      hr1 = 0.0_wp; hr2 = 0.0_wp; hr3 = 0.0_wp
      vj = 0.0_wp; vjr = 0.0_wp; nder = 0.0_wp

      call one%prepare(point, err)
      if (.not. allocated(err)) call scaled%prepare(point, err)
      if (allocated(err)) then
         call test_failed(error, "prepare: "//err%message)
         return
      end if
      if (one%unreached() .or. scaled%unreached()) then
         call test_failed(error, "test point unreached")
         return
      end if

      ! Spatial jet
      call one%f3_rrr(lsf0=s0(1), lsf1_r=g1(:, 1), lsf2_rr=g2(:, 1), lsf3_rrr=g3(:, 1))
      call scaled%f3_rrr(lsf0=s0(2), lsf1_r=g1(:, 2), lsf2_rr=g2(:, 2), lsf3_rrr=g3(:, 2))
      call one%f4_rrrr(g4(:, 1))
      call scaled%f4_rrrr(g4(:, 2))
      call compare(error, "f0", s0(1:1), s0(2:2), alpha)
      if (allocated(error)) return
      call compare(error, "f1_r", g1(:, 1), g1(:, 2), alpha)
      if (allocated(error)) return
      call compare(error, "f2_rr", g2(:, 1), g2(:, 2), alpha)
      if (allocated(error)) return
      call compare(error, "f3_rrr", g3(:, 1), g3(:, 2), alpha)
      if (allocated(error)) return
      call compare(error, "f4_rrrr", g4(:, 1), g4(:, 2), alpha)
      if (allocated(error)) return

      ! Nuclear and radius families
      call one%f3_rr_rA(lsf1_rA=t1(:, :, 1), lsf2_r_rA=t2(:, :, :, 1), lsf3_rr_rA=t3(:, :, :, 1))
      call scaled%f3_rr_rA(lsf1_rA=t1(:, :, 2), lsf2_r_rA=t2(:, :, :, 2), lsf3_rr_rA=t3(:, :, :, 2))
      call one%f4_rrr_rA(t4(:, :, :, 1))
      call scaled%f4_rrr_rA(t4(:, :, :, 2))
      call one%f3_rr_rad(lsf1_rad=q1(:, 1), lsf2_r_rad=q2(:, :, 1), lsf3_rr_rad=q3(:, :, 1))
      call scaled%f3_rr_rad(lsf1_rad=q1(:, 2), lsf2_r_rad=q2(:, :, 2), lsf3_rr_rad=q3(:, :, 2))
      call compare(error, "f1_rA", pack(t1(:, :, 1), .true.), pack(t1(:, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f2_r_rA", pack(t2(:, :, :, 1), .true.), pack(t2(:, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f3_rr_rA", pack(t3(:, :, :, 1), .true.), pack(t3(:, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f4_rrr_rA", pack(t4(:, :, :, 1), .true.), pack(t4(:, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f1_rad", q1(:, 1), q1(:, 2), alpha)
      if (allocated(error)) return
      call compare(error, "f2_r_rad", pack(q2(:, :, 1), .true.), pack(q2(:, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f3_rr_rad", pack(q3(:, :, 1), .true.), pack(q3(:, :, 2), .true.), alpha)
      if (allocated(error)) return

      ! Two-atom blocks
      call one%f2_rArB(p2(:, :, :, :, 1))
      call scaled%f2_rArB(p2(:, :, :, :, 2))
      call one%f3_r_rArB(p3(:, :, :, :, :, 1))
      call scaled%f3_r_rArB(p3(:, :, :, :, :, 2))
      call one%f4_rr_rArB(p4(:, :, :, :, :, 1))
      call scaled%f4_rr_rArB(p4(:, :, :, :, :, 2))
      call one%f2_radrad(rr2(:, :, 1))
      call scaled%f2_radrad(rr2(:, :, 2))
      call one%f3_r_radrad(rr3(:, :, :, 1))
      call scaled%f3_r_radrad(rr3(:, :, :, 2))
      call one%f4_rr_radrad(rr4(:, :, :, 1))
      call scaled%f4_rr_radrad(rr4(:, :, :, 2))
      call one%f2_rA_rad(ar2(:, :, :, 1))
      call scaled%f2_rA_rad(ar2(:, :, :, 2))
      call one%f3_r_rA_rad(ar3(:, :, :, :, 1))
      call scaled%f3_r_rA_rad(ar3(:, :, :, :, 2))
      call one%f4_rr_rA_rad(ar4(:, :, :, :, 1))
      call scaled%f4_rr_rA_rad(ar4(:, :, :, :, 2))
      call compare(error, "f2_rArB", pack(p2(:, :, :, :, 1), .true.), pack(p2(:, :, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f3_r_rArB", pack(p3(:, :, :, :, :, 1), .true.), pack(p3(:, :, :, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f4_rr_rArB", pack(p4(:, :, :, :, :, 1), .true.), pack(p4(:, :, :, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f2_radrad", pack(rr2(:, :, 1), .true.), pack(rr2(:, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f3_r_radrad", pack(rr3(:, :, :, 1), .true.), pack(rr3(:, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f4_rr_radrad", pack(rr4(:, :, :, 1), .true.), pack(rr4(:, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f2_rA_rad", pack(ar2(:, :, :, 1), .true.), pack(ar2(:, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f3_r_rA_rad", pack(ar3(:, :, :, :, 1), .true.), pack(ar3(:, :, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "f4_rr_rA_rad", pack(ar4(:, :, :, :, 1), .true.), pack(ar4(:, :, :, :, 2), .true.), alpha)
      if (allocated(error)) return

      ! Directional and adjoint products
      call one%tangent_f0(v, tan0(1))
      call scaled%tangent_f0(v, tan0(2))
      call one%tangent_f1_r(v, tan1(:, 1))
      call scaled%tangent_f1_r(v, tan1(:, 2))
      call one%tangent_f2_rr(v, tan2(:, 1))
      call scaled%tangent_f2_rr(v, tan2(:, 2))
      call one%tangent_f3_rrr(v, tan3(:, 1))
      call scaled%tangent_f3_rrr(v, tan3(:, 2))
      call one%hvp_f1_rA(v, hvp1(:, :, 1), vrad)
      call scaled%hvp_f1_rA(v, hvp1(:, :, 2), vrad)
      call one%hvp_f2_r_rA(v, hvp2(:, :, :, 1), vrad)
      call scaled%hvp_f2_r_rA(v, hvp2(:, :, :, 2), vrad)
      call one%hvp_f3_rr_rA(v, hvp3(:, :, :, 1), vrad)
      call scaled%hvp_f3_rr_rA(v, hvp3(:, :, :, 2), vrad)
      call one%hvp_f1_rad(v, vrad, hr1(:, 1))
      call scaled%hvp_f1_rad(v, vrad, hr1(:, 2))
      call one%hvp_f2_r_rad(v, vrad, hr2(:, :, 1))
      call scaled%hvp_f2_r_rad(v, vrad, hr2(:, :, 2))
      call one%hvp_f3_rr_rad(v, vrad, hr3(:, :, 1))
      call scaled%hvp_f3_rr_rad(v, vrad, hr3(:, :, 2))
      call one%vjp_f1_rA(0.8_wp, w1, w2, vj(:, :, 1))
      call scaled%vjp_f1_rA(0.8_wp, w1, w2, vj(:, :, 2))
      call one%vjp_f1_rad(0.8_wp, w1, w2, vjr(:, 1))
      call scaled%vjp_f1_rad(0.8_wp, w1, w2, vjr(:, 2))
      call compare(error, "tangent_f0", tan0(1:1), tan0(2:2), alpha)
      if (allocated(error)) return
      call compare(error, "tangent_f1_r", tan1(:, 1), tan1(:, 2), alpha)
      if (allocated(error)) return
      call compare(error, "tangent_f2_rr", tan2(:, 1), tan2(:, 2), alpha)
      if (allocated(error)) return
      call compare(error, "tangent_f3_rrr", tan3(:, 1), tan3(:, 2), alpha)
      if (allocated(error)) return
      call compare(error, "hvp_f1_rA", pack(hvp1(:, :, 1), .true.), pack(hvp1(:, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "hvp_f2_r_rA", pack(hvp2(:, :, :, 1), .true.), pack(hvp2(:, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "hvp_f3_rr_rA", pack(hvp3(:, :, :, 1), .true.), pack(hvp3(:, :, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "hvp_f1_rad", hr1(:, 1), hr1(:, 2), alpha)
      if (allocated(error)) return
      call compare(error, "hvp_f2_r_rad", pack(hr2(:, :, 1), .true.), pack(hr2(:, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "hvp_f3_rr_rad", pack(hr3(:, :, 1), .true.), pack(hr3(:, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "vjp_f1_rA", pack(vj(:, :, 1), .true.), pack(vj(:, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "vjp_f1_rad", vjr(:, 1), vjr(:, 2), alpha)
      if (allocated(error)) return

      ! Scale-invariant quantities
      call one%normalized_f01_rA(nval(1), deriv_rA=nder(:, :, 1))
      call scaled%normalized_f01_rA(nval(2), deriv_rA=nder(:, :, 2))
      call compare(error, "normalized_f0", nval(1:1), nval(2:2), 1.0_wp)
      if (allocated(error)) return
      call compare(error, "normalized_f1_rA", pack(nder(:, :, 1), .true.), &
                   pack(nder(:, :, 2), .true.), 1.0_wp)
      if (allocated(error)) return
      r0(1) = one%exclusion_radius(s0(1))
      r0(2) = scaled%exclusion_radius(s0(2))
      call compare(error, "exclusion_radius", r0(1:1), r0(2:2), 1.0_wp, 1.0e-12_wp)
   end subroutine check_point

   !> Compare the batched jet and the batched nuclear vjp
   !>
   !> @param[out]   error  Test failure
   !> @param[inout] one    Instance at `scale = 1`
   !> @param[inout] scaled Instance at `scale = alpha`
   !> @param[in]    nat    Atom count
   !> @param[in]    points Evaluation points (3, np)
   subroutine check_batch(error, one, scaled, nat, points)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Instance at `scale = 1`
      class(moist_cavity_drop_lsf_type), intent(inout) :: one
      !> Instance at `scale = alpha`
      class(moist_cavity_drop_lsf_type), intent(inout) :: scaled
      !> Atom count
      integer, intent(in) :: nat
      !> Evaluation points
      real(wp), intent(in) :: points(:, :)

      type(mctc_error), allocatable :: err
      real(wp), allocatable :: l0(:, :), l1(:, :, :), l2(:, :, :), l3(:, :, :), l4(:, :, :)
      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :), res(:, :, :, :)
      integer :: np, k

      np = size(points, 2)
      allocate (l0(np, 2), l1(ndim, np, 2), l2(6, np, 2), l3(10, np, 2), l4(15, np, 2))
      call one%evaluate_batch(points, 4, l0(:, 1), l1(:, :, 1), l2(:, :, 1), l3(:, :, 1), &
                              l4(:, :, 1), error=err)
      if (allocated(err)) then
         call test_failed(error, "evaluate_batch: "//err%message)
         return
      end if
      call scaled%evaluate_batch(points, 4, l0(:, 2), l1(:, :, 2), l2(:, :, 2), l3(:, :, 2), &
                                    l4(:, :, 2), error=err)
      if (allocated(err)) then
         call test_failed(error, "evaluate_batch: "//err%message)
         return
      end if
      call compare(error, "batch f0", l0(:, 1), l0(:, 2), alpha)
      if (allocated(error)) return
      call compare(error, "batch f1_r", pack(l1(:, :, 1), .true.), pack(l1(:, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "batch f2_rr", pack(l2(:, :, 1), .true.), pack(l2(:, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "batch f3_rrr", pack(l3(:, :, 1), .true.), pack(l3(:, :, 2), .true.), alpha)
      if (allocated(error)) return
      call compare(error, "batch f4_rrrr", pack(l4(:, :, 1), .true.), pack(l4(:, :, 2), .true.), alpha)
      if (allocated(error)) return

      allocate (w0(np), w1(ndim, np), w2(6, np), res(ndim, nat, np, 2))
      do k = 1, np
         w0(k) = 0.5_wp + 0.1_wp*k
         w1(:, k) = [0.2_wp, -0.4_wp, 0.1_wp*k]
         w2(:, k) = [0.3_wp, 0.1_wp, -0.2_wp, 0.05_wp*k, 0.4_wp, -0.1_wp]
      end do
      call one%vjp_f1_rA_batch(points, w0, w1, w2, res(:, :, :, 1), err)
      if (.not. allocated(err)) call scaled%vjp_f1_rA_batch(points, w0, w1, w2, res(:, :, :, 2), err)
      if (allocated(err)) then
         call test_failed(error, "vjp_f1_rA_batch: "//err%message)
         return
      end if
      call compare(error, "batch vjp_f1_rA", pack(res(:, :, :, 1), .true.), &
                   pack(res(:, :, :, 2), .true.), alpha)
   end subroutine check_batch

   !> Check `b = factor a` against the largest entry of `factor a`
   !>
   !> @param[out] error  Test failure naming the output
   !> @param[in]  name   Output name
   !> @param[in]  a      Output at `scale = 1`
   !> @param[in]  b      Output at the other scale
   !> @param[in]  factor Expected ratio `b/a`
   !> @param[in]  tol    Relative tolerance (default `rel_tol`)
   subroutine compare(error, name, a, b, factor, tol)
      !> Test failure naming the output
      type(error_type), allocatable, intent(out) :: error
      !> Output name
      character(len=*), intent(in) :: name
      !> Output at `scale = 1`
      real(wp), intent(in) :: a(:)
      !> Output at the other scale
      real(wp), intent(in) :: b(:)
      !> Expected ratio `b/a`
      real(wp), intent(in) :: factor
      !> Relative tolerance
      real(wp), intent(in), optional :: tol

      real(wp) :: bound, dev

      bound = rel_tol
      if (present(tol)) bound = tol
      bound = bound*max(maxval(abs(factor*a)), tiny(1.0_wp))
      dev = maxval(abs(b - factor*a))
      if (.not. (dev <= bound)) then
         call test_failed(error, name//": deviation "//to_string(dev)//" above "//to_string(bound))
      end if
   end subroutine compare

end module test_cavity_drop_lsf_scale
