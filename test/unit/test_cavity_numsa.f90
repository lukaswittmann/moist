module test_cavity_numsa
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_io, only: structure_type
   use mctc_io_convert, only: aatoau
   use mstore, only: get_structure
   use moist_data_radii_legacy, only: get_radius_func
   use moist_cavity_numsa, only: cavity_type_numsa, new_cavity_numsa
   use moist_radii, only : static_radius_type
   use moist_radii, only : new_d3_radii, new_bondi_radii, new_cosmo_radii, new_cpcm_radii

   implicit none
   private

   public :: collect_cavity_numsa

   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

   subroutine collect_cavity_numsa(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("mb01_area", test_mb01), &
         & new_unittest("mb02_area", test_mb02), &
         & new_unittest("mb03_area", test_mb03), &
         & new_unittest("mb01_gradient", test_mb01_grad), &
         & new_unittest("mb05_gradient", test_mb05_grad), &
         & new_unittest("mb03_gradient", test_mb03_grad) &
         & ]
   end subroutine collect_cavity_numsa


   subroutine test_mb01(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_numsa) :: cav
      type(mctc_error), allocatable :: cavity_error
      type(static_radius_type) :: radii
      integer :: i

      real(wp), parameter :: probe = 1.4_wp * aatoau
      integer, parameter :: nleb = 110
      real(wp), parameter :: ref(16) = [&
         & 1.98249964603498E+02_wp, &
         & 9.34967918541344E+01_wp, &
         & 7.26746425976157E+01_wp, &
         & 3.72308705072405E+01_wp, &
         & 1.00057039498616E+02_wp, &
         & 8.72703799995796E+01_wp, &
         & 1.75563553107864E+01_wp, &
         & 5.79324044295481E+01_wp, &
         & 9.81701754804677E-03_wp, &
         & 1.05256238904348E+02_wp, &
         & 6.62363240313345E+01_wp, &
         & 1.44944528018566E+02_wp, &
         & 3.33346853562456E+01_wp, &
         & 5.79746582175529E+01_wp, &
         & 6.69252984752073E+00_wp, &
         & 4.86484694486336E+01_wp]

      call get_structure(mol, "MB16-43", "01")

      call new_d3_radii(radii)
      call new_cavity_numsa(cav, nleb=nleb, probe_r=probe, radii=radii, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      do i = 1, mol%nat
         call check(error, cav%asph(i), ref(i), thr=thr2)
      end do

   end subroutine test_mb01


   subroutine test_mb02(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_numsa) :: cav
      type(mctc_error), allocatable :: cavity_error
      type(static_radius_type) :: radii
      integer :: i
      real(wp), parameter :: probe = 1.2_wp * aatoau
      integer, parameter :: nleb = 230
      real(wp), parameter :: ref(16) = [&
         & 2.86084867868854E+01_wp, &
         & 7.50937555534059E+01_wp, &
         & 8.05879869880977E+01_wp, &
         & 8.24020440962820E+01_wp, &
         & 6.48136730299052E+01_wp, &
         & 1.97586791688521E+01_wp, &
         & 4.90632288004349E+01_wp, &
         & 5.29220735596789E+01_wp, &
         & 9.14599031786151E+01_wp, &
         & 1.38294851260743E+01_wp, &
         & 9.02032751808618E+01_wp, &
         & 1.13713659875286E+02_wp, &
         & 9.83820274680035E+01_wp, &
         & 5.95926059359978E+01_wp, &
         & 2.96614646358023E+00_wp, &
         & 1.44874751490690E+02_wp]

      call get_structure(mol, "MB16-43", "02")

      call new_bondi_radii(radii)
      call new_cavity_numsa(cav, nleb=nleb, probe_r=probe, radii=radii, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      do i = 1, mol%nat
         call check(error, cav%asph(i), ref(i), thr=thr2)
      end do

   end subroutine test_mb02


   subroutine test_mb03(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_numsa) :: cav
      type(mctc_error), allocatable :: cavity_error
      type(static_radius_type) :: radii
      integer :: i
      real(wp), parameter :: probe = 0.2_wp * aatoau
      integer, parameter :: nleb = 110
      real(wp), parameter :: ref(16) = [&
         & 4.93447390726497E+01_wp, &
         & 5.42387849176901E+01_wp, &
         & 2.58043997374119E+01_wp, &
         & 3.26892803192176E+01_wp, &
         & 1.27988010759842E+01_wp, &
         & 9.45810634518707E+01_wp, &
         & 3.43532470377123E+01_wp, &
         & 2.76341416140764E+01_wp, &
         & 2.74903764017798E+01_wp, &
         & 2.85813017859723E+01_wp, &
         & 7.99313005786035E+01_wp, &
         & 1.26258175473983E+02_wp, &
         & 5.38016574162998E+01_wp, &
         & 4.16287245622076E+01_wp, &
         & 9.95930646536509E+01_wp, &
         & 2.36024718294637E+01_wp]

      call get_structure(mol, "MB16-43", "03")

      call new_cosmo_radii(radii)
      call new_cavity_numsa(cav, nleb=nleb, probe_r=probe, radii=radii, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      do i = 1, mol%nat
         call check(error, cav%asph(i), ref(i), thr=thr2)
      end do

   end subroutine test_mb03


   subroutine test_mb01_grad(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_numsa) :: cav
      type(mctc_error), allocatable :: cavity_error
      type(static_radius_type) :: radii
      real(wp), allocatable :: grad_analytic(:, :)
      real(wp), allocatable :: grad_numeric(:, :)
      integer :: i, j

      real(wp), parameter :: probe = 1.4_wp * aatoau
      integer, parameter :: nleb = 110

      call get_structure(mol, "MB16-43", "01")

      call new_d3_radii(radii)
      call new_cavity_numsa(cav, nleb=nleb, probe_r=probe, radii=radii, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%get_gradient()

      allocate(grad_analytic(3, mol%nat))
      grad_analytic = cav%area_grad

      allocate(grad_numeric(3, mol%nat))
      call compute_total_area_gradient_fd(cav, mol, step, grad_numeric, error)
      if (allocated(error)) return

      do i = 1, mol%nat
         do j = 1, 3
            call check(error, grad_analytic(j, i), grad_numeric(j, i), thr=thr2)
            if (allocated(error)) return
         end do
      end do

   end subroutine test_mb01_grad


   subroutine test_mb05_grad(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_numsa) :: cav
      type(mctc_error), allocatable :: cavity_error
      type(static_radius_type) :: radii
      real(wp), allocatable :: grad_analytic(:, :)
      real(wp), allocatable :: grad_numeric(:, :)
      integer :: i, j

      real(wp), parameter :: probe = 1.2_wp * aatoau
      integer, parameter :: nleb = 110

      call get_structure(mol, "MB16-43", "05")

      call new_bondi_radii(radii)
      call new_cavity_numsa(cav, nleb=nleb, probe_r=probe, radii=radii, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%get_gradient()

      allocate(grad_analytic(3, mol%nat))
      grad_analytic = cav%area_grad

      allocate(grad_numeric(3, mol%nat))
      call compute_total_area_gradient_fd(cav, mol, step, grad_numeric, error)
      if (allocated(error)) return

      do i = 1, mol%nat
         do j = 1, 3
            call check(error, grad_analytic(j, i), grad_numeric(j, i), thr=thr2)
            if (allocated(error)) return
         end do
      end do

   end subroutine test_mb05_grad


   subroutine test_mb03_grad(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_numsa) :: cav
      type(mctc_error), allocatable :: cavity_error
      type(static_radius_type) :: radii
      real(wp), allocatable :: grad_analytic(:, :)
      real(wp), allocatable :: grad_numeric(:, :)
      integer :: i, j

      real(wp), parameter :: probe = 0.2_wp * aatoau
      integer, parameter :: nleb = 110

      call get_structure(mol, "MB16-43", "03")

      call new_cosmo_radii(radii)
      call new_cavity_numsa(cav, nleb=nleb, probe_r=probe, radii=radii, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%get_gradient()

      allocate(grad_analytic(3, mol%nat))
      grad_analytic = cav%area_grad

      allocate(grad_numeric(3, mol%nat))
      call compute_total_area_gradient_fd(cav, mol, step, grad_numeric, error)
      if (allocated(error)) return

      do i = 1, mol%nat
         do j = 1, 3
            call check(error, grad_analytic(j, i), grad_numeric(j, i), thr=thr2)
            if (allocated(error)) return
         end do
      end do

   end subroutine test_mb03_grad


   !> Compute gradient of total cavity area via 4-point central finite differences
   !>
   !> Uses the 4-point central difference stencil (4th-order accurate):
   !>
   !>   f'(x) = (f(x-2h) - 8 f(x-h) + 8 f(x+h) - f(x+2h)) / (12 h)
   !>
   !> Each atomic coordinate is perturbed in turn, the cavity is re-evaluated,
   !> and the corresponding total-area derivative is assembled.
   subroutine compute_total_area_gradient_fd(cav, mol, stepsize, grad, error)
      type(cavity_type_numsa), intent(inout) :: cav
      type(structure_type), intent(inout) :: mol
      real(wp), intent(in) :: stepsize
      real(wp), intent(out) :: grad(:, :)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: cavity_error

      integer :: iat, ic
      real(wp) :: x0
      real(wp) :: am2, am1, ap1, ap2

      do iat = 1, mol%nat
         do ic = 1, 3
            x0 = mol%xyz(ic, iat)

            mol%xyz(ic, iat) = x0 - 2.0_wp * stepsize
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            am2 = cav%total_area
            mol%xyz(ic, iat) = x0 - stepsize
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            am1 = cav%total_area
            mol%xyz(ic, iat) = x0 + stepsize
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            ap1 = cav%total_area
            mol%xyz(ic, iat) = x0 + 2.0_wp * stepsize
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            ap2 = cav%total_area

            mol%xyz(ic, iat) = x0

            grad(ic, iat) = (am2 - 8.0_wp*am1 + 8.0_wp*ap1 - ap2) / (12.0_wp*stepsize)
         end do
      end do

   end subroutine compute_total_area_gradient_fd

end module test_cavity_numsa
