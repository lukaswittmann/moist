module test_cavity_iswig
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use moist_cavity, only: cavity_type_iswig, new_cavity_iswig
   use moist_radii, only: default_cpcm_radii, new_radii_custom_atoms, radius_type

   implicit none
   private

   public :: collect_cavity_iswig

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: STEP_SIZE = 1.0E-4_wp
   real(wp), parameter :: ABS_THR = 5.0E-9_wp
   real(wp), parameter :: REL_THR = 5.0E-8_wp

contains

   !> Collect all exported unit tests
   subroutine collect_cavity_iswig(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("Spherical cavity", test_spherical_cavity), &
         & new_unittest("Molecular cavity", test_molecular_cavity), &
         & new_unittest("AreaSum", test_area_summation), &
         & new_unittest("AreaVariants", test_area_variants), &
         & new_unittest("GradientSwitch", test_gradient_switch), &
         & new_unittest("GradientArea", test_gradient_area), &
         & new_unittest("GradientVolume", test_gradient_volume), &
         & new_unittest("AmatProperties", test_amat_properties), &
         & new_unittest("AmatGradient", test_amat_gradient) &
         & ]

   end subroutine collect_cavity_iswig

   !> Smoke test for spherical cavity
   subroutine test_spherical_cavity(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp) :: radii(1)
      real(wp) :: area_ref, volume_ref, swi_ref
      real(wp) :: xyz(3, 1)

      xyz(:, 1) = 0.0_wp
      call new(mol, [1], xyz)

      radii = 6.9_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate(cav)
      call new_cavity_iswig(cav, nleb=1202, &
         & radius_model=radius_model, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      area_ref = 4.0_wp*pi*radii(1)**2
      call check(error, cav%total_area, area_ref, thr=1.0E-11_wp, &
         & more="Single-atom total area does not match")

      volume_ref = 4.0_wp/3.0_wp*pi*radii(1)**3
      call check(error, cav%total_volume, volume_ref, thr=1.0E-11_wp, &
         & more="Single-atom total volume does not match")

      swi_ref = 1.0_wp
      call check(error, sum(cav%f)/cav%ngrid, swi_ref, thr=1.0E-11_wp, &
         & more="Single-atom switching function does not match")

   end subroutine test_spherical_cavity

   !> Smoke test for molecular cavity
   subroutine test_molecular_cavity(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp), allocatable :: radii(:)
      real(wp) :: volume_ref, area_ref, switch_ref
      integer :: ngrid_ref

      call get_structure(mol, "MB16-43", "01")
      allocate (radii(mol%nat))
      radii = 2.0_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate(cav)
      call new_cavity_iswig(cav, radius_model=radius_model, &
         & error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ngrid_ref = 1544
      call check(error, cav%ngrid, ngrid_ref, &
         & more="Number of grid points does not match")

      switch_ref = 1.240536285050911E3_wp
      call check(error, sum(cav%f), switch_ref, thr_abs=ABS_THR, thr_rel=REL_THR, &
         & more="Switching function does not match")

      area_ref = 5.650168713524450e2_wp
      call check(error, cav%total_area, area_ref, thr_abs=ABS_THR, thr_rel=REL_THR, &
         & more="Cavity total area does not match")

      volume_ref = 454.41275406590046_wp
      call check(error, cav%total_volume, volume_ref, thr_abs=ABS_THR, thr_rel=REL_THR, &
         & more="Cavity total volume does not match")

   end subroutine test_molecular_cavity

   !> Test of surface area summations (from gridpoints and from spheres)
   subroutine test_area_summation(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "MB16-43", "01")
      allocate (radii(mol%nat))
      radii = 3.0_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate(cav)
      call new_cavity_iswig(cav, radius_model=radius_model, &
         & error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call check(error, sum(cav%a), sum(cav%asph), thr=thr2, &
         & more="Cavity area summation does not match")

      call check(error, sum(cav%a), cav%total_area, thr=thr2, &
         & more="Cavity total area does not match")

      call check(error, sum(cav%asph), cav%total_area, thr=thr2, &
         & more="Cavity total area from spheres does not match")

   end subroutine test_area_summation

   !> Test of cavity creation routines
   subroutine test_area_variants(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      integer :: nsph, ngrid
      real(wp), allocatable :: radii(:)
      real(wp), allocatable :: asph_full(:), asph_eff(:)
      integer :: num_leb
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model

      call get_structure(mol, "MB16-43", "01")
      nsph = mol%nat
      allocate (radii(nsph))
      radii = 2.0_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      num_leb = 110

      allocate (asph_full(nsph), asph_eff(nsph))

      allocate(cav)
      call new_cavity_iswig(cav, num_leb, 0.0_wp, 0.0_wp, &
         & radius_model=radius_model, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      asph_full = cav%asph
      deallocate (cav)

      allocate(cav)
      call new_cavity_iswig(cav, num_leb, 0.0_wp, 0.0_wp, &
         & radius_model=radius_model, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      asph_eff = cav%asph

      call check(error, sum(asph_full), sum(asph_eff), thr=thr2, &
         & more="Cavity total areas of regular and efficient routine do not match")
      call check(error, maxval(abs(asph_full - asph_eff)), 0.0_wp, thr=thr2, &
         & more="Cavity atomic areas of regular and efficient routine do not match")

   end subroutine test_area_variants

   subroutine test_gradient_switch(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp), allocatable :: radii(:)
      real(wp), allocatable :: num2d(:, :), ana2d(:, :)
      real(wp) :: fwd, bwd
      integer :: i, j

      call get_structure(mol, "MB16-43", "01")
      allocate (radii(mol%nat))
      radii = 2.0_wp
      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      allocate (num2d(3, mol%nat))
      allocate (ana2d(3, mol%nat), source=0.0_wp)
      do i = 1, mol%nat
         do j = 1, 3
            mol%xyz(j, i) = mol%xyz(j, i) + STEP_SIZE
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, radius_model=radius_model, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            fwd = cav%f(1)
            mol%xyz(j, i) = mol%xyz(j, i) - 2*STEP_SIZE
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, radius_model=radius_model, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            bwd = cav%f(1)
            mol%xyz(j, i) = mol%xyz(j, i) + STEP_SIZE
            num2d(j, i) = (fwd - bwd)/(2*STEP_SIZE)
         end do
      end do

      ! simple structural check to avoid unused warnings
      do i = 1, mol%nat
         do j = 1, 3
            call check(error, ana2d(j, i), num2d(j, i), thr_abs=ABS_THR, thr_rel=REL_THR, &
               more="Analytical and numerical gradients do not match for switching function")
         end do
      end do

   end subroutine test_gradient_switch

   subroutine test_gradient_area(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp), allocatable :: radii(:)
      real(wp), allocatable :: num2d(:, :), ana2d(:, :)
      real(wp) :: h, fwd, bwd, ffwd, bbwd, cut_a, cut_f, s
      integer :: i, j, nleb, nlebs(5)

      nlebs = [14, 26, 50, 110, 194]
      cut_a = 0.0_wp
      cut_f = 1.0E-7_wp

      call get_structure(mol, "MB16-43", "06")
      allocate (radii(mol%nat))
      radii = 2.0_wp
      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      allocate (num2d(3, mol%nat))
      h = 5.0e-3_wp

      allocate (ana2d(3, mol%nat), source=0.0_wp)

      do nleb = 1, size(nlebs)

         num2d = 0.0_wp

         do i = 1, mol%nat
            do j = 1, 3

               mol%xyz(j, i) = mol%xyz(j, i) + 2.0_wp*STEP_SIZE
               if (allocated(cav)) deallocate(cav)
               allocate(cav)
               call new_cavity_iswig(cav, nleb=nlebs(nleb), &
                  & cut_a=cut_a, cut_f=cut_f, &
                  & radius_model=radius_model, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, cavity_error%message)
                  return
               end if
               call cav%update(mol, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, cavity_error%message)
                  return
               end if
               ffwd = cav%total_area

               mol%xyz(j, i) = mol%xyz(j, i) - STEP_SIZE
               if (allocated(cav)) deallocate(cav)
               allocate(cav)
               call new_cavity_iswig(cav, nleb=nlebs(nleb), &
                  & cut_a=cut_a, cut_f=cut_f, &
                  & radius_model=radius_model, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, cavity_error%message)
                  return
               end if
               call cav%update(mol, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, cavity_error%message)
                  return
               end if
               fwd = cav%total_area

               mol%xyz(j, i) = mol%xyz(j, i) - 2.0_wp*STEP_SIZE
               if (allocated(cav)) deallocate(cav)
               allocate(cav)
               call new_cavity_iswig(cav, nleb=nlebs(nleb), &
                  & cut_a=cut_a, cut_f=cut_f, &
                  & radius_model=radius_model, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, cavity_error%message)
                  return
               end if
               call cav%update(mol, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, cavity_error%message)
                  return
               end if
               bwd = cav%total_area

               mol%xyz(j, i) = mol%xyz(j, i) - STEP_SIZE
               if (allocated(cav)) deallocate(cav)
               allocate(cav)
               call new_cavity_iswig(cav, nleb=nlebs(nleb), &
                  & cut_a=cut_a, cut_f=cut_f, &
                  & radius_model=radius_model, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, cavity_error%message)
                  return
               end if
               call cav%update(mol, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, cavity_error%message)
                  return
               end if
               bbwd = cav%total_area
               mol%xyz(j, i) = mol%xyz(j, i) + 2.0_wp*STEP_SIZE

               num2d(j, i) = (-ffwd + 8.0_wp*fwd - 8.0_wp*bwd + bbwd)/(12.0_wp*STEP_SIZE)
            end do
         end do

         ana2d = 0.0_wp
         ! Use type-bound gradient on the cavity (3, nat)
         if (allocated(cav)) deallocate(cav)
         allocate(cav)
         call new_cavity_iswig(cav, nleb=nlebs(nleb), &
            & cut_a=cut_a, cut_f=cut_f, &
            & radius_model=radius_model, error=cavity_error)
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
         if (allocated(cav%area_grad)) ana2d = cav%area_grad

         do i = 1, mol%nat
            do j = 1, 3
               call check(error, ana2d(j, i), num2d(j, i), thr_abs=ABS_THR, thr_rel=REL_THR, &
                  more="Analytical and numerical gradients do not match for area")
            end do
         end do

      end do

   end subroutine test_gradient_area

   subroutine test_gradient_volume(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp), allocatable :: radii(:)
      real(wp), allocatable :: num(:), ana(:)
      real(wp) :: ffwd, fwd, bwd, bbwd, h
      integer :: i, j

      h = 1.0e-3_wp
      call get_structure(mol, "MB16-43", "03")
      allocate (radii(mol%nat))
      radii = 2.0_wp
      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      allocate (num(3*mol%nat))
      do i = 1, mol%nat
         do j = 1, 3
            mol%xyz(j, i) = mol%xyz(j, i) + 2.0_wp*h
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, &
               & radius_model=radius_model, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            ffwd = cav%total_volume

            mol%xyz(j, i) = mol%xyz(j, i) - h
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, &
               & radius_model=radius_model, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            fwd = cav%total_volume

            mol%xyz(j, i) = mol%xyz(j, i) - 2.0_wp*h
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, &
               & radius_model=radius_model, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            bwd = cav%total_volume

            mol%xyz(j, i) = mol%xyz(j, i) - h
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, &
               & radius_model=radius_model, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            bbwd = cav%total_volume

            mol%xyz(j, i) = mol%xyz(j, i) + 2.0_wp*h
            num(3*(i - 1) + j) = (-ffwd + 8.0_wp*fwd &
               & - 8.0_wp*bwd + bbwd) / (12.0_wp*h)
         end do
      end do

      ! Compute analytical volume gradient
      if (allocated(cav)) deallocate(cav)
      allocate(cav)
      call new_cavity_iswig(cav, radius_model=radius_model, &
         & error=cavity_error)
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

      allocate (ana(3*mol%nat))
      do i = 1, mol%nat
         do j = 1, 3
            ana(3*(i - 1) + j) = cav%volume_grad(j, i)
         end do
      end do

      do i = 1, size(num)
         call check(error, ana(i), num(i), &
            & thr_abs=ABS_THR, thr_rel=REL_THR, &
            & more="Volume gradient mismatch")
      end do

   end subroutine test_gradient_volume

   !> Smoke test for A-matrix properties (symmetry, sign, diagonal dominance)
   subroutine test_amat_properties(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp), allocatable :: radii(:)
      real(wp), allocatable :: amat(:, :)
      real(wp) :: offdiag_sum
      integer :: i, j, n

      call get_structure(mol, "MB16-43", "01")
      allocate (radii(mol%nat))
      radii = 2.0_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate(cav)
      call new_cavity_iswig(cav, radius_model=radius_model, &
         & error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      n = cav%ngrid
      allocate(amat(n, n))
      call cav%get_amat(amat, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ! Test symmetry: amat(i,j) == amat(j,i)
      do i = 1, n
         do j = i + 1, n
            call check(error, amat(i, j), amat(j, i), thr=1.0E-11_wp, &
               & more="A-matrix is not symmetric")
            if (allocated(error)) return
         end do
      end do

      ! Test positive diagonal
      do i = 1, n
         if (amat(i, i) <= 0.0_wp) then
            call test_failed(error, "A-matrix diagonal element is not positive")
            return
         end if
      end do

      ! Test negative off-diagonal
      do i = 1, n
         do j = 1, n
            if (i /= j .and. amat(i, j) > 0.0_wp) then
               call test_failed(error, "A-matrix off-diagonal element is not <= 0")
               return
            end if
         end do
      end do

      ! Test diagonal dominance: amat(i,i) >= sum(|amat(i,j)|, j/=i)
      do i = 1, n
         offdiag_sum = 0.0_wp
         do j = 1, n
            if (j /= i) offdiag_sum = offdiag_sum + abs(amat(i, j))
         end do
         if (amat(i, i) < offdiag_sum) then
            call test_failed(error, "A-matrix is not diagonally dominant")
            return
         end if
      end do

   end subroutine test_amat_properties

   !> Numerical vs analytical test for contracted A-matrix derivative.
   !> Uses 5-point finite differences on q1^T A q2 to verify
   !> contract_amat1_q1q2_rA.
   subroutine test_amat_gradient(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp), allocatable :: radii(:)
      real(wp), allocatable :: amat(:, :)
      real(wp), allocatable :: q1(:), q2(:)
      real(wp), allocatable :: ana_grad(:, :), num_grad(:, :)
      real(wp) :: ffwd, fwd, bwd, bbwd, energy
      integer :: i, j, n, iat, iax

      call get_structure(mol, "MB16-43", "13")
      allocate(radii(mol%nat))
      radii = 4.0_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ! Build reference cavity to set up charge vectors
      allocate(cav)
      call new_cavity_iswig(cav, cut_f=0.01_wp, radius_model=radius_model, &
         & error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      n = cav%ngrid

      ! Set up test charge vectors
      allocate(q1(n), q2(n))
      do i = 1, n
         q1(i) = real(i, wp) / real(n + 1, wp)
         if (mod(i, 2) == 0) then
            q2(i) = real(i, wp) / real(n + 1, wp)
         else
            q2(i) = -real(i, wp) / real(n + 1, wp)
         end if
      end do

      ! Compute analytical gradient
      allocate(ana_grad(3, mol%nat))
      call cav%contract_amat1_q1q2_rA(q1, q2, ana_grad, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ! Compute numerical gradient via 5-point stencil
      allocate(num_grad(3, mol%nat))
      num_grad = 0.0_wp

      do iat = 1, mol%nat
         do iax = 1, 3

            ! +2h
            mol%xyz(iax, iat) = mol%xyz(iax, iat) + 2.0_wp * STEP_SIZE
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, cut_f=0.01_wp, radius_model=radius_model, &
               & error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            allocate(amat(cav%ngrid, cav%ngrid))
            call cav%get_amat(amat, cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            ffwd = dot_product(q1, matmul(amat, q2))
            deallocate(amat)

            ! +h
            mol%xyz(iax, iat) = mol%xyz(iax, iat) - STEP_SIZE
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, cut_f=0.01_wp, radius_model=radius_model, &
               & error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            allocate(amat(cav%ngrid, cav%ngrid))
            call cav%get_amat(amat, cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            fwd = dot_product(q1, matmul(amat, q2))
            deallocate(amat)

            ! -h
            mol%xyz(iax, iat) = mol%xyz(iax, iat) - 2.0_wp * STEP_SIZE
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, cut_f=0.01_wp, radius_model=radius_model, &
               & error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            allocate(amat(cav%ngrid, cav%ngrid))
            call cav%get_amat(amat, cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            bwd = dot_product(q1, matmul(amat, q2))
            deallocate(amat)

            ! -2h
            mol%xyz(iax, iat) = mol%xyz(iax, iat) - STEP_SIZE
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, cut_f=0.01_wp, radius_model=radius_model, &
               & error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            allocate(amat(cav%ngrid, cav%ngrid))
            call cav%get_amat(amat, cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            bbwd = dot_product(q1, matmul(amat, q2))
            deallocate(amat)

            ! Restore position
            mol%xyz(iax, iat) = mol%xyz(iax, iat) + 2.0_wp * STEP_SIZE

            ! 5-point stencil
            num_grad(iax, iat) = (-ffwd + 8.0_wp*fwd &
               & - 8.0_wp*bwd + bbwd) / (12.0_wp * STEP_SIZE)
         end do
      end do

      ! Compare analytical vs numerical
      do iat = 1, mol%nat
         do iax = 1, 3
            call check(error, ana_grad(iax, iat), num_grad(iax, iat), &
               & thr_abs=ABS_THR, thr_rel=REL_THR, &
               & more="A-matrix gradient mismatch")
            if (allocated(error)) return
         end do
      end do

   end subroutine test_amat_gradient


end module test_cavity_iswig
