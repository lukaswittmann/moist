module test_cavity_iswig
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use moist_cavity, only: cavity_type_iswig, new_cavity_iswig
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, &
      & pcm_amat_surface_weights, pcm_amat_nuclear_gradient
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_radii, only: default_cpcm_radii, new_radii_custom_atoms, radius_type
   use moist_context, only: moist_context_type, new_context
   implicit none (type, external)
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
         & new_unittest("spherical_cavity", test_spherical_cavity), &
         & new_unittest("molecular_cavity", test_molecular_cavity), &
         & new_unittest("area_sum", test_area_summation), &
         & new_unittest("area_variants", test_area_variants), &
         & new_unittest("gradient_switch", test_gradient_switch), &
         & new_unittest("gradient_area", test_gradient_area), &
         & new_unittest("gradient_volume", test_gradient_volume), &
         & new_unittest("amat_properties", test_amat_properties), &
         & new_unittest("amat_gradient", test_amat_gradient), &
         & new_unittest("amat_orca_reference", test_amat_orca_reference), &
         & new_unittest("surface_gradient", test_surface_gradient) &
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
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      xyz(:, 1) = 0.0_wp
      call new(mol, [1], xyz)

      radii = 6.9_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate (cav)
      call new_cavity_iswig(cav, ctx, nleb=1202, &
         & radius_model=radius_model, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      !> The cavity must *borrow* the caller-owned run context, not copy it
      call check(error, associated(cav%ctx, ctx), &
         & more="iSwiG cavity does not borrow the caller-owned run context")
      if (allocated(error)) return

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
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      allocate (radii(mol%nat))
      radii = 2.0_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate (cav)
      call new_cavity_iswig(cav, ctx, radius_model=radius_model, &
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

   !> Test of reface area summations (from gridpoints and from spheres)
   subroutine test_area_summation(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp), allocatable :: radii(:)
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      allocate (radii(mol%nat))
      radii = 3.0_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate (cav)
      call new_cavity_iswig(cav, ctx, radius_model=radius_model, &
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
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

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

      allocate (cav)
      call new_cavity_iswig(cav, ctx, num_leb, 0.0_wp, 0.0_wp, &
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

      allocate (cav)
      call new_cavity_iswig(cav, ctx, num_leb, 0.0_wp, 0.0_wp, &
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
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

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
            if (allocated(cav)) deallocate (cav)
            allocate (cav)
            call new_cavity_iswig(cav, ctx, radius_model=radius_model, error=cavity_error)
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
            if (allocated(cav)) deallocate (cav)
            allocate (cav)
            call new_cavity_iswig(cav, ctx, radius_model=radius_model, error=cavity_error)
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
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

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
               if (allocated(cav)) deallocate (cav)
               allocate (cav)
               call new_cavity_iswig(cav, ctx, nleb=nlebs(nleb), &
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
               if (allocated(cav)) deallocate (cav)
               allocate (cav)
               call new_cavity_iswig(cav, ctx, nleb=nlebs(nleb), &
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
               if (allocated(cav)) deallocate (cav)
               allocate (cav)
               call new_cavity_iswig(cav, ctx, nleb=nlebs(nleb), &
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
               if (allocated(cav)) deallocate (cav)
               allocate (cav)
               call new_cavity_iswig(cav, ctx, nleb=nlebs(nleb), &
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
         if (allocated(cav)) deallocate (cav)
         allocate (cav)
         call new_cavity_iswig(cav, ctx, nleb=nlebs(nleb), &
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
         call cav%get_gradient(cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if
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
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

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
            if (allocated(cav)) deallocate (cav)
            allocate (cav)
            call new_cavity_iswig(cav, ctx, &
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
            if (allocated(cav)) deallocate (cav)
            allocate (cav)
            call new_cavity_iswig(cav, ctx, &
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
            if (allocated(cav)) deallocate (cav)
            allocate (cav)
            call new_cavity_iswig(cav, ctx, &
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
            if (allocated(cav)) deallocate (cav)
            allocate (cav)
            call new_cavity_iswig(cav, ctx, &
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
               & - 8.0_wp*bwd + bbwd)/(12.0_wp*h)
         end do
      end do

      ! Compute analytical volume gradient
      if (allocated(cav)) deallocate (cav)
      allocate (cav)
      call new_cavity_iswig(cav, ctx, radius_model=radius_model, &
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
      call cav%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

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

   !> Test the genuine iSwiG Amat
   subroutine test_amat_properties(error)
      use moist_model_component_pcm_solvers, only: solve_pcm_cholesky
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp), allocatable :: radii(:)
      real(wp), allocatable :: amat(:, :)
      real(wp), allocatable :: rhs(:), q(:)
      real(wp) :: zeta_ij, r_dist, diag_ref
      integer :: i, j, n
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      allocate (radii(mol%nat))
      radii = 2.0_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate (cav)
      call new_cavity_iswig(cav, ctx, radius_model=radius_model, &
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
      allocate (amat(n, n))
      call assemble_pcm_amat(cav%xi0, cav%f, cav%xyz, amat, cavity_error)
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

      ! Test diagonal: A_ii = xi_i * sqrt(2/pi) / F_i
      do i = 1, n
         diag_ref = cav%xi0(i)*sqrt(2.0_wp/pi)/cav%f(i)
         call check(error, amat(i, i), diag_ref, thr=thr, &
            & more="A-matrix diagonal does not match xi*sqrt(2/pi)/F")
         if (allocated(error)) return
      end do

      ! Test off-diagonal: 0 < erf(zeta_ij r)/r <= 2*zeta_ij/sqrt(pi)
      do i = 1, n
         do j = i + 1, n
            if (amat(i, j) <= 0.0_wp) then
               call test_failed(error, "A-matrix off-diagonal element is not > 0")
               return
            end if
            zeta_ij = cav%xi0(i)*cav%xi0(j)/sqrt(cav%xi0(i)**2 + cav%xi0(j)**2)
            if (amat(i, j) > 2.0_wp*zeta_ij/sqrt(pi) + thr) then
               call test_failed(error, &
                  & "A-matrix off-diagonal exceeds Gaussian short-range bound")
               return
            end if
            r_dist = sqrt(sum((cav%xyz(:, i) - cav%xyz(:, j))**2))
            if (amat(i, j) > 1.0_wp/r_dist + thr) then
               call test_failed(error, &
                  & "A-matrix off-diagonal exceeds bare Coulomb bound")
               return
            end if
         end do
      end do

      ! Test positive definiteness: Cholesky factorization must succeed
      allocate (rhs(n), q(n))
      rhs = 1.0_wp
      call solve_pcm_cholesky(amat, rhs, q, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, &
            & "Cholesky failed: "//cavity_error%message)
         return
      end if

   end subroutine test_amat_properties

   !> Cross-validation of the iSwiG CPCM electrostatics against ORCA 6.1.1
   subroutine test_amat_orca_reference(error)
      use moist_model_component_pcm_solvers, only: solve_pcm_cholesky
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_iswig) :: cav
      type(mctc_error), allocatable :: cavity_error
      integer, parameter :: npts = 272
      !> Dielectric constant of the ORCA reference run (water)
      real(wp), parameter :: eps_ref = 80.1510_wp
      !> ORCA CPCM dielectric energy (Eh)
      real(wp), parameter :: e_ref = -0.037719806_wp
      real(wp) :: feps, e_ours
      real(wp) :: ref(7, npts)
      real(wp), allocatable :: phi(:), q_ref(:), q(:), rhs(:), amat(:, :)

      ! ORCA reface table: x, y, z, zeta, F, phi, q_conductor per point
      ref(:, 1) = [3.292331345_wp, 3.586705171_wp, -2.417637305_wp, &
         & 1.21923656069597_wp, 0.09173595197406_wp, 0.048856705_wp, -0.002369459_wp]
      ref(:, 2) = [-7.002896633_wp, 3.586705171_wp, -2.417637305_wp, &
         & 1.21923656069597_wp, 0.99940960491722_wp, 0.017900905_wp, -0.015293723_wp]
      ref(:, 3) = [-1.855282644_wp, 8.734319160_wp, -2.417637305_wp, &
         & 1.21923656069597_wp, 0.99984325044015_wp, 0.011488300_wp, -0.016164838_wp]
      ref(:, 4) = [-1.855282644_wp, -1.560908818_wp, -2.417637305_wp, &
         & 1.21923656069597_wp, 0.05703557303492_wp, 0.040338558_wp, -0.000886460_wp]
      ref(:, 5) = [-1.855282644_wp, 3.586705171_wp, -7.565251294_wp, &
         & 1.21923656069597_wp, 0.99872854994434_wp, 0.055046052_wp, -0.032765883_wp]
      ref(:, 6) = [-1.855282644_wp, 7.226617930_wp, 1.222275453_wp, &
         & 1.36314791518462_wp, 0.06867469897137_wp, -0.131999395_wp, 0.004136603_wp]
      ref(:, 7) = [-1.855282644_wp, 7.226617930_wp, -6.057550063_wp, &
         & 1.36314791518462_wp, 0.99999999999999_wp, 0.052494718_wp, -0.028708067_wp]
      ref(:, 8) = [-1.855282644_wp, -0.053207587_wp, -6.057550063_wp, &
         & 1.36314791518462_wp, 0.71710606819880_wp, 0.040181169_wp, -0.009819519_wp]
      ref(:, 9) = [1.784630114_wp, 3.586705171_wp, 1.222275453_wp, &
         & 1.36314791518462_wp, 0.85765981530706_wp, 0.014930619_wp, -0.010302509_wp]
      ref(:, 10) = [-5.495195402_wp, 3.586705171_wp, 1.222275453_wp, &
         & 1.36314791518462_wp, 0.00021422587559_wp, -0.055676625_wp, 0.000004090_wp]
      ref(:, 11) = [1.784630114_wp, 3.586705171_wp, -6.057550063_wp, &
         & 1.36314791518462_wp, 0.14271878209459_wp, 0.038090841_wp, -0.001469478_wp]
      ref(:, 12) = [-5.495195402_wp, 3.586705171_wp, -6.057550063_wp, &
         & 1.36314791518462_wp, 1.00000000000000_wp, 0.053783875_wp, -0.028434877_wp]
      ref(:, 13) = [1.784630114_wp, 7.226617930_wp, -2.417637305_wp, &
         & 1.36314791518462_wp, 0.99999996267404_wp, 0.032706069_wp, -0.024470370_wp]
      ref(:, 14) = [-5.495195402_wp, 7.226617930_wp, -2.417637305_wp, &
         & 1.36314791518462_wp, 0.99925915341925_wp, -0.000682273_wp, -0.003254562_wp]
      ref(:, 15) = [-5.495195402_wp, -0.053207587_wp, -2.417637305_wp, &
         & 1.36314791518462_wp, 0.99786496413726_wp, 0.032608403_wp, -0.011253502_wp]
      ref(:, 16) = [1.116693678_wp, 6.558681493_wp, 0.554339017_wp, &
         & 1.48400748272934_wp, 0.98869208792631_wp, -0.017787546_wp, 0.003494812_wp]
      ref(:, 17) = [-4.827258966_wp, 6.558681493_wp, 0.554339017_wp, &
         & 1.48400748272934_wp, 0.04083163153654_wp, -0.136776854_wp, 0.002442609_wp]
      ref(:, 18) = [1.116693678_wp, 0.614728849_wp, 0.554339017_wp, &
         & 1.48400748272934_wp, 0.00000010539424_wp, 0.412353647_wp, -0.000000035_wp]
      ref(:, 19) = [-4.827258966_wp, 0.614728849_wp, 0.554339017_wp, &
         & 1.48400748272934_wp, 0.00123849315439_wp, 0.014451239_wp, -0.000008522_wp]
      ref(:, 20) = [1.116693678_wp, 6.558681493_wp, -5.389613627_wp, &
         & 1.48400748272934_wp, 0.99999867601251_wp, 0.050123202_wp, -0.024546366_wp]
      ref(:, 21) = [-4.827258966_wp, 6.558681493_wp, -5.389613627_wp, &
         & 1.48400748272934_wp, 1.00000000000000_wp, 0.045504073_wp, -0.022901534_wp]
      ref(:, 22) = [-4.827258966_wp, 0.614728849_wp, -5.389613627_wp, &
         & 1.48400748272934_wp, 0.99999999991195_wp, 0.050787783_wp, -0.022560133_wp]
      ref(:, 23) = [6.896218754_wp, 0.023388444_wp, -4.954577518_wp, &
         & 2.51606090252714_wp, 0.80916564394798_wp, 0.008842478_wp, 0.000514243_wp]
      ref(:, 24) = [4.401780257_wp, 2.517826941_wp, -4.954577518_wp, &
         & 2.51606090252714_wp, 0.00001314543483_wp, 0.027257608_wp, -0.000000066_wp]
      ref(:, 25) = [4.401780257_wp, -2.471050053_wp, -4.954577518_wp, &
         & 2.51606090252714_wp, 0.00051948530839_wp, 0.008690474_wp, -0.000000184_wp]
      ref(:, 26) = [4.401780257_wp, 0.023388444_wp, -7.449016015_wp, &
         & 2.51606090252714_wp, 0.14397471642388_wp, 0.002131712_wp, 0.000507628_wp]
      ref(:, 27) = [4.401780257_wp, 1.787222821_wp, -6.718411894_wp, &
         & 2.81304160679008_wp, 0.00828695747367_wp, 0.005395145_wp, 0.000021188_wp]
      ref(:, 28) = [4.401780257_wp, -1.740445932_wp, -6.718411894_wp, &
         & 2.81304160679008_wp, 0.06506483903895_wp, 0.003120453_wp, 0.000069775_wp]
      ref(:, 29) = [6.165614633_wp, 0.023388444_wp, -3.190743142_wp, &
         & 2.81304160679008_wp, 0.01720651122405_wp, 0.013516290_wp, 0.000038860_wp]
      ref(:, 30) = [6.165614633_wp, 0.023388444_wp, -6.718411894_wp, &
         & 2.81304160679008_wp, 0.90388797684435_wp, 0.006616614_wp, 0.000738448_wp]
      ref(:, 31) = [6.165614633_wp, 1.787222821_wp, -4.954577518_wp, &
         & 2.81304160679008_wp, 0.23484527439512_wp, 0.008848682_wp, 0.000459895_wp]
      ref(:, 32) = [6.165614633_wp, -1.740445932_wp, -4.954577518_wp, &
         & 2.81304160679008_wp, 0.52981343038882_wp, 0.006131547_wp, 0.000707639_wp]
      ref(:, 33) = [5.841944995_wp, 1.463553182_wp, -3.514412781_wp, &
         & 3.06245180526872_wp, 0.00109574948087_wp, 0.013307272_wp, 0.000002372_wp]
      ref(:, 34) = [5.841944995_wp, -1.416776293_wp, -3.514412781_wp, &
         & 3.06245180526872_wp, 0.01213211266558_wp, 0.012853710_wp, 0.000013306_wp]
      ref(:, 35) = [5.841944995_wp, 1.463553182_wp, -6.394742256_wp, &
         & 3.06245180526872_wp, 0.56092475519694_wp, 0.006531664_wp, 0.000838058_wp]
      ref(:, 36) = [5.841944995_wp, -1.416776293_wp, -6.394742256_wp, &
         & 3.06245180526872_wp, 0.79714814725963_wp, 0.004888750_wp, 0.000530867_wp]
      ref(:, 37) = [0.459800119_wp, 4.762520678_wp, 1.270433022_wp, &
         & 1.82083354788149_wp, 0.02066797649315_wp, -0.039606026_wp, 0.000325987_wp]
      ref(:, 38) = [-6.433920818_wp, 4.762520678_wp, 1.270433022_wp, &
         & 1.82083354788149_wp, 0.34400885820192_wp, -0.083994147_wp, 0.009792229_wp]
      ref(:, 39) = [-2.987060349_wp, 8.209381146_wp, 1.270433022_wp, &
         & 1.82083354788149_wp, 0.98760721519801_wp, -0.101161090_wp, 0.042345201_wp]
      ref(:, 40) = [-2.987060349_wp, 4.762520678_wp, 4.717293490_wp, &
         & 1.82083354788149_wp, 0.00090939945548_wp, -0.088173248_wp, 0.000010541_wp]
      ref(:, 41) = [-2.987060349_wp, 7.199819089_wp, 3.707731433_wp, &
         & 2.03575379438756_wp, 0.92993125957553_wp, -0.118770500_wp, 0.032200460_wp]
      ref(:, 42) = [-2.987060349_wp, 7.199819089_wp, -1.166865389_wp, &
         & 2.03575379438756_wp, 0.00041883537985_wp, -0.032372173_wp, 0.000003581_wp]
      ref(:, 43) = [-0.549761938_wp, 4.762520678_wp, 3.707731433_wp, &
         & 2.03575379438756_wp, 0.01587074852488_wp, -0.062608351_wp, 0.000153759_wp]
      ref(:, 44) = [-5.424358760_wp, 4.762520678_wp, 3.707731433_wp, &
         & 2.03575379438756_wp, 0.00466568420055_wp, -0.107376706_wp, 0.000089183_wp]
      ref(:, 45) = [-5.424358760_wp, 4.762520678_wp, -1.166865389_wp, &
         & 2.03575379438756_wp, 0.00028087859576_wp, -0.025155703_wp, 0.000002667_wp]
      ref(:, 46) = [-0.549761938_wp, 7.199819089_wp, 1.270433022_wp, &
         & 2.03575379438756_wp, 0.69464093314709_wp, -0.081510062_wp, 0.019654355_wp]
      ref(:, 47) = [-5.424358760_wp, 7.199819089_wp, 1.270433022_wp, &
         & 2.03575379438756_wp, 0.99730830385762_wp, -0.106580933_wp, 0.037797957_wp]
      ref(:, 48) = [-5.424358760_wp, 2.325222267_wp, 1.270433022_wp, &
         & 2.03575379438756_wp, 0.00000014805526_wp, 0.035799218_wp, -0.000000005_wp]
      ref(:, 49) = [-0.997014530_wp, 6.752566497_wp, 3.260478841_wp, &
         & 2.21624801697078_wp, 0.90530524183282_wp, -0.101003907_wp, 0.026006092_wp]
      ref(:, 50) = [-4.977106169_wp, 6.752566497_wp, 3.260478841_wp, &
         & 2.21624801697078_wp, 0.88728801664143_wp, -0.121499698_wp, 0.026878033_wp]
      ref(:, 51) = [-0.997014530_wp, 6.752566497_wp, -0.719612797_wp, &
         & 2.21624801697078_wp, 0.00000259086140_wp, -0.015726676_wp, 0.000000000_wp]
      ref(:, 52) = [-4.977106169_wp, 6.752566497_wp, -0.719612797_wp, &
         & 2.21624801697078_wp, 0.11182919899493_wp, -0.066151606_wp, 0.002232047_wp]
      ref(:, 53) = [0.799808865_wp, 1.411034563_wp, -7.540991738_wp, &
         & 2.51606090252714_wp, 0.21325862342442_wp, 0.022677579_wp, 0.000092053_wp]
      ref(:, 54) = [0.799808865_wp, 3.174868939_wp, -6.810387617_wp, &
         & 2.81304160679008_wp, 0.29598010459624_wp, 0.041186129_wp, -0.001446753_wp]
      ref(:, 55) = [0.799808865_wp, -0.352799813_wp, -6.810387617_wp, &
         & 2.81304160679008_wp, 0.00004642077621_wp, 0.026085587_wp, -0.000000141_wp]
      ref(:, 56) = [2.563643241_wp, 1.411034563_wp, -6.810387617_wp, &
         & 2.81304160679008_wp, 0.00000038719954_wp, 0.014949535_wp, 0.000000000_wp]
      ref(:, 57) = [-0.964025512_wp, 1.411034563_wp, -6.810387617_wp, &
         & 2.81304160679008_wp, 0.21287871479541_wp, 0.046942629_wp, -0.001208230_wp]
      ref(:, 58) = [2.563643241_wp, 3.174868939_wp, -5.046553241_wp, &
         & 2.81304160679008_wp, 0.00000025106991_wp, 0.048399410_wp, -0.000000003_wp]
      ref(:, 59) = [-0.964025512_wp, -0.352799813_wp, -5.046553241_wp, &
         & 2.81304160679008_wp, 0.00005951541747_wp, 0.042408570_wp, -0.000000360_wp]
      ref(:, 60) = [2.239973602_wp, 2.851199301_wp, -6.486717979_wp, &
         & 3.06245180526872_wp, 0.00111876462997_wp, 0.024402347_wp, -0.000000694_wp]
      ref(:, 61) = [-0.640355873_wp, 2.851199301_wp, -6.486717979_wp, &
         & 3.06245180526872_wp, 0.00013688507745_wp, 0.068578556_wp, -0.000001555_wp]
      ref(:, 62) = [-0.640355873_wp, -0.029130175_wp, -6.486717979_wp, &
         & 3.06245180526872_wp, 0.06433666934552_wp, 0.030779089_wp, -0.000039943_wp]
      ref(:, 63) = [-0.872997815_wp, 1.842757685_wp, 4.550380871_wp, &
         & 1.88276666175501_wp, 0.00074119098760_wp, 0.005120134_wp, -0.000009651_wp]
      ref(:, 64) = [-7.539951615_wp, 1.842757685_wp, 4.550380871_wp, &
         & 1.88276666175501_wp, 0.97414874596005_wp, -0.058132258_wp, 0.018424234_wp]
      ref(:, 65) = [-4.206474715_wp, 5.176234585_wp, 4.550380871_wp, &
         & 1.88276666175501_wp, 0.01812435451429_wp, -0.109974136_wp, 0.000361813_wp]
      ref(:, 66) = [-4.206474715_wp, -1.490719216_wp, 4.550380871_wp, &
         & 1.88276666175501_wp, 0.04150178914482_wp, 0.036284419_wp, -0.000551300_wp]
      ref(:, 67) = [-4.206474715_wp, 1.842757685_wp, 7.883857771_wp, &
         & 1.88276666175501_wp, 0.99988001328340_wp, -0.059606369_wp, 0.017807824_wp]
      ref(:, 68) = [-4.206474715_wp, 4.199881806_wp, 6.907504992_wp, &
         & 2.10499712072727_wp, 0.97691422625890_wp, -0.078195181_wp, 0.015606125_wp]
      ref(:, 69) = [-4.206474715_wp, -0.514366437_wp, 6.907504992_wp, &
         & 2.10499712072727_wp, 0.99987533126932_wp, -0.029636522_wp, 0.012157395_wp]
      ref(:, 70) = [-4.206474715_wp, -0.514366437_wp, 2.193256750_wp, &
         & 2.10499712072727_wp, 0.00000025876704_wp, 0.039312708_wp, -0.000000004_wp]
      ref(:, 71) = [-1.849350594_wp, 1.842757685_wp, 6.907504992_wp, &
         & 2.10499712072727_wp, 0.90267543029733_wp, -0.045946165_wp, 0.011409450_wp]
      ref(:, 72) = [-6.563598836_wp, 1.842757685_wp, 6.907504992_wp, &
         & 2.10499712072727_wp, 0.99999498895838_wp, -0.062713177_wp, 0.013834278_wp]
      ref(:, 73) = [-6.563598836_wp, 1.842757685_wp, 2.193256750_wp, &
         & 2.10499712072727_wp, 0.01317392677995_wp, -0.021698749_wp, -0.000038531_wp]
      ref(:, 74) = [-1.849350594_wp, 4.199881806_wp, 4.550380871_wp, &
         & 2.10499712072727_wp, 0.00002164841567_wp, -0.043718167_wp, -0.000000107_wp]
      ref(:, 75) = [-6.563598836_wp, 4.199881806_wp, 4.550380871_wp, &
         & 2.10499712072727_wp, 0.67122548632883_wp, -0.083081967_wp, 0.010009210_wp]
      ref(:, 76) = [-1.849350594_wp, -0.514366437_wp, 4.550380871_wp, &
         & 2.10499712072727_wp, 0.00003827228073_wp, 0.078908173_wp, -0.000001556_wp]
      ref(:, 77) = [-6.563598836_wp, -0.514366437_wp, 4.550380871_wp, &
         & 2.10499712072727_wp, 0.99295736652091_wp, -0.031813082_wp, 0.014097628_wp]
      ref(:, 78) = [-2.281890929_wp, 3.767341470_wp, 6.474964657_wp, &
         & 2.29163060258203_wp, 0.55915801006906_wp, -0.063302393_wp, 0.004613806_wp]
      ref(:, 79) = [-6.131058500_wp, 3.767341470_wp, 6.474964657_wp, &
         & 2.29163060258203_wp, 0.99925450214966_wp, -0.077624619_wp, 0.013840312_wp]
      ref(:, 80) = [-2.281890929_wp, -0.081826101_wp, 6.474964657_wp, &
         & 2.29163060258203_wp, 0.95855957846650_wp, -0.014100288_wp, 0.005556128_wp]
      ref(:, 81) = [-6.131058500_wp, -0.081826101_wp, 6.474964657_wp, &
         & 2.29163060258203_wp, 0.99999887473438_wp, -0.044680373_wp, 0.011330476_wp]
      ref(:, 82) = [-6.131058500_wp, 3.767341470_wp, 2.625797085_wp, &
         & 2.29163060258203_wp, 0.00039183373975_wp, -0.069690088_wp, 0.000003585_wp]
      ref(:, 83) = [-6.131058500_wp, -0.081826101_wp, 2.625797085_wp, &
         & 2.29163060258203_wp, 0.18218862681761_wp, -0.009037967_wp, 0.000444250_wp]
      ref(:, 84) = [-6.037999733_wp, -3.188356667_wp, 1.462400225_wp, &
         & 2.51606090252714_wp, 0.99238031813194_wp, 0.065395710_wp, -0.014951829_wp]
      ref(:, 85) = [-3.543561236_wp, -5.682795164_wp, 1.462400225_wp, &
         & 2.51606090252714_wp, 0.84088436340670_wp, 0.063336446_wp, -0.010274919_wp]
      ref(:, 86) = [-3.543561236_wp, -3.188356667_wp, 3.956838722_wp, &
         & 2.51606090252714_wp, 0.00000031568341_wp, 0.122657678_wp, -0.000000013_wp]
      ref(:, 87) = [-3.543561236_wp, -3.188356667_wp, -1.032038272_wp, &
         & 2.51606090252714_wp, 0.73117960101129_wp, 0.063055056_wp, -0.009875376_wp]
      ref(:, 88) = [-3.543561236_wp, -4.952191044_wp, 3.226234601_wp, &
         & 2.81304160679008_wp, 0.13733680697170_wp, 0.066002160_wp, -0.001269059_wp]
      ref(:, 89) = [-3.543561236_wp, -1.424522291_wp, -0.301434151_wp, &
         & 2.81304160679008_wp, 0.00016482125737_wp, 0.032209471_wp, -0.000000089_wp]
      ref(:, 90) = [-3.543561236_wp, -4.952191044_wp, -0.301434151_wp, &
         & 2.81304160679008_wp, 0.96205319524920_wp, 0.066652924_wp, -0.010391929_wp]
      ref(:, 91) = [-5.307395612_wp, -3.188356667_wp, 3.226234601_wp, &
         & 2.81304160679008_wp, 0.59210928267325_wp, 0.058896968_wp, -0.005935022_wp]
      ref(:, 92) = [-1.779726860_wp, -3.188356667_wp, -0.301434151_wp, &
         & 2.81304160679008_wp, 0.00000103877230_wp, 0.092973057_wp, -0.000000027_wp]
      ref(:, 93) = [-5.307395612_wp, -3.188356667_wp, -0.301434151_wp, &
         & 2.81304160679008_wp, 0.99798577881884_wp, 0.069167776_wp, -0.012306592_wp]
      ref(:, 94) = [-5.307395612_wp, -1.424522291_wp, 1.462400225_wp, &
         & 2.81304160679008_wp, 0.39932460000139_wp, 0.043689934_wp, -0.002936288_wp]
      ref(:, 95) = [-1.779726860_wp, -4.952191044_wp, 1.462400225_wp, &
         & 2.81304160679008_wp, 0.00147389248811_wp, 0.040730740_wp, -0.000003135_wp]
      ref(:, 96) = [-5.307395612_wp, -4.952191044_wp, 1.462400225_wp, &
         & 2.81304160679008_wp, 0.99915452789980_wp, 0.073123429_wp, -0.012667997_wp]
      ref(:, 97) = [-4.983725973_wp, -1.748191930_wp, 2.902564963_wp, &
         & 3.06245180526872_wp, 0.03736921480784_wp, 0.042399433_wp, -0.000243497_wp]
      ref(:, 98) = [-2.103396498_wp, -4.628521405_wp, 2.902564963_wp, &
         & 3.06245180526872_wp, 0.00000149567612_wp, 0.066306270_wp, -0.000000019_wp]
      ref(:, 99) = [-4.983725973_wp, -4.628521405_wp, 2.902564963_wp, &
         & 3.06245180526872_wp, 0.92730479368583_wp, 0.067784082_wp, -0.009327840_wp]
      ref(:, 100) = [-4.983725973_wp, -1.748191930_wp, 0.022235487_wp, &
         & 3.06245180526872_wp, 0.80728031981703_wp, 0.050235498_wp, -0.004833778_wp]
      ref(:, 101) = [-2.103396498_wp, -4.628521405_wp, 0.022235487_wp, &
         & 3.06245180526872_wp, 0.07977377106631_wp, 0.055625039_wp, -0.000271736_wp]
      ref(:, 102) = [-4.983725973_wp, -4.628521405_wp, 0.022235487_wp, &
         & 3.06245180526872_wp, 0.99988733905701_wp, 0.074186817_wp, -0.010972724_wp]
      ref(:, 103) = [5.194760111_wp, 1.068184530_wp, -1.732346512_wp, &
         & 2.51606090252714_wp, 0.00457873283989_wp, 0.035365432_wp, -0.000004956_wp]
      ref(:, 104) = [2.700321614_wp, 3.562623027_wp, -1.732346512_wp, &
         & 2.51606090252714_wp, 0.00226606072985_wp, 0.055858241_wp, -0.000034009_wp]
      ref(:, 105) = [2.700321614_wp, 1.068184530_wp, 0.762091985_wp, &
         & 2.51606090252714_wp, 0.00000265536500_wp, 0.151465592_wp, -0.000000166_wp]
      ref(:, 106) = [2.700321614_wp, 2.832018907_wp, 0.031487864_wp, &
         & 2.81304160679008_wp, 0.52651434544447_wp, 0.052968143_wp, -0.007458658_wp]
      ref(:, 107) = [4.464155991_wp, 1.068184530_wp, 0.031487864_wp, &
         & 2.81304160679008_wp, 0.34586832633569_wp, 0.062348032_wp, -0.003033722_wp]
      ref(:, 108) = [4.464155991_wp, 2.832018907_wp, -1.732346512_wp, &
         & 2.81304160679008_wp, 0.05048764336415_wp, 0.036624918_wp, -0.000162845_wp]
      ref(:, 109) = [4.464155991_wp, -0.695649846_wp, -1.732346512_wp, &
         & 2.81304160679008_wp, 0.00000024284675_wp, 0.069475961_wp, -0.000000004_wp]
      ref(:, 110) = [4.140486352_wp, 2.508349268_wp, -0.292181775_wp, &
         & 3.06245180526872_wp, 0.85372216615598_wp, 0.044027831_wp, -0.004177346_wp]
      ref(:, 111) = [4.140486352_wp, -0.371980207_wp, -0.292181775_wp, &
         & 3.06245180526872_wp, 0.00000020844318_wp, 0.105513234_wp, -0.000000006_wp]
      ref(:, 112) = [7.178001375_wp, -2.070015444_wp, 2.231609387_wp, &
         & 1.82083354788149_wp, 0.05951711176281_wp, 0.086908668_wp, -0.002186276_wp]
      ref(:, 113) = [3.731140907_wp, 1.376845024_wp, 2.231609387_wp, &
         & 1.82083354788149_wp, 0.49501152488093_wp, 0.022762328_wp, -0.003150914_wp]
      ref(:, 114) = [3.731140907_wp, -5.516875912_wp, 2.231609387_wp, &
         & 1.82083354788149_wp, 0.99598267145560_wp, -0.039993324_wp, 0.030107873_wp]
      ref(:, 115) = [3.731140907_wp, -2.070015444_wp, 5.678469855_wp, &
         & 1.82083354788149_wp, 0.99999090384328_wp, -0.055361585_wp, 0.034129587_wp]
      ref(:, 116) = [3.731140907_wp, -2.070015444_wp, -1.215251081_wp, &
         & 1.82083354788149_wp, 0.00009080803088_wp, 0.109694459_wp, -0.000005350_wp]
      ref(:, 117) = [3.731140907_wp, 0.367282967_wp, 4.668907798_wp, &
         & 2.03575379438756_wp, 0.99878176791358_wp, -0.035774051_wp, 0.022096532_wp]
      ref(:, 118) = [3.731140907_wp, -4.507313855_wp, 4.668907798_wp, &
         & 2.03575379438756_wp, 0.99999925465023_wp, -0.056475253_wp, 0.027318490_wp]
      ref(:, 119) = [3.731140907_wp, 0.367282967_wp, -0.205689024_wp, &
         & 2.03575379438756_wp, 0.00002082181200_wp, 0.123409413_wp, -0.000001128_wp]
      ref(:, 120) = [3.731140907_wp, -4.507313855_wp, -0.205689024_wp, &
         & 2.03575379438756_wp, 0.40058007182303_wp, 0.002767608_wp, 0.000559975_wp]
      ref(:, 121) = [6.168439318_wp, -2.070015444_wp, 4.668907798_wp, &
         & 2.03575379438756_wp, 0.94512123061249_wp, -0.020562029_wp, 0.009477859_wp]
      ref(:, 122) = [1.293842496_wp, -2.070015444_wp, 4.668907798_wp, &
         & 2.03575379438756_wp, 0.56057710548863_wp, -0.008949540_wp, 0.004099773_wp]
      ref(:, 123) = [6.168439318_wp, -2.070015444_wp, -0.205689024_wp, &
         & 2.03575379438756_wp, 0.16398765347291_wp, 0.068892198_wp, -0.002988968_wp]
      ref(:, 124) = [6.168439318_wp, 0.367282967_wp, 2.231609387_wp, &
         & 2.03575379438756_wp, 0.15371327562668_wp, 0.061698242_wp, -0.002671959_wp]
      ref(:, 125) = [1.293842496_wp, 0.367282967_wp, 2.231609387_wp, &
         & 2.03575379438756_wp, 0.00000010065117_wp, 0.193692413_wp, -0.000000011_wp]
      ref(:, 126) = [6.168439318_wp, -4.507313855_wp, 2.231609387_wp, &
         & 2.03575379438756_wp, 0.94848644892200_wp, -0.016325878_wp, 0.010465634_wp]
      ref(:, 127) = [1.293842496_wp, -4.507313855_wp, 2.231609387_wp, &
         & 2.03575379438756_wp, 0.04382754164728_wp, 0.026066015_wp, -0.000526681_wp]
      ref(:, 128) = [5.721186726_wp, -0.079969625_wp, 4.221655206_wp, &
         & 2.21624801697078_wp, 0.89991094131715_wp, -0.008171203_wp, 0.005572075_wp]
      ref(:, 129) = [1.741095087_wp, -0.079969625_wp, 4.221655206_wp, &
         & 2.21624801697078_wp, 0.29589214327814_wp, -0.004010245_wp, -0.000018100_wp]
      ref(:, 130) = [5.721186726_wp, -4.060061263_wp, 4.221655206_wp, &
         & 2.21624801697078_wp, 0.99929153022470_wp, -0.039408497_wp, 0.015529631_wp]
      ref(:, 131) = [1.741095087_wp, -4.060061263_wp, 4.221655206_wp, &
         & 2.21624801697078_wp, 0.93095855387275_wp, -0.029768990_wp, 0.014023894_wp]
      ref(:, 132) = [5.721186726_wp, -0.079969625_wp, 0.241563568_wp, &
         & 2.21624801697078_wp, 0.13663779876487_wp, 0.081275156_wp, -0.002222278_wp]
      ref(:, 133) = [5.721186726_wp, -4.060061263_wp, 0.241563568_wp, &
         & 2.21624801697078_wp, 0.90958835848927_wp, 0.003767250_wp, 0.005301832_wp]
      ref(:, 134) = [1.741095087_wp, -4.060061263_wp, 0.241563568_wp, &
         & 2.21624801697078_wp, 0.00000150204216_wp, 0.153755518_wp, -0.000000116_wp]
      ref(:, 135) = [-5.267958810_wp, 0.359514173_wp, 1.053234067_wp, &
         & 1.78559160824507_wp, 0.00579820213809_wp, 0.017607324_wp, -0.000051713_wp]
      ref(:, 136) = [-1.753068201_wp, 0.359514173_wp, 4.568124676_wp, &
         & 1.78559160824507_wp, 0.00027060946415_wp, 0.014837157_wp, -0.000002945_wp]
      ref(:, 137) = [-1.753068201_wp, 0.359514173_wp, -2.461656542_wp, &
         & 1.78559160824507_wp, 0.00000026204437_wp, 0.076834805_wp, -0.000000010_wp]
      ref(:, 138) = [-1.753068201_wp, -2.125888812_wp, -1.432168918_wp, &
         & 1.99635210804458_wp, 0.00263762180929_wp, 0.068546804_wp, -0.000066496_wp]
      ref(:, 139) = [0.732334784_wp, 0.359514173_wp, 3.538637052_wp, &
         & 1.99635210804458_wp, 0.02109623588977_wp, 0.035563222_wp, -0.000409192_wp]
      ref(:, 140) = [-4.238471186_wp, 0.359514173_wp, -1.432168918_wp, &
         & 1.99635210804458_wp, 0.00135790707632_wp, 0.016454335_wp, 0.000005292_wp]
      ref(:, 141) = [0.732334784_wp, 2.844917158_wp, 1.053234067_wp, &
         & 1.99635210804458_wp, 0.00141173219169_wp, 0.019087479_wp, -0.000009395_wp]
      ref(:, 142) = [-4.238471186_wp, -2.125888812_wp, 1.053234067_wp, &
         & 1.99635210804458_wp, 0.00000120856905_wp, 0.175309905_wp, -0.000000102_wp]
      ref(:, 143) = [0.276254838_wp, 2.388837213_wp, 3.082557106_wp, &
         & 2.17335289406167_wp, 0.00144118165262_wp, -0.001197987_wp, -0.000003302_wp]
      ref(:, 144) = [-3.782391240_wp, -1.669808866_wp, -0.976088972_wp, &
         & 2.17335289406167_wp, 0.60605783661059_wp, 0.031175417_wp, 0.002203581_wp]
      ref(:, 145) = [7.911996409_wp, -1.578818309_wp, 1.753940036_wp, &
         & 2.51606090252714_wp, 0.99752430466304_wp, 0.059537675_wp, -0.017728158_wp]
      ref(:, 146) = [5.417557912_wp, 0.915620188_wp, 1.753940036_wp, &
         & 2.51606090252714_wp, 0.52102576436912_wp, 0.043224590_wp, -0.003623807_wp]
      ref(:, 147) = [5.417557912_wp, -4.073256805_wp, 1.753940036_wp, &
         & 2.51606090252714_wp, 0.00260714873077_wp, -0.028994912_wp, 0.000026820_wp]
      ref(:, 148) = [5.417557912_wp, -1.578818309_wp, 4.248378533_wp, &
         & 2.51606090252714_wp, 0.00299380097580_wp, -0.035287424_wp, 0.000029107_wp]
      ref(:, 149) = [5.417557912_wp, -1.578818309_wp, -0.740498461_wp, &
         & 2.51606090252714_wp, 0.39834641118546_wp, 0.059360519_wp, -0.004829838_wp]
      ref(:, 150) = [5.417557912_wp, 0.185016068_wp, 3.517774412_wp, &
         & 2.81304160679008_wp, 0.08120963602580_wp, 0.001772463_wp, 0.000072143_wp]
      ref(:, 151) = [5.417557912_wp, -3.342652685_wp, 3.517774412_wp, &
         & 2.81304160679008_wp, 0.00005384823718_wp, -0.048537393_wp, 0.000000641_wp]
      ref(:, 152) = [5.417557912_wp, 0.185016068_wp, -0.009894340_wp, &
         & 2.81304160679008_wp, 0.68551453121954_wp, 0.066386781_wp, -0.007036174_wp]
      ref(:, 153) = [5.417557912_wp, -3.342652685_wp, -0.009894340_wp, &
         & 2.81304160679008_wp, 0.07159690383542_wp, 0.018219452_wp, -0.000077333_wp]
      ref(:, 154) = [7.181392289_wp, -1.578818309_wp, 3.517774412_wp, &
         & 2.81304160679008_wp, 0.85677081688370_wp, 0.030450327_wp, -0.008180940_wp]
      ref(:, 155) = [7.181392289_wp, -1.578818309_wp, -0.009894340_wp, &
         & 2.81304160679008_wp, 0.99721437616018_wp, 0.063476909_wp, -0.013237884_wp]
      ref(:, 156) = [7.181392289_wp, 0.185016068_wp, 1.753940036_wp, &
         & 2.81304160679008_wp, 0.99740417651070_wp, 0.060619435_wp, -0.013263902_wp]
      ref(:, 157) = [7.181392289_wp, -3.342652685_wp, 1.753940036_wp, &
         & 2.81304160679008_wp, 0.85091090905656_wp, 0.032968640_wp, -0.007916556_wp]
      ref(:, 158) = [6.857722650_wp, -0.138653571_wp, 3.194104774_wp, &
         & 3.06245180526872_wp, 0.93635966880702_wp, 0.038421905_wp, -0.008110117_wp]
      ref(:, 159) = [6.857722650_wp, -3.018983046_wp, 3.194104774_wp, &
         & 3.06245180526872_wp, 0.43017101421689_wp, 0.012784064_wp, -0.001678410_wp]
      ref(:, 160) = [6.857722650_wp, -0.138653571_wp, 0.313775299_wp, &
         & 3.06245180526872_wp, 0.99875701365517_wp, 0.067561387_wp, -0.011312015_wp]
      ref(:, 161) = [6.857722650_wp, -3.018983046_wp, 0.313775299_wp, &
         & 3.06245180526872_wp, 0.93062105192083_wp, 0.044002985_wp, -0.008030601_wp]
      ref(:, 162) = [0.259809803_wp, -2.138565061_wp, 4.109222878_wp, &
         & 2.51606090252714_wp, 0.00234630726808_wp, 0.042129077_wp, -0.000033225_wp]
      ref(:, 163) = [-4.729067190_wp, -2.138565061_wp, 4.109222878_wp, &
         & 2.51606090252714_wp, 0.25114424002685_wp, 0.040347328_wp, -0.001857885_wp]
      ref(:, 164) = [-2.234628694_wp, -4.633003558_wp, 4.109222878_wp, &
         & 2.51606090252714_wp, 0.04235175721088_wp, 0.053470619_wp, -0.000315408_wp]
      ref(:, 165) = [-2.234628694_wp, -2.138565061_wp, 6.603661374_wp, &
         & 2.51606090252714_wp, 0.99931551422751_wp, 0.042445556_wp, -0.016810004_wp]
      ref(:, 166) = [-2.234628694_wp, -0.374730685_wp, 5.873057254_wp, &
         & 2.81304160679008_wp, 0.32619030980113_wp, 0.010026434_wp, -0.002032035_wp]
      ref(:, 167) = [-2.234628694_wp, -3.902399437_wp, 5.873057254_wp, &
         & 2.81304160679008_wp, 0.97754916782438_wp, 0.050657193_wp, -0.010593427_wp]
      ref(:, 168) = [-0.470794317_wp, -2.138565061_wp, 5.873057254_wp, &
         & 2.81304160679008_wp, 0.96936423037214_wp, 0.036675263_wp, -0.010575866_wp]
      ref(:, 169) = [-3.998463070_wp, -2.138565061_wp, 5.873057254_wp, &
         & 2.81304160679008_wp, 0.99584229804347_wp, 0.035296323_wp, -0.010787311_wp]
      ref(:, 170) = [-0.470794317_wp, -0.374730685_wp, 4.109222878_wp, &
         & 2.81304160679008_wp, 0.00226243279906_wp, 0.044965858_wp, -0.000030685_wp]
      ref(:, 171) = [-0.470794317_wp, -3.902399437_wp, 4.109222878_wp, &
         & 2.81304160679008_wp, 0.01003588252964_wp, 0.023684925_wp, -0.000027253_wp]
      ref(:, 172) = [-3.998463070_wp, -3.902399437_wp, 4.109222878_wp, &
         & 2.81304160679008_wp, 0.08984362851396_wp, 0.068774124_wp, -0.000668784_wp]
      ref(:, 173) = [-0.794463956_wp, -0.698400323_wp, 5.549387615_wp, &
         & 3.06245180526872_wp, 0.95312776440853_wp, 0.031305571_wp, -0.010186792_wp]
      ref(:, 174) = [-3.674793431_wp, -0.698400323_wp, 5.549387615_wp, &
         & 3.06245180526872_wp, 0.00732036793236_wp, -0.007336751_wp, 0.000014354_wp]
      ref(:, 175) = [-0.794463956_wp, -3.578729799_wp, 5.549387615_wp, &
         & 3.06245180526872_wp, 0.88186619807911_wp, 0.038194040_wp, -0.006613116_wp]
      ref(:, 176) = [-3.674793431_wp, -3.578729799_wp, 5.549387615_wp, &
         & 3.06245180526872_wp, 0.97293321531912_wp, 0.053533400_wp, -0.009529674_wp]
      ref(:, 177) = [4.984083548_wp, -3.219521561_wp, -3.360509647_wp, &
         & 1.58152399587420_wp, 0.50514758266420_wp, 0.009904481_wp, 0.000219482_wp]
      ref(:, 178) = [-2.952766214_wp, -3.219521561_wp, -3.360509647_wp, &
         & 1.58152399587420_wp, 0.99991631808354_wp, 0.020851445_wp, 0.006802889_wp]
      ref(:, 179) = [1.015658667_wp, -7.187946443_wp, -3.360509647_wp, &
         & 1.58152399587420_wp, 0.99999999981714_wp, 0.006625792_wp, 0.002989027_wp]
      ref(:, 180) = [1.015658667_wp, -3.219521561_wp, -7.328934528_wp, &
         & 1.58152399587420_wp, 0.94662084379630_wp, 0.010783696_wp, 0.000539601_wp]
      ref(:, 181) = [1.015658667_wp, -6.025621706_wp, -0.554409503_wp, &
         & 1.76819758141091_wp, 0.93303999504620_wp, 0.024985456_wp, -0.000450564_wp]
      ref(:, 182) = [1.015658667_wp, -0.413421417_wp, -6.166609791_wp, &
         & 1.76819758141091_wp, 0.00000511847737_wp, 0.053886717_wp, -0.000000125_wp]
      ref(:, 183) = [1.015658667_wp, -6.025621706_wp, -6.166609791_wp, &
         & 1.76819758141091_wp, 0.99999999969107_wp, 0.003894287_wp, 0.002934583_wp]
      ref(:, 184) = [3.821758811_wp, -3.219521561_wp, -0.554409503_wp, &
         & 1.76819758141091_wp, 0.00421173770371_wp, 0.044664887_wp, -0.000088517_wp]
      ref(:, 185) = [-1.790441477_wp, -3.219521561_wp, -0.554409503_wp, &
         & 1.76819758141091_wp, 0.00290087377644_wp, 0.083896171_wp, -0.000101521_wp]
      ref(:, 186) = [3.821758811_wp, -3.219521561_wp, -6.166609791_wp, &
         & 1.76819758141091_wp, 0.52096477959209_wp, 0.003862755_wp, 0.001151007_wp]
      ref(:, 187) = [-1.790441477_wp, -3.219521561_wp, -6.166609791_wp, &
         & 1.76819758141091_wp, 0.99995071102360_wp, 0.011461842_wp, 0.005762863_wp]
      ref(:, 188) = [-1.790441477_wp, -0.413421417_wp, -3.360509647_wp, &
         & 1.76819758141091_wp, 0.00203974880639_wp, 0.056214787_wp, -0.000043029_wp]
      ref(:, 189) = [3.821758811_wp, -6.025621706_wp, -3.360509647_wp, &
         & 1.76819758141091_wp, 0.99999993882157_wp, -0.003594745_wp, 0.008788908_wp]
      ref(:, 190) = [-1.790441477_wp, -6.025621706_wp, -3.360509647_wp, &
         & 1.76819758141091_wp, 0.99999999987835_wp, 0.015814609_wp, 0.002290434_wp]
      ref(:, 191) = [3.306829841_wp, -5.510692735_wp, -1.069338473_wp, &
         & 1.92496970616891_wp, 0.98069021869565_wp, 0.004581507_wp, 0.002249596_wp]
      ref(:, 192) = [-1.275512506_wp, -5.510692735_wp, -1.069338473_wp, &
         & 1.92496970616891_wp, 0.98924892254265_wp, 0.033760663_wp, -0.003549800_wp]
      ref(:, 193) = [-1.275512506_wp, -0.928350388_wp, -5.651680820_wp, &
         & 1.92496970616891_wp, 0.33623529207419_wp, 0.029170410_wp, -0.000226570_wp]
      ref(:, 194) = [3.306829841_wp, -5.510692735_wp, -5.651680820_wp, &
         & 1.92496970616891_wp, 0.99999894419809_wp, -0.001178342_wp, 0.005133414_wp]
      ref(:, 195) = [-1.275512506_wp, -5.510692735_wp, -5.651680820_wp, &
         & 1.92496970616891_wp, 0.99999999999993_wp, 0.007810554_wp, 0.003154720_wp]
      ref(:, 196) = [6.775121582_wp, 0.266264352_wp, -3.918624763_wp, &
         & 1.44149322540618_wp, 0.58594438229965_wp, 0.009433236_wp, 0.004194215_wp]
      ref(:, 197) = [-1.932736443_wp, 0.266264352_wp, -3.918624763_wp, &
         & 1.44149322540618_wp, 0.00089896130279_wp, 0.078652101_wp, -0.000038832_wp]
      ref(:, 198) = [2.421192569_wp, 4.620193365_wp, -3.918624763_wp, &
         & 1.44149322540618_wp, 0.15372860903672_wp, 0.055768995_wp, -0.004022466_wp]
      ref(:, 199) = [2.421192569_wp, -4.087664660_wp, -3.918624763_wp, &
         & 1.44149322540618_wp, 0.00000285238487_wp, 0.386648592_wp, -0.000000935_wp]
      ref(:, 200) = [2.421192569_wp, 0.266264352_wp, 0.435304250_wp, &
         & 1.44149322540618_wp, 0.00000116123699_wp, 0.398091819_wp, -0.000000379_wp]
      ref(:, 201) = [2.421192569_wp, 0.266264352_wp, -8.272553775_wp, &
         & 1.44149322540618_wp, 0.99332341695172_wp, 0.008569672_wp, 0.007632124_wp]
      ref(:, 202) = [2.421192569_wp, 3.344957082_wp, -0.839932033_wp, &
         & 1.61163842055682_wp, 0.04231484084804_wp, 0.055137975_wp, -0.000996353_wp]
      ref(:, 203) = [2.421192569_wp, 3.344957082_wp, -6.997317492_wp, &
         & 1.61163842055682_wp, 0.93853145674244_wp, 0.021960392_wp, 0.001374117_wp]
      ref(:, 204) = [2.421192569_wp, -2.812428377_wp, -6.997317492_wp, &
         & 1.61163842055682_wp, 0.45608030236155_wp, 0.010953172_wp, -0.000766659_wp]
      ref(:, 205) = [5.499885299_wp, 0.266264352_wp, -0.839932033_wp, &
         & 1.61163842055682_wp, 0.76630113831579_wp, 0.049726846_wp, -0.006681247_wp]
      ref(:, 206) = [5.499885299_wp, 0.266264352_wp, -6.997317492_wp, &
         & 1.61163842055682_wp, 0.35549380208991_wp, 0.006801355_wp, 0.000448569_wp]
      ref(:, 207) = [-0.657500160_wp, 0.266264352_wp, -6.997317492_wp, &
         & 1.61163842055682_wp, 0.62214257441247_wp, 0.029108931_wp, -0.000162972_wp]
      ref(:, 208) = [5.499885299_wp, 3.344957082_wp, -3.918624763_wp, &
         & 1.61163842055682_wp, 0.99569400384091_wp, 0.018905385_wp, -0.000843105_wp]
      ref(:, 209) = [5.499885299_wp, -2.812428377_wp, -3.918624763_wp, &
         & 1.61163842055682_wp, 0.85666792029421_wp, 0.007097048_wp, 0.004299430_wp]
      ref(:, 210) = [-0.657500160_wp, -2.812428377_wp, -3.918624763_wp, &
         & 1.61163842055682_wp, 0.00000041906065_wp, 0.333810726_wp, -0.000000103_wp]
      ref(:, 211) = [4.934934657_wp, 2.780006440_wp, -1.404882675_wp, &
         & 1.75452968010187_wp, 0.80018922293848_wp, 0.030863125_wp, -0.002161393_wp]
      ref(:, 212) = [4.934934657_wp, -2.247477735_wp, -1.404882675_wp, &
         & 1.75452968010187_wp, 0.48476451058773_wp, 0.043917148_wp, -0.007781672_wp]
      ref(:, 213) = [4.934934657_wp, 2.780006440_wp, -6.432366850_wp, &
         & 1.75452968010187_wp, 0.95384718356785_wp, 0.009217608_wp, 0.004977683_wp]
      ref(:, 214) = [-0.092549518_wp, 2.780006440_wp, -6.432366850_wp, &
         & 1.75452968010187_wp, 0.00834047467973_wp, 0.064887920_wp, -0.000172177_wp]
      ref(:, 215) = [4.934934657_wp, -2.247477735_wp, -6.432366850_wp, &
         & 1.75452968010187_wp, 0.74381405233053_wp, 0.001633283_wp, 0.003868090_wp]
      ref(:, 216) = [-0.092549518_wp, -2.247477735_wp, -6.432366850_wp, &
         & 1.75452968010187_wp, 0.08188624502532_wp, 0.020609459_wp, -0.000267163_wp]
      ref(:, 217) = [1.328668009_wp, 2.536678903_wp, 2.316649859_wp, &
         & 1.44149322540618_wp, 0.46415591980070_wp, 0.011601962_wp, -0.003818037_wp]
      ref(:, 218) = [-7.379190016_wp, 2.536678903_wp, 2.316649859_wp, &
         & 1.44149322540618_wp, 0.89196001392740_wp, -0.036287465_wp, 0.007271805_wp]
      ref(:, 219) = [-3.025261003_wp, 6.890607916_wp, 2.316649859_wp, &
         & 1.44149322540618_wp, 0.01327222513480_wp, -0.177266874_wp, 0.001013433_wp]
      ref(:, 220) = [-3.025261003_wp, 2.536678903_wp, 6.670578871_wp, &
         & 1.44149322540618_wp, 0.04951288636591_wp, -0.082600159_wp, 0.001423059_wp]
      ref(:, 221) = [-3.025261003_wp, 5.615371633_wp, 5.395342588_wp, &
         & 1.61163842055682_wp, 0.90858128150109_wp, -0.084478934_wp, 0.018601585_wp]
      ref(:, 222) = [-3.025261003_wp, -0.542013826_wp, 5.395342588_wp, &
         & 1.61163842055682_wp, 0.01967563462713_wp, 0.012193924_wp, -0.000182741_wp]
      ref(:, 223) = [-3.025261003_wp, -0.542013826_wp, -0.762042871_wp, &
         & 1.61163842055682_wp, 0.00022011458150_wp, -0.016538986_wp, 0.000006810_wp]
      ref(:, 224) = [0.053431727_wp, 2.536678903_wp, 5.395342588_wp, &
         & 1.61163842055682_wp, 0.99229931289734_wp, -0.017881687_wp, -0.003339297_wp]
      ref(:, 225) = [-6.103953733_wp, 2.536678903_wp, 5.395342588_wp, &
         & 1.61163842055682_wp, 0.00457616978210_wp, -0.098313364_wp, 0.000140785_wp]
      ref(:, 226) = [-6.103953733_wp, 2.536678903_wp, -0.762042871_wp, &
         & 1.61163842055682_wp, 0.13977901235252_wp, 0.004569194_wp, -0.000352845_wp]
      ref(:, 227) = [0.053431727_wp, 5.615371633_wp, 2.316649859_wp, &
         & 1.61163842055682_wp, 0.30750017064463_wp, -0.079004273_wp, 0.010497218_wp]
      ref(:, 228) = [-6.103953733_wp, 5.615371633_wp, 2.316649859_wp, &
         & 1.61163842055682_wp, 0.45406844621640_wp, -0.108765536_wp, 0.016932121_wp]
      ref(:, 229) = [-6.103953733_wp, -0.542013826_wp, 2.316649859_wp, &
         & 1.61163842055682_wp, 0.82521236977836_wp, 0.004774797_wp, 0.000238172_wp]
      ref(:, 230) = [-0.511518916_wp, 5.050420991_wp, 4.830391946_wp, &
         & 1.75452968010187_wp, 0.98707997693871_wp, -0.052487397_wp, 0.008186952_wp]
      ref(:, 231) = [-5.539003091_wp, 5.050420991_wp, 4.830391946_wp, &
         & 1.75452968010187_wp, 0.64001076605148_wp, -0.093808289_wp, 0.013066520_wp]
      ref(:, 232) = [-0.511518916_wp, 0.022936816_wp, 4.830391946_wp, &
         & 1.75452968010187_wp, 0.51330041073732_wp, 0.019480321_wp, -0.005583580_wp]
      ref(:, 233) = [-5.539003091_wp, 0.022936816_wp, 4.830391946_wp, &
         & 1.75452968010187_wp, 0.00424642075355_wp, -0.060675874_wp, 0.000104404_wp]
      ref(:, 234) = [-5.539003091_wp, 5.050420991_wp, -0.197092229_wp, &
         & 1.75452968010187_wp, 0.00752937315771_wp, -0.085959269_wp, 0.000301792_wp]
      ref(:, 235) = [-5.539003091_wp, 0.022936816_wp, -0.197092229_wp, &
         & 1.75452968010187_wp, 0.76317353512283_wp, 0.013259712_wp, 0.003089472_wp]
      ref(:, 236) = [-5.519280106_wp, -2.292351381_wp, 2.197828084_wp, &
         & 1.78559160824507_wp, 0.30204039362416_wp, 0.060388022_wp, -0.005477750_wp]
      ref(:, 237) = [-2.004389497_wp, -5.807241990_wp, 2.197828084_wp, &
         & 1.78559160824507_wp, 0.94439663195767_wp, 0.033083869_wp, -0.000252831_wp]
      ref(:, 238) = [-2.004389497_wp, -2.292351381_wp, 5.712718693_wp, &
         & 1.78559160824507_wp, 0.01426445530323_wp, 0.122687739_wp, -0.000883774_wp]
      ref(:, 239) = [-2.004389497_wp, -2.292351381_wp, -1.317062525_wp, &
         & 1.78559160824507_wp, 0.02062412554830_wp, 0.061518019_wp, -0.000459361_wp]
      ref(:, 240) = [-2.004389497_wp, 0.193051604_wp, 4.683231069_wp, &
         & 1.99635210804458_wp, 0.00009361384673_wp, 0.014839289_wp, -0.000000902_wp]
      ref(:, 241) = [-2.004389497_wp, -4.777754366_wp, 4.683231069_wp, &
         & 1.99635210804458_wp, 0.72921822280349_wp, 0.041730630_wp, -0.006225216_wp]
      ref(:, 242) = [-2.004389497_wp, -4.777754366_wp, -0.287574900_wp, &
         & 1.99635210804458_wp, 0.73048247336994_wp, 0.049614366_wp, -0.006330013_wp]
      ref(:, 243) = [0.481013488_wp, -2.292351381_wp, 4.683231069_wp, &
         & 1.99635210804458_wp, 0.52325801852375_wp, 0.019294261_wp, -0.004992345_wp]
      ref(:, 244) = [-4.489792481_wp, -2.292351381_wp, 4.683231069_wp, &
         & 1.99635210804458_wp, 0.31890841588113_wp, 0.045558304_wp, -0.004330102_wp]
      ref(:, 245) = [-4.489792481_wp, -2.292351381_wp, -0.287574900_wp, &
         & 1.99635210804458_wp, 0.17569196758670_wp, 0.071319126_wp, -0.002865723_wp]
      ref(:, 246) = [0.481013488_wp, -4.777754366_wp, 2.197828084_wp, &
         & 1.99635210804458_wp, 0.19252442913454_wp, 0.023815136_wp, -0.001085929_wp]
      ref(:, 247) = [-4.489792481_wp, -4.777754366_wp, 2.197828084_wp, &
         & 1.99635210804458_wp, 0.07743360505699_wp, 0.098927139_wp, -0.001984214_wp]
      ref(:, 248) = [0.024933543_wp, -0.263028342_wp, 4.227151124_wp, &
         & 2.17335289406167_wp, 0.17000271296144_wp, 0.032343029_wp, -0.002346633_wp]
      ref(:, 249) = [-4.033712536_wp, -0.263028342_wp, 4.227151124_wp, &
         & 2.17335289406167_wp, 0.00000040702456_wp, -0.013451320_wp, 0.000000002_wp]
      ref(:, 250) = [0.024933543_wp, -4.321674420_wp, 4.227151124_wp, &
         & 2.17335289406167_wp, 0.96381335373141_wp, 0.008794766_wp, 0.001945100_wp]
      ref(:, 251) = [-4.033712536_wp, -4.321674420_wp, 4.227151124_wp, &
         & 2.17335289406167_wp, 0.80685866422043_wp, 0.057718403_wp, -0.010301980_wp]
      ref(:, 252) = [-4.033712536_wp, -0.263028342_wp, 0.168505045_wp, &
         & 2.17335289406167_wp, 0.00000963480489_wp, -0.004649642_wp, 0.000000131_wp]
      ref(:, 253) = [0.024933543_wp, -4.321674420_wp, 0.168505045_wp, &
         & 2.17335289406167_wp, 0.00032675416448_wp, 0.104956015_wp, -0.000014458_wp]
      ref(:, 254) = [-4.033712536_wp, -4.321674420_wp, 0.168505045_wp, &
         & 2.17335289406167_wp, 0.01501536484590_wp, 0.121960186_wp, -0.000518863_wp]
      ref(:, 255) = [5.294780850_wp, -1.369420077_wp, 0.484550557_wp, &
         & 1.50416684390210_wp, 0.00008760541353_wp, 0.185482086_wp, -0.000010859_wp]
      ref(:, 256) = [-3.050249757_wp, -1.369420077_wp, 0.484550557_wp, &
         & 1.50416684390210_wp, 0.00000083126587_wp, 0.107677811_wp, -0.000000056_wp]
      ref(:, 257) = [1.122265547_wp, 2.803095227_wp, 0.484550557_wp, &
         & 1.50416684390210_wp, 0.01126036312408_wp, 0.043613979_wp, -0.000271261_wp]
      ref(:, 258) = [1.122265547_wp, -5.541935381_wp, 0.484550557_wp, &
         & 1.50416684390210_wp, 0.86181694145560_wp, 0.036005728_wp, -0.023338431_wp]
      ref(:, 259) = [1.122265547_wp, -1.369420077_wp, 4.657065861_wp, &
         & 1.50416684390210_wp, 0.56726977849291_wp, 0.003281321_wp, -0.002276729_wp]
      ref(:, 260) = [1.122265547_wp, 1.580993789_wp, 3.434964423_wp, &
         & 1.68170965623320_wp, 0.45467895972205_wp, 0.008327237_wp, -0.003745642_wp]
      ref(:, 261) = [1.122265547_wp, -4.319833943_wp, 3.434964423_wp, &
         & 1.68170965623320_wp, 0.57250749202404_wp, -0.004934376_wp, 0.000060528_wp]
      ref(:, 262) = [4.072679412_wp, -1.369420077_wp, 3.434964423_wp, &
         & 1.68170965623320_wp, 0.00000018142518_wp, 0.124188160_wp, -0.000000018_wp]
      ref(:, 263) = [4.072679412_wp, -1.369420077_wp, -2.465863309_wp, &
         & 1.68170965623320_wp, 0.00001049687484_wp, 0.086473470_wp, -0.000000496_wp]
      ref(:, 264) = [-1.828148319_wp, -1.369420077_wp, -2.465863309_wp, &
         & 1.68170965623320_wp, 0.03126043446735_wp, 0.040817782_wp, -0.000354453_wp]
      ref(:, 265) = [4.072679412_wp, 1.580993789_wp, 0.484550557_wp, &
         & 1.68170965623320_wp, 0.53959323121149_wp, 0.057423079_wp, -0.012242749_wp]
      ref(:, 266) = [4.072679412_wp, -4.319833943_wp, 0.484550557_wp, &
         & 1.68170965623320_wp, 0.08190584357246_wp, -0.023163473_wp, 0.001498317_wp]
      ref(:, 267) = [-1.828148319_wp, -4.319833943_wp, 0.484550557_wp, &
         & 1.68170965623320_wp, 0.00625191135976_wp, 0.059141046_wp, -0.000088846_wp]
      ref(:, 268) = [3.531268380_wp, 1.039582757_wp, 2.893553391_wp, &
         & 1.83081357923673_wp, 0.24730744975798_wp, 0.000303486_wp, 0.001237826_wp]
      ref(:, 269) = [3.531268380_wp, -3.778422911_wp, 2.893553391_wp, &
         & 1.83081357923673_wp, 0.00001554372045_wp, -0.069550717_wp, 0.000000523_wp]
      ref(:, 270) = [-1.286737287_wp, -3.778422911_wp, 2.893553391_wp, &
         & 1.83081357923673_wp, 0.00000077880333_wp, 0.130640225_wp, -0.000000055_wp]
      ref(:, 271) = [3.531268380_wp, -3.778422911_wp, -1.924452277_wp, &
         & 1.83081357923673_wp, 0.00322879967805_wp, 0.035091494_wp, -0.000049599_wp]
      ref(:, 272) = [-1.286737287_wp, -3.778422911_wp, -1.924452277_wp, &
         & 1.83081357923673_wp, 0.00095424479095_wp, 0.065792319_wp, -0.000025192_wp]

      ! Populate iSwiG cavity with ORCA equiv. inputs
      cav%ngrid = npts
      cav%xyz = ref(1:3, :)
      cav%xi0 = ref(4, :)
      cav%f = ref(5, :)

      ! Scale ORCA's unscaled conductor charges
      feps = (eps_ref - 1.0_wp)/eps_ref
      phi = ref(6, :)
      q_ref = ref(7, :)*feps

      allocate (amat(npts, npts), q(npts))
      call assemble_pcm_amat(cav%xi0, cav%f, cav%xyz, amat, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      rhs = -feps*phi
      call solve_pcm_cholesky(amat, rhs, q, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ! Check charges; thr is number of reference decimals
      call check(error, maxval(abs(q - q_ref)), 0.0_wp, thr=5.0E-9_wp, &
         & more="iSwiG CPCM charges deviate from ORCA reference")
      if (allocated(error)) return

      ! Dielectric energy: E = 1/2 q^T phi (observed deviation 2.3e-10)
      e_ours = 0.5_wp*dot_product(q, phi)
      call check(error, e_ours, e_ref, thr=5.0E-9_wp, &
         & more="iSwiG CPCM dielectric energy deviates from ORCA reference")
      if (allocated(error)) return

   end subroutine test_amat_orca_reference

   !> Numerical vs analytical test for the contracted A-matrix derivative.
   !> Uses 5-point finite differences on q1^T A q2 to verify the
   !> `pcm_amat_surface_weights` + `pcm_amat_nuclear_gradient` chain.
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
      real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :)
      real(wp), allocatable :: ana_grad(:, :), num_grad(:, :)
      real(wp) :: ffwd, fwd, bwd, bbwd
      integer :: i, n, iat, iax
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "13")
      allocate (radii(mol%nat))
      radii = 4.0_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ! Build reference cavity to set up charge vectors
      allocate (cav)
      call new_cavity_iswig(cav, ctx, cut_f=0.01_wp, radius_model=radius_model, &
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
      allocate (q1(n), q2(n))
      do i = 1, n
         q1(i) = real(i, wp)/real(n + 1, wp)
         if (mod(i, 2) == 0) then
            q2(i) = real(i, wp)/real(n + 1, wp)
         else
            q2(i) = -real(i, wp)/real(n + 1, wp)
         end if
      end do

      ! Compute analytical gradient through the generic PCM A-matrix adjoints
      call cav%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate (w_xi(n), w_f(n), w_xyz(3, n))
      call pcm_amat_surface_weights(cav%xi0, cav%f, cav%xyz, q1, q2, &
         & w_xi, w_f, w_xyz, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate (ana_grad(3, mol%nat))
      call pcm_amat_nuclear_gradient(cav%xi1_rA, cav%f1_rA, cav%xyz1_rA, &
         & w_xi, w_f, w_xyz, ana_grad, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ! Compute numerical gradient via 5-point stencil
      allocate (num_grad(3, mol%nat))
      num_grad = 0.0_wp

      do iat = 1, mol%nat
         do iax = 1, 3

            ! +2h
            mol%xyz(iax, iat) = mol%xyz(iax, iat) + 2.0_wp*STEP_SIZE
            call amat_energy(ffwd)
            if (allocated(error)) return

            ! +h
            mol%xyz(iax, iat) = mol%xyz(iax, iat) - STEP_SIZE
            call amat_energy(fwd)
            if (allocated(error)) return

            ! -h
            mol%xyz(iax, iat) = mol%xyz(iax, iat) - 2.0_wp*STEP_SIZE
            call amat_energy(bwd)
            if (allocated(error)) return

            ! -2h
            mol%xyz(iax, iat) = mol%xyz(iax, iat) - STEP_SIZE
            call amat_energy(bbwd)
            if (allocated(error)) return

            ! Restore position
            mol%xyz(iax, iat) = mol%xyz(iax, iat) + 2.0_wp*STEP_SIZE

            ! 5-point stencil
            num_grad(iax, iat) = (-ffwd + 8.0_wp*fwd &
               & - 8.0_wp*bwd + bbwd)/(12.0_wp*STEP_SIZE)
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

   contains

      !> Rebuild the cavity at the current geometry and return q1^T A q2
      subroutine amat_energy(value)
         !> Contracted A-matrix energy
         real(wp), intent(out) :: value

         value = 0.0_wp
         if (allocated(cav)) deallocate (cav)
         allocate (cav)
         call new_cavity_iswig(cav, ctx, cut_f=0.01_wp, radius_model=radius_model, &
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
         if (cav%ngrid /= n) then
            call test_failed(error, "A-matrix gradient FD: grid size changed")
            return
         end if
         if (allocated(amat)) deallocate (amat)
         allocate (amat(cav%ngrid, cav%ngrid))
         call assemble_pcm_amat(cav%xi0, cav%f, cav%xyz, amat, cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if
         value = dot_product(q1, matmul(amat, q2))
         deallocate (amat)
      end subroutine amat_energy

   end subroutine test_amat_gradient

   !> Reverse-mode surface gradient against finite differences
   !>
   !> [[get_surface_gradient_iswig]] contracts an accumulated surface adjoint
   !> straight into `dE/dR_A`. The reference is the numerical derivative of the
   !> very quantity that adjoint defines,
   !>
   !>   E = sum_i [ w_xi_i xi_i + w_f_i f_i + w_a_i a_i + w_w_i wleb_i
   !>             + w_xyz_i . r_i + w_n_i . n_i + (w_k1_i + w_k2_i)/R_owner(i) ]
   !>
   !> re-evaluated on a cavity rebuilt at each displaced geometry. Every channel
   !> is driven at once, including the four the contraction drops as
   !> geometry-independent: if any of them did move with the nuclei, the
   !> numerical derivative would see it and the analytic one would not.
   !>
   !> Weights are pinned to the *raw* Lebedev index rather than the filtered
   !> grid index, so a rebuilt grid keeps assigning the same weight to the same
   !> point; the energy helper additionally refuses to run if the surviving set
   !> changes at all, which would put a step into the finite difference.
   subroutine test_surface_gradient(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      type(cavity_surface_adjoint_type) :: acc
      real(wp), allocatable :: radii(:)
      !> Adjoint weights, indexed by raw (pre-filter) grid point
      real(wp), allocatable :: rw_xi(:), rw_f(:), rw_a(:), rw_w(:)
      real(wp), allocatable :: rw_xyz(:, :), rw_n(:, :), rw_k1(:), rw_k2(:)
      !> Channel weights on the filtered grid
      real(wp), allocatable :: w_xi(:), w_f(:), w_a(:), w_w(:)
      real(wp), allocatable :: w_xyz(:, :), w_n(:, :), w_k1(:), w_k2(:)
      !> Reverse-mode and numerical gradients
      real(wp), allocatable :: ana_grad(:, :), num_grad(:, :)
      !> Reference grid extent and point identities
      integer, allocatable :: numbering_ref(:)
      integer :: nraw, n
      real(wp) :: ffwd, fwd, bwd, bbwd
      integer :: i, iraw, iat, iax
      !> Number of Lebedev points per sphere of the fixture
      integer, parameter :: NLEB = 50
      !> Switching cutoff, low enough that the dropped points carry no weight
      real(wp), parameter :: CUT_F = 1.0E-7_wp
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx)

      call get_structure(mol, "MB16-43", "01")
      allocate (radii(mol%nat))
      radii = 2.0_wp

      call new_radii_custom_atoms(radii, radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate (cav)
      call new_cavity_iswig(cav, ctx, nleb=NLEB, cut_f=CUT_F, &
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

      n = cav%ngrid
      nraw = mol%nat*NLEB
      numbering_ref = cav%numbering

      ! Deterministic weights of comparable size in every channel, scaled so the
      ! grid sum stays of order one
      allocate (rw_xi(nraw), rw_f(nraw), rw_a(nraw), rw_w(nraw))
      allocate (rw_k1(nraw), rw_k2(nraw))
      allocate (rw_xyz(3, nraw), rw_n(3, nraw))
      do iraw = 1, nraw
         rw_xi(iraw) = channel_weight(iraw, 1)
         rw_f(iraw) = channel_weight(iraw, 2)
         rw_a(iraw) = channel_weight(iraw, 3)
         rw_w(iraw) = channel_weight(iraw, 4)
         rw_k1(iraw) = channel_weight(iraw, 5)
         rw_k2(iraw) = channel_weight(iraw, 6)
         do iax = 1, 3
            rw_xyz(iax, iraw) = channel_weight(iraw, 6 + iax)
            rw_n(iax, iraw) = channel_weight(iraw, 9 + iax)
         end do
      end do

      ! Analytical gradient through the reverse-mode contraction
      allocate (w_xi(n), w_f(n), w_a(n), w_w(n), w_k1(n), w_k2(n))
      allocate (w_xyz(3, n), w_n(3, n))
      do i = 1, n
         iraw = cav%numbering(i)
         w_xi(i) = rw_xi(iraw)
         w_f(i) = rw_f(iraw)
         w_a(i) = rw_a(iraw)
         w_w(i) = rw_w(iraw)
         w_k1(i) = rw_k1(iraw)
         w_k2(i) = rw_k2(iraw)
         w_xyz(:, i) = rw_xyz(:, iraw)
         w_n(:, i) = rw_n(:, iraw)
      end do

      call acc%init(n)
      call acc%add_surface_weights(cavity_error, w_xi=w_xi, w_f=w_f, w_a=w_a, &
         & w_w=w_w, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, w_k2=w_k2)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate (ana_grad(3, mol%nat), source=0.0_wp)
      call cav%get_surface_gradient(acc, ana_grad, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ! Numerical gradient via 5-point stencil
      allocate (num_grad(3, mol%nat), source=0.0_wp)
      do iat = 1, mol%nat
         do iax = 1, 3

            mol%xyz(iax, iat) = mol%xyz(iax, iat) + 2.0_wp*STEP_SIZE
            call adjoint_energy(ffwd)
            if (allocated(error)) return

            mol%xyz(iax, iat) = mol%xyz(iax, iat) - STEP_SIZE
            call adjoint_energy(fwd)
            if (allocated(error)) return

            mol%xyz(iax, iat) = mol%xyz(iax, iat) - 2.0_wp*STEP_SIZE
            call adjoint_energy(bwd)
            if (allocated(error)) return

            mol%xyz(iax, iat) = mol%xyz(iax, iat) - STEP_SIZE
            call adjoint_energy(bbwd)
            if (allocated(error)) return

            mol%xyz(iax, iat) = mol%xyz(iax, iat) + 2.0_wp*STEP_SIZE

            num_grad(iax, iat) = (-ffwd + 8.0_wp*fwd &
               & - 8.0_wp*bwd + bbwd)/(12.0_wp*STEP_SIZE)
         end do
      end do

      ! A vanishing reference would let any implementation pass
      if (maxval(abs(num_grad)) <= 1.0E-6_wp) then
         call test_failed(error, "Surface-gradient reference is vacuous")
         return
      end if

      do iat = 1, mol%nat
         do iax = 1, 3
            call check(error, ana_grad(iax, iat), num_grad(iax, iat), &
               & thr_abs=ABS_THR, thr_rel=REL_THR, &
               & more="Reverse-mode surface gradient mismatch")
            if (allocated(error)) return
         end do
      end do

   contains

      !> Deterministic pseudo-random weight for one raw point and channel
      pure function channel_weight(ipoint, ichannel) result(weight)
         !> Raw grid point index
         integer, intent(in) :: ipoint
         !> Channel identifier
         integer, intent(in) :: ichannel
         !> Weight of this point in this channel
         real(wp) :: weight

         weight = sin(0.37_wp*real(ipoint, wp) + 1.7_wp*real(ichannel, wp)) &
            & /real(nraw, wp)

      end function channel_weight

      !> Rebuild the cavity at the current geometry and return the adjoint energy
      subroutine adjoint_energy(value)
         !> Surface energy whose adjoints `acc` holds
         real(wp), intent(out) :: value

         integer :: igrid, jraw, jat

         value = 0.0_wp
         if (allocated(cav)) deallocate (cav)
         allocate (cav)
         call new_cavity_iswig(cav, ctx, nleb=NLEB, cut_f=CUT_F, &
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
         if (cav%ngrid /= n) then
            call test_failed(error, "Surface-gradient FD: grid size changed")
            return
         end if
         if (any(cav%numbering /= numbering_ref)) then
            call test_failed(error, "Surface-gradient FD: surviving point set changed")
            return
         end if

         do igrid = 1, cav%ngrid
            jraw = cav%numbering(igrid)
            jat = cav%owner(igrid)
            value = value &
               & + rw_xi(jraw)*cav%xi0(igrid) &
               & + rw_f(jraw)*cav%f(igrid) &
               & + rw_a(jraw)*cav%a(igrid) &
               & + rw_w(jraw)*cav%wleb(igrid) &
               & + dot_product(rw_xyz(:, jraw), cav%xyz(:, igrid)) &
               & + dot_product(rw_n(:, jraw), cav%normal0(:, igrid)) &
               & + (rw_k1(jraw) + rw_k2(jraw))/cav%radii(jat)
         end do

      end subroutine adjoint_energy

   end subroutine test_surface_gradient

end module test_cavity_iswig
