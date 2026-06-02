module test_cavity_drop_convergence
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_marchingcubes, only: integrate_surface_marching_cubes
   use moist_cavity, only: cavity_type_iswig, new_cavity_iswig
   use moist_radii, only: default_cpcm_radii, new_radii_custom_atoms, radius_type
   use moist_data_radii_legacy, only: get_radius_func
   use mstore, only: get_structure
   implicit none
   private

   public :: collect_cavity_drop_convergence

   integer, parameter :: ndim = 3

   real(wp), parameter :: blend_k = 5.0_wp
   real(wp), parameter :: blend_2b = 1.0_wp
   real(wp), parameter :: blend_3b = 1.0_wp
   real(wp), parameter :: proj_tol = 1.0e-15_wp
   integer, parameter :: proj_maxiter = 1000
   integer, parameter :: proj_level = 2

contains

   subroutine collect_cavity_drop_convergence(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("convergence_drop_lebedev", test_convergence_drop_lebedev), &
         new_unittest("convergence_marching_cubes", test_convergence_marching_cubes), &
         new_unittest("convergence_iswig_lebedev", test_convergence_iswig_lebedev), &
         new_unittest("convergence_drop_gradient", test_convergence_drop_gradient), &
         new_unittest("convergence_drop_nleb_blendk", test_convergence_drop_nleb_blendk), &
         new_unittest("convergence_drop_nleb_projtol", test_convergence_drop_nleb_projtol) &
         ]
   end subroutine collect_cavity_drop_convergence

   !> Test convergence of DROP cavity area and volume w.r.t. Lebedev grid size.
   subroutine test_convergence_drop_lebedev(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error
      real(wp) :: areas(29), volumes(29)
      real(wp) :: ref_area, ref_volume
      integer :: igrid, imol

      integer, parameter :: n_grids = 29
      integer, parameter :: nleb_values(n_grids) = [ &
         6, 14, 26, 38, 50, 86, 110, 146, &
         170, 194, 302, 350, 434, 590, 770, 974, &
         1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, &
         3890, 4334, 4802, 5294, 5810]

      integer, parameter :: n_mols = 3
      character(len=12), parameter :: dataset_names(n_mols) = [ &
         'MB16-43     ', 'Amino20x4   ', 'UPU23       ']
      character(len=7), parameter :: mol_ids(n_mols) = [ &
         'CH4    ', 'THR_xab', '4b     ']

      write (*, '(a)') ''
      write (*, '(a)') '========================================================================'
      write (*, '(a)') 'Convergence: DROP cavity area & volume vs. Lebedev grid size'
      write (*, '(a)') '========================================================================'

      do imol = 1, n_mols
         call get_structure(mol, trim(dataset_names(imol)), trim(mol_ids(imol)))

         do igrid = 1, n_grids
            if (allocated(cavity)) deallocate(cavity)
            allocate(cavity)
            block
               type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
               call svdw_template%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
               call new_cavity_drop(cavity, nleb=nleb_values(igrid), &
                  tolerance=proj_tol, proj_maxiter=proj_maxiter, proj_level=proj_level, &
                  debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
                  lsf_model=svdw_template, error=cavity_error)
            end block
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if

            call cavity%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if

            areas(igrid) = cavity%total_area
            volumes(igrid) = cavity%total_volume
         end do

         ref_area = areas(n_grids)
         ref_volume = volumes(n_grids)

         write (*, '(a)') ''
         write (*, '(a,a,a,a,a,i0,a)') '  Molecule: ', trim(dataset_names(imol)), &
            '/', trim(mol_ids(imol)), ' (', mol%nat, ' atoms)'
         write (*, '(a8, 4a20)') &
            'nleb', 'Area (bohr^2)', 'Volume (bohr^3)', 'dA_ref (%)', 'dV_ref (%)'
         write (*, '(a8, 4a20)') &
            '-------', '-------------------', '-------------------', &
            '-------------------', '-------------------'

         do igrid = n_grids, 1, -1
            if (igrid == n_grids) then
               write (*, '(i8, 2f20.12, 2a20)') nleb_values(igrid), &
                  areas(igrid), volumes(igrid), '         ref        ', '         ref        '
            else
               write (*, '(i8, 4f20.12)') nleb_values(igrid), &
                  areas(igrid), volumes(igrid), &
                  100.0_wp * (areas(igrid) - ref_area) / ref_area, &
                  100.0_wp * (volumes(igrid) - ref_volume) / ref_volume
            end if
         end do

      end do

      write (*, '(a)') ''

      if (allocated(cavity)) deallocate(cavity)

   end subroutine test_convergence_drop_lebedev

   !> Test convergence of marching cubes area and volume w.r.t. grid spacing.
   subroutine test_convergence_marching_cubes(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_parameters_type) :: param
      real(wp), allocatable :: radii(:)
      real(wp) :: areas(7), volumes(7)
      real(wp) :: ref_area, ref_volume
      integer :: igrid, imol, iat

      integer, parameter :: n_spacings = 8
      real(wp), parameter :: spacings(n_spacings) = [ &
         1.0_wp, 0.5_wp, 0.3_wp, 0.2_wp, 0.1_wp, 0.05_wp, 0.03_wp, 0.01_wp]

      integer, parameter :: n_mols = 3
      character(len=12), parameter :: dataset_names(n_mols) = [ &
         'MB16-43     ', 'Amino20x4   ', 'UPU23       ']
      character(len=7), parameter :: mol_ids(n_mols) = [ &
         'CH4    ', 'THR_xab', '4b     ']

      write (*, '(a)') ''
      write (*, '(a)') '========================================================================'
      write (*, '(a)') 'Convergence: Marching cubes area & volume vs. grid spacing'
      write (*, '(a)') '========================================================================'

      do imol = 1, n_mols
         call get_structure(mol, trim(dataset_names(imol)), trim(mol_ids(imol)))

         if (allocated(radii)) deallocate(radii)
         allocate(radii(mol%nat))
         do iat = 1, mol%nat
            radii(iat) = get_radius_func(mol%num(mol%id(iat)))
         end do

         call lsf%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
         !> Without a cavity to set this, the direct user owns the
         !> screening threshold; lsf_update reads it when sizing SSD.
         lsf%screening_threshold = proj_tol * 0.1_wp
         call lsf%update(mol, radii)

         do igrid = 1, n_spacings
            call integrate_surface_marching_cubes(lsf, mol%xyz, areas(igrid), volumes(igrid), &
               target_spacing=spacings(igrid))
         end do

         ref_area = areas(n_spacings)
         ref_volume = volumes(n_spacings)

         write (*, '(a)') ''
         write (*, '(a,a,a,a,a,i0,a)') '  Molecule: ', trim(dataset_names(imol)), &
            '/', trim(mol_ids(imol)), ' (', mol%nat, ' atoms)'
         write (*, '(a12, 4a20)') &
            'spacing', 'Area (bohr^2)', 'Volume (bohr^3)', 'dA_ref (%)', 'dV_ref (%)'
         write (*, '(a12, 4a20)') &
            '-----------', '-------------------', '-------------------', &
            '-------------------', '-------------------'

         do igrid = n_spacings, 1, -1
            if (igrid == n_spacings) then
               write (*, '(f12.4, 2f20.12, 2a20)') spacings(igrid), &
                  areas(igrid), volumes(igrid), '         ref        ', '         ref        '
            else
               write (*, '(f12.4, 4f20.12)') spacings(igrid), &
                  areas(igrid), volumes(igrid), &
                  100.0_wp * (areas(igrid) - ref_area) / ref_area, &
                  100.0_wp * (volumes(igrid) - ref_volume) / ref_volume
            end if
         end do

      end do

      write (*, '(a)') ''

      if (allocated(radii)) deallocate(radii)

   end subroutine test_convergence_marching_cubes

   !> Test convergence of iSWiG cavity area and volume w.r.t. Lebedev grid size.
   subroutine test_convergence_iswig_lebedev(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cav
      type(mctc_error), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model
      real(wp), allocatable :: radii(:)
      real(wp) :: areas(11), volumes(11)
      real(wp) :: ref_area, ref_volume
      integer :: igrid, imol, iat

      integer, parameter :: n_grids = 11
      integer, parameter :: nleb_values(n_grids) = [ &
         14, 26, 50, 110, 194, 302, 434, 590, 770, 974, 1202]

      integer, parameter :: n_mols = 3
      character(len=12), parameter :: dataset_names(n_mols) = [ &
         'MB16-43     ', 'Amino20x4   ', 'UPU23       ']
      character(len=7), parameter :: mol_ids(n_mols) = [ &
         'CH4    ', 'THR_xab', '4b     ']

      write (*, '(a)') ''
      write (*, '(a)') '========================================================================'
      write (*, '(a)') 'Convergence: iSWiG cavity area & volume vs. Lebedev grid size'
      write (*, '(a)') '========================================================================'

      do imol = 1, n_mols
         call get_structure(mol, trim(dataset_names(imol)), trim(mol_ids(imol)))

         if (allocated(radii)) deallocate(radii)
         allocate(radii(mol%nat))
         do iat = 1, mol%nat
            radii(iat) = get_radius_func(mol%num(mol%id(iat)))
         end do

         call new_radii_custom_atoms(radii, radius_model, cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if

         do igrid = 1, n_grids
            if (allocated(cav)) deallocate(cav)
            allocate(cav)
            call new_cavity_iswig(cav, nleb=nleb_values(igrid), &
               radius_model=radius_model, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if

            call cav%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if

            areas(igrid) = cav%total_area
            volumes(igrid) = cav%total_volume
         end do

         ref_area = areas(n_grids)
         ref_volume = volumes(n_grids)

         write (*, '(a)') ''
         write (*, '(a,a,a,a,a,i0,a)') '  Molecule: ', trim(dataset_names(imol)), &
            '/', trim(mol_ids(imol)), ' (', mol%nat, ' atoms)'
         write (*, '(a8, 4a20)') &
            'nleb', 'Area (bohr^2)', 'Volume (bohr^3)', 'dA_ref (%)', 'dV_ref (%)'
         write (*, '(a8, 4a20)') &
            '-------', '-------------------', '-------------------', &
            '-------------------', '-------------------'

         do igrid = n_grids, 1, -1
            if (igrid == n_grids) then
               write (*, '(i8, 2f20.12, 2a20)') nleb_values(igrid), &
                  areas(igrid), volumes(igrid), '         ref        ', '         ref        '
            else
               write (*, '(i8, 4f20.12)') nleb_values(igrid), &
                  areas(igrid), volumes(igrid), &
                  100.0_wp * (areas(igrid) - ref_area) / ref_area, &
                  100.0_wp * (volumes(igrid) - ref_volume) / ref_volume
            end if
         end do

      end do

      write (*, '(a)') ''

      if (allocated(cav)) deallocate(cav)
      if (allocated(radii)) deallocate(radii)

   end subroutine test_convergence_iswig_lebedev

   !> Test convergence of DROP total-area and total-volume gradients w.r.t. Lebedev grid size.
   !> Computes the gradient at the finest level as reference, then reports the deviation
   !> of each coarser level from that reference (max abs, RMS, mean abs).
   subroutine test_convergence_drop_gradient(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error
      real(wp), allocatable :: ref_dA(:, :), ref_dV(:, :)
      real(wp), allocatable :: cur_dA(:, :), cur_dV(:, :)
      real(wp) :: max_a, rms_a, mad_a
      real(wp) :: max_v, rms_v, mad_v
      real(wp) :: diff
      integer :: igrid, iat, idir, n_comp

      integer, parameter :: n_grids = 29
      integer, parameter :: nleb_values(n_grids) = [ &
         6, 14, 26, 38, 50, 86, 110, 146, &
         170, 194, 302, 350, 434, 590, 770, 974, &
         1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, &
         3890, 4334, 4802, 5294, 5810]

      call get_structure(mol, "MB16-43", "CH4")
      n_comp = ndim * mol%nat

      allocate(ref_dA(ndim, mol%nat), ref_dV(ndim, mol%nat))
      allocate(cur_dA(ndim, mol%nat), cur_dV(ndim, mol%nat))

      ! Compute reference gradient at finest level
      allocate(cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
         call new_cavity_drop(cavity, nleb=nleb_values(n_grids), &
            tolerance=proj_tol, proj_maxiter=proj_maxiter, proj_level=proj_level, &
            do_fine=.true., &
            debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
            lsf_model=svdw_template, error=cavity_error)
      end block
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call cavity%get_gradient()

      do iat = 1, mol%nat
         do idir = 1, ndim
            ref_dA(idir, iat) = sum(cavity%asph1_rA(idir, :, iat))
            ref_dV(idir, iat) = sum(cavity%vsph1_rA(idir, :, iat))
         end do
      end do

      write (*, '(a)') ''
      write (*, '(a)') '================================================================================================'
      write (*, '(a)') 'Convergence: DROP gradient vs. Lebedev grid size (reference: finest level)'
      write (*, '(a,i0,a)') '  Molecule: MB16-43/CH4 (', mol%nat, ' atoms)'
      write (*, '(a)') '================================================================================================'
      write (*, '(a8, 3a20, 3a20)') &
         'nleb', 'max|dA_err|', 'rms(dA_err)', 'mad(dA_err)', &
         'max|dV_err|', 'rms(dV_err)', 'mad(dV_err)'
      write (*, '(a8, 6a20)') &
         '-------', '-------------------', '-------------------', '-------------------', &
         '-------------------', '-------------------', '-------------------'

      ! Print finest level as reference
      write (*, '(i8, 6a20)') nleb_values(n_grids), &
         '         ref        ', '         ref        ', '         ref        ', &
         '         ref        ', '         ref        ', '         ref        '

      ! Loop from second-finest to coarsest
      do igrid = n_grids - 1, 1, -1
         deallocate(cavity)
         allocate(cavity)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
            call new_cavity_drop(cavity, nleb=nleb_values(igrid), &
               tolerance=proj_tol, proj_maxiter=proj_maxiter, proj_level=proj_level, &
               do_fine=.true., &
               debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
               lsf_model=svdw_template, error=cavity_error)
         end block
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if

         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if
         call cavity%get_gradient()

         do iat = 1, mol%nat
            do idir = 1, ndim
               cur_dA(idir, iat) = sum(cavity%asph1_rA(idir, :, iat))
               cur_dV(idir, iat) = sum(cavity%vsph1_rA(idir, :, iat))
            end do
         end do

         max_a = 0.0_wp; rms_a = 0.0_wp; mad_a = 0.0_wp
         max_v = 0.0_wp; rms_v = 0.0_wp; mad_v = 0.0_wp

         do iat = 1, mol%nat
            do idir = 1, ndim
               diff = cur_dA(idir, iat) - ref_dA(idir, iat)
               max_a = max(max_a, abs(diff))
               rms_a = rms_a + diff**2
               mad_a = mad_a + abs(diff)

               diff = cur_dV(idir, iat) - ref_dV(idir, iat)
               max_v = max(max_v, abs(diff))
               rms_v = rms_v + diff**2
               mad_v = mad_v + abs(diff)
            end do
         end do

         rms_a = sqrt(rms_a / real(n_comp, wp))
         rms_v = sqrt(rms_v / real(n_comp, wp))
         mad_a = mad_a / real(n_comp, wp)
         mad_v = mad_v / real(n_comp, wp)

         write (*, '(i8, 6es20.8)') nleb_values(igrid), &
            max_a, rms_a, mad_a, max_v, rms_v, mad_v
      end do

      write (*, '(a)') ''

      deallocate(cavity)
      deallocate(ref_dA, ref_dV, cur_dA, cur_dV)

   end subroutine test_convergence_drop_gradient

   !> Test convergence of DROP cavity area and volume w.r.t. Lebedev grid size
   !> for various values of the blending parameter k.
   subroutine test_convergence_drop_nleb_blendk(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error
      integer :: igrid, ik

      integer, parameter :: n_grids = 29
      integer, parameter :: nleb_values(n_grids) = [ &
         6, 14, 26, 38, 50, 86, 110, 146, &
         170, 194, 302, 350, 434, 590, 770, 974, &
         1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, &
         3890, 4334, 4802, 5294, 5810]

      integer, parameter :: n_blendk = 16
      real(wp), parameter :: blendk_values(n_blendk) = [ &
         1.0_wp, 1.25_wp, 1.5_wp, 1.75_wp, 2.0_wp, 2.5_wp, 3.0_wp, 3.5_wp, &
         4.0_wp, 4.5_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp, 10.0_wp]

      call get_structure(mol, "UPU23", "4b")

      write (*, '(a)') ''
      write (*, '(a)') 'blend_k, nleb, area, volume'

      do ik = 1, n_blendk
         do igrid = 1, n_grids
            if (allocated(cavity)) deallocate(cavity)
            allocate(cavity)
            block
               type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
               call svdw_template%new(blend_k=blendk_values(ik), blend_2b=blend_2b, &
                  blend_3b=blend_3b)
               call new_cavity_drop(cavity, nleb=nleb_values(igrid), &
                  tolerance=proj_tol, proj_maxiter=proj_maxiter, proj_level=proj_level, &
                  debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
                  lsf_model=svdw_template, error=cavity_error)
            end block
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if

            call cavity%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if

            write (*, '(f8.2, a, i8, a, es24.15, a, es24.15)') &
               blendk_values(ik), ',', nleb_values(igrid), ',', &
               cavity%total_area, ',', cavity%total_volume
         end do
      end do

      if (allocated(cavity)) deallocate(cavity)

   end subroutine test_convergence_drop_nleb_blendk

   !> Test convergence of DROP cavity area and volume w.r.t. Lebedev grid size
   !> for various values of the projection tolerance.
   subroutine test_convergence_drop_nleb_projtol(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error
      integer :: igrid, itol

      integer, parameter :: n_grids = 29
      integer, parameter :: nleb_values(n_grids) = [ &
         6, 14, 26, 38, 50, 86, 110, 146, &
         170, 194, 302, 350, 434, 590, 770, 974, &
         1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, &
         3890, 4334, 4802, 5294, 5810]

      integer, parameter :: n_tol = 37
      real(wp), parameter :: tol_values(n_tol) = [ &
         1.0e-6_wp, 7.5e-7_wp, 5.0e-7_wp, 2.5e-7_wp, &
         1.0e-7_wp, 7.5e-8_wp, 5.0e-8_wp, 2.5e-8_wp, &
         1.0e-8_wp, 7.5e-9_wp, 5.0e-9_wp, 2.5e-9_wp, &
         1.0e-9_wp, 7.5e-10_wp, 5.0e-10_wp, 2.5e-10_wp, &
         1.0e-10_wp, 7.5e-11_wp, 5.0e-11_wp, 2.5e-11_wp, &
         1.0e-11_wp, 7.5e-12_wp, 5.0e-12_wp, 2.5e-12_wp, &
         1.0e-12_wp, 7.5e-13_wp, 5.0e-13_wp, 2.5e-13_wp, &
         1.0e-13_wp, 7.5e-14_wp, 5.0e-14_wp, 2.5e-14_wp, &
         1.0e-14_wp, 7.5e-15_wp, 5.0e-15_wp, 2.5e-15_wp, &
         1.0e-15_wp]

      call get_structure(mol, "UPU23", "4b")

      write (*, '(a)') ''
      write (*, '(a)') 'proj_tol, nleb, area, volume'

      do itol = 1, n_tol
         do igrid = 1, n_grids
            if (allocated(cavity)) deallocate(cavity)
            allocate(cavity)
            block
               type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
               call svdw_template%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
               call new_cavity_drop(cavity, nleb=nleb_values(igrid), &
                  tolerance=tol_values(itol), proj_maxiter=proj_maxiter, &
                  proj_level=proj_level, &
                  debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
                  lsf_model=svdw_template, error=cavity_error)
            end block
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if

            call cavity%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if

            write (*, '(es12.1, a, i8, a, es24.15, a, es24.15)') &
               tol_values(itol), ',', nleb_values(igrid), ',', &
               cavity%total_area, ',', cavity%total_volume
         end do
      end do

      if (allocated(cavity)) deallocate(cavity)

   end subroutine test_convergence_drop_nleb_projtol

end module test_cavity_drop_convergence
