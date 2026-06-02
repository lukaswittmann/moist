module test_cavity_drop_integration
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_radii, only : default_cpcm_radii
   use moist_data_radii_legacy, only: get_radius_func
   use mstore, only: get_structure
   use, intrinsic :: iso_fortran_env, only: error_unit

   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_marchingcubes, only: integrate_surface_marching_cubes
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_utils_env, only: resolve_dir, ensure_dir

   implicit none
   private

   public :: collect_cavity_drop_integration

   integer, parameter :: ndim = 3

   real(wp), parameter :: k = 2.5
   real(wp), parameter :: beta = 1.0_wp
   real(wp), parameter :: gamma = 1.0
   integer, parameter :: NUM_LEB = 194

   real(wp), parameter :: PROJ_TOL = 1E-12_wp
   integer, parameter :: PROJ_MAXITER = 1000
   integer, parameter :: PROJ_LEVEL = 1

   ! real, parameter :: k_array(7) = [1.5_wp, 2.0_wp, 2.5_wp, 3.0_wp, 5.0_wp, 10.0_wp, 20.0_wp]
   real, parameter :: k_array(7) = [1.75_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 10.0_wp, 15.0_wp]
   ! real, parameter :: gamma_array(4) = [0.0_wp, 0.5_wp, 1.0_wp, 1.5_wp]
   real, parameter :: gamma_array(2) = [0.0_wp, 1.0_wp]

contains
   subroutine collect_cavity_drop_integration(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      ! TODO: re-enable the per-dataset marching-cubes tests (upu23, amino20x4, mb16_43,
      ! but14diol, il16, dimer_pes) once the comparison harness is finalized.
      testsuite = [ &
                  new_unittest("mc_mixed", test_mc_mixed) &
                  ]
   end subroutine collect_cavity_drop_integration

   !> Test ray-casting on UPU23 molecules (sample)
   subroutine test_mc_upu23(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp) :: mc_area, mc_volume, cavity_area, cavity_volume
      real(wp) :: grid_bounds(3, 2)
      integer :: grid_steps(3)
      integer :: iid, k_idx, gamma_idx
      real(wp) :: k_val, gamma_val

      ! Sample UPU23 IDs to test
      character(len=2), parameter :: test_ids(4) = [ &
                                     '1c', '2h', '4b', '7p' &
                                     ]

      ! Print header
      call print_comparison_header("UPU23")

      do k_idx = 1, size(k_array)
         do gamma_idx = 1, size(gamma_array)
            k_val = k_array(k_idx)
            gamma_val = gamma_array(gamma_idx)

            do iid = 1, size(test_ids)
               call get_structure(mol, "UPU23", trim(test_ids(iid)))

               ! Auto-compute grid bounds and steps
               grid_bounds = 0.0_wp
               grid_steps = 0

               ! Compare (radii will be auto-computed from atomic numbers)
               call compare_mc_cavity_cached(k_val, beta, gamma_val, mol, "UPU23", trim(test_ids(iid)), &
                                             grid_bounds=grid_bounds, grid_steps=grid_steps, &
                                             mc_area=mc_area, mc_volume=mc_volume, &
                                             cavity_area=cavity_area, cavity_volume=cavity_volume, error=error)
               if (allocated(error)) return

               call print_comparison_row(trim(test_ids(iid)), k_val, beta, gamma_val, &
                                         mc_area, cavity_area, mc_volume, cavity_volume)

            end do
         end do
      end do

   end subroutine test_mc_upu23


   !> Test ray-casting on amino20x4 molecules (sample)
   subroutine test_mc_amino20x4(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp) :: mc_area, mc_volume, cavity_area, cavity_volume
      real(wp) :: grid_bounds(3, 2)
      integer :: grid_steps(3)
      integer :: iid, k_idx, gamma_idx
      real(wp) :: k_val, gamma_val

      ! Sample a few amino20x4 IDs to test
      character(len=7), parameter :: test_ids(4) = [ &
                                     'GLN_xai', 'PHE_xab', 'THR_xab', 'VAL_xad' &
                                     ]

      call print_comparison_header("amino20x4")

      do k_idx = 1, size(k_array)
         do gamma_idx = 1, size(gamma_array)
            k_val = k_array(k_idx)
            gamma_val = gamma_array(gamma_idx)

            do iid = 1, size(test_ids)
               call get_structure(mol, "Amino20x4", test_ids(iid))

               ! Auto-compute grid bounds and steps
               grid_bounds = 0.0_wp
               grid_steps = 0

               ! Compare (radii will be auto-computed from atomic numbers)
               call compare_mc_cavity_cached(k_val, beta, gamma_val, mol, "Amino20x4", test_ids(iid), &
                                             grid_bounds=grid_bounds, grid_steps=grid_steps, &
                                             mc_area=mc_area, mc_volume=mc_volume, &
                                             cavity_area=cavity_area, cavity_volume=cavity_volume, error=error)
               if (allocated(error)) return

               call print_comparison_row(trim(test_ids(iid)), k_val, beta, gamma_val, &
                                         mc_area, cavity_area, mc_volume, cavity_volume)

            end do
         end do
      end do

   end subroutine test_mc_amino20x4

   !> Test ray-casting on MB16-43 molecules (sample)
   subroutine test_mc_mb16_43(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp) :: mc_area, mc_volume, cavity_area, cavity_volume
      real(wp) :: grid_bounds(3, 2)
      integer :: grid_steps(3)
      integer :: iid, k_idx, gamma_idx
      real(wp) :: k_val, gamma_val

      ! Sample MB16-43 IDs to test
      character(len=4), parameter :: test_ids(4) = [ &
                                     'O2  ', 'CH4 ', '26  ', '40  ' &
                                     ]

      ! Print header
      call print_comparison_header("MB16-43")

      do k_idx = 1, size(k_array)
         do gamma_idx = 1, size(gamma_array)
            k_val = k_array(k_idx)
            gamma_val = gamma_array(gamma_idx)

            do iid = 1, size(test_ids)
               call get_structure(mol, "MB16-43", trim(test_ids(iid)))

               ! Auto-compute grid bounds and steps
               grid_bounds = 0.0_wp
               grid_steps = 0

               ! Compare (radii will be auto-computed from atomic numbers)
               call compare_mc_cavity_cached(k_val, beta, gamma_val, mol, "MB16-43", &
                                             trim(test_ids(iid)), &
                                             grid_bounds=grid_bounds, grid_steps=grid_steps, &
                                             mc_area=mc_area, mc_volume=mc_volume, &
                                             cavity_area=cavity_area, cavity_volume=cavity_volume, error=error)
               if (allocated(error)) return

               call print_comparison_row(trim(test_ids(iid)), k_val, beta, gamma_val, &
                                         mc_area, cavity_area, mc_volume, cavity_volume)

            end do
         end do
      end do

   end subroutine test_mc_mb16_43

   !> Test ray-casting on But14diol molecules (sample)
   subroutine test_mc_but14diol(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp) :: mc_area, mc_volume, cavity_area, cavity_volume
      real(wp) :: grid_bounds(3, 2)
      integer :: grid_steps(3)
      integer :: iid, k_idx, gamma_idx
      real(wp) :: k_val, gamma_val

      ! Sample But14diol IDs to test
      character(len=2), parameter :: test_ids(4) = [ &
                                     '20', '30', '40', '50']

      ! Print header
      call print_comparison_header("But14diol")

      do k_idx = 1, size(k_array)
         do gamma_idx = 1, size(gamma_array)
            k_val = k_array(k_idx)
            gamma_val = gamma_array(gamma_idx)

            do iid = 1, size(test_ids)
               call get_structure(mol, "But14diol", trim(test_ids(iid)))

               ! Auto-compute grid bounds and steps
               grid_bounds = 0.0_wp
               grid_steps = 0

               ! Compare (radii will be auto-computed from atomic numbers)
               call compare_mc_cavity_cached(k_val, beta, gamma_val, mol, "But14diol", &
                                             trim(test_ids(iid)), &
                                             grid_bounds=grid_bounds, grid_steps=grid_steps, &
                                             mc_area=mc_area, mc_volume=mc_volume, &
                                             cavity_area=cavity_area, cavity_volume=cavity_volume, error=error)
               if (allocated(error)) return

               call print_comparison_row(trim(test_ids(iid)), k_val, beta, gamma_val, &
                                         mc_area, cavity_area, mc_volume, cavity_volume)

            end do
         end do
      end do

   end subroutine test_mc_but14diol

   !> Test ray-casting on IL16 molecules (sample)
   subroutine test_mc_il16(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp) :: mc_area, mc_volume, cavity_area, cavity_volume
      real(wp) :: grid_bounds(3, 2)
      integer :: grid_steps(3)
      integer :: iid, k_idx, gamma_idx
      real(wp) :: k_val, gamma_val

      ! Sample IL16 IDs to test
      character(len=4), parameter :: test_ids(4) = [ &
                                     '008 ', '144 ', '187 ', '212 ']

      ! Print header
      call print_comparison_header("IL16")

      do k_idx = 1, size(k_array)
         do gamma_idx = 1, size(gamma_array)
            k_val = k_array(k_idx)
            gamma_val = gamma_array(gamma_idx)

            do iid = 1, size(test_ids)
               call get_structure(mol, "IL16", trim(test_ids(iid)))

               ! Auto-compute grid bounds and steps
               grid_bounds = 0.0_wp
               grid_steps = 0

               ! Compare (radii will be auto-computed from atomic numbers)
               call compare_mc_cavity_cached(k_val, beta, gamma_val, mol, "IL16", trim(test_ids(iid)), &
                                             grid_bounds=grid_bounds, grid_steps=grid_steps, &
                                             mc_area=mc_area, mc_volume=mc_volume, &
                                             cavity_area=cavity_area, cavity_volume=cavity_volume, error=error)
               if (allocated(error)) return

               call print_comparison_row(trim(test_ids(iid)), k_val, beta, gamma_val, &
                                         mc_area, cavity_area, mc_volume, cavity_volume)

            end do
         end do
      end do

   end subroutine test_mc_il16

   !> Test ray-casting on mixed molecules from different datasets
   subroutine test_mc_mixed(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp) :: mc_area, mc_volume, cavity_area, cavity_volume
      real(wp) :: grid_bounds(3, 2)
      integer :: grid_steps(3)
      integer :: iid, k_idx, gamma_idx
      real(wp) :: k_val, gamma_val

      integer :: n_samples
      real(wp) :: err_A, err_V, pct_err_A, pct_err_V
      real(wp) :: sum_err_A, sum_abs_err_A, sum_sq_err_A, max_abs_err_A
      real(wp) :: sum_err_V, sum_abs_err_V, sum_sq_err_V, max_abs_err_V
      real(wp) :: sum_pct_err_A, sum_abs_pct_err_A, max_abs_pct_err_A
      real(wp) :: sum_pct_err_V, sum_abs_pct_err_V, max_abs_pct_err_V
      real(wp) :: mean_err_A, mean_abs_err_A, std_err_A, mpe_A, mape_A
      real(wp) :: mean_err_V, mean_abs_err_V, std_err_V, mpe_V, mape_V

      ! Mixed molecules from different datasets
      integer, parameter :: num_mols = 8
      integer, parameter :: do_num_mols = 8
      character(len=12), parameter :: dataset_names(num_mols) = [ &
                                     'MB16-43     ', 'MB16-43     ', 'MB16-43     ', &
                                     'Amino20x4   ', 'MB16-43     ', &
                                     'But14diol   ', 'IL16        ', &
                                     'UPU23       ' &
                                     ]
      character(len=7), parameter :: mol_ids(num_mols) = [ &
                                     'Ar     ', 'O2     ', 'CH4    ', &
                                     'THR_xab', '16     ', &
                                     '30     ', '144    ', &
                                     '4b     ' &
                                     ]

                                     ! Print header
      call print_comparison_header("mixed datasets")

      n_samples = 0
      sum_err_A = 0.0_wp; sum_abs_err_A = 0.0_wp; sum_sq_err_A = 0.0_wp; max_abs_err_A = 0.0_wp
      sum_err_V = 0.0_wp; sum_abs_err_V = 0.0_wp; sum_sq_err_V = 0.0_wp; max_abs_err_V = 0.0_wp
      sum_pct_err_A = 0.0_wp; sum_abs_pct_err_A = 0.0_wp; max_abs_pct_err_A = 0.0_wp
      sum_pct_err_V = 0.0_wp; sum_abs_pct_err_V = 0.0_wp; max_abs_pct_err_V = 0.0_wp

      do iid = 1, do_num_mols
         do k_idx = 1, size(k_array)
            do gamma_idx = 1, size(gamma_array)
               k_val = k_array(k_idx)
               gamma_val = gamma_array(gamma_idx)


               if (trim(mol_ids(iid)) == 'Ar') then
                  call new(mol, [18], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
               else
                  call get_structure(mol, trim(dataset_names(iid)), trim(mol_ids(iid)))
               end if

               ! Auto-compute grid bounds and steps
               grid_bounds = 0.0_wp
               grid_steps = 0

               ! Compare (radii will be auto-computed from atomic numbers)
               call compare_mc_cavity_cached(k_val, beta, gamma_val, mol, &
                                             trim(dataset_names(iid)), trim(mol_ids(iid)), &
                                             grid_bounds=grid_bounds, grid_steps=grid_steps, &
                                             mc_area=mc_area, mc_volume=mc_volume, &
                                             cavity_area=cavity_area, cavity_volume=cavity_volume, error=error)
               if (allocated(error)) return

               call print_comparison_row(trim(mol_ids(iid)), k_val, beta, gamma_val, &
                                         mc_area, cavity_area, mc_volume, cavity_volume)

               err_A = cavity_area - mc_area
               err_V = cavity_volume - mc_volume
               if (mc_area > 1.0e-12_wp) then
                  pct_err_A = 100.0_wp * err_A / mc_area
               else
                  pct_err_A = 0.0_wp
               end if
               if (mc_volume > 1.0e-12_wp) then
                  pct_err_V = 100.0_wp * err_V / mc_volume
               else
                  pct_err_V = 0.0_wp
               end if

               n_samples = n_samples + 1
               sum_err_A = sum_err_A + err_A
               sum_err_V = sum_err_V + err_V
               sum_abs_err_A = sum_abs_err_A + abs(err_A)
               sum_abs_err_V = sum_abs_err_V + abs(err_V)
               sum_sq_err_A = sum_sq_err_A + err_A**2
               sum_sq_err_V = sum_sq_err_V + err_V**2
               max_abs_err_A = max(max_abs_err_A, abs(err_A))
               max_abs_err_V = max(max_abs_err_V, abs(err_V))

               sum_pct_err_A = sum_pct_err_A + pct_err_A
               sum_pct_err_V = sum_pct_err_V + pct_err_V
               sum_abs_pct_err_A = sum_abs_pct_err_A + abs(pct_err_A)
               sum_abs_pct_err_V = sum_abs_pct_err_V + abs(pct_err_V)
               max_abs_pct_err_A = max(max_abs_pct_err_A, abs(pct_err_A))
               max_abs_pct_err_V = max(max_abs_pct_err_V, abs(pct_err_V))

               ! ! Area
               ! call check(error, cavity_area, &
               !    mc_area, thr=5.0_wp, &
               !    more="Absolute difference in area too large")
               ! if (allocated(error)) return
               ! call check(error, cavity_area / mc_area, &
               !    1.0_wp, thr=0.008_wp, &
               !    more="Relative difference in area too large")
               ! if (allocated(error)) return

               ! ! Volume
               ! call check(error, cavity_volume, &
               !    mc_volume, thr=10.0_wp, &
               !    more="Absolute difference in volume too large")
               ! if (allocated(error)) return
               ! call check(error, cavity_volume / mc_volume, &
               !    1.0_wp, thr=0.007_wp, &
               !    more="Relative difference in volume too large")
               ! if (allocated(error)) return

            end do
         end do
      end do

      if (n_samples > 0) then
         mean_err_A = sum_err_A / n_samples
         mean_abs_err_A = sum_abs_err_A / n_samples
         std_err_A = sqrt(max(sum_sq_err_A / n_samples - mean_err_A**2, 0.0_wp))
         mpe_A = sum_pct_err_A / n_samples
         mape_A = sum_abs_pct_err_A / n_samples

         mean_err_V = sum_err_V / n_samples
         mean_abs_err_V = sum_abs_err_V / n_samples
         std_err_V = sqrt(max(sum_sq_err_V / n_samples - mean_err_V**2, 0.0_wp))
         mpe_V = sum_pct_err_V / n_samples
         mape_V = sum_abs_pct_err_V / n_samples

         write(*,*) ' '
         write(*,*) '--- Summary Statistics across whole set ---'
         write(*,*) 'Area:'
         write(*,'(A,F10.4)')   '  Mean Error:           ', mean_err_A
         write(*,'(A,F10.4)')   '  Mean Absolute Error:  ', mean_abs_err_A
         write(*,'(A,F10.4)')   '  Standard Deviation:   ', std_err_A
         write(*,'(A,F10.4)')   '  Max Abs Error:        ', max_abs_err_A
         write(*,'(A,F10.4,A)') '  Mean % Error (MPE):   ', mpe_A, ' %'
         write(*,'(A,F10.4,A)') '  Mean Abs % Error:     ', mape_A, ' %'
         write(*,'(A,F10.4,A)') '  Max Abs % Error:      ', max_abs_pct_err_A, ' %'
         write(*,*) 'Volume:'
         write(*,'(A,F10.4)')   '  Mean Error:           ', mean_err_V
         write(*,'(A,F10.4)')   '  Mean Absolute Error:  ', mean_abs_err_V
         write(*,'(A,F10.4)')   '  Standard Deviation:   ', std_err_V
         write(*,'(A,F10.4)')   '  Max Abs Error:        ', max_abs_err_V
         write(*,'(A,F10.4,A)') '  Mean % Error (MPE):   ', mpe_V, ' %'
         write(*,'(A,F10.4,A)') '  Mean Abs % Error:     ', mape_V, ' %'
         write(*,'(A,F10.4,A)') '  Max Abs % Error:      ', max_abs_pct_err_V, ' %'
         write(*,*) ' '
      end if

   end subroutine test_mc_mixed

   !> Test dimer PES scan
   subroutine test_mc_dimer_pes(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp) :: mc_area, mc_volume, cavity_area, cavity_volume
      real(wp) :: grid_bounds(3, 2)
      integer :: grid_steps(3)
      integer :: k_idx, gamma_idx, i_dist
      real(wp) :: k_val, gamma_val, distance
      real(wp), allocatable :: radii(:)
      character(len=12) :: dist_str
      real(wp), parameter :: step_size = 0.2_wp !bohr

      ! Allocate radii for 2 atoms
      allocate(radii(2))
      radii = 3.0_wp

      do k_idx = 1, size(k_array)
         do gamma_idx = 1, size(gamma_array)

            ! Print header
            call print_comparison_header("Dimer scan with k=" // &
                                         trim(to_string(k_array(k_idx))) // &
                                         ", gamma=" // trim(to_string(gamma_array(gamma_idx))))

            do i_dist = 0, 100
               distance = (step_size * i_dist + 0.01_wp) ** (1.0_wp / 3.0_wp)

               ! Create dimer with specified distance
               call new(mol, [8, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                  distance, distance, distance], [3, 2]))


               k_val = k_array(k_idx)
               gamma_val = gamma_array(gamma_idx)

               ! Auto-compute grid bounds and steps
               grid_bounds = 0.0_wp
               grid_steps = 0

               ! Compare with explicit radii
               call compare_mc_cavity_cached(k_val, beta, gamma_val, mol, "dimer_scan", &
                                            to_string(distance ** 3), radii, &
                                            grid_bounds, grid_steps, &
                                            mc_area, mc_volume, &
                                            cavity_area, cavity_volume, error=error)
               if (allocated(error)) return

               ! Format distance as string for display
               write(dist_str, '(F5.1)') distance ** 3
               call print_comparison_row(trim(adjustl(dist_str)), k_val, beta, gamma_val, &
                                        mc_area, cavity_area, mc_volume, cavity_volume)

               ! Exit if the are is that of two separated spheres
               if (mc_area >= 2.0_wp * 4.0_wp * pi * radii(1)**2) exit

            end do
         end do
      end do

      deallocate(radii)

   end subroutine test_mc_dimer_pes

   !* ================================================================================= *!
   !*                               Helper print routines                               *!
   !* ================================================================================= *!

   !> Print table header for MC vs DROP comparison
   subroutine print_comparison_header(title)
      character(len=*), intent(in) :: title

      write (*, *)
      write (*, '(A)') "  Comparing MC vs DROP for "//trim(title)//" molecules"
      write (*, '(1x, A12, 3A6, 2x, 3A12, 2A8, 2x, 3A12, 2A8)') "struct", "k", "b", "g", &
         "A_cav", "A_exact", "dA", "%", "x", "V_cav", "V_exact", "dV", "%", "x"
      write (*, '(1x, A12, 3A6, 2x, 3A12, 2A8, 2x, 3A12, 2A8)') "-----------", "-----", &
         "-----", "-----", "-----------", "-----------", "-----------", "-------", "-------", &
         "-----------", "-----------", "-----------", "-------", "-------"
   end subroutine print_comparison_header

   !> Print comparison data row
   subroutine print_comparison_row(struct_name, k_val, beta_val, gamma_val, &
                                   mc_area, cavity_area, mc_volume, cavity_volume)
      character(len=*), intent(in) :: struct_name
      real(wp), intent(in) :: k_val, beta_val, gamma_val
      real(wp), intent(in) :: mc_area, cavity_area, mc_volume, cavity_volume
      real(wp) :: ratio_area, ratio_volume

      ! Calculate ratios: cavity/exact, negative if cavity < exact
      ! Handle zero values to avoid division by zero
      if (abs(cavity_area) < 1.0e-14_wp .or. abs(mc_area) < 1.0e-14_wp) then
         ratio_area = 0.0_wp
      else if (cavity_area >= mc_area) then
         ratio_area = cavity_area / mc_area
      else
         ratio_area = -(mc_area / cavity_area)
      end if

      if (abs(cavity_volume) < 1.0e-14_wp .or. abs(mc_volume) < 1.0e-14_wp) then
         ratio_volume = 0.0_wp
      else if (cavity_volume >= mc_volume) then
         ratio_volume = cavity_volume / mc_volume
      else
         ratio_volume = -(mc_volume / cavity_volume)
      end if

      write (*, '(1x, A12, 3F6.1, 2x, 3F12.1, F8.1, F8.2, 2x, 3F12.1, F8.1, F8.2)') &
         trim(struct_name), k_val, beta_val, gamma_val, &
         cavity_area, mc_area, cavity_area - mc_area, &
         100.0_wp*(cavity_area - mc_area)/cavity_area, ratio_area, &
         cavity_volume, mc_volume, cavity_volume - mc_volume, &
         100.0_wp*(cavity_volume - mc_volume)/cavity_volume, ratio_volume
   end subroutine print_comparison_row

   !* ================================================================================= *!
   !*                                  Working routines                                 *!
   !* ================================================================================= *!

   !> Compare raycast integration vs cavity implementation for a given molecule
   subroutine compare_mc_cavity(k, beta, gamma, mol, radii, grid_bounds, grid_steps, &
                                mc_area, mc_volume, cavity_area, cavity_volume, &
                                target_spacing, error)
      real(wp), intent(in) :: k, beta, gamma
      type(structure_type), intent(in) :: mol
      real(wp), intent(in), optional :: radii(:)
      real(wp), intent(inout) :: grid_bounds(3, 2)
      integer, intent(inout) :: grid_steps(3)
      real(wp), intent(out) :: mc_area, mc_volume
      real(wp), intent(out) :: cavity_area, cavity_volume
      real(wp), intent(in), optional :: target_spacing  ! Target grid spacing (default: 0.3 Bohr)
      !> Error handle: set (and caller should return) if the cavity build fails
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_parameters_type) :: param
      real(wp), allocatable :: radii_local(:)
      integer :: i
      type(mctc_error), allocatable :: cavity_error

      ! Get radii (either from argument or compute from atomic numbers)
      if (present(radii)) then
         allocate (radii_local(size(radii)))
         radii_local = radii
      else
         allocate (radii_local(mol%nat))
         do i = 1, mol%nat
            radii_local(i) = get_radius_func(mol%num(mol%id(i)))
         end do
      end if

      ! Initialize cavity with Lebedev grid
      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=k, blend_2b=beta, blend_3b=gamma)
         call new_cavity_drop(cavity, nleb=NUM_LEB, &
                             tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL, &
                             debug=.false., verbose=0, &
                             radius_model=default_cpcm_radii(), &
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

      ! Set up LSF primitive for raycast integration
      call lsf%new(blend_k=k, blend_2b=beta, blend_3b=gamma)
      !> Direct LSF use (no cavity to set screening); we own this.
      lsf%screening_threshold = PROJ_TOL * 0.1_wp
      call lsf%update(mol, radii_local)

      ! Integrate using marching cubes
      call integrate_surface_marching_cubes(lsf, mol%xyz, &
                                            mc_area, mc_volume, debug=.false., &
                                            target_spacing=target_spacing)

      ! Get cavity results
      cavity_area = cavity%total_area
      cavity_volume = cavity%total_volume

      deallocate (radii_local)

   end subroutine compare_mc_cavity

   !> Compare with caching - reads from cache if available, otherwise computes and caches
   subroutine compare_mc_cavity_cached(k, beta, gamma, mol, benchmark_name, structure_id, &
                                       radii, grid_bounds, grid_steps, &
                                       mc_area, mc_volume, cavity_area, cavity_volume, &
                                       target_spacing, error)
      real(wp), intent(in) :: k, beta, gamma
      type(structure_type), intent(in) :: mol
      character(len=*), intent(in) :: benchmark_name, structure_id
      real(wp), intent(in), optional :: radii(:)
      real(wp), intent(inout) :: grid_bounds(3, 2)
      integer, intent(inout) :: grid_steps(3)
      real(wp), intent(out) :: mc_area, mc_volume
      real(wp), intent(out) :: cavity_area, cavity_volume
      real(wp), intent(in), optional :: target_spacing
      !> Error handle: set (and caller should return) if the cavity build fails
      type(error_type), allocatable, intent(out) :: error

      character(len=512) :: cache_file
      character(len=:), allocatable :: cache_dir
      character(len=256) :: sanitized_benchmark, sanitized_structure
      character(len=32) :: k_str, beta_str, gamma_str
      integer :: io_unit, io_stat, i
      logical :: file_exists, need_compute
      type(cavity_type_drop), allocatable :: cavity
      real(wp), allocatable :: radii_local(:)
      type(mctc_error), allocatable :: cavity_error

      ! Sanitize names (replace spaces and special chars with underscores)
      sanitized_benchmark = trim(benchmark_name)
      sanitized_structure = trim(structure_id)
      do i = 1, len_trim(sanitized_benchmark)
         if (sanitized_benchmark(i:i) == ' ' .or. sanitized_benchmark(i:i) == '-') then
            sanitized_benchmark(i:i) = '_'
         end if
      end do
      do i = 1, len_trim(sanitized_structure)
         if (sanitized_structure(i:i) == ' ' .or. sanitized_structure(i:i) == '-') then
            sanitized_structure(i:i) = '_'
         end if
      end do

      ! Format k, beta and gamma as strings
      write (k_str, '(F0.1)') k
      write (beta_str, '(F0.1)') beta
      write (gamma_str, '(F0.1)') gamma

      ! Resolve cache directory (MOIST_DROP_CACHE_DIR override, repo-relative default)
      cache_dir = resolve_dir("MOIST_DROP_CACHE_DIR", "test/outputs/drop_cache")
      call ensure_dir(cache_dir)

      ! Construct cache file path
      cache_file = trim(cache_dir)//"/benchmark_"//trim(sanitized_benchmark)//"_"// &
                   trim(sanitized_structure)//"_"//trim(adjustl(k_str))//"_"// &
                   trim(adjustl(beta_str))//"_"// &
                   trim(adjustl(gamma_str))//".txt"

      ! Check if cache file exists
      inquire (file=trim(cache_file), exist=file_exists)
      need_compute = .true.

      if (file_exists) then
         ! Read from cache
         open (newunit=io_unit, file=trim(cache_file), status='old', action='read', iostat=io_stat)
         if (io_stat == 0) then
            read (io_unit, *, iostat=io_stat) mc_area
            if (io_stat == 0) then
               read (io_unit, *, iostat=io_stat) mc_volume
               if (io_stat == 0) then
                  need_compute = .false.
               end if
            end if
            close (io_unit)

            if (io_stat /= 0) then
               write (error_unit, '(A)') "Warning: Error reading cache file, recomputing..."
            end if
         else
            write (error_unit, '(A)') "Warning: Could not open cache file, recomputing..."
         end if
      end if

      ! If cache doesn't exist or read failed, compute using compare_mc_cavity
      if (need_compute) then
         call compare_mc_cavity(k, beta, gamma, mol, radii, grid_bounds, grid_steps, &
                                mc_area, mc_volume, cavity_area, cavity_volume, 0.01_wp, error)
         if (allocated(error)) return

         ! Write mc results to cache (cavity results are not cached as they're fast)
         open (newunit=io_unit, file=trim(cache_file), status='replace', &
               action='write', iostat=io_stat)
         if (io_stat == 0) then
            write (io_unit, '(F20.10)') mc_area
            write (io_unit, '(F20.10)') mc_volume
            close (io_unit)
         else
            write (error_unit, '(A)') "Warning: Could not write cache file"
         end if
      else
         ! Cache hit! Only need to compute cavity results (mc results were read from cache)

         ! Get radii (either from argument or compute from atomic numbers)
         if (present(radii)) then
            allocate (radii_local(size(radii)))
            radii_local = radii
         else
            allocate (radii_local(mol%nat))
            do i = 1, mol%nat
               radii_local(i) = get_radius_func(mol%num(mol%id(i)))
            end do
         end if

         ! Initialize and compute cavity
         allocate (cavity)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=k, blend_2b=beta, blend_3b=gamma)
            call new_cavity_drop(cavity, nleb=NUM_LEB, &
                                tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                                proj_level=PROJ_LEVEL, &
                                debug=.false., verbose=0, &
                                radius_model=default_cpcm_radii(), &
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

         cavity_area = cavity%total_area
         cavity_volume = cavity%total_volume

         deallocate (radii_local)
      end if

   end subroutine compare_mc_cavity_cached

end module test_cavity_drop_integration
