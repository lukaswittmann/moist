module test_cavity_drop_timings
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use mctc_io_convert, only: aatoau
   use mstore, only: get_structure
   use moist_data_radii_legacy, only: get_radius_func
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_math_cell_grid, only: moist_cell_grid_type
   use moist_radii, only: default_cpcm_radii
   use moist_radii_static, only: static_radius_type
   use moist_cavity_drop_marchingcubes, only: integrate_surface_marching_cubes
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_cavity_drop_timings

   integer, parameter :: ndim = 3

contains

   subroutine collect_cavity_drop_timings(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("timing_drop_proj_levels", test_timing_drop_proj_levels), &
                  new_unittest("timing_cell_fraction_benchmark", test_timing_cell_fraction_benchmark), &
                  new_unittest("timing_drop_scaling", test_timing_drop_scaling), &
                  new_unittest("timing_mc_scaling", test_timing_mc_scaling), &
                  new_unittest("timing_iswig_scaling", test_timing_iswig_scaling) &
                  ]
   end subroutine collect_cavity_drop_timings

   subroutine test_timing_drop_proj_levels(error)
      use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity_drop
      real(wp), allocatable :: radii(:)
      character(len=20) :: struct_names(25)
      integer :: iat, iter, istruct, ilevel
      real :: t0, t1
      real(wp) :: time_update, time_gradient, time_total
      integer :: n_iter, nleb, proj_level
      real(wp) :: blend_k, blend_3b
      type(mctc_error), allocatable :: cavity_error
      integer, parameter :: n_proj_levels = 1
      integer, parameter :: proj_levels(n_proj_levels) = [1]

      nleb = 110
      blend_k = 3.0_wp
      blend_3b = 1.0_wp
      n_iter = 1

      struct_names = [character(len=20) :: 'polyala_04', 'polyala_08', 'polyala_12', &
                      'polyala_16', 'polyala_20', 'polyala_24', 'polyala_28', 'polyala_32', 'polyala_36', &
                      'polyala_40', 'polyala_44', 'polyala_48', 'polyala_52', 'polyala_56', 'polyala_60', &
                      'polyala_64', 'polyala_68', 'polyala_72', 'polyala_76', 'polyala_80', 'polyala_84', &
                      'polyala_88', 'polyala_92', 'polyala_96', 'polyala_100']

      write (*, '(a)') ''
      write (*, '(a)') '================================================================'
      write (*, '(a)') 'Benchmark: DROP Cavity - Projection Level Comparison'
      write (*, '(a, i0)') 'Parameters: nleb = ', nleb
      write (*, '(a, f6.3)') '            blend_k = ', blend_k
      write (*, '(a, f6.3)') '            blend_3b = ', blend_3b
      write (*, '(a, i0)') '            iterations = ', n_iter
      write (*, '(a)') '================================================================'
      write (*, '(a)') ''
      write (*, '(a14, a9, a10, a12, 3a13)') 'Structure', 'N_atoms', 'N_grid', &
         'proj_level', 'Update(s)', 'Gradient(s)', 'Total(s)'
      write (*, '(a14, a9, a10, a12, 3a13)') '-------------', '--------', '--------', &
         '-----------', '------------', '------------', '------------'

      do istruct = 1, size(struct_names)
         call get_structure(mol, 'POLYALANINE', trim(struct_names(istruct)))

         if (allocated(radii)) deallocate (radii)
         allocate (radii(mol%nat))
         do iat = 1, mol%nat
            radii(iat) = get_radius_func(mol%num(mol%id(iat)))*aatoau
         end do

         do ilevel = 1, n_proj_levels
            proj_level = proj_levels(ilevel)

            time_update = 0.0_wp
            time_gradient = 0.0_wp

            do iter = 1, n_iter
               if (allocated(cavity_drop)) deallocate (cavity_drop)
               allocate (cavity_drop)

               call cpu_time(t0)
               block
                  type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
                  call svdw_template%new(blend_k=blend_k, blend_3b=blend_3b)
                  call new_cavity_drop(cavity_drop, nleb=nleb, &
                                      verbose=0, debug=.false., &
                                      tolerance=1.0E-10_wp, proj_level=proj_level, &
                                      radius_model=default_cpcm_radii(), &
                                      lsf_model=svdw_template, error=cavity_error)
               end block
               if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
               call cavity_drop%update(mol, error=cavity_error)
               if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
               call cpu_time(t1)
               time_update = time_update + real(t1 - t0, wp)

               call cpu_time(t0)
               call cavity_drop%get_gradient()
               call cpu_time(t1)
               time_gradient = time_gradient + real(t1 - t0, wp)
            end do

            time_update = time_update/real(n_iter, wp)
            time_gradient = time_gradient/real(n_iter, wp)
            time_total = time_update + time_gradient

            write (*, '(a14, i9, i10, i12, 3f13.6)') &
               trim(struct_names(istruct)), mol%nat, cavity_drop%ngrid, proj_level, &
               time_update, time_gradient, time_total
         end do
      end do

      write (*, '(a)') ''

      if (allocated(radii)) deallocate (radii)

   end subroutine test_timing_drop_proj_levels

   !> Generate categorized test points: surface, interior, and exterior.
   !>
   !> Surface points are placed near atomic sphere surfaces where screening is
   !> least effective (many active atoms). Interior points sit deep inside where
   !> screening rapidly culls. Exterior points lie well outside the molecule.
   !>
   !> @param[in]  mol           Molecular structure
   !> @param[in]  radii         Atomic radii [n_atoms]
   !> @param[out] points        Output points [ndim, n_total]
   !> @param[in]  n_surface     Number of surface-proximate points
   !> @param[in]  n_interior    Number of interior points
   !> @param[in]  n_exterior    Number of exterior points
   subroutine generate_categorized_points(mol, radii, points, n_surface, n_interior, n_exterior)
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)
      real(wp), intent(out) :: points(:, :)
      integer, intent(in) :: n_surface, n_interior, n_exterior

      real(wp) :: rmin(ndim), rmax(ndim), random_vec(ndim), dir(ndim), dnorm
      real(wp) :: offset, r_max
      integer :: ipt, iat, idim, idx

      if (n_surface + n_interior + n_exterior > size(points, 2)) &
         error stop "generate_categorized_points: output array too small"

      points = 0.0_wp

      ! Compute bounding box
      rmin = mol%xyz(:, 1) - radii(1)
      rmax = mol%xyz(:, 1) + radii(1)
      do iat = 2, mol%nat
         do idim = 1, ndim
            rmin(idim) = min(rmin(idim), mol%xyz(idim, iat) - radii(iat))
            rmax(idim) = max(rmax(idim), mol%xyz(idim, iat) + radii(iat))
         end do
      end do
      r_max = maxval(radii)

      ! --- Surface points: atom center + radius * (1 + small_offset) * random_direction ---
      idx = 0
      do ipt = 1, n_surface
         ! Pick a random atom
         call random_number(offset)
         iat = 1 + int(offset*real(mol%nat, wp))
         iat = min(iat, mol%nat)
         ! Random unit direction
         call random_number(dir)
         dir = dir - 0.5_wp
         dnorm = norm2(dir)
         if (dnorm > 0.0_wp) dir = dir/dnorm
         ! Place at radius * (1 + small offset in [-0.1, 0.15])
         call random_number(offset)
         offset = -0.1_wp + offset*0.25_wp
         idx = idx + 1
         points(:, idx) = mol%xyz(:, iat) + radii(iat)*(1.0_wp + offset)*dir
      end do

      ! --- Interior points: atom center + radius * 0.3..0.6 * random_direction ---
      do ipt = 1, n_interior
         call random_number(offset)
         iat = 1 + int(offset*real(mol%nat, wp))
         iat = min(iat, mol%nat)
         call random_number(dir)
         dir = dir - 0.5_wp
         dnorm = norm2(dir)
         if (dnorm > 0.0_wp) dir = dir/dnorm
         call random_number(offset)
         offset = 0.3_wp + offset*0.3_wp
         idx = idx + 1
         points(:, idx) = mol%xyz(:, iat) + radii(iat)*offset*dir
      end do

      ! --- Exterior points: bounding box expanded by 2 * r_max ---
      do ipt = 1, n_exterior
         call random_number(random_vec)
         idx = idx + 1
         points(:, idx) = (rmin - 2.0_wp*r_max) &
                          + random_vec*((rmax + 2.0_wp*r_max) - (rmin - 2.0_wp*r_max))
      end do
   end subroutine generate_categorized_points

   !> Benchmark cell_fraction impact on cell-grid-screened SSD + LSF evaluation.
   !>
   !> For each structure, compares:
   !>   - Full SSD compute (no cell grid, brute-force baseline)
   !>   - Cell-grid-screened SSD + LSF with cell_fraction = 1.0, 0.5, 0.25
   !>
   !> Reports per-call timings, average candidates per query, and speedup
   !> relative to the unscreened baseline.
   subroutine test_timing_cell_fraction_benchmark(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf_prim
      type(moist_cell_grid_type) :: cell_grid
      real(wp), allocatable :: radii(:), r_eff(:)
      real(wp) :: point(ndim)
      real(wp) :: lsf_val, lsf_grad(ndim), lsf_hess(ndim, ndim)

      integer, parameter :: n_structs = 6
      character(len=20), parameter :: struct_sets(n_structs) = &
                                      [character(len=20) :: 'AMYLOSE', 'AMYLOSE', 'AMYLOSE', &
                                                             'POLYALANINE', 'POLYALANINE', 'POLYALANINE']
      character(len=20), parameter :: struct_names(n_structs) = &
                                      [character(len=20) :: 'Amylose16', 'Amylose32', 'Amylose64', &
                                                             'polyala_28', 'polyala_52', 'polyala_100']

      integer, parameter :: n_fractions = 6
      real(wp), parameter :: fractions(n_fractions) = [1.0_wp, 0.5_wp, 0.25_wp, 0.125_wp, 0.0625_wp, 0.03125_wp]

      integer, parameter :: n_surface = 1500
      integer, parameter :: n_interior = 500
      integer, parameter :: n_pts = n_surface + n_interior
      integer, parameter :: n_iter = 50
      real(wp), parameter :: blend_k = 4.0_wp
      real(wp), parameter :: blend_3b = 1.0_wp
      real(wp), parameter :: threshold = 1.0e-12_wp

      real(wp) :: all_points(ndim, n_pts)
      real :: t0, t1
      real(wp) :: t_full, t_screened, t_build, cand_sum
      integer :: istruct, ipt, iter, iat, ifrac, n_calls
      integer :: start, n_cand
      real(wp) :: delta
      character(len=8) :: frac_str

      ! Screening shell: exp(-k/3 * delta) = 0.1 * threshold
      delta = -3.0_wp/blend_k*log(0.1_wp*threshold)

      write (*, '(a)') ''
      write (*, '(a)') '========================================================================================='
      write (*, '(a)') 'Benchmark: cell_fraction impact on cell-grid-screened SSD + LSF'
      write (*, '(a, i0, a, i0)') 'Points: surface=', n_surface, '  interior=', n_interior
      write (*, '(a, i0, a, f5.2, a, es8.1)') 'Parameters: iterations=', n_iter, &
         '  k=', blend_k, '  threshold=', threshold
      write (*, '(a, f8.2, a)') 'Screening shell delta=', delta, ' bohr'
      write (*, '(a)') '========================================================================================='
      write (*, '(a)') ''
      write (*, '(a14, a9, a10, a12, a14, a14, a12, a9)') &
         'Structure', 'N_atoms', 'fraction', 'N_cells', 't_build(us)', 't_eval(us)', 'N_cand', 'speedup'
      write (*, '(a14, a9, a10, a12, a14, a14, a12, a9)') &
         '-------------', '--------', '---------', '-----------', &
         '-------------', '-------------', '-----------', '--------'

      do istruct = 1, n_structs
         call get_structure(mol, trim(struct_sets(istruct)), trim(struct_names(istruct)))

         if (allocated(radii)) deallocate (radii)
         if (allocated(r_eff)) deallocate (r_eff)
         allocate (radii(mol%nat), r_eff(mol%nat))
         do iat = 1, mol%nat
            radii(iat) = get_radius_func(mol%num(mol%id(iat)))*aatoau
            r_eff(iat) = radii(iat) + delta
         end do

         ! Setup LSF primitive (owns its internal SSD system)
         call lsf_prim%new(blend_k=blend_k, blend_1b=1.0_wp, blend_2b=1.0_wp, &
                           blend_3b=blend_3b)
         lsf_prim%screening_threshold = threshold
         call lsf_prim%set_max_deriv(2)
         call lsf_prim%update(mol, radii)

         ! Generate test points
         call generate_categorized_points(mol, radii, all_points, n_surface, n_interior, 0)

         ! --- Baseline: full SSD compute (no cell grid) ---
         t_full = 0.0_wp
         n_calls = 0
         do iter = 1, n_iter
            do ipt = 1, n_pts
               point = all_points(:, ipt)
               call cpu_time(t0)
               call lsf_prim%prepare(point)
               call lsf_prim%f012_r_screened(lsf_val, lsf_grad, lsf_hess)
               call cpu_time(t1)
               t_full = t_full + real(t1 - t0, wp)
               n_calls = n_calls + 1
            end do
         end do
         t_full = t_full/real(n_calls, wp)

         write (*, '(a14, i9, a10, a12, a14, f14.4, a12, a9)') &
            trim(struct_names(istruct)), mol%nat, 'full', &
            '-', '-', 1.0e6_wp*t_full, '-', '1.00x'

         ! --- Cell-grid-screened path for each fraction ---
         do ifrac = 1, n_fractions
            ! Build cell grid
            call cpu_time(t0)
            call cell_grid%build(mol%xyz, r_eff, cell_fraction=fractions(ifrac))
            call cpu_time(t1)
            t_build = real(t1 - t0, wp)

            ! Benchmark screened evaluation
            t_screened = 0.0_wp
            cand_sum = 0.0_wp
            n_calls = 0
            do iter = 1, n_iter
               do ipt = 1, n_pts
                  point = all_points(:, ipt)
                  call cpu_time(t0)
                  call cell_grid%query(point, start, n_cand)
                  call lsf_prim%prepare_subset(point, &
                                               cell_grid%cell_nlat(start + 1:start + n_cand))
                  call lsf_prim%f012_r_screened(lsf_val, lsf_grad, lsf_hess)
                  call cpu_time(t1)
                  t_screened = t_screened + real(t1 - t0, wp)
                  cand_sum = cand_sum + real(n_cand, wp)
                  n_calls = n_calls + 1
               end do
            end do
            t_screened = t_screened/real(n_calls, wp)

            write (frac_str, '(f5.2)') fractions(ifrac)
            write (*, '(a14, i9, a10, i12, f14.4, f14.4, f12.1, f7.2, a)') &
               '', mol%nat, adjustr(frac_str), cell_grid%ncells, &
               1.0e6_wp*t_build, &
               1.0e6_wp*t_screened, &
               cand_sum/real(n_calls, wp), &
               t_full/max(t_screened, 1.0e-30_wp), 'x'

            call cell_grid%destroy()
         end do

         write (*, '(a)') ''
      end do

      if (allocated(radii)) deallocate (radii)
      if (allocated(r_eff)) deallocate (r_eff)
   end subroutine test_timing_cell_fraction_benchmark

   !> Benchmark DROP cavity component scaling with system size.
   !>
   !> Runs full DROP cavity update and gradient computation for polyalanine
   !> chains of increasing size. Reads per-component wall times from the
   !> internal timer and fits t(N) = A * N^X via grid search to
   !> determine the formal scaling exponent of each component.
   !>
   !> @param[out] error  Test error (set on failure)
   subroutine test_timing_drop_scaling(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error

      !> Number of polyalanine structures to benchmark
      integer, parameter :: n_struct = 25
      !> Number of timer slots registered in new_cavity_drop
      integer, parameter :: n_timers = 33
      !> Minimum time (s) for a data point to be included in the fit
      real(wp), parameter :: t_min = 1.0e-6_wp

      !> Cavity parameters
      integer, parameter :: nleb = 194
      real(wp), parameter :: blend_k = 4.5_wp
      real(wp), parameter :: blend_3b = 1.0_wp
      integer, parameter :: proj_level = 3

      character(len=20) :: struct_names(n_struct)
      character(len=30) :: timer_labels(n_timers)
      logical :: is_parent(n_timers)

      !> Collected benchmark data
      integer :: n_atoms_arr(n_struct)
      integer :: n_grid_arr(n_struct)
      real(wp) :: times(n_timers, n_struct)
      real(wp) :: total_times(n_struct)

      !> Number of repetitions per structure for stable timings
      integer, parameter :: n_iter = 10

      !> Fit workspace
      real(wp) :: real_n(n_struct), raw_t(n_struct)
      real(wp) :: exponent, r_sq, coeff_a
      integer :: n_valid

      integer :: istruct, itimer, iter, iprint
      real(wp) :: rn_iter

      !> Print order (groups children under their parent)
      integer :: print_order(n_timers), n_print

      !> CSV output
      integer :: csv_unit
      character(len=*), parameter :: csv_file = 'drop_timings.csv'

      !> Polyalanine structures (increasing size)
      struct_names = [character(len=20) :: &
                      'polyala_04', 'polyala_08', 'polyala_12', 'polyala_16', 'polyala_20', &
                      'polyala_24', 'polyala_28', 'polyala_32', 'polyala_36', 'polyala_40', &
                      'polyala_44', 'polyala_48', 'polyala_52', 'polyala_56', 'polyala_60', &
                      'polyala_64', 'polyala_68', 'polyala_72', 'polyala_76', 'polyala_80', &
                      'polyala_84', 'polyala_88', 'polyala_92', 'polyala_96', 'polyala_100']

      !> Timer labels matching IDs registered in new_cavity_drop
      timer_labels = ' '
      timer_labels(1) = 'Setup'
      timer_labels(2) = '  Lebedev cache'
      timer_labels(3) = '  Array setup'
      timer_labels(4) = '  Adj. lists'
      timer_labels(5) = '  Switching func.'
      timer_labels(6) = '  Pre-filter'
      timer_labels(7) = 'Projector'
      timer_labels(8) = 'Post processing'
      timer_labels(9) = '  Filter'
      timer_labels(10) = '  Grid adj. list'
      timer_labels(11) = '  CP Jacobian'
      timer_labels(12) = '  Disconnected cav.'
      timer_labels(13) = 'Properties'
      timer_labels(14) = '  Grid density'
      timer_labels(15) = '  Curvatures'
      timer_labels(16) = '  Area & Volume'
      timer_labels(17) = '  Gaussians'
      timer_labels(18) = '  CPCM energy'
      timer_labels(19) = 'Gradients'
      timer_labels(20) = '  Primitives'
      timer_labels(21) = '  Positions'
      timer_labels(22) = '  Displacement'
      timer_labels(23) = '  Distances'
      timer_labels(25) = '  CP Jacobian (grad)'
      timer_labels(26) = '  Gaussian widths'
      timer_labels(27) = '  Switching (grad)'
      timer_labels(28) = '  Area (grad)'
      timer_labels(29) = '  Volume (grad)'
      timer_labels(30) = '  Surface normal'
      timer_labels(31) = '  Branch weights'
      timer_labels(32) = '  Branch weights (pp)'
      timer_labels(33) = '  CPCM (grad)'

      is_parent = .false.
      is_parent([1, 7, 8, 13, 19]) = .true.

      !> Print order: group children under their parent
      n_print = 0
      call add_order(1); call add_order(2); call add_order(3)
      call add_order(4); call add_order(5); call add_order(6)
      call add_order(7)
      call add_order(8); call add_order(9); call add_order(10)
      call add_order(11); call add_order(12); call add_order(32)
      call add_order(13); call add_order(14); call add_order(15)
      call add_order(16); call add_order(17); call add_order(18)
      call add_order(19); call add_order(20); call add_order(21)
      call add_order(22); call add_order(23); call add_order(30)
      call add_order(25); call add_order(26); call add_order(27)
      call add_order(28); call add_order(29); call add_order(31)
      call add_order(33)

      !> ====================== Benchmark loop ======================
      rn_iter = real(n_iter, wp)

      write (*, '(a)') ''
      write (*, '(a)') '=================================================================='
      write (*, '(a)') 'Benchmark: DROP cavity component scaling (update + gradient)'
      write (*, '(a,i0,a,f5.2,a,f5.2,a,i0,a,i0)') &
         'nleb=', nleb, '  k=', blend_k, '  b3=', blend_3b, &
         '  proj=', proj_level, '  iter=', n_iter
      write (*, '(a)') '=================================================================='

      !> Open CSV file and write header
      open (newunit=csv_unit, file=csv_file, status='replace', action='write')
      write (csv_unit, '(a)', advance='no') 'structure,n_atoms,n_grid,iter,total'
      do itimer = 1, n_timers
         if (len_trim(timer_labels(itimer)) > 0) then
            write (csv_unit, '(a,a)', advance='no') ',', trim(adjustl(timer_labels(itimer)))
         else
            write (csv_unit, '(a,i0)', advance='no') ',timer_', itimer
         end if
      end do
      write (csv_unit, *)

      do istruct = 1, n_struct
         call get_structure(mol, 'POLYALANINE', trim(struct_names(istruct)))
         n_atoms_arr(istruct) = mol%nat

         !> Zero accumulators
         times(:, istruct) = 0.0_wp
         total_times(istruct) = 0.0_wp

         do iter = 1, n_iter
            if (allocated(cavity)) deallocate (cavity)
            allocate (cavity)
            block
               type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
               call svdw_template%new(blend_k=blend_k, blend_3b=blend_3b)
               call new_cavity_drop(cavity, nleb=nleb, &
                                   verbose=0, debug=.false., do_fine=.true., &
                                   tolerance=1.0E-10_wp, proj_level=proj_level, &
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

            call cavity%get_gradient()

            !> Write per-iteration row to CSV
            write (csv_unit, '(a,a,i0,a,i0,a,i0,a,es14.6)', advance='no') &
               trim(struct_names(istruct)), ',', mol%nat, ',', cavity%ngrid, &
               ',', iter, ',', cavity%timer%get()
            do itimer = 1, n_timers
               write (csv_unit, '(a,es14.6)', advance='no') ',', cavity%timer%get(itimer)
            end do
            write (csv_unit, *)

            !> Accumulate timer values
            do itimer = 1, n_timers
               times(itimer, istruct) = times(itimer, istruct) &
                                        + cavity%timer%get(itimer)
            end do
            total_times(istruct) = total_times(istruct) + cavity%timer%get()
            n_grid_arr(istruct) = cavity%ngrid
         end do

         !> Average over iterations
         times(:, istruct) = times(:, istruct)/rn_iter
         total_times(istruct) = total_times(istruct)/rn_iter

         write (*, '(2x,a14,a,i5,a,i7,a,f10.3,a)') &
            trim(struct_names(istruct)), &
            '  N_at=', n_atoms_arr(istruct), &
            '  N_grid=', n_grid_arr(istruct), &
            '  avg=', total_times(istruct), ' s'
      end do

      if (allocated(cavity)) deallocate (cavity)
      close (csv_unit)
      write (*, '(a,a,a)') 'Per-iteration timings written to: ', csv_file, ''

      !> =================== Raw timing table ====================
      write (*, '(a)') ''
      write (*, '(a,i0,a)') '--- Average timings over ', n_iter, ' iterations (seconds) ---'

      !> Header: component name + N_atoms for each structure
      write (*, '(a30)', advance='no') 'Component'
      do istruct = 1, n_struct
         write (*, '(i10)', advance='no') n_atoms_arr(istruct)
      end do
      write (*, *)

      write (*, '(a30)', advance='no') repeat('-', 30)
      do istruct = 1, n_struct
         write (*, '(a10)', advance='no') '----------'
      end do
      write (*, *)

      !> Total row
      write (*, '(a30)', advance='no') 'TOTAL'
      do istruct = 1, n_struct
         write (*, '(f10.4)', advance='no') total_times(istruct)
      end do
      write (*, *)

      !> Per-timer rows (skip zero timers)
      do iprint = 1, n_print
         itimer = print_order(iprint)
         if (maxval(times(itimer, :)) < t_min) cycle
         if (is_parent(itimer)) then
            write (*, '(a30)', advance='no') repeat('-', 30)
            do istruct = 1, n_struct
               write (*, '(a10)', advance='no') '----------'
            end do
            write (*, *)
         end if
         write (*, '(a30)', advance='no') timer_labels(itimer)
         do istruct = 1, n_struct
            write (*, '(f10.4)', advance='no') times(itimer, istruct)
         end do
         write (*, *)
      end do

      !> ================== Scaling exponents ===================
      write (*, '(a)') ''
      write (*, '(a)') '--- Scaling fit: t(N) = A * N^X ---'
      write (*, '(a30, a10, a10, a12)') 'Component', 'X', 'R^2', 'A'
      write (*, '(a30, a10, a10, a12)') repeat('-', 30), &
         repeat('-', 10), repeat('-', 10), repeat('-', 12)

      !> Total
      call collect_valid_points(n_struct, n_atoms_arr, total_times, t_min, &
                                real_n, raw_t, n_valid)
      if (n_valid >= 4) then
         call fit_power_law(n_valid, real_n, raw_t, exponent, r_sq, &
                            coeff_a)
         write (*, '(a30, f10.3, f10.4, es12.3)') &
            'TOTAL', exponent, r_sq, coeff_a
      end if

      !> Per-timer
      do iprint = 1, n_print
         itimer = print_order(iprint)
         if (maxval(times(itimer, :)) < t_min) cycle
         if (is_parent(itimer)) then
            write (*, '(a30, a10, a10, a12)') repeat('-', 30), &
               repeat('-', 10), repeat('-', 10), repeat('-', 12)
         end if

         call collect_valid_points(n_struct, n_atoms_arr, times(itimer, :), t_min, &
                                   real_n, raw_t, n_valid)
         if (n_valid >= 4) then
            call fit_power_law(n_valid, real_n, raw_t, exponent, r_sq, &
                               coeff_a)
            write (*, '(a30, f10.3, f10.4, es12.3)') &
               timer_labels(itimer), exponent, r_sq, coeff_a
         else
            write (*, '(a30, a22)') timer_labels(itimer), '  (insufficient data)'
         end if
      end do

      write (*, '(a)') ''

   contains

      subroutine add_order(id)
         integer, intent(in) :: id
         n_print = n_print + 1
         print_order(n_print) = id
      end subroutine add_order

   end subroutine test_timing_drop_scaling

   !> Benchmark marching cubes integration scaling with system size.
   !> Uses the same polyalanine series as test_timing_drop_scaling.
   subroutine test_timing_mc_scaling(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(static_radius_type) :: radius_model
      type(mctc_error), allocatable :: radii_error

      !> Number of polyalanine structures to benchmark
      integer, parameter :: n_struct = 25
      !> Minimum time (s) for a data point to be included in the fit
      real(wp), parameter :: t_min = 1.0e-6_wp

      !> LSF parameters (match the DROP scaling benchmark)
      real(wp), parameter :: blend_k = 4.5_wp
      real(wp), parameter :: blend_3b = 1.0_wp

      !> Marching cubes target spacing
      real(wp), parameter :: mc_spacing = 0.2_wp

      character(len=20) :: struct_names(n_struct)

      !> Collected benchmark data
      integer :: n_atoms_arr(n_struct)
      real(wp) :: mc_times(n_struct)
      real(wp) :: mc_areas(n_struct), mc_volumes(n_struct)

      !> Number of repetitions per structure for stable timings
      integer, parameter :: n_iter = 3

      !> Fit workspace
      real(wp) :: real_n(n_struct), raw_t(n_struct)
      real(wp) :: exponent, r_sq, coeff_a
      integer :: n_valid

      integer :: istruct, iter
      real(wp) :: area, volume, rn_iter
      real :: t0, t1

      !> Polyalanine structures (increasing size)
      struct_names = [character(len=20) :: &
                      'polyala_04', 'polyala_08', 'polyala_12', 'polyala_16', 'polyala_20', &
                      'polyala_24', 'polyala_28', 'polyala_32', 'polyala_36', 'polyala_40', &
                      'polyala_44', 'polyala_48', 'polyala_52', 'polyala_56', 'polyala_60', &
                      'polyala_64', 'polyala_68', 'polyala_72', 'polyala_76', 'polyala_80', &
                      'polyala_84', 'polyala_88', 'polyala_92', 'polyala_96', 'polyala_100']

      rn_iter = real(n_iter, wp)

      write (*, '(a)') ''
      write (*, '(a)') '=================================================================='
      write (*, '(a)') 'Benchmark: Marching cubes scaling'
      write (*, '(a,f5.2,a,f5.2,a,f5.2,a,i0)') &
         'k=', blend_k, '  b3=', blend_3b, &
         '  spacing=', mc_spacing, '  iter=', n_iter
      write (*, '(a)') '=================================================================='
      write (*, '(a)') ''
      write (*, '(a14, a8, a12, a16, a16)') &
         'Structure', 'N_at', 'Time (s)', 'Area', 'Volume'
      write (*, '(a14, a8, a12, a16, a16)') &
         '-------------', '-------', '-----------', &
         '---------------', '---------------'

      radius_model = default_cpcm_radii()

      do istruct = 1, n_struct
         call get_structure(mol, 'POLYALANINE', trim(struct_names(istruct)))
         n_atoms_arr(istruct) = mol%nat

         call radius_model%update(mol, radii_error)
         if (allocated(radii_error)) then
            call test_failed(error, radii_error%message)
            return
         end if

         call lsf%new(blend_k=blend_k, blend_1b=1.0_wp, blend_2b=1.0_wp, &
                      blend_3b=blend_3b)
         lsf%screening_threshold = 0.0_wp
         call lsf%update(mol, radius_model%f0)

         mc_times(istruct) = 0.0_wp
         do iter = 1, n_iter
            call cpu_time(t0)
            call integrate_surface_marching_cubes(lsf, mol%xyz, area, volume, &
                                                  target_spacing=mc_spacing)
            call cpu_time(t1)
            mc_times(istruct) = mc_times(istruct) + real(t1 - t0, wp)
         end do
         mc_times(istruct) = mc_times(istruct)/rn_iter
         mc_areas(istruct) = area
         mc_volumes(istruct) = volume

         write (*, '(2x,a14, i6, f12.4, f16.4, f16.4)') &
            trim(struct_names(istruct)), n_atoms_arr(istruct), &
            mc_times(istruct), mc_areas(istruct), mc_volumes(istruct)
      end do

      !> Scaling fit
      write (*, '(a)') ''
      write (*, '(a)') '--- Scaling fit: t(N) = A * N^X ---'
      write (*, '(a10, a10, a10, a12)') 'Component', 'X', 'R^2', 'A'
      write (*, '(a10, a10, a10, a12)') repeat('-', 10), &
         repeat('-', 10), repeat('-', 10), repeat('-', 12)

      call collect_valid_points(n_struct, n_atoms_arr, mc_times, t_min, &
                                real_n, raw_t, n_valid)
      if (n_valid >= 4) then
         call fit_power_law(n_valid, real_n, raw_t, exponent, r_sq, &
                            coeff_a)
         write (*, '(a10, f10.3, f10.4, es12.3)') &
            'MC total', exponent, r_sq, coeff_a
      else
         write (*, '(a10, a22)') 'MC total', '  (insufficient data)'
      end if

      write (*, '(a)') ''

   end subroutine test_timing_mc_scaling

   !> Benchmark iSwiG cavity scaling with system size.
   !> Uses the same polyalanine series as test_timing_drop_scaling.
   subroutine test_timing_iswig_scaling(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_iswig), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error

      !> Number of polyalanine structures to benchmark
      integer, parameter :: n_struct = 25
      !> Minimum time (s) for a data point to be included in the fit
      real(wp), parameter :: t_min = 1.0e-6_wp

      !> iSwiG parameters
      integer, parameter :: nleb = 194

      character(len=20) :: struct_names(n_struct)

      !> Collected benchmark data
      integer :: n_atoms_arr(n_struct)
      integer :: n_grid_arr(n_struct)
      real(wp) :: total_times(n_struct)
      real(wp) :: areas(n_struct), volumes(n_struct)

      !> Number of repetitions per structure for stable timings
      integer, parameter :: n_iter = 10

      !> Fit workspace
      real(wp) :: real_n(n_struct), raw_t(n_struct)
      real(wp) :: exponent, r_sq, coeff_a
      integer :: n_valid

      integer :: istruct, iter
      real(wp) :: rn_iter
      real :: t0, t1

      !> Polyalanine structures (increasing size)
      struct_names = [character(len=20) :: &
                      'polyala_04', 'polyala_08', 'polyala_12', 'polyala_16', 'polyala_20', &
                      'polyala_24', 'polyala_28', 'polyala_32', 'polyala_36', 'polyala_40', &
                      'polyala_44', 'polyala_48', 'polyala_52', 'polyala_56', 'polyala_60', &
                      'polyala_64', 'polyala_68', 'polyala_72', 'polyala_76', 'polyala_80', &
                      'polyala_84', 'polyala_88', 'polyala_92', 'polyala_96', 'polyala_100']

      rn_iter = real(n_iter, wp)

      write (*, '(a)') ''
      write (*, '(a)') '=================================================================='
      write (*, '(a)') 'Benchmark: iSwiG cavity scaling (update + gradient)'
      write (*, '(a,i0,a,i0)') 'nleb=', nleb, '  iter=', n_iter
      write (*, '(a)') '=================================================================='
      write (*, '(a)') ''
      write (*, '(a14, a8, a8, a12, a16, a16)') &
         'Structure', 'N_at', 'N_grid', 'Time (s)', 'Area', 'Volume'
      write (*, '(a14, a8, a8, a12, a16, a16)') &
         '-------------', '-------', '-------', '-----------', &
         '---------------', '---------------'

      do istruct = 1, n_struct
         call get_structure(mol, 'POLYALANINE', trim(struct_names(istruct)))
         n_atoms_arr(istruct) = mol%nat

         total_times(istruct) = 0.0_wp
         do iter = 1, n_iter
            if (allocated(cavity)) deallocate (cavity)
            allocate (cavity)
            call new_cavity_iswig(cavity, nleb=nleb, &
                                  radius_model=default_cpcm_radii(), error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if

            call cpu_time(t0)
            call cavity%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call cavity%get_gradient()
            call cpu_time(t1)
            total_times(istruct) = total_times(istruct) + real(t1 - t0, wp)
         end do
         total_times(istruct) = total_times(istruct)/rn_iter
         n_grid_arr(istruct) = cavity%ngrid
         areas(istruct) = cavity%total_area
         volumes(istruct) = cavity%total_volume

         write (*, '(2x,a14, i6, i8, f12.4, f16.4, f16.4)') &
            trim(struct_names(istruct)), n_atoms_arr(istruct), &
            n_grid_arr(istruct), total_times(istruct), &
            areas(istruct), volumes(istruct)
      end do

      if (allocated(cavity)) deallocate (cavity)

      !> Scaling fit
      write (*, '(a)') ''
      write (*, '(a)') '--- Scaling fit: t(N) = A * N^X ---'
      write (*, '(a14, a10, a10, a12)') 'Component', 'X', 'R^2', 'A'
      write (*, '(a14, a10, a10, a12)') repeat('-', 14), &
         repeat('-', 10), repeat('-', 10), repeat('-', 12)

      call collect_valid_points(n_struct, n_atoms_arr, total_times, t_min, &
                                real_n, raw_t, n_valid)
      if (n_valid >= 4) then
         call fit_power_law(n_valid, real_n, raw_t, exponent, r_sq, &
                            coeff_a)
         write (*, '(a14, f10.3, f10.4, es12.3)') &
            'iSwiG total', exponent, r_sq, coeff_a
      else
         write (*, '(a14, a22)') 'iSwiG total', '  (insufficient data)'
      end if

      write (*, '(a)') ''

   end subroutine test_timing_iswig_scaling

   ! !> Benchmark GEPOL SES cavity scaling with system size.
   ! !> Uses the same polyalanine series as test_timing_drop_scaling.
   ! subroutine test_timing_gepol_scaling(error)
   !    type(error_type), allocatable, intent(out) :: error
   !    type(structure_type) :: mol
   !    type(cavity_type_gepol), allocatable :: cavity
   !    type(mctc_error), allocatable :: cavity_error

   !    !> Number of polyalanine structures to benchmark
   !    integer, parameter :: n_struct = 25
   !    !> Minimum time (s) for a data point to be included in the fit
   !    real(wp), parameter :: t_min = 1.0e-6_wp

   !    !> GEPOL parameters
   !    integer, parameter :: ndiv = 3

   !    character(len=20) :: struct_names(n_struct)

   !    !> Collected benchmark data
   !    integer :: n_atoms_arr(n_struct)
   !    integer :: n_grid_arr(n_struct)
   !    real(wp) :: total_times(n_struct)
   !    real(wp) :: areas(n_struct), volumes(n_struct)

   !    !> Number of repetitions per structure for stable timings
   !    integer, parameter :: n_iter = 10

   !    !> Fit workspace
   !    real(wp) :: real_n(n_struct), raw_t(n_struct)
   !    real(wp) :: exponent, r_sq, coeff_a
   !    integer :: n_valid

   !    integer :: istruct, iter
   !    real(wp) :: rn_iter
   !    real :: t0, t1

   !    !> Polyalanine structures (increasing size)
   !    struct_names = [character(len=20) :: &
   !                    'polyala_04', 'polyala_08', 'polyala_12', 'polyala_16', 'polyala_20', &
   !                    'polyala_24', 'polyala_28', 'polyala_32', 'polyala_36', 'polyala_40', &
   !                    'polyala_44', 'polyala_48', 'polyala_52', 'polyala_56', 'polyala_60', &
   !                    'polyala_64', 'polyala_68', 'polyala_72', 'polyala_76', 'polyala_80', &
   !                    'polyala_84', 'polyala_88', 'polyala_92', 'polyala_96', 'polyala_100']

   !    rn_iter = real(n_iter, wp)

   !    write (*, '(a)') ''
   !    write (*, '(a)') '=================================================================='
   !    write (*, '(a)') 'Benchmark: GEPOL SES cavity scaling (update + gradient)'
   !    write (*, '(a,i0,a,i0)') 'ndiv=', ndiv, '  iter=', n_iter
   !    write (*, '(a)') '=================================================================='
   !    write (*, '(a)') ''
   !    write (*, '(a14, a8, a8, a12, a16, a16)') &
   !       'Structure', 'N_at', 'N_grid', 'Time (s)', 'Area', 'Volume'
   !    write (*, '(a14, a8, a8, a12, a16, a16)') &
   !       '-------------', '-------', '-------', '-----------', &
   !       '---------------', '---------------'

   !    do istruct = 1, n_struct
   !       call get_structure(mol, 'POLYALANINE', trim(struct_names(istruct)))
   !       n_atoms_arr(istruct) = mol%nat

   !       total_times(istruct) = 0.0_wp
   !       do iter = 1, n_iter
   !          if (allocated(cavity)) deallocate (cavity)
   !          allocate (cavity)
   !          call new_cavity_gepol(cavity, ndiv=ndiv, verbosity=0, &
   !                                radius_model=default_cpcm_radii(), error=cavity_error)
   !          if (allocated(cavity_error)) then
   !             call test_failed(error, cavity_error%message)
   !             return
   !          end if

   !          call cpu_time(t0)
   !          call cavity%update(mol, error=cavity_error)
   !          if (allocated(cavity_error)) then
   !             call test_failed(error, cavity_error%message)
   !             return
   !          end if
   !          call cavity%get_gradient()
   !          call cpu_time(t1)
   !          total_times(istruct) = total_times(istruct) + real(t1 - t0, wp)
   !       end do
   !       total_times(istruct) = total_times(istruct)/rn_iter
   !       n_grid_arr(istruct) = cavity%ngrid
   !       areas(istruct) = cavity%total_area
   !       volumes(istruct) = cavity%total_volume

   !       write (*, '(2x,a14, i6, i8, f12.4, f16.4, f16.4)') &
   !          trim(struct_names(istruct)), n_atoms_arr(istruct), &
   !          n_grid_arr(istruct), total_times(istruct), &
   !          areas(istruct), volumes(istruct)
   !    end do

   !    if (allocated(cavity)) deallocate (cavity)

   !    !> Scaling fit
   !    write (*, '(a)') ''
   !    write (*, '(a)') '--- Scaling fit: t(N) = A * N^X ---'
   !    write (*, '(a14, a10, a10, a12)') 'Component', 'X', 'R^2', 'A'
   !    write (*, '(a14, a10, a10, a12)') repeat('-', 14), &
   !       repeat('-', 10), repeat('-', 10), repeat('-', 12)

   !    call collect_valid_points(n_struct, n_atoms_arr, total_times, t_min, &
   !                              real_n, raw_t, n_valid)
   !    if (n_valid >= 4) then
   !       call fit_power_law(n_valid, real_n, raw_t, exponent, r_sq, &
   !                          coeff_a)
   !       write (*, '(a14, f10.3, f10.4, es12.3)') &
   !          'GEPOL total', exponent, r_sq, coeff_a
   !    else
   !       write (*, '(a14, a22)') 'GEPOL total', '  (insufficient data)'
   !    end if

   !    write (*, '(a)') ''

   ! end subroutine test_timing_gepol_scaling

   !> Collect valid data points for power-law fitting.
   !> Only includes points where the measured time exceeds t_min.
   !>
   !> @param[in]  n_struct  Number of structures
   !> @param[in]  n_atoms   Atom counts per structure
   !> @param[in]  times     Measured times per structure (s)
   !> @param[in]  t_min     Minimum time threshold
   !> @param[out] real_n    Valid atom counts (as real)
   !> @param[out] raw_t     Valid measured times
   !> @param[out] n_valid   Number of valid data points
   pure subroutine collect_valid_points(n_struct, n_atoms, times, t_min, &
                                        real_n, raw_t, n_valid)
      integer, intent(in) :: n_struct
      integer, intent(in) :: n_atoms(n_struct)
      real(wp), intent(in) :: times(n_struct)
      real(wp), intent(in) :: t_min
      real(wp), intent(out) :: real_n(n_struct), raw_t(n_struct)
      integer, intent(out) :: n_valid
      integer :: i

      n_valid = 0
      do i = 1, n_struct
         if (times(i) > t_min) then
            n_valid = n_valid + 1
            real_n(n_valid) = real(n_atoms(i), wp)
            raw_t(n_valid) = times(i)
         end if
      end do
   end subroutine collect_valid_points

   !> Fit power law: t(N) = A * N^X  (constrained through origin).
   !>
   !> Uses a two-pass grid search over X. For each candidate X the
   !> model is linear in A, so the optimal prefactor is obtained as
   !> A = sum(N_i^X * t_i) / sum(N_i^{2X}).
   !>
   !> Pass 1: X in [0, 4], step 0.05  (coarse scan)
   !> Pass 2: X in [best-0.05, best+0.05], step 0.001  (refinement)
   !>
   !> @param[in]  n          Number of data points (must be >= 4)
   !> @param[in]  real_n     System sizes (atom counts)
   !> @param[in]  raw_t      Measured times (s)
   !> @param[out] exponent   Scaling exponent X
   !> @param[out] r_squared  Coefficient of determination
   !> @param[out] coeff_a    Prefactor A
   pure subroutine fit_power_law(n, real_n, raw_t, exponent, r_squared, &
                                 coeff_a)
      integer, intent(in) :: n
      real(wp), intent(in) :: real_n(n), raw_t(n)
      real(wp), intent(out) :: exponent, r_squared, coeff_a

      real(wp) :: log_n_pre(n), nix(n)
      real(wp) :: s_n2x, s_nxt
      real(wp) :: a_cand, sse, best_sse
      real(wp) :: best_x, best_a
      real(wp) :: x_lo, x_hi, dx, x_cand
      real(wp) :: t_mean, sst
      integer :: i, pass

      t_mean = sum(raw_t(1:n))/real(n, wp)

      do i = 1, n
         log_n_pre(i) = log(real_n(i))
      end do

      best_sse = huge(1.0_wp)
      best_x = 1.0_wp
      best_a = 0.0_wp

      !> Two-pass grid search
      do pass = 1, 2
         if (pass == 1) then
            x_lo = 0.0_wp; x_hi = 4.0_wp; dx = 0.05_wp
         else
            x_lo = max(0.0_wp, best_x - 0.05_wp)
            x_hi = min(4.0_wp, best_x + 0.05_wp)
            dx = 0.001_wp
         end if

         x_cand = x_lo
         do while (x_cand <= x_hi + 0.5_wp*dx)
            do i = 1, n
               nix(i) = exp(x_cand*log_n_pre(i))
            end do

            s_n2x = sum(nix(1:n)**2)
            s_nxt = sum(nix(1:n)*raw_t(1:n))

            if (s_n2x < 1.0e-30_wp) then
               x_cand = x_cand + dx
               cycle
            end if

            a_cand = s_nxt/s_n2x

            !> Enforce A >= 0 (time contribution must be non-negative)
            if (a_cand < 0.0_wp) then
               x_cand = x_cand + dx
               cycle
            end if

            sse = 0.0_wp
            do i = 1, n
               sse = sse + (raw_t(i) - a_cand*nix(i))**2
            end do

            if (sse < best_sse) then
               best_sse = sse
               best_x = x_cand
               best_a = a_cand
            end if

            x_cand = x_cand + dx
         end do
      end do

      !> Compute R^2
      sst = 0.0_wp
      do i = 1, n
         sst = sst + (raw_t(i) - t_mean)**2
      end do

      exponent = best_x
      coeff_a = best_a
      if (sst > 1.0e-30_wp) then
         r_squared = 1.0_wp - best_sse/sst
      else
         r_squared = 0.0_wp
      end if
   end subroutine fit_power_law

end module test_cavity_drop_timings
