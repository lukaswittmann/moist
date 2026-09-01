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
   use moist_radii_static, only: radius_type_static
   use moist_cavity_marchingcubes, only: integrate_surface_marching_cubes
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   implicit none
   private

   public :: collect_cavity_drop_timings
   !> Shared with the other scaling benchmarks in test/dev
   public :: collect_valid_points, fit_power_law

   integer, parameter :: ndim = 3

contains

   subroutine collect_cavity_drop_timings(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  ! new_unittest("timing_drop_proj_levels", test_timing_drop_proj_levels), &
                  ! new_unittest("timing_cell_fraction_benchmark", test_timing_cell_fraction_benchmark), &
                  new_unittest("timing_drop_scaling", test_timing_drop_scaling), &
                  new_unittest("timing_mc_scaling", test_timing_mc_scaling), &
                  new_unittest("timing_iswig_scaling", test_timing_iswig_scaling), &
                  new_unittest("timing_lsf_accessors", test_timing_lsf_accessors), &
                  new_unittest("timing_prepare_split", test_timing_prepare_split), &
                  new_unittest("timing_cavity_build", test_timing_cavity_build) &
                  ]
   end subroutine collect_cavity_drop_timings

   !> Time a serial cavity build and report the geometry it produced.
   !>
   !> Area and volume are printed alongside the timing so that a change aimed
   !> at the projection can be checked for having left the cavity alone --
   !> which is the whole premise of relaxing anything in the *seed* stage.
   subroutine test_timing_cavity_build(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n_struct = 3
      character(len=20), parameter :: struct_names(n_struct) = &
                                      [character(len=20) :: 'polyala_04', 'polyala_16', &
                                       'polyala_32']
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(radius_type_static) :: radius_model
      type(mctc_error), allocatable :: cavity_error
      type(moist_context_type), target :: ctx
      integer :: istruct
      real(wp) :: t_build
      real :: c0, c1

      call new_context(ctx, verbosity=0)
      radius_model = default_cpcm_radii()

      write (*, '(a)') ''
      write (*, '(a)') '=================================================================='
      write (*, '(a)') 'Cavity build time and resulting geometry (OMP_NUM_THREADS=1)'
      write (*, '(a)') '=================================================================='
      write (*, '(a)') ''
      write (*, '(2x,a12,a8,a12,a15,a15)') &
         'Structure', 'N_at', 'build(s)', 'area', 'volume'
      write (*, '(2x,a12,a8,a12,a15,a15)') &
         '-----------', '-------', '-----------', '--------------', '--------------'

      do istruct = 1, n_struct
         call get_structure(mol, 'POLYALANINE', trim(struct_names(istruct)))


         allocate (cavity)
         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
            call svdw_template%new(blend_k=5.5_wp, blend_2b=0.0_wp, blend_3b=3.0_wp)
            call new_cavity_drop(cavity, ctx, radius_model=radius_model, &
                                 lsf_model=svdw_template, error=cavity_error)
         end block
         if (allocated(cavity_error)) then
            call test_failed(error, "cavity build failed: "//cavity_error%message)
            return
         end if
         call cpu_time(c0)
         call cavity%update(mol, error=cavity_error)
         call cpu_time(c1)
         if (allocated(cavity_error)) then
            call test_failed(error, "cavity update failed: "//cavity_error%message)
            return
         end if
         t_build = real(c1 - c0, wp)

         write (*, '(2x,a12,i8,f12.4,2f15.6)') &
            trim(struct_names(istruct)), mol%nat, t_build, &
            cavity%total_area, cavity%total_volume

         deallocate (cavity)
      end do

      write (*, '(a)') ''
      write (*, '(a)') 'Serial cavity build; area and volume pin the geometry down'
      write (*, '(a)') ''

   end subroutine test_timing_cavity_build

   !> Split `prepare` into its screening and accumulation halves.
   !>
   !> `prepare` does two things: reject candidates that cannot contribute
   !> (`screen_candidates`, one squared-distance compare per candidate) and
   !> accumulate the power sums over the survivors (one `exp` and the kind
   !> tensors per active atom). Only the second half depends on the evaluation
   !> point continuously; the first is a set membership that barely moves
   !> between two nearby points. Knowing the split therefore decides whether
   !> freezing the active set across solver iterations is worth anything.
   !>
   !> Both halves are measured through the public API, by timing the same
   !> points twice with different candidate-list lengths:
   !>
   !>   `prepare`         screens all `N` atoms       t_full = a*N     + acc
   !>   `prepare_subset`  screens the grid candidates t_grid = a*ncand + acc
   !>
   !> `acc` is identical in the two because the surviving set is the same, so
   !> the pair determines both unknowns. `a` must come out independent of the
   !> derivative order (screening does no derivative work), which is printed as
   !> a consistency check on the model rather than assumed.
   !>
   !> The cell grid is built exactly as `setup_mol_cell_grid` builds it, and the
   !> query is timed separately so the grid lookup is not charged to screening.
   !> SvdW only: CFC does not remap the cell grid into its candidate space.
   subroutine test_timing_prepare_split(error)
      type(error_type), allocatable, intent(out) :: error

      !> Points per structure
      integer, parameter :: npts = 2000
      !> Highest derivative level measured
      integer, parameter :: max_lvl = 3
      integer, parameter :: n_struct = 6
      !> Mirrors moist_cavity_drop_parameters defaults
      integer, parameter :: full_scan_below = 200
      real(wp), parameter :: cell_fraction = 0.25_wp
      real(wp), parameter :: screen_thr = 1.0e-11_wp
      real(wp), parameter :: blend_k = 5.5_wp, blend_2b = 0.0_wp, blend_3b = 3.0_wp

      character(len=20), parameter :: struct_names(n_struct) = &
                                      [character(len=20) :: 'polyala_04', 'polyala_16', &
                                       'polyala_32', 'polyala_52', 'polyala_76', 'polyala_100']

      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cell_grid_type) :: cell_grid
      type(mctc_error), allocatable :: lsf_err
      real(wp), allocatable :: radii(:), r_eff(:), points(:, :)
      real(wp) :: t_full(0:max_lvl), t_grid(0:max_lvl), t_query
      real(wp) :: a_screen(0:max_lvl), t_acc(0:max_lvl)
      real(wp) :: ncand_mean, nact_mean
      integer :: istruct, ilvl, ipt, start, ncand, cand_sum, act_sum
      integer :: nat
      real :: c0, c1

      write (*, '(a)') ''
      write (*, '(a)') '=================================================================='
      write (*, '(a)') 'Benchmark: SvdW prepare split (screening vs accumulation)'
      write (*, '(a,i0)') 'points/system: ', npts
      write (*, '(a,i0,a,f5.2)') 'cell grid: full_scan_below=', full_scan_below, &
         '  cell_fraction=', cell_fraction
      write (*, '(a)') '=================================================================='
      write (*, '(a)') ''
      write (*, '(2x,a12,a6,a9,a9,a4,a11,a11,a11,a11)') &
         'Structure', 'N_at', 'ncand', 'n_act', 'd', 'full(us)', 'grid(us)', &
         'screen/cand', 'accum(us)'
      write (*, '(2x,a12,a6,a9,a9,a4,a11,a11,a11,a11)') &
         '-----------', '-----', '--------', '--------', '---', &
         '----------', '----------', '----------', '----------'

      do istruct = 1, n_struct
         call get_structure(mol, 'POLYALANINE', trim(struct_names(istruct)))
         nat = mol%nat
         call fill_cpcm_radii(mol, radii, error)
         if (allocated(error)) return
         radii = radii*aatoau
         allocate (r_eff(nat), points(ndim, npts))
         call build_lsf_shell_points(mol, radii, points)

         call lsf%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
         lsf%screening_threshold = screen_thr
         call lsf%update(mol, radii)

         ! Same construction as setup_mol_cell_grid, then relabel the candidate
         ! ids into the LSF's internal (sorted) atom ordering
         do ipt = 1, nat
            r_eff(ipt) = radii(ipt) + lsf%neighbor_cutoff(radii(ipt))
         end do
         call cell_grid%build(mol%xyz, r_eff, full_scan_below=full_scan_below, &
                              cell_fraction=cell_fraction)
         if (allocated(cell_grid%cell_nlat)) then
            call lsf%remap_candidate_grid(cell_grid%cell_nlat)
         end if

         !> Candidate and active counts (both fixed across derivative levels)
         call lsf%set_max_deriv(0)
         cand_sum = 0
         act_sum = 0
         do ipt = 1, npts
            call cell_grid%query(points(:, ipt), start, ncand)
            cand_sum = cand_sum + ncand
            call lsf%prepare_subset(points(:, ipt), &
                                    cell_grid%cell_nlat(start + 1:start + ncand), lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "prepare_subset failed: "//lsf_err%message)
               return
            end if
            act_sum = act_sum + lsf%active_count()
         end do
         ncand_mean = real(cand_sum, wp)/real(npts, wp)
         nact_mean = real(act_sum, wp)/real(npts, wp)

         !> Grid query alone, so the lookup is not charged to screening
         call cpu_time(c0)
         do ipt = 1, npts
            call cell_grid%query(points(:, ipt), start, ncand)
         end do
         call cpu_time(c1)
         t_query = real(c1 - c0, wp)/real(npts, wp)*1.0e6_wp

         do ilvl = 0, max_lvl
            call lsf%set_max_deriv(ilvl)

            call cpu_time(c0)
            do ipt = 1, npts
               call lsf%prepare(points(:, ipt), lsf_err)
            end do
            call cpu_time(c1)
            t_full(ilvl) = real(c1 - c0, wp)/real(npts, wp)*1.0e6_wp

            call cpu_time(c0)
            do ipt = 1, npts
               call cell_grid%query(points(:, ipt), start, ncand)
               call lsf%prepare_subset(points(:, ipt), &
                                       cell_grid%cell_nlat(start + 1:start + ncand), lsf_err)
            end do
            call cpu_time(c1)
            t_grid(ilvl) = real(c1 - c0, wp)/real(npts, wp)*1.0e6_wp - t_query

            ! t_full = a*N + acc, t_grid = a*ncand + acc. Below
            ! `full_scan_below` the grid is a single cell, so the two
            ! measurements coincide and the split is not determined.
            if (real(nat, wp) - ncand_mean > 1.0_wp) then
               a_screen(ilvl) = (t_full(ilvl) - t_grid(ilvl))/(real(nat, wp) - ncand_mean)
               t_acc(ilvl) = t_grid(ilvl) - a_screen(ilvl)*ncand_mean
            else
               a_screen(ilvl) = -1.0_wp
               t_acc(ilvl) = -1.0_wp
            end if
         end do

         do ilvl = 0, max_lvl
            if (ilvl == 0) then
               write (*, '(2x,a12,i6,f9.1,f9.1,i4,4f11.4)') &
                  trim(struct_names(istruct)), nat, ncand_mean, nact_mean, ilvl, &
                  t_full(ilvl), t_grid(ilvl), a_screen(ilvl), t_acc(ilvl)
            else
               write (*, '(2x,a12,a6,a9,a9,i4,4f11.4)') &
                  '', '', '', '', ilvl, &
                  t_full(ilvl), t_grid(ilvl), a_screen(ilvl), t_acc(ilvl)
            end if
         end do

         call cell_grid%destroy()
         deallocate (radii, r_eff, points)
      end do

      write (*, '(a)') ''
      write (*, '(a)') 'full  = prepare (screens every atom), grid = prepare_subset'
      write (*, '(a)') '        (screens the cell-grid candidates), query excluded'
      write (*, '(a)') 'accum = grid - screen/cand * ncand: the point-dependent half'
      write (*, '(a)') '  "-1" = grid is a single cell (N_at < full_scan_below),'
      write (*, '(a)') '        so full and grid coincide and the split is undetermined'
      write (*, '(a)') ''

      call sweep_screening_threshold(error)

   end subroutine test_timing_prepare_split

   !> How much does the active set actually cost in accuracy?
   !>
   !> The SvdW gate keeps every atom whose one-body weight
   !> `u = exp(-(k/3)(x-R))` exceeds `screening_threshold`, a reach of
   !> `-3 ln(thr)/k` Bohr that is the same for every atom. The gate is therefore
   !> phrased in a *one-body* quantity, while the error it controls is the shift
   !> in `S = -(1/k) ln Z` -- and `Z` is not 1 away from the surface. This sweep
   !> measures what the absolute gate buys and costs: active count, `prepare`
   !> cost, and the resulting error in `S` against an effectively unscreened
   !> reference, over the shell points where a cavity surface would sit.
   !>
   !> `dS` is reported in absolute terms on purpose. The projection solves
   !> `S = 0`, so a shift in `S` displaces the surface by `dS/|grad S|`; a
   !> relative measure against `S ~ 0` would be meaningless there.
   subroutine sweep_screening_threshold(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: npts = 2000
      integer, parameter :: n_thr = 6
      real(wp), parameter :: thr_list(n_thr) = &
                             [1.0e-11_wp, 1.0e-9_wp, 1.0e-7_wp, 1.0e-5_wp, &
                              1.0e-4_wp, 1.0e-3_wp]
      !> Effectively unscreened reference
      real(wp), parameter :: thr_ref = 1.0e-15_wp
      real(wp), parameter :: blend_k = 5.5_wp, blend_2b = 0.0_wp, blend_3b = 3.0_wp
      integer, parameter :: full_scan_below = 200
      real(wp), parameter :: cell_fraction = 0.25_wp

      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cell_grid_type) :: cell_grid
      type(mctc_error), allocatable :: lsf_err
      real(wp), allocatable :: radii(:), r_eff(:), points(:, :), s_ref(:)
      real(wp) :: t_prep, nact_mean, ncand_mean, ds_max, ds_mean, offset, val
      integer :: ithr, ipt, start, ncand, act_sum, cand_sum, nat
      real :: c0, c1

      write (*, '(a)') '=================================================================='
      write (*, '(a)') 'Sweep: SvdW screening_threshold (polyala_100, max_deriv=2)'
      write (*, '(a)') '=================================================================='
      write (*, '(a)') ''

      call get_structure(mol, 'POLYALANINE', 'polyala_100')
      nat = mol%nat
      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return
      radii = radii*aatoau
      allocate (r_eff(nat), points(ndim, npts), s_ref(npts))
      call build_lsf_shell_points(mol, radii, points)

      !> Reference S at a threshold far below any of the sweep values
      call lsf%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
      lsf%screening_threshold = thr_ref
      call lsf%update(mol, radii)
      call lsf%set_max_deriv(0)
      do ipt = 1, npts
         call lsf%prepare(points(:, ipt), lsf_err)
         call lsf%f0(s_ref(ipt))
      end do

      write (*, '(2x,a10,a10,a9,a9,a12,a13,a13)') &
         'threshold', 'reach', 'ncand', 'n_act', 'prep_d2(us)', 'max|dS|', 'mean|dS|'
      write (*, '(2x,a10,a10,a9,a9,a12,a13,a13)') &
         '---------', '--------', '--------', '--------', '-----------', &
         '------------', '------------'

      do ithr = 1, n_thr
         call lsf%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
         lsf%screening_threshold = thr_list(ithr)
         call lsf%update(mol, radii)
         offset = lsf%screening_offset(radii(1))

         do ipt = 1, nat
            r_eff(ipt) = radii(ipt) + lsf%neighbor_cutoff(radii(ipt))
         end do
         call cell_grid%build(mol%xyz, r_eff, full_scan_below=full_scan_below, &
                              cell_fraction=cell_fraction)
         if (allocated(cell_grid%cell_nlat)) then
            call lsf%remap_candidate_grid(cell_grid%cell_nlat)
         end if

         !> Counts and accuracy
         call lsf%set_max_deriv(0)
         act_sum = 0
         cand_sum = 0
         ds_max = 0.0_wp
         ds_mean = 0.0_wp
         do ipt = 1, npts
            call cell_grid%query(points(:, ipt), start, ncand)
            cand_sum = cand_sum + ncand
            call lsf%prepare_subset(points(:, ipt), &
                                    cell_grid%cell_nlat(start + 1:start + ncand), lsf_err)
            act_sum = act_sum + lsf%active_count()
            call lsf%f0(val)
            ds_max = max(ds_max, abs(val - s_ref(ipt)))
            ds_mean = ds_mean + abs(val - s_ref(ipt))
         end do
         nact_mean = real(act_sum, wp)/real(npts, wp)
         ncand_mean = real(cand_sum, wp)/real(npts, wp)
         ds_mean = ds_mean/real(npts, wp)

         !> prepare cost at the order the projection's Newton runs at
         call lsf%set_max_deriv(2)
         call cpu_time(c0)
         do ipt = 1, npts
            call cell_grid%query(points(:, ipt), start, ncand)
            call lsf%prepare_subset(points(:, ipt), &
                                    cell_grid%cell_nlat(start + 1:start + ncand), lsf_err)
         end do
         call cpu_time(c1)
         t_prep = real(c1 - c0, wp)/real(npts, wp)*1.0e6_wp

         write (*, '(2x,es10.1,f10.2,f9.1,f9.1,f12.4,es13.3,es13.3)') &
            thr_list(ithr), offset, ncand_mean, nact_mean, t_prep, ds_max, ds_mean

         call cell_grid%destroy()
      end do

      write (*, '(a)') ''
      write (*, '(a)') 'reach = -3*ln(thr)/k, the radial offset beyond each atom radius'
      write (*, '(a)') 'dS    = shift in the level-set value against thr=1e-15;'
      write (*, '(a)') '        the surface moves by dS/|grad S|'
      write (*, '(a)') ''

   end subroutine sweep_screening_threshold

   subroutine test_timing_drop_proj_levels(error)
      use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity_drop
      real(wp), allocatable :: radii(:)
      character(len=20) :: struct_names(25)
      integer :: iter, istruct, ilevel
      real :: t0, t1
      real(wp) :: time_update, time_gradient, time_total
      integer :: n_iter, nleb, proj_level
      real(wp) :: blend_k, blend_2b, blend_3b
      type(mctc_error), allocatable :: cavity_error
      integer, parameter :: n_proj_levels = 1
      integer, parameter :: proj_levels(n_proj_levels) = [1]
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, verbosity=0)

      nleb = 194
      blend_k = 5.5_wp
      blend_2b = 0.0_wp
      blend_3b = 3.0_wp
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
      write (*, '(a, f6.3)') '            blend_k  = ', blend_k
      write (*, '(a, f6.3)') '            blend_2b = ', blend_2b
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
         call fill_cpcm_radii(mol, radii, error)
         if (allocated(error)) return
         radii = radii*aatoau

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
                  call svdw_template%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
                  call new_cavity_drop(cavity_drop, ctx, nleb=nleb, &
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
               call cavity_drop%get_gradient(cavity_error)
               if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
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
      type(mctc_error), allocatable :: lsf_err

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
      integer :: istruct, ipt, iter, ifrac, n_calls
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
         call fill_cpcm_radii(mol, radii, error)
         if (allocated(error)) return
         allocate (r_eff(mol%nat))
         radii = radii*aatoau
         r_eff = radii + delta

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
               call lsf_prim%prepare(point, lsf_err)
               call lsf_prim%f012_r(lsf_val, lsf_grad, lsf_hess)
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
                                               cell_grid%cell_nlat(start + 1:start + n_cand), &
                                               lsf_err)
                  call lsf_prim%f012_r(lsf_val, lsf_grad, lsf_hess)
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
      !> Upper bound on the number of timer nodes; the actual count is read from
      !> the cavity timer tree at runtime into n_timers.
      integer, parameter :: max_timers = 64
      !> Minimum time (s) for a data point to be included in the fit
      real(wp), parameter :: t_min = 1.0e-6_wp

      !> Cavity parameters
      integer, parameter :: nleb = 194
      real(wp), parameter :: blend_k = 5.5_wp
      real(wp), parameter :: blend_2b = 0.0_wp
      real(wp), parameter :: blend_3b = 3.0_wp
      integer, parameter :: proj_level = 3

      character(len=20) :: struct_names(n_struct)
      !> Display labels (indented by tree depth) and top-level flags, enumerated
      !> from the cavity timer itself rather than a hardcoded table.
      integer :: n_timers
      character(len=40) :: timer_labels(max_timers)
      logical :: is_parent(max_timers)

      !> Collected benchmark data
      integer :: n_atoms_arr(n_struct)
      integer :: n_grid_arr(n_struct)
      real(wp) :: times(max_timers, n_struct)
      real(wp) :: total_times(n_struct)

      !> Number of repetitions per structure for stable timings
      integer, parameter :: n_iter = 3

      !> Fit workspace
      real(wp) :: real_n(n_struct), raw_t(n_struct)
      real(wp) :: exponent, r_sq, coeff_a
      integer :: n_valid

      integer :: istruct, itimer, iter, iprint
      real(wp) :: rn_iter

      !> CSV output
      integer :: csv_unit
      character(len=*), parameter :: csv_file = 'drop_timings.csv'
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, verbosity=0, do_profile=.true.)

      !> Polyalanine structures (increasing size)
      struct_names = [character(len=20) :: &
                      'polyala_04', 'polyala_08', 'polyala_12', 'polyala_16', 'polyala_20', &
                      'polyala_24', 'polyala_28', 'polyala_32', 'polyala_36', 'polyala_40', &
                      'polyala_44', 'polyala_48', 'polyala_52', 'polyala_56', 'polyala_60', &
                      'polyala_64', 'polyala_68', 'polyala_72', 'polyala_76', 'polyala_80', &
                      'polyala_84', 'polyala_88', 'polyala_92', 'polyala_96', 'polyala_100']

      !> The timer node set (names, nesting, order) is not declared here: it is
      !> enumerated from the cavity's own timer after the first build, so this
      !> benchmark automatically tracks whatever the DROP source measures.
      is_parent = .false.
      n_timers = 0

      !> ====================== Benchmark loop ======================
      rn_iter = real(n_iter, wp)

      write (*, '(a)') ''
      write (*, '(a)') '=================================================================='
      write (*, '(a)') 'Benchmark: DROP cavity component scaling (update + gradient)'
      write (*, '(a,i0,a,f5.2,a,f5.2,a,f5.2,a,i0,a,i0)') &
         'nleb=', nleb, '  k=', blend_k, '  b2=', blend_2b, '  b3=', blend_3b, &
         '  proj=', proj_level, '  iter=', n_iter
      write (*, '(a)') '=================================================================='

      !> Open the CSV file; its header is written once the timer tree is known
      !> (after the first cavity builds, so columns match what was measured).
      open (newunit=csv_unit, file=csv_file, status='replace', action='write')

      do istruct = 1, n_struct
         call get_structure(mol, 'POLYALANINE', trim(struct_names(istruct)))
         n_atoms_arr(istruct) = mol%nat

         !> Zero accumulators
         times(:, istruct) = 0.0_wp
         total_times(istruct) = 0.0_wp

         do iter = 1, n_iter
            if (allocated(cavity)) deallocate (cavity)
            allocate (cavity)
            !> Zero the shared timer so every node reports this build alone. The
            !> node tree (and hence the itimer -> label mapping snapshotted below)
            !> survives the reset, unlike a freshly constructed context.
            call ctx%timer%reset()
            block
               type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
               call svdw_template%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
               call new_cavity_drop(cavity, ctx, nleb=nleb, &
                                    do_fine=.true., &
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

            call cavity%get_gradient(cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if

            !> On the first build, enumerate the timer tree (labels + structure)
            !> and write the now-known CSV header.
            if (n_timers == 0) then
               call snapshot_timer_tree()
               write (csv_unit, '(a)', advance='no') 'structure,n_atoms,n_grid,iter,total'
               do itimer = 1, n_timers
                  write (csv_unit, '(a,a)', advance='no') ',', trim(adjustl(timer_labels(itimer)))
               end do
               write (csv_unit, *)
            end if

            !> Write per-iteration row to CSV
            write (csv_unit, '(a,a,i0,a,i0,a,i0,a,es14.6)', advance='no') &
               trim(struct_names(istruct)), ',', mol%nat, ',', cavity%ngrid, &
               ',', iter, ',', cavity%ctx%timer%get()
            do itimer = 1, n_timers
               write (csv_unit, '(a,es14.6)', advance='no') ',', &
                  cavity%ctx%timer%node_time(itimer)
            end do
            write (csv_unit, *)

            !> Accumulate timer values
            do itimer = 1, n_timers
               times(itimer, istruct) = times(itimer, istruct) &
                                        + cavity%ctx%timer%node_time(itimer)
            end do
            total_times(istruct) = total_times(istruct) + cavity%ctx%timer%get()
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
      do iprint = 1, n_timers
         itimer = iprint
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
      do iprint = 1, n_timers
         itimer = iprint
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

      !> Enumerate the cavity timer tree into the display arrays: label indented
      !> by nesting depth, top-level nodes flagged as parents. Uses the timer's
      !> own introspection so nothing about the node set is duplicated here.
      subroutine snapshot_timer_tree()
         integer :: id, depth
         n_timers = min(cavity%ctx%timer%num_nodes(), max_timers)
         do id = 1, n_timers
            depth = cavity%ctx%timer%node_depth(id)
            timer_labels(id) = repeat('  ', depth)//cavity%ctx%timer%node_name(id)
            is_parent(id) = depth == 0
         end do
      end subroutine snapshot_timer_tree

   end subroutine test_timing_drop_scaling

   !> Benchmark marching cubes integration scaling with system size.
   !> Uses the same polyalanine series as test_timing_drop_scaling.
   subroutine test_timing_mc_scaling(error)
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: mc_error
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(radius_type_static) :: radius_model
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
            call integrate_surface_marching_cubes(lsf, mol%xyz, area, volume, mc_error, &
                                                  target_spacing=mc_spacing)
            if (allocated(mc_error)) then
               call test_failed(error, mc_error%message)
               return
            end if
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
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      !> Polyalanine structures (increasing size)
      call new_context(ctx, verbosity=0)

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
            call new_cavity_iswig(cavity, ctx, nleb=nleb, &
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
            call cavity%get_gradient(cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
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

   !> Benchmark LSF evaluation cost per derivative order and per LSF type.
   !>
   !> Baseline for the code-generation refactor of the SvdW and CFC level set
   !> functions. Deliberately independent of the DROP projection: the evaluation
   !> points come from a fixed deterministic shell sample (see
   !> [[build_lsf_shell_points]]) rather than from a cavity, so the series stays
   !> comparable while the projection code changes underneath it.
   !>
   !> `prepare` and the accessors are timed separately, because the refactor
   !> deliberately moves work across exactly that boundary. `prepare` gets one
   !> sweep per `set_max_deriv` level; every accessor then gets its own sweep in
   !> which `prepare` (and, for the two accessors consuming its output,
   !> `f3_rr_rA`) runs outside the timer window, so an accessor row is
   !> never a difference of two much larger numbers.
   subroutine test_timing_lsf_accessors(error)
      type(error_type), allocatable, intent(out) :: error

      !> Number of polyalanine structures to benchmark
      integer, parameter :: n_struct = 25
      !> Fixed number of evaluation points per system (not per atom), so per-call
      !> costs are directly comparable across system sizes and the growth of the
      !> active-atom count is the scaling signal
      integer, parameter :: npts = 2000
      !> Reduced point budget for the accessors whose result array is quadratic
      !> in the atom count (tens of MB per call at the top of the series). The
      !> point order is quasi-uniform over the molecule, so a prefix of the
      !> sequence is still a fair sample
      integer, parameter :: npts_heavy = 200
      !> Point budget for CFC. Its pair term makes `prepare` up to three orders
      !> of magnitude more expensive than SvdW's, so the same 2000 points would
      !> dominate the runtime of the whole benchmark; CFC uses the leading
      !> `npts_cfc` points of the very same sequence
      integer, parameter :: npts_cfc = 250
      !> Minimum per-call time (s) for a data point to enter the power-law fit
      real(wp), parameter :: t_min = 1.0e-9_wp
      !> Highest derivative level supported by any LSF benchmarked here
      integer, parameter :: max_lvl = 4
      !> Largest accessor-row count of the LSFs benchmarked here
      integer, parameter :: max_acc = 16

      !> SvdW blend parameters, matching the other benchmarks in this file
      real(wp), parameter :: blend_k = 5.5_wp
      real(wp), parameter :: blend_2b = 0.0_wp
      real(wp), parameter :: blend_3b = 3.0_wp
      !> LSF screening threshold as the DROP cavity derives it
      !> (`screening_threshold = tolerance*0.1`, at the tolerance the other
      !> benchmarks in this file use)
      real(wp), parameter :: screen_thr = 1.0E-10_wp*0.1_wp

      !> SvdW accessor rows (orders 0-4)
      integer, parameter :: n_acc_svdw = 12
      character(len=24), parameter :: acc_labels_svdw(n_acc_svdw) = [character(len=24) :: &
                                      'f0', 'f012_r val', 'f012_r val+grad', &
                                      'f012_r val+grad+hess', 'f3_rrr', 'f3_rr_rA', &
                                      'f2_rArB', 'f3_r_rArB', 'f4_rrrr', 'f4_rrr_rA', &
                                      'f4_rr_rArB', 'normalized_f01_rA']
      !> `set_max_deriv` level each SvdW row is measured at
      integer, parameter :: acc_deriv_svdw(n_acc_svdw) = [0, 0, 1, 2, 3, 3, 3, 3, 4, 4, 4, 2]

      !> CFC accessor rows (orders 0-4, plus the direction-contracted families)
      integer, parameter :: n_acc_cfc = 16
      character(len=24), parameter :: acc_labels_cfc(n_acc_cfc) = [character(len=24) :: &
                                     'f0', 'f012_r val', 'f012_r val+grad', &
                                     'f012_r val+grad+hess', 'f3_rrr', 'f3_rr_rA', &
                                     'f4_rrrr', 'f4_rrr_rA', 'normalized_f01_rA', &
                                     'f2_rArB', 'f3_r_rArB', 'f4_rr_rArB', &
                                     'tangent_f2_rr', 'hvp_f1_rA', 'hvp_f2_r_rA', &
                                     'hvp_f3_rr_rA']
      !> `set_max_deriv` level each CFC row is measured at.
      !>
      !> CFC's `max_deriv` is the highest *total* order `prepare` provisions, so a
      !> row reading the one-nuclear-index family sits one level above the tensor
      !> rank it returns: `f3_rr_rA` reads `qn2_rr`, a total order 3.
      integer, parameter :: acc_deriv_cfc(n_acc_cfc) = &
                            [0, 0, 1, 2, 3, 3, 4, 4, 2, 1, 2, 3, 2, 1, 2, 3]

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      character(len=20) :: struct_names(n_struct)

      !> Collected benchmark data
      integer :: n_atoms_arr(n_struct)
      real(wp) :: prep_times(0:max_lvl, n_struct)
      real(wp) :: acc_times(max_acc, n_struct)
      integer :: acc_ncalls(max_acc), acc_reps(max_acc, n_struct)
      real(wp) :: act_mean(n_struct), checksums(n_struct)
      integer :: act_max(n_struct)

      integer :: ilsf, istruct, n_acc, n_lvl
      integer(kind=8) :: clock0, clock1, clock_rate
      real(wp) :: elapsed

      struct_names = [character(len=20) :: &
                      'polyala_04', 'polyala_08', 'polyala_12', 'polyala_16', 'polyala_20', &
                      'polyala_24', 'polyala_28', 'polyala_32', 'polyala_36', 'polyala_40', &
                      'polyala_44', 'polyala_48', 'polyala_52', 'polyala_56', 'polyala_60', &
                      'polyala_64', 'polyala_68', 'polyala_72', 'polyala_76', 'polyala_80', &
                      'polyala_84', 'polyala_88', 'polyala_92', 'polyala_96', 'polyala_100']

      write (*, '(a)') ''
      write (*, '(a)') '=================================================================='
      write (*, '(a)') 'Benchmark: LSF cost per derivative order (prepare vs accessors)'
      write (*, '(a,i0,a,i0,a,i0)') &
         'points/system: SvdW=', npts, '  SvdW heavy accessors=', npts_heavy, &
         '  CFC=', npts_cfc
      write (*, '(a,es9.2)') 'screening_threshold=', screen_thr
      write (*, '(a,f5.2,a,f5.2,a,f5.2)') &
         'SvdW: k=', blend_k, '  b2=', blend_2b, '  b3=', blend_3b
      write (*, '(a)') 'CFC : Diedenhofen-Klamt defaults'
      write (*, '(a)') 'Points: deterministic Fibonacci shells, no cavity, no random_number'
      write (*, '(a)') '=================================================================='

      do ilsf = 1, 2
         if (ilsf == 1) then
            n_acc = n_acc_svdw
            n_lvl = 4
         else
            n_acc = n_acc_cfc
            n_lvl = 4
         end if
         prep_times = 0.0_wp
         acc_times = 0.0_wp
         acc_ncalls = 0

         call system_clock(clock0, clock_rate)

         do istruct = 1, n_struct
            call get_structure(mol, 'POLYALANINE', trim(struct_names(istruct)))
            n_atoms_arr(istruct) = mol%nat

            if (allocated(radii)) deallocate (radii)
            call fill_cpcm_radii(mol, radii, error)
            if (allocated(error)) return
            radii = radii*aatoau

            if (allocated(points)) deallocate (points)
            allocate (points(ndim, npts))
            call build_lsf_shell_points(mol, radii, points)

            if (ilsf == 1) then
               call bench_lsf_svdw(mol, radii, points, npts_heavy, screen_thr, &
                                   blend_k, blend_2b, blend_3b, &
                                   prep_times(0:n_lvl, istruct), &
                                   acc_times(1:n_acc, istruct), acc_ncalls(1:n_acc), &
                                   acc_reps(1:n_acc, istruct), act_mean(istruct), &
                                   act_max(istruct), &
                                   checksums(istruct), error)
            else
               call bench_lsf_cfc(mol, radii, points, npts_cfc, screen_thr, &
                                  prep_times(0:n_lvl, istruct), &
                                  acc_times(1:n_acc, istruct), acc_ncalls(1:n_acc), &
                                  acc_reps(1:n_acc, istruct), act_mean(istruct), &
                                  act_max(istruct), &
                                  checksums(istruct), error)
            end if
            if (allocated(error)) return

            write (*, '(2x,a14,a,i5,a,f7.1,a,i5)') &
               trim(struct_names(istruct)), '  N_at=', n_atoms_arr(istruct), &
               '  act_mean=', act_mean(istruct), '  act_max=', act_max(istruct)
         end do

         call system_clock(clock1)
         elapsed = real(clock1 - clock0, wp)/real(clock_rate, wp)

         if (ilsf == 1) then
            call report_lsf_timings('SvdW', n_struct, struct_names, n_atoms_arr, &
                                    n_lvl, prep_times(0:n_lvl, :), n_acc, &
                                    acc_labels_svdw, acc_deriv_svdw, &
                                    acc_times(1:n_acc, :), acc_ncalls(1:n_acc), &
                                    acc_reps(1:n_acc, :), act_mean, act_max, &
                                    checksums, t_min, elapsed)
         else
            call report_lsf_timings('CFC', n_struct, struct_names, n_atoms_arr, &
                                    n_lvl, prep_times(0:n_lvl, :), n_acc, &
                                    acc_labels_cfc, acc_deriv_cfc, &
                                    acc_times(1:n_acc, :), acc_ncalls(1:n_acc), &
                                    acc_reps(1:n_acc, :), act_mean, act_max, &
                                    checksums, t_min, elapsed)
         end if
      end do

      write (*, '(a)') ''

   end subroutine test_timing_lsf_accessors

   !> Deterministic shell sample of evaluation points around a molecule.
   !>
   !> Point `p` is assigned to atom `a = 1 + mod((p-1)*stride, nat)` with a
   !> stride prime larger than any structure benchmarked here, so consecutive
   !> points hop across the molecule and any *prefix* of the sequence is still
   !> spread over all atoms (the O(ncenters^2) accessors are timed on a prefix).
   !> The point is placed at `r_a + (R_a + delta)*u_p`, with `u_p` the p-th
   !> direction of a fixed Fibonacci sphere and `delta` cycling through
   !> {0, 0.25, 0.5} Bohr, so points sit near where a cavity surface would lie.
   !>
   !> No `random_number` anywhere: the sample depends only on the geometry.
   !>
   !> @param[in]  mol    Molecular structure
   !> @param[in]  radii  Atomic radii in Bohr [nat]
   !> @param[out] points Output points [ndim, npts]
   subroutine build_lsf_shell_points(mol, radii, points)
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)
      real(wp), intent(out) :: points(:, :)

      !> Prime larger than any structure in the series, hence coprime to `nat`
      integer, parameter :: stride = 7919
      !> Radial offsets from the atomic sphere (Bohr), cycled over the points
      real(wp), parameter :: shell_offsets(3) = [0.0_wp, 0.25_wp, 0.5_wp]

      real(wp) :: golden, theta, uz, sr, delta
      integer :: npts, ipt, iat

      npts = size(points, 2)
      golden = acos(-1.0_wp)*(3.0_wp - sqrt(5.0_wp))

      do ipt = 1, npts
         iat = 1 + mod((ipt - 1)*stride, mol%nat)
         uz = 1.0_wp - 2.0_wp*(real(ipt, wp) - 0.5_wp)/real(npts, wp)
         sr = sqrt(max(0.0_wp, 1.0_wp - uz*uz))
         theta = golden*real(ipt - 1, wp)
         delta = shell_offsets(1 + mod(ipt - 1, size(shell_offsets)))
         points(:, ipt) = mol%xyz(:, iat) &
                          + (radii(iat) + delta)*[sr*cos(theta), sr*sin(theta), uz]
      end do
   end subroutine build_lsf_shell_points

   !> Time SvdW `prepare` and every SvdW accessor on a fixed point set.
   !>
   !> @param[in]  mol         Molecular structure
   !> @param[in]  radii       Atomic radii in Bohr [nat]
   !> @param[in]  points      Evaluation points [ndim, npts]
   !> @param[in]  npts_heavy  Point budget for O(ncenters^2) accessors
   !> @param[in]  thr         LSF screening threshold
   !> @param[in]  blend_k     SvdW blending sharpness
   !> @param[in]  blend_2b    SvdW two-body weight
   !> @param[in]  blend_3b    SvdW three-body weight
   !> @param[out] prep_time   Seconds per `prepare` call, indexed 0:max_deriv
   !> @param[out] acc_time    Seconds per accessor call (prepare excluded)
   !> @param[out] acc_ncalls  Points actually used per accessor row
   !> @param[out] acc_reps    Accessor calls per point used per row
   !> @param[out] act_mean    Mean active-atom count over the point set
   !> @param[out] act_max     Maximum active-atom count over the point set
   !> @param[out] chk         Value checksum, so a refactor can be checked
   !> @param[out] error       Test failure
   subroutine bench_lsf_svdw(mol, radii, points, npts_heavy, thr, &
                             blend_k, blend_2b, blend_3b, &
                             prep_time, acc_time, acc_ncalls, acc_reps, &
                             act_mean, act_max, chk, error)
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)
      real(wp), intent(in), target :: points(:, :)
      integer, intent(in) :: npts_heavy
      real(wp), intent(in) :: thr, blend_k, blend_2b, blend_3b
      real(wp), intent(out) :: prep_time(0:)
      real(wp), intent(out) :: acc_time(:)
      integer, intent(out) :: acc_ncalls(:), acc_reps(:)
      real(wp), intent(out) :: act_mean
      integer, intent(out) :: act_max
      real(wp), intent(out) :: chk
      type(error_type), allocatable, intent(out) :: error

      !> Accessor rows whose baseline pass has to run `f3_rr_rA`
      !> first, because they take its output as `intent(in)`
      logical, parameter :: needs_rr_rA(12) = [.false., .false., .false., .false., &
                                               .false., .false., .false., .true., &
                                               .false., .false., .true., .false.]
      !> Accessor rows whose result array is quadratic in the atom count
      !> (`f2_rArB` in ncenters, `f3_r_rArB` and `f4_rr_rArB` in the active
      !> count). At the cavity screening threshold the active count reaches a
      !> few hundred atoms, so these allocate tens of MB per call
      logical, parameter :: is_heavy(12) = [.false., .false., .false., .false., &
                                            .false., .false., .true., .true., &
                                            .false., .false., .true., .false.]
      integer, parameter :: acc_deriv(12) = [0, 0, 1, 2, 3, 3, 3, 3, 4, 4, 4, 2]

      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(mctc_error), allocatable :: lsf_err
      real(wp) :: val, lsf0, lsf1_r(ndim), lsf2_rr(ndim, ndim)
      real(wp) :: lsf3_rrr(ndim, ndim, ndim), res_rrrr(ndim, ndim, ndim, ndim)
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :), lsf3_rr_rA(:, :, :, :)
      real(wp), allocatable :: res_rArB(:, :, :, :), res_r_rArB(:, :, :, :, :)
      real(wp), allocatable :: res_rrr_rA(:, :, :, :, :)
      real(wp), allocatable :: res_rr_rArB(:, :, :, :, :, :), deriv_rA(:, :)
      real :: c0, c1
      integer :: npts, ipt, ilvl, iacc, nc, nrep, n_act, act_sum
      !> Set for the first accessor call at each point, so that repeated calls
      !> do not multiply the checksum
      logical :: chk_on

      npts = size(points, 2)
      chk = 0.0_wp

      call lsf%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
      lsf%screening_threshold = thr
      call lsf%update(mol, radii)


      !> ================= prepare, one pass per derivative level =================
      do ilvl = 0, ubound(prep_time, 1)
         call lsf%set_max_deriv(ilvl)
         call cpu_time(c0)
         do ipt = 1, npts
            call lsf%prepare(points(:, ipt), lsf_err)
            if (allocated(lsf_err)) exit
         end do
         call cpu_time(c1)
         if (allocated(lsf_err)) then
            call test_failed(error, "SvdW prepare failed: "//lsf_err%message)
            return
         end if
         prep_time(ilvl) = real(c1 - c0, wp)/real(npts, wp)
      end do

      !> ======================== Active-atom statistics =========================
      call lsf%set_max_deriv(2)
      act_sum = 0
      act_max = 0
      do ipt = 1, npts
         call lsf%prepare(points(:, ipt), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "SvdW prepare failed: "//lsf_err%message)
            return
         end if
         n_act = lsf%active_count()
         act_sum = act_sum + n_act
         act_max = max(act_max, n_act)
      end do
      act_mean = real(act_sum, wp)/real(npts, wp)

      !> ==================== Caller-owned result buffers ====================
      !> Every nuclear output is active-indexed, so `act_max` (measured just
      !> above over the very same point set) is the exact capacity needed --
      !> no accessor allocates anything of its own any more
      allocate (lsf1_rA(ndim, act_max), source=0.0_wp)
      allocate (lsf2_r_rA(ndim, ndim, act_max), source=0.0_wp)
      allocate (lsf3_rr_rA(ndim, ndim, ndim, act_max), source=0.0_wp)
      allocate (res_rArB(ndim, act_max, ndim, act_max), source=0.0_wp)
      allocate (res_r_rArB(ndim, ndim, act_max, ndim, act_max), source=0.0_wp)
      allocate (res_rrr_rA(ndim, ndim, ndim, ndim, act_max), source=0.0_wp)
      allocate (res_rr_rArB(ndim, ndim, ndim, act_max, ndim, act_max), source=0.0_wp)
      allocate (deriv_rA(ndim, act_max), source=0.0_wp)

      !> ======================= Accessors, one row at a time =======================
      do iacc = 1, size(acc_time)
         nc = npts
         if (is_heavy(iacc)) nc = min(npts, npts_heavy)
         acc_ncalls(iacc) = nc
         call lsf%set_max_deriv(acc_deriv(iacc))

         nrep = calibrate_reps(iacc)
         acc_reps(iacc) = nrep
         call time_accessor(iacc, nc, nrep, acc_time(iacc))
         if (allocated(error)) return
      end do

   contains

      !> Time accessor row `iacc` over the first `nc` points.
      !>
      !> `prepare` (and, for the two accessors consuming it, the
      !> `f3_rr_rA` that produces their input) runs outside the timer
      !> window, so the row measures the accessor alone rather than a difference
      !> of two much larger numbers - the CFC accessors are five orders of
      !> magnitude cheaper than their `prepare`, which no subtraction survives.
      !>
      !> @param[in]  iacc     Accessor row
      !> @param[in]  nc       Points to sweep
      !> @param[in]  nrep     Accessor calls per prepared point
      !> @param[out] per_call Seconds per accessor call
      subroutine time_accessor(iacc, nc, nrep, per_call)
         integer, intent(in) :: iacc, nc, nrep
         real(wp), intent(out) :: per_call
         integer :: ip, ir
         real :: p0, p1
         real(wp) :: total

         total = 0.0_wp
         do ip = 1, nc
            call lsf%prepare(points(:, ip), lsf_err)
            if (allocated(lsf_err)) exit
            if (needs_rr_rA(iacc)) &
               call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
            chk_on = .true.
            call cpu_time(p0)
            do ir = 1, nrep
               call call_accessor(iacc)
               chk_on = .false.
            end do
            call cpu_time(p1)
            total = total + real(p1 - p0, wp)
         end do
         per_call = total/real(nc, wp)/real(nrep, wp)
         if (allocated(lsf_err)) &
            call test_failed(error, "SvdW prepare failed: "//lsf_err%message)
      end subroutine time_accessor

      !> Accessor calls per prepared point, so that one timer window sits well
      !> above the clock resolution. Calibrated on the first point by growing
      !> the count until the window is long enough; the repeated calls hit a
      !> warm cache, which is also how the projector uses these accessors
      !> (several of them back-to-back after a single `prepare`)
      !>
      !> @param[in] iacc  Accessor row
      !> @returns   nrep  Calls per prepared point
      integer function calibrate_reps(iacc) result(nrep)
         integer, intent(in) :: iacc

         !> Target length of one timer window (s)
         real(wp), parameter :: window = 5.0e-5_wp
         !> Cap, so a near-free accessor cannot spin here forever
         integer, parameter :: max_reps = 4096
         integer :: ir
         real :: p0, p1

         chk_on = .false.
         call lsf%prepare(points(:, 1), lsf_err)
         if (allocated(lsf_err)) then
            nrep = 1
            return
         end if
         if (needs_rr_rA(iacc)) &
            call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)

         nrep = 1
         do
            call cpu_time(p0)
            do ir = 1, nrep
               call call_accessor(iacc)
            end do
            call cpu_time(p1)
            if (real(p1 - p0, wp) >= window .or. nrep >= max_reps) exit
            nrep = min(max_reps, nrep*4)
         end do
      end function calibrate_reps

      !> Invoke accessor row `iacc` once at the prepared point. On the first
      !> call at each point one element of the result feeds the checksum: enough
      !> to pin the values down across a refactor without adding a traversal of
      !> the (large) result arrays
      subroutine call_accessor(iacc)
         integer, intent(in) :: iacc

         select case (iacc)
         case (1)
            call lsf%f0(val)
            if (chk_on) chk = chk + val
         case (2)
            call lsf%f012_r(lsf0=lsf0)
            if (chk_on) chk = chk + lsf0
         case (3)
            call lsf%f012_r(lsf0=lsf0, lsf1_r=lsf1_r)
            if (chk_on) chk = chk + lsf1_r(1)
         case (4)
            call lsf%f012_r(lsf0=lsf0, lsf1_r=lsf1_r, lsf2_rr=lsf2_rr)
            if (chk_on) chk = chk + lsf2_rr(1, 1)
         case (5)
            call lsf%f3_rrr(lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
            if (chk_on) chk = chk + lsf3_rrr(1, 1, 1)
         case (6)
            call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
            if (chk_on) chk = chk + lsf3_rr_rA(1, 1, 1, 1)
         case (7)
            call lsf%f2_rArB(res_rArB)
            if (chk_on .and. size(res_rArB) > 0) chk = chk + res_rArB(1, 1, 1, 1)
         case (8)
            call lsf%f3_r_rArB(res_r_rArB)
            if (chk_on .and. size(res_r_rArB) > 0) chk = chk + res_r_rArB(1, 1, 1, 1, 1)
         case (9)
            call lsf%f4_rrrr(res_rrrr)
            if (chk_on) chk = chk + res_rrrr(1, 1, 1, 1)
         case (10)
            call lsf%f4_rrr_rA(res_rrr_rA)
            if (chk_on .and. size(res_rrr_rA) > 0) chk = chk + res_rrr_rA(1, 1, 1, 1, 1)
         case (11)
            call lsf%f4_rr_rArB(res_rr_rArB)
            if (chk_on .and. size(res_rr_rArB) > 0) chk = chk + res_rr_rArB(1, 1, 1, 1, 1, 1)
         case (12)
            call lsf%normalized_f01_rA(val, deriv_rA=deriv_rA)
            if (chk_on) chk = chk + val + deriv_rA(1, 1)
         end select
      end subroutine call_accessor

   end subroutine bench_lsf_svdw

   !> Time CFC `prepare` and every CFC accessor on a fixed point set.
   !>
   !> Same protocol as [[bench_lsf_svdw]]; the CFC LSF stops at order 3 and has
   !> no nuclear-pair accessors, so no row needs its own point budget - but the
   !> pair term makes every `prepare` expensive, hence the smaller `npts_use`.
   !>
   !> @param[in]  mol         Molecular structure
   !> @param[in]  radii       Atomic radii in Bohr [nat]
   !> @param[in]  points      Evaluation points [ndim, npts]
   !> @param[in]  npts_use    Leading points of `points` to evaluate on
   !> @param[in]  thr         LSF screening threshold
   !> @param[out] prep_time   Seconds per `prepare` call, indexed 0:max_deriv
   !> @param[out] acc_time    Seconds per accessor call (prepare excluded)
   !> @param[out] acc_ncalls  Points actually used per accessor row
   !> @param[out] acc_reps    Accessor calls per point used per row
   !> @param[out] act_mean    Mean active-atom count over the point set
   !> @param[out] act_max     Maximum active-atom count over the point set
   !> @param[out] chk         Value checksum, so a refactor can be checked
   !> @param[out] error       Test failure
   subroutine bench_lsf_cfc(mol, radii, points, npts_use, thr, &
                            prep_time, acc_time, acc_ncalls, acc_reps, &
                            act_mean, act_max, chk, error)
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)
      real(wp), intent(in), target :: points(:, :)
      integer, intent(in) :: npts_use
      real(wp), intent(in) :: thr
      real(wp), intent(out) :: prep_time(0:)
      real(wp), intent(out) :: acc_time(:)
      integer, intent(out) :: acc_ncalls(:), acc_reps(:)
      real(wp), intent(out) :: act_mean
      integer, intent(out) :: act_max
      real(wp), intent(out) :: chk
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: acc_deriv(16) = &
                            [0, 0, 1, 2, 3, 3, 4, 4, 2, 1, 2, 3, 2, 1, 2, 3]
      !> Per-row point cap; 0 means "use the full `npts_use` prefix".
      !>
      !> Two costs drive these, and both are quadratic in the active-atom count.
      !> First, `time_accessor` re-`prepare`s at the row's own `set_max_deriv`
      !> level for every point, and a CFC `prepare` at level 4 accumulates the
      !> two largest pair branches in the module (1101 + 2484 CSE temporaries),
      !> so a row at that level costs tens of milliseconds per point before the
      !> accessor is even called. Second, the two-nucleus and direction-
      !> contracted rows run their own O(n_active**2) sweep on top of that --
      !> the `f*_rArB` rows twice over, once for the diagonal blocks and once
      !> for the ordered pairs.
      !>
      !> The caps keep every row near or below a second per structure. The point
      !> sequence is quasi-uniform over the molecule, so a shorter prefix stays
      !> a fair sample; `points` in the report names the budget each row used.
      integer, parameter :: acc_npts_cap(16) = &
                            [0, 0, 0, 0, 50, 50, 10, 10, 0, 10, 4, 2, 10, 10, 10, 4]

      !> Point budget of the `prepare` pass at each derivative level.
      !>
      !> Level 4 accumulates the order-4 spatial ladder and the order-3 nuclear
      !> one, i.e. the two largest pair branches in the module (1101 + 2484 CSE
      !> temporaries), so at the top of the series a single point costs tens of
      !> milliseconds. Same prefix argument as above.
      integer, parameter :: prep_npts_cap(0:4) = [0, 0, 0, 50, 10]

      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(mctc_error), allocatable :: lsf_err
      real(wp) :: val, lsf0, lsf1_r(ndim), lsf2_rr(ndim, ndim)
      real(wp) :: lsf3_rrr(ndim, ndim, ndim), lsf4_rrrr(ndim, ndim, ndim, ndim)
      real(wp) :: norm0, tng2(ndim, ndim)
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :), lsf3_rr_rA(:, :, :, :)
      real(wp), allocatable :: lsf4_rrr_rA(:, :, :, :, :)
      real(wp), allocatable :: lsf2_rArB(:, :, :, :), lsf3_r_rArB(:, :, :, :, :)
      real(wp), allocatable :: lsf4_rr_rArB(:, :, :, :, :, :)
      real(wp), allocatable :: hvp1(:, :), hvp2(:, :, :), hvp3(:, :, :, :)
      real(wp), allocatable :: vdir(:, :)
      real :: c0, c1
      integer :: npts, ipt, ilvl, iacc, nrep, n_act, act_sum, iat, npts_row
      !> Set for the first accessor call at each point, so that repeated calls
      !> do not multiply the checksum
      logical :: chk_on

      npts = min(npts_use, size(points, 2))
      chk = 0.0_wp

      call lsf%new()
      lsf%screening_threshold = thr
      call lsf%update(mol, radii)

      !> ================= prepare, one pass per derivative level =================
      do ilvl = 0, ubound(prep_time, 1)
         npts_row = npts
         if (prep_npts_cap(ilvl) > 0) npts_row = min(npts, prep_npts_cap(ilvl))
         call lsf%set_max_deriv(ilvl)
         call cpu_time(c0)
         do ipt = 1, npts_row
            call lsf%prepare(points(:, ipt), lsf_err)
            if (allocated(lsf_err)) exit
         end do
         call cpu_time(c1)
         if (allocated(lsf_err)) then
            call test_failed(error, "CFC prepare failed: "//lsf_err%message)
            return
         end if
         prep_time(ilvl) = real(c1 - c0, wp)/real(npts_row, wp)
      end do

      !> ======================== Active-atom statistics =========================
      call lsf%set_max_deriv(2)
      act_sum = 0
      act_max = 0
      do ipt = 1, npts
         call lsf%prepare(points(:, ipt), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "CFC prepare failed: "//lsf_err%message)
            return
         end if
         n_act = lsf%active_count()
         act_sum = act_sum + n_act
         act_max = max(act_max, n_act)
      end do
      act_mean = real(act_sum, wp)/real(npts, wp)

      !> The nuclear outputs are active-indexed and caller-owned; `act_max`
      !> (measured just above on this point set) is the capacity needed
      allocate (lsf1_rA(ndim, act_max), source=0.0_wp)
      allocate (lsf2_r_rA(ndim, ndim, act_max), source=0.0_wp)
      allocate (lsf3_rr_rA(ndim, ndim, ndim, act_max), source=0.0_wp)
      allocate (lsf4_rrr_rA(ndim, ndim, ndim, ndim, act_max), source=0.0_wp)
      allocate (lsf2_rArB(ndim, act_max, ndim, act_max), source=0.0_wp)
      allocate (lsf3_r_rArB(ndim, ndim, act_max, ndim, act_max), source=0.0_wp)
      allocate (lsf4_rr_rArB(ndim, ndim, ndim, act_max, ndim, act_max), source=0.0_wp)
      allocate (hvp1(ndim, act_max), source=0.0_wp)
      allocate (hvp2(ndim, ndim, act_max), source=0.0_wp)
      allocate (hvp3(ndim, ndim, ndim, act_max), source=0.0_wp)

      !> Deterministic displacement field for the contracted rows
      allocate (vdir(ndim, mol%nat))
      do iat = 1, mol%nat
         vdir(1, iat) = 0.31_wp*real(iat, wp) - 0.7_wp
         vdir(2, iat) = -0.17_wp*real(iat, wp) + 0.4_wp
         vdir(3, iat) = 0.23_wp*real(mod(iat, 3) + 1, wp)
      end do

      !> ======================= Accessors, one row at a time =======================
      do iacc = 1, size(acc_time)
         npts_row = npts
         if (acc_npts_cap(iacc) > 0) npts_row = min(npts, acc_npts_cap(iacc))
         acc_ncalls(iacc) = npts_row
         call lsf%set_max_deriv(acc_deriv(iacc))

         nrep = calibrate_reps(iacc)
         acc_reps(iacc) = nrep
         call time_accessor(iacc, npts_row, nrep, acc_time(iacc))
         if (allocated(error)) return
      end do

   contains

      !> Time accessor row `iacc` over the first `nc` points, with `prepare`
      !> outside the timer window. CFC's pair term makes `prepare` five orders
      !> of magnitude more expensive than the accessors that read its caches, a
      !> gap no difference of two passes could resolve
      !>
      !> @param[in]  iacc     Accessor row
      !> @param[in]  nc       Points to sweep
      !> @param[in]  nrep     Accessor calls per prepared point
      !> @param[out] per_call Seconds per accessor call
      subroutine time_accessor(iacc, nc, nrep, per_call)
         integer, intent(in) :: iacc, nc, nrep
         real(wp), intent(out) :: per_call
         integer :: ip, ir
         real :: p0, p1
         real(wp) :: total

         total = 0.0_wp
         do ip = 1, nc
            call lsf%prepare(points(:, ip), lsf_err)
            if (allocated(lsf_err)) exit
            chk_on = .true.
            call cpu_time(p0)
            do ir = 1, nrep
               call call_accessor(iacc)
               chk_on = .false.
            end do
            call cpu_time(p1)
            total = total + real(p1 - p0, wp)
         end do
         per_call = total/real(nc, wp)/real(nrep, wp)
         if (allocated(lsf_err)) &
            call test_failed(error, "CFC prepare failed: "//lsf_err%message)
      end subroutine time_accessor

      !> Accessor calls per prepared point, so that one timer window sits well
      !> above the clock resolution; calibrated on the first point
      !>
      !> @param[in] iacc  Accessor row
      !> @returns   nrep  Calls per prepared point
      integer function calibrate_reps(iacc) result(nrep)
         integer, intent(in) :: iacc

         !> Target length of one timer window (s)
         real(wp), parameter :: window = 5.0e-5_wp
         !> Cap, so a near-free accessor cannot spin here forever
         integer, parameter :: max_reps = 4096
         integer :: ir
         real :: p0, p1

         chk_on = .false.
         call lsf%prepare(points(:, 1), lsf_err)
         if (allocated(lsf_err)) then
            nrep = 1
            return
         end if

         nrep = 1
         do
            call cpu_time(p0)
            do ir = 1, nrep
               call call_accessor(iacc)
            end do
            call cpu_time(p1)
            if (real(p1 - p0, wp) >= window .or. nrep >= max_reps) exit
            nrep = min(max_reps, nrep*4)
         end do
      end function calibrate_reps

      !> Invoke accessor row `iacc` once at the prepared point
      subroutine call_accessor(iacc)
         integer, intent(in) :: iacc

         select case (iacc)
         case (1)
            call lsf%f0(val)
            if (chk_on) chk = chk + val
         case (2)
            call lsf%f012_r(lsf0=lsf0)
            if (chk_on) chk = chk + lsf0
         case (3)
            call lsf%f012_r(lsf0=lsf0, lsf1_r=lsf1_r)
            if (chk_on) chk = chk + lsf1_r(1)
         case (4)
            call lsf%f012_r(lsf0=lsf0, lsf1_r=lsf1_r, lsf2_rr=lsf2_rr)
            if (chk_on) chk = chk + lsf2_rr(1, 1)
         case (5)
            call lsf%f3_rrr(lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
            if (chk_on) chk = chk + lsf3_rrr(1, 1, 1)
         case (6)
            call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
            if (chk_on) chk = chk + lsf3_rr_rA(1, 1, 1, 1)
         case (7)
            call lsf%f4_rrrr(lsf4_rrrr)
            if (chk_on) chk = chk + lsf4_rrrr(1, 1, 1, 1)
         case (8)
            call lsf%f4_rrr_rA(lsf4_rrr_rA)
            if (chk_on) chk = chk + lsf4_rrr_rA(1, 1, 1, 1, 1)
         case (9)
            call lsf%normalized_f01_rA(norm0, deriv_rA=lsf1_rA)
            if (chk_on) chk = chk + norm0
         case (10)
            call lsf%f2_rArB(lsf2_rArB)
            if (chk_on) chk = chk + lsf2_rArB(1, 1, 1, 1)
         case (11)
            call lsf%f3_r_rArB(lsf3_r_rArB)
            if (chk_on) chk = chk + lsf3_r_rArB(1, 1, 1, 1, 1)
         case (12)
            call lsf%f4_rr_rArB(lsf4_rr_rArB)
            if (chk_on) chk = chk + lsf4_rr_rArB(1, 1, 1, 1, 1, 1)
         case (13)
            call lsf%tangent_f2_rr(vdir, tng2)
            if (chk_on) chk = chk + tng2(1, 1)
         case (14)
            call lsf%hvp_f1_rA(vdir, hvp1)
            if (chk_on) chk = chk + hvp1(1, 1)
         case (15)
            call lsf%hvp_f2_r_rA(vdir, hvp2)
            if (chk_on) chk = chk + hvp2(1, 1, 1)
         case (16)
            call lsf%hvp_f3_rr_rA(vdir, hvp3)
            if (chk_on) chk = chk + hvp3(1, 1, 1, 1)
         end select
      end subroutine call_accessor

   end subroutine bench_lsf_cfc

   !> Print the per-system, per-accessor and scaling tables for one LSF.
   !>
   !> @param[in] title        LSF name used in the table headers
   !> @param[in] n_struct     Number of structures
   !> @param[in] struct_names Structure labels
   !> @param[in] n_atoms      Atom counts per structure
   !> @param[in] n_lvl        Highest derivative level measured
   !> @param[in] prep_times   Seconds per prepare call [0:n_lvl, n_struct]
   !> @param[in] n_acc        Number of accessor rows
   !> @param[in] acc_labels   Accessor row labels
   !> @param[in] acc_deriv    `set_max_deriv` level of each accessor row
   !> @param[in] acc_times    Seconds per accessor call [n_acc, n_struct]
   !> @param[in] acc_ncalls   Points used per accessor row
   !> @param[in] acc_reps     Accessor calls per point [n_acc, n_struct]
   !> @param[in] act_mean     Mean active-atom count per structure
   !> @param[in] act_max      Maximum active-atom count per structure
   !> @param[in] checksums    Value checksum per structure
   !> @param[in] t_min        Minimum per-call time entering the fit
   !> @param[in] elapsed      Wall time spent benchmarking this LSF (s)
   subroutine report_lsf_timings(title, n_struct, struct_names, n_atoms, n_lvl, &
                                 prep_times, n_acc, acc_labels, acc_deriv, &
                                 acc_times, acc_ncalls, acc_reps, act_mean, act_max, &
                                 checksums, t_min, elapsed)
      character(len=*), intent(in) :: title
      integer, intent(in) :: n_struct
      character(len=*), intent(in) :: struct_names(:)
      integer, intent(in) :: n_atoms(:)
      integer, intent(in) :: n_lvl
      real(wp), intent(in) :: prep_times(0:, :)
      integer, intent(in) :: n_acc
      character(len=*), intent(in) :: acc_labels(:)
      integer, intent(in) :: acc_deriv(:)
      real(wp), intent(in) :: acc_times(:, :)
      integer, intent(in) :: acc_ncalls(:)
      integer, intent(in) :: acc_reps(:, :)
      real(wp), intent(in) :: act_mean(:)
      integer, intent(in) :: act_max(:)
      real(wp), intent(in) :: checksums(:)
      real(wp), intent(in) :: t_min, elapsed

      !> Structures per column block of the accessor table
      integer, parameter :: cols = 7

      real(wp) :: real_n(n_struct), raw_t(n_struct), row(n_struct)
      real(wp) :: exponent, r_sq, coeff_a
      integer :: n_valid, istruct, ilvl, iacc, iblock, ifirst, ilast
      character(len=12) :: hdr

      !> ==================== Per-system summary (prepare) ====================
      write (*, '(a)') ''
      write (*, '(a,a,a,f8.1,a)') '=== ', title, ': per-system summary   (benchmarked in ', &
         elapsed, ' s wall) ==='
      write (*, '(a14, a7, a10, a9)', advance='no') 'Structure', 'N_at', 'act_mean', 'act_max'
      do ilvl = 0, n_lvl
         write (hdr, '(a,i0,a)') 'prep_d', ilvl, '(us)'
         write (*, '(a13)', advance='no') hdr
      end do
      write (*, '(a18)') 'checksum'

      write (*, '(a14, a7, a10, a9)', advance='no') repeat('-', 13), repeat('-', 6), &
         repeat('-', 9), repeat('-', 8)
      do ilvl = 0, n_lvl
         write (*, '(a13)', advance='no') repeat('-', 12)
      end do
      write (*, '(a18)') repeat('-', 17)

      do istruct = 1, n_struct
         write (*, '(a14, i7, f10.2, i9)', advance='no') &
            trim(struct_names(istruct)), n_atoms(istruct), &
            act_mean(istruct), act_max(istruct)
         do ilvl = 0, n_lvl
            write (*, '(f13.4)', advance='no') prep_times(ilvl, istruct)*1.0e6_wp
         end do
         write (*, '(es18.9)') checksums(istruct)
      end do

      !> ================== Per-accessor cost, blocked columns ==================
      write (*, '(a)') ''
      write (*, '(a,a,a)') '=== ', title, ': accessor cost in us/call (prepare excluded) ==='
      do iblock = 1, (n_struct + cols - 1)/cols
         ifirst = (iblock - 1)*cols + 1
         ilast = min(n_struct, ifirst + cols - 1)

         write (*, '(a24, a4, a8)', advance='no') 'N_at ->', 'd', 'points'
         do istruct = ifirst, ilast
            write (*, '(i12)', advance='no') n_atoms(istruct)
         end do
         write (*, *)
         write (*, '(a24, a4, a8)', advance='no') repeat('-', 23), repeat('-', 3), &
            repeat('-', 7)
         do istruct = ifirst, ilast
            write (*, '(a12)', advance='no') repeat('-', 11)
         end do
         write (*, *)

         do iacc = 1, n_acc
            write (*, '(a24, i4, i8)', advance='no') &
               acc_labels(iacc), acc_deriv(iacc), acc_ncalls(iacc)
            do istruct = ifirst, ilast
               write (*, '(f12.4)', advance='no') acc_times(iacc, istruct)*1.0e6_wp
            end do
            write (*, *)
         end do
         write (*, '(a)') ''
      end do

      !> ==================== Repeat counts, blocked columns ====================
      !> Accessor calls timed per prepared point behind each timing above; a
      !> row at 1 was measured call-for-call, a higher count means the accessor
      !> is too cheap to resolve in a single timer window
      write (*, '(a,a,a)') '=== ', title, ': accessor calls per prepared point (repeat count) ==='
      do iblock = 1, (n_struct + cols - 1)/cols
         ifirst = (iblock - 1)*cols + 1
         ilast = min(n_struct, ifirst + cols - 1)

         write (*, '(a24, a4)', advance='no') 'N_at ->', 'd'
         do istruct = ifirst, ilast
            write (*, '(i12)', advance='no') n_atoms(istruct)
         end do
         write (*, *)
         write (*, '(a24, a4)', advance='no') repeat('-', 23), repeat('-', 3)
         do istruct = ifirst, ilast
            write (*, '(a12)', advance='no') repeat('-', 11)
         end do
         write (*, *)

         do iacc = 1, n_acc
            write (*, '(a24, i4)', advance='no') acc_labels(iacc), acc_deriv(iacc)
            do istruct = ifirst, ilast
               write (*, '(i12)', advance='no') acc_reps(iacc, istruct)
            end do
            write (*, *)
         end do
         write (*, '(a)') ''
      end do

      !> ========================== Scaling exponents ==========================
      write (*, '(a)') ''
      write (*, '(a,a,a)') '=== ', title, ': per-call scaling fit  t(N) = A * N^X ==='
      write (*, '(a24, a4, a10, a10, a12, a12, a12)') &
         'Row', 'd', 'X', 'R^2', 'A', 'us@N_min', 'us@N_max'
      write (*, '(a24, a4, a10, a10, a12, a12, a12)') repeat('-', 23), repeat('-', 3), &
         repeat('-', 9), repeat('-', 9), repeat('-', 11), repeat('-', 11), repeat('-', 11)

      do ilvl = 0, n_lvl
         write (hdr, '(a,i0,a)') 'prepare d', ilvl
         row = prep_times(ilvl, 1:n_struct)
         call collect_valid_points(n_struct, n_atoms, row, t_min, real_n, raw_t, n_valid)
         if (n_valid >= 4) then
            call fit_power_law(n_valid, real_n, raw_t, exponent, r_sq, coeff_a)
            write (*, '(a24, i4, f10.3, f10.4, es12.3, f12.4, f12.4)') &
               hdr, ilvl, exponent, r_sq, coeff_a, &
               row(1)*1.0e6_wp, row(n_struct)*1.0e6_wp
         else
            write (*, '(a24, i4, a22)') hdr, ilvl, '  (insufficient data)'
         end if
      end do

      do iacc = 1, n_acc
         row = acc_times(iacc, 1:n_struct)
         call collect_valid_points(n_struct, n_atoms, row, t_min, real_n, raw_t, n_valid)
         if (n_valid >= 4) then
            call fit_power_law(n_valid, real_n, raw_t, exponent, r_sq, coeff_a)
            write (*, '(a24, i4, f10.3, f10.4, es12.3, f12.4, f12.4)') &
               acc_labels(iacc), acc_deriv(iacc), exponent, r_sq, coeff_a, &
               row(1)*1.0e6_wp, row(n_struct)*1.0e6_wp
         else
            write (*, '(a24, i4, a22)') acc_labels(iacc), acc_deriv(iacc), &
               '  (insufficient data)'
         end if
      end do

   end subroutine report_lsf_timings

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

   !> Fill per-atom CPCM radii, turning a failed lookup into a test failure.
   subroutine fill_cpcm_radii(mol, radii, error)
      !> Structure whose per-atom radii are filled
      type(structure_type), intent(in) :: mol
      !> Allocated on exit to mol%nat
      real(wp), allocatable, intent(out) :: radii(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: err
      integer :: iat

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)), err)
         if (allocated(err)) then
            call test_failed(error, "radius lookup failed: "//trim(err%message))
            return
         end if
      end do
   end subroutine fill_cpcm_radii

end module test_cavity_drop_timings
