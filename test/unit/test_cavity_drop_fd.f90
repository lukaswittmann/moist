!> Finite-difference checks of the DROP nuclear derivatives (not a suite)
!>
!> - [[check_drop_totals_gradient]]: total area and volume, 4-point stencil,
!>   optionally on the anchors where branch weights and gates act
!> - [[check_drop_row_gradient]]: every row quantity, below
!>
!> Per-point check:
!>
!> - every atom and direction: 6-point stencil of the rebuilt cavity
!> - rows matched across the stencil by `cavity%numbering`, kept in reference
!>   order; a row is compared only if present at every stencil geometry
!> - per row: xyz, normal, cpjac, wleb, w_f, xi, f, a, v, compared in
!>   that order, each over all atoms, directions and points
!> - a compared row must keep its anchor branch count across the stencil
!> - optional per-point floor on the absolute threshold of the scalar
!>   quantities, from the stencil's own value noise
module test_cavity_drop_fd
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use testdrive, only: error_type, check, test_failed
   use test_helpers, only: build_numbering_map, fd_probe_type, fd_derivative, fd_compare
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type, drop_solve_parameters_type
   use moist_radii, only: default_cpcm_radii
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_grid_model_symmetry_type, only: grid_model_symmetry_type
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: check_drop_row_gradient, check_drop_totals_gradient, drop_fd_param, drop_fd_svdw

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Default SvdW blend steepness
   real(wp), parameter :: DEFAULT_BLEND_K = 2.5_wp
   !> Default SvdW three-body weight
   real(wp), parameter :: DEFAULT_BLEND_3B = 1.0_wp
   !> Default Lebedev grid size
   integer, parameter :: DEFAULT_NUM_LEB = 50
   !> Default projection level
   integer, parameter :: DEFAULT_PROJ_LEVEL = 1
   !> Default weight prune level
   integer, parameter :: DEFAULT_PRUNE_LEVEL = 4
   !> Projection tolerance
   real(wp), parameter :: PROJ_TOL = 1E-14_wp
   !> Projection iteration cap
   integer, parameter :: PROJ_MAXITER = 50

   !> Default FD step, bohr
   real(wp), parameter :: DEFAULT_STEP_SIZE = 2.5E-4_wp
   !> Default FD step of the totals check, bohr
   real(wp), parameter :: DEFAULT_TOTALS_STEP_SIZE = 1.0E-3_wp
   !> Stencil points of the totals check
   integer, parameter :: TOTALS_STENCIL_POINTS = 4
   !> Stencil points
   !>
   !> - O(h^6) central difference, as hessian_fixed: with the 79 branched anchors of the
   !>   level-3 cross the O(h^4) truncation error exceeds REL_THR (cpjac 8.6e-8 relative)
   integer, parameter :: STENCIL_POINTS = 6
   !> Absolute threshold
   real(wp), parameter :: ABS_THR = 5.0E-9_wp
   !> Relative threshold
   real(wp), parameter :: REL_THR = 5.0E-8_wp
   !> Standard deviations of the stencil noise the noise floor admits
   real(wp), parameter :: NOISE_SIGMAS = 5.0_wp

   !> Most branched anchors an anchor subset keeps; more are thinned evenly
   integer, parameter :: SUBSET_MAX_BRANCHED = 48
   !> An anchor subset keeps every this many single-row anchors with a partly
   !> open gate or branch weight
   integer, parameter :: SUBSET_ACTIVE_STRIDE = 8
   !> Pseudo-random share of the other anchors an anchor subset keeps
   real(wp), parameter :: SUBSET_CONTROL_SHARE = 0.1_wp
   !> Largest change of an anchor's area when projected without the others (bohr^2)
   !>
   !> - measured: 1e-15 on SvdW and CFC, 2.5e-12 on SSES, whose level-set
   !>   batches depend on the strip
   real(wp), parameter :: SUBSET_AREA_THR = 1.0E-10_wp

   !> Per-row quantities
   integer, parameter :: NQUANT = 9
   !> Values per row
   integer, parameter :: NVAL = 13
   !> First row value of each quantity, and one past the last
   integer, parameter :: QFIRST(NQUANT + 1) = [1, 4, 7, 8, 9, 10, 11, 12, 13, 14]
   !> Failure label of each quantity
   character(len=*), parameter :: QLABEL(NQUANT) = [character(len=40) :: &
                                                    "Grid point position gradient mismatch", &
                                                    "Surface normal gradient mismatch", &
                                                    "Jacobian scaling gradient mismatch", &
                                                    "Lebedev weight gradient mismatch", &
                                                    "w_f switching function gradient mismatch", &
                                                    "Gaussian width (xi) gradient mismatch", &
                                                    "Switching function gradient mismatch", &
                                                    "Surface area gradient mismatch", &
                                                    "Volume gradient mismatch"]

   !> Cavity rows at a displacement of one atom along one axis
   !>
   !> - component `(m - 1)*nrow + irow`: row value `m` of reference row `irow`
   type, extends(fd_probe_type) :: drop_row_probe_type
      !> Cavity rebuilt at every sample
      type(cavity_type_drop), pointer :: cavity => null()
      !> Reference geometry
      type(structure_type) :: mol
      !> Displaced atom
      integer :: iat = 0
      !> Displaced axis
      integer :: idir = 0
      !> Persistent numbering -> reference row, 0 if absent
      integer, allocatable :: numbering_to_idx(:)
      !> Anchor branch count per reference row
      integer, allocatable :: bc_ref(:)
      !> Branch count equal to `bc_ref` at every sample so far
      logical, allocatable :: bc_kept(:)
   contains
      !> Rows at one displacement
      procedure :: sample => sample_rows
   end type drop_row_probe_type

   !> Total area and volume at a displacement of one atom along one axis
   !>
   !> - components: `sum(a)`, `sum(v)` over the rebuilt grid
   type, extends(fd_probe_type) :: drop_totals_probe_type
      !> Cavity rebuilt at every sample
      type(cavity_type_drop), pointer :: cavity => null()
      !> Reference geometry
      type(structure_type) :: mol
      !> Displaced atom
      integer :: iat = 0
      !> Displaced axis
      integer :: idir = 0
   contains
      !> Totals at one displacement
      procedure :: sample => sample_totals
   end type drop_totals_probe_type

contains

   !* ================================================================================= *!
   !*                                      Driver                                       *!
   !* ================================================================================= *!

   !> DROP settings of the per-point check: tight projection, fine grid
   !>
   !> - 50 Lebedev points, projection level 1, prune level 4, tolerance 1e-14,
   !>   50 Newton iterations; callers change single fields
   function drop_fd_param() result(param)
      !> DROP settings
      type(moist_cavity_drop_parameters_type) :: param

      param = moist_cavity_drop_parameters_type(num_leb=DEFAULT_NUM_LEB, do_fine=.true., tolerance=PROJ_TOL, &
                                                solve=drop_solve_parameters_type(max_iter=PROJ_MAXITER), &
                                                proj_level=DEFAULT_PROJ_LEVEL, &
                                                wleb_prune_level=DEFAULT_PRUNE_LEVEL)
   end function drop_fd_param

   !> SvdW level set of the per-point check
   !>
   !> @param[in] blend_k   blend steepness, default 2.5
   !> @param[in] blend_3b  three-body weight, default 1
   function drop_fd_svdw(blend_k, blend_3b) result(lsf)
      !> Blend steepness
      real(wp), intent(in), optional :: blend_k
      !> Three-body weight
      real(wp), intent(in), optional :: blend_3b
      !> SvdW level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf

      real(wp) :: k, w3

      k = DEFAULT_BLEND_K
      if (present(blend_k)) k = blend_k
      w3 = DEFAULT_BLEND_3B
      if (present(blend_3b)) w3 = blend_3b
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=k, blend_3b=w3))
   end function drop_fd_svdw

   !> Check every per-point nuclear derivative of a DROP cavity against finite differences
   !>
   !> @param[out] error              test failure
   !> @param[in]  mol                molecular structure
   !> @param[in]  param              DROP settings, default [[drop_fd_param]]
   !> @param[in]  lsf                level set, default [[drop_fd_svdw]]
   !> @param[in]  require_branching  reference build and every stencil keep a branched anchor
   !> @param[in]  step               FD step in bohr, default 2.5e-4
   !> @param[in]  noise_floor        absolute threshold of the scalar quantities raised per point
   !>                                to [[stencil_noise_floor]], default off
   !> @param[in]  symmetry           symmetry model of the cavity, updated at every
   !>                                stencil geometry; default none (lab frame)
   subroutine check_drop_row_gradient(error, mol, param, lsf, require_branching, step, noise_floor, symmetry)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> DROP settings
      type(moist_cavity_drop_parameters_type), intent(in), optional :: param
      !> Level set
      class(moist_cavity_drop_lsf_type), intent(in), optional :: lsf
      !> Require branched anchors
      logical, intent(in), optional :: require_branching
      !> FD step
      real(wp), intent(in), optional :: step
      !> Per-point noise floor
      logical, intent(in), optional :: noise_floor
      !> Symmetry model of the cavity
      type(grid_model_symmetry_type), intent(in), optional :: symmetry

      !> Run context borrowed by the cavity
      type(moist_context_type), target :: ctx
      type(cavity_type_drop), allocatable, target :: cavity
      type(drop_row_probe_type) :: probe
      type(mctc_error), allocatable :: cavity_error
      !> Analytic and numeric derivatives (NVAL*nrow, ndim, nat)
      real(wp), allocatable :: analytic(:, :, :), numeric(:, :, :)
      !> Component compared (NVAL*nrow, ndim, nat)
      logical, allocatable :: valid(:, :, :)
      !> Branch count kept across the stencil (nrow, ndim, nat)
      logical, allocatable :: bc_kept(:, :, :)
      !> Reference rows (NVAL*nrow), stencil samples (NVAL*nrow, STENCIL_POINTS)
      real(wp), allocatable :: centre(:), samples(:, :)
      !> Sixth differences (NVAL*nrow, ndim, nat), absolute floor (NVAL*nrow)
      real(wp), allocatable :: d6(:, :, :), thr_floor(:)
      character(len=:), allocatable :: label
      real(wp) :: h
      logical :: branching, use_floor
      integer :: nrow, jgrid, irow, num_idn, iat, idir, iq, m, lo, hi, nflipped, nbranched

      call new_context(ctx, nthreads=0, verbosity=0)
      branching = .false.
      if (present(require_branching)) branching = require_branching
      h = DEFAULT_STEP_SIZE
      if (present(step)) h = step
      use_floor = .false.
      if (present(noise_floor)) use_floor = noise_floor

      allocate (cavity)
      block
         type(moist_cavity_drop_parameters_type) :: drop_param

         if (present(param)) then
            drop_param = param
         else
            drop_param = drop_fd_param()
         end if
         if (present(lsf)) then
            call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=lsf, &
                                 error=cavity_error, param=drop_param, ctx=ctx, grid_model=symmetry)
         else
            call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=drop_fd_svdw(), &
                                 error=cavity_error, param=drop_param, ctx=ctx, grid_model=symmetry)
         end if
      end block
      if (allocated(cavity_error)) call test_failed(error, "Failed to initialize cavity: "//cavity_error%message)
      if (allocated(error)) return
      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) call test_failed(error, "Failed to build cavity: "//cavity_error%message)
      if (allocated(error)) return
      ! Empty reference grid: fail before the 18*nat stencil rebuilds
      if (cavity%ngrid == 0) then
         call test_failed(error, "Reference cavity has no grid points")
         return
      end if

      if (branching) then
         if (.not. allocated(cavity%offset)) then
            call test_failed(error, "Fixture requested branching but has no row layout")
            return
         end if
         if (.not. any(cavity%offset(2:cavity%anchor%ngrid + 1) - cavity%offset(1:cavity%anchor%ngrid) > 1)) then
            call test_failed(error, "Fixture requested branching but no anchor branched; "// &
                             "the branch post-pass would go untested")
            return
         end if
      end if

      ! Reference rows
      nrow = cavity%ngrid
      call build_numbering_map(cavity%numbering(1:nrow), probe%numbering_to_idx)
      allocate (probe%bc_ref(nrow), source=0)
      allocate (centre(NVAL*nrow), source=ieee_value(0.0_wp, ieee_quiet_nan))
      do jgrid = 1, cavity%ngrid
         irow = reference_row(probe%numbering_to_idx, cavity%numbering(jgrid))
         if (irow <= 0) cycle
         probe%bc_ref(irow) = branch_count(cavity, jgrid)
         centre(irow::nrow) = sample_row(cavity, jgrid)
      end do

      call cavity%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      if (.not. allocated(cavity%xyz1_rA)) then
         call test_failed(error, "xyz1_rA not allocated after get_gradient")
         return
      end if

      allocate (analytic(NVAL*nrow, ndim, mol%nat), source=ieee_value(0.0_wp, ieee_quiet_nan))
      do jgrid = 1, cavity%ngrid
         num_idn = cavity%numbering(jgrid)
         irow = reference_row(probe%numbering_to_idx, num_idn)
         if (irow <= 0) cycle
         do iat = 1, mol%nat
            do idir = 1, ndim
               analytic(irow::nrow, idir, iat) = analytic_row(cavity, iat, idir, jgrid)
            end do
         end do
      end do

      ! Finite differences of every atom and direction
      allocate (numeric(NVAL*nrow, ndim, mol%nat), valid(NVAL*nrow, ndim, mol%nat))
      allocate (bc_kept(nrow, ndim, mol%nat))
      allocate (samples(NVAL*nrow, STENCIL_POINTS), d6(NVAL*nrow, ndim, mol%nat))
      allocate (probe%bc_kept(nrow))
      probe%cavity => cavity
      probe%mol = mol
      do iat = 1, mol%nat
         do idir = 1, ndim
            probe%iat = iat
            probe%idir = idir
            probe%bc_kept(:) = .true.
            call fd_derivative(probe, h, STENCIL_POINTS, numeric(:, idir, iat), &
                               valid(:, idir, iat), error, samples)
            if (allocated(error)) return
            bc_kept(:, idir, iat) = probe%bc_kept
            ! Samples at +3h, +2h, +h, -h, -2h, -3h
            if (use_floor) d6(:, idir, iat) = sixth_difference(samples(:, 6), samples(:, 5), samples(:, 4), &
                                                               centre, samples(:, 3), samples(:, 2), samples(:, 1))
         end do
      end do

      allocate (thr_floor(NVAL*nrow), source=0.0_wp)
      if (use_floor) then
         do iq = 1, NQUANT
            if (QFIRST(iq + 1) - QFIRST(iq) > 1) cycle
            lo = (QFIRST(iq) - 1)*nrow + 1
            hi = QFIRST(iq)*nrow
            call stencil_noise_floor(d6(lo:hi, :, :), valid(lo:hi, :, :), h, thr_floor(lo:hi))
         end do
      end if

      ! Every nuclear direction must keep at least one point; none left would pass
      do iat = 1, mol%nat
         do idir = 1, ndim
            call check(error, count(valid(1:nrow, idir, iat)) > 0, &
                       "No valid FD comparisons for atom "//to_string(iat)// &
                       " direction "//to_string(idir))
            if (allocated(error)) return
            nflipped = count(valid(1:nrow, idir, iat) .and. .not. bc_kept(:, idir, iat))
            if (nflipped > 0) then
               call test_failed(error, to_string(nflipped)//" compared points change their branch "// &
                                "count across the stencil for atom "//to_string(iat)// &
                                " direction "//to_string(idir))
               return
            end if
            if (branching) then
               nbranched = count(valid(1:nrow, idir, iat) .and. probe%bc_ref > 1)
               if (nbranched < 1) then
                  call test_failed(error, "No branched point survives the stencil for atom "// &
                                   to_string(iat)//" direction "//to_string(idir))
                  return
               end if
            end if
         end do
      end do

      do iq = 1, NQUANT
         do iat = 1, mol%nat
            do idir = 1, ndim
               do m = QFIRST(iq), QFIRST(iq + 1) - 1
                  label = trim(QLABEL(iq))//" for atom "//to_string(iat)//" direction "//to_string(idir)
                  if (QFIRST(iq + 1) - QFIRST(iq) > 1) then
                     label = label//" component "//to_string(m - QFIRST(iq) + 1)
                  end if
                  lo = (m - 1)*nrow + 1
                  hi = m*nrow
                  call fd_compare(error, analytic(lo:hi, idir, iat), numeric(lo:hi, idir, iat), &
                                  valid(lo:hi, idir, iat), max(ABS_THR, thr_floor(lo:hi)), REL_THR, &
                                  label, "point")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do

   end subroutine check_drop_row_gradient

   !> Check the total area and volume gradients of a built DROP cavity against finite differences
   !>
   !> - 4-point central difference of `sum(a)` and `sum(v)` over every atom and axis
   !> - both accumulators compared: `A_tot1_rA`/`V_tot1_rA` and the per-sphere
   !>   sums `sum(asph1_rA)`/`sum(vsph1_rA)`
   !> - with `subset`, both sides on the anchors of [[restrict_to_switching_anchors]];
   !>   the cavity stays restricted
   !>
   !> @param[in,out] cavity  built cavity, rebuilt at the displaced geometries
   !> @param[in]     mol     reference structure
   !> @param[out]    error   test failure
   !> @param[in]     step    FD step in bohr, default 1e-3
   !> @param[in]     subset  check the anchors where branch weights and gates act, default off
   subroutine check_drop_totals_gradient(cavity, mol, error, step, subset)
      !> Built cavity
      type(cavity_type_drop), intent(inout), target :: cavity
      !> Reference structure
      type(structure_type), intent(in) :: mol
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> FD step
      real(wp), intent(in), optional :: step
      !> Check the switching anchors only
      logical, intent(in), optional :: subset

      !> Number of compared totals per direction
      integer, parameter :: NTOT = 4
      !> Failure label of each total
      character(len=*), parameter :: TLABEL(NTOT) = [character(len=36) :: &
                                                     "Total area gradient mismatch", &
                                                     "Total volume gradient mismatch", &
                                                     "Per-sphere area gradient mismatch", &
                                                     "Per-sphere volume gradient mismatch"]
      !> Sampled total each compares against: area, volume, area, volume
      integer, parameter :: TSAMPLE(NTOT) = [1, 2, 1, 2]
      type(drop_totals_probe_type) :: probe
      type(mctc_error), allocatable :: cavity_error
      !> Analytic derivatives (NTOT, ndim, nat)
      real(wp), allocatable :: analytic(:, :, :)
      real(wp) :: numeric(2), h
      logical :: valid(2)
      integer :: iat, idir, it

      h = DEFAULT_TOTALS_STEP_SIZE
      if (present(step)) h = step
      if (present(subset)) then
         if (subset) then
            call restrict_to_switching_anchors(cavity, mol, error)
            if (allocated(error)) return
         end if
      end if

      call cavity%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      allocate (analytic(NTOT, ndim, mol%nat))
      do iat = 1, mol%nat
         do idir = 1, ndim
            analytic(:, idir, iat) = [cavity%A_tot1_rA(idir, iat), cavity%V_tot1_rA(idir, iat), &
                                      sum(cavity%asph1_rA(idir, :, iat)),sum(cavity%vsph1_rA(idir, :, iat))]
         end do
      end do

      probe%cavity => cavity
      probe%mol = mol
      do iat = 1, mol%nat
         do idir = 1, ndim
            probe%iat = iat
            probe%idir = idir
            call fd_derivative(probe, h, TOTALS_STENCIL_POINTS, numeric, valid, error)
            if (allocated(error)) return
            do it = 1, NTOT
               call fd_compare(error, analytic(it:it, idir, iat), numeric(TSAMPLE(it):TSAMPLE(it)), &
                               valid(TSAMPLE(it):TSAMPLE(it)), ABS_THR, REL_THR, &
                               trim(TLABEL(it))//", atom "//to_string(iat)//" axis "//to_string(idir))
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine check_drop_totals_gradient

   !> Rebuild a cavity on the anchors where branch weights and gates act
   !>
   !> - every branched anchor, at most [[SUBSET_MAX_BRANCHED]], beyond that
   !>   evenly thinned in anchor order
   !> - every [[SUBSET_ACTIVE_STRIDE]]-th single-row anchor with a partly open
   !>   gate or branch weight, the first included
   !> - [[SUBSET_CONTROL_SHARE]] of the other anchors, picked by a
   !>   multiplicative hash of the anchor id
   !> - every anchor projects on its own: each kept anchor must keep its rows
   !>   and area of the full build, else the subset checks something else
   !>
   !> @param[in,out] cavity  built cavity, rebuilt on the subset
   !> @param[in]     mol     reference structure
   !> @param[out]    error   test failure
   subroutine restrict_to_switching_anchors(cavity, mol, error)
      !> Built cavity
      type(cavity_type_drop), intent(inout) :: cavity
      !> Reference structure
      type(structure_type), intent(in) :: mol
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: cavity_error
      !> Kept anchor ids
      integer, allocatable :: ids(:)
      !> Row count and area of the full build per anchor id, -1 rows if absent
      integer, allocatable :: ref_rows(:)
      real(wp), allocatable :: ref_area(:)
      integer :: nanch, stride_branched, ia, r0, r1, id, nkeep, nbranched, nactive
      logical :: keep

      nanch = cavity%anchor%ngrid
      stride_branched = count(cavity%offset(2:nanch + 1) - cavity%offset(1:nanch) > 1)
      stride_branched = max(1, (stride_branched + SUBSET_MAX_BRANCHED - 1)/SUBSET_MAX_BRANCHED)
      allocate (ref_rows(cavity%nsph*cavity%param%num_leb), source=-1)
      allocate (ref_area(size(ref_rows)), source=0.0_wp)
      allocate (ids(nanch))
      nkeep = 0
      nbranched = 0
      nactive = 0
      do ia = 1, nanch
         r0 = cavity%offset(ia) + 1
         r1 = cavity%offset(ia + 1)
         id = cavity%anchor%id(ia)
         ref_rows(id) = r1 - r0 + 1
         ref_area(id) = sum(cavity%a(r0:r1))
         if (r1 > r0) then
            nbranched = nbranched + 1
            keep = modulo(nbranched - 1, stride_branched) == 0
         else if (any(partly_open(cavity%w_f0(r0:r1))) .or. any(partly_open(cavity%wbranch(r0:r1)))) then
            nactive = nactive + 1
            keep = modulo(nactive - 1, SUBSET_ACTIVE_STRIDE) == 0
         else
            keep = control_anchor(id)
         end if
         if (keep) then
            nkeep = nkeep + 1
            ids(nkeep) = id
         end if
      end do

      cavity%anchor_subset = ids(1:nkeep)
      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call check(error, cavity%anchor%ngrid == nkeep, "anchor subset: "//to_string(cavity%anchor%ngrid)// &
                 " anchors projected of "//to_string(nkeep)//" kept")
      if (allocated(error)) return
      do ia = 1, cavity%anchor%ngrid
         r0 = cavity%offset(ia) + 1
         r1 = cavity%offset(ia + 1)
         id = cavity%anchor%id(ia)
         call check(error, &
                    r1 - r0 + 1 == ref_rows(id) .and. abs(sum(cavity%a(r0:r1)) - ref_area(id)) <= SUBSET_AREA_THR, &
                    "anchor "//to_string(id)//" projects differently without the others: "// &
                    to_string(r1 - r0 + 1)//" rows and area "//to_string(sum(cavity%a(r0:r1)))//", full build "// &
                    to_string(ref_rows(id))//" rows and area "//to_string(ref_area(id)))
         if (allocated(error)) return
      end do
   end subroutine restrict_to_switching_anchors

   !> Whether a gate or branch weight is partly open
   !>
   !> @param[in] w  gate or branch weight
   !> @returns      `0 < w < 1`
   elemental function partly_open(w) result(partly)
      !> Gate or branch weight
      real(wp), intent(in) :: w
      !> Partly open
      logical :: partly

      partly = w > 0.0_wp .and. w < 1.0_wp
   end function partly_open

   !> Whether an anchor without switching activity joins an anchor subset
   !>
   !> - Fibonacci hash of the id: `id*2654435761 mod 2^32` below the share of `2^32`
   !>
   !> @param[in] id  anchor id
   !> @returns       anchor kept
   pure function control_anchor(id) result(kept)
      !> Anchor id
      integer, intent(in) :: id
      !> Anchor kept
      logical :: kept

      integer(int64), parameter :: golden = 2654435761_int64, word = 4294967296_int64

      kept = real(modulo(int(id, int64)*golden, word), wp) < SUBSET_CONTROL_SHARE*real(word, wp)
   end function control_anchor

   !* ================================================================================= *!
   !*                                     Row probe                                     *!
   !* ================================================================================= *!

   !> Rebuild the cavity with one atom displaced and scatter its rows into reference order
   !>
   !> @param[in,out] self   probe; cavity rebuilt, branch counts tracked
   !> @param[in]     t      displacement of atom `iat` along `idir`
   !> @param[out]    vals   row values, NaN for absent rows
   !> @param[out]    valid  row present
   !> @param[out]    error  test failure, set if the rebuild fails
   subroutine sample_rows(self, t, vals, valid, error)
      !> Probe
      class(drop_row_probe_type), intent(inout) :: self
      !> Displacement
      real(wp), intent(in) :: t
      !> Row values
      real(wp), intent(out) :: vals(:)
      !> Row present
      logical, intent(out) :: valid(:)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol_fd
      type(mctc_error), allocatable :: cavity_error
      integer :: nrow, jgrid, irow

      nrow = size(self%bc_ref)
      if (size(vals) /= NVAL*nrow) then
         call test_failed(error, "sample_rows: "//to_string(size(vals))//" components for "// &
                          to_string(nrow)//" rows")
         return
      end if

      mol_fd = self%mol
      mol_fd%xyz(self%idir, self%iat) = mol_fd%xyz(self%idir, self%iat) + t
      call self%cavity%update(mol_fd, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      vals = ieee_value(0.0_wp, ieee_quiet_nan)
      valid = .false.
      do jgrid = 1, self%cavity%ngrid
         irow = reference_row(self%numbering_to_idx, self%cavity%numbering(jgrid))
         if (irow <= 0) cycle
         vals(irow::nrow) = sample_row(self%cavity, jgrid)
         valid(irow::nrow) = .true.
         self%bc_kept(irow) = self%bc_kept(irow) .and. &
                              branch_count(self%cavity, jgrid) == self%bc_ref(irow)
      end do
   end subroutine sample_rows

   !> Rebuild the cavity with one atom displaced and total its area and volume
   !>
   !> @param[in,out] self   probe; cavity rebuilt
   !> @param[in]     t      displacement of atom `iat` along `idir`
   !> @param[out]    vals   `sum(a)`, `sum(v)`
   !> @param[out]    valid  always set
   !> @param[out]    error  test failure, set if the rebuild fails
   subroutine sample_totals(self, t, vals, valid, error)
      !> Probe
      class(drop_totals_probe_type), intent(inout) :: self
      !> Displacement
      real(wp), intent(in) :: t
      !> Totals
      real(wp), intent(out) :: vals(:)
      !> Totals present
      logical, intent(out) :: valid(:)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol_fd
      type(mctc_error), allocatable :: cavity_error

      vals = ieee_value(0.0_wp, ieee_quiet_nan)
      valid = .false.
      mol_fd = self%mol
      mol_fd%xyz(self%idir, self%iat) = mol_fd%xyz(self%idir, self%iat) + t
      call self%cavity%update(mol_fd, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      vals(1:2) = [sum(self%cavity%a(1:self%cavity%ngrid)), sum(self%cavity%v(1:self%cavity%ngrid))]
      valid(1:2) = .true.
   end subroutine sample_totals

   !* ================================================================================= *!
   !*                                   Row columns                                     *!
   !* ================================================================================= *!

   !> Values of one cavity row, in [[QFIRST]] order
   !>
   !> @param[in] cavity  built cavity
   !> @param[in] jgrid   row
   pure function sample_row(cavity, jgrid) result(row)
      !> Built cavity
      type(cavity_type_drop), intent(in) :: cavity
      !> Row
      integer, intent(in) :: jgrid
      !> Row values
      real(wp) :: row(NVAL)

      row = [cavity%xyz(:, jgrid), cavity%normal0(:, jgrid), &
             cavity%cpjac_scal0(jgrid), cavity%wleb(jgrid), cavity%w_f0(jgrid), cavity%xi0(jgrid), &
             cavity%f(jgrid), cavity%a(jgrid), cavity%v(jgrid)]
   end function sample_row

   !> Analytic derivatives of one cavity row along atom `iat`, axis `idir`, in [[QFIRST]] order
   !>
   !> - `normal1_rA` holds the normal component first and the axis third
   !>
   !> @param[in] cavity  cavity after `get_gradient`
   !> @param[in] iat     atom
   !> @param[in] idir    axis
   !> @param[in] jgrid   row
   pure function analytic_row(cavity, iat, idir, jgrid) result(row)
      !> Cavity after `get_gradient`
      type(cavity_type_drop), intent(in) :: cavity
      !> Atom
      integer, intent(in) :: iat
      !> Axis
      integer, intent(in) :: idir
      !> Row
      integer, intent(in) :: jgrid
      !> Row derivatives
      real(wp) :: row(NVAL)

      row = [cavity%xyz1_rA(:, idir, iat, jgrid), &
             cavity%normal1_rA(:, iat, idir, jgrid), cavity%cpjac_scal1_rA(idir, iat, jgrid), &
             cavity%wleb1_rA(idir, iat, jgrid), cavity%w_f1_rA(idir, iat, jgrid), &
             cavity%xi1_rA(idir, iat, jgrid), cavity%f1_rA(idir, iat, jgrid), &
             cavity%a_i1_rA(idir, iat, jgrid), cavity%v1_rA(idir, iat, jgrid)]
   end function analytic_row

   !> Reference row of a persistent number, 0 if absent from the reference
   !>
   !> @param[in] numbering_to_idx  persistent numbering -> reference row
   !> @param[in] num_idn           persistent number
   pure function reference_row(numbering_to_idx, num_idn) result(irow)
      !> Persistent numbering -> reference row
      integer, intent(in) :: numbering_to_idx(:)
      !> Persistent number
      integer, intent(in) :: num_idn
      !> Reference row
      integer :: irow

      irow = 0
      if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) return
      irow = numbering_to_idx(num_idn)
   end function reference_row

   !> Branch count of the anchor of a cavity row
   !>
   !> @param[in] cavity  built cavity
   !> @param[in] jgrid   row
   pure function branch_count(cavity, jgrid) result(nbranch)
      !> Built cavity
      type(cavity_type_drop), intent(in) :: cavity
      !> Row
      integer, intent(in) :: jgrid
      !> Branch count
      integer :: nbranch

      nbranch = cavity%offset(cavity%anchor_of_row(jgrid) + 1) - cavity%offset(cavity%anchor_of_row(jgrid))
   end function branch_count

   !* ================================================================================= *!
   !*                               Stencil round-off floor                             *!
   !* ================================================================================= *!

   !> Sixth difference of the seven stencil samples, `h^6 f^(6)` plus value noise
   !>
   !> - offsets from the centre value first: each is exact for small `h`
   !>   (Sterbenz), so no `30 c` cancellation adds rounding of `c` itself
   !>
   !> @param[in] m3  value at -3h
   !> @param[in] m2  value at -2h
   !> @param[in] m1  value at -h
   !> @param[in] c   value at the reference geometry
   !> @param[in] p1  value at +h
   !> @param[in] p2  value at +2h
   !> @param[in] p3  value at +3h
   elemental function sixth_difference(m3, m2, m1, c, p1, p2, p3) result(d6)
      !> Value at -3h
      real(wp), intent(in) :: m3
      !> Value at -2h
      real(wp), intent(in) :: m2
      !> Value at -h
      real(wp), intent(in) :: m1
      !> Value at the reference geometry
      real(wp), intent(in) :: c
      !> Value at +h
      real(wp), intent(in) :: p1
      !> Value at +2h
      real(wp), intent(in) :: p2
      !> Value at +3h
      real(wp), intent(in) :: p3
      !> Sixth difference
      real(wp) :: d6

      d6 = ((m3 - c) + (p3 - c)) - 6.0_wp*((m2 - c) + (p2 - c)) + 15.0_wp*((m1 - c) + (p1 - c))
   end function sixth_difference

   !> Round-off floor of the 6-point derivative per point, from the samples
   !>
   !> - value noise far above a few ulp on some points: on the cross, minor
   !>   branches at the concave seams on the ramp of the focal switch, whose
   !>   input `beta2 = lambda_min(w - lambda Hess S)` cancels (`w_f` 0.06-0.4,
   !>   relative noise 1e-13 to 4e-13; `xi` half of it, values up to 18)
   !> - worst point: white noise, independent of the projection tolerance
   !>   (1e-12 to 1e-16), no topology change; analytic derivative right to 1e-8
   !>   relative against least-squares slopes over 400 samples
   !> - stencil: value noise `sigma` turns into derivative noise
   !>   `sqrt(2 (45^2 + 9^2 + 1)) / 60 * sigma / h`, beyond a fixed `ABS_THR`
   !> - `sigma` per point: the sixth difference has noise `sqrt(924) sigma`; its
   !>   lower quartile over the point's coordinates, most of which move the point
   !>   slowly, is `0.3186 sqrt(924) sigma`
   !> - floor: [[NOISE_SIGMAS]] times the derivative noise; above `ABS_THR` only
   !>   where `sigma` exceeds ~1e-13
   !>
   !> @param[in]  d6         sixth differences (npoint, ndim, nat)
   !> @param[in]  valid      compared entries (npoint, ndim, nat)
   !> @param[in]  h          FD step
   !> @param[out] thr_floor  absolute floor per point (npoint), zero without a compared entry
   pure subroutine stencil_noise_floor(d6, valid, h, thr_floor)
      !> Sixth differences
      real(wp), intent(in) :: d6(:, :, :)
      !> Compared entries
      logical, intent(in) :: valid(:, :, :)
      !> FD step
      real(wp), intent(in) :: h
      !> Absolute floor per point
      real(wp), intent(out) :: thr_floor(:)

      real(wp), parameter :: d6_noise = 0.3186_wp*sqrt(924.0_wp)
      real(wp), parameter :: fd_noise = sqrt(2.0_wp*(45.0_wp**2 + 9.0_wp**2 + 1.0_wp))/60.0_wp
      real(wp) :: v(size(d6, 2)*size(d6, 3)), x
      integer :: ipt, idir, iat, m, i, j

      do ipt = 1, size(d6, 1)
         m = 0
         do iat = 1, size(d6, 3)
            do idir = 1, size(d6, 2)
               if (.not. valid(ipt, idir, iat)) cycle
               m = m + 1
               v(m) = abs(d6(ipt, idir, iat))
            end do
         end do
         thr_floor(ipt) = 0.0_wp
         if (m == 0) cycle
         ! Insertion sort, a handful of entries
         do i = 2, m
            x = v(i)
            j = i - 1
            do while (j >= 1)
               if (v(j) <= x) exit
               v(j + 1) = v(j)
               j = j - 1
            end do
            v(j + 1) = x
         end do
         ! Lower quartile, interpolated
         x = 0.25_wp*real(m - 1, wp) + 1.0_wp
         i = int(x)
         x = v(i) + (x - real(i, wp))*(v(min(i + 1, m)) - v(i))
         thr_floor(ipt) = NOISE_SIGMAS*fd_noise*(x/d6_noise)/h
      end do
   end subroutine stencil_noise_floor

end module test_cavity_drop_fd
