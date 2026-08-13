!> Shared support routines for the unit-test suite; not a test itself
!>
!> Provides:
!>   * `center_at_origin(mol)` - centroid-shift a structure_type
!>   * `get_test_structures(mols, n)` - sample `n` MB16-43 + Heavy28 +
!>                                         Amino20x4 + But14diol + UPU23
!>                                         records (20% each, no replacement)
!>   * `get_test_radii(mol, radii)` - CPCM-table radii
!>   * `get_test_points(mol, points, n)` - `n` deterministic random sampling
!>                                         points inside/near mol's box
!>   * `get_test_cavity_iswig(mol, cavity, error, ...)` - ready-to-use iSwiG surface
!>   * `build_test_cavity(mol, nleb, ctx, ...)` - iSwiG surface on a caller-owned
!>                                         context and COSMO radii
!>   * `make_charge_coupling(qat, coupling)` - single-column charge coupling
!>   * `get_test_cross(mol)` - five-carbon cross with concave seams
!>   * `check_moist_error(error, err, context)` - moist error -> testdrive failure
!>   * `fd4_scalar(fpp, fp, fm, fmm, h)` - 4-point central FD formula
!>   * `fd4_offsets` - the matching stencil offsets, in units of h
!>   * `rel_deviation(a, b)` - |a - b| / (1 + |b|)
!>   * `fill_legacy_radii(mol, radii, error)` - legacy per-element radius table
!>   * `build_numbering_map(numbering, map)` - persistent grid numbering ->
!>                                             current array index
!>
!> No global Fortran RNG state is touched (self-contained LCG), so the
!> point and structure samplers are safe under parallel test execution.
module test_helpers
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use mctc_env_error, only: moist_error_type => error_type
   use mstore, only: get_structure
   use mstore_data_record, only: record_type
   use mstore_mb16_43, only: get_mb16_43_records
   use mstore_heavy28, only: get_heavy28_records
   use mstore_amino20x4, only: get_amino20x4_records
   use mstore_but14diol, only: get_but14diol_records
   use mstore_upu23, only: get_upu23_records
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_context, only: moist_context_type, new_context
   use moist_radii, only: default_cpcm_radii, radius_type, new_radii_custom_atoms, &
                          radius_type_static, new_cosmo_radii
   use moist_channels, only: coupling_type
   use moist_data_radii_legacy, only: get_radius_func
   use testdrive, only: error_type, test_failed
   implicit none
   private

   public :: center_at_origin
   public :: get_test_structures
   public :: get_test_radii
   public :: get_test_points
   public :: get_test_cavity_iswig
   public :: build_test_cavity
   public :: make_charge_coupling
   public :: get_test_cross
   public :: fd4_scalar
   public :: fd4_offsets
   public :: rel_deviation
   public :: check_moist_error
   public :: fill_legacy_radii
   public :: build_numbering_map

   !> Default n for get_test_structures (must be a multiple of 5).
   integer, parameter :: default_n_structures = 5
   !> Default n for get_test_points.
   integer, parameter :: default_n_points = 7
   !> Default Lebedev order for get_test_cavity_iswig.
   integer, parameter :: default_nleb = 26

   !> Stencil offsets, in units of h, matching `fd4_scalar`'s argument order
   real(wp), parameter :: fd4_offsets(4) = [2.0_wp, 1.0_wp, -1.0_wp, -2.0_wp]

   !> The 5 mstore collections that get_test_structures samples from.
   integer, parameter :: n_datasets = 5
   character(len=*), parameter :: datasets(n_datasets) = [character(len=10):: &
                                                          "MB16-43", "Heavy28", "Amino20x4", "But14diol", "UPU23"]

contains

   !> Translate `mol` so its arithmetic centroid sits at the origin.
   !> Pure positional shift; atomic identities and ordering preserved.
   !> @param[inout] mol  structure whose %xyz is shifted in place
   subroutine center_at_origin(mol)
      !> Structure whose Cartesian coordinates are shifted to put the
      !> arithmetic centroid at the origin.
      type(structure_type), intent(inout) :: mol
      real(wp) :: centroid(3)
      integer :: iat

      centroid = 0.0_wp
      do iat = 1, mol%nat
         centroid = centroid + mol%xyz(:, iat)
      end do
      centroid = centroid/real(mol%nat, wp)
      do iat = 1, mol%nat
         mol%xyz(:, iat) = mol%xyz(:, iat) - centroid
      end do
   end subroutine center_at_origin

   !> Populate `structures` with `n` mstore records sampled evenly across
   !> 5 collections (MB16-43, Heavy28, Amino20x4, But14diol, UPU23). Each
   !> collection contributes 20% of the total via Fisher-Yates shuffle on
   !> its record list (no duplicates within a collection). Different `n`
   !> values produce different but reproducible samples; the LCG state is
   !> self-contained.
   !>
   !> `n` must be a multiple of 5 and >= 5. Default = 10 (2 per set).
   !> @param[out] structures  allocated array of mstore-sourced structures
   !> @param[in]  n           optional total count; default 10
   subroutine get_test_structures(structures, n)
      !> Output array of mstore-sourced structures.
      type(structure_type), allocatable, intent(out) :: structures(:)
      !> Optional total count (must be multiple of 5).
      integer, optional, intent(in) :: n

      integer :: total, per_set, ids, k
      integer(int64) :: rng
      type(record_type), allocatable :: records(:)
      integer, allocatable :: order(:)

      total = default_n_structures
      if (present(n)) total = n
      if (total < n_datasets .or. mod(total, n_datasets) /= 0) then
         error stop "get_test_structures: n must be a positive multiple of 5"
      end if
      per_set = total/n_datasets

      allocate (structures(total))
      !* LCG seed mixes the total count so different n values yield
      !* different draws while each n is fully reproducible.
      rng = int(total, int64)*1009_int64 + 12345_int64

      k = 0
      do ids = 1, n_datasets
         call load_dataset(trim(datasets(ids)), records)
         call shuffle_indices(size(records), rng, order)
         do k = 1, per_set
            call get_structure(structures((ids - 1)*per_set + k), &
                               trim(datasets(ids)), &
                               trim(records(order(k))%id))
         end do
         deallocate (records, order)
      end do
   end subroutine get_test_structures

   !> Return CPCM-table radii for `mol` using the project's standard
   !> `default_cpcm_radii()` model.
   !> @param[in]  mol    structure to look up radii for
   !> @param[out] radii  allocated to size mol%nat; filled with CPCM radii
   subroutine get_test_radii(mol, radii)
      !> Structure whose per-atom radii are to be filled.
      type(structure_type), intent(in) :: mol
      !> Allocated on exit to mol%nat; filled with CPCM-table radii.
      real(wp), allocatable, intent(out) :: radii(:)
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      model = default_cpcm_radii()
      call model%update(mol, err)
      if (allocated(err)) error stop "get_test_radii: "//trim(err%message)
      allocate (radii, source=model%f0)
   end subroutine get_test_radii

   !> Generate `n` deterministic random sampling points inside and near
   !> `mol`. Points are drawn uniformly from the atom bounding box padded
   !> by +/-2 bohr in each axis; candidates within 0.2 bohr of any nucleus
   !> are rejected. A self-contained 64-bit LCG seeded from mol%nat keeps
   !> results reproducible without touching the global Fortran RNG state,
   !> so the routine is safe under parallel test execution.
   !> @param[in]  mol     structure used to anchor the bounding box
   !> @param[out] points  (3, n) allocatable array of sampling points
   !> @param[in]  n       optional point count; default 12
   subroutine get_test_points(mol, points, n)
      !> Structure used to seed the bounding-box random sampler.
      type(structure_type), intent(in) :: mol
      !> Allocated on exit to (3, n); filled with sampling points.
      real(wp), allocatable, intent(out) :: points(:, :)
      !> Optional point count; default `default_n_points`.
      integer, optional, intent(in) :: n

      real(wp), parameter :: pad = 2.0_wp
      real(wp), parameter :: min_dist = 2.0e-1_wp
      integer, parameter :: max_attempts = 10000

      integer :: total
      integer(int64) :: rng
      real(wp) :: box_min(3), box_max(3), span(3), cand(3)
      integer :: ax, attempt, accepted

      total = default_n_points
      if (present(n)) total = n
      if (total < 1) error stop "get_test_points: n must be >= 1"
      if (mol%nat < 1) error stop "get_test_points: empty mol"

      do ax = 1, 3
         box_min(ax) = minval(mol%xyz(ax, :)) - pad
         box_max(ax) = maxval(mol%xyz(ax, :)) + pad
      end do
      span = box_max - box_min

      allocate (points(3, total))
      !* Seed mixes mol%nat with a fixed offset so different-sized
      !* molecules get distinct draws while each size is reproducible.
      rng = int(mol%nat, int64)*97_int64 + 12345_int64
      accepted = 0
      do attempt = 1, max_attempts
         do ax = 1, 3
            cand(ax) = box_min(ax) + lcg_uniform(rng)*span(ax)
         end do
         if (.not. far_from_all_atoms(cand, mol%xyz, min_dist)) cycle
         accepted = accepted + 1
         points(:, accepted) = cand
         if (accepted == total) return
      end do
      error stop "get_test_points: not enough valid points"
   end subroutine get_test_points

   !> Build an iSwiG surface for `mol`
   !>
   !> @param[in]  mol           Structure to wrap
   !> @param[out] cavity        Constructed iSwiG cavity
   !> @param[out] error         Error handling
   !> @param[in]  nleb          Lebedev order; default `default_nleb`
   !> @param[in]  radius_model  Radius model; default CPCM-table per-atom radii
   !> @param[in]  cut_f         Switching-factor cutoff; cavity default if absent
   subroutine get_test_cavity_iswig(mol, cavity, error, nleb, radius_model, cut_f)
      !> Structure to wrap.
      type(structure_type), intent(in) :: mol
      !> Constructed iSwiG cavity.
      type(cavity_type_iswig), intent(out) :: cavity
      !> Error handling.
      type(moist_error_type), allocatable, intent(out) :: error
      !> Optional Lebedev order.
      integer, intent(in), optional :: nleb
      !> Optional radius model; CPCM-table per-atom radii if absent.
      class(radius_type), intent(in), optional :: radius_model
      !> Optional switching-factor cutoff;  grid points with `f <= cut_f` are
      !> dropped at construction. Raise it to keep only the exposed surface.
      real(wp), intent(in), optional :: cut_f

      !> Per-atom radii of the default model.
      real(wp), allocatable :: radii(:)
      !> Default radius model, built only when none was supplied.
      class(radius_type), allocatable :: default_model
      !> Resolved Lebedev order.
      integer :: num_leb
      !> Run context
      type(moist_context_type), pointer :: ctx

      num_leb = default_nleb
      if (present(nleb)) num_leb = nleb

      allocate (ctx)
      call new_context(ctx)
      if (present(radius_model)) then
         call new_cavity_iswig(cavity, ctx, nleb=num_leb, cut_f=cut_f, &
                               radius_model=radius_model, error=error)
      else
         call get_test_radii(mol, radii)
         call new_radii_custom_atoms(radii, default_model, error)
         if (allocated(error)) return
         call new_cavity_iswig(cavity, ctx, nleb=num_leb, cut_f=cut_f, &
                               radius_model=default_model, error=error)
      end if
      if (allocated(error)) return

      call cavity%update(mol, error=error)
   end subroutine get_test_cavity_iswig

   !> Build a COSMO-radii iSwiG test cavity for a given molecule and Lebedev
   !> grid size.
   !>
   !> Unlike `get_test_cavity_iswig`, the run context is created and owned by
   !> the caller so that it outlives the cavity borrowing it, so several
   !> cavities can share one context, and so the same context can be handed to
   !> a solvation-model component.
   !>
   !> @param[in]  mol          Molecular structure
   !> @param[in]  nleb         Lebedev grid size
   !> @param[in]  ctx          Run context owned by the caller, borrowed by the cavity
   !> @param[out] radius_model Radius model storage
   !> @param[out] cavity       Constructed cavity
   !> @param[out] error        Error handling
   subroutine build_test_cavity(mol, nleb, ctx, radius_model, cavity, error)

      !> Molecular structure
      type(structure_type), intent(in) :: mol

      !> Lebedev grid size
      integer, intent(in) :: nleb

      !> Run context owned by the caller, borrowed by the cavity
      type(moist_context_type), intent(in), target :: ctx

      !> Radius model storage
      type(radius_type_static), intent(out) :: radius_model

      !> Constructed cavity
      type(cavity_type_iswig), intent(out) :: cavity

      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      call new_cosmo_radii(radius_model)
      call new_cavity_iswig(cavity, ctx, nleb=nleb, radius_model=radius_model, error=error)
      if (allocated(error)) return

      call cavity%update(mol, error=error)

   end subroutine build_test_cavity

   !> Build a single-column wavefunction charge array.
   !>
   !> @param[in]  qat      Atomic charges
   !> @param[out] coupling Wavefunction to populate
   subroutine make_charge_coupling(qat, coupling)

      !> Atomic charges
      real(wp), intent(in) :: qat(:)

      !> Wavefunction to populate
      type(coupling_type), intent(out) :: coupling

      coupling%electrostatics%qat = reshape(qat, [size(qat), 1])

   end subroutine make_charge_coupling

   !> Five-carbon cross, converted to bohr.
   !>
   !> Deliberately concave seams between the four outer atoms: at the
   !> unconditional-multistart projection level this geometry produces branched
   !> anchors, which is what the warm-start and gradient FD suites need.
   !>
   !> @param[out] mol  Five-carbon cross structure, coordinates in bohr
   subroutine get_test_cross(mol)
      !> Resulting structure; coordinates in bohr.
      type(structure_type), intent(out) :: mol

      call new(mol, [6, 6, 6, 6, 6], reshape([ &
                                             0.00_wp, 4.21_wp, 0.00_wp, &
                                             0.00_wp, 0.00_wp, 4.22_wp, &
                                             0.00_wp, -4.18_wp, 0.00_wp, &
                                             0.00_wp, 0.00_wp, -4.15_wp, &
                                             0.02_wp, 0.10_wp, -0.20_wp], &
                                             [3, 5])*aatoau)
   end subroutine get_test_cross

   !> Turn a moist error into a testdrive test failure.
   !>
   !> Safe to call unconditionally: an unallocated `err` is a no-op, so the
   !> idiom at every call site collapses to
   !> `call check_moist_error(error, err); if (allocated(error)) return`.
   !> This is the one assertion in an otherwise fixture-only module; it lives
   !> here because the alternative is the same four lines copied into every
   !> suite that touches a moist routine.
   !>
   !> @param[out] error    Test failure, allocated only when `err` was
   !> @param[in]  err      Moist error to translate
   !> @param[in]  context  Optional prefix, e.g. the operation that failed
   subroutine check_moist_error(error, err, context)
      !> Test failure; allocated only when `err` is.
      type(error_type), allocatable, intent(out) :: error
      !> Moist error to translate; unallocated means success.
      type(moist_error_type), allocatable, intent(in) :: err
      !> Optional prefix describing what was being attempted.
      character(len=*), intent(in), optional :: context

      if (.not. allocated(err)) return
      if (present(context)) then
         call test_failed(error, context//": "//trim(err%message))
      else
         call test_failed(error, trim(err%message))
      end if
   end subroutine check_moist_error

   !> 4-point central finite-difference formula:
   !>   f'(x) ~ (-f(x+2h) + 8 f(x+h) - 8 f(x-h) + f(x-2h)) / (12 h).
   !> Truncation O(h^4 f^(5)); useful for FD-checking analytic derivatives.
   pure real(wp) function fd4_scalar(fpp, fp, fm, fmm, h) result(df)
      !> Value at x + 2h.
      real(wp), intent(in) :: fpp
      !> Value at x + h.
      real(wp), intent(in) :: fp
      !> Value at x - h.
      real(wp), intent(in) :: fm
      !> Value at x - 2h.
      real(wp), intent(in) :: fmm
      !> Step size h.
      real(wp), intent(in) :: h

      df = (-fpp + 8.0_wp*fp - 8.0_wp*fm + fmm)/(12.0_wp*h)
   end function fd4_scalar

   !> Deviation of `a` from reference `b`, relative but safe near zero
   !>   |a - b| / (1 + |b|).
   !>
   !> @param[in] a  Value under test
   !> @param[in] b  Reference value
   elemental pure real(wp) function rel_deviation(a, b) result(dev)
      !> Value under test.
      real(wp), intent(in) :: a
      !> Reference value.
      real(wp), intent(in) :: b

      dev = abs(a - b)/(1.0_wp + abs(b))
   end function rel_deviation

   !* ===================================================================
   !*                          Private helpers
   !* ===================================================================

   !> Dispatch to the per-dataset records getter. Caller frees `records`.
   subroutine load_dataset(name, records)
      character(len=*), intent(in) :: name
      type(record_type), allocatable, intent(out) :: records(:)

      select case (name)
      case ("MB16-43")
         call get_mb16_43_records(records)
      case ("Heavy28")
         call get_heavy28_records(records)
      case ("Amino20x4")
         call get_amino20x4_records(records)
      case ("But14diol")
         call get_but14diol_records(records)
      case ("UPU23")
         call get_upu23_records(records)
      case default
         error stop "load_dataset: unknown collection '"//trim(name)//"'"
      end select
   end subroutine load_dataset

   !> Fisher-Yates: produce a random permutation of [1..n] using `rng`.
   !> Allocates `order(n)` on output. Stateful in `rng`, no globals touched.
   subroutine shuffle_indices(n, rng, order)
      integer, intent(in) :: n
      integer(int64), intent(inout) :: rng
      integer, allocatable, intent(out) :: order(:)
      integer :: i, j, tmp

      allocate (order(n))
      do i = 1, n
         order(i) = i
      end do
      !* Standard Fisher-Yates: for i = n downto 2, swap order(i) with
      !* order(j) where j uniform in [1, i].
      do i = n, 2, -1
         j = 1 + int(lcg_uniform(rng)*real(i, wp))
         if (j > i) j = i   ! guard against rounding to exactly i
         tmp = order(i)
         order(i) = order(j)
         order(j) = tmp
      end do
   end subroutine shuffle_indices

   !> True iff `point` is at least `min_dist` from every column of `centers`.
   logical function far_from_all_atoms(point, centers, min_dist) result(ok)
      real(wp), intent(in) :: point(3)
      real(wp), intent(in) :: centers(:, :)
      real(wp), intent(in) :: min_dist
      real(wp) :: min_dist_sq
      integer :: iat

      min_dist_sq = min_dist*min_dist
      ok = .true.
      do iat = 1, size(centers, dim=2)
         if (sum((point - centers(:, iat))**2) < min_dist_sq) then
            ok = .false.
            return
         end if
      end do
   end function far_from_all_atoms

   !> 64-bit LCG (Knuth MMIX constants from TAOCP vol 2). Self-contained;
   !> does not touch the global Fortran RNG state, so safe under parallel
   !> test execution.
   !> @param[inout] state  LCG state, advanced by one step
   real(wp) function lcg_uniform(state) result(u)
      integer(int64), intent(inout) :: state
      state = state*6364136223846793005_int64 + 1442695040888963407_int64
      !* Top 31 bits as a nonneg integer; divide by 2^31 to land in [0, 1).
      u = real(ishft(state, -33), wp)/real(2_int64**31, wp)
   end function lcg_uniform

   !> Fill per-atom radii from the legacy per-element table
   !>
   !> The CPCM-flavoured cavity tests want the same radii the legacy code used,
   !> which is not what `get_test_radii` returns.
   !>
   !> @param[in]  mol   Structure whose per-atom radii are filled
   !> @param[out] radii Allocated on exit to `mol%nat`
   !> @param[out] error Test failure, allocated if a lookup fails
   subroutine fill_legacy_radii(mol, radii, error)
      !> Structure whose per-atom radii are filled
      type(structure_type), intent(in) :: mol
      !> Allocated on exit to mol%nat
      real(wp), allocatable, intent(out) :: radii(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Lookup failure from the radius table
      type(moist_error_type), allocatable :: err
      !> Atom index
      integer :: iat

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)), err)
         if (allocated(err)) then
            call test_failed(error, "radius lookup failed: "//trim(err%message))
            return
         end if
      end do
   end subroutine fill_legacy_radii

   !> Invert a cavity's persistent grid numbering into an index lookup
   !>
   !> DROP filters and reorders its grid on every rebuild, so a finite-difference
   !> comparison between two builds cannot be keyed on the array index. It has to
   !> go through `cavity%numbering`, which survives the rebuild. `map(n)` is the
   !> current array index of the point numbered `n`, or 0 if that point is absent
   !> from this build.
   !>
   !> @param[in]  numbering Persistent point numbers, one per current grid point
   !> @param[out] map       Numbering -> index, sized to the largest number seen
   subroutine build_numbering_map(numbering, map)
      !> Persistent point numbers
      integer, intent(in) :: numbering(:)
      !> Numbering -> current index
      integer, allocatable, intent(out) :: map(:)

      !> Grid index and the number it carries
      integer :: igrid, inum
      !> Largest number present
      integer :: max_num

      max_num = 0
      if (size(numbering) > 0) max_num = maxval(numbering)

      allocate (map(max(1, max_num)), source=0)
      do igrid = 1, size(numbering)
         inum = numbering(igrid)
         if (inum > 0 .and. inum <= max_num) map(inum) = igrid
      end do
   end subroutine build_numbering_map

end module test_helpers
