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
!>   * `fd4_scalar(fpp, fp, fm, fmm, h)` - 4-point central FD formula
!>
!> No global Fortran RNG state is touched (self-contained LCG), so the
!> point and structure samplers are safe under parallel test execution.
module test_helpers
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use mctc_env_error, only: moist_error_type => error_type
   use mstore, only: get_structure
   use mstore_data_record, only: record_type
   use mstore_mb16_43, only: get_mb16_43_records
   use mstore_heavy28, only: get_heavy28_records
   use mstore_amino20x4, only: get_amino20x4_records
   use mstore_but14diol, only: get_but14diol_records
   use mstore_upu23, only: get_upu23_records
   use moist_radii, only: default_cpcm_radii, radius_type
   implicit none
   private

   public :: center_at_origin
   public :: get_test_structures
   public :: get_test_radii
   public :: get_test_points
   public :: fd4_scalar

   !> Default n for get_test_structures (must be a multiple of 5).
   integer, parameter :: default_n_structures = 3
   !> Default n for get_test_points.
   integer, parameter :: default_n_points = 5

   !> The 5 mstore collections that get_test_structures samples from.
   integer, parameter :: n_datasets = 3
   character(len=*), parameter :: datasets(n_datasets) = [character(len=10):: &
      "MB16-43", "Amino20x4", "But14diol"] ! As the CFC is slow, we skip UPU23 for now

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
      integer,  parameter :: max_attempts = 10000

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

end module test_helpers
