!> Test suite for batched molecular partition schemes
!>
!> - Partition of unity over all owners
!> - SSF hard zeros and ones on a homonuclear bond axis and near nuclei
!> - Closed-form weights of a homonuclear pair, Becke stiffness, and
!>   heteronuclear midpoint weights with the size adjustment
!> - Closed-form switch values, tails and endpoints; explicit SSF and power widths
!> - Nuclear derivatives: one-sided continuity and the analytic pair switch
!> - Coincident sites share equal weights
!> - Single atom, batch against single points, threaded against serial, and a
!>   call from inside an enclosing parallel region
!> - Every test runs as a selected test, outside test-drive's OpenMP team, so
!>   the routine's own parallel region is active
module test_math_grid_3d_partition
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads, omp_in_parallel, &
!$    & omp_get_thread_num, omp_get_num_procs
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use moist_data_atomicrad, only: covalent_rad
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test, test_failed
   use test_helpers, only: get_test_structures, get_test_points, fd4_scalar
   use moist_math_grid_3d_partition, only: becke_partition_weights, ssf_partition_weights, &
      & pvoronoi_partition_weights
   use moist_math_grid_3d_partition_becke, only: becke_cell
   use moist_math_grid_3d_partition_ssf, only: ssf_cell
   use moist_math_grid_3d_partition_pvoronoi, only: bump_cell
   implicit none(type, external)
   private

   public :: collect_math_grid_3d_partition

   !> SSF cutoff parameter (Stratmann-Scuseria-Frisch recommendation)
   real(wp), parameter :: ssf_a = 0.64_wp
   !> Becke cell with one polynomial iteration (generic iterate path)
   integer, parameter :: cell_becke_k1 = 1
   !> Becke cell with three iterations passed explicitly (unrolled path)
   integer, parameter :: cell_becke_k3 = 2
   !> Becke cell with the stiffness left at its default
   integer, parameter :: cell_becke_default = 3
   !> SSF cell with `a = ssf_a`
   integer, parameter :: cell_ssf = 4
   !> C-infinity power-Voronoi partition
   integer, parameter :: cell_power = 6
   !> Number of structures drawn from `get_test_structures`
   integer, parameter :: nmol = 5
   !> Random box points per structure, on top of the atom-centered points
   integer, parameter :: nbox = 30
   !> Radii of the atom-centered sample shells (bohr)
   real(wp), parameter :: shell_radii(2) = [0.45_wp, 1.6_wp]

contains

   !> Collect all math_grid_3d_partition tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_grid_3d_partition(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("partition_of_unity", test_partition_of_unity), &
                  new_unittest("compact_switches", test_compact_switches), &
                  new_unittest("nuclear_derivative_continuity", test_nuclear_derivative_continuity), &
                  new_unittest("power_geometry", test_power_geometry), &
                  new_unittest("partition_invariance", test_partition_invariance), &
                  new_unittest("large_coincident_partition", test_large_coincident), &
                  new_unittest("ssf_bond_axis_hard_zeros", test_ssf_bond_axis), &
                  new_unittest("ssf_owned_near_nucleus", test_ssf_near_nucleus), &
                  new_unittest("single_atom_weight_one", test_single_atom), &
                  new_unittest("homonuclear_analytic_weights", test_homonuclear_analytic), &
                  new_unittest("batch_matches_single_points", test_batch_matches_single), &
                  new_unittest("threaded_matches_serial", test_threaded_matches_serial), &
                  new_unittest("nested_call_inside_parallel_region", test_nested_call) &
                  ]
   end subroutine collect_math_grid_3d_partition

   !> Batched weights of one owner under a cell setting
   !>
   !> @param[in]  cell     Cell-function setting (`cell_*`)
   !> @param[in]  owner    Atom whose weight is returned
   !> @param[in]  points   Sample points, shape (3, npts)
   !> @param[in]  xyz      Atom positions, shape (3, nat)
   !> @param[in]  numbers  Atomic numbers, shape (nat)
   !> @param[out] w        Owner weights, shape (npts)
   subroutine batched(cell, owner, points, xyz, numbers, w)
      !> Cell-function setting
      integer, intent(in) :: cell
      !> Atom whose weight is returned
      integer, intent(in) :: owner
      !> Sample points
      real(wp), intent(in) :: points(:, :)
      !> Atom positions
      real(wp), intent(in) :: xyz(:, :)
      !> Atomic numbers
      integer, intent(in) :: numbers(:)
      !> Owner weights
      real(wp), intent(out) :: w(:)

      select case (cell)
      case (cell_becke_k1)
         call becke_partition_weights(owner, points, xyz, numbers, w, stiffness=1)
      case (cell_becke_k3)
         call becke_partition_weights(owner, points, xyz, numbers, w, stiffness=3)
      case (cell_becke_default)
         call becke_partition_weights(owner, points, xyz, numbers, w)
      case (cell_power)
         call pvoronoi_partition_weights(owner, points, xyz, numbers, w)
      case default
         call ssf_partition_weights(owner, points, xyz, numbers, w, a=ssf_a)
      end select
   end subroutine batched

   !> Atomic numbers and sample points for one structure
   !>
   !> - Six axis directions at each of `shell_radii` around every atom
   !> - `nbox` deterministic random points in the padded bounding box
   !>
   !> @param[in]  mol      Structure
   !> @param[out] numbers  Atomic numbers, shape (nat)
   !> @param[out] points   Sample points, shape (3, npts)
   subroutine make_fixture(mol, numbers, points)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Atomic numbers
      integer, allocatable, intent(out) :: numbers(:)
      !> Sample points
      real(wp), allocatable, intent(out) :: points(:, :)

      real(wp), allocatable :: box(:, :)
      real(wp) :: dirs(3, 6)
      integer :: iat, ir, id, ip

      dirs = reshape([1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, &
         & 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp], [3, 6])
      allocate (numbers(mol%nat))
      numbers(:) = mol%num(mol%id)
      call get_test_points(mol, box, nbox)
      allocate (points(3, mol%nat*size(shell_radii)*6 + nbox))
      ip = 0
      do iat = 1, mol%nat
         do ir = 1, size(shell_radii)
            do id = 1, 6
               ip = ip + 1
               points(:, ip) = mol%xyz(:, iat) + shell_radii(ir)*dirs(:, id)
            end do
         end do
      end do
      points(:, ip + 1:) = box
   end subroutine make_fixture

   !> Weights of all owners at one point sum to one within round-off
   !>
   !> Bound `4*nat*eps`: one rounding per normalized weight plus the
   !> accumulated error of the normalization sum
   subroutine test_partition_of_unity(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: cells(3) = [cell_becke_k3, cell_ssf, cell_power]
      type(structure_type), allocatable :: mols(:)
      integer, allocatable :: numbers(:)
      real(wp), allocatable :: points(:, :), w(:, :)
      integer :: im, ic, owner, nat
      real(wp) :: bound

      call get_test_structures(mols, nmol)
      do im = 1, size(mols)
         call make_fixture(mols(im), numbers, points)
         nat = mols(im)%nat
         bound = 4.0_wp*real(nat, wp)*epsilon(1.0_wp)
         allocate (w(size(points, 2), nat))
         do ic = 1, size(cells)
            do owner = 1, nat
               call batched(cells(ic), owner, points, mols(im)%xyz, numbers, w(:, owner))
            end do
            call check(error, all(w >= 0.0_wp .and. w <= 1.0_wp), &
               & "partition weight outside [0, 1]")
            if (allocated(error)) return
            call check(error, maxval(abs(sum(w, dim=2) - 1.0_wp)) <= bound, &
               & "partition weights of all owners do not sum to one")
            if (allocated(error)) return
         end do
         deallocate (w)
      end do
   end subroutine test_partition_of_unity

   !> SSF weights on a homonuclear bond axis follow the cell window exactly
   !>
   !> - Equal elements give a zero size adjustment, so on the axis
   !>   `nu = (2*t - R)/R` at distance `t` from atom 1
   !> - Exactly 1 (atom 1) and 0 (atom 2) for `nu <= -a`, the reverse for
   !>   `nu >= a`, strictly inside (0, 1) in between
   !> - Points beyond either nucleus lie outside the window as well
   subroutine test_ssf_bond_axis(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: npts = 81
      real(wp), parameter :: bond = 2.1_wp
      !> Clearance from the window edges; 1 - g(z) ~ (1 - |z|)**4 must stay above round-off
      real(wp), parameter :: margin = 1.0e-3_wp
      real(wp) :: xyz(3, 2), points(3, npts), t(npts), nu, w1(npts), w2(npts)
      integer :: numbers(2), ip

      numbers = [7, 7]
      xyz = 0.0_wp
      xyz(3, 2) = bond
      points = 0.0_wp
      do ip = 1, npts
         t(ip) = -1.0_wp + (bond + 2.0_wp)*real(ip - 1, wp)/real(npts - 1, wp)
         points(3, ip) = t(ip)
      end do
      call becke_partition_weights(1, points, xyz, numbers, w1, ssf_a=ssf_a)
      call becke_partition_weights(2, points, xyz, numbers, w2, ssf_a=ssf_a)

      do ip = 1, npts
         nu = (2.0_wp*min(max(t(ip), 0.0_wp), bond) - bond)/bond
         if (nu <= -ssf_a - margin) then
            call check(error, w1(ip) == 1.0_wp .and. w2(ip) == 0.0_wp, &
               & "SSF weights are not exactly (1, 0) inside the first cell")
         else if (nu >= ssf_a + margin) then
            call check(error, w1(ip) == 0.0_wp .and. w2(ip) == 1.0_wp, &
               & "SSF weights are not exactly (0, 1) inside the second cell")
         else if (abs(nu) <= ssf_a - margin) then
            call check(error, w1(ip) > 0.0_wp .and. w1(ip) < 1.0_wp .and. &
               & w2(ip) > 0.0_wp .and. w2(ip) < 1.0_wp, &
               & "SSF weights inside the switching window are not strictly fractional")
         end if
         if (allocated(error)) return
      end do
   end subroutine test_ssf_bond_axis

   !> Points 0.05 bohr from a nucleus are owned outright under SSF
   subroutine test_ssf_near_nucleus(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: offset = 0.05_wp
      type(structure_type), allocatable :: mols(:)
      integer, allocatable :: numbers(:)
      real(wp), allocatable :: points(:, :), w(:)
      integer :: im, iat, owner, nat

      call get_test_structures(mols, nmol)
      do im = 1, size(mols)
         nat = mols(im)%nat
         allocate (numbers(nat), points(3, 2*nat), w(2*nat))
         numbers(:) = mols(im)%num(mols(im)%id)
         do iat = 1, nat
            points(:, 2*iat - 1) = mols(im)%xyz(:, iat) + [offset, 0.0_wp, 0.0_wp]
            points(:, 2*iat) = mols(im)%xyz(:, iat) - [0.0_wp, offset, offset]/sqrt(2.0_wp)
         end do
         do owner = 1, nat
            call becke_partition_weights(owner, points, mols(im)%xyz, numbers, w, ssf_a=ssf_a)
            do iat = 1, nat
               if (iat == owner) then
                  call check(error, all(w(2*iat - 1:2*iat) == 1.0_wp), &
                     & "SSF weight next to the owner's nucleus is not exactly 1")
               else
                  call check(error, all(w(2*iat - 1:2*iat) == 0.0_wp), &
                     & "SSF weight next to another nucleus is not exactly 0")
               end if
               if (allocated(error)) return
            end do
         end do
         deallocate (numbers, points, w)
      end do
   end subroutine test_ssf_near_nucleus

   !> A single atom owns every point with weight exactly 1
   subroutine test_single_atom(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: cells(5) = [cell_becke_k1, cell_becke_k3, cell_becke_default, cell_ssf, cell_power]
      real(wp) :: xyz(3, 1), points(3, 4), w(4)
      integer :: ic

      xyz(:, 1) = [0.3_wp, -0.2_wp, 1.1_wp]
      points(:, 1) = xyz(:, 1)
      points(:, 2) = xyz(:, 1) + [0.5_wp, 0.0_wp, 0.0_wp]
      points(:, 3) = [-4.0_wp, 7.0_wp, 2.5_wp]
      points(:, 4) = [100.0_wp, 0.0_wp, -50.0_wp]
      do ic = 1, size(cells)
         call batched(cells(ic), 1, points, xyz, [8], w)
         call check(error, all(w == 1.0_wp), "single-atom partition weight is not exactly 1")
         if (allocated(error)) return
      end do
   end subroutine test_single_atom

   !> Closed-form weights of a homonuclear pair
   !>
   !> - Equal radii switch the size adjustment off, so `nu = mu`
   !> - On the bisecting plane `mu = 0`: both atoms get 1/2 for every cell function
   !> - On the axis at `mu = 1/2` with one Becke iteration: `p = 11/16`, so
   !>   the far atom gets 5/32 and the near atom 27/32
   !> - Three iterations, explicit or by default, at the same point
   !> - Heteronuclear midpoint `mu = 0`, so `nu = a`: O-H gives Becke and SSF
   !>   an unclamped size adjustment, H-Cs one clamped to `|a| = 1/2` (Becke
   !>   1988, App. A); both power-cell pairs saturate and pin only the sign of
   !>   the squared-radius offset
   subroutine test_homonuclear_analytic(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: cells(5) = [cell_becke_k1, cell_becke_k3, cell_becke_default, cell_ssf, cell_power]
      integer, parameter :: stiff_cells(2) = [cell_becke_k3, cell_becke_default]
      real(wp), parameter :: tol = 4.0_wp*epsilon(1.0_wp)
      real(wp) :: xyz(3, 2), plane(3, 3), axis(3, 1), w(3), w1(1), w2(1)
      integer :: ic, owner, k, iter, ipair
      integer, parameter :: pairs(2, 2) = reshape([8, 1, 1, 55], [2, 2])
      real(wp) :: chi, u, nu, expected, x

      xyz(:, 1) = [0.0_wp, 0.0_wp, -1.0_wp]
      xyz(:, 2) = [0.0_wp, 0.0_wp, 1.0_wp]
      plane(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      plane(:, 2) = [0.7_wp, -0.4_wp, 0.0_wp]
      plane(:, 3) = [-6.0_wp, 2.5_wp, 0.0_wp]
      do ic = 1, size(cells)
         do owner = 1, 2
            call batched(cells(ic), owner, plane, xyz, [8, 8], w)
            call check(error, all(abs(w - 0.5_wp) <= tol), &
               & "weight on the bisecting plane of a homonuclear pair is not 1/2")
            if (allocated(error)) return
         end do
      end do

      ! r1 = 1.5, r2 = 0.5, R = 2
      axis(:, 1) = [0.0_wp, 0.0_wp, 0.5_wp]
      call batched(cell_becke_k1, 1, axis, xyz, [8, 8], w1)
      call batched(cell_becke_k1, 2, axis, xyz, [8, 8], w2)
      call check(error, abs(w1(1) - 5.0_wp/32.0_wp) <= tol .and. abs(w2(1) - 27.0_wp/32.0_wp) <= tol, &
         & "one-iteration Becke weights at mu = 1/2 deviate from 5/32 and 27/32")
      if (allocated(error)) return
      x = 0.5_wp
      do iter = 1, 3
         x = (3.0_wp*x - x**3)/2.0_wp
      end do
      do ic = 1, size(stiff_cells)
         call batched(stiff_cells(ic), 1, axis, xyz, [8, 8], w1)
         call check(error, abs(w1(1) - (1.0_wp - x)/2.0_wp) < 2.0e-15_wp, &
            & "explicit and default Becke stiffness must follow repeated cubic switches")
         if (allocated(error)) return
      end do

      ! Midpoint weights: unclamped (O-H) and clamped (H-Cs) Becke/SSF size adjustment
      axis = 0.0_wp
      do ipair = 1, size(pairs, 2)
         chi = covalent_rad(pairs(2, ipair))/covalent_rad(pairs(1, ipair))
         u = (chi - 1.0_wp)/(chi + 1.0_wp)
         nu = -max(-0.5_wp, min(0.5_wp, u/(u*u - 1.0_wp)))
         do ic = 1, size(cells)
            x = nu
            select case (cells(ic))
            case (cell_becke_k1, cell_becke_k3, cell_becke_default)
               k = 3
               if (cells(ic) == cell_becke_k1) k = 1
               do iter = 1, k
                  x = (3.0_wp*x - x**3)/2.0_wp
               end do
               expected = (1.0_wp - x)/2.0_wp
            case (cell_ssf)
               x = x/ssf_a
               expected = (1.0_wp - (35.0_wp*x - 35.0_wp*x**3 + &
                  & 21.0_wp*x**5 - 5.0_wp*x**7)/16.0_wp)/2.0_wp
            case (cell_power)
               x = covalent_rad(pairs(2, ipair))**2 - covalent_rad(pairs(1, ipair))**2
               expected = 0.0_wp
               if (abs(x) < 1.0_wp) expected = 1.0_wp/(1.0_wp + exp(2.0_wp*x/(1.0_wp - x*x)))
               if (x <= -1.0_wp) expected = 1.0_wp
            case default
               call check(error, .false., "unknown analytical pair cell")
               return
            end select
            call batched(cells(ic), 1, axis, xyz, pairs(:, ipair), w1)
            call check(error, abs(w1(1) - expected) < 2.0e-15_wp, &
               & "heteronuclear midpoint weights must follow element radii and size correction")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_homonuclear_analytic

   !> A batch gives the same numbers as its points passed one at a time
   subroutine test_batch_matches_single(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: cells(3) = [cell_becke_k3, cell_ssf, cell_power]
      type(structure_type), allocatable :: mols(:)
      integer, allocatable :: numbers(:)
      real(wp), allocatable :: points(:, :), w(:), w1(:)
      integer :: ic, ip, owner

      call get_test_structures(mols, nmol)
      call make_fixture(mols(1), numbers, points)
      allocate (w(size(points, 2)), w1(size(points, 2)))
      do ic = 1, size(cells)
         do owner = 1, mols(1)%nat
            call batched(cells(ic), owner, points, mols(1)%xyz, numbers, w)
            do ip = 1, size(points, 2)
               call batched(cells(ic), owner, points(:, ip:ip), mols(1)%xyz, numbers, w1(ip:ip))
            end do
            call check(error, all(w == w1), "batched weights differ from single-point calls")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_batch_matches_single

   !> Weights computed by a thread team equal the single-thread result
   !>
   !> - Forces a four-thread team even when `OMP_NUM_THREADS=1`, so the
   !>   parallel path is exercised in every environment; skips without
   !>   OpenMP or on a single processor
   !> - Fails inside an enclosing team, where the routine cannot thread
   subroutine test_threaded_matches_serial(error)
      !> Test failure, or set by `skip_test`
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: cells(3) = [cell_becke_k3, cell_ssf, cell_power]
      integer, parameter :: nteam = 4
      type(structure_type), allocatable :: mols(:)
      integer, allocatable :: numbers(:)
      real(wp), allocatable :: points(:, :), w_serial(:), w_team(:)
      integer :: max_threads, nprocs, im, ic, owner
      logical :: has_omp, in_parallel

      has_omp = .false.
      max_threads = 1
      nprocs = 1
      in_parallel = .false.
!$    has_omp = .true.
!$    max_threads = omp_get_max_threads()
!$    nprocs = omp_get_num_procs()
!$    in_parallel = omp_in_parallel()
      if (.not. has_omp .or. nprocs < 2) then
         call skip_test(error, "OpenMP disabled or only one processor available")
         return
      end if
      call check(error, .not. in_parallel, &
         & "threaded partition check must run as a selected test, not inside the suite-level team")
      if (allocated(error)) return

      call get_test_structures(mols, nmol)
      do im = 1, size(mols)
         call make_fixture(mols(im), numbers, points)
         allocate (w_serial(size(points, 2)), w_team(size(points, 2)))
         do ic = 1, size(cells)
            do owner = 1, mols(im)%nat
!$             call omp_set_num_threads(1)
               call batched(cells(ic), owner, points, mols(im)%xyz, numbers, w_serial)
!$             call omp_set_num_threads(nteam)
               call batched(cells(ic), owner, points, mols(im)%xyz, numbers, w_team)
               call check(error, all(w_team == w_serial), &
                  & "threaded partition weights differ from the single-thread result")
               if (allocated(error)) exit
            end do
            if (allocated(error)) exit
         end do
!$       call omp_set_num_threads(max_threads)
         if (allocated(error)) return
         deallocate (w_serial, w_team)
      end do
   end subroutine test_threaded_matches_serial

   !> Calls from inside an enclosing parallel region run and agree
   !>
   !> Each thread of a two-thread team computes the same owner weights into
   !> its own column; every column written equals the result of a call made
   !> outside the region
   !>
   subroutine test_nested_call(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nteam = 2
      type(structure_type), allocatable :: mols(:)
      integer, allocatable :: numbers(:)
      real(wp), allocatable :: points(:, :), w_ref(:), w_team(:, :)
      logical :: ran(nteam)
      integer :: owner, it, ic
      integer, parameter :: cells(3) = [cell_becke_k3, cell_ssf, cell_power]

      call get_test_structures(mols, nmol)
      call make_fixture(mols(2), numbers, points)
      allocate (w_ref(size(points, 2)), w_team(size(points, 2), nteam))
      do ic = 1, size(cells)
         do owner = 1, mols(2)%nat
            call batched(cells(ic), owner, points, mols(2)%xyz, numbers, w_ref)
            ran = .false.
            !$omp parallel num_threads(nteam) default(shared) private(it)
            it = 1
!$          it = omp_get_thread_num() + 1
            ran(it) = .true.
            call batched(cells(ic), owner, points, mols(2)%xyz, numbers, w_team(:, it))
            !$omp end parallel
            do it = 1, nteam
               if (.not. ran(it)) cycle
               call check(error, all(w_team(:, it) == w_ref), &
                  & "partition weights computed inside a parallel region differ")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_nested_call

   !> Published SSF polynomial and bump switch values, tails and exact endpoints
   !>
   !> - Becke switch: positive tail after three iterations, every iteration
   !>   kept beyond the default stiffness
   !> - Bump switch: logistic interior value, complement, super-polynomial tail
   !> - Legacy SSF option and an explicit SSF width `a`
   subroutine test_compact_switches(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: x, expected, h, points(3, 1), xyz(3, 2), legacy(1), direct(1)
      real(wp) :: tail
      integer :: i

      tail = becke_cell(1.0_wp - 1.0e-3_wp, 3)
      call check(error, tail > 0.0_wp .and. tail < 1.0e-20_wp, "Becke tail remains positive after three iterations")
      if (allocated(error)) return
      x = 0.2_wp
      do i = 1, 5
         x = (3.0_wp*x - x**3)/2.0_wp
      end do
      call check(error, abs(becke_cell(0.2_wp, 5) - (1.0_wp - x)/2.0_wp) < 1.0e-15_wp, &
         & "Becke stiffness beyond the default must retain every iteration")
      if (allocated(error)) return
      do i = -10, 10
         x = real(i, wp)/10.0_wp
         expected = (1.0_wp - (35.0_wp*x - 35.0_wp*x**3 + 21.0_wp*x**5 - 5.0_wp*x**7)/16.0_wp)/2.0_wp
         call check(error, abs(ssf_cell(x) - expected) < 8.0_wp*epsilon(x), "SSF polynomial")
         if (allocated(error)) return
         call check(error, abs(bump_cell(x) + bump_cell(-x) - 1.0_wp) < epsilon(x), "bump complement")
         if (allocated(error)) return
      end do
      call check(error, ssf_cell(1.0_wp) == 0.0_wp .and. ssf_cell(-1.0_wp) == 1.0_wp &
         & .and. bump_cell(1.0_wp) == 0.0_wp .and. bump_cell(-1.0_wp) == 1.0_wp, "exact endpoints")
      if (allocated(error)) return
      h = 1.0e-4_wp
      call check(error, ssf_cell(1.0_wp - h) > 0.0_wp .and. ssf_cell(1.0_wp - h) < 3.0_wp*h**4, &
         & "SSF tail remains positive with fourth-order approach to zero")
      if (allocated(error)) return
      ! Relative bound: rounding x*x at x = 0.95 (0.5 ulp, the subtraction is
      ! exact) moves the exponent 2x/(1 - x*x) = 19.5 by up to 1.1e-14, which
      ! carries over to the logistic tail as relative error
      x = 0.95_wp
      expected = 1.0_wp/(1.0_wp + exp(2.0_wp*x/(1.0_wp - x*x)))
      call check(error, abs(bump_cell(x)/expected - 1.0_wp) < 2.0e-14_wp, &
         & "bump switch interior tail must follow its logistic formula")
      if (allocated(error)) return
      h = 0.02_wp
      call check(error, bump_cell(1.0_wp - h) > 0.0_wp .and. bump_cell(1.0_wp - h) < h**10, &
         & "bump tail approaches zero faster than a finite-order switch")
      if (allocated(error)) return
      xyz = 0.0_wp
      xyz(1, 2) = 2.0_wp
      points(:, 1) = [0.7_wp, 0.2_wp, 0.0_wp]
      call becke_partition_weights(1, points, xyz, [8, 1], legacy, stiffness=1, ssf_a=ssf_a)
      call ssf_partition_weights(1, points, xyz, [8, 1], direct)
      call check(error, all(legacy == direct), "legacy SSF option selects the separate SSF routine")
      if (allocated(error)) return
      x = -0.3_wp/0.4_wp
      expected = (1.0_wp - (35.0_wp*x - 35.0_wp*x**3 + 21.0_wp*x**5 - 5.0_wp*x**7)/16.0_wp)/2.0_wp
      points(:, 1) = [0.7_wp, 0.0_wp, 0.0_wp]
      call ssf_partition_weights(1, points, xyz, [1, 1], direct, a=0.4_wp)
      call check(error, abs(direct(1) - expected) < 1.0e-15_wp, "explicit SSF width must set its polynomial argument")
   end subroutine test_compact_switches

   !> Nuclear derivatives agree from both sides at a pair distance, nucleus and SSF boundary
   !>
   !> - Central difference (`fd4_scalar`) of the weight at the off-axis point
   !>   matches the analytic derivative of each pair switch
   !>
   subroutine test_nuclear_derivative_continuity(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: cells(4) = [cell_becke_k1, cell_becke_k3, cell_ssf, cell_power]
      real(wp), parameter :: h = 1.0e-5_wp
      real(wp) :: xyz(3, 2), points(3, 3), w(3, -2:2), left(3), right(3)
      integer :: ic, i, iter, k
      real(wp) :: r1, r2, mu, dmu, x, dx, sw, expected, fd

      points = 0.0_wp
      points(:, 1) = [2.45_wp, 0.3_wp, 0.2_wp]
      points(:, 2) = [5.0_wp, 0.0_wp, 0.0_wp]
      points(:, 3) = [0.9_wp, 0.0_wp, 0.0_wp]
      do ic = 1, size(cells)
         do i = -2, 2
            xyz = 0.0_wp
            xyz(1, 2) = 5.0_wp + real(i, wp)*h
            call batched(cells(ic), 2, points, xyz, [1, 1], w(:, i))
         end do
         left = (3.0_wp*w(:, 0) - 4.0_wp*w(:, -1) + w(:, -2))/(2.0_wp*h)
         right = (-3.0_wp*w(:, 0) + 4.0_wp*w(:, 1) - w(:, 2))/(2.0_wp*h)
         call check(error, maxval(abs(left - right)) < 2.0e-8_wp, &
            & "partition nuclear derivative must be continuous")
         if (allocated(error)) return
         r1 = norm2(points(:, 1))
         r2 = norm2(points(:, 1) - [5.0_wp, 0.0_wp, 0.0_wp])
         mu = (r2 - r1)/5.0_wp
         dmu = ((5.0_wp - points(1, 1))/r2 - mu)/5.0_wp
         select case (cells(ic))
         case (cell_becke_k1, cell_becke_k3)
            x = mu
            dx = dmu
            k = 3
            if (cells(ic) == cell_becke_k1) k = 1
            do iter = 1, k
               dx = 1.5_wp*(1.0_wp - x*x)*dx
               x = (3.0_wp*x - x**3)/2.0_wp
            end do
            expected = -0.5_wp*dx
         case (cell_ssf)
            x = mu/ssf_a
            expected = -35.0_wp/32.0_wp*(1.0_wp - x*x)**3*dmu/ssf_a
         case (cell_power)
            x = (25.0_wp - 10.0_wp*points(1, 1))
            dx = 2.0_wp*(5.0_wp - points(1, 1))
            sw = 1.0_wp/(1.0_wp + exp(2.0_wp*x/(1.0_wp - x*x)))
            expected = -2.0_wp*(1.0_wp + x*x)/(1.0_wp - x*x)**2*sw*(1.0_wp - sw)*dx
         case default
            call check(error, .false., "unknown analytical derivative cell")
            return
         end select
         if (.not. ieee_is_finite(expected)) then
            call test_failed(error, "analytic pair-switch derivative is not finite")
            return
         end if
         call fd4_scalar(w(1, 2), w(1, 1), w(1, -1), w(1, -2), h, fd, error)
         if (allocated(error)) return
         call check(error, fd, expected, thr=1.0e-9_wp, &
            & more="nuclear derivative must match the analytic pair switch")
         if (allocated(error)) return
      end do
   end subroutine test_nuclear_derivative_continuity

   !> Power-radius boundary shift, bump values, all-space coverage and regularity at a nucleus
   subroutine test_power_geometry(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: xyz(3, 2), points(3, 5), w(5), w2(5), expected, h, derivative, fd

      xyz = 0.0_wp
      xyz(1, 2) = 2.0_wp
      points = 0.0_wp
      points(1, :) = [-1.0e6_wp, 0.25_wp, 0.375_wp, 0.5_wp, 1.0e6_wp]
      call pvoronoi_partition_weights(1, points, xyz, [1, 1], w, radii=[1.0_wp, 2.0_wp])
      call pvoronoi_partition_weights(2, points, xyz, [1, 1], w2, radii=[1.0_wp, 2.0_wp])
      expected = 1.0_wp/(1.0_wp + exp(4.0_wp/3.0_wp))
      call check(error, w(1) == 1.0_wp .and. w(2) == 0.5_wp .and. abs(w(3) - expected) < 1.0e-15_wp &
         & .and. all(w(4:5) == 0.0_wp) .and. maxval(abs(w + w2 - 1.0_wp)) < 1.0e-15_wp, &
         & "power cells must shift with squared radii, have exact zeros, and cover the far field")
      if (allocated(error)) return

      h = 1.0e-4_wp
      points(1, :) = [-2.0_wp*h, -h, 0.0_wp, h, 2.0_wp*h]
      call pvoronoi_partition_weights(1, points, xyz, [1, 1], w, width=8.0_wp, radii=[1.0_wp, 1.0_wp])
      expected = 1.0_wp/(1.0_wp + exp(-4.0_wp/3.0_wp))
      call check(error, abs(w(3) - expected) < 1.0e-15_wp, "explicit power width must set the transition scale")
      if (allocated(error)) return
      derivative = -0.5_wp*(40.0_wp/9.0_wp)*expected*(1.0_wp - expected)
      call fd4_scalar(w(5), w(4), w(2), w(1), h, fd, error)
      if (allocated(error)) return
      call check(error, fd, derivative, thr=1.0e-11_wp, &
         & more="power partition derivative through a nucleus must follow the smooth bump")
   end subroutine test_power_geometry

   !> Rotation, translation and relabeling preserve all partition schemes
   subroutine test_partition_invariance(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: cells(3) = [cell_becke_k3, cell_ssf, cell_power]
      integer, parameter :: perm(3) = [3, 1, 2], numbers(3) = [8, 1, 6]
      real(wp) :: xyz(3, 3), points(3, 4), moved_xyz(3, 3), moved_points(3, 4), w(4), ref(4)
      real(wp), parameter :: shift(3) = [3.0_wp, -4.0_wp, 2.0_wp]
      integer :: ic, i

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.5_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [0.0_wp, 2.0_wp, 0.0_wp]
      points(:, :3) = xyz
      points(:, 4) = [0.7_wp, 1.1_wp, 0.4_wp]
      do i = 1, 3
         moved_xyz(:, i) = xyz([2, 3, 1], perm(i)) + shift
      end do
      do i = 1, 4
         moved_points(:, i) = points([2, 3, 1], i) + shift
      end do
      do ic = 1, size(cells)
         call batched(cells(ic), 1, points, xyz, numbers, ref)
         call batched(cells(ic), 2, moved_points, moved_xyz, numbers(perm), w)
         call check(error, maxval(abs(w - ref)) < 1.0e-14_wp, "partition changed under rigid motion or relabeling")
         if (allocated(error)) return
      end do
   end subroutine test_partition_invariance

   !> Log normalization retains equal weights when every direct cell product underflows
   subroutine test_large_coincident(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nat = 1100
      real(wp) :: xyz(3, nat), points(3, 1), w(1)
      integer :: numbers(nat)

      xyz = 0.0_wp
      points = 0.0_wp
      numbers = 1
      call becke_partition_weights(1, points, xyz(:, :3), numbers(:3), w)
      call check(error, abs(w(1) - 1.0_wp/3.0_wp) < epsilon(1.0_wp), &
         & "coincident Becke sites share equal weights")
      if (allocated(error)) return
      call ssf_partition_weights(1, points, xyz(:, :3), numbers(:3), w)
      call check(error, abs(w(1) - 1.0_wp/3.0_wp) < epsilon(1.0_wp), &
         & "coincident SSF sites share equal weights")
      if (allocated(error)) return
      call pvoronoi_partition_weights(1, points, xyz, numbers, w)
      call check(error, abs(w(1)*real(nat, wp) - 1.0_wp) < 4.0_wp*epsilon(1.0_wp), &
         & "underflowing coincident power-cell products must still normalize")
   end subroutine test_large_coincident

end module test_math_grid_3d_partition
