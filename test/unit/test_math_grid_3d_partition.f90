!> Test suite for molecular partition schemes and the molecular grid they assemble
!>
!> Partition weights (batched owner weights of every scheme)
!>
!> - Partition of unity over all owners
!> - SSF hard zeros and ones on a homonuclear bond axis and near nuclei
!> - Closed-form weights of a homonuclear pair, Becke stiffness, and
!>   heteronuclear midpoint weights with the size adjustment
!> - Closed-form switch values, tails and endpoints; explicit SSF and power widths
!> - Nuclear derivatives: one-sided continuity and the analytic pair switch
!> - Coincident sites share equal weights
!>
!> Molecular grid assembly
!>
!> - Every scheme selects, prunes and repartitions the molecular grid
!> - Updates translate fixed local grids, as an independent reconstruction
!> - Atom-count changes match a fresh grid
!> - Element overrides assign each atom its local grid; sector policies
!>   report the largest retained angular rule
!>
!> Owner-level partition derivatives
!>
!> - Switch values, first and second derivatives against finite differences
!> - Owner gradients (nuclear and point) against finite differences
!> - Partition-unity gradients and Hessian-vector products cancel at the
!>   tiny-tail scale
!>
!> Threaded weights are compared in `math_grid_3d_threaded`
module test_math_grid_3d_partition
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use moist_data_atomicrad, only: covalent_rad
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: get_test_structures, get_test_points, fd4_scalar, check_moist_error, &
      & get_becke_recipe
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_shell_sector_type, new_sector_shell_policy, &
      & moist_math_grid_atomic_recipe_type, moist_math_grid_atomic_recipe_override_type, &
      & element_override_index
   use moist_math_grid_atomic_grid, only: moist_math_grid_atomic_type, new_atomic_grid
   use moist_math_grid_3d_kernel_base, only: moist_math_grid_3d_partition_type
   use moist_math_grid_3d_kernel_becke, only: becke_partition_type
   use moist_math_grid_3d_kernel_ssf, only: ssf_partition_type
   use moist_math_grid_3d_kernel_pvoronoi, only: pvoronoi_partition_type
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid, &
      & partition_becke, partition_ssf, partition_pvoronoi
   implicit none(type, external)
   private

   public :: collect_math_grid_3d_partition

   !> SSF cutoff parameter (Stratmann-Scuseria-Frisch recommendation)
   real(wp), parameter :: ssf_a = 0.64_wp
   !> Becke cell with one polynomial iteration (generic iterate path)
   integer, parameter :: cell_becke_k1 = 1
   !> Becke cell with three iterations passed explicitly (unrolled path)
   integer, parameter :: cell_becke_k3 = 2
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
   !> Weight-pruning threshold of the molecular grid (bohr^3)
   real(wp), parameter :: wthr = 1.0e-14_wp

contains

   !> Collect all math_grid_3d_partition tests
   !>
   !> @param[out] testsuite  collected unit tests
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
                  new_unittest("homonuclear_analytic_weights", test_homonuclear_analytic), &
                  new_unittest("partition_assembly", test_partition_assembly), &
                  new_unittest("fixed_local_quadrature", test_fixed_local_quadrature), &
                  new_unittest("update_changes_atom_count", test_update_atom_count), &
                  new_unittest("element_overrides", test_element_overrides), &
                  new_unittest("switch_jets", test_switch_jets), &
                  new_unittest("partition_finite_difference", test_partition_fd), &
                  new_unittest("gradient_tails", test_gradient_tails), &
                  new_unittest("hessian_vector_tails", test_hessian_vector_tails) &
                  ]
   end subroutine collect_math_grid_3d_partition

   !* ------------------------------- Partition weights ------------------------------- *!

   !> Batched weights of one owner under a cell setting
   !>
   !> @param[in]  cell     cell-function setting (`cell_*`)
   !> @param[in]  owner    atom whose weight is returned
   !> @param[in]  points   sample points, shape (3, npts)
   !> @param[in]  xyz      atom positions, shape (3, nat)
   !> @param[in]  numbers  atomic numbers, shape (nat)
   !> @param[out] w        owner weights, shape (npts)
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

      type(becke_partition_type) :: becke
      type(ssf_partition_type) :: ssf
      type(pvoronoi_partition_type) :: pvoronoi

      select case (cell)
      case (cell_becke_k1)
         becke%k = 1
         call becke%owner_weights(owner, points, xyz, numbers, w)
      case (cell_becke_k3)
         becke%k = 3
         call becke%owner_weights(owner, points, xyz, numbers, w)
      case (cell_power)
         call pvoronoi%owner_weights(owner, points, xyz, numbers, w)
      case default
         ssf%a = ssf_a
         call ssf%owner_weights(owner, points, xyz, numbers, w)
      end select
   end subroutine batched

   !> Atomic numbers and sample points for one structure
   !>
   !> - Six axis directions at each of `shell_radii` around every atom
   !> - `nbox` deterministic random points in the padded bounding box
   !>
   !> @param[in]  mol      structure
   !> @param[out] numbers  atomic numbers, shape (nat)
   !> @param[out] points   sample points, shape (3, npts)
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
      integer :: im, ic, owner, nat, ip
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
            do ip = 1, size(w, 1)
               call check(error, sum(w(ip, :)), 1.0_wp, thr=bound, &
                  & more="partition weights of all owners do not sum to one")
               if (allocated(error)) return
            end do
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
      !> Clearance from the window edges
      !>
      !> - 1 - g(z) ~ (1 - |z|)**4 above round-off
      real(wp), parameter :: margin = 1.0e-3_wp
      real(wp) :: xyz(3, 2), points(3, npts), t(npts), nu, w1(npts), w2(npts)
      integer :: numbers(2), ip
      type(ssf_partition_type) :: ssf

      numbers = [7, 7]
      xyz = 0.0_wp
      xyz(3, 2) = bond
      points = 0.0_wp
      do ip = 1, npts
         t(ip) = -1.0_wp + (bond + 2.0_wp)*real(ip - 1, wp)/real(npts - 1, wp)
         points(3, ip) = t(ip)
      end do
      ssf%a = ssf_a
      call ssf%owner_weights(1, points, xyz, numbers, w1)
      call ssf%owner_weights(2, points, xyz, numbers, w2)

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
      type(ssf_partition_type) :: ssf

      ssf%a = ssf_a
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
            call ssf%owner_weights(owner, points, mols(im)%xyz, numbers, w)
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

   !> Closed-form weights of a homonuclear pair
   !>
   !> - Equal radii switch the size adjustment off, so `nu = mu`
   !> - On the bisecting plane `mu = 0`: both atoms get 1/2 for every cell function
   !> - On the axis at `mu = 1/2` with one Becke iteration: `p = 11/16`, so
   !>   the far atom gets 5/32 and the near atom 27/32
   !> - Three iterations at the same point
   !> - Heteronuclear midpoint `mu = 0`, so `nu = a`: O-H gives Becke and SSF
   !>   an unclamped size adjustment, H-Cs one clamped to `|a| = 1/2` (Becke
   !>   1988, Appendix A); both power-cell pairs saturate and pin only the sign of
   !>   the squared-radius offset
   subroutine test_homonuclear_analytic(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: cells(4) = [cell_becke_k1, cell_becke_k3, cell_ssf, cell_power]
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
      call batched(cell_becke_k3, 1, axis, xyz, [8, 8], w1)
      call check(error, abs(w1(1) - (1.0_wp - x)/2.0_wp) < 2.0e-15_wp, &
         & "three-iteration Becke weight must follow repeated cubic switches")
      if (allocated(error)) return

      ! Midpoint weights: unclamped (O-H) and clamped (H-Cs) Becke/SSF size adjustment
      axis = 0.0_wp
      do ipair = 1, size(pairs, 2)
         chi = covalent_rad(pairs(2, ipair))/covalent_rad(pairs(1, ipair))
         u = (chi - 1.0_wp)/(chi + 1.0_wp)
         nu = -max(-0.5_wp, min(0.5_wp, u/(u*u - 1.0_wp)))
         do ic = 1, size(cells)
            x = nu
            select case (cells(ic))
            case (cell_becke_k1, cell_becke_k3)
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

   !> Published SSF polynomial and bump switch values, tails and exact endpoints
   !>
   !> - Becke switch: positive tail after three iterations, every iteration
   !>   kept beyond the default stiffness
   !> - Bump switch: logistic interior value, complement, super-polynomial tail
   !> - Legacy SSF option and an explicit SSF width `a`
   subroutine test_compact_switches(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: x, expected, h, points(3, 1), xyz(3, 2), direct(1)
      real(wp) :: tail
      integer :: i
      type(becke_partition_type) :: becke
      type(ssf_partition_type) :: ssf
      type(pvoronoi_partition_type) :: bump

      ssf%a = 1.0_wp
      becke%k = 3
      tail = becke%eval(1.0_wp - 1.0e-3_wp)
      call check(error, tail > 0.0_wp .and. tail < 1.0e-20_wp, "Becke tail remains positive after three iterations")
      if (allocated(error)) return
      x = 0.2_wp
      do i = 1, 5
         x = (3.0_wp*x - x**3)/2.0_wp
      end do
      becke%k = 5
      call check(error, abs(becke%eval(0.2_wp) - (1.0_wp - x)/2.0_wp) < 1.0e-15_wp, &
         & "Becke stiffness beyond the default must retain every iteration")
      if (allocated(error)) return
      do i = -10, 10
         x = real(i, wp)/10.0_wp
         expected = (1.0_wp - (35.0_wp*x - 35.0_wp*x**3 + 21.0_wp*x**5 - 5.0_wp*x**7)/16.0_wp)/2.0_wp
         call check(error, abs(ssf%eval(x) - expected) < 8.0_wp*epsilon(x), "SSF polynomial")
         if (allocated(error)) return
         call check(error, abs(bump%eval(x) + bump%eval(-x) - 1.0_wp) < epsilon(x), "bump complement")
         if (allocated(error)) return
      end do
      call check(error, ssf%eval(1.0_wp) == 0.0_wp .and. ssf%eval(-1.0_wp) == 1.0_wp &
         & .and. bump%eval(1.0_wp) == 0.0_wp .and. bump%eval(-1.0_wp) == 1.0_wp, "exact endpoints")
      if (allocated(error)) return
      h = 1.0e-4_wp
      call check(error, ssf%eval(1.0_wp - h) > 0.0_wp .and. ssf%eval(1.0_wp - h) < 3.0_wp*h**4, &
         & "SSF tail remains positive with fourth-order approach to zero")
      if (allocated(error)) return
      ! Relative bound: rounding x*x at x = 0.95 (0.5 ulp, the subtraction is
      ! exact) moves the exponent 2x/(1 - x*x) = 19.5 by up to 1.1e-14, which
      ! carries over to the logistic tail as relative error
      x = 0.95_wp
      expected = 1.0_wp/(1.0_wp + exp(2.0_wp*x/(1.0_wp - x*x)))
      call check(error, abs(bump%eval(x)/expected - 1.0_wp) < 2.0e-14_wp, &
         & "bump switch interior tail must follow its logistic formula")
      if (allocated(error)) return
      h = 0.02_wp
      call check(error, bump%eval(1.0_wp - h) > 0.0_wp .and. bump%eval(1.0_wp - h) < h**10, &
         & "bump tail approaches zero faster than a finite-order switch")
      if (allocated(error)) return
      xyz = 0.0_wp
      xyz(1, 2) = 2.0_wp
      x = -0.3_wp/0.4_wp
      expected = (1.0_wp - (35.0_wp*x - 35.0_wp*x**3 + 21.0_wp*x**5 - 5.0_wp*x**7)/16.0_wp)/2.0_wp
      points(:, 1) = [0.7_wp, 0.0_wp, 0.0_wp]
      ssf%a = 0.4_wp
      call ssf%owner_weights(1, points, xyz, [1, 1], direct)
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
         do i = 1, size(left)
            call check(error, ieee_is_finite(right(i)), "right one-sided derivative is not finite")
            if (allocated(error)) return
            call check(error, left(i), right(i), thr=2.0e-8_wp, more="partition nuclear derivative must be continuous")
            if (allocated(error)) return
         end do
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
      real(wp) :: far(3, 1), w_far(1)
      type(pvoronoi_partition_type) :: pvoronoi
      integer :: i

      xyz = 0.0_wp
      xyz(1, 2) = 2.0_wp
      points = 0.0_wp
      points(1, :) = [-1.0e6_wp, 0.25_wp, 0.375_wp, 0.5_wp, 1.0e6_wp]
      pvoronoi%radii = [1.0_wp, 2.0_wp]
      call pvoronoi%owner_weights(1, points, xyz, [1, 1], w)
      call pvoronoi%owner_weights(2, points, xyz, [1, 1], w2)
      expected = 1.0_wp/(1.0_wp + exp(4.0_wp/3.0_wp))
      call check(error, w(1) == 1.0_wp .and. w(2) == 0.5_wp .and. abs(w(3) - expected) < 1.0e-15_wp &
         & .and. all(w(4:5) == 0.0_wp), &
         & "power cells must shift with squared radii, have exact zeros, and cover the far field")
      if (allocated(error)) return
      do i = 1, size(w)
         call check(error, w(i) + w2(i), 1.0_wp, thr=1.0e-15_wp, &
            & more="power cells must shift with squared radii, have exact zeros, and cover the far field")
         if (allocated(error)) return
      end do
      ! Far transverse offset: subtracting squared distances cancels to zero there
      far(:, 1) = [0.375_wp, 0.0_wp, 1.0e12_wp]
      call pvoronoi%owner_weights(1, far, xyz, [1, 1], w_far)
      call check(error, abs(w_far(1) - expected) < 1.0e-15_wp, &
         & "power difference must stay exact far off the bond axis")
      if (allocated(error)) return

      h = 1.0e-4_wp
      points(1, :) = [-2.0_wp*h, -h, 0.0_wp, h, 2.0_wp*h]
      pvoronoi%width = 8.0_wp
      pvoronoi%radii = [1.0_wp, 1.0_wp]
      call pvoronoi%owner_weights(1, points, xyz, [1, 1], w)
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
         do i = 1, size(w)
            call check(error, ieee_is_finite(ref(i)), "reference partition weight is not finite")
            if (allocated(error)) return
            call check(error, w(i), ref(i), thr=1.0e-14_wp, more="partition changed under rigid motion or relabeling")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_partition_invariance

   !> Log normalization retains equal weights when every direct cell product underflows
   subroutine test_large_coincident(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nat = 1100
      real(wp) :: xyz(3, nat), points(3, 1), w(1)
      integer :: numbers(nat)
      type(becke_partition_type) :: becke
      type(ssf_partition_type) :: ssf
      type(pvoronoi_partition_type) :: pvoronoi

      xyz = 0.0_wp
      points = 0.0_wp
      numbers = 1
      call becke%owner_weights(1, points, xyz(:, :3), numbers(:3), w)
      call check(error, abs(w(1) - 1.0_wp/3.0_wp) < epsilon(1.0_wp), &
         & "coincident Becke sites share equal weights")
      if (allocated(error)) return
      call ssf%owner_weights(1, points, xyz(:, :3), numbers(:3), w)
      call check(error, abs(w(1) - 1.0_wp/3.0_wp) < epsilon(1.0_wp), &
         & "coincident SSF sites share equal weights")
      if (allocated(error)) return
      call pvoronoi%owner_weights(1, points, xyz, numbers, w)
      call check(error, abs(w(1)*real(nat, wp) - 1.0_wp) < 4.0_wp*epsilon(1.0_wp), &
         & "underflowing coincident power-cell products must still normalize")
   end subroutine test_large_coincident

   !* ---------------------------- Molecular grid assembly ---------------------------- *!

   !> All schemes select, prune and repartition the molecular grid correctly
   subroutine test_partition_assembly(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_type) :: local(2)
      type(moist_math_grid_3d_molecular_type) :: grid, pruned
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      class(moist_math_grid_3d_partition_type), allocatable :: partition
      real(wp), allocatable :: points(:, :), bw(:)
      real(wp) :: weight
      integer :: scheme, ischeme, step, iat, i, j, zeros
      integer, parameter :: schemes(3) = [partition_becke, partition_ssf, partition_pvoronoi]
      integer, parameter :: numbers(3) = [8, 1, 1]

      call get_becke_recipe(recipe, 14, 0.5_wp, 11, merr, rcut_upper=5.0_wp)
      if (.not. allocated(merr)) call new_atomic_grid(local(1), recipe, 8, merr)
      if (.not. allocated(merr)) call new_atomic_grid(local(2), recipe, 1, merr)
      call check_moist_error(error, merr, "local grids")
      if (allocated(error)) return
      do ischeme = 1, size(schemes)
         scheme = schemes(ischeme)
         if (allocated(partition)) deallocate (partition)
         select case (scheme)
         case (partition_becke)
            allocate (becke_partition_type :: partition)
         case (partition_ssf)
            allocate (ssf_partition_type :: partition)
         case default
            allocate (pvoronoi_partition_type :: partition)
         end select
         call make_water(mol)
         call new_molecular_point_grid(grid, merr, recipe=recipe, reciprocal=.false., partition=scheme, weight_threshold=0.0_wp)
         call check_moist_error(error, merr, "partition construction")
         if (allocated(error)) return
         do step = 1, 2
            call grid%update(mol, merr)
            call check_moist_error(error, merr, "partition update")
            if (allocated(error)) return
            j = 0
            zeros = 0
            do iat = 1, 3
               associate (atom => local(merge(1, 2, iat == 1)))
                  allocate (points(3, atom%npts), bw(atom%npts))
                  do i = 1, atom%npts
                     points(:, i) = atom%shell_r(atom%shell(i))*atom%u(:, i) + mol%xyz(:, iat)
                  end do
                  call partition%owner_weights(iat, points, mol%xyz, numbers, bw)
                  do i = 1, atom%npts
                     weight = atom%w(i)*bw(i)
                     if (weight == 0.0_wp) then
                        zeros = zeros + 1
                        cycle
                     end if
                     j = j + 1
                     call check(error, j <= grid%ngrid, "nonzero point must be retained")
                     if (allocated(error)) return
                     call check(error, grid%owner(j) == iat .and. all(grid%xyz(:, j) == points(:, i)) &
                        & .and. abs(grid%w(j) - weight) <= 1.0e-14_wp*abs(weight), "partitioned point or weight mismatch")
                     if (allocated(error)) return
                  end do
                  deallocate (points, bw)
               end associate
            end do
            call check(error, j == grid%ngrid, "exact-zero compaction must preserve every nonzero point")
            if (allocated(error)) return
            if (scheme == partition_ssf .or. scheme == partition_pvoronoi) then
               call check(error, zeros > 0, "compact partition must remove redundant points")
               if (allocated(error)) return
            end if
            if (step == 1) mol%xyz(:, 2) = mol%xyz(:, 2) + [7.0_wp, 0.35_wp, -0.2_wp]
         end do
      end do
      ! A nonzero threshold remains an explicit, separate pruning option
      call new_molecular_point_grid(pruned, merr, recipe=recipe, reciprocal=.false., &
         & partition=partition_pvoronoi, weight_threshold=1.0e-6_wp)
      if (.not. allocated(merr)) call pruned%update(mol, merr)
      call check_moist_error(error, merr, "threshold pruning")
      if (allocated(error)) return
      call check(error, pruned%ngrid < grid%ngrid, "positive threshold must prune additional nonzero weights")
   end subroutine test_partition_assembly

   !> Updates translate fixed local grids, repartition, and prune, as an independent reconstruction
   !>
   !> Reference per atom: the element's atomic grid, points `r*u + R`,
   !> batched partition weights, weight `w_local*partition`, kept if
   !> `|w| >= 1e-14` and the partition weight reaches the threshold; at two
   !> geometries, the second crossing the threshold
   subroutine test_fixed_local_quadrature(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: threshold = 0.1_wp
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_type) :: local(2)
      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      type(becke_partition_type) :: becke
      real(wp), allocatable :: points(:, :), bw(:)
      real(wp) :: x(3), w, scale
      integer :: step, iat, i, j, k, count_before, numbers(3)
      character(len=120) :: msg

      call make_water(mol)
      numbers = [8, 1, 1]
      call get_becke_recipe(recipe, 14, 0.5_wp, 11, merr)
      if (.not. allocated(merr)) call new_atomic_grid(local(1), recipe, 8, merr)
      if (.not. allocated(merr)) call new_atomic_grid(local(2), recipe, 1, merr)
      if (.not. allocated(merr)) call new_molecular_point_grid(grid, merr, recipe=recipe, becke_k=2, &
         & pruning_threshold=threshold, reciprocal=.false.)
      call check_moist_error(error, merr, "construction")
      if (allocated(error)) return

      count_before = 0
      do step = 1, 2
         call grid%update(mol, merr)
         call check_moist_error(error, merr, "update")
         if (allocated(error)) return
         j = 0
         do iat = 1, mol%nat
            associate (atom => local(merge(1, 2, iat == 1)))
               allocate (points(3, atom%npts), bw(atom%npts))
               do i = 1, atom%npts
                  do k = 1, 3
                     points(k, i) = atom%shell_r(atom%shell(i))*atom%u(k, i) + mol%xyz(k, iat)
                  end do
               end do
               becke%k = 2
               call becke%owner_weights(iat, points, mol%xyz, numbers, bw)
               do i = 1, atom%npts
                  w = atom%w(i)*bw(i)
                  if (abs(w) < wthr .or. bw(i) < threshold) cycle
                  j = j + 1
                  if (j > grid%ngrid) exit
                  x = points(:, i)
                  write (msg, "(a,i0,a,i0,a,i0)") "step ", step, ", atom ", iat, ", local point ", i
                  call check(error, grid%owner(j) == iat .and. abs(grid%w(j) - w) <= 1.0e-14_wp*abs(w), trim(msg))
                  if (allocated(error)) return
                  scale = 1.0_wp + maxval(abs(x))
                  do k = 1, 3
                     call check(error, ieee_is_finite(x(k)), "reference point is not finite")
                     if (allocated(error)) return
                     call check(error, grid%xyz(k, j), x(k), thr=4.0_wp*epsilon(1.0_wp)*scale, more=trim(msg))
                     if (allocated(error)) return
                  end do
               end do
               deallocate (points, bw)
            end associate
         end do
         call check(error, j == grid%ngrid, "retained points differ from the reconstruction")
         if (allocated(error)) return
         if (step == 1) then
            count_before = grid%ngrid
            mol%xyz(:, 2) = mol%xyz(:, 2) + [0.0_wp, 0.35_wp, -0.2_wp]
         else
            call check(error, grid%ngrid /= count_before, "deformation must cross the pruning threshold")
            if (allocated(error)) return
         end if
      end do
   end subroutine test_fixed_local_quadrature

   !> An update to a molecule with a different atom count equals a fresh grid
   !>
   !> H2 to water to H2 on one grid, each against a grid constructed and
   !> updated once; the same code builds both, so the results must be equal
   subroutine test_update_atom_count(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid, fresh
      type(structure_type) :: mols(3)
      type(mctc_error), allocatable :: merr
      character(len=32) :: label
      integer :: step

      call make_h2(mols(1))
      call make_water(mols(2))
      call make_h2(mols(3))
      call get_becke_recipe(recipe, 12, 0.5_wp, 7, merr, rcut_upper=5.0_wp)
      if (.not. allocated(merr)) call new_molecular_point_grid(grid, merr, recipe=recipe, reciprocal=.false.)
      call check_moist_error(error, merr, "construction")
      if (allocated(error)) return

      do step = 1, size(mols)
         write (label, "(a,i0)") "step ", step
         call grid%update(mols(step), merr)
         if (.not. allocated(merr)) call new_molecular_point_grid(fresh, merr, recipe=recipe, reciprocal=.false.)
         if (.not. allocated(merr)) call fresh%update(mols(step), merr)
         call check_moist_error(error, merr, trim(label))
         if (allocated(error)) return
         call check(error, grid%ngrid == fresh%ngrid, trim(label)//": sizes")
         if (allocated(error)) return
         call check(error, all(grid%owner == fresh%owner), trim(label)//": ownership")
         if (allocated(error)) return
         call check(error, all(grid%xyz == fresh%xyz) .and. all(grid%w == fresh%w), &
            & trim(label)//": points and weights")
         if (allocated(error)) return
      end do
   end subroutine test_update_atom_count

   !> Per-element overrides select the recipe per atom; per-shell results follow it
   !>
   !> - First override listing an element wins
   !> - Requested shells reported before the cutoff, retained shells in the CSR arrays
   !> - A variable shell policy reports its largest angular rule per atom
   subroutine test_element_overrides(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      type(moist_math_grid_atomic_type) :: local
      type(moist_math_grid_atomic_shell_sector_type) :: sectors
      type(moist_math_grid_3d_molecular_type) :: grid
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      integer :: iat, idx, lo, hi, z, order

      call make_water(mol)
      call get_becke_recipe(recipe, 20, 0.5_wp, 11, merr)
      allocate (overrides(2))
      overrides(1)%elements = [2, 1]
      if (.not. allocated(merr)) call get_becke_recipe(overrides(1)%recipe, 16, 1.0_wp, 5, merr, rcut_upper=4.0_wp)
      overrides(2)%elements = [1, 8]
      if (.not. allocated(merr)) call get_becke_recipe(overrides(2)%recipe, 30, 0.5_wp, 29, merr)
      if (.not. allocated(merr)) call new_molecular_point_grid(grid, merr, recipe=recipe, overrides=overrides, &
         & reciprocal=.false.)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      call check_moist_error(error, merr, "overrides")
      if (allocated(error)) return
      call check(error, all(grid%nrad_per_atom == [30, 16, 16]) .and. all(grid%nang_per_atom == [302, 14, 14]), &
         & "per-atom sizes do not follow the first matching override")
      if (allocated(error)) return

      do iat = 1, mol%nat
         z = mol%num(mol%id(iat))
         idx = element_override_index(z, overrides)
         call new_atomic_grid(local, overrides(idx)%recipe, z, merr)
         call check_moist_error(error, merr, "atomic grid")
         if (allocated(error)) return
         lo = grid%atom_shell_offset(iat)
         hi = grid%atom_shell_offset(iat + 1) - 1
         call check(error, hi - lo + 1 == local%nshell .and. all(grid%shell_r(lo:hi) == local%shell_r) &
            & .and. all(grid%shell_npts(lo:hi) == local%shell_npts) &
            & .and. all(grid%shell_degree(lo:hi) == local%shell_degree) &
            & .and. all(grid%shell_spacing(lo:hi) == local%shell_spacing), "per-shell results")
         if (allocated(error)) return
         call check(error, grid%atom_offset(iat + 1) - grid%atom_offset(iat) <= local%npts &
            & .and. all(grid%owner(grid%atom_offset(iat):grid%atom_offset(iat + 1) - 1) == iat), &
            & "atom points")
         if (allocated(error)) return
      end do
      ! Fewer retained hydrogen shells than requested after cutoff
      call check(error, grid%atom_shell_offset(3) - grid%atom_shell_offset(2) < grid%nrad_per_atom(2) &
         & .and. all(grid%shell_r(grid%atom_shell_offset(2):grid%atom_shell_offset(3) - 1) <= 4.0_wp), &
         & "rcut_upper of the hydrogen override")
      if (allocated(error)) return

      ! Largest retained angular rule for a variable shell policy
      ! Both sector orders: neither first nor last shell sufficient for the maximum
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      do order = 1, 2
         call get_becke_recipe(recipe, 14, 0.5_wp, 5, merr, rcut_upper=5.0_wp)
         if (.not. allocated(merr)) call new_sector_shell_policy(sectors, [1.0_wp], &
            & merge([5, 17], [17, 5], order == 1), merr)
         call check_moist_error(error, merr, "sector recipe")
         if (allocated(error)) return
         deallocate (recipe%shells)
         allocate (recipe%shells, source=sectors)
         call new_atomic_grid(local, recipe, 1, merr)
         if (.not. allocated(merr)) call new_molecular_point_grid(grid, merr, recipe=recipe, reciprocal=.false.)
         if (.not. allocated(merr)) call grid%update(mol, merr)
         call check_moist_error(error, merr, "variable angular shells")
         if (allocated(error)) return
         ! Lebedev orders 5 and 17 are the 14- and 110-point rules
         call check(error, minval(local%shell_npts) == 14 .and. maxval(local%shell_npts) == 110, &
            & "sector fixture must mix the 14- and 110-point rules")
         if (allocated(error)) return
         call check(error, grid%nang_per_atom(1) == 110, "largest per-atom angular size")
         if (allocated(error)) return
      end do
   end subroutine test_element_overrides

   !* ----------------------- Owner-level partition derivatives ----------------------- *!

   !> Fused switch values and first two derivatives against independent differences
   subroutine test_switch_jets(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      class(moist_math_grid_3d_partition_type), allocatable :: cell
      type(becke_partition_type) :: becke
      type(ssf_partition_type) :: ssf
      type(pvoronoi_partition_type) :: bump
      real(wp), parameter :: xs(9) = [-1.2_wp, -0.9_wp, -0.7_wp, -0.3_wp, 0.0_wp, &
         & 0.3_wp, 0.7_wp, 0.9_wp, 1.2_wp], step = 1.0e-5_wp
      integer, parameter :: offsets(4) = [-2, -1, 1, 2]
      real(wp) :: x, s, ds, d2s, discarded, discarded_ds, vals(4), derivs(4), first, second, negative
      integer :: scheme, stiffness, i, k

      ssf%a = 0.53_wp
      do scheme = partition_becke, partition_pvoronoi
         do stiffness = 1, 5, 2
            if (scheme /= partition_becke .and. stiffness /= 3) cycle
            if (allocated(cell)) deallocate (cell)
            select case (scheme)
            case (partition_becke)
               becke%k = stiffness
               allocate (cell, source=becke)
            case (partition_ssf)
               allocate (cell, source=ssf)
            case (partition_pvoronoi)
               allocate (cell, source=bump)
            case default
               call test_failed(error, "unknown partition in switch-jet test")
               return
            end select
            do i = 1, size(xs)
               x = xs(i)
               call production_switch_jet(cell, x, s, ds, d2s, error)
               if (allocated(error)) return
               call production_switch_jet(cell, -x, discarded, discarded_ds, negative, error)
               if (allocated(error)) return
               call check(error, abs(d2s + negative) < 2.0e-14_wp, "switch second derivative must be odd")
               if (allocated(error)) return
               do k = 1, 4
                  call production_switch_jet(cell, x + real(offsets(k), wp)*step, &
                     & vals(k), derivs(k), discarded, error)
                  if (allocated(error)) return
               end do
               first = (vals(1) - 8.0_wp*vals(2) + 8.0_wp*vals(3) - vals(4))/(12.0_wp*step)
               second = (derivs(1) - 8.0_wp*derivs(2) + 8.0_wp*derivs(3) - derivs(4))/(12.0_wp*step)
               call check(error, abs(first - ds) < 2.0e-9_wp .and. abs(second - d2s) < 2.0e-9_wp, &
                  & "generated first and second switch derivatives")
               if (allocated(error)) return
            end do
         end do
      end do
      call bump%eval_with_second_derivative(0.99_wp, s, ds, d2s)
      call check(error, s > 0.0_wp .and. ds < 0.0_wp .and. d2s > 0.0_wp .and. d2s < 1.0e-30_wp, &
         & "second derivatives must retain tiny positive switch tails")
   end subroutine test_switch_jets

   !> Production switch jets and agreement of the value, first- and second-order bindings
   !>
   !> @param[in]  cell   configured partition scheme
   !> @param[in]  x      pair coordinate
   !> @param[out] s      switch value
   !> @param[out] ds     first derivative
   !> @param[out] d2s    second derivative
   !> @param[out] error  test failure
   subroutine production_switch_jet(cell, x, s, ds, d2s, error)
      !> Configured partition scheme
      class(moist_math_grid_3d_partition_type), intent(in) :: cell
      !> Pair coordinate
      real(wp), intent(in) :: x
      !> Switch value
      real(wp), intent(out) :: s
      !> First derivative
      real(wp), intent(out) :: ds
      !> Second derivative
      real(wp), intent(out) :: d2s
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: value, s1, ds1

      call cell%eval_with_second_derivative(x, s, ds, d2s)
      call cell%eval_with_derivative(x, s1, ds1)
      value = cell%eval(x)
      call check(error, abs(s - value) < 2.0e-15_wp .and. abs(s - s1) < 2.0e-15_wp .and. &
         & abs(ds - ds1) < 2.0e-14_wp, "switch bindings must agree")
   end subroutine production_switch_jet

   !> Fixed-point partition partials and point partials against finite differences
   subroutine test_partition_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      class(moist_math_grid_3d_partition_type), allocatable :: partition
      real(wp) :: points(3, 4), moved_points(3, 4), xyz(3, 3), weights(4), seed(4)
      real(wp) :: gradient(3, 3), point_gradient(3, 4), energies(4), fd
      real(wp), parameter :: h = 1.0e-5_wp
      integer, parameter :: offsets(4) = [-2, -1, 1, 2], numbers(3) = [8, 1, 6]
      integer :: scheme, owner, channel, a, c, k, count, stiffness

      call make_molecule(mol)
      points(:, 1) = [0.1_wp, 0.4_wp, 0.3_wp]
      points(:, 2) = mol%xyz(:, 2)
      points(:, 3) = [1.6_wp, 0.9_wp, 0.2_wp]
      points(:, 4) = [-1.2_wp, -0.4_wp, -0.3_wp]
      seed = [0.3_wp, -0.2_wp, 0.7_wp, -0.4_wp]
      do scheme = partition_becke, partition_pvoronoi
         do stiffness = 1, 5, 2
            if (scheme /= partition_becke .and. stiffness /= 3) cycle
            points(:, 2) = mol%xyz(:, 2)
            ! One-iteration switch only once differentiable at a nucleus
            if (scheme == partition_becke .and. stiffness == 1) then
               points(:, 2) = points(:, 2) + [0.05_wp, -0.03_wp, 0.02_wp]
            end if
            if (allocated(partition)) deallocate (partition)
            select case (scheme)
            case (partition_becke)
               allocate (partition, source=becke_partition_type(k=stiffness))
            case (partition_ssf)
               allocate (partition, source=ssf_partition_type(a=0.53_wp))
            case default
               allocate (partition, source=pvoronoi_partition_type(width=2.3_wp, radii=[1.0_wp, 0.7_wp, 1.3_wp]))
            end select
            do owner = 1, 3
               gradient = 0.0_wp
               call partition%owner_gradient(owner, points, mol%xyz, numbers, seed, gradient, point_gradient, merr)
               call check_moist_error(error, merr, "partition reverse")
               if (allocated(error)) return
               do c = 1, 3
                  call check(error, sum(gradient(c, :)) + sum(point_gradient(c, :)), 0.0_wp, thr=1.0e-12_wp, &
                     & more="partition reverse is translation invariant")
                  if (allocated(error)) return
               end do
               do channel = 1, 2
                  count = merge(3, 4, channel == 1)
                  do a = 1, count
                     do c = 1, 3
                        do k = 1, 4
                           xyz = mol%xyz
                           moved_points = points
                           if (channel == 1) xyz(c, a) = xyz(c, a) + real(offsets(k), wp)*h
                           if (channel == 2) moved_points(c, a) = moved_points(c, a) + real(offsets(k), wp)*h
                           call partition%owner_weights(owner, moved_points, xyz, numbers, weights)
                           energies(k) = sum(seed*weights)
                        end do
                        fd = (energies(1) - 8.0_wp*energies(2) + 8.0_wp*energies(3) - energies(4))/(12.0_wp*h)
                        call check(error, ieee_is_finite(fd), "finite-difference partition gradient is not finite")
                        if (allocated(error)) return
                        if (channel == 1) then
                           call check(error, gradient(c, a), fd, thr=1.0e-10_wp)
                        else
                           call check(error, point_gradient(c, a), fd, thr=1.0e-10_wp)
                        end if
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_partition_fd

   !> Partition-unity derivatives retain contributions from tiny nonzero tails
   subroutine test_gradient_tails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: merr
      type(becke_partition_type) :: becke
      type(pvoronoi_partition_type) :: pvoronoi
      real(wp) :: xyz(3, 2), points(3, 1), gradient(3, 2, 2), point_gradient(3, 1, 2), scale
      integer :: owner, scheme, a, c

      xyz = 0.0_wp
      xyz(1, :) = [-1.0_wp, 1.0_wp]
      do scheme = 1, 2
         if (scheme == 1) then
            points(:, 1) = [2.0_wp, 0.001_wp, 0.0_wp]
         else
            points(:, 1) = [0.2475_wp, 0.1_wp, 0.0_wp]
         end if
         gradient = 0.0_wp
         do owner = 1, 2
            if (scheme == 1) then
               call becke%owner_gradient(owner, points, xyz, [1, 1], [1.0_wp], &
                  & gradient(:, :, owner), point_gradient(:, :, owner), merr)
            else
               call pvoronoi%owner_gradient(owner, points, xyz, [1, 1], [1.0_wp], &
                  & gradient(:, :, owner), point_gradient(:, :, owner), merr)
            end if
            call check_moist_error(error, merr, "tiny partition tail")
            if (allocated(error)) return
         end do
         scale = maxval(abs(gradient))
         call check(error, scale > 0.0_wp .and. scale < 1.0e-30_wp, "tail derivative must remain nonzero")
         if (allocated(error)) return
         do a = 1, size(gradient, 2)
            do c = 1, 3
               call check(error, sum(gradient(c, a, :)), 0.0_wp, thr=2.0e-12_wp*scale, &
                  & more="partition-unity nuclear derivatives must cancel at the tail scale")
               if (allocated(error)) return
            end do
         end do
         do a = 1, size(point_gradient, 2)
            do c = 1, 3
               call check(error, sum(point_gradient(c, a, :)), 0.0_wp, thr=2.0e-12_wp*scale, &
                  & more="partition-unity point derivatives must cancel at the tail scale")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_gradient_tails

   !> Partition-unity Hessian-vector products cancel at the tiny-tail scale
   subroutine test_hessian_vector_tails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: merr
      type(becke_partition_type) :: becke
      type(pvoronoi_partition_type) :: pvoronoi
      real(wp) :: xyz(3, 2), points(3, 1), direction(3, 2), hvp(3, 2, 2), point_hvp(3, 1, 2), scale
      integer :: owner, scheme, a, c

      xyz = 0.0_wp
      xyz(1, :) = [-1.0_wp, 1.0_wp]
      direction = reshape([0.1_wp, 0.3_wp, -0.2_wp, -0.4_wp, 0.2_wp, 0.5_wp], [3, 2])
      do scheme = 1, 2
         if (scheme == 1) then
            points(:, 1) = [2.0_wp, 0.001_wp, 0.0_wp]
         else
            points(:, 1) = [0.2475_wp, 0.1_wp, 0.0_wp]
         end if
         hvp = 0.0_wp
         do owner = 1, 2
            if (scheme == 1) then
               call becke%owner_hessian_vector(owner, points, xyz, [1, 1], [1.0_wp], &
                  & direction, [0.0_wp, 0.0_wp, 0.0_wp], hvp(:, :, owner), point_hvp(:, :, owner), merr)
            else
               call pvoronoi%owner_hessian_vector(owner, points, xyz, [1, 1], [1.0_wp], &
                  & direction, [0.0_wp, 0.0_wp, 0.0_wp], hvp(:, :, owner), point_hvp(:, :, owner), merr)
            end if
            call check_moist_error(error, merr, "partition Hessian tail")
            if (allocated(error)) return
         end do
         scale = maxval(abs(hvp))
         call check(error, scale > 0.0_wp .and. scale < 1.0e-25_wp, "Hessian tails must remain nonzero")
         if (allocated(error)) return
         do a = 1, size(hvp, 2)
            do c = 1, 3
               call check(error, sum(hvp(c, a, :)), 0.0_wp, thr=2.0e-10_wp*scale, &
                  & more="partition-unity nuclear curvature")
               if (allocated(error)) return
            end do
         end do
         do a = 1, size(point_hvp, 2)
            do c = 1, 3
               call check(error, sum(point_hvp(c, a, :)), 0.0_wp, thr=2.0e-10_wp*scale, &
                  & more="partition-unity point curvature")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_hessian_vector_tails

   !* ------------------------------------ Helpers ------------------------------------ *!

   !> Hydrogen molecule on the x axis, 1.4 bohr bond
   !>
   !> @param[out] mol  structure
   subroutine make_h2(mol)
      !> Structure
      type(structure_type), intent(out) :: mol

      call new(mol, [1, 1], reshape([-0.7_wp, 0.0_wp, 0.0_wp, 0.7_wp, 0.0_wp, 0.0_wp], [3, 2]))
   end subroutine make_h2

   !> Water, oxygen near the origin
   !>
   !> @param[out] mol  structure
   subroutine make_water(mol)
      !> Structure
      type(structure_type), intent(out) :: mol

      call new(mol, [8, 1, 1], reshape([0.05_wp, -0.02_wp, 0.12_wp, &
         & 0.0_wp, 1.43_wp, -0.98_wp, &
         & 0.0_wp, -1.43_wp, -0.98_wp], [3, 3]))
   end subroutine make_water

   !> Three distinct nuclei away from partition switching boundaries
   !>
   !> @param[out] mol  test molecule
   subroutine make_molecule(mol)
      !> Test molecule
      type(structure_type), intent(out) :: mol

      call new(mol, [8, 1, 6], reshape([-0.6_wp, -0.2_wp, 0.1_wp, &
         & 0.7_wp, 0.3_wp, -0.25_wp, 0.1_wp, 1.2_wp, 0.4_wp], [3, 3]))
   end subroutine make_molecule

end module test_math_grid_3d_partition
