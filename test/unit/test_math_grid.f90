!> Test suite for the moist math/grid submodule:
!>   - Lebedev angular grid normalisation
!>   - Chebyshev-2 radial quadrature
!>   - Becke partition-of-unity
!>   - Atom-centered molecular grid (auto + uniform modes)
module test_math_grid
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: center_at_origin
   use moist_math_grid_lebedev, only: get_angular_grid
   use moist_math_grid_radial, only: chebyshev2_radii
   use moist_math_grid_becke, only: becke_weights
   use moist_math_grid, only: molecular_grid_type, &
      & new_molecular_grid, new_molecular_grid_uniform
   implicit none (type, external)
   private

   public :: collect_math_grid

contains

   !> Collect all math_grid tests
   subroutine collect_math_grid(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("angular_weights_sum_to_one", test_angular_weights_sum), &
         new_unittest("angular_points_on_unit_sphere", test_angular_unit_sphere), &
         new_unittest("radial_chebyshev_gauss", test_radial_chebyshev_gauss), &
         new_unittest("becke_partition_of_unity", test_becke_partition_of_unity), &
         new_unittest("molecular_grid_single_h_gaussian", test_mol_grid_single_h), &
         new_unittest("molecular_grid_uniform_h_gaussian", test_mol_grid_uniform_h), &
         new_unittest("molecular_grid_integrate_constant", test_mol_grid_integrate_const), &
         new_unittest("molecular_grid_destroy_idempotent", test_mol_grid_destroy), &
         new_unittest("molecular_grid_pruning_threshold", test_mol_grid_pruning) &
         ]
   end subroutine collect_math_grid


   !> Lebedev weights should sum to 1 on the unit sphere.
   subroutine test_angular_weights_sum(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: npts_list(3) = [74, 230, 434]
      integer, parameter :: orders(3)    = [ 6,  12,  16]
      real(wp), allocatable :: xyz(:, :), w(:)
      type(mctc_error), allocatable :: merr
      integer :: k

      do k = 1, 3
         allocate(xyz(3, npts_list(k)), w(npts_list(k)))
         call get_angular_grid(orders(k), xyz, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, abs(sum(w) - 1.0_wp) < 1.0e-12_wp, &
            & "Lebedev weights do not sum to 1")
         deallocate(xyz, w)
         if (allocated(error)) return
      end do
   end subroutine test_angular_weights_sum


   !> Lebedev points should lie exactly on the unit sphere.
   subroutine test_angular_unit_sphere(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: npts_list(3) = [74, 230, 434]
      integer, parameter :: orders(3)    = [ 6,  12,  16]
      real(wp), allocatable :: xyz(:, :), w(:)
      type(mctc_error), allocatable :: merr
      integer :: k, i
      real(wp) :: r2

      do k = 1, 3
         allocate(xyz(3, npts_list(k)), w(npts_list(k)))
         call get_angular_grid(orders(k), xyz, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         do i = 1, npts_list(k)
            r2 = xyz(1, i)**2 + xyz(2, i)**2 + xyz(3, i)**2
            call check(error, abs(r2 - 1.0_wp) < 1.0e-12_wp, &
               & "Lebedev point off unit sphere")
            if (allocated(error)) exit
         end do
         deallocate(xyz, w)
         if (allocated(error)) return
      end do
   end subroutine test_angular_unit_sphere


   !> Chebyshev-2 radial quadrature (with r^2 dr Jacobian folded into w):
   !> integral_0^inf exp(-r^2) dr = sqrt(pi)/4 when r^2 included.
   !> Here we want integral_0^inf exp(-r^2) r^2 dr = sqrt(pi)/4; the
   !> chebyshev2_radii weights already include r^2 dr, so the sum
   !> reduces to sum_i w_i exp(-r_i^2).
   subroutine test_radial_chebyshev_gauss(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nr = 80
      real(wp) :: radii(nr), weights(nr)
      real(wp) :: result, expected
      integer :: ir

      call chebyshev2_radii(nr, 1.0_wp, radii, weights)
      result = 0.0_wp
      do ir = 1, nr
         result = result + weights(ir) * exp(-radii(ir)**2)
      end do
      expected = sqrt(pi) / 4.0_wp
      call check(error, abs(result - expected) < 1.0e-6_wp, &
         & "Chebyshev-2 quadrature of exp(-r^2) r^2 deviates from sqrt(pi)/4")
   end subroutine test_radial_chebyshev_gauss


   !> Becke partition weights must sum to 1 at every sample point.
   !> The molecule is MB16-43/H2 (any 2-atom system would do; the
   !> partition-of-unity property is geometry-independent).
   subroutine test_becke_partition_of_unity(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :)
      integer,  allocatable :: numbers(:)
      real(wp) :: samples(3, 5)
      real(wp) :: weights(2)
      integer  :: k

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      xyz = mol%xyz
      allocate(numbers(mol%nat))
      do k = 1, mol%nat
         numbers(k) = mol%num(mol%id(k))
      end do

      samples(:, 1) = [ 0.1_wp,  0.2_wp, -0.3_wp]
      samples(:, 2) = [ 1.0_wp,  0.0_wp,  0.0_wp]
      samples(:, 3) = [-0.5_wp,  0.4_wp,  0.7_wp]
      samples(:, 4) = [ 0.0_wp,  0.0_wp,  3.0_wp]
      samples(:, 5) = [ 2.0_wp, -1.0_wp,  0.5_wp]

      do k = 1, size(samples, 2)
         call becke_weights(samples(:, k), mol%nat, xyz, numbers, weights)
         call check(error, abs(sum(weights) - 1.0_wp) < 1.0e-12_wp, &
            & "Becke weights do not sum to 1")
         if (allocated(error)) return
         call check(error, weights(1) >= -1.0e-14_wp .and. weights(2) >= -1.0e-14_wp, &
            & "Becke weights unexpectedly negative")
         if (allocated(error)) return
      end do
   end subroutine test_becke_partition_of_unity


   !> Integrate exp(-|r|^2) over a centered MB16-43/H2 grid using default
   !> grid sizes. Analytic value is pi^(3/2). Reference is integrand-only
   !> and so does not depend on the carrier molecule, provided the grid
   !> samples r ~ O(1) densely enough - H2 (nat=2) satisfies this.
   subroutine test_mol_grid_single_h(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(molecular_grid_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: result, expected
      integer :: i

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid(grid, mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      result = 0.0_wp
      do i = 1, grid%npts
         result = result + grid%weights(i) * exp(-sum(grid%xyz(:, i)**2))
      end do
      expected = pi**1.5_wp

      call check(error, abs(result - expected) < 1.0e-3_wp, &
         & "Single-H Gaussian integral deviates from pi^(3/2)")
      call grid%destroy()
   end subroutine test_mol_grid_single_h


   !> Same integrand as above, but with the uniform constructor and
   !> Pople "fine" sizes (75, 302).
   subroutine test_mol_grid_uniform_h(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(molecular_grid_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: result, expected
      integer :: i

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid_uniform(grid, mol, nrad=75, nang=302, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      result = 0.0_wp
      do i = 1, grid%npts
         result = result + grid%weights(i) * exp(-sum(grid%xyz(:, i)**2))
      end do
      expected = pi**1.5_wp

      call check(error, abs(result - expected) < 1.0e-3_wp, &
         & "Uniform-grid H Gaussian integral deviates from pi^(3/2)")
      call grid%destroy()
   end subroutine test_mol_grid_uniform_h


   !> %integrate(f=1, result) should return sum(weights).
   subroutine test_mol_grid_integrate_const(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(molecular_grid_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: result, expected

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid(grid, mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call grid%integrate(one_function, result)
      expected = sum(grid%weights)

      call check(error, abs(result - expected) < 1.0e-10_wp, &
         & "integrate(f=1) != sum(weights)")
      call grid%destroy()
   end subroutine test_mol_grid_integrate_const


   !> destroy() must be safely callable twice.
   subroutine test_mol_grid_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(molecular_grid_type) :: grid
      type(mctc_error), allocatable :: merr

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid(grid, mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call grid%destroy()
      call grid%destroy()

      call check(error, grid%npts == 0, "npts should be 0 after destroy")
      if (allocated(error)) return
      call check(error, .not. allocated(grid%xyz), "xyz should be deallocated")
   end subroutine test_mol_grid_destroy


   !> All retained weights should have |w| >= default threshold (1e-14).
   subroutine test_mol_grid_pruning(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(molecular_grid_type) :: grid
      type(mctc_error), allocatable :: merr
      integer :: i

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call new_molecular_grid(grid, mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      do i = 1, grid%npts
         call check(error, abs(grid%weights(i)) >= 1.0e-14_wp, &
            & "Point with |weight| below threshold survived pruning")
         if (allocated(error)) exit
      end do
      call grid%destroy()
   end subroutine test_mol_grid_pruning


   !> Trivial `pure` integrand (f == 1) used by the integrate test.
   pure function one_function(r) result(val)
      real(wp), intent(in) :: r(3)
      real(wp) :: val
      val = 1.0_wp + 0.0_wp * r(1)   ! reference r to silence unused warning
   end function one_function

end module test_math_grid
