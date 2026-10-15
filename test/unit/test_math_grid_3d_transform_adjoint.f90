!> Euclidean transpose identities for the packed FFT and both NUFFT routes
module test_math_grid_3d_transform_adjoint
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type, moist_math_grid_3d_trafo_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, &
      & moist_math_grid_3d_molecular_trafo_type, new_molecular_point_grid, new_molecular_grid_trafo, molecular_grid_set_kgrid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   use test_helpers, only: check_moist_error, get_qc_handymod_recipe
   implicit none(type, external)
   private
   public :: collect_math_grid_3d_transform_adjoint

   !> A one-point transform implementing only the original energy contract
   type, extends(moist_math_grid_3d_trafo_type) :: energy_only_trafo_type
   contains
      procedure :: fft_r2k => energy_only_r2k
      procedure :: fft_k2r => energy_only_k2r
      procedure :: destroy => energy_only_destroy
   end type energy_only_trafo_type

contains

   subroutine collect_math_grid_3d_transform_adjoint(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("optional_derivatives", optional_derivatives), &
         & new_unittest("cartesian_even", cartesian_even), &
         & new_unittest("cartesian_odd", cartesian_odd), &
         & new_unittest("molecular_type12", molecular_type12), &
         & new_unittest("molecular_type3", molecular_type3), &
         & new_unittest("geometry_cartesian", geometry_cartesian), &
         & new_unittest("geometry_type12", geometry_type12), &
         & new_unittest("geometry_type3", geometry_type3)]
   end subroutine collect_math_grid_3d_transform_adjoint

   !> Energy-only backends remain concrete; unsupported derivatives return errors
   subroutine optional_derivatives(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(energy_only_trafo_type), target :: backend
      class(moist_math_grid_3d_trafo_type), pointer :: trafo
      type(volume_adjoint_type) :: acc
      real(wp) :: real_block(1, 1)
      complex(wp) :: complex_block(1, 1)
      integer :: operation
      character(len=32), parameter :: names(4) = [character(len=32) :: &
         & "fft_r2k_adj", "fft_k2r_adj", "accumulate_r2k_geometry", "accumulate_k2r_geometry"]

      trafo => backend
      real_block = 2.0_wp
      call trafo%prepare(1, err)
      if (.not. allocated(err)) call trafo%fft_r2k(real_block, complex_block, err)
      if (.not. allocated(err)) call trafo%fft_k2r(complex_block, real_block, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(real_block == 2.0_wp), more="the energy-only transform failed")
      if (allocated(error)) return
      call acc%init(1)
      acc%w_xyz = 3.0_wp
      acc%w_w = 4.0_wp
      acc%w_xi = 5.0_wp
      acc%w_kref = 6.0_wp
      do operation = 1, 4
         real_block = 2.0_wp
         complex_block = cmplx(2.0_wp, 0.0_wp, wp)
         select case (operation)
         case (1)
            call trafo%fft_r2k_adj(complex_block, real_block, err)
         case (2)
            call trafo%fft_k2r_adj(real_block, complex_block, err)
         case (3)
            call trafo%accumulate_r2k_geometry(real_block, complex_block, acc, err)
         case default
            call trafo%accumulate_k2r_geometry(complex_block, real_block, acc, err)
         end select
         call check(error, allocated(err), more="unsupported "//trim(names(operation))//" succeeded")
         if (allocated(error)) return
         call check(error, index(err%message, trim(names(operation))) > 0, more="error omitted the unsupported method")
         if (allocated(error)) return
         deallocate (err)
      end do
      call check(error, all(acc%w_xyz == 3.0_wp) .and. all(acc%w_w == 4.0_wp) .and. &
         & all(acc%w_xi == 5.0_wp) .and. all(acc%w_kref == 6.0_wp), more="an unsupported derivative changed the accumulator")
      call trafo%destroy()
   end subroutine optional_derivatives

   subroutine energy_only_r2k(self, f_r, f_k, error)
      class(energy_only_trafo_type), intent(inout) :: self
      real(wp), intent(inout), contiguous, target :: f_r(:, :)
      complex(wp), intent(out), contiguous, target :: f_k(:, :)
      type(moist_error), allocatable, intent(out) :: error
      f_k = cmplx(f_r, 0.0_wp, wp)
   end subroutine energy_only_r2k

   subroutine energy_only_k2r(self, f_k, f_r, error)
      class(energy_only_trafo_type), intent(inout) :: self
      complex(wp), intent(inout), contiguous, target :: f_k(:, :)
      real(wp), intent(out), contiguous, target :: f_r(:, :)
      type(moist_error), allocatable, intent(out) :: error
      f_r = real(f_k, wp)
   end subroutine energy_only_k2r

   subroutine energy_only_destroy(self)
      class(energy_only_trafo_type), intent(inout) :: self
   end subroutine energy_only_destroy

   subroutine cartesian_even(error)
      type(error_type), allocatable, intent(out) :: error
      call cartesian_case(8, error)
   end subroutine cartesian_even

   subroutine cartesian_odd(error)
      type(error_type), allocatable, intent(out) :: error
      call cartesian_case(9, error)
   end subroutine cartesian_odd

   subroutine cartesian_case(nx, error)
      integer, intent(in) :: nx
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      integer :: nt

      call new(mol, [1], reshape([0.3_wp, -0.4_wp, 0.7_wp], [3, 1]))
      call new_cartesian_point_grid(grid, err, nx, 7, 6, 0.73_wp, margin=0.0_wp)
      if (.not. allocated(err)) call grid%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do nt = 1, 4, 3
         grid%nthreads = nt
         call grid%new_trafo(trafo, err)
         if (.not. allocated(err)) call trafo%prepare(3, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check_pairing(grid, trafo, 1.0e-11_wp, error)
         if (allocated(error)) return
         deallocate (trafo)
      end do
   end subroutine cartesian_case

   subroutine molecular_type12(error)
      type(error_type), allocatable, intent(out) :: error
      call molecular_case(.true., error)
   end subroutine molecular_type12

   subroutine molecular_type3(error)
      type(error_type), allocatable, intent(out) :: error
      call molecular_case(.false., error)
   end subroutine molecular_type3

   subroutine molecular_case(type12, error)
      logical, intent(in) :: type12
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: grid
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      type(moist_math_grid_atomic_recipe_type) :: recipe

      call new(mol, [1, 8], reshape([0.3_wp, -0.4_wp, 0.7_wp, 1.7_wp, 0.2_wp, -0.1_wp], [3, 2]))
      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(grid, err, recipe=recipe, dr=1.1_wp)
      if (.not. allocated(err)) call grid%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      ! The transpose must work even when a quadrature weight cannot be inverted.
      grid%w(1) = 0.0_wp
      grid%w(2) = -abs(grid%w(2))
      grid%nthreads = 2
      call new_molecular_grid_trafo(trafo, grid, use_type12=type12, nufft_tol=1.0e-12_wp)
      call trafo%prepare(3, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, trafo%uses_type12() .eqv. type12, more="requested NUFFT route not exercised")
      if (allocated(error)) return
      call check_pairing(grid, trafo, 2.0e-9_wp, error)
   end subroutine molecular_case

   !> Arbitrary complex seeds include imaginary boundary/Nyquist entries;
   !> they deliberately need not satisfy Hermitian symmetry.
   subroutine check_pairing(grid, trafo, tolerance, error)
      class(moist_math_grid_3d_type), intent(in) :: grid
      class(moist_math_grid_3d_trafo_type), intent(inout) :: trafo
      real(wp), intent(in) :: tolerance
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp), allocatable :: x(:, :), work_r(:, :), adj_r(:, :), back(:, :)
      complex(wp), allocatable :: y(:, :), work_k(:, :), adj_k(:, :), forward(:, :)
      real(wp) :: lhs, rhs
      integer :: i, col

      allocate (x(grid%ngrid, 3), work_r(grid%ngrid, 3), adj_r(grid%ngrid, 3), back(grid%ngrid, 3))
      allocate (y(grid%npts_k, 3), work_k(grid%npts_k, 3), adj_k(grid%npts_k, 3), forward(grid%npts_k, 3))
      do col = 1, 3
         do i = 1, grid%ngrid
            x(i, col) = sin(0.13_wp*i + col) + 0.2_wp*cos(0.07_wp*i*col)
         end do
         do i = 1, grid%npts_k
            y(i, col) = cmplx(cos(0.19_wp*i + col), sin(0.11_wp*i*col), wp)
         end do
      end do
      work_r = x
      call trafo%fft_r2k(work_r, forward, err)
      if (.not. allocated(err)) then
         work_k = y
         call trafo%fft_r2k_adj(work_k, adj_r, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) return
      lhs = real(sum(conjg(y)*forward), wp)
      rhs = sum(x*adj_r)
      call check(error, lhs, rhs, thr=tolerance*max(1.0_wp, abs(lhs), abs(rhs)), more="forward transform transpose")
      if (allocated(error)) return
      work_k = y
      call trafo%fft_k2r(work_k, back, err)
      if (.not. allocated(err)) then
         work_r = x
         call trafo%fft_k2r_adj(work_r, adj_k, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) return
      lhs = sum(x*back)
      rhs = real(sum(conjg(adj_k)*y), wp)
      call check(error, lhs, rhs, thr=tolerance*max(1.0_wp, abs(lhs), abs(rhs)), more="backward transform transpose")
   end subroutine check_pairing

   !> Both molecular plans differentiate points, measures and reference motion.
   subroutine geometry_type12(error)
      type(error_type), allocatable, intent(out) :: error
      call molecular_geometry_case(.true., error)
   end subroutine geometry_type12

   subroutine geometry_type3(error)
      type(error_type), allocatable, intent(out) :: error
      call molecular_geometry_case(.false., error)
   end subroutine geometry_type3

   subroutine molecular_geometry_case(type12, error)
      logical, intent(in) :: type12
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: grid, moved
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(volume_adjoint_type) :: acc
      real(wp), allocatable :: f(:, :), seed_r(:, :), out_r(:, :), dx(:, :), dw(:)
      complex(wp), allocatable :: fk(:, :), seed_k(:, :), work_k(:, :), out_k(:, :)
      real(wp), parameter :: eps = 2.0e-5_wp, reference(3) = [0.1_wp, 0.3_wp, -0.2_wp]
      real(wp), parameter :: dref(3) = [0.13_wp, -0.21_wp, 0.17_wp]
      real(wp) :: ep, em, dot, sign
      integer :: i, col, route, channel, side

      call new(mol, [1, 8], reshape([0.3_wp, -0.4_wp, 0.7_wp, 1.7_wp, 0.2_wp, -0.1_wp], [3, 2]))
      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(grid, err, recipe=recipe, dr=1.1_wp)
      if (.not. allocated(err)) call grid%update(mol, err)
      if (.not. allocated(err)) call molecular_grid_set_kgrid(grid, 1.1_wp, err, buffer=2.0_wp, reference=reference)
      call check_moist_error(error, err)
      if (allocated(error)) return
      grid%w(1) = 0.0_wp
      grid%w(2) = -abs(grid%w(2))
      grid%nthreads = 2
      allocate (f(grid%ngrid, 3), seed_r(grid%ngrid, 3), out_r(grid%ngrid, 3), dx(3, grid%ngrid), dw(grid%ngrid))
      allocate (fk(grid%npts_k, 3), seed_k(grid%npts_k, 3), work_k(grid%npts_k, 3), out_k(grid%npts_k, 3))
      do i = 1, grid%ngrid
         dx(:, i) = [sin(0.13_wp*i), cos(0.17_wp*i), sin(0.23_wp*i)]
         dw(i) = cos(0.11_wp*i)
         do col = 1, 3
            f(i, col) = sin(0.19_wp*i + col)
            seed_r(i, col) = cos(0.31_wp*i + col)
         end do
      end do
      do i = 1, grid%npts_k
         do col = 1, 3
            fk(i, col) = cmplx(sin(0.17_wp*i + col), cos(0.31_wp*i*col), wp)
            seed_k(i, col) = cmplx(cos(0.29_wp*i*col), sin(0.13_wp*i + col), wp)
         end do
      end do
      do route = 1, 2
         call new_molecular_grid_trafo(trafo, grid, nufft_tol=1.0e-12_wp, use_type12=type12)
         call trafo%prepare(3, err)
         call acc%init(grid%ngrid)
         if (.not. allocated(err)) then
            if (route == 1) then
               call trafo%accumulate_r2k_geometry(f, seed_k, acc, err)
            else
               call trafo%accumulate_k2r_geometry(fk, seed_r, acc, err)
            end if
         end if
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, maxval(abs(sum(acc%w_xyz, dim=2) + acc%w_kref)), 0.0_wp, &
            & thr=2.0e-9_wp, more="phase translation identity")
         if (allocated(error)) return
         if (route == 2) then
            call check(error, maxval(abs(acc%w_w)), 0.0_wp, thr=0.0_wp, more="inverse has no measure dependency")
            if (allocated(error)) return
         end if
         do channel = 1, 3
            select case (channel)
            case (1)
               dot = sum(acc%w_xyz*dx)
            case (2)
               dot = sum(acc%w_w*dw)
            case default
               dot = sum(acc%w_kref*dref)
            end select
            do side = 1, 2
               sign = real(3 - 2*side, wp)
               moved = grid
               select case (channel)
               case (1)
                  moved%xyz = grid%xyz + sign*eps*dx
               case (2)
                  moved%w = grid%w + sign*eps*dw
               case default
                  call molecular_grid_set_kgrid(moved, 1.1_wp, err, buffer=2.0_wp, reference=reference + sign*eps*dref)
               end select
               call new_molecular_grid_trafo(trafo, moved, nufft_tol=1.0e-12_wp, use_type12=type12)
               if (.not. allocated(err)) call trafo%prepare(3, err)
               call check_moist_error(error, err)
               if (allocated(error)) return
               call check(error, moved%npts_k == grid%npts_k, more="FD changed reciprocal topology")
               if (allocated(error)) return
               if (route == 1) then
                  out_r = f
                  call trafo%fft_r2k(out_r, out_k, err)
                  ep = real(sum(conjg(seed_k)*out_k), wp)
               else
                  work_k = fk
                  call trafo%fft_k2r(work_k, out_r, err)
                  ep = sum(seed_r*out_r)
               end if
               call check_moist_error(error, err)
               if (allocated(error)) return
               if (side == 1) em = ep
            end do
            call check(error, dot, (em - ep)/(2.0_wp*eps), &
               & thr=2.0e-7_wp*max(1.0_wp, abs(dot)), more="Fourier geometry finite difference")
            if (allocated(error)) return
         end do
      end do
   end subroutine molecular_geometry_case

   !> Cartesian centroid motion leaves the discrete FFT unchanged.
   subroutine geometry_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      type(volume_adjoint_type) :: acc
      real(wp), allocatable :: f(:, :), out_r(:, :), old_r(:, :)
      complex(wp), allocatable :: fk(:, :), old_k(:, :), work_k(:, :)
      call new(mol, [1], reshape([0.3_wp, -0.4_wp, 0.7_wp], [3, 1]))
      call new_cartesian_point_grid(grid, err, 9, 7, 6, 0.73_wp, margin=0.0_wp)
      if (.not. allocated(err)) call grid%update(mol, err)
      if (.not. allocated(err)) call grid%new_trafo(trafo, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (f(grid%ngrid, 2), out_r(grid%ngrid, 2), old_r(grid%ngrid, 2))
      allocate (fk(grid%npts_k, 2), old_k(grid%npts_k, 2), work_k(grid%npts_k, 2))
      f = 0.31_wp
      fk = cmplx(0.17_wp, -0.23_wp, wp)
      call acc%init(grid%ngrid)
      call trafo%accumulate_r2k_geometry(f, fk, acc, err)
      if (.not. allocated(err)) call trafo%accumulate_k2r_geometry(fk, f, acc, err)
      if (.not. allocated(err)) call trafo%fft_r2k(f, old_k, err)
      work_k = fk
      if (.not. allocated(err)) call trafo%fft_k2r(work_k, old_r, err)
      mol%xyz = mol%xyz + 0.13_wp
      deallocate (trafo)
      if (.not. allocated(err)) call grid%update(mol, err)
      if (.not. allocated(err)) call grid%new_trafo(trafo, err)
      if (.not. allocated(err)) call trafo%fft_r2k(f, work_k, err)
      if (.not. allocated(err)) call trafo%fft_k2r(fk, out_r, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(work_k - old_k)) + maxval(abs(out_r - old_r)) + &
         & maxval(abs(acc%w_xyz)) + maxval(abs(acc%w_w)) + maxval(abs(acc%w_kref)), 0.0_wp, thr=0.0_wp)
   end subroutine geometry_cartesian

end module test_math_grid_3d_transform_adjoint
