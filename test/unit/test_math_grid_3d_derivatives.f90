!> Finite-difference validation of the volume-grid nuclear gradient and Hessian
!>
!> - One test per molecule; each runs every 3D grid: Cartesian point, Cartesian
!>   Gaussian, molecular point and molecular Gaussian
!> - Gradient: every adjoint channel of a linear grid observable, see
!>   `check_grid_gradient_fd`
!> - Hessian: fixed position and weight adjoints, see `check_grid_hessian_fd`
!> - Partition schemes besides the default Becke on a small molecule, and the
!>   Gaussian width term for tiny retained weights
module test_math_grid_3d_derivatives
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: get_qc_handymod_recipe, check_moist_error
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid, &
      & new_cartesian_gaussian_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid, &
      & new_molecular_gaussian_grid, partition_becke, partition_ssf, partition_pvoronoi
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   implicit none(type, external)
   private

   public :: collect_math_grid_3d_derivatives

   !> Gradient finite-difference step (bohr)
   real(wp), parameter :: GRADIENT_STEP = 4.0E-4_wp
   !> Gradient thresholds against the finite differences of the grid observable
   real(wp), parameter :: GRADIENT_ABS_THR = 5.0E-9_wp
   real(wp), parameter :: GRADIENT_REL_THR = 3.0E-8_wp
   !> Hessian finite-difference step (bohr), applied to the reverse gradient
   real(wp), parameter :: HESSIAN_STEP = 2.0E-4_wp
   !> Hessian thresholds against the finite differences of the reverse gradient
   real(wp), parameter :: HESSIAN_ABS_THR = 1.0E-8_wp
   real(wp), parameter :: HESSIAN_REL_THR = 1.0E-7_wp
   !> Exact identities of the gradient: translation and additive accumulation
   real(wp), parameter :: IDENTITY_ABS_THR = 1.0E-11_wp
   real(wp), parameter :: IDENTITY_REL_THR = 1.0E-12_wp
   !> Exact identities of the Hessian: symmetry, directional product
   real(wp), parameter :: CURVATURE_ABS_THR = 1.0E-10_wp
   real(wp), parameter :: CURVATURE_REL_THR = 1.0E-9_wp

   !> Six-point O(h**6) central stencil, `f' = sum_k STENCIL_WEIGHTS(k)*(f(k*h) - f(-k*h))/h`;
   !> differencing before weighting cancels identical values exactly
   real(wp), parameter :: STENCIL_WEIGHTS(3) = [45.0_wp, -9.0_wp, 1.0_wp]/60.0_wp

   !> Cartesian box: points per axis and spacing (bohr)
   integer, parameter :: CARTESIAN_POINTS = 64
   real(wp), parameter :: CARTESIAN_DR = 0.6_wp
   !> Nucleus-to-face distance of the Cartesian box (the production default)
   real(wp), parameter :: CARTESIAN_MARGIN = 10.0_wp
   !> Molecular test recipe: radial nodes, Lebedev degree and outer radius (bohr)
   integer, parameter :: MOLECULAR_NRAD = 5
   integer, parameter :: MOLECULAR_DEGREE = 5
   real(wp), parameter :: MOLECULAR_RMAX = 5.0_wp

contains

   !> Collect the finite-difference derivative tests
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_math_grid_3d_derivatives(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("single_atom", test_single_atom), &
                  new_unittest("dimer", test_dimer), &
                  new_unittest("amino20x4_gly_xab", test_amino20x4_gly_xab), &
                  new_unittest("mb16_43_01", test_mb16_43_01), &
                  new_unittest("but14diol_1", test_but14diol_1), &
                  new_unittest("il16_008", test_il16_008), &
                  new_unittest("partition_schemes", test_partition_schemes), &
                  new_unittest("gaussian_tiny_weights", test_gaussian_tiny_weights) &
                  ]
   end subroutine collect_math_grid_3d_derivatives

   !> Single oxygen atom
   subroutine test_single_atom(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [8], reshape([0.3_wp, -0.4_wp, 0.2_wp], [3, 1]))
      call do_test(error, mol)
   end subroutine test_single_atom

   !> Off-axis heteronuclear dimer
   subroutine test_dimer(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [6, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, 2.1_wp, 0.4_wp, -0.3_wp], [3, 2]))
      call do_test(error, mol)
   end subroutine test_dimer

   !> Amino20x4 GLY_xab
   subroutine test_amino20x4_gly_xab(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "Amino20x4", "GLY_xab")
      call do_test(error, mol)
   end subroutine test_amino20x4_gly_xab

   !> MB16-43 01
   subroutine test_mb16_43_01(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "01")
      call do_test(error, mol)
   end subroutine test_mb16_43_01

   !> But14diol 1
   subroutine test_but14diol_1(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "But14diol", "1")
      call do_test(error, mol)
   end subroutine test_but14diol_1

   !> IL16 008
   subroutine test_il16_008(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "IL16", "008")
      call do_test(error, mol)
   end subroutine test_il16_008

   !> Partition schemes besides the default Becke k=3, on both molecular grids
   !>
   !> Becke k=1 and k=5, SSF and power Voronoi on an oxygen, hydrogen, carbon triple
   subroutine test_partition_schemes(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: schemes(4) = [partition_becke, partition_becke, partition_ssf, partition_pvoronoi]
      integer, parameter :: stiffness(2) = [1, 5]
      class(moist_math_grid_3d_type), allocatable :: grid
      type(structure_type) :: mol
      ! Unallocated for SSF and power Voronoi, which take no stiffness
      integer, allocatable :: becke_k
      integer :: is, kind

      call new(mol, [8, 1, 6], reshape([-0.6_wp, -0.2_wp, 0.1_wp, 0.7_wp, 0.3_wp, -0.25_wp, &
         & 0.1_wp, 1.2_wp, 0.4_wp], [3, 3]))
      do is = 1, size(schemes)
         if (allocated(becke_k)) deallocate (becke_k)
         if (is <= size(stiffness)) becke_k = stiffness(is)
         do kind = 3, 4
            call make_grid(kind, mol, grid, error, scheme=schemes(is), becke_k=becke_k)
            if (allocated(error)) return
            call check_grid_gradient_fd(grid, mol, error)
            if (allocated(error)) return
            call check_grid_hessian_fd(grid, mol, error)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_partition_schemes

   !> Gaussian width term stays finite for retained weights far below the xi0/w overflow
   !>
   !> Regression: `xi0/w` overflows for weights below ~1e-240, so a zero width adjoint
   !> must not touch it and a width adjoint on ordinary points must not reach it
   subroutine test_gaussian_tiny_weights(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(volume_adjoint_type) :: acc
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      real(wp) :: xyz(3, 6), gradient(3, 6)
      integer :: a

      xyz = 0.0_wp
      do a = 1, 6
         xyz(1, a) = 1.4_wp*real(a - 1, wp)
      end do
      call new(mol, [1, 1, 1, 1, 1, 1], xyz)
      call get_qc_handymod_recipe(recipe, merr, nrad=40, degree=11, rmax=10.0_wp)
      if (.not. allocated(merr)) call new_molecular_gaussian_grid(grid, merr, recipe=recipe, &
         & partition=partition_becke, reciprocal=.false., weight_threshold=0.0_wp, pruning_threshold=0.0_wp)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      call check_moist_error(error, merr, "tiny-weight grid")
      if (allocated(error)) return
      call check(error, minval(grid%w) < 1.0e-240_wp, "fixture must retain weights that overflow xi0/w")
      if (allocated(error)) return

      call acc%init(grid%ngrid)
      acc%w_w = 1.0_wp
      gradient = 0.0_wp
      call grid%get_volume_gradient(acc, gradient, merr)
      call check_moist_error(error, merr, "zero width adjoint")
      if (allocated(error)) return
      call check(error, all(ieee_is_finite(gradient)), "zero width adjoint gave a non-finite gradient")
      if (allocated(error)) return

      where (grid%w > 1.0e-8_wp) acc%w_xi = 0.1_wp
      gradient = 0.0_wp
      call grid%get_volume_gradient(acc, gradient, merr)
      call check_moist_error(error, merr, "width adjoint on ordinary points")
      if (allocated(error)) return
      call check(error, all(ieee_is_finite(gradient)), "width adjoint gave a non-finite gradient")
   end subroutine test_gaussian_tiny_weights

   !> Gradient and Hessian finite differences of `mol` on every 3D grid
   !> @param[in]  mol   Reference geometry
   subroutine do_test(error, mol)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Reference geometry
      type(structure_type), intent(in) :: mol

      class(moist_math_grid_3d_type), allocatable :: grid
      integer :: kind

      do kind = 1, 4
         call make_grid(kind, mol, grid, error)
         if (allocated(error)) return
         call check_grid_gradient_fd(grid, mol, error)
         if (allocated(error)) return
         call check_grid_hessian_fd(grid, mol, error)
         if (allocated(error)) return
      end do
   end subroutine do_test

   !> Construct and update one of the four 3D grids for `mol`
   !>
   !> - 1, 2: Cartesian point and Gaussian grids; the box holds the molecule
   !>   with the default margin
   !> - 3, 4: molecular point and Gaussian grids, Becke partition unless `scheme`
   !>   is given, small recipe, thresholds that keep point membership under the
   !>   stencil steps
   !>
   !> @param[in]  kind     Grid kind, 1..4
   !> @param[in]  mol      Geometry to update with
   !> @param[out] grid     Updated grid
   !> @param[out] error    Test failure
   !> @param[in]  scheme   Molecular partition scheme
   !> @param[in]  becke_k  Becke stiffness
   subroutine make_grid(kind, mol, grid, error, scheme, becke_k)
      !> Grid kind
      integer, intent(in) :: kind
      !> Geometry
      type(structure_type), intent(in) :: mol
      !> Updated grid
      class(moist_math_grid_3d_type), allocatable, intent(out) :: grid
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Molecular partition scheme
      integer, intent(in), optional :: scheme
      !> Becke stiffness
      integer, intent(in), optional :: becke_k

      type(moist_math_grid_3d_cartesian_type), allocatable :: cart
      type(moist_math_grid_3d_molecular_type), allocatable :: molecular
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(mctc_error), allocatable :: merr

      select case (kind)
      case (1, 2)
         allocate (cart)
         if (kind == 1) then
            call new_cartesian_point_grid(cart, merr, nx=CARTESIAN_POINTS, ny=CARTESIAN_POINTS, &
               & nz=CARTESIAN_POINTS, dr=CARTESIAN_DR, margin=CARTESIAN_MARGIN)
         else
            call new_cartesian_gaussian_grid(cart, merr, nx=CARTESIAN_POINTS, ny=CARTESIAN_POINTS, &
               & nz=CARTESIAN_POINTS, dr=CARTESIAN_DR, margin=CARTESIAN_MARGIN)
         end if
         if (.not. allocated(merr)) call cart%update(mol, merr)
         call check_moist_error(error, merr, "Cartesian grid setup")
         if (allocated(error)) return
         call move_alloc(cart, grid)
      case default
         call get_qc_handymod_recipe(recipe, merr, nrad=MOLECULAR_NRAD, degree=MOLECULAR_DEGREE, &
            & rmax=MOLECULAR_RMAX)
         allocate (molecular)
         if (.not. allocated(merr)) then
            if (kind == 3) then
               call new_molecular_point_grid(molecular, merr, recipe=recipe, reciprocal=.false., &
                  & weight_threshold=1.0e-5_wp, pruning_threshold=0.01_wp, partition=scheme, becke_k=becke_k)
            else
               call new_molecular_gaussian_grid(molecular, merr, recipe=recipe, reciprocal=.false., &
                  & weight_threshold=1.0e-5_wp, pruning_threshold=0.01_wp, xi0_factor=1.2_wp, &
                  & partition=scheme, becke_k=becke_k)
            end if
         end if
         molecular%nthreads = 1
         if (.not. allocated(merr)) call molecular%update(mol, merr)
         call check_moist_error(error, merr, "molecular grid setup")
         if (allocated(error)) return
         call move_alloc(molecular, grid)
      end select
   end subroutine make_grid

   !> Reverse grid gradient of every adjoint channel against finite differences
   !>
   !> - Linear synthetic energy `E = sum(w_xyz*xyz) + sum(w_w*w) + sum(w_xi*xi0)`
   !>   with fixed seeds; channels: positions, weights, widths (Gaussian grids), all
   !> - Six-point O(h**6) stencil with `GRADIENT_STEP`; every step must keep point
   !>   membership and move each point rigidly with its owner (molecular) or the
   !>   centroid (Cartesian); each displaced grid serves every channel
   !> - Also checks the translation identity and additive accumulation
   !> - Every element compared one by one with absolute and relative thresholds
   !>
   !> @param[in]  grid   Successfully updated grid
   !> @param[in]  mol    Geometry of that update
   !> @param[out] error  Test failure
   subroutine check_grid_gradient_fd(grid, mol, error)
      !> Successfully updated grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Geometry of that update
      type(structure_type), intent(in) :: mol
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nchannel = 4
      class(moist_math_grid_3d_type), allocatable :: moved
      type(structure_type) :: displaced
      type(volume_adjoint_type) :: acc(nchannel)
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: gradient(:, :, :), fd(:, :, :), again(:, :)
      real(wp) :: energies(-3:3, nchannel), expected(3), total(3)
      logical :: active(nchannel)
      integer :: channel, i, a, c, k

      allocate (gradient(3, mol%nat, nchannel), fd(3, mol%nat, nchannel), again(3, mol%nat))
      active = [.true., .true., allocated(grid%xi0), .true.]
      do channel = 1, nchannel
         if (.not. active(channel)) cycle
         call acc(channel)%init(grid%ngrid)
         do i = 1, grid%ngrid
            if (channel == 1 .or. channel == 4) then
               acc(channel)%w_xyz(:, i) = [sin(real(i, wp)), cos(0.3_wp*real(i, wp)), &
                  & sin(0.7_wp*real(i, wp))]/real(grid%ngrid, wp)
            end if
            if (channel == 2 .or. channel == 4) acc(channel)%w_w(i) = cos(0.2_wp*real(i, wp))
            if (allocated(grid%xi0) .and. (channel == 3 .or. channel == 4)) then
               acc(channel)%w_xi(i) = 0.04_wp*sin(0.37_wp*real(i, wp))*grid%w(i)
            end if
         end do
         gradient(:, :, channel) = 0.0_wp
         call grid%get_volume_gradient(acc(channel), gradient(:, :, channel), merr)
         call check_moist_error(error, merr, "reverse contraction")
         if (allocated(error)) return
         again = 0.25_wp
         call grid%get_volume_gradient(acc(channel), again, merr)
         call check_moist_error(error, merr, "additive contraction")
         if (allocated(error)) return
         expected = sum(acc(channel)%w_xyz, dim=2)
         total = sum(gradient(:, :, channel), dim=2)
         do c = 1, 3
            call check(error, total(c), expected(c), thr_abs=IDENTITY_ABS_THR, thr_rel=IDENTITY_REL_THR, &
               & more="translation identity")
            if (allocated(error)) return
         end do
         do a = 1, mol%nat
            do c = 1, 3
               call check(error, again(c, a), gradient(c, a, channel) + 0.25_wp, thr_abs=IDENTITY_ABS_THR, &
                  & thr_rel=IDENTITY_REL_THR, more="additive accumulation")
               if (allocated(error)) return
            end do
         end do
      end do

      allocate (moved, source=grid)
      energies = 0.0_wp
      do a = 1, mol%nat
         do c = 1, 3
            do k = -3, 3
               if (k == 0) cycle
               displaced = mol
               displaced%xyz(c, a) = displaced%xyz(c, a) + real(k, wp)*GRADIENT_STEP
               call moved%update(displaced, merr)
               call check_moist_error(error, merr, "displaced grid update")
               if (allocated(error)) return
               call check_rigid_motion(grid, moved, a, c, real(k, wp)*GRADIENT_STEP, error)
               if (allocated(error)) return
               do channel = 1, nchannel
                  if (active(channel)) energies(k, channel) = grid_observable_energy(moved, acc(channel))
               end do
            end do
            do channel = 1, nchannel
               fd(c, a, channel) = sum(STENCIL_WEIGHTS*(energies(1:3, channel) - energies(-1:-3:-1, channel))) &
                  & /GRADIENT_STEP
            end do
         end do
      end do
      do channel = 1, nchannel
         if (.not. active(channel)) cycle
         do a = 1, mol%nat
            do c = 1, 3
               ! test-drive accepts a non-finite expected value
               call check(error, ieee_is_finite(fd(c, a, channel)), "finite-difference gradient is not finite")
               if (allocated(error)) return
               call check(error, gradient(c, a, channel), fd(c, a, channel), thr_abs=GRADIENT_ABS_THR, &
                  & thr_rel=GRADIENT_REL_THR, more="gradient against finite differences")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine check_grid_gradient_fd

   !> Nuclear volume Hessian against finite differences of the reverse gradient
   !>
   !> - Fixed position and weight adjoints; width adjoints stay zero because
   !>   geometry-dependent width curvature is unsupported
   !> - Full Hessian symmetric, Hessian-vector product equal to the full Hessian
   !>   times the direction, rigid translation in the nullspace
   !> - Six-point O(h**6) stencil of `get_volume_gradient` with `HESSIAN_STEP`
   !> - Every element compared one by one with absolute and relative thresholds
   !>
   !> @param[in]  grid   Successfully updated grid
   !> @param[in]  mol    Geometry of that update
   !> @param[out] error  Test failure
   subroutine check_grid_hessian_fd(grid, mol, error)
      !> Successfully updated grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Geometry of that update
      type(structure_type), intent(in) :: mol
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      class(moist_math_grid_3d_type), allocatable :: moved
      type(structure_type) :: displaced
      type(volume_adjoint_type) :: acc
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: hessian(:, :, :, :), fd(:, :, :, :), grads(:, :, :)
      real(wp), allocatable :: direction(:, :), hvp(:, :), expected(:, :), translation(:, :)
      integer :: n, i, a, b, c, d, k

      n = 3*mol%nat
      allocate (hessian(3, mol%nat, 3, mol%nat), fd(3, mol%nat, 3, mol%nat), grads(3, mol%nat, -3:3))
      allocate (direction(3, mol%nat), hvp(3, mol%nat), translation(3, mol%nat))
      call acc%init(grid%ngrid)
      acc%w_xyz = 0.3_wp
      do i = 1, grid%ngrid
         acc%w_w(i) = cos(0.2_wp*real(i, wp))
      end do
      hessian = 0.0_wp
      call grid%get_volume_hessian(acc, hessian, merr)
      call check_moist_error(error, merr, "full nuclear Hessian")
      if (allocated(error)) return
      direction = reshape(sin(real([(i, i=1, n)], wp)), [3, mol%nat])
      hvp = 0.0_wp
      call grid%get_volume_hessian_vector(acc, direction, hvp, merr)
      call check_moist_error(error, merr, "Hessian-vector contraction")
      if (allocated(error)) return
      expected = reshape(matmul(reshape(hessian, [n, n]), reshape(direction, [n])), [3, mol%nat])
      direction = spread([0.3_wp, -0.2_wp, 0.7_wp], 2, mol%nat)
      translation = 0.0_wp
      call grid%get_volume_hessian_vector(acc, direction, translation, merr)
      call check_moist_error(error, merr, "translation Hessian")
      if (allocated(error)) return
      do b = 1, mol%nat
         do d = 1, 3
            call check(error, hvp(d, b), expected(d, b), thr_abs=CURVATURE_ABS_THR, &
               & thr_rel=CURVATURE_REL_THR, more="directional product")
            if (allocated(error)) return
            call check(error, translation(d, b), 0.0_wp, thr_abs=CURVATURE_ABS_THR, &
               & thr_rel=CURVATURE_REL_THR, more="translation nullspace")
            if (allocated(error)) return
            do a = 1, mol%nat
               do c = 1, 3
                  call check(error, hessian(c, a, d, b), hessian(d, b, c, a), thr_abs=CURVATURE_ABS_THR, &
                     & thr_rel=CURVATURE_REL_THR, more="Hessian symmetry")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do

      allocate (moved, source=grid)
      do a = 1, mol%nat
         do c = 1, 3
            do k = -3, 3
               if (k == 0) cycle
               displaced = mol
               displaced%xyz(c, a) = displaced%xyz(c, a) + real(k, wp)*HESSIAN_STEP
               call moved%update(displaced, merr)
               call check_moist_error(error, merr, "displaced Hessian geometry")
               if (allocated(error)) return
               call check_rigid_motion(grid, moved, a, c, real(k, wp)*HESSIAN_STEP, error)
               if (allocated(error)) return
               grads(:, :, k) = 0.0_wp
               call moved%get_volume_gradient(acc, grads(:, :, k), merr)
               call check_moist_error(error, merr, "displaced reverse gradient")
               if (allocated(error)) return
            end do
            fd(:, :, c, a) = 0.0_wp
            do k = 1, 3
               fd(:, :, c, a) = fd(:, :, c, a) + STENCIL_WEIGHTS(k)*(grads(:, :, k) - grads(:, :, -k))
            end do
            fd(:, :, c, a) = fd(:, :, c, a)/HESSIAN_STEP
         end do
      end do
      do b = 1, mol%nat
         do d = 1, 3
            do a = 1, mol%nat
               do c = 1, 3
                  call check(error, ieee_is_finite(fd(c, a, d, b)), "finite-difference Hessian is not finite")
                  if (allocated(error)) return
                  call check(error, hessian(c, a, d, b), fd(c, a, d, b), thr_abs=HESSIAN_ABS_THR, &
                     & thr_rel=HESSIAN_REL_THR, more="Hessian against finite differences")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine check_grid_hessian_fd

   !> Displaced update keeps membership and moves points rigidly
   !>
   !> Molecular points follow their owner; Cartesian points follow the centroid
   !>
   !> @param[in]  grid   Reference grid
   !> @param[in]  moved  Grid updated with atom `a` displaced along `c`
   !> @param[in]  a      Displaced atom
   !> @param[in]  c      Displaced Cartesian component
   !> @param[in]  shift  Displacement (bohr)
   !> @param[out] error  Test failure
   subroutine check_rigid_motion(grid, moved, a, c, shift, error)
      !> Reference grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Displaced grid
      class(moist_math_grid_3d_type), intent(in) :: moved
      !> Displaced atom and component
      integer, intent(in) :: a, c
      !> Displacement
      real(wp), intent(in) :: shift
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: motion(3)
      integer :: i

      if (moved%ngrid /= grid%ngrid) then
         call test_failed(error, grid%kind_name()//": finite-difference step changed point membership")
         return
      end if
      do i = 1, grid%ngrid
         motion = 0.0_wp
         select type (grid)
         type is (moist_math_grid_3d_cartesian_type)
            motion(c) = shift/real(grid%natom, wp)
         type is (moist_math_grid_3d_molecular_type)
            if (grid%owner(i) == a) motion(c) = shift
         end select
         ! Per-element check calls cost ~30x here (millions of points); `.not. all(... < thr)` also fails on NaN
         if (.not. all(abs(moved%xyz(:, i) - grid%xyz(:, i) - motion) < 1.0e-12_wp)) then
            call test_failed(error, grid%kind_name()//": retained point identities must stay fixed")
            return
         end if
      end do
   end subroutine check_rigid_motion

   !> Linear synthetic energy in grid observables with fixed seeds
   !>
   !> @param[in] grid Realized grid
   !> @param[in] acc Fixed energy coefficients
   pure function grid_observable_energy(grid, acc) result(energy)
      !> Realized grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Fixed coefficients
      type(volume_adjoint_type), intent(in) :: acc
      !> Synthetic energy
      real(wp) :: energy

      energy = sum(acc%w_xyz*grid%xyz) + sum(acc%w_w*grid%w)
      if (allocated(grid%xi0)) energy = energy + sum(acc%w_xi*grid%xi0)
   end function grid_observable_energy

end module test_math_grid_3d_derivatives
