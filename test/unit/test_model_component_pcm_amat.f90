!> Interface-level unit tests and shared harness for the Gaussian PCM interaction matrix
module test_model_component_pcm_amat
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use moist_model_continuum_component_pcm_amat, only: assemble_pcm_amat, &
                                             pcm_amat_surface_weights, &
                                             assemble_pcm_amat_with_gradient, &
                                             pcm_amat_nuclear_gradient
   use moist_math_boys, only: boys_x_far
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_amat
   public :: count_branches
   public :: nmol, nleb_survey

   !> Number of mstore structures drawn per test
   integer, parameter :: nmol = 10

   !> Lebedev grid used for the survey tests
   integer, parameter :: nleb_survey = 50

   !> Routes of the interface-contract test, in call order
   integer, parameter :: nroute = 4
   character(len=*), parameter :: route_names(nroute) = [character(len=31) :: &
                                  "assemble_pcm_amat", "assemble_pcm_amat_with_gradient", &
                                  "pcm_amat_surface_weights", "pcm_amat_nuclear_gradient"]

   !> Malformed-argument groups, each shared by every route taking those arguments
   integer, parameter :: group_surface = 1, group_matrix = 2, group_derivative = 3, &
                         group_matrix_derivative = 4, group_contraction = 5, &
                         group_weights = 6, group_gradient = 7

   !> Argument extents and values of one interface-contract case, valid by default
   type :: interface_case_type
      !> Case description for failure messages
      character(len=64) :: label = "valid arguments"
      !> Width and switching-factor lengths
      integer :: nx = 3, nf = 3
      !> Position and nuclear-derivative extents (3 points, 2 atoms)
      integer :: zs(2) = [3, 3], dxs(3) = [3, 2, 3], dfs(3) = [3, 2, 3]
      integer :: dzs(4) = [3, 3, 2, 3]
      !> Matrix and matrix-derivative extents
      integer :: as(2) = [3, 3], das(4) = [3, 2, 3, 3]
      !> Contraction-vector, weight, and gradient extents
      integer :: nq(2) = [3, 3], nw(2) = [3, 3], wzs(2) = [3, 3], gs(2) = [3, 2]
      !> Width and switching factor written to every point
      real(wp) :: xi = 1.5_wp, f = 0.8_wp
   end type interface_case_type

contains

   !> Collect the interface-level test suite
   subroutine collect_model_component_pcm_amat(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rejects_invalid_surface", test_rejects_invalid_surface), &
                  new_unittest("interface_contracts", test_interface_contracts) &
                  ]

   end subroutine collect_model_component_pcm_amat

   !* --------------------------------- Shared helpers -------------------------------- *!

   !> Count the pairs on either side of the erf-saturation threshold
   !>
   !> @param[in]  xi    Gaussian widths (ngrid)
   !> @param[in]  xyz   Surface positions (3, ngrid)
   !> @param[out] nfar  Number of saturated pairs
   !> @param[out] nnear Number of unsaturated pairs
   subroutine count_branches(xi, xyz, nfar, nnear)
      !> Gaussian widths and surface positions
      real(wp), intent(in) :: xi(:), xyz(:, :)
      !> Saturated and unsaturated pair counts
      integer, intent(out) :: nfar, nnear

      !> Pair indices and grid size
      integer :: i, j, ngrid
      !> Squared separation
      real(wp) :: r2
      !> Per-point saturation bounds
      real(wp), allocatable :: bound(:)

      ngrid = size(xi)
      allocate (bound(ngrid))
      bound = boys_x_far/(xi*xi)

      nfar = 0
      nnear = 0
      do i = 1, ngrid
         do j = 1, i - 1
            r2 = sum((xyz(:, i) - xyz(:, j))**2)
            if (r2 >= bound(i) + bound(j)) then
               nfar = nfar + 1
            else
               nnear = nnear + 1
            end if
         end do
      end do

   end subroutine count_branches

   !* ------------------------------------- Tests ------------------------------------- *!

   !> Malformed surfaces are rejected and leave a defined output behind
   subroutine test_rejects_invalid_surface(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Synthetic surface
      real(wp) :: xi(3), f(3), xyz(3, 3), amat(3, 3), amat_small(2, 2)
      !> Contraction vectors and adjoints
      real(wp) :: q1(3), q2(3), w_xi(3), w_f(3), w_xyz(3, 3)

      xi = [1.5_wp, 2.0_wp, 0.9_wp]
      f = [0.8_wp, 0.95_wp, 0.6_wp]
      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.7_wp, -0.4_wp, 0.9_wp]
      xyz(:, 3) = [-3.1_wp, 2.2_wp, -1.4_wp]
      q1 = [0.3_wp, -0.2_wp, 0.1_wp]
      q2 = [0.1_wp, 0.4_wp, -0.3_wp]

      ! A non-positive width has no pair width, so the kernel is undefined
      amat = 1.0_wp
      call assemble_pcm_amat([1.5_wp, 0.0_wp, 0.9_wp], f, xyz, amat, err)
      call check(error, allocated(err), more="a zero Gaussian width was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(amat)), 0.0_wp, thr=0.0_wp, &
                 more="rejected assembly left the matrix untouched")
      if (allocated(error)) return
      deallocate (err)

      ! A non-positive switching factor divides the diagonal by zero
      amat = 1.0_wp
      call assemble_pcm_amat(xi, [0.8_wp, -0.1_wp, 0.6_wp], xyz, amat, err)
      call check(error, allocated(err), more="a negative switching factor was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(amat)), 0.0_wp, thr=0.0_wp, &
                 more="rejected assembly left the matrix untouched")
      if (allocated(error)) return
      deallocate (err)

      ! The output matrix must match the surface it is assembled on
      amat_small = 1.0_wp
      call assemble_pcm_amat(xi, f, xyz, amat_small, err)
      call check(error, allocated(err), more="a mis-shaped output matrix was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(amat_small)), 0.0_wp, thr=0.0_wp, &
                 more="rejected assembly left the matrix untouched")
      if (allocated(error)) return
      deallocate (err)

      ! The adjoint route validates the same surface plus its own vectors
      w_xi = 1.0_wp
      w_f = 1.0_wp
      w_xyz = 1.0_wp
      call pcm_amat_surface_weights(xi, f, xyz, q1(1:2), q2, w_xi, w_f, w_xyz, err)
      call check(error, allocated(err), &
                 more="a mis-shaped contraction vector was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(w_xi)) + maxval(abs(w_f)) + maxval(abs(w_xyz)), &
                 0.0_wp, thr=0.0_wp, more="rejected contraction left the weights untouched")

   end subroutine test_rejects_invalid_surface

   !> Reject each independently malformed argument and zero every output
   !>
   !> Each route first accepts the valid default case, then walks its argument
   !> groups until a group runs out of cases, so a new case needs no count update
   subroutine test_interface_contracts(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Argument groups of the current route
      integer, allocatable :: groups(:)
      !> Route, group, and in-group case indices
      integer :: route, igroup, k
      !> Current case
      type(interface_case_type) :: this_case
      !> Whether the group has a case k
      logical :: found

      do route = 1, nroute
         select case (route)
         case (1)
            groups = [group_surface, group_matrix]
         case (2)
            groups = [group_surface, group_derivative, group_matrix, &
                      group_matrix_derivative]
         case (3)
            groups = [group_surface, group_contraction, group_weights]
         case (4)
            groups = [group_derivative, group_weights, group_gradient]
         case default
            call check(error, .false., more="invalid interface test route")
            return
         end select

         call run_interface_case(error, route, interface_case_type(), .true.)
         if (allocated(error)) return
         do igroup = 1, size(groups)
            k = 0
            do
               k = k + 1
               call get_interface_case(groups(igroup), k, this_case, found)
               if (.not. found) exit
               call run_interface_case(error, route, this_case, .false.)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_interface_contracts

   !> Malformed-argument case k of one argument group
   !>
   !> @param[in]  group     Argument group
   !> @param[in]  k         Case index within the group
   !> @param[out] this_case Case extents and values, valid except for one argument
   !> @param[out] found     Whether the group has a case k
   subroutine get_interface_case(group, k, this_case, found)
      !> Argument group and case index within it
      integer, intent(in) :: group, k
      !> Case extents and values
      type(interface_case_type), intent(out) :: this_case
      !> Whether the group has a case k
      logical, intent(out) :: found

      found = .true.
      select case (group)
      case (group_surface)
         select case (k)
         case (1)
            this_case%label = "empty surface"
            this_case%nx = 0
            this_case%nf = 0
            this_case%zs(2) = 0
            this_case%as = 0
            this_case%das(3:4) = 0
            this_case%nq = 0
            this_case%nw = 0
            this_case%wzs(2) = 0
            this_case%dxs(3) = 0
            this_case%dfs(3) = 0
            this_case%dzs(4) = 0
         case (2)
            this_case%label = "switching factors longer than the widths"
            this_case%nf = 4
         case (3)
            this_case%label = "positions with 4 Cartesian rows"
            this_case%zs(1) = 4
         case (4)
            this_case%label = "positions longer than the widths"
            this_case%zs(2) = 4
         case (5)
            this_case%label = "zero Gaussian width"
            this_case%xi = 0.0_wp
         case (6)
            this_case%label = "negative Gaussian width"
            this_case%xi = -0.1_wp
         case (7)
            this_case%label = "zero switching factor"
            this_case%f = 0.0_wp
         case (8)
            this_case%label = "negative switching factor"
            this_case%f = -0.1_wp
         case default
            found = .false.
         end select
      case (group_matrix)
         select case (k)
         case (1)
            this_case%label = "matrix with an extra row"
            this_case%as(1) = 4
         case (2)
            this_case%label = "matrix with an extra column"
            this_case%as(2) = 4
         case default
            found = .false.
         end select
      case (group_derivative)
         select case (k)
         case (1)
            this_case%label = "width derivative with 4 Cartesian rows"
            this_case%dxs(1) = 4
         case (2)
            this_case%label = "switching derivative with 4 Cartesian rows"
            this_case%dfs(1) = 4
         case (3)
            this_case%label = "position derivative with 4 position rows"
            this_case%dzs(1) = 4
         case (4)
            this_case%label = "position derivative with 4 Cartesian axes"
            this_case%dzs(2) = 4
         case (5)
            this_case%label = "switching derivative with an extra atom"
            this_case%dfs(2) = 3
         case (6)
            this_case%label = "position derivative with an extra atom"
            this_case%dzs(3) = 3
         case (7)
            this_case%label = "width derivative with an extra point"
            this_case%dxs(3) = 4
         case (8)
            this_case%label = "switching derivative with an extra point"
            this_case%dfs(3) = 4
         case (9)
            this_case%label = "position derivative with an extra point"
            this_case%dzs(4) = 4
         case default
            found = .false.
         end select
      case (group_matrix_derivative)
         select case (k)
         case (1)
            this_case%label = "matrix derivative with 4 Cartesian rows"
            this_case%das(1) = 4
         case (2)
            this_case%label = "matrix derivative with an extra atom"
            this_case%das(2) = 3
         case (3)
            this_case%label = "matrix derivative with an extra row"
            this_case%das(3) = 4
         case (4)
            this_case%label = "matrix derivative with an extra column"
            this_case%das(4) = 4
         case default
            found = .false.
         end select
      case (group_contraction)
         select case (k)
         case (1)
            this_case%label = "left contraction vector with an extra point"
            this_case%nq(1) = 4
         case (2)
            this_case%label = "right contraction vector with an extra point"
            this_case%nq(2) = 4
         case (3)
            this_case%label = "width weights with an extra point"
            this_case%nw(1) = 4
         case default
            found = .false.
         end select
      case (group_weights)
         select case (k)
         case (1)
            this_case%label = "switching weights with an extra point"
            this_case%nw(2) = 4
         case (2)
            this_case%label = "position weights with 4 Cartesian rows"
            this_case%wzs(1) = 4
         case (3)
            this_case%label = "position weights with an extra point"
            this_case%wzs(2) = 4
         case default
            found = .false.
         end select
      case (group_gradient)
         select case (k)
         case (1)
            this_case%label = "gradient with 4 Cartesian rows"
            this_case%gs(1) = 4
         case (2)
            this_case%label = "gradient with an extra atom"
            this_case%gs(2) = 3
         case default
            found = .false.
         end select
      case default
         found = .false.
      end select

   end subroutine get_interface_case

   !> Call one route on one case and check acceptance, or rejection with zeroed outputs
   !>
   !> @param[out] error     Test failure
   !> @param[in]  route     Route index into route_names
   !> @param[in]  this_case Case extents and values
   !> @param[in]  valid     Whether the case must be accepted
   subroutine run_interface_case(error, route, this_case, valid)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Route index
      integer, intent(in) :: route
      !> Case extents and values
      type(interface_case_type), intent(in) :: this_case
      !> Whether the case must be accepted
      logical, intent(in) :: valid

      !> Library error
      type(moist_error_type), allocatable :: err
      !> Sum of absolute outputs after the call
      real(wp) :: residual
      !> Surface, derivative, contraction, and output arrays
      real(wp), allocatable :: xi(:), f(:), xyz(:, :), dx(:, :, :), df(:, :, :)
      real(wp), allocatable :: dz(:, :, :, :), a(:, :), da(:, :, :, :)
      real(wp), allocatable :: q1(:), q2(:), wx(:), wf(:), wz(:, :), grad(:, :)
      !> Failure-message prefix
      character(len=128) :: label

      allocate (xi(this_case%nx), f(this_case%nf), xyz(this_case%zs(1), this_case%zs(2)))
      allocate (dx(this_case%dxs(1), this_case%dxs(2), this_case%dxs(3)))
      allocate (df(this_case%dfs(1), this_case%dfs(2), this_case%dfs(3)))
      allocate (dz(this_case%dzs(1), this_case%dzs(2), this_case%dzs(3), this_case%dzs(4)))
      allocate (a(this_case%as(1), this_case%as(2)))
      allocate (da(this_case%das(1), this_case%das(2), this_case%das(3), this_case%das(4)))
      allocate (q1(this_case%nq(1)), q2(this_case%nq(2)))
      allocate (wx(this_case%nw(1)), wf(this_case%nw(2)))
      allocate (wz(this_case%wzs(1), this_case%wzs(2)))
      allocate (grad(this_case%gs(1), this_case%gs(2)))
      xi = this_case%xi
      f = this_case%f
      xyz = 0.0_wp
      dx = 0.1_wp
      df = 0.2_wp
      dz = 0.3_wp
      q1 = 0.3_wp
      q2 = -0.2_wp
      a = 1.0_wp
      da = 1.0_wp
      wx = 1.0_wp
      wf = 1.0_wp
      wz = 1.0_wp
      grad = 1.0_wp
      select case (route)
      case (1)
         call assemble_pcm_amat(xi, f, xyz, a, err)
         residual = sum(abs(a))
      case (2)
         call assemble_pcm_amat_with_gradient(xi, f, xyz, dx, df, dz, a, da, err)
         residual = sum(abs(a)) + sum(abs(da))
      case (3)
         call pcm_amat_surface_weights(xi, f, xyz, q1, q2, wx, wf, wz, err)
         residual = sum(abs(wx)) + sum(abs(wf)) + sum(abs(wz))
      case (4)
         call pcm_amat_nuclear_gradient(dx, df, dz, wx, wf, wz, grad, err)
         residual = sum(abs(grad))
      case default
         call check(error, .false., more="invalid interface test route")
         return
      end select

      label = trim(route_names(route))//": "//trim(this_case%label)
      if (valid) then
         call check(error, .not. allocated(err), more=trim(label)//" rejected")
      else
         call check(error, allocated(err), more=trim(label)//" accepted")
         if (allocated(error)) return
         call check(error, residual, 0.0_wp, thr=0.0_wp, &
                    more=trim(label)//" left a nonzero output")
      end if

   end subroutine run_interface_case

end module test_model_component_pcm_amat
