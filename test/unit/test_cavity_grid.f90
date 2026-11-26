!> Grid components of the Lebedev-based cavities: container and split spheres
!>
!> A mock grid model gives every sphere several unrotated components whose
!> weights sum to one; a cavity that loops over the components correctly then
!> reproduces the plain cavity (areas, volume and their gradients), since
!> every node appears once per component with its share of the weight
!>
!> A second mock gives every sphere two components whose weights, and the
!> rotation of the heavier one, depend smoothly on the coordinates: the
!> gradients of the cavities must then carry the weight and rotation terms
!> of the model (central differences)
module test_cavity_grid
   use mctc_env, only: wp
   use moist_math_linalg_geometry, only: identity3
   use mctc_env_error, only: mctc_error => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_io, only: structure_type
   use moist_cavity_adjoint, only: cavity_adjoint_type
   use moist_cavity_type, only: cavity_type
   use mstore, only: get_structure
   use moist_context, only: moist_context_type, new_context
   use moist_channels_fields, only: field_query_type
   use moist_radii, only: default_cpcm_radii
   use moist_cavity_grid, only: grid_model_type, grid_components_type, grid_size_default, cavity_grid_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_cavity_numsa, only: cavity_type_numsa, new_cavity_numsa
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use test_cavity_drop_fd, only: drop_fd_param, drop_fd_svdw
   implicit none(type, external)
   private

   public :: collect_cavity_grid

   !> Mock model: `weight` gives the components of every sphere, all unrotated
   type, extends(grid_model_type) :: split_model_type
      !> Weights of the components of every sphere
      real(wp), allocatable :: weight(:)
      !> Lebedev size asked for
      integer :: size = grid_size_default
   contains
      procedure :: update => split_update
      procedure :: is_current => split_is_current
      procedure :: is_trivial => split_is_trivial
      procedure :: is_plain => split_is_plain
      procedure :: components => split_components
      procedure :: add_gradient => split_add_gradient
      procedure :: jacobian => split_jacobian
      procedure :: list_fields => split_list_fields
      procedure :: print => split_print
   end type split_model_type

   !> Mock model with weights and a rotation that follow the geometry
   !>
   !> - with s = sum(coef*xyz): component 1 of every sphere unrotated with
   !>   weight 0.3 + 0.1 sin(s); component 2 turned by 0.4 + 0.3 cos(s) about
   !>   z with weight 0.8 - half the first weight, always the heavier
   !> - `single`: component 2 alone, one per sphere with a weight below one
   type, extends(grid_model_type) :: fade_model_type
      !> Coordinates of the latest update (3, nat)
      real(wp), allocatable :: xyz(:, :)
      !> Whether every sphere carries the second component alone
      logical :: single = .false.
   contains
      procedure :: update => fade_update
      procedure :: is_current => fade_is_current
      procedure :: is_trivial => fade_is_trivial
      procedure :: is_plain => fade_is_plain
      procedure :: components => fade_components
      procedure :: add_gradient => fade_add_gradient
      procedure :: jacobian => fade_jacobian
      procedure :: list_fields => fade_list_fields
      procedure :: print => fade_print
      procedure :: state => fade_state
   end type fade_model_type

   !> Cavities of the gradient checks
   integer, parameter :: kind_iswig = 1, kind_numsa = 2, kind_drop = 3

contains

   !> Collect the grid component tests
   !>
   !> @param[out] testsuite Tests
   subroutine collect_cavity_grid(testsuite)
      !> Tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("container_plain", test_container_plain), &
         & new_unittest("container_split", test_container_split), &
         & new_unittest("other_size_rejected", test_other_size), &
         & new_unittest("iswig_split", test_iswig_split), &
         & new_unittest("numsa_split", test_numsa_split), &
         & new_unittest("drop_heaviest_component", test_drop_split), &
         & new_unittest("heaviest_contract", test_heaviest_contract), &
         & new_unittest("weight_gradient_iswig", test_fade_iswig), &
         & new_unittest("weight_gradient_numsa", test_fade_numsa), &
         & new_unittest("weight_gradient_drop_reduced", test_fade_drop), &
         & new_unittest("drop_single_weighted_component", test_fade_drop_single) &
         & ]
   end subroutine collect_cavity_grid

   !* ------------------------------------ Mock model ------------------------------------ *!

   !> Mock update, no work
   subroutine split_update(self, mol, radii, error, ctx)
      class(split_model_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)
      type(mctc_error), allocatable, intent(out) :: error
      type(moist_context_type), intent(inout), optional :: ctx
   end subroutine split_update

   !> Mock always current
   function split_is_current(self, mol, radii) result(current)
      class(split_model_type), intent(in) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)
      logical :: current
      current = .true.
   end function split_is_current

   !> Mock asks for nothing
   pure function split_is_trivial(self) result(yes)
      class(split_model_type), intent(in) :: self
      logical :: yes
      yes = .true.
   end function split_is_trivial

   !> Mock components are not plain
   pure function split_is_plain(self) result(yes)
      class(split_model_type), intent(in) :: self
      logical :: yes
      yes = .false.
   end function split_is_plain

   !> Fixed weighted components of fixed sizes on every sphere
   subroutine split_components(self, nsph, comps)
      class(split_model_type), intent(in) :: self
      integer, intent(in) :: nsph
      type(grid_components_type), intent(out) :: comps
      integer :: isph, k, nk

      nk = size(self%weight)
      allocate (comps%sph_first(nsph + 1), comps%weight(nsph*nk))
      allocate (comps%size(nsph*nk), source=self%size)
      allocate (comps%rot(nsph*nk), source=1)
      do isph = 1, nsph + 1
         comps%sph_first(isph) = (isph - 1)*nk + 1
      end do
      do isph = 1, nsph
         do k = 1, nk
            comps%weight((isph - 1)*nk + k) = self%weight(k)
         end do
      end do
      allocate (comps%rotation(3, 3, 1), source=0.0_wp)
      do k = 1, 3
         comps%rotation(k, k, 1) = 1.0_wp
      end do
   end subroutine split_components

   !> No gradient: fixed rotation and weights
   subroutine split_add_gradient(self, torque, wbar, gradient)
      class(split_model_type), intent(in) :: self
      real(wp), intent(in) :: torque(:, :)
      real(wp), intent(in) :: wbar(:)
      real(wp), intent(inout) :: gradient(:, :)
   end subroutine split_add_gradient

   !> Zero rotation Jacobian, no weight Jacobian
   subroutine split_jacobian(self, jac_rot, jac_weight, weight_class)
      class(split_model_type), intent(in) :: self
      real(wp), allocatable, intent(out) :: jac_rot(:, :, :, :)
      real(wp), allocatable, intent(out) :: jac_weight(:, :, :)
      integer, allocatable, intent(out) :: weight_class(:)
      ! Fixed rotation and weights; the tests size it for their structure
      allocate (jac_rot(3, 3, 16, 1), source=0.0_wp)
   end subroutine split_jacobian

   !> Mock declares no fields
   subroutine split_list_fields(self, query)
      class(split_model_type), intent(in) :: self
      type(field_query_type), intent(inout) :: query
   end subroutine split_list_fields

   !> Mock prints nothing
   subroutine split_print(self, unit, level)
      class(split_model_type), intent(in) :: self
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: level
   end subroutine split_print

   !* ----------------------------- Geometry-dependent mock ------------------------------ *!

   !> Weights, angle and their derivatives of the fading mock
   !>
   !> @param[in]  self  Mock model after an update
   !> @param[out] w     Weights of the two components
   !> @param[out] theta Angle of the second component about z
   !> @param[out] ds    ds/dx (3, nat)
   !> @param[out] dw    dw/ds of the two weights
   !> @param[out] dth   dtheta/ds
   pure subroutine fade_state(self, w, theta, ds, dw, dth)
      !> Mock model
      class(fade_model_type), intent(in) :: self
      !> Weights
      real(wp), intent(out) :: w(2)
      !> Angle
      real(wp), intent(out) :: theta
      !> Derivative of s
      real(wp), allocatable, intent(out) :: ds(:, :)
      !> Weight derivatives
      real(wp), intent(out) :: dw(2)
      !> Angle derivative
      real(wp), intent(out) :: dth

      real(wp) :: s
      integer :: ia, nat

      nat = size(self%xyz, 2)
      allocate (ds(3, nat))
      do ia = 1, nat
         ds(:, ia) = 0.2_wp*sin(1.3_wp*real(ia, wp) + [0.4_wp, 2.1_wp, 3.3_wp])
      end do
      s = sum(ds*self%xyz)
      w(1) = 0.3_wp + 0.1_wp*sin(s)
      w(2) = 0.8_wp - 0.5_wp*w(1)
      dw(1) = 0.1_wp*cos(s)
      dw(2) = -0.5_wp*dw(1)
      theta = 0.4_wp + 0.3_wp*cos(s)
      dth = -0.3_wp*sin(s)
   end subroutine fade_state

   !> Mock update, keep the positions
   subroutine fade_update(self, mol, radii, error, ctx)
      class(fade_model_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)
      type(mctc_error), allocatable, intent(out) :: error
      type(moist_context_type), intent(inout), optional :: ctx
      self%xyz = mol%xyz
   end subroutine fade_update

   !> Mock never current
   function fade_is_current(self, mol, radii) result(current)
      class(fade_model_type), intent(in) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)
      logical :: current
      current = .false.
   end function fade_is_current

   !> Mock always does work
   pure function fade_is_trivial(self) result(yes)
      class(fade_model_type), intent(in) :: self
      logical :: yes
      yes = .false.
   end function fade_is_trivial

   !> Mock components are not plain
   pure function fade_is_plain(self) result(yes)
      class(fade_model_type), intent(in) :: self
      logical :: yes
      yes = .false.
   end function fade_is_plain

   !> Two components per sphere, or the turned one alone
   subroutine fade_components(self, nsph, comps)
      class(fade_model_type), intent(in) :: self
      integer, intent(in) :: nsph
      type(grid_components_type), intent(out) :: comps
      real(wp), allocatable :: ds(:, :)
      real(wp) :: w(2), theta, dw(2), dth
      integer :: isph, nk

      call self%state(w, theta, ds, dw, dth)
      nk = merge(1, 2, self%single)
      allocate (comps%sph_first(nsph + 1), comps%weight(nsph*nk), comps%rot(nsph*nk))
      allocate (comps%size(nsph*nk), source=grid_size_default)
      comps%sph_first = [((isph - 1)*nk + 1, isph=1, nsph + 1)]
      if (self%single) then
         comps%weight = w(2)
         comps%rot = 2
      else
         comps%weight = [(w, isph=1, nsph)]
         comps%rot = [([1, 2], isph=1, nsph)]
      end if
      allocate (comps%rotation(3, 3, 2), source=0.0_wp)
      comps%rotation(:, :, 1) = identity3
      comps%rotation(:, :, 2) = reshape([cos(theta), sin(theta), 0.0_wp, -sin(theta), cos(theta), 0.0_wp, &
         & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
   end subroutine fade_components

   !> Map weight adjoints and z torques to the atoms
   subroutine fade_add_gradient(self, torque, wbar, gradient)
      class(fade_model_type), intent(in) :: self
      real(wp), intent(in) :: torque(:, :)
      real(wp), intent(in) :: wbar(:)
      real(wp), intent(inout) :: gradient(:, :)
      real(wp), allocatable :: ds(:, :)
      real(wp) :: w(2), theta, dw(2), dth
      integer :: icomp, k

      call self%state(w, theta, ds, dw, dth)
      do icomp = 1, size(wbar)
         k = merge(2, 2 - mod(icomp, 2), self%single)
         gradient = gradient + wbar(icomp)*dw(k)*ds
         if (k == 2) gradient = gradient + torque(3, icomp)*dth*ds
      end do
   end subroutine fade_add_gradient

   !> Rotation and weight Jacobians of the fading mock
   subroutine fade_jacobian(self, jac_rot, jac_weight, weight_class)
      class(fade_model_type), intent(in) :: self
      real(wp), allocatable, intent(out) :: jac_rot(:, :, :, :)
      real(wp), allocatable, intent(out) :: jac_weight(:, :, :)
      integer, allocatable, intent(out) :: weight_class(:)
      real(wp), allocatable :: ds(:, :)
      real(wp) :: w(2), theta, dw(2), dth
      integer :: nat, isph

      call self%state(w, theta, ds, dw, dth)
      nat = size(ds, 2)
      allocate (jac_rot(3, 3, nat, 2), source=0.0_wp)
      jac_rot(3, :, :, 2) = dth*ds
      allocate (jac_weight(3, nat, 2))
      jac_weight(:, :, 1) = dw(1)*ds
      jac_weight(:, :, 2) = dw(2)*ds
      if (self%single) then
         weight_class = [(2, isph=1, nat)]
      else
         weight_class = [([1, 2], isph=1, nat)]
      end if
   end subroutine fade_jacobian

   !> Mock declares no fields
   subroutine fade_list_fields(self, query)
      class(fade_model_type), intent(in) :: self
      type(field_query_type), intent(inout) :: query
   end subroutine fade_list_fields

   !> Mock prints nothing
   subroutine fade_print(self, unit, level)
      class(fade_model_type), intent(in) :: self
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: level
   end subroutine fade_print

   !> Updated cavity of one kind with the fading mock
   !>
   !> - DROP: the settings of the per-point DROP checks
   !>
   !> @param[in]  kind   Cavity kind
   !> @param[in]  mol    Structure
   !> @param[in]  ctx    Run context
   !> @param[out] cav    Updated cavity
   !> @param[out] err    Construction or update error
   !> @param[in]  single One weighted component per sphere (default two)
   !> @param[in]  plain  No model at all (default the mock)
   subroutine fade_cavity(kind, mol, ctx, cav, err, single, plain)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Updated cavity
      class(cavity_type), allocatable, intent(out) :: cav
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: err
      !> One weighted component per sphere
      logical, intent(in), optional :: single
      !> No model
      logical, intent(in), optional :: plain

      type(fade_model_type), allocatable :: model

      allocate (model)
      if (present(single)) model%single = single
      if (present(plain)) then
         if (plain) deallocate (model)
      end if
      select case (kind)
      case (kind_iswig)
         block
            type(cavity_type_iswig), allocatable :: iswig
            allocate (iswig)
            call new_cavity_iswig(iswig, default_cpcm_radii(), err, ctx=ctx, grid_model=model)
            call move_alloc(iswig, cav)
         end block
      case (kind_numsa)
         block
            type(cavity_type_numsa), allocatable :: numsa
            allocate (numsa)
            call new_cavity_numsa(numsa, default_cpcm_radii(), err, ctx=ctx, grid_model=model)
            call move_alloc(numsa, cav)
         end block
      case default
         block
            type(cavity_type_drop), allocatable :: drop
            allocate (drop)
            call new_cavity_drop(drop, default_cpcm_radii(), drop_fd_svdw(), err, param=drop_fd_param(), &
               & ctx=ctx, grid_model=model)
            call move_alloc(drop, cav)
         end block
      end select
      if (.not. allocated(err)) call cav%update(mol, err)
   end subroutine fade_cavity

   !> Area and volume of a cavity with the fading mock, and their gradients
   !>
   !> - gradients: forward area, reverse area (unit area adjoint), forward
   !>   volume; NumSA: forward area and the sum of `dsdr` over the spheres
   !>
   !> @param[in]  kind  Cavity kind
   !> @param[in]  mol   Structure
   !> @param[in]  ctx   Run context
   !> @param[out] vals  Area, area, volume (NumSA: area twice, zero)
   !> @param[out] error Test failure
   !> @param[out] grads Gradients (3, nat, 3)
   subroutine fade_values(kind, mol, ctx, vals, error, grads)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Values
      real(wp), intent(out) :: vals(3)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Gradients
      real(wp), intent(out), optional :: grads(:, :, :)

      class(cavity_type), allocatable :: cav
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: err

      vals = 0.0_wp
      call fade_cavity(kind, mol, ctx, cav, err)
      if (.not. allocated(err) .and. present(grads)) call cav%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, "cavity failed: "//err%message)
         return
      end if
      vals(1:2) = cav%total_area
      if (kind /= kind_numsa) vals(3) = cav%total_volume
      if (.not. present(grads)) return
      grads = 0.0_wp
      select type (cav)
      type is (cavity_type_iswig)
         grads(:, :, 1) = cav%area_grad
         grads(:, :, 3) = cav%volume_grad
      type is (cavity_type_numsa)
         grads(:, :, 1) = cav%area_grad
         grads(:, :, 2) = sum(cav%dsdr, dim=3)
      type is (cavity_type_drop)
         grads(:, :, 1) = cav%A_tot1_rA
         grads(:, :, 3) = cav%V_tot1_rA
      end select
      if (kind == kind_numsa) return
      call acc%init(cav%ngrid)
      acc%w_a = 1.0_wp
      call cav%get_surface_gradient(acc, grads(:, :, 2), err)
      if (allocated(err)) call test_failed(error, "reverse gradient failed: "//err%message)
   end subroutine fade_values

   !> Gradients of a cavity with the fading mock against central differences
   !>
   !> @param[in]  kind  Cavity kind
   !> @param[in]  name  MB16-43 record
   !> @param[in]  step  Step, bohr
   !> @param[in]  thr   Tolerance relative to the largest gradient component
   !> @param[out] error Test failure
   subroutine check_fade_gradient(kind, name, step, thr, error)
      !> Cavity kind
      integer, intent(in) :: kind
      !> Record
      character(len=*), intent(in) :: name
      !> Step
      real(wp), intent(in) :: step
      !> Tolerance
      real(wp), intent(in) :: thr
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: label(3) = [character(len=14) :: "forward area", "second area", &
         & "forward volume"]
      type(structure_type) :: mol, moved
      type(moist_context_type), target :: ctx
      real(wp), allocatable :: grads(:, :, :), dir(:, :), plain_grad(:, :, :)
      real(wp) :: vals(3), vp(3), vm(3), analytic, numeric, scale
      character(len=160) :: text
      integer :: idir, iq, ia

      call get_structure(mol, "MB16-43", name)
      call new_context(ctx, nthreads=0, verbosity=0)
      allocate (grads(3, mol%nat, 3), plain_grad(3, mol%nat, 3), dir(3, mol%nat))
      call fade_values(kind, mol, ctx, vals, error, grads)
      if (allocated(error)) return
      do idir = 1, 3
         do ia = 1, mol%nat
            dir(:, ia) = sin(real(idir, wp)*2.1_wp + 0.7_wp*real(ia, wp) + [0.3_wp, 1.7_wp, 4.0_wp])
         end do
         moved = mol
         moved%xyz = mol%xyz + step*dir
         call fade_values(kind, moved, ctx, vp, error)
         if (allocated(error)) return
         moved%xyz = mol%xyz - step*dir
         call fade_values(kind, moved, ctx, vm, error)
         if (allocated(error)) return
         do iq = 1, 3
            if (kind == kind_numsa .and. iq == 3) cycle
            analytic = sum(grads(:, :, iq)*dir)
            numeric = (vp(iq) - vm(iq))/(2.0_wp*step)
            scale = max(1.0_wp, maxval(abs(grads(:, :, iq))))
            write (text, '(a,", direction ",i0,": ",es14.6," vs FD ",es14.6)') trim(label(iq)), idir, &
               & analytic, numeric
            call check(error, abs(analytic - numeric) <= thr*scale, &
               & "gradient with geometry-dependent weights differs: "//trim(text))
            if (allocated(error)) return
         end do
      end do
   end subroutine check_fade_gradient

   !* --------------------------------------- Tests -------------------------------------- *!

   !> Without a model: one plain component per sphere on the unit grid itself
   subroutine test_container_plain(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_grid_type) :: grid
      type(mctc_error), allocatable :: err
      real(wp) :: unit_grid(3, 4), unit_weight(4)

      call random_number(unit_grid)
      unit_weight = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      call grid%assemble(3, unit_grid, unit_weight, 0.5_wp, err)
      call check(error,.not. allocated(err), "plain assembly failed")
      if (allocated(error)) return
      call check(error, grid%plain .and. grid%ncomp == 3 .and. grid%nnode == 12 .and. grid%nset == 1, &
         & "plain counts")
      if (allocated(error)) return
      call check(error, all(grid%sph_first == [1, 2, 3, 4]) .and. all(grid%comp_first == [1, 5, 9, 13]) &
         & .and. all(grid%comp_sphere == [1, 2, 3]) .and. all(grid%comp_weight == 1.0_wp) &
         & .and. all(grid%comp_zeta == 0.5_wp), "plain ranges")
      if (allocated(error)) return
      call check(error, all(grid%node == unit_grid) .and. all(grid%wang == unit_weight) .and. &
         & grid%set_first(grid%comp_set(3)) == 1 .and. grid%comp_first(3) - grid%comp_first(2) == 4, "plain nodes")
   end subroutine test_container_plain

   !> A model with three components per sphere: contiguous ranges per sphere
   subroutine test_container_split(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_grid_type) :: grid
      type(split_model_type) :: model
      type(mctc_error), allocatable :: err
      real(wp) :: unit_grid(3, 4), unit_weight(4)

      call random_number(unit_grid)
      unit_weight = 1.0_wp
      model%weight = [0.2_wp, 0.3_wp, 0.5_wp]
      call grid%assemble(2, unit_grid, unit_weight, 0.5_wp, err, model)
      call check(error,.not. allocated(err), "split assembly failed")
      if (allocated(error)) return
      call check(error,.not. grid%plain .and. grid%ncomp == 6 .and. grid%nnode == 24 .and. grid%nset == 1, &
         & "split counts")
      if (allocated(error)) return
      call check(error, all(grid%sph_first == [1, 4, 7]) .and. all(grid%comp_sphere == [1, 1, 1, 2, 2, 2]) &
         & .and. all(grid%comp_first == [1, 5, 9, 13, 17, 21, 25]) .and. all(grid%comp_set == 1), "split ranges")
      if (allocated(error)) return
      call check(error, all(grid%comp_weight == [0.2_wp, 0.3_wp, 0.5_wp, 0.2_wp, 0.3_wp, 0.5_wp]) .and. &
         & maxval(abs(grid%node - unit_grid)) < 1.0e-15_wp, "split weights or nodes")
   end subroutine test_container_split

   !> A Lebedev size other than the cavity's is refused for now
   subroutine test_other_size(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_grid_type) :: grid
      type(split_model_type) :: model
      type(mctc_error), allocatable :: err
      real(wp) :: unit_grid(3, 4), unit_weight(4)

      unit_grid = 0.0_wp
      unit_weight = 1.0_wp
      model%weight = [1.0_wp]
      model%size = 50
      call grid%assemble(2, unit_grid, unit_weight, 0.5_wp, err, model)
      call check(error, allocated(err), "another Lebedev size was accepted")
   end subroutine test_other_size

   !> iSwiG with split spheres reproduces the plain cavity
   subroutine test_iswig_split(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_context_type), target :: ctx
      type(cavity_type_iswig) :: plain, split
      type(split_model_type) :: model
      type(mctc_error), allocatable :: err

      call get_structure(mol, "MB16-43", "01")
      call new_context(ctx, verbosity=0)
      model%weight = [0.3_wp, 0.7_wp]
      call new_cavity_iswig(plain, default_cpcm_radii(), err, ctx=ctx)
      if (.not. allocated(err)) call plain%update(mol, err)
      if (.not. allocated(err)) call plain%get_gradient(err)
      if (.not. allocated(err)) call new_cavity_iswig(split, default_cpcm_radii(), err, ctx=ctx, grid_model=model)
      if (.not. allocated(err)) call split%update(mol, err)
      if (.not. allocated(err)) call split%get_gradient(err)
      call check(error,.not. allocated(err), "iSwiG failed")
      if (allocated(error)) return

      call check(error, split%grid%ncomp == 2*mol%nat .and. split%ngrid > plain%ngrid .and. &
         & all(split%grid%comp_sphere(split%comp) == split%owner), "split points or components")
      if (allocated(error)) return
      ! The switching cutoff drops a few more of the down-weighted points
      call check(error, abs(split%total_area - plain%total_area) < 1.0e-7_wp*plain%total_area .and. &
         & abs(split%total_volume - plain%total_volume) < 1.0e-7_wp*plain%total_volume .and. &
         & maxval(abs(split%asph - plain%asph)) < 1.0e-7_wp*plain%total_area, "split area or volume")
      if (allocated(error)) return
      call check(error, maxval(abs(split%area_grad - plain%area_grad)) < 1.0e-6_wp .and. &
         & maxval(abs(split%volume_grad - plain%volume_grad)) < 1.0e-6_wp, "split gradients")
   end subroutine test_iswig_split

   !> NumSA with split spheres reproduces the plain areas and gradient
   subroutine test_numsa_split(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_context_type), target :: ctx
      type(cavity_type_numsa) :: plain, split
      type(split_model_type) :: model
      type(mctc_error), allocatable :: err

      call get_structure(mol, "MB16-43", "01")
      call new_context(ctx, verbosity=0)
      model%weight = [0.25_wp, 0.25_wp, 0.5_wp]
      call new_cavity_numsa(plain, default_cpcm_radii(), err, ctx=ctx)
      if (.not. allocated(err)) call plain%update(mol, err)
      if (.not. allocated(err)) call new_cavity_numsa(split, default_cpcm_radii(), err, ctx=ctx, grid_model=model)
      if (.not. allocated(err)) call split%update(mol, err)
      call check(error,.not. allocated(err), "NumSA failed")
      if (allocated(error)) return

      call check(error, maxval(abs(split%asph - plain%asph)) < 1.0e-11_wp*maxval(plain%asph) .and. &
         & maxval(abs(split%area_grad - plain%area_grad)) < 1.0e-10_wp, "split areas or gradient")
   end subroutine test_numsa_split

   !> DROP keeps the heaviest component of every sphere, with weight one
   subroutine test_drop_split(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_context_type), target :: ctx
      type(cavity_type_drop) :: cav
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_parameters_type) :: param
      type(split_model_type) :: model
      type(mctc_error), allocatable :: err

      call get_structure(mol, "MB16-43", "01")
      call new_context(ctx, verbosity=0)
      model%weight = [0.5_wp, 0.5_wp]
      call param%init_defaults()
      call lsf%new()
      call new_cavity_drop(cav, default_cpcm_radii(), lsf, err, param=param, ctx=ctx, grid_model=model)
      if (.not. allocated(err)) call cav%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "DROP failed: "//err%message)
         return
      end if
      call check(error, cav%grid%reduced .and. cav%grid%ncomp == cav%nsph .and. &
         & all(cav%grid%comp_weight == 1.0_wp), "DROP did not keep one component of weight one per sphere")
   end subroutine test_drop_split

   !> `heaviest`: one component of weight one per sphere, whatever the model gives
   !>
   !> - a sphere without components is refused
   !> - one weighted component per sphere is kept with weight one, reduced
   !> - one component of weight one per sphere needs no reduction
   subroutine test_heaviest_contract(error)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_grid_type) :: grid
      type(split_model_type) :: split
      type(fade_model_type) :: fade
      type(mctc_error), allocatable :: err
      real(wp) :: unit_grid(3, 4), unit_weight(4)

      call random_number(unit_grid)
      unit_weight = 1.0_wp
      fade%single = .true.
      allocate (fade%xyz(3, 2), source=0.3_wp)
      call grid%assemble(2, unit_grid, unit_weight, 0.5_wp, err, fade, heaviest=.true.)
      call check(error,.not. allocated(err), "single weighted component refused")
      if (allocated(error)) return
      call check(error, grid%reduced .and. grid%ncomp == 2 .and. all(grid%comp_weight == 1.0_wp) &
         & .and. all(grid%comp_model == [1, 2]) .and. grid%ncomp_model == 2, &
         & "a weighted component per sphere was not taken with weight one")
      if (allocated(error)) return

      split%weight = [1.0_wp]
      call grid%assemble(2, unit_grid, unit_weight, 0.5_wp, err, split, heaviest=.true.)
      call check(error,.not. allocated(err) .and. .not. grid%reduced .and. grid%ncomp == 2, &
         & "one component of weight one per sphere was reduced")
      if (allocated(error)) return

      ! Two components on two spheres, the third sphere empty: as many
      ! components as spheres, yet no grid for the third
      split%weight = [0.4_wp, 0.6_wp]
      call grid%assemble(4, unit_grid, unit_weight, 0.5_wp, err, split, heaviest=.true.)
      call check(error,.not. allocated(err) .and. grid%reduced .and. grid%ncomp == 4 .and. &
         & all(grid%comp_model == [2, 4, 6, 8]), "heaviest components not kept")
   end subroutine test_heaviest_contract

   !> iSwiG with weights and a rotation that follow the geometry
   subroutine test_fade_iswig(error)
      type(error_type), allocatable, intent(out) :: error

      call check_fade_gradient(kind_iswig, "01", 1.0e-5_wp, 1.0e-6_wp, error)
   end subroutine test_fade_iswig

   !> NumSA with weights and a rotation that follow the geometry
   subroutine test_fade_numsa(error)
      type(error_type), allocatable, intent(out) :: error

      call check_fade_gradient(kind_numsa, "01", 1.0e-5_wp, 1.0e-6_wp, error)
   end subroutine test_fade_numsa

   !> DROP on the reduced grid
   !>
   !> - the rotation of the kept component reaches the gradient, its weight
   !>   does not
   subroutine test_fade_drop(error)
      type(error_type), allocatable, intent(out) :: error

      call check_fade_gradient(kind_drop, "01", 1.0e-5_wp, 1.0e-6_wp, error)
   end subroutine test_fade_drop

   !> DROP with one weighted component per sphere
   !>
   !> - the weight is dropped, the gradients are those of the same grid with
   !>   weight one
   subroutine test_fade_drop_single(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_context_type), target :: ctx
      class(cavity_type), allocatable :: single, double
      type(mctc_error), allocatable :: err

      call get_structure(mol, "MB16-43", "01")
      call new_context(ctx, nthreads=0, verbosity=0)
      call fade_cavity(kind_drop, mol, ctx, single, err, single=.true.)
      if (.not. allocated(err)) call single%get_gradient(err)
      if (.not. allocated(err)) call fade_cavity(kind_drop, mol, ctx, double, err)
      if (.not. allocated(err)) call double%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, "DROP failed: "//err%message)
         return
      end if
      call check(error, single%grid%reduced .and. all(single%grid%comp_weight == 1.0_wp), &
         & "a weighted component per sphere kept its weight")
      if (allocated(error)) return
      ! The heavier component of the two-component mock is the same grid
      call check(error, single%ngrid == double%ngrid .and. &
         & abs(single%total_area - double%total_area) <= 1.0e-12_wp*double%total_area, &
         & "single weighted component differs from the heaviest of two")
      if (allocated(error)) return
      select type (single)
      type is (cavity_type_drop)
         select type (double)
         type is (cavity_type_drop)
            call check(error, maxval(abs(single%A_tot1_rA - double%A_tot1_rA)) <= 1.0e-10_wp, &
               & "gradient of a single weighted component differs from the heaviest of two")
         end select
      end select
   end subroutine test_fade_drop_single

end module test_cavity_grid
