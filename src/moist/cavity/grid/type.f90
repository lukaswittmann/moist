!> Angular grid of every atomic sphere as a weighted sum of components
!>
!> Flat container a Lebedev-based cavity loops over, rebuilt once per update
!> from the cavity's own unit grid and, optionally, a grid model
!> ([[moist_cavity_grid_model_type]])
!>
!> - sphere a owns the components `sph_first(a):sph_first(a + 1) - 1`
!> - component c owns the global nodes `comp_first(c):comp_first(c + 1) - 1`;
!>   its unit nodes and angular weights are the node set `comp_set(c)`, stored
!>   once per distinct (size, rotation): node j of the component is
!>   `node(:, set_first(comp_set(c)) + j - 1)` and sits at c_a + r_a node
!> - `comp_weight(c)` scales the switching value of the component's points,
!>   `comp_zeta(c)` is the width parameter of its Lebedev size
!> - `heaviest`: a cavity that cannot carry several components per sphere
!>   keeps the heaviest component of each sphere with weight one, a constant;
!>   `comp_model` maps the kept components to those of the model
!> - without a model, or with a plain one, every sphere has one component of
!>   weight one on the cavity's own nodes (`plain`): bitwise the grid of a
!>   cavity without any grid handling
!>
!> A loop over the nodes of a cavity:
!>
!>     do isph = 1, grid%nsph
!>        do icomp = grid%sph_first(isph), grid%sph_first(isph + 1) - 1
!>           n0 = grid%set_first(grid%comp_set(icomp)) - 1
!>           do j = 1, grid%comp_first(icomp + 1) - grid%comp_first(icomp)
!>              p = center(:, isph) + radius(isph)*grid%node(:, n0 + j)
module moist_cavity_grid_type
   use mctc_env, only: wp, error_type, fatal_error
   use moist_cavity_grid_model_type, only: grid_model_type, grid_components_type, grid_size_default
   implicit none(type, external)
   private

   public :: cavity_grid_type

   !> Angular grid components of every sphere of a cavity
   type :: cavity_grid_type
      !> Number of spheres
      integer :: nsph = 0
      !> Number of components over all spheres
      integer :: ncomp = 0
      !> Number of nodes over all components
      integer :: nnode = 0
      !> Number of distinct node sets
      integer :: nset = 0
      !> Whether every sphere has one unrotated default component of weight one
      logical :: plain = .true.
      !> Component range of each sphere (nsph + 1)
      integer, allocatable :: sph_first(:)
      !> Global node range of each component (ncomp + 1)
      integer, allocatable :: comp_first(:)
      !> Sphere of each component (ncomp)
      integer, allocatable :: comp_sphere(:)
      !> Smooth weight of each component, sums to one per sphere (ncomp)
      real(wp), allocatable :: comp_weight(:)
      !> Width parameter of each component, fitted per Lebedev size (ncomp)
      real(wp), allocatable :: comp_zeta(:)
      !> Distinct rotation of each component, for the frame Jacobian (ncomp)
      integer, allocatable :: comp_rot(:)
      !> Node set of each component (ncomp)
      integer, allocatable :: comp_set(:)
      !> Component of the grid model behind each component (ncomp)
      integer, allocatable :: comp_model(:)
      !> Number of components of the grid model
      integer :: ncomp_model = 0
      !> Whether only the heaviest component of each sphere was kept
      logical :: reduced = .false.
      !> Node range of each node set (nset + 1)
      integer, allocatable :: set_first(:)
      !> Rotated unit nodes of every node set (3, sum of the set sizes)
      real(wp), allocatable :: node(:, :)
      !> Plain angular (Lebedev) weights of every node set
      real(wp), allocatable :: wang(:)
   contains
      !> Build the container from the cavity's unit grid and a grid model
      procedure :: assemble => assemble_grid
      !> Forget the components
      procedure :: clear => clear_grid
   end type cavity_grid_type

contains

   !> Build the components of every sphere
   !>
   !> - no model, or a plain one: one component per sphere on `unit_grid`
   !>   itself (copied, so c + r u is bitwise unchanged)
   !> - otherwise the components of the model; a node set per distinct
   !>   rotation holds R u for every unit node u
   !> - only the cavity's own Lebedev size is available here: a model asking
   !>   for another size is an error until the cavities supply further sizes
   !>
   !> @param[in,out] self        Container
   !> @param[in]     nsph        Number of spheres
   !> @param[in]     unit_grid   Unit nodes of the cavity's Lebedev size (3, npt)
   !> @param[in]     unit_weight Angular weights of the nodes (npt)
   !> @param[in]     zeta        Width parameter of the cavity's Lebedev size
   !> @param[out]    error       Inconsistent or unsupported components
   !> @param[in]     model       Updated grid model
   !> @param[in]     heaviest    Keep the heaviest component of each sphere only
   subroutine assemble_grid(self, nsph, unit_grid, unit_weight, zeta, error, model, heaviest)
      !> Container
      class(cavity_grid_type), intent(inout) :: self
      !> Number of spheres
      integer, intent(in) :: nsph
      !> Unit nodes
      real(wp), intent(in) :: unit_grid(:, :)
      !> Angular weights
      real(wp), intent(in) :: unit_weight(:)
      !> Width parameter
      real(wp), intent(in) :: zeta
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Grid model
      class(grid_model_type), intent(in), optional :: model
      !> Keep the heaviest component only
      logical, intent(in), optional :: heaviest

      !> Components of the model
      type(grid_components_type) :: comps
      logical :: plain
      integer :: npt, isph, icomp, iset, ncomp, nset, jcomp
      integer, allocatable :: keep(:)

      npt = size(unit_weight)
      plain = .true.
      if (present(model)) plain = model%is_plain()

      if (plain) then
         self%nsph = nsph
         self%ncomp = nsph
         self%nnode = nsph*npt
         self%nset = 1
         self%plain = .true.
         self%sph_first = [(isph, isph=1, nsph + 1)]
         self%comp_first = [((icomp - 1)*npt + 1, icomp=1, nsph + 1)]
         self%comp_sphere = [(isph, isph=1, nsph)]
         self%comp_weight = spread(1.0_wp, 1, nsph)
         self%comp_zeta = spread(zeta, 1, nsph)
         self%comp_rot = spread(1, 1, nsph)
         self%comp_set = spread(1, 1, nsph)
         self%comp_model = [(isph, isph=1, nsph)]
         self%ncomp_model = nsph
         self%reduced = .false.
         self%set_first = [1, npt + 1]
         self%node = unit_grid
         self%wang = unit_weight
         return
      end if

      call model%components(nsph, comps)
      if (.not. (allocated(comps%sph_first) .and. allocated(comps%size) .and. allocated(comps%rot) &
          & .and. allocated(comps%weight) .and. allocated(comps%rotation))) then
         call fatal_error(error, "Grid model returned no components")
         return
      end if
      ncomp = size(comps%rot)
      nset = size(comps%rotation, 3)
      if (size(comps%sph_first) /= nsph + 1 .or. size(comps%size) /= ncomp .or. size(comps%weight) /= ncomp) then
         call fatal_error(error, "Grid model returned inconsistent components")
         return
      end if
      if (comps%sph_first(1) /= 1 .or. comps%sph_first(nsph + 1) /= ncomp + 1 .or. &
          & any(comps%sph_first(2:) < comps%sph_first(:nsph)) .or. any(comps%rot < 1) .or. any(comps%rot > nset)) then
         call fatal_error(error, "Grid model returned inconsistent components")
         return
      end if
      if (any(comps%size /= grid_size_default .and. comps%size /= npt)) then
         call fatal_error(error, "Grid model asks for a Lebedev size other than the cavity's, "// &
            & "which is not supported yet")
         return
      end if

      self%ncomp_model = ncomp
      self%comp_model = [(icomp, icomp=1, ncomp)]
      self%reduced = .false.
      if (present(heaviest)) then
         if (heaviest) then
            ! Every sphere needs a component to keep
            if (any(comps%sph_first(2:) <= comps%sph_first(:nsph))) then
               call fatal_error(error, "Grid model returned a sphere without components")
               return
            end if
            ! One component per sphere with weight one needs no reduction
            if (ncomp /= nsph .or. any(comps%weight /= 1.0_wp)) then
               ! The heaviest component of every sphere, the first among equals
               allocate (keep(nsph))
               do isph = 1, nsph
                  jcomp = comps%sph_first(isph)
                  do icomp = comps%sph_first(isph) + 1, comps%sph_first(isph + 1) - 1
                     if (comps%weight(icomp) > comps%weight(jcomp)) jcomp = icomp
                  end do
                  keep(isph) = jcomp
               end do
               comps%size = comps%size(keep)
               comps%rot = comps%rot(keep)
               comps%weight = spread(1.0_wp, 1, nsph)
               comps%sph_first = [(isph, isph=1, nsph + 1)]
               self%comp_model = keep
               self%reduced = .true.
               ncomp = nsph
            end if
         end if
      end if

      self%nsph = nsph
      self%ncomp = ncomp
      self%nnode = ncomp*npt
      self%nset = nset
      self%plain = .false.
      self%sph_first = comps%sph_first
      self%comp_first = [((icomp - 1)*npt + 1, icomp=1, ncomp + 1)]
      if (allocated(self%comp_sphere)) deallocate (self%comp_sphere)
      allocate (self%comp_sphere(ncomp))
      do isph = 1, nsph
         self%comp_sphere(comps%sph_first(isph):comps%sph_first(isph + 1) - 1) = isph
      end do
      self%comp_weight = comps%weight
      self%comp_zeta = spread(zeta, 1, ncomp)
      self%comp_rot = comps%rot
      ! One node set per distinct rotation while every component has the cavity's size
      self%comp_set = comps%rot
      self%set_first = [((iset - 1)*npt + 1, iset=1, nset + 1)]
      if (allocated(self%node)) then
         if (any(shape(self%node) /= [3, nset*npt])) deallocate (self%node)
      end if
      if (.not. allocated(self%node)) allocate (self%node(3, nset*npt))
      if (allocated(self%wang)) then
         if (size(self%wang) /= nset*npt) deallocate (self%wang)
      end if
      if (.not. allocated(self%wang)) allocate (self%wang(nset*npt))
      do iset = 1, nset
         self%node(:, (iset - 1)*npt + 1:iset*npt) = matmul(comps%rotation(:, :, iset), unit_grid)
         self%wang((iset - 1)*npt + 1:iset*npt) = unit_weight
      end do

   end subroutine assemble_grid

   !> Forget the components
   !>
   !> @param[in,out] self Container
   subroutine clear_grid(self)
      !> Container
      class(cavity_grid_type), intent(inout) :: self

      self%nsph = 0
      self%ncomp = 0
      self%nnode = 0
      self%nset = 0
      self%plain = .true.
      if (allocated(self%sph_first)) deallocate (self%sph_first)
      if (allocated(self%comp_first)) deallocate (self%comp_first)
      if (allocated(self%comp_sphere)) deallocate (self%comp_sphere)
      if (allocated(self%comp_weight)) deallocate (self%comp_weight)
      if (allocated(self%comp_zeta)) deallocate (self%comp_zeta)
      if (allocated(self%comp_rot)) deallocate (self%comp_rot)
      if (allocated(self%comp_set)) deallocate (self%comp_set)
      if (allocated(self%comp_model)) deallocate (self%comp_model)
      self%ncomp_model = 0
      self%reduced = .false.
      if (allocated(self%set_first)) deallocate (self%set_first)
      if (allocated(self%node)) deallocate (self%node)
      if (allocated(self%wang)) deallocate (self%wang)
   end subroutine clear_grid

end module moist_cavity_grid_type
