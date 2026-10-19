!> Frozen Gaussian surface quadrature on a retained marching-tetrahedra mesh
module moist_cavity_mt
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_structure, only: structure_type
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_cavity_type, only: cavity_type, update_cavity_base, build_cavity_base, &
      & cavity_declare_pass, list_cavity_fields_base
   use moist_cavity_mt_parameters, only: moist_cavity_mt_parameters_type
   use moist_cavity_mt_mesh, only: mt_mesh_type, build_mt_mesh
   use moist_cavity_mt_patches, only: mt_patch_type, build_mt_patches
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_switching, only: moist_cavity_drop_swif_sigmoid_bump_type, new_swif_sigmoid_bump
   use moist_context, only: moist_context_type
   use moist_radius_type, only: radius_type
   use moist_channels_coupling, only: coupling_type, coupling_view_type, coupling_make_view, coupling_close_view
   use moist_channels_fields, only: field_query_type
   use moist_utils_timer, only: cat_properties
!$ use omp_lib, only: omp_get_thread_num
   implicit none(type, external)
   private
   public :: cavity_type_mt, new_cavity_mt
   public :: moist_cavity_mt_parameters_type

   !> Frozen cavity with separate geometric and electrostatic discretizations
   type, extends(cavity_type) :: cavity_type_mt
      !> Copied LSF template, including build-phase inputs when required
      class(moist_cavity_drop_lsf_type), allocatable :: lsf_model
      !> Geometry and quadrature controls
      type(moist_cavity_mt_parameters_type) :: param
      !> Complete indexed planar mesh, including suppressed surface regions
      type(mt_mesh_type) :: mesh
      !> Complete unswitched integration patches
      type(mt_patch_type) :: patches
      !> Raw planar-mesh area, bohr**2
      real(wp), allocatable :: mesh_area
      !> Raw oriented planar-mesh volume, bohr**3
      real(wp), allocatable :: mesh_volume
      !> Geometric area of each retained grid point before switching, bohr**2
      real(wp), allocatable :: patch_area(:)
      !> Area-averaged critical-gradient switch of each retained patch
      real(wp), allocatable :: critical_weight(:)
      !> One-based index into the complete patch collection for each grid point
      integer, allocatable :: patch_index(:)
   contains
      !> Accept normals; reject unavailable curvature
      procedure :: configure => configure_mt
      !> Update geometry, deferring density-dependent LSFs until build
      procedure :: update => update_mt
      !> Whether the LSF requires build-phase host data
      procedure :: defers_build => mt_defers_build
      !> Consume host data and construct a frozen surface
      procedure :: build => build_mt
      !> Declare LSF host requests
      procedure :: declare_coupling => declare_mt_coupling
      !> Release both the grid and the retained geometry
      procedure :: clear_surface => clear_mt_surface
      !> Reject nuclear derivatives
      procedure :: get_gradient => get_mt_gradient
      !> Expose raw geometry alongside effective quadrature
      procedure :: list_fields => list_mt_fields
      !> Construction settings
      procedure :: parameters => mt_parameters
   end type cavity_type_mt

   !> Private gradient evaluator and first error for one OpenMP thread
   type :: mt_gradient_slot
      !> Spatial-derivative LSF copy
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      !> First evaluation error on this thread
      type(error_type), allocatable :: error
   end type mt_gradient_slot

contains

   !> Construct a cavity by copying its radius model, LSF, and settings
   !>
   !> @param[in,out] self Cavity to initialize
   !> @param[in] radius_model Atomic radius model to copy
   !> @param[in] lsf_model Level set function to copy
   !> @param[out] error Invalid settings or allocation failure
   !> @param[in] param Construction settings; omitted uses compiled defaults
   !> @param[in] ctx Borrowed run context; a model may supply its own
   subroutine new_cavity_mt(self, radius_model, lsf_model, error, param, ctx)
      !> Cavity to initialize
      type(cavity_type_mt), intent(inout) :: self
      !> Atomic radius model to copy
      class(radius_type), intent(in) :: radius_model
      !> Level set function to copy
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf_model
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Optional construction settings
      type(moist_cavity_mt_parameters_type), intent(in), optional :: param
      !> Borrowed run context
      type(moist_context_type), intent(in), optional, target :: ctx
      type(moist_cavity_mt_parameters_type) :: settings
      integer :: stat

      if (present(param)) settings = param
      call settings%validate(error)
      if (allocated(error)) return
      call self%clear_surface()
      call self%clear_couplings()
      self%updated = .false.
      if (allocated(self%mol)) deallocate (self%mol)
      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%sphxyz)) deallocate (self%sphxyz)
      self%nsph = 0
      nullify (self%ctx)
      if (present(ctx)) self%ctx => ctx
      self%label = "Marching tetrahedra"
      self%param = settings
      if (allocated(self%radius_model)) deallocate (self%radius_model)
      allocate (self%radius_model, source=radius_model, stat=stat)
      if (stat == 0) then
         if (allocated(self%lsf_model)) deallocate (self%lsf_model)
         allocate (self%lsf_model, source=lsf_model, stat=stat)
      end if
      if (stat /= 0) call fatal_error(error, "Cannot allocate marching tetrahedra models")
   end subroutine new_cavity_mt

   !> Accept the always-present LSF normals and reject curvatures
   !>
   !> @param[in,out] self Cavity instance
   !> @param[in] normal Optional normal request
   !> @param[in] curvature Optional curvature request
   !> @param[out] error Unsupported property
   subroutine configure_mt(self, normal, curvature, error)
      !> Cavity instance
      class(cavity_type_mt), intent(inout) :: self
      !> Optional normal request
      logical, intent(in), optional :: normal
      !> Optional curvature request
      logical, intent(in), optional :: curvature
      !> Unsupported property
      type(error_type), allocatable, intent(out) :: error

      if (present(curvature)) then
         if (curvature) call fatal_error(error, "Marching tetrahedra does not provide surface curvatures")
      end if
   end subroutine configure_mt

   !> Update sphere geometry and build geometric level sets immediately
   !>
   !> @param[in,out] self Cavity instance
   !> @param[in] mol Molecular structure
   !> @param[out] error Construction failure
   subroutine update_mt(self, mol, error)
      !> Cavity instance
      class(cavity_type_mt), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      integer :: stat

      call update_cavity_base(self, mol, error)
      if (allocated(error)) return
      call self%clear_surface()
      call self%require_context(error)
      if (allocated(error)) return
      call self%param%validate(error)
      if (allocated(error)) return
      if (.not. allocated(self%radius_model) .or. .not. allocated(self%lsf_model)) then
         call fatal_error(error, "Construct the marching tetrahedra cavity before updating it")
         return
      end if
      if (mol%nat < 1) then
         call fatal_error(error, "Marching tetrahedra requires at least one atom")
         return
      end if
      self%nsph = mol%nat
      call self%radius_model%update(mol, error)
      if (allocated(error)) return
      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%sphxyz)) deallocate (self%sphxyz)
      allocate (self%radii(self%nsph), stat=stat)
      if (stat == 0) allocate (self%sphxyz(3, self%nsph), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate marching tetrahedra atom data")
         return
      end if
      self%radii = self%radius_model%f0
      self%sphxyz = mol%xyz
      if (.not. self%defers_build()) then
         call construct_mt_surface(self, error)
         if (allocated(error)) return
      end if
      self%updated = .true.
      self%built = .not. self%defers_build()
   end subroutine update_mt

   !> Whether construction waits for host data
   !>
   !> @param[in] self Cavity instance
   function mt_defers_build(self) result(defers)
      !> Cavity instance
      class(cavity_type_mt), intent(in) :: self
      !> Build-phase dependency
      logical :: defers

      defers = .false.
      if (allocated(self%lsf_model)) defers = self%lsf_model%requires_build_inputs()
   end function mt_defers_build

   !> Build a frozen surface after reading the LSF's host inputs
   !>
   !> @param[in,out] self Updated cavity
   !> @param[in,out] coupling Staged build coupling
   !> @param[out] error Build failure
   subroutine build_mt(self, coupling, error)
      !> Updated cavity
      class(cavity_type_mt), intent(inout) :: self
      !> Staged build coupling
      type(coupling_type), intent(inout), target :: coupling
      !> Build failure
      type(error_type), allocatable, intent(out) :: error
      type(coupling_view_type) :: view

      call build_cavity_base(self, coupling, error)
      if (allocated(error)) return
      if (.not. self%defers_build()) return
      call self%clear_surface()
      call self%param%validate(error)
      if (allocated(error)) return
      call coupling_make_view(coupling, 0, view)
      call self%lsf_model%read_build(view, self%mol, error)
      call coupling_close_view(coupling)
      if (allocated(error)) return
      call construct_mt_surface(self, error)
      if (allocated(error)) return
      self%built = .true.
      call cavity_declare_pass(self, coupling, error)
   end subroutine build_mt

   !> Declare build-phase LSF requests without density-induced geometry response
   !>
   !> @param[in] self Cavity instance
   !> @param[in,out] coupling Host coupling
   !> @param[out] error Declaration failure
   subroutine declare_mt_coupling(self, coupling, error)
      !> Cavity instance
      class(cavity_type_mt), intent(in) :: self
      !> Host coupling
      type(coupling_type), intent(inout) :: coupling
      !> Declaration failure
      type(error_type), allocatable, intent(out) :: error

      if (allocated(self%lsf_model)) call self%lsf_model%declare_coupling(coupling, error)
   end subroutine declare_mt_coupling

   !> Extract geometry, partition it, then evaluate the LSF at patch points
   !>
   !> @param[in,out] self Updated cavity geometry
   !> @param[out] error Construction failure
   subroutine construct_mt_surface(self, error)
      !> Updated cavity geometry
      class(cavity_type_mt), intent(inout) :: self
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error

      ! Screened exterior points can return a zero sentinel instead of the LSF
      ! Full-scan evaluation preserves the zero set and its spatial normals
      call self%lsf_model%set_screening_threshold(0.0_wp)
      call self%lsf_model%update(self%mol, self%radii)
      call self%ctx%timer%start("Marching tetrahedra", category=cat_properties)
      call build_mt_mesh(self%lsf_model, self%sphxyz, self%param%spacing, self%mesh, error, &
         & nthreads=self%ctx%get_num_threads())
      if (.not. allocated(error)) then
         call build_mt_patches(self%mesh, self%param%target_area, self%param%area_tolerance, self%patches, error)
      end if
      if (.not. allocated(error) .and. .not. self%patches%balanced) call warn_unbalanced(self)
      if (.not. allocated(error)) call finalize_mt_grid(self, error)
      call self%ctx%timer%stop("Marching tetrahedra")
      if (allocated(error)) call self%clear_surface()
   end subroutine construct_mt_surface

   !> Integrate critical switching over patches and populate the common cavity fields
   !>
   !> - The switch is sampled at every mesh face centroid, so seams converge
   !>   with the mesh spacing rather than the patch size
   !> - Normals come from the LSF gradient at each patch point, or from the
   !>   switched mesh normal where that gradient vanishes
   !>
   !> @param[in,out] self Cavity with geometric patches
   !> @param[out] error Evaluation or allocation failure
   subroutine finalize_mt_grid(self, error)
      !> Cavity with geometric patches
      class(cavity_type_mt), intent(inout) :: self
      !> Evaluation or allocation failure
      type(error_type), allocatable, intent(out) :: error
      type(mt_gradient_slot), allocatable :: slots(:)
      type(moist_cavity_drop_swif_sigmoid_bump_type) :: critical
      real(wp), allocatable :: gradients(:, :), weights(:), effective(:), face_weight(:)
      logical, allocatable :: active(:)
      real(wp) :: grad(3), magnitude, d2, best, origin(3), mean_normal(3)
      integer :: nf, np, nt, stat, i, j, k, atom, tid

      nf = size(self%mesh%area)
      np = size(self%patches%area)
      nt = min(self%ctx%get_num_threads(), max(1, nf))
      allocate (slots(nt), stat=stat)
      if (stat == 0) allocate (face_weight(nf), stat=stat)
      if (stat == 0) allocate (gradients(3, np), stat=stat)
      if (stat == 0) allocate (weights(np), stat=stat)
      if (stat == 0) allocate (effective(np), stat=stat)
      if (stat == 0) allocate (active(np), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate marching tetrahedra gradient workspace")
         return
      end if
      do i = 1, nt
         allocate (slots(i)%lsf, source=self%lsf_model, stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot clone the marching tetrahedra level set")
            return
         end if
         call slots(i)%lsf%set_max_deriv(1)
      end do
      call new_swif_sigmoid_bump(critical, self%param%w_0ls_from, self%param%w_0ls_to, &
         & p_hi=self%param%w_0ls_p, a_hi=self%param%w_0ls_a, &
         & p_lo=self%param%w_0ls_p, a_lo=self%param%w_0ls_a)
      face_weight = 0.0_wp
!$omp parallel do default(none) schedule(static) num_threads(nt) &
!$omp shared(self, slots, face_weight, critical, nf) private(i, tid, grad, magnitude)
      do i = 1, nf
         tid = 1
!$       tid = omp_get_thread_num() + 1
         if (allocated(slots(tid)%error)) cycle
         call point_gradient(slots(tid), self%mesh%centroid(:, i), grad, magnitude)
         if (allocated(slots(tid)%error)) cycle
         face_weight(i) = critical%f0(magnitude)
         if (.not. ieee_is_finite(face_weight(i))) then
            call fatal_error(slots(tid)%error, "Nonfinite MT critical-gradient switch")
         end if
      end do
!$omp end parallel do
      call collect_slot_error(slots, error)
      if (allocated(error)) return
      weights = 0.0_wp
      effective = 0.0_wp
      gradients = 0.0_wp
!$omp parallel do default(none) schedule(dynamic, 16) num_threads(nt) &
!$omp shared(self, slots, face_weight, gradients, weights, effective, np) &
!$omp private(i, j, k, tid, grad, magnitude, mean_normal)
      do i = 1, np
         tid = 1
!$       tid = omp_get_thread_num() + 1
         if (allocated(slots(tid)%error)) cycle
         mean_normal = 0.0_wp
         do j = self%patches%member_offset(i), self%patches%member_offset(i + 1) - 1
            k = self%patches%member_face(j)
            effective(i) = effective(i) + face_weight(k)*self%patches%member_area(j)
            mean_normal = mean_normal + face_weight(k)*self%patches%member_area(j)*self%mesh%normal(:, k)
         end do
         weights(i) = effective(i)/self%patches%area(i)
         if (weights(i) <= self%param%weight_cutoff) then
            weights(i) = 0.0_wp
            cycle
         end if
         call point_gradient(slots(tid), self%patches%xyz(:, i), grad, magnitude)
         if (allocated(slots(tid)%error)) cycle
         if (magnitude > tiny(1.0_wp)) then
            gradients(:, i) = grad/magnitude
            cycle
         end if
         magnitude = norm2(mean_normal)
         if (magnitude > tiny(1.0_wp)) then
            gradients(:, i) = mean_normal/magnitude
         else
            weights(i) = 0.0_wp
         end if
      end do
!$omp end parallel do
      call collect_slot_error(slots, error)
      if (allocated(error)) return
      active = weights > 0.0_wp
      call self%lsf_model%check_surface(self%patches%xyz, active, error)
      if (allocated(error)) return
      self%ngrid = count(active)
      allocate (self%xyz(3, self%ngrid), stat=stat)
      if (stat == 0) allocate (self%a(self%ngrid), stat=stat)
      if (stat == 0) allocate (self%f(self%ngrid), stat=stat)
      if (stat == 0) allocate (self%xi0(self%ngrid), stat=stat)
      if (stat == 0) allocate (self%normal0(3, self%ngrid), stat=stat)
      if (stat == 0) allocate (self%owner(self%ngrid), stat=stat)
      if (stat == 0) allocate (self%asph(self%nsph), stat=stat)
      if (stat == 0) allocate (self%v(self%ngrid), stat=stat)
      if (stat == 0) allocate (self%patch_area(self%ngrid), stat=stat)
      if (stat == 0) allocate (self%critical_weight(self%ngrid), stat=stat)
      if (stat == 0) allocate (self%patch_index(self%ngrid), stat=stat)
      if (stat == 0) allocate (self%total_area, stat=stat)
      if (stat == 0) allocate (self%total_volume, stat=stat)
      if (stat == 0) allocate (self%mesh_area, stat=stat)
      if (stat == 0) allocate (self%mesh_volume, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate the marching tetrahedra integration grid")
         return
      end if
      self%asph = 0.0_wp
      self%f = 1.0_wp
      origin = self%get_volume_origin()
      k = 0
      do i = 1, np
         if (.not. active(i)) cycle
         k = k + 1
         self%xyz(:, k) = self%patches%xyz(:, i)
         self%normal0(:, k) = gradients(:, i)
         self%patch_index(k) = i
         self%patch_area(k) = self%patches%area(i)
         self%critical_weight(k) = weights(i)
         self%a(k) = effective(i)
         if (self%a(k) <= tiny(1.0_wp)) then
            call fatal_error(error, "Marching tetrahedra effective patch area underflow")
            return
         end if
         self%xi0(k) = self%param%gaussian_scale/sqrt(self%a(k))
         self%v(k) = self%a(k)*dot_product(self%xyz(:, k) - origin, self%normal0(:, k))/3.0_wp
         best = huge(1.0_wp)
         atom = 1
         do j = 1, self%nsph
            d2 = sum((self%xyz(:, k) - self%sphxyz(:, j))**2)
            if (d2 >= best) cycle
            best = d2
            atom = j
         end do
         self%owner(k) = atom
         self%asph(atom) = self%asph(atom) + self%a(k)
      end do
      if (.not. all(ieee_is_finite(self%xi0))) then
         call fatal_error(error, "Nonfinite marching tetrahedra Gaussian widths")
         return
      end if
      if (any(self%xi0 <= 0.0_wp)) then
         call fatal_error(error, "Marching tetrahedra Gaussian widths underflow")
         return
      end if
      self%total_area = sum(self%a)
      self%total_volume = sum(self%v)
      self%mesh_area = self%mesh%total_area
      self%mesh_volume = self%mesh%total_volume
   end subroutine finalize_mt_grid

   !> Finite spatial LSF gradient at one point; failures latch on the slot
   !>
   !> @param[in,out] slot Thread-local evaluator
   !> @param[in] xyz Evaluation point, bohr
   !> @param[out] grad Spatial gradient
   !> @param[out] magnitude Gradient norm
   subroutine point_gradient(slot, xyz, grad, magnitude)
      !> Thread-local evaluator
      type(mt_gradient_slot), intent(inout) :: slot
      !> Evaluation point
      real(wp), intent(in) :: xyz(3)
      !> Spatial gradient
      real(wp), intent(out) :: grad(3)
      !> Gradient norm
      real(wp), intent(out) :: magnitude

      grad = 0.0_wp
      magnitude = 0.0_wp
      call slot%lsf%prepare(xyz, slot%error)
      if (allocated(slot%error)) return
      call slot%lsf%f012_r(lsf1_r=grad)
      if (.not. all(ieee_is_finite(grad))) then
         call fatal_error(slot%error, "Nonfinite marching tetrahedra LSF gradient")
         return
      end if
      magnitude = norm2(grad)
      if (.not. ieee_is_finite(magnitude)) then
         call fatal_error(slot%error, "Nonfinite MT LSF gradient magnitude")
      end if
   end subroutine point_gradient

   !> Move the first latched thread error to the caller
   !>
   !> @param[in,out] slots Thread-local evaluators
   !> @param[out] error First failure
   subroutine collect_slot_error(slots, error)
      !> Thread-local evaluators
      type(mt_gradient_slot), intent(inout) :: slots(:)
      !> First failure
      type(error_type), allocatable, intent(out) :: error
      integer :: i

      do i = 1, size(slots)
         if (.not. allocated(slots(i)%error)) cycle
         call move_alloc(slots(i)%error, error)
         return
      end do
   end subroutine collect_slot_error

   !> Report patches kept outside the requested area tolerance
   !>
   !> @param[in] self Cavity with fallback patches
   subroutine warn_unbalanced(self)
      !> Cavity with fallback patches
      class(cavity_type_mt), intent(in) :: self
      character(len=32) :: deviation, tolerance

      write (deviation, "(f12.1)") 100.0_wp*self%patches%max_deviation
      write (tolerance, "(f12.1)") 100.0_wp*self%param%area_tolerance
      call self%ctx%warning("Marching tetrahedra patches missed the area tolerance ("// &
         & trim(adjustl(deviation))//"% > "//trim(adjustl(tolerance))//"%); using the best connected partition found")
   end subroutine warn_unbalanced

   !> Release the grid and its retained construction geometry
   !>
   !> @param[in,out] self Cavity instance
   subroutine clear_mt_surface(self)
      !> Cavity instance
      class(cavity_type_mt), intent(inout) :: self

      call self%clear_surface()
      call self%mesh%clear()
      call self%patches%clear()
      if (allocated(self%mesh_area)) deallocate (self%mesh_area)
      if (allocated(self%mesh_volume)) deallocate (self%mesh_volume)
      if (allocated(self%patch_area)) deallocate (self%patch_area)
      if (allocated(self%critical_weight)) deallocate (self%critical_weight)
      if (allocated(self%patch_index)) deallocate (self%patch_index)
   end subroutine clear_mt_surface

   !> Reject nuclear derivatives of the frozen discretization
   !>
   !> @param[in,out] self Cavity instance
   !> @param[out] error Unsupported derivative
   subroutine get_mt_gradient(self, error)
      !> Cavity instance
      class(cavity_type_mt), intent(inout) :: self
      !> Unsupported derivative
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "Marching tetrahedra is frozen and provides no nuclear geometry derivatives")
   end subroutine get_mt_gradient

   !> Expose common quadrature and independent mesh diagnostics
   !>
   !> @param[in] self Cavity instance
   !> @param[in,out] query Named-field collector
   subroutine list_mt_fields(self, query)
      !> Cavity instance
      class(cavity_type_mt), intent(in) :: self
      !> Named-field collector
      type(field_query_type), intent(inout) :: query

      call list_cavity_fields_base(self, query)
      if (.not. self%built) return
      call query%add_real_scalar("mesh_area", "Unswitched planar-mesh area, bohr**2", self%mesh_area)
      call query%add_real_scalar("mesh_volume", "Oriented planar-mesh volume, bohr**3", self%mesh_volume)
      call query%add_real("patch_area", "Geometric point area before switching, bohr**2 (ngrid)", self%patch_area)
      call query%add_real("critical_weight", "Area-averaged critical-gradient switch (ngrid)", &
         & self%critical_weight)
      call query%add_int("patch_index", "Index in unswitched patches, 0-based (ngrid)", &
         & self%patch_index, zero_based=.true.)
      call query%add_real("all_patch_area", "All geometric patch areas before suppression, bohr**2", self%patches%area)
      call query%add_real("component_target", "Attainable area per patch within its component, bohr**2", &
         & self%patches%component_target)
      call query%add_bool("area_exception", "Component rounding prevents the requested area tolerance", &
         & self%patches%area_exception)
      call query%add_int("patch_component", "Connected component of every unswitched patch, 0-based", &
         & self%patches%component, zero_based=.true.)
   end subroutine list_mt_fields

   !> Return construction settings borrowed from the cavity
   !>
   !> @param[in] self Cavity instance
   function mt_parameters(self) result(param)
      !> Cavity instance
      class(cavity_type_mt), intent(in), target :: self
      !> Borrowed parameter reference
      class(moist_model_parameters_type), pointer :: param

      param => self%param
   end function mt_parameters

end module moist_cavity_mt
