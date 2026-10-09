!> Host potential with constrained ESP-fitted Ng tail charges
!>
!> - Own point-potential request: phi in every phase, dphi_dr for gradients
!> - Exact total charge, quadratic reference restraint and linear excess-charge penalty
!> - Smooth nuclear exclusion and extra weight on distant samples
!> - Fit built once per geometry in `update_grid` from the model's volume grid
!> - Matrix-free: mask, weights and Hessian kept; 1/r basis recomputed in grid blocks
!> - Charges fitted to each current host potential; any other grid refused by name
!> - Fit adjoints returned through the potential, nuclei and volume grid
!> - Volume grids only; no host charge request
module moist_model_moz_potential_electrostatic_ec
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use moist_channels_coupling, only: coupling_type, coupling_view_type, point_potential_request_type, &
      & request_require, coupling_register, moist_phase_energy, moist_phase_response, moist_phase_gradient
   use moist_channels_response, only: response_type, response_accumulate, potential_adjoint_response_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_lapack, only: potrf, potrs
   use moist_model_moz_potential_term, only: potential_term_type, potential_name_len, charge_label_len
   !$ use omp_lib, only: omp_get_thread_num
   implicit none(type, external)
   private
   public :: ec_charges_type

   !> Phases requiring the host potential
   integer, parameter :: all_phases(3) = [moist_phase_energy, moist_phase_response, moist_phase_gradient]

   !> Grid points per block of the basis products
   integer, parameter :: block_size = 256

   !> Weighted ESP fit on the volume grid of the update
   !>
   !> - No stored basis: A_ga = 1/|r_g - R_a| recomputed blockwise where the mask is positive
   type :: ec_fit_type
      !> Grid coordinates of the update, bohr (3, ngrid)
      real(wp), allocatable :: xyz(:, :)
      !> Volume weights of the update, bohr**3 (ngrid)
      real(wp), allocatable :: volume(:)
      !> Exclusion mask times long-range weight (ngrid)
      real(wp), allocatable :: mask(:)
      !> Normalized fitting weights (ngrid)
      real(wp), allocatable :: weight(:)
      !> Restrained Hessian, 1/bohr**2 (natom, natom)
      real(wp), allocatable :: h(:, :)
      !> Sum of unnormalized fitting weights, bohr**3
      real(wp) :: normalization = 0.0_wp
      !> Threads of the grid's team for the basis products
      integer :: nthreads = 1
   end type ec_fit_type

   !> Potential term with fitted, restrained atom-centered tail charges
   !>
   !> Configure before `add` or `update`; references absent by default
   type, extends(potential_term_type) :: ec_charges_type

      !* ---------------------------- EC LR charge fitting --------------------------- *!

      !> Nuclear exclusion radius; full weight beyond twice this radius, bohr
      real(wp) :: exclusion_radius = 1.0_wp
      !> Distance scale of the far-field weighting, bohr
      real(wp) :: range_scale = 5.0_wp
      !> Far-field weight increases from 1 to 1 + long_range_bias
      real(wp) :: long_range_bias = 1.0_wp
      !> Quadratic reference-charge restraint, 1/bohr**2; positive
      real(wp) :: restraint_strength = 0.01_wp
      !> Slope of the linear excess-charge penalty, e/bohr**2; nonnegative
      real(wp) :: charge_penalty = 0.1_wp

      !* ------------------------ Reference charges and bounds ------------------------ *!

      !> Optional reference charges, e (natom); default total charge / natom
      real(wp), allocatable :: reference_charges(:)
      !> Optional thresholds for the excess-charge penalty, e (natom)
      real(wp), allocatable :: charge_thresholds(:)
      !> Optional lower charge bounds, e (natom); unbounded by default
      real(wp), allocatable :: lower_bounds(:)
      !> Optional upper charge bounds, e (natom); unbounded by default
      real(wp), allocatable :: upper_bounds(:)

      !* ----------------------------- Updated solute data ---------------------------- *!

      !> Updated nuclear coordinates, bohr (3, natom)
      real(wp), allocatable, private :: atom_xyz(:, :)
      !> Updated reference charges, e (natom)
      real(wp), allocatable, private :: reference(:)
      !> Updated excess-charge penalty thresholds, e (natom)
      real(wp), allocatable, private :: threshold(:)
      !> Updated lower charge bounds, e (natom)
      real(wp), allocatable, private :: lower(:)
      !> Updated upper charge bounds, e (natom)
      real(wp), allocatable, private :: upper(:)
      !> Structure's total charge, e
      real(wp), private :: total_charge = 0.0_wp

      !* ------------------------------ Fit of the update ----------------------------- *!

      !> Fit on the volume grid, absent until a successful `update_grid`
      type(ec_fit_type), allocatable, private :: fit
   contains

      !  Term lifecycle
      procedure :: name => ec_charges_name
      procedure :: update => ec_charges_update
      procedure :: update_grid_1d => ec_charges_update_grid_1d
      procedure :: update_grid_3d => ec_charges_update_grid_3d

      ! Field and charge coverage
      procedure :: is_host_fed => ec_charges_is_host_fed
      procedure :: has_field => ec_charges_has_field
      procedure :: charges_covered => ec_charges_covered
      procedure :: charge_label => ec_charges_label

      ! Host requirements
      procedure :: declare_pass_1d => ec_charges_declare_1d
      procedure :: declare_pass_3d => ec_charges_declare_3d

      ! Potential evaluation and adjoint
      procedure :: tail_charges => ec_charges_tail
      procedure :: field_3d => ec_charges_field
      procedure :: accumulate_adjoint_3d => ec_charges_adjoint
   end type ec_charges_type

contains

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function ec_charges_name(self) result(name)
      !> Term
      class(ec_charges_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "ec_charges"
   end function ec_charges_name

   !> Host-fed by construction
   !>
   !> @param[in] self  term
   pure function ec_charges_is_host_fed(self) result(host_fed)
      !> Term
      class(ec_charges_type), intent(in) :: self
      !> Host-fed flag
      logical :: host_fed
      host_fed = .true.
   end function ec_charges_is_host_fed

   !> Supplies the full host field
   !>
   !> @param[in] self  term
   pure function ec_charges_has_field(self) result(has)
      !> Term
      class(ec_charges_type), intent(in) :: self
      !> Field presence
      logical :: has
      has = .true.
   end function ec_charges_has_field

   !> Covers every atom after update
   !>
   !> @param[in] self   term
   !> @param[in] natom  number of atoms
   pure function ec_charges_covered(self, natom) result(covered)
      !> Term
      class(ec_charges_type), intent(in) :: self
      !> Number of atoms
      integer, intent(in) :: natom
      !> Charge coverage (natom)
      logical :: covered(natom)
      covered = .false.
      if (allocated(self%atom_xyz)) covered = size(self%atom_xyz, 2) == natom
   end function ec_charges_covered

   !> Charge-source label
   !>
   !> @param[in] self  term
   pure function ec_charges_label(self) result(label)
      !> Term
      class(ec_charges_type), intent(in) :: self
      !> Label
      character(len=charge_label_len) :: label
      label = "EC"
   end function ec_charges_label

   !> Record the structure and validate charge restraints
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         solute structure
   !> @param[out]    error       invalid settings, structure or infeasible bounds
   !> @param[in]     solvent_id  unused; EC serves a solute
   subroutine ec_charges_update(self, mol, error, solvent_id)
      !> Term
      class(ec_charges_type), intent(inout) :: self
      !> Solute structure
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Unused solvent id
      integer, intent(in), optional :: solvent_id
      real(wp) :: settings(5)
      integer :: i, n
      if (allocated(self%fit)) deallocate (self%fit)
      if (allocated(self%atom_xyz)) deallocate (self%atom_xyz)
      if (allocated(self%reference)) deallocate (self%reference, self%threshold, self%lower, self%upper)
      settings = [self%exclusion_radius, self%range_scale, self%long_range_bias, &
         & self%restraint_strength, self%charge_penalty]
      if (.not. all(ieee_is_finite(settings))) then
         call fatal_error(error, "EC fit settings must be finite")
         return
      end if
      if (self%exclusion_radius <= 0.0_wp .or. self%range_scale <= 0.0_wp .or. &
         & self%long_range_bias < 0.0_wp .or. self%restraint_strength <= 0.0_wp .or. self%charge_penalty < 0.0_wp) then
         call fatal_error(error, "EC fit requires positive radii and restraint, nonnegative bias and penalty")
         return
      end if
      n = mol%nat
      if (n < 1 .or. .not. ieee_is_finite(mol%charge)) then
         call fatal_error(error, "EC fit needs atoms and a finite total charge")
         return
      end if
      if (.not. all(ieee_is_finite(mol%xyz))) then
         call fatal_error(error, "EC fit needs finite nuclear coordinates")
         return
      end if
      call check_setting(self%reference_charges, n, "reference charges", error)
      if (.not. allocated(error)) call check_setting(self%charge_thresholds, n, "charge thresholds", error)
      if (.not. allocated(error)) call check_setting(self%lower_bounds, n, "lower bounds", error)
      if (.not. allocated(error)) call check_setting(self%upper_bounds, n, "upper bounds", error)
      if (allocated(error)) return
      allocate (self%reference(n), self%threshold(n), self%lower(n), self%upper(n))
      self%reference = mol%charge/real(n, wp)
      if (allocated(self%reference_charges)) self%reference = self%reference_charges
      do i = 1, n
         self%threshold(i) = 0.5_wp*real(mol%num(mol%id(i)), wp)
         if (mol%num(mol%id(i)) == 1 .or. mol%num(mol%id(i)) == 3) then
            self%threshold(i) = real(mol%num(mol%id(i)), wp)
         end if
      end do
      if (allocated(self%charge_thresholds)) self%threshold = self%charge_thresholds
      self%lower = -huge(1.0_wp)/(4.0_wp*real(n, wp))
      self%upper = -self%lower
      if (allocated(self%lower_bounds)) self%lower = self%lower_bounds
      if (allocated(self%upper_bounds)) self%upper = self%upper_bounds
      if (any(self%threshold <= 0.0_wp) .or. any(self%lower > self%upper) .or. &
         & sum(self%lower) > mol%charge .or. sum(self%upper) < mol%charge) then
         call fatal_error(error, "EC fit has invalid thresholds or bounds incompatible with the total charge")
         return
      end if
      self%total_charge = mol%charge
      self%atom_xyz = mol%xyz
   end subroutine ec_charges_update

   !> Validate an optional per-atom setting
   !>
   !> @param[in]  values  optional allocated values (natom)
   !> @param[in]  natom   number of atoms
   !> @param[in]  label   diagnostic label
   !> @param[out] error   invalid length or nonfinite value
   subroutine check_setting(values, natom, label, error)
      !> Optional allocated values
      real(wp), allocatable, intent(in) :: values(:)
      !> Number of atoms
      integer, intent(in) :: natom
      !> Diagnostic label
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(values)) return
      if (size(values) /= natom .or. .not. all(ieee_is_finite(values))) then
         call fatal_error(error, "EC "//label//" must be finite with one entry per atom")
      end if
   end subroutine check_setting

   !> Refuse a radial EC fit
   !>
   !> @param[in,out] self   updated term
   !> @param[in]     grid   radial grid
   !> @param[out]    error  always set
   subroutine ec_charges_update_grid_1d(self, grid, error)
      !> Updated term
      class(ec_charges_type), intent(inout) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "EC charge fitting supports volume grids only")
   end subroutine ec_charges_update_grid_1d

   !> Build the fit on the volume grid of the update: mask, weights and restrained Hessian
   !>
   !> - Once per geometry; the grid is kept to refuse any other grid later
   !> - Hessian from per-thread block sums, added up serially
   !>
   !> @param[in,out] self   updated term
   !> @param[in]     grid   volume grid of the owning model, updated to the structure
   !> @param[out]    error  term not updated or unusable grid
   subroutine ec_charges_update_grid_3d(self, grid, error)
      !> Updated term
      class(ec_charges_type), intent(inout) :: self
      !> Volume grid of the owning model
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(ec_fit_type), allocatable :: fit
      !> Fit mask and weights; plain arrays inside the parallel regions, moved into the fit after
      real(wp), allocatable :: mask(:), weight(:)
      !> Per-thread Hessian sums (natom, natom, nthreads)
      real(wp), allocatable :: h_t(:, :, :)
      real(wp) :: center(3)
      integer :: n, a, ib, nblock, g0, g1, it, nt
      if (allocated(self%fit)) deallocate (self%fit)
      if (.not. allocated(self%atom_xyz)) then
         call fatal_error(error, "EC term is not updated")
         return
      end if
      if (grid%ngrid < 1 .or. .not. allocated(grid%xyz) .or. .not. allocated(grid%w)) then
         call fatal_error(error, "EC fit needs volume grid coordinates and weights")
         return
      end if
      if (any(shape(grid%xyz) /= [3, grid%ngrid]) .or. size(grid%w) /= grid%ngrid) then
         call fatal_error(error, "EC fit grid storage does not match its point count")
         return
      end if
      if (.not. all(ieee_is_finite(grid%xyz)) .or. .not. all(ieee_is_finite(grid%w)) .or. any(grid%w < 0.0_wp)) then
         call fatal_error(error, "EC fit needs finite coordinates and nonnegative volume weights")
         return
      end if
      n = size(self%atom_xyz, 2)
      nt = max(1, grid%team_size())
      allocate (mask(grid%ngrid))
      center = sum(self%atom_xyz, dim=2)/real(n, wp)
      nblock = (grid%ngrid + block_size - 1)/block_size
      !$omp parallel do default(shared) schedule(static) num_threads(nt) private(ib, g0, g1)
      do ib = 1, nblock
         g0 = (ib - 1)*block_size + 1
         g1 = min(ib*block_size, grid%ngrid)
         call mask_block(self, center, grid%xyz(:, g0:g1), mask(g0:g1))
      end do
      !$omp end parallel do
      allocate (fit)
      fit%normalization = sum(grid%w*mask)
      if (.not. ieee_is_finite(fit%normalization) .or. fit%normalization <= tiny(1.0_wp)) then
         call fatal_error(error, "EC fit has no positively weighted samples outside the nuclear exclusion")
         return
      end if
      weight = grid%w*mask/fit%normalization
      allocate (h_t(n, n, nt), source=0.0_wp)
      !$omp parallel do default(shared) schedule(static) num_threads(nt) private(ib, g0, g1, it)
      do ib = 1, nblock
         g0 = (ib - 1)*block_size + 1
         g1 = min(ib*block_size, grid%ngrid)
         it = 1
         !$ it = omp_get_thread_num() + 1
         call hessian_block(self, grid%xyz(:, g0:g1), weight(g0:g1), h_t(:, :, it))
      end do
      !$omp end parallel do
      fit%xyz = grid%xyz
      fit%volume = grid%w
      fit%nthreads = nt
      call move_alloc(mask, fit%mask)
      call move_alloc(weight, fit%weight)
      fit%h = sum(h_t, dim=3)
      do a = 1, n
         fit%h(a, a) = fit%h(a, a) + self%restraint_strength
      end do
      call move_alloc(fit, self%fit)
   end subroutine ec_charges_update_grid_3d

   !> Refuse a radial EC fit
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      radial grid
   !> @param[in,out] coupling  coupling in registration
   !> @param[out]    error     always set
   subroutine ec_charges_declare_1d(self, grid, coupling, error)
      !> Term
      class(ec_charges_type), intent(inout) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Coupling in registration
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "EC charge fitting supports volume grids only")
   end subroutine ec_charges_declare_1d

   !> Declare the host potential on the grid of the update
   !>
   !> - phi in every phase, dphi_dr in the gradient phase
   !> - No rebuild: the fit belongs to `update_grid`
   !>
   !> @param[in,out] self      term updated with a volume grid
   !> @param[in]     grid      volume grid
   !> @param[in,out] coupling  coupling in registration, scope already set by the potential
   !> @param[out]    error     fit not built, another grid or failed declaration
   subroutine ec_charges_declare_3d(self, grid, coupling, error)
      !> Term updated with a volume grid
      class(ec_charges_type), intent(inout) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Coupling in registration
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(point_potential_request_type) :: request
      integer :: ip
      call require_fit_grid(self, grid, error)
      if (allocated(error)) return
      do ip = 1, size(all_phases)
         call request_require(request, all_phases(ip), "phi", error)
         if (allocated(error)) return
      end do
      call request_require(request, moist_phase_gradient, "dphi_dr", error)
      if (.not. allocated(error)) call coupling_register(coupling, "potential", request, error)
   end subroutine ec_charges_declare_3d

   !> Require the fit and the very grid it was built on
   !>
   !> @param[in]  self   term
   !> @param[in]  grid   volume grid
   !> @param[out] error  fit not built or a grid differing from the grid of the update
   subroutine require_fit_grid(self, grid, error)
      !> Term
      class(ec_charges_type), intent(in) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(self%fit)) then
         call fatal_error(error, "EC fit is not built; update the potential with a volume grid")
         return
      end if
      if (.not. allocated(grid%xyz) .or. .not. allocated(grid%w)) then
         call fatal_error(error, "EC grid has no coordinates or weights")
         return
      end if
      if (any(shape(grid%xyz) /= shape(self%fit%xyz)) .or. size(grid%w) /= size(self%fit%volume)) then
         call fatal_error(error, "EC grid differs from the grid of its update")
         return
      end if
      if (any(abs(grid%xyz - self%fit%xyz) > 0.0_wp) .or. any(abs(grid%w - self%fit%volume) > 0.0_wp)) then
         call fatal_error(error, "EC grid differs from the grid of its update")
      end if
   end subroutine require_fit_grid

   !> Fit mask of a block of grid points: product of the exclusion switches times the far-field weight
   !>
   !> @param[in]  self    updated term
   !> @param[in]  center  unweighted nuclear centroid, bohr (3)
   !> @param[in]  xyz     grid points of the block, bohr (3, nb)
   !> @param[out] mask    fit mask (nb)
   subroutine mask_block(self, center, xyz, mask)
      !> Updated term
      class(ec_charges_type), intent(in) :: self
      !> Unweighted nuclear centroid
      real(wp), intent(in) :: center(3)
      !> Grid points of the block
      real(wp), intent(in) :: xyz(:, :)
      !> Fit mask
      real(wp), intent(out) :: mask(:)
      real(wp) :: r, switch, deriv, product_switch
      integer :: i, a
      do i = 1, size(xyz, 2)
         product_switch = 1.0_wp
         do a = 1, size(self%atom_xyz, 2)
            call exclusion_switch(self%exclusion_radius, norm2(xyz(:, i) - self%atom_xyz(:, a)), switch, deriv)
            product_switch = product_switch*switch
         end do
         r = sum((xyz(:, i) - center)**2)
         mask(i) = product_switch*(1.0_wp + self%long_range_bias*r/(r + self%range_scale**2))
      end do
   end subroutine mask_block

   !> Add the weighted outer products of a block's basis rows to a Hessian sum
   !>
   !> @param[in]     self    updated term
   !> @param[in]     xyz     grid points of the block, bohr (3, nb)
   !> @param[in]     weight  normalized fitting weights of the block (nb)
   !> @param[in,out] h       Hessian sum, 1/bohr**2 (natom, natom)
   subroutine hessian_block(self, xyz, weight, h)
      !> Updated term
      class(ec_charges_type), intent(in) :: self
      !> Grid points of the block
      real(wp), intent(in) :: xyz(:, :)
      !> Normalized fitting weights of the block
      real(wp), intent(in) :: weight(:)
      !> Hessian sum
      real(wp), intent(inout) :: h(:, :)
      ! TODO: check ifx stack/heap handling of this automatic array inside the OpenMP region
      real(wp) :: row(size(self%atom_xyz, 2))
      integer :: i, a, b
      do i = 1, size(xyz, 2)
         if (weight(i) <= 0.0_wp) cycle
         do a = 1, size(row)
            row(a) = 1.0_wp/norm2(xyz(:, i) - self%atom_xyz(:, a))
         end do
         do b = 1, size(row)
            do a = 1, size(row)
               h(a, b) = h(a, b) + weight(i)*row(a)*row(b)
            end do
         end do
      end do
   end subroutine hessian_block

   !> Basis product y = A x on the fit grid; zero where the mask vanishes
   !>
   !> @param[in]  self  term with a built fit
   !> @param[in]  x     atom values (natom)
   !> @param[out] y     grid values (ngrid)
   subroutine apply_basis(self, x, y)
      !> Term with a built fit
      class(ec_charges_type), intent(in) :: self
      !> Atom values
      real(wp), intent(in) :: x(:)
      !> Grid values
      real(wp), allocatable, intent(out) :: y(:)
      integer :: g, a
      allocate (y(size(self%fit%mask)))
      !$omp parallel do default(shared) schedule(static) num_threads(self%fit%nthreads) private(g, a)
      do g = 1, size(y)
         y(g) = 0.0_wp
         if (self%fit%mask(g) <= 0.0_wp) cycle
         do a = 1, size(x)
            y(g) = y(g) + x(a)/norm2(self%fit%xyz(:, g) - self%atom_xyz(:, a))
         end do
      end do
      !$omp end parallel do
   end subroutine apply_basis

   !> Transposed basis product x = A^T y on the fit grid, from per-thread block sums
   !>
   !> @param[in]  self  term with a built fit
   !> @param[in]  y     grid values (ngrid)
   !> @param[out] x     atom values (natom)
   subroutine apply_basis_transpose(self, y, x)
      !> Term with a built fit
      class(ec_charges_type), intent(in) :: self
      !> Grid values
      real(wp), intent(in) :: y(:)
      !> Atom values
      real(wp), allocatable, intent(out) :: x(:)
      !> Per-thread sums (natom, nthreads)
      real(wp), allocatable :: x_t(:, :)
      integer :: ib, nblock, g, a, it
      allocate (x_t(size(self%atom_xyz, 2), self%fit%nthreads), source=0.0_wp)
      nblock = (size(y) + block_size - 1)/block_size
      !$omp parallel do default(shared) schedule(static) num_threads(self%fit%nthreads) private(ib, g, a, it)
      do ib = 1, nblock
         it = 1
         !$ it = omp_get_thread_num() + 1
         do g = (ib - 1)*block_size + 1, min(ib*block_size, size(y))
            if (self%fit%mask(g) <= 0.0_wp) cycle
            do a = 1, size(x_t, 1)
               x_t(a, it) = x_t(a, it) + y(g)/norm2(self%fit%xyz(:, g) - self%atom_xyz(:, a))
            end do
         end do
      end do
      !$omp end parallel do
      x = sum(x_t, dim=2)
   end subroutine apply_basis_transpose

   !> Smooth exclusion switch of one nucleus and its radial derivative
   !>
   !> - Quintic smoothstep in t = (r - radius)/radius: zero inside the radius, one beyond twice it
   !>
   !> @param[in]  radius  exclusion radius, bohr
   !> @param[in]  r       distance to the nucleus, bohr
   !> @param[out] switch  switch value
   !> @param[out] deriv   radial derivative, 1/bohr
   pure subroutine exclusion_switch(radius, r, switch, deriv)
      !> Exclusion radius
      real(wp), intent(in) :: radius
      !> Distance to the nucleus
      real(wp), intent(in) :: r
      !> Switch value
      real(wp), intent(out) :: switch
      !> Radial derivative
      real(wp), intent(out) :: deriv
      real(wp) :: t
      t = (r - radius)/radius
      switch = 0.0_wp
      deriv = 0.0_wp
      if (t <= 0.0_wp) return
      switch = 1.0_wp
      if (t >= 1.0_wp) return
      switch = t**3*(10.0_wp + t*(-15.0_wp + 6.0_wp*t))
      deriv = 30.0_wp*t**2*(1.0_wp - t)**2/radius
   end subroutine exclusion_switch

   !> Smooth exclusion switches and their radial derivatives
   !>
   !> @param[in]  self    term
   !> @param[in]  point   volume point, bohr (3)
   !> @param[out] switch  switches per nucleus (natom)
   !> @param[out] deriv   radial derivatives, 1/bohr (natom)
   subroutine exclusion(self, point, switch, deriv)
      !> Term
      class(ec_charges_type), intent(in) :: self
      !> Volume point
      real(wp), intent(in) :: point(3)
      !> Switches per nucleus
      real(wp), intent(out) :: switch(:)
      !> Radial derivatives
      real(wp), intent(out) :: deriv(:)
      integer :: a
      do a = 1, size(switch)
         call exclusion_switch(self%exclusion_radius, norm2(point - self%atom_xyz(:, a)), switch(a), deriv(a))
      end do
   end subroutine exclusion

   !> Read and validate the current host potential
   !>
   !> @param[in]  self   declared term
   !> @param[in]  view   scoped coupling view
   !> @param[out] phi    potential, Hartree/e (ngrid)
   !> @param[out] error  missing fit, missing potential or invalid values
   subroutine read_potential(self, view, phi, error)
      !> Declared term
      class(ec_charges_type), intent(in) :: self
      !> Scoped coupling view
      type(coupling_view_type), intent(in) :: view
      !> Potential
      real(wp), allocatable, intent(out) :: phi(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(self%fit)) then
         call fatal_error(error, "EC fit is not built; update the potential with a volume grid")
         return
      end if
      call view%read("potential", "phi", phi, error)
      if (allocated(error)) return
      if (size(phi) /= size(self%fit%weight)) then
         call fatal_error(error, "EC host potential does not match the fit grid")
      else if (.not. all(ieee_is_finite(phi))) then
         call fatal_error(error, "EC host potential must be finite")
      end if
      if (allocated(error)) deallocate (phi)
   end subroutine read_potential

   !> Fit tail charges from the current host potential
   !>
   !> @param[in]  self   declared term
   !> @param[in]  view   scoped coupling view
   !> @param[out] q      fitted charges, e (natom)
   !> @param[out] error  invalid host potential or failed constrained solve
   subroutine ec_charges_tail(self, view, q, error)
      !> Declared term
      class(ec_charges_type), intent(in) :: self
      !> Scoped coupling view
      type(coupling_view_type), intent(in) :: view
      !> Fitted charges
      real(wp), allocatable, intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: phi(:)
      integer, allocatable :: active(:)
      ! TODO: cache q and the active set per host answer; tail_charges and the adjoint refit every call
      call read_potential(self, view, phi, error)
      if (.not. allocated(error)) call solve_fit(self, phi, q, active, error)
   end subroutine ec_charges_tail

   !> Supply the full host field, independently of the fitted tail charges
   !>
   !> @param[in]  self     declared term
   !> @param[in]  grid     volume grid
   !> @param[in]  view     scoped coupling view
   !> @param[out] phi      host potential, Hartree/e (ngrid)
   !> @param[out] error    fit not built, another grid, or missing or invalid host outputs
   !> @param[out] dphi_dr  optional point derivative, Hartree/(e bohr) (3, ngrid)
   subroutine ec_charges_field(self, grid, view, phi, error, dphi_dr)
      !> Declared term
      class(ec_charges_type), intent(in) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Scoped coupling view
      type(coupling_view_type), intent(in) :: view
      !> Host potential
      real(wp), allocatable, intent(out) :: phi(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional point derivative
      real(wp), allocatable, intent(out), optional :: dphi_dr(:, :)
      call require_fit_grid(self, grid, error)
      if (allocated(error)) return
      call read_potential(self, view, phi, error)
      if (allocated(error)) return
      if (.not. present(dphi_dr)) return
      call view%read("potential", "dphi_dr", dphi_dr, error)
      if (allocated(error)) return
      if (any(shape(dphi_dr) /= [3, grid%ngrid])) then
         call fatal_error(error, "EC potential derivative does not match the volume grid")
         return
      end if
      if (.not. all(ieee_is_finite(dphi_dr))) call fatal_error(error, "EC potential derivative must be finite")
   end subroutine ec_charges_field

   !> Reverse the constrained fit and the full-field evaluation
   !>
   !> Derivatives use the current active constraints; charge-penalty kinks are piecewise differentiable
   !>
   !> @param[in,out] self          declared term
   !> @param[in]     grid          volume grid
   !> @param[in]     view          scoped coupling view
   !> @param[in]     w_phi         adjoint of the full field (ngrid)
   !> @param[in]     w_qinf        adjoint of fitted tail charges (natom)
   !> @param[in,out] response      potential adjoint for the host
   !> @param[out]    error         invalid shape, stale geometry or failed solve
   !> @param[in,out] grid_adjoint  optional position and quadrature-weight adjoints
   !> @param[in,out] gradient      optional explicit nuclear gradient, Hartree/bohr (3, natom)
   subroutine ec_charges_adjoint(self, grid, view, w_phi, w_qinf, response, error, grid_adjoint, gradient)
      !> Declared term
      class(ec_charges_type), intent(inout) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Scoped coupling view
      type(coupling_view_type), intent(in) :: view
      !> Adjoint of the full field
      real(wp), intent(in) :: w_phi(:)
      !> Adjoint of fitted tail charges
      real(wp), intent(in) :: w_qinf(:)
      !> Potential adjoint for the host
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional volume adjoints
      type(volume_adjoint_type), intent(inout), optional :: grid_adjoint
      !> Optional explicit nuclear gradient
      real(wp), intent(inout), optional :: gradient(:, :)
      type(potential_adjoint_response_type) :: item
      real(wp), allocatable :: phi(:), dphi(:, :), q(:), v(:), residual(:), t(:), w_raw(:), w_xyz(:, :), w_w(:)
      real(wp), allocatable :: switches(:), deriv(:), zeros(:)
      integer, allocatable :: active(:)
      real(wp) :: multiplier, center(3), d(3), dr(3), r, far, radial, other, w_a, mask_deriv, mean_weight
      integer :: n, g, a, b
      ! Grid motion needs dphi_dr, a gradient-phase output; refused elsewhere instead of zeros
      if (present(grid_adjoint)) then
         call ec_charges_field(self, grid, view, phi, error, dphi)
      else
         call ec_charges_field(self, grid, view, phi, error)
      end if
      if (allocated(error)) return
      n = size(self%atom_xyz, 2)
      if (size(w_phi) /= grid%ngrid .or. size(w_qinf) /= n) then
         call fatal_error(error, "EC field or tail-charge adjoint has the wrong length")
         return
      end if
      if (present(gradient)) then
         if (any(shape(gradient) /= [3, n])) then
            call fatal_error(error, "EC nuclear gradient has the wrong shape")
            return
         end if
      end if
      ! TODO: reuse the charges and active set of tail_charges instead of refitting
      call solve_fit(self, phi, q, active, error)
      if (allocated(error)) return
      allocate (zeros(n), source=0.0_wp)
      call solve_face(self%fit%h, w_qinf, zeros, active, 0.0_wp, v, multiplier, error)
      if (allocated(error)) return
      call apply_basis(self, q, residual)
      residual = phi - residual
      call apply_basis(self, v, t)
      item%w_phi = w_phi + self%fit%weight*t
      if (present(grid_adjoint) .or. present(gradient)) then
         allocate (w_xyz(3, grid%ngrid), source=0.0_wp)
         allocate (switches(n), deriv(n))
         mean_weight = sum(self%fit%weight*residual*t)
         w_raw = (residual*t - mean_weight)/self%fit%normalization
         w_w = self%fit%mask*w_raw
         center = sum(self%atom_xyz, dim=2)/real(n, wp)
         ! TODO: parallelize this per-point loop (per-thread gradient buffers); last serial O(ngrid natom) part
         do g = 1, grid%ngrid
            call exclusion(self, grid%xyz(:, g), switches, deriv)
            d = grid%xyz(:, g) - center
            r = sum(d**2)
            far = 1.0_wp + self%long_range_bias*r/(r + self%range_scale**2)
            dr = w_raw(g)*grid%w(g)*product(switches)* &
               & 2.0_wp*self%long_range_bias*self%range_scale**2*d/(r + self%range_scale**2)**2
            w_xyz(:, g) = dr
            if (present(gradient)) then
               do a = 1, n
                  gradient(:, a) = gradient(:, a) - dr/real(n, wp)
               end do
            end if
            do a = 1, n
               d = grid%xyz(:, g) - self%atom_xyz(:, a)
               r = norm2(d)
               dr = 0.0_wp
               if (self%fit%mask(g) > 0.0_wp) then
                  w_a = self%fit%weight(g)*(residual(g)*v(a) - t(g)*q(a))
                  dr = -w_a*d/r**3
               end if
               if (deriv(a) > 0.0_wp) then
                  other = 1.0_wp
                  do b = 1, n
                     if (b /= a) other = other*switches(b)
                  end do
                  radial = w_raw(g)*grid%w(g)*far*other*deriv(a)
                  mask_deriv = radial/r
                  dr = dr + mask_deriv*d
               end if
               w_xyz(:, g) = w_xyz(:, g) + dr
               if (present(gradient)) gradient(:, a) = gradient(:, a) - dr
            end do
         end do
         if (present(grid_adjoint)) then
            do g = 1, grid%ngrid
               w_xyz(:, g) = w_xyz(:, g) + item%w_phi(g)*dphi(:, g)
            end do
            call grid_adjoint%add_weights(error, w_xyz=w_xyz, w_w=w_w)
            if (allocated(error)) return
         end if
      end if
      call response_accumulate(response, item, error)
   end subroutine ec_charges_adjoint

   !> Minimize the restrained fit with exact total charge and piecewise linear penalties
   !>
   !> Active codes: -1/1 hard lower/upper, -2/2 negative/positive penalty kink, 3 fixed
   !>
   !> @param[in]  self    declared term
   !> @param[in]  phi     host potential, Hartree/e (ngrid)
   !> @param[out] q       fitted charges, e (natom)
   !> @param[out] active  active constraints (natom)
   !> @param[out] error   failed constrained solve
   subroutine solve_fit(self, phi, q, active, error)
      !> Declared term
      class(ec_charges_type), intent(in) :: self
      !> Host potential
      real(wp), intent(in) :: phi(:)
      !> Fitted charges
      real(wp), allocatable, intent(out) :: q(:)
      !> Active constraints
      integer, allocatable, intent(out) :: active(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: b(:), candidate(:), direction(:), rhs(:), lo(:), hi(:), grad(:)
      integer, allocatable :: region(:)
      real(wp) :: left, right, shift, step, alpha, ratio, nu, tol, low_nu, high_nu, gmin, gmax, violation, worst
      integer :: n, i, iter, hit, release, low_i, high_i, hit_code, release_dir
      n = size(self%reference)
      call apply_basis_transpose(self, self%fit%weight*phi, b)
      b = b + self%restraint_strength*self%reference
      allocate (q(n), active(n), region(n), lo(n), hi(n))
      q = min(self%upper, max(self%lower, self%reference))
      left = 0.0_wp
      right = 0.0_wp
      step = max(1.0_wp, abs(self%total_charge), maxval(abs(self%reference)))
      if (sum(q) < self%total_charge) then
         right = step
         do while (sum(min(self%upper, max(self%lower, self%reference + right))) < self%total_charge)
            right = 2.0_wp*right
         end do
      else if (sum(q) > self%total_charge) then
         left = -step
         do while (sum(min(self%upper, max(self%lower, self%reference + left))) > self%total_charge)
            left = 2.0_wp*left
         end do
      end if
      do iter = 1, 100
         shift = 0.5_wp*(left + right)
         q = min(self%upper, max(self%lower, self%reference + shift))
         if (sum(q) < self%total_charge) then
            left = shift
         else
            right = shift
         end if
      end do
      active = 0
      region = 0
      where (q < -self%threshold) region = -1
      where (q > self%threshold) region = 1
      where (q <= self%lower) active = -1
      where (q >= self%upper) active = 1
      where (self%lower >= self%upper) active = 3
      tol = 256.0_wp*epsilon(1.0_wp)*max(1.0_wp, maxval(abs(b)), self%charge_penalty)
      do iter = 1, 1000*n + 100
         if (any(active == 0)) then
            rhs = b - self%charge_penalty*real(region, wp)
            call solve_face(self%fit%h, rhs, q, active, self%total_charge, candidate, nu, error)
            if (allocated(error)) return
            direction = candidate - q
            lo = self%lower
            hi = self%upper
            where (region == -1) hi = min(hi, -self%threshold)
            where (region == 0) lo = max(lo, -self%threshold)
            where (region == 0) hi = min(hi, self%threshold)
            where (region == 1) lo = max(lo, self%threshold)
            alpha = 1.0_wp
            hit = 0
            hit_code = 0
            do i = 1, n
               if (active(i) /= 0) cycle
               if (candidate(i) < lo(i) .and. direction(i) < 0.0_wp) then
                  ratio = max(0.0_wp, (lo(i) - q(i))/direction(i))
                  if (ratio >= alpha) cycle
                  alpha = ratio
                  hit = i
                  hit_code = -1
                  if (lo(i) > self%lower(i)) hit_code = 2*region(i) - 2
                  if (region(i) == 1 .and. lo(i) > self%lower(i)) hit_code = 2
               else if (candidate(i) > hi(i) .and. direction(i) > 0.0_wp) then
                  ratio = max(0.0_wp, (hi(i) - q(i))/direction(i))
                  if (ratio >= alpha) cycle
                  alpha = ratio
                  hit = i
                  hit_code = 1
                  if (hi(i) < self%upper(i)) hit_code = 2
                  if (region(i) == -1 .and. hi(i) < self%upper(i)) hit_code = -2
               end if
            end do
            q = q + alpha*direction
            if (hit > 0) then
               active(hit) = hit_code
               select case (hit_code)
               case (-1)
                  q(hit) = self%lower(hit)
               case (1)
                  q(hit) = self%upper(hit)
               case (-2)
                  q(hit) = -self%threshold(hit)
               case (2)
                  q(hit) = self%threshold(hit)
               case default
                  call fatal_error(error, "EC fit encountered an invalid active constraint")
                  return
               end select
               cycle
            end if
         else
            ! All fixed: intersect the multiplier intervals of the active constraints
            grad = matmul(self%fit%h, q) - b
            low_nu = -huge(1.0_wp)
            high_nu = huge(1.0_wp)
            low_i = 0
            high_i = 0
            do i = 1, n
               if (active(i) == 3) cycle
               call penalty_interval(q(i), self%threshold(i), self%charge_penalty, gmin, gmax)
               gmin = gmin + grad(i)
               gmax = gmax + grad(i)
               if (active(i) /= 1 .and. -gmax > low_nu) then
                  low_nu = -gmax
                  low_i = i
               end if
               if (active(i) /= -1 .and. -gmin < high_nu) then
                  high_nu = -gmin
                  high_i = i
               end if
            end do
            if (low_nu <= high_nu + tol) return
            active(low_i) = 0
            active(high_i) = 0
            call release_region(q(low_i), self%threshold(low_i), 1, region(low_i))
            call release_region(q(high_i), self%threshold(high_i), -1, region(high_i))
            cycle
         end if
         grad = matmul(self%fit%h, q) - b + nu
         release = 0
         release_dir = 0
         worst = tol
         do i = 1, n
            if (active(i) == 0 .or. active(i) == 3) cycle
            call penalty_interval(q(i), self%threshold(i), self%charge_penalty, gmin, gmax)
            gmin = gmin + grad(i)
            gmax = gmax + grad(i)
            violation = 0.0_wp
            if (active(i) /= 1) violation = max(violation, -gmax)
            if (active(i) /= -1) violation = max(violation, gmin)
            if (violation <= worst) cycle
            worst = violation
            release = i
            release_dir = 1
            if (gmin > 0.0_wp) release_dir = -1
         end do
         if (release == 0) return
         active(release) = 0
         call release_region(q(release), self%threshold(release), release_dir, region(release))
      end do
      call fatal_error(error, "EC constrained ESP fit did not converge")
   end subroutine solve_fit

   !> Subgradient interval of the linear excess-charge penalty
   !>
   !> @param[in]  q          charge, e
   !> @param[in]  threshold  excess-charge threshold, e
   !> @param[in]  slope      penalty slope, e/bohr**2
   !> @param[out] lower      lower subgradient
   !> @param[out] upper      upper subgradient
   pure subroutine penalty_interval(q, threshold, slope, lower, upper)
      !> Charge
      real(wp), intent(in) :: q
      !> Excess-charge threshold
      real(wp), intent(in) :: threshold
      !> Penalty slope
      real(wp), intent(in) :: slope
      !> Lower subgradient
      real(wp), intent(out) :: lower
      !> Upper subgradient
      real(wp), intent(out) :: upper
      lower = 0.0_wp
      upper = 0.0_wp
      if (q <= -threshold) lower = -slope
      if (q < -threshold) upper = -slope
      if (q >= threshold) upper = slope
      if (q > threshold) lower = slope
   end subroutine penalty_interval

   !> Choose the smooth charge region on release from a constraint
   !>
   !> @param[in]  q          charge, e
   !> @param[in]  threshold  excess-charge threshold, e
   !> @param[in]  direction  release direction, -1 or 1
   !> @param[out] region     smooth penalty region, -1, 0 or 1
   pure subroutine release_region(q, threshold, direction, region)
      !> Charge
      real(wp), intent(in) :: q
      !> Excess-charge threshold
      real(wp), intent(in) :: threshold
      !> Release direction
      integer, intent(in) :: direction
      !> Smooth penalty region
      integer, intent(out) :: region
      region = 0
      if (direction < 0) then
         if (q <= -threshold) region = -1
         if (q > threshold) region = 1
      else
         if (q < -threshold) region = -1
         if (q >= threshold) region = 1
      end if
   end subroutine release_region

   !> Solve an equality-constrained quadratic on the free charge face
   !>
   !> @param[in]  h           positive definite Hessian (natom, natom)
   !> @param[in]  b           right-hand side (natom)
   !> @param[in]  fixed       fixed charge values (natom)
   !> @param[in]  active      active constraints; zero for free charges (natom)
   !> @param[in]  total       required charge sum
   !> @param[out] q           solved charges (natom)
   !> @param[out] multiplier  equality multiplier
   !> @param[out] error       factorization or solve failure
   subroutine solve_face(h, b, fixed, active, total, q, multiplier, error)
      !> Positive definite Hessian
      real(wp), intent(in) :: h(:, :)
      !> Right-hand side
      real(wp), intent(in) :: b(:)
      !> Fixed charge values
      real(wp), intent(in) :: fixed(:)
      !> Active constraints
      integer, intent(in) :: active(:)
      !> Required charge sum
      real(wp), intent(in) :: total
      !> Solved charges
      real(wp), allocatable, intent(out) :: q(:)
      !> Equality multiplier
      real(wp), intent(out) :: multiplier
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: factor(:, :), rhs(:, :)
      integer, allocatable :: free(:)
      integer :: i, nf, info
      q = fixed
      multiplier = 0.0_wp
      nf = count(active == 0)
      if (nf == 0) return
      free = pack([(i, i=1, size(active))], active == 0)
      q(free) = 0.0_wp
      factor = h(free, free)
      allocate (rhs(nf, 2))
      rhs(:, 1) = b(free) - matmul(h(free, :), q)
      rhs(:, 2) = 1.0_wp
      call potrf(factor, info)
      if (info == 0) call potrs(factor, rhs, info)
      if (info /= 0) then
         call fatal_error(error, "EC charge-fit factorization failed")
         return
      end if
      multiplier = (sum(rhs(:, 1)) + sum(q) - total)/sum(rhs(:, 2))
      q(free) = rhs(:, 1) - multiplier*rhs(:, 2)
   end subroutine solve_face

end module moist_model_moz_potential_electrostatic_ec
