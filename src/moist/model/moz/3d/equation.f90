!> 3D UV equation on any 3D grid: the map c_s(r) -> t_s(r) and its derivative
!>
!> - Per k-point: h_a = sum_g c_g P_ga(|k|) with the solvent response table
!> - Ng split: the map takes short-range c_s(r); c(k) = c_s(k) - beta u_lr(k)
!>   and t_s(k) = h(k) - c(k) - beta u_lr(k)
!>   The nonlinear driver iterates on t_s and supplies c_s through the closure
!> - The map is affine in c: dt = k2r(P^T dc(k) - dc(k)) needs no linearization
!> - Columns: one per solvent site, (ngrid, nv) in real and (npts_k, nv) in k space
!> - Transforms through the engine of the grid (`grid%new_trafo`), so the
!>   equation never sees the grid's layout; the equation owns a copy of the grid
!>   and prepares its engine during setup
!> - Assignment copies everything but the engine, which the copy rebuilds from
!>   its own grid; finalization releases the engine and its plans
module moist_model_moz_3d_equation
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: resolve_num_threads
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type, moist_math_grid_3d_trafo_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type
   use moist_model_moz_3d_kernel, only: uv3d_point_map, uv3d_point_jvp, uv3d_point_adjoint, uv3d_tail_adjoint
   use moist_model_moz_fields, only: copy_field
   implicit none(type, external)
   private

   public :: uv3d_equation_type

   !> Fixed 3D UV data, the owned grid and its transform engine
   !>
   !> Copy with assignment (=), allocating the destination first if allocatable.
   !> allocate(source=...) bypasses this assignment and may share native plans
   !> and the original grid binding; do not use it for prepared equations.
   type :: uv3d_equation_type
      !> OpenMP team size of the k-point loops and the transforms; zero takes omp_get_max_threads
      integer :: nthreads = 0
      !> Solvent sites, real-space points and k-points; 0 before setup
      integer, private :: nv = 0, ngrid = 0, nk = 0
      !> Inverse thermal energy, 1/Hartree
      real(wp), private :: beta = 0.0_wp
      !> Copy of the updated grid
      class(moist_math_grid_3d_type), allocatable, private :: grid
      !> Transform engine, prepared by setup; copies rebuild it when needed
      class(moist_math_grid_3d_trafo_type), allocatable, private :: trafo
      !> Solvent response per distinct |k|
      type(uv3d_susceptibility_type), private :: response
      !> Long-range reciprocal potential, Hartree bohr^3, (npts_k, nv)
      complex(wp), allocatable, private :: uk_lr(:, :)
      !> Real-space and k-space scratch
      real(wp), allocatable, private :: r_work(:, :)
      complex(wp), allocatable, private :: ck_work(:, :), tk_work(:, :)
   contains
      !> Bind the grid and the fixed data
      procedure :: setup => equation_setup
      !> Whether setup succeeded
      procedure :: ready => equation_ready
      !> Replace the tail table on the same grid, retaining transform plans and scratch
      procedure :: update_potential => equation_update_potential
      !> Map c_s(r) -> t_s(r)
      procedure :: map => equation_map
      !> Derivative dc(r) -> dt(r)
      procedure :: jvp => equation_jvp
      !> Euclidean transpose; optionally also return the u_lr(k) sensitivity
      procedure :: adjoint => equation_adjoint
      !> Volume integral of a real-space field by the owned grid's quadrature
      procedure :: integrate => equation_integrate
      !> Compensated contribution of a bounded field slice
      procedure :: integrate_block => equation_integrate_block
      !> Quadrature weights of the owned grid, bohr^3
      procedure :: integration_weights => equation_integration_weights
      !> Copy without the transform engine
      generic :: assignment(=) => assign
      procedure, private :: assign => equation_assign
      !> Build the transform engine of the owned grid unless present
      procedure, private :: ensure_trafo => equation_ensure_trafo
      !> Check the shapes of a pair of real-space fields
      procedure, private :: check_fields => equation_check_fields
      !> Primal map before the inverse transform, shared with its geometry pullback
      procedure, private :: map_reciprocal => equation_map_reciprocal
      final :: destroy_equation
   end type uv3d_equation_type

contains

   !> Copy an updated grid and bind the fixed UV data
   !>
   !> @param[inout] self      equation; a previous setup is dropped
   !> @param[in]    grid      updated grid
   !> @param[in]    response  solvent response on the grid's k-points
   !> @param[in]    uk_lr     long-range reciprocal potential, Hartree bohr^3, (npts_k, nv)
   !> @param[in]    beta      inverse thermal energy, 1/Hartree
   !> @param[out]   error     grid not updated, inconsistent shapes, invalid beta or engine failure
   subroutine equation_setup(self, grid, response, uk_lr, beta, error)
      !> Equation
      class(uv3d_equation_type), intent(inout), target :: self
      !> Updated grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Solvent response
      type(uv3d_susceptibility_type), intent(in) :: response
      !> Long-range reciprocal potential
      complex(wp), intent(in) :: uk_lr(:, :)
      !> Inverse thermal energy
      real(wp), intent(in) :: beta
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer :: nv, stat

      call release(self)
      if (grid%ngrid < 1 .or. grid%npts_k < 1) then
         call fatal_error(error, "3D UV equation needs an updated grid with k-points")
         return
      end if
      nv = response%nv
      if (nv < 1 .or. .not. allocated(response%values) .or. .not. allocated(response%index)) then
         call fatal_error(error, "3D UV equation needs a solvent response table")
         return
      end if
      if (size(response%index) /= grid%npts_k .or. size(response%values, 1) /= nv .or. &
         & size(response%values, 2) /= nv .or. minval(response%index) < 1 .or. &
         & maxval(response%index) > size(response%values, 3)) then
         call fatal_error(error, "3D UV solvent response does not match the grid's k-points")
         return
      end if
      if (any(shape(uk_lr) /= [grid%npts_k, nv])) then
         call fatal_error(error, "3D UV long-range potential must have shape (npts_k, nv)")
         return
      end if
      if (allocated(response%charge_response)) then
         if (.not. allocated(response%charge) .or. .not. allocated(response%radius)) then
            call fatal_error(error, "3D UV OZ charge response needs solvent charges and shell radii")
            return
         end if
         if (any(shape(response%charge_response) /= [nv, size(response%values, 3)]) .or. &
            & size(response%charge) /= nv .or. size(response%radius) /= size(response%values, 3)) then
            call fatal_error(error, "3D UV OZ charge response does not match the susceptibility table")
            return
         end if
      end if
      if (.not. ieee_is_finite(beta) .or. beta <= 0.0_wp) then
         call fatal_error(error, "3D UV equation needs a positive inverse thermal energy")
         return
      end if
      allocate (self%grid, source=grid, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to copy the 3D UV grid")
         return
      end if
      self%grid%nthreads = self%nthreads
      self%response = response
      self%uk_lr = uk_lr
      self%beta = beta
      self%ngrid = grid%ngrid
      self%nk = grid%npts_k
      allocate (self%r_work(self%ngrid, nv), self%ck_work(self%nk, nv), self%tk_work(self%nk, nv))
      self%nv = nv
      call self%ensure_trafo(error)
      if (allocated(error)) call release(self)
   end subroutine equation_setup

   !> Whether setup succeeded
   !>
   !> @param[in] self Equation
   pure logical function equation_ready(self) result(ready)
      !> Equation
      class(uv3d_equation_type), intent(in) :: self
      ready = self%nv > 0
   end function equation_ready

   !> Replace only the reciprocal potential; geometry, response and beta stay fixed
   !>
   !> @param[inout] self   equation after setup
   !> @param[in]    uk_lr  long-range potential, Hartree bohr^3, (npts_k, nv)
   !> @param[out]   error  unset equation or shape mismatch; leaves the old table intact
   subroutine equation_update_potential(self, uk_lr, error)
      class(uv3d_equation_type), intent(inout) :: self
      complex(wp), intent(in) :: uk_lr(:, :)
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%ready()) then
         call fatal_error(error, "Update a set-up 3D UV equation's potential")
      else if (any(shape(uk_lr) /= [self%nk, self%nv])) then
         call fatal_error(error, "3D UV long-range potential must have shape (npts_k, nv)")
      else
         self%uk_lr(:, :) = uk_lr
      end if
   end subroutine equation_update_potential

   !> Map the short-range direct correlation c_s(r) to t_s(r)
   !>
   !> @param[inout] self   equation after setup
   !> @param[in]    cr     c_s(r), (ngrid, nv)
   !> @param[out]   tr     t_s(r), (ngrid, nv)
   !> @param[out]   error  not set up, shape mismatch or transform failure
   subroutine equation_map(self, cr, tr, error)
      !> Equation
      class(uv3d_equation_type), intent(inout), target :: self
      !> Short-range direct correlation
      real(wp), intent(in) :: cr(:, :)
      !> Short-range indirect correlation
      real(wp), intent(out), contiguous, target :: tr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%check_fields(cr, tr, error)
      if (allocated(error)) return
      call self%ensure_trafo(error)
      if (allocated(error)) return
      call self%map_reciprocal(cr, error)
      if (allocated(error)) return
      call self%trafo%fft_k2r(self%tk_work, tr, error)
   end subroutine equation_map

   !> Transform c_s and apply the same per-k map as the primal equation.
   subroutine equation_map_reciprocal(self, cr, error)
      !> Equation with a prepared engine
      class(uv3d_equation_type), intent(inout) :: self
      !> Primal short-range direct correlation
      real(wp), intent(in) :: cr(:, :)
      !> Transform failure
      type(error_type), allocatable, intent(out) :: error
      integer :: ik

      call copy_field(cr, self%r_work, resolve_num_threads(self%nthreads))
      call self%trafo%fft_r2k(self%r_work, self%ck_work, error)
      if (allocated(error)) return
      !$omp parallel do schedule(static) num_threads(resolve_num_threads(self%nthreads))
      do ik = 1, self%nk
         if (allocated(self%response%charge_response)) then
            call uv3d_point_map(self%response%values(:, :, self%response%index(ik)), self%ck_work(ik, :), &
               & self%uk_lr(ik, :), self%beta, self%tk_work(ik, :), self%response%charge, &
               & self%response%charge_response(:, self%response%index(ik)), &
               & self%response%radius(self%response%index(ik))**2)
         else
            call uv3d_point_map(self%response%values(:, :, self%response%index(ik)), self%ck_work(ik, :), &
               & self%uk_lr(ik, :), self%beta, self%tk_work(ik, :))
         end if
      end do
      !$omp end parallel do
   end subroutine equation_map_reciprocal

   !> Derivative of the map: dt(r) = k2r(P^T dc(k) - dc(k)), dc(k) = r2k(dc(r))
   !>
   !> @param[inout] self   equation after setup
   !> @param[in]    dcr    direction dc(r), (ngrid, nv)
   !> @param[out]   dtr    response dt(r), (ngrid, nv)
   !> @param[out]   error  not set up, shape mismatch or transform failure
   subroutine equation_jvp(self, dcr, dtr, error)
      !> Equation
      class(uv3d_equation_type), intent(inout), target :: self
      !> Direction
      real(wp), intent(in) :: dcr(:, :)
      !> Response
      real(wp), intent(out), contiguous, target :: dtr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer :: ik

      call self%check_fields(dcr, dtr, error)
      if (allocated(error)) return
      call self%ensure_trafo(error)
      if (allocated(error)) return
      call copy_field(dcr, self%r_work, resolve_num_threads(self%nthreads))
      call self%trafo%fft_r2k(self%r_work, self%ck_work, error)
      if (allocated(error)) return
      !$omp parallel do schedule(static) num_threads(resolve_num_threads(self%nthreads))
      do ik = 1, self%nk
         call uv3d_point_jvp(self%response%values(:, :, self%response%index(ik)), self%ck_work(ik, :), &
            & self%tk_work(ik, :))
      end do
      !$omp end parallel do
      call self%trafo%fft_k2r(self%tk_work, dtr, error)
   end subroutine equation_jvp

   !> Transpose of c_s(r) -> t_s(r), optionally differentiating u_lr(k)
   !>
   !> Real fields use sum(w*direction); complex fields use
   !> Re(sum(conjg(wk)*direction)). The optional wk includes -beta.
   !> cr and grid_adjoint together request the geometry pullback of the primal
   !> affine map, with the bulk response and reciprocal lattice held fixed.
   subroutine equation_adjoint(self, wt, wc, error, wk, cr, grid_adjoint)
      class(uv3d_equation_type), intent(inout), target :: self
      real(wp), intent(in) :: wt(:, :)
      real(wp), intent(out), contiguous, target :: wc(:, :)
      type(error_type), allocatable, intent(out) :: error
      complex(wp), allocatable, optional, intent(out) :: wk(:, :)
      real(wp), optional, intent(in) :: cr(:, :)
      type(volume_adjoint_type), optional, intent(inout) :: grid_adjoint
      integer :: ik, shell

      call self%check_fields(wt, wc, error)
      if (allocated(error)) return
      if (present(cr) .neqv. present(grid_adjoint)) then
         call fatal_error(error, "3D UV geometry pullback needs both primal correlations and a grid accumulator")
         return
      end if
      if (present(grid_adjoint)) then
         call self%check_fields(cr, wc, error)
         if (allocated(error)) return
         call self%grid%check_volume_adjoint(grid_adjoint, error)
         if (allocated(error)) return
      end if
      call self%ensure_trafo(error)
      if (allocated(error)) return
      if (present(grid_adjoint)) then
         call self%map_reciprocal(cr, error)
         if (allocated(error)) return
         call self%trafo%accumulate_k2r_geometry(self%tk_work, wt, grid_adjoint, error)
         if (allocated(error)) return
      end if
      call copy_field(wt, self%r_work, resolve_num_threads(self%nthreads))
      call self%trafo%fft_k2r_adj(self%r_work, self%tk_work, error)
      if (allocated(error)) return
      if (present(wk)) allocate (wk(self%nk, self%nv))
      !$omp parallel do schedule(static) num_threads(resolve_num_threads(self%nthreads)) private(shell)
      do ik = 1, self%nk
         shell = self%response%index(ik)
         call uv3d_point_adjoint(self%response%values(:, :, shell), self%tk_work(ik, :), self%ck_work(ik, :))
         if (present(wk)) then
            if (allocated(self%response%charge_response)) then
               call uv3d_tail_adjoint(self%response%values(:, :, shell), self%tk_work(ik, :), self%beta, wk(ik, :), &
                  & self%response%charge, self%response%charge_response(:, shell), self%response%radius(shell)**2)
            else
               call uv3d_tail_adjoint(self%response%values(:, :, shell), self%tk_work(ik, :), self%beta, wk(ik, :))
            end if
         end if
      end do
      !$omp end parallel do
      if (present(grid_adjoint)) then
         call self%trafo%accumulate_r2k_geometry(cr, self%ck_work, grid_adjoint, error)
         if (allocated(error)) return
      end if
      call self%trafo%fft_r2k_adj(self%ck_work, wc, error)
   end subroutine equation_adjoint

   !> Copy the owned grid's final quadrature weights for response evaluation
   !>
   !> @param[in] self Prepared equation
   !> @param[out] weights Integration weights, bohr**3 (ngrid)
   !> @param[out] error Unprepared equation or invalid weight storage
   subroutine equation_integration_weights(self, weights, error)
      !> Prepared equation
      class(uv3d_equation_type), intent(in) :: self
      !> Integration weights, bohr**3 (ngrid)
      real(wp), allocatable, intent(out) :: weights(:)
      !> Unprepared equation or invalid weight storage
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%ready()) then
         call fatal_error(error, "3D UV equation is not set up")
         return
      end if
      if (.not. allocated(self%grid%w)) then
         call fatal_error(error, "3D UV response needs grid integration weights")
         return
      end if
      if (size(self%grid%w) /= self%ngrid) then
         call fatal_error(error, "3D UV integration weights must have ngrid entries")
         return
      end if
      weights = self%grid%w
   end subroutine equation_integration_weights

   !> Volume integral of a real-space field by the owned grid's quadrature
   !>
   !> @param[in]  self    equation after setup
   !> @param[in]  f       field on the grid points, (ngrid)
   !> @param[out] result  integral, field units times bohr^3
   !> @param[out] error   not set up or size mismatch
   subroutine equation_integrate(self, f, result, error)
      !> Equation
      class(uv3d_equation_type), intent(in) :: self
      !> Field
      real(wp), intent(in) :: f(:)
      !> Integral
      real(wp), intent(out) :: result
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      result = 0.0_wp
      if (.not. self%ready()) then
         call fatal_error(error, "3D UV equation is not set up")
         return
      end if
      call self%grid%integrate_field(f, result, error)
   end subroutine equation_integrate

   !> Integrate a slice through the owned grid without a whole-field temporary
   !>
   !> @param[in] self Prepared equation
   !> @param[in] f Field values at first:first+size(f)-1
   !> @param[in] first Absolute grid index of the first slice value
   !> @param[out] partial Compensated total and correction, in integral units
   !> @param[out] error Unprepared equation or invalid slice
   subroutine equation_integrate_block(self, f, first, partial, error)
      class(uv3d_equation_type), intent(in) :: self
      real(wp), intent(in) :: f(:)
      integer, intent(in) :: first
      real(wp), intent(out) :: partial(2)
      type(error_type), allocatable, intent(out) :: error

      partial = 0.0_wp
      if (.not. self%ready()) then
         call fatal_error(error, "3D UV equation is not set up")
         return
      end if
      call self%grid%integrate_field_block(f, first, partial, error)
   end subroutine equation_integrate_block

   !> Check that the equation is set up and both fields are (ngrid, nv)
   !>
   !> @param[in]  self   equation
   !> @param[in]  fin    input field
   !> @param[in]  fout   output field
   !> @param[out] error  not set up or shape mismatch
   subroutine equation_check_fields(self, fin, fout, error)
      !> Equation
      class(uv3d_equation_type), intent(in) :: self
      !> Input and output fields
      real(wp), intent(in) :: fin(:, :), fout(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%ready()) then
         call fatal_error(error, "3D UV equation is not set up")
      else if (any(shape(fin) /= [self%ngrid, self%nv]) .or. any(shape(fout) /= shape(fin))) then
         call fatal_error(error, "3D UV fields must have shape (ngrid, nv)")
      end if
   end subroutine equation_check_fields

   !> Build and prepare the transform engine of the owned grid unless present
   !>
   !> @param[inout] self   equation with a grid
   !> @param[out]   error  engine failure
   subroutine equation_ensure_trafo(self, error)
      !> Equation
      class(uv3d_equation_type), intent(inout), target :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (allocated(self%trafo)) return
      call self%grid%new_trafo(self%trafo, error)
      if (allocated(error)) return
      call self%trafo%prepare(self%nv, error)
      if (allocated(error)) then
         call self%trafo%destroy()
         deallocate (self%trafo)
      end if
   end subroutine equation_ensure_trafo

   !> Copy an equation without its transform engine
   !>
   !> @param[inout] lhs  copy; its own engine is released
   !> @param[in]    rhs  source
   subroutine equation_assign(lhs, rhs)
      !> Copy
      class(uv3d_equation_type), intent(inout) :: lhs
      !> Source
      class(uv3d_equation_type), intent(in) :: rhs
      call release(lhs)
      lhs%nthreads = rhs%nthreads
      if (.not. rhs%ready()) return
      allocate (lhs%grid, source=rhs%grid)
      lhs%response = rhs%response
      lhs%uk_lr = rhs%uk_lr
      lhs%r_work = rhs%r_work
      lhs%ck_work = rhs%ck_work
      lhs%tk_work = rhs%tk_work
      lhs%beta = rhs%beta
      lhs%ngrid = rhs%ngrid
      lhs%nk = rhs%nk
      lhs%nv = rhs%nv
   end subroutine equation_assign

   !> Release the engine and the setup
   !>
   !> @param[inout] self  equation
   subroutine release(self)
      !> Equation
      class(uv3d_equation_type), intent(inout) :: self
      if (allocated(self%trafo)) then
         call self%trafo%destroy()
         deallocate (self%trafo)
      end if
      if (allocated(self%grid)) deallocate (self%grid)
      if (allocated(self%uk_lr)) deallocate (self%uk_lr)
      if (allocated(self%r_work)) deallocate (self%r_work)
      if (allocated(self%ck_work)) deallocate (self%ck_work)
      if (allocated(self%tk_work)) deallocate (self%tk_work)
      self%nv = 0
      self%ngrid = 0
      self%nk = 0
   end subroutine release

   !> Release the engine and its plans
   !>
   !> @param[inout] self  equation being destroyed
   subroutine destroy_equation(self)
      !> Equation being destroyed
      type(uv3d_equation_type), intent(inout) :: self
      call release(self)
   end subroutine destroy_equation

end module moist_model_moz_3d_equation
