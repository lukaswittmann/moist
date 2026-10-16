!> Closure objects for molecular Ornstein-Zernike models
!>
!> Wraps closure values, derivatives and local free-energy integrands from kernel.f90,
!> and the partial derivatives of the integrands for response calculations
!> Built-in HNC, KH and PSE2-4 values and derivatives accept full u and t = h - c
!> or matching Ng short-range u_sr and t_s, returning c or c_s respectively
!> At full coupling: u_sr = u - u_lr, c_s = c + beta*u_lr, t_s = t - beta*u_lr
!> The caller owns the split and the lambda-stage convention; kernels scale the supplied u
!> Local energy integrands and their derivatives always require physical h and c
!> Arrays have shape (npts, npair); the caller handles splitting and grid layout
!>
!> TODO: Add closures with bridge functions evaluated from physical quantities
module moist_model_moz_closure
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: resolve_num_threads
!$ use omp_lib, only: omp_in_parallel
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use moist_model_moz_closure_kernel, only: &
      & closure_hnc_c, closure_hnc_phi, &
      & closure_kh_c, closure_kh_phi, &
      & closure_pse2_c, closure_pse2_phi, &
      & closure_pse3_c, closure_pse3_phi, &
      & closure_pse4_c, closure_pse4_phi
   use moist_model_moz_closure_kernel, only: closure_du, closure_phi_dc, &
      & closure_hnc_phi_dh, closure_kh_phi_dh, closure_pse2_phi_dh, closure_pse3_phi_dh, closure_pse4_phi_dh
   use moist_model_moz_closure_kernel, only: closure_dc_sweep, closure_d2c_sweep, closure_c_dc_sweep
   implicit none(type, external)
   private
   public :: moist_model_moz_closure_type
   public :: moist_model_moz_closure_hnc_type, moist_model_moz_closure_kh_type
   public :: moist_model_moz_closure_pse_type
   public :: new_closure

   !> Minimum number of values for parallel closure evaluation
   integer, parameter :: closure_parallel_threshold = 65536
   !> Points per block, allowing a single pair column to use several threads
   integer, parameter :: closure_sweep_block = 4096
   !> Lower bound on g = 1 + h; allows hard-core underflow and roundoff near zero
   real(wp), parameter :: g_floor = -1.0e-12_wp

   !> Local closure relation and its pointwise energy integrand
   type, abstract :: moist_model_moz_closure_type
      !> OpenMP team size; models should set ctx%get_num_threads(), zero takes omp_get_max_threads
      integer :: nthreads = 0
   contains
      !> OpenMP team size: one for small sweeps or inside an active parallel region
      procedure :: team_size => closure_team_size
      !> Evaluate c or c_s from the matching potential and indirect correlation
      procedure(closure_array_operation), deferred :: apply
      !> Evaluate dc/dt at fixed potential, beta and lambda
      procedure(closure_array_operation), deferred :: derivative
      !> Evaluate c and dc/dt together; HNC, KH and PSE share one exponential
      procedure :: apply_and_derivative => default_apply_and_derivative
      !> Evaluate d2c/dt2; optional for custom closures
      procedure :: derivative2 => unsupported_second_derivative
      !> Evaluate the dimensionless local energy integrand Phi(h, c)
      procedure(closure_point_integrand), deferred :: free_energy_integrand
      !> Pure local energy integrand without input checks, for parallel reductions
      procedure(closure_local_integrand), deferred :: local_integrand
      !> Local energy integrands over (npts, npair) arrays with input checks
      procedure :: free_energy_integrands => closure_free_energy_integrands
      !> Evaluate dc/du at fixed supplied t and split; optional for custom closures
      procedure :: potential_derivative => unsupported_potential_derivative
      !> Partial derivatives of the local energy integrand at independent h and c
      procedure :: free_energy_derivatives => default_energy_derivatives
      !> Optional pure derivatives; the default returns NaNs, checked calls report errors
      procedure :: local_energy_derivatives => unsupported_local_energy_derivatives
      !> Local energy derivatives over (npts, npair) arrays with input checks
      procedure :: free_energy_derivative_fields => closure_free_energy_derivative_fields
   end type moist_model_moz_closure_type

   !> Array and pointwise operations shared by MOZ models
   abstract interface
      !> Evaluate a closure value or its derivative with respect to t
      !>
      !> @param[in] self Closure
      !> @param[in] u Full u or Ng short-range u_sr (npts, npair), energy units
      !> @param[in] t Matching full t = h - c or Ng t_s (npts, npair), dimensionless
      !> @param[in] beta Inverse thermal energy, reciprocal units of u
      !> @param[in] lambda Potential scaling parameter, between zero and one
      !> @param[out] values c or dc/dt (npts, npair), defined only on success
      !> @param[out] error Invalid input or non-finite result
      subroutine closure_array_operation(self, u, t, beta, lambda, values, error)
         import :: moist_model_moz_closure_type, wp, error_type
         implicit none(type, external)
         !> Closure
         class(moist_model_moz_closure_type), intent(in) :: self
         !> Full u or Ng short-range u_sr (npts, npair), energy units
         real(wp), intent(in) :: u(:, :)
         !> Matching full t = h - c or Ng t_s (npts, npair), dimensionless
         real(wp), intent(in) :: t(:, :)
         !> Inverse thermal energy, reciprocal units of u
         real(wp), intent(in) :: beta
         !> Potential scaling parameter, between zero and one
         real(wp), intent(in) :: lambda
         !> c or dc/dt (npts, npair), defined only on success
         real(wp), intent(out) :: values(:, :)
         !> Invalid input or non-finite result
         type(error_type), allocatable, intent(out) :: error
      end subroutine closure_array_operation

      !> Evaluate a local energy integrand without density or grid weights
      !>
      !> @param[in] self Closure
      !> @param[in] h Physical total correlation, dimensionless
      !> @param[in] c Physical direct correlation, dimensionless
      !> @param[out] phi Dimensionless integrand, defined only on success
      !> @param[out] error Invalid input or non-finite result
      subroutine closure_point_integrand(self, h, c, phi, error)
         import :: moist_model_moz_closure_type, wp, error_type
         implicit none(type, external)
         !> Closure
         class(moist_model_moz_closure_type), intent(in) :: self
         !> Physical total correlation, dimensionless
         real(wp), intent(in) :: h
         !> Physical direct correlation, dimensionless
         real(wp), intent(in) :: c
         !> Dimensionless integrand, defined only on success
         real(wp), intent(out) :: phi
         !> Invalid input or non-finite result
         type(error_type), allocatable, intent(out) :: error
      end subroutine closure_point_integrand

      !> Evaluate the local energy integrand without checks; callers validate h and c
      !>
      !> @param[in] self Closure
      !> @param[in] h Physical total correlation, dimensionless; 1 + h >= g_floor
      !> @param[in] c Physical direct correlation, dimensionless
      !> @return Dimensionless integrand
      pure function closure_local_integrand(self, h, c) result(phi)
         import :: moist_model_moz_closure_type, wp
         implicit none(type, external)
         !> Closure
         class(moist_model_moz_closure_type), intent(in) :: self
         !> Physical total correlation, dimensionless
         real(wp), intent(in) :: h
         !> Physical direct correlation, dimensionless
         real(wp), intent(in) :: c
         !> Dimensionless integrand
         real(wp) :: phi
      end function closure_local_integrand

      !> Evaluate c or c_s over a block of matching full or Ng inputs
      !>
      !> @param[in] u Full u or Ng short-range u_sr, energy units
      !> @param[in] t Matching full t = h - c or Ng t_s, dimensionless
      !> @param[in] beta Inverse thermal energy, reciprocal units of u
      !> @param[in] lambda Potential scaling parameter, between zero and one
      !> @param[out] values Direct correlation
      pure subroutine closure_c_sweep(u, t, beta, lambda, values)
         import :: wp
         implicit none(type, external)
         !> Full u or Ng short-range u_sr, energy units
         real(wp), intent(in) :: u(:)
         !> Matching full t = h - c or Ng t_s, dimensionless
         real(wp), intent(in) :: t(:)
         !> Inverse thermal energy and potential scaling
         real(wp), intent(in) :: beta, lambda
         !> Direct correlation
         real(wp), intent(out) :: values(:)
      end subroutine closure_c_sweep
   end interface

   !> Hypernetted-chain closure
   type, extends(moist_model_moz_closure_type) :: moist_model_moz_closure_hnc_type
   contains
      !> Evaluate the direct correlation
      procedure :: apply => hnc_apply
      !> Evaluate dc/dt
      procedure :: derivative => hnc_derivative
      !> Evaluate c and dc/dt in one fused sweep
      procedure :: apply_and_derivative => hnc_apply_and_derivative
      !> Evaluate d2c/dt2
      procedure :: derivative2 => hnc_derivative2
      !> Evaluate the local energy integrand
      procedure :: free_energy_integrand => hnc_free_energy_integrand
      !> Pure local energy integrand
      procedure :: local_integrand => hnc_local_integrand
      !> Potential and local energy derivatives used by response calculations
      procedure :: potential_derivative => hnc_potential_derivative
      procedure :: free_energy_derivatives => hnc_energy_derivatives
      procedure :: local_energy_derivatives => hnc_local_derivatives
      procedure :: free_energy_derivative_fields => hnc_energy_derivative_fields
   end type moist_model_moz_closure_hnc_type

   !> Kovalenko-Hirata closure
   type, extends(moist_model_moz_closure_type) :: moist_model_moz_closure_kh_type
   contains
      !> Evaluate the direct correlation
      procedure :: apply => kh_apply
      !> Evaluate dc/dt
      procedure :: derivative => kh_derivative
      !> Evaluate c and dc/dt in one fused sweep
      procedure :: apply_and_derivative => kh_apply_and_derivative
      !> Evaluate d2c/dt2
      procedure :: derivative2 => kh_derivative2
      !> Evaluate the local energy integrand
      procedure :: free_energy_integrand => kh_free_energy_integrand
      !> Pure local energy integrand
      procedure :: local_integrand => kh_local_integrand
      !> Potential and local energy derivatives used by response calculations
      procedure :: potential_derivative => kh_potential_derivative
      procedure :: free_energy_derivatives => kh_energy_derivatives
      procedure :: local_energy_derivatives => kh_local_derivatives
      procedure :: free_energy_derivative_fields => kh_energy_derivative_fields
   end type moist_model_moz_closure_kh_type

   !> Partial series expansion closure of order two, three or four
   type, extends(moist_model_moz_closure_type) :: moist_model_moz_closure_pse_type
      !> Polynomial expansion order; valid values are 2, 3 and 4
      integer :: order = 2
   contains
      !> Evaluate the direct correlation
      procedure :: apply => pse_apply
      !> Evaluate dc/dt
      procedure :: derivative => pse_derivative
      !> Evaluate c and dc/dt in one fused sweep
      procedure :: apply_and_derivative => pse_apply_and_derivative
      !> Evaluate d2c/dt2
      procedure :: derivative2 => pse_derivative2
      !> Evaluate the local energy integrand
      procedure :: free_energy_integrand => pse_free_energy_integrand
      !> Pure local energy integrand
      procedure :: local_integrand => pse_local_integrand
      !> Potential and local energy derivatives used by response calculations
      procedure :: potential_derivative => pse_potential_derivative
      procedure :: free_energy_derivatives => pse_energy_derivatives
      procedure :: local_energy_derivatives => pse_local_derivatives
      procedure :: free_energy_derivative_fields => pse_energy_derivative_fields
   end type moist_model_moz_closure_pse_type

contains

   !> Construct an HNC, KH or PSE2-4 closure
   !>
   !> @param[out] closure Allocated closure; unallocated on failure
   !> @param[in] name Closure identifier: HNC, KH or PSE2-4
   !> @param[out] error Unknown closure identifier or negative thread count
   !> @param[in] nthreads OpenMP team size; absent or zero takes omp_get_max_threads
   subroutine new_closure(closure, name, error, nthreads)
      !> Allocated closure; unallocated on failure
      class(moist_model_moz_closure_type), allocatable, intent(out) :: closure
      !> Closure identifier: HNC, KH or PSE2-4
      character(len=*), intent(in) :: name
      !> Unknown closure identifier or negative thread count
      type(error_type), allocatable, intent(out) :: error
      !> OpenMP team size, e.g. ctx%get_num_threads(); absent or zero takes omp_get_max_threads
      integer, intent(in), optional :: nthreads

      select case (trim(adjustl(name)))
      case ("HNC", "hnc")
         allocate (moist_model_moz_closure_hnc_type :: closure)
      case ("KH", "kh")
         allocate (moist_model_moz_closure_kh_type :: closure)
      case ("PSE2", "pse2", "PSE3", "pse3", "PSE4", "pse4")
         allocate (moist_model_moz_closure_pse_type :: closure)
         select type (closure)
         type is (moist_model_moz_closure_pse_type)
            select case (trim(adjustl(name)))
            case ("PSE2", "pse2")
               closure%order = 2
            case ("PSE3", "pse3")
               closure%order = 3
            case ("PSE4", "pse4")
               closure%order = 4
            case default
               call fatal_error(error, "Unknown PSE closure: "//trim(name))
            end select
         end select
      case default
         call fatal_error(error, "Unknown MOZ closure: "//trim(name))
      end select
      if (.not. allocated(closure) .or. .not. present(nthreads)) return
      if (nthreads < 0) then
         deallocate (closure)
         call fatal_error(error, "Closure thread count must be nonnegative")
         return
      end if
      closure%nthreads = nthreads
   end subroutine new_closure

   !> OpenMP team size for n values from the closure's thread count
   !>
   !> @param[in] self Closure
   !> @param[in] n Number of values in the sweep
   function closure_team_size(self, n) result(nt)
      !> Closure
      class(moist_model_moz_closure_type), intent(in) :: self
      !> Number of values in the sweep
      integer, intent(in) :: n
      !> Team size (>= 1)
      integer :: nt
      nt = 1
      ! Block producers already distribute work; closures must not create nested teams.
!$    if (omp_in_parallel()) return
      if (n >= closure_parallel_threshold) nt = resolve_num_threads(self%nthreads)
   end function closure_team_size

   !> Evaluate a closure kernel and check finite values within each block
   !>
   !> @param[in] sweep Closure block kernel
   !> @param[in] u Full u or Ng short-range u_sr (npts, npair), energy units
   !> @param[in] t Matching full t = h - c or Ng t_s (npts, npair), dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of u
   !> @param[in] lambda Potential scaling parameter, between zero and one
   !> @param[out] values Direct correlation (npts, npair)
   !> @param[in] nt OpenMP team size
   !> @param[out] error Non-finite input or result
   subroutine sweep_closure_c(sweep, u, t, beta, lambda, values, nt, error)
      !> Closure block kernel
      procedure(closure_c_sweep) :: sweep
      !> Full u or Ng short-range u_sr (npts, npair), energy units
      real(wp), intent(in) :: u(:, :)
      !> Matching full t = h - c or Ng t_s (npts, npair), dimensionless
      real(wp), intent(in) :: t(:, :)
      !> Inverse thermal energy and potential scaling
      real(wp), intent(in) :: beta, lambda
      !> Direct correlation (npts, npair)
      real(wp), intent(out) :: values(:, :)
      !> OpenMP team size
      integer, intent(in) :: nt
      !> Non-finite input or result
      type(error_type), allocatable, intent(out) :: error
      integer :: j, ib, i0, i1, nblk
      logical :: bad_input, bad_result

      if (size(t) == 0) return
      nblk = (size(t, 1) - 1)/closure_sweep_block + 1
      bad_input = .false.
      bad_result = .false.
      ! Keep the checks close to the kernel while the block is in cache.
      ! Only flags are reduced; errors are allocated after the parallel region.
      !$omp parallel do collapse(2) schedule(static) private(i0, i1) &
      !$omp reduction(.or.:bad_input, bad_result) &
      !$omp num_threads(nt)
      do j = 1, size(t, 2)
         do ib = 1, nblk
            i0 = (ib - 1)*closure_sweep_block + 1
            i1 = min(ib*closure_sweep_block, size(t, 1))
            if (.not. all(ieee_is_finite(u(i0:i1, j))) .or. &
               & .not. all(ieee_is_finite(t(i0:i1, j)))) then
               bad_input = .true.
               cycle
            end if
            call sweep(u(i0:i1, j), t(i0:i1, j), beta, lambda, values(i0:i1, j))
            if (.not. all(ieee_is_finite(values(i0:i1, j)))) bad_result = .true.
         end do
      end do
      !$omp end parallel do
      if (bad_input) then
         call fatal_error(error, "Closure potential and correlation must be finite")
      else if (bad_result) then
         call fatal_error(error, "Closure kernel returned a non-finite result")
      end if
   end subroutine sweep_closure_c

   !> Evaluate dc/dt, and optionally c, over blocks with finite checks
   !>
   !> @param[in] kind Closure kind: 0 HNC, 1 KH, 2-4 PSE order
   !> @param[in] u Full u or Ng short-range u_sr (npts, npair), energy units
   !> @param[in] t Matching full t = h - c or Ng t_s (npts, npair), dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of u
   !> @param[in] lambda Potential scaling parameter, between zero and one
   !> @param[out] derivs Derivative dc/dt (npts, npair)
   !> @param[in] nt OpenMP team size
   !> @param[out] error Non-finite input or result
   !> @param[out] values Direct correlation (npts, npair), from the same exponential
   !> @param[in] second Evaluate d2c/dt2 instead of dc/dt; not combined with values
   subroutine sweep_closure_dc(kind, u, t, beta, lambda, derivs, nt, error, values, second)
      !> Closure kind: 0 HNC, 1 KH, 2-4 PSE order
      integer, intent(in) :: kind
      !> Full u or Ng short-range u_sr (npts, npair), energy units
      real(wp), intent(in) :: u(:, :)
      !> Matching full t = h - c or Ng t_s (npts, npair), dimensionless
      real(wp), intent(in) :: t(:, :)
      !> Inverse thermal energy and potential scaling
      real(wp), intent(in) :: beta, lambda
      !> Derivative dc/dt (npts, npair)
      real(wp), intent(out) :: derivs(:, :)
      !> OpenMP team size
      integer, intent(in) :: nt
      !> Non-finite input or result
      type(error_type), allocatable, intent(out) :: error
      !> Direct correlation (npts, npair)
      real(wp), intent(out), optional :: values(:, :)
      !> Evaluate d2c/dt2 instead of dc/dt
      logical, intent(in), optional :: second
      integer :: j, ib, i0, i1, nblk
      logical :: bad_input, bad_result, fused, second_order

      if (size(t) == 0) return
      fused = present(values)
      second_order = .false.
      if (present(second)) second_order = second
      nblk = (size(t, 1) - 1)/closure_sweep_block + 1
      bad_input = .false.
      bad_result = .false.
      !$omp parallel do collapse(2) schedule(static) private(i0, i1) &
      !$omp reduction(.or.:bad_input, bad_result) &
      !$omp num_threads(nt)
      do j = 1, size(t, 2)
         do ib = 1, nblk
            i0 = (ib - 1)*closure_sweep_block + 1
            i1 = min(ib*closure_sweep_block, size(t, 1))
            if (.not. all(ieee_is_finite(u(i0:i1, j))) .or. &
               & .not. all(ieee_is_finite(t(i0:i1, j)))) then
               bad_input = .true.
               cycle
            end if
            if (fused) then
               call closure_c_dc_sweep(u(i0:i1, j), t(i0:i1, j), beta, lambda, kind, &
                  & values(i0:i1, j), derivs(i0:i1, j))
               if (.not. all(ieee_is_finite(values(i0:i1, j)))) bad_result = .true.
            else if (second_order) then
               call closure_d2c_sweep(u(i0:i1, j), t(i0:i1, j), beta, lambda, kind, derivs(i0:i1, j))
            else
               call closure_dc_sweep(u(i0:i1, j), t(i0:i1, j), beta, lambda, kind, derivs(i0:i1, j))
            end if
            if (.not. all(ieee_is_finite(derivs(i0:i1, j)))) bad_result = .true.
         end do
      end do
      !$omp end parallel do
      if (bad_input) then
         call fatal_error(error, "Closure potential and correlation must be finite")
      else if (bad_result) then
         call fatal_error(error, "Closure kernel returned a non-finite result")
      end if
   end subroutine sweep_closure_dc

   !> Evaluate c and dc/dt with two calls; closures without a fused kernel use this
   !>
   !> @param[in] self Closure
   !> @param[in] u Full u or Ng short-range u_sr (npts, npair), energy units
   !> @param[in] t Matching full t = h - c or Ng t_s (npts, npair), dimensionless
   !> @param[in] beta Inverse thermal energy, reciprocal units of u
   !> @param[in] lambda Potential scaling parameter, between zero and one
   !> @param[out] values Direct correlation (npts, npair), defined only on success
   !> @param[out] derivs Derivative dc/dt (npts, npair), defined only on success
   !> @param[out] error Invalid input or non-finite result
   subroutine default_apply_and_derivative(self, u, t, beta, lambda, values, derivs, error)
      !> Closure
      class(moist_model_moz_closure_type), intent(in) :: self
      !> Full u or Ng short-range u_sr (npts, npair), energy units
      real(wp), intent(in) :: u(:, :)
      !> Matching full t = h - c or Ng t_s (npts, npair), dimensionless
      real(wp), intent(in) :: t(:, :)
      !> Inverse thermal energy and potential scaling
      real(wp), intent(in) :: beta, lambda
      !> Direct correlation (npts, npair), defined only on success
      real(wp), intent(out) :: values(:, :)
      !> Derivative dc/dt (npts, npair), defined only on success
      real(wp), intent(out) :: derivs(:, :)
      !> Invalid input or non-finite result
      type(error_type), allocatable, intent(out) :: error

      call self%apply(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      call self%derivative(u, t, beta, lambda, derivs, error)
   end subroutine default_apply_and_derivative

   !> HNC direct correlation and derivative in one fused sweep
   subroutine hnc_apply_and_derivative(self, u, t, beta, lambda, values, derivs, error)
      class(moist_model_moz_closure_hnc_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta, lambda
      real(wp), intent(out) :: values(:, :)
      real(wp), intent(out) :: derivs(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      if (any(shape(derivs) /= shape(values))) then
         call fatal_error(error, "Closure arrays must have matching shapes")
         return
      end if
      call sweep_closure_dc(0, u, t, beta, lambda, derivs, self%team_size(size(t)), error, values)
   end subroutine hnc_apply_and_derivative

   !> KH direct correlation and derivative in one fused sweep
   subroutine kh_apply_and_derivative(self, u, t, beta, lambda, values, derivs, error)
      class(moist_model_moz_closure_kh_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta, lambda
      real(wp), intent(out) :: values(:, :)
      real(wp), intent(out) :: derivs(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      if (any(shape(derivs) /= shape(values))) then
         call fatal_error(error, "Closure arrays must have matching shapes")
         return
      end if
      call sweep_closure_dc(1, u, t, beta, lambda, derivs, self%team_size(size(t)), error, values)
   end subroutine kh_apply_and_derivative

   !> PSE direct correlation and derivative in one fused sweep
   subroutine pse_apply_and_derivative(self, u, t, beta, lambda, values, derivs, error)
      class(moist_model_moz_closure_pse_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta, lambda
      real(wp), intent(out) :: values(:, :)
      real(wp), intent(out) :: derivs(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      if (any(shape(derivs) /= shape(values))) then
         call fatal_error(error, "Closure arrays must have matching shapes")
         return
      end if
      if (self%order < 2 .or. self%order > 4) then
         call fatal_error(error, "PSE order must be 2, 3 or 4")
         return
      end if
      call sweep_closure_dc(self%order, u, t, beta, lambda, derivs, self%team_size(size(t)), error, values)
   end subroutine pse_apply_and_derivative

   !> HNC direct correlation
   subroutine hnc_apply(self, u, t, beta, lambda, values, error)
      class(moist_model_moz_closure_hnc_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta
      real(wp), intent(in) :: lambda
      real(wp), intent(out) :: values(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      call sweep_closure_c(closure_hnc_c, u, t, beta, lambda, values, self%team_size(size(t)), error)
   end subroutine hnc_apply

   !> HNC derivative with respect to t
   subroutine hnc_derivative(self, u, t, beta, lambda, values, error)
      class(moist_model_moz_closure_hnc_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta
      real(wp), intent(in) :: lambda
      real(wp), intent(out) :: values(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      call sweep_closure_dc(0, u, t, beta, lambda, values, self%team_size(size(t)), error)
   end subroutine hnc_derivative

   !> HNC local energy integrand
   !>
   !> @param[in] self Closure
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   !> @param[out] phi Dimensionless integrand, defined only on success
   !> @param[out] error Invalid input or non-finite result
   subroutine hnc_free_energy_integrand(self, h, c, phi, error)
      !> Closure
      class(moist_model_moz_closure_hnc_type), intent(in) :: self
      !> Physical total correlation, dimensionless
      real(wp), intent(in) :: h
      !> Physical direct correlation, dimensionless
      real(wp), intent(in) :: c
      !> Dimensionless integrand, defined only on success
      real(wp), intent(out) :: phi
      !> Invalid input or non-finite result
      type(error_type), allocatable, intent(out) :: error

      call validate_integrand_inputs(h, c, error)
      if (allocated(error)) return
      phi = closure_hnc_phi(h, c)
      if (.not. ieee_is_finite(phi)) call fatal_error(error, "Closure energy integrand is non-finite")
   end subroutine hnc_free_energy_integrand

   !> KH direct correlation
   subroutine kh_apply(self, u, t, beta, lambda, values, error)
      class(moist_model_moz_closure_kh_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta
      real(wp), intent(in) :: lambda
      real(wp), intent(out) :: values(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      call sweep_closure_c(closure_kh_c, u, t, beta, lambda, values, self%team_size(size(t)), error)
   end subroutine kh_apply

   !> KH derivative with respect to t
   subroutine kh_derivative(self, u, t, beta, lambda, values, error)
      class(moist_model_moz_closure_kh_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta
      real(wp), intent(in) :: lambda
      real(wp), intent(out) :: values(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      call sweep_closure_dc(1, u, t, beta, lambda, values, self%team_size(size(t)), error)
   end subroutine kh_derivative

   !> KH local energy integrand
   !>
   !> @param[in] self Closure
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   !> @param[out] phi Dimensionless integrand, defined only on success
   !> @param[out] error Invalid input or non-finite result
   subroutine kh_free_energy_integrand(self, h, c, phi, error)
      !> Closure
      class(moist_model_moz_closure_kh_type), intent(in) :: self
      !> Physical total correlation, dimensionless
      real(wp), intent(in) :: h
      !> Physical direct correlation, dimensionless
      real(wp), intent(in) :: c
      !> Dimensionless integrand, defined only on success
      real(wp), intent(out) :: phi
      !> Invalid input or non-finite result
      type(error_type), allocatable, intent(out) :: error

      call validate_integrand_inputs(h, c, error)
      if (allocated(error)) return
      phi = closure_kh_phi(h, c)
      if (.not. ieee_is_finite(phi)) call fatal_error(error, "Closure energy integrand is non-finite")
   end subroutine kh_free_energy_integrand

   !> PSE direct correlation
   subroutine pse_apply(self, u, t, beta, lambda, values, error)
      class(moist_model_moz_closure_pse_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta
      real(wp), intent(in) :: lambda
      real(wp), intent(out) :: values(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      select case (self%order)
      case (2)
         call sweep_closure_c(closure_pse2_c, u, t, beta, lambda, values, self%team_size(size(t)), error)
      case (3)
         call sweep_closure_c(closure_pse3_c, u, t, beta, lambda, values, self%team_size(size(t)), error)
      case (4)
         call sweep_closure_c(closure_pse4_c, u, t, beta, lambda, values, self%team_size(size(t)), error)
      case default
         call validate_arrays(u, t, beta, lambda, values, error)
         if (allocated(error)) return
         call fatal_error(error, "PSE order must be 2, 3 or 4")
         return
      end select
   end subroutine pse_apply

   !> PSE derivative with respect to t
   subroutine pse_derivative(self, u, t, beta, lambda, values, error)
      class(moist_model_moz_closure_pse_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta
      real(wp), intent(in) :: lambda
      real(wp), intent(out) :: values(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      if (self%order < 2 .or. self%order > 4) then
         call fatal_error(error, "PSE order must be 2, 3 or 4")
         return
      end if
      call sweep_closure_dc(self%order, u, t, beta, lambda, values, self%team_size(size(t)), error)
   end subroutine pse_derivative

   !> PSE local energy integrand
   !>
   !> @param[in] self Closure
   !> @param[in] h Physical total correlation, dimensionless
   !> @param[in] c Physical direct correlation, dimensionless
   !> @param[out] phi Dimensionless integrand, defined only on success
   !> @param[out] error Invalid input or non-finite result
   subroutine pse_free_energy_integrand(self, h, c, phi, error)
      !> Closure
      class(moist_model_moz_closure_pse_type), intent(in) :: self
      !> Physical total correlation, dimensionless
      real(wp), intent(in) :: h
      !> Physical direct correlation, dimensionless
      real(wp), intent(in) :: c
      !> Dimensionless integrand, defined only on success
      real(wp), intent(out) :: phi
      !> Invalid input or non-finite result
      type(error_type), allocatable, intent(out) :: error

      call validate_integrand_inputs(h, c, error)
      if (allocated(error)) return
      select case (self%order)
      case (2)
         phi = closure_pse2_phi(h, c)
      case (3)
         phi = closure_pse3_phi(h, c)
      case (4)
         phi = closure_pse4_phi(h, c)
      case default
         call fatal_error(error, "PSE order must be 2, 3 or 4")
         return
      end select
      if (.not. ieee_is_finite(phi)) call fatal_error(error, "Closure energy integrand is non-finite")
   end subroutine pse_free_energy_integrand

   !> Check array shapes and scalar parameters without scanning the fields
   !>
   !> @param[in] u Full u or Ng short-range u_sr (npts, npair), energy units
   !> @param[in] t Matching full t or Ng t_s (npts, npair), dimensionless
   !> @param[in] beta Positive inverse thermal energy
   !> @param[in] lambda Potential scaling in [0, 1]
   !> @param[out] values Output array (npts, npair); only its shape is checked
   !> @param[out] error Invalid shape or scalar parameter
   subroutine validate_array_parameters(u, t, beta, lambda, values, error)
      !> Full u or Ng short-range u_sr (npts, npair), energy units
      real(wp), intent(in) :: u(:, :)
      !> Matching full t or Ng t_s (npts, npair), dimensionless
      real(wp), intent(in) :: t(:, :)
      !> Positive inverse thermal energy
      real(wp), intent(in) :: beta
      !> Potential scaling in [0, 1]
      real(wp), intent(in) :: lambda
      !> Output array (npts, npair); only its shape is checked
      real(wp), intent(out) :: values(:, :)
      !> Invalid shape or scalar parameter
      type(error_type), allocatable, intent(out) :: error

      if (any(shape(u) /= shape(t)) .or. any(shape(u) /= shape(values))) then
         call fatal_error(error, "Closure arrays must have matching shapes")
      else if (.not. ieee_is_finite(beta) .or. .not. ieee_is_finite(lambda)) then
         call fatal_error(error, "Closure beta and lambda must be finite")
      else if (beta <= 0.0_wp .or. lambda < 0.0_wp .or. lambda > 1.0_wp) then
         call fatal_error(error, "Closure requires beta > 0 and lambda in [0, 1]")
      end if
   end subroutine validate_array_parameters

   !> Check shapes, scalar parameters and finite inputs for derivative evaluation
   !>
   !> @param[in] u Full u or Ng short-range u_sr (npts, npair), energy units
   !> @param[in] t Matching full t or Ng t_s (npts, npair), dimensionless
   !> @param[in] beta Positive inverse thermal energy
   !> @param[in] lambda Potential scaling in [0, 1]
   !> @param[out] values Output array (npts, npair); only its shape is checked
   !> @param[out] error Invalid shape, domain or non-finite input
   subroutine validate_arrays(u, t, beta, lambda, values, error)
      real(wp), intent(in) :: u(:, :), t(:, :), beta, lambda
      real(wp), intent(out) :: values(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      if (.not. all(ieee_is_finite(u)) .or. .not. all(ieee_is_finite(t))) then
         call fatal_error(error, "Closure potential and correlation must be finite")
      end if
   end subroutine validate_arrays

   !> Report a non-finite kernel output, including exponential overflow
   !>
   !> @param[in] values Computed closure values (npts, npair)
   !> @param[out] error Non-finite result
   subroutine validate_result(values, error)
      !> Computed closure values (npts, npair)
      real(wp), intent(in) :: values(:, :)
      !> Non-finite result
      type(error_type), allocatable, intent(out) :: error

      if (.not. all(ieee_is_finite(values))) call fatal_error(error, "Closure kernel returned a non-finite result")
   end subroutine validate_result

   !> Local energy integrands over arrays, validated and evaluated in parallel
   !>
   !> @param[in] self Closure
   !> @param[in] h Physical total correlation (npts, npair); 1 + h >= g_floor
   !> @param[in] c Physical direct correlation (npts, npair)
   !> @param[out] phi Dimensionless integrands (npts, npair), defined only on success
   !> @param[out] error Shape mismatch, invalid input or non-finite result
   subroutine closure_free_energy_integrands(self, h, c, phi, error)
      !> Closure
      class(moist_model_moz_closure_type), intent(in) :: self
      !> Physical total correlation (npts, npair)
      real(wp), intent(in) :: h(:, :)
      !> Physical direct correlation (npts, npair)
      real(wp), intent(in) :: c(:, :)
      !> Dimensionless integrands (npts, npair)
      real(wp), intent(out) :: phi(:, :)
      !> Shape mismatch, invalid input or non-finite result
      type(error_type), allocatable, intent(out) :: error
      integer :: i, j, nt
      logical :: bad_input, bad_result

      if (any(shape(h) /= shape(c)) .or. any(shape(h) /= shape(phi))) then
         call fatal_error(error, "Closure energy arrays must have matching shapes")
         return
      end if
      bad_input = .false.
      bad_result = .false.
      nt = self%team_size(size(h))
      !$omp parallel do collapse(2) schedule(static) reduction(.or.:bad_input, bad_result) &
      !$omp num_threads(nt)
      do j = 1, size(h, 2)
         do i = 1, size(h, 1)
            if (.not. ieee_is_finite(h(i, j)) .or. .not. ieee_is_finite(c(i, j)) &
               & .or. 1.0_wp + h(i, j) < g_floor) then
               bad_input = .true.
               cycle
            end if
            phi(i, j) = self%local_integrand(h(i, j), c(i, j))
            if (.not. ieee_is_finite(phi(i, j))) bad_result = .true.
         end do
      end do
      !$omp end parallel do
      if (bad_input) then
         call fatal_error(error, "Closure energy correlations must be finite with 1 + h >= g_floor")
      else if (bad_result) then
         call fatal_error(error, "Closure energy integrand is non-finite")
      end if
   end subroutine closure_free_energy_integrands

   !> HNC local energy integrand without checks
   pure function hnc_local_integrand(self, h, c) result(phi)
      !> Closure
      class(moist_model_moz_closure_hnc_type), intent(in) :: self
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Dimensionless integrand
      real(wp) :: phi
      phi = closure_hnc_phi(h, c)
   end function hnc_local_integrand

   !> KH local energy integrand without checks
   pure function kh_local_integrand(self, h, c) result(phi)
      !> Closure
      class(moist_model_moz_closure_kh_type), intent(in) :: self
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Dimensionless integrand
      real(wp) :: phi
      phi = closure_kh_phi(h, c)
   end function kh_local_integrand

   !> PSE local energy integrand without checks; NaN for an invalid order
   pure function pse_local_integrand(self, h, c) result(phi)
      !> Closure
      class(moist_model_moz_closure_pse_type), intent(in) :: self
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Dimensionless integrand
      real(wp) :: phi
      select case (self%order)
      case (2)
         phi = closure_pse2_phi(h, c)
      case (3)
         phi = closure_pse3_phi(h, c)
      case (4)
         phi = closure_pse4_phi(h, c)
      case default
         phi = ieee_value(0.0_wp, ieee_quiet_nan)
      end select
   end function pse_local_integrand

   !> Checked derivative fields for custom closures through their scalar operation
   !>
   !> Energy-only closures inherit an unsupported error. A custom closure may
   !> override the checked scalar operation or supply the optional pure hook.
   subroutine closure_free_energy_derivative_fields(self, h, c, dh, dc, error)
      class(moist_model_moz_closure_type), intent(in) :: self
      !> Physical correlations and unweighted integrand derivatives
      real(wp), intent(in) :: h(:, :), c(:, :)
      real(wp), intent(out) :: dh(:, :), dc(:, :)
      !> Invalid input, missing derivatives or non-finite result
      type(error_type), allocatable, intent(out) :: error
      integer :: i, j

      if (any(shape(h) /= shape(c)) .or. any(shape(h) /= shape(dh)) .or. any(shape(h) /= shape(dc))) then
         call fatal_error(error, "Closure energy derivative arrays must have matching shapes")
         return
      end if
      do j = 1, size(h, 2)
         do i = 1, size(h, 1)
            call validate_integrand_inputs(h(i, j), c(i, j), error)
            if (allocated(error)) return
            call self%free_energy_derivatives(h(i, j), c(i, j), dh(i, j), dc(i, j), error)
            if (allocated(error)) return
            if (.not. ieee_is_finite(dh(i, j)) .or. .not. ieee_is_finite(dc(i, j))) then
               call fatal_error(error, "Closure energy derivative is non-finite")
               return
            end if
         end do
      end do
   end subroutine closure_free_energy_derivative_fields

   !> Parallel pure derivative sweep retained for the built-in closures
   !>
   !> dh = dPhi/dh at fixed c and dc = dPhi/dc at fixed h, pointwise, without
   !> density or grid weights; the response contracts them with its own weights
   !>
   !> @param[in]  self   Closure
   !> @param[in]  h      Physical total correlation (npts, npair)
   !> @param[in]  c      Physical direct correlation (npts, npair)
   !> @param[out] dh     dPhi/dh (npts, npair), defined only on success
   !> @param[out] dc     dPhi/dc (npts, npair), defined only on success
   !> @param[out] error  Shape mismatch, invalid input or non-finite result
   subroutine local_energy_derivative_fields(self, h, c, dh, dc, error)
      !> Closure
      class(moist_model_moz_closure_type), intent(in) :: self
      !> Physical total and direct correlations (npts, npair)
      real(wp), intent(in) :: h(:, :), c(:, :)
      !> Integrand partial derivatives (npts, npair)
      real(wp), intent(out) :: dh(:, :), dc(:, :)
      !> Shape mismatch, invalid input or non-finite result
      type(error_type), allocatable, intent(out) :: error
      integer :: i, j, nt
      logical :: bad_input, bad_result

      if (any(shape(h) /= shape(c)) .or. any(shape(h) /= shape(dh)) .or. any(shape(h) /= shape(dc))) then
         call fatal_error(error, "Closure energy derivative arrays must have matching shapes")
         return
      end if
      bad_input = .false.
      bad_result = .false.
      nt = self%team_size(size(h))
      !$omp parallel do collapse(2) schedule(static) reduction(.or.:bad_input, bad_result) &
      !$omp num_threads(nt)
      do j = 1, size(h, 2)
         do i = 1, size(h, 1)
            if (.not. ieee_is_finite(h(i, j)) .or. .not. ieee_is_finite(c(i, j)) &
               & .or. 1.0_wp + h(i, j) < g_floor) then
               bad_input = .true.
               cycle
            end if
            call self%local_energy_derivatives(h(i, j), c(i, j), dh(i, j), dc(i, j))
            if (.not. ieee_is_finite(dh(i, j)) .or. .not. ieee_is_finite(dc(i, j))) bad_result = .true.
         end do
      end do
      !$omp end parallel do
      if (bad_input) then
         call fatal_error(error, "Closure energy correlations must be finite with 1 + h >= g_floor")
      else if (bad_result) then
         call fatal_error(error, "Closure energy derivative is non-finite")
      end if
   end subroutine local_energy_derivative_fields

   !> HNC array derivatives through the shared parallel sweep
   subroutine hnc_energy_derivative_fields(self, h, c, dh, dc, error)
      class(moist_model_moz_closure_hnc_type), intent(in) :: self
      real(wp), intent(in) :: h(:, :), c(:, :)
      real(wp), intent(out) :: dh(:, :), dc(:, :)
      type(error_type), allocatable, intent(out) :: error
      call local_energy_derivative_fields(self, h, c, dh, dc, error)
   end subroutine hnc_energy_derivative_fields

   !> KH array derivatives through the shared parallel sweep
   subroutine kh_energy_derivative_fields(self, h, c, dh, dc, error)
      class(moist_model_moz_closure_kh_type), intent(in) :: self
      real(wp), intent(in) :: h(:, :), c(:, :)
      real(wp), intent(out) :: dh(:, :), dc(:, :)
      type(error_type), allocatable, intent(out) :: error
      call local_energy_derivative_fields(self, h, c, dh, dc, error)
   end subroutine kh_energy_derivative_fields

   !> PSE array derivatives through the shared parallel sweep
   subroutine pse_energy_derivative_fields(self, h, c, dh, dc, error)
      class(moist_model_moz_closure_pse_type), intent(in) :: self
      real(wp), intent(in) :: h(:, :), c(:, :)
      real(wp), intent(out) :: dh(:, :), dc(:, :)
      type(error_type), allocatable, intent(out) :: error
      call local_energy_derivative_fields(self, h, c, dh, dc, error)
   end subroutine pse_energy_derivative_fields

   !> HNC local energy derivatives without checks
   pure subroutine hnc_local_derivatives(self, h, c, dh, dc)
      !> Closure
      class(moist_model_moz_closure_hnc_type), intent(in) :: self
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Integrand partial derivatives
      real(wp), intent(out) :: dh, dc
      dh = closure_hnc_phi_dh(h, c)
      dc = closure_phi_dc(h)
   end subroutine hnc_local_derivatives

   !> KH local energy derivatives without checks
   pure subroutine kh_local_derivatives(self, h, c, dh, dc)
      !> Closure
      class(moist_model_moz_closure_kh_type), intent(in) :: self
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Integrand partial derivatives
      real(wp), intent(out) :: dh, dc
      dh = closure_kh_phi_dh(h, c)
      dc = closure_phi_dc(h)
   end subroutine kh_local_derivatives

   !> PSE local energy derivatives without checks; NaN for an invalid order
   pure subroutine pse_local_derivatives(self, h, c, dh, dc)
      !> Closure
      class(moist_model_moz_closure_pse_type), intent(in) :: self
      !> Physical total and direct correlations
      real(wp), intent(in) :: h, c
      !> Integrand partial derivatives
      real(wp), intent(out) :: dh, dc
      dc = closure_phi_dc(h)
      select case (self%order)
      case (2)
         dh = closure_pse2_phi_dh(h, c)
      case (3)
         dh = closure_pse3_phi_dh(h, c)
      case (4)
         dh = closure_pse4_phi_dh(h, c)
      case default
         dh = ieee_value(0.0_wp, ieee_quiet_nan)
      end select
   end subroutine pse_local_derivatives

   !> Check physical correlations before evaluating a local integrand
   !>
   !> @param[in] h Physical total correlation, dimensionless; 1 + h >= g_floor
   !> @param[in] c Physical direct correlation, dimensionless
   !> @param[out] error Invalid correlation input
   subroutine validate_integrand_inputs(h, c, error)
      !> Physical total correlation, dimensionless; 1 + h >= g_floor
      real(wp), intent(in) :: h
      !> Physical direct correlation, dimensionless
      real(wp), intent(in) :: c
      !> Invalid correlation input
      type(error_type), allocatable, intent(out) :: error

      if (.not. ieee_is_finite(h) .or. .not. ieee_is_finite(c)) then
         call fatal_error(error, "Closure energy correlations must be finite")
      else if (1.0_wp + h < g_floor) then
         call fatal_error(error, "Closure total correlation must be at least -1 within roundoff")
      end if
   end subroutine validate_integrand_inputs


   !> HNC second derivative with respect to t
   subroutine hnc_derivative2(self, u, t, beta, lambda, values, error)
      class(moist_model_moz_closure_hnc_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta, lambda
      real(wp), intent(out) :: values(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      call sweep_closure_dc(0, u, t, beta, lambda, values, self%team_size(size(t)), error, second=.true.)
   end subroutine hnc_derivative2

   !> KH second derivative with respect to t
   subroutine kh_derivative2(self, u, t, beta, lambda, values, error)
      class(moist_model_moz_closure_kh_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta, lambda
      real(wp), intent(out) :: values(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      call sweep_closure_dc(1, u, t, beta, lambda, values, self%team_size(size(t)), error, second=.true.)
   end subroutine kh_derivative2

   !> PSE second derivative with respect to t
   subroutine pse_derivative2(self, u, t, beta, lambda, values, error)
      class(moist_model_moz_closure_pse_type), intent(in) :: self
      real(wp), intent(in) :: u(:, :)
      real(wp), intent(in) :: t(:, :)
      real(wp), intent(in) :: beta, lambda
      real(wp), intent(out) :: values(:, :)
      type(error_type), allocatable, intent(out) :: error

      call validate_array_parameters(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      if (self%order < 2 .or. self%order > 4) then
         call fatal_error(error, "PSE order must be 2, 3 or 4")
         return
      end if
      call sweep_closure_dc(self%order, u, t, beta, lambda, values, self%team_size(size(t)), error, second=.true.)
   end subroutine pse_derivative2

   !> Report a missing second derivative in a custom closure
   subroutine unsupported_second_derivative(self, u, t, beta, lambda, values, error)
      !> Closure
      class(moist_model_moz_closure_type), intent(in) :: self
      !> Matching full or Ng potential/correlation, inverse temperature and scaling
      real(wp), intent(in) :: u(:, :), t(:, :), beta, lambda
      !> Second derivative; defined only on success
      real(wp), intent(out) :: values(:, :)
      !> Unsupported operation
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "This closure does not provide a second derivative")
   end subroutine unsupported_second_derivative

   !> Report a missing potential derivative in a custom closure
   subroutine unsupported_potential_derivative(self, u, t, beta, lambda, values, error)
      !> Closure
      class(moist_model_moz_closure_type), intent(in) :: self
      !> Matching full or Ng potential/correlation, inverse temperature and scaling
      real(wp), intent(in) :: u(:, :), t(:, :), beta, lambda
      !> Potential derivative; defined only on success
      real(wp), intent(out) :: values(:, :)
      !> Unsupported operation
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "This closure does not provide a potential derivative")
   end subroutine unsupported_potential_derivative

   !> Mark an absent pure hook without inventing response derivatives
   pure subroutine unsupported_local_energy_derivatives(self, h, c, dh, dc)
      class(moist_model_moz_closure_type), intent(in) :: self
      real(wp), intent(in) :: h, c
      real(wp), intent(out) :: dh, dc
      dh = ieee_value(0.0_wp, ieee_quiet_nan)
      dc = dh
   end subroutine unsupported_local_energy_derivatives

   !> Checked scalar derivatives through the optional pure hook
   subroutine default_energy_derivatives(self, h, c, dh, dc, error)
      !> Closure
      class(moist_model_moz_closure_type), intent(in) :: self
      !> Independent physical correlations
      real(wp), intent(in) :: h, c
      !> Local partial derivatives; defined only on success
      real(wp), intent(out) :: dh, dc
      !> Invalid input, missing derivatives or non-finite result
      type(error_type), allocatable, intent(out) :: error
      call validate_integrand_inputs(h, c, error)
      if (allocated(error)) return
      call self%local_energy_derivatives(h, c, dh, dc)
      if (.not. ieee_is_finite(dh) .or. .not. ieee_is_finite(dc)) then
         call fatal_error(error, "Closure energy integrand derivatives are unsupported or non-finite")
      end if
   end subroutine default_energy_derivatives

   !> Evaluate the HNC potential derivative using its closure kernel
   subroutine hnc_potential_derivative(self, u, t, beta, lambda, values, error)
      !> Closure
      class(moist_model_moz_closure_hnc_type), intent(in) :: self
      !> Matching full or Ng potential/correlation, inverse temperature and scaling
      real(wp), intent(in) :: u(:, :), t(:, :), beta, lambda
      !> Potential derivative at fixed t
      real(wp), intent(out) :: values(:, :)
      !> Invalid input or derivative
      type(error_type), allocatable, intent(out) :: error
      call self%derivative(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      values = closure_du(values, beta, lambda)
      call validate_result(values, error)
   end subroutine hnc_potential_derivative

   !> Evaluate HNC local energy derivatives without density or grid weights
   subroutine hnc_energy_derivatives(self, h, c, dh, dc, error)
      !> Closure
      class(moist_model_moz_closure_hnc_type), intent(in) :: self
      !> Independent physical correlations
      real(wp), intent(in) :: h, c
      !> Integrand partial derivatives with respect to h and c
      real(wp), intent(out) :: dh, dc
      !> Invalid input or derivative
      type(error_type), allocatable, intent(out) :: error
      call validate_integrand_inputs(h, c, error)
      if (allocated(error)) return
      dc = closure_phi_dc(h)
      dh = closure_hnc_phi_dh(h, c)
      if (.not. ieee_is_finite(dh) .or. .not. ieee_is_finite(dc)) then
         call fatal_error(error, "Closure energy derivative is non-finite")
      end if
   end subroutine hnc_energy_derivatives

   !> Evaluate the KH potential derivative using its closure kernel
   subroutine kh_potential_derivative(self, u, t, beta, lambda, values, error)
      !> Closure
      class(moist_model_moz_closure_kh_type), intent(in) :: self
      !> Matching full or Ng potential/correlation, inverse temperature and scaling
      real(wp), intent(in) :: u(:, :), t(:, :), beta, lambda
      !> Potential derivative at fixed t
      real(wp), intent(out) :: values(:, :)
      !> Invalid input or derivative
      type(error_type), allocatable, intent(out) :: error
      call self%derivative(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      values = closure_du(values, beta, lambda)
      call validate_result(values, error)
   end subroutine kh_potential_derivative

   !> Evaluate KH local energy derivatives without density or grid weights
   subroutine kh_energy_derivatives(self, h, c, dh, dc, error)
      !> Closure
      class(moist_model_moz_closure_kh_type), intent(in) :: self
      !> Independent physical correlations
      real(wp), intent(in) :: h, c
      !> Integrand partial derivatives with respect to h and c
      real(wp), intent(out) :: dh, dc
      !> Invalid input or derivative
      type(error_type), allocatable, intent(out) :: error
      call validate_integrand_inputs(h, c, error)
      if (allocated(error)) return
      dc = closure_phi_dc(h)
      dh = closure_kh_phi_dh(h, c)
      if (.not. ieee_is_finite(dh) .or. .not. ieee_is_finite(dc)) then
         call fatal_error(error, "Closure energy derivative is non-finite")
      end if
   end subroutine kh_energy_derivatives

   !> Evaluate the PSE potential derivative using its closure kernel
   subroutine pse_potential_derivative(self, u, t, beta, lambda, values, error)
      !> Closure
      class(moist_model_moz_closure_pse_type), intent(in) :: self
      !> Matching full or Ng potential/correlation, inverse temperature and scaling
      real(wp), intent(in) :: u(:, :), t(:, :), beta, lambda
      !> Potential derivative at fixed t
      real(wp), intent(out) :: values(:, :)
      !> Invalid input or derivative
      type(error_type), allocatable, intent(out) :: error
      call self%derivative(u, t, beta, lambda, values, error)
      if (allocated(error)) return
      values = closure_du(values, beta, lambda)
      call validate_result(values, error)
   end subroutine pse_potential_derivative

   !> Evaluate PSE local energy derivatives without density or grid weights
   subroutine pse_energy_derivatives(self, h, c, dh, dc, error)
      !> Closure
      class(moist_model_moz_closure_pse_type), intent(in) :: self
      !> Independent physical correlations
      real(wp), intent(in) :: h, c
      !> Integrand partial derivatives with respect to h and c
      real(wp), intent(out) :: dh, dc
      !> Invalid input or derivative
      type(error_type), allocatable, intent(out) :: error
      call validate_integrand_inputs(h, c, error)
      if (allocated(error)) return
      dc = closure_phi_dc(h)
      select case (self%order)
      case (2)
         dh = closure_pse2_phi_dh(h, c)
      case (3)
         dh = closure_pse3_phi_dh(h, c)
      case (4)
         dh = closure_pse4_phi_dh(h, c)
      case default
         call fatal_error(error, "PSE closure order must be 2, 3 or 4")
         return
      end select
      if (.not. ieee_is_finite(dh) .or. .not. ieee_is_finite(dc)) then
         call fatal_error(error, "Closure energy derivative is non-finite")
      end if
   end subroutine pse_energy_derivatives

end module moist_model_moz_closure
