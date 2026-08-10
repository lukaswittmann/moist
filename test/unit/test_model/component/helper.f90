!> Reusable finite-difference harness for model component surface weights
!>
!> NOT a test; only provides helper funcitons
module test_model_component_helper
   use mctc_env, only: wp
   use test_helpers, only: fd4_scalar, fd4_offsets
   use testdrive, only: error_type, check, test_failed
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   implicit none
   private

   public :: surface_fixture
   public :: new_surface_fixture
   public :: check_surface_weights

   public :: fixture_ngrid_param
   public :: fixture_areas_param
   public :: fixture_xis_param
   public :: fixture_fs_param
   public :: fixture_xyz_param
   public :: fixture_normals_param
   public :: fixture_radial_normals

   !* --------------------- The shared 7-point synthetic surface --------------------- *!

   !> Number of points in the non-symmetric test surface
   integer, parameter :: fixture_ngrid_param = 7
   !> Synthetic surface areas
   real(wp), parameter :: fixture_areas_param(fixture_ngrid_param) = &
      [1.2_wp, 0.9_wp, 1.4_wp, 1.1_wp, 0.73_wp, 1.31_wp, 0.84_wp]
   !> Gaussian widths defining the DROP area channel
   real(wp), parameter :: fixture_xis_param(fixture_ngrid_param) = &
      [0.65_wp, 0.72_wp, 0.81_wp, 0.76_wp, 0.59_wp, 0.87_wp, 0.69_wp]
   !> Anchor switch factors defining the DROP area channel
   real(wp), parameter :: fixture_fs_param(fixture_ngrid_param) = &
      [0.91_wp, 0.83_wp, 0.88_wp, 0.79_wp, 0.94_wp, 0.74_wp, 0.86_wp]
   !> Surface positions
   real(wp), parameter :: fixture_xyz_param(3, fixture_ngrid_param) = reshape([ &
                                                  2.2_wp, 0.2_wp, 0.1_wp, &
                                                  -0.3_wp, 2.0_wp, 0.4_wp, &
                                                  0.1_wp, -0.2_wp, 2.4_wp, &
                                                  -1.5_wp, -1.2_wp, -1.0_wp, &
                                                  1.1_wp, -2.3_wp, 0.7_wp, &
                                                  -2.0_wp, 0.8_wp, 1.4_wp, &
                                                  0.6_wp, 1.3_wp, -2.1_wp], [3, fixture_ngrid_param])
   !> Tilted normal field: deliberately NOT the normalised positions, so a
   !> kernel that silently assumes a radial surface fails here.
   real(wp), parameter :: fixture_normals_param(3, fixture_ngrid_param) = reshape([ &
                                          0.95_wp, 0.28_wp, 0.14_wp, &
                                          -0.11_wp, 0.97_wp, 0.21_wp, &
                                          0.08_wp, -0.16_wp, 0.98_wp, &
                                          -0.69_wp, -0.55_wp, -0.46_wp, &
                                          0.42_wp, -0.88_wp, 0.22_wp, &
                                          -0.78_wp, 0.31_wp, 0.54_wp, &
                                          0.23_wp, 0.49_wp, -0.84_wp], [3, fixture_ngrid_param])

   !> Channel selectors for `sample_channel`.
   integer, parameter :: channel_xi = 1
   integer, parameter :: channel_f = 2
   integer, parameter :: channel_xyz = 3
   integer, parameter :: channel_normal = 4

   !> Independent DROP surface variables used by component energy callbacks.
   type :: surface_fixture
      !> Area prefactor held fixed when xi or f is perturbed
      real(wp), allocatable :: area_base(:)
      !> Gaussian widths
      real(wp), allocatable :: xi(:)
      !> Anchor switch factors
      real(wp), allocatable :: f(:)
      !> Projected surface positions
      real(wp), allocatable :: xyz(:, :)
      !> Outward surface normals
      real(wp), allocatable :: normal(:, :)
   contains
      !> Reconstruct quadrature areas from the independent DROP variables.
      procedure :: areas => fixture_areas
      !> Return the number of surface points.
      procedure :: ngrid => fixture_ngrid
   end type surface_fixture

   abstract interface
      !> Evaluate the cavity-dependent component energy on a surface fixture.
      function surface_energy_callback(surface) result(energy)
         import :: surface_fixture, wp
         !> Surface variables at which to evaluate the energy
         type(surface_fixture), intent(in) :: surface
         !> Component energy
         real(wp) :: energy
      end function surface_energy_callback
   end interface

contains

   !> Construct a surface fixture while preserving the DROP area invariant.
   !> @param[out] surface Fixture to initialize
   !> @param[in]  areas   Surface quadrature areas
   !> @param[in]  xi      Gaussian widths
   !> @param[in]  f       Anchor switch factors
   !> @param[in]  xyz     Projected surface positions
   !> @param[in]  normal  Outward surface normals
   subroutine new_surface_fixture(surface, areas, xi, f, xyz, normal)
      !> Fixture to initialize
      type(surface_fixture), intent(out) :: surface
      !> Surface quadrature areas
      real(wp), intent(in) :: areas(:)
      !> Gaussian widths
      real(wp), intent(in) :: xi(:)
      !> Anchor switch factors
      real(wp), intent(in) :: f(:)
      !> Projected surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Outward surface normals
      real(wp), intent(in) :: normal(:, :)

      integer :: ngrid

      ngrid = size(areas)
      if (size(xi) /= ngrid .or. size(f) /= ngrid) &
         error stop "new_surface_fixture: scalar channel size mismatch"
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= ngrid) &
         error stop "new_surface_fixture: xyz shape mismatch"
      if (size(normal, 1) /= 3 .or. size(normal, 2) /= ngrid) &
         error stop "new_surface_fixture: normal shape mismatch"
      if (any(abs(f) <= tiny(1.0_wp))) &
         error stop "new_surface_fixture: f must be nonzero"

      allocate (surface%area_base, source=areas*xi**2/f)
      allocate (surface%xi, source=xi)
      allocate (surface%f, source=f)
      allocate (surface%xyz, source=xyz)
      allocate (surface%normal, source=normal)
   end subroutine new_surface_fixture

   !> Reconstruct quadrature areas for the current fixture variables.
   !> @param[in] self Surface fixture
   !> @return Surface quadrature areas
   pure function fixture_areas(self) result(areas)
      !> Surface fixture
      class(surface_fixture), intent(in) :: self
      !> Reconstructed surface quadrature areas
      real(wp), allocatable :: areas(:)

      allocate (areas, source=self%area_base*self%f/self%xi**2)
   end function fixture_areas

   !> Return the number of points in a surface fixture.
   !> @param[in] self Surface fixture
   !> @return Number of surface points
   pure integer function fixture_ngrid(self) result(ngrid)
      !> Surface fixture
      class(surface_fixture), intent(in) :: self

      ngrid = size(self%xi)
   end function fixture_ngrid

   !> Check DROP surface weights against fourth-order central differences.
   !>
   !> The energy callback receives the complete perturbed surface, so coupled
   !> point terms are tested without special treatment. Individual channels can
   !> be disabled only when a component does not depend on that surface variable.
   !>
   !> @param[out] error       Test failure information
   !> @param[in]  surface     Reference surface fixture
   !> @param[in]  weights     Analytic component surface weights
   !> @param[in]  evaluate    Independent component energy callback
   !> @param[in]  label       Component label used in failure messages
   !> @param[in]  step        Finite-difference displacement
   !> @param[in]  thr_abs     Absolute comparison tolerance
   !> @param[in]  thr_rel     Relative comparison tolerance
   !> @param[in]  check_xi    Check the Gaussian-width channel
   !> @param[in]  check_f     Check the anchor-switch channel
   !> @param[in]  check_xyz   Check the position channel
   !> @param[in]  check_normal Check the normal channel
   subroutine check_surface_weights(error, surface, weights, evaluate, label, &
                                         step, thr_abs, thr_rel, check_xi, check_f, check_xyz, check_normal)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Reference surface fixture
      type(surface_fixture), intent(in) :: surface
      !> Analytic component surface weights
      type(cavity_surface_adjoint_type), intent(in) :: weights
      !> Independent component energy callback
      procedure(surface_energy_callback) :: evaluate
      !> Component label used in failure messages
      character(len=*), intent(in) :: label
      !> Finite-difference displacement
      real(wp), intent(in), optional :: step
      !> Absolute comparison tolerance
      real(wp), intent(in), optional :: thr_abs
      !> Relative comparison tolerance
      real(wp), intent(in), optional :: thr_rel
      !> Check the Gaussian-width channel
      logical, intent(in), optional :: check_xi
      !> Check the anchor-switch channel
      logical, intent(in), optional :: check_f
      !> Check the position channel
      logical, intent(in), optional :: check_xyz
      !> Check the normal channel
      logical, intent(in), optional :: check_normal

      real(wp) :: h, atol, rtol, fd, analytic
      !> Stencil samples, in `fd4_offsets` order
      real(wp) :: vals(4)
      real(wp), allocatable :: areas(:)
      logical :: do_xi, do_f, do_xyz, do_normal
      integer :: axis, igrid, ngrid
      character(len=160) :: context

      ngrid = surface%ngrid()
      areas = surface%areas()
      if (.not. weights_are_valid(weights, ngrid)) then
         call test_failed(error, trim(label)//": invalid surface-weight shapes")
         return
      end if

      h = 2.0e-5_wp
      if (present(step)) h = step
      atol = 5.0e-11_wp
      if (present(thr_abs)) atol = thr_abs
      rtol = 5.0e-10_wp
      if (present(thr_rel)) rtol = thr_rel
      do_xi = .true.
      if (present(check_xi)) do_xi = check_xi
      do_f = .true.
      if (present(check_f)) do_f = check_f
      do_xyz = .true.
      if (present(check_xyz)) do_xyz = check_xyz
      do_normal = .true.
      if (present(check_normal)) do_normal = check_normal

      do igrid = 1, ngrid
         if (do_xi) then
            call sample_channel(surface, evaluate, igrid, channel_xi, 1, h, vals)
            fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), h)
            write (context, '(A," xi point ",I0)') trim(label), igrid
            analytic = weights%w_xi(igrid) &
                       - 2.0_wp*areas(igrid)*weights%w_a(igrid)/surface%xi(igrid)
            call check(error, analytic, fd, thr_abs=atol, &
                       thr_rel=rtol, more=trim(context))
            if (allocated(error)) return
         end if

         if (do_f) then
            call sample_channel(surface, evaluate, igrid, channel_f, 1, h, vals)
            fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), h)
            write (context, '(A," f point ",I0)') trim(label), igrid
            analytic = weights%w_f(igrid) &
                       + areas(igrid)*weights%w_a(igrid)/surface%f(igrid)
            call check(error, analytic, fd, thr_abs=atol, &
                       thr_rel=rtol, more=trim(context))
            if (allocated(error)) return
         end if

         do axis = 1, 3
            if (do_xyz) then
               call sample_channel(surface, evaluate, igrid, channel_xyz, axis, h, vals)
               fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), h)
               write (context, '(A," xyz(",I0,") point ",I0)') &
                  trim(label), axis, igrid
               call check(error, weights%w_xyz(axis, igrid), fd, &
                          thr_abs=atol, thr_rel=rtol, more=trim(context))
               if (allocated(error)) return
            end if

            if (do_normal) then
               call sample_channel(surface, evaluate, igrid, channel_normal, axis, h, vals)
               fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), h)
               write (context, '(A," normal(",I0,") point ",I0)') &
                  trim(label), axis, igrid
               call check(error, weights%w_n(axis, igrid), fd, &
                          thr_abs=atol, thr_rel=rtol, more=trim(context))
               if (allocated(error)) return
            end if
         end do
      end do
   end subroutine check_surface_weights

   !> Radial unit normals of the shared fixture: `xyz` normalised point by point.
   !>
   !> A genuinely different surface from `fixture_normals_param` on the same
   !> point set, for components whose energy is only defined on a closed radial
   !> surface; swapping in the tilted field would silently change what such a
   !> suite tests. Kept as a function because `norm2` is not permitted in a
   !> constant expression.
   !>
   !> @return  Unit normals, (3, fixture_ngrid_param)
   pure function fixture_radial_normals() result(normals)
      !> Radially outward unit normals
      real(wp) :: normals(3, fixture_ngrid_param)
      !> Point index
      integer :: i

      do i = 1, fixture_ngrid_param
         normals(:, i) = fixture_xyz_param(:, i)/norm2(fixture_xyz_param(:, i))
      end do
   end function fixture_radial_normals

   !> Sample one surface channel at the four `fd4_offsets` displacements.
   !>
   !> Every channel perturbs a single scalar slot of the fixture, so the stencil
   !> order lives in exactly one place here rather than being re-derived per
   !> channel -- a reversed stencil would otherwise silently flip the sign of
   !> whichever channel got it wrong.
   !>
   !> @param[in]  surface   Reference surface fixture
   !> @param[in]  evaluate  Component energy callback
   !> @param[in]  igrid     Surface point index
   !> @param[in]  channel   One of `channel_xi`, `channel_f`, `channel_xyz`, `channel_normal`
   !> @param[in]  axis      Cartesian axis; ignored by the scalar channels
   !> @param[in]  h         Finite-difference displacement
   !> @param[out] vals      Energies, in `fd4_offsets` order
   subroutine sample_channel(surface, evaluate, igrid, channel, axis, h, vals)
      !> Reference surface fixture
      type(surface_fixture), intent(in) :: surface
      !> Component energy callback
      procedure(surface_energy_callback) :: evaluate
      !> Surface point index
      integer, intent(in) :: igrid
      !> Channel selector
      integer, intent(in) :: channel
      !> Cartesian axis; ignored by the scalar channels
      integer, intent(in) :: axis
      !> Finite-difference displacement
      real(wp), intent(in) :: h
      !> Energies, in `fd4_offsets` order
      real(wp), intent(out) :: vals(4)

      !> Perturbed copy of the fixture
      type(surface_fixture) :: trial
      !> Stencil index
      integer :: k

      trial = surface
      do k = 1, size(fd4_offsets)
         select case (channel)
         case (channel_xi)
            trial%xi(igrid) = surface%xi(igrid) + fd4_offsets(k)*h
         case (channel_f)
            trial%f(igrid) = surface%f(igrid) + fd4_offsets(k)*h
         case (channel_xyz)
            trial%xyz(axis, igrid) = surface%xyz(axis, igrid) + fd4_offsets(k)*h
         case (channel_normal)
            trial%normal(axis, igrid) = surface%normal(axis, igrid) + fd4_offsets(k)*h
         case default
            error stop "sample_channel: invalid channel"
         end select
         vals(k) = evaluate(trial)
      end do
   end subroutine sample_channel

   !> Check that every DROP weight channel has the fixture shape.
   pure logical function weights_are_valid(weights, ngrid) result(valid)
      type(cavity_surface_adjoint_type), intent(in) :: weights
      integer, intent(in) :: ngrid

      valid = allocated(weights%w_xi) .and. allocated(weights%w_f) &
              .and. allocated(weights%w_a) .and. allocated(weights%w_w) &
              .and. allocated(weights%w_xyz) .and. allocated(weights%w_n)
      if (.not. valid) return
      valid = size(weights%w_xi) == ngrid .and. size(weights%w_f) == ngrid &
              .and. size(weights%w_a) == ngrid .and. size(weights%w_w) == ngrid &
              .and. size(weights%w_xyz, 1) == 3 .and. size(weights%w_xyz, 2) == ngrid &
              .and. size(weights%w_n, 1) == 3 .and. size(weights%w_n, 2) == ngrid
   end function weights_are_valid

end module test_model_component_helper
