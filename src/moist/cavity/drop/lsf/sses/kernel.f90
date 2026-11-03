!> Implicit Taylor coefficients for the optimized moment field
module moist_cavity_drop_lsf_sses_kernel
   use mctc_env_accuracy, only: wp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use moist_cavity_drop_lsf_sses_algebra, only: sses_jet_type, jet_constant, jet_seed, &
                                                 operator(+), operator(-), operator(*), &
                                                 operator(/), operator(**), exp, log
   use moist_cavity_drop_lsf_sses_native, only: moment_state_type, moment_init, moment_eval_native
   implicit none(type, external)
   private
   public :: sses_evaluate, sses_jets, sses_moment
contains

   !> Converge the moment variable from a cold start
   !>
   !> @param[in] point Evaluation point in Bohr
   !> @param[in] centers Atom centers in Bohr
   !> @param[in] radii Bare radii in Bohr
   !> @param[in] probe Probe radius in Bohr
   !> @param[in] epsilon Smoothing length in Bohr
   !> @param[in] tau Boltzmann temperature in Bohr^2
   !> @param[out] y0 Converged moment optimizer variable
   !> @param[out] valid Successful finite evaluation
   subroutine sses_moment(point, centers, radii, probe, epsilon, tau, y0, valid)
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Atom centers in Bohr
      real(wp), intent(in) :: centers(:, :)
      !> Bare radii in Bohr
      real(wp), intent(in) :: radii(:)
      !> Probe radius in Bohr
      real(wp), intent(in) :: probe
      !> Smoothing length in Bohr
      real(wp), intent(in) :: epsilon
      !> Boltzmann temperature in Bohr^2
      real(wp), intent(in) :: tau
      !> Converged moment optimizer variable
      real(wp), intent(out) :: y0(3)
      !> Successful finite evaluation
      logical, intent(out) :: valid
      type(moment_state_type) :: state
      real(wp) :: f, g(3), h(3, 3), stats(3)
      valid = .false.
      y0 = 0.0_wp
      if (size(radii) == 1) then
         valid = .true.
         return
      end if
      if (moment_init(state, size(radii), centers, radii, probe, epsilon, tau) /= 0) return
      if (moment_eval_native(state, point, 0, f, g, h, stats) /= 0) return
      y0 = state%y
      valid = .true.
   end subroutine sses_moment

   !> Solve the moment equation and differentiate it by Taylor degree
   !>
   !> @param[in] point Evaluation point in Bohr
   !> @param[in] centers Atom centers in Bohr
   !> @param[in] radii Bare radii in Bohr
   !> @param[in] probe Probe radius in Bohr
   !> @param[in] epsilon Smoothing length in Bohr
   !> @param[in] tau Boltzmann temperature in Bohr^2
   !> @param[in] ps Spatial seed directions
   !> @param[in] cs Atom-center seed directions
   !> @param[in] rs Radius seed directions
   !> @param[in] order Highest total derivative order
   !> @param[out] result Normalized Taylor coefficients
   !> @param[out] valid Successful finite evaluation
   subroutine sses_evaluate(point, centers, radii, probe, epsilon, tau, ps, cs, rs, order, result, valid)
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Atom centers in Bohr
      real(wp), intent(in) :: centers(:, :)
      !> Bare radii in Bohr
      real(wp), intent(in) :: radii(:)
      !> Probe radius in Bohr
      real(wp), intent(in) :: probe
      !> Smoothing length in Bohr
      real(wp), intent(in) :: epsilon
      !> Boltzmann temperature in Bohr^2
      real(wp), intent(in) :: tau
      !> Spatial seed directions
      real(wp), intent(in) :: ps(:, :)
      !> Atom-center seed directions
      real(wp), intent(in) :: cs(:, :, :)
      !> Radius seed directions
      real(wp), intent(in) :: rs(:, :)
      !> Highest total derivative order
      integer, intent(in) :: order
      !> Normalized Taylor coefficients
      type(sses_jet_type), intent(out) :: result
      !> Successful finite evaluation
      logical, intent(out) :: valid
      real(wp) :: y0(3)
      result = jet_constant(0.0_wp, max(0, min(order, 4)))
      call sses_moment(point, centers, radii, probe, epsilon, tau, y0, valid)
      if (.not. valid) return
      call sses_jets(point, centers, radii, probe, epsilon, tau, ps, cs, rs, order, y0, result, valid)
   end subroutine sses_evaluate

   !> Field jet and envelope atom-derivative jets from a converged moment
   !>
   !> gc(k, a) is the jet of dF/dc(k, a) and gr(a) the jet of dF/dr(a), both as
   !> functions of the five seeded variables
   !>
   !> @param[in] point Evaluation point in Bohr
   !> @param[in] centers Atom centers in Bohr
   !> @param[in] radii Bare radii in Bohr
   !> @param[in] probe Probe radius in Bohr
   !> @param[in] epsilon Smoothing length in Bohr
   !> @param[in] tau Boltzmann temperature in Bohr^2
   !> @param[in] ps Spatial seed directions
   !> @param[in] cs Atom-center seed directions
   !> @param[in] rs Radius seed directions
   !> @param[in] order Highest total derivative order
   !> @param[in] y0 Converged moment optimizer variable
   !> @param[out] result Normalized Taylor coefficients
   !> @param[out] valid Successful finite evaluation
   !> @param[out] gc Atom-center derivative jets dF/dc(k, a)
   !> @param[out] gr Radius derivative jets dF/dr(a)
   subroutine sses_jets(point, centers, radii, probe, epsilon, tau, ps, cs, rs, order, y0, result, valid, gc, gr)
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Atom centers in Bohr
      real(wp), intent(in) :: centers(:, :)
      !> Bare radii in Bohr
      real(wp), intent(in) :: radii(:)
      !> Probe radius in Bohr
      real(wp), intent(in) :: probe
      !> Smoothing length in Bohr
      real(wp), intent(in) :: epsilon
      !> Boltzmann temperature in Bohr^2
      real(wp), intent(in) :: tau
      !> Spatial seed directions
      real(wp), intent(in) :: ps(:, :)
      !> Atom-center seed directions
      real(wp), intent(in) :: cs(:, :, :)
      !> Radius seed directions
      real(wp), intent(in) :: rs(:, :)
      !> Highest total derivative order
      integer, intent(in) :: order
      !> Converged moment optimizer variable
      real(wp), intent(in) :: y0(3)
      !> Normalized Taylor coefficients
      type(sses_jet_type), intent(out) :: result
      !> Successful finite evaluation
      logical, intent(out) :: valid
      !> Atom-center derivative jets dF/dc(k, a)
      type(sses_jet_type), intent(out), optional :: gc(:, :)
      !> Radius derivative jets dF/dr(a)
      type(sses_jet_type), intent(out), optional :: gr(:)
      type(sses_jet_type), allocatable :: v(:, :), q(:), w(:), exponent(:), slope(:)
      type(sses_jet_type) :: rad, root, y(3), u(3), d, rd, z, m(3), residual(3), shift
      real(wp) :: mean(3), covariance(3, 3), dv(3), t(3, 3), jac(3, 3), inverse(3, 3), scalar_d, correction(3), rtau
      integer :: n, i, j, k, degree, coefficient
      integer, parameter :: counts(0:4) = [1, 6, 21, 56, 126]
      valid = .false.
      result = jet_constant(0.0_wp, max(0, min(order, 4)))
      n = size(radii)
      if (n < 1 .or. order < 0 .or. order > 4) return
      if (.not. ieee_is_finite(tau) .or. tau <= 0.0_wp) return
      allocate (v(3, n), q(n))
      if (present(gr)) allocate (slope(n))
      do i = 1, n
         rad = jet_seed(radii(i), rs(i, :), order)
         root = (rad*rad + epsilon**2)**0.5_wp
         q(i) = -rad*rad - 2.0_wp*probe*root
         if (present(gr)) slope(i) = -2.0_wp*rad - 2.0_wp*probe*rad/root
         do j = 1, 3
            v(j, i) = jet_seed(point(j) - centers(j, i), ps(j, :) - cs(j, i, :), order)
            q(i) = q(i) + v(j, i)*v(j, i)
         end do
      end do
      if (n == 1) then
         ! Exact squared-sphere override and its center and radius derivatives
         result = -rad*rad
         do j = 1, 3
            result = result + v(j, 1)*v(j, 1)
            if (present(gc)) gc(j, 1) = -2.0_wp*v(j, 1)
         end do
         if (present(gr)) gr(1) = -2.0_wp*rad
         valid = all(ieee_is_finite(result%c))
         return
      end if
      allocate (w(n), exponent(n))
      rtau = 1.0_wp/tau
      do j = 1, 3
         y(j) = jet_constant(y0(j), order)
         y(j)%constant = order == 0
      end do
      call weights()
      ! Scalar moment Jacobian J = I + 2/tau C T at the converged moment
      mean = 0.0_wp
      do i = 1, n
         mean = mean + w(i)%c(0)*v(:, i)%c(0)
      end do
      covariance = 0.0_wp
      do i = 1, n
         dv = v(:, i)%c(0) - mean
         do k = 1, 3
            covariance(:, k) = covariance(:, k) + w(i)%c(0)*dv*dv(k)
         end do
      end do
      scalar_d = sqrt(sum(y0**2) + epsilon**2)
      do k = 1, 3
         t(:, k) = -probe*y0*y0(k)/scalar_d**3
         t(k, k) = t(k, k) + probe/scalar_d
      end do
      jac = 2.0_wp*matmul(covariance, t)/tau
      do j = 1, 3
         jac(j, j) = jac(j, j) + 1.0_wp
      end do
      call inverse3(jac, inverse, valid)
      if (.not. valid) return
      ! Degree d of the residual is linear in degree d of y, with the scalar Jacobian
      do degree = 1, order
         do j = 1, 3
            m(j) = jet_constant(0.0_wp, order)
            do i = 1, n
               m(j) = m(j) + w(i)*v(j, i)
            end do
            residual(j) = y(j) - m(j)
         end do
         do coefficient = counts(degree - 1), counts(degree) - 1
            correction = matmul(inverse, residual(:)%c(coefficient))
            do j = 1, 3
               y(j)%c(coefficient) = y(j)%c(coefficient) - correction(j)
            end do
         end do
         call weights()
      end do
      result = shift - tau*log(z) + 2.0_wp*probe*epsilon**2*rd
      valid = all(ieee_is_finite(result%c))
      if (present(gc)) then
         do i = 1, n
            do j = 1, 3
               gc(j, i) = -2.0_wp*w(i)*(v(j, i) + u(j))
               valid = valid .and. all(ieee_is_finite(gc(j, i)%c))
            end do
         end do
      end if
      if (present(gr)) then
         ! Constant tau: the radius enters only through q_a
         do i = 1, n
            gr(i) = w(i)*slope(i)
            valid = valid .and. all(ieee_is_finite(gr(i)%c))
         end do
      end if
   contains
      !> Stable shifted Boltzmann weights and regularized probe displacement
      subroutine weights()
         type(sses_jet_type) :: rz
         d = jet_constant(epsilon**2, order)
         do j = 1, 3
            d = d + y(j)*y(j)
         end do
         rd = d**(-0.5_wp)
         do j = 1, 3
            u(j) = probe*y(j)*rd
         end do
         do i = 1, n
            exponent(i) = q(i)
            do j = 1, 3
               exponent(i) = exponent(i) + 2.0_wp*u(j)*v(j, i)
            end do
         end do
         ! Remove the full common Taylor jet, avoiding large far-field cancellations
         shift = exponent(minloc(exponent%c(0), dim=1))
         z = jet_constant(0.0_wp, order)
         do i = 1, n
            w(i) = exp(-((exponent(i) - shift)*rtau))
            z = z + w(i)
         end do
         rz = 1.0_wp/z
         do i = 1, n
            w(i) = w(i)*rz
         end do
      end subroutine weights
   end subroutine sses_jets

   !> Pivoted inverse of a scalar three-dimensional response matrix
   !>
   !> @param[in] a Scalar response matrix
   !> @param[out] inverse Inverse response matrix
   !> @param[out] valid Successful finite evaluation
   subroutine inverse3(a, inverse, valid)
      !> Scalar response matrix
      real(wp), intent(in) :: a(3, 3)
      !> Inverse response matrix
      real(wp), intent(out) :: inverse(3, 3)
      !> Successful finite evaluation
      logical, intent(out) :: valid
      real(wp) :: augmented(3, 6), row(6), factor
      integer :: i, k, pivot
      augmented(:, 1:3) = a
      augmented(:, 4:6) = 0.0_wp
      do i = 1, 3
         augmented(i, i + 3) = 1.0_wp
      end do
      valid = .false.
      do k = 1, 3
         pivot = maxloc(abs(augmented(k:3, k)), dim=1) + k - 1
         if (abs(augmented(pivot, k)) <= tiny(1.0_wp)) return
         row = augmented(k, :)
         augmented(k, :) = augmented(pivot, :)
         augmented(pivot, :) = row
         augmented(k, :) = augmented(k, :)/augmented(k, k)
         do i = 1, 3
            if (i == k) cycle
            factor = augmented(i, k)
            augmented(i, :) = augmented(i, :) - factor*augmented(k, :)
         end do
      end do
      inverse = augmented(:, 4:6)
      valid = .true.
   end subroutine inverse3

end module moist_cavity_drop_lsf_sses_kernel
