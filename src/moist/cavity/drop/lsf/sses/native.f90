!> Converged moment optimizer and analytic spatial response
module moist_cavity_drop_lsf_sses_native
   use mctc_env_accuracy, only: wp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use moist_cavity_drop_lsf_sses_spatial, only: moment_response_matrices, moment_hessian_product, &
                                                 moment_third_rhs, moment_third_product
   implicit none(type, external)
   private
   public :: moment_state_type, moment_init, moment_eval_native, moment_third_native, moment_atom_native

   !> Value-owned geometry and point cache, copied independently for each DROP thread
   type :: moment_state_type
      !> Atom count
      integer :: n = 0
      !> Highest total derivative order
      integer :: order = -1
      !> Successful finite evaluation
      logical :: valid = .false.
      !> Available moment sensitivity
      logical :: sensitivity_valid = .false.
      !> Available third derivative
      logical :: third_valid = .false.
      !> Probe radius
      real(wp) :: p = 0.0_wp
      !> Smoothing length
      real(wp) :: eps = 0.0_wp
      !> Boltzmann temperature in Bohr^2
      real(wp) :: tau = 0.0_wp
      !> Atom centers in Bohr
      real(wp), allocatable :: xyz(:, :)
      !> Geometry coefficients
      real(wp), allocatable :: offset(:)
      !> Slope
      real(wp), allocatable :: slope(:)
      !> Relative coordinates
      real(wp), allocatable :: v(:, :)
      !> Power distances
      real(wp), allocatable :: q(:)
      !> Shifted Boltzmann weights
      real(wp), allocatable :: w(:)
      !> Cached evaluation point
      real(wp) :: x(3) = 0.0_wp
      !> Moment optimizer variable
      real(wp) :: y(3) = 0.0_wp
      !> Optimized first moment
      real(wp) :: m(3) = 0.0_wp
      !> Weighted covariance
      real(wp) :: covariance(3, 3) = 0.0_wp
      !> Cached field value
      real(wp) :: value = 0.0_wp
      !> Cached spatial gradient
      real(wp) :: gradient(3) = 0.0_wp
      !> Cached spatial Hessian
      real(wp) :: hessian(3, 3) = 0.0_wp
      !> Moment sensitivity
      real(wp) :: sensitivity(3, 3) = 0.0_wp
      !> Softmax shift
      real(wp) :: lo = 0.0_wp
      !> Shifted Boltzmann sum
      real(wp) :: total = 0.0_wp
      !> Stationarity residual
      real(wp) :: residual = 0.0_wp
      !> Third spatial derivative
      real(wp) :: third(3, 3, 3) = 0.0_wp
   end type moment_state_type
   
contains

   !> Initialize immutable geometry and per-thread workspace
   !>
   !> @param[out] state Native moment state
   !> @param[in] n Atom count
   !> @param[in] xyz Atom centers in Bohr
   !> @param[in] radii Bare radii in Bohr
   !> @param[in] probe Probe radius in Bohr
   !> @param[in] epsilon Smoothing length in Bohr
   !> @param[in] tau Boltzmann temperature in Bohr^2
   function moment_init(state, n, xyz, radii, probe, epsilon, tau) result(status)
      !> Native moment state
      type(moment_state_type), intent(out) :: state
      !> Atom count
      integer, intent(in) :: n
      !> Atom centers in Bohr
      real(wp), intent(in) :: xyz(3, n)
      !> Bare radii in Bohr
      real(wp), intent(in) :: radii(n)
      !> Probe radius in Bohr
      real(wp), intent(in) :: probe
      !> Smoothing length in Bohr
      real(wp), intent(in) :: epsilon
      !> Boltzmann temperature in Bohr^2
      real(wp), intent(in) :: tau
      !> Native status code
      integer :: status
      integer :: atom, allocation_status
      status = 1
      if (n < 1) return
      if (.not. ieee_is_finite(probe) .or. .not. ieee_is_finite(epsilon)) return
      if (.not. ieee_is_finite(tau)) return
      if (probe < 0.0_wp .or. epsilon <= 0.0_wp .or. tau <= 0.0_wp) return
      if (.not. all(ieee_is_finite(radii)) .or. .not. all(ieee_is_finite(xyz))) return
      if (any(radii <= 0.0_wp)) return
      status = 2
      allocate (state%xyz(3, n), stat=allocation_status)
      if (allocation_status /= 0) return
      allocate (state%offset(n), stat=allocation_status)
      if (allocation_status /= 0) return
      allocate (state%v(3, n), stat=allocation_status)
      if (allocation_status /= 0) return
      allocate (state%q(n), stat=allocation_status)
      if (allocation_status /= 0) return
      allocate (state%w(n), stat=allocation_status)
      if (allocation_status /= 0) return
      allocate (state%slope(n), stat=allocation_status)
      if (allocation_status /= 0) return
      state%n = n
      state%p = probe
      state%eps = epsilon
      state%xyz = xyz
      state%tau = tau
      do atom = 1, n
         state%offset(atom) = -radii(atom)**2 - 2.0_wp*probe*norm2([radii(atom), epsilon])
         state%slope(atom) = -2.0_wp*radii(atom) - 2.0_wp*probe*radii(atom)/norm2([radii(atom), epsilon])
      end do
      status = 0
   end function moment_init

   !> Factor a three-dimensional response matrix with partial pivoting
   !>
   !> @param[in,out] a Scalar response matrix
   !> @param[out] pivots Row pivot indices
   !> @param[out] valid Successful finite evaluation
   pure subroutine factor3(a, pivots, valid)
      !> Scalar response matrix
      real(wp), intent(inout) :: a(3, 3)
      !> Row pivot indices
      integer, intent(out) :: pivots(3)
      !> Successful finite evaluation
      logical, intent(out) :: valid
      real(wp) :: row(3)
      integer :: k, pivot, i, j
      valid = .false.
      do k = 1, 3
         pivot = maxloc(abs(a(k:3, k)), dim=1) + k - 1
         if (abs(a(pivot, k)) <= tiny(1.0_wp)) return
         pivots(k) = pivot
         row = a(k, :)
         a(k, :) = a(pivot, :)
         a(pivot, :) = row
         do i = k + 1, 3
            a(i, k) = a(i, k)/a(k, k)
            do j = k + 1, 3
               a(i, j) = a(i, j) - a(i, k)*a(k, j)
            end do
         end do
      end do
      valid = .true.
   end subroutine factor3

   !> Solve one right-hand side with an existing pivoted factorization
   !>
   !> @param[in] a Scalar response matrix
   !> @param[in] pivots Row pivot indices
   !> @param[in,out] x Cached evaluation point
   pure subroutine solve3(a, pivots, x)
      !> Scalar response matrix
      real(wp), intent(in) :: a(3, 3)
      !> Row pivot indices
      integer, intent(in) :: pivots(3)
      !> Cached evaluation point
      real(wp), intent(inout) :: x(3)
      real(wp) :: temporary
      integer :: k, i, j
      do k = 1, 3
         temporary = x(k)
         x(k) = x(pivots(k))
         x(pivots(k)) = temporary
      end do
      do i = 1, 3
         do j = 1, i - 1
            x(i) = x(i) - a(i, j)*x(j)
         end do
      end do
      do i = 3, 1, -1
         do j = i + 1, 3
            x(i) = x(i) - a(i, j)*x(j)
         end do
         x(i) = x(i)/a(i, i)
      end do
   end subroutine solve3

   !> Accumulate shifted Boltzmann weights and centered moments without logarithms
   !>
   !> @param[in,out] state Native moment state
   !> @param[in] u Regularized probe displacement
   !> @param[out] m Optimized first moment
   !> @param[out] c Weighted covariance
   !> @param[out] lo Softmax shift
   !> @param[out] total Shifted Boltzmann sum
   subroutine moments(state, u, m, c, lo, total)
      !> Native moment state
      type(moment_state_type), intent(inout) :: state
      !> Regularized probe displacement
      real(wp), intent(in) :: u(3)
      !> Optimized first moment
      real(wp), intent(out) :: m(3)
      !> Weighted covariance
      real(wp), intent(out) :: c(3, 3)
      !> Softmax shift
      real(wp), intent(out) :: lo
      !> Shifted Boltzmann sum
      real(wp), intent(out) :: total
      real(wp) :: dv(3), weight, inverse
      integer :: atom, i, j
      lo = huge(1.0_wp)
      do atom = 1, state%n
         state%w(atom) = state%q(atom) + 2.0_wp*dot_product(u, state%v(:, atom))
         lo = min(lo, state%w(atom))
      end do
      total = 0.0_wp
      do atom = 1, state%n
         state%w(atom) = exp(-(state%w(atom) - lo)/state%tau)
         total = total + state%w(atom)
      end do
      m = 0.0_wp
      do atom = 1, state%n
         m = m + state%w(atom)*state%v(:, atom)
      end do
      inverse = 1.0_wp/total
      m = m*inverse
      c = 0.0_wp
      do atom = 1, state%n
         dv = state%v(:, atom) - m
         weight = state%w(atom)*inverse
         do j = 1, 3
            do i = 1, j
               c(i, j) = c(i, j) + weight*dv(i)*dv(j)
            end do
         end do
      end do
      do j = 1, 3
         do i = j + 1, 3
            c(i, j) = c(j, i)
         end do
      end do
   end subroutine moments

   !> Complete the Hessian and moment sensitivity from the converged first moment
   !>
   !> @param[in,out] state Native moment state
   !> @param[in] order Highest total derivative order
   function finish(state, order) result(status)
      !> Native moment state
      type(moment_state_type), intent(inout) :: state
      !> Highest total derivative order
      integer, intent(in) :: order
      !> Native status code
      integer :: status
      real(wp) :: d, t(3, 3), jac(3, 3), z(3, 3), rhs(3)
      integer :: pivots(3), k
      logical :: valid
      status = 0
      if (order < 2) return
      d = sqrt(dot_product(state%y, state%y) + state%eps**2)
      call moment_response_matrices(state%y, state%covariance, state%p, state%tau, 1.0_wp/d, t, jac)
      call factor3(jac, pivots, valid)
      if (.not. valid) then
         status = 2
         return
      end if
      do k = 1, 3
         rhs = -4.0_wp*state%covariance(:, k)/state%tau
         rhs(k) = rhs(k) + 2.0_wp
         call solve3(jac, pivots, rhs)
         z(:, k) = rhs
      end do
      state%sensitivity = 0.5_wp*z
      call moment_hessian_product(t, z, state%hessian)
      state%hessian = 0.5_wp*(state%hessian + transpose(state%hessian))
      state%sensitivity_valid = .true.
   end function finish

   !> Evaluate the field through order two with safeguarded warm-started Newton steps
   !>
   !> @param[in,out] state Native moment state
   !> @param[in] point Evaluation point in Bohr
   !> @param[in] order Highest total derivative order
   !> @param[out] f Field value
   !> @param[out] g Spatial gradient
   !> @param[out] h Spatial Hessian
   !> @param[out] stats Iterations, atom passes and residual
   function moment_eval_native(state, point, order, f, g, h, stats) result(status)
      !> Native moment state
      type(moment_state_type), intent(inout) :: state
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Highest total derivative order
      integer, intent(in) :: order
      !> Field value
      real(wp), intent(out) :: f
      !> Spatial gradient
      real(wp), intent(out) :: g(3)
      !> Spatial Hessian
      real(wp), intent(out) :: h(3, 3)
      !> Iterations, atom passes and residual
      real(wp), intent(out) :: stats(3)
      !> Native status code
      integer :: status
      real(wp) :: initial(3), y(3), m(3), c(3, 3), lo, total, d, u(3), residual(3), res
      real(wp) :: t(3, 3), jac(3, 3), step(3), alpha, trial(3), tm(3), tc(3, 3), td, tl, ts
      integer :: atom, attempt, iteration, passes, ls, pivots(3)
      logical :: same, warm, converged, accepted, valid
      status = 1
      f = 0.0_wp
      g = 0.0_wp
      h = 0.0_wp
      stats = 0.0_wp
      if (state%n < 1 .or. order < 0 .or. order > 2) return
      if (.not. all(ieee_is_finite(point))) return
      stats(3) = state%residual
      same = state%valid .and. all(point == state%x)
      if (same) then
         if (order > state%order) then
            status = finish(state, order)
            if (status /= 0) return
            state%order = order
         end if
      else
         do atom = 1, state%n
            state%v(:, atom) = point - state%xyz(:, atom)
            state%q(atom) = dot_product(state%v(:, atom), state%v(:, atom)) + state%offset(atom)
         end do
         warm = state%valid
         initial = state%y
         if (warm) then
            if (state%sensitivity_valid) then
               initial = initial + matmul(state%sensitivity, point - state%x)
            else
               initial = initial + point - state%x
            end if
         end if
         state%valid = .false.
         state%sensitivity_valid = .false.
         state%third_valid = .false.
         converged = .false.
         passes = 0
         iteration = 0
         res = huge(1.0_wp)
         do attempt = 1, merge(2, 1, warm)
            if (warm .and. attempt == 1) then
               y = initial
            else
               call moments(state, [0.0_wp, 0.0_wp, 0.0_wp], y, c, lo, total)
               passes = passes + 1
            end if
            d = sqrt(dot_product(y, y) + state%eps**2)
            u = state%p*y/d
            call moments(state, u, m, c, lo, total)
            passes = passes + 1
            do iteration = 0, 99
               residual = y - m
               res = sqrt(dot_product(residual, residual))
               if (res <= 1.0e-11_wp*(1.0_wp + sqrt(dot_product(m, m)))) then
                  converged = .true.
                  exit
               end if
               call moment_response_matrices(y, c, state%p, state%tau, 1.0_wp/d, t, jac)
               call factor3(jac, pivots, valid)
               if (.not. valid) exit
               step = -residual
               call solve3(jac, pivots, step)
               alpha = 1.0_wp
               accepted = .false.
               do ls = 1, 60
                  trial = y + alpha*step
                  td = sqrt(dot_product(trial, trial) + state%eps**2)
                  call moments(state, state%p*trial/td, tm, tc, tl, ts)
                  passes = passes + 1
                  residual = trial - tm
                  if (sqrt(dot_product(residual, residual)) <= (1.0_wp - 1.0e-4_wp*alpha)*res) then
                     y = trial
                     m = tm
                     c = tc
                     lo = tl
                     total = ts
                     d = td
                     accepted = .true.
                     exit
                  end if
                  alpha = 0.5_wp*alpha
               end do
               if (.not. accepted) exit
            end do
            if (converged) exit
         end do
         stats = [real(iteration, wp), real(passes, wp), res]
         status = 2
         if (.not. converged) return
         state%x = point
         state%y = y
         state%m = m
         state%covariance = c
         d = sqrt(dot_product(y, y) + state%eps**2)
         state%value = lo - state%tau*log(total) + 2.0_wp*state%p*state%eps**2/d
         state%gradient = 2.0_wp*(m + state%p*y/d)
         state%residual = res
         state%lo = lo
         state%total = total
         status = finish(state, order)
         if (status /= 0) return
         state%order = order
         state%valid = .true.
      end if
      f = state%value
      if (order >= 1) g = state%gradient
      if (order >= 2) h = state%hessian
      status = 0
   end function moment_eval_native

   !> Spatial derivatives of the probe Jacobian and of twice the moment sensitivity
   !>
   !> @param[in] state Native moment state
   !> @param[out] t Probe displacement Jacobian
   !> @param[out] z Twice the moment sensitivity
   !> @param[out] dt Spatial derivative of the probe Jacobian
   !> @param[out] dz Spatial derivative of twice the moment sensitivity
   !> @param[out] valid Successful finite evaluation
   pure subroutine third_response(state, t, z, dt, dz, valid)
      !> Native moment state
      type(moment_state_type), intent(in) :: state
      !> Probe displacement Jacobian
      real(wp), intent(out) :: t(3, 3)
      !> Twice the moment sensitivity
      real(wp), intent(out) :: z(3, 3)
      !> Spatial derivative of the probe Jacobian
      real(wp), intent(out) :: dt(3, 3, 3)
      !> Spatial derivative of twice the moment sensitivity
      real(wp), intent(out) :: dz(3, 3, 3)
      !> Successful finite evaluation
      logical, intent(out) :: valid
      real(wp) :: cm3(3, 3, 3), dv(3), weight, d, jac(3, 3), rhs(3, 3, 3)
      integer :: atom, i, j, k, column, pivots(3)
      cm3 = 0.0_wp
      do atom = 1, state%n
         dv = state%v(:, atom) - state%m
         weight = state%w(atom)/state%total
         do k = 1, 3
            do j = 1, k
               do i = 1, j
                  cm3(i, j, k) = cm3(i, j, k) + weight*dv(i)*dv(j)*dv(k)
               end do
            end do
         end do
      end do
      do k = 1, 3
         do j = 1, k
            do i = 1, j
               weight = cm3(i, j, k)
               cm3(i, k, j) = weight
               cm3(j, i, k) = weight
               cm3(j, k, i) = weight
               cm3(k, i, j) = weight
               cm3(k, j, i) = weight
            end do
         end do
      end do
      d = sqrt(dot_product(state%y, state%y) + state%eps**2)
      call moment_response_matrices(state%y, state%covariance, state%p, state%tau, 1.0_wp/d, t, jac)
      call factor3(jac, pivots, valid)
      if (.not. valid) return
      z = 2.0_wp*state%sensitivity
      call moment_third_rhs(state%y, state%covariance, cm3, state%sensitivity, t, z, &
                            state%p, state%tau, 1.0_wp/d, dt, rhs)
      dz = rhs
      do k = 1, 3
         do column = 1, 3
            call solve3(jac, pivots, dz(:, column, k))
         end do
      end do
   end subroutine third_response

   !> Evaluate the third spatial derivative from centered third moments
   !>
   !> @param[in,out] state Native moment state
   !> @param[in] point Evaluation point in Bohr
   !> @param[out] third Third spatial derivative
   !> @param[out] stats Iterations, atom passes and residual
   function moment_third_native(state, point, third, stats) result(status)
      !> Native moment state
      type(moment_state_type), intent(inout) :: state
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Third spatial derivative
      real(wp), intent(out) :: third(3, 3, 3)
      !> Iterations, atom passes and residual
      real(wp), intent(out) :: stats(3)
      !> Native status code
      integer :: status
      real(wp) :: f, g(3), h(3, 3), t(3, 3), z(3, 3), dt(3, 3, 3), dz(3, 3, 3), raw(3, 3, 3)
      integer :: i, j, k
      logical :: valid
      third = 0.0_wp
      status = moment_eval_native(state, point, 2, f, g, h, stats)
      if (status /= 0) return
      if (.not. state%third_valid) then
         call third_response(state, t, z, dt, dz, valid)
         if (.not. valid) then
            status = 2
            return
         end if
         call moment_third_product(t, z, dt, dz, raw)
         do k = 1, 3
            do j = 1, 3
               do i = 1, 3
                  state%third(i, j, k) = (raw(i, j, k) + raw(i, k, j) + raw(j, i, k) &
                                          + raw(j, k, i) + raw(k, i, j) + raw(k, j, i))/6.0_wp
               end do
            end do
         end do
         state%third_valid = .true.
         stats(2) = stats(2) + 1.0_wp
      end if
      third = state%third
   end function moment_third_native

   !> Envelope atom derivatives and their spatial derivatives through order two
   !>
   !> With weights w, exponents e = q + 2 u.v and u(x) the converged probe
   !> displacement, dF/dc(s, a) = -2 w_a (v_a + u)_s and dF/dr(a) = w_a q'_a
   !> with tau constant. Spatial derivatives use dw_a = -2 w_a/tau B d_a
   !> with B = I + du^T/dx and d_a = v_a - m, so every atom costs O(1).
   !> Requires the converged state of the point with its moment sensitivity.
   !>
   !> @param[in] state Native moment state
   !> @param[in] order Highest total derivative order
   !> @param[out] valid Successful finite evaluation
   !> @param[out] gc0 gc0
   !> @param[out] gc1 gc1
   !> @param[out] gc2 gc2
   !> @param[out] gr0 gr0
   !> @param[out] gr1 gr1
   !> @param[out] gr2 gr2
   pure subroutine moment_atom_native(state, order, valid, gc0, gc1, gc2, gr0, gr1, gr2)
      !> Native moment state
      type(moment_state_type), intent(in) :: state
      !> Highest total derivative order
      integer, intent(in) :: order
      !> Successful finite evaluation
      logical, intent(out) :: valid
      !> Gc0
      real(wp), intent(out), optional :: gc0(:, :)
      !> Gc1
      real(wp), intent(out), optional :: gc1(:, :, :)
      !> Gc2
      real(wp), intent(out), optional :: gc2(:, :, :, :)
      !> Gr0
      real(wp), intent(out), optional :: gr0(:)
      !> Gr1
      real(wp), intent(out), optional :: gr1(:, :)
      !> Gr2
      real(wp), intent(out), optional :: gr2(:, :, :)
      real(wp) :: d, u(3), t(3, 3), jac(3, 3), z(3, 3), dt(3, 3, 3), dz(3, 3, 3)
      real(wp) :: du(3, 3), du2(3, 3, 3), b(3, 3), bs(3, 3), unit(3, 3), lift(3, 3)
      real(wp) :: wa, h(3), da(3), g(3), dw(3), gmat(3, 3), d2w(3, 3)
      integer :: a, i, k, l, s
      logical :: radius
      valid = .false.
      if (.not. state%valid .or. order < 0 .or. order > 2) return
      if (order >= 1 .and. .not. state%sensitivity_valid) return
      radius = present(gr0) .or. present(gr1) .or. present(gr2)
      d = sqrt(dot_product(state%y, state%y) + state%eps**2)
      u = state%p*state%y/d
      unit = 0.0_wp
      do k = 1, 3
         unit(k, k) = 1.0_wp
      end do
      du = 0.0_wp
      du2 = 0.0_wp
      if (order >= 1) then
         call moment_response_matrices(state%y, state%covariance, state%p, state%tau, 1.0_wp/d, t, jac)
         ! du(i, k) = du_i/dx_k
         du = matmul(t, state%sensitivity)
      end if
      if (order >= 2) then
         call third_response(state, t, z, dt, dz, valid)
         if (.not. valid) return
         valid = .false.
         ! du2(i, k, l) = d2u_i/dx_k dx_l, symmetrized in k and l
         do l = 1, 3
            do k = 1, 3
               du2(:, k, l) = matmul(dt(:, :, l), state%sensitivity(:, k)) + 0.5_wp*matmul(t, dz(:, k, l))
            end do
         end do
         do l = 1, 3
            do k = 1, l
               du2(:, k, l) = 0.5_wp*(du2(:, k, l) + du2(:, l, k))
               du2(:, l, k) = du2(:, k, l)
            end do
         end do
      end if
      b = unit + transpose(du)
      bs = matmul(b, unit - state%sensitivity)
      lift = unit + du
      do a = 1, state%n
         wa = state%w(a)/state%total
         h = state%v(:, a) + u
         da = state%v(:, a) - state%m
         g = matmul(b, da)
         dw = -2.0_wp*wa/state%tau*g
         d2w = 0.0_wp
         if (order >= 2) then
            do l = 1, 3
               do k = 1, 3
                  gmat(k, l) = bs(k, l)
                  do i = 1, 3
                     gmat(k, l) = gmat(k, l) + du2(i, k, l)*da(i)
                  end do
               end do
            end do
            gmat = 0.5_wp*(gmat + transpose(gmat))
            do l = 1, 3
               d2w(:, l) = 4.0_wp*wa/state%tau**2*g*g(l) - 2.0_wp*wa/state%tau*gmat(:, l)
            end do
         end if
         if (present(gc0)) gc0(:, a) = -2.0_wp*wa*h
         if (present(gc1) .and. order >= 1) then
            do s = 1, 3
               gc1(:, s, a) = -2.0_wp*(dw*h(s) + wa*lift(s, :))
            end do
         end if
         if (present(gc2) .and. order >= 2) then
            do s = 1, 3
               do l = 1, 3
                  gc2(:, l, s, a) = -2.0_wp*(d2w(:, l)*h(s) + dw*lift(s, l) + dw(l)*lift(s, :) &
                                             + wa*du2(s, :, l))
               end do
            end do
         end if
         if (radius) then
            if (present(gr0)) gr0(a) = wa*state%slope(a)
            if (present(gr1) .and. order >= 1) gr1(:, a) = dw*state%slope(a)
            if (present(gr2) .and. order >= 2) gr2(:, :, a) = d2w*state%slope(a)
         end if
      end do
      valid = .true.
   end subroutine moment_atom_native

end module moist_cavity_drop_lsf_sses_native
