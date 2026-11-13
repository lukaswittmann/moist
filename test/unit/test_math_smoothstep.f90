!> Test suite for smoothstep functions in moist_math_smoothstep
!>
!> Polynomial type:
!>
!> - polynomial_plateaus: value exactly 0/1, derivatives exactly 0 at and beyond the edges
!> - polynomial_range_monotone: value in [0, 1] and monotone up to both edges
!> - polynomial_continuity_class: orders 1..n vanish at both edges, order n+1 does not
!> - polynomial_leading_power: leading power t^(n+1) at the lower edge
!> - polynomial_point_symmetry: point symmetry about the band center
!> - polynomial_fd_derivatives: derivatives 1..4 against central differences
!> - polynomial_default_order: absent order equals order 4
!> - polynomial_reference_values: exact rationals for S_n(1/4) and S_n'(1/2), pins the table
!>
!> Erf and sigmoid types:
!>
!> - erf_pinned_values, sigmoid_pinned_values: values pinned from the pre-move `switching.f90`
!> - sigmoid_endpoint: no divide by zero next to the lower edge
!>
!> Nuclear chain rule of the base type:
!>
!> - f1_ra_chain_rule: `f1_rA` against central differences of `f0` over a synthetic
!>   argument, and the `active` column mask
!>
!> Error stops, not testable in process:
!>
!> - Polynomial order outside 1..5
!> - Polynomial evaluation before new_smoothstep_polynomial
!> - f3 or f4 requested from the erf or sigmoid type
module test_math_smoothstep
   use, intrinsic :: ieee_exceptions, only: ieee_get_flag, ieee_set_flag, ieee_divide_by_zero
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_math_smoothstep, only: moist_math_smoothstep_type, &
      & moist_math_smoothstep_polynomial_type, new_smoothstep_polynomial, &
      & moist_math_smoothstep_erf_type, new_smoothstep_erf, &
      & moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   implicit none(type, external)
   private

   public :: collect_math_smoothstep

   !> Band ends, non-unit width so the 1/width chain rule is exercised
   real(wp), parameter :: band_a = 0.3_wp, band_b = 1.7_wp
   !> Orders of the polynomial smoothstep
   integer, parameter :: min_order = 1, max_order = 5
   !> FD step in units of the band width
   real(wp), parameter :: h_t = 1.0e-2_wp
   !> FD tolerance in t-units, absolute with 10x relative fallback; max abs 8e-9 at order 5, rel 1e-11
   real(wp), parameter :: fd_thr = 1.0e-9_wp
   !> Interior points in t, at least 4 FD steps from either edge
   real(wp), parameter :: t_interior(7) = [0.05_wp, 0.2_wp, 0.37_wp, 0.5_wp, 0.71_wp, 0.9_wp, 0.95_wp]
   !> Value symmetry tolerance, absolute on an O(1) sum; measured 1.1e-16
   real(wp), parameter :: sym_thr = 1.0e-14_wp
   !> Derivative symmetry tolerance in t-units, 10x relative fallback; measured abs 9e-13, rel 7e-15
   real(wp), parameter :: dsym_thr = 5.0e-12_wp
   !> Pinned erf/sigmoid values, relative tolerance
   real(wp), parameter :: pin_thr = 1.0e-13_wp
   !> Central first-derivative weights, offsets 1..4, O(h^8)
   real(wp), parameter :: w8(4) = [4.0_wp/5.0_wp, -1.0_wp/5.0_wp, 4.0_wp/105.0_wp, -1.0_wp/280.0_wp]
   !> Centre count of the synthetic chain-rule argument
   integer, parameter :: ncent = 4

contains

   !> Collect all smoothstep tests
   subroutine collect_math_smoothstep(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("polynomial_plateaus", test_polynomial_plateaus), &
                  new_unittest("polynomial_range_monotone", test_polynomial_range_monotone), &
                  new_unittest("polynomial_continuity_class", test_polynomial_continuity_class), &
                  new_unittest("polynomial_leading_power", test_polynomial_leading_power), &
                  new_unittest("polynomial_point_symmetry", test_polynomial_point_symmetry), &
                  new_unittest("polynomial_fd_derivatives", test_polynomial_fd_derivatives), &
                  new_unittest("polynomial_default_order", test_polynomial_default_order), &
                  new_unittest("polynomial_reference_values", test_polynomial_reference_values), &
                  new_unittest("erf_pinned_values", test_erf_pinned_values), &
                  new_unittest("sigmoid_pinned_values", test_sigmoid_pinned_values), &
                  new_unittest("sigmoid_endpoint", test_sigmoid_endpoint), &
                  new_unittest("f1_ra_chain_rule", test_f1_rA_chain_rule) &
                  ]
   end subroutine collect_math_smoothstep

   !> k-th derivative (k = 0..4) of a smoothstep at x
   function deriv(sw, x, k) result(val)
      !> Smoothstep instance
      class(moist_math_smoothstep_type), intent(in) :: sw
      !> Input value
      real(wp), intent(in) :: x
      !> Derivative order, 0..4
      integer, intent(in) :: k
      !> Derivative value
      real(wp) :: val
      real(wp) :: f0, f1, f2, f3, f4

      call sw%eval(x, f0, f1, f2, f3, f4)
      select case (k)
      case (0)
         val = f0
      case (1)
         val = f1
      case (2)
         val = f2
      case (3)
         val = f3
      case default
         val = f4
      end select
   end function deriv

   !> `x^k` as a multiply chain; `x**k` traps under -ffpe-trap on arm64
   pure function wpow(x, k) result(p)
      !> Base
      real(wp), intent(in) :: x
      !> Exponent, >= 0
      integer, intent(in) :: k
      !> Power
      real(wp) :: p
      integer :: j

      p = 1.0_wp
      do j = 1, k
         p = p*x
      end do
   end function wpow

   !> Absolute tolerance with 10x relative fallback
   pure function within(val, ref, thr) result(ok)
      !> Computed value
      real(wp), intent(in) :: val
      !> Reference value
      real(wp), intent(in) :: ref
      !> Absolute tolerance
      real(wp), intent(in) :: thr
      !> Whether val matches ref
      logical :: ok

      ok = abs(val - ref) <= thr .or. abs(val - ref) <= 10.0_wp*thr*abs(ref)
   end function within

   !> Polynomial smoothstep on the test band, rising or falling
   subroutine make_polynomial(sw, order, rising)
      !> Smoothstep instance
      type(moist_math_smoothstep_polynomial_type), intent(out) :: sw
      !> Continuity order n
      integer, intent(in) :: order
      !> Rising (from < to) or falling orientation
      logical, intent(in) :: rising

      if (rising) then
         call new_smoothstep_polynomial(sw, band_a, band_b, order)
      else
         call new_smoothstep_polynomial(sw, band_b, band_a, order)
      end if
   end subroutine make_polynomial

   !> Value exactly 0/1 and derivatives exactly 0 at the edges and beyond
   subroutine test_polynomial_plateaus(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_smoothstep_polynomial_type) :: sw
      real(wp) :: x_lo(3), x_hi(3), f(0:4), v_lo, v_hi
      character(len=96) :: msg
      integer :: n, i, idir
      logical :: rising

      x_lo = [band_a, band_a - 0.25_wp, -1.0e6_wp]
      x_hi = [band_b, band_b + 0.25_wp, 1.0e6_wp]
      do idir = 1, 2
         rising = idir == 1
         v_lo = merge(0.0_wp, 1.0_wp, rising)
         v_hi = 1.0_wp - v_lo
         do n = min_order, max_order
            call make_polynomial(sw, n, rising)
            do i = 1, size(x_lo)
               call sw%eval(x_lo(i), f(0), f(1), f(2), f(3), f(4))
               write (msg, "(a,i0,a,l1,a,i0)") "lower plateau, order ", n, " rising ", rising, " point ", i
               call check(error, f(0) == v_lo .and. all(f(1:4) == 0.0_wp), msg)
               if (allocated(error)) return
               call sw%eval(x_hi(i), f(0), f(1), f(2), f(3), f(4))
               write (msg, "(a,i0,a,l1,a,i0)") "upper plateau, order ", n, " rising ", rising, " point ", i
               call check(error, f(0) == v_hi .and. all(f(1:4) == 0.0_wp), msg)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_polynomial_plateaus

   !> Value in [0, 1] and first derivative of the orientation's sign, edges included
   !>
   !> - Points placed by offset from each band end, so `hi - x` is not rounded through `x - lo`
   subroutine test_polynomial_range_monotone(error)
      type(error_type), allocatable, intent(out) :: error
      !> Offsets from a band end in t-units, edge layer through the band center
      real(wp), parameter :: offsets(8) = [1.0e-10_wp, 1.0e-7_wp, 1.0e-5_wp, 1.0e-3_wp, &
                                           0.05_wp, 0.2_wp, 0.37_wp, 0.5_wp]
      type(moist_math_smoothstep_polynomial_type) :: sw
      real(wp) :: width, x, f0, f1, direction
      character(len=128) :: msg
      integer :: n, i, iedge, idir
      logical :: rising

      width = band_b - band_a
      do idir = 1, 2
         rising = idir == 1
         direction = merge(1.0_wp, -1.0_wp, rising)
         do n = min_order, max_order
            call make_polynomial(sw, n, rising)
            do iedge = 1, 2
               do i = 1, size(offsets)
                  if (iedge == 1) then
                     x = band_a + offsets(i)*width
                  else
                     x = band_b - offsets(i)*width
                  end if
                  call sw%eval(x, f0, f1)
                  write (msg, "(a,i0,a,i0,a,es8.1,a,l1,a,es10.3,a,es10.3)") "order ", n, " edge ", iedge, &
                     & " offset ", offsets(i), " rising ", rising, " f0 ", f0, " f1 ", f1
                  call check(error, f0 >= 0.0_wp .and. f0 <= 1.0_wp .and. direction*f1 >= 0.0_wp, msg)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_polynomial_range_monotone

   !> One-sided limits at both edges: orders 1..n vanish, order n+1 does not
   subroutine test_polynomial_continuity_class(error)
      type(error_type), allocatable, intent(out) :: error
      !> Edge offset in t-units
      real(wp), parameter :: delta_t = 1.0e-10_wp
      !> Bound on |f_k| width^k, k <= n, at delta_t from the edge; measured 1.5e-6 (n = k = 4) at both edges
      real(wp), parameter :: vanish_thr = 1.0e-5_wp
      type(moist_math_smoothstep_polynomial_type) :: sw
      real(wp) :: width, x_edge(2), scaled
      character(len=96) :: msg
      integer :: n, k, iedge, idir
      logical :: rising

      width = band_b - band_a
      x_edge = [band_a + delta_t*width, band_b - delta_t*width]
      do idir = 1, 2
         rising = idir == 1
         do n = min_order, max_order
            call make_polynomial(sw, n, rising)
            do iedge = 1, 2
               do k = 1, min(n + 1, 4)
                  scaled = abs(deriv(sw, x_edge(iedge), k))*wpow(width, k)
                  write (msg, "(a,i0,a,i0,a,i0,a,l1,a,es10.3)") "edge ", iedge, " order ", n, " deriv ", k, &
                     & " rising ", rising, " |f_k| w^k ", scaled
                  if (k <= n) then
                     call check(error, scaled <= vanish_thr, msg)
                  else
                     call check(error, scaled >= 1.0_wp, msg)
                  end if
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_polynomial_continuity_class

   !> Leading power t^(n+1) of the value at the lower edge
   !>
   !> - Only check that order n+1 is nonzero for n = 4, 5, beyond the order-4 continuity test
   !> - Upper edge follows by point symmetry
   subroutine test_polynomial_leading_power(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: delta_t = 1.0e-5_wp
      type(moist_math_smoothstep_polynomial_type) :: sw
      real(wp) :: width, power
      character(len=96) :: msg
      integer :: n

      width = band_b - band_a
      do n = min_order, max_order
         call make_polynomial(sw, n, .true.)
         power = log(sw%f0(band_a + 2.0_wp*delta_t*width)/sw%f0(band_a + delta_t*width))/log(2.0_wp)
         write (msg, "(a,i0,a,f12.8)") "order ", n, " leading power ", power
         call check(error, abs(power - real(n + 1, wp)) < 1.0e-3_wp, msg)
         if (allocated(error)) return
      end do
   end subroutine test_polynomial_leading_power

   !> Point symmetry about the band center
   !>
   !> - `S(lo + d) + S(hi - d) = 1`, `f_k(lo + d) = (-1)^(k+1) f_k(hi - d)`
   subroutine test_polynomial_point_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_smoothstep_polynomial_type) :: sw
      real(wp) :: width, d, lhs, rhs, sgn, diff
      character(len=96) :: msg
      integer :: n, k, i

      width = band_b - band_a
      do n = min_order, max_order
         call make_polynomial(sw, n, .true.)
         do i = 1, size(t_interior)
            d = t_interior(i)*width
            write (msg, "(a,i0,a,i0)") "value symmetry, order ", n, " point ", i
            call check(error, abs(sw%f0(band_a + d) + sw%f0(band_b - d) - 1.0_wp) < sym_thr, msg)
            if (allocated(error)) return
            sgn = 1.0_wp
            do k = 1, 4
               lhs = deriv(sw, band_a + d, k)*wpow(width, k)
               rhs = sgn*deriv(sw, band_b - d, k)*wpow(width, k)
               diff = abs(lhs - rhs)
               write (msg, "(a,i0,a,i0,a,i0,a,es10.3)") "deriv symmetry, order ", n, " deriv ", k, &
                  & " point ", i, " diff ", diff
               call check(error, diff <= dsym_thr .or. diff <= 10.0_wp*dsym_thr*abs(lhs), msg)
               if (allocated(error)) return
               sgn = -sgn
            end do
         end do
      end do
   end subroutine test_polynomial_point_symmetry

   !> Derivatives 1..4 against an 8th-order central difference of the next lower one
   subroutine test_polynomial_fd_derivatives(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_smoothstep_polynomial_type) :: sw
      real(wp) :: width, h, x, fd, an, diff
      character(len=128) :: msg
      integer :: n, k, i, j, idir
      logical :: rising

      width = band_b - band_a
      h = h_t*width
      do idir = 1, 2
         rising = idir == 1
         do n = min_order, max_order
            call make_polynomial(sw, n, rising)
            do i = 1, size(t_interior)
               x = band_a + t_interior(i)*width
               do k = 1, 4
                  fd = 0.0_wp
                  do j = 1, 4
                     fd = fd + w8(j)*(deriv(sw, x + j*h, k - 1) - deriv(sw, x - j*h, k - 1))
                  end do
                  fd = fd/h
                  an = deriv(sw, x, k)
                  diff = abs(fd - an)*wpow(width, k)
                  write (msg, "(a,i0,a,i0,a,f5.2,a,l1,a,es10.3,a,es10.3)") "order ", n, " deriv ", k, " t ", &
                     & t_interior(i), " rising ", rising, " diff ", diff, " ref ", an*wpow(width, k)
                  call check(error, diff <= fd_thr .or. diff <= 10.0_wp*fd_thr*abs(an)*wpow(width, k), msg)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_polynomial_fd_derivatives

   !> Absent order equals order 4, bit for bit
   subroutine test_polynomial_default_order(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_smoothstep_polynomial_type) :: sw_default, sw_four
      real(wp) :: x
      integer :: i, k

      call new_smoothstep_polynomial(sw_default, band_a, band_b)
      call new_smoothstep_polynomial(sw_four, band_a, band_b, 4)
      do i = 1, size(t_interior)
         x = band_a + t_interior(i)*(band_b - band_a)
         do k = 0, 4
            call check(error, deriv(sw_default, x, k) == deriv(sw_four, x, k), "default order differs from 4")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_polynomial_default_order

   !> Values and midpoint slopes against exact rationals, both orientations
   !>
   !> - Dyadic band, t = 1/4, 1/2, 3/4 exact in binary
   !> - References from the normalized integral of t^n (1-t)^n, independent of the coefficient table
   !> - S_n(1/4) = 5/32, 53/512, 289/4096, 6413/131072, 35995/1048576
   !> - S_n'(1/2) = (2n+1)!/((n!)^2 4^n) = 3/2, 15/8, 35/16, 315/128, 693/256
   !> - t = 3/4 through the mirrored branch, S_n(3/4) = 1 - S_n(1/4)
   !> - Negative control: no order matches the reference of another order
   subroutine test_polynomial_reference_values(error)
      type(error_type), allocatable, intent(out) :: error
      !> Dyadic band, width 2 so the 1/width scaling enters
      real(wp), parameter :: ref_lo = 0.25_wp, ref_hi = 2.25_wp
      !> S_n(1/4), n = 1..5
      real(wp), parameter :: s_quarter(max_order) = [5.0_wp/32.0_wp, 53.0_wp/512.0_wp, &
                                                     289.0_wp/4096.0_wp, 6413.0_wp/131072.0_wp, &
                                                     35995.0_wp/1048576.0_wp]
      !> S_n'(1/2), n = 1..5
      real(wp), parameter :: ds_half(max_order) = [3.0_wp/2.0_wp, 15.0_wp/8.0_wp, 35.0_wp/16.0_wp, &
                                                   315.0_wp/128.0_wp, 693.0_wp/256.0_wp]
      !> Absolute tolerance, 10x relative fallback
      real(wp), parameter :: ref_thr = 1.0e-15_wp
      type(moist_math_smoothstep_polynomial_type) :: sw
      real(wp) :: width, x_q, x_mid, x_3q, f_q, f_mid, f_3q, f1_mid, v_q, v_3q, d_mid, direction
      character(len=128) :: msg
      integer :: n, m, idir
      logical :: rising

      width = ref_hi - ref_lo
      x_q = ref_lo + 0.25_wp*width
      x_mid = ref_lo + 0.5_wp*width
      x_3q = ref_lo + 0.75_wp*width
      do idir = 1, 2
         rising = idir == 1
         direction = merge(1.0_wp, -1.0_wp, rising)
         do n = min_order, max_order
            if (rising) then
               call new_smoothstep_polynomial(sw, ref_lo, ref_hi, n)
            else
               call new_smoothstep_polynomial(sw, ref_hi, ref_lo, n)
            end if
            v_q = merge(s_quarter(n), 1.0_wp - s_quarter(n), rising)
            v_3q = 1.0_wp - v_q
            d_mid = direction*ds_half(n)/width
            call sw%eval(x_q, f_q)
            call sw%eval(x_mid, f_mid, f1_mid)
            call sw%eval(x_3q, f_3q)

            write (msg, "(a,i0,a,l1,a,es23.16,a,es23.16)") "S(1/4), order ", n, " rising ", rising, &
               & " f0 ", f_q, " ref ", v_q
            call check(error, within(f_q, v_q, ref_thr), msg)
            if (allocated(error)) return
            write (msg, "(a,i0,a,l1,a,es23.16,a,es23.16)") "S(3/4), order ", n, " rising ", rising, &
               & " f0 ", f_3q, " ref ", v_3q
            call check(error, within(f_3q, v_3q, ref_thr), msg)
            if (allocated(error)) return
            write (msg, "(a,i0,a,l1,a,es23.16,a,es23.16)") "S'(1/2), order ", n, " rising ", rising, &
               & " f1 ", f1_mid, " ref ", d_mid
            call check(error, within(f1_mid, d_mid, ref_thr), msg)
            if (allocated(error)) return

            if (.not. rising) cycle
            do m = min_order, max_order
               if (m == n) cycle
               write (msg, "(a,i0,a,i0)") "negative control: order ", n, " matches reference of order ", m
               call check(error, .not. within(f_q, s_quarter(m), ref_thr) &
                  & .and. .not. within(f1_mid, ds_half(m)/width, ref_thr), msg)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_polynomial_reference_values

   !> Compare f0..f2 against pinned values
   subroutine check_pinned(error, sw, x, ref, label)
      type(error_type), allocatable, intent(out) :: error
      !> Smoothstep instance
      class(moist_math_smoothstep_type), intent(in) :: sw
      !> Input value
      real(wp), intent(in) :: x
      !> Pinned f0, f1, f2
      real(wp), intent(in) :: ref(3)
      !> Case label
      character(len=*), intent(in) :: label
      real(wp) :: f(3)
      integer :: k

      call sw%eval(x, f(1), f(2), f(3))
      do k = 1, 3
         if (abs(f(k) - ref(k)) > pin_thr*max(abs(ref(k)), 1.0e-300_wp)) then
            call test_failed(error, label//": f0..f2 differ from pre-move switching.f90")
            return
         end if
      end do
   end subroutine check_pinned

   !> Erf type unchanged by the move, both orientations, plateaus
   subroutine test_erf_pinned_values(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_smoothstep_erf_type) :: e1, e2

      call new_smoothstep_erf(e1, -0.5_wp, 0.5_wp)
      call new_smoothstep_erf(e2, 2.0_wp, 0.5_wp, 0.7_wp)
      call check_pinned(error, e1, -0.15_wp, [3.10979211464817752e-2_wp, 1.04389116265920179_wp, &
         & 2.50883470675889022e1_wp], "erf rising 1")
      if (allocated(error)) return
      call check_pinned(error, e1, 0.3_wp, [9.99999943136371749e-1_wp, 1.17065212897972775e-5_wp, &
         & -2.22370106154632860e-3_wp], "erf rising 2")
      if (allocated(error)) return
      call check_pinned(error, e1, 0.7_wp, [1.0_wp, 0.0_wp, 0.0_wp], "erf rising plateau")
      if (allocated(error)) return
      call check_pinned(error, e2, 0.85_wp, [9.98418326555498448e-1_wp, -6.78312268828850096e-2_wp, &
         & -2.30679470544622722_wp], "erf falling 1")
      if (allocated(error)) return
      call check_pinned(error, e2, 1.3_wp, [3.95441161470311942e-1_wp, -2.06086045970365683_wp, &
         & 1.82318512416276191_wp], "erf falling 2")
      if (allocated(error)) return
      call check_pinned(error, e2, 1.7_wp, [1.02687868070661104e-4_wp, -7.11443437396035855e-3_wp, &
         & 4.19051954830192308e-1_wp], "erf falling 3")
   end subroutine test_erf_pinned_values

   !> Sigmoid type unchanged by the move, both orientations, plateaus
   subroutine test_sigmoid_pinned_values(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_smoothstep_sigmoid_type) :: s1, s2

      call new_smoothstep_sigmoid(s1, 0.0_wp, 1.0_wp)
      call new_smoothstep_sigmoid(s2, 1.5_wp, 0.25_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=1.3_wp, a_lo=0.9_wp)
      call check_pinned(error, s1, 0.05_wp, [0.0_wp, 0.0_wp, 0.0_wp], "sigmoid rising lower saturation")
      if (allocated(error)) return
      call check_pinned(error, s1, 0.3_wp, [1.15019392821386643e-4_wp, 9.18956493757581289e-3_wp, &
         & 6.51807273246784447e-1_wp], "sigmoid rising 1")
      if (allocated(error)) return
      call check_pinned(error, s1, 0.75_wp, [9.99999334164296427e-1_wp, 8.83834656549667316e-5_wp, &
         & -1.07219743266472201e-2_wp], "sigmoid rising 2")
      if (allocated(error)) return
      call check_pinned(error, s1, 1.15_wp, [1.0_wp, 0.0_wp, 0.0_wp], "sigmoid rising plateau")
      if (allocated(error)) return
      call check_pinned(error, s2, 0.55_wp, [9.77271597446259355e-1_wp, -5.91102248137361896e-1_wp, &
         & -1.08400091289437412e1_wp], "sigmoid falling 1")
      if (allocated(error)) return
      call check_pinned(error, s2, 1.0_wp, [1.70541790808148241e-1_wp, -1.18243242314795860_wp, &
         & 5.11389089415691611_wp], "sigmoid falling 2")
      if (allocated(error)) return
      call check_pinned(error, s2, 1.4_wp, [1.56454857869192665e-5_wp, -1.52823798204837621e-3_wp, &
         & 1.22121184095400168e-1_wp], "sigmoid falling 3")
   end subroutine test_sigmoid_pinned_values

   !> Next to the lower edge the sigmoid stays finite without a divide by zero
   subroutine test_sigmoid_endpoint(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_smoothstep_sigmoid_type) :: sw
      real(wp) :: value
      logical :: raised, previous

      call new_smoothstep_sigmoid(sw, 0.0_wp, 1.0_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=0.8_wp, a_lo=1.6_wp)
      call ieee_get_flag(ieee_divide_by_zero, previous)
      call ieee_set_flag(ieee_divide_by_zero, .false.)
      value = sw%f0(1.0e-20_wp)
      call ieee_get_flag(ieee_divide_by_zero, raised)
      call ieee_set_flag(ieee_divide_by_zero, previous)
      call check(error,.not. raised .and. value == 0.0_wp, "finite endpoint switching without divide by zero")
   end subroutine test_sigmoid_endpoint

   !> Synthetic smoothstep argument: linear plus isotropic quadratic in the centres
   pure subroutine synthetic_arg(centers, offset, x0, x1)
      !> Centre coordinates (3, ncent)
      real(wp), intent(in) :: centers(3, ncent)
      !> Constant shift of the argument
      real(wp), intent(in) :: offset
      !> Argument value
      real(wp), intent(out) :: x0
      !> Argument gradient w.r.t. the centres (3, ncent)
      real(wp), intent(out), optional :: x1(3, ncent)
      !> Curvature of the quadratic part
      real(wp), parameter :: curv = 0.1_wp
      real(wp) :: lin
      integer :: ic, ia

      x0 = offset
      do ia = 1, ncent
         do ic = 1, 3
            lin = 0.05_wp*real(ic - 2, wp) + 0.02_wp*real(ia, wp)
            x0 = x0 + lin*centers(ic, ia) + 0.5_wp*curv*centers(ic, ia)*centers(ic, ia)
            if (present(x1)) x1(ic, ia) = lin + curv*centers(ic, ia)
         end do
      end do
   end subroutine synthetic_arg

   !> `f1_rA` against central differences of `f0` over the synthetic argument
   !>
   !> - `active`: listed columns bit-identical to the full result, the others exactly 0
   !> - against `fd_thr`: measured abs 1.4e-13, rel 3.5e-12
   subroutine test_f1_rA_chain_rule(error)
      type(error_type), allocatable, intent(out) :: error
      !> Argument values inside the erf transition
      real(wp), parameter :: x_target(3) = [-0.2_wp, -0.1_wp, 0.0_wp]
      !> Listed centres, unsorted
      integer, parameter :: active(2) = [3, 1]
      !> Centre step
      real(wp), parameter :: h = 1.0e-3_wp
      type(moist_math_smoothstep_erf_type) :: sw
      real(wp) :: centers(3, ncent), shifted(3, ncent), x1(3, ncent)
      real(wp) :: grad(3, ncent), masked(3, ncent)
      real(wp) :: x0, x_plus, x_minus, offset, fd
      character(len=128) :: msg
      integer :: it, ia, ic, j

      call new_smoothstep_erf(sw, -0.5_wp, 0.5_wp)
      do ia = 1, ncent
         do ic = 1, 3
            centers(ic, ia) = 0.3_wp*sin(real(ic + 3*ia, wp))
         end do
      end do

      do it = 1, size(x_target)
         ! Shift that puts the argument on the target at the reference centres
         call synthetic_arg(centers, 0.0_wp, x0)
         offset = x_target(it) - x0
         call synthetic_arg(centers, offset, x0, x1)
         grad = sw%f1_rA(x0, x1)

         ! Negative control: a vanishing gradient would pass every comparison below
         call check(error, maxval(abs(grad)) > 1.0e-3_wp, "chain-rule gradient vanishes at the test point")
         if (allocated(error)) return

         do ia = 1, ncent
            do ic = 1, 3
               fd = 0.0_wp
               do j = 1, 4
                  shifted = centers
                  shifted(ic, ia) = centers(ic, ia) + real(j, wp)*h
                  call synthetic_arg(shifted, offset, x_plus)
                  shifted(ic, ia) = centers(ic, ia) - real(j, wp)*h
                  call synthetic_arg(shifted, offset, x_minus)
                  fd = fd + w8(j)*(sw%f0(x_plus) - sw%f0(x_minus))
               end do
               fd = fd/h
               write (msg, "(a,i0,a,i0,a,f5.2,a,es10.3,a,es10.3)") "centre ", ia, " axis ", ic, " x0 ", &
                  & x_target(it), " diff ", abs(grad(ic, ia) - fd), " ref ", fd
               call check(error, within(grad(ic, ia), fd, fd_thr), msg)
               if (allocated(error)) return
            end do
         end do

         masked = sw%f1_rA(x0, x1, active)
         call check(error, all(masked(:, active) == grad(:, active)), "active columns differ from the full gradient")
         if (allocated(error)) return
         call check(error, all(masked(:, [2, 4]) == 0.0_wp), "columns outside active are not zero")
         if (allocated(error)) return
      end do
   end subroutine test_f1_rA_chain_rule

end module test_math_smoothstep
