!> Implementation of various smoothing Kernels
!> References:
!> https://pysph.readthedocs.io/en/main/reference/kernels.html
!> https://ludwigboess.github.io/SPHKernels.jl/stable/kernels/
module moist_math_smoothing_kernels

   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi

   implicit none (type, external)
   private

   public :: wendland_kernel_type

   !> Abstract base type for smoothing kernels
   type, abstract :: smoothing_kernel_type
      real(wp) :: h              !< Smoothing length
      integer :: dimension       !< Spatial dimension (1, 2, or 3)
      integer :: order           !< Kernel order (2=C2, 4=C4, 6=C6, etc.)
   contains
      !> Initialize the kernel
      procedure(init_interface), deferred :: init
      !> Evaluate kernel at distance r
      procedure(f0_interface), deferred :: f0
      !> Evaluate first derivative at distance r
      procedure(f1_interface), deferred :: f1
   end type smoothing_kernel_type

   !> Wendland smoothing kernel type (supports C2, C4, C6, etc.)
   type, extends(smoothing_kernel_type) :: wendland_kernel_type
      real(wp) :: prefactor      !< Precomputed normalization factor (includes h^n)
      procedure(wendland_compute_interface), pointer, nopass :: compute => null()
      procedure(wendland_compute_interface), pointer, nopass :: compute_deriv => null()
   contains
      procedure :: init => wendland_init
      procedure :: f0 => wendland_f0
      procedure :: f1 => wendland_f1
      procedure :: gradient => wendland_gradient
      procedure :: gradient_h => wendland_gradient_h
   end type wendland_kernel_type

   abstract interface
      !> Function signature for dimension-specific kernel computation
      pure function wendland_compute_interface(q) result(val)
         import :: wp
         implicit none (type, external)
         real(wp), intent(in) :: q
         real(wp) :: val
      end function wendland_compute_interface
   end interface

   abstract interface
      !> Initialize kernel interface
      !>
      !> @param[inout] self      Kernel instance
      !> @param[in]    order     Kernel order
      !> @param[in]    dimension Spatial dimension
      !> @param[in]    h         Smoothing length
      !> @param[out]   error     Unsupported order/dimension, or invalid `h`
      subroutine init_interface(self, order, dimension, h, error)
         import :: smoothing_kernel_type, wp, error_type
         implicit none (type, external)
         class(smoothing_kernel_type), intent(inout) :: self
         integer, intent(in) :: order
         integer, intent(in) :: dimension
         real(wp), intent(in) :: h
         type(error_type), allocatable, intent(out) :: error
      end subroutine init_interface

      !> Kernel evaluation interface
      pure function f0_interface(self, r) result(kernel_val)
         import :: smoothing_kernel_type, wp
         implicit none (type, external)
         class(smoothing_kernel_type), intent(in) :: self
         real(wp), intent(in) :: r
         real(wp) :: kernel_val
      end function f0_interface

      !> Kernel derivative interface (dW/dr)
      pure function f1_interface(self, r) result(derivative)
         import :: smoothing_kernel_type, wp
         implicit none (type, external)
         class(smoothing_kernel_type), intent(in) :: self
         real(wp), intent(in) :: r
         real(wp) :: derivative
      end function f1_interface
   end interface

contains

   !> Wendland C2 kernel for 1D: (1-q/2)^3 (1.5q+1), q < 2
   pure function wendland_c2_1d(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         val = tmp*tmp*tmp*(1.5_wp*q + 1.0_wp)
      else
         val = 0.0_wp
      end if
   end function wendland_c2_1d

   !> Wendland C2 kernel for 2D/3D: (1-q/2)^4 (2q+1), q < 2
   pure function wendland_c2_23d(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         val = tmp*tmp*tmp*tmp*(2.0_wp*q + 1.0_wp)
      else
         val = 0.0_wp
      end if
   end function wendland_c2_23d

   !> Derivative of Wendland C2 kernel for 1D: dW/dq = -3q(1-q/2)^2 , q < 2
   pure function wendland_c2_1d_deriv(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         val = -3.0_wp*q*tmp*tmp
      else
         val = 0.0_wp
      end if
   end function wendland_c2_1d_deriv

   !> Derivative of Wendland C2 kernel for 2D/3D: dW/dq = -5q(1-q/2)^3 , q < 2
   pure function wendland_c2_23d_deriv(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         val = -5.0_wp*q*tmp*tmp*tmp
      else
         val = 0.0_wp
      end if
   end function wendland_c2_23d_deriv

   !> Wendland C4 kernel for 1D: (1-q/2)^5 (2q^2 +2.5q+1), q < 2
   pure function wendland_c4_1d(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         val = tmp*tmp*tmp*tmp*tmp*(2.0_wp*q*q + 2.5_wp*q + 1.0_wp)
      else
         val = 0.0_wp
      end if
   end function wendland_c4_1d

   !> Wendland C4 kernel for 2D/3D: (1-q/2)^6 (35q^2 /12+3q+1), q < 2
   pure function wendland_c4_23d(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         val = tmp**6*((35.0_wp/12.0_wp)*q*q + 3.0_wp*q + 1.0_wp)
      else
         val = 0.0_wp
      end if
   end function wendland_c4_23d

   !> Derivative of Wendland C4 kernel for 1D: dW/dq = -3.5q(2q+1)(1-q/2)^4 , q < 2
   pure function wendland_c4_1d_deriv(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         val = -3.5_wp*q*(2.0_wp*q + 1.0_wp)*tmp**4
      else
         val = 0.0_wp
      end if
   end function wendland_c4_1d_deriv

   !> Derivative of Wendland C4 kernel for 2D/3D: dW/dq = -(14/3)q(1+2.5q)(1-q/2)^5 , q < 2
   pure function wendland_c4_23d_deriv(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         val = -(14.0_wp/3.0_wp)*q*(1.0_wp + 2.5_wp*q)*tmp**5
      else
         val = 0.0_wp
      end if
   end function wendland_c4_23d_deriv

   !> Wendland C6 kernel for 1D: (1-q/2)^7 (21q^3 /8+19q^2 /4+3.5q+1), q < 2
   pure function wendland_c6_1d(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         val = tmp**7*(2.625_wp*q**3 + 4.75_wp*q*q + 3.5_wp*q + 1.0_wp)
      else
         val = 0.0_wp
      end if
   end function wendland_c6_1d

   !> Wendland C6 kernel for 2D/3D: (1-q/2)^8 (4q^3 +6.25q^2 +4q+1), q < 2
   pure function wendland_c6_23d(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         val = tmp**8*(4.0_wp*q**3 + 6.25_wp*q*q + 4.0_wp*q + 1.0_wp)
      else
         val = 0.0_wp
      end if
   end function wendland_c6_23d

   !> Derivative of Wendland C6 kernel for 1D: dW/dq, q < 2
   pure function wendland_c6_1d_deriv(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         ! Using the chain rule from Python reference
         val = -(63.0_wp/16.0_wp)*q*(q*q + q + 4.0_wp/9.0_wp)*tmp**6
      else
         val = 0.0_wp
      end if
   end function wendland_c6_1d_deriv

   !> Derivative of Wendland C6 kernel for 2D/3D: dW/dq, q < 2
   pure function wendland_c6_23d_deriv(q) result(val)
      real(wp), intent(in) :: q
      real(wp) :: val
      real(wp) :: tmp

      if (q < 2.0_wp) then
         tmp = 1.0_wp - 0.5_wp*q
         ! dW/dq = -22q(q^2 + q + 2/11)(1-q/2)^7
         val = -22.0_wp*q*(q*q + q + 2.0_wp/11.0_wp)*tmp**7
      else
         val = 0.0_wp
      end if
   end function wendland_c6_23d_deriv

   !> Initialize Wendland kernel with order and dimension-specific normalization
   !>
   !> @param[inout] self      Kernel instance
   !> @param[in]    order     Kernel order (2, 4 or 6)
   !> @param[in]    dimension Spatial dimension (1, 2 or 3)
   !> @param[in]    h         Smoothing length, must be positive
   !> @param[out]   error     Unsupported order/dimension, or non-positive `h`
   subroutine wendland_init(self, order, dimension, h, error)
      class(wendland_kernel_type), intent(inout) :: self
      integer, intent(in) :: order
      integer, intent(in) :: dimension
      real(wp), intent(in) :: h
      !> Unsupported order/dimension, or non-positive smoothing length
      type(error_type), allocatable, intent(out) :: error
      !> Reciprocal smoothing length shared by every normalization below
      real(wp) :: h_inv

      self%compute => null()
      self%compute_deriv => null()

      !> Guard the reciprocal below; a zero or negative smoothing length has no
      !> kernel and would otherwise divide by zero.
      if (h <= 0.0_wp) then
         call fatal_error(error, "wendland_init: smoothing length must be positive")
         return
      end if

      self%order = order
      self%dimension = dimension
      self%h = h
      h_inv = 1.0_wp/h

      ! Set order and dimension-specific normalization and kernel function
      select case (order)
      case (2)  ! Wendland C2
         select case (dimension)
         case (1)
            ! C2-1D: alpha = 5/(8h)
            self%prefactor = (5.0_wp/8.0_wp)*h_inv
            self%compute => wendland_c2_1d
            self%compute_deriv => wendland_c2_1d_deriv
         case (2)
            ! C2-2D: alpha = 7/(4 pi h^2 )
            self%prefactor = (7.0_wp/(4.0_wp*pi))*h_inv*h_inv
            self%compute => wendland_c2_23d
            self%compute_deriv => wendland_c2_23d_deriv
         case (3)
            ! C2-3D: alpha = 21/(16 pi h^3 )
            self%prefactor = (21.0_wp/(16.0_wp*pi))*h_inv*h_inv*h_inv
            self%compute => wendland_c2_23d
            self%compute_deriv => wendland_c2_23d_deriv
         case default
            call unsupported_dimension(error, "C2", dimension)
            return
         end select
      case (4)  ! Wendland C4
         select case (dimension)
         case (1)
            ! C4-1D: alpha = 3/(4h) = 0.75/h
            self%prefactor = 0.75_wp*h_inv
            self%compute => wendland_c4_1d
            self%compute_deriv => wendland_c4_1d_deriv
         case (2)
            ! C4-2D: alpha = 9/(4 pi h^2 )
            self%prefactor = (9.0_wp/(4.0_wp*pi))*h_inv*h_inv
            self%compute => wendland_c4_23d
            self%compute_deriv => wendland_c4_23d_deriv
         case (3)
            ! C4-3D: alpha = 495/(256 pi h^3 )
            self%prefactor = (495.0_wp/(256.0_wp*pi))*h_inv*h_inv*h_inv
            self%compute => wendland_c4_23d
            self%compute_deriv => wendland_c4_23d_deriv
         case default
            call unsupported_dimension(error, "C4", dimension)
            return
         end select
      case (6)  ! Wendland C6
         select case (dimension)
         case (1)
            ! C6-1D: alpha = 55/(64h)
            self%prefactor = (55.0_wp/64.0_wp)*h_inv
            self%compute => wendland_c6_1d
            self%compute_deriv => wendland_c6_1d_deriv
         case (2)
            ! C6-2D: alpha = 78/(28 pi h^2 ) = 39/(14 pi h^2 )
            self%prefactor = (78.0_wp/(28.0_wp*pi))*h_inv*h_inv
            self%compute => wendland_c6_23d
            self%compute_deriv => wendland_c6_23d_deriv
         case (3)
            ! C6-3D: alpha = 1365/(512 pi h^3 )
            self%prefactor = (1365.0_wp/(512.0_wp*pi))*h_inv*h_inv*h_inv
            self%compute => wendland_c6_23d
            self%compute_deriv => wendland_c6_23d_deriv
         case default
            call unsupported_dimension(error, "C6", dimension)
            return
         end select
      case default
         call fatal_error(error, "wendland_init: unsupported order "// &
                          trim(as_string(order))//" (must be 2, 4, or 6)")
         return
      end select
   end subroutine wendland_init

   !> Report a spatial dimension the requested Wendland order does not cover
   !>
   !> @param[out] error     Populated diagnostic
   !> @param[in]  kernel    Kernel label used in the message, e.g. "C2"
   !> @param[in]  dimension Offending spatial dimension
   subroutine unsupported_dimension(error, kernel, dimension)
      !> Populated diagnostic
      type(error_type), allocatable, intent(out) :: error
      !> Kernel label used in the message
      character(*), intent(in) :: kernel
      !> Offending spatial dimension
      integer, intent(in) :: dimension

      call fatal_error(error, "wendland_init: unsupported dimension "// &
                       trim(as_string(dimension))//" for "//kernel// &
                       " (must be 1, 2, or 3)")
   end subroutine unsupported_dimension

   !> Render an integer for use in a diagnostic message
   !>
   !> @param[in] value Integer to render
   !> @returns         Left-justified decimal representation
   pure function as_string(value) result(text)
      !> Integer to render
      integer, intent(in) :: value
      !> Left-justified decimal representation
      character(len=16) :: text

      write (text, "(i0)") value
   end function as_string

   !> Evaluate Wendland kernel at distance r
   pure function wendland_f0(self, r) result(kernel_val)
      class(wendland_kernel_type), intent(in) :: self
      real(wp), intent(in) :: r
      real(wp) :: kernel_val
      real(wp) :: q

      q = r/self%h
      kernel_val = self%prefactor*self%compute(q)
   end function wendland_f0

   !> Evaluate derivative dW/dr at distance r
   pure function wendland_f1(self, r) result(derivative)
      class(wendland_kernel_type), intent(in) :: self
      real(wp), intent(in) :: r
      real(wp) :: derivative
      real(wp) :: q

      q = r/self%h
      ! dW/dr = (prefactor * dW/dq) / h = prefactor * dW/dq * (1/h)
      derivative = self%prefactor*self%compute_deriv(q)/self%h
   end function wendland_f1

   !> Compute gradient vector dW = ( dW/ dx, dW/ dy, dW/ dz)
   pure subroutine wendland_gradient(self, r, xij, grad)
      class(wendland_kernel_type), intent(in) :: self
      real(wp), intent(in) :: r
      real(wp), intent(in) :: xij(3)
      real(wp), intent(out) :: grad(3)
      real(wp) :: tmp, dwdr

      ! Compute dW/dr and convert to gradient
      if (r > 1.0e-12_wp) then
         dwdr = self%f1(r)
         tmp = dwdr/r
      else
         tmp = 0.0_wp
      end if

      grad(1) = tmp*xij(1)
      grad(2) = tmp*xij(2)
      grad(3) = tmp*xij(3)
   end subroutine wendland_gradient

   !> Compute derivative with respect to smoothing length h: dW/dh
   pure function wendland_gradient_h(self, r) result(dwdh)
      class(wendland_kernel_type), intent(in) :: self
      real(wp), intent(in) :: r
      real(wp) :: dwdh
      real(wp) :: q, w, dwdq

      q = r/self%h

      ! Compute kernel value and derivative at q
      w = self%compute(q)
      dwdq = self%compute_deriv(q)

      ! dW/dh = -prefactor/h * (dwdq * q + w * dim)
      dwdh = -self%prefactor/self%h*(dwdq*q + w*real(self%dimension, wp))
   end function wendland_gradient_h

end module moist_math_smoothing_kernels
