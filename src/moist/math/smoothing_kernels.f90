!> Implementation of various smoothing Kernels
!> References:
!> https://pysph.readthedocs.io/en/main/reference/kernels.html
!> https://ludwigboess.github.io/SPHKernels.jl/stable/kernels/
module moist_math_smoothing_kernels

   use mctc_env, only: wp
   use mctc_io_constants, only: pi

   implicit none
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
         real(wp), intent(in) :: q
         real(wp) :: val
      end function wendland_compute_interface
   end interface

   abstract interface
      !> Initialize kernel interface
      subroutine init_interface(self, order, dimension, h)
         import :: smoothing_kernel_type, wp
         class(smoothing_kernel_type), intent(inout) :: self
         integer, intent(in) :: order
         integer, intent(in) :: dimension
         real(wp), intent(in) :: h
      end subroutine init_interface

      !> Kernel evaluation interface
      pure function f0_interface(self, r) result(kernel_val)
         import :: smoothing_kernel_type, wp
         class(smoothing_kernel_type), intent(in) :: self
         real(wp), intent(in) :: r
         real(wp) :: kernel_val
      end function f0_interface

      !> Kernel derivative interface (dW/dr)
      pure function f1_interface(self, r) result(derivative)
         import :: smoothing_kernel_type, wp
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
   subroutine wendland_init(self, order, dimension, h)
      class(wendland_kernel_type), intent(inout) :: self
      integer, intent(in) :: order
      integer, intent(in) :: dimension
      real(wp), intent(in) :: h
      real(wp) :: h_inv

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
            ! TODO: Proper error propagration
            error stop "wendland_init: unsupported dimension for C2 (must be 1, 2, or 3)"
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
            ! TODO: Proper error propagration
            error stop "wendland_init: unsupported dimension for C4 (must be 1, 2, or 3)"
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
            ! TODO: Proper error propagration
            error stop "wendland_init: unsupported dimension for C6 (must be 1, 2, or 3)"
         end select
      case default
         ! TODO: Proper error propagration
         error stop "wendland_init: unsupported order (must be 2, 4, or 6)"
      end select
   end subroutine wendland_init

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
