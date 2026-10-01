!> Radial mappings from the reference interval [-1, 1] to radii
!>
!> A mapping turns reference nodes and `dx` weights into radii and `dr`
!> weights, applying abs(dr/dx) internally; parameters may depend on the
!> element z but are fixed under nuclear motion
!>
!> - Linear onto a finite interval [lower, upper]:
!>   r = (lower + upper)/2 + (upper - lower)/2 * x, dr/dx = (upper - lower)/2;
!>   the same affine arithmetic as the rule interval scaling; ignores z
!> - Becke rational onto [0, inf): r = (1 + x)/(1 - x)*p, dr/dx = 2*p/(1 - x)^2,
!>   with p a fixed scale or radius_factor*covalent_rad(z); diverges at x = 1,
!>   which is rejected
!> - HandyMod finite onto [rmin, rmax]: with d = rmax - rmin, c = 2^m,
!>   a = c*(1 - c + d), u = (1 + x)^m:
!>
!>       r     = rmin + d*u/(a - u*(d - c))
!>       dr/dx = d*m*(1 + x)^(m - 1)*a/(a - u*(d - c))^2
!>
!>   Requires d > 2^m - 1; ignores z. For m < 1 the derivative diverges at
!>   x = -1, which is rejected
!> - Knowles logarithmic onto [0, inf): with t = (1 + x)/2 and u = t^k:
!>
!>       r     = -R*log(1 - u)
!>       dr/dx = R*k*t^(k - 1)/(2*(1 - u))
!>
!>   R is a fixed scale or comes from a caller-supplied table indexed by atomic
!>   number. Midpoint nodes with k = 3 give the Mura-Knowles recipe. Diverges at
!>   x = 1, and for k < 1 the derivative diverges at x = -1; both are rejected
module moist_math_grid_radial_mapping
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, &
      & ieee_positive_inf
   use mctc_env, only: wp, error_type, fatal_error
   use moist_data_atomicrad, only: covalent_rad, max_elem
   implicit none(type, external)
   private

   public :: moist_math_grid_radial_mapping_type
   public :: check_mapping_input
   public :: check_mapping_output
   public :: moist_math_grid_radial_mapping_linear_type
   public :: new_linear_mapping
   public :: moist_math_grid_radial_mapping_becke_type
   public :: new_becke_mapping
   public :: moist_math_grid_radial_mapping_handymod_type
   public :: new_handymod_mapping
   public :: moist_math_grid_radial_mapping_knowles_type
   public :: new_knowles_mapping

   !> Abstract radial mapping
   type, abstract :: moist_math_grid_radial_mapping_type
   contains
      !> Forward map r = f(x) for element z
      procedure(mapping_map_i), deferred :: map
      !> Map reference nodes and dx weights to radii and dr weights for element z
      procedure(mapping_transform_i), deferred :: transform
   end type moist_math_grid_radial_mapping_type

   !> Deferred operations shared by every concrete mapping
   abstract interface
      !> Forward map r = f(x) for element z
      !>
      !> Outside the domain the result is not finite: +inf at a divergent
      !> endpoint, NaN for x outside [-1, 1] or an element without a scale
      !>
      !> @param[in] self  Mapping instance
      !> @param[in] z     Atomic number
      !> @param[in] x     Reference node on [-1, 1]
      pure function mapping_map_i(self, z, x) result(r)
         import :: moist_math_grid_radial_mapping_type, wp
         implicit none(type, external)
         !> Mapping instance
         class(moist_math_grid_radial_mapping_type), intent(in) :: self
         !> Atomic number
         integer, intent(in) :: z
         !> Reference node
         real(wp), intent(in) :: x
         !> Radius (bohr)
         real(wp) :: r
      end function mapping_map_i

      !> Map reference nodes and dx weights to radii and dr weights
      !>
      !> A node outside the domain, a divergent endpoint, or a non-finite
      !> generated radius or weight is an error; no node is dropped
      !>
      !> @param[in]  self   Mapping instance
      !> @param[in]  z      Atomic number
      !> @param[in]  x_ref  Reference nodes on [-1, 1], shape (n)
      !> @param[in]  w_ref  Reference dx weights, shape (n)
      !> @param[out] r      Radii in bohr, shape (n)
      !> @param[out] w      Radial dr weights in bohr, shape (n)
      !> @param[out] error  Set on invalid input or non-finite output
      subroutine mapping_transform_i(self, z, x_ref, w_ref, r, w, error)
         import :: moist_math_grid_radial_mapping_type, wp, error_type
         implicit none(type, external)
         !> Mapping instance
         class(moist_math_grid_radial_mapping_type), intent(in) :: self
         !> Atomic number
         integer, intent(in) :: z
         !> Reference nodes
         real(wp), intent(in) :: x_ref(:)
         !> Reference dx weights
         real(wp), intent(in) :: w_ref(:)
         !> Radii (bohr)
         real(wp), allocatable, intent(out) :: r(:)
         !> Radial dr weights (bohr)
         real(wp), allocatable, intent(out) :: w(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine mapping_transform_i
   end interface

   !> Linear mapping onto [lower, upper]
   type, extends(moist_math_grid_radial_mapping_type) :: moist_math_grid_radial_mapping_linear_type
      !> Lower radius in bohr
      real(wp) :: lower = 0.0_wp
      !> Upper radius in bohr
      real(wp) :: upper = 1.0_wp
   contains
      !> Forward map
      procedure :: map => linear_map
      !> Nodes and dr weights
      procedure :: transform => linear_transform
   end type moist_math_grid_radial_mapping_linear_type

   !> Becke rational mapping
   type, extends(moist_math_grid_radial_mapping_type) :: moist_math_grid_radial_mapping_becke_type
      !> Whether p = radius_factor*covalent_rad(z) instead of the fixed scale
      logical :: per_element = .false.
      !> Fixed radial scale p in bohr
      real(wp) :: scale = 1.0_wp
      !> Factor multiplying covalent_rad(z)
      real(wp) :: radius_factor = 0.0_wp
   contains
      !> Forward map
      procedure :: map => becke_map
      !> Nodes and dr weights
      procedure :: transform => becke_transform
   end type moist_math_grid_radial_mapping_becke_type

   !> HandyMod finite mapping
   type, extends(moist_math_grid_radial_mapping_type) :: moist_math_grid_radial_mapping_handymod_type
      !> Lower radius in bohr
      real(wp) :: rmin = 0.0_wp
      !> Upper radius in bohr
      real(wp) :: rmax = 2.0_wp
      !> Map exponent
      real(wp) :: m = 1.0_wp
   contains
      !> Forward map
      procedure :: map => handymod_map
      !> Nodes and dr weights
      procedure :: transform => handymod_transform
   end type moist_math_grid_radial_mapping_handymod_type

   !> Knowles logarithmic mapping
   type, extends(moist_math_grid_radial_mapping_type) :: moist_math_grid_radial_mapping_knowles_type
      !> Map exponent
      real(wp) :: k = 3.0_wp
      !> Whether R comes from the element table instead of the fixed scale
      logical :: per_element = .false.
      !> Fixed radial scale R in bohr
      real(wp) :: scale = 1.0_wp
      !> Radial scale R in bohr per atomic number, shape (max z)
      real(wp), allocatable :: element_scale(:)
   contains
      !> Forward map
      procedure :: map => knowles_map
      !> Nodes and dr weights
      procedure :: transform => knowles_transform
   end type moist_math_grid_radial_mapping_knowles_type

contains

   !> Validate reference nodes and weights and allocate the outputs
   !>
   !> Requires matching sizes, finite values, and nodes inside [-1, 1]
   !>
   !> @param[in]  x_ref  Reference nodes, shape (n)
   !> @param[in]  w_ref  Reference dx weights, shape (n)
   !> @param[in]  label  Error-message label
   !> @param[out] r      Radii, allocated to shape (n) and zeroed
   !> @param[out] w      Radial weights, allocated to shape (n) and zeroed
   !> @param[out] error  Set on invalid input or allocation failure
   subroutine check_mapping_input(x_ref, w_ref, label, r, w, error)
      !> Reference nodes
      real(wp), intent(in) :: x_ref(:)
      !> Reference dx weights
      real(wp), intent(in) :: w_ref(:)
      !> Error-message label
      character(len=*), intent(in) :: label
      !> Radii
      real(wp), allocatable, intent(out) :: r(:)
      !> Radial weights
      real(wp), allocatable, intent(out) :: w(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: stat

      if (size(x_ref) /= size(w_ref)) then
         call fatal_error(error, label//": reference nodes and weights differ in size")
         return
      end if
      if (.not. (all(ieee_is_finite(x_ref)) .and. all(ieee_is_finite(w_ref)))) then
         call fatal_error(error, label//": reference nodes and weights must be finite")
         return
      end if
      if (any(x_ref < -1.0_wp) .or. any(x_ref > 1.0_wp)) then
         call fatal_error(error, label//": reference node outside [-1, 1]")
         return
      end if
      allocate (r(size(x_ref)), stat=stat)
      if (stat == 0) allocate (w(size(x_ref)), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, label//": failed to allocate radial nodes")
         return
      end if
      r(:) = 0.0_wp
      w(:) = 0.0_wp
   end subroutine check_mapping_input

   !> Require finite generated radii and weights
   !>
   !> @param[in]  r      Radii, shape (n)
   !> @param[in]  w      Radial weights, shape (n)
   !> @param[in]  label  Error-message label
   !> @param[out] error  Set on a non-finite radius or weight
   subroutine check_mapping_output(r, w, label, error)
      !> Radii
      real(wp), intent(in) :: r(:)
      !> Radial weights
      real(wp), intent(in) :: w(:)
      !> Error-message label
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. (all(ieee_is_finite(r)) .and. all(ieee_is_finite(w)))) then
         call fatal_error(error, label//": generated non-finite radius or weight")
         return
      end if
   end subroutine check_mapping_output

   !> Create a linear mapping onto [lower, upper]
   !>
   !> @param[out] mapping  New mapping
   !> @param[in]  lower    Lower radius in bohr, finite, >= 0
   !> @param[in]  upper    Upper radius in bohr, finite, > lower
   !> @param[out] error    Set on invalid bounds
   subroutine new_linear_mapping(mapping, lower, upper, error)
      !> New mapping
      type(moist_math_grid_radial_mapping_linear_type), intent(out) :: mapping
      !> Lower radius in bohr
      real(wp), intent(in) :: lower
      !> Upper radius in bohr
      real(wp), intent(in) :: upper
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. (ieee_is_finite(lower) .and. ieee_is_finite(upper))) then
         call fatal_error(error, "Linear mapping: bounds must be finite")
         return
      end if
      if (lower < 0.0_wp) then
         call fatal_error(error, "Linear mapping: lower bound must be >= 0")
         return
      end if
      if (.not. lower < upper) then
         call fatal_error(error, "Linear mapping: lower bound must be below upper bound")
         return
      end if
      mapping%lower = lower
      mapping%upper = upper
   end subroutine new_linear_mapping

   !> Forward map r = (lower + upper)/2 + (upper - lower)/2 * x
   !>
   !> @param[in] self  Mapping instance
   !> @param[in] z     Atomic number (unused)
   !> @param[in] x     Reference node on [-1, 1]
   pure function linear_map(self, z, x) result(r)
      !> Mapping instance
      class(moist_math_grid_radial_mapping_linear_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Reference node
      real(wp), intent(in) :: x
      !> Radius (bohr); NaN outside [-1, 1]
      real(wp) :: r

      if (.not. ieee_is_finite(x)) then
         r = ieee_value(r, ieee_quiet_nan)
         return
      end if
      if (x < -1.0_wp .or. x > 1.0_wp) then
         r = ieee_value(r, ieee_quiet_nan)
         return
      end if
      r = 0.5_wp*(self%lower + self%upper) + 0.5_wp*(self%upper - self%lower)*x
   end function linear_map

   !> Map reference nodes and dx weights to radii and dr weights
   !>
   !> @param[in]  self   Mapping instance
   !> @param[in]  z      Atomic number (unused)
   !> @param[in]  x_ref  Reference nodes on [-1, 1], shape (n)
   !> @param[in]  w_ref  Reference dx weights, shape (n)
   !> @param[out] r      Radii in bohr, shape (n)
   !> @param[out] w      Radial dr weights in bohr, shape (n)
   !> @param[out] error  Set on invalid input or non-finite output
   subroutine linear_transform(self, z, x_ref, w_ref, r, w, error)
      !> Mapping instance
      class(moist_math_grid_radial_mapping_linear_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Reference nodes
      real(wp), intent(in) :: x_ref(:)
      !> Reference dx weights
      real(wp), intent(in) :: w_ref(:)
      !> Radii (bohr)
      real(wp), allocatable, intent(out) :: r(:)
      !> Radial dr weights (bohr)
      real(wp), allocatable, intent(out) :: w(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i
      real(wp) :: mid, half

      call check_mapping_input(x_ref, w_ref, "Linear mapping", r, w, error)
      if (allocated(error)) return

      mid = 0.5_wp*(self%lower + self%upper)
      half = 0.5_wp*(self%upper - self%lower)
      do i = 1, size(x_ref)
         r(i) = mid + half*x_ref(i)
         w(i) = w_ref(i)*abs(half)
      end do

      call check_mapping_output(r, w, "Linear mapping", error)
   end subroutine linear_transform

   !> Create a Becke mapping with a fixed scale or an element-dependent scale
   !>
   !> Exactly one of `scale` and `radius_factor` must be present
   !>
   !> @param[out] mapping        New mapping
   !> @param[out] error          Set on a missing, duplicated, or invalid scale
   !> @param[in]  scale          Fixed radial scale p in bohr, finite, > 0
   !> @param[in]  radius_factor  Factor on covalent_rad(z) in bohr, finite, > 0
   subroutine new_becke_mapping(mapping, error, scale, radius_factor)
      !> New mapping
      type(moist_math_grid_radial_mapping_becke_type), intent(out) :: mapping
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Fixed radial scale p in bohr
      real(wp), intent(in), optional :: scale
      !> Factor multiplying covalent_rad(z)
      real(wp), intent(in), optional :: radius_factor

      real(wp) :: value

      if (present(scale) .eqv. present(radius_factor)) then
         call fatal_error(error, "Becke mapping: give exactly one of scale and radius_factor")
         return
      end if
      if (present(scale)) then
         value = scale
      else
         value = radius_factor
      end if
      if (.not. ieee_is_finite(value)) then
         call fatal_error(error, "Becke mapping: scale must be finite")
         return
      end if
      if (value <= 0.0_wp) then
         call fatal_error(error, "Becke mapping: scale must be > 0")
         return
      end if
      mapping%per_element = present(radius_factor)
      if (mapping%per_element) then
         mapping%radius_factor = value
      else
         mapping%scale = value
      end if
   end subroutine new_becke_mapping

   !> Radial scale p for element z; NaN for an element without a covalent radius
   !>
   !> @param[in] self  Mapping instance
   !> @param[in] z     Atomic number
   pure function becke_scale(self, z) result(p)
      !> Mapping instance
      class(moist_math_grid_radial_mapping_becke_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Radial scale in bohr
      real(wp) :: p

      if (.not. self%per_element) then
         p = self%scale
      else if (z >= 1 .and. z <= max_elem) then
         p = self%radius_factor*covalent_rad(z)
      else
         p = ieee_value(p, ieee_quiet_nan)
      end if
   end function becke_scale

   !> Forward map r = (1 + x)/(1 - x)*p
   !>
   !> @param[in] self  Mapping instance
   !> @param[in] z     Atomic number
   !> @param[in] x     Reference node on [-1, 1)
   pure function becke_map(self, z, x) result(r)
      !> Mapping instance
      class(moist_math_grid_radial_mapping_becke_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Reference node
      real(wp), intent(in) :: x
      !> Radius (bohr); +inf at x = 1, NaN outside [-1, 1]
      real(wp) :: r

      if (.not. ieee_is_finite(x)) then
         r = ieee_value(r, ieee_quiet_nan)
         return
      end if
      if (x < -1.0_wp .or. x > 1.0_wp) then
         r = ieee_value(r, ieee_quiet_nan)
         return
      end if
      if (x >= 1.0_wp) then
         r = ieee_value(r, ieee_positive_inf)
         return
      end if
      r = (1.0_wp + x)/(1.0_wp - x)*becke_scale(self, z)
   end function becke_map

   !> Map reference nodes and dx weights to radii and dr weights
   !>
   !> @param[in]  self   Mapping instance
   !> @param[in]  z      Atomic number, 1..max_elem for an element-dependent scale
   !> @param[in]  x_ref  Reference nodes on [-1, 1), shape (n)
   !> @param[in]  w_ref  Reference dx weights, shape (n)
   !> @param[out] r      Radii in bohr, shape (n)
   !> @param[out] w      Radial dr weights in bohr, shape (n)
   !> @param[out] error  Set on invalid input, a node at x = 1, or non-finite output
   subroutine becke_transform(self, z, x_ref, w_ref, r, w, error)
      !> Mapping instance
      class(moist_math_grid_radial_mapping_becke_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Reference nodes
      real(wp), intent(in) :: x_ref(:)
      !> Reference dx weights
      real(wp), intent(in) :: w_ref(:)
      !> Radii (bohr)
      real(wp), allocatable, intent(out) :: r(:)
      !> Radial dr weights (bohr)
      real(wp), allocatable, intent(out) :: w(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i
      real(wp) :: p, one_minus_x

      call check_mapping_input(x_ref, w_ref, "Becke mapping", r, w, error)
      if (allocated(error)) return
      if (self%per_element .and. (z < 1 .or. z > max_elem)) then
         call fatal_error(error, "Becke mapping: atomic number has no covalent radius")
         return
      end if
      if (any(x_ref >= 1.0_wp)) then
         call fatal_error(error, "Becke mapping: reference node at the divergent endpoint x = 1")
         return
      end if

      p = becke_scale(self, z)
      do i = 1, size(x_ref)
         one_minus_x = 1.0_wp - x_ref(i)
         r(i) = (1.0_wp + x_ref(i))/one_minus_x*p
         w(i) = w_ref(i)*abs(2.0_wp*p/(one_minus_x*one_minus_x))
      end do

      call check_mapping_output(r, w, "Becke mapping", error)
   end subroutine becke_transform

   !> Create a HandyMod mapping onto [rmin, rmax]
   !>
   !> @param[out] mapping  New mapping
   !> @param[in]  rmin     Lower radius in bohr, finite, >= 0
   !> @param[in]  rmax     Upper radius in bohr, finite, rmax - rmin > 2^m - 1
   !> @param[in]  m        Map exponent, finite, > 0
   !> @param[out] error    Set on invalid parameters
   subroutine new_handymod_mapping(mapping, rmin, rmax, m, error)
      !> New mapping
      type(moist_math_grid_radial_mapping_handymod_type), intent(out) :: mapping
      !> Lower radius in bohr
      real(wp), intent(in) :: rmin
      !> Upper radius in bohr
      real(wp), intent(in) :: rmax
      !> Map exponent
      real(wp), intent(in) :: m
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: a, c, d

      if (.not. (ieee_is_finite(rmin) .and. ieee_is_finite(rmax))) then
         call fatal_error(error, "HandyMod mapping: radial bounds must be finite")
         return
      end if
      if (.not. ieee_is_finite(m)) then
         call fatal_error(error, "HandyMod mapping: m must be finite")
         return
      end if
      if (m <= 0.0_wp) then
         call fatal_error(error, "HandyMod mapping: m must be > 0")
         return
      end if
      ! 2**m must not overflow
      if (m >= real(maxexponent(1.0_wp), wp)) then
         call fatal_error(error, "HandyMod mapping: m is too large, 2^m overflows")
         return
      end if
      if (rmin < 0.0_wp) then
         call fatal_error(error, "HandyMod mapping: rmin must be >= 0")
         return
      end if
      if (rmax <= rmin) then
         call fatal_error(error, "HandyMod mapping: rmax must be greater than rmin")
         return
      end if
      d = rmax - rmin
      c = 2.0_wp**m
      if (d <= c - 1.0_wp) then
         call fatal_error(error, "HandyMod mapping: rmax - rmin must be greater than 2^m - 1")
         return
      end if
      a = c*(1.0_wp - c + d)
      if (a <= 0.0_wp .or. .not. ieee_is_finite(a)) then
         call fatal_error(error, "HandyMod mapping: invalid denominator scale")
         return
      end if
      mapping%rmin = rmin
      mapping%rmax = rmax
      mapping%m = m
   end subroutine new_handymod_mapping

   !> Forward map r = rmin + d*u/(a - u*(d - c))
   !>
   !> @param[in] self  Mapping instance
   !> @param[in] z     Atomic number (unused)
   !> @param[in] x     Reference node on [-1, 1]
   pure function handymod_map(self, z, x) result(r)
      !> Mapping instance
      class(moist_math_grid_radial_mapping_handymod_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Reference node
      real(wp), intent(in) :: x
      !> Radius (bohr); NaN outside [-1, 1]
      real(wp) :: r

      real(wp) :: a, c, d, u

      if (.not. ieee_is_finite(x)) then
         r = ieee_value(r, ieee_quiet_nan)
         return
      end if
      if (x < -1.0_wp .or. x > 1.0_wp) then
         r = ieee_value(r, ieee_quiet_nan)
         return
      end if
      d = self%rmax - self%rmin
      c = 2.0_wp**self%m
      a = c*(1.0_wp - c + d)
      u = (1.0_wp + x)**self%m
      r = self%rmin + d*u/(a - u*(d - c))
   end function handymod_map

   !> Map reference nodes and dx weights to radii and dr weights
   !>
   !> @param[in]  self   Mapping instance
   !> @param[in]  z      Atomic number (unused)
   !> @param[in]  x_ref  Reference nodes on [-1, 1], shape (n); x = -1 only for m >= 1
   !> @param[in]  w_ref  Reference dx weights, shape (n)
   !> @param[out] r      Radii in bohr, shape (n)
   !> @param[out] w      Radial dr weights in bohr, shape (n)
   !> @param[out] error  Set on invalid input, a divergent node, or non-finite output
   subroutine handymod_transform(self, z, x_ref, w_ref, r, w, error)
      !> Mapping instance
      class(moist_math_grid_radial_mapping_handymod_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Reference nodes
      real(wp), intent(in) :: x_ref(:)
      !> Reference dx weights
      real(wp), intent(in) :: w_ref(:)
      !> Radii (bohr)
      real(wp), allocatable, intent(out) :: r(:)
      !> Radial dr weights (bohr)
      real(wp), allocatable, intent(out) :: w(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i
      real(wp) :: a, c, d, den, drdx, u, x

      call check_mapping_input(x_ref, w_ref, "HandyMod mapping", r, w, error)
      if (allocated(error)) return
      if (self%m < 1.0_wp .and. any(x_ref <= -1.0_wp)) then
         call fatal_error(error, &
            & "HandyMod mapping: reference node at x = -1, where dr/dx diverges for m < 1")
         return
      end if

      d = self%rmax - self%rmin
      c = 2.0_wp**self%m
      a = c*(1.0_wp - c + d)
      do i = 1, size(x_ref)
         x = x_ref(i)
         u = (1.0_wp + x)**self%m
         den = a - u*(d - c)
         if (den <= 0.0_wp) then
            call fatal_error(error, "HandyMod mapping: invalid denominator")
            return
         end if
         r(i) = self%rmin + d*u/den
         drdx = d*self%m*(1.0_wp + x)**(self%m - 1.0_wp)*a/den**2
         w(i) = w_ref(i)*abs(drdx)
      end do

      call check_mapping_output(r, w, "HandyMod mapping", error)
   end subroutine handymod_transform

   !> Create a Knowles mapping with a fixed scale or an element table
   !>
   !> Exactly one of `scale` and `element_scale` must be present
   !>
   !> @param[out] mapping        New mapping
   !> @param[in]  k              Map exponent, finite, > 0
   !> @param[out] error          Set on invalid parameters
   !> @param[in]  scale          Fixed radial scale R in bohr, finite, > 0
   !> @param[in]  element_scale  R in bohr indexed by atomic number, shape (max z), finite, > 0
   subroutine new_knowles_mapping(mapping, k, error, scale, element_scale)
      !> New mapping
      type(moist_math_grid_radial_mapping_knowles_type), intent(out) :: mapping
      !> Map exponent
      real(wp), intent(in) :: k
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Fixed radial scale R in bohr
      real(wp), intent(in), optional :: scale
      !> Radial scale R in bohr indexed by atomic number
      real(wp), intent(in), optional :: element_scale(:)

      if (.not. ieee_is_finite(k)) then
         call fatal_error(error, "Knowles mapping: k must be finite")
         return
      end if
      if (k <= 0.0_wp) then
         call fatal_error(error, "Knowles mapping: k must be > 0")
         return
      end if
      if (present(scale) .eqv. present(element_scale)) then
         call fatal_error(error, "Knowles mapping: give exactly one of scale and element_scale")
         return
      end if
      if (present(scale)) then
         if (.not. ieee_is_finite(scale)) then
            call fatal_error(error, "Knowles mapping: scale must be finite")
            return
         end if
         if (scale <= 0.0_wp) then
            call fatal_error(error, "Knowles mapping: scale must be > 0")
            return
         end if
         mapping%scale = scale
      else
         if (size(element_scale) < 1) then
            call fatal_error(error, "Knowles mapping: element table is empty")
            return
         end if
         if (.not. all(ieee_is_finite(element_scale))) then
            call fatal_error(error, "Knowles mapping: element scales must be finite")
            return
         end if
         if (any(element_scale <= 0.0_wp)) then
            call fatal_error(error, "Knowles mapping: element scales must be > 0")
            return
         end if
         mapping%per_element = .true.
         mapping%element_scale = element_scale
      end if
      mapping%k = k
   end subroutine new_knowles_mapping

   !> Radial scale R for element z; NaN for an element outside the table
   !>
   !> @param[in] self  Mapping instance
   !> @param[in] z     Atomic number
   pure function knowles_scale(self, z) result(scale)
      !> Mapping instance
      class(moist_math_grid_radial_mapping_knowles_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Radial scale in bohr
      real(wp) :: scale

      if (.not. self%per_element) then
         scale = self%scale
      else if (z >= 1 .and. z <= size(self%element_scale)) then
         scale = self%element_scale(z)
      else
         scale = ieee_value(scale, ieee_quiet_nan)
      end if
   end function knowles_scale

   !> -log(1 - u) for 0 <= u < 1, accurate for small u
   !>
   !> Uses log1p(y) = y*log(1 + y)/((1 + y) - 1) with the rounded 1 + y
   !>
   !> @param[in] u            Argument, 0 <= u < 1
   !> @param[in] one_minus_u  Rounded 1 - u, > 0
   pure function neg_log1m(u, one_minus_u) result(val)
      !> Argument
      real(wp), intent(in) :: u
      !> Rounded 1 - u
      real(wp), intent(in) :: one_minus_u
      !> -log(1 - u)
      real(wp) :: val

      if (one_minus_u >= 1.0_wp) then
         val = u
      else
         val = u*log(one_minus_u)/(one_minus_u - 1.0_wp)
      end if
   end function neg_log1m

   !> Forward map r = -R*log(1 - ((1 + x)/2)^k)
   !>
   !> @param[in] self  Mapping instance
   !> @param[in] z     Atomic number
   !> @param[in] x     Reference node on [-1, 1)
   pure function knowles_map(self, z, x) result(r)
      !> Mapping instance
      class(moist_math_grid_radial_mapping_knowles_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Reference node
      real(wp), intent(in) :: x
      !> Radius (bohr); +inf at x = 1, NaN outside [-1, 1]
      real(wp) :: r

      real(wp) :: u, one_minus_u

      if (.not. ieee_is_finite(x)) then
         r = ieee_value(r, ieee_quiet_nan)
         return
      end if
      if (x < -1.0_wp .or. x > 1.0_wp) then
         r = ieee_value(r, ieee_quiet_nan)
         return
      end if
      u = (0.5_wp*(1.0_wp + x))**self%k
      one_minus_u = 1.0_wp - u
      if (one_minus_u <= 0.0_wp) then
         r = ieee_value(r, ieee_positive_inf)
         return
      end if
      r = knowles_scale(self, z)*neg_log1m(u, one_minus_u)
   end function knowles_map

   !> Map reference nodes and dx weights to radii and dr weights
   !>
   !> @param[in]  self   Mapping instance
   !> @param[in]  z      Atomic number, inside the element table if one is set
   !> @param[in]  x_ref  Reference nodes on [-1, 1), shape (n); x = -1 only for k >= 1
   !> @param[in]  w_ref  Reference dx weights, shape (n)
   !> @param[out] r      Radii in bohr, shape (n)
   !> @param[out] w      Radial dr weights in bohr, shape (n)
   !> @param[out] error  Set on invalid input, a divergent node, or non-finite output
   subroutine knowles_transform(self, z, x_ref, w_ref, r, w, error)
      !> Mapping instance
      class(moist_math_grid_radial_mapping_knowles_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Reference nodes
      real(wp), intent(in) :: x_ref(:)
      !> Reference dx weights
      real(wp), intent(in) :: w_ref(:)
      !> Radii (bohr)
      real(wp), allocatable, intent(out) :: r(:)
      !> Radial dr weights (bohr)
      real(wp), allocatable, intent(out) :: w(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i
      real(wp) :: scale, t, u, one_minus_u, drdx

      call check_mapping_input(x_ref, w_ref, "Knowles mapping", r, w, error)
      if (allocated(error)) return
      if (self%per_element) then
         if (z < 1 .or. z > size(self%element_scale)) then
            call fatal_error(error, "Knowles mapping: atomic number outside the element table")
            return
         end if
      end if
      if (self%k < 1.0_wp .and. any(x_ref <= -1.0_wp)) then
         call fatal_error(error, &
            & "Knowles mapping: reference node at x = -1, where dr/dx diverges for k < 1")
         return
      end if

      scale = knowles_scale(self, z)
      do i = 1, size(x_ref)
         t = 0.5_wp*(1.0_wp + x_ref(i))
         u = t**self%k
         one_minus_u = 1.0_wp - u
         if (one_minus_u <= 0.0_wp) then
            call fatal_error(error, &
               & "Knowles mapping: reference node at the divergent endpoint x = 1")
            return
         end if
         r(i) = scale*neg_log1m(u, one_minus_u)
         drdx = scale*self%k*t**(self%k - 1.0_wp)/(2.0_wp*one_minus_u)
         w(i) = w_ref(i)*abs(drdx)
      end do

      call check_mapping_output(r, w, "Knowles mapping", error)
   end subroutine knowles_transform

end module moist_math_grid_radial_mapping
