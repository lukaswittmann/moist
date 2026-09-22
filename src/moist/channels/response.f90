!> MOIST -> host response items
!>
!> A response is a list of *items*, one per host contraction. An item that is
!> absent from the list is correct physics, not an error: a cavity with
!> field-independent geometry has no density response, a model without
!> GOSTSHYP has no Gaussian amplitudes. The host selects on the dynamic type of
!> `item(i)` and contracts what it finds
!>
!> Within one model call the items are accumulators over components:
!> `accumulate` finds the item of the same dynamic type and adds to it, or
!> appends a copy when there is none yet. Every `get_*` clears the list on
!> entry, so the sum never leaks across calls
module moist_channels_response
   use mctc_env, only: wp, error_type, fatal_error

   implicit none
   private

   public :: response_channel_type, response_slot, response_type
   public :: surface_charge_response_type, density_response_type
   public :: gostshyp_amplitude_response_type
   public :: find_surface_charge, find_density, find_gostshyp_amplitude
   public :: response_name_len

   !> Length of the fixed-size item names returned by `name()`
   integer, parameter :: response_name_len = 32

   !* ============================================================================== *!
   !*                              Response base type                                *!
   !* ============================================================================== *!

   !> One contraction handed back to the host, accumulated over components
   type, abstract :: response_channel_type
   contains
      !> Short fixed-length name for diagnostics, e.g. "surface_charge"
      procedure(response_item_name), deferred :: name
      !> Add another item of the same dynamic type into this one
      procedure(response_item_add), deferred :: add
      !> Deallocate every array of this item
      procedure(response_item_clear), deferred :: clear
   end type response_channel_type

   abstract interface

      !> Short fixed-length name of a response item for diagnostics
      function response_item_name(self) result(name)
         import :: response_channel_type, response_name_len
         !> Item to name
         class(response_channel_type), intent(in) :: self
         !> Name, blank padded
         character(len=response_name_len) :: name
      end function response_item_name

      !> Add another item of the same dynamic type into this one
      !>
      !> Arrays absent on `self` are copied from `other`; arrays absent on
      !> `other` contribute nothing; a shape mismatch is an error
      subroutine response_item_add(self, other, error)
         import :: response_channel_type, error_type
         !> Accumulator
         class(response_channel_type), intent(inout) :: self
         !> Item to add
         class(response_channel_type), intent(in) :: other
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine response_item_add

      !> Deallocate every array of a response item
      subroutine response_item_clear(self)
         import :: response_channel_type
         !> Item to clear
         class(response_channel_type), intent(inout) :: self
      end subroutine response_item_clear

   end interface

   !> Owning box that makes heterogeneous response items storable in one array
   type :: response_slot
      !> Concrete item owned by this slot
      class(response_channel_type), allocatable :: item
   end type response_slot

   !* ============================================================================== *!
   !*                              Concrete items                                    *!
   !* ============================================================================== *!

   !> Surface charge, the adjoint `dE/dphi_i` by stationarity
   !>
   !> The host contracts it with its own potential integrals,
   !> `F_uv += sum_i q_i V_uv(r_i)`. Only charge-like contributions may be
   !> accumulated here, so that the name stays true of the sum
   type, extends(response_channel_type) :: surface_charge_response_type
      !> Surface charge (ngrid)
      real(wp), allocatable :: q(:)
   contains
      procedure :: name => surface_charge_name
      procedure :: add => surface_charge_add
      procedure :: clear => surface_charge_clear
   end type surface_charge_response_type

   !> Weights conjugate to the solute density on the cavity grid
   !>
   !> The host contracts them with its own `d rho/dP`, `d grad rho/dP` and
   !> `d hess rho/dP`. The `dS/drho = -scale` factor of the level set is folded
   !> in by MOIST, so no convention crosses the boundary: the weights are
   !> conjugate to the density itself, never to the level set built from it
   !>
   !> The item is **present only for a cavity whose level set is a function of
   !> a density** (both isodensity variants). A geometric cavity -- iSwiG, or
   !> DROP on SvdW or CFC -- produces no density item at all, because its
   !> surface does not move with the host density; that absence is the physical
   !> answer and is distinguishable from an item of zeros, which zeros would
   !> not be
   type, extends(response_channel_type) :: density_response_type
      !> Weights for the density values (ngrid)
      real(wp), allocatable :: w_rho(:)
      !> Weights for the density gradients (3, ngrid)
      real(wp), allocatable :: w_grad_rho(:, :)
      !> Weights for the density Hessians (3, 3, ngrid)
      real(wp), allocatable :: w_hess_rho(:, :, :)
   contains
      procedure :: name => density_name
      procedure :: add => density_add
      procedure :: clear => density_clear
   end type density_response_type

   !> Amplitudes conjugate to the host's Gaussian integral blocks
   !>
   !> Counterpart of the Gaussian moment request. The host builds its Fock
   !> contribution as a plain sum over grid points,
   !>
   !>    F_uv += sum_i [ w_overlap(i) g_uv,i + w_normal_deriv(i) f_uv,i ]
   !>
   !> with `g_uv,i = <u|G_i|v>` and `f_uv,i = n_i . grad_r g_uv,i`. Both signs
   !> are folded in here. Grid points the model has switched off carry exactly
   !> zero, so the mask propagates without the host repeating it
   type, extends(response_channel_type) :: gostshyp_amplitude_response_type
      !> Overlap amplitudes (ngrid)
      real(wp), allocatable :: w_overlap(:)
      !> Normal derivative amplitudes (ngrid)
      real(wp), allocatable :: w_normal_deriv(:)
   contains
      procedure :: name => gostshyp_amplitude_name
      procedure :: add => gostshyp_amplitude_add
      procedure :: clear => gostshyp_amplitude_clear
   end type gostshyp_amplitude_response_type

   !* ============================================================================== *!
   !*                              Response container                                *!
   !* ============================================================================== *!

   !> Response item list handed back to the host for one model call
   !>
   !> `accumulate` is the only binding that grows `items(:)`; `clear` empties
   !> it. A variable of this type must be declared `target` wherever `item` or
   !> a finder is called on it
   type :: response_type
      !> Accumulated items, one per dynamic type
      type(response_slot), allocatable :: items(:)
   contains
      !> Number of items
      procedure :: n => response_n
      !> Pointer to the i-th item
      procedure :: item => response_item
      !> Add an item to the one of the same type, appending a copy if absent
      procedure :: accumulate => response_accumulate
      !> Drop every item
      procedure :: clear => response_clear
   end type response_type

contains

   !* ============================================================================== *!
   !*                              Array accumulation                                *!
   !* ============================================================================== *!

   !> Accumulate one rank-1 array: copy when absent, add when present
   !>
   !> @param[inout] acc   Accumulator, allocated on first use
   !> @param[in]    add   Array to add, absent contributes nothing
   !> @param[in]    item  Item name for diagnostics
   !> @param[in]    array Array name for diagnostics
   !> @param[out]   error Error handling
   subroutine accumulate_vector(acc, add, item, array, error)
      !> Accumulator
      real(wp), allocatable, intent(inout) :: acc(:)
      !> Array to add
      real(wp), allocatable, intent(in) :: add(:)
      !> Item name for diagnostics
      character(len=*), intent(in) :: item
      !> Array name for diagnostics
      character(len=*), intent(in) :: array
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(add)) return
      if (.not. allocated(acc)) then
         acc = add
         return
      end if
      if (any(shape(acc) /= shape(add))) then
         call shape_mismatch(item, array, error)
         return
      end if
      acc = acc + add

   end subroutine accumulate_vector

   !> Accumulate one rank-2 array: copy when absent, add when present
   !>
   !> @param[inout] acc   Accumulator, allocated on first use
   !> @param[in]    add   Array to add, absent contributes nothing
   !> @param[in]    item  Item name for diagnostics
   !> @param[in]    array Array name for diagnostics
   !> @param[out]   error Error handling
   subroutine accumulate_matrix(acc, add, item, array, error)
      !> Accumulator
      real(wp), allocatable, intent(inout) :: acc(:, :)
      !> Array to add
      real(wp), allocatable, intent(in) :: add(:, :)
      !> Item name for diagnostics
      character(len=*), intent(in) :: item
      !> Array name for diagnostics
      character(len=*), intent(in) :: array
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(add)) return
      if (.not. allocated(acc)) then
         acc = add
         return
      end if
      if (any(shape(acc) /= shape(add))) then
         call shape_mismatch(item, array, error)
         return
      end if
      acc = acc + add

   end subroutine accumulate_matrix

   !> Accumulate one rank-3 array: copy when absent, add when present
   !>
   !> @param[inout] acc   Accumulator, allocated on first use
   !> @param[in]    add   Array to add, absent contributes nothing
   !> @param[in]    item  Item name for diagnostics
   !> @param[in]    array Array name for diagnostics
   !> @param[out]   error Error handling
   subroutine accumulate_tensor3(acc, add, item, array, error)
      !> Accumulator
      real(wp), allocatable, intent(inout) :: acc(:, :, :)
      !> Array to add
      real(wp), allocatable, intent(in) :: add(:, :, :)
      !> Item name for diagnostics
      character(len=*), intent(in) :: item
      !> Array name for diagnostics
      character(len=*), intent(in) :: array
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(add)) return
      if (.not. allocated(acc)) then
         acc = add
         return
      end if
      if (any(shape(acc) /= shape(add))) then
         call shape_mismatch(item, array, error)
         return
      end if
      acc = acc + add

   end subroutine accumulate_tensor3

   !> Report an accumulation with mismatched shapes
   subroutine shape_mismatch(item, array, error)
      !> Item name
      character(len=*), intent(in) :: item
      !> Array name
      character(len=*), intent(in) :: array
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "Response item '"//item//"': accumulated '"//array// &
         & "' has a different shape than the stored one")

   end subroutine shape_mismatch

   !> Report an accumulation across different dynamic types
   subroutine type_mismatch(item, error)
      !> Item name
      character(len=*), intent(in) :: item
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "Response item '"//item// &
         & "' cannot accumulate an item of a different type")

   end subroutine type_mismatch

   !* ============================================================================== *!
   !*                              Surface charge item                               *!
   !* ============================================================================== *!

   !> Name of the surface charge item
   function surface_charge_name(self) result(name)
      !> Item
      class(surface_charge_response_type), intent(in) :: self
      !> Name
      character(len=response_name_len) :: name

      name = "surface_charge"

   end function surface_charge_name

   !> Add another surface charge item into this one
   subroutine surface_charge_add(self, other, error)
      !> Accumulator
      class(surface_charge_response_type), intent(inout) :: self
      !> Item to add
      class(response_channel_type), intent(in) :: other
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select type (other)
      type is (surface_charge_response_type)
         call accumulate_vector(self%q, other%q, "surface_charge", "q", error)
      class default
         call type_mismatch("surface_charge", error)
      end select

   end subroutine surface_charge_add

   !> Deallocate the surface charge
   subroutine surface_charge_clear(self)
      !> Item
      class(surface_charge_response_type), intent(inout) :: self

      if (allocated(self%q)) deallocate (self%q)

   end subroutine surface_charge_clear

   !* ============================================================================== *!
   !*                              Density item                                      *!
   !* ============================================================================== *!

   !> Name of the density item
   function density_name(self) result(name)
      !> Item
      class(density_response_type), intent(in) :: self
      !> Name
      character(len=response_name_len) :: name

      name = "density"

   end function density_name

   !> Add another density item into this one
   subroutine density_add(self, other, error)
      !> Accumulator
      class(density_response_type), intent(inout) :: self
      !> Item to add
      class(response_channel_type), intent(in) :: other
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select type (other)
      type is (density_response_type)
         call accumulate_vector(self%w_rho, other%w_rho, "density", "w_rho", error)
         if (allocated(error)) return
         call accumulate_matrix(self%w_grad_rho, other%w_grad_rho, "density", "w_grad_rho", error)
         if (allocated(error)) return
         call accumulate_tensor3(self%w_hess_rho, other%w_hess_rho, "density", "w_hess_rho", error)
      class default
         call type_mismatch("density", error)
      end select

   end subroutine density_add

   !> Deallocate the density weights
   subroutine density_clear(self)
      !> Item
      class(density_response_type), intent(inout) :: self

      if (allocated(self%w_rho)) deallocate (self%w_rho)
      if (allocated(self%w_grad_rho)) deallocate (self%w_grad_rho)
      if (allocated(self%w_hess_rho)) deallocate (self%w_hess_rho)

   end subroutine density_clear

   !* ============================================================================== *!
   !*                           GOSTSHYP amplitude item                              *!
   !* ============================================================================== *!

   !> Name of the GOSTSHYP amplitude item
   function gostshyp_amplitude_name(self) result(name)
      !> Item
      class(gostshyp_amplitude_response_type), intent(in) :: self
      !> Name
      character(len=response_name_len) :: name

      name = "gostshyp_amplitude"

   end function gostshyp_amplitude_name

   !> Add another GOSTSHYP amplitude item into this one
   subroutine gostshyp_amplitude_add(self, other, error)
      !> Accumulator
      class(gostshyp_amplitude_response_type), intent(inout) :: self
      !> Item to add
      class(response_channel_type), intent(in) :: other
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select type (other)
      type is (gostshyp_amplitude_response_type)
         call accumulate_vector(self%w_overlap, other%w_overlap, &
            & "gostshyp_amplitude", "w_overlap", error)
         if (allocated(error)) return
         call accumulate_vector(self%w_normal_deriv, other%w_normal_deriv, &
            & "gostshyp_amplitude", "w_normal_deriv", error)
      class default
         call type_mismatch("gostshyp_amplitude", error)
      end select

   end subroutine gostshyp_amplitude_add

   !> Deallocate the amplitudes
   subroutine gostshyp_amplitude_clear(self)
      !> Item
      class(gostshyp_amplitude_response_type), intent(inout) :: self

      if (allocated(self%w_overlap)) deallocate (self%w_overlap)
      if (allocated(self%w_normal_deriv)) deallocate (self%w_normal_deriv)

   end subroutine gostshyp_amplitude_clear

   !* ============================================================================== *!
   !*                              Container bindings                                *!
   !* ============================================================================== *!

   !> Number of items in the response
   function response_n(self) result(n)
      !> Response
      class(response_type), intent(in) :: self
      !> Number of items
      integer :: n

      n = 0
      if (allocated(self%items)) n = size(self%items)

   end function response_n

   !> Pointer to the i-th item, null when `i` is out of range
   !>
   !> @param[in] self Response, which must be a target
   !> @param[in] i    Item index
   function response_item(self, i) result(item)
      !> Response
      class(response_type), intent(in), target :: self
      !> Item index
      integer, intent(in) :: i
      !> Item, or null
      class(response_channel_type), pointer :: item

      item => null()
      if (.not. allocated(self%items)) return
      if (i < 1 .or. i > size(self%items)) return
      item => self%items(i)%item

   end function response_item

   !> Add an item into the stored one of the same dynamic type, or append a copy
   !>
   !> Appending moves the existing slots and copies only the new item (the
   !> gfortran safe pattern for polymorphic slot arrays)
   !>
   !> @param[inout] self  Response
   !> @param[in]    item  Item to accumulate
   !> @param[out]   error Error handling
   subroutine response_accumulate(self, item, error)
      !> Response
      class(response_type), intent(inout) :: self
      !> Item to accumulate
      class(response_channel_type), intent(in) :: item
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Grown slot array
      type(response_slot), allocatable :: grown(:)
      integer :: i, n, stat

      n = self%n()
      do i = 1, n
         if (same_type_as(self%items(i)%item, item)) then
            call self%items(i)%item%add(item, error)
            return
         end if
      end do

      allocate (grown(n + 1), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to grow response items")
         return
      end if
      allocate (grown(n + 1)%item, source=item, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate response item '"//trim(item%name())//"'")
         return
      end if
      ! Keep existing items intact until all allocations succeed
      do i = 1, n
         call move_alloc(self%items(i)%item, grown(i)%item)
      end do
      call move_alloc(grown, self%items)

   end subroutine response_accumulate

   !> Drop every item
   subroutine response_clear(self)
      !> Response
      class(response_type), intent(inout) :: self

      if (allocated(self%items)) deallocate (self%items)
      allocate (self%items(0))

   end subroutine response_clear

   !* ============================================================================== *!
   !*                              Per-type finders                                  *!
   !* ============================================================================== *!

   !> Find the surface charge item, null when absent
   !>
   !> @param[in] response Response, which must be a target
   function find_surface_charge(response) result(item)
      !> Response
      class(response_type), intent(in), target :: response
      !> Item, or null
      type(surface_charge_response_type), pointer :: item

      integer :: i

      item => null()
      if (.not. allocated(response%items)) return
      do i = 1, size(response%items)
         select type (stored => response%items(i)%item)
         type is (surface_charge_response_type)
            item => stored
            return
         end select
      end do

   end function find_surface_charge

   !> Find the density item, null when absent
   !>
   !> @param[in] response Response, which must be a target
   function find_density(response) result(item)
      !> Response
      class(response_type), intent(in), target :: response
      !> Item, or null
      type(density_response_type), pointer :: item

      integer :: i

      item => null()
      if (.not. allocated(response%items)) return
      do i = 1, size(response%items)
         select type (stored => response%items(i)%item)
         type is (density_response_type)
            item => stored
            return
         end select
      end do

   end function find_density

   !> Find the GOSTSHYP amplitude item, null when absent
   !>
   !> @param[in] response Response, which must be a target
   function find_gostshyp_amplitude(response) result(item)
      !> Response
      class(response_type), intent(in), target :: response
      !> Item, or null
      type(gostshyp_amplitude_response_type), pointer :: item

      integer :: i

      item => null()
      if (.not. allocated(response%items)) return
      do i = 1, size(response%items)
         select type (stored => response%items(i)%item)
         type is (gostshyp_amplitude_response_type)
            item => stored
            return
         end select
      end do

   end function find_gostshyp_amplitude

end module moist_channels_response
