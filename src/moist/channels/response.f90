!> MOIST -> host response items
!>
!> A response is a list of items, one per host contraction
!>
!> - an item absent from the list is correct physics, not an error: a cavity
!>   with field-independent geometry has no density response, a model without
!>   GOSTSHYP no Gaussian amplitudes
!> - the host walks the items with a cursor, `do while (response%next())`,
!>   selects on the dynamic type of the copy `response%item()` and contracts
!>   what it finds
!> - a default branch that stops the host is the one check that no item goes
!>   uncontracted: unlike an unanswered request, a skipped item fails nowhere
!> - every item also declares its arrays through `list_fields`, so they can be
!>   listed and read by name without selecting on the type
!>
!> Within one model call the items are accumulators over components
!>
!> - `response_accumulate` finds the item of the same dynamic type and adds to
!>   it, or appends a copy when there is none yet
!> - a `get_*` clears the list once it accepts the request, so the sum never
!>   leaks across calls and every walk starts afresh
!> - a rejected request, including an unimplemented theory, leaves the list
!>   untouched
module moist_channels_response
   use mctc_env, only: wp, error_type, fatal_error
   use moist_channels_fields, only: field_query_type

   implicit none(type, external)
   private

   public :: response_item_type, response_type
   public :: potential_adjoint_response_type, density_response_type
   public :: gaussian_amplitude_response_type, atomic_multipole_adjoint_response_type
   public :: atomic_charge_adjoint_response_type, radial_potential_adjoint_response_type
   public :: current_response_item, response_accumulate, response_clear
   public :: response_name_len

   !> Length of the fixed-size item names returned by `name()`
   integer, parameter :: response_name_len = 32

   !> Error of a read outside a `next()` window
   character(len=*), parameter :: no_current_item = &
      & "No current response item - call next() first"

   !* ============================================================================== *!
   !*                              Response base type                                *!
   !* ============================================================================== *!

   !> One contraction handed back to the host, accumulated over components
   type, abstract :: response_item_type
   contains
      !> Short fixed-length name for diagnostics, e.g. "potential_adjoint"
      procedure(response_item_name), deferred :: name
      !> Declare the item's arrays, named like its components
      procedure(response_item_list_fields), deferred :: list_fields
      procedure(response_item_add), deferred, private :: add
      procedure(response_item_clear), deferred, private :: clear
   end type response_item_type

   abstract interface

      !> Short fixed-length name of a response item for diagnostics
      function response_item_name(self) result(name)
         import :: response_item_type, response_name_len
         implicit none(type, external)
         !> Item to name
         class(response_item_type), intent(in) :: self
         !> Name, blank padded
         character(len=response_name_len) :: name
      end function response_item_name

      !> Declare every allocated array of a response item
      !>
      !> Field names are the component names; an array the item was
      !> accumulated without is not declared
      subroutine response_item_list_fields(self, query)
         import :: response_item_type, field_query_type
         implicit none(type, external)
         !> Item to describe
         class(response_item_type), intent(in) :: self
         !> Field walker
         type(field_query_type), intent(inout) :: query
      end subroutine response_item_list_fields

      !> Add another item of the same dynamic type into this one
      !>
      !> Arrays absent on `self` are copied from `other`; arrays absent on
      !> `other` contribute nothing; a shape mismatch is an error
      subroutine response_item_add(self, other, error)
         import :: response_item_type, error_type
         implicit none(type, external)
         !> Accumulator
         class(response_item_type), intent(inout) :: self
         !> Item to add
         class(response_item_type), intent(in) :: other
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine response_item_add

      !> Deallocate every array of a response item
      subroutine response_item_clear(self)
         import :: response_item_type
         implicit none(type, external)
         !> Item to clear
         class(response_item_type), intent(inout) :: self
      end subroutine response_item_clear

   end interface

   !> Owning box that makes heterogeneous response items storable in one array
   type :: response_slot
      !> Concrete item owned by this slot
      class(response_item_type), allocatable :: item
   end type response_slot

   !* ============================================================================== *!
   !*                              Concrete items                                    *!
   !* ============================================================================== *!

   !> Weights conjugate to the host potential on the cavity grid, `dE/dphi_i`
   !>
   !> The host contracts them with its own potential integrals,
   !> `F_uv += sum_i w_phi(i) V_uv(r_i)`, and with their basis-center
   !> derivative in the gradient phase
   !>
   !> - for a stationary PCM (CPCM, COSMO) the weights are the induced surface
   !>   charges `q_i`; for a non-symmetric response matrix (IEF-PCM, SS(V)PE)
   !>   they are the symmetrized adjoint, which is not the apparent charge
   !> - any component whose energy depends on the potential accumulates here,
   !>   so the sum is the adjoint of the model energy, never a charge of one
   !>   component
   type, extends(response_item_type) :: potential_adjoint_response_type
      !> Weights for the potential values (ngrid)
      real(wp), allocatable :: w_phi(:)
   contains
      procedure :: name => potential_adjoint_name
      procedure :: list_fields => potential_adjoint_list_fields
      procedure, private :: add => potential_adjoint_add
      procedure, private :: clear => potential_adjoint_clear
   end type potential_adjoint_response_type

   !> Weights conjugate to the solute density on the cavity grid
   !>
   !> The host contracts them with its own `d rho/dP`, `d grad rho/dP` and
   !> `d hess rho/dP`. The `dS/drho = -scale` factor of the level set is folded
   !> in by MOIST, so no convention crosses the boundary: the weights are
   !> conjugate to the density itself, never to the level set built from it
   !>
   !> The item is **present only for a cavity whose level set is a function of
   !> a density** (both isodensity variants); a geometric cavity -- iSwiG, or
   !> DROP on SvdW or CFC -- produces no density item at all, because its
   !> surface does not move with the host density; that absence is the physical
   !> answer and is distinguishable from an item of zeros, which zeros would
   !> not be
   !>
   !> The **response phase alone** produces it. The weights are `dE/drho` at
   !> fixed nuclei, the same object in either phase, so they are formed once and
   !> the gradient phase leaves the contraction out -- a host that does not need
   !> them would otherwise pay for it. A host on a density-backed cavity runs
   !> the response phase first and carries the weights over: against its own
   !> `d rho/dP` they complete the Fock matrix, against the basis-center
   !> `d rho/dR` the nuclear gradient
   type, extends(response_item_type) :: density_response_type
      !> Weights for the density values (ngrid)
      real(wp), allocatable :: w_rho(:)
      !> Weights for the density gradients (3, ngrid)
      real(wp), allocatable :: w_grad_rho(:, :)
      !> Weights for the density Hessians (3, 3, ngrid)
      real(wp), allocatable :: w_hess_rho(:, :, :)
   contains
      procedure :: name => density_name
      procedure :: list_fields => density_list_fields
      procedure, private :: add => density_add
      procedure, private :: clear => density_clear
   end type density_response_type

   !> Amplitudes conjugate to the host's Gaussian integral blocks
   !>
   !> Counterpart of the Gaussian moment request; the host builds its Fock
   !> contribution as a plain sum over grid points,
   !>
   !>    F_uv += sum_i [ w_overlap(i) g_uv,i + w_normal_deriv(i) f_uv,i ]
   !>
   !> with `g_uv,i = <u|G_i|v>` and `f_uv,i = n_i . grad_r g_uv,i`
   !>
   !> - both signs are folded in here
   !> - grid points the model has switched off carry exactly zero, so the mask
   !>   propagates without the host repeating it
   type, extends(response_item_type) :: gaussian_amplitude_response_type
      !> Overlap amplitudes (ngrid)
      real(wp), allocatable :: w_overlap(:)
      !> Normal derivative amplitudes (ngrid)
      real(wp), allocatable :: w_normal_deriv(:)
   contains
      procedure :: name => gaussian_amplitude_name
      procedure :: list_fields => gaussian_amplitude_list_fields
      procedure, private :: add => gaussian_amplitude_add
      procedure, private :: clear => gaussian_amplitude_clear
   end type gaussian_amplitude_response_type

   !> Weights conjugate to the solute's per-atom multipoles, `dE/dq_A`,...
   !>
   !> Counterpart of the `atomic_multipoles` request; host contracts them
   !> with its own derivatives of the charge model, `d q_A/dP` for the Fock
   !> matrix and `d q_A/dR` at fixed density for the nuclear gradient
   !>
   !> Every array runs over the solute atoms; an array of a multipole
   !> order the model never consumed is absent
   type, extends(response_item_type) :: atomic_multipole_adjoint_response_type
      !> Weights for the atomic charges (natom)
      real(wp), allocatable :: dg_dq(:)
      !> Weights for the atomic dipoles (3, natom)
      real(wp), allocatable :: dg_dmu(:, :)
      !> Weights for the atomic quadrupoles (3, 3, natom)
      real(wp), allocatable :: dg_dtheta(:, :, :)
   contains
      procedure :: name => atomic_multipole_adjoint_name
      procedure :: list_fields => atomic_multipole_adjoint_list_fields
      procedure, private :: add => atomic_multipole_adjoint_add
      procedure, private :: clear => atomic_multipole_adjoint_clear
   end type atomic_multipole_adjoint_response_type

   !> Weights conjugate to the solute's per-atom partial charges, `dE/dq_A`
   !>
   !> Counterpart of the `atomic_charges` request; the host contracts them
   !> with its own `d q_A/dR` and its charge operator, `d q_A/dP` for the
   !> Fock matrix
   type, extends(response_item_type) :: atomic_charge_adjoint_response_type
      !> Weights for the atomic charges (natom)
      real(wp), allocatable :: dg_dq(:)
   contains
      procedure :: name => atomic_charge_adjoint_name
      procedure :: list_fields => atomic_charge_adjoint_list_fields
      procedure, private :: add => atomic_charge_adjoint_add
      procedure, private :: clear => atomic_charge_adjoint_clear
   end type atomic_charge_adjoint_response_type

   !> Weights conjugate to the site-resolved radial potential, `dE/dphi_a(r_i)`
   !>
   !> Counterpart of the `radial_potential` request; the host contracts them
   !> with its own site-decomposed potential operator for the Fock matrix
   !> and with its own d phi / dR for the nuclear gradient
   type, extends(response_item_type) :: radial_potential_adjoint_response_type
      !> Weights for the radial potential values (ngrid, natom)
      real(wp), allocatable :: dg_dphi(:, :)
   contains
      procedure :: name => radial_potential_adjoint_name
      procedure :: list_fields => radial_potential_adjoint_list_fields
      procedure, private :: add => radial_potential_adjoint_add
      procedure, private :: clear => radial_potential_adjoint_clear
   end type radial_potential_adjoint_response_type

   !> Placeholder `response%item()` returns outside a `next()` window; no arrays
   type, extends(response_item_type) :: no_response_item_type
   contains
      procedure :: name => no_item_name
      procedure :: list_fields => no_item_list_fields
      procedure, private :: add => no_item_add
      procedure, private :: clear => no_item_clear
   end type no_response_item_type

   !* ============================================================================== *!
   !*                              Response container                                *!
   !* ============================================================================== *!

   !> Response item list handed back to the host for one model call
   !>
   !> Hosts use `next` and `item`; the model and its components grow and empty
   !> the list through the module procedures `response_accumulate` and
   !> `response_clear`, which the `moist` umbrella does not re-export
   !>
   !> - `response_accumulate` is the only route that grows `items(:)`,
   !>   `response_clear` empties it; both reset the cursor
   type :: response_type
      private
      !> Accumulated items, one per dynamic type
      type(response_slot), allocatable :: items(:)
      !> Item of the host walk; zero before the first and after the last
      integer :: cursor = 0
   contains
      !> Advance to the next item
      procedure :: next => response_next
      !> Copy of the current item
      procedure :: item => response_item
      !> Number of items
      procedure, private :: n => response_n
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
   !*                            Potential adjoint item                             *!
   !* ============================================================================== *!

   !> Name of the potential adjoint item
   function potential_adjoint_name(self) result(name)
      !> Item
      class(potential_adjoint_response_type), intent(in) :: self
      !> Name
      character(len=response_name_len) :: name

      name = "potential_adjoint"

   end function potential_adjoint_name

   !> Declare the potential weights
   !>
   !> @param[in]    self   Item
   !> @param[inout] query  Field walker
   subroutine potential_adjoint_list_fields(self, query)
      !> Item
      class(potential_adjoint_response_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query

      call query%add_real("w_phi", "dE/dphi on the cavity grid (ngrid)", self%w_phi)

   end subroutine potential_adjoint_list_fields

   !> Add another potential adjoint item into this one
   subroutine potential_adjoint_add(self, other, error)
      !> Accumulator
      class(potential_adjoint_response_type), intent(inout) :: self
      !> Item to add
      class(response_item_type), intent(in) :: other
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select type (other)
      type is (potential_adjoint_response_type)
         call accumulate_vector(self%w_phi, other%w_phi, "potential_adjoint", "w_phi", error)
      class default
         call type_mismatch("potential_adjoint", error)
      end select

   end subroutine potential_adjoint_add

   !> Deallocate the potential weights
   subroutine potential_adjoint_clear(self)
      !> Item
      class(potential_adjoint_response_type), intent(inout) :: self

      if (allocated(self%w_phi)) deallocate (self%w_phi)

   end subroutine potential_adjoint_clear

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

   !> Declare the density weights
   !>
   !> @param[in]    self   Item
   !> @param[inout] query  Field walker
   subroutine density_list_fields(self, query)
      !> Item
      class(density_response_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query

      call query%add_real("w_rho", "dE/drho on the cavity grid (ngrid)", self%w_rho)
      call query%add_real2("w_grad_rho", "dE/d(grad rho) on the cavity grid (3, ngrid)", self%w_grad_rho)
      call query%add_real3("w_hess_rho", "dE/d(hess rho) on the cavity grid (3, 3, ngrid)", self%w_hess_rho)

   end subroutine density_list_fields

   !> Add another density item into this one
   subroutine density_add(self, other, error)
      !> Accumulator
      class(density_response_type), intent(inout) :: self
      !> Item to add
      class(response_item_type), intent(in) :: other
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
   !*                           Gaussian amplitude item                              *!
   !* ============================================================================== *!

   !> Name of the Gaussian amplitude item
   function gaussian_amplitude_name(self) result(name)
      !> Item
      class(gaussian_amplitude_response_type), intent(in) :: self
      !> Name
      character(len=response_name_len) :: name

      name = "gaussian_amplitude"

   end function gaussian_amplitude_name

   !> Declare the amplitudes
   !>
   !> @param[in]    self   Item
   !> @param[inout] query  Field walker
   subroutine gaussian_amplitude_list_fields(self, query)
      !> Item
      class(gaussian_amplitude_response_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query

      call query%add_real("w_overlap", "dE/dg_i, conjugate to the overlap blocks <u|G_i|v> (ngrid)", &
         & self%w_overlap)
      call query%add_real("w_normal_deriv", "dE/df_i, conjugate to the blocks n_i . grad <u|G_i|v> (ngrid)", &
         & self%w_normal_deriv)

   end subroutine gaussian_amplitude_list_fields

   !> Add another Gaussian amplitude item into this one
   subroutine gaussian_amplitude_add(self, other, error)
      !> Accumulator
      class(gaussian_amplitude_response_type), intent(inout) :: self
      !> Item to add
      class(response_item_type), intent(in) :: other
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select type (other)
      type is (gaussian_amplitude_response_type)
         call accumulate_vector(self%w_overlap, other%w_overlap, &
            & "gaussian_amplitude", "w_overlap", error)
         if (allocated(error)) return
         call accumulate_vector(self%w_normal_deriv, other%w_normal_deriv, &
            & "gaussian_amplitude", "w_normal_deriv", error)
      class default
         call type_mismatch("gaussian_amplitude", error)
      end select

   end subroutine gaussian_amplitude_add

   !> Deallocate the amplitudes
   subroutine gaussian_amplitude_clear(self)
      !> Item
      class(gaussian_amplitude_response_type), intent(inout) :: self

      if (allocated(self%w_overlap)) deallocate (self%w_overlap)
      if (allocated(self%w_normal_deriv)) deallocate (self%w_normal_deriv)

   end subroutine gaussian_amplitude_clear

   !* ============================================================================== *!
   !*                         Atomic multipole adjoint item                          *!
   !* ============================================================================== *!

   !> Name of the atomic multipole adjoint item
   !>
   !> @param[in] self Item
   function atomic_multipole_adjoint_name(self) result(name)
      !> Item
      class(atomic_multipole_adjoint_response_type), intent(in) :: self
      !> Name
      character(len=response_name_len) :: name

      name = "atomic_multipole_adjoint"

   end function atomic_multipole_adjoint_name

   !> Declare the multipole weights the model consumed
   !>
   !> @param[in]    self   Item
   !> @param[inout] query  Field walker
   subroutine atomic_multipole_adjoint_list_fields(self, query)
      !> Item
      class(atomic_multipole_adjoint_response_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query

      call query%add_real("dg_dq", "dE/dq of the atomic charges (natom)", self%dg_dq)
      call query%add_real2("dg_dmu", "dE/dmu of the atomic dipoles (3, natom)", self%dg_dmu)
      call query%add_real3("dg_dtheta", "dE/dtheta of the atomic quadrupoles (3, 3, natom)", self%dg_dtheta)

   end subroutine atomic_multipole_adjoint_list_fields

   !> Add another atomic multipole adjoint item into this one
   !>
   !> @param[in,out] self Accumulator
   !> @param[in] other Item to add
   !> @param[out] error Error handling
   subroutine atomic_multipole_adjoint_add(self, other, error)
      !> Accumulator
      class(atomic_multipole_adjoint_response_type), intent(inout) :: self
      !> Item to add
      class(response_item_type), intent(in) :: other
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select type (other)
      type is (atomic_multipole_adjoint_response_type)
         call accumulate_vector(self%dg_dq, other%dg_dq, "atomic_multipole_adjoint", "dg_dq", error)
         if (allocated(error)) return
         call accumulate_matrix(self%dg_dmu, other%dg_dmu, "atomic_multipole_adjoint", "dg_dmu", error)
         if (allocated(error)) return
         call accumulate_tensor3(self%dg_dtheta, other%dg_dtheta, "atomic_multipole_adjoint", &
            & "dg_dtheta", error)
      class default
         call type_mismatch("atomic_multipole_adjoint", error)
      end select

   end subroutine atomic_multipole_adjoint_add

   !> Deallocate the multipole weights
   !>
   !> @param[in,out] self Item
   subroutine atomic_multipole_adjoint_clear(self)
      !> Item
      class(atomic_multipole_adjoint_response_type), intent(inout) :: self

      if (allocated(self%dg_dq)) deallocate (self%dg_dq)
      if (allocated(self%dg_dmu)) deallocate (self%dg_dmu)
      if (allocated(self%dg_dtheta)) deallocate (self%dg_dtheta)

   end subroutine atomic_multipole_adjoint_clear

   !* ============================================================================== *!
   !*                          Atomic charge adjoint item                            *!
   !* ============================================================================== *!

   !> Name of the atomic charge adjoint item
   !>
   !> @param[in] self Item
   function atomic_charge_adjoint_name(self) result(name)
      !> Item
      class(atomic_charge_adjoint_response_type), intent(in) :: self
      !> Name
      character(len=response_name_len) :: name

      name = "atomic_charge_adjoint"

   end function atomic_charge_adjoint_name

   !> Declare the charge weights
   !>
   !> @param[in]    self   Item
   !> @param[inout] query  Field walker
   subroutine atomic_charge_adjoint_list_fields(self, query)
      !> Item
      class(atomic_charge_adjoint_response_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query

      call query%add_real("dg_dq", "dE/dq of the atomic charges (natom)", self%dg_dq)

   end subroutine atomic_charge_adjoint_list_fields

   !> Add another atomic charge adjoint item into this one
   !>
   !> @param[in,out] self Accumulator
   !> @param[in] other Item to add
   !> @param[out] error Error handling
   subroutine atomic_charge_adjoint_add(self, other, error)
      !> Accumulator
      class(atomic_charge_adjoint_response_type), intent(inout) :: self
      !> Item to add
      class(response_item_type), intent(in) :: other
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select type (other)
      type is (atomic_charge_adjoint_response_type)
         call accumulate_vector(self%dg_dq, other%dg_dq, "atomic_charge_adjoint", "dg_dq", error)
      class default
         call type_mismatch("atomic_charge_adjoint", error)
      end select

   end subroutine atomic_charge_adjoint_add

   !> Deallocate the charge weights
   !>
   !> @param[in,out] self Item
   subroutine atomic_charge_adjoint_clear(self)
      !> Item
      class(atomic_charge_adjoint_response_type), intent(inout) :: self

      if (allocated(self%dg_dq)) deallocate (self%dg_dq)

   end subroutine atomic_charge_adjoint_clear

   !* ============================================================================== *!
   !*                         Radial potential adjoint item                          *!
   !* ============================================================================== *!

   !> Name of the radial potential adjoint item
   !>
   !> @param[in] self Item
   function radial_potential_adjoint_name(self) result(name)
      !> Item
      class(radial_potential_adjoint_response_type), intent(in) :: self
      !> Name
      character(len=response_name_len) :: name

      name = "radial_potential_adjoint"

   end function radial_potential_adjoint_name

   !> Declare the radial potential weights
   !>
   !> @param[in]    self   Item
   !> @param[inout] query  Field walker
   subroutine radial_potential_adjoint_list_fields(self, query)
      !> Item
      class(radial_potential_adjoint_response_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query

      call query%add_real2("dg_dphi", "dE/dphi of the site-resolved radial potential (ngrid, natom)", &
         & self%dg_dphi)

   end subroutine radial_potential_adjoint_list_fields

   !> Add another radial potential adjoint item into this one
   !>
   !> @param[in,out] self Accumulator
   !> @param[in] other Item to add
   !> @param[out] error Error handling
   subroutine radial_potential_adjoint_add(self, other, error)
      !> Accumulator
      class(radial_potential_adjoint_response_type), intent(inout) :: self
      !> Item to add
      class(response_item_type), intent(in) :: other
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select type (other)
      type is (radial_potential_adjoint_response_type)
         call accumulate_matrix(self%dg_dphi, other%dg_dphi, "radial_potential_adjoint", "dg_dphi", error)
      class default
         call type_mismatch("radial_potential_adjoint", error)
      end select

   end subroutine radial_potential_adjoint_add

   !> Deallocate the radial potential weights
   !>
   !> @param[in,out] self Item
   subroutine radial_potential_adjoint_clear(self)
      !> Item
      class(radial_potential_adjoint_response_type), intent(inout) :: self

      if (allocated(self%dg_dphi)) deallocate (self%dg_dphi)

   end subroutine radial_potential_adjoint_clear

   !* ============================================================================== *!
   !*                              Placeholder item                                  *!
   !* ============================================================================== *!

   !> Name of the placeholder returned outside a `next()` window
   !>
   !> @param[in] self Placeholder to describe
   function no_item_name(self) result(name)
      !> Placeholder to describe
      class(no_response_item_type), intent(in) :: self
      !> Diagnostic name
      character(len=response_name_len) :: name

      name = "no_current_item"

   end function no_item_name

   !> The placeholder holds no arrays, so it declares none
   !>
   !> @param[in]    self   Placeholder
   !> @param[inout] query  Field walker
   subroutine no_item_list_fields(self, query)
      !> Placeholder
      class(no_response_item_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query

   end subroutine no_item_list_fields

   !> The placeholder stands for no item, so nothing accumulates into it
   !>
   !> @param[inout] self  Placeholder
   !> @param[in]    other Item to add
   !> @param[out]   error Always set
   subroutine no_item_add(self, other, error)
      !> Placeholder
      class(no_response_item_type), intent(inout) :: self
      !> Item to add
      class(response_item_type), intent(in) :: other
      !> Always set
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "The placeholder response item '"//trim(self%name())// &
         & "' cannot accumulate '"//trim(other%name())//"'")

   end subroutine no_item_add

   !> The placeholder holds no arrays, so there is nothing to deallocate
   subroutine no_item_clear(self)
      !> Placeholder
      class(no_response_item_type), intent(inout) :: self

   end subroutine no_item_clear

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

   !> Advance to the next item of the response
   !>
   !> - each item is visited once per pass, in accumulation order
   !> - false ends the pass and rewinds, so the next call starts a new one; an
   !>   empty response returns false at once
   !> - a pass left early resumes here; `accumulate` and `clear` reset the
   !>   cursor
   !> - moves the cursor, so it must stand alone in a loop condition
   !>
   !> @param[inout] self Response to walk
   function response_next(self) result(more)
      !> Response to walk
      class(response_type), intent(inout) :: self
      !> Whether an item is now current
      logical :: more

      more = self%cursor < self%n()
      if (more) then
         self%cursor = self%cursor + 1
      else
         self%cursor = 0
      end if

   end function response_next

   !> Copy of the current item; select on its dynamic type to contract it
   !>
   !> Outside a `next()` window it is a placeholder named "no_current_item" with no
   !> arrays, which falls through to `class default`
   !>
   !> @param[in] self Response to inspect
   function response_item(self) result(item)
      !> Response to inspect
      class(response_type), intent(in) :: self
      !> Independent item value
      class(response_item_type), allocatable :: item

      if (self%cursor < 1 .or. self%cursor > self%n()) then
         allocate (no_response_item_type :: item)
         return
      end if
      allocate (item, source=self%items(self%cursor)%item)

   end function response_item

   !> Copy the current item, or report that none is current
   !>
   !> The checked form of `response%item()` for the C layer
   !>
   !> @param[in]  response Response to inspect
   !> @param[out] item     Independent item value
   !> @param[out] error    No current item
   subroutine current_response_item(response, item, error)
      !> Response to inspect
      class(response_type), intent(in) :: response
      !> Independent item value
      class(response_item_type), allocatable, intent(out) :: item
      !> No current item
      type(error_type), allocatable, intent(out) :: error

      if (response%cursor < 1 .or. response%cursor > response%n()) then
         call fatal_error(error, no_current_item)
         return
      end if
      allocate (item, source=response%items(response%cursor)%item)

   end subroutine current_response_item

   !> Add an item into the stored one of the same type, or append a copy
   !>
   !> Appending moves the existing slots and copies only the new item (the
   !> gfortran safe pattern for polymorphic slot arrays)
   !>
   !> @param[inout] response  Response
   !> @param[in]    item  Item to accumulate
   !> @param[out]   error Error handling
   subroutine response_accumulate(response, item, error)
      !> Response
      class(response_type), intent(inout) :: response
      !> Item to accumulate
      class(response_item_type), intent(in) :: item
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Grown slot array
      type(response_slot), allocatable :: grown(:)
      integer :: i, n, stat

      response%cursor = 0
      n = response%n()
      do i = 1, n
         if (same_type_as(response%items(i)%item, item)) then
            call response%items(i)%item%add(item, error)
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
         call move_alloc(response%items(i)%item, grown(i)%item)
      end do
      call move_alloc(grown, response%items)

   end subroutine response_accumulate

   !> Drop every item
   subroutine response_clear(response)
      !> Response
      class(response_type), intent(inout) :: response

      response%cursor = 0
      if (allocated(response%items)) deallocate (response%items)
      allocate (response%items(0))

   end subroutine response_clear

end module moist_channels_response
