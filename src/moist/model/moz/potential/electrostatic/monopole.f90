!> Host-fed point monopoles of the MOZ potential layer
!>
!> - Own "charges" request: q in every phase on radial and volume grids
!> - Tail charges read from the scoped coupling view at every call
!> - Full point-charge field from kernels
!> - Charge adjoint returned as `atomic_charge_adjoint`
module moist_model_moz_potential_electrostatic_monopole
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_channels_coupling, only: coupling_type, coupling_view_type, atomic_charge_request_type, &
      & request_require, coupling_register, moist_phase_energy, moist_phase_response, moist_phase_gradient
   use moist_channels_response, only: response_type, response_accumulate, atomic_charge_adjoint_response_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_model_moz_potential_term, only: potential_term_type, potential_name_len, charge_label_len
   implicit none(type, external)
   private

   public :: monopole_type

   !> Phases requiring host charges
   integer, parameter :: all_phases(3) = [moist_phase_energy, moist_phase_response, moist_phase_gradient]

   !> Point-monopole potential with host charges
   type, extends(potential_term_type) :: monopole_type
      !> Number of atoms of the latest update; 0 until built
      integer :: natom = 0
   contains
      procedure :: name => monopole_name
      procedure :: update => monopole_update
      !> Host-fed by construction, so a solvent potential refuses it at `add`
      procedure :: is_host_fed => monopole_is_host_fed
      procedure :: declare_pass_1d => monopole_declare_pass_1d
      procedure :: declare_pass_3d => monopole_declare_pass_3d
      procedure :: tail_charges => monopole_tail_charges
      procedure :: charges_covered => monopole_covered
      procedure :: charge_label => monopole_label
      procedure :: accumulate_adjoint_1d => monopole_accumulate_adjoint_1d
      procedure :: accumulate_adjoint_3d => monopole_accumulate_adjoint_3d
   end type monopole_type

contains

   !> Host-fed by construction
   !>
   !> @param[in] self  term
   pure function monopole_is_host_fed(self) result(host_fed)
      !> Term
      class(monopole_type), intent(in) :: self
      !> Host-fed flag
      logical :: host_fed
      host_fed = .true.
   end function monopole_is_host_fed

   !> Charge coverage, every atom
   !>
   !> @param[in] self   built term
   !> @param[in] natom  number of atoms of the side
   pure function monopole_covered(self, natom) result(covered)
      !> Built term
      class(monopole_type), intent(in) :: self
      !> Number of atoms of the side
      integer, intent(in) :: natom
      !> Whether each atom has a charge from this term (natom)
      logical :: covered(natom)
      covered = self%natom == natom
   end function monopole_covered

   !> Charge-source label
   !>
   !> @param[in] self  term
   pure function monopole_label(self) result(label)
      !> Term
      class(monopole_type), intent(in) :: self
      !> Label
      character(len=charge_label_len) :: label
      label = "host"
   end function monopole_label

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function monopole_name(self) result(name)
      !> Term
      class(monopole_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "monopole"
   end function monopole_name

   !> Record the number of atoms; the charges come from the coupling
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         structure of this side
   !> @param[out]    error       never set
   !> @param[in]     solvent_id  unused; host terms serve a solute
   subroutine monopole_update(self, mol, error, solvent_id)
      !> Term
      class(monopole_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Unused; host terms serve a solute
      integer, intent(in), optional :: solvent_id
      self%natom = mol%nat
   end subroutine monopole_update

   !> Declare the partial charges for a radial grid, q in every phase
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      radial grid
   !> @param[in,out] coupling  coupling, scope already set by the potential
   !> @param[out]    error     declaration failure
   subroutine monopole_declare_pass_1d(self, grid, coupling, error)
      !> Term
      class(monopole_type), intent(inout) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Coupling, scope already set by the potential
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call declare_charges(coupling, error)
   end subroutine monopole_declare_pass_1d

   !> Declare the partial charges for a volume grid, q in every phase
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      volume grid
   !> @param[in,out] coupling  coupling, scope already set by the potential
   !> @param[out]    error     declaration failure
   subroutine monopole_declare_pass_3d(self, grid, coupling, error)
      !> Term
      class(monopole_type), intent(inout) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Coupling, scope already set by the potential
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call declare_charges(coupling, error)
   end subroutine monopole_declare_pass_3d

   !> Read the host charges, the tail charges of the term
   !>
   !> @param[in]  self   built term
   !> @param[in]  view   scoped coupling view of this term
   !> @param[out] q      tail charges per atom, e (natom)
   !> @param[out] error  unbuilt term, missing output or wrong length
   subroutine monopole_tail_charges(self, view, q, error)
      !> Built term
      class(monopole_type), intent(in) :: self
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Tail charges per atom, e
      real(wp), allocatable, intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (self%natom < 1) then
         call fatal_error(error, "Monopole term '"//trim(self%name())//"' is not built")
         return
      end if
      call view%read("charges", "q", q, error)
      if (allocated(error)) return
      if (size(q) /= self%natom) then
         deallocate (q)
         call fatal_error(error, "Host charges do not match the number of solute atoms")
      end if
   end subroutine monopole_tail_charges

   !> Push the adjoint of the charges into the response (1D)
   !>
   !> @param[in,out] self      built term
   !> @param[in]     grid      radial grid
   !> @param[in]     view      scoped coupling view of this term
   !> @param[in]     w_phi     adjoint of the radial field (npts, natom); not an input of this term
   !> @param[in]     w_qinf    adjoint of the charges (natom)
   !> @param[in,out] response  host response items
   !> @param[out]    error     wrong adjoint length or accumulation failure
   subroutine monopole_accumulate_adjoint_1d(self, grid, view, w_phi, w_qinf, response, error)
      !> Built term
      class(monopole_type), intent(inout) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Adjoint of the radial field (npts, natom)
      real(wp), intent(in) :: w_phi(:, :)
      !> Adjoint of the charges (natom)
      real(wp), intent(in) :: w_qinf(:)
      !> Host response items
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call push_charge_adjoint(self%natom, w_qinf, response, error)
   end subroutine monopole_accumulate_adjoint_1d

   !> Push the adjoint of the charges into the response (3D)
   !>
   !> @param[in,out] self          built term
   !> @param[in]     grid          volume grid
   !> @param[in]     view          scoped coupling view of this term
   !> @param[in]     w_phi         adjoint of the full field (ngrid); not an input of this term
   !> @param[in]     w_qinf        adjoint of the charges (natom)
   !> @param[in,out] response      host response items
   !> @param[out]    error         wrong adjoint length or accumulation failure
   !> @param[in,out] grid_adjoint  optional volume adjoint accumulator, untouched
   !> @param[in,out] gradient      optional explicit nuclear gradient, Hartree/bohr (3, natom)
   subroutine monopole_accumulate_adjoint_3d(self, grid, view, w_phi, w_qinf, response, error, grid_adjoint, gradient)
      !> Built term
      class(monopole_type), intent(inout) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Adjoint of the full field (ngrid)
      real(wp), intent(in) :: w_phi(:)
      !> Adjoint of the charges (natom)
      real(wp), intent(in) :: w_qinf(:)
      !> Host response items
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional volume adjoint accumulator
      type(volume_adjoint_type), intent(inout), optional :: grid_adjoint
      !> Optional explicit nuclear gradient, Hartree/bohr (3, natom)
      real(wp), intent(inout), optional :: gradient(:, :)
      call push_charge_adjoint(self%natom, w_qinf, response, error)
   end subroutine monopole_accumulate_adjoint_3d

   !> Declare the host charges, q in every phase
   !>
   !> @param[in,out] coupling  coupling, scope already set by the potential
   !> @param[out]    error     declaration failure
   subroutine declare_charges(coupling, error)
      !> Coupling, scope already set by the potential
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(atomic_charge_request_type) :: charges
      integer :: ip
      do ip = 1, size(all_phases)
         call request_require(charges, all_phases(ip), "q", error)
         if (allocated(error)) return
      end do
      call coupling_register(coupling, "charges", charges, error)
   end subroutine declare_charges

   !> Push the adjoint of the host charges as an `atomic_charge_adjoint` item
   !>
   !> @param[in]     natom     number of solute atoms of the built term
   !> @param[in]     w_qinf    adjoint of the charges (natom)
   !> @param[in,out] response  host response items
   !> @param[out]    error     wrong adjoint length or accumulation failure
   subroutine push_charge_adjoint(natom, w_qinf, response, error)
      !> Number of solute atoms of the built term
      integer, intent(in) :: natom
      !> Adjoint of the charges (natom)
      real(wp), intent(in) :: w_qinf(:)
      !> Host response items
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(atomic_charge_adjoint_response_type) :: item
      integer :: stat
      if (size(w_qinf) /= natom) then
         call fatal_error(error, "Host charge adjoint does not match the number of solute atoms")
         return
      end if
      allocate (item%dg_dq, source=w_qinf, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate the host charge adjoint")
         return
      end if
      call response_accumulate(response, item, error)
   end subroutine push_charge_adjoint

end module moist_model_moz_potential_electrostatic_monopole
