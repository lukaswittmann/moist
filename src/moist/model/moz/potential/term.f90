!> Abstract term of a MOZ potential
!>
!> - Updated from a structure; solvent id for solvent-specific parameters
!> - Grid-dependent data from the owning model's grid in `update_grid`, after `update`; none by default
!> - Short range: plain-data 12-6 Lennard-Jones slot `pair`, partial atom coverage allowed
!> - Uncovered atoms filled by later terms of the potential
!> - Electrostatics: tail charges (`tail_charges`) and/or grid potential (`has_field`, `field`)
!> - Common term type for `solvent%potential%add` and `model%potential%add`
!> - Host data through coupling: `is_host_fed`; refused by a solvent potential
!> - Requests in `declare_pass`, one coupling scope per term; radial (1D) and volume (3D) generics
!> - Scoped views for `tail_charges` and `field`
!> - Adjoint dispatch to every field supplier or charge owner, including parameter terms
!> - Input adjoints through `accumulate_adjoint`; explicit no-op for constant data
!> - `field` and `accumulate_adjoint`: radial (1D) and volume (3D) generics
!> - Default field and adjoint implementations: errors
module moist_model_moz_potential_term
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_channels_coupling, only: coupling_type, coupling_view_type
   use moist_channels_response, only: response_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_model_moz_potential_lj_base, only: lj_12_6_type
   implicit none(type, external)
   private

   public :: potential_term_type, potential_name_len, charge_label_len, atom_label, atom_label_len

   !> Length of a term name
   integer, parameter :: potential_name_len = 32
   !> Length of the charge-source label of a term
   integer, parameter :: charge_label_len = 16
   !> Length of an atom label
   integer, parameter :: atom_label_len = 16

   !> Abstract potential term
   type, abstract :: potential_term_type
      !> Lennard-Jones parameters of this side; absent when the term supplies none
      type(lj_12_6_type), allocatable :: pair
   contains
      !> Diagnostic name
      procedure(potential_name_i), deferred :: name
      !> Fill the term from the structure
      procedure(potential_update_i), deferred :: update
      !> Prepare grid-dependent data for a radial grid after `update`; nothing by default
      procedure :: update_grid_1d => potential_update_grid_1d
      !> Prepare grid-dependent data for a volume grid after `update`; nothing by default
      procedure :: update_grid_3d => potential_update_grid_3d
      generic :: update_grid => update_grid_1d, update_grid_3d
      !> Whether evaluation reads host data from the coupling; false by default
      procedure :: is_host_fed => potential_is_host_fed
      !> Whether `field` supplies a grid potential; false by default
      procedure :: has_field => potential_has_field
      !> Declare coupling requests for a radial grid in the current scope; none by default
      procedure :: declare_pass_1d => potential_declare_pass_1d
      !> Declare coupling requests for a volume grid in the current scope; none by default
      procedure :: declare_pass_3d => potential_declare_pass_3d
      generic :: declare_pass => declare_pass_1d, declare_pass_3d
      !> Tail charges per atom; none (unallocated) by default
      procedure :: tail_charges => potential_tail_charges
      !> Atoms whose charges this built term supplies; none by default
      procedure :: charges_covered => potential_charges_covered
      !> Short charge-source label for the site table; blank by default
      procedure :: charge_label => potential_charge_label
      !> Radial potential per solute atom on a 1D grid; errors by default
      procedure :: field_1d => potential_field_1d
      !> Potential on a volume grid; errors by default
      procedure :: field_3d => potential_field_3d
      generic :: field => field_1d, field_3d
      !> Push field and tail-charge adjoints into host inputs (1D); errors by default
      procedure :: accumulate_adjoint_1d => potential_accumulate_adjoint_1d
      !> Push field and tail-charge adjoints into host inputs (3D); errors by default
      procedure :: accumulate_adjoint_3d => potential_accumulate_adjoint_3d
      generic :: accumulate_adjoint => accumulate_adjoint_1d, accumulate_adjoint_3d
   end type potential_term_type

   abstract interface

      !> Diagnostic name
      !>
      !> @param[in] self  term
      pure function potential_name_i(self) result(name)
         import potential_term_type, potential_name_len
         implicit none(type, external)
         !> Term
         class(potential_term_type), intent(in) :: self
         !> Name
         character(len=potential_name_len) :: name
      end function potential_name_i

      !> Fill the term from the structure
      !>
      !> - `solvent_id` required for solvent-keyed parameters, ignored by structure-based typing
      !> - Present `pair`: `mol%nat` entries
      !>
      !> @param[in,out] self        term
      !> @param[in]     mol         structure of this side
      !> @param[out]    error       typing or parameter failure
      !> @param[in]     solvent_id  optional solvent id of this side; absent for a solute
      subroutine potential_update_i(self, mol, error, solvent_id)
         import potential_term_type, structure_type, error_type
         implicit none(type, external)
         !> Term
         class(potential_term_type), intent(inout) :: self
         !> Structure of this side
         class(structure_type), intent(in) :: mol
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Optional solvent id of this side; absent for a solute
         integer, intent(in), optional :: solvent_id
      end subroutine potential_update_i

   end interface

contains

   !> Parameter terms read no host data
   !>
   !> @param[in] self  term
   pure function potential_is_host_fed(self) result(host_fed)
      !> Term
      class(potential_term_type), intent(in) :: self
      !> Host-fed flag
      logical :: host_fed
      host_fed = .false.
   end function potential_is_host_fed

   !> Grid-field presence, false by default; point-charge field from tail charges
   !>
   !> @param[in] self  term
   pure function potential_has_field(self) result(has)
      !> Term
      class(potential_term_type), intent(in) :: self
      !> Whether `field` supplies a grid potential
      logical :: has
      has = .false.
   end function potential_has_field

   !> Prepare nothing for a radial grid
   !>
   !> @param[in,out] self   updated term
   !> @param[in]     grid   radial grid of the owning model
   !> @param[out]    error  never set
   subroutine potential_update_grid_1d(self, grid, error)
      !> Updated term
      class(potential_term_type), intent(inout) :: self
      !> Radial grid of the owning model
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   end subroutine potential_update_grid_1d

   !> Prepare nothing for a volume grid
   !>
   !> @param[in,out] self   updated term
   !> @param[in]     grid   volume grid of the owning model
   !> @param[out]    error  never set
   subroutine potential_update_grid_3d(self, grid, error)
      !> Updated term
      class(potential_term_type), intent(inout) :: self
      !> Volume grid of the owning model
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   end subroutine potential_update_grid_3d

   !> Declare nothing for a radial grid
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      radial grid
   !> @param[in,out] coupling  coupling, scope already set by the potential
   !> @param[out]    error     never set
   subroutine potential_declare_pass_1d(self, grid, coupling, error)
      !> Term
      class(potential_term_type), intent(inout) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Coupling, scope already set by the potential
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   end subroutine potential_declare_pass_1d

   !> Declare nothing for a volume grid
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      volume grid
   !> @param[in,out] coupling  coupling, scope already set by the potential
   !> @param[out]    error     never set
   subroutine potential_declare_pass_3d(self, grid, coupling, error)
      !> Term
      class(potential_term_type), intent(inout) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Coupling, scope already set by the potential
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   end subroutine potential_declare_pass_3d

   !> No tail charges: `q` stays unallocated
   !>
   !> @param[in]  self   term
   !> @param[in]  view   scoped coupling view of this term
   !> @param[out] q      tail charges per atom, e; unallocated
   !> @param[out] error  never set
   subroutine potential_tail_charges(self, view, q, error)
      !> Term
      class(potential_term_type), intent(in) :: self
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Tail charges per atom, e
      real(wp), allocatable, intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   end subroutine potential_tail_charges

   !> Charge coverage, none by default
   !>
   !> @param[in] self   built term
   !> @param[in] natom  number of atoms of the side
   pure function potential_charges_covered(self, natom) result(covered)
      !> Built term
      class(potential_term_type), intent(in) :: self
      !> Number of atoms of the side
      integer, intent(in) :: natom
      !> Whether each atom has a charge from this term (natom)
      logical :: covered(natom)
      covered = .false.
   end function potential_charges_covered

   !> No charge-source label
   !>
   !> @param[in] self  term
   pure function potential_charge_label(self) result(label)
      !> Term
      class(potential_term_type), intent(in) :: self
      !> Label
      character(len=charge_label_len) :: label
      label = ""
   end function potential_charge_label

   !> Refuse a radial field
   !>
   !> @param[in]  self     term
   !> @param[in]  grid     radial grid
   !> @param[in]  view     scoped coupling view of this term
   !> @param[out] phi      radial potential per solute atom, Hartree/e (npts, natom)
   !> @param[out] error    always set
   subroutine potential_field_1d(self, grid, view, phi, error)
      !> Term
      class(potential_term_type), intent(in) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Radial potential per solute atom, Hartree/e (npts, natom)
      real(wp), allocatable, intent(out) :: phi(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "Potential term '"//trim(self%name())//"' supplies no radial field")
   end subroutine potential_field_1d

   !> Refuse a volume field
   !>
   !> @param[in]  self     term
   !> @param[in]  grid     volume grid
   !> @param[in]  view     scoped coupling view of this term
   !> @param[out] phi      potential on the grid, Hartree/e (ngrid)
   !> @param[out] error    always set
   !> @param[out] dphi_dr  optional gradient with respect to the grid point, Hartree/(e bohr) (3, ngrid)
   subroutine potential_field_3d(self, grid, view, phi, error, dphi_dr)
      !> Term
      class(potential_term_type), intent(in) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Potential on the grid, Hartree/e (ngrid)
      real(wp), allocatable, intent(out) :: phi(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional gradient with respect to the grid point, Hartree/(e bohr) (3, ngrid)
      real(wp), allocatable, intent(out), optional :: dphi_dr(:, :)
      call fatal_error(error, "Potential term '"//trim(self%name())//"' supplies no volume field")
   end subroutine potential_field_3d

   !> Refuse a radial adjoint
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      radial grid
   !> @param[in]     view      scoped coupling view of this term
   !> @param[in]     w_phi     adjoint of the radial potential (npts, natom)
   !> @param[in]     w_qinf    adjoint of the tail charges (natom); zero-sized when the term has none
   !> @param[in,out] response  host response items
   !> @param[out]    error     always set
   subroutine potential_accumulate_adjoint_1d(self, grid, view, w_phi, w_qinf, response, error)
      !> Term
      class(potential_term_type), intent(inout) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Adjoint of the radial potential (npts, natom)
      real(wp), intent(in) :: w_phi(:, :)
      !> Adjoint of the tail charges (natom)
      real(wp), intent(in) :: w_qinf(:)
      !> Host response items
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "Potential term '"//trim(self%name())//"' has no radial adjoint")
   end subroutine potential_accumulate_adjoint_1d

   !> Refuse a volume adjoint
   !>
   !> @param[in,out] self          term
   !> @param[in]     grid          volume grid
   !> @param[in]     view          scoped coupling view of this term
   !> @param[in]     w_phi         adjoint of the grid potential (ngrid)
   !> @param[in]     w_qinf        adjoint of the tail charges (natom); zero-sized when the term has none
   !> @param[in,out] response      host response items
   !> @param[out]    error         always set
   !> @param[in,out] grid_adjoint  optional volume adjoint accumulator; present when the geometry derivative is requested
   !> @param[in,out] gradient      optional explicit nuclear gradient, Hartree/bohr (3, natom)
   subroutine potential_accumulate_adjoint_3d(self, grid, view, w_phi, w_qinf, response, error, grid_adjoint, gradient)
      !> Term
      class(potential_term_type), intent(inout) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Adjoint of the grid potential (ngrid)
      real(wp), intent(in) :: w_phi(:)
      !> Adjoint of the tail charges (natom)
      real(wp), intent(in) :: w_qinf(:)
      !> Host response items
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional volume adjoint accumulator
      type(volume_adjoint_type), intent(inout), optional :: grid_adjoint
      !> Optional explicit nuclear gradient, Hartree/bohr (3, natom)
      real(wp), intent(inout), optional :: gradient(:, :)
      call fatal_error(error, "Potential term '"//trim(self%name())//"' has no volume adjoint")
   end subroutine potential_accumulate_adjoint_3d

   !> Form the diagnostic atom label "Sym<index>", e.g. "O1"
   !>
   !> - Fixed-length subroutine output, safe inside OpenMP regions
   !>
   !> @param[in]  mol    structure
   !> @param[in]  iat    atom index
   !> @param[out] label  element symbol followed by the atom index
   subroutine atom_label(mol, iat, label)
      !> Structure
      class(structure_type), intent(in) :: mol
      !> Atom index
      integer, intent(in) :: iat
      !> Element symbol followed by the atom index
      character(len=atom_label_len), intent(out) :: label
      write (label, "(a, i0)") trim(mol%sym(mol%id(iat))), iat
   end subroutine atom_label

end module moist_model_moz_potential_term
