!> Solvent electrostatics: table charges and solvent-model charges
!>
!> - Solvent tail charges with kernel-built point-charge field; no host data
!> - Missing table entry, custom solvent or non-water SPC/E structure: uncovered atoms
!> - Later charge terms fill uncovered atoms
!> - `coulomb_<model>_<environment>`: `moist_data_solvents` charges keyed by solvent id
!> - Charge models: hirshfeld, resp, mbis, chelpg
!> - Environments: gas, solvent, conductor
!> - `coulomb_spce`: SPC/E charges for neutral water only
!> - `multipoles_mbis_<environment>`: stub for atom-centred multipoles in 6D MOZ, no kernel yet
module moist_model_moz_potential_electrostatic_solvent
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_utils, only: to_lower
   use moist_channels_coupling, only: coupling_view_type
   use moist_channels_response, only: response_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_data_solvents, only: get_solvent_charges, get_solvent_id
   use moist_model_moz_potential_term, only: potential_term_type, potential_name_len, charge_label_len
   implicit none(type, external)
   private

   public :: solvent_charges_type, solvent_multipoles_type, solvent_model_charges_type
   public :: solvent_choice_len, charges_set_spce
   public :: coulomb_hirshfeld_gas, coulomb_hirshfeld_solvent, coulomb_hirshfeld_conductor
   public :: coulomb_resp_gas, coulomb_resp_solvent, coulomb_resp_conductor
   public :: coulomb_mbis_gas, coulomb_mbis_solvent, coulomb_mbis_conductor
   public :: coulomb_chelpg_gas, coulomb_chelpg_solvent, coulomb_chelpg_conductor
   public :: coulomb_spce, is_water
   public :: multipoles_mbis_gas, multipoles_mbis_solvent, multipoles_mbis_conductor

   !> Length of an environment or charge-model choice
   integer, parameter :: solvent_choice_len = 16

   !> Atomic partial charges of a table solvent
   type, extends(potential_term_type) :: solvent_charges_type
      !> Environment of the table: gas, solvent or conductor; blank until constructed
      character(len=solvent_choice_len) :: environment = ""
      !> Charge model of the table: hirshfeld, resp, mbis or chelpg; blank until constructed
      character(len=solvent_choice_len) :: model = ""
      !> Charges per atom, e, after update; zero when the table has none for the solvent
      real(wp), allocatable :: q(:)
      !> Whether the table had charges for the solvent of the latest update
      logical :: found = .false.
   contains
      procedure :: name => solvent_charges_name
      procedure :: update => solvent_charges_update
      procedure :: tail_charges => solvent_charges_tail_charges
      procedure :: charges_covered => solvent_charges_covered
      procedure :: charge_label => solvent_charges_label
      !> Constant charges: no adjoint to push
      procedure :: accumulate_adjoint_1d => solvent_charges_adjoint_1d
      !> Constant charges: no adjoint to push
      procedure :: accumulate_adjoint_3d => solvent_charges_adjoint_3d
   end type solvent_charges_type

   !> SPC/E water charges, Berendsen et al., J. Phys. Chem. 91, 6269 (1987)
   integer, parameter :: charges_set_spce = 1
   !> SPC/E oxygen charge, e
   real(wp), parameter :: spce_q_ow = -0.8476_wp
   !> SPC/E hydrogen charge, e
   real(wp), parameter :: spce_q_hw = 0.4238_wp

   !> Fixed charges of a solvent model
   type, extends(potential_term_type) :: solvent_model_charges_type
      !> Solvent model, `charges_set_spce`; 0 until chosen
      integer :: set = 0
      !> Charges per atom, e, after update; zero when the model does not describe the structure
      real(wp), allocatable :: q(:)
      !> Whether the model described the structure of the latest update
      logical :: found = .false.
   contains
      procedure :: charges_covered => solvent_model_charges_covered
      procedure :: charge_label => solvent_model_charges_label
      procedure :: name => solvent_model_charges_name
      procedure :: update => solvent_model_charges_update
      procedure :: tail_charges => solvent_model_charges_tail_charges
      !> Constant charges: no adjoint to push
      procedure :: accumulate_adjoint_1d => solvent_model_charges_adjoint_1d
      !> Constant charges: no adjoint to push
      procedure :: accumulate_adjoint_3d => solvent_model_charges_adjoint_3d
   end type solvent_model_charges_type

   !> Atom-centred multipoles of a table solvent; TODO for the 6D MOZ, errors at update
   type, extends(potential_term_type) :: solvent_multipoles_type
      !> Environment of the table: gas, solvent or conductor; blank until constructed
      character(len=solvent_choice_len) :: environment = ""
      !> Multipole model of the table; blank until constructed
      character(len=solvent_choice_len) :: model = ""
   contains
      procedure :: name => solvent_multipoles_name
      procedure :: update => solvent_multipoles_update
   end type solvent_multipoles_type

   !> Environments of the tables
   character(len=solvent_choice_len), parameter :: environments(3) = [character(len=solvent_choice_len) :: &
      & "gas", "solvent", "conductor"]
   !> Charge models of the table
   character(len=solvent_choice_len), parameter :: charge_models(4) = [character(len=solvent_choice_len) :: &
      & "hirshfeld", "resp", "mbis", "chelpg"]
   !> Multipole models of the table
   character(len=solvent_choice_len), parameter :: multipole_models(1) = [character(len=solvent_choice_len) :: "mbis"]

   !> Hirshfeld table charges, gas phase
   type(solvent_charges_type), parameter :: coulomb_hirshfeld_gas = &
      & solvent_charges_type(environment="gas", model="hirshfeld")
   !> Hirshfeld table charges, CPCM with the solvent's dielectric
   type(solvent_charges_type), parameter :: coulomb_hirshfeld_solvent = &
      & solvent_charges_type(environment="solvent", model="hirshfeld")
   !> Hirshfeld table charges, CPCM conductor
   type(solvent_charges_type), parameter :: coulomb_hirshfeld_conductor = &
      & solvent_charges_type(environment="conductor", model="hirshfeld")
   !> RESP table charges, gas phase
   type(solvent_charges_type), parameter :: coulomb_resp_gas = &
      & solvent_charges_type(environment="gas", model="resp")
   !> RESP table charges, CPCM with the solvent's dielectric
   type(solvent_charges_type), parameter :: coulomb_resp_solvent = &
      & solvent_charges_type(environment="solvent", model="resp")
   !> RESP table charges, CPCM conductor
   type(solvent_charges_type), parameter :: coulomb_resp_conductor = &
      & solvent_charges_type(environment="conductor", model="resp")
   !> MBIS table charges, gas phase
   type(solvent_charges_type), parameter :: coulomb_mbis_gas = &
      & solvent_charges_type(environment="gas", model="mbis")
   !> MBIS table charges, CPCM with the solvent's dielectric
   type(solvent_charges_type), parameter :: coulomb_mbis_solvent = &
      & solvent_charges_type(environment="solvent", model="mbis")
   !> MBIS table charges, CPCM conductor
   type(solvent_charges_type), parameter :: coulomb_mbis_conductor = &
      & solvent_charges_type(environment="conductor", model="mbis")
   !> CHELPG table charges, gas phase
   type(solvent_charges_type), parameter :: coulomb_chelpg_gas = &
      & solvent_charges_type(environment="gas", model="chelpg")
   !> CHELPG table charges, CPCM with the solvent's dielectric
   type(solvent_charges_type), parameter :: coulomb_chelpg_solvent = &
      & solvent_charges_type(environment="solvent", model="chelpg")
   !> CHELPG table charges, CPCM conductor
   type(solvent_charges_type), parameter :: coulomb_chelpg_conductor = &
      & solvent_charges_type(environment="conductor", model="chelpg")
   !> SPC/E water charges
   type(solvent_model_charges_type), parameter :: coulomb_spce = solvent_model_charges_type(set=charges_set_spce)
   !> MBIS table multipoles, gas phase; stub
   type(solvent_multipoles_type), parameter :: multipoles_mbis_gas = &
      & solvent_multipoles_type(environment="gas", model="mbis")
   !> MBIS table multipoles, CPCM with the solvent's dielectric; stub
   type(solvent_multipoles_type), parameter :: multipoles_mbis_solvent = &
      & solvent_multipoles_type(environment="solvent", model="mbis")
   !> MBIS table multipoles, CPCM conductor; stub
   type(solvent_multipoles_type), parameter :: multipoles_mbis_conductor = &
      & solvent_multipoles_type(environment="conductor", model="mbis")

contains

   !> Refuse a choice outside the allowed list, case-insensitively
   !>
   !> @param[in]  choice   given choice
   !> @param[in]  allowed  allowed choices, lowercase
   !> @param[in]  what     name of the choice for the message
   !> @param[out] error    unknown choice
   subroutine check_choice(choice, allowed, what, error)
      !> Given choice
      character(len=*), intent(in) :: choice
      !> Allowed choices, lowercase
      character(len=*), intent(in) :: allowed(:)
      !> Name of the choice for the message
      character(len=*), intent(in) :: what
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: list
      integer :: i
      if (any(to_lower(trim(adjustl(choice))) == allowed)) return
      list = ""
      do i = 1, size(allowed)
         if (i > 1) list = list//", "
         list = list//trim(allowed(i))
      end do
      call fatal_error(error, "Unknown solvent "//what//" '"//trim(choice)//"'; use one of "//list)
   end subroutine check_choice

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function solvent_charges_name(self) result(name)
      !> Term
      class(solvent_charges_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "solvent_charges"
   end function solvent_charges_name

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function solvent_multipoles_name(self) result(name)
      !> Term
      class(solvent_multipoles_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "solvent_multipoles"
   end function solvent_multipoles_name

   !> Read the charges of the solvent from the table
   !>
   !> - No table charges: atoms left uncovered
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         structure of this side
   !> @param[out]    error       unknown choice or a structure of another size than the table entry
   !> @param[in]     solvent_id  solvent id of the table; required
   subroutine solvent_charges_update(self, mol, error, solvent_id)
      !> Term
      class(solvent_charges_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Solvent id of the table
      integer, intent(in), optional :: solvent_id
      real(wp), allocatable :: q(:)
      logical :: table_solvent
      type(error_type), allocatable :: missing
      if (allocated(self%q)) deallocate (self%q)
      self%found = .false.
      call check_choice(self%environment, environments, "environment", error)
      if (allocated(error)) return
      call check_choice(self%model, charge_models, "charge model", error)
      if (allocated(error)) return
      allocate (self%q(mol%nat), source=0.0_wp)
      ! Solute, custom solvent or missing table charges: uncovered atoms
      table_solvent = present(solvent_id)
      if (table_solvent) table_solvent = solvent_id > 0
      if (.not. table_solvent) return
      call get_solvent_charges(solvent_id, self%environment, self%model, q, missing)
      if (allocated(missing)) return
      if (size(q) /= mol%nat) then
         call fatal_error(error, "Solvent charges of the table do not match the number of atoms of the solvent")
         deallocate (self%q)
         return
      end if
      self%q = q
      self%found = .true.
   end subroutine solvent_charges_update

   !> Refuse construction of solvent multipoles without a MOZ kernel
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         structure of this side
   !> @param[out]    error       always: not implemented
   !> @param[in]     solvent_id  solvent id of the table; required
   subroutine solvent_multipoles_update(self, mol, error, solvent_id)
      !> Term
      class(solvent_multipoles_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Solvent id of the table
      integer, intent(in), optional :: solvent_id
      call check_choice(self%environment, environments, "environment", error)
      if (allocated(error)) return
      call check_choice(self%model, multipole_models, "multipole model", error)
      if (allocated(error)) return
      call require_table_solvent("Solvent multipoles", solvent_id, error)
      if (allocated(error)) return
      ! TODO: atom-centred multipoles for the 6D MOZ (moist_data_solvents: get_solvent_multipoles)
      call fatal_error(error, "Solvent multipoles are not implemented yet; they are needed for the 6D MOZ")
   end subroutine solvent_multipoles_update

   !> Refuse a side without a solvent id of the table
   !>
   !> @param[in]  what        term name for the message
   !> @param[in]  solvent_id  optional solvent id of the table
   !> @param[out] error       absent or non-positive id
   subroutine require_table_solvent(what, solvent_id, error)
      !> Term name for the message
      character(len=*), intent(in) :: what
      !> Optional solvent id of the table
      integer, intent(in), optional :: solvent_id
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      logical :: ok
      ok = present(solvent_id)
      if (ok) ok = solvent_id > 0
      if (.not. ok) call fatal_error(error, what//" need a solvent of the table; a custom solvent or a solute "// &
         & "takes fixed charges instead")
   end subroutine require_table_solvent

   !> Read the table tail charges
   !>
   !> @param[in]  self   built term
   !> @param[in]  view   scoped coupling view of this term; unused
   !> @param[out] q      tail charges per atom, e
   !> @param[out] error  unbuilt term
   subroutine solvent_charges_tail_charges(self, view, q, error)
      !> Built term
      class(solvent_charges_type), intent(in) :: self
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Tail charges per atom, e
      real(wp), allocatable, intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(self%q)) then
         call fatal_error(error, "Solvent charges are not built")
         return
      end if
      q = self%q
   end subroutine solvent_charges_tail_charges

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function solvent_model_charges_name(self) result(name)
      !> Term
      class(solvent_model_charges_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "solvent_model_charges"
   end function solvent_model_charges_name

   !> Assign the solvent-model charges by element
   !>
   !> - Non-neutral or non-water structure: atoms left for later charge terms
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         structure of this side
   !> @param[out]    error       unknown model or missing atom identities
   !> @param[in]     solvent_id  optional solvent id of the table
   subroutine solvent_model_charges_update(self, mol, error, solvent_id)
      !> Term
      class(solvent_model_charges_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional solvent id of the table
      integer, intent(in), optional :: solvent_id
      character(len=16) :: label
      if (allocated(self%q)) deallocate (self%q)
      self%found = .false.
      if (self%set /= charges_set_spce) then
         write (label, "(i0)") self%set
         call fatal_error(error, "Unknown solvent charge model "//trim(label)//"; use coulomb_spce")
         return
      end if
      if (.not. allocated(mol%id) .or. .not. allocated(mol%num)) then
         call fatal_error(error, "SPC/E charges need initialized atom identities")
         return
      end if
      self%found = is_water(mol, solvent_id)
      if (self%found) then
         self%q = merge(spce_q_ow, spce_q_hw, mol%num(mol%id) == 8)
      else
         allocate (self%q(mol%nat), source=0.0_wp)
      end if
   end subroutine solvent_model_charges_update

   !> Read the solvent-model tail charges
   !>
   !> @param[in]  self   built term
   !> @param[in]  view   scoped coupling view of this term; unused
   !> @param[out] q      tail charges per atom, e
   !> @param[out] error  unbuilt term
   subroutine solvent_model_charges_tail_charges(self, view, q, error)
      !> Built term
      class(solvent_model_charges_type), intent(in) :: self
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Tail charges per atom, e
      real(wp), allocatable, intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(self%q)) then
         call fatal_error(error, "Solvent-model charges are not built")
         return
      end if
      q = self%q
   end subroutine solvent_model_charges_tail_charges

   !> Charge coverage, every atom with table charges, none otherwise
   !>
   !> @param[in] self   built term
   !> @param[in] natom  number of atoms of the side
   pure function solvent_charges_covered(self, natom) result(covered)
      !> Built term
      class(solvent_charges_type), intent(in) :: self
      !> Number of atoms of the side
      integer, intent(in) :: natom
      !> Whether each atom has a charge from this term (natom)
      logical :: covered(natom)
      covered = self%found
   end function solvent_charges_covered

   !> Charge-source label, model/environment
   !>
   !> @param[in] self  term
   pure function solvent_charges_label(self) result(label)
      !> Term
      class(solvent_charges_type), intent(in) :: self
      !> Label
      character(len=charge_label_len) :: label
      label = trim(self%model)//"/"//trim(self%environment)
   end function solvent_charges_label

   !> Charge coverage, every atom for modelled structures, none otherwise
   !>
   !> @param[in] self   built term
   !> @param[in] natom  number of atoms of the side
   pure function solvent_model_charges_covered(self, natom) result(covered)
      !> Built term
      class(solvent_model_charges_type), intent(in) :: self
      !> Number of atoms of the side
      integer, intent(in) :: natom
      !> Whether each atom has a charge from this term (natom)
      logical :: covered(natom)
      covered = self%found
   end function solvent_model_charges_covered

   !> Charge-source label
   !>
   !> @param[in] self  term
   pure function solvent_model_charges_label(self) result(label)
      !> Term
      class(solvent_model_charges_type), intent(in) :: self
      !> Label
      character(len=charge_label_len) :: label
      label = "SPC/E"
   end function solvent_model_charges_label

   !> Leave radial adjoints unchanged; table charges independent of geometry and host data
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      radial grid
   !> @param[in]     view      scoped coupling view of this term; unused
   !> @param[in]     w_phi     adjoint of the radial potential (npts, natom); unused
   !> @param[in]     w_qinf    adjoint of the tail charges (natom); unused
   !> @param[in,out] response  host response items; untouched
   !> @param[out]    error     never set
   subroutine solvent_charges_adjoint_1d(self, grid, view, w_phi, w_qinf, response, error)
      !> Term
      class(solvent_charges_type), intent(inout) :: self
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
   end subroutine solvent_charges_adjoint_1d

   !> Leave volume adjoints unchanged; table charges independent of geometry and host data
   !>
   !> @param[in,out] self          term
   !> @param[in]     grid          volume grid
   !> @param[in]     view          scoped coupling view of this term; unused
   !> @param[in]     w_phi         adjoint of the grid potential (ngrid); unused
   !> @param[in]     w_qinf        adjoint of the tail charges (natom); unused
   !> @param[in,out] response      host response items; untouched
   !> @param[out]    error         never set
   !> @param[in,out] grid_adjoint  optional volume adjoint accumulator; untouched
   !> @param[in,out] gradient      optional explicit nuclear gradient, Hartree/bohr (3, natom)
   subroutine solvent_charges_adjoint_3d(self, grid, view, w_phi, w_qinf, response, error, grid_adjoint, gradient)
      !> Term
      class(solvent_charges_type), intent(inout) :: self
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
   end subroutine solvent_charges_adjoint_3d

   !> Leave radial adjoints unchanged; solvent-model charges independent of geometry and host data
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      radial grid
   !> @param[in]     view      scoped coupling view of this term; unused
   !> @param[in]     w_phi     adjoint of the radial potential (npts, natom); unused
   !> @param[in]     w_qinf    adjoint of the tail charges (natom); unused
   !> @param[in,out] response  host response items; untouched
   !> @param[out]    error     never set
   subroutine solvent_model_charges_adjoint_1d(self, grid, view, w_phi, w_qinf, response, error)
      !> Term
      class(solvent_model_charges_type), intent(inout) :: self
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
   end subroutine solvent_model_charges_adjoint_1d

   !> Leave volume adjoints unchanged; solvent-model charges independent of geometry and host data
   !>
   !> @param[in,out] self          term
   !> @param[in]     grid          volume grid
   !> @param[in]     view          scoped coupling view of this term; unused
   !> @param[in]     w_phi         adjoint of the grid potential (ngrid); unused
   !> @param[in]     w_qinf        adjoint of the tail charges (natom); unused
   !> @param[in,out] response      host response items; untouched
   !> @param[out]    error         never set
   !> @param[in,out] grid_adjoint  optional volume adjoint accumulator; untouched
   !> @param[in,out] gradient      optional explicit nuclear gradient, Hartree/bohr (3, natom)
   subroutine solvent_model_charges_adjoint_3d(self, grid, view, w_phi, w_qinf, response, error, grid_adjoint, gradient)
      !> Term
      class(solvent_model_charges_type), intent(inout) :: self
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
   end subroutine solvent_model_charges_adjoint_3d

   !> Whether the structure is one neutral, closed-shell water molecule
   !>
   !> - Non-water table solvent id: false
   !>
   !> @param[in] mol         structure with initialized atom identities
   !> @param[in] solvent_id  optional solvent id of the table; 0 for a custom solvent
   function is_water(mol, solvent_id) result(water)
      !> Structure
      class(structure_type), intent(in) :: mol
      !> Optional solvent id of the table
      integer, intent(in), optional :: solvent_id
      !> Whether the structure is neutral H2O
      logical :: water
      type(error_type), allocatable :: error
      integer :: water_id
      water = .false.
      if (present(solvent_id)) then
         if (solvent_id /= 0) then
            call get_solvent_id("water", water_id, error)
            if (allocated(error)) return
            if (solvent_id /= water_id) return
         end if
      end if
      if (mol%nat /= 3) return
      if (count(mol%num(mol%id) == 8) /= 1 .or. count(mol%num(mol%id) == 1) /= 2) return
      if (.not. ieee_is_finite(mol%charge) .or. abs(mol%charge) > 1.0e-8_wp .or. mol%uhf /= 0) return
      water = .true.
   end function is_water

end module moist_model_moz_potential_electrostatic_solvent
