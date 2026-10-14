!> Fixed per-atom point charges of the MOZ potential layer
!>
!> - Explicit charges at construction, for either side
!> - Tail charges and kernel-built point-charge field
!> - Optional mask for atoms left to later charge terms
!> - No host data; all derivatives from kernels
!> - Explicit no-op adjoint for constant charges
module moist_model_moz_potential_electrostatic_fixed
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use moist_channels_coupling, only: coupling_view_type
   use moist_channels_response, only: response_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_model_moz_potential_term, only: potential_term_type, potential_name_len, charge_label_len
   implicit none(type, external)
   private

   public :: fixed_charges_type, new_fixed_charges

   !> Potential term of fixed per-atom point charges
   type, extends(potential_term_type) :: fixed_charges_type
      !> Charges per atom, e
      real(wp), allocatable :: q(:)
      !> Whether each atom takes its charge from this term (natom)
      logical, allocatable :: covered(:)
   contains
      procedure :: name => fixed_charges_name
      procedure :: update => fixed_charges_update
      procedure :: tail_charges => fixed_charges_tail_charges
      procedure :: charges_covered => fixed_charges_covered
      procedure :: charge_label => fixed_charges_label
      !> Constant charges: no adjoint to push
      procedure :: accumulate_adjoint_1d => fixed_charges_adjoint_1d
      !> Constant charges: no adjoint to push
      procedure :: accumulate_adjoint_3d => fixed_charges_adjoint_3d
   end type fixed_charges_type

contains

   !> Store per-atom charges
   !>
   !> @param[out] self     term
   !> @param[in]  q        charges per atom, e
   !> @param[out] error    empty, non-finite or unmasked charges
   !> @param[in]  covered  optional coverage mask; every atom when absent
   subroutine new_fixed_charges(self, q, error, covered)
      !> Term
      type(fixed_charges_type), intent(out) :: self
      !> Charges per atom, e
      real(wp), intent(in) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional coverage mask; every atom when absent
      logical, intent(in), optional :: covered(:)
      if (size(q) < 1) then
         call fatal_error(error, "Fixed charges cover no atom")
         return
      end if
      if (present(covered)) then
         if (size(covered) /= size(q)) then
            call fatal_error(error, "Fixed charges coverage mask does not match the charges")
            return
         end if
         self%covered = covered
      else
         allocate (self%covered(size(q)), source=.true.)
      end if
      if (.not. all(ieee_is_finite(q) .or. .not. self%covered)) then
         call fatal_error(error, "Fixed charges must be finite for every covered atom")
         return
      end if
      self%q = merge(q, 0.0_wp, self%covered)
   end subroutine new_fixed_charges

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function fixed_charges_name(self) result(name)
      !> Term
      class(fixed_charges_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "fixed_charges"
   end function fixed_charges_name

   !> Check the stored charges against the structure
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         structure of this side
   !> @param[out]    error       unconstructed term or charges not matching the structure
   !> @param[in]     solvent_id  unused; the charges are explicit
   subroutine fixed_charges_update(self, mol, error, solvent_id)
      !> Term
      class(fixed_charges_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Unused; the charges are explicit
      integer, intent(in), optional :: solvent_id
      if (.not. allocated(self%q)) then
         call fatal_error(error, "Fixed charges are not set; construct the term with new_fixed_charges")
         return
      end if
      if (size(self%q) /= mol%nat) then
         call fatal_error(error, "Fixed charges do not match the number of atoms of the structure")
      end if
   end subroutine fixed_charges_update

   !> Read the stored tail charges
   !>
   !> @param[in]  self   term
   !> @param[in]  view   scoped coupling view of this term; unused
   !> @param[out] q      tail charges per atom, e
   !> @param[out] error  unconstructed term
   subroutine fixed_charges_tail_charges(self, view, q, error)
      !> Term
      class(fixed_charges_type), intent(in) :: self
      !> Scoped coupling view of this term
      type(coupling_view_type), intent(in) :: view
      !> Tail charges per atom, e
      real(wp), allocatable, intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(self%q)) then
         call fatal_error(error, "Fixed charges are not set; construct the term with new_fixed_charges")
         return
      end if
      q = self%q
   end subroutine fixed_charges_tail_charges

   !> Atoms whose charges this term supplies
   !>
   !> @param[in] self   built term
   !> @param[in] natom  number of atoms of the side
   pure function fixed_charges_covered(self, natom) result(covered)
      !> Built term
      class(fixed_charges_type), intent(in) :: self
      !> Number of atoms of the side
      integer, intent(in) :: natom
      !> Whether each atom has a charge from this term (natom)
      logical :: covered(natom)
      covered = .false.
      if (allocated(self%covered)) then
         if (size(self%covered) == natom) covered = self%covered
      end if
   end function fixed_charges_covered

   !> Charge-source label
   !>
   !> @param[in] self  term
   pure function fixed_charges_label(self) result(label)
      !> Term
      class(fixed_charges_type), intent(in) :: self
      !> Label
      character(len=charge_label_len) :: label
      label = "fixed"
   end function fixed_charges_label

   !> Leave radial adjoints unchanged; fixed charges independent of geometry and host data
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      radial grid
   !> @param[in]     view      scoped coupling view of this term; unused
   !> @param[in]     w_phi     adjoint of the radial potential (npts, natom); unused
   !> @param[in]     w_qinf    adjoint of the tail charges (natom); unused
   !> @param[in,out] response  host response items; untouched
   !> @param[out]    error     never set
   subroutine fixed_charges_adjoint_1d(self, grid, view, w_phi, w_qinf, response, error)
      !> Term
      class(fixed_charges_type), intent(inout) :: self
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
   end subroutine fixed_charges_adjoint_1d

   !> Leave volume adjoints unchanged; fixed charges independent of geometry and host data
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
   subroutine fixed_charges_adjoint_3d(self, grid, view, w_phi, w_qinf, response, error, grid_adjoint, gradient)
      !> Term
      class(fixed_charges_type), intent(inout) :: self
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
   end subroutine fixed_charges_adjoint_3d

end module moist_model_moz_potential_electrostatic_fixed
