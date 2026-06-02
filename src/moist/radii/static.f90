module moist_radii_static
   use mctc_env, only: wp
   use mctc_env, only: error_type, fatal_error
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_io, only: structure_type
   use mctc_io_convert, only: autoaa
   use mctc_io_symbols, only: to_symbol
   use moist_data_radii_legacy, only: get_radius, get_radius_func, rad_type
   use moist_radius_type, only: radius_type
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   implicit none
   private

   public :: static_radius_type
   public :: new_cpcm_radii
   public :: new_smd_radii
   public :: new_d3_radii
   public :: new_cosmo_radii
   public :: new_bondi_radii
   public :: new_rahm_radii
   public :: new_gauss_radii

   !> Static, table-based radii model.
   type, extends(radius_type) :: static_radius_type
      !> Selected legacy radius-model tag.
      integer :: model_tag = rad_type%cpcm
      !> Atomic numbers for the current structure.
      integer, allocatable :: atomic_numbers(:)
   contains
      !> Refresh cached radii values and derivatives for a structure.
      procedure :: update => update_static_radii
      !> Print unique elemental radii for the current molecule.
      procedure :: print => print_static_radii
   end type static_radius_type

contains

   !> Initialize a static radii object with a legacy model tag.
   !> @param[out] self       static radii model
   !> @param[in]  model_tag  legacy radii selector tag
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_static_radii(self, model_tag, verbosity)
      !> Static radii model
      type(static_radius_type), intent(out) :: self
      !> Legacy radii selector tag
      integer, intent(in) :: model_tag
      !> Optional print level
      integer, intent(in), optional :: verbosity

      self%nat = 0
      self%model_tag = model_tag
      if (present(verbosity)) then
         self%verbosity = verbosity
      else
         self%verbosity = 0
      end if
   end subroutine new_static_radii

   !> Constructor for CPCM static radii.
   !> @param[out] self       static radii model
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_cpcm_radii(self, verbosity)
      !> Static radii model
      type(static_radius_type), intent(out) :: self
      !> Optional print level
      integer, intent(in), optional :: verbosity

      call new_static_radii(self, rad_type%cpcm, verbosity)
   end subroutine new_cpcm_radii

   !> Constructor for SMD static radii.
   !> @param[out] self       static radii model
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_smd_radii(self, verbosity)
      !> Static radii model
      type(static_radius_type), intent(out) :: self
      !> Optional print level
      integer, intent(in), optional :: verbosity

      call new_static_radii(self, rad_type%smd, verbosity)
   end subroutine new_smd_radii

   !> Constructor for D3 static radii.
   !> @param[out] self       static radii model
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_d3_radii(self, verbosity)
      !> Static radii model
      type(static_radius_type), intent(out) :: self
      !> Optional print level
      integer, intent(in), optional :: verbosity

      call new_static_radii(self, rad_type%d3, verbosity)
   end subroutine new_d3_radii

   !> Constructor for COSMO static radii.
   !> @param[out] self       static radii model
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_cosmo_radii(self, verbosity)
      !> Static radii model
      type(static_radius_type), intent(out) :: self
      !> Optional print level
      integer, intent(in), optional :: verbosity

      call new_static_radii(self, rad_type%cosmo, verbosity)
   end subroutine new_cosmo_radii

   !> Constructor for Bondi static radii.
   !> @param[out] self       static radii model
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_bondi_radii(self, verbosity)
      !> Static radii model
      type(static_radius_type), intent(out) :: self
      !> Optional print level
      integer, intent(in), optional :: verbosity

      call new_static_radii(self, rad_type%bondi, verbosity)
   end subroutine new_bondi_radii

   !> Constructor for Rahm (2016) atomic radii.
   !> @param[out] self       static radii model
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_rahm_radii(self, verbosity)
      !> Static radii model
      type(static_radius_type), intent(out) :: self
      !> Optional print level
      integer, intent(in), optional :: verbosity

      call new_static_radii(self, rad_type%rahm, verbosity)
   end subroutine new_rahm_radii

   !> Constructor for Gaussian charge scheme (Bondi-based, scaled) radii.
   !> @param[out] self       static radii model
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_gauss_radii(self, verbosity)
      !> Static radii model
      type(static_radius_type), intent(out) :: self
      !> Optional print level
      integer, intent(in), optional :: verbosity

      call new_static_radii(self, rad_type%gauss, verbosity)
   end subroutine new_gauss_radii

   !> Print unique elemental radii present in the current molecule.
   !> @param[in] self  static radii model
   !> @param[in] unit  optional output unit
   subroutine print_static_radii(self, unit)
      !> Static radii model
      class(static_radius_type), intent(in) :: self
      !> Optional output unit
      integer, intent(in), optional :: unit

      !> Atom index
      integer :: iat
      !> Atomic number loop index
      integer :: iz
      !> Maximum atomic number in current molecule
      integer :: zmax
      !> Present-element mask
      logical, allocatable :: has_element(:)
      !> Pretty table printer
      type(prettylistprinter) :: plp
      !> Radius in bohr
      real(wp) :: r_bohr
      !> Radius in angstrom
      real(wp) :: r_ang
      !> Element symbol
      character(len=4) :: sym

      zmax = maxval(self%atomic_numbers)
      allocate (has_element(zmax), source=.false.)
      do iat = 1, size(self%atomic_numbers)
         if (self%atomic_numbers(iat) >= 1) then
            has_element(self%atomic_numbers(iat)) = .true.
         end if
      end do

      plp = new_prettylistprinter([6, 6, 10, 10], &
                                  [character(len=9) :: "Num", "Sym", "R (A)", "R (bohr)"], unit=unit)
      call plp%header("RADII")
      call plp%blank()
      call plp%print_header()
      call plp%separator()

      do iz = 1, zmax
         if (.not. has_element(iz)) cycle

         r_bohr = get_radius_func(iz, self%model_tag)
         if (r_bohr <= 0.0_wp) cycle
         r_ang = r_bohr*autoaa
         sym = to_symbol(iz)

         call plp%begin_row()
         call plp%add(iz)
         call plp%add(trim(sym))
         call plp%add(r_ang, fmt='f10.4')
         call plp%add(r_bohr, fmt='f10.4')
         call plp%end_row()
      end do
   end subroutine print_static_radii

   !> Update cached static radii and derivatives for a molecular structure.
   !> @param[inout] self   static radii model
   !> @param[in]    mol    molecular structure
   !> @param[out]   error  error handle on invalid input/lookup
   subroutine update_static_radii(self, mol, error)
      !> Static radii model
      class(static_radius_type), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Atom index
      integer :: iat
      !> Error message buffer
      character(len=128) :: msg

      self%nat = mol%nat
      if (self%nat < 1) then
         write (msg, '(a,i0)') "Invalid atom count in structure: ", self%nat
         call fatal_error(error, trim(msg))
         return
      end if

      if (allocated(self%atomic_numbers)) deallocate (self%atomic_numbers)
      if (allocated(self%f0)) deallocate (self%f0)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      allocate (self%atomic_numbers(self%nat), self%f0(self%nat), self%f1_rA(3, self%nat, self%nat))

      do iat = 1, self%nat
         self%atomic_numbers(iat) = mol%num(mol%id(iat))
         call get_radius(self%atomic_numbers(iat), self%model_tag, self%f0(iat), error)
         if (allocated(error)) then
            deallocate (self%atomic_numbers, self%f0, self%f1_rA)
            self%nat = 0
            return
         end if
      end do

      self%f1_rA = 0.0_wp
   end subroutine update_static_radii

end module moist_radii_static
