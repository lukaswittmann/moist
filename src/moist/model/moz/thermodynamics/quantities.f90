!> Concrete thermodynamic quantities and the existing bulk/UV dependency graphs
module moist_model_moz_thermodynamics_quantities
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_moz_thermodynamics_type, only: thermodynamic_quantity_type
   use moist_model_moz_thermodynamics_list, only: thermo_list_type
   use moist_model_moz_thermodynamics_kernel, only: &
      & bulk_compressibility_denominator, denominator_adjoint, bulk_compressibility, &
      & compressibility_adjoint, bulk_pressure, pressure_adjoint, &
      & bulk_ideal_pressure, ideal_pressure_adjoint, partial_molar_volume, &
      & volume_adjoint, pc_correction, pc_adjoint, &
      & pc_plus_correction, pc_plus_adjoint, coulomb_box_mean, electrostatic_finite_size_correction
   implicit none(type, external)
   private
   public :: new_bulk_thermo_list, new_uv_thermo_list
   public :: compressibility_denominator_type
   public :: compressibility_item_type
   public :: pressure_item_type
   public :: ideal_pressure_item_type
   public :: partial_molar_volume_type
   public :: pc_item_type
   public :: pc_plus_item_type
   public :: hunenberger_correction_type, new_hunenberger_correction

   type, extends(thermodynamic_quantity_type) :: compressibility_denominator_type
   contains
      procedure :: evaluate => denominator_evaluate
      procedure :: accumulate_adjoint => denominator_accumulate
   end type compressibility_denominator_type

   type, extends(thermodynamic_quantity_type) :: compressibility_item_type
   contains
      procedure :: evaluate => compressibility_evaluate
      procedure :: accumulate_adjoint => compressibility_accumulate
   end type compressibility_item_type

   type, extends(thermodynamic_quantity_type) :: pressure_item_type
      integer :: nsite = 1
   contains
      procedure :: evaluate => pressure_evaluate
      procedure :: accumulate_adjoint => pressure_accumulate
   end type pressure_item_type

   type, extends(thermodynamic_quantity_type) :: ideal_pressure_item_type
   contains
      procedure :: evaluate => ideal_pressure_evaluate
      procedure :: accumulate_adjoint => ideal_pressure_accumulate
   end type ideal_pressure_item_type

   type, extends(thermodynamic_quantity_type) :: partial_molar_volume_type
   contains
      procedure :: evaluate => volume_evaluate
      procedure :: accumulate_adjoint => volume_accumulate
   end type partial_molar_volume_type

   type, extends(thermodynamic_quantity_type) :: pc_item_type
   contains
      procedure :: evaluate => pc_evaluate
      procedure :: accumulate_adjoint => pc_accumulate
   end type pc_item_type

   type, extends(thermodynamic_quantity_type) :: pc_plus_item_type
   contains
      procedure :: evaluate => pc_plus_evaluate
      procedure :: accumulate_adjoint => pc_plus_accumulate
   end type pc_plus_item_type

   !> Leading net-charge finite-size term at fixed charge, box and bulk solvent
   type, extends(thermodynamic_quantity_type) :: hunenberger_correction_type
      private
      real(wp) :: charge = 0.0_wp, screening = 0.0_wp, coulomb_mean = 0.0_wp
   contains
      procedure :: evaluate => hunenberger_evaluate
      procedure :: accumulate_adjoint => hunenberger_accumulate
   end type hunenberger_correction_type

contains

   !> Own the fixed parameters; evaluate the Coulomb box average once at setup
   subroutine new_hunenberger_correction(self, charge, screening, lengths, error)
      type(hunenberger_correction_type), intent(out) :: self
      real(wp), intent(in) :: charge, screening, lengths(3)
      type(error_type), allocatable, intent(out) :: error

      if (.not. ieee_is_finite(charge) .or. .not. ieee_is_finite(screening) .or. &
         & .not. all(ieee_is_finite(lengths))) then
         call fatal_error(error, "Hunenberger correction needs finite charge, screening and box lengths")
         return
      end if
      if (screening < 0.0_wp .or. screening > 1.0_wp .or. any(lengths <= 0.0_wp)) then
         call fatal_error(error, "Hunenberger correction needs screening in [0,1] and positive box lengths")
         return
      end if
      self%name = "hunenberger"
      self%unit = "Hartree"
      self%correction = .true.
      self%exclusive_group = "electrostatic_finite_size"
      self%charge = charge
      self%screening = screening
      self%coulomb_mean = coulomb_box_mean(lengths)
   end subroutine new_hunenberger_correction

   subroutine hunenberger_evaluate(self, values, value, error)
      class(hunenberger_correction_type), intent(inout) :: self
      real(wp), intent(in) :: values(:)
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 0) then
         call fatal_error(error, "Hunenberger correction has no variable dependencies")
         return
      end if
      value = electrostatic_finite_size_correction(self%charge, self%screening, self%coulomb_mean)
   end subroutine hunenberger_evaluate

   !> The fixed total charge, reciprocal lattice and bulk response have no host adjoints
   pure subroutine hunenberger_accumulate(self, values, value, seed, adjoints)
      class(hunenberger_correction_type), intent(in) :: self
      real(wp), intent(in) :: values(:), value, seed
      real(wp), intent(inout) :: adjoints(:)
   end subroutine hunenberger_accumulate

   !> Bulk scalar graph; grid reductions supply C and the uncorrected mu separately
   subroutine new_bulk_thermo_list(self, nsite, rho, beta, moment, error)
      type(thermo_list_type), intent(out) :: self
      integer, intent(in) :: nsite
      real(wp), intent(in) :: rho, beta, moment
      type(error_type), allocatable, intent(out) :: error
      type(compressibility_denominator_type) :: denominator
      type(compressibility_item_type) :: compressibility
      type(pressure_item_type) :: pressure
      type(ideal_pressure_item_type) :: ideal_pressure

      call self%add_input("rho", rho, error, unit="bohr^-3")
      if (.not. allocated(error)) call self%add_input("beta", beta, error, unit="Hartree^-1")
      if (.not. allocated(error)) call self%add_input("bulk_moment", moment, error, unit="bohr^3")
      pressure%nsite = nsite
      denominator%name = "denominator"
      denominator%unit = "1"
      denominator%dependencies = [character(len=32) :: "rho", "bulk_moment"]
      if (.not. allocated(error)) call self%add(denominator, error)
      compressibility%name = "compressibility"
      compressibility%unit = "bohr^3/Hartree"
      compressibility%dependencies = [character(len=32) :: "rho", "beta", "denominator"]
      if (.not. allocated(error)) call self%add(compressibility, error)
      pressure%name = "pressure"
      pressure%unit = "Hartree/bohr^3"
      pressure%dependencies = [character(len=32) :: "rho", "beta", "bulk_moment"]
      if (.not. allocated(error)) call self%add(pressure, error)
      ideal_pressure%name = "pressure_ideal"
      ideal_pressure%unit = "Hartree/bohr^3"
      ideal_pressure%dependencies = [character(len=32) :: "rho", "beta"]
      if (.not. allocated(error)) call self%add(ideal_pressure, error)
   end subroutine new_bulk_thermo_list

   !> UV graph at fixed bulk state; both PC terms are diagnostic until selected
   subroutine new_uv_thermo_list(self, free_energy, moment, pressure, ideal, kappa, beta, error, corrections)
      type(thermo_list_type), intent(out) :: self
      real(wp), intent(in) :: free_energy, moment, pressure, ideal, kappa, beta
      type(error_type), allocatable, intent(out) :: error
      character(len=*), optional, intent(in) :: corrections(:)
      type(partial_molar_volume_type) :: volume
      type(pc_item_type) :: pc
      type(pc_plus_item_type) :: pc_plus
      integer :: i

      call self%add_input("free_energy", free_energy, error, unit="Hartree")
      if (.not. allocated(error)) call self%set_base_energy("free_energy", error)
      if (.not. allocated(error)) call self%add_input("uv_moment", moment, error, unit="1")
      if (.not. allocated(error)) call self%add_input("pressure", pressure, error, unit="Hartree/bohr^3")
      if (.not. allocated(error)) call self%add_input("pressure_ideal", ideal, error, unit="Hartree/bohr^3")
      if (.not. allocated(error)) call self%add_input("compressibility", kappa, error, unit="bohr^3/Hartree")
      if (.not. allocated(error)) call self%add_input("beta", beta, error, unit="Hartree^-1")
      volume%name = "partial_molar_volume"
      volume%unit = "bohr^3"
      volume%dependencies = [character(len=32) :: "uv_moment", "compressibility", "beta"]
      if (.not. allocated(error)) call self%add(volume, error)
      pc%name = "pc"
      pc%unit = "Hartree"
      pc%dependencies = [character(len=32) :: "pressure", "partial_molar_volume"]
      pc%correction = .true.
      pc%exclusive_group = "pressure_correction"
      if (.not. allocated(error)) call self%add(pc, error)
      pc_plus%name = "pc+"
      pc_plus%unit = "Hartree"
      pc_plus%dependencies = [character(len=32) :: "pressure", "pressure_ideal", "partial_molar_volume"]
      pc_plus%correction = .true.
      pc_plus%exclusive_group = "pressure_correction"
      if (.not. allocated(error)) call self%add(pc_plus, error)
      if (allocated(error) .or. .not. present(corrections)) return
      do i = 1, size(corrections)
         if (to_lower(trim(adjustl(corrections(i)))) == "none") cycle
         call self%select(corrections(i), .true., error)
         if (allocated(error)) return
      end do
   end subroutine new_uv_thermo_list

   subroutine denominator_evaluate(self, values, value, error)
      class(compressibility_denominator_type), intent(inout) :: self
      real(wp), intent(in) :: values(:)
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 2) then
         call fatal_error(error, "Invalid dependency count for the denominator thermodynamic item")
         return
      end if
      value = bulk_compressibility_denominator(values(1), values(2))
      if (.not. ieee_is_finite(value) .or. abs(value) <= epsilon(1.0_wp)) then
         call fatal_error(error, "1D VV compressibility is undefined: 1 - rho C vanishes")
      end if
   end subroutine denominator_evaluate

   pure subroutine denominator_accumulate(self, values, value, seed, adjoints)
      class(compressibility_denominator_type), intent(in) :: self
      real(wp), intent(in) :: values(:), value, seed
      real(wp), intent(inout) :: adjoints(:)
      call denominator_adjoint(values(1), values(2), seed, adjoints)
   end subroutine denominator_accumulate

   subroutine compressibility_evaluate(self, values, value, error)
      class(compressibility_item_type), intent(inout) :: self
      real(wp), intent(in) :: values(:)
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 3) then
         call fatal_error(error, "Invalid dependency count for the compressibility thermodynamic item")
         return
      end if
      value = bulk_compressibility(values(1), values(2), values(3))
   end subroutine compressibility_evaluate

   pure subroutine compressibility_accumulate(self, values, value, seed, adjoints)
      class(compressibility_item_type), intent(in) :: self
      real(wp), intent(in) :: values(:), value, seed
      real(wp), intent(inout) :: adjoints(:)
      call compressibility_adjoint(values(1), values(2), values(3), value, seed, adjoints)
   end subroutine compressibility_accumulate

   subroutine pressure_evaluate(self, values, value, error)
      class(pressure_item_type), intent(inout) :: self
      real(wp), intent(in) :: values(:)
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 3) then
         call fatal_error(error, "Invalid dependency count for the pressure thermodynamic item")
         return
      end if
      value = bulk_pressure(self%nsite, values(1), values(2), values(3))
   end subroutine pressure_evaluate

   pure subroutine pressure_accumulate(self, values, value, seed, adjoints)
      class(pressure_item_type), intent(in) :: self
      real(wp), intent(in) :: values(:), value, seed
      real(wp), intent(inout) :: adjoints(:)
      call pressure_adjoint(self%nsite, values(1), values(2), values(3), value, seed, adjoints)
   end subroutine pressure_accumulate

   subroutine ideal_pressure_evaluate(self, values, value, error)
      class(ideal_pressure_item_type), intent(inout) :: self
      real(wp), intent(in) :: values(:)
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 2) then
         call fatal_error(error, "Invalid dependency count for the ideal_pressure thermodynamic item")
         return
      end if
      value = bulk_ideal_pressure(values(1), values(2))
   end subroutine ideal_pressure_evaluate

   pure subroutine ideal_pressure_accumulate(self, values, value, seed, adjoints)
      class(ideal_pressure_item_type), intent(in) :: self
      real(wp), intent(in) :: values(:), value, seed
      real(wp), intent(inout) :: adjoints(:)
      call ideal_pressure_adjoint(values(2), value, seed, adjoints)
   end subroutine ideal_pressure_accumulate

   subroutine volume_evaluate(self, values, value, error)
      class(partial_molar_volume_type), intent(inout) :: self
      real(wp), intent(in) :: values(:)
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 3) then
         call fatal_error(error, "Invalid dependency count for the volume thermodynamic item")
         return
      end if
      value = partial_molar_volume(values(1), values(2), values(3))
   end subroutine volume_evaluate

   pure subroutine volume_accumulate(self, values, value, seed, adjoints)
      class(partial_molar_volume_type), intent(in) :: self
      real(wp), intent(in) :: values(:), value, seed
      real(wp), intent(inout) :: adjoints(:)
      call volume_adjoint(values(1), values(2), values(3), value, seed, adjoints)
   end subroutine volume_accumulate

   subroutine pc_evaluate(self, values, value, error)
      class(pc_item_type), intent(inout) :: self
      real(wp), intent(in) :: values(:)
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 2) then
         call fatal_error(error, "Invalid dependency count for the pc thermodynamic item")
         return
      end if
      value = pc_correction(values(1), values(2))
   end subroutine pc_evaluate

   pure subroutine pc_accumulate(self, values, value, seed, adjoints)
      class(pc_item_type), intent(in) :: self
      real(wp), intent(in) :: values(:), value, seed
      real(wp), intent(inout) :: adjoints(:)
      call pc_adjoint(values(1), values(2), seed, adjoints)
   end subroutine pc_accumulate

   subroutine pc_plus_evaluate(self, values, value, error)
      class(pc_plus_item_type), intent(inout) :: self
      real(wp), intent(in) :: values(:)
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 3) then
         call fatal_error(error, "Invalid dependency count for the pc_plus thermodynamic item")
         return
      end if
      value = pc_plus_correction(values(1), values(2), values(3))
   end subroutine pc_plus_evaluate

   pure subroutine pc_plus_accumulate(self, values, value, seed, adjoints)
      class(pc_plus_item_type), intent(in) :: self
      real(wp), intent(in) :: values(:), value, seed
      real(wp), intent(inout) :: adjoints(:)
      call pc_plus_adjoint(values(1), values(2), values(3), seed, adjoints)
   end subroutine pc_plus_accumulate

end module moist_model_moz_thermodynamics_quantities
