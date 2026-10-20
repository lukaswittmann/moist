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

   !> Bulk compressibility denominator 1 - rho C, dimensionless
   type, extends(thermodynamic_quantity_type) :: compressibility_denominator_type
   contains
      procedure :: evaluate => denominator_evaluate
      procedure :: accumulate_adjoint => denominator_accumulate
   end type compressibility_denominator_type

   !> Bulk isothermal compressibility, bohr^3/Hartree
   type, extends(thermodynamic_quantity_type) :: compressibility_item_type
   contains
      procedure :: evaluate => compressibility_evaluate
      procedure :: accumulate_adjoint => compressibility_accumulate
   end type compressibility_item_type

   !> Bulk 3D-RISM pressure, Hartree/bohr^3
   type, extends(thermodynamic_quantity_type) :: pressure_item_type
      !> Sites per solvent molecule; default 1
      integer :: nsite = 1
   contains
      procedure :: evaluate => pressure_evaluate
      procedure :: accumulate_adjoint => pressure_accumulate
   end type pressure_item_type

   !> Molecular ideal-gas pressure rho/beta, Hartree/bohr^3
   type, extends(thermodynamic_quantity_type) :: ideal_pressure_item_type
   contains
      procedure :: evaluate => ideal_pressure_evaluate
      procedure :: accumulate_adjoint => ideal_pressure_accumulate
   end type ideal_pressure_item_type

   !> Partial molar volume kappa/beta (1 - D), bohr^3
   type, extends(thermodynamic_quantity_type) :: partial_molar_volume_type
   contains
      procedure :: evaluate => volume_evaluate
      procedure :: accumulate_adjoint => volume_accumulate
   end type partial_molar_volume_type

   !> PC energy -P Vbar, Hartree; diagnostic until selected
   type, extends(thermodynamic_quantity_type) :: pc_item_type
   contains
      procedure :: evaluate => pc_evaluate
      procedure :: accumulate_adjoint => pc_accumulate
   end type pc_item_type

   !> PC+ energy -(P - P_ideal) Vbar, Hartree; diagnostic until selected
   type, extends(thermodynamic_quantity_type) :: pc_plus_item_type
   contains
      procedure :: evaluate => pc_plus_evaluate
      procedure :: accumulate_adjoint => pc_plus_accumulate
   end type pc_plus_item_type

   !> Leading net-charge finite-size term at fixed charge, box and bulk solvent
   type, extends(thermodynamic_quantity_type) :: hunenberger_correction_type
      private
      !> Total charge (e), screening S = 1 - 1/epsilon and box mean of 1/r (bohr^-1)
      real(wp) :: charge = 0.0_wp, screening = 0.0_wp, coulomb_mean = 0.0_wp
   contains
      procedure :: evaluate => hunenberger_evaluate
      procedure :: accumulate_adjoint => hunenberger_accumulate
   end type hunenberger_correction_type

contains

   !> Set up the leading net-charge finite-size correction
   !>
   !> - Fixes charge, screening and box lengths; evaluates the Coulomb box average once
   !> - Item hunenberger in Hartree, exclusive group electrostatic_finite_size, a correction
   !> - Error for non-finite input, screening outside [0, 1] or a nonpositive box length
   !>
   !> @param[out]   self       correction item to set up
   !> @param[in]    charge     total charge, e
   !> @param[in]    screening  screening S = 1 - 1/epsilon, dimensionless, in [0, 1]
   !> @param[in]    lengths    box lengths (3), bohr, positive
   !> @param[out]   error      error handler, set for invalid input
   subroutine new_hunenberger_correction(self, charge, screening, lengths, error)
      !> Correction item to set up
      type(hunenberger_correction_type), intent(out) :: self
      !> Charge (e), screening S in [0, 1], box lengths (3) in bohr
      real(wp), intent(in) :: charge, screening, lengths(3)
      !> Error handler, set for invalid input
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

   !> Evaluate -Q^2 S <1/r>_box / 2 in Hartree
   !>
   !> @param[inout] self    finite-size correction
   !> @param[in]    values  dependency values, must be empty
   !> @param[out]   value   finite-size energy, Hartree
   !> @param[out]   error   error handler, set if any dependency is passed
   subroutine hunenberger_evaluate(self, values, value, error)
      !> Finite-size correction
      class(hunenberger_correction_type), intent(inout) :: self
      !> Dependency values, must be empty
      real(wp), intent(in) :: values(:)
      !> Finite-size energy, Hartree
      real(wp), intent(out) :: value
      !> Error handler, set if any dependency is passed
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 0) then
         call fatal_error(error, "Hunenberger correction has no variable dependencies")
         return
      end if
      value = electrostatic_finite_size_correction(self%charge, self%screening, self%coulomb_mean)
   end subroutine hunenberger_evaluate

   !> No-op; the fixed charge, lattice and bulk response have no host adjoints
   !>
   !> @param[in]    self      finite-size correction
   !> @param[in]    values    dependency values, empty, unused
   !> @param[in]    value     forward value, unused
   !> @param[in]    seed      dE/dvalue, unused
   !> @param[inout] adjoints  dependency sensitivities, empty, unused
   pure subroutine hunenberger_accumulate(self, values, value, seed, adjoints)
      !> Finite-size correction
      class(hunenberger_correction_type), intent(in) :: self
      !> Dependency values, forward value and seed, all unused
      real(wp), intent(in) :: values(:), value, seed
      !> Dependency sensitivities, empty, unused
      real(wp), intent(inout) :: adjoints(:)
   end subroutine hunenberger_accumulate

   !> Build the bulk scalar graph
   !>
   !> - Inputs rho, beta, bulk_moment; items denominator, compressibility, pressure, pressure_ideal
   !> - Grid reductions supply C and the uncorrected mu separately
   !>
   !> @param[out]   self    item list, set up with the bulk inputs and items
   !> @param[in]    nsite   sites per solvent molecule
   !> @param[in]    rho     molecular number density, bohr^-3
   !> @param[in]    beta    inverse thermal energy, Hartree^-1
   !> @param[in]    moment  bulk C = sum_ab int c_ab d^3r, bohr^3
   !> @param[out]   error   error handler, allocated on failure
   subroutine new_bulk_thermo_list(self, nsite, rho, beta, moment, error)
      !> Bulk thermodynamic item list
      type(thermo_list_type), intent(out) :: self
      !> Sites per solvent molecule
      integer, intent(in) :: nsite
      !> Number density (bohr^-3), beta (Hartree^-1) and bulk C (bohr^3)
      real(wp), intent(in) :: rho, beta, moment
      !> Error handler, allocated on failure
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

   !> Build the UV scalar graph at fixed bulk state
   !>
   !> - Inputs free_energy (base energy), uv_moment, pressure, pressure_ideal, compressibility, beta
   !> - Items partial_molar_volume, pc, pc+; both PC terms are diagnostic until selected
   !>
   !> @param[out]   self         item list, set up with the UV inputs and items
   !> @param[in]    free_energy  uncorrected free energy, Hartree
   !> @param[in]    moment       UV D = sum_ia rho_a int c_ia d^3r, dimensionless
   !> @param[in]    pressure     bulk pressure, Hartree/bohr^3
   !> @param[in]    ideal        molecular ideal-gas pressure, Hartree/bohr^3
   !> @param[in]    kappa        bulk isothermal compressibility, bohr^3/Hartree
   !> @param[in]    beta         inverse thermal energy, Hartree^-1
   !> @param[out]   error        error handler, allocated on failure
   !> @param[in]    corrections  names of corrections to select (n_sel), "none" entries skipped; default none selected
   subroutine new_uv_thermo_list(self, free_energy, moment, pressure, ideal, kappa, beta, error, corrections)
      !> UV thermodynamic item list
      type(thermo_list_type), intent(out) :: self
      !> Free energy (Hartree), UV D (dimensionless), pressures (Hartree/bohr^3), kappa (bohr^3/Hartree), beta (Hartree^-1)
      real(wp), intent(in) :: free_energy, moment, pressure, ideal, kappa, beta
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      !> Names of corrections to select (n_sel), "none" entries skipped
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

   !> Evaluate 1 - rho C, dimensionless
   !>
   !> @param[inout] self    denominator item
   !> @param[in]    values  dependency values (2): rho [bohr^-3], bulk_moment [bohr^3]
   !> @param[out]   value   denominator 1 - rho C, dimensionless
   !> @param[out]   error   error handler, set for an invalid dependency count or a vanishing denominator
   subroutine denominator_evaluate(self, values, value, error)
      !> Denominator item
      class(compressibility_denominator_type), intent(inout) :: self
      !> Dependency values (2): rho [bohr^-3], bulk_moment [bohr^3]
      real(wp), intent(in) :: values(:)
      !> Denominator 1 - rho C, dimensionless
      real(wp), intent(out) :: value
      !> Error handler, set for an invalid dependency count or a vanishing denominator
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

   !> Add reverse sensitivities to rho and bulk_moment
   !>
   !> @param[in]    self      denominator item
   !> @param[in]    values    dependency values (2): rho, bulk_moment
   !> @param[in]    value     forward value of the item
   !> @param[in]    seed      dE/dvalue, Hartree per item unit
   !> @param[inout] adjoints  sensitivities (2) to rho and bulk_moment, accumulated in place
   pure subroutine denominator_accumulate(self, values, value, seed, adjoints)
      !> Denominator item
      class(compressibility_denominator_type), intent(in) :: self
      !> Dependency values, item value and seed dE/dvalue [Hartree per item unit]
      real(wp), intent(in) :: values(:), value, seed
      !> Sensitivities (2) to rho and bulk_moment, accumulated in place
      real(wp), intent(inout) :: adjoints(:)
      call denominator_adjoint(values(1), values(2), seed, adjoints)
   end subroutine denominator_accumulate

   !> Evaluate beta/[rho (1 - rho C)], bohr^3/Hartree
   !>
   !> @param[inout] self    compressibility item
   !> @param[in]    values  dependency values (3): rho [bohr^-3], beta [Hartree^-1], denominator
   !> @param[out]   value   isothermal compressibility, bohr^3/Hartree
   !> @param[out]   error   error handler, set if the dependency count is invalid
   subroutine compressibility_evaluate(self, values, value, error)
      !> Compressibility item
      class(compressibility_item_type), intent(inout) :: self
      !> Dependency values (3): rho [bohr^-3], beta [Hartree^-1], denominator
      real(wp), intent(in) :: values(:)
      !> Isothermal compressibility, bohr^3/Hartree
      real(wp), intent(out) :: value
      !> Error handler, set if the dependency count is invalid
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 3) then
         call fatal_error(error, "Invalid dependency count for the compressibility thermodynamic item")
         return
      end if
      value = bulk_compressibility(values(1), values(2), values(3))
   end subroutine compressibility_evaluate

   !> Add reverse sensitivities to rho, beta and denominator
   !>
   !> @param[in]    self      compressibility item
   !> @param[in]    values    dependency values (3): rho, beta, denominator
   !> @param[in]    value     forward value of the item
   !> @param[in]    seed      dE/dvalue, Hartree per item unit
   !> @param[inout] adjoints  sensitivities (3) to rho, beta and denominator, accumulated in place
   pure subroutine compressibility_accumulate(self, values, value, seed, adjoints)
      !> Compressibility item
      class(compressibility_item_type), intent(in) :: self
      !> Dependency values, item value and seed dE/dvalue [Hartree per item unit]
      real(wp), intent(in) :: values(:), value, seed
      !> Sensitivities (3) to rho, beta and denominator, accumulated in place
      real(wp), intent(inout) :: adjoints(:)
      call compressibility_adjoint(values(1), values(2), values(3), value, seed, adjoints)
   end subroutine compressibility_accumulate

   !> Evaluate the 3D-RISM pressure, Hartree/bohr^3
   !>
   !> @param[inout] self    pressure item
   !> @param[in]    values  dependency values (3): rho [bohr^-3], beta [Hartree^-1], bulk_moment [bohr^3]
   !> @param[out]   value   bulk pressure, Hartree/bohr^3
   !> @param[out]   error   error handler, set if the dependency count is invalid
   subroutine pressure_evaluate(self, values, value, error)
      !> Pressure item
      class(pressure_item_type), intent(inout) :: self
      !> Dependency values (3): rho [bohr^-3], beta [Hartree^-1], bulk_moment [bohr^3]
      real(wp), intent(in) :: values(:)
      !> Bulk pressure, Hartree/bohr^3
      real(wp), intent(out) :: value
      !> Error handler, set if the dependency count is invalid
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 3) then
         call fatal_error(error, "Invalid dependency count for the pressure thermodynamic item")
         return
      end if
      value = bulk_pressure(self%nsite, values(1), values(2), values(3))
   end subroutine pressure_evaluate

   !> Add reverse sensitivities to rho, beta and bulk_moment
   !>
   !> @param[in]    self      pressure item
   !> @param[in]    values    dependency values (3): rho, beta, bulk_moment
   !> @param[in]    value     forward value of the item
   !> @param[in]    seed      dE/dvalue, Hartree per item unit
   !> @param[inout] adjoints  sensitivities (3) to rho, beta and bulk_moment, accumulated in place
   pure subroutine pressure_accumulate(self, values, value, seed, adjoints)
      !> Pressure item
      class(pressure_item_type), intent(in) :: self
      !> Dependency values, item value and seed dE/dvalue [Hartree per item unit]
      real(wp), intent(in) :: values(:), value, seed
      !> Sensitivities (3) to rho, beta and bulk_moment, accumulated in place
      real(wp), intent(inout) :: adjoints(:)
      call pressure_adjoint(self%nsite, values(1), values(2), values(3), value, seed, adjoints)
   end subroutine pressure_accumulate

   !> Evaluate the molecular ideal-gas pressure rho/beta, Hartree/bohr^3
   !>
   !> @param[inout] self    ideal pressure item
   !> @param[in]    values  dependency values (2): rho [bohr^-3], beta [Hartree^-1]
   !> @param[out]   value   ideal-gas pressure, Hartree/bohr^3
   !> @param[out]   error   error handler, set if the dependency count is invalid
   subroutine ideal_pressure_evaluate(self, values, value, error)
      !> Ideal pressure item
      class(ideal_pressure_item_type), intent(inout) :: self
      !> Dependency values (2): rho [bohr^-3], beta [Hartree^-1]
      real(wp), intent(in) :: values(:)
      !> Ideal-gas pressure, Hartree/bohr^3
      real(wp), intent(out) :: value
      !> Error handler, set if the dependency count is invalid
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 2) then
         call fatal_error(error, "Invalid dependency count for the ideal_pressure thermodynamic item")
         return
      end if
      value = bulk_ideal_pressure(values(1), values(2))
   end subroutine ideal_pressure_evaluate

   !> Add reverse sensitivities to rho and beta
   !>
   !> @param[in]    self      ideal pressure item
   !> @param[in]    values    dependency values (2): rho, beta
   !> @param[in]    value     forward value of the item
   !> @param[in]    seed      dE/dvalue, Hartree per item unit
   !> @param[inout] adjoints  sensitivities (2) to rho and beta, accumulated in place
   pure subroutine ideal_pressure_accumulate(self, values, value, seed, adjoints)
      !> Ideal pressure item
      class(ideal_pressure_item_type), intent(in) :: self
      !> Dependency values, item value and seed dE/dvalue [Hartree per item unit]
      real(wp), intent(in) :: values(:), value, seed
      !> Sensitivities (2) to rho and beta, accumulated in place
      real(wp), intent(inout) :: adjoints(:)
      call ideal_pressure_adjoint(values(2), value, seed, adjoints)
   end subroutine ideal_pressure_accumulate

   !> Evaluate the partial molar volume, bohr^3
   !>
   !> @param[inout] self    partial molar volume item
   !> @param[in]    values  dependency values (3): uv_moment [1], compressibility [bohr^3/Hartree], beta [Hartree^-1]
   !> @param[out]   value   partial molar volume, bohr^3
   !> @param[out]   error   error handler, set if the dependency count is invalid
   subroutine volume_evaluate(self, values, value, error)
      !> Partial molar volume item
      class(partial_molar_volume_type), intent(inout) :: self
      !> Dependency values (3): uv_moment [1], compressibility [bohr^3/Hartree], beta [Hartree^-1]
      real(wp), intent(in) :: values(:)
      !> Partial molar volume, bohr^3
      real(wp), intent(out) :: value
      !> Error handler, set if the dependency count is invalid
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 3) then
         call fatal_error(error, "Invalid dependency count for the volume thermodynamic item")
         return
      end if
      value = partial_molar_volume(values(1), values(2), values(3))
   end subroutine volume_evaluate

   !> Add reverse sensitivities to uv_moment, compressibility and beta
   !>
   !> @param[in]    self      partial molar volume item
   !> @param[in]    values    dependency values (3): uv_moment, compressibility, beta
   !> @param[in]    value     forward value of the item
   !> @param[in]    seed      dE/dvalue, Hartree per item unit
   !> @param[inout] adjoints  sensitivities (3) to uv_moment, compressibility and beta, accumulated in place
   pure subroutine volume_accumulate(self, values, value, seed, adjoints)
      !> Partial molar volume item
      class(partial_molar_volume_type), intent(in) :: self
      !> Dependency values, item value and seed dE/dvalue [Hartree per item unit]
      real(wp), intent(in) :: values(:), value, seed
      !> Sensitivities (3) to uv_moment, compressibility and beta, accumulated in place
      real(wp), intent(inout) :: adjoints(:)
      call volume_adjoint(values(1), values(2), values(3), value, seed, adjoints)
   end subroutine volume_accumulate

   !> Evaluate the PC energy -P Vbar, Hartree
   !>
   !> @param[inout] self    PC item
   !> @param[in]    values  dependency values (2): pressure [Hartree/bohr^3], partial_molar_volume [bohr^3]
   !> @param[out]   value   PC energy, Hartree
   !> @param[out]   error   error handler, set if the dependency count is invalid
   subroutine pc_evaluate(self, values, value, error)
      !> PC item
      class(pc_item_type), intent(inout) :: self
      !> Dependency values (2): pressure [Hartree/bohr^3], partial_molar_volume [bohr^3]
      real(wp), intent(in) :: values(:)
      !> PC energy, Hartree
      real(wp), intent(out) :: value
      !> Error handler, set if the dependency count is invalid
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 2) then
         call fatal_error(error, "Invalid dependency count for the pc thermodynamic item")
         return
      end if
      value = pc_correction(values(1), values(2))
   end subroutine pc_evaluate

   !> Add reverse sensitivities to pressure and partial_molar_volume
   !>
   !> @param[in]    self      PC item
   !> @param[in]    values    dependency values (2): pressure, partial_molar_volume
   !> @param[in]    value     forward value of the item
   !> @param[in]    seed      dE/dvalue, Hartree per item unit
   !> @param[inout] adjoints  sensitivities (2) to pressure and partial_molar_volume, accumulated in place
   pure subroutine pc_accumulate(self, values, value, seed, adjoints)
      !> PC item
      class(pc_item_type), intent(in) :: self
      !> Dependency values, item value and seed dE/dvalue [Hartree per item unit]
      real(wp), intent(in) :: values(:), value, seed
      !> Sensitivities (2) to pressure and partial_molar_volume, accumulated in place
      real(wp), intent(inout) :: adjoints(:)
      call pc_adjoint(values(1), values(2), seed, adjoints)
   end subroutine pc_accumulate

   !> Evaluate the PC+ energy -(P - P_ideal) Vbar, Hartree
   !>
   !> @param[inout] self    PC+ item
   !> @param[in]    values  dependency values (3): pressure, pressure_ideal [Hartree/bohr^3], partial_molar_volume [bohr^3]
   !> @param[out]   value   PC+ energy, Hartree
   !> @param[out]   error   error handler, set if the dependency count is invalid
   subroutine pc_plus_evaluate(self, values, value, error)
      !> PC+ item
      class(pc_plus_item_type), intent(inout) :: self
      !> Dependency values (3): pressure, pressure_ideal [Hartree/bohr^3], partial_molar_volume [bohr^3]
      real(wp), intent(in) :: values(:)
      !> PC+ energy, Hartree
      real(wp), intent(out) :: value
      !> Error handler, set if the dependency count is invalid
      type(error_type), allocatable, intent(out) :: error
      if (size(values) /= 3) then
         call fatal_error(error, "Invalid dependency count for the pc_plus thermodynamic item")
         return
      end if
      value = pc_plus_correction(values(1), values(2), values(3))
   end subroutine pc_plus_evaluate

   !> Add reverse sensitivities to pressure, pressure_ideal and partial_molar_volume
   !>
   !> @param[in]    self      PC+ item
   !> @param[in]    values    dependency values (3): pressure, pressure_ideal, partial_molar_volume
   !> @param[in]    value     forward value of the item
   !> @param[in]    seed      dE/dvalue, Hartree per item unit
   !> @param[inout] adjoints  sensitivities (3) to pressure, pressure_ideal and partial_molar_volume, accumulated in place
   pure subroutine pc_plus_accumulate(self, values, value, seed, adjoints)
      !> PC+ item
      class(pc_plus_item_type), intent(in) :: self
      !> Dependency values, item value and seed dE/dvalue [Hartree per item unit]
      real(wp), intent(in) :: values(:), value, seed
      !> Sensitivities (3) to pressure, pressure_ideal and partial_molar_volume, accumulated in place
      real(wp), intent(inout) :: adjoints(:)
      call pc_plus_adjoint(values(1), values(2), values(3), seed, adjoints)
   end subroutine pc_plus_accumulate

end module moist_model_moz_thermodynamics_quantities
