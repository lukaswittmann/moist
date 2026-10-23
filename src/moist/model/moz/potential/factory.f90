!> RISM potential settings and assembly of parameter-only native potential terms
module moist_model_moz_potential_factory
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_parameters, only: moist_model_parameters_type, print_hidden
   use moist_model_moz_potential, only: moz_potential_type, lj_uff, lj_dreiding, lj_tm, lj_gaff, lj_oplsaa, &
      & lj_spce, lj_cspce, coulomb_spce, coulomb_cspce, custom_lj_type, new_custom_lj, &
      & fixed_charges_type, new_fixed_charges, &
      & lj_mixing_lorentz_berthelot, lj_mixing_geometric
   use moist_model_moz_potential_electrostatic_solvent, only: solvent_charges_type
   implicit none(type, external)
   private
   public :: rism_potential_parameters_type, configure_potential

   !> Choices for assembling native, host-independent potential terms
   type, extends(moist_model_parameters_type) :: rism_potential_parameters_type
      !> Whether the settings describe the bulk solvent side
      logical, private :: solvent_side = .false.
      !> Names of primary LJ, fallback LJ (optional), charge model and charge environment
      character(len=16) :: lj = "uff", lj_fallback = "", charges = "fixed", environment = "gas"
      !> LJ mixing rule, lorentz-berthelot (default) or geometric; solvent side only
      character(len=24) :: mixing = "lorentz-berthelot"
      !> Whether the Ng split is applied, default true
      logical :: ng_split = .true.
      !> Ng split parameter, 1/bohr, positive, default 1
      real(wp) :: alpha = 1.0_wp
      !> Optional explicit per-atom charges, e, LJ sigma, bohr, and LJ epsilon, Hartree (natom)
      real(wp), allocatable :: q(:), sigma(:), epsilon(:)
   contains
      procedure :: init_defaults => potential_defaults
      procedure :: register_entries => potential_entries
      procedure :: validate => potential_validate
      procedure :: set_solvent_side => potential_set_solvent_side
   end type rism_potential_parameters_type

contains

   !> Assemble parameter-only potential terms; the native owner updates them
   !>
   !> @param[inout] potential  empty potential, receives the parameter-only terms
   !> @param[in]    params     potential settings
   !> @param[in]    natom      atom count of this side
   !> @param[out]   error      array mismatch or native term/settings failure
   !> @param[in]    mixing     optional common solvent mixing rule for UV tables
   subroutine configure_potential(potential, params, natom, error, mixing)
      !> Potential being assembled
      type(moz_potential_type), intent(inout) :: potential
      !> Potential settings
      type(rism_potential_parameters_type), intent(in) :: params
      !> Atom count of this side
      integer, intent(in) :: natom
      !> Assembly failure
      type(error_type), allocatable, intent(out) :: error
      !> Common solvent mixing rule for UV tables, integer rule constant
      integer, intent(in), optional :: mixing
      type(fixed_charges_type) :: charges
      type(solvent_charges_type) :: table_charges
      real(wp), allocatable :: zero(:)
      integer :: rule

      call potential%set_ng_split(params%ng_split, params%alpha, error)
      if (allocated(error)) return
      rule = lj_mixing_lorentz_berthelot
      if (params%mixing == "geometric") rule = lj_mixing_geometric
      if (present(mixing)) rule = mixing
      call potential%set_mixing(rule, error)
      if (.not. allocated(error)) call add_lj(potential, params%lj, params, natom, error)
      if (.not. allocated(error) .and. len_trim(params%lj_fallback) > 0) then
         call add_lj(potential, params%lj_fallback, params, natom, error)
      end if
      if (allocated(error)) return
      select case (params%charges)
      case ("fixed")
         if (size(params%q) /= natom) then
            call fatal_error(error, "Potential q must contain one charge per atom, in structure order")
            return
         end if
         call new_fixed_charges(charges, params%q, error)
         if (.not. allocated(error)) call potential%add(charges, error)
      case ("none")
         allocate(zero(natom), source=0.0_wp)
         call new_fixed_charges(charges, zero, error)
         if (.not. allocated(error)) call potential%add(charges, error)
      case ("spce")
         call potential%add(coulomb_spce, error)
      case ("cspce")
         call potential%add(coulomb_cspce, error)
      case ("hirshfeld", "resp", "mbis", "chelpg")
         table_charges%model = params%charges
         table_charges%environment = params%environment
         call potential%add(table_charges, error)
      case default
         call fatal_error(error, "Unknown charge factory: "//trim(params%charges))
      end select
   end subroutine configure_potential

   !> Add one native LJ implementation to the potential's owned term list
   !>
   !> @param[inout] potential  potential being assembled
   !> @param[in]    name       native LJ choice
   !> @param[in]    params     explicit arrays for custom LJ
   !> @param[in]    natom      atom count
   !> @param[out]   error      unknown choice, array mismatch or native constructor failure
   subroutine add_lj(potential, name, params, natom, error)
      !> Potential being assembled
      type(moz_potential_type), intent(inout) :: potential
      !> Native LJ choice
      character(len=*), intent(in) :: name
      !> Explicit arrays for custom LJ
      type(rism_potential_parameters_type), intent(in) :: params
      !> Atom count
      integer, intent(in) :: natom
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      type(custom_lj_type) :: custom

      select case (name)
      case ("uff")
         call potential%add(lj_uff, error)
      case ("dreiding")
         call potential%add(lj_dreiding, error)
      case ("tm")
         call potential%add(lj_tm, error)
      case ("gaff")
         call potential%add(lj_gaff, error)
      case ("oplsaa")
         call potential%add(lj_oplsaa, error)
      case ("spce")
         call potential%add(lj_spce, error)
      case ("cspce")
         call potential%add(lj_cspce, error)
      case ("custom")
         if (size(params%sigma) /= natom .or. size(params%epsilon) /= natom) then
            call fatal_error(error, "Custom LJ sigma and epsilon must contain one value per atom")
            return
         end if
         call new_custom_lj(custom, params%sigma, params%epsilon, error)
         if (.not. allocated(error)) call potential%add(custom, error)
      case default
         call fatal_error(error, "Unknown LJ factory: "//trim(name))
      end select
   end subroutine add_lj

   !> Mark settings as the bulk-solvent side before defaults are restored
   !>
   !> @param[inout] self          potential settings
   !> @param[in]    solvent_side  whether the settings describe the solvent
   subroutine potential_set_solvent_side(self, solvent_side)
      !> Potential settings
      class(rism_potential_parameters_type), intent(inout) :: self
      !> Whether the settings describe the solvent
      logical, intent(in) :: solvent_side
      self%solvent_side = solvent_side
   end subroutine potential_set_solvent_side

   !> Restore potential choices; explicit solute charges have no guessed default
   !>
   !> @param[inout] self  potential settings
   subroutine potential_defaults(self)
      !> Potential settings
      class(rism_potential_parameters_type), intent(inout) :: self
      self%lj = "uff"
      self%charges = "fixed"
      if (self%solvent_side) then
         self%lj = "cspce"
         self%charges = "cspce"
      end if
      self%lj_fallback = ""
      self%environment = "gas"
      self%mixing = "lorentz-berthelot"
      self%ng_split = .true.
      self%alpha = 1.0_wp
      if (allocated(self%q)) deallocate(self%q)
      if (allocated(self%sigma)) deallocate(self%sigma)
      if (allocated(self%epsilon)) deallocate(self%epsilon)
   end subroutine potential_defaults

   !> Register potential selection and optional per-atom arrays
   !>
   !> @param[inout] self  potential settings
   subroutine potential_entries(self)
      !> Potential settings
      class(rism_potential_parameters_type), intent(inout), target :: self
      call self%register_string("lj", self%lj, level=2)
      ! An empty fallback selects none and is left out of printed tables
      call self%register_string("lj_fallback", self%lj_fallback, &
         & level=merge(2, print_hidden, len_trim(self%lj_fallback) > 0))
      call self%register_string("charges", self%charges, level=2)
      ! Per-atom values print in the potential table of the update instead
      call self%register_alloc_real_vector("q", self%q, level=print_hidden)
      call self%register_alloc_real_vector("sigma", self%sigma, level=print_hidden)
      call self%register_alloc_real_vector("epsilon", self%epsilon, level=print_hidden)
      call self%register_logical("ng_split", self%ng_split)
      call self%register_real_scalar("alpha", self%alpha)
      if (self%solvent_side) then
         call self%register_string("mixing", self%mixing)
         select case (to_lower(trim(adjustl(self%charges))))
         case ("hirshfeld", "resp", "mbis", "chelpg")
            call self%register_string("environment", self%environment, level=2)
         case default
            ! Other charge models do not use an environment setting
            continue
         end select
      end if
   end subroutine potential_entries

   !> Check choices and array values before native potential construction
   !>
   !> @param[inout] self   potential settings
   !> @param[out]   error  unknown choice, missing data or an unused array
   subroutine potential_validate(self, error)
      !> Potential settings
      class(rism_potential_parameters_type), intent(inout) :: self
      !> Invalid setting
      type(error_type), allocatable, intent(out) :: error
      self%lj = to_lower(trim(adjustl(self%lj)))
      self%lj_fallback = to_lower(trim(adjustl(self%lj_fallback)))
      self%charges = to_lower(trim(adjustl(self%charges)))
      self%environment = to_lower(trim(adjustl(self%environment)))
      self%mixing = to_lower(trim(adjustl(self%mixing)))
      if (.not. valid_lj(self%lj) .or. &
         & (len_trim(self%lj_fallback) > 0 .and. .not. valid_lj(self%lj_fallback)) &
         & .or. self%lj_fallback == "custom") then
         call fatal_error(error, "LJ choices are uff, dreiding, tm, gaff, oplsaa, spce, cspce or custom (primary only)")
         return
      end if
      if (self%lj == "custom") then
         if (.not. allocated(self%sigma) .or. .not. allocated(self%epsilon)) then
            call fatal_error(error, "Custom LJ needs per-atom sigma (bohr) and epsilon (Hartree)")
            return
         end if
      else if (allocated(self%sigma) .or. allocated(self%epsilon)) then
         call fatal_error(error, "Explicit sigma/epsilon require lj = 'custom'")
         return
      end if
      select case (self%charges)
      case ("fixed")
         if (.not. allocated(self%q)) then
            call fatal_error(error, "Fixed charges need potential.q: one charge per atom in electrons")
            return
         end if
      case ("spce", "cspce", "none")
         continue
      case ("hirshfeld", "resp", "mbis", "chelpg")
         if (.not. self%solvent_side) then
            call fatal_error(error, "Table charge models belong to solvent.potential; supply explicit solute q")
            return
         end if
         select case (self%environment)
         case ("gas", "solvent", "conductor")
            continue
         case default
            call fatal_error(error, "Charge environment must be gas, solvent or conductor")
            return
         end select
      case default
         call fatal_error(error, "Unknown potential.charges: "//trim(self%charges))
         return
      end select
      if (self%charges /= "fixed" .and. allocated(self%q)) then
         call fatal_error(error, "Explicit q requires charges = 'fixed'")
      else if (.not. positive(self%alpha)) then
         call fatal_error(error, "Potential alpha must be positive (1/bohr)")
      else if (self%mixing /= "lorentz-berthelot" .and. self%mixing /= "geometric") then
         call fatal_error(error, "Solvent mixing must be lorentz-berthelot or geometric")
      end if
   end subroutine potential_validate

   !> Whether an LJ factory is available
   !>
   !> @param[in] name  normalized choice
   pure logical function valid_lj(name)
      !> Normalized choice
      character(len=*), intent(in) :: name
      select case (name)
      case ("uff", "dreiding", "tm", "gaff", "oplsaa", "spce", "cspce", "custom")
         valid_lj = .true.
      case default
         valid_lj = .false.
      end select
   end function valid_lj

   !> Finite positive scalar
   !>
   !> @param[in] value  scalar to check
   pure logical function positive(value)
      !> Scalar to check
      real(wp), intent(in) :: value
      positive = ieee_is_finite(value) .and. value > 0.0_wp
   end function positive
end module moist_model_moz_potential_factory
