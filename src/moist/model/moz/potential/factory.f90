!> RISM potential settings and assembly of parameter-only native potential terms
!>
!> - `lj` and `charges` list choices in order of precedence; each atom takes the first choice covering it
!> - Charge choices name the source of the site charges
!>   - `custom`: explicit `q`, one charge per atom
!>   - `none`: zero charges
!>   - `spce`: SPC/E water charges; other structures left to later choices
!>   - `<model>-<environment>`: solvent table charges, solvent side only;
!>     models hirshfeld, resp, mbis, chelpg; environments gas, solvent, conductor
!>   - `ec`: charges fitted to the host potential, solute side only; fit settings in `ec`
!>   - `custom`, `none` and `ec` cover every atom and stand alone
!> - Solvent-side defaults: cSPC/E for water, else GAFF2, UFF and TM LJ with RESP charges
module moist_model_moz_potential_factory
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_parameters, only: moist_model_parameters_type, print_hidden
   use moist_model_moz_potential, only: moz_potential_type, lj_uff, lj_dreiding, lj_tm, lj_gaff, lj_oplsaa, &
      & lj_spce, lj_cspce, coulomb_spce, custom_lj_type, new_custom_lj, &
      & fixed_charges_type, new_fixed_charges, new_neutral_charges, ec_charges_type, &
      & lj_mixing_lorentz_berthelot, lj_mixing_geometric
   use moist_model_moz_potential_electrostatic_solvent, only: solvent_charges_type
   implicit none(type, external)
   private
   public :: rism_potential_parameters_type, configure_potential

   !> Choices for assembling native, host-independent potential terms
   type, extends(moist_model_parameters_type) :: rism_potential_parameters_type
      !> Whether the settings describe the bulk solvent side
      logical, private :: solvent_side = .false.
      !> LJ tables in order of precedence; each atom takes the first table covering it
      character(len=16), allocatable :: lj(:)
      !> Sources of the site charges in order of precedence; each atom takes the first covering it
      character(len=24), allocatable :: charges(:)
      !> LJ mixing rule, lorentz-berthelot (default) or geometric; solvent side only
      character(len=24) :: mixing = "lorentz-berthelot"
      !> Whether the Ng split is applied, default true
      logical :: ng_split = .true.
      !> Ng split parameter, 1/bohr, positive, default 1
      real(wp) :: alpha = 1.0_wp
      !> Optional explicit per-atom charges, e, LJ sigma, bohr, and LJ epsilon, Hartree (natom)
      real(wp), allocatable :: q(:), sigma(:), epsilon(:)
      !> EC term carrying the fit settings of `charges = "ec"`; per-atom settings unused
      type(ec_charges_type) :: ec
   contains
      procedure :: init_defaults => potential_defaults
      procedure :: register_entries => potential_entries
      procedure :: validate => potential_validate
      procedure :: set_solvent_side => potential_set_solvent_side
   end type rism_potential_parameters_type

contains

   !> Assemble parameter-only potential terms; the native owner updates them
   !>
   !> - Atom counts of explicit arrays are checked by the terms at update
   !>
   !> @param[inout] potential  empty potential, receives the parameter-only terms
   !> @param[in]    params     validated potential settings
   !> @param[out]   error      native term or settings failure
   !> @param[in]    mixing     optional common solvent mixing rule for UV tables
   subroutine configure_potential(potential, params, error, mixing)
      !> Potential being assembled
      type(moz_potential_type), intent(inout) :: potential
      !> Validated potential settings
      type(rism_potential_parameters_type), intent(in) :: params
      !> Assembly failure
      type(error_type), allocatable, intent(out) :: error
      !> Common solvent mixing rule for UV tables, integer rule constant
      integer, intent(in), optional :: mixing
      integer :: rule, i

      call potential%set_ng_split(params%ng_split, params%alpha, error)
      if (allocated(error)) return
      rule = lj_mixing_lorentz_berthelot
      if (params%mixing == "geometric") rule = lj_mixing_geometric
      if (present(mixing)) rule = mixing
      call potential%set_mixing(rule, error)
      do i = 1, size(params%lj)
         if (allocated(error)) exit
         call add_lj(potential, trim(params%lj(i)), params, error)
      end do
      do i = 1, size(params%charges)
         if (allocated(error)) exit
         call add_charges(potential, trim(params%charges(i)), params, error)
      end do
   end subroutine configure_potential

   !> Add one native charge term to the potential's owned term list
   !>
   !> @param[inout] potential  potential being assembled
   !> @param[in]    name       native charge choice
   !> @param[in]    params     explicit charges and EC fit settings
   !> @param[out]   error      unknown choice or native constructor failure
   subroutine add_charges(potential, name, params, error)
      !> Potential being assembled
      type(moz_potential_type), intent(inout) :: potential
      !> Native charge choice
      character(len=*), intent(in) :: name
      !> Explicit charges and EC fit settings
      type(rism_potential_parameters_type), intent(in) :: params
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      type(fixed_charges_type) :: charges
      type(solvent_charges_type) :: table_charges
      integer :: dash

      select case (name)
      case ("custom")
         call new_fixed_charges(charges, params%q, error)
         if (.not. allocated(error)) call potential%add(charges, error)
      case ("none")
         call new_neutral_charges(charges)
         call potential%add(charges, error)
      case ("spce")
         call potential%add(coulomb_spce, error)
      case ("ec")
         call potential%add(params%ec, error)
      case default
         if (.not. is_table_charges(name)) then
            call fatal_error(error, "Unknown charge factory: "//name)
            return
         end if
         dash = index(name, "-")
         table_charges%model = name(:dash - 1)
         table_charges%environment = name(dash + 1:)
         call potential%add(table_charges, error)
      end select
   end subroutine add_charges

   !> Add one native LJ implementation to the potential's owned term list
   !>
   !> @param[inout] potential  potential being assembled
   !> @param[in]    name       native LJ choice
   !> @param[in]    params     explicit arrays for custom LJ
   !> @param[out]   error      unknown choice or native constructor failure
   subroutine add_lj(potential, name, params, error)
      !> Potential being assembled
      type(moz_potential_type), intent(inout) :: potential
      !> Native LJ choice
      character(len=*), intent(in) :: name
      !> Explicit arrays for custom LJ
      type(rism_potential_parameters_type), intent(in) :: params
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
      case ("gaff2")
         call potential%add(lj_gaff, error)
      case ("oplsaa")
         call potential%add(lj_oplsaa, error)
      case ("spce")
         call potential%add(lj_spce, error)
      case ("cspce")
         call potential%add(lj_cspce, error)
      case ("custom")
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
      type(ec_charges_type) :: ec
      self%lj = [character(len=16) :: "uff"]
      self%charges = [character(len=24) :: "custom"]
      if (self%solvent_side) then
         ! SPC/E models for water; GAFF2, UFF, TM and RESP charges for any other solvent
         self%lj = [character(len=16) :: "cspce", "gaff2", "uff", "tm"]
         self%charges = [character(len=24) :: "spce", "resp-solvent"]
      end if
      self%mixing = "lorentz-berthelot"
      self%ng_split = .true.
      self%alpha = 1.0_wp
      if (allocated(self%q)) deallocate(self%q)
      if (allocated(self%sigma)) deallocate(self%sigma)
      if (allocated(self%epsilon)) deallocate(self%epsilon)
      self%ec = ec
   end subroutine potential_defaults

   !> Register potential selection and optional per-atom arrays
   !>
   !> - The EC fit settings are read only with `charges = "ec"`; validate refuses them on the solvent side
   !>
   !> @param[inout] self  potential settings
   subroutine potential_entries(self)
      !> Potential settings
      class(rism_potential_parameters_type), intent(inout), target :: self
      call self%register_alloc_string_vector("lj", self%lj, level=2)
      call self%register_alloc_string_vector("charges", self%charges, level=2)
      ! Per-atom values print in the potential table of the update instead
      call self%register_alloc_real_vector("q", self%q, level=print_hidden)
      call self%register_alloc_real_vector("sigma", self%sigma, level=print_hidden)
      call self%register_alloc_real_vector("epsilon", self%epsilon, level=print_hidden)
      call self%register_logical("ng_split", self%ng_split)
      call self%register_real_scalar("alpha", self%alpha)
      if (self%solvent_side) call self%register_string("mixing", self%mixing)
      if (has_choice(self%charges, "ec")) then
         call self%register_real_scalar("ec.exclusion_radius", self%ec%exclusion_radius)
         call self%register_real_scalar("ec.range_scale", self%ec%range_scale)
         call self%register_real_scalar("ec.long_range_bias", self%ec%long_range_bias)
         call self%register_real_scalar("ec.restraint_strength", self%ec%restraint_strength)
         call self%register_real_scalar("ec.charge_penalty", self%ec%charge_penalty)
      end if
   end subroutine potential_entries

   !> Check choices and array values before native potential construction; normalizes the names
   !>
   !> @param[inout] self   potential settings
   !> @param[out]   error  unknown choice, missing data or an unused array
   subroutine potential_validate(self, error)
      !> Potential settings
      class(rism_potential_parameters_type), intent(inout) :: self
      !> Invalid setting
      type(error_type), allocatable, intent(out) :: error
      self%mixing = to_lower(trim(adjustl(self%mixing)))
      call validate_lj(self, error)
      if (.not. allocated(error)) call validate_charges(self, error)
      if (allocated(error)) return
      if (.not. positive(self%alpha)) then
         call fatal_error(error, "Potential alpha must be positive (1/bohr)")
      else if (self%mixing /= "lorentz-berthelot" .and. self%mixing /= "geometric") then
         call fatal_error(error, "Solvent mixing must be lorentz-berthelot or geometric")
      end if
   end subroutine potential_validate

   !> Check the LJ tables and the custom LJ arrays
   !>
   !> @param[inout] self   potential settings
   !> @param[out]   error  empty, unknown or repeated table, custom LJ beside another table or missing arrays
   subroutine validate_lj(self, error)
      !> Potential settings
      class(rism_potential_parameters_type), intent(inout) :: self
      !> Invalid setting
      type(error_type), allocatable, intent(out) :: error
      integer :: i

      call normalize_choices(self%lj, "lj", error)
      if (allocated(error)) return
      do i = 1, size(self%lj)
         if (.not. valid_lj(self%lj(i))) then
            call fatal_error(error, "Unknown LJ table: "//trim(self%lj(i))// &
               & "; choices are uff, dreiding, tm, gaff2, oplsaa, spce, cspce or custom")
            return
         end if
      end do
      if (any(self%lj == "custom") .and. size(self%lj) > 1) then
         call fatal_error(error, "Custom LJ covers every atom and must be the only LJ table")
      else if (self%lj(1) == "custom") then
         if (.not. allocated(self%sigma) .or. .not. allocated(self%epsilon)) then
            call fatal_error(error, "Custom LJ needs per-atom sigma (bohr) and epsilon (Hartree)")
         end if
      else if (allocated(self%sigma) .or. allocated(self%epsilon)) then
         call fatal_error(error, "Explicit sigma/epsilon require lj = 'custom'")
      end if
   end subroutine validate_lj

   !> Check the charge choices, their side and the explicit charges
   !>
   !> @param[inout] self   potential settings
   !> @param[out]   error  empty, unknown or repeated choice, a whole-side choice beside another,
   !>                      a choice of the other side or q without custom charges
   subroutine validate_charges(self, error)
      !> Potential settings
      class(rism_potential_parameters_type), intent(inout) :: self
      !> Invalid setting
      type(error_type), allocatable, intent(out) :: error
      integer :: i

      call normalize_choices(self%charges, "charges", error)
      if (allocated(error)) return
      do i = 1, size(self%charges)
         select case (self%charges(i))
         case ("custom", "none", "ec")
            if (size(self%charges) > 1) then
               call fatal_error(error, "Charges "//trim(self%charges(i))// &
                  & " cover every atom and must be the only charge choice")
            else if (self%charges(i) == "custom" .and. .not. allocated(self%q)) then
               call fatal_error(error, "Custom charges need potential.q: one charge per atom in electrons")
            else if (self%charges(i) == "ec" .and. self%solvent_side) then
               call fatal_error(error, "EC charges belong to the solute potential")
            end if
         case ("spce")
            continue
         case default
            if (.not. is_table_charges(self%charges(i))) then
               call fatal_error(error, "Unknown potential charges: "//trim(self%charges(i))// &
                  & "; choices are custom, none, spce, ec or <hirshfeld|resp|mbis|chelpg>-<gas|solvent|conductor>")
            else if (.not. self%solvent_side) then
               call fatal_error(error, "Table charge models belong to solvent.potential; supply explicit solute q")
            end if
         end select
         if (allocated(error)) return
      end do
      if (allocated(self%q) .and. self%charges(1) /= "custom") then
         call fatal_error(error, "Explicit q requires charges = 'custom'")
      end if
   end subroutine validate_charges

   !> Normalize a choice list and require it to be nonempty without repeats
   !>
   !> @param[inout] choices  choice list
   !> @param[in]    key      setting name for the error
   !> @param[out]   error    empty list or repeated choice
   subroutine normalize_choices(choices, key, error)
      !> Choice list
      ! allow(assumed-size-character-intent): fixed-length entries of the settings
      character(len=*), allocatable, intent(inout) :: choices(:)
      !> Setting name for the error
      character(len=*), intent(in) :: key
      !> Empty list or repeated choice
      type(error_type), allocatable, intent(out) :: error
      integer :: i

      if (.not. allocated(choices)) then
         call fatal_error(error, "Potential "//key//" needs at least one choice")
         return
      else if (size(choices) == 0) then
         call fatal_error(error, "Potential "//key//" needs at least one choice")
         return
      end if
      choices = to_lower(adjustl(choices))
      do i = 2, size(choices)
         if (any(choices(:i - 1) == choices(i))) then
            call fatal_error(error, "Potential "//key//" lists "//trim(choices(i))//" twice")
            return
         end if
      end do
   end subroutine normalize_choices

   !> Whether a choice list, not yet normalized, names a choice
   !>
   !> @param[in] choices  optional choice list
   !> @param[in] name     normalized choice
   pure function has_choice(choices, name) result(found)
      !> Optional choice list
      character(len=*), allocatable, intent(in) :: choices(:)
      !> Normalized choice
      character(len=*), intent(in) :: name
      !> Whether the list names the choice
      logical :: found
      found = .false.
      if (allocated(choices)) found = any(to_lower(adjustl(choices)) == name)
   end function has_choice

   !> Whether an LJ factory is available
   !>
   !> @param[in] name  normalized choice
   pure function valid_lj(name) result(valid)
      !> Normalized choice
      character(len=*), intent(in) :: name
      !> Whether the factory exists
      logical :: valid
      select case (name)
      case ("uff", "dreiding", "tm", "gaff2", "oplsaa", "spce", "cspce", "custom")
         valid = .true.
      case default
         valid = .false.
      end select
   end function valid_lj

   !> Whether a charge choice names solvent table charges, `<model>-<environment>`
   !>
   !> @param[in] name  normalized choice
   pure function is_table_charges(name) result(table)
      !> Normalized choice
      character(len=*), intent(in) :: name
      !> Whether the choice names table charges
      logical :: table
      integer :: dash
      table = .false.
      dash = index(name, "-")
      if (dash < 2) return
      select case (name(:dash - 1))
      case ("hirshfeld", "resp", "mbis", "chelpg")
         continue
      case default
         return
      end select
      select case (trim(name(dash + 1:)))
      case ("gas", "solvent", "conductor")
         table = .true.
      case default
         continue
      end select
   end function is_table_charges

   !> Finite positive scalar
   !>
   !> @param[in] value  scalar to check
   pure logical function positive(value)
      !> Scalar to check
      real(wp), intent(in) :: value
      positive = ieee_is_finite(value) .and. value > 0.0_wp
   end function positive
end module moist_model_moz_potential_factory
