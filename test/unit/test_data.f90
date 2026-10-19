!> Tests for the parameter tables under src/moist/data/ and the accessors
module test_data
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type, fatal_error
   use mctc_io, only: structure_type, new_structure
   use mctc_io_codata2018, only: Avogadro_constant, atomic_unit_of_mass
   use mctc_io_symbols, only: to_symbol
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

   use moist_data_en, only: get_electronegativity, en_max_elem => max_elem
   use moist_data_hardness, only: get_hardness, hardness_max_elem => max_elem
   use moist_data_mass, only: get_mass, mass_max_elem => max_elem
   use moist_data_atomicrad, only: get_atomic_rad, get_covalent_rad, &
      & arad_max_elem => max_elem
   use moist_data_radii_legacy, only: get_radius, get_radius_func, &
      & get_upper_bound, rad_type
   use moist_data_solvents, only: get_solvent_id, max_solvents, &
      & solvation_system_type, new_solvation_system, solvent_multipole_data_type, &
      & get_solvent_charges, get_solvent_multipoles

   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none(type, external)
   private

   public :: collect_data

   !> Anchor comparisons are between identical literal expressions
   real(wp), parameter :: thr = 1.0e-10_wp

   !> Tags selecting one of the element-indexed accessors
   integer, parameter :: acc_en = 1
   integer, parameter :: acc_hardness = 2
   integer, parameter :: acc_mass = 3
   integer, parameter :: acc_arad = 4
   integer, parameter :: acc_crad = 5
   integer, parameter :: n_accessors = 5

   character(len=*), parameter :: acc_name(n_accessors) = [character(len=17) :: &
      & "electronegativity", "hardness", "mass", "atomic_rad", "covalent_rad"]

   !> Every radius model tag known to moist_data_radii_legacy
   integer, parameter :: n_models = 7
   character(len=*), parameter :: model_name(n_models) = [character(len=8) :: &
      & "cpcm", "smd", "d3", "cosmo", "bondi", "rahm", "gauss"]

contains

   !> Collect all data-table tests
   subroutine collect_data(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("element_tables_complete", test_element_tables_complete), &
                  new_unittest("element_tables_physical", test_element_tables_physical), &
                  new_unittest("hardness_superheavy_are_zero", test_hardness_superheavy), &
                  new_unittest("accessor_rejects_out_of_range", test_accessor_out_of_range), &
                  new_unittest("accessor_rejects_bad_symbol", test_accessor_bad_symbol), &
                  new_unittest("accessor_symbol_case_insensitive", test_accessor_symbol_case), &
                  new_unittest("accessor_symbol_matches_number", test_accessor_symbol_consistency), &
                  new_unittest("radius_all_models_complete", test_radius_models_complete), &
                  new_unittest("radius_model_upper_bounds", test_radius_upper_bounds), &
                  new_unittest("radius_keyword_normalisation", test_radius_keyword_normalisation), &
                  new_unittest("radius_rejects_bad_keyword", test_radius_bad_keyword), &
                  new_unittest("radius_models_are_distinct", test_radius_models_distinct), &
                  new_unittest("radius_model_error_wins_over_symbol", test_radius_error_precedence), &
                  new_unittest("radius_func_reports_error", test_radius_func_reports_error), &
                  new_unittest("solvent_alias_round_trip", test_solvent_alias_round_trip), &
                  new_unittest("solvent_alias_case_and_blanks", test_solvent_alias_normalisation), &
                  new_unittest("solvent_rejects_blank_alias", test_solvent_blank_alias), &
                  new_unittest("solvent_system_constructs", test_solvent_system_constructs), &
                  new_unittest("solvent_system_validates_input", test_solvent_system_validation), &
                  new_unittest("solvent_system_all_ids", test_solvent_system_all_ids), &
                  new_unittest("solvent_charge_data", test_solvent_charge_data), &
                  new_unittest("solvent_system_charge_accessors", test_solvent_system_charge_accessors), &
                  new_unittest("solute_update_ownership", test_solute_update_ownership), &
                  new_unittest("water_higher_multipoles", test_water_higher_multipoles) &
                  ]
   end subroutine collect_data

    !* -------------------------------- Private helpers -------------------------------- *!

   !> Dispatch to one of the element-indexed accessors by tag
   subroutine lookup_num(tag, num, val, err)
      !> Accessor tag, one of the acc_* parameters
      integer, intent(in) :: tag
      !> Atomic number to look up
      integer, intent(in) :: num
      !> Table value
      real(wp), intent(out) :: val
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: err

      select case (tag)
         case (acc_en); call get_electronegativity(num, val, err)
         case (acc_hardness); call get_hardness(num, val, err)
         case (acc_mass); call get_mass(num, val, err)
         case (acc_arad); call get_atomic_rad(num, val, err)
         case (acc_crad); call get_covalent_rad(num, val, err)
         case default
            call fatal_error(err, "lookup_num: unknown accessor tag")
            return
      end select
   end subroutine lookup_num

   !> Dispatch to the symbol overload of one of the element-indexed accessors
   subroutine lookup_sym(tag, sym, val, err)
      !> Accessor tag, one of the acc_* parameters
      integer, intent(in) :: tag
      !> Element symbol to look up
      character(len=*), intent(in) :: sym
      !> Table value
      real(wp), intent(out) :: val
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: err

      select case (tag)
         case (acc_en); call get_electronegativity(sym, val, err)
         case (acc_hardness); call get_hardness(sym, val, err)
         case (acc_mass); call get_mass(sym, val, err)
         case (acc_arad); call get_atomic_rad(sym, val, err)
         case (acc_crad); call get_covalent_rad(sym, val, err)
         case default
            call fatal_error(err, "lookup_num: unknown accessor tag")
            return
      end select
   end subroutine lookup_sym

   !> Upper bound of the table behind a given accessor tag
   pure function accessor_max_elem(tag) result(upper)
      !> Accessor tag, one of the acc_* parameters
      integer, intent(in) :: tag
      !> Highest atomic number the table covers
      integer :: upper

      select case (tag)
      case (acc_en); upper = en_max_elem
      case (acc_hardness); upper = hardness_max_elem
      case (acc_mass); upper = mass_max_elem
      case default; upper = arad_max_elem
      end select
   end function accessor_max_elem

   !* -------------- Group A: structural invariants over the whole table -------------- *!

   !> Every element table must answer for every atomic number in its declared
   !> range, with a finite value and no error. A constructor that lost an entry
   !> would either trip the bound here or leave a trailing element unreachable
   subroutine test_element_tables_complete(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: tag, iz, upper
      real(wp) :: val

      do tag = 1, n_accessors
         upper = accessor_max_elem(tag)
         call check(error, upper, 118, more="unexpected range for "//trim(acc_name(tag)))
         if (allocated(error)) return

         do iz = 1, upper
            call lookup_num(tag, iz, val, err)
            if (allocated(err)) then
               call test_failed(error, trim(acc_name(tag))//" rejected a valid Z: "//trim(err%message))
               return
            end if
            if (.not. ieee_is_finite(val)) then
               call test_failed(error, trim(acc_name(tag))//" returned a non-finite value")
               return
            end if
         end do
      end do
   end subroutine test_element_tables_complete

   !> Masses, radii and electronegativities are strictly positive for every
   !> element. A shifted table tends to survive the completeness check above but
   !> not this one, because the shifted-in filler is usually zero
   subroutine test_element_tables_physical(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: tag, iz
      real(wp) :: val

      do tag = 1, n_accessors
         if (tag == acc_hardness) cycle  ! legitimately zero for Rf-Og, see below
         do iz = 1, accessor_max_elem(tag)
            call lookup_num(tag, iz, val, err)
            if (allocated(err)) then
               call test_failed(error, "unexpected error from "//trim(acc_name(tag)))
               return
            end if
            if (val <= 0.0_wp) then
               call test_failed(error, trim(acc_name(tag))//" is not positive for all elements")
               return
            end if
         end do
      end do
   end subroutine test_element_tables_physical

   !> DFT-D4 does not parametrise Rf-Og, so those hardnesses are exactly zero
   !> while every lighter element is positive. Pinning this matters because a
   !> zero used to be indistinguishable from the old out-of-range sentinel:
   !> callers must branch on the error, not on the value
   subroutine test_hardness_superheavy(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: iz, nzero
      real(wp) :: eta

      nzero = 0
      do iz = 1, hardness_max_elem
         call get_hardness(iz, eta, err)
         if (allocated(err)) then
            call test_failed(error, "hardness rejected a valid Z: "//trim(err%message))
            return
         end if
         if (eta == 0.0_wp) then
            nzero = nzero + 1
            if (iz < 104) then
               call test_failed(error, "unexpected zero hardness below Rf")
               return
            end if
         else if (eta < 0.0_wp) then
            call test_failed(error, "negative chemical hardness")
            return
         end if
      end do

      call check(error, nzero, 15, more="Rf-Og (15 elements) must be the only zero hardnesses")
   end subroutine test_hardness_superheavy

   !* ------------------------- Group B: the accessor contract ------------------------ *!

   !> Out-of-range atomic numbers must be rejected, and the output left at zero
   !> rather than carrying a sentinel the caller might mistake for data
   subroutine test_accessor_out_of_range(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: bad_z(5) = [0, -1, 119, huge(1), -huge(1)]

      type(moist_error_type), allocatable :: err
      integer :: tag, i
      real(wp) :: val

      do tag = 1, n_accessors
         do i = 1, size(bad_z)
            val = 1.0_wp
            call lookup_num(tag, bad_z(i), val, err)
            call check(error, allocated(err), &
                       more=trim(acc_name(tag))//" accepted an out-of-range atomic number")
            if (allocated(error)) return
            call check(error, val, 0.0_wp, thr=0.0_wp, &
                       more="rejected lookup left a non-zero value behind")
            if (allocated(error)) return
            deallocate (err)
         end do
      end do
   end subroutine test_accessor_out_of_range

   !> Unknown, empty and blank symbols must all be rejected. Previously these
   !> resolved to atomic number zero and fell through to a silent sentinel
   subroutine test_accessor_bad_symbol(error)
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: bad_sym(4) = [character(len=4) :: "Xx", "", "  ", "1234"]

      type(moist_error_type), allocatable :: err
      integer :: tag, i
      real(wp) :: val

      do tag = 1, n_accessors
         do i = 1, size(bad_sym)
            val = 1.0_wp
            call lookup_sym(tag, bad_sym(i), val, err)
            call check(error, allocated(err), &
                       more=trim(acc_name(tag))//" accepted an invalid element symbol")
            if (allocated(error)) return
            call check(error, val, 0.0_wp, thr=0.0_wp, &
                       more="rejected symbol lookup left a non-zero value behind")
            if (allocated(error)) return
            deallocate (err)
         end do
      end do
   end subroutine test_accessor_bad_symbol

   !> Symbol lookup is case-insensitive and tolerates padding, and the
   !> deuterium/tritium aliases resolve to hydrogen
   subroutine test_accessor_symbol_case(error)
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: hydrogen(5) = [character(len=4) :: "H", "h", " H  ", "D", "T"]

      type(moist_error_type), allocatable :: err
      integer :: i
      real(wp) :: reference, val

      call get_electronegativity(1, reference, err)
      if (allocated(err)) then
         call test_failed(error, "hydrogen lookup failed: "//trim(err%message))
         return
      end if

      do i = 1, size(hydrogen)
         call get_electronegativity(hydrogen(i), val, err)
         if (allocated(err)) then
            call test_failed(error, "symbol '"//trim(hydrogen(i))//"' rejected: "//trim(err%message))
            return
         end if
         call check(error, val, reference, thr=0.0_wp, &
                    more="symbol '"//trim(hydrogen(i))//"' did not resolve to hydrogen")
         if (allocated(error)) return
      end do

      ! Mixed case on a two-letter symbol
      call get_electronegativity("hE", val, err)
      if (allocated(err)) then
         call test_failed(error, "mixed-case symbol rejected: "//trim(err%message))
         return
      end if
      call get_electronegativity(2, reference, err)
      call check(error, val, reference, thr=0.0_wp, more="'hE' did not resolve to helium")
   end subroutine test_accessor_symbol_case

   !> The symbol and atomic-number overloads must agree for every element
   !> This walks every entry of every table through both paths without
   !> restating a single tabulated value
   subroutine test_accessor_symbol_consistency(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err_num, err_sym
      integer :: tag, iz
      real(wp) :: by_num, by_sym

      do tag = 1, n_accessors
         do iz = 1, accessor_max_elem(tag)
            call lookup_num(tag, iz, by_num, err_num)
            call lookup_sym(tag, to_symbol(iz), by_sym, err_sym)

            if (allocated(err_num) .or. allocated(err_sym)) then
               call test_failed(error, trim(acc_name(tag))//" failed on a valid element")
               return
            end if
            if (by_num /= by_sym) then
               call test_failed(error, trim(acc_name(tag))//" symbol and number overloads disagree")
               return
            end if
         end do
      end do
   end subroutine test_accessor_symbol_consistency

   !* --------------------- Group C: radius-model keyword dispatch -------------------- *!

   !> Every model must answer for every atomic number up to its own bound. The
   !> d3, cosmo, bondi and rahm tables declare literal extents decoupled from the
   !> max_elem_* constants used to guard them, so a mismatch would otherwise be
   !> an unchecked out-of-bounds read
   subroutine test_radius_models_complete(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: imodel, iz, upper, nmissing
      real(wp) :: rad

      do imodel = 1, n_models
         call get_upper_bound(imodel, upper, err)
         if (allocated(err)) then
            call test_failed(error, "model "//trim(model_name(imodel))//" has no upper bound")
            return
         end if
         nmissing = 0

         do iz = 1, upper
            call get_radius(iz, imodel, rad, err)
            if (allocated(err)) then
               ! Bondi has genuine gaps, flagged by the negative `missing`
               ! sentinel; every other model must be complete
               nmissing = nmissing + 1
               deallocate (err)
               cycle
            end if
            if (.not. ieee_is_finite(rad) .or. rad <= 0.0_wp) then
               call test_failed(error, "model "//trim(model_name(imodel))//" gave an unusable radius")
               return
            end if
         end do

         if (imodel /= rad_type%bondi .and. nmissing /= 0) then
            call test_failed(error, "model "//trim(model_name(imodel))//" has unexpected gaps")
            return
         end if
      end do
   end subroutine test_radius_models_complete

   !> Each model's documented upper bound is accepted and one past it rejected
   subroutine test_radius_upper_bounds(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: expected_upper(n_models) = [118, 118, 94, 94, 88, 96, 118]

      type(moist_error_type), allocatable :: err
      integer :: imodel, upper
      real(wp) :: rad

      do imodel = 1, n_models
         call get_upper_bound(imodel, upper, err)
         call check(error, .not. allocated(err), &
                    "model "//trim(model_name(imodel))//" was not recognised")
         if (allocated(error)) return
         call check(error, upper, expected_upper(imodel), &
                    more="unexpected upper bound for "//trim(model_name(imodel)))
         if (allocated(error)) return

         rad = 1.0_wp
         call get_radius(upper + 1, imodel, rad, err)
         call check(error, allocated(err), &
                    more="model "//trim(model_name(imodel))//" accepted Z past its bound")
         if (allocated(error)) return
         call check(error, rad, 0.0_wp, thr=0.0_wp, more="rejected radius lookup was not zeroed")
         if (allocated(error)) return
         deallocate (err)

         call get_radius(0, imodel, rad, err)
         call check(error, allocated(err), &
                    more="model "//trim(model_name(imodel))//" accepted Z = 0")
         if (allocated(error)) return
         deallocate (err)
      end do

      ! An unknown tag has no bound and must be reported as an error, not as a
      ! sentinel the caller could mistake for a real bound
      call get_upper_bound(n_models + 1, upper, err)
      call check(error, allocated(err), more="an unknown model tag was accepted")
      if (allocated(error)) return
      call check(error, upper, -1, more="rejected bound lookup left a usable value")
   end subroutine test_radius_upper_bounds

   !> Model names are matched case-insensitively after trimming and adjusting,
   !> and every name must agree with its integer tag
   subroutine test_radius_keyword_normalisation(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: imodel
      real(wp) :: by_tag, by_name

      do imodel = 1, n_models
         call get_radius(6, imodel, by_tag, err)
         if (allocated(err)) then
            call test_failed(error, "carbon rejected by tag: "//trim(err%message))
            return
         end if

         call get_radius(6, trim(model_name(imodel)), by_name, err)
         if (allocated(err)) then
            call test_failed(error, "carbon rejected by name: "//trim(err%message))
            return
         end if
         call check(error, by_name, by_tag, thr=0.0_wp, &
                    more="name and tag disagree for "//trim(model_name(imodel)))
         if (allocated(error)) return
      end do

      ! Upper case, and leading/trailing blanks
      call get_radius(6, "CPCM", by_name, err)
      if (allocated(err)) then
         call test_failed(error, "upper-case model name rejected")
         return
      end if
      call get_radius(6, rad_type%cpcm, by_tag, err)
      call check(error, by_name, by_tag, thr=0.0_wp, more="'CPCM' did not resolve to cpcm")
      if (allocated(error)) return

      call get_radius(6, "  smd  ", by_name, err)
      if (allocated(err)) then
         call test_failed(error, "padded model name rejected")
         return
      end if
      call get_radius(6, rad_type%smd, by_tag, err)
      call check(error, by_name, by_tag, thr=0.0_wp, more="'  smd  ' did not resolve to smd")
   end subroutine test_radius_keyword_normalisation

   !> Unknown model names, and names with interior blanks, are rejected
   subroutine test_radius_bad_keyword(error)
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: bad(5) = [character(len=8) :: "", "   ", "xyz", "s md", "cpcm2"]

      type(moist_error_type), allocatable :: err
      integer :: i
      real(wp) :: rad

      do i = 1, size(bad)
         rad = 1.0_wp
         call get_radius(6, trim(bad(i)), rad, err)
         call check(error, allocated(err), more="an invalid model name was accepted")
         if (allocated(error)) return
         call check(error, rad, 0.0_wp, thr=0.0_wp, more="rejected model name left a radius behind")
         if (allocated(error)) return
         deallocate (err)
      end do

      ! Integer tags outside the known set are rejected too
      rad = 1.0_wp
      call get_radius(6, 99, rad, err)
      call check(error, allocated(err), more="an unknown model tag was accepted")
      if (allocated(error)) return
      call check(error, rad, 0.0_wp, thr=0.0_wp, more="rejected model tag left a radius behind")
   end subroutine test_radius_bad_keyword

   !> The seven models must not be aliases of one another. A branch of
   !> fetch_radius wired to the wrong array would otherwise go unnoticed
   subroutine test_radius_models_distinct(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: i, j, iz
      real(wp) :: ri, rj
      logical :: differs

      do i = 1, n_models
         do j = i + 1, n_models
            differs = .false.
            ! 1-88 is inside every model's range
            do iz = 1, 88
               call get_radius(iz, i, ri, err)
               if (allocated(err)) then
                  deallocate (err)
                  cycle
               end if
               call get_radius(iz, j, rj, err)
               if (allocated(err)) then
                  deallocate (err)
                  cycle
               end if
               if (ri /= rj) then
                  differs = .true.
                  exit
               end if
            end do

            if (.not. differs) then
               call test_failed(error, "models "//trim(model_name(i))//" and " &
                                //trim(model_name(j))//" return identical radii")
               return
            end if
         end do
      end do
   end subroutine test_radius_models_distinct

   !> With both a bad symbol and a bad model name, the model name is resolved
   !> first, so its error is the one reported
   subroutine test_radius_error_precedence(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp) :: rad

      call get_radius("Xx", "nosuchmodel", rad, err)
      call check(error, allocated(err), more="a doubly invalid lookup was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "radius type") > 0, &
                 "the model-name error must take precedence over the symbol error")
   end subroutine test_radius_error_precedence

   !> Every overload reports a failed lookup through the error *and* the
   !> negative `missing` sentinel, and leaves the error unallocated on success
   subroutine test_radius_func_reports_error(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp) :: rad

      rad = get_radius_func(0, err)
      call check(error, allocated(err) .and. rad < 0.0_wp, "default overload accepted Z = 0")
      if (allocated(error)) return
      deallocate (err)

      rad = get_radius_func(95, rad_type%d3, err)
      call check(error, allocated(err) .and. rad < 0.0_wp, "tag overload accepted Z past the d3 table")
      if (allocated(error)) return
      deallocate (err)

      rad = get_radius_func(119, "cpcm", err)
      call check(error, allocated(err) .and. rad < 0.0_wp, "name overload accepted Z past the cpcm table")
      if (allocated(error)) return
      deallocate (err)

      rad = get_radius_func(6, "nosuchmodel", err)
      call check(error, allocated(err) .and. rad < 0.0_wp, "name overload accepted an unknown model")
      if (allocated(error)) return
      call check(error, index(err%message, "radius type") > 0, &
                 "the propagated message must name the offending model")
      if (allocated(error)) return
      deallocate (err)

      ! Bondi does not parametrise Tc, and that surfaces through the error too
      rad = get_radius_func(43, "bondi", err)
      call check(error, allocated(err) .and. rad < 0.0_wp, "bondi accepted Tc")
      if (allocated(error)) return
      deallocate (err)

      ! On success the error must be left unallocated
      rad = get_radius_func(6, err)
      call check(error, .not. allocated(err) .and. rad > 0.0_wp, "default overload failed for carbon")
      if (allocated(error)) return
      rad = get_radius_func(6, "cpcm", err)
      call check(error, .not. allocated(err) .and. rad > 0.0_wp, "name overload failed for carbon")
   end subroutine test_radius_func_reports_error

   !* ---------------------------- Group D: solvent tables ---------------------------- *!

   !> Every stored alias resolves to its own solvent, so no alias is shadowed by
   !> an earlier entry or unreachable through normalisation
   subroutine test_solvent_alias_round_trip(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: i, j, id
      real(wp), dimension(max_solvents) :: eps, refr, A, B, g, rho
      integer :: id_list(max_solvents)
      character(len=64) :: name_list(max_solvents)
      character(len=64) :: alias_list(10, max_solvents)

      include "../src/moist/data/solvents.inc"

      do i = 1, max_solvents
         do j = 1, 10
            if (len_trim(alias_list(j, i)) == 0) cycle
            call get_solvent_id(alias_list(j, i), id, err)
            if (allocated(err)) then
               call test_failed(error, "alias '"//trim(alias_list(j, i))//"' does not resolve: "//trim(err%message))
               return
            end if
            call check(error, id, id_list(i), more="alias '"//trim(alias_list(j, i))//"' resolved elsewhere")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_solvent_alias_round_trip

   !> Alias matching ignores case and surrounding blanks
   subroutine test_solvent_alias_normalisation(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: id, reference

      call get_solvent_id("water", reference, err)
      if (allocated(err)) then
         call test_failed(error, "'water' did not resolve: "//trim(err%message))
         return
      end if

      call get_solvent_id("wAtEr", id, err)
      if (allocated(err)) then
         call test_failed(error, "mixed-case alias rejected: "//trim(err%message))
         return
      end if
      call check(error, id, reference, more="mixed-case alias resolved elsewhere")
      if (allocated(error)) return

      call get_solvent_id("  WATER  ", id, err)
      if (allocated(err)) then
         call test_failed(error, "padded alias rejected: "//trim(err%message))
         return
      end if
      call check(error, id, reference, more="padded alias resolved elsewhere")
      if (allocated(error)) return

      ! Stored aliases are normalised like the query
      call get_solvent_id("furan", reference, err)
      if (allocated(err)) then
         call test_failed(error, "'furan' did not resolve: "//trim(err%message))
         return
      end if
      call get_solvent_id(" Tetrole ", id, err)
      if (allocated(err)) then
         call test_failed(error, "'Tetrole' did not resolve: "//trim(err%message))
         return
      end if
      call check(error, id, reference, more="'Tetrole' must map to furan")
   end subroutine test_solvent_alias_normalisation

   !> Solvents with fewer than ten aliases have their remaining alias slots
   !> blank-padded. A blank query must be rejected rather than matching that
   !> padding and silently resolving to whichever solvent comes first
   subroutine test_solvent_blank_alias(error)
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: blank(3) = [character(len=8) :: "", " ", "        "]

      type(moist_error_type), allocatable :: err
      integer :: i, id

      do i = 1, size(blank)
         id = -1
         call get_solvent_id(blank(i), id, err)
         call check(error, allocated(err), more="a blank solvent alias was accepted")
         if (allocated(error)) return
         call check(error, id, 0, more="rejected alias lookup left an id behind")
         if (allocated(error)) return
         deallocate (err)
      end do

      ! And an ordinary unknown alias is still rejected
      call get_solvent_id("definitely-not-a-solvent", id, err)
      call check(error, allocated(err), more="an unknown solvent alias was accepted")
   end subroutine test_solvent_blank_alias

   !> A full solvation system builds for a real solvent and carries the table
   !> values through into the derived type
   subroutine test_solvent_system_constructs(error)
      type(error_type), allocatable, intent(out) :: error

      type(solvation_system_type) :: system
      type(moist_error_type), allocatable :: err
      integer :: water_id

      call get_solvent_id("water", water_id, err)
      if (allocated(err)) then
         call test_failed(error, "could not resolve water: "//trim(err%message))
         return
      end if

      call new_solvation_system(system, water_id, error=err)
      if (allocated(err)) then
         call test_failed(error, "water system failed to build: "//trim(err%message))
         return
      end if

      call check(error, system%solvent_id, water_id, more="constructor stored the wrong id")
      if (allocated(error)) return
      call check(error, trim(system%solvent_name), "water", more="constructor stored the wrong name")
      if (allocated(error)) return
      call check(error, system%solvent_epsilon > 1.0_wp, "water permittivity must exceed vacuum")
      if (allocated(error)) return
      call check(error, system%temperature, 298.15_wp, thr=thr, more="default temperature")
      if (allocated(error)) return
      call check(error, system%pressure_si, 101325.0_wp, thr=thr, more="default pressure")
      if (allocated(error)) return
      call check(error, system%solvent_molar_mass_si > 0.0_wp, "solvent molar mass must be positive")
      if (allocated(error)) return

      ! Mass density in atomic units is number density times molecular mass
      call check(error, system%solvent_mass_density_au, &
                 system%solvent_number_density_au*system%solvent_mass_au, thr=thr, rel=.true., &
                 more="mass density and number density disagree in atomic units")
   end subroutine test_solvent_system_constructs

   !> The constructor validates its inputs before doing any work, and an
   !> unmatched solvent id must be reported instead of leaving the object
   !> half-initialised
   subroutine test_solvent_system_validation(error)
      type(error_type), allocatable, intent(out) :: error

      type(solvation_system_type) :: system
      type(moist_error_type), allocatable :: err

      call new_solvation_system(system, max_solvents + 1, error=err)
      call check(error, allocated(err), more="an unknown solvent id was accepted")
      if (allocated(error)) return
      deallocate (err)

      call new_solvation_system(system, 0, error=err)
      call check(error, allocated(err), more="solvent id 0 was accepted")
      if (allocated(error)) return
      deallocate (err)

      call new_solvation_system(system, 175, temperature=-1.0_wp, error=err)
      call check(error, allocated(err), more="a negative temperature was accepted")
      if (allocated(error)) return
      deallocate (err)

      call new_solvation_system(system, 175, temperature=0.0_wp, error=err)
      call check(error, allocated(err), more="a zero temperature was accepted")
      if (allocated(error)) return
      deallocate (err)

      call new_solvation_system(system, 175, pressure_si=-1.0_wp, error=err)
      call check(error, allocated(err), more="a negative pressure was accepted")
   end subroutine test_solvent_system_validation

   !> Every table entry either builds a solvation system or reports an error
   !> that names its id. Entries without a geometry must not fail with a
   !> garbled message, and the geometry and charge tables must cover the same
   !> ids with the same atom counts
   subroutine test_solvent_system_all_ids(error)
      type(error_type), allocatable, intent(out) :: error

      type(solvation_system_type) :: system
      type(moist_error_type), allocatable :: err, err_charges
      real(wp), allocatable :: charges(:)
      character(len=16) :: id_str
      integer :: id
      real(wp), dimension(max_solvents) :: eps, refr, A, B, g, rho
      integer :: id_list(max_solvents)
      character(len=64) :: name_list(max_solvents)
      character(len=64) :: alias_list(10, max_solvents)

      include "../src/moist/data/solvents.inc"

      do id = 1, max_solvents
         call new_solvation_system(system, id, error=err)
         call get_solvent_charges(id, "gas", "mbis", charges, err_charges)
         if (.not. allocated(err)) then
            call check(error, system%solvent_id, id, more="constructor stored the wrong id")
            if (allocated(error)) return
            call check(error, system%solvent_epsilon, eps(id), thr=thr, more="constructor permittivity")
            if (allocated(error)) return
            call check(error, system%solvent_refractive_index, refr(id), thr=thr, more="constructor refractive index")
            if (allocated(error)) return
            call check(error, system%solvent_alpha, A(id), thr=thr, more="constructor HB acidity")
            if (allocated(error)) return
            call check(error, system%solvent_beta, B(id), thr=thr, more="constructor HB basicity")
            if (allocated(error)) return
            call check(error, system%solvent_surface_tension_si, g(id)*0.001_wp, thr=thr, more="constructor tension")
            if (allocated(error)) return
            call check(error, system%solvent_mass_density_si, rho(id), thr=thr, more="constructor density")
            if (allocated(error)) return
            call check(error, .not. allocated(err_charges), more="solvent builds but has no charges")
            if (allocated(error)) return
            call check(error, size(charges), system%solv_mol%nat, more="geometry and charge atom counts disagree")
            if (allocated(error)) return
            cycle
         end if

         write (id_str, "(i0)") id
         if (index(err%message, "(ID "//trim(id_str)//")") == 0) then
            call test_failed(error, "solvent error does not name its id: "//err%message)
            return
         end if
         call check(error, allocated(err_charges), more="solvent without a geometry has charges")
         if (allocated(error)) return
         deallocate (err)
      end do
   end subroutine test_solvent_system_all_ids

   !> Every charge scheme and MBIS moment set agrees on the atom count, and
   !> failed lookups clear the output
   subroutine test_solvent_charge_data(error)
      type(error_type), allocatable, intent(out) :: error
      type(solvent_multipole_data_type) :: multipoles
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: charges(:)
      integer :: id, ienv, imodel, nat
      character(len=9), parameter :: environments(3) = [character(len=9) :: &
         "gas", "solvent", "conductor"]
      character(len=9), parameter :: models(4) = [character(len=9) :: &
         "hirshfeld", "resp", "mbis", "chelpg"]

      do id = 1, 187
         if (id == 94) cycle
         do ienv = 1, size(environments)
            call get_solvent_multipoles(id, environments(ienv), "mbis", multipoles, err)
            if (allocated(err)) then
               call test_failed(error, "multipole data unavailable: "//err%message)
               return
            end if
            nat = size(multipoles%monopole)
            call check(error, all(shape(multipoles%dipole) == [3, nat]), &
               "dipole shape does not match atom count")
            if (allocated(error)) return
            call check(error, all(shape(multipoles%quadrupole) == [6, nat]), &
               "quadrupole shape does not match atom count")
            if (allocated(error)) return
            call check(error, all(shape(multipoles%octupole) == [10, nat]), &
               "octupole shape does not match atom count")
            if (allocated(error)) return
            do imodel = 1, size(models)
               call get_solvent_charges(id, environments(ienv), models(imodel), charges, err)
               if (allocated(err)) then
                  call test_failed(error, "charge data unavailable: "//err%message)
                  return
               end if
               call check(error, size(charges), nat)
               if (allocated(error)) return
            end do
         end do
      end do

      call get_solvent_charges(94, "gas", "mbis", charges, err)
      call check(error, allocated(err), "excluded ID 94 should return an error")
      if (allocated(error)) return
      deallocate (err)
      call get_solvent_multipoles(94, "gas", "mbis", multipoles, err)
      call check(error, allocated(err), "excluded multipole ID 94 should return an error")
      if (allocated(error)) return
      call check(error, .not. allocated(multipoles%monopole), "excluded multipole ID should clear output")
      if (allocated(error)) return
      deallocate (err)
      call get_solvent_charges(175, "unknown", "mbis", charges, err)
      call check(error, allocated(err), "unknown charge environment should return an error")
      if (allocated(error)) return
      call check(error, .not. allocated(charges), "failed charge environment should clear output")
      if (allocated(error)) return
      deallocate (err)
      call get_solvent_multipoles(175, "unknown", "mbis", multipoles, err)
      call check(error, allocated(err), "unknown environment should return an error")
   end subroutine test_solvent_charge_data

   !> Type-bound accessors return the table entries for their own solvent and
   !> select charge and multipole models independently
   subroutine test_solvent_system_charge_accessors(error)
      type(error_type), allocatable, intent(out) :: error
      type(solvation_system_type) :: system
      type(solvent_multipole_data_type) :: multipoles, ref_multipoles
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: charges(:), ref_charges(:)
      character(len=9), parameter :: models(4) = [character(len=9) :: &
         "hirshfeld", "resp", "mbis", "chelpg"]
      integer :: i

      call new_solvation_system(system, 175, error=err)
      if (allocated(err)) then
         call test_failed(error, "water system unavailable: "//err%message)
         return
      end if

      do i = 1, size(models)
         call get_solvent_charges(175, "gas", models(i), ref_charges, err)
         if (.not. allocated(err)) call system%get_charges("gas", models(i), charges, err)
         if (allocated(err)) then
            call test_failed(error, "charge lookup failed: "//err%message)
            return
         end if
         call check(error, size(charges), size(ref_charges))
         if (allocated(error)) return
         call check(error, all(charges == ref_charges), &
                    "type-bound "//trim(models(i))//" charges differ from the table")
         if (allocated(error)) return
      end do

      call get_solvent_charges(175, "solvent", "mbis", ref_charges, err)
      if (.not. allocated(err)) call system%get_charges(" SoLvEnT ", "MBIS", charges, err)
      if (allocated(err)) then
         call test_failed(error, "case-insensitive charge lookup failed: "//err%message)
         return
      end if
      call check(error, size(charges), size(ref_charges))
      if (allocated(error)) return
      call check(error, all(charges == ref_charges), "case-insensitive lookup returned other charges")
      if (allocated(error)) return

      call get_solvent_multipoles(175, "conductor", "mbis", ref_multipoles, err)
      if (.not. allocated(err)) call system%get_multipoles(" CoNdUcToR ", "MBIS", multipoles, err)
      if (allocated(err)) then
         call test_failed(error, "multipole lookup failed: "//err%message)
         return
      end if
      call check(error, size(multipoles%monopole), size(ref_multipoles%monopole))
      if (allocated(error)) return
      call check(error, all(multipoles%monopole == ref_multipoles%monopole) &
                 .and. all(multipoles%dipole == ref_multipoles%dipole), &
                 "case-insensitive lookup returned other moments")
      if (allocated(error)) return

      call system%get_charges("gas", "unknown", charges, err)
      call check(error, allocated(err), "unknown charge model should return an error")
      if (allocated(error)) return
      call check(error, .not. allocated(charges), "failed charge lookup should clear output")
      if (allocated(error)) return
      deallocate (err)
      call system%get_multipoles("gas", "unknown", multipoles, err)
      call check(error, allocated(err), "unknown multipole model should return an error")
      if (allocated(error)) return
      call check(error, .not. allocated(multipoles%monopole), "failed multipole lookup should clear output")
      if (allocated(error)) return
      deallocate (err)
      call system%get_multipoles("unknown", "mbis", multipoles, err)
      call check(error, allocated(err), "unknown environment should return an error")
   end subroutine test_solvent_system_charge_accessors

   !> Repeated solute updates replace the geometry and recompute the molar
   !> mass from zero
   subroutine test_solute_update_ownership(error)
      type(error_type), allocatable, intent(out) :: error

      type(solvation_system_type) :: system
      type(structure_type) :: solute
      type(moist_error_type), allocatable :: err
      real(wp) :: carbon_mass, hydrogen_mass

      ! Molar masses in kg/mol
      call get_mass(6, carbon_mass, err)
      if (.not. allocated(err)) call get_mass(1, hydrogen_mass, err)
      if (allocated(err)) then
         call test_failed(error, "mass lookup failed: "//trim(err%message))
         return
      end if
      carbon_mass = carbon_mass*0.001_wp
      hydrogen_mass = hydrogen_mass*0.001_wp

      call new_solvation_system(system, 175, error=err)
      if (allocated(err)) then
         call test_failed(error, "water system failed to build: "//trim(err%message))
         return
      end if
      call new_structure(solute, num=[6, 1], sym=["C", "H"], &
                         xyz=reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp], [3, 2]))
      call system%update(solute, err)
      if (allocated(err)) then
         call test_failed(error, "first solute update failed: "//trim(err%message))
         return
      end if
      call check(error, allocated(system%solu_mol), more="update did not allocate the solute")
      if (allocated(error)) return
      call check(error, system%solu_mol%nat, 2, more="first solute atom count")
      if (allocated(error)) return
      call check(error, system%solu_mol%xyz(1, 1), 1.0_wp, thr=thr, more="first solute geometry")
      if (allocated(error)) return
      call check(error, system%solute_molar_mass_si, carbon_mass + hydrogen_mass, thr=thr, rel=.true., &
                 more="first solute molar mass")
      if (allocated(error)) return

      call new_structure(solute, num=[6], sym=["C"], xyz=reshape([7.0_wp, 8.0_wp, 9.0_wp], [3, 1]))
      call system%update(solute, err)
      if (allocated(err)) then
         call test_failed(error, "replacement solute update failed: "//trim(err%message))
         return
      end if
      call check(error, system%solu_mol%nat, 1, more="replacement solute atom count")
      if (allocated(error)) return
      call check(error, system%solu_mol%xyz(1, 1), 7.0_wp, thr=thr, more="replacement solute geometry")
      if (allocated(error)) return
      call check(error, system%solute_molar_mass_si, carbon_mass, thr=thr, rel=.true., &
                 more="molar mass not reset by the replacement")
      if (allocated(error)) return
      call check(error, system%solute_mass_au, carbon_mass/Avogadro_constant/atomic_unit_of_mass, &
                 thr=thr, rel=.true., more="replacement solute mass in atomic units")
   end subroutine test_solute_update_ownership

   !> Water MBIS moments are stored component-major, (component, atom). The
   !> second quadrupole and third octupole component of oxygen differ from
   !> the entries an atom-major reading would return
   subroutine test_water_higher_multipoles(error)
      type(error_type), allocatable, intent(out) :: error

      type(solvent_multipole_data_type) :: multipoles
      type(moist_error_type), allocatable :: err

      call get_solvent_multipoles(175, "gas", "mbis", multipoles, err)
      if (allocated(err)) then
         call test_failed(error, "water gas multipoles unavailable: "//trim(err%message))
         return
      end if
      call check(error, multipoles%quadrupole(2, 1), -5.014708_wp, thr=thr, &
                 more="oxygen quadrupole component 2")
      if (allocated(error)) return
      call check(error, multipoles%octupole(3, 1), 0.653404_wp, thr=thr, &
                 more="oxygen octupole component 3")
   end subroutine test_water_higher_multipoles

end module test_data
