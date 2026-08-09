!> Tests for the parameter tables under src/moist/data/ and the accessors
module test_data
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type, fatal_error
   use mctc_io_convert, only: aatoau
   use mctc_io_symbols, only: to_symbol
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

   use moist_data_en, only: get_electronegativity, en_max_elem => max_elem
   use moist_data_hardness, only: get_hardness, hardness_max_elem => max_elem
   use moist_data_mass, only: get_mass, mass_max_elem => max_elem
   use moist_data_atomicrad, only: get_atomic_rad, get_covalent_rad, &
      & arad_max_elem => max_elem
   use moist_data_radii_legacy, only: get_radius, get_radius_func, &
      & get_upper_bound, rad_type
   use moist_data_solvents, only: get_solvent_id, get_solvent_for_alpb, max_solvents, &
      & solvation_system_parameters, new_solvation_system_parameters

   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none
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

   !> Atomic numbers sitting on the array-constructor row boundaries of the elem. tab
   integer, parameter :: anchor_z(8) = [1, 10, 18, 30, 36, 54, 70, 118]

contains

   !> Collect all data-table tests.
   subroutine collect_data(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ElementTablesComplete", test_element_tables_complete), &
                  new_unittest("ElementTablesPhysical", test_element_tables_physical), &
                  new_unittest("HardnessSuperheavyAreZero", test_hardness_superheavy), &
                  new_unittest("AnchorElectronegativity", test_anchor_en), &
                  new_unittest("AnchorHardness", test_anchor_hardness), &
                  new_unittest("AnchorMass", test_anchor_mass), &
                  new_unittest("AnchorAtomicRad", test_anchor_atomic_rad), &
                  new_unittest("AnchorCovalentRad", test_anchor_covalent_rad), &
                  new_unittest("AccessorRejectsOutOfRange", test_accessor_out_of_range), &
                  new_unittest("AccessorAcceptsBoundaries", test_accessor_boundaries), &
                  new_unittest("AccessorRejectsBadSymbol", test_accessor_bad_symbol), &
                  new_unittest("AccessorSymbolCaseInsensitive", test_accessor_symbol_case), &
                  new_unittest("AccessorSymbolMatchesNumber", test_accessor_symbol_consistency), &
                  new_unittest("RadiusAllModelsComplete", test_radius_models_complete), &
                  new_unittest("RadiusModelUpperBounds", test_radius_upper_bounds), &
                  new_unittest("RadiusKeywordNormalisation", test_radius_keyword_normalisation), &
                  new_unittest("RadiusRejectsBadKeyword", test_radius_bad_keyword), &
                  new_unittest("RadiusModelsAreDistinct", test_radius_models_distinct), &
                  new_unittest("RadiusBondiMissingRejected", test_radius_bondi_missing), &
                  new_unittest("RadiusModelErrorWinsOverSymbol", test_radius_error_precedence), &
                  new_unittest("RadiusFuncSentinel", test_radius_func_sentinel), &
                  new_unittest("RadiusFuncReportsError", test_radius_func_reports_error), &
                  new_unittest("SolventIdsAreContiguous", test_solvent_ids_contiguous), &
                  new_unittest("SolventNameRoundTrip", test_solvent_name_round_trip), &
                  new_unittest("SolventAliasCaseAndBlanks", test_solvent_alias_normalisation), &
                  new_unittest("SolventRejectsBlankAlias", test_solvent_blank_alias), &
                  new_unittest("SolventRejectsBadId", test_solvent_bad_id), &
                  new_unittest("SolventSystemConstructs", test_solvent_system_constructs), &
                  new_unittest("SolventSystemValidatesInput", test_solvent_system_validation) &
                  ]
   end subroutine collect_data

    !* -------------------------------- Private helpers ------------------------------- *!

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
   !> not this one, because the shifted-in filler is usually zero.
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
   !> callers must branch on the error, not on the value.
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

   !* ------------ Group B: anchor values on the constructor row boundaries ----------- *!

   !> Pauling electronegativities at H, Ne, Ar, Zn, Kr, Xe, Yb, Og.
   subroutine test_anchor_en(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: expected(8) = [ &
         & 2.20_wp, 4.50_wp, 3.50_wp, 1.65_wp, 3.00_wp, 2.60_wp, 1.26_wp, 1.50_wp]

      type(moist_error_type), allocatable :: err
      integer :: i
      real(wp) :: en

      do i = 1, size(anchor_z)
         call get_electronegativity(anchor_z(i), en, err)
         if (allocated(err)) then
            call test_failed(error, "anchor lookup failed: "//trim(err%message))
            return
         end if
         call check(error, en, expected(i), thr=thr, more="electronegativity anchor mismatch")
         if (allocated(error)) return
      end do
   end subroutine test_anchor_en

   !> DFT-D4 chemical hardnesses on the same anchors.
   subroutine test_anchor_hardness(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: expected(8) = [ &
         & 0.47259288_wp, 0.75191607_wp, 0.47308269_wp, 0.27592565_wp, &
         & 0.46611708_wp, 0.44105777_wp, 0.31159587_wp, 0.00000000_wp]

      type(moist_error_type), allocatable :: err
      integer :: i
      real(wp) :: eta

      do i = 1, size(anchor_z)
         call get_hardness(anchor_z(i), eta, err)
         if (allocated(err)) then
            call test_failed(error, "anchor lookup failed: "//trim(err%message))
            return
         end if
         call check(error, eta, expected(i), thr=thr, more="hardness anchor mismatch")
         if (allocated(error)) return
      end do
   end subroutine test_anchor_hardness

   !> NIST atomic masses in u on the same anchors.
   subroutine test_anchor_mass(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: expected(8) = [ &
         &   1.00794075_wp, 20.18004638_wp, 39.94779856_wp, 65.37778253_wp, &
         &  83.79800000_wp, 131.29276145_wp, 173.05415017_wp, 294.21392000_wp]

      type(moist_error_type), allocatable :: err
      integer :: i
      real(wp) :: mass

      do i = 1, size(anchor_z)
         call get_mass(anchor_z(i), mass, err)
         if (allocated(err)) then
            call test_failed(error, "anchor lookup failed: "//trim(err%message))
            return
         end if
         call check(error, mass, expected(i), thr=thr, more="mass anchor mismatch")
         if (allocated(error)) return
      end do
   end subroutine test_anchor_mass

   !> Mantina/Truhlar atomic radii. The table is stored in bohr, so the
   !> Angstrom source values are converted here the same way.
   subroutine test_anchor_atomic_rad(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: expected(8) = aatoau*[ &
         & 0.32_wp, 0.62_wp, 1.01_wp, 1.20_wp, 1.16_wp, 1.36_wp, 1.78_wp, 1.57_wp]

      type(moist_error_type), allocatable :: err
      integer :: i
      real(wp) :: rad

      do i = 1, size(anchor_z)
         call get_atomic_rad(anchor_z(i), rad, err)
         if (allocated(err)) then
            call test_failed(error, "anchor lookup failed: "//trim(err%message))
            return
         end if
         call check(error, rad, expected(i), thr=thr, more="atomic radius anchor mismatch")
         if (allocated(error)) return
      end do
   end subroutine test_anchor_atomic_rad

   !> Alvarez 2008 covalent radii, likewise stored in bohr.
   subroutine test_anchor_covalent_rad(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: expected(8) = aatoau*[ &
         & 0.31_wp, 0.58_wp, 1.06_wp, 1.22_wp, 1.16_wp, 1.40_wp, 1.87_wp, 1.50_wp]

      type(moist_error_type), allocatable :: err
      integer :: i
      real(wp) :: rad

      do i = 1, size(anchor_z)
         call get_covalent_rad(anchor_z(i), rad, err)
         if (allocated(err)) then
            call test_failed(error, "anchor lookup failed: "//trim(err%message))
            return
         end if
         call check(error, rad, expected(i), thr=thr, more="covalent radius anchor mismatch")
         if (allocated(error)) return
      end do
   end subroutine test_anchor_covalent_rad

   !* ------------------------- Group C: the accessor contract ------------------------ *!

   !> Out-of-range atomic numbers must be rejected, and the output left at zero
   !> rather than carrying a sentinel the caller might mistake for data.
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

   !> Both ends of the valid range must be accepted. This is the off-by-one
   !> guard: Z = max_elem is data, Z = max_elem + 1 is not.
   subroutine test_accessor_boundaries(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: tag, upper
      real(wp) :: val

      do tag = 1, n_accessors
         upper = accessor_max_elem(tag)

         call lookup_num(tag, 1, val, err)
         call check(error, .not. allocated(err), trim(acc_name(tag))//" rejected Z = 1")
         if (allocated(error)) return

         call lookup_num(tag, upper, val, err)
         call check(error, .not. allocated(err), trim(acc_name(tag))//" rejected its last element")
         if (allocated(error)) return

         call lookup_num(tag, upper + 1, val, err)
         call check(error, allocated(err), trim(acc_name(tag))//" accepted max_elem + 1")
         if (allocated(error)) return
         deallocate (err)
      end do
   end subroutine test_accessor_boundaries

   !> Unknown, empty and blank symbols must all be rejected. Previously these
   !> resolved to atomic number zero and fell through to a silent sentinel.
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
   !> deuterium/tritium aliases resolve to hydrogen.
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

      ! Mixed case on a two-letter symbol.
      call get_electronegativity("hE", val, err)
      if (allocated(err)) then
         call test_failed(error, "mixed-case symbol rejected: "//trim(err%message))
         return
      end if
      call get_electronegativity(2, reference, err)
      call check(error, val, reference, thr=0.0_wp, more="'hE' did not resolve to helium")
   end subroutine test_accessor_symbol_case

   !> The symbol and atomic-number overloads must agree for every element.
   !> This walks every entry of every table through both paths without
   !> restating a single tabulated value.
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

   !* --------------------- Group D: radius-model keyword dispatch -------------------- *!

   !> Every model must answer for every atomic number up to its own bound. The
   !> radius tables have four different lengths (88, 94, 96, 118) whose declared
   !> extents are decoupled from the max_elem_* constants used to guard them, so
   !> a mismatch would otherwise be an unchecked out-of-bounds read.
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
               ! sentinel; every other model must be complete.
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

   !> Each model's documented upper bound is accepted and one past it rejected.
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
      ! sentinel the caller could mistake for a real bound.
      call get_upper_bound(n_models + 1, upper, err)
      call check(error, allocated(err), more="an unknown model tag was accepted")
      if (allocated(error)) return
      call check(error, upper, -1, more="rejected bound lookup left a usable value")
   end subroutine test_radius_upper_bounds

   !> Model names are matched case-insensitively after trimming and adjusting,
   !> and every name must agree with its integer tag.
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

      ! Upper case, and leading/trailing blanks.
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

   !> Unknown model names, and names with interior blanks, are rejected.
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

      ! Integer tags outside the known set are rejected too.
      rad = 1.0_wp
      call get_radius(6, 99, rad, err)
      call check(error, allocated(err), more="an unknown model tag was accepted")
      if (allocated(error)) return
      call check(error, rad, 0.0_wp, thr=0.0_wp, more="rejected model tag left a radius behind")
   end subroutine test_radius_bad_keyword

   !> The seven models must not be aliases of one another. A branch of
   !> fetch_radius wired to the wrong array would otherwise go unnoticed, since
   !> its default case silently falls through to the CPCM table.
   subroutine test_radius_models_distinct(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: i, j, iz
      real(wp) :: ri, rj
      logical :: differs

      do i = 1, n_models
         do j = i + 1, n_models
            differs = .false.
            ! 1-88 is inside every model's range.
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

   !> Bondi has no radius for the mid-row transition metals. Those entries carry
   !> a negative sentinel in the table and must surface as an error, never as a
   !> negative radius handed back to the caller.
   subroutine test_radius_bondi_missing(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: iz, nmissing
      real(wp) :: rad

      nmissing = 0
      do iz = 1, 88
         call get_radius(iz, rad_type%bondi, rad, err)
         if (allocated(err)) then
            nmissing = nmissing + 1
            deallocate (err)
            cycle
         end if
         if (rad <= 0.0_wp) then
            call test_failed(error, "bondi returned a non-positive radius without an error")
            return
         end if
      end do

      call check(error, nmissing > 0, "bondi is expected to have unparametrised elements")
      if (allocated(error)) return
      ! Technetium is one of the documented gaps.
      call get_radius(43, rad_type%bondi, rad, err)
      call check(error, allocated(err), more="bondi accepted Tc, which it does not parametrise")
   end subroutine test_radius_bondi_missing

   !> With both a bad symbol and a bad model name, the model name is resolved
   !> first, so its error is the one reported.
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

   !> A rejected lookup sets the error *and* returns the negative sentinel, so
   !> code that only inspects the value (print_static_radii skips unparametrised
   !> rows this way) still sees the failure.
   subroutine test_radius_func_sentinel(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp) :: rad

      rad = get_radius_func(0, err)
      call check(error, rad < 0.0_wp, "get_radius_func must report failure as a negative value")
      if (allocated(error)) return

      rad = get_radius_func(119, "cpcm", err)
      call check(error, rad < 0.0_wp, "get_radius_func must reject Z past the table")
      if (allocated(error)) return

      rad = get_radius_func(6, "nosuchmodel", err)
      call check(error, rad < 0.0_wp, "get_radius_func must reject an unknown model name")
      if (allocated(error)) return

      rad = get_radius_func(6, err)
      call check(error, rad > 0.0_wp, "get_radius_func must succeed for carbon")
   end subroutine test_radius_func_sentinel

   !> Passing the optional error reports the same failures the subroutine form
   !> would, across all three overloads, and leaves it unallocated on success.
   subroutine test_radius_func_reports_error(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp) :: rad

      rad = get_radius_func(0, err)
      call check(error, allocated(err), more="default overload accepted Z = 0")
      if (allocated(error)) return
      call check(error, rad < 0.0_wp, "the sentinel must still be returned alongside the error")
      if (allocated(error)) return
      deallocate (err)

      rad = get_radius_func(95, rad_type%d3, err)
      call check(error, allocated(err), more="tag overload accepted Z past the d3 table")
      if (allocated(error)) return
      deallocate (err)

      rad = get_radius_func(6, "nosuchmodel", err)
      call check(error, allocated(err), more="name overload accepted an unknown model")
      if (allocated(error)) return
      call check(error, index(err%message, "radius type") > 0, &
                 "the propagated message must name the offending model")
      if (allocated(error)) return
      deallocate (err)

      ! Bondi does not parametrise Tc, and that surfaces through the error too.
      rad = get_radius_func(43, "bondi", err)
      call check(error, allocated(err), more="bondi accepted Tc")
      if (allocated(error)) return
      deallocate (err)

      ! On success the error must be left unallocated.
      rad = get_radius_func(6, "cpcm", err)
      call check(error, .not. allocated(err), "a valid lookup must not raise an error")
      if (allocated(error)) return
      call check(error, rad > 0.0_wp, "carbon must have a positive CPCM radius")
   end subroutine test_radius_func_reports_error

   !* ---------------------------- Group E: solvent tables ---------------------------- *!

   !> Solvent ids run 1..max_solvents without gaps, every entry has a name, and
   !> every permittivity is physical. The lookups assume this identity mapping.
   subroutine test_solvent_ids_contiguous(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      character(:), allocatable :: name
      integer :: id
      real(wp) :: eps_val

      do id = 1, max_solvents
         call get_solvent_for_alpb(id, eps_val, name, err)
         if (allocated(err)) then
            call test_failed(error, "solvent id gap in the table: "//trim(err%message))
            return
         end if
         if (.not. allocated(name)) then
            call test_failed(error, "solvent name was not returned")
            return
         end if
         if (len_trim(name) == 0) then
            call test_failed(error, "solvent has a blank name")
            return
         end if
         if (eps_val < 1.0_wp) then
            call test_failed(error, "solvent permittivity below the vacuum limit")
            return
         end if
      end do
   end subroutine test_solvent_ids_contiguous

   !> Every solvent's own name must resolve back to its own id. This walks all
   !> 180 entries and would fail if a name were misspelt relative to its alias
   !> list, or if two solvents shared an alias and the wrong one won.
   subroutine test_solvent_name_round_trip(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      character(:), allocatable :: name
      integer :: id, resolved
      real(wp) :: eps_val

      do id = 1, max_solvents
         call get_solvent_for_alpb(id, eps_val, name, err)
         if (allocated(err)) then
            call test_failed(error, "solvent lookup failed: "//trim(err%message))
            return
         end if

         call get_solvent_id(name, resolved, err)
         if (allocated(err)) then
            call test_failed(error, "solvent name '"//name//"' does not resolve: "//trim(err%message))
            return
         end if
         if (resolved /= id) then
            call test_failed(error, "solvent name '"//name//"' resolved to the wrong id")
            return
         end if
      end do
   end subroutine test_solvent_name_round_trip

   !> Alias matching ignores case and surrounding blanks.
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

      ! A secondary alias reaches the same solvent as the primary name.
      call get_solvent_id("methyl chloroform", id, err)
      if (allocated(err)) then
         call test_failed(error, "secondary alias rejected: "//trim(err%message))
         return
      end if
      call check(error, id, 1, more="'methyl chloroform' must map to 1,1,1-trichloroethane")
   end subroutine test_solvent_alias_normalisation

   !> Solvents with fewer than ten aliases have their remaining alias slots
   !> blank-padded. A blank query must be rejected rather than matching that
   !> padding and silently resolving to whichever solvent comes first.
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

      ! And an ordinary unknown alias is still rejected.
      call get_solvent_id("definitely-not-a-solvent", id, err)
      call check(error, allocated(err), more="an unknown solvent alias was accepted")
   end subroutine test_solvent_blank_alias

   !> Ids outside the table are rejected by the permittivity lookup.
   subroutine test_solvent_bad_id(error)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      character(:), allocatable :: name
      integer :: i
      integer :: bad_id(4)
      real(wp) :: eps_val

      bad_id = [0, -1, max_solvents + 1, huge(1)]

      do i = 1, size(bad_id)
         call get_solvent_for_alpb(bad_id(i), eps_val, name, err)
         call check(error, allocated(err), more="an out-of-range solvent id was accepted")
         if (allocated(error)) return
         deallocate (err)
      end do
   end subroutine test_solvent_bad_id

   !> A full solvation system builds for a real solvent and carries the table
   !> values through into the derived type.
   subroutine test_solvent_system_constructs(error)
      type(error_type), allocatable, intent(out) :: error

      type(solvation_system_parameters) :: system
      type(moist_error_type), allocatable :: err
      integer :: water_id

      call get_solvent_id("water", water_id, err)
      if (allocated(err)) then
         call test_failed(error, "could not resolve water: "//trim(err%message))
         return
      end if

      call new_solvation_system_parameters(system, water_id, error=err)
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
   end subroutine test_solvent_system_constructs

   !> The constructor validates its inputs before doing any work, and an
   !> unmatched solvent id must be reported instead of leaving the object
   !> half-initialised.
   subroutine test_solvent_system_validation(error)
      type(error_type), allocatable, intent(out) :: error

      type(solvation_system_parameters) :: system
      type(moist_error_type), allocatable :: err

      call new_solvation_system_parameters(system, max_solvents + 1, error=err)
      call check(error, allocated(err), more="an unknown solvent id was accepted")
      if (allocated(error)) return
      deallocate (err)

      call new_solvation_system_parameters(system, 0, error=err)
      call check(error, allocated(err), more="solvent id 0 was accepted")
      if (allocated(error)) return
      deallocate (err)

      call new_solvation_system_parameters(system, 175, temperature=-1.0_wp, error=err)
      call check(error, allocated(err), more="a negative temperature was accepted")
      if (allocated(error)) return
      deallocate (err)

      call new_solvation_system_parameters(system, 175, temperature=0.0_wp, error=err)
      call check(error, allocated(err), more="a zero temperature was accepted")
      if (allocated(error)) return
      deallocate (err)

      call new_solvation_system_parameters(system, 175, pressure_si=-1.0_wp, error=err)
      call check(error, allocated(err), more="a negative pressure was accepted")
      if (allocated(error)) return
      deallocate (err)

      ! The error argument is optional. Omitting it on a failing call must still
      ! return cleanly: the routine has to route through its own local error
      ! rather than probing an absent optional.
      call new_solvation_system_parameters(system, max_solvents + 1, error=err)
      call check(error, allocated(err), more="an unknown solvent id was accepted without an error argument")
      if (allocated(error)) return
   end subroutine test_solvent_system_validation

end module test_data
