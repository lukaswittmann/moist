!> Fixed 1D UV operators, reused independently of host potential tables
module moist_model_moz_1d_operators
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_packed_sym, only: packed_index_lut, npair_from_ns
   use moist_model_moz_solvent_vv, only: solvent_vv_type
   use moist_model_moz_solvent_vv_kernel, only: vv_intramolecular
   use moist_model_moz_potential_sites, only: potential_sites_type
   use moist_model_moz_1d_kernel, only: uv_charge_moment_correction
   implicit none(type, external)
   private
   public :: uv_operators_type, same_radial_grid

   !> Owned response and intramolecular tables on one reciprocal radial grid
   type :: uv_operators_type
      !> Charge-moment-corrected solvent response P_ba, (npts, nv, nv)
      real(wp), allocatable :: response(:, :, :)
      !> Packed symmetric solute W(k), (npts, nu(nu+1)/2)
      real(wp), allocatable :: wu(:, :)
      !> Exact inputs of the cached tables
      type(moist_math_grid_radial_type), private :: kgrid
      !> Solvent solution generation of the cached response; -1 before the first build
      integer(int64), private :: bulk_generation = -1_int64
      !> Solute coordinates of the cached W, bohr, (3, nu)
      real(wp), allocatable, private :: xyz(:, :)
   contains
      !> Rebuild only tables whose bulk, grid or geometry inputs changed
      procedure :: prepare => operators_prepare
      !> Whether the cached response matches a bulk solution and grid
      procedure :: response_matches
      !> Whether the cached W matches solute coordinates and a grid
      procedure :: intramolecular_matches
   end type uv_operators_type

contains

   !> Prepare fixed UV operators without changing the sampling mathematics
   !>
   !> - New tables replace cached data only after preparation succeeds
   !>
   !> @param[inout] self     operators; cached tables rebuilt when inputs changed
   !> @param[in]    kgrid    reciprocal radial grid
   !> @param[in]    solvent  solved bulk solvent
   !> @param[in]    xyz      solute site coordinates, bohr, (3, nu)
   !> @param[out]   error    unsolved solvent, bad grid or coordinates, or table failure
   subroutine operators_prepare(self, kgrid, solvent, xyz, error)
      !> Operators
      class(uv_operators_type), intent(inout) :: self
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Solved bulk solvent
      type(solvent_vv_type), intent(in) :: solvent
      !> Solute site coordinates, bohr, (3, nu)
      real(wp), intent(in) :: xyz(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(packed_index_lut) :: lut
      type(potential_sites_type) :: sites
      real(wp), allocatable :: chi(:, :), charge_response(:, :), response(:, :, :), wu(:, :)
      integer :: nu, nv, a, b
      logical :: rebuild_response, rebuild_w, applied

      if (.not. solvent%is_solved()) then
         call fatal_error(error, "Solve the bulk solvent before preparing 1D UV operators")
         return
      end if
      if (.not. allocated(kgrid%r) .or. .not. allocated(kgrid%w)) then
         call fatal_error(error, "1D UV operators need a built reciprocal radial grid")
         return
      end if
      if (kgrid%npts < 1 .or. size(kgrid%r) /= kgrid%npts .or. size(kgrid%w) /= kgrid%npts) then
         call fatal_error(error, "1D UV operator grid has inconsistent sizes")
         return
      end if
      if (.not. all(ieee_is_finite(kgrid%r)) .or. .not. all(ieee_is_finite(kgrid%w))) then
         call fatal_error(error, "1D UV operator grid nodes and weights must be finite")
         return
      end if
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) < 1 .or. .not. all(ieee_is_finite(xyz))) then
         call fatal_error(error, "1D UV operators need finite solute coordinates (3, nu)")
         return
      end if
      rebuild_response = .not. self%response_matches(kgrid, solvent)
      rebuild_w = .not. self%intramolecular_matches(kgrid, xyz)
      if (.not. rebuild_response .and. .not. rebuild_w) return

      if (rebuild_response) then
         if (solvent%theory_name() == "DRISM") then
            call solvent%susceptibility_oz_at(kgrid%r, chi, error, charge_response)
         else
            call solvent%susceptibility_at(kgrid%r, chi, error)
         end if
         if (allocated(error)) return
         call solvent%potential%sites(sites, error)
         if (allocated(error)) return
         nv = solvent%nsite()
         call lut%init(nv)
         allocate (response(size(chi, 1), nv, nv))
         do a = 1, nv
            do b = 1, nv
               response(:, b, a) = chi(:, lut%idx(a, b))
            end do
         end do
         if (allocated(sites%q)) then
            call uv_charge_moment_correction(kgrid%r, sites%q, response, applied, charge_response)
         end if
      end if
      if (rebuild_w) then
         nu = size(xyz, 2)
         call lut%init(nu)
         allocate (wu(kgrid%npts, npair_from_ns(nu)))
         call vv_intramolecular(kgrid%r, xyz, lut, wu)
      end if

      if (rebuild_response) then
         call move_alloc(response, self%response)
         self%bulk_generation = solvent%solution_generation()
      end if
      if (rebuild_w) then
         call move_alloc(wu, self%wu)
         self%xyz = xyz
      end if
      if (.not. same_radial_grid(self%kgrid, kgrid)) self%kgrid = kgrid
   end subroutine operators_prepare

   !> Whether the corrected response still represents this bulk solution and grid
   !>
   !> @param[in] self     operators
   !> @param[in] kgrid    reciprocal radial grid
   !> @param[in] solvent  bulk solvent
   pure logical function response_matches(self, kgrid, solvent)
      !> Operators
      class(uv_operators_type), intent(in) :: self
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Bulk solvent
      type(solvent_vv_type), intent(in) :: solvent

      response_matches = .false.
      if (.not. allocated(self%response) .or. .not. solvent%is_solved()) return
      if (self%bulk_generation /= solvent%solution_generation()) return
      response_matches = same_radial_grid(self%kgrid, kgrid)
   end function response_matches

   !> Whether W still represents these solute coordinates and this grid
   !>
   !> @param[in] self   operators
   !> @param[in] kgrid  reciprocal radial grid
   !> @param[in] xyz    solute site coordinates, bohr, (3, nu)
   pure logical function intramolecular_matches(self, kgrid, xyz)
      !> Operators
      class(uv_operators_type), intent(in) :: self
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Solute site coordinates, bohr, (3, nu)
      real(wp), intent(in) :: xyz(:, :)

      intramolecular_matches = .false.
      if (.not. allocated(self%wu) .or. .not. allocated(self%xyz)) return
      if (any(shape(xyz) /= shape(self%xyz))) return
      if (.not. all(xyz == self%xyz)) return
      intramolecular_matches = same_radial_grid(self%kgrid, kgrid)
   end function intramolecular_matches

   !> Exact equality of radial nodes, weights and transform geometry
   !>
   !> @param[in] a  first radial grid
   !> @param[in] b  second radial grid
   pure logical function same_radial_grid(a, b)
      !> First and second radial grid
      type(moist_math_grid_radial_type), intent(in) :: a, b

      same_radial_grid = .false.
      if (a%npts /= b%npts .or. a%transform /= b%transform .or. a%spacing /= b%spacing) return
      if (.not. allocated(a%r) .or. .not. allocated(b%r) .or. &
         & .not. allocated(a%w) .or. .not. allocated(b%w)) return
      if (any(shape(a%r) /= shape(b%r)) .or. any(shape(a%w) /= shape(b%w))) return
      same_radial_grid = all(a%r == b%r) .and. all(a%w == b%w)
   end function same_radial_grid

end module moist_model_moz_1d_operators
