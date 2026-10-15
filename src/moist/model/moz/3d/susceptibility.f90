!> Solvent response of the 3D UV equation on the k-points of any 3D grid
!>
!> - P_ga(|k|) = w_ga(k) + rho_g h_ga(k) of the solvent's VV solution; isotropic, so
!>   one matrix per distinct |k| serves every k-point of that length
!> - Distinct means floating-point equal, without merging; index(j) names the matrix
!>   of k-point j in the grid's k-point order
!> - direct: radial quadrature of stored h(r) at each distinct |k| (susceptibility_at);
!>   avoids k interpolation while retaining radial discretization and cutoff errors
!> - oz: cubic interpolation of native cbar(k), analytical Coulomb and form
!>   factors, then an OZ solve per |k|, including the analytical zero limit
!> - linear: interpolation of the response at the VV grid's own k nodes, at
!>   x = |k|/dk + 1/2 between nodes int(x) and int(x) + 1, held constant outside
!>   them; the tabulation of moist_dev/pysofft (44874e7e), model/rism/3d/uv
!> - The oz charge response is evaluated without low-k cancellation
!> - Every |k| must lie in the VV grid's band, |k| <= pi/dr
module moist_model_moz_3d_susceptibility
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_math_packed_sym, only: packed_index_lut
   use moist_math_sorter, only: qsort
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_model_moz_solvent_vv, only: solvent_vv_type
   use moist_model_moz_potential_sites, only: potential_sites_type
   implicit none(type, external)
   private

   public :: uv3d_susceptibility_type, new_uv3d_susceptibility

   !> Solvent response per distinct |k| of a 3D grid
   type :: uv3d_susceptibility_type
      !> Solvent sites
      integer :: nv = 0
      !> Response matrices P(g, a, shell), (nv, nv, nshell)
      real(wp), allocatable :: values(:, :, :)
      !> Matrix of every k-point, (npts_k)
      integer, allocatable :: index(:)
      !> Wave number of every matrix, 1/bohr, ascending, (nshell)
      real(wp), allocatable :: radius(:)
      !> Stable chi q / k^2 per shell, including its zero limit (oz mode)
      real(wp), allocatable :: charge_response(:, :), charge(:)
      !> Inputs of the cached table
      integer(int64), private :: bulk_generation = -1_int64
      character(len=8), private :: mode = ""
   contains
      !> Whether this table still matches the bulk solution and ordered grid radii
      procedure :: matches => response_matches
   end type uv3d_susceptibility_type

contains

   !> Tabulate the solvent response on the k-points of an updated 3D grid
   !>
   !> @param[out] self     response table
   !> @param[in]  grid     updated grid with its k-points
   !> @param[in]  solvent  solved bulk solvent
   !> @param[in]  mode     oz, direct or linear
   !> @param[out] error    grid without k-points, unknown mode, unsolved solvent or |k| out of band
   subroutine new_uv3d_susceptibility(self, grid, solvent, mode, error)
      !> Response table
      type(uv3d_susceptibility_type), intent(out) :: self
      !> Updated grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Solved solvent
      type(solvent_vv_type), intent(in) :: solvent
      !> Mode
      character(len=*), intent(in) :: mode
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(packed_index_lut) :: lut
      type(potential_sites_type) :: sites
      real(wp), allocatable :: radii(:), chi(:, :), z(:, :)
      integer, allocatable :: order(:)
      integer :: nk, j, nshell, g, a

      nk = grid%npts_k
      if (nk < 1) then
         call fatal_error(error, "3D UV susceptibility needs an updated grid with k-points")
         return
      end if
      allocate (radii(nk), order(nk), self%index(nk))
      do j = 1, nk
         radii(j) = norm2(grid%kpoint(j))
         order(j) = j
      end do
      if (.not. all(ieee_is_finite(radii))) then
         call fatal_error(error, "3D UV susceptibility wave numbers must be finite")
         return
      end if
      call qsort(radii, order, error)
      if (allocated(error)) return
      nshell = 0
      do j = 1, nk
         if (j == 1) then
            nshell = 1
         else if (radii(j) /= radii(nshell)) then
            nshell = nshell + 1
            radii(nshell) = radii(j)
         end if
         self%index(order(j)) = nshell
      end do
      self%radius = radii(:nshell)

      select case (to_lower(trim(adjustl(mode))))
      case ("oz")
         call solvent%susceptibility_oz_at(self%radius, chi, error, z)
         if (allocated(error)) return
         self%charge_response = transpose(z)
         call solvent%potential%sites(sites, error)
         if (allocated(error)) return
         allocate (self%charge(solvent%nsite()), source=0.0_wp)
         if (allocated(sites%q)) self%charge = sites%q
      case ("direct")
         call solvent%susceptibility_at(self%radius, chi, error)
      case ("linear")
         call interpolated_response(solvent, self%radius, chi, error)
      case default
         call fatal_error(error, "Unknown 3D UV susceptibility mode: "//trim(mode)//" (choose oz, direct or linear)")
      end select
      if (allocated(error)) return

      self%nv = solvent%nsite()
      call lut%init(self%nv)
      allocate (self%values(self%nv, self%nv, nshell))
      do j = 1, nshell
         do a = 1, self%nv
            do g = 1, self%nv
               self%values(g, a, j) = chi(j, lut%idx(g, a))
            end do
         end do
      end do
      self%bulk_generation = solvent%solution_generation()
      self%mode = to_lower(trim(adjustl(mode)))
   end subroutine new_uv3d_susceptibility

   !> Cache identity includes every successful bulk solve, mode and ordered |k|
   logical function response_matches(self, grid, solvent, mode) result(matches)
      class(uv3d_susceptibility_type), intent(in) :: self
      class(moist_math_grid_3d_type), intent(in) :: grid
      type(solvent_vv_type), intent(in) :: solvent
      character(len=*), intent(in) :: mode
      integer :: j
      matches = .false.
      if (.not. solvent%is_solved()) return
      if (.not. allocated(self%index) .or. .not. allocated(self%radius)) return
      if (self%bulk_generation /= solvent%solution_generation() .or. &
         & self%mode /= to_lower(trim(adjustl(mode))) .or. size(self%index) /= grid%npts_k) return
      do j = 1, grid%npts_k
         if (self%radius(self%index(j)) /= norm2(grid%kpoint(j))) return
      end do
      matches = .true.
   end function response_matches

   !> Response at the given wave numbers, interpolated between the VV grid's k nodes
   !>
   !> @param[in]  solvent  solved bulk solvent
   !> @param[in]  k        wave numbers, 1/bohr (nk)
   !> @param[out] chi      response, packed (nk, npair)
   !> @param[out] error    unsolved solvent or wave number out of band
   subroutine interpolated_response(solvent, k, chi, error)
      !> Solved solvent
      type(solvent_vv_type), intent(in) :: solvent
      !> Wave numbers
      real(wp), intent(in) :: k(:)
      !> Response
      real(wp), allocatable, intent(out) :: chi(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: nodes(:, :)
      real(wp) :: dk, x, frac
      integer :: j, ir, nr

      call solvent%radial_grid(rgrid, kgrid, error)
      if (allocated(error)) return
      call solvent%susceptibility(nodes, error)
      if (allocated(error)) return
      nr = size(nodes, 1)
      dk = kgrid%w(1)
      if (any(k < 0.0_wp) .or. any(k > real(nr, wp)*dk)) then
         call fatal_error(error, "3D UV susceptibility: a wave number lies outside the VV grid's band")
         return
      end if
      allocate (chi(size(k), size(nodes, 2)))
      do j = 1, size(k)
         x = k(j)/dk + 0.5_wp
         ir = int(x)
         if (ir < 1) then
            chi(j, :) = nodes(1, :)
         else if (ir >= nr) then
            chi(j, :) = nodes(nr, :)
         else
            frac = x - real(ir, wp)
            chi(j, :) = (1.0_wp - frac)*nodes(ir, :) + frac*nodes(ir + 1, :)
         end if
      end do
   end subroutine interpolated_response

end module moist_model_moz_3d_susceptibility
