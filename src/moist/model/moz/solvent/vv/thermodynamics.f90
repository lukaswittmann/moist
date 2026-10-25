!> Bulk thermodynamics of a 1D VV solution of one rigid molecular species
!>
!> - Excess chemical potential per molecule from the closure's local energy integrand
!> - mu = (rho/beta) sum_ab int Phi_ab(h, c) d^3r, with physical h and c
!> - Compressibility route from the direct correlation at k = 0,
!>   C = sum_ab c_ab(k = 0) = sum_ab int c_ab d^3r
!> - rho kT kappa_T = 1/(1 - rho C)
!> - P = (n + 1)/2 rho kT - kT/2 rho^2 C for n sites
!> - Reference: Sergiievskyi, Jeanmairet, Levesque, Borgis, JCP 143, 184116 (2015), Eq. 21
!> - Derived there for HNC, evaluated unchanged for the other closures
!> - Ideal-gas pressure rho kT of the same density, for the PC+ correction
!> - Every site carries the molecular density rho
!> - Sums run over all ordered site pairs, so off-diagonal packed pairs count twice
!> - Coulomb tails of c cancel pointwise in these sums for a neutral molecule
module moist_model_moz_solvent_vv_thermodynamics
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_packed_sym, only: packed_index_lut
   use moist_model_moz_closure, only: moist_model_moz_closure_type
   use moist_model_moz_thermodynamics_list, only: thermo_list_type
   use moist_model_moz_thermodynamics_quantities, only: new_bulk_thermo_list
   implicit none(type, external)
   private

   public :: solvent_vv_thermodynamics_type, new_vv_thermodynamics, vv_chemical_potential, vv_direct_moment

   !> Bulk thermodynamics of a VV solution, atomic units
   type :: solvent_vv_thermodynamics_type
      !> Excess chemical potential per molecule, Hartree
      real(wp) :: chemical_potential = 0.0_wp
      !> Pressure, Hartree/bohr^3
      real(wp) :: pressure = 0.0_wp
      !> Ideal-gas pressure rho kT, Hartree/bohr^3
      real(wp) :: pressure_ideal = 0.0_wp
      !> Isothermal compressibility, bohr^3/Hartree
      real(wp) :: compressibility = 0.0_wp
      !> Direct-correlation moment C = sum_ab c_ab(k = 0), bohr^3
      real(wp) :: direct_moment = 0.0_wp
   end type solvent_vv_thermodynamics_type

contains

   !> Thermodynamics of a VV solution
   !>
   !> @param[out] self     thermodynamics
   !> @param[in]  closure  closure of the solve
   !> @param[in]  rgrid    real-space radial grid
   !> @param[in]  lut      packed-index lookup for the sites
   !> @param[in]  hr       physical h(r), packed (npts, npair)
   !> @param[in]  cr       physical c(r), packed (npts, npair)
   !> @param[in]  rho      molecular number density, 1/bohr^3
   !> @param[in]  beta     inverse thermal energy, 1/Hartree
   !> @param[out] error    shape mismatch, invalid correlations or a vanishing 1 - rho C
   !> @param[out] items    optional owned forward cache of the evaluated bulk list
   subroutine new_vv_thermodynamics(self, closure, rgrid, lut, hr, cr, rho, beta, error, items)
      !> Thermodynamics
      type(solvent_vv_thermodynamics_type), intent(out) :: self
      !> Closure
      class(moist_model_moz_closure_type), intent(in) :: closure
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Packed-index lookup
      type(packed_index_lut), intent(in) :: lut
      !> Physical correlations
      real(wp), intent(in) :: hr(:, :), cr(:, :)
      !> Molecular density and inverse thermal energy
      real(wp), intent(in) :: rho, beta
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional owned forward cache of the evaluated bulk list
      type(thermo_list_type), optional, intent(out), target :: items
      !> Owned cache when no output list is requested
      type(thermo_list_type), target :: local_items
      !> Build in the output directly; ifort loses cache state in the optional-output copy
      type(thermo_list_type), pointer :: evaluated

      if (present(items)) then
         evaluated => items
      else
         evaluated => local_items
      end if

      call vv_chemical_potential(closure, rgrid, lut, hr, cr, rho, beta, self%chemical_potential, error)
      if (allocated(error)) return
      call vv_direct_moment(rgrid, lut, cr, self%direct_moment, error)
      if (allocated(error)) return
      call new_bulk_thermo_list(evaluated, lut%ns, rho, beta, self%direct_moment, error)
      if (.not. allocated(error)) then
         call evaluated%add_input("chemical_potential", self%chemical_potential, error, unit="Hartree")
      end if
      if (.not. allocated(error)) call evaluated%set_base_energy("chemical_potential", error)
      if (.not. allocated(error)) call evaluated%evaluate(error)
      if (.not. allocated(error)) call evaluated%get_value("compressibility", self%compressibility, error)
      if (.not. allocated(error)) call evaluated%get_value("pressure", self%pressure, error)
      if (.not. allocated(error)) call evaluated%get_value("pressure_ideal", self%pressure_ideal, error)
      if (allocated(error)) return
   end subroutine new_vv_thermodynamics

   !> Excess chemical potential mu = (rho/beta) sum_ab int Phi_ab(h, c) d^3r
   !>
   !> - Also valid for a stage iterate; then mu of the partially coupled fluid
   !>
   !> @param[in]  closure  closure of the solve
   !> @param[in]  rgrid    real-space radial grid
   !> @param[in]  lut      packed-index lookup for the sites
   !> @param[in]  hr       physical h(r), packed (npts, npair)
   !> @param[in]  cr       physical c(r), packed (npts, npair)
   !> @param[in]  rho      molecular number density, 1/bohr^3
   !> @param[in]  beta     inverse thermal energy, 1/Hartree
   !> @param[out] mu       excess chemical potential per molecule, Hartree
   !> @param[out] error    shape mismatch, invalid correlations or quadrature failure
   subroutine vv_chemical_potential(closure, rgrid, lut, hr, cr, rho, beta, mu, error)
      !> Closure
      class(moist_model_moz_closure_type), intent(in) :: closure
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Packed-index lookup
      type(packed_index_lut), intent(in) :: lut
      !> Physical correlations
      real(wp), intent(in) :: hr(:, :), cr(:, :)
      !> Molecular density and inverse thermal energy
      real(wp), intent(in) :: rho, beta
      !> Excess chemical potential, Hartree
      real(wp), intent(out) :: mu
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: phi(:, :)
      real(wp) :: integral
      integer :: ip

      mu = 0.0_wp
      if (any(shape(hr) /= shape(cr)) .or. size(hr, 1) /= rgrid%npts .or. &
         & size(hr, 2) /= lut%ns*(lut%ns + 1)/2) then
         call fatal_error(error, "1D VV thermodynamics needs packed (npts, npair) correlations")
         return
      end if
      allocate (phi, mold=hr)
      call closure%free_energy_integrands(hr, cr, phi, error)
      if (allocated(error)) return
      do ip = 1, size(phi, 2)
         call rgrid%integrate_field(phi(:, ip), integral, error)
         if (allocated(error)) return
         mu = mu + pair_weight(lut, ip)*integral
      end do
      mu = rho*mu/beta
   end subroutine vv_chemical_potential

   !> Direct-correlation moment C = sum_ab c_ab(k = 0) = sum_ab int c_ab d^3r
   !>
   !> @param[in]  rgrid   real-space radial grid
   !> @param[in]  lut     packed-index lookup for the sites
   !> @param[in]  cr      physical c(r), packed (npts, npair)
   !> @param[out] moment  C, bohr^3
   !> @param[out] error   shape mismatch
   subroutine vv_direct_moment(rgrid, lut, cr, moment, error)
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Packed-index lookup
      type(packed_index_lut), intent(in) :: lut
      !> Physical direct correlation
      real(wp), intent(in) :: cr(:, :)
      !> Moment
      real(wp), intent(out) :: moment
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: integral
      integer :: ip

      moment = 0.0_wp
      if (size(cr, 2) /= lut%ns*(lut%ns + 1)/2) then
         call fatal_error(error, "1D VV direct-correlation moment needs ns*(ns+1)/2 site pairs")
         return
      end if
      do ip = 1, size(cr, 2)
         call rgrid%integrate_field(cr(:, ip), integral, error)
         if (allocated(error)) return
         moment = moment + pair_weight(lut, ip)*integral
      end do
   end subroutine vv_direct_moment

   !> Ordered site pairs of a packed pair: 1 on the diagonal, 2 off it
   !>
   !> @param[in] lut  packed-index lookup
   !> @param[in] ip   packed pair
   pure real(wp) function pair_weight(lut, ip) result(weight)
      !> Packed-index lookup
      type(packed_index_lut), intent(in) :: lut
      !> Packed pair
      integer, intent(in) :: ip
      integer :: a

      weight = 2.0_wp
      do a = 1, lut%ns
         if (lut%idx(a, a) == ip) weight = 1.0_wp
      end do
   end function pair_weight

end module moist_model_moz_solvent_vv_thermodynamics
