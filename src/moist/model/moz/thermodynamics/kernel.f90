!> Scalar thermodynamic kernels shared by the MOZ models
!>
!> - Atomic units: rho in bohr^-3, beta in Hartree^-1,
!>   pressure in Hartree/bohr^3, compressibility in bohr^3/Hartree
!> - Bulk C = sum_ab int c_ab d^3r in bohr^3
!> - UV D = sum_ia rho_a int c_ia d^3r, dimensionless
!> - Grid reductions and input checks belong to the callers
!> - Pressure is the existing 3D-RISM/HNC expression, also used for KH/PSE
!> - Response shifts are beta dE/dD, dimensionless, at fixed bulk state
!> - No grid, model, closure or correction-selector dependencies
module moist_model_moz_thermodynamics_kernel
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: bulk_compressibility_denominator, bulk_compressibility, bulk_pressure, bulk_ideal_pressure
   public :: partial_molar_volume, pc_correction, pc_plus_correction
   public :: pc_moment_shift, pc_plus_moment_shift
   public :: denominator_adjoint, compressibility_adjoint, pressure_adjoint, ideal_pressure_adjoint
   public :: volume_adjoint, pc_adjoint, pc_plus_adjoint, moment_response_shift
   public :: scale_energy_derivative
   public :: coulomb_box_mean, electrostatic_finite_size_correction

contains

   !> Average 1/r over a centred orthorhombic box, bohr^-1
   !>
   !> - Integrated analytically over an octant
   !> - Cube: mean = [6 asinh(1/sqrt(2)) - pi/2]/L = 2.380077363979553/L
   !> - Positive lengths are checked by the caller
   !> - Lengths are normalized by their geometric mean to avoid squared units
   !>
   !> @param[in] lengths  box edge lengths, (3), bohr
   pure function coulomb_box_mean(lengths) result(mean)
      !> Box edge lengths, bohr
      real(wp), intent(in) :: lengths(3)
      !> Mean 1/r, bohr^-1, and work values: scaled box, scale, radius, octant integral
      real(wp) :: mean, scale, box(3), radius, integral
      integer :: i, j, k

      scale = exp(sum(log(lengths))/3.0_wp)
      box = 0.5_wp*(lengths/scale)
      radius = norm2(box)
      integral = 0.0_wp
      do i = 1, 3
         j = mod(i, 3) + 1
         k = mod(i + 1, 3) + 1
         integral = integral + box(j)*box(k)*asinh(box(i)/norm2(box([j, k]))) - &
            & 0.5_wp*box(i)**2*atan(box(j)*box(k)/(box(i)*radius))
      end do
      mean = (integral/product(box))/scale
   end function coulomb_box_mean

   !> Leading missing-zero-mode correction for the UV energy using physical c(r), Hartree
   !>
   !> - delta h_a = -4 pi beta Q z_a(0)/V; dPhi/dh_a ~ beta q_a Q/(2r)
   !> - Hence delta G = -Q^2 S <1/r>_box / 2, with S = 1 - 1/epsilon
   !> - Physical c already restores the isolated tail
   !> - The periodic Ewald Madelung coefficient cannot be added again
   !> - Follows from the background response of Kovalenko and Hirata,
   !>   JCP 112, 10391 (2000), Eq. (20), applied to this energy
   !> - Cavity/multipole and higher finite-volume terms are outside the leading approximation
   !>
   !> @param[in] charge        total solute charge Q, elementary charges
   !> @param[in] screening     bulk screening S = 1 - 1/epsilon, dimensionless
   !> @param[in] coulomb_mean  box average of 1/r, bohr^-1
   pure elemental function electrostatic_finite_size_correction(charge, screening, coulomb_mean) result(energy)
      !> Total solute charge Q, elementary charges; screening S, dimensionless; box average of 1/r, bohr^-1
      real(wp), intent(in) :: charge, screening, coulomb_mean
      !> Correction energy, Hartree
      real(wp) :: energy
      energy = -0.5_wp*charge**2*screening*coulomb_mean
   end function electrostatic_finite_size_correction

   !> Dimensionless denominator 1 - rho C of the bulk compressibility
   !>
   !> @param[in] rho     molecular number density, bohr^-3
   !> @param[in] moment  bulk C = sum_ab int c_ab d^3r, bohr^3
   pure elemental function bulk_compressibility_denominator(rho, moment) result(denominator)
      !> Molecular number density, bohr^-3; bulk C, bohr^3
      real(wp), intent(in) :: rho, moment
      !> Denominator 1 - rho C, dimensionless
      real(wp) :: denominator

      denominator = 1.0_wp - rho*moment
   end function bulk_compressibility_denominator

   !> Isothermal compressibility beta/[rho (1 - rho C)], bohr^3/Hartree
   !>
   !> @param[in] rho          molecular number density, bohr^-3
   !> @param[in] beta         inverse thermal energy, Hartree^-1
   !> @param[in] denominator  precomputed, nonzero 1 - rho C, dimensionless
   pure elemental function bulk_compressibility(rho, beta, denominator) result(kappa)
      !> Molecular number density, bohr^-3; inverse thermal energy, Hartree^-1; nonzero 1 - rho C, dimensionless
      real(wp), intent(in) :: rho, beta, denominator
      !> Isothermal compressibility, bohr^3/Hartree
      real(wp) :: kappa

      kappa = beta/(rho*denominator)
   end function bulk_compressibility

   !> Existing 3D-RISM pressure [(n+1) rho/2 - rho^2 C/2]/beta, Hartree/bohr^3
   !>
   !> @param[in] nsite   sites per solvent molecule
   !> @param[in] rho     molecular number density, bohr^-3
   !> @param[in] beta    inverse thermal energy, Hartree^-1
   !> @param[in] moment  bulk C = sum_ab int c_ab d^3r, bohr^3
   pure elemental function bulk_pressure(nsite, rho, beta, moment) result(pressure)
      !> Sites per solvent molecule
      integer, intent(in) :: nsite
      !> Molecular number density, bohr^-3; inverse thermal energy, Hartree^-1; bulk C, bohr^3
      real(wp), intent(in) :: rho, beta, moment
      !> Bulk pressure, Hartree/bohr^3
      real(wp) :: pressure

      pressure = (0.5_wp*real(nsite + 1, wp)*rho - 0.5_wp*rho**2*moment)/beta
   end function bulk_pressure

   !> Molecular ideal-gas pressure rho/beta, Hartree/bohr^3
   !>
   !> @param[in] rho   molecular number density, bohr^-3
   !> @param[in] beta  inverse thermal energy, Hartree^-1
   pure elemental function bulk_ideal_pressure(rho, beta) result(pressure)
      !> Molecular number density, bohr^-3; inverse thermal energy, Hartree^-1
      real(wp), intent(in) :: rho, beta
      !> Ideal-gas pressure, Hartree/bohr^3
      real(wp) :: pressure

      pressure = rho/beta
   end function bulk_ideal_pressure

   !> Partial molar volume kappa/beta (1 - D), bohr^3
   !>
   !> @param[in] moment  UV D = sum_ia rho_a int c_ia d^3r, dimensionless
   !> @param[in] kappa   bulk isothermal compressibility, bohr^3/Hartree
   !> @param[in] beta    inverse thermal energy, Hartree^-1
   pure elemental function partial_molar_volume(moment, kappa, beta) result(volume)
      !> UV D, dimensionless; bulk compressibility, bohr^3/Hartree; inverse thermal energy, Hartree^-1
      real(wp), intent(in) :: moment, kappa, beta
      !> Partial molar volume, bohr^3
      real(wp) :: volume

      volume = kappa/beta*(1.0_wp - moment)
   end function partial_molar_volume

   !> PC energy -P Vbar, Hartree
   !>
   !> @param[in] pressure  bulk RISM pressure, Hartree/bohr^3
   !> @param[in] volume    partial molar volume, bohr^3
   pure elemental function pc_correction(pressure, volume) result(energy)
      !> Bulk pressure, Hartree/bohr^3; partial molar volume, bohr^3
      real(wp), intent(in) :: pressure, volume
      !> PC energy, Hartree
      real(wp) :: energy

      energy = -pressure*volume
   end function pc_correction

   !> PC+ energy -(P - P_ideal) Vbar, Hartree
   !>
   !> @param[in] pressure        bulk RISM pressure, Hartree/bohr^3
   !> @param[in] pressure_ideal  molecular ideal-gas pressure, Hartree/bohr^3
   !> @param[in] volume          partial molar volume, bohr^3
   pure elemental function pc_plus_correction(pressure, pressure_ideal, volume) result(energy)
      !> Bulk and ideal-gas pressure, Hartree/bohr^3; partial molar volume, bohr^3
      real(wp), intent(in) :: pressure, pressure_ideal, volume
      !> PC+ energy, Hartree
      real(wp) :: energy

      energy = -(pressure - pressure_ideal)*volume
   end function pc_plus_correction

   !> PC response shift beta dE_PC/dD = P kappa, dimensionless
   !>
   !> - The UV energy reduction already carries 1/beta
   !> - Its dPhi/dc receives this shift directly, without another beta conversion
   !>
   !> @param[in] pressure  bulk RISM pressure, Hartree/bohr^3
   !> @param[in] kappa     bulk isothermal compressibility, bohr^3/Hartree
   pure elemental function pc_moment_shift(pressure, kappa) result(shift)
      !> Bulk pressure, Hartree/bohr^3; bulk compressibility, bohr^3/Hartree
      real(wp), intent(in) :: pressure, kappa
      !> PC response shift, dimensionless
      real(wp) :: shift

      shift = pressure*kappa
   end function pc_moment_shift

   !> PC+ response shift beta dE_PC+/dD = (P - P_ideal) kappa, dimensionless
   !>
   !> - As for PC, the shift already includes beta and is added once to dPhi/dc
   !>
   !> @param[in] pressure        bulk RISM pressure, Hartree/bohr^3
   !> @param[in] pressure_ideal  molecular ideal-gas pressure, Hartree/bohr^3
   !> @param[in] kappa           bulk isothermal compressibility, bohr^3/Hartree
   pure elemental function pc_plus_moment_shift(pressure, pressure_ideal, kappa) result(shift)
      !> Bulk and ideal-gas pressure, Hartree/bohr^3; bulk compressibility, bohr^3/Hartree
      real(wp), intent(in) :: pressure, pressure_ideal, kappa
      !> PC+ response shift, dimensionless
      real(wp) :: shift

      shift = (pressure - pressure_ideal)*kappa
   end function pc_plus_moment_shift

   !> Reverse kernel of the denominator 1 - rho C
   !>
   !> - Reverse kernels add ordinary Hartree-energy sensitivities, without beta scaling
   !> - Adds to adjoints(1:2) for rho and C
   !>
   !> @param[in]    rho       molecular number density, bohr^-3
   !> @param[in]    moment    bulk C, bohr^3
   !> @param[in]    seed      energy sensitivity of the denominator
   !> @param[inout] adjoints  accumulated sensitivities, (2): rho, C
   pure subroutine denominator_adjoint(rho, moment, seed, adjoints)
      !> Molecular number density, bohr^-3; bulk C, bohr^3; energy sensitivity of the denominator
      real(wp), intent(in) :: rho, moment, seed
      !> Accumulated sensitivities, (2): rho, C
      real(wp), intent(inout) :: adjoints(:)
      adjoints(1) = adjoints(1) - seed*moment
      adjoints(2) = adjoints(2) - seed*rho
   end subroutine denominator_adjoint

   !> Reverse kernel of the bulk compressibility
   !>
   !> - Dependencies: rho [bohr^-3], beta [Hartree^-1], dimensionless denominator
   !>
   !> @param[in]    rho          molecular number density, bohr^-3
   !> @param[in]    beta         inverse thermal energy, Hartree^-1
   !> @param[in]    denominator  1 - rho C, dimensionless
   !> @param[in]    kappa        bulk isothermal compressibility, bohr^3/Hartree
   !> @param[in]    seed         energy sensitivity of kappa
   !> @param[inout] adjoints     accumulated sensitivities, (3): rho, beta, denominator
   pure subroutine compressibility_adjoint(rho, beta, denominator, kappa, seed, adjoints)
      !> Density rho, bohr^-3; beta, Hartree^-1; 1 - rho C, dimensionless; kappa, bohr^3/Hartree; energy sensitivity seed
      real(wp), intent(in) :: rho, beta, denominator, kappa, seed
      !> Accumulated sensitivities, (3): rho, beta, denominator
      real(wp), intent(inout) :: adjoints(:)
      adjoints(1) = adjoints(1) - seed*kappa/rho
      adjoints(2) = adjoints(2) + seed*kappa/beta
      adjoints(3) = adjoints(3) - seed*kappa/denominator
   end subroutine compressibility_adjoint

   !> Reverse kernel of the bulk pressure
   !>
   !> - Dependencies: rho [bohr^-3], beta [Hartree^-1], bulk C [bohr^3]
   !>
   !> @param[in]    nsite     sites per solvent molecule
   !> @param[in]    rho       molecular number density, bohr^-3
   !> @param[in]    beta      inverse thermal energy, Hartree^-1
   !> @param[in]    moment    bulk C, bohr^3
   !> @param[in]    pressure  bulk pressure, Hartree/bohr^3
   !> @param[in]    seed      energy sensitivity of the pressure
   !> @param[inout] adjoints  accumulated sensitivities, (3): rho, beta, C
   pure subroutine pressure_adjoint(nsite, rho, beta, moment, pressure, seed, adjoints)
      !> Sites per solvent molecule
      integer, intent(in) :: nsite
      !> Density rho, bohr^-3; beta, Hartree^-1; bulk C, bohr^3; pressure, Hartree/bohr^3; energy sensitivity seed
      real(wp), intent(in) :: rho, beta, moment, pressure, seed
      !> Accumulated sensitivities, (3): rho, beta, C
      real(wp), intent(inout) :: adjoints(:)
      adjoints(1) = adjoints(1) + seed*(0.5_wp*real(nsite + 1, wp) - rho*moment)/beta
      adjoints(2) = adjoints(2) - seed*pressure/beta
      adjoints(3) = adjoints(3) - seed*0.5_wp*rho**2/beta
   end subroutine pressure_adjoint

   !> Reverse kernel of the ideal-gas pressure
   !>
   !> - Dependencies: rho [bohr^-3], beta [Hartree^-1]
   !>
   !> @param[in]    beta      inverse thermal energy, Hartree^-1
   !> @param[in]    pressure  ideal-gas pressure, Hartree/bohr^3
   !> @param[in]    seed      energy sensitivity of the ideal-gas pressure
   !> @param[inout] adjoints  accumulated sensitivities, (2): rho, beta
   pure subroutine ideal_pressure_adjoint(beta, pressure, seed, adjoints)
      !> Inverse thermal energy, Hartree^-1; ideal-gas pressure, Hartree/bohr^3; energy sensitivity seed
      real(wp), intent(in) :: beta, pressure, seed
      !> Accumulated sensitivities, (2): rho, beta
      real(wp), intent(inout) :: adjoints(:)
      adjoints(1) = adjoints(1) + seed/beta
      adjoints(2) = adjoints(2) - seed*pressure/beta
   end subroutine ideal_pressure_adjoint

   !> Reverse kernel of the partial molar volume
   !>
   !> - Dependencies: dimensionless UV D, kappa [bohr^3/Hartree], beta [Hartree^-1]
   !>
   !> @param[in]    moment    UV D, dimensionless
   !> @param[in]    kappa     bulk isothermal compressibility, bohr^3/Hartree
   !> @param[in]    beta      inverse thermal energy, Hartree^-1
   !> @param[in]    volume    partial molar volume, bohr^3
   !> @param[in]    seed      energy sensitivity of the volume
   !> @param[inout] adjoints  accumulated sensitivities, (3): D, kappa, beta
   pure subroutine volume_adjoint(moment, kappa, beta, volume, seed, adjoints)
      !> UV D, dimensionless; kappa, bohr^3/Hartree; beta, Hartree^-1; volume, bohr^3; energy sensitivity seed
      real(wp), intent(in) :: moment, kappa, beta, volume, seed
      !> Accumulated sensitivities, (3): D, kappa, beta
      real(wp), intent(inout) :: adjoints(:)
      adjoints(1) = adjoints(1) - seed*kappa/beta
      adjoints(2) = adjoints(2) + seed*(1.0_wp - moment)/beta
      adjoints(3) = adjoints(3) - seed*volume/beta
   end subroutine volume_adjoint

   !> Reverse kernel of the PC energy
   !>
   !> - Dependencies: pressure [Hartree/bohr^3], Vbar [bohr^3]
   !>
   !> @param[in]    pressure  bulk pressure, Hartree/bohr^3
   !> @param[in]    volume    partial molar volume, bohr^3
   !> @param[in]    seed      energy sensitivity of the PC energy
   !> @param[inout] adjoints  accumulated sensitivities, (2): pressure, Vbar
   pure subroutine pc_adjoint(pressure, volume, seed, adjoints)
      !> Bulk pressure, Hartree/bohr^3; partial molar volume, bohr^3; energy sensitivity seed
      real(wp), intent(in) :: pressure, volume, seed
      !> Accumulated sensitivities, (2): pressure, Vbar
      real(wp), intent(inout) :: adjoints(:)
      adjoints(1) = adjoints(1) - seed*volume
      adjoints(2) = adjoints(2) - seed*pressure
   end subroutine pc_adjoint

   !> Reverse kernel of the PC+ energy
   !>
   !> - Dependencies: pressure, ideal pressure [Hartree/bohr^3], Vbar [bohr^3]
   !>
   !> @param[in]    pressure  bulk pressure, Hartree/bohr^3
   !> @param[in]    ideal     ideal-gas pressure, Hartree/bohr^3
   !> @param[in]    volume    partial molar volume, bohr^3
   !> @param[in]    seed      energy sensitivity of the PC+ energy
   !> @param[inout] adjoints  accumulated sensitivities, (3): pressure, ideal pressure, Vbar
   pure subroutine pc_plus_adjoint(pressure, ideal, volume, seed, adjoints)
      !> Pressures, Hartree/bohr^3; partial molar volume, bohr^3; energy sensitivity seed
      real(wp), intent(in) :: pressure, ideal, volume, seed
      !> Accumulated sensitivities, (3): pressure, ideal pressure, Vbar
      real(wp), intent(inout) :: adjoints(:)
      adjoints(1) = adjoints(1) - seed*volume
      adjoints(2) = adjoints(2) + seed*volume
      adjoints(3) = adjoints(3) - seed*(pressure - ideal)
   end subroutine pc_plus_adjoint

   !> Convert accumulated dE/dD [Hartree] once to the dimensionless UV response shift
   !>
   !> @param[in] beta        inverse thermal energy, Hartree^-1
   !> @param[in] derivative  accumulated dE/dD, Hartree
   pure elemental function moment_response_shift(beta, derivative) result(shift)
      !> Inverse thermal energy, Hartree^-1; accumulated dE/dD, Hartree
      real(wp), intent(in) :: beta, derivative
      !> Response shift beta dE/dD, dimensionless
      real(wp) :: shift
      shift = beta*derivative
   end function moment_response_shift

   !> Apply the accumulated dimensionless sensitivity of the raw energy root
   !>
   !> @param[in] seed        sensitivity of the energy to the raw energy root
   !> @param[in] derivative  accumulated derivative to scale
   pure elemental function scale_energy_derivative(seed, derivative) result(scaled)
      !> Sensitivity of the energy to the raw energy root; accumulated derivative to scale
      real(wp), intent(in) :: seed, derivative
      !> Scaled derivative
      real(wp) :: scaled
      scaled = seed*derivative
   end function scale_energy_derivative

end module moist_model_moz_thermodynamics_kernel
