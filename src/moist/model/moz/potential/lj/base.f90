!> Lennard-Jones 12-6 pair data of the MOZ potential layer
!>
!> Plain per-atom data of one side: sigma (bohr), epsilon (Hartree) and a
!> coverage mask. The mixing rule belongs to the interaction, not to a side:
!> the kernels mix the parameters of the two atoms of a pair once with
!> `lj_12_6_mix` and evaluate the potential with `lj_12_6_evaluate`
!>
!> Every LJ term (custom, element, typed) fills this one type, so any two of
!> them combine; the potential set merges partial terms by coverage
module moist_model_moz_potential_lj_base
   use mctc_env, only: wp, error_type, fatal_error
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none(type, external)
   private

   public :: lj_12_6_type, new_lj_12_6, lj_12_6_mix, lj_12_6_evaluate, check_lj_mixing
   public :: lj_mixing_lorentz_berthelot, lj_mixing_geometric, lj_label_len

   !> Length of a per-atom "Type/Src" label
   integer, parameter :: lj_label_len = 16

   !> Lorentz-Berthelot rule: sigma_ij = (sigma_i + sigma_j)/2, eps_ij = sqrt(eps_i eps_j)
   integer, parameter :: lj_mixing_lorentz_berthelot = 1
   !> Geometric rule: sigma_ij = sqrt(sigma_i sigma_j), eps_ij = sqrt(eps_i eps_j)
   integer, parameter :: lj_mixing_geometric = 2

   !> 12-6 Lennard-Jones parameters of one side
   type :: lj_12_6_type
      !> Sigma per atom, bohr
      real(wp), allocatable :: sigma(:)
      !> Epsilon per atom, Hartree
      real(wp), allocatable :: epsilon(:)
      !> Whether each atom has parameters
      logical, allocatable :: covered(:)
      !> "Type/Src" label per atom, e.g. "c3/GAFF" or "Fe/UFF"; blank when uncovered
      character(len=lj_label_len), allocatable :: label(:)
   end type lj_12_6_type

contains

   !> Validate and store per-atom parameters
   !>
   !> Parameters of uncovered atoms are ignored; a term typed from the
   !> structure passes the atoms it could type as `covered`. A covered atom
   !> needs finite parameters, epsilon >= 0 and, where epsilon > 0, sigma > 0;
   !> epsilon = 0 (a polar hydrogen) switches the atom's pair potential off
   !>
   !> @param[out] self Parameters
   !> @param[in] sigma Sigma per atom, bohr
   !> @param[in] epsilon Epsilon per atom, Hartree
   !> @param[out] error Size mismatch or invalid covered parameter
   !> @param[in] covered Optional coverage mask; every atom when absent
   !> @param[in] label Optional "Type/Src" label per atom, kept for covered atoms
   subroutine new_lj_12_6(self, sigma, epsilon, error, covered, label)
      !> Parameters
      type(lj_12_6_type), intent(out) :: self
      !> Sigma per atom, bohr
      real(wp), intent(in) :: sigma(:)
      !> Epsilon per atom, Hartree
      real(wp), intent(in) :: epsilon(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional coverage mask; every atom when absent
      logical, intent(in), optional :: covered(:)
      !> Optional "Type/Src" label per atom; blank when absent
      character(len=*), intent(in), optional :: label(:)
      logical, allocatable :: mask(:)
      integer :: iat
      if (size(sigma) /= size(epsilon)) then
         call fatal_error(error, "Lennard-Jones sigma and epsilon arrays differ in length")
         return
      end if
      if (size(sigma) < 1) then
         call fatal_error(error, "Lennard-Jones parameters cover no atom")
         return
      end if
      if (present(covered)) then
         if (size(covered) /= size(sigma)) then
            call fatal_error(error, "Lennard-Jones coverage mask does not match the parameter arrays")
            return
         end if
         mask = covered
      else
         allocate (mask(size(sigma)), source=.true.)
      end if
      if (any(mask .and. .not. (ieee_is_finite(sigma) .and. ieee_is_finite(epsilon)))) then
         call fatal_error(error, "Lennard-Jones sigma and epsilon must be finite for every covered atom")
         return
      end if
      if (any(mask .and. epsilon < 0.0_wp)) then
         call fatal_error(error, "Lennard-Jones epsilon must not be negative for any covered atom")
         return
      end if
      if (any(mask .and. epsilon > 0.0_wp .and. sigma <= 0.0_wp)) then
         call fatal_error(error, "Lennard-Jones sigma must be positive for every covered atom with epsilon > 0")
         return
      end if
      if (present(label)) then
         if (size(label) /= size(sigma)) then
            call fatal_error(error, "Lennard-Jones labels do not match the parameter arrays")
            return
         end if
      end if
      self%sigma = sigma
      self%epsilon = epsilon
      allocate (self%label(size(sigma)))
      do iat = 1, size(sigma)
         self%label(iat) = ""
         if (present(label) .and. mask(iat)) self%label(iat) = label(iat)
      end do
      call move_alloc(mask, self%covered)
   end subroutine new_lj_12_6

   !> Refuse an unknown mixing rule
   !>
   !> @param[in] mixing Mixing rule, `lj_mixing_lorentz_berthelot` or `lj_mixing_geometric`
   !> @param[out] error Unknown rule
   subroutine check_lj_mixing(mixing, error)
      !> Mixing rule
      integer, intent(in) :: mixing
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      character(len=16) :: label
      select case (mixing)
      case (lj_mixing_lorentz_berthelot, lj_mixing_geometric)
      case default
         write (label, "(i0)") mixing
         call fatal_error(error, "Unknown Lennard-Jones mixing rule "//trim(label)// &
            & "; use lj_mixing_lorentz_berthelot or lj_mixing_geometric")
      end select
   end subroutine check_lj_mixing

   !> Cross parameters of one atom pair
   !>
   !> A pair with an inactive atom (epsilon = 0, whose sigma is not
   !> validated) is inactive: zero cross parameters under every rule. An
   !> unknown rule yields NaN parameters, never a silent zero; callers
   !> validate the rule with `check_lj_mixing`
   !>
   !> @param[in] mixing Mixing rule
   !> @param[in] sigma_i Sigma of the first atom, bohr
   !> @param[in] eps_i Epsilon of the first atom, Hartree
   !> @param[in] sigma_j Sigma of the second atom, bohr
   !> @param[in] eps_j Epsilon of the second atom, Hartree
   !> @param[out] sigma_ij Cross sigma, bohr
   !> @param[out] eps_ij Cross epsilon, Hartree
   pure subroutine lj_12_6_mix(mixing, sigma_i, eps_i, sigma_j, eps_j, sigma_ij, eps_ij)
      !> Mixing rule
      integer, intent(in) :: mixing
      !> Sigma of the first atom, bohr
      real(wp), intent(in) :: sigma_i
      !> Epsilon of the first atom, Hartree
      real(wp), intent(in) :: eps_i
      !> Sigma of the second atom, bohr
      real(wp), intent(in) :: sigma_j
      !> Epsilon of the second atom, Hartree
      real(wp), intent(in) :: eps_j
      !> Cross sigma, bohr
      real(wp), intent(out) :: sigma_ij
      !> Cross epsilon, Hartree
      real(wp), intent(out) :: eps_ij
      if (eps_i == 0.0_wp .or. eps_j == 0.0_wp) then
         sigma_ij = 0.0_wp
         eps_ij = 0.0_wp
         return
      end if
      select case (mixing)
      case (lj_mixing_lorentz_berthelot)
         sigma_ij = 0.5_wp*(sigma_i + sigma_j)
         eps_ij = sqrt(eps_i*eps_j)
      case (lj_mixing_geometric)
         sigma_ij = sqrt(sigma_i*sigma_j)
         eps_ij = sqrt(eps_i*eps_j)
      case default
         sigma_ij = ieee_value(sigma_ij, ieee_quiet_nan)
         eps_ij = ieee_value(eps_ij, ieee_quiet_nan)
      end select
   end subroutine lj_12_6_mix

   !> 12-6 potential of one mixed pair over a distance vector
   !>
   !> sr6 = (sigma/r)^6, u = 4 eps (sr6^2 - sr6), du/dr = -24 eps (2 sr6^2 - sr6)/r.
   !> A distance of zero is a caller error; the kernels guard it
   !>
   !> @param[in] sigma_ij Cross sigma, bohr
   !> @param[in] eps_ij Cross epsilon, Hartree
   !> @param[in] r Distances, bohr (n)
   !> @param[out] u Pair potential, Hartree (n)
   !> @param[out] du_dr Optional radial derivative, Hartree/bohr (n)
   pure subroutine lj_12_6_evaluate(sigma_ij, eps_ij, r, u, du_dr)
      !> Cross sigma, bohr
      real(wp), intent(in) :: sigma_ij
      !> Cross epsilon, Hartree
      real(wp), intent(in) :: eps_ij
      !> Distances, bohr
      real(wp), intent(in) :: r(:)
      !> Pair potential, Hartree
      real(wp), intent(out) :: u(:)
      !> Optional radial derivative, Hartree/bohr
      real(wp), intent(out), optional :: du_dr(:)
      real(wp) :: sr6
      integer :: ir
      if (present(du_dr)) then
         do ir = 1, size(r)
            sr6 = (sigma_ij/r(ir))**6
            u(ir) = 4.0_wp*eps_ij*(sr6*sr6 - sr6)
            du_dr(ir) = -24.0_wp*eps_ij*(2.0_wp*sr6*sr6 - sr6)/r(ir)
         end do
      else
         do ir = 1, size(r)
            sr6 = (sigma_ij/r(ir))**6
            u(ir) = 4.0_wp*eps_ij*(sr6*sr6 - sr6)
         end do
      end if
   end subroutine lj_12_6_evaluate

end module moist_model_moz_potential_lj_base
