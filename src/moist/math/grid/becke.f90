! This file is part of moist.
!> Becke fuzzy-cell partitioning weights for atom-centered molecular grids.
!>
!> A. D. Becke, "A multicenter numerical integration scheme for polyatomic
!> molecules", J. Chem. Phys. 88, 2547 (1988).
!>
!> For a point `r` and a set of atoms A_1, ..., A_n, the Becke partition
!> assigns non-negative weights `w_A(r)` that sum to unity:
!>    sum_A w_A(r) = 1  for all r (provided atoms are distinct).
!> Each atomic grid contribution is multiplied by `w_A(r)` to avoid
!> double counting in overlapping atomic spheres. Size-dependence is
!> handled via the covalent-radius ratio (chi) adjustment from Becke's
!> Eq. (A4)-(A6).
!>
!> Stiffness: the three-fold iteration of the cutoff polynomial
!>    p(mu) = 1.5*mu - 0.5*mu**3
!> (Becke k=3) is used, matching the original recommendation.
module moist_math_grid_becke
   use mctc_env, only: wp
   use moist_data_atomicrad, only: covalent_rad
   implicit none
   private

   public :: becke_weights

contains

   !> Compute Becke partition weights for a single sample point.
   !>
   !> Uses the size-adjusted (covalent-radius ratio) variant of Becke's
   !> smoothed Voronoi construction with `k = 3` iterations of the
   !> cutoff polynomial.
   !>
   !> @param[in]  point     Sample point in bohr, shape (3).
   !> @param[in]  nat       Number of atoms.
   !> @param[in]  xyz       Atom positions in bohr, shape (3, nat).
   !> @param[in]  numbers   Atomic numbers, shape (nat).
   !> @param[out] weights   Per-atom partition weights, shape (nat);
   !>                       sum(weights) = 1 for distinct atoms.
   pure subroutine becke_weights(point, nat, xyz, numbers, weights)
      !> Sample point in bohr
      real(wp), intent(in)  :: point(3)
      !> Number of atoms
      integer, intent(in)  :: nat
      !> Atom positions in bohr
      real(wp), intent(in)  :: xyz(3, nat)
      !> Atomic numbers
      integer, intent(in)  :: numbers(nat)
      !> Per-atom partition weights, summing to 1
      real(wp), intent(out) :: weights(nat)

      integer  :: ii, jj
      real(wp) :: ri, rj, rij, chi, mu_prime, a_adj, mu, nu, s, total

      weights = 1.0_wp
      do ii = 1, nat
         ri = norm2(point - xyz(:, ii))
         do jj = 1, nat
            if (jj == ii) cycle
            rj = norm2(point - xyz(:, jj))
            rij = norm2(xyz(:, ii) - xyz(:, jj))
            if (rij <= 0.0_wp) cycle
            mu = (ri - rj)/rij
            ! Size adjustment based on covalent-radius ratio (Becke Eq. A6)
            chi = covalent_rad(numbers(ii))/covalent_rad(numbers(jj))
            mu_prime = (chi - 1.0_wp)/(chi + 1.0_wp)
            a_adj = mu_prime/(mu_prime*mu_prime - 1.0_wp)
            if (a_adj > 0.5_wp) a_adj = 0.5_wp
            if (a_adj < -0.5_wp) a_adj = -0.5_wp
            nu = mu + a_adj*(1.0_wp - mu*mu)
            s = 0.5_wp*(1.0_wp - becke_k3(nu))
            weights(ii) = weights(ii)*s
         end do
      end do
      total = sum(weights)
      if (total > 0.0_wp) then
         weights = weights/total
      end if
   end subroutine becke_weights

   !> Three-fold iterated Becke cutoff polynomial p(x) = 1.5*x - 0.5*x**3.
   !>
   !> Hardcoded for k = 3 (Becke's recommendation). Avoids the overhead
   !> and pitfalls of a recursive implementation.
   pure function becke_k3(x) result(y)
      !> Input value in [-1, 1]
      real(wp), intent(in) :: x
      !> Three-fold iterate of p
      real(wp) :: y

      real(wp) :: y1, y2

      y1 = 1.5_wp*x - 0.5_wp*x*x*x
      y2 = 1.5_wp*y1 - 0.5_wp*y1*y1*y1
      y = 1.5_wp*y2 - 0.5_wp*y2*y2*y2
   end function becke_k3

end module moist_math_grid_becke
