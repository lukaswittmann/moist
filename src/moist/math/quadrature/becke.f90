!> Becke fuzzy-cell partitioning weights for atom-centered molecular grids
module moist_math_quadrature_becke
   use mctc_env, only: wp
   use moist_data_atomicrad, only: covalent_rad
   implicit none(type, external)
   private

   public :: becke_weights

   !> Default number of iterations of the Becke cutoff polynomial
   integer, parameter :: default_stiffness = 3

contains

   !> Compute Becke partition weights for a single sample point
   !>
   !> Uses the size-adjusted (covalent-radius ratio) variant of Becke's
   !> smoothed Voronoi construction; the cell (cutoff) function defaults to
   !> `k = 3` iterations of `p(mu) = 1.5*mu - 0.5*mu**3`; `stiffness` selects a
   !> different iteration count and `ssf_a` switches to the compactly supported
   !> Stratmann-Scuseria-Frisch cell function instead (it takes precedence over
   !> `stiffness` when both are given)
   !>
   !> @param[in]  point      Sample point in bohr, shape (3)
   !> @param[in]  nat        Number of atoms
   !> @param[in]  xyz        Atom positions in bohr, shape (3, nat)
   !> @param[in]  numbers    Atomic numbers, shape (nat)
   !> @param[out] weights    Per-atom partition weights, shape (nat);
   !>                        sum(weights) = 1 for distinct atoms
   !> @param[in]  stiffness  Optional iteration count `k` of the Becke cutoff
   !>                        polynomial (default 3); larger is stiffer
   !> @param[in]  ssf_a      Optional SSF cutoff parameter `a` in (0, 1];
   !>                        when present the SSF cell function replaces the
   !>                        iterated Becke polynomial
   pure subroutine becke_weights(point, nat, xyz, numbers, weights, stiffness, ssf_a)
      !> Sample point in bohr
      real(wp), intent(in)  :: point(3)
      !> Number of atoms
      integer,  intent(in)  :: nat
      !> Atom positions in bohr
      real(wp), intent(in)  :: xyz(3, nat)
      !> Atomic numbers
      integer,  intent(in)  :: numbers(nat)
      !> Per-atom partition weights, summing to 1
      real(wp), intent(out) :: weights(nat)
      !> Iterations of the Becke cutoff polynomial (default 3)
      integer,  intent(in), optional :: stiffness
      !> SSF cutoff parameter; selects the SSF cell function when present
      real(wp), intent(in), optional :: ssf_a

      integer  :: ii, jj, kiter
      real(wp) :: ri, rj, rij, chi, mu_prime, a_adj, mu, nu, s, total
      logical  :: use_ssf

      use_ssf = present(ssf_a)
      kiter = default_stiffness
      if (present(stiffness)) kiter = stiffness

      weights = 1.0_wp
      do ii = 1, nat
         ri = norm2(point - xyz(:, ii))
         do jj = 1, nat
            if (jj == ii) cycle
            rj = norm2(point - xyz(:, jj))
            rij = norm2(xyz(:, ii) - xyz(:, jj))
            if (rij <= 0.0_wp) cycle
            mu = (ri - rj) / rij
            ! Size adjustment based on covalent-radius ratio (Becke Eq. A6)
            chi = covalent_rad(numbers(ii)) / covalent_rad(numbers(jj))
            mu_prime = (chi - 1.0_wp) / (chi + 1.0_wp)
            a_adj = mu_prime / (mu_prime * mu_prime - 1.0_wp)
            if (a_adj >  0.5_wp) a_adj =  0.5_wp
            if (a_adj < -0.5_wp) a_adj = -0.5_wp
            nu = mu + a_adj * (1.0_wp - mu * mu)
            if (use_ssf) then
               s = ssf_cell(nu, ssf_a)
            else if (kiter == default_stiffness) then
               s = 0.5_wp * (1.0_wp - becke_k3(nu))
            else
               s = 0.5_wp * (1.0_wp - becke_iterate(nu, kiter))
            end if
            weights(ii) = weights(ii) * s
         end do
      end do
      total = sum(weights)
      if (total > 0.0_wp) then
         weights = weights / total
      end if
   end subroutine becke_weights

   !> Three-fold iterated Becke cutoff polynomial p(x) = 1.5*x - 0.5*x**3
   !>
   !> Hardcoded for k = 3 (Becke's recommendation); avoids the overhead
   !> and pitfalls of a recursive implementation
   !>
   !> @param[in] x  Input value in [-1, 1]
   pure function becke_k3(x) result(y)
      !> Input value in [-1, 1]
      real(wp), intent(in) :: x
      !> Three-fold iterate of p
      real(wp) :: y

      real(wp) :: y1, y2

      y1 = 1.5_wp * x  - 0.5_wp * x * x * x
      y2 = 1.5_wp * y1 - 0.5_wp * y1 * y1 * y1
      y  = 1.5_wp * y2 - 0.5_wp * y2 * y2 * y2
   end function becke_k3

   !> `k`-fold iterated Becke cutoff polynomial p(x) = 1.5*x - 0.5*x**3
   !>
   !> Generic counterpart of [[becke_k3]] for a caller-chosen stiffness; a
   !> non-positive `k` returns the argument unchanged (no smoothing)
   !>
   !> @param[in] x  Input value in [-1, 1]
   !> @param[in] k  Number of iterations of p
   pure function becke_iterate(x, k) result(y)
      !> Input value in [-1, 1]
      real(wp), intent(in) :: x
      !> Number of iterations
      integer,  intent(in) :: k
      !> `k`-fold iterate of p
      real(wp) :: y

      integer :: i

      y = x
      do i = 1, k
         y = 1.5_wp * y - 0.5_wp * y * y * y
      end do
   end function becke_iterate

   !> Stratmann-Scuseria-Frisch compactly supported cell function
   !>
   !> R. E. Stratmann, G. E. Scuseria, M. J. Frisch, Chem. Phys. Lett. 257,
   !> 213 (1996)
   !>
   !> With `z = nu/a` the switching polynomial is
   !>    g(z) = (35 z - 35 z**3 + 21 z**5 - 5 z**7) / 16,
   !> C^3-continuous at `z = +/-1`, so `s = 0.5*(1 - g)` is exactly 1 for
   !> `nu < -a` and exactly 0 for `nu > a`
   !>
   !> The hard zero is what makes SSF cheap: a point outside every
   !> neighbour's `a`-window is owned outright, and points another atom
   !> owns get weight identically zero
   !>
   !> @param[in] nu  Size-adjusted elliptic coordinate
   !> @param[in] a   Cutoff parameter in (0, 1]; 0.64 is the SSF recommendation
   pure function ssf_cell(nu, a) result(s)
      !> Size-adjusted elliptic coordinate
      real(wp), intent(in) :: nu
      !> Cutoff parameter
      real(wp), intent(in) :: a
      !> Cell function value in [0, 1]
      real(wp) :: s

      real(wp) :: z, z2, g

      if (a <= 0.0_wp) then
         s = 0.5_wp
         if (nu < 0.0_wp) s = 1.0_wp
         if (nu > 0.0_wp) s = 0.0_wp
         return
      end if

      z = nu / a
      if (z <= -1.0_wp) then
         s = 1.0_wp
      else if (z >= 1.0_wp) then
         s = 0.0_wp
      else
         z2 = z * z
         g = (35.0_wp * z * (1.0_wp + z2 * (-1.0_wp + z2 * (0.6_wp &
            & - z2 * (1.0_wp / 7.0_wp))))) / 16.0_wp
         s = 0.5_wp * (1.0_wp - g)
      end if
   end function ssf_cell

end module moist_math_quadrature_becke
