!> Core fill for ECP nuclei of the isodensity level set
!>
!> An ECP pseudo-density can fall below `rho_iso` at its own nucleus, which
!> gives the level set a spurious inner sheet. Every atom whose nuclear charge
!> `z`, answered by the host through the `nuclear_charges` coupling request, is
!> below its atomic number `Z` receives the compact bump
!>
!>    g(r) = A (1 - |r - R| / X)^5,   |r - R| < X,   A = 100 rho_iso
!>
!> zero with its first four derivatives at `X = min(f R_A, 2 Bohr)` and beyond,
!> with `R_A` the atom's cavity radius and `f` the parameter
!> `core_radius_fraction` (default 1/2); `g >= rho_iso` within `0.60 X`
!>
!> - a cusp at the nucleus, `|grad g| = 5 A / X` there instead of zero, so the
!>   nucleus is no critical point of the level set a projector could settle on;
!>   the higher derivatives grow like `1 / |r - R|^k` towards it, and within
!>   `1e-8 X` of it only the value is added
!> - with every built-in radius model, `rho + g > rho_iso` throughout each ball
!>   for def2, cc-pVnZ-PP and LANL2DZ ECPs (Na to Rn, cations included) at
!>   `rho_iso <= 1e-3`, narrowest for LANL2DZ Na (1.2 `rho_iso`); at `1e-2` only
!>   LANL2DZ Na and Mg trip the check; these data hold for `f = 1/2`
!> - a smaller `f` admits a surface that comes closer to an ECP nucleus
!> - the ball must stay inside the cavity, so the level set carries no explicit
!>   nuclear derivative at any surface point: [[iso_core_fill_check_surface]]
!>   rejects every converged surface point inside a ball, and
!>   [[iso_core_fill_check]] already every evaluated point with `rho + g <= rho_iso`
module moist_cavity_drop_lsf_isodensity_core
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_isodensity_param, only: moist_cavity_drop_lsf_isodensity_param_type
   implicit none(type, external)
   private

   public :: moist_iso_core_fill_type, iso_check_nuclear_charges

   integer, parameter :: ndim = 3

   !> Fill at the nucleus in units of `rho_iso`
   real(wp), parameter :: peak_ratio = 100.0_wp
   !> Largest fill radius in Bohr
   real(wp), parameter :: max_radius = 2.0_wp
   !> Distance from the nucleus, in units of the fill radius, inside which only
   !> the value of the fill is added
   real(wp), parameter :: tip_radius = 1.0e-8_wp

   !> Core fill around the ECP nuclei of one structure
   type :: moist_iso_core_fill_type
      !> Number of filled nuclei
      integer :: n = 0
      !> Atom index of each filled nucleus
      integer, allocatable :: atom(:)
      !> Position of each filled nucleus in Bohr (3, n)
      real(wp), allocatable :: center(:, :)
      !> Fill radius of each filled nucleus in Bohr
      real(wp), allocatable :: radius(:)
      !> Squared fill radius of each filled nucleus in Bohr^2
      real(wp), allocatable :: radius2(:)
      !> Most fill balls that can hold one point
      integer :: overlap = 1
      !> Fill at the nucleus in Bohr^-3
      real(wp) :: amplitude = 0.0_wp
      !> Density isovalue the fill is guarded against
      real(wp) :: rho_iso = 0.0_wp
   contains
      !> Select the ECP nuclei of a structure and size their fill
      procedure :: update => iso_core_fill_update
      !> Fill whose ball holds a point
      procedure :: holder => iso_core_fill_holder
      !> Add the fill and its spatial derivatives at one point
      procedure :: add => iso_core_fill_add
      !> Reject a point inside a fill ball that is not inside the cavity
      procedure :: check => iso_core_fill_check
      !> Reject converged surface points inside a fill ball
      procedure :: check_surface => iso_core_fill_check_surface
      !> Bound on `|grad g| / max(g, rho_iso)` (zero without fill)
      procedure :: log_grad_bound => iso_core_fill_log_grad_bound
   end type moist_iso_core_fill_type

contains

   !> Refuse nuclear charges outside (0, Z] by name
   !>
   !> An ECP removes core electrons, so `0 < z <= Z`; `z < Z` then marks an
   !> ECP nucleus
   !>
   !> @param[in]  mol   Molecular structure supplying `Z`
   !> @param[in]  z     Nuclear charges of the host Hamiltonian (nat)
   !> @param[out] error First charge out of range, e.g. "z(3) = 80.5 exceeds Z = 79"
   subroutine iso_check_nuclear_charges(mol, z, error)
      !> Molecular structure supplying `Z`
      type(structure_type), intent(in) :: mol
      !> Nuclear charges of the host Hamiltonian
      real(wp), intent(in) :: z(:)
      !> First charge out of range
      type(error_type), allocatable, intent(out) :: error

      !> Atom index and its atomic number
      integer :: iat, number
      !> Formatted index, charge and atomic number
      character(len=32) :: sat, sz, snum

      do iat = 1, min(size(z), mol%nat)
         number = mol%num(mol%id(iat))
         if (z(iat) > 0.0_wp .and. z(iat) <= real(number, wp)) cycle
         write (sat, "(i0)") iat
         write (snum, "(i0)") number
         call format_charge(z(iat), sz)
         if (z(iat) > 0.0_wp) then
            call fatal_error(error, "nuclear_charges: z("//trim(sat)//") = "//trim(sz)// &
               & " exceeds Z = "//trim(snum))
         else
            call fatal_error(error, "nuclear_charges: z("//trim(sat)//") = "//trim(sz)// &
               & " is not positive")
         end if
         return
      end do

   end subroutine iso_check_nuclear_charges

   !> Shortest decimal form of a charge: no trailing zeros, no trailing point
   !>
   !> @param[in]  value Charge
   !> @param[out] text  Formatted charge, left adjusted
   subroutine format_charge(value, text)
      !> Charge
      real(wp), intent(in) :: value
      !> Formatted charge
      character(len=*), intent(out) :: text

      !> Last character kept
      integer :: last

      write (text, "(f0.6)") value
      if (index(text, ".") == 0) return
      last = len_trim(text)
      do while (last > 1 .and. text(last:last) == "0")
         last = last - 1
      end do
      if (text(last:last) == ".") last = last - 1
      if (last < 1) then
         text = "0"
         return
      end if
      text = text(:last)
      if (text(1:1) == ".") text = "0"//text
      if (text(1:2) == "-.") text = "-0"//text(2:)

   end subroutine format_charge

   !> Select the ECP nuclei of a structure and size their fill
   !>
   !> ECP nuclei are those with `z < Z`; without `z`, or with a fill radius
   !> fraction <= 0, there are none
   !>
   !> @param[inout] self  Core fill
   !> @param[in]    mol   Molecular structure
   !> @param[in]    radii Cavity radius of each atom in Bohr
   !> @param[in]    param Isodensity parameters
   !> @param[in]    z     Nuclear charges of the host Hamiltonian (nat)
   subroutine iso_core_fill_update(self, mol, radii, param, z)
      class(moist_iso_core_fill_type), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Cavity radius of each atom in Bohr
      real(wp), intent(in) :: radii(:)
      !> Isodensity parameters
      type(moist_cavity_drop_lsf_isodensity_param_type), intent(in) :: param
      !> Nuclear charges of the host Hamiltonian
      real(wp), intent(in), optional :: z(:)

      !> Atom and fill indices, overlapping balls of one fill
      integer :: iat, ic, jc, nb
      !> Whether an atom is a filled ECP nucleus
      logical, allocatable :: ecp(:)

      self%n = 0
      self%overlap = 1
      self%rho_iso = param%rho_iso
      if (.not. param%core_fill .or. param%core_radius_fraction <= 0.0_wp) return
      if (.not. present(z)) return
      if (size(z) /= mol%nat .or. size(radii) < mol%nat) return

      allocate (ecp(mol%nat))
      do iat = 1, mol%nat
         ! A ball too small to square is no ball
         ecp(iat) = z(iat) < real(mol%num(mol%id(iat)), wp) .and. &
                    min(param%core_radius_fraction*radii(iat), max_radius)**2 > 0.0_wp
      end do
      if (.not. any(ecp)) return

      self%atom = pack([(iat, iat=1, mol%nat)], ecp)
      self%center = mol%xyz(:, self%atom)
      self%radius = min(param%core_radius_fraction*radii(self%atom), max_radius)
      self%radius2 = self%radius**2
      self%n = size(self%atom)
      self%amplitude = peak_ratio*param%rho_iso

      ! A point inside p balls makes each of them overlap the other p - 1
      do ic = 1, self%n
         nb = 1
         do jc = 1, self%n
            if (jc == ic) cycle
            if (norm2(self%center(:, ic) - self%center(:, jc)) < &
                sqrt(self%radius2(ic)) + sqrt(self%radius2(jc))) nb = nb + 1
         end do
         self%overlap = max(self%overlap, nb)
      end do
   end subroutine iso_core_fill_update

   !> Add the fill and its spatial derivatives at one point
   !>
   !> As a function of `u = |r - R|^2`, every derivative of `g` is `2^k d^k g/du^k`
   !> times products of `r - R` and Kronecker deltas; with `s = |r - R|` and
   !> `g' = dg/ds`, `2^k d^k g/du^k` is
   !>
   !>   k = 1: g'/s
   !>   k = 2: (g'' - g'/s)/s^2
   !>   k = 3: (g''' - 3 g''/s + 3 g'/s^2)/s^3
   !>   k = 4: (g'''' - 6 g'''/s + 15 g''/s^2 - 15 g'/s^3)/s^4
   !>
   !> whose terms share one sign at each order, so they never cancel
   !>
   !> @param[in]    self   Core fill
   !> @param[in]    point  Evaluation point in Bohr
   !> @param[in]    nderiv Highest derivative order to update (1 to 4)
   !> @param[inout] rho    Density
   !> @param[inout] drho   Density gradient
   !> @param[inout] d2rho  Density Hessian (nderiv >= 2)
   !> @param[inout] d3rho  Third density derivative (nderiv >= 3)
   !> @param[inout] d4rho  Fourth density derivative (nderiv >= 4)
   !> @param[out]   inside Filled atom whose ball holds the point (0: none)
   subroutine iso_core_fill_add(self, point, nderiv, rho, drho, d2rho, d3rho, d4rho, inside)
      class(moist_iso_core_fill_type), intent(in) :: self
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(ndim)
      !> Highest derivative order to update
      integer, intent(in) :: nderiv
      !> Density and its spatial derivatives
      real(wp), intent(inout) :: rho, drho(ndim)
      real(wp), intent(inout), optional :: d2rho(ndim, ndim), d3rho(ndim, ndim, ndim)
      real(wp), intent(inout), optional :: d4rho(ndim, ndim, ndim, ndim)
      !> Filled atom whose ball holds the point
      integer, intent(out) :: inside

      !> Offset from the nucleus, its square and length, 1/|r - R|, 1/X and 1 - |r - R|/X
      real(wp) :: d(ndim), u, s, si, q, w
      !> Radial derivatives `d^k g/ds^k` for k = 1..4
      real(wp) :: g1, g2, g3, g4
      !> `2^k d^k g/du^k` for k = 0..4
      real(wp) :: p0, p1, p2, p3, p4
      !> Kronecker deltas
      real(wp) :: delta(ndim, ndim)
      integer :: ic, i, j, k, l

      inside = 0
      if (self%n == 0) return
      delta = 0.0_wp
      do i = 1, ndim
         delta(i, i) = 1.0_wp
      end do

      do ic = 1, self%n
         d = point - self%center(:, ic)
         u = sum(d*d)
         ! Written so that a NaN point lies in no ball
         if (.not. (u < self%radius2(ic))) cycle
         if (inside == 0) inside = self%atom(ic)
         s = sqrt(u)
         q = 1.0_wp/self%radius(ic)
         w = 1.0_wp - q*s
         p0 = self%amplitude*w**5
         rho = rho + p0
         ! The tip of the cusp has no derivative; its neighbourhood would overflow
         if (s <= tip_radius*self%radius(ic)) cycle

         si = 1.0_wp/s
         g1 = -5.0_wp*self%amplitude*q*w**4
         g2 = 20.0_wp*self%amplitude*q**2*w**3
         p1 = g1*si
         p2 = (g2 - g1*si)*si**2

         drho = drho + p1*d
         if (nderiv >= 2 .and. present(d2rho)) then
            do j = 1, ndim
               do i = 1, ndim
                  d2rho(i, j) = d2rho(i, j) + p2*d(i)*d(j) + p1*delta(i, j)
               end do
            end do
         end if
         if ((nderiv >= 3 .and. present(d3rho)) .or. (nderiv >= 4 .and. present(d4rho))) then
            g3 = -60.0_wp*self%amplitude*q**3*w**2
            p3 = (g3 - 3.0_wp*g2*si + 3.0_wp*g1*si**2)*si**3
         end if
         if (nderiv >= 3 .and. present(d3rho)) then
            do k = 1, ndim
               do j = 1, ndim
                  do i = 1, ndim
                     d3rho(i, j, k) = d3rho(i, j, k) + p3*d(i)*d(j)*d(k) &
                        & + p2*(delta(i, j)*d(k) + delta(i, k)*d(j) + delta(j, k)*d(i))
                  end do
               end do
            end do
         end if
         if (nderiv >= 4 .and. present(d4rho)) then
            g4 = 120.0_wp*self%amplitude*q**4*w
            p4 = (g4 - 6.0_wp*g3*si + 15.0_wp*g2*si**2 - 15.0_wp*g1*si**3)*si**4
            do l = 1, ndim
               do k = 1, ndim
                  do j = 1, ndim
                     do i = 1, ndim
                        d4rho(i, j, k, l) = d4rho(i, j, k, l) + p4*d(i)*d(j)*d(k)*d(l) &
                           & + p3*(delta(i, j)*d(k)*d(l) + delta(i, k)*d(j)*d(l) &
                           & + delta(i, l)*d(j)*d(k) + delta(j, k)*d(i)*d(l) &
                           & + delta(j, l)*d(i)*d(k) + delta(k, l)*d(i)*d(j)) &
                           & + p2*(delta(i, j)*delta(k, l) + delta(i, k)*delta(j, l) &
                           & + delta(i, l)*delta(j, k))
                     end do
                  end do
               end do
            end do
         end if
      end do
   end subroutine iso_core_fill_add

   !> Reject a point inside a fill ball that is not inside the cavity
   !>
   !> The surface must never enter a ball, or the fill would move it and its
   !> nuclear motion would be missing from the level-set derivatives
   !>
   !> @param[in]  self   Core fill
   !> @param[in]  inside Filled atom whose ball holds the point (0: none)
   !> @param[in]  rho    Filled density at the point
   !> @param[out] error  Set when the point is not inside the cavity
   subroutine iso_core_fill_check(self, inside, rho, error)
      class(moist_iso_core_fill_type), intent(in) :: self
      !> Filled atom whose ball holds the point
      integer, intent(in) :: inside
      !> Filled density at the point
      real(wp), intent(in) :: rho
      !> Set when the point is not inside the cavity
      type(error_type), allocatable, intent(out) :: error

      if (inside == 0) return
      ! A NaN density is the solver's to reject, as without the fill
      if (.not. (rho <= self%rho_iso)) return
      call fill_error(self, findloc(self%atom, inside, 1), error)
   end subroutine iso_core_fill_check

   !> Fill whose ball holds a point, by the same test [[iso_core_fill_add]] uses
   !>
   !> @param[in] self  Core fill
   !> @param[in] point Point in Bohr
   !> @returns   ic    First fill whose ball holds the point (0: none)
   pure function iso_core_fill_holder(self, point) result(ic)
      class(moist_iso_core_fill_type), intent(in) :: self
      !> Point in Bohr
      real(wp), intent(in) :: point(ndim)
      !> First fill whose ball holds the point
      integer :: ic

      !> Offset from the nucleus
      real(wp) :: d(ndim)

      do ic = 1, self%n
         d = point - self%center(:, ic)
         if (sum(d*d) < self%radius2(ic)) return
      end do
      ic = 0
   end function iso_core_fill_holder

   !> Reject projected surface points inside a fill ball
   !>
   !> The fill is nonzero there, but none of its nuclear derivatives reach the
   !> level set; the per-point [[iso_core_fill_check]] cannot see such a point,
   !> whose filled density only matches `rho_iso` to the projection tolerance
   !>
   !> @param[in]  self      Core fill
   !> @param[in]  xyz       Cavity points in Bohr (3, n)
   !> @param[in]  projected Whether each point is a converged surface projection (n)
   !> @param[out] error     Set at the first projected point inside a ball
   subroutine iso_core_fill_check_surface(self, xyz, projected, error)
      class(moist_iso_core_fill_type), intent(in) :: self
      !> Cavity points in Bohr
      real(wp), intent(in) :: xyz(:, :)
      !> Whether each point is a converged surface projection
      logical, intent(in) :: projected(:)
      !> Set at the first projected point inside a ball
      type(error_type), allocatable, intent(out) :: error

      !> Surface point and the fill holding it
      integer :: i, ic

      if (self%n == 0) return
      do i = 1, size(xyz, 2)
         if (.not. projected(i)) cycle
         ic = self%holder(xyz(:, i))
         if (ic == 0) cycle
         call fill_error(self, ic, error, i, norm2(xyz(:, i) - self%center(:, ic)))
         return
      end do
   end subroutine iso_core_fill_check_surface

   !> Report a surface inside a fill ball, with the fill radius and how to resolve it
   !>
   !> Formats into fixed-length buffers: the per-point guard runs inside the
   !> threaded projection, where deferred-length function results are unsafe
   !>
   !> @param[in]  self   Core fill
   !> @param[in]  ic     Fill index
   !> @param[out] error  Error
   !> @param[in]  ipoint Index of the offending surface point (optional)
   !> @param[in]  dist   Its distance from the nucleus in Bohr (with `ipoint`)
   subroutine fill_error(self, ic, error, ipoint, dist)
      type(moist_iso_core_fill_type), intent(in) :: self
      !> Fill index
      integer, intent(in) :: ic
      !> Error
      type(error_type), allocatable, intent(out) :: error
      !> Index of the offending surface point
      integer, intent(in), optional :: ipoint
      !> Its distance from the nucleus in Bohr
      real(wp), intent(in), optional :: dist

      character(len=16) :: atom, radius, point, distance

      write (atom, "(i0)") self%atom(ic)
      write (radius, "(f16.3)") sqrt(self%radius2(ic))
      if (present(ipoint) .and. present(dist)) then
         write (point, "(i0)") ipoint
         write (distance, "(f16.3)") dist
         call fatal_error(error, "Isodensity surface point "//trim(point)//" lies "// &
                          trim(adjustl(distance))//" Bohr from the nucleus of atom "//trim(atom)// &
                          ", inside its ECP core fill (X = "//trim(adjustl(radius))//" Bohr); "// &
                          "lower rho_iso, lower core_radius_fraction, or disable core_fill")
      else
         call fatal_error(error, "Isodensity surface enters the ECP core fill of atom "//trim(atom)// &
                          ", where the filled density falls to rho_iso (X = "//trim(adjustl(radius))// &
                          " Bohr); lower rho_iso, lower core_radius_fraction, or disable core_fill")
      end if
   end subroutine fill_error

   !> Bound on `|grad g| / max(g, rho_iso)`, which adds to the interior bound on
   !> `|grad ln rho|` wherever `rho + g >= rho_iso`
   !>
   !> One ball gives `5 / ((1 - s) X)` where `g >= rho_iso` and `100 (1 - s)^4 5 / X`
   !> beyond, with `s = |r - R| / X`; both peak where `g = rho_iso`, at
   !> `(1 - s)^5 = 1/100`: `5 100^(1/5) / X`. With `g` summed over `p`
   !> overlapping balls each term is bounded alone, so the bound takes `p` times
   !> the steepest ball
   !>
   !> @param[in] self Core fill
   !> @returns   k    Bound in 1/Bohr (zero without fill)
   pure function iso_core_fill_log_grad_bound(self) result(k)
      class(moist_iso_core_fill_type), intent(in) :: self
      !> Bound in 1/Bohr
      real(wp) :: k

      k = 0.0_wp
      if (self%n == 0) return
      k = real(self%overlap, wp)*5.0_wp*peak_ratio**0.2_wp/minval(self%radius)
   end function iso_core_fill_log_grad_bound

end module moist_cavity_drop_lsf_isodensity_core
