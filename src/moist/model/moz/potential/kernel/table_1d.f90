!> Site-site potential tables of the 1D models
!>
!> - Ng-split arrays for listed (solute atom, solvent atom) pairs, from resolved sites
!> - Short range: u_sr = (u_pair + q_v phi_u) - ur_lr, shape (npts, npair)
!> - Real tail: ur_lr = q_u q_v erf(alpha r)/r, shape (npts, npair)
!> - Reciprocal tail: uk_lr = q_u q_v 4 pi exp(-k^2/(4 alpha^2))/k^2
!> - Reciprocal shape (nk, npair), zero at k = 0
!> - q_u, q_v: tail charges of resolved sites
!> - Full first-side radial field phi_u: supplied field (npts, nu), otherwise q_u/r
!> - u_pair: 12-6 potential with the interaction's mixing rule
!> - One-sided category or no categories: error; category absent on both sides: no contribution
!> - No Ng split or no first-side tail charges: zero tails, full Coulomb potential in u_sr
!> - No second-side field; radial grid excludes r = 0
!> - Caller validation: alpha positive and finite, valid mixing rule
!> - Reverse functional: L = sum w_sr u_sr + sum w_lr ur_lr + sum w_k uk_lr
!> - Field adjoint: w_phi(:, a) = sum q_v w_sr over pairs of atom a
!> - w_q: first-side tail-charge adjoint, fixed second-side sites
!> - No geometry in radial grids; no gradient
!> - Further adjoint propagation: `moz_potential_type%compute_adjoint`
module moist_model_moz_potential_kernel_table_1d
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_model_moz_potential_sites, only: potential_sites_type
   use moist_model_moz_potential_lj_base, only: lj_12_6_mix, lj_12_6_evaluate
   use moist_model_moz_potential_kernel_category, only: resolve_categories
   use moist_model_moz_potential_kernel_ng, only: ng_real_tail, ng_fourier_envelope
   implicit none(type, external)
   private

   public :: evaluate_table_1d, evaluate_table_1d_adjoint, require_radial_grids

contains

   !> Assemble the 1D Ng-split tables for the listed site pairs
   !>
   !> @param[in]  rgrid     real-space radial grid
   !> @param[in]  kgrid     reciprocal radial grid
   !> @param[in]  site_u    resolved sites of the first side
   !> @param[in]  site_v    resolved sites of the second side; the same sites for VV
   !> @param[in]  pairs     site pairs, (2, npair): atom of the first side, atom of the second
   !> @param[in]  mixing    lennard-jones mixing rule of the interaction
   !> @param[in]  ng_split  whether the Coulomb tail is split off
   !> @param[in]  alpha     ng split parameter, 1/bohr; validated by the caller
   !> @param[out] u_sr      short-range potential, Hartree (npts, npair)
   !> @param[out] ur_lr     long-range real-space potential, Hartree (npts, npair)
   !> @param[out] uk_lr     long-range reciprocal potential, Hartree bohr^3 (nk, npair)
   !> @param[out] error     category or input failure
   !> @param[in]  phi_u     optional first-side radial field, Hartree/e (npts, nu); required for field-supplying side
   subroutine evaluate_table_1d(rgrid, kgrid, site_u, site_v, pairs, mixing, ng_split, alpha, u_sr, ur_lr, uk_lr, &
         & error, phi_u)
      !> Real-space radial grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Resolved sites of the first side
      type(potential_sites_type), intent(in) :: site_u
      !> Resolved sites of the second side
      type(potential_sites_type), intent(in) :: site_v
      !> Site pairs, (2, npair)
      integer, intent(in) :: pairs(:, :)
      !> Lennard-Jones mixing rule of the interaction
      integer, intent(in) :: mixing
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Ng split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> Short-range potential (npts, npair)
      real(wp), allocatable, intent(out) :: u_sr(:, :)
      !> Long-range real-space potential (npts, npair)
      real(wp), allocatable, intent(out) :: ur_lr(:, :)
      !> Long-range reciprocal potential (nk, npair)
      real(wp), allocatable, intent(out) :: uk_lr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Radial field of the first side (npts, nu)
      real(wp), intent(in), optional :: phi_u(:, :)
      !> Pair potential, Coulomb part and the Ng real tail per node
      real(wp), allocatable :: u_pair(:), u_coul(:), f(:)
      !> Reciprocal envelope per node
      real(wp), allocatable :: env(:)
      logical :: do_pair, do_elec, tails
      real(wp) :: q_prod, sigma_ij, eps_ij
      integer :: ip, iu, iv, npts, nk, npair

      call check_inputs(rgrid, kgrid, site_u, site_v, pairs, do_pair, do_elec, error)
      if (allocated(error)) return
      if (do_elec .and. site_u%has_field) then
         if (.not. present(phi_u)) then
            call fatal_error(error, "1D tables need the radial field of a side with a field")
            return
         end if
         if (any(shape(phi_u) /= [rgrid%npts, site_u%natom])) then
            call fatal_error(error, "1D radial field must be (npts, natom) of the first side")
            return
         end if
      end if
      tails = do_elec .and. ng_split .and. allocated(site_u%q)

      npts = rgrid%npts
      nk = kgrid%npts
      npair = size(pairs, 2)
      allocate (u_sr(npts, npair), ur_lr(npts, npair), uk_lr(nk, npair), source=0.0_wp)
      allocate (u_pair(npts), u_coul(npts), f(npts), env(nk))
      if (tails) then
         call ng_real_tail(rgrid%r, alpha, f)
         env = ng_fourier_envelope(kgrid%r*kgrid%r, alpha)
      end if
      do ip = 1, npair
         iu = pairs(1, ip)
         iv = pairs(2, ip)
         if (do_pair) then
            call lj_12_6_mix(mixing, site_u%pair%sigma(iu), site_u%pair%epsilon(iu), &
               & site_v%pair%sigma(iv), site_v%pair%epsilon(iv), sigma_ij, eps_ij)
            call lj_12_6_evaluate(sigma_ij, eps_ij, rgrid%r, u_pair)
         else
            u_pair = 0.0_wp
         end if
         if (do_elec) then
            if (site_u%has_field) then
               u_coul = site_v%q(iv)*phi_u(:, iu)
            else
               u_coul = site_u%q(iu)*site_v%q(iv)/rgrid%r
            end if
            if (tails) then
               q_prod = site_u%q(iu)*site_v%q(iv)
               ur_lr(:, ip) = q_prod*f
               uk_lr(:, ip) = q_prod*env
            end if
            u_sr(:, ip) = (u_pair + u_coul) - ur_lr(:, ip)
         else
            u_sr(:, ip) = u_pair
         end if
      end do
   end subroutine evaluate_table_1d

   !> Form adjoints of `evaluate_table_1d` in the first-side field and tail charges
   !>
   !> - L = sum w_sr u_sr + sum w_lr ur_lr + sum w_k uk_lr
   !> - First-side atoms absent from pair list: zero adjoints
   !>
   !> @param[in]  rgrid     real-space radial grid
   !> @param[in]  kgrid     reciprocal radial grid
   !> @param[in]  site_u    resolved sites of the first side
   !> @param[in]  site_v    resolved sites of the second side
   !> @param[in]  pairs     site pairs, (2, npair): atom of the first side, atom of the second
   !> @param[in]  ng_split  whether the Coulomb tail is split off
   !> @param[in]  alpha     ng split parameter, 1/bohr; validated by the caller
   !> @param[in]  w_sr      adjoint of u_sr (npts, npair)
   !> @param[in]  w_lr      adjoint of ur_lr (npts, npair)
   !> @param[in]  w_k       adjoint of uk_lr (nk, npair)
   !> @param[out] w_phi     adjoint of the first-side field (npts, nu); unallocated without electrostatics
   !> @param[out] w_q       adjoint of the first-side tail charges (nu); zero-sized when the side has none
   !> @param[out] error     category or input failure
   subroutine evaluate_table_1d_adjoint(rgrid, kgrid, site_u, site_v, pairs, ng_split, alpha, w_sr, w_lr, w_k, &
         & w_phi, w_q, error)
      !> Real-space radial grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Resolved sites of the first side
      type(potential_sites_type), intent(in) :: site_u
      !> Resolved sites of the second side
      type(potential_sites_type), intent(in) :: site_v
      !> Site pairs, (2, npair)
      integer, intent(in) :: pairs(:, :)
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Ng split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> Adjoint of u_sr (npts, npair)
      real(wp), intent(in) :: w_sr(:, :)
      !> Adjoint of ur_lr (npts, npair)
      real(wp), intent(in) :: w_lr(:, :)
      !> Adjoint of uk_lr (nk, npair)
      real(wp), intent(in) :: w_k(:, :)
      !> Adjoint of the first-side field (npts, nu)
      real(wp), allocatable, intent(out) :: w_phi(:, :)
      !> Adjoint of the first-side tail charges (nu)
      real(wp), allocatable, intent(out) :: w_q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Ng real tail per node and reciprocal envelope per node
      real(wp), allocatable :: f(:), env(:)
      logical :: do_pair, do_elec, tails, point
      real(wp) :: q_v
      integer :: ip, iu, iv, npts, nk, npair

      call check_inputs(rgrid, kgrid, site_u, site_v, pairs, do_pair, do_elec, error)
      if (allocated(error)) return
      npts = rgrid%npts
      nk = kgrid%npts
      npair = size(pairs, 2)
      if (any(shape(w_sr) /= [npts, npair]) .or. any(shape(w_lr) /= [npts, npair])) then
         call fatal_error(error, "1D table adjoints w_sr and w_lr must be (npts, npair)")
         return
      end if
      if (any(shape(w_k) /= [nk, npair])) then
         call fatal_error(error, "1D table adjoint w_k must be (nk, npair)")
         return
      end if
      if (.not. do_elec) return

      tails = ng_split .and. allocated(site_u%q)
      point = allocated(site_u%q) .and. .not. site_u%has_field
      allocate (w_phi(npts, site_u%natom), source=0.0_wp)
      if (allocated(site_u%q)) then
         allocate (w_q(site_u%natom), source=0.0_wp)
      else
         allocate (w_q(0))
      end if
      if (tails) then
         allocate (f(npts))
         call ng_real_tail(rgrid%r, alpha, f)
         env = ng_fourier_envelope(kgrid%r*kgrid%r, alpha)
      end if
      do ip = 1, npair
         iu = pairs(1, ip)
         iv = pairs(2, ip)
         q_v = site_v%q(iv)
         w_phi(:, iu) = w_phi(:, iu) + q_v*w_sr(:, ip)
         if (point) w_q(iu) = w_q(iu) + q_v*sum(w_sr(:, ip)/rgrid%r)
         if (tails) w_q(iu) = w_q(iu) + q_v*(sum((w_lr(:, ip) - w_sr(:, ip))*f) + sum(w_k(:, ip)*env))
      end do
   end subroutine evaluate_table_1d_adjoint

   !> Resolve the categories and check the sites, grids and pair list
   !>
   !> @param[in]  rgrid    real-space radial grid
   !> @param[in]  kgrid    reciprocal radial grid
   !> @param[in]  site_u   resolved sites of the first side
   !> @param[in]  site_v   resolved sites of the second side
   !> @param[in]  pairs    site pairs, (2, npair)
   !> @param[out] do_pair  whether the pair category enters
   !> @param[out] do_elec  whether the electrostatic category enters
   !> @param[out] error    category or input failure
   subroutine check_inputs(rgrid, kgrid, site_u, site_v, pairs, do_pair, do_elec, error)
      !> Real-space radial grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Resolved sites of the first side
      type(potential_sites_type), intent(in) :: site_u
      !> Resolved sites of the second side
      type(potential_sites_type), intent(in) :: site_v
      !> Site pairs, (2, npair)
      integer, intent(in) :: pairs(:, :)
      !> Whether the pair category enters
      logical, intent(out) :: do_pair
      !> Whether the electrostatic category enters
      logical, intent(out) :: do_elec
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      do_pair = .false.
      do_elec = .false.
      if (site_u%natom < 1 .or. site_v%natom < 1) then
         call fatal_error(error, "1D MOZ tables need resolved sites on both sides")
         return
      end if
      if (site_v%has_field) then
         call fatal_error(error, "1D MOZ tables take no field on the second side")
         return
      end if
      call resolve_categories(allocated(site_u%pair), allocated(site_v%pair), &
         & site_u%has_field .or. allocated(site_u%q), allocated(site_v%q), do_pair, do_elec, error)
      if (allocated(error)) return
      call require_radial_grids(rgrid, kgrid, error)
      if (allocated(error)) return
      if (size(pairs, 1) /= 2) then
         call fatal_error(error, "1D pair list must have two rows, (2, npair)")
         return
      end if
      if (size(pairs, 2) < 1) then
         call fatal_error(error, "1D pair list is empty")
         return
      end if
      if (any(pairs(1, :) < 1) .or. any(pairs(1, :) > site_u%natom) .or. &
         & any(pairs(2, :) < 1) .or. any(pairs(2, :) > site_v%natom)) then
         call fatal_error(error, "1D pair list refers to atoms outside the two sides")
      end if
   end subroutine check_inputs

   !> Require built radial grids whose real-space nodes exclude r = 0
   !>
   !> @param[in]  rgrid  real-space radial grid
   !> @param[in]  kgrid  reciprocal radial grid
   !> @param[out] error  unbuilt grid or a node at r <= 0
   subroutine require_radial_grids(rgrid, kgrid, error)
      !> Real-space radial grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(rgrid%r) .or. rgrid%npts < 1) then
         call fatal_error(error, "1D MOZ needs a built real-space radial grid")
         return
      end if
      if (.not. allocated(kgrid%r) .or. kgrid%npts < 1) then
         call fatal_error(error, "1D MOZ needs a built reciprocal radial grid")
         return
      end if
      if (size(rgrid%r) /= rgrid%npts .or. size(kgrid%r) /= kgrid%npts) then
         call fatal_error(error, "1D MOZ radial grid node count does not match its nodes")
         return
      end if
      if (minval(rgrid%r) <= 0.0_wp) then
         call fatal_error(error, "1D MOZ real-space radial grid contains r <= 0, "// &
            & "where the pair and Coulomb potentials are singular")
      end if
   end subroutine require_radial_grids

end module moist_model_moz_potential_kernel_table_1d
