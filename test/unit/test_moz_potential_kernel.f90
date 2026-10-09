!> Tests of the MOZ potential tables: moz_potential_kernel
!>
!> - 1D and 3D Ng-split LJ and fixed-charge tables against dev-fork formulas
!> - Split identity, category rules and r = 0 guards
!> - 1D and 3D reverse modes against finite differences
module test_moz_potential_kernel
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_channels_coupling, only: coupling_type, coupling_view_type
   use moist_channels_response, only: response_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_model_moz_potential_term, only: potential_term_type, potential_name_len
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_potential_lj_base, only: lj_mixing_lorentz_berthelot, lj_mixing_geometric
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use moist_model_moz_potential_electrostatic_fixed, only: fixed_charges_type, new_fixed_charges
   use moist_model_moz_potential_kernel_ng, only: ng_real_tail, ng_fourier_envelope
   use moist_model_moz_potential_kernel_table_1d, only: evaluate_table_1d, evaluate_table_1d_adjoint
   use moist_model_moz_potential_sites, only: potential_sites_type
   use test_helpers, only: check_moist_error, fd4_scalar, fd4_offsets
   use test_moz_fixtures, only: ref_lj
   implicit none(type, external)
   private
   public :: collect_moz_potential_kernel
   public :: make_potential, make_solute, make_solvent, make_grid
   public :: check_close, solute_xyz, sigma_u, eps_u, q_u, sigma_v, eps_v, q_v, weights_3d, eval_functional
   public :: weights_1d, functional_1d, tables_1d, tables_3d, adjoint_3d
   public :: lb

   !> Pi
   real(wp), parameter :: pi = acos(-1.0_wp)
   !> Solute coordinates on integer offsets from their centroid, away from the half-integer grid nodes, bohr
   real(wp), parameter :: solute_xyz(3, 3) = reshape([-1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, &
      & 0.0_wp, -1.0_wp, 0.0_wp], [3, 3])
   !> Solute Lennard-Jones sigma, bohr
   real(wp), parameter :: sigma_u(3) = [0.9_wp, 0.8_wp, 0.6_wp]
   !> Solute Lennard-Jones epsilon, Hartree
   real(wp), parameter :: eps_u(3) = [0.004_wp, 0.006_wp, 0.002_wp]
   !> Solute charges, e
   real(wp), parameter :: q_u(3) = [0.35_wp, -0.55_wp, 0.12_wp]
   !> Solvent Lennard-Jones sigma, bohr
   real(wp), parameter :: sigma_v(2) = [1.0_wp, 0.4_wp]
   !> Solvent Lennard-Jones epsilon, Hartree
   real(wp), parameter :: eps_v(2) = [0.003_wp, 0.001_wp]
   !> Solvent charges, e
   real(wp), parameter :: q_v(2) = [-0.8_wp, 0.4_wp]
   !> Ng split parameter, 1/bohr
   real(wp), parameter :: alpha = 1.0_wp
   !> Mixing rule of the kernel tests
   integer, parameter :: lb = lj_mixing_lorentz_berthelot

   !> Field term from parameters alone
   !>
   !> - No host data
   !> - Volume field: phi(r) = x; radial field: phi(r, a) = a r
   type, extends(potential_term_type) :: linear_field_type
      !> Atoms of the latest build
      integer :: natom = 0
      !> Receives the radial field adjoint when associated (npts, natom)
      real(wp), pointer :: sink(:, :) => null()
   contains
      procedure :: name => linear_field_name
      procedure :: update => linear_field_update
      procedure :: has_field => linear_field_has_field
      procedure :: field_1d => linear_field_1d
      procedure :: field_3d => linear_field_3d
      procedure :: accumulate_adjoint_1d => linear_field_adjoint_1d
      procedure :: accumulate_adjoint_3d => linear_field_adjoint_3d
   end type linear_field_type

contains

   !> Register tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_potential_kernel(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("ng_tail_series", check_ng_tail_series), &
         & new_unittest("table_1d_reference", check_table_1d_reference), &
         & new_unittest("table_3d_reference", check_table_3d_reference), &
         & new_unittest("split_identity", check_split_identity), &
         & new_unittest("single_categories", check_single_categories), &
         & new_unittest("category_mismatch", check_category_mismatch), &
         & new_unittest("r0_guards", check_r0_guards), &
         & new_unittest("adjoint_3d_fd", check_adjoint_3d_fd), &
         & new_unittest("custom_field_adjoint_fd", check_custom_field_adjoint_fd), &
         & new_unittest("adjoint_1d_fd", check_adjoint_1d_fd), &
         & new_unittest("adjoint_1d_field_fd", check_adjoint_1d_field_fd), &
         & new_unittest("custom_radial_field_adjoint", check_custom_radial_field_adjoint)]
   end subroutine collect_moz_potential_kernel

   !* ================================================================================= *!
   !*                                   Shared helpers                                  *!
   !* ================================================================================= *!

   !> Update a potential from optional LJ parameters and optional fixed charges
   !>
   !> @param[out] pot    updated potential
   !> @param[in]  mol    structure of the side
   !> @param[out] err    build failure
   !> @param[in]  sigma  optional LJ sigma, bohr
   !> @param[in]  eps    optional LJ epsilon, Hartree
   !> @param[in]  q      optional fixed charges, e
   subroutine make_potential(pot, mol, err, sigma, eps, q)
      !> Updated potential
      type(moz_potential_type), intent(out) :: pot
      !> Structure of the side
      type(structure_type), intent(in) :: mol
      !> Build failure
      type(moist_error), allocatable, intent(out) :: err
      !> Optional LJ sigma, bohr
      real(wp), intent(in), optional :: sigma(:)
      !> Optional LJ epsilon, Hartree
      real(wp), intent(in), optional :: eps(:)
      !> Optional fixed charges, e
      real(wp), intent(in), optional :: q(:)
      type(custom_lj_type) :: lj
      type(fixed_charges_type) :: charges
      if (present(sigma)) then
         call new_custom_lj(lj, sigma, eps, err)
         if (allocated(err)) return
         call pot%add(lj, err)
         if (allocated(err)) return
      end if
      if (present(q)) then
         call new_fixed_charges(charges, q, err)
         if (allocated(err)) return
         call pot%add(charges, err)
         if (allocated(err)) return
      end if
      call pot%update(mol, err)
   end subroutine make_potential

   !> Construct a three-atom solute
   !>
   !> @param[out] mol  solute structure
   subroutine make_solute(mol)
      !> Solute structure
      type(structure_type), intent(out) :: mol
      call new(mol, [6, 8, 1], solute_xyz)
   end subroutine make_solute

   !> Construct a two-site solvent
   !>
   !> @param[out] mol  solvent structure
   subroutine make_solvent(mol)
      !> Solvent structure
      type(structure_type), intent(out) :: mol
      call new(mol, [8, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp, 0.0_wp, 0.0_wp], [3, 2]))
   end subroutine make_solvent

   !> Construct an 8^3 Cartesian point grid with unit spacing around the solute
   !>
   !> @param[out] grid  updated grid
   !> @param[in]  mol   solute structure
   !> @param[out] err   grid failure
   subroutine make_grid(grid, mol, err)
      !> Updated grid
      type(moist_math_grid_3d_cartesian_type), intent(out) :: grid
      !> Solute structure
      type(structure_type), intent(in) :: mol
      !> Grid failure
      type(moist_error), allocatable, intent(out) :: err
      call new_cartesian_point_grid(grid, err, nx=8, ny=8, nz=8, dr=1.0_wp, margin=0.0_wp)
      if (allocated(err)) return
      call grid%update(mol, err)
   end subroutine make_grid

   !> Form deterministic adjoint weights for the three 3D tables
   !>
   !> @param[in]  ngrid  number of grid points
   !> @param[in]  nk     number of reciprocal points
   !> @param[in]  nv     number of solvent atoms
   !> @param[out] w_sr   adjoint of u_sr (ngrid, nv)
   !> @param[out] w_lr   adjoint of ur_lr (ngrid, nv)
   !> @param[out] w_k    adjoint of uk_lr (nk, nv)
   subroutine weights_3d(ngrid, nk, nv, w_sr, w_lr, w_k)
      !> Number of grid points
      integer, intent(in) :: ngrid
      !> Number of reciprocal points
      integer, intent(in) :: nk
      !> Number of solvent atoms
      integer, intent(in) :: nv
      !> Adjoint of u_sr (ngrid, nv)
      real(wp), allocatable, intent(out) :: w_sr(:, :)
      !> Adjoint of ur_lr (ngrid, nv)
      real(wp), allocatable, intent(out) :: w_lr(:, :)
      !> Adjoint of uk_lr (nk, nv)
      complex(wp), allocatable, intent(out) :: w_k(:, :)
      integer :: g, v
      allocate (w_sr(ngrid, nv), w_lr(ngrid, nv), w_k(nk, nv))
      do v = 1, nv
         do g = 1, ngrid
            w_sr(g, v) = sin(1.37_wp*real(g, wp) + 0.71_wp*real(v, wp))
            w_lr(g, v) = cos(0.53_wp*real(g, wp) - 1.1_wp*real(v, wp))
         end do
         do g = 1, nk
            w_k(g, v) = cmplx(sin(0.91_wp*real(g, wp) + real(v, wp)), cos(1.21_wp*real(g, wp) - real(v, wp)), &
               & kind=wp)*1.0e-3_wp
         end do
      end do
   end subroutine weights_3d

   !> Form deterministic adjoint weights for the three 1D tables
   !>
   !> @param[in]  npts   real-space nodes
   !> @param[in]  nk     reciprocal nodes
   !> @param[in]  npair  site pairs
   !> @param[out] w_sr   adjoint of u_sr (npts, npair)
   !> @param[out] w_lr   adjoint of ur_lr (npts, npair)
   !> @param[out] w_k    adjoint of uk_lr (nk, npair)
   subroutine weights_1d(npts, nk, npair, w_sr, w_lr, w_k)
      !> Real-space nodes
      integer, intent(in) :: npts
      !> Reciprocal nodes
      integer, intent(in) :: nk
      !> Site pairs
      integer, intent(in) :: npair
      !> Adjoint of u_sr
      real(wp), allocatable, intent(out) :: w_sr(:, :)
      !> Adjoint of ur_lr
      real(wp), allocatable, intent(out) :: w_lr(:, :)
      !> Adjoint of uk_lr
      real(wp), allocatable, intent(out) :: w_k(:, :)
      integer :: i, p
      allocate (w_sr(npts, npair), w_lr(npts, npair), w_k(nk, npair))
      do p = 1, npair
         do i = 1, npts
            w_sr(i, p) = sin(1.37_wp*real(i, wp) + 0.71_wp*real(p, wp))
            w_lr(i, p) = cos(0.53_wp*real(i, wp) - 1.1_wp*real(p, wp))
         end do
         do i = 1, nk
            w_k(i, p) = sin(0.91_wp*real(i, wp) + real(p, wp))*1.0e-3_wp
         end do
      end do
   end subroutine weights_1d

   !> Evaluate the scalar functional of the three 1D tables of plain sites
   !>
   !> - L = sum w_sr u_sr + sum w_lr ur_lr + sum w_k uk_lr
   !>
   !> @param[in]  rgrid     real-space radial grid
   !> @param[in]  kgrid     reciprocal radial grid
   !> @param[in]  site_u    first-side sites
   !> @param[in]  site_v    second-side sites
   !> @param[in]  pairs     site pairs (2, npair)
   !> @param[in]  mixing    mixing rule
   !> @param[in]  ng_split  whether the Coulomb tail is split off
   !> @param[in]  w_sr      adjoint of u_sr
   !> @param[in]  w_lr      adjoint of ur_lr
   !> @param[in]  w_k       adjoint of uk_lr
   !> @param[out] val       functional value
   !> @param[out] err       assembly failure
   !> @param[in]  phi_u     optional first-side radial field (npts, nu)
   subroutine functional_1d(rgrid, kgrid, site_u, site_v, pairs, mixing, ng_split, w_sr, w_lr, w_k, val, err, phi_u)
      !> Real-space radial grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> First-side sites
      type(potential_sites_type), intent(in) :: site_u
      !> Second-side sites
      type(potential_sites_type), intent(in) :: site_v
      !> Site pairs
      integer, intent(in) :: pairs(:, :)
      !> Mixing rule
      integer, intent(in) :: mixing
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Adjoint of u_sr
      real(wp), intent(in) :: w_sr(:, :)
      !> Adjoint of ur_lr
      real(wp), intent(in) :: w_lr(:, :)
      !> Adjoint of uk_lr
      real(wp), intent(in) :: w_k(:, :)
      !> Functional value
      real(wp), intent(out) :: val
      !> Assembly failure
      type(moist_error), allocatable, intent(out) :: err
      !> First-side radial field
      real(wp), intent(in), optional :: phi_u(:, :)
      real(wp), allocatable :: sr(:, :), lr(:, :), uk(:, :)
      val = 0.0_wp
      call evaluate_table_1d(rgrid, kgrid, site_u, site_v, pairs, mixing, ng_split, alpha, sr, lr, uk, err, phi_u)
      if (allocated(err)) return
      val = sum(w_sr*sr) + sum(w_lr*lr) + sum(w_k*uk)
   end subroutine functional_1d

   !> Compare arrays with an absolute tolerance relaxed to relative above 1
   !>
   !> @param[out] error   test error
   !> @param[in]  actual  actual values
   !> @param[in]  ref     reference values
   !> @param[in]  tol     tolerance
   !> @param[in]  what    diagnostic label
   subroutine check_close(error, actual, ref, tol, what)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Actual values
      real(wp), intent(in) :: actual(:)
      !> Reference values
      real(wp), intent(in) :: ref(:)
      !> Tolerance
      real(wp), intent(in) :: tol
      !> Diagnostic label
      character(len=*), intent(in) :: what
      character(len=32) :: buffer
      real(wp) :: dev
      if (size(actual) /= size(ref)) then
         call test_failed(error, what//": size mismatch")
         return
      end if
      dev = maxval(abs(actual - ref)/max(1.0_wp, abs(ref)))
      if (dev > tol) then
         write (buffer, "(es12.4)") dev
         call test_failed(error, what//": deviation "//trim(buffer))
      end if
   end subroutine check_close

   !* ================================================================================= *!
   !*                                       Tests                                       *!
   !* ================================================================================= *!

   !> Check real Ng tail, radial factor and reciprocal envelope
   !>
   !> - Both sides of series switch and r = 0 limit
   !> - Dropped k = 0 mode
   subroutine check_ng_tail_series(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: radii(5) = [0.03_wp, 0.09_wp, 0.099_wp, 0.11_wp, 0.5_wp]
      real(wp), parameter :: alphas(2) = [1.0_wp, 0.7_wp]
      real(wp) :: f, g, x, ref_f, ref_g
      character(len=64) :: label
      integer :: i, j
      do j = 1, size(alphas)
         do i = 1, size(radii)
            call ng_real_tail(radii(i), alphas(j), f, g=g)
            x = alphas(j)*radii(i)
            ref_f = erf(x)/radii(i)
            ref_g = (2.0_wp/sqrt(pi)*x*exp(-x*x) - erf(x))/radii(i)**3
            write (label, "(a, f6.3, a, f4.2)") "Ng tail at r =", radii(i), ", alpha =", alphas(j)
            call check_close(error, [f], [ref_f], 1.0e-14_wp, trim(label)//" f")
            if (allocated(error)) return
            call check_close(error, [g], [ref_g], 1.0e-10_wp, trim(label)//" g")
            if (allocated(error)) return
         end do
         call ng_real_tail(0.0_wp, alphas(j), f, g)
         call check_close(error, [f, g], [2.0_wp*alphas(j)/sqrt(pi), -4.0_wp*alphas(j)**3/(3.0_wp*sqrt(pi))], &
            & 1.0e-15_wp, "Ng tail at r = 0")
         if (allocated(error)) return
      end do
      call check(error, ng_fourier_envelope(0.0_wp, 1.0_wp) == 0.0_wp, more="k = 0 mode must be dropped")
      if (allocated(error)) return
      call check_close(error, [ng_fourier_envelope(0.25_wp, 0.7_wp)], &
         & [4.0_wp*pi*exp(-0.25_wp/(4.0_wp*0.49_wp))/0.25_wp], 1.0e-14_wp, "Ng envelope")
   end subroutine check_ng_tail_series

   !> Compare 1D UV rectangle and VV triangle with dev formulas
   subroutine check_table_1d_reference(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :), uk_lr(:, :), ref_sr(:), ref_lr(:), ref_k(:)
      integer, allocatable :: pairs(:, :)
      real(wp) :: qq
      integer :: iu, iv, ip, side

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_potential(pot_u, mol_u, err, sigma_u, eps_u, q_u)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 64, 0.05_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      do side = 1, 2
         if (side == 1) then
            ! UV rectangle in the dev flat order ip = (iu - 1)*nv + iv
            call pot_u%site_pairs(pairs, err, solvent=pot_v)
            if (.not. allocated(err)) call check(error, all(pairs(:, 4) == [2, 2]), more="UV layout is not solute major")
            if (allocated(error)) return
            if (.not. allocated(err)) call tables_1d(rgrid, kgrid, pot_u, lb, .true., u_sr, ur_lr, uk_lr, err, solvent=pot_v)
         else
            ! VV upper triangle
            call pot_v%site_pairs(pairs, err)
            if (.not. allocated(err)) call check(error, all(pairs == reshape([1, 1, 1, 2, 2, 2], [2, 3])), &
               & more="VV layout is not the packed upper triangle")
            if (allocated(error)) return
            if (.not. allocated(err)) call tables_1d(rgrid, kgrid, pot_v, lb, .true., u_sr, ur_lr, uk_lr, err)
         end if
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, all(shape(u_sr) == [64, size(pairs, 2)]) .and. all(shape(uk_lr) == [64, size(pairs, 2)]))
         if (allocated(error)) return
         do ip = 1, size(pairs, 2)
            iu = pairs(1, ip)
            iv = pairs(2, ip)
            if (side == 1) then
               qq = q_u(iu)*q_v(iv)
               ref_sr = ref_lj(rgrid%r, eps_u(iu), sigma_u(iu), eps_v(iv), sigma_v(iv)) + qq/rgrid%r
            else
               qq = q_v(iu)*q_v(iv)
               ref_sr = ref_lj(rgrid%r, eps_v(iu), sigma_v(iu), eps_v(iv), sigma_v(iv)) + qq/rgrid%r
            end if
            ref_lr = qq*erf(alpha*rgrid%r)/rgrid%r
            ref_k = qq*4.0_wp*pi*exp(-0.25_wp*kgrid%r**2/alpha**2)/kgrid%r**2
            ref_sr = ref_sr - ref_lr
            call check_close(error, u_sr(:, ip), ref_sr, 1.0e-12_wp, "1D u_sr")
            if (allocated(error)) return
            call check_close(error, ur_lr(:, ip), ref_lr, 1.0e-12_wp, "1D ur_lr")
            if (allocated(error)) return
            call check_close(error, uk_lr(:, ip), ref_k, 1.0e-12_wp, "1D uk_lr")
            if (allocated(error)) return
         end do
      end do
   end subroutine check_table_1d_reference

   !> Compare 3D tables with dev formulas, including phased k-space sum
   subroutine check_table_3d_reference(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(coupling_type), target :: coupling
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(moist_math_grid_3d_cartesian_type) :: grid
      real(wp), allocatable :: u_sr(:, :), ur_lr(:, :), ref_sr(:, :), ref_lr(:, :)
      complex(wp), allocatable :: uk_lr(:, :), ref_k(:, :)
      real(wp) :: rv(3), kv(3), r0(3), r, lj, coul, lr, k2, env, re_sum, im_sum, kdotr
      integer :: g, j, a, v

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_potential(pot_u, mol_u, err, sigma_u, eps_u, q_u)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      if (.not. allocated(err)) call make_grid(grid, mol_u, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, grid%npts_k > 0, more="the grid has no reciprocal points")
      if (allocated(error)) return

      call tables_3d(grid, pot_u, pot_v, lb, .true., u_sr, ur_lr, uk_lr, err, coupling=coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return

      allocate (ref_sr(grid%ngrid, 2), ref_lr(grid%ngrid, 2), ref_k(grid%npts_k, 2))
      do g = 1, grid%ngrid
         rv = grid%point(g)
         do v = 1, 2
            lj = 0.0_wp
            coul = 0.0_wp
            lr = 0.0_wp
            do a = 1, 3
               r = norm2(rv - mol_u%xyz(:, a))
               lj = lj + ref_lj(r, eps_u(a), sigma_u(a), eps_v(v), sigma_v(v))
               coul = coul + q_u(a)*q_v(v)/r
               lr = lr + q_u(a)*q_v(v)*erf(alpha*r)/r
            end do
            ref_lr(g, v) = lr
            ref_sr(g, v) = (lj + coul) - lr
         end do
      end do
      r0 = grid%point(1)
      do j = 1, grid%npts_k
         kv = grid%kpoint(j)
         k2 = sum(kv**2)
         if (k2 < 1.0e-20_wp) then
            ref_k(j, :) = (0.0_wp, 0.0_wp)
            cycle
         end if
         env = 4.0_wp*pi*exp(-0.25_wp*k2/alpha**2)/k2
         re_sum = 0.0_wp
         im_sum = 0.0_wp
         do a = 1, 3
            kdotr = sum(kv*(mol_u%xyz(:, a) - r0))
            re_sum = re_sum + q_u(a)*cos(kdotr)
            im_sum = im_sum - q_u(a)*sin(kdotr)
         end do
         do v = 1, 2
            ref_k(j, v) = cmplx(q_v(v)*env*re_sum, q_v(v)*env*im_sum, kind=wp)
         end do
      end do
      call check_close(error, reshape(u_sr, [size(u_sr)]), reshape(ref_sr, [size(ref_sr)]), 1.0e-12_wp, "3D u_sr")
      if (allocated(error)) return
      call check_close(error, reshape(ur_lr, [size(ur_lr)]), reshape(ref_lr, [size(ref_lr)]), 1.0e-12_wp, "3D ur_lr")
      if (allocated(error)) return
      call check_close(error, reshape(real(uk_lr), [size(uk_lr)]), reshape(real(ref_k), [size(ref_k)]), &
         & 1.0e-12_wp, "3D Re uk_lr")
      if (allocated(error)) return
      call check_close(error, reshape(aimag(uk_lr), [size(uk_lr)]), reshape(aimag(ref_k), [size(ref_k)]), &
         & 1.0e-12_wp, "3D Im uk_lr")
   end subroutine check_table_3d_reference

   !> Check unsplit identity u_sr + ur_lr and zero tails without split
   subroutine check_split_identity(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(coupling_type), target :: coupling
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: sr(:, :), lr(:, :), sr0(:, :), lr0(:, :), k1(:, :), k0(:, :)
      complex(wp), allocatable :: uk(:, :), uk0(:, :)

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_potential(pot_u, mol_u, err, sigma_u, eps_u, q_u)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      if (.not. allocated(err)) call make_grid(grid, mol_u, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 32, 0.1_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      call tables_1d(rgrid, kgrid, pot_u, lb, .true., sr, lr, k1, err, solvent=pot_v)
      if (.not. allocated(err)) call tables_1d(rgrid, kgrid, pot_u, lb, .false., sr0, lr0, k0, err, solvent=pot_v)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_close(error, reshape(sr + lr, [size(sr)]), reshape(sr0, [size(sr0)]), 1.0e-12_wp, "1D split sum")
      if (allocated(error)) return
      call check(error, all(lr0 == 0.0_wp) .and. all(k0 == 0.0_wp), more="1D tails without the split")
      if (allocated(error)) return

      call tables_3d(grid, pot_u, pot_v, lb, .true., sr, lr, uk, err, coupling=coupling)
      if (.not. allocated(err)) call tables_3d(grid, pot_u, pot_v, lb, .false., sr0, lr0, uk0, err, coupling=coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_close(error, reshape(sr + lr, [size(sr)]), reshape(sr0, [size(sr0)]), 1.0e-12_wp, "3D split sum")
      if (allocated(error)) return
      call check(error, all(lr0 == 0.0_wp) .and. all(uk0 == (0.0_wp, 0.0_wp)), more="3D tails without the split")
      if (allocated(error)) return
   end subroutine check_split_identity

   !> Check pure pair tables for LJ-only potentials and Coulomb tables for charge-only potentials
   subroutine check_single_categories(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(coupling_type), target :: coupling
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(moist_math_grid_3d_cartesian_type) :: grid
      real(wp), allocatable :: sr(:, :), lr(:, :), ref(:, :)
      complex(wp), allocatable :: uk(:, :)
      real(wp) :: r
      integer :: g, v, a

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_grid(grid, mol_u, err)
      if (.not. allocated(err)) call make_potential(pot_u, mol_u, err, sigma=sigma_u, eps=eps_u)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma=sigma_v, eps=eps_v)
      if (.not. allocated(err)) call tables_3d(grid, pot_u, pot_v, lb, .true., sr, lr, uk, err, coupling=coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (ref(grid%ngrid, 2), source=0.0_wp)
      do g = 1, grid%ngrid
         do v = 1, 2
            do a = 1, 3
               r = norm2(grid%point(g) - mol_u%xyz(:, a))
               ref(g, v) = ref(g, v) + ref_lj(r, eps_u(a), sigma_u(a), eps_v(v), sigma_v(v))
            end do
         end do
      end do
      call check_close(error, reshape(sr, [size(sr)]), reshape(ref, [size(ref)]), 1.0e-12_wp, "LJ-only u_sr")
      if (allocated(error)) return
      call check(error, all(lr == 0.0_wp) .and. all(uk == (0.0_wp, 0.0_wp)), more="LJ-only tails")
      if (allocated(error)) return

      call make_potential(pot_u, mol_u, err, q=q_u)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, q=q_v)
      if (.not. allocated(err)) call tables_3d(grid, pot_u, pot_v, lb, .false., sr, lr, uk, err, coupling=coupling)
      call check_moist_error(error, err)
      if (allocated(error)) return
      ref = 0.0_wp
      do g = 1, grid%ngrid
         do v = 1, 2
            do a = 1, 3
               r = norm2(grid%point(g) - mol_u%xyz(:, a))
               ref(g, v) = ref(g, v) + q_u(a)*q_v(v)/r
            end do
         end do
      end do
      call check_close(error, reshape(sr, [size(sr)]), reshape(ref, [size(ref)]), 1.0e-12_wp, "charge-only u_sr")
   end subroutine check_single_categories

   !> Check one-sided category and empty-potential refusal in 1D and 3D
   subroutine check_category_mismatch(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(coupling_type), target :: coupling
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: sr(:, :), lr(:, :), k1(:, :)
      complex(wp), allocatable :: uk(:, :)
      !> Expected message fragment per case
      character(len=40), parameter :: expect(3) = [character(len=40) :: &
         & "Electrostatic potential on one side", "Pair potential on one side", "No potential terms"]
      type(potential_sites_type) :: site_u, site_v
      integer :: icase

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_grid(grid, mol_u, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 16, 0.1_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do icase = 1, 3
         select case (icase)
         case (1)
            call make_potential(pot_u, mol_u, err, sigma_u, eps_u, q_u)
            if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma=sigma_v, eps=eps_v)
         case (2)
            call make_potential(pot_u, mol_u, err, sigma=sigma_u, eps=eps_u)
            if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, q=q_v)
         case default
            ! A potential without terms is refused at update; the kernels refuse empty sites
            call make_potential(pot_u, mol_u, err)
            call check(error, allocated(err), more="a potential without terms was updated")
            if (allocated(error)) return
            call check(error, index(err%message, "has no terms") > 0, more=err%message)
            if (allocated(error)) return
            deallocate (err)
            site_u%natom = 3
            site_v%natom = 2
            call evaluate_table_1d(rgrid, kgrid, site_u, site_v, reshape([1, 1], [2, 1]), lb, .true., alpha, &
               & sr, lr, k1, err)
            call check(error, allocated(err), more="1D accepted: "//trim(expect(icase)))
            if (allocated(error)) return
            call check(error, index(err%message, trim(expect(icase))) > 0, more=err%message)
            return
         end select
         call check_moist_error(error, err)
         if (allocated(error)) return
         call tables_3d(grid, pot_u, pot_v, lb, .true., sr, lr, uk, err, coupling=coupling)
         call check(error, allocated(err), more="3D accepted: "//trim(expect(icase)))
         if (allocated(error)) return
         call check(error, index(err%message, trim(expect(icase))) > 0, more=err%message)
         if (allocated(error)) return
         call tables_1d(rgrid, kgrid, pot_u, lb, .true., sr, lr, k1, err, solvent=pot_v)
         call check(error, allocated(err), more="1D accepted: "//trim(expect(icase)))
         if (allocated(error)) return
         call check(error, index(err%message, trim(expect(icase))) > 0, more=err%message)
         if (allocated(error)) return
      end do
   end subroutine check_category_mismatch

   !> Check refusal of atom-coincident 3D grid points and r = 0 radial nodes
   subroutine check_r0_guards(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(coupling_type), target :: coupling
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: sr(:, :), lr(:, :), k1(:, :)
      type(structure_type) :: moved_mol
      complex(wp), allocatable :: uk(:, :)

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_grid(grid, mol_u, err)
      if (.not. allocated(err)) call make_potential(pot_u, mol_u, err, sigma_u, eps_u, q_u)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      call check_moist_error(error, err)
      if (allocated(error)) return
      moved_mol = mol_u
      moved_mol%xyz(:, 2) = grid%point(5)
      call tables_3d(grid, pot_u, pot_v, lb, .true., sr, lr, uk, err, mol=moved_mol)
      call check(error, allocated(err), more="a grid point on an atom was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "solute atom 2") > 0, more=err%message)
      if (allocated(error)) return

      rgrid%npts = 4
      rgrid%r = [0.0_wp, 0.1_wp, 0.2_wp, 0.3_wp]
      rgrid%w = [0.1_wp, 0.1_wp, 0.1_wp, 0.1_wp]
      kgrid = rgrid
      kgrid%r = kgrid%r + 0.05_wp
      call tables_1d(rgrid, kgrid, pot_u, lb, .true., sr, lr, k1, err, solvent=pot_v)
      call check(error, allocated(err), more="r = 0 on the radial grid was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "r <= 0") > 0, more=err%message)
   end subroutine check_r0_guards

   !> Evaluate the scalar functional of the three 3D tables
   !>
   !> - L = sum w_sr u_sr + sum w_lr ur_lr + Re sum conj(w_k) uk_lr
   !>
   !> @param[in]     grid      grid
   !> @param[in]     mol       solute structure at which the solute potential is updated
   !> @param[in]     pot_u     solute potential
   !> @param[in]     pot_v     solvent potential
   !> @param[in]     w_sr      adjoint of u_sr
   !> @param[in]     w_lr      adjoint of ur_lr
   !> @param[in]     w_k       adjoint of uk_lr
   !> @param[out]    val       functional value
   !> @param[out]    err       assembly failure
   !> @param[in]     mixing    optional mixing rule; Lorentz-Berthelot when absent
   !> @param[in]     ng_split  optional split flag; split when absent
   subroutine eval_functional(grid, mol, pot_u, pot_v, w_sr, w_lr, w_k, val, err, mixing, ng_split)
      !> Grid
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
      !> Solute structure at which the solute potential is updated
      type(structure_type), intent(in) :: mol
      !> Solute potential
      type(moz_potential_type), intent(in) :: pot_u
      !> Solvent potential
      type(moz_potential_type), intent(in) :: pot_v
      !> Adjoint of u_sr
      real(wp), intent(in) :: w_sr(:, :)
      !> Adjoint of ur_lr
      real(wp), intent(in) :: w_lr(:, :)
      !> Adjoint of uk_lr
      complex(wp), intent(in) :: w_k(:, :)
      !> Functional value
      real(wp), intent(out) :: val
      !> Assembly failure
      type(moist_error), allocatable, intent(out) :: err
      !> Mixing rule
      integer, intent(in), optional :: mixing
      !> Whether the Coulomb tail is split off
      logical, intent(in), optional :: ng_split
      type(coupling_type), target :: coupling
      real(wp), allocatable :: sr(:, :), lr(:, :)
      complex(wp), allocatable :: uk(:, :)
      integer :: mix
      logical :: split
      mix = lb
      if (present(mixing)) mix = mixing
      split = .true.
      if (present(ng_split)) split = ng_split
      val = 0.0_wp
      call tables_3d(grid, pot_u, pot_v, mix, split, sr, lr, uk, err, mol=mol, coupling=coupling)
      if (allocated(err)) return
      val = sum(w_sr*sr) + sum(w_lr*lr) + sum(real(conjg(w_k)*uk, wp))
   end subroutine eval_functional

   !> 1D tables of a potential at the given settings, computed on a copy
   !>
   !> @param[in]  rgrid     real-space radial grid
   !> @param[in]  kgrid     reciprocal radial grid
   !> @param[in]  pot_u     updated potential, first side
   !> @param[in]  mixing    mixing rule
   !> @param[in]  ng_split  whether the Coulomb tail is split off
   !> @param[out] u_sr      short-range table
   !> @param[out] ur_lr     long-range real-space table
   !> @param[out] uk_lr     long-range reciprocal table
   !> @param[out] err       setting or assembly failure
   !> @param[in]  solvent   optional solvent potential; UV when present
   subroutine tables_1d(rgrid, kgrid, pot_u, mixing, ng_split, u_sr, ur_lr, uk_lr, err, solvent)
      !> Real-space radial grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Reciprocal radial grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Updated potential, first side
      type(moz_potential_type), intent(in) :: pot_u
      !> Mixing rule
      integer, intent(in) :: mixing
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Short-range table
      real(wp), allocatable, intent(out) :: u_sr(:, :)
      !> Long-range real-space table
      real(wp), allocatable, intent(out) :: ur_lr(:, :)
      !> Long-range reciprocal table
      real(wp), allocatable, intent(out) :: uk_lr(:, :)
      !> Setting or assembly failure
      type(moist_error), allocatable, intent(out) :: err
      !> Optional solvent potential
      type(moz_potential_type), intent(in), optional :: solvent
      type(moz_potential_type) :: work
      work = pot_u
      call work%set_mixing(mixing, err)
      if (.not. allocated(err)) call work%set_ng_split(ng_split, alpha, err)
      if (.not. allocated(err)) call work%compute(rgrid, kgrid, u_sr, ur_lr, uk_lr, err, solvent=solvent)
   end subroutine tables_1d

   !> 3D UV tables of a solute potential at the given settings, computed on a copy
   !>
   !> @param[in]     grid      volume grid
   !> @param[in]     pot_u     updated solute potential
   !> @param[in]     pot_v     updated solvent potential
   !> @param[in]     mixing    mixing rule
   !> @param[in]     ng_split  whether the Coulomb tail is split off
   !> @param[out]    u_sr      short-range table
   !> @param[out]    ur_lr     long-range real-space table
   !> @param[out]    uk_lr     long-range reciprocal table
   !> @param[out]    err       setting, update or assembly failure
   !> @param[in]     mol       optional structure to update the copy at; the potential's own when absent
   !> @param[in,out] coupling  optional coupling for host-data terms
   subroutine tables_3d(grid, pot_u, pot_v, mixing, ng_split, u_sr, ur_lr, uk_lr, err, mol, coupling)
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Updated solute potential
      type(moz_potential_type), intent(in) :: pot_u
      !> Updated solvent potential
      type(moz_potential_type), intent(in) :: pot_v
      !> Mixing rule
      integer, intent(in) :: mixing
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Short-range table
      real(wp), allocatable, intent(out) :: u_sr(:, :)
      !> Long-range real-space table
      real(wp), allocatable, intent(out) :: ur_lr(:, :)
      !> Long-range reciprocal table
      complex(wp), allocatable, intent(out) :: uk_lr(:, :)
      !> Setting, update or assembly failure
      type(moist_error), allocatable, intent(out) :: err
      !> Optional structure to update the copy at
      type(structure_type), intent(in), optional :: mol
      !> Optional coupling for host-data terms
      type(coupling_type), intent(inout), target, optional :: coupling
      type(moz_potential_type) :: work
      work = pot_u
      if (present(mol)) call work%update(mol, err)
      if (.not. allocated(err)) call work%set_mixing(mixing, err)
      if (.not. allocated(err)) call work%set_ng_split(ng_split, alpha, err)
      if (.not. allocated(err)) call work%compute(grid, u_sr, ur_lr, uk_lr, err, solvent=pot_v, coupling=coupling)
   end subroutine tables_3d

   !> 3D reverse mode of a solute potential after setting its mixing rule and split
   !>
   !> @param[in]     grid          volume grid
   !> @param[in,out] pot_u         updated solute potential; settings changed in place
   !> @param[in]     pot_v         updated solvent potential
   !> @param[in]     mixing        mixing rule
   !> @param[in]     ng_split      whether the Coulomb tail is split off
   !> @param[in,out] coupling      coupling of the model
   !> @param[in]     w_sr          adjoint of u_sr
   !> @param[in]     w_lr          adjoint of ur_lr
   !> @param[in]     w_k           adjoint of uk_lr
   !> @param[in,out] response      host response items
   !> @param[out]    err           setting or adjoint failure
   !> @param[in,out] gradient      optional solute gradient
   !> @param[in,out] grid_adjoint  optional grid adjoint accumulator
   subroutine adjoint_3d(grid, pot_u, pot_v, mixing, ng_split, coupling, w_sr, w_lr, w_k, response, err, &
         & gradient, grid_adjoint)
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Updated solute potential
      type(moz_potential_type), intent(inout) :: pot_u
      !> Updated solvent potential
      type(moz_potential_type), intent(in) :: pot_v
      !> Mixing rule
      integer, intent(in) :: mixing
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Coupling of the model
      type(coupling_type), intent(inout), target :: coupling
      !> Adjoint of u_sr
      real(wp), intent(in) :: w_sr(:, :)
      !> Adjoint of ur_lr
      real(wp), intent(in) :: w_lr(:, :)
      !> Adjoint of uk_lr
      complex(wp), intent(in) :: w_k(:, :)
      !> Host response items
      type(response_type), intent(inout) :: response
      !> Setting or adjoint failure
      type(moist_error), allocatable, intent(out) :: err
      !> Optional solute gradient
      real(wp), intent(inout), optional :: gradient(:, :)
      !> Optional grid adjoint accumulator
      type(volume_adjoint_type), intent(inout), optional :: grid_adjoint
      call pot_u%set_mixing(mixing, err)
      if (.not. allocated(err)) call pot_u%set_ng_split(ng_split, alpha, err)
      if (.not. allocated(err)) call pot_u%compute_adjoint(grid, pot_v, coupling, w_sr, w_lr, w_k, response, err, &
         & gradient, grid_adjoint)
   end subroutine adjoint_3d

   !> Compare 3D reverse mode with four-point central differences
   !>
   !> - Solute and grid coordinates; total translation invariance
   !> - Both mixing rules, with and without Ng split
   !> - No geometry derivatives: no adjoint dispatch for parameter-only potential
   subroutine check_adjoint_3d_fd(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(coupling_type), target :: coupling
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(response_type) :: response
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), gradient(:, :)
      complex(wp), allocatable :: w_k(:, :)
      !> Mixing rules probed
      integer, parameter :: mixings(2) = [lj_mixing_lorentz_berthelot, lj_mixing_geometric]
      !> Split flags probed
      logical, parameter :: splits(2) = [.true., .false.]
      integer :: im, is

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_grid(grid, mol_u, err)
      if (.not. allocated(err)) call make_potential(pot_u, mol_u, err, sigma_u, eps_u, q_u)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call weights_3d(grid%ngrid, grid%npts_k, 2, w_sr, w_lr, w_k)
      do im = 1, size(mixings)
         do is = 1, size(splits)
            call check_adjoint_3d_case(error, grid, mol_u, pot_u, pot_v, mixings(im), splits(is), w_sr, w_lr, w_k)
            if (allocated(error)) return
         end do
      end do

      ! Response phase: no geometry, no parameter-only adjoint dispatch
      call adjoint_3d(grid, pot_u, pot_v, lb, .true., coupling, w_sr, w_lr, w_k, response, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. response%next(), more="a parameter-only potential pushed response items")
      if (allocated(error)) return
      allocate (gradient(3, 3), source=0.0_wp)
      call adjoint_3d(grid, pot_u, pot_v, lb, .true., coupling, w_sr, w_lr, w_k, response, err, gradient=gradient)
      call check(error, allocated(err), more="a gradient without a grid adjoint was accepted")
   end subroutine check_adjoint_3d_fd

   !> Check one mixing rule and split flag of `check_adjoint_3d_fd`
   !>
   !> @param[out]    error     test error
   !> @param[in]     grid      grid
   !> @param[in]     mol_u     solute
   !> @param[in,out] pot_u     solute potential
   !> @param[in]     pot_v     solvent potential
   !> @param[in]     mixing    mixing rule
   !> @param[in]     ng_split  whether the Coulomb tail is split off
   !> @param[in]     w_sr      adjoint of u_sr
   !> @param[in]     w_lr      adjoint of ur_lr
   !> @param[in]     w_k       adjoint of uk_lr
   subroutine check_adjoint_3d_case(error, grid, mol_u, pot_u, pot_v, mixing, ng_split, w_sr, w_lr, w_k)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Grid
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
      !> Solute
      type(structure_type), intent(in) :: mol_u
      !> Solute potential
      type(moz_potential_type), intent(inout) :: pot_u
      !> Solvent potential
      type(moz_potential_type), intent(in) :: pot_v
      !> Mixing rule
      integer, intent(in) :: mixing
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Adjoint of u_sr
      real(wp), intent(in) :: w_sr(:, :)
      !> Adjoint of ur_lr
      real(wp), intent(in) :: w_lr(:, :)
      !> Adjoint of uk_lr
      complex(wp), intent(in) :: w_k(:, :)
      type(moist_error), allocatable :: err
      type(coupling_type), target :: coupling
      type(moist_math_grid_3d_cartesian_type) :: moved
      type(volume_adjoint_type) :: acc
      type(response_type) :: response
      real(wp), allocatable :: gradient(:, :)
      type(structure_type) :: moved_mol
      real(wp), parameter :: h = 1.0e-4_wp
      !> Grid points probed: the phase reference, points near atoms, the last point
      integer :: probe(4)
      real(wp) :: val(4), fd, total(3)
      character(len=32) :: case_label
      character(len=128) :: label
      integer :: a, c, s, ip, g

      write (case_label, "(a, i0, a, l1, a)") "mixing ", mixing, ", split ", ng_split, ": "
      allocate (gradient(3, 3), source=0.0_wp)
      call acc%init(grid%ngrid)
      call adjoint_3d(grid, pot_u, pot_v, mixing, ng_split, coupling, w_sr, w_lr, w_k, response, err, &
         & gradient=gradient, grid_adjoint=acc)
      call check_moist_error(error, err)
      if (allocated(error)) return

      do a = 1, 3
         do c = 1, 3
            do s = 1, 4
               moved_mol = mol_u
               moved_mol%xyz(c, a) = moved_mol%xyz(c, a) + fd4_offsets(s)*h
               call eval_functional(grid, moved_mol, pot_u, pot_v, w_sr, w_lr, w_k, val(s), err, mixing, ng_split)
               call check_moist_error(error, err)
               if (allocated(error)) return
            end do
            call fd4_scalar(val(1), val(2), val(3), val(4), h, fd, error)
            if (allocated(error)) return
            write (label, "(a, a, i0, a, i0)") trim(case_label), " solute gradient atom ", a, " component ", c
            call check_close(error, [gradient(c, a)], [fd], 1.0e-7_wp, trim(label))
            if (allocated(error)) return
         end do
      end do

      probe = [1, 2, grid%ngrid/2 + 4, grid%ngrid]
      do ip = 1, size(probe)
         g = probe(ip)
         do c = 1, 3
            do s = 1, 4
               moved = grid
               moved%xyz(c, g) = moved%xyz(c, g) + fd4_offsets(s)*h
               call eval_functional(moved, mol_u, pot_u, pot_v, w_sr, w_lr, w_k, val(s), err, mixing, ng_split)
               call check_moist_error(error, err)
               if (allocated(error)) return
            end do
            call fd4_scalar(val(1), val(2), val(3), val(4), h, fd, error)
            if (allocated(error)) return
            write (label, "(a, a, i0, a, i0)") trim(case_label), " grid adjoint point ", g, " component ", c
            call check_close(error, [acc%w_xyz(c, g)], [fd], 1.0e-7_wp, trim(label))
            if (allocated(error)) return
         end do
      end do

      ! Table invariance under rigid translation of grid and solute
      total = sum(acc%w_xyz, dim=2) + sum(gradient, dim=2)
      call check(error, maxval(abs(total)) < 1.0e-10_wp, more=trim(case_label)//" translation invariance")
      if (allocated(error)) return
      call check(error, all(acc%w_w == 0.0_wp) .and. all(acc%w_xi == 0.0_wp), &
         & more=trim(case_label)//" parameter-built tables touched the weight or width adjoints")
   end subroutine check_adjoint_3d_case

   !> Compare parameter-field adjoint with four-point central differences in grid coordinates
   !>
   !> - No host data; volume field phi(r) = x
   !> - LJ on both sides, charges on solvent only
   subroutine check_custom_field_adjoint_fd(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(coupling_type), target :: coupling
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(custom_lj_type) :: lj
      type(linear_field_type) :: field
      type(moist_math_grid_3d_cartesian_type) :: grid, moved
      type(volume_adjoint_type) :: acc
      type(response_type) :: response
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), gradient(:, :)
      complex(wp), allocatable :: w_k(:, :)
      real(wp), parameter :: h = 1.0e-4_wp
      !> Grid points probed: the phase reference, points near atoms, the last point
      integer :: probe(4)
      real(wp) :: val(4), fd
      character(len=96) :: label
      integer :: c, s, ip, g

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_grid(grid, mol_u, err)
      if (.not. allocated(err)) call new_custom_lj(lj, sigma_u, eps_u, err)
      if (.not. allocated(err)) call pot_u%add(lj, err)
      if (.not. allocated(err)) call pot_u%add(field, err)
      if (.not. allocated(err)) call pot_u%update(mol_u, err)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. field%is_host_fed() .and. field%has_field(), &
         & more="the test term must supply a field without host data")
      if (allocated(error)) return
      call weights_3d(grid%ngrid, grid%npts_k, 2, w_sr, w_lr, w_k)
      allocate (gradient(3, 3), source=0.0_wp)
      call acc%init(grid%ngrid)
      call adjoint_3d(grid, pot_u, pot_v, lb, .true., coupling, w_sr, w_lr, w_k, response, err, gradient=gradient, grid_adjoint=acc)
      call check_moist_error(error, err)
      if (allocated(error)) return

      probe = [1, 2, grid%ngrid/2 + 4, grid%ngrid]
      do ip = 1, size(probe)
         g = probe(ip)
         do c = 1, 3
            do s = 1, 4
               moved = grid
               moved%xyz(c, g) = moved%xyz(c, g) + fd4_offsets(s)*h
               call eval_functional(moved, mol_u, pot_u, pot_v, w_sr, w_lr, w_k, val(s), err)
               call check_moist_error(error, err)
               if (allocated(error)) return
            end do
            call fd4_scalar(val(1), val(2), val(3), val(4), h, fd, error)
            if (allocated(error)) return
            write (label, "(a, i0, a, i0)") "grid adjoint point ", g, " component ", c
            call check_close(error, [acc%w_xyz(c, g)], [fd], 1.0e-7_wp, trim(label))
            if (allocated(error)) return
         end do
      end do
   end subroutine check_custom_field_adjoint_fd

   !> Compare 1D charge reverse mode with four-point central differences
   !>
   !> - First-side charges; both mixing rules, with and without Ng split
   !> - Repeated pair; solute atom 3 absent, zero adjoints
   !> - w_phi(:, a) = sum q_v w_sr over pairs of atom a
   subroutine check_adjoint_1d_fd(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(potential_sites_type) :: site_u, site_v, moved
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), w_k(:, :), w_phi(:, :), w_q(:), ref_phi(:, :)
      !> Site pairs: (1, 1) twice, solute atom 3 absent
      integer, parameter :: pairs(2, 5) = reshape([1, 1, 1, 2, 2, 1, 2, 2, 1, 1], [2, 5])
      !> Mixing rules probed
      integer, parameter :: mixings(2) = [lj_mixing_lorentz_berthelot, lj_mixing_geometric]
      !> Split flags probed
      logical, parameter :: splits(2) = [.true., .false.]
      !> - L linear in charges and field; exact central differences for any step
      !> - Step large enough to exceed LJ-wall roundoff at first node
      real(wp), parameter :: h = 0.1_wp
      real(wp) :: val(4), fd
      character(len=32) :: case_label
      character(len=96) :: label
      integer :: im, is, a, s, ip

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_potential(pot_u, mol_u, err, sigma_u, eps_u, q_u)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      if (.not. allocated(err)) call pot_u%sites(site_u, err)
      if (.not. allocated(err)) call pot_v%sites(site_v, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 32, 0.6_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call weights_1d(rgrid%npts, kgrid%npts, size(pairs, 2), w_sr, w_lr, w_k)
      allocate (ref_phi(rgrid%npts, 3), source=0.0_wp)
      do ip = 1, size(pairs, 2)
         ref_phi(:, pairs(1, ip)) = ref_phi(:, pairs(1, ip)) + q_v(pairs(2, ip))*w_sr(:, ip)
      end do

      do im = 1, size(mixings)
         do is = 1, size(splits)
            write (case_label, "(a, i0, a, l1, a)") "mixing ", mixings(im), ", split ", splits(is), ":"
            call evaluate_table_1d_adjoint(rgrid, kgrid, site_u, site_v, pairs, splits(is), alpha, w_sr, w_lr, w_k, &
               & w_phi, w_q, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, all(shape(w_phi) == [rgrid%npts, 3]) .and. size(w_q) == 3, &
               & more=trim(case_label)//" adjoint shapes")
            if (allocated(error)) return
            call check_close(error, reshape(w_phi, [size(w_phi)]), reshape(ref_phi, [size(ref_phi)]), 1.0e-14_wp, &
               & trim(case_label)//" w_phi")
            if (allocated(error)) return
            call check(error, w_q(3) == 0.0_wp, more=trim(case_label)//" atom outside the pair list has a charge adjoint")
            if (allocated(error)) return
            do a = 1, 3
               do s = 1, 4
                  moved = site_u
                  moved%q(a) = moved%q(a) + fd4_offsets(s)*h
                  call functional_1d(rgrid, kgrid, moved, site_v, pairs, mixings(im), splits(is), w_sr, w_lr, w_k, &
                     & val(s), err)
                  call check_moist_error(error, err)
                  if (allocated(error)) return
               end do
               call fd4_scalar(val(1), val(2), val(3), val(4), h, fd, error)
               if (allocated(error)) return
               write (label, "(a, a, i0)") trim(case_label), " charge adjoint atom ", a
               call check_close(error, [w_q(a)], [fd], 1.0e-9_wp, trim(label))
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine check_adjoint_1d_fd

   !> Compare 1D radial-field reverse mode with four-point central differences
   !>
   !> - With and without Ng split
   !> - w_phi against differences in single field entries
   !> - w_q (tail only) against differences in charges
   !> - Field without charges: no tails, empty w_q
   subroutine check_adjoint_1d_field_fd(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(potential_sites_type) :: site_u, site_v, moved
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), w_k(:, :), w_phi(:, :), w_q(:), phi(:, :), phi_moved(:, :)
      real(wp), allocatable :: sr(:, :), lr(:, :), uk(:, :)
      !> Site pairs, solute major
      integer, parameter :: pairs(2, 6) = reshape([1, 1, 1, 2, 2, 1, 2, 2, 3, 1, 3, 2], [2, 6])
      !> Split flags probed
      logical, parameter :: splits(2) = [.true., .false.]
      !> - L linear in charges and field; exact central differences for any step
      !> - Step large enough to exceed LJ-wall roundoff at first node
      real(wp), parameter :: h = 0.1_wp
      integer :: probe(3)
      real(wp) :: val(4), fd
      character(len=32) :: case_label
      character(len=96) :: label
      integer :: is, a, s, ip, i

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call make_potential(pot_u, mol_u, err, sigma_u, eps_u, q_u)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      if (.not. allocated(err)) call pot_u%sites(site_u, err)
      if (.not. allocated(err)) call pot_v%sites(site_v, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, 32, 0.6_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      site_u%has_field = .true.
      allocate (phi(rgrid%npts, 3))
      do a = 1, 3
         phi(:, a) = q_u(a)*exp(-0.3_wp*rgrid%r)/rgrid%r + 0.01_wp*real(a, wp)
      end do
      call weights_1d(rgrid%npts, kgrid%npts, size(pairs, 2), w_sr, w_lr, w_k)
      probe = [1, rgrid%npts/2, rgrid%npts]

      call evaluate_table_1d(rgrid, kgrid, site_u, site_v, pairs, lb, .true., alpha, sr, lr, uk, err)
      call check(error, allocated(err), more="a side with a field was assembled without its field")
      if (allocated(error)) return
      deallocate (err)

      do is = 1, size(splits)
         write (case_label, "(a, l1, a)") "field, split ", splits(is), ":"
         call evaluate_table_1d_adjoint(rgrid, kgrid, site_u, site_v, pairs, splits(is), alpha, w_sr, w_lr, w_k, &
            & w_phi, w_q, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         do a = 1, 3
            do ip = 1, size(probe)
               i = probe(ip)
               do s = 1, 4
                  phi_moved = phi
                  phi_moved(i, a) = phi_moved(i, a) + fd4_offsets(s)*h
                  call functional_1d(rgrid, kgrid, site_u, site_v, pairs, lb, splits(is), w_sr, w_lr, w_k, &
                     & val(s), err, phi_moved)
                  call check_moist_error(error, err)
                  if (allocated(error)) return
               end do
               call fd4_scalar(val(1), val(2), val(3), val(4), h, fd, error)
               if (allocated(error)) return
               write (label, "(a, a, i0, a, i0)") trim(case_label), " field adjoint atom ", a, " node ", i
               call check_close(error, [w_phi(i, a)], [fd], 1.0e-9_wp, trim(label))
               if (allocated(error)) return
            end do
            do s = 1, 4
               moved = site_u
               moved%q(a) = moved%q(a) + fd4_offsets(s)*h
               call functional_1d(rgrid, kgrid, moved, site_v, pairs, lb, splits(is), w_sr, w_lr, w_k, val(s), err, phi)
               call check_moist_error(error, err)
               if (allocated(error)) return
            end do
            call fd4_scalar(val(1), val(2), val(3), val(4), h, fd, error)
            if (allocated(error)) return
            write (label, "(a, a, i0)") trim(case_label), " tail charge adjoint atom ", a
            call check_close(error, [w_q(a)], [fd], 1.0e-9_wp, trim(label))
            if (allocated(error)) return
         end do
         if (.not. splits(is)) then
            call check(error, all(w_q == 0.0_wp), more="without the split a field side has tail charge adjoints")
            if (allocated(error)) return
         end if
      end do

      deallocate (site_u%q)
      call evaluate_table_1d_adjoint(rgrid, kgrid, site_u, site_v, pairs, .true., alpha, w_sr, w_lr, w_k, &
         & w_phi, w_q, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, allocated(w_phi) .and. size(w_q) == 0, more="a field without charges has charge adjoints")
   end subroutine check_adjoint_1d_field_fd

   !> Check radial parameter-field assembly and adjoint dispatch
   !>
   !> - No host data; field sum matches kernel input
   !> - Kernel field adjoint returned to term; empty response
   subroutine check_custom_radial_field_adjoint(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(coupling_type), target :: coupling
      type(structure_type) :: mol_u, mol_v
      type(moz_potential_type) :: pot_u, pot_v
      type(potential_sites_type) :: site_u, site_v
      type(custom_lj_type) :: lj
      type(linear_field_type) :: field
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(response_type) :: response
      real(wp), allocatable :: w_sr(:, :), w_lr(:, :), w_k(:, :), w_phi(:, :), w_q(:), phi(:, :)
      real(wp), allocatable :: sr(:, :), lr(:, :), uk(:, :), sr_ref(:, :), lr_ref(:, :), uk_ref(:, :)
      real(wp), allocatable, target :: sink(:, :)
      !> Site pairs, solute major
      integer, parameter :: pairs(2, 6) = reshape([1, 1, 1, 2, 2, 1, 2, 2, 3, 1, 3, 2], [2, 6])
      integer :: a

      call make_solute(mol_u)
      call make_solvent(mol_v)
      call new_uniform_radial_pair(rgrid, kgrid, 32, 0.6_wp, err)
      if (.not. allocated(err)) call new_custom_lj(lj, sigma_u, eps_u, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (sink(rgrid%npts, 3), source=0.0_wp)
      field%sink => sink
      call pot_u%add(lj, err)
      if (.not. allocated(err)) call pot_u%add(field, err)
      if (.not. allocated(err)) call pot_u%update(mol_u, err)
      if (.not. allocated(err)) call make_potential(pot_v, mol_v, err, sigma_v, eps_v, q_v)
      if (.not. allocated(err)) call pot_u%sites(site_u, err)
      if (.not. allocated(err)) call pot_v%sites(site_v, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      allocate (phi(rgrid%npts, 3))
      do a = 1, 3
         phi(:, a) = real(a, wp)*rgrid%r
      end do
      call evaluate_table_1d(rgrid, kgrid, site_u, site_v, pairs, lb, .true., alpha, sr_ref, lr_ref, uk_ref, err, phi)
      if (.not. allocated(err)) call pot_u%compute(rgrid, kgrid, sr, lr, uk, err, solvent=pot_v)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(sr == sr_ref) .and. all(lr == lr_ref) .and. all(uk == uk_ref), &
         & more="the potential's radial field does not reach the tables")
      if (allocated(error)) return

      call weights_1d(rgrid%npts, kgrid%npts, size(pairs, 2), w_sr, w_lr, w_k)
      call evaluate_table_1d_adjoint(rgrid, kgrid, site_u, site_v, pairs, .true., alpha, w_sr, w_lr, w_k, &
         & w_phi, w_q, err)
      if (.not. allocated(err)) call pot_u%compute_adjoint(rgrid, kgrid, pot_v, coupling, w_sr, w_lr, w_k, response, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(sink == w_phi), more="the field term did not receive the kernel's field adjoint")
      if (allocated(error)) return
      call check(error, .not. response%next(), more="a parameter-built field pushed response items")
   end subroutine check_custom_radial_field_adjoint

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function linear_field_name(self) result(name)
      !> Term
      class(linear_field_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "linear_field"
   end function linear_field_name

   !> Record the number of atoms
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         structure of this side
   !> @param[out]    error       never set
   !> @param[in]     solvent_id  unused
   subroutine linear_field_update(self, mol, error, solvent_id)
      !> Term
      class(linear_field_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(moist_error), allocatable, intent(out) :: error
      !> Unused
      integer, intent(in), optional :: solvent_id
      self%natom = mol%nat
   end subroutine linear_field_update

   !> Grid-field presence
   !>
   !> @param[in] self  term
   pure function linear_field_has_field(self) result(has)
      !> Term
      class(linear_field_type), intent(in) :: self
      !> Presence
      logical :: has
      has = .true.
   end function linear_field_has_field

   !> Form phi(r, a) = a r on the radial grid
   !>
   !> @param[in]  self     built term
   !> @param[in]  grid     radial grid
   !> @param[in]  view     unbound view; unused
   !> @param[out] phi      field per atom, Hartree/e (npts, natom)
   !> @param[out] error    never set
   subroutine linear_field_1d(self, grid, view, phi, error)
      !> Built term
      class(linear_field_type), intent(in) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Unbound view
      type(coupling_view_type), intent(in) :: view
      !> Field per atom (npts, natom)
      real(wp), allocatable, intent(out) :: phi(:, :)
      !> Error handling
      type(moist_error), allocatable, intent(out) :: error
      integer :: a
      allocate (phi(grid%npts, self%natom))
      do a = 1, self%natom
         phi(:, a) = real(a, wp)*grid%r
      end do
   end subroutine linear_field_1d

   !> Add the radial field adjoint to the sink, when one is attached
   !>
   !> @param[in,out] self      term
   !> @param[in]     grid      radial grid
   !> @param[in]     view      unbound view; unused
   !> @param[in]     w_phi     adjoint of the field (npts, natom)
   !> @param[in]     w_qinf    adjoint of the tail charges; zero-sized, the term has none
   !> @param[in,out] response  host response items; untouched
   !> @param[out]    error     never set
   subroutine linear_field_adjoint_1d(self, grid, view, w_phi, w_qinf, response, error)
      !> Term
      class(linear_field_type), intent(inout) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Unbound view
      type(coupling_view_type), intent(in) :: view
      !> Adjoint of the field (npts, natom)
      real(wp), intent(in) :: w_phi(:, :)
      !> Adjoint of the tail charges
      real(wp), intent(in) :: w_qinf(:)
      !> Host response items
      type(response_type), intent(inout) :: response
      !> Error handling
      type(moist_error), allocatable, intent(out) :: error
      if (associated(self%sink)) self%sink = self%sink + w_phi
   end subroutine linear_field_adjoint_1d

   !> Form phi(r) = x at every grid point
   !>
   !> @param[in]  self     term
   !> @param[in]  grid     volume grid
   !> @param[in]  view     unbound view; unused
   !> @param[out] phi      field, Hartree/e (ngrid)
   !> @param[out] error    never set
   !> @param[out] dphi_dr  optional gradient, (3, ngrid)
   subroutine linear_field_3d(self, grid, view, phi, error, dphi_dr)
      !> Term
      class(linear_field_type), intent(in) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Unbound view
      type(coupling_view_type), intent(in) :: view
      !> Field (ngrid)
      real(wp), allocatable, intent(out) :: phi(:)
      !> Error handling
      type(moist_error), allocatable, intent(out) :: error
      !> Optional gradient (3, ngrid)
      real(wp), allocatable, intent(out), optional :: dphi_dr(:, :)
      phi = grid%xyz(1, :)
      if (present(dphi_dr)) then
         allocate (dphi_dr(3, grid%ngrid), source=0.0_wp)
         dphi_dr(1, :) = 1.0_wp
      end if
   end subroutine linear_field_3d

   !> Push w_phi d phi/d x into the grid adjoint
   !>
   !> @param[in,out] self          term
   !> @param[in]     grid          volume grid
   !> @param[in]     view          unbound view; unused
   !> @param[in]     w_phi         adjoint of the field (ngrid)
   !> @param[in]     w_qinf        adjoint of the tail charges; zero-sized, the term has none
   !> @param[in,out] response      host response items; untouched
   !> @param[out]    error         grid adjoint failure
   !> @param[in,out] grid_adjoint  optional volume adjoint accumulator
   !> @param[in,out] gradient      optional explicit nuclear gradient, Hartree/bohr (3, natom)
   subroutine linear_field_adjoint_3d(self, grid, view, w_phi, w_qinf, response, error, grid_adjoint, gradient)
      !> Term
      class(linear_field_type), intent(inout) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Unbound view
      type(coupling_view_type), intent(in) :: view
      !> Adjoint of the field (ngrid)
      real(wp), intent(in) :: w_phi(:)
      !> Adjoint of the tail charges
      real(wp), intent(in) :: w_qinf(:)
      !> Host response items
      type(response_type), intent(inout) :: response
      !> Error handling
      type(moist_error), allocatable, intent(out) :: error
      !> Optional volume adjoint accumulator
      type(volume_adjoint_type), intent(inout), optional :: grid_adjoint
      !> Optional explicit nuclear gradient, Hartree/bohr (3, natom)
      real(wp), intent(inout), optional :: gradient(:, :)
      real(wp), allocatable :: w_xyz(:, :)
      if (.not. present(grid_adjoint)) return
      allocate (w_xyz(3, grid%ngrid), source=0.0_wp)
      w_xyz(1, :) = w_phi
      call grid_adjoint%add_weights(error, w_xyz=w_xyz)
   end subroutine linear_field_adjoint_3d

end module test_moz_potential_kernel
