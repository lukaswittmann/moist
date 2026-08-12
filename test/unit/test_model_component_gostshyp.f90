!> Unit tests for the GOSTSHYP hydrostatic-pressure solvation-model component
!>
!> We are testing
!>   * the energy is `sum_i p_i gtilde_i` with `p_i = p a_i / ftilde_i`
!>   * the amplitudes handed back to the host reproduce that energy
!>   * the surface adjoints are the derivatives of that energy with respect to
!>     the cavity area, position and normal channels
!>   * a grid point that has left the density is switched off in *every*
!>     quantity, not just in the derivative
!>
module test_model_component_gostshyp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_type, only: coupling_type, potential_type
   use moist_model_components, only: solvation_model_component_gostshyp, new_component_gostshyp
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_drop, only: cavity_type_drop
   use test_model_component_helper, only: surface_fixture, &
      & new_surface_fixture, check_surface_weights, fixture_radial_normals, &
      & ngrid_sw => fixture_ngrid_param, sw_areas => fixture_areas_param, &
      & sw_xis => fixture_xis_param, sw_fs => fixture_fs_param, &
      & sw_xyz => fixture_xyz_param
   implicit none (type, external)
   private

   public :: collect_model_component_gostshyp

   !> Tolerance for values that must agree to roundoff
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)

   !* ------------------------- The model solute density ------------------------- *!

   !> Number of s-Gaussians in the model density
   integer, parameter :: nprim = 2
   !> Contraction coefficients of the model density
   real(wp), parameter :: rho_coeff(nprim) = [1.0_wp, 0.6_wp]
   !> Exponents of the model density, bohr**-2
   real(wp), parameter :: rho_expo(nprim) = [0.55_wp, 0.90_wp]
   !> Centers of the model density, bohr
   real(wp), parameter :: rho_center(3, nprim) = reshape([ &
                                                 0.10_wp, -0.05_wp, 0.08_wp, &
                                                 -0.25_wp, 0.18_wp, -0.12_wp], [3, nprim])

   !> A lot of pressure (in a.u.)
   real(wp), parameter :: test_pressure = 1.67_wp

contains

   !> Collect the GOSTSHYP component test suite.
   !>
   !> @param[out] testsuite Collected unit tests
   subroutine collect_model_component_gostshyp(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("gostshyp_energy_matches_amplitudes", test_gostshyp_energy), &
         & new_unittest("gostshyp_surface_weights", test_gostshyp_surface_weights), &
         & new_unittest("gostshyp_switching_factor_is_inert", test_gostshyp_w_f_zero), &
         & new_unittest("gostshyp_zero_pressure_short_circuit", test_gostshyp_short_circuit), &
         & new_unittest("gostshyp_lifecycle_guards", test_gostshyp_guards), &
         & new_unittest("gostshyp_unrepresentable_amplitudes", test_gostshyp_nonfinite) &
         & ]

   end subroutine collect_model_component_gostshyp

   !* ================================================================================= *!
   !*                        Closed-form moments of the model density                   *!
   !* ================================================================================= *!

   !> Gaussian moments of the model density about one grid point
   !>
   !> For a single primitive `c exp(-eta |r - A|^2)` against the point Gaussian
   !> `exp(-w |r - C|^2)`, the Gaussian product theorem gives an effective
   !> center `P = (eta A + w C)/(eta + w)` and, with `s = eta + w`,
   !> `d = P - C = eta (A - C)/s`,
   !>
   !>    S    = c (pi/s)^(3/2) exp(-eta w |A - C|^2 / s)
   !>    <G>  = S
   !>    <u>  = S d                        with u = r - C
   !>    <uu> = S (d d + I/(2 s))
   !>    <u |u|^2> = S d (5/(2 s) + |d|^2)
   !>
   !> The last follows from `u = v + d` with `v` isotropic about `P`: the odd
   !> moments of `v` vanish, `<v_a v_b> = delta_ab/(2 s)`, hence
   !> `<u|u|^2> = d/s + 3 d/(2 s) + d |d|^2`.
   !>
   !> @param[in]  center Grid-point center, bohr
   !> @param[in]  omega  Gaussian width of the grid point, bohr**-2
   !> @param[out] gt     Zeroth moment
   !> @param[out] pt     First moment (3)
   !> @param[out] mt     Second moment (3, 3)
   !> @param[out] rt     Third moment contracted to a vector (3)
   pure subroutine model_moments(center, omega, gt, pt, mt, rt)
      !> Grid-point center
      real(wp), intent(in) :: center(3)
      !> Gaussian width of the grid point
      real(wp), intent(in) :: omega
      !> Moments of the model density about `center`
      real(wp), intent(out) :: gt, pt(3), mt(3, 3), rt(3)

      !> Primitive index and cartesian indices
      integer :: iprim, a, b
      !> Combined exponent, overlap prefactor and displaced center
      real(wp) :: s, prefactor, offset(3), disp(3)

      gt = 0.0_wp
      pt = 0.0_wp
      mt = 0.0_wp
      rt = 0.0_wp

      do iprim = 1, nprim
         s = rho_expo(iprim) + omega
         offset = rho_center(:, iprim) - center
         disp = rho_expo(iprim)*offset/s
         prefactor = rho_coeff(iprim)*(pi/s)**1.5_wp &
            & *exp(-rho_expo(iprim)*omega*dot_product(offset, offset)/s)

         gt = gt + prefactor
         pt = pt + prefactor*disp
         do b = 1, 3
            do a = 1, 3
               mt(a, b) = mt(a, b) + prefactor*disp(a)*disp(b)
            end do
            mt(b, b) = mt(b, b) + prefactor/(2.0_wp*s)
         end do
         rt = rt + prefactor*disp*(2.5_wp/s + dot_product(disp, disp))
      end do

   end subroutine model_moments

   !> Fill the coupling moments for a whole grid
   !>
   !> @param[inout] coupling   Coupling data to fill
   !> @param[in]    centers    Grid-point centers (3, ngrid)
   !> @param[in]    areas      Grid-point areas (ngrid)
   subroutine set_model_moments(coupling, centers, areas)
      !> Coupling data to fill
      type(coupling_type), intent(inout) :: coupling
      !> Grid-point centers and areas
      real(wp), intent(in) :: centers(:, :), areas(:)

      !> Grid-point index and grid size
      integer :: igrid, ngrid

      ngrid = size(areas)
      call coupling%clear()
      allocate (coupling%gauss_gt(ngrid))
      allocate (coupling%gauss_pt(3, ngrid))
      allocate (coupling%gauss_mt(3, 3, ngrid))
      allocate (coupling%gauss_rt(3, ngrid))

      do igrid = 1, ngrid
         call model_moments(centers(:, igrid), gaussian_width(areas(igrid)), &
            & coupling%gauss_gt(igrid), coupling%gauss_pt(:, igrid), &
            & coupling%gauss_mt(:, :, igrid), coupling%gauss_rt(:, igrid))
      end do

   end subroutine set_model_moments

   !> The component's own width convention, `w = pi ln2 / a`
   !>
   !> @param[in] area Grid-point area, bohr**2
   !> @return Gaussian width, bohr**-2
   pure real(wp) function gaussian_width(area) result(omega)
      !> Grid-point area
      real(wp), intent(in) :: area

      omega = pi*log(2.0_wp)/area

   end function gaussian_width

   !> The GOSTSHYP energy of an arbitrary surface, moments rebuilt from scratch
   !>
   !> The finite-difference reference for the surface weights. It shares no
   !> state with the component and rebuilds the moments at the displaced
   !> surface, which is the whole point: the Gaussians live *on* the grid
   !> points, so every surface channel moves them.
   !>
   !> @param[in] areas   Grid-point areas
   !> @param[in] centers Grid-point centers
   !> @param[in] normals Grid-point outward normals
   !> @return GOSTSHYP energy
   function surface_energy(areas, centers, normals) result(energy)
      !> Grid-point areas, centers and normals
      real(wp), intent(in) :: areas(:), centers(:, :), normals(:, :)
      !> GOSTSHYP energy
      real(wp) :: energy

      !> Moments of the model density at one grid point
      real(wp) :: gt, pt(3), mt(3, 3), rt(3)
      !> Width and normal-projected gradient trace at one grid point
      real(wp) :: omega, ftilde
      !> Grid-point index
      integer :: igrid

      energy = 0.0_wp
      do igrid = 1, size(areas)
         omega = gaussian_width(areas(igrid))
         call model_moments(centers(:, igrid), omega, gt, pt, mt, rt)
         ftilde = -2.0_wp*omega*dot_product(normals(:, igrid), pt)
         energy = energy + test_pressure*areas(igrid)*gt/ftilde
      end do

   end function surface_energy

   !* ================================================================================= *!
   !*                                       Tests                                       *!
   !* ================================================================================= *!

   !> Build the synthetic surface every test in this module runs on.
   !>
   !> @param[out] cavity   Synthetic DROP surface
   !> @param[out] coupling Coupling data carrying the model-density moments
   !> @param[out] normals  Radial normal field of the fixture
   subroutine build_fixture_surface(cavity, coupling, normals)
      !> Synthetic DROP surface
      type(cavity_type_drop), intent(out) :: cavity
      !> Coupling data carrying the model-density moments
      type(coupling_type), intent(out) :: coupling
      !> Radial normal field of the fixture
      real(wp), intent(out) :: normals(3, ngrid_sw)

      normals = fixture_radial_normals()

      cavity%ngrid = ngrid_sw
      cavity%nsph = 1
      allocate (cavity%a, source=sw_areas)
      allocate (cavity%xi0, source=sw_xis)
      allocate (cavity%f, source=sw_fs)
      allocate (cavity%xyz, source=sw_xyz)
      allocate (cavity%normal0, source=normals)

      call set_model_moments(coupling, sw_xyz, sw_areas)

   end subroutine build_fixture_surface

   !> The energy, and the amplitudes the host would contract, agree
   !>
   !> The component reports an energy and, separately, the two amplitudes the
   !> host multiplies into its own integral blocks. Those are independent code
   !> paths through the same masking, so a mask applied in one and not the other
   !> shows up as a mismatch here rather than as an unclosable finite difference
   !> three layers up.
   !>
   !> @param[out] error Error handling
   subroutine test_gostshyp_energy(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Dummy structure; the component only stores it
      type(structure_type) :: mol
      !> Synthetic DROP surface and its coupling moments
      type(cavity_type_drop) :: cavity
      type(coupling_type) :: coupling
      !> Component under test
      type(solvation_model_component_gostshyp) :: component
      !> Potential accumulator receiving the amplitudes
      type(potential_type) :: potential
      !> Radial normal field
      real(wp) :: normals(3, ngrid_sw)
      !> Dummy molecular geometry
      real(wp) :: xyz_mol(3, 1)
      !> Reported and reconstructed energies
      real(wp) :: energy, rebuilt
      !> Independent reference energy
      real(wp) :: reference
      !> Grid-point index
      integer :: igrid

      call build_fixture_surface(cavity, coupling, normals)
      xyz_mol(:, 1) = 0.0_wp
      call new (mol, [1], xyz_mol)

      call new_component_gostshyp(component, test_pressure)
      call component%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP update failed: "//err%message)
         return
      end if

      energy = 0.0_wp
      call component%get_energy(coupling, cavity, energy, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP energy failed: "//err%message)
         return
      end if

      reference = surface_energy(sw_areas, sw_xyz, normals)
      call check(error, energy, reference, thr=thr, &
         & more="GOSTSHYP energy does not match the independent surface sum")
      if (allocated(error)) return

      call component%get_potential(coupling, cavity, potential, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP potential failed: "//err%message)
         return
      end if
      if (.not. allocated(potential%w_gauss_g) .or. .not. allocated(potential%w_gauss_f)) then
         call test_failed(error, "GOSTSHYP wrote no host amplitudes")
         return
      end if

      ! E = sum_i alpha_i gtilde_i, and beta_i = gtilde_i alpha_i / ftilde_i, so
      ! contracting the amplitudes against the traces the host would supply must
      ! return the energy: sum_i [w_gauss_g gtilde + w_gauss_f ftilde] = E - E = 0
      ! for the f channel alone, hence the two are checked separately.
      rebuilt = dot_product(potential%w_gauss_g, coupling%gauss_gt)
      call check(error, rebuilt, energy, thr=thr, &
         & more="GOSTSHYP amplitudes do not reproduce the reported energy")
      if (allocated(error)) return

      do igrid = 1, ngrid_sw
         call check(error, potential%w_gauss_f(igrid) < 0.0_wp, &
            & more="GOSTSHYP f-amplitude lost its sign fold")
         if (allocated(error)) return
      end do

   end subroutine test_gostshyp_energy

   !> The GOSTSHYP surface adjoints against fourth-order central differences
   !>
   !> Driven on the shared synthetic seven-point DROP surface with the radial
   !> normal field: the amplitude `p a / ftilde` inverts the normal-projected
   !> trace, so a field tilted away from the density would put `n . Pt` near
   !> zero at some point and cost the comparison its conditioning, not its
   !> correctness. All four channels stay enabled -- the area reaches the energy
   !> twice over, explicitly and through the Gaussian width, and the harness
   !> folds `w_a` into both scalar channels.
   !>
   !> @param[out] error Error handling
   subroutine test_gostshyp_surface_weights(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Dummy structure; the component only stores it
      type(structure_type) :: mol
      !> Synthetic DROP surface and its coupling moments
      type(cavity_type_drop) :: cavity
      type(coupling_type) :: coupling
      !> Component under test
      type(solvation_model_component_gostshyp) :: component
      !> Analytic surface weights
      type(cavity_surface_adjoint_type) :: weights
      !> Fixture mirroring the synthetic surface for the harness
      type(surface_fixture) :: surface
      !> Radial normal field
      real(wp) :: normals(3, ngrid_sw)
      !> Dummy molecular geometry
      real(wp) :: xyz_mol(3, 1)

      !> Finite-difference displacement.  Measured to be the optimum: the worst
      !> stencil deviation over all four channels is 3.8e-12 here, against
      !> 6.4e-10 at 1e-3 (truncation) and 8.2e-11 at 1e-5 (roundoff).
      real(wp), parameter :: step = 3.0e-4_wp
      !> Tolerances of the finite-difference comparison.
      !>
      !> Set from what the stencil achieves rather than from what the energy
      !> looks like it should need. Every channel agrees to 4e-12 absolute at
      !> the step above, so `fd_atol` sits ~50x clear of the worst observed
      !> deviation -- room for a compiler that contracts differently, not for a
      !> wrong derivative.
      !>
      !> The pair is what fixes the sensitivity of this test, so it is chosen
      !> against a target: scaling any single weight by 1 + 1e-8 must fail here.
      !> `check` thresholds at `max(fd_atol, fd_rtol*expected)`, so detection is
      !> driven by the largest entry of a channel -- 1.18 for `w_n`, giving a
      !> 1.2e-8 discrepancy against a 1.2e-9 threshold. Loosening `fd_rtol` by
      !> one decade is enough to stop seeing that, which is where these values
      !> started; anything relying on the tolerances should measure first.
      real(wp), parameter :: fd_atol = 2.0e-10_wp
      real(wp), parameter :: fd_rtol = 1.0e-9_wp

      call build_fixture_surface(cavity, coupling, normals)
      xyz_mol(:, 1) = 0.0_wp
      call new (mol, [1], xyz_mol)

      call new_component_gostshyp(component, test_pressure)
      call component%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP update failed: "//err%message)
         return
      end if

      call weights%init(ngrid_sw)
      call component%get_surface_weights(coupling, cavity, weights, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP surface-weight assembly failed: "//err%message)
         return
      end if

      ! Vacuous unless every channel the harness checks actually carries signal.
      call check(error, maxval(abs(weights%w_a)) > 0.0_wp, &
         & more="GOSTSHYP wrote no area weights")
      if (allocated(error)) return
      call check(error, maxval(abs(weights%w_xyz)) > 0.0_wp, &
         & more="GOSTSHYP wrote no position weights")
      if (allocated(error)) return
      call check(error, maxval(abs(weights%w_n)) > 0.0_wp, &
         & more="GOSTSHYP wrote no normal weights")
      if (allocated(error)) return

      call new_surface_fixture(surface, sw_areas, sw_xis, sw_fs, sw_xyz, normals)
      call check_surface_weights(error, surface, weights, gostshyp_surface_energy, &
         & "gostshyp", step=step, thr_abs=fd_atol, thr_rel=fd_rtol)

   contains

      !> The energy of a perturbed surface fixture
      !>
      !> @param[in] trial Perturbed surface fixture
      !> @return GOSTSHYP energy of the perturbed surface
      function gostshyp_surface_energy(trial) result(energy)

         !> Perturbed surface fixture
         type(surface_fixture), intent(in) :: trial

         !> GOSTSHYP energy of the perturbed surface
         real(wp) :: energy

         energy = surface_energy(trial%areas(), trial%xyz, trial%normal)

      end function gostshyp_surface_energy

   end subroutine test_gostshyp_surface_weights

   !> The switching factor carries no GOSTSHYP dependence of its own
   !>
   !> `w_f` must be *exactly* zero, not merely small: the area is the only route
   !> by which the switching factor reaches this energy, and the cavity folds
   !> that route itself. A nonzero `w_f` would double-count it. Pinned here
   !> rather than left to a comment because the harness cannot tell a genuine
   !> `w_f` from an area contribution that leaked into it.
   !>
   !> @param[out] error Error handling
   subroutine test_gostshyp_w_f_zero(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Dummy structure; the component only stores it
      type(structure_type) :: mol
      !> Synthetic DROP surface and its coupling moments
      type(cavity_type_drop) :: cavity
      type(coupling_type) :: coupling
      !> Component under test
      type(solvation_model_component_gostshyp) :: component
      !> Analytic surface weights
      type(cavity_surface_adjoint_type) :: weights
      !> Radial normal field
      real(wp) :: normals(3, ngrid_sw)
      !> Dummy molecular geometry
      real(wp) :: xyz_mol(3, 1)

      call build_fixture_surface(cavity, coupling, normals)
      xyz_mol(:, 1) = 0.0_wp
      call new (mol, [1], xyz_mol)

      call new_component_gostshyp(component, test_pressure)
      call component%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP update failed: "//err%message)
         return
      end if

      call weights%init(ngrid_sw)
      call component%get_surface_weights(coupling, cavity, weights, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP surface-weight assembly failed: "//err%message)
         return
      end if

      call check(error, maxval(abs(weights%w_f)), 0.0_wp, thr=0.0_wp, &
         & more="GOSTSHYP wrote a switching-factor weight")
      if (allocated(error)) return
      call check(error, maxval(abs(weights%w_k1)) + maxval(abs(weights%w_k2)), 0.0_wp, &
         & thr=0.0_wp, more="GOSTSHYP wrote a curvature weight")

   end subroutine test_gostshyp_w_f_zero

   !> A zero pressure must short-circuit before the moments are ever read
   !>
   !> Driven with no moments supplied at all, so a missing short circuit is
   !> observable as the error the component raises when it looks for them.
   !>
   !> @param[out] error Error handling
   subroutine test_gostshyp_short_circuit(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Dummy structure; the component only stores it
      type(structure_type) :: mol
      !> Synthetic DROP surface, carrying no coupling moments
      type(cavity_type_drop) :: cavity
      type(coupling_type) :: coupling
      !> Component under test
      type(solvation_model_component_gostshyp) :: component
      !> Radial normal field
      real(wp) :: normals(3, ngrid_sw)
      !> Dummy molecular geometry
      real(wp) :: xyz_mol(3, 1)

      call build_fixture_surface(cavity, coupling, normals)
      call coupling%clear()
      xyz_mol(:, 1) = 0.0_wp
      call new (mol, [1], xyz_mol)

      call new_component_gostshyp(component, 0.0_wp)
      call check_inert(error, component, coupling, cavity, mol, "zero pressure")
      if (allocated(error)) return

      ! A component scaled to zero contributes nothing either, and must reach
      ! that conclusion without asking for moments it will not use.
      call new_component_gostshyp(component, test_pressure)
      component%scale = 0.0_wp
      call check_inert(error, component, coupling, cavity, mol, "zero scale")
      if (allocated(error)) return

   end subroutine test_gostshyp_short_circuit

   !> A component that cannot contribute must touch nothing at all
   !>
   !> Driven with the coupling cleared, so a missing short circuit shows up as
   !> the error the component raises when it looks for its moments.
   !>
   !> @param[out]   error     Test failure information
   !> @param[inout] component Component under test
   !> @param[in]    coupling  Coupling data carrying no moments
   !> @param[in]    cavity    Synthetic DROP surface
   !> @param[in]    mol       Molecular structure
   !> @param[in]    label     Case name used in failure messages
   subroutine check_inert(error, component, coupling, cavity, mol, label)
      !> Test failure information
      type(error_type), allocatable, intent(out) :: error
      !> Component under test
      type(solvation_model_component_gostshyp), intent(inout) :: component
      !> Coupling data carrying no moments
      type(coupling_type), intent(in) :: coupling
      !> Synthetic DROP surface.  Not `intent(in)`: the component's own energy
      !> and potential hooks take the cavity as `intent(inout)`.
      type(cavity_type_drop), intent(inout) :: cavity
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Case name used in failure messages
      character(len=*), intent(in) :: label

      type(moist_error_type), allocatable :: err
      !> Prefilled accumulator the component must not touch
      type(cavity_surface_adjoint_type) :: prefilled
      !> Potential accumulator the component must not touch
      type(potential_type) :: potential
      !> Energy accumulator carrying a sentinel
      real(wp) :: energy

      !> Value prefilled into every channel of the untouched accumulator
      real(wp), parameter :: prefill = 0.75_wp
      !> Sentinel the energy must keep
      real(wp), parameter :: sentinel = 1.5_wp

      call component%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP update failed: "//err%message)
         return
      end if

      energy = sentinel
      call component%get_energy(coupling, cavity, energy, err)
      call check(error, .not. allocated(err), &
         & more="GOSTSHYP read the moments at "//label)
      if (allocated(error)) return
      call check(error, energy, sentinel, thr=0.0_wp, &
         & more="GOSTSHYP moved the energy at "//label)
      if (allocated(error)) return

      call component%get_potential(coupling, cavity, potential, err)
      call check(error, .not. allocated(err), &
         & more="GOSTSHYP potential failed at "//label)
      if (allocated(error)) return
      call check(error, .not. allocated(potential%w_gauss_g), &
         & more="GOSTSHYP allocated host amplitudes at "//label)
      if (allocated(error)) return

      call prefilled%init(ngrid_sw)
      prefilled%w_a = prefill
      prefilled%w_xyz = prefill
      prefilled%w_n = prefill
      call component%get_surface_weights(coupling, cavity, prefilled, err)
      call check(error, .not. allocated(err), &
         & more="GOSTSHYP surface weights failed at "//label)
      if (allocated(error)) return
      call check(error, maxval(abs(prefilled%w_a - prefill)), 0.0_wp, thr=0.0_wp, &
         & more="GOSTSHYP wrote area weights at "//label)
      if (allocated(error)) return
      call check(error, maxval(abs(prefilled%w_xyz - prefill)), 0.0_wp, thr=0.0_wp, &
         & more="GOSTSHYP wrote position weights at "//label)

   end subroutine check_inert

   !> Missing or mis-shaped host moments are refused rather than assumed
   !>
   !> The moments are rebuilt by the host after every cavity update, so a stale
   !> set is the natural failure. It must be loud: silently reusing moments from
   !> the previous geometry would produce a plausible energy for a surface that
   !> no longer exists.
   !>
   !> @param[out] error Error handling
   subroutine test_gostshyp_guards(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Dummy structure; the component only stores it
      type(structure_type) :: mol
      !> Synthetic DROP surface and its coupling moments
      type(cavity_type_drop) :: cavity
      type(coupling_type) :: coupling
      !> Component under test
      type(solvation_model_component_gostshyp) :: component
      !> Cavity missing its surface arrays
      type(cavity_type_drop) :: bare
      !> Radial normal field
      real(wp) :: normals(3, ngrid_sw)
      !> Dummy molecular geometry
      real(wp) :: xyz_mol(3, 1)
      !> Energy accumulator
      real(wp) :: energy

      call build_fixture_surface(cavity, coupling, normals)
      xyz_mol(:, 1) = 0.0_wp
      call new (mol, [1], xyz_mol)

      call new_component_gostshyp(component, test_pressure)

      ! An un-updated cavity carries no surface arrays at all.
      bare%ngrid = ngrid_sw
      bare%nsph = 1
      call component%update(mol, bare, err)
      call check(error, allocated(err), &
         & more="GOSTSHYP accepted a cavity without surface data")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)

      call component%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP update failed: "//err%message)
         return
      end if

      ! Moments never supplied.
      call coupling%clear()
      energy = 0.0_wp
      call component%get_energy(coupling, cavity, energy, err)
      call check(error, allocated(err), &
         & more="GOSTSHYP accepted a coupling without Gaussian moments")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)

      ! Moments supplied for a different grid size: the stale-cavity case.
      call set_model_moments(coupling, sw_xyz(:, 1:ngrid_sw - 1), sw_areas(1:ngrid_sw - 1))
      energy = 0.0_wp
      call component%get_energy(coupling, cavity, energy, err)
      call check(error, allocated(err), &
         & more="GOSTSHYP accepted Gaussian moments of the wrong grid size")
      if (allocated(error)) return
      if (allocated(err)) deallocate (err)

      ! And the forward-mode gradient is refused rather than silently zero.
      block
         !> Nuclear-gradient accumulator
         real(wp) :: gradient(3, 1)

         gradient = 0.0_wp
         call set_model_moments(coupling, sw_xyz, sw_areas)
         call component%get_gradient(coupling, cavity, gradient, err)
         call check(error, allocated(err), &
            & more="GOSTSHYP returned a forward-mode gradient instead of refusing")
         if (allocated(error)) return
         if (allocated(err)) deallocate (err)
         call check(error, maxval(abs(gradient)), 0.0_wp, thr=0.0_wp, &
            & more="GOSTSHYP wrote a forward-mode gradient")
      end block

   end subroutine test_gostshyp_guards

   !> Moments that divide their way out of the reals switch their point off
   !>
   !> The activity floor is *relative*, so it can only compare grid points with
   !> one another -- a moment supply that is wrong by a uniform factor passes it
   !> untouched and only fails at the divisions. Here `gauss_pt` alone is scaled
   !> into the far underflow, which leaves `alpha = p a / ftilde` large but
   !> perfectly finite and sends `beta = gtilde alpha / ftilde` past the largest
   !> double: the energy would look plausible while the amplitude the host folds
   !> into its Fock matrix is an infinity. Both are pinned, because only one of
   !> them shows the damage.
   !>
   !> @param[out] error Error handling
   subroutine test_gostshyp_nonfinite(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: err

      !> Dummy structure; the component only stores it
      type(structure_type) :: mol
      !> Synthetic DROP surface and its coupling moments
      type(cavity_type_drop) :: cavity
      type(coupling_type) :: coupling
      !> Component under test
      type(solvation_model_component_gostshyp) :: component
      !> Potential accumulator receiving the amplitudes
      type(potential_type) :: potential
      !> Radial normal field
      real(wp) :: normals(3, ngrid_sw)
      !> Dummy molecular geometry
      real(wp) :: xyz_mol(3, 1)
      !> Energy accumulator
      real(wp) :: energy

      call build_fixture_surface(cavity, coupling, normals)
      xyz_mol(:, 1) = 0.0_wp
      call new (mol, [1], xyz_mol)

      call new_component_gostshyp(component, test_pressure)
      call component%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP update failed: "//err%message)
         return
      end if

      !> Uniform, so every point keeps its share of the total and the relative
      !> floor has nothing to bite on.
      coupling%gauss_pt = coupling%gauss_pt*1.0e-300_wp

      energy = 0.0_wp
      call component%get_energy(coupling, cavity, energy, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP energy failed: "//err%message)
         return
      end if
      call check(error, ieee_is_finite(energy), &
         & more="GOSTSHYP reported a non-finite energy")
      if (allocated(error)) return
      call check(error, energy, 0.0_wp, thr=0.0_wp, &
         & more="GOSTSHYP kept an unrepresentable grid point in the energy")
      if (allocated(error)) return

      call component%get_potential(coupling, cavity, potential, err)
      if (allocated(err)) then
         call test_failed(error, "GOSTSHYP potential failed: "//err%message)
         return
      end if
      call check(error, all(ieee_is_finite(potential%w_gauss_g)) &
         & .and. all(ieee_is_finite(potential%w_gauss_f)), &
         & more="GOSTSHYP handed the host a non-finite amplitude")
      if (allocated(error)) return
      call check(error, maxval(abs(potential%w_gauss_g)) + maxval(abs(potential%w_gauss_f)), &
         & 0.0_wp, thr=0.0_wp, &
         & more="GOSTSHYP kept an unrepresentable grid point in the amplitudes")

   end subroutine test_gostshyp_nonfinite

end module test_model_component_gostshyp
