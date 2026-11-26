! Shared fixture modules of the unit tests, none of which registers a suite
!
! - test_moz_fixtures                  MOZ tests
! - test_cavity_symmetry_fixtures      cavity symmetry tests
! - test_model_component_pcm_fixtures  PCM component tests

! ============================================================================================
! MOZ fixtures
! ============================================================================================

!> Shared fixtures of the MOZ tests (not a suite)
!>
!> - Table water solvent with custom LJ and fixed charges
!> - Pair-only term with switchable declaration failure
!> - Reference 12-6 potential with Lorentz-Berthelot mixing
!> - Answers to pending atomic_charges requests
module test_moz_fixtures
   use mctc_env, only: wp, moist_error => error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_channels_coupling, only: coupling_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_solvent_system, only: solvent_system_type, new_solvent_system
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent
   use moist_model_moz_potential_term, only: potential_term_type, potential_name_len
   use moist_model_moz_potential_lj_base, only: new_lj_12_6
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use moist_model_moz_potential_electrostatic_fixed, only: fixed_charges_type, new_fixed_charges
   implicit none(type, external)
   private

   public :: water_id, water_sigma, water_epsilon, water_charges
   public :: new_water_system, new_water_solvent, ref_lj, answer_charges

   !> Water solvent id of the table
   integer, parameter :: water_id = 175
   !> Distinct per-site sigma, bohr; swapped pair columns detectable
   real(wp), parameter :: water_sigma(3) = [6.0_wp, 1.5_wp, 2.0_wp]
   !> Epsilon per site, Hartree
   real(wp), parameter :: water_epsilon(3) = [2.4e-4_wp, 1.0e-5_wp, 2.0e-5_wp]
   !> Charges per site, e
   real(wp), parameter :: water_charges(3) = [-0.8_wp, 0.3_wp, 0.5_wp]

contains

   !> Construct table water system
   !>
   !> @param[out] system  water solvent system, structure in bohr
   !> @param[out] err     table failure
   subroutine new_water_system(system, err)
      !> Water solvent system
      type(solvent_system_type), intent(out) :: system
      !> Table failure
      type(moist_error), allocatable, intent(out) :: err
      call new_solvent_system(system, water_id, error=err)
   end subroutine new_water_system

   !> Updated table water solvent with custom LJ and fixed charges
   !>
   !> @param[in]  ctx      run context, outlives the solvent
   !> @param[out] solvent  updated solvent
   !> @param[out] err      construction or update error
   !> @param[in]  mixing   optional Lennard-Jones mixing rule; Lorentz-Berthelot when absent
   subroutine new_water_solvent(ctx, solvent, err, mixing)
      !> Run context, outlives the solvent
      type(moist_context_type), intent(in), target :: ctx
      !> Updated solvent
      type(solvent_vv_type), intent(out) :: solvent
      !> Construction or update error
      type(moist_error), allocatable, intent(out) :: err
      !> Optional Lennard-Jones mixing rule
      integer, intent(in), optional :: mixing
      type(solvent_system_type) :: system
      type(custom_lj_type) :: lj
      type(fixed_charges_type) :: charges
      call new_water_system(system, err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      if (present(mixing) .and. .not. allocated(err)) call solvent%potential%set_mixing(mixing, err)
      if (.not. allocated(err)) call new_custom_lj(lj, water_sigma, water_epsilon, err)
      if (.not. allocated(err)) call new_fixed_charges(charges, water_charges, err)
      if (.not. allocated(err)) call solvent%potential%add(lj, err)
      if (.not. allocated(err)) call solvent%potential%add(charges, err)
      if (.not. allocated(err)) call solvent%update(err)
   end subroutine new_water_solvent

   !> Reference 12-6 potential with Lorentz-Berthelot mixing
   !>
   !> @param[in] r   distance, bohr
   !> @param[in] e1  epsilon of site 1, Hartree
   !> @param[in] s1  sigma of site 1, bohr
   !> @param[in] e2  epsilon of site 2, Hartree
   !> @param[in] s2  sigma of site 2, bohr
   elemental function ref_lj(r, e1, s1, e2, s2) result(u)
      !> Distance, bohr
      real(wp), intent(in) :: r
      !> Epsilon of site 1
      real(wp), intent(in) :: e1
      !> Sigma of site 1
      real(wp), intent(in) :: s1
      !> Epsilon of site 2
      real(wp), intent(in) :: e2
      !> Sigma of site 2
      real(wp), intent(in) :: s2
      !> Potential, Hartree
      real(wp) :: u
      real(wp) :: sr6
      sr6 = (0.5_wp*(s1 + s2)/r)**6
      u = 4.0_wp*sqrt(e1*e2)*(sr6*sr6 - sr6)
   end function ref_lj

   !> Answer the pending atomic charges of a staged walk
   !>
   !> @param[in,out] coupling  staged coupling
   !> @param[in]     q         charges per atom
   !> @param[out]    err       missing request or rejected answer
   subroutine answer_charges(coupling, q, err)
      !> Staged coupling
      type(coupling_type), intent(inout) :: coupling
      !> Charges per atom
      real(wp), intent(in) :: q(:)
      !> Missing request or rejected answer
      type(moist_error), allocatable, intent(out) :: err
      do while (coupling%next())
         associate (item => coupling%request())
            if (item%name() /= "atomic_charges") cycle
         end associate
         call coupling%answer("q", q, err)
         return
      end do
      call fatal_error(err, "atomic_charges is not pending")
   end subroutine answer_charges

end module test_moz_fixtures

! ============================================================================================
! Cavity symmetry fixtures
! ============================================================================================

!> Shared fixtures of the cavity symmetry tests
!>
!> - reference structures in bohr, exactly symmetric unless distorted
!> - rigid motions: rotation, translation and atom permutation of a structure,
!>   deterministic rotations (large ones included) and permutations
!> - surface points of any cavity in a comparable form, compared as sets
!>   after mapping the moved cavity back, optionally modulo a symmetry of the
!>   molecule (tied frame choices give the image of the cavity under one);
!>   NumSA, which stores no surface, gives its quadrature nodes
!> - covariance of nuclear gradients, total and per sphere
!>
!> Used by `test_cavity_symmetry`; meant for every Lebedev-based cavity
module test_cavity_symmetry_fixtures
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use moist_math_linalg_geometry, only: identity3
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use mctc_io_convert, only: aatoau
   use mstore, only: get_structure
   use moist_cavity_type, only: cavity_type
   use moist_cavity_numsa, only: cavity_type_numsa
   use moist_cavity_grid_model_symmetry_type, only: grid_model_symmetry_type
   use moist_cavity_grid_model_type, only: grid_components_type
   use moist_math_symmetry_msym, only: msym_group_type, find_point_group
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: symmetry_of, model_nrot, model_rot_index, model_rotation, model_rotate_grid
   public :: model_add_torque, model_frame_jacobian
   public :: structure_of, new_tilted, rotation_matrix, test_rotation, test_permutation, rigid_motion
   public :: surface_points_type, surface_points_of, match_surface_points, molecule_symmetry
   public :: same_bits, gradient_covariance_error, sphere_gradient_covariance_error

   !> Whether two real arrays agree in shape and bit for bit
   interface same_bits
      module procedure same_bits_1, same_bits_2, same_bits_3, same_bits_4
   end interface same_bits

   !> Number of scalar values compared per point: area, switching, width, volume
   integer, parameter :: nval = 4
   !> Their names in messages
   character(len=*), parameter :: val_name(nval) = [character(len=9) :: "area", "switching", &
      & "width", "volume"]

   !> Surface points of a cavity in a comparable form
   type :: surface_points_type
      !> Number of points
      integer :: n = 0
      !> Positions (3, n), bohr
      real(wp), allocatable :: xyz(:, :)
      !> Outward unit normals (3, n)
      real(wp), allocatable :: normal(:, :)
      !> Owner sphere (atom) of every point (n)
      integer, allocatable :: owner(:)
      !> Area, switching, Gaussian width and volume element of every point
      !>
      !> - shape (nval, n); zero where the cavity has no such array
      real(wp), allocatable :: val(:, :)
   end type surface_points_type

contains

   !> Symmetry model a cavity holds as its grid model
   !>
   !> Null unless the cavity's grid model is a [[grid_model_symmetry_type]]; valid while
   !> the cavity keeps its model
   !>
   !> @param[in] cav Cavity constructed with a symmetry model
   function symmetry_of(cav) result(sym)
      !> Cavity
      class(cavity_type), intent(in), target :: cav
      !> Its symmetry model
      type(grid_model_symmetry_type), pointer :: sym

      sym => null()
      if (.not. allocated(cav%grid_model)) return
      select type (model => cav%grid_model)
      type is (grid_model_symmetry_type)
         sym => model
      end select
   end function symmetry_of

   !> Components of a symmetry model, for one sphere when nothing is computed
   !>
   !> @param[in]  model Symmetry model
   !> @param[out] comps Its components
   subroutine model_components(model, comps)
      !> Symmetry model
      type(grid_model_symmetry_type), intent(in) :: model
      !> Components
      type(grid_components_type), intent(out) :: comps

      integer :: nsph

      nsph = 1
      if (model%computed) nsph = size(model%xyz, 2)
      call model%components(nsph, comps)
   end subroutine model_components

   !> Number of distinct rotations of a symmetry model
   !>
   !> @param[in] model Symmetry model
   function model_nrot(model) result(nrot)
      !> Symmetry model
      type(grid_model_symmetry_type), intent(in) :: model
      !> Number of distinct rotations
      integer :: nrot

      type(grid_components_type) :: comps

      call model_components(model, comps)
      nrot = size(comps%rotation, 3)
   end function model_nrot

   !> Distinct rotation of the first component of a sphere
   !>
   !> @param[in] model Symmetry model
   !> @param[in] isph  Sphere
   function model_rot_index(model, isph) result(irot)
      !> Symmetry model
      type(grid_model_symmetry_type), intent(in) :: model
      !> Sphere
      integer, intent(in) :: isph
      !> Rotation index
      integer :: irot

      type(grid_components_type) :: comps

      call model_components(model, comps)
      irot = 1
      if (model%computed) irot = comps%rot(comps%sph_first(isph))
   end function model_rot_index

   !> Rotation of the first component of a sphere
   !>
   !> @param[in] model Symmetry model
   !> @param[in] isph  Sphere
   function model_rotation(model, isph) result(rot)
      !> Symmetry model
      type(grid_model_symmetry_type), intent(in) :: model
      !> Sphere
      integer, intent(in) :: isph
      !> Rotation
      real(wp) :: rot(3, 3)

      type(grid_components_type) :: comps

      call model_components(model, comps)
      rot = comps%rotation(:, :, model_rot_index(model, isph))
   end function model_rotation

   !> Copies of a unit grid under every distinct rotation of a symmetry model
   !>
   !> A plain model copies, as the grid container does
   !>
   !> @param[in]  model     Symmetry model
   !> @param[in]  unit_grid Unit grid (3, npt)
   !> @param[out] rotated   Rotated copies (3, npt, nrot)
   subroutine model_rotate_grid(model, unit_grid, rotated)
      !> Symmetry model
      type(grid_model_symmetry_type), intent(in) :: model
      !> Unit grid
      real(wp), intent(in) :: unit_grid(:, :)
      !> Rotated copies
      real(wp), allocatable, intent(out) :: rotated(:, :, :)

      type(grid_components_type) :: comps
      integer :: irot

      call model_components(model, comps)
      allocate (rotated(3, size(unit_grid, 2), size(comps%rotation, 3)))
      if (model%is_plain()) then
         rotated(:, :, 1) = unit_grid
         return
      end if
      do irot = 1, size(rotated, 3)
         rotated(:, :, irot) = matmul(comps%rotation(:, :, irot), unit_grid)
      end do
   end subroutine model_rotate_grid

   !> Add the gradient of component torques at constant weights
   !>
   !> @param[in]     model    Symmetry model
   !> @param[in]     torque   Torque per component (3, ncomp)
   !> @param[in,out] gradient Nuclear gradient (3, nat), accumulated
   subroutine model_add_torque(model, torque, gradient)
      !> Symmetry model
      type(grid_model_symmetry_type), intent(in) :: model
      !> Torque per component
      real(wp), intent(in) :: torque(:, :)
      !> Nuclear gradient
      real(wp), intent(inout) :: gradient(:, :)

      real(wp), allocatable :: wbar(:)

      allocate (wbar(size(torque, 2)), source=0.0_wp)
      call model%add_gradient(torque, wbar, gradient)
   end subroutine model_add_torque

   !> Angular velocity of every distinct rotation per sphere displacement
   !>
   !> @param[in]  model Symmetry model
   !> @param[out] jac   Jacobians (3, 3, nat, nrot)
   subroutine model_frame_jacobian(model, jac)
      !> Symmetry model
      type(grid_model_symmetry_type), intent(in) :: model
      !> Jacobians
      real(wp), allocatable, intent(out) :: jac(:, :, :, :)

      real(wp), allocatable :: jac_weight(:, :, :)
      integer, allocatable :: weight_class(:)

      call model%jacobian(jac, jac_weight, weight_class)
   end subroutine model_frame_jacobian

   !* ================================================================================= *!
   !*                                    Structures                                     *!
   !* ================================================================================= *!

   !> Reference structures in bohr, exactly symmetric unless distorted
   !>
   !> water, water_distorted, ammonia, ammonia_distorted, methane, benzene,
   !> benzene_noisy, sf6, hocl, c1, c1_plane, ch3f, ethane, mb16, bf3 (default);
   !> at a generic pose: co2, hcn, co2_bent (nearly a line), ring13,
   !> h2o2 (C2), h2o2_trans (C2h), ci,
   !> ethene_twisted (D2), s4, allene (D2d), pyramid (C4v), c6v, pyrite (Th),
   !> icosahedron (Ih), t12 (T), o24 (O), cube (Oh); c1_sphere and ci_sphere:
   !> c1 and ci for a sphere-level frame
   !>
   !> @param[in]  name Molecule
   !> @param[out] mol  Structure
   recursive subroutine structure_of(name, mol)
      !> Molecule
      character(len=*), intent(in) :: name
      !> Structure
      type(structure_type), intent(out) :: mol

      real(wp), allocatable :: xyz(:, :)
      real(wp) :: phi, a
      integer :: k

      select case (name)
      case ("water", "water_distorted")
         xyz = reshape([0.0_wp, 0.0_wp, 0.1173_wp, 0.0_wp, 0.7572_wp, -0.4692_wp, &
            & 0.0_wp, -0.7572_wp, -0.4692_wp], [3, 3])*aatoau
         ! An in-plane stretch of one O-H bond keeps only the molecular plane
         if (name == "water_distorted") xyz(2, 2) = xyz(2, 2) + 1.0e-5_wp
         call new(mol, [8, 1, 1], xyz)
      case ("ammonia", "ammonia_distorted")
         allocate (xyz(3, 4))
         xyz(:, 1) = 0.0_wp
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k + 1) = [0.9377_wp*cos(phi), 0.9377_wp*sin(phi), -0.3816_wp]
         end do
         xyz = xyz*aatoau
         ! Off C3v by 1e-4 bohr: a plane-level frame without ties
         if (name == "ammonia_distorted") xyz(:, 2) = xyz(:, 2) + [1.0e-4_wp, 0.4e-4_wp, 0.0_wp]
         call new(mol, [7, 1, 1, 1], xyz)
      case ("methane")
         a = 0.629_wp*aatoau
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, a, a, a, a, -a, -a, -a, a, -a, -a, -a, a], [3, 5])
         call new(mol, [6, 1, 1, 1, 1], xyz)
      case ("benzene", "benzene_noisy")
         allocate (xyz(3, 12))
         do k = 1, 6
            phi = 2.0_wp*pi*real(k - 1, wp)/6.0_wp
            xyz(:, k) = [1.397_wp*cos(phi), 1.397_wp*sin(phi), 0.0_wp]*aatoau
            xyz(:, k + 6) = [2.481_wp*cos(phi), 2.481_wp*sin(phi), 0.0_wp]*aatoau
         end do
         ! Every atom off by up to 2e-5 bohr: D6h only to libmsym's thresholds
         if (name == "benzene_noisy") then
            do k = 1, 12
               xyz(:, k) = xyz(:, k) + 2.0e-5_wp*sin(1.7_wp*real(k, wp) + [0.3_wp, 1.9_wp, 4.1_wp])
            end do
         end if
         call new(mol, [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1], xyz)
      case ("sf6")
         a = 1.56_wp*aatoau
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, a, 0.0_wp, 0.0_wp, -a, 0.0_wp, 0.0_wp, &
            & 0.0_wp, a, 0.0_wp, 0.0_wp, -a, 0.0_wp, 0.0_wp, 0.0_wp, a, 0.0_wp, 0.0_wp, -a], [3, 7])
         call new(mol, [16, 9, 9, 9, 9, 9, 9], xyz)
      case ("hocl")
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.966_wp, 0.0_wp, 0.0_wp, -0.383_wp, 1.643_wp, 0.0_wp], &
            & [3, 3])*aatoau
         call new(mol, [8, 1, 17], xyz)
      case ("c1")
         ! Bromochlorofluoromethane-like, no symmetry
         xyz = reshape([0.02_wp, -0.01_wp, 0.03_wp, 0.65_wp, 0.71_wp, 0.48_wp, -0.81_wp, 0.52_wp, &
            & -0.69_wp, 0.93_wp, -1.24_wp, -0.61_wp, -1.12_wp, -0.88_wp, 1.37_wp], [3, 5])*aatoau
         call new(mol, [6, 1, 9, 17, 35], xyz)
      case ("c1_plane")
         ! Three-bladed propeller without a mirror plane, one atom moved by
         ! about 1e-5 bohr: no symmetry, two eigenvalues of the frame tensor
         ! agree to far below the gap threshold, a plane-level frame
         allocate (xyz(3, 7), source=0.0_wp)
         xyz(:, 1) = [0.0_wp, 0.0_wp, 0.21_wp]
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k + 1) = [1.31_wp*cos(phi), 1.31_wp*sin(phi), -0.42_wp]
            xyz(:, k + 4) = [0.83_wp*cos(phi + 0.7_wp), 0.83_wp*sin(phi + 0.7_wp), 0.95_wp]
         end do
         xyz = xyz*aatoau
         xyz(:, 2) = xyz(:, 2) + [1.0e-5_wp, -0.6e-5_wp, 0.8e-5_wp]
         call new(mol, [6, 9, 9, 9, 1, 1, 1], xyz)
      case ("ch3f")
         ! C3v: carbon and fluorine on the axis
         allocate (xyz(3, 5), source=0.0_wp)
         xyz(:, 2) = [0.0_wp, 0.0_wp, 1.383_wp]
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k + 2) = [1.029_wp*cos(phi), 1.029_wp*sin(phi), -0.350_wp]
         end do
         call new(mol, [6, 9, 1, 1, 1], xyz*aatoau)
      case ("ethane")
         ! Staggered, D3d
         allocate (xyz(3, 8), source=0.0_wp)
         xyz(:, 1) = [0.0_wp, 0.0_wp, 0.7655_wp]
         xyz(:, 2) = [0.0_wp, 0.0_wp, -0.7655_wp]
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k + 2) = [1.020_wp*cos(phi), 1.020_wp*sin(phi), 1.1606_wp]
            phi = phi + pi/3.0_wp
            xyz(:, k + 5) = [1.020_wp*cos(phi), 1.020_wp*sin(phi), -1.1606_wp]
         end do
         call new(mol, [6, 6, 1, 1, 1, 1, 1, 1], xyz*aatoau)
      case ("mb16")
         call get_structure(mol, "MB16-43", "01")
      case ("c1_sphere")
         ! The c1 molecule; the fixtures of the tests give it a gap
         ! threshold that no pair of eigenvalues passes, so that the frame
         ! comes from the cubic moment
         call structure_of("c1", mol)
      case ("ci_sphere")
         ! The same for the ci structure: the axis from the quartic moment
         call structure_of("ci", mol)
      case ("co2")
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.16_wp, 0.0_wp, 0.0_wp, -1.16_wp], [3, 3])
         call new_tilted(mol, [6, 8, 8], xyz*aatoau)
      case ("hcn")
         xyz = reshape([0.0_wp, 0.0_wp, -1.064_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.156_wp], [3, 3])
         call new_tilted(mol, [1, 6, 7], xyz*aatoau)
      case ("co2_bent")
         ! Bent by 0.005 bohr: C2v, yet a line to the frame
         xyz = reshape([-2.2_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.005_wp, 0.0_wp, 2.2_wp, 0.0_wp, 0.0_wp], [3, 3])
         call new_tilted(mol, [8, 6, 8], xyz)
      case ("ring13")
         ! Regular 13-ring: no in-plane harmonic up to the largest order used
         allocate (xyz(3, 13))
         do k = 1, 13
            phi = 2.0_wp*pi*real(k, wp)/13.0_wp
            xyz(:, k) = 6.0_wp*[cos(phi), sin(phi), 0.0_wp]
         end do
         call new_tilted(mol, [(6, k=1, 13)], xyz)
      case ("h2o2")
         ! Twisted, C2
         xyz = reshape([0.0_wp, 0.73_wp, 0.0_wp, 0.0_wp, -0.73_wp, 0.0_wp, 0.6_wp, 0.95_wp, 0.5_wp, &
            & -0.6_wp, -0.95_wp, 0.5_wp], [3, 4])
         call new_tilted(mol, [8, 8, 1, 1], xyz*aatoau)
      case ("h2o2_trans")
         ! Planar trans, C2h
         xyz = reshape([0.0_wp, 0.73_wp, 0.0_wp, 0.0_wp, -0.73_wp, 0.0_wp, 0.9_wp, 0.95_wp, 0.0_wp, &
            & -0.9_wp, -0.95_wp, 0.0_wp], [3, 4])
         call new_tilted(mol, [8, 8, 1, 1], xyz*aatoau)
      case ("ci")
         ! Four pairs of atoms through the origin, Ci
         allocate (xyz(3, 8))
         xyz(:, 1) = [0.25_wp, -0.10_wp, 0.70_wp]
         xyz(:, 3) = [1.35_wp, 0.40_wp, 1.20_wp]
         xyz(:, 5) = [-0.90_wp, 1.20_wp, 1.40_wp]
         xyz(:, 7) = [0.30_wp, -1.10_wp, 1.10_wp]
         xyz(:, 2:8:2) = -xyz(:, 1:7:2)
         call new_tilted(mol, [6, 6, 9, 9, 17, 17, 1, 1], xyz*aatoau)
      case ("ethene_twisted")
         ! Planes of the two CH2 groups 40 degrees apart, D2
         phi = 20.0_wp*pi/180.0_wp
         allocate (xyz(3, 6))
         xyz(:, 1) = [0.0_wp, 0.0_wp, 0.67_wp]
         xyz(:, 2) = [0.0_wp, 0.0_wp, -0.67_wp]
         xyz(:, 3) = [0.92_wp*cos(phi), 0.92_wp*sin(phi), 1.24_wp]
         xyz(:, 4) = [-0.92_wp*cos(phi), -0.92_wp*sin(phi), 1.24_wp]
         xyz(:, 5) = [0.92_wp*cos(phi), -0.92_wp*sin(phi), -1.24_wp]
         xyz(:, 6) = [-0.92_wp*cos(phi), 0.92_wp*sin(phi), -1.24_wp]
         call new_tilted(mol, [6, 6, 1, 1, 1, 1], xyz*aatoau)
      case ("s4")
         ! Two orbits of the improper fourfold axis at different angles, S4
         allocate (xyz(3, 9), source=0.0_wp)
         xyz(:, 2) = [1.0_wp, 0.4_wp, 0.8_wp]
         xyz(:, 6) = [0.3_wp, 1.7_wp, 1.2_wp]
         do k = 1, 3
            xyz(:, k + 2) = [-xyz(2, k + 1), xyz(1, k + 1), -xyz(3, k + 1)]
            xyz(:, k + 6) = [-xyz(2, k + 5), xyz(1, k + 5), -xyz(3, k + 5)]
         end do
         call new_tilted(mol, [6, 9, 9, 9, 9, 17, 17, 17, 17], xyz*aatoau)
      case ("allene")
         ! D2d
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.31_wp, 0.0_wp, 0.0_wp, -1.31_wp, &
            & 0.93_wp, 0.0_wp, 1.87_wp, -0.93_wp, 0.0_wp, 1.87_wp, 0.0_wp, 0.93_wp, -1.87_wp, &
            & 0.0_wp, -0.93_wp, -1.87_wp], [3, 7])
         call new_tilted(mol, [6, 6, 6, 1, 1, 1, 1], xyz*aatoau)
      case ("pyramid")
         ! Square pyramid, C4v
         xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.7_wp, 1.9_wp, 0.0_wp, -0.1_wp, &
            & -1.9_wp, 0.0_wp, -0.1_wp, 0.0_wp, 1.9_wp, -0.1_wp, 0.0_wp, -1.9_wp, -0.1_wp], [3, 6])
         call new_tilted(mol, [54, 8, 9, 9, 9, 9], xyz*aatoau)
      case ("c6v")
         ! Hexagon with one atom on its axis, C6v
         allocate (xyz(3, 7))
         do k = 1, 6
            phi = pi*real(k - 1, wp)/3.0_wp
            xyz(:, k) = [1.40_wp*cos(phi), 1.40_wp*sin(phi), 0.0_wp]
         end do
         xyz(:, 7) = [0.0_wp, 0.0_wp, 1.8_wp]
         call new_tilted(mol, [6, 6, 6, 6, 6, 6, 3], xyz*aatoau)
      case ("pyrite")
         ! Twelve atoms at the cyclic images of (0, a, b) about a centre, Th
         allocate (xyz(3, 13), source=0.0_wp)
         call add_orbit(xyz(:, 2:13), [0.0_wp, 1.0_wp, 1.9_wp], .false.)
         call new_tilted(mol, [26, (16, k=1, 12)], xyz*aatoau)
      case ("icosahedron")
         ! Twelve atoms at the corners of an icosahedron, Ih
         allocate (xyz(3, 12))
         call add_orbit(xyz, [0.0_wp, 1.0_wp, 0.5_wp*(1.0_wp + sqrt(5.0_wp))], .false.)
         call new_tilted(mol, [(5, k=1, 12)], 0.9_wp*xyz*aatoau)
      case ("t12")
         ! Orbit of a generic point under the rotations of a tetrahedron, T
         allocate (xyz(3, 13), source=0.0_wp)
         call add_orbit(xyz(:, 2:13), [1.9_wp, 0.7_wp, 0.3_wp], .false.)
         call new_tilted(mol, [6, (9, k=1, 12)], xyz*aatoau)
      case ("o24")
         ! Orbit of a generic point under the rotations of a cube, O
         allocate (xyz(3, 25), source=0.0_wp)
         call add_orbit(xyz(:, 2:25), [2.4_wp, 1.1_wp, 0.5_wp], .true.)
         call new_tilted(mol, [16, (9, k=1, 24)], xyz*aatoau)
      case ("cube")
         ! Centred cube, Oh, its atoms on the threefold axes
         allocate (xyz(3, 9), source=0.0_wp)
         do k = 1, 8
            xyz(:, k + 1) = 0.95_wp*[merge(1.0_wp, -1.0_wp, btest(k - 1, 0)), &
               & merge(1.0_wp, -1.0_wp, btest(k - 1, 1)), merge(1.0_wp, -1.0_wp, btest(k - 1, 2))]
         end do
         call new_tilted(mol, [14, (6, k=1, 8)], xyz*aatoau)
      case default
         ! BF3
         allocate (xyz(3, 4))
         xyz(:, 1) = 0.0_wp
         do k = 1, 3
            phi = 2.0_wp*pi*real(k - 1, wp)/3.0_wp
            xyz(:, k + 1) = [1.30_wp*cos(phi), 1.30_wp*sin(phi), 0.0_wp]
         end do
         call new(mol, [5, 9, 9, 9], xyz*aatoau)
      end select
   end subroutine structure_of

   !> Structure at a generic pose, off the lab axes
   !>
   !> @param[out] mol Structure
   !> @param[in]  num Atomic numbers (nat)
   !> @param[in]  xyz Coordinates (3, nat), bohr
   subroutine new_tilted(mol, num, xyz)
      !> Structure
      type(structure_type), intent(out) :: mol
      !> Atomic numbers
      integer, intent(in) :: num(:)
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)

      call new(mol, num, matmul(rotation_matrix(0.5_wp, 1.2_wp, -0.4_wp), xyz) &
         & + spread([0.4_wp, -1.1_wp, 2.2_wp], 2, size(num)))
   end subroutine new_tilted

   !> Images of a point under the rotations of a tetrahedron or a cube
   !>
   !> - generators: the threefold axis (1, 1, 1), the twofold axis z, and with
   !>   `fourfold` the fourfold axis z
   !> - the orbit must have exactly `size(orbit, 2)` points
   !>
   !> @param[out] orbit    Images (3, n)
   !> @param[in]  point    Starting point
   !> @param[in]  fourfold Whether the fourfold axis joins (O, else T)
   pure subroutine add_orbit(orbit, point, fourfold)
      !> Images
      real(wp), intent(out) :: orbit(:, :)
      !> Starting point
      real(wp), intent(in) :: point(3)
      !> Whether the fourfold axis joins
      logical, intent(in) :: fourfold

      real(wp) :: image(3)
      integer :: n, k, ig, j
      logical :: known

      orbit = 0.0_wp
      orbit(:, 1) = point
      n = 1
      k = 0
      do while (k < n)
         k = k + 1
         do ig = 1, merge(3, 2, fourfold)
            select case (ig)
            case (1)
               image = [orbit(3, k), orbit(1, k), orbit(2, k)]
            case (2)
               image = [-orbit(1, k), -orbit(2, k), orbit(3, k)]
            case default
               image = [-orbit(2, k), orbit(1, k), orbit(3, k)]
            end select
            known = .false.
            do j = 1, n
               if (maxval(abs(orbit(:, j) - image)) < 1.0e-10_wp) known = .true.
            end do
            if (known .or. n >= size(orbit, 2)) cycle
            n = n + 1
            orbit(:, n) = image
         end do
      end do
   end subroutine add_orbit

   !* ================================================================================= *!
   !*                                   Rigid motions                                   *!
   !* ================================================================================= *!

   !> Rotation from z-y-z Euler angles
   !>
   !> @param[in] a First angle
   !> @param[in] b Second angle
   !> @param[in] c Third angle
   pure function rotation_matrix(a, b, c) result(q)
      !> First angle
      real(wp), intent(in) :: a
      !> Second angle
      real(wp), intent(in) :: b
      !> Third angle
      real(wp), intent(in) :: c
      !> Rotation
      real(wp) :: q(3, 3)

      real(wp) :: ra(3, 3), rb(3, 3), rc(3, 3)

      ra = reshape([cos(a), sin(a), 0.0_wp, -sin(a), cos(a), 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      rb = reshape([cos(b), 0.0_wp, -sin(b), 0.0_wp, 1.0_wp, 0.0_wp, sin(b), 0.0_wp, cos(b)], [3, 3])
      rc = reshape([cos(c), sin(c), 0.0_wp, -sin(c), cos(c), 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      q = matmul(ra, matmul(rb, rc))
   end function rotation_matrix

   !> Deterministic test rotation number k
   !>
   !> - 1: moderate Euler angles; 2: a half turn about (1, 1, 1); 3: close to
   !>   a half turn about every Euler axis; from 4 on, pseudo-random unit
   !>   quaternions covering every angle up to a half turn
   !>
   !> @param[in] k Rotation number, positive
   pure function test_rotation(k) result(q)
      !> Rotation number
      integer, intent(in) :: k
      !> Rotation
      real(wp) :: q(3, 3)

      real(wp) :: u(3), quat(4), w, x, y, z

      select case (k)
      case (1)
         q = rotation_matrix(0.81_wp, -1.27_wp, 2.33_wp)
      case (2)
         ! R = 2 u u^T - I, a rotation by pi about u
         u = [1.0_wp, 1.0_wp, 1.0_wp]/sqrt(3.0_wp)
         q = 2.0_wp*spread(u, 2, 3)*spread(u, 1, 3) - identity3
      case (3)
         q = rotation_matrix(3.05_wp, 3.02_wp, -3.09_wp)
      case default
         quat = sin(real(k, wp)*[1.31_wp, 2.17_wp, 0.73_wp, 1.93_wp] + [0.2_wp, 0.7_wp, 1.1_wp, 2.3_wp])
         quat = quat/norm2(quat)
         w = quat(1)
         x = quat(2)
         y = quat(3)
         z = quat(4)
         q = reshape([1.0_wp - 2.0_wp*(y*y + z*z), 2.0_wp*(x*y + w*z), 2.0_wp*(x*z - w*y), &
            & 2.0_wp*(x*y - w*z), 1.0_wp - 2.0_wp*(x*x + z*z), 2.0_wp*(y*z + w*x), &
            & 2.0_wp*(x*z + w*y), 2.0_wp*(y*z - w*x), 1.0_wp - 2.0_wp*(x*x + y*y)], [3, 3])
      end select
   end function test_rotation

   !> Deterministic test permutation number k of nat atoms
   !>
   !> - 1: reversed and cycled by one; otherwise a pseudo-random shuffle
   !>
   !> @param[in] nat Number of atoms
   !> @param[in] k   Permutation number, positive
   pure function test_permutation(nat, k) result(perm)
      !> Number of atoms
      integer, intent(in) :: nat
      !> Permutation number
      integer, intent(in) :: k
      !> Permutation: entry i is the original atom placed at i
      integer :: perm(nat)

      integer(int64) :: state
      integer :: i, j, tmp

      if (k == 1) then
         perm = cshift([(nat + 1 - i, i=1, nat)], 1)
         return
      end if
      perm = [(i, i=1, nat)]
      state = 12345_int64 + 7919_int64*int(k, int64)
      do i = nat, 2, -1
         state = modulo(1103515245_int64*state + 12345_int64, 2147483648_int64)
         j = 1 + int(modulo(state, int(i, int64)))
         tmp = perm(i)
         perm(i) = perm(j)
         perm(j) = tmp
      end do
   end function test_permutation

   !> Rotated, translated and permuted copy of a structure
   !>
   !> Atom i of the copy is atom perm(i) of the original at q x + shift
   !>
   !> @param[in]  mol   Structure
   !> @param[in]  q     Rotation
   !> @param[in]  shift Translation, bohr
   !> @param[in]  perm  Original atom of every atom of the copy (nat)
   !> @param[out] moved Moved copy
   subroutine rigid_motion(mol, q, shift, perm, moved)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Rotation
      real(wp), intent(in) :: q(3, 3)
      !> Translation
      real(wp), intent(in) :: shift(3)
      !> Original atom of every atom of the copy
      integer, intent(in) :: perm(:)
      !> Moved copy
      type(structure_type), intent(out) :: moved

      call new(moved, mol%num(mol%id(perm)), matmul(q, mol%xyz(:, perm)) + spread(shift, 2, mol%nat))
   end subroutine rigid_motion

   !> Largest deviation of a moved gradient from the rotated, permuted original
   !>
   !> grad_moved(:, i) should be q grad(:, perm(i)); relative to the largest
   !> component of the original. With an operation g of the molecule (the one
   !> `match_surface_points` matched, atom images `gperm`), the moved cavity
   !> is the image of g(original): grad_moved(:, i) = q g grad(:, a) with
   !> gperm(a) = perm(i)
   !>
   !> @param[in] grad       Gradient of the original (3, nat)
   !> @param[in] grad_moved Gradient of the moved copy (3, nat)
   !> @param[in] q          Rotation of the copy
   !> @param[in] perm       Original atom of every atom of the copy (nat)
   !> @param[in] g          Operation of the molecule
   !> @param[in] gperm      Image atom of every atom under `g` (nat)
   pure function gradient_covariance_error(grad, grad_moved, q, perm, g, gperm) result(err)
      !> Gradient of the original
      real(wp), intent(in) :: grad(:, :)
      !> Gradient of the moved copy
      real(wp), intent(in) :: grad_moved(:, :)
      !> Rotation
      real(wp), intent(in) :: q(3, 3)
      !> Original atom of every atom of the copy
      integer, intent(in) :: perm(:)
      !> Operation of the molecule
      real(wp), intent(in), optional :: g(3, 3)
      !> Image atoms of the operation
      integer, intent(in), optional :: gperm(:)
      !> Largest relative deviation
      real(wp) :: err

      real(wp) :: qg(3, 3)
      integer :: i, a

      qg = q
      if (present(g)) qg = matmul(q, g)
      err = 0.0_wp
      do i = 1, size(perm)
         a = perm(i)
         if (present(gperm)) a = findloc(gperm, perm(i), dim=1)
         err = max(err, maxval(abs(grad_moved(:, i) - matmul(qg, grad(:, a)))))
      end do
      err = err/max(maxval(abs(grad)), tiny(1.0_wp))
   end function gradient_covariance_error

   !> Largest deviation of moved per-sphere gradients from the original ones
   !>
   !> grad(:, j, i) = dX_i/dR_j of a per-sphere quantity X_i; as
   !> [[gradient_covariance_error]] for both atom indices: grad_moved(:, j, i)
   !> should be q g grad(:, b, a) with gperm(a) = perm(i), gperm(b) = perm(j)
   !>
   !> @param[in] grad       Per-sphere gradients of the original (3, nat, nat)
   !> @param[in] grad_moved Per-sphere gradients of the moved copy (3, nat, nat)
   !> @param[in] q          Rotation of the copy
   !> @param[in] perm       Original atom of every atom of the copy (nat)
   !> @param[in] g          Operation of the molecule
   !> @param[in] gperm      Image atom of every atom under `g` (nat)
   pure function sphere_gradient_covariance_error(grad, grad_moved, q, perm, g, gperm) result(err)
      !> Per-sphere gradients of the original
      real(wp), intent(in) :: grad(:, :, :)
      !> Per-sphere gradients of the moved copy
      real(wp), intent(in) :: grad_moved(:, :, :)
      !> Rotation
      real(wp), intent(in) :: q(3, 3)
      !> Original atom of every atom of the copy
      integer, intent(in) :: perm(:)
      !> Operation of the molecule
      real(wp), intent(in) :: g(3, 3)
      !> Image atoms of the operation
      integer, intent(in) :: gperm(:)
      !> Largest relative deviation
      real(wp) :: err

      real(wp) :: qg(3, 3)
      integer, allocatable :: back(:)
      integer :: i, j

      qg = matmul(q, g)
      back = [(findloc(gperm, perm(i), dim=1), i=1, size(perm))]
      err = 0.0_wp
      do i = 1, size(perm)
         do j = 1, size(perm)
            err = max(err, maxval(abs(grad_moved(:, j, i) - matmul(qg, grad(:, back(j), back(i))))))
         end do
      end do
      err = err/max(maxval(abs(grad)), tiny(1.0_wp))
   end function sphere_gradient_covariance_error

   !> Point-group operations of a molecule and the atom each one maps to
   !>
   !> libmsym's group of the exact reference; the identity alone when the
   !> detection fails or the group is linear
   !>
   !> @param[in]  mol    Structure
   !> @param[out] ops    Operations (3, 3, order) about `center`
   !> @param[out] perm   Image atom of every atom under every operation (nat, order)
   !> @param[out] center Fixed point
   subroutine molecule_symmetry(mol, ops, perm, center)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Operations
      real(wp), allocatable, intent(out) :: ops(:, :, :)
      !> Image atoms
      integer, allocatable, intent(out) :: perm(:, :)
      !> Fixed point
      real(wp), intent(out) :: center(3)

      type(msym_group_type) :: group
      type(mctc_error), allocatable :: err
      real(wp) :: image(3), dist, best
      integer :: ig, ia, ib, iat

      center = sum(mol%xyz, dim=2)/real(mol%nat, wp)
      ops = reshape(identity3, [3, 3, 1])
      call find_point_group(mol%xyz, mol%num(mol%id), group, err)
      if (.not. allocated(err) .and. .not. group%linear) then
         if (allocated(group%ops)) then
            ops = group%ops
            center = group%center
         end if
      end if
      allocate (perm(mol%nat, size(ops, 3)))
      do ig = 1, size(ops, 3)
         do ia = 1, mol%nat
            image = center + matmul(ops(:, :, ig), mol%xyz(:, ia) - center)
            best = huge(1.0_wp)
            iat = ia
            do ib = 1, mol%nat
               if (mol%num(mol%id(ib)) /= mol%num(mol%id(ia))) cycle
               dist = norm2(mol%xyz(:, ib) - image)
               if (dist < best) then
                  best = dist
                  iat = ib
               end if
            end do
            perm(ia, ig) = iat
         end do
      end do
   end subroutine molecule_symmetry

   !* ================================================================================= *!
   !*                                  Surface points                                   *!
   !* ================================================================================= *!

   !> Surface points of a built cavity
   !>
   !> - NumSA stores no surface: its quadrature nodes c + r R u of every
   !>   sphere (normal R u, buried ones included) with their area before
   !>   switching (solid-angle times radial weight), from the cavity's own
   !>   cached nodes
   !>
   !> @param[in]  cav Built cavity
   !> @param[out] pts Points, normals, owners and values
   subroutine surface_points_of(cav, pts)
      !> Built cavity
      class(cavity_type), intent(in) :: cav
      !> Surface points
      type(surface_points_type), intent(out) :: pts

      select type (cav)
      type is (cavity_type_numsa)
         call numsa_nodes_of(cav, pts)
         return
      end select
      pts%n = cav%ngrid
      pts%xyz = cav%xyz(:, 1:pts%n)
      allocate (pts%normal(3, pts%n), source=0.0_wp)
      if (allocated(cav%normal0)) pts%normal = cav%normal0(:, 1:pts%n)
      pts%owner = cav%owner(1:pts%n)
      allocate (pts%val(nval, pts%n), source=0.0_wp)
      if (allocated(cav%a)) pts%val(1, :) = cav%a(1:pts%n)
      if (allocated(cav%f)) pts%val(2, :) = cav%f(1:pts%n)
      if (allocated(cav%xi0)) pts%val(3, :) = cav%xi0(1:pts%n)
      if (allocated(cav%v)) pts%val(4, :) = cav%v(1:pts%n)
   end subroutine surface_points_of

   !> Quadrature nodes of an updated NumSA cavity
   !>
   !> Node ip of sphere a at c_a + r_a u, u the node of the sphere's (single)
   !> grid component
   !>
   !> @param[in]  cav Updated NumSA cavity
   !> @param[out] pts Nodes, unit normals u, owners and pre-switching areas
   subroutine numsa_nodes_of(cav, pts)
      !> Updated NumSA cavity
      type(cavity_type_numsa), intent(in) :: cav
      !> Nodes
      type(surface_points_type), intent(out) :: pts

      real(wp) :: u(3)
      integer :: nleb, iat, ip, k

      nleb = size(cav%ang_grid, 2)
      pts%n = cav%nsph*nleb
      allocate (pts%xyz(3, pts%n), pts%normal(3, pts%n), pts%owner(pts%n))
      allocate (pts%val(nval, pts%n), source=0.0_wp)
      k = 0
      do iat = 1, cav%nsph
         do ip = 1, nleb
            k = k + 1
            u = cav%grid%node(:, cav%grid%set_first(cav%grid%comp_set(cav%grid%sph_first(iat))) - 1 + ip)
            pts%xyz(:, k) = cav%sphxyz(:, iat) + cav%vdwsa(iat)*u
            pts%normal(:, k) = u
            pts%owner(k) = iat
            pts%val(1, k) = cav%ang_weight(ip)*cav%wrp(iat)
         end do
      end do
   end subroutine numsa_nodes_of

   !> Compare the points of a moved cavity with those of the original as sets
   !>
   !> The moved points map back by x -> q^T (x - shift), normals by q^T; the
   !> point of moved atom i belongs to original atom perm(i). With `ops`, the
   !> back-mapped set may also be the image of the original under one
   !> operation g of the molecule (about `center`, atom images `atom_perm`);
   !> the identity is tried first
   !>
   !> - every moved point is matched to the nearest unmatched original point
   !>   of its owner, within tol * (1 + largest coordinate) bohr
   !> - normals within tol, every value within tol * (largest of its kind)
   !>
   !> @param[in]  ref        Points of the original
   !> @param[in]  moved      Points of the moved cavity
   !> @param[in]  q          Rotation of the copy
   !> @param[in]  shift      Translation of the copy, bohr
   !> @param[in]  perm       Original atom of every moved atom (nat)
   !> @param[in]  tol        Relative tolerance
   !> @param[out] ok         Whether the sets match
   !> @param[out] message    Worst mismatch, empty on success
   !> @param[in]  ops        Symmetry operations of the molecule (3, 3, order)
   !> @param[in]  atom_perm  Image atom of every atom under every operation (nat, order)
   !> @param[in]  center     Fixed point of the operations
   !> @param[out] matched_op Operation the sets match under, 0 without `ops`
   !> @param[in]  allowed    Operations tried, from the identity at 0 (0:order); all if absent
   subroutine match_surface_points(ref, moved, q, shift, perm, tol, ok, message, ops, atom_perm, &
      & center, matched_op, allowed)
      !> Points of the original
      type(surface_points_type), intent(in) :: ref
      !> Points of the moved cavity
      type(surface_points_type), intent(in) :: moved
      !> Rotation
      real(wp), intent(in) :: q(3, 3)
      !> Translation
      real(wp), intent(in) :: shift(3)
      !> Original atom of every moved atom
      integer, intent(in) :: perm(:)
      !> Relative tolerance
      real(wp), intent(in) :: tol
      !> Whether the sets match
      logical, intent(out) :: ok
      !> Worst mismatch
      character(len=:), allocatable, intent(out) :: message
      !> Symmetry operations of the molecule
      real(wp), intent(in), optional :: ops(:, :, :)
      !> Image atoms of the operations
      integer, intent(in), optional :: atom_perm(:, :)
      !> Fixed point
      real(wp), intent(in), optional :: center(3)
      !> Operation of the match
      integer, intent(out), optional :: matched_op
      !> Operations tried
      logical, intent(in), optional :: allowed(0:)

      !> Moved points mapped back to the original frame
      type(surface_points_type) :: back
      !> Back-mapped points under the inverse of one operation
      type(surface_points_type) :: trial
      real(wp) :: g(3, 3), c(3)
      integer :: nop, iop, ip, ia
      character(len=:), allocatable :: first_message

      ok = .false.
      message = ""
      first_message = ""
      if (present(matched_op)) matched_op = 0
      if (ref%n /= moved%n) then
         message = "point counts differ: "//to_string(ref%n)//" vs "//to_string(moved%n)
         return
      end if

      back = moved
      do ip = 1, moved%n
         back%xyz(:, ip) = matmul(transpose(q), moved%xyz(:, ip) - shift)
         back%normal(:, ip) = matmul(transpose(q), moved%normal(:, ip))
         back%owner(ip) = perm(moved%owner(ip))
      end do

      nop = 1
      if (present(ops) .and. present(atom_perm) .and. present(center)) nop = size(ops, 3)
      do iop = 0, nop
         if (present(allowed)) then
            if (.not. allowed(iop)) cycle
         end if
         if (iop == 0) then
            trial = back
         else
            if (.not. present(ops)) exit
            ! back = c + g (y - c) for a point y of the original, owner(back) = perm_g(owner(y))
            g = ops(:, :, iop)
            c = center
            trial = back
            do ip = 1, back%n
               trial%xyz(:, ip) = c + matmul(transpose(g), back%xyz(:, ip) - c)
               trial%normal(:, ip) = matmul(transpose(g), back%normal(:, ip))
               ia = findloc(atom_perm(:, iop), back%owner(ip), dim=1)
               trial%owner(ip) = ia
            end do
         end if
         call match_once(ref, trial, tol, ok, message)
         if (len(first_message) == 0) first_message = message
         if (ok) then
            if (present(matched_op)) matched_op = iop
            message = ""
            return
         end if
      end do
      message = first_message
   end subroutine match_surface_points

   !> Match two point sets in one frame
   !>
   !> @param[in]  ref     Original points
   !> @param[in]  pts     Points in the original frame
   !> @param[in]  tol     Relative tolerance
   !> @param[out] ok      Whether they match
   !> @param[out] message Worst mismatch, empty on success
   subroutine match_once(ref, pts, tol, ok, message)
      !> Original points
      type(surface_points_type), intent(in) :: ref
      !> Points in the original frame
      type(surface_points_type), intent(in) :: pts
      !> Relative tolerance
      real(wp), intent(in) :: tol
      !> Whether they match
      logical, intent(out) :: ok
      !> Worst mismatch
      character(len=:), allocatable, intent(out) :: message

      logical, allocatable :: used(:)
      real(wp) :: scale, dist, best, vscale(nval)
      character(len=160) :: buf
      integer :: ip, jp, jbest, iv

      ok = .false.
      message = ""
      scale = tol*(1.0_wp + maxval(abs(ref%xyz)))
      do iv = 1, nval
         vscale(iv) = tol*max(maxval(abs(ref%val(iv, :))), tiny(1.0_wp))
      end do
      allocate (used(ref%n), source=.false.)
      do ip = 1, pts%n
         best = huge(1.0_wp)
         jbest = 0
         do jp = 1, ref%n
            if (used(jp) .or. ref%owner(jp) /= pts%owner(ip)) cycle
            dist = norm2(ref%xyz(:, jp) - pts%xyz(:, ip))
            if (dist < best) then
               best = dist
               jbest = jp
            end if
         end do
         if (jbest == 0) then
            message = "no unmatched point of sphere "//to_string(pts%owner(ip))
            return
         end if
         if (best > scale) then
            write (buf, '("nearest point ",es10.3," bohr away (sphere ",i0,")")') best, pts%owner(ip)
            message = trim(buf)
            return
         end if
         used(jbest) = .true.
         if (maxval(abs(ref%normal(:, jbest) - pts%normal(:, ip))) > tol) then
            write (buf, '("normal differs by ",es10.3)') maxval(abs(ref%normal(:, jbest) - pts%normal(:, ip)))
            message = trim(buf)
            return
         end if
         do iv = 1, nval
            if (abs(ref%val(iv, jbest) - pts%val(iv, ip)) > vscale(iv)) then
               write (buf, '(a," differs by ",es10.3," (",es10.3," vs ",es10.3,")")') trim(val_name(iv)), &
                  & abs(ref%val(iv, jbest) - pts%val(iv, ip)), ref%val(iv, jbest), pts%val(iv, ip)
               message = trim(buf)
               return
            end if
         end do
      end do
      ok = .true.
   end subroutine match_once

   !* ================================================================================= *!
   !*                                   Small helpers                                   *!
   !* ================================================================================= *!

   !> Rank-1 `same_bits`
   !>
   !> @param[in] a First array
   !> @param[in] b Second array
   pure function same_bits_1(a, b) result(same)
      !> First array
      real(wp), intent(in) :: a(:)
      !> Second array
      real(wp), intent(in) :: b(:)
      !> Whether shapes and bit patterns agree
      logical :: same

      same = all(shape(a) == shape(b))
      if (same) same = same_bits_flat(size(a), a, b)
   end function same_bits_1

   !> Rank-2 `same_bits`
   !>
   !> @param[in] a First array
   !> @param[in] b Second array
   pure function same_bits_2(a, b) result(same)
      !> First array
      real(wp), intent(in) :: a(:, :)
      !> Second array
      real(wp), intent(in) :: b(:, :)
      !> Whether shapes and bit patterns agree
      logical :: same

      same = all(shape(a) == shape(b))
      if (same) same = same_bits_flat(size(a), a, b)
   end function same_bits_2

   !> Rank-3 `same_bits`
   !>
   !> @param[in] a First array
   !> @param[in] b Second array
   pure function same_bits_3(a, b) result(same)
      !> First array
      real(wp), intent(in) :: a(:, :, :)
      !> Second array
      real(wp), intent(in) :: b(:, :, :)
      !> Whether shapes and bit patterns agree
      logical :: same

      same = all(shape(a) == shape(b))
      if (same) same = same_bits_flat(size(a), a, b)
   end function same_bits_3

   !> Rank-4 `same_bits`
   !>
   !> @param[in] a First array
   !> @param[in] b Second array
   pure function same_bits_4(a, b) result(same)
      !> First array
      real(wp), intent(in) :: a(:, :, :, :)
      !> Second array
      real(wp), intent(in) :: b(:, :, :, :)
      !> Whether shapes and bit patterns agree
      logical :: same

      same = all(shape(a) == shape(b))
      if (same) same = same_bits_flat(size(a), a, b)
   end function same_bits_4

   !> Whether two real sequences agree bit for bit
   !>
   !> - element by element: no array temporaries, which ifx would place on the
   !>   stack and overflow it for whole-cavity Jacobians
   !>
   !> @param[in] n Number of elements
   !> @param[in] a First sequence (n)
   !> @param[in] b Second sequence (n)
   pure function same_bits_flat(n, a, b) result(same)
      !> Number of elements
      integer, intent(in) :: n
      !> First sequence
      real(wp), intent(in) :: a(n)
      !> Second sequence
      real(wp), intent(in) :: b(n)
      !> Whether every bit pattern agrees
      logical :: same

      integer :: i

      same = .true.
      do i = 1, n
         if (transfer(a(i), 0_int64) /= transfer(b(i), 0_int64)) then
            same = .false.
            return
         end if
      end do
   end function same_bits_flat

end module test_cavity_symmetry_fixtures

! ============================================================================================
! PCM component fixtures
! ============================================================================================

!> Shared fixtures of the PCM component tests, NOT a test suite
!>
!> - Component factories of all nine models: `new_pcm` (difference-check physics,
!>   dense solves) and `make_component` (constructor defaults, optional context)
!> - Solutes of the difference checks and the iSwiG/SvdW-DROP seed cavities
!> - `fd_cavity`, a synthetic cavity whose surface fields the test sets directly
!> - Host of point charges answering every request of a staged phase
module test_model_component_pcm_fixtures
   use mctc_env, only: wp, fatal_error
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: error_type, check, test_failed
   use test_helpers, only: check_moist_error, component_view, point_charge_potential, &
      & point_charge_field, submit, get_test_cross
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_type, only: cavity_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, &
      & moist_cavity_iswig_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type, drop_solve_parameters_type, &
                                           drop_branch_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_radii, only: default_cpcm_radii
   use moist_channels_coupling, only: coupling_type, gaussian_potential_request_type, &
      & point_potential_request_type, nuclear_charge_request_type
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, &
      & moist_pcm_parameters_type, solver_type
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, &
      & new_component_cpcm, model_continuum_component_cosmo, new_component_cosmo
   use moist_model_continuum_component_pcm_iefpcm, only: model_continuum_component_iefpcm, &
      & new_component_iefpcm
   use moist_model_continuum_component_pcm_ssvpe, only: model_continuum_component_ssvpe, &
      & new_component_ssvpe
   use moist_model_continuum_component_pcm_dpcm, only: model_continuum_component_dpcm, &
      & new_component_dpcm
   use moist_model_continuum_component_pcm_dpcm_original, only: &
      & model_continuum_component_dpcm_original, new_component_dpcm_original
   use moist_model_continuum_component_pcm_lpb, only: model_continuum_component_lpb, &
      & new_component_lpb
   use moist_model_continuum_component_pcm_slic, only: model_continuum_component_slic, &
      & new_component_slic
   use moist_model_continuum_component_pcm_nonlocal, only: model_continuum_component_nonlocal, &
      & new_component_nonlocal
   implicit none(type, external)
   private

   public :: pcm_cpcm, pcm_cosmo, pcm_iefpcm, pcm_ssvpe, pcm_dpcm, pcm_dpcm_original, pcm_lpb, &
      & pcm_slic, pcm_nonlocal, npcm
   public :: pcm_physics, new_pcm, make_component
   public :: seed_iswig, seed_drop, seed_name, first_seed, new_seed_cavity, new_sibling_cavity
   public :: fd_cavity, get_seed
   public :: nmol, mol_triatomic, mol_label, get_solute
   public :: thr_vacuous
   public :: host_energy, stage_energy, answer_sources, host_charges

   !> Model identifiers of the component factory
   integer, parameter :: pcm_cpcm = 1, pcm_cosmo = 2, pcm_iefpcm = 3, pcm_ssvpe = 4, &
      & pcm_dpcm = 5, pcm_dpcm_original = 6, pcm_lpb = 7, pcm_slic = 8, pcm_nonlocal = 9
   !> Number of models
   integer, parameter :: npcm = 9

   !> Seeds of the synthetic cavity and cavities of the nuclear check
   integer, parameter :: seed_iswig = 1, seed_drop = 2
   !> Seed names for messages
   character(len=5), parameter :: seed_name(2) = ["iSwiG", "DROP "]

   !> Physical parameters of a component, defaults those of the difference checks
   type :: pcm_physics
      !> Static permittivity of every model, bulk permittivity of the nonlocal one
      real(wp) :: eps = 78.4_wp
      !> High-frequency permittivity of the nonlocal model
      real(wp) :: eps_inf = 1.8_wp
      !> Correlation length of the nonlocal model, bohr
      real(wp) :: lambda = 3.0_wp
      !> Relative permittivity inside the cavity of the nonlocal model
      real(wp) :: eps_solute = 1.0_wp
      !> Inverse Debye length of the LPB model, 1/bohr
      real(wp) :: kappa = 0.3_wp
      !> SLIC h(E_n) = alpha tanh(beta E_n - gamma) + mu, unsaturated at the test fields
      real(wp) :: slic_alpha = 0.3_wp
      !> SLIC beta, 1/(atomic field unit)
      real(wp) :: slic_beta = -2.0_wp
      !> SLIC gamma
      real(wp) :: slic_gamma = -0.4_wp
      !> SLIC mu
      real(wp) :: slic_mu = 0.05_wp
      !> SLIC static potential, hartree per elementary charge
      real(wp) :: slic_phi_static = 0.0_wp
   end type pcm_physics
   !> Newton tolerance of SLIC, far below the finite-difference tolerances
   real(wp), parameter :: newton_tol = 1.0e-13_wp

   !> Lebedev order of the seeding and nuclear-check cavities
   integer, parameter :: nleb = 26
   !> Switching cutoff of the seeds, keeping f well away from zero
   real(wp), parameter :: cut_f = 1.0e-2_wp
   !> Minimum point separation of the point seeds, bohr
   !>
   !> A point surface keeps seam points whatever their f; collocation of
   !> full-area points closer than their own extent is nearly singular, which
   !> spoils the differences (truncation or roundoff) on larger solutes
   real(wp), parameter :: min_sep = 1.0_wp
   !> Magnitude below which an expected-active quantity counts as vacuous
   real(wp), parameter :: thr_vacuous = 1.0e-7_wp

   !> Solutes of the difference checks
   integer, parameter :: nmol = 9
   !> Hand-placed bent triatomic without symmetry, the solute of the lifecycle check
   integer, parameter :: mol_triatomic = 1
   !> Five-carbon cross of `get_test_cross`
   integer, parameter :: mol_cross = 9
   !> Solute names for messages
   !>
   !> - triatomic: bent O-H-H, three spheres and no symmetry
   !> - IL16 008B: chloride, a single charged sphere (k1 = k2 on DROP too)
   !> - MB16-43 LiH: two spheres of very different radii
   !> - Amino20x4 GLY_xab: glycine, polar organic
   !> - But14diol 1: intramolecular hydrogen bond, concave crevices
   !> - MB16-43 01: mindless molecule with unusual elements and bonding
   !> - IL16 008A: dimethylimidazolium, net charge +1
   !> - Heavy28 bih3_h2o: bismuth (large radius) and a weakly bound water
   !> - C5 cross: the five-carbon cross of the DROP gradient tests, concave seams
   character(len=*), parameter :: mol_label(nmol) = [character(len=17) :: "triatomic", &
      & "IL16 008B", "MB16-43 LiH", "Amino20x4 GLY_xab", "But14diol 1", "MB16-43 01", &
      & "IL16 008A", "Heavy28 bih3_h2o", "C5 cross"]
   !> mstore collection of each solute (blank for the hand-placed ones)
   character(len=*), parameter :: mol_collection(nmol) = [character(len=9) :: "", "IL16", &
      & "MB16-43", "Amino20x4", "But14diol", "MB16-43", "IL16", "Heavy28", ""]
   !> mstore record of each solute (blank for the hand-placed ones)
   character(len=*), parameter :: mol_record(nmol) = [character(len=8) :: "", "008B", "LiH", &
      & "GLY_xab", "1", "01", "008A", "bih3_h2o", ""]
   !> Atomic numbers of the triatomic
   integer, parameter :: solute_num(3) = [8, 1, 1]
   !> Coordinates of the triatomic, bohr
   real(wp), parameter :: solute_xyz(3, 3) = reshape([0.0_wp, 0.0_wp, 0.1_wp, &
      & 1.45_wp, 1.05_wp, 0.0_wp, -1.40_wp, 1.15_wp, -0.2_wp], [3, 3])

   !> Synthetic cavity whose surface fields are set directly by the test
   type, extends(cavity_type) :: fd_cavity
   contains
      !> No-op: the fields are set by the test
      procedure :: update => fd_cavity_update
      !> No-op: the synthetic surface has no nuclear derivatives
      procedure :: get_gradient => fd_cavity_get_gradient
   end type fd_cavity

contains

   !> Component factory
   !>
   !> Dense solves; SLIC with a Newton tolerance far below the differences
   !>
   !> @param[in]  kind    Model identifier
   !> @param[in]  ctx     Run context borrowed by the component
   !> @param[out] pcm     Constructed component
   !> @param[out] error   Test failure
   !> @param[in]  physics Physical parameters, those of the difference checks by default
   subroutine new_pcm(kind, ctx, pcm, error, physics)
      !> Model identifier
      integer, intent(in) :: kind
      !> Run context borrowed by the component
      type(moist_context_type), target, intent(in) :: ctx
      !> Constructed component
      class(model_continuum_component_pcm), allocatable, intent(out) :: pcm
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Physical parameters
      type(pcm_physics), intent(in), optional :: physics

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Solver controls
      type(moist_pcm_parameters_type) :: param
      !> Physical parameters in effect
      type(pcm_physics) :: phys

      param%solver = solver_type%dense
      param%newton_tol = newton_tol
      if (present(physics)) phys = physics
      select case (kind)
      case (pcm_cpcm)
         allocate (model_continuum_component_cpcm :: pcm)
      case (pcm_cosmo)
         allocate (model_continuum_component_cosmo :: pcm)
      case (pcm_iefpcm)
         allocate (model_continuum_component_iefpcm :: pcm)
      case (pcm_ssvpe)
         allocate (model_continuum_component_ssvpe :: pcm)
      case (pcm_dpcm)
         allocate (model_continuum_component_dpcm :: pcm)
      case (pcm_dpcm_original)
         allocate (model_continuum_component_dpcm_original :: pcm)
      case (pcm_lpb)
         allocate (model_continuum_component_lpb :: pcm)
      case (pcm_slic)
         allocate (model_continuum_component_slic :: pcm)
      case (pcm_nonlocal)
         allocate (model_continuum_component_nonlocal :: pcm)
      case default
         call test_failed(error, "unknown PCM model")
         return
      end select
      associate (eps => phys%eps)
         select type (pcm)
         type is (model_continuum_component_cpcm)
            call new_component_cpcm(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_cosmo)
            call new_component_cosmo(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_iefpcm)
            call new_component_iefpcm(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_ssvpe)
            call new_component_ssvpe(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_dpcm)
            call new_component_dpcm(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_dpcm_original)
            call new_component_dpcm_original(pcm, eps, err, param, ctx=ctx)
         type is (model_continuum_component_lpb)
            call new_component_lpb(pcm, eps, phys%kappa, err, param, ctx=ctx)
         type is (model_continuum_component_slic)
            call new_component_slic(pcm, eps, phys%slic_alpha, phys%slic_beta, phys%slic_gamma, &
               & err, param, mu=phys%slic_mu, phi_static=phys%slic_phi_static, ctx=ctx)
         type is (model_continuum_component_nonlocal)
            call new_component_nonlocal(pcm, eps, phys%eps_inf, phys%lambda, err, param, &
               & epsilon_solute=phys%eps_solute, ctx=ctx)
         end select
      end associate
      call check_moist_error(error, err, "component setup")
   end subroutine new_pcm

   !> Construct each PCM family with optional context forwarding
   !>
   !> @param[in] kind Constructor index
   !> @param[out] component Constructed component
   !> @param[out] err Library diagnostic
   !> @param[in] ctx Optional borrowed context
   subroutine make_component(kind, component, err, ctx)
      !> Constructor index
      integer, intent(in) :: kind
      !> Constructed component
      class(model_continuum_component_pcm), allocatable, intent(out) :: component
      !> Library diagnostic
      type(moist_error_type), allocatable, intent(out) :: err
      !> Optional borrowed context
      type(moist_context_type), target, optional, intent(in) :: ctx

      select case (kind)
      case (1)
         allocate (model_continuum_component_cpcm :: component)
         select type (component)
         type is (model_continuum_component_cpcm)
            call new_component_cpcm(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (2)
         allocate (model_continuum_component_cosmo :: component)
         select type (component)
         type is (model_continuum_component_cosmo)
            call new_component_cosmo(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (3)
         allocate (model_continuum_component_iefpcm :: component)
         select type (component)
         type is (model_continuum_component_iefpcm)
            call new_component_iefpcm(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (4)
         allocate (model_continuum_component_ssvpe :: component)
         select type (component)
         type is (model_continuum_component_ssvpe)
            call new_component_ssvpe(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (5)
         allocate (model_continuum_component_dpcm :: component)
         select type (component)
         type is (model_continuum_component_dpcm)
            call new_component_dpcm(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (6)
         allocate (model_continuum_component_dpcm_original :: component)
         select type (component)
         type is (model_continuum_component_dpcm_original)
            call new_component_dpcm_original(component, 80.0_wp, error=err, ctx=ctx)
         end select
      case (7)
         allocate (model_continuum_component_lpb :: component)
         select type (component)
         type is (model_continuum_component_lpb)
            call new_component_lpb(component, 80.0_wp, 0.1_wp, error=err, ctx=ctx)
         end select
      case (8)
         allocate (model_continuum_component_slic :: component)
         select type (component)
         type is (model_continuum_component_slic)
            call new_component_slic(component, 80.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, error=err, ctx=ctx)
         end select
      case (9)
         allocate (model_continuum_component_nonlocal :: component)
         select type (component)
         type is (model_continuum_component_nonlocal)
            call new_component_nonlocal(component, 80.0_wp, 2.0_wp, 1.0_wp, error=err, ctx=ctx)
         end select
      case default
         call fatal_error(err, "Unknown PCM component fixture")
      end select
   end subroutine make_component

   !> Energy with the traces of host point charges
   !>
   !> @param[in,out] pcm    Component under test
   !> @param[in]     mol    Solute
   !> @param[in,out] cavity Cavity, current for `mol`
   !> @param[in]     src    Structure carrying the host charge positions
   !> @param[in]     qsrc   Host charges
   !> @param[out]    energy Polarization energy
   !> @param[out]    error  Test failure
   subroutine host_energy(pcm, mol, cavity, src, qsrc, energy, error)
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Cavity
      class(cavity_type), intent(inout) :: cavity
      !> Host charge positions
      type(structure_type), intent(in) :: src
      !> Host charges
      real(wp), intent(in) :: qsrc(:)
      !> Polarization energy
      real(wp), intent(out) :: energy
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Host coupling
      type(coupling_type), target :: coupling

      energy = 0.0_wp
      call stage_energy(pcm, mol, cavity, coupling, error)
      if (allocated(error)) return
      call answer_sources(cavity, coupling, src, qsrc, real(mol%num(mol%id), wp))
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err, "energy")
   end subroutine host_energy

   !> Update the component and stage a new coupling for the energy phase
   !>
   !> With `keep` the operators of the last update stay, for new host traces on
   !> an unchanged cavity; only the charge cache sees the change
   !>
   !> @param[in,out] pcm      Component under test
   !> @param[in]     mol      Solute
   !> @param[in,out] cavity   Cavity, current for `mol`
   !> @param[out]    coupling Coupling staged for the energy phase
   !> @param[out]    error    Test failure
   !> @param[in]     keep     Skip the update (default update)
   subroutine stage_energy(pcm, mol, cavity, coupling, error, keep)
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Cavity
      class(cavity_type), intent(inout) :: cavity
      !> Coupling staged for the energy phase
      type(coupling_type), intent(out), target :: coupling
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Skip the update
      logical, intent(in), optional :: keep

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Update skipped
      logical :: skip

      skip = .false.
      if (present(keep)) skip = keep
      if (.not. skip) then
         call pcm%update(mol, cavity, err)
         call check_moist_error(error, err, "update")
         if (allocated(error)) return
      end if
      call pcm%new_coupling(cavity, coupling, err)
      call check_moist_error(error, err, "coupling setup")
      if (allocated(error)) return
      call pcm%prepare_energy(cavity, coupling, err)
      call check_moist_error(error, err, "energy staging")
   end subroutine stage_energy

   !> Answer every request of the staged phase from host point charges
   !>
   !> Potential and all its derivatives, Gaussian-probed on a Gaussian request;
   !> `z` as the nuclear charges
   !>
   !> @param[in]     cavity   Cavity the coupling was staged on
   !> @param[in,out] coupling Staged coupling
   !> @param[in]     src      Structure carrying the host charge positions
   !> @param[in]     qsrc     Host charges
   !> @param[in]     z        Nuclear charges of the solute
   subroutine answer_sources(cavity, coupling, src, qsrc, z)
      !> Cavity the coupling was staged on
      class(cavity_type), intent(in) :: cavity
      !> Staged coupling
      type(coupling_type), intent(inout), target :: coupling
      !> Host charge positions
      type(structure_type), intent(in) :: src
      !> Host charges
      real(wp), intent(in) :: qsrc(:)
      !> Nuclear charges of the solute
      real(wp), intent(in) :: z(:)

      !> Raw potential derivatives
      real(wp), allocatable :: g(:, :), gxi(:), hrr(:, :, :), hrxi(:, :)

      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            call point_charge_field(cavity, qsrc, src, .true., g, gxi, hrr, hrxi)
            call submit(coupling, "phi", point_charge_potential(cavity, qsrc, src, .true.))
            call submit(coupling, "dphi_dr", g)
            call submit(coupling, "dphi_dxi", gxi)
            call submit(coupling, "d2phi_dr2", hrr)
            call submit(coupling, "d2phi_dr_dxi", hrxi)
         type is (point_potential_request_type)
            call point_charge_field(cavity, qsrc, src, .false., g, gxi, hrr, hrxi)
            call submit(coupling, "phi", point_charge_potential(cavity, qsrc, src, .false.))
            call submit(coupling, "dphi_dr", g)
            call submit(coupling, "d2phi_dr2", hrr)
         type is (nuclear_charge_request_type)
            call submit(coupling, "z", z)
         end select
      end do
   end subroutine answer_sources

   !> Host point charges on the solute atoms, summing to the solute charge
   !>
   !> @param[in] mol Solute
   pure function host_charges(mol) result(q)
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Host charges (nat)
      real(wp), allocatable :: q(:)
      !> Atom index
      integer :: iatom

      allocate (q(mol%nat))
      do iatom = 1, mol%nat
         q(iatom) = 0.6_wp*sin(1.7_wp*real(iatom, wp) + 0.3_wp)
      end do
      q = q + (mol%charge - sum(q))/real(mol%nat, wp)
   end function host_charges

   !> Solute of the difference checks
   !>
   !> @param[in]  imol Solute index
   !> @param[out] mol  Solute
   subroutine get_solute(imol, mol)
      !> Solute index
      integer, intent(in) :: imol
      !> Solute
      type(structure_type), intent(out) :: mol

      select case (imol)
      case (mol_triatomic)
         call new(mol, solute_num, solute_xyz)
      case (mol_cross)
         call get_test_cross(mol)
      case default
         call get_structure(mol, trim(mol_collection(imol)), trim(mol_record(imol)))
      end select
   end subroutine get_solute

   !> First cavity seed of a model's checks: SvdW-DROP for SLIC, else iSwiG
   !>
   !> Concave iSwiG seams can fold the SLIC interface condition (the component
   !> warns there), so SLIC is checked on the smooth cavity only
   !>
   !> @param[in] kind Model identifier
   pure function first_seed(kind) result(seed)
      !> Model identifier
      integer, intent(in) :: kind
      !> Cavity seed
      integer :: seed

      seed = merge(seed_drop, seed_iswig, kind == pcm_slic)
   end function first_seed

   !> Empty iSwiG or SvdW-DROP cavity with CPCM radii
   !>
   !> DROP projects to 1e-14 so that its surface is smooth at the difference steps
   !>
   !> @param[in]  seed   Cavity kind
   !> @param[in]  ctx    Run context borrowed by the cavity
   !> @param[out] cavity Cavity
   !> @param[out] error  Test failure
   subroutine new_seed_cavity(seed, ctx, cavity, error)
      !> Cavity kind
      integer, intent(in) :: seed
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(in) :: ctx
      !> Cavity
      class(cavity_type), allocatable, intent(out) :: cavity
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Concrete cavities
      type(cavity_type_iswig), allocatable :: iswig
      type(cavity_type_drop), allocatable :: drop
      !> Level set of the DROP cavity
      type(moist_cavity_drop_lsf_svdw_type) :: svdw

      select case (seed)
      case (seed_iswig)
         allocate (iswig)
         call new_cavity_iswig(iswig, radius_model=default_cpcm_radii(), error=err, &
            & param=moist_cavity_iswig_parameters_type(num_leb=nleb), ctx=ctx)
         call check_moist_error(error, err, "iSwiG setup")
         if (allocated(error)) return
         call move_alloc(iswig, cavity)
      case (seed_drop)
         allocate (drop)
         call svdw%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=1.5_wp, blend_3b=1.0_wp))
         call new_cavity_drop(drop, radius_model=default_cpcm_radii(), lsf_model=svdw, error=err, &
            & param=moist_cavity_drop_parameters_type(num_leb=nleb, tolerance=1.0e-14_wp, &
            & solve=drop_solve_parameters_type(max_iter=150), proj_level=1), ctx=ctx)
         call check_moist_error(error, err, "DROP setup")
         if (allocated(error)) return
         call move_alloc(drop, cavity)
      case default
         call test_failed(error, "unknown cavity seed")
      end select
   end subroutine new_seed_cavity

   !> DROP cavity whose anchors on `get_test_sibling_gate` hold partly open sibling points
   !>
   !> - SvdW `blend_k = 1.5`, `blend_3b = 3`, a 0.2 bohr branch window and the
   !>   default switches, 302 Lebedev points per atom
   !>
   !> @param[in]  ctx    Run context borrowed by the cavity
   !> @param[out] cavity Cavity, not updated
   !> @param[out] error  Test failure
   subroutine new_sibling_cavity(ctx, cavity, error)
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(in) :: ctx
      !> Cavity
      class(cavity_type), allocatable, intent(out) :: cavity
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Concrete cavity
      type(cavity_type_drop), allocatable :: drop
      !> Level set
      type(moist_cavity_drop_lsf_svdw_type) :: svdw

      allocate (drop)
      call svdw%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=1.5_wp, blend_3b=3.0_wp))
      call new_cavity_drop(drop, radius_model=default_cpcm_radii(), lsf_model=svdw, error=err, &
         & param=moist_cavity_drop_parameters_type(num_leb=302, tolerance=1.0e-14_wp, &
         & solve=drop_solve_parameters_type(max_iter=1000), proj_level=3, &
         & branch=drop_branch_parameters_type(dpdist_max=0.2_wp), do_fine=.true.), ctx=ctx)
      call check_moist_error(error, err, "DROP setup")
      if (allocated(error)) return
      call move_alloc(drop, cavity)
   end subroutine new_sibling_cavity

   !> Synthetic cavity seeded from the surface of the solute
   !>
   !> - iSwiG: k1 = k2, modulated around 1/R so that the curvature is not
   !>   piecewise constant
   !> - DROP: its principal curvatures, k1 /= k2
   !> - Only switched-on points (f > cut_f); widths and switching factors only
   !>   for a Gaussian surface
   !> - Point surface: no two points closer than `min_sep`, the most
   !>   switched-on point of a close group kept
   !>
   !> @param[in]  seed    Seed cavity
   !> @param[in]  smeared Keep the widths and switching factors
   !> @param[in]  mol     Solute
   !> @param[out] cavity  Synthetic cavity
   !> @param[out] error   Test failure
   subroutine get_seed(seed, smeared, mol, cavity, error)
      !> Seed cavity
      integer, intent(in) :: seed
      !> Keep the widths and switching factors
      logical, intent(in) :: smeared
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Synthetic cavity
      type(fd_cavity), intent(out) :: cavity
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the seeding cavity
      type(moist_context_type), target :: ctx
      !> Seeding cavity
      class(cavity_type), allocatable :: source
      !> Seed points kept
      integer, allocatable :: keep(:)
      !> Grid-point index
      integer :: i

      call new_context(ctx, verbosity=0)
      call new_seed_cavity(seed, ctx, source, error)
      if (allocated(error)) return
      call source%configure(normal=.true., curvature=.true., error=err)
      call check_moist_error(error, err, "seed configuration")
      if (allocated(error)) return
      call source%update(mol, err)
      call check_moist_error(error, err, "seed update")
      if (allocated(error)) return

      keep = pack([(i, i=1, source%ngrid)], source%f(1:source%ngrid) > cut_f)
      if (.not. smeared) call separate_points(source%xyz, source%f, keep)
      cavity%ngrid = size(keep)
      cavity%nsph = source%nsph
      cavity%xyz = source%xyz(:, keep)
      cavity%normal0 = source%normal0(:, keep)
      cavity%a = source%a(keep)
      cavity%k1 = source%k1(keep)
      cavity%k2 = source%k2(keep)
      if (seed == seed_iswig) then
         do i = 1, cavity%ngrid
            cavity%k1(i) = cavity%k1(i)*(1.0_wp + 0.2_wp*sin(0.37_wp*real(i, wp)))
         end do
         cavity%k2 = cavity%k1
      end if
      if (smeared) then
         cavity%xi0 = source%xi0(keep)
         cavity%f = source%f(keep)
      end if
      cavity%updated = .true.
      cavity%built = .true.
      call check(error, cavity%ngrid > 10, message="unexpected seed size")
      if (allocated(error)) return
      if (seed == seed_drop .and. mol%nat > 1) then
         call check(error, maxval(abs(cavity%k1 - cavity%k2)) > thr_vacuous, &
            & message="DROP seed: principal curvatures coincide")
      end if
   end subroutine get_seed

   !> Thin a point selection so that no two points lie closer than `min_sep`
   !>
   !> Greedy in order of decreasing switching factor; the kept indices stay in
   !> ascending order
   !>
   !> @param[in]     xyz  Point positions (3, n)
   !> @param[in]     f    Switching factors (n)
   !> @param[in,out] keep Selected point indices
   subroutine separate_points(xyz, f, keep)
      !> Point positions
      real(wp), intent(in) :: xyz(:, :)
      !> Switching factors
      real(wp), intent(in) :: f(:)
      !> Selected point indices
      integer, allocatable, intent(inout) :: keep(:)

      !> Candidates by decreasing f and the accepted mask
      integer, allocatable :: order(:)
      logical, allocatable :: accepted(:)
      !> Candidate indices
      integer :: i, j, ibest

      order = keep
      do i = 1, size(order) - 1
         ibest = i - 1 + maxloc(f(order(i:)), dim=1)
         j = order(i)
         order(i) = order(ibest)
         order(ibest) = j
      end do
      allocate (accepted(size(xyz, 2)), source=.false.)
      do i = 1, size(order)
         accepted(order(i)) = all(norm2(xyz(:, order(:i - 1)) &
            & - spread(xyz(:, order(i)), 2, i - 1), dim=1) >= min_sep .or. &
            & .not. accepted(order(:i - 1)))
      end do
      keep = pack(keep, accepted(keep))
   end subroutine separate_points

   !> No-op update: the fields are set by the test
   !>
   !> @param[in,out] self  Synthetic cavity
   !> @param[in]     mol   Solute, unused
   !> @param[out]    error Never allocated
   subroutine fd_cavity_update(self, mol, error)
      !> Synthetic cavity
      class(fd_cavity), intent(inout) :: self
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      associate (unused_self => self, unused_mol => mol)
      end associate
   end subroutine fd_cavity_update

   !> No-op: the synthetic surface has no nuclear derivatives
   !>
   !> @param[in,out] self  Synthetic cavity
   !> @param[out]    error Never allocated
   subroutine fd_cavity_get_gradient(self, error)
      !> Synthetic cavity
      class(fd_cavity), intent(inout) :: self
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      associate (unused_self => self)
      end associate
   end subroutine fd_cavity_get_gradient

end module test_model_component_pcm_fixtures
