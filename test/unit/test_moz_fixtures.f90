!> Shared fixtures of the MOZ tests (not a suite)
!>
!> - Table water solvent with custom LJ and fixed charges
!> - Table water solvent with cSPC/E LJ and charges, for physical VV solves
!> - Synthetic UV grids, solvent response and smooth potential tables
!> - Reference 12-6 potential with Lorentz-Berthelot mixing
!> - Answers to pending atomic_charges requests
module test_moz_fixtures
   use mctc_io, only: structure_type, new
   use mctc_env, only: wp, moist_error => error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_channels_coupling, only: coupling_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   use moist_model_moz_3d_susceptibility, only: uv3d_susceptibility_type
   use test_helpers, only: get_qc_handymod_recipe
   use moist_solvent_system, only: solvent_system_type, new_solvent_system
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use moist_model_moz_potential_electrostatic_fixed, only: fixed_charges_type, new_fixed_charges
   use moist_model_moz_potential, only: lj_cspce, coulomb_cspce
   implicit none(type, external)
   private

   public :: make_uv3d_grid, make_uv3d_data
   public :: water_id, water_sigma, water_epsilon, water_charges
   public :: new_water_system, new_water_solvent, new_cspce_water_solvent, ref_lj, answer_charges

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

   !> Updated table water solvent with cSPC/E LJ and charges
   !>
   !> - Table geometry, not the SPC/E geometry
   !>
   !> @param[in]  ctx      run context, outlives the solvent
   !> @param[out] solvent  updated solvent
   !> @param[out] err      construction or update error
   subroutine new_cspce_water_solvent(ctx, solvent, err)
      !> Run context, outlives the solvent
      type(moist_context_type), intent(in), target :: ctx
      !> Updated solvent
      type(solvent_vv_type), intent(out) :: solvent
      !> Construction or update error
      type(moist_error), allocatable, intent(out) :: err
      type(solvent_system_type) :: system
      call new_water_system(system, err)
      if (.not. allocated(err)) call new_vv_solvent(solvent, ctx, system, err)
      if (.not. allocated(err)) call solvent%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call solvent%potential%add(coulomb_cspce, err)
      if (.not. allocated(err)) call solvent%update(err)
   end subroutine new_cspce_water_solvent

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

   !> Small UV grids with reproducible geometry and optional unequal weights
   !>
   !> @param[in]  kind  Cartesian (1), molecular (2), weighted molecular (3), weighted Cartesian (4)
   !> @param[out] grid  updated grid with two threads
   !> @param[out] err   construction failure
   subroutine make_uv3d_grid(kind, grid, err)
      !> Cartesian (1), molecular (2), weighted molecular (3), weighted Cartesian (4)
      integer, intent(in) :: kind
      !> Updated grid
      class(moist_math_grid_3d_type), allocatable, intent(out) :: grid
      !> Construction failure
      type(moist_error), allocatable, intent(out) :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(moist_math_grid_3d_molecular_type) :: mole
      type(moist_math_grid_atomic_recipe_type) :: recipe
      integer :: i

      call new(mol, [1, 8], reshape([-0.4_wp, 0.2_wp, 0.1_wp, 0.9_wp, -0.1_wp, 0.3_wp], [3, 2]))
      if (kind == 1 .or. kind == 4) then
         call new_cartesian_point_grid(cart, err, 9, 8, 7, 0.9_wp, margin=0.0_wp)
         if (.not. allocated(err)) call cart%update(mol, err)
         if (.not. allocated(err)) allocate (grid, source=cart)
      else
         call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=4.0_wp)
         if (.not. allocated(err)) call new_molecular_point_grid(mole, err, recipe=recipe, &
            & dr=1.1_wp, nufft_tol=1.0e-12_wp)
         if (.not. allocated(err)) call mole%update(mol, err)
         if (.not. allocated(err)) allocate (grid, source=mole)
      end if
      if (allocated(err)) return
      grid%nthreads = 2
      if (kind >= 3) then
         do i = 1, grid%ngrid
            grid%w(i) = grid%w(i)*(1.0_wp + 0.125_wp*real(mod(i, 5), wp))
         end do
      end if
   end subroutine make_uv3d_grid

   !> Nonsymmetric solvent response and smooth real and reciprocal UV potentials
   !>
   !> @param[in]  grid      updated grid
   !> @param[out] response  two-site response on the grid modes
   !> @param[out] sr        short-range potential (ngrid, 2)
   !> @param[out] lr        long-range potential (ngrid, 2)
   !> @param[out] uk        reciprocal long-range potential (npts_k, 2)
   subroutine make_uv3d_data(grid, response, sr, lr, uk)
      !> Updated grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Solvent response
      type(uv3d_susceptibility_type), intent(out) :: response
      !> Real-space potential tables
      real(wp), allocatable, intent(out) :: sr(:, :), lr(:, :)
      !> Reciprocal potential table
      complex(wp), allocatable, intent(out) :: uk(:, :)
      real(wp) :: r(3), k(3), k2
      integer :: i, col

      response%nv = 2
      allocate (response%values(2, 2, grid%npts_k), response%index(grid%npts_k), response%radius(grid%npts_k))
      allocate (response%charge_response(2, grid%npts_k))
      response%charge = [-0.7_wp, 0.7_wp]
      allocate (sr(grid%ngrid, 2), lr(grid%ngrid, 2), uk(grid%npts_k, 2))
      do i = 1, grid%npts_k
         k = grid%kpoint(i)
         k2 = sum(k**2)
         response%index(i) = i
         response%radius(i) = sqrt(k2)
         response%values(:, :, i) = reshape([0.04_wp, -0.02_wp, 0.03_wp, 0.06_wp], [2, 2])*exp(-k2)
         response%values(1, 1, i) = response%values(1, 1, i) + 1.0_wp
         response%values(2, 2, i) = response%values(2, 2, i) + 1.0_wp
         response%charge_response(:, i) = [0.13_wp, -0.09_wp]/(1.0_wp + k2)
         do col = 1, 2
            uk(i, col) = cmplx(0.03_wp*col, 0.02_wp*k(1), wp)*exp(-k2)
         end do
      end do
      do i = 1, grid%ngrid
         r = grid%point(i)
         do col = 1, 2
            sr(i, col) = 0.08_wp*cos(r(1) + col)*exp(-sum(r**2)/8.0_wp)
            lr(i, col) = 0.03_wp*sin(r(2) + col)/(1.0_wp + sum(r**2))
         end do
      end do
   end subroutine make_uv3d_data

end module test_moz_fixtures
