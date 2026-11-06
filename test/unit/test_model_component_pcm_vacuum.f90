!> PCM models in a medium without dielectric contrast
!>
!> - Vacuum (eps = 1; LPB kappa = 0; nonlocal eps = eps_inf = eps_solute = 1) and,
!>   for the nonlocal model, a uniform medium (eps = eps_inf = eps_solute = 4)
!> - Every model: energy, charges, trace adjoints, direct nuclear term and all
!>   surface weights vanish at any discretization, dense and iterative, Gaussian
!>   and point surfaces; the full nuclear gradient on a real cavity too
!> - Exact zeros except LPB, whose right-hand side cancels through S_kappa S^-1 at
!>   roundoff (about 1e-16)
module test_model_component_pcm_vacuum
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use test_helpers, only: check_moist_error, component_view
   use test_model_component_pcm_fixtures, only: npcm, pcm_physics, new_pcm, seed_drop, &
      & fd_cavity, get_seed, new_seed_cavity, mol_triatomic, get_solute, host_charges, &
      & stage_energy, answer_sources
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, &
      & moist_pcm_parameters_type, solver_type
   use moist_model_continuum_component_pcm_nonlocal, only: model_continuum_component_nonlocal, &
      & new_component_nonlocal
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_vacuum

   !> Bound on every quantity that vanishes, LPB roundoff included
   real(wp), parameter :: thr = 1.0e-12_wp

   !> Vacuum of every model
   type(pcm_physics), parameter :: vacuum = pcm_physics(eps=1.0_wp, eps_inf=1.0_wp, kappa=0.0_wp)

contains

   !> Collect all exported unit tests
   !>
   !> @param[out] testsuite Unit tests of this suite
   subroutine collect_model_component_pcm_vacuum(testsuite)
      !> Unit tests of this suite
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("vacuum_response", test_vacuum_response), &
         & new_unittest("vacuum_nuclear_gradient", test_vacuum_nuclear_gradient), &
         & new_unittest("nonlocal_uniform_medium", test_nonlocal_uniform_medium) &
         & ]
   end subroutine collect_model_component_pcm_vacuum

   !> Every model, solver and surface kind responds to nothing in vacuum
   subroutine test_vacuum_response(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Run context borrowed by the components
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Solute and synthetic cavity
      type(structure_type) :: mol
      type(fd_cavity) :: cavity
      !> Solver kinds
      integer, parameter :: solvers(2) = [solver_type%dense, solver_type%iterative]
      !> Model, solver and surface indices
      integer :: kind, isolver, isurface
      !> Gaussian surface
      logical :: smeared
      !> Failure context
      character(len=:), allocatable :: label

      call new_context(ctx, verbosity=0)
      call get_solute(mol_triatomic, mol)
      do isurface = 1, 2
         smeared = isurface == 1
         call get_seed(seed_drop, smeared, mol, cavity, error)
         if (allocated(error)) return
         do kind = 1, npcm
            do isolver = 1, size(solvers)
               call new_pcm(kind, ctx, pcm, error, physics=vacuum)
               if (allocated(error)) return
               pcm%param%solver = solvers(isolver)
               label = trim(pcm%name)//merge(", Gaussian", ", point   ", smeared) &
                  & //merge(", dense    ", ", iterative", isolver == 1)
               call check_response(error, pcm, mol, cavity, label)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_vacuum_response

   !> The full nuclear gradient of every model vanishes in vacuum on a real cavity
   subroutine test_vacuum_nuclear_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the cavity and the components
      type(moist_context_type), target :: ctx
      !> Component under test
      class(model_continuum_component_pcm), allocatable :: pcm
      !> Solute
      type(structure_type) :: mol
      !> SvdW-DROP cavity of the solute
      class(cavity_type), allocatable :: cavity
      !> Host coupling and the host part of the gradient phase
      type(coupling_type), target :: coupling
      type(response_type) :: response
      !> Nuclear gradient
      real(wp), allocatable :: gradient(:, :)
      !> Energy
      real(wp) :: energy
      !> Model index
      integer :: kind

      call new_context(ctx, verbosity=0)
      call get_solute(mol_triatomic, mol)
      do kind = 1, npcm
         call new_pcm(kind, ctx, pcm, error, physics=vacuum)
         if (allocated(error)) return
         call new_seed_cavity(seed_drop, ctx, cavity, error)
         if (allocated(error)) return
         call pcm%configure_cavity(cavity, err)
         call check_moist_error(error, err, trim(pcm%name)//": cavity configuration")
         if (allocated(error)) return
         call cavity%update(mol, err)
         call check_moist_error(error, err, trim(pcm%name)//": cavity update")
         if (allocated(error)) return
         call stage_gradient(error, pcm, mol, cavity, coupling, energy)
         if (allocated(error)) return
         allocate (gradient(3, mol%nat), source=0.0_wp)
         call pcm%get_gradient(component_view(coupling), cavity, response, gradient, err)
         call check_moist_error(error, err, trim(pcm%name)//": nuclear gradient")
         if (allocated(error)) return
         call check_zero(error, trim(pcm%name)//": energy", [energy])
         if (allocated(error)) return
         call check_zero(error, trim(pcm%name)//": nuclear gradient", reshape(gradient, [size(gradient)]))
         if (allocated(error)) return
         deallocate (gradient, cavity)
      end do
   end subroutine test_vacuum_nuclear_gradient

   !> Nonlocal response vanishes in any uniform medium, not only in vacuum
   subroutine test_nonlocal_uniform_medium(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Run context borrowed by the component
      type(moist_context_type), target :: ctx
      !> Component under test
      type(model_continuum_component_nonlocal) :: pcm
      !> Solute and synthetic cavity
      type(structure_type) :: mol
      type(fd_cavity) :: cavity
      !> Surface index
      integer :: isurface
      !> Gaussian surface
      logical :: smeared

      call new_context(ctx, verbosity=0)
      call get_solute(mol_triatomic, mol)
      do isurface = 1, 2
         smeared = isurface == 1
         call get_seed(seed_drop, smeared, mol, cavity, error)
         if (allocated(error)) return
         call new_component_nonlocal(pcm, 4.0_wp, 4.0_wp, 3.0_wp, err, &
            & moist_pcm_parameters_type(solver=solver_type%dense), epsilon_solute=4.0_wp, ctx=ctx)
         call check_moist_error(error, err, "nonlocal setup")
         if (allocated(error)) return
         call check_response(error, pcm, mol, cavity, &
            & "nonlocal, eps = eps_inf = eps_solute = 4"//merge(", Gaussian", ", point   ", smeared))
         if (allocated(error)) return
      end do
   end subroutine test_nonlocal_uniform_medium

   !> Energy, charges, trace adjoints, direct term and surface weights all vanish
   !>
   !> @param[out]    error  Test failure
   !> @param[in,out] pcm    Component in a medium without contrast
   !> @param[in]     mol    Solute
   !> @param[in,out] cavity Synthetic cavity of the solute
   !> @param[in]     label  Failure context
   subroutine check_response(error, pcm, mol, cavity, label)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Component in a medium without contrast
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Synthetic cavity of the solute
      class(cavity_type), intent(inout) :: cavity
      !> Failure context
      character(len=*), intent(in) :: label
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Host coupling of the gradient phase
      type(coupling_type), target :: coupling
      !> Surface weights of the gradient
      type(cavity_surface_adjoint_type) :: acc
      !> Trace adjoints and the direct nuclear term
      real(wp), allocatable :: w_phi(:), w_vn(:), direct(:, :)
      !> Energy
      real(wp) :: energy

      call stage_gradient(error, pcm, mol, cavity, coupling, energy)
      if (allocated(error)) return
      call check_zero(error, label//": energy", [energy])
      if (allocated(error)) return
      call check_zero(error, label//": charges", pcm%q)
      if (allocated(error)) return
      call pcm%potential_adjoint(w_phi, err)
      call check_moist_error(error, err, label//": potential adjoint")
      if (allocated(error)) return
      call check_zero(error, label//": potential adjoint", w_phi)
      if (allocated(error)) return
      if (pcm%uses_normal_trace()) then
         call pcm%normal_trace_adjoint(w_vn, err)
         call check_moist_error(error, err, label//": normal-trace adjoint")
         if (allocated(error)) return
         call check_zero(error, label//": normal-trace adjoint", w_vn)
         if (allocated(error)) return
      end if
      allocate (direct(3, mol%nat), source=0.0_wp)
      call pcm%get_direct_gradient(component_view(coupling), cavity, direct, err)
      call check_moist_error(error, err, label//": direct nuclear term")
      if (allocated(error)) return
      call check_zero(error, label//": direct nuclear term", reshape(direct, [size(direct)]))
      if (allocated(error)) return
      call acc%init(cavity%ngrid)
      call pcm%get_gradient_surface_weights(component_view(coupling), cavity, acc, err)
      call check_moist_error(error, err, label//": surface weights")
      if (allocated(error)) return
      call check_zero(error, label//": w_xyz", reshape(acc%w_xyz, [size(acc%w_xyz)]))
      if (allocated(error)) return
      call check_zero(error, label//": w_n", reshape(acc%w_n, [size(acc%w_n)]))
      if (allocated(error)) return
      call check_zero(error, label//": w_k", [acc%w_k1, acc%w_k2])
      if (allocated(error)) return
      call check_zero(error, label//": w_xi, w_f, w_a, w_w", [acc%w_xi, acc%w_f, acc%w_a, acc%w_w])
   end subroutine check_response

   !> Solve the energy with host point charges on the nuclei, then stage the gradient
   !>
   !> @param[out]    error    Test failure
   !> @param[in,out] pcm      Component under test
   !> @param[in]     mol      Solute
   !> @param[in,out] cavity   Cavity, current for `mol`
   !> @param[out]    coupling Coupling answered for the gradient phase
   !> @param[out]    energy   Polarization energy
   subroutine stage_gradient(error, pcm, mol, cavity, coupling, energy)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Component under test
      class(model_continuum_component_pcm), intent(inout) :: pcm
      !> Solute
      type(structure_type), intent(in) :: mol
      !> Cavity
      class(cavity_type), intent(inout) :: cavity
      !> Coupling answered for the gradient phase
      type(coupling_type), intent(out), target :: coupling
      !> Polarization energy
      real(wp), intent(out) :: energy
      !> Library error handling
      type(moist_error_type), allocatable :: err

      energy = 0.0_wp
      call stage_energy(pcm, mol, cavity, coupling, error)
      if (allocated(error)) return
      call answer_sources(cavity, coupling, mol, host_charges(mol), real(mol%num(mol%id), wp))
      call pcm%get_energy(component_view(coupling), cavity, energy, err)
      call check_moist_error(error, err, trim(pcm%name)//": energy")
      if (allocated(error)) return
      call pcm%prepare_gradient(cavity, coupling, err)
      call check_moist_error(error, err, trim(pcm%name)//": gradient staging")
      if (allocated(error)) return
      call answer_sources(cavity, coupling, mol, host_charges(mol), real(mol%num(mol%id), wp))
   end subroutine stage_gradient

   !> Every value finite and within `thr` of zero
   !>
   !> @param[out] error  Test failure
   !> @param[in]  label  Failure context
   !> @param[in]  values Values that must vanish
   subroutine check_zero(error, label, values)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Failure context
      character(len=*), intent(in) :: label
      !> Values that must vanish
      real(wp), intent(in) :: values(:)
      !> Largest magnitude
      character(len=24) :: largest

      if (size(values) == 0) return
      write (largest, "(es12.4)") maxval(abs(values))
      call check(error, all(ieee_is_finite(values)), message=label//" is not finite")
      if (allocated(error)) return
      call check(error, maxval(abs(values)) <= thr, &
         & message=label//" does not vanish", more="largest magnitude "//trim(largest))
   end subroutine check_zero

end module test_model_component_pcm_vacuum
