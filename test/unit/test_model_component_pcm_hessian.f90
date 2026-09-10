!> Second-order unit tests for the PCM solvation-model component
!>
!> The component's `get_hessian_surface_weights` must return the directional
!> derivative of its gradient-path surface adjoints, and `get_direct_hessian`
!> that of its direct gradient, along a batch of nuclear directions with a
!> prescribed surface tangent. The reference differences the shipped first-order
!> hooks on a component rebuilt at every stencil point: the surface displaced
!> along the tangent, the nuclei along the direction, and the charges re-solved
!> on the displaced matrix, so every moving part of the response -- the moved
!> matrix, the moved potential, the charge response through the linear solve,
!> the second derivatives of the kernel and of the nuclear field -- is
!> differentiated together.
!>
!> The width, switching, position and charge responses are all live in this
!> fixture; the area, weight, normal and curvature tangents are set nonzero on
!> purpose because the PCM adjoints do not read them, so their response must
!> be exactly zero.
module test_model_component_pcm_hessian
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type
   use moist_cavity_iswig, only: cavity_type_iswig
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_surface_tangent, only: cavity_surface_tangent_type
   use moist_channels, only: coupling_type
   use moist_context, only: moist_context_type, new_context
   use moist_model_component_pcm_type, only: solver_type, potential_source
   use moist_model_component_pcm_cpcm, only: solvation_model_component_cpcm, new_component_cpcm
   use test_helpers, only: get_test_structures, center_at_origin, get_test_cavity_iswig
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_hessian

   !> Number of directions in the batch
   integer, parameter :: ndir = 2

   !> Dielectric constant of the component under test
   real(wp), parameter :: epsilon = 78.4_wp

   !> Lebedev grid size of the fixture
   integer, parameter :: nleb = 26

   !> Central-difference step and 4-point stencil
   real(wp), parameter :: step = 1.0e-3_wp
   integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
   real(wp), parameter :: COEFF(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp

   !> Agreement bounds of the differenced response
   real(wp), parameter :: fd_atol = 1.0e-9_wp
   real(wp), parameter :: fd_rtol = 1.0e-8_wp

   !> Below this a reference carries no information
   real(wp), parameter :: vacuity_thr = 1.0e-5_wp

contains

   !> Collect the suite
   subroutine collect_model_component_pcm_hessian(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("cpcm_hessian_surface_weights", test_hessian_surface_weights), &
                  new_unittest("cpcm_direct_hessian_stash", test_direct_hessian_stash), &
                  new_unittest("cpcm_factor_cache_follows_update", test_factor_cache), &
                  new_unittest("cpcm_solvers_agree", test_solvers_agree), &
                  new_unittest("cpcm_external_source_refused", test_external_refused), &
                  new_unittest("cpcm_vacuum_short_circuit", test_vacuum_short_circuit) &
                  ]

   end subroutine collect_model_component_pcm_hessian

   !* --------------------------------- Local helpers --------------------------------- *!

   !> Sampled structure with its iSwiG surface and a neutral point-charge coupling
   !>
   !> @param[out] mol      Structure
   !> @param[out] cavity   iSwiG surface of the structure
   !> @param[out] coupling Point-charge coupling
   !> @param[out] error    Test failure
   subroutine build_fixture(mol, cavity, coupling, error)
      type(structure_type), intent(out) :: mol
      type(cavity_type_iswig), intent(out) :: cavity
      type(coupling_type), intent(out) :: coupling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: qat(:)
      integer :: iat

      call get_test_structures(mols, 5)
      call center_at_origin(mols(1))
      mol = mols(1)
      allocate (qat(mol%nat))
      do iat = 1, mol%nat
         qat(iat) = 0.2_wp*sin(1.3_wp*real(iat, wp))
      end do
      qat = qat - sum(qat)/real(mol%nat, wp)
      coupling%electrostatics%qat = reshape(qat, [mol%nat, 1])

      call get_test_cavity_iswig(mol, cavity, err, nleb=nleb)
      if (allocated(err)) then
         call test_failed(error, "cavity setup failed: "//err%message)
         return
      end if

   end subroutine build_fixture

   !> A dense, reproducible tangent batch on every surface channel and a
   !> matching batch of nuclear directions
   !>
   !> @param[in]  cavity  Surface the tangent is prescribed on
   !> @param[in]  nat     Number of atoms
   !> @param[out] tangent Surface tangent batch
   !> @param[out] dirs    Nuclear directions (3, nat, ndir)
   subroutine prescribe(cavity, nat, tangent, dirs)
      type(cavity_type_iswig), intent(in) :: cavity
      integer, intent(in) :: nat
      type(cavity_surface_tangent_type), intent(out) :: tangent
      real(wp), allocatable, intent(out) :: dirs(:, :, :)

      integer :: igrid, idir, iaxis, iatom

      call tangent%init(cavity%ngrid, ndir, .false.)
      do idir = 1, ndir
         do igrid = 1, cavity%ngrid
            ! Relative width and switching tangents keep both positive at
            ! every stencil point
            tangent%d_xi(igrid, idir) = 0.1_wp*cavity%xi0(igrid)*cos(0.5_wp*igrid + idir)
            tangent%d_f(igrid, idir) = 0.1_wp*cavity%f(igrid)*sin(0.9_wp*igrid - idir)
            tangent%d_a(igrid, idir) = 0.3_wp*sin(0.7_wp*igrid + 1.1_wp*idir)
            tangent%d_w(igrid, idir) = 0.2_wp*cos(0.4_wp*igrid + idir)
            do iaxis = 1, 3
               tangent%d_xyz(iaxis, igrid, idir) = 0.3_wp*cos(0.3_wp*igrid + 0.8_wp*iaxis + idir)
               tangent%d_n(iaxis, igrid, idir) = 0.25_wp*sin(0.6_wp*igrid + 1.3_wp*iaxis - idir)
            end do
         end do
      end do
      allocate (dirs(3, nat, ndir))
      do idir = 1, ndir
         do iatom = 1, nat
            do iaxis = 1, 3
               dirs(iaxis, iatom, idir) = 0.5_wp + 0.4_wp*sin(0.37_wp*iaxis + 0.61_wp*iatom &
                                                                + 1.1_wp*idir)
            end do
         end do
      end do

   end subroutine prescribe

   !> A CPCM component on the point-charge potential source
   subroutine new_cpcm(ctx, cpcm, error, phi_source, eps, solver)
      type(moist_context_type), intent(in), target :: ctx
      type(solvation_model_component_cpcm), intent(out) :: cpcm
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: phi_source
      real(wp), intent(in), optional :: eps
      integer, intent(in), optional :: solver

      type(moist_error_type), allocatable :: err
      real(wp) :: eps_loc
      integer :: solver_loc

      eps_loc = epsilon
      if (present(eps)) eps_loc = eps
      solver_loc = solver_type%cholesky
      if (present(solver)) solver_loc = solver
      call new_component_cpcm(cpcm, ctx, eps_loc, solver=solver_loc, &
                              phi_source=phi_source, error=err)
      if (allocated(err)) then
         call test_failed(error, "CPCM construction failed: "//err%message)
         return
      end if

   end subroutine new_cpcm

   !> The second-order response of a freshly built component
   !>
   !> @param[in]  ctx      Run context
   !> @param[in]  mol      Structure
   !> @param[in]  cavity   Surface
   !> @param[in]  coupling Coupling data
   !> @param[in]  tangent  Surface tangent batch
   !> @param[in]  dirs     Nuclear directions
   !> @param[out] dacc     Adjoint response per direction
   !> @param[out] hvp      Direct columns per direction
   !> @param[out] error    Test failure
   subroutine analytic_response(cpcm, mol, cavity, coupling, tangent, dirs, dacc, hvp, error)
      type(solvation_model_component_cpcm), intent(inout) :: cpcm
      type(structure_type), intent(in) :: mol
      type(cavity_type_iswig), intent(inout) :: cavity
      type(coupling_type), intent(in) :: coupling
      type(cavity_surface_tangent_type), intent(in) :: tangent
      real(wp), intent(in) :: dirs(:, :, :)
      type(cavity_surface_adjoint_type), allocatable, intent(out) :: dacc(:)
      real(wp), allocatable, intent(out) :: hvp(:, :, :)
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      integer :: idir

      call cpcm%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if
      allocate (dacc(ndir))
      do idir = 1, ndir
         call dacc(idir)%init(cavity%ngrid)
      end do
      allocate (hvp(3, mol%nat, ndir), source=0.0_wp)
      call cpcm%get_hessian_surface_weights(coupling, cavity, dirs, tangent, dacc, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM second-order surface weights failed: "//err%message)
         return
      end if
      call cpcm%get_direct_hessian(coupling, cavity, dirs, tangent, hvp, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM direct Hessian failed: "//err%message)
         return
      end if

   end subroutine analytic_response

   !* ------------------------------------- Tests ------------------------------------- *!

   !> The adjoint and direct responses against the differenced first-order hooks
   subroutine test_hessian_surface_weights(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, trial_mol
      type(cavity_type_iswig) :: cavity, trial
      type(coupling_type) :: coupling
      type(solvation_model_component_cpcm) :: cpcm, probe
      type(cavity_surface_tangent_type) :: tangent
      type(cavity_surface_adjoint_type), allocatable :: dacc(:)
      type(cavity_surface_adjoint_type) :: acc, fd
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :), grad(:, :), fd_grad(:, :)
      real(wp) :: s, dev, ref
      integer :: idir, ioff
      character(len=64) :: context

      call new_context(ctx)
      call build_fixture(mol, cavity, coupling, error)
      if (allocated(error)) return
      call prescribe(cavity, mol%nat, tangent, dirs)
      call new_cpcm(ctx, cpcm, error)
      if (allocated(error)) return
      call analytic_response(cpcm, mol, cavity, coupling, tangent, dirs, dacc, hvp, error)
      if (allocated(error)) return

      call new_cpcm(ctx, probe, error)
      if (allocated(error)) return
      allocate (grad(3, mol%nat), fd_grad(3, mol%nat))
      do idir = 1, ndir
         call fd%init(cavity%ngrid)
         fd_grad = 0.0_wp
         do ioff = 1, size(OFFSET)
            s = real(OFFSET(ioff), wp)*step
            trial_mol = mol
            trial_mol%xyz = mol%xyz + s*dirs(:, :, idir)
            trial = cavity
            trial%xi0 = cavity%xi0 + s*tangent%d_xi(:, idir)
            trial%f = cavity%f + s*tangent%d_f(:, idir)
            trial%xyz = cavity%xyz + s*tangent%d_xyz(:, :, idir)
            call probe%update(trial_mol, trial, err)
            if (allocated(err)) then
               call test_failed(error, "displaced CPCM update failed: "//err%message)
               return
            end if
            grad = 0.0_wp
            call probe%get_direct_gradient(coupling, trial, grad, err)
            if (allocated(err)) then
               call test_failed(error, "displaced direct gradient failed: "//err%message)
               return
            end if
            call acc%init(cavity%ngrid)
            call probe%get_gradient_surface_weights(coupling, trial, acc, err)
            if (allocated(err)) then
               call test_failed(error, "displaced surface weights failed: "//err%message)
               return
            end if
            fd%w_xi = fd%w_xi + COEFF(ioff)*acc%w_xi/step
            fd%w_f = fd%w_f + COEFF(ioff)*acc%w_f/step
            fd%w_xyz = fd%w_xyz + COEFF(ioff)*acc%w_xyz/step
            fd_grad = fd_grad + COEFF(ioff)*grad/step
         end do

         if (min(maxval(abs(fd%w_xi)), maxval(abs(fd%w_f)), maxval(abs(fd%w_xyz)), &
                 maxval(abs(fd_grad))) <= vacuity_thr) then
            call test_failed(error, "CPCM adjoint response reference is vacuous")
            return
         end if

         write (context, "(a,i0)") "width-weight response, direction ", idir
         dev = maxval(abs(dacc(idir)%w_xi - fd%w_xi))
         ref = maxval(abs(fd%w_xi))
         call check(error, dev, 0.0_wp, thr=fd_atol + fd_rtol*ref, more=trim(context))
         if (allocated(error)) return
         write (context, "(a,i0)") "switching-weight response, direction ", idir
         dev = maxval(abs(dacc(idir)%w_f - fd%w_f))
         ref = maxval(abs(fd%w_f))
         call check(error, dev, 0.0_wp, thr=fd_atol + fd_rtol*ref, more=trim(context))
         if (allocated(error)) return
         write (context, "(a,i0)") "position-weight response, direction ", idir
         dev = maxval(abs(dacc(idir)%w_xyz - fd%w_xyz))
         ref = maxval(abs(fd%w_xyz))
         call check(error, dev, 0.0_wp, thr=fd_atol + fd_rtol*ref, more=trim(context))
         if (allocated(error)) return
         write (context, "(a,i0)") "direct Hessian columns, direction ", idir
         dev = maxval(abs(hvp(:, :, idir) - fd_grad))
         ref = maxval(abs(fd_grad))
         call check(error, dev, 0.0_wp, thr=fd_atol + fd_rtol*ref, more=trim(context))
         if (allocated(error)) return

         ! Channels the PCM adjoints never read must stay exactly zero
         call check(error, max(maxval(abs(dacc(idir)%w_a)), maxval(abs(dacc(idir)%w_w)), &
                               maxval(abs(dacc(idir)%w_n)), maxval(abs(dacc(idir)%w_k1)), &
                               maxval(abs(dacc(idir)%w_k2))), 0.0_wp, thr=0.0_wp, &
                    more="CPCM adjoint response wrote to a channel the PCM does not read")
         if (allocated(error)) return
      end do

   end subroutine test_hessian_surface_weights

   !> The direct columns are handed out for the block they were formed on only
   subroutine test_direct_hessian_stash(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(solvation_model_component_cpcm) :: cpcm
      type(cavity_surface_tangent_type) :: tangent
      type(cavity_surface_adjoint_type), allocatable :: dacc(:)
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :), other(:, :, :)
      integer :: idir

      call new_context(ctx)
      call build_fixture(mol, cavity, coupling, error)
      if (allocated(error)) return
      call prescribe(cavity, mol%nat, tangent, dirs)
      call new_cpcm(ctx, cpcm, error)
      if (allocated(error)) return
      call cpcm%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if
      allocate (hvp(3, mol%nat, ndir), source=0.0_wp)

      ! Before any surface response was formed
      call cpcm%get_direct_hessian(coupling, cavity, dirs, tangent, hvp, err)
      call check(error, allocated(err), more="direct Hessian without a stashed block was accepted")
      if (allocated(error)) return
      deallocate (err)
      call check(error, maxval(abs(hvp)), 0.0_wp, thr=0.0_wp, more="hvp was touched on failure")
      if (allocated(error)) return

      allocate (dacc(ndir))
      do idir = 1, ndir
         call dacc(idir)%init(cavity%ngrid)
      end do
      call cpcm%get_hessian_surface_weights(coupling, cavity, dirs, tangent, dacc, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM second-order surface weights failed: "//err%message)
         return
      end if

      ! Another block of the same shape, and another shape
      other = 2.0_wp*dirs
      call cpcm%get_direct_hessian(coupling, cavity, other, tangent, hvp, err)
      call check(error, allocated(err), more="direct Hessian for other directions was accepted")
      if (allocated(error)) return
      deallocate (err)
      call cpcm%get_direct_hessian(coupling, cavity, dirs(:, :, 1:1), tangent, hvp(:, :, 1:1), err)
      call check(error, allocated(err), more="direct Hessian for another block shape was accepted")
      if (allocated(error)) return
      deallocate (err)
      call check(error, maxval(abs(hvp)), 0.0_wp, thr=0.0_wp, more="hvp was touched on failure")
      if (allocated(error)) return

      ! The block it was formed on
      call cpcm%get_direct_hessian(coupling, cavity, dirs, tangent, hvp, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM direct Hessian failed: "//err%message)
         return
      end if
      call check(error, maxval(abs(hvp)) > vacuity_thr, more="direct Hessian columns are vacuous")
      if (allocated(error)) return

      ! A new geometry drops the stash
      call cpcm%update(mol, cavity, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM update failed: "//err%message)
         return
      end if
      call cpcm%get_direct_hessian(coupling, cavity, dirs, tangent, hvp, err)
      call check(error, allocated(err), more="a stale direct block survived the update")
      if (allocated(error)) return

   end subroutine test_direct_hessian_stash

   !> The cached factorization follows the matrix through an update
   subroutine test_factor_cache(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol, mol2
      type(cavity_type_iswig) :: cavity, cavity2
      type(coupling_type) :: coupling
      type(solvation_model_component_cpcm) :: reused, fresh
      type(cavity_surface_tangent_type) :: tangent, tangent2
      type(cavity_surface_adjoint_type), allocatable :: dacc(:), dacc_reused(:), dacc_fresh(:)
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: dirs(:, :, :), dirs2(:, :, :), hvp(:, :, :)
      real(wp), allocatable :: hvp_reused(:, :, :), hvp_fresh(:, :, :)
      real(wp) :: scale
      integer :: idir

      call new_context(ctx)
      call build_fixture(mol, cavity, coupling, error)
      if (allocated(error)) return
      call prescribe(cavity, mol%nat, tangent, dirs)

      ! First geometry: the factor is built and cached here
      call new_cpcm(ctx, reused, error)
      if (allocated(error)) return
      call analytic_response(reused, mol, cavity, coupling, tangent, dirs, dacc, hvp, error)
      if (allocated(error)) return

      ! Second geometry, with its own surface, on the same component and on a fresh one
      mol2 = mol
      mol2%xyz(1, 1) = mol%xyz(1, 1) + 0.3_wp
      mol2%xyz(3, 2) = mol%xyz(3, 2) - 0.2_wp
      call get_test_cavity_iswig(mol2, cavity2, err, nleb=nleb)
      if (allocated(err)) then
         call test_failed(error, "second cavity setup failed: "//err%message)
         return
      end if
      call prescribe(cavity2, mol2%nat, tangent2, dirs2)
      call analytic_response(reused, mol2, cavity2, coupling, tangent2, dirs2, dacc_reused, &
                             hvp_reused, error)
      if (allocated(error)) return
      call new_cpcm(ctx, fresh, error)
      if (allocated(error)) return
      call analytic_response(fresh, mol2, cavity2, coupling, tangent2, dirs2, dacc_fresh, &
                             hvp_fresh, error)
      if (allocated(error)) return

      scale = max(maxval(abs(hvp_fresh)), maxval(abs(dacc_fresh(1)%w_xyz)))
      call check(error, scale > vacuity_thr, more="second-geometry response is vacuous")
      if (allocated(error)) return
      do idir = 1, ndir
         call check(error, max(maxval(abs(dacc_reused(idir)%w_xi - dacc_fresh(idir)%w_xi)), &
                               maxval(abs(dacc_reused(idir)%w_f - dacc_fresh(idir)%w_f)), &
                               maxval(abs(dacc_reused(idir)%w_xyz - dacc_fresh(idir)%w_xyz)), &
                               maxval(abs(hvp_reused(:, :, idir) - hvp_fresh(:, :, idir)))), &
                    0.0_wp, thr=0.0_wp, &
                    more="the reused component answers with the factor of the old matrix")
         if (allocated(error)) return
      end do

   end subroutine test_factor_cache

   !> Every solver reproduces the Cholesky response through its own factored path
   !>
   !> LU and the explicit inverse factorize the same matrix another way and
   !> must agree to round-off; the iterative solver keeps no factor and solves
   !> the response columns to its own tolerance, tightened here.
   subroutine test_solvers_agree(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(solvation_model_component_cpcm) :: cpcm
      type(cavity_surface_tangent_type) :: tangent
      type(cavity_surface_adjoint_type), allocatable :: dacc_ref(:), dacc(:)
      real(wp), allocatable :: dirs(:, :, :), hvp_ref(:, :, :), hvp(:, :, :)
      real(wp) :: scale, dev, bound
      integer :: idir, isolver
      integer, parameter :: solvers(3) = [solver_type%lu, solver_type%inversion, &
                                          solver_type%iterative]
      character(len=*), parameter :: names(3) = [character(len=9) :: "LU", "inversion", &
                                                 "iterative"]
      !> Agreement of the direct factorizations, relative to the response scale
      real(wp), parameter :: direct_rtol = 1.0e-11_wp
      !> Tolerance of the iterative solver and the agreement it buys
      real(wp), parameter :: cg_tol = 1.0e-14_wp, cg_rtol = 1.0e-8_wp

      call new_context(ctx)
      call build_fixture(mol, cavity, coupling, error)
      if (allocated(error)) return
      call prescribe(cavity, mol%nat, tangent, dirs)
      call new_cpcm(ctx, cpcm, error)
      if (allocated(error)) return
      call analytic_response(cpcm, mol, cavity, coupling, tangent, dirs, dacc_ref, hvp_ref, error)
      if (allocated(error)) return
      scale = max(maxval(abs(hvp_ref)), maxval(abs(dacc_ref(1)%w_xyz)), &
                  maxval(abs(dacc_ref(1)%w_xi)))
      call check(error, scale > vacuity_thr, more="Cholesky response is vacuous")
      if (allocated(error)) return

      do isolver = 1, size(solvers)
         call new_cpcm(ctx, cpcm, error, solver=solvers(isolver))
         if (allocated(error)) return
         cpcm%solver_tol = cg_tol
         call analytic_response(cpcm, mol, cavity, coupling, tangent, dirs, dacc, hvp, error)
         if (allocated(error)) return
         bound = direct_rtol*scale
         if (solvers(isolver) == solver_type%iterative) bound = cg_rtol*scale
         dev = maxval(abs(hvp - hvp_ref))
         do idir = 1, ndir
            dev = max(dev, maxval(abs(dacc(idir)%w_xi - dacc_ref(idir)%w_xi)), &
                      maxval(abs(dacc(idir)%w_f - dacc_ref(idir)%w_f)), &
                      maxval(abs(dacc(idir)%w_xyz - dacc_ref(idir)%w_xyz)))
         end do
         call check(error, dev, 0.0_wp, thr=bound, &
                    more="the "//trim(names(isolver))//" solver disagrees with Cholesky")
         if (allocated(error)) return
      end do

   end subroutine test_solvers_agree

   !> An external potential source is refused by name
   subroutine test_external_refused(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(solvation_model_component_cpcm) :: cpcm
      type(cavity_surface_tangent_type) :: tangent
      type(cavity_surface_adjoint_type), allocatable :: dacc(:)
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: dirs(:, :, :)
      real(wp) :: energy
      integer :: idir, igrid

      call new_context(ctx)
      call build_fixture(mol, cavity, coupling, error)
      if (allocated(error)) return
      call prescribe(cavity, mol%nat, tangent, dirs)
      call new_cpcm(ctx, cpcm, error, phi_source=potential_source%external)
      if (allocated(error)) return

      allocate (coupling%electrostatics%phi(cavity%ngrid))
      do igrid = 1, cavity%ngrid
         coupling%electrostatics%phi(igrid) = 0.1_wp*sin(0.3_wp*igrid)
      end do
      allocate (coupling%electrostatics%qefield(3, cavity%ngrid), source=0.0_wp)
      call cpcm%update(mol, cavity, err)
      energy = 0.0_wp
      if (.not. allocated(err)) call cpcm%get_energy(coupling, cavity, energy, err)
      if (allocated(err)) then
         call test_failed(error, "external-source CPCM setup failed: "//err%message)
         return
      end if

      allocate (dacc(ndir))
      do idir = 1, ndir
         call dacc(idir)%init(cavity%ngrid)
      end do
      call cpcm%get_hessian_surface_weights(coupling, cavity, dirs, tangent, dacc, err)
      call check(error, allocated(err), &
                 more="second-order surface weights with an external potential were accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "second-order") > 0 .and. &
                 index(err%message, "CPCM") > 0, &
                 more="refusal does not name the channel and the component: "//err%message)
      if (allocated(error)) return
      do idir = 1, ndir
         call check(error, max(maxval(abs(dacc(idir)%w_xi)), maxval(abs(dacc(idir)%w_f)), &
                               maxval(abs(dacc(idir)%w_xyz))), 0.0_wp, thr=0.0_wp, &
                    more="the response accumulator was touched on failure")
         if (allocated(error)) return
      end do

   end subroutine test_external_refused

   !> At eps = 1 the response is exactly zero and no factor is built
   subroutine test_vacuum_short_circuit(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_type_iswig) :: cavity
      type(coupling_type) :: coupling
      type(solvation_model_component_cpcm) :: cpcm
      type(cavity_surface_tangent_type) :: tangent
      type(cavity_surface_adjoint_type), allocatable :: dacc(:)
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :)
      integer :: idir

      call new_context(ctx)
      call build_fixture(mol, cavity, coupling, error)
      if (allocated(error)) return
      call prescribe(cavity, mol%nat, tangent, dirs)
      call new_cpcm(ctx, cpcm, error, eps=1.0_wp)
      if (allocated(error)) return
      call analytic_response(cpcm, mol, cavity, coupling, tangent, dirs, dacc, hvp, error)
      if (allocated(error)) return

      call check(error, maxval(abs(hvp)), 0.0_wp, thr=0.0_wp, &
                 more="direct Hessian columns are not zero at eps = 1")
      if (allocated(error)) return
      do idir = 1, ndir
         call check(error, max(maxval(abs(dacc(idir)%w_xi)), maxval(abs(dacc(idir)%w_f)), &
                               maxval(abs(dacc(idir)%w_xyz))), 0.0_wp, thr=0.0_wp, &
                    more="adjoint response is not zero at eps = 1")
         if (allocated(error)) return
      end do
      call check(error, .not. cpcm%factor_valid, more="a factorization was built at eps = 1")

   end subroutine test_vacuum_short_circuit

end module test_model_component_pcm_hessian
