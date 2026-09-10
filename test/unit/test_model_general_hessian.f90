!> The nuclear Hessian of the general solvation model
!>
!> `solvation_model_general%get_hessian` is, by definition, the derivative of
!> `get_gradient` at fixed host coupling data. The reference here is exactly
!> that: the shipped model gradient, central-differenced over rebuilt models,
!> so every moving part -- the cavity's second derivatives, the geometry
!> dependence of the weight fold, and the components' own adjoint response
!> through [[surface_adjoint_response_type]] -- is differentiated together.
!>
!> Two components on a DROP cavity carry a complete second-order surface
!> channel and are driven here: PV, whose adjoints move with every observable
!> they read, and CPCM on the point-charge potential source, whose adjoints
!> move with the surface, the nuclei and the re-solved surface charges. The
!> assertions, each for both components:
!>
!>   * `*_hessian_fd`: the analytic model Hessian against the differenced
!>     model gradient, at two steps, absolute *or* relative;
!>   * `*_frozen_misses_response`: the cavity's frozen-adjoint block, built
!>     from the same surface weights without the model's response, must miss
!>     that reference by a wide margin. Both blocks are symmetric and both
!>     pass every structural check, so only this separates a working response
!>     channel from an absent one;
!>   * `*_hvp_matches_dense`: the per-direction entry reproduces the dense
!>     block, to the bit on unit directions;
!>   * `pv_linearity`: two PV components at `p1` and `p2` give the Hessian of
!>     one at `p1 + p2`, so the response of a second component is added and
!>     not overwritten;
!>   * `hessian_guards`: an un-updated model and a component without a
!>     second-order channel (CPCM on an external potential, whose host data
!>     carry no response) are refused by name.
module test_model_general_hessian
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use testdrive, only: new_unittest, unittest_type, error_type, to_string, test_failed
   use moist_channels, only: coupling_type
   use moist_cavity_drop, only: cavity_type_drop, drop_hvp_per_dir_max
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_component_pcm_type, only: solver_type, potential_source
   use moist_model_component_pcm_cpcm, only: solvation_model_component_cpcm, new_component_cpcm
   use moist_model_components, only: solvation_model_component_pv, new_component_pv
   use moist_model_general, only: solvation_model_general, new_model_general
   use moist_context, only: moist_context_type
   use test_helpers, only: drop_fixture_geometry, build_drop_test_cavity, &
                           make_charge_coupling, LSF_SVDW, FIX_PLAIN

   implicit none(type, external)
   private

   public :: collect_model_general_hessian

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Component selectors of the shared fixture and reference routines
   integer, parameter :: MODEL_PV = 1, MODEL_CPCM = 2

   !> Pressure of the component under test; large enough that the block
   !> clears the vacuity threshold by orders of magnitude
   real(wp), parameter :: PRESSURE = 2.5_wp

   !> Dielectric constant and point charges of the CPCM component under test;
   !> the charges are large so that the block clears the vacuity threshold by
   !> an order of magnitude
   real(wp), parameter :: EPSILON = 32.0_wp
   real(wp), parameter :: QAT(3) = [1.2_wp, -1.8_wp, 0.6_wp]

   !> Central-difference steps of the reference; two, so that a value that
   !> agrees at one step only still fails
   real(wp), parameter :: FD_STEPS(*) = [3.0E-4_wp, 2.5E-4_wp]

   !> Agreement bound, absolute *and* relative: a component fails only when
   !> it misses both
   !>
   !> Measured, worst over the block and both steps: `2.1e-10` at `3e-4` and
   !> `3.3e-10` at `2.5e-4` on a block whose largest entry is `87`, the
   !> reference's round-off floor at this pressure (the volume Hessian of the
   !> cavity suite measures `9.4e-11` at unit pressure).
   real(wp), parameter :: HESS_TOL = 8.0E-10_wp

   !> The same bound for the CPCM block, measured separately: worst `2.2e-13`
   !> at `2.5e-4` on a block whose largest entry is `1.4e-2`, the reference's
   !> round-off floor for a gradient of that size (the frozen block misses by
   !> `2.2e-2`, ten decades over the bound)
   real(wp), parameter :: HESS_TOL_CPCM = 2.0E-12_wp

   !> Symmetry bound of the analytic block, relative to its largest entry
   real(wp), parameter :: SYM_TOL = 1.0E-11_wp

   !> Below this the reference carries no information
   real(wp), parameter :: VACUITY_THR = 1.0E-3_wp

   !> How far the frozen block must miss the reference, in units of `HESS_TOL`
   !>
   !> Measured: `28.3` absolute on a block whose largest entry is `92`, nine
   !> decades over the bound.
   real(wp), parameter :: TEETH_FACTOR = 1.0E2_wp

   !> Agreement of the per-direction entry with the dense block
   real(wp), parameter :: HVP_UNIT_TOL = 0.0_wp
   real(wp), parameter :: HVP_GEN_TOL = 1.0E-13_wp

   !> Agreement of the two-component sum with the single component
   real(wp), parameter :: LIN_TOL = 1.0E-13_wp

contains

   !> Collect the suite
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_model_general_hessian(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("pv_hessian_fd", test_pv_hessian_fd), &
                  new_unittest("pv_frozen_misses_response", test_pv_teeth), &
                  new_unittest("pv_hvp_matches_dense", test_pv_hvp), &
                  new_unittest("pv_linearity", test_pv_linearity), &
                  new_unittest("cpcm_hessian_fd", test_cpcm_hessian_fd), &
                  new_unittest("cpcm_frozen_misses_response", test_cpcm_teeth), &
                  new_unittest("cpcm_hvp_matches_dense", test_cpcm_hvp), &
                  new_unittest("hessian_guards", test_guards) &
                  ]
   end subroutine collect_model_general_hessian

   !* ================================================================================= *!
   !*                                    Fixture                                       *!
   !* ================================================================================= *!

   !> A general model with PV components on the DROP fixture
   !>
   !> @param[in]  mol       Structure
   !> @param[in]  pressures Pressure of each PV component
   !> @param[out] model     Updated model
   !> @param[out] ctx       Run context the model borrows; must outlive it
   !> @param[out] error     Error handle
   subroutine build_pv_model(mol, pressures, model, ctx, error)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Pressure of each PV component
      real(wp), intent(in) :: pressures(:)
      !> Updated model
      type(solvation_model_general), intent(out) :: model
      !> Run context
      type(moist_context_type), intent(inout), target :: ctx
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(solvation_model_component_pv) :: pv
      type(mctc_error), allocatable :: err
      integer :: ip

      call build_drop_test_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error, &
                                  want_fine=.true.)
      if (allocated(error)) return

      call new_model_general(model, cavity, ctx, err)
      if (.not. allocated(err)) then
         do ip = 1, size(pressures)
            call new_component_pv(pv, pressures(ip))
            call model%add_component(pv, err)
            if (allocated(err)) exit
         end do
      end if
      if (.not. allocated(err)) call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "model setup failed: "//err%message)
         return
      end if

   end subroutine build_pv_model

   !> A general model with one CPCM component on the point-charge potential
   !> source, on the DROP fixture
   !>
   !> @param[in]  mol   Structure
   !> @param[out] model Updated model
   !> @param[out] ctx   Run context the model borrows; must outlive it
   !> @param[out] error Error handle
   subroutine build_cpcm_model(mol, model, ctx, error)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Updated model
      type(solvation_model_general), intent(out) :: model
      !> Run context
      type(moist_context_type), intent(inout), target :: ctx
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(solvation_model_component_cpcm) :: cpcm
      type(mctc_error), allocatable :: err

      call build_drop_test_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error, &
                                  want_fine=.true.)
      if (allocated(error)) return

      call new_component_cpcm(cpcm, ctx, EPSILON, solver=solver_type%cholesky, &
                              phi_source=potential_source%charges, error=err)
      if (.not. allocated(err)) call new_model_general(model, cavity, ctx, err)
      if (.not. allocated(err)) call model%add_component(cpcm, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM model setup failed: "//err%message)
         return
      end if

   end subroutine build_cpcm_model

   !> The model of one selector, with the coupling it is driven by
   !>
   !> @param[in]  kind     Component selector
   !> @param[in]  mol      Structure
   !> @param[out] model    Updated model
   !> @param[out] ctx      Run context the model borrows; must outlive it
   !> @param[out] coupling Host coupling data of the component
   !> @param[out] error    Error handle
   subroutine build_model(kind, mol, model, ctx, coupling, error)
      !> Component selector
      integer, intent(in) :: kind
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Updated model
      type(solvation_model_general), intent(out) :: model
      !> Run context
      type(moist_context_type), intent(inout), target :: ctx
      !> Host coupling data
      type(coupling_type), intent(out) :: coupling
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      select case (kind)
      case (MODEL_PV)
         call build_pv_model(mol, [PRESSURE], model, ctx, error)
      case (MODEL_CPCM)
         if (mol%nat /= size(QAT)) then
            call test_failed(error, "the CPCM point charges do not match the fixture")
            return
         end if
         call make_charge_coupling(QAT, coupling)
         call build_cpcm_model(mol, model, ctx, error)
      case default
         call test_failed(error, "unknown component selector")
      end select

   end subroutine build_model

   !> The agreement bound of one selector
   pure real(wp) function hessian_tolerance(kind) result(tol)
      integer, intent(in) :: kind

      tol = HESS_TOL
      if (kind == MODEL_CPCM) tol = HESS_TOL_CPCM
   end function hessian_tolerance

   !> The persistent grid identity of the model's DROP cavity
   !>
   !> @param[in]  model     Updated model
   !> @param[out] numbering Persistent point ids
   !> @param[out] owner     Owner sphere per point
   !> @param[out] error     Error handle
   subroutine grid_identity(model, numbering, owner, error)
      !> Updated model
      type(solvation_model_general), intent(in) :: model
      !> Persistent point ids and owners
      integer, allocatable, intent(out) :: numbering(:), owner(:)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      select type (cav => model%cavity)
      type is (cavity_type_drop)
         numbering = cav%numbering(1:cav%ngrid)
         owner = cav%owner(1:cav%ngrid)
      class default
         call test_failed(error, "the model cavity is not a DROP cavity")
      end select

   end subroutine grid_identity

   !* ================================================================================= *!
   !*                          Finite-difference reference                             *!
   !* ================================================================================= *!

   !> The analytic model Hessian against the differenced model gradient
   subroutine test_pv_hessian_fd(error)
      type(error_type), allocatable, intent(out) :: error

      call run_fd(MODEL_PV, .true., error)
   end subroutine test_pv_hessian_fd

   !> The frozen-adjoint block must miss the differenced model gradient
   subroutine test_pv_teeth(error)
      type(error_type), allocatable, intent(out) :: error

      call run_fd(MODEL_PV, .false., error)
   end subroutine test_pv_teeth

   !> The analytic CPCM model Hessian against the differenced model gradient
   subroutine test_cpcm_hessian_fd(error)
      type(error_type), allocatable, intent(out) :: error

      call run_fd(MODEL_CPCM, .true., error)
   end subroutine test_cpcm_hessian_fd

   !> The frozen-adjoint CPCM block must miss the differenced model gradient
   subroutine test_cpcm_teeth(error)
      type(error_type), allocatable, intent(out) :: error

      call run_fd(MODEL_CPCM, .false., error)
   end subroutine test_cpcm_teeth

   !> Compare the model Hessian, or the cavity's frozen block built from the
   !> model's surface weights, to the differenced model gradient
   !>
   !> @param[in]  kind          Component selector
   !> @param[in]  with_response Whether the model Hessian (true) or the frozen
   !>                           cavity block (false) is compared
   !> @param[out] error         Error handle
   subroutine run_fd(kind, with_response, error)
      !> Component selector
      integer, intent(in) :: kind
      !> Whether the model Hessian or the frozen block is compared
      logical, intent(in) :: with_response
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(solvation_model_general) :: model
      type(moist_context_type), target :: ctx
      type(coupling_type) :: coupling
      type(mctc_error), allocatable :: err
      type(cavity_surface_adjoint_type) :: acc

      real(wp), allocatable :: hess(:, :, :, :), ref(:, :, :, :)
      real(wp) :: scale, sym, worst, worst_ref, tol
      integer :: nsph, istep, i1, i2, i3, i4

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_model(kind, mol, model, ctx, coupling, error)
      if (allocated(error)) return
      nsph = model%cavity%nsph
      tol = hessian_tolerance(kind)

      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      if (with_response) then
         call model%get_hessian(coupling, hess, err)
         if (allocated(err)) then
            call test_failed(error, "model Hessian failed: "//err%message)
            return
         end if
      else
         ! The same surface weights the model differentiates, contracted by
         ! the cavity with the adjoints held frozen
         call acc%init(model%cavity%ngrid)
         call model%components(1)%item%get_gradient_surface_weights(coupling, model%cavity, &
                                                                    acc, err)
         if (.not. allocated(err)) call model%cavity%get_hessian(acc, hess, err)
         if (allocated(err)) then
            call test_failed(error, "frozen cavity Hessian failed: "//err%message)
            return
         end if
      end if

      scale = maxval(abs(hess))
      if (scale <= VACUITY_THR) then
         call test_failed(error, "analytic Hessian is vacuous")
         return
      end if
      sym = 0.0_wp
      do i4 = 1, nsph
         do i3 = 1, ndim
            do i2 = 1, nsph
               do i1 = 1, ndim
                  sym = max(sym, abs(hess(i1, i2, i3, i4) - hess(i3, i4, i1, i2)))
               end do
            end do
         end do
      end do
      if (sym > SYM_TOL*scale) then
         call test_failed(error, "analytic Hessian is asymmetric: defect "//to_string(sym)// &
                          " against "//to_string(scale))
         return
      end if

      do istep = 1, size(FD_STEPS)
         call numerical_model_hessian(kind, mol, coupling, FD_STEPS(istep), ref, error)
         if (allocated(error)) return
         if (maxval(abs(ref)) <= VACUITY_THR) then
            call test_failed(error, "reference Hessian is vacuous")
            return
         end if

         call worst_deviation(hess, ref, tol, worst, worst_ref)
         if (with_response) then
            if (worst > tol .and. worst > tol*worst_ref) then
               call test_failed(error, "model Hessian mismatch (h = "// &
                                to_string(FD_STEPS(istep))//"): worst deviation "// &
                                to_string(worst)//" against reference "//to_string(worst_ref))
               return
            end if
         else
            if (worst <= TEETH_FACTOR*tol) then
               call test_failed(error, "the frozen-adjoint block reproduces the model"// &
                                " Hessian (h = "//to_string(FD_STEPS(istep))// &
                                "): worst deviation "//to_string(worst)// &
                                "; the model's adjoint response has no teeth")
               return
            end if
         end if
      end do

   end subroutine run_fd

   !> Five-point central difference of the model gradient over rebuilt models
   !>
   !> @param[in]  kind     Component selector
   !> @param[in]  mol      Base structure
   !> @param[in]  coupling Host coupling data, held fixed
   !> @param[in]  step     Central-difference step
   !> @param[out] hess     Differenced gradient `(3, nsph, 3, nsph)`
   !> @param[out] error    Error handle
   subroutine numerical_model_hessian(kind, mol, coupling, step, hess, error)
      !> Component selector
      integer, intent(in) :: kind
      !> Base structure
      type(structure_type), intent(in) :: mol
      !> Host coupling data
      type(coupling_type), intent(in) :: coupling
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Differenced gradient
      real(wp), allocatable, intent(out) :: hess(:, :, :, :)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
      real(wp), parameter :: COEFF(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp

      type(solvation_model_general) :: model
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol_disp
      type(mctc_error), allocatable :: err
      type(coupling_type) :: unused
      integer, allocatable :: ref_numbering(:), ref_owner(:), numbering(:), owner(:)
      real(wp), allocatable :: grad(:, :)
      integer :: nsph, batom, baxis, ioff

      call build_model(kind, mol, model, ctx, unused, error)
      if (allocated(error)) return
      call grid_identity(model, ref_numbering, ref_owner, error)
      if (allocated(error)) return
      nsph = model%cavity%nsph
      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      allocate (grad(ndim, nsph))

      do batom = 1, nsph
         do baxis = 1, ndim
            do ioff = 1, size(OFFSET)
               mol_disp = mol
               mol_disp%xyz(baxis, batom) = mol%xyz(baxis, batom) + real(OFFSET(ioff), wp)*step
               call model%update(mol_disp, err)
               if (allocated(err)) then
                  call test_failed(error, "model update on the displaced geometry failed: "// &
                                   err%message)
                  return
               end if

               ! A point that appears, vanishes, is reordered or changes owner
               ! puts a step into the differenced gradient; report it rather
               ! than absorb it
               call grid_identity(model, numbering, owner, error)
               if (allocated(error)) return
               if (size(numbering) /= size(ref_numbering)) then
                  call test_failed(error, "grid changed under the step "//to_string(step))
                  return
               end if
               if (any(numbering /= ref_numbering) .or. any(owner /= ref_owner)) then
                  call test_failed(error, "grid changed under the step "//to_string(step))
                  return
               end if

               grad = 0.0_wp
               call model%get_gradient(coupling, grad, err)
               if (allocated(err)) then
                  call test_failed(error, "model gradient failed: "//err%message)
                  return
               end if
               hess(:, :, baxis, batom) = hess(:, :, baxis, batom) + COEFF(ioff)*grad/step
            end do
         end do
      end do

   end subroutine numerical_model_hessian

   !> Worst deviation between two blocks among the entries that miss both bounds
   subroutine worst_deviation(hess, ref, tol, worst, worst_ref)
      real(wp), intent(in) :: hess(:, :, :, :), ref(:, :, :, :)
      real(wp), intent(in) :: tol
      real(wp), intent(out) :: worst, worst_ref

      real(wp) :: diff
      integer :: i1, i2, i3, i4

      worst = 0.0_wp
      worst_ref = 0.0_wp
      do i4 = 1, size(hess, 4)
         do i3 = 1, size(hess, 3)
            do i2 = 1, size(hess, 2)
               do i1 = 1, size(hess, 1)
                  diff = abs(hess(i1, i2, i3, i4) - ref(i1, i2, i3, i4))
                  if (diff > tol .and. diff > tol*abs(ref(i1, i2, i3, i4))) then
                     if (diff > worst) then
                        worst = diff
                        worst_ref = abs(ref(i1, i2, i3, i4))
                     end if
                  end if
               end do
            end do
         end do
      end do
   end subroutine worst_deviation

   !* ================================================================================= *!
   !*                        Products, linearity and the guards                        *!
   !* ================================================================================= *!

   !> Per-direction products reproduce the dense block
   subroutine test_pv_hvp(error)
      type(error_type), allocatable, intent(out) :: error

      call run_hvp(MODEL_PV, error)
   end subroutine test_pv_hvp

   !> Per-direction CPCM products reproduce the dense block
   subroutine test_cpcm_hvp(error)
      type(error_type), allocatable, intent(out) :: error

      call run_hvp(MODEL_CPCM, error)
   end subroutine test_cpcm_hvp

   !> Compare the per-direction entry of one selector to its dense block
   !>
   !> @param[in]  kind  Component selector
   !> @param[out] error Error handle
   subroutine run_hvp(kind, error)
      !> Component selector
      integer, intent(in) :: kind
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(solvation_model_general) :: model
      type(moist_context_type), target :: ctx
      type(coupling_type) :: coupling
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: hess(:, :, :, :), dirs(:, :, :), hvp(:, :, :), expect(:, :, :)
      integer :: nsph, ndir, idir, iatom, iaxis

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_model(kind, mol, model, ctx, coupling, error)
      if (allocated(error)) return
      nsph = model%cavity%nsph

      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      call model%get_hessian(coupling, hess, err)
      if (allocated(err)) then
         call test_failed(error, "model Hessian failed: "//err%message)
         return
      end if

      ! Unit directions: the dense columns, to the bit
      ndir = ndim*nsph
      allocate (dirs(ndim, nsph, ndir), source=0.0_wp)
      do iatom = 1, nsph
         do iaxis = 1, ndim
            dirs(iaxis, iatom, ndim*(iatom - 1) + iaxis) = 1.0_wp
         end do
      end do
      allocate (hvp(ndim, nsph, ndir), source=0.0_wp)
      call model%get_hvp(coupling, dirs, hvp, err)
      if (allocated(err)) then
         call test_failed(error, "model Hessian-vector products failed: "//err%message)
         return
      end if
      do iatom = 1, nsph
         do iaxis = 1, ndim
            idir = ndim*(iatom - 1) + iaxis
            if (maxval(abs(hvp(:, :, idir) - hess(:, :, iaxis, iatom))) > HVP_UNIT_TOL) then
               call test_failed(error, "unit-direction product differs from the dense column: "// &
                                to_string(maxval(abs(hvp(:, :, idir) - hess(:, :, iaxis, iatom)))))
               return
            end if
         end do
      end do
      deallocate (dirs, hvp)

      ! General directions, few enough for the per-direction fixed channel
      ndir = min(drop_hvp_per_dir_max, ndim*nsph - 1)
      allocate (dirs(ndim, nsph, ndir))
      do idir = 1, ndir
         do iatom = 1, nsph
            do iaxis = 1, ndim
               dirs(iaxis, iatom, idir) = 0.5_wp + 0.4_wp*sin(0.37_wp*iaxis + 0.61_wp*iatom &
                                                                + 1.1_wp*idir)
            end do
         end do
      end do
      allocate (hvp(ndim, nsph, ndir), source=0.0_wp)
      allocate (expect(ndim, nsph, ndir), source=0.0_wp)
      call model%get_hvp(coupling, dirs, hvp, err)
      if (allocated(err)) then
         call test_failed(error, "model Hessian-vector products failed: "//err%message)
         return
      end if
      do idir = 1, ndir
         do iatom = 1, nsph
            do iaxis = 1, ndim
               expect(:, :, idir) = expect(:, :, idir) &
                                    + hess(:, :, iaxis, iatom)*dirs(iaxis, iatom, idir)
            end do
         end do
      end do
      if (maxval(abs(hvp - expect)) > HVP_GEN_TOL*maxval(abs(hess))) then
         call test_failed(error, "general-direction product differs from the contracted"// &
                          " block: "//to_string(maxval(abs(hvp - expect))))
         return
      end if

   end subroutine run_hvp

   !> Two PV components sum to one at the summed pressure
   subroutine test_pv_linearity(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(solvation_model_general) :: model_one, model_two
      type(moist_context_type), target :: ctx_one, ctx_two
      type(coupling_type) :: coupling
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: hess_one(:, :, :, :), hess_two(:, :, :, :)
      real(wp), parameter :: P1 = 1.7_wp, P2 = 0.8_wp
      integer :: nsph

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_pv_model(mol, [P1 + P2], model_one, ctx_one, error)
      if (allocated(error)) return
      call build_pv_model(mol, [P1, P2], model_two, ctx_two, error)
      if (allocated(error)) return
      nsph = model_one%cavity%nsph

      allocate (hess_one(ndim, nsph, ndim, nsph), source=0.0_wp)
      allocate (hess_two(ndim, nsph, ndim, nsph), source=0.0_wp)
      call model_one%get_hessian(coupling, hess_one, err)
      if (.not. allocated(err)) call model_two%get_hessian(coupling, hess_two, err)
      if (allocated(err)) then
         call test_failed(error, "model Hessian failed: "//err%message)
         return
      end if
      if (maxval(abs(hess_one)) <= VACUITY_THR) then
         call test_failed(error, "model Hessian is vacuous")
         return
      end if
      if (maxval(abs(hess_one - hess_two)) > LIN_TOL*maxval(abs(hess_one))) then
         call test_failed(error, "two PV components do not sum to one: "// &
                          to_string(maxval(abs(hess_one - hess_two)))//" against "// &
                          to_string(maxval(abs(hess_one))))
         return
      end if

   end subroutine test_pv_linearity

   !> An un-updated model and a component without a second-order channel are refused
   subroutine test_guards(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(solvation_model_general) :: model
      type(cavity_type_drop), allocatable :: cavity
      type(solvation_model_component_cpcm) :: cpcm
      type(solvation_model_component_pv) :: pv
      type(moist_context_type), target :: ctx
      type(coupling_type) :: coupling
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: hess(:, :, :, :), hvp(:, :, :), dirs(:, :, :)
      integer :: nsph

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_drop_test_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error, &
                                  want_fine=.true.)
      if (allocated(error)) return
      nsph = cavity%nsph
      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      allocate (dirs(ndim, nsph, 1), source=1.0_wp)
      allocate (hvp(ndim, nsph, 1), source=0.0_wp)

      ! Before the first update
      call new_model_general(model, cavity, ctx, err)
      if (.not. allocated(err)) then
         call new_component_pv(pv, PRESSURE)
         call model%add_component(pv, err)
      end if
      if (allocated(err)) then
         call test_failed(error, "model setup failed: "//err%message)
         return
      end if
      call model%get_hessian(coupling, hess, err)
      if (.not. allocated(err)) then
         call test_failed(error, "the Hessian of an un-updated model was accepted")
         return
      end if
      deallocate (err)
      call model%get_hvp(coupling, dirs, hvp, err)
      if (.not. allocated(err)) then
         call test_failed(error, "Hessian-vector products of an un-updated model were accepted")
         return
      end if
      deallocate (err)
      if (maxval(abs(hess)) /= 0.0_wp .or. maxval(abs(hvp)) /= 0.0_wp) then
         call test_failed(error, "an accumulator was touched on failure")
         return
      end if

      ! A component without a second-order surface channel, named in the
      ! refusal: CPCM on an external potential, whose host data carry no
      ! response along a direction
      call new_component_cpcm(cpcm, ctx, 32.0_wp, solver=solver_type%cholesky, &
                              phi_source=potential_source%external, error=err)
      if (.not. allocated(err)) call new_model_general(model, cavity, ctx, err)
      if (.not. allocated(err)) call model%add_component(cpcm, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM model setup failed: "//err%message)
         return
      end if
      ! The surface charges exist once the energy has been evaluated, so the
      ! refusal below is the second-order hook's and not an earlier guard's
      block
         real(wp) :: energy
         integer :: igrid
         allocate (coupling%electrostatics%phi(model%cavity%ngrid))
         do igrid = 1, model%cavity%ngrid
            coupling%electrostatics%phi(igrid) = 0.1_wp*sin(0.3_wp*igrid)
         end do
         allocate (coupling%electrostatics%qefield(3, model%cavity%ngrid), source=0.0_wp)
         energy = 0.0_wp
         call model%get_energy(coupling, energy, err)
         if (allocated(err)) then
            call test_failed(error, "CPCM energy failed: "//err%message)
            return
         end if
      end block
      call model%get_hessian(coupling, hess, err)
      if (.not. allocated(err)) then
         call test_failed(error, "a CPCM Hessian was accepted although CPCM has no"// &
                          " second-order surface channel")
         return
      end if
      if (index(err%message, "second-order") == 0 .or. .not. allocated(cpcm%name)) then
         call test_failed(error, "refusal does not name the missing channel: "//err%message)
         return
      end if
      if (index(err%message, cpcm%name) == 0) then
         call test_failed(error, "refusal does not name the component: "//err%message)
         return
      end if
      if (maxval(abs(hess)) /= 0.0_wp) then
         call test_failed(error, "the Hessian accumulator was touched on failure")
         return
      end if

   end subroutine test_guards

end module test_model_general_hessian
