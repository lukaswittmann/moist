!> Reverse-mode nuclear gradient of the DROP cavity
!>
!> [[get_surface_gradient_drop]] contracts a surface adjoint straight into
!> `dE/dR_A` without ever forming the forward Jacobian. These tests pin it
!> against the forward `*_rA` arrays, which are independently validated by the
!> finite-difference suite in gradient.f90.
!>
!> The identity under test is, for every atom `A` and axis `beta`,
!>
!>   dE/dR_A,beta = sum_i [ w_xi_i  * xi1_rA(beta,A,i)
!>                        + w_f_i   * f1_rA(beta,A,i)
!>                        + w_a_i   * a_i1_rA(beta,A,i)
!>                        + w_w_i   * wleb1_rA(beta,A,i)
!>                        + w_xyz_i . xyz1_rA(:,beta,A,i)
!>                        + w_n_i   . normal1_rA(:,A,beta,i)
!>                        + w_k1_i  * k1_rA(beta,A,i)
!>                        + w_k2_i  * k2_rA(beta,A,i) ]
!>
!> `test_all_channels` drives every channel at once; `test_single_channels`
!> drives them one at a time so a bug in one cannot hide behind another.
module test_cavity_drop_nuclear_adjoint
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, test_failed
   use test_helpers, only: get_test_cross
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_radii, only: default_cpcm_radii
   use moist_data_radii_legacy, only: get_radius_func
   use moist_context, only: moist_context_type, new_context
   use moist_type, only: coupling_type
   use moist_model_general, only: general_solvation_model, new_general_model
   use moist_model_component_pcm_cpcm, only: cpcm, new_cpcm
   use moist_model_component_pcm_type, only: solver_type
   use moist_model_components, only: pv, new_pv
   use test_helpers, only: make_charge_coupling
   implicit none(type, external)
   private

   public :: collect_cavity_drop_nuclear_adjoint

   integer, parameter :: ndim = 3

   !> Level-set smoothing of the shared fixture
   real(wp), parameter :: k = 2.5_wp
   real(wp), parameter :: gamma = 1.0_wp
   integer, parameter :: NUM_LEB = 50

   real(wp), parameter :: PROJ_TOL = 1E-14_wp
   integer, parameter :: PROJ_MAXITER = 1000
   integer, parameter :: PROJ_LEVEL = 2

   !> Forward-versus-reverse agreement bounds. Both paths evaluate the same
   !> analytic derivatives, so the only difference is summation order and the
   !> reverse path's seed decomposition; the residual is pure round-off.
   real(wp), parameter :: EQ_ABS = 1.0E-9_wp
   real(wp), parameter :: EQ_REL = 1.0E-9_wp

   !> Below this the reference gradient is too small to carry a relative test
   real(wp), parameter :: VACUITY_THR = 1.0E-6_wp

   !> Softmax scale for the branching fixture. The production value set by
   !> `new_cavity_drop` (0.05) is peaked enough that the branch prune in
   !> filter.f90 keeps a single sibling and `branch_count` collapses to 1,
   !> which would silently retire the softmax reverse pass from the test.
   !> Measured branching onset for this fixture is s ~ 0.2.
   real(wp), parameter :: BRANCH_SOFTMAX_S = 1.0_wp

   !> Channel identifiers, in the order [[apply_channel]] understands
   integer, parameter :: CH_XI = 1, CH_F = 2, CH_A = 3, CH_W = 4
   integer, parameter :: CH_XYZ = 5, CH_N = 6, CH_K1 = 7, CH_K2 = 8
   integer, parameter :: NCHAN = 8

contains

   !> Register the reverse-mode nuclear gradient tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_cavity_drop_nuclear_adjoint(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("all_channels", test_all_channels), &
                  new_unittest("single_channels", test_single_channels), &
                  new_unittest("branching_channels", test_branching_channels), &
                  new_unittest("branching_xi_fd", test_branching_xi_fd), &
                  new_unittest("model_forward_reverse", test_model_forward_reverse), &
                  new_unittest("shape_guard", test_shape_guard) &
                  ]

   end subroutine collect_cavity_drop_nuclear_adjoint

   !> Every surface-adjoint channel driven simultaneously
   !>
   !> @param[out] error  Error handle
   subroutine test_all_channels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_equivalence([CH_XI, CH_F, CH_A, CH_W, CH_XYZ, CH_N, CH_K1, CH_K2], &
                           "all channels", .false., error)

   end subroutine test_all_channels

   !> Every channel again, on a geometry whose anchors branch
   !>
   !> The default fixture is asymmetric, so its multistart seeds all refine to
   !> one minimum and `branch_count` never exceeds one -- which leaves the
   !> softmax reverse pass and its explicit owner term completely untested.
   !> A near-symmetric dimer under multistart projection does branch, and
   !> [[run_equivalence]] asserts that it actually did.
   !>
   !> @param[out] error  Error handle
   subroutine test_branching_channels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Channel index
      integer :: ich
      !> Channel label
      character(len=16) :: label(NCHAN)

      label = [character(len=16) :: "branch w_xi", "branch w_f", "branch w_a", &
               "branch w_w", "branch w_xyz", "branch w_n", "branch w_k1", "branch w_k2"]

      do ich = 1, NCHAN
         call run_equivalence([ich], trim(label(ich)), .true., error)
         if (allocated(error)) return
      end do

   end subroutine test_branching_channels

   !> Each surface-adjoint channel driven on its own
   !>
   !> A channel that the reverse path drops entirely would still pass the
   !> combined test if another channel dominated the sum, so each one is also
   !> checked in isolation with its own non-vacuity guard.
   !>
   !> @param[out] error  Error handle
   subroutine test_single_channels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Channel index
      integer :: ich
      !> Channel label
      character(len=8) :: label(NCHAN)

      label = [character(len=8) :: "w_xi", "w_f", "w_a", "w_w", "w_xyz", "w_n", "w_k1", "w_k2"]

      do ich = 1, NCHAN
         call run_equivalence([ich], trim(label(ich)), .false., error)
         if (allocated(error)) return
      end do

   end subroutine test_single_channels

   !> Finite-difference check of the Gaussian-width channel under branching
   !>
   !> `test_branching_channels` only proves the reverse path matches the
   !> forward one. Both apply a branch-weight correction to `xi` that no
   !> shipped fixture ever exercised, so this pins the pair against an
   !> independent 4-point stencil on the scalar functional
   !>
   !>   L(R) = sum_i w_i * xi_i(R)
   !>
   !> Grid points are filtered and reordered on every rebuild, so the weights
   !> are keyed on the persistent `cavity%numbering` and restricted to points
   !> that survive at every stencil geometry -- otherwise a point appearing or
   !> vanishing would put a step discontinuity into L.
   !>
   !> @param[out] error  Error handle
   subroutine test_branching_xi_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Central-difference stencil
      integer, parameter :: NSTEP = 4
      real(wp), parameter :: FD_STEP = 1.0E-4_wp
      !> Displaced atom and axis. One coordinate is enough: the branch
      !> correction is the same code path for every atom and axis.
      integer, parameter :: FD_ATOM = 5, FD_AXIS = 3
      !> FD-versus-analytic bounds, limited by the stencil round-off floor
      real(wp), parameter :: FD_ABS = 5.0E-8_wp
      real(wp), parameter :: FD_REL = 5.0E-6_wp

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(moist_context_type), target :: ctx_fd
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol, mol_fd

      real(wp) :: fd_coeff(NSTEP), fd_delta(NSTEP)
      real(wp), allocatable :: grad_ana(:, :), wmap(:), xi_store(:, :), w_xi(:)
      logical, allocatable :: have(:, :), usable(:)
      real(wp) :: lvals(NSTEP), num_deriv, ana_deriv, diff
      integer :: max_num, istep, igrid, inum, ngrid, ncommon

      fd_coeff = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/(12.0_wp*FD_STEP)
      fd_delta = [-2.0_wp, -1.0_wp, 1.0_wp, 2.0_wp]*FD_STEP

      call fixture_geometry(.true., mol)
      call build_cavity(cavity, ctx, mol, .true., error)
      if (allocated(error)) return

      ngrid = cavity%ngrid
      if (.not. any(cavity%branch_count(1:ngrid) > 1)) then
         call test_failed(error, "FD fixture did not branch")
         return
      end if
      max_num = maxval(cavity%numbering(1:ngrid))

      allocate (wmap(max_num))
      do inum = 1, max_num
         wmap(inum) = ramp(inum, 0)
      end do
      allocate (xi_store(max_num, NSTEP), source=0.0_wp)
      allocate (have(max_num, NSTEP), source=.false.)

      ! Collect xi at each stencil geometry, keyed on the persistent numbering
      do istep = 1, NSTEP
         mol_fd = mol
         mol_fd%xyz(FD_AXIS, FD_ATOM) = mol_fd%xyz(FD_AXIS, FD_ATOM) + fd_delta(istep)
         block
            type(cavity_type_drop), allocatable :: cavity_fd
            call build_cavity(cavity_fd, ctx_fd, mol_fd, .true., error)
            if (allocated(error)) return
            do igrid = 1, cavity_fd%ngrid
               inum = cavity_fd%numbering(igrid)
               if (inum < 1 .or. inum > max_num) cycle
               have(inum, istep) = .true.
               xi_store(inum, istep) = cavity_fd%xi0(igrid)
            end do
         end block
      end do

      ! Only points present at every geometry carry weight
      allocate (usable(max_num))
      usable = all(have, dim=2)
      ncommon = count(usable)
      if (ncommon < 10) then
         call test_failed(error, "too few grid points survive the stencil ("// &
                          to_string(ncommon)//" of "//to_string(max_num)//")")
         return
      end if

      do istep = 1, NSTEP
         lvals(istep) = sum(wmap*xi_store(:, istep), mask=usable)
      end do
      num_deriv = sum(fd_coeff*lvals)

      ! Analytic side: same restricted weights through the reverse path
      allocate (w_xi(ngrid), source=0.0_wp)
      do igrid = 1, ngrid
         inum = cavity%numbering(igrid)
         if (inum >= 1 .and. inum <= max_num) then
            if (usable(inum)) w_xi(igrid) = wmap(inum)
         end if
      end do
      call acc%init(ngrid)
      call acc%add_surface_weights(cav_error, w_xi=w_xi)
      if (allocated(cav_error)) then
         call test_failed(error, "add_surface_weights failed: "//cav_error%message)
         return
      end if

      allocate (grad_ana(3, cavity%nsph), source=0.0_wp)
      call cavity%get_surface_gradient(acc, grad_ana, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "reverse gradient failed: "//cav_error%message)
         return
      end if
      ana_deriv = grad_ana(FD_AXIS, FD_ATOM)

      if (abs(num_deriv) <= VACUITY_THR) then
         call test_failed(error, "FD reference is vacuous ("//to_string(num_deriv)//")")
         return
      end if

      diff = abs(ana_deriv - num_deriv)
      if (diff > FD_ABS .and. diff > FD_REL*abs(num_deriv)) then
         call test_failed(error, "branching xi gradient disagrees with finite differences: "// &
                          "analytic "//to_string(ana_deriv)//" numeric "//to_string(num_deriv)// &
                          " (common points "//to_string(ncommon)//")")
         return
      end if

      ! Same stencil against the forward path. The branch post-pass has to
      ! reach xi1_rA as well as wleb1_rA and a_i1_rA -- it originally did not,
      ! and nothing caught it because no shipped fixture branches.
      call cavity%get_gradient(cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "forward gradient failed: "//cav_error%message)
         return
      end if
      ana_deriv = 0.0_wp
      do igrid = 1, ngrid
         ana_deriv = ana_deriv + w_xi(igrid)*cavity%xi1_rA(FD_AXIS, FD_ATOM, igrid)
      end do

      diff = abs(ana_deriv - num_deriv)
      if (diff > FD_ABS .and. diff > FD_REL*abs(num_deriv)) then
         call test_failed(error, "forward xi1_rA disagrees with finite differences under "// &
                          "branching: analytic "//to_string(ana_deriv)// &
                          " numeric "//to_string(num_deriv))
         return
      end if

   end subroutine test_branching_xi_fd

   !> Model-level gradient: reverse path must equal the legacy forward path
   !>
   !> This is the only test that drives `general_get_gradient` over a cavity
   !> that supports the surface contraction -- the model suite uses an iSwiG
   !> cavity, which falls back to the forward path -- so it is what actually
   !> exercises the component hooks `get_gradient_surface_weights` and
   !> `get_direct_gradient`.
   !>
   !> The two paths are required to agree, not merely to be close: the
   !> gradient-side surface weights were chosen to reproduce exactly the set
   !> of terms the forward path assembles.
   !>
   !> @param[out] error  Error handle
   subroutine test_model_forward_reverse(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Dielectric constant and pressure of the probe model
      real(wp), parameter :: epsilon_r = 32.0_wp
      real(wp), parameter :: pressure = 0.75_wp
      !> Both paths sum the same analytic terms, so only round-off separates them
      real(wp), parameter :: AB_ABS = 1.0E-9_wp
      real(wp), parameter :: AB_REL = 1.0E-9_wp
      !> Atomic charges driving the electrostatics
      real(wp), parameter :: qat_vals(*) = [0.20_wp, -0.15_wp, -0.05_wp]

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(general_solvation_model) :: model_rev, model_fwd
      type(cpcm) :: pcm_component
      type(pv) :: pv_component
      type(coupling_type) :: coupling
      type(structure_type) :: mol
      type(mctc_error), allocatable :: err

      real(wp), allocatable :: grad_rev(:, :), grad_fwd(:, :)
      real(wp) :: diff, scale
      integer :: nat, iatom, iaxis

      call fixture_geometry(.false., mol)
      call build_cavity(cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      nat = mol%nat

      call make_charge_coupling(qat_vals, coupling)
      call new_cpcm(pcm_component, ctx, epsilon_r, solver=solver_type%cholesky, error=err)
      if (allocated(err)) then
         call test_failed(error, "CPCM construction failed: "//err%message)
         return
      end if
      call new_pv(pv_component, pressure)

      call build_model(model_rev, cavity, ctx, pcm_component, pv_component, mol, error)
      if (allocated(error)) return
      call build_model(model_fwd, cavity, ctx, pcm_component, pv_component, mol, error)
      if (allocated(error)) return
      model_fwd%force_forward_gradient = .true.

      allocate (grad_rev(3, nat), source=0.0_wp)
      allocate (grad_fwd(3, nat), source=0.0_wp)

      call model_rev%get_gradient(coupling, grad_rev, err)
      if (allocated(err)) then
         call test_failed(error, "reverse-mode model gradient failed: "//err%message)
         return
      end if
      call model_fwd%get_gradient(coupling, grad_fwd, err)
      if (allocated(err)) then
         call test_failed(error, "forward model gradient failed: "//err%message)
         return
      end if

      if (maxval(abs(grad_fwd)) <= VACUITY_THR) then
         call test_failed(error, "model gradient is vacuous")
         return
      end if

      do iatom = 1, nat
         do iaxis = 1, 3
            diff = abs(grad_rev(iaxis, iatom) - grad_fwd(iaxis, iatom))
            scale = abs(grad_fwd(iaxis, iatom))
            if (diff > AB_ABS .and. diff > AB_REL*scale) then
               call test_failed(error, "model gradient differs between paths at atom "// &
                                to_string(iatom)//" axis "//to_string(iaxis)// &
                                ": reverse "//to_string(grad_rev(iaxis, iatom))// &
                                " forward "//to_string(grad_fwd(iaxis, iatom)))
               return
            end if
         end do
      end do

   end subroutine test_model_forward_reverse

   !> Assemble a CPCM + PV model on a DROP cavity
   !>
   !> @param[out]   model  Model to build
   !> @param[in]    cavity Cavity template copied into the model
   !> @param[in]    ctx    Run context owned by the caller
   !> @param[in]    pcmc   CPCM component template
   !> @param[in]    pvc    PV component template
   !> @param[in]    mol    Molecular structure
   !> @param[out]   error  Error handle
   subroutine build_model(model, cavity, ctx, pcmc, pvc, mol, error)
      !> Model to build
      type(general_solvation_model), intent(out) :: model
      !> Cavity template copied into the model
      type(cavity_type_drop), intent(in) :: cavity
      !> Run context owned by the caller
      type(moist_context_type), intent(in), target :: ctx
      !> Component templates
      type(cpcm), intent(in) :: pcmc
      type(pv), intent(in) :: pvc
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: err

      call new_general_model(model, cavity, ctx, err)
      if (.not. allocated(err)) call model%add_component(pcmc, err)
      if (.not. allocated(err)) call model%add_component(pvc, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "model setup failed: "//err%message)
         return
      end if

   end subroutine build_model

   !> A mis-shaped gradient accumulator must be rejected, not silently ignored
   !>
   !> @param[out] error  Error handle
   subroutine test_shape_guard(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      real(wp), allocatable :: bad_grad(:, :)

      call build_fixture(cavity, ctx, .false., error)
      if (allocated(error)) return

      call acc%init(cavity%ngrid)
      allocate (bad_grad(3, cavity%nsph + 1), source=0.0_wp)
      call cavity%get_surface_gradient(acc, bad_grad, cav_error)
      call check(error, allocated(cav_error), &
                 "mis-shaped gradient accumulator was accepted")

   end subroutine test_shape_guard

   !> Compare the reverse contraction against the forward arrays
   !>
   !> @param[in]  channels   Channel identifiers to populate
   !> @param[in]  label      Human-readable channel description
   !> @param[in]  branching  Use the branching fixture and require it to branch
   !> @param[out] error      Error handle
   subroutine run_equivalence(channels, label, branching, error)
      !> Channel identifiers to populate
      integer, intent(in) :: channels(:)
      !> Human-readable channel description
      character(len=*), intent(in) :: label
      !> Use the branching fixture and require it to branch
      logical, intent(in) :: branching
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error

      !> Reverse-mode and forward-mode gradients
      real(wp), allocatable :: grad_rev(:, :), grad_fwd(:, :)
      !> Grid extents and loop indices
      integer :: ngrid, nsph, igrid, iatom, iaxis, ich
      !> Comparison bookkeeping
      real(wp) :: diff, scale, worst
      integer :: worst_atom, worst_axis

      call build_fixture(cavity, ctx, branching, error)
      if (allocated(error)) return

      ngrid = cavity%ngrid
      nsph = cavity%nsph

      ! Without real branching the softmax reverse pass and its explicit owner
      ! term are dead code here, so the branching variant must prove it branched
      if (branching) then
         if (.not. allocated(cavity%branch_count)) then
            call test_failed(error, "branching fixture has no branch_count array")
            return
         end if
         if (.not. any(cavity%branch_count(1:ngrid) > 1)) then
            call test_failed(error, "branching fixture did not branch; "// &
                             "the softmax channel would go untested (ngrid "// &
                             to_string(ngrid)//", max branch_count "// &
                             to_string(maxval(cavity%branch_count(1:ngrid)))//")")
            return
         end if
      end if

      ! Forward path: fills every *_rA array the reference contraction reads
      call cavity%get_gradient(cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "forward gradient failed: "//cav_error%message)
         return
      end if
      if (.not. allocated(cavity%xyz1_rA) .or. .not. allocated(cavity%normal1_rA) &
          .or. .not. allocated(cavity%k1_rA)) then
         call test_failed(error, "forward derivative arrays are incomplete")
         return
      end if

      call acc%init(ngrid)
      do ich = 1, size(channels)
         call apply_channel(acc, channels(ich), cavity, error)
         if (allocated(error)) return
      end do

      allocate (grad_rev(3, nsph), source=0.0_wp)
      call cavity%get_surface_gradient(acc, grad_rev, cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "reverse gradient failed ("//label//"): "//cav_error%message)
         return
      end if

      ! Reference: contract the forward Jacobian with the same weights
      allocate (grad_fwd(3, nsph), source=0.0_wp)
      do igrid = 1, ngrid
         do iatom = 1, nsph
            do iaxis = 1, 3
               grad_fwd(iaxis, iatom) = grad_fwd(iaxis, iatom) &
                  & + acc%w_xi(igrid)*cavity%xi1_rA(iaxis, iatom, igrid) &
                  & + acc%w_f(igrid)*cavity%f1_rA(iaxis, iatom, igrid) &
                  & + acc%w_a(igrid)*cavity%a_i1_rA(iaxis, iatom, igrid) &
                  & + acc%w_w(igrid)*cavity%wleb1_rA(iaxis, iatom, igrid) &
                  & + acc%w_k1(igrid)*cavity%k1_rA(iaxis, iatom, igrid) &
                  & + acc%w_k2(igrid)*cavity%k2_rA(iaxis, iatom, igrid) &
                  & + dot_product(acc%w_xyz(:, igrid), cavity%xyz1_rA(:, iaxis, iatom, igrid)) &
                  & + dot_product(acc%w_n(:, igrid), cavity%normal1_rA(:, iatom, iaxis, igrid))
            end do
         end do
      end do

      ! Non-vacuity: a channel that produced nothing would pass trivially
      if (maxval(abs(grad_fwd)) <= VACUITY_THR) then
         call test_failed(error, "reference gradient is vacuous for "//label// &
                          " (max |g| = "//to_string(maxval(abs(grad_fwd)))//")")
         return
      end if

      worst = 0.0_wp
      worst_atom = 0
      worst_axis = 0
      do iatom = 1, nsph
         do iaxis = 1, 3
            diff = abs(grad_rev(iaxis, iatom) - grad_fwd(iaxis, iatom))
            scale = max(abs(grad_fwd(iaxis, iatom)), 1.0_wp)
            if (diff/scale > worst) then
               worst = diff/scale
               worst_atom = iatom
               worst_axis = iaxis
            end if
            if (diff > EQ_ABS .and. diff > EQ_REL*abs(grad_fwd(iaxis, iatom))) then
               call test_failed(error, "reverse/forward gradient mismatch for "//label// &
                                " at atom "//to_string(iatom)//" axis "//to_string(iaxis)// &
                                ": reverse "//to_string(grad_rev(iaxis, iatom))// &
                                " forward "//to_string(grad_fwd(iaxis, iatom)))
               return
            end if
         end do
      end do

   end subroutine run_equivalence

   !> Populate one surface-adjoint channel with a reproducible weight pattern
   !>
   !> The pattern is deterministic and varies across the grid so that a bug
   !> that happens to cancel for uniform weights still shows up.
   !>
   !> @param[inout] acc      Surface-adjoint accumulator
   !> @param[in]    channel  Channel identifier
   !> @param[in]    cavity   Cavity supplying the grid extent
   !> @param[out]   error    Error handle
   subroutine apply_channel(acc, channel, cavity, error)
      !> Surface-adjoint accumulator
      type(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Channel identifier
      integer, intent(in) :: channel
      !> Cavity supplying the grid extent
      type(cavity_type_drop), intent(in) :: cavity
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: add_error
      real(wp), allocatable :: ws(:), wv(:, :)
      integer :: ngrid, igrid, iaxis

      ngrid = cavity%ngrid
      allocate (ws(ngrid), wv(3, ngrid))
      do igrid = 1, ngrid
         ws(igrid) = ramp(igrid, 0)
         do iaxis = 1, 3
            wv(iaxis, igrid) = ramp(igrid, iaxis)
         end do
      end do

      select case (channel)
      case (CH_XI)
         call acc%add_surface_weights(add_error, w_xi=ws)
      case (CH_F)
         call acc%add_surface_weights(add_error, w_f=ws)
      case (CH_A)
         call acc%add_surface_weights(add_error, w_a=ws)
      case (CH_W)
         call acc%add_surface_weights(add_error, w_w=ws)
      case (CH_XYZ)
         call acc%add_surface_weights(add_error, w_xyz=wv)
      case (CH_N)
         call acc%add_surface_weights(add_error, w_n=wv)
      case (CH_K1)
         call acc%add_surface_weights(add_error, w_k1=ws)
      case (CH_K2)
         call acc%add_surface_weights(add_error, w_k2=ws)
      case default
         call test_failed(error, "unknown adjoint channel "//to_string(channel))
         return
      end select

      if (allocated(add_error)) then
         call test_failed(error, "add_surface_weights failed: "//add_error%message)
         return
      end if

   end subroutine apply_channel

   !> Deterministic pseudo-random weight in roughly [-1, 1]
   !>
   !> @param[in] igrid   Grid index
   !> @param[in] stride  Component offset, 0 for scalar channels
   !> @returns           Weight value
   pure real(wp) function ramp(igrid, stride)
      !> Grid index
      integer, intent(in) :: igrid
      !> Component offset
      integer, intent(in) :: stride

      ramp = sin(0.7_wp*real(igrid, wp) + 1.3_wp*real(stride, wp)) &
             + 0.3_wp*cos(0.21_wp*real(igrid, wp)*real(stride + 1, wp))

   end function ramp

   !> Build the shared DROP fixture with every optional property enabled
   !>
   !> The two heavy centres give an elongated cavity with well-separated
   !> principal curvatures and the off-axis hydrogen removes the residual
   !> rotational symmetry, so no channel is accidentally degenerate.
   !>
   !> The `branching` variant instead uses the five-carbon cross with concave
   !> seams. Branching needs the projector to find sibling minima, which takes
   !> a concave seam and multistart projection -- a plain dimer never branches
   !> at any Lebedev order or softmax scale -- *and* a softmax flat enough that
   !> the prune keeps the siblings. The caller asserts `branch_count > 1`
   !> actually occurred rather than trusting the configuration.
   !>
   !> @param[out]   cavity     Constructed cavity
   !> @param[inout] ctx        Run context borrowed by the cavity
   !> @param[in]    branching  Build the branching variant
   !> @param[out]   error      Error handle
   subroutine build_fixture(cavity, ctx, branching, error)
      !> Constructed cavity
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      !> Run context borrowed by the cavity; must outlive it
      type(moist_context_type), target, intent(inout) :: ctx
      !> Build the branching variant
      logical, intent(in) :: branching
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      call fixture_geometry(branching, mol)
      call build_cavity(cavity, ctx, mol, branching, error)

   end subroutine build_fixture

   !> Geometry of the requested fixture variant
   !>
   !> @param[in]  branching  Select the branching variant
   !> @param[out] mol        Fixture structure
   subroutine fixture_geometry(branching, mol)
      !> Select the branching variant
      logical, intent(in) :: branching
      !> Fixture structure
      type(structure_type), intent(out) :: mol

      if (branching) then
         call get_test_cross(mol)
      else
         call new(mol, [8, 6, 1], reshape([ &
                                          0.00_wp, 0.00_wp, 0.00_wp, &
                                          0.00_wp, 0.00_wp, 4.60_wp, &
                                          2.60_wp, 0.40_wp, -1.10_wp], [3, 3]))
      end if

   end subroutine fixture_geometry

   !> Build the DROP cavity for a given structure and fixture variant
   !>
   !> @param[out]   cavity     Constructed cavity
   !> @param[inout] ctx        Run context borrowed by the cavity
   !> @param[in]    mol        Structure to build on
   !> @param[in]    branching  Use the branching configuration
   !> @param[out]   error      Error handle
   subroutine build_cavity(cavity, ctx, mol, branching, error)
      !> Constructed cavity
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      !> Run context borrowed by the cavity; must outlive it
      type(moist_context_type), target, intent(inout) :: ctx
      !> Structure to build on
      type(structure_type), intent(in) :: mol
      !> Use the branching configuration
      logical, intent(in) :: branching
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: radii(:)
      type(mctc_error), allocatable :: cav_error
      integer :: iat, proj_level_loc, nleb_loc, prune_loc
      real(wp) :: blend_k_loc, gamma_loc

      nleb_loc = NUM_LEB
      prune_loc = 4
      gamma_loc = gamma
      blend_k_loc = k
      proj_level_loc = PROJ_LEVEL
      if (branching) then
         proj_level_loc = 7
         blend_k_loc = 1.0_wp
         gamma_loc = 1.0_wp
         nleb_loc = 110
         prune_loc = 0
      end if

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)), cav_error)
         if (allocated(cav_error)) then
            call test_failed(error, "radius lookup failed: "//trim(cav_error%message))
            return
         end if
      end do

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=blend_k_loc, blend_3b=gamma_loc)
         call new_context(ctx, verbosity=0)
         call new_cavity_drop(cavity, ctx, nleb=nleb_loc, &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                              proj_level=proj_level_loc, wleb_prune_level=prune_loc, &
                              radius_model=default_cpcm_radii(), &
                              lsf_model=svdw_template, error=cav_error)
      end block
      if (allocated(cav_error)) then
         call test_failed(error, "failed to initialize cavity: "//cav_error%message)
         return
      end if

      ! Multistart projection finds sibling branches, but the default softmax
      ! scale (0.0025) is so peaked that the prune in filter.f90 keeps only the
      ! strongest one and branch_count collapses to 1. Widening the softmax lets
      ! siblings survive, which is what puts the branch channel under test.
      if (branching) then
         cavity%param%branch_weight_s = BRANCH_SOFTMAX_S
         cavity%param%branch_rho_cut = log(1.0_wp/cavity%param%wleb_cut) &
                                       *sqrt(BRANCH_SOFTMAX_S)
         call cavity%branch_weight%init(BRANCH_SOFTMAX_S)
      end if

      ! do_fine turns on curvature, normals, r_iI and rho so that every
      ! optional *_rA array the reference contraction reads is allocated
      call cavity%properties(do_fine=.true.)

      call cavity%update(mol, error=cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "failed to build cavity: "//cav_error%message)
         return
      end if
      if (cavity%ngrid == 0) then
         call test_failed(error, "empty grid")
         return
      end if

   end subroutine build_cavity

end module test_cavity_drop_nuclear_adjoint
